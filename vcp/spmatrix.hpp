// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATRIX_HPP
#define VCP_SPMATRIX_HPP

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

// spmats.hpp already includes spmats_eigs_types.hpp which defines all type traits,
// enums, and result structs. No redefinition needed here.
#include <vcp/spmats.hpp>
// SPC-1: scalar conversion layer for the conversion constructors/assignments.
#include <vcp/vcp_converter.hpp>

namespace vcp {

	// SPC-3: forward declaration only.  spmatrix.hpp must NOT include
	// matrix.hpp (design SPC-3_design.md v1.1 S3-2, mirror of the SPC-2
	// forward declaration in matrix.hpp): every use of matrix in the
	// dense-to-sparse operator= below is inside a template, so name
	// resolution happens at instantiation time, and any translation unit
	// that instantiates it necessarily includes matrix.hpp itself.  No
	// default argument here -- matrix.hpp declares it.
	template <typename _T, class _P> class matrix;

	// SUB-1: block-write proxy (defined after spmatrix below).
	template <typename _T, class _P> class spmatrix_block;

	// SLU-C3-ELEMENT-ACCESSOR: lightweight proxy returned by the non-const
	// spmatrix::operator()(i,j), so that `T c = A(i,j)`, `A(i,j) = v`,
	// `A(i,j) += v`, `A(i,j) *= v` etc. read naturally while still routing
	// through get()/set()/add() (Option A tagged merge, see
	// sandbox/docs/design/spmats_finalize_policy.md and
	// sandbox/docs/issues/VCP_task_list.md C-3). += routes straight to
	// add(v); -= routes to add(-v) (subtraction is just signed addition, so
	// it rides the same tagged-merge fast path -- SLU-C3-FIX-COMPOUND-MINUS).
	// Neither needs a get() round trip. *= and /= have no such fast path
	// (the current value must be known before the new value can be
	// computed), so they stay on get -> compute -> set.
	template <typename _T, class _Owner>
	class spmats_element_proxy {
	public:
		typedef typename _Owner::index_type index_type;

		spmats_element_proxy(_Owner& owner, const index_type i, const index_type j)
			: owner_(owner), i_(i), j_(j) {}

		operator _T() const { return owner_.get(i_, j_); }

		// Explicit proxy-to-proxy assignment: without this, the compiler's
		// implicitly-deleted copy-assignment (exact-match, no conversion
		// needed since _Owner& has no copy assignment) beats the operator=
		// (const _T&) template below in overload resolution and makes
		// `A(i,j) = A(k,l)` ill-formed.
		spmats_element_proxy& operator=(const spmats_element_proxy& other) {
			owner_.set(i_, j_, static_cast<_T>(other));
			return *this;
		}
		spmats_element_proxy& operator=(const _T& v) {
			owner_.set(i_, j_, v);
			return *this;
		}
		spmats_element_proxy& operator+=(const _T& v) {
			owner_.add(i_, j_, v);
			return *this;
		}
		spmats_element_proxy& operator-=(const _T& v) {
			owner_.add(i_, j_, -v);
			return *this;
		}
		spmats_element_proxy& operator*=(const _T& v) {
			owner_.set(i_, j_, owner_.get(i_, j_) * v);
			return *this;
		}
		spmats_element_proxy& operator/=(const _T& v) {
			owner_.set(i_, j_, owner_.get(i_, j_) / v);
			return *this;
		}

	private:
		_Owner& owner_;
		index_type i_;
		index_type j_;
	};

	template <typename _T, class _P = spmats<_T> > class spmatrix : protected _P {
	public:
		typedef _T value_type;
		typedef _P policy_type;
		typedef typename _P::index_type index_type;
		typedef typename _P::format_type format_type;
		typedef typename spmatrix_real_type<_T>::type scalar_real_type;
		typedef typename _P::dense_matrix_type dense_matrix_type;
		typedef vcp::linear_solve_options<_T> linear_solve_options_type;
		typedef vcp::eig_options<_T> eig_options_type;
		typedef vcp::ldl_options<_T> ldl_options_type;
		typedef vcp::ldl_result<_T, typename _P::index_type> ldl_result_type;
		typedef vcp::inertia_options<_T> inertia_options_type;
		typedef vcp::inertia_result<typename _P::index_type> inertia_result_type;
		typedef vcp::lu_extract_options<_T> lu_extract_options_type;
		typedef vcp::lu_extract_result<_T, typename _P::index_type> lu_extract_result_type;
		typedef vcp::chol_options<_T> chol_options_type;
		typedef vcp::chol_result<_T, typename _P::index_type> chol_result_type;
		typedef vcp::ainv_options<_T> ainv_options_type;
		typedef vcp::ainv_result<_T, typename _P::index_type> ainv_result_type;
		typedef vcp::fsai_options<_T> fsai_options_type;
		typedef vcp::fsai_result<_T, typename _P::index_type> fsai_result_type;
		typedef vcp::fsai_adaptive_options<_T> fsai_adaptive_options_type;
		typedef vcp::fsai_adaptive_result<_T, typename _P::index_type> fsai_adaptive_result_type;
		// LSS-1 P-4
		typedef vcp::sparse_lu_options<_T> sparse_lu_options_type;
		typedef vcp::lu_factor_handle<_T, typename _P::index_type> lu_factor_handle_type;
		// SLDL-SH (B2)
		typedef vcp::ldl_shift_handle<_T, typename _P::index_type> ldl_shift_handle_type;
		typedef typename ldl_shift_handle_type::workspace_type ldl_shift_workspace_type;
		typedef typename ldl_shift_handle_type::inertia_result_type ldl_shift_inertia_result_type;

		spmatrix() : _P() {}
		spmatrix(const index_type rows, const index_type cols) : _P() { this->resize(rows, cols); }
		~spmatrix() = default;
		spmatrix(const spmatrix&) = default;
		spmatrix(spmatrix&&) = default;
		spmatrix& operator=(const spmatrix&) = default;
		spmatrix& operator=(spmatrix&&) = default;

		// ---------------------------------------------------------------
		// SPC-1: conversion constructors / conversion assignments
		// (sparse mirror of matrix.hpp L57-99; design SPC-1_design.md v2.2).
		// Signatures and SFINAE conditions are identical in form to the
		// dense members, so the identity instantiation (same _T AND same
		// _P) is excluded and falls back to the defaulted copy members
		// above -- existing paths stay bit-identical.  All four members
		// delegate to convert_from_ below (finalize policy category 4:
		// read through an as_csr()/as_csc() copy, source A is never
		// modified, result is born-finalized, format is preserved).
		// Note: in expressions these create an O(nnz) temporary, same
		// trade-off as the dense matrix conversion members.
		// ---------------------------------------------------------------

		// (1) Policy conversion constructor: T same, P different
		template <class _P2,
		          typename std::enable_if<!std::is_same<_P, _P2>::value, int>::type = 0>
		spmatrix(const spmatrix<_T, _P2>& A) : _P() {
			convert_from_(A);
		}

		// (2) Type conversion constructor: T different, P arbitrary
		template <typename _T2, class _P2,
		          typename std::enable_if<!std::is_same<_T, _T2>::value, int>::type = 0>
		spmatrix(const spmatrix<_T2, _P2>& A) : _P() {
			convert_from_(A);
		}

		// (3) Policy conversion assignment: T same, P different
		template <class _P2,
		          typename std::enable_if<!std::is_same<_P, _P2>::value, int>::type = 0>
		spmatrix<_T, _P>& operator=(const spmatrix<_T, _P2>& A) {
			this->clear();
			convert_from_(A);
			return *this;
		}

		// (4) Type conversion assignment: T different, P arbitrary
		template <typename _T2, class _P2,
		          typename std::enable_if<!std::is_same<_T, _T2>::value, int>::type = 0>
		spmatrix<_T, _P>& operator=(const spmatrix<_T2, _P2>& A) {
			this->clear();
			convert_from_(A);
			return *this;
		}

		// ---------------------------------------------------------------
		// SPC-3: dense-to-sparse conversion assignment
		// (design SPC-3_design.md v1.1 S3-1..S3-7).
		// Stores the components that are strictly nonzero in the source
		// AND remain strictly nonzero after value conversion (zero rule
		// S3-3': no thresholding; an interval is zero only if both
		// endpoints are exactly zero, so [-eps, eps] is kept).  Values
		// that become exactly zero through conversion (interval -> point
		// mid = 0, underflow to 0 on narrowing) are dropped, per the
		// spmats invariant "explicit zero is not allowed" (assign_* would
		// throw on them).  Because -0.0 == 0.0, a -0.0 source component
		// is skipped and the sign of zero is not preserved (unstored
		// elements read back as +0.0).
		// The dense side is read through the public 2-arg
		// operator()(i, j) const and rowsize()/columnsize() only (S3-5;
		// correct for all matstype states, and 'N' has
		// row == column == 0, so it naturally yields a 0x0 empty result,
		// S3-6).  A is never copied or modified.  Row-major single pass:
		// the column indices of each row are ascending and unique by
		// construction, no zero value is ever pushed, so the assign_csr
		// preconditions (sorted / unique / no explicit zero) hold
		// constructively and the result is born-finalized CSR.
		// Complexity O(rows*cols) (the lower bound for a dense input;
		// push_back is amortized O(1)).  Only fires on an explicit
		// assignment statement; there is no converting constructor, so
		// no implicit conversion path exists (S3-1, negative check
		// SPC3-T7).
		// ---------------------------------------------------------------
		template <typename _T2, class _P2>
		spmatrix<_T, _P>& operator=(const matrix<_T2, _P2>& A) {
			const index_type rows = static_cast<index_type>(A.rowsize());
			const index_type cols = static_cast<index_type>(A.columnsize());
			std::vector<index_type> row_ptr(static_cast<std::size_t>(rows) + 1, 0);
			std::vector<index_type> col_idx;
			std::vector<_T> val;
			for (index_type i = 0; i < rows; ++i) {
				for (index_type j = 0; j < cols; ++j) {
					const _T2& x = A(static_cast<int>(i), static_cast<int>(j));
					if (is_strict_zero_(x)) continue;    // (i) source-side test
					_T y;
					convert_value_(x, y);
					if (is_strict_zero_(y)) continue;    // (ii) post-conversion test
					col_idx.push_back(j);
					val.push_back(y);
				}
				row_ptr[static_cast<std::size_t>(i) + 1] = static_cast<index_type>(col_idx.size());
			}
			this->assign_csr(rows, cols, row_ptr, col_idx, val);
			return *this;
		}

	private:
		// SPC-1: strict-zero test on the DESTINATION value type (design
		// v2.2 §2): point types compare against _Tv(0); interval types
		// compare endpoints only (no abs / three-way certified compare
		// needed -- this mirrors the [0,0]-only semantics of the spmats
		// invariant "explicit zero is not allowed").
		template <typename _Tv>
		static bool is_strict_zero_(const _Tv& x) {
			return x == _Tv(0);
		}
#if defined(INTERVAL_HPP)
		template <typename _Tv>
		static bool is_strict_zero_(const kv::interval<_Tv>& x) {
			return x.lower() == _Tv(0) && x.upper() == _Tv(0);
		}
#endif

		// SPC-3: value-conversion helper for the dense-to-sparse
		// operator= above (directive §4, same tag dispatch as the SPC-2
		// helper in matrix.hpp).  Same scalar type: plain assignment
		// (does not rely on an identity vcp::convert overload).
		// Different scalar type: delegate to the vcp::convert scalar
		// layer (point->point widening exact / narrowing nearest,
		// ->interval inclusion-preserving, interval->point mid).
		template <typename _S>
		static void convert_value_(const _S& x, _S& y) { y = x; }
		template <typename _S2, typename _S,
		          typename std::enable_if<!std::is_same<_S, _S2>::value, int>::type = 0>
		static void convert_value_(const _S2& x, _S& y) { vcp::convert(x, y); }

		// SPC-1 common core for the four conversion members above.
		// Complexity O(nnz + max(rows, cols)).  Access to A is via the
		// public spmatrix API only; the destructive helpers (to_csr etc.)
		// are never called on A itself -- as_csr()/as_csc() convert on a
		// copy, which also normalizes an unfinalized (COO) source.
		template <typename _T2, class _P2>
		void convert_from_(const spmatrix<_T2, _P2>& A) {
			typedef typename spmatrix<_T2, _P2>::index_type src_index_type;
			// Format preservation (design v2.2 §2.1): finalized CSC stays
			// CSC; finalized CSR and unfinalized input become CSR.
			const bool use_csc = A.is_finalized() && A.format() == vcp::sparse_csc;
			const spmatrix<_T2, _P2> src = use_csc ? A.as_csc() : A.as_csr();
			const index_type rows = static_cast<index_type>(src.rowsize());
			const index_type cols = static_cast<index_type>(src.columnsize());
			const std::size_t nouter = static_cast<std::size_t>(use_csc ? cols : rows);
			const std::vector<src_index_type>& src_outer = src.outer_index();
			const std::vector<src_index_type>& src_inner = src.inner_index();
			const std::vector<_T2>& src_value = src.values();
			// Index transfer is a per-element static_cast: if the
			// destination index_type differs (future policies), every
			// stored index is < rows/cols, so whenever assign_* succeeds
			// for (rows, cols) all indices are representable in it.
			std::vector<index_type> outer(nouter + 1);
			std::vector<index_type> inner(src_value.size());
			std::vector<_T> val(src_value.size());
			// Single pass: convert each stored value (vcp::convert; for
			// the same-_T policy conversions this is the identity
			// overload, i.e. a plain copy) and drop values that became
			// strictly zero (v2.2 zero-drop pass: the spmats invariant
			// forbids explicit zeros, assign_* would reject them).  The
			// outer pointers are rebuilt on the fly with a second write
			// cursor, keeping the whole pass O(nnz + max(rows, cols)).
			std::size_t out = 0;
			outer[0] = 0;
			for (std::size_t i = 0; i < nouter; i++) {
				for (src_index_type k = src_outer[i];
				     k < src_outer[i + 1]; k++) {
					_T y;
					vcp::convert(src_value[static_cast<std::size_t>(k)], y);
					if (!is_strict_zero_(y)) {
						inner[out] = static_cast<index_type>(src_inner[static_cast<std::size_t>(k)]);
						val[out] = y;
						out++;
					}
				}
				outer[i + 1] = static_cast<index_type>(out);
			}
			inner.resize(out);
			val.resize(out);
			// assign_csr/assign_csc postconditions: born-finalized,
			// sorted, unique (arrays are already zero-free and sorted).
			if (use_csc) {
				this->assign_csc(rows, cols, outer, inner, val);
			}
			else {
				this->assign_csr(rows, cols, outer, inner, val);
			}
		}

	public:
		index_type rowsize() const { return _P::rowsize(); }
		index_type columnsize() const { return _P::columnsize(); }
		index_type nnz() const { return _P::nnz(); }
		index_type stored_nnz() const { return _P::stored_nnz(); }
		bool is_finalized() const { return _P::is_finalized(); }
		bool is_sorted() const { return _P::is_sorted(); }
		bool is_unique() const { return _P::is_unique(); }
		format_type format() const { return _P::format(); }

		void resize(const index_type rows, const index_type cols) { _P::resize(rows, cols); }
		void clear() { _P::clear(); }
		void reserve(const index_type n) { _P::reserve(n); }

		void add(const index_type i, const index_type j, const _T& value) { _P::add(i, j, value); }
		void set(const index_type i, const index_type j, const _T& value) { _P::set(i, j, value); }
		_T get(const index_type i, const index_type j) const { return _P::get(i, j); }

		// SLU-C3-ELEMENT-ACCESSOR: `T c = A(i,j)`, `A(i,j) = v`,
		// `A(i,j) += v`, `A(i,j) *= v` (also -=, /=) via spmats_element_proxy.
		// finalize() is never forced (same as get()).
		typedef spmats_element_proxy<_T, spmatrix> element_proxy_type;
		element_proxy_type operator()(const index_type i, const index_type j) {
			return element_proxy_type(*this, i, j);
		}
		_T operator()(const index_type i, const index_type j) const { return get(i, j); }

		// -----------------------------------------------------------------
		// SUB-1: submatrix extraction (dense matrix::submatrix parity).
		// Selector grammar (0-based, exact mirror of mats::submat):
		//   {}      whole axis (MATLAB ':')
		//   {i}     single index
		//   {a,b}   CLOSED range a..b (both ends included -- NOT half-open)
		//   {a,s,b} stride a, a+s, ... while <= b (b need not be hit)
		// Invalid selectors throw vcp::index_error with the dense-side
		// conditions and message form.  The result is a new value,
		// born-finalized, source format preserved (finalized CSC -> CSC;
		// finalized CSR and unfinalized -> CSR); *this is never modified.
		// Return-type rule (S11): bare ints only = element (T, A(i,j));
		// ANY braced selector = submatrix: A({1},{2}), A(1,{2}), A({1},2)
		// are 1x1 spmatrix, A(1,2) is the element.
		// -----------------------------------------------------------------
		spmatrix submatrix(const std::initializer_list<int>& list1,
		                   const std::initializer_list<int>& list2) const {
			spmatrix B;
			_P::submat(static_cast<_P&>(B), list1, list2);
			return B;
		}

		// const sugar -- pure forwarding to submatrix (extraction)
		spmatrix operator()(const std::initializer_list<int>& list1,
		                    const std::initializer_list<int>& list2) const {
			return this->submatrix(list1, list2);
		}
		// mixed forms (S11): the int is wrapped as the single selector {i}
		spmatrix operator()(const int i, const std::initializer_list<int>& list2) const {
			return this->submatrix({i}, list2);
		}
		spmatrix operator()(const std::initializer_list<int>& list1, const int j) const {
			return this->submatrix(list1, {j});
		}

		// SUB-1 Phase W: the non-const forms return the block-write proxy
		// (spmatrix_block, defined after this class): A({..},{..}) = B
		// (same-type spmatrix, block REPLACEMENT) or = scalar (fill).
		// Selector validation happens HERE, at proxy construction.
		// See spmatrix_block for the semantics (zeros in the right-hand
		// side delete stored entries; nonzero scalar fill DENSIFIES the
		// block) and the lifetime warning (`auto x = A({..},{..});`
		// captures the proxy, not an spmatrix).
		spmatrix_block< _T, _P > operator()(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2);
		spmatrix_block< _T, _P > operator()(const int i, const std::initializer_list<int>& list2);
		spmatrix_block< _T, _P > operator()(const std::initializer_list<int>& list1, const int j);

		void sort_coo() { _P::sort_coo(); }
		void normalize_coo() { _P::normalize_coo(); }
		void finalize() { _P::finalize(); }
		void to_csr() { _P::to_csr(); }
		void to_csc() { _P::to_csc(); }

		spmatrix as_csr() const {
			spmatrix B;
			static_cast<vcp::spmats<_T, typename _P::index_type>&>(static_cast<_P&>(B)) = _P::as_csr();
			// TPX-1: assign into the base spmats subobject so derived policies
			// (spumar, spcmodumar, ...) compile.  The result's policy-specific state is
			// default-constructed -- the same convention as every other
			// result-returning operation of spmatrix (T-2).
			return B;
		}

		spmatrix as_csc() const {
			spmatrix B;
			static_cast<vcp::spmats<_T, typename _P::index_type>&>(static_cast<_P&>(B)) = _P::as_csc();
			// TPX-1: assign into the base spmats subobject so derived policies
			// (spumar, spcmodumar, ...) compile.  The result's policy-specific state is
			// default-constructed -- the same convention as every other
			// result-returning operation of spmatrix (T-2).
			return B;
		}

		std::vector<_T> mul_vec(const std::vector<_T>& x) const { return _P::mul_vec(x); }
		void mul_vec(const _T* x, _T* y) const { _P::mul_vec(x, y); }
		std::vector<_T> trans_mul_vec(const std::vector<_T>& x) const { return _P::trans_mul_vec(x); }
		void trans_mul_vec(const _T* x, _T* y) const { _P::trans_mul_vec(x, y); }

		// Phase 7.7: to_dense / is_symmetric are public API; implementation is
		// delegated to policy_to_dense / policy_is_symmetric so that a custom
		// policy P can override the semantics (e.g. verified dense enclosure).
		dense_matrix_type to_dense() const {
			return this->policy_to_dense(*this);
		}

		bool is_symmetric() const {
			return this->policy_is_symmetric(*this);
		}

		bool is_symmetric(const scalar_real_type& tol) const {
			return this->policy_is_symmetric(*this, tol);
		}

		// ---------------------------------------------------------------
		// Named arithmetic methods — all delegate to policy
		// These remain for backward compatibility and as building blocks
		// for the friend operators below.
		// ---------------------------------------------------------------
		spmatrix add(const spmatrix& rhs) const {
			if (rowsize() != rhs.rowsize() || columnsize() != rhs.columnsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::add: dimension mismatch");
			spmatrix C;
			static_cast<_P&>(C) = this->policy_add(
				static_cast<const _P&>(*this), static_cast<const _P&>(rhs));
			return C;
		}

		spmatrix sub(const spmatrix& rhs) const {
			if (rowsize() != rhs.rowsize() || columnsize() != rhs.columnsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::sub: dimension mismatch");
			spmatrix C;
			static_cast<_P&>(C) = this->policy_sub(
				static_cast<const _P&>(*this), static_cast<const _P&>(rhs));
			return C;
		}

		spmatrix matmul(const spmatrix& rhs) const {
			if (columnsize() != rhs.rowsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::matmul: dimension mismatch");
			spmatrix C;
			static_cast<_P&>(C) = this->policy_mul(
				static_cast<const _P&>(*this), static_cast<const _P&>(rhs));
			return C;
		}

		// ---------------------------------------------------------------
		// In-place arithmetic — delegate to policy; aligned with matrix<T,P>
		// naming (addmm, subsmmA, subsmmB, mulmm, mulsm, mulms, divms, minusm)
		// ---------------------------------------------------------------

		// *this = *this + B
		void addmm(const spmatrix& B) {
			if (rowsize() != B.rowsize() || columnsize() != B.columnsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::addmm: dimension mismatch");
			static_cast<_P&>(*this) = this->policy_add(
				static_cast<const _P&>(*this), static_cast<const _P&>(B));
		}

		// *this = *this - B
		void subsmmA(const spmatrix& B) {
			if (rowsize() != B.rowsize() || columnsize() != B.columnsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::subsmmA: dimension mismatch");
			static_cast<_P&>(*this) = this->policy_sub(
				static_cast<const _P&>(*this), static_cast<const _P&>(B));
		}

		// *this = A - *this
		void subsmmB(const spmatrix& A) {
			if (A.rowsize() != rowsize() || A.columnsize() != columnsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::subsmmB: dimension mismatch");
			static_cast<_P&>(*this) = this->policy_sub(
				static_cast<const _P&>(A), static_cast<const _P&>(*this));
		}

		// C = *this * B  (out-of-place; avoids aliasing issues)
		void mulmm(const spmatrix& B, spmatrix& C) const {
			if (columnsize() != B.rowsize())
				vcp::throw_error<vcp::dimension_error>("spmatrix::mulmm: dimension mismatch");
			static_cast<_P&>(C) = this->policy_mul(
				static_cast<const _P&>(*this), static_cast<const _P&>(B));
		}

		// *this = alpha * *this
		// SPC-P1: calls the destructive policy directly (in-place, no copy).
		void mulsm(const _T& alpha) {
			static_cast<_P&>(*this).policy_mulsm(alpha);
		}

		// *this = *this * alpha  (commutativity: same as mulsm)
		void mulms(const _T& alpha) { mulsm(alpha); }

		// *this = *this / alpha
		// SPC-P1: calls the destructive policy directly (in-place, no copy).
		void divms(const _T& alpha) {
			static_cast<_P&>(*this).policy_divms(alpha);
		}

		// *this = -*this
		// SPC-P1: calls the destructive policy directly (in-place, no copy).
		void minusm() {
			static_cast<_P&>(*this).policy_minusm();
		}

		// ---------------------------------------------------------------
		// MATLAB-like friend operators — A+B, A-B, A*B with rvalue overloads
		// All computation is delegated to the named methods above, which
		// in turn call policy methods.  Verified policy can override any
		// policy_* method to change the arithmetic semantics.
		// ---------------------------------------------------------------

		// --- A + B ---
		friend spmatrix operator+(const spmatrix& A, const spmatrix& B) {
			return A.add(B);
		}
		friend spmatrix operator+(spmatrix&& A, const spmatrix& B) {
			A.addmm(B);
			return std::move(A);
		}
		friend spmatrix operator+(const spmatrix& A, spmatrix&& B) {
			B.addmm(A);
			return std::move(B);
		}
		friend spmatrix operator+(spmatrix&& A, spmatrix&& B) {
			A.addmm(B);
			return std::move(A);
		}

		// --- A - B ---
		friend spmatrix operator-(const spmatrix& A, const spmatrix& B) {
			return A.sub(B);
		}
		friend spmatrix operator-(spmatrix&& A, const spmatrix& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		friend spmatrix operator-(const spmatrix& A, spmatrix&& B) {
			B.subsmmB(A);
			return std::move(B);
		}
		friend spmatrix operator-(spmatrix&& A, spmatrix&& B) {
			A.subsmmA(B);
			return std::move(A);
		}

		// --- A * B (matrix-matrix) ---
		// Sparse matrix product cannot safely reuse either operand's storage
		// (output pattern differs), so rvalue overloads just call matmul.
		friend spmatrix operator*(const spmatrix& A, const spmatrix& B) {
			return A.matmul(B);
		}
		friend spmatrix operator*(spmatrix&& A, const spmatrix& B) {
			return A.matmul(B);
		}
		friend spmatrix operator*(const spmatrix& A, spmatrix&& B) {
			return A.matmul(B);
		}
		friend spmatrix operator*(spmatrix&& A, spmatrix&& B) {
			return A.matmul(B);
		}

		// --- scalar * A, A * scalar (Phase 7.8: constructible scalar template) ---
		// _Sm must be constructible to _T (std::is_constructible, not is_convertible).
		// This allows explicit constructors such as kv::dd(int) or kv::mpfr<N>(int).
		// The !std::is_same exclusion prevents collision with spmatrix * spmatrix overloads.
		// Internal conversion _T(_Sm) is explicit and performed exactly once.
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value &&
			!std::is_same<typename std::decay<_Sm>::type, spmatrix>::value,
			spmatrix
		>::type
		operator*(const _Sm& alpha, const spmatrix& A) {
			_T Ta = _T(alpha);
			spmatrix C = A;
			C.mulsm(Ta);
			return C;
		}
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value &&
			!std::is_same<typename std::decay<_Sm>::type, spmatrix>::value,
			spmatrix
		>::type
		operator*(const _Sm& alpha, spmatrix&& A) {
			_T Ta = _T(alpha);
			A.mulsm(Ta);
			return std::move(A);
		}
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value &&
			!std::is_same<typename std::decay<_Sm>::type, spmatrix>::value,
			spmatrix
		>::type
		operator*(const spmatrix& A, const _Sm& alpha) {
			_T Ta = _T(alpha);
			spmatrix C = A;
			C.mulsm(Ta);
			return C;
		}
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value &&
			!std::is_same<typename std::decay<_Sm>::type, spmatrix>::value,
			spmatrix
		>::type
		operator*(spmatrix&& A, const _Sm& alpha) {
			_T Ta = _T(alpha);
			A.mulsm(Ta);
			return std::move(A);
		}

		// --- A / scalar ---
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value,
			spmatrix
		>::type
		operator/(const spmatrix& A, const _Sm& alpha) {
			_T Ta = _T(alpha);
			spmatrix C = A;
			C.divms(Ta);
			return C;
		}
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value,
			spmatrix
		>::type
		operator/(spmatrix&& A, const _Sm& alpha) {
			_T Ta = _T(alpha);
			A.divms(Ta);
			return std::move(A);
		}

		// --- unary - ---
		friend spmatrix operator-(const spmatrix& A) {
			spmatrix C = A;
			C.minusm();
			return C;
		}
		friend spmatrix operator-(spmatrix&& A) {
			A.minusm();
			return std::move(A);
		}

		// --- compound assignment ---
		friend spmatrix& operator+=(spmatrix& A, const spmatrix& B) {
			A.addmm(B);
			return A;
		}
		friend spmatrix& operator-=(spmatrix& A, const spmatrix& B) {
			A.subsmmA(B);
			return A;
		}
		friend spmatrix& operator*=(spmatrix& A, const spmatrix& B) {
			spmatrix C;
			A.mulmm(B, C);
			A = std::move(C);
			return A;
		}
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value &&
			!std::is_same<typename std::decay<_Sm>::type, spmatrix>::value,
			spmatrix&
		>::type
		operator*=(spmatrix& A, const _Sm& alpha) {
			_T Ta = _T(alpha);
			A.mulsm(Ta);
			return A;
		}
		template <typename _Sm>
		friend typename std::enable_if<
			std::is_constructible<_T, _Sm>::value,
			spmatrix&
		>::type
		operator/=(spmatrix& A, const _Sm& alpha) {
			_T Ta = _T(alpha);
			A.divms(Ta);
			return A;
		}

		// ---------------------------------------------------------------
		// Linear system solve — all delegate to policy_lss / policy_lss_with_info
		// ---------------------------------------------------------------

		// strict solve: policy decides convergence checking
		std::vector<_T> solve(const std::vector<_T>& b, const linear_solve_options_type& options = linear_solve_options_type()) const {
			return this->policy_lss(b, options);
		}

		// non-strict: return full diagnostic result
		linear_solve_result<_T> solve_with_info(const std::vector<_T>& b, const linear_solve_options_type& options = linear_solve_options_type()) const {
			return this->policy_lss_with_info(b, options);
		}

		// LSS-1 P-4: reusable LU factorization handle (thin forwarding only,
		// WFIX layer discipline) — delegates to policy_lu_factorize_with_info.
		lu_factor_handle_type lu_factorize_with_info(
			const sparse_lu_options_type& opt = sparse_lu_options_type()) const {
			return this->policy_lu_factorize_with_info(opt);
		}

		// ---------------------------------------------------------------
		// LDL^T factorization (LDL-3) — delegates to policy_ldl_with_info.
		// Convention P^T A P = L D L^T with perm p new->old (A(p,p) =
		// L D L^T, MATLAB 'vector' form) and P(p[k],k) = 1 (design v2
		// SS5.3-5.4).  strict ldl: ANY status != success throws
		// (zero_pivot / not_symmetric / inconclusive_pivot_test /
		// structural_singularity included); use ldl_with_info to inspect
		// such factorizations (decision 2).
		// ---------------------------------------------------------------

		// non-strict, permutation-vector form
		ldl_result_type ldl_with_info(spmatrix& L, spmatrix& D, std::vector<index_type>& p,
		                              const ldl_options_type& options = ldl_options_type()) const {
			return this->policy_ldl_with_info(
				static_cast<_P&>(L), static_cast<_P&>(D), p, options);
		}

		// non-strict, permutation-matrix form (P finalized, P(p[k],k) = 1;
		// WFIX: materialization lives in the policy-layer matrix-form
		// overload -- this wrapper only forwards)
		ldl_result_type ldl_with_info(spmatrix& L, spmatrix& D, spmatrix& P,
		                              const ldl_options_type& options = ldl_options_type()) const {
			return this->policy_ldl_with_info(
				static_cast<_P&>(L), static_cast<_P&>(D),
				static_cast<_P&>(P), options);
		}

		// strict, permutation-vector form
		void ldl(spmatrix& L, spmatrix& D, std::vector<index_type>& p,
		         const ldl_options_type& options = ldl_options_type()) const {
			const ldl_result_type result = ldl_with_info(L, D, p, options);
			if (result.status != sparse_ldl_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::ldl: factorization failed with status ",
					sparse_ldl_status_to_string(result.status));
			}
		}

		// strict, permutation-matrix form
		void ldl(spmatrix& L, spmatrix& D, spmatrix& P,
		         const ldl_options_type& options = ldl_options_type()) const {
			const ldl_result_type result = ldl_with_info(L, D, P, options);
			if (result.status != sparse_ldl_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::ldl: factorization failed with status ",
					sparse_ldl_status_to_string(result.status));
			}
		}

		// ---------------------------------------------------------------
		// LL^T Cholesky factorization (CHOL-3) — delegates to
		// policy_chol_with_info.  Convention P^T A P = L L^T with perm p
		// new->old (A(p,p) = L L^T, MATLAB chol(A,'lower','vector') form)
		// and P(p[k],k) = 1 (chol design v1 SS1.1; same orientation as
		// ldl).  L is NON-unit lower triangular with positive diagonal;
		// there is no D factor and no pivoting (p is the ordering output).
		// L / p (P) are valid outputs only when the status is success
		// (D-3: no partial factor).  strict chol: ANY status != success
		// throws (not_positive_definite / inconclusive_pivot_test /
		// not_symmetric included); use chol_with_info to inspect such
		// inputs (D-4).
		// ---------------------------------------------------------------

		// non-strict, permutation-vector form
		chol_result_type chol_with_info(spmatrix& L, std::vector<index_type>& p,
		                                const chol_options_type& options = chol_options_type()) const {
			return this->policy_chol_with_info(
				static_cast<_P&>(L), p, options);
		}

		// non-strict, permutation-matrix form (P finalized, P(p[k],k) = 1;
		// materialization lives in the policy-layer matrix-form overload --
		// this wrapper only forwards)
		chol_result_type chol_with_info(spmatrix& L, spmatrix& P,
		                                const chol_options_type& options = chol_options_type()) const {
			return this->policy_chol_with_info(
				static_cast<_P&>(L), static_cast<_P&>(P), options);
		}

		// strict, permutation-vector form
		void chol(spmatrix& L, std::vector<index_type>& p,
		          const chol_options_type& options = chol_options_type()) const {
			const chol_result_type result = chol_with_info(L, p, options);
			if (result.status != sparse_chol_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::chol: factorization failed with status ",
					sparse_chol_status_to_string(result.status));
			}
		}

		// strict, permutation-matrix form
		void chol(spmatrix& L, spmatrix& P,
		          const chol_options_type& options = chol_options_type()) const {
			const chol_result_type result = chol_with_info(L, P, options);
			if (result.status != sparse_chol_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::chol: factorization failed with status ",
					sparse_chol_status_to_string(result.status));
			}
		}

		// ---------------------------------------------------------------
		// inertia (LDL-4) — delegates to policy_inertia_with_info (default
		// implementation: LDL through the public policy API + certified scan
		// of the block diagonal D; design v2 SS7).  strict inertia throws on
		// any status != success.
		// ---------------------------------------------------------------

		// non-strict
		inertia_result_type inertia_with_info(const inertia_options_type& options = inertia_options_type()) const {
			return this->policy_inertia_with_info(options);
		}

		// strict
		inertia_result_type inertia(const inertia_options_type& options = inertia_options_type()) const {
			const inertia_result_type result = inertia_with_info(options);
			if (result.status != inertia_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::inertia: computation failed with status ",
					inertia_status_to_string(result.status));
			}
			return result;
		}

		// ---------------------------------------------------------------
		// A - sigma*B LDL^T shift handle (SLDL-SH / B2) — delegates to
		// policy_ldl_shift_setup_with_info.  Setup runs every
		// sigma-independent stage once; the returned read-only handle then
		// serves inertia_at(sigma [, ws]) (primary, H-5) without
		// constructing L / D / perm.  The overload without B iterates
		// A - sigma*I.  Defaults (H-1): method = supernodal,
		// pivoting = none (ldl_shift_default_options; the one-shot ldl
		// defaults are unchanged, H-6).  Info-only reporting through the
		// handle status (H-2): no strict variant exists.
		// ---------------------------------------------------------------

		ldl_shift_handle_type ldl_shift_setup_with_info(
			const ldl_options_type& options = vcp::ldl_shift_default_options<_T>()) const {
			return this->policy_ldl_shift_setup_with_info(options);
		}

		ldl_shift_handle_type ldl_shift_setup_with_info(const spmatrix& B,
			const ldl_options_type& options = vcp::ldl_shift_default_options<_T>()) const {
			return this->policy_ldl_shift_setup_with_info(B, options);
		}

		// ---------------------------------------------------------------
		// LU factor extraction (LUX-1) — delegates to policy_lu_with_info.
		// Convention (SSC, shared with spumar): P A Q = L U with p / q
		// new->old (A(p,q) = L U, MATLAB [L,U,P,Q] = lu(A) orientation)
		// and P(k,p[k]) = 1, Q(q[k],k) = 1.  NOTE: the row side is
		// TRANSPOSED relative to the LDL convention (LDL: P(p[k],k) = 1).
		// L unit lower (explicit unit diagonal), U upper triangular.
		// equilibration == true is rejected (unsupported_options); strict
		// lu throws on ANY status != success.  L / U / p / q (P / Q) are
		// valid outputs only on success.
		// ---------------------------------------------------------------

		// non-strict, permutation-vector form
		lu_extract_result_type lu_with_info(spmatrix& L, spmatrix& U,
		                                    std::vector<index_type>& p, std::vector<index_type>& q,
		                                    const lu_extract_options_type& options = lu_extract_options_type()) const {
			return this->policy_lu_with_info(
				static_cast<_P&>(L), static_cast<_P&>(U), p, q, options);
		}

		// non-strict, permutation-matrix form (P, Q finalized;
		// P(k,p[k]) = 1, Q(q[k],k) = 1; WFIX: materialization lives in the
		// policy-layer NVI (policy_lu_matrices_with_info_impl is the
		// backend replacement point) -- this wrapper only forwards)
		lu_extract_result_type lu_with_info(spmatrix& L, spmatrix& U, spmatrix& P, spmatrix& Q,
		                                    const lu_extract_options_type& options = lu_extract_options_type()) const {
			return this->policy_lu_with_info(
				static_cast<_P&>(L), static_cast<_P&>(U),
				static_cast<_P&>(P), static_cast<_P&>(Q), options);
		}

		// strict, permutation-vector form
		void lu(spmatrix& L, spmatrix& U, std::vector<index_type>& p, std::vector<index_type>& q,
		        const lu_extract_options_type& options = lu_extract_options_type()) const {
			const lu_extract_result_type result = lu_with_info(L, U, p, q, options);
			if (result.status != sparse_lu_extract_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::lu: factor extraction failed with status ",
					sparse_lu_extract_status_to_string(result.status));
			}
		}

		// strict, permutation-matrix form
		void lu(spmatrix& L, spmatrix& U, spmatrix& P, spmatrix& Q,
		        const lu_extract_options_type& options = lu_extract_options_type()) const {
			const lu_extract_result_type result = lu_with_info(L, U, P, Q, options);
			if (result.status != sparse_lu_extract_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::lu: factor extraction failed with status ",
					sparse_lu_extract_status_to_string(result.status));
			}
		}

		// ---------------------------------------------------------------
		// LU factor consumers (LUX-2) — thin forwarding wrappers, LDL
		// style: policy_* call + strict throw decision ONLY (the entire
		// implementation lives in the policy NVI pair; see
		// spmats_base/spmats_lu_extract_impl.hpp).  The factors are
		// ARGUMENTS in the SSC convention (any §C-conformant source works:
		// A.lu(...) or an external backend).  Vector permutation form
		// (p, q) only.  strict forms throw on ANY status != success.
		// ---------------------------------------------------------------

		// non-strict solve: x = Q U^{-1} L^{-1} P b
		lu_apply_result lu_solve_with_info(const spmatrix& L, const spmatrix& U,
		                                   const std::vector<index_type>& p, const std::vector<index_type>& q,
		                                   const std::vector<_T>& b, std::vector<_T>& x) const {
			return this->policy_lu_solve_with_info(
				static_cast<const _P&>(L), static_cast<const _P&>(U), p, q, b, x);
		}

		// strict solve
		std::vector<_T> lu_solve(const spmatrix& L, const spmatrix& U,
		                         const std::vector<index_type>& p, const std::vector<index_type>& q,
		                         const std::vector<_T>& b) const {
			std::vector<_T> x;
			const lu_apply_result result = this->policy_lu_solve_with_info(
				static_cast<const _P&>(L), static_cast<const _P&>(U), p, q, b, x);
			if (result.status != lu_apply_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::lu_solve: failed with status ",
					lu_apply_status_to_string(result.status));
			}
			return x;
		}

		// non-strict inverse row: row_i(A^{-1})
		lu_apply_result lu_inverse_row_with_info(const spmatrix& L, const spmatrix& U,
		                                         const std::vector<index_type>& p, const std::vector<index_type>& q,
		                                         const index_type i, std::vector<_T>& row) const {
			return this->policy_lu_inverse_row_with_info(
				static_cast<const _P&>(L), static_cast<const _P&>(U), p, q, i, row);
		}

		// strict inverse row
		std::vector<_T> lu_inverse_row(const spmatrix& L, const spmatrix& U,
		                               const std::vector<index_type>& p, const std::vector<index_type>& q,
		                               const index_type i) const {
			std::vector<_T> row;
			const lu_apply_result result = this->policy_lu_inverse_row_with_info(
				static_cast<const _P&>(L), static_cast<const _P&>(U), p, q, i, row);
			if (result.status != lu_apply_status::success) {
				vcp::throw_error<vcp::numerical_error>(
					"spmatrix::lu_inverse_row: failed with status ",
					lu_apply_status_to_string(result.status));
			}
			return row;
		}

		// ---------------------------------------------------------------
		// AINV approximate inverse (AINV-1) — thin forwarding wrappers
		// only, _with_info style exclusively (no strict throwing sugar,
		// design D-17).  Construction outputs Z / W unit upper triangular
		// and D diagonal, born-finalized; R = Z D^{-1} W^T is never
		// materialized (D-13) — use ainv_apply for its action on a vector
		// and ainv_residual_norm_estimate for a NON-GUARANTEED
		// ||I - R A||_inf value (D-11).  The entire implementation lives
		// in the policy layer (spmats_base/spmats_ainv_impl.hpp).
		// ---------------------------------------------------------------

		// non-strict construction (factor outputs; ldl_with_info(L,D,P)
		// forwarding form)
		ainv_result_type ainv_with_info(spmatrix& Z, spmatrix& W, spmatrix& D,
		                                const ainv_options_type& options = ainv_options_type()) const {
			return this->policy_ainv_with_info(
				static_cast<_P&>(Z), static_cast<_P&>(W), static_cast<_P&>(D), options);
		}

		// non-guaranteed residual norm estimate (const factors;
		// lu_solve_with_info forwarding form)
		ainv_status ainv_residual_norm_estimate(const spmatrix& Z, const spmatrix& W,
		                                        const spmatrix& D,
		                                        typename vcp::tsparse_scalar::real_type<_T>::type& est) const {
			return this->policy_ainv_residual_norm_estimate(
				static_cast<const _P&>(Z), static_cast<const _P&>(W),
				static_cast<const _P&>(D), est);
		}

		// apply: z = Z (D^{-1} (W^T r))
		ainv_status ainv_apply(const spmatrix& Z, const spmatrix& W, const spmatrix& D,
		                       const std::vector<_T>& r, std::vector<_T>& z) const {
			return this->policy_ainv_apply(
				static_cast<const _P&>(Z), static_cast<const _P&>(W),
				static_cast<const _P&>(D), r, z);
		}

		// ---------------------------------------------------------------
		// FSAI static factored approximate inverse (FSAI-1) — thin
		// forwarding wrappers only, _with_info style exclusively (F-D10).
		// Construction outputs the sqrt-free triple: U unit upper
		// triangular, D diagonal, perm new->old (P(perm[k],k) = 1),
		// born-finalized; R = P U D^{-1} U^T P^T is never materialized —
		// use fsai_apply for its action on a vector and
		// fsai_residual_norm_estimate for a NON-GUARANTEED
		// ||I - R A||_inf value.  The entire implementation lives in the
		// policy layer (spmats_base/spmats_fsai_impl.hpp).
		// ---------------------------------------------------------------

		// non-strict construction (factor outputs; ainv_with_info
		// forwarding form)
		fsai_result_type fsai_with_info(spmatrix& U, spmatrix& D,
		                                std::vector<index_type>& perm,
		                                const fsai_options_type& options = fsai_options_type()) const {
			return this->policy_fsai_with_info(
				static_cast<_P&>(U), static_cast<_P&>(D), perm, options);
		}

		// non-guaranteed residual norm estimate (const factors)
		fsai_status fsai_residual_norm_estimate(const spmatrix& U, const spmatrix& D,
		                                        const std::vector<index_type>& perm,
		                                        typename vcp::tsparse_scalar::real_type<_T>::type& est) const {
			return this->policy_fsai_residual_norm_estimate(
				static_cast<const _P&>(U), static_cast<const _P&>(D), perm, est);
		}

		// apply: z = P (U (D^{-1} (U^T (P^T r))))
		fsai_status fsai_apply(const spmatrix& U, const spmatrix& D,
		                       const std::vector<index_type>& perm,
		                       const std::vector<_T>& r, std::vector<_T>& z) const {
			return this->policy_fsai_apply(
				static_cast<const _P&>(U), static_cast<const _P&>(D), perm, r, z);
		}

		// ---------------------------------------------------------------
		// ADAPTIVE FSAI (FSAI-2) — thin forwarding wrappers only, same
		// output contract as fsai_with_info (sqrt-free triple U / D /
		// perm, born-finalized; R never materialized).  Estimate / apply
		// on the returned triple: use fsai_residual_norm_estimate /
		// fsai_apply above (no adaptive-specific helpers, F2-D6).  The
		// implementation lives in the policy layer
		// (spmats_base/spmats_fsai_adaptive_impl.hpp).
		// ---------------------------------------------------------------

		// diagonal-initialized adaptive construction
		fsai_adaptive_result_type fsai_adaptive_with_info(spmatrix& U, spmatrix& D,
		                                                  std::vector<index_type>& perm,
		                                                  const fsai_adaptive_options_type& options = fsai_adaptive_options_type()) const {
			return this->policy_fsai_adaptive_with_info(
				static_cast<_P&>(U), static_cast<_P&>(D), perm, options);
		}

		// initial-value adaptive construction (combined strategy, F2-D3):
		// U0 unit upper triangular in the same permuted frame as perm0;
		// options.ordering is ignored and perm0 is reused (design SS3.6)
		fsai_adaptive_result_type fsai_adaptive_with_info(spmatrix& U, spmatrix& D,
		                                                  std::vector<index_type>& perm,
		                                                  const spmatrix& U0,
		                                                  const std::vector<index_type>& perm0,
		                                                  const fsai_adaptive_options_type& options = fsai_adaptive_options_type()) const {
			return this->policy_fsai_adaptive_with_info(
				static_cast<_P&>(U), static_cast<_P&>(D), perm,
				static_cast<const _P&>(U0), perm0, options);
		}

		// Convenience overloads — build options and delegate to solve / solve_with_info
		std::vector<_T> solve_jacobi(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::jacobi;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = true;
			return this->policy_lss(b, opt_);
		}

		std::vector<_T> solve_gauss_seidel(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::gauss_seidel;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = true;
			return this->policy_lss(b, opt_);
		}

		std::vector<_T> solve_cg(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::conjugate_gradient;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.check_symmetric = true;
			opt_.preconditioner = preconditioner_type::none;
			opt_.use_relative_residual = true;
			return this->policy_lss(b, opt_);
		}

		std::vector<_T> solve_bicgstab(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::bicgstab;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = true;
			return this->policy_lss(b, opt_);
		}

		std::vector<_T> solve_gmres(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::gmres;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.restart = 30;
			opt_.use_relative_residual = true;
			return this->policy_lss(b, opt_);
		}

		linear_solve_result<_T> solve_jacobi_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                               const bool use_relative_residual = true) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::jacobi;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = use_relative_residual;
			return solve_with_info(b, opt_);
		}

		linear_solve_result<_T> solve_gauss_seidel_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                                     const bool use_relative_residual = true) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::gauss_seidel;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = use_relative_residual;
			return solve_with_info(b, opt_);
		}

		linear_solve_result<_T> solve_cg_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                           const bool check_symmetric = true,
		                                           const preconditioner_type preconditioner = preconditioner_type::none,
		                                           const bool use_relative_residual = true) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::conjugate_gradient;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.check_symmetric = check_symmetric;
			opt_.preconditioner = preconditioner;
			opt_.use_relative_residual = use_relative_residual;
			return solve_with_info(b, opt_);
		}

		linear_solve_result<_T> solve_bicgstab_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                                 const bool use_relative_residual = true,
		                                                 const preconditioner_type preconditioner = preconditioner_type::none) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::bicgstab;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = use_relative_residual;
			opt_.preconditioner = preconditioner;
			return solve_with_info(b, opt_);
		}

		linear_solve_result<_T> solve_gmres_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                             const std::size_t restart, const bool use_relative_residual = true,
		                                             const preconditioner_type preconditioner = preconditioner_type::none) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::gmres;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.restart = restart;
			opt_.use_relative_residual = use_relative_residual;
			opt_.preconditioner = preconditioner;
			return solve_with_info(b, opt_);
		}

		// ---------------------------------------------------------------
		// Eigenvalue API — all delegate to policy (no result checking here)
		// ---------------------------------------------------------------

		// Full dense eig (strict): policy_eig owns convergence checking
		eig_result<_T> eig(const eig_options_type& options) const {
			validate_eig_input("spmatrix::eig");
			return this->policy_eig(options);
		}

		// Full dense eig (non-strict / diagnostic)
		eig_result<_T> eig_with_info(const eig_options_type& options) const {
			validate_eig_input("spmatrix::eig");
			if (options.max_iter == 0 || options.tol <= scalar_real_type(0))
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eig: invalid iteration option");
			if (options.method != eig_solver_method::dense_fallback_explicit)
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eig: full dense eig requires dense_fallback_explicit");
			return this->policy_eigs_with_info(
				static_cast<std::size_t>(rowsize()), options);
		}

		eig_options_type default_eigs_options() const {
			eig_options_type options;
			options.structure = matrix_structure_hint::auto_detect;
			return options;
		}

		static eig_options_type dense_eigs_options() {
			eig_options_type options;
			options.method = eig_solver_method::dense_fallback_explicit;
			return options;
		}

		std::vector<_T> eigvals(const eig_options_type& options) const {
			return eig(options).eigenvalues;
		}

		// eigs(k): shorthand for eigs(k, default)
		std::vector<_T> eigs(const std::size_t k) const {
			return eigs(k, default_eigs_options());
		}

		// eigs(k, opt): strict — policy_eigs decides convergence
		std::vector<_T> eigs(const std::size_t k, const eig_options_type& options) const {
			return this->policy_eigs(k, options);
		}

		// eigs_with_info(k): shorthand
		eig_result<_T> eigs_with_info(const std::size_t k) const {
			return eigs_with_info(k, default_eigs_options());
		}

		// eigs_with_info(k, opt): non-strict, return full result
		eig_result<_T> eigs_with_info(const std::size_t k, const eig_options_type& options) const {
			return this->policy_eigs_with_info(
				k, options);
		}

		// Full generalized eig (strict, returns eigenvalue vector)
		std::vector<_T> eig(const spmatrix& B) const {
			validate_eig_input("spmatrix::eig(A,B)");
			return this->policy_generalized_eigs(
				static_cast<const _P&>(B),
				static_cast<std::size_t>(rowsize()), eig_options_type());
		}

		// Full generalized eig (strict, returns eig_result)
		eig_result<_T> eig(const spmatrix& B, const eig_options_type& options) const {
			validate_eig_input("spmatrix::eig(A,B)");
			return this->policy_generalized_eig(
				static_cast<const _P&>(B),
				static_cast<std::size_t>(rowsize()), options);
		}

		// Partial generalized eigs (strict): policy_generalized_eigs decides convergence
		std::vector<_T> eigs(const spmatrix& B, const std::size_t k, const eig_options_type& options = eig_options_type()) const {
			return this->policy_generalized_eigs(
				static_cast<const _P&>(B), k, options);
		}

		// Partial generalized eigs_with_info (non-strict)
		eig_result<_T> eigs_with_info(const spmatrix& B, const std::size_t k, const eig_options_type& options = eig_options_type()) const {
			return this->policy_generalized_eigs_with_info(
				static_cast<const _P&>(B), k, options);
		}

		// ------------------------------------------------------------------
		// Phase 6: preconditioner overloads — all delegate to policy
		// ------------------------------------------------------------------

		template <class Preconditioner>
		std::vector<_T> eigs(const std::size_t k, const eig_options_type& options,
		                     const Preconditioner& M) const {
			return this->policy_eigs(k, options, M);
		}

		template <class Preconditioner>
		eig_result<_T> eigs_with_info(const std::size_t k, const eig_options_type& options,
		                              const Preconditioner& M) const {
			return this->policy_eigs_with_info(
				k, options, M);
		}

		template <class Preconditioner>
		std::vector<_T> eigs(const spmatrix& B, const std::size_t k,
		                     const eig_options_type& options, const Preconditioner& M) const {
			return this->policy_generalized_eigs(
				static_cast<const _P&>(B), k, options, M);
		}

		template <class Preconditioner>
		eig_result<_T> eigs_with_info(const spmatrix& B, const std::size_t k,
		                              const eig_options_type& options,
		                              const Preconditioner& M) const {
			return this->policy_generalized_eigs_with_info(
				static_cast<const _P&>(B), k, options, M);
		}

		spmatrix transpose() const {
			spmatrix B;
			static_cast<vcp::spmats<_T, typename _P::index_type>&>(static_cast<_P&>(B)) = _P::transpose();
			// TPX-1: assign into the base spmats subobject so derived policies
			// (spumar, spcmodumar, ...) compile.  The result's policy-specific state is
			// default-constructed -- the same convention as every other
			// result-returning operation of spmatrix (T-2).
			return B;
		}

		// Matlab C = [A,B] -- thin forwarding to the policy's horzcat (CSR
		// direct-merge, see spmats::horzcat).
		void horzcat(const spmatrix& B, spmatrix& C) const {
			_P::horzcat(static_cast<const _P&>(B), static_cast<_P&>(C));
		}

		// Matlab C = [A;B] -- thin forwarding to the policy's vercat.
		void vercat(const spmatrix& B, spmatrix& C) const {
			_P::vercat(static_cast<const _P&>(B), static_cast<_P&>(C));
		}

		void assign_csr(const index_type rows, const index_type cols,
		                const std::vector<index_type>& row_ptr,
		                const std::vector<index_type>& col_ind,
		                const std::vector<_T>& value) {
			_P::assign_csr(rows, cols, row_ptr, col_ind, value);
		}

		void assign_csc(const index_type rows, const index_type cols,
		                const std::vector<index_type>& col_ptr,
		                const std::vector<index_type>& row_ind,
		                const std::vector<_T>& value) {
			_P::assign_csc(rows, cols, col_ptr, row_ind, value);
		}

		const std::vector<index_type>& outer_index() const { return _P::outer_index(); }
		const std::vector<index_type>& inner_index() const { return _P::inner_index(); }
		const std::vector<_T>& values() const { return _P::values(); }
		const std::vector<index_type>& coo_rows() const { return _P::coo_rows(); }
		const std::vector<index_type>& coo_columns() const { return _P::coo_columns(); }
		const std::vector<_T>& coo_values() const { return _P::coo_values(); }

		// A * x — delegates to policy_mul_vec
		friend std::vector<_T> operator*(const spmatrix& A, const std::vector<_T>& x) {
			return A.policy_mul_vec(static_cast<const _P&>(A), x);
		}

		// x * A — delegates to policy_left_mul_vec
		friend std::vector<_T> operator*(const std::vector<_T>& x, const spmatrix& A) {
			return A.policy_left_mul_vec(x, static_cast<const _P&>(A));
		}

		friend spmatrix transpose(const spmatrix& A) {
			return A.transpose();
		}

		//***************** Math functions (SPFN; mirrors matrix.hpp 627-686) *****************//
		// abs/sqrt/sin: sparsity-preserving.  cos/exp/log: MATLAB sparse
		// semantics -- applied at every position, result structurally
		// (near-)full (see the SPFN block in spmats.hpp).  const version
		// copies then applies; rvalue version applies in place.
		friend spmatrix abs(const spmatrix& A) {
			spmatrix C;
			C = A;
			C.abs();
			return C;
		}
		friend spmatrix abs(spmatrix&& A) {
			A.abs();
			return std::move(A);
		}
		friend spmatrix sqrt(const spmatrix& A) {
			spmatrix C;
			C = A;
			C.sqrt();
			return C;
		}
		friend spmatrix sqrt(spmatrix&& A) {
			A.sqrt();
			return std::move(A);
		}
		friend spmatrix sin(const spmatrix& A) {
			spmatrix C;
			C = A;
			C.sin();
			return C;
		}
		friend spmatrix sin(spmatrix&& A) {
			A.sin();
			return std::move(A);
		}
		friend spmatrix cos(const spmatrix& A) {
			spmatrix C;
			C = A;
			C.cos();
			return C;
		}
		friend spmatrix cos(spmatrix&& A) {
			A.cos();
			return std::move(A);
		}
		friend spmatrix exp(const spmatrix& A) {
			spmatrix C;
			C = A;
			C.exp();
			return C;
		}
		friend spmatrix exp(spmatrix&& A) {
			A.exp();
			return std::move(A);
		}
		friend spmatrix log(const spmatrix& A) {
			spmatrix C;
			C = A;
			C.log();
			return C;
		}
		friend spmatrix log(spmatrix&& A) {
			A.log();
			return std::move(A);
		}

		//************* matlab like reductions (SPFN; mirrors matrix.hpp 704-732) *************//
		// max/min: 1x1 for scalars/vectors, 1 x column (column-wise) for
		// matrices; implicit zeros participate (判断B).  normone/norminf/
		// normtwo: 1x1 spmatrix (operator norms for matrices, as in
		// mats<T>); normtwo's matrix branch keeps the dynamic policy type
		// (判断C, see spmats::normtwo).
		friend spmatrix max(const spmatrix& A) {
			spmatrix c;
			A.max(c);
			return c;
		}
		friend spmatrix min(const spmatrix& A) {
			spmatrix c;
			A.min(c);
			return c;
		}
		friend spmatrix normone(const spmatrix& A) {
			spmatrix c;
			A.normone(c);
			return c;
		}
		friend spmatrix normtwo(const spmatrix& A) {
			spmatrix c = A;
			c.normtwo();
			return c;
		}
		friend spmatrix normtwo(spmatrix&& A) {
			A.normtwo();
			return std::move(A);
		}
		friend spmatrix norminf(const spmatrix& A) {
			spmatrix c;
			A.norminf(c);
			return c;
		}

	private:
		void validate_eig_input(const char* routine) const {
			if (rowsize() != columnsize()) vcp::throw_error<vcp::dimension_error>(routine, ": matrix must be square");
		}
	};

	// -----------------------------------------------------------------------
	// SUB-1: block-write proxy returned by the non-const
	// spmatrix::operator()({..},{..}) and the mixed int/list forms.
	// Holds a reference to the parent and the NORMALIZED selectors only;
	// selector validation is completed at proxy construction
	// (vcp::index_error), and operator= validates the right-hand-side
	// dimensions (and an index-width preflight) BEFORE the first set(), so
	// a throwing assignment leaves the parent untouched.
	// Semantics (SUB-1_design.md v1.1 §3.3, MATLAB block replacement):
	//   proxy = spmatrix<_T,_P> -- the block is REPLACED: positions where
	//                              the right-hand side is zero lose any
	//                              previously stored entry.  The right-hand
	//                              side is materialized first (as_csr()
	//                              copy), so aliased assignments such as
	//                              A({0,2},{0,2}) = A({1,3},{1,3}) follow
	//                              the MATLAB semantics; B is unchanged.
	//   proxy = scalar _T       -- fill: a nonzero scalar DENSIFIES the
	//                              block (nnz grows by the block area, the
	//                              MATLAB consequence); scalar 0 clears it.
	//   spmatrix<_T,_P>(proxy)  -- read; delegates to spmatrix::submatrix.
	// Mechanism: public tagged set()/to_csr()/to_csc() path only (C-3) --
	// set() pushes tagged COO entries, the normalize inside the final
	// to_csr()/to_csc() applies the replace-then-drop-zeros semantics, and
	// the parent is re-finalized into its pre-write format (finalized CSC
	// stays CSC; CSR and unfinalized parents finalize to CSR).  Complexity
	// O(nnz + block area + sort).
	// LIFETIME WARNING: `auto x = A({..},{..});` captures the PROXY, not an
	// spmatrix -- it must not outlive the parent A.  Use an explicit
	// spmatrix<_T,_P> variable to take a copy of the block.
	// -----------------------------------------------------------------------
	template <typename _T, class _P> class spmatrix_block {
	public:
		typedef typename spmatrix< _T, _P >::index_type index_type;

		spmatrix_block(spmatrix< _T, _P >& A,
		               const index_type r0, const index_type rs, const index_type rn, const bool rfull,
		               const index_type c0, const index_type cs, const index_type cn, const bool cfull)
			: A_(A), r0_(r0), rs_(rs), rn_(rn), rfull_(rfull),
			  c0_(c0), cs_(cs), cn_(cn), cfull_(cfull) {}

		spmatrix_block(const spmatrix_block&) = default;

		// block replacement (right-hand side materialized first, S6)
		spmatrix_block& operator=(const spmatrix< _T, _P >& B) {
			if (B.rowsize() != rn_ || B.columnsize() != cn_) {
				vcp::throw_error<vcp::dimension_error>(
					"spmatrix_block: block assignment size mismatch: ",
					B.rowsize(), "x", B.columnsize(), " != ", rn_, "x", cn_);
			}
			preflight_();
			// materialize the right-hand side FIRST: an as_csr() copy is
			// aliasing-safe (B may be the parent or overlap the block),
			// leaves B untouched, and gives O(log) finalized reads.
			const spmatrix< _T, _P > Bc = B.as_csr();
			const bool was_csc = A_.is_finalized() && A_.format() == vcp::sparse_csc;
			for (index_type i = 0; i < rn_; i++) {
				for (index_type j = 0; j < cn_; j++) {
					// tagged set() on EVERY block position (zeros included):
					// the tagged merge discards older entries at the
					// position and the zero-drop pass removes the zeros, so
					// "zeros in B delete old entries" holds automatically.
					A_.set(r0_ + i * rs_, c0_ + j * cs_, Bc.get(i, j));
				}
			}
			if (was_csc) {
				A_.to_csc();
			}
			else {
				A_.to_csr();
			}
			return *this;
		}

		// scalar fill (nonzero: densifies the block; zero: clears it)
		spmatrix_block& operator=(const _T& s) {
			preflight_();
			const bool was_csc = A_.is_finalized() && A_.format() == vcp::sparse_csc;
			for (index_type i = 0; i < rn_; i++) {
				for (index_type j = 0; j < cn_; j++) {
					A_.set(r0_ + i * rs_, c0_ + j * cs_, s);
				}
			}
			if (was_csc) {
				A_.to_csc();
			}
			else {
				A_.to_csr();
			}
			return *this;
		}

		// proxy = proxy (e.g. A({0,2},{0,2}) = A({1,3},{1,3})): materialize
		// the right-hand block first, then block-assign.  Without this
		// overload the implicitly-deleted copy assignment would win the
		// overload resolution (same reason as spmats_element_proxy).
		spmatrix_block& operator=(const spmatrix_block& other) {
			return (*this) = static_cast<spmatrix< _T, _P > >(other);
		}

		// read conversion -- delegates to spmatrix::submatrix, rebuilding
		// the normalized selectors as {} / {start, stride, last} (identical
		// selections by construction; the values came from int selectors).
		operator spmatrix< _T, _P >() const {
			const int r0 = static_cast<int>(r0_), rs = static_cast<int>(rs_);
			const int rl = static_cast<int>(r0_ + (rn_ - 1) * rs_);
			const int c0 = static_cast<int>(c0_), cs = static_cast<int>(cs_);
			const int cl = static_cast<int>(c0_ + (cn_ - 1) * cs_);
			if (rfull_ && cfull_) {
				return A_.submatrix({}, {});
			}
			if (rfull_) {
				return A_.submatrix({}, {c0, cs, cl});
			}
			if (cfull_) {
				return A_.submatrix({r0, rs, rl}, {});
			}
			return A_.submatrix({r0, rs, rl}, {c0, cs, cl});
		}

	private:
		// no-throw-after-start guard: once validation is done, the only
		// failure mode inside the set()/normalize/finalize sequence (apart
		// from std::bad_alloc, which no path can exclude) would be the
		// size_to_index overflow check firing when stored entries + block
		// area exceed the index_type range; reject that case here, before
		// the first set().
		void preflight_() const {
			const std::size_t limit = static_cast<std::size_t>((std::numeric_limits<index_type>::max)());
			const std::size_t area = static_cast<std::size_t>(rn_) * static_cast<std::size_t>(cn_);
			if (area > limit - static_cast<std::size_t>(A_.stored_nnz())) {
				vcp::throw_error<vcp::invalid_argument>(
					"spmatrix_block: stored entries + block area exceed index_type range");
			}
		}

		spmatrix< _T, _P >& A_;
		index_type r0_, rs_, rn_;
		bool rfull_;
		index_type c0_, cs_, cn_;
		bool cfull_;
	};

	// SUB-1: out-of-line definitions of the proxy-returning operator()
	// overloads (declared inside spmatrix; spmatrix_block must be complete
	// here).  Validation and normalization use the same protected policy
	// helpers as spmats::submat (single validation implementation on the
	// sparse side); the mixed forms wrap the int as the single selector {i}.
	template <typename _T, class _P>
	spmatrix_block< _T, _P > spmatrix< _T, _P >::operator()(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) {
		_P::submat_check_sizes_(list1, list2);
		const typename _P::submat_axis_ ra = _P::submat_normalize_(list1, this->rowsize(), true);
		const typename _P::submat_axis_ ca = _P::submat_normalize_(list2, this->columnsize(), false);
		return spmatrix_block< _T, _P >(*this,
			ra.start, ra.stride, ra.count, ra.full,
			ca.start, ca.stride, ca.count, ca.full);
	}
	template <typename _T, class _P>
	spmatrix_block< _T, _P > spmatrix< _T, _P >::operator()(const int i, const std::initializer_list<int>& list2) {
		return (*this)({i}, list2);
	}
	template <typename _T, class _P>
	spmatrix_block< _T, _P > spmatrix< _T, _P >::operator()(const std::initializer_list<int>& list1, const int j) {
		return (*this)(list1, {j});
	}


	// -----------------------------------------------------------------------
	// MATLAB-like free functions: lss / lss_with_info
	// MATLAB: x = A \ b  →  C++: x = lss(A, b)
	// These delegate to A.solve() / A.solve_with_info(), which delegate to
	// the policy methods policy_lss / policy_lss_with_info.
	// Verified policy can override those methods for enclosure semantics.
	// -----------------------------------------------------------------------

	template <typename _T, class _P>
	std::vector<_T> lss(const spmatrix<_T, _P>& A, const std::vector<_T>& b) {
		return A.solve(b);
	}

	template <typename _T, class _P>
	std::vector<_T> lss(const spmatrix<_T, _P>& A, const std::vector<_T>& b,
	                    const linear_solve_options<_T>& opt) {
		return A.solve(b, opt);
	}

	template <typename _T, class _P>
	linear_solve_result<_T> lss_with_info(const spmatrix<_T, _P>& A,
	                                      const std::vector<_T>& b) {
		return A.solve_with_info(b);
	}

	template <typename _T, class _P>
	linear_solve_result<_T> lss_with_info(const spmatrix<_T, _P>& A,
	                                      const std::vector<_T>& b,
	                                      const linear_solve_options<_T>& opt) {
		return A.solve_with_info(b, opt);
	}

	// -----------------------------------------------------------------------
	// MATLAB-like free functions: horzcat / vercat (variadic)
	// MATLAB: C = [A,B]  ->  C++: C = horzcat(A, B)
	// MATLAB: C = [A;B]  ->  C++: C = vercat(A, B)
	// Ordinary namespace-scope templates (not friends, matching the lss/
	// lss_with_info convention just above) so the recursive unqualified
	// call in the variadic overload resolves to this overload set rather
	// than being hidden by spmatrix::horzcat/vercat's own 2-arg member
	// (a friend defined only inside the class is found solely via ADL,
	// which the class's own member of the same name would shadow).
	// -----------------------------------------------------------------------

	template <typename _T, class _P>
	spmatrix<_T, _P> horzcat(const spmatrix<_T, _P>& A) {
		return A;
	}

	template <typename _T, class _P>
	spmatrix<_T, _P> horzcat(const spmatrix<_T, _P>& A, const spmatrix<_T, _P>& B) {
		spmatrix<_T, _P> C;
		A.horzcat(B, C);
		return C;
	}

	template <typename _T, class _P, typename... Args>
	spmatrix<_T, _P> horzcat(const spmatrix<_T, _P>& A, const spmatrix<_T, _P>& B, const Args&... args) {
		spmatrix<_T, _P> C;
		A.horzcat(B, C);
		return horzcat(C, args...);
	}

	template <typename _T, class _P>
	spmatrix<_T, _P> vercat(const spmatrix<_T, _P>& A) {
		return A;
	}

	template <typename _T, class _P>
	spmatrix<_T, _P> vercat(const spmatrix<_T, _P>& A, const spmatrix<_T, _P>& B) {
		spmatrix<_T, _P> C;
		A.vercat(B, C);
		return C;
	}

	template <typename _T, class _P, typename... Args>
	spmatrix<_T, _P> vercat(const spmatrix<_T, _P>& A, const spmatrix<_T, _P>& B, const Args&... args) {
		spmatrix<_T, _P> C;
		A.vercat(B, C);
		return vercat(C, args...);
	}

	template <typename _T, typename _Index = int> class spmatrix_builder {
	public:
		typedef _T value_type;
		typedef _Index index_type;
		typedef vcp::spmatrix<value_type, vcp::spmats<value_type, index_type> > matrix_type;

		spmatrix_builder(const index_type rows, const index_type cols)
			: matrix_(rows, cols) {}

		void reserve(const std::size_t n) {
			if (n > static_cast<std::size_t>(std::numeric_limits<index_type>::max())) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix_builder::reserve: size exceeds index range");
			}
			matrix_.reserve(static_cast<index_type>(n));
		}

		void add(const index_type i, const index_type j, const value_type& value) {
			matrix_.add(i, j, value);
		}

		void set(const index_type i, const index_type j, const value_type& value) {
			matrix_.set(i, j, value);
		}

		matrix_type finalize_csr() const {
			matrix_type result = matrix_;
			result.finalize();
			return result;
		}

		matrix_type finalize_csc() const {
			matrix_type result = matrix_;
			result.to_csc();
			return result;
		}

	private:
		matrix_type matrix_;
	};
}

#endif
