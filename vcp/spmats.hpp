// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATS_HPP
#define VCP_SPMATS_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse.hpp>
#include <vcp/spmats_base/spmats_eigs_types.hpp>
#include <vcp/spmats_base/spmats_ldl.hpp>
#include <vcp/spmats_base/spmats_chol.hpp>
#include <vcp/spmats_base/spmats_lu_extract.hpp>
#include <vcp/spmats_base/spmats_lu_factor.hpp>
#include <vcp/spmats_base/spmats_ldl_shift.hpp>
#include <vcp/spmats_base/spmats_ainv.hpp>
#include <vcp/spmats_base/spmats_fsai.hpp>
#include <vcp/spmats_base/spmats_fsai_adaptive.hpp>
#include <vcp/spmats_base/spmats_policy_traits.hpp>

namespace vcp {

	template <typename _T, typename _Index = int> class spmats {
	public:
		typedef _Index index_type;
		typedef _T value_type;
		typedef vcp::sparse_format format_type;

		spmats()
			: row(0), column(0), fmt(vcp::sparse_coo), finalized(false), sorted(true), unique(true) {}

		spmats(const spmats&) = default;
		spmats(spmats&&) = default;
		spmats& operator=(const spmats&) = default;
		spmats& operator=(spmats&&) = default;
		virtual ~spmats() = default;

		index_type rowsize() const { return row; }
		index_type columnsize() const { return column; }
		index_type stored_nnz() const {
			return finalized ? size_to_index(value.size(), "spmats::stored_nnz") : size_to_index(coo_value.size(), "spmats::stored_nnz");
		}
		index_type nnz() const {
			if (finalized) return size_to_index(value.size(), "spmats::nnz");
			index_type count = 0;
			for (std::size_t k = 0; k < coo_value.size(); k++) {
				if (!(coo_value[k] == _T(0))) count++;
			}
			return count;
		}
		bool is_finalized() const { return finalized; }
		bool is_sorted() const { return sorted; }
		bool is_unique() const { return unique; }
		format_type format() const { return fmt; }

		void resize(const index_type rows, const index_type cols) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::resize: negative size");
			}
			row = rows;
			column = cols;
			clear_storage();
			fmt = vcp::sparse_coo;
			finalized = false;
			sorted = true;
			unique = true;
		}

		void clear() {
			row = 0;
			column = 0;
			clear_storage();
			fmt = vcp::sparse_coo;
			finalized = false;
			sorted = true;
			unique = true;
		}

		void reserve(const index_type n) {
			if (n < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::reserve: negative size");
			}
			const std::size_t m = index_to_size(n, "spmats::reserve");
			coo_row.reserve(m);
			coo_col.reserve(m);
			coo_value.reserve(m);
			coo_is_set.reserve(m);
			inner.reserve(m);
			value.reserve(m);
		}

		void add(const index_type i, const index_type j, const _T& a) {
			check_index(i, j, "spmats::add");
			ensure_coo_buffer();
			coo_row.push_back(i);
			coo_col.push_back(j);
			coo_value.push_back(a);
			coo_is_set.push_back(0);
			finalized = false;
			fmt = vcp::sparse_coo;
			sorted = false;
			unique = false;
		}

		// Option A (see sandbox/docs/design/spmats_finalize_policy.md and
		// SLU-C3-ELEMENT-ACCESSOR): set() is just as light as add() (O(1)
		// amortized push, no COO rescan). It tags the pushed entry as
		// is_set=true; the actual "discard everything before this at (i,j)"
		// semantics are applied later by the tagged merge in normalize_coo()
		// (or, before finalize, by get()'s own tagged scan below).
		void set(const index_type i, const index_type j, const _T& a) {
			check_index(i, j, "spmats::set");
			ensure_coo_buffer();
			coo_row.push_back(i);
			coo_col.push_back(j);
			coo_value.push_back(a);
			coo_is_set.push_back(1);
			finalized = false;
			fmt = vcp::sparse_coo;
			sorted = false;
			unique = false;
		}

		_T get(const index_type i, const index_type j) const {
			check_index(i, j, "spmats::get");
			if (!finalized) {
				// Scan in stored order (insertion order, unless sort_coo()/
				// normalize_coo() reordered it -- either way entries sharing
				// (i,j) keep their relative insertion order, see the tagged
				// tcoo_sort overload). A set() resets the accumulator to its
				// value; an add() accumulates. Byte-identical to the old
				// unconditional sum when no set() ever touched (i,j).
				_T sum = _T(0);
				for (std::size_t k = 0; k < coo_value.size(); k++) {
					if (coo_row[k] == i && coo_col[k] == j) {
						if (coo_is_set[k]) sum = coo_value[k];
						else sum += coo_value[k];
					}
				}
				return sum;
			}
			if (fmt == vcp::sparse_csr) {
				const index_type first = outer[i];
				const index_type last = outer[i + 1];
				const typename std::vector<index_type>::const_iterator begin = inner.begin() + first;
				const typename std::vector<index_type>::const_iterator end = inner.begin() + last;
				typename std::vector<index_type>::const_iterator it = std::lower_bound(begin, end, j);
				if (it != end && *it == j) {
					return value[static_cast<std::size_t>(it - inner.begin())];
				}
				return _T(0);
			}
			if (fmt == vcp::sparse_csc) {
				const index_type first = outer[j];
				const index_type last = outer[j + 1];
				const typename std::vector<index_type>::const_iterator begin = inner.begin() + first;
				const typename std::vector<index_type>::const_iterator end = inner.begin() + last;
				typename std::vector<index_type>::const_iterator it = std::lower_bound(begin, end, i);
				if (it != end && *it == i) {
					return value[static_cast<std::size_t>(it - inner.begin())];
				}
				return _T(0);
			}
			return _T(0);
		}

		void finalize() const { to_csr(); }

		// Calls ensure_coo_buffer() (see the warning above its declaration):
		// on a finalized (CSR/CSC) matrix this silently reverts it to
		// unfinalized COO, and sort_coo() has no path back to finalized
		// state. Safe only from to_csr()/to_csc() (which always re-finalize
		// afterward) or on a matrix that is already unfinalized. Do not call
		// directly on a finalized matrix in new code -- use as_csr()/
		// as_csc() instead (see
		// sandbox/docs/misc/spmats_destructive_helpers_misc.md).
		void sort_coo() const {
			ensure_coo_buffer();
			const index_type n = size_to_index(coo_value.size(), "spmats::sort_coo");
			if (n > 0) {
				vcp::tcoo_sort(n, coo_row.data(), coo_col.data(), coo_value.data(), coo_is_set.data());
			}
			sorted = true;
		}

		// Calls ensure_coo_buffer()/sort_coo() (same caveat as sort_coo()
		// above): reverts a finalized matrix to unfinalized COO with no way
		// back. Only safe from to_csr()/to_csc() or on an already-
		// unfinalized matrix; do not call directly on a finalized matrix in
		// new code -- use as_csr()/as_csc() instead.
		void normalize_coo() const {
			ensure_coo_buffer();
			sort_coo();
			index_type n = size_to_index(coo_value.size(), "spmats::normalize_coo");
			if (n > 0) {
				n = vcp::tcoo_sum_duplicates_tagged(n, coo_row.data(), coo_col.data(), coo_value.data(), coo_is_set.data());
				n = vcp::tcoo_remove_zeros(n, coo_row.data(), coo_col.data(), coo_value.data());
			}
			coo_row.resize(static_cast<std::size_t>(n));
			coo_col.resize(static_cast<std::size_t>(n));
			coo_value.resize(static_cast<std::size_t>(n));
			// Each (row,col) is now unique, so the tag no longer matters
			// (see get()'s reset-from-zero equivalence); keep the array in
			// lockstep with the other three so indices always line up.
			coo_is_set.resize(static_cast<std::size_t>(n));
			sorted = true;
			unique = true;
		}

		void to_csr() const {
			if (finalized && fmt == vcp::sparse_csr) return;
			ensure_coo_buffer();
			normalize_coo();
			const index_type n = size_to_index(coo_value.size(), "spmats::to_csr");
			outer.assign(index_to_size(checked_plus_one(row, "spmats::to_csr"), "spmats::to_csr"), 0);
			inner.assign(index_to_size(n, "spmats::to_csr"), 0);
			value.assign(index_to_size(n, "spmats::to_csr"), _T(0));
			if (n > 0 || row >= 0) {
				vcp::tcoo_to_csr(row, column, n, coo_row.data(), coo_col.data(), coo_value.data(), outer.data(), inner.data(), value.data());
			}
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			coo_is_set.clear();
			fmt = vcp::sparse_csr;
			finalized = true;
			sorted = true;
			unique = true;
		}

		void to_csc() const {
			if (finalized && fmt == vcp::sparse_csc) return;
			ensure_coo_buffer();
			normalize_coo();
			const index_type n = size_to_index(coo_value.size(), "spmats::to_csc");
			outer.assign(index_to_size(checked_plus_one(column, "spmats::to_csc"), "spmats::to_csc"), 0);
			inner.assign(index_to_size(n, "spmats::to_csc"), 0);
			value.assign(index_to_size(n, "spmats::to_csc"), _T(0));
			if (n > 0 || column >= 0) {
				vcp::tcoo_to_csc(row, column, n, coo_row.data(), coo_col.data(), coo_value.data(), outer.data(), inner.data(), value.data());
			}
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			coo_is_set.clear();
			fmt = vcp::sparse_csc;
			finalized = true;
			sorted = true;
			unique = true;
		}

		spmats as_csr() const {
			spmats tmp(*this);
			tmp.to_csr();
			return tmp;
		}

		spmats as_csc() const {
			spmats tmp(*this);
			tmp.to_csc();
			return tmp;
		}

		std::vector<_T> mul_vec(const std::vector<_T>& x) const {
			if (column != size_to_index(x.size(), "spmats::mul_vec")) {
				vcp::throw_error<vcp::dimension_error>("spmats::mul_vec: dimension mismatch");
			}
			std::vector<_T> y(index_to_size(row, "spmats::mul_vec"), _T(0));
			mul_vec(x.data(), y.data());
			return y;
		}

		void mul_vec(const _T* x, _T* y) const {
			require_finalized("spmats::mul_vec");
			if (fmt == vcp::sparse_csr) {
				vcp::tcsrmv('N', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else if (fmt == vcp::sparse_csc) {
				vcp::tcscmv('N', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else {
				vcp::throw_error<vcp::state_error>("spmats::mul_vec: unsupported sparse format");
			}
		}

		std::vector<_T> trans_mul_vec(const std::vector<_T>& x) const {
			if (row != size_to_index(x.size(), "spmats::trans_mul_vec")) {
				vcp::throw_error<vcp::dimension_error>("spmats::trans_mul_vec: dimension mismatch");
			}
			std::vector<_T> y(index_to_size(column, "spmats::trans_mul_vec"), _T(0));
			trans_mul_vec(x.data(), y.data());
			return y;
		}

		void trans_mul_vec(const _T* x, _T* y) const {
			require_finalized("spmats::trans_mul_vec");
			if (fmt == vcp::sparse_csr) {
				vcp::tcsrmv('T', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else if (fmt == vcp::sparse_csc) {
				vcp::tcscmv('T', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else {
				vcp::throw_error<vcp::state_error>("spmats::trans_mul_vec: unsupported sparse format");
			}
		}

		spmats transpose() const {
			spmats tmp(*this);
			tmp.ensure_coo_buffer();
			spmats out;
			out.resize(column, row);
			out.reserve(size_to_index(tmp.coo_value.size(), "spmats::transpose"));
			for (std::size_t k = 0; k < tmp.coo_value.size(); k++) {
				out.coo_row.push_back(tmp.coo_col[k]);
				out.coo_col.push_back(tmp.coo_row[k]);
				out.coo_value.push_back(tmp.coo_value[k]);
				out.coo_is_set.push_back(tmp.coo_is_set[k]);
			}
			out.finalize();
			return out;
		}

		// ------------------------------------------------------------------
		// SPFN: elementwise math functions / min / max / norms
		// (pure addition; mirrors mats<T> semantics in mats.hpp 600-670 /
		//  1068-1251 and the shape contracts of matrix.hpp 627-731.  MATLAB
		//  sparse semantics for the elementwise maps: f with f(0)=0 (abs,
		//  sqrt, sin) preserves the sparsity pattern; f with f(0)!=0 (cos,
		//  exp, log) is applied at EVERY position, so the result is
		//  structurally (near-)full (MathWorks spfun doc: "exp(S) returns 1
		//  for the elements of S that are 0s").  Requirements on _T are the
		//  spmats<T> baseline: arithmetic, comparison, abs and the used
		//  elementary function found by ADL (using std::xxx).
		// ------------------------------------------------------------------
	private:
		// f(0)=0 maps: finalize first (COO duplicates are summed by the
		// finalize/normalize pass; mapping before summation would be wrong
		// since f(a+b) != f(a)+f(b)), then map the stored values in place.
		// If f maps a stored nonzero to exactly _T(0), the entry is
		// compacted away to keep the no-explicit-zero invariant of
		// finalized storage (same invariant to_csr enforces via
		// tcoo_remove_zeros; the zero test `== _T(0)` is the established
		// idiom of nnz()/tcoo_remove_zeros/validate_csr).  A finalized CSC
		// matrix keeps its CSC format (value-array map is format-agnostic).
		template <class _Fn> void spfn_map_stored_(_Fn f) {
			if (!finalized) to_csr();
			bool has_zero = false;
			for (std::size_t k = 0; k < value.size(); k++) {
				value[k] = f(value[k]);
				if (value[k] == _T(0)) has_zero = true;
			}
			if (!has_zero) return;
			// compact within the current (CSR or CSC) format
			const std::size_t nouter = (fmt == vcp::sparse_csr)
				? index_to_size(row, "spmats::spfn_map_stored_")
				: index_to_size(column, "spmats::spfn_map_stored_");
			std::size_t w = 0;
			index_type base = 0;
			for (std::size_t o = 0; o < nouter; o++) {
				const index_type first = outer[o];
				const index_type last = outer[o + 1];
				outer[o] = base;
				for (index_type p = first; p < last; p++) {
					const std::size_t sp = static_cast<std::size_t>(p);
					if (!(value[sp] == _T(0))) {
						inner[w] = inner[sp];
						value[w] = value[sp];
						w++;
						base++;
					}
				}
			}
			outer[nouter] = base;
			inner.resize(w);
			value.resize(w);
		}
		struct spfn_abs_ { _T operator()(const _T& x) const { using std::abs; return abs(x); } };
		// 判断E: no complex _T exists, so a certifiably negative component
		// is an error (double: explicit throw below; kv::interval: kv's own
		// throw also remains possible for cases the certified `<` cannot
		// decide, e.g. an interval straddling zero -- both are acceptable).
		// log(_T(0)) is not negative, so it passes this guard and inherits
		// the type's own behavior (-inf for double, no exception).
		struct spfn_sqrt_ {
			_T operator()(const _T& x) const {
				if (x < _T(0)) vcp::throw_error<vcp::invalid_argument>("spmats::sqrt: negative component (no complex support)");
				using std::sqrt; return sqrt(x);
			}
		};
		struct spfn_sin_ { _T operator()(const _T& x) const { using std::sin; return sin(x); } };
		// f(0)!=0 maps (SPFN 判断A): build the structurally (near-)full CSR
		// result.  INVARIANT (measured): finalized storage never contains
		// explicit zeros (to_csr drops them via normalize_coo; assign_csr
		// rejects them).  Therefore positions where f evaluates to exactly
		// _T(0) (e.g. log(1)) are NOT stored -- they become implicit zeros,
		// which represents the same mathematical matrix and matches MATLAB
		// (sparse storage never keeps explicit zeros).  f(_T(0)) is
		// evaluated once and reused for every implicit-zero position.
		template <class _Fn> void spfn_map_full_(_Fn f, const char* routine) {
			if (row <= 0 || column <= 0) {
				vcp::throw_error<vcp::state_error>(routine, ": empty matrix");
			}
			to_csr();
			const std::size_t m  = index_to_size(row, routine);
			const std::size_t nc = index_to_size(column, routine);
			const _T f0 = f(_T(0));
			// dense row buffer of mapped values, flushed row by row
			std::vector<index_type> new_outer(m + 1);
			std::vector<index_type> new_inner;
			std::vector<_T> new_value;
			new_inner.reserve(m * nc);
			new_value.reserve(m * nc);
			std::vector<_T> rowbuf(nc);
			new_outer[0] = 0;
			for (std::size_t i = 0; i < m; i++) {
				for (std::size_t j = 0; j < nc; j++) rowbuf[j] = f0;
				for (index_type p = outer[i]; p < outer[i + 1]; p++) {
					const std::size_t j = static_cast<std::size_t>(inner[static_cast<std::size_t>(p)]);
					rowbuf[j] = f(value[static_cast<std::size_t>(p)]);
				}
				for (std::size_t j = 0; j < nc; j++) {
					if (!(rowbuf[j] == _T(0))) {   // drop exact zeros (invariant)
						new_inner.push_back(static_cast<index_type>(j));
						new_value.push_back(rowbuf[j]);
					}
				}
				new_outer[i + 1] = size_to_index(new_value.size(), routine);
			}
			assign_csr(row, column, new_outer, new_inner, new_value);
		}
		struct spfn_cos_ { _T operator()(const _T& x) const { using std::cos; return cos(x); } };
		struct spfn_exp_ { _T operator()(const _T& x) const { using std::exp; return exp(x); } };
		// 判断E guard as in spfn_sqrt_; log(_T(0)) is NOT negative, so
		// implicit zeros inherit the type's own log(0) behavior (SPFN 判断A
		// [実測]: double/kv::dd -> -inf value, no exception; the -inf is a
		// storable nonzero).
		struct spfn_log_ {
			_T operator()(const _T& x) const {
				if (x < _T(0)) vcp::throw_error<vcp::invalid_argument>("spmats::log: negative component (no complex support)");
				using std::log; return log(x);
			}
		};
		// SPFN 判断F: NaN must PROPAGATE through max/min folds, never be
		// silently dropped (MATLAB's omitnan default is deliberately NOT
		// followed -- silent NaN omission is a bug incubator for verified
		// computation).  Detection is x != x (true only for NaN-like
		// values; for kv::interval this is always false [実測], so the
		// guard is a no-op there).  std::max/std::min alone would be
		// order-dependent for NaN, hence the explicit guard.  The 2-ary
		// max/min and the comparison are resolved by ADL (using std::xxx),
		// the spmats<T> baseline requirement -- NOT a total order (kv's
		// endpointwise interval max/min is a verified enclosure of the true
		// max and is in scope; complex T has no order and is out of scope,
		// as in dense mats).
		struct spfn_takemax_ {
			_T operator()(const _T& a, const _T& b) const {
				if (a != a) return a;
				if (b != b) return b;
				using std::max; return max(a, b);
			}
		};
		struct spfn_takemin_ {
			_T operator()(const _T& a, const _T& b) const {
				if (a != a) return a;
				if (b != b) return b;
				using std::min; return min(a, b);
			}
		};
		// SPFN 判断G: empty matrices are errors (dense mats style; MATLAB's
		// max([]) = [] is deliberately not followed).
		void spfn_require_nonempty_(const char* routine) const {
			if (row <= 0 || column <= 0) {
				vcp::throw_error<vcp::state_error>(routine, ": empty matrix");
			}
		}
		// SPFN 判断B: shared column-wise reduction core for max/min over a
		// CSC copy (as_csc(); *this stays untouched, matching the const
		// contract without relying on the mutable-format loophole).
		// Implicit zeros participate: a column with fewer stored entries
		// than `row` folds a _T(0) candidate in; a FULLY stored column must
		// NOT (a phantom zero would corrupt e.g. an all-negative column).
		template <class _Cmp> void spfn_colreduce_(spmats& B, _Cmp better) const {
			const spmats C = as_csc();
			B.resize(1, column);
			for (index_type j = 0; j < column; j++) {
				const index_type first = C.outer[static_cast<std::size_t>(j)];
				const index_type last  = C.outer[static_cast<std::size_t>(j) + 1];
				_T best;
				if (first == last) {
					best = _T(0);
				}
				else {
					best = C.value[static_cast<std::size_t>(first)];
					for (index_type p = first + 1; p < last; p++) {
						best = better(C.value[static_cast<std::size_t>(p)], best);
					}
					if (last - first < row) best = better(_T(0), best);
				}
				B.set(0, j, best);
			}
			B.finalize();
		}
		// Whole-vector reduction (row==1 || column==1) incl. implicit zeros.
		template <class _Cmp> void spfn_vecreduce_(spmats& B, _Cmp better) const {
			const spmats C = as_csr();
			const index_type len = (row == 1) ? column : row;
			_T best;
			if (C.value.empty()) {
				best = _T(0);
			}
			else {
				best = C.value[0];
				for (std::size_t k = 1; k < C.value.size(); k++) {
					best = better(C.value[k], best);
				}
				if (size_to_index(C.value.size(), "spmats::spfn_vecreduce_") < len) best = better(_T(0), best);
			}
			B.resize(1, 1);
			B.set(0, 0, best);
			B.finalize();
		}
	public:
		// A = abs(A) / sqrt(A) / sin(A): sparsity-preserving (f(0)=0).
		void abs()  { spfn_map_stored_(spfn_abs_());  }
		void sqrt() { spfn_map_stored_(spfn_sqrt_()); }
		void sin()  { spfn_map_stored_(spfn_sin_());  }
		// A = cos(A) / exp(A) / log(A): MATLAB sparse semantics (SPFN 判断A),
		// applied at every position; the result is structurally (near-)full
		// (row*column minus the exactly-zero images, e.g. log(1)=0).  log
		// evaluates log(_T(0)) at implicit zeros, exactly as dense
		// mats<T>::log does on zero entries.
		void cos()  { spfn_map_full_(spfn_cos_(), "spmats::cos"); }
		void exp()  { spfn_map_full_(spfn_exp_(), "spmats::exp"); }
		void log()  { spfn_map_full_(spfn_log_(), "spmats::log"); }

		// max/min (SPFN): same shape contract as mats<T> --
		//   1x1        -> 1x1 (the value itself)
		//   vector     -> 1x1 (whole-vector reduction)
		//   m x n      -> 1 x n (column-wise reduction, MATLAB style)
		// Implicit zeros participate in the comparisons (判断B), NaN
		// propagates (判断F), empty matrices throw (判断G).
		void max(spmats& B) const {
			spfn_require_nonempty_("spmats::max");
			if (row == 1 && column == 1) { B.resize(1, 1); B.set(0, 0, get(0, 0)); B.finalize(); return; }
			if (row == 1 || column == 1) { spfn_vecreduce_(B, spfn_takemax_()); return; }
			spfn_colreduce_(B, spfn_takemax_());
		}
		void min(spmats& B) const {
			spfn_require_nonempty_("spmats::min");
			if (row == 1 && column == 1) { B.resize(1, 1); B.set(0, 0, get(0, 0)); B.finalize(); return; }
			if (row == 1 || column == 1) { spfn_vecreduce_(B, spfn_takemin_()); return; }
			spfn_colreduce_(B, spfn_takemin_());
		}

		// normone (SPFN): 1x1 |a| / vector sum_i |v_i| / matrix operator
		// 1-norm (max column abs-sum).  B is 1x1, as in mats<T>.  Implicit
		// zeros contribute |0| = 0, so only stored entries are summed; an
		// empty column's sum is 0 and the max fold starts from the first
		// column's sum so empty axes are included correctly.  The max fold
		// uses spfn_takemax_ (判断F: NaN propagates; the sums themselves
		// propagate NaN through arithmetic).
		void normone(spmats& B) const {
			spfn_require_nonempty_("spmats::normone");
			using std::abs;
			_T res;
			if (row == 1 && column == 1) {
				res = abs(get(0, 0));
			}
			else if (row == 1 || column == 1) {
				const spmats C = as_csr();
				res = _T(0);
				for (std::size_t k = 0; k < C.value.size(); k++) res += abs(C.value[k]);
			}
			else {
				const spmats C = as_csc();
				res = _T(0);
				for (index_type j = 0; j < column; j++) {
					_T s = _T(0);
					for (index_type p = C.outer[static_cast<std::size_t>(j)];
					     p < C.outer[static_cast<std::size_t>(j) + 1]; p++) {
						s += abs(C.value[static_cast<std::size_t>(p)]);
					}
					res = (j == 0) ? s : spfn_takemax_()(s, res);   // 判断F: NaN propagates
				}
			}
			B.resize(1, 1);
			B.set(0, 0, res);
			B.finalize();
		}
		// norminf (SPFN): 1x1 |a| / vector max_i |v_i| / matrix operator
		// inf-norm (max row abs-sum).  B is 1x1, as in mats<T>.  The vector
		// fold starts from 0 (implicit zeros give |0| = 0 and abs >= 0
		// makes the 0-init exact even when nothing is stored).
		void norminf(spmats& B) const {
			spfn_require_nonempty_("spmats::norminf");
			using std::abs;
			_T res;
			if (row == 1 && column == 1) {
				res = abs(get(0, 0));
			}
			else if (row == 1 || column == 1) {
				const spmats C = as_csr();
				res = _T(0);
				for (std::size_t k = 0; k < C.value.size(); k++) res = spfn_takemax_()(abs(C.value[k]), res);   // 判断F
			}
			else {
				const spmats C = as_csr();
				res = _T(0);
				for (index_type i = 0; i < row; i++) {
					_T s = _T(0);
					for (index_type p = C.outer[static_cast<std::size_t>(i)];
					     p < C.outer[static_cast<std::size_t>(i) + 1]; p++) {
						s += abs(C.value[static_cast<std::size_t>(p)]);
					}
					res = (i == 0) ? s : spfn_takemax_()(s, res);   // 判断F: NaN propagates
				}
			}
			B.resize(1, 1);
			B.set(0, 0, res);
			B.finalize();
		}
		// normtwo (SPFN): destructive (self becomes 1x1), as mats<T>::normtwo.
		//   1x1    -> |a|
		//   vector -> sqrt(sum v_i^2)
		//   matrix -> sqrt(lambda_max(A^T A)) via policy_eigs on *this.
		// 判断C (most important): the A^T A product is assigned into the
		// BASE subobject of *this (`*this = AtA;` inside this spmats member
		// resolves to the base operator=, mirroring dense mats::normtwo's
		// "mulltmm -> (*this)=A -> eigsym").  This keeps the dynamic type
		// of *this, so the NVI chain policy_eigs -> policy_eigs_impl ->
		// policy_eigs_with_info -> policy_eigs_with_info_impl dispatches to
		// a derived policy's override (e.g. spumar's ARPACK path) when one
		// is installed.  Calling eigs on a fresh spmats<T> temporary
		// instead would lose that override and silently fall back to the
		// base Lanczos -- do NOT "simplify" this into a temporary.
		void normtwo() {
			spfn_require_nonempty_("spmats::normtwo");
			using std::abs;
			using std::sqrt;
			if (row == 1 && column == 1) {
				const _T res = abs(get(0, 0));
				this->resize(1, 1);
				this->set(0, 0, res);
				this->finalize();
				return;
			}
			if (row == 1 || column == 1) {
				const spmats C = as_csr();
				_T s = _T(0);
				for (std::size_t k = 0; k < C.value.size(); k++) {
					s += C.value[k] * C.value[k];
				}
				const _T res = sqrt(s);
				this->resize(1, 1);
				this->set(0, 0, res);
				this->finalize();
				return;
			}
			{
				const spmats At = this->transpose();
				const spmats AtA = this->policy_mul(At, *this);   // (n x m)(m x n) = n x n, symmetric PSD
				*this = AtA;   // base-subobject assignment; dynamic type preserved (判断C)
			}
			eig_options<_T> opt;
			opt.structure = matrix_structure_hint::symmetric;
			opt.target = eig_target::largest_magnitude;
			const std::vector<_T> ev = this->policy_eigs(1, opt);
			const _T res = sqrt(abs(ev[0]));   // abs guards a tiny negative rounding of lambda_max >= 0
			this->resize(1, 1);
			this->set(0, 0, res);
			this->finalize();
		}

		// matlab C = [A,B] -- concatenate columns (row counts must match).
		// Built directly from as_csr() copies of *this and B (see
		// sandbox/docs/misc/spmats_destructive_helpers_misc.md): as_csr() never
		// mutates *this/B, and is near-free when they are already CSR-finalized
		// (SLU-C2-FIX-IDEMPOTENT-FINALIZE). Within each output row, B's shifted
		// column indices (+= this->column) are all greater than A's, so the
		// merged row stays sorted with no re-sort/dedup pass needed.
		void horzcat(const spmats& B, spmats& C) const {
			if (row != B.row) {
				vcp::throw_error<vcp::dimension_error>("spmats::horzcat: row size mismatch: ", row, " != ", B.row);
			}
			const spmats Ac = this->as_csr();
			const spmats Bc = B.as_csr();
			const index_type out_rows = row;
			const index_type out_cols = column + B.column;
			const index_type nnzA = size_to_index(Ac.value.size(), "spmats::horzcat");
			const index_type nnzB = size_to_index(Bc.value.size(), "spmats::horzcat");

			C.row = out_rows;
			C.column = out_cols;
			C.outer.assign(index_to_size(checked_plus_one(out_rows, "spmats::horzcat"), "spmats::horzcat"), 0);
			C.inner.assign(index_to_size(nnzA + nnzB, "spmats::horzcat"), 0);
			C.value.assign(index_to_size(nnzA + nnzB, "spmats::horzcat"), _T(0));

			index_type pos = 0;
			for (index_type i = 0; i < out_rows; i++) {
				C.outer[i] = pos;
				for (index_type k = Ac.outer[i]; k < Ac.outer[i + 1]; k++) {
					C.inner[pos] = Ac.inner[k];
					C.value[pos] = Ac.value[k];
					pos++;
				}
				for (index_type k = Bc.outer[i]; k < Bc.outer[i + 1]; k++) {
					C.inner[pos] = Bc.inner[k] + column;
					C.value[pos] = Bc.value[k];
					pos++;
				}
			}
			C.outer[out_rows] = pos;

			C.coo_row.clear();
			C.coo_col.clear();
			C.coo_value.clear();
			C.coo_is_set.clear();
			C.fmt = vcp::sparse_csr;
			C.finalized = true;
			C.sorted = true;
			C.unique = true;
		}

		// matlab C = [A;B] -- concatenate rows (column counts must match).
		// Same as_csr()-copy approach as horzcat() above. Row pointers are
		// concatenated directly (B's row pointers shifted by A's total nnz,
		// i.e. Ac.outer[row]); column indices are untouched since rows don't
		// interleave columns, so no re-sort/dedup pass is needed either.
		void vercat(const spmats& B, spmats& C) const {
			if (column != B.column) {
				vcp::throw_error<vcp::dimension_error>("spmats::vercat: column size mismatch: ", column, " != ", B.column);
			}
			const spmats Ac = this->as_csr();
			const spmats Bc = B.as_csr();
			const index_type out_rows = row + B.row;
			const index_type out_cols = column;
			const index_type nnzA = size_to_index(Ac.value.size(), "spmats::vercat");
			const index_type nnzB = size_to_index(Bc.value.size(), "spmats::vercat");

			C.row = out_rows;
			C.column = out_cols;
			C.outer.assign(index_to_size(checked_plus_one(out_rows, "spmats::vercat"), "spmats::vercat"), 0);
			C.inner.assign(index_to_size(nnzA + nnzB, "spmats::vercat"), 0);
			C.value.assign(index_to_size(nnzA + nnzB, "spmats::vercat"), _T(0));

			for (index_type i = 0; i <= row; i++) {
				C.outer[i] = Ac.outer[i];
			}
			for (index_type i = 1; i <= B.row; i++) {
				C.outer[row + i] = Ac.outer[row] + Bc.outer[i];
			}
			for (index_type k = 0; k < nnzA; k++) {
				C.inner[k] = Ac.inner[k];
				C.value[k] = Ac.value[k];
			}
			for (index_type k = 0; k < nnzB; k++) {
				C.inner[nnzA + k] = Bc.inner[k];
				C.value[nnzA + k] = Bc.value[k];
			}

			C.coo_row.clear();
			C.coo_col.clear();
			C.coo_value.clear();
			C.coo_is_set.clear();
			C.fmt = vcp::sparse_csr;
			C.finalized = true;
			C.sorted = true;
			C.unique = true;
		}

		// SUB-1: sparse mirror of mats::submat -- same name, (a)-style out
		// parameter, non-virtual.  Selector grammar / validation order /
		// throw messages are an exact mirror of the dense implementation
		// (SUB-1_design.md v1.1 §1): {} whole axis, {i} single index,
		// {a,b} CLOSED range a..b (both ends included), {a,s,b} stride
		// a, a+s, ... while <= b; 0-based; vcp::index_error on reversed
		// range / negative index / stride < 1 / out-of-range / selector
		// with more than 3 elements.
		// Finalize policy: category 4 (self-contained) -- reads through an
		// as_csr()/as_csc() copy, *this is never modified, B is
		// born-finalized with the source format preserved (finalized CSC ->
		// CSC; finalized CSR and unfinalized COO -> CSR).  Single gather
		// pass: O(nnz of the touched outer lines + max(row, column) +
		// output nnz).
		void submat(spmats& B, const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
			submat_check_sizes_(list1, list2);
			const submat_axis_ ra = submat_normalize_(list1, row, true);
			const submat_axis_ ca = submat_normalize_(list2, column, false);
			const bool use_csc = finalized && fmt == vcp::sparse_csc;
			const spmats src = use_csc ? as_csc() : as_csr();
			const submat_axis_& oa = use_csc ? ca : ra;   // outer axis of the storage
			const submat_axis_& ia = use_csc ? ra : ca;   // inner axis of the storage
			const index_type inner_dim = use_csc ? row : column;
			// old inner index -> new inner index (-1 = not selected)
			std::vector<index_type> remap(index_to_size(inner_dim, "spmats::submat"), static_cast<index_type>(-1));
			for (index_type k = 0; k < ia.count; k++) {
				remap[static_cast<std::size_t>(ia.start + k * ia.stride)] = k;
			}
			std::vector<index_type> new_outer(index_to_size(checked_plus_one(oa.count, "spmats::submat"), "spmats::submat"), 0);
			std::vector<index_type> new_inner;
			std::vector<_T> new_value;
			for (index_type r = 0; r < oa.count; r++) {
				const std::size_t o = static_cast<std::size_t>(oa.start + r * oa.stride);
				for (index_type k = src.outer[o]; k < src.outer[o + 1]; k++) {
					const index_type m = remap[static_cast<std::size_t>(src.inner[static_cast<std::size_t>(k)])];
					if (m >= 0) {
						new_inner.push_back(m);
						new_value.push_back(src.value[static_cast<std::size_t>(k)]);
					}
				}
				new_outer[static_cast<std::size_t>(r) + 1] = size_to_index(new_inner.size(), "spmats::submat");
			}
			// The selected old indices are increasing and remap is monotone,
			// so each gathered line stays sorted and unique; stored values
			// are nonzero by the spmats invariant.  The assign_*
			// preconditions therefore hold constructively and B is
			// born-finalized in the preserved format.
			if (use_csc) {
				B.assign_csc(ra.count, ca.count, new_outer, new_inner, new_value);
			}
			else {
				B.assign_csr(ra.count, ca.count, new_outer, new_inner, new_value);
			}
		}

		void assign_csr(const index_type rows, const index_type cols,
		                const std::vector<index_type>& row_ptr,
		                const std::vector<index_type>& col_ind,
		                const std::vector<_T>& val) {
			validate_csr(rows, cols, row_ptr, col_ind, val);
			row = rows;
			column = cols;
			outer = row_ptr;
			inner = col_ind;
			value = val;
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			fmt = vcp::sparse_csr;
			finalized = true;
			sorted = true;
			unique = true;
		}

		void assign_csc(const index_type rows, const index_type cols,
		                const std::vector<index_type>& col_ptr,
		                const std::vector<index_type>& row_ind,
		                const std::vector<_T>& val) {
			validate_csc(rows, cols, col_ptr, row_ind, val);
			row = rows;
			column = cols;
			outer = col_ptr;
			inner = row_ind;
			value = val;
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			fmt = vcp::sparse_csc;
			finalized = true;
			sorted = true;
			unique = true;
		}

		const std::vector<index_type>& outer_index() const { return outer; }
		const std::vector<index_type>& inner_index() const { return inner; }
		const std::vector<_T>& values() const { return value; }
		const std::vector<index_type>& coo_rows() const { return coo_row; }
		const std::vector<index_type>& coo_columns() const { return coo_col; }
		const std::vector<_T>& coo_values() const { return coo_value; }

		// 非 const アクセサ(SPI-R9、オーナー裁定 2026-07-25): ポリシー層
		// (spmats)はアルゴリズムの作業層であり、ユーザ向けカプセル化は
		// spmatrix 層(protected 継承 + const 転送窓)が担う。格納配列への
		// 直接書き込みを許すことで、値のみの操作(符号反転・スケーリング等)を
		// 一時配列や assign_* 再構成なしの O(nnz) 1 パスで行える。
		// 【契約(呼び出し側の責務)】
		//   - values() の書き換えは値のみ(パターン不変)が原則。
		//   - outer/inner を書き換える場合は CSR/CSC 不変条件(ポインタ単調
		//     増加・行内昇順・重複なし・配列長整合)を保つこと。破れば以後の
		//     全操作が未定義。
		//   - 厳密ゼロの書き込みは invariant(格納値は非零)を破る。ゼロ化は
		//     パターンから除く操作(assign_* 再構成)で行うこと。
		//   - spmatrix 側でこれらを using 宣言により再公開することを禁ずる
		//     (const 転送窓 L1329-1334 がユーザ境界である)。
		std::vector<index_type>& outer_index() { return outer; }
		std::vector<index_type>& inner_index() { return inner; }
		std::vector<_T>&         values()      { return value; }

	protected:
		index_type row;
		index_type column;
		// finalize()/to_csr()/to_csc() are logically const (they change the
		// internal COO<->CSR/CSC representation but not the mathematical
		// value of the matrix); mutable lets them run from const methods
		// (see sandbox/docs/design/spmats_finalize_policy.md §4).
		mutable format_type fmt;
		mutable bool finalized;
		mutable bool sorted;
		mutable bool unique;

		mutable std::vector<index_type> outer;
		mutable std::vector<index_type> inner;
		mutable std::vector<_T> value;

		mutable std::vector<index_type> coo_row;
		mutable std::vector<index_type> coo_col;
		mutable std::vector<_T> coo_value;
		// Option A tag array (SLU-C3-ELEMENT-ACCESSOR): parallel to
		// coo_row/coo_col/coo_value (always same size). 1 = pushed by
		// set(), 0 = pushed by add(). unsigned char rather than
		// std::vector<bool>, which has no contiguous .data() to hand to
		// the tcoo_sort/tcoo_sum_duplicates_tagged C-style array API.
		mutable std::vector<unsigned char> coo_is_set;

		void clear_storage() {
			outer.clear();
			inner.clear();
			value.clear();
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			coo_is_set.clear();
		}

		void check_index(const index_type i, const index_type j, const char* routine) const {
			if (i < 0 || i >= row || j < 0 || j >= column) {
				vcp::throw_error<vcp::index_error>(routine, ": index out of range");
			}
		}

		// SUB-1: one selector axis normalized to start/stride/count in the
		// closed-interval grammar of mats::submat; full marks {} (whole
		// axis).  Shared by submat() above and by the spmatrix_block write
		// proxy (spmatrix.hpp) so that validation is complete before any
		// write begins.
		struct submat_axis_ {
			index_type start;
			index_type stride;
			index_type count;
			bool full;
		};

		// Mirrors the leading combined size check of mats::submat.
		static void submat_check_sizes_(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) {
			if (list1.size() > 3 || list2.size() > 3) {
				vcp::throw_error<vcp::index_error>(
					"submat: invalid selector size: ", list1.size(), ", ", list2.size());
			}
		}

		// Validation conditions, order and message text mirror mats::submat
		// exactly (the row selector is validated before the column selector
		// -- callers must keep that call order).
		static submat_axis_ submat_normalize_(const std::initializer_list<int>& list, const index_type dim, const bool is_row) {
			const std::vector<int> l = list;
			submat_axis_ a;
			a.full = false;
			if (l.size() == 0) {
				a.start = 0;
				a.stride = 1;
				a.count = dim;
				a.full = true;
			}
			else if (l.size() == 1) {
				if (l[0] < 0 || static_cast<index_type>(l[0]) >= dim) {
					if (is_row) {
						vcp::throw_error<vcp::index_error>("submat: row index out of range: ", l[0]);
					}
					vcp::throw_error<vcp::index_error>("submat: column index out of range: ", l[0]);
				}
				a.start = static_cast<index_type>(l[0]);
				a.stride = 1;
				a.count = 1;
			}
			else if (l.size() == 2) {
				if (l[0] > l[1] || l[0] < 0 || static_cast<index_type>(l[1]) >= dim) {
					if (is_row) {
						vcp::throw_error<vcp::index_error>(
							"submat: invalid row range: ", l[0], ":", l[1]);
					}
					vcp::throw_error<vcp::index_error>(
						"submat: invalid column range: ", l[0], ":", l[1]);
				}
				a.start = static_cast<index_type>(l[0]);
				a.stride = 1;
				a.count = static_cast<index_type>(l[1] - l[0] + 1);
			}
			else {
				if (l[0] > l[2] || l[0] < 0 || l[1] < 1 || static_cast<index_type>(l[2]) >= dim) {
					if (is_row) {
						vcp::throw_error<vcp::index_error>(
							"submat: invalid row range: ", l[0], ":", l[1], ":", l[2]);
					}
					vcp::throw_error<vcp::index_error>(
						"submat: invalid column range: ", l[0], ":", l[1], ":", l[2]);
				}
				index_type k = 0;
				for (int i = l[0]; i <= l[2]; i += l[1]) {
					k++;
				}
				a.start = static_cast<index_type>(l[0]);
				a.stride = static_cast<index_type>(l[1]);
				a.count = k;
			}
			return a;
		}

		void require_finalized(const char* routine) const {
			if (!finalized) {
				vcp::throw_error<vcp::state_error>(routine, ": matrix is not finalized");
			}
		}

		static std::size_t index_to_size(const index_type n, const char* routine) {
			if (n < 0) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": negative size");
			}
			return static_cast<std::size_t>(n);
		}

		static index_type size_to_index(const std::size_t n, const char* routine) {
			if (n > static_cast<std::size_t>((std::numeric_limits<index_type>::max)())) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": size exceeds index_type range");
			}
			return static_cast<index_type>(n);
		}

		static index_type checked_plus_one(const index_type n, const char* routine) {
			if (n < 0) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": negative size");
			}
			if (n == (std::numeric_limits<index_type>::max)()) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": index_type overflow");
			}
			return static_cast<index_type>(n + 1);
		}

		// NOT read-only: this transitions a finalized (CSR/CSC) matrix into
		// unfinalized COO as a side effect (clears outer/inner/value, sets
		// finalized=false, fmt=sparse_coo). Safe to call from add()/set()/
		// to_csr()/to_csc() (each restores its own target state
		// afterward). For any other read-only need, do not call this
		// directly on the target -- use as_csr()/as_csc() instead (a
		// copy-based operation, near-free once the source is already
		// finalized thanks to SLU-C2-FIX-IDEMPOTENT-FINALIZE). See
		// sandbox/docs/misc/spmats_destructive_helpers_misc.md.
		void ensure_coo_buffer() const {
			if (!finalized && fmt == vcp::sparse_coo) return;

			std::vector<index_type> nr;
			std::vector<index_type> nc;
			std::vector<_T> nv;
			const index_type n = size_to_index(value.size(), "spmats::ensure_coo_buffer");
			nr.resize(index_to_size(n, "spmats::ensure_coo_buffer"));
			nc.resize(index_to_size(n, "spmats::ensure_coo_buffer"));
			nv.resize(index_to_size(n, "spmats::ensure_coo_buffer"));
			// CSR/CSC values are already merged and unique, so the tag is
			// irrelevant here (see get()'s reset-from-zero equivalence);
			// zero-fill just to keep the array in lockstep.
			std::vector<unsigned char> ni(index_to_size(n, "spmats::ensure_coo_buffer"), 0);

			if (fmt == vcp::sparse_csr) {
				vcp::tcsr_to_coo(row, column, n, outer.data(), inner.data(), value.data(), nr.data(), nc.data(), nv.data());
			}
			else if (fmt == vcp::sparse_csc) {
				index_type p = 0;
				for (index_type j = 0; j < column; j++) {
					for (index_type k = outer[j]; k < outer[j + 1]; k++) {
						nr[static_cast<std::size_t>(p)] = inner[k];
						nc[static_cast<std::size_t>(p)] = j;
						nv[static_cast<std::size_t>(p)] = value[k];
						p++;
					}
				}
			}
			else {
				vcp::throw_error<vcp::state_error>("spmats::ensure_coo_buffer: unsupported sparse format");
			}
			coo_row.swap(nr);
			coo_col.swap(nc);
			coo_value.swap(nv);
			coo_is_set.swap(ni);
			outer.clear();
			inner.clear();
			value.clear();
			fmt = vcp::sparse_coo;
			finalized = false;
			sorted = true;
			unique = true;
		}

		static void validate_csr(const index_type rows, const index_type cols,
		                         const std::vector<index_type>& row_ptr,
		                         const std::vector<index_type>& col_ind,
		                         const std::vector<_T>& val) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: negative size");
			}
			if (row_ptr.size() != index_to_size(checked_plus_one(rows, "spmats::assign_csr"), "spmats::assign_csr") || col_ind.size() != val.size()) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: invalid array size");
			}
			if (row_ptr.empty() || row_ptr[0] != 0 || row_ptr.back() != size_to_index(col_ind.size(), "spmats::assign_csr")) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: invalid row_ptr boundary");
			}
			for (index_type i = 0; i < rows; i++) {
				if (row_ptr[i] > row_ptr[i + 1]) {
					vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: row_ptr must be monotone");
				}
				index_type previous = -1;
				for (index_type p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
					if (col_ind[p] < 0 || col_ind[p] >= cols) {
						vcp::throw_error<vcp::index_error>("spmats::assign_csr: col_ind out of range");
					}
					if (col_ind[p] <= previous) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: columns must be sorted and unique in each row");
					}
					if (val[p] == _T(0)) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: explicit zero is not allowed");
					}
					previous = col_ind[p];
				}
			}
		}

		static void validate_csc(const index_type rows, const index_type cols,
		                         const std::vector<index_type>& col_ptr,
		                         const std::vector<index_type>& row_ind,
		                         const std::vector<_T>& val) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: negative size");
			}
			if (col_ptr.size() != index_to_size(checked_plus_one(cols, "spmats::assign_csc"), "spmats::assign_csc") || row_ind.size() != val.size()) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: invalid array size");
			}
			if (col_ptr.empty() || col_ptr[0] != 0 || col_ptr.back() != size_to_index(row_ind.size(), "spmats::assign_csc")) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: invalid col_ptr boundary");
			}
			for (index_type j = 0; j < cols; j++) {
				if (col_ptr[j] > col_ptr[j + 1]) {
					vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: col_ptr must be monotone");
				}
				index_type previous = -1;
				for (index_type p = col_ptr[j]; p < col_ptr[j + 1]; p++) {
					if (row_ind[p] < 0 || row_ind[p] >= rows) {
						vcp::throw_error<vcp::index_error>("spmats::assign_csc: row_ind out of range");
					}
					if (row_ind[p] <= previous) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: rows must be sorted and unique in each column");
					}
					if (val[p] == _T(0)) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: explicit zero is not allowed");
					}
					previous = row_ind[p];
				}
			}
		}
		// ------------------------------------------------------------------
		// Symmetry check (used by eigs dispatch and lss)
		// ------------------------------------------------------------------
	public:
		typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;

		bool is_symmetric() const {
			return is_symmetric(vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(12));
		}

		bool is_symmetric(const scalar_real_type& tol) const {
			if (tol <= scalar_real_type(0))
				vcp::throw_error<vcp::invalid_argument>("spmats::is_symmetric: tol must be positive");
			if (row != column) return false;
			spmats A = this->as_csr();
			const std::vector<_Index>& outer = A.outer;
			const std::vector<_Index>& inner = A.inner;
			const std::vector<_T>& val = A.value;
			for (_Index i = 0; i < A.row; i++) {
				for (_Index p = outer[static_cast<std::size_t>(i)];
				     p < outer[static_cast<std::size_t>(i + 1)]; p++) {
					const _Index j = inner[static_cast<std::size_t>(p)];
					if (i == j) continue;
					const _Index first = outer[static_cast<std::size_t>(j)];
					const _Index last = outer[static_cast<std::size_t>(j + 1)];
					const typename std::vector<_Index>::const_iterator begin = inner.begin() + first;
					const typename std::vector<_Index>::const_iterator end = inner.begin() + last;
					typename std::vector<_Index>::const_iterator it = std::lower_bound(begin, end, i);
					_T mirrored = _T(0);
					if (it != end && *it == i)
						mirrored = val[static_cast<std::size_t>(it - inner.begin())];
					const scalar_real_type diff = vcp::tsparse_scalar::abs_value(
						val[static_cast<std::size_t>(p)] - mirrored);
					if (diff > tol) return false;
				}
			}
			return true;
		}

		// ------------------------------------------------------------------
		// Phase 7.7: policy methods for to_dense / is_symmetric forwarding
		// spmatrix<T,P> delegates A.to_dense() / A.is_symmetric() here.
		// Custom policy P can override these to change dense conversion or
		// symmetry-check semantics (e.g. verified / interval enclosure).
		// ------------------------------------------------------------------
		typedef std::vector<std::vector<_T> > dense_matrix_type;

		// policy_to_dense: CSR traversal → row-major dense 2-D vector.
		// Missing entries become T(0).  Behaviour matches the pre-7.7
		// spmatrix::to_dense() implementation.
		template <class Matrix>
		dense_matrix_type policy_to_dense(const Matrix& A) const {
			dense_matrix_type dense(
				static_cast<std::size_t>(A.rowsize()),
				std::vector<_T>(static_cast<std::size_t>(A.columnsize()), _T(0)));
			Matrix Acsr = A.as_csr();
			const std::vector<_Index>& outerv = Acsr.outer_index();
			const std::vector<_Index>& innerv = Acsr.inner_index();
			const std::vector<_T>& val = Acsr.values();
			for (_Index i = 0; i < Acsr.rowsize(); i++) {
				for (_Index p = outerv[static_cast<std::size_t>(i)];
				     p < outerv[static_cast<std::size_t>(i + 1)]; p++) {
					dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(innerv[static_cast<std::size_t>(p)])] =
						val[static_cast<std::size_t>(p)];
				}
			}
			return dense;
		}

		// policy_is_symmetric (no-tol): delegates to tolerance overload.
		template <class Matrix>
		bool policy_is_symmetric(const Matrix& A) const {
			return policy_is_symmetric(A, vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(12));
		}

		// policy_is_symmetric (with tol): complex-symmetric check, NOT Hermitian.
		// For complex T, compares A(i,j) with A(j,i), not conj(A(j,i)).
		// Behaviour matches the pre-7.7 spmatrix::is_symmetric(tol) implementation.
		template <class Matrix>
		bool policy_is_symmetric(const Matrix& A, const scalar_real_type& tol) const {
			if (tol <= scalar_real_type(0))
				vcp::throw_error<vcp::invalid_argument>("spmats::policy_is_symmetric: tol must be positive");
			if (A.rowsize() != A.columnsize()) return false;
			Matrix Acsr = A.as_csr();
			const std::vector<_Index>& outerv = Acsr.outer_index();
			const std::vector<_Index>& innerv = Acsr.inner_index();
			const std::vector<_T>& val = Acsr.values();
			for (_Index i = 0; i < Acsr.rowsize(); i++) {
				for (_Index p = outerv[static_cast<std::size_t>(i)];
				     p < outerv[static_cast<std::size_t>(i + 1)]; p++) {
					const _Index j = innerv[static_cast<std::size_t>(p)];
					if (i == j) continue;
					const _Index first = outerv[static_cast<std::size_t>(j)];
					const _Index last  = outerv[static_cast<std::size_t>(j + 1)];
					const typename std::vector<_Index>::const_iterator begin = innerv.begin() + first;
					const typename std::vector<_Index>::const_iterator end   = innerv.begin() + last;
					typename std::vector<_Index>::const_iterator it = std::lower_bound(begin, end, i);
					_T mirrored = _T(0);
					if (it != end && *it == i)
						mirrored = val[static_cast<std::size_t>(it - innerv.begin())];
					const scalar_real_type diff = vcp::tsparse_scalar::abs_value(
						val[static_cast<std::size_t>(p)] - mirrored);
					if (diff > tol) return false;
				}
			}
			return true;
		}

		// ------------------------------------------------------------------
		// Private helpers for policy methods
		// ------------------------------------------------------------------
	private:
		// solve_* private helpers (called from policy_lss_with_info)
		linear_solve_result<_T> policy_solve_jacobi_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol, bool use_relative) const;
		linear_solve_result<_T> policy_solve_gauss_seidel_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol, bool use_relative) const;
		linear_solve_result<_T> policy_solve_cg_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol,
			bool check_symmetric, preconditioner_type prec, bool use_relative) const;
		linear_solve_result<_T> policy_solve_bicgstab_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol,
			bool use_relative, preconditioner_type prec) const;
		linear_solve_result<_T> policy_solve_gmres_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol,
			std::size_t restart, bool use_relative, preconditioner_type prec) const;

	public:
		// ------------------------------------------------------------------
		// Policy methods: arithmetic
		//
		// finalize safety (see sandbox/docs/design/spmats_finalize_policy.md):
		// policy_mul is NVI-split (non-virtual outer + virtual _impl) because
		// spgemm is an algorithm future policies may want to replace wholesale.
		// policy_add/policy_sub/policy_mul_vec/policy_left_mul_vec get a plain
		// auto-finalize check at the entry point (no virtual split needed).
		// policy_scalar_mul/policy_scalar_div/policy_neg are unchanged: they
		// operate on an as_csr() copy and are correct regardless of the
		// input's finalize state (self-contained, category 4).
		// ------------------------------------------------------------------
		spmats<_T,_Index> policy_add(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;
		spmats<_T,_Index> policy_sub(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;

		// policy_mul: non-virtual outer. Must never be overridden; override
		// policy_mul_impl instead.
		spmats<_T,_Index> policy_mul(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;
		virtual spmats<_T,_Index> policy_mul_impl(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;

		std::vector<_T> policy_mul_vec(const spmats<_T,_Index>& A, const std::vector<_T>& x) const;
		std::vector<_T> policy_left_mul_vec(const std::vector<_T>& x, const spmats<_T,_Index>& A) const;
		spmats<_T,_Index> policy_scalar_mul(const spmats<_T,_Index>& A, const _T& alpha) const;
		spmats<_T,_Index> policy_scalar_div(const spmats<_T,_Index>& A, const _T& alpha) const;
		spmats<_T,_Index> policy_neg(const spmats<_T,_Index>& A) const;

		// SPC-P1(オーナー裁定 2026-07-26): パターン不変演算の破壊的
		// (in-place)ポリシー。追加確保ゼロで *this の値配列を直接更新する。
		// 厳密ゼロが生じた要素(alpha==0・アンダーフロー・div のゼロ化)は
		// その場で前詰め圧縮しパターンから除く(invariant: 格納値は非零)。
		// 未 finalize(COO)の場合は先に finalize() する。
		void policy_mulsm(const _T& alpha);   // *this *= alpha
		void policy_divms(const _T& alpha);   // *this /= alpha
		void policy_minusm();                 // *this = -*this

		// ------------------------------------------------------------------
		// Policy methods: linear system solve
		//
		// policy_lss_with_info: non-virtual outer, NVI pattern (see design doc
		// §3). Must never be overridden; override policy_lss_with_info_impl
		// instead. solve_jacobi/gauss_seidel/cg/bicgstab/gmres (+_with_info)
		// and policy_lss all route through policy_lss_with_info, so they
		// inherit the finalize guarantee automatically.
		// WFIX-2: the subject matrix is *this (dense-side mats precedent);
		// the former leading `const spmats& A` argument is removed from every
		// solve/factorize/scan policy method below.
		// ------------------------------------------------------------------
		linear_solve_result<_T> policy_lss_with_info(
			const std::vector<_T>& b,
			const linear_solve_options<_T>& opt) const;
		virtual linear_solve_result<_T> policy_lss_with_info_impl(
			const std::vector<_T>& b,
			const linear_solve_options<_T>& opt) const;
		std::vector<_T> policy_lss(
			const std::vector<_T>& b,
			const linear_solve_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: reusable LU factorization handle (LSS-1 P-4;
		// definitions at the end of spmats_base/spmats_lss.hpp).
		//
		// policy_lu_factorize_with_info: non-virtual outer, NVI pattern
		// (finalize guarantee + squareness check).  Must never be
		// overridden; override policy_lu_factorize_with_info_impl instead.
		// Takes sparse_lu_options directly (no linear_solve_options: the
		// method concept does not apply to an explicit factorization).
		// unsigned Index throws vcp::state_error (signed-only, same SFINAE
		// split as the sparse_lu solve dispatch).
		// ------------------------------------------------------------------
		lu_factor_handle<_T,_Index> policy_lu_factorize_with_info(
			const sparse_lu_options<_T>& opt = sparse_lu_options<_T>()) const;
		virtual lu_factor_handle<_T,_Index> policy_lu_factorize_with_info_impl(
			const sparse_lu_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: eigenvalue _with_info (non-throwing, definitions in spmats_eigs.hpp)
		//
		// policy_eigs_with_info / policy_generalized_eigs_with_info (no
		// Preconditioner): NVI pattern, non-virtual outer + virtual _impl.
		// Must never override the outer; override the _impl instead.
		// Preconditioner overloads (templates) cannot be virtual in C++; they
		// get a plain auto-finalize check on *this (and B) at the entry point
		// instead of an _impl split.
		// WFIX-2: the subject matrix is *this; the generalized forms keep B
		// (operand) only.
		// ------------------------------------------------------------------
		eig_result<_T> policy_eigs_with_info(
			std::size_t k,
			const eig_options<_T>& opt) const;
		virtual eig_result<_T> policy_eigs_with_info_impl(
			std::size_t k,
			const eig_options<_T>& opt) const;

		template <class Prec>
		eig_result<_T> policy_eigs_with_info(
			std::size_t k,
			const eig_options<_T>& opt, const Prec& M) const;

		eig_result<_T> policy_generalized_eigs_with_info(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;
		virtual eig_result<_T> policy_generalized_eigs_with_info_impl(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;

		template <class Prec>
		eig_result<_T> policy_generalized_eigs_with_info(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt, const Prec& M) const;

		// ------------------------------------------------------------------
		// Policy methods: eigenvalue strict (throwing on failure)
		// These own all convergence/count checking; spmatrix.hpp does none.
		//
		// policy_eig / policy_eigs / policy_generalized_eig /
		// policy_generalized_eigs (no Preconditioner): NVI pattern, same rule
		// as above (outer never overridden, override the _impl).
		// ------------------------------------------------------------------
		eig_result<_T> policy_eig(
			const eig_options<_T>& opt) const;
		virtual eig_result<_T> policy_eig_impl(
			const eig_options<_T>& opt) const;

		std::vector<_T> policy_eigs(
			std::size_t k,
			const eig_options<_T>& opt) const;
		virtual std::vector<_T> policy_eigs_impl(
			std::size_t k,
			const eig_options<_T>& opt) const;

		template <class Prec>
		std::vector<_T> policy_eigs(
			std::size_t k,
			const eig_options<_T>& opt, const Prec& M) const;

		eig_result<_T> policy_generalized_eig(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;
		virtual eig_result<_T> policy_generalized_eig_impl(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;

		std::vector<_T> policy_generalized_eigs(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;
		virtual std::vector<_T> policy_generalized_eigs_impl(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;

		template <class Prec>
		std::vector<_T> policy_generalized_eigs(
			const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt, const Prec& M) const;

		// ------------------------------------------------------------------
		// Policy methods: LDL^T factorization (LDL-3, design v2 SS5.1)
		//
		// policy_ldl_with_info: non-virtual outer, NVI pattern (finalize
		// guarantee + squareness entry check).  Must never be overridden;
		// override policy_ldl_with_info_impl instead.  Convention:
		// P^T A P = L D L^T with perm new->old (design v2 SS5.4); L unit
		// lower (explicit unit diagonal), D 1x1/2x2 block diagonal with
		// structural zeros not stored (design v2 SS5.2).  L / D / perm are
		// valid outputs only when the returned status is success or
		// zero_pivot.  Definitions in spmats_base/spmats_ldl_impl.hpp.
		// ------------------------------------------------------------------
		ldl_result<_T,_Index> policy_ldl_with_info(
			spmats<_T,_Index>& L, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const ldl_options<_T>& opt) const;
		virtual ldl_result<_T,_Index> policy_ldl_with_info_impl(
			spmats<_T,_Index>& L, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const ldl_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: LL^T Cholesky factorization (CHOL-2, chol design
		// v1 SS5)
		//
		// policy_chol_with_info: non-virtual outer, NVI pattern (finalize
		// guarantee + squareness entry check).  Must never be overridden;
		// override policy_chol_with_info_impl instead (the designated
		// replacement point for external backends, e.g. a future CHOLMOD
		// delegation; design SS7.6).  Convention: P^T A P = L L^T with perm
		// new->old and P(p[k],k) = 1 (same orientation as LDL); L NON-unit
		// lower triangular with positive diagonal (l_kk = sqrt of the
		// certified-positive pivot); perm is the ordering output itself
		// (no pivoting; design SS5.6).  Exact zeros are not stored (the
		// certified-zero drop happens at this boundary; kernel nnz_L is the
		// symbolic count, so the materialized L stores <= nnz_L entries).
		// L / perm are valid outputs ONLY when the returned status is
		// success (D-3: chol has no LDL-style "completed failure").
		// Definitions in spmats_base/spmats_chol_impl.hpp.
		// ------------------------------------------------------------------
		chol_result<_T,_Index> policy_chol_with_info(
			spmats<_T,_Index>& L,
			std::vector<_Index>& perm,
			const chol_options<_T>& opt) const;
		virtual chol_result<_T,_Index> policy_chol_with_info_impl(
			spmats<_T,_Index>& L,
			std::vector<_Index>& perm,
			const chol_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: inertia (LDL-4, design v2 SS7; decision 4)
		//
		// policy_inertia_with_info: non-virtual outer (finalize + squareness),
		// NVI as above.  The default _impl calls policy_ldl_with_info and
		// scans the returned D; derived policies may override the _impl to
		// replace the computation.  policy_inertia_from_block_diagonal is the
		// non-virtual shared D-consumer (does NOT assume BK origin): certified
		// sign counting over a general 1x1/2x2 block diagonal, structurally
		// empty columns counted as 1x1 zero blocks (SS7.2 rules 1-5).
		// Definitions in spmats_base/spmats_ldl_impl.hpp.
		// ------------------------------------------------------------------
		inertia_result<_Index> policy_inertia_with_info(
			const inertia_options<_T>& opt) const;
		virtual inertia_result<_Index> policy_inertia_with_info_impl(
			const inertia_options<_T>& opt) const;
		// WFIX-2 (W2-3): the scanned block diagonal D is the SUBJECT --
		// call as D.policy_inertia_from_block_diagonal(tol).
		inertia_result<_Index> policy_inertia_from_block_diagonal(
			const scalar_real_type& tol) const;

		// ------------------------------------------------------------------
		// Policy methods: A - sigma*B LDL^T shift handle (SLDL-SH / B2)
		//
		// policy_ldl_shift_setup_with_info: non-virtual outers, NVI pattern
		// (finalize guarantee + squareness / dimension entry checks).  Must
		// never be overridden; override the _impl instead (derived policies
		// inherit the default unchanged).  The overload without B sets up
		// the A - sigma*I iteration; the overload with B requires B to be
		// certified symmetric only (never SPD-tested, H-4).  Every sigma-
		// independent stage (conversion, symmetry checks, pattern merge,
		// ordering, symbolic) runs once here; the returned handle is
		// read-only afterwards and const-shareable across threads.  Failure
		// reporting is info-only through the handle status (H-2; no strict
		// variant exists).  Definitions in spmats_base/spmats_ldl_shift.hpp
		// (phase 2 of its two-phase include).
		// ------------------------------------------------------------------
		ldl_shift_handle<_T,_Index> policy_ldl_shift_setup_with_info(
			const ldl_options<_T>& opt) const;
		virtual ldl_shift_handle<_T,_Index> policy_ldl_shift_setup_with_info_impl(
			const ldl_options<_T>& opt) const;
		ldl_shift_handle<_T,_Index> policy_ldl_shift_setup_with_info(
			const spmats<_T,_Index>& B, const ldl_options<_T>& opt) const;
		virtual ldl_shift_handle<_T,_Index> policy_ldl_shift_setup_with_info_impl(
			const spmats<_T,_Index>& B, const ldl_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: LU factor extraction (LUX-1, lux_design_v0 SS2)
		//
		// policy_lu_with_info: non-virtual outer, NVI pattern (finalize
		// guarantee + squareness entry check).  Must never be overridden;
		// override policy_lu_with_info_impl instead.  Convention (SSC):
		// P A Q = L U with p / q new->old (A(p,q) = L U, MATLAB
		// [L,U,P,Q] = lu(A) orientation) and P(k,p[k]) = 1, Q(q[k],k) = 1
		// -- note the row side is TRANSPOSED relative to the LDL convention
		// (LDL's left factor is P^T, here it is P itself).  L unit lower
		// (explicit unit diagonal), U upper; exact zeros not stored.
		// equilibration == true in the options is rejected with
		// unsupported_options (the SSC form has no scaling).  L / U / p / q
		// are valid outputs only when the returned status is success.
		// The virtual _impl is the designated replacement point for
		// external-backend policies (spumar / UMFPACK delegation).
		// Definitions in spmats_base/spmats_lu_extract_impl.hpp.
		// ------------------------------------------------------------------
		lu_extract_result<_T,_Index> policy_lu_with_info(
			spmats<_T,_Index>& L, spmats<_T,_Index>& U,
			std::vector<_Index>& p, std::vector<_Index>& q,
			const lu_extract_options<_T>& opt) const;
		virtual lu_extract_result<_T,_Index> policy_lu_with_info_impl(
			spmats<_T,_Index>& L, spmats<_T,_Index>& U,
			std::vector<_Index>& p, std::vector<_Index>& q,
			const lu_extract_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: LU factor consumers (LUX-2, lux_design_v0 SS2a)
		//
		// Both take SSC-convention factors (P A Q = L U, p/q new->old,
		// L unit lower with explicit unit diagonal, U upper) as ARGUMENTS
		// and consume nothing else -- the default _impl reads no policy
		// state (P-7), so any derived policy (spumar included) inherits it
		// unchanged, and a derived policy may override the _impl to replace
		// the computation (same extension-point design as inertia).
		//   policy_lu_solve_with_info:        x = Q U^{-1} L^{-1} P b
		//   policy_lu_inverse_row_with_info:  row_i(A^{-1})^T =
		//       P^T L^{-T} U^{-T} Q^T e_i  (transposed triangular solves,
		//       CSC arrays read row-wise; no transpose is materialized)
		// Non-virtual outers own the finalize guarantee and ALL dimension /
		// index-range checks (reported as dimension_mismatch, non-throwing);
		// the default _impl validates the structural contract (p/q
		// bijections, L unit lower, U upper -> invalid_input) and gates
		// every U-diagonal division through the certified three-branch
		// (P-9: nonzero certified -> divide / zero certified ->
		// singular_factor / undecidable -> inconclusive_division).
		// x / row are valid outputs only when the status is success.
		// Definitions in spmats_base/spmats_lu_extract_impl.hpp.
		// ------------------------------------------------------------------
		lu_apply_result policy_lu_solve_with_info(
			const spmats<_T,_Index>& L, const spmats<_T,_Index>& U,
			const std::vector<_Index>& p, const std::vector<_Index>& q,
			const std::vector<_T>& b, std::vector<_T>& x) const;
		virtual lu_apply_result policy_lu_solve_with_info_impl(
			const spmats<_T,_Index>& L, const spmats<_T,_Index>& U,
			const std::vector<_Index>& p, const std::vector<_Index>& q,
			const std::vector<_T>& b, std::vector<_T>& x) const;

		lu_apply_result policy_lu_inverse_row_with_info(
			const spmats<_T,_Index>& L, const spmats<_T,_Index>& U,
			const std::vector<_Index>& p, const std::vector<_Index>& q,
			const _Index i, std::vector<_T>& row) const;
		virtual lu_apply_result policy_lu_inverse_row_with_info_impl(
			const spmats<_T,_Index>& L, const spmats<_T,_Index>& U,
			const std::vector<_Index>& p, const std::vector<_Index>& q,
			const _Index i, std::vector<_T>& row) const;

		// ------------------------------------------------------------------
		// WFIX: matrix-form overloads -- factor materialization moved INTO
		// the policy layer (spmatrix is a pure forwarding wrapper; the P/Q
		// construction loops formerly in spmatrix::ldl_with_info /
		// lu_with_info were a layer-discipline violation).
		//
		// W-1 (LU, NVI): materialization itself is a replacement point --
		// backends differ in their natural factor representation (own
		// baseline / supernodal, UMFPACK, future SuperLU), so a derived
		// policy may either (a) override only the vector-form
		// policy_lu_with_info_impl (the matrix form then works through the
		// default materialization below), or (b) override
		// policy_lu_matrices_with_info_impl to build the matrices directly
		// from its internal representation.  Never override the outer.
		// Convention (SSC): P(k, p[k]) = 1, Q(q[k], k) = 1; P, Q returned
		// finalized; valid outputs only when status == success.
		//
		// W-2 (LDL, non-virtual helper): the LDL factorization is own-code
		// only, so materialization has no replacement demand -- a plain
		// overload calling the existing vector-form outer.  Convention
		// (LDL design v2 SS5.4): P(p[k], k) = 1 -- NOTE the row side is
		// TRANSPOSED relative to the LU convention above.  P returned
		// finalized; valid when status is success or zero_pivot (same
		// validity rule as the vector-form L/D/perm outputs).
		// ------------------------------------------------------------------

		// W-1 outer (non-virtual, matrix-form overload).  Must never be
		// overridden -- override policy_lu_matrices_with_info_impl instead.
		lu_extract_result<_T,_Index> policy_lu_with_info(
			spmats<_T,_Index>& L, spmats<_T,_Index>& U,
			spmats<_T,_Index>& P, spmats<_T,_Index>& Q,
			const lu_extract_options<_T>& opt) const {
			return policy_lu_matrices_with_info_impl(L, U, P, Q, opt);
		}

		// W-1 replacement point.  Default: vector-form outer (-> existing
		// vector-form virtual _impl; contract unchanged), then SSC
		// materialization P(k,p[k]) = 1, Q(q[k],k) = 1, finalized.
		virtual lu_extract_result<_T,_Index> policy_lu_matrices_with_info_impl(
			spmats<_T,_Index>& L, spmats<_T,_Index>& U,
			spmats<_T,_Index>& P, spmats<_T,_Index>& Q,
			const lu_extract_options<_T>& opt) const {
			std::vector<_Index> p, q;
			lu_extract_result<_T,_Index> result =
			    policy_lu_with_info(L, U, p, q, opt);
			P.resize(_Index(0), _Index(0));
			Q.resize(_Index(0), _Index(0));
			if (result.status == sparse_lu_extract_status::success) {
				const _Index n = static_cast<_Index>(p.size());
				P.resize(n, n);
				Q.resize(n, n);
				for (_Index k = 0; k < n; k++) {
					P.add(k, p[static_cast<std::size_t>(k)], _T(1));
					Q.add(q[static_cast<std::size_t>(k)], k, _T(1));
				}
				P.finalize();
				Q.finalize();
			}
			return result;
		}

		// W-2 (non-virtual, matrix-form overload): vector-form outer
		// (-> existing virtual _impl), then LDL materialization
		// P(p[k], k) = 1, finalized.
		ldl_result<_T,_Index> policy_ldl_with_info(
			spmats<_T,_Index>& L, spmats<_T,_Index>& D,
			spmats<_T,_Index>& P,
			const ldl_options<_T>& opt) const {
			std::vector<_Index> perm;
			ldl_result<_T,_Index> result =
			    policy_ldl_with_info(L, D, perm, opt);
			P.resize(_Index(0), _Index(0));
			if (result.status == sparse_ldl_status::success ||
			    result.status == sparse_ldl_status::zero_pivot) {
				const _Index n = static_cast<_Index>(perm.size());
				P.resize(n, n);
				for (_Index k = 0; k < n; k++) {
					P.add(perm[static_cast<std::size_t>(k)], k, _T(1));
				}
				P.finalize();
			}
			return result;
		}

		// CHOL-2 (non-virtual, matrix-form overload, W-2 pattern): the
		// Cholesky factorization is own-code only, so materialization has
		// no replacement demand -- a plain overload calling the vector-form
		// outer (-> virtual _impl).  Convention (chol design v1 SS1.1, same
		// orientation as LDL): P(p[k], k) = 1.  P returned finalized; valid
		// ONLY when status is success (simpler than LDL: chol has no
		// zero_pivot-style completed failure, design SS5).
		chol_result<_T,_Index> policy_chol_with_info(
			spmats<_T,_Index>& L,
			spmats<_T,_Index>& P,
			const chol_options<_T>& opt) const {
			std::vector<_Index> perm;
			chol_result<_T,_Index> result =
			    policy_chol_with_info(L, perm, opt);
			P.resize(_Index(0), _Index(0));
			if (result.status == sparse_chol_status::success) {
				const _Index n = static_cast<_Index>(perm.size());
				P.resize(n, n);
				for (_Index k = 0; k < n; k++) {
					P.add(perm[static_cast<std::size_t>(k)], k, _T(1));
				}
				P.finalize();
			}
			return result;
		}

		// ------------------------------------------------------------------
		// Policy methods: AINV approximate inverse (AINV-1, ainv design v1.2)
		//
		// policy_ainv_with_info: non-virtual outer, NVI pattern (finalize
		// guarantee + squareness entry check -- NON-throwing: a non-square
		// input is reported as invalid_input, design §1.3, unlike the
		// throwing ldl / lu outers).  Must never be overridden; override
		// policy_ainv_with_info_impl instead.  Incomplete biconjugation
		// after [BT98] (Benzi/Tuma, SIAM J. Sci. Comput. 19(3), 1998):
		// outputs Z / W unit upper triangular (explicit unit diagonal) and
		// D diagonal with every pivot lifted to nonzero (design §4.3), all
		// born-finalized.  R = Z D^{-1} W^T is NEVER materialized (D-13);
		// the apply / estimate methods below consume the factors directly.
		// Numerical events (tiny pivots, fill growth) are never a failure
		// status (D-4); a singular input completes with success (D-3).
		// API is _with_info only -- no strict throwing sugar (D-17).
		// Definitions in spmats_base/spmats_ainv_impl.hpp.
		// ------------------------------------------------------------------
		ainv_result<_T,_Index> policy_ainv_with_info(
			spmats<_T,_Index>& Z, spmats<_T,_Index>& W,
			spmats<_T,_Index>& D,
			const ainv_options<_T>& opt) const;
		virtual ainv_result<_T,_Index> policy_ainv_with_info_impl(
			spmats<_T,_Index>& Z, spmats<_T,_Index>& W,
			spmats<_T,_Index>& D,
			const ainv_options<_T>& opt) const;

		// policy_ainv_apply: z = Z (D^{-1} (W^T r)) as three sparse stages
		// ([BT98] eq. (6) structure); W^T is not materialized -- the W^T
		// product is a scatter scan over W's stored lines (design §2.4).
		// Consumer-type NVI (policy_lu_solve_with_info precedent): the
		// factors are ARGUMENTS, the non-virtual outer owns the finalize
		// guarantee and ALL dimension / structure checks (invalid_input,
		// non-throwing); override the _impl only.  z is a valid output only
		// when the returned status is success.
		ainv_status policy_ainv_apply(
			const spmats<_T,_Index>& Z, const spmats<_T,_Index>& W,
			const spmats<_T,_Index>& D,
			const std::vector<_T>& r, std::vector<_T>& z) const;
		virtual ainv_status policy_ainv_apply_impl(
			const spmats<_T,_Index>& Z, const spmats<_T,_Index>& W,
			const spmats<_T,_Index>& D,
			const std::vector<_T>& r, std::vector<_T>& z) const;

		// policy_ainv_residual_norm_estimate: NON-GUARANTEED estimate of
		// ||I - Z D^{-1} W^T A||_inf in plain T point arithmetic (D-11),
		// row by row and factor-based (R never materialized, working
		// memory O(n); design §2.3).  (a)-type: A is *this, the factors
		// are arguments.  Single non-virtual method (design §2.3;
		// consumer-type checks inline, non-throwing).
		ainv_status policy_ainv_residual_norm_estimate(
			const spmats<_T,_Index>& Z, const spmats<_T,_Index>& W,
			const spmats<_T,_Index>& D,
			typename vcp::tsparse_scalar::real_type<_T>::type& est) const;

		// ------------------------------------------------------------------
		// Policy methods: static FSAI factored approximate inverse (FSAI-1,
		// fsai design v0)
		//
		// policy_fsai_with_info: non-virtual outer, NVI pattern (finalize
		// guarantee + squareness entry check -- NON-throwing: a non-square
		// input is reported as invalid_input, like the ainv outer).  Must
		// never be overridden; override policy_fsai_with_info_impl instead.
		// Static FSAI after [JFSG15] (Janna et al., ACM TOMS 41(2), 2015),
		// computed from the entries of A only (no Cholesky factorization is
		// performed).  Outputs the SQRT-FREE triple: U unit upper triangular
		// (explicit unit diagonal, = Ghat^T), D diagonal (no structural
		// zeros barring overflow), perm new->old with P(perm[k],k) = 1
		// (chol/ldl orientation; ordering enum reused from chol, F-D8), all
		// born-finalized.  R = P U D^{-1} U^T P^T is NEVER materialized;
		// the estimate / apply methods below delegate to the AINV helpers
		// with Z = W = U (F-D10).  Numerical events (tiny pivots, non-SPD /
		// singular input) are lifted silently and never a failure status
		// (F-D3/F-D4).  API is _with_info only -- no strict sugar.
		// Definitions in spmats_base/spmats_fsai_impl.hpp.
		// ------------------------------------------------------------------
		fsai_result<_T,_Index> policy_fsai_with_info(
			spmats<_T,_Index>& U, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const fsai_options<_T>& opt) const;
		virtual fsai_result<_T,_Index> policy_fsai_with_info_impl(
			spmats<_T,_Index>& U, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const fsai_options<_T>& opt) const;

		// policy_fsai_residual_norm_estimate: NON-GUARANTEED estimate of
		// ||I - R A||_inf, delegated to policy_ainv_residual_norm_estimate
		// with Z = W = U against A (natural) or a permuted copy P^T A P
		// (the estimate is permutation invariant).  Single non-virtual
		// method (design SS2.3); consumer-type checks inline, non-throwing.
		fsai_status policy_fsai_residual_norm_estimate(
			const spmats<_T,_Index>& U, const spmats<_T,_Index>& D,
			const std::vector<_Index>& perm,
			typename vcp::tsparse_scalar::real_type<_T>::type& est) const;

		// policy_fsai_apply: z = R r = P (U D^{-1} U^T) P^T r, delegated to
		// policy_ainv_apply with Z = W = U plus O(n) permutation copies.
		// Single non-virtual method; checks inline, non-throwing.
		fsai_status policy_fsai_apply(
			const spmats<_T,_Index>& U, const spmats<_T,_Index>& D,
			const std::vector<_Index>& perm,
			const std::vector<_T>& r, std::vector<_T>& z) const;

		// ------------------------------------------------------------------
		// Policy methods: ADAPTIVE FSAI factored approximate inverse
		// (FSAI-2, fsai_adaptive design v0) -- [JFSG15] Algorithm 3,
		// per-row adaptive pattern generation.  Pure addition over FSAI-1:
		// the static path above is untouched.
		//
		// policy_fsai_adaptive_with_info: non-virtual outers, NVI pattern
		// (finalize guarantee + entry checks -- NON-throwing, invalid_input
		// like the fsai outer).  Must never be overridden; override the
		// _impl overloads instead.  Outputs the SAME sqrt-free triple as
		// FSAI-1 (F2-D6): U unit upper triangular (explicit unit diagonal),
		// D diagonal with D_ii = psi_k,i (design SS3.3 derivation), perm
		// new->old with P(perm[k],k) = 1, all born-finalized.
		// R = P U D^{-1} U^T P^T is NEVER materialized; use the EXISTING
		// policy_fsai_residual_norm_estimate / policy_fsai_apply on the
		// returned triple (no adaptive-specific helpers).  Numerical events
		// are lifted silently and never a failure status (F-D3/F-D4).
		// The second overload starts from an initial factor (U0, perm0)
		// ([JFSG15] Table VII combined strategy, F2-D3): U0 must be unit
		// upper triangular in the same permuted frame as perm0;
		// opt.ordering is ignored and perm0 is reused (design SS3.6).
		// Definitions in spmats_base/spmats_fsai_adaptive_impl.hpp.
		// ------------------------------------------------------------------
		fsai_adaptive_result<_T,_Index> policy_fsai_adaptive_with_info(
			spmats<_T,_Index>& U, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const fsai_adaptive_options<_T>& opt) const;
		virtual fsai_adaptive_result<_T,_Index> policy_fsai_adaptive_with_info_impl(
			spmats<_T,_Index>& U, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const fsai_adaptive_options<_T>& opt) const;

		fsai_adaptive_result<_T,_Index> policy_fsai_adaptive_with_info(
			spmats<_T,_Index>& U, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const spmats<_T,_Index>& U0, const std::vector<_Index>& perm0,
			const fsai_adaptive_options<_T>& opt) const;
		virtual fsai_adaptive_result<_T,_Index> policy_fsai_adaptive_with_info_impl(
			spmats<_T,_Index>& U, spmats<_T,_Index>& D,
			std::vector<_Index>& perm,
			const spmats<_T,_Index>& U0, const std::vector<_Index>& perm0,
			const fsai_adaptive_options<_T>& opt) const;
	};
}

// Include policy method implementations (out-of-line definitions)
#include <vcp/spmats_base/spmats_product.hpp>
#include <vcp/spmats_base/spmats_lss.hpp>
#include <vcp/spmats_base/spmats_eigs.hpp>
#include <vcp/spmats_base/spmats_ldl_impl.hpp>
#include <vcp/spmats_base/spmats_ldl_shift.hpp>
#include <vcp/spmats_base/spmats_chol_impl.hpp>
#include <vcp/spmats_base/spmats_lu_extract_impl.hpp>
#include <vcp/spmats_base/spmats_ainv_impl.hpp>
#include <vcp/spmats_base/spmats_fsai_impl.hpp>
#include <vcp/spmats_base/spmats_fsai_adaptive_impl.hpp>

#endif
