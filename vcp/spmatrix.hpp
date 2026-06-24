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

namespace vcp {


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

		spmatrix() : _P() {}
		spmatrix(const index_type rows, const index_type cols) : _P() { this->resize(rows, cols); }
		~spmatrix() = default;
		spmatrix(const spmatrix&) = default;
		spmatrix(spmatrix&&) = default;
		spmatrix& operator=(const spmatrix&) = default;
		spmatrix& operator=(spmatrix&&) = default;

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

		void sort_coo() { _P::sort_coo(); }
		void normalize_coo() { _P::normalize_coo(); }
		void finalize() { _P::finalize(); }
		void to_csr() { _P::to_csr(); }
		void to_csc() { _P::to_csc(); }

		spmatrix as_csr() const {
			spmatrix B;
			static_cast<_P&>(B) = _P::as_csr();
			return B;
		}

		spmatrix as_csc() const {
			spmatrix B;
			static_cast<_P&>(B) = _P::as_csc();
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
		void mulsm(const _T& alpha) {
			static_cast<_P&>(*this) = this->policy_scalar_mul(
				static_cast<const _P&>(*this), alpha);
		}

		// *this = *this * alpha  (commutativity: same as mulsm)
		void mulms(const _T& alpha) { mulsm(alpha); }

		// *this = *this / alpha
		void divms(const _T& alpha) {
			static_cast<_P&>(*this) = this->policy_scalar_div(
				static_cast<const _P&>(*this), alpha);
		}

		// *this = -*this
		void minusm() {
			static_cast<_P&>(*this) = this->policy_neg(
				static_cast<const _P&>(*this));
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

		// --- scalar * A, A * scalar ---
		friend spmatrix operator*(const _T& alpha, const spmatrix& A) {
			spmatrix C = A;
			C.mulsm(alpha);
			return C;
		}
		friend spmatrix operator*(const _T& alpha, spmatrix&& A) {
			A.mulsm(alpha);
			return std::move(A);
		}
		friend spmatrix operator*(const spmatrix& A, const _T& alpha) {
			spmatrix C = A;
			C.mulsm(alpha);
			return C;
		}
		friend spmatrix operator*(spmatrix&& A, const _T& alpha) {
			A.mulsm(alpha);
			return std::move(A);
		}

		// --- A / scalar ---
		friend spmatrix operator/(const spmatrix& A, const _T& alpha) {
			spmatrix C = A;
			C.divms(alpha);
			return C;
		}
		friend spmatrix operator/(spmatrix&& A, const _T& alpha) {
			A.divms(alpha);
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
		friend spmatrix& operator*=(spmatrix& A, const _T& alpha) {
			A.mulsm(alpha);
			return A;
		}
		friend spmatrix& operator/=(spmatrix& A, const _T& alpha) {
			A.divms(alpha);
			return A;
		}

		// ---------------------------------------------------------------
		// Linear system solve — all delegate to policy_lss / policy_lss_with_info
		// ---------------------------------------------------------------

		// strict solve: policy decides convergence checking
		std::vector<_T> solve(const std::vector<_T>& b, const linear_solve_options_type& options = linear_solve_options_type()) const {
			return this->policy_lss(static_cast<const _P&>(*this), b, options);
		}

		// non-strict: return full diagnostic result
		linear_solve_result<_T> solve_with_info(const std::vector<_T>& b, const linear_solve_options_type& options = linear_solve_options_type()) const {
			return this->policy_lss_with_info(static_cast<const _P&>(*this), b, options);
		}

		// Convenience overloads — build options and delegate to solve / solve_with_info
		std::vector<_T> solve_jacobi(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::jacobi;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = true;
			return this->policy_lss(static_cast<const _P&>(*this), b, opt_);
		}

		std::vector<_T> solve_gauss_seidel(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::gauss_seidel;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = true;
			return this->policy_lss(static_cast<const _P&>(*this), b, opt_);
		}

		std::vector<_T> solve_cg(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::conjugate_gradient;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.check_symmetric = true;
			opt_.preconditioner = preconditioner_type::none;
			opt_.use_relative_residual = true;
			return this->policy_lss(static_cast<const _P&>(*this), b, opt_);
		}

		std::vector<_T> solve_bicgstab(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::bicgstab;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.use_relative_residual = true;
			return this->policy_lss(static_cast<const _P&>(*this), b, opt_);
		}

		std::vector<_T> solve_gmres(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			linear_solve_options_type opt_;
			opt_.method = linear_solver_method::gmres;
			opt_.max_iter = max_iter;
			opt_.tol = tol;
			opt_.restart = 30;
			opt_.use_relative_residual = true;
			return this->policy_lss(static_cast<const _P&>(*this), b, opt_);
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
			return this->policy_eig(static_cast<const _P&>(*this), options);
		}

		// Full dense eig (non-strict / diagnostic)
		eig_result<_T> eig_with_info(const eig_options_type& options) const {
			validate_eig_input("spmatrix::eig");
			if (options.max_iter == 0 || options.tol <= scalar_real_type(0))
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eig: invalid iteration option");
			if (options.method != eig_solver_method::dense_fallback_explicit)
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eig: full dense eig requires dense_fallback_explicit");
			return this->policy_eigs_with_info(
				static_cast<const _P&>(*this), static_cast<std::size_t>(rowsize()), options);
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
			return this->policy_eigs(static_cast<const _P&>(*this), k, options);
		}

		// eigs_with_info(k): shorthand
		eig_result<_T> eigs_with_info(const std::size_t k) const {
			return eigs_with_info(k, default_eigs_options());
		}

		// eigs_with_info(k, opt): non-strict, return full result
		eig_result<_T> eigs_with_info(const std::size_t k, const eig_options_type& options) const {
			return this->policy_eigs_with_info(
				static_cast<const _P&>(*this), k, options);
		}

		// Full generalized eig (strict, returns eigenvalue vector)
		std::vector<_T> eig(const spmatrix& B) const {
			validate_eig_input("spmatrix::eig(A,B)");
			return this->policy_generalized_eigs(
				static_cast<const _P&>(*this), static_cast<const _P&>(B),
				static_cast<std::size_t>(rowsize()), eig_options_type());
		}

		// Full generalized eig (strict, returns eig_result)
		eig_result<_T> eig(const spmatrix& B, const eig_options_type& options) const {
			validate_eig_input("spmatrix::eig(A,B)");
			return this->policy_generalized_eig(
				static_cast<const _P&>(*this), static_cast<const _P&>(B),
				static_cast<std::size_t>(rowsize()), options);
		}

		// Partial generalized eigs (strict): policy_generalized_eigs decides convergence
		std::vector<_T> eigs(const spmatrix& B, const std::size_t k, const eig_options_type& options = eig_options_type()) const {
			return this->policy_generalized_eigs(
				static_cast<const _P&>(*this), static_cast<const _P&>(B), k, options);
		}

		// Partial generalized eigs_with_info (non-strict)
		eig_result<_T> eigs_with_info(const spmatrix& B, const std::size_t k, const eig_options_type& options = eig_options_type()) const {
			return this->policy_generalized_eigs_with_info(
				static_cast<const _P&>(*this), static_cast<const _P&>(B), k, options);
		}

		// ------------------------------------------------------------------
		// Phase 6: preconditioner overloads — all delegate to policy
		// ------------------------------------------------------------------

		template <class Preconditioner>
		std::vector<_T> eigs(const std::size_t k, const eig_options_type& options,
		                     const Preconditioner& M) const {
			return this->policy_eigs(static_cast<const _P&>(*this), k, options, M);
		}

		template <class Preconditioner>
		eig_result<_T> eigs_with_info(const std::size_t k, const eig_options_type& options,
		                              const Preconditioner& M) const {
			return this->policy_eigs_with_info(
				static_cast<const _P&>(*this), k, options, M);
		}

		template <class Preconditioner>
		std::vector<_T> eigs(const spmatrix& B, const std::size_t k,
		                     const eig_options_type& options, const Preconditioner& M) const {
			return this->policy_generalized_eigs(
				static_cast<const _P&>(*this), static_cast<const _P&>(B), k, options, M);
		}

		template <class Preconditioner>
		eig_result<_T> eigs_with_info(const spmatrix& B, const std::size_t k,
		                              const eig_options_type& options,
		                              const Preconditioner& M) const {
			return this->policy_generalized_eigs_with_info(
				static_cast<const _P&>(*this), static_cast<const _P&>(B), k, options, M);
		}

		spmatrix transpose() const {
			spmatrix B;
			static_cast<_P&>(B) = _P::transpose();
			return B;
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

	private:
		void validate_eig_input(const char* routine) const {
			if (rowsize() != columnsize()) vcp::throw_error<vcp::dimension_error>(routine, ": matrix must be square");
		}
	};


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
