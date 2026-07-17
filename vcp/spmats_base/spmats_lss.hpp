// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spmats_lss.hpp
// Policy method implementations for linear system solving on spmats<_T,_Index>.
// This file is included inside the namespace vcp {} block, AFTER the closing
// brace of spmats<_T,_Index>, via spmats.hpp.

#ifndef VCP_SPMATS_LSS_HPP
#define VCP_SPMATS_LSS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_solvers.hpp>
#include <vcp/spmats_base/spmats_eigs_types.hpp>

namespace vcp {

// ---------------------------------------------------------------------------
// Private helper methods (implemented as free functions operating on spmats)
// ---------------------------------------------------------------------------

namespace spmats_lss_detail {

	template <typename _T, typename _Index>
	static typename vcp::tsparse_scalar::real_type<_T>::type
	dot_value_lss(const std::vector<_T>& a, const std::vector<_T>& b) {
		return vcp::tsparse_scalar::real_dot_value(a, b);
	}

	template <typename _T, typename _Index>
	static typename vcp::tsparse_scalar::real_type<_T>::type
	norm_value_lss(const std::vector<_T>& a) {
		return vcp::tsparse_scalar::real_norm_value(a);
	}

	template <typename _T, typename _Index>
	static typename vcp::tsparse_scalar::real_type<_T>::type
	residual_norm_lss(const spmats<_T,_Index>& A,
	                  const std::vector<_T>& x,
	                  const std::vector<_T>& b) {
		return vcp::tsparse_solvers::residual_norm_value(A, x, b);
	}

	template <typename _T, typename _Index>
	static void set_linear_residual_lss(linear_solve_result<_T>& result,
	                                    const spmats<_T,_Index>& A,
	                                    const std::vector<_T>& b) {
		vcp::tsparse_solvers::set_linear_residual_fields(result, A, b);
	}

	template <typename _T, typename _Index>
	static std::vector<_T> make_jacobi_inv_diag_lss(const spmats<_T,_Index>& A,
	                                                  const preconditioner_type prec,
	                                                  const char* routine) {
		std::vector<_T> inv_diag;
		if (prec == preconditioner_type::none) return inv_diag;
		if (prec != preconditioner_type::jacobi) {
			vcp::throw_error<vcp::state_error>(routine, ": unknown preconditioner");
			return inv_diag;
		}
		const std::size_t n = static_cast<std::size_t>(A.rowsize());
		inv_diag.assign(n, _T(0));
		for (std::size_t i = 0; i < n; i++) {
			const _T diag = A.get(static_cast<_Index>(i), static_cast<_Index>(i));
			if (diag == _T(0)) {
				vcp::throw_error<vcp::numerical_error>(routine, ": Jacobi preconditioner has zero diagonal");
			}
			inv_diag[i] = _T(1) / diag;
		}
		return inv_diag;
	}

	template <typename _T>
	static std::vector<_T> apply_left_prec_lss(const std::vector<_T>& r,
	                                             const std::vector<_T>& inv_diag,
	                                             const preconditioner_type prec) {
		if (prec == preconditioner_type::none) return r;
		std::vector<_T> z(r.size(), _T(0));
		for (std::size_t i = 0; i < r.size(); i++) z[i] = inv_diag[i] * r[i];
		return z;
	}

	template <typename _T, typename _Index>
	static std::vector<_T> apply_prec_op_lss(const spmats<_T,_Index>& A,
	                                           const std::vector<_T>& x,
	                                           const std::vector<_T>& inv_diag,
	                                           const preconditioner_type prec) {
		const std::vector<_T> ax = A.mul_vec(x);
		return apply_left_prec_lss(ax, inv_diag, prec);
	}

	template <typename _T, typename _Index>
	static std::vector<_T> solve_upper_tri_lss(const std::vector<std::vector<_T> >& R,
	                                             const std::vector<_T>& rhs,
	                                             const std::size_t n) {
		typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
		std::vector<_T> y(n, _T(0));
		for (std::size_t i = n; i-- > 0; ) {
			_T sum = rhs[i];
			for (std::size_t j = i + 1; j < n; j++) sum -= R[i][j] * y[j];
			const real_type denom = vcp::tsparse_scalar::abs_value(R[i][i]);
			if (!(denom > vcp::tsparse_scalar::epsilon<real_type>())) { y[i] = _T(0); continue; }
			y[i] = sum / R[i][i];
		}
		return y;
	}

	// ---------------------------------------------------------------------------
	// dispatch_sparse_lu_: SFINAE-guarded helper to avoid static_assert
	// instantiation for unsigned Index when sparse_lu case is compiled.
	// ---------------------------------------------------------------------------

	// signed Index path: calls sparse_lu_factorize_with_info
	template <typename _T, typename _Index>
	inline typename std::enable_if<std::is_signed<_Index>::value,
	                               linear_solve_result<_T> >::type
	dispatch_sparse_lu_(
	    const spmats<_T,_Index>& A,
	    const std::vector<_T>& b,
	    const linear_solve_options<_T>& opt)
	{
	    typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	    const std::size_t n = static_cast<std::size_t>(A.rowsize());
	    vcp::sparse_lu_factorization<_T, _Index> fac =
	        vcp::sparse_lu_factorize_with_info(A, opt.sparse_lu);

	    linear_solve_result<_T> result;
	    result.method = linear_solver_method::sparse_lu;
	    result.iterations = 0;

	    if (fac.info().success) {
	        if (opt.sparse_lu.iterative_refinement) {
	            // IR path: sparse_lu_solve_refined runs fac.solve(b) then refines.
	            // A is in scope here so residual r = b - A*x can be computed.
	            vcp::sparse_lu_refinement_info<_T> ir_info;
	            result.x = vcp::sparse_lu_solve_refined(A, fac, b, opt.sparse_lu, &ir_info);
	            result.iterations             = ir_info.iterations;
	            result.initial_residual_norm  = ir_info.initial_residual;
	            result.residual_norm          = ir_info.final_residual;
	            result.absolute_residual_norm = ir_info.final_residual;
	            result.relative_residual_norm = ir_info.final_relative_residual;
	        } else {
	            // Plain solve path (LSS-1 P-6): x is unchanged (byte-identical to
	            // the pre-P-6 path); the residual fields now report the ACTUAL
	            // residual of the returned x (one extra matvec, report-only)
	            // instead of the former fabricated real_type(0).
	            // initial_residual_norm = absolute_residual_norm is the
	            // 0-iteration reading, consistent with the IR path above.
	            // converged stays fac.info().success (D-7: unchanged).
	            result.x = fac.solve(b);
	            set_linear_residual_lss(result, A, b);
	            result.initial_residual_norm = result.absolute_residual_norm;
	        }
	        result.solution  = result.x;
	        result.converged = true;
	    } else {
	        result.converged = false;
	        result.x.assign(n, _T(0));
	        result.solution = result.x;
	        // SLU-GT1 D5: residual fields remain at the default real_type(0);
	        // converged == false marks them as undefined (do not read).
	    }
	    return result;
	}

	// unsigned Index path: sparse_lu cannot be used (Index must be signed)
	template <typename _T, typename _Index>
	inline typename std::enable_if<!std::is_signed<_Index>::value,
	                               linear_solve_result<_T> >::type
	dispatch_sparse_lu_(
	    const spmats<_T,_Index>& A,
	    const std::vector<_T>& b,
	    const linear_solve_options<_T>& opt)
	{
	    (void)A; (void)b; (void)opt;
	    vcp::throw_error<vcp::state_error>(
	        "spmats::policy_lss_with_info: sparse_lu requires a signed Index type");
	    return linear_solve_result<_T>();
	}

	// ---------------------------------------------------------------------------
	// resolve_auto_nonsymmetric_: SFINAE-guarded nonsymmetric branch of the
	// LSS-1 P-1 auto_select resolution (same split pattern as
	// dispatch_sparse_lu_ / dispatch_sparse_lu_extract_ -- the signed-only
	// body is never instantiated for unsigned Index).
	// ---------------------------------------------------------------------------

	// signed Index path: direct method.  D-3: amd is set ONLY when the user
	// left opt.sparse_lu.ordering == auto_select; the sparse_lu_options
	// default itself is unchanged (no propagation to the eig shift-invert LU).
	template <typename _T, typename _Index>
	inline typename std::enable_if<std::is_signed<_Index>::value, void>::type
	resolve_auto_nonsymmetric_(linear_solve_options<_T>& resolved)
	{
	    resolved.method = linear_solver_method::sparse_lu;
	    // SLU-L1 L-1 makes this redundant (auto_select now resolves to amd inside
	    // sparse_lu_symbolic); RETAINED as defense for the auto path (design §4.4).
	    if (resolved.sparse_lu.ordering == sparse_lu_ordering::auto_select)
	        resolved.sparse_lu.ordering = sparse_lu_ordering::amd;   // D-3
	}

	// unsigned Index path: sparse_lu is impossible; no silent fallback (D-4).
	template <typename _T, typename _Index>
	inline typename std::enable_if<!std::is_signed<_Index>::value, void>::type
	resolve_auto_nonsymmetric_(linear_solve_options<_T>& resolved)
	{
	    (void)resolved;
	    vcp::throw_error<vcp::state_error>(
	        "spmats::policy_lss_with_info: auto_select cannot solve a nonsymmetric system with unsigned Index (sparse_lu requires a signed Index type); specify an iterative method explicitly");
	}

} // namespace spmats_lss_detail

// ---------------------------------------------------------------------------
// solve_jacobi_with_info (private policy helper)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_solve_jacobi_with_info_(
	const spmats<_T, _Index>& A_in,
	const std::vector<_T>& b,
	const std::size_t max_iter,
	const scalar_real_type& tol,
	const bool use_relative_residual) const
{
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	const _Index n = A_in.rowsize();
	std::vector<_T> diag(static_cast<std::size_t>(n), _T(0));
	for (_Index i = 0; i < n; i++) {
		diag[static_cast<std::size_t>(i)] = A_in.get(i, i);
		if (diag[static_cast<std::size_t>(i)] == _T(0)) {
			vcp::throw_error<vcp::numerical_error>("spmats::solve_jacobi: zero diagonal");
		}
	}
	spmats<_T,_Index> A = A_in.as_csr();
	std::vector<_T> x(static_cast<std::size_t>(n), _T(0));
	std::vector<_T> next(static_cast<std::size_t>(n), _T(0));
	linear_solve_result<_T> result;
	result.converged = false;
	result.iterations = 0;
	result.method = linear_solver_method::jacobi;
	const vcp::tsparse_solvers::residual_control<_T> control =
		vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
	real_type current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
	result.residual_norm = current_residual;
	result.converged = current_residual <= control.threshold;
	for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
		const std::vector<_Index>& outer = A.outer_index();
		const std::vector<_Index>& inner = A.inner_index();
		const std::vector<_T>& val = A.values();
		for (_Index i = 0; i < n; i++) {
			_T sigma = _T(0);
			for (_Index p = outer[static_cast<std::size_t>(i)];
			     p < outer[static_cast<std::size_t>(i + 1)]; p++) {
				const _Index j = inner[static_cast<std::size_t>(p)];
				if (j != i) sigma += val[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(j)];
			}
			next[static_cast<std::size_t>(i)] =
				(b[static_cast<std::size_t>(i)] - sigma) / diag[static_cast<std::size_t>(i)];
		}
		x.swap(next);
		current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
		result.residual_norm = current_residual;
		result.iterations = iter;
		if (current_residual <= control.threshold) { result.converged = true; break; }
	}
	result.x = x;
	spmats_lss_detail::set_linear_residual_lss(result, A, b);
	return result;
}

// ---------------------------------------------------------------------------
// solve_gauss_seidel_with_info (private policy helper)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_solve_gauss_seidel_with_info_(
	const spmats<_T, _Index>& A_in,
	const std::vector<_T>& b,
	const std::size_t max_iter,
	const scalar_real_type& tol,
	const bool use_relative_residual) const
{
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	const _Index n = A_in.rowsize();
	for (_Index i = 0; i < n; i++) {
		if (A_in.get(i, i) == _T(0)) {
			vcp::throw_error<vcp::numerical_error>("spmats::solve_gauss_seidel: zero diagonal");
		}
	}
	spmats<_T,_Index> A = A_in.as_csr();
	const std::vector<_Index>& outer = A.outer_index();
	const std::vector<_Index>& inner = A.inner_index();
	const std::vector<_T>& val = A.values();
	std::vector<_T> x(static_cast<std::size_t>(n), _T(0));
	linear_solve_result<_T> result;
	result.converged = false;
	result.iterations = 0;
	result.method = linear_solver_method::gauss_seidel;
	const vcp::tsparse_solvers::residual_control<_T> control =
		vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
	real_type current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
	result.residual_norm = current_residual;
	result.converged = current_residual <= control.threshold;
	for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
		for (_Index i = 0; i < n; i++) {
			_T sigma = _T(0);
			_T diag = _T(0);
			for (_Index p = outer[static_cast<std::size_t>(i)];
			     p < outer[static_cast<std::size_t>(i + 1)]; p++) {
				const _Index j = inner[static_cast<std::size_t>(p)];
				if (j == i) diag = val[static_cast<std::size_t>(p)];
				else sigma += val[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(j)];
			}
			x[static_cast<std::size_t>(i)] = (b[static_cast<std::size_t>(i)] - sigma) / diag;
		}
		current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
		result.residual_norm = current_residual;
		result.iterations = iter;
		if (current_residual <= control.threshold) { result.converged = true; break; }
	}
	result.x = x;
	spmats_lss_detail::set_linear_residual_lss(result, A, b);
	return result;
}

// ---------------------------------------------------------------------------
// solve_cg_with_info (private policy helper)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_solve_cg_with_info_(
	const spmats<_T, _Index>& A_in,
	const std::vector<_T>& b,
	const std::size_t max_iter,
	const scalar_real_type& tol,
	const bool check_symmetric,
	const preconditioner_type prec_type,
	const bool use_relative_residual) const
{
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	if (check_symmetric) {
		const real_type sym_tol =
			vcp::tsparse_scalar::decimal_power_negative<real_type>(10);
		if (!A_in.is_symmetric(sym_tol)) {
			vcp::throw_error<vcp::domain_error>("spmats::solve_cg: matrix must be symmetric");
		}
	}
	spmats<_T,_Index> A = A_in.as_csr();
	const std::size_t n = b.size();
	std::vector<_T> x(n, _T(0));
	std::vector<_T> r = b;
	const std::vector<_T> inv_diag =
		spmats_lss_detail::make_jacobi_inv_diag_lss<_T,_Index>(A, prec_type, "spmats::solve_cg");
	std::vector<_T> z = spmats_lss_detail::apply_left_prec_lss(r, inv_diag, prec_type);
	std::vector<_T> p = z;
	real_type rzold = spmats_lss_detail::dot_value_lss<_T,_Index>(r, z);
	const vcp::tsparse_solvers::residual_control<_T> control =
		vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
	linear_solve_result<_T> result;
	real_type current_residual = spmats_lss_detail::norm_value_lss<_T,_Index>(r);
	result.converged = current_residual <= control.threshold;
	result.iterations = 0;
	result.residual_norm = current_residual;
	result.method = linear_solver_method::conjugate_gradient;
	for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
		const std::vector<_T> Ap = A.mul_vec(p);
		const real_type denom = spmats_lss_detail::dot_value_lss<_T,_Index>(p, Ap);
		if (!(denom > real_type(0)) || !vcp::tsparse_scalar::is_finite(denom)) {
			result.converged = false; break;
		}
		const _T alpha = _T(rzold / denom);
		for (std::size_t i = 0; i < n; i++) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
		current_residual = spmats_lss_detail::norm_value_lss<_T,_Index>(r);
		result.residual_norm = current_residual;
		result.iterations = iter;
		if (current_residual <= control.threshold) { result.converged = true; break; }
		z = spmats_lss_detail::apply_left_prec_lss(r, inv_diag, prec_type);
		const real_type rznew = spmats_lss_detail::dot_value_lss<_T,_Index>(r, z);
		if (!(rznew > real_type(0)) || !vcp::tsparse_scalar::is_finite(rznew)) {
			result.converged = false; break;
		}
		const _T beta = _T(rznew / rzold);
		for (std::size_t i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
		rzold = rznew;
	}
	result.x = x;
	spmats_lss_detail::set_linear_residual_lss(result, A, b);
	return result;
}

// ---------------------------------------------------------------------------
// solve_bicgstab_with_info (private policy helper)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_solve_bicgstab_with_info_(
	const spmats<_T, _Index>& A_in,
	const std::vector<_T>& b,
	const std::size_t max_iter,
	const scalar_real_type& tol,
	const bool use_relative_residual,
	const preconditioner_type prec_type) const
{
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	spmats<_T,_Index> A = A_in.as_csr();
	const std::size_t n = b.size();
	const std::vector<_T> inv_diag =
		spmats_lss_detail::make_jacobi_inv_diag_lss<_T,_Index>(A, prec_type, "spmats::solve_bicgstab");
	const std::vector<_T> pb =
		spmats_lss_detail::apply_left_prec_lss(b, inv_diag, prec_type);
	const vcp::tsparse_solvers::residual_control<_T> control =
		vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
	std::vector<_T> x(n, _T(0));
	std::vector<_T> r = pb;
	std::vector<_T> r_hat = r;
	std::vector<_T> p(n, _T(0));
	std::vector<_T> v(n, _T(0));
	real_type rho_old = real_type(1);
	real_type alpha = real_type(1);
	real_type omega = real_type(1);
	linear_solve_result<_T> result;
	result.x = x;
	real_type current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
	result.residual_norm = current_residual;
	result.converged = current_residual <= control.threshold;
	result.iterations = 0;
	result.method = linear_solver_method::bicgstab;
	const real_type eps = vcp::tsparse_scalar::epsilon<real_type>();
	for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
		const real_type rho = spmats_lss_detail::dot_value_lss<_T,_Index>(r_hat, r);
		if (!(vcp::tsparse_scalar::abs_value(rho) > eps) || !vcp::tsparse_scalar::is_finite(rho)) break;
		const real_type beta = (rho / rho_old) * (alpha / omega);
		for (std::size_t i = 0; i < n; i++)
			p[i] = r[i] + _T(beta) * (p[i] - _T(omega) * v[i]);
		v = spmats_lss_detail::apply_prec_op_lss<_T,_Index>(A, p, inv_diag, prec_type);
		const real_type denom = spmats_lss_detail::dot_value_lss<_T,_Index>(r_hat, v);
		if (!(vcp::tsparse_scalar::abs_value(denom) > eps) || !vcp::tsparse_scalar::is_finite(denom)) break;
		alpha = rho / denom;
		std::vector<_T> s(n);
		for (std::size_t i = 0; i < n; i++) s[i] = r[i] - _T(alpha) * v[i];
		const real_type snorm = spmats_lss_detail::norm_value_lss<_T,_Index>(s);
		if (snorm <= control.threshold) {
			for (std::size_t i = 0; i < n; i++) x[i] += _T(alpha) * p[i];
			current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
			result.residual_norm = current_residual;
			result.iterations = iter;
			result.converged = current_residual <= control.threshold;
			break;
		}
		const std::vector<_T> t = spmats_lss_detail::apply_prec_op_lss<_T,_Index>(A, s, inv_diag, prec_type);
		const real_type tt = spmats_lss_detail::dot_value_lss<_T,_Index>(t, t);
		if (!(tt > eps) || !vcp::tsparse_scalar::is_finite(tt)) break;
		omega = spmats_lss_detail::dot_value_lss<_T,_Index>(t, s) / tt;
		if (!(vcp::tsparse_scalar::abs_value(omega) > eps) || !vcp::tsparse_scalar::is_finite(omega)) break;
		for (std::size_t i = 0; i < n; i++) {
			x[i] += _T(alpha) * p[i] + _T(omega) * s[i];
			r[i] = s[i] - _T(omega) * t[i];
		}
		current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
		result.residual_norm = current_residual;
		result.iterations = iter;
		result.converged = current_residual <= control.threshold;
		rho_old = rho;
	}
	result.x = x;
	spmats_lss_detail::set_linear_residual_lss(result, A, b);
	return result;
}

// ---------------------------------------------------------------------------
// solve_gmres_with_info (private policy helper)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_solve_gmres_with_info_(
	const spmats<_T, _Index>& A_in,
	const std::vector<_T>& b,
	const std::size_t max_iter,
	const scalar_real_type& tol,
	const std::size_t restart,
	const bool use_relative_residual,
	const preconditioner_type prec_type) const
{
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	if (restart == 0) vcp::throw_error<vcp::invalid_argument>("spmats::solve_gmres: restart must be positive");
	spmats<_T,_Index> A = A_in.as_csr();
	const std::size_t n = b.size();
	const std::vector<_T> inv_diag =
		spmats_lss_detail::make_jacobi_inv_diag_lss<_T,_Index>(A, prec_type, "spmats::solve_gmres");
	const std::vector<_T> pb =
		spmats_lss_detail::apply_left_prec_lss(b, inv_diag, prec_type);
	const vcp::tsparse_solvers::residual_control<_T> control =
		vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
	std::vector<_T> x(n, _T(0));
	linear_solve_result<_T> result;
	result.x = x;
	result.converged = false;
	result.iterations = 0;
	real_type current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
	result.residual_norm = current_residual;
	result.method = linear_solver_method::gmres;
	if (current_residual <= control.threshold) {
		result.converged = true;
		spmats_lss_detail::set_linear_residual_lss(result, A, b);
		return result;
	}
	const std::size_t mmax = std::min(restart, n);
	while (result.iterations < max_iter && !result.converged) {
		std::vector<_T> Ax = A.mul_vec(x);
		const std::vector<_T> pAx = spmats_lss_detail::apply_left_prec_lss(Ax, inv_diag, prec_type);
		std::vector<_T> r(n);
		for (std::size_t i = 0; i < n; i++) r[i] = pb[i] - pAx[i];
		const real_type beta = spmats_lss_detail::norm_value_lss<_T,_Index>(r);
		current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
		result.residual_norm = current_residual;
		if (current_residual <= control.threshold) { result.converged = true; break; }
		std::vector<std::vector<_T> > V(mmax + 1, std::vector<_T>(n, _T(0)));
		std::vector<std::vector<_T> > H(mmax + 1, std::vector<_T>(mmax, _T(0)));
		std::vector<real_type> cs(mmax, real_type(0));
		std::vector<real_type> sn(mmax, real_type(0));
		std::vector<_T> g(mmax + 1, _T(0));
		g[0] = _T(beta);
		for (std::size_t i = 0; i < n; i++) V[0][i] = r[i] / _T(beta);
		std::size_t m = 0;
		for (; m < mmax && result.iterations < max_iter; m++) {
			std::vector<_T> w = spmats_lss_detail::apply_prec_op_lss<_T,_Index>(A, V[m], inv_diag, prec_type);
			for (std::size_t j = 0; j <= m; j++) {
				H[j][m] = _T(spmats_lss_detail::dot_value_lss<_T,_Index>(w, V[j]));
				for (std::size_t i = 0; i < n; i++) w[i] -= H[j][m] * V[j][i];
			}
			const real_type hnext = spmats_lss_detail::norm_value_lss<_T,_Index>(w);
			H[m + 1][m] = _T(hnext);
			if (hnext > vcp::tsparse_scalar::epsilon<real_type>()) {
				for (std::size_t i = 0; i < n; i++) V[m + 1][i] = w[i] / _T(hnext);
			}
			for (std::size_t j = 0; j < m; j++) {
				const _T hij = H[j][m];
				const _T hip1j = H[j + 1][m];
				H[j][m] = _T(cs[j]) * hij + _T(sn[j]) * hip1j;
				H[j + 1][m] = _T(-sn[j]) * hij + _T(cs[j]) * hip1j;
			}
			const real_type h0 = vcp::tsparse_scalar::abs_value(H[m][m]);
			const real_type h1 = vcp::tsparse_scalar::abs_value(H[m + 1][m]);
			const real_type rho = vcp::tsparse_scalar::hypot_value(h0, h1);
			if (!(rho > vcp::tsparse_scalar::epsilon<real_type>())) {
				cs[m] = real_type(1); sn[m] = real_type(0);
			} else {
				cs[m] = vcp::tsparse_scalar::real_part(H[m][m]) / rho;
				sn[m] = vcp::tsparse_scalar::real_part(H[m + 1][m]) / rho;
			}
			const _T hmm = H[m][m];
			const _T hp1m = H[m + 1][m];
			H[m][m] = _T(cs[m]) * hmm + _T(sn[m]) * hp1m;
			H[m + 1][m] = _T(0);
			const _T gm = g[m];
			g[m] = _T(cs[m]) * gm;
			g[m + 1] = _T(-sn[m]) * gm;
			result.iterations++;
			const real_type projected_residual = vcp::tsparse_scalar::abs_value(g[m + 1]);
			result.residual_norm = projected_residual;
			if (projected_residual <= control.threshold) {
				m++;
				std::vector<_T> y = spmats_lss_detail::solve_upper_tri_lss<_T,_Index>(H, g, m);
				for (std::size_t j = 0; j < m; j++)
					for (std::size_t i = 0; i < n; i++) x[i] += V[j][i] * y[j];
				current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
				result.residual_norm = current_residual;
				result.converged = current_residual <= control.threshold;
				break;
			}
			if (!(hnext > vcp::tsparse_scalar::epsilon<real_type>())) { m++; break; }
		}
		if (m == 0) break;
		if (!result.converged) {
			std::vector<_T> y = spmats_lss_detail::solve_upper_tri_lss<_T,_Index>(H, g, m);
			for (std::size_t j = 0; j < m; j++)
				for (std::size_t i = 0; i < n; i++) x[i] += V[j][i] * y[j];
			current_residual = spmats_lss_detail::residual_norm_lss<_T,_Index>(A, x, b);
			result.residual_norm = current_residual;
			result.converged = current_residual <= control.threshold;
		}
	}
	result.x = x;
	spmats_lss_detail::set_linear_residual_lss(result, A, b);
	return result;
}

// ---------------------------------------------------------------------------
// policy_lss_with_info: NVI outer (non-virtual). Finalizes A, then delegates
// to the virtual policy_lss_with_info_impl. Must never be overridden —
// override policy_lss_with_info_impl instead (see
// sandbox/docs/design/spmats_finalize_policy.md §3).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_lss_with_info(
	const std::vector<_T>& b,
	const linear_solve_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	return policy_lss_with_info_impl(b, opt);
}

// ---------------------------------------------------------------------------
// policy_lss_with_info_impl: virtual algorithm body (input validation +
// method dispatch). Custom policies (SuperLU, etc.) override this, not
// policy_lss_with_info.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
linear_solve_result<_T> spmats<_T, _Index>::policy_lss_with_info_impl(
	const std::vector<_T>& b,
	const linear_solve_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	// Validate inputs (misuse contract: invalid input THROWS vcp::error;
	// the SLU-GT1 D6 net below rethrows these unchanged)
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>("spmats::policy_lss_with_info: matrix must be square");
	if (b.size() != static_cast<std::size_t>(A.rowsize()))
		vcp::throw_error<vcp::dimension_error>("spmats::policy_lss_with_info: rhs dimension mismatch");
	if (opt.max_iter == 0)
		vcp::throw_error<vcp::invalid_argument>("spmats::policy_lss_with_info: max_iter must be positive");
	if (opt.tol <= scalar_real_type(0))
		vcp::throw_error<vcp::invalid_argument>("spmats::policy_lss_with_info: tol must be positive");

	// -------------------------------------------------------------------
	// LSS-1 P-1: auto_select resolution (design §2.2).  Executed ONLY when
	// the caller left method == auto_select; every explicitly requested
	// method reaches the switch below through the unchanged path.
	//   symmetric (policy_is_symmetric, default tol 1e-12; complex-symmetric
	//   check, not Hermitian) -> conjugate_gradient (D-2; the CG-internal
	//   check_symmetric with tol 1e-10 stays active, D-5),
	//   nonsymmetric -> sparse_lu (signed Index; +amd only when the user's
	//   sparse_lu.ordering == auto_select, D-3) or vcp::state_error for
	//   unsigned Index (D-4, no silent fallback).
	// The returned result.method is the RESOLVED method (each solve helper
	// stamps its own value; auto_select is never returned).  Re-entry depth
	// is exactly 1: resolved.method != auto_select.
	// NOTE: auto is not a universal best pick — for a huge nonsymmetric
	// system sparse_lu can be expensive; choose an explicit iterative
	// method there.
	// -------------------------------------------------------------------
	if (opt.method == linear_solver_method::auto_select) {
		linear_solve_options<_T> resolved = opt;
		if (A.policy_is_symmetric(A)) {
			resolved.method = linear_solver_method::conjugate_gradient;
		} else {
			spmats_lss_detail::resolve_auto_nonsymmetric_<_T, _Index>(resolved);
		}
		return policy_lss_with_info_impl(b, resolved);
	}

	try {
		switch (opt.method) {
		case linear_solver_method::jacobi:
			return policy_solve_jacobi_with_info_(A, b, opt.max_iter, opt.tol, opt.use_relative_residual);
		case linear_solver_method::gauss_seidel:
			return policy_solve_gauss_seidel_with_info_(A, b, opt.max_iter, opt.tol, opt.use_relative_residual);
		case linear_solver_method::conjugate_gradient:
			return policy_solve_cg_with_info_(A, b, opt.max_iter, opt.tol,
			                                  opt.check_symmetric, opt.preconditioner, opt.use_relative_residual);
		case linear_solver_method::bicgstab:
			return policy_solve_bicgstab_with_info_(A, b, opt.max_iter, opt.tol,
			                                        opt.use_relative_residual, opt.preconditioner);
		case linear_solver_method::gmres:
			return policy_solve_gmres_with_info_(A, b, opt.max_iter, opt.tol,
			                                     opt.restart, opt.use_relative_residual, opt.preconditioner);
		case linear_solver_method::sparse_lu:
			return spmats_lss_detail::dispatch_sparse_lu_<_T, _Index>(A, b, opt);
		}
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		// SLU-GT1 D6: certified ゲート(D1/D3)が正しければ到達しない最終防護網。
		// 発火は「ゲートの取りこぼし」を意味する(調査対象)。
		// linear_solve_result にはステータス欄がないため converged=false のまま
		// 返す(residual フィールドは D5 契約により未定義)。
		linear_solve_result<_T> result;
		result.method = opt.method;
		result.converged = false;
		result.x.assign(b.size(), _T(0));
		result.solution = result.x;
		return result;
	}
	vcp::throw_error<vcp::state_error>("spmats::policy_lss_with_info: unknown method");
	return linear_solve_result<_T>();
}

// ---------------------------------------------------------------------------
// policy_lss: thin wrapper (no-throw version)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
std::vector<_T> spmats<_T, _Index>::policy_lss(
	const std::vector<_T>& b,
	const linear_solve_options<_T>& opt) const
{
	linear_solve_result<_T> result = policy_lss_with_info(b, opt);
	if (!result.converged)
		vcp::throw_error<vcp::state_error>("spmats::policy_lss: iterative solver did not converge");
	return result.x;
}

// ===========================================================================
// LSS-1 P-4: policy_lu_factorize_with_info (reusable LU factorization handle)
// Type: spmats_base/spmats_lu_factor.hpp.  Defined here so the handle's
// solve semantics stay next to dispatch_sparse_lu_, which they mirror.
// ===========================================================================

namespace spmats_lss_detail {

	// signed Index path: factorize once and pack the handle.  On success the
	// handle keeps a finalized CSR copy of A for the IR residual (§3.3).
	template <typename _T, typename _Index>
	inline typename std::enable_if<std::is_signed<_Index>::value,
	                               lu_factor_handle<_T,_Index> >::type
	dispatch_lu_factorize_(
	    const spmats<_T,_Index>& A,
	    const sparse_lu_options<_T>& opt)
	{
	    lu_factor_handle<_T,_Index> handle;
	    const _Index n = A.rowsize();
	    vcp::sparse_lu_factorization<_T, _Index> fac =
	        vcp::sparse_lu_factorize_with_info(A, opt);
	    const bool ok = fac.info().success;
	    spmats_lu_factor_detail::handle_access::assign(
	        handle, std::move(fac),
	        ok ? A.as_csr() : spmats<_T,_Index>(),
	        opt, n, ok);
	    return handle;
	}

	// unsigned Index path: sparse LU cannot be used (Index must be signed);
	// same reporting convention as dispatch_sparse_lu_.
	template <typename _T, typename _Index>
	inline typename std::enable_if<!std::is_signed<_Index>::value,
	                               lu_factor_handle<_T,_Index> >::type
	dispatch_lu_factorize_(
	    const spmats<_T,_Index>& A,
	    const sparse_lu_options<_T>& opt)
	{
	    (void)A; (void)opt;
	    vcp::throw_error<vcp::state_error>(
	        "spmats::policy_lu_factorize_with_info: sparse_lu requires a signed Index type");
	    return lu_factor_handle<_T,_Index>();
	}

} // namespace spmats_lss_detail

// ---------------------------------------------------------------------------
// policy_lu_factorize_with_info: NVI outer (non-virtual).  Finalize guarantee
// + squareness entry check, then delegates to the virtual _impl.  Must never
// be overridden -- override policy_lu_factorize_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_factor_handle<_T,_Index> spmats<_T, _Index>::policy_lu_factorize_with_info(
	const sparse_lu_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_lu_factorize_with_info: matrix must be square");
	return policy_lu_factorize_with_info_impl(opt);
}

// ---------------------------------------------------------------------------
// policy_lu_factorize_with_info_impl: virtual algorithm body (default:
// signed-Index guard -> SLU factorization -> handle assembly).  Runtime
// failure is a non-valid handle; misuse (vcp::error) keeps its throwing
// contract; the final std::exception net returns a non-valid handle with
// n recorded (same P3 pattern as the other policies).  Designated
// replacement point for external-backend policies (spumar override is out
// of scope for LSS-1; the default implementation runs its own SLU there).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_factor_handle<_T,_Index> spmats<_T, _Index>::policy_lu_factorize_with_info_impl(
	const sparse_lu_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	try {
		return spmats_lss_detail::dispatch_lu_factorize_<_T, _Index>(A, opt);
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		// runtime failure net: non-valid handle, n recorded so the D5
		// solve_with_info contract (converged=false, x = 0 of size n) holds.
		lu_factor_handle<_T,_Index> handle;
		spmats_lu_factor_detail::handle_access::assign(
		    handle,
		    vcp::sparse_lu_factorization<_T,
		        typename lu_factor_handle<_T,_Index>::factor_index_type>(),
		    spmats<_T,_Index>(), opt, A.rowsize(), false);
		return handle;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_LSS_HPP
