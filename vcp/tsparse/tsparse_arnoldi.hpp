// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_ARNOLDI_HPP
#define VCP_TSPARSE_ARNOLDI_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_arnoldi {

// ---------------------------------------------------------------------------
// Result package (analogous to tsparse_lanczos::lanczos_result_package)
// ---------------------------------------------------------------------------
template <typename T>
struct arnoldi_result_package {
	typedef typename tsparse_scalar::real_type<T>::type R;
	typedef std::pair<R, R> complex_pair;

	std::vector<T> eigenvalues;                    // real eigenvalues
	std::vector<complex_pair> complex_eigenvalues; // complex pairs (re, im)
	std::vector<std::vector<T> > eigenvectors;
		std::vector<R> residuals_abs;
		std::vector<R> residuals_rel;
		std::vector<R> history_abs;
		std::vector<R> history_rel;
	bool converged;
	bool has_complex;
	std::size_t converged_count;
	std::size_t returned_count;
	std::size_t iterations;
	std::size_t restarts;
	std::size_t mv_count;
	std::string breakdown_reason;
	std::string failure_reason;
	std::string used_method;

	arnoldi_result_package()
		: converged(false), has_complex(false),
		  converged_count(0), returned_count(0),
		  iterations(0), restarts(0), mv_count(0),
		  used_method("arnoldi") {}
};

// ---------------------------------------------------------------------------
// Build one Arnoldi factorization: A V_m = V_m H_m + h_{m+1,m} v_{m+1} e_m^T
//
// Returns (V, H, m_actual, beta_last, breakdown_reason).
// Orthogonalization is dispatched to tsparse_eigensolvers::orthogonalize.
//
// V_start: starting vector (already normalized)
// ---------------------------------------------------------------------------
template <typename T, class ApplyA>
struct arnoldi_factorization_result {
	std::vector<std::vector<T> > V;       // basis V[0..m_actual]
	std::vector<std::vector<T> > H;       // (m_limit+1) × m_limit Hessenberg
	std::size_t m_actual;
	typename tsparse_scalar::real_type<T>::type beta_last;
	std::size_t mv_count;
	std::string breakdown_reason;
};

template <typename T, class ApplyA>
arnoldi_factorization_result<T, ApplyA> build_arnoldi(
	const std::size_t n,
	const std::size_t m_limit,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const orthogonalization_method orth,
	const bool full_reorthogonalization,
	const std::vector<T>& v_start,
	ApplyA apply)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	arnoldi_factorization_result<T, ApplyA> res;
	res.m_actual = 0;
	res.beta_last = R(0);
	res.mv_count = 0;

	const R small_tol = tol / R(10);
	res.V.reserve(m_limit + 1);
	res.V.push_back(v_start);
	res.H.assign(m_limit + 1, std::vector<T>(m_limit, T(0)));

	for (std::size_t j = 0; j < m_limit; j++) {
		std::vector<T> w;
		apply(res.V[j], w);
		res.mv_count++;

		// First orthogonalization pass
		std::vector<T> h_col;
		tsparse_eigensolvers::orthogonalize(res.V, j + 1, w, h_col, orth);
		h_col.resize(j + 2, T(0));

		// Optional second reorthogonalization pass
		if (full_reorthogonalization) {
			std::vector<T> h2;
			tsparse_eigensolvers::orthogonalize(res.V, j + 1, w, h2, orth);
			for (std::size_t i = 0; i <= j; i++) h_col[i] += h2[i];
		}

		res.beta_last = tsparse_scalar::real_norm_value(w);
		h_col[j + 1] = T(res.beta_last);

		for (std::size_t i = 0; i <= j + 1; i++) res.H[i][j] = h_col[i];

		res.m_actual++;

		if (res.beta_last <= small_tol) {
			res.breakdown_reason = "happy breakdown";
			break;
		}

		std::vector<T> vnew(n);
		for (std::size_t i = 0; i < n; i++) vnew[i] = w[i] / T(res.beta_last);
		res.V.push_back(vnew);
	}

	return res;
}

// ---------------------------------------------------------------------------
// Extract Hessenberg square matrix Hm (m × m portion)
// ---------------------------------------------------------------------------
template <typename T>
std::vector<std::vector<T> > extract_hm(const std::vector<std::vector<T> >& H, const std::size_t m) {
	std::vector<std::vector<T> > Hm(m, std::vector<T>(m, T(0)));
	for (std::size_t i = 0; i < m; i++)
		for (std::size_t j = 0; j < m; j++)
			Hm[i][j] = H[i][j];
	return Hm;
}

// ---------------------------------------------------------------------------
// Standard Arnoldi eigensolver with deterministic restart vectors
//
// Design: Arnoldi is for GENERAL (non-symmetric) matrices.
// Unlike Lanczos, eigenvectors are NOT mutually orthogonal, so
// orthogonal deflation does not work. We use Ritz-restart:
//   1. Build Krylov subspace of size m_limit
//   2. Extract Ritz values via Francis QR on the Hessenberg
//   3. Check convergence (Arnoldi residual bound + exact residual)
//   4. If not converged: restart from the leading Ritz vector
//   5. Repeat until k eigenvalues converge or max_restarts reached
//
// ApplyA: void(const vector<T>& x, vector<T>& y)  →  y = A x
// ApplyNorm: R(const vector<T>& v)                →  ||v||_B
// ---------------------------------------------------------------------------
template <typename T, class ApplyA, class ApplyNorm>
arnoldi_result_package<T> arnoldi_eigs(
	const std::size_t n,            // matrix dimension
	const std::size_t k,            // number of eigenvalues requested
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const orthogonalization_method orth,
	const bool full_reorthogonalization,
	const bool use_restart,
	const unsigned int random_seed,
	const bool random_start,
	const eig_target target,
	const typename tsparse_scalar::real_type<T>::type& shift_val,
	const bool compute_residual_history,
	ApplyA apply,
	ApplyNorm apply_norm)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	typedef std::complex<R> C;
	typedef std::pair<R, R> cpair;

	arnoldi_result_package<T> res;
	if (n == 0 || k == 0) return res;

	const R small_tol = tol / R(10);
	const std::size_t m_limit = std::min(n, std::max(subspace_dim, k + 2));

	// Track found eigenvalues (to avoid duplicates on restart)
	std::vector<T> found_real;
	std::vector<cpair> found_complex;
	std::size_t total_found = 0;
	bool any_complex = false;

	// Residuals for found eigenvalues
	std::vector<R> found_res_abs;
	std::vector<std::vector<T> > found_vecs;

	unsigned int seed_counter = random_start ? random_seed : 0u;

	// Starting vector for first iteration
	std::vector<T> v_start = tsparse_lanczos::deterministic_start_vector<T>(n, seed_counter++);

	for (std::size_t restart = 0; restart <= max_restarts; restart++) {
		if (total_found >= k) break;
		res.restarts = restart;

		// Build Arnoldi factorization
		auto fac = build_arnoldi<T, ApplyA>(n, m_limit, tol, orth, full_reorthogonalization,
		                                     v_start, apply);
		res.mv_count += fac.mv_count;
		res.iterations += fac.m_actual;
		if (!fac.breakdown_reason.empty()) res.breakdown_reason = fac.breakdown_reason;

		const std::size_t m = fac.m_actual;
		if (m == 0) { res.failure_reason = "empty Arnoldi basis"; break; }

		// Extract square Hessenberg and compute complex eigenvalues
		std::vector<std::vector<T> > Hm = extract_hm(fac.H, m);
		const R hess_tol = small_tol;
		const std::size_t hess_iter = m * m * 50 + 200;
		std::vector<C> ceigs = tsparse_eigensolvers::hessenberg_complex_eigenvalues<T>(
			Hm, hess_iter, hess_tol);

		if (ceigs.empty()) { res.failure_reason = "Hessenberg eigensolver failed"; break; }

		// Select target eigenvalues
		std::vector<std::size_t> sel = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs, ceigs.size(), target, shift_val);

		// The "best" Ritz vector for next restart (from top-1 selected)
		bool found_best_restart_vec = false;
		std::vector<T> best_ritz_vec;

		// Check each selected Ritz pair for convergence
		bool found_new = false;

		for (std::size_t si = 0; si < sel.size() && total_found < k; si++) {
			const std::size_t idx = sel[si];
			const C& cv = ceigs[idx];
			const R im_abs = tsparse_scalar::abs_value(cv.imag());
			const R re_abs_plus = tsparse_scalar::abs_value(cv.real()) + R(1);
			const bool is_real = im_abs <= hess_tol * re_abs_plus;

			if (is_real) {
				// Check if we already found this eigenvalue
				const R theta = cv.real();
				bool already_found = false;
				for (std::size_t fi = 0; fi < found_real.size(); fi++) {
					if (tsparse_scalar::abs_value(tsparse_scalar::real_part(found_real[fi]) - theta)
					    <= tol * (R(1) + tsparse_scalar::abs_value(theta))) {
						already_found = true; break;
					}
				}
				if (already_found) continue;

				// Compute Ritz vector
				std::vector<T> ysmall = tsparse_dense_linalg::dense_eigenvector_inverse_iteration(
					Hm, T(theta));
				if (ysmall.empty()) continue;

				std::vector<T> ritz(n, T(0));
				for (std::size_t jj = 0; jj < ysmall.size() && jj < fac.V.size(); jj++)
					for (std::size_t i = 0; i < n; i++)
						ritz[i] += fac.V[jj][i] * ysmall[jj];

				const R ritz_nrm = apply_norm(ritz);
				if (ritz_nrm <= small_tol) continue;
				for (std::size_t i = 0; i < n; i++) ritz[i] /= T(ritz_nrm);

				// Save best Ritz vector for restart (use first real pair)
				if (!found_best_restart_vec) {
					best_ritz_vec = ritz;
					found_best_restart_vec = true;
				}

				// Arnoldi residual bound: |h_{m+1,m}| * |y[m-1]|
				const R arnoldi_bound = fac.beta_last * tsparse_scalar::abs_value(ysmall.back());

				// Skip exact residual check if Arnoldi bound is too large (fast exit)
				if (arnoldi_bound > tol * R(1000) && !fac.breakdown_reason.empty()) {
					// Happy breakdown: Krylov subspace invariant → compute exact residual
				} else if (arnoldi_bound > tol * R(100) && fac.breakdown_reason.empty()) {
					continue; // definitely not converged
				}

				// Compute exact residual ||A*ritz - theta*ritz||
				std::vector<T> Ar;
				apply(ritz, Ar);
				res.mv_count++;
				const T mu = T(theta);
				std::vector<T> r_vec(n);
				for (std::size_t i = 0; i < n; i++) r_vec[i] = Ar[i] - mu * ritz[i];
				const R res_abs = tsparse_scalar::real_norm_value(r_vec);
				const R res_rel = res_abs / (R(1) + tsparse_scalar::abs_value(mu));

					if (compute_residual_history) {
						res.history_abs.push_back(res_abs);
						res.history_rel.push_back(res_rel);
					}

				if (res_abs <= tol || res_rel <= tol) {
					found_real.push_back(mu);
					found_res_abs.push_back(res_abs);
					found_vecs.push_back(ritz);
					total_found++;
					found_new = true;
					res.converged_count++;
				}
			} else {
				// Complex eigenvalue pair
				const R re_part = cv.real();
				const R im_part = tsparse_scalar::abs_value(cv.imag());

				// Check if we already found this complex pair
				bool already_found = false;
				for (std::size_t fi = 0; fi < found_complex.size(); fi++) {
					if (tsparse_scalar::abs_value(found_complex[fi].first - re_part) <= tol * R(10)
					 && tsparse_scalar::abs_value(found_complex[fi].second - im_part) <= tol * R(10)) {
						already_found = true; break;
					}
				}
				if (already_found) continue;

				// Arnoldi residual bound for the complex pair
				// Use |h_{m+1,m}| as conservative upper bound
				// (exact bound requires complex eigenvector computation)
				const R arnoldi_bound = fac.beta_last;

				if (arnoldi_bound <= tol || !fac.breakdown_reason.empty()) {
					found_complex.push_back(cpair(re_part, im_part));
					any_complex = true;
					total_found += 2; // conjugate pair
					res.converged_count += 2;
					found_new = true;
				}
			}
		}

		// Check if we found all requested eigenvalues
		if (total_found >= k) {
			res.converged = !any_complex;
			break;
		}

		// If no new eigenvalues found on this restart, try a different start
		if (!found_new || !use_restart) {
			// Generate fresh start (orthogonal to already-found eigenvectors if possible)
			v_start = tsparse_lanczos::deterministic_start_vector<T>(n, seed_counter++);
			// Soft deflation: try to remove found components
			for (std::size_t fi = 0; fi < found_vecs.size(); fi++) {
				const R c = tsparse_scalar::real_dot_value(found_vecs[fi], v_start);
				for (std::size_t i = 0; i < n; i++) v_start[i] -= T(c) * found_vecs[fi][i];
			}
			const R nv = tsparse_scalar::real_norm_value(v_start);
			if (nv > small_tol) {
				for (std::size_t i = 0; i < n; i++) v_start[i] /= T(nv);
			} else {
				v_start = tsparse_lanczos::deterministic_start_vector<T>(n, seed_counter++);
			}
		} else {
			// Ritz restart: start from the best unconverged Ritz vector
			if (found_best_restart_vec) {
				v_start = best_ritz_vec;
				// Soft deflation for already-found eigenvalues
				for (std::size_t fi = 0; fi < found_vecs.size(); fi++) {
					const R c = tsparse_scalar::real_dot_value(found_vecs[fi], v_start);
					for (std::size_t i = 0; i < n; i++) v_start[i] -= T(c) * found_vecs[fi][i];
				}
				const R nv = tsparse_scalar::real_norm_value(v_start);
				if (nv > small_tol) {
					for (std::size_t i = 0; i < n; i++) v_start[i] /= T(nv);
				} else {
					v_start = tsparse_lanczos::deterministic_start_vector<T>(n, seed_counter++);
				}
			} else {
				v_start = tsparse_lanczos::deterministic_start_vector<T>(n, seed_counter++);
			}
		}
	}

	// Sort found real eigenvalues by target
	if (!found_real.empty()) {
		std::vector<C> ceigs_found;
		ceigs_found.reserve(found_real.size());
		for (std::size_t i = 0; i < found_real.size(); i++)
			ceigs_found.push_back(C(tsparse_scalar::real_part(found_real[i]), R(0)));

		const std::size_t ksorted = std::min(k, found_real.size());
		std::vector<std::size_t> order = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs_found, ksorted, target, shift_val);

		for (std::size_t i = 0; i < order.size(); i++) {
			res.eigenvalues.push_back(found_real[order[i]]);
			res.residuals_abs.push_back(found_res_abs[order[i]]);
			if (!found_vecs.empty() && order[i] < found_vecs.size())
				res.eigenvectors.push_back(found_vecs[order[i]]);
		}
		res.residuals_rel.resize(res.residuals_abs.size());
		for (std::size_t i = 0; i < res.eigenvalues.size(); i++)
			res.residuals_rel[i] = res.residuals_abs[i] / (R(1) + tsparse_scalar::abs_value(res.eigenvalues[i]));
		res.returned_count = res.eigenvalues.size();
	}

	res.complex_eigenvalues = found_complex;
	res.has_complex = any_complex;
	res.converged = res.converged && (total_found >= k);

	if (!res.converged && res.failure_reason.empty()) {
		if (any_complex)
			res.failure_reason = "complex eigenvalues in requested subset";
		else
			res.failure_reason = "maximum restarts reached without full convergence";
	}

	return res;
}

// ---------------------------------------------------------------------------
// Convenience overload with standard L2 norm
// ---------------------------------------------------------------------------
template <typename T, class ApplyA>
arnoldi_result_package<T> arnoldi_eigs_standard(
	const std::size_t n,
	const std::size_t k,
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const orthogonalization_method orth,
	const bool full_reorthogonalization,
	const bool use_restart,
	const unsigned int random_seed,
	const bool random_start,
	const eig_target target,
	const typename tsparse_scalar::real_type<T>::type& shift_val,
	const bool compute_residual_history,
	ApplyA apply)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	struct norm_fn {
		R operator()(const std::vector<T>& v) const {
			return tsparse_scalar::real_norm_value(v);
		}
	} norm_func;
	return arnoldi_eigs<T, ApplyA, norm_fn>(
		n, k, subspace_dim, max_restarts, tol,
		orth, full_reorthogonalization, use_restart,
		random_seed, random_start, target, shift_val,
		compute_residual_history, apply, norm_func);
}

} // namespace tsparse_arnoldi
} // namespace vcp

#endif
