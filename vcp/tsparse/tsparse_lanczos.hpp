// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_LANCZOS_HPP
#define VCP_TSPARSE_LANCZOS_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_lanczos {

// ---------------------------------------------------------------------------
// Generate a deterministic starting vector (based on indices)
// ---------------------------------------------------------------------------
template <typename T>
std::vector<T> deterministic_start_vector(const std::size_t n, const unsigned int seed) {
	typedef typename tsparse_scalar::real_type<T>::type R;
	std::vector<T> v(n, T(0));
	// LCG-based pseudorandom (simple, deterministic, no <random> dependency issue)
	unsigned int state = seed + 1u;
	for (std::size_t i = 0; i < n; i++) {
		state = state * 1664525u + 1013904223u;
		// Map to [-1, 1]
		const R val = R(static_cast<int>(state >> 16)) / R(32768) - R(1);
		v[i] = T(val);
	}
	const R nv = tsparse_scalar::real_norm_value(v);
	if (nv > R(0)) {
		for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	}
	return v;
}

// ---------------------------------------------------------------------------
// Orthogonalize v against a list of vectors (Gram-Schmidt, full)
// Returns false if resulting vector is (near) zero.
// ---------------------------------------------------------------------------
template <typename T>
bool orthogonalize_against(std::vector<T>& v,
                            const std::vector<std::vector<T> >& vecs,
                            const typename tsparse_scalar::real_type<T>::type& tol)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	const std::size_t n = v.size();
	// Two passes for numerical stability
	for (int pass = 0; pass < 2; pass++) {
		for (std::size_t j = 0; j < vecs.size(); j++) {
			const R c = tsparse_scalar::real_dot_value(vecs[j], v);
			for (std::size_t i = 0; i < n; i++) v[i] -= T(c) * vecs[j][i];
		}
	}
	const R nv = tsparse_scalar::real_norm_value(v);
	if (nv <= tol) return false;
	for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	return true;
}

// ---------------------------------------------------------------------------
// Symmetric tridiagonal eigensolver (Jacobi-based, small m)
// Returns eigenvalues and eigenvectors of tridiagonal T(alpha, beta)
// ---------------------------------------------------------------------------
template <typename T>
tsparse_dense_linalg::dense_eigen_result<T> solve_tridiag(
	const std::vector<T>& alpha,
	const std::vector<T>& beta,
	const std::size_t max_iter,
	const typename tsparse_scalar::real_type<T>::type& tol)
{
	const std::size_t m = alpha.size();
	std::vector<std::vector<T> > Tm(m, std::vector<T>(m, T(0)));
	for (std::size_t i = 0; i < m; i++) {
		Tm[i][i] = alpha[i];
		if (i + 1 < m) { Tm[i][i+1] = beta[i]; Tm[i+1][i] = beta[i]; }
	}
	return tsparse_dense_linalg::jacobi_eig_dense(Tm, max_iter, tol);
}

// ---------------------------------------------------------------------------
// Full-reorthogonalized standard Lanczos eigensolver
//
// ApplyA: void(const vector<T>& x, vector<T>& y)  →  y = A x
// Returns eig_result<T> with all required fields populated.
// ---------------------------------------------------------------------------
template <typename T, class ApplyA>
struct lanczos_result_package {
	std::vector<T> eigenvalues;
	std::vector<std::vector<T> > eigenvectors;
		std::vector<typename tsparse_scalar::real_type<T>::type> residuals_abs;
		std::vector<typename tsparse_scalar::real_type<T>::type> residuals_rel;
		std::vector<typename tsparse_scalar::real_type<T>::type> history_abs;
		std::vector<typename tsparse_scalar::real_type<T>::type> history_rel;
	bool converged;
	std::size_t converged_count;
	std::size_t returned_count;
	std::size_t iterations;
	std::size_t restarts;
	std::size_t mv_count;
	std::string breakdown_reason;
	std::string failure_reason;
	std::string used_method;
};

template <typename T, class ApplyA, class ApplyNorm>
lanczos_result_package<T, ApplyA> lanczos_eigs(
	const std::size_t n,           // matrix size
	const std::size_t k,           // number requested
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const unsigned int random_seed,
	const bool random_start,
	const eig_target target,
	const typename tsparse_scalar::real_type<T>::type& shift_val,
	const bool compute_residual_history,
	ApplyA apply,                  // y = A x
	ApplyNorm apply_norm)          // ||v||_B (for standard: standard norm)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	lanczos_result_package<T, ApplyA> res;
	res.converged = false;
	res.converged_count = 0;
	res.returned_count = 0;
	res.iterations = 0;
	res.restarts = 0;
	res.mv_count = 0;
	res.used_method = "lanczos";

	if (n == 0 || k == 0) return res;

	const std::size_t m_limit = std::min(n, std::max(subspace_dim, k + std::size_t(3)));
	const R small_tol = tol / R(10);

	// locked eigenvalues and their eigenvectors
	std::vector<T> locked_vals;
	std::vector<std::vector<T> > locked_vecs;
	std::vector<R> locked_res;
	std::vector<T> best_vals;
	std::vector<std::vector<T> > best_vecs;
	std::vector<R> best_res;

	unsigned int seed_counter = random_start ? random_seed : 0u;

	for (std::size_t restart = 0; restart <= max_restarts; restart++) {
		if (locked_vals.size() >= k) break;
		res.restarts = restart;

		// Build starting vector orthogonal to locked vectors
		std::vector<T> v0 = deterministic_start_vector<T>(n, seed_counter);
		seed_counter++;
		if (!orthogonalize_against(v0, locked_vecs, small_tol)) {
			// Try a few more seeds
			bool ok = false;
			for (int attempt = 0; attempt < 20 && !ok; attempt++) {
				v0 = deterministic_start_vector<T>(n, seed_counter++);
				ok = orthogonalize_against(v0, locked_vecs, small_tol);
			}
			if (!ok) { res.failure_reason = "exhausted orthogonal starting vectors"; break; }
		}

		// Run Lanczos from v0
		const std::size_t m_active = std::min(n - locked_vals.size(), m_limit);
		if (m_active == 0) break;

		std::vector<std::vector<T> > basis;
		basis.reserve(m_active + 1);
		basis.push_back(v0);

		std::vector<T> alpha_vec, beta_vec;
		std::vector<T> v_prev(n, T(0));
		R beta_prev = R(0);
		std::string bd_reason;

		for (std::size_t j = 0; j < m_active; j++) {
			// z = A * basis[j]
			std::vector<T> z;
			apply(basis[j], z);
			res.mv_count++;

			if (j > 0) {
				for (std::size_t i = 0; i < n; i++) z[i] -= T(beta_prev) * v_prev[i];
			}

			// Orthogonalize z against locked vectors
			for (std::size_t li = 0; li < locked_vecs.size(); li++) {
				const R c = tsparse_scalar::real_dot_value(locked_vecs[li], z);
				for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * locked_vecs[li][i];
			}

			const R a = tsparse_scalar::real_dot_value(basis[j], z);
			alpha_vec.push_back(T(a));
			for (std::size_t i = 0; i < n; i++) z[i] -= T(a) * basis[j][i];

			for (int pass = 0; pass < 2; pass++) {
				for (std::size_t bi = 0; bi < basis.size(); bi++) {
					const R c = tsparse_scalar::real_dot_value(basis[bi], z);
					for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * basis[bi][i];
				}
				for (std::size_t li = 0; li < locked_vecs.size(); li++) {
					const R c = tsparse_scalar::real_dot_value(locked_vecs[li], z);
					for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * locked_vecs[li][i];
				}
			}

			const R bn = apply_norm(z);
			if (compute_residual_history) {
				res.history_abs.push_back(bn);
				res.history_rel.push_back(bn / (R(1) + tsparse_scalar::abs_value(a)));
			}

			if (bn <= small_tol) {
				bd_reason = "happy breakdown";
				break;
			}

			beta_vec.push_back(T(bn));
			v_prev = basis[j];
			beta_prev = bn;
			std::vector<T> vnew(n);
			for (std::size_t i = 0; i < n; i++) vnew[i] = z[i] / T(bn);
			basis.push_back(vnew);
		}

		res.iterations += alpha_vec.size();

		if (alpha_vec.empty()) { res.breakdown_reason = "empty Lanczos basis"; break; }

		// Solve small tridiagonal eigenproblem
		const std::size_t m = alpha_vec.size();
		const std::size_t proj_iter = std::max(m * m * 100, std::size_t(1000));
		tsparse_dense_linalg::dense_eigen_result<T> small =
			solve_tridiag(alpha_vec, beta_vec, proj_iter, small_tol);

		if (small.eigenvalues.empty()) { res.breakdown_reason = "projected eigensolver failed"; break; }

		// Select target eigenvalues
		std::vector<std::complex<R> > ceigs;
		ceigs.reserve(small.eigenvalues.size());
		for (std::size_t i = 0; i < small.eigenvalues.size(); i++) {
			ceigs.push_back(std::complex<R>(tsparse_scalar::real_part(small.eigenvalues[i]), R(0)));
		}
		const std::size_t k_remaining = k - locked_vals.size();
		std::vector<std::size_t> sel = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs, std::min(m, k_remaining + m), target, shift_val);

		// Check convergence for each selected Ritz pair.
		// Lock at most ONE new pair per restart so that repeated eigenvalues
		// are found on subsequent restarts (each with a new deflated start vector).
		for (std::size_t si = 0; si < sel.size() && locked_vals.size() < k; si++) {
			const std::size_t idx = sel[si];
			if (idx >= small.eigenvectors.size()) continue;

			// Lift Ritz vector
			const std::vector<T>& y = small.eigenvectors[idx];
			std::vector<T> ritz_vec(n, T(0));
			for (std::size_t j = 0; j < y.size() && j < basis.size(); j++) {
				for (std::size_t i = 0; i < n; i++) ritz_vec[i] += basis[j][i] * y[j];
			}
			// Normalize and orthogonalize against already-locked vectors
			for (std::size_t li = 0; li < locked_vecs.size(); li++) {
				const R c = tsparse_scalar::real_dot_value(locked_vecs[li], ritz_vec);
				for (std::size_t i = 0; i < n; i++) ritz_vec[i] -= T(c) * locked_vecs[li][i];
			}
			const R ritz_nrm = apply_norm(ritz_vec);
			if (ritz_nrm <= small_tol) continue;
			for (std::size_t i = 0; i < n; i++) ritz_vec[i] /= T(ritz_nrm);

			// Compute residual ||A v - mu v||
			std::vector<T> Av;
			apply(ritz_vec, Av);
			res.mv_count++;
			const T mu = small.eigenvalues[idx];
			std::vector<T> r(n);
			for (std::size_t i = 0; i < n; i++) r[i] = Av[i] - mu * ritz_vec[i];
			const R res_abs = tsparse_scalar::real_norm_value(r);
			const R res_rel = res_abs / (R(1) + tsparse_scalar::abs_value(mu));

				if (compute_residual_history) {
					res.history_abs.push_back(res_abs);
					res.history_rel.push_back(res_rel);
				}
			if (best_vals.size() < k) {
				best_vals.push_back(mu);
				best_vecs.push_back(ritz_vec);
				best_res.push_back(res_abs);
			}

			if (res_abs <= tol || res_rel <= tol) {
				locked_vals.push_back(mu);
				locked_vecs.push_back(ritz_vec);
				locked_res.push_back(res_abs);
				res.converged_count++;
				break;  // one lock per restart → next restart deflates this vector out
			}
		}

		if (!bd_reason.empty() && locked_vals.size() < k) {
			// Happy breakdown with insufficient locked pairs → need restart
			// Keep going to next restart
		}

		if (locked_vals.size() >= k) {
			res.converged = true;
			res.breakdown_reason = bd_reason;
			break;
		}
	}

	// Sort locked eigenvalues by target
	if (!locked_vals.empty()) {
		std::vector<std::complex<R> > ceigs_locked;
		ceigs_locked.reserve(locked_vals.size());
		for (std::size_t i = 0; i < locked_vals.size(); i++) {
			ceigs_locked.push_back(std::complex<R>(tsparse_scalar::real_part(locked_vals[i]), R(0)));
		}
		const std::size_t ksorted = std::min(k, locked_vals.size());
		std::vector<std::size_t> order = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs_locked, ksorted, target, shift_val);

		std::vector<T> sorted_vals;
		std::vector<std::vector<T> > sorted_vecs;
		std::vector<R> sorted_res;
		for (std::size_t i = 0; i < order.size(); i++) {
			sorted_vals.push_back(locked_vals[order[i]]);
			sorted_vecs.push_back(locked_vecs[order[i]]);
			sorted_res.push_back(locked_res[order[i]]);
		}
		res.eigenvalues = sorted_vals;
		res.eigenvectors = sorted_vecs;
		res.residuals_abs = sorted_res;
		res.residuals_rel.resize(sorted_res.size());
		for (std::size_t i = 0; i < sorted_vals.size(); i++) {
			res.residuals_rel[i] = sorted_res[i] / (R(1) + tsparse_scalar::abs_value(sorted_vals[i]));
		}
		res.returned_count = sorted_vals.size();
	}

	if (!res.converged && res.failure_reason.empty()) {
		res.failure_reason = "maximum restarts reached without full convergence";
	}
	if (res.eigenvalues.empty() && !best_vals.empty()) {
		std::vector<std::complex<R> > ceigs_best;
		ceigs_best.reserve(best_vals.size());
		for (std::size_t i = 0; i < best_vals.size(); i++) {
			ceigs_best.push_back(std::complex<R>(tsparse_scalar::real_part(best_vals[i]), R(0)));
		}
		const std::vector<std::size_t> order = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs_best, std::min(k, best_vals.size()), target, shift_val);
		for (std::size_t i = 0; i < order.size(); i++) {
			res.eigenvalues.push_back(best_vals[order[i]]);
			res.eigenvectors.push_back(best_vecs[order[i]]);
			res.residuals_abs.push_back(best_res[order[i]]);
		}
		res.residuals_rel.resize(res.residuals_abs.size());
		for (std::size_t i = 0; i < res.eigenvalues.size(); i++) {
			res.residuals_rel[i] = res.residuals_abs[i] / (R(1) + tsparse_scalar::abs_value(res.eigenvalues[i]));
		}
		res.returned_count = res.eigenvalues.size();
	}

	return res;
}

// Convenience overload with standard L2 norm
template <typename T, class ApplyA>
lanczos_result_package<T, ApplyA> lanczos_eigs_standard(
	const std::size_t n,
	const std::size_t k,
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
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
	return lanczos_eigs<T, ApplyA, norm_fn>(
		n, k, subspace_dim, max_restarts, tol,
		random_seed, random_start, target, shift_val,
		compute_residual_history, apply, norm_func);
}

} // namespace tsparse_lanczos
} // namespace vcp

#endif
