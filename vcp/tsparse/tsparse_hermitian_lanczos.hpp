// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_HERMITIAN_LANCZOS_HPP
#define VCP_TSPARSE_HERMITIAN_LANCZOS_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_projected_eigensolver.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_hermitian_lanczos {

namespace detail {
// Compile-time T construction from two real values.
// For real T: ignore im, return re. For complex T: return T(re, im).
template <typename T, bool IsComplex = vcp::tsparse_scalar::is_complex<T>::value>
struct make_from_reals {
	typedef typename vcp::tsparse_scalar::real_type<T>::type R;
	static T apply(R re, R im) { return T(re, im); }
};
template <typename T>
struct make_from_reals<T, false> {
	typedef typename vcp::tsparse_scalar::real_type<T>::type R;
	static T apply(R re, R /*im*/) { return re; }
};
} // namespace detail

// Hermitian inner product: <x, y> = sum_i conj(x[i]) * y[i]
template <typename T>
T hermitian_inner_product(const std::vector<T>& x, const std::vector<T>& y) {
	T sum = T(0);
	for (std::size_t i = 0; i < x.size(); i++) {
		sum += vcp::tsparse_scalar::conjugate_if_needed(x[i]) * y[i];
	}
	return sum;
}

// Hermitian norm: sqrt(real(<x, x>))
template <typename T>
typename vcp::tsparse_scalar::real_type<T>::type hermitian_norm(const std::vector<T>& x) {
	typedef typename vcp::tsparse_scalar::real_type<T>::type R;
	const T v = hermitian_inner_product(x, x);
	const R r = vcp::tsparse_scalar::real_part(v);
	return vcp::tsparse_scalar::sqrt_value(r > R(0) ? r : R(0));
}

// Deterministic LCG-based starting vector; works for both real and complex T.
template <typename T>
std::vector<T> complex_start_vector(const std::size_t n, const unsigned int seed) {
	typedef typename vcp::tsparse_scalar::real_type<T>::type R;
	std::vector<T> v(n, T(0));
	unsigned int state = seed + 1u;
	for (std::size_t i = 0; i < n; i++) {
		state = state * 1664525u + 1013904223u;
		const R re = R(static_cast<int>(state >> 16)) / R(32768) - R(1);
		state = state * 1664525u + 1013904223u;
		const R im = R(static_cast<int>(state >> 16)) / R(32768) - R(1);
		v[i] = detail::make_from_reals<T>::apply(re, im);
	}
	const R nv = hermitian_norm(v);
	if (nv > R(0)) {
		for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	}
	return v;
}

// Orthogonalize v against a set of Hermitian-orthonormal vectors (2-pass CGS)
// Returns false if resulting vector is near zero.
template <typename T>
bool hermitian_orthogonalize_against(
	std::vector<T>& v,
	const std::vector<std::vector<T> >& vecs,
	const typename vcp::tsparse_scalar::real_type<T>::type& tol)
{
	const std::size_t n = v.size();
	for (int pass = 0; pass < 2; pass++) {
		for (std::size_t j = 0; j < vecs.size(); j++) {
			const T c = hermitian_inner_product(vecs[j], v);
			for (std::size_t i = 0; i < n; i++) v[i] -= c * vecs[j][i];
		}
	}
	const typename vcp::tsparse_scalar::real_type<T>::type nv = hermitian_norm(v);
	if (nv <= tol) return false;
	for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	return true;
}

// Result package from Hermitian Lanczos
template <typename T>
struct hermitian_lanczos_result {
	typedef typename vcp::tsparse_scalar::real_type<T>::type R;
	std::vector<T> eigenvalues;              // stored as T(lambda_real, 0)
	std::vector<std::vector<T> > eigenvectors;
	std::vector<R> residuals_abs;
	std::vector<R> residuals_rel;
	std::vector<R> history_abs;
	std::vector<R> history_rel;
	bool converged;
	std::size_t converged_count;
	std::size_t returned_count;
	std::size_t iterations;
	std::size_t restarts;
	std::size_t mv_count;
	std::string breakdown_reason;
	std::string failure_reason;
	std::string used_method;

	hermitian_lanczos_result()
		: converged(false), converged_count(0), returned_count(0),
		  iterations(0), restarts(0), mv_count(0),
		  used_method("hermitian_lanczos") {}
};

// Full-reorthogonalized Hermitian Lanczos eigensolver.
//
// T is typically std::complex<R> but also compiles for real T.
// ApplyA: void(const vector<T>& x, vector<T>& y) — compute y = A*x.
//
// Eigenvalues returned as T(lambda_real, 0).
template <typename T, class ApplyA>
hermitian_lanczos_result<T> hermitian_lanczos_eigs(
	const std::size_t n,
	const std::size_t k,
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename vcp::tsparse_scalar::real_type<T>::type& tol,
	const unsigned int random_seed,
	const bool random_start,
	const vcp::eig_target target,
	const typename vcp::tsparse_scalar::real_type<T>::type& shift_val,
	const typename vcp::tsparse_scalar::real_type<T>::type& normA_fro,
	const bool compute_residual_history,
	ApplyA apply)
{
	typedef typename vcp::tsparse_scalar::real_type<T>::type R;
	hermitian_lanczos_result<T> res;

	if (n == 0 || k == 0) {
		res.converged = (k == 0);
		return res;
	}

	const std::size_t m_limit = std::min(n, std::max(subspace_dim, k + std::size_t(3)));
	const R small_tol = tol / R(10);
	const R eps_R = std::numeric_limits<R>::epsilon();

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

		std::vector<T> v0 = complex_start_vector<T>(n, seed_counter);
		seed_counter++;
		if (!hermitian_orthogonalize_against(v0, locked_vecs, small_tol)) {
			bool ok = false;
			for (int attempt = 0; attempt < 20 && !ok; attempt++) {
				v0 = complex_start_vector<T>(n, seed_counter++);
				ok = hermitian_orthogonalize_against(v0, locked_vecs, small_tol);
			}
			if (!ok) {
				res.failure_reason = "exhausted orthogonal starting vectors";
				break;
			}
		}

		const std::size_t m_active = std::min(n - locked_vals.size(), m_limit);
		if (m_active == 0) break;

		std::vector<std::vector<T> > basis;
		basis.reserve(m_active + 1);
		basis.push_back(v0);

		std::vector<R> alpha_vec;
		std::vector<R> beta_vec;
		std::vector<T> v_prev(n, T(0));
		R beta_prev(0);
		std::string bd_reason;

		for (std::size_t j = 0; j < m_active; j++) {
			std::vector<T> z;
			apply(basis[j], z);
			res.mv_count++;

			if (j > 0) {
				for (std::size_t i = 0; i < n; i++) z[i] -= T(beta_prev) * v_prev[i];
			}

			for (std::size_t li = 0; li < locked_vecs.size(); li++) {
				const T c = hermitian_inner_product(locked_vecs[li], z);
				for (std::size_t i = 0; i < n; i++) z[i] -= c * locked_vecs[li][i];
			}

			const T hdot = hermitian_inner_product(basis[j], z);
			const R a = vcp::tsparse_scalar::real_part(hdot);
			alpha_vec.push_back(a);

			for (std::size_t i = 0; i < n; i++) z[i] -= T(a) * basis[j][i];

			// Full reorthogonalization (2-pass)
			for (int pass = 0; pass < 2; pass++) {
				for (std::size_t bi = 0; bi < basis.size(); bi++) {
					const T c = hermitian_inner_product(basis[bi], z);
					for (std::size_t i = 0; i < n; i++) z[i] -= c * basis[bi][i];
				}
				for (std::size_t li = 0; li < locked_vecs.size(); li++) {
					const T c = hermitian_inner_product(locked_vecs[li], z);
					for (std::size_t i = 0; i < n; i++) z[i] -= c * locked_vecs[li][i];
				}
			}

			const R bn = hermitian_norm(z);
			if (compute_residual_history) {
				res.history_abs.push_back(bn);
				// Standard definition: ||v||=1 (basis is normalized), denom = ||A||_F + |a| + eps
				const R step_denom = normA_fro + vcp::tsparse_scalar::abs_value(a) + eps_R;
				res.history_rel.push_back(bn / step_denom);
			}

			if (bn <= small_tol) {
				bd_reason = "happy breakdown";
				break;
			}

			beta_vec.push_back(bn);
			v_prev = basis[j];
			beta_prev = bn;
			std::vector<T> vnew(n);
			for (std::size_t i = 0; i < n; i++) vnew[i] = z[i] / T(bn);
			basis.push_back(vnew);
		}

		res.iterations += alpha_vec.size();

		if (alpha_vec.empty()) {
			res.breakdown_reason = "empty Lanczos basis";
			break;
		}

		// Solve real symmetric tridiagonal projected problem via solve_real_symmetric_projected
		const std::size_t m = alpha_vec.size();
		std::vector<std::vector<R> > Tm(m, std::vector<R>(m, R(0)));
		for (std::size_t i = 0; i < m; i++) {
			Tm[i][i] = alpha_vec[i];
			if (i + 1 < m) {
				Tm[i][i + 1] = beta_vec[i];
				Tm[i + 1][i] = beta_vec[i];
			}
		}
		const std::size_t proj_iter = std::max(m * m * std::size_t(100), std::size_t(1000));
		vcp::tsparse_projected::projected_eigensolver_result<R> proj =
			vcp::tsparse_projected::solve_real_symmetric_projected<R>(Tm, proj_iter, small_tol);

		if (!proj.success) {
			res.breakdown_reason = proj.status.empty() ? "projected_eigensolver_failed" : proj.status;
			res.failure_reason = proj.message.empty() ? "projected eigensolver failed" : proj.message;
			break;
		}
		if (proj.pairs.empty()) {
			res.breakdown_reason = "projected_eigensolver_empty_result";
			res.failure_reason = "projected eigensolver returned empty result";
			break;
		}

		const std::size_t k_remaining = k - locked_vals.size();
		std::vector<std::size_t> sel = vcp::tsparse_projected::select_projected_indices(
			proj.pairs, std::min(m, k_remaining + m), target, shift_val);

		for (std::size_t si = 0; si < sel.size() && locked_vals.size() < k; si++) {
			const std::size_t idx = sel[si];
			if (idx >= proj.pairs.size()) continue;

			const std::vector<R>& y = proj.pairs[idx].vector;
			std::vector<T> ritz_vec(n, T(0));
			for (std::size_t j = 0; j < y.size() && j < basis.size(); j++) {
				for (std::size_t i = 0; i < n; i++) ritz_vec[i] += basis[j][i] * T(y[j]);
			}

			for (std::size_t li = 0; li < locked_vecs.size(); li++) {
				const T c = hermitian_inner_product(locked_vecs[li], ritz_vec);
				for (std::size_t i = 0; i < n; i++) ritz_vec[i] -= c * locked_vecs[li][i];
			}

			const R ritz_nrm = hermitian_norm(ritz_vec);
			if (ritz_nrm <= small_tol) continue;
			for (std::size_t i = 0; i < n; i++) ritz_vec[i] /= T(ritz_nrm);

			std::vector<T> Av;
			apply(ritz_vec, Av);
			res.mv_count++;

			const R lambda_r = vcp::tsparse_scalar::real_part(proj.pairs[idx].value);
			std::vector<T> r(n);
			for (std::size_t i = 0; i < n; i++) r[i] = Av[i] - T(lambda_r) * ritz_vec[i];
			const R res_abs = hermitian_norm(r);
			// ritz_vec is unit-norm: ||v||=1, so denom = ||A||_F + |lambda_r| + eps
			const R ritz_denom = normA_fro + vcp::tsparse_scalar::abs_value(lambda_r) + eps_R;
			const R res_rel = res_abs / ritz_denom;

			if (compute_residual_history) {
				res.history_abs.push_back(res_abs);
				res.history_rel.push_back(res_rel);
			}

			if (best_vals.size() < k) {
				best_vals.push_back(T(lambda_r));
				best_vecs.push_back(ritz_vec);
				best_res.push_back(res_abs);
			}

			if (res_abs <= tol || res_rel <= tol) {
				locked_vals.push_back(T(lambda_r));
				locked_vecs.push_back(ritz_vec);
				locked_res.push_back(res_abs);
				res.converged_count++;
				break; // one lock per restart
			}
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
			ceigs_locked.push_back(std::complex<R>(
				vcp::tsparse_scalar::real_part(locked_vals[i]), R(0)));
		}
		const std::size_t ksorted = std::min(k, locked_vals.size());
		const std::vector<std::size_t> order = vcp::tsparse_eigensolvers::select_ritz_indices<R>(
			ceigs_locked, ksorted, target, shift_val);

		std::vector<T> sorted_vals;
		std::vector<std::vector<T> > sorted_vecs;
		std::vector<R> sorted_res;
		sorted_vals.reserve(order.size());
		sorted_vecs.reserve(order.size());
		sorted_res.reserve(order.size());
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
			const R lr = vcp::tsparse_scalar::real_part(sorted_vals[i]);
			// eigenvectors are unit-norm: denom = ||A||_F + |lambda| + eps
			const R denom = normA_fro + vcp::tsparse_scalar::abs_value(lr) + eps_R;
			res.residuals_rel[i] = sorted_res[i] / denom;
		}
		res.returned_count = sorted_vals.size();
	}

	if (!res.converged && res.failure_reason.empty()) {
		res.failure_reason = "maximum restarts reached without full convergence";
	}

	// Fallback to best-seen if nothing locked
	if (res.eigenvalues.empty() && !best_vals.empty()) {
		std::vector<std::complex<R> > ceigs_best;
		ceigs_best.reserve(best_vals.size());
		for (std::size_t i = 0; i < best_vals.size(); i++) {
			ceigs_best.push_back(std::complex<R>(
				vcp::tsparse_scalar::real_part(best_vals[i]), R(0)));
		}
		const std::vector<std::size_t> order =
			vcp::tsparse_eigensolvers::select_ritz_indices<R>(
				ceigs_best, std::min(k, best_vals.size()), target, shift_val);
		for (std::size_t i = 0; i < order.size(); i++) {
			res.eigenvalues.push_back(best_vals[order[i]]);
			res.eigenvectors.push_back(best_vecs[order[i]]);
			res.residuals_abs.push_back(best_res[order[i]]);
		}
		res.residuals_rel.resize(res.residuals_abs.size());
		for (std::size_t i = 0; i < res.eigenvalues.size(); i++) {
			const R lr = vcp::tsparse_scalar::real_part(res.eigenvalues[i]);
			// eigenvectors are unit-norm: denom = ||A||_F + |lambda| + eps
			const R denom = normA_fro + vcp::tsparse_scalar::abs_value(lr) + eps_R;
			res.residuals_rel[i] = res.residuals_abs[i] / denom;
		}
		res.returned_count = res.eigenvalues.size();
	}

	return res;
}

} // namespace tsparse_hermitian_lanczos
} // namespace vcp

#endif
