// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_B_INNER_LANCZOS_HPP
#define VCP_TSPARSE_B_INNER_LANCZOS_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_b_inner_lanczos {

// ---------------------------------------------------------------------------
// b_inner_product: compute x^T B y via sparse SpMV (B.mul_vec(y))
//
// Template parameters:
//   SparseMatrix - must support mul_vec, rowsize, columnsize
//   T            - scalar type (real only; real_type<T>::type == T required)
//
// Throws vcp::dimension_error if sizes mismatch.
// ---------------------------------------------------------------------------
template <class SparseMatrix, class T>
typename vcp::tsparse_scalar::real_type<T>::type
b_inner_product(
    const SparseMatrix& B,
    const std::vector<T>& x,
    const std::vector<T>& y)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    if (x.size() != static_cast<std::size_t>(B.rowsize()) ||
        y.size() != static_cast<std::size_t>(B.columnsize())) {
        vcp::throw_error<vcp::dimension_error>(
            "tsparse_b_inner_lanczos::b_inner_product: dimension mismatch");
    }
    const std::vector<T> By = B.mul_vec(y);
    R sum = R(0);
    for (std::size_t i = 0; i < x.size(); i++) {
        sum += x[i] * By[i];
    }
    return sum;
}

// ---------------------------------------------------------------------------
// b_norm: compute sqrt(x^T B x) via sparse SpMV
//
// Returns 0 when x^T B x <= 0 (caller checks for SPD failure).
// Throws vcp::dimension_error if sizes mismatch.
// ---------------------------------------------------------------------------
template <class SparseMatrix, class T>
typename vcp::tsparse_scalar::real_type<T>::type
b_norm(
    const SparseMatrix& B,
    const std::vector<T>& x)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    if (x.size() != static_cast<std::size_t>(B.rowsize())) {
        vcp::throw_error<vcp::dimension_error>(
            "tsparse_b_inner_lanczos::b_norm: dimension mismatch");
    }
    const R v = b_inner_product<SparseMatrix, T>(B, x, x);
    if (v <= R(0)) return R(0);
    return vcp::tsparse_scalar::sqrt_value(v);
}

// ---------------------------------------------------------------------------
// b_lanczos_result: result package returned by b_inner_lanczos_eigs
// ---------------------------------------------------------------------------
template <typename T>
struct b_lanczos_result {
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    std::vector<T> eigenvalues;
    std::vector<std::vector<T> > eigenvectors;
    std::vector<R> residuals_abs;     // generalized absolute residuals
    std::vector<R> residuals_rel;     // generalized relative residuals
    std::vector<R> history_abs;       // per-step absolute residual history
    std::vector<R> history_rel;       // per-step relative residual history (NOT a copy of history_abs)
    bool converged;
    bool spd_check_failed;
    bool budget_exhausted;            // true when mv_budget was reached before convergence
    std::size_t converged_count;
    std::size_t returned_count;
    std::size_t iterations;
    std::size_t restarts;
    std::size_t mv_count;             // total sparse matrix-vector products (A and B combined)
    std::string breakdown_reason;
    std::string failure_reason;

    b_lanczos_result()
        : converged(false), spd_check_failed(false), budget_exhausted(false),
          converged_count(0), returned_count(0),
          iterations(0), restarts(0), mv_count(0) {}
};

// ---------------------------------------------------------------------------
// Internal helper: generate a deterministic or seeded starting vector
// ---------------------------------------------------------------------------
template <typename T>
std::vector<T> b_lanczos_start_vector(
    const std::size_t n, const unsigned int seed, const bool use_random)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    std::vector<T> v(n, T(0));
    unsigned int state = (use_random ? seed : 0u) + 1u;
    for (std::size_t i = 0; i < n; i++) {
        state = state * 1664525u + 1013904223u;
        const R val = R(static_cast<int>(state >> 16)) / R(32768) - R(1);
        v[i] = T(val);
    }
    const R nv = vcp::tsparse_scalar::real_norm_value(v);
    if (nv > R(0)) {
        for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
    }
    return v;
}

// ---------------------------------------------------------------------------
// b_inner_lanczos_eigs: B-inner product Lanczos for generalized symmetric SPD
//
// Solves:  A x = lambda B x
//   A: n x n symmetric sparse matrix
//   B: n x n SPD sparse matrix
//
// Basis is B-orthonormal:  V^T B V = I
// Projected matrix:        T_m = V^T A V  (standard inner product, dense)
// Ritz vector:             u = V y
// Generalized residual:    r = A u - lambda B u
//
// Returns b_lanczos_result<T> with eigenvalues sorted by target.
// Sets spd_check_failed=true and returns early on SPD violation.
// Sets budget_exhausted=true when mv_budget is reached.
//
// mv_budget: maximum total sparse matrix-vector products (A and B combined).
// norm_A, norm_B: Frobenius norms of A and B for relative residual denominator.
// ---------------------------------------------------------------------------
template <class SparseMatrix>
b_lanczos_result<typename SparseMatrix::value_type>
b_inner_lanczos_eigs(
    const SparseMatrix& A,
    const SparseMatrix& B,
    const std::size_t n,
    const std::size_t k,
    const std::size_t subspace_dim,
    const std::size_t mv_budget,
    const typename vcp::tsparse_scalar::real_type<typename SparseMatrix::value_type>::type& tol,
    const unsigned int random_seed,
    const bool random_start,
    const eig_target target,
    const typename vcp::tsparse_scalar::real_type<typename SparseMatrix::value_type>::type& norm_A,
    const typename vcp::tsparse_scalar::real_type<typename SparseMatrix::value_type>::type& norm_B,
    const bool compute_residual_history)
{
    typedef typename SparseMatrix::value_type T;
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    typedef b_lanczos_result<T> Result;

    Result res;

    // Edge cases
    if (n == 0) {
        res.failure_reason = "matrix is empty (n == 0)";
        return res;
    }
    if (k == 0) {
        res.converged = true;
        return res;
    }

    const std::size_t k_actual = (k < n) ? k : n;
    const std::size_t m_limit  = std::min(n, std::max(subspace_dim, k_actual + std::size_t(3)));
    const R eps   = std::numeric_limits<R>::epsilon();
    const R small_tol = tol / R(10);

    // Threshold for SPD violation detection
    const R spd_eps = eps * R(static_cast<int>(n) + 1);

    // Running total of ALL sparse matrix-vector products (A and B)
    std::size_t mv_total = 0;

    // Locked eigenpairs (converged across restarts)
    std::vector<T> locked_vals;
    std::vector<std::vector<T> > locked_vecs;
    std::vector<R> locked_res_abs;

    // Best (not yet converged) pairs for fallback
    std::vector<T> best_vals;
    std::vector<std::vector<T> > best_vecs;
    std::vector<R> best_res_abs;

    unsigned int seed_counter = random_seed;

    // Outer restart loop: runs until convergence, budget exhausted, or SPD failure.
    // No restart count limit — the mv_budget is the termination criterion.
    for (std::size_t restart = 0; locked_vals.size() < k_actual; restart++) {
        if (res.budget_exhausted) break;
        res.restarts = restart;

        // ---- Build B-orthonormal starting vector ----
        std::vector<T> v0;
        bool v0_ok = false;

        for (int attempt = 0; attempt < 5 && !v0_ok; attempt++) {
            v0 = b_lanczos_start_vector<T>(n, seed_counter, random_start);
            seed_counter++;

            // B-orthogonalize against locked vectors (two CGS passes)
            for (int pass = 0; pass < 2 && !locked_vecs.empty(); pass++) {
                if (mv_total >= mv_budget) { res.budget_exhausted = true; break; }
                const std::vector<T> Bv0 = B.mul_vec(v0);
                mv_total++;
                for (std::size_t li = 0; li < locked_vecs.size(); li++) {
                    R c = R(0);
                    for (std::size_t i = 0; i < n; i++) c += locked_vecs[li][i] * Bv0[i];
                    for (std::size_t i = 0; i < n; i++) v0[i] -= T(c) * locked_vecs[li][i];
                }
            }
            if (res.budget_exhausted) break;

            // B-normalize
            if (mv_total >= mv_budget) { res.budget_exhausted = true; break; }
            const std::vector<T> Bv0 = B.mul_vec(v0);
            mv_total++;
            R bsq = R(0);
            for (std::size_t i = 0; i < n; i++) bsq += v0[i] * Bv0[i];

            if (bsq <= spd_eps) {
                // Possibly SPD failure; try another seed
                if (attempt == 4) {
                    res.spd_check_failed = true;
                    res.failure_reason = "starting vector has non-positive B-norm; B may not be SPD";
                    res.mv_count = mv_total;
                    return res;
                }
                continue;
            }

            const R bn0 = vcp::tsparse_scalar::sqrt_value(bsq);
            for (std::size_t i = 0; i < n; i++) v0[i] /= T(bn0);
            v0_ok = true;
        }

        res.mv_count = mv_total;
        if (res.budget_exhausted) {
            res.failure_reason = "mv_budget exhausted during starting vector computation";
            break;
        }
        if (!v0_ok) {
            res.failure_reason = "exhausted B-orthogonal starting vectors";
            break;
        }

        // ---- Run B-inner Lanczos from v0 ----
        const std::size_t m_active = std::min(n - locked_vals.size(), m_limit);
        if (m_active == 0) break;

        std::vector<std::vector<T> > basis;
        basis.reserve(m_active + 1);
        basis.push_back(v0);

        // Tm_dense[i][j] = v_i^T A v_j  (standard inner product, dense m x m)
        // V^T B V = I (B-orthonormal), so T_m y = mu y gives generalized Ritz values.
        std::vector<std::vector<T> > Tm_dense(m_active, std::vector<T>(m_active, T(0)));

        std::size_t num_steps = 0;
        std::vector<T> v_prev(n, T(0));
        R beta_prev = R(0);
        std::string bd_reason;

        for (std::size_t j = 0; j < m_active; j++) {
            // Budget check before A.mul_vec
            if (mv_total >= mv_budget) { res.budget_exhausted = true; break; }

            // z = A * basis[j]  (w_j = A v_j before any modification)
            std::vector<T> z = A.mul_vec(basis[j]);
            mv_total++;

            // ---- Fill T_m column j using standard inner product v_i^T (A v_j) ----
            // Must be done BEFORE the three-term recurrence subtracts v_prev from z.
            for (std::size_t bi = 0; bi < basis.size(); bi++) {
                T dot_val = T(0);
                for (std::size_t i = 0; i < n; i++) dot_val += basis[bi][i] * z[i];
                Tm_dense[j][bi] = dot_val;
                Tm_dense[bi][j] = dot_val;   // symmetric
            }

            // Subtract previous beta term (Lanczos three-term recurrence)
            if (j > 0) {
                for (std::size_t i = 0; i < n; i++) z[i] -= T(beta_prev) * v_prev[i];
            }

            // ---- CGS Pass 1: B-reorthogonalization ----
            if (mv_total >= mv_budget) { res.budget_exhausted = true; break; }
            std::vector<T> Bz = B.mul_vec(z);
            mv_total++;

            for (std::size_t bi = 0; bi < basis.size(); bi++) {
                R c = R(0);
                for (std::size_t i = 0; i < n; i++) c += basis[bi][i] * Bz[i];
                for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * basis[bi][i];
            }
            for (std::size_t li = 0; li < locked_vecs.size(); li++) {
                R c = R(0);
                for (std::size_t i = 0; i < n; i++) c += locked_vecs[li][i] * Bz[i];
                for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * locked_vecs[li][i];
            }

            // ---- CGS Pass 2: second B-reorthogonalization for stability ----
            if (mv_total >= mv_budget) { res.budget_exhausted = true; break; }
            Bz = B.mul_vec(z);
            mv_total++;
            for (std::size_t bi = 0; bi < basis.size(); bi++) {
                R c = R(0);
                for (std::size_t i = 0; i < n; i++) c += basis[bi][i] * Bz[i];
                for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * basis[bi][i];
            }
            for (std::size_t li = 0; li < locked_vecs.size(); li++) {
                R c = R(0);
                for (std::size_t i = 0; i < n; i++) c += locked_vecs[li][i] * Bz[i];
                for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * locked_vecs[li][i];
            }

            // ---- Compute B-norm of z ----
            if (mv_total >= mv_budget) { res.budget_exhausted = true; break; }
            Bz = B.mul_vec(z);
            mv_total++;
            R beta_sq = R(0);
            for (std::size_t i = 0; i < n; i++) beta_sq += z[i] * Bz[i];

            // Check for SPD failure (x^T B x < 0 with significant margin)
            if (beta_sq < -spd_eps) {
                res.spd_check_failed = true;
                res.failure_reason = "B-inner product became negative; B may not be SPD";
                res.mv_count = mv_total;
                return res;
            }

            // Record per-step history (B-norm approximates local residual)
            if (compute_residual_history) {
                const R beta_approx = (beta_sq > R(0))
                    ? vcp::tsparse_scalar::sqrt_value(beta_sq) : R(0);
                res.history_abs.push_back(beta_approx);
                const T alpha_diag = Tm_dense[j][j];
                const R denom_h = R(1) + vcp::tsparse_scalar::abs_value(alpha_diag);
                res.history_rel.push_back(beta_approx / denom_h);
            }

            num_steps++;

            // Check for happy breakdown (B-norm too small to continue)
            if (beta_sq <= small_tol * small_tol) {
                bd_reason = "happy breakdown (zero B-norm residual)";
                break;
            }

            const R beta_j = vcp::tsparse_scalar::sqrt_value(beta_sq);
            v_prev = basis[j];
            beta_prev = beta_j;

            // Build B-normalized next basis vector
            std::vector<T> vnew(n);
            for (std::size_t i = 0; i < n; i++) vnew[i] = z[i] / T(beta_j);
            basis.push_back(vnew);
        }

        res.iterations += num_steps;
        res.mv_count = mv_total;

        if (res.budget_exhausted) {
            res.failure_reason = "mv_budget exhausted during Lanczos steps";
            break;
        }

        if (num_steps == 0) {
            res.breakdown_reason = "empty Lanczos basis";
            break;
        }

        // ---- Solve small dense symmetric eigenproblem T_m y = mu y ----
        // T_m = V^T A V  (standard projection, V B-orthonormal => V^T B V = I)
        // Eigenvalues mu of T_m approximate generalized eigenvalues lambda of Ax = lambda Bx.
        const std::size_t m = num_steps;
        std::vector<std::vector<T> > Tm(m, std::vector<T>(m, T(0)));
        for (std::size_t i = 0; i < m; i++) {
            for (std::size_t jj = 0; jj < m; jj++) {
                Tm[i][jj] = Tm_dense[i][jj];
            }
        }
        const std::size_t proj_max_iter = std::max(m * m * 100, std::size_t(1000));
        tsparse_dense_linalg::dense_eigen_result<T> small =
            tsparse_dense_linalg::jacobi_eig_dense(Tm, proj_max_iter, small_tol);

        if (small.eigenvalues.empty()) {
            res.breakdown_reason = "projected eigensolver failed";
            break;
        }

        // ---- Select target Ritz pairs ----
        std::vector<std::complex<R> > ceigs;
        ceigs.reserve(small.eigenvalues.size());
        for (std::size_t i = 0; i < small.eigenvalues.size(); i++) {
            ceigs.push_back(std::complex<R>(
                vcp::tsparse_scalar::real_part(small.eigenvalues[i]), R(0)));
        }
        const std::size_t k_remaining = k_actual - locked_vals.size();
        const std::vector<std::size_t> sel = tsparse_eigen_selection::select_eigen_indices(
            ceigs, std::min(m, k_remaining + m), target, R(0));

        // ---- Check convergence; lock one pair per restart ----
        bool ritz_budget_exhausted = false;
        for (std::size_t si = 0; si < sel.size() && locked_vals.size() < k_actual; si++) {
            const std::size_t idx = sel[si];
            if (idx >= small.eigenvectors.size()) continue;

            // Lift Ritz vector: u = V * y
            const std::vector<T>& y = small.eigenvectors[idx];
            std::vector<T> u(n, T(0));
            for (std::size_t ji = 0; ji < y.size() && ji < basis.size(); ji++) {
                for (std::size_t i = 0; i < n; i++) u[i] += basis[ji][i] * y[ji];
            }

            // B-orthogonalize u against locked vectors (two passes)
            bool orth_budget_ok = true;
            if (!locked_vecs.empty()) {
                for (int pass = 0; pass < 2 && orth_budget_ok; pass++) {
                    if (mv_total >= mv_budget) { orth_budget_ok = false; ritz_budget_exhausted = true; break; }
                    const std::vector<T> Bu_orth = B.mul_vec(u);
                    mv_total++;
                    for (std::size_t li = 0; li < locked_vecs.size(); li++) {
                        R c = R(0);
                        for (std::size_t i = 0; i < n; i++) c += locked_vecs[li][i] * Bu_orth[i];
                        for (std::size_t i = 0; i < n; i++) u[i] -= T(c) * locked_vecs[li][i];
                    }
                }
            }
            if (ritz_budget_exhausted) break;

            // B-norm of u (inline: avoids an extra B.mul_vec from b_norm helper)
            if (mv_total >= mv_budget) { ritz_budget_exhausted = true; break; }
            const std::vector<T> Bu_norm_vec = B.mul_vec(u);
            mv_total++;
            R bsq_u = R(0);
            for (std::size_t i = 0; i < n; i++) bsq_u += u[i] * Bu_norm_vec[i];
            if (bsq_u <= small_tol * small_tol) continue;
            const R u_bnorm = vcp::tsparse_scalar::sqrt_value(bsq_u);
            for (std::size_t i = 0; i < n; i++) u[i] /= T(u_bnorm);

            // Compute generalized residual:  r = A u - mu * B u
            const T mu = small.eigenvalues[idx];
            if (mv_total >= mv_budget) { ritz_budget_exhausted = true; break; }
            const std::vector<T> Au = A.mul_vec(u);
            mv_total++;
            if (mv_total >= mv_budget) { ritz_budget_exhausted = true; break; }
            const std::vector<T> Bu = B.mul_vec(u);
            mv_total++;

            std::vector<T> r(n);
            for (std::size_t i = 0; i < n; i++) r[i] = Au[i] - mu * Bu[i];

            const R u_l2norm = vcp::tsparse_scalar::real_norm_value(u);
            const R res_abs  = vcp::tsparse_scalar::real_norm_value(r);
            // Generalized relative residual denominator:
            //   ||A||_F * ||u|| + |mu| * ||B||_F * ||u|| + eps
            const R denom    = norm_A * u_l2norm
                + vcp::tsparse_scalar::abs_value(mu) * norm_B * u_l2norm
                + std::numeric_limits<R>::epsilon();
            const R res_rel  = res_abs / denom;

            if (compute_residual_history) {
                res.history_abs.push_back(res_abs);
                res.history_rel.push_back(res_rel);
            }

            // Keep best unconverged pairs for fallback
            if (best_vals.size() < k_actual) {
                best_vals.push_back(mu);
                best_vecs.push_back(u);
                best_res_abs.push_back(res_abs);
            }

            // Convergence check
            if (res_abs <= tol || res_rel <= tol) {
                locked_vals.push_back(mu);
                locked_vecs.push_back(u);
                locked_res_abs.push_back(res_abs);
                res.converged_count++;
                break;  // One lock per restart; next restart deflates this pair
            }
        }

        res.mv_count = mv_total;
        if (ritz_budget_exhausted) {
            res.budget_exhausted = true;
            res.failure_reason = "mv_budget exhausted during Ritz pair evaluation";
            break;
        }

        if (locked_vals.size() >= k_actual) {
            res.converged = true;
            res.breakdown_reason = bd_reason;
            break;
        }
    }

    // ---- Sort and assemble locked eigenpairs ----
    if (!locked_vals.empty()) {
        std::vector<std::complex<R> > ceigs_locked;
        ceigs_locked.reserve(locked_vals.size());
        for (std::size_t i = 0; i < locked_vals.size(); i++) {
            ceigs_locked.push_back(std::complex<R>(
                vcp::tsparse_scalar::real_part(locked_vals[i]), R(0)));
        }
        const std::size_t ksorted = std::min(k_actual, locked_vals.size());
        const std::vector<std::size_t> order = tsparse_eigen_selection::select_eigen_indices(
            ceigs_locked, ksorted, target, R(0));

        for (std::size_t i = 0; i < order.size(); i++) {
            res.eigenvalues.push_back(locked_vals[order[i]]);
            res.eigenvectors.push_back(locked_vecs[order[i]]);
            res.residuals_abs.push_back(locked_res_abs[order[i]]);
        }
        res.residuals_rel.resize(res.residuals_abs.size());
        for (std::size_t i = 0; i < res.eigenvalues.size(); i++) {
            const R u_l2norm = vcp::tsparse_scalar::real_norm_value(res.eigenvectors[i]);
            const R denom    = norm_A * u_l2norm
                + vcp::tsparse_scalar::abs_value(res.eigenvalues[i]) * norm_B * u_l2norm
                + std::numeric_limits<R>::epsilon();
            res.residuals_rel[i] = res.residuals_abs[i] / denom;
        }
        res.returned_count = res.eigenvalues.size();
    }

    if (!res.converged && res.failure_reason.empty()) {
        res.failure_reason = "maximum restarts reached without full convergence";
    }

    // ---- Fallback: return best unconverged pairs if no locked pairs ----
    if (res.eigenvalues.empty() && !best_vals.empty()) {
        std::vector<std::complex<R> > ceigs_best;
        ceigs_best.reserve(best_vals.size());
        for (std::size_t i = 0; i < best_vals.size(); i++) {
            ceigs_best.push_back(std::complex<R>(
                vcp::tsparse_scalar::real_part(best_vals[i]), R(0)));
        }
        const std::vector<std::size_t> order = tsparse_eigen_selection::select_eigen_indices(
            ceigs_best, std::min(k_actual, best_vals.size()), target, R(0));
        for (std::size_t i = 0; i < order.size(); i++) {
            res.eigenvalues.push_back(best_vals[order[i]]);
            res.eigenvectors.push_back(best_vecs[order[i]]);
            res.residuals_abs.push_back(best_res_abs[order[i]]);
        }
        res.residuals_rel.resize(res.residuals_abs.size());
        for (std::size_t i = 0; i < res.eigenvalues.size(); i++) {
            const R u_l2norm = vcp::tsparse_scalar::real_norm_value(res.eigenvectors[i]);
            const R denom    = norm_A * u_l2norm
                + vcp::tsparse_scalar::abs_value(res.eigenvalues[i]) * norm_B * u_l2norm
                + std::numeric_limits<R>::epsilon();
            res.residuals_rel[i] = res.residuals_abs[i] / denom;
        }
        res.returned_count = res.eigenvalues.size();
    }

    return res;
}

} // namespace tsparse_b_inner_lanczos
} // namespace vcp

#endif // VCP_TSPARSE_B_INNER_LANCZOS_HPP
