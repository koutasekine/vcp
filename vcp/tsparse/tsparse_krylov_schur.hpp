// VCP Library
// http://verified.computation.jp
//
// vcp/tsparse/tsparse_krylov_schur.hpp
//
// Experimental Krylov-Schur Arnoldi eigensolver for REAL GENERAL operators.
// Phase 3 of the spmatrix next-stage design.
//
// Namespace: vcp::tsparse_experimental
//
// NOT connected to vcp::spmatrix::eigs / eig_method enum.
// NOT connected to the existing spmatrix::eigs dispatch.
// Complex types T are rejected via static_assert.
//
// Apply functor convention:
//   void operator()(const std::vector<T>& x, std::vector<T>& y) const;
//
// options.max_iter = total matrix-vector product budget (main loop).
// Final residual re-evaluation apply() calls are counted in
// matrix_vector_products but are NOT subject to the budget check.
//
// Design notes:
//   - Arnoldi basis is orthonormal; full reorthogonalization is performed.
//   - Restart uses approximate Ritz vectors (via inverse iteration on the
//     projected Hessenberg) to form the retained subspace.  This is a
//     Ritz-based variant of Krylov-Schur; it provides the same invariant-
//     subspace property when Ritz vectors are well-converged.
//   - Complex Ritz pairs are tracked via complex_eigenvalues but are NOT
//     used as restart vectors in Phase 3.  This limitation is documented
//     in the Phase 3 spec (tsparse_krylov_schur complex pair restart
//     limitation).
//   - The projected Hessenberg after restart may be a dense p x p block.
//     It is reduced to upper Hessenberg form before eigenvalue extraction.
//   - n==0, k==0, k>n are handled per spec.
//   - used_dense_fallback / used_shift_invert / used_generalized_operator
//     are always false.

#pragma once

#ifndef VCP_TSPARSE_KRYLOV_SCHUR_HPP
#define VCP_TSPARSE_KRYLOV_SCHUR_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/spmatrix.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>
#include <vcp/tsparse/tsparse_projected_hessenberg.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_experimental {

// ---------------------------------------------------------------------------
// Type guard: only real floating-point T is supported.
// ---------------------------------------------------------------------------
template <class T>
struct ks_is_real_floating {
    static const bool value =
        std::is_floating_point<T>::value &&
        !vcp::tsparse_scalar::is_complex<T>::value;
};

// ---------------------------------------------------------------------------
// Diagnostic result wrapper
// ---------------------------------------------------------------------------
template <class T>
struct krylov_schur_result {
    vcp::eig_result<T> eigs;
    std::size_t restart_count;           // Krylov-Schur restarts performed
    std::size_t matrix_vector_products;  // total apply() calls (including final eval)
    std::size_t locked_real_count;       // converged real Ritz pairs returned
    std::size_t locked_complex_count;    // conjugate complex pairs in complex_eigenvalues
    std::size_t projected_dimension;     // effective Krylov subspace dimension m

    krylov_schur_result()
        : restart_count(0), matrix_vector_products(0),
          locked_real_count(0), locked_complex_count(0), projected_dimension(0) {}
};

// ===========================================================================
// Internal helpers
// ===========================================================================
namespace ks_detail {

// Effective subspace dimension m.  Always satisfies k <= m <= n.
inline std::size_t effective_m(std::size_t n, std::size_t k, std::size_t requested)
{
    if (n == 0 || k == 0) return 0;
    std::size_t m;
    if (requested == 0) {
        m = std::max<std::size_t>(
            std::max<std::size_t>(2 * k + 6, std::size_t(20)), k + 3);
    } else if (requested <= k) {
        // subspace_dim too small: pad to k + max(requested, 2)
        m = k + std::max<std::size_t>(requested, std::size_t(2));
    } else {
        m = requested;
    }
    if (m > n) m = n;
    if (k < n && m < k + 1) m = std::min(k + 1, n);
    return m;
}

// Orthogonalize w against V[0..j] and accumulate Gram-Schmidt coefficients h.
template <typename T>
void arnoldi_orthogonalize(
    const std::vector<std::vector<T> >& V,
    const std::size_t j,          // last basis index (inclusive)
    std::vector<T>& w,
    std::vector<T>& h,            // output: length j+1
    const vcp::orthogonalization_method orth)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = w.size();
    h.assign(j + 1, T(0));

    const bool twice =
        (orth == vcp::orthogonalization_method::classical_gram_schmidt_twice);
    const int passes = twice ? 2 : 1;

    for (int pass = 0; pass < passes; pass++) {
        for (std::size_t i = 0; i <= j; i++) {
            const R c = vcp::tsparse_scalar::real_dot_value(V[i], w);
            h[i] += T(c);
            for (std::size_t ii = 0; ii < n; ii++) w[ii] -= T(c) * V[i][ii];
        }
    }
}

// Extract the square m x m matrix H_m from H_mat.
template <typename T>
std::vector<std::vector<T> > extract_hm(
    const std::vector<std::vector<T> >& H_mat,
    const std::size_t m)
{
    std::vector<std::vector<T> > Hm(m, std::vector<T>(m, T(0)));
    for (std::size_t i = 0; i < m; i++)
        for (std::size_t j = 0; j < m; j++)
            Hm[i][j] = H_mat[i][j];
    return Hm;
}

// Lift a projected Ritz vector y (length m) to full space via V_basis.
template <typename T>
std::vector<T> lift_ritz(
    const std::vector<std::vector<T> >& V,
    const std::vector<T>& y,
    const std::size_t n,
    const std::size_t m)
{
    std::vector<T> u(n, T(0));
    const std::size_t len = std::min(y.size(), std::min(m, V.size()));
    for (std::size_t j = 0; j < len; j++)
        for (std::size_t i = 0; i < n; i++)
            u[i] += V[j][i] * y[j];
    return u;
}

// Compute H_proj = Q_p^T * Hm * Q_p  (p_keep x p_keep).
// Q_p is stored as p_keep vectors each of length m.
template <typename T>
std::vector<std::vector<T> > project_hessenberg(
    const std::vector<std::vector<T> >& Hm,   // m x m
    const std::vector<std::vector<T> >& Qp)   // p_keep vectors of length m
{
    const std::size_t m  = Hm.size();
    const std::size_t p  = Qp.size();

    // tmp[i][j] = (Hm * Qp[j])[i]
    std::vector<std::vector<T> > tmp(m, std::vector<T>(p, T(0)));
    for (std::size_t i = 0; i < m; i++) {
        for (std::size_t j = 0; j < p; j++) {
            T s = T(0);
            const std::size_t cols = std::min(m, Hm[i].size());
            for (std::size_t kk = 0; kk < cols; kk++)
                s += Hm[i][kk] * Qp[j][kk];
            tmp[i][j] = s;
        }
    }

    // H_proj[i][j] = Qp[i]^T * tmp[:,j]
    std::vector<std::vector<T> > Hp(p, std::vector<T>(p, T(0)));
    for (std::size_t i = 0; i < p; i++) {
        for (std::size_t j = 0; j < p; j++) {
            T s = T(0);
            for (std::size_t kk = 0; kk < m; kk++)
                s += Qp[i][kk] * tmp[kk][j];
            Hp[i][j] = s;
        }
    }
    return Hp;
}

// Orthogonalize a set of m-vectors in-place (CGS2 within the set).
// Returns the number of linearly independent vectors retained.
template <typename T>
std::size_t orthonormalize_set(
    std::vector<std::vector<T> >& vecs,
    const std::size_t m,
    const typename vcp::tsparse_scalar::real_type<T>::type& tol)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    std::size_t kept = 0;

    for (std::size_t j = 0; j < vecs.size(); j++) {
        std::vector<T>& v = vecs[j];
        if (v.size() != m) v.assign(m, T(0));

        for (int pass = 0; pass < 2; pass++) {
            for (std::size_t i = 0; i < kept; i++) {
                R c = R(0);
                for (std::size_t kk = 0; kk < m; kk++)
                    c += vcp::tsparse_scalar::real_part(vecs[i][kk] * v[kk]);
                for (std::size_t kk = 0; kk < m; kk++)
                    v[kk] -= T(c) * vecs[i][kk];
            }
        }

        R nrm = vcp::tsparse_scalar::real_norm_value(v);
        if (nrm <= tol) continue;
        for (std::size_t kk = 0; kk < m; kk++) v[kk] /= T(nrm);

        if (j != kept) vecs[kept] = vecs[j];
        kept++;
    }
    vecs.resize(kept);
    return kept;
}

// Check whether a complex eigenvalue is effectively real.
template <typename R>
bool is_effectively_real(const std::complex<R>& z, const R& tol)
{
    return vcp::tsparse_scalar::abs_value(z.imag()) <=
           tol * (vcp::tsparse_scalar::abs_value(z.real()) + R(1));
}

} // namespace ks_detail

// ===========================================================================
// krylov_schur_eigs_with_diagnostics
// ===========================================================================
template <class Apply, class T>
krylov_schur_result<T> krylov_schur_eigs_with_diagnostics(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    static_assert(
        ks_is_real_floating<T>::value,
        "Krylov-Schur Arnoldi experimental solver supports real scalar types only in Phase 3");

    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    typedef std::complex<R> C;

    krylov_schur_result<T> diag;
    vcp::eig_result<T>& result = diag.eigs;

    result.used_dense_fallback       = false;
    result.used_shift_invert         = false;
    result.used_generalized_operator = false;
    result.used_method               = "krylov_schur_experimental";
    result.used_orthogonalization    =
        (options.orthogonalization ==
         vcp::orthogonalization_method::classical_gram_schmidt_twice)
        ? "classical_gram_schmidt_twice"
        : "modified_gram_schmidt";

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------
    if (k == 0) {
        result.requested_count = 0;
        result.returned_count  = 0;
        result.converged       = true;
        result.status          = "success";
        result.message         = "k=0: nothing to compute";
        return diag;
    }
    if (n == 0) {
        result.requested_count = k;
        result.returned_count  = 0;
        result.converged       = false;
        result.status          = "failed";
        result.failure_reason  = "dimension n=0 with k>0";
        return diag;
    }

    const std::size_t k_original = k;
    if (k > n) k = n;
    result.requested_count = k_original;

    // -----------------------------------------------------------------------
    // Parameters
    // -----------------------------------------------------------------------
    const R tol = (options.tol > R(0))
        ? options.tol
        : vcp::tsparse_scalar::decimal_power_negative<R>(12);

    const std::size_t max_mv = (options.max_iter > 0)
        ? options.max_iter
        : std::size_t(300) * n;

    const std::size_t m_limit =
        ks_detail::effective_m(n, k, options.subspace_dim);

    result.used_subspace_dim = m_limit;
    diag.projected_dimension  = m_limit;

    const vcp::eig_target target   = options.target;
    const R shift_val              = options.shift;
    const bool compute_hist        = options.compute_residual_history;
    const vcp::orthogonalization_method orth = options.orthogonalization;

    unsigned int seed = options.random_start ? options.random_seed : 0u;

    const R small_tol = std::max(
        tol * R(1e-4),
        std::numeric_limits<R>::epsilon() * R(n + 1));

    // -----------------------------------------------------------------------
    // State
    // -----------------------------------------------------------------------
    // V[j] = j-th Arnoldi basis vector (length n).
    std::vector<std::vector<T> > V;
    V.reserve(m_limit + 2);

    // H_mat[row][col]: Hessenberg matrix, shape (m_limit+1) x m_limit.
    std::vector<std::vector<T> > H_mat(m_limit + 1, std::vector<T>(m_limit, T(0)));

    std::size_t m_current = 0;   // completed Arnoldi steps
    std::size_t mv_count  = 0;
    std::size_t restart_count = 0;
    bool budget_exhausted = false;
    bool happy_breakdown  = false;

    // Converged real Ritz pairs.
    struct real_ritz_t {
        T value;
        std::vector<T> vector;
        R res_abs;
        R res_rel;
    };
    std::vector<real_ritz_t> locked;
    locked.reserve(k + 1);

    // -----------------------------------------------------------------------
    // Initial vector
    // -----------------------------------------------------------------------
    V.push_back(vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++));

    // -----------------------------------------------------------------------
    // Main loop
    // -----------------------------------------------------------------------
    while (locked.size() < k && mv_count < max_mv) {

        // -------------------------------------------------------------------
        // PHASE 1: Expand Arnoldi from m_current to m_limit.
        // -------------------------------------------------------------------
        happy_breakdown = false;

        while (m_current < m_limit && mv_count < max_mv) {
            const std::size_t j = m_current;

            std::vector<T> w;
            apply(V[j], w);
            mv_count++;

            std::vector<T> h;
            ks_detail::arnoldi_orthogonalize(V, j, w, h, orth);

            const R beta = vcp::tsparse_scalar::real_norm_value(w);

            for (std::size_t i = 0; i <= j; i++) H_mat[i][j] = h[i];
            H_mat[j + 1][j] = T(beta);

            m_current = j + 1;

            if (beta <= small_tol) {
                // When happy breakdown occurs before filling the subspace, try
                // extending into the orthogonal complement of the current basis.
                // This handles repeated eigenvalues: a start vector with equal
                // components in both copies of an eigenspace collapses them into
                // one Krylov direction; a complement vector breaks that symmetry
                // and lets the second copy be discovered in the same pass.
                bool extended = false;
                if (m_current < m_limit) {
                    std::vector<T> v_ext =
                        vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
                    for (int xpass = 0; xpass < 2; xpass++) {
                        for (std::size_t jj = 0; jj < m_current; jj++) {
                            const R xc = vcp::tsparse_scalar::real_dot_value(
                                V[jj], v_ext);
                            for (std::size_t ii = 0; ii < n; ii++)
                                v_ext[ii] -= T(xc) * V[jj][ii];
                        }
                        for (std::size_t li = 0; li < locked.size(); li++) {
                            const R xc = vcp::tsparse_scalar::real_dot_value(
                                locked[li].vector, v_ext);
                            for (std::size_t ii = 0; ii < n; ii++)
                                v_ext[ii] -= T(xc) * locked[li].vector[ii];
                        }
                    }
                    const R vxn = vcp::tsparse_scalar::real_norm_value(v_ext);
                    if (vxn > small_tol) {
                        for (std::size_t ii = 0; ii < n; ii++) v_ext[ii] /= T(vxn);
                        if (V.size() <= m_current) V.push_back(v_ext);
                        else                       V[m_current] = v_ext;
                        extended = true;
                    }
                }
                if (extended) continue;
                happy_breakdown = true;
                std::vector<T> dummy(n, T(0));
                if (V.size() <= m_current) V.push_back(dummy);
                else                       V[m_current] = dummy;
                break;
            }

            std::vector<T> vnew(n);
            for (std::size_t i = 0; i < n; i++) vnew[i] = w[i] / T(beta);
            if (V.size() <= m_current) V.push_back(vnew);
            else                       V[m_current] = vnew;
        }

        const std::size_t m = m_current;
        if (m == 0) { result.failure_reason = "Arnoldi produced empty basis"; break; }

        // -------------------------------------------------------------------
        // PHASE 2: Extract eigenvalues of H[0:m][0:m].
        //
        // Reduce Hm to upper Hessenberg form via Householder, accumulating the
        // orthogonal transformation Q (so that Q^T * Hm * Q = Hm_hess).
        // Eigenvalues are extracted from Hm_hess.
        // Eigenvectors of Hm_hess are then mapped back to the Krylov basis via Q.
        // -------------------------------------------------------------------
        std::vector<std::vector<T> > Hm      = ks_detail::extract_hm(H_mat, m);
        std::vector<std::vector<T> > Hm_hess = Hm;  // modified in-place
        std::vector<std::vector<T> > Hess_Q;         // Hm = Hm_hess after Q-transform
        vcp::tsparse_proj_hess::reduce_to_hessenberg_with_q(Hm_hess, Hess_Q);

        const std::size_t hess_iter = m * m * 80 + 300;
        const R hess_tol = small_tol;

        std::vector<C> all_eigs =
            vcp::tsparse_eigensolvers::hessenberg_complex_eigenvalues<T>(
                Hm_hess, hess_iter, hess_tol);

        if (all_eigs.empty()) {
            result.failure_reason = "projected Hessenberg eigensolver returned no eigenvalues";
            break;
        }

        // Separate real and complex eigenvalues (positive-imag complex pairs stored once).
        std::vector<R> real_eig_vals;
        std::vector<C> complex_eig_vals;
        const R real_tol = hess_tol * R(100);

        for (std::size_t ei = 0; ei < all_eigs.size(); ei++) {
            if (ks_detail::is_effectively_real(all_eigs[ei], real_tol)) {
                real_eig_vals.push_back(all_eigs[ei].real());
            } else if (all_eigs[ei].imag() > R(0)) {
                complex_eig_vals.push_back(all_eigs[ei]);
            }
        }

        // -------------------------------------------------------------------
        // PHASE 3: Compute Ritz residuals for real eigenvalues.
        //
        // Inverse iteration is performed on Hm_hess (not Hm) for consistency:
        // eigenvalues came from Hm_hess.  The resulting eigenvector y_h is in
        // the Hm_hess basis.  To get back to the Krylov (V) basis:
        //   y_krylov = Hess_Q * y_h   (since Q^T * Hm * Q = Hm_hess)
        // The Arnoldi residual bound uses the last element of y_krylov.
        //
        // For repeated Ritz values, inverse iteration from [1,...,1] always
        // converges to the same direction.  We deflate each y_h against
        // previously extracted eigenvectors in the Hm_hess basis; when the
        // result is near-zero, we re-run inverse iteration from a fresh start
        // vector orthogonal to all previous ones, so that the independent
        // direction of the repeated eigenspace is found without restarting.
        // -------------------------------------------------------------------
        const R beta_overflow = vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(H_mat[m][m - 1]));

        struct ritz_candidate_t {
            R value;
            std::vector<T> y_proj;   // Ritz vector in Krylov basis (length m)
            R res_bound;             // Arnoldi residual bound
            bool converged;
        };
        std::vector<ritz_candidate_t> candidates;
        candidates.reserve(real_eig_vals.size());

        // Track normalized eigenvectors of Hm_hess extracted in this sweep.
        // Used to detect repeated Ritz values: when inverse iteration from
        // [1,...,1] converges to the same direction for the same eigenvalue,
        // we regenerate with a start vector orthogonal to all prior directions
        // so that the independent eigenspace direction is found.
        // For distinct eigenvalues the eigenvectors differ by shift, so the
        // near-parallelism check does not fire and y_h is used unmodified.
        std::vector<std::vector<T> > prev_y_h_vecs;
        prev_y_h_vecs.reserve(real_eig_vals.size());
        const R parallel_thresh = R(1) - std::sqrt(small_tol);

        for (std::size_t ei = 0; ei < real_eig_vals.size(); ei++) {
            const R theta = real_eig_vals[ei];

            // Eigenvector of Hm_hess for eigenvalue theta (standard start).
            std::vector<T> y_h =
                vcp::tsparse_dense_linalg::dense_eigenvector_inverse_iteration(
                    Hm_hess, T(theta));
            if (y_h.empty() || y_h.size() != m) continue;
            const R yhn0 = vcp::tsparse_scalar::real_norm_value(y_h);
            if (yhn0 <= small_tol) continue;
            for (std::size_t i = 0; i < m; i++) y_h[i] /= T(yhn0);

            // Check whether y_h is nearly parallel to any already-extracted
            // eigenvector.  For a repeated Ritz value the same shift produces
            // the same converged direction, giving |dot| close to 1; for a
            // distinct eigenvalue the directions differ and the check is safe.
            bool need_fresh = false;
            for (std::size_t ji = 0; ji < prev_y_h_vecs.size(); ji++) {
                R c = R(0);
                for (std::size_t ii = 0; ii < m; ii++)
                    c += vcp::tsparse_scalar::real_part(
                        prev_y_h_vecs[ji][ii] * y_h[ii]);
                if (c > parallel_thresh || c < -parallel_thresh) {
                    need_fresh = true;
                    break;
                }
            }

            if (need_fresh) {
                // Build a fresh start vector orthogonal ONLY to the parallel
                // prev vectors (those for the same repeated eigenvalue).
                // Deflating against all prev vectors can collapse the start
                // below small_tol when the subspace is nearly full.
                std::vector<T> y_h_start(m, T(1));
                for (int pass = 0; pass < 2; pass++) {
                    for (std::size_t ji = 0; ji < prev_y_h_vecs.size(); ji++) {
                        R cp = R(0);
                        for (std::size_t ii = 0; ii < m; ii++)
                            cp += vcp::tsparse_scalar::real_part(
                                prev_y_h_vecs[ji][ii] * y_h[ii]);
                        if (cp <= parallel_thresh && cp >= -parallel_thresh) continue;
                        R c = R(0);
                        for (std::size_t ii = 0; ii < m; ii++)
                            c += vcp::tsparse_scalar::real_part(
                                prev_y_h_vecs[ji][ii] * y_h_start[ii]);
                        for (std::size_t ii = 0; ii < m; ii++)
                            y_h_start[ii] -= T(c) * prev_y_h_vecs[ji][ii];
                    }
                }
                const R y_h_start_n = vcp::tsparse_scalar::real_norm_value(y_h_start);
                if (y_h_start_n <= small_tol) continue;
                for (std::size_t ii = 0; ii < m; ii++)
                    y_h_start[ii] /= T(y_h_start_n);

                y_h = vcp::tsparse_dense_linalg::dense_eigenvector_inverse_iteration_from(
                    Hm_hess, T(theta), y_h_start);
                if (y_h.empty() || y_h.size() != m) continue;
                const R yhn_fresh = vcp::tsparse_scalar::real_norm_value(y_h);
                if (yhn_fresh <= small_tol) continue;
                for (std::size_t i = 0; i < m; i++) y_h[i] /= T(yhn_fresh);
            }
            prev_y_h_vecs.push_back(y_h);

            // Map to Krylov basis: y = Q * y_h  (Q = Hess_Q).
            std::vector<T> y(m, T(0));
            for (std::size_t i = 0; i < m; i++)
                for (std::size_t j = 0; j < m; j++)
                    y[i] += Hess_Q[i][j] * y_h[j];

            const R yn = vcp::tsparse_scalar::real_norm_value(y);
            if (yn <= small_tol) continue;
            for (std::size_t i = 0; i < m; i++) y[i] /= T(yn);

            // Arnoldi residual bound: |h_{m+1,m}| * |y[m-1]|.
            const R res_bound =
                beta_overflow * vcp::tsparse_scalar::abs_value(y.back());

            ritz_candidate_t rc;
            rc.value     = theta;
            rc.y_proj    = y;
            rc.res_bound = res_bound;
            rc.converged = false;
            candidates.push_back(rc);
        }

        // Sort candidates by target.
        if (!candidates.empty()) {
            std::vector<C> cvals;
            cvals.reserve(candidates.size());
            for (std::size_t i = 0; i < candidates.size(); i++)
                cvals.push_back(C(candidates[i].value, R(0)));
            const std::vector<std::size_t> sel =
                vcp::tsparse_eigensolvers::select_ritz_indices<T>(
                    cvals, candidates.size(), target, shift_val);
            std::vector<ritz_candidate_t> sorted;
            sorted.reserve(sel.size());
            for (std::size_t i = 0; i < sel.size(); i++)
                sorted.push_back(candidates[sel[i]]);
            candidates.swap(sorted);
        }

        // Check convergence: compute exact residual in target-priority order.
        // Only lock candidates from the top of the sorted list.  If a higher-
        // priority Ritz value does not converge (poor Ritz vector approximation),
        // we break immediately and restart so that the subspace improves for that
        // target eigenvalue rather than locking out-of-target ones.
        //
        // Collinearity skip: if the lifted Ritz vector is nearly collinear with
        // an already-locked eigenvector (|cos angle| near 1), it is a re-discovery
        // of a direction already in the locked set.  We continue to the next
        // candidate rather than breaking, so that genuinely new eigenvalues in
        // the same iteration can still be locked.  For repeated eigenvalues the
        // Phase 3 extraction deflation above has already produced independent
        // directions, so the collinearity check fires only for re-discovered
        // distinct locked eigenvalues (e.g. after restart the locked eigenvalue
        // appears again as the first candidate).
        for (std::size_t ci = 0; ci < candidates.size() && locked.size() < k; ci++) {
            ritz_candidate_t& rc = candidates[ci];

            const bool candidate_ok =
                (rc.res_bound <= tol * R(100)) || happy_breakdown;
            if (!candidate_ok) break;  // higher-priority candidate not ready: restart

            // Lift to full space.
            std::vector<T> vfull = ks_detail::lift_ritz(V, rc.y_proj, n, m);
            const R vn = vcp::tsparse_scalar::real_norm_value(vfull);
            if (vn <= small_tol) { break; }
            for (std::size_t i = 0; i < n; i++) vfull[i] /= T(vn);

            // Exact residual: r = A v - theta v.
            std::vector<T> Av;
            apply(vfull, Av);
            mv_count++;

            const T mu = T(rc.value);
            std::vector<T> r(n);
            for (std::size_t i = 0; i < n; i++) r[i] = Av[i] - mu * vfull[i];
            const R res_abs = vcp::tsparse_scalar::real_norm_value(r);
            const R res_rel = res_abs / (R(1) + vcp::tsparse_scalar::abs_value(rc.value));

            if (res_abs <= tol || res_rel <= tol) {
                // Check whether this Ritz vector is collinear with an already-
                // locked eigenvector (|cos angle| > 1 - eps_col).  If so, it
                // is a re-discovery of a locked direction (can happen after
                // restart) and we skip it rather than duplicating the entry.
                const R eps_col = small_tol * R(100);
                bool collinear = false;
                for (std::size_t li = 0; li < locked.size(); li++) {
                    const R c = vcp::tsparse_scalar::real_dot_value(
                        locked[li].vector, vfull);
                    if (c > R(1) - eps_col || c < -(R(1) - eps_col)) {
                        collinear = true;
                        break;
                    }
                }
                if (collinear) { continue; }

                real_ritz_t rr;
                rr.value   = T(rc.value);
                rr.vector  = vfull;
                rr.res_abs = res_abs;
                rr.res_rel = res_rel;
                locked.push_back(rr);
                rc.converged = true;
            } else {
                // Target eigenvalue did not pass exact residual: restart.
                break;
            }
        }

        // -------------------------------------------------------------------
        // PHASE 4: Residual history.
        // -------------------------------------------------------------------
        if (compute_hist && !candidates.empty()) {
            R best_abs   = candidates[0].res_bound;
            R best_theta = vcp::tsparse_scalar::abs_value(candidates[0].value);
            for (std::size_t ci = 1; ci < candidates.size(); ci++) {
                if (candidates[ci].res_bound < best_abs) {
                    best_abs   = candidates[ci].res_bound;
                    best_theta = vcp::tsparse_scalar::abs_value(candidates[ci].value);
                }
            }
            result.residual_history_absolute.push_back(best_abs);
            result.residual_history_relative.push_back(best_abs / (R(1) + best_theta));
        }

        if (locked.size() >= k) break;

        if (mv_count >= max_mv) { budget_exhausted = true; break; }

        // If happy breakdown and no new convergence, signal failure.
        if (happy_breakdown && locked.size() < k && candidates.empty()) {
            result.breakdown_reason = "happy breakdown: Krylov subspace exhausted";
            break;
        }

        // -------------------------------------------------------------------
        // PHASE 5: Krylov-Schur restart.
        // -------------------------------------------------------------------
        const std::size_t padding     = std::max<std::size_t>(k / 2 + 1, std::size_t(2));
        const std::size_t p_keep_max  = (m > 1) ? (m - 1) : 0;
        std::size_t p_keep = std::min(k + padding, p_keep_max);

        if (p_keep == 0) {
            // Degenerate: full fresh restart deflated against locked vectors.
            for (std::size_t i = 0; i <= m_limit; i++)
                std::fill(H_mat[i].begin(), H_mat[i].end(), T(0));
            std::vector<T> v0 = vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
            for (int pass = 0; pass < 2; pass++)
                for (std::size_t li = 0; li < locked.size(); li++) {
                    const R lc = vcp::tsparse_scalar::real_dot_value(locked[li].vector, v0);
                    for (std::size_t ii = 0; ii < n; ii++) v0[ii] -= T(lc) * locked[li].vector[ii];
                }
            const R v0n = vcp::tsparse_scalar::real_norm_value(v0);
            if (v0n > small_tol) for (std::size_t ii = 0; ii < n; ii++) v0[ii] /= T(v0n);
            V.clear();
            V.push_back(v0);
            m_current = 0;
            restart_count++;
            continue;
        }

        // Collect up to p_keep Ritz vectors in projected space.
        // Exclude already-locked candidates so the restart subspace focuses
        // on the remaining target eigenvalues and does not re-introduce
        // already-converged directions (which would persist in the Hessenberg
        // coupling and slow convergence of subsequent eigenvalues).
        std::vector<std::vector<T> > Qp;
        Qp.reserve(p_keep);

        for (std::size_t ci = 0; ci < candidates.size() && Qp.size() < p_keep; ci++) {
            if (!candidates[ci].converged)
                Qp.push_back(candidates[ci].y_proj);
        }

        // Pad with random directions in R^m if not enough real Ritz vectors.
        while (Qp.size() < p_keep) {
            std::vector<T> rv(m);
            unsigned int s2 = seed++;
            for (std::size_t i = 0; i < m; i++) {
                s2 = s2 * 1664525u + 1013904223u;
                rv[i] = T(static_cast<int>(s2 >> 16)) / T(32768) - T(1);
            }
            Qp.push_back(rv);
        }

        // Orthonormalize Qp in R^m.
        const std::size_t n_kept =
            ks_detail::orthonormalize_set(Qp, m, small_tol);

        if (n_kept == 0) {
            // Fallback: fresh restart deflated against locked vectors.
            for (std::size_t i = 0; i <= m_limit; i++)
                std::fill(H_mat[i].begin(), H_mat[i].end(), T(0));
            std::vector<T> v0b = vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
            for (int pass = 0; pass < 2; pass++)
                for (std::size_t li = 0; li < locked.size(); li++) {
                    const R lc = vcp::tsparse_scalar::real_dot_value(locked[li].vector, v0b);
                    for (std::size_t ii = 0; ii < n; ii++) v0b[ii] -= T(lc) * locked[li].vector[ii];
                }
            const R v0bn = vcp::tsparse_scalar::real_norm_value(v0b);
            if (v0bn > small_tol) for (std::size_t ii = 0; ii < n; ii++) v0b[ii] /= T(v0bn);
            V.clear();
            V.push_back(v0b);
            m_current = 0;
            restart_count++;
            continue;
        }
        p_keep = n_kept;

        // Compute H_proj = Qp^T Hm Qp  (p_keep x p_keep dense restart block).
        const std::vector<std::vector<T> > H_proj =
            ks_detail::project_hessenberg(Hm, Qp);

        // Coupling row: c[j] = beta_overflow * Qp[j][m-1]
        std::vector<T> coupling(p_keep);
        for (std::size_t j = 0; j < p_keep; j++)
            coupling[j] = T(beta_overflow) * Qp[j][m - 1];

        // Lift Ritz vectors to full space.
        std::vector<std::vector<T> > V_new(p_keep, std::vector<T>(n, T(0)));
        for (std::size_t j = 0; j < p_keep; j++) {
            for (std::size_t i = 0; i < m && i < V.size(); i++)
                for (std::size_t ii = 0; ii < n; ii++)
                    V_new[j][ii] += V[i][ii] * Qp[j][i];
            // Re-normalize for floating-point hygiene.
            const R nrm = vcp::tsparse_scalar::real_norm_value(V_new[j]);
            if (nrm > small_tol)
                for (std::size_t ii = 0; ii < n; ii++) V_new[j][ii] /= T(nrm);
        }

        // Overflow vector: V[m] is already ⊥ span(V[0..m-1]).
        std::vector<T> v_overflow;
        if (m < V.size() && !happy_breakdown) {
            v_overflow = V[m];
            // Re-deflate against locked in case of numerical drift.
            for (int pass = 0; pass < 2; pass++)
                for (std::size_t li = 0; li < locked.size(); li++) {
                    const R lc = vcp::tsparse_scalar::real_dot_value(locked[li].vector, v_overflow);
                    for (std::size_t i = 0; i < n; i++) v_overflow[i] -= T(lc) * locked[li].vector[i];
                }
            const R vovn = vcp::tsparse_scalar::real_norm_value(v_overflow);
            if (vovn > small_tol) for (std::size_t i = 0; i < n; i++) v_overflow[i] /= T(vovn);
        } else {
            // Happy breakdown or missing overflow: generate fresh vector ⊥ V_new and locked.
            v_overflow = vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
            for (int pass = 0; pass < 2; pass++) {
                for (std::size_t j = 0; j < p_keep; j++) {
                    const R c = vcp::tsparse_scalar::real_dot_value(V_new[j], v_overflow);
                    for (std::size_t i = 0; i < n; i++) v_overflow[i] -= T(c) * V_new[j][i];
                }
                for (std::size_t li = 0; li < locked.size(); li++) {
                    const R lc = vcp::tsparse_scalar::real_dot_value(locked[li].vector, v_overflow);
                    for (std::size_t i = 0; i < n; i++) v_overflow[i] -= T(lc) * locked[li].vector[i];
                }
            }
            const R vn = vcp::tsparse_scalar::real_norm_value(v_overflow);
            if (vn > small_tol)
                for (std::size_t i = 0; i < n; i++) v_overflow[i] /= T(vn);
        }

        // Reset H_mat.
        for (std::size_t i = 0; i <= m_limit; i++)
            std::fill(H_mat[i].begin(), H_mat[i].end(), T(0));

        // Set dense restart block.
        for (std::size_t i = 0; i < p_keep; i++)
            for (std::size_t j = 0; j < p_keep; j++)
                H_mat[i][j] = H_proj[i][j];

        // Set coupling row.
        for (std::size_t j = 0; j < p_keep; j++)
            H_mat[p_keep][j] = coupling[j];

        // Reset basis.
        V.clear();
        V.reserve(m_limit + 2);
        for (std::size_t j = 0; j < p_keep; j++) V.push_back(V_new[j]);
        V.push_back(v_overflow);

        m_current = p_keep;
        restart_count++;

    } // end main loop

    // -----------------------------------------------------------------------
    // Collect complex eigenvalues from the final projected Hessenberg.
    // -----------------------------------------------------------------------
    {
        const std::size_t m = m_current;
        if (m > 0) {
            std::vector<std::vector<T> > Hm      = ks_detail::extract_hm(H_mat, m);
            std::vector<std::vector<T> > Hm_hess = Hm;
            std::vector<std::vector<T> > dummy_Q;
            vcp::tsparse_proj_hess::reduce_to_hessenberg_with_q(Hm_hess, dummy_Q);

            const std::size_t hess_iter = m * m * 80 + 300;
            std::vector<C> feigs =
                vcp::tsparse_eigensolvers::hessenberg_complex_eigenvalues<T>(
                    Hm_hess, hess_iter, small_tol);

            const R real_tol = small_tol * R(100);
            for (std::size_t ei = 0; ei < feigs.size(); ei++) {
                if (!ks_detail::is_effectively_real(feigs[ei], real_tol) &&
                    feigs[ei].imag() > R(0))
                {
                    typedef typename vcp::eig_result<T>::eigenvalue_type EV;
                    result.complex_eigenvalues.push_back(
                        EV(feigs[ei].real(),  feigs[ei].imag()));
                    result.complex_eigenvalues.push_back(
                        EV(feigs[ei].real(), -feigs[ei].imag()));
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Final exact residual re-evaluation for all locked pairs.
    // -----------------------------------------------------------------------
    for (std::size_t i = 0; i < locked.size(); i++) {
        std::vector<T> Av;
        apply(locked[i].vector, Av);
        mv_count++;
        const T lam = locked[i].value;
        std::vector<T> r(n);
        for (std::size_t j = 0; j < n; j++) r[j] = Av[j] - lam * locked[i].vector[j];
        locked[i].res_abs = vcp::tsparse_scalar::real_norm_value(r);
        locked[i].res_rel = locked[i].res_abs /
            (R(1) + vcp::tsparse_scalar::abs_value(
                vcp::tsparse_scalar::real_part(lam)));
    }

    // -----------------------------------------------------------------------
    // Sort locked pairs by target and fill result.
    // -----------------------------------------------------------------------
    if (!locked.empty()) {
        std::vector<C> ceigs;
        ceigs.reserve(locked.size());
        for (std::size_t i = 0; i < locked.size(); i++)
            ceigs.push_back(C(vcp::tsparse_scalar::real_part(locked[i].value), R(0)));

        const std::vector<std::size_t> order =
            vcp::tsparse_eigensolvers::select_ritz_indices<T>(
                ceigs, locked.size(), target, shift_val);

        const std::size_t take = std::min(k, order.size());
        for (std::size_t i = 0; i < take; i++) {
            const std::size_t idx = order[i];
            result.eigenvalues.push_back(locked[idx].value);
            result.eigenvectors.push_back(locked[idx].vector);
            result.residuals_absolute.push_back(locked[idx].res_abs);
            result.residuals_relative.push_back(locked[idx].res_rel);
        }
    }

    // -----------------------------------------------------------------------
    // Fill result counts, norms, history, and status.
    // -----------------------------------------------------------------------
    result.returned_real_count    = result.eigenvalues.size();
    result.returned_complex_count = result.complex_eigenvalues.size();
    result.returned_count         =
        result.returned_real_count + result.returned_complex_count;
    result.converged_count        = result.returned_real_count;

    if (!result.residuals_absolute.empty()) {
        result.residual_norm_absolute =
            *std::max_element(result.residuals_absolute.begin(),
                              result.residuals_absolute.end());
    }
    if (!result.residuals_relative.empty()) {
        result.residual_norm_relative =
            *std::max_element(result.residuals_relative.begin(),
                              result.residuals_relative.end());
    }

    if (!compute_hist) {
        result.residual_history_absolute.clear();
        result.residual_history_relative.clear();
    }

    result.matrix_vector_products = mv_count;

    if (result.returned_real_count >= k) {
        result.converged = true;
        result.status    = "converged";
        result.message   = "krylov_schur_experimental converged";
    } else if (budget_exhausted || mv_count >= max_mv) {
        result.converged      = false;
        result.status         = "max_iter_exhausted";
        result.failure_reason =
            "matrix-vector product budget exhausted before full convergence";
        result.message        = "krylov_schur: budget exhausted";
    } else {
        result.converged = false;
        result.status    = "failed";
        if (result.failure_reason.empty()) {
            if (!result.complex_eigenvalues.empty() &&
                result.returned_real_count < k) {
                result.failure_reason =
                    "insufficient real eigenvalues converged; complex eigenvalues detected "
                    "(complex restart is outside Phase 3 scope)";
            } else {
                result.failure_reason =
                    "eigensolver terminated without full convergence";
            }
        }
        result.message = "krylov_schur: not all eigenvalues converged";
    }

    diag.restart_count          = restart_count;
    diag.matrix_vector_products = mv_count;
    diag.locked_real_count      = result.returned_real_count;
    diag.locked_complex_count   = result.returned_complex_count / 2;

    return diag;
}

// ===========================================================================
// Simplified wrapper: returns only eig_result<T>.
// ===========================================================================
template <class Apply, class T>
vcp::eig_result<T> krylov_schur_eigs(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    return krylov_schur_eigs_with_diagnostics<Apply, T>(
        apply, n, k, options).eigs;
}

} // namespace tsparse_experimental
} // namespace vcp

#endif // VCP_TSPARSE_KRYLOV_SCHUR_HPP
