// VCP Library
// http://verified.computation.jp
//
// vcp/tsparse/tsparse_thick_restart_lanczos.hpp
//
// Experimental thick-restart Lanczos eigensolver for REAL SYMMETRIC operators.
// Implements Wu & Simon (2000) thick-restart with arrowhead projected matrix
// and progressive locking for correct repeated-eigenvalue discovery.
//
// Namespace: vcp::tsparse_experimental
//
// PUBLIC API IS NOT CONNECTED to vcp::spmatrix::eigs / eig_method enum.
// Complex types are rejected via static_assert.
//
// Apply functor convention:
//   void operator()(const std::vector<T>& x, std::vector<T>& y) const;
//
// options.max_iter = total matrix-vector product budget (main loop only).
//
// PRECONDITION (EIG-1 F-6, D-12): Apply must represent a REAL SYMMETRIC
// operator.  The Apply abstraction cannot verify symmetry at runtime; a
// violation is made honest by the end-of-run exact-residual acceptance check
// (EIG-0 C-1): the returned pairs then carry large true residuals and the
// result is demoted to converged=false / status="residual_check_failed"
// instead of silent garbage.
//
// Termination contract (EIG-1 F-3; EIG-0 C-1/C-2):
//   converged=true requires (i) k pairs returned, (ii) no unconverged Ritz
//   candidate certainly inside the returned set at termination
//   (honest_termination_check_), and (iii) all end-of-run EXACT residuals
//   pass the acceptance test (residual_acceptance_check_).
//   Note (EIG-0 C-4): converged=true is a residual claim plus absence of a
//   visible contradiction; it is NOT a completeness guarantee (Krylov methods
//   cannot see eigenspaces orthogonal to the start vector).

#pragma once

#ifndef VCP_TSPARSE_THICK_RESTART_LANCZOS_HPP
#define VCP_TSPARSE_THICK_RESTART_LANCZOS_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <vcp/tsparse/tsparse_restart.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>
#include <vcp/tsparse/tsparse_projected_eigensolver.hpp>
#include <vcp/spmatrix.hpp>

namespace vcp {
namespace tsparse_experimental {

// ---------------------------------------------------------------------------
// Diagnostic result wrapper
// ---------------------------------------------------------------------------
template <class T>
struct thick_restart_lanczos_result {
    vcp::eig_result<T> eigs;
    std::size_t restart_count;           // thick restarts performed
    std::size_t matrix_vector_products;  // total apply() calls
    std::size_t locked_count;            // converged + locked pairs
};

// ---------------------------------------------------------------------------
// Type guard: only real floating-point T is supported.
// ---------------------------------------------------------------------------
template <class T>
struct trl_is_real_floating {
    static const bool value =
        std::is_floating_point<T>::value &&
        !vcp::tsparse_scalar::is_complex<T>::value;
};

// ===========================================================================
// Internal helpers
// ===========================================================================
namespace trl_detail {

// Effective subspace dimension m, always satisfying k <= m <= n.
// Fix: m = min(k+1, n) when m < k+1, preventing m > n when k == n.
inline std::size_t effective_m(std::size_t n, std::size_t k, std::size_t requested)
{
    if (n == 0 || k == 0) return 0;
    std::size_t m;
    if (requested == 0) {
        m = std::max<std::size_t>(
                std::max<std::size_t>(2 * k + 6, std::size_t(20)), k + 3);
    } else if (requested <= k) {
        m = k + 2;
    } else {
        m = requested;
    }
    if (m > n) m = n;
    if (m < k + 1) m = std::min(k + 1, n);
    return m;
}

// Orthogonalize v against locked (weighted) and optional extra unit-norm vectors.
// Returns ||v|| after orthogonalization.
template <class T>
typename vcp::tsparse_scalar::real_type<T>::type
orth_and_norm(
    std::vector<T>& v,
    const std::vector<vcp::tsparse::locked_pair<T> >& locked,
    const std::vector<std::vector<T> >* extra_vecs)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = v.size();

    for (int pass = 0; pass < 2; pass++) {
        for (std::size_t li = 0; li < locked.size(); li++) {
            const std::vector<T>& q = locked[li].vector;
            if (q.size() != n) continue;
            R q2 = R(0);
            T dot = T(0);
            for (std::size_t i = 0; i < n; i++) {
                q2 += vcp::tsparse_scalar::real_part(
                    vcp::tsparse_scalar::conjugate_if_needed(q[i]) * q[i]);
                dot += vcp::tsparse_scalar::conjugate_if_needed(q[i]) * v[i];
            }
            if (q2 > std::numeric_limits<R>::epsilon()) {
                const T c = dot / T(q2);
                for (std::size_t i = 0; i < n; i++) v[i] -= c * q[i];
            }
        }
        if (extra_vecs) {
            for (std::size_t j = 0; j < extra_vecs->size(); j++) {
                const R c = vcp::tsparse_scalar::real_dot_value((*extra_vecs)[j], v);
                for (std::size_t i = 0; i < n; i++) v[i] -= T(c) * (*extra_vecs)[j][i];
            }
        }
    }
    return vcp::tsparse_scalar::real_norm_value(v);
}

// Generate a start vector (tries multiple seeds until orthogonal to locked).
template <class T>
std::vector<T> new_start_vector(
    std::size_t n,
    unsigned int& seed,
    const std::vector<vcp::tsparse::locked_pair<T> >& locked,
    const std::vector<std::vector<T> >& extra_vecs)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const R eps_n = std::numeric_limits<R>::epsilon() * R(n + 1);

    for (int attempt = 0; attempt < 60; attempt++) {
        std::vector<T> v =
            vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
        const R nrm = orth_and_norm(v, locked, &extra_vecs);
        if (nrm > eps_n) {
            for (std::size_t i = 0; i < n; i++) v[i] /= T(nrm);
            return v;
        }
    }
    std::vector<T> fb(n, T(0));
    if (n > 0) fb[0] = T(1);
    return fb;
}

// Full reorthogonalization of z against basis (unit-norm) and locked.
// Returns ||z|| after orthogonalization.
template <class T>
typename vcp::tsparse_scalar::real_type<T>::type
reorthogonalize(
    std::vector<T>& z,
    const std::vector<std::vector<T> >& basis,
    const std::vector<vcp::tsparse::locked_pair<T> >& locked)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = z.size();

    for (int pass = 0; pass < 2; pass++) {
        for (std::size_t li = 0; li < locked.size(); li++) {
            const std::vector<T>& q = locked[li].vector;
            if (q.size() != n) continue;
            R q2 = R(0);
            T dot = T(0);
            for (std::size_t i = 0; i < n; i++) {
                q2 += vcp::tsparse_scalar::real_part(
                    vcp::tsparse_scalar::conjugate_if_needed(q[i]) * q[i]);
                dot += vcp::tsparse_scalar::conjugate_if_needed(q[i]) * z[i];
            }
            if (q2 > std::numeric_limits<R>::epsilon()) {
                const T c = dot / T(q2);
                for (std::size_t i = 0; i < n; i++) z[i] -= c * q[i];
            }
        }
        for (std::size_t j = 0; j < basis.size(); j++) {
            const R c = vcp::tsparse_scalar::real_dot_value(basis[j], z);
            for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * basis[j][i];
        }
    }
    return vcp::tsparse_scalar::real_norm_value(z);
}

// Build projected symmetric matrix with ARROWHEAD structure.
//
// Standard tridiagonal (k_restart == 0):
//   M[i][i] = alpha[i],  M[i][i+1] = M[i+1][i] = beta[i]
//
// Arrowhead (k_restart > 0):
//   M[i][i]             = alpha[i]
//   M[j][k_restart]     = M[k_restart][j] = beta[j]  (j < k_restart)
//   M[i][i+1]           = M[i+1][i]       = beta[i]  (i >= k_restart)
template <class T>
std::vector<std::vector<T> > build_projected_matrix(
    const std::vector<T>& alpha,
    const std::vector<T>& beta,
    std::size_t k_restart)
{
    const std::size_t m = alpha.size();
    std::vector<std::vector<T> > M(m, std::vector<T>(m, T(0)));
    for (std::size_t i = 0; i < m; i++) M[i][i] = alpha[i];

    if (k_restart == 0 || k_restart >= m) {
        for (std::size_t i = 0; i + 1 < m && i < beta.size(); i++) {
            M[i][i + 1] = M[i + 1][i] = beta[i];
        }
    } else {
        for (std::size_t j = 0; j < k_restart && j < beta.size(); j++) {
            M[j][k_restart] = M[k_restart][j] = beta[j];
        }
        for (std::size_t i = k_restart; i + 1 < m && i < beta.size(); i++) {
            M[i][i + 1] = M[i + 1][i] = beta[i];
        }
    }
    return M;
}

// EIG-1 F-3-1 freshness guard: indices (ascending) of the locked pairs that
// would form the returned target-order prefix k right now.  Used to decide
// whether the current subspace postdates the last change of the returned set
// (see the termination check in the main loop).
template <class T>
std::vector<std::size_t> locked_prefix_indices(
    const std::vector<vcp::tsparse::locked_pair<T> >& locked,
    std::size_t k,
    eig_target target,
    const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
    std::vector<T> vals;
    vals.reserve(locked.size());
    for (std::size_t i = 0; i < locked.size(); i++) vals.push_back(locked[i].value);
    // コア実装は tsparse_honest_termination.hpp(1 箇所)。
    return vcp::tsparse::locked_prefix_indices_(vals, k, target, shift);
}

// Lift projected eigenvector to full space: u = sum_j basis[j] * y[j]
template <class T>
std::vector<T> lift_ritz(
    const std::vector<std::vector<T> >& basis,
    const std::vector<T>& y,
    std::size_t n)
{
    std::vector<T> u(n, T(0));
    const std::size_t len = std::min(y.size(), basis.size());
    for (std::size_t j = 0; j < len; j++) {
        for (std::size_t i = 0; i < n; i++) u[i] += basis[j][i] * y[j];
    }
    return u;
}

} // namespace trl_detail

// ===========================================================================
// thick_restart_lanczos_eigs_with_diagnostics
// ===========================================================================
template <class Apply, class T>
thick_restart_lanczos_result<T> thick_restart_lanczos_eigs_with_diagnostics(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    static_assert(
        trl_is_real_floating<T>::value,
        "thick_restart_lanczos_eigs: T must be real floating-point. "
        "Complex types are not supported.");

    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    thick_restart_lanczos_result<T> diag;
    vcp::eig_result<T>& result = diag.eigs;
    diag.restart_count          = 0;
    diag.matrix_vector_products = 0;
    diag.locked_count           = 0;

    result.used_dense_fallback       = false;
    result.used_shift_invert         = false;
    result.used_generalized_operator = false;
    result.used_method               = "thick_restart_lanczos_experimental";
    result.used_orthogonalization    = "full_reorthogonalization";

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

    // Clamp k to n; preserve original k in requested_count.
    const std::size_t k_original = k;
    if (k > n) k = n;
    result.requested_count = k_original;

    // -----------------------------------------------------------------------
    // Options
    // -----------------------------------------------------------------------
    const R tol = (options.tol > R(0))
        ? options.tol
        : vcp::tsparse_scalar::decimal_power_negative<R>(12);

    const std::size_t max_mv = (options.max_iter > 0)
        ? options.max_iter
        : std::size_t(200) * n;

    const std::size_t m_limit =
        trl_detail::effective_m(n, k, options.subspace_dim);
    result.used_subspace_dim = m_limit;

    const eig_target target       = options.target;
    const R          shift        = options.shift;
    const bool       compute_hist = options.compute_residual_history;
    unsigned int     seed         =
        options.random_start ? options.random_seed : 0u;

    const R small_tol =
        (tol * R(1e-4) > std::numeric_limits<R>::epsilon() * R(n + 1))
        ? tol * R(1e-4)
        : std::numeric_limits<R>::epsilon() * R(n + 1);

    // -----------------------------------------------------------------------
    // State
    // -----------------------------------------------------------------------
    std::vector<vcp::tsparse::locked_pair<T> > locked;
    locked.reserve(k + 1);

    std::vector<std::vector<T> > basis;
    basis.reserve(m_limit + 2);

    // alpha[i] = projected diagonal
    // beta[i]  = coupling (arrowhead) or off-diagonal (tridiagonal)
    std::vector<T> alpha;
    alpha.reserve(m_limit + 1);
    std::vector<T> beta;
    beta.reserve(m_limit + 1);

    // k_restart = number of retained Ritz vectors at start of current cycle
    std::size_t k_restart = 0;

    // EIG-4 T-4: running projected-operator infinity-norm estimate (shared
    // C-1 scale input; updated once per restart, zero extra applies)
    R anorm_est_run = R(0);
    // beta_overflow = ||z|| from last expansion step (Wu & Simon f_m)
    // z_overflow   = z / ||z|| (unit-norm overflow Lanczos vector)
    R               beta_overflow = R(0);
    std::vector<T>  z_overflow(n, T(0));

    std::size_t mv_count       = 0;
    std::size_t total_restarts = 0;

    // EIG-1 F-3-1 (EIG-0 C-2): true while termination is held back because an
    // unconverged Ritz candidate certainly lies inside the would-be returned
    // prefix k.  Re-evaluated every cycle once locked.size() >= k.
    bool honest_hold = false;

    // EIG-1 F-3-1 freshness guard: the returned-prefix indices as of the last
    // basis (re)build.  C-2 may only clear termination on a subspace that was
    // built orthogonal to a locked set whose returned prefix equals the
    // current one; otherwise a multi-lock cycle could terminate on evidence
    // from a spent subspace in which a not-yet-surfaced inner candidate (e.g.
    // a further multiplicity copy) is invisible (t1 mechanism, D-1).
    std::vector<std::size_t> locked_prefix_at_restart;

    // -----------------------------------------------------------------------
    // Initial starting vector
    // -----------------------------------------------------------------------
    {
        std::vector<std::vector<T> > empty;
        basis.push_back(
            trl_detail::new_start_vector<T>(n, seed, locked, empty));
    }

    // -----------------------------------------------------------------------
    // Main loop
    // -----------------------------------------------------------------------
    // EIG-1 F-3-1: keep iterating while the C-2 check holds termination back,
    // even though locked.size() >= k (the lock set is kept, not discarded).
    while ((locked.size() < k || honest_hold) && mv_count < max_mv) {

        // -------------------------------------------------------------------
        // PHASE 1: Expand Lanczos basis to m_limit steps.
        //
        // j == k_restart: multi-term recurrence subtracts all retained betas.
        // j >  k_restart: standard 3-term recurrence.
        // Break on happy breakdown (nrm <= small_tol) or full subspace.
        // -------------------------------------------------------------------
        beta_overflow = R(0);
        std::fill(z_overflow.begin(), z_overflow.end(), T(0));

        while (alpha.size() < m_limit && mv_count < max_mv) {
            const std::size_t j = alpha.size();
            const std::vector<T>& vj = basis[j];

            std::vector<T> z;
            apply(vj, z);
            mv_count++;

            if (j > 0) {
                if (j == k_restart && k_restart > 0) {
                    for (std::size_t j2 = 0; j2 < k_restart && j2 < beta.size(); j2++) {
                        const T cj2 = beta[j2];
                        for (std::size_t i = 0; i < n; i++)
                            z[i] -= cj2 * basis[j2][i];
                    }
                } else {
                    if (j - 1 < beta.size()) {
                        const T bprev = beta[j - 1];
                        const std::vector<T>& vprev = basis[j - 1];
                        for (std::size_t i = 0; i < n; i++)
                            z[i] -= bprev * vprev[i];
                    }
                }
            }

            const R alpha_j_r = vcp::tsparse_scalar::real_dot_value(vj, z);
            const T alpha_j   = T(alpha_j_r);
            for (std::size_t i = 0; i < n; i++) z[i] -= alpha_j * vj[i];

            const R nrm = trl_detail::reorthogonalize(z, basis, locked);

            alpha.push_back(alpha_j);

            if (j + 1 >= m_limit || nrm <= small_tol) {
                beta_overflow = nrm;
                if (nrm > small_tol) {
                    for (std::size_t i = 0; i < n; i++) z_overflow[i] = z[i] / T(nrm);
                } else {
                    std::fill(z_overflow.begin(), z_overflow.end(), T(0));
                }
                break;
            }

            beta.push_back(T(nrm));
            std::vector<T> vnew(n);
            for (std::size_t i = 0; i < n; i++) vnew[i] = z[i] / T(nrm);
            basis.push_back(vnew);
        }

        if (alpha.empty()) {
            result.failure_reason = "Lanczos produced empty basis";
            break;
        }

        // -------------------------------------------------------------------
        // PHASE 2: Projected eigenproblem (arrowhead or tridiagonal).
        // -------------------------------------------------------------------
        const std::size_t m_actual = alpha.size();

        while (beta.size() + 1 < m_actual) beta.push_back(T(0));

        const std::vector<std::vector<T> > T_proj =
            trl_detail::build_projected_matrix(alpha, beta, k_restart);

        // EIG-4 T-4: running ||A|| estimate = max over restarts of the
        // projected operator's infinity norm (definition documented at the
        // shared helper, tsparse_honest_termination.hpp).  Zero extra applies;
        // the run trajectory is untouched (final-acceptance input only).
        for (std::size_t i2 = 0; i2 < T_proj.size(); i2++) {
            R rowsum(0);
            for (std::size_t j2 = 0; j2 < T_proj[i2].size(); j2++)
                rowsum += vcp::tsparse_scalar::abs_value(T_proj[i2][j2]);
            if (rowsum > anorm_est_run) anorm_est_run = rowsum;
        }

        const std::size_t proj_iter =
            std::max<std::size_t>(m_actual * m_actual * 300, std::size_t(3000));

        vcp::tsparse_projected::projected_eigensolver_result<T> proj_sym =
            vcp::tsparse_projected::solve_real_symmetric_projected(
                T_proj, proj_iter, small_tol);

        if (!proj_sym.success) {
            result.failure_reason = proj_sym.message.empty()
                ? "projected eigensolver failed" : proj_sym.message;
            break;
        }
        if (proj_sym.pairs.empty()) {
            result.failure_reason = "projected eigensolver returned empty result";
            break;
        }

        // -------------------------------------------------------------------
        // PHASE 3: Select Ritz pairs and compute Wu & Simon residual estimates.
        //
        //   ||A*u_j - theta_j*u_j|| ≈ |beta_overflow * y_j[m-1]|
        //
        // Ritz vectors are lifted to full space and orthogonalized against
        // locked vectors using the Phase 1 helper orthogonalize_against_locked.
        // -------------------------------------------------------------------
        const std::size_t n_ritz = proj_sym.pairs.size();
        const std::size_t k_extra =
            std::max<std::size_t>(k / 2, std::size_t(2));
        const std::size_t k_want = std::min(n_ritz, k + k_extra);

        const std::vector<std::size_t> sel_order =
            vcp::tsparse_projected::select_projected_indices(
                proj_sym.pairs, k_want, target, shift);

        struct ritz_data_t {
            T              value;
            std::vector<T> vector;    // lifted, unit-norm, orth vs locked
            R              res_abs;   // |beta_overflow * y_j[m-1]|
            R              res_rel;   // res_abs / (1 + |theta|)
            bool           converged;
            T              proj_last; // y_j[m_actual-1] (Wu & Simon coupling)
        };
        std::vector<ritz_data_t> ritz_sel;
        ritz_sel.reserve(sel_order.size());

        for (std::size_t si = 0; si < sel_order.size(); si++) {
            const std::size_t idx = sel_order[si];
            if (idx >= proj_sym.pairs.size()) continue;

            const T proj_last = proj_sym.pairs[idx].vector.empty()
                ? T(0) : proj_sym.pairs[idx].vector.back();
            const R res_abs_val =
                beta_overflow * vcp::tsparse_scalar::abs_value(
                    vcp::tsparse_scalar::real_part(proj_last));

            const T theta = proj_sym.pairs[idx].value;
            const R abs_theta = vcp::tsparse_scalar::abs_value(
                vcp::tsparse_scalar::real_part(theta));
            const R res_rel_val = res_abs_val / (R(1) + abs_theta);

            // Lift and orthogonalize against locked using Phase 1 helper
            std::vector<T> u = trl_detail::lift_ritz(basis, proj_sym.pairs[idx].vector, n);
            const R unrm = vcp::tsparse::orthogonalize_against_locked(u, locked, false);
            if (unrm <= small_tol) continue;
            for (std::size_t i = 0; i < n; i++) u[i] /= T(unrm);

            ritz_data_t rd;
            rd.value     = theta;
            rd.vector    = u;
            rd.res_abs   = res_abs_val;
            rd.res_rel   = res_rel_val;
            rd.converged = (res_rel_val <= tol || res_abs_val <= tol);
            rd.proj_last = proj_last;
            ritz_sel.push_back(rd);
        }

        // -------------------------------------------------------------------
        // PHASE 4: Record residual history.
        // -------------------------------------------------------------------
        if (compute_hist && !ritz_sel.empty()) {
            R best_abs = ritz_sel[0].res_abs;
            R best_rel = ritz_sel[0].res_rel;
            for (std::size_t i = 1; i < ritz_sel.size(); i++) {
                if (ritz_sel[i].res_abs < best_abs) best_abs = ritz_sel[i].res_abs;
                if (ritz_sel[i].res_rel < best_rel) best_rel = ritz_sel[i].res_rel;
            }
            result.residual_history_absolute.push_back(best_abs);
            result.residual_history_relative.push_back(best_rel);
        }

        // -------------------------------------------------------------------
        // PHASE 5: Lock converged pairs via lock_converged_pairs.
        //
        // EIG-1 F-1 (D-9 root fix): every candidate that converged in this
        // cycle is locked, in target order (active_rp preserves the
        // target-sorted sel_order).  Locking may exceed k; what is returned
        // is decided by the final target-ordered prefix-k selection
        // (take_locked_prefix).  Repeated eigenvalues are still discovered
        // progressively because each new cycle orthogonalizes against all
        // locked vectors.
        // -------------------------------------------------------------------
        // EIG-1 F-3-1: set when this cycle changed the returned prefix while
        // locked.size() >= k; the following restart must then start from a
        // fresh random direction (see PHASE 6).
        bool fresh_restart_needed = false;

        std::vector<vcp::tsparse::ritz_pair<T> > active_rp;
        std::size_t n_locked_now = 0;
        {
            active_rp.reserve(ritz_sel.size());
            for (std::size_t i = 0; i < ritz_sel.size(); i++) {
                vcp::tsparse::ritz_pair<T> rp;
                rp.value             = ritz_sel[i].value;
                rp.vector            = ritz_sel[i].vector;
                rp.residual_absolute = ritz_sel[i].res_abs;
                rp.residual_relative = ritz_sel[i].res_rel;
                rp.converged         = ritz_sel[i].converged;
                active_rp.push_back(rp);
            }
            n_locked_now = vcp::tsparse::lock_converged_pairs(
                active_rp, locked, active_rp.size());
        }
        if (locked.size() >= k) {
            // ---------------------------------------------------------------
            // EIG-1 F-3-1 (EIG-0 C-2, honest stopping rule): do not terminate
            // while an unconverged active candidate certainly lies inside
            // (more target-preferred than) the worst value of the would-be
            // returned prefix k.  The lock set is kept and iteration
            // continues; if the budget runs out first, the result is an
            // honest not_converged (max_iter_exhausted).
            //
            // Freshness guard: if this cycle changed the returned prefix,
            // the current candidates were computed in a subspace that did not
            // yet know the final locked set, so an inner candidate (e.g. a
            // further multiplicity copy) may simply not have surfaced yet.
            // Termination then waits for at least one verification cycle
            // whose subspace was built orthogonal to the current prefix.
            // Locks OUTSIDE the prefix do not retrigger verification, so
            // outward-converging cycles cannot stall termination forever.
            // ---------------------------------------------------------------
            const std::vector<std::size_t> cur_prefix =
                trl_detail::locked_prefix_indices(locked, k, target, shift);
            const bool prefix_fresh = (cur_prefix == locked_prefix_at_restart);
            const bool no_inner_unconverged =
                vcp::tsparse::honest_termination_check_(
                    active_rp, locked, k, target, shift);
            honest_hold = !(no_inner_unconverged && prefix_fresh);
            if (!honest_hold) break;
            // Blind-spot probe: when the value check sees no inner
            // unconverged candidate BUT the subspace predates the current
            // prefix, the evidence may be blind — a single Krylov sequence
            // sees at most one direction per eigenspace, so a further
            // multiplicity copy of a locked value can stay invisible forever
            // (t1 mechanism).  The next restart then starts from a fresh
            // random direction orthogonal to the current locked set (the
            // progressive-locking multiplicity discovery, promoted to a
            // termination requirement).
            // If an inner unconverged candidate IS visible, keep the normal
            // thick restart instead so its convergence progress is retained
            // (lap2d second-copy mechanism) — probing here would destroy it.
            if (no_inner_unconverged && !prefix_fresh) fresh_restart_needed = true;
        }

        // -------------------------------------------------------------------
        // PHASE 6: Thick restart (Wu & Simon 2000).
        //
        // EIG-1 F-2 (invariant, B-15): a Ritz pair judged converged is never
        // silently discarded.  Retention takes everything NOT locked this
        // cycle, regardless of its converged flag.  With F-1's uncapped
        // locking the converged-but-unlocked set is empty by construction;
        // this code pins the invariant against future locking policy changes.
        // Coupling:  c_j = beta_overflow * y_j[m_actual-1]
        // v_new = z_overflow, orthogonalized against locked via Phase 1 helper.
        // -------------------------------------------------------------------
        const std::size_t max_retained =
            (m_limit >= std::size_t(2)) ? (m_limit - 2) : std::size_t(0);

        std::vector<ritz_data_t> retained;
        retained.reserve(std::min(k_want, max_retained));
        {
            std::size_t conv_seen = 0;
            for (std::size_t i = 0;
                 i < ritz_sel.size() && retained.size() < max_retained; i++) {
                bool locked_this_cycle = false;
                if (ritz_sel[i].converged) {
                    // lock_converged_pairs took the first n_locked_now
                    // converged entries in active (= sel) order.
                    if (conv_seen < n_locked_now) locked_this_cycle = true;
                    conv_seen++;
                }
                if (locked_this_cycle) continue;
                // Verification cycle: drop UNCONVERGED retained pairs so the
                // restart is a fresh random probe orthogonal to locked
                // (multiplicity discovery).  B-15 only protects converged
                // pairs; those are still retained even here.
                if (fresh_restart_needed && !ritz_sel[i].converged) continue;
                retained.push_back(ritz_sel[i]);
            }
        }

        const std::size_t k_ret = retained.size();

        std::vector<std::vector<T> > new_basis;
        std::vector<T>               new_alpha;
        std::vector<T>               new_beta;
        new_basis.reserve(k_ret + 2);
        new_alpha.reserve(k_ret + 1);
        new_beta.reserve(k_ret + 1);

        for (std::size_t i = 0; i < k_ret; i++) {
            new_basis.push_back(retained[i].vector);
            new_alpha.push_back(retained[i].value);
        }

        if (k_ret > 0 && beta_overflow > small_tol && !fresh_restart_needed) {
            std::vector<T> v_new = z_overflow;

            // Orthogonalize v_new against locked using Phase 1 helper
            vcp::tsparse::orthogonalize_against_locked(v_new, locked, false);

            const R vnew_nrm = vcp::tsparse_scalar::real_norm_value(v_new);
            if (vnew_nrm > small_tol) {
                for (std::size_t i = 0; i < n; i++) v_new[i] /= T(vnew_nrm);
                new_basis.push_back(v_new);

                // Wu & Simon coupling: c_j = beta_overflow * y_j[m-1]
                for (std::size_t j = 0; j < k_ret; j++) {
                    new_beta.push_back(T(beta_overflow) * retained[j].proj_last);
                }
            } else {
                // v_new collapsed; use fresh random direction orthogonal to all
                std::vector<std::vector<T> > extra;
                extra.reserve(k_ret);
                for (std::size_t i = 0; i < k_ret; i++)
                    extra.push_back(retained[i].vector);
                new_basis.push_back(
                    trl_detail::new_start_vector<T>(n, seed, locked, extra));
                for (std::size_t j = 0; j < k_ret; j++) new_beta.push_back(T(0));
            }
        } else if (k_ret > 0) {
            // Happy breakdown with no retained unconverged pairs will fall
            // into the else branch below; this branch handles retained > 0
            // but beta_overflow ≈ 0 (invariant subspace + retained pairs).
            std::vector<std::vector<T> > extra;
            extra.reserve(k_ret);
            for (std::size_t i = 0; i < k_ret; i++)
                extra.push_back(retained[i].vector);
            new_basis.push_back(
                trl_detail::new_start_vector<T>(n, seed, locked, extra));
            for (std::size_t j = 0; j < k_ret; j++) new_beta.push_back(T(0));
        } else {
            // No retained pairs: fresh random restart orthogonal to all locked.
            // This is the key step for repeated-eigenvalue discovery:
            // the new start vector is orthogonal to all previously locked
            // eigenvectors, so the next Krylov space can find additional
            // eigenvectors of the same eigenvalue.
            std::vector<std::vector<T> > empty;
            new_basis.push_back(
                trl_detail::new_start_vector<T>(n, seed, locked, empty));
        }

        basis.clear();
        alpha.clear();
        beta.clear();
        for (std::size_t i = 0; i < new_basis.size(); i++)
            basis.push_back(new_basis[i]);
        basis.reserve(m_limit + 2);
        alpha     = new_alpha;
        beta      = new_beta;
        k_restart = k_ret;

        // EIG-1 F-3-1 freshness guard: the basis just built is (and will be
        // kept, via reorthogonalize) orthogonal to the current locked set;
        // record which returned prefix that subspace knows about.
        locked_prefix_at_restart =
            trl_detail::locked_prefix_indices(locked, k, target, shift);

        result.iterations++;
        total_restarts++;

    } // end main while

    // -----------------------------------------------------------------------
    // Budget flag (captured before re-evaluation calls)
    // -----------------------------------------------------------------------
    const bool budget_exhausted =
        (mv_count >= max_mv) && (locked.size() < k || honest_hold);

    // -----------------------------------------------------------------------
    // Re-evaluate actual residuals via apply(): ||Av - lambda*v||
    // Replaces Wu & Simon estimates with exact operator residuals.
    // -----------------------------------------------------------------------
    for (std::size_t i = 0; i < locked.size(); i++) {
        std::vector<T> w;
        apply(locked[i].vector, w);
        mv_count++;
        const T lam = locked[i].value;
        std::vector<T> diff(n);
        for (std::size_t j = 0; j < n; j++) diff[j] = w[j] - lam * locked[i].vector[j];
        const R res_abs = vcp::tsparse_scalar::real_norm_value(diff);
        const R abs_lam = vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(lam));
        const R res_rel = res_abs / (R(1) + abs_lam);
        locked[i].residual_absolute = res_abs;
        locked[i].residual_relative = res_rel;
    }

    // -----------------------------------------------------------------------
    // Clear or keep residual history
    // -----------------------------------------------------------------------
    if (!compute_hist) {
        result.residual_history_absolute.clear();
        result.residual_history_relative.clear();
    }

    diag.restart_count            = total_restarts;
    diag.matrix_vector_products   = mv_count;
    diag.locked_count             = locked.size();
    result.matrix_vector_products = mv_count;

    // -----------------------------------------------------------------------
    // Sort locked pairs by target, take first k via take_locked_prefix
    // -----------------------------------------------------------------------
    if (!locked.empty()) {
        std::vector<std::complex<R> > ceigs_locked;
        ceigs_locked.reserve(locked.size());
        for (std::size_t i = 0; i < locked.size(); i++) {
            ceigs_locked.push_back(std::complex<R>(
                vcp::tsparse_scalar::real_part(locked[i].value), R(0)));
        }
        const std::vector<std::size_t> order =
            vcp::tsparse_eigensolvers::select_ritz_indices<T>(
                ceigs_locked, locked.size(), target, shift);

        std::vector<vcp::tsparse::locked_pair<T> > sorted_locked;
        sorted_locked.reserve(locked.size());
        for (std::size_t i = 0; i < order.size(); i++) {
            sorted_locked.push_back(locked[order[i]]);
        }

        const std::vector<vcp::tsparse::locked_pair<T> > final_locked =
            vcp::tsparse::take_locked_prefix(sorted_locked, k);

        for (std::size_t i = 0; i < final_locked.size(); i++) {
            result.eigenvalues.push_back(final_locked[i].value);
            result.eigenvectors.push_back(final_locked[i].vector);
            result.residuals_absolute.push_back(final_locked[i].residual_absolute);
            result.residuals_relative.push_back(final_locked[i].residual_relative);
        }
    }

    result.returned_count         = result.eigenvalues.size();
    result.returned_real_count    = result.returned_count;
    result.returned_complex_count = 0;
    result.converged_count        = result.returned_count;

    // Convergence is based on the clamped k (not k_original).
    //
    // EIG-1 F-3 (EIG-0 C-1/C-2): converged=true requires
    //   (i)  k pairs returned,
    //   (ii) no honest-termination hold at exit (C-2), and
    //   (iii) the end-of-run EXACT residuals of all returned pairs pass the
    //        acceptance test.  EIG-4 T-4 (D4-4 / R-1): the acceptance is the
    //        legacy test (res_abs <= tol or res_rel <= tol, scale = 1+|theta|)
    //        OR res_abs <= tol * max(1+|theta|, anorm_est_run) via the shared
    //        helper (strictly a relaxation -- B-29; anorm_est_run = running
    //        projected-operator infinity norm, see the helper's docs).  This
    //        still demotes non-symmetric misuse (D-12) to an honest failure.
    std::vector<R> theta_abs_c1;
    theta_abs_c1.reserve(result.eigenvalues.size());
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++)
        theta_abs_c1.push_back(vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(result.eigenvalues[i])));
    if (result.returned_count >= k && !honest_hold &&
        vcp::tsparse::residual_acceptance_check_scaled_(
            result.residuals_absolute, result.residuals_relative, tol,
            theta_abs_c1, anorm_est_run)) {
        result.converged = true;
        result.status    = "converged";
        result.message   = "thick_restart_lanczos converged";
    } else if (budget_exhausted) {
        result.converged = false;
        result.status    = "max_iter_exhausted";
        if (result.failure_reason.empty()) {
            result.failure_reason =
                "matrix-vector product budget exhausted before full convergence";
        }
        result.message = "thick_restart_lanczos: budget exhausted";
    } else if (result.returned_count >= k) {
        // k pairs are present but the converged claim was refused by C-1/C-2.
        // Values and residuals are still returned as diagnostics (B-15: the
        // locked pairs are not discarded).
        result.converged = false;
        result.status    = "residual_check_failed";
        if (result.failure_reason.empty()) {
            result.failure_reason = honest_hold
                ? "unconverged candidate inside returned set at termination (C-2)"
                : "end-of-run exact residual failed acceptance (C-1)";
        }
        result.message =
            "thick_restart_lanczos: converged claim rejected by exact residual "
            "/ honest termination check";
    } else {
        result.converged = false;
        result.status    = "failed";
        if (result.failure_reason.empty()) {
            result.failure_reason = "eigensolver terminated without convergence";
        }
        result.message = "thick_restart_lanczos: not all eigenvalues converged";
    }

    if (!result.residuals_absolute.empty()) {
        R max_abs = result.residuals_absolute[0];
        R max_rel = result.residuals_relative.empty() ? R(0)
                                                       : result.residuals_relative[0];
        for (std::size_t i = 1; i < result.residuals_absolute.size(); i++) {
            if (result.residuals_absolute[i] > max_abs)
                max_abs = result.residuals_absolute[i];
        }
        for (std::size_t i = 1; i < result.residuals_relative.size(); i++) {
            if (result.residuals_relative[i] > max_rel)
                max_rel = result.residuals_relative[i];
        }
        result.residual_norm_absolute = max_abs;
        result.residual_norm_relative = max_rel;
    }

    return diag;
}

// ===========================================================================
// Simplified wrapper: returns only eig_result<T>.
// ===========================================================================
template <class Apply, class T>
vcp::eig_result<T> thick_restart_lanczos_eigs(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    return thick_restart_lanczos_eigs_with_diagnostics<Apply, T>(
        apply, n, k, options).eigs;
}

} // namespace tsparse_experimental
} // namespace vcp

#endif // VCP_TSPARSE_THICK_RESTART_LANCZOS_HPP
