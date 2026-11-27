// VCP Library
// http://verified.computation.jp
//
// vcp/tsparse/tsparse_thick_restart_lanczos.hpp
//
// Thick-restart Lanczos eigensolver for REAL SYMMETRIC operators and, since
// EIG-T2, COMPLEX HERMITIAN operators (T = std::complex of a real type).
// Implements Wu & Simon (2000) thick-restart with arrowhead projected matrix
// and progressive locking for correct repeated-eigenvalue discovery.
//
// Namespace: vcp::tsparse_experimental
//
// Scalar types (static_assert, trl_is_supported_scalar_): a real point scalar
// with a numeric_limits specialization, or std::complex of such a type.
// EIG-T2 D-1: every inner product of the run is the Hermitian one,
// <x, y> = sum conj(x_i) y_i (tsparse_scalar::hermitian_dot); orthogonalization
// coefficients are of type T, norms are hermitian_norm_value, the Lanczos
// alpha is the real part of <v_j, A v_j>.  The projected matrix stays real
// symmetric (stored in T with zero imaginary parts) and the Ritz vectors are
// real combinations of the complex basis.  For a real T the conjugation does
// nothing and every formula is the one used before EIG-T2.
// EIG-T2 D-2: the B-inner-product slot (trl_run_b_ with a functor bmul) is for
// real scalars only; a complex T with bmul is refused at compile time.
//
// Apply functor convention:
//   void operator()(const std::vector<T>& x, std::vector<T>& y) const;
//
// options.max_iter = total matrix-vector product budget (main loop only).
//
// PRECONDITION (EIG-1 F-6, D-12): Apply must represent a REAL SYMMETRIC
// (complex T: HERMITIAN) operator.  The Apply abstraction cannot verify symmetry and nothing after
// the run checks it (EIG-T1 T-9): with a nonsymmetric operator the returned
// values are meaningless, as with ARPACK's dsaupd.  The dispatch layer checks
// the symmetry of the matrix before it names this kernel.
//
// Termination contract (EIG-T1 T-7a / T-8 / T-9a):
//   converged=true requires (i) k pairs locked by the estimate test
//   estimate <= tol * max(eps^(2/3) * ||A||_est, |theta|), (ii) no unconverged
//   Ritz candidate certainly inside the returned set at termination
//   (honest_termination_check_) in a subspace built after the last change of
//   the returned set, and (iii) the repeated-eigenvalue conditions of T-8.
//   The exact residuals in the result are diagnostics; they do not decide.
//   Note (EIG-0 C-4): converged=true is NOT a completeness guarantee (Krylov
//   methods cannot see eigenspaces orthogonal to the start vector).

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
#include <vcp/tsparse/tsparse_eigs_types.hpp>

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
    // EIG-T1 P5a: B-inner-product run only -- a fresh value z^T B z was
    // negative or not finite (B not positive semidefinite on the Krylov
    // space); the run stopped without a verdict.
    bool b_inner_invalid;
    thick_restart_lanczos_result()
        : restart_count(0), matrix_vector_products(0), locked_count(0),
          b_inner_invalid(false) {}
};

// ---------------------------------------------------------------------------
// Type guard.
// ---------------------------------------------------------------------------
// EIG-T1 T-1f: every real point scalar whose real type has a numeric_limits
// specialization (float / double / long double / kv::dd / kv::mpfr<N>).
// Complex scalars and types without numeric_limits (kv::interval) stay out.
template <class T>
struct trl_is_real_floating {
    static const bool value =
        !vcp::tsparse_scalar::is_complex<T>::value &&
        !vcp::tsparse_scalar::numeric_limits_is_unspecialized<
            typename vcp::tsparse_scalar::real_type<T>::type>::value;
};

// EIG-T2 D-1: the scalar types the kernel itself accepts -- the real types of
// trl_is_real_floating and their std::complex.  trl_is_real_floating keeps
// its meaning (the dispatch layer uses it to select the real-only routes).
template <class T>
struct trl_is_supported_scalar_ {
    static const bool value =
        !vcp::tsparse_scalar::numeric_limits_is_unspecialized<
            typename vcp::tsparse_scalar::real_type<T>::type>::value;
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
                // EIG-T2 D-1: Hermitian coefficient <e_j, v> (type T)
                const T c = vcp::tsparse_scalar::hermitian_dot((*extra_vecs)[j], v);
                for (std::size_t i = 0; i < n; i++) v[i] -= c * (*extra_vecs)[j][i];
            }
        }
    }
    return vcp::tsparse_scalar::hermitian_norm_value(v);
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
            // EIG-T2 D-1: Hermitian coefficient <v_j, z> (type T)
            const T c = vcp::tsparse_scalar::hermitian_dot(basis[j], z);
            for (std::size_t i = 0; i < n; i++) z[i] -= c * basis[j][i];
        }
    }
    return vcp::tsparse_scalar::hermitian_norm_value(z);
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

// ---------------------------------------------------------------------------
// EIG-T1 P5a: B-inner-product slot (ARPACK mode 3 shape: the operator
// OP = (A - sigma B)^{-1} B is symmetric in the B inner product).
//
// trl_run_b_ takes a pointer to a functor  bmul(x, y): y = B x.  A null
// pointer selects the Euclidean inner product: every formula of the run is
// then the one used before this slot existed (the B branches are not
// entered).  With B the run keeps B * v beside every basis / Ritz / locked
// vector, so a B inner product against a kept vector is one Euclidean dot;
// only the norm of a new Lanczos vector needs a multiplication by B.
// ---------------------------------------------------------------------------
template <class T>
struct trl_no_b_ {
    void operator()(const std::vector<T>&, std::vector<T>&) const {}
};

// z <- z - sum_j <v_j, z>_B v_j over locked (any B-norm) and basis
// (B-orthonormal); two passes.  bbasis[j] = B basis[j], blocked[i] = B q_i.
template <class T>
void b_project_(
    std::vector<T>& z,
    const std::vector<std::vector<T> >& basis,
    const std::vector<std::vector<T> >& bbasis,
    const std::vector<vcp::tsparse::locked_pair<T> >& locked,
    const std::vector<std::vector<T> >& blocked)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = z.size();
    for (int pass = 0; pass < 2; pass++) {
        for (std::size_t li = 0; li < locked.size() && li < blocked.size(); li++) {
            const std::vector<T>& q = locked[li].vector;
            if (q.size() != n) continue;
            const R q2 = vcp::tsparse_scalar::real_dot_value(q, blocked[li]);
            if (q2 > std::numeric_limits<R>::epsilon()) {
                const T c = T(vcp::tsparse_scalar::real_dot_value(blocked[li], z) / q2);
                for (std::size_t i = 0; i < n; i++) z[i] -= c * q[i];
            }
        }
        for (std::size_t j = 0; j < basis.size() && j < bbasis.size(); j++) {
            const T c = T(vcp::tsparse_scalar::real_dot_value(bbasis[j], z));
            for (std::size_t i = 0; i < n; i++) z[i] -= c * basis[j][i];
        }
    }
}

// v <- v - sum_i <q_i, v>_B q_i over locked, with bv = B v updated by the
// same coefficients; two passes.  Returns v^T B v (from the updated pair).
template <class T>
typename vcp::tsparse_scalar::real_type<T>::type b_deflate_locked_(
    std::vector<T>& v,
    std::vector<T>& bv,
    const std::vector<vcp::tsparse::locked_pair<T> >& locked,
    const std::vector<std::vector<T> >& blocked)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = v.size();
    for (int pass = 0; pass < 2; pass++) {
        for (std::size_t li = 0; li < locked.size() && li < blocked.size(); li++) {
            const std::vector<T>& q = locked[li].vector;
            if (q.size() != n) continue;
            const R q2 = vcp::tsparse_scalar::real_dot_value(q, blocked[li]);
            if (q2 > std::numeric_limits<R>::epsilon()) {
                const T c = T(vcp::tsparse_scalar::real_dot_value(blocked[li], v) / q2);
                for (std::size_t i = 0; i < n; i++) {
                    v[i]  -= c * q[i];
                    bv[i] -= c * blocked[li][i];
                }
            }
        }
    }
    return vcp::tsparse_scalar::real_dot_value(v, bv);
}

// true when x is a usable squared B-norm (not negative, finite)
template <class R>
bool b_value_ok_(const R& x)
{
    if (!(x >= R(0))) return false;
    if (!(x <= (std::numeric_limits<R>::max)())) return false;
    return true;
}

// bz = B z, nrm = sqrt(z^T B z).  false: the value is negative or not finite.
template <class BMul, class T>
bool b_norm_(
    const BMul& bmul,
    const std::vector<T>& z,
    std::vector<T>& bz,
    typename vcp::tsparse_scalar::real_type<T>::type& nrm)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    bz.clear();
    bmul(z, bz);
    if (bz.size() != z.size()) return false;
    const R s2 = vcp::tsparse_scalar::real_dot_value(z, bz);
    if (!b_value_ok_<R>(s2)) return false;
    nrm = vcp::tsparse_scalar::sqrt_value(s2);
    return true;
}

// Start vector of a B run (ARPACK dgetv0 with bmat = 'G'): a deterministic
// vector goes through OP once (it then lies in the range of OP; the
// application is counted), is B-orthogonalized against locked and extra and
// B-normalized.  Returns false when no vector was produced: `invalid` set =
// the B inner product failed; otherwise the budget ended or every attempt
// collapsed.
template <class Apply, class BMul, class T>
bool b_start_vector_(
    const Apply& apply,
    const BMul& bmul,
    std::size_t n,
    unsigned int& seed,
    const std::vector<vcp::tsparse::locked_pair<T> >& locked,
    const std::vector<std::vector<T> >& blocked,
    const std::vector<std::vector<T> >& extra,
    const std::vector<std::vector<T> >& bextra,
    std::size_t& mv_count,
    const std::size_t max_mv,
    std::vector<T>& v,
    std::vector<T>& bv,
    bool& invalid)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const R eps_n = std::numeric_limits<R>::epsilon() * R(n + 1);
    for (int attempt = 0; attempt < 60; attempt++) {
        if (mv_count >= max_mv) return false;
        const std::vector<T> s0 =
            vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
        std::vector<T> w;
        apply(s0, w);
        mv_count++;
        if (w.size() != n) { invalid = true; return false; }
        R nrm0(0);
        std::vector<T> bw;
        if (!b_norm_<BMul, T>(bmul, w, bw, nrm0)) { invalid = true; return false; }
        b_project_<T>(w, extra, bextra, locked, blocked);
        R nrm(0);
        if (!b_norm_<BMul, T>(bmul, w, bw, nrm)) { invalid = true; return false; }
        if (nrm > R(0) && nrm > eps_n * nrm0) {
            v.resize(n);
            bv.resize(n);
            for (std::size_t i = 0; i < n; i++) { v[i] = w[i] / T(nrm); bv[i] = bw[i] / T(nrm); }
            return true;
        }
    }
    return false;
}

} // namespace trl_detail

// ===========================================================================
// thick_restart_lanczos_eigs_with_diagnostics
//
// EIG-T1 T-2: the body is trl_run_; evaluate_residuals = false (used by the
// dispatch layer's routes that work on a transformed operator and compute
// the residuals of the ORIGINAL problem themselves) skips the end-of-run
// applications, so that the number of applications never exceeds max_iter.
// The residual fields then hold the estimates.  The public entries below
// always pass true (the behaviour before EIG-T1).
// ===========================================================================
//
// EIG-T1 P5a: trl_run_b_ is the body with the B-inner-product slot (see
// trl_detail above); trl_run_ passes a null pointer (Euclidean).
template <class Apply, class BMul, class T>
thick_restart_lanczos_result<T> trl_run_b_(
    const Apply& apply,
    const BMul* bmul,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options,
    const bool evaluate_residuals)
{
    static_assert(
        trl_is_supported_scalar_<T>::value,
        "thick_restart_lanczos_eigs: T must be a real floating-point type "
        "or std::complex of one.");
    // EIG-T2 D-2: the B-inner-product slot is not used with complex scalars
    static_assert(
        !vcp::tsparse_scalar::is_complex<T>::value
        || (std::is_same<BMul, trl_detail::trl_no_b_<T> >::value),
        "thick_restart_lanczos_eigs: the B-inner-product slot (bmul) is "
        "for real scalars only.");

    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    thick_restart_lanczos_result<T> diag;
    vcp::eig_result<T>& result = diag.eigs;
    diag.restart_count          = 0;
    diag.matrix_vector_products = 0;
    diag.locked_count           = 0;
    diag.b_inner_invalid        = false;

    // EIG-T1 P5a: B-inner-product state (used only when bmul != 0)
    const bool useB = (bmul != 0);
    std::vector<std::vector<T> > bbasis;    // B * basis[j]
    std::vector<std::vector<T> > blocked;   // B * locked[i].vector
    std::vector<T> bz_overflow;             // B * z_overflow
    bool b_invalid = false;                 // a fresh z^T B z was negative / not finite
    bool b_no_start = false;                // no start vector could be produced

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

    // EIG-T1 T-7a: eps_R^(2/3) by arithmetic only (eps_R = 2^-p is a power of
    // two: p = number of doublings that bring eps_R to 1; the value is
    // 0.5^floor(2p/3)).  No pow / cbrt, so the same code serves every real
    // scalar type with a numeric_limits specialization.
    R eps23 = R(1);
    {
        std::size_t p_bits = 0;
        R e = std::numeric_limits<R>::epsilon();
        while (e < R(1) && p_bits < std::size_t(100000)) { e = e * R(2); p_bits++; }
        const R half = R(1) / R(2);
        for (std::size_t i = 0; i < (2 * p_bits) / 3; i++) eps23 = eps23 * half;
    }

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
    // EIG-T1 T-8: pairs locked outside the returned prefix since the prefix
    // last changed (see the termination conditions in the main loop).
    std::size_t t8_outside_locks = 0;
    // EIG-T1 T-8: the returned values (ascending) as of the last basis
    // (re)build, to tell a real change of the returned prefix from an
    // exchange among copies of one repeated eigenvalue.
    std::vector<R> t8_prefix_values_at_restart;
    // EIG-T1 T-9a: the unlocked Ritz pairs of the last completed cycle (target
    // order).  When the run ends without convergence the best of them fill
    // the returned set up to k as diagnostics (approximations with their
    // exact residuals; converged_count counts the locked pairs only).
    std::vector<vcp::tsparse::ritz_pair<T> > t9_last_active;

    // -----------------------------------------------------------------------
    // Initial starting vector
    // -----------------------------------------------------------------------
    if (useB) {
        std::vector<std::vector<T> > empty;
        std::vector<T> v0, bv0;
        if (trl_detail::b_start_vector_<Apply, BMul, T>(apply, *bmul, n, seed,
                locked, blocked, empty, empty, mv_count, max_mv, v0, bv0, b_invalid)) {
            basis.push_back(v0);
            bbasis.push_back(bv0);
        } else {
            b_no_start = true;
        }
    } else {
        std::vector<std::vector<T> > empty;
        basis.push_back(
            trl_detail::new_start_vector<T>(n, seed, locked, empty));
    }

    // -----------------------------------------------------------------------
    // Main loop
    // -----------------------------------------------------------------------
    // EIG-1 F-3-1: keep iterating while the C-2 check holds termination back,
    // even though locked.size() >= k (the lock set is kept, not discarded).
    while ((locked.size() < k || honest_hold) && mv_count < max_mv
           && !b_invalid && !b_no_start) {

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

            // EIG-T2 D-1: alpha = Re <v_j, A v_j> (Hermitian inner product)
            const R alpha_j_r = useB
                ? vcp::tsparse_scalar::real_dot_value(bbasis[j], z)
                : vcp::tsparse_scalar::real_part(
                      vcp::tsparse_scalar::hermitian_dot(vj, z));
            const T alpha_j   = T(alpha_j_r);
            for (std::size_t i = 0; i < n; i++) z[i] -= alpha_j * vj[i];

            R nrm(0);
            std::vector<T> bz;   // B * z (B run only)
            if (useB) {
                trl_detail::b_project_<T>(z, basis, bbasis, locked, blocked);
                if (!trl_detail::b_norm_<BMul, T>(*bmul, z, bz, nrm)) { b_invalid = true; break; }
            } else {
                nrm = trl_detail::reorthogonalize(z, basis, locked);
            }

            alpha.push_back(alpha_j);

            // EIG-T1 T-7a (iii): the breakdown threshold is relative to the
            // operator scale (running norm estimate; inside a cycle also the
            // largest |alpha| / |beta| seen so far, because the estimate of
            // the first cycle is not available before its projected matrix).
            R brk_tol = small_tol;
            {
                R sc = anorm_est_run;
                for (std::size_t i2 = 0; i2 < alpha.size(); i2++) {
                    const R a2 = vcp::tsparse_scalar::abs_value(
                        vcp::tsparse_scalar::real_part(alpha[i2]));
                    if (a2 > sc) sc = a2;
                }
                for (std::size_t i2 = 0; i2 < beta.size(); i2++) {
                    const R b2 = vcp::tsparse_scalar::abs_value(
                        vcp::tsparse_scalar::real_part(beta[i2]));
                    if (b2 > sc) sc = b2;
                }
                brk_tol = small_tol * sc;
            }

            if (j + 1 >= m_limit || nrm <= brk_tol) {
                beta_overflow = nrm;
                if (nrm > brk_tol) {
                    for (std::size_t i = 0; i < n; i++) z_overflow[i] = z[i] / T(nrm);
                    if (useB) {
                        bz_overflow.resize(n);
                        for (std::size_t i = 0; i < n; i++) bz_overflow[i] = bz[i] / T(nrm);
                    }
                } else {
                    std::fill(z_overflow.begin(), z_overflow.end(), T(0));
                    if (useB) bz_overflow.assign(n, T(0));
                }
                break;
            }

            beta.push_back(T(nrm));
            std::vector<T> vnew(n);
            for (std::size_t i = 0; i < n; i++) vnew[i] = z[i] / T(nrm);
            basis.push_back(vnew);
            if (useB) {
                for (std::size_t i = 0; i < n; i++) bz[i] = bz[i] / T(nrm);
                bbasis.push_back(bz);
            }
        }
        if (b_invalid) break;

        // EIG-T1 T-9a: when the budget ends the expansion before the subspace
        // is full, the last step has already produced its residual (beta and
        // the next basis vector).  They are the overflow of this shorter
        // cycle; without this the estimates of the cycle would be computed
        // with beta_overflow = 0 and every Ritz pair would count as converged.
        if (!alpha.empty()
            && beta.size() == alpha.size() && basis.size() == alpha.size() + 1) {
            beta_overflow = vcp::tsparse_scalar::abs_value(
                vcp::tsparse_scalar::real_part(beta.back()));
            z_overflow = basis.back();
            beta.pop_back();
            basis.pop_back();
            if (useB && bbasis.size() == basis.size() + 1) {
                bz_overflow = bbasis.back();
                bbasis.pop_back();
            }
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

        // EIG-T1 T-1a: the projected problem is solved to
        // eps * (subspace dimension + 1) * ||T_proj||_inf (relative to the
        // projected matrix itself), not to an absolute value tied to n.
        R proj_tol = small_tol;
        {
            R tnorm = R(0);
            for (std::size_t i2 = 0; i2 < T_proj.size(); i2++) {
                R rowsum(0);
                for (std::size_t j2 = 0; j2 < T_proj[i2].size(); j2++)
                    rowsum += vcp::tsparse_scalar::abs_value(T_proj[i2][j2]);
                if (rowsum > tnorm) tnorm = rowsum;
            }
            proj_tol = std::numeric_limits<R>::epsilon()
                * R(static_cast<int>(m_actual + 1)) * tnorm;
        }

        vcp::tsparse_projected::projected_eigensolver_result<T> proj_sym =
            vcp::tsparse_projected::solve_real_symmetric_projected(
                T_proj, proj_iter, proj_tol);

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
        // EIG-T1 T-1c: keep max(k + max(k/2, 2), subspace / 2) Ritz pairs
        // over a restart (the upper limit m_limit - 2 below is unchanged).
        const std::size_t k_want =
            std::min(n_ritz, std::max<std::size_t>(k + k_extra, m_limit / 2));

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
            std::vector<T> bvector;   // B * vector (B run only)
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
            std::vector<T> bu;
            R unrm(0);
            if (useB) {
                // B run: the same lift for B u; B-orthogonalize against locked.
                // u^T B u comes from the lifted pair (no fresh product): a
                // value that is not positive means a collapsed vector here.
                bu = trl_detail::lift_ritz(bbasis, proj_sym.pairs[idx].vector, n);
                const R u2 = trl_detail::b_deflate_locked_<T>(u, bu, locked, blocked);
                if (!(u2 > R(0)) || !trl_detail::b_value_ok_<R>(u2)) continue;
                unrm = vcp::tsparse_scalar::sqrt_value(u2);
            } else {
                unrm = vcp::tsparse::orthogonalize_against_locked(u, locked, false);
            }
            if (unrm <= small_tol) continue;
            for (std::size_t i = 0; i < n; i++) u[i] /= T(unrm);
            if (useB) for (std::size_t i = 0; i < n; i++) bu[i] /= T(unrm);

            ritz_data_t rd;
            rd.value     = theta;
            rd.vector    = u;
            rd.bvector   = bu;
            rd.res_abs   = res_abs_val;
            rd.res_rel   = res_rel_val;
            // EIG-T1 T-7a (i): ARPACK's test (dsconv.f: bound <= tol *
            // max(eps23, |theta|)) with the floor made relative to the
            // operator scale: estimate <= tol * max(eps_R^(2/3) * anorm_est_run,
            // |theta|).  No absolute term, so the verdict does not depend on
            // the scale of the matrix.
            {
                const R floor_t7 = eps23 * anorm_est_run;
                rd.converged = (res_abs_val
                    <= tol * ((floor_t7 > abs_theta) ? floor_t7 : abs_theta));
            }
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
            if (useB) {
                // lock_converged_pairs appended the first n_locked_now
                // converged entries in this order
                std::size_t taken = 0;
                for (std::size_t i = 0; i < ritz_sel.size() && taken < n_locked_now; i++) {
                    if (!ritz_sel[i].converged) continue;
                    blocked.push_back(ritz_sel[i].bvector);
                    taken++;
                }
            }
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
            bool prefix_fresh = (cur_prefix == locked_prefix_at_restart);
            // EIG-T1 T-8: when a newly locked pair only takes the place of an
            // equal value in the returned prefix (another copy of a repeated
            // eigenvalue at the edge of the prefix), the returned VALUES are
            // the same and the prefix counts as unchanged.  Equal = within
            // the convergence tolerance of the values themselves.
            if (!prefix_fresh
                && cur_prefix.size() == t8_prefix_values_at_restart.size()
                && !cur_prefix.empty()) {
                std::vector<R> cur_vals;
                cur_vals.reserve(cur_prefix.size());
                for (std::size_t pi = 0; pi < cur_prefix.size(); pi++)
                    cur_vals.push_back(vcp::tsparse_scalar::real_part(locked[cur_prefix[pi]].value));
                std::sort(cur_vals.begin(), cur_vals.end());
                bool same_values = true;
                for (std::size_t pi = 0; pi < cur_vals.size() && same_values; pi++) {
                    const R a8 = vcp::tsparse_scalar::abs_value(cur_vals[pi]);
                    const R f8 = eps23 * anorm_est_run;
                    const R d8 = vcp::tsparse_scalar::abs_value(
                        R(cur_vals[pi] - t8_prefix_values_at_restart[pi]));
                    if (!(d8 <= tol * ((f8 > a8) ? f8 : a8))) same_values = false;
                }
                if (same_values) prefix_fresh = true;
            }
            const bool no_inner_unconverged =
                vcp::tsparse::honest_termination_check_(
                    active_rp, locked, k, target, shift);
            honest_hold = !(no_inner_unconverged && prefix_fresh);
            // ---------------------------------------------------------------
            // EIG-T1 T-8: a hidden copy of a repeated eigenvalue needs more
            // than one verification cycle to surface.  Besides the two
            // conditions above, termination needs one of
            //   (i)  a (k+1)-th pair: at least one pair was locked OUTSIDE the
            //        returned prefix after the prefix last changed (or all n
            //        pairs are locked), or
            //   (ii) the best unlocked Ritz value (first in target order) lies
            //        outside the worst returned value and its residual
            //        estimate is smaller than the distance between the two.
            // A hidden copy converges inward, fails (ii), gets locked and
            // changes the returned prefix (the verification then restarts).
            // (i) counts only locks made after the last change of the prefix:
            // pairs locked earlier come from a subspace that cannot contain a
            // hidden copy, so a bare count of k + 1 locked pairs would let a
            // verification end after a single cycle (triple eigenvalue /
            // several pairs locked in one cycle).
            // Without any unlocked candidate neither holds and the iteration
            // continues (budget exhaustion ends as max_iter_exhausted).
            // ---------------------------------------------------------------
            if (!prefix_fresh) t8_outside_locks = 0;
            else t8_outside_locks += n_locked_now;
            if (!honest_hold) {
                bool t8_ok = (locked.size() >= n || t8_outside_locks >= 1);
                if (!t8_ok && !active_rp.empty() && cur_prefix.size() >= k) {
                    typedef std::complex<R> C8;
                    // worst returned value = the last of the target-ordered prefix
                    R worst_key = vcp::tsparse_eigen_selection::target_distance(
                        C8(vcp::tsparse_scalar::real_part(locked[cur_prefix[0]].value), R(0)),
                        target, shift);
                    const bool prefer_large =
                        (target == eig_target::largest_magnitude ||
                         target == eig_target::largest_algebraic);
                    for (std::size_t pi = 1; pi < cur_prefix.size(); pi++) {
                        const R key = vcp::tsparse_eigen_selection::target_distance(
                            C8(vcp::tsparse_scalar::real_part(locked[cur_prefix[pi]].value), R(0)),
                            target, shift);
                        if (prefer_large ? (key < worst_key) : (key > worst_key)) worst_key = key;
                    }
                    // active_rp keeps the target order of the selection; its
                    // first entry is the best unlocked candidate.
                    const R best_key = vcp::tsparse_eigen_selection::target_distance(
                        C8(vcp::tsparse_scalar::real_part(active_rp[0].value), R(0)),
                        target, shift);
                    const bool outside = prefer_large ? (best_key < worst_key)
                                                      : (best_key > worst_key);
                    if (outside) {
                        const R gap = prefer_large ? R(worst_key - best_key)
                                                   : R(best_key - worst_key);
                        if (active_rp[0].residual_absolute < gap) t8_ok = true;
                    }
                }
                if (!t8_ok) honest_hold = true;
            }
            // EIG-T1 T-8: with all n pairs locked nothing can be hidden; no
            // further verification cycle is needed.
            if (locked.size() >= n) honest_hold = false;
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

        t9_last_active.swap(active_rp);   // not read below

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
        std::vector<std::vector<T> > new_bbasis;   // B run only
        std::vector<T>               new_alpha;
        std::vector<T>               new_beta;
        new_basis.reserve(k_ret + 2);
        new_alpha.reserve(k_ret + 1);
        new_beta.reserve(k_ret + 1);

        for (std::size_t i = 0; i < k_ret; i++) {
            new_basis.push_back(retained[i].vector);
            new_alpha.push_back(retained[i].value);
            if (useB) new_bbasis.push_back(retained[i].bvector);
        }
        // B run: a fresh start vector (through OP, B-orthogonal to locked and
        // the retained vectors).  false = none was produced; the run ends.
        struct b_fresh_t {
            static bool make(const Apply& ap, const BMul* bm, std::size_t nn, unsigned int& sd,
                             const std::vector<vcp::tsparse::locked_pair<T> >& lk,
                             const std::vector<std::vector<T> >& blk,
                             std::vector<std::vector<T> >& nb,
                             std::vector<std::vector<T> >& nbb,
                             std::size_t& mv, const std::size_t mx, bool& inv) {
                const std::vector<std::vector<T> > extra(nb), bextra(nbb);
                std::vector<T> v1, bv1;
                if (!trl_detail::b_start_vector_<Apply, BMul, T>(ap, *bm, nn, sd,
                        lk, blk, extra, bextra, mv, mx, v1, bv1, inv)) return false;
                nb.push_back(v1);
                nbb.push_back(bv1);
                return true;
            }
        };

        // EIG-T1 T-7a (iii): coupling threshold relative to the operator scale
        const R couple_tol = small_tol * anorm_est_run;
        if (k_ret > 0 && beta_overflow > couple_tol && !fresh_restart_needed) {
            std::vector<T> v_new = z_overflow;
            std::vector<T> bv_new;   // B run only

            // Orthogonalize v_new against locked using Phase 1 helper
            R vnew_nrm(0);
            if (useB) {
                bv_new = bz_overflow;
                const R v2 = trl_detail::b_deflate_locked_<T>(v_new, bv_new, locked, blocked);
                vnew_nrm = (v2 > R(0) && trl_detail::b_value_ok_<R>(v2))
                    ? vcp::tsparse_scalar::sqrt_value(v2) : R(0);
            } else {
                vcp::tsparse::orthogonalize_against_locked(v_new, locked, false);
                vnew_nrm = vcp::tsparse_scalar::hermitian_norm_value(v_new);
            }
            if (vnew_nrm > small_tol) {
                for (std::size_t i = 0; i < n; i++) v_new[i] /= T(vnew_nrm);
                new_basis.push_back(v_new);
                if (useB) {
                    for (std::size_t i = 0; i < n; i++) bv_new[i] /= T(vnew_nrm);
                    new_bbasis.push_back(bv_new);
                }

                // Wu & Simon coupling: c_j = beta_overflow * y_j[m-1]
                for (std::size_t j = 0; j < k_ret; j++) {
                    new_beta.push_back(T(beta_overflow) * retained[j].proj_last);
                }
            } else if (useB) {
                if (!b_fresh_t::make(apply, bmul, n, seed, locked, blocked,
                        new_basis, new_bbasis, mv_count, max_mv, b_invalid)) { b_no_start = true; break; }
                for (std::size_t j = 0; j < k_ret; j++) new_beta.push_back(T(0));
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
            if (useB) {
                if (!b_fresh_t::make(apply, bmul, n, seed, locked, blocked,
                        new_basis, new_bbasis, mv_count, max_mv, b_invalid)) { b_no_start = true; break; }
                for (std::size_t j = 0; j < k_ret; j++) new_beta.push_back(T(0));
            } else {
            std::vector<std::vector<T> > extra;
            extra.reserve(k_ret);
            for (std::size_t i = 0; i < k_ret; i++)
                extra.push_back(retained[i].vector);
            new_basis.push_back(
                trl_detail::new_start_vector<T>(n, seed, locked, extra));
            for (std::size_t j = 0; j < k_ret; j++) new_beta.push_back(T(0));
            }
        } else if (useB) {
            if (!b_fresh_t::make(apply, bmul, n, seed, locked, blocked,
                    new_basis, new_bbasis, mv_count, max_mv, b_invalid)) { b_no_start = true; break; }
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
        if (useB) bbasis.swap(new_bbasis);
        alpha     = new_alpha;
        beta      = new_beta;
        k_restart = k_ret;

        // EIG-1 F-3-1 freshness guard: the basis just built is (and will be
        // kept, via reorthogonalize) orthogonal to the current locked set;
        // record which returned prefix that subspace knows about.
        locked_prefix_at_restart =
            trl_detail::locked_prefix_indices(locked, k, target, shift);
        t8_prefix_values_at_restart.clear();
        for (std::size_t pi = 0; pi < locked_prefix_at_restart.size(); pi++)
            t8_prefix_values_at_restart.push_back(vcp::tsparse_scalar::real_part(
                locked[locked_prefix_at_restart[pi]].value));
        std::sort(t8_prefix_values_at_restart.begin(), t8_prefix_values_at_restart.end());

        result.iterations++;
        total_restarts++;

    } // end main while

    // -----------------------------------------------------------------------
    // Budget flag (captured before re-evaluation calls)
    // -----------------------------------------------------------------------
    const bool budget_exhausted =
        (mv_count >= max_mv) && (locked.size() < k || honest_hold);
    // EIG-T1 P5a: a B run that stopped on its inner product has no verdict
    diag.b_inner_invalid = b_invalid;
    if (b_invalid && result.failure_reason.empty())
        result.failure_reason = "B inner product negative or not finite";
    if (b_no_start && !b_invalid && !budget_exhausted && result.failure_reason.empty())
        result.failure_reason = "no start vector in the range of the operator";

    // EIG-T1 T-9a: number of pairs the solver itself judged converged.  A run
    // that ends short of k locked pairs returns the best unlocked Ritz pairs
    // of its last cycle after them (diagnostics; never counted as converged).
    const std::size_t t9_true_locked = locked.size();
    if (locked.size() < k) {
        for (std::size_t i = 0; i < t9_last_active.size() && locked.size() < k; i++) {
            vcp::tsparse::locked_pair<T> lp;
            lp.value = t9_last_active[i].value;
            lp.vector = t9_last_active[i].vector;
            lp.residual_absolute = t9_last_active[i].residual_absolute;
            lp.residual_relative = t9_last_active[i].residual_relative;
            if (lp.vector.size() != n) continue;
            locked.push_back(lp);
        }
    }

    // -----------------------------------------------------------------------
    // Re-evaluate actual residuals via apply(): ||Av - lambda*v||
    // Replaces Wu & Simon estimates with exact operator residuals.
    // -----------------------------------------------------------------------
    for (std::size_t i = 0; evaluate_residuals && i < locked.size(); i++) {
        std::vector<T> w;
        apply(locked[i].vector, w);
        mv_count++;
        const T lam = locked[i].value;
        std::vector<T> diff(n);
        for (std::size_t j = 0; j < n; j++) diff[j] = w[j] - lam * locked[i].vector[j];
        const R res_abs = vcp::tsparse_scalar::hermitian_norm_value(diff);
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
    diag.locked_count             = t9_true_locked;
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
    result.converged_count        = std::min(result.returned_count, t9_true_locked);

    // Convergence is based on the clamped k (not k_original).
    //
    // EIG-T1 T-9a: the verdict is the solver's own -- k pairs locked by the
    // estimate test and the termination conditions of the main loop.  The
    // exact residuals computed above stay in the result as diagnostics and
    // do not decide `converged`.
    if (result.returned_count >= k && !honest_hold && t9_true_locked >= k
        && !b_invalid && !b_no_start) {
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

template <class Apply, class T>
thick_restart_lanczos_result<T> trl_run_(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options,
    const bool evaluate_residuals)
{
    return trl_run_b_<Apply, trl_detail::trl_no_b_<T>, T>(
        apply, static_cast<const trl_detail::trl_no_b_<T>*>(0), n, k, options,
        evaluate_residuals);
}

// EIG-T1 P5a: the kernel in the B inner product (bmul(x, y): y = B x), without
// the end-of-run residual applications.  The returned vectors are
// B-orthonormal.  diag.b_inner_invalid tells the caller that the inner
// product failed (the caller then takes another route).
template <class Apply, class BMul, class T>
thick_restart_lanczos_result<T> thick_restart_lanczos_eigs_b_inner_(
    const Apply& apply,
    const BMul& bmul,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    return trl_run_b_<Apply, BMul, T>(apply, &bmul, n, k, options, false);
}

template <class Apply, class T>
thick_restart_lanczos_result<T> thick_restart_lanczos_eigs_with_diagnostics(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    return trl_run_<Apply, T>(apply, n, k, options, true);
}

// EIG-T1 T-2: the kernel without the end-of-run residual applications (see
// trl_run_); for callers that evaluate the residuals of their own problem.
template <class Apply, class T>
vcp::eig_result<T> thick_restart_lanczos_eigs_no_residual_(
    const Apply& apply,
    std::size_t n,
    std::size_t k,
    const vcp::eig_options<T>& options)
{
    return trl_run_<Apply, T>(apply, n, k, options, false).eigs;
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
