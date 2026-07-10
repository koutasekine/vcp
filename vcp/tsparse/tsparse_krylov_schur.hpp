// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// vcp/tsparse/tsparse_krylov_schur.hpp
//
// EIG-3 T-1: Krylov-Schur eigensolver for REAL GENERAL operators, rebuilt on
// the EIG-2 real Schur QR core (vcp::tsparse_real_schur).  This file is a
// full replacement of the Phase-3 experimental driver (whose projected-QR
// dependencies carried the D-QR/D-4/D-7/D-8/D-15 defect family); it shares
// no code with the legacy core (tsparse_eigensolvers.hpp projected QR).
//
// Contract (EIG-0 + EIG-3 design D3-1..D3-5; frozen unit set ks_units.cpp):
//  * options.max_iter > 0  = total matrix-vector product budget, STRICT:
//    every apply() call is counted (B-24, 1:1) and mv <= max_iter always
//    (the final C-1 re-evaluation applies are inside the budget too).
//    options.max_iter == 0 = default budget 300*n (TRL convention).
//    Budget reached => honest not_converged, status "max_iter_exhausted".
//  * converged=true requires, at the terminating analysis:
//      C-1  exact lambda-space residuals of all k returned pairs, re-evaluated
//           against the operator (certified acceptance res<=tol or rel<=tol),
//      C-2  no certainly-inner unconverged Ritz candidate in the exported
//           candidate set (shared check honest_termination_check_complex_,
//           target metric per D3-3: algebraic=real part, magnitude=hypot),
//      D3-2 no complex pair inside the k return window (real-only eig_result
//           cannot return complex pairs as converged; they are reported as
//           diagnostics with an explicit reason),
//      freshness: the declaring analysis happened on a subspace enriched by
//           deterministic probe directions after the previous analysis first
//           satisfied C-1/C-2 (verification cycle; blind-spot probe for
//           multiplicity discovery), or the basis spans the whole space.
//  * C-3: no value-proximity dedup anywhere; multiplicities are discovered
//    via probe-enriched verification cycles and returned as independent pairs.
//  * C-2 evidence (D3-4, EIG-1 a-5): the final candidate set (values, flags)
//    and the freshness stamp are exported on EVERY termination path with
//    n>0 && k>0 (ks_c2_evidence).  A path that could not export evidence
//    could not declare converged.
//  * Contraction (G-2.1 method A, approved 2026-07-06): Schur-form
//    reordering by adjacent block swaps (dtrexc/dlaexc style) with certified
//    swap rejection; a rejected swap never discards the target candidate --
//    the keep window is extended to the candidate's current position
//    (keep-more, condition c-1).  Rejections are counted in
//    contraction_swap_rejections (condition c-2).  2x2 complex blocks are
//    never split across the contraction boundary (D3-5 keep-together;
//    boundary extension is the conservative direction).
//  * T-generic (GT1 P1-P6): module-scalar ops only, certified gates with
//    honest-failure third branches, no SWO-dependent std algorithms on T
//    (orderings below are hand-rolled insertion sorts, partial-order
//    tolerant), no infinity sentinels.  kv::interval instantiates and fails
//    honestly by default.
//  * Determinism: no randomness; identical input -> byte-identical output.
//
// Exact zero assignments after certified-small annihilation (swap kernel,
// deflation-style flushes) follow the EIG-2 real Schur core precedent
// (backward-stability semantics, approved at EIG-2 G-1.1).

#pragma once

#ifndef VCP_TSPARSE_KRYLOV_SCHUR_HPP
#define VCP_TSPARSE_KRYLOV_SCHUR_HPP

#include <complex>
#include <cstddef>
#include <string>
#include <vector>

// eig_result / eig_options only (NOT vcp/spmatrix.hpp: this header is also
// included from the spmats_eigs dispatch layer, which spmatrix.hpp includes
// after the spmats class definition — a spmatrix.hpp include here would be
// circular; spmats_eigs_types.hpp is the self-contained type surface)
#include <vcp/spmats_base/spmats_eigs_types.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>            // deterministic_start_vector
#include <vcp/tsparse/tsparse_real_schur.hpp>          // EIG-2 core
#include <vcp/tsparse/tsparse_dense_schur_driver.hpp>  // hessenberg_reduce_ (reuse)
#include <vcp/tsparse/tsparse_honest_termination.hpp>  // shared C-2 check

namespace vcp {
namespace tsparse_experimental {

// ---------------------------------------------------------------------------
// C-2 evidence package (D3-4)
// ---------------------------------------------------------------------------
template <class T>
struct ks_c2_evidence {
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    std::vector<R> candidate_real;         // final-analysis Ritz candidates
    std::vector<R> candidate_imag;         // complex pairs adjacent (+s, -s)
    std::vector<bool> candidate_converged; // residual-bound convergence flags
    bool fresh;                            // declaring analysis was verification-fresh
    bool exported;                         // fields populated (true for n>0 && k>0)

    ks_c2_evidence() : fresh(false), exported(false) {}
};

// ---------------------------------------------------------------------------
// Diagnostic result wrapper (field-compatible superset of the old driver)
// ---------------------------------------------------------------------------
template <class T>
struct krylov_schur_result {
    vcp::eig_result<T> eigs;
    std::size_t restart_count;
    std::size_t matrix_vector_products;
    std::size_t locked_real_count;
    std::size_t locked_complex_count;
    std::size_t projected_dimension;
    std::size_t contraction_swap_rejections;   // G-2.1 (c-2)

    ks_c2_evidence<T> c2_evidence;
    // terminating Krylov-Schur factorization  A V = V S + v_next c^T
    std::vector<std::vector<T> > final_basis;       // V: M vectors (length n)
    std::vector<std::vector<T> > final_compressed;  // S: M x M
    std::vector<T> final_residual_vector;           // v_next (length n; may be empty)
    std::vector<T> final_coupling;                  // c (length M)

    krylov_schur_result()
        : restart_count(0), matrix_vector_products(0),
          locked_real_count(0), locked_complex_count(0),
          projected_dimension(0), contraction_swap_rejections(0) {}
};

// ===========================================================================
// Internal helpers
// ===========================================================================
namespace ks_detail {

// Effective subspace dimension m (same policy as the old driver).
inline std::size_t effective_m(std::size_t n, std::size_t k, std::size_t requested)
{
    if (n == 0 || k == 0) return 0;
    std::size_t m;
    if (requested == 0) {
        std::size_t a = 2 * k + 6;
        if (a < 20) a = 20;
        if (a < k + 3) a = k + 3;
        m = a;
    } else if (requested <= k) {
        std::size_t pad = (requested > 2) ? requested : std::size_t(2);
        m = k + pad;
    } else {
        m = requested;
    }
    if (m > n) m = n;
    if (k < n && m < k + 1) m = (k + 1 < n) ? k + 1 : n;
    return m;
}

// Orthogonalize w against V[0..count-1]; coefficients accumulated into h
// (length count).  passes = 1 (MGS-like single CGS pass) or 2 (CGS2).
template <typename T>
void orthogonalize_block(const std::vector<std::vector<T> >& V,
                         const std::size_t count,
                         std::vector<T>& w,
                         std::vector<T>& h,
                         const int passes)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = w.size();
    h.assign(count, T(0));
    for (int pass = 0; pass < passes; pass++) {
        for (std::size_t i = 0; i < count; i++) {
            const R c = vcp::tsparse_scalar::real_dot_value(V[i], w);
            h[i] += T(c);
            for (std::size_t ii = 0; ii < n; ii++) w[ii] -= T(c) * V[i][ii];
        }
    }
}

// Deterministic direction orthogonal to V[0..count-1]; returns false when no
// certifiably nonzero direction is found within the attempt budget.
template <typename T>
bool fresh_orthogonal_direction(const std::vector<std::vector<T> >& V,
                                const std::size_t count,
                                const std::size_t n,
                                unsigned int& seed,
                                const typename vcp::tsparse_scalar::real_type<T>::type& floor_tol,
                                std::vector<T>& out)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    for (int attempt = 0; attempt < 24; attempt++) {
        out = vcp::tsparse_lanczos::deterministic_start_vector<T>(n, seed++);
        std::vector<T> h;
        orthogonalize_block(V, count, out, h, 2);
        const R nv = vcp::tsparse_scalar::real_norm_value(out);
        if (nv > floor_tol) {
            for (std::size_t i = 0; i < n; i++) out[i] /= T(nv);
            return true;
        }
    }
    return false;
}

// Small dense linear solve (q <= 4) with certified partial pivoting.
// Returns false when no pivot is certifiably nonzero (honest failure).
template <typename T>
bool solve_small(std::vector<std::vector<T> > A, std::vector<T> b, std::vector<T>& x)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    using vcp::tsparse_scalar::abs_value;
    const std::size_t q = A.size();
    x.assign(q, T(0));
    for (std::size_t col = 0; col < q; col++) {
        std::size_t piv = col;
        R best = abs_value(A[col][col]);
        for (std::size_t r = col + 1; r < q; r++) {
            const R a = abs_value(A[r][col]);
            if (a > best) { best = a; piv = r; }
        }
        if (!(best > R(0))) return false;   // certified-only gate (P1)
        if (piv != col) { A[piv].swap(A[col]); const T tb = b[piv]; b[piv] = b[col]; b[col] = tb; }
        for (std::size_t r = col + 1; r < q; r++) {
            const T f = A[r][col] / A[col][col];
            for (std::size_t cc = col; cc < q; cc++) A[r][cc] -= f * A[col][cc];
            b[r] -= f * b[col];
        }
    }
    for (std::size_t ri = q; ri-- > 0;) {
        T s = b[ri];
        for (std::size_t cc = ri + 1; cc < q; cc++) s -= A[ri][cc] * x[cc];
        x[ri] = s / A[ri][ri];
    }
    return true;
}

// Block boundaries of a quasi-upper-triangular T (structural subdiagonal).
template <typename T>
void enumerate_blocks(const std::vector<std::vector<T> >& Tm,
                      std::vector<std::size_t>& pos,
                      std::vector<std::size_t>& size)
{
    using vcp::tsparse_scalar::abs_value;
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t m = Tm.size();
    pos.clear(); size.clear();
    std::size_t p = 0;
    while (p < m) {
        if (p + 1 < m && !(abs_value(Tm[p + 1][p]) <= R(0))) {
            pos.push_back(p); size.push_back(2); p += 2;
        } else {
            pos.push_back(p); size.push_back(1); p += 1;
        }
    }
}

// Certified eigenvalues of the 2x2 block at (p,p).  Returns:
//   +1 certified complex pair (re +/- im, im > 0),
//    0 certified real pair (values in re1/re2),
//   -1 not certifiable (interval 0-straddle) -> honest failure upstream.
template <typename T>
int eig_2x2(const std::vector<std::vector<T> >& Tm, const std::size_t p,
            T& re, T& im, T& re1, T& re2)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    using vcp::tsparse_scalar::abs_value;
    using vcp::tsparse_scalar::sqrt_value;
    R bs = abs_value(Tm[p][p]);
    {
        const R t1 = abs_value(Tm[p][p + 1]);
        const R t2 = abs_value(Tm[p + 1][p]);
        const R t3 = abs_value(Tm[p + 1][p + 1]);
        if (t1 > bs) bs = t1;
        if (t2 > bs) bs = t2;
        if (t3 > bs) bs = t3;
    }
    if (!(bs > R(0))) return -1;
    const T a = Tm[p][p] / T(bs);
    const T b = Tm[p][p + 1] / T(bs);
    const T c = Tm[p + 1][p] / T(bs);
    const T d = Tm[p + 1][p + 1] / T(bs);
    const T tr = a + d;
    const T amd = a - d;
    const T disc = amd * amd + T(4) * (b * c);
    if (disc < T(0)) {
        const T s = sqrt_value(-disc);
        re = T(bs) * tr / T(2);
        im = T(bs) * s / T(2);
        return 1;
    }
    if (disc >= T(0)) {
        const T sd = sqrt_value(disc);
        re1 = T(bs) * (tr + sd) / T(2);
        re2 = T(bs) * (tr - sd) / T(2);
        return 0;
    }
    return -1;   // 0-straddling discriminant (interval): not certifiable
}

// Split a certified-real 2x2 block at p into two 1x1 blocks by a
// deterministic rotation (real_schur standardization pattern).  Applies the
// similarity to the whole matrix and accumulates into U.  Returns false when
// the rotation is not certifiable.
template <typename T>
bool split_real_2x2(std::vector<std::vector<T> >& Tm,
                    std::vector<std::vector<T> >& U,
                    const std::size_t p)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    using vcp::tsparse_scalar::abs_value;
    using vcp::tsparse_scalar::sqrt_value;
    const std::size_t m = Tm.size();
    R bs = abs_value(Tm[p][p]);
    {
        const R t1 = abs_value(Tm[p][p + 1]);
        const R t2 = abs_value(Tm[p + 1][p]);
        const R t3 = abs_value(Tm[p + 1][p + 1]);
        if (t1 > bs) bs = t1;
        if (t2 > bs) bs = t2;
        if (t3 > bs) bs = t3;
    }
    if (!(bs > R(0))) return false;
    const T a = Tm[p][p] / T(bs);
    const T b = Tm[p][p + 1] / T(bs);
    const T c = Tm[p + 1][p] / T(bs);
    const T d = Tm[p + 1][p + 1] / T(bs);
    const T tr = a + d;
    const T amd = a - d;
    const T disc = amd * amd + T(4) * (b * c);
    if (!(disc >= T(0))) return false;
    const T sd = sqrt_value(disc);
    const T lam1 = (amd >= T(0)) ? (tr + sd) / T(2) : (tr - sd) / T(2);
    const T v0 = lam1 - d;
    const R av0 = abs_value(v0);
    const R avc = abs_value(c);
    R vs = av0;
    if (avc > vs) vs = avc;
    if (!(vs > R(0))) return false;
    const T v0s = v0 / T(vs);
    const T cs0 = c / T(vs);
    const R vv = abs_value(v0s) * abs_value(v0s) + abs_value(cs0) * abs_value(cs0);
    if (!(vv > R(0))) return false;
    const T vn = T(sqrt_value(vv));
    const T gc = v0s / vn;
    const T gs = cs0 / vn;
    // Givens similarity on (p, p+1): rows p..m-1, cols 0..p+1, U cols
    for (std::size_t j = p; j < m; j++) {
        const T hp = Tm[p][j];
        const T hq = Tm[p + 1][j];
        Tm[p][j] = gc * hp + gs * hq;
        Tm[p + 1][j] = gc * hq - gs * hp;
    }
    for (std::size_t r = 0; r <= p + 1; r++) {
        const T hp = Tm[r][p];
        const T hq = Tm[r][p + 1];
        Tm[r][p] = gc * hp + gs * hq;
        Tm[r][p + 1] = gc * hq - gs * hp;
    }
    for (std::size_t r = 0; r < U.size(); r++) {
        const T zp = U[r][p];
        const T zq = U[r][p + 1];
        U[r][p] = gc * zp + gs * zq;
        U[r][p + 1] = gc * zq - gs * zp;
    }
    Tm[p + 1][p] = T(0);   // annihilated analytically (EIG-2 precedent)
    return true;
}

// Swap the adjacent diagonal blocks (pos p, sizes s1 then s2) of the
// quasi-triangular Tm; accumulate the orthogonal transform into U.
// Certified, deterministic; returns false = swap REJECTED (matrix unchanged).
template <typename T>
bool swap_adjacent_blocks(std::vector<std::vector<T> >& Tm,
                          std::vector<std::vector<T> >& U,
                          const std::size_t p,
                          const std::size_t s1,
                          const std::size_t s2)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    using vcp::tsparse_scalar::abs_value;
    using vcp::tsparse_scalar::sqrt_value;
    const std::size_t m = Tm.size();
    const std::size_t sz = s1 + s2;
    const R eps = vcp::tsparse_scalar::epsilon<R>();

    if (s1 == 1 && s2 == 1) {
        // dlaexc n1=n2=1: rotation from (T[p][p+1], T[p+1][p+1]-T[p][p])
        const T f = Tm[p][p + 1];
        const T g = Tm[p + 1][p + 1] - Tm[p][p];
        R sc = abs_value(f);
        {
            const R ag = abs_value(g);
            if (ag > sc) sc = ag;
        }
        if (!(sc > R(0))) {
            // equal eigenvalues, no coupling: blocks are interchangeable as-is
            return true;
        }
        const T fs = f / T(sc);
        const T gs = g / T(sc);
        const R rr = abs_value(fs) * abs_value(fs) + abs_value(gs) * abs_value(gs);
        if (!(rr > R(0))) return false;
        const T rn = T(sqrt_value(rr));
        const T cs = fs / rn;
        const T sn = gs / rn;
        for (std::size_t j = p; j < m; j++) {
            const T hp = Tm[p][j];
            const T hq = Tm[p + 1][j];
            Tm[p][j] = cs * hp + sn * hq;
            Tm[p + 1][j] = cs * hq - sn * hp;
        }
        for (std::size_t r = 0; r <= p + 1; r++) {
            const T hp = Tm[r][p];
            const T hq = Tm[r][p + 1];
            Tm[r][p] = cs * hp + sn * hq;
            Tm[r][p + 1] = cs * hq - sn * hp;
        }
        for (std::size_t r = 0; r < U.size(); r++) {
            const T zp = U[r][p];
            const T zq = U[r][p + 1];
            U[r][p] = cs * zp + sn * zq;
            U[r][p + 1] = cs * zq - sn * zp;
        }
        Tm[p + 1][p] = T(0);
        return true;
    }

    // Direct swap (dlaexc style): solve A11 X - X A22 = A12, then QR of
    // [-X; I] gives the orthogonal Q with Q^T [A11 A12; 0 A22] Q = [A22' ...].
    // Local copies of the blocks:
    std::vector<std::vector<T> > A11(s1, std::vector<T>(s1)),
                                 A22(s2, std::vector<T>(s2)),
                                 A12(s1, std::vector<T>(s2));
    for (std::size_t i = 0; i < s1; i++)
        for (std::size_t j = 0; j < s1; j++) A11[i][j] = Tm[p + i][p + j];
    for (std::size_t i = 0; i < s2; i++)
        for (std::size_t j = 0; j < s2; j++) A22[i][j] = Tm[p + s1 + i][p + s1 + j];
    for (std::size_t i = 0; i < s1; i++)
        for (std::size_t j = 0; j < s2; j++) A12[i][j] = Tm[p + i][p + s1 + j];

    // Sylvester as a small dense system over vec(X) (row-major):
    // (A11 X)_{ij} - (X A22)_{ij} = A12_{ij}
    const std::size_t q = s1 * s2;
    std::vector<std::vector<T> > M(q, std::vector<T>(q, T(0)));
    std::vector<T> rhs(q, T(0));
    for (std::size_t i = 0; i < s1; i++) {
        for (std::size_t j = 0; j < s2; j++) {
            const std::size_t row = i * s2 + j;
            rhs[row] = A12[i][j];
            for (std::size_t t = 0; t < s1; t++) M[row][t * s2 + j] += A11[i][t];
            for (std::size_t t = 0; t < s2; t++) M[row][i * s2 + t] -= A22[t][j];
        }
    }
    std::vector<T> xv;
    if (!solve_small(M, rhs, xv)) return false;   // rejected (certified)

    // Q from Householder QR of the (sz x s2) matrix [-X; I].
    std::vector<std::vector<T> > W(sz, std::vector<T>(s2, T(0)));
    for (std::size_t i = 0; i < s1; i++)
        for (std::size_t j = 0; j < s2; j++) W[i][j] = -xv[i * s2 + j];
    for (std::size_t j = 0; j < s2; j++) W[s1 + j][j] = T(1);

    std::vector<std::vector<T> > Q(sz, std::vector<T>(sz, T(0)));
    for (std::size_t i = 0; i < sz; i++) Q[i][i] = T(1);
    for (std::size_t col = 0; col < s2; col++) {
        R sc(0);
        for (std::size_t i = col; i < sz; i++) {
            const R a = abs_value(W[i][col]);
            if (a > sc) sc = a;
        }
        if (!(sc > R(0))) return false;
        std::vector<T> u(sz, T(0));
        R n2(0);
        for (std::size_t i = col; i < sz; i++) {
            u[i] = W[i][col] / T(sc);
            const R au = abs_value(u[i]);
            n2 += au * au;
        }
        if (!(n2 > R(0))) return false;
        const T nrm = T(sqrt_value(n2));
        const T sgn = (u[col] >= T(0)) ? T(1) : T(-1);
        u[col] += sgn * nrm;
        R uu(0);
        for (std::size_t i = col; i < sz; i++) {
            const R au = abs_value(u[i]);
            uu += au * au;
        }
        if (!(uu > R(0))) return false;
        const T two = T(2);
        for (std::size_t j = 0; j < s2; j++) {          // W <- P W
            T s(0);
            for (std::size_t i = col; i < sz; i++) s += u[i] * W[i][j];
            const T fct = two * s / T(uu);
            for (std::size_t i = col; i < sz; i++) W[i][j] -= fct * u[i];
        }
        for (std::size_t j = 0; j < sz; j++) {          // Q <- Q P
            T s(0);
            for (std::size_t i = col; i < sz; i++) s += Q[j][i] * u[i];
            const T fct = two * s / T(uu);
            for (std::size_t i = col; i < sz; i++) Q[j][i] -= fct * u[i];
        }
    }

    // Candidate transform of the local block D = Q^T [A11 A12; 0 A22] Q.
    std::vector<std::vector<T> > B(sz, std::vector<T>(sz, T(0)));
    for (std::size_t i = 0; i < s1; i++) {
        for (std::size_t j = 0; j < s1; j++) B[i][j] = A11[i][j];
        for (std::size_t j = 0; j < s2; j++) B[i][s1 + j] = A12[i][j];
    }
    for (std::size_t i = 0; i < s2; i++)
        for (std::size_t j = 0; j < s2; j++) B[s1 + i][s1 + j] = A22[i][j];
    std::vector<std::vector<T> > QB(sz, std::vector<T>(sz, T(0)));
    for (std::size_t i = 0; i < sz; i++)
        for (std::size_t j = 0; j < sz; j++) {
            T s(0);
            for (std::size_t t = 0; t < sz; t++) s += Q[t][i] * B[t][j];
            QB[i][j] = s;
        }
    std::vector<std::vector<T> > D(sz, std::vector<T>(sz, T(0)));
    R bnorm(0);
    for (std::size_t i = 0; i < sz; i++)
        for (std::size_t j = 0; j < sz; j++) {
            T s(0);
            for (std::size_t t = 0; t < sz; t++) s += QB[i][t] * Q[t][j];
            D[i][j] = s;
            const R a = abs_value(B[i][j]);
            if (a > bnorm) bnorm = a;
        }
    // Stability acceptance: the would-be-zero lower-left block (new trailing
    // block has size s1, leading has size s2) must be certifiably negligible.
    const R swap_tol = R(20) * eps * bnorm;
    for (std::size_t i = s2; i < sz; i++)
        for (std::size_t j = 0; j < s2; j++)
            if (!(abs_value(D[i][j]) <= swap_tol)) return false;   // rejected

    // Accept: apply Q to the full matrix and U, then install D with the
    // certified-negligible entries flushed to exact zero (EIG-2 precedent).
    // rows p..p+sz-1, columns p..m-1  (T <- Q^T T)
    std::vector<T> tmp(sz);
    for (std::size_t j = p; j < m; j++) {
        for (std::size_t i = 0; i < sz; i++) {
            T s(0);
            for (std::size_t t = 0; t < sz; t++) s += Q[t][i] * Tm[p + t][j];
            tmp[i] = s;
        }
        for (std::size_t i = 0; i < sz; i++) Tm[p + i][j] = tmp[i];
    }
    // columns p..p+sz-1, rows 0..p+sz-1  (T <- T Q)
    for (std::size_t r = 0; r < p + sz; r++) {
        for (std::size_t j = 0; j < sz; j++) {
            T s(0);
            for (std::size_t t = 0; t < sz; t++) s += Tm[r][p + t] * Q[t][j];
            tmp[j] = s;
        }
        for (std::size_t j = 0; j < sz; j++) Tm[r][p + j] = tmp[j];
    }
    for (std::size_t r = 0; r < U.size(); r++) {
        for (std::size_t j = 0; j < sz; j++) {
            T s(0);
            for (std::size_t t = 0; t < sz; t++) s += U[r][p + t] * Q[t][j];
            tmp[j] = s;
        }
        for (std::size_t j = 0; j < sz; j++) U[r][p + j] = tmp[j];
    }
    for (std::size_t i = s2; i < sz; i++)
        for (std::size_t j = 0; j < s2; j++)
            Tm[p + i][p + j] = T(0);
    return true;
}

// ---------------------------------------------------------------------------
// EIG-4 T-3 helpers (allow_complex_pairs opt-in only; unreachable when the
// option is false).  Ordering loops are hand-rolled (P6: module scalars are
// not passed to SWO-based std algorithms).
// ---------------------------------------------------------------------------

// target-order selection over (re, im) items where an item with im>0 is a
// complex pair occupying 2 slots.  Fills at least k_slots slots (keep-together:
// a straddling pair extends the fill to k_slots+1).  Returns selected indices
// in target order; empty result = not enough slots available.
template <class R>
std::vector<std::size_t> select_slot_indices_pairs(
    const std::vector<R>& re,
    const std::vector<R>& im,
    const std::size_t k_slots,
    const vcp::eig_target target,
    const R& shift)
{
    typedef std::complex<R> C;
    const std::size_t m = re.size();
    std::vector<std::size_t> out;
    if (k_slots == 0) return out;
    const bool prefer_large =
        (target == vcp::eig_target::largest_magnitude ||
         target == vcp::eig_target::largest_algebraic);
    std::vector<R> key(m, R(0));
    std::size_t total = 0;
    for (std::size_t i = 0; i < m; i++) {
        key[i] = vcp::tsparse_eigen_selection::target_distance(
            C(re[i], im[i]), target, shift);
        total += (im[i] > R(0)) ? 2 : 1;
    }
    if (total < k_slots) return out;
    std::vector<bool> used(m, false);
    std::size_t filled = 0;
    while (filled < k_slots) {
        std::size_t best = m;
        for (std::size_t i = 0; i < m; i++) {
            if (used[i]) continue;
            if (best == m) { best = i; continue; }
            const bool better = prefer_large ? (key[i] > key[best])
                                             : (key[i] < key[best]);
            if (better) best = i;
        }
        if (best == m) break;
        used[best] = true;
        out.push_back(best);
        filled += (im[best] > R(0)) ? 2 : 1;
    }
    if (filled < k_slots) { out.clear(); }
    return out;
}

// worst (k_slots-th, slot-weighted) target key over pooled (re, im) values.
// Returns false when the pool holds fewer than k_slots slots.
template <class R>
bool pool_worst_key_pairs(
    const std::vector<R>& re,
    const std::vector<R>& im,
    const std::size_t k_slots,
    const vcp::eig_target target,
    const R& shift,
    R& worst_key_out)
{
    const std::vector<std::size_t> sel =
        select_slot_indices_pairs<R>(re, im, k_slots, target, shift);
    if (sel.empty()) return false;
    typedef std::complex<R> C;
    worst_key_out = vcp::tsparse_eigen_selection::target_distance(
        C(re[sel.back()], im[sel.back()]), target, shift);
    return true;
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
    static_assert(!vcp::tsparse_scalar::is_complex<T>::value,
        "krylov_schur: complex scalar types are rejected (real general operators only)");

    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    using vcp::tsparse_scalar::abs_value;
    using vcp::tsparse_scalar::sqrt_value;
    namespace kd = ks_detail;

    krylov_schur_result<T> diag;
    vcp::eig_result<T>& result = diag.eigs;
    result.used_dense_fallback       = false;
    result.used_shift_invert         = false;
    result.used_generalized_operator = false;
    result.used_method               = "krylov_schur";
    result.used_orthogonalization    =
        (options.orthogonalization ==
         vcp::orthogonalization_method::classical_gram_schmidt_twice)
        ? "classical_gram_schmidt_twice" : "modified_gram_schmidt";

    // ---- edge cases -------------------------------------------------------
    if (k == 0) {
        result.requested_count = 0;
        result.converged = true;
        result.status = "success";
        result.message = "k=0: nothing to compute";
        return diag;
    }
    if (n == 0) {
        result.requested_count = k;
        result.converged = false;
        result.status = "failed";
        result.failure_reason = "dimension n=0 with k>0";
        return diag;
    }
    const std::size_t k_original = k;
    if (k > n) k = n;
    result.requested_count = k_original;

    // ---- parameters ---------------------------------------------------------
    const R tol = (options.tol > R(0))
        ? options.tol
        : vcp::tsparse_scalar::decimal_power_negative<R>(12);
    const std::size_t max_mv = (options.max_iter > 0)
        ? options.max_iter : std::size_t(300) * n;
    const std::size_t m_limit = kd::effective_m(n, k, options.subspace_dim);
    result.used_subspace_dim = m_limit;
    diag.projected_dimension = m_limit;

    const vcp::eig_target target = options.target;
    const R shift_val = vcp::tsparse_scalar::real_part(options.shift);
    const bool prefer_large =
        (target == vcp::eig_target::largest_magnitude ||
         target == vcp::eig_target::largest_algebraic);
    // CGS2 always: full reorthogonalization is required for the frozen
    // relation/orthonormality invariants; the option only labels the result.
    const int orth_passes = 2;
    const R eps = vcp::tsparse_scalar::epsilon<R>();
    R floor_tol = eps * R(static_cast<int>(n) + 1);
    {
        const R t4 = tol * vcp::tsparse_scalar::decimal_power_negative<R>(4);
        if (t4 > floor_tol) floor_tol = t4;
    }
    unsigned int seed = options.random_start ? options.random_seed : 0u;

    // ---- state --------------------------------------------------------------
    // Factorization: A V_M = V_M S + v_next c^T  (V: M vectors + v_next).
    std::vector<std::vector<T> > V;      // applied basis columns (M)
    std::vector<T> v_next;               // residual direction (unit or empty)
    std::vector<std::vector<T> > S;      // M x M (grows)
    std::vector<T> coup;                 // c (length M)
    std::size_t mv_count = 0;
    std::size_t restart_count = 0;
    bool budget_exhausted = false;
    bool pending_verification = false;   // in-solve probe cycle flag
    bool subspace_exhausted = false;     // no fresh orthogonal direction exists

    // ---- cross-solve pool (EIG-3 blind-spot protocol, TRL F-3-1 transplant) --
    // A single-vector Krylov factorization can represent at most one direction
    // per eigenspace, so hidden multiplicity copies of already-verified values
    // are structurally invisible to the solve that verified them.  The honest
    // freshness protocol therefore runs SEQUENTIAL SOLVES: each new solve
    // starts from a deterministic direction orthogonalized against every
    // pooled (exact-C-1-verified) eigenvector.  For (near-)normal operators
    // the new Krylov space then carries fresh projections of every eigenspace,
    // so hidden copies surface as certainly-inner candidates and block C-2.
    // Declaration requires one confirming solve that started orthogonal to a
    // complete pool and added nothing certainly-inner (C-4: still a detection
    // mechanism, not a completeness guarantee).  Pool admission is by
    // certified linear independence (residual after CGS2 projection >= 1e-3),
    // NOT by value proximity (C-3: no value dedup).  A dependent candidate
    // whose Rayleigh value is inconsistent with the matched pool member (an
    // extremely non-normal configuration) blocks declaration and degrades to
    // honest not_converged after two consecutive occurrences.
    struct pool_pair_t {
        T theta;
        std::vector<T> vec;
        R res_abs;
        R res_rel;
        // EIG-4 T-3 (allow_complex_pairs only): im > 0 marks a complex pair
        // occupying 2 slots; vec holds u and vec2 holds v of x = u + i v.
        // Always im == 0 / vec2 empty when the opt-in is off.
        R im = R(0);
        std::vector<T> vec2;
    };
    std::vector<pool_pair_t> pool;
    // slot count of the pool (pairs count 2).  Equal to pool.size() whenever
    // allow_complex_pairs is off (no pair is ever pooled then).
    std::size_t pool_slots = 0;
    // EIG-4 T-3 opt-in gate (D4-3 / B-27): all pair-returning branches below
    // are reachable only when this is true.
    const bool allow_pairs = options.allow_complex_pairs;
    std::size_t returned_pair_count = 0;
    // Orthonormal basis of span(pool vectors).  Eigenvectors of a nonnormal
    // operator are mutually non-orthogonal, so the independence test and the
    // pool-orthogonal starts must project onto this basis, NOT Gram-Schmidt
    // against the raw (oblique) pool vectors.
    std::vector<std::vector<T> > pool_onb;
    // Independence admission threshold: a residual-verified eigenvector with
    // exact residual res has direction error ~ res*kappa_v/gap; content
    // outside span(pool) exceeding tau_add therefore certifies a NEW
    // eigenspace dimension whenever res*kappa_v/gap << tau_add (documented
    // limitation: extremely nonnormal pairs with vector error > 1e-6 fall to
    // the dependent/consistency path and can only degrade to honest
    // not_converged, never to a lie).
    const R tau_add = vcp::tsparse_scalar::decimal_power_negative<R>(6);
    // Complex blocks whose invariant subspace lies within span(pool) up to
    // this tolerance are rediscovery artifacts (two pooled copies of a real
    // eigenvalue projecting to a 2x2 with an eps-level discriminant): they are
    // excluded from the return prefix and the D3-2 veto (their subspace is
    // already returned as verified real pairs).
    const R cover_tol = vcp::tsparse_scalar::decimal_power_negative<R>(2);
    const std::size_t k_eff = k;   // k already clamped to n
    std::size_t stable_confirms = 0;
    std::size_t dependent_block_solves = 0;
    bool start_new_solve = true;
    bool pool_complete_at_start = false;
    bool declare_from_pool = false;
    // ---- EIG-7 α/β(G-0.1 承認: additive route + per-solve cap)------------
    // α: 確認ソルブ(pool_complete_at_start)中、prefix_ready に至らない解析でも
    //    契約 C-2 の certainly 判定(共有 honest_termination の呼び出しのみ —
    //    B-43)で confirm を評価する**追加**経路。prefix_ready == true の既存
    //    経路は完全不変(G-0.1 c-2: fast-path 無改変)。最小走査(設計 §2.3):
    //    全展開到達 + probe cycle 1 回経由 + 候補全景の evidence 輸出
    //    (evidence は最終解析の final_cands から従来経路で輸出される)。
    // β: 確認ソルブ 1 本あたりのリスタート上限(cap 到達で当該ソルブを放棄し
    //    新しい pool 直交ソルブへ)。放棄→再試行ループは B-1 で必ず停止する:
    //    各リスタートの展開は apply を >= 1 消費し、mv_count >= max_mv で
    //    budget_exhausted 停止に到達するため(G-0.1 c-4 設計メモ)。
    bool alpha_pending = false;               // α の probe cycle 待ち
    std::size_t solve_start_restart = 0;      // 現在ソルブ開始時点の restart_count
    // EIG-7 (iii)(G-0.1rev 裁定・r-1): α のエンゲージは「停滞」検知による。
    // 進捗イベント :=
    //   (E1) pool 追加(pool_slots の増加を外側検知 — 多重度コピー浮上を含む)
    //   (E2) 確認ソルブ内の解析で target prefix の bound 収束ブロック数が
    //        当該チェイス区間(pool 追加で区切る)の過去最大を更新
    // stall = 最後の進捗イベントからのリスタート数。stall >= S' で α をエンゲージ。
    // 実測(r-2): diagmult チェイス(浮上 33)の全域でイベント間隔 <= 15
    // (E2 が 5〜8 リスタートごと・チェイスが構造的に α を遅延)。west 空転は
    // rc=11 以降イベントゼロ。S' = 24 = 実測最大チェイス間隔 15 の 1.6 倍。
    const std::size_t alpha_stall_engage = 24;
    std::size_t stall_last_event_rc = 0;
    std::size_t solve_max_prefix_conv = 0;
    std::size_t solve_seen_pool_slots = 0;
    const std::size_t confirm_solve_cap = 128; // β cap(G-0.1 の 8 は実測で健全多重度追跡と衝突 — 改訂承認待ち。健全最大 70 の ~2 倍)
    std::size_t confirm_restarts_total = 0;   // 診断(累計 — B-44)
    R val_scale_seen = R(1);             // running max |Ritz value| (consistency scale)

    // final-analysis bookkeeping (exported on all exits)
    struct cand_t {
        std::size_t pos, size;
        T re, im;        // for complex block: re +/- im (im > 0)
        T re2;           // second value of a (split) real pair — unused after split
        R bound;
        bool converged;
        bool is_complex;
    };
    std::vector<cand_t> final_cands;
    bool have_final_analysis = false;
    std::string analysis_fail_reason;

    // declared results (filled at converged declaration)
    bool declared_converged = false;
    bool complex_window_reject = false;
    std::size_t complex_window_count = 0;

    struct lifted_pair_t {
        std::vector<T> v;
        R res_abs;
        R res_rel;
        // EIG-4 T-3: second column (imaginary part v of x = u + i v) for a
        // complex-pair block; empty for real blocks / when opt-in is off.
        std::vector<T> v2;
    };

    // =========================================================================
    // main loop (sequential solves; per-solve phases 1..7)
    // =========================================================================
    while (true) {
        // ---- PHASE 0: (re)start a solve orthogonal to the pool ---------------
        if (start_new_solve) {
            // EIG-7 β: 直前ソルブが確認ソルブなら消費リスタートを累計に計上
            if (pool_complete_at_start)
                confirm_restarts_total += restart_count - solve_start_restart;
            solve_start_restart = restart_count;
            alpha_pending = false;
            stall_last_event_rc = restart_count;   // EIG-7 (iii): 停滞計測リセット
            solve_max_prefix_conv = 0;
            solve_seen_pool_slots = pool_slots;
            start_new_solve = false;
            V.clear();
            S.clear();
            coup.clear();
            v_next.clear();
            pending_verification = false;
            subspace_exhausted = false;
            pool_complete_at_start = (pool_slots >= k_eff);
            // deterministic start, orthogonalized against span(pool)
            std::vector<T> v0;
            if (!kd::fresh_orthogonal_direction(pool_onb, pool_onb.size(), n, seed,
                                                floor_tol, v0)) {
                // the pool spans (numerically) the whole space: nothing new can
                // be scanned — declare from the pool if it is complete
                if (pool_slots >= k_eff) { declare_from_pool = true; break; }
                result.converged = false;
                result.status = "failed";
                result.failure_reason = "no start direction available (pool-orthogonal)";
                break;
            }
            v_next = v0;
        }

        // ---- PHASE 1: expand to m_limit (budget-gated) ----------------------
        bool did_expand = false;
        while (V.size() < m_limit && !v_next.empty()) {
            if (mv_count >= max_mv) { budget_exhausted = true; break; }
            const std::size_t M = V.size();
            // absorb v_next as applied column M
            V.push_back(v_next);
            std::vector<T> w;
            apply(V[M], w);
            mv_count++;
            did_expand = true;

            std::vector<T> h;
            kd::orthogonalize_block(V, M + 1, w, h, orth_passes);

            // probe fold (verification cycle): orthogonalize against a fresh
            // deterministic direction p too, then fold p and the remainder
            // into a single new residual direction (rank-one tail preserved:
            //   A V[M] = sum h_i V_i + h_p p + beta u = sum h_i V_i + gamma v')
            std::vector<T> pdir;
            T h_p = T(0);
            bool have_probe = false;
            if (pending_verification) {
                if (kd::fresh_orthogonal_direction(V, M + 1, n, seed, floor_tol, pdir)) {
                    R c1 = vcp::tsparse_scalar::real_dot_value(pdir, w);
                    for (std::size_t i = 0; i < n; i++) w[i] -= T(c1) * pdir[i];
                    const R c2 = vcp::tsparse_scalar::real_dot_value(pdir, w);
                    for (std::size_t i = 0; i < n; i++) w[i] -= T(c2) * pdir[i];
                    h_p = T(c1) + T(c2);
                    have_probe = true;
                }
            }

            // grow S by one row/column; install column M
            S.push_back(std::vector<T>(M, T(0)));
            for (std::size_t i = 0; i <= M; i++) {
                if (S[i].size() < M + 1) S[i].resize(M + 1, T(0));
                S[i][M] = h[i];
            }
            // previous coupling row becomes S row M for old columns
            for (std::size_t j = 0; j < coup.size() && j < M; j++) S[M][j] = coup[j];

            const R beta = vcp::tsparse_scalar::real_norm_value(w);
            std::vector<T> q;                 // new residual direction
            T gamma = T(0);
            if (have_probe) {
                // q_vec = h_p * p + w  (w = beta * u already)
                q.assign(n, T(0));
                for (std::size_t i = 0; i < n; i++) q[i] = h_p * pdir[i] + w[i];
                const R gq = vcp::tsparse_scalar::real_norm_value(q);
                if (gq > floor_tol) {
                    for (std::size_t i = 0; i < n; i++) q[i] /= T(gq);
                    gamma = T(gq);
                } else if (beta > floor_tol) {
                    q = w;
                    for (std::size_t i = 0; i < n; i++) q[i] /= T(beta);
                    gamma = T(beta);
                } else {
                    // both negligible: keep the probe as the next direction,
                    // coupling recorded honestly as the computed fold norm
                    q = pdir;
                    gamma = T(gq);
                }
            } else if (beta > floor_tol) {
                q = w;
                for (std::size_t i = 0; i < n; i++) q[i] /= T(beta);
                gamma = T(beta);
            } else {
                // near-breakdown: continue with a fresh orthogonal direction;
                // the tiny coupling is recorded honestly (never forced to 0)
                std::vector<T> fd;
                if (kd::fresh_orthogonal_direction(V, M + 1, n, seed, floor_tol, fd)) {
                    // fold the (tiny) remainder into the fresh direction so the
                    // rank-one tail stays exact: gamma q = w + 0*fd is wrong;
                    // instead q = normalize(w + floor-scaled fd)?  Honest and
                    // simple: q = normalized(w) if certifiable else fd with
                    // coupling = ||w|| (residual bounded by the recorded norm).
                    if (beta > R(0)) {
                        q = w;
                        for (std::size_t i = 0; i < n; i++) q[i] /= T(beta);
                        gamma = T(beta);
                    } else {
                        q = fd;
                        gamma = T(beta);   // == 0 exactly (or interval enclosure)
                    }
                } else {
                    // no direction left: invariant subspace spans everything
                    subspace_exhausted = true;
                    q.clear();
                    gamma = T(beta);
                }
            }
            // coupling row for the new factorization: gamma * e_M
            coup.assign(M + 1, T(0));
            coup[M] = gamma;
            v_next = q;   // may be empty (subspace exhausted)
            if (subspace_exhausted) break;
        }
        (void)did_expand;

        const std::size_t M = V.size();
        if (M == 0) {
            result.converged = false;
            result.status = budget_exhausted ? "max_iter_exhausted" : "failed";
            result.failure_reason = budget_exhausted
                ? "matrix-vector product budget exhausted before any expansion"
                : "empty Krylov basis";
            break;
        }

        // ---- PHASE 2: projected analysis  S = U T U^T -----------------------
        std::vector<std::vector<T> > Tm(M, std::vector<T>(M, T(0)));
        for (std::size_t i = 0; i < M; i++)
            for (std::size_t j = 0; j < M; j++) Tm[i][j] = S[i][j];

        std::vector<std::vector<T> > Q0;
        bool analysis_ok =
            vcp::tsparse_dense_schur::dense_schur_detail::hessenberg_reduce_(Tm, Q0);
        std::vector<std::vector<T> > U;
        if (analysis_ok) {
            vcp::tsparse_real_schur::real_schur_result<T> sr =
                vcp::tsparse_real_schur::real_schur_decompose<T>(Tm, true);
            if (!sr.success) {
                analysis_ok = false;
                analysis_fail_reason = "projected real Schur core: " + sr.failure_reason;
            } else {
                // U = Q0 * Z
                U.assign(M, std::vector<T>(M, T(0)));
                for (std::size_t i = 0; i < M; i++)
                    for (std::size_t j = 0; j < M; j++) {
                        T s(0);
                        for (std::size_t t = 0; t < M; t++)
                            s += Q0[i][t] * sr.schur_vectors[t][j];
                        U[i][j] = s;
                    }
                Tm.swap(sr.schur_form);
            }
        } else {
            analysis_fail_reason = "projected Hessenberg reduction not certifiable";
        }
        if (!analysis_ok) {
            result.converged = false;
            result.status = "failed";
            result.failure_reason = analysis_fail_reason;
            break;
        }

        // ---- PHASE 3: blocks, values, target order --------------------------
        std::vector<std::size_t> bpos, bsize;
        kd::enumerate_blocks(Tm, bpos, bsize);
        // certified re-standardization (post-swap real pairs are split)
        {
            bool changed = true;
            int guard = 0;
            while (changed && guard++ < 64) {
                changed = false;
                kd::enumerate_blocks(Tm, bpos, bsize);
                for (std::size_t bi = 0; bi < bpos.size(); bi++) {
                    if (bsize[bi] != 2) continue;
                    T re, im, r1, r2;
                    const int cls = kd::eig_2x2(Tm, bpos[bi], re, im, r1, r2);
                    if (cls == 0) {
                        if (!kd::split_real_2x2(Tm, U, bpos[bi])) {
                            analysis_ok = false;
                            analysis_fail_reason =
                                "2x2 real-pair standardization not certifiable";
                        }
                        changed = analysis_ok;
                        break;
                    }
                    if (cls < 0) {
                        analysis_ok = false;
                        analysis_fail_reason =
                            "2x2 discriminant sign not certifiable (0-straddle)";
                        break;
                    }
                }
                if (!analysis_ok) break;
            }
        }
        if (!analysis_ok) {
            result.converged = false;
            result.status = "failed";
            result.failure_reason = analysis_fail_reason;
            break;
        }
        kd::enumerate_blocks(Tm, bpos, bsize);
        const std::size_t nblocks = bpos.size();

        // block values + target keys
        std::vector<T> bre(nblocks, T(0)), bim(nblocks, T(0));
        std::vector<R> bkey(nblocks, R(0));
        bool value_fail = false;
        for (std::size_t bi = 0; bi < nblocks; bi++) {
            if (bsize[bi] == 1) {
                bre[bi] = Tm[bpos[bi]][bpos[bi]];
                bim[bi] = T(0);
            } else {
                T re, im, r1, r2;
                const int cls = kd::eig_2x2(Tm, bpos[bi], re, im, r1, r2);
                if (cls != 1) { value_fail = true; break; }
                bre[bi] = re;
                bim[bi] = im;
            }
            bkey[bi] = vcp::tsparse_eigen_selection::target_distance(
                std::complex<R>(vcp::tsparse_scalar::real_part(bre[bi]),
                                vcp::tsparse_scalar::real_part(bim[bi])),
                target, shift_val);
        }
        if (value_fail) {
            result.converged = false;
            result.status = "failed";
            result.failure_reason = "block eigenvalues not certifiable";
            break;
        }

        // target order of blocks: hand-rolled insertion sort with certainly
        // comparisons (partial-order tolerant — P6; exact for double)
        std::vector<std::size_t> order(nblocks);
        for (std::size_t i = 0; i < nblocks; i++) order[i] = i;
        for (std::size_t i = 1; i < nblocks; i++) {
            const std::size_t oi = order[i];
            std::size_t j = i;
            while (j > 0) {
                const bool before = prefer_large ? (bkey[oi] > bkey[order[j - 1]])
                                                 : (bkey[oi] < bkey[order[j - 1]]);
                if (!before) break;
                order[j] = order[j - 1];
                j--;
            }
            order[j] = oi;
        }

        // ---- PHASE 4: reorder selected blocks to the front (G-2.1 method A) --
        // Selection: target-order prefix filling k return slots (2x2 blocks
        // keep-together, counted as 2 slots, boundary widened conservatively),
        // plus padding, plus every bound-converged block (B-15), capped at
        // M-1 columns (at least one expansion slot must remain — the mv
        // budget then guarantees overall termination, G-2.1 c-2).
        //
        // Bounds need the reordered T; but selection needs bounds… resolve by
        // two passes: (pass A) reorder the target-prefix; (pass B) compute
        // bounds on the reordered form; keep set final = prefix ∪ converged.
        // Movement: selection-sort by target rank with adjacent block swaps.
        std::size_t rejections_this = 0;
        std::size_t keep_more_cols = 0;   // G-2.1 c-1: window extension demands
        {
            // current block id order by position
            std::vector<std::size_t> ids(order.size());
            // ids[slot] = block id currently at slot (position order)
            for (std::size_t i = 0; i < nblocks; i++) ids[i] = i;
            // positions recomputed from sizes as we go
            std::vector<std::size_t> cursize(nblocks);
            for (std::size_t i = 0; i < nblocks; i++) cursize[i] = bsize[i];

            // how many leading slots we want ordered: enough blocks to cover
            // k return slots + padding
            std::size_t want_cols = k + ((k / 2 + 1 > 2) ? k / 2 + 1 : 2);
            if (want_cols > M - 1) want_cols = M - 1;

            std::size_t placed_slots = 0;   // block slots already固定 at front
            std::size_t placed_cols = 0;
            for (std::size_t r = 0; r < order.size() && placed_cols < want_cols; r++) {
                const std::size_t id = order[r];
                // find current slot of block id
                std::size_t slot = placed_slots;
                bool found = false;
                for (std::size_t s2 = 0; s2 < ids.size(); s2++)
                    if (ids[s2] == id) { slot = s2; found = true; break; }
                if (!found || slot < placed_slots) continue;   // already placed
                // bubble the block from `slot` down to `placed_slots`
                bool stuck = false;
                while (slot > placed_slots) {
                    // positions of the two adjacent blocks
                    std::size_t p = 0;
                    for (std::size_t s2 = 0; s2 + 1 <= slot - 1; s2++) p += cursize[ids[s2]];
                    const std::size_t s_up = cursize[ids[slot - 1]];
                    const std::size_t s_dn = cursize[ids[slot]];
                    if (kd::swap_adjacent_blocks(Tm, U, p, s_up, s_dn)) {
                        const std::size_t tid = ids[slot - 1];
                        ids[slot - 1] = ids[slot];
                        ids[slot] = tid;
                        slot--;
                    } else {
                        rejections_this++;
                        stuck = true;
                        break;
                    }
                }
                if (stuck) {
                    // keep-more (G-2.1 c-1): the candidate stays where it is;
                    // the keep window is extended up to its current end
                    std::size_t endpos = 0;
                    for (std::size_t s2 = 0; s2 <= slot; s2++) endpos += cursize[ids[s2]];
                    if (endpos > want_cols) want_cols = (endpos <= M - 1) ? endpos : M - 1;
                    if (endpos > keep_more_cols) keep_more_cols = endpos;
                    // (an m_limit-1 cap-forced narrowing can only delay
                    //  convergence, never fabricate it — C-1/C-2 still gate)
                } else {
                    placed_slots++;
                    placed_cols += cursize[id];
                }
            }
            diag.contraction_swap_rejections += rejections_this;
        }
        kd::enumerate_blocks(Tm, bpos, bsize);
        // refresh values/keys after reordering
        {
            const std::size_t nb2 = bpos.size();
            bre.assign(nb2, T(0)); bim.assign(nb2, T(0)); bkey.assign(nb2, R(0));
            bool vf = false;
            for (std::size_t bi = 0; bi < nb2; bi++) {
                if (bsize[bi] == 1) { bre[bi] = Tm[bpos[bi]][bpos[bi]]; }
                else {
                    T re, im, r1, r2;
                    if (kd::eig_2x2(Tm, bpos[bi], re, im, r1, r2) != 1) { vf = true; break; }
                    bre[bi] = re; bim[bi] = im;
                }
                bkey[bi] = vcp::tsparse_eigen_selection::target_distance(
                    std::complex<R>(vcp::tsparse_scalar::real_part(bre[bi]),
                                    vcp::tsparse_scalar::real_part(bim[bi])),
                    target, shift_val);
            }
            if (vf) {
                result.converged = false;
                result.status = "failed";
                result.failure_reason = "block eigenvalues not certifiable (post-reorder)";
                break;
            }
        }

        // transformed coupling row  c_T = c * U
        std::vector<T> cT(M, T(0));
        for (std::size_t j = 0; j < M; j++) {
            T s(0);
            for (std::size_t i = 0; i < M && i < coup.size(); i++) s += coup[i] * U[i][j];
            cT[j] = s;
        }

        // ---- PHASE 5: candidate bounds via quasi-triangular back-substitution
        // Invariant columns Y for block bi: rows 0..pos+size-1; Y[pos..] = I.
        // Upward solve with pivot blocks of Tm; a pivot that is not certifiably
        // nonsingular contributes X = 0 (deterministic; final honesty is
        // guaranteed by the exact C-1 verification, not by these bounds).
        const std::size_t nb = bpos.size();
        std::vector<R> bound(nb, R(0));
        std::vector<bool> bconv(nb, false);
        std::vector<std::vector<std::vector<T> > > yof(nb);   // per block: Y (rows x size)
        R tnorm(0);
        for (std::size_t i = 0; i < M; i++)
            for (std::size_t j = 0; j < M; j++) {
                const R a = abs_value(Tm[i][j]);
                if (a > tnorm) tnorm = a;
            }
        for (std::size_t bi = 0; bi < nb; bi++) {
            const std::size_t p = bpos[bi];
            const std::size_t sz = bsize[bi];
            std::vector<std::vector<T> > Y(p + sz, std::vector<T>(sz, T(0)));
            for (std::size_t j = 0; j < sz; j++) Y[p + j][j] = T(1);
            // pivot blocks above
            std::vector<std::size_t> ppos, psize;
            for (std::size_t bj = 0; bj < bi; bj++) { ppos.push_back(bpos[bj]); psize.push_back(bsize[bj]); }
            // B = the block's own matrix (sz x sz)
            std::vector<std::vector<T> > Bblk(sz, std::vector<T>(sz, T(0)));
            for (std::size_t i = 0; i < sz; i++)
                for (std::size_t j = 0; j < sz; j++) Bblk[i][j] = Tm[p + i][p + j];
            for (std::size_t bj = ppos.size(); bj-- > 0;) {
                const std::size_t q = ppos[bj];
                const std::size_t qs = psize[bj];
                // RHS = - T[q..q+qs, q+qs..p+sz-1] * Y[q+qs..p+sz-1]
                std::vector<std::vector<T> > Rhs(qs, std::vector<T>(sz, T(0)));
                for (std::size_t i = 0; i < qs; i++)
                    for (std::size_t c2 = q + qs; c2 < p + sz; c2++)
                        for (std::size_t j = 0; j < sz; j++)
                            Rhs[i][j] -= Tm[q + i][c2] * Y[c2][j];
                // solve P X - X B = RHS  (X qs x sz)
                const std::size_t nun = qs * sz;
                std::vector<std::vector<T> > Msm(nun, std::vector<T>(nun, T(0)));
                std::vector<T> rv(nun, T(0));
                for (std::size_t i = 0; i < qs; i++)
                    for (std::size_t j = 0; j < sz; j++) {
                        const std::size_t row = i * sz + j;
                        rv[row] = Rhs[i][j];
                        for (std::size_t t = 0; t < qs; t++)
                            Msm[row][t * sz + j] += Tm[q + i][q + t];
                        for (std::size_t t = 0; t < sz; t++)
                            Msm[row][i * sz + t] -= Bblk[t][j];
                    }
                std::vector<T> xv;
                if (kd::solve_small(Msm, rv, xv)) {
                    // magnitude guard: reject explosive solutions (near-singular
                    // pivots for multiplicities) — deterministic X = 0 instead
                    R xmax(0);
                    for (std::size_t t = 0; t < nun; t++) {
                        const R a = abs_value(xv[t]);
                        if (a > xmax) xmax = a;
                    }
                    R ymax(0);
                    for (std::size_t r2 = q + qs; r2 < p + sz; r2++)
                        for (std::size_t j = 0; j < sz; j++) {
                            const R a = abs_value(Y[r2][j]);
                            if (a > ymax) ymax = a;
                        }
                    const R lim = (ymax > R(1) ? ymax : R(1)) /
                                  (sqrt_value(eps) > R(0) ? sqrt_value(eps) : R(1));
                    if (!(xmax > lim)) {
                        for (std::size_t i = 0; i < qs; i++)
                            for (std::size_t j = 0; j < sz; j++)
                                Y[q + i][j] = xv[i * sz + j];
                    }
                    // else: leave zero (deterministic complement choice)
                }
                // else: pivot not certifiably nonsingular -> X = 0
            }
            // bound = ||cT * Y||_F / ||Y||_F  (relative residual of the pair)
            R num(0), den(0);
            for (std::size_t j = 0; j < sz; j++) {
                T s(0);
                for (std::size_t r2 = 0; r2 < p + sz; r2++) s += cT[r2] * Y[r2][j];
                const R a = abs_value(s);
                num += a * a;
            }
            for (std::size_t r2 = 0; r2 < p + sz; r2++)
                for (std::size_t j = 0; j < sz; j++) {
                    const R a = abs_value(Y[r2][j]);
                    den += a * a;
                }
            const R nn = sqrt_value(num);
            const R dd = sqrt_value(den);
            bound[bi] = (dd > R(0)) ? nn / dd : nn;
            const R scale = R(1) + vcp::tsparse_eigen_selection::target_distance(
                std::complex<R>(vcp::tsparse_scalar::real_part(bre[bi]),
                                vcp::tsparse_scalar::real_part(bim[bi])),
                vcp::eig_target::largest_magnitude, R(0));
            bconv[bi] = (bound[bi] <= tol * scale);
            yof[bi].swap(Y);
        }

        // record final analysis for evidence export
        final_cands.clear();
        for (std::size_t bi = 0; bi < nb; bi++) {
            cand_t c;
            c.pos = bpos[bi]; c.size = bsize[bi];
            c.re = bre[bi]; c.im = bim[bi]; c.re2 = T(0);
            c.bound = bound[bi];
            c.converged = bconv[bi];
            c.is_complex = (bsize[bi] == 2);
            final_cands.push_back(c);
        }
        have_final_analysis = true;

        // residual history (best bound this analysis)
        if (options.compute_residual_history && nb > 0) {
            R best = bound[0];
            R besttheta = abs_value(vcp::tsparse_scalar::real_part(bre[0]));
            for (std::size_t bi = 1; bi < nb; bi++)
                if (bound[bi] < best) {
                    best = bound[bi];
                    besttheta = abs_value(vcp::tsparse_scalar::real_part(bre[bi]));
                }
            result.residual_history_absolute.push_back(best);
            result.residual_history_relative.push_back(best / (R(1) + besttheta));
        }

        // ---- PHASE 6: termination logic -------------------------------------
        // coverage of complex blocks by the pool span (rediscovery artifacts)
        std::vector<bool> covered(nb, false);
        if (!pool_onb.empty()) {
            for (std::size_t bi = 0; bi < nb; bi++) {
                if (bsize[bi] != 2) continue;
                const std::vector<std::vector<T> >& Y = yof[bi];
                bool all_cov = true;
                for (std::size_t col = 0; col < 2 && all_cov; col++) {
                    std::vector<T> z(M, T(0));
                    for (std::size_t r2 = 0; r2 < M; r2++) {
                        T s(0);
                        for (std::size_t t = 0; t < Y.size(); t++) s += U[r2][t] * Y[t][col];
                        z[r2] = s;
                    }
                    std::vector<T> v(n, T(0));
                    for (std::size_t r2 = 0; r2 < M; r2++) {
                        const T zr = z[r2];
                        for (std::size_t ii = 0; ii < n; ii++) v[ii] += V[r2][ii] * zr;
                    }
                    const R vn2 = vcp::tsparse_scalar::real_norm_value(v);
                    if (!(vn2 > R(0))) { all_cov = false; break; }
                    for (std::size_t ii = 0; ii < n; ii++) v[ii] /= T(vn2);
                    for (int pass = 0; pass < 2; pass++)
                        for (std::size_t j2 = 0; j2 < pool_onb.size(); j2++) {
                            const R c2v = vcp::tsparse_scalar::real_dot_value(pool_onb[j2], v);
                            for (std::size_t ii = 0; ii < n; ii++)
                                v[ii] -= T(c2v) * pool_onb[j2][ii];
                        }
                    if (!(vcp::tsparse_scalar::real_norm_value(v) <= cover_tol))
                        all_cov = false;
                }
                if (all_cov) {
                    covered[bi] = true;
                    final_cands[bi].converged = true;   // accounted-for in C-2 terms
                }
            }
        }

        // target-order over blocks of the FINAL (reordered) analysis
        std::vector<std::size_t> ord2(nb);
        for (std::size_t i = 0; i < nb; i++) ord2[i] = i;
        for (std::size_t i = 1; i < nb; i++) {
            const std::size_t oi = ord2[i];
            std::size_t j = i;
            while (j > 0) {
                const bool before = prefer_large ? (bkey[oi] > bkey[ord2[j - 1]])
                                                 : (bkey[oi] < bkey[ord2[j - 1]]);
                if (!before) break;
                ord2[j] = ord2[j - 1];
                j--;
            }
            ord2[j] = oi;
        }
        // return prefix: first non-covered blocks filling k slots
        std::vector<std::size_t> prefix_blocks;
        std::size_t slots = 0;
        for (std::size_t i = 0; i < nb && slots < k; i++) {
            if (covered[ord2[i]]) continue;   // pool-covered rediscovery artifact
            prefix_blocks.push_back(ord2[i]);
            slots += bsize[ord2[i]];
        }

        // D3-2: a converged complex pair inside the return window means the
        // real-only result can never honestly return k converged real pairs.
        // Reference window: the pool prefix once the pool is complete (later
        // scanning solves legitimately meet outer complex pairs in their own
        // complement prefix — those must not fire the rejection).
        {
            R pool_worst_key = R(0);
            bool have_pool_ref = false;
            if (allow_pairs) {
                // EIG-4 T-3: slot-aware worst key over pooled (re, im) values
                if (pool_slots >= k_eff) {
                    std::vector<R> pre, pim;
                    for (std::size_t i = 0; i < pool.size(); i++) {
                        pre.push_back(vcp::tsparse_scalar::real_part(pool[i].theta));
                        pim.push_back(pool[i].im);
                    }
                    have_pool_ref = kd::pool_worst_key_pairs<R>(
                        pre, pim, k_eff, target, shift_val, pool_worst_key);
                }
            } else if (pool.size() >= k_eff) {
                std::vector<R> pv2;
                for (std::size_t i = 0; i < pool.size(); i++)
                    pv2.push_back(vcp::tsparse_scalar::real_part(pool[i].theta));
                const std::vector<std::size_t> po =
                    vcp::tsparse_eigen_selection::select_eigen_indices_from_real(
                        pv2, pv2.size(), target, shift_val);
                if (po.size() >= k_eff) {
                    pool_worst_key = vcp::tsparse_eigen_selection::target_distance(
                        std::complex<R>(pv2[po[k_eff - 1]], R(0)), target, shift_val);
                    have_pool_ref = true;
                }
            }
            bool cx_conv_inner = false;
            std::size_t cx_count = 0;
            for (std::size_t i = 0; i < prefix_blocks.size(); i++) {
                const std::size_t bi = prefix_blocks[i];
                if (bsize[bi] != 2 || covered[bi]) continue;
                bool counts = true;
                if (have_pool_ref) {
                    const R key = bkey[bi];
                    counts = prefer_large ? (key > pool_worst_key)
                                          : (key < pool_worst_key);
                }
                if (counts) {
                    cx_count++;
                    if (bconv[bi]) cx_conv_inner = true;
                }
            }
            // EIG-4 T-3 (B-27): with the opt-in ON the D3-2 veto is lifted --
            // converged pairs inside the window become returnable below.
            if (cx_conv_inner && !allow_pairs) {
                complex_window_reject = true;
                complex_window_count = cx_count;
                result.converged = false;
                break;
            }
        }

        // prefix all real and bound-converged?  attempt declaration
        // (EIG-4 T-3: with the opt-in ON, bound-converged 2x2 pair blocks are
        // also declaration-eligible)
        bool prefix_ready = (slots >= k);
        for (std::size_t i = 0; prefix_ready && i < prefix_blocks.size(); i++) {
            const std::size_t bi = prefix_blocks[i];
            const bool size_ok = (bsize[bi] == 1)
                              || (allow_pairs && bsize[bi] == 2);
            if (!size_ok || !bconv[bi]) prefix_ready = false;
        }

        if (prefix_ready) {
            // exact C-1 verification of the k prefix pairs (budget-gated;
            // pair blocks need one apply per column = 2)
            std::size_t prefix_cols = 0;
            for (std::size_t i = 0; i < prefix_blocks.size(); i++)
                prefix_cols += bsize[prefix_blocks[i]];
            if (mv_count + prefix_cols > max_mv) {
                budget_exhausted = true;
                result.converged = false;
                break;
            }
            bool all_pass = true;
            std::vector<lifted_pair_t> lifted(prefix_blocks.size());
            for (std::size_t i = 0; i < prefix_blocks.size(); i++) {
                const std::size_t bi = prefix_blocks[i];
                const std::vector<std::vector<T> >& Y = yof[bi];
                if (bsize[bi] == 1) {
                // T-basis -> V-basis: z = U * y (Y lives in the reordered
                // Schur coordinates; the state basis is V = (basis) with
                // S = U T U^T), then v = V z.
                std::vector<T> z(M, T(0));
                for (std::size_t r2 = 0; r2 < M; r2++) {
                    T s(0);
                    for (std::size_t t = 0; t < Y.size(); t++) s += U[r2][t] * Y[t][0];
                    z[r2] = s;
                }
                std::vector<T> v(n, T(0));
                for (std::size_t r2 = 0; r2 < M; r2++) {
                    const T zr = z[r2];
                    for (std::size_t ii = 0; ii < n; ii++) v[ii] += V[r2][ii] * zr;
                }
                const R vn = vcp::tsparse_scalar::real_norm_value(v);
                if (!(vn > R(0))) { all_pass = false; bconv[bi] = false; break; }
                for (std::size_t ii = 0; ii < n; ii++) v[ii] /= T(vn);
                std::vector<T> Av;
                apply(v, Av);
                mv_count++;
                const T theta = bre[bi];
                std::vector<T> rr(n);
                for (std::size_t ii = 0; ii < n; ii++) rr[ii] = Av[ii] - theta * v[ii];
                const R ra = vcp::tsparse_scalar::real_norm_value(rr);
                const R rrel = ra / (R(1) + abs_value(vcp::tsparse_scalar::real_part(theta)));
                lifted[i].v.swap(v);
                lifted[i].res_abs = ra;
                lifted[i].res_rel = rrel;
                const bool acc = (ra <= tol) || (rrel <= tol);
                if (!acc) {
                    all_pass = false;
                    bconv[bi] = false;
                    final_cands[bi].converged = false;
                    break;
                }
                } else {
                // EIG-4 T-3 (allow_pairs only; prefix_ready forbids 2x2
                // otherwise): lift BOTH columns of Y -- Tm Y = Y Bblk, so
                // W = V (U Y) spans the invariant plane with A W ~= W Bblk.
                // Transform Bblk to the standard rotation form via the real
                // 2x2 eigenbasis E: Bblk E = E [[p, q], [-q, p]] (real
                // arithmetic only, D4-3), then verify the exact pair residual
                // || A [u v] - [u v] [[p, q], [-q, p]] ||_F.
                const std::size_t p0 = bpos[bi];
                std::vector<std::vector<T> > w(2, std::vector<T>(n, T(0)));
                for (std::size_t j = 0; j < 2; j++) {
                    std::vector<T> z(M, T(0));
                    for (std::size_t r2 = 0; r2 < M; r2++) {
                        T s(0);
                        for (std::size_t t = 0; t < Y.size(); t++) s += U[r2][t] * Y[t][j];
                        z[r2] = s;
                    }
                    for (std::size_t r2 = 0; r2 < M; r2++) {
                        const T zr = z[r2];
                        for (std::size_t ii = 0; ii < n; ii++) w[j][ii] += V[r2][ii] * zr;
                    }
                }
                const T s11 = Tm[p0][p0],     s12 = Tm[p0][p0 + 1];
                const T s21 = Tm[p0 + 1][p0], s22 = Tm[p0 + 1][p0 + 1];
                const T pr = bre[bi];
                const T qi = bim[bi];   // > 0 (certified pair)
                // real 2x2 eigenbasis of Bblk for lambda = p + i q
                T e_r0, e_r1, e_i0, e_i1;
                if (!(abs_value(s12) < abs_value(s21))) {
                    e_r0 = s12; e_r1 = pr - s11; e_i0 = T(0); e_i1 = qi;
                } else {
                    e_r0 = pr - s22; e_r1 = s21; e_i0 = qi; e_i1 = T(0);
                }
                std::vector<T> u(n), v2c(n);
                for (std::size_t ii = 0; ii < n; ii++) {
                    u[ii]   = w[0][ii] * e_r0 + w[1][ii] * e_r1;
                    v2c[ii] = w[0][ii] * e_i0 + w[1][ii] * e_i1;
                }
                const R un = vcp::tsparse_scalar::real_norm_value(u);
                const R vn2 = vcp::tsparse_scalar::real_norm_value(v2c);
                const R pn = sqrt_value(un * un + vn2 * vn2);
                if (!(pn > R(0))) { all_pass = false; bconv[bi] = false; break; }
                for (std::size_t ii = 0; ii < n; ii++) { u[ii] /= T(pn); v2c[ii] /= T(pn); }
                std::vector<T> Au, Av2;
                apply(u, Au);
                apply(v2c, Av2);
                mv_count += 2;
                R rs(0);
                for (std::size_t ii = 0; ii < n; ii++) {
                    const T d1 = Au[ii]  - (pr * u[ii]   - qi * v2c[ii]);
                    const T d2 = Av2[ii] - (qi * u[ii]   + pr * v2c[ii]);
                    const R a1 = abs_value(d1);
                    const R a2 = abs_value(d2);
                    rs += a1 * a1 + a2 * a2;
                }
                const R ra = sqrt_value(rs);
                const R mag = sqrt_value(
                    vcp::tsparse_scalar::real_part(pr) * vcp::tsparse_scalar::real_part(pr)
                  + vcp::tsparse_scalar::real_part(qi) * vcp::tsparse_scalar::real_part(qi));
                const R rrel = ra / (R(1) + mag);
                lifted[i].v.swap(u);
                lifted[i].v2.swap(v2c);
                lifted[i].res_abs = ra;
                lifted[i].res_rel = rrel;
                const bool acc = (ra <= tol) || (rrel <= tol);
                if (!acc) {
                    all_pass = false;
                    bconv[bi] = false;
                    final_cands[bi].converged = false;
                    break;
                }
                }
            }
            if (all_pass) {
                // C-2 on the remaining candidates (shared complex check)
                std::vector<R> ar, ai;
                std::vector<bool> ac;
                std::vector<R> locked_vals;
                for (std::size_t i = 0; i < prefix_blocks.size(); i++)
                    locked_vals.push_back(
                        vcp::tsparse_scalar::real_part(bre[prefix_blocks[i]]));
                for (std::size_t bi = 0; bi < nb; bi++) {
                    bool in_prefix = false;
                    for (std::size_t i = 0; i < prefix_blocks.size(); i++)
                        if (prefix_blocks[i] == bi) { in_prefix = true; break; }
                    if (in_prefix) continue;
                    const bool acct = bconv[bi] || covered[bi];
                    ar.push_back(vcp::tsparse_scalar::real_part(bre[bi]));
                    ai.push_back(vcp::tsparse_scalar::real_part(bim[bi]));
                    ac.push_back(acct);
                    if (bsize[bi] == 2) {
                        ar.push_back(vcp::tsparse_scalar::real_part(bre[bi]));
                        ai.push_back(-vcp::tsparse_scalar::real_part(bim[bi]));
                        ac.push_back(acct);
                    }
                }
                bool c2_ok;
                if (allow_pairs) {
                    // EIG-4 T-3: pair-aware locked worst key (shared additive
                    // check; magnitude targets weigh a pair by hypot(re, im))
                    std::vector<R> locked_im2;
                    for (std::size_t i = 0; i < prefix_blocks.size(); i++) {
                        const std::size_t bi = prefix_blocks[i];
                        locked_im2.push_back((bsize[bi] == 2)
                            ? vcp::tsparse_scalar::real_part(bim[bi]) : R(0));
                    }
                    c2_ok = vcp::tsparse::honest_termination_check_complex_pairs_<R>(
                        ar, ai, ac, locked_vals, locked_im2, k, target, shift_val);
                } else {
                    c2_ok = vcp::tsparse::honest_termination_check_complex_<R>(
                        ar, ai, ac, locked_vals, k, target, shift_val);
                }
                if (c2_ok) {
                    const bool structurally_complete = (M >= n) || subspace_exhausted;
                    if (structurally_complete && pool.empty()) {
                        // whole space analyzed in a single factorization:
                        // the prefix is complete by construction — declare
                        declared_converged = true;
                        diag.c2_evidence.fresh = true;
                        for (std::size_t i = 0; i < prefix_blocks.size(); i++) {
                            const std::size_t bi = prefix_blocks[i];
                            if (bsize[bi] == 1) {
                            result.eigenvalues.push_back(bre[bi]);
                            result.eigenvectors.push_back(lifted[i].v);
                            result.residuals_absolute.push_back(lifted[i].res_abs);
                            result.residuals_relative.push_back(lifted[i].res_rel);
                            if (allow_pairs) result.eigenvalues_imag.push_back(T(0));
                            } else {
                            // EIG-4 T-3: adjacent (re, re) / (+im, -im) / (u, v)
                            result.eigenvalues.push_back(bre[bi]);
                            result.eigenvalues.push_back(bre[bi]);
                            result.eigenvalues_imag.push_back(bim[bi]);
                            result.eigenvalues_imag.push_back(-bim[bi]);
                            result.eigenvectors.push_back(lifted[i].v);
                            result.eigenvectors.push_back(lifted[i].v2);
                            result.residuals_absolute.push_back(lifted[i].res_abs);
                            result.residuals_absolute.push_back(lifted[i].res_abs);
                            result.residuals_relative.push_back(lifted[i].res_rel);
                            result.residuals_relative.push_back(lifted[i].res_rel);
                            returned_pair_count++;
                            }
                        }
                        break;
                    }
                    if (!pending_verification && !structurally_complete) {
                        // one cheap in-solve probe cycle before ending the solve
                        pending_verification = true;
                    } else {
                        // ---- SOLVE VERIFIED: pool update + confirm protocol --
                        bool added_inner = false;
                        bool dependent_inconsistent = false;
                        // current pool prefix worst key (if complete)
                        R pool_worst_key = R(0);
                        bool have_pool_ref = false;
                        if (allow_pairs) {
                            // EIG-4 T-3: slot-aware worst key over (re, im)
                            if (pool_slots >= k_eff) {
                                std::vector<R> pre2, pim2;
                                for (std::size_t i2 = 0; i2 < pool.size(); i2++) {
                                    pre2.push_back(vcp::tsparse_scalar::real_part(pool[i2].theta));
                                    pim2.push_back(pool[i2].im);
                                }
                                have_pool_ref = kd::pool_worst_key_pairs<R>(
                                    pre2, pim2, k_eff, target, shift_val, pool_worst_key);
                            }
                        } else {
                            if (pool.size() >= k_eff) {
                                std::vector<R> pv2;
                                for (std::size_t i2 = 0; i2 < pool.size(); i2++)
                                    pv2.push_back(vcp::tsparse_scalar::real_part(pool[i2].theta));
                                const std::vector<std::size_t> po =
                                    vcp::tsparse_eigen_selection::select_eigen_indices_from_real(
                                        pv2, pv2.size(), target, shift_val);
                                if (po.size() >= k_eff) {
                                    pool_worst_key = vcp::tsparse_eigen_selection::target_distance(
                                        std::complex<R>(pv2[po[k_eff - 1]], R(0)), target, shift_val);
                                    have_pool_ref = true;
                                }
                            }
                        }
                        for (std::size_t i = 0; i < prefix_blocks.size(); i++) {
                            const std::size_t bi = prefix_blocks[i];
                            const bool is_pair = (bsize[bi] == 2);   // allow_pairs only
                            const R th_r = vcp::tsparse_scalar::real_part(bre[bi]);
                            const R th_i = is_pair ? vcp::tsparse_scalar::real_part(bim[bi]) : R(0);
                            {
                                R a2 = th_r; if (a2 < R(0)) a2 = -a2;
                                if (is_pair) {
                                    const R m2 = sqrt_value(th_r * th_r + th_i * th_i);
                                    if (m2 > a2) a2 = m2;
                                }
                                if (a2 > val_scale_seen) val_scale_seen = a2;
                            }
                            bool needed = !have_pool_ref;
                            if (have_pool_ref) {
                                const R key = vcp::tsparse_eigen_selection::target_distance(
                                    std::complex<R>(th_r, th_i), target, shift_val);
                                needed = prefer_large ? (key > pool_worst_key)
                                                      : (key < pool_worst_key);
                            }
                            if (!needed) continue;
                            // certified independence: residual after CGS2
                            // projection onto span(pool) (orthonormal basis)
                            std::vector<T> w = lifted[i].v;
                            for (int pass = 0; pass < 2; pass++)
                                for (std::size_t j2 = 0; j2 < pool_onb.size(); j2++) {
                                    const R c2v = vcp::tsparse_scalar::real_dot_value(pool_onb[j2], w);
                                    for (std::size_t ii = 0; ii < n; ii++)
                                        w[ii] -= T(c2v) * pool_onb[j2][ii];
                                }
                            const R rind = vcp::tsparse_scalar::real_norm_value(w);
                            // pair: the second column may carry the new content
                            // even when the first is covered (2D subspace)
                            std::vector<T> w2;
                            R rind2 = R(0);
                            if (is_pair) {
                                w2 = lifted[i].v2;
                                for (int pass = 0; pass < 2; pass++) {
                                    for (std::size_t j2 = 0; j2 < pool_onb.size(); j2++) {
                                        const R c2v = vcp::tsparse_scalar::real_dot_value(pool_onb[j2], w2);
                                        for (std::size_t ii = 0; ii < n; ii++)
                                            w2[ii] -= T(c2v) * pool_onb[j2][ii];
                                    }
                                    if (rind > R(0)) {
                                        const R cw = vcp::tsparse_scalar::real_dot_value(w, w2) / (rind * rind);
                                        for (std::size_t ii = 0; ii < n; ii++)
                                            w2[ii] -= T(cw) * w[ii];
                                    }
                                }
                                rind2 = vcp::tsparse_scalar::real_norm_value(w2);
                            }
                            const bool fresh_dir = is_pair
                                ? (rind >= tau_add || rind2 >= tau_add || pool.empty())
                                : (rind >= tau_add || pool.empty());
                            if (fresh_dir) {
                                pool_pair_t pp;
                                pp.theta = bre[bi];
                                pp.vec = lifted[i].v;
                                pp.res_abs = lifted[i].res_abs;
                                pp.res_rel = lifted[i].res_rel;
                                if (is_pair) {
                                    pp.im = th_i;
                                    pp.vec2 = lifted[i].v2;
                                }
                                pool.push_back(pp);
                                pool_slots += is_pair ? 2 : 1;
                                if (rind > R(0)) {
                                    std::vector<T> q2 = w;
                                    for (std::size_t ii = 0; ii < n; ii++) q2[ii] /= T(rind);
                                    pool_onb.push_back(q2);
                                }
                                if (is_pair && rind2 > R(0)) {
                                    std::vector<T> q3 = w2;
                                    for (std::size_t ii = 0; ii < n; ii++) q3[ii] /= T(rind2);
                                    pool_onb.push_back(q3);
                                }
                                if (have_pool_ref) added_inner = true;
                            } else {
                                // dependent: rediscovery iff the Rayleigh value is
                                // consistent with the matched pool member
                                std::size_t jbest = 0;
                                R best = R(0);
                                for (std::size_t j2 = 0; j2 < pool.size(); j2++) {
                                    R d2 = vcp::tsparse_scalar::real_dot_value(pool[j2].vec, lifted[i].v);
                                    if (d2 < R(0)) d2 = -d2;
                                    if (d2 > best) { best = d2; jbest = j2; }
                                }
                                R dv = th_r - vcp::tsparse_scalar::real_part(pool[jbest].theta);
                                if (dv < R(0)) dv = -dv;
                                if (is_pair) {
                                    // value consistency must include the imaginary part
                                    R dvi = th_i - pool[jbest].im;
                                    if (dvi < R(0)) dvi = -dvi;
                                    dv += dvi;
                                }
                                const R eps_match = R(10) * (lifted[i].res_abs + pool[jbest].res_abs
                                                             + rind * R(4) * val_scale_seen);
                                if (!(dv <= eps_match)) dependent_inconsistent = true;
                                // else: rediscovery of a pooled direction — ignored
                            }
                        }
                        if (dependent_inconsistent) {
                            dependent_block_solves++;
                            stable_confirms = 0;
                            if (dependent_block_solves >= 2) {
                                result.converged = false;
                                result.status = "not_converged";
                                result.failure_reason =
                                    "inner candidate not independently representable"
                                    " (dependent eigendirection with inconsistent value;"
                                    " honest not_converged)";
                                break;
                            }
                        } else if (added_inner || !pool_complete_at_start) {
                            stable_confirms = 0;
                        } else if (pool_slots >= k_eff) {
                            // this solve started orthogonal to a complete pool and
                            // added nothing certainly-inner: fresh confirmation
                            stable_confirms++;
                        }
                        if (stable_confirms >= 1 && pool_slots >= k_eff) {
                            declare_from_pool = true;
                            break;
                        }
                        // continue scanning with a new pool-orthogonal solve
                        start_new_solve = true;
                        restart_count++;
                        continue;
                    }
                } else {
                    pending_verification = false;
                }
            } else {
                pending_verification = false;
            }
        } else if (!prefix_ready && pending_verification) {
            // verification analysis contradicted readiness: reset the protocol
            // EIG-7 α: additive route の probe cycle 待ち(alpha_pending)中は
            // リセットしない(prefix_ready 経路の規律自体は不変)
            if (!alpha_pending) pending_verification = false;
        }

        // ---- EIG-7 α(additive route; G-0.1 承認)---------------------------
        // 確認ソルブ中で prefix が bound 収束に至らない解析でも、契約 C-2 の
        // certainly 判定で confirm を評価する。prefix_ready == true の解析は
        // 上の既存経路のみが処理し、本ブロックは一切関与しない(c-2)。
        // 発動条件(最小走査 §2.3): (i) 全展開到達(M >= m_limit)、
        // (ii) probe cycle 1 回経由(alpha_pending の 2 段階)。
        // EIG-7 (iii): 進捗イベント追跡(確認ソルブのみ・fast-path 無編集の
        // 読み取り検知。E1 = pool_slots 増加の外側検知 / E2 = prefix 転換数の
        // チェイス区間内新最大)
        if (pool_complete_at_start) {
            if (pool_slots > solve_seen_pool_slots) {                    // E1
                solve_seen_pool_slots = pool_slots;
                stall_last_event_rc = restart_count;
                solve_max_prefix_conv = 0;   // チェイス区切り
            }
            std::size_t nconv_pfx = 0;
            for (std::size_t i = 0; i < prefix_blocks.size(); i++)
                if (bconv[prefix_blocks[i]]) nconv_pfx++;
            if (nconv_pfx > solve_max_prefix_conv) {                     // E2
                solve_max_prefix_conv = nconv_pfx;
                stall_last_event_rc = restart_count;
            }
        }
        if (!declared_converged && !declare_from_pool && pool_complete_at_start
            && !prefix_ready && have_final_analysis && M >= m_limit
            && restart_count - stall_last_event_rc >= alpha_stall_engage) {
            if (!alpha_pending) {
                alpha_pending = true;
                pending_verification = true;   // 次の展開に probe 折込(既存機構を使用)
            } else {
                // 候補全景(全ブロックの値 + bound 収束/被覆フラグ)を共有
                // certainly 検査へ(比較は共有ヘルパの呼び出しのみ — B-43)。
                // locked = pool(エントリ単位、複素対は im > 0 で 2 スロット)。
                std::vector<R> ar3, ai3;
                std::vector<bool> ac3;
                for (std::size_t bi = 0; bi < nb; bi++) {
                    const bool acct = bconv[bi] || covered[bi];
                    ar3.push_back(vcp::tsparse_scalar::real_part(bre[bi]));
                    ai3.push_back(vcp::tsparse_scalar::real_part(bim[bi]));
                    ac3.push_back(acct);
                    if (bsize[bi] == 2) {
                        ar3.push_back(vcp::tsparse_scalar::real_part(bre[bi]));
                        ai3.push_back(-vcp::tsparse_scalar::real_part(bim[bi]));
                        ac3.push_back(acct);
                    }
                }
                std::vector<R> pv3, pim3;
                for (std::size_t i = 0; i < pool.size(); i++) {
                    pv3.push_back(vcp::tsparse_scalar::real_part(pool[i].theta));
                    pim3.push_back(pool[i].im);
                }
                bool c2_alpha;
                if (allow_pairs) {
                    c2_alpha = vcp::tsparse::honest_termination_check_complex_pairs_<R>(
                        ar3, ai3, ac3, pv3, pim3, k_eff, target, shift_val);
                } else {
                    c2_alpha = vcp::tsparse::honest_termination_check_complex_<R>(
                        ar3, ai3, ac3, pv3, k_eff, target, shift_val);
                }
                if (c2_alpha) {
                    // certainly-inner な未収束候補が不存在: pool からの宣言へ。
                    // 返却対の C-1 厳密再検査は従来の declare_from_pool 経路が行う。
                    declare_from_pool = true;
                    break;
                }
                // 不成立(certainly-inner が見えている): 従来どおり継続。
                // 次の全展開解析で probe cycle からやり直す。
                alpha_pending = false;
                pending_verification = false;
            }
        }

        // ---- EIG-7 β(per-solve cap; G-0.1 承認 cap=8)-----------------------
        // 確認ソルブが cap を超えて空転する場合は当該ソルブを放棄し、新しい
        // pool 直交ソルブで再試行する(正直さ不変・総予算 B-1 で必ず停止)。
        if (!declared_converged && !declare_from_pool && pool_complete_at_start
            && restart_count - solve_start_restart >= confirm_solve_cap) {
            start_new_solve = true;
            restart_count++;
            continue;
        }

        if (budget_exhausted || mv_count >= max_mv) {
            budget_exhausted = true;
            result.converged = false;
            break;
        }
        if (subspace_exhausted) {
            // no expansion direction remains: analyzing the same matrix again
            // cannot change anything — honest failure (never a livelock)
            result.converged = false;
            result.status = "failed";
            result.failure_reason =
                "invariant subspace exhausted without an honest converged declaration";
            break;
        }

        // ---- PHASE 7: contraction (window = prefix ∪ converged ∪ padding) ---
        {
            // choose keep columns: leading blocks after reorder up to
            // want_cols, then extend for keep-together and keep-more handled
            // above; additionally include all bound-converged blocks that are
            // already inside the leading region.  We contract simply to the
            // leading `Kcols` columns of the reordered form: the reorder pass
            // has already moved the selected prefix to the front, converged
            // blocks among them (B-15: nothing converged is dropped unless the
            // m_limit-1 cap forces it — counted above and reported).
            std::size_t Kcols = 0;
            {
                std::size_t want_cols = k + ((k / 2 + 1 > 2) ? k / 2 + 1 : 2);
                if (keep_more_cols > want_cols) want_cols = keep_more_cols;   // G-2.1 c-1
                // B-15: extend to cover every bound-converged block
                for (std::size_t bi = 0; bi < nb; bi++) {
                    const std::size_t bend = bpos[bi] + bsize[bi];
                    if (bconv[bi] && bend > want_cols) want_cols = bend;
                }
                // keep-together (D3-5): widen to the boundary of a straddled block
                for (std::size_t bi = 0; bi < nb; bi++) {
                    const std::size_t bend = bpos[bi] + bsize[bi];
                    if (bpos[bi] < want_cols && bend > want_cols) want_cols = bend;
                }
                if (want_cols > M - 1) want_cols = M - 1;   // >=1 expansion slot
                // final alignment DOWN to a block boundary (never split a 2x2;
                // an m_limit-1 cap clip is the reported rare corner)
                std::size_t aligned = 0;
                for (std::size_t bi = 0; bi < nb; bi++) {
                    const std::size_t bend = bpos[bi] + bsize[bi];
                    if (bend <= want_cols) aligned = bend;
                    else break;
                }
                Kcols = aligned;   // 0 => degenerate fresh restart below
            }

            if (Kcols == 0) {
                // degenerate: fresh restart from current v_next (or a fresh
                // direction) — factorization restarts empty
                V.clear();
                S.clear();
                coup.clear();
                if (v_next.empty()) {
                    std::vector<T> fd;
                    if (!kd::fresh_orthogonal_direction(V, 0, n, seed, floor_tol, fd)) {
                        result.converged = false;
                        result.status = "failed";
                        result.failure_reason = "no restart direction available";
                        break;
                    }
                    v_next = fd;
                }
                restart_count++;
                continue;
            }

            // V_new = V * U(:, 0..Kcols-1)
            std::vector<std::vector<T> > Vn(Kcols, std::vector<T>(n, T(0)));
            for (std::size_t j = 0; j < Kcols; j++)
                for (std::size_t i = 0; i < M; i++) {
                    const T u = U[i][j];
                    const std::vector<T>& vi = V[i];
                    for (std::size_t ii = 0; ii < n; ii++) Vn[j][ii] += vi[ii] * u;
                }
            // re-orthonormalize Vn (MGS QR) and push the R-correction into
            // S and c so the relation stays exact bookkeeping:
            //   V = Q Rc  =>  A Q = Q (Rc Tblk Rc^{-1}) + v (cT Rc^{-1})
            std::vector<std::vector<T> > Rc(Kcols, std::vector<T>(Kcols, T(0)));
            bool orth_ok = true;
            for (std::size_t j = 0; j < Kcols && orth_ok; j++) {
                for (std::size_t i = 0; i < j; i++) {
                    const R c2 = vcp::tsparse_scalar::real_dot_value(Vn[i], Vn[j]);
                    Rc[i][j] = T(c2);
                    for (std::size_t ii = 0; ii < n; ii++) Vn[j][ii] -= T(c2) * Vn[i][ii];
                }
                const R nv = vcp::tsparse_scalar::real_norm_value(Vn[j]);
                if (!(nv > R(0))) { orth_ok = false; break; }
                Rc[j][j] = T(nv);
                for (std::size_t ii = 0; ii < n; ii++) Vn[j][ii] /= T(nv);
            }
            if (!orth_ok) {
                result.converged = false;
                result.status = "failed";
                result.failure_reason = "contracted basis not certifiably independent";
                break;
            }
            // Tblk = leading Kcols block of Tm ; S_new = Rc Tblk Rc^{-1}
            // c_new = cT(0..Kcols-1) Rc^{-1}
            std::vector<std::vector<T> > Tb(Kcols, std::vector<T>(Kcols, T(0)));
            for (std::size_t i = 0; i < Kcols; i++)
                for (std::size_t j = 0; j < Kcols; j++) Tb[i][j] = Tm[i][j];
            // X = Tblk Rc^{-1}: solve X Rc = Tblk column-wise (Rc upper tri)
            std::vector<std::vector<T> > X(Kcols, std::vector<T>(Kcols, T(0)));
            for (std::size_t i = 0; i < Kcols; i++) {
                for (std::size_t j = 0; j < Kcols; j++) {
                    T s = Tb[i][j];
                    for (std::size_t t = 0; t < j; t++) s -= X[i][t] * Rc[t][j];
                    X[i][j] = s / Rc[j][j];
                }
            }
            std::vector<std::vector<T> > Sn(Kcols, std::vector<T>(Kcols, T(0)));
            for (std::size_t i = 0; i < Kcols; i++)
                for (std::size_t j = 0; j < Kcols; j++) {
                    T s(0);
                    for (std::size_t t = i; t < Kcols; t++) s += Rc[i][t] * X[t][j];
                    Sn[i][j] = s;
                }
            std::vector<T> cn(Kcols, T(0));
            for (std::size_t j = 0; j < Kcols; j++) {
                T s = cT[j];
                for (std::size_t t = 0; t < j; t++) s -= cn[t] * Rc[t][j];
                cn[j] = s / Rc[j][j];
            }

            V.swap(Vn);
            S.swap(Sn);
            coup.swap(cn);
            // v_next unchanged (still orthogonal to span(V) up to eps)
            restart_count++;
        }
    }   // main loop

    // ---- declaration from the pool (confirmed by a fresh orthogonal solve) --
    if (declare_from_pool && !declared_converged && allow_pairs) {
        // EIG-4 T-3: slot-aware pool declaration (pairs occupy 2 slots;
        // keep-together on a straddled pair returns k_eff + 1 values).
        std::vector<R> pre3, pim3;
        for (std::size_t i = 0; i < pool.size(); i++) {
            pre3.push_back(vcp::tsparse_scalar::real_part(pool[i].theta));
            pim3.push_back(pool[i].im);
        }
        const std::vector<std::size_t> sel =
            kd::select_slot_indices_pairs<R>(pre3, pim3, k_eff, target, shift_val);
        std::size_t need_mv = 0;
        for (std::size_t i = 0; i < sel.size(); i++)
            need_mv += (pool[sel[i]].im > R(0)) ? 2 : 1;
        if (!sel.empty() && mv_count + need_mv <= max_mv) {
            bool all_ok = true;
            struct outp_t { std::vector<T> u, v; R ra, rrel; bool pair; T re, im; };
            std::vector<outp_t> outp;
            for (std::size_t i = 0; i < sel.size() && all_ok; i++) {
                const pool_pair_t& pp = pool[sel[i]];
                outp_t o;
                o.pair = (pp.im > R(0));
                o.re = pp.theta;
                o.im = T(pp.im);
                if (!o.pair) {
                    std::vector<T> Av;
                    apply(pp.vec, Av);
                    mv_count++;
                    std::vector<T> rr(n);
                    for (std::size_t ii = 0; ii < n; ii++)
                        rr[ii] = Av[ii] - pp.theta * pp.vec[ii];
                    o.ra = vcp::tsparse_scalar::real_norm_value(rr);
                    o.rrel = o.ra / (R(1) + abs_value(vcp::tsparse_scalar::real_part(pp.theta)));
                    o.u = pp.vec;
                } else {
                    std::vector<T> Au, Av2;
                    apply(pp.vec, Au);
                    apply(pp.vec2, Av2);
                    mv_count += 2;
                    R rs(0);
                    for (std::size_t ii = 0; ii < n; ii++) {
                        const T d1 = Au[ii]  - (pp.theta * pp.vec[ii]  - T(pp.im) * pp.vec2[ii]);
                        const T d2 = Av2[ii] - (T(pp.im) * pp.vec[ii]  + pp.theta * pp.vec2[ii]);
                        const R a1 = abs_value(d1);
                        const R a2 = abs_value(d2);
                        rs += a1 * a1 + a2 * a2;
                    }
                    o.ra = sqrt_value(rs);
                    const R thr = vcp::tsparse_scalar::real_part(pp.theta);
                    o.rrel = o.ra / (R(1) + sqrt_value(thr * thr + pp.im * pp.im));
                    o.u = pp.vec;
                    o.v = pp.vec2;
                }
                if (!((o.ra <= tol) || (o.rrel <= tol))) { all_ok = false; break; }
                outp.push_back(o);
            }
            if (all_ok) {
                declared_converged = true;
                diag.c2_evidence.fresh = true;
                for (std::size_t i = 0; i < outp.size(); i++) {
                    const outp_t& o = outp[i];
                    if (!o.pair) {
                        result.eigenvalues.push_back(o.re);
                        result.eigenvalues_imag.push_back(T(0));
                        result.eigenvectors.push_back(o.u);
                        result.residuals_absolute.push_back(o.ra);
                        result.residuals_relative.push_back(o.rrel);
                    } else {
                        result.eigenvalues.push_back(o.re);
                        result.eigenvalues.push_back(o.re);
                        result.eigenvalues_imag.push_back(o.im);
                        result.eigenvalues_imag.push_back(-o.im);
                        result.eigenvectors.push_back(o.u);
                        result.eigenvectors.push_back(o.v);
                        result.residuals_absolute.push_back(o.ra);
                        result.residuals_absolute.push_back(o.ra);
                        result.residuals_relative.push_back(o.rrel);
                        result.residuals_relative.push_back(o.rrel);
                        returned_pair_count++;
                    }
                }
            } else {
                result.converged = false;
                result.status = "residual_check_failed";
                result.failure_reason =
                    "end-of-run exact residual re-check failed on a pooled pair (C-1)";
            }
        } else if (!sel.empty()) {
            budget_exhausted = true;
            result.converged = false;
        } else {
            result.converged = false;
            result.status = "failed";
            result.failure_reason = "pool ordering shorter than k (internal)";
        }
    } else if (declare_from_pool && !declared_converged) {
        std::vector<R> pv2;
        for (std::size_t i = 0; i < pool.size(); i++)
            pv2.push_back(vcp::tsparse_scalar::real_part(pool[i].theta));
        const std::vector<std::size_t> po =
            vcp::tsparse_eigen_selection::select_eigen_indices_from_real(
                pv2, pv2.size(), target, shift_val);
        if (po.size() >= k_eff && mv_count + k_eff <= max_mv) {
            bool all_ok = true;
            std::vector<lifted_pair_t> outp(k_eff);
            for (std::size_t i = 0; i < k_eff; i++) {
                const pool_pair_t& pp = pool[po[i]];
                std::vector<T> Av;
                apply(pp.vec, Av);
                mv_count++;
                std::vector<T> rr(n);
                for (std::size_t ii = 0; ii < n; ii++)
                    rr[ii] = Av[ii] - pp.theta * pp.vec[ii];
                const R ra = vcp::tsparse_scalar::real_norm_value(rr);
                const R rrel = ra / (R(1) + abs_value(vcp::tsparse_scalar::real_part(pp.theta)));
                if (!((ra <= tol) || (rrel <= tol))) { all_ok = false; break; }
                outp[i].v = pp.vec;
                outp[i].res_abs = ra;
                outp[i].res_rel = rrel;
            }
            if (all_ok) {
                declared_converged = true;
                diag.c2_evidence.fresh = true;
                for (std::size_t i = 0; i < k_eff; i++) {
                    result.eigenvalues.push_back(pool[po[i]].theta);
                    result.eigenvectors.push_back(outp[i].v);
                    result.residuals_absolute.push_back(outp[i].res_abs);
                    result.residuals_relative.push_back(outp[i].res_rel);
                }
            } else {
                result.converged = false;
                result.status = "residual_check_failed";
                result.failure_reason =
                    "end-of-run exact residual re-check failed on a pooled pair (C-1)";
            }
        } else if (po.size() >= k_eff) {
            budget_exhausted = true;
            result.converged = false;
        } else {
            result.converged = false;
            result.status = "failed";
            result.failure_reason = "pool ordering shorter than k (internal)";
        }
    }

    // =========================================================================
    // termination: evidence + factorization export, result assembly
    // =========================================================================
    diag.restart_count = restart_count;
    diag.matrix_vector_products = mv_count;
    result.matrix_vector_products = mv_count;
    result.iterations = restart_count;
    // EIG-7 β(B-44): 最終ソルブが確認ソルブならその消費分も累計に含めて公開
    if (pool_complete_at_start)
        confirm_restarts_total += restart_count - solve_start_restart;
    result.confirm_restarts = confirm_restarts_total;

    // evidence (exported on every n>0 && k>0 path)
    diag.c2_evidence.exported = true;
    diag.c2_evidence.candidate_real.clear();
    diag.c2_evidence.candidate_imag.clear();
    diag.c2_evidence.candidate_converged.clear();
    if (have_final_analysis) {
        for (std::size_t i = 0; i < final_cands.size(); i++) {
            const cand_t& c = final_cands[i];
            diag.c2_evidence.candidate_real.push_back(
                vcp::tsparse_scalar::real_part(c.re));
            diag.c2_evidence.candidate_imag.push_back(
                vcp::tsparse_scalar::real_part(c.im));
            diag.c2_evidence.candidate_converged.push_back(c.converged);
            if (c.is_complex) {
                diag.c2_evidence.candidate_real.push_back(
                    vcp::tsparse_scalar::real_part(c.re));
                diag.c2_evidence.candidate_imag.push_back(
                    -vcp::tsparse_scalar::real_part(c.im));
                diag.c2_evidence.candidate_converged.push_back(c.converged);
            }
            if (c.is_complex) {
                typedef typename vcp::eig_result<T>::eigenvalue_type EV;
                result.complex_eigenvalues.push_back(
                    EV(vcp::tsparse_scalar::real_part(c.re),
                       vcp::tsparse_scalar::real_part(c.im)));
                result.complex_eigenvalues.push_back(
                    EV(vcp::tsparse_scalar::real_part(c.re),
                       -vcp::tsparse_scalar::real_part(c.im)));
            }
        }
    }

    // terminating factorization export
    diag.final_basis = V;
    diag.final_compressed = S;
    diag.final_coupling = coup;
    diag.final_residual_vector = v_next;

    // EIG-4 T-3: complex-pair API surface -- populated only when the opt-in
    // actually returned pairs (all-real returns keep the legacy empty shape)
    if (returned_pair_count == 0) result.eigenvalues_imag.clear();
    result.complex_pair_count = returned_pair_count;

    // counts and status
    result.returned_real_count = result.eigenvalues.size();
    result.returned_complex_count = result.complex_eigenvalues.size();
    result.returned_count = result.returned_real_count + result.returned_complex_count;
    result.converged_count = result.returned_real_count;
    diag.locked_real_count = result.returned_real_count;
    diag.locked_complex_count = result.returned_complex_count / 2;

    if (!result.residuals_absolute.empty()) {
        R mx = result.residuals_absolute[0];
        for (std::size_t i = 1; i < result.residuals_absolute.size(); i++)
            if (result.residuals_absolute[i] > mx) mx = result.residuals_absolute[i];
        result.residual_norm_absolute = mx;
    }
    if (!result.residuals_relative.empty()) {
        R mx = result.residuals_relative[0];
        for (std::size_t i = 1; i < result.residuals_relative.size(); i++)
            if (result.residuals_relative[i] > mx) mx = result.residuals_relative[i];
        result.residual_norm_relative = mx;
    }
    if (!options.compute_residual_history) {
        result.residual_history_absolute.clear();
        result.residual_history_relative.clear();
    }

    if (declared_converged) {
        result.converged = true;
        result.status = "converged";
        result.message = "krylov_schur converged (C-1/C-2 verified, fresh)";
    } else if (complex_window_reject) {
        result.converged = false;
        result.status = "not_converged";
        result.failure_reason =
            "complex conjugate pairs present (" +
            std::to_string(complex_window_count) +
            ") in the target window: real-only eig_result cannot return them as"
            " converged pairs (complex eigenpair API pending); honest not_converged";
        result.message = "krylov_schur: complex pair in return window";
    } else if (budget_exhausted) {
        result.converged = false;
        result.status = "max_iter_exhausted";
        result.failure_reason =
            "matrix-vector product budget exhausted before full convergence";
        result.message = "krylov_schur: budget exhausted";
    } else {
        result.converged = false;
        if (result.status.empty()) result.status = "failed";
        if (result.failure_reason.empty())
            result.failure_reason = "eigensolver terminated without full convergence";
        result.message = "krylov_schur: not all eigenvalues converged";
    }

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
    return krylov_schur_eigs_with_diagnostics<Apply, T>(apply, n, k, options).eigs;
}

} // namespace tsparse_experimental
} // namespace vcp

#endif // VCP_TSPARSE_KRYLOV_SCHUR_HPP
