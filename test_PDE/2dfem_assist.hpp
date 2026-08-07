// test_PDE/2dfem_assist.hpp
//
// CM-1: verified upper bound of the H^1_0 projection error constant C_M
// (Liu and Oishi 2010, Theorem 3.6; notation follows section 3.2 of the
// Takayasu dissertation, equations (31)-(47)).
//
//     || u - P_h u ||_V <= C_M || f ||_X                        ... (38)
//     C_M = sqrt( (C_{h,0})^2 + kappa^2 )                       ... Thm 3.6
//
// The single public entry point is
//
//     vcp::fem2d_assist::projection_error_constant<T, DP, SP>(Th, k)
//
// whose .upper() is the guaranteed bound.  Everything else lives in
// namespace detail.
//
// Authority: sandbox/docs/design/CM-1_design_v1.0.md.
//
// Lexical policy (design 6.1): no decimal literals, no `double` / `float`
// tokens; the underlying point type is reached through typename T::base_type
// only.  sqrt and kv::constants<T>::pi() are allowed.
//
// Scope note (design 8): this is the "make it work" pass.  Memory
// reduction, blocking, sparse paths, parallelism and time optimisation are
// explicitly out of scope; the dense path below is intentionally literal.

#ifndef VCP_TEST_PDE_2DFEM_ASSIST_HPP
#define VCP_TEST_PDE_2DFEM_ASSIST_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <cmath>
#include <algorithm>

#include <kv/interval.hpp>
#include <kv/constants.hpp>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>
#include <vcp/imats.hpp>
#include <vcp/spimats.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_assemble.hpp>

namespace vcp {
namespace fem2d_assist {

namespace detail {

// ---------------------------------------------------------------------------
// scalar contract (design 6): T must be a kv::interval-like type, i.e. it must
// expose T::base_type.  lss / eigsymge are only available on imats<B, _P>, so
// a point type would not compile further down anyway; this turns that into an
// immediate, readable failure.
// ---------------------------------------------------------------------------
template <typename T>
struct interval_scalar_contract {
    typedef typename T::base_type base_type;
    static void require() {
        static_assert(sizeof(base_type) > 0,
                      "vcp::fem2d_assist: T must be a kv::interval-like type "
                      "exposing T::base_type");
    }
};

// ---------------------------------------------------------------------------
// result bundle.  The public function returns only c_m; the remaining fields
// exist so that the CM-1 audit gates (design 9) can inspect the intermediate
// quantities without re-deriving them.  keep_debug controls whether the two
// dense matrices are retained (they are not needed by the public path).
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct cm_parts {
    T c_h0;                          // (39)
    T kappa2;                        // (47)
    T c_m;                           // Theorem 3.6
    std::vector<T> lambda;           // diagonal of E from eigsymge(Q, Md, E)
    int sym_pairs_checked;           // number of (i, j), i < j, intersected
    vcp::matrix<T, DP> q_raw;        // Q BEFORE symmetrisation (debug only)
    vcp::matrix<T, DP> md;           // M_h mass matrix (debug only)

    cm_parts()
        : c_h0(), kappa2(), c_m(), lambda(), sym_pairs_checked(0),
          q_raw(), md() {}
};

// ---------------------------------------------------------------------------
// c0_element_bound: the Kikuchi-Liu bound for ONE choice of the origin
// vertex O (design 5, 5.1, 5.2).
//
//     C_0(K_h) <= (h / pi) * sqrt( nu_plus(alpha, theta) / 2 )
//     nu_plus   = 1 + alpha^2 + sqrt( 1 + 2 alpha^2 cos 2theta + alpha^4 )
//
// with h = |OA|, alpha = |OB| / |OA| in (0, 1), theta = angle AOB in (0, pi).
//
// Source: Lemma 3.2 of the dissertation, attributed to Kikuchi and Liu
// (2007).  [出典未逐語確認: the lemma is used as quoted by the dissertation;
// no verbatim comparison against the Kikuchi-Liu original was performed.]
//
// Evaluated in the squared form of design 5.2 so that the only square roots
// are the three written below.  a and b are the two other vertices; they are
// re-ordered internally so that |u| >= |v| (hence alpha <= 1).  The ordering
// test is made on the upper ends so that it is a total, deterministic order
// even when the two squared lengths overlap as intervals.
// ---------------------------------------------------------------------------
template <typename T>
T c0_element_bound(const std::array<T, 2>& o,
                   const std::array<T, 2>& a,
                   const std::array<T, 2>& b) {
    using std::sqrt;
    typedef typename T::base_type B;

    T ux = a[0] - o[0], uy = a[1] - o[1];
    T vx = b[0] - o[0], vy = b[1] - o[1];
    T su = ux * ux + uy * uy;                    // |u|^2
    T sv = vx * vx + vy * vy;                    // |v|^2

    // degenerate element (design 5.2): a vanishing edge cannot be certified
    // away from zero.  Failure-side gate: NOT certainly positive.
    if (!(su > T(0)) || !(sv > T(0)))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::fem2d_assist::c0_element_bound: degenerate element "
            "(an edge from the origin vertex is not certainly nonzero)");

    if (su.upper() < sv.upper()) {               // enforce |u| >= |v|
        T t = su; su = sv; sv = t;
        t = ux; ux = vx; vx = t;
        t = uy; uy = vy; vy = t;
    }

    const T d = ux * vx + uy * vy;               // u . v
    const T alpha2 = sv / su;                    // alpha^2
    T cos2t = (d * d) / (su * sv);               // cos^2 theta
    cos2t = cos2t + cos2t - T(1);                // cos 2theta

    T rad = T(1) + T(2) * alpha2 * cos2t + alpha2 * alpha2;
    // The exact radicand equals (1 - alpha^2)^2 + 2 alpha^2 (1 + cos 2theta)
    // and is therefore nonnegative; only outward rounding can push the lower
    // end below zero, so it is clipped back to the rigorous bound 0.
    if (rad.lower() < B(0)) {
        if (rad.upper() < B(0))
            vcp::throw_error<vcp::verification_error>(
                "vcp::fem2d_assist::c0_element_bound: nu_plus radicand "
                "certainly negative (inclusion broken)");
        rad.lower() = B(0);
    }

    const T nu_plus = T(1) + alpha2 + sqrt(rad);
    return sqrt(su) * sqrt(nu_plus / T(2)) / kv::constants<T>::pi();
}

// ---------------------------------------------------------------------------
// c_h0: C_{h,0} = max_{K_h} C_0(K_h)   ... (39)
//
// Per element the three vertex labellings are candidates and the SMALLEST
// upper bound is kept (design 5.1); across elements the LARGEST is kept.
// Both comparisons are made on the upper ends, and the selected interval is
// returned unchanged (design 5.2).
// ---------------------------------------------------------------------------
template <typename T>
T c_h0(const vcp::bfem::mesh<2, T>& Th) {
    if (Th.num_elements() <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::fem2d_assist::c_h0: mesh has no element");

    T best;
    bool have_best = false;
    for (int e = 0; e < Th.num_elements(); ++e) {
        const std::array<int, 3>& el = Th.element(e);
        const std::array<T, 2>& p0 = Th.vertex(el[0]);
        const std::array<T, 2>& p1 = Th.vertex(el[1]);
        const std::array<T, 2>& p2 = Th.vertex(el[2]);

        T c = c0_element_bound(p0, p1, p2);
        T c1 = c0_element_bound(p1, p2, p0);
        if (c1.upper() < c.upper()) c = c1;
        T c2 = c0_element_bound(p2, p0, p1);
        if (c2.upper() < c.upper()) c = c2;

        if (!have_best || best.upper() < c.upper()) {
            best = c;
            have_best = true;
        }
    }
    return best;
}

// ---------------------------------------------------------------------------
// dense material bundle after the homogeneous Dirichlet reduction
// (design 4, steps P1-P3).  Sizes:
//
//     Sr (ni x ni)   Br (ni x nb)   Gr (ni x nr)
//     Pd (nr x nr)   Md (nb x nb)   Nd (nb x nr)
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct dense_materials {
    int ni, nr, nb;
    vcp::matrix<T, DP> Sr, Br, Gr, Pd, Md, Nd;
    dense_materials() : ni(0), nr(0), nb(0), Sr(), Br(), Gr(), Pd(), Md(), Nd() {}
};

// ---------------------------------------------------------------------------
// assemble_dense: build the six materials of (46) and densify them.
//
// Space correspondence (design 3):  V_h = P^k,  W_h = RT_{k-1},
// M_h = P^{k-1}.  The KKT of (46) reduces to the pointwise per-element
// condition div p_h + f_h = 0 only when M_h = div(W_h), and
// div(RT_j) = P_j.  [出典未逐語確認: the dissertation attributes
// div(RT_j) = P_j to its reference [28]; no verbatim check was made.]
//
// Densification walks the sparse entries (COO before finalize, CSR/CSC
// after) and accumulates into matrix<T, DP> with +=; spmatrix::to_dense()
// and the policy dependent dense_matrix_type are deliberately not used
// (design 4, P3).
// ---------------------------------------------------------------------------
template <typename T, class DP, class SP>
void assemble_dense(const vcp::bfem::mesh<2, T>& Th, int k,
                    dense_materials<T, DP>& out) {
    typedef vcp::bfem::detail::spm_adapter<T, SP> adapter;

    vcp::bfem::fe_space<2, T, DP, SP> fs(Th, k);
    vcp::bfem::rt_space<2, T, DP, SP> rs(Th, k - 1);
    vcp::bfem::broken_space<2, T, DP, SP> bs(Th, k - 1);

    const int nd = fs.ndof(k);
    const int nr = rs.ndof();
    const int nb = bs.ndof();

    // interior (non-Dirichlet) index list of V_h
    const std::vector<int> bdry = fs.dofs(k).boundary_dofs();
    std::vector<int> pos(static_cast<std::size_t>(nd), -1);
    for (std::size_t t = 0; t < bdry.size(); ++t)
        pos[static_cast<std::size_t>(bdry[t])] = -2;      // mark as boundary
    int ni = 0;
    for (int i = 0; i < nd; ++i)
        if (pos[static_cast<std::size_t>(i)] != -2)
            pos[static_cast<std::size_t>(i)] = ni++;
    for (std::size_t t = 0; t < bdry.size(); ++t)
        pos[static_cast<std::size_t>(bdry[t])] = -1;

    // ni == 0 is legitimate (e.g. the 2-element unit square with k = 1): the
    // reduced V_h is then {0}, the Galerkin solution operator K is the zero
    // map and the two K-terms of Q drop out.  It is handled in compute().
    if (nr <= 0 || nb <= 0)
        vcp::throw_error<vcp::dimension_error>(
            "vcp::fem2d_assist::assemble_dense: empty W_h or M_h (nr = ", nr,
            ", nb = ", nb, ")");

    // ---- P2: sparse assembly (all six, then finalize) ----
    vcp::spmatrix<T, SP> S  = fs.stiffness(k);                       // S
    vcp::spmatrix<T, SP> Pm = rs.mass();                             // P
    vcp::spmatrix<T, SP> Mm = bs.mass();                             // M
    vcp::spmatrix<T, SP> Nm = vcp::bfem::assemble_div_mass(bs, rs);  // N
    vcp::spmatrix<T, SP> Bm = vcp::bfem::assemble_mixed_mass(fs, k, bs);  // B
    vcp::spmatrix<T, SP> Xm = vcp::bfem::assemble_cross_grad(rs, fs, k);  // G^T
    S.finalize();
    Pm.finalize();
    Mm.finalize();
    Nm.finalize();
    Bm.finalize();
    Xm.finalize();

    // ---- P3: densify with the Dirichlet reduction folded in ----
    out.ni = ni;
    out.nr = nr;
    out.nb = nb;
    out.Sr.zeros(ni, ni);
    out.Br.zeros(ni, nb);
    out.Gr.zeros(ni, nr);
    out.Pd.zeros(nr, nr);
    out.Md.zeros(nb, nb);
    out.Nd.zeros(nb, nr);

    vcp::matrix<T, DP>& Sr = out.Sr;
    vcp::matrix<T, DP>& Br = out.Br;
    vcp::matrix<T, DP>& Gr = out.Gr;
    vcp::matrix<T, DP>& Pd = out.Pd;
    vcp::matrix<T, DP>& Md = out.Md;
    vcp::matrix<T, DP>& Nd = out.Nd;

    adapter::for_each_entry(S, [&](int i, int j, const T& v) {
        const int a = pos[static_cast<std::size_t>(i)];
        const int b = pos[static_cast<std::size_t>(j)];
        if (a >= 0 && b >= 0) Sr(a, b) += v;
    });
    adapter::for_each_entry(Bm, [&](int i, int j, const T& v) {
        const int a = pos[static_cast<std::size_t>(i)];
        if (a >= 0) Br(a, j) += v;
    });
    // Xm has rows = RT, columns = fe; Gr = Xm[:, interior]^T (design 4 P3;
    // same orientation as bfem_rt_e2e_tests.cpp L368 Gr.at(a,j) = X.at(j,in[a]))
    adapter::for_each_entry(Xm, [&](int i, int j, const T& v) {
        const int a = pos[static_cast<std::size_t>(j)];
        if (a >= 0) Gr(a, i) += v;
    });
    adapter::for_each_entry(Pm, [&](int i, int j, const T& v) { Pd(i, j) += v; });
    adapter::for_each_entry(Mm, [&](int i, int j, const T& v) { Md(i, j) += v; });
    adapter::for_each_entry(Nm, [&](int i, int j, const T& v) { Nd(i, j) += v; });
}

// ---------------------------------------------------------------------------
// compute: the whole pipeline of design 4 (P1-P7).
//
// keep_debug retains q_raw (Q BEFORE the symmetry intersection) and md, which
// the CM-1 audit gates need; the public entry point passes false.
// ---------------------------------------------------------------------------
template <typename T, class DP, class SP>
void compute(const vcp::bfem::mesh<2, T>& Th, int k,
             cm_parts<T, DP>& out, bool keep_debug) {
    interval_scalar_contract<T>::require();
    if (k < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::fem2d_assist: k must be >= 1 (got ", k, ")");

    // C_{h,0} first: a degenerate element is reported as invalid_argument
    // before any space is built (design 5.2).
    out.c_h0 = c_h0(Th);

    dense_materials<T, DP> mat;
    assemble_dense<T, DP, SP>(Th, k, mat);
    const int ni = mat.ni;
    const int nb = mat.nb;

    // ---- P4: K = Sr^{-1} Br, the Galerkin solution operator (ni x nb) ----
    // Sr is SPD; lss certifies through its ||RA - I|| < 1 test and throws
    // vcp::verification_error otherwise (design 4 P4, design 7).
    vcp::matrix<T, DP> K;
    if (ni > 0) K = lss(mat.Sr, mat.Br);

    // ---- P5: H, the equilibrated flux operator (nr x nb) ----
    // Design 4 P5 solves the saddle system of the dissertation's variant a)
    //
    //     [ Pd   Nd^T ] [ H ]   [   0  ]
    //     [ Nd    0   ] [ L ] = [ -Md  ]      (nb right hand sides)
    //
    // which is equivalent to the explicit form
    // H = -P^{-1} N^T (N P^{-1} N^T)^{-1} M of (46).  The nonsingularity of
    // N P^{-1} N^T assumed just after (46) is discharged by the ||RA - I|| < 1
    // test inside lss; if it fails, vcp::verification_error propagates.
    // Layout matches bfem_rt_e2e_tests.cpp L374-387.
    const int nr = mat.nr;
    vcp::matrix<T, DP> H;
    {
        vcp::matrix<T, DP> A, rhs;
        A.zeros(nr + nb, nr + nb);
        rhs.zeros(nr + nb, nb);
        for (int i = 0; i < nr; ++i)
            for (int j = 0; j < nr; ++j) A(i, j) = mat.Pd(i, j);
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nr; ++j) {
                A(nr + i, j) = mat.Nd(i, j);
                A(j, nr + i) = mat.Nd(i, j);
            }
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nb; ++j) rhs(nr + i, j) = -mat.Md(i, j);

        vcp::matrix<T, DP> sol = lss(A, rhs);
        H.zeros(nr, nb);
        for (int i = 0; i < nr; ++i)
            for (int j = 0; j < nb; ++j) H(i, j) = sol(i, j);
    }

    // ---- P6: Q = K^T (Sr K) + H^T (Pd H) - 2 K^T (Gr H)   ... (46) ----
    vcp::matrix<T, DP> Q = transpose(H) * (mat.Pd * H);
    if (ni > 0) {
        vcp::matrix<T, DP> Kt = transpose(K);
        Q = Q + Kt * (mat.Sr * K);
        vcp::matrix<T, DP> C = Kt * (mat.Gr * H);
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nb; ++j) Q(i, j) -= C(i, j) + C(i, j);
    }
    if (keep_debug) out.q_raw = Q;

    // The exact Q is symmetric, so Q(i, j) and Q(j, i) enclose the same real
    // number; the intersection is the sharper valid enclosure of both.  An
    // empty intersection means the computed enclosure does not contain the
    // exact value, i.e. the inclusion is broken -- no claim is made
    // (design 4 P6, design 7).
    typedef typename T::base_type B;
    out.sym_pairs_checked = 0;
    for (int i = 0; i < nb; ++i)
        for (int j = i + 1; j < nb; ++j) {
            using std::min;
            using std::max;
            const B lo = max(Q(i, j).lower(), Q(j, i).lower());
            const B up = min(Q(i, j).upper(), Q(j, i).upper());
            if (up < lo)
                vcp::throw_error<vcp::verification_error>(
                    "vcp::fem2d_assist: Q and Q^T do not overlap at (", i, ", ",
                    j, "): [", Q(i, j).lower(), ", ", Q(i, j).upper(),
                    "] vs [", Q(j, i).lower(), ", ", Q(j, i).upper(), "]");
            Q(i, j).lower() = lo;
            Q(i, j).upper() = up;
            Q(j, i) = Q(i, j);
            ++out.sym_pairs_checked;
        }

    if (keep_debug) out.md = mat.Md;

    // ---- P7: kappa^2 = max eigenvalue of Q f_v = lambda Md f_v ... (47) ----
    // eigsymge returns E diagonal with an enclosure of every generalized
    // eigenvalue; Md (the M_h mass matrix) is SPD.  Its own certification
    // (||YBX - I|| < 1) throws vcp::verification_error on failure.
    vcp::matrix<T, DP> E;
    eigsymge(Q, mat.Md, E);

    out.lambda.resize(static_cast<std::size_t>(nb));
    B kappa2u = E(0, 0).upper();
    for (int i = 0; i < nb; ++i) {
        out.lambda[static_cast<std::size_t>(i)] = E(i, i);
        if (kappa2u < E(i, i).upper()) kappa2u = E(i, i).upper();
    }
    if (kappa2u < B(0))
        vcp::throw_error<vcp::verification_error>(
            "vcp::fem2d_assist: kappa^2 upper bound is negative (", kappa2u,
            "); the enclosure of (47) is not usable");

    out.kappa2 = T(kappa2u);

    // ---- C_M = sqrt( C_{h,0}^2 + kappa^2 )  ... Theorem 3.6 ----
    // C_{h,0} comes from Lemma 3.2, which bounds the P^0 projection error.
    // For k >= 2 the space M_h = P^{k-1} contains P^0, so its projection
    // error is not larger and the bound stays VALID, though not sharp.
    // [出典未逐語確認: this reuse of Lemma 3.2 for P^{k-1} is the design's own
    // argument (design 3.1), with no source; C_M therefore remains O(h)
    // regardless of k.]
    {
        using std::sqrt;
        out.c_m = sqrt(out.c_h0 * out.c_h0 + out.kappa2);
    }
}

} // namespace detail

// ---------------------------------------------------------------------------
// projection_error_constant (design 6)
//
// Returns an enclosure of C_M for the pair (Th, P^k).  The GUARANTEED UPPER
// BOUND IS THE RETURN VALUE'S .upper(); the lower end carries no claim beyond
// being a valid enclosure end of the computed quantity.
//
// Throws vcp::invalid_argument (k < 1, degenerate element),
// vcp::verification_error (inclusion broken, or lss / eigsymge could not
// certify), vcp::dimension_error (internal size mismatch).
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T projection_error_constant(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::cm_parts<T, DP> parts;
    detail::compute<T, DP, SP>(Th, k, parts, false);
    return parts.c_m;
}

} // namespace fem2d_assist
} // namespace vcp

#endif // VCP_TEST_PDE_2DFEM_ASSIST_HPP
