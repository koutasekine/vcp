// test_PDE/2dfem_assist.hpp
//
// CM-1: verified upper bound of the H^1_0 projection error constant C_M
// (Liu and Oishi 2010, Theorem 3.6; notation follows section 3.2 of the
// Takayasu dissertation, equations (31)-(47)).
//
//     || u - P_h u ||_V <= C_M || f ||_X                        ... (38)
//     C_M = sqrt( (C_{h,0})^2 + kappa^2 )                       ... Thm 3.6
//
// The public entry points (CM-1R design 2) are
//
//     c0_element_bound(o, a, b)          Lemma 3.2, one vertex labelling
//     c0_element(Th, e)                  C_0(K_h) of element e
//     c_h0(Th)                           (39), k independent
//     kappa_squared(Th, k)               (47)
//     projection_error_constant(Th, k)   Theorem 3.6; .upper() is the bound
//     projection_constants(Th, k)        all three in ONE dense pass
//
// all in namespace vcp::fem2d_assist.  The pipeline itself stays in
// namespace detail.
//
// Authority: sandbox/docs/design/CM-1_design_v1.0.md and
// sandbox/docs/design/CM-1R_design_v1.0.md.
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
// vcp::compsym.  imats_assist.hpp does not pull in its own dependencies: it is
// normally reached only through vcp/matrix_assist.hpp, which includes
// vcp/vcp_metafunction.hpp (vcp::is_interval) first.  matrix.hpp does not
// include matrix_assist.hpp, so the metafunction header has to precede it here
// or vcp::is_interval is undeclared inside imats_assist.hpp.
#include <vcp/vcp_metafunction.hpp>
#include <vcp/imats_assist.hpp>
#include <vcp/spmatrix.hpp>
#include <vcp/imats.hpp>
#include <vcp/spimats.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_assemble.hpp>

namespace vcp {
namespace fem2d_assist {

// ---------------------------------------------------------------------------
// projection_constant_set (CM-1R design 2): the three constants of Theorem 3.6
// together with the generalized spectrum they came from.  Returned by
// projection_constants; kappa_squared and projection_error_constant are thin
// projections of it.
//
// The two dense debug matrices of the CM-1 bundle (Q before symmetrisation and
// M_h) are deliberately NOT here: Q before the symmetry intersection is not yet
// a valid enclosure, so it must not be reachable from a public type.  They live
// in detail::core_result instead (CM-1R design 2.2).
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct projection_constant_set {
    T c_h0;                          // (39)
    T kappa2;                        // (47)
    T c_m;                           // Theorem 3.6
    std::vector<T> lambda;           // diagonal of E from eigsymge(Q, Md, E)
    int sym_pairs_checked;           // number of (i, j), i < j, intersected

    projection_constant_set()
        : c_h0(), kappa2(), c_m(), lambda(), sym_pairs_checked(0) {}
};

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
// core_result: the public constant set plus the two dense matrices that the
// CM-1 audit gates (design 9) need in order to inspect the intermediate
// quantities without re-deriving them.  keep_debug controls whether those two
// are retained (the public entry points pass false).
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct core_result : public projection_constant_set<T, DP> {
    vcp::matrix<T, DP> q_raw;        // Q BEFORE symmetrisation (debug only)
    vcp::matrix<T, DP> md;           // M_h mass matrix (debug only)

    core_result() : projection_constant_set<T, DP>(), q_raw(), md() {}
};

} // namespace detail

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
// c0_element: C_0(K_h) for the single element e.
//
// The three vertex labellings are candidates and the SMALLEST upper bound is
// kept (design 5.1); the comparison is made on the upper ends and the selected
// interval is returned unchanged (design 5.2).
//
// e outside [0, num_elements) is vcp::invalid_argument (this is the contract of
// the function itself; c_h0 below never produces such an index).
// ---------------------------------------------------------------------------
template <typename T>
T c0_element(const vcp::bfem::mesh<2, T>& Th, int e) {
    if (e < 0 || e >= Th.num_elements())
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::fem2d_assist::c0_element: element index out of range (e = ",
            e, ", num_elements = ", Th.num_elements(), ")");

    const std::array<int, 3>& el = Th.element(e);
    const std::array<T, 2>& p0 = Th.vertex(el[0]);
    const std::array<T, 2>& p1 = Th.vertex(el[1]);
    const std::array<T, 2>& p2 = Th.vertex(el[2]);

    T c = c0_element_bound(p0, p1, p2);
    T c1 = c0_element_bound(p1, p2, p0);
    if (c1.upper() < c.upper()) c = c1;
    T c2 = c0_element_bound(p2, p0, p1);
    if (c2.upper() < c.upper()) c = c2;
    return c;
}

// ---------------------------------------------------------------------------
// c_h0: C_{h,0} = max_{K_h} C_0(K_h)   ... (39)
//
// Across elements the LARGEST upper bound is kept; the comparison is made on
// the upper ends and the selected interval is returned unchanged (design 5.2).
// ---------------------------------------------------------------------------
template <typename T>
T c_h0(const vcp::bfem::mesh<2, T>& Th) {
    if (Th.num_elements() <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::fem2d_assist::c_h0: mesh has no element");

    T best;
    bool have_best = false;
    for (int e = 0; e < Th.num_elements(); ++e) {
        T c = c0_element(Th, e);
        if (!have_best || best.upper() < c.upper()) {
            best = c;
            have_best = true;
        }
    }
    return best;
}

namespace detail {

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
// assemble_dense_materials: build the six materials of (46) and densify them.
//
// Space correspondence (design 3):  V_h = P^k,  W_h = RT_{k-1},
// M_h = P^{k-1}.  The KKT of (46) reduces to the pointwise per-element
// condition div p_h + f_h = 0 only when M_h = div(W_h), and
// div(RT_j) = P_j.  [出典未逐語確認: the dissertation attributes
// div(RT_j) = P_j to its reference [28]; no verbatim check was made.]
//
// Densification is matrix<T, DP>::operator=(const spmatrix&) (CM-1R R2, SPC-2):
// every source is finalized first, so each (i, j) carries exactly one value and
// the assignment reproduces the CM-1 `+=` accumulation onto a zeroed matrix
// entry for entry.
// ---------------------------------------------------------------------------
template <typename T, class DP, class SP>
void assemble_dense_materials(const vcp::bfem::mesh<2, T>& Th, int k,
                    dense_materials<T, DP>& out) {
    vcp::bfem::fe_space<2, T, DP, SP> fs(Th, k);
    vcp::bfem::rt_space<2, T, DP, SP> rs(Th, k - 1);
    vcp::bfem::broken_space<2, T, DP, SP> bs(Th, k - 1);

    const int nd = fs.ndof(k);
    const int nr = rs.ndof();
    const int nb = bs.ndof();

    // homogeneous Dirichlet reduction of V_h (CM-1R R1): dirichlet_reduction
    // numbers the unconstrained dofs in increasing order of their global index
    // and maps the constrained ones to -1, which is exactly the index map the
    // CM-1 `pos` array carried.  The reduction itself stays SPARSE.
    vcp::bfem::dirichlet_reduction<T, DP, SP> dr(nd, fs.dofs(k).boundary_dofs());
    const int ni = dr.reduced_size();

    // ni == 0 is legitimate (e.g. the 2-element unit square with k = 1): the
    // reduced V_h is then {0}, the Galerkin solution operator K is the zero
    // map and the two K-terms of Q drop out.  It is handled in
    // projection_constants_core().
    if (nr <= 0 || nb <= 0)
        vcp::throw_error<vcp::dimension_error>(
            "vcp::fem2d_assist::assemble_dense_materials: empty W_h or M_h (nr = ", nr,
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

    // ---- P3a: the Dirichlet reduction, still sparse (CM-1R R1 / R3) ----
    // Xm has rows = RT, columns = fe, and Gr = Xm[:, interior]^T (design 4 P3;
    // same orientation as bfem_rt_e2e_tests.cpp L368 Gr.at(a,j) = X.at(j,in[a])).
    // Transposing FIRST turns the fe index into the row index, so the interior
    // selection is the same row filter the other two use:
    //     transpose(Xm)(j, i) = Xm(i, j)  ==>  Gr(a, i) = Xm(i, interior(a)).
    vcp::spmatrix<T, SP> Ssp = dr.reduce(S);                    // ni x ni
    vcp::spmatrix<T, SP> Bsp = dr.reduce_rows(Bm);              // ni x nb
    vcp::spmatrix<T, SP> Gsp = dr.reduce_rows(transpose(Xm));   // ni x nr
    Ssp.finalize();
    Bsp.finalize();
    Gsp.finalize();

    // ---- P3b: densify (every source is finalized) ----
    out.ni = ni;
    out.nr = nr;
    out.nb = nb;
    out.Sr = Ssp;
    out.Br = Bsp;
    out.Gr = Gsp;
    out.Pd = Pm;
    out.Md = Mm;
    out.Nd = Nm;
}

// ---------------------------------------------------------------------------
// projection_constants_core: the whole pipeline of design 4 (P1-P7).
//
// keep_debug retains q_raw (Q BEFORE the symmetry intersection) and md, which
// the CM-1 audit gates need; the public entry points pass false.
// ---------------------------------------------------------------------------
template <typename T, class DP, class SP>
void projection_constants_core(const vcp::bfem::mesh<2, T>& Th, int k,
                               core_result<T, DP>& out, bool keep_debug) {
    interval_scalar_contract<T>::require();
    if (k < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::fem2d_assist: k must be >= 1 (got ", k, ")");

    // C_{h,0} first: a degenerate element is reported as invalid_argument
    // before any space is built (design 5.2).
    out.c_h0 = c_h0(Th);

    dense_materials<T, DP> mat;
    assemble_dense_materials<T, DP, SP>(Th, k, mat);
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
        A({0, nr - 1}, {0, nr - 1}) = mat.Pd;
        A({nr, nr + nb - 1}, {0, nr - 1}) = mat.Nd;
        A({0, nr - 1}, {nr, nr + nb - 1}) = transpose(mat.Nd);
        rhs({nr, nr + nb - 1}, {}) = -mat.Md;

        vcp::matrix<T, DP> sol = lss(A, rhs);
        H = sol({0, nr - 1}, {});
    }

    // ---- P6: Q = K^T (Sr K) + H^T (Pd H) - 2 K^T (Gr H)   ... (46) ----
    vcp::matrix<T, DP> Q = transpose(H) * (mat.Pd * H);
    if (ni > 0) {
        vcp::matrix<T, DP> Kt = transpose(K);
        Q = Q + Kt * (mat.Sr * K);
        vcp::matrix<T, DP> C = Kt * (mat.Gr * H);
        // componentwise Q(i, j) - (C(i, j) + C(i, j)): the same operation
        // sequence as the CM-1 loop, so the rounding order is unchanged.
        Q = Q - (C + C);
    }
    if (keep_debug) out.q_raw = Q;

    // The exact Q is symmetric, so Q(i, j) and Q(j, i) enclose the same real
    // number; the intersection is the sharper valid enclosure of both.  An
    // empty intersection means the computed enclosure does not contain the
    // exact value, i.e. the inclusion is broken -- no claim is made
    // (design 4 P6, design 7).
    //
    // vcp::compsym (CM-1R R7) is that operation verbatim: kv's overlap is
    // max(lower) <= min(upper) and its intersect is [max(lower), min(upper)],
    // and a non-overlapping pair raises vcp::verification_error.  It visits the
    // same (i, j), i < j, pairs, so the count below is the number of pairs it
    // examined.
    typedef typename T::base_type B;
    vcp::compsym(Q);
    out.sym_pairs_checked = nb * (nb - 1) / 2;

    if (keep_debug) out.md = mat.Md;

    // ---- P7: kappa^2 = max eigenvalue of Q f_v = lambda Md f_v ... (47) ----
    // eigsymge returns E diagonal with an enclosure of every generalized
    // eigenvalue; Md (the M_h mass matrix) is SPD.  Its own certification
    // (||YBX - I|| < 1) throws vcp::verification_error on failure.
    vcp::matrix<T, DP> E;
    eigsymge(Q, mat.Md, E);

    // diag(E) extracts the diagonal exactly into an nb x 1 column (CM-1R R8),
    // leaving a single scan for the maximum.
    vcp::matrix<T, DP> ev = diag(E);

    out.lambda.resize(static_cast<std::size_t>(nb));
    B kappa2u = ev(0, 0).upper();
    for (int i = 0; i < nb; ++i) {
        out.lambda[static_cast<std::size_t>(i)] = ev(i, 0);
        if (kappa2u < ev(i, 0).upper()) kappa2u = ev(i, 0).upper();
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
// projection_constants (design 6, CM-1R design 2)
//
// The three constants of Theorem 3.6 for the pair (Th, P^k), obtained in ONE
// pass through the dense path.  kappa_squared and projection_error_constant
// below are thin projections of this function; calling both of them costs two
// dense passes, so ask for the set when both are wanted.
//
// Throws vcp::invalid_argument (k < 1, degenerate element),
// vcp::verification_error (inclusion broken, or lss / eigsymge could not
// certify), vcp::dimension_error (internal size mismatch).
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
projection_constant_set<T, DP>
projection_constants(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::core_result<T, DP> r;
    detail::projection_constants_core<T, DP, SP>(Th, k, r, false);
    return static_cast<const projection_constant_set<T, DP>&>(r);
}

// ---------------------------------------------------------------------------
// kappa_squared: kappa^2 = lambda_max(Q, M)   ... (47)
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T kappa_squared(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::core_result<T, DP> r;
    detail::projection_constants_core<T, DP, SP>(Th, k, r, false);
    return r.kappa2;
}

// ---------------------------------------------------------------------------
// projection_error_constant (design 6)
//
// Returns an enclosure of C_M for the pair (Th, P^k).  The GUARANTEED UPPER
// BOUND IS THE RETURN VALUE'S .upper(); the lower end carries no claim beyond
// being a valid enclosure end of the computed quantity.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T projection_error_constant(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::core_result<T, DP> r;
    detail::projection_constants_core<T, DP, SP>(Th, k, r, false);
    return r.c_m;
}

} // namespace fem2d_assist
} // namespace vcp

#endif // VCP_TEST_PDE_2DFEM_ASSIST_HPP
