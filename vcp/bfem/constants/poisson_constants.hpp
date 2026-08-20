// vcp/bfem/constants/poisson_constants.hpp
//
// Verified constants of the Laplacian (Poisson) problem: the H^1_0 projection
// error constant C_h of Liu and Oishi 2010 (CM-1), the piecewise constant
// L^2 projection constant C_0 h and the hypercircle quantity kappa_h^2 it is
// built from.  Symbols and equation numbers follow chapter 8 of
//
//     劉 雪峰・関根 晃太「偏微分方程式の精度保証付き数値計算法」,
//     『精度保証付き数値計算の基礎』第 8 章
//
// which is the primary source of this file (CM-1S design 1).
//
//     || u - u_h ||_V <= C_h || f ||_X               ... Theorem 8.3 / (8.25)
//     C_h := sqrt( C_0^2 h^2 + kappa_h^2 )           ... (8.25)
//
// with u the solution of the variational problem (8.9), u_h = P_h u its
// Galerkin approximation (8.10) and P_h the projection (8.11).
//
// VER-0: this header is the verbatim move of the CM-1 header formerly under
// test_PDE/ into namespace vcp::bfem::constants, with the public symbols
// renamed as below.  The computation code is unchanged; the detail layer and
// the struct member names (c_h0 / kappa2 / c_m / lambda / sym_pairs_checked)
// keep their old names.  Old names are kept in the table for cross-reference
// with the CM-1 / CM-1R / CM-1S / CM-1T documents:
//
//     old name                    new public entry point
//     projection_error_constant   ritz_projection_error_constant_h01(Th, k)  C_h
//     projection_constants        ritz_projection_constants_h01(Th, k)  the set
//     projection_constant_set     ritz_projection_constant_set_h01<T, DP>
//     kappa_squared               hypercircle_kappa_h01_squared(Th, k)  (8.21)
//     c_h0                        l2_projection_error_constant(Th)  C_0 h,
//                                 (8.19)(8.20), k independent
//     c0_element                  l2_projection_element_constant(Th, e)
//     c0_element_bound            l2_projection_element_bound(o, a, b)
//
// VER-1 boundary-condition ledger (owner ruling R23): the BC assumption of
// an entry point is part of its NAME -- the BC tag sits right after the
// concept name and before the form suffixes (_sq / _bound / _squared /
// _measure / ...).  _h01 marks an entry that assumes the homogeneous
// Dirichlet space H^1_0; future additions to the layer follow the same
// convention.  BC column, spelled out for every entry point of this header
// and its detail layer:
//
//     entry point                              BC assumption
//     ritz_projection_error_constant_h01       H^1_0: the Ritz projection
//                                              maps into V_h in H^1_0 and
//                                              kappa is the Dirichlet dual
//                                              bound
//     ritz_projection_error_constant_h01_sq    H^1_0: same frame (the
//                                              CONST-C1 squared composition
//                                              of detail/poisson_dict_impl)
//     hypercircle_kappa_h01_squared            H^1_0: the defining equation
//                                              of kappa_h (see the function
//                                              comment below) frames
//                                              -Delta phi = psi with
//                                              phi in H^1_0
//     poincare_constant_h01_sq_bound           H^1_0: the supplied lambda_1
//                                              is a DIRICHLET eigenvalue
//                                              (Friedrichs form
//                                              ||u|| <= C ||grad u||,
//                                              u in H^1_0)
//     poincare_constant_h01_bound              H^1_0: same
//     ritz_projection_constants_h01 /          H^1_0: the set the _h01
//     ritz_projection_constant_set_h01         entries above are thin
//                                              projections of (renamed by
//                                              CONST-D under R23: _set is
//                                              part of the carrier noun,
//                                              not a form suffix, so the
//                                              tag sits at the end)
//     l2_projection_error_constant(_sq)        BC-free: the L^2 projection
//     l2_projection_element_constant           carries no boundary
//     l2_projection_element_constants_sq       condition -- these bounds
//     l2_projection_element_bound              are valid for Dirichlet,
//     dictionary layer (dict/ registry,        Neumann and Robin problems
//     coverage tables and the resolve chain)   alike
//     cr_interpolation_constant                BC-free: element local
//     cr_projection_error_constant_sq          BC-free: element local
//     eigenvalue_lower_bound                   BC-free: general discrete
//                                              lower bound formula
//
// l2_projection_error_constant is the constant of the projection onto the
// piecewise constant space P^0: the target space is FIXED at P^0 and the
// value does not depend on k.  A (Th, d) overload with a k dependent target
// degree is reserved for the 定数層-B track.
//
// static_assert convention of this layer: the entry points below take a mesh
// and build the conforming spaces themselves, so there is nothing to assert
// here.  Whoever adds to this layer a function that takes a SPACE as an
// argument and assumes V_h in H^1_0 must place
// static_assert(vcp::bfem::space_traits<S>::h1_conforming, ...) in it (see
// space_traits.hpp).
//
// C_0 h against C_h -- the two are easy to confuse.  (8.19)(8.20) normalise by
// h, i.e. C_0 = max_K C_0^{(1)}(K) / h, so the quantity
// l2_projection_error_constant returns is C_0 h = max_K C_0^{(1)}(K) itself;
// it is the constant of the piecewise constant projection pi_{0,h} of
// (8.18)(8.19) and only the FIRST term under the square root of (8.25).  C_h
// is the projection error constant itself.
//
// Authority: sandbox/docs/design/CM-1_design_v1.0.md,
// sandbox/docs/design/CM-1R_design_v1.0.md,
// sandbox/docs/design/CM-1S_design_v1.0.md and, for the move and renaming,
// sandbox/docs/design/VER-0_design_v1.0.md.
//
// Lexical policy (design 6.1): no decimal literals, no `double` / `float`
// tokens; the underlying point type is reached through typename T::base_type
// only.  sqrt and kv::constants<T>::pi() are allowed.
//
// Scope note (design 8): this is the "make it work" pass.  Memory
// reduction, blocking, sparse paths, parallelism and time optimisation are
// explicitly out of scope; the dense path below is intentionally literal.
//
// Dense policy DP (CM-1T design 2).  The default stays
// vcp::imats<typename T::base_type>, so that the header needs no external
// BLAS.  When one is available, pass DP = vcp::pidblas.  It derives from
// vcp::imats<double, vcp::pdblas> and overrides mulmm / mul_im_m / mul_m_im /
// vmulmm / mulltmm, so that the INTERVAL products themselves are turned into
// BLAS3 calls through a mid/rad split with directed rounding.
//
//     vcp::pidblas                      <- recommended when BLAS is available
//     vcp::imats<double, vcp::pdblas>   <- NOT the same thing
//
// The second one accelerates the POINT operations only; interval x interval
// still falls back to the scalar loops of imats, which is what dominates the
// pipeline below.  Measured by the CM-1T author on the book mesh h = 1/4, one
// core, reference BLAS: k = 1/2/3 took 1.26 / 38.41 / (over 300, unfinished)
// seconds with vcp::imats<double> against 0.17 / 4.51 / 40.42 seconds with
// vcp::pidblas, the two agreeing to nine digits.  Neither the timings nor the
// agreement is contracted here: the rounding ORDER differs between policies,
// so the results are not required to be bit-identical.
#ifndef VCP_BFEM_CONSTANTS_POISSON_CONSTANTS_HPP
#define VCP_BFEM_CONSTANTS_POISSON_CONSTANTS_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <cmath>
#include <algorithm>

#include <kv/interval.hpp>
#include <kv/constants.hpp>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
// vcp::compsym.  matrix_assist.hpp is the proper entry point: it pulls in
// vcp/vcp_metafunction.hpp (vcp::is_interval) and then takes in
// vcp/imats_assist.hpp under #if defined(INTERVAL_HPP), which kv/interval.hpp
// above has already defined.  It must follow matrix.hpp (it #errors otherwise),
// which is also the order the existing test_PDE/*.cpp use.
#include <vcp/matrix_assist.hpp>
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
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// ritz_projection_constant_set_h01 (CM-1R design 2): the three constants of Theorem 8.3
// together with the generalized spectrum they came from.  Returned by
// ritz_projection_constants_h01; hypercircle_kappa_h01_squared and ritz_projection_error_constant_h01 are thin
// projections of it.
//
// The two dense debug matrices of the CM-1 bundle (Q before symmetrisation and
// the X_h mass matrix) are deliberately NOT here: Q before the symmetry
// intersection is not yet a valid enclosure, so it must not be reachable from a
// public type.  They live in detail::core_result instead (CM-1R design 2.2).
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct ritz_projection_constant_set_h01 {
    T c_h0;                          // C_0 h, (8.19)(8.20)
    T kappa2;                        // kappa_h^2, kappa_h is (8.21)
    T c_m;                           // C_h, Theorem 8.3 / (8.25)
    std::vector<T> lambda;           // diagonal of E from eigsymge(Q, Md, E)
    int sym_pairs_checked;           // number of (i, j), i < j, intersected

    ritz_projection_constant_set_h01()
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
                      "vcp::bfem::constants: T must be a kv::interval-like type "
                      "exposing T::base_type");
    }
};

// ---------------------------------------------------------------------------
// core_result: the public constant set plus the three dense matrices that the
// CM-1 audit gates (design 9) need in order to inspect the intermediate
// quantities without re-deriving them.  keep_debug controls whether those are
// retained (the public entry points pass false).
//
// q_raw and q_sym are Q on either side of vcp::compsym, so gate G2 can check
// the post-conditions of the symmetrisation itself (CM-1S design 4): q_sym must
// be exactly symmetric and must be contained in q_raw entry by entry.
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct core_result : public ritz_projection_constant_set_h01<T, DP> {
    vcp::matrix<T, DP> q_raw;        // Q BEFORE symmetrisation (debug only)
    vcp::matrix<T, DP> q_sym;        // Q AFTER  symmetrisation (debug only)
    vcp::matrix<T, DP> md;           // X_h mass matrix (debug only)

    core_result() : ritz_projection_constant_set_h01<T, DP>(), q_raw(), q_sym(), md() {}
};

} // namespace detail

// ---------------------------------------------------------------------------
// l2_projection_element_bound: an upper bound of C_0^{(1)}(K) (table 8.1, section 8.3) for
// ONE choice of the origin vertex O (design 5, 5.1, 5.2).  The symbols O, A, B,
// L, alpha, theta are those of figure 8.4 of the book.
//
//     C_0^{(1)}(K) <= (L / pi) * sqrt( nu_plus(alpha, theta) / 2 )
//     nu_plus       = 1 + alpha^2 + sqrt( 1 + 2 alpha^2 cos 2theta + alpha^4 )
//
// with L = |OA|, alpha = |OB| / |OA| in (0, 1], theta = angle AOB in (0, pi).
//
// Which norm C_0^{(1)}(K) is taken over (CM-1T design 6.1).  It is the
// SEMINORM.  Table 8.1 of the book writes the projection and its estimate as
//
//     Pi_0^{(1)} u := ( 1 / |K| ) int_K u dxdy,
//     || u - Pi_0^{(1)} u ||_{0,K} <= C_0^{(1)} | u |_{1,K}
//
// so the denominator of the supremum defining C_0^{(1)}(K) is | u |_{1,K} and
// not the full H^1 norm.  That is the reading this implementation uses; it was
// left as an open question by CM-1S, and the book settles it.  It is also what
// makes C_0^{(1)} = 1 / pi on the unit right isosceles triangle an equality
// rather than an estimate (section 8.3, and (b) at the end of this file).
//
// Which of the two bounds is implemented (CM-1S design 1.2).  Section 8.3 of
// the book gives C_0^{(1)}(K) <= (L / pi) sqrt( 1 + |cos theta| ).  For
// alpha = 1 the identity nu_plus / 2 = 1 + |cos theta| holds, so that estimate
// is exactly the alpha = 1 specialisation of the one above; for alpha < 1 the
// nu_plus form is sharper (up to 26 percent on the sampled alpha-theta grid),
// so it is the one used here.  Its source is Lemma 3.2 of the Takayasu
// dissertation, attributed there to Kikuchi and Liu (2007), and this is the
// ONLY dissertation reference left in the file (CM-1S gate S-G9 exception).
// [出典未逐語確認: the lemma is used as quoted by the dissertation; no verbatim
// comparison against the Kikuchi-Liu original was performed.]
//
// The book's other estimate, C_0^{(1)}(K) <= |AB| / j_{1,1}, is NOT used: a
// rigorous enclosure of j_{1,1} would need a decimal literal, which the lexical
// policy below forbids.
//
// theta means the same angle in both forms.  The book takes the largest
// interior angle; the reordering below makes |OB| <= |OA| within one call, and
// l2_projection_element_constant keeps the smallest of the three vertex labellings, among which the
// one placing O opposite the longest edge AB carries that largest angle.
//
// Evaluated in the squared form of design 5.2 so that the only square roots
// are the three written below.  a and b are the two other vertices; they are
// re-ordered internally so that |u| >= |v| (hence alpha <= 1).  The ordering
// test is made on the upper ends so that it is a total, deterministic order
// even when the two squared lengths overlap as intervals.
// ---------------------------------------------------------------------------
template <typename T>
T l2_projection_element_bound(const std::array<T, 2>& o,
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
            "vcp::bfem::constants::l2_projection_element_bound: degenerate element "
            "(an edge from the origin vertex is not certainly nonzero)");

    // collinear degeneracy (CM-1S design 3.2): three DISTINCT vertices can
    // still span zero area, and the test above does not see it.  Twice the
    // signed area is the cross product below; the same failure-side gate is
    // applied to it, so a sign that cannot be certified either way -- an
    // enclosure containing 0 -- is rejected.  This matches the philosophy of
    // vcp/bfem/geometry.hpp's geometry_traits<T>::sign(), with the exception
    // type kept inside the vcp::error hierarchy.  Either sign is accepted: the
    // element orientation is free, only the area has to be certainly nonzero.
    const T cross = ux * vy - uy * vx;           // 2 * signed area
    if (!(cross > T(0)) && !(cross < T(0)))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_bound: degenerate element "
            "(zero or sign-indefinite area)");

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
                "vcp::bfem::constants::l2_projection_element_bound: nu_plus radicand "
                "certainly negative (inclusion broken)");
        rad.lower() = B(0);
    }

    const T nu_plus = T(1) + alpha2 + sqrt(rad);
    return sqrt(su) * sqrt(nu_plus / T(2)) / kv::constants<T>::pi();
}

// ---------------------------------------------------------------------------
// l2_projection_element_constant: C_0^{(1)}(K) for the single element e.
//
// The three vertex labellings are candidates and the SMALLEST upper bound is
// kept (design 5.1); the comparison is made on the upper ends and the selected
// interval is returned unchanged (design 5.2).
//
// e outside [0, num_elements) is vcp::invalid_argument (this is the contract of
// the function itself; l2_projection_error_constant below never produces such an index).
// ---------------------------------------------------------------------------
template <typename T>
T l2_projection_element_constant(const vcp::bfem::mesh<2, T>& Th, int e) {
    if (e < 0 || e >= Th.num_elements())
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_constant: element index out of range (e = ",
            e, ", num_elements = ", Th.num_elements(), ")");

    const std::array<int, 3>& el = Th.element(e);
    const std::array<T, 2>& p0 = Th.vertex(el[0]);
    const std::array<T, 2>& p1 = Th.vertex(el[1]);
    const std::array<T, 2>& p2 = Th.vertex(el[2]);

    T c = l2_projection_element_bound(p0, p1, p2);
    T c1 = l2_projection_element_bound(p1, p2, p0);
    if (c1.upper() < c.upper()) c = c1;
    T c2 = l2_projection_element_bound(p2, p0, p1);
    if (c2.upper() < c.upper()) c = c2;
    return c;
}

// ---------------------------------------------------------------------------
// l2_projection_error_constant: C_0 h = max_K C_0^{(1)}(K)   ... (8.19)(8.20)
//
// This is the constant of the piecewise constant projection pi_{0,h} of (8.18):
// (8.19) reads || v - pi_{0,h} v || <= C_0 h | v |_1, and (8.20) normalises it
// as C_0 = max_K C_0^{(1)}(K) / h, so the value returned here is the
// unnormalised max_K C_0^{(1)}(K) = C_0 h.  It does not depend on k.
//
// Across elements the LARGEST upper bound is kept; the comparison is made on
// the upper ends and the selected interval is returned unchanged (design 5.2).
// ---------------------------------------------------------------------------
template <typename T>
T l2_projection_error_constant(const vcp::bfem::mesh<2, T>& Th) {
    if (Th.num_elements() <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_error_constant: mesh has no element");

    T best;
    bool have_best = false;
    for (int e = 0; e < Th.num_elements(); ++e) {
        T c = l2_projection_element_constant(Th, e);
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
// assemble_dense_materials: build the six materials of the mixed (saddle point)
// minimisation of section 8.4.3 and densify them.
//
// Space correspondence (design 3):  V_h = P^k,  W_h = RT_{k-1} contained in
// H(div, Omega) (8.7),  X_h = P^{k-1} (the piecewise polynomial space of
// section 8.4.3).  The KKT system reduces to the pointwise per-element
// condition div p_h + f_h = 0 only when X_h = div(W_h), and div(RT_j) = P_j.
// [出典未逐語確認: section 8.4.3 attributes div(RT_j) = P_j to Boffi, Brezzi
// and Fortin, "Mixed Finite Element Methods and Applications", Springer, 2013;
// no verbatim check against that book was made.]
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
    // reduced V_h is then {0}, the Galerkin solution operator of (8.10) is the
    // zero map and the two K-terms of Q drop out.  It is handled in
    // projection_constants_core().
    if (nr <= 0 || nb <= 0)
        vcp::throw_error<vcp::dimension_error>(
            "vcp::bfem::constants::assemble_dense_materials: empty W_h or X_h (nr = ", nr,
            ", nb = ", nb, ")");

    // ---- P2: sparse assembly (all six, then finalize) ----
    vcp::spmatrix<T, SP> S  = fs.stiffness(k);                       // S
    vcp::spmatrix<T, SP> Pm = rs.mass();                             // P
    vcp::spmatrix<T, SP> Mm = bs.mass();                             // M (X_h)
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
            "vcp::bfem::constants: k must be >= 1 (got ", k, ")");

    // C_0 h first: a degenerate element is reported as invalid_argument
    // before any space is built (design 5.2).
    out.c_h0 = l2_projection_error_constant(Th);

    dense_materials<T, DP> mat;
    assemble_dense_materials<T, DP, SP>(Th, k, mat);
    const int ni = mat.ni;
    const int nb = mat.nb;

    // ---- P4: K = Sr^{-1} Br, the Galerkin solution operator of (8.10) ----
    // (ni x nb; column j is the Galerkin approximation P_h of (8.11) applied to
    // the j-th basis function of X_h.)
    // Sr is SPD; lss certifies through its ||RA - I|| < 1 test and throws
    // vcp::verification_error otherwise (design 4 P4, design 7).
    vcp::matrix<T, DP> K;
    if (ni > 0) K = lss(mat.Sr, mat.Br);

    // ---- P5: H, the equilibrated flux operator (nr x nb) ----
    // The equilibrated flux p_h in W_h (div p_h + f_h = 0) is the ingredient of
    // the hypercircle equation (8.17) / (8.24), which chapter 8 attributes to
    // Prager and Synge; kappa_h of (8.21) is the bound it yields.
    //
    // Design 4 P5 solves the saddle system of section 8.4.3
    //
    //     [ Pd   Nd^T ] [ H ]   [   0  ]
    //     [ Nd    0   ] [ L ] = [ -Md  ]      (nb right hand sides)
    //
    // which is equivalent to the explicit form
    // H = -P^{-1} N^T (N P^{-1} N^T)^{-1} M.  The nonsingularity of
    // N P^{-1} N^T, assumed in the proof of Theorem 8.2, is discharged by the
    // ||RA - I|| < 1 test inside lss; if it fails, vcp::verification_error
    // propagates.  Layout matches bfem_rt_e2e_tests.cpp L374-387.
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

    // ---- P6: Q = K^T (Sr K) + H^T (Pd H) - 2 K^T (Gr H)  ... section 8.4.3 --
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
    if (keep_debug) out.q_sym = Q;

    if (keep_debug) out.md = mat.Md;

    // ---- P7: kappa_h^2 = max eigenvalue of Q f_v = lambda Md f_v ----
    // kappa_h itself is defined by (8.21); Theorem 8.2 / (8.23) is the error
    // estimate it carries.  NOTE that the book's kappa_h is NOT squared, while
    // the quantity computed here (and returned by hypercircle_kappa_h01_squared) is kappa_h^2.
    //
    // eigsymge returns E diagonal with an enclosure of every generalized
    // eigenvalue; Md (the X_h mass matrix) is SPD.  Its own certification
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
            "vcp::bfem::constants: kappa_h^2 upper bound is negative (", kappa2u,
            "); the enclosure is not usable");

    out.kappa2 = T(kappa2u);

    // ---- C_h = sqrt( (C_0 h)^2 + kappa_h^2 )  ... Theorem 8.3 / (8.25) ----
    // C_0 h is the constant of the piecewise constant projection pi_{0,h},
    // (8.18)(8.19).  For k >= 2 the space X_h = P^{k-1} contains P^0, so its
    // projection error is not larger and the bound stays VALID, though not
    // sharp.
    //
    // WHERE THE SLACK IS, AND WHERE IT IS NOT (CM-1T design 3).  Three points,
    // because the first two are the ones easy to get wrong.
    //
    // (a) The O(h) is ESSENTIAL; it is not an artefact of this implementation
    //     and it cannot be removed by raising k.  In the proof of Theorem 8.3,
    //     (8.27), the function that C_0 h is applied to is u - u~ in H^1_0:
    //     the regularity available there is H^1 and nothing better, so that
    //     term cannot decay faster than O(h) whatever the degree of X_h is.
    //     What a larger k buys is a constant factor.
    //
    // (b) On a uniform partition there is NO slack in C_0 h itself, so
    //     l2_projection_element_bound is NOT the place to sharpen.  Section 8.3 quotes
    //     Kikuchi and Liu for C_0^{(1)} = 1 / pi on the unit right isosceles
    //     triangle, and that is an EQUALITY, not an estimate.  Every element of
    //     a partition into right isosceles triangles attains it, so on such a
    //     mesh C_0 h = h / pi is the optimal value rather than a bound with
    //     room in it; the nu_plus form above reproduces exactly that at
    //     alpha = 1, theta = pi/2.
    //
    // (c) The one slack that does remain is the DEGREE of the projection: X_h
    //     is P^{k-1}, yet the constant carried for pi_{0,h} is the one of the
    //     P^0 projection.  Sharpening it needs an upper bound of
    //
    //         sup_{v in H^1(K)} ||(I - pi_K^{(m)}) v||_{0,K} / |v|_{1,K},
    //         m = k - 1 >= 1,
    //
    //     and table 8.1 of the book treats the 0-th order projection only
    //     (C_0^{(1)}, C_0^{(2)}).  Whether values for m >= 1 are published
    //     anywhere is UNKNOWN here.  C_0^{(1)} = 1 / sqrt(lambda_1) comes from
    //     the Neumann eigenvalue problem (8.8), and the analogous quantity for
    //     m >= 1 is not simply the (m+1)-st eigenvalue, because P^m is not
    //     spanned by the eigenfunctions.  This may be an open research question
    //     rather than a gap in the material at hand, so NO claim is made that
    //     it is a matter of substituting a published value later on.
    //
    // Measured on the book mesh h = 1/4 (96 elements) by
    // sandbox/probes/cm1t_u3_probe.cpp: kappa_h falls with k as 0.146500307 /
    // 0.060494965 / 0.037385599 for k = 1/2/3, while C_0 h stays at
    // 0.079577472 for every k, so C_h = 0.166718067 / 0.099961066 /
    // 0.087921880 and the share (C_0 h / C_h)^2 that C_0 h contributes rises
    // from 23 to 63 to 82 percent.  From k = 2 on, C_h is C_0 h limited:
    // refining h is what helps, not raising k.
    // [出典未逐語確認: this reuse of the pi_{0,h} constant for P^{k-1} is the
    // design's own argument (CM-1T design 3), with no source; C_h therefore
    // remains O(h) regardless of k.]
    {
        using std::sqrt;
        out.c_m = sqrt(out.c_h0 * out.c_h0 + out.kappa2);
    }
}

} // namespace detail

// ---------------------------------------------------------------------------
// ritz_projection_constants_h01 (design 6, CM-1R design 2)
//
// The three constants of Theorem 8.3 for the pair (Th, P^k), obtained in ONE
// pass through the dense path.  hypercircle_kappa_h01_squared and ritz_projection_error_constant_h01
// below are thin projections of this function; calling both of them costs two
// dense passes, so ask for the set when both are wanted.
//
// Throws vcp::invalid_argument (k < 1, degenerate element),
// vcp::verification_error (inclusion broken, or lss / eigsymge could not
// certify), vcp::dimension_error (internal size mismatch).  A degenerate mesh
// may also let vcp::bfem::degenerate_element through from the space
// construction; see ritz_projection_error_constant_h01 below.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
ritz_projection_constant_set_h01<T, DP>
ritz_projection_constants_h01(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::core_result<T, DP> r;
    detail::projection_constants_core<T, DP, SP>(Th, k, r, false);
    return static_cast<const ritz_projection_constant_set_h01<T, DP>&>(r);
}

// ---------------------------------------------------------------------------
// hypercircle_kappa_h01_squared: lambda_max(Q, M), i.e. kappa_h^2 with kappa_h of (8.21).
//
// The book's kappa_h is NOT squared: this function returns kappa_h SQUARED, so
// the quantity the book prints (for instance in table 8.2) is sqrt() of what is
// returned here.  The squared form is what (8.25) needs directly.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T hypercircle_kappa_h01_squared(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::core_result<T, DP> r;
    detail::projection_constants_core<T, DP, SP>(Th, k, r, false);
    return r.kappa2;
}

// ---------------------------------------------------------------------------
// ritz_projection_error_constant_h01 (design 6)
//
// Returns an enclosure of C_h (Theorem 8.3 / (8.25)) for the pair (Th, P^k).
// The GUARANTEED UPPER BOUND IS THE RETURN VALUE'S .upper(); the lower end
// carries NO claim beyond being a valid enclosure end of the computed
// quantity.  In particular hypercircle_kappa_h01_squared is reduced to the point interval
// T(kappa2u) built from the largest upper end of the spectrum, so the lower end
// of the value returned here is NOT a lower bound of C_h (CM-1S U2, kept as is
// by the owner's decision).
//
// Exceptions: as for ritz_projection_constants_h01 above.  Note that on a degenerate
// mesh vcp::bfem::degenerate_element may propagate from the space construction;
// it derives from std::runtime_error and is NOT part of the vcp::error
// hierarchy, so a catch on vcp::error alone does not see it.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T ritz_projection_error_constant_h01(const vcp::bfem::mesh<2, T>& Th, int k) {
    detail::core_result<T, DP> r;
    detail::projection_constants_core<T, DP, SP>(Th, k, r, false);
    return r.c_m;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

// ===========================================================================
// CONST-A additions (everything below this line; the VER-0 moved code above
// is untouched).  The CR (Crouzeix-Raviart) interpolation error constants,
// their conversion to a working scalar type, the nonconforming eigenvalue
// lower bound formula and the Poincare constant bounds derived from it.
//
// Authority: sandbox/docs/design/CONST-A_design_v1.0.md.
//
// Lexical regime of the addition (owner ruling R16): decimal digits may
// appear ONLY inside string literals of the VCP_CONSTANTS_TABLE block below.
// Bare decimal literals and `double` / `float` tokens stay forbidden in this
// whole file, comments keep carrying no values (R14: a value written twice is
// a divergence accident waiting to happen), and every other header keeps the
// traditional regime unchanged.
// ===========================================================================

#include <vcp/bfem/constants/literal.hpp>
#include <vcp/bfem/cr1/cr1_space.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// VCP_CONSTANTS_TABLE_BEGIN  (authorized decimal strings; the single source of
// truth for every certified constant of this header.  Values appear ONLY here.)
static const char cr_interpolation_constant_2d_str[] = "0.1893";
    // certified: Liu 2015 Sec. 4 (verified computation + perturbation argument)
static const char cr_interpolation_constant_3d_str[] = "0.3804";
    // certified: Liu 2015 Thm. 3.4 (proof)
// VCP_CONSTANTS_TABLE_END

// ---------------------------------------------------------------------------
// cr_interpolation_constant<D, T>: the dimensionless constant c_D of the CR
// interpolation error estimate C(K, R^D) <= c_D h, with h the MAXIMUM EDGE
// LENGTH of the element (Liu 2015 Thm. 3.4 and Sec. 5 fix that convention).
// D = 2: Liu 2015 Sec. 4; D = 3: Liu 2015 Thm. 3.4.  The value is converted
// from the table string by constant_from_string<T> (literal.hpp): a kv
// interval T receives an outward rounded enclosure, rational receives the
// exact value, double receives the upward point value.
// ---------------------------------------------------------------------------
template <int D, typename T>
T cr_interpolation_constant() {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::cr_interpolation_constant: "
                  "D must be 2 or 3");
    return constant_from_string<T>::get(
        D == 2 ? cr_interpolation_constant_2d_str
               : cr_interpolation_constant_3d_str);
}

// ---------------------------------------------------------------------------
// cr_projection_error_constant_sq: C_h^2 = c_D^2 * h^2 for the CR space V,
// with h^2 = V.max_edge_length_sq() (the exact-for-every-T master API of
// cr1_space, ruling R4 there).  The SQUARED form is deliberate (ruling R17):
// no square root appears on this path, so an exact rational T goes through
// exactly, and the lower bound formula below consumes C_h^2 directly.
//
// NO h1_conforming static_assert here, deliberately: the static_assert
// convention of this layer (see the file head) guards functions that ASSUME
// V_h in H^1_0, and this function is the opposite case -- it exists precisely
// because the CR space is NONCONFORMING (space_traits h1_conforming == false),
// so the convention of VER-0 design section two does not apply to it.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
T cr_projection_error_constant_sq(const vcp::bfem::cr1_space<D, T, P, SP>& V) {
    const T c = cr_interpolation_constant<D, T>();
    return c * c * V.max_edge_length_sq();
}

// ---------------------------------------------------------------------------
// eigenvalue_lower_bound: lambda_k >= lambda_{h,k} / (1 + C_h^2 lambda_{h,k})
// (Liu 2015 Thm. 2.1).  lambda_h is the k-th DISCRETE nonconforming
// eigenvalue, ch_sq is C_h^2 from cr_projection_error_constant_sq.  The
// formula is evaluated in whatever scalar T the two arguments carry: a kv
// interval T yields an enclosure of the bound, rational yields it exactly,
// a point T yields the point evaluation.
// ---------------------------------------------------------------------------
template <typename T>
T eigenvalue_lower_bound(const T& lambda_h, const T& ch_sq) {
    return lambda_h / (T(1) + ch_sq * lambda_h);
}

// ---------------------------------------------------------------------------
// poincare_constant_h01_sq_bound: C_P^2 <= 1 / lambda_1^{lower}, with
// lambda_1^{lower} a POSITIVE lower bound of the first Dirichlet eigenvalue
// (for instance eigenvalue_lower_bound applied to lambda_{h,1}).  Square-root
// free, so rational T goes through exactly.
//
// poincare_constant_h01_bound: the convenience square root C_P <= 1 / sqrt(...),
// for floating point / interval T.  Same lazy policy as cr1_space's
// max_edge_length(): an ordinary non-virtual template, instantiated only when
// called, so a rational T that has no square root still compiles as long as
// only the _sq form is used.
// ---------------------------------------------------------------------------
template <typename T>
T poincare_constant_h01_sq_bound(const T& lambda1_lower) {
    return T(1) / lambda1_lower;
}

template <typename T>
T poincare_constant_h01_bound(const T& lambda1_lower) {
    using std::sqrt;
    return sqrt(T(1) / lambda1_lower);
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#include <vcp/bfem/constants/detail/poisson_dict_impl.hpp>
#endif // VCP_BFEM_CONSTANTS_POISSON_CONSTANTS_HPP
