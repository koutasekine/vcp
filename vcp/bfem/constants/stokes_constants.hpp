// vcp/bfem/constants/stokes_constants.hpp
//
// Verified constants of the Stokes problem -- the VECTOR-FIELD ENTRY of the
// constants layer (CONST-E).  Every entry point below is a THIN WRAPPER of
// the corresponding SCALAR entry point of poisson_constants.hpp / its detail
// layer / sobolev_constants.hpp: it forwards its arguments unchanged and
// returns the scalar result unchanged.  No arithmetic, no coefficient, no
// rounding and no numeric literal is added anywhere on this file (gate
// G-H1 of CONST-E nails the value AND every report field down bit for bit,
// so a later coefficient smuggled into a wrapper turns that gate red).
//
// What a wrapper buys, then, is not a number but a THEOREM: the statement
// that the scalar constant is a valid constant for the vector-valued space
// as well.  The three lemmas that carry those statements are section one of
// the CONST-E design and are transcribed here in full.
//
// SCOPE OF THE WRAPPER STATEMENTS ABOVE (CONST-G, owner ruling of
// 2026-09-07): the "thin wrapper, no arithmetic, no coefficient" statements
// describe the FOUR CONST-E entry points (the two L^2 projection entries,
// the Sobolev embedding entry and stokes_poincare_constant_h01_sq_bound),
// and they keep holding for those four bit for bit (gate G-H1).  CONST-G
// adds, at the END of this file, the first entries of the Stokes side that
// COMPUTE something: the verified lower bounds of the Stokes eigenvalues
// (stokes_eigenvalue_lower_bounds_h01) and the div-free Poincare constant
// built on them (stokes_divfree_poincare_constant_h01_sq_bound).  Those two
// consume the vector CR space of vcp/bfem/cr1/vcr1_space.hpp, the exact
// divergence-free null space of vcp/bfem/cr1/vcr1_divfree.hpp and the
// EXISTING lower bound formula eigenvalue_lower_bound of
// poisson_constants.hpp; they add no numeric literal and no new formula of
// their own (the lexical regime of this file, gate G-H4, is unchanged: no
// decimal, no `double`, no `float` in the code).
//
// ---------------------------------------------------------------------------
// BUNDLING CONVENTION (CONST-E design section one).  For a vector field
// u = (u_1, ..., u_D) in (H^1_0(Omega))^D the norms of this file are the
// COMPONENTWISE l^2 bundles
//
//     ||u||^2 := sum_i ||u_i||^2,        |u|_1^2 := sum_i ||grad u_i||^2.
//
// Every statement below is made with respect to THESE norms; a different
// bundling (for instance a max over components) would change the constants.
//
// LEMMA 1 (L^2 projection, factor one).  Let Pi_h u := (pi_h u_i)_i be the
// componentwise L^2 projection.  Then, summing the scalar inequality
// ||u_i - pi_h u_i||^2 <= C_{0,h}^2 |u_i|_1^2 over i = 1 .. D,
//
//     ||u - Pi_h u||^2 = sum_i ||u_i - pi_h u_i||^2
//                     <= C_{0,h}^2 sum_i |u_i|_1^2 = C_{0,h}^2 |u|_1^2,
//
// i.e. THE SCALAR CONSTANT CARRIES OVER WITH FACTOR ONE -- the derivation is
// the sum of the scalar inequalities and nothing else.  The same argument
// applies element by element, which is why the per-element wrapper is
// legitimate too.  [self-contained, elementary]
//
// LEMMA 2 (Sobolev embedding, p >= 2 ONLY).  Let ||u(x)|| denote the
// pointwise Euclidean norm of the vector field.  Minkowski's integral
// inequality in the form  ||  ||u||  ||_p <= ( sum_i ||u_i||_p^2 )^{1/2}
// REQUIRES p >= 2 (for p < 2 the inequality points the other way and the
// argument collapses); with p >= 2 and the scalar bound ||u_i||_p <= C_p
// |u_i|_1 applied componentwise,
//
//     ||  ||u||  ||_p <= ( sum_i ||u_i||_p^2 )^{1/2}
//                     <= C_p ( sum_i |u_i|_1^2 )^{1/2} = C_p |u|_1,
//
// and the sigma-weighted variant of the scalar route follows the same
// pattern.  p < 2 IS NOT COVERED by this header.  This costs nothing in
// practice: the scalar entry point's own domain is p in (2, infinity) for
// n = 2 and p in [2, 2n/(n-2)] for n >= 3, both inside p >= 2, so an
// argument that would break Lemma 2 is already refused by the scalar side
// (which is where the throw comes from -- this file adds no check of its
// own).  [elementary and standard; the SOURCE HAS NOT BEEN VERIFIED
// VERBATIM on this track, unlike the scalar Sobolev constants themselves,
// whose Theorem A.1 / Corollary A.2 / Theorem A.3 were checked against the
// original in CONST-D]
//
// LEMMA 3 (Poincare on the divergence-free subspace -- a SUBSTITUTE, not a
// sharp constant).  The best constant of a Friedrichs inequality is a
// supremum over the space; restricting the space to the divergence-free
// subspace V_sigma of (H^1_0)^D can only DECREASE that supremum, so the
// scalar (full-space) Poincare constant is a VALID UPPER BOUND on
// V_sigma.  It is, however, NOT SHARP: the first divergence-free (Stokes)
// eigenvalue is in general strictly larger than the first Dirichlet
// Laplacian eigenvalue, so the true div-free constant is strictly smaller.
// Sharpening it is a future track (see the planned contents below); until
// then this wrapper is the honest, conservative substitute.  [self-contained]
// ---------------------------------------------------------------------------
//
// STILL PLANNED (future tracks; VER-0 design, section 2 -- these names are
// kept here deliberately as the forward reference of this file, and CONST-E
// gives NONE of them a body):
//
//     stokes_ritz_projection_error_constant   the Stokes counterpart of the
//                                             Ritz projection error constant
//                                             (GIVEN A BODY by CONST-H in
//                                             the SQUARED, dictionary served
//                                             form, see the end of this file)
//     a div-free Poincare constant            (the SHARP one of Lemma 3 --
//                                             the wrapper below is only the
//                                             full-space substitute; GIVEN A
//                                             BODY by CONST-G, see
//                                             stokes_divfree_poincare_constant_h01_sq_bound
//                                             at the end of this file)
//     the hypercircle kappa_h of the mixed    Stokes formulation (GIVEN A
//                                             BODY by CONST-H:
//                                             stokes_hypercircle_kappa_h01_squared)
//     C_{4,P}                                 (appears on the Stokes side
//                                             when it is needed; R10)
//
// The remaining ones are still to be given bodies only after the vcr1 and
// the local div-free constraint tracks are complete (the vcr1 space and the
// exact div-free constraint now exist: CONST-G); the C_h of the Stokes
// projection itself remains out of scope (CONST-H).
//
// VER-1 boundary-condition ledger (owner ruling R23; same format as the one
// at the head of poisson_constants.hpp -- the BC assumption of an entry
// point is part of its NAME, the BC tag sitting right after the concept
// name and before the form suffixes):
//
//     entry point                                    BC assumption
//     stokes_l2_projection_element_constants_sq      BC-free: the
//     stokes_l2_projection_error_constant_sq         componentwise L^2
//                                                    projection of Lemma 1
//                                                    carries no boundary
//                                                    condition -- valid for
//                                                    Dirichlet, Neumann and
//                                                    Robin velocity alike
//                                                    (the scalar entries it
//                                                    forwards to are BC-free
//                                                    for the same reason)
//     stokes_sobolev_embedding_constant_h01          H^1_0 VECTOR fields:
//                                                    (H^1_0(Omega))^D, the
//                                                    componentwise
//                                                    homogeneous Dirichlet
//                                                    velocity space of
//                                                    Lemma 2
//     stokes_poincare_constant_h01_sq_bound          H^1_0 VECTOR fields:
//                                                    the supplied lambda_1
//                                                    is a DIRICHLET
//                                                    eigenvalue and the
//                                                    inequality is the
//                                                    Friedrichs form on
//                                                    (H^1_0)^D (Lemma 3;
//                                                    on the div-free
//                                                    subspace it stays
//                                                    valid but not sharp)
//     stokes_eigenvalue_lower_bounds_h01             H^1_0 VECTOR fields:
//                                                    the Stokes eigenvalue
//                                                    problem (5) of
//                                                    Xie-Xie-Liu 2018 on
//                                                    V_0 = div-free
//                                                    (H^1_0)^2, CR space
//                                                    with vanishing
//                                                    boundary facet means
//                                                    (CONST-G)
//     stokes_divfree_poincare_constant_h01_sq_bound  H^1_0 VECTOR fields,
//                                                    DIV-FREE subspace: the
//                                                    supplied lambda_1 is a
//                                                    verified lower bound of
//                                                    the first STOKES
//                                                    eigenvalue (CONST-G)
//     stokes_hypercircle_kappa_h01_squared           H^1_0 VECTOR fields:
//                                                    the SV velocity space
//                                                    carries the homogeneous
//                                                    Dirichlet condition
//                                                    (no-slip), the flux
//                                                    RT_h carries none
//                                                    (CONST-H)
//     stokes_ritz_projection_error_constant_h01_sq   H^1_0 VECTOR fields:
//                                                    the Stokes projection
//                                                    P_h of (11) onto the
//                                                    no-slip SV space; C_0,h
//                                                    is the BC-free L^2
//                                                    projection constant
//                                                    (CONST-H)
//
// The planned entries above are likewise all framed on H^1_0 vector fields
// -- the div-free test space and the mixed hypercircle assume homogeneous
// Dirichlet velocity -- so, when they are given bodies, their names carry
// the BC tag of R23 as well.
//
// UMBRELLA (the Stokes side of the owner ruling of 2026-08-20, "one entry
// point"): this header includes poisson_constants.hpp, which is itself the
// umbrella over the dictionary detail layer and over sobolev_constants.hpp.
// So a translation unit that includes THIS file alone can call the scalar
// Poisson, dictionary and Sobolev entry points as well as the Stokes ones.
// No cycle: poisson_constants.hpp reaches no Stokes header.
//
// Lexical policy (inherited from poisson_constants.hpp, design 6.1): no
// decimal literals, no `double` / `float` tokens in code; sqrt and
// kv::constants<T>::pi() are allowed.  This file uses none of the three.
//
// static_assert convention: as in poisson_constants.hpp -- whoever adds to
// this layer a function that takes a SPACE as an argument and assumes a
// conforming space must place the corresponding
// static_assert(vcp::bfem::space_traits<S>::..., ...) in it.  The entry
// points below take a mesh or plain scalars, so there is nothing to assert.
//
// Authority: sandbox/docs/design/CONST-E_design_v1.0.md and, for the file's
// place in the layer, sandbox/docs/design/VER-0_design_v1.0.md.

#ifndef VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP
#define VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <stdexcept>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/cr1/vcr1_space.hpp>
#include <vcp/bfem/cr1/vcr1_divfree.hpp>
#include <vcp/bfem/constants/poisson_constants.hpp>
// CONST-G: the exact rationalization of point-interval vertices
// (detail::rational_of_point_interval) lives in element_projection.hpp,
// which itself includes poisson_constants.hpp (no cycle: neither reaches
// this file).
#include <vcp/bfem/constants/element_projection.hpp>
// CONST-H: the Stokes hypercircle core (SV velocity on the barycentric
// refinement, RT flux with the multiplier constraint, the prolongation
// between the two meshes).  Includes the sv/ and rt/ layers; no cycle
// (none of them reaches this file).
#include <vcp/bfem/constants/detail/stokes_hypercircle_impl.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// stokes_l2_projection_element_constants_sq(Th, d, rep): the per-element
// vector of squared L^2 projection error constants for VECTOR fields, i.e.
// exactly l2_projection_element_constants_sq(Th, d, rep) forwarded.
//
// Lemma 1 is what makes the forwarding legitimate: with the componentwise
// projection Pi_h and the l^2 bundling, the element-local scalar inequality
// summed over the D components gives the vector inequality WITH THE SAME
// CONSTANT (factor one).  Entry e is therefore an upper bound of the
// element-local vector constant C_d(K_e)^2 as well.
//
// Everything else -- the dictionary resolution order, the report census,
// the maximum rule, the exceptions (empty mesh, unservable element,
// negative d, non-dyadic vertex, degenerate element) -- is the scalar
// entry point's, unchanged, because this IS the scalar entry point.
//
// D is deduced from the mesh; the two dimensional and three dimensional
// scalar overloads are picked by ordinary overload resolution at
// instantiation.
// ---------------------------------------------------------------------------
template <int D, typename T>
std::vector<T> stokes_l2_projection_element_constants_sq(
        const vcp::bfem::mesh<D, T>& Th, int d,
        mesh_resolution_report* rep = nullptr) {
    return l2_projection_element_constants_sq(Th, d, rep);
}

// ---------------------------------------------------------------------------
// stokes_l2_projection_error_constant_sq(Th, d, rep): the mesh-level squared
// L^2 projection error constant for VECTOR fields,
//
//     ||u - Pi_h u||^2 <= C_{0,h}(d)^2 |u|_1^2,   u in (H^1(Omega))^D,
//
// i.e. exactly l2_projection_error_constant_sq(Th, d, rep) forwarded.  Same
// justification as above (Lemma 1: the sum over components of the scalar
// inequality, factor one), same report, same maximum rule, same exceptions.
// SQUARED deliberately (R17): no square root on this path, so an exact
// fraction scalar T closes end to end here too.
// ---------------------------------------------------------------------------
template <int D, typename T>
T stokes_l2_projection_error_constant_sq(const vcp::bfem::mesh<D, T>& Th,
                                         int d,
                                         mesh_resolution_report* rep = nullptr) {
    return l2_projection_error_constant_sq(Th, d, rep);
}

// ---------------------------------------------------------------------------
// stokes_sobolev_embedding_constant_h01<D, T>(pnum, pden, mea, rho_lb,
// sigma, which): the constant of the embedding (H^1_0(Omega))^D -> L^p in
// the sense
//
//     ||  ||u||  ||_p <= C_p |u|_1,      p >= 2,
//
// with ||u(x)|| the pointwise Euclidean norm -- exactly
// sobolev_embedding_constant_h01<D, T>(pnum, pden, mea, rho_lb, sigma,
// which) forwarded, including the min selection between the measure route
// (A.2) and the spectrum route (A.3) and the `which` output (zero for
// (A.2), one for (A.3)).
//
// Lemma 2 is what makes the forwarding legitimate, and it is the ONE place
// on this file where the vector statement is weaker than the scalar one:
// the Minkowski step needs p >= 2.  No check is added here -- the scalar
// domain (n = 2: p in (2, infinity); n >= 3: p in [2, 2n/(n-2)]) already
// lies inside p >= 2, so out-of-domain input throws from the scalar route
// with the scalar message, bit-identically to a direct scalar call.
//
// D is NOT deducible from the arguments and must be given explicitly, as on
// the scalar side.
// ---------------------------------------------------------------------------
template <int D, typename T>
T stokes_sobolev_embedding_constant_h01(long long pnum, long long pden,
                                        const T& mea, const T& rho_lb,
                                        const T& sigma,
                                        int* which = nullptr) {
    return sobolev_embedding_constant_h01<D, T>(pnum, pden, mea, rho_lb,
                                                sigma, which);
}

// ---------------------------------------------------------------------------
// stokes_poincare_constant_h01_sq_bound(lambda1_lower): the squared
// Poincare (Friedrichs) constant for VECTOR fields,
//
//     ||u||^2 <= C_P^2 |u|_1^2,      u in (H^1_0(Omega))^D,
//
// i.e. exactly poincare_constant_h01_sq_bound(lambda1_lower) forwarded, the
// supplied argument being a verified LOWER bound of the first DIRICHLET
// eigenvalue of the scalar Laplacian.
//
// Two things are worth being explicit about.  (a) On the whole space
// (H^1_0)^D the constant is the scalar one with factor one -- Lemma 1's
// summation argument applied to the Friedrichs inequality.  (b) On the
// DIVERGENCE-FREE subspace V_sigma it remains a valid upper bound (Lemma 3:
// the best constant is a supremum, and restricting the space can only
// decrease it) but it is NOT SHARP, because the first Stokes eigenvalue is
// in general strictly larger than lambda_1.  A user who needs the sharp
// div-free constant is waiting for the future track named at the head of
// this file; what this entry point promises is a guaranteed bound, not the
// best one.
//
// Square-root free, so a rational T goes through exactly.  There is
// deliberately no stokes_poincare_constant_h01_bound: the square root form
// adds nothing that the scalar poincare_constant_h01_bound does not already
// give, and CONST-E keeps the vector surface to the four entry points of
// its design section two.
// ---------------------------------------------------------------------------
template <typename T>
T stokes_poincare_constant_h01_sq_bound(const T& lambda1_lower) {
    return poincare_constant_h01_sq_bound(lambda1_lower);
}

// ===========================================================================
// CONST-G: verified lower bounds of the Stokes eigenvalues and the div-free
// Poincare constant (the SHARP constant of Lemma 3).
//
// Primary source (verbatim-verified, CONST-G P0, correspondence table in
// sandbox/docs/reports/CONST-G_P0_reconnaissance_report.md):
//     M. Xie, H. Xie, X. Liu, "Explicit lower bounds for Stokes eigenvalue
//     problems by using nonconforming finite elements", Japan J. Indust.
//     Appl. Math. 35 (2018) 335-354.
//
//   (5)   continuous problem on V_0 = div-free (H^1_0(Omega))^2:
//         a(u, v) = lambda r(u, v) for all v in V_0, eigenvalues
//         0 < lambda_1 <= lambda_2 <= ...;  the div-free Poincare constant
//         is C_P^sigma^2 = 1 / lambda_1 (the Friedrichs inequality on V_0
//         with its best constant, i.e. Lemma 3 made sharp).
//   (9)   discrete space V_{0,h}^(1) = { v_h in (CR)^2 : b_h(v_h, q_h) = 0
//         for all piecewise constant q_h } -- the element-mean divergence
//         constraint (vcr1_divfree.hpp, section (i): with vanishing
//         boundary facet means this is "int_K div v_h = 0 on every K", by
//         the discrete divergence theorem).
//   Prop. 1 (item 1)  the CR interpolation maps V_0 into V_{0,h}^(1)
//         (Fortin property; the proof uses only grad q_h|_K = 0 and the
//         preservation of facet means).
//   Thm. 4 / (38)     lambda_i >= lambda_{i,h} / (1 + (alpha_1 h)^2
//         lambda_{i,h}), 1 <= i <= dim V_{0,h}^(1), alpha_1 = 0.1893, h =
//         the largest element diameter (longest edge).
//
// alpha_1 is the EXISTING cr_interpolation_constant_2d_str of
// poisson_constants.hpp (Liu 2015, the same constant the scalar CR lower
// bound uses -- Lemma 1 of the paper cites it) and the bound is the EXISTING
// eigenvalue_lower_bound(lambda_h, C_h^2) with C_h^2 =
// cr_projection_error_constant_sq(V) = alpha_1^2 h^2 of the scalar CR space:
// no new decimal string and no new formula enter this file.
//
// The paper's sharper constant alpha~_1 = 0.1761 of its equation (40) holds
// ONLY for uniform meshes of isosceles right triangles and is what its
// Tables 1 / 3 are computed with; it is NOT a constant of this library (the
// CONST-G gates keep it test-side, as a literature value).  The entry points
// below use alpha_1 for every admissible mesh.
//
// Two dimensions only: the primary source states Theorem 4 for triangular
// partitions of Omega in R^2.  The three dimensional analogue would be a
// self-derived extension (dimension-free proof steps plus the 3D CR
// constant) and is deferred by the CONST-G ruling until its mathematics is
// established; D == 3 is rejected at compile time.
// ===========================================================================

// ---------------------------------------------------------------------------
// stokes_eigenvalue_report<T>: the diagnostics of one call of
// stokes_eigenvalue_lower_bounds_h01 -- the sizes, the null space
// certificate (vcr1_divfree.hpp, section (ii)), C_h^2 and the discrete
// eigenvalue enclosures the bounds were made from.
// ---------------------------------------------------------------------------
template <typename T>
struct stokes_eigenvalue_report {
    int n_elements;
    int n_scalar_dofs;            // facets of the scalar CR space
    int n_reduced_scalar;         // interior facets n_r
    int columns;                  // 2 n_r, the reduced vector dofs
    int rank;                     // exact rank of the constraint matrix
    int rank_deficiency;          // n_elements - rank (1 on a connected domain)
    int dim;                      // dim V_{0,h}^(1) after the constraint
    long long pairs_checked;      // B z == 0 products verified (all of them)
    T ch_sq;                      // alpha_1^2 h^2
    std::vector<T> lambda_h;      // discrete eigenvalue enclosures, first nev
};

// ---------------------------------------------------------------------------
// stokes_eigenvalue_lower_bounds_h01<D, T, DP, SP>(Th, nev, rep): verified
// lower bounds of the first nev Stokes eigenvalues on the mesh Th, i.e. the
// vector ( lambda_{k,h} / (1 + C_h^2 lambda_{k,h}) )_{k = 1..nev} with every
// lambda_{k,h} a verified enclosure of the k-th eigenvalue of the reduced
// pencil Z^T A Z x = lambda Z^T M Z x (A, M the block-diagonal broken
// stiffness / mass of the vector CR space after the homogeneous Dirichlet
// reduction, Z the certified exact null space basis of the element-mean
// divergence constraint) and C_h^2 = cr_projection_error_constant_sq of the
// scalar CR space.  Entry k of the result is a T whose LOWER end is the
// guaranteed lower bound of lambda_{k+1}; the upper end carries no claim
// beyond enclosing the formula's value.
//
// Pipeline: vcr1_space (block-diagonal forms) -> dirichlet_reduction of the
// scalar space (the reduced vector numbering is d * n_r + r) -> exact
// rationalization of the point-interval vertices -> vcr1_divfree_basis
// (constraint, null space, certificate) -> reduce_pencil (interval,
// symmetric by mirroring) -> eigsymge (Miyajima, all eigenvalues enclosed)
// -> eigenvalue_lower_bound, k = 1 .. nev.
//
// T: a kv interval scalar (T::base_type required by the policy defaults, as
// in element_projection.hpp).  Exceptions: std::invalid_argument (nev < 1
// or nev > dim V_{0,h}^(1) -- Theorem 4 is stated for 1 <= i <= dim only;
// empty mesh), vcp::invalid_argument (a vertex that is not a point interval
// or not exactly rationalizable, from rational_of_point_interval),
// vcp::verification_error (null space certificate failed, or eigsymge could
// not certify the enclosures), vcp::bfem::degenerate_element (flat element).
// ---------------------------------------------------------------------------
template <int D, typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
std::vector<T> stokes_eigenvalue_lower_bounds_h01(
        const vcp::bfem::mesh<D, T>& Th, int nev,
        stokes_eigenvalue_report<T>* rep = nullptr) {
    static_assert(D == 2,
                  "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01: "
                  "two dimensions only -- the primary source (Xie-Xie-Liu "
                  "2018, Theorem 4) is stated for triangulations of a plane "
                  "domain; the 3D extension is deferred by the CONST-G ruling");
    typedef ::vcp::bfem::detail::rational exact_t;
    if (Th.num_elements() < 1)
        throw std::invalid_argument(
            "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01: empty mesh");
    if (nev < 1)
        throw std::invalid_argument(
            "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01: nev < 1");

    // the vector CR space and the scalar Dirichlet reduction
    vcp::bfem::vcr1_space<D, T, DP, SP> V(Th);
    vcp::bfem::dirichlet_reduction<T, DP, SP> red(V.scalar_ndof(),
                                                  V.scalar().boundary_dofs());
    std::vector<int> s2r(static_cast<std::size_t>(V.scalar_ndof()));
    for (int i = 0; i < V.scalar_ndof(); ++i)
        s2r[static_cast<std::size_t>(i)] = red.to_reduced(i);

    // exact rationalization of the mesh (point-interval vertices only)
    std::vector<std::array<exact_t, D> > xv;
    xv.reserve(static_cast<std::size_t>(Th.num_vertices()));
    for (int v = 0; v < Th.num_vertices(); ++v) {
        std::array<exact_t, D> q;
        for (int d = 0; d < D; ++d)
            q[static_cast<std::size_t>(d)] = detail::rational_of_point_interval<T>(
                Th.vertex(v)[static_cast<std::size_t>(d)]);
        xv.push_back(q);
    }
    std::vector<std::array<int, D + 1> > xe;
    xe.reserve(static_cast<std::size_t>(Th.num_elements()));
    for (int e = 0; e < Th.num_elements(); ++e) xe.push_back(Th.element(e));
    const vcp::bfem::mesh<D, exact_t> Tx =
        vcp::bfem::mesh<D, exact_t>::from_lists(xv, xe);

    // constraint, certified null space, reduced pencil
    const vcp::bfem::vcr1_divfree_basis<D, exact_t> Z(Tx, s2r);
    const int dim = Z.dim();
    if (nev > dim)
        throw std::invalid_argument(
            "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01: nev exceeds "
            "dim V_{0,h} (Theorem 4 covers 1 <= i <= dim only)");
    vcp::matrix<T, DP> A = Z.template reduce_pencil<T, DP, SP>(
        red.reduce(V.scalar().stiffness()));
    vcp::matrix<T, DP> M = Z.template reduce_pencil<T, DP, SP>(
        red.reduce(V.scalar().mass()));
    vcp::matrix<T, DP> E;
    eigsymge(A, M, E);                   // verified enclosures, ascending

    // the existing lower bound formula with the existing CR constant
    const T ch_sq = cr_projection_error_constant_sq(V.scalar());
    std::vector<T> out;
    out.reserve(static_cast<std::size_t>(nev));
    std::vector<T> lam;
    lam.reserve(static_cast<std::size_t>(nev));
    for (int k = 0; k < nev; ++k) {
        const T lk = (E.columnsize() == 1) ? E(k, 0) : E(k, k);
        lam.push_back(lk);
        out.push_back(eigenvalue_lower_bound(lk, ch_sq));
    }
    if (rep != nullptr) {
        const vcp::bfem::vcr1_divfree_certificate& c = Z.certificate();
        rep->n_elements = Th.num_elements();
        rep->n_scalar_dofs = V.scalar_ndof();
        rep->n_reduced_scalar = red.reduced_size();
        rep->columns = c.columns;
        rep->rank = c.rank;
        rep->rank_deficiency = c.rank_deficiency;
        rep->dim = dim;
        rep->pairs_checked = c.pairs_checked;
        rep->ch_sq = ch_sq;
        rep->lambda_h = lam;
    }
    return out;
}

// ---------------------------------------------------------------------------
// stokes_divfree_poincare_constant_h01_sq_bound(stokes_lambda1_lower): the
// squared Poincare (Friedrichs) constant on the DIVERGENCE-FREE subspace
//
//     ||u||^2 <= C_P^sigma^2 |u|_1^2,      u in V_sigma = div-free (H^1_0)^2,
//
// namely C_P^sigma^2 <= 1 / lambda_1^S,lower, with the argument a verified
// LOWER bound of the first STOKES eigenvalue (entry 0 of
// stokes_eigenvalue_lower_bounds_h01, or any other certified bound of
// lambda_1 of (5)).  This is the SHARP constant Lemma 3 at the head of this
// file could not give: on V_sigma the best Friedrichs constant is exactly
// 1 / lambda_1^S (the min-max principle (6) of the paper), so the bound is
// tight up to the discretization gap of lambda_1^S,lower.  Square-root free,
// hence exact on a rational T (R17); no sqrt form is added (same reasoning
// as for the transferred entry).
//
// The transferred entry stokes_poincare_constant_h01_sq_bound(lambda_1)
// (argument: a Dirichlet LAPLACIAN eigenvalue bound, Lemma 3) is KEPT as a
// separate entry point, and the two are deliberately NOT combined into a
// min: the two arguments are numbers of different theorems (a Laplacian
// bound versus a Stokes bound), and an automatic min would hide which
// theorem the returned number rests on.  A caller who holds both bounds may
// take the min explicitly, knowing what each one is.
// ---------------------------------------------------------------------------
template <typename T>
T stokes_divfree_poincare_constant_h01_sq_bound(const T& stokes_lambda1_lower) {
    return T(1) / stokes_lambda1_lower;
}

// ===========================================================================
// CONST-H: the Stokes projection error constant C_h^S of Theorem 7 of
// Liu-Nakao-Oishi, Commun. Nonlinear Sci. Numer. Simul. 108 (2022) 106223,
//
//     ||grad(u - P_h u)|| <= C_h ||f||,  ||u - P_h u|| <= C_h ||grad(u - P_h u)||,
//     C_h := sqrt( kappa_h^2 + C_{0,h}^2 ),
//
// for the Stokes projection P_h (11) onto the Scott-Vogelius space on the
// barycentric refinement of the mesh, with kappa_h the hypercircle quantity
// of the paper (its section 3-1; a verified generalized eigenvalue: no literature constant
// enters it -- detail/stokes_hypercircle_impl.hpp, whose head carries the
// derivation record) and C_{0,h} the L^2 projection error constant of pi_h
// onto X_h = broken P^{k-1} on the COARSE mesh (paper (9), bounded by (8)).
//
// Two things the paper's example section makes explicit and this layer
// mirrors: the two meshes are independent (kappa lives on the pair T_h /
// T_SV, C_{0,h} on T_h alone), and d = m = k - 1 (X_h and RT_h of degree
// k - 1 for the velocity degree k).  Ownership of C_{0,h}: by DEFAULT the
// dictionary of this layer serves it (stokes_l2_projection_error_constant_sq,
// i.e. the CONST-C2 dictionary composition, Lemma 1 factor one); the
// paper itself uses the closed form of its (31) for its
// cube-based mesh family (Liu, JCAM 371 (2020) 112666, Table 3), and the
// owner ruling of CONST-H gate 1 (item 2) added the EXPLICIT overload below
// for a caller who holds such a mesh-family specific certified value: the
// value is the caller's responsibility and the report says nothing about it.
//
// Three dimensions only: the SV construction class this layer is gated on
// is "barycentric (Alfeld) split, k >= 3" in 3D; D != 3 is rejected at
// compile time (the core itself is dimension generic).
// ===========================================================================

// ---------------------------------------------------------------------------
// stokes_hypercircle_kappa_h01_squared<T, DP, SP>(Th, k, rep): kappa_h^2 of
// the Stokes hypercircle for the coarse mesh Th and the velocity degree k
// (SV velocity: vector P^k on alfeld_refine(Th) with no-slip; flux: RT_{k-1}
// on Th; f_h and the multiplier: broken P^{k-1} on Th; phi_h: P^k on Th
// modulo constants).  Returned SQUARED, as the point interval of the largest
// upper end of the certified spectrum (the convention of
// hypercircle_kappa_h01_squared); rep, when non-null, receives the sizes,
// the certificate counts and every eigenvalue enclosure.
//
// Exceptions: as detail::stokes_hypercircle_kappa_core.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T stokes_hypercircle_kappa_h01_squared(const vcp::bfem::mesh<3, T>& Th, int k,
                                       stokes_hypercircle_report<T>* rep = nullptr) {
    stokes_hypercircle_report<T> local;
    stokes_hypercircle_report<T>& r = (rep != nullptr) ? *rep : local;
    detail::stokes_hypercircle_kappa_core<3, T, DP, SP>(Th, k, r,
        static_cast<detail::stokes_hypercircle_materials<T, DP>*>(nullptr));
    return r.kappa2;
}

// ---------------------------------------------------------------------------
// stokes_ritz_projection_error_constant_h01_sq<T, DP, SP>(Th, k, rep, drep):
//
//     C_h(k)^2 = stokes_hypercircle_kappa_h01_squared(Th, k)
//                + stokes_l2_projection_error_constant_sq(Th, k - 1)
//
// the Stokes twin of the three dimensional ritz_projection_error_constant_h01_sq
// of the Poisson side: the SAME expression in the SAME order (kappa first,
// then the dictionary served C_{0,h}^2 of the coarse mesh with d = k - 1 --
// the only k -> d conversion point of this path, pinned by gate G-T4), the
// sum returned SQUARED (R17: no square root, an exact fraction closes).
// rep / drep, when non-null, receive the kappa report and the dictionary
// resolution census.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T stokes_ritz_projection_error_constant_h01_sq(const vcp::bfem::mesh<3, T>& Th, int k,
                                               stokes_hypercircle_report<T>* rep = nullptr,
                                               mesh_resolution_report* drep = nullptr) {
    const T kappa2 = stokes_hypercircle_kappa_h01_squared<T, DP, SP>(Th, k, rep);
    const T c0h_sq = stokes_l2_projection_error_constant_sq(Th, k - 1, drep);
    return kappa2 + c0h_sq;
}

// ---------------------------------------------------------------------------
// stokes_ritz_projection_error_constant_h01_sq<T, DP, SP>(Th, k, c0h_sq, rep):
// the EXPLICIT C_{0,h}^2 overload (owner ruling, CONST-H gate 1 item 2).
// The caller supplies a certified upper bound of C_{0,h}^2 of (9) for ITS
// mesh family (for instance the square of the closed form (31) of the paper for the
// cube / five-tetrahedra family, or a sharper dictionary of its own); this
// function composes it with the same kappa in the same order and vouches
// for nothing about the supplied value.  The default entry point above is
// the dictionary served one.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T stokes_ritz_projection_error_constant_h01_sq(const vcp::bfem::mesh<3, T>& Th, int k,
                                               const T& c0h_sq,
                                               stokes_hypercircle_report<T>* rep = nullptr) {
    const T kappa2 = stokes_hypercircle_kappa_h01_squared<T, DP, SP>(Th, k, rep);
    return kappa2 + c0h_sq;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP
