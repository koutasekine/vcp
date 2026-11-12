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
// CONST-J adds four more THIN WRAPPERS in the CONST-E sense, each right after
// the entry point it belongs to: the square root forms
// stokes_poincare_constant_h01_bound and
// stokes_divfree_poincare_constant_h01_bound, and the two Sobolev routes
// stokes_sobolev_embedding_constant_h01_measure / _spectrum.  They forward to
// the scalar entry points of the same name (the div-free one takes the square
// root of the squared entry above it) and add no formula, no coefficient and
// no check.  CONST-J also adds the squared constant set
// stokes_ritz_projection_constants_h01_sq after the explicit composition at
// the end of this file, built from the same calls as the compositions.
//
// CONST-K extends stokes_eigenvalue_lower_bounds_h01 to D = 3 (source: Liu,
// Nakao, Oishi 2022, section 3-2, cited in the CONST-G section below) and
// routes both dimensions through detail::stokes_divfree_eigen_enclosures
// (dense now; the replacement point of the planned sparse path) and the
// common index-certifying function of poisson_constants.hpp.  The report
// stokes_eigenvalue_report gains four trailing fields with the names of
// eigenvalue_report: index_lo / index_hi (k..k on the dense path) and
// nnz_stiffness / nnz_mass (0 on the dense path).  The policy defaults of
// the entry follow the eigenvalue backend guideline of spimats.hpp (ruling
// M-8).  The entry accepts interval vertices: the divergence-free basis is
// an interval enclosure built on the given mesh (supplementary ruling 2), so
// the point-interval / rationalization requirement of CONST-G is gone; the
// exact null space of vcr1_divfree.hpp stays in the library, unused by the
// entry.  No new decimal string, no new formula.
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
//     stokes_sobolev_embedding_constant_h01_measure  H^1_0 VECTOR fields:
//     stokes_sobolev_embedding_constant_h01_spectrum the two routes of the
//                                                    entry above, same frame
//                                                    (Lemma 2; CONST-J)
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
//     stokes_poincare_constant_h01_bound             H^1_0 VECTOR fields:
//                                                    same (square root form;
//                                                    CONST-J)
//     stokes_eigenvalue_lower_bounds_h01             H^1_0 VECTOR fields:
//                                                    the Stokes eigenvalue
//                                                    problem (5) of
//                                                    Xie-Xie-Liu 2018 on
//                                                    V_0 = div-free
//                                                    (H^1_0)^2, CR space
//                                                    with vanishing
//                                                    boundary facet means
//                                                    (CONST-G)
//     stokes_eigenvalue_lower_bounds_h01 (D = 3)     H^1_0 VECTOR fields:
//                                                    problem (21) of
//                                                    Liu-Nakao-Oishi 2022,
//                                                    no-slip (CONST-K)
//     stokes_divfree_poincare_constant_h01_sq_bound  H^1_0 VECTOR fields,
//                                                    DIV-FREE subspace: the
//                                                    supplied lambda_1 is a
//                                                    verified lower bound of
//                                                    the first STOKES
//                                                    eigenvalue (CONST-G)
//     stokes_divfree_poincare_constant_h01_bound     H^1_0 VECTOR fields,
//                                                    DIV-FREE subspace: same
//                                                    (square root form;
//                                                    CONST-J)
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
//     stokes_ritz_projection_constants_h01_sq /      H^1_0 VECTOR fields:
//     stokes_ritz_projection_constant_set_h01_sq     the same frame as the
//                                                    composition above
//                                                    (CONST-J squared set)
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
// stokes_sobolev_embedding_constant_h01_measure<D, T>(pnum, pden, mea) and
// stokes_sobolev_embedding_constant_h01_spectrum<D, T>(pnum, pden, rho_lb,
// sigma) (CONST-J): the two routes of the entry above, each on its own --
// exactly sobolev_embedding_constant_h01_measure / _spectrum forwarded.
// Lemma 2 covers them for the same reason as the min selection: each scalar
// route's domain lies inside p >= 2, and out-of-domain input throws from the
// scalar route with the scalar message.  rho_lb of the spectrum route is a
// SCALAR Dirichlet lower bound, as on the scalar side.  D is given
// explicitly.
// ---------------------------------------------------------------------------
template <int D, typename T>
T stokes_sobolev_embedding_constant_h01_measure(long long pnum, long long pden,
                                                const T& mea) {
    return sobolev_embedding_constant_h01_measure<D, T>(pnum, pden, mea);
}

template <int D, typename T>
T stokes_sobolev_embedding_constant_h01_spectrum(long long pnum, long long pden,
                                                 const T& rho_lb,
                                                 const T& sigma) {
    return sobolev_embedding_constant_h01_spectrum<D, T>(pnum, pden, rho_lb,
                                                         sigma);
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
// Square-root free, so a rational T goes through exactly.  The square root
// form is stokes_poincare_constant_h01_bound right below (CONST-J), a thin
// wrapper of the scalar poincare_constant_h01_bound.
// ---------------------------------------------------------------------------
template <typename T>
T stokes_poincare_constant_h01_sq_bound(const T& lambda1_lower) {
    return poincare_constant_h01_sq_bound(lambda1_lower);
}

// ---------------------------------------------------------------------------
// stokes_poincare_constant_h01_bound(lambda1_lower) (CONST-J): the square
// root form C_P <= 1 / sqrt(lambda1_lower) of the entry above, exactly
// poincare_constant_h01_bound(lambda1_lower) forwarded (same argument, same
// Lemma 1 / Lemma 3 reading).  Same lazy policy as the scalar side: a
// rational T that has no square root still compiles as long as only the
// squared form is called.
// ---------------------------------------------------------------------------
template <typename T>
T stokes_poincare_constant_h01_bound(const T& lambda1_lower) {
    return poincare_constant_h01_bound(lambda1_lower);
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
// Three dimensions (CONST-K; this replaces the CONST-G deferral).  The three
// dimensional statement is in the literature: Liu, Nakao, Oishi, Commun.
// Nonlinear Sci. Numer. Simul. 108 (2022) 106223, section 3-2 (checked
// against the original by the CONST-K design side).  For problem (21),
// -Delta u + grad p = lambda u, div u = 0 in Omega in R^3, u = 0 on the
// boundary, with the discrete space "the CR space (vanishing boundary facet
// integrals) of three components whose restriction to every element K has
// div = 0", that section states
//     lambda_k >= lambda_{h,k} / (1 + (alpha_3 h)^2 lambda_{h,k}),
//     k = 1 .. Dim, h the largest tetrahedron diameter,
// resting on Theorem 2.1 of Liu 2015 (its Theorem 9).  alpha_3 is the
// EXISTING cr_interpolation_constant_3d_str, the diameter of a tetrahedron
// is its longest edge, and "div = 0 on K" is the constraint "int_K div v_h =
// 0" of vcr1_divfree.hpp because div of an elementwise linear field is
// elementwise constant [elementary].  So D == 3 uses the same pipeline with
// the same existing formula and C_h^2 = cr_projection_error_constant_sq of
// the scalar CR space (which switches to alpha_3 by itself); no new decimal
// string and no new formula.
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
    // CONST-K: the common fields of eigenvalue_report (same names, same
    // meaning).  DENSE PATH (to be replaced by the sparse path): index_lo /
    // index_hi are k..k for k = 1 .. nev (the k-th diagonal enclosure of
    // eigsymge, see detail::stokes_divfree_eigen_enclosures) and
    // nnz_stiffness / nnz_mass are 0 in two and three dimensions.
    std::vector<std::size_t> index_lo;
    std::vector<std::size_t> index_hi;
    long long nnz_stiffness;
    long long nnz_mass;
};

namespace detail {

// ---------------------------------------------------------------------------
// CONST-K (G-0 ruling M-1 and supplementary ruling 2): the discrete
// divergence-free space and its eigenvalue enclosures on the INTERVAL mesh.
//
// stokes_divfree_interval_basis<T, DP>: an interval enclosure of a basis
// matrix T* of the discrete divergence-free space (the space with Z), built
// by stokes_divfree_build_basis below:
//
//   1. C (n_elements x D n_r): row e = int_K div phi over the reduced vector
//      dofs, |K| (-D) d/dx_d lambda_i (cr1_element_op: phi_i = 1 - D
//      lambda_i), assembled in the interval scalar from the interval
//      geometry of the input mesh (no rationalization, no point-interval
//      requirement).
//   2. the LAST row is dropped.  The rows of the exact C sum to zero on the
//      reduced space: summing int_K div v over K gives the facet integrals,
//      the interior facet contributions cancel pairwise (same facet mean,
//      opposite normals) and the boundary facets are removed by the
//      Dirichlet reduction.  So rank C <= n_elements - 1, and the dropped
//      row is implied by the others.
//   3. rank = n_elements - 1 columns are chosen on the MIDPOINT of C by the
//      existing LU with partial pivoting of tlapack applied to mid(C)^T
//      (C = (C1, C2) after the split).  The choice only decides whether step 4 succeeds; it plays no
//      part in the correctness of what is returned.
//   4. X = an interval enclosure of C1^{-1} C2 by the verified dense linear
//      solver (lss: ||RA - I|| < 1, which also certifies that C1 is
//      nonsingular, hence rank C = n_elements - 1 exactly); a failure is
//      vcp::verification_error.
//   5. T = [-X; I] (pivot rows -X, free rows the identity) contains the exact
//      T* = [-C1^{-1} C2; I], whose columns are a basis of ker C (the
//      identity block is the independence witness, dim = D n_r - rank).  The
//      interval product C T over ALL rows (the dropped one included) must
//      contain zero entrywise; otherwise vcp::verification_error.
//
// Report meaning on this path: columns = D n_r, rank = the size of C1
// (n_elements - 1), rank_deficiency = 1, dim = columns - rank, pairs_checked
// = the number of entries of C T checked (n_elements * dim).
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct stokes_divfree_interval_basis {
    int n_elements;
    int n_reduced;                     // n_r
    int columns;                       // D n_r
    int rank;                          // size of C1 = n_elements - 1
    int dim;                           // columns - rank
    long long pairs_checked;           // entries of C T checked
    std::vector<int> pivot_cols;       // the columns of C1 (row p of X)
    std::vector<int> free_cols;        // the columns of C2 (column j of X)
    vcp::matrix<T, DP> X;              // enclosure of C1^{-1} C2, rank x dim
};

template <int D, typename T, class DP, class SP>
stokes_divfree_interval_basis<T, DP> stokes_divfree_build_basis(
        const vcp::bfem::vcr1_space<D, T, DP, SP>& V,
        const vcp::bfem::dirichlet_reduction<T, DP, SP>& red) {
    typedef typename T::base_type B;
    static const char who[] = "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01";
    stokes_divfree_interval_basis<T, DP> out;
    const int nt = V.num_elements();
    const int nr = red.reduced_size();
    const int nc = D * nr;
    out.n_elements = nt;
    out.n_reduced = nr;
    out.columns = nc;
    out.rank = nt - 1;
    out.dim = nc - out.rank;
    out.pairs_checked = 0;
    if (out.rank < 0 || out.rank > nc)
        vcp::throw_error<vcp::verification_error>(
            who, ": the divergence constraint has more independent rows than columns (rows - 1 = ",
            out.rank, ", columns = ", nc, ")");

    // 1. C, sparse rows in the interval scalar
    std::vector<std::vector<std::pair<int, T> > > crow(static_cast<std::size_t>(nt));
    for (int e = 0; e < nt; ++e) {
        const vcp::bfem::element_geometry<D, T>& g = V.scalar().geometry(e);
        const T w = g.measure() * T(-D);
        std::vector<std::pair<int, T> >& row = crow[static_cast<std::size_t>(e)];
        for (int i = 0; i <= D; ++i) {
            const int r = red.to_reduced(V.scalar().dofs().l2g(e, i));
            if (r < 0) continue;
            for (int d = 0; d < D; ++d)
                row.push_back(std::make_pair(d * nr + r, w * g.grad_lambda(i, d)));
        }
    }

    // 3. column choice on the midpoint of the first nt - 1 rows (approximation
    //    only): the existing LU with partial pivoting of tlapack (tgetrf, the
    //    LAPACK dgetrf semantics, ipiv 0-based) applied to mid(C)^T, nc x rank; its pivot
    //    rows are the chosen columns of C.  No new elimination code.
    std::vector<int> is_pivot(static_cast<std::size_t>(nc), 0);
    if (out.rank > 0) {
        const int m = out.rank;
        std::vector<B> Mt(static_cast<std::size_t>(nc) * static_cast<std::size_t>(m), B(0));
        for (int e = 0; e < m; ++e) {            // column e of mid(C)^T = row e of mid(C)
            const std::vector<std::pair<int, T> >& row = crow[static_cast<std::size_t>(e)];
            for (std::size_t q = 0; q < row.size(); ++q)
                Mt[static_cast<std::size_t>(e) * static_cast<std::size_t>(nc) +
                   static_cast<std::size_t>(row[q].first)] += mid(row[q].second);
        }
        std::vector<int> ipiv(static_cast<std::size_t>(m), 0);
        const int info = vcp::tgetrf(nc, m, Mt.data(), nc, ipiv.data());
        if (info != 0)
            vcp::throw_error<vcp::verification_error>(
                who, ": no independent columns found on the midpoint of the constraint "
                "(tgetrf info = ", info, "; the domain may be disconnected)");
        std::vector<int> perm(static_cast<std::size_t>(nc));
        for (int c = 0; c < nc; ++c) perm[static_cast<std::size_t>(c)] = c;
        for (int k = 0; k < m; ++k)
            std::swap(perm[static_cast<std::size_t>(k)],
                      perm[static_cast<std::size_t>(ipiv[static_cast<std::size_t>(k)])]);   // tlapack ipiv: 0-based
        for (int k = 0; k < m; ++k) {
            is_pivot[static_cast<std::size_t>(perm[static_cast<std::size_t>(k)])] = 1;
            out.pivot_cols.push_back(perm[static_cast<std::size_t>(k)]);
        }
    }
    std::vector<int> pos(static_cast<std::size_t>(nc), -1);   // column -> row of X / column of X
    for (int p = 0; p < out.rank; ++p) pos[static_cast<std::size_t>(out.pivot_cols[static_cast<std::size_t>(p)])] = p;
    for (int c = 0; c < nc; ++c)
        if (!is_pivot[static_cast<std::size_t>(c)]) {
            pos[static_cast<std::size_t>(c)] = static_cast<int>(out.free_cols.size());
            out.free_cols.push_back(c);
        }

    // 4. X = C1^{-1} C2 (verified)
    if (out.rank > 0 && out.dim > 0) {
        vcp::matrix<T, DP> C1, C2;
        C1.zeros(out.rank, out.rank);
        C2.zeros(out.rank, out.dim);
        for (int e = 0; e < out.rank; ++e) {
            const std::vector<std::pair<int, T> >& row = crow[static_cast<std::size_t>(e)];
            for (std::size_t q = 0; q < row.size(); ++q) {
                const int c = row[q].first;
                if (is_pivot[static_cast<std::size_t>(c)]) C1(e, pos[static_cast<std::size_t>(c)]) += row[q].second;
                else C2(e, pos[static_cast<std::size_t>(c)]) += row[q].second;
            }
        }
        out.X = lss(C1, C2);                     // throws vcp::verification_error
    } else {
        out.X.zeros(out.rank > 0 ? out.rank : 1, out.dim > 0 ? out.dim : 1);
    }

    // 5. C T contains zero (all rows, the dropped one included)
    for (int e = 0; e < nt; ++e) {
        const std::vector<std::pair<int, T> >& row = crow[static_cast<std::size_t>(e)];
        for (int j = 0; j < out.dim; ++j) {
            T sum(0);
            for (std::size_t q = 0; q < row.size(); ++q) {
                const int c = row[q].first;
                if (is_pivot[static_cast<std::size_t>(c)])
                    sum -= row[q].second * out.X(pos[static_cast<std::size_t>(c)], j);
                else if (pos[static_cast<std::size_t>(c)] == j)
                    sum += row[q].second;
            }
            ++out.pairs_checked;
            if (!(sum.lower() <= B(0) && B(0) <= sum.upper()))
                vcp::throw_error<vcp::verification_error>(
                    who, ": the interval product C T does not contain zero (row ", e,
                    ", column ", j, ")");
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// stokes_divfree_eigen_enclosures<D, T, DP, SP>(V, red, Zb, nev): the
// discrete eigenvalue enclosures of T^T blockdiag(A) T x = lambda T^T
// blockdiag(M) T x as an eig_result, k = 1 .. nev -- "(the space with Z,
// nev) -> eig_result".  The Stokes entry below calls only this function,
// stokes_divfree_build_basis (the space) and the common function
// indexed_eigenvalue_lower_bounds (poisson_constants.hpp); the report, the
// bound and the exception checks stay outside.
//
// What it does now (DENSE PATH): with T = [-X; I] of the basis above,
// W = blockdiag(K) T (through the sparsity of K), then (T^T W)(i, j) =
// W(free_i, j) - sum_p X(p, i) W(pivot_p, j) for j >= i, MIRRORED (exact
// symmetry for eigsymge), in the interval scalar for K = the scalar
// Dirichlet-reduced stiffness and mass; eigsymge (Miyajima, every
// eigenvalue enclosed for every matrix in the interval pencil); the k-th
// diagonal enclosure E(k, k) is returned with the index range k..k,
// converged = true and rank_certified = 1.  The reading "E(k, k) contains
// the k-th eigenvalue in ascending order" rests on the Miyajima enclosure
// named at imats::eigsymge (JCAM 236 (2012), its Theorem 3, enclosures in
// the order of the approximate eigenvalues, overlapping ones merged) and on
// the point side returning the approximate eigenvalues ascending (mats2:
// tsygv -> tsyev -> tsteqr, which sorts ascending).  The exact pencil of the
// true mesh (A*, M*, T*) lies inside the interval pencil, so its eigenvalues
// are enclosed.
//
// SPARSE-PATH REPLACEMENT POINT.  The ruling builds this function to be
// swapped for the sparse path (a) -- a sparse pencil and
// spmatrix::eigs_with_info -- with the API, the common function, the report
// fields and the tests unchanged.  The work that swap needs:
//   (1) a new entry that builds the reduced pencil as a spmatrix with a
//       pattern fixed symbolically (not from floating point cancellation);
//       the dense X above is not sparse, so a sparse basis is needed;
//   (2) spimats problem 1 resolved (below);
//   (3) the report filled from the sparse result: index_lo / index_hi from
//       the eig_result (no longer k..k) and nnz_stiffness / nnz_mass from
//       the finalized pencil (no longer 0);
//   (4) interval vertices accepted, as now (supplementary ruling 2).
// spimats problems recorded by CONST-K (sandbox/docs/reports/CONST-K_report.md,
// section "the current problems of spimats"):
//   1. lambda_min(B) is approximated without shift-invert, contrary to
//      veigs.m 'smallestabs'; fails for N > 300 on Z^T M Z.
//   2. the approximation stage needs the caller's options (target alone
//      does not converge).
//   3. mats<double>::eigsym inside the cluster refinement can hang,
//      contrary to the spimats note that it no longer reproduces.
//   4. P1 on the unit cube: a cluster cut by the approximation stage sends
//      the whole result to the Rayleigh fallback.
//   5. path (a) at n = 5: the approximation stage itself did not converge;
//      whether resolving 1 suffices is unconfirmed.
//
// Exceptions: vcp::verification_error from eigsymge (enclosures not
// certified).  nev <= dim is the caller's check.
// ---------------------------------------------------------------------------
template <typename T, class SP>
struct stokes_row_gather {
    std::vector<std::vector<std::pair<int, T> > >* rows;
    template <typename I, typename J>
    void operator()(I i, J j, const T& v) const {
        (*rows)[static_cast<std::size_t>(i)].push_back(std::make_pair(static_cast<int>(j), v));
    }
};

template <int D, typename T, class DP, class SP>
vcp::matrix<T, DP> stokes_divfree_reduce_dense(
        const stokes_divfree_interval_basis<T, DP>& Zb,
        const vcp::spmatrix<T, SP>& K) {
    const int nr = Zb.n_reduced, nc = Zb.columns, N = Zb.dim, m = Zb.rank;
    std::vector<std::vector<std::pair<int, T> > > krow(static_cast<std::size_t>(nr));
    stokes_row_gather<T, SP> g = { &krow };
    ::vcp::bfem::detail::spm_adapter<T, SP>::for_each_entry(K, g);
    // the rows of T: pivot column -> (-X row p), free column -> unit row
    std::vector<int> prow(static_cast<std::size_t>(nc), -1), frow(static_cast<std::size_t>(nc), -1);
    for (int p = 0; p < m; ++p) prow[static_cast<std::size_t>(Zb.pivot_cols[static_cast<std::size_t>(p)])] = p;
    for (int j = 0; j < N; ++j) frow[static_cast<std::size_t>(Zb.free_cols[static_cast<std::size_t>(j)])] = j;
    // W = blockdiag(K) T, nc x N
    vcp::matrix<T, DP> W;
    W.zeros(nc, N);
    for (int d = 0; d < D; ++d)
        for (int r = 0; r < nr; ++r) {
            const int gi = d * nr + r;
            const std::vector<std::pair<int, T> >& kr = krow[static_cast<std::size_t>(r)];
            for (std::size_t q = 0; q < kr.size(); ++q) {
                const int gj = d * nr + kr[q].first;
                const T& kv_ = kr[q].second;
                if (frow[static_cast<std::size_t>(gj)] >= 0) {
                    W(gi, frow[static_cast<std::size_t>(gj)]) += kv_;
                } else {
                    const int p = prow[static_cast<std::size_t>(gj)];
                    for (int j = 0; j < N; ++j) W(gi, j) -= kv_ * Zb.X(p, j);
                }
            }
        }
    // T^T W, upper triangle, mirrored
    vcp::matrix<T, DP> out;
    out.zeros(N, N);
    for (int i = 0; i < N; ++i) {
        const int fi = Zb.free_cols[static_cast<std::size_t>(i)];
        for (int j = i; j < N; ++j) {
            T s = W(fi, j);
            for (int p = 0; p < m; ++p)
                s -= Zb.X(p, i) * W(Zb.pivot_cols[static_cast<std::size_t>(p)], j);
            out(i, j) = s;
            if (j != i) out(j, i) = s;
        }
    }
    return out;
}

template <int D, typename T, class DP, class SP>
vcp::eig_result<T> stokes_divfree_eigen_enclosures(
        const vcp::bfem::vcr1_space<D, T, DP, SP>& V,
        const vcp::bfem::dirichlet_reduction<T, DP, SP>& red,
        const stokes_divfree_interval_basis<T, DP>& Zb,
        int nev) {
    vcp::matrix<T, DP> A = stokes_divfree_reduce_dense<D, T, DP, SP>(
        Zb, red.reduce(V.scalar().stiffness()));
    vcp::matrix<T, DP> M = stokes_divfree_reduce_dense<D, T, DP, SP>(
        Zb, red.reduce(V.scalar().mass()));
    vcp::matrix<T, DP> E;
    eigsymge(A, M, E);                   // verified enclosures, ascending
    vcp::eig_result<T> r;
    r.eigenvalues.reserve(static_cast<std::size_t>(nev));
    for (int k = 0; k < nev; ++k) {
        r.eigenvalues.push_back((E.columnsize() == 1) ? E(k, 0) : E(k, k));
        r.index_lo.push_back(static_cast<std::size_t>(k) + 1);
        r.index_hi.push_back(static_cast<std::size_t>(k) + 1);
    }
    r.converged = true;
    r.rank_certified = 1;
    r.requested_count = static_cast<std::size_t>(nev);
    r.returned_count = static_cast<std::size_t>(nev);
    r.returned_real_count = static_cast<std::size_t>(nev);
    r.converged_count = static_cast<std::size_t>(nev);
    r.used_method = "dense_eigsymge";
    return r;
}

} // namespace detail

// ---------------------------------------------------------------------------
// stokes_eigenvalue_lower_bounds_h01<D, T, DP, SP>(Th, nev, rep): verified
// lower bounds of the first nev Stokes eigenvalues on the mesh Th, i.e. the
// vector ( lambda_{k,h} / (1 + C_h^2 lambda_{k,h}) )_{k = 1..nev} with every
// lambda_{k,h} a verified enclosure of the k-th eigenvalue of the reduced
// pencil T^T A T x = lambda T^T M T x (A, M the block-diagonal broken
// stiffness / mass of the vector CR space after the homogeneous Dirichlet
// reduction, T an interval basis of the kernel of the element-mean
// divergence constraint) and C_h^2 = cr_projection_error_constant_sq of the
// scalar CR space.  Entry k of the result is a T whose LOWER end is the
// guaranteed lower bound of lambda_{k+1}; the upper end carries no claim
// beyond enclosing the formula's value.  D = 2 or 3 (the sources of both
// dimensions are cited at the head of the CONST-G section above).
//
// Pipeline (CONST-K, G-0 ruling M-1 and supplementary ruling 2):
// vcr1_space (interval vertices accepted) -> dirichlet_reduction of the
// scalar space (reduced vector numbering d * n_r + r) ->
// detail::stokes_divfree_build_basis (interval constraint C, last row
// dropped, midpoint column choice, verified C1^{-1} C2, C T contains zero)
// -> detail::stokes_divfree_eigen_enclosures (T^T A T and T^T M T in the
// interval scalar, eigsymge, E(k, k) numbered k..k: the dense path and the
// replacement point of the planned sparse path) -> the common function
// indexed_eigenvalue_lower_bounds of poisson_constants.hpp (rules 1..3; the
// bound at the LOWER END of each enclosure as a point, G-0 supplementary
// ruling -- the two dimensional values may move in the last digits).
//
// The domain the bound speaks about.  The bound is a bound for the
// polyhedron that the given mesh fills.  If boundary vertices are given as
// points that do not lie exactly on the boundary of the intended domain,
// the bound is NOT a bound for that domain.  Vertices given as intervals
// that contain the true positions give a correct bound for the intended
// domain (every enclosure of this pipeline holds for every mesh inside the
// vertex boxes).
//
// Report: as stokes_eigenvalue_report, with columns = D n_r, rank = the size
// of C1 (n_elements - 1), rank_deficiency = 1, dim = columns - rank,
// pairs_checked = the entries of C T checked (n_elements * dim); index_lo /
// index_hi = k..k and nnz_stiffness / nnz_mass = 0 on the dense path.
//
// Policy defaults (ruling M-8): DP = imats<B, mats2<B> >, SP = spimats<B,
// spmats<B>, imats<B, mats2<B> > >, B = T::base_type, the eigenvalue
// backend guideline of spimats.hpp, the same as eigenvalue_lower_bounds_h01
// (DP = vcp::pidblas is the other choice of that guideline).
//
// T: a kv interval scalar (T::base_type required by the policy defaults).
// Exceptions: std::invalid_argument (empty mesh, nev < 1, nev > dim -- the
// lower bound theorem covers 1 <= k <= dim only), vcp::verification_error
// (no independent columns, lss on C1 not certified, C T not containing
// zero, eigsymge not certified, rules 1 and 3 of the common function),
// vcp::bfem::degenerate_element (flat element).
// ---------------------------------------------------------------------------
template <int D, typename T,
          class DP = vcp::imats<typename T::base_type,
                                vcp::mats2<typename T::base_type> >,
          class SP = vcp::spimats<typename T::base_type,
                                  vcp::spmats<typename T::base_type>,
                                  vcp::imats<typename T::base_type,
                                             vcp::mats2<typename T::base_type> > > >
std::vector<T> stokes_eigenvalue_lower_bounds_h01(
        const vcp::bfem::mesh<D, T>& Th, int nev,
        stokes_eigenvalue_report<T>* rep = nullptr) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01: "
                  "two or three dimensions only");
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

    // CONST-K (supplementary ruling 2): the interval basis of the discrete
    // divergence-free space on the interval mesh as given (no
    // rationalization; the exact null space path of vcr1_divfree.hpp stays in
    // the library but is not used here)
    const detail::stokes_divfree_interval_basis<T, DP> Zb =
        detail::stokes_divfree_build_basis<D, T, DP, SP>(V, red);
    const int dim = Zb.dim;
    if (nev > dim)
        throw std::invalid_argument(
            "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01: nev exceeds "
            "dim V_{0,h} (the lower bound theorem covers 1 <= k <= dim only)");
    // the enclosures (replacement point) and the common function
    const vcp::eig_result<T> er =
        detail::stokes_divfree_eigen_enclosures<D, T, DP, SP>(V, red, Zb, nev);

    // the existing lower bound formula with the existing CR constant
    const T ch_sq = cr_projection_error_constant_sq(V.scalar());
    const detail::indexed_eigenvalue_bounds<T> b = detail::indexed_eigenvalue_lower_bounds(
        er, nev, ch_sq, "vcp::bfem::constants::stokes_eigenvalue_lower_bounds_h01");
    if (rep != nullptr) {
        rep->n_elements = Th.num_elements();
        rep->n_scalar_dofs = V.scalar_ndof();
        rep->n_reduced_scalar = red.reduced_size();
        rep->columns = Zb.columns;            // D n_r
        rep->rank = Zb.rank;                  // size of C1 = n_elements - 1
        rep->rank_deficiency = Zb.n_elements - Zb.rank;   // 1
        rep->dim = dim;                       // columns - rank
        rep->pairs_checked = Zb.pairs_checked; // entries of C T checked
        rep->ch_sq = ch_sq;
        rep->lambda_h = b.lambda_h;
        rep->index_lo = b.index_lo;           // dense path: k..k
        rep->index_hi = b.index_hi;
        rep->nnz_stiffness = 0;               // dense path (sparse path planned)
        rep->nnz_mass = 0;
    }
    return b.bounds;
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
// hence exact on a rational T (R17); the square root form is
// stokes_divfree_poincare_constant_h01_bound right below (CONST-J).
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

// ---------------------------------------------------------------------------
// stokes_divfree_poincare_constant_h01_bound(stokes_lambda1_lower)
// (CONST-J): the square root form C_P^sigma <= 1 / sqrt(lambda_1^S,lower),
// defined as the square root of the squared entry above (same argument: a
// verified lower bound of the first STOKES eigenvalue).  Instantiated only
// when called, so a rational T keeps using the squared form.
// ---------------------------------------------------------------------------
template <typename T>
T stokes_divfree_poincare_constant_h01_bound(const T& stokes_lambda1_lower) {
    using std::sqrt;
    return sqrt(stokes_divfree_poincare_constant_h01_sq_bound(stokes_lambda1_lower));
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

// ---------------------------------------------------------------------------
// stokes_ritz_projection_constant_set_h01_sq<T> (CONST-J): the squared
// constant set of one Stokes hypercircle pass, the Stokes twin of
// ritz_projection_constant_set_h01_sq of the Poisson side.
//
//   kappa2       kappa_h^2 (as stokes_hypercircle_kappa_h01_squared)
//   c0h_sq       C_{0,h}(k - 1)^2 of the coarse mesh: dictionary served, or
//                the value the caller supplied (explicit form)
//   ch_sq        C_h(k)^2 = kappa2 + c0h_sq
//   hypercircle  the kappa report (sizes, certificate counts, spectrum)
//   resolution   the dictionary census of the c0h_sq pass; left default
//                constructed by the explicit form (no dictionary call)
// ---------------------------------------------------------------------------
template <typename T>
struct stokes_ritz_projection_constant_set_h01_sq {
    T kappa2;
    T c0h_sq;
    T ch_sq;
    stokes_hypercircle_report<T> hypercircle;
    mesh_resolution_report resolution;

    stokes_ritz_projection_constant_set_h01_sq()
        : kappa2(), c0h_sq(), ch_sq(), hypercircle(), resolution() {}
};

// ---------------------------------------------------------------------------
// stokes_ritz_projection_constants_h01_sq<T, DP, SP>(Th, k [, c0h_sq])
// (CONST-J): the set above, from the same calls in the same order as the
// two compositions above -- stokes_hypercircle_kappa_h01_squared with its
// report, then stokes_l2_projection_error_constant_sq(Th, k - 1) with its
// census (or the supplied value, which is the caller's responsibility),
// then ch_sq = kappa2 + c0h_sq.  Exceptions as the matching composition.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
stokes_ritz_projection_constant_set_h01_sq<T>
stokes_ritz_projection_constants_h01_sq(const vcp::bfem::mesh<3, T>& Th, int k) {
    stokes_ritz_projection_constant_set_h01_sq<T> s;
    s.kappa2 = stokes_hypercircle_kappa_h01_squared<T, DP, SP>(Th, k, &s.hypercircle);
    s.c0h_sq = stokes_l2_projection_error_constant_sq(Th, k - 1, &s.resolution);
    s.ch_sq = s.kappa2 + s.c0h_sq;
    return s;
}

template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
stokes_ritz_projection_constant_set_h01_sq<T>
stokes_ritz_projection_constants_h01_sq(const vcp::bfem::mesh<3, T>& Th, int k,
                                        const T& c0h_sq) {
    stokes_ritz_projection_constant_set_h01_sq<T> s;
    s.kappa2 = stokes_hypercircle_kappa_h01_squared<T, DP, SP>(Th, k, &s.hypercircle);
    s.c0h_sq = c0h_sq;
    s.ch_sq = s.kappa2 + s.c0h_sq;
    return s;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP
