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
//     a div-free Poincare constant            (the SHARP one of Lemma 3 --
//                                             the wrapper below is only the
//                                             full-space substitute)
//     the hypercircle kappa_h of the mixed    Stokes formulation
//     C_{4,P}                                 (appears on the Stokes side
//                                             when it is needed; R10)
//
// These are still to be given bodies only after the vcr1 and the local
// div-free constraint tracks are complete; the C_h of the Stokes projection
// itself remains out of scope.
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

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/constants/poisson_constants.hpp>

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

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP
