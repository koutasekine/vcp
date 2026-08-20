// vcp/bfem/constants/stokes_constants.hpp
//
// Verified constants of the Stokes problem -- SKELETON ONLY (VER-0, R10).
// This header deliberately contains no function and no type: it fixes the
// file name and the namespace so that the Stokes side of the constants layer
// has its place next to poisson_constants.hpp.
//
// Planned contents (future tracks; see the VER-0 design, section 2):
//
//     stokes_ritz_projection_error_constant   the Stokes counterpart of the
//                                             Ritz projection error constant
//     a div-free Poincare constant            (divergence-free test space)
//     the hypercircle kappa_h of the mixed    Stokes formulation
//     C_{4,P}                                 (appears on the Stokes side
//                                             when it is needed; R10)
//
// These are to be given bodies only after the vcr1 and the local div-free
// constraint tracks are complete; nothing here is usable until then.
//
// Boundary-condition note (VER-1, same ledger format as the one at the head
// of poisson_constants.hpp): the planned entries above are all framed on
// H^1_0 VECTOR fields -- the div-free test space and the mixed hypercircle
// assume homogeneous Dirichlet velocity -- so, when they are given bodies,
// their names carry the BC tag of ruling R23 (concept name, then BC tag,
// then form suffix).
//
// Lexical policy (inherited from poisson_constants.hpp, design 6.1): no
// decimal literals, no `double` / `float` tokens in code; sqrt and
// kv::constants<T>::pi() are allowed.
//
// static_assert convention: as in poisson_constants.hpp -- whoever adds to
// this layer a function that takes a SPACE as an argument and assumes a
// conforming space must place the corresponding
// static_assert(vcp::bfem::space_traits<S>::..., ...) in it.
//
// Authority: sandbox/docs/design/VER-0_design_v1.0.md.

#ifndef VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP
#define VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP

namespace vcp {
namespace bfem {
namespace constants {

// intentionally empty (VER-0: skeleton only)

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_STOKES_CONSTANTS_HPP
