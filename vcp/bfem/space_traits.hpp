// vcp/bfem/space_traits.hpp
// CR1-1 (ruling R2): the centrally placed trait that records, per finite
// element space, whether the space is CONFORMING in H^1.
//
// Conforms to: CR1-1 design v1.0 (sections 0 R2, 3.4) and
//              CR1-1 implementation directive v1.0 (phase P4).
//
// Placement rationale (R2).  The property belongs to the space, but writing it
// into each space header would mean editing five frozen headers for a purely
// declarative fact.  This header therefore FORWARD DECLARES the five class
// templates and specializes space_traits for them.  It includes none of them
// and modifies none of them; bfem stays a pure addition.  Default template
// arguments are deliberately absent from the forward declarations -- they live
// on the definitions, and repeating them here would be an error.
//
// The primary template is left UNDEFINED on purpose: asking for the trait of a
// space that has not been classified is a compile error rather than a silent
// default.
//
// Convention for later tracks (the assert is NOT planted here -- R2 defers it
// to the constants track): a function whose correctness rests on
// V_h subset of H^1_0 -- a projection error constant, a conforming-space
// eigenvalue upper bound, anything that integrates a global gradient rather
// than a broken one -- should open with
//
//     static_assert(space_traits<S>::h1_conforming,
//                   "<function>: requires an H^1 conforming space");
//
// so that handing it a nonconforming space fails at compile time instead of
// returning a number that means nothing.
//
// Classification and why:
//  - fe_space   : the Lagrange P^n family, continuous across facets      -> true
//  - vfe_space  : a componentwise composition of conforming scalars      -> true
//  - c1_space   : the Argyris C^1 family, and C^1 is contained in H^1    -> true
//  - rt_space   : Raviart-Thomas, conforming in H(div) but NOT in H^1
//                 (only the normal trace is continuous)                  -> false
//  - cr1_space  : Crouzeix-Raviart, only the facet MEAN is continuous    -> false

#ifndef VCP_BFEM_SPACE_TRAITS_HPP
#define VCP_BFEM_SPACE_TRAITS_HPP

namespace vcp {
namespace bfem {

// ---- forward declarations (argument lists copied from the definitions;
//      no default arguments here) ----
template <int D, typename T, typename P, class SP> class fe_space;    // vcp/bfem/fe_space.hpp
template <int D, typename T, typename P, class SP> class vfe_space;   // vcp/bfem/sv/vfe_space.hpp
template <int D, typename T, typename P, class SP> class c1_space;    // vcp/bfem/c1/c1_space.hpp
template <int D, typename T, typename P, class SP> class rt_space;    // vcp/bfem/rt/rt_space.hpp
template <int D, typename T, typename P, class SP> class cr1_space;   // vcp/bfem/cr1/cr1_space.hpp

// ---- the trait ----
template <class S>
struct space_traits;                 // intentionally undefined for unclassified S

template <int D, typename T, typename P, class SP>
struct space_traits<fe_space<D, T, P, SP> > {
    static const bool h1_conforming = true;
};

template <int D, typename T, typename P, class SP>
struct space_traits<vfe_space<D, T, P, SP> > {
    static const bool h1_conforming = true;
};

template <int D, typename T, typename P, class SP>
struct space_traits<c1_space<D, T, P, SP> > {
    static const bool h1_conforming = true;
};

template <int D, typename T, typename P, class SP>
struct space_traits<rt_space<D, T, P, SP> > {
    static const bool h1_conforming = false;
};

template <int D, typename T, typename P, class SP>
struct space_traits<cr1_space<D, T, P, SP> > {
    static const bool h1_conforming = false;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SPACE_TRAITS_HPP
