// vcp/bfem/geometry_traits_kv.hpp
// Layer 2: opt-in geometry_traits specialization for kv::interval<F>
// (L2 external design v0.3, section 3.2 / internal design 2.3).
//
// sup(det) < 0 -> -1, inf(det) > 0 -> +1, otherwise (0 enclosed) the element
// is not verifiable and degenerate_element is thrown -- an unverifiable
// element never passes silently.

#ifndef VCP_BFEM_GEOMETRY_TRAITS_KV_HPP
#define VCP_BFEM_GEOMETRY_TRAITS_KV_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/bfem/geometry.hpp>

namespace vcp {
namespace bfem {

template <typename F>
struct geometry_traits<kv::interval<F> > {
    static int sign(const kv::interval<F>& det) {
        if (det.upper() < F(0)) return -1;
        if (det.lower() > F(0)) return +1;
        throw degenerate_element(
            "bfem::element_geometry: interval det encloses zero (unverifiable element)");
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_GEOMETRY_TRAITS_KV_HPP
