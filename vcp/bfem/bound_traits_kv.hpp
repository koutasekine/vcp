// vcp/bfem/bound_traits_kv.hpp
// Layer 1: opt-in bound_traits specialization for kv::interval<F>
// (L1 external design v0.3, section 8.3).
//
// Specialization contract: after the scan, acc encloses the true minimum
// (resp. maximum) coefficient value. With intervals [a,b], [c,d] the running
// enclosure of the min is [min(a,c), min(b,d)] and of the max is
// [max(a,c), max(b,d)]; both preserve the enclosure inductively, which is the
// premise of the range/range_refined inclusion property (external 8.1, 8.2).

#ifndef VCP_BFEM_BOUND_TRAITS_KV_HPP
#define VCP_BFEM_BOUND_TRAITS_KV_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/bfem/refine.hpp>

namespace vcp {
namespace bfem {

template <typename F>
struct bound_traits<kv::interval<F> > {
    static void min_update(kv::interval<F>& acc, const kv::interval<F>& x) {
        if (x.lower() < acc.lower()) acc.lower() = x.lower();
        if (x.upper() < acc.upper()) acc.upper() = x.upper();
    }
    static void max_update(kv::interval<F>& acc, const kv::interval<F>& x) {
        if (x.lower() > acc.lower()) acc.lower() = x.lower();
        if (x.upper() > acc.upper()) acc.upper() = x.upper();
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_BOUND_TRAITS_KV_HPP
