// vcp/bfem/detail/table_cache.hpp
// L4 seam alignment (design L4_concept_design_v0.1 section 10, item 1): the
// single shared cache vessel replacing the four hand-rolled copies of the
// "magic static + mutex + std::map" pattern (L0 coeff/typed registries,
// L2 ref_stiffness, rt rational/typed registries, c1 rational/typed
// registries).
//
// Contract (identical to every replaced implementation):
//  - lazy generation: the factory runs only on a lookup miss;
//  - returned references stay valid until program termination;
//  - concurrency safe: one mutex per registry x stage. The LOCK GRANULARITY
//    IS PRESERVED -- table_cache itself does NOT lock. Each registry-stage
//    state owns exactly one mutex (table_cache_state) covering all of its
//    keyed maps, exactly as before; the rational-stage builders may therefore
//    keep filling sibling maps of the same stage while holding that lock.
//  - the factory is injected per call; generated values, generation order
//    and exception behavior (a throwing factory inserts nothing) are those
//    of the caller's unchanged generation code.
//
// This header adds NO new capability -- it is the pure-refactor vessel only.

#ifndef VCP_BFEM_DETAIL_TABLE_CACHE_HPP
#define VCP_BFEM_DETAIL_TABLE_CACHE_HPP

#include <map>
#include <mutex>
#include <utility>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// table_cache<Key, Table>: one keyed table map of a registry-stage cache.
// get_or_build must be called with the owning stage's mutex held for the
// whole call (the replaced find/insert idiom, verbatim).
// ---------------------------------------------------------------------------
template <typename Key, typename Table>
class table_cache {
public:
    template <typename Factory>
    const Table& get_or_build(const Key& key, const Factory& make) {
        typename std::map<Key, Table>::iterator it = map_.find(key);
        if (it != map_.end()) return it->second;
        Table t = make();
        return map_.insert(std::make_pair(key, std::move(t))).first->second;
    }

private:
    std::map<Key, Table> map_;
};

// ---------------------------------------------------------------------------
// table_cache_state<Maps>: a registry-stage state = the stage's table_cache
// members (Maps) + the one stage mutex. Inheriting Maps keeps the member
// access syntax (s.mass, s.imaps, ...) of the replaced hand-rolled states.
// ---------------------------------------------------------------------------
template <typename Maps>
struct table_cache_state : Maps {
    std::mutex mtx;
};

// ---------------------------------------------------------------------------
// table_cache_instance<State>(): the magic-static holder (C++11 thread-safe
// initialization). State is a registry-private type, so each registry x
// stage keeps its own independent instance -- cache instance granularity is
// unchanged.
// ---------------------------------------------------------------------------
template <typename State>
State& table_cache_instance() {
    static State s;
    return s;
}

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_DETAIL_TABLE_CACHE_HPP
