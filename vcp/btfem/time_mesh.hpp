// vcp/btfem/time_mesh.hpp
// btfem T1: uniform time mesh time_mesh<T> and degenerate_interval.
//
// Conforms to: BT-1 external design v1.1 (section 3) and
//              BT-1 internal design v1.0 (section 2).
//
// Division-once contract (external 8, internal 2): the ONLY T division of
// this header is inv_k_ = T(m) / (t1 - t0) in uniform(). The width k itself
// is formed as (t1 - t0) * (1/m) where 1/m is an enclose-once rational
// constant (convert_traits<T>, shared per (T, m) like bfem's 1/D! constant);
// its one-time conversion division belongs to the L0 conversion budget, not
// to the per-mesh division count (same convention as bfem geometry.hpp).

#ifndef VCP_BTFEM_TIME_MESH_HPP
#define VCP_BTFEM_TIME_MESH_HPP

#include <string>
#include <stdexcept>
#include <mutex>

#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/detail/table_cache.hpp>

namespace vcp {
namespace btfem {

class degenerate_interval : public std::runtime_error {
public:
    explicit degenerate_interval(const std::string& msg) : std::runtime_error(msg) {}
};

namespace detail {

// enclose-once 1/m per (T, m): process shared, generated lazily
template <typename T>
class inverse_int_cache {
public:
    static const T& get(int m) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.inv.get_or_build(m, [&]() -> T {
            return bfem::rational_to<T>(1, m);
        });
    }
private:
    struct maps { bfem::detail::table_cache<int, T> inv; };
    typedef bfem::detail::table_cache_state<maps> state;
    static state& st() { return bfem::detail::table_cache_instance<state>(); }
};

} // namespace detail

// ---------------------------------------------------------------------------
// time_mesh<T>: [t0, t1] divided into m equal sub-intervals J_l = [tau_{l-1},
// tau_l], l = 1..m, width k = (t1 - t0)/m (R14: uniform only).
// ---------------------------------------------------------------------------
template <typename T>
class time_mesh {
public:
    static time_mesh uniform(const T& t0, const T& t1, int m) {
        if (m < 1)
            throw std::invalid_argument("btfem::time_mesh::uniform: m < 1");
        T len = t1 - t0;
        int sg;
        try {
            sg = bfem::geometry_traits<T>::sign(len);
        } catch (const bfem::degenerate_element&) {
            throw degenerate_interval(
                "btfem::time_mesh::uniform: interval length is zero or sign-indefinite");
        }
        if (sg != +1)
            throw degenerate_interval("btfem::time_mesh::uniform: t1 - t0 is not positive");
        time_mesh tm;
        tm.t0_ = t0;
        tm.t1_ = t1;
        tm.m_ = m;
        tm.k_ = len * detail::inverse_int_cache<T>::get(m);   // multiplication only
        tm.inv_k_ = T(m) / len;                                // the single division
        return tm;
    }

    const T& t_begin() const { return t0_; }
    const T& t_end() const { return t1_; }
    int num_intervals() const { return m_; }
    const T& width() const { return k_; }          // k (common to all sub-intervals)
    const T& inv_width() const { return inv_k_; }  // 1/k (independent value; k*inv_k is
                                                   // only an enclosure of 1 for interval T)
    // tau_l = t0 + l k, l = 0..m (one multiplication, one addition)
    T node(int l) const {
        if (l < 0 || l > m_)
            throw std::invalid_argument("btfem::time_mesh::node: l out of range");
        return t0_ + T(l) * k_;
    }

    // value identity of two meshes (internal 3.2): equal m and non-separable
    // end points. For point types this is exact equality; for interval types
    // (whose == requires degenerate intervals) it is the overlap test
    // !(a < b) && !(b < a), which identical intervals always pass.
    bool same_as(const time_mesh& o) const {
        return m_ == o.m_ && same_value(t0_, o.t0_) && same_value(t1_, o.t1_);
    }

private:
    time_mesh() : t0_(), t1_(), k_(), inv_k_(), m_(0) {}
    T t0_, t1_, k_, inv_k_;
    int m_;

    static bool same_value(const T& a, const T& b) {
        return !(a < b) && !(b < a);
    }
};

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_TIME_MESH_HPP
