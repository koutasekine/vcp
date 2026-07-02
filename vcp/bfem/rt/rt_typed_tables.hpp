// vcp/bfem/rt/rt_typed_tables.hpp
// RT Layer 0: T stage of the two-stage cache -- typed_rt_registry<D, T>.
//
// Conforms to: RT-L0 external design v0.3 (section 3) and
//              RT-L0 internal design v0.3 (section 7).
//
// Every entry of every table is the image of the RATIONAL-STAGE value under
// convert_traits<T> exactly once (enclose-once). All geometry-free
// contractions already happened at the rational stage (external design
// 1.3-3); no T arithmetic beyond the conversion occurs here, so for interval
// T every stored entry is an enclosure of the exact rational value and the
// endpoints are bit-reproducible (structural check S-RT0-1).
//
// Accessor surface is identical to the rational stage by construction (the
// table shells are scalar-generic, rt_tables.hpp).

#ifndef VCP_BFEM_RT_RT_TYPED_TABLES_HPP
#define VCP_BFEM_RT_RT_TYPED_TABLES_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/rt/rt_tables.hpp>

namespace vcp {
namespace bfem {

namespace detail {

template <typename T>
inline T rt_conv(const rational& q) {
    return convert_traits<T>::from_rational(q.num(), q.den());
}

template <typename T>
inline std::vector<T> rt_conv_vec_mat(const rt_mat_tbl<rational>& s) {
    std::vector<T> v;
    v.reserve(static_cast<std::size_t>(s.rows()) * static_cast<std::size_t>(s.cols()));
    for (int i = 0; i < s.rows(); ++i)
        for (int j = 0; j < s.cols(); ++j)
            v.push_back(rt_conv<T>(s.at(i, j)));
    return v;
}

} // namespace detail

// ---------------------------------------------------------------------------
// typed_rt_registry<D, T> (external design section 3: same shape as
// rt_registry, values of type T). One mutex per (D, T); different T
// initialize concurrently. References stay valid until program termination.
// ---------------------------------------------------------------------------
template <int D, typename T>
class typed_rt_registry {
    static_assert(D == 2,
                  "bfem::typed_rt_registry: initial version supports D == 2 only");
public:
    typedef rt_basis_tbl<T> basis_table;
    typedef rt_div_tbl<T>   div_table;
    typedef rt_flux_tbl<T>  flux_table;
    typedef rt_block_tbl<T> block_table;
    typedef rt_cross_tbl<T> cross_table;
    typedef rt_mat_tbl<T>   mat_table;

    static int dim(int k) { return rt_registry<D>::dim(k); }

    static const basis_table& basis(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, basis_table>::iterator it = s.basis.find(k);
        if (it != s.basis.end()) return it->second;
        const rt_basis_table& src = rt_registry<D>::basis(k);
        std::vector<T> v;
        v.reserve(static_cast<std::size_t>(src.rows())
                  * static_cast<std::size_t>(src.cols()));
        for (int i = 0; i < src.rows(); ++i)
            for (int j = 0; j < src.cols(); ++j)
                v.push_back(detail::rt_conv<T>(src.at(i, j)));
        basis_table t = detail::rt_table_access::make_basis(
            src.order(), src.comp_size(), src.dim(), std::move(v));
        return s.basis.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static const div_table& divergence(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, div_table>::iterator it = s.div.find(k);
        if (it != s.div.end()) return it->second;
        const rt_div_table& src = rt_registry<D>::divergence(k);
        std::vector<T> v;
        v.reserve(static_cast<std::size_t>(src.rows())
                  * static_cast<std::size_t>(src.cols()));
        for (int i = 0; i < src.rows(); ++i)
            for (int j = 0; j < src.cols(); ++j)
                v.push_back(detail::rt_conv<T>(src.at(i, j)));
        div_table t = detail::rt_table_access::make_div(
            src.order(), src.rows(), src.cols(), std::move(v));
        return s.div.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static const flux_table& edge_flux(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, flux_table>::iterator it = s.flux.find(k);
        if (it != s.flux.end()) return it->second;
        const rt_flux_table& src = rt_registry<D>::edge_flux(k);
        const int kk = src.order();
        const int dm = src.dim();
        std::vector<T> v;
        v.reserve(static_cast<std::size_t>(3) * static_cast<std::size_t>(kk + 1)
                  * static_cast<std::size_t>(dm));
        for (int e = 0; e < 3; ++e)
            for (int j = 0; j <= kk; ++j)
                for (int c = 0; c < dm; ++c)
                    v.push_back(detail::rt_conv<T>(src.at(e, j, c)));
        flux_table t = detail::rt_table_access::make_flux(kk, dm, std::move(v));
        return s.flux.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static const block_table& comp_mass(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, block_table>::iterator it = s.cmass.find(k);
        if (it != s.cmass.end()) return it->second;
        block_table t = conv_block(rt_registry<D>::comp_mass(k));
        return s.cmass.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static const block_table& div_mass(int k, int l) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, l);
        typename std::map<std::pair<int, int>, block_table>::iterator it =
            s.dmass.find(key);
        if (it != s.dmass.end()) return it->second;
        block_table t = conv_block(rt_registry<D>::div_mass(k, l));
        return s.dmass.insert(std::make_pair(key, std::move(t))).first->second;
    }

    static const cross_table& cross_grad(int k, int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, n);
        typename std::map<std::pair<int, int>, cross_table>::iterator it =
            s.cross.find(key);
        if (it != s.cross.end()) return it->second;
        const rt_cross_table& src = rt_registry<D>::cross_grad(k, n);
        std::vector<std::vector<T> > blk;
        blk.reserve(6);
        for (int d = 0; d < 2; ++d) {
            for (int i = 0; i <= 2; ++i) {
                rt_block_view<detail::rational> B = src.block(d, i);
                std::vector<T> b;
                b.reserve(static_cast<std::size_t>(B.rows())
                          * static_cast<std::size_t>(B.cols()));
                for (int r = 0; r < B.rows(); ++r)
                    for (int c = 0; c < B.cols(); ++c)
                        b.push_back(detail::rt_conv<T>(B.at(r, c)));
                blk.push_back(std::move(b));
            }
        }
        cross_table t = detail::rt_table_access::make_cross(
            src.dim(), src.nn(), std::move(blk));
        return s.cross.insert(std::make_pair(key, std::move(t))).first->second;
    }

    static const mat_table& cross_grad_contracted(int k, int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, n);
        typename std::map<std::pair<int, int>, mat_table>::iterator it =
            s.crossc.find(key);
        if (it != s.crossc.end()) return it->second;
        const rt_mat_table& src = rt_registry<D>::cross_grad_contracted(k, n);
        mat_table t = detail::rt_table_access::make_mat(
            src.rows(), src.cols(), detail::rt_conv_vec_mat<T>(src));
        return s.crossc.insert(std::make_pair(key, std::move(t))).first->second;
    }

    static const mat_table& inv_mass(int l) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, mat_table>::iterator it = s.invm.find(l);
        if (it != s.invm.end()) return it->second;
        const rt_mat_table& src = rt_registry<D>::inv_mass(l);
        mat_table t = detail::rt_table_access::make_mat(
            src.rows(), src.cols(), detail::rt_conv_vec_mat<T>(src));
        return s.invm.insert(std::make_pair(l, std::move(t))).first->second;
    }

private:
    struct state {
        std::mutex mtx;
        std::map<int, basis_table> basis;
        std::map<int, div_table> div;
        std::map<int, flux_table> flux;
        std::map<int, block_table> cmass;
        std::map<std::pair<int, int>, block_table> dmass;
        std::map<std::pair<int, int>, cross_table> cross;
        std::map<std::pair<int, int>, mat_table> crossc;
        std::map<int, mat_table> invm;
    };
    static state& st() {
        static state s;
        return s;
    }

    static block_table conv_block(const rt_block_table& src) {
        std::vector<std::vector<T> > blk;
        blk.reserve(static_cast<std::size_t>(src.num_blocks()));
        // convert the STORED blocks (block(1,0) stays a transposed view)
        for (int b = 0; b < src.num_blocks(); ++b) {
            // stored order: single block, or (0,0), (0,1), (1,1)
            rt_block_view<detail::rational> V =
                src.num_blocks() == 1
                    ? rt_block_view<detail::rational>(&src.at(0, 0),
                                                      src.rows(), src.cols(), false)
                    : src.block(b == 0 ? 0 : (b == 1 ? 0 : 1),
                                b == 0 ? 0 : 1);
            std::vector<T> w;
            w.reserve(static_cast<std::size_t>(V.rows())
                      * static_cast<std::size_t>(V.cols()));
            for (int i = 0; i < V.rows(); ++i)
                for (int j = 0; j < V.cols(); ++j)
                    w.push_back(detail::rt_conv<T>(V.at(i, j)));
            blk.push_back(std::move(w));
        }
        return detail::rt_table_access::make_block(src.block_rows(),
                                                   src.block_cols(),
                                                   std::move(blk));
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_TYPED_TABLES_HPP
