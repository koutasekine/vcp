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
#include <vcp/bfem/detail/table_cache.hpp>

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
    static_assert(D == 2 || D == 3,
                  "bfem::typed_rt_registry: only D == 2 or D == 3");
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
        return s.basis.get_or_build(k, [&]() -> basis_table {
            const rt_basis_table& src = rt_registry<D>::basis(k);
            std::vector<T> v;
            v.reserve(static_cast<std::size_t>(src.rows())
                      * static_cast<std::size_t>(src.cols()));
            for (int i = 0; i < src.rows(); ++i)
                for (int j = 0; j < src.cols(); ++j)
                    v.push_back(detail::rt_conv<T>(src.at(i, j)));
            basis_table t = detail::rt_table_access::make_basis(
                src.num_comp(), src.order(), src.comp_size(), src.dim(), std::move(v));
            return t;
        });
    }

    static const div_table& divergence(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.div.get_or_build(k, [&]() -> div_table {
            const rt_div_table& src = rt_registry<D>::divergence(k);
            std::vector<T> v;
            v.reserve(static_cast<std::size_t>(src.rows())
                      * static_cast<std::size_t>(src.cols()));
            for (int i = 0; i < src.rows(); ++i)
                for (int j = 0; j < src.cols(); ++j)
                    v.push_back(detail::rt_conv<T>(src.at(i, j)));
            div_table t = detail::rt_table_access::make_div(
                src.order(), src.rows(), src.cols(), std::move(v));
            return t;
        });
    }

    static const flux_table& edge_flux(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.flux.get_or_build(k, [&]() -> flux_table {
            const rt_flux_table& src = rt_registry<D>::edge_flux(k);
            const int kk = src.order();
            const int dm = src.dim();
            const int nfac = src.num_facets();
            const int pf = src.per_facet();
            std::vector<T> v;
            v.reserve(static_cast<std::size_t>(nfac) * static_cast<std::size_t>(pf)
                      * static_cast<std::size_t>(dm));
            for (int e = 0; e < nfac; ++e)
                for (int j = 0; j < pf; ++j)
                    for (int c = 0; c < dm; ++c)
                        v.push_back(detail::rt_conv<T>(src.at(e, j, c)));
            flux_table t = detail::rt_table_access::make_flux(kk, dm, nfac, pf,
                                                              std::move(v));
            return t;
        });
    }

    static const block_table& comp_mass(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.cmass.get_or_build(k, [&]() -> block_table {
            return conv_block(rt_registry<D>::comp_mass(k));
        });
    }

    static const block_table& div_mass(int k, int l) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, l);
        return s.dmass.get_or_build(key, [&]() -> block_table {
            return conv_block(rt_registry<D>::div_mass(k, l));
        });
    }

    static const cross_table& cross_grad(int k, int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, n);
        return s.cross.get_or_build(key, [&]() -> cross_table {
            const rt_cross_table& src = rt_registry<D>::cross_grad(k, n);
            const int nc = src.n_comp();
            const int nv = src.n_vert();
            std::vector<std::vector<T> > blk;
            blk.reserve(static_cast<std::size_t>(nc) * static_cast<std::size_t>(nv));
            for (int d = 0; d < nc; ++d) {
                for (int i = 0; i < nv; ++i) {
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
                src.dim(), src.nn(), nv, std::move(blk));
            return t;
        });
    }

    static const mat_table& cross_grad_contracted(int k, int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, n);
        return s.crossc.get_or_build(key, [&]() -> mat_table {
            const rt_mat_table& src = rt_registry<D>::cross_grad_contracted(k, n);
            mat_table t = detail::rt_table_access::make_mat(
                src.rows(), src.cols(), detail::rt_conv_vec_mat<T>(src));
            return t;
        });
    }

    static const mat_table& inv_mass(int l) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.invm.get_or_build(l, [&]() -> mat_table {
            const rt_mat_table& src = rt_registry<D>::inv_mass(l);
            mat_table t = detail::rt_table_access::make_mat(
                src.rows(), src.cols(), detail::rt_conv_vec_mat<T>(src));
            return t;
        });
    }

private:
    // L4: shared detail::table_cache vessel (one mutex per D x T, unchanged)
    struct maps {
        detail::table_cache<int, basis_table> basis;
        detail::table_cache<int, div_table> div;
        detail::table_cache<int, flux_table> flux;
        detail::table_cache<int, block_table> cmass;
        detail::table_cache<std::pair<int, int>, block_table> dmass;
        detail::table_cache<std::pair<int, int>, cross_table> cross;
        detail::table_cache<std::pair<int, int>, mat_table> crossc;
        detail::table_cache<int, mat_table> invm;
    };
    typedef detail::table_cache_state<maps> state;
    static state& st() {
        return detail::table_cache_instance<state>();
    }

    static block_table conv_block(const rt_block_table& src) {
        std::vector<std::vector<T> > blk;
        blk.reserve(static_cast<std::size_t>(src.num_blocks()));
        // convert the STORED blocks only (the lower triangle stays a
        // transposed view); stored order: single block, or the upper
        // triangle (d <= d') in row-major order
        const int nblk = src.num_blocks();
        int d = 0, dp = 0;
        for (int b = 0; b < nblk; ++b) {
            rt_block_view<detail::rational> V =
                nblk == 1
                    ? rt_block_view<detail::rational>(&src.at(0, 0),
                                                      src.rows(), src.cols(), false)
                    : src.block(d, dp);
            std::vector<T> w;
            w.reserve(static_cast<std::size_t>(V.rows())
                      * static_cast<std::size_t>(V.cols()));
            for (int i = 0; i < V.rows(); ++i)
                for (int j = 0; j < V.cols(); ++j)
                    w.push_back(detail::rt_conv<T>(V.at(i, j)));
            blk.push_back(std::move(w));
            if (nblk > 1) {                       // next upper-triangle pair
                ++dp;
                if (dp >= src.num_comp()) { ++d; dp = d; }
            }
        }
        return detail::rt_table_access::make_block(src.block_rows(),
                                                   src.block_cols(),
                                                   std::move(blk));
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_TYPED_TABLES_HPP
