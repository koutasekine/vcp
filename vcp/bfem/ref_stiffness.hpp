// vcp/bfem/ref_stiffness.hpp
// Layer 2: reference stiffness tensor R^{(ij)} (V3) and its typed cache.
//
// Conforms to: L2 external design v0.3 (section 5.1) and
//              L2 internal design v0.3 (section 4).
//
//   R^{(ij)}_{ab} = n^2 * M^{n-1,n-1}[dm(a,i), dm(b,j)]   (0 on vanishing)
//
// Only the frozen L0 API (typed mass + derivative_map) is used. The cache
// vessel is the shared detail::table_cache (the post-stabilization
// unification announced in internal design 10.3, executed as L4 seam
// alignment section 10-1; behavior and lock granularity unchanged).

#ifndef VCP_BFEM_REF_STIFFNESS_HPP
#define VCP_BFEM_REF_STIFFNESS_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/typed_tables.hpp>
#include <vcp/bfem/bpoly.hpp>   // detail::deriv_cache
#include <vcp/bfem/detail/table_cache.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// upper-triangle block index for (i, j), i <= j, over 0..D
inline int ref_stiff_block_index(int D, int i, int j) {
    assert(0 <= i && i <= j && j <= D);
    return i * (D + 1) - i * (i - 1) / 2 + (j - i);
}

template <int D, typename T>
class ref_stiffness_cache {
public:
    struct tensor {
        int n;
        int N;                                  // N(D, n)
        // blocks for i <= j (packed upper triangle), each N*N row-major
        std::vector<std::vector<T> > blocks;
        const T& at(int i, int j, int a, int b) const {   // requires i <= j
            return blocks[static_cast<std::size_t>(ref_stiff_block_index(D, i, j))]
                         [static_cast<std::size_t>(a) * static_cast<std::size_t>(N)
                          + static_cast<std::size_t>(b)];
        }
    };

    // n >= 1 (enforced by the element_op entry checks, C-3)
    static const tensor& get(int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.tensors.get_or_build(n, [&]() -> tensor {
            return build(n);
        });
    }

private:
    // L4: shared detail::table_cache vessel (one mutex per D x T, unchanged)
    struct maps {
        table_cache<int, tensor> tensors;
    };
    typedef table_cache_state<maps> state;
    static state& st() {
        return table_cache_instance<state>();
    }

    static tensor build(int n) {
        assert(n >= 1);
        const typed_mass_table<D, T>& M = typed_registry<D, T>::mass(n - 1, n - 1);
        const derivative_map<D>& dm = deriv_cache<D>::get(n);
        const int N = coeff_registry<D>::indices(n).size();
        const T n2 = T(n) * T(n);
        tensor t;
        t.n = n;
        t.N = N;
        t.blocks.resize(static_cast<std::size_t>((D + 1) * (D + 2) / 2));
        for (int i = 0; i <= D; ++i) {
            for (int j = i; j <= D; ++j) {
                std::vector<T>& blk =
                    t.blocks[static_cast<std::size_t>(ref_stiff_block_index(D, i, j))];
                blk.assign(static_cast<std::size_t>(N) * static_cast<std::size_t>(N), T(0));
                for (int a = 0; a < N; ++a) {
                    int ta = dm.target(a, i);
                    if (ta < 0) continue;
                    for (int b = 0; b < N; ++b) {
                        int tb = dm.target(b, j);
                        if (tb < 0) continue;
                        blk[static_cast<std::size_t>(a) * static_cast<std::size_t>(N)
                            + static_cast<std::size_t>(b)] = n2 * M.at(ta, tb);
                    }
                }
            }
        }
        return t;
    }
};

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_REF_STIFFNESS_HPP
