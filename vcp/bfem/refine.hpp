// vcp/bfem/refine.hpp
// Layer 1: restriction to sub-simplices (F9, blossom) and rigorous range
// enclosure (F10: range / range_refined) with bound_traits.
//
// Conforms to: L1 external design v0.3 (section 8) and
//              L1 internal design v0.3 (sections 7, 8).

#ifndef VCP_BFEM_REFINE_HPP
#define VCP_BFEM_REFINE_HPP

#include <array>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <cassert>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/convert_traits.hpp>

namespace vcp {
namespace bfem {

template <typename T>
struct range_pair {
    T lo;
    T hi;
};

// ---------------------------------------------------------------------------
// bound_traits (C-2): the primary template is valid for builtin arithmetic
// types only; any other T requires an explicit specialization (point types
// can simply inherit point_bound_traits<T>; interval types need an enclosing
// specialization, see bound_traits_kv.hpp). Forgetting the specialization is
// a compile error, never a silently unsound result.
// ---------------------------------------------------------------------------
template <typename T>
struct point_bound_traits {
    static void min_update(T& acc, const T& x) { if (x < acc) acc = x; }
    static void max_update(T& acc, const T& x) { if (acc < x) acc = x; }
};

template <typename T>
struct bound_traits {
    static_assert(std::is_arithmetic<T>::value,
        "bound_traits<T>: specialize for non-arithmetic T "
        "(point types may inherit point_bound_traits<T>; "
        "interval types need an enclosing specialization, e.g. bound_traits_kv.hpp)");
    static void min_update(T& acc, const T& x) { point_bound_traits<T>::min_update(acc, x); }
    static void max_update(T& acc, const T& x) { point_bound_traits<T>::max_update(acc, x); }
};

// ---------------------------------------------------------------------------
// F10: coefficient bound (convex hull property). For every x in the reference
// simplex, lo <= u(x) <= hi in the sense of bound_traits<T>; this is a range
// ENCLOSURE, not the min/max themselves.
// ---------------------------------------------------------------------------
namespace detail {

template <int D, typename T>
void range_merge(range_pair<T>& acc, bool& started, const bpoly<D, T>& u) {
    const std::vector<T>& c = bpoly_access::vec(u);
    std::size_t i0 = 0;
    if (!started) {
        acc.lo = c[0];
        acc.hi = c[0];
        started = true;
        i0 = 1;
    }
    for (std::size_t i = i0; i < c.size(); ++i) {
        bound_traits<T>::min_update(acc.lo, c[i]);
        bound_traits<T>::max_update(acc.hi, c[i]);
    }
}

} // namespace detail

template <int D, typename T>
range_pair<T> range(const bpoly<D, T>& u) {
    range_pair<T> acc;
    bool started = false;
    detail::range_merge(acc, started, u);
    return acc;
}

// ---------------------------------------------------------------------------
// F9: restriction to the sub-simplex spanned by V = (v_0, ..., v_D) given in
// the parent's barycentric coordinates. Exact coefficient generation via the
// blossom (iterated de Casteljau); additions/multiplications only, so
// interval T automatically yields enclosures.
//
// Implementation: prefix-sharing depth-first traversal (internal design 7).
// dfs(j, k): buf[j] holds the blossom intermediate of degree k after
// consuming beta_0 copies of w_0, ..., beta_{j-1} copies of w_{j-1} and t
// copies of w_j so far; siblings share the reductions of common prefixes.
// ---------------------------------------------------------------------------
namespace detail {

// one de Casteljau reduction step (degree k -> k-1) in direction w,
// in place with forward packing (write rank in degree k-1 is never above any
// read rank in degree k: rank_k(alpha + e_i) >= rank_{k-1}(alpha))
template <int D, typename T>
void decasteljau_step(std::vector<T>& A, int k, const bary_point<D, T>& w) {
    const index_map<D>& hi = coeff_registry<D>::indices(k);
    const index_map<D>& lo = coeff_registry<D>::indices(k - 1);
    for (int r = 0; r < lo.size(); ++r) {
        multi_index<D> alpha = lo.unrank(r);
        T acc(0);
        for (int i = 0; i <= D; ++i) {
            multi_index<D> ai = alpha;
            ai.a[static_cast<std::size_t>(i)] += 1;
            acc += w[static_cast<std::size_t>(i)]
                   * A[static_cast<std::size_t>(hi.rank(ai))];
        }
        A[static_cast<std::size_t>(r)] = acc;
    }
}

template <int D, typename T>
struct blossom_dfs {
    const std::array<bary_point<D, T>, D + 1>* V;
    const index_map<D>* im_n;
    std::array<std::vector<T>, static_cast<std::size_t>(D + 1)> buf;
    std::array<int, static_cast<std::size_t>(D + 1)> beta;
    std::vector<T>* out;

    void run(int j, int k) {
        if (j == D) {
            beta[static_cast<std::size_t>(D)] = k;
            std::vector<T>& A = buf[static_cast<std::size_t>(D)];
            for (int t = k; t >= 1; --t)
                decasteljau_step<D, T>(A, t, (*V)[static_cast<std::size_t>(D)]);
            multi_index<D> b;
            for (int i = 0; i <= D; ++i)
                b.a[static_cast<std::size_t>(i)] = beta[static_cast<std::size_t>(i)];
            (*out)[static_cast<std::size_t>(im_n->rank(b))] =
                A[0];
            return;
        }
        std::vector<T>& A = buf[static_cast<std::size_t>(j)];
        std::vector<T>& C = buf[static_cast<std::size_t>(j + 1)];
        for (int t = 0; t <= k; ++t) {
            beta[static_cast<std::size_t>(j)] = t;
            int sz = coeff_registry<D>::indices(k - t).size();
            for (int i = 0; i < sz; ++i)
                C[static_cast<std::size_t>(i)] = A[static_cast<std::size_t>(i)];
            run(j + 1, k - t);
            if (t < k)
                decasteljau_step<D, T>(A, k - t, (*V)[static_cast<std::size_t>(j)]);
        }
    }
};

} // namespace detail

template <int D, typename T>
bpoly<D, T> restrict_to(const bpoly<D, T>& u,
                        const std::array<bary_point<D, T>, D + 1>& V) {
    const int n = u.degree();
    bpoly<D, T> out = bpoly<D, T>::zero(n);
    detail::blossom_dfs<D, T> dfs;
    dfs.V = &V;
    dfs.im_n = &coeff_registry<D>::indices(n);
    dfs.out = &detail::bpoly_access::vec(out);
    std::size_t N = detail::bpoly_access::vec(u).size();
    for (int j = 0; j <= D; ++j)
        dfs.buf[static_cast<std::size_t>(j)].assign(N, T(0));
    dfs.buf[0] = detail::bpoly_access::vec(u);
    dfs.run(0, n);
    return out;
}

// ---------------------------------------------------------------------------
// F10 refinement: uniform red refinement (D == 2, four children) applied
// depth times; union of the leaf coefficient bounds. depth == 0 equals range.
// The midpoint constant 1/2 goes through convert_traits (enclose-once).
// ---------------------------------------------------------------------------
namespace detail {

template <typename T>
std::array<std::array<bary_point<2, T>, 3>, 4> red_children_2d() {
    const T o(1);
    const T z(0);
    const T h = ::vcp::bfem::rational_to<T>(1, 2);
    std::array<std::array<bary_point<2, T>, 3>, 4> kids;
    // C0: (1,0,0) (h,h,0) (h,0,h)
    kids[0][0][0] = o; kids[0][0][1] = z; kids[0][0][2] = z;
    kids[0][1][0] = h; kids[0][1][1] = h; kids[0][1][2] = z;
    kids[0][2][0] = h; kids[0][2][1] = z; kids[0][2][2] = h;
    // C1: (h,h,0) (0,1,0) (0,h,h)
    kids[1][0][0] = h; kids[1][0][1] = h; kids[1][0][2] = z;
    kids[1][1][0] = z; kids[1][1][1] = o; kids[1][1][2] = z;
    kids[1][2][0] = z; kids[1][2][1] = h; kids[1][2][2] = h;
    // C2: (h,0,h) (0,h,h) (0,0,1)
    kids[2][0][0] = h; kids[2][0][1] = z; kids[2][0][2] = h;
    kids[2][1][0] = z; kids[2][1][1] = h; kids[2][1][2] = h;
    kids[2][2][0] = z; kids[2][2][1] = z; kids[2][2][2] = o;
    // C3 (central, inverted): (0,h,h) (h,0,h) (h,h,0)
    kids[3][0][0] = z; kids[3][0][1] = h; kids[3][0][2] = h;
    kids[3][1][0] = h; kids[3][1][1] = z; kids[3][1][2] = h;
    kids[3][2][0] = h; kids[3][2][1] = h; kids[3][2][2] = z;
    return kids;
}

template <typename T>
void range_refined_rec(const bpoly<2, T>& u, int d,
                       const std::array<std::array<bary_point<2, T>, 3>, 4>& kids,
                       range_pair<T>& acc, bool& started) {
    if (d == 0) {
        range_merge(acc, started, u);
        return;
    }
    for (int c = 0; c < 4; ++c)
        range_refined_rec(restrict_to(u, kids[static_cast<std::size_t>(c)]),
                          d - 1, kids, acc, started);
}

} // namespace detail

template <typename T>
range_pair<T> range_refined(const bpoly<2, T>& u, int depth) {
    if (depth < 0)
        throw std::invalid_argument("bfem::range_refined: negative depth");
    if (depth == 0) return range(u);
    std::array<std::array<bary_point<2, T>, 3>, 4> kids =
        detail::red_children_2d<T>();
    range_pair<T> acc;
    bool started = false;
    detail::range_refined_rec(u, depth, kids, acc, started);
    return acc;
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_REFINE_HPP
