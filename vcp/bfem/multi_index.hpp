// vcp/bfem/multi_index.hpp
// Layer 0: multi_index<D> and index_map<D> (U1: the single index authority).
//
// Conforms to: L0 external design v0.3 (sections 3.1, 3.2) and
//              L0 internal design v0.3 (section 6).
//
// Canonical order (normative): descending in alpha_0, ties broken descending
// in alpha_1, and so on (alpha_D is dependent). Example D=2, n=2:
//   (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2)
//
// rank/unrank at runtime touch only an int table baked at construction time
// (no bigint, no binomial_cache, no allocation): S-L0-1 structural check.

#ifndef VCP_BFEM_MULTI_INDEX_HPP
#define VCP_BFEM_MULTI_INDEX_HPP

#include <array>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/rational.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// multi_index<D>: multi index alpha = (alpha_0, ..., alpha_D) over the
// barycentric coordinates of a D dimensional simplex.
// ---------------------------------------------------------------------------
template <int D>
struct multi_index {
    std::array<int, D + 1> a;

    int degree() const {
        int s = 0;
        for (int i = 0; i <= D; ++i) s += a[static_cast<std::size_t>(i)];
        return s;
    }
    bool operator==(const multi_index& b) const {
        for (int i = 0; i <= D; ++i)
            if (a[static_cast<std::size_t>(i)] != b.a[static_cast<std::size_t>(i)])
                return false;
        return true;
    }
    bool operator!=(const multi_index& b) const { return !(*this == b); }
    // componentwise alpha_i <= b_i
    bool leq(const multi_index& b) const {
        for (int i = 0; i <= D; ++i)
            if (a[static_cast<std::size_t>(i)] > b.a[static_cast<std::size_t>(i)])
                return false;
        return true;
    }
    multi_index operator+(const multi_index& b) const {
        multi_index r;
        for (int i = 0; i <= D; ++i)
            r.a[static_cast<std::size_t>(i)] =
                a[static_cast<std::size_t>(i)] + b.a[static_cast<std::size_t>(i)];
        return r;
    }
    // alpha - e_i (for derivatives). second == false when alpha_i == 0.
    std::pair<multi_index, bool> minus_e(int i) const {
        assert(i >= 0 && i <= D);
        multi_index r = *this;
        if (r.a[static_cast<std::size_t>(i)] == 0)
            return std::pair<multi_index, bool>(r, false);
        r.a[static_cast<std::size_t>(i)] -= 1;
        return std::pair<multi_index, bool>(r, true);
    }
};

// ---------------------------------------------------------------------------
// index_map<D>: bijection between { alpha : |alpha| = n } and [0, N) with
// N = C(n + D, D), in the canonical order above.
//
// Construction bakes tab[d][s] = C(s + d, d) (d = 0..D, s = 0..n) into an int
// table via binomial_cache; a bigint value that does not fit int raises
// std::overflow_error (such degrees are outside practical use).
// rank/unrank afterwards are pure int table arithmetic, allocation free.
// ---------------------------------------------------------------------------
template <int D>
class index_map {
public:
    explicit index_map(int n) : n_(n), size_(0), tab_() {
        if (n < 0)
            throw std::invalid_argument("bfem::index_map: negative degree");
        tab_.resize(static_cast<std::size_t>(D + 1) * static_cast<std::size_t>(n + 1));
        for (int d = 0; d <= D; ++d) {
            for (int s = 0; s <= n; ++s) {
                const detail::bigint& v = detail::binomial_cache::binom(s + d, d);
                if (!v.fits_int64())
                    throw std::overflow_error("bfem::index_map: table exceeds int");
                long long x = v.to_int64();
                if (x > 2147483647LL)
                    throw std::overflow_error("bfem::index_map: table exceeds int");
                tab_[idx(d, s)] = static_cast<int>(x);
            }
        }
        size_ = tab_[idx(D, n)];
    }

    int size() const { return size_; }     // N(D, n)
    int degree() const { return n_; }

    // O(D) table lookups; no allocation, no bigint (S-L0-1)
    int rank(const multi_index<D>& alpha) const {
        assert(alpha.degree() == n_);
        int rem = n_;
        int r = 0;
        for (int k = 0; k < D; ++k) {
            int ak = alpha.a[static_cast<std::size_t>(k)];
            assert(ak >= 0 && ak <= rem);
            // number of predecessors at position k:
            //   sum_{v > ak} C(rem - v + d - 1, d - 1) = C(rem - ak - 1 + d, d)
            int d = D - k;
            if (rem - ak - 1 >= 0) r += tab_[idx(d, rem - ak - 1)];
            rem -= ak;
        }
        return r;
    }

    // O(D * n) table lookups; no allocation (multi_index is a std::array)
    multi_index<D> unrank(int r) const {
        assert(r >= 0 && r < size_);
        multi_index<D> alpha;
        int rem = n_;
        for (int k = 0; k < D; ++k) {
            int d = D - k;
            int v = rem;
            for (;; --v) {
                assert(v >= 0);
                // count of indices with value v at position k:
                //   C(rem - v + d - 1, d - 1) = tab[d-1][rem - v]
                int cnt = tab_[idx(d - 1, rem - v)];
                if (r < cnt) break;
                r -= cnt;
            }
            alpha.a[static_cast<std::size_t>(k)] = v;
            rem -= v;
        }
        alpha.a[static_cast<std::size_t>(D)] = rem;
        return alpha;
    }

private:
    int n_;
    int size_;
    std::vector<int> tab_;   // immutable after construction

    std::size_t idx(int d, int s) const {
        return static_cast<std::size_t>(d) * static_cast<std::size_t>(n_ + 1)
               + static_cast<std::size_t>(s);
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MULTI_INDEX_HPP
