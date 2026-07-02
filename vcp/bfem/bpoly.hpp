// vcp/bfem/bpoly.hpp
// Layer 1: local polynomial algebra bpoly<D,T> (F1-F5, F7, F8).
//
// Conforms to: L1 external design v0.3 (sections 3.1, 4, 5, 7) and
//              L1 internal design v0.3 (sections 2, 3, 4, 6).
//
// Every algorithm here is "L0 table reads + ring operations of T".
// No floating point literal, no division, no numerical integration.

#ifndef VCP_BFEM_BPOLY_HPP
#define VCP_BFEM_BPOLY_HPP

#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/typed_tables.hpp>

namespace vcp {
namespace bfem {

// barycentric point (sum of coordinates == 1 is the caller's responsibility)
template <int D, typename T>
using bary_point = std::array<T, D + 1>;

template <int D, typename T> class bpoly;

namespace detail {
// access key for the internal prepare/coefficient plumbing (bpoly_core)
struct bpoly_access;

// B-4: L1-owned static cache of derivative_map<D> keyed by degree
// (same magic static + mutex + std::map pattern as the L0 registry).
template <int D>
class deriv_cache {
public:
    static const derivative_map<D>& get(int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, derivative_map<D> >::iterator it = s.maps.find(n);
        if (it != s.maps.end()) return it->second;
        return s.maps.insert(std::make_pair(n, derivative_map<D>(n))).first->second;
    }
private:
    struct state {
        std::mutex mtx;
        std::map<int, derivative_map<D> > maps;
    };
    static state& st() {
        static state s;
        return s;
    }
};
} // namespace detail

// ---------------------------------------------------------------------------
// bpoly<D, T>: u in P^n on the reference simplex, stored as degree + the
// coefficient vector in the L0 canonical order.
// Invariant: c_.size() == N(D, deg_); degree 0 is the valid constant state
// (there is no "empty" state).
// ---------------------------------------------------------------------------
template <int D, typename T>
class bpoly {
public:
    // ---- construction ----
    bpoly() : deg_(0), c_(1) {}                      // degree 0 zero polynomial

    static bpoly constant(const T& s) {
        bpoly r;
        r.c_[0] = s;
        return r;
    }
    static bpoly zero(int n) {
        bpoly r;
        r.deg_ = check_degree(n);
        r.c_.assign(static_cast<std::size_t>(coeff_registry<D>::indices(n).size()), T(0));
        return r;
    }
    static bpoly from_coeffs(int n, std::vector<T> c) {
        check_degree(n);
        if (static_cast<int>(c.size()) != coeff_registry<D>::indices(n).size())
            throw std::invalid_argument("bfem::bpoly::from_coeffs: size != N(D, n)");
        bpoly r;
        r.deg_ = n;
        r.c_ = std::move(c);
        return r;
    }

    // ---- observers ----
    int degree() const { return deg_; }              // stored degree (U6)
    int size() const { return static_cast<int>(c_.size()); }
    const T& coeff(int rank) const {
        assert(rank >= 0 && rank < size());
        return c_[static_cast<std::size_t>(rank)];
    }
    T& coeff(int rank) {                             // degree stays invariant
        assert(rank >= 0 && rank < size());
        return c_[static_cast<std::size_t>(rank)];
    }
    const std::vector<T>& coeffs() const { return c_; }

    // ---- in-place (degree preserving only, U1) ----
    bpoly& operator+=(const bpoly& v) {
        if (v.deg_ != deg_)
            throw std::invalid_argument("bfem::bpoly::operator+=: degree mismatch");
        for (std::size_t i = 0; i < c_.size(); ++i) c_[i] += v.c_[i];
        return *this;
    }
    bpoly& operator-=(const bpoly& v) {
        if (v.deg_ != deg_)
            throw std::invalid_argument("bfem::bpoly::operator-=: degree mismatch");
        for (std::size_t i = 0; i < c_.size(); ++i) c_[i] -= v.c_[i];
        return *this;
    }
    bpoly& operator*=(const T& s) {
        for (std::size_t i = 0; i < c_.size(); ++i) c_[i] *= s;
        return *this;
    }
    bpoly& add_scalar(const T& s) {                  // u + s via partition of unity
        for (std::size_t i = 0; i < c_.size(); ++i) c_[i] += s;
        return *this;
    }

private:
    int deg_;
    std::vector<T> c_;

    static int check_degree(int n) {
        if (n < 0)
            throw std::invalid_argument("bfem::bpoly: negative degree");
        return n;
    }
    friend struct detail::bpoly_access;
};

namespace detail {

// bpoly_core plumbing: prepare(dst, n, zero_fill) with the reuse contract of
// internal design 4.1 (resize only when the size changes; zero by assignment,
// never by reconstruction, so mantissa memory of mpfr-like T is reused).
struct bpoly_access {
    template <int D, typename T>
    static void prepare(bpoly<D, T>& dst, int n, bool zero_fill) {
        dst.deg_ = n;
        std::size_t N = static_cast<std::size_t>(coeff_registry<D>::indices(n).size());
        if (dst.c_.size() != N) dst.c_.resize(N);
        if (zero_fill) {
            const T z(0);
            for (std::size_t i = 0; i < dst.c_.size(); ++i) dst.c_[i] = z;
        }
    }
    template <int D, typename T>
    static std::vector<T>& vec(bpoly<D, T>& u) { return u.c_; }
    template <int D, typename T>
    static const std::vector<T>& vec(const bpoly<D, T>& u) { return u.c_; }
};

} // namespace detail

// ---------------------------------------------------------------------------
// F1: add / sub. Result degree max(n, m); the lower-degree side is lifted by
// the L0 elevation scatter in one pass (no intermediate bpoly). Identity
// elevation is skipped by degree comparison (U3; no L0 flag API involved).
// Aliasing dst == u or dst == v is forbidden (debug assert).
// ---------------------------------------------------------------------------
template <int D, typename T>
void add_into(bpoly<D, T>& dst, const bpoly<D, T>& u, const bpoly<D, T>& v) {
    assert(&dst != &u && &dst != &v);
    const bpoly<D, T>& hi = (u.degree() >= v.degree()) ? u : v;
    const bpoly<D, T>& lo = (u.degree() >= v.degree()) ? v : u;
    detail::bpoly_access::prepare(dst, hi.degree(), false);
    std::vector<T>& d = detail::bpoly_access::vec(dst);
    const std::vector<T>& h = detail::bpoly_access::vec(hi);
    const std::vector<T>& l = detail::bpoly_access::vec(lo);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = h[i];
    if (lo.degree() == hi.degree()) {
        for (std::size_t i = 0; i < d.size(); ++i) d[i] += l[i];
        return;
    }
    const typed_elevation_table<D, T>& E =
        typed_registry<D, T>::elevation(lo.degree(), hi.degree());
    for (int i = 0; i < E.source_size(); ++i) {
        typename typed_elevation_table<D, T>::entry_range r = E.row(i);
        const T& li = l[static_cast<std::size_t>(i)];
        for (typename typed_elevation_table<D, T>::entry_iterator p = r.begin();
             p != r.end(); ++p) {
            typename typed_elevation_table<D, T>::entry e = *p;
            d[static_cast<std::size_t>(e.target_rank)] += *e.coeff * li;
        }
    }
}

template <int D, typename T>
void sub_into(bpoly<D, T>& dst, const bpoly<D, T>& u, const bpoly<D, T>& v) {
    assert(&dst != &u && &dst != &v);
    if (u.degree() >= v.degree()) {
        // copy u, then subtract (elevated) v
        detail::bpoly_access::prepare(dst, u.degree(), false);
        std::vector<T>& d = detail::bpoly_access::vec(dst);
        const std::vector<T>& a = detail::bpoly_access::vec(u);
        const std::vector<T>& b = detail::bpoly_access::vec(v);
        for (std::size_t i = 0; i < d.size(); ++i) d[i] = a[i];
        if (v.degree() == u.degree()) {
            for (std::size_t i = 0; i < d.size(); ++i) d[i] -= b[i];
            return;
        }
        const typed_elevation_table<D, T>& E =
            typed_registry<D, T>::elevation(v.degree(), u.degree());
        for (int i = 0; i < E.source_size(); ++i) {
            typename typed_elevation_table<D, T>::entry_range r = E.row(i);
            const T& bi = b[static_cast<std::size_t>(i)];
            for (typename typed_elevation_table<D, T>::entry_iterator p = r.begin();
                 p != r.end(); ++p) {
                typename typed_elevation_table<D, T>::entry e = *p;
                d[static_cast<std::size_t>(e.target_rank)] -= *e.coeff * bi;
            }
        }
    } else {
        // store -v while copying (one pass, internal design 3.1), then add u
        detail::bpoly_access::prepare(dst, v.degree(), false);
        std::vector<T>& d = detail::bpoly_access::vec(dst);
        const std::vector<T>& a = detail::bpoly_access::vec(u);
        const std::vector<T>& b = detail::bpoly_access::vec(v);
        for (std::size_t i = 0; i < d.size(); ++i) d[i] = -b[i];
        const typed_elevation_table<D, T>& E =
            typed_registry<D, T>::elevation(u.degree(), v.degree());
        for (int i = 0; i < E.source_size(); ++i) {
            typename typed_elevation_table<D, T>::entry_range r = E.row(i);
            const T& ai = a[static_cast<std::size_t>(i)];
            for (typename typed_elevation_table<D, T>::entry_iterator p = r.begin();
                 p != r.end(); ++p) {
                typename typed_elevation_table<D, T>::entry e = *p;
                d[static_cast<std::size_t>(e.target_rank)] += *e.coeff * ai;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// F2: mul. Result degree a + b. The hot loop contains only the precomputed
// target lookup + coefficient multiplication (no rank computation, no
// allocation): structural check S-L1-1. Exactly 2 T-multiplications per
// coefficient pair (2 * N_a * N_b in total).
// ---------------------------------------------------------------------------
template <int D, typename T>
void mul_into(bpoly<D, T>& dst, const bpoly<D, T>& u, const bpoly<D, T>& v) {
    assert(&dst != &u && &dst != &v);
    int a = u.degree(), b = v.degree();
    const typed_product_table<D, T>& P = typed_registry<D, T>::product(a, b);
    detail::bpoly_access::prepare(dst, a + b, true);
    std::vector<T>& d = detail::bpoly_access::vec(dst);
    const std::vector<T>& ua = detail::bpoly_access::vec(u);
    const std::vector<T>& vb = detail::bpoly_access::vec(v);
    const int Na = P.rows(), Nb = P.cols();
    for (int i = 0; i < Na; ++i) {
        const T& ui = ua[static_cast<std::size_t>(i)];
        for (int j = 0; j < Nb; ++j) {
            d[static_cast<std::size_t>(P.target_rank(i, j))] +=
                P.coeff(i, j) * ui * vb[static_cast<std::size_t>(j)];
        }
    }
}

// ---------------------------------------------------------------------------
// F3: scale
// ---------------------------------------------------------------------------
template <int D, typename T>
void scale_into(bpoly<D, T>& dst, const bpoly<D, T>& u, const T& s) {
    assert(&dst != &u);
    detail::bpoly_access::prepare(dst, u.degree(), false);
    std::vector<T>& d = detail::bpoly_access::vec(dst);
    const std::vector<T>& a = detail::bpoly_access::vec(u);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = a[i] * s;
}

// ---------------------------------------------------------------------------
// F4: elevate. Precondition m >= u.degree(); m == degree copies (identity by
// degree comparison, U3 -- the L0 table is not touched in that case).
// ---------------------------------------------------------------------------
template <int D, typename T>
void elevate_into(bpoly<D, T>& dst, const bpoly<D, T>& u, int m) {
    assert(&dst != &u);
    if (m < u.degree())
        throw std::invalid_argument("bfem::elevate: m < degree(u)");
    if (m == u.degree()) {
        detail::bpoly_access::prepare(dst, m, false);
        std::vector<T>& d = detail::bpoly_access::vec(dst);
        const std::vector<T>& a = detail::bpoly_access::vec(u);
        for (std::size_t i = 0; i < d.size(); ++i) d[i] = a[i];
        return;
    }
    const typed_elevation_table<D, T>& E =
        typed_registry<D, T>::elevation(u.degree(), m);
    detail::bpoly_access::prepare(dst, m, true);
    std::vector<T>& d = detail::bpoly_access::vec(dst);
    const std::vector<T>& a = detail::bpoly_access::vec(u);
    for (int i = 0; i < E.source_size(); ++i) {
        typename typed_elevation_table<D, T>::entry_range r = E.row(i);
        const T& ai = a[static_cast<std::size_t>(i)];
        for (typename typed_elevation_table<D, T>::entry_iterator p = r.begin();
             p != r.end(); ++p) {
            typename typed_elevation_table<D, T>::entry e = *p;
            d[static_cast<std::size_t>(e.target_rank)] += *e.coeff * ai;
        }
    }
}

// ---------------------------------------------------------------------------
// F5: dlambda. Degree 0 input yields the degree 0 zero. The integer factor n
// is exact in T; the index shift comes from the L1-cached derivative_map.
// ---------------------------------------------------------------------------
template <int D, typename T>
void dlambda_into(bpoly<D, T>& dst, const bpoly<D, T>& u, int i) {
    assert(&dst != &u);
    if (i < 0 || i > D)
        throw std::invalid_argument("bfem::dlambda: direction out of range");
    int n = u.degree();
    if (n == 0) {
        detail::bpoly_access::prepare(dst, 0, true);
        return;
    }
    const derivative_map<D>& dm = detail::deriv_cache<D>::get(n);
    detail::bpoly_access::prepare(dst, n - 1, true);
    std::vector<T>& d = detail::bpoly_access::vec(dst);
    const std::vector<T>& a = detail::bpoly_access::vec(u);
    const T fn(n);
    for (int r = 0; r < dm.source_size(); ++r) {
        int t = dm.target(r, i);
        if (t >= 0) d[static_cast<std::size_t>(t)] = fn * a[static_cast<std::size_t>(r)];
    }
}

// ---------------------------------------------------------------------------
// pure-function forms (the default style, U1)
// ---------------------------------------------------------------------------
template <int D, typename T>
bpoly<D, T> add(const bpoly<D, T>& u, const bpoly<D, T>& v) {
    bpoly<D, T> r;
    add_into(r, u, v);
    return r;
}
template <int D, typename T>
bpoly<D, T> sub(const bpoly<D, T>& u, const bpoly<D, T>& v) {
    bpoly<D, T> r;
    sub_into(r, u, v);
    return r;
}
template <int D, typename T>
bpoly<D, T> mul(const bpoly<D, T>& u, const bpoly<D, T>& v) {
    bpoly<D, T> r;
    mul_into(r, u, v);
    return r;
}
template <int D, typename T>
bpoly<D, T> scale(const bpoly<D, T>& u, const T& s) {
    bpoly<D, T> r;
    scale_into(r, u, s);
    return r;
}
template <int D, typename T>
bpoly<D, T> elevate(const bpoly<D, T>& u, int m) {
    bpoly<D, T> r;
    elevate_into(r, u, m);
    return r;
}
template <int D, typename T>
bpoly<D, T> dlambda(const bpoly<D, T>& u, int i) {
    bpoly<D, T> r;
    dlambda_into(r, u, i);
    return r;
}

// operator sugar (same semantics as the named functions)
template <int D, typename T>
bpoly<D, T> operator+(const bpoly<D, T>& u, const bpoly<D, T>& v) { return add(u, v); }
template <int D, typename T>
bpoly<D, T> operator-(const bpoly<D, T>& u, const bpoly<D, T>& v) { return sub(u, v); }
template <int D, typename T>
bpoly<D, T> operator*(const bpoly<D, T>& u, const bpoly<D, T>& v) { return mul(u, v); }
template <int D, typename T>
bpoly<D, T> operator*(const bpoly<D, T>& u, const T& s) { return scale(u, s); }
template <int D, typename T>
bpoly<D, T> operator*(const T& s, const bpoly<D, T>& u) { return scale(u, s); }

// ---------------------------------------------------------------------------
// F7: reference-element L2 inner product (coefficient of |T|; the physical
// value is Layer 2's measure multiplication). c_u^T M^{a,b} c_v with the row
// accumulation of internal design 3.3 (N_a*N_b + N_a multiplications).
// ---------------------------------------------------------------------------
template <int D, typename T>
T inner(const bpoly<D, T>& u, const bpoly<D, T>& v) {
    const typed_mass_table<D, T>& M =
        typed_registry<D, T>::mass(u.degree(), v.degree());
    const std::vector<T>& a = detail::bpoly_access::vec(u);
    const std::vector<T>& b = detail::bpoly_access::vec(v);
    T acc(0);
    for (int i = 0; i < M.rows(); ++i) {
        T row(0);
        for (int j = 0; j < M.cols(); ++j)
            row += M.at(i, j) * b[static_cast<std::size_t>(j)];
        acc += a[static_cast<std::size_t>(i)] * row;
    }
    return acc;
}

// ---------------------------------------------------------------------------
// F8: point evaluation (de Casteljau; additions and multiplications only, so
// interval T yields an enclosure). Off the integral path by design (U8).
// ---------------------------------------------------------------------------
template <int D, typename T>
T eval(const bpoly<D, T>& u, const bary_point<D, T>& lam) {
    int n = u.degree();
    std::vector<T> work(detail::bpoly_access::vec(u));
    for (int k = n; k >= 1; --k) {
        const index_map<D>& im_hi = coeff_registry<D>::indices(k);
        const index_map<D>& im_lo = coeff_registry<D>::indices(k - 1);
        for (int r = 0; r < im_lo.size(); ++r) {
            multi_index<D> alpha = im_lo.unrank(r);
            T acc(0);
            for (int i = 0; i <= D; ++i) {
                multi_index<D> ai = alpha;
                ai.a[static_cast<std::size_t>(i)] += 1;
                acc += lam[static_cast<std::size_t>(i)]
                       * work[static_cast<std::size_t>(im_hi.rank(ai))];
            }
            work[static_cast<std::size_t>(r)] = acc;   // forward packing
        }
    }
    return work[0];
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_BPOLY_HPP
