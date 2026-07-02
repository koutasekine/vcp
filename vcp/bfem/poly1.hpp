// vcp/bfem/poly1.hpp
// Layer 1: univariate polynomial poly1<T> (U4) and composition F6.
//
// Conforms to: L1 external design v0.3 (sections 3.2, 6) and
//              L1 internal design v0.3 (sections 4.3, 5).

#ifndef VCP_BFEM_POLY1_HPP
#define VCP_BFEM_POLY1_HPP

#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/convert_traits.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// poly1<T>: f(x) = sum_{k=0}^{K} a_k x^k, monomial coefficients in T.
// ---------------------------------------------------------------------------
template <typename T>
class poly1 {
public:
    static poly1 from_coeffs(std::vector<T> a) {
        if (a.empty())
            throw std::invalid_argument("bfem::poly1::from_coeffs: empty");
        poly1 r;
        r.a_ = std::move(a);
        return r;
    }

    // enclose-once construction from exact rational coefficients (the
    // recommended path when f is exactly representable): each a_k passes
    // through convert_traits<T> exactly once.
    static poly1 from_rational(
        const std::vector<std::pair<long long, long long> >& num_den) {
        if (num_den.empty())
            throw std::invalid_argument("bfem::poly1::from_rational: empty");
        std::vector<T> a;
        a.reserve(num_den.size());
        for (std::size_t k = 0; k < num_den.size(); ++k) {
            if (num_den[k].second == 0)
                throw std::invalid_argument("bfem::poly1::from_rational: zero denominator");
            a.push_back(convert_traits<T>::from_rational(
                detail::bigint(num_den[k].first), detail::bigint(num_den[k].second)));
        }
        return from_coeffs(std::move(a));
    }

    int degree() const { return static_cast<int>(a_.size()) - 1; }   // K
    const T& coeff(int k) const {
        assert(k >= 0 && k < static_cast<int>(a_.size()));
        return a_[static_cast<std::size_t>(k)];
    }

    // formal derivative f' (degree K-1; K == 0 gives the zero constant)
    poly1 derivative() const {
        poly1 r;
        if (a_.size() <= 1) {
            r.a_.assign(1, T(0));
            return r;
        }
        r.a_.reserve(a_.size() - 1);
        for (std::size_t k = 1; k < a_.size(); ++k)
            r.a_.push_back(T(static_cast<int>(k)) * a_[k]);
        return r;
    }

private:
    poly1() : a_() {}
    std::vector<T> a_;
};

// ---------------------------------------------------------------------------
// compose_workspace<D,T>: the two ping-pong buffers of the Horner loop.
// Reused across element loops (same (K, n) => no allocation from the second
// call on; internal design 4.3).
// ---------------------------------------------------------------------------
template <int D, typename T>
struct compose_workspace {
    bpoly<D, T> buf_a, buf_b;
};

// ---------------------------------------------------------------------------
// F6: composition f(u) via Horner on the Bernstein algebra:
//   r = a_K; for k = K-1 .. 0: r = r * u + a_k
// Intermediate degrees increase monotonically (n, 2n, ..., K n); the scalar
// addition is exact thanks to the partition of unity. K == 0 returns the
// degree 0 constant a_0.
// ---------------------------------------------------------------------------
template <int D, typename T>
void compose_into(bpoly<D, T>& dst, const poly1<T>& f, const bpoly<D, T>& u,
                  compose_workspace<D, T>& ws) {
    assert(&dst != &u);
    const int K = f.degree();
    if (K == 0) {
        detail::bpoly_access::prepare(dst, 0, false);
        detail::bpoly_access::vec(dst)[0] = f.coeff(0);
        return;
    }
    bpoly<D, T>* r = &ws.buf_a;
    bpoly<D, T>* tmp = &ws.buf_b;
    detail::bpoly_access::prepare(*r, 0, false);
    detail::bpoly_access::vec(*r)[0] = f.coeff(K);
    for (int k = K - 1; k >= 0; --k) {
        mul_into(*tmp, *r, u);
        tmp->add_scalar(f.coeff(k));
        bpoly<D, T>* sw = r;
        r = tmp;
        tmp = sw;
    }
    // copy the final buffer into dst (assignment; buffer reuse preserved)
    detail::bpoly_access::prepare(dst, r->degree(), false);
    std::vector<T>& d = detail::bpoly_access::vec(dst);
    const std::vector<T>& s = detail::bpoly_access::vec(*r);
    for (std::size_t i = 0; i < d.size(); ++i) d[i] = s[i];
}

template <int D, typename T>
bpoly<D, T> compose(const poly1<T>& f, const bpoly<D, T>& u) {
    compose_workspace<D, T> ws;
    bpoly<D, T> r;
    compose_into(r, f, u, ws);
    return r;
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_POLY1_HPP
