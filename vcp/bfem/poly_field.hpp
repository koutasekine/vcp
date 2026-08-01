// vcp/bfem/poly_field.hpp
// PF-1 Layer 1: poly_field<D,T> -- a coordinate-dependent multivariate
// polynomial f(x) = sum_alpha c_alpha x^alpha (alpha a D-dimensional
// exponent), together with its exact restriction to a mesh element as a
// Bernstein polynomial (bpoly<D,T>).
//
// Conforms to: PF-1 design v1.0 (sections 2, 3) and
//             PF-1 implementation directive v1.0 (section 2).
//
// Construction of restrict_to (design 3.2): the coordinate function x_d is
// linear on every element, so its degree-1 Bernstein coefficients are the
// vertex coordinates themselves (linear precision; the canonical rank i of
// index_map<D>(1) is alpha = e_i, i.e. vertex i -- multi_index.hpp header
// note). Monomials x^alpha are exact degree products (mul_into) of the
// coordinate bpolys, raised to the common degree n (elevate_into), scaled by
// c_alpha and accumulated. The operation sequence is deterministic (term
// insertion order); for interval T the enclosure property of the bpoly
// product is inherited unchanged, and with dyadic-rational vertices and
// coefficients in an exact scalar the result is exact.
//
// This header must NOT include fe_space.hpp (the space headers include this
// one); its dependencies are fixed by the implementation directive section 1.

#ifndef VCP_BFEM_POLY_FIELD_HPP
#define VCP_BFEM_POLY_FIELD_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>

#include <vcp/error.hpp>
#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/bpoly.hpp>

namespace vcp {
namespace bfem {

namespace detail {

// minimal local enable_if so that only the D-matching add_term overload
// exists (directive 2: "SFINAE or static_assert, never both signatures");
// self-contained to keep this header's dependency list closed
template <bool B, typename U = void> struct pf_enable_if {};
template <typename U> struct pf_enable_if<true, U> { typedef U type; };

} // namespace detail

template <int D, typename T>
class poly_field {
    static_assert(D == 2 || D == 3, "bfem::poly_field: only D == 2 or D == 3");
public:
    poly_field() : terms_() {}                       // zero polynomial

    // add the monomial c * x0^i0 * x1^i1 (D == 2). Re-adding an existing
    // exponent pair merges by coefficient addition.
    template <int DD = D>
    typename detail::pf_enable_if<DD == 2, void>::type
    add_term(const T& c, int i0, int i1) {
        std::array<int, D> ex;
        ex[0] = i0;
        ex[1] = i1;
        push_term(c, ex);
    }

    // add the monomial c * x0^i0 * x1^i1 * x2^i2 (D == 3)
    template <int DD = D>
    typename detail::pf_enable_if<DD == 3, void>::type
    add_term(const T& c, int i0, int i1, int i2) {
        std::array<int, D> ex;
        ex[0] = i0;
        ex[1] = i1;
        ex[2] = i2;
        push_term(c, ex);
    }

    // total degree max_alpha |alpha| over the stored terms (0 for the zero
    // polynomial)
    int total_degree() const {
        int n = 0;
        for (std::size_t k = 0; k < terms_.size(); ++k) {
            int s = 0;
            for (int d = 0; d < D; ++d)
                s += terms_[k].first[static_cast<std::size_t>(d)];
            if (s > n) n = s;
        }
        return n;
    }

    // no terms, or every stored coefficient compares == to T(0) (for interval
    // T this is the "certainly zero point interval" comparison)
    bool is_zero() const {
        for (std::size_t k = 0; k < terms_.size(); ++k)
            if (!(terms_[k].second == T(0))) return false;
        return true;
    }

    // Bernstein coefficients (degree n) of f on element e of msh.
    // Requires n >= total_degree().
    bpoly<D, T> restrict_to(const mesh<D, T>& msh, int e, int n) const {
        if (e < 0 || e >= msh.num_elements())
            vcp::throw_error<vcp::invalid_argument>(
                "bfem::poly_field::restrict_to: element index out of range");
        if (n < total_degree())
            vcp::throw_error<vcp::invalid_argument>(
                "bfem::poly_field::restrict_to: n < total_degree()");

        // degree-1 coordinate bpolys: coefficient i = vertex i coordinate
        const std::array<int, D + 1>& el = msh.element(e);
        std::array<bpoly<D, T>, D> coord;
        for (int d = 0; d < D; ++d) {
            std::vector<T> c(static_cast<std::size_t>(D + 1));
            for (int k = 0; k <= D; ++k)
                c[static_cast<std::size_t>(k)] =
                    msh.vertex(el[static_cast<std::size_t>(k)])
                       [static_cast<std::size_t>(d)];
            coord[static_cast<std::size_t>(d)] =
                bpoly<D, T>::from_coeffs(1, std::move(c));
        }

        bpoly<D, T> acc = bpoly<D, T>::zero(n);
        bpoly<D, T> buf0, buf1, raised;
        for (std::size_t k = 0; k < terms_.size(); ++k) {
            const std::array<int, D>& ex = terms_[k].first;
            // monomial x^alpha as the exact degree product of the coordinate
            // bpolys (existing mul_into; dst never aliases an operand)
            bpoly<D, T>* cur = &buf0;
            bpoly<D, T>* nxt = &buf1;
            bool first = true;
            for (int d = 0; d < D; ++d) {
                const bpoly<D, T>& xd = coord[static_cast<std::size_t>(d)];
                for (int j = 0; j < ex[static_cast<std::size_t>(d)]; ++j) {
                    if (first) {
                        *cur = xd;
                        first = false;
                    } else {
                        mul_into(*nxt, *cur, xd);
                        bpoly<D, T>* t = cur;
                        cur = nxt;
                        nxt = t;
                    }
                }
            }
            if (first) *cur = bpoly<D, T>::constant(T(1));   // alpha = 0
            elevate_into(raised, *cur, n);                   // raise to n
            raised *= terms_[k].second;                      // scale by c_alpha
            acc += raised;                                   // same degree n
        }
        return acc;
    }

private:
    // (exponent, coefficient) list in insertion order; merging is a linear
    // search (directive 2.2: at most a few dozen terms, no map -- determinism
    // and dependency reduction)
    std::vector<std::pair<std::array<int, D>, T> > terms_;

    void push_term(const T& c, const std::array<int, D>& ex) {
        for (int d = 0; d < D; ++d)
            if (ex[static_cast<std::size_t>(d)] < 0)
                vcp::throw_error<vcp::invalid_argument>(
                    "bfem::poly_field::add_term: negative exponent");
        for (std::size_t k = 0; k < terms_.size(); ++k) {
            if (terms_[k].first == ex) {
                terms_[k].second += c;
                return;
            }
        }
        terms_.push_back(std::pair<std::array<int, D>, T>(ex, c));
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_POLY_FIELD_HPP
