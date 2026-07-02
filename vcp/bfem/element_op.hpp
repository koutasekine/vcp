// vcp/bfem/element_op.hpp
// Layer 2: element operations element_op<D,T,P> (G2-G8).
//
// Conforms to: L2 external design v0.3 (sections 4, 5) and
//              L2 internal design v0.3 (sections 3, 5, 6, 7).
//
// matrix.hpp reconciliation (external design 11.2, resolved 2026-07-02, see
// sandbox/docs/reviews/bfem_api_reconciliation_review.md): vcp::matrix<T,P>
// with the mats policy satisfies "size-preserving re-setting does not
// reallocate" -- zeros(r, c) calls std::vector::resize (no-op for an equal
// size) and zero-fills BY ASSIGNMENT. The main plan (direct matrix output) is
// therefore adopted; the local_dense fallback is not needed. The kernels use
// only zeros(r, c) and operator()(i, j) so a future adapter swap stays a
// one-file change (internal design 10.2).
//
// Thread model (external 4.1): one element_op instance is NOT thread safe;
// clone the instance per thread (the shared L0/L2 caches are safe).

#ifndef VCP_BFEM_ELEMENT_OP_HPP
#define VCP_BFEM_ELEMENT_OP_HPP

#include <array>
#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>

#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/ref_stiffness.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/poly1.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P = vcp::mats<T> >
class element_op {
    static_assert(D == 2 || D == 3, "bfem::element_op: only D == 2 or D == 3");
public:
    element_op()
        : geom_(detail::geometry_access::make_empty<D, T>()),
          geometry_set_(false) {}

    // element switch; the single division of the element happened inside
    // element_geometry::from_vertices (V6). Here: copy + multiplications.
    void set_geometry(const element_geometry<D, T>& g) {
        geom_ = g;
        const T c = g.inv_absdet() * dfact_inv();
        for (int i = 0; i <= D; ++i) {
            for (int j = i; j <= D; ++j) {
                T dot = detail::geometry_access::cof(g, i, 0)
                        * detail::geometry_access::cof(g, j, 0);
                for (int d = 1; d < D; ++d)
                    dot += detail::geometry_access::cof(g, i, d)
                           * detail::geometry_access::cof(g, j, d);
                gprime_[static_cast<std::size_t>(
                    detail::ref_stiff_block_index(D, i, j))] = dot * c;
            }
        }
        geometry_set_ = true;
    }

    // ---- G2: local stiffness (geometry separated, V3). n >= 1 (C-3) ----
    void local_stiffness(int n, vcp::matrix<T, P>& out) {
        require_geometry();
        if (n < 1)
            throw std::invalid_argument("bfem::element_op::local_stiffness: n < 1");
        const typename detail::ref_stiffness_cache<D, T>::tensor& R =
            detail::ref_stiffness_cache<D, T>::get(n);
        const int N = R.N;
        out.zeros(N, N);
        for (int i = 0; i <= D; ++i) {
            for (int j = i; j <= D; ++j) {
                const T& w = gprime_[static_cast<std::size_t>(
                    detail::ref_stiff_block_index(D, i, j))];
                if (i == j) {
                    for (int a = 0; a < N; ++a)
                        for (int b = 0; b < N; ++b)
                            out(a, b) += w * R.at(i, i, a, b);
                } else {
                    // R^{(ji)}[a][b] = R^{(ij)}[b][a]: fold the symmetric pair
                    for (int a = 0; a < N; ++a)
                        for (int b = 0; b < N; ++b)
                            out(a, b) += w * (R.at(i, j, a, b) + R.at(i, j, b, a));
                }
            }
        }
    }

    // ---- G3: local (mixed-degree) mass ----
    void local_mass(int a, int b, vcp::matrix<T, P>& out) {
        require_geometry();
        if (a < 0 || b < 0)
            throw std::invalid_argument("bfem::element_op::local_mass: negative degree");
        const typed_mass_table<D, T>& M = typed_registry<D, T>::mass(a, b);
        out.zeros(M.rows(), M.cols());
        const T& mea = geom_.measure();
        for (int i = 0; i < M.rows(); ++i)
            for (int j = 0; j < M.cols(); ++j)
                out(i, j) = mea * M.at(i, j);
    }

    // ---- G4: local load F_alpha = (w, phi^n_alpha)_T (N x 1) ----
    void local_load(const bpoly<D, T>& w, int n, vcp::matrix<T, P>& out) {
        require_geometry();
        if (n < 0)
            throw std::invalid_argument("bfem::element_op::local_load: n < 0");
        moment_vec(0, w, n);
        const int N = coeff_registry<D>::indices(n).size();
        out.zeros(N, 1);
        const T& mea = geom_.measure();
        for (int a = 0; a < N; ++a)
            out(a, 0) = mea * mom_[0][static_cast<std::size_t>(a)];
    }
    // convenience overload (external 5.3): w = f(u_h) composed internally
    void local_load(const poly1<T>& f, const bpoly<D, T>& uh, int n,
                    vcp::matrix<T, P>& out) {
        compose_into(wbuf_, f, uh, cws_);
        local_load(wbuf_, n, out);
    }

    // ---- G5: local weighted mass (moment-vector method, V4) ----
    // J_{ab} = |T| c(a,b) m[target(a,b)],  m = M^{2n,c} coeffs(w)
    void local_weighted_mass(const bpoly<D, T>& w, int n, vcp::matrix<T, P>& out) {
        require_geometry();
        if (n < 0)
            throw std::invalid_argument("bfem::element_op::local_weighted_mass: n < 0");
        moment_vec(0, w, 2 * n);
        const int Nq = coeff_registry<D>::indices(2 * n).size();
        const T& mea = geom_.measure();
        for (int g = 0; g < Nq; ++g)
            mom_[0][static_cast<std::size_t>(g)] *= mea;   // fold |T| into m
        const typed_product_table<D, T>& Pt = typed_registry<D, T>::product(n, n);
        const int N = Pt.rows();
        out.zeros(N, N);
        for (int a = 0; a < N; ++a)
            for (int b = 0; b < N; ++b)
                out(a, b) = Pt.coeff(a, b)
                            * mom_[0][static_cast<std::size_t>(Pt.target_rank(a, b))];
    }
    void local_weighted_mass(const poly1<T>& fprime, const bpoly<D, T>& uh, int n,
                             vcp::matrix<T, P>& out) {
        compose_into(wbuf_, fprime, uh, cws_);
        local_weighted_mass(wbuf_, n, out);
    }

    // ---- G6: local convection C_{ab} = (b.grad phi_b, phi_a)_T. n >= 1 (C-3) ----
    void local_convection(const std::array<bpoly<D, T>, D>& b, int n,
                          vcp::matrix<T, P>& out) {
        require_geometry();
        if (n < 1)
            throw std::invalid_argument("bfem::element_op::local_convection: n < 1");
        // w_k = b . grad_lambda_k (bpoly), then m_k = M^{2n-1, deg w_k} coeffs(w_k)
        for (int k = 0; k <= D; ++k) {
            scale_into(scratch_[0], b[0], geom_.grad_lambda(k, 0));
            for (int d = 1; d < D; ++d) {
                scale_into(scratch_[1], b[static_cast<std::size_t>(d)],
                           geom_.grad_lambda(k, d));
                add_into(scratch_[2], scratch_[0], scratch_[1]);
                std::swap(scratch_[0], scratch_[2]);
            }
            moment_vec(k, scratch_[0], 2 * n - 1);
        }
        const typed_product_table<D, T>& Pt = typed_registry<D, T>::product(n - 1, n);
        const derivative_map<D>& dm = detail::deriv_cache<D>::get(n);
        const int N = coeff_registry<D>::indices(n).size();
        out.zeros(N, N);
        const T mn = geom_.measure() * T(n);
        for (int a = 0; a < N; ++a) {
            for (int bi = 0; bi < N; ++bi) {
                T acc(0);
                for (int k = 0; k <= D; ++k) {
                    int r = dm.target(bi, k);
                    if (r < 0) continue;
                    acc += Pt.coeff(r, a)
                           * mom_[static_cast<std::size_t>(k)]
                                 [static_cast<std::size_t>(Pt.target_rank(r, a))];
                }
                out(a, bi) = mn * acc;
            }
        }
    }

    // ---- G7: local scalar (physical value of the L1 inner product) ----
    T local_inner(const bpoly<D, T>& w1, const bpoly<D, T>& w2) const {
        require_geometry();
        return geom_.measure() * ::vcp::bfem::inner(w1, w2);
    }

    // ---- G8: gradient assembly ----
    // (grad u)_d = sum_i dlambda(u, i) * grad_lambda(i, d), degree n-1
    void grad_component(const bpoly<D, T>& u, int d, bpoly<D, T>& out) {
        require_geometry();
        assert(&out != &u);
        if (d < 0 || d >= D)
            throw std::invalid_argument("bfem::element_op::grad_component: bad d");
        dlambda_into(scratch_[3], u, 0);
        scale_into(out, scratch_[3], geom_.grad_lambda(0, d));
        for (int i = 1; i <= D; ++i) {
            dlambda_into(scratch_[3], u, i);
            scale_into(scratch_[4], scratch_[3], geom_.grad_lambda(i, d));
            out += scratch_[4];                       // same degree n-1
        }
    }
    // b . grad u  (degree deg b + n - 1; components of b may differ in degree)
    void b_dot_grad(const std::array<bpoly<D, T>, D>& b, const bpoly<D, T>& u,
                    bpoly<D, T>& out) {
        require_geometry();
        assert(&out != &u);
        grad_component(u, 0, scratch_[5]);
        mul_into(scratch_[6], b[0], scratch_[5]);
        for (int d = 1; d < D; ++d) {
            grad_component(u, d, scratch_[5]);
            mul_into(scratch_[7], b[static_cast<std::size_t>(d)], scratch_[5]);
            add_into(scratch_[2], scratch_[6], scratch_[7]);
            std::swap(scratch_[6], scratch_[2]);
        }
        detail::bpoly_access::prepare(out, scratch_[6].degree(), false);
        std::vector<T>& o = detail::bpoly_access::vec(out);
        const std::vector<T>& s = detail::bpoly_access::vec(scratch_[6]);
        for (std::size_t i = 0; i < o.size(); ++i) o[i] = s[i];
    }

private:
    element_geometry<D, T> geom_;                            // copy (V1)
    std::array<T, static_cast<std::size_t>((D + 1) * (D + 2) / 2)> gprime_;
    std::array<std::vector<T>, static_cast<std::size_t>(D + 1)> mom_;
    std::array<bpoly<D, T>, 8> scratch_;
    compose_workspace<D, T> cws_;
    bpoly<D, T> wbuf_;
    bool geometry_set_;

    void require_geometry() const {
        if (!geometry_set_)
            throw std::logic_error("bfem::element_op: set_geometry not called");
    }

    static const T& dfact_inv() {
        static const T c = rational_to<T>(1, detail::factorial_of<D>::value);
        return c;
    }

    // m = M^{q, deg w} coeffs(w) into mom_[slot] (grow-only buffer)
    void moment_vec(int slot, const bpoly<D, T>& w, int q) {
        const typed_mass_table<D, T>& M = typed_registry<D, T>::mass(q, w.degree());
        std::vector<T>& m = mom_[static_cast<std::size_t>(slot)];
        if (m.size() < static_cast<std::size_t>(M.rows()))
            m.resize(static_cast<std::size_t>(M.rows()));
        const std::vector<T>& wc = detail::bpoly_access::vec(w);
        for (int g = 0; g < M.rows(); ++g) {
            T& acc = m[static_cast<std::size_t>(g)];
            acc = M.at(g, 0) * wc[0];
            for (int j = 1; j < M.cols(); ++j)
                acc += M.at(g, j) * wc[static_cast<std::size_t>(j)];
        }
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_ELEMENT_OP_HPP
