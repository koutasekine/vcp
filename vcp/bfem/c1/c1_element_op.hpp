// vcp/bfem/c1/c1_element_op.hpp
// Phase 6 (2D C1 Argyris family): local element kernels -- the M_T sandwich
// around the rational-stage contracted reference tables.
//
// Conforms to: C1 external design v0.2 (sections 3, 4, 5) and
//              C1 internal design v0.2 (section 4).
//
// Kernel structure (internal design 4):
//   local matrix = P^T (metric-weighted contracted reference table) P
//     stiffness:        sum_{i<=j} w'_{ij} S^{(ij)} fold   (w' as in P^n L2)
//     laplacian:        |T| sum_{p<=q} c_p c_q R2^{(p)(q)} fold,
//                       c_p = g''_{ii} or 2 g''_{ij} (i < j),
//                       g''_{ij} = grad lambda_i . grad lambda_j
//     hessian (D2:D2):  same blocks, weights W_{pq} = sum over the pair
//                       expansions of g''_{ik} g''_{jl}
//     mass/mixed:       |T| x contracted mass
//   field gather:       hat_c = P (l2g gather), Bernstein coeffs = C hat_c
//   laplacian_field:    Delta coeffs = sum_p c_p (dlambda_i dlambda_j) --
//                       derivative_map twice, broken P_{k-2}
//
// Division count: ZERO in every kernel and in the P / P^{-1} applications
// (S-C1-2: the per-element division happened in element_geometry, the
// per-edge division in the c1_geometry edge cache).
//
// Thread model: one instance is NOT thread safe (owned buffers), the shared
// registries are safe -- the fe_space / element_op contract.

#ifndef VCP_BFEM_C1_C1_ELEMENT_OP_HPP
#define VCP_BFEM_C1_C1_ELEMENT_OP_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>

#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/c1/c1_tables.hpp>
#include <vcp/bfem/c1/c1_geometry.hpp>

namespace vcp {
namespace bfem {

template <typename T, typename P = vcp::mats<T> >
class c1_element_op {
public:
    c1_element_op()
        : geom_(detail::geometry_access::make_empty<2, T>()),
          geometry_set_(false) {}

    // element switch: metric weights (multiplications only) + pullback slots
    // are marked stale and lazily rebuilt per requested degree.
    void set_geometry(const element_geometry<2, T>& g,
                      const std::array<T, 3>& inv_tsq) {
        geom_ = g;
        inv_tsq_ = inv_tsq;
        // first-derivative metric (the P^n L2 formula: cof dot cof x
        // inv_absdet x 1/2)
        const T c = g.inv_absdet() * half();
        for (int i = 0; i <= 2; ++i)
            for (int j = i; j <= 2; ++j) {
                T dot = detail::geometry_access::cof(g, i, 0)
                            * detail::geometry_access::cof(g, j, 0)
                        + detail::geometry_access::cof(g, i, 1)
                              * detail::geometry_access::cof(g, j, 1);
                w1_[static_cast<std::size_t>(
                    detail::ref_stiff_block_index(2, i, j))] = dot * c;
            }
        // second-derivative metric g''_{ij} and the Delta weights c_p
        for (int i = 0; i <= 2; ++i)
            for (int j = i; j <= 2; ++j) {
                T dot = g.grad_lambda(i, 0) * g.grad_lambda(j, 0)
                        + g.grad_lambda(i, 1) * g.grad_lambda(j, 1);
                int p = detail::ref_stiff_block_index(2, i, j);
                gpp_[static_cast<std::size_t>(p)] = dot;
                cp_[static_cast<std::size_t>(p)] = (i == j) ? dot : dot + dot;
            }
        for (std::size_t s = 0; s < pulls_.size(); ++s)
            pulls_[s].second = false;                    // stale
        geometry_set_ = true;
    }

    const element_geometry<2, T>& geometry() const {
        require_geometry();
        return geom_;
    }

    // pullback of the current element at degree k (lazy slot reuse: no
    // allocation from the second element on for a fixed degree set)
    const detail::c1_pullback<T>& pull(int k) {
        require_geometry();
        detail::c1_check_k(k, "c1_element_op::pull");
        for (std::size_t s = 0; s < pulls_.size(); ++s) {
            if (pulls_[s].first.order() == k) {
                if (!pulls_[s].second) {
                    pulls_[s].first.build(k, geom_, inv_tsq_,
                                          c1_typed_registry<T>::hermite_w(k));
                    pulls_[s].second = true;
                }
                return pulls_[s].first;
            }
        }
        pulls_.push_back(std::make_pair(detail::c1_pullback<T>(), true));
        pulls_.back().first.build(k, geom_, inv_tsq_,
                                  c1_typed_registry<T>::hermite_w(k));
        return pulls_.back().first;
    }

    // ---- local stiffness (grad, grad): P^T [sum w' S] P ----
    void local_stiffness(int k, vcp::matrix<T, P>& out) {
        require_geometry();
        const typename c1_typed_registry<T>::block_table& S =
            c1_typed_registry<T>::stiffness(k);
        const int dim = detail::c1_dim(k);
        out.zeros(dim, dim);
        for (int i = 0; i <= 2; ++i) {
            for (int j = i; j <= 2; ++j) {
                const T& w = w1_[static_cast<std::size_t>(
                    detail::ref_stiff_block_index(2, i, j))];
                rt_block_view<T> B = S.block(i, j);
                if (i == j) {
                    for (int a = 0; a < dim; ++a)
                        for (int b = 0; b < dim; ++b)
                            out(a, b) += w * B.at(a, b);
                } else {
                    for (int a = 0; a < dim; ++a)
                        for (int b = 0; b < dim; ++b)
                            out(a, b) += w * (B.at(a, b) + B.at(b, a));
                }
            }
        }
        sandwich(pull(k), pull(k), out);
    }

    // ---- local (mixed-degree) mass: P_a^T [|T| CM] P_b ----
    void local_mass(int a, int b, vcp::matrix<T, P>& out) {
        require_geometry();
        const typename c1_typed_registry<T>::mat_table& M =
            c1_typed_registry<T>::mass_contracted(a, b);
        out.zeros(M.rows(), M.cols());
        const T& mea = geom_.measure();
        for (int i = 0; i < M.rows(); ++i)
            for (int j = 0; j < M.cols(); ++j)
                out(i, j) = mea * M.at(i, j);
        sandwich(pull(a), pull(b), out);
    }

    // ---- local laplacian stiffness (Delta, Delta) ----
    void local_laplacian(int k, vcp::matrix<T, P>& out) {
        require_geometry();
        T w[6];
        const T& mea = geom_.measure();
        for (int p = 0; p < 6; ++p)
            w[p] = mea * cp_[static_cast<std::size_t>(p)];
        r2_weighted(k, w, cp_.data(), out);
        sandwich(pull(k), pull(k), out);
    }

    // ---- local full-Hessian stiffness (D^2 u : D^2 v) ----
    void local_hessian(int k, vcp::matrix<T, P>& out) {
        require_geometry();
        // W_{pq} = sum over expansions of p, q of g''_{ik} g''_{jl}
        T wfull[6][6];
        for (int p = 0; p < 6; ++p) {
            int pi, pj;
            detail::c1_pair_dirs(p, pi, pj);
            for (int q = p; q < 6; ++q) {
                int qi, qj;
                detail::c1_pair_dirs(q, qi, qj);
                T acc = gpp(pi, qi) * gpp(pj, qj);
                if (qi != qj) acc += gpp(pi, qj) * gpp(pj, qi);
                if (pi != pj) {
                    acc += gpp(pj, qi) * gpp(pi, qj);
                    if (qi != qj) acc += gpp(pj, qj) * gpp(pi, qi);
                }
                wfull[p][q] = geom_.measure() * acc;
            }
        }
        r2_weighted_full(k, wfull, out);
        sandwich(pull(k), pull(k), out);
    }

    // ---- local load (w, Phi_a): P^T [|T| CM^{(k, deg w)} coeffs(w)] ----
    void local_load(const bpoly<2, T>& w, int k, vcp::matrix<T, P>& out) {
        require_geometry();
        const typename c1_typed_registry<T>::mat_table& CM =
            c1_typed_registry<T>::load_contracted(k, w.degree());
        const int dim = CM.rows();
        grow(fhat_, dim);
        const T& mea = geom_.measure();
        const std::vector<T>& wc = detail::bpoly_access::vec(w);
        for (int a = 0; a < dim; ++a) {
            T acc = CM.at(a, 0) * wc[0];
            for (int q = 1; q < CM.cols(); ++q)
                acc += CM.at(a, q) * wc[static_cast<std::size_t>(q)];
            fhat_[static_cast<std::size_t>(a)] = mea * acc;
        }
        grow(fout_, dim);
        pull(k).apply_transpose(fhat_.data(), fout_.data());
        out.zeros(dim, 1);
        for (int a = 0; a < dim; ++a)
            out(a, 0) = fout_[static_cast<std::size_t>(a)];
    }

    // ---- local weighted mass (w Phi_b, Phi_a) ----
    void local_weighted_mass(const bpoly<2, T>& w, int k, vcp::matrix<T, P>& out) {
        require_geometry();
        // Bernstein-level J via the moment-vector method (the L2 V4 recipe)
        const typed_mass_table<2, T>& M2 =
            typed_registry<2, T>::mass(2 * k, w.degree());
        grow(mom_, M2.rows());
        const std::vector<T>& wc = detail::bpoly_access::vec(w);
        const T& mea = geom_.measure();
        for (int g = 0; g < M2.rows(); ++g) {
            T acc = M2.at(g, 0) * wc[0];
            for (int j = 1; j < M2.cols(); ++j)
                acc += M2.at(g, j) * wc[static_cast<std::size_t>(j)];
            mom_[static_cast<std::size_t>(g)] = mea * acc;
        }
        const typed_product_table<2, T>& Pt = typed_registry<2, T>::product(k, k);
        const int N = Pt.rows();
        // J_B then the C-contraction: out = P^T C^T J_B C P
        jb_.zeros(N, N);
        for (int a = 0; a < N; ++a)
            for (int b = 0; b < N; ++b)
                jb_(a, b) = Pt.coeff(a, b)
                            * mom_[static_cast<std::size_t>(Pt.target_rank(a, b))];
        const typename c1_typed_registry<T>::mat_table& C =
            c1_typed_registry<T>::basis(k);
        const int dim = C.cols();
        // tmp = J_B C (N x dim), out = C^T tmp (dim x dim)
        tmp_.zeros(N, dim);
        for (int a = 0; a < N; ++a)
            for (int c = 0; c < dim; ++c) {
                T acc = jb_(a, 0) * C.at(0, c);
                for (int b = 1; b < N; ++b)
                    acc += jb_(a, b) * C.at(b, c);
                tmp_(a, c) = acc;
            }
        out.zeros(dim, dim);
        for (int r = 0; r < dim; ++r)
            for (int c = 0; c < dim; ++c) {
                T acc = C.at(0, r) * tmp_(0, c);
                for (int b = 1; b < N; ++b)
                    acc += C.at(b, r) * tmp_(b, c);
                out(r, c) = acc;
            }
        sandwich(pull(k), pull(k), out);
    }

    // ---- local laplacian mixed (Delta Phi_c, q_g): rows broken P_l ----
    void local_laplacian_mixed(int k, int l, vcp::matrix<T, P>& out) {
        require_geometry();
        const typename c1_typed_registry<T>::block_table& LM =
            c1_typed_registry<T>::lap_mixed(k, l);
        const int rows = LM.block_rows();
        const int dim = LM.block_cols();
        out.zeros(rows, dim);
        const T& mea = geom_.measure();
        for (int p = 0; p < 6; ++p) {
            int i, j;
            detail::c1_pair_dirs(p, i, j);
            rt_block_view<T> B = LM.block(i, j);
            const T w = mea * cp_[static_cast<std::size_t>(p)];
            for (int g = 0; g < rows; ++g)
                for (int c = 0; c < dim; ++c)
                    out(g, c) += w * B.at(g, c);
        }
        // columns: out := out P  (row transform by P^T)
        const detail::c1_pullback<T>& Pk = pull(k);
        grow(xbuf_, dim);
        grow(ybuf_, dim);
        for (int g = 0; g < rows; ++g) {
            for (int c = 0; c < dim; ++c)
                xbuf_[static_cast<std::size_t>(c)] = out(g, c);
            Pk.apply_transpose(xbuf_.data(), ybuf_.data());
            for (int c = 0; c < dim; ++c)
                out(g, c) = ybuf_[static_cast<std::size_t>(c)];
        }
    }

    // ---- local laplacian load (w, Delta Phi_a) ----
    void local_laplacian_load(const bpoly<2, T>& w, int k, vcp::matrix<T, P>& out) {
        require_geometry();
        const typed_mass_table<2, T>& M2 =
            typed_registry<2, T>::mass(k - 2, w.degree());
        grow(mom_, M2.rows());
        const std::vector<T>& wc = detail::bpoly_access::vec(w);
        for (int g = 0; g < M2.rows(); ++g) {
            T acc = M2.at(g, 0) * wc[0];
            for (int j = 1; j < M2.cols(); ++j)
                acc += M2.at(g, j) * wc[static_cast<std::size_t>(j)];
            mom_[static_cast<std::size_t>(g)] = acc;
        }
        const typename c1_typed_registry<T>::block_table& D2 =
            c1_typed_registry<T>::d2c(k);
        const int dim = D2.block_cols();
        const int rows = D2.block_rows();
        grow(fhat_, dim);
        const T& mea = geom_.measure();
        for (int a = 0; a < dim; ++a) fhat_[static_cast<std::size_t>(a)] = T(0);
        for (int p = 0; p < 6; ++p) {
            int i, j;
            detail::c1_pair_dirs(p, i, j);
            rt_block_view<T> B = D2.block(i, j);
            const T w6 = mea * cp_[static_cast<std::size_t>(p)];
            for (int a = 0; a < dim; ++a) {
                T acc = B.at(0, a) * mom_[0];
                for (int g = 1; g < rows; ++g)
                    acc += B.at(g, a) * mom_[static_cast<std::size_t>(g)];
                fhat_[static_cast<std::size_t>(a)] += w6 * acc;
            }
        }
        grow(fout_, dim);
        pull(k).apply_transpose(fhat_.data(), fout_.data());
        out.zeros(dim, 1);
        for (int a = 0; a < dim; ++a)
            out(a, 0) = fout_[static_cast<std::size_t>(a)];
    }

    // ---- field gather: local Bernstein coefficients of u o F ----
    // (hat_c = P x, then C hat_c)
    void gather_hat(int k, const T* x, bpoly<2, T>& out) {
        require_geometry();
        const detail::c1_pullback<T>& Pk = pull(k);
        const int dim = Pk.dim();
        grow(xbuf_, dim);
        Pk.apply(x, xbuf_.data());
        const typename c1_typed_registry<T>::mat_table& C =
            c1_typed_registry<T>::basis(k);
        detail::bpoly_access::prepare(out, k, false);
        std::vector<T>& oc = detail::bpoly_access::vec(out);
        for (int r = 0; r < C.rows(); ++r) {
            T acc = C.at(r, 0) * xbuf_[0];
            for (int a = 1; a < dim; ++a)
                acc += C.at(r, a) * xbuf_[static_cast<std::size_t>(a)];
            oc[static_cast<std::size_t>(r)] = acc;
        }
    }

    // ---- Delta of a local field: broken P_{k-2} coefficients ----
    // Delta u = sum_p c_p dlambda_i dlambda_j u (derivative_map twice)
    void laplacian_coeffs(const bpoly<2, T>& uhat, bpoly<2, T>& out) {
        require_geometry();
        const int k = uhat.degree();
        detail::c1_check_k(k, "c1_element_op::laplacian_coeffs");
        const derivative_map<2>& dmk = detail::deriv_cache<2>::get(k);
        const derivative_map<2>& dm1 = detail::deriv_cache<2>::get(k - 1);
        detail::bpoly_access::prepare(out, k - 2, true);
        std::vector<T>& oc = detail::bpoly_access::vec(out);
        const std::vector<T>& uc = detail::bpoly_access::vec(uhat);
        const T fac = T(k) * T(k - 1);
        for (int p = 0; p < 6; ++p) {
            int i, j;
            detail::c1_pair_dirs(p, i, j);
            const T w = fac * cp_[static_cast<std::size_t>(p)];
            for (int r = 0; r < dmk.source_size(); ++r) {
                int t1 = dmk.target(r, i);
                if (t1 < 0) continue;
                int t2 = dm1.target(t1, j);
                if (t2 < 0) continue;
                oc[static_cast<std::size_t>(t2)] += w * uc[static_cast<std::size_t>(r)];
            }
        }
    }

    // ---- physical gradient component of a local field (L2 G8 formula) ----
    void grad_component(const bpoly<2, T>& u, int d, bpoly<2, T>& out) {
        require_geometry();
        assert(&out != &u);
        if (d < 0 || d >= 2)
            throw std::invalid_argument("bfem::c1_element_op::grad_component: bad d");
        dlambda_into(scr_[0], u, 0);
        scale_into(out, scr_[0], geom_.grad_lambda(0, d));
        for (int i = 1; i <= 2; ++i) {
            dlambda_into(scr_[0], u, i);
            scale_into(scr_[1], scr_[0], geom_.grad_lambda(i, d));
            out += scr_[1];
        }
    }

    // ---- local scalar ----
    T local_inner(const bpoly<2, T>& a, const bpoly<2, T>& b) const {
        require_geometry();
        return geom_.measure() * ::vcp::bfem::inner(a, b);
    }

private:
    element_geometry<2, T> geom_;
    std::array<T, 3> inv_tsq_;
    std::array<T, 6> w1_;              // first-derivative metric weights
    std::array<T, 6> gpp_;             // g''_{ij}, packed i <= j
    std::array<T, 6> cp_;              // Delta weights per direction pair
    std::vector<std::pair<detail::c1_pullback<T>, bool> > pulls_;
    std::vector<T> xbuf_, ybuf_, fhat_, fout_, mom_;
    vcp::matrix<T, P> jb_, tmp_;
    bpoly<2, T> scr_[2];
    bool geometry_set_;

    void require_geometry() const {
        if (!geometry_set_)
            throw std::logic_error("bfem::c1_element_op: set_geometry not called");
    }
    static const T& half() {
        static const T c = rational_to<T>(1, 2);
        return c;
    }
    static void grow(std::vector<T>& v, int n) {
        if (v.size() < static_cast<std::size_t>(n))
            v.resize(static_cast<std::size_t>(n));
    }
    const T& gpp(int i, int j) const {
        return gpp_[static_cast<std::size_t>(
            detail::ref_stiff_block_index(2, i < j ? i : j, i < j ? j : i))];
    }

    // out += sum_{p<=q} weights (diag: wd[p] applied as full c_p c_q fold)
    // laplacian form: out = sum_p |T| c_p c_p R2^{pp}
    //                     + sum_{p<q} |T| c_p c_q (R2^{pq} + R2^{pq T})
    void r2_weighted(int k, const T* mcp, const T* cp, vcp::matrix<T, P>& out) {
        const typename c1_typed_registry<T>::pair_table& R2 =
            c1_typed_registry<T>::r2(k);
        const int dim = R2.rows();
        out.zeros(dim, dim);
        for (int p = 0; p < 6; ++p) {
            for (int q = p; q < 6; ++q) {
                rt_block_view<T> B = R2.block(p, q);
                const T w = mcp[p] * cp[q];
                if (p == q) {
                    for (int a = 0; a < dim; ++a)
                        for (int b = 0; b < dim; ++b)
                            out(a, b) += w * B.at(a, b);
                } else {
                    for (int a = 0; a < dim; ++a)
                        for (int b = 0; b < dim; ++b)
                            out(a, b) += w * (B.at(a, b) + B.at(b, a));
                }
            }
        }
    }

    // hessian form with a full symmetric weight table (upper triangle given)
    void r2_weighted_full(int k, const T wfull[6][6], vcp::matrix<T, P>& out) {
        const typename c1_typed_registry<T>::pair_table& R2 =
            c1_typed_registry<T>::r2(k);
        const int dim = R2.rows();
        out.zeros(dim, dim);
        for (int p = 0; p < 6; ++p) {
            for (int q = p; q < 6; ++q) {
                rt_block_view<T> B = R2.block(p, q);
                const T& w = wfull[p][q];
                if (p == q) {
                    for (int a = 0; a < dim; ++a)
                        for (int b = 0; b < dim; ++b)
                            out(a, b) += w * B.at(a, b);
                } else {
                    for (int a = 0; a < dim; ++a)
                        for (int b = 0; b < dim; ++b)
                            out(a, b) += w * (B.at(a, b) + B.at(b, a));
                }
            }
        }
    }

    // A := Pr^T A Pc (structured application; scratch vectors only)
    void sandwich(const detail::c1_pullback<T>& Pr,
                  const detail::c1_pullback<T>& Pc, vcp::matrix<T, P>& A) {
        const int rows = Pr.dim();
        const int cols = Pc.dim();
        grow(xbuf_, rows > cols ? rows : cols);
        grow(ybuf_, rows > cols ? rows : cols);
        for (int a = 0; a < rows; ++a) {                 // rows: A := A Pc
            for (int c = 0; c < cols; ++c)
                xbuf_[static_cast<std::size_t>(c)] = A(a, c);
            Pc.apply_transpose(xbuf_.data(), ybuf_.data());
            for (int c = 0; c < cols; ++c)
                A(a, c) = ybuf_[static_cast<std::size_t>(c)];
        }
        for (int c = 0; c < cols; ++c) {                 // cols: A := Pr^T A
            for (int a = 0; a < rows; ++a)
                xbuf_[static_cast<std::size_t>(a)] = A(a, c);
            Pr.apply_transpose(xbuf_.data(), ybuf_.data());
            for (int a = 0; a < rows; ++a)
                A(a, c) = ybuf_[static_cast<std::size_t>(a)];
        }
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_C1_C1_ELEMENT_OP_HPP
