// vcp/bfem/c1/c1_geometry.hpp
// Phase 6 (2D C1 Argyris family): the edge cache 1/|t_e|^2 and the
// per-element pullback map P_T (the basis-change of the C1 element).
//
// Conforms to: C1 external design v0.2 (section 3) and
//              C1 internal design v0.2 (section 3).
//
// M_T direction note (deviation record, see
// sandbox/docs/issues/bfem_c1_mt_direction_issue.md): the design section 3
// defines M_T by "the reference DOFs of the pullback are linear combinations
// of the physical DOFs". This file implements exactly that map and names it
// P (c1_pullback): (P c)_a = ell_hat_a(u o F) for c_j = ell_j(u). Under this
// definition the field gather is hat_c = P c and the local matrices are
// A = P^T (contracted reference table) P; the inverse P^{-1} (also provided,
// closed form) drives interpolation and elevation. The design text writes
// the local-matrix rule with the inverse convention; the two clauses of the
// design contradict each other and the tests (C-T1..C-T4) pin the math
// implemented here.
//
// Division contract (S-C1-2, the revised budget of external design 3):
//  - c1_edge_inv_tsq performs the ONE T division per mesh edge (the single
//    division operator of vcp/bfem/c1/ -- G-C1-1), after certifying
//    0 < |t_e|^2; an uncertified (degenerate or interval zero-straddling)
//    edge throws degenerate_element BEFORE the division (B-4 fail-fast);
//  - c1_pullback::build/apply/apply_transpose/apply_inverse contain NO
//    division: the inverse uses the closed forms adj(B) inv_det (gradient),
//    Sym^2(B^{-1}) (Hessian), and 1/alpha_e = |t_e|^2 inv_det (1/|e_hat|^2)
//    with the rational reference constant 1/|e_hat|^2 (enclose-once).
//
// The ND-row coupling uses the T-C3 Hermite table W: the tangential
// parameter derivative of the edge trace at the ND points is the W-weighted
// combination of the Hermite data (endpoint values / first / second
// parameter derivatives + interior trace values), whose geometric
// coefficients are the edge-vector components (coordinate differences).

#ifndef VCP_BFEM_C1_C1_GEOMETRY_HPP
#define VCP_BFEM_C1_C1_GEOMETRY_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/dofmap.hpp>            // detail::mesh_topology2
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/c1/c1_tables.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// edge cache: inv_tsq[edge id] = 1 / |t_e|^2 of the canonical physical edge
// vector. Exactly ONE division per edge; fail-fast before dividing.
// ---------------------------------------------------------------------------
template <typename T>
std::vector<T> c1_edge_inv_tsq(const mesh<2, T>& msh, const mesh_topology2& tp) {
    std::vector<T> out;
    out.reserve(tp.edges.size());
    for (std::size_t ed = 0; ed < tp.edges.size(); ++ed) {
        const std::array<T, 2>& a = msh.vertex(tp.edges[ed][0]);
        const std::array<T, 2>& b = msh.vertex(tp.edges[ed][1]);
        T t0 = b[0] - a[0];
        T t1 = b[1] - a[1];
        T tsq = t0 * t0 + t1 * t1;
        // certification BEFORE the division (B-4): degenerate edges and
        // interval inputs straddling zero both fail here, counter untouched
        if (!(T(0) < tsq))
            throw degenerate_element(
                "bfem::c1: |t_e|^2 not certified positive "
                "(degenerate edge or interval straddling zero)");
        out.push_back(T(1) / tsq);        // the single division of vcp/bfem/c1/
    }
    return out;
}

// ---------------------------------------------------------------------------
// c1_pullback<T>: P_T and its closed-form inverse as structured applications
// (never a dense matrix). Fields are open for the AXN-C1 in-test replication.
// ---------------------------------------------------------------------------
template <typename T>
class c1_pullback {
public:
    c1_pullback() : k_(0), dim_(0), W_(0) {}

    int order() const { return k_; }
    int dim() const { return dim_; }

    // build from the element geometry, the cached inv |t_e|^2 of the three
    // LOCAL edges (local edge e = global edge tri_edge[e]) and the typed
    // Hermite table W of degree k. No division, no allocation after the
    // first build at a given k (buffers are grow-only).
    void build(int k, const element_geometry<2, T>& g,
               const std::array<T, 3>& inv_tsq,
               const rt_mat_tbl<T>& W) {
        c1_check_k(k, "c1_pullback::build");
        assert(W.rows() == k - 4 && W.cols() == k + 1);
        k_ = k;
        dim_ = c1_dim(k);
        W_ = &W;
        // B and Sym^2(B^T) rows (physical Hessian dofs -> reference Hessian)
        for (int r = 0; r < 2; ++r)
            for (int c = 0; c < 2; ++c)
                b_[r][c] = g.edge_matrix(r, c);
        sym2_rows(b_, sym2_);
        // closed-form inverse materials: B^{-1} = adj(B) inv_det
        T inv_det = (g.orientation() > 0) ? g.inv_absdet() : -g.inv_absdet();
        binv_[0][0] = b_[1][1] * inv_det;
        binv_[0][1] = -b_[0][1] * inv_det;
        binv_[1][0] = -b_[1][0] * inv_det;
        binv_[1][1] = b_[0][0] * inv_det;
        sym2_rows(binv_, sym2inv_);
        // per local edge: t_e, (t0^2, 2 t0 t1, t1^2), alpha, beta and the
        // inverse coefficients (no division: |e_hat|^2 constants are exact)
        const std::array<std::array<T, 2>, 3>& v = g.vertices();
        for (int e = 0; e < 3; ++e) {
            const int s = (e + 1) % 3, q = (e + 2) % 3;
            T t0 = v[static_cast<std::size_t>(q)][0] - v[static_cast<std::size_t>(s)][0];
            T t1 = v[static_cast<std::size_t>(q)][1] - v[static_cast<std::size_t>(s)][1];
            tvec_[static_cast<std::size_t>(e)][0] = t0;
            tvec_[static_cast<std::size_t>(e)][1] = t1;
            T t01 = t0 * t1;
            ttvec_[static_cast<std::size_t>(e)][0] = t0 * t0;
            ttvec_[static_cast<std::size_t>(e)][1] = t01 + t01;
            ttvec_[static_cast<std::size_t>(e)][2] = t1 * t1;
            // B nu_hat_e: integer combination of the columns of B
            T bn0(0), bn1(0);
            for (int c = 0; c < 2; ++c) {
                int nu = rt_edge_normal(e, c);
                if (nu > 0) { bn0 += b_[0][c]; bn1 += b_[1][c]; }
                else if (nu < 0) { bn0 -= b_[0][c]; bn1 -= b_[1][c]; }
            }
            // nu_e = R_{-90} t_e = (t1, -t0)
            alpha_[static_cast<std::size_t>(e)] =
                (bn0 * t1 - bn1 * t0) * inv_tsq[static_cast<std::size_t>(e)];
            T bnt = bn0 * t0 + bn1 * t1;
            beta_[static_cast<std::size_t>(e)] =
                bnt * inv_tsq[static_cast<std::size_t>(e)];
            // 1/alpha = |t|^2 inv_det (1/|e_hat|^2); beta/alpha = (B nu . t)
            // inv_det (1/|e_hat|^2) -- multiplications only
            T tsq = ttvec_[static_cast<std::size_t>(e)][0]
                    + ttvec_[static_cast<std::size_t>(e)][2];
            if (c1_edge_sq(e) == 2) {
                inv_alpha_[static_cast<std::size_t>(e)] =
                    tsq * inv_det * ehalf();
                beta_over_alpha_[static_cast<std::size_t>(e)] =
                    bnt * inv_det * ehalf();
            } else {
                inv_alpha_[static_cast<std::size_t>(e)] = tsq * inv_det;
                beta_over_alpha_[static_cast<std::size_t>(e)] = bnt * inv_det;
            }
        }
        if (dbuf_.size() < static_cast<std::size_t>(k + 1))
            dbuf_.resize(static_cast<std::size_t>(k + 1));
    }

    // y = P x (reference DOFs of the pullback from local physical DOFs)
    void apply(const T* x, T* y) const {
        require_built();
        const int k = k_;
        for (int i = 0; i < 3; ++i) vertex_forward(i, x, y);
        for (int r = 18 + 3 * c1_nedge(k); r < dim_; ++r) y[r] = x[r];  // interior
        for (int e = 0; e < 3; ++e) {
            const int s = (e + 1) % 3, q = (e + 2) % 3;
            for (int i = 0; i < c1_ntrace(k); ++i) {
                int idx = c1_edge_trace_dof(k, e, i);
                y[idx] = x[idx];
            }
            hermite_data(e, s, q, x);
            const T& be = beta_[static_cast<std::size_t>(e)];
            for (int j = 0; j < c1_nnd(k); ++j) {
                T acc = W_->at(j, 0) * dbuf_[0];
                for (int r = 1; r <= k; ++r)
                    acc += W_->at(j, r) * dbuf_[static_cast<std::size_t>(r)];
                int idx = c1_edge_nd_dof(k, e, j);
                y[idx] = alpha_[static_cast<std::size_t>(e)] * x[idx] + be * acc;
            }
        }
    }

    // z = P^T y
    void apply_transpose(const T* y, T* z) const {
        require_built();
        const int k = k_;
        for (int i = 0; i < 3; ++i) {
            const int vb = c1_vertex_dof(i, 0);
            z[vb] = y[vb];
            // gradient block transpose: z_g = B y_g
            z[vb + 1] = b_[0][0] * y[vb + 1] + b_[0][1] * y[vb + 2];
            z[vb + 2] = b_[1][0] * y[vb + 1] + b_[1][1] * y[vb + 2];
            // Hessian block transpose: z_h = Sym2^T y_h
            for (int c = 0; c < 3; ++c)
                z[vb + 3 + c] = sym2_[0][c] * y[vb + 3]
                                + sym2_[1][c] * y[vb + 4]
                                + sym2_[2][c] * y[vb + 5];
        }
        for (int r = 18 + 3 * c1_nedge(k); r < dim_; ++r) z[r] = y[r];
        for (int e = 0; e < 3; ++e)
            for (int i = 0; i < c1_ntrace(k); ++i) {
                int idx = c1_edge_trace_dof(k, e, i);
                z[idx] = y[idx];
            }
        for (int e = 0; e < 3; ++e) {
            const int s = (e + 1) % 3, q = (e + 2) % 3;
            const T& be = beta_[static_cast<std::size_t>(e)];
            // s_r = beta sum_j W[j][r] y_nd(e, j); diagonal alpha to z_nd
            for (int r = 0; r <= k; ++r) {
                T acc = W_->at(0, r) * y[c1_edge_nd_dof(k, e, 0)];
                for (int j = 1; j < c1_nnd(k); ++j)
                    acc += W_->at(j, r) * y[c1_edge_nd_dof(k, e, j)];
                dbuf_[static_cast<std::size_t>(r)] = be * acc;
            }
            for (int j = 0; j < c1_nnd(k); ++j) {
                int idx = c1_edge_nd_dof(k, e, j);
                z[idx] = alpha_[static_cast<std::size_t>(e)] * y[idx];
            }
            scatter_data(e, s, q, z);
        }
    }

    // x = P^{-1} y (closed form; no division)
    void apply_inverse(const T* y, T* x) const {
        require_built();
        const int k = k_;
        for (int i = 0; i < 3; ++i) {
            const int vb = c1_vertex_dof(i, 0);
            x[vb] = y[vb];
            // x_g = B^{-T} y_g
            x[vb + 1] = binv_[0][0] * y[vb + 1] + binv_[1][0] * y[vb + 2];
            x[vb + 2] = binv_[0][1] * y[vb + 1] + binv_[1][1] * y[vb + 2];
            // x_h = Sym2(B^{-1}) y_h
            for (int r = 0; r < 3; ++r)
                x[vb + 3 + r] = sym2inv_[r][0] * y[vb + 3]
                                + sym2inv_[r][1] * y[vb + 4]
                                + sym2inv_[r][2] * y[vb + 5];
        }
        for (int r = 18 + 3 * c1_nedge(k); r < dim_; ++r) x[r] = y[r];
        for (int e = 0; e < 3; ++e)
            for (int i = 0; i < c1_ntrace(k); ++i) {
                int idx = c1_edge_trace_dof(k, e, i);
                x[idx] = y[idx];
            }
        for (int e = 0; e < 3; ++e) {
            const int s = (e + 1) % 3, q = (e + 2) % 3;
            hermite_data(e, s, q, x);          // from the recovered entries
            const T& boa = beta_over_alpha_[static_cast<std::size_t>(e)];
            for (int j = 0; j < c1_nnd(k); ++j) {
                T acc = W_->at(j, 0) * dbuf_[0];
                for (int r = 1; r <= k; ++r)
                    acc += W_->at(j, r) * dbuf_[static_cast<std::size_t>(r)];
                int idx = c1_edge_nd_dof(k, e, j);
                x[idx] = inv_alpha_[static_cast<std::size_t>(e)] * y[idx]
                         - boa * acc;
            }
        }
    }

    // partial inverse (the c1_elevate path): the physical ND entries of
    // local edge e from their reference values, using the ALREADY recovered
    // physical vertex / trace entries in x (vertex DOFs are copied there,
    // never re-evaluated -- C1-10)
    void invert_nd_edge(int e, const T* yhat_nd, T* x) const {
        require_built();
        assert(e >= 0 && e < 3);
        const int k = k_;
        hermite_data(e, (e + 1) % 3, (e + 2) % 3, x);
        const T& boa = beta_over_alpha_[static_cast<std::size_t>(e)];
        for (int j = 0; j < c1_nnd(k); ++j) {
            T acc = W_->at(j, 0) * dbuf_[0];
            for (int r = 1; r <= k; ++r)
                acc += W_->at(j, r) * dbuf_[static_cast<std::size_t>(r)];
            x[c1_edge_nd_dof(k, e, j)] =
                inv_alpha_[static_cast<std::size_t>(e)] * yhat_nd[j] - boa * acc;
        }
    }

    // ---- open materials (AXN-C1 in-test replication) ----
    T b_[2][2];
    T sym2_[3][3];
    T binv_[2][2];
    T sym2inv_[3][3];
    std::array<T, 3> alpha_, beta_, inv_alpha_, beta_over_alpha_;
    std::array<std::array<T, 2>, 3> tvec_;
    std::array<std::array<T, 3>, 3> ttvec_;

private:
    int k_, dim_;
    const rt_mat_tbl<T>* W_;
    mutable std::vector<T> dbuf_;              // Hermite data / scatter buffer

    void require_built() const {
        if (!W_)
            throw std::logic_error("bfem::c1_pullback: build not called");
    }

    static const T& ehalf() {                  // 1/|e_hat_0|^2 = 1/2, exact
        static const T c = rational_to<T>(1, 2);
        return c;
    }

    // Sym^2 rows of the map H_phys -> H_ref for the matrix M (reference
    // Hessian components (11, 12, 22) from physical (xx, xy, yy)):
    //   H_ref[cd] = sum_ab M[a][c] M[b][d] H_phys[ab]
    static void sym2_rows(const T m[2][2], T out[3][3]) {
        out[0][0] = m[0][0] * m[0][0];
        out[0][1] = m[0][0] * m[1][0] + m[0][0] * m[1][0];
        out[0][2] = m[1][0] * m[1][0];
        out[1][0] = m[0][0] * m[0][1];
        out[1][1] = m[0][0] * m[1][1] + m[0][1] * m[1][0];
        out[1][2] = m[1][0] * m[1][1];
        out[2][0] = m[0][1] * m[0][1];
        out[2][1] = m[0][1] * m[1][1] + m[0][1] * m[1][1];
        out[2][2] = m[1][1] * m[1][1];
    }

    // forward vertex blocks: value identity, gradient B^T, Hessian Sym2
    void vertex_forward(int i, const T* x, T* y) const {
        const int vb = c1_vertex_dof(i, 0);
        y[vb] = x[vb];
        y[vb + 1] = b_[0][0] * x[vb + 1] + b_[1][0] * x[vb + 2];
        y[vb + 2] = b_[0][1] * x[vb + 1] + b_[1][1] * x[vb + 2];
        for (int r = 0; r < 3; ++r)
            y[vb + 3 + r] = sym2_[r][0] * x[vb + 3]
                            + sym2_[r][1] * x[vb + 4]
                            + sym2_[r][2] * x[vb + 5];
    }

    // Hermite data of edge e into dbuf_: (p(0), p'(0), p''(0), p(1), p'(1),
    // p''(1), trace values), parameter derivatives via t_e / tt_e weights
    void hermite_data(int e, int s, int q, const T* x) const {
        const int k = k_;
        const T* t = tvec_[static_cast<std::size_t>(e)].data();
        const T* tt = ttvec_[static_cast<std::size_t>(e)].data();
        const int sb = c1_vertex_dof(s, 0);
        const int qb = c1_vertex_dof(q, 0);
        dbuf_[0] = x[sb];
        dbuf_[1] = t[0] * x[sb + 1] + t[1] * x[sb + 2];
        dbuf_[2] = tt[0] * x[sb + 3] + tt[1] * x[sb + 4] + tt[2] * x[sb + 5];
        dbuf_[3] = x[qb];
        dbuf_[4] = t[0] * x[qb + 1] + t[1] * x[qb + 2];
        dbuf_[5] = tt[0] * x[qb + 3] + tt[1] * x[qb + 4] + tt[2] * x[qb + 5];
        for (int i = 0; i < c1_ntrace(k); ++i)
            dbuf_[static_cast<std::size_t>(6 + i)] =
                x[c1_edge_trace_dof(k, e, i)];
    }

    // transpose scatter of the Hermite data weights (dbuf_ holds s_r)
    void scatter_data(int e, int s, int q, T* z) const {
        const int k = k_;
        const T* t = tvec_[static_cast<std::size_t>(e)].data();
        const T* tt = ttvec_[static_cast<std::size_t>(e)].data();
        const int sb = c1_vertex_dof(s, 0);
        const int qb = c1_vertex_dof(q, 0);
        z[sb] += dbuf_[0];
        z[sb + 1] += t[0] * dbuf_[1];
        z[sb + 2] += t[1] * dbuf_[1];
        z[sb + 3] += tt[0] * dbuf_[2];
        z[sb + 4] += tt[1] * dbuf_[2];
        z[sb + 5] += tt[2] * dbuf_[2];
        z[qb] += dbuf_[3];
        z[qb + 1] += t[0] * dbuf_[4];
        z[qb + 2] += t[1] * dbuf_[4];
        z[qb + 3] += tt[0] * dbuf_[5];
        z[qb + 4] += tt[1] * dbuf_[5];
        z[qb + 5] += tt[2] * dbuf_[5];
        for (int i = 0; i < c1_ntrace(k); ++i)
            z[c1_edge_trace_dof(k, e, i)] += dbuf_[static_cast<std::size_t>(6 + i)];
    }
};

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_C1_C1_GEOMETRY_HPP
