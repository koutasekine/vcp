// vcp/bfem/c1/c1_tables.hpp
// Phase 6 (2D C1 Argyris family): exact rational reference tables of the
// C1 element for general k >= 5 -- DOF application matrix, dual basis,
// contracted stiffness / second-derivative (R2) / mass / load tables, the
// 1D Hermite inverse table (T-C3), and the two-stage registry
// c1_registry / c1_typed_registry<T>.
//
// Conforms to: C1 external design v0.2 (sections 2, 3, 5) and
//              C1 internal design v0.2 (section 2).
//
// Normative conventions implemented here:
//  - reference vertices (0,0), (1,0), (0,1); edge e runs (e+1)%3 -> (e+2)%3;
//    nu_hat_e = R_{-90}(edge vector) and grad-hat lambda are the RT constants
//    (rt_edge_normal / rt_grad_lambda -- rt_tables machinery reuse, C1-2);
//  - reference DOFs are built from REFERENCE CARTESIAN derivatives (v0.2,
//    A-1): d/dx_hat_d = sum_i (grad-hat lambda_i)_d d/dlambda_i with integer
//    grad-hat lambda, realized as derivative_map index shifts with integer
//    factors (the RT div construction transplanted; Hessian = two shifts);
//  - local DOF order (normative): vertex blocks (value, d1, d2, d11, d12,
//    d22 for local vertices 0, 1, 2), edge blocks (trace point values
//    t = i/(k-4), i = 1..k-5, then normal-derivative point values
//    t = j/(k-3), j = 1..k-4, local edges 0, 1, 2), interior point values
//    (degree k-3 lattice interior points, canonical L0 order);
//  - point evaluation is the closed form B^k_alpha(p) = multinomial(k,
//    alpha) * lambda(p)^alpha (v0.2, B-1 -- the rational-stage point
//    evaluation helper that the rt machinery does not have; the 1D version
//    B^k_j(t) = C(k,j) t^j (1-t)^{k-j} serves the edge rows and T-C3);
//  - all geometry-free contractions happen at the RATIONAL stage; the typed
//    stage is a pure convert_traits image (enclose-once, S-C1-3);
//  - generation-time duality checks: L C == I (element Vandermonde) and
//    H H^{-1} == I (T-C3 Hermite table) raise std::logic_error on failure.
//
// No floating point appears in this header. No division operator appears in
// this header (G-C1-1: the single c1/ division is the edge cache of
// c1_geometry.hpp); rational divisions run inside the frozen rational_la.

#ifndef VCP_BFEM_C1_C1_TABLES_HPP
#define VCP_BFEM_C1_C1_TABLES_HPP

#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/ref_stiffness.hpp>       // detail::ref_stiff_block_index
#include <vcp/bfem/rt/rational_la.hpp>
#include <vcp/bfem/rt/rt_tables.hpp>        // table shells + rt constants

namespace vcp {
namespace bfem {

namespace detail {

// ---------------------------------------------------------------------------
// C1 reference constants and dimension bookkeeping (k >= 5 throughout)
// ---------------------------------------------------------------------------

inline void c1_check_k(int k, const char* where) {
    if (k < 5) {
        std::string msg("bfem::c1: ");
        msg += where;
        msg += ": k must be >= 5";
        throw std::invalid_argument(msg);
    }
}

// |e_hat|^2 of reference edge e (edge 0: (-1,1), edges 1, 2: unit axes)
inline int c1_edge_sq(int e) {
    assert(e >= 0 && e < 3);
    return e == 0 ? 2 : 1;
}

inline int c1_dim(int k) { return coeff_registry<2>::indices(k).size(); }
inline int c1_ntrace(int k) { return k - 5; }        // trace points per edge
inline int c1_nnd(int k) { return k - 4; }           // ND points per edge
inline int c1_nint(int k) { return ((k - 4) * (k - 5)) >> 1; }
inline int c1_nedge(int k) { return 2 * k - 9; }     // per-edge dof count

// local DOF indices (the normative local order above)
inline int c1_vertex_dof(int i, int c) {
    assert(i >= 0 && i < 3 && c >= 0 && c < 6);
    return 6 * i + c;
}
inline int c1_edge_trace_dof(int k, int e, int i) {  // i in 0..k-6
    assert(e >= 0 && e < 3 && i >= 0 && i < c1_ntrace(k));
    return 18 + e * c1_nedge(k) + i;
}
inline int c1_edge_nd_dof(int k, int e, int j) {     // j in 0..k-5
    assert(e >= 0 && e < 3 && j >= 0 && j < c1_nnd(k));
    return 18 + e * c1_nedge(k) + c1_ntrace(k) + j;
}
inline int c1_interior_dof(int k, int r) {
    assert(r >= 0 && r < c1_nint(k));
    return 18 + 3 * c1_nedge(k) + r;
}

// interior lattice indices of degree k-3 (all components >= 1), canonical
// L0 order; the DOF point of entry alpha is alpha / (k-3)
inline std::vector<multi_index<2> > c1_interior_indices(int k) {
    const index_map<2>& im = coeff_registry<2>::indices(k - 3);
    std::vector<multi_index<2> > out;
    for (int r = 0; r < im.size(); ++r) {
        multi_index<2> al = im.unrank(r);
        if (al.a[0] >= 1 && al.a[1] >= 1 && al.a[2] >= 1) out.push_back(al);
    }
    assert(static_cast<int>(out.size()) == c1_nint(k));
    return out;
}

// ---------------------------------------------------------------------------
// rational-stage point evaluation helpers (v0.2, B-1)
// ---------------------------------------------------------------------------

inline rational c1_pow(const rational& x, int p) {
    rational r(1);
    for (int i = 0; i < p; ++i) r *= x;
    return r;
}

// 2D: B^n_alpha(p) = multinomial(n, alpha) * prod_i lambda_i(p)^alpha_i
inline rational c1_eval_b2(int n, const multi_index<2>& al,
                           const std::array<rational, 3>& lam) {
    rational r(multinomial(n, al.a.data(), 3), bigint(1));
    for (int i = 0; i < 3; ++i)
        r *= c1_pow(lam[static_cast<std::size_t>(i)],
                    al.a[static_cast<std::size_t>(i)]);
    return r;
}

// 1D: B^n_j(t) = C(n, j) t^j (1-t)^(n-j)
inline rational c1_eval_b1(int n, int j, const rational& t) {
    assert(j >= 0 && j <= n);
    rational r(binomial_cache::binom(n, j), bigint(1));
    r *= c1_pow(t, j);
    r *= c1_pow(rational(1) - t, n - j);
    return r;
}

// barycentric point of edge e at parameter t (local direction (e+1)%3 ->
// (e+2)%3): lambda_e = 0, lambda_start = 1-t, lambda_end = t
inline std::array<rational, 3> c1_edge_point(int e, const rational& t) {
    std::array<rational, 3> lam;
    lam[static_cast<std::size_t>(e)] = rational(0);
    lam[static_cast<std::size_t>((e + 1) % 3)] = rational(1) - t;
    lam[static_cast<std::size_t>((e + 2) % 3)] = t;
    return lam;
}

// ---------------------------------------------------------------------------
// DOF row builders on degree-k Bernstein coefficient vectors
// ---------------------------------------------------------------------------

// row += (value at lam)
inline void c1_row_value(rmat& L, int row, int k,
                         const std::array<rational, 3>& lam) {
    const index_map<2>& im = coeff_registry<2>::indices(k);
    for (int r = 0; r < im.size(); ++r)
        L.at(row, r) += c1_eval_b2(k, im.unrank(r), lam);
}

// row += (w . grad-hat) u at lam  (w integer Cartesian weight vector):
// (w . grad-hat u)(p) = k sum_i g_i sum_{|gam|=k-1} b_{gam+e_i} B^{k-1}_gam(p)
// with g_i = sum_d w_d (grad-hat lambda_i)_d (integer)
inline void c1_row_cart1(rmat& L, int row, int k, int w0, int w1,
                         const std::array<rational, 3>& lam) {
    const index_map<2>& imk = coeff_registry<2>::indices(k);
    const index_map<2>& imm = coeff_registry<2>::indices(k - 1);
    for (int r = 0; r < imm.size(); ++r) {
        multi_index<2> gam = imm.unrank(r);
        rational v = c1_eval_b2(k - 1, gam, lam);
        if (v.is_zero()) continue;
        for (int i = 0; i < 3; ++i) {
            int g = w0 * rt_grad_lambda(i, 0) + w1 * rt_grad_lambda(i, 1);
            if (g == 0) continue;
            multi_index<2> be = gam;
            be.a[static_cast<std::size_t>(i)] += 1;
            L.at(row, imk.rank(be)) += rational(g * k) * v;
        }
    }
}

// row += (d/dx_hat_c d/dx_hat_d) u at lam
inline void c1_row_cart2(rmat& L, int row, int k, int c, int d,
                         const std::array<rational, 3>& lam) {
    const index_map<2>& imk = coeff_registry<2>::indices(k);
    const index_map<2>& imm = coeff_registry<2>::indices(k - 2);
    const int kk1 = k * (k - 1);
    for (int r = 0; r < imm.size(); ++r) {
        multi_index<2> gam = imm.unrank(r);
        rational v = c1_eval_b2(k - 2, gam, lam);
        if (v.is_zero()) continue;
        for (int i = 0; i < 3; ++i) {
            int gi = rt_grad_lambda(i, c);
            if (gi == 0) continue;
            for (int j = 0; j < 3; ++j) {
                int gj = rt_grad_lambda(j, d);
                if (gj == 0) continue;
                multi_index<2> be = gam;
                be.a[static_cast<std::size_t>(i)] += 1;
                be.a[static_cast<std::size_t>(j)] += 1;
                L.at(row, imk.rank(be)) += rational(gi * gj * kk1) * v;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// c1_dof_matrix: the DOF application matrix L (dim x N_k) acting on degree-k
// Bernstein coefficients, rows in the normative local DOF order. The element
// Vandermonde IS L (the spanning set is the P_k identity, C1-2); the dual
// basis is C = L^{-1}.
// ---------------------------------------------------------------------------
inline rmat c1_dof_matrix(int k) {
    c1_check_k(k, "c1_dof_matrix");
    const int dim = c1_dim(k);
    rmat L(dim, dim);
    int row = 0;
    // vertex blocks: value, d1, d2, d11, d12, d22 (reference Cartesian)
    for (int i = 0; i < 3; ++i) {
        std::array<rational, 3> lam;
        lam[0] = rational(0);
        lam[1] = rational(0);
        lam[2] = rational(0);
        lam[static_cast<std::size_t>(i)] = rational(1);
        c1_row_value(L, row++, k, lam);
        c1_row_cart1(L, row++, k, 1, 0, lam);
        c1_row_cart1(L, row++, k, 0, 1, lam);
        c1_row_cart2(L, row++, k, 0, 0, lam);
        c1_row_cart2(L, row++, k, 0, 1, lam);
        c1_row_cart2(L, row++, k, 1, 1, lam);
    }
    // edge blocks: trace point values then normal-derivative point values
    for (int e = 0; e < 3; ++e) {
        for (int i = 1; i <= c1_ntrace(k); ++i) {
            std::array<rational, 3> lam = c1_edge_point(e, rational(i, k - 4));
            c1_row_value(L, row++, k, lam);
        }
        for (int j = 1; j <= c1_nnd(k); ++j) {
            std::array<rational, 3> lam = c1_edge_point(e, rational(j, k - 3));
            c1_row_cart1(L, row++, k, rt_edge_normal(e, 0),
                         rt_edge_normal(e, 1), lam);
        }
    }
    // interior point values
    std::vector<multi_index<2> > ints = c1_interior_indices(k);
    for (std::size_t r = 0; r < ints.size(); ++r) {
        std::array<rational, 3> lam;
        for (int i = 0; i < 3; ++i)
            lam[static_cast<std::size_t>(i)] =
                rational(ints[r].a[static_cast<std::size_t>(i)], k - 3);
        c1_row_value(L, row++, k, lam);
    }
    assert(row == dim);
    return L;
}

// ---------------------------------------------------------------------------
// 1D Hermite table (T-C3, v0.2 B-2): the data functionals of the edge trace
//   d = (p(0), p'(0), p''(0), p(1), p'(1), p''(1), p(i/(k-4)) i = 1..k-5)
// determine the degree-k trace polynomial; W (k-4 x k+1) gives the tangential
// parameter derivative at the ND points: p'(j/(k-3)) = sum_r W[j-1][r] d_r.
// Geometry free, exact, one per k; generation-time duality check H H^{-1} = I.
// ---------------------------------------------------------------------------

// 1D derivative rows over degree-k 1D Bernstein coefficients
inline void c1_row1_value(rmat& H, int row, int k, const rational& t) {
    for (int j = 0; j <= k; ++j) H.at(row, j) += c1_eval_b1(k, j, t);
}
inline void c1_row1_d1(rmat& H, int row, int k, const rational& t) {
    for (int g = 0; g <= k - 1; ++g) {
        rational b = rational(k) * c1_eval_b1(k - 1, g, t);
        if (b.is_zero()) continue;
        H.at(row, g + 1) += b;
        H.at(row, g) -= b;
    }
}
inline void c1_row1_d2(rmat& H, int row, int k, const rational& t) {
    const rational kk1(k * (k - 1));
    for (int g = 0; g <= k - 2; ++g) {
        rational b = kk1 * c1_eval_b1(k - 2, g, t);
        if (b.is_zero()) continue;
        H.at(row, g + 2) += b;
        H.at(row, g + 1) -= b + b;
        H.at(row, g) += b;
    }
}

inline rmat c1_hermite_matrix(int k) {
    c1_check_k(k, "c1_hermite_matrix");
    rmat H(k + 1, k + 1);
    int row = 0;
    c1_row1_value(H, row++, k, rational(0));
    c1_row1_d1(H, row++, k, rational(0));
    c1_row1_d2(H, row++, k, rational(0));
    c1_row1_value(H, row++, k, rational(1));
    c1_row1_d1(H, row++, k, rational(1));
    c1_row1_d2(H, row++, k, rational(1));
    for (int i = 1; i <= k - 5; ++i)
        c1_row1_value(H, row++, k, rational(i, k - 4));
    assert(row == k + 1);
    return H;
}

inline rmat c1_hermite_w(int k) {
    rmat H = c1_hermite_matrix(k);
    rmat Hkeep = H;                                  // solve_exact consumes H
    rmat Hinv = solve_exact(std::move(H), rmat::identity(k + 1));
    // generation-time duality check (T-C3)
    rmat P = mul(Hkeep, Hinv);
    for (int i = 0; i <= k; ++i)
        for (int j = 0; j <= k; ++j)
            if (!(P.at(i, j) == rational(i == j ? 1 : 0)))
                throw std::logic_error(
                    "bfem::c1_hermite_w: Hermite duality check failed");
    rmat G(k - 4, k + 1);
    for (int j = 1; j <= k - 4; ++j)
        c1_row1_d1(G, j - 1, k, rational(j, k - 3));
    return mul(G, Hinv);
}

// ---------------------------------------------------------------------------
// c1_pair_tbl<S>: 21 stored blocks over pairs (p <= q) of the 6 second-
// derivative directions (p, q index the barycentric pairs (i <= j) through
// ref_stiff_block_index); block(q, p) is the transposed view of (p, q).
// ---------------------------------------------------------------------------
template <typename S>
class c1_pair_tbl {
public:
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    rt_block_view<S> block(int p, int q) const {
        assert(p >= 0 && p < 6 && q >= 0 && q < 6);
        const bool tr = p > q;
        const int lo = tr ? q : p;
        const int hi = tr ? p : q;
        static const int off[6] = { 0, 6, 11, 15, 18, 20 };
        return rt_block_view<S>(
            blk_[static_cast<std::size_t>(off[lo] + (hi - lo))].data(),
            rows_, cols_, tr);
    }

private:
    friend struct c1_table_access;
    c1_pair_tbl() : rows_(0), cols_(0), blk_() {}
    int rows_, cols_;
    std::vector<std::vector<S> > blk_;               // 21 blocks, upper pairs
};

struct c1_table_access {
    template <typename S>
    static c1_pair_tbl<S> make_pair(int rows, int cols,
                                    std::vector<std::vector<S> > blk) {
        assert(blk.size() == 21u);
        c1_pair_tbl<S> t;
        t.rows_ = rows;
        t.cols_ = cols;
        t.blk_ = std::move(blk);
        return t;
    }
};

// second-derivative direction pair p <-> (i <= j) over 0..2
inline void c1_pair_dirs(int p, int& i, int& j) {
    assert(p >= 0 && p < 6);
    static const int pi[6] = { 0, 0, 0, 1, 1, 2 };
    static const int pj[6] = { 0, 1, 2, 1, 2, 2 };
    i = pi[p];
    j = pj[p];
}

// D2C^{(p)} (N_{k-2} x dim): degree k-2 Bernstein coefficients of
// d/dlambda_i d/dlambda_j Phi_a, entries k (k-1) C[gamma + e_i + e_j, a]
inline rmat c1_d2c_block(int k, const rmat& C, int p) {
    const index_map<2>& imk = coeff_registry<2>::indices(k);
    const index_map<2>& im2 = coeff_registry<2>::indices(k - 2);
    const int dim = C.cols;
    int di, dj;
    c1_pair_dirs(p, di, dj);
    const rational kk1(k * (k - 1));
    rmat D(im2.size(), dim);
    for (int g = 0; g < im2.size(); ++g) {
        multi_index<2> be = im2.unrank(g);
        be.a[static_cast<std::size_t>(di)] += 1;
        be.a[static_cast<std::size_t>(dj)] += 1;
        int src = imk.rank(be);
        for (int a = 0; a < dim; ++a) {
            const rational& c = C.at(src, a);
            if (c.is_zero()) continue;
            D.at(g, a) = kk1 * c;
        }
    }
    return D;
}

} // namespace detail

// ---------------------------------------------------------------------------
// c1_registry: rational stage. Lazy generation under one mutex; returned
// references stay valid until program termination (the L0 cache contract).
// ---------------------------------------------------------------------------
class c1_registry {
public:
    static int dim(int k) {
        detail::c1_check_k(k, "dim");
        return detail::c1_dim(k);
    }

    // DOF application matrix L (dim x N_k)
    static const rt_mat_table& dof_matrix(int k) {
        detail::c1_check_k(k, "dof_matrix");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return dofm_locked(s, k);
    }

    // dual basis C = L^{-1} (N_k x dim; generation checks L C == I)
    static const rt_mat_table& basis(int k) {
        detail::c1_check_k(k, "basis");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return basis_locked(s, k);
    }

    // contracted first-derivative stiffness blocks S^{(ij)} = C^T R^{(ij)} C
    // (6 stored blocks over barycentric i <= j; block(j, i) is the transpose)
    static const rt_block_table& stiffness(int k) {
        detail::c1_check_k(k, "stiffness");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, rt_block_table>::iterator it = s.stiff.find(k);
        if (it != s.stiff.end()) return it->second;
        rt_block_table t = build_stiffness(s, k);
        return s.stiff.insert(std::make_pair(k, std::move(t))).first->second;
    }

    // contracted second-derivative blocks R2^{(p),(q)} =
    // (D2C^{(p)})^T M^{(k-2,k-2)} D2C^{(q)} over direction pairs p <= q
    static const detail::c1_pair_tbl<detail::rational>& r2(int k) {
        detail::c1_check_k(k, "r2");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, detail::c1_pair_tbl<detail::rational> >::iterator
            it = s.r2.find(k);
        if (it != s.r2.end()) return it->second;
        detail::c1_pair_tbl<detail::rational> t = build_r2(s, k);
        return s.r2.insert(std::make_pair(k, std::move(t))).first->second;
    }

    // D2C^{(p)} blocks (N_{k-2} x dim), stored as a 6-block table indexed by
    // the barycentric pair (i <= j) through block(i, j)
    static const rt_block_table& d2c(int k) {
        detail::c1_check_k(k, "d2c");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return d2c_locked(s, k);
    }

    // contracted mass C_a^T M^{(a,b)} C_b (dim_a x dim_b), a, b >= 5
    static const rt_mat_table& mass_contracted(int a, int b) {
        detail::c1_check_k(a, "mass_contracted");
        detail::c1_check_k(b, "mass_contracted");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        typename std::map<std::pair<int, int>, rt_mat_table>::iterator it =
            s.massc.find(key);
        if (it != s.massc.end()) return it->second;
        rt_mat_table t = build_mass_contracted(s, a, b);
        return s.massc.insert(std::make_pair(key, std::move(t))).first->second;
    }

    // contracted load C^T M^{(k,q)} (dim x N_q), q >= 0
    static const rt_mat_table& load_contracted(int k, int q) {
        detail::c1_check_k(k, "load_contracted");
        if (q < 0)
            throw std::invalid_argument("bfem::c1_registry::load_contracted: q < 0");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, q);
        typename std::map<std::pair<int, int>, rt_mat_table>::iterator it =
            s.loadc.find(key);
        if (it != s.loadc.end()) return it->second;
        rt_mat_table t = build_load_contracted(s, k, q);
        return s.loadc.insert(std::make_pair(key, std::move(t))).first->second;
    }

    // laplacian mixed blocks LM^{(p)} = M^{(l,k-2)} D2C^{(p)} (N_l x dim),
    // stored as a 6-block table indexed through block(i, j)
    static const rt_block_table& lap_mixed(int k, int l) {
        detail::c1_check_k(k, "lap_mixed");
        if (l < 0)
            throw std::invalid_argument("bfem::c1_registry::lap_mixed: l < 0");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, l);
        typename std::map<std::pair<int, int>, rt_block_table>::iterator it =
            s.lapm.find(key);
        if (it != s.lapm.end()) return it->second;
        rt_block_table t = build_lap_mixed(s, k, l);
        return s.lapm.insert(std::make_pair(key, std::move(t))).first->second;
    }

    // T-C3: tangential-derivative combination table W ((k-4) x (k+1))
    static const rt_mat_table& hermite_w(int k) {
        detail::c1_check_k(k, "hermite_w");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, rt_mat_table>::iterator it = s.herm.find(k);
        if (it != s.herm.end()) return it->second;
        detail::rmat W = detail::c1_hermite_w(k);
        rt_mat_table t = detail::rt_table_access::make_mat(
            W.rows, W.cols, detail::rmat_flat(W));
        return s.herm.insert(std::make_pair(k, std::move(t))).first->second;
    }

private:
    struct state {
        std::mutex mtx;
        std::map<int, rt_mat_table> dofm;
        std::map<int, rt_mat_table> basis;
        std::map<int, rt_block_table> stiff;
        std::map<int, detail::c1_pair_tbl<detail::rational> > r2;
        std::map<int, rt_block_table> d2c;
        std::map<std::pair<int, int>, rt_mat_table> massc;
        std::map<std::pair<int, int>, rt_mat_table> loadc;
        std::map<std::pair<int, int>, rt_block_table> lapm;
        std::map<int, rt_mat_table> herm;
    };
    static state& st() {
        static state s;
        return s;
    }

    static const rt_mat_table& dofm_locked(state& s, int k) {
        typename std::map<int, rt_mat_table>::iterator it = s.dofm.find(k);
        if (it != s.dofm.end()) return it->second;
        detail::rmat L = detail::c1_dof_matrix(k);
        rt_mat_table t = detail::rt_table_access::make_mat(
            L.rows, L.cols, detail::rmat_flat(L));
        return s.dofm.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static const rt_mat_table& basis_locked(state& s, int k) {
        typename std::map<int, rt_mat_table>::iterator it = s.basis.find(k);
        if (it != s.basis.end()) return it->second;
        const rt_mat_table& Lt = dofm_locked(s, k);
        const int dim = Lt.rows();
        detail::rmat L(dim, dim);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                L.at(i, j) = Lt.at(i, j);
        detail::rmat Lkeep = L;
        detail::rmat C = detail::solve_exact(std::move(L),
                                             detail::rmat::identity(dim));
        // generation-time duality check (C-T1 at the source)
        detail::rmat P = detail::mul(Lkeep, C);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                if (!(P.at(i, j) == detail::rational(i == j ? 1 : 0)))
                    throw std::logic_error(
                        "bfem::c1_registry::basis: duality check failed");
        rt_mat_table t = detail::rt_table_access::make_mat(
            C.rows, C.cols, detail::rmat_flat(C));
        return s.basis.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static detail::rmat basis_rmat_locked(state& s, int k) {
        const rt_mat_table& Ct = basis_locked(s, k);
        detail::rmat C(Ct.rows(), Ct.cols());
        for (int i = 0; i < Ct.rows(); ++i)
            for (int j = 0; j < Ct.cols(); ++j)
                C.at(i, j) = Ct.at(i, j);
        return C;
    }

    static rt_block_table build_stiffness(state& s, int k) {
        using detail::rmat;
        const rmat C = basis_rmat_locked(s, k);
        const int dim = C.cols;
        const int N = C.rows;
        const mass_table<2>& M = coeff_registry<2>::mass(k - 1, k - 1);
        derivative_map<2> dm(k);
        const detail::rational k2(k * k);
        rmat Ct = detail::transpose(C);
        std::vector<std::vector<detail::rational> > blk;
        blk.reserve(6u);
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                // R^{(ij)}[a][b] = k^2 M^{k-1,k-1}[dm(a,i), dm(b,j)]
                rmat R(N, N);
                for (int a = 0; a < N; ++a) {
                    int ta = dm.target(a, i);
                    if (ta < 0) continue;
                    for (int b = 0; b < N; ++b) {
                        int tb = dm.target(b, j);
                        if (tb < 0) continue;
                        R.at(a, b) = k2 * M.at(ta, tb);
                    }
                }
                blk.push_back(detail::rmat_flat(
                    detail::mul(Ct, detail::mul(R, C))));
            }
        }
        return detail::rt_table_access::make_block(dim, dim, std::move(blk));
    }

    static const rt_block_table& d2c_locked(state& s, int k) {
        typename std::map<int, rt_block_table>::iterator it = s.d2c.find(k);
        if (it != s.d2c.end()) return it->second;
        using detail::rmat;
        const rmat C = basis_rmat_locked(s, k);
        std::vector<std::vector<detail::rational> > blk;
        blk.reserve(6u);
        int rows = 0;
        for (int p = 0; p < 6; ++p) {
            rmat D = detail::c1_d2c_block(k, C, p);
            rows = D.rows;
            blk.push_back(detail::rmat_flat(D));
        }
        rt_block_table t = detail::rt_table_access::make_block(
            rows, C.cols, std::move(blk));
        return s.d2c.insert(std::make_pair(k, std::move(t))).first->second;
    }

    static detail::c1_pair_tbl<detail::rational> build_r2(state& s, int k) {
        using detail::rmat;
        const rmat C = basis_rmat_locked(s, k);
        const int dim = C.cols;
        const mass_table<2>& Mt = coeff_registry<2>::mass(k - 2, k - 2);
        rmat M(Mt.rows(), Mt.cols());
        for (int i = 0; i < Mt.rows(); ++i)
            for (int j = 0; j < Mt.cols(); ++j)
                M.at(i, j) = Mt.at(i, j);
        std::vector<rmat> D;
        D.reserve(6u);
        for (int p = 0; p < 6; ++p)
            D.push_back(detail::c1_d2c_block(k, C, p));
        std::vector<std::vector<detail::rational> > blk;
        blk.reserve(21u);
        for (int p = 0; p < 6; ++p) {
            rmat DtM = detail::mul(detail::transpose(D[static_cast<std::size_t>(p)]), M);
            for (int q = p; q < 6; ++q)
                blk.push_back(detail::rmat_flat(
                    detail::mul(DtM, D[static_cast<std::size_t>(q)])));
        }
        return detail::c1_table_access::make_pair<detail::rational>(
            dim, dim, std::move(blk));
    }

    static rt_mat_table build_mass_contracted(state& s, int a, int b) {
        using detail::rmat;
        const rmat Ca = basis_rmat_locked(s, a);
        const rmat Cb = basis_rmat_locked(s, b);
        rmat M = detail::mass_rmat<2>(a, b);
        rmat R = detail::mul(detail::transpose(Ca), detail::mul(M, Cb));
        return detail::rt_table_access::make_mat(R.rows, R.cols,
                                                 detail::rmat_flat(R));
    }

    static rt_mat_table build_load_contracted(state& s, int k, int q) {
        using detail::rmat;
        const rmat C = basis_rmat_locked(s, k);
        rmat M = detail::mass_rmat<2>(k, q);
        rmat R = detail::mul(detail::transpose(C), M);
        return detail::rt_table_access::make_mat(R.rows, R.cols,
                                                 detail::rmat_flat(R));
    }

    static rt_block_table build_lap_mixed(state& s, int k, int l) {
        using detail::rmat;
        const rt_block_table& D2 = d2c_locked(s, k);
        rmat M = detail::mass_rmat<2>(l, k - 2);
        std::vector<std::vector<detail::rational> > blk;
        blk.reserve(6u);
        int rows = 0, cols = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                rt_block_view<detail::rational> Dv = D2.block(i, j);
                rmat D(Dv.rows(), Dv.cols());
                for (int r = 0; r < Dv.rows(); ++r)
                    for (int c = 0; c < Dv.cols(); ++c)
                        D.at(r, c) = Dv.at(r, c);
                rmat R = detail::mul(M, D);
                rows = R.rows;
                cols = R.cols;
                blk.push_back(detail::rmat_flat(R));
            }
        }
        return detail::rt_table_access::make_block(rows, cols, std::move(blk));
    }
};

// ---------------------------------------------------------------------------
// c1_typed_registry<T>: T stage. Every entry is the convert_traits image of
// the rational-stage value exactly once (enclose-once, S-C1-3); no T
// arithmetic beyond the conversion occurs here.
// ---------------------------------------------------------------------------
template <typename T>
class c1_typed_registry {
public:
    typedef rt_mat_tbl<T> mat_table;
    typedef rt_block_tbl<T> block_table;
    typedef detail::c1_pair_tbl<T> pair_table;

    static const mat_table& dof_matrix(int k) {
        return mat_entry(st().dofm, k, c1_registry::dof_matrix(k));
    }
    static const mat_table& basis(int k) {
        return mat_entry(st().basis, k, c1_registry::basis(k));
    }
    static const mat_table& hermite_w(int k) {
        return mat_entry(st().herm, k, c1_registry::hermite_w(k));
    }
    static const block_table& stiffness(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, block_table>::iterator it = s.stiff.find(k);
        if (it != s.stiff.end()) return it->second;
        block_table t = conv_block6(c1_registry::stiffness(k));
        return s.stiff.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const block_table& d2c(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, block_table>::iterator it = s.d2c.find(k);
        if (it != s.d2c.end()) return it->second;
        block_table t = conv_block6(c1_registry::d2c(k));
        return s.d2c.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const pair_table& r2(int k) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, pair_table>::iterator it = s.r2.find(k);
        if (it != s.r2.end()) return it->second;
        const detail::c1_pair_tbl<detail::rational>& src = c1_registry::r2(k);
        std::vector<std::vector<T> > blk;
        blk.reserve(21u);
        for (int p = 0; p < 6; ++p) {
            for (int q = p; q < 6; ++q) {
                rt_block_view<detail::rational> V = src.block(p, q);
                std::vector<T> w;
                w.reserve(static_cast<std::size_t>(V.rows())
                          * static_cast<std::size_t>(V.cols()));
                for (int i = 0; i < V.rows(); ++i)
                    for (int j = 0; j < V.cols(); ++j)
                        w.push_back(conv(V.at(i, j)));
                blk.push_back(std::move(w));
            }
        }
        pair_table t = detail::c1_table_access::make_pair<T>(
            src.rows(), src.cols(), std::move(blk));
        return s.r2.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const mat_table& mass_contracted(int a, int b) {
        return mat_entry2(st().massc, a, b, c1_registry::mass_contracted(a, b));
    }
    static const mat_table& load_contracted(int k, int q) {
        return mat_entry2(st().loadc, k, q, c1_registry::load_contracted(k, q));
    }
    static const block_table& lap_mixed(int k, int l) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, l);
        typename std::map<std::pair<int, int>, block_table>::iterator it =
            s.lapm.find(key);
        if (it != s.lapm.end()) return it->second;
        block_table t = conv_block6(c1_registry::lap_mixed(k, l));
        return s.lapm.insert(std::make_pair(key, std::move(t))).first->second;
    }

private:
    struct state {
        std::mutex mtx;
        std::map<int, mat_table> dofm;
        std::map<int, mat_table> basis;
        std::map<int, mat_table> herm;
        std::map<int, block_table> stiff;
        std::map<int, block_table> d2c;
        std::map<int, pair_table> r2;
        std::map<std::pair<int, int>, mat_table> massc;
        std::map<std::pair<int, int>, mat_table> loadc;
        std::map<std::pair<int, int>, block_table> lapm;
    };
    static state& st() {
        static state s;
        return s;
    }
    static T conv(const detail::rational& r) {
        return convert_traits<T>::from_rational(r.num(), r.den());
    }
    static mat_table conv_mat(const rt_mat_table& src) {
        std::vector<T> v;
        v.reserve(static_cast<std::size_t>(src.rows())
                  * static_cast<std::size_t>(src.cols()));
        for (int i = 0; i < src.rows(); ++i)
            for (int j = 0; j < src.cols(); ++j)
                v.push_back(conv(src.at(i, j)));
        return detail::rt_table_access::make_mat(src.rows(), src.cols(),
                                                 std::move(v));
    }
    static block_table conv_block6(const rt_block_table& src) {
        std::vector<std::vector<T> > blk;
        blk.reserve(6u);
        for (int i = 0; i < 3; ++i) {
            for (int j = i; j < 3; ++j) {
                rt_block_view<detail::rational> V = src.block(i, j);
                std::vector<T> w;
                w.reserve(static_cast<std::size_t>(V.rows())
                          * static_cast<std::size_t>(V.cols()));
                for (int r = 0; r < V.rows(); ++r)
                    for (int c = 0; c < V.cols(); ++c)
                        w.push_back(conv(V.at(r, c)));
                blk.push_back(std::move(w));
            }
        }
        return detail::rt_table_access::make_block(src.block_rows(),
                                                   src.block_cols(),
                                                   std::move(blk));
    }
    static const mat_table& mat_entry(std::map<int, mat_table>& m, int k,
                                      const rt_mat_table& src) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, mat_table>::iterator it = m.find(k);
        if (it != m.end()) return it->second;
        mat_table t = conv_mat(src);
        return m.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const mat_table& mat_entry2(std::map<std::pair<int, int>, mat_table>& m,
                                       int a, int b, const rt_mat_table& src) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        typename std::map<std::pair<int, int>, mat_table>::iterator it =
            m.find(key);
        if (it != m.end()) return it->second;
        mat_table t = conv_mat(src);
        return m.insert(std::make_pair(key, std::move(t))).first->second;
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_C1_C1_TABLES_HPP
