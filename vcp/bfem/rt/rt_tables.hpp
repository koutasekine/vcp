// vcp/bfem/rt/rt_tables.hpp
// RT Layer 0: exact rational coefficient tables of the reference Raviart-
// Thomas space RT_k (T-R1..T-R7) and the rational-stage registry
// rt_registry<D> (D == 2).
//
// Conforms to: RT-L0 external design v0.3 (sections 1, 3, 5) and
//              RT-L0 internal design v0.3 (sections 3, 4, 5, 7).
//
// The frozen P^n L0 files are reused read-only (bigint / rational /
// binomial_cache / index_map / coeff_registry); nothing outside vcp/bfem/rt/
// is modified. No floating point appears in this header. Rational division
// appears only inside detail::rational operations (generation time).
//
// Normative conventions implemented here (external design 1.1):
//  - embedding RT_k subset (P_{k+1})^2, component-major degree k+1 Bernstein
//    coefficients (vector length 2 N_{k+1});
//  - DOF order: edge block (e = 0,1,2; j = 0..k along the local positive
//    direction (e+1)%3 -> (e+2)%3) then interior block (component major,
//    canonical L0 order of index_map(k-1));
//  - nu_e = R_{-90}(end - start): edge-length absorbed rational normal;
//  - area tables are "coefficients of |T|" (the L0 mass convention);
//    the 1D edge tables are plain [0,1] parameter integrals.
//  - interior moments are the plain integrals int_T sigma_d B^{k-1}_alpha dx
//    of the external design 1.1; the |T_hat| = 1/D! measure factor is
//    therefore applied on top of the L0 "coefficient of |T|" mass value.
//    (The internal design 4.2 recipe omits this constant; the external
//    normative DOF definition prevails -- recorded in the gate report.)

#ifndef VCP_BFEM_RT_RT_TABLES_HPP
#define VCP_BFEM_RT_RT_TABLES_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/rt/rational_la.hpp>

namespace vcp {
namespace bfem {

namespace detail {
struct rt_table_access;   // the single construction/fill point of all tables
}

// ---------------------------------------------------------------------------
// Table shells, generic in the scalar S. The rational stage instantiates
// S = detail::rational; the typed stage (rt_typed_tables.hpp) instantiates
// S = T with the identical accessor surface (external design section 3).
// All tables are immutable after construction; accessors return const
// references only (S-RT0-3: no arithmetic, no allocation on the read path).
// ---------------------------------------------------------------------------

// generic dense matrix table (T-R6' / T-R7)
template <typename S>
class rt_mat_tbl {
public:
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    const S& at(int i, int j) const {
        assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
        return v_[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)
                  + static_cast<std::size_t>(j)];
    }

private:
    friend struct detail::rt_table_access;
    rt_mat_tbl() : rows_(0), cols_(0), v_() {}
    int rows_, cols_;
    std::vector<S> v_;
};

// T-R1: basis coefficient matrix C, 2 N_{k+1} x dim, component major
template <typename S>
class rt_basis_tbl {
public:
    int order() const { return k_; }            // k
    int comp_size() const { return nc_; }       // N_{k+1}
    int dim() const { return dim_; }            // (k+1)(k+3)
    int rows() const { return 2 * nc_; }
    int cols() const { return dim_; }
    const S& at(int i, int j) const {
        assert(i >= 0 && i < 2 * nc_ && j >= 0 && j < dim_);
        return v_[static_cast<std::size_t>(i) * static_cast<std::size_t>(dim_)
                  + static_cast<std::size_t>(j)];
    }
    // component d (0/1), coefficient rank r in index_map(k+1), basis column j
    const S& comp_at(int d, int r, int j) const {
        assert(d == 0 || d == 1);
        assert(r >= 0 && r < nc_);
        return at(d * nc_ + r, j);
    }

private:
    friend struct detail::rt_table_access;
    rt_basis_tbl() : k_(0), nc_(0), dim_(0), v_() {}
    int k_, nc_, dim_;
    std::vector<S> v_;
};

// T-R2: divergence coefficients, N_k x dim (degree k Bernstein rows)
template <typename S>
class rt_div_tbl {
public:
    int order() const { return k_; }
    int rows() const { return rows_; }          // N_k
    int cols() const { return cols_; }          // dim
    const S& at(int i, int j) const {
        assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
        return v_[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)
                  + static_cast<std::size_t>(j)];
    }

private:
    friend struct detail::rt_table_access;
    rt_div_tbl() : k_(0), rows_(0), cols_(0), v_() {}
    int k_, rows_, cols_;
    std::vector<S> v_;
};

// T-R3: edge normal flux, 3 edges x (k+1) 1D coefficients x dim
template <typename S>
class rt_flux_tbl {
public:
    int order() const { return k_; }
    int coeffs_per_edge() const { return k_ + 1; }
    int dim() const { return dim_; }
    // edge e, 1D Bernstein coefficient j (degree k), basis column c
    const S& at(int e, int j, int c) const {
        assert(e >= 0 && e < 3);
        assert(j >= 0 && j <= k_);
        assert(c >= 0 && c < dim_);
        return v_[(static_cast<std::size_t>(e) * static_cast<std::size_t>(k_ + 1)
                   + static_cast<std::size_t>(j)) * static_cast<std::size_t>(dim_)
                  + static_cast<std::size_t>(c)];
    }

private:
    friend struct detail::rt_table_access;
    rt_flux_tbl() : k_(0), dim_(0), v_() {}
    int k_, dim_;
    std::vector<S> v_;
};

// read-only view of one (transposable) block of a block table
template <typename S>
class rt_block_view {
public:
    rt_block_view(const S* base, int r, int c, bool trans)
        : base_(base), r_(r), c_(c), trans_(trans) {}
    int rows() const { return trans_ ? c_ : r_; }
    int cols() const { return trans_ ? r_ : c_; }
    const S& at(int i, int j) const {
        assert(i >= 0 && i < rows() && j >= 0 && j < cols());
        return trans_
            ? base_[static_cast<std::size_t>(j) * static_cast<std::size_t>(c_)
                    + static_cast<std::size_t>(i)]
            : base_[static_cast<std::size_t>(i) * static_cast<std::size_t>(c_)
                    + static_cast<std::size_t>(j)];
    }

private:
    const S* base_;
    int r_, c_;
    bool trans_;
};

// T-R4 (3 stored blocks (0,0), (0,1), (1,1); block(1,0) is the transposed
// view) and T-R5 (single block; the plain rows/cols/at surface).
template <typename S>
class rt_block_tbl {
public:
    int block_rows() const { return brows_; }
    int block_cols() const { return bcols_; }
    int num_blocks() const { return nblk_; }    // 1 (T-R5) or 3 (T-R4)

    // single-block surface (T-R5 usage; also block (0,0) of T-R4)
    int rows() const { return brows_; }
    int cols() const { return bcols_; }
    const S& at(int i, int j) const {
        assert(i >= 0 && i < brows_ && j >= 0 && j < bcols_);
        return blk_[0][static_cast<std::size_t>(i) * static_cast<std::size_t>(bcols_)
                      + static_cast<std::size_t>(j)];
    }

    // component pair surface (T-R4): symmetry R^{(d'd)} = R^{(dd')^T} is
    // realized as a transposed view of the stored (0,1) block
    rt_block_view<S> block(int d, int dp) const {
        assert(nblk_ == 3);
        assert((d == 0 || d == 1) && (dp == 0 || dp == 1));
        if (d == 0 && dp == 0) return rt_block_view<S>(blk_[0].data(), brows_, bcols_, false);
        if (d == 0 && dp == 1) return rt_block_view<S>(blk_[1].data(), brows_, bcols_, false);
        if (d == 1 && dp == 0) return rt_block_view<S>(blk_[1].data(), brows_, bcols_, true);
        return rt_block_view<S>(blk_[2].data(), brows_, bcols_, false);
    }

private:
    friend struct detail::rt_table_access;
    rt_block_tbl() : brows_(0), bcols_(0), nblk_(0), blk_() {}
    int brows_, bcols_, nblk_;
    std::vector<std::vector<S> > blk_;
};

// T-R6: cross table, 6 blocks (d = 0,1; i = 0,1,2), each dim x N_n
template <typename S>
class rt_cross_tbl {
public:
    int dim() const { return dim_; }            // rows (RT side)
    int nn() const { return nn_; }              // cols (N_n)
    rt_block_view<S> block(int d, int i) const {
        assert((d == 0 || d == 1) && i >= 0 && i < 3);
        return rt_block_view<S>(blk_[static_cast<std::size_t>(d * 3 + i)].data(),
                                dim_, nn_, false);
    }

private:
    friend struct detail::rt_table_access;
    rt_cross_tbl() : dim_(0), nn_(0), blk_() {}
    int dim_, nn_;
    std::vector<std::vector<S> > blk_;
};

// rational-stage table names of the external design section 3
typedef rt_mat_tbl<detail::rational>   rt_mat_table;
typedef rt_basis_tbl<detail::rational> rt_basis_table;
typedef rt_div_tbl<detail::rational>   rt_div_table;
typedef rt_flux_tbl<detail::rational>  rt_flux_table;
typedef rt_block_tbl<detail::rational> rt_block_table;
typedef rt_cross_tbl<detail::rational> rt_cross_table;

namespace detail {

// ---------------------------------------------------------------------------
// Reference geometry constants (D = 2, normative; RD-1 hardcode table).
// Edge e runs (e+1)%3 -> (e+2)%3; nu_e = R_{-90}(end - start),
// R_{-90}(x, y) = (y, -x). Reference vertices (0,0), (1,0), (0,1).
// ---------------------------------------------------------------------------
inline int rt_edge_vertex(int e, int end) {
    assert(e >= 0 && e < 3 && (end == 0 || end == 1));
    return (e + 1 + end) % 3;
}

inline int rt_edge_normal(int e, int d) {
    assert(e >= 0 && e < 3 && (d == 0 || d == 1));
    static const int nu[3][2] = { { 1, 1 }, { -1, 0 }, { 0, -1 } };
    return nu[e][d];
}

// reference gradients grad lambda_i (integer components)
inline int rt_grad_lambda(int i, int d) {
    assert(i >= 0 && i <= 2 && (d == 0 || d == 1));
    static const int gl[3][2] = { { -1, -1 }, { 1, 0 }, { 0, 1 } };
    return gl[i][d];
}

inline int rt_dim2(int k) { return (k + 1) * (k + 3); }

// trace slice (internal design 4.1): rank in index_map<2>(m) of the
// multi-index with alpha_e = 0, alpha_{(e+2)%3} = t, alpha_{(e+1)%3} = m - t.
// The 2D Bernstein B^m_alpha restricted to edge e equals the 1D B^m_t along
// the local positive direction (RD-2).
inline int rt_trace_index(int e, int m, int t) {
    assert(e >= 0 && e < 3 && t >= 0 && t <= m);
    multi_index<2> al;
    al.a[static_cast<std::size_t>(e)] = 0;
    al.a[static_cast<std::size_t>((e + 2) % 3)] = t;
    al.a[static_cast<std::size_t>((e + 1) % 3)] = m - t;
    return coeff_registry<2>::indices(m).rank(al);
}

// ---------------------------------------------------------------------------
// rt_span (internal design 3): spanning set S, 2 N_{k+1} x dim.
// Column order is the GENERATION order (P_k)^2 block then x P~_k block; the
// canonical DOF order is produced by C = S V^{-1} automatically.
// ---------------------------------------------------------------------------
inline rmat rt_span_matrix(int k) {
    assert(k >= 0);
    const int m = k + 1;
    const int n1 = coeff_registry<2>::indices(m).size();
    const int nk = coeff_registry<2>::indices(k).size();
    const int dim = rt_dim2(k);
    rmat S(2 * n1, dim);
    // (P_k)^2 block: e_d B^k_alpha elevated to degree k+1
    const elevation_table<2>& E = coeff_registry<2>::elevation(k, m);
    int col = 0;
    for (int d = 0; d < 2; ++d) {
        for (int ar = 0; ar < nk; ++ar) {
            elevation_table<2>::entry_range rr = E.row(ar);
            for (elevation_table<2>::entry_iterator p = rr.begin(); p != rr.end(); ++p) {
                elevation_table<2>::entry en = *p;
                S.at(d * n1 + en.target_rank, col) = *en.coeff;
            }
            ++col;
        }
    }
    // x P~_k block: x lambda_1^a lambda_2^b, a + b = k; closed form
    // lambda^gamma = B^{k+1}_gamma / multinomial(k+1, gamma): one nonzero
    // per component (RTX-3)
    const index_map<2>& im1 = coeff_registry<2>::indices(m);
    for (int a = k; a >= 0; --a) {
        int b = k - a;
        multi_index<2> g0, g1;
        g0.a[0] = 0; g0.a[1] = a + 1; g0.a[2] = b;
        g1.a[0] = 0; g1.a[1] = a;     g1.a[2] = b + 1;
        S.at(0 * n1 + im1.rank(g0), col) =
            rational(bigint(1), multinomial(m, g0.a.data(), 3));
        S.at(1 * n1 + im1.rank(g1), col) =
            rational(bigint(1), multinomial(m, g1.a.data(), 3));
        ++col;
    }
    assert(col == dim);
    return S;
}

// ---------------------------------------------------------------------------
// rt_dof (internal design 4): the DOF application matrix L (dim x 2 N_{k+1})
// acting on embedded degree k+1 component-major coefficients, and the
// Vandermonde V = L S. Row order is the canonical DOF order.
// ---------------------------------------------------------------------------
inline rmat rt_dof_apply_matrix(int k) {
    assert(k >= 0);
    const int m = k + 1;
    const int n1 = coeff_registry<2>::indices(m).size();
    const int dim = rt_dim2(k);
    rmat L(dim, 2 * n1);
    // edge moment rows: ell^e_j(sigma) = sum_d (nu_e)_d [M^{(k,m)}_{1D} c^e_d]_j
    const mass_table<1>& M1 = coeff_registry<1>::mass(k, m);
    for (int e = 0; e < 3; ++e) {
        for (int j = 0; j <= k; ++j) {
            int row = e * (k + 1) + j;
            for (int t = 0; t <= m; ++t) {
                int tr = rt_trace_index(e, m, t);
                for (int d = 0; d < 2; ++d) {
                    int nu = rt_edge_normal(e, d);
                    if (nu == 0) continue;
                    L.at(row, d * n1 + tr) += rational(nu) * M1.at(j, t);
                }
            }
        }
    }
    // interior moment rows (k >= 1):
    // ell^int_{d,alpha}(sigma) = int_T sigma_d B^{k-1}_alpha dx
    //                          = |T_hat| [M^{(k-1,m)} c^{(d)}]_{rank(alpha)}
    // (|T_hat| = 1/D! on top of the L0 "coefficient of |T|" convention)
    if (k >= 1) {
        const int nkm1 = coeff_registry<2>::indices(k - 1).size();
        const mass_table<2>& M2 = coeff_registry<2>::mass(k - 1, m);
        const rational that(1, 2);              // |T_hat|, D = 2
        for (int d = 0; d < 2; ++d) {
            for (int ar = 0; ar < nkm1; ++ar) {
                int row = 3 * (k + 1) + d * nkm1 + ar;
                for (int b = 0; b < n1; ++b) {
                    const rational& v = M2.at(ar, b);
                    if (v.is_zero()) continue;
                    L.at(row, d * n1 + b) = that * v;
                }
            }
        }
    }
    return L;
}

inline rmat rt_dof_matrix(int k) {
    return mul(rt_dof_apply_matrix(k), rt_span_matrix(k));
}

} // namespace detail

// ---------------------------------------------------------------------------
// rt_registry<D>: rational stage of the two-stage cache (external design
// section 3). Lazy generation under one mutex per D; returned references
// stay valid until program termination (L0 section 4.1 contract inherited).
// ---------------------------------------------------------------------------
template <int D>
class rt_registry {
    static_assert(D == 2, "bfem::rt_registry: initial version supports D == 2 only");
public:
    // dim RT_k = (k+1)(k+3)
    static int dim(int k) {
        check_k(k, "dim");
        return detail::rt_dim2(k);
    }

    // T-R1
    static const rt_basis_table& basis(int k) {
        check_k(k, "basis");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return basis_locked(s, k);
    }

    // T-R2
    static const rt_div_table& divergence(int k) {
        check_k(k, "divergence");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return div_locked(s, k);
    }

    // T-R3
    static const rt_flux_table& edge_flux(int k) {
        check_k(k, "edge_flux");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return flux_locked(s, k);
    }

    // T-R4
    static const rt_block_table& comp_mass(int k) {
        check_k(k, "comp_mass");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return cmass_locked(s, k);
    }

    // T-R5
    static const rt_block_table& div_mass(int k, int l) {
        check_k(k, "div_mass");
        if (l < 0)
            throw std::invalid_argument("bfem::rt_registry::div_mass: l < 0");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, l);
        typename std::map<std::pair<int, int>, rt_block_table>::iterator it =
            s.dmass.find(key);
        if (it != s.dmass.end()) return it->second;
        rt_block_table t = build_div_mass(s, k, l);
        return s.dmass.insert(std::make_pair(key, std::move(t))).first->second;
    }

    // T-R6 (n >= 1)
    static const rt_cross_table& cross_grad(int k, int n) {
        check_k(k, "cross_grad");
        check_n(n, "cross_grad");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return cross_locked(s, k, n);
    }

    // T-R6' (n >= 1)
    static const rt_mat_table& cross_grad_contracted(int k, int n) {
        check_k(k, "cross_grad_contracted");
        check_n(n, "cross_grad_contracted");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(k, n);
        typename std::map<std::pair<int, int>, rt_mat_table>::iterator it =
            s.crossc.find(key);
        if (it != s.crossc.end()) return it->second;
        rt_mat_table t = build_cross_contracted(s, k, n);
        return s.crossc.insert(std::make_pair(key, std::move(t))).first->second;
    }

    // T-R7 (l >= 0; v0.3)
    static const rt_mat_table& inv_mass(int l) {
        if (l < 0)
            throw std::invalid_argument("bfem::rt_registry::inv_mass: l < 0");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, rt_mat_table>::iterator it = s.invm.find(l);
        if (it != s.invm.end()) return it->second;
        rt_mat_table t = build_inv_mass(l);
        return s.invm.insert(std::make_pair(l, std::move(t))).first->second;
    }

private:
    struct state {
        std::mutex mtx;
        std::map<int, rt_basis_table> basis;
        std::map<int, rt_div_table> div;
        std::map<int, rt_flux_table> flux;
        std::map<int, rt_block_table> cmass;
        std::map<std::pair<int, int>, rt_block_table> dmass;
        std::map<std::pair<int, int>, rt_cross_table> cross;
        std::map<std::pair<int, int>, rt_mat_table> crossc;
        std::map<int, rt_mat_table> invm;
    };
    static state& st() {
        static state s;
        return s;
    }
    static void check_k(int k, const char* where) {
        if (k < 0) {
            std::string msg("bfem::rt_registry::");
            msg += where;
            msg += ": k < 0";
            throw std::invalid_argument(msg);
        }
    }
    static void check_n(int n, const char* where) {
        if (n < 1) {
            std::string msg("bfem::rt_registry::");
            msg += where;
            msg += ": n < 1";
            throw std::invalid_argument(msg);
        }
    }

    // ---- locked generation (callers hold s.mtx) ----

    static const rt_basis_table& basis_locked(state& s, int k) {
        typename std::map<int, rt_basis_table>::iterator it = s.basis.find(k);
        if (it != s.basis.end()) return it->second;
        rt_basis_table t = build_basis(k);
        return s.basis.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const rt_div_table& div_locked(state& s, int k) {
        typename std::map<int, rt_div_table>::iterator it = s.div.find(k);
        if (it != s.div.end()) return it->second;
        rt_div_table t = build_div(basis_locked(s, k));
        return s.div.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const rt_flux_table& flux_locked(state& s, int k) {
        typename std::map<int, rt_flux_table>::iterator it = s.flux.find(k);
        if (it != s.flux.end()) return it->second;
        rt_flux_table t = build_flux(basis_locked(s, k));
        return s.flux.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const rt_block_table& cmass_locked(state& s, int k) {
        typename std::map<int, rt_block_table>::iterator it = s.cmass.find(k);
        if (it != s.cmass.end()) return it->second;
        rt_block_table t = build_comp_mass(basis_locked(s, k));
        return s.cmass.insert(std::make_pair(k, std::move(t))).first->second;
    }
    static const rt_cross_table& cross_locked(state& s, int k, int n) {
        std::pair<int, int> key(k, n);
        typename std::map<std::pair<int, int>, rt_cross_table>::iterator it =
            s.cross.find(key);
        if (it != s.cross.end()) return it->second;
        rt_cross_table t = build_cross(basis_locked(s, k), n);
        return s.cross.insert(std::make_pair(key, std::move(t))).first->second;
    }

    // ---- builders (rational stage; all contractions happen HERE, external
    //      design 1.3-3: never at the typed stage) ----

    static rt_basis_table build_basis(int k);
    static rt_div_table build_div(const rt_basis_table& C);
    static rt_flux_table build_flux(const rt_basis_table& C);
    static rt_block_table build_comp_mass(const rt_basis_table& C);
    static rt_block_table build_div_mass(state& s, int k, int l);
    static rt_cross_table build_cross(const rt_basis_table& C, int n);
    static rt_mat_table build_cross_contracted(state& s, int k, int n);
    static rt_mat_table build_inv_mass(int l);
};

namespace detail {

// the single fill point of the immutable tables (rational and typed stages)
struct rt_table_access {
    template <typename S>
    static rt_mat_tbl<S> make_mat(int rows, int cols, std::vector<S> v) {
        rt_mat_tbl<S> t;
        t.rows_ = rows;
        t.cols_ = cols;
        t.v_ = std::move(v);
        return t;
    }
    template <typename S>
    static rt_basis_tbl<S> make_basis(int k, int nc, int dim, std::vector<S> v) {
        rt_basis_tbl<S> t;
        t.k_ = k;
        t.nc_ = nc;
        t.dim_ = dim;
        t.v_ = std::move(v);
        return t;
    }
    template <typename S>
    static rt_div_tbl<S> make_div(int k, int rows, int cols, std::vector<S> v) {
        rt_div_tbl<S> t;
        t.k_ = k;
        t.rows_ = rows;
        t.cols_ = cols;
        t.v_ = std::move(v);
        return t;
    }
    template <typename S>
    static rt_flux_tbl<S> make_flux(int k, int dim, std::vector<S> v) {
        rt_flux_tbl<S> t;
        t.k_ = k;
        t.dim_ = dim;
        t.v_ = std::move(v);
        return t;
    }
    template <typename S>
    static rt_block_tbl<S> make_block(int brows, int bcols,
                                      std::vector<std::vector<S> > blk) {
        rt_block_tbl<S> t;
        t.brows_ = brows;
        t.bcols_ = bcols;
        t.nblk_ = static_cast<int>(blk.size());
        t.blk_ = std::move(blk);
        return t;
    }
    template <typename S>
    static rt_cross_tbl<S> make_cross(int dim, int nn,
                                      std::vector<std::vector<S> > blk) {
        rt_cross_tbl<S> t;
        t.dim_ = dim;
        t.nn_ = nn;
        t.blk_ = std::move(blk);
        return t;
    }
};

// rmat -> flat row-major rational vector
inline std::vector<rational> rmat_flat(const rmat& A) {
    return A.a;
}

// basis table block (component d) as an rmat (N_{k+1} x dim)
inline rmat basis_comp_rmat(const rt_basis_tbl<rational>& C, int d) {
    const int n1 = C.comp_size();
    const int dim = C.dim();
    rmat B(n1, dim);
    for (int r = 0; r < n1; ++r)
        for (int j = 0; j < dim; ++j)
            B.at(r, j) = C.comp_at(d, r, j);
    return B;
}

// L0 mass table as rmat
template <int D>
inline rmat mass_rmat(int a, int b) {
    const mass_table<D>& M = coeff_registry<D>::mass(a, b);
    rmat R(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); ++i)
        for (int j = 0; j < M.cols(); ++j)
            R.at(i, j) = M.at(i, j);
    return R;
}

} // namespace detail

// ---------------------------------------------------------------------------
// builder definitions
// ---------------------------------------------------------------------------

template <int D>
rt_basis_table rt_registry<D>::build_basis(int k) {
    using detail::rmat;
    const int m = k + 1;
    const int n1 = coeff_registry<2>::indices(m).size();
    const int dim = detail::rt_dim2(k);
    rmat S = detail::rt_span_matrix(k);
    rmat V = detail::rt_dof_matrix(k);
    // exact inversion; singular V would throw std::logic_error (RD-5 wiring)
    rmat Vinv = detail::solve_exact(std::move(V), rmat::identity(dim));
    rmat C = detail::mul(S, Vinv);
    return detail::rt_table_access::make_basis(k, n1, dim, detail::rmat_flat(C));
}

template <int D>
rt_div_table rt_registry<D>::build_div(const rt_basis_table& C) {
    const int k = C.order();
    const int m = k + 1;
    const int n1 = C.comp_size();
    const int nk = coeff_registry<2>::indices(k).size();
    const int dim = C.dim();
    derivative_map<2> dm(m);
    std::vector<detail::rational> v(static_cast<std::size_t>(nk)
                                    * static_cast<std::size_t>(dim));
    for (int j = 0; j < dim; ++j) {
        for (int d = 0; d < 2; ++d) {
            for (int b = 0; b < n1; ++b) {
                const detail::rational& c = C.comp_at(d, b, j);
                if (c.is_zero()) continue;
                for (int i = 0; i <= 2; ++i) {
                    int g = detail::rt_grad_lambda(i, d);
                    if (g == 0) continue;
                    int tr = dm.target(b, i);
                    if (tr < 0) continue;
                    v[static_cast<std::size_t>(tr) * static_cast<std::size_t>(dim)
                      + static_cast<std::size_t>(j)] +=
                        detail::rational(g * m) * c;
                }
            }
        }
    }
    return detail::rt_table_access::make_div(k, nk, dim, std::move(v));
}

template <int D>
rt_flux_table rt_registry<D>::build_flux(const rt_basis_table& C) {
    using detail::rmat;
    const int k = C.order();
    const int m = k + 1;
    const int dim = C.dim();
    // 1D elevation matrix E^{k -> k+1} as a tall (m+1) x (k+1) rmat
    const elevation_table<1>& E = coeff_registry<1>::elevation(k, m);
    rmat E1(m + 1, k + 1);
    for (int src = 0; src <= k; ++src) {
        elevation_table<1>::entry_range rr = E.row(src);
        for (elevation_table<1>::entry_iterator p = rr.begin(); p != rr.end(); ++p) {
            elevation_table<1>::entry en = *p;
            E1.at(en.target_rank, src) = *en.coeff;
        }
    }
    std::vector<detail::rational> v(static_cast<std::size_t>(3)
                                    * static_cast<std::size_t>(k + 1)
                                    * static_cast<std::size_t>(dim));
    for (int e = 0; e < 3; ++e) {
        // degree k+1 1D flux coefficients of all basis columns on edge e
        rmat F(m + 1, dim);
        for (int t = 0; t <= m; ++t) {
            int tr = detail::rt_trace_index(e, m, t);
            for (int j = 0; j < dim; ++j) {
                detail::rational acc;
                for (int d = 0; d < 2; ++d) {
                    int nu = detail::rt_edge_normal(e, d);
                    if (nu == 0) continue;
                    acc += detail::rational(nu) * C.comp_at(d, tr, j);
                }
                F.at(t, j) = acc;
            }
        }
        // exact degree reduction k+1 -> k; a nonzero residual raises
        // std::logic_error inside solve_consistent (RA-4 generation-time
        // version, internal design 5.2)
        rmat G = detail::solve_consistent(E1, std::move(F));
        for (int j2 = 0; j2 <= k; ++j2)
            for (int c = 0; c < dim; ++c)
                v[(static_cast<std::size_t>(e) * static_cast<std::size_t>(k + 1)
                   + static_cast<std::size_t>(j2)) * static_cast<std::size_t>(dim)
                  + static_cast<std::size_t>(c)] = G.at(j2, c);
    }
    return detail::rt_table_access::make_flux(k, dim, std::move(v));
}

template <int D>
rt_block_table rt_registry<D>::build_comp_mass(const rt_basis_table& C) {
    using detail::rmat;
    const int k = C.order();
    const int dim = C.dim();
    rmat M = detail::mass_rmat<2>(k + 1, k + 1);
    rmat C0 = detail::basis_comp_rmat(C, 0);
    rmat C1 = detail::basis_comp_rmat(C, 1);
    rmat MC0 = detail::mul(M, C0);
    rmat MC1 = detail::mul(M, C1);
    rmat C0t = detail::transpose(C0);
    rmat C1t = detail::transpose(C1);
    std::vector<std::vector<detail::rational> > blk;
    blk.reserve(3);
    blk.push_back(detail::rmat_flat(detail::mul(C0t, MC0)));   // (0,0)
    blk.push_back(detail::rmat_flat(detail::mul(C0t, MC1)));   // (0,1)
    blk.push_back(detail::rmat_flat(detail::mul(C1t, MC1)));   // (1,1)
    return detail::rt_table_access::make_block(dim, dim, std::move(blk));
}

template <int D>
rt_block_table rt_registry<D>::build_div_mass(state& s, int k, int l) {
    using detail::rmat;
    const rt_div_table& Dv = div_locked(s, k);
    const int nk = Dv.rows();
    const int dim = Dv.cols();
    rmat Dm(nk, dim);
    for (int i = 0; i < nk; ++i)
        for (int j = 0; j < dim; ++j)
            Dm.at(i, j) = Dv.at(i, j);
    rmat M = detail::mass_rmat<2>(l, k);
    rmat R = detail::mul(M, Dm);                              // N_l x dim
    std::vector<std::vector<detail::rational> > blk;
    blk.push_back(detail::rmat_flat(R));
    return detail::rt_table_access::make_block(R.rows, R.cols, std::move(blk));
}

template <int D>
rt_cross_table rt_registry<D>::build_cross(const rt_basis_table& C, int n) {
    using detail::rmat;
    const int k = C.order();
    const int dim = C.dim();
    const int nn = coeff_registry<2>::indices(n).size();
    const int n1 = C.comp_size();
    derivative_map<2> dmn(n);
    rmat M = detail::mass_rmat<2>(k + 1, n - 1);              // N_{k+1} x N_{n-1}
    const detail::rational fn(n);
    std::vector<std::vector<detail::rational> > blk;
    blk.reserve(6);
    for (int d = 0; d < 2; ++d) {
        rmat Cdt = detail::transpose(detail::basis_comp_rmat(C, d));  // dim x N_{k+1}
        for (int i = 0; i <= 2; ++i) {
            // W(b, ar) = n * M(b, dmap_n(ar, i)) (0 on vanishing)
            rmat W(n1, nn);
            for (int ar = 0; ar < nn; ++ar) {
                int t = dmn.target(ar, i);
                if (t < 0) continue;
                for (int b = 0; b < n1; ++b) {
                    const detail::rational& mv = M.at(b, t);
                    if (mv.is_zero()) continue;
                    W.at(b, ar) = fn * mv;
                }
            }
            blk.push_back(detail::rmat_flat(detail::mul(Cdt, W)));
        }
    }
    return detail::rt_table_access::make_cross(dim, nn, std::move(blk));
}

template <int D>
rt_mat_table rt_registry<D>::build_cross_contracted(state& s, int k, int n) {
    const rt_cross_table& X = cross_locked(s, k, n);
    const int dim = X.dim();
    const int nn = X.nn();
    std::vector<detail::rational> v(static_cast<std::size_t>(dim)
                                    * static_cast<std::size_t>(nn));
    // X' = sum_{d,i} (grad lambda_i)_d X^{(d,i)}: signed rational additions
    // only (the gradient components are -1/0/+1)
    for (int d = 0; d < 2; ++d) {
        for (int i = 0; i <= 2; ++i) {
            int g = detail::rt_grad_lambda(i, d);
            if (g == 0) continue;
            rt_block_view<detail::rational> B = X.block(d, i);
            for (int j = 0; j < dim; ++j) {
                for (int a = 0; a < nn; ++a) {
                    const detail::rational& x = B.at(j, a);
                    if (x.is_zero()) continue;
                    std::size_t idx = static_cast<std::size_t>(j)
                                      * static_cast<std::size_t>(nn)
                                      + static_cast<std::size_t>(a);
                    if (g > 0) v[idx] += x;
                    else       v[idx] -= x;
                }
            }
        }
    }
    return detail::rt_table_access::make_mat(dim, nn, std::move(v));
}

template <int D>
rt_mat_table rt_registry<D>::build_inv_mass(int l) {
    using detail::rmat;
    rmat M = detail::mass_rmat<2>(l, l);
    const int n = M.rows;
    rmat Minv = detail::solve_exact(M, rmat::identity(n));
    // generation-time exactness check (RA-9 generation version)
    rmat P = detail::mul(detail::mass_rmat<2>(l, l), Minv);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (!(P.at(i, j) == detail::rational(i == j ? 1 : 0)))
                throw std::logic_error(
                    "bfem::rt_registry::inv_mass: exact inverse check failed");
    return detail::rt_table_access::make_mat(n, n, detail::rmat_flat(Minv));
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_TABLES_HPP
