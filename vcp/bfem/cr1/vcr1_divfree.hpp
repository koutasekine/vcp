// vcp/bfem/cr1/vcr1_divfree.hpp
// CONST-G (P1): the element-mean divergence constraint of the vector CR space
// and its EXACT null space -- the discrete divergence-free subspace
// V_{0,h}^(1) of Xie-Xie-Liu 2018, equation (9), as an explicit basis Z over
// the exact rational scalar, together with the certificate that Z IS a basis
// of the kernel, and the reduced pencil Z^T A Z / Z^T M Z in the interval
// scalar of the caller.
//
// Conforms to: CONST-G design v1.1 (sections 1-2, 3, 4 G-S5) and the P0
//              ruling of 2026-09-07 (N-1 = option (i) exact null space;
//              cr1 layer; NO constant here).
//
// ---------------------------------------------------------------------------
// (i) Mathematics (verbatim-verified against the primary source, CONST-G P0)
// ---------------------------------------------------------------------------
// Constraint space, Xie-Xie-Liu equation (9):
//
//     V_{0,h}^(1) = { v_h in V_h^(1) : b_h(v_h, q_h) = 0  for all q_h in Q_h },
//     b_h(v_h, q_h) = sum_K int_K q_h div v_h,   Q_h = piecewise P0, mean zero.
//
// Reading as "element-mean divergence zero" (design section 1, item 2):
// Q_h is the MEAN-ZERO piecewise constants, so (9) literally says that the
// vector (int_K div v_h)_K is orthogonal to every mean-zero piecewise
// constant, i.e. div v_h has one and the same element mean c on every
// element.  The discrete divergence theorem then forces c = 0: summing over
// the elements, sum_K int_K div v_h = sum_K int_{dK} v_h . n; the interior
// facet contributions cancel pairwise (a CR function has the same facet mean
// seen from both sides and the normals are opposite) and the boundary facet
// contributions vanish because the boundary facet means are zero.  Hence
// c |Omega| = 0.  So, on the Dirichlet-reduced space,
//
//     (9)  <=>  int_K div v_h = 0  for every element K,
//
// and the constraint matrix B of this header has one row per element.  The
// same identity says that the rows of B sum to zero on the reduced space,
// which is the ONE linear dependency among them (rank B = nt - 1 on a
// connected domain, an observation the certificate below reports as
// rank_deficiency and the CONST-G gates check; the certificate itself does
// not depend on it).  [self-contained; note (f) of the P0 ruling]
//
// Fortin property, Xie-Xie-Liu Proposition 1 item 1 (with proof, p. 340-341):
// the CR interpolation Pi_h^(1) maps V_0 (the continuous divergence-free
// space) INTO V_{0,h}^(1).  The proof uses only that grad(q_h|_K) = 0 and
// that Pi_h preserves facet means:
//
//     b_h(Pi_h u, q_h) = sum_K int_{dK} q_h n . Pi_h u = sum_K int_{dK} q_h n . u
//                      = b_h(u, q_h) = 0.
//
// That inclusion is what makes the lower bound theorem (Theorem 4 of the
// paper, consumed in stokes_constants.hpp through eigenvalue_lower_bound)
// applicable to the eigenvalues of the reduced pencil built here.
//
// Element integral of the basis divergence (the entries of B):
//
//     int_K d/dx_d phi_i = |K| (-D) d/dx_d lambda_i,
//
// because phi_i = 1 - D lambda_i on K (cr1_element_op.hpp) and grad lambda_i
// is constant on K.  |K| and grad lambda come from element_geometry<D, R>
// on the EXACT rational scalar R, so every entry of B is exact.  (The
// equivalent division-free form |f_i| n_{i,d} -- facet measure times outward
// normal -- is used by the CONST-G gates as an independent cross-check of B;
// the library keeps the geometry form because it is D-generic and already
// verified.)
//
// ---------------------------------------------------------------------------
// (ii) The null space and its certificate (design G-S5)
// ---------------------------------------------------------------------------
// B is reduced to REDUCED ROW ECHELON FORM by exact rational elimination on
// sparse rows (pivot = the last column of the reduced row; every pivot
// column is eliminated from all other pivot rows, so the pivot rows are
// mutually reduced).  The exact rank is the number of pivot rows.  For every
// FREE (non-pivot) column f the vector
//
//     z_f = e_f - sum_{pivot rows p}  R_p[f] e_{pivot(p)}
//
// lies in ker B, and the z_f are linearly independent because z_f is the
// only one with a nonzero entry in column f.  Rank-nullity gives
// dim ker B = ncols - rank, which is exactly the number of free columns, so
// {z_f} is a basis of ker B.  Each z_f is then scaled by the least common
// multiple of its denominators, which keeps it in the kernel and makes every
// entry an INTEGER (on the uniform meshes of the CONST-G anchors every entry
// is 0 or +-1, measured in P0), so that its conversion to an interval scalar
// through convert_traits<T>::from_rational is exact.
//
// The certificate is computed from the FINISHED Z, independently of how it
// was built:
//   (a) B z_f == 0 exactly, for EVERY row of B and EVERY column of Z (all
//       nt * dim products; count reported);
//   (b) dim == ncols - rank (the rank-nullity identity, with the exact rank);
//   (c) every column carries the unit entry in its own free column and every
//       other column is zero there (the independence witness).
// A failure of (a), (b) or (c) throws vcp::verification_error from the
// constructor: an uncertified basis is never handed out (the certified-only
// policy of the constants layer, applied at the source).
//
// ---------------------------------------------------------------------------
// (iii) The reduced pencil
// ---------------------------------------------------------------------------
// reduce_pencil(K) takes the scalar Dirichlet-REDUCED matrix K (n_r x n_r,
// the stiffness or the mass of cr1_space after dirichlet_reduction) and
// returns Z^T blockdiag(K, ..., K) Z (dim x dim, dense, in the interval
// scalar T of K).  Z is converted entry by entry through
// convert_traits<T>::from_rational (enclose-once); W = blockdiag(K) Z is
// accumulated through the sparsity of Z, then Z^T W is formed for j >= i
// and MIRRORED, so that the returned matrix is bit-for-bit symmetric (the
// symmetry precondition of eigsymge is checked on the midpoint matrix) while
// every entry still encloses the exact value.  The reduced vector dof of
// component d and reduced scalar dof r is d * n_r + r (vcr1_space.hpp,
// section (ii)); B is assembled directly in that numbering.
//
// Lexical policy: no decimal literal, no `double` / `float` token in the
// code.  No OpenMP pragma.  R must be an exact field type (the canonical one
// is vcp::bfem::detail::rational); T must be a scalar with inclusion
// arithmetic or a point type -- the reduced pencil claims exactness only for
// the Z side, which is integer.

#ifndef VCP_BFEM_CR1_VCR1_DIVFREE_HPP
#define VCP_BFEM_CR1_VCR1_DIVFREE_HPP

#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <cstddef>
#include <stdexcept>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/fe_space.hpp>        // detail::spm_adapter
#include <vcp/bfem/cr1/cr1_dofmap.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// vcr1_divfree_certificate: what the constructor verified.
// ---------------------------------------------------------------------------
struct vcr1_divfree_certificate {
    long long pairs_checked;      // rows of B times columns of Z (all of them)
    long long violations;         // pairs with B z != 0 (must be 0)
    int rows;                     // nt (elements)
    int columns;                  // D * n_r
    int rank;                     // exact rank of B
    int rank_deficiency;          // rows - rank (1 on a connected domain)
    int dim;                      // columns of Z
    bool dim_identity;            // dim == columns - rank
    bool independence_witness;    // unit entries in the free columns
    bool ok() const {
        return violations == 0 && dim_identity && independence_witness;
    }
};

// ---------------------------------------------------------------------------
// vcr1_divfree_basis<D, R>
// ---------------------------------------------------------------------------
template <int D, typename R = detail::rational>
class vcr1_divfree_basis {
    static_assert(D == 2 || D == 3,
                  "bfem::vcr1_divfree_basis: only D == 2 or D == 3");
public:
    typedef std::vector<std::pair<int, R> > sparse_vec;   // (index, value), ascending
    typedef detail::cr1_dofmap<D> dofmap_type;

    // exact_mesh: the mesh on the exact scalar R.
    // scalar_to_reduced: size = number of facets; the reduced scalar dof of
    // facet i, or -1 for a Dirichlet facet (dirichlet_reduction::to_reduced
    // of the scalar space, tabulated).  The unconstrained facets must map
    // bijectively onto 0 .. n_r - 1.
    vcr1_divfree_basis(const mesh<D, R>& exact_mesh,
                       const std::vector<int>& scalar_to_reduced)
        : nt_(0), n_reduced_(0), ncols_(0), rank_(0), rows_(), cols_(),
          free_col_(), cert_() {
        const dofmap_type dm = dofmap_type::from_mesh(exact_mesh);
        nt_ = dm.num_elements();
        if (static_cast<int>(scalar_to_reduced.size()) != dm.ndof())
            throw std::invalid_argument(
                "bfem::vcr1_divfree_basis: scalar_to_reduced.size() != number of facets");
        n_reduced_ = check_numbering_(scalar_to_reduced);
        ncols_ = D * n_reduced_;

        build_constraint_(exact_mesh, dm, scalar_to_reduced);
        rref_result rr = rref_();
        rank_ = static_cast<int>(rr.rows.size());
        build_nullspace_(rr);
        certify_();
        if (!cert_.ok())
            vcp::throw_error<vcp::verification_error>(
                "bfem::vcr1_divfree_basis: null space certificate failed "
                "(violations=", cert_.violations,
                ", dim_identity=", cert_.dim_identity ? 1 : 0,
                ", independence=", cert_.independence_witness ? 1 : 0, ")");
    }

    // ---- observers ----
    int num_elements() const { return nt_; }
    int num_reduced_scalar() const { return n_reduced_; }
    int num_columns() const { return ncols_; }          // D * n_r
    int rank() const { return rank_; }
    int dim() const { return static_cast<int>(cols_.size()); }
    const vcr1_divfree_certificate& certificate() const { return cert_; }

    // row e of B: int_K div phi over the reduced vector dofs (exact)
    const sparse_vec& constraint_row(int e) const {
        if (e < 0 || e >= nt_)
            throw std::invalid_argument("bfem::vcr1_divfree_basis::constraint_row: out of range");
        return rows_[static_cast<std::size_t>(e)];
    }
    // column c of Z (integer entries, exact)
    const sparse_vec& basis_column(int c) const {
        if (c < 0 || c >= dim())
            throw std::invalid_argument("bfem::vcr1_divfree_basis::basis_column: out of range");
        return cols_[static_cast<std::size_t>(c)];
    }
    // the free column that carries the unit entry of column c
    int free_column(int c) const {
        if (c < 0 || c >= dim())
            throw std::invalid_argument("bfem::vcr1_divfree_basis::free_column: out of range");
        return free_col_[static_cast<std::size_t>(c)];
    }

    // ---- the reduced pencil Z^T blockdiag(K) Z (section (iii)) ----
    template <typename T, class DP, class SP>
    vcp::matrix<T, DP> reduce_pencil(const vcp::spmatrix<T, SP>& K) const {
        if (static_cast<int>(K.rowsize()) != n_reduced_ ||
            static_cast<int>(K.columnsize()) != n_reduced_)
            throw std::invalid_argument(
                "bfem::vcr1_divfree_basis::reduce_pencil: K must be n_r x n_r "
                "(the scalar Dirichlet-reduced matrix)");
        const int N = dim();
        // row lists of K
        std::vector<std::vector<std::pair<int, T> > > krow(
            static_cast<std::size_t>(n_reduced_));
        row_gather<T> g = { &krow };
        detail::spm_adapter<T, SP>::for_each_entry(K, g);

        // Z in T (enclose-once; exact for integer entries)
        std::vector<std::vector<std::pair<int, T> > > zt(static_cast<std::size_t>(N));
        for (int c = 0; c < N; ++c) {
            const sparse_vec& z = cols_[static_cast<std::size_t>(c)];
            std::vector<std::pair<int, T> >& o = zt[static_cast<std::size_t>(c)];
            o.reserve(z.size());
            for (std::size_t k = 0; k < z.size(); ++k)
                o.push_back(std::make_pair(z[k].first, entry_to<T>(z[k].second)));
        }

        // W = blockdiag(K) Z, dense (D n_r) x N, accumulated over nnz(Z)
        vcp::matrix<T, DP> W;
        W.zeros(ncols_, N);
        for (int c = 0; c < N; ++c) {
            const std::vector<std::pair<int, T> >& z = zt[static_cast<std::size_t>(c)];
            for (std::size_t k = 0; k < z.size(); ++k) {
                const int col = z[k].first;
                const int d = col / n_reduced_;
                const int j = col - d * n_reduced_;
                const T& zv = z[k].second;
                const std::vector<std::pair<int, T> >& kr = krow[static_cast<std::size_t>(j)];
                // K symmetric: column j of K is row j of K
                for (std::size_t q = 0; q < kr.size(); ++q)
                    W(d * n_reduced_ + kr[q].first, c) += kr[q].second * zv;
            }
        }
        // Z^T W, upper triangle, mirrored
        vcp::matrix<T, DP> out;
        out.zeros(N, N);
        for (int i = 0; i < N; ++i) {
            const std::vector<std::pair<int, T> >& z = zt[static_cast<std::size_t>(i)];
            for (int j = i; j < N; ++j) {
                T s(0);
                for (std::size_t k = 0; k < z.size(); ++k)
                    s += z[k].second * W(z[k].first, j);
                out(i, j) = s;
                if (j != i) out(j, i) = s;
            }
        }
        return out;
    }

private:
    // ---- numbering check: bijection of the unconstrained facets ----
    static int check_numbering_(const std::vector<int>& s2r) {
        int nr = 0;
        for (std::size_t i = 0; i < s2r.size(); ++i)
            if (s2r[i] >= 0) ++nr;
        std::vector<int> seen(static_cast<std::size_t>(nr), 0);
        for (std::size_t i = 0; i < s2r.size(); ++i) {
            const int r = s2r[i];
            if (r < -1 || r >= nr)
                throw std::invalid_argument(
                    "bfem::vcr1_divfree_basis: scalar_to_reduced entry out of range");
            if (r >= 0) {
                if (seen[static_cast<std::size_t>(r)])
                    throw std::invalid_argument(
                        "bfem::vcr1_divfree_basis: scalar_to_reduced is not injective");
                seen[static_cast<std::size_t>(r)] = 1;
            }
        }
        return nr;
    }

    template <typename T>
    static T entry_to(const R& q) { return convert_traits<T>::from_rational(q.num(), q.den()); }

    template <typename T>
    struct row_gather {
        std::vector<std::vector<std::pair<int, T> > >* rows;
        template <typename I>
        void operator()(I i, I j, const T& v) const {
            (*rows)[static_cast<std::size_t>(i)].push_back(
                std::make_pair(static_cast<int>(j), v));
        }
    };

    struct by_index {
        bool operator()(const std::pair<int, R>& a, const std::pair<int, R>& b) const {
            return a.first < b.first;
        }
    };

    // ---- B: one row per element, int_K d/dx_d phi_i = |K| (-D) grad lambda_i[d]
    void build_constraint_(const mesh<D, R>& msh, const dofmap_type& dm,
                           const std::vector<int>& s2r) {
        rows_.assign(static_cast<std::size_t>(nt_), sparse_vec());
        const R minus_d = R(-D);
        for (int e = 0; e < nt_; ++e) {
            std::array<std::array<R, D>, D + 1> vv;
            for (int k = 0; k <= D; ++k)
                vv[static_cast<std::size_t>(k)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(k)]);
            const element_geometry<D, R> g = element_geometry<D, R>::from_vertices(vv);
            const R w = g.measure() * minus_d;
            sparse_vec row;
            for (int i = 0; i <= D; ++i) {
                const int r = s2r[static_cast<std::size_t>(dm.l2g(e, i))];
                if (r < 0) continue;
                for (int d = 0; d < D; ++d) {
                    const R v = w * g.grad_lambda(i, d);
                    if (!v.is_zero()) row.push_back(std::make_pair(d * n_reduced_ + r, v));
                }
            }
            std::sort(row.begin(), row.end(), by_index());
            rows_[static_cast<std::size_t>(e)].swap(row);
        }
    }

    // y += a * x on sorted sparse rows, exact zeros dropped
    static void axpy_(sparse_vec& y, const R& a, const sparse_vec& x) {
        sparse_vec out;
        out.reserve(y.size() + x.size());
        std::size_t i = 0, j = 0;
        while (i < y.size() || j < x.size()) {
            if (j >= x.size() || (i < y.size() && y[i].first < x[j].first)) {
                out.push_back(y[i]); ++i;
            } else if (i >= y.size() || x[j].first < y[i].first) {
                out.push_back(std::make_pair(x[j].first, a * x[j].second)); ++j;
            } else {
                const R s = y[i].second + a * x[j].second;
                if (!s.is_zero()) out.push_back(std::make_pair(y[i].first, s));
                ++i; ++j;
            }
        }
        y.swap(out);
    }

    struct rref_result {
        std::vector<sparse_vec> rows;     // pivot rows, mutually reduced, pivot entry 1
        std::vector<int> pivot;           // pivot column of each row
        std::vector<int> pivot_row_of;    // column -> row index, or -1
    };

    rref_result rref_() const {
        rref_result rr;
        rr.pivot_row_of.assign(static_cast<std::size_t>(ncols_), -1);
        for (int e = 0; e < nt_; ++e) {
            sparse_vec row = rows_[static_cast<std::size_t>(e)];
            for (;;) {                    // eliminate every existing pivot column
                bool changed = false;
                for (std::size_t k = 0; k < row.size(); ++k) {
                    const int pr = rr.pivot_row_of[static_cast<std::size_t>(row[k].first)];
                    if (pr < 0) continue;
                    const R a = -row[k].second;
                    axpy_(row, a, rr.rows[static_cast<std::size_t>(pr)]);
                    changed = true;
                    break;
                }
                if (!changed) break;
            }
            if (row.empty()) continue;    // dependent row
            const std::size_t pk = row.size() - 1;      // pivot: the last column
            const int pc = row[pk].first;
            const R inv = R(1) / row[pk].second;
            for (std::size_t k = 0; k < row.size(); ++k) row[k].second = row[k].second * inv;
            for (std::size_t q = 0; q < rr.rows.size(); ++q) {   // keep pivot rows mutually reduced
                sparse_vec& pr = rr.rows[q];
                for (std::size_t k = 0; k < pr.size(); ++k)
                    if (pr[k].first == pc) {
                        const R a = -pr[k].second;
                        axpy_(pr, a, row);
                        break;
                    }
            }
            rr.pivot_row_of[static_cast<std::size_t>(pc)] = static_cast<int>(rr.rows.size());
            rr.rows.push_back(row);
            rr.pivot.push_back(pc);
        }
        return rr;
    }

    static detail::bigint lcm_(const detail::bigint& a, const detail::bigint& b) {
        detail::bigint g = detail::bigint::gcd(a, b), q, r;
        detail::bigint::divmod(a, g, q, r);
        return q * b;
    }

    // z_f = e_f - sum_p R_p[f] e_{pivot(p)}, scaled to integers
    void build_nullspace_(const rref_result& rr) {
        std::vector<sparse_vec> per_free(static_cast<std::size_t>(ncols_));
        for (std::size_t q = 0; q < rr.rows.size(); ++q) {
            const sparse_vec& pr = rr.rows[q];
            const int pc = rr.pivot[q];
            for (std::size_t k = 0; k < pr.size(); ++k) {
                if (pr[k].first == pc) continue;
                per_free[static_cast<std::size_t>(pr[k].first)].push_back(
                    std::make_pair(pc, -pr[k].second));
            }
        }
        cols_.clear();
        free_col_.clear();
        for (int f = 0; f < ncols_; ++f) {
            if (rr.pivot_row_of[static_cast<std::size_t>(f)] >= 0) continue;
            sparse_vec z = per_free[static_cast<std::size_t>(f)];
            z.push_back(std::make_pair(f, R(1)));
            std::sort(z.begin(), z.end(), by_index());
            detail::bigint L(1);
            for (std::size_t k = 0; k < z.size(); ++k) L = lcm_(L, z[k].second.den());
            const R Lr(L, detail::bigint(1));
            for (std::size_t k = 0; k < z.size(); ++k) z[k].second = z[k].second * Lr;
            cols_.push_back(z);
            free_col_.push_back(f);
        }
    }

    // certificate (section (ii)): computed from the finished B and Z only
    void certify_() {
        const int N = dim();
        cert_.rows = nt_;
        cert_.columns = ncols_;
        cert_.rank = rank_;
        cert_.rank_deficiency = nt_ - rank_;
        cert_.dim = N;
        cert_.dim_identity = (N == ncols_ - rank_);
        cert_.pairs_checked = 0;
        cert_.violations = 0;
        // (a) every row of B against every column of Z
        for (int e = 0; e < nt_; ++e) {
            const sparse_vec& row = rows_[static_cast<std::size_t>(e)];
            for (int c = 0; c < N; ++c) {
                const sparse_vec& z = cols_[static_cast<std::size_t>(c)];
                R s(0);
                std::size_t i = 0, j = 0;
                while (i < row.size() && j < z.size()) {
                    if (row[i].first < z[j].first) ++i;
                    else if (z[j].first < row[i].first) ++j;
                    else { s += row[i].second * z[j].second; ++i; ++j; }
                }
                ++cert_.pairs_checked;
                if (!s.is_zero()) ++cert_.violations;
            }
        }
        // (c) independence witness: the unit entry of column c sits in its
        // free column, and no other column has an entry there
        bool wit = true;
        std::vector<int> owner(static_cast<std::size_t>(ncols_), -1);
        for (int c = 0; c < N && wit; ++c) {
            const int f = free_col_[static_cast<std::size_t>(c)];
            if (owner[static_cast<std::size_t>(f)] != -1) { wit = false; break; }
            owner[static_cast<std::size_t>(f)] = c;
            const sparse_vec& z = cols_[static_cast<std::size_t>(c)];
            bool found = false;
            for (std::size_t k = 0; k < z.size(); ++k)
                if (z[k].first == f) { found = !z[k].second.is_zero(); break; }
            if (!found) wit = false;
        }
        for (int c = 0; c < N && wit; ++c) {
            const sparse_vec& z = cols_[static_cast<std::size_t>(c)];
            for (std::size_t k = 0; k < z.size(); ++k) {
                const int o = owner[static_cast<std::size_t>(z[k].first)];
                if (o != -1 && o != c) { wit = false; break; }
            }
        }
        cert_.independence_witness = wit;
    }

    int nt_, n_reduced_, ncols_, rank_;
    std::vector<sparse_vec> rows_;        // B, nt rows over D * n_r columns
    std::vector<sparse_vec> cols_;        // Z, dim columns (integer entries)
    std::vector<int> free_col_;
    vcr1_divfree_certificate cert_;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CR1_VCR1_DIVFREE_HPP
