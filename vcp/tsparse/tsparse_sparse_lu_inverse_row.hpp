// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-INVROW: single-row / single-column extraction of the approximate
// inverse from a baseline-CSC sparse LU factor.
//
// Factor pipeline (see detail/tsparse_sparse_lu_solve_impl.hpp,
// solve_baseline_storage):
//
//     x = Dc . Q . U^{-1} . L^{-1} . P . Dr . b,
//     (P v)[new] = v[row_perm[new]],   (Q z)[col_perm[new]] = z[new],
//
// hence   A^{-1} ~= Dc Q U^{-1} L^{-1} P Dr   on original coordinates.
//
// ROW i (adjoint direction -- two TRANSPOSE triangular solves):
//     row_i(A^{-1}) = e_i^T A^{-1} = Dc[i] . (e_jc^T U^{-1} L^{-1}) P Dr,
//     jc = inv_col_perm[i]
//   stage 1:  U^T y = e_jc   (U^T lower triangular, explicit diagonal)
//   stage 2:  L^T z = y      (L^T upper triangular, unit diagonal)
//   exit:     row_i(A^{-1})[row_perm[k]] = Dc[i] * Dr[row_perm[k]] * z[k]
// U^T, L^T are built ONCE at extractor construction (CSC transpose ==
// storing U, L row-wise), after which both stages are plain Gilbert-Peierls
// scatter solves.
//
// COLUMN j (forward direction -- NO transposes needed; the original CSC
// L/U already have the scatter orientation):
//     col_j(A^{-1}) = A^{-1} e_j = Dc . Q . U^{-1} L^{-1} . (Dr[j] e_ir),
//     ir = inv_row_perm[j]
//   stage 1:  L y = Dr[j] e_ir   (unit diagonal, forward)
//   stage 2:  U z = y            (explicit diagonal, backward)
//   exit:     col_j(A^{-1})[col_perm[k]] = Dc[col_perm[k]] * z[k]
//
// CONTRACT (slice locality): the nonzero pattern of the result is confined
// to the graph reach of the seed(s) on the corresponding triangular factor
// graphs (symbolic DFS computed before any arithmetic).  Components outside
// the reach are NEVER touched: no other row/column of A^{-1} is formed, not
// even as an intermediate.  Cost is O(flops actually performed) = O(sum of
// column lengths over the reach), independent of n and of nnz(L)+nnz(U)
// (Gilbert-Peierls).  The pattern can still grow toward [0, n) when the
// requested slice of A^{-1} is itself dense (a property of A, not of the
// algorithm) -- this is the documented "not completely forbidden" caveat.
//
// AUTOMATIC FAST PATHS (no special casing; certified by the reach counters):
//   row i with jc == n-1:    column n-1 of U^T holds only the diagonal
//                            -> reach_u == 1, y = e_{n-1}/u_{nn}
//   column j with ir == n-1: column n-1 of L holds only the diagonal
//                            -> reach_l == 1, y = Dr[j] e_{n-1}
// Under natural ordering these are exactly the LAST row / LAST column.
//
// STORAGE REQUIREMENTS: baseline CSC L/U must be present (baseline_storage()).
// method=supernodal factors are accepted iff the CSC fallback is retained
// (same condition as solve() dispatch priority 3); otherwise state_error.
// U's diagonal is validated (present and nonzero) ONCE at extractor
// construction; per-extract solves therefore cannot hit a pivot error.
//
// SCALAR GATE: same as the SLU-8/14R module scope -- point arithmetic types
// {float, double, long double, kv::dd, kv::mpfr<N>}.  interval<...> scalars
// are outside the scope of this module (verified enclosures of inverse
// slices belong to the verified layer, on top of this approximate kernel).
//
// COORDINATES: input indices and output patterns are in ORIGINAL
// (unpermuted, unscaled) coordinates.
//
// LIFETIME: extractors copy everything they need (permutations, scalings,
// factor storages / transposes); after construction they are independent of
// the factorization object's lifetime.

#pragma once

#ifndef VCP_TSPARSE_SPARSE_LU_INVERSE_ROW_HPP
#define VCP_TSPARSE_SPARSE_LU_INVERSE_ROW_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_sparse_lu.hpp>

namespace vcp {

// ===========================================================================
// Result structs.  index[] is sorted ascending; index/value are in original
// coordinates.  The pattern is structural: it may contain exact numerical
// zeros (cancellation), never omits a structurally reachable component.
// ===========================================================================

template <class T, class Index>
struct sparse_lu_inverse_row_result {
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_inverse_row: Index must be signed");

    Index row;                   // requested row index (original coordinates)
    std::vector<Index> index;    // column indices, sorted ascending
    std::vector<T>     value;    // matching values

    Index reach_u;               // |Reach_{G(U^T)}(jc)|      (stage-1 pattern)
    Index reach_l;               // |Reach_{G(L^T)}(pat(y))|  (== index.size())
    std::size_t flops;           // off-diagonal multiply-subtract count

    sparse_lu_inverse_row_result()
        : row(Index(0)), reach_u(Index(0)), reach_l(Index(0)), flops(0) {}

    // Densify into a length-n vector (zeros outside the pattern).
    std::vector<T> to_dense(Index n) const {
        std::vector<T> v(static_cast<std::size_t>(n), T(0));
        for (std::size_t k = 0; k < index.size(); ++k) {
            v[static_cast<std::size_t>(index[k])] = value[k];
        }
        return v;
    }

    // Build a 1 x n sparse matrix.  Matrix must provide resize/add/finalize
    // (vcp::spmatrix does); template keeps this header spmatrix-independent.
    template <class Matrix>
    Matrix to_matrix(Index n) const {
        Matrix M;
        M.resize(typename Matrix::index_type(1),
                 static_cast<typename Matrix::index_type>(n));
        for (std::size_t k = 0; k < index.size(); ++k) {
            M.add(typename Matrix::index_type(0),
                  static_cast<typename Matrix::index_type>(index[k]),
                  value[k]);
        }
        M.finalize();
        return M;
    }
};

template <class T, class Index>
struct sparse_lu_inverse_col_result {
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_inverse_col: Index must be signed");

    Index col;                   // requested column index (original coordinates)
    std::vector<Index> index;    // row indices, sorted ascending
    std::vector<T>     value;    // matching values

    Index reach_l;               // |Reach_{G(L)}(ir)|        (stage-1 pattern)
    Index reach_u;               // |Reach_{G(U)}(pat(y))|    (== index.size())
    std::size_t flops;           // off-diagonal multiply-subtract count

    sparse_lu_inverse_col_result()
        : col(Index(0)), reach_l(Index(0)), reach_u(Index(0)), flops(0) {}

    std::vector<T> to_dense(Index n) const {
        std::vector<T> v(static_cast<std::size_t>(n), T(0));
        for (std::size_t k = 0; k < index.size(); ++k) {
            v[static_cast<std::size_t>(index[k])] = value[k];
        }
        return v;
    }

    // Build an n x 1 sparse matrix.
    template <class Matrix>
    Matrix to_matrix(Index n) const {
        Matrix M;
        M.resize(static_cast<typename Matrix::index_type>(n),
                 typename Matrix::index_type(1));
        for (std::size_t k = 0; k < index.size(); ++k) {
            M.add(static_cast<typename Matrix::index_type>(index[k]),
                  typename Matrix::index_type(0),
                  value[k]);
        }
        M.finalize();
        return M;
    }
};

namespace sparse_lu_detail {

// ===========================================================================
// inverse_slice_extractor_base<T, Index>
// Shared machinery for the row and column extractors:
//   - factor validation / baseline availability check
//   - CSC structural transpose, triangularity check, diagonal precompute
//   - symbolic reach (iterative DFS, postorder) with persistent mark scratch
//   - the two numeric scatter kernels (unit diagonal / explicit diagonal)
//   - touched-list scratch reset (O(reach) per extract, never O(n))
// ===========================================================================
template <class T, class Index>
class inverse_slice_extractor_base {
protected:
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu inverse slice: Index must be signed");

    inverse_slice_extractor_base() : n_(Index(0)) {}

    // ---- factor validation common to both extractors ----
    static const baseline_lu_storage<T, Index>&
    checked_baseline_(const sparse_lu_factorization<T, Index>& fac,
                      const char* who, Index& n_out)
    {
        if (!fac.valid()) {
            vcp::throw_error<vcp::state_error>(
                who, ": factor is not valid (status=",
                sparse_lu_status_to_string(fac.info().status), ")");
        }
        n_out = fac.info().n;
        const baseline_lu_storage<T, Index>& lu = fac.baseline_storage();

        // Same availability condition as solve() dispatch priority 3.
        if (n_out > Index(0) && lu.L.col_ptr.empty()) {
            vcp::throw_error<vcp::state_error>(
                who, ": baseline CSC L/U not available (supernodal-only "
                "storage; inverse-slice extraction requires the CSC factors)");
        }
        validate_baseline_lu_storage_for_solve(lu, n_out);
        return lu;
    }

    void init_scratch_(const Index n) {
        n_ = n;
        const std::size_t un = static_cast<std::size_t>(n);
        x_.assign(un, T(0));
        mark_.assign(un, static_cast<unsigned char>(0));
        dfs_node_.reserve(un);
        dfs_edge_.reserve(un);
    }

    // ---- structural transpose: CSC(M) -> CSC(M^T), counting sort.
    // Rows within each output column come out sorted ascending. ----
    static void transpose_csc_(const csc_storage<T, Index>& M,
                               const Index n,
                               csc_storage<T, Index>& MT)
    {
        const std::size_t un  = static_cast<std::size_t>(n);
        const std::size_t nnz = M.row_ind.size();

        MT.col_ptr.assign(un + 1u, Index(0));
        MT.row_ind.assign(nnz, Index(0));
        MT.values.assign(nnz, T(0));

        for (std::size_t p = 0; p < nnz; ++p) {
            ++MT.col_ptr[static_cast<std::size_t>(M.row_ind[p]) + 1u];
        }
        for (std::size_t c = 0; c < un; ++c) {
            MT.col_ptr[c + 1u] += MT.col_ptr[c];
        }
        std::vector<Index> next(MT.col_ptr.begin(), MT.col_ptr.end() - 1);
        for (Index j = Index(0); j < n; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            for (Index p = M.col_ptr[sj]; p < M.col_ptr[sj + 1u]; ++p) {
                const std::size_t sp = static_cast<std::size_t>(p);
                const std::size_t r = static_cast<std::size_t>(M.row_ind[sp]);
                const std::size_t q = static_cast<std::size_t>(next[r]++);
                MT.row_ind[q] = j;          // column of M = row of M^T
                MT.values[q]  = M.values[sp];
            }
        }
    }

    static void check_triangular_(const csc_storage<T, Index>& M,
                                  const bool lower,
                                  const char* msg)
    {
        const Index n = static_cast<Index>(M.col_ptr.empty()
            ? 0 : M.col_ptr.size() - 1u);
        for (Index j = Index(0); j < n; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            for (Index p = M.col_ptr[sj]; p < M.col_ptr[sj + 1u]; ++p) {
                const Index r = M.row_ind[static_cast<std::size_t>(p)];
                if ((lower && r < j) || (!lower && r > j)) {
                    vcp::throw_error<vcp::invalid_argument>(msg);
                }
            }
        }
    }

    // ---- one-time explicit-diagonal validation + position table.
    // diag_pos[k] = index into M.values of the diagonal of column k.
    // Throws (zero_pivot / numerical_singularity semantics) if a diagonal
    // is missing or exactly zero, so per-extract solves cannot fail. ----
    static std::vector<Index>
    build_explicit_diag_positions_(const csc_storage<T, Index>& M,
                                   const Index n, const char* who)
    {
        std::vector<Index> diag_pos(static_cast<std::size_t>(n), Index(-1));
        for (Index k = Index(0); k < n; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            for (Index p = M.col_ptr[sk]; p < M.col_ptr[sk + 1u]; ++p) {
                if (M.row_ind[static_cast<std::size_t>(p)] == k) {
                    diag_pos[sk] = p;
                    break;
                }
            }
            if (diag_pos[sk] < Index(0)) {
                vcp::throw_error<vcp::invalid_argument>(
                    who, ": missing diagonal in U (zero_pivot)");
            }
            if (sparse_lu_scalar_policy<T>::is_exact_zero(
                    M.values[static_cast<std::size_t>(diag_pos[sk])])) {
                vcp::throw_error<vcp::state_error>(
                    who, ": zero diagonal in U (numerical_singularity)");
            }
        }
        return diag_pos;
    }

    // ---- symbolic reach: iterative DFS from `seed` over the off-diagonal
    // entries of G(M) (out-edges of node k = off-diagonal row indices of
    // column k).  Appends the reach to `topo` in DFS POSTORDER; consuming
    // `topo` back-to-front yields a topological order of the reach DAG
    // (triangularity of M guarantees acyclicity).  mark_ persists across
    // the two stages of one extract and is reset via the touched lists. ----
    void reach_(const csc_storage<T, Index>& M,
                const Index seed,
                std::vector<Index>& topo)
    {
        if (mark_[static_cast<std::size_t>(seed)]) return;

        dfs_node_.clear();
        dfs_edge_.clear();
        dfs_node_.push_back(seed);
        dfs_edge_.push_back(M.col_ptr[static_cast<std::size_t>(seed)]);
        mark_[static_cast<std::size_t>(seed)] = 1;

        while (!dfs_node_.empty()) {
            const Index k = dfs_node_.back();
            const std::size_t sk = static_cast<std::size_t>(k);
            Index& p = dfs_edge_.back();

            bool descended = false;
            while (p < M.col_ptr[sk + 1u]) {
                const Index r = M.row_ind[static_cast<std::size_t>(p)];
                ++p;
                if (r == k) continue;                       // diagonal: no edge
                if (mark_[static_cast<std::size_t>(r)]) continue;
                mark_[static_cast<std::size_t>(r)] = 1;
                dfs_node_.push_back(r);
                dfs_edge_.push_back(M.col_ptr[static_cast<std::size_t>(r)]);
                descended = true;
                break;
            }
            if (!descended) {
                topo.push_back(k);                          // postorder finish
                dfs_node_.pop_back();
                dfs_edge_.pop_back();
            }
        }
    }

    // ---- numeric kernels: process the reach in reverse postorder,
    // scattering column k's off-diagonal contributions out of x_[k].
    // Only components inside the reach are read or written. ----

    // unit-diagonal factor (L-type): stored diagonal entries silently
    // ignored (unit contract, matching csc_forward_solve_L).
    void solve_unit_diag_(const csc_storage<T, Index>& M,
                          const std::vector<Index>& topo,
                          std::size_t& flops)
    {
        for (std::size_t t = topo.size(); t-- > 0; ) {
            const Index k = topo[t];
            const std::size_t sk = static_cast<std::size_t>(k);
            for (Index p = M.col_ptr[sk]; p < M.col_ptr[sk + 1u]; ++p) {
                const std::size_t sp = static_cast<std::size_t>(p);
                const Index r = M.row_ind[sp];
                if (r == k) continue;
                x_[static_cast<std::size_t>(r)] -= M.values[sp] * x_[sk];
                ++flops;
            }
        }
    }

    // explicit-diagonal factor (U-type): divide by the (pre-validated)
    // diagonal, then scatter.
    void solve_explicit_diag_(const csc_storage<T, Index>& M,
                              const std::vector<Index>& diag_pos,
                              const std::vector<Index>& topo,
                              std::size_t& flops)
    {
        for (std::size_t t = topo.size(); t-- > 0; ) {
            const Index k = topo[t];
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index dp = diag_pos[sk];

            x_[sk] /= M.values[static_cast<std::size_t>(dp)];

            for (Index p = M.col_ptr[sk]; p < M.col_ptr[sk + 1u]; ++p) {
                if (p == dp) continue;
                const std::size_t sp = static_cast<std::size_t>(p);
                x_[static_cast<std::size_t>(M.row_ind[sp])]
                    -= M.values[sp] * x_[sk];
                ++flops;
            }
        }
    }

    void clear_marks_(const std::vector<Index>& touched) {
        for (std::size_t s = 0; s < touched.size(); ++s) {
            mark_[static_cast<std::size_t>(touched[s])] = 0;
        }
    }

    void clear_scratch_(const std::vector<Index>& touched) {
        for (std::size_t s = 0; s < touched.size(); ++s) {
            const std::size_t k = static_cast<std::size_t>(touched[s]);
            x_[k]    = T(0);
            mark_[k] = 0;
        }
    }

    static void sort_pattern_(std::vector<Index>& index, std::vector<T>& value) {
        const std::size_t m = index.size();
        std::vector<std::size_t> ord(m);
        for (std::size_t s = 0; s < m; ++s) ord[s] = s;
        const std::vector<Index>& idx = index;
        std::sort(ord.begin(), ord.end(),
                  [&idx](const std::size_t a, const std::size_t b) {
                      return idx[a] < idx[b];
                  });
        std::vector<Index> si(m);
        std::vector<T>     sv(m);
        for (std::size_t s = 0; s < m; ++s) {
            si[s] = index[ord[s]];
            sv[s] = value[ord[s]];
        }
        index.swap(si);
        value.swap(sv);
    }

    Index n_;

    // per-extract scratch (allocated once; reset via touched lists)
    std::vector<T> x_;
    std::vector<unsigned char> mark_;
    std::vector<Index> topo1_, topo2_;
    std::vector<Index> dfs_node_;
    std::vector<Index> dfs_edge_;
};

} // namespace sparse_lu_detail

// ===========================================================================
// sparse_lu_inverse_row_extractor<T, Index>
//
// Construction cost: O(n + nnz(L) + nnz(U)) -- one structural/numerical
// transpose of each factor, U-diagonal validation, O(n) scratch.  Each
// extract(i) then costs O(flops of the two GP solves) + O(|pattern| log)
// for the final sort.
// ===========================================================================
template <class T, class Index>
class sparse_lu_inverse_row_extractor
    : private sparse_lu_detail::inverse_slice_extractor_base<T, Index> {

    typedef sparse_lu_detail::inverse_slice_extractor_base<T, Index> base;

public:
    typedef T     value_type;
    typedef Index index_type;

    explicit sparse_lu_inverse_row_extractor(
        const sparse_lu_factorization<T, Index>& fac)
    {
        Index n = Index(0);
        const baseline_lu_storage<T, Index>& lu =
            base::checked_baseline_(fac, "sparse_lu_inverse_row", n);

        // Triangularity validated ONCE (the per-call solves then trust the
        // transposed structures).
        base::check_triangular_(lu.L, /*lower=*/true,
            "sparse_lu_inverse_row: upper entry in L (invalid factor storage)");
        base::check_triangular_(lu.U, /*lower=*/false,
            "sparse_lu_inverse_row: lower entry in U (invalid factor storage)");

        // U^T in CSC == U row-wise: column k of UT_ = row k of U (lower tri).
        // L^T in CSC == L row-wise: column k of LT_ = row k of L (upper tri).
        base::transpose_csc_(lu.U, n, UT_);
        base::transpose_csc_(lu.L, n, LT_);
        udiag_pos_ = base::build_explicit_diag_positions_(
            UT_, n, "sparse_lu_inverse_row");

        row_perm_ = lu.row_perm;
        Dr_ = lu.Dr;
        Dc_ = lu.Dc;
        if (!lu.inv_col_perm.empty()) {
            inv_col_perm_ = lu.inv_col_perm;
        } else if (!lu.col_perm.empty()) {
            inv_col_perm_ = sparse_lu_inverse_permutation(lu.col_perm);
        }

        base::init_scratch_(n);
    }

    Index size() const { return this->n_; }

    // Extract row i (original coordinates) of the approximate inverse.
    sparse_lu_inverse_row_result<T, Index> extract(const Index i) {
        if (i < Index(0) || i >= this->n_) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_inverse_row: row index out of range");
        }

        sparse_lu_inverse_row_result<T, Index> result;
        result.row = i;

        // seed:  jc = inv_col_perm[i]  (e_i^T Q = e_jc^T)
        const Index jc = inv_col_perm_.empty()
            ? i : inv_col_perm_[static_cast<std::size_t>(i)];

        // Stage 1:  U^T y = e_jc  (lower triangular, explicit diagonal)
        this->topo1_.clear();
        this->reach_(UT_, jc, this->topo1_);
        result.reach_u = static_cast<Index>(this->topo1_.size());

        this->x_[static_cast<std::size_t>(jc)] = T(1);
        this->solve_explicit_diag_(UT_, udiag_pos_, this->topo1_, result.flops);

        // Stage 2:  L^T z = y  (upper triangular, unit diagonal);
        // seeds = pattern(y); z overwrites y in-place in x_.
        this->clear_marks_(this->topo1_);
        this->topo2_.clear();
        for (std::size_t s = 0; s < this->topo1_.size(); ++s) {
            this->reach_(LT_, this->topo1_[s], this->topo2_);
        }
        result.reach_l = static_cast<Index>(this->topo2_.size());
        this->solve_unit_diag_(LT_, this->topo2_, result.flops);

        // Exit:  row_i(A^{-1})[row_perm[k]] = Dc[i] * Dr[row_perm[k]] * z[k]
        const T dci = Dc_.empty() ? T(1) : Dc_[static_cast<std::size_t>(i)];

        result.index.reserve(this->topo2_.size());
        result.value.reserve(this->topo2_.size());
        for (std::size_t s = 0; s < this->topo2_.size(); ++s) {
            const Index k = this->topo2_[s];
            const Index col = row_perm_.empty()
                ? k : row_perm_[static_cast<std::size_t>(k)];
            T v = dci * this->x_[static_cast<std::size_t>(k)];
            if (!Dr_.empty()) v *= Dr_[static_cast<std::size_t>(col)];
            result.index.push_back(col);
            result.value.push_back(v);
        }

        this->clear_scratch_(this->topo2_);
        base::sort_pattern_(result.index, result.value);
        return result;
    }

    // Last row: for natural column ordering this hits the automatic fast
    // path (reach_u == 1).
    sparse_lu_inverse_row_result<T, Index> extract_last() {
        if (this->n_ == Index(0)) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_inverse_row: extract_last on empty factor");
        }
        return extract(this->n_ - Index(1));
    }

private:
    csc_storage<T, Index> UT_;          // U^T (CSC) == U row-wise; lower tri
    csc_storage<T, Index> LT_;          // L^T (CSC) == L row-wise; upper tri
    std::vector<Index> udiag_pos_;      // diag position per UT_ column
    std::vector<Index> row_perm_;       // empty == identity
    std::vector<Index> inv_col_perm_;   // empty == identity
    std::vector<T> Dr_, Dc_;            // empty == identity
};

// ===========================================================================
// sparse_lu_inverse_col_extractor<T, Index>
//
// col_j(A^{-1}) = A^{-1} e_j: the FORWARD sparse-RHS solve.  The original
// CSC L/U already have the scatter orientation, so no transposes are built;
// construction copies L/U (lifetime independence) and validates once.
// ===========================================================================
template <class T, class Index>
class sparse_lu_inverse_col_extractor
    : private sparse_lu_detail::inverse_slice_extractor_base<T, Index> {

    typedef sparse_lu_detail::inverse_slice_extractor_base<T, Index> base;

public:
    typedef T     value_type;
    typedef Index index_type;

    explicit sparse_lu_inverse_col_extractor(
        const sparse_lu_factorization<T, Index>& fac)
    {
        Index n = Index(0);
        const baseline_lu_storage<T, Index>& lu =
            base::checked_baseline_(fac, "sparse_lu_inverse_col", n);

        base::check_triangular_(lu.L, /*lower=*/true,
            "sparse_lu_inverse_col: upper entry in L (invalid factor storage)");
        base::check_triangular_(lu.U, /*lower=*/false,
            "sparse_lu_inverse_col: lower entry in U (invalid factor storage)");

        L_ = lu.L;                       // owned copies: lifetime independence
        U_ = lu.U;
        udiag_pos_ = base::build_explicit_diag_positions_(
            U_, n, "sparse_lu_inverse_col");

        col_perm_ = lu.col_perm;
        Dr_ = lu.Dr;
        Dc_ = lu.Dc;
        if (!lu.inv_row_perm.empty()) {
            inv_row_perm_ = lu.inv_row_perm;
        } else if (!lu.row_perm.empty()) {
            inv_row_perm_ = sparse_lu_inverse_permutation(lu.row_perm);
        }

        base::init_scratch_(n);
    }

    Index size() const { return this->n_; }

    // Extract column j (original coordinates) of the approximate inverse.
    sparse_lu_inverse_col_result<T, Index> extract(const Index j) {
        if (j < Index(0) || j >= this->n_) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_inverse_col: column index out of range");
        }

        sparse_lu_inverse_col_result<T, Index> result;
        result.col = j;

        // seed:  P Dr e_j = Dr[j] e_ir,  ir = inv_row_perm[j]
        const Index ir = inv_row_perm_.empty()
            ? j : inv_row_perm_[static_cast<std::size_t>(j)];

        // Stage 1:  L y = Dr[j] e_ir  (unit diagonal, forward)
        this->topo1_.clear();
        this->reach_(L_, ir, this->topo1_);
        result.reach_l = static_cast<Index>(this->topo1_.size());

        this->x_[static_cast<std::size_t>(ir)] =
            Dr_.empty() ? T(1) : Dr_[static_cast<std::size_t>(j)];
        this->solve_unit_diag_(L_, this->topo1_, result.flops);

        // Stage 2:  U z = y  (explicit diagonal, backward);
        // seeds = pattern(y); z overwrites y in-place in x_.
        this->clear_marks_(this->topo1_);
        this->topo2_.clear();
        for (std::size_t s = 0; s < this->topo1_.size(); ++s) {
            this->reach_(U_, this->topo1_[s], this->topo2_);
        }
        result.reach_u = static_cast<Index>(this->topo2_.size());
        this->solve_explicit_diag_(U_, udiag_pos_, this->topo2_, result.flops);

        // Exit:  col_j(A^{-1})[col_perm[k]] = Dc[col_perm[k]] * z[k]
        result.index.reserve(this->topo2_.size());
        result.value.reserve(this->topo2_.size());
        for (std::size_t s = 0; s < this->topo2_.size(); ++s) {
            const Index k = this->topo2_[s];
            const Index row = col_perm_.empty()
                ? k : col_perm_[static_cast<std::size_t>(k)];
            T v = this->x_[static_cast<std::size_t>(k)];
            if (!Dc_.empty()) v *= Dc_[static_cast<std::size_t>(row)];
            result.index.push_back(row);
            result.value.push_back(v);
        }

        this->clear_scratch_(this->topo2_);
        base::sort_pattern_(result.index, result.value);
        return result;
    }

    // Last column: for natural row ordering this hits the automatic fast
    // path (reach_l == 1).
    sparse_lu_inverse_col_result<T, Index> extract_last() {
        if (this->n_ == Index(0)) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_inverse_col: extract_last on empty factor");
        }
        return extract(this->n_ - Index(1));
    }

private:
    csc_storage<T, Index> L_;           // owned copy; lower tri, unit diag
    csc_storage<T, Index> U_;           // owned copy; upper tri, explicit diag
    std::vector<Index> udiag_pos_;      // diag position per U_ column
    std::vector<Index> col_perm_;       // empty == identity
    std::vector<Index> inv_row_perm_;   // empty == identity
    std::vector<T> Dr_, Dc_;            // empty == identity
};

// ===========================================================================
// One-shot conveniences.  Prefer the extractor classes when extracting
// several slices from the same factor (amortizes construction).
// ===========================================================================

template <class T, class Index>
sparse_lu_inverse_row_result<T, Index>
sparse_lu_inverse_row(const sparse_lu_factorization<T, Index>& fac,
                      const Index i)
{
    sparse_lu_inverse_row_extractor<T, Index> ex(fac);
    return ex.extract(i);
}

template <class T, class Index>
sparse_lu_inverse_row_result<T, Index>
sparse_lu_inverse_last_row(const sparse_lu_factorization<T, Index>& fac)
{
    sparse_lu_inverse_row_extractor<T, Index> ex(fac);
    return ex.extract_last();
}

template <class T, class Index>
sparse_lu_inverse_col_result<T, Index>
sparse_lu_inverse_col(const sparse_lu_factorization<T, Index>& fac,
                      const Index j)
{
    sparse_lu_inverse_col_extractor<T, Index> ex(fac);
    return ex.extract(j);
}

template <class T, class Index>
sparse_lu_inverse_col_result<T, Index>
sparse_lu_inverse_last_col(const sparse_lu_factorization<T, Index>& fac)
{
    sparse_lu_inverse_col_extractor<T, Index> ex(fac);
    return ex.extract_last();
}

// From a matrix: factorize (strict) then extract.  Matrix follows the same
// concept as sparse_lu_factorize (vcp::spmatrix<T, P> qualifies).
template <class Matrix>
sparse_lu_inverse_row_result<typename Matrix::value_type,
                             typename Matrix::index_type>
sparse_lu_inverse_row(
    const Matrix& A,
    const typename Matrix::index_type i,
    const sparse_lu_options<typename Matrix::value_type>& opt
        = sparse_lu_options<typename Matrix::value_type>())
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    sparse_lu_factorization<T, Index> fac = sparse_lu_factorize(A, opt);
    sparse_lu_inverse_row_extractor<T, Index> ex(fac);
    return ex.extract(i);
}

template <class Matrix>
sparse_lu_inverse_col_result<typename Matrix::value_type,
                             typename Matrix::index_type>
sparse_lu_inverse_col(
    const Matrix& A,
    const typename Matrix::index_type j,
    const sparse_lu_options<typename Matrix::value_type>& opt
        = sparse_lu_options<typename Matrix::value_type>())
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    sparse_lu_factorization<T, Index> fac = sparse_lu_factorize(A, opt);
    sparse_lu_inverse_col_extractor<T, Index> ex(fac);
    return ex.extract(j);
}

} // namespace vcp

#endif // VCP_TSPARSE_SPARSE_LU_INVERSE_ROW_HPP
