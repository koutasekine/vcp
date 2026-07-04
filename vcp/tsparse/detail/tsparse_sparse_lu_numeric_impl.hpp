// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-4.1 dense-reference + SLU-5 sparse GP baseline numeric factorization.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER all SLU-1/2/3
// helpers (csc_storage, baseline_lu_storage, sparse_lu_scalar_policy,
// sparse_lu_is_acceptable_pivot, sparse_lu_identity_permutation,
// sparse_lu_inverse_permutation, baseline_lu_storage) are in scope.
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse.hpp>  (umbrella)

#ifndef VCP_TSPARSE_SPARSE_LU_NUMERIC_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_NUMERIC_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <vector>

#include <vcp/error.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// sparse_lu_column_reach_from_L  (SLU-6: DFS-based column reach helper)
//
// Computes the set of columns k < j needed for the GP update of column j,
// using the L structure built so far (columns 0..j-1).
//
// Seeds: pattern entries with index i < j (active rows from scatter of A[:,j]).
// DFS:   from each seed k, follow L[:,k] rows m with m < j.  Those m can
//        become active in the SPA through the update x[m] -= L[m,k]*x[k],
//        so they also need their L column applied.  Repeat transitively.
// Result: sorted in increasing order — the correct topological order for the
//         GP update (processing k1 < k2 never requires result of k2 first).
//
// reach_mark: persistent generation-marker array, size n, initialized to -1.
//             reach_mark[k] == j means column k is already in reach for col j.
//             Reused across columns; no re-allocation per column.
//
// Correctness note:
//   If m < j and L[m,k] != 0 for some k already in reach, then m is also in
//   reach (found by DFS from k).  Hence every column that can become active
//   during the update is included.  Extra nodes in reach that are not
//   numerically active are harmless — the guard (mark[sk] != j) in the caller
//   skips them.
//
// Conservative scan fallback (for k=0..j-1 unconditionally) has been removed
// from the production default (SLU-5) and replaced by this helper.  See
// baseline_reference_lu_factorize (TEST-ONLY) for the naive dense reference.
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index>
sparse_lu_column_reach_from_L(
    Index                      n,
    Index                      j,
    const std::vector<Index>&  L_col_ptr_buf,
    const std::vector<Index>&  L_row_ind_buf,
    const std::vector<Index>&  pattern,
    std::vector<Index>&        reach_mark)
{
    (void)n;

    if (j == Index(0)) {
        return std::vector<Index>();
    }

    std::vector<Index> reach;
    std::vector<Index> dfs_stack;

    for (std::size_t ai = 0u; ai < pattern.size(); ++ai) {
        const Index i = pattern[ai];
        if (i < Index(0) || i >= j) continue;
        const std::size_t si = static_cast<std::size_t>(i);
        if (reach_mark[si] == j) continue;
        dfs_stack.push_back(i);

        while (!dfs_stack.empty()) {
            const Index k = dfs_stack.back();
            dfs_stack.pop_back();
            const std::size_t sk = static_cast<std::size_t>(k);
            if (reach_mark[sk] == j) continue;
            reach_mark[sk] = j;
            reach.push_back(k);

            // Follow L[:,k]: rows m with m < j can activate fill in the SPA
            for (Index lp = L_col_ptr_buf[sk]; lp < L_col_ptr_buf[sk + 1u]; ++lp) {
                const std::size_t slp = static_cast<std::size_t>(lp);
                const Index m = L_row_ind_buf[slp];
                if (m < j) {
                    const std::size_t sm = static_cast<std::size_t>(m);
                    if (reach_mark[sm] != j) {
                        dfs_stack.push_back(m);
                    }
                }
            }
        }
    }

    std::sort(reach.begin(), reach.end());
    return reach;
}

// ---------------------------------------------------------------------------
// baseline_reference_factorize_result: shared output type for both helpers
// ---------------------------------------------------------------------------
template <class T, class Index>
struct baseline_reference_factorize_result {
    bool                          success;
    sparse_lu_status              status;
    baseline_lu_storage<T, Index> storage;

    baseline_reference_factorize_result()
        : success(false), status(sparse_lu_status::numerical_singularity) {}
};

// ---------------------------------------------------------------------------
// baseline_sparse_gp_lu_factorize  (SLU-5 production path, SLU-6 updated)
//
// Sparse Gilbert-Peierls / SPA baseline numeric factorization.
//
// SPA workspace:
//   x[n]         : numeric accumulator, zero between columns
//   mark[n]      : generation marker; mark[i] == j iff row i is active in col j
//   pattern      : list of active row indices for current column
//   reach_mark[n]: DFS reach generation marker (SLU-6); reused across columns
//
// Algorithm (per column j):
//   1. Scatter A[:,j] into SPA via inv_row_perm (original row → pivot position).
//   2. SLU-6 DFS reach-based update: compute reach set from L-structure DFS
//      starting at seeds (pattern rows < j), then apply L[:,k] for each active
//      k in reach in increasing order.
//      SLU-5 conservative scan (for k=0..j-1) removed from production default;
//      retained only in baseline_reference_lu_factorize (TEST-ONLY reference).
//   3. Pivot selection: argmax_{i>=j, marked} |x[i]| (partial) or use j (none).
//   4. Pivot acceptability check via sparse_lu_is_acceptable_pivot (SLU-3).
//   5. Row swap: swap x[j] <-> x[pivot_pos], update row_perm + inv_row_perm.
//   6. Store U column j: x[0..j] (diagonal always; zero strict-upper skipped).
//   7. Store L column j: x[j+1..n-1] / x[j] (zero entries skipped).
//   8. Clear SPA.
//
// Preconditions (caller ensures):
//   - A_csc is a valid n×n CSC storage.
//   - opt pivot parameters have been validated.
//   - col_perm / inv_col_perm are valid permutations of size n.
// ---------------------------------------------------------------------------
template <class T, class Index>
baseline_reference_factorize_result<T, Index>
baseline_sparse_gp_lu_factorize(
    const csc_storage<T, Index>& A_csc,
    Index n,
    const std::vector<Index>& col_perm,
    const std::vector<Index>& inv_col_perm,
    const sparse_lu_options<T>& opt)
{
    typedef typename sparse_lu_scalar_policy<T>::real_type real_type;

    baseline_reference_factorize_result<T, Index> result;
    const std::size_t un = static_cast<std::size_t>(n);

    // Trivial 0×0
    if (n == Index(0)) {
        baseline_lu_storage<T, Index>& lu = result.storage;
        lu.L.col_ptr.resize(1u, Index(0));
        lu.U.col_ptr.resize(1u, Index(0));
        lu.row_perm     = col_perm;
        lu.inv_row_perm = inv_col_perm;
        lu.col_perm     = col_perm;
        lu.inv_col_perm = inv_col_perm;
        result.success = true;
        result.status  = sparse_lu_status::success;
        return result;
    }

    // ------------------------------------------------------------------
    // SPA workspace
    // ------------------------------------------------------------------
    std::vector<T>     x(un, T(0));
    std::vector<Index> mark(un, Index(-1));       // mark[i]==j: row i active in col j
    std::vector<Index> reach_mark(un, Index(-1)); // reach_mark[k]==j: k in reach for col j
    std::vector<Index> pattern;
    pattern.reserve(un);

    // Row permutation (pivot-position space):
    //   row_perm[new_row]    = old_row
    //   inv_row_perm[old_row] = new_row
    std::vector<Index> row_perm     = sparse_lu_identity_permutation(n);
    std::vector<Index> inv_row_perm = sparse_lu_identity_permutation(n);

    // L and U built column by column
    std::vector<Index> L_col_ptr(un + 1u, Index(0));
    std::vector<Index> U_col_ptr(un + 1u, Index(0));
    std::vector<Index> L_row_ind_buf, U_row_ind_buf;
    std::vector<T>     L_val_buf,     U_val_buf;

    const real_type threshold        = opt.pivot_threshold;
    const real_type abs_tol          = opt.absolute_pivot_tolerance;
    const bool      do_partial_pivot =
        (opt.pivoting == sparse_lu_pivoting::threshold_partial);

    // SLU-NQ2: reverse index  row -> stored L entries currently carrying that
    // row index.  Maintained only when do_partial_pivot (F1: swaps cannot
    // occur otherwise).  Each per-row list is sorted by col automatically
    // (entries are appended while j increases; merge rebuild preserves it).
    struct l_row_rec { Index col; Index pos; };   // pos: offset into L_*_buf
    std::vector<std::vector<l_row_rec> > l_rows;  // size n when active
    if (do_partial_pivot) {
        l_rows.assign(static_cast<std::size_t>(n), std::vector<l_row_rec>());
    }

    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        pattern.clear();

        // ------------------------------------------------------------------
        // Step 1: Scatter A[:,j] into SPA
        //   orig_row → pivot_row via inv_row_perm
        // ------------------------------------------------------------------
        for (Index p = A_csc.col_ptr[sj]; p < A_csc.col_ptr[sj + 1u]; ++p) {
            const std::size_t sp       = static_cast<std::size_t>(p);
            const Index       orig_row = A_csc.row_ind[sp];
            const Index       piv_row  = inv_row_perm[static_cast<std::size_t>(orig_row)];
            const std::size_t spiv     = static_cast<std::size_t>(piv_row);
            x[spiv] += A_csc.values[sp];
            if (mark[spiv] != j) {
                mark[spiv] = j;
                pattern.push_back(piv_row);
            }
        }

        // ------------------------------------------------------------------
        // Step 2: SLU-6 DFS reach-based column update
        //   Compute reach set via L-structure DFS from pattern seeds (rows < j),
        //   then apply L[:,k] for each active k in reach in increasing order.
        //
        //   Reach note: the DFS may include structural k not numerically active
        //   (x[k]==0).  The guard (mark[sk] != j) below skips those safely.
        //   All k that do become active during this update are guaranteed to be
        //   in reach (proven by transitive closure of the L-column DFS).
        //
        //   Row-pivoting correctness: reach and L structure are in pivot-space
        //   coordinates (inv_row_perm applied at scatter, L rows renamed at
        //   swap steps); no additional coordinate mapping is required here.
        // ------------------------------------------------------------------
        {
            const std::vector<Index> reach =
                sparse_lu_column_reach_from_L(n, j, L_col_ptr, L_row_ind_buf,
                                              pattern, reach_mark);

            for (std::size_t rr = 0u; rr < reach.size(); ++rr) {
                const Index       k   = reach[rr];
                const std::size_t sk  = static_cast<std::size_t>(k);
                if (mark[sk] != j) continue;  // not active in SPA; skip

                const T u_kj = x[sk];  // U[k,j]: reduced value at pivot row k

                for (Index lp = L_col_ptr[sk]; lp < L_col_ptr[sk + 1u]; ++lp) {
                    const std::size_t slp = static_cast<std::size_t>(lp);
                    const Index       i   = L_row_ind_buf[slp];
                    const std::size_t si  = static_cast<std::size_t>(i);
                    x[si] -= L_val_buf[slp] * u_kj;
                    if (mark[si] != j) {
                        mark[si] = j;
                        pattern.push_back(i);
                    }
                }
            }
        }

        // ------------------------------------------------------------------
        // Step 3: Pivot selection
        // ------------------------------------------------------------------
        Index     pivot_pos = j;
        real_type col_max   = real_type(0);
        T         pivot_val = T(0);

        if (do_partial_pivot) {
            for (std::size_t pi = 0u; pi < pattern.size(); ++pi) {
                const Index i = pattern[pi];
                if (i < j) continue;
                const real_type ai =
                    sparse_lu_scalar_policy<T>::abs_value(x[static_cast<std::size_t>(i)]);
                if (ai > col_max) {
                    col_max   = ai;
                    pivot_pos = i;
                    pivot_val = x[static_cast<std::size_t>(i)];
                }
            }
        } else {
            // none / diagonal: use position j directly
            pivot_val = x[sj];
            col_max   = sparse_lu_scalar_policy<T>::abs_value(x[sj]);
        }

        // ------------------------------------------------------------------
        // Step 4: Pivot acceptability check (SLU-3 policy)
        //   pivot is the argmax, so threshold condition is trivially satisfied;
        //   only the absolute_tol rejection can fire here.
        // ------------------------------------------------------------------
        if (!sparse_lu_is_acceptable_pivot(pivot_val, pivot_val, threshold, abs_tol)) {
            result.success = false;
            result.status  = sparse_lu_status::numerical_singularity;
            for (std::size_t pi = 0u; pi < pattern.size(); ++pi)
                x[static_cast<std::size_t>(pattern[pi])] = T(0);
            return result;
        }

        // ------------------------------------------------------------------
        // Step 5: Row swap (bring pivot to position j)
        // ------------------------------------------------------------------
        if (pivot_pos != j) {
            const std::size_t sp = static_cast<std::size_t>(pivot_pos);
            std::swap(x[sj], x[sp]);

            // Ensure j is marked so it gets stored in U
            if (mark[sj] != j) {
                mark[sj] = j;
                pattern.push_back(j);
            }

            // Update row permutation
            const Index r_j = row_perm[sj];
            const Index r_p = row_perm[sp];
            std::swap(row_perm[sj], row_perm[sp]);
            inv_row_perm[static_cast<std::size_t>(r_j)] = pivot_pos;
            inv_row_perm[static_cast<std::size_t>(r_p)] = j;

            // SLU-NQ2: rewire rows j / pivot_pos using the reverse index
            // instead of scanning every stored L column.  Two-pointer merge
            // over the (col-sorted) lists reproduces the old per-column
            // three-case semantics exactly:
            //   col in both lists  -> swap VALUES (row indices unchanged;
            //                          list membership unchanged)
            //   col only in row j  -> rename index j -> pivot_pos; record
            //                          migrates to pivot_pos's list
            //   col only in pivot  -> rename pivot_pos -> j; record migrates
            // Written values and indices are identical to the old full scan,
            // so the final L buffers remain byte-identical.
            {
                std::vector<l_row_rec>& lj = l_rows[sj];
                std::vector<l_row_rec>& lp = l_rows[sp];
                std::vector<l_row_rec> nj, np;
                nj.reserve(lj.size() + lp.size());
                np.reserve(lj.size() + lp.size());
                std::size_t a = 0u, b = 0u;
                while (a < lj.size() || b < lp.size()) {
                    if (b >= lp.size() ||
                        (a < lj.size() && lj[a].col < lp[b].col)) {
                        // col only in row j: rename j -> pivot_pos
                        L_row_ind_buf[static_cast<std::size_t>(lj[a].pos)] = pivot_pos;
                        np.push_back(lj[a]); ++a;
                    } else if (a >= lj.size() || lp[b].col < lj[a].col) {
                        // col only in pivot row: rename pivot_pos -> j
                        L_row_ind_buf[static_cast<std::size_t>(lp[b].pos)] = j;
                        nj.push_back(lp[b]); ++b;
                    } else {
                        // both rows in this column: swap values only
                        std::swap(L_val_buf[static_cast<std::size_t>(lj[a].pos)],
                                  L_val_buf[static_cast<std::size_t>(lp[b].pos)]);
                        nj.push_back(lj[a]); ++a;
                        np.push_back(lp[b]); ++b;
                    }
                }
                lj.swap(nj);
                lp.swap(np);
            }
        } else {
            if (mark[sj] != j) {
                mark[sj] = j;
                pattern.push_back(j);
            }
        }

        const T pivot_diag = x[sj];

        // ------------------------------------------------------------------
        // Steps 6+7 (SLU-NQ1): Store U column j (rows 0..j) and L column j
        //   (rows j+1..n-1, normalized by pivot) by a single ascending sweep
        //   of the sorted SPA pattern (pattern is built in sync with mark and
        //   holds no duplicates -- construction invariant of Steps 1/2/5).
        //   The previous implementation scanned the full row range
        //   (i=0..j then i=j+1..n-1), costing n iterations per column and
        //   n^2 total; replaced by SLU-NQ1
        //   (issue_SLU_numeric_store_quadratic.md, plan A).
        //   Store order (ascending row indices), store conditions, and stored
        //   values are identical to the old scans, so the resulting L/U
        //   storage is byte-identical.
        //   - i <  j: U entry.  The mark guard is always true for pattern
        //     members (invariant above); kept defensively.  Zero strict-upper
        //     entries skipped, as before.
        //   - i == j: diagonal, stored unconditionally (same semantics as the
        //     old `if (i == j)` branch -- no is_exact_zero skip).  j is
        //     guaranteed to be in pattern by both branches of Step 5.
        //   - i >  j: L entry; zero test applies to the value AFTER division
        //     by the pivot (same semantics as the old Step 7).
        // ------------------------------------------------------------------
        std::sort(pattern.begin(), pattern.end());
        Index u_cnt = Index(0);
        Index l_cnt = Index(0);
        for (std::size_t pi = 0u; pi < pattern.size(); ++pi) {
            const Index       i  = pattern[pi];
            const std::size_t si = static_cast<std::size_t>(i);
            if (i < j) {
                if (mark[si] == j &&
                    !sparse_lu_scalar_policy<T>::is_exact_zero(x[si])) {
                    U_row_ind_buf.push_back(i);
                    U_val_buf.push_back(x[si]);
                    ++u_cnt;
                }
            } else if (i == j) {
                U_row_ind_buf.push_back(i);
                U_val_buf.push_back(x[si]);
                ++u_cnt;
            } else if (mark[si] == j) {
                const T l_val = x[si] / pivot_diag;
                if (!sparse_lu_scalar_policy<T>::is_exact_zero(l_val)) {
                    L_row_ind_buf.push_back(i);
                    L_val_buf.push_back(l_val);
                    // SLU-NQ2: O(1) reverse-index record per stored L entry.
                    if (do_partial_pivot) {
                        l_rows[si].push_back(l_row_rec{
                            j, static_cast<Index>(L_row_ind_buf.size() - 1u) });
                    }
                    ++l_cnt;
                }
            }
        }
        U_col_ptr[sj + 1u] = U_col_ptr[sj] + u_cnt;
        L_col_ptr[sj + 1u] = L_col_ptr[sj] + l_cnt;

        // ------------------------------------------------------------------
        // Step 8: Clear SPA for this column
        // ------------------------------------------------------------------
        for (std::size_t pi = 0u; pi < pattern.size(); ++pi)
            x[static_cast<std::size_t>(pattern[pi])] = T(0);
        x[sj] = T(0);  // safety: clear diagonal even if already in pattern
    }

    // Assemble result
    baseline_lu_storage<T, Index>& lu = result.storage;
    lu.L.col_ptr  = L_col_ptr;
    lu.L.row_ind  = L_row_ind_buf;
    lu.L.values   = L_val_buf;
    lu.U.col_ptr  = U_col_ptr;
    lu.U.row_ind  = U_row_ind_buf;
    lu.U.values   = U_val_buf;
    lu.row_perm     = row_perm;
    lu.inv_row_perm = sparse_lu_inverse_permutation(row_perm);
    lu.col_perm     = col_perm;
    lu.inv_col_perm = inv_col_perm;

    result.success = true;
    result.status  = sparse_lu_status::success;
    return result;
}

// ---------------------------------------------------------------------------
// baseline_reference_lu_factorize  (SLU-4.1 dense-reference: TEST-ONLY)
//
// Correctness-first dense-reference baseline numeric factorization retained as
// a test/diagnostic reference path.  It is NOT the production default as of
// SLU-5; the production path uses baseline_sparse_gp_lu_factorize above.
//
// Uses an internal dense workspace (W[row][col] = std::vector<std::vector<T>>)
// for robustness and deterministic testing, then emits owning CSC L/U baseline
// storage consumed by the SLU-2 solve path.
// ---------------------------------------------------------------------------
template <class T, class Index>
baseline_reference_factorize_result<T, Index>
baseline_reference_lu_factorize(
    const csc_storage<T, Index>& A_csc,
    Index n,
    const std::vector<Index>& col_perm,
    const std::vector<Index>& inv_col_perm,
    const sparse_lu_options<T>& opt)
{
    typedef typename sparse_lu_scalar_policy<T>::real_type real_type;

    baseline_reference_factorize_result<T, Index> result;

    const std::size_t un = static_cast<std::size_t>(n);

    if (n == Index(0)) {
        baseline_lu_storage<T, Index>& lu = result.storage;
        lu.L.col_ptr.resize(1u, Index(0));
        lu.U.col_ptr.resize(1u, Index(0));
        lu.row_perm     = col_perm;
        lu.inv_row_perm = inv_col_perm;
        lu.col_perm     = col_perm;
        lu.inv_col_perm = inv_col_perm;
        result.success = true;
        result.status  = sparse_lu_status::success;
        return result;
    }

    // Dense workspace W[row][col]  (test-only; production uses sparse GP above)
    std::vector<std::vector<T> > W(un, std::vector<T>(un, T(0)));
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index kp = A_csc.col_ptr[sj]; kp < A_csc.col_ptr[sj + 1u]; ++kp) {
            const std::size_t skp = static_cast<std::size_t>(kp);
            const std::size_t sr  = static_cast<std::size_t>(A_csc.row_ind[skp]);
            W[sr][sj] = A_csc.values[skp];
        }
    }

    std::vector<Index> row_perm = sparse_lu_identity_permutation(n);

    const real_type threshold = opt.pivot_threshold;
    const real_type abs_tol   = opt.absolute_pivot_tolerance;
    const bool do_partial_pivot =
        (opt.pivoting == sparse_lu_pivoting::threshold_partial);

    for (Index k = Index(0); k < n; ++k) {
        const std::size_t sk = static_cast<std::size_t>(k);

        Index pivot_row  = k;
        T     col_max_val = W[sk][sk];

        if (do_partial_pivot) {
            real_type col_max_abs =
                sparse_lu_scalar_policy<T>::abs_value(W[sk][sk]);
            for (Index i = k + Index(1); i < n; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                const real_type ai =
                    sparse_lu_scalar_policy<T>::abs_value(W[si][sk]);
                if (ai > col_max_abs) {
                    col_max_abs = ai;
                    col_max_val = W[si][sk];
                    pivot_row   = i;
                }
            }
        }

        const T pivot_val = col_max_val;

        if (!sparse_lu_is_acceptable_pivot(pivot_val, col_max_val,
                                           threshold, abs_tol)) {
            result.success = false;
            result.status  = sparse_lu_status::numerical_singularity;
            return result;
        }

        if (pivot_row != k) {
            W[sk].swap(W[static_cast<std::size_t>(pivot_row)]);
            std::swap(row_perm[sk],
                      row_perm[static_cast<std::size_t>(pivot_row)]);
        }

        const T pivot_diag = W[sk][sk];
        for (Index i = k + Index(1); i < n; ++i) {
            const std::size_t si = static_cast<std::size_t>(i);
            W[si][sk] /= pivot_diag;
        }

        for (Index i = k + Index(1); i < n; ++i) {
            const std::size_t si = static_cast<std::size_t>(i);
            const T lik = W[si][sk];
            if (sparse_lu_scalar_policy<T>::is_exact_zero(lik)) continue;
            for (Index j = k + Index(1); j < n; ++j) {
                const std::size_t sj = static_cast<std::size_t>(j);
                W[si][sj] -= lik * W[sk][sj];
            }
        }
    }

    baseline_lu_storage<T, Index>& lu = result.storage;
    csc_storage<T, Index>& L = lu.L;

    L.col_ptr.resize(un + 1u, Index(0));
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        Index cnt = Index(0);
        for (Index i = j + Index(1); i < n; ++i) {
            if (!sparse_lu_scalar_policy<T>::is_exact_zero(
                    W[static_cast<std::size_t>(i)][sj])) {
                ++cnt;
            }
        }
        L.col_ptr[sj + 1u] = cnt;
    }
    for (std::size_t j = 0u; j < un; ++j) {
        L.col_ptr[j + 1u] += L.col_ptr[j];
    }
    {
        const Index nnz_L = L.col_ptr[un];
        L.row_ind.resize(static_cast<std::size_t>(nnz_L));
        L.values.resize(static_cast<std::size_t>(nnz_L));
        for (Index j = Index(0); j < n; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            Index out = L.col_ptr[sj];
            for (Index i = j + Index(1); i < n; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                if (!sparse_lu_scalar_policy<T>::is_exact_zero(W[si][sj])) {
                    L.row_ind[static_cast<std::size_t>(out)] = i;
                    L.values[static_cast<std::size_t>(out)]  = W[si][sj];
                    ++out;
                }
            }
        }
    }

    csc_storage<T, Index>& U = lu.U;

    U.col_ptr.resize(un + 1u, Index(0));
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        Index cnt = Index(0);
        for (Index i = Index(0); i <= j; ++i) {
            const std::size_t si = static_cast<std::size_t>(i);
            if (i == j || !sparse_lu_scalar_policy<T>::is_exact_zero(W[si][sj])) {
                ++cnt;
            }
        }
        U.col_ptr[sj + 1u] = cnt;
    }
    for (std::size_t j = 0u; j < un; ++j) {
        U.col_ptr[j + 1u] += U.col_ptr[j];
    }
    {
        const Index nnz_U = U.col_ptr[un];
        U.row_ind.resize(static_cast<std::size_t>(nnz_U));
        U.values.resize(static_cast<std::size_t>(nnz_U));
        for (Index j = Index(0); j < n; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            Index out = U.col_ptr[sj];
            for (Index i = Index(0); i <= j; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                if (i == j ||
                    !sparse_lu_scalar_policy<T>::is_exact_zero(W[si][sj])) {
                    U.row_ind[static_cast<std::size_t>(out)] = i;
                    U.values[static_cast<std::size_t>(out)]  = W[si][sj];
                    ++out;
                }
            }
        }
    }

    lu.row_perm     = row_perm;
    lu.inv_row_perm = sparse_lu_inverse_permutation(row_perm);
    lu.col_perm     = col_perm;
    lu.inv_col_perm = inv_col_perm;

    result.success = true;
    result.status  = sparse_lu_status::success;
    return result;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_NUMERIC_IMPL_HPP
