// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-8R.4: §17.2(B) Within-panel factorization -- internal implementation.
//
// Implements §17.2(B): within-panel factorization with active panel-height pivoting.
//
// Algorithm:
//   For each supernode j:
//     1. Gather panel workspace from panel_values (post §17.2(A) update).
//     2. For each pivot column c in [0, w_j):
//        a. Search ALL active rows (c..R-1) for pivot candidate (active panel height).
//        b. Threshold partial pivoting decision.
//        c. Row swap if needed (updates workspace, row_ind, storage.row_perm).
//        d. Check zero / near-zero pivot status.
//        e. Scale column below pivot (L multipliers).
//        f. Update trailing columns via dense_kernel::ger (rank-1 update).
//     3. Scatter factored workspace back to panel_values.
//     4. Update desc_j.row_indices to reflect within-panel row ordering.
//
// SCOPE:
//   This file implements §17.2(B) ONLY.
//   §18.2 storage-native supernodal solve is NOT implemented here.
//   true_supernodal_numeric remains false (CSC bootstrap is still source of truth).
//   Gate 6 remains PENDING.
//   supernodal_solve_native remains false.
//   SLU-8 full conformance is NOT claimed.
//
// BOUNDARY:
//   within_panel_factorization_is_numeric_source = false:
//     The within-panel factorization is applied to transitional storage data
//     (bootstrapped from CSC L/U + §17.2(A) update). It is NOT the final
//     numeric source of truth for solve operations.
//   Solve continues to use the baseline CSC L/U path (Gate 4 maintained).
//   storage.row_perm / inv_row_perm are updated to reflect within-panel row swaps,
//   but these updates do not affect solve (which uses baseline_.row_perm).
//
// DENSE KERNEL ADAPTER USAGE:
//   Trailing column update uses sparse_lu_dense_kernel<T>::ger.
//   Direct vcp::tger calls are prohibited here.
//   dense_kernel::getrf is NOT used (§25 prohibition, §17.2(B) pivot search).
//   within_panel_used_getrf == false in normal path.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - supernodal_lu_storage<T, Index> is defined (with panel_values, U_segments)
//   - supernode_desc<Index> is defined (with row_indices, u_seg_col_ptr)
//   - sparse_lu_dense_kernel<T> adapter is defined and implemented
//   - sparse_lu_scalar_policy<T> is defined
//   - sparse_lu_options<T> is defined
//   - <chrono>, <algorithm>, <vector> are available
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//
// SLU-8 Gate status after SLU-8R.4:
//   Gate 1: PASS (dense kernel adapter connected)
//   Gate 2: PARTIAL / PENDING-REQUIRED (§17.2(A)+(B) calls exist; full PASS requires
//            these to be main factorization component, not transitional calls)
//   Gate 3: PENDING-REQUIRED
//   Gate 4: PASS (solve uses CSC baseline; panel_values not used in solve)
//   Gate 5: PENDING-REQUIRED
//   Gate 6: PENDING (baseline GP still numeric source of truth)
//
// §18.2 storage-native supernodal solve: PENDING.
// issue_SLU8_contract_violation.md: OPEN.

#ifndef VCP_TSPARSE_SPARSE_LU_WITHIN_PANEL_FACTOR_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_WITHIN_PANEL_FACTOR_IMPL_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <vector>

// within_panel_factor_stats is defined in tsparse_sparse_lu.hpp
// (in namespace sparse_lu_detail, before sparse_lu_factorization class body).
// Do NOT redefine it here.

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// within_panel_workspace<T, Index>
//
// Dense workspace for the current supernode panel during §17.2(B) factorization.
//
// Layout:
//   row_ind:     current global row ids for local rows 0..rows-1.
//                Starts as a copy of desc_j.row_indices.
//                Updated in-place when rows are swapped.
//   values:      panel data, column-major, rows x cols, leading dim = ld.
//                values[c * ld + r] = (local row r, panel col c) entry.
//                Starts as a copy of panel_values[values_offset..].
//   local_perm:  local_perm[r] = original local row now at position r.
//                Initially identity. Updated on each row swap.
//   inv_local_perm: inv_local_perm[orig_r] = current position of original row orig_r.
//                Initially identity. Maintained as inverse of local_perm.
//   rows:        total panel rows (= desc_j.row_indices.size()).
//   cols:        panel width (= desc_j.num_cols).
//   ld:          leading dimension of values (= desc_j.leading_dimension >= rows).
//   valid:       true iff workspace was successfully gathered.
//
// Invariant: inv_local_perm[local_perm[r]] == r for all r in [0, rows).
// ---------------------------------------------------------------------------
template <class T, class Index>
struct within_panel_workspace {
    static_assert(std::is_signed<Index>::value,
                  "within_panel_workspace: Index must be signed");

    std::vector<Index> row_ind;
    std::vector<T>     values;
    std::vector<Index> local_perm;
    std::vector<Index> inv_local_perm;
    Index rows;
    Index cols;
    Index ld;
    bool  valid;

    within_panel_workspace()
        : rows(Index(0)), cols(Index(0)), ld(Index(0)), valid(false) {}
};

// ---------------------------------------------------------------------------
// gather_within_panel_workspace<T, Index>
//
// Gathers the panel data for supernode j from panel_values into a dense
// column-major workspace covering all panel rows (diagonal block + off-diagonal L).
//
// Initializes identity local permutation (no swaps yet).
// ---------------------------------------------------------------------------
template <class T, class Index>
void gather_within_panel_workspace(
    const supernodal_lu_storage<T, Index>& storage,
    Index supernode_j,
    within_panel_workspace<T, Index>& work)
{
    work.valid = false;
    work.rows  = Index(0);
    work.cols  = Index(0);
    work.ld    = Index(0);
    work.row_ind.clear();
    work.values.clear();
    work.local_perm.clear();
    work.inv_local_perm.clear();

    if (!storage.valid) return;
    const std::size_t j = static_cast<std::size_t>(supernode_j);
    if (j >= storage.supernodes.size()) return;

    const supernode_desc<Index>& desc_j = storage.supernodes[j];
    const Index w   = desc_j.num_cols;
    const Index R   = static_cast<Index>(desc_j.row_indices.size());
    const Index ld  = desc_j.leading_dimension;

    if (w <= Index(0) || R <= Index(0) || ld < R) return;

    const std::size_t offset = static_cast<std::size_t>(desc_j.values_offset);
    const std::size_t total  =
        static_cast<std::size_t>(ld) * static_cast<std::size_t>(w);
    if (offset + total > storage.panel_values.size()) return;

    // Copy row indices and values
    work.row_ind.assign(desc_j.row_indices.begin(), desc_j.row_indices.end());
    work.values.assign(
        storage.panel_values.begin() + static_cast<std::ptrdiff_t>(offset),
        storage.panel_values.begin() + static_cast<std::ptrdiff_t>(offset + total));
    work.rows = R;
    work.cols = w;
    work.ld   = ld;

    // Identity local permutation
    const std::size_t R_sz = static_cast<std::size_t>(R);
    work.local_perm.resize(R_sz);
    work.inv_local_perm.resize(R_sz);
    for (std::size_t r = 0u; r < R_sz; ++r) {
        work.local_perm[r]     = static_cast<Index>(r);
        work.inv_local_perm[r] = static_cast<Index>(r);
    }

    work.valid = true;
}

// ---------------------------------------------------------------------------
// select_threshold_pivot<T, Index>
//
// Selects the pivot for column pivot_col in the within-panel workspace.
// Searches ALL active rows (pivot_col .. R-1), i.e. active panel height.
//
// Returns: pivot_decision (acceptable / reject / inconclusive)
//   - acceptable: diagonal entry at (pivot_col, pivot_col) passes threshold test.
//     pivot_row_out = pivot_col (no swap needed).
//   - reject: diagonal rejected; max-abs row is the pivot.
//     pivot_row_out = max-abs row index (>= pivot_col).
//   - inconclusive: non-finite values encountered; max-abs row used as fallback.
//     pivot_row_out = max-abs row index.
//
// Also sets:
//   pivot_value_out: the value at (pivot_col, pivot_row_out) BEFORE any swap.
//   max_abs_out: maximum abs value in the active column.
//
// IMPORTANT: getrf is NOT used here. This is an explicit sparse-LU pivot search
//            over the active panel height (§17.2(B) requirement, §25 prohibition).
//
// Complexity: O(R - pivot_col) per column.
// ---------------------------------------------------------------------------
template <class T, class Index>
pivot_decision select_threshold_pivot(
    const within_panel_workspace<T, Index>& work,
    Index pivot_col,
    const sparse_lu_options<T>& opt,
    Index& pivot_row_out,
    T& pivot_value_out)
{
    typedef sparse_lu_scalar_policy<T> scalar_pol;
    typedef typename scalar_pol::real_type real_type;

    const Index R  = work.rows;
    const Index ld = work.ld;
    const Index c  = pivot_col;

    // Default: diagonal position
    pivot_row_out   = c;
    pivot_value_out = T(0);

    if (c >= R || c >= work.cols) return pivot_decision::reject;

    // Scan active rows c..R-1 for max abs
    real_type max_abs = real_type(0);
    Index     max_row = c;
    bool      has_nonfinite = false;

    for (Index r = c; r < R; ++r) {
        const T& val = work.values[static_cast<std::size_t>(c * ld + r)];
        const real_type av = scalar_pol::abs_value(val);
        if (!vcp::tsparse_scalar::is_finite(av)) {
            has_nonfinite = true;
            continue;
        }
        if (av > max_abs) {
            max_abs = av;
            max_row = r;
        }
    }

    if (has_nonfinite) {
        // Non-finite value in active column: inconclusive
        pivot_row_out   = max_row;
        pivot_value_out = work.values[static_cast<std::size_t>(c * ld + max_row)];
        return pivot_decision::inconclusive;
    }

    if (max_abs == real_type(0)) {
        // All active entries are zero: zero pivot (reject)
        pivot_row_out   = c;
        pivot_value_out = T(0);
        return pivot_decision::reject;
    }

    // Threshold partial pivoting decision on the diagonal entry
    const T& diag_val = work.values[static_cast<std::size_t>(c * ld + c)];
    const real_type diag_abs = scalar_pol::abs_value(diag_val);

    pivot_decision decision;
    if (!vcp::tsparse_scalar::is_finite(diag_abs)) {
        // Diagonal non-finite: use max row
        decision = pivot_decision::inconclusive;
    } else {
        // acceptable if abs(diag) >= threshold * max_abs
        decision = scalar_pol::acceptable_pivot(diag_val, max_abs, opt.pivot_threshold);
    }

    if (decision == pivot_decision::acceptable) {
        pivot_row_out   = c;
        pivot_value_out = diag_val;
    } else {
        // reject or inconclusive: use max-abs row
        pivot_row_out   = max_row;
        pivot_value_out = work.values[static_cast<std::size_t>(c * ld + max_row)];
    }

    return decision;
}

// ---------------------------------------------------------------------------
// apply_within_panel_row_swap<T, Index>
//
// Swaps local rows r1 and r2 in the panel workspace.
// Updates:
//   - work.values: all column entries swapped
//   - work.row_ind: global row ids swapped
//   - work.local_perm / inv_local_perm: permutation tracking
//   - storage.row_perm / inv_row_perm: global permutation consistency
//
// Maintains invariant: storage.inv_row_perm[storage.row_perm[i]] == i.
// Does nothing if r1 == r2.
// ---------------------------------------------------------------------------
template <class T, class Index>
void apply_within_panel_row_swap(
    within_panel_workspace<T, Index>& work,
    Index r1, Index r2,
    supernodal_lu_storage<T, Index>& storage,
    within_panel_factor_stats& stats)
{
    if (r1 == r2) return;

    const Index ld = work.ld;
    const Index w  = work.cols;
    const std::size_t sr1 = static_cast<std::size_t>(r1);
    const std::size_t sr2 = static_cast<std::size_t>(r2);

    // Swap values in all columns
    for (Index c = Index(0); c < w; ++c) {
        std::swap(
            work.values[static_cast<std::size_t>(c * ld) + sr1],
            work.values[static_cast<std::size_t>(c * ld) + sr2]);
    }

    // Swap global row ids and update local permutation
    const Index g_r1 = work.row_ind[sr1];
    const Index g_r2 = work.row_ind[sr2];
    std::swap(work.row_ind[sr1], work.row_ind[sr2]);

    // Update local permutation
    std::swap(work.local_perm[sr1], work.local_perm[sr2]);
    // Repair inverse permutation
    work.inv_local_perm[static_cast<std::size_t>(work.local_perm[sr1])] = r1;
    work.inv_local_perm[static_cast<std::size_t>(work.local_perm[sr2])] = r2;

    // Update storage.row_perm / inv_row_perm for the two global rows swapped.
    // storage.row_perm[global_r] = original_r
    // storage.inv_row_perm[original_r] = global_r
    const Index n_perm     = static_cast<Index>(storage.row_perm.size());
    const Index n_inv_perm = static_cast<Index>(storage.inv_row_perm.size());
    if (g_r1 >= Index(0) && g_r1 < n_perm &&
        g_r2 >= Index(0) && g_r2 < n_perm) {
        const Index old_orig_r1 = storage.row_perm[static_cast<std::size_t>(g_r1)];
        const Index old_orig_r2 = storage.row_perm[static_cast<std::size_t>(g_r2)];
        std::swap(storage.row_perm[static_cast<std::size_t>(g_r1)],
                  storage.row_perm[static_cast<std::size_t>(g_r2)]);
        if (old_orig_r1 >= Index(0) && old_orig_r1 < n_inv_perm)
            storage.inv_row_perm[static_cast<std::size_t>(old_orig_r1)] = g_r2;
        if (old_orig_r2 >= Index(0) && old_orig_r2 < n_inv_perm)
            storage.inv_row_perm[static_cast<std::size_t>(old_orig_r2)] = g_r1;
    }

    stats.row_swap_count++;
}

// ---------------------------------------------------------------------------
// scale_pivot_column<T, Index>
//
// Scales rows below the pivot in column pivot_col by 1/pivot_value.
// Assumes pivot is at (pivot_col, pivot_col) (after any row swap).
// Increments stats.scale_count (one per column where scaling is applied).
//
// Precondition: abs(pivot_value) > 0 must be CERTIFIED by the caller via the
// Step 4 gate !(abs_pv > opt.zero_tolerance) -> skip (SLU-GT1 D3).  For
// interval scalars this means the pivot magnitude is certainly positive (an
// interval containing 0 never reaches this division).
// ---------------------------------------------------------------------------
template <class T, class Index>
void scale_pivot_column(
    within_panel_workspace<T, Index>& work,
    Index pivot_col,
    const T& pivot_value,
    within_panel_factor_stats& stats)
{
    const Index R  = work.rows;
    const Index ld = work.ld;
    const Index c  = pivot_col;
    const std::size_t sc = static_cast<std::size_t>(c);

    // Scale L multipliers: rows c+1..R-1 in column c
    bool scaled = false;
    for (Index r = c + Index(1); r < R; ++r) {
        work.values[sc * static_cast<std::size_t>(ld) + static_cast<std::size_t>(r)]
            /= pivot_value;
        scaled = true;
    }
    if (scaled) {
        stats.scale_count++;
    }
}

// ---------------------------------------------------------------------------
// update_trailing_columns<T, Index>
//
// Applies rank-1 update (Schur complement) to trailing columns c+1..w-1.
// Uses sparse_lu_dense_kernel<T>::ger via the adapter.
// Direct vcp::tger calls are prohibited here.
//
// Update: panel[c+1:R, c+1:w] -= L_col * U_row^T
//   where L_col = panel[c+1:R, c]  (L multipliers, already scaled)
//         U_row = panel[c, c+1:w]  (U row entries, non-contiguous in column-major)
//
// y_copy is allocated internally to provide contiguous U row for ger.
// Increments stats.ger_count on successful call.
// Timing accumulated into stats.dense_kernel_ticks.
//
// Postcondition: panel[r, q] -= L[r,c] * U[c,q] for all r > c, q > c.
// ---------------------------------------------------------------------------
template <class T, class Index>
void update_trailing_columns(
    within_panel_workspace<T, Index>& work,
    Index pivot_col,
    within_panel_factor_stats& stats)
{
    const Index R  = work.rows;
    const Index ld = work.ld;
    const Index w  = work.cols;
    const Index c  = pivot_col;

    const Index m_idx = R - c - Index(1);   // rows below pivot
    const Index n_idx = w - c - Index(1);   // trailing columns

    if (m_idx <= Index(0) || n_idx <= Index(0)) return;

    const std::size_t m = static_cast<std::size_t>(m_idx);
    const std::size_t n = static_cast<std::size_t>(n_idx);
    const std::size_t sc = static_cast<std::size_t>(c);
    const std::size_t sld = static_cast<std::size_t>(ld);

    // x = L column c, rows c+1..R-1 (contiguous in column-major)
    const T* x = &work.values[sc * sld + sc + 1u];

    // y = U row c in cols c+1..w-1 (NOT contiguous; stride = ld)
    // Copy to contiguous temp buffer
    std::vector<T> y_copy(n);
    for (std::size_t q = 0u; q < n; ++q) {
        y_copy[q] = work.values[(sc + 1u + q) * sld + sc];
    }

    // A = trailing submatrix at (c+1, c+1), leading dim ld
    T* A = &work.values[(sc + 1u) * sld + sc + 1u];

    // Rank-1 update: A -= x * y^T
    // This goes through the adapter (no direct vcp::tger call).
    auto t0 = std::chrono::steady_clock::now();
    sparse_lu_dense_kernel<T>::ger(m, n, T(-1), x, y_copy.data(), A, sld);
    auto t1 = std::chrono::steady_clock::now();

    stats.dense_kernel_ticks += static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    stats.ger_count++;
    // SLU-PERF: ger(m,n) FLOP = 2*m*n (rank-1 trailing update, BLAS-2).
    stats.flop_ger += 2.0 * static_cast<double>(m) * static_cast<double>(n);
}

// ---------------------------------------------------------------------------
// scatter_within_panel_factor<T, Index>
//
// Writes the factored workspace back to storage.panel_values for supernode j.
// Also updates desc_j.row_indices to reflect any within-panel row swaps
// (so that the supernodal storage reflects the within-panel pivot sequence).
// ---------------------------------------------------------------------------
template <class T, class Index>
void scatter_within_panel_factor(
    supernodal_lu_storage<T, Index>& storage,
    Index supernode_j,
    const within_panel_workspace<T, Index>& work,
    within_panel_factor_stats& stats)
{
    (void)stats;
    if (!work.valid) return;

    const std::size_t j = static_cast<std::size_t>(supernode_j);
    if (j >= storage.supernodes.size()) return;

    supernode_desc<Index>& desc_j = storage.supernodes[j];
    const Index ld = desc_j.leading_dimension;
    const Index w  = desc_j.num_cols;
    const std::size_t offset = static_cast<std::size_t>(desc_j.values_offset);
    const std::size_t total  =
        static_cast<std::size_t>(ld) * static_cast<std::size_t>(w);

    if (offset + total > storage.panel_values.size()) return;
    if (work.values.size() < total) return;

    // Write back factored values
    for (std::size_t i = 0u; i < total; ++i) {
        storage.panel_values[offset + i] = work.values[i];
    }

    // Update row_indices to reflect within-panel row ordering after swaps
    const std::size_t R_sz = static_cast<std::size_t>(work.rows);
    if (desc_j.row_indices.size() == R_sz && work.row_ind.size() == R_sz) {
        for (std::size_t r = 0u; r < R_sz; ++r) {
            desc_j.row_indices[r] = work.row_ind[r];
        }
    }
}

// ---------------------------------------------------------------------------
// run_within_panel_factorization<T, Index>
//
// Main driver: applies §17.2(B) within-panel factorization to all supernodes
// in the given supernodal_lu_storage.
//
// For each supernode j with width w and total panel rows R:
//   - Gathers panel workspace from panel_values (post §17.2(A) update)
//   - Performs w-step LU factorization with active panel-height pivot search
//   - Each step: pivot search (rows c..R-1), threshold decision, row swap,
//                scale, rank-1 ger update
//   - Scatters factored workspace back to panel_values
//   - Updates desc_j.row_indices and storage.row_perm/inv_row_perm
//
// PROHIBITIONS (§25):
//   - getrf is NOT called as a pivot-search substitute
//   - within_panel_used_getrf == false in normal path
//   - Direct vcp::tger / vcp::tgemm calls are prohibited (use adapter only)
//
// BOUNDARY:
//   - within_panel_factorization_is_numeric_source = false (CSC source of truth)
//   - true_supernodal_numeric remains false (set in caller via set_within_panel_factor_info_)
//   - Gate 6 remains PENDING (baseline CSC L/U is still numeric source of truth)
//   - Solve uses baseline CSC L/U (Gate 4 maintained; panel_values not used in solve)
//
// Returns: within_panel_factor_stats with evidence of §17.2(B) execution.
// ---------------------------------------------------------------------------
template <class T, class Index>
within_panel_factor_stats
run_within_panel_factorization(
    supernodal_lu_storage<T, Index>& storage,
    const sparse_lu_options<T>& opt)
{
    typedef sparse_lu_scalar_policy<T> scalar_pol;
    typedef typename scalar_pol::real_type real_type;

    within_panel_factor_stats stats;
    stats.executed   = true;
    stats.used_getrf = false;  // §17.2(B) does NOT use getrf (§25)

    if (!storage.valid) {
        stats.completed = false;
        return stats;
    }

    const std::size_t nsup = storage.supernodes.size();
    bool all_completed = true;

    for (std::size_t j = 0u; j < nsup; ++j) {
        const supernode_desc<Index>& desc_j_const = storage.supernodes[j];
        const Index w = desc_j_const.num_cols;
        const Index R = static_cast<Index>(desc_j_const.row_indices.size());

        if (w <= Index(0) || R <= Index(0)) continue;

        stats.panel_count++;

        // Gather panel workspace (copies panel_values + row_indices for this supernode)
        within_panel_workspace<T, Index> work;
        gather_within_panel_workspace(storage, static_cast<Index>(j), work);

        if (!work.valid) {
            all_completed = false;
            continue;
        }

        // Within-panel factorization loop: one step per pivot column
        bool panel_completed = true;
        for (Index c = Index(0); c < w; ++c) {
            // ----------------------------------------------------------------
            // Step 1: Active panel-height pivot search (rows c..R-1).
            //         NOT delegated to getrf. This is the explicit Sparse LU
            //         pivot search as required by §17.2(B).
            // ----------------------------------------------------------------
            stats.pivot_search_count++;

            Index pivot_row;
            T     pivot_value;
            pivot_decision decision = select_threshold_pivot(
                work, c, opt, pivot_row, pivot_value);

            // ----------------------------------------------------------------
            // Step 2: Count pivot decisions
            // ----------------------------------------------------------------
            if (decision == pivot_decision::acceptable) {
                stats.pivot_accept_count++;
            } else if (decision == pivot_decision::reject) {
                stats.pivot_reject_count++;
            } else {
                // inconclusive
                stats.inconclusive_pivot_count++;
            }

            // ----------------------------------------------------------------
            // Step 3: Row swap if pivot row != current diagonal position
            // ----------------------------------------------------------------
            if (pivot_row != c) {
                apply_within_panel_row_swap(work, c, pivot_row, storage, stats);
            }

            // After potential swap, pivot is at (c, c)
            const T& actual_pivot =
                work.values[static_cast<std::size_t>(c * work.ld + c)];
            const real_type abs_pv = scalar_pol::abs_value(actual_pivot);

            // ----------------------------------------------------------------
            // Step 4: Zero / near-zero pivot status tracking
            //
            // SLU-GT1 D3 (certified-only): !(abs_pv > tol) means "cannot
            // certify abs_pv > tol".  For interval scalars a pivot containing
            // 0 is skipped HERE, so only certifiably nonzero pivots reach the
            // division in scale_pivot_column (Step 5).  For totally ordered
            // scalars this is identical to abs_pv <= tol.
            // ----------------------------------------------------------------
            if (!(abs_pv > opt.zero_tolerance)) {
                stats.zero_pivot_count++;
                // Zero pivot: skip scale and trailing update for this column.
                // Status propagation: not silently consumed.
                panel_completed = false;
                continue;
            }
            if (!(abs_pv > opt.near_zero_tolerance)) {
                stats.near_zero_pivot_count++;
                // Near-zero pivot: warning only, continue factorization.
            }

            // ----------------------------------------------------------------
            // Step 5: Scale L multipliers in column c below pivot
            // ----------------------------------------------------------------
            scale_pivot_column(work, c, actual_pivot, stats);

            // ----------------------------------------------------------------
            // Step 6: Trailing columns update via dense_kernel::ger
            //         A[c+1:R, c+1:w] -= L_col * U_row^T
            // ----------------------------------------------------------------
            update_trailing_columns(work, c, stats);
        }

        if (!panel_completed) all_completed = false;

        // ----------------------------------------------------------------
        // Step 7: Scatter factored workspace back to panel_values.
        //         Update desc_j.row_indices to reflect within-panel row order.
        // ----------------------------------------------------------------
        scatter_within_panel_factor(storage, static_cast<Index>(j), work, stats);
    }

    stats.completed         = all_completed;
    stats.is_numeric_source = false;  // CSC still source of truth in SLU-8R.4

    // SLU-8R.4.1: Determine worst pivot event status.
    // Priority (highest wins, set last): zero_pivot > rejected_pivot > inconclusive_pivot > near_zero_warning > success.
    // Any != success means a non-silent event occurred and must NOT be masked as "success".
    if (stats.pivot_search_count > 0u) {
        stats.worst_status = within_panel_factor_status::success;
        if (stats.near_zero_pivot_count > 0u)
            stats.worst_status = within_panel_factor_status::near_zero_warning;
        if (stats.inconclusive_pivot_count > 0u)
            stats.worst_status = within_panel_factor_status::inconclusive_pivot;
        if (stats.pivot_reject_count > 0u)
            stats.worst_status = within_panel_factor_status::rejected_pivot;
        if (stats.zero_pivot_count > 0u)
            stats.worst_status = within_panel_factor_status::zero_pivot;
    } else {
        stats.worst_status = within_panel_factor_status::not_run;
    }

    return stats;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_WITHIN_PANEL_FACTOR_IMPL_HPP
