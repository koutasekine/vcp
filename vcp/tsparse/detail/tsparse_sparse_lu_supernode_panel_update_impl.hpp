// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-8R.3 / SLU-8R.3.1: Supernode-panel left-looking update -- internal implementation.
//
// Implements §17.2(A): supernode-panel left-looking update.
//
// SLU-8R.3.1 REPAIR: Fixed update set criterion and workspace gather/scatter.
//   - compute_panel_update_set: now uses U_segments rows (rows < col_begin_j)
//     instead of row_indices rows (rows >= col_begin_j). This is the correct
//     criterion: K updates J iff J's U_segments contain a row in K's column range.
//   - gather_panel_workspace: extended to include U_segment rows at the head of
//     the workspace (before panel rows). Uses u_seg_col_ptr for per-column lookup.
//   - scatter_panel_workspace: writes U rows back to U_segments and panel rows
//     back to panel_values.
//   - run_supernode_panel_leftlooking_update: sets stats.update_applied only when
//     actual updates occurred (update_count > 0), and sets driver_called/applied/
//     dense_kernel_called in stats.
//     SLU-8R.7.1: the §17.2(A) "applied" diagnostic is authoritative in
//     sparse_lu_info (info.supernode_panel_update_applied), not in storage.
//
// SCOPE:
//   This file implements §17.2(A) ONLY.
//   §17.2(B) within-panel factorization is NOT implemented here.
//   true_supernodal_numeric remains false until §17.2(A)/(B) are both complete.
//   Gate 6 remains PENDING.
//   supernodal_solve_native remains false.
//   SLU-8 full conformance is NOT claimed.
//
// BOUNDARY:
//   stats.update_applied = true only when update_count > 0 (actual updates);
//   surfaced as info.supernode_panel_update_applied (SLU-8R.7.1 authoritative).
//   supernode_panel_update_is_numeric_source = false (§17.2(B) PENDING).
//   The update is applied in-place to panel_values and U_segments (bootstrapped
//   from CSC L/U), so the result is a §17.2(A)-updated transitional state, NOT
//   final factor values. Values may differ from the correct factored values.
//   Solve continues to use the baseline CSC L/U path (Gate 4 maintained).
//   Modifying U_segments is safe: solve uses baseline_ CSC, not U_segments.
//
// DENSE KERNEL ADAPTER USAGE:
//   All dense operations go through sparse_lu_dense_kernel<T> adapter.
//   Direct vcp::tgemm / vcp::ttrsm / vcp::tgemv calls are prohibited here.
//   getrf is NOT used in §17.2(A) (no pivot search; §25 prohibition).
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - supernodal_lu_storage<T, Index> is defined
//   - sparse_lu_dense_kernel<T> adapter is defined and implemented
//   - <chrono> and <algorithm> are available
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//
// SLU-8 Gate status after SLU-8R.3.1:
//   Gate 1: PASS
//   Gate 2: PARTIAL / PENDING-REQUIRED (real §17.2(A) dense calls exist, but full
//           PASS requires §17.2(A)/(B) as main factorization component)
//   Gate 3: PENDING-REQUIRED
//   Gate 4: PASS (solve uses CSC baseline; panel_values/U_segments not used in solve)
//   Gate 5: PENDING-REQUIRED
//   Gate 6: PENDING (baseline GP still numeric source of truth)
//
// §17.2(B) within-panel factorization: PENDING.
// §18.2 storage-native supernodal solve: PENDING.
// issue_SLU8_contract_violation.md: OPEN.

#ifndef VCP_TSPARSE_SPARSE_LU_SUPERNODE_PANEL_UPDATE_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SUPERNODE_PANEL_UPDATE_IMPL_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <vector>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// supernode_panel_workspace<T, Index>
//
// Dense workspace for the current supernode panel during §17.2(A) update.
//
// SLU-8R.3.1 EXTENDED layout:
//   row_ind: sorted, global row ids, covering BOTH U rows AND panel rows.
//     First u_row_count entries: U_segment rows (global row < col_begin_j), sorted.
//     Next (rows - u_row_count) entries: panel rows = desc_j.row_indices, sorted.
//     Since U rows < col_begin_j <= panel rows, the combined row_ind is sorted.
//
//   values: column-major, rows × cols, leading dimension = ld = rows.
//     values[r + c * ld] = (workspace row r, panel column c) entry.
//     Rows 0..u_row_count-1: U_segment values (gathered from U_segments).
//     Rows u_row_count..rows-1: panel values (gathered from panel_values).
//
//   u_row_count: number of U_segment rows at the start of the workspace.
//
//   valid: true iff workspace was successfully gathered.
//
// NOTE: The gather/scatter includes U_segment rows so that apply_supernode_panel_update
//   can find inter_j positions in the U part (where k's column range intersects).
//   After scatter, both panel_values and U_segments are updated (transitional state).
// ---------------------------------------------------------------------------
template <class T, class Index>
struct supernode_panel_workspace {
    static_assert(std::is_signed<Index>::value,
                  "supernode_panel_workspace: Index must be signed");

    std::vector<Index> row_ind;  // global row ids: U rows ++ panel rows (sorted)
    std::vector<T>     values;   // dense panel: rows x num_cols, column-major
    Index rows;                  // total rows = u_row_count + panel_row_count
    Index cols;                  // num_cols of this supernode
    Index ld;                    // leading dimension (== rows, no padding)
    bool  valid;
    std::size_t u_row_count;     // count of U rows at start of row_ind / values

    // SLU-SN-OPT: reusable inverse map  global_row -> position in row_ind  (else -1).
    // Sized to global n and reused across supernodes (gather clears its own previous
    // footprint, then refills). Replaces the per-call lower_bound scans in gather /
    // apply (off-diagonal search) / scatter with O(1) lookups. The positions found
    // are identical to the lower_bound results, so numerics are byte-identical.
    std::vector<Index> row_pos;

    supernode_panel_workspace()
        : rows(Index(0)), cols(Index(0)), ld(Index(0)),
          valid(false), u_row_count(0u) {}
};

// ---------------------------------------------------------------------------
// supernode_update_set<Index>
//
// Set of updating supernodes for a given current supernode.
// SLU-8R.3.1: criterion uses U_segments rows (not row_indices).
// ---------------------------------------------------------------------------
template <class Index>
struct supernode_update_set {
    static_assert(std::is_signed<Index>::value,
                  "supernode_update_set: Index must be signed");

    std::vector<Index> updating_supernodes;  // ascending order, k < current
    bool conservative;                       // false = exact U-segment criterion
    bool symmetric_pruning_active;           // ALWAYS false in SLU-8R.3.1 (hook)
    std::size_t pruned_edges;               // ALWAYS 0 in SLU-8R.3.1

    supernode_update_set()
        : conservative(false),
          symmetric_pruning_active(false),
          pruned_edges(0u) {}
};

// supernode_panel_update_stats is defined in tsparse_sparse_lu.hpp
// (in namespace sparse_lu_detail, before sparse_lu_factorization class body)
// so that set_panel_update_info_ can use it inline in the class body.
// Do NOT redefine it here.

// ---------------------------------------------------------------------------
// compute_panel_update_set<T, Index>
//
// SLU-8R.3.1 FIXED panel symbolic reachability: determines the set of prior
// supernodes that update current supernode j's panel.
//
// FIXED CRITERION (U-segment rows):
//   Supernode k < j updates j iff j's U_segments contain at least one
//   global row r such that:
//     supernodes[k].first_col <= r < supernodes[k].first_col + supernodes[k].num_cols
//
//   U_segment rows are rows < col_begin_j (off-diagonal U entries for j's columns).
//   k < j means col_end_k <= col_begin_j, so U_segment rows CAN be in k's range.
//
//   This is the correct criterion for §17.2(A): K has a non-zero U connection
//   to J iff J's U_segment contains a row in K's column range.
//
// PREVIOUS BROKEN CRITERION (row_indices):
//   row_indices contains rows >= col_begin_j. For k < j, k's range is < col_begin_j.
//   Intersection was always empty → update_count was always 0.
//
// Returns: supernode_update_set with updating_supernodes in ascending order.
// ---------------------------------------------------------------------------
template <class T, class Index>
supernode_update_set<Index>
compute_panel_update_set(
    const supernodal_lu_storage<T, Index>& storage,
    Index current_supernode)
{
    supernode_update_set<Index> result;

    if (!storage.valid) return result;
    const std::size_t nsup = storage.supernodes.size();
    const std::size_t j    = static_cast<std::size_t>(current_supernode);
    if (j >= nsup) return result;

    // Get J's U_segment entry range.
    if (j + 1u >= storage.U_segments.seg_ptr.size()) return result;
    const Index j_seg_begin = storage.U_segments.seg_ptr[j];
    const Index j_seg_end   = storage.U_segments.seg_ptr[j + 1u];

    if (j_seg_begin >= j_seg_end) {
        // No U_segment rows for J → no prior supernodes can contribute.
        result.conservative = false;
        return result;
    }

    // Collect unique U_segment rows for J (may have duplicates: one per column).
    std::vector<Index> u_rows_j(
        storage.U_segments.row_ind.begin() +
            static_cast<std::ptrdiff_t>(j_seg_begin),
        storage.U_segments.row_ind.begin() +
            static_cast<std::ptrdiff_t>(j_seg_end));
    std::sort(u_rows_j.begin(), u_rows_j.end());
    u_rows_j.erase(
        std::unique(u_rows_j.begin(), u_rows_j.end()),
        u_rows_j.end());

    // For each prior supernode k, check if any U_row_j is in k's column range.
    for (std::size_t k = 0u; k < j; ++k) {
        const supernode_desc<Index>& desc_k = storage.supernodes[k];
        const Index col_begin_k = desc_k.first_col;
        const Index col_end_k   = col_begin_k + desc_k.num_cols;

        if (col_begin_k >= col_end_k) continue;

        // Binary search in sorted u_rows_j for first row >= col_begin_k.
        typename std::vector<Index>::const_iterator it =
            std::lower_bound(u_rows_j.begin(), u_rows_j.end(), col_begin_k);

        if (it != u_rows_j.end() && *it < col_end_k) {
            result.updating_supernodes.push_back(static_cast<Index>(k));
        }
    }

    result.conservative             = false;  // exact U-segment criterion
    result.symmetric_pruning_active = false;  // SLU-8R.3.1: hook only
    result.pruned_edges             = 0u;

    return result;
}

// ---------------------------------------------------------------------------
// compute_panel_update_set<T, Index>  (SLU-SN-OPT fast overload)
//
// Identical result to the range-scan overload above, but O(#U-rows of j)
// instead of O(nsup): instead of scanning every prior supernode k and binary
// searching j's U rows, map each U-segment row directly to its owning supernode
// via a precomputed col_to_supernode array (built once per factorization, O(n)).
//
// Each U-segment row r is, by construction, a row < col_begin_j, i.e. a COLUMN
// owned by exactly one earlier supernode k < j; col_to_supernode[r] = k. So the
// updating-supernode set is { col_to_supernode[r] : r in U_segment(j) }, which
// equals the set produced by the range-scan criterion. Sorting + unique restores
// ascending order, so the §17.2(A) update sequence is byte-identical and the
// numeric factorization is unchanged (this is a pure symbolic speedup).
// ---------------------------------------------------------------------------
template <class T, class Index>
supernode_update_set<Index>
compute_panel_update_set(
    const supernodal_lu_storage<T, Index>& storage,
    Index current_supernode,
    const std::vector<Index>& col_to_supernode)
{
    supernode_update_set<Index> result;

    if (!storage.valid) return result;
    const std::size_t nsup = storage.supernodes.size();
    const std::size_t j    = static_cast<std::size_t>(current_supernode);
    if (j >= nsup) return result;

    if (j + 1u >= storage.U_segments.seg_ptr.size()) return result;
    const Index j_seg_begin = storage.U_segments.seg_ptr[j];
    const Index j_seg_end   = storage.U_segments.seg_ptr[j + 1u];

    if (j_seg_begin >= j_seg_end) {
        result.conservative = false;
        return result;
    }

    for (Index idx = j_seg_begin; idx < j_seg_end; ++idx) {
        const std::size_t sidx = static_cast<std::size_t>(idx);
        if (sidx >= storage.U_segments.row_ind.size()) continue;
        const Index r = storage.U_segments.row_ind[sidx];
        if (r < Index(0) ||
            static_cast<std::size_t>(r) >= col_to_supernode.size()) continue;
        const Index k = col_to_supernode[static_cast<std::size_t>(r)];
        if (k < Index(0) || static_cast<std::size_t>(k) >= j) continue;
        result.updating_supernodes.push_back(k);
    }

    std::sort(result.updating_supernodes.begin(),
              result.updating_supernodes.end());
    result.updating_supernodes.erase(
        std::unique(result.updating_supernodes.begin(),
                    result.updating_supernodes.end()),
        result.updating_supernodes.end());

    result.conservative             = false;
    result.symmetric_pruning_active = false;
    result.pruned_edges             = 0u;

    return result;
}

// ---------------------------------------------------------------------------
// gather_panel_workspace<T, Index>
//
// SLU-8R.3.1 EXTENDED: Gathers the dense panel for current supernode j into
// a workspace that covers BOTH U_segment rows AND panel rows.
//
// workspace layout:
//   row_ind: [U_seg_rows (sorted)] ++ [desc_j.row_indices (sorted)]
//            All U rows < col_begin_j <= panel rows, so combined is sorted.
//   values: column-major, (u_row_count + panel_row_count) x w_j.
//     U part (rows 0..u_row_count-1): from U_segments via u_seg_col_ptr.
//     Panel part (rows u_row_count..total-1): from panel_values.
//   u_row_count: number of U_segment unique rows.
//
// Requires: desc_j.u_seg_col_ptr populated (done by bootstrap).
// ---------------------------------------------------------------------------
template <class T, class Index>
void gather_panel_workspace(
    const supernodal_lu_storage<T, Index>& storage,
    Index current_supernode,
    supernode_panel_workspace<T, Index>& work)
{
    // SLU-SN-OPT: clear this workspace's previous global_row->position footprint
    // before the row_ind it indexes is discarded, so row_pos is reusable across
    // supernodes without an O(n) wipe each call.
    for (std::size_t i = 0u; i < work.row_ind.size(); ++i) {
        const Index gr = work.row_ind[i];
        if (gr >= Index(0) &&
            static_cast<std::size_t>(gr) < work.row_pos.size())
            work.row_pos[static_cast<std::size_t>(gr)] = Index(-1);
    }

    work.valid       = false;
    work.rows        = Index(0);
    work.cols        = Index(0);
    work.ld          = Index(0);
    work.u_row_count = 0u;
    work.row_ind.clear();
    work.values.clear();

    if (!storage.valid) return;
    const std::size_t j = static_cast<std::size_t>(current_supernode);
    if (j >= storage.supernodes.size()) return;

    const supernode_desc<Index>& desc_j = storage.supernodes[j];
    const Index w_j       = desc_j.num_cols;
    const Index ld_j      = desc_j.leading_dimension;
    const Index panel_row_count = static_cast<Index>(desc_j.row_indices.size());

    if (w_j <= Index(0) || panel_row_count <= Index(0)) return;
    if (ld_j < panel_row_count) return;

    // Check panel_values bounds for panel rows.
    const std::size_t panel_start  = static_cast<std::size_t>(desc_j.values_offset);
    const std::size_t panel_end    =
        panel_start +
        static_cast<std::size_t>(ld_j) * static_cast<std::size_t>(w_j);
    if (panel_end > storage.panel_values.size()) return;

    // -----------------------------------------------------------------------
    // Part 1: U_segment rows for J.
    // -----------------------------------------------------------------------
    std::vector<Index> u_rows_sorted;
    if (j + 1u < storage.U_segments.seg_ptr.size()) {
        const Index j_seg_begin = storage.U_segments.seg_ptr[j];
        const Index j_seg_end   = storage.U_segments.seg_ptr[j + 1u];
        if (j_seg_begin < j_seg_end) {
            std::vector<Index> u_tmp(
                storage.U_segments.row_ind.begin() +
                    static_cast<std::ptrdiff_t>(j_seg_begin),
                storage.U_segments.row_ind.begin() +
                    static_cast<std::ptrdiff_t>(j_seg_end));
            std::sort(u_tmp.begin(), u_tmp.end());
            u_tmp.erase(std::unique(u_tmp.begin(), u_tmp.end()), u_tmp.end());
            u_rows_sorted = u_tmp;
        }
    }
    const std::size_t u_row_cnt = u_rows_sorted.size();

    // -----------------------------------------------------------------------
    // Combined row_ind: U rows ++ panel rows (sorted since U < col_begin_j <= panel).
    // -----------------------------------------------------------------------
    work.row_ind.clear();
    work.row_ind.reserve(u_row_cnt + static_cast<std::size_t>(panel_row_count));
    work.row_ind.insert(work.row_ind.end(),
                        u_rows_sorted.begin(), u_rows_sorted.end());
    work.row_ind.insert(work.row_ind.end(),
                        desc_j.row_indices.begin(), desc_j.row_indices.end());

    const Index total_rows = static_cast<Index>(u_row_cnt) + panel_row_count;
    work.rows        = total_rows;
    work.cols        = w_j;
    work.ld          = total_rows;
    work.u_row_count = u_row_cnt;

    const std::size_t total_values =
        static_cast<std::size_t>(total_rows) * static_cast<std::size_t>(w_j);
    work.values.assign(total_values, T(0));

    // SLU-SN-OPT: build the global_row -> position map for this workspace.
    const std::size_t gn = storage.row_perm.size();
    if (work.row_pos.size() < gn) work.row_pos.assign(gn, Index(-1));
    for (std::size_t i = 0u; i < work.row_ind.size(); ++i) {
        const Index gr = work.row_ind[i];
        if (gr >= Index(0) && static_cast<std::size_t>(gr) < work.row_pos.size())
            work.row_pos[static_cast<std::size_t>(gr)] = static_cast<Index>(i);
    }

    // -----------------------------------------------------------------------
    // Fill U part from U_segments using u_seg_col_ptr.  Position via row_pos.
    // -----------------------------------------------------------------------
    if (u_row_cnt > 0u &&
        !desc_j.u_seg_col_ptr.empty() &&
        desc_j.u_seg_col_ptr.size() == static_cast<std::size_t>(w_j) + 1u)
    {
        for (Index c = Index(0); c < w_j; ++c) {
            const std::size_t sc = static_cast<std::size_t>(c);
            const Index rel_start = desc_j.u_seg_col_ptr[sc];
            const Index rel_end   = desc_j.u_seg_col_ptr[sc + 1u];
            for (Index ki = rel_start; ki < rel_end; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(desc_j.u_segment_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                const T     u_val = storage.U_segments.values[abs_idx];
                const Index pos = (static_cast<std::size_t>(u_row) < work.row_pos.size())
                    ? work.row_pos[static_cast<std::size_t>(u_row)] : Index(-1);
                if (pos >= Index(0) &&
                    static_cast<std::size_t>(pos) < u_row_cnt) {
                    work.values[static_cast<std::size_t>(pos) +
                                sc * static_cast<std::size_t>(total_rows)] = u_val;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Fill panel part from panel_values.
    // -----------------------------------------------------------------------
    const T* panel_j = &storage.panel_values[panel_start];
    for (Index c = Index(0); c < w_j; ++c) {
        const std::size_t sc = static_cast<std::size_t>(c);
        for (Index r = Index(0); r < panel_row_count; ++r) {
            const std::size_t sr = static_cast<std::size_t>(r);
            work.values[(u_row_cnt + sr) + sc * static_cast<std::size_t>(total_rows)] =
                panel_j[sr + sc * static_cast<std::size_t>(ld_j)];
        }
    }

    work.valid = true;
}

// ---------------------------------------------------------------------------
// apply_supernode_panel_update<T, Index>
//
// Applies the §17.2(A) left-looking update from updating supernode k to
// current supernode j's workspace.
//
// SLU-8R.3.1 NOTE: No change needed here. The function already uses work.rows
// for stride and work.row_ind for lookups. With the extended workspace (U rows
// at the front), inter_j will find U_segment rows that fall in k's column range,
// which is the correct criterion.
//
// Algorithm:
//   1. Find inter_j/inter_k: positions in j's workspace where global row
//      falls within k's column range [col_begin_k, col_end_k).
//      With extended workspace, these will be U rows (rows < col_begin_j).
//
//   2. Build Z (w_k x w_j): Z[inter_k[i], c] = work.values[inter_j[i] + c*rows]
//
//   3. trsm: L_k_diag^{-1} * Z (unit lower triangular, in-place).
//
//   4. Scatter Z[inter_k] back to workspace[inter_j] (updates U workspace).
//
//   5. Find off_j/off_k: k's off-diagonal L rows that appear in j's workspace.
//      With extended workspace, these can be in either U or panel part of j.
//
//   6. gemm/gemv: workspace[off_j, :] -= L_off_sub * Z_inter.
//      Scatter C back to workspace (updates affected rows).
//
// Dense kernel calls (via adapter, NOT direct tblas/tlapack):
//   trsm: sparse_lu_dense_kernel<T>::trsm ('L','L','N','U', ...)
//   gemm: sparse_lu_dense_kernel<T>::gemm (...)
//   gemv: sparse_lu_dense_kernel<T>::gemv (...) for w_j==1 case
//
// getrf is NOT called in §17.2(A) (no pivot search; §25 prohibition).
//
// SLU-SN-OPT: this function is called once per (k,j) update pair — the hottest
// non-kernel path in the factorization (panelNK dominated total numeric time).
// To eliminate the per-call heap churn (previously 4 vector allocations per
// call) and the linear inter-row scan, the caller threads a reusable
// panel_apply_scratch buffer through, and the inter-row search uses binary
// search over the sorted workspace row_ind. Both are pure speedups: the values
// built, the dense-kernel calls, and their order are byte-identical.
// ---------------------------------------------------------------------------
template <class T, class Index>
struct panel_apply_scratch {
    std::vector<std::size_t> inter_j;
    std::vector<std::size_t> inter_k;
    std::vector<std::size_t> off_j;
    std::vector<std::size_t> off_k;
    std::vector<T>           Z;
    std::vector<T>           Z_inter;
    std::vector<T>           L_off_sub;
    std::vector<T>           C;
};

template <class T, class Index>
void apply_supernode_panel_update(
    const supernodal_lu_storage<T, Index>& storage,
    Index updating_supernode,
    Index current_supernode,
    supernode_panel_workspace<T, Index>& work,
    supernode_panel_update_stats& stats,
    panel_apply_scratch<T, Index>& scr);

// Convenience wrapper: allocates a transient scratch (non-hot callers/tests).
template <class T, class Index>
void apply_supernode_panel_update(
    const supernodal_lu_storage<T, Index>& storage,
    Index updating_supernode,
    Index current_supernode,
    supernode_panel_workspace<T, Index>& work,
    supernode_panel_update_stats& stats)
{
    panel_apply_scratch<T, Index> scr;
    apply_supernode_panel_update(storage, updating_supernode,
                                 current_supernode, work, stats, scr);
}

template <class T, class Index>
void apply_supernode_panel_update(
    const supernodal_lu_storage<T, Index>& storage,
    Index updating_supernode,
    Index /*current_supernode*/,
    supernode_panel_workspace<T, Index>& work,
    supernode_panel_update_stats& stats,
    panel_apply_scratch<T, Index>& scr)
{
    if (!work.valid) return;

    const std::size_t k = static_cast<std::size_t>(updating_supernode);
    if (k >= storage.supernodes.size()) return;

    const supernode_desc<Index>& desc_k    = storage.supernodes[k];
    const Index col_begin_k  = desc_k.first_col;
    const Index col_end_k    = col_begin_k + desc_k.num_cols;
    const std::size_t w_k    = static_cast<std::size_t>(desc_k.num_cols);
    const std::size_t ld_k   = static_cast<std::size_t>(desc_k.leading_dimension);
    const std::size_t row_count_k = desc_k.row_indices.size();

    if (w_k == 0u || ld_k < w_k) return;

    // Check panel_k bounds.
    const std::size_t panel_k_start = static_cast<std::size_t>(desc_k.values_offset);
    const std::size_t panel_k_end   = panel_k_start + ld_k * w_k;
    if (panel_k_end > storage.panel_values.size()) return;

    const T* panel_k = &storage.panel_values[panel_k_start];

    const std::size_t w_j          = static_cast<std::size_t>(work.cols);
    const std::size_t row_count_j  = static_cast<std::size_t>(work.rows);  // U + panel
    const std::vector<Index>& row_ind_j = work.row_ind;

    if (w_j == 0u || row_count_j == 0u) return;

    // ------------------------------------------------------------------
    // Step 1: Find connection rows (U rows of j in k's column range).
    // With extended workspace, inter_j finds positions in U part of workspace.
    // ------------------------------------------------------------------
    std::vector<std::size_t>& inter_j = scr.inter_j;   // positions in workspace rows
    std::vector<std::size_t>& inter_k = scr.inter_k;   // local index in k's diagonal block
    inter_j.clear();
    inter_k.clear();

    // row_ind_j is sorted ascending; the rows in [col_begin_k, col_end_k) form a
    // contiguous block. Binary search to the start, then walk while < col_end_k.
    // Identical rows and order as the prior O(row_count_j) linear scan.
    {
        typename std::vector<Index>::const_iterator lo =
            std::lower_bound(row_ind_j.begin(), row_ind_j.end(), col_begin_k);
        for (; lo != row_ind_j.end() && *lo < col_end_k; ++lo) {
            inter_j.push_back(
                static_cast<std::size_t>(lo - row_ind_j.begin()));
            inter_k.push_back(
                static_cast<std::size_t>(*lo - col_begin_k));
        }
    }

    if (inter_j.empty()) return;  // no connection

    const std::size_t n_inter = inter_j.size();

    // ------------------------------------------------------------------
    // Step 2: Build Z (w_k x w_j, column-major).
    // ------------------------------------------------------------------
    std::vector<T>& Z = scr.Z;
    Z.assign(w_k * w_j, T(0));
    for (std::size_t c = 0u; c < w_j; ++c) {
        for (std::size_t i = 0u; i < n_inter; ++i) {
            Z[inter_k[i] + c * w_k] =
                work.values[inter_j[i] + c * row_count_j];
        }
    }

    // ------------------------------------------------------------------
    // Step 3: trsm -- solve L_k_diag * Z = Z (unit lower tri, in-place).
    // ADAPTER CALL (not direct tblas/tlapack).
    // ------------------------------------------------------------------
    {
        const auto t0 = std::chrono::steady_clock::now();
        sparse_lu_dense_kernel<T>::trsm(
            'L', 'L', 'N', 'U',
            w_k, w_j,
            T(1), panel_k, ld_k,
            Z.data(), w_k);
        const auto t1 = std::chrono::steady_clock::now();
        stats.trsm_count++;
        // SLU-PERF: trsm(m,n) FLOP = m*m*n, m=w_k (triangular order), n=w_j (rhs).
        stats.flop_trsm +=
            static_cast<double>(w_k) * static_cast<double>(w_k) *
            static_cast<double>(w_j);
        stats.dense_kernel_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                t1 - t0).count());
    }

    // ------------------------------------------------------------------
    // Step 4: Scatter trsm result back to workspace at inter_j positions.
    // ------------------------------------------------------------------
    for (std::size_t c = 0u; c < w_j; ++c) {
        for (std::size_t i = 0u; i < n_inter; ++i) {
            work.values[inter_j[i] + c * row_count_j] =
                Z[inter_k[i] + c * w_k];
        }
    }

    // ------------------------------------------------------------------
    // Step 5: Off-diagonal update (gemm/gemv).
    // k's off-diagonal rows (positions w_k..row_count_k-1 in k's panel)
    // that also appear in j's workspace (could be U or panel part).
    // ADAPTER CALLS (not direct tblas/tlapack).
    // ------------------------------------------------------------------
    std::vector<std::size_t>& off_j = scr.off_j;
    std::vector<std::size_t>& off_k = scr.off_k;
    off_j.clear();
    off_k.clear();

    for (std::size_t r_k = w_k; r_k < row_count_k; ++r_k) {
        if (r_k >= desc_k.row_indices.size()) break;
        const Index global_row = desc_k.row_indices[r_k];
        // SLU-SN-OPT: O(1) position lookup instead of lower_bound over row_ind_j.
        const Index pos = (global_row >= Index(0) &&
                           static_cast<std::size_t>(global_row) < work.row_pos.size())
            ? work.row_pos[static_cast<std::size_t>(global_row)] : Index(-1);
        if (pos >= Index(0)) {
            off_j.push_back(static_cast<std::size_t>(pos));
            off_k.push_back(r_k);
        }
    }

    if (!off_j.empty()) {
        const std::size_t n_off = off_j.size();

        // Build Z_inter (n_inter x w_j, column-major).
        std::vector<T>& Z_inter = scr.Z_inter;
        Z_inter.resize(n_inter * w_j);
        for (std::size_t c = 0u; c < w_j; ++c) {
            for (std::size_t i = 0u; i < n_inter; ++i) {
                Z_inter[i + c * n_inter] = Z[inter_k[i] + c * w_k];
            }
        }

        // Build L_off_sub (n_off x n_inter, column-major).
        std::vector<T>& L_off_sub = scr.L_off_sub;
        L_off_sub.resize(n_off * n_inter);
        for (std::size_t i_inter = 0u; i_inter < n_inter; ++i_inter) {
            for (std::size_t i_off = 0u; i_off < n_off; ++i_off) {
                L_off_sub[i_off + i_inter * n_off] =
                    panel_k[off_k[i_off] + inter_k[i_inter] * ld_k];
            }
        }

        // Build C (n_off x w_j) from workspace at off_j rows.
        std::vector<T>& C = scr.C;
        C.resize(n_off * w_j);
        for (std::size_t c = 0u; c < w_j; ++c) {
            for (std::size_t i = 0u; i < n_off; ++i) {
                C[i + c * n_off] = work.values[off_j[i] + c * row_count_j];
            }
        }

        {
            const auto t0 = std::chrono::steady_clock::now();
            if (w_j == 1u) {
                sparse_lu_dense_kernel<T>::gemv(
                    n_off, n_inter,
                    T(-1), L_off_sub.data(), n_off,
                    Z_inter.data(), T(1), C.data());
                stats.gemv_count++;
                // SLU-PERF: gemv(m,n) FLOP = 2*m*n, m=n_off, n=n_inter.
                stats.flop_gemv +=
                    2.0 * static_cast<double>(n_off) *
                    static_cast<double>(n_inter);
            } else {
                sparse_lu_dense_kernel<T>::gemm(
                    n_off, w_j, n_inter,
                    T(-1), L_off_sub.data(), n_off,
                    Z_inter.data(), n_inter,
                    T(1), C.data(), n_off);
                stats.gemm_count++;
                // SLU-PERF: gemm(m,n,k) FLOP = 2*m*n*k; record shape (m=n_off,
                // n=w_j, k=n_inter) for the line-894 BLAS-3 concentration verdict.
                stats.flop_gemm +=
                    2.0 * static_cast<double>(n_off) *
                    static_cast<double>(w_j) * static_cast<double>(n_inter);
                stats.gemm_m_sum += n_off;
                stats.gemm_n_sum += w_j;
                stats.gemm_k_sum += n_inter;
                if (n_off  > stats.gemm_max_m) stats.gemm_max_m = n_off;
                if (w_j    > stats.gemm_max_n) stats.gemm_max_n = w_j;
                if (n_inter > stats.gemm_max_k) stats.gemm_max_k = n_inter;
                std::size_t mind = n_off;
                if (w_j    < mind) mind = w_j;
                if (n_inter < mind) mind = n_inter;
                ++stats.gemm_dim_hist[gemm_dim_bucket(mind)];
            }
            const auto t1 = std::chrono::steady_clock::now();
            stats.dense_kernel_ticks += static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    t1 - t0).count());
        }

        // Scatter C back to workspace.
        for (std::size_t c = 0u; c < w_j; ++c) {
            for (std::size_t i = 0u; i < n_off; ++i) {
                work.values[off_j[i] + c * row_count_j] = C[i + c * n_off];
            }
        }
    }

    stats.update_count++;
}

// ---------------------------------------------------------------------------
// scatter_panel_workspace<T, Index>
//
// SLU-8R.3.1 EXTENDED: Scatters the updated workspace back to both:
//   (a) U_segments.values: for U rows (work rows 0..u_row_count-1)
//   (b) panel_values: for panel rows (work rows u_row_count..total-1)
//
// Both panel_values and U_segments are transitional state after §17.2(A).
//   true_numeric_source remains false (§17.2(B) PENDING).
//   Solve uses baseline CSC (Gate 4 maintained).
//
// Modifying U_segments is safe: production solve does not use U_segments.
// ---------------------------------------------------------------------------
template <class T, class Index>
void scatter_panel_workspace(
    supernodal_lu_storage<T, Index>& storage,
    Index current_supernode,
    const supernode_panel_workspace<T, Index>& work,
    supernode_panel_update_stats& stats)
{
    if (!work.valid || !storage.valid) return;
    const std::size_t j = static_cast<std::size_t>(current_supernode);
    if (j >= storage.supernodes.size()) return;

    const supernode_desc<Index>& desc_j = storage.supernodes[j];
    const Index w_j             = desc_j.num_cols;
    const Index ld_j            = desc_j.leading_dimension;
    const Index panel_row_count = static_cast<Index>(desc_j.row_indices.size());

    if (w_j <= Index(0) || panel_row_count <= Index(0)) return;

    // Expected workspace total rows: u_row_count + panel_row_count.
    const Index expected_rows =
        static_cast<Index>(work.u_row_count) + panel_row_count;
    if (work.rows != expected_rows || work.cols != w_j) return;

    const std::size_t total_rows = static_cast<std::size_t>(work.rows);
    const std::size_t u_row_cnt  = work.u_row_count;

    // -----------------------------------------------------------------------
    // (a) Scatter U rows back to U_segments.
    // -----------------------------------------------------------------------
    if (u_row_cnt > 0u &&
        j + 1u < storage.U_segments.seg_ptr.size() &&
        !desc_j.u_seg_col_ptr.empty() &&
        desc_j.u_seg_col_ptr.size() == static_cast<std::size_t>(w_j) + 1u)
    {
        for (Index c = Index(0); c < w_j; ++c) {
            const std::size_t sc = static_cast<std::size_t>(c);
            const Index rel_start = desc_j.u_seg_col_ptr[sc];
            const Index rel_end   = desc_j.u_seg_col_ptr[sc + 1u];
            for (Index ki = rel_start; ki < rel_end; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(desc_j.u_segment_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                // SLU-SN-OPT: O(1) position lookup; U rows live in [0,u_row_cnt).
                const Index pos = (static_cast<std::size_t>(u_row) < work.row_pos.size())
                    ? work.row_pos[static_cast<std::size_t>(u_row)] : Index(-1);
                if (pos >= Index(0) &&
                    static_cast<std::size_t>(pos) < u_row_cnt)
                {
                    storage.U_segments.values[abs_idx] =
                        work.values[static_cast<std::size_t>(pos) + sc * total_rows];
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // (b) Scatter panel rows back to panel_values.
    // -----------------------------------------------------------------------
    const std::size_t panel_start  = static_cast<std::size_t>(desc_j.values_offset);
    const std::size_t panel_end    =
        panel_start +
        static_cast<std::size_t>(ld_j) * static_cast<std::size_t>(w_j);
    if (panel_end > storage.panel_values.size()) return;

    T* panel_j = &storage.panel_values[panel_start];
    for (Index c = Index(0); c < w_j; ++c) {
        const std::size_t sc = static_cast<std::size_t>(c);
        for (Index r = Index(0); r < panel_row_count; ++r) {
            const std::size_t sr = static_cast<std::size_t>(r);
            panel_j[sr + sc * static_cast<std::size_t>(ld_j)] =
                work.values[(u_row_cnt + sr) + sc * total_rows];
        }
    }

    stats.scatter_count++;
}

// ---------------------------------------------------------------------------
// run_supernode_panel_leftlooking_update<T, Index>
//
// SLU-8R.3.1 FIXED top-level §17.2(A) driver.
//
// Changes from SLU-8R.3:
//   - stats.driver_called = true (always set at entry)
//   - compute_panel_update_set uses U_segment rows → update_count > 0 for
//     matrices with off-diagonal U structure (previously always 0)
//   - stats.update_applied = true only when update_count > 0 (SLU-8R.7.1:
//     authoritative as info.supernode_panel_update_applied; no storage flag)
//   - stats.update_applied, stats.dense_kernel_called computed at end
//
// For each supernode j in forward order:
//   1. Compute update set via U-segment criterion.
//   2. If non-empty: gather workspace (U + panel rows).
//   3. Apply updates from each k in update set (trsm + gemm/gemv via adapter).
//   4. Scatter updated workspace back to U_segments and panel_values.
//
// After completing all supernodes:
//   stats.update_applied        = (update_count > 0)  [SLU-8R.7.1: info-side]
//   storage.true_numeric_source = false (unchanged; §17.2(B) PENDING)
// ---------------------------------------------------------------------------
template <class T, class Index>
supernode_panel_update_stats
run_supernode_panel_leftlooking_update(
    supernodal_lu_storage<T, Index>& storage)
{
    supernode_panel_update_stats stats;
    stats.driver_called = true;

    if (!storage.valid) {
        return stats;
    }

    const std::size_t nsup = storage.supernodes.size();

    for (std::size_t j = 0u; j < nsup; ++j) {
        const Index j_idx = static_cast<Index>(j);

        supernode_update_set<Index> update_set =
            compute_panel_update_set(storage, j_idx);

        if (update_set.updating_supernodes.empty()) continue;

        supernode_panel_workspace<T, Index> work;
        gather_panel_workspace(storage, j_idx, work);

        if (!work.valid) continue;

        for (std::size_t ui = 0u;
             ui < update_set.updating_supernodes.size(); ++ui)
        {
            apply_supernode_panel_update(
                storage,
                update_set.updating_supernodes[ui],
                j_idx,
                work,
                stats);
        }

        scatter_panel_workspace(storage, j_idx, work, stats);
        stats.panel_count++;
    }

    stats.used_conservative_reach   = false;
    stats.symmetric_pruning_active  = false;

    // Compute derived diagnostics.
    stats.update_applied    = (stats.update_count > 0u && stats.scatter_count > 0u);
    stats.dense_kernel_called =
        (stats.trsm_count > 0u &&
         (stats.gemm_count > 0u || stats.gemv_count > 0u));

    // SLU-8R.7.1: the §17.2(A) "panel update applied" diagnostic is authoritative
    // in sparse_lu_info (info.supernode_panel_update_applied), fed from
    // stats.update_applied below. No storage-resident flag is written: after the
    // SLU-8R.7 re-bootstrap, final storage is the A_eff-origin true-numeric source
    // storage, and a storage flag here would be reset by that re-bootstrap.
    // true_numeric_source remains false (§17.2(B) PENDING).

    return stats;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_SUPERNODE_PANEL_UPDATE_IMPL_HPP
