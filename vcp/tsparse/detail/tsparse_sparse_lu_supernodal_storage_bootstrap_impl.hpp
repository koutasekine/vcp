// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-8R.2: Supernodal storage bootstrap from CSC L/U -- internal implementation.
//
// Provides bootstrap_supernodal_storage_from_csc in namespace sparse_lu_detail.
//
// This is a TRANSITIONAL bootstrap helper: it populates supernodal_lu_storage
// from the baseline CSC L/U result produced by baseline_sparse_gp_lu_factorize.
// The resulting storage is flagged as:
//   bootstrapped_from_csc = true  (values came from CSC L/U, not §17.2)
//   source_of_truth_storage = true (factor owns this as supernodal storage)
//   true_numeric_source = false    (§17.2(A)/(B) NOT yet implemented)
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - all storage type definitions (csc_storage, baseline_lu_storage,
//     supernode_desc, supernodal_lu_storage, u_segment_storage)
//   - sparse_lu_supernode_numeric<T, Index>
//   - sparse_lu_is_valid_csc_storage (from SLU-10.1 numeric section)
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//
// SLU-8 Gate status after SLU-8R.2:
//   Gate 1: PASS
//   Gate 2: PENDING
//   Gate 3: PENDING-REQUIRED
//   Gate 4: PASS
//   Gate 5: PENDING-REQUIRED
//   Gate 6: PENDING (baseline_sparse_gp still used as numeric source)
//
// §17.2(A)/(B) supernode-panel update / within-panel factorization: PENDING.
// §18.2 storage-native supernodal solve: PENDING.
// issue_SLU8_contract_violation.md: OPEN.

#ifndef VCP_TSPARSE_SPARSE_LU_SUPERNODAL_STORAGE_BOOTSTRAP_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SUPERNODAL_STORAGE_BOOTSTRAP_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// csc_lookup_value
//
// Looks up the value at (row, col) in a CSC matrix stored as col_ptr/row_ind/values.
// Uses binary search on the sorted row_ind for column col.
// Returns T(0) if the entry is structurally absent.
// Precondition: col is a valid column index (col_ptr has at least col+2 entries).
// ---------------------------------------------------------------------------
template <class T, class Index>
T csc_lookup_value(
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     values,
    Index col,
    Index row)
{
    const std::size_t scol = static_cast<std::size_t>(col);
    const Index kb = col_ptr[scol];
    const Index ke = col_ptr[scol + 1u];
    const typename std::vector<Index>::const_iterator begin_it =
        row_ind.begin() + static_cast<std::ptrdiff_t>(kb);
    const typename std::vector<Index>::const_iterator end_it =
        row_ind.begin() + static_cast<std::ptrdiff_t>(ke);
    const typename std::vector<Index>::const_iterator it =
        std::lower_bound(begin_it, end_it, row);
    if (it == end_it || *it != row) return T(0);
    return values[static_cast<std::size_t>(it - row_ind.begin())];
}

// ---------------------------------------------------------------------------
// bootstrap_supernodal_storage_from_csc
//
// SLU-8R.2 transitional bootstrap: builds a supernodal_lu_storage<T, Index>
// from the baseline CSC L/U factors produced by baseline_sparse_gp_lu_factorize.
//
// This bootstrap:
//   - copies permutations and scaling from baseline_lu_storage
//   - builds supernode descriptors from sparse_lu_supernode_numeric partition
//   - constructs combined dense panels (row_count x num_cols, column-major):
//       diagonal block rows: combined L\U (U on/above diag, L below diag)
//       off-diagonal rows:   L entries below the diagonal block
//   - populates U_segments with off-diagonal U entries (rows < first_col)
//     stored column-by-column within each supernode's segment
//
// SLU-8R.2 status flags:
//   bootstrapped_from_csc = true
//   source_of_truth_storage = true
//   true_numeric_source = false  (§17.2 NOT implemented, Gate 6 PENDING)
//   valid = true on success
//
// §18.2 storage-native solve PENDING: U_segments column association is implicit
// (entries ordered col_begin first, then col_begin+1, etc. within each segment).
//
// Returns invalid storage (valid == false) if sn_num.valid is false or n <= 0.
// ---------------------------------------------------------------------------
template <class T, class Index>
supernodal_lu_storage<T, Index>
bootstrap_supernodal_storage_from_csc(
    Index n,
    const baseline_lu_storage<T, Index>& csc_lu,
    const sparse_lu_supernode_numeric<T, Index>& sn_num,
    bool add_alignment_padding = false)
{
    static_assert(std::is_signed<Index>::value,
                  "bootstrap_supernodal_storage_from_csc: Index must be signed");

    supernodal_lu_storage<T, Index> result;
    // SLU-8R.2 transition flags.
    result.bootstrapped_from_csc   = true;
    result.source_of_truth_storage = false; // set true on success
    result.true_numeric_source     = false; // §17.2 NOT implemented
    result.valid                   = false;

    if (!sn_num.valid || n <= Index(0)) return result;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.L, n, n)) return result;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.U, n, n)) return result;

    const std::size_t nsup = sn_num.supernode_ptr.size() - 1u;

    // ---- Permutations and equilibration scaling (copy from baseline CSC L/U) ----
    // Convention: row_perm[new_row] = old_row, inv_row_perm[old_row] = new_row
    //             col_perm[new_col] = old_col, inv_col_perm[old_col] = new_col
    result.row_perm     = csc_lu.row_perm;
    result.inv_row_perm = csc_lu.inv_row_perm;
    result.col_perm     = csc_lu.col_perm;
    result.inv_col_perm = csc_lu.inv_col_perm;
    result.Dr = csc_lu.Dr;
    result.Dc = csc_lu.Dc;

    // ---- Allocate supernode descriptors ----
    result.supernodes.resize(nsup);

    // ---- U_segments: off-diagonal U entries (rows < first_col per supernode) ----
    // Layout: seg_ptr[s]..seg_ptr[s+1] indexes the entries for supernode s.
    // Within each segment, entries are stored column-by-column:
    //   first all (row, val) for col_begin, then for col_begin+1, etc.
    // §18.2 storage-native solve PENDING: explicit column sub-partition not stored.
    result.U_segments.seg_ptr.resize(nsup + 1u, Index(0));

    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index col_begin = sn_num.supernode_ptr[s];
        const Index col_end   = sn_num.supernode_ptr[s + 1u];

        result.U_segments.seg_ptr[s] =
            static_cast<Index>(result.U_segments.row_ind.size());

        // Collect off-diagonal U entries (rows < col_begin) column-by-column.
        for (Index j = col_begin; j < col_end; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            for (Index k = csc_lu.U.col_ptr[sj]; k < csc_lu.U.col_ptr[sj + 1u]; ++k) {
                const std::size_t sk = static_cast<std::size_t>(k);
                const Index r = csc_lu.U.row_ind[sk];
                if (r < col_begin) {
                    result.U_segments.row_ind.push_back(r);
                    result.U_segments.values.push_back(csc_lu.U.values[sk]);
                }
            }
        }
    }
    result.U_segments.seg_ptr[nsup] =
        static_cast<Index>(result.U_segments.row_ind.size());

    // ---- Panel values and supernode descriptors ----
    //
    // For each supernode s at columns [col_begin, col_end), width w = col_end - col_begin:
    //
    //   row_indices = sorted union of:
    //     [col_begin, col_begin+1, ..., col_end-1]  (diagonal block rows)
    //     sn_num.l_row_ind[l_row_ptr[s]..l_row_ptr[s+1]) (L rows, includes
    //       within-block rows col_begin+1..col_end-1 and off-diagonal rows >=col_end)
    //   row_count = row_indices.size()
    //   leading_dimension >= row_count
    //
    //   panel_values: column-major, row_count x w block at offset values_offset.
    //   For column c (local, 0..w-1), global_col = col_begin + c:
    //     For each row r_pos in [0, row_count), global_row = row_indices[r_pos]:
    //       - If col_begin <= global_row < col_end:  diagonal block row
    //           b = global_row - col_begin  (local block row)
    //           If b < c:   above diagonal  -> U[global_row, global_col]
    //           If b == c:  diagonal        -> U[global_row, global_col]
    //           If b > c:   below diagonal  -> L[global_row, global_col]
    //       - If global_row >= col_end:  off-diagonal L
    //           -> L[global_row, global_col]
    //     Values of 0 are stored explicitly (dense panel).
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index col_begin = sn_num.supernode_ptr[s];
        const Index col_end   = sn_num.supernode_ptr[s + 1u];
        const Index w         = col_end - col_begin;

        supernode_desc<Index>& desc = result.supernodes[s];
        desc.first_col = col_begin;
        desc.num_cols  = w;

        // Build row_indices: sorted union of diagonal block rows and L rows.
        {
            const Index l_begin = sn_num.l_row_ptr[s];
            const Index l_end   = sn_num.l_row_ptr[s + 1u];
            const std::size_t total =
                static_cast<std::size_t>(w) +
                static_cast<std::size_t>(l_end - l_begin);
            desc.row_indices.clear();
            desc.row_indices.reserve(total);
            // Add diagonal block rows [col_begin, col_end)
            for (Index b = col_begin; b < col_end; ++b) {
                desc.row_indices.push_back(b);
            }
            // Add L rows for this supernode
            for (Index k = l_begin; k < l_end; ++k) {
                desc.row_indices.push_back(
                    sn_num.l_row_ind[static_cast<std::size_t>(k)]);
            }
            // Sort and remove duplicates
            // (some l_row_ind entries may overlap with [col_begin+1, col_end-1])
            std::sort(desc.row_indices.begin(), desc.row_indices.end());
            desc.row_indices.erase(
                std::unique(desc.row_indices.begin(), desc.row_indices.end()),
                desc.row_indices.end());
        }

        const Index row_count      = static_cast<Index>(desc.row_indices.size());
        desc.values_offset   = static_cast<Index>(result.panel_values.size());
        desc.leading_dimension =
            row_count + (add_alignment_padding ? Index(1) : Index(0));

        // U_segment info for this supernode.
        desc.u_segment_start = result.U_segments.seg_ptr[s];
        desc.u_segment_count =
            result.U_segments.seg_ptr[s + 1u] - result.U_segments.seg_ptr[s];

        // SLU-8R.3.1: build per-column entry offset within this supernode's U_segment.
        // u_seg_col_ptr[c] = cumulative U entries for local cols 0..c-1.
        // Mirrors the column-by-column fill order used above in U_segments population.
        {
            desc.u_seg_col_ptr.resize(static_cast<std::size_t>(w) + 1u, Index(0));
            Index col_offset = Index(0);
            for (Index c = Index(0); c < w; ++c) {
                desc.u_seg_col_ptr[static_cast<std::size_t>(c)] = col_offset;
                const Index global_col = col_begin + c;
                const std::size_t sgc = static_cast<std::size_t>(global_col);
                for (Index k = csc_lu.U.col_ptr[sgc]; k < csc_lu.U.col_ptr[sgc + 1u]; ++k) {
                    const Index r = csc_lu.U.row_ind[static_cast<std::size_t>(k)];
                    if (r < col_begin) {
                        ++col_offset;
                    }
                }
            }
            desc.u_seg_col_ptr[static_cast<std::size_t>(w)] = col_offset;
        }

        // Allocate panel (initialized to T(0)).
        const std::size_t panel_size =
            static_cast<std::size_t>(desc.leading_dimension) *
            static_cast<std::size_t>(w);
        result.panel_values.resize(result.panel_values.size() + panel_size, T(0));
        T* panel = &result.panel_values[static_cast<std::size_t>(desc.values_offset)];

        // Fill panel column by column.
        for (Index c = Index(0); c < w; ++c) {
            const Index global_col = col_begin + c;
            const std::size_t sc = static_cast<std::size_t>(c);
            T* col_data =
                panel + sc * static_cast<std::size_t>(desc.leading_dimension);

            // U entries: rows in [col_begin, col_end) with row <= global_col
            // (upper triangular part of diagonal block, including diagonal).
            {
                const std::size_t sgc = static_cast<std::size_t>(global_col);
                for (Index k = csc_lu.U.col_ptr[sgc]; k < csc_lu.U.col_ptr[sgc + 1u]; ++k) {
                    const std::size_t sk = static_cast<std::size_t>(k);
                    const Index r = csc_lu.U.row_ind[sk];
                    if (r >= col_begin && r < col_end) {
                        // Diagonal block row: find position in row_indices.
                        const typename std::vector<Index>::const_iterator it =
                            std::lower_bound(desc.row_indices.begin(),
                                             desc.row_indices.end(), r);
                        if (it != desc.row_indices.end() && *it == r) {
                            const std::size_t r_pos =
                                static_cast<std::size_t>(it - desc.row_indices.begin());
                            col_data[r_pos] = csc_lu.U.values[sk];
                        }
                    }
                }
            }

            // L entries: all rows r > global_col (strictly lower triangular).
            // Covers both within-block L (col_begin < r < col_end) and
            // off-diagonal L (r >= col_end).
            {
                const std::size_t sgc = static_cast<std::size_t>(global_col);
                for (Index k = csc_lu.L.col_ptr[sgc]; k < csc_lu.L.col_ptr[sgc + 1u]; ++k) {
                    const std::size_t sk = static_cast<std::size_t>(k);
                    const Index r = csc_lu.L.row_ind[sk];
                    // L is strictly lower triangular: r > global_col always holds.
                    // Find position in row_indices.
                    const typename std::vector<Index>::const_iterator it =
                        std::lower_bound(desc.row_indices.begin(),
                                         desc.row_indices.end(), r);
                    if (it != desc.row_indices.end() && *it == r) {
                        const std::size_t r_pos =
                            static_cast<std::size_t>(it - desc.row_indices.begin());
                        col_data[r_pos] = csc_lu.L.values[sk];
                    }
                }
            }
        }
    }

    // All supernode invariants met on success.
    result.source_of_truth_storage = true;
    result.valid                   = true;
    return result;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_SUPERNODAL_STORAGE_BOOTSTRAP_IMPL_HPP
