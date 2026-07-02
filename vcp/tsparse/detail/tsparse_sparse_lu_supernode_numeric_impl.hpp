// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-10 supernode numeric metadata builder -- internal implementation.
//
// Builds sparse_lu_supernode_numeric<T, Index> from actual CSC L/U factors
// produced by the baseline sparse GP factorization.  This is a CSC-backed
// prototype: the numeric source of truth remains the baseline CSC L/U storage.
//
// IMPORTANT: All row patterns and values are derived from the actual emitted
// CSC L/U factors.  Symbolic panel_row_ind is NEVER used as the final numeric
// row structure.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - all storage type definitions (csc_storage, baseline_lu_storage, etc.)
//   - sparse_lu_supernode_symbolic<Index>
//   - sparse_lu_supernode_numeric<T, Index>
//   - sparse_lu_is_valid_supernode_symbolic
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse.hpp>  (umbrella)

#ifndef VCP_TSPARSE_SPARSE_LU_SUPERNODE_NUMERIC_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SUPERNODE_NUMERIC_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// build_supernode_numeric_from_csc
//
// SLU-10 prototype: Build supernode numeric block descriptors from actual
// CSC L/U factors (baseline_lu_storage).
//
// For each supernode s with columns [col_begin, col_end):
//   l_col_start[s]  = L.col_ptr[col_begin]
//   u_col_start[s]  = U.col_ptr[col_begin]
//   l_row_ind slice = sorted unique row indices from actual L CSC columns
//   u_row_ind slice = sorted unique row indices from actual U CSC columns
//   diag_block      = column-major (width x width) block from actual U diagonal
//
// Preconditions (returns invalid result if violated):
//   n >= 0
//   sym_info is valid per sparse_lu_is_valid_supernode_symbolic
//   csc_lu.L.col_ptr.size() == n+1
//   csc_lu.U.col_ptr.size() == n+1
//
// Invariants guaranteed on success:
//   result.valid == true
//   result.supernode_ptr  == sym_info.supernode_ptr
//   result.column_to_supernode.size() == n
//   result.l_row_ptr.size() == nsup + 1, monotone, front==0
//   result.u_row_ptr.size() == nsup + 1, monotone, front==0
//   result.l_row_ind per supernode: sorted unique, entries in [0,n)
//   result.u_row_ind per supernode: sorted unique, entries in [0,n)
//   result.diag_block_ptr.size() == nsup + 1
//   each block size == width * width (column-major)
//   values derived from actual U columns (not symbolic pattern)
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_supernode_numeric<T, Index>
build_supernode_numeric_from_csc(
    Index n,
    const baseline_lu_storage<T, Index>& csc_lu,
    const sparse_lu_supernode_symbolic<Index>& sym_info)
{
    static_assert(std::is_signed<Index>::value,
                  "build_supernode_numeric_from_csc: Index must be signed");

    sparse_lu_supernode_numeric<T, Index> result;
    result.valid = false;

    if (n < Index(0)) return result;

    // SLU-10.1: validate CSC L/U storage before any row/value access.
    // Checks col_ptr size, col_ptr[0]==0, monotone, row_ind/values size, row range [0,n).
    if (!sparse_lu_is_valid_csc_storage(csc_lu.L, n, n)) return result;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.U, n, n)) return result;

    // Validate symbolic supernode metadata
    if (!sym_info.valid ||
        !sparse_lu_is_valid_supernode_symbolic(n, sym_info)) {
        return result;
    }

    const std::size_t un   = static_cast<std::size_t>(n);
    const std::size_t nsup = sym_info.supernode_ptr.size() - 1u;

    // Copy partition from symbolic metadata
    result.supernode_ptr       = sym_info.supernode_ptr;
    result.column_to_supernode = sym_info.column_to_supernode;

    result.l_col_start.resize(nsup, Index(0));
    result.u_col_start.resize(nsup, Index(0));
    result.l_row_ptr.resize(nsup + 1u, Index(0));
    result.u_row_ptr.resize(nsup + 1u, Index(0));
    result.diag_block_ptr.resize(nsup + 1u, Index(0));

    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index col_begin = sym_info.supernode_ptr[s];
        const Index col_end   = sym_info.supernode_ptr[s + 1u];
        const std::size_t scol_begin = static_cast<std::size_t>(col_begin);

        // Column start in baseline CSC for first column of this supernode
        result.l_col_start[s] = csc_lu.L.col_ptr[scol_begin];
        result.u_col_start[s] = csc_lu.U.col_ptr[scol_begin];

        // --- L row indices for this supernode ---
        // Derived from actual emitted CSC L row_ind; NOT from symbolic panel_row_ind.
        result.l_row_ptr[s] = static_cast<Index>(result.l_row_ind.size());
        {
            std::vector<Index> tmp_rows;
            for (Index j = col_begin; j < col_end; ++j) {
                const std::size_t sj = static_cast<std::size_t>(j);
                const Index kbegin = csc_lu.L.col_ptr[sj];
                const Index kend   = csc_lu.L.col_ptr[sj + 1u];
                for (Index k = kbegin; k < kend; ++k) {
                    tmp_rows.push_back(
                        csc_lu.L.row_ind[static_cast<std::size_t>(k)]);
                }
            }
            std::sort(tmp_rows.begin(), tmp_rows.end());
            tmp_rows.erase(
                std::unique(tmp_rows.begin(), tmp_rows.end()), tmp_rows.end());
            result.l_row_ind.insert(
                result.l_row_ind.end(), tmp_rows.begin(), tmp_rows.end());
        }

        // --- U row indices for this supernode ---
        // Derived from actual emitted CSC U row_ind; NOT from symbolic panel_row_ind.
        result.u_row_ptr[s] = static_cast<Index>(result.u_row_ind.size());
        {
            std::vector<Index> tmp_rows;
            for (Index j = col_begin; j < col_end; ++j) {
                const std::size_t sj = static_cast<std::size_t>(j);
                const Index kbegin = csc_lu.U.col_ptr[sj];
                const Index kend   = csc_lu.U.col_ptr[sj + 1u];
                for (Index k = kbegin; k < kend; ++k) {
                    tmp_rows.push_back(
                        csc_lu.U.row_ind[static_cast<std::size_t>(k)]);
                }
            }
            std::sort(tmp_rows.begin(), tmp_rows.end());
            tmp_rows.erase(
                std::unique(tmp_rows.begin(), tmp_rows.end()), tmp_rows.end());
            result.u_row_ind.insert(
                result.u_row_ind.end(), tmp_rows.begin(), tmp_rows.end());
        }

        // --- Diagonal block values from actual U (column-major, width x width) ---
        // Each entry is looked up in the actual CSC U column.
        // Required diagonal U(j,j): missing entry is invalid (not silently zero).
        // Non-diagonal in-block entries: absent entry is treated as zero.
        result.diag_block_ptr[s] = static_cast<Index>(result.diag_block_values.size());
        {
            const Index width = col_end - col_begin;
            for (Index j = col_begin; j < col_end; ++j) {
                const std::size_t sj = static_cast<std::size_t>(j);
                const Index kbegin = csc_lu.U.col_ptr[sj];
                const Index kend   = csc_lu.U.col_ptr[sj + 1u];
                for (Index i = col_begin; i < col_end; ++i) {
                    bool found = false;
                    T val = T(0);
                    for (Index k = kbegin; k < kend; ++k) {
                        const std::size_t sk = static_cast<std::size_t>(k);
                        if (csc_lu.U.row_ind[sk] == i) {
                            val = csc_lu.U.values[sk];
                            found = true;
                            break;
                        }
                    }
                    // SLU-11.1: required diagonal entry must be present in U CSC.
                    if (i == j && !found) {
                        result = sparse_lu_supernode_numeric<T, Index>();
                        return result; // valid == false
                    }
                    result.diag_block_values.push_back(val);
                }
            }
            (void)width;
        }
    }

    result.l_row_ptr[nsup] = static_cast<Index>(result.l_row_ind.size());
    result.u_row_ptr[nsup] = static_cast<Index>(result.u_row_ind.size());
    result.diag_block_ptr[nsup] = static_cast<Index>(result.diag_block_values.size());

    // SLU-10.1: postcondition check — only mark valid if all invariants hold.
    // Tentatively set valid = true so the public validator can inspect the payload.
    result.valid = true;
    if (!sparse_lu_is_valid_supernode_numeric(n, sym_info, result)) {
        result = sparse_lu_supernode_numeric<T, Index>();
        return result;
    }
    return result;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_SUPERNODE_NUMERIC_IMPL_HPP
