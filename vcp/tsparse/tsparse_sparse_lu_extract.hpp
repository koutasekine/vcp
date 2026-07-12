// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// LUX-0: extraction of the L / U / p / q factors from a baseline-CSC sparse
// LU factorization object, in the COMMON OUTPUT CONVENTION (design
// lux_design_v0 SSC, shared with spumar):
//
//     P . A . Q = L . U
//     L: unit lower triangular (unit diagonal EXPLICITLY stored, value T(1)),
//     U: upper triangular (explicit diagonal),
//     vector form (new->old):  A(p, q) = L . U,  i.e.
//         (P A Q)(i, j) = A(p[i], q[j]),
//     permutation matrices:  P(k, p[k]) = 1,  Q(q[k], k) = 1.
//
// Mapping from the SLU internal convention (P-1, frozen decision; see
// tsparse_sparse_lu_inverse_row.hpp header for the operator convention):
//     solve pipeline (equilibration disabled, Dr/Dc empty):
//         x = Q_op U^{-1} L^{-1} P_op b,
//         (P_op v)[new] = v[row_perm[new]],  (Q_op z)[col_perm[new]] = z[new]
//     hence  A^{-1} = Q_op U^{-1} L^{-1} P_op  =>  P_op A Q_op = L U,  and
//     (P_op A Q_op)(i, j) = A(row_perm[i], col_perm[j]),  therefore
//         p = row_perm,  q = col_perm   (both copied new->old, unchanged).
//     The correctness of this mapping is MACHINE-VERIFIED by the A(p,q) = LU
//     residual gate (G-X0.1); it must never be adjusted by guesswork.
//
// Scope (P-3/P-4/P-5, frozen decisions):
//   - baseline_csc storage only; supernodal factors -> unsupported_storage.
//   - equilibrated factors (Dr or Dc non-empty; produced by
//     equilibration=true or pivoting=static_mc64) -> unsupported_options
//     (the scaled extended form P Dr A Dc Q = L U is a future addition).
//   - status priority: invalid_factorization > unsupported_options >
//     unsupported_storage > extraction.
//   - output CSC columns are sorted ascending by row index; exact zeros in
//     the baseline storage are NOT copied (spmats invariant: explicit zeros
//     are never stored).  The zero-drop test uses the same  v == T(0)  rule
//     as the spmats validators and the LDL D assembly (certified-zero drop);
//     it is the single sanctioned T comparison of this layer (P-3 overrides
//     the general B-2 no-T-comparison rule for exactly this drop).
//   - otherwise the extraction performs ONLY index manipulation and value
//     copies (B-2).  The unit-diagonal constant T(1) of L is the single
//     sanctioned T construction.
//
// Non-throwing contract: every runtime failure is reported through the
// status field.  Malformed factor storage (a vcp::error thrown by the
// structural validator) and any std::exception map to internal_error.

#pragma once

#ifndef VCP_TSPARSE_SPARSE_LU_EXTRACT_HPP
#define VCP_TSPARSE_SPARSE_LU_EXTRACT_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <type_traits>
#include <vector>

#include <vcp/tsparse/tsparse_sparse_lu.hpp>

namespace vcp {

// ===========================================================================
// Status / result types (design lux_design_v0 SS1)
// ===========================================================================

enum class sparse_lu_extract_status {
    success,
    invalid_factorization,   // fac.valid() == false (the factorization failed)
    unsupported_storage,     // non-baseline_csc storage (supernodal; v0 scope)
    unsupported_options,     // equilibrated factor (Dr or Dc non-identity)
    internal_error
};

inline const char* sparse_lu_extract_status_to_string(sparse_lu_extract_status s) {
    switch (s) {
    case sparse_lu_extract_status::success:               return "success";
    case sparse_lu_extract_status::invalid_factorization: return "invalid_factorization";
    case sparse_lu_extract_status::unsupported_storage:   return "unsupported_storage";
    case sparse_lu_extract_status::unsupported_options:   return "unsupported_options";
    case sparse_lu_extract_status::internal_error:        return "internal_error";
    }
    return "unknown";
}

template <class T, class Index>
struct sparse_lu_extracted {
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_extract: Index must be signed");

    sparse_lu_extract_status status;

    // L: CSC, unit diagonal explicitly stored (T(1)), columns sorted ascending.
    std::vector<Index> L_col_ptr, L_row_ind;
    std::vector<T>     L_val;
    // U: CSC, upper triangular, explicit diagonal, columns sorted ascending.
    std::vector<Index> U_col_ptr, U_row_ind;
    std::vector<T>     U_val;

    std::vector<Index> p, q;   // SSC new->old convention: A(p,q) = L U

    Index nnz_L, nnz_U;

    sparse_lu_extracted()
        : status(sparse_lu_extract_status::internal_error),
          nnz_L(Index(0)), nnz_U(Index(0)) {}
};

namespace sparse_lu_extract_detail {

// Copy one factor from the baseline CSC storage into (col_ptr, row_ind, val),
// sorting each column ascending by row index and dropping exact zeros
// (v == T(0), the spmats invariant rule).  Index-only comparator: the sort
// key is the row index, never a T value.
//
// unit_lower == true (L): stored diagonal entries are IGNORED (unit-diagonal
// contract of the SLU baseline storage) and a fresh explicit T(1) diagonal is
// emitted as the first entry of every column; any stored row < col entry is
// a malformed factor -> returns false.
// unit_lower == false (U): entries are copied verbatim; any stored row > col
// entry is malformed -> returns false.  (A missing / zero U diagonal is NOT
// rejected here: extraction is a faithful copy, and the factor-consuming
// layer re-gates divisions itself.)
template <class T, class Index>
bool copy_factor_sorted_(const csc_storage<T, Index>& M,
                         const Index n,
                         const bool unit_lower,
                         std::vector<Index>& col_ptr,
                         std::vector<Index>& row_ind,
                         std::vector<T>& val)
{
    const std::size_t un = static_cast<std::size_t>(n);
    col_ptr.assign(un + 1u, Index(0));
    row_ind.clear();
    val.clear();
    row_ind.reserve(M.row_ind.size() + (unit_lower ? un : 0u));
    val.reserve(M.row_ind.size() + (unit_lower ? un : 0u));

    std::vector<Index> pos;   // per-column scratch: positions into M storage
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);

        pos.clear();
        for (Index k = M.col_ptr[sj]; k < M.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index r = M.row_ind[sk];
            if (unit_lower) {
                if (r < j) return false;          // upper entry in L: malformed
                if (r == j) continue;             // stored diagonal: unit contract, ignore
            } else {
                if (r > j) return false;          // lower entry in U: malformed
            }
            if (M.values[sk] == T(0)) continue;   // exact zero: not stored (P-3)
            pos.push_back(k);
        }

        // sort by row index only (integer comparator; no T comparison)
        std::sort(pos.begin(), pos.end(),
                  [&M](const Index a, const Index b) {
                      return M.row_ind[static_cast<std::size_t>(a)]
                           < M.row_ind[static_cast<std::size_t>(b)];
                  });

        if (unit_lower) {
            row_ind.push_back(j);                 // explicit unit diagonal first
            val.push_back(T(1));
        }
        for (std::size_t s = 0; s < pos.size(); ++s) {
            const std::size_t sk = static_cast<std::size_t>(pos[s]);
            row_ind.push_back(M.row_ind[sk]);
            val.push_back(M.values[sk]);
        }
        col_ptr[sj + 1u] = static_cast<Index>(row_ind.size());
    }
    return true;
}

} // namespace sparse_lu_extract_detail

// ===========================================================================
// sparse_lu_extract_factors
// Reads a CONSTRUCTED factorization object (the SLU call sites are
// unchanged; this function only reads) and returns the SSC-convention
// factors.  Non-throwing: all outcomes are status values.
// ===========================================================================
template <class T, class Index>
sparse_lu_extracted<T, Index>
sparse_lu_extract_factors(const sparse_lu_factorization<T, Index>& fac)
{
    sparse_lu_extracted<T, Index> out;

    try {
        // ---- entry gates, P-4 priority order ----
        if (!fac.valid()) {
            out.status = sparse_lu_extract_status::invalid_factorization;
            return out;
        }

        const Index n = fac.info().n;
        const baseline_lu_storage<T, Index>& lu = fac.baseline_storage();

        if (!lu.Dr.empty() || !lu.Dc.empty()) {
            out.status = sparse_lu_extract_status::unsupported_options;
            return out;
        }

        if (fac.storage_kind() != sparse_lu_storage_kind::baseline_csc) {
            out.status = sparse_lu_extract_status::unsupported_storage;
            return out;
        }

        // Structural validation (index-level only; same validator as the
        // solve path).  A malformed storage throws vcp::error, mapped to
        // internal_error by the net below.
        sparse_lu_detail::validate_baseline_lu_storage_for_solve(lu, n);

        // ---- extraction (P-3) ----
        if (!sparse_lu_extract_detail::copy_factor_sorted_<T, Index>(
                lu.L, n, /*unit_lower=*/true,
                out.L_col_ptr, out.L_row_ind, out.L_val)) {
            out.status = sparse_lu_extract_status::internal_error;
            return out;
        }
        if (!sparse_lu_extract_detail::copy_factor_sorted_<T, Index>(
                lu.U, n, /*unit_lower=*/false,
                out.U_col_ptr, out.U_row_ind, out.U_val)) {
            out = sparse_lu_extracted<T, Index>();
            out.status = sparse_lu_extract_status::internal_error;
            return out;
        }

        // ---- permutations (P-1): p = row_perm, q = col_perm, new->old;
        // empty internal permutation means identity -> materialize it. ----
        const std::size_t un = static_cast<std::size_t>(n);
        if (lu.row_perm.empty()) {
            out.p.resize(un);
            for (Index k = Index(0); k < n; ++k) {
                out.p[static_cast<std::size_t>(k)] = k;
            }
        } else {
            out.p = lu.row_perm;
        }
        if (lu.col_perm.empty()) {
            out.q.resize(un);
            for (Index k = Index(0); k < n; ++k) {
                out.q[static_cast<std::size_t>(k)] = k;
            }
        } else {
            out.q = lu.col_perm;
        }

        out.nnz_L = static_cast<Index>(out.L_row_ind.size());
        out.nnz_U = static_cast<Index>(out.U_row_ind.size());
        out.status = sparse_lu_extract_status::success;
        return out;
    } catch (const std::exception&) {
        out = sparse_lu_extracted<T, Index>();
        out.status = sparse_lu_extract_status::internal_error;
        return out;
    }
}

} // namespace vcp

#endif // VCP_TSPARSE_SPARSE_LU_EXTRACT_HPP
