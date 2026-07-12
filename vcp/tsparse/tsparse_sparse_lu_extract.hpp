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
// Scope (P-3/P-4/P-5 frozen decisions; P-11..P-14 added by LUX-3):
//   - baseline_csc AND supernodal storage kinds (LUX-3, P-11).  Unknown
//     storage kinds -> unsupported_storage.
//   - supernodal dispatch mirrors fac.solve() (SLU-8R.5.5.1 priority):
//       * accepted A_eff-origin storage (valid && source_of_truth &&
//         true_numeric_source && !supernodes.empty(), the same predicate as
//         can_use_supernodal_storage_solve): flatten the supernodal panels /
//         U_segments (survey: LUX-3_S0_survey.md SS1.4 mapping).
//       * transitional storage (true_numeric_source == false): the panel
//         values are NOT the numeric source (R.5.5.2 in-place may-modify
//         contract); extraction falls back to the baseline CSC factors,
//         exactly like the solve dispatch.  storage_kind_extracted records
//         which route produced the output.
//   - equilibrated factors (Dr or Dc non-empty; produced by
//     equilibration=true or pivoting=static_mc64 on the GP path) ->
//     unsupported_options (the scaled extended form P Dr A Dc Q = L U is a
//     future addition).  The MF7 native-MC64 path is matching-only (empty
//     Dr/Dc) and extracts fine.
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
    unsupported_storage,     // unknown storage kind (supernodal supported since LUX-3)
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

    // LUX-3 (P-11): storage kind the factors were actually read from.
    // Meaningful only when status == success.  For a supernodal factor this
    // is `supernodal` when the accepted A_eff-origin storage was flattened,
    // and `baseline_csc` when the transitional CSC fallback route was taken
    // (mirroring the fac.solve() dispatch).  Existing fields are unchanged.
    sparse_lu_storage_kind storage_kind_extracted;

    sparse_lu_extracted()
        : status(sparse_lu_extract_status::internal_error),
          nnz_L(Index(0)), nnz_U(Index(0)),
          storage_kind_extracted(sparse_lu_storage_kind::baseline_csc) {}
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

// ===========================================================================
// LUX-3 (P-11..P-14): supernodal storage extraction
// ===========================================================================

// Index-level structural validation of a supernodal_lu_storage before
// flattening (no T arithmetic; mirrors the invariants the SS18.2 native solve
// relies on and the self-symbolic builder self-checks, SLU-MF2 SS12.5):
//   - supernodes form a contiguous partition of [0, n)
//   - row_indices: size >= num_cols, strictly ascending, row_indices[c] ==
//     first_col + c for the pivot block, all < n
//   - panel bounds: leading_dimension >= row_count and the panel slab fits
//     inside panel_values (rows in [row_count, ld) are alignment padding and
//     are never read)
//   - U_segments: per-supernode u_seg_col_ptr present (size num_cols+1,
//     monotone, [0] == 0, terminal == u_segment_count) whenever entries
//     exist; entry rows in [0, first_col); segment range inside the arrays
//   - permutations: size n (or empty == identity) and bijective
// Returns false on any violation (mapped to internal_error by the caller,
// consistent with the malformed-factor contract of the baseline route).
template <class T, class Index>
bool validate_supernodal_storage_for_extract_(
    const supernodal_lu_storage<T, Index>& st,
    const Index n)
{
    if (n < Index(0)) return false;
    const std::size_t un = static_cast<std::size_t>(n);

    if (st.U_segments.values.size() != st.U_segments.row_ind.size()) return false;

    Index next_col = Index(0);
    for (std::size_t s = 0; s < st.supernodes.size(); ++s) {
        const supernode_desc<Index>& d = st.supernodes[s];
        const Index fc = d.first_col;
        const Index w  = d.num_cols;
        if (fc != next_col || w <= Index(0)) return false;
        next_col = fc + w;
        if (next_col > n) return false;

        const std::size_t sw = static_cast<std::size_t>(w);
        const std::size_t rcnt = d.row_indices.size();
        if (rcnt < sw) return false;
        for (std::size_t c = 0; c < sw; ++c) {
            if (d.row_indices[c] != fc + static_cast<Index>(c)) return false;
        }
        for (std::size_t r = 1; r < rcnt; ++r) {
            if (!(d.row_indices[r - 1u] < d.row_indices[r])) return false;
        }
        if (d.row_indices[rcnt - 1u] >= n) return false;

        if (d.leading_dimension < static_cast<Index>(rcnt)) return false;
        if (d.values_offset < Index(0)) return false;
        const std::size_t slab =
            static_cast<std::size_t>(d.leading_dimension) * sw;
        if (static_cast<std::size_t>(d.values_offset) + slab >
            st.panel_values.size()) return false;

        if (d.u_segment_count < Index(0) || d.u_segment_start < Index(0)) return false;
        if (d.u_segment_count > Index(0)) {
            if (d.u_seg_col_ptr.size() != sw + 1u) return false;
            if (d.u_seg_col_ptr[0] != Index(0)) return false;
            for (std::size_t c = 1; c <= sw; ++c) {
                if (d.u_seg_col_ptr[c] < d.u_seg_col_ptr[c - 1u]) return false;
            }
            if (d.u_seg_col_ptr[sw] != d.u_segment_count) return false;
            const std::size_t seg_end =
                static_cast<std::size_t>(d.u_segment_start) +
                static_cast<std::size_t>(d.u_segment_count);
            if (seg_end > st.U_segments.row_ind.size()) return false;
            for (Index k = d.u_segment_start;
                 k < d.u_segment_start + d.u_segment_count; ++k) {
                const Index r = st.U_segments.row_ind[static_cast<std::size_t>(k)];
                if (r < Index(0) || r >= fc) return false;
            }
        }
    }
    if (next_col != n) return false;

    // permutations: size n or empty; bijective when present (index-only)
    if (!st.row_perm.empty()) {
        if (st.row_perm.size() != un) return false;
        std::vector<char> seen(un, 0);
        for (std::size_t k = 0; k < un; ++k) {
            const Index v = st.row_perm[k];
            if (v < Index(0) || v >= n) return false;
            if (seen[static_cast<std::size_t>(v)]) return false;
            seen[static_cast<std::size_t>(v)] = 1;
        }
    }
    if (!st.col_perm.empty()) {
        if (st.col_perm.size() != un) return false;
        std::vector<char> seen(un, 0);
        for (std::size_t k = 0; k < un; ++k) {
            const Index v = st.col_perm[k];
            if (v < Index(0) || v >= n) return false;
            if (seen[static_cast<std::size_t>(v)]) return false;
            seen[static_cast<std::size_t>(v)] = 1;
        }
    }
    return true;
}

// Copy a new->old permutation (empty == identity is materialized).
template <class Index>
void copy_or_identity_perm_(const std::vector<Index>& perm,
                            const Index n,
                            std::vector<Index>& out)
{
    if (perm.empty()) {
        out.resize(static_cast<std::size_t>(n));
        for (Index k = Index(0); k < n; ++k) {
            out[static_cast<std::size_t>(k)] = k;
        }
    } else {
        out = perm;
    }
}

// Transitional-storage fallback: extract the baseline CSC factors, mirroring
// the fac.solve() CSC fallback (the panel values of a transitional supernodal
// storage are NOT the numeric source, R.5.5.2 may-modify contract).  Same
// gates and steps as the baseline route of sparse_lu_extract_factors (P-4
// order after the storage dispatch); storage_kind_extracted records the
// baseline_csc origin.  May throw (vcp::error from the structural validator);
// the caller's net maps it to internal_error.
template <class T, class Index>
sparse_lu_extracted<T, Index>
extract_baseline_csc_factors_(const baseline_lu_storage<T, Index>& lu,
                              const Index n)
{
    sparse_lu_extracted<T, Index> out;

    if (!lu.Dr.empty() || !lu.Dc.empty()) {
        out.status = sparse_lu_extract_status::unsupported_options;
        return out;
    }

    sparse_lu_detail::validate_baseline_lu_storage_for_solve(lu, n);

    if (!copy_factor_sorted_<T, Index>(
            lu.L, n, /*unit_lower=*/true,
            out.L_col_ptr, out.L_row_ind, out.L_val)) {
        out.status = sparse_lu_extract_status::internal_error;
        return out;
    }
    if (!copy_factor_sorted_<T, Index>(
            lu.U, n, /*unit_lower=*/false,
            out.U_col_ptr, out.U_row_ind, out.U_val)) {
        out = sparse_lu_extracted<T, Index>();
        out.status = sparse_lu_extract_status::internal_error;
        return out;
    }

    copy_or_identity_perm_(lu.row_perm, n, out.p);
    copy_or_identity_perm_(lu.col_perm, n, out.q);

    out.nnz_L = static_cast<Index>(out.L_row_ind.size());
    out.nnz_U = static_cast<Index>(out.U_row_ind.size());
    out.storage_kind_extracted = sparse_lu_storage_kind::baseline_csc;
    out.status = sparse_lu_extract_status::success;
    return out;
}

// Supernodal-storage extraction (LUX-3 core).  Flattens the accepted
// A_eff-origin supernodal storage into the SSC-convention CSC factors:
//   panel diag block col c: rows r < c -> U, r == c -> U diagonal,
//                           r > c -> L; off-diagonal panel rows -> L;
//   U_segments (per u_seg_col_ptr) -> U rows < first_col;
//   L unit diagonal: fresh explicit T(1) (not stored in the panel);
//   alignment padding rows (>= row_count) are never read;
//   exact zeros (dense-panel fill / relaxed-amalgamation padding /
//   dense-front-model superset) are dropped by the sanctioned  v == T(0)
//   rule (P-12 == P-3).  All other work is index manipulation and value
//   copies (P-14).
// L columns are emitted directly ascending (row_indices is strictly
// ascending); U columns sort the U_segments part by row index (ascending
// order there is NOT code-guaranteed; survey SS1.6) with an index-only
// comparator, then append the diag-block part (all rows >= first_col).
// May throw only through std::vector; caller's net maps to internal_error.
template <class T, class Index>
sparse_lu_extracted<T, Index>
extract_supernodal_factors_(const supernodal_lu_storage<T, Index>& st,
                            const Index n)
{
    sparse_lu_extracted<T, Index> out;

    if (!st.Dr.empty() || !st.Dc.empty()) {
        out.status = sparse_lu_extract_status::unsupported_options;
        return out;
    }

    if (!validate_supernodal_storage_for_extract_(st, n)) {
        out.status = sparse_lu_extract_status::internal_error;
        return out;
    }

    const std::size_t un = static_cast<std::size_t>(n);
    out.L_col_ptr.assign(un + 1u, Index(0));
    out.U_col_ptr.assign(un + 1u, Index(0));

    std::vector<Index> useg_rows;   // per-column scratch (U_segments part)
    std::vector<T>     useg_vals;
    std::vector<std::size_t> order;

    for (std::size_t s = 0; s < st.supernodes.size(); ++s) {
        const supernode_desc<Index>& d = st.supernodes[s];
        const Index fc = d.first_col;
        const std::size_t w    = static_cast<std::size_t>(d.num_cols);
        const std::size_t rcnt = d.row_indices.size();
        const std::size_t ld   = static_cast<std::size_t>(d.leading_dimension);
        const std::size_t off  = static_cast<std::size_t>(d.values_offset);
        const T* panel = st.panel_values.data() + off;

        for (std::size_t c = 0; c < w; ++c) {
            const std::size_t gc = static_cast<std::size_t>(fc) + c;
            const T* col = panel + c * ld;

            // ---- L column gc: explicit unit diagonal first, then the
            // strictly-lower diag-block rows, then the off-diagonal rows
            // (already ascending; padding rows >= rcnt never read) ----
            out.L_row_ind.push_back(static_cast<Index>(gc));
            out.L_val.push_back(T(1));
            for (std::size_t r = c + 1u; r < rcnt; ++r) {
                if (col[r] == T(0)) continue;      // P-12 drop
                out.L_row_ind.push_back(d.row_indices[r]);
                out.L_val.push_back(col[r]);
            }
            out.L_col_ptr[gc + 1u] = static_cast<Index>(out.L_row_ind.size());

            // ---- U column gc: U_segments part (rows < fc, sorted by an
            // index-only comparator), then diag-block rows 0..c ----
            useg_rows.clear();
            useg_vals.clear();
            if (d.u_segment_count > Index(0)) {
                const std::size_t base =
                    static_cast<std::size_t>(d.u_segment_start);
                const std::size_t kb =
                    base + static_cast<std::size_t>(d.u_seg_col_ptr[c]);
                const std::size_t ke =
                    base + static_cast<std::size_t>(d.u_seg_col_ptr[c + 1u]);
                for (std::size_t k = kb; k < ke; ++k) {
                    if (st.U_segments.values[k] == T(0)) continue;   // P-12 drop
                    useg_rows.push_back(st.U_segments.row_ind[k]);
                    useg_vals.push_back(st.U_segments.values[k]);
                }
            }
            order.resize(useg_rows.size());
            for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
            std::sort(order.begin(), order.end(),
                      [&useg_rows](const std::size_t a, const std::size_t b) {
                          return useg_rows[a] < useg_rows[b];
                      });
            Index prev_row = Index(-1);
            for (std::size_t i = 0; i < order.size(); ++i) {
                const Index r = useg_rows[order[i]];
                if (r == prev_row) {               // duplicate row: malformed
                    out = sparse_lu_extracted<T, Index>();
                    out.status = sparse_lu_extract_status::internal_error;
                    return out;
                }
                prev_row = r;
                out.U_row_ind.push_back(r);
                out.U_val.push_back(useg_vals[order[i]]);
            }
            for (std::size_t r = 0; r <= c; ++r) {
                if (col[r] == T(0)) continue;      // P-12 drop
                out.U_row_ind.push_back(d.row_indices[r]);
                out.U_val.push_back(col[r]);
            }
            out.U_col_ptr[gc + 1u] = static_cast<Index>(out.U_row_ind.size());
        }
    }

    copy_or_identity_perm_(st.row_perm, n, out.p);   // P-13 == P-1
    copy_or_identity_perm_(st.col_perm, n, out.q);

    out.nnz_L = static_cast<Index>(out.L_row_ind.size());
    out.nnz_U = static_cast<Index>(out.U_row_ind.size());
    out.storage_kind_extracted = sparse_lu_storage_kind::supernodal;
    out.status = sparse_lu_extract_status::success;
    return out;
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

        // ---- P-11 dispatch (LUX-3): supernodal storage kind ----
        // Mirrors the fac.solve() dispatch: the accepted A_eff-origin
        // storage (same predicate as can_use_supernodal_storage_solve) is
        // flattened; transitional storage (true_numeric_source == false)
        // falls back to the baseline CSC factors that solve() actually
        // consumes.  The baseline route below is unchanged.
        if (fac.storage_kind() == sparse_lu_storage_kind::supernodal) {
            const supernodal_lu_storage<T, Index>& st = fac.supernodal_storage();
            if (st.valid && st.source_of_truth_storage &&
                st.true_numeric_source && !st.supernodes.empty()) {
                return sparse_lu_extract_detail::extract_supernodal_factors_<T, Index>(
                    st, fac.info().n);
            }
            return sparse_lu_extract_detail::extract_baseline_csc_factors_<T, Index>(
                fac.baseline_storage(), fac.info().n);
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
