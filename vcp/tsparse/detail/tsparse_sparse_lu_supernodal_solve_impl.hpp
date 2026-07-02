// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-8R.5: §18.2 storage-native supernodal solve -- internal implementation.
//
// ALGORITHM (§18.2 solve order):
//   Step 1: rhs1 = Dr * b                           (row scaling)
//   Step 2: rhs2[new_i] = rhs1[row_perm[new_i]]    (row permutation)
//   Step 3: z1 = L^{-1} * rhs2                     (supernodal L forward solve)
//   Step 4: z2 = U^{-1} * z1                       (supernodal U backward solve)
//   Step 5: x_orig[col_perm[new_j]] = z2[new_j]    (column permutation exit)
//   Step 6: x = Dc * x_orig                        (column scaling)
//
// L SOLVE (forward supernode order):
//   For each supernode s (s = 0 .. nsup-1):
//     Gather: work_block[r] = work[desc.row_indices[r]] for r in [0, R)
//     trsm('L','L','N','U', w, 1, ...): unit-lower diagonal block
//     gemv(-1, L_offdiag, work_block[0..w-1], work_block[w..R-1]):  off-diag L update
//     Scatter: work[desc.row_indices[r]] = work_block[r]
//
// U SOLVE (reverse supernode order):
//   For each supernode s (s = nsup-1 .. 0):
//     Gather: work_block[c] = work[desc.row_indices[c]] for c in [0, w)
//     trsm('L','U','N','N', w, 1, ...): non-unit upper diagonal block
//     Apply U_segments: work[row_r] -= val * work_block[c] for each (row_r,val) in col c
//     Scatter: work[desc.row_indices[c]] = work_block[c]
//
// SCOPE:
//   dense kernel adapter (sparse_lu_dense_kernel<T>) used for ALL dense ops.
//   No direct vcp::tblas / vcp::tlapack calls.
//   No getrf / getrs (U diagonal is pre-factored in panel_values).
//   Gate 6: PENDING (baseline CSC still numeric source-of-truth for §17.2).
//   Solve parity guaranteed only for matrices where §17.2(A) is a no-op.
//   SLU-8 full conformance: NOT claimed.
//   issue_SLU8_contract_violation.md: OPEN.
//
// supernodal_storage_solve_stats is defined in tsparse_sparse_lu.hpp
// (in namespace sparse_lu_detail, before sparse_lu_factorization class body).
// Do NOT redefine it here.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - supernodal_lu_storage<T, Index> defined (with panel_values, U_segments)
//   - supernode_desc<Index> defined (with row_indices, u_seg_col_ptr)
//   - sparse_lu_dense_kernel<T> adapter defined and implemented
//   - within_panel_factor_status enum defined
//   - sparse_lu_factorization<T, Index> class defined
//   - supernodal_storage_solve_stats defined (in sparse_lu_detail)
//   - <chrono>, <string>, <vector> are available
// It has no "namespace vcp { }" wrapper.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//
// Gate status after SLU-8R.5:
//   Gate 1: PASS  Gate 2: PARTIAL  Gate 3: PENDING-REQUIRED
//   Gate 4: PASS  Gate 5: PENDING-REQUIRED  Gate 6: PENDING

#ifndef VCP_TSPARSE_SPARSE_LU_SUPERNODAL_SOLVE_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SUPERNODAL_SOLVE_IMPL_HPP

#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

#include <vcp/error.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// can_use_supernodal_storage_solve
//
// Returns true iff storage-native §18.2 solve can be used.
//
// Eligibility is based on the accepted A_eff-origin supernodal storage contract:
//   1. storage.valid == true
//   2. storage.source_of_truth_storage == true
//   3. storage.true_numeric_source == true
//      SLU-8R.5 contract: production transitional storage (bootstrapped from
//      CSC, §17.2(A)/(B) not fully replacing values) has true_numeric_source==false.
//      Using native solve on transitional storage risks false-pass or false-fail.
//      Production path falls back to CSC while true_numeric_source == false.
//      SLU-8R.5.5 (Numeric Source Switch) sets true_numeric_source == true when
//      the A_eff-origin factorization is accepted by the residual gate.
//      Gate 6 is tied to true_numeric_source, not merely to supernodal_solve_native.
//   4. !storage.supernodes.empty()
//
// NOTE (B2+ Option A): info_.within_panel_status is NOT a parameter and is NOT
// consulted here.  That status may describe the TRANSITIONAL CSC-bootstrapped
// §17.2(B) run; it must not gate native solve once storage.true_numeric_source
// has been accepted.  A_eff-origin pivot failure is already reflected in
// storage.true_numeric_source == false (condition 3 above).
// ---------------------------------------------------------------------------
template <class T, class Index>
bool can_use_supernodal_storage_solve(
    const supernodal_lu_storage<T, Index>& storage,
    std::string*                           reason_out)
{
    if (!storage.valid) {
        if (reason_out) *reason_out = "storage.valid == false";
        return false;
    }
    if (!storage.source_of_truth_storage) {
        if (reason_out) *reason_out = "storage.source_of_truth_storage == false";
        return false;
    }
    if (!storage.true_numeric_source) {
        if (reason_out) *reason_out =
            "true_numeric_source_false: transitional storage; "
            "R.5.5 Numeric Source Switch required for production native solve";
        return false;
    }
    if (storage.supernodes.empty()) {
        if (reason_out) *reason_out = "storage.supernodes is empty";
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// supernodal_l_solve_single_rhs
//
// Forward L solve for one RHS reading from supernodal_lu_storage.
// Processes supernodes in forward order (s = 0 .. nsup-1).
//
// For each supernode s:
//   Gather work_block[r] = work[row_indices[r]] for r in [0, R)
//   trsm('L','L','N','U', w, 1): solve unit-lower diagonal block
//   If R > w: gemv(-1, L_offdiag, x, y): update off-diagonal rows
//   Scatter work_block back to work
//
// L_offdiag = panel[w..R-1, 0..w-1] = off-diagonal L rows (rows w..R-1, lda=ld).
// Adapter usage: trsm + gemv for each supernode (even w==1 trivial cases).
// No direct tblas/tlapack calls.
// ---------------------------------------------------------------------------
template <class T, class Index>
void supernodal_l_solve_single_rhs(
    const supernodal_lu_storage<T, Index>& storage,
    std::vector<T>&                        work,
    supernodal_storage_solve_stats&        stats)
{
    const std::size_t nsup = storage.supernodes.size();
    const T one  = T(1);
    const T mone = T(-1);

    std::vector<T> work_block;

    for (std::size_t sj = 0; sj < nsup; ++sj) {
        const supernode_desc<Index>& desc = storage.supernodes[sj];
        const std::size_t w   = static_cast<std::size_t>(desc.num_cols);
        const std::size_t R   = desc.row_indices.size();
        const std::size_t ld  = static_cast<std::size_t>(desc.leading_dimension);
        const std::size_t off = static_cast<std::size_t>(desc.values_offset);

        if (w == 0 || R == 0) continue;

        const T* panel = storage.panel_values.data() + off;

        // Gather: R rows for this supernode panel
        work_block.resize(R);
        for (std::size_t r = 0; r < R; ++r) {
            work_block[r] = work[static_cast<std::size_t>(desc.row_indices[r])];
        }

        // L diagonal block trsm (unit-lower, w×w), lda=ld, ldb=R >= w
        sparse_lu_dense_kernel<T>::trsm('L', 'L', 'N', 'U',
            w, 1, one, panel, ld, work_block.data(), R);
        ++stats.trsm_count;
        ++stats.l_block_count;

        // Off-diagonal L update: rows w..R-1
        if (R > w) {
            const std::size_t m = R - w;
            // L_offdiag: m×w at panel+w (column-major, lda=ld)
            sparse_lu_dense_kernel<T>::gemv(m, w, mone, panel + w, ld,
                work_block.data(), one, work_block.data() + w);
            ++stats.gemv_count;
            ++stats.l_update_count;
        }

        // Scatter back
        for (std::size_t r = 0; r < R; ++r) {
            work[static_cast<std::size_t>(desc.row_indices[r])] = work_block[r];
        }
    }
}

// ---------------------------------------------------------------------------
// supernodal_u_solve_single_rhs
//
// Backward U solve for one RHS reading from supernodal_lu_storage.
// Processes supernodes in reverse order (s = nsup-1 .. 0).
//
// For each supernode s:
//   Gather work_block[c] = work[row_indices[c]] for c in [0, w)
//   trsm('L','U','N','N', w, 1): solve non-unit upper diagonal block
//   Apply U_segments: work[row_r] -= val * work_block[c]
//   Scatter work_block back to work
//
// U_segment apply updates work[row_r] for row_r < first_col of this supernode.
// These rows belong to earlier supernodes (not yet processed in backward pass).
//
// u_seg_col_ptr[c] = cumulative count for local cols 0..c-1 (relative to seg start).
// Absolute U_segment index for local col c:
//   [u_segment_start + u_seg_col_ptr[c], u_segment_start + u_seg_col_ptr[c+1])
//
// Adapter usage: trsm for diagonal block.
// U_segment apply is sparse (element-wise loop; no dense adapter call).
// No direct tblas/tlapack calls.
// ---------------------------------------------------------------------------
template <class T, class Index>
void supernodal_u_solve_single_rhs(
    const supernodal_lu_storage<T, Index>& storage,
    std::vector<T>&                        work,
    supernodal_storage_solve_stats&        stats)
{
    const std::size_t nsup = storage.supernodes.size();
    const T one = T(1);

    std::vector<T> work_block;

    for (std::size_t sj_idx = 0; sj_idx < nsup; ++sj_idx) {
        const std::size_t sj = nsup - 1u - sj_idx; // reverse order
        const supernode_desc<Index>& desc = storage.supernodes[sj];
        const std::size_t w   = static_cast<std::size_t>(desc.num_cols);
        const std::size_t ld  = static_cast<std::size_t>(desc.leading_dimension);
        const std::size_t off = static_cast<std::size_t>(desc.values_offset);

        if (w == 0) continue;

        const T* panel = storage.panel_values.data() + off;

        // Gather diagonal block rows (indices 0..w-1)
        work_block.resize(w);
        for (std::size_t c = 0; c < w; ++c) {
            work_block[c] = work[static_cast<std::size_t>(desc.row_indices[c])];
        }

        // U diagonal block backward solve (non-unit upper, w×w), lda=ld, ldb=w
        sparse_lu_dense_kernel<T>::trsm('L', 'U', 'N', 'N',
            w, 1, one, panel, ld, work_block.data(), w);
        ++stats.trsm_count;
        ++stats.u_block_count;

        // Apply U_segments: update rows < first_col of this supernode
        // u_seg_col_ptr.size() == w+1 when off-diagonal U exists
        if (!desc.u_seg_col_ptr.empty() &&
            desc.u_seg_col_ptr.size() == w + 1u) {
            const std::size_t seg_base =
                static_cast<std::size_t>(desc.u_segment_start);
            for (std::size_t c = 0; c < w; ++c) {
                const std::size_t u_begin =
                    seg_base + static_cast<std::size_t>(desc.u_seg_col_ptr[c]);
                const std::size_t u_end =
                    seg_base + static_cast<std::size_t>(desc.u_seg_col_ptr[c + 1u]);
                const T x_c = work_block[c];
                for (std::size_t k = u_begin; k < u_end; ++k) {
                    const std::size_t row_r = static_cast<std::size_t>(
                        storage.U_segments.row_ind[k]);
                    work[row_r] -= storage.U_segments.values[k] * x_c;
                    ++stats.u_seg_apply_count;
                }
            }
        }

        // Scatter back
        for (std::size_t c = 0; c < w; ++c) {
            work[static_cast<std::size_t>(desc.row_indices[c])] = work_block[c];
        }
    }
}

// ---------------------------------------------------------------------------
// supernodal_storage_solve_single_rhs
//
// §18.2 full solve pipeline for one RHS.
// Reads directly from supernodal_lu_storage (not CSC baseline L/U).
//
// Permutation convention matches solve_baseline_storage exactly:
//   row_perm[new_i] = old_i, col_perm[new_j] = old_j.
//
// Does NOT call can_use_supernodal_storage_solve (caller must check).
// Does NOT validate b.size() == n (caller must check).
// Updates stats with adapter call counts and timing.
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
supernodal_storage_solve_single_rhs(
    const supernodal_lu_storage<T, Index>& storage,
    Index                                  n,
    const std::vector<T>&                  b,
    supernodal_storage_solve_stats&        stats)
{
    const std::size_t un = static_cast<std::size_t>(n);
    stats.nrhs      = 1;
    stats.attempted = true;

    auto t0 = std::chrono::steady_clock::now();

    // Step 1: Dr row scaling
    std::vector<T> work(un);
    if (storage.Dr.empty()) {
        work = b;
    } else {
        for (std::size_t i = 0; i < un; ++i) {
            work[i] = storage.Dr[i] * b[i];
        }
    }

    // Step 2: Row permutation P: rhs2[new_i] = work[row_perm[new_i]]
    if (!storage.row_perm.empty()) {
        std::vector<T> rhs2(un);
        for (std::size_t ni = 0; ni < un; ++ni) {
            rhs2[ni] = work[static_cast<std::size_t>(storage.row_perm[ni])];
        }
        work = rhs2;
    }

    // Step 3: Supernodal L forward solve
    supernodal_l_solve_single_rhs(storage, work, stats);

    // Step 4: Supernodal U backward solve
    supernodal_u_solve_single_rhs(storage, work, stats);

    // Step 5: Column permutation Q exit: x_orig[col_perm[new_j]] = work[new_j]
    std::vector<T> x_orig(un);
    if (storage.col_perm.empty()) {
        x_orig = work;
    } else {
        for (std::size_t nj = 0; nj < un; ++nj) {
            x_orig[static_cast<std::size_t>(storage.col_perm[nj])] = work[nj];
        }
    }

    // Step 6: Dc column scaling
    std::vector<T> x(un);
    if (storage.Dc.empty()) {
        x = x_orig;
    } else {
        for (std::size_t i = 0; i < un; ++i) {
            x[i] = storage.Dc[i] * x_orig[i];
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    stats.ticks = static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    stats.succeeded = true;
    return x;
}

} // namespace sparse_lu_detail

// ===========================================================================
// sparse_lu_make_supernodal_factor_for_testing
//
// SLU-8R.5 testing hook: constructs a factor from manually provided
// supernodal_lu_storage.  NOT a numeric factorization; use only in tests.
//
// Allows injecting synthetic L/U values in panel_values with explicit
// permutations and scaling to test the native §18.2 solve path.
//
// The returned factor:
//   - has_supernodal_storage() == true
//   - supernodal_solve_native() == true
//   - solve() dispatches to storage-native path
//   - baseline_ L/U CSC is empty (CSC fallback is not available)
//
// If native solve cannot be used, solve() throws (no CSC fallback).
//
// storage.valid must be true, source_of_truth_storage must be true,
// supernodes must be non-empty.
// ===========================================================================
template <class T, class Index>
sparse_lu_factorization<T, Index>
sparse_lu_make_supernodal_factor_for_testing(
    Index n,
    const supernodal_lu_storage<T, Index>& storage)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_make_supernodal_factor_for_testing: Index must be signed");

    if (!storage.valid) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_make_supernodal_factor_for_testing: storage.valid == false");
    }
    if (!storage.source_of_truth_storage) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_make_supernodal_factor_for_testing: "
            "storage.source_of_truth_storage == false");
    }
    if (storage.supernodes.empty()) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_make_supernodal_factor_for_testing: "
            "storage.supernodes is empty");
    }

    sparse_lu_factorization<T, Index> fac;
    fac.set_supernodal_factor_for_testing_(n, storage);
    return fac;
}

// ===========================================================================
// sparse_lu_supernodal_storage_solve
//
// SLU-8R.5 free function: explicitly invokes §18.2 storage-native solve.
// Bypasses the solve() dispatch; directly calls the native pipeline.
//
// Useful in tests to:
//   - verify native path is used (not CSC fallback)
//   - inspect per-call stats (l_block_count, trsm_count, etc.)
//   - verify multi-RHS by calling this function nrhs times
//
// Throws vcp::state_error if:
//   - fac.has_supernodal_storage() == false
//   - can_use_supernodal_storage_solve() returns false
//   - b.size() != fac.info().n
//
// Gate 6 PENDING: correctness depends on panel_values being valid L/U.
// For diagonal matrices (empty U_segments), result matches baseline CSC.
// ===========================================================================
template <class T, class Index>
std::vector<T>
sparse_lu_supernodal_storage_solve(
    const sparse_lu_factorization<T, Index>&          fac,
    const std::vector<T>&                             b,
    sparse_lu_detail::supernodal_storage_solve_stats* stats_out)
{
    sparse_lu_detail::supernodal_storage_solve_stats stats;
    stats.attempted = true;

    if (!fac.has_supernodal_storage()) {
        stats.fallback_csc    = true;
        stats.fallback_reason = "has_supernodal_storage == false";
        if (stats_out) *stats_out = stats;
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_supernodal_storage_solve: "
            "factor does not have supernodal storage");
    }

    const Index n = fac.info().n;
    if (static_cast<Index>(b.size()) != n) {
        stats.fallback_csc    = true;
        stats.fallback_reason = "b.size() != n";
        if (stats_out) *stats_out = stats;
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_supernodal_storage_solve: b.size() != n");
    }

    std::string reason;
    if (!sparse_lu_detail::can_use_supernodal_storage_solve(
            fac.supernodal_storage(), &reason)) {
        stats.fallback_csc    = true;
        stats.fallback_reason = reason;
        if (stats_out) *stats_out = stats;
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_supernodal_storage_solve: cannot use native solve: ", reason);
    }

    std::vector<T> x = sparse_lu_detail::supernodal_storage_solve_single_rhs(
        fac.supernodal_storage(), n, b, stats);

    if (stats_out) *stats_out = stats;
    return x;
}

// ===========================================================================
// sparse_lu_supernodal_storage_solve_all_rhs
//
// SLU-8R.5.1 multi-RHS free function: invokes §18.2 storage-native solve
// for every column in B (vector of column vectors).
//
// Design:
//   - Each column b_j is solved with supernodal_storage_solve_single_rhs.
//   - stats.nrhs == B.size() after completion.
//   - stats counts (l_block_count, trsm_count, ...) accumulate across all RHS.
//   - Throws vcp::state_error if any precondition fails (same as single-RHS).
//
// Returns vector of solution columns X: X[j] = A^{-1} * B[j].
//
// This is the true multi-RHS entry point for R.5 testing.
// Sequential single-RHS calls to sparse_lu_supernodal_storage_solve are NOT
// a substitute for this function as multi-RHS evidence.
//
// Gate 6 PENDING: correctness depends on panel_values being valid L/U.
// ===========================================================================
template <class T, class Index>
std::vector<std::vector<T>>
sparse_lu_supernodal_storage_solve_all_rhs(
    const sparse_lu_factorization<T, Index>&                fac,
    const std::vector<std::vector<T>>&                      B,
    sparse_lu_detail::supernodal_storage_solve_stats*       stats_out)
{
    sparse_lu_detail::supernodal_storage_solve_stats accum;
    accum.attempted = true;

    if (!fac.has_supernodal_storage()) {
        accum.fallback_csc    = true;
        accum.fallback_reason = "has_supernodal_storage == false";
        if (stats_out) *stats_out = accum;
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_supernodal_storage_solve_all_rhs: "
            "factor does not have supernodal storage");
    }

    const Index n = fac.info().n;

    // Check all RHS columns have correct size before entering solve loop
    for (std::size_t j = 0; j < B.size(); ++j) {
        if (static_cast<Index>(B[j].size()) != n) {
            accum.fallback_csc    = true;
            accum.fallback_reason = "B[j].size() != n";
            if (stats_out) *stats_out = accum;
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_supernodal_storage_solve_all_rhs: B[j].size() != n");
        }
    }

    std::string reason;
    if (!sparse_lu_detail::can_use_supernodal_storage_solve(
            fac.supernodal_storage(), &reason)) {
        accum.fallback_csc    = true;
        accum.fallback_reason = reason;
        if (stats_out) *stats_out = accum;
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_supernodal_storage_solve_all_rhs: cannot use native solve: ",
            reason);
    }

    std::vector<std::vector<T>> X;
    X.reserve(B.size());

    for (std::size_t j = 0; j < B.size(); ++j) {
        sparse_lu_detail::supernodal_storage_solve_stats col_stats;
        std::vector<T> xj = sparse_lu_detail::supernodal_storage_solve_single_rhs(
            fac.supernodal_storage(), n, B[j], col_stats);
        X.push_back(xj);

        // Accumulate per-column stats
        accum.l_block_count      += col_stats.l_block_count;
        accum.l_update_count     += col_stats.l_update_count;
        accum.u_block_count      += col_stats.u_block_count;
        accum.u_seg_apply_count  += col_stats.u_seg_apply_count;
        accum.trsm_count         += col_stats.trsm_count;
        accum.gemv_count         += col_stats.gemv_count;
        accum.gemm_count         += col_stats.gemm_count;
        accum.ticks              += col_stats.ticks;
    }

    accum.nrhs      = B.size();
    accum.succeeded = true;
    if (stats_out) *stats_out = accum;
    return X;
}

#endif // VCP_TSPARSE_SPARSE_LU_SUPERNODAL_SOLVE_IMPL_HPP
