// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-9 supernode-aware triangular solve -- internal implementation.
//
// SLU-9 supernode-aware solve helper.
// This helper uses supernode metadata only to group column traversal and
// gather/scatter work vectors. Numeric updates are still performed from the
// actual emitted CSC L/U factors. Symbolic panel rows are not treated as the
// final numeric row structure under threshold partial pivoting.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - all storage type definitions (csc_storage, baseline_lu_storage, etc.)
//   - sparse_lu_scalar_policy
//   - validate_baseline_lu_storage_for_solve (from SLU-2 detail)
//   - sparse_lu_is_valid_supernode_symbolic (from SLU-7 section)
//   - sparse_lu_supernode_symbolic type
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse.hpp>  (umbrella)

#ifndef VCP_TSPARSE_SPARSE_LU_SUPERNODE_SOLVE_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SUPERNODE_SOLVE_IMPL_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// supernode_aware_forward_solve_L
//
// SLU-9: L forward solve grouped by supernode column ranges.
// Uses supernode metadata to group column traversal only.
// Numeric updates use actual CSC L.row_ind / L.values (not symbolic panel rows).
//
// For supernode s with columns [supernode_ptr[s], supernode_ptr[s+1]):
//   for each column j in increasing order:
//     for each stored entry (i, v) in L column j:
//       if i > j: y[i] -= v * y[j]   (strict lower, unit diagonal contract)
//       if i == j: skip              (unit diagonal: ignore stored value)
//       if i < j: error              (upper entry in L is invalid)
//
// Precondition: y holds the permuted/scaled RHS.
// Postcondition: y = L^{-1} * rhs (same result as csc_forward_solve_L).
// ---------------------------------------------------------------------------
template <class T, class Index>
void supernode_aware_forward_solve_L(
    const csc_storage<T, Index>& L,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    std::vector<T>& y)
{
    (void)n; // n used only for bounds in outer caller
    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index col_begin = supernodes.supernode_ptr[s];
        const Index col_end   = supernodes.supernode_ptr[s + 1u];

        for (Index j = col_begin; j < col_end; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            for (Index k = L.col_ptr[sj]; k < L.col_ptr[sj + 1u]; ++k) {
                const std::size_t sk = static_cast<std::size_t>(k);
                const Index i = L.row_ind[sk];
                if (i < j) {
                    vcp::throw_error<vcp::invalid_argument>(
                        "sparse_lu supernode-aware L solve: "
                        "upper entry in L (invalid factor storage)");
                }
                if (i == j) {
                    continue; // unit diagonal contract: skip stored diagonal
                }
                // i > j: strict lower entry -- forward substitution scatter
                y[static_cast<std::size_t>(i)] -= L.values[sk] * y[sj];
            }
        }
    }
}

// ---------------------------------------------------------------------------
// supernode_aware_backward_solve_U
//
// SLU-9: U backward solve grouped by supernode column ranges.
// Uses supernode metadata to group column traversal only.
// Numeric updates use actual CSC U.row_ind / U.values (not symbolic panel rows).
//
// For supernode s (decreasing) with columns [supernode_ptr[s], supernode_ptr[s+1]):
//   for each column j in decreasing order:
//     locate diagonal U[j,j]; error if missing or zero
//     x[j] /= U[j,j]
//     for each stored entry (i, v) in U column j with i < j:
//       x[i] -= v * x[j]   (upper triangle scatter)
//
// Precondition: x holds the result after L forward solve.
// Postcondition: x = U^{-1} * L^{-1} * rhs (same as csc_backward_solve_U).
// ---------------------------------------------------------------------------
template <class T, class Index>
void supernode_aware_backward_solve_U(
    const csc_storage<T, Index>& U,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    std::vector<T>& x)
{
    (void)n;
    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    for (std::size_t si = nsup; si-- > 0u; ) {
        const Index col_begin = supernodes.supernode_ptr[si];
        const Index col_end   = supernodes.supernode_ptr[si + 1u];

        // Within supernode: traverse columns in decreasing order
        for (Index jj = col_end; jj-- > col_begin; ) {
            const Index j = jj;
            const std::size_t sj = static_cast<std::size_t>(j);

            // First pass: locate diagonal, validate no lower entries
            bool found_diag = false;
            T    diag_val   = T(0);

            for (Index k = U.col_ptr[sj]; k < U.col_ptr[sj + 1u]; ++k) {
                const std::size_t sk = static_cast<std::size_t>(k);
                const Index i = U.row_ind[sk];
                if (i > j) {
                    vcp::throw_error<vcp::invalid_argument>(
                        "sparse_lu supernode-aware U solve: "
                        "lower entry in U (invalid factor storage)");
                }
                if (i == j) {
                    diag_val   = U.values[sk];
                    found_diag = true;
                }
            }

            if (!found_diag) {
                vcp::throw_error<vcp::invalid_argument>(
                    "sparse_lu supernode-aware U solve: "
                    "missing diagonal in U (zero_pivot)");
            }
            if (sparse_lu_scalar_policy<T>::is_exact_zero(diag_val)) {
                vcp::throw_error<vcp::state_error>(
                    "sparse_lu supernode-aware U solve: "
                    "zero diagonal in U (numerical_singularity)");
            }

            x[sj] /= diag_val;

            // Second pass: scatter strict upper-triangle contributions
            for (Index k = U.col_ptr[sj]; k < U.col_ptr[sj + 1u]; ++k) {
                const std::size_t sk = static_cast<std::size_t>(k);
                const Index i = U.row_ind[sk];
                if (i < j) {
                    x[static_cast<std::size_t>(i)] -= U.values[sk] * x[sj];
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// make_singleton_supernodes
//
// Builds trivial singleton supernode metadata: each column is its own
// supernode.  Used as fallback when stored supernode_info is invalid.
// Singleton supernodes preserve the sequential column traversal order
// and satisfy all sparse_lu_is_valid_supernode_symbolic invariants.
// ---------------------------------------------------------------------------
template <class Index>
sparse_lu_supernode_symbolic<Index>
make_singleton_supernodes(Index n)
{
    static_assert(std::is_signed<Index>::value,
                  "make_singleton_supernodes: Index must be signed");
    sparse_lu_supernode_symbolic<Index> sn;
    if (n < Index(0)) { sn.valid = false; return sn; }

    const std::size_t un = static_cast<std::size_t>(n);
    sn.supernode_ptr.resize(un + 1u);
    sn.column_to_supernode.resize(un);
    sn.parent.resize(un, Index(-1));
    sn.row_ptr.resize(un + 1u, Index(0));
    // row_ind empty: no structural rows stored (not needed for solve)
    for (Index i = Index(0); i <= n; ++i) {
        sn.supernode_ptr[static_cast<std::size_t>(i)] = i;
    }
    for (Index i = Index(0); i < n; ++i) {
        sn.column_to_supernode[static_cast<std::size_t>(i)] = i;
    }
    sn.valid = true;
    return sn;
}

// ---------------------------------------------------------------------------
// solve_baseline_storage_supernode_aware
//
// SLU-9: full solve pipeline for baseline_csc storage with supernode-grouped
// L/U traversal.  Mathematically equivalent to solve_baseline_storage.
//
// Pipeline:
//   1. Dr scaling        (row equilibration, identity if Dr empty)
//   2. Row permutation P (rhs2[new_i] = rhs1[row_perm[new_i]])
//   3. L^{-1}            (supernode-aware CSC forward solve)
//   4. U^{-1}            (supernode-aware CSC backward solve)
//   5. Column perm Q     (x_orig[col_perm[new_j]] = z2[new_j])
//   6. Dc scaling        (column equilibration, identity if Dc empty)
//
// If supernode_info.valid == false or fails is_valid_supernode_symbolic,
// falls back to singleton supernodes (sequential column traversal).
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
solve_baseline_storage_supernode_aware(
    const baseline_lu_storage<T, Index>& lu,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernode_info,
    const std::vector<T>& b)
{
    if (static_cast<Index>(b.size()) != n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu supernode-aware solve: b.size() != n");
    }

    if (n == Index(0)) {
        return std::vector<T>();
    }

    // Structural validation via SLU-2 validator
    validate_baseline_lu_storage_for_solve(lu, n);

    const std::size_t un = static_cast<std::size_t>(n);

    // Step 1: Dr scaling (row equilibration)
    std::vector<T> work(un);
    if (lu.Dr.empty()) {
        work = b;
    } else {
        for (std::size_t i = 0u; i < un; ++i) {
            work[i] = lu.Dr[i] * b[i];
        }
    }

    // Step 2: Row permutation P  (work[new_i] = scaled[row_perm[new_i]])
    if (!lu.row_perm.empty()) {
        std::vector<T> tmp(un);
        for (Index new_i = Index(0); new_i < n; ++new_i) {
            const std::size_t sni = static_cast<std::size_t>(new_i);
            tmp[sni] = work[static_cast<std::size_t>(lu.row_perm[sni])];
        }
        work = tmp;
    }

    // Select supernode partition: use stored metadata if valid, else singleton.
    // Singleton ensures correctness when metadata is unavailable.
    sparse_lu_supernode_symbolic<Index> singleton;
    const sparse_lu_supernode_symbolic<Index>* sn =
        &supernode_info;

    if (!supernode_info.valid ||
        !sparse_lu_is_valid_supernode_symbolic(n, supernode_info)) {
        singleton = make_singleton_supernodes(n);
        sn = &singleton;
    }

    // Step 3: L^{-1} (supernode-aware forward solve)
    // Numeric updates use actual CSC L.row_ind / L.values.
    // Symbolic panel rows (panel_row_ind) are never accessed here.
    supernode_aware_forward_solve_L(lu.L, n, *sn, work);

    // Step 4: U^{-1} (supernode-aware backward solve)
    // Numeric updates use actual CSC U.row_ind / U.values.
    supernode_aware_backward_solve_U(lu.U, n, *sn, work);

    // Step 5: Column permutation Q exit
    std::vector<T> x_orig(un);
    if (lu.col_perm.empty()) {
        x_orig = work;
    } else {
        for (Index new_j = Index(0); new_j < n; ++new_j) {
            const std::size_t snj = static_cast<std::size_t>(new_j);
            x_orig[static_cast<std::size_t>(lu.col_perm[snj])] = work[snj];
        }
    }

    // Step 6: Dc scaling (column equilibration)
    if (lu.Dc.empty()) {
        return x_orig;
    }
    std::vector<T> x(un);
    for (std::size_t i = 0u; i < un; ++i) {
        x[i] = lu.Dc[i] * x_orig[i];
    }
    return x;
}

// ---------------------------------------------------------------------------
// sparse_lu_solve_upper_diag_block_column_major
//
// SLU-14: Manual back-substitution for a small dense upper triangular block
// stored in column-major layout.
//
// Layout: diag_block_values[block_offset + lc*width + lr] == U_block(lr, lc)
//   lc = local column index (0..width-1)
//   lr = local row index    (0..width-1)
//
// Solves U_block * y = y in-place by back substitution:
//   for i = width-1 downto 0:
//     y[i] -= sum_{j=i+1..width-1} U_block(i,j) * y[j]
//     require U_block(i,i) != 0  (throws state_error if zero)
//     y[i] /= U_block(i,i)
//
// w == 0: no-op. w == 1: scalar divide.
// Zero diagonal throws state_error (numerical_singularity).
// Does not call BLAS/LAPACK/tblas/tlapack.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_lu_solve_upper_diag_block_column_major(
    Index width,
    const std::vector<T>& block_values,
    Index block_offset,
    std::vector<T>& y)
{
    // SLU-14.1: validate inputs before any vector access
    if (width < Index(0) || block_offset < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_solve_upper_diag_block_column_major: negative width or offset");
    }

    const std::size_t sw = static_cast<std::size_t>(width);
    const std::size_t so = static_cast<std::size_t>(block_offset);

    if (sw > y.size()) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_solve_upper_diag_block_column_major: y.size() < width");
    }

    if (so > block_values.size()) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_solve_upper_diag_block_column_major: block_offset out of range");
    }

    // overflow-safe: need sw*sw <= block_values.size() - so
    if (sw != 0u && sw > (block_values.size() - so) / sw) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_solve_upper_diag_block_column_major: block is out of range");
    }

    if (sw == 0u) return;

    for (Index i = width - Index(1); i >= Index(0); --i) {
        const std::size_t si = static_cast<std::size_t>(i);
        // Subtract off-diagonal contributions: U_block(i, j) for j > i
        for (Index j = i + Index(1); j < width; ++j) {
            const std::size_t idx =
                static_cast<std::size_t>(block_offset + j * width + i);
            y[si] -= block_values[idx] * y[static_cast<std::size_t>(j)];
        }
        // Divide by diagonal U_block(i, i)
        const std::size_t diag_idx =
            static_cast<std::size_t>(block_offset + i * width + i);
        const T diag_val = block_values[diag_idx];
        if (sparse_lu_scalar_policy<T>::is_exact_zero(diag_val)) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_solve_upper_diag_block_column_major: "
                "zero diagonal pivot (numerical_singularity)");
        }
        y[si] /= diag_val;

        if (i == Index(0)) break;
    }
}

// ---------------------------------------------------------------------------
// supernode_aware_backward_solve_U_diag_block
//
// SLU-14: U backward solve using validated diag_block_values for within-supernode
// diagonal block, and actual CSC U data for off-diagonal contributions.
//
// For each supernode s in reverse order:
//   1. Gather x[col_begin..col_end) into local vector.
//   2. Solve U_diag_block * local = local using diag_block_values (back-sub).
//   3. Write local back to x[col_begin..col_end).
//   4. Scatter off-diagonal: for each col j in [col_begin,col_end), for each
//      U CSC entry (i,j) with i < col_begin: x[i] -= U[i,j] * x[j].
//
// Preconditions (caller ensures):
//   sn_numeric is valid, invariant-checked, and CSC-consistent.
//   diag_block_ptr sizes match supernodes.supernode_ptr.
//
// Result is numerically equivalent to supernode_aware_backward_solve_U
// when diag_block_values is consistent with the actual U CSC factors.
// ---------------------------------------------------------------------------
template <class T, class Index>
void supernode_aware_backward_solve_U_diag_block(
    const csc_storage<T, Index>& U,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& sn_numeric,
    std::vector<T>& x)
{
    (void)n;
    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    for (std::size_t si = nsup; si-- > 0u; ) {
        const Index col_begin = supernodes.supernode_ptr[si];
        const Index col_end   = supernodes.supernode_ptr[si + 1u];
        const Index width     = col_end - col_begin;
        const std::size_t uwidth = static_cast<std::size_t>(width);

        // Step 1: gather local RHS from x
        std::vector<T> local(uwidth);
        for (std::size_t lc = 0u; lc < uwidth; ++lc) {
            local[lc] = x[static_cast<std::size_t>(col_begin) + lc];
        }

        // Step 2: solve diagonal block using diag_block_values (back-substitution)
        sparse_lu_solve_upper_diag_block_column_major(
            width, sn_numeric.diag_block_values, sn_numeric.diag_block_ptr[si], local);

        // Step 3: write solved values back
        for (std::size_t lc = 0u; lc < uwidth; ++lc) {
            x[static_cast<std::size_t>(col_begin) + lc] = local[lc];
        }

        // Step 4: scatter off-diagonal U contributions (rows < col_begin)
        for (Index j = col_begin; j < col_end; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const T xj = x[sj];
            for (Index k = U.col_ptr[sj]; k < U.col_ptr[sj + 1u]; ++k) {
                const std::size_t sk = static_cast<std::size_t>(k);
                const Index i = U.row_ind[sk];
                if (i < col_begin) {
                    x[static_cast<std::size_t>(i)] -= U.values[sk] * xj;
                }
                // i >= col_begin: diagonal block already handled above
            }
        }
    }
}

// ---------------------------------------------------------------------------
// solve_baseline_storage_supernode_aware_diag_block
//
// SLU-14: Full solve pipeline using diag_block_values for U diagonal block solve.
//
// Pipeline (identical to solve_baseline_storage_supernode_aware except step 4):
//   1. Dr scaling
//   2. Row permutation P
//   3. L^{-1} (supernode-aware CSC forward solve, unchanged from SLU-9)
//   4. U^{-1} with diag-block: diagonal block via diag_block_values,
//      off-diagonal contributions via actual CSC U
//   5. Column permutation Q
//   6. Dc scaling
//
// Preconditions (caller must ensure):
//   sn_numeric.valid == true
//   sparse_lu_is_valid_supernode_numeric passes
//   sparse_lu_verify_supernode_numeric_against_csc passes
//
// Result is numerically equivalent to solve_baseline_storage_supernode_aware
// when diag_block_values is consistent with CSC U factors.
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
solve_baseline_storage_supernode_aware_diag_block(
    const baseline_lu_storage<T, Index>& lu,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernode_info,
    const sparse_lu_supernode_numeric<T, Index>& sn_numeric,
    const std::vector<T>& b)
{
    if (static_cast<Index>(b.size()) != n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu supernode-aware diag-block solve: b.size() != n");
    }

    if (n == Index(0)) {
        return std::vector<T>();
    }

    validate_baseline_lu_storage_for_solve(lu, n);

    const std::size_t un = static_cast<std::size_t>(n);

    // Step 1: Dr scaling
    std::vector<T> work(un);
    if (lu.Dr.empty()) {
        work = b;
    } else {
        for (std::size_t i = 0u; i < un; ++i) {
            work[i] = lu.Dr[i] * b[i];
        }
    }

    // Step 2: Row permutation P
    if (!lu.row_perm.empty()) {
        std::vector<T> tmp(un);
        for (Index new_i = Index(0); new_i < n; ++new_i) {
            const std::size_t sni = static_cast<std::size_t>(new_i);
            tmp[sni] = work[static_cast<std::size_t>(lu.row_perm[sni])];
        }
        work = tmp;
    }

    // Step 3: L^{-1} (supernode-aware CSC forward solve, same as SLU-9)
    supernode_aware_forward_solve_L(lu.L, n, supernode_info, work);

    // Step 4: U^{-1} with diag-block (SLU-14: uses diag_block_values for block solve)
    supernode_aware_backward_solve_U_diag_block(lu.U, n, supernode_info, sn_numeric, work);

    // Step 5: Column permutation Q
    std::vector<T> x_orig(un);
    if (lu.col_perm.empty()) {
        x_orig = work;
    } else {
        for (Index new_j = Index(0); new_j < n; ++new_j) {
            const std::size_t snj = static_cast<std::size_t>(new_j);
            x_orig[static_cast<std::size_t>(lu.col_perm[snj])] = work[snj];
        }
    }

    // Step 6: Dc scaling
    if (lu.Dc.empty()) {
        return x_orig;
    }
    std::vector<T> x(un);
    for (std::size_t i = 0u; i < un; ++i) {
        x[i] = lu.Dc[i] * x_orig[i];
    }
    return x;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_SUPERNODE_SOLVE_IMPL_HPP
