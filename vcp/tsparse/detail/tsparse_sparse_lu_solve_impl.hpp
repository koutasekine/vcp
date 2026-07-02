// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-2 baseline triangular solve -- internal implementation.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER all storage
// type definitions, sparse_lu_scalar_policy, sparse_lu_is_valid_csc_storage,
// and sparse_lu_inverse_permutation are in scope.
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse.hpp>  (umbrella)

#ifndef VCP_TSPARSE_SPARSE_LU_SOLVE_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SOLVE_IMPL_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// validate_optional_inverse_permutation  (SLU-2.1)
// Validates inv_row_perm or inv_col_perm against its forward permutation.
//
// Semantics:
//   - inverse empty  → accepted as "unspecified" (identity implied)
//   - inverse non-empty:
//       1. size must equal n
//       2. inverse must itself be a valid permutation (range + no duplicates)
//       3. inverse must equal the exact inverse of effective_forward
//          where effective_forward = forward if non-empty, else identity(n)
// ---------------------------------------------------------------------------
template <class Index>
void validate_optional_inverse_permutation(
    const std::vector<Index>& forward,
    const std::vector<Index>& inverse,
    Index n,
    const char* inverse_name)
{
    if (inverse.empty()) return;

    if (static_cast<Index>(inverse.size()) != n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: ", inverse_name, " size mismatch");
    }

    // Validate that inverse is itself a valid permutation (throws on range/duplicate)
    (void)sparse_lu_inverse_permutation(inverse);

    // Determine effective forward permutation
    std::vector<Index> effective_forward;
    if (forward.empty()) {
        effective_forward = sparse_lu_identity_permutation(n);
    } else {
        effective_forward = forward;
    }

    // Compute the expected inverse from the effective forward
    std::vector<Index> expected_inverse = sparse_lu_inverse_permutation(effective_forward);

    if (inverse != expected_inverse) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: ", inverse_name,
            " is inconsistent with forward permutation");
    }
}

// ---------------------------------------------------------------------------
// validate_baseline_lu_storage_for_solve
// Checks structural validity of baseline_lu_storage.
// Does NOT check U diagonal (that is detected during backward solve).
// ---------------------------------------------------------------------------
template <class T, class Index>
void validate_baseline_lu_storage_for_solve(
    const baseline_lu_storage<T, Index>& lu,
    Index n)
{
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: negative n");
    }
    if (!sparse_lu_is_valid_csc_storage(lu.L, n, n)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: L is not a valid n x n CSC storage");
    }
    if (!sparse_lu_is_valid_csc_storage(lu.U, n, n)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: U is not a valid n x n CSC storage");
    }
    if (!lu.row_perm.empty()) {
        if (static_cast<Index>(lu.row_perm.size()) != n) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu validate: row_perm size mismatch");
        }
        // Validate by attempting to compute the inverse (throws on range/duplicate)
        const std::vector<Index> tmp = sparse_lu_inverse_permutation(lu.row_perm);
        (void)tmp;
    }
    // SLU-2.1: validate inv_row_perm consistency (empty = unspecified, accepted)
    validate_optional_inverse_permutation(lu.row_perm, lu.inv_row_perm, n, "inv_row_perm");
    if (!lu.col_perm.empty()) {
        if (static_cast<Index>(lu.col_perm.size()) != n) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu validate: col_perm size mismatch");
        }
        const std::vector<Index> tmp = sparse_lu_inverse_permutation(lu.col_perm);
        (void)tmp;
    }
    // SLU-2.1: validate inv_col_perm consistency (empty = unspecified, accepted)
    validate_optional_inverse_permutation(lu.col_perm, lu.inv_col_perm, n, "inv_col_perm");
    if (!lu.Dr.empty() && static_cast<Index>(lu.Dr.size()) != n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: Dr size mismatch");
    }
    if (!lu.Dc.empty() && static_cast<Index>(lu.Dc.size()) != n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu validate: Dc size mismatch");
    }
}

// ---------------------------------------------------------------------------
// csc_forward_solve_L
// Solve L * y = rhs, L is unit lower triangular CSC.
// Unit diagonal contract: stored diagonal entries are silently ignored.
// Throws if any upper entry (row < col) is found in L.
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
csc_forward_solve_L(
    const csc_storage<T, Index>& L,
    Index n,
    const std::vector<T>& rhs)
{
    std::vector<T> y(rhs);
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index k = L.col_ptr[sj]; k < L.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index i = L.row_ind[sk];
            if (i < j) {
                vcp::throw_error<vcp::invalid_argument>(
                    "sparse_lu L solve: upper entry in L (invalid factor storage)");
            }
            if (i == j) {
                continue; // unit diagonal contract: skip stored diagonal value
            }
            // i > j: strict lower entry, forward substitution update
            y[static_cast<std::size_t>(i)] -= L.values[sk] * y[sj];
        }
    }
    return y;
}

// ---------------------------------------------------------------------------
// csc_backward_solve_U
// Solve U * x = y, U is upper triangular CSC with explicit diagonal.
// Throws if any lower entry (row > col) is found in U.
// Throws if the diagonal entry of any column is missing or zero.
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
csc_backward_solve_U(
    const csc_storage<T, Index>& U,
    Index n,
    const std::vector<T>& y)
{
    std::vector<T> x(y);
    for (Index jj = n; jj-- > 0; ) {
        const Index j = jj;
        const std::size_t sj = static_cast<std::size_t>(j);

        // First pass: locate diagonal and reject any lower-triangle entries
        bool found_diag = false;
        T diag_val = T(0);
        for (Index k = U.col_ptr[sj]; k < U.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index i = U.row_ind[sk];
            if (i > j) {
                vcp::throw_error<vcp::invalid_argument>(
                    "sparse_lu U solve: lower entry in U (invalid factor storage)");
            }
            if (i == j) {
                diag_val   = U.values[sk];
                found_diag = true;
            }
        }

        if (!found_diag) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu U solve: missing diagonal in U (zero_pivot)");
        }
        if (sparse_lu_scalar_policy<T>::is_exact_zero(diag_val)) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu U solve: zero diagonal in U (numerical_singularity)");
        }

        x[sj] /= diag_val;

        // Second pass: scatter the column's strict upper-triangle contributions
        for (Index k = U.col_ptr[sj]; k < U.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index i = U.row_ind[sk];
            if (i < j) {
                x[static_cast<std::size_t>(i)] -= U.values[sk] * x[sj];
            }
        }
    }
    return x;
}

// ---------------------------------------------------------------------------
// solve_baseline_storage
// Full solve pipeline for baseline_csc storage:
//   rhs1       = Dr * b                          (Dr scaling, identity if Dr empty)
//   rhs2[new]  = rhs1[row_perm[new]]             (row permutation, identity if empty)
//   z1         = L^{-1} * rhs2                   (CSC forward solve, unit diagonal L)
//   z2         = U^{-1} * z1                     (CSC backward solve, explicit diagonal U)
//   x_orig[col_perm[new]] = z2[new]              (col permutation exit, identity if empty)
//   x          = Dc * x_orig                     (Dc scaling, identity if Dc empty)
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
solve_baseline_storage(
    const baseline_lu_storage<T, Index>& lu,
    Index n,
    const std::vector<T>& b)
{
    if (static_cast<Index>(b.size()) != n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu solve: b.size() != n");
    }

    if (n == Index(0)) {
        return std::vector<T>();
    }

    // Validate CSC structure, permutation sizes, scaling sizes
    validate_baseline_lu_storage_for_solve(lu, n);

    const std::size_t un = static_cast<std::size_t>(n);

    // Step 1: Dr scaling (row equilibration entry)
    std::vector<T> rhs1(un);
    if (lu.Dr.empty()) {
        rhs1 = b;
    } else {
        for (std::size_t i = 0; i < un; ++i) {
            rhs1[i] = lu.Dr[i] * b[i];
        }
    }

    // Step 2: Row permutation P  (rhs2[new_i] = rhs1[row_perm[new_i]])
    std::vector<T> rhs2(un);
    if (lu.row_perm.empty()) {
        rhs2 = rhs1;
    } else {
        for (Index new_i = Index(0); new_i < n; ++new_i) {
            const std::size_t sni = static_cast<std::size_t>(new_i);
            rhs2[sni] = rhs1[static_cast<std::size_t>(lu.row_perm[sni])];
        }
    }

    // Step 3: L^{-1} rhs2  (CSC forward solve, unit diagonal contract)
    std::vector<T> z1 = csc_forward_solve_L(lu.L, n, rhs2);

    // Step 4: U^{-1} z1  (CSC backward solve, explicit diagonal)
    std::vector<T> z2 = csc_backward_solve_U(lu.U, n, z1);

    // Step 5: Column permutation Q exit  (x_orig[col_perm[new_j]] = z2[new_j])
    std::vector<T> x_orig(un);
    if (lu.col_perm.empty()) {
        x_orig = z2;
    } else {
        for (Index new_j = Index(0); new_j < n; ++new_j) {
            const std::size_t snj = static_cast<std::size_t>(new_j);
            x_orig[static_cast<std::size_t>(lu.col_perm[snj])] = z2[snj];
        }
    }

    // Step 6: Dc scaling (column equilibration exit)
    if (lu.Dc.empty()) {
        return x_orig;
    }
    std::vector<T> x(un);
    for (std::size_t i = 0; i < un; ++i) {
        x[i] = lu.Dc[i] * x_orig[i];
    }
    return x;
}

} // namespace sparse_lu_detail

// ===========================================================================
// sparse_lu_make_baseline_factor_for_testing
// SLU-2 testing hook: constructs a factor object from a manually provided
// baseline_lu_storage.  NOT a numeric factorization; use only in tests.
// ===========================================================================
template <class T, class Index>
sparse_lu_factorization<T, Index>
sparse_lu_make_baseline_factor_for_testing(
    Index n,
    const baseline_lu_storage<T, Index>& storage)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_make_baseline_factor_for_testing: Index must be signed");

    // Validate structural properties before setting success=true
    sparse_lu_detail::validate_baseline_lu_storage_for_solve(storage, n);

    sparse_lu_factorization<T, Index> fac;
    fac.set_baseline_storage_for_internal_use_(n, storage);
    return fac;
}

#endif // VCP_TSPARSE_SPARSE_LU_SOLVE_IMPL_HPP
