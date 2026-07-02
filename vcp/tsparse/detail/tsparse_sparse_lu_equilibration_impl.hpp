// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// O4.2a: Two-sided max-norm equilibration helper.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER all storage
// type definitions, sparse_lu_scalar_policy, and tsparse_sparse_lu_convert_impl.hpp
// are in scope.  It has no "namespace vcp { }" wrapper.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>

#ifndef VCP_TSPARSE_SPARSE_LU_EQUILIBRATION_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_EQUILIBRATION_IMPL_HPP

#include <cstddef>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// sparse_lu_compute_equilibration
//
// Two-sided max-norm scaling for a CSC matrix A (n x n).
// On exit:
//   Dr[i] = 1/max_j|A(i,j)| if row i is nonzero, else 1.
//   Dc[j] = 1/max_i|Dr[i]*A(i,j)| if col j is nonzero after row-scaling, else 1.
//
// Invariants:
//   - No NaN/Inf is produced: zero rows/cols yield scale 1.
//   - Purely T-arithmetic; no external BLAS; no mixed precision.
//   - Dr, Dc are in ORIGINAL (unpermuted) coordinates so the existing
//     solve_baseline_storage un-apply path (Step 1: rhs1=Dr*b, Step 6: x=Dc*x_orig)
//     works without modification.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_lu_compute_equilibration(
    const csc_storage<T, Index>& A,
    Index n,
    std::vector<T>& Dr,
    std::vector<T>& Dc)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    const std::size_t un = static_cast<std::size_t>(n);

    Dr.assign(un, T(1));
    Dc.assign(un, T(1));

    if (n <= Index(0)) return;

    // --- Row max-norm ---
    std::vector<real_type> row_max(un, real_type(0));
    for (std::size_t j = 0u; j < un; ++j) {
        for (Index k = A.col_ptr[j]; k < A.col_ptr[j + 1u]; ++k) {
            const std::size_t i = static_cast<std::size_t>(A.row_ind[k]);
            const real_type v = vcp::tsparse_scalar::abs_value(A.values[k]);
            if (v > row_max[i]) row_max[i] = v;
        }
    }
    for (std::size_t i = 0u; i < un; ++i) {
        if (row_max[i] > real_type(0)) {
            Dr[i] = T(real_type(1) / row_max[i]);
        }
    }

    // --- Column max-norm (after row scaling) ---
    std::vector<real_type> col_max(un, real_type(0));
    for (std::size_t j = 0u; j < un; ++j) {
        for (Index k = A.col_ptr[j]; k < A.col_ptr[j + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const std::size_t i  = static_cast<std::size_t>(A.row_ind[sk]);
            const real_type v = vcp::tsparse_scalar::abs_value(Dr[i] * A.values[sk]);
            if (v > col_max[j]) col_max[j] = v;
        }
    }
    for (std::size_t j = 0u; j < un; ++j) {
        if (col_max[j] > real_type(0)) {
            Dc[j] = T(real_type(1) / col_max[j]);
        }
    }
}

// ---------------------------------------------------------------------------
// sparse_lu_apply_equilibration
//
// Returns Aeq = diag(Dr) * A * diag(Dc) in CSC, sharing A's sparsity pattern.
// Dr[i] is the row-i scaling, Dc[j] the original-col-j scaling.
// Aeq[i,j] = Dr[i] * A[i,j] * Dc[j].
// The column permutation (ordering) is applied AFTER this by the caller.
// ---------------------------------------------------------------------------
template <class T, class Index>
csc_storage<T, Index>
sparse_lu_apply_equilibration(
    const csc_storage<T, Index>& A,
    const std::vector<T>& Dr,
    const std::vector<T>& Dc,
    Index n)
{
    csc_storage<T, Index> Aeq;
    Aeq.col_ptr = A.col_ptr;
    Aeq.row_ind = A.row_ind;
    Aeq.values.resize(A.values.size());

    for (std::size_t j = 0u; j < static_cast<std::size_t>(n); ++j) {
        const T dcj = Dc[j];
        for (Index k = A.col_ptr[j]; k < A.col_ptr[j + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const std::size_t i  = static_cast<std::size_t>(A.row_ind[sk]);
            Aeq.values[sk] = Dr[i] * A.values[sk] * dcj;
        }
    }
    return Aeq;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_EQUILIBRATION_IMPL_HPP
