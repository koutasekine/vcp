// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-1 conversion helpers — internal implementation.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER csc_storage<>,
// sparse_lu_symbolic_result<>, and all storage type definitions are in scope.
// It has no "namespace vcp { }" wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse_sparse_lu_convert.hpp>  (public wrapper)
//   <vcp/tsparse/tsparse.hpp>                    (umbrella)

#ifndef VCP_TSPARSE_SPARSE_LU_CONVERT_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_CONVERT_IMPL_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_convert.hpp>
#include <vcp/tsparse/tsparse_format.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

// ---------------------------------------------------------------------------
// 1.1  owning CSC validation helper
// ---------------------------------------------------------------------------
template <class T, class Index>
bool sparse_lu_is_valid_csc_storage(
    const csc_storage<T, Index>& A,
    Index nrow,
    Index ncol)
{
    if (nrow < Index(0) || ncol < Index(0)) return false;

    const std::size_t expected_cp = static_cast<std::size_t>(ncol) + 1u;
    if (A.col_ptr.size() != expected_cp) return false;
    if (A.col_ptr[0] != Index(0)) return false;

    for (Index j = Index(0); j < ncol; ++j) {
        if (A.col_ptr[static_cast<std::size_t>(j)] >
            A.col_ptr[static_cast<std::size_t>(j) + 1u]) {
            return false;
        }
    }

    const Index nnz = A.col_ptr[static_cast<std::size_t>(ncol)];
    if (nnz < Index(0)) return false;
    const std::size_t unnz = static_cast<std::size_t>(nnz);
    if (A.row_ind.size() != unnz) return false;
    if (A.values.size()  != unnz) return false;

    for (std::size_t k = 0; k < unnz; ++k) {
        if (A.row_ind[k] < Index(0) || A.row_ind[k] >= nrow) return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// 1.2  max_abs helper
// ---------------------------------------------------------------------------
template <class T, class Index>
typename vcp::tsparse_scalar::real_type<T>::type
sparse_lu_max_abs_csc(const csc_storage<T, Index>& A)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
    real_type m = real_type(0);
    const std::size_t n = A.values.size();
    for (std::size_t k = 0; k < n; ++k) {
        const real_type v = vcp::tsparse_scalar::abs_value(A.values[k]);
        if (v > m) m = v;
    }
    return m;
}

// ---------------------------------------------------------------------------
// 1.3  permutation helpers
//   Convention: perm[new] = old,  inv_perm[old] = new
// ---------------------------------------------------------------------------

template <class Index>
std::vector<Index> sparse_lu_identity_permutation(Index n)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_identity_permutation: negative n");
    }
    std::vector<Index> p(static_cast<std::size_t>(n));
    for (Index i = Index(0); i < n; ++i) {
        p[static_cast<std::size_t>(i)] = i;
    }
    return p;
}

template <class Index>
std::vector<Index> sparse_lu_inverse_permutation(
    const std::vector<Index>& perm)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");
    const Index n = static_cast<Index>(perm.size());
    std::vector<Index> inv(static_cast<std::size_t>(n), Index(-1));

    for (Index new_j = Index(0); new_j < n; ++new_j) {
        const std::size_t nj    = static_cast<std::size_t>(new_j);
        const Index        old_j = perm[nj];
        if (old_j < Index(0) || old_j >= n) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_inverse_permutation: perm value out of range");
        }
        const std::size_t oj = static_cast<std::size_t>(old_j);
        if (inv[oj] != Index(-1)) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_inverse_permutation: duplicate in permutation");
        }
        inv[oj] = new_j;
    }
    return inv;
}

// ---------------------------------------------------------------------------
// 1.4  column permutation application
//   col_perm[new_col] = old_col
//   A_perm[:, new_col] = A[:, col_perm[new_col]]
//
// SLU-1.1: validates that col_perm is a proper permutation (no duplicates)
// by calling sparse_lu_inverse_permutation, which throws invalid_argument
// for out-of-range or duplicate entries.
// ---------------------------------------------------------------------------
template <class T, class Index>
csc_storage<T, Index>
sparse_lu_apply_column_permutation_csc(
    const csc_storage<T, Index>& A,
    Index nrow,
    Index ncol,
    const std::vector<Index>& col_perm)
{
    (void)nrow;

    if (static_cast<Index>(col_perm.size()) != ncol) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_apply_column_permutation_csc: col_perm size mismatch");
    }

    // Validate: size, range, and uniqueness (rejects duplicates such as {0,0}).
    const std::vector<Index> inv_perm = sparse_lu_inverse_permutation(col_perm);
    (void)inv_perm;

    csc_storage<T, Index> P;
    P.col_ptr.resize(static_cast<std::size_t>(ncol) + 1u, Index(0));

    // Step 1: compute column sizes for P
    for (Index new_j = Index(0); new_j < ncol; ++new_j) {
        const std::size_t nj    = static_cast<std::size_t>(new_j);
        const Index        old_j = col_perm[nj];
        const std::size_t oj   = static_cast<std::size_t>(old_j);
        P.col_ptr[nj + 1u] = A.col_ptr[oj + 1u] - A.col_ptr[oj];
    }

    // Step 2: prefix sum → col_ptr
    for (std::size_t j = 0; j < static_cast<std::size_t>(ncol); ++j) {
        P.col_ptr[j + 1u] += P.col_ptr[j];
    }

    // Step 3: copy entries
    const Index nnz = P.col_ptr[static_cast<std::size_t>(ncol)];
    P.row_ind.resize(static_cast<std::size_t>(nnz));
    P.values.resize(static_cast<std::size_t>(nnz));

    for (Index new_j = Index(0); new_j < ncol; ++new_j) {
        const std::size_t nj    = static_cast<std::size_t>(new_j);
        const Index        old_j = col_perm[nj];
        const std::size_t oj   = static_cast<std::size_t>(old_j);
        Index out = P.col_ptr[nj];
        for (Index k = A.col_ptr[oj]; k < A.col_ptr[oj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const std::size_t so = static_cast<std::size_t>(out);
            P.row_ind[so] = A.row_ind[sk];
            P.values[so]  = A.values[sk];
            ++out;
        }
    }
    return P;
}

// ---------------------------------------------------------------------------
// 1.6  CSC conversion from Matrix
//
// Finalized path (CSR or CSC):
//   After finalize(), spmats stores data in CSR with deduplication and zero
//   removal already performed by normalize_coo() inside finalize().
//   tcsr_to_csc expands rows 0..nrow-1 in order; within each output CSC
//   column row indices are ascending.
//
// Unfinalized COO path (SLU-1.1 addition):
//   Copies COO buffers, sorts by (row,col), sums duplicates in-place.
//   Zero policy: does NOT call tcoo_remove_zeros so that structural zeros
//   arising from cancellation are preserved before the numeric phase.
//   Row order within each output CSC column is ascending (same as finalized
//   path) because tcoo_sort sorts by (row,col) lexicographic order and
//   tcoo_to_csc scatters in that order.
// ---------------------------------------------------------------------------
template <class Matrix>
csc_storage<
    typename Matrix::value_type,
    typename Matrix::index_type
>
sparse_lu_make_csc_storage(const Matrix& A)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");

    const Index nrow = static_cast<Index>(A.rowsize());
    const Index ncol = static_cast<Index>(A.columnsize());

    csc_storage<T, Index> C;

    if (A.is_finalized()) {
        const Index nnz = static_cast<Index>(A.stored_nnz());
        C.col_ptr.resize(static_cast<std::size_t>(ncol) + 1u, Index(0));
        C.row_ind.resize(static_cast<std::size_t>(nnz));
        C.values.resize(static_cast<std::size_t>(nnz));

        if (A.format() == vcp::sparse_csr) {
            vcp::tcsr_to_csc(
                nrow, ncol, nnz,
                A.outer_index().data(), A.inner_index().data(), A.values().data(),
                C.col_ptr.data(), C.row_ind.data(), C.values.data());
        }
        else if (A.format() == vcp::sparse_csc) {
            C.col_ptr = A.outer_index();
            C.row_ind = A.inner_index();
            C.values  = A.values();
        }
        else {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_make_csc_storage: unsupported finalized format");
        }
        return C;
    }

    // Unfinalized COO path
    std::vector<Index> row_v = A.coo_rows();
    std::vector<Index> col_v = A.coo_columns();
    std::vector<T>     val_v = A.coo_values();

    if (row_v.size() != col_v.size() || row_v.size() != val_v.size()) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_make_csc_storage: COO buffers have inconsistent sizes");
    }

    // Range validation (tcoo_to_csc also validates, but we give a clearer message)
    const std::size_t coo_sz = row_v.size();
    for (std::size_t k = 0u; k < coo_sz; ++k) {
        if (row_v[k] < Index(0) || row_v[k] >= nrow ||
            col_v[k] < Index(0) || col_v[k] >= ncol) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_make_csc_storage: COO entry index out of range");
        }
    }

    Index nnz_coo = static_cast<Index>(coo_sz);

    if (nnz_coo > Index(0)) {
        // Sort by (row,col): ensures ascending row order within each CSC column.
        vcp::tcoo_sort(nnz_coo, row_v.data(), col_v.data(), val_v.data());
        // Sum duplicates in-place; no zero removal (SLU-1 policy).
        nnz_coo = vcp::tcoo_sum_duplicates(
            nnz_coo, row_v.data(), col_v.data(), val_v.data());
    }

    C.col_ptr.resize(static_cast<std::size_t>(ncol) + 1u, Index(0));
    C.row_ind.resize(static_cast<std::size_t>(nnz_coo));
    C.values.resize(static_cast<std::size_t>(nnz_coo));

    if (nnz_coo > Index(0)) {
        vcp::tcoo_to_csc(nrow, ncol, nnz_coo,
                         row_v.data(), col_v.data(), val_v.data(),
                         C.col_ptr.data(), C.row_ind.data(), C.values.data());
    }
    return C;
}

// ---------------------------------------------------------------------------
// 1.7  conversion with column permutation from symbolic result
// ---------------------------------------------------------------------------
template <class Matrix>
csc_storage<
    typename Matrix::value_type,
    typename Matrix::index_type
>
sparse_lu_make_permuted_csc_storage(
    const Matrix& A,
    const sparse_lu_symbolic_result<typename Matrix::index_type>& sym)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;

    csc_storage<T, Index> C = sparse_lu_make_csc_storage(A);

    if (sym.col_perm.empty()) {
        return C;
    }

    const Index nrow = static_cast<Index>(A.rowsize());
    const Index ncol = static_cast<Index>(A.columnsize());
    return sparse_lu_apply_column_permutation_csc(C, nrow, ncol, sym.col_perm);
}

#endif // VCP_TSPARSE_SPARSE_LU_CONVERT_IMPL_HPP
