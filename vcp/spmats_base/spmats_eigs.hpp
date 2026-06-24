// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATS_EIGS_HPP
#define VCP_SPMATS_EIGS_HPP

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include <vcp/spmats_base/spmats_eigs_types.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>
#include <vcp/tsparse/tsparse_arnoldi.hpp>
#include <vcp/tsparse/tsparse_factorization.hpp>
#include <vcp/tsparse/tsparse_preconditioner.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_generalized_shift_invert.hpp>
#include <vcp/tsparse/tsparse_b_inner_lanczos.hpp>
#include <vcp/tsparse/tsparse_hermitian_lanczos.hpp>
#include <vcp/tsparse/tsparse_dense_fallback.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_solvers.hpp>
// NOTE: spmats.hpp includes this file after spmats<_T,_Index> is defined.
// Do NOT #include <vcp/spmats.hpp> here to avoid circular dependency.

namespace vcp {

// tag dispatch helpers for b_inner_lanczos generalized eigs
struct spmats_b_inner_real_tag_ {};
struct spmats_b_inner_complex_tag_ {};

// ---------------------------------------------------------------------------
// 1. is_symmetric_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static bool is_symmetric_value_(const spmats<_T,_Index>& A, const typename vcp::tsparse_scalar::real_type<_T>::type& tol)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (tol <= scalar_real_type(0)) vcp::throw_error<vcp::invalid_argument>("spmats::is_symmetric: tol must be positive");
    if (A.rowsize() != A.columnsize()) return false;
    spmats<_T,_Index> C = A.as_csr();
    const std::vector<_Index>& outer = C.outer_index();
    const std::vector<_Index>& inner = C.inner_index();
    const std::vector<_T>& val = C.values();
    for (_Index i = 0; i < C.rowsize(); i++) {
        for (_Index p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
            const _Index j = inner[static_cast<std::size_t>(p)];
            if (i == j) continue;
            const _Index first = outer[static_cast<std::size_t>(j)];
            const _Index last  = outer[static_cast<std::size_t>(j + 1)];
            const typename std::vector<_Index>::const_iterator begin = inner.begin() + first;
            const typename std::vector<_Index>::const_iterator end   = inner.begin() + last;
            typename std::vector<_Index>::const_iterator it = std::lower_bound(begin, end, i);
            _T mirrored = _T(0);
            if (it != end && *it == i) {
                mirrored = val[static_cast<std::size_t>(it - inner.begin())];
            }
            if (vcp::tsparse_scalar::abs_value(val[static_cast<std::size_t>(p)] - mirrored) > tol) return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// 2. to_dense_impl_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::vector<std::vector<_T> > to_dense_impl_(const spmats<_T,_Index>& A)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type); // suppress unused typedef warning
    std::vector<std::vector<_T> > dense(
        static_cast<std::size_t>(A.rowsize()),
        std::vector<_T>(static_cast<std::size_t>(A.columnsize()), _T(0)));
    spmats<_T,_Index> C = A.as_csr();
    const std::vector<_Index>& outer = C.outer_index();
    const std::vector<_Index>& inner = C.inner_index();
    const std::vector<_T>& val = C.values();
    for (_Index i = 0; i < C.rowsize(); i++) {
        for (_Index p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
            dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(inner[static_cast<std::size_t>(p)])] =
                val[static_cast<std::size_t>(p)];
        }
    }
    return dense;
}

// ---------------------------------------------------------------------------
// 3. eig_method_to_string_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::string eig_method_to_string_(const eig_solver_method m)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    switch (m) {
    case eig_solver_method::lanczos:                  return "lanczos";
    case eig_solver_method::arnoldi:                  return "arnoldi";
    case eig_solver_method::shift_invert_lanczos:     return "shift_invert_lanczos";
    case eig_solver_method::shift_invert_arnoldi:     return "shift_invert_arnoldi";
    case eig_solver_method::dense_fallback_explicit:  return "dense_fallback_explicit";
    }
    return "unknown";
}

// ---------------------------------------------------------------------------
// 6. set_result_counts_  (declared before set_eig_diagnostics_ which calls it)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void set_result_counts_(eig_result<_T>& result, const std::size_t requested)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    result.requested_count = requested;
    result.returned_real_count = result.eigenvalues.size();
    result.returned_complex_count = 0;
    for (std::size_t i = 0; i < result.complex_eigenvalues.size(); i++) {
        if (vcp::tsparse_scalar::abs_value(result.complex_eigenvalues[i].imag())
                > std::numeric_limits<scalar_real_type>::epsilon()) {
            result.returned_complex_count++;
        }
    }
    result.returned_count = result.returned_real_count + result.returned_complex_count;
    if (result.converged_count == 0 && result.converged) result.converged_count = result.returned_count;
    if (!result.residuals_absolute.empty()) {
        result.residual_norm_absolute = *std::max_element(result.residuals_absolute.begin(), result.residuals_absolute.end());
    }
    if (!result.residuals_relative.empty()) {
        result.residual_norm_relative = *std::max_element(result.residuals_relative.begin(), result.residuals_relative.end());
    }
}

// ---------------------------------------------------------------------------
// 4. set_eig_diagnostics_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void set_eig_diagnostics_(eig_result<_T>& result,
                                  const eig_solver_method method,
                                  const std::size_t subspace_dim,
                                  const std::size_t matvec_count,
                                  const std::string& breakdown,
                                  const std::string& failure)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    result.method           = method;
    result.used_method      = eig_method_to_string_<_T,_Index>(method);
    result.used_subspace_dim    = subspace_dim;
    result.matrix_vector_products = matvec_count;
    result.breakdown_reason = breakdown;
    if (result.converged) {
        result.status   = "converged";
        result.message  = breakdown.empty() ? "converged" : breakdown;
        result.failure_reason.clear();
    }
    else {
        result.status         = "not_converged";
        result.failure_reason = failure.empty() ? "residual tolerance not reached" : failure;
        result.message        = result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, result.requested_count);
}

// ---------------------------------------------------------------------------
// 5. convert_dense_result_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> convert_dense_result_(
    const vcp::tsparse_dense_linalg::dense_eigen_result<_T>& source,
    const eig_solver_method method)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    eig_result<_T> result;
    result.eigenvalues         = source.eigenvalues;
    result.complex_eigenvalues = source.complex_eigenvalues;
    result.eigenvectors        = source.eigenvectors;
    result.converged           = source.converged;
    result.iterations          = source.iterations;
    result.residuals_absolute  = source.residuals;
    result.residual_norm_absolute = source.residual_norm;
    set_eig_diagnostics_<_T,_Index>(result, method, source.eigenvalues.size(), 0, std::string(), std::string());
    return result;
}

// ---------------------------------------------------------------------------
// 7. orthogonalization_to_string_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::string orthogonalization_to_string_(const orthogonalization_method method)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    switch (method) {
    case orthogonalization_method::modified_gram_schmidt:         return "modified_gram_schmidt";
    case orthogonalization_method::classical_gram_schmidt_twice:  return "classical_gram_schmidt_twice";
    }
    return "modified_gram_schmidt";
}

// ---------------------------------------------------------------------------
// 8. sort_eigenpairs_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void sort_eigenpairs_(eig_result<_T>& result)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (result.eigenvalues.empty()) return;  // no real eigenvalues to sort; preserve complex_eigenvalues
    std::vector<std::size_t> order(result.eigenvalues.size());
    for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
        return vcp::tsparse_scalar::real_part(result.eigenvalues[a]) <
               vcp::tsparse_scalar::real_part(result.eigenvalues[b]);
    });
    std::vector<_T> values(order.size());
    std::vector<typename eig_result<_T>::eigenvalue_type> cvalues;
    if (result.complex_eigenvalues.size() == result.eigenvalues.size()) cvalues.resize(order.size());
    std::vector<scalar_real_type> residuals_abs;
    if (result.residuals_absolute.size() == result.eigenvalues.size()) residuals_abs.resize(order.size());
    std::vector<scalar_real_type> residuals_rel;
    if (result.residuals_relative.size() == result.eigenvalues.size()) residuals_rel.resize(order.size());
    std::vector<std::vector<_T> > vectors;
    if (result.eigenvectors.size() == result.eigenvalues.size()) vectors.resize(order.size());
    for (std::size_t i = 0; i < order.size(); i++) {
        values[i] = result.eigenvalues[order[i]];
        if (!cvalues.empty())       cvalues[i]       = result.complex_eigenvalues[order[i]];
        if (!residuals_abs.empty()) residuals_abs[i] = result.residuals_absolute[order[i]];
        if (!residuals_rel.empty()) residuals_rel[i] = result.residuals_relative[order[i]];
        if (!vectors.empty())       vectors[i]       = result.eigenvectors[order[i]];
    }
    result.eigenvalues.swap(values);
    if (!cvalues.empty()) result.complex_eigenvalues.swap(cvalues);
    else if (!result.complex_eigenvalues.empty() && result.complex_eigenvalues.size() == result.eigenvalues.size())
        result.complex_eigenvalues.clear();
    if (!residuals_abs.empty()) result.residuals_absolute.swap(residuals_abs);
    else result.residuals_absolute.clear();
    if (!residuals_rel.empty()) result.residuals_relative.swap(residuals_rel);
    else result.residuals_relative.clear();
    if (!vectors.empty()) result.eigenvectors.swap(vectors);
}

// ---------------------------------------------------------------------------
// 10. populate_real_complex_eigenvalues_  (declared before select_eigenpairs_ which calls it)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void populate_real_complex_eigenvalues_(eig_result<_T>& result)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    if (!result.complex_eigenvalues.empty() && result.eigenvalues.empty()) return;
    result.complex_eigenvalues.clear();
    result.complex_eigenvalues.reserve(result.eigenvalues.size());
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
        result.complex_eigenvalues.push_back(
            typename eig_result<_T>::eigenvalue_type(result.eigenvalues[i]));
    }
}

// ---------------------------------------------------------------------------
// 9. select_eigenpairs_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void select_eigenpairs_(eig_result<_T>& result,
                                const std::size_t k,
                                const eig_target target,
                                const typename vcp::tsparse_scalar::real_type<_T>::type& shift
                                    = typename vcp::tsparse_scalar::real_type<_T>::type(0))
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::vector<std::size_t> order =
        vcp::tsparse_eigen_selection::select_real_eigenpairs(result.eigenvalues, k, target, shift);
    std::vector<_T> values;
    std::vector<std::vector<_T> > vectors;
    std::vector<scalar_real_type> residuals_abs;
    std::vector<scalar_real_type> residuals_rel;
    values.reserve(order.size());
    if (result.eigenvectors.size() == result.eigenvalues.size())       vectors.reserve(order.size());
    if (result.residuals_absolute.size() == result.eigenvalues.size()) residuals_abs.reserve(order.size());
    if (result.residuals_relative.size() == result.eigenvalues.size()) residuals_rel.reserve(order.size());
    for (std::size_t i = 0; i < order.size(); i++) {
        const std::size_t j = order[i];
        values.push_back(result.eigenvalues[j]);
        if (result.eigenvectors.size()       == result.eigenvalues.size()) vectors.push_back(result.eigenvectors[j]);
        if (result.residuals_absolute.size() == result.eigenvalues.size()) residuals_abs.push_back(result.residuals_absolute[j]);
        if (result.residuals_relative.size() == result.eigenvalues.size()) residuals_rel.push_back(result.residuals_relative[j]);
    }
    result.eigenvalues.swap(values);
    if (result.eigenvectors.size()       == vectors.size())       result.eigenvectors.swap(vectors);
    if (result.residuals_absolute.size() == residuals_abs.size()) result.residuals_absolute.swap(residuals_abs);
    if (result.residuals_relative.size() == residuals_rel.size()) result.residuals_relative.swap(residuals_rel);
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    set_result_counts_<_T,_Index>(result, k);
}

// ---------------------------------------------------------------------------
// 11. dot_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
dot_value_(const std::vector<_T>& a, const std::vector<_T>& b)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    return vcp::tsparse_scalar::real_dot_value(a, b);
}

// ---------------------------------------------------------------------------
// 12. norm_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
norm_value_(const std::vector<_T>& a)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    return vcp::tsparse_scalar::real_norm_value(a);
}

// ---------------------------------------------------------------------------
// 13. frobenius_norm_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
frobenius_norm_value_(const spmats<_T,_Index>& A)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> C = A.as_csr();
    scalar_real_type s(0);
    const std::vector<_T>& val = C.values();
    for (std::size_t i = 0; i < val.size(); i++) {
        const scalar_real_type a = vcp::tsparse_scalar::abs_value(val[i]);
        s += a * a;
    }
    return vcp::tsparse_scalar::sqrt_value(s);
}

// ---------------------------------------------------------------------------
// 14. eigenpair_residual_norm_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
eigenpair_residual_norm_value_(const spmats<_T,_Index>& A, const _T& lambda, const std::vector<_T>& v)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    std::vector<_T> r = A.mul_vec(v);
    for (std::size_t i = 0; i < r.size(); i++) r[i] -= lambda * v[i];
    return norm_value_<_T,_Index>(r);
}

// ---------------------------------------------------------------------------
// 15. eigenpair_residual_norm_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
eigenpair_residual_norm_(const spmats<_T,_Index>& A, const _T& lambda, const std::vector<_T>& v)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    return eigenpair_residual_norm_value_<_T,_Index>(A, lambda, v);
}

// ---------------------------------------------------------------------------
// 16. eigenpair_relative_residual_norm_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
eigenpair_relative_residual_norm_value_(const spmats<_T,_Index>& A, const _T& lambda, const std::vector<_T>& v)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const scalar_real_type abs_res = eigenpair_residual_norm_value_<_T,_Index>(A, lambda, v);
    const scalar_real_type vn      = norm_value_<_T,_Index>(v);
    const scalar_real_type denom   = frobenius_norm_value_<_T,_Index>(A) * vn
        + vcp::tsparse_scalar::abs_value(lambda) * vn
        + std::numeric_limits<scalar_real_type>::epsilon();
    return abs_res / denom;
}

// ---------------------------------------------------------------------------
// 17. max_eigenpair_residual_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
max_eigenpair_residual_value_(const spmats<_T,_Index>& A,
                               const std::vector<_T>& eigenvalues,
                               const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    scalar_real_type maximum(0);
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
        return std::numeric_limits<scalar_real_type>::infinity();
    }
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        const scalar_real_type residual = eigenpair_residual_norm_value_<_T,_Index>(A, eigenvalues[i], eigenvectors[i]);
        if (residual > maximum) maximum = residual;
    }
    return maximum;
}

// ---------------------------------------------------------------------------
// 18. eigenpair_residuals_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::vector<typename vcp::tsparse_scalar::real_type<_T>::type>
eigenpair_residuals_(const spmats<_T,_Index>& A,
                     const std::vector<_T>& eigenvalues,
                     const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::vector<scalar_real_type> residuals;
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
    residuals.reserve(eigenvalues.size());
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        residuals.push_back(eigenpair_residual_norm_<_T,_Index>(A, eigenvalues[i], eigenvectors[i]));
    }
    return residuals;
}

// ---------------------------------------------------------------------------
// 19. eigenpair_relative_residuals_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::vector<typename vcp::tsparse_scalar::real_type<_T>::type>
eigenpair_relative_residuals_(const spmats<_T,_Index>& A,
                               const std::vector<_T>& eigenvalues,
                               const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::vector<scalar_real_type> residuals;
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
    residuals.reserve(eigenvalues.size());
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        residuals.push_back(eigenpair_relative_residual_norm_value_<_T,_Index>(A, eigenvalues[i], eigenvectors[i]));
    }
    return residuals;
}

// ---------------------------------------------------------------------------
// 20. krylov_subspace_dim_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::size_t krylov_subspace_dim_(const std::size_t n, const std::size_t k, const std::size_t requested)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    return vcp::tsparse_eigensolvers::krylov_subspace_dim(n, k, requested);
}

// ---------------------------------------------------------------------------
// 21. check_dense_allowed_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void check_dense_allowed_(const spmats<_T,_Index>& A, const eig_options<_T>& options, const char*)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    vcp::tsparse_dense_fallback::check_dense_size(
        static_cast<std::size_t>(A.rowsize()),
        static_cast<std::size_t>(A.columnsize()),
        options.max_dense_size,
        options.allow_dense_conversion);
}

// ---------------------------------------------------------------------------
// 22. dense_eig_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> dense_eig_(std::vector<std::vector<_T> > dense, const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    eig_result<_T> result = vcp::tsparse_dense_linalg::is_dense_symmetric(dense, options.tol * scalar_real_type(10))
        ? convert_dense_result_<_T,_Index>(
              vcp::tsparse_dense_linalg::jacobi_eig_dense(dense, options.max_iter, options.tol),
              eig_solver_method::dense_fallback_explicit)
        : convert_dense_result_<_T,_Index>(
              vcp::tsparse_dense_linalg::qr_eig_dense(dense, options.max_iter, options.tol),
              eig_solver_method::dense_fallback_explicit);
    result.method           = eig_solver_method::dense_fallback_explicit;
    result.used_method      = eig_method_to_string_<_T,_Index>(eig_solver_method::dense_fallback_explicit);
    result.used_dense_fallback = true;
    return result;
}

// ---------------------------------------------------------------------------
// 23. is_diagonal_matrix_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static bool is_diagonal_matrix_(const spmats<_T,_Index>& A)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> C = A.as_csr();
    const std::vector<_Index>& outer = C.outer_index();
    const std::vector<_Index>& inner = C.inner_index();
    const std::vector<_T>& val       = C.values();
    for (std::size_t i = 0; i < static_cast<std::size_t>(C.rowsize()); i++) {
        for (_Index p = outer[i]; p < outer[i + 1]; p++) {
            const std::size_t j = static_cast<std::size_t>(inner[static_cast<std::size_t>(p)]);
            if (i != j && vcp::tsparse_scalar::abs_value(val[static_cast<std::size_t>(p)]) > scalar_real_type(0)) return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// 24. generalized_eigenpair_residual_norm_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
generalized_eigenpair_residual_norm_value_(const spmats<_T,_Index>& A,
                                            const spmats<_T,_Index>& B,
                                            const _T& eigenvalue,
                                            const std::vector<_T>& eigenvector)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    std::vector<_T> r  = A.mul_vec(eigenvector);
    std::vector<_T> bv = B.mul_vec(eigenvector);
    for (std::size_t i = 0; i < r.size(); i++) r[i] -= eigenvalue * bv[i];
    return norm_value_<_T,_Index>(r);
}

// ---------------------------------------------------------------------------
// 25. generalized_eigenpair_relative_residual_norm_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
generalized_eigenpair_relative_residual_norm_value_(const spmats<_T,_Index>& A,
                                                     const spmats<_T,_Index>& B,
                                                     const _T& eigenvalue,
                                                     const std::vector<_T>& eigenvector)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const scalar_real_type abs_res = generalized_eigenpair_residual_norm_value_<_T,_Index>(A, B, eigenvalue, eigenvector);
    const scalar_real_type vn      = norm_value_<_T,_Index>(eigenvector);
    const scalar_real_type denom   = frobenius_norm_value_<_T,_Index>(A) * vn
        + vcp::tsparse_scalar::abs_value(eigenvalue) * frobenius_norm_value_<_T,_Index>(B) * vn
        + std::numeric_limits<scalar_real_type>::epsilon();
    return abs_res / denom;
}

// ---------------------------------------------------------------------------
// 26. max_generalized_eigenpair_residual_value_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
max_generalized_eigenpair_residual_value_(const spmats<_T,_Index>& A,
                                           const spmats<_T,_Index>& B,
                                           const std::vector<_T>& eigenvalues,
                                           const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    scalar_real_type maximum(0);
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
        return std::numeric_limits<scalar_real_type>::infinity();
    }
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        const scalar_real_type residual =
            generalized_eigenpair_residual_norm_value_<_T,_Index>(A, B, eigenvalues[i], eigenvectors[i]);
        if (residual > maximum) maximum = residual;
    }
    return maximum;
}

// ---------------------------------------------------------------------------
// 27. generalized_eigenpair_residuals_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::vector<typename vcp::tsparse_scalar::real_type<_T>::type>
generalized_eigenpair_residuals_(const spmats<_T,_Index>& A,
                                  const spmats<_T,_Index>& B,
                                  const std::vector<_T>& eigenvalues,
                                  const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::vector<scalar_real_type> residuals;
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
    residuals.reserve(eigenvalues.size());
    for (std::size_t p = 0; p < eigenvalues.size(); p++) {
        residuals.push_back(generalized_eigenpair_residual_norm_value_<_T,_Index>(A, B, eigenvalues[p], eigenvectors[p]));
    }
    return residuals;
}

// ---------------------------------------------------------------------------
// 28. generalized_eigenpair_relative_residuals_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::vector<typename vcp::tsparse_scalar::real_type<_T>::type>
generalized_eigenpair_relative_residuals_(const spmats<_T,_Index>& A,
                                           const spmats<_T,_Index>& B,
                                           const std::vector<_T>& eigenvalues,
                                           const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::vector<scalar_real_type> residuals;
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
    residuals.reserve(eigenvalues.size());
    for (std::size_t p = 0; p < eigenvalues.size(); p++) {
        residuals.push_back(
            generalized_eigenpair_relative_residual_norm_value_<_T,_Index>(A, B, eigenvalues[p], eigenvectors[p]));
    }
    return residuals;
}

// ---------------------------------------------------------------------------
// lanczos_package_to_result_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class ApplyA>
static eig_result<_T> lanczos_package_to_result_(
    const vcp::tsparse_lanczos::lanczos_result_package<_T, ApplyA>& pkg,
    const spmats<_T,_Index>& self,
    const std::size_t k,
    const eig_solver_method method)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    eig_result<_T> result;
    result.eigenvalues               = pkg.eigenvalues;
    result.eigenvectors              = pkg.eigenvectors;
    result.converged                 = pkg.converged;
    result.iterations                = pkg.iterations;
    result.converged_count           = pkg.converged_count;
    result.returned_count            = pkg.returned_count;
    result.matrix_vector_products    = pkg.mv_count;
    result.residuals_absolute        = pkg.residuals_abs;
    result.residuals_relative        = pkg.residuals_rel;
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.breakdown_reason          = pkg.breakdown_reason;
    result.failure_reason            = pkg.failure_reason;
    if (!pkg.residuals_abs.empty()) {
        result.residual_norm_absolute =
            *std::max_element(pkg.residuals_abs.begin(), pkg.residuals_abs.end());
    }
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    if (!result.eigenvectors.empty()) {
        spmats<_T,_Index> A = self.as_csr();
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    set_eig_diagnostics_<_T,_Index>(result, method, result.returned_count, result.matrix_vector_products,
        result.breakdown_reason, result.failure_reason);
    result.used_shift_invert = (method == eig_solver_method::shift_invert_lanczos
                             || method == eig_solver_method::shift_invert_arnoldi);
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// validate_eig_input_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void validate_eig_input_(const spmats<_T,_Index>& A, const char* routine)
{
    if (A.rowsize() != A.columnsize())
        vcp::throw_error<vcp::dimension_error>(routine, ": matrix must be square");
}

// ---------------------------------------------------------------------------
// validate_eigs_input_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void validate_eigs_input_(const spmats<_T,_Index>& A, const std::size_t k,
                                  const char* routine, const eig_solver_method = eig_solver_method::lanczos)
{
    validate_eig_input_<_T,_Index>(A, routine);
    if (k == 0) vcp::throw_error<vcp::invalid_argument>(routine, ": k must be positive");
}

// ---------------------------------------------------------------------------
// validate_generalized_eig_input_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void validate_generalized_eig_input_(const spmats<_T,_Index>& A,
                                             const spmats<_T,_Index>& B,
                                             const char* routine)
{
    if (A.rowsize() != A.columnsize() || B.rowsize() != B.columnsize())
        vcp::throw_error<vcp::dimension_error>(routine, ": matrices must be square");
    if (A.rowsize() != B.rowsize())
        vcp::throw_error<vcp::dimension_error>(routine, ": matrix sizes must match");
}

// ---------------------------------------------------------------------------
// resolve_eigs_options_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_options<_T> resolve_eigs_options_(const spmats<_T,_Index>& A, const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    eig_options<_T> active = options;
    if (active.structure == matrix_structure_hint::hermitian) {
        vcp::throw_error<vcp::state_error>("spmats::eigs: hermitian structure hint is reserved for future complex support");
        return active;
    }
    if (active.structure == matrix_structure_hint::symmetric) {
        if (A.rowsize() != A.columnsize()) {
            vcp::throw_error<vcp::dimension_error>("spmats::eigs: symmetric structure hint requires a square matrix");
            return active;
        }
        if (active.method == eig_solver_method::arnoldi) active.method = eig_solver_method::lanczos;
    }
    else if (active.structure == matrix_structure_hint::general) {
        if (active.method == eig_solver_method::lanczos) active.method = eig_solver_method::arnoldi;
    }
    if (active.use_shift) {
        if (active.method == eig_solver_method::lanczos) active.method = eig_solver_method::shift_invert_lanczos;
        else if (active.method == eig_solver_method::arnoldi) active.method = eig_solver_method::shift_invert_arnoldi;
        else if (active.method == eig_solver_method::dense_fallback_explicit) {
            vcp::throw_error<vcp::invalid_argument>("spmats::eigs: shift is not supported by dense_fallback_explicit");
        }
    }
    if (active.structure == matrix_structure_hint::auto_detect
     && active.method == eig_solver_method::lanczos
     && !is_symmetric_value_<_T,_Index>(A, vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(10))) {
        vcp::throw_error<vcp::domain_error>("spmats::eigs: Lanczos requires a symmetric matrix");
    }
    return active;
}

// ---------------------------------------------------------------------------
// generalized_diagonal_eigs_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> generalized_diagonal_eigs_(const spmats<_T,_Index>& A,
                                                   const spmats<_T,_Index>& B,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(A.rowsize());
    eig_result<_T> result;
    result.requested_count = k;
    result.method = options.method;
    result.used_method = eig_method_to_string_<_T,_Index>(options.method);
    result.used_generalized_operator = true;
    result.used_dense_fallback = false;
    std::vector<_T> values;
    std::vector<std::vector<_T> > vectors;
    values.reserve(n);
    vectors.reserve(n);
    const scalar_real_type pivot_tol = std::numeric_limits<scalar_real_type>::epsilon();
    for (std::size_t i = 0; i < n; i++) {
        const _T bdiag = B.get(static_cast<_Index>(i), static_cast<_Index>(i));
        if (vcp::tsparse_scalar::abs_value(bdiag) <= pivot_tol) {
            vcp::throw_error<vcp::numerical_error>("spmats::eigs(A,B): zero diagonal in B");
        }
        values.push_back(A.get(static_cast<_Index>(i), static_cast<_Index>(i)) / bdiag);
        std::vector<_T> e(n, _T(0));
        e[i] = _T(1);
        vectors.push_back(e);
    }
    const std::vector<std::size_t> order =
        vcp::tsparse_eigen_selection::select_real_eigenpairs(values, k, options.target, options.shift);
    for (std::size_t i = 0; i < order.size(); i++) {
        result.eigenvalues.push_back(values[order[i]]);
        result.eigenvectors.push_back(vectors[order[i]]);
    }
    result.residuals_absolute = generalized_eigenpair_residuals_<_T,_Index>(A, B, result.eigenvalues, result.eigenvectors);
    result.residuals_relative = generalized_eigenpair_relative_residuals_<_T,_Index>(A, B, result.eigenvalues, result.eigenvectors);
    result.converged = true;
    result.status = "converged";
    result.message = "converged";
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// hermitian residual helpers (Phase 7)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
hermitian_eigenpair_residual_norm_value_(const spmats<_T,_Index>& A, const _T& lambda, const std::vector<_T>& v)
{
    std::vector<_T> r = A.mul_vec(v);
    for (std::size_t i = 0; i < r.size(); i++) r[i] -= lambda * v[i];
    return vcp::tsparse_hermitian_lanczos::hermitian_norm(r);
}

template <typename _T, typename _Index>
static std::vector<typename vcp::tsparse_scalar::real_type<_T>::type>
hermitian_eigenpair_residuals_(const spmats<_T,_Index>& A,
                                const std::vector<_T>& eigenvalues,
                                const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::vector<scalar_real_type> residuals;
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
    residuals.reserve(eigenvalues.size());
    for (std::size_t i = 0; i < eigenvalues.size(); i++)
        residuals.push_back(hermitian_eigenpair_residual_norm_value_<_T,_Index>(A, eigenvalues[i], eigenvectors[i]));
    return residuals;
}

template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
hermitian_eigenpair_relative_residual_norm_value_(const spmats<_T,_Index>& A, const _T& lambda, const std::vector<_T>& v)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const scalar_real_type abs_res = hermitian_eigenpair_residual_norm_value_<_T,_Index>(A, lambda, v);
    const scalar_real_type vn = vcp::tsparse_hermitian_lanczos::hermitian_norm(v);
    const scalar_real_type denom = frobenius_norm_value_<_T,_Index>(A) * vn
        + vcp::tsparse_scalar::abs_value(lambda) * vn
        + std::numeric_limits<scalar_real_type>::epsilon();
    return abs_res / denom;
}

template <typename _T, typename _Index>
static std::vector<typename vcp::tsparse_scalar::real_type<_T>::type>
hermitian_eigenpair_relative_residuals_(const spmats<_T,_Index>& A,
                                         const std::vector<_T>& eigenvalues,
                                         const std::vector<std::vector<_T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::vector<scalar_real_type> residuals;
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
    residuals.reserve(eigenvalues.size());
    for (std::size_t i = 0; i < eigenvalues.size(); i++)
        residuals.push_back(hermitian_eigenpair_relative_residual_norm_value_<_T,_Index>(A, eigenvalues[i], eigenvectors[i]));
    return residuals;
}

// ---------------------------------------------------------------------------
// check_hermitian_sparse_  (uses *this -> takes A by const ref)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static std::string check_hermitian_sparse_(const spmats<_T,_Index>& A,
                                            const typename vcp::tsparse_scalar::real_type<_T>::type& tol)
{
    if (A.rowsize() != A.columnsize()) return "non-square matrix";
    spmats<_T,_Index> C = A.as_csr();
    const std::vector<_Index>& outer = C.outer_index();
    const std::vector<_Index>& inner = C.inner_index();
    const std::vector<_T>& val = C.values();
    for (_Index i = 0; i < C.rowsize(); i++) {
        for (_Index p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
            const _Index j = inner[static_cast<std::size_t>(p)];
            const _T v = val[static_cast<std::size_t>(p)];
            if (i == j) {
                const _T diff = v - vcp::tsparse_scalar::conjugate_if_needed(v);
                if (vcp::tsparse_scalar::abs_value(diff) > tol)
                    return "non-real diagonal element";
                continue;
            }
            const _Index first = outer[static_cast<std::size_t>(j)];
            const _Index last  = outer[static_cast<std::size_t>(j + 1)];
            const typename std::vector<_Index>::const_iterator begin = inner.begin() + first;
            const typename std::vector<_Index>::const_iterator end   = inner.begin() + last;
            typename std::vector<_Index>::const_iterator it = std::lower_bound(begin, end, i);
            _T mirrored = _T(0);
            if (it != end && *it == i)
                mirrored = val[static_cast<std::size_t>(it - inner.begin())];
            const _T expected = vcp::tsparse_scalar::conjugate_if_needed(mirrored);
            if (vcp::tsparse_scalar::abs_value(v - expected) > tol)
                return "non-Hermitian entry found";
        }
    }
    return "";
}

// ---------------------------------------------------------------------------
// hermitian_package_to_result_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> hermitian_package_to_result_(
    const vcp::tsparse_hermitian_lanczos::hermitian_lanczos_result<_T>& pkg,
    const spmats<_T,_Index>& self,
    const std::size_t k)
{
    eig_result<_T> result;
    result.eigenvalues               = pkg.eigenvalues;
    result.eigenvectors              = pkg.eigenvectors;
    result.converged                 = pkg.converged;
    result.iterations                = pkg.iterations;
    result.converged_count           = pkg.converged_count;
    result.returned_count            = pkg.returned_count;
    result.matrix_vector_products    = pkg.mv_count;
    result.residuals_absolute        = pkg.residuals_abs;
    result.residuals_relative        = pkg.residuals_rel;
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.breakdown_reason          = pkg.breakdown_reason;
    result.failure_reason            = pkg.failure_reason;
    result.used_method               = pkg.used_method;
    result.method                    = eig_solver_method::lanczos;
    result.used_shift_invert         = false;
    result.used_dense_fallback       = false;
    result.used_generalized_operator = false;
    result.requested_count           = k;
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    if (!result.eigenvectors.empty()) {
        spmats<_T,_Index> A = self.as_csr();
        result.residuals_absolute = hermitian_eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = hermitian_eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    (void)sizeof(scalar_real_type);
    if (!result.residuals_absolute.empty())
        result.residual_norm_absolute = *std::max_element(result.residuals_absolute.begin(), result.residuals_absolute.end());
    if (!result.residuals_relative.empty())
        result.residual_norm_relative = *std::max_element(result.residuals_relative.begin(), result.residuals_relative.end());
    if (result.converged) {
        result.status  = "converged";
        result.message = result.breakdown_reason.empty() ? "converged" : result.breakdown_reason;
    } else {
        result.status  = "not_converged";
        result.message = result.failure_reason.empty() ? "not converged" : result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// complex_generalized_eigs_unsupported_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> complex_generalized_eigs_unsupported_(const std::size_t k, const eig_options<_T>& options)
{
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported_complex_generalized";
    result.failure_reason = "generalized eigenproblem (A,B) is not supported for complex scalar types in Phase 7";
    result.message = result.failure_reason;
    result.method = options.method;
    result.used_method = eig_method_to_string_<_T,_Index>(options.method);
    result.used_generalized_operator = true;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// complex_preconditioned_eigs_unsupported_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> complex_preconditioned_eigs_unsupported_(const std::size_t k, const eig_options<_T>& options)
{
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.method = options.method;
    result.used_method = eig_method_to_string_<_T,_Index>(options.method);
    const bool wants_shift_invert =
        (options.method == eig_solver_method::shift_invert_lanczos)
        || (options.method == eig_solver_method::shift_invert_arnoldi)
        || options.use_shift;
    if (wants_shift_invert) {
        result.status = "unsupported_complex_shift_invert";
        result.failure_reason = "shift-invert with preconditioner is not supported for complex scalar types";
        result.used_shift_invert = true;
    } else {
        result.status = "unsupported_preconditioner_for_method";
        result.failure_reason = "preconditioner overload is not supported for complex scalar types in Phase 7";
    }
    result.message = result.failure_reason;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// hermitian_lanczos_eigs_impl_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> hermitian_lanczos_eigs_impl_(const spmats<_T,_Index>& self,
                                                     const std::size_t k,
                                                     const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    const R normA_fro = frobenius_norm_value_<_T,_Index>(A);
    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + std::size_t(5), std::min(n, std::size_t(30)))
        : options.subspace_dim;
    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    } apply_op = { &A };
    const std::size_t max_restarts_l = (sdim > 0) ? (options.max_iter / sdim + k + 1) : options.max_iter;
    vcp::tsparse_hermitian_lanczos::hermitian_lanczos_result<_T> pkg =
        vcp::tsparse_hermitian_lanczos::hermitian_lanczos_eigs<_T, apply_fn>(
            n, k, sdim, max_restarts_l, options.tol,
            options.random_seed, options.random_start,
            options.target, R(options.shift),
            normA_fro, options.compute_residual_history, apply_op);
    return hermitian_package_to_result_<_T,_Index>(pkg, self, k);
}

// ---------------------------------------------------------------------------
// complex_standard_eigs_with_info_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> complex_standard_eigs_with_info_(const spmats<_T,_Index>& self,
                                                         const std::size_t k,
                                                         const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (k == 0) {
        eig_result<_T> result;
        result.requested_count = 0;
        result.converged = true;
        result.status = "converged";
        result.message = "converged";
        result.method = eig_solver_method::lanczos;
        result.used_method = "hermitian_lanczos";
        set_result_counts_<_T,_Index>(result, 0);
        return result;
    }
    if (self.rowsize() != self.columnsize())
        vcp::throw_error<vcp::dimension_error>("spmats::eigs(complex): matrix must be square");
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    if (n == 0) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "not_converged";
        result.failure_reason = "matrix is empty (n == 0) but k > 0";
        result.message = result.failure_reason;
        result.method = options.method;
        result.used_method = eig_method_to_string_<_T,_Index>(options.method);
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    const bool wants_shift_invert =
        (options.method == eig_solver_method::shift_invert_lanczos)
        || (options.method == eig_solver_method::shift_invert_arnoldi)
        || ((options.use_shift)
            && (options.method == eig_solver_method::lanczos
                || options.method == eig_solver_method::arnoldi));
    if (wants_shift_invert) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_complex_shift_invert";
        result.failure_reason = "shift-invert is not supported for complex scalar types in Phase 7";
        result.message = result.failure_reason;
        result.method = options.method;
        result.used_method = eig_method_to_string_<_T,_Index>(options.method);
        result.used_shift_invert = true;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    if (options.method == eig_solver_method::arnoldi) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_complex_non_hermitian";
        result.failure_reason = "Arnoldi method is not supported for complex scalar types; use Lanczos with hermitian structure hint";
        result.message = result.failure_reason;
        result.method = options.method;
        result.used_method = "arnoldi";
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    if (options.structure == matrix_structure_hint::general) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_complex_non_hermitian";
        result.failure_reason = "general matrix structure is not supported for complex scalar types; Phase 7 supports only complex Hermitian";
        result.message = result.failure_reason;
        result.method = options.method;
        result.used_method = eig_method_to_string_<_T,_Index>(options.method);
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    if (options.method == eig_solver_method::dense_fallback_explicit) {
        if (options.max_iter == 0 || options.tol <= scalar_real_type(0))
            vcp::throw_error<vcp::invalid_argument>("spmats::eigs(complex, dense): invalid iteration option");
        check_dense_allowed_<_T,_Index>(self, options, "spmats::eigs(complex)");
        std::vector<std::vector<_T> > dense = to_dense_impl_<_T,_Index>(self);
        eig_result<_T> result = dense_eig_<_T,_Index>(dense, options);
        result.used_dense_fallback = true;
        result.requested_count = k;
        select_eigenpairs_<_T,_Index>(result, k, options.target, options.shift);
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    if (options.max_iter == 0)
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs(complex): max_iter must be positive");
    if (options.tol <= scalar_real_type(0))
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs(complex): tol must be positive");
    const scalar_real_type herm_tol_candidate = options.tol * scalar_real_type(100);
    const scalar_real_type herm_tol = (herm_tol_candidate > scalar_real_type(0))
        ? herm_tol_candidate
        : vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(10);
    const std::string herm_check = check_hermitian_sparse_<_T,_Index>(self, herm_tol);
    const bool is_herm = herm_check.empty();
    if (!is_herm) {
        if (options.structure == matrix_structure_hint::hermitian
         || options.structure == matrix_structure_hint::symmetric) {
            eig_result<_T> result;
            result.requested_count = k;
            result.converged = false;
            result.status = "non_hermitian";
            result.failure_reason = "Hermitian structure asserted but check failed: " + herm_check;
            result.message = result.failure_reason;
            result.method = options.method;
            result.used_method = "hermitian_lanczos";
            set_result_counts_<_T,_Index>(result, k);
            return result;
        }
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_complex_non_hermitian";
        result.failure_reason = "complex matrix is not Hermitian; non-Hermitian complex eigs not supported in Phase 7";
        result.message = result.failure_reason;
        result.method = options.method;
        result.used_method = "hermitian_lanczos";
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    const std::size_t k_eff = (k < n) ? k : n;
    eig_result<_T> result = hermitian_lanczos_eigs_impl_<_T,_Index>(self, k_eff, options);
    result.requested_count = k;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// shift_invert_lanczos_eigs_  (standard, no user preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> shift_invert_lanczos_eigs_(const spmats<_T,_Index>& self,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options,
                                                   const typename vcp::tsparse_scalar::real_type<_T>::type& sigma)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());

    spmats<_T,_Index> shifted = A;
    for (std::size_t i = 0; i < n; i++) {
        const _T cur = shifted.get(static_cast<_Index>(i), static_cast<_Index>(i));
        shifted.set(static_cast<_Index>(i), static_cast<_Index>(i), cur - _T(sigma));
    }
    shifted.finalize();
    spmats<_T,_Index> shiftedCSR = shifted.as_csr();

    typedef vcp::tsparse_factorization::ilu0_data<_T, _Index> ILU;
    ILU ilu = vcp::tsparse_factorization::ilu0_factorize<_T, _Index>(
        shiftedCSR.outer_index(), shiftedCSR.inner_index(), shiftedCSR.values(), n,
        vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(14));
    if (ilu.singular_or_unstable) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = eig_solver_method::shift_invert_lanczos;
        result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_lanczos);
        result.used_shift_invert = true;
        result.used_dense_fallback = false;
        result.status = "factorization_failed";
        result.failure_reason = "ILU zero or near-zero pivot detected";
        result.message = result.failure_reason;
        result.factorization_diagnostics = ilu.diagnostics;
        result.factorization_zero_pivots = ilu.zero_pivots;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const std::size_t inner_max = std::min(n, std::size_t(100));
    const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
    std::size_t linear_solve_count = 0;
    std::size_t inner_failure_count = 0;
    std::size_t inner_iteration_count = 0;
    scalar_real_type inner_residual_norm = scalar_real_type(0);

    struct apply_fn {
        const spmats<_T,_Index>* mat;
        const ILU* ilu_ptr;
        std::size_t inner_max;
        scalar_real_type inner_tol;
        scalar_real_type sigma;
        std::size_t* solve_count;
        std::size_t* failure_count;
        std::size_t* iteration_count;
        scalar_real_type* max_inner_residual;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
            const std::size_t nn = x.size();
            struct av {
                const spmats<_T,_Index>* m;
                scalar_real_type s;
                void operator()(const std::vector<_T>& u, std::vector<_T>& v) const {
                    v = m->mul_vec(u);
                    for (std::size_t i = 0; i < v.size(); i++) v[i] -= _T(s) * u[i];
                }
            } av_fn = {mat, sigma};
            struct pv {
                const ILU* p;
                void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
                    z = vcp::tsparse_factorization::ilu0_solve(*p, r);
                }
            } pv_fn = {ilu_ptr};
            auto gr = vcp::tsparse_factorization::gmres_solve<_T, av, pv>(
                av_fn, pv_fn, x, inner_max / 10 + 1, inner_tol, std::min(nn, inner_max));
            (*solve_count)++;
            (*iteration_count) += gr.iterations;
            if (gr.residual_norm > *max_inner_residual) *max_inner_residual = gr.residual_norm;
            if (!gr.converged) (*failure_count)++;
            y = gr.x;
        }
    } apply_si = { &A, &ilu, inner_max, inner_tol, sigma, &linear_solve_count,
        &inner_failure_count, &inner_iteration_count, &inner_residual_norm };

    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const std::size_t max_restarts_si = (sdim > 0) ? (options.max_iter / sdim + k + 1) : options.max_iter;
    auto pkg = vcp::tsparse_lanczos::lanczos_eigs_standard<_T, apply_fn>(
        n, k, sdim, max_restarts_si, options.tol,
        options.random_seed, options.random_start,
        eig_target::largest_magnitude, scalar_real_type(0), options.compute_residual_history, apply_si);

    for (std::size_t i = 0; i < pkg.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(pkg.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0)) {
            pkg.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
        }
    }
    eig_result<_T> result = lanczos_package_to_result_<_T,_Index>(pkg, self, k, eig_solver_method::shift_invert_lanczos);
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_failure_count = inner_failure_count;
    result.inner_residual_norm = inner_residual_norm;
    result.factorization_diagnostics = ilu.diagnostics;
    result.factorization_zero_pivots = ilu.zero_pivots;
    result.used_shift_invert = true;
    if (inner_failure_count != 0) {
        result.converged = false;
        result.status = "inner_solve_failed";
        result.failure_reason = "inner GMRES solve failed";
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// lanczos_eigs_new_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> lanczos_eigs_new_(const spmats<_T,_Index>& self,
                                          const std::size_t k,
                                          const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const scalar_real_type shift_val = options.shift;

    if (options.method == eig_solver_method::shift_invert_lanczos) {
        return shift_invert_lanczos_eigs_<_T,_Index>(self, k, options, shift_val);
    }

    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    } apply = { &A };

    const std::size_t max_restarts_l = (sdim > 0) ? (options.max_iter / sdim + k + 1) : options.max_iter;
    auto pkg = vcp::tsparse_lanczos::lanczos_eigs_standard<_T, apply_fn>(
        n, k, sdim, max_restarts_l, options.tol,
        options.random_seed, options.random_start,
        options.target, shift_val, options.compute_residual_history, apply);

    return lanczos_package_to_result_<_T,_Index>(pkg, self, k, eig_solver_method::lanczos);
}

// ---------------------------------------------------------------------------
// shift_invert_arnoldi_eigs_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> shift_invert_arnoldi_eigs_(const spmats<_T,_Index>& self,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options,
                                                   const typename vcp::tsparse_scalar::real_type<_T>::type& sigma)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    spmats<_T,_Index> shifted = A;
    for (std::size_t i = 0; i < n; i++) {
        const _T cur = shifted.get(static_cast<_Index>(i), static_cast<_Index>(i));
        shifted.set(static_cast<_Index>(i), static_cast<_Index>(i), cur - _T(sigma));
    }
    shifted.finalize();
    spmats<_T,_Index> shiftedCSR = shifted.as_csr();
    typedef vcp::tsparse_factorization::ilu0_data<_T, _Index> ILU;
    ILU ilu = vcp::tsparse_factorization::ilu0_factorize<_T, _Index>(
        shiftedCSR.outer_index(), shiftedCSR.inner_index(), shiftedCSR.values(), n,
        vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(14));
    if (ilu.singular_or_unstable) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = eig_solver_method::shift_invert_arnoldi;
        result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_arnoldi);
        result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
        result.used_shift_invert = true;
        result.used_dense_fallback = false;
        result.status = "factorization_failed";
        result.failure_reason = "ILU zero or near-zero pivot detected";
        result.message = result.failure_reason;
        result.factorization_diagnostics = ilu.diagnostics;
        result.factorization_zero_pivots = ilu.zero_pivots;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    const std::size_t inner_max = std::min(n, std::size_t(100));
    const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
    std::size_t linear_solve_count = 0;
    std::size_t inner_failure_count = 0;
    std::size_t inner_iteration_count = 0;
    scalar_real_type inner_residual_norm = scalar_real_type(0);
    struct apply_fn {
        const spmats<_T,_Index>* mat;
        const ILU* ilu_ptr;
        std::size_t inner_max;
        scalar_real_type inner_tol;
        scalar_real_type sigma;
        std::size_t* solve_count;
        std::size_t* failure_count;
        std::size_t* iteration_count;
        scalar_real_type* max_inner_residual;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
            const std::size_t nn = x.size();
            struct av {
                const spmats<_T,_Index>* m;
                scalar_real_type s;
                void operator()(const std::vector<_T>& u, std::vector<_T>& v) const {
                    v = m->mul_vec(u);
                    for (std::size_t i = 0; i < v.size(); i++) v[i] -= _T(s) * u[i];
                }
            } av_fn = {mat, sigma};
            struct pv {
                const ILU* p;
                void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
                    z = vcp::tsparse_factorization::ilu0_solve(*p, r);
                }
            } pv_fn = {ilu_ptr};
            auto gr = vcp::tsparse_factorization::gmres_solve<_T, av, pv>(
                av_fn, pv_fn, x, inner_max / 10 + 1, inner_tol, std::min(nn, inner_max));
            (*solve_count)++;
            (*iteration_count) += gr.iterations;
            if (gr.residual_norm > *max_inner_residual) *max_inner_residual = gr.residual_norm;
            if (!gr.converged) (*failure_count)++;
            y = gr.x;
        }
    } apply_si = { &A, &ilu, inner_max, inner_tol, sigma, &linear_solve_count,
        &inner_failure_count, &inner_iteration_count, &inner_residual_norm };
    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
        vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
            n, k, sdim, options.max_iter + options.max_iter * k,
            options.tol, options.orthogonalization, true, true,
            options.random_seed, options.random_start, eig_target::largest_magnitude, scalar_real_type(0),
            options.compute_residual_history, apply_si);
    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::shift_invert_arnoldi;
    result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_arnoldi);
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_failure_count = inner_failure_count;
    result.inner_residual_norm = inner_residual_norm;
    result.factorization_diagnostics = ilu.diagnostics;
    result.factorization_zero_pivots = ilu.zero_pivots;
    result.used_subspace_dim = sdim;
    result.breakdown_reason = pkg.breakdown_reason;
    result.failure_reason = pkg.failure_reason;
    result.converged_count = pkg.converged_count;
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.eigenvalues = pkg.eigenvalues;
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(result.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            result.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }
    result.eigenvectors = pkg.eigenvectors;
    result.used_shift_invert = true;
    result.used_dense_fallback = false;
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    result.converged = pkg.converged && result.eigenvalues.size() >= k && inner_failure_count == 0;
    if (inner_failure_count != 0) {
        result.status = "inner_solve_failed";
        result.failure_reason = "inner GMRES solve failed";
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
    } else {
        set_eig_diagnostics_<_T,_Index>(result, eig_solver_method::shift_invert_arnoldi, sdim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// arnoldi_eigs_new_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> arnoldi_eigs_new_(const spmats<_T,_Index>& self,
                                          const std::size_t k,
                                          const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    if (options.method == eig_solver_method::shift_invert_arnoldi) {
        return shift_invert_arnoldi_eigs_<_T,_Index>(self, k, options, options.shift);
    }
    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const std::size_t max_restarts = options.max_iter + static_cast<std::size_t>(options.max_iter * k);

    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    } apply_op = { &A };

    const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
        vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
            n, k, sdim, max_restarts, options.tol, options.orthogonalization,
            true, true, options.random_seed, options.random_start,
            options.target, vcp::tsparse_scalar::real_part(options.shift),
            options.compute_residual_history, apply_op);

    eig_result<_T> result;
    result.method = eig_solver_method::arnoldi;
    result.used_method = "arnoldi";
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.used_subspace_dim = sdim;
    result.breakdown_reason = pkg.breakdown_reason;
    result.failure_reason = pkg.failure_reason;
    result.converged_count = pkg.converged_count;
    result.returned_count = pkg.returned_count;
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.eigenvalues = pkg.eigenvalues;
    result.eigenvectors = pkg.eigenvectors;
    for (std::size_t i = 0; i < pkg.complex_eigenvalues.size(); i++) {
        result.complex_eigenvalues.push_back(
            typename eig_result<_T>::eigenvalue_type(
                pkg.complex_eigenvalues[i].first,
                pkg.complex_eigenvalues[i].second));
    }
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        const scalar_real_type res_val = max_eigenpair_residual_value_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residual_norm_absolute = res_val;
    } else if (!pkg.residuals_abs.empty()) {
        result.residual_norm_absolute = pkg.residuals_abs[0];
    }
    const bool has_complex = pkg.has_complex;
    result.converged = pkg.converged && !has_complex && (result.eigenvalues.size() >= k);
    if (has_complex) {
        result.status = "complex_ritz_values";
        result.message = "complex Ritz values detected";
        if (result.failure_reason.empty())
            result.failure_reason = "complex Ritz values in requested subset";
    } else {
        set_eig_diagnostics_<_T,_Index>(result, eig_solver_method::arnoldi, sdim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// generalized_shift_invert_arnoldi_eigs_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> generalized_shift_invert_arnoldi_eigs_(const spmats<_T,_Index>& self,
                                                               const spmats<_T,_Index>& B,
                                                               const std::size_t k,
                                                               const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    const scalar_real_type sigma = options.shift;

    eig_solver_method actual_method = options.method;
    bool promoted_from_lanczos = false;
    if (options.method == eig_solver_method::shift_invert_lanczos) {
        actual_method = eig_solver_method::shift_invert_arnoldi;
        promoted_from_lanczos = true;
    } else if (options.use_shift && options.method == eig_solver_method::lanczos) {
        actual_method = eig_solver_method::shift_invert_arnoldi;
        promoted_from_lanczos = true;
    } else if (options.use_shift && options.method == eig_solver_method::arnoldi) {
        actual_method = eig_solver_method::shift_invert_arnoldi;
    }

    const std::size_t inner_max_iter = std::min(n, std::size_t(100));
    const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
    const std::size_t inner_restart = std::min(n, std::size_t(30));

    typedef vcp::tsparse::generalized_shift_invert_operator<spmats<_T,_Index> > GSIOperator;
    GSIOperator gsi_op(self, B, _T(sigma), inner_max_iter, inner_tol, inner_restart);

    if (!gsi_op.factorization_ok()) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = actual_method;
        result.used_method = promoted_from_lanczos
            ? "shift_invert_arnoldi(promoted_from_lanczos)"
            : eig_method_to_string_<_T,_Index>(actual_method);
        result.used_shift_invert = true;
        result.used_generalized_operator = true;
        result.used_dense_fallback = false;
        result.status = "factorization_failed";
        result.failure_reason = "ILU zero or near-zero pivot in (A - sigma*B)";
        result.message = result.failure_reason;
        result.factorization_diagnostics = gsi_op.factorization_diagnostics();
        result.factorization_zero_pivots = gsi_op.factorization_zero_pivots();
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    struct ApplyFn {
        GSIOperator* op;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
    } apply_fn = { &gsi_op };

    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const std::size_t max_restarts = options.max_iter + options.max_iter * k;

    const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
        vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
            n, k, sdim, max_restarts, options.tol, options.orthogonalization,
            true, true, options.random_seed, options.random_start,
            eig_target::largest_magnitude, scalar_real_type(0),
            options.compute_residual_history, apply_fn);

    std::vector<_T> eigenvalues = pkg.eigenvalues;
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0)) {
            eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
        }
    }

    eig_result<_T> result;
    result.requested_count = k;
    result.method = actual_method;
    result.used_method = promoted_from_lanczos
        ? "shift_invert_arnoldi(promoted_from_lanczos)"
        : eig_method_to_string_<_T,_Index>(actual_method);
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.linear_solves = gsi_op.linear_solves();
    result.inner_iterations = gsi_op.inner_iterations();
    result.inner_failure_count = gsi_op.inner_failure_count();
    result.inner_residual_norm = gsi_op.inner_residual_norm();
    result.factorization_diagnostics = gsi_op.factorization_diagnostics();
    result.factorization_zero_pivots = gsi_op.factorization_zero_pivots();
    result.used_subspace_dim = sdim;
    result.breakdown_reason = pkg.breakdown_reason;
    result.failure_reason = pkg.failure_reason;
    result.converged_count = pkg.converged_count;
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.used_shift_invert = true;
    result.used_generalized_operator = true;
    result.used_dense_fallback = false;
    result.eigenvalues = eigenvalues;
    result.eigenvectors = pkg.eigenvectors;

    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = generalized_eigenpair_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = generalized_eigenpair_relative_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        const scalar_real_type res_val = max_generalized_eigenpair_residual_value_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residual_norm_absolute = res_val;
    }

    for (std::size_t i = 0; i < pkg.complex_eigenvalues.size(); i++) {
        result.complex_eigenvalues.push_back(
            typename eig_result<_T>::eigenvalue_type(
                pkg.complex_eigenvalues[i].first,
                pkg.complex_eigenvalues[i].second));
    }

    const bool has_complex = pkg.has_complex;
    const bool inner_ok = (gsi_op.inner_failure_count() == 0);
    result.converged = pkg.converged && !has_complex && (result.eigenvalues.size() >= k) && inner_ok;

    if (!inner_ok) {
        result.status = "inner_solve_failed";
        result.failure_reason = "inner GMRES solve failed in generalized shift-invert";
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
    } else if (has_complex) {
        result.status = "complex_ritz_values";
        result.message = "complex Ritz values detected in generalized shift-invert";
        if (result.failure_reason.empty())
            result.failure_reason = "complex Ritz values in requested subset";
    } else {
        set_eig_diagnostics_<_T,_Index>(result, actual_method, sdim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
        if (promoted_from_lanczos)
            result.used_method = "shift_invert_arnoldi(promoted_from_lanczos)";
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// b_inner_lanczos_generalized_eigs_ (tag dispatch)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> b_inner_lanczos_generalized_eigs_impl_(const spmats<_T,_Index>&,
                                                               const spmats<_T,_Index>&,
                                                               const std::size_t k,
                                                               const eig_options<_T>& options,
                                                               spmats_b_inner_complex_tag_)
{
    return complex_generalized_eigs_unsupported_<_T,_Index>(k, options);
}

template <typename _T, typename _Index>
static eig_result<_T> b_inner_lanczos_generalized_eigs_impl_(const spmats<_T,_Index>& self,
                                                               const spmats<_T,_Index>& B,
                                                               const std::size_t k,
                                                               const eig_options<_T>& options,
                                                               spmats_b_inner_real_tag_)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());

    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::lanczos;
    result.used_method = "b_inner_lanczos";
    result.used_generalized_operator = true;
    result.used_shift_invert = false;
    result.used_dense_fallback = false;

    const scalar_real_type sym_tol = vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(10);
    if (!is_symmetric_value_<_T,_Index>(B, sym_tol)) {
        result.converged = false;
        result.status = "spd_check_failed";
        result.failure_reason = "B is not symmetric; cannot be SPD";
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const scalar_real_type norm_A_val = frobenius_norm_value_<_T,_Index>(self);
    const scalar_real_type norm_B_val = frobenius_norm_value_<_T,_Index>(B);

    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const std::size_t mv_budget = options.max_iter;

    typedef vcp::tsparse_b_inner_lanczos::b_lanczos_result<_T> BLResult;
    BLResult pkg = vcp::tsparse_b_inner_lanczos::b_inner_lanczos_eigs(
        self, B, n, k, sdim, mv_budget, options.tol,
        options.random_seed, options.random_start,
        options.target, norm_A_val, norm_B_val,
        options.compute_residual_history);

    if (pkg.spd_check_failed) {
        result.converged = false;
        result.status = "spd_check_failed";
        result.failure_reason = pkg.failure_reason.empty()
            ? "SPD check failed during B-inner Lanczos"
            : pkg.failure_reason;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.used_subspace_dim = sdim;
    result.breakdown_reason = pkg.breakdown_reason;
    result.failure_reason = pkg.failure_reason;
    result.converged_count = pkg.converged_count;
    result.eigenvalues = pkg.eigenvalues;
    result.eigenvectors = pkg.eigenvectors;
    result.residuals_absolute = pkg.residuals_abs;
    result.residuals_relative = pkg.residuals_rel;
    result.converged = pkg.converged;

    if (options.compute_residual_history) {
        result.residual_history_absolute = pkg.history_abs;
        result.residual_history_relative = pkg.history_rel;
    }

    populate_real_complex_eigenvalues_<_T,_Index>(result);
    if (result.converged) {
        result.status = "converged";
        result.message = pkg.breakdown_reason.empty() ? "converged" : pkg.breakdown_reason;
        result.failure_reason.clear();
    } else if (pkg.budget_exhausted) {
        result.status = "max_iter_exhausted";
        if (result.failure_reason.empty())
            result.failure_reason = "MV budget (max_iter) exhausted before convergence";
        result.message = result.failure_reason;
    } else {
        result.status = "not_converged";
        if (result.failure_reason.empty())
            result.failure_reason = "B-inner Lanczos did not converge";
        result.message = result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> b_inner_lanczos_generalized_eigs_(const spmats<_T,_Index>& self,
                                                          const spmats<_T,_Index>& B,
                                                          const std::size_t k,
                                                          const eig_options<_T>& options)
{
    typedef typename std::conditional<
        spmatrix_is_complex<_T>::value,
        spmats_b_inner_complex_tag_,
        spmats_b_inner_real_tag_>::type dispatch_tag;
    return b_inner_lanczos_generalized_eigs_impl_<_T,_Index>(self, B, k, options, dispatch_tag{});
}

// ---------------------------------------------------------------------------
// shift_invert_lanczos_eigs_with_prec_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class Preconditioner>
static eig_result<_T> shift_invert_lanczos_eigs_with_prec_(
    const spmats<_T,_Index>& self,
    const std::size_t k,
    const eig_options<_T>& options,
    const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
    const Preconditioner& M)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());

    spmats<_T,_Index> K_shifted = A;
    for (std::size_t i = 0; i < n; i++) {
        const _T cur = K_shifted.get(static_cast<_Index>(i), static_cast<_Index>(i));
        K_shifted.set(static_cast<_Index>(i), static_cast<_Index>(i), cur - _T(sigma));
    }
    K_shifted.finalize();
    spmats<_T,_Index> K_csr = K_shifted.as_csr();

    if (!vcp::tsparse_prec_traits::get_valid(M)) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = eig_solver_method::shift_invert_lanczos;
        result.used_method = "shift_invert_lanczos+user_preconditioner";
        result.used_shift_invert = true;
        result.used_dense_fallback = false;
        result.converged = false;
        result.status = "preconditioner_failed";
        result.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
        result.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
        result.failure_reason = "preconditioner is invalid: " + result.factorization_diagnostics;
        result.inner_failure_reason = result.failure_reason;
        result.inner_failure_count = 1;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const std::size_t inner_max = std::min(n, std::size_t(100));
    const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
    std::size_t linear_solve_count = 0;
    std::size_t inner_failure_count = 0;
    std::size_t inner_iteration_count = 0;
    scalar_real_type inner_residual_norm = scalar_real_type(0);
    std::string prec_failure_reason;

    struct SIApply {
        const spmats<_T,_Index>* K_ptr;
        const Preconditioner* prec_ptr;
        std::size_t imax;
        scalar_real_type itol;
        std::size_t* sc;
        std::size_t* fc;
        std::size_t* ic;
        scalar_real_type* mr;
        std::string* fm;

        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
            const std::size_t nn = x.size();
            struct AV {
                const spmats<_T,_Index>* K;
                void operator()(const std::vector<_T>& u, std::vector<_T>& v) const { v = K->mul_vec(u); }
            } av_fn = { K_ptr };
            bool pf = false;
            std::string pm;
            struct PV {
                const Preconditioner* M;
                bool* failed;
                std::string* msg;
                void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
                    if (*failed) { z = r; return; }
                    try {
                        M->apply(r, z);
                        if (z.size() != r.size()) { *failed = true; *msg = "preconditioner output size mismatch"; z = r; }
                    } catch (const vcp::error& e) {
                        *failed = true; *msg = std::string("preconditioner failed: ") + e.what(); z = r;
                    }
                }
            } pv_fn = { prec_ptr, &pf, &pm };
            typedef vcp::tsparse_factorization::gmres_result<_T, AV, PV> GR;
            GR gr = vcp::tsparse_factorization::gmres_solve<_T, AV, PV>(
                av_fn, pv_fn, x, imax / 10 + 1, itol, std::min(nn, imax));
            (*sc)++;
            (*ic) += gr.iterations;
            if (gr.residual_norm > *mr) *mr = gr.residual_norm;
            if (pf) { (*fc)++; if (fm->empty()) *fm = pm; } else if (!gr.converged) { (*fc)++; }
            y = gr.x;
        }
    } apply_si = { &K_csr, &M, inner_max, inner_tol,
                   &linear_solve_count, &inner_failure_count,
                   &inner_iteration_count, &inner_residual_norm,
                   &prec_failure_reason };

    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const std::size_t max_restarts_si = (sdim > 0) ? (options.max_iter / sdim + k + 1) : options.max_iter;

    typedef vcp::tsparse_lanczos::lanczos_result_package<_T, SIApply> LPkg;
    LPkg pkg = vcp::tsparse_lanczos::lanczos_eigs_standard<_T, SIApply>(
        n, k, sdim, max_restarts_si, options.tol,
        options.random_seed, options.random_start,
        eig_target::largest_magnitude, scalar_real_type(0),
        options.compute_residual_history, apply_si);

    for (std::size_t i = 0; i < pkg.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(pkg.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            pkg.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }

    eig_result<_T> result = lanczos_package_to_result_<_T,_Index>(pkg, self, k, eig_solver_method::shift_invert_lanczos);
    result.used_method = "shift_invert_lanczos+user_preconditioner";
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_failure_count = inner_failure_count;
    result.inner_residual_norm = inner_residual_norm;
    result.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
    result.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
    result.used_shift_invert = true;

    if (inner_failure_count != 0) {
        result.converged = false;
        result.eigenvalues.clear();
        result.eigenvectors.clear();
        result.status = "inner_solve_failed";
        result.failure_reason = prec_failure_reason.empty() ? "inner GMRES solve failed" : prec_failure_reason;
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// shift_invert_arnoldi_eigs_with_prec_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class Preconditioner>
static eig_result<_T> shift_invert_arnoldi_eigs_with_prec_(
    const spmats<_T,_Index>& self,
    const std::size_t k,
    const eig_options<_T>& options,
    const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
    const Preconditioner& M)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());

    spmats<_T,_Index> K_shifted = A;
    for (std::size_t i = 0; i < n; i++) {
        const _T cur = K_shifted.get(static_cast<_Index>(i), static_cast<_Index>(i));
        K_shifted.set(static_cast<_Index>(i), static_cast<_Index>(i), cur - _T(sigma));
    }
    K_shifted.finalize();
    spmats<_T,_Index> K_csr = K_shifted.as_csr();

    if (!vcp::tsparse_prec_traits::get_valid(M)) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = eig_solver_method::shift_invert_arnoldi;
        result.used_method = "shift_invert_arnoldi+user_preconditioner";
        result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
        result.used_shift_invert = true;
        result.used_dense_fallback = false;
        result.converged = false;
        result.status = "preconditioner_failed";
        result.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
        result.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
        result.failure_reason = "preconditioner is invalid: " + result.factorization_diagnostics;
        result.inner_failure_reason = result.failure_reason;
        result.inner_failure_count = 1;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const std::size_t inner_max = std::min(n, std::size_t(100));
    const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
    std::size_t linear_solve_count = 0;
    std::size_t inner_failure_count = 0;
    std::size_t inner_iteration_count = 0;
    scalar_real_type inner_residual_norm = scalar_real_type(0);
    std::string prec_failure_reason;

    struct SIApply {
        const spmats<_T,_Index>* K_ptr;
        const Preconditioner* prec_ptr;
        std::size_t imax;
        scalar_real_type itol;
        std::size_t* sc;
        std::size_t* fc;
        std::size_t* ic;
        scalar_real_type* mr;
        std::string* fm;

        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
            const std::size_t nn = x.size();
            struct AV {
                const spmats<_T,_Index>* K;
                void operator()(const std::vector<_T>& u, std::vector<_T>& v) const { v = K->mul_vec(u); }
            } av_fn = { K_ptr };
            bool pf = false;
            std::string pm;
            struct PV {
                const Preconditioner* M;
                bool* failed;
                std::string* msg;
                void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
                    if (*failed) { z = r; return; }
                    try {
                        M->apply(r, z);
                        if (z.size() != r.size()) { *failed = true; *msg = "preconditioner output size mismatch"; z = r; }
                    } catch (const vcp::error& e) {
                        *failed = true; *msg = std::string("preconditioner failed: ") + e.what(); z = r;
                    }
                }
            } pv_fn = { prec_ptr, &pf, &pm };
            typedef vcp::tsparse_factorization::gmres_result<_T, AV, PV> GR;
            GR gr = vcp::tsparse_factorization::gmres_solve<_T, AV, PV>(
                av_fn, pv_fn, x, imax / 10 + 1, itol, std::min(nn, imax));
            (*sc)++;
            (*ic) += gr.iterations;
            if (gr.residual_norm > *mr) *mr = gr.residual_norm;
            if (pf) { (*fc)++; if (fm->empty()) *fm = pm; } else if (!gr.converged) { (*fc)++; }
            y = gr.x;
        }
    } apply_si = { &K_csr, &M, inner_max, inner_tol,
                   &linear_solve_count, &inner_failure_count,
                   &inner_iteration_count, &inner_residual_norm,
                   &prec_failure_reason };

    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
        vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
            n, k, sdim, options.max_iter + options.max_iter * k,
            options.tol, options.orthogonalization, true, true,
            options.random_seed, options.random_start,
            eig_target::largest_magnitude, scalar_real_type(0),
            options.compute_residual_history, apply_si);

    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::shift_invert_arnoldi;
    result.used_method = "shift_invert_arnoldi+user_preconditioner";
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_failure_count = inner_failure_count;
    result.inner_residual_norm = inner_residual_norm;
    result.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
    result.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
    result.used_subspace_dim = sdim;
    result.breakdown_reason = pkg.breakdown_reason;
    result.failure_reason = pkg.failure_reason;
    result.converged_count = pkg.converged_count;
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.used_shift_invert = true;
    result.used_dense_fallback = false;
    result.eigenvalues = pkg.eigenvalues;
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(result.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            result.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }
    result.eigenvectors = pkg.eigenvectors;

    if (inner_failure_count != 0) {
        result.converged = false;
        result.eigenvalues.clear();
        result.eigenvectors.clear();
        result.status = "inner_solve_failed";
        result.failure_reason = prec_failure_reason.empty() ? "inner GMRES solve failed" : prec_failure_reason;
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    result.converged = pkg.converged && result.eigenvalues.size() >= k;
    if (result.converged) {
        result.status = "converged";
        result.message = "converged";
    } else {
        result.status = "not_converged";
        if (result.failure_reason.empty()) result.failure_reason = "Arnoldi shift-invert did not converge";
        result.message = result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// generalized_shift_invert_arnoldi_eigs_with_prec_
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class Preconditioner>
static eig_result<_T> generalized_shift_invert_arnoldi_eigs_with_prec_(
    const spmats<_T,_Index>& self,
    const spmats<_T,_Index>& B,
    const std::size_t k,
    const eig_options<_T>& options,
    const Preconditioner& M)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    const scalar_real_type sigma = options.shift;

    eig_solver_method actual_method = options.method;
    bool promoted_from_lanczos = false;
    if (options.method == eig_solver_method::shift_invert_lanczos
        || (options.use_shift && options.method == eig_solver_method::lanczos)) {
        actual_method = eig_solver_method::shift_invert_arnoldi;
        promoted_from_lanczos = true;
    } else if (options.use_shift && options.method == eig_solver_method::arnoldi) {
        actual_method = eig_solver_method::shift_invert_arnoldi;
    }

    const std::string base_name = promoted_from_lanczos
        ? "shift_invert_arnoldi(promoted_from_lanczos)"
        : "shift_invert_arnoldi";
    const std::string used_method_str = base_name + "+user_preconditioner";

    if (!vcp::tsparse_prec_traits::get_valid(M)) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = actual_method;
        result.used_method = used_method_str;
        result.used_shift_invert = true;
        result.used_generalized_operator = true;
        result.used_dense_fallback = false;
        result.converged = false;
        result.status = "preconditioner_failed";
        result.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
        result.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
        result.failure_reason = "preconditioner is invalid: " + result.factorization_diagnostics;
        result.inner_failure_reason = result.failure_reason;
        result.inner_failure_count = 1;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    spmats<_T,_Index> K_csr = vcp::tsparse::subtract_scaled_sparse(self, B, _T(sigma)).as_csr();
    spmats<_T,_Index> B_csr = B.as_csr();

    const std::size_t inner_max = std::min(n, std::size_t(100));
    const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
    const std::size_t inner_restart = std::min(n, std::size_t(30));
    std::size_t linear_solve_count = 0;
    std::size_t inner_failure_count = 0;
    std::size_t inner_iteration_count = 0;
    scalar_real_type inner_residual_norm = scalar_real_type(0);
    std::string prec_failure_reason;

    struct GSIApply {
        const spmats<_T,_Index>* K_ptr;
        const spmats<_T,_Index>* B_ptr;
        const Preconditioner* prec_ptr;
        std::size_t imax;
        scalar_real_type itol;
        std::size_t irestart;
        std::size_t* sc;
        std::size_t* fc;
        std::size_t* ic;
        scalar_real_type* mr;
        std::string* fm;

        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
            const std::size_t nn = x.size();
            std::vector<_T> rhs = B_ptr->mul_vec(x);
            struct AV {
                const spmats<_T,_Index>* K;
                void operator()(const std::vector<_T>& u, std::vector<_T>& v) const { v = K->mul_vec(u); }
            } av_fn = { K_ptr };
            bool pf = false;
            std::string pm;
            struct PV {
                const Preconditioner* M;
                bool* failed;
                std::string* msg;
                void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
                    if (*failed) { z = r; return; }
                    try {
                        M->apply(r, z);
                        if (z.size() != r.size()) { *failed = true; *msg = "preconditioner output size mismatch"; z = r; }
                    } catch (const vcp::error& e) {
                        *failed = true; *msg = std::string("preconditioner failed: ") + e.what(); z = r;
                    }
                }
            } pv_fn = { prec_ptr, &pf, &pm };
            const std::size_t act_restart = (irestart == 0) ? std::min(nn, std::size_t(30)) : irestart;
            const std::size_t act_max = (imax == 0) ? (act_restart + 1) : imax;
            typedef vcp::tsparse_factorization::gmres_result<_T, AV, PV> GR;
            GR gr = vcp::tsparse_factorization::gmres_solve<_T, AV, PV>(av_fn, pv_fn, rhs, act_max, itol, act_restart);
            (*sc)++;
            (*ic) += gr.iterations;
            if (gr.residual_norm > *mr) *mr = gr.residual_norm;
            if (pf) { (*fc)++; if (fm->empty()) *fm = pm; } else if (!gr.converged) { (*fc)++; }
            y = gr.x;
        }
    } apply_gsi = { &K_csr, &B_csr, &M, inner_max, inner_tol, inner_restart,
                    &linear_solve_count, &inner_failure_count,
                    &inner_iteration_count, &inner_residual_norm,
                    &prec_failure_reason };

    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    const std::size_t max_restarts = options.max_iter + options.max_iter * k;

    const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
        vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
            n, k, sdim, max_restarts, options.tol, options.orthogonalization,
            true, true, options.random_seed, options.random_start,
            eig_target::largest_magnitude, scalar_real_type(0),
            options.compute_residual_history, apply_gsi);

    std::vector<_T> eigenvalues = pkg.eigenvalues;
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }

    eig_result<_T> result;
    result.requested_count = k;
    result.method = actual_method;
    result.used_method = used_method_str;
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_failure_count = inner_failure_count;
    result.inner_residual_norm = inner_residual_norm;
    result.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
    result.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
    result.used_subspace_dim = sdim;
    result.breakdown_reason = pkg.breakdown_reason;
    result.failure_reason = pkg.failure_reason;
    result.converged_count = pkg.converged_count;
    result.residual_history_absolute = pkg.history_abs;
    result.residual_history_relative = pkg.history_rel;
    result.used_shift_invert = true;
    result.used_generalized_operator = true;
    result.used_dense_fallback = false;
    result.eigenvalues = eigenvalues;
    result.eigenvectors = pkg.eigenvectors;

    if (inner_failure_count != 0) {
        result.converged = false;
        result.eigenvalues.clear();
        result.eigenvectors.clear();
        result.status = "inner_solve_failed";
        result.failure_reason = prec_failure_reason.empty()
            ? "inner GMRES solve failed in generalized shift-invert with user preconditioner"
            : prec_failure_reason;
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = generalized_eigenpair_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = generalized_eigenpair_relative_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        const scalar_real_type res_val = max_generalized_eigenpair_residual_value_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residual_norm_absolute = res_val;
    }

    for (std::size_t i = 0; i < pkg.complex_eigenvalues.size(); i++) {
        result.complex_eigenvalues.push_back(
            typename eig_result<_T>::eigenvalue_type(
                pkg.complex_eigenvalues[i].first,
                pkg.complex_eigenvalues[i].second));
    }

    const bool has_complex = pkg.has_complex;
    const bool inner_ok = (inner_failure_count == 0);
    result.converged = pkg.converged && !has_complex && (result.eigenvalues.size() >= k) && inner_ok;

    if (has_complex) {
        result.status = "complex_ritz_values";
        result.message = "complex Ritz values detected";
        if (result.failure_reason.empty())
            result.failure_reason = "complex Ritz values in requested subset";
    } else {
        set_eig_diagnostics_<_T,_Index>(result, actual_method, sdim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
        result.used_method = used_method_str;
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ===========================================================================
// PUBLIC POLICY METHODS  (non-member free functions operating on spmats)
// These mirror spmatrix::eigs_with_info / spmatrix::eigs_with_info(B,...).
// ===========================================================================

// ---------------------------------------------------------------------------
// policy_eigs_with_info  (standard, no preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
eig_result<_T> policy_eigs_with_info(const spmats<_T,_Index>& A,
                                      const std::size_t k,
                                      const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (spmatrix_is_complex<_T>::value)
        return complex_standard_eigs_with_info_<_T,_Index>(A, k, options);
    eig_options<_T> active = resolve_eigs_options_<_T,_Index>(A, options);
    validate_eigs_input_<_T,_Index>(A, k, "spmats::eigs", active.method);
    if (active.max_iter == 0 || active.tol <= scalar_real_type(0)) {
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs: invalid iteration option");
    }
    if (active.method == eig_solver_method::lanczos
     || active.method == eig_solver_method::shift_invert_lanczos)
        return lanczos_eigs_new_<_T,_Index>(A, k, active);
    if (active.method == eig_solver_method::arnoldi
     || active.method == eig_solver_method::shift_invert_arnoldi)
        return arnoldi_eigs_new_<_T,_Index>(A, k, active);
    if (active.method != eig_solver_method::dense_fallback_explicit) {
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs: unknown eigensolver method");
    }
    // dense fallback
    if (active.max_iter == 0 || active.tol <= scalar_real_type(0))
        vcp::throw_error<vcp::invalid_argument>("spmats::eig: invalid iteration option");
    if (A.rowsize() != A.columnsize())
        vcp::throw_error<vcp::dimension_error>("spmats::eig: matrix must be square");
    check_dense_allowed_<_T,_Index>(A, active, "spmats::eig");
    std::vector<std::vector<_T> > dense = to_dense_impl_<_T,_Index>(A);
    eig_result<_T> result = dense_eig_<_T,_Index>(dense, active);
    result.used_dense_fallback = true;
    select_eigenpairs_<_T,_Index>(result, k, active.target, active.shift);
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// policy_eigs_with_info  (standard, with preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class Preconditioner>
eig_result<_T> policy_eigs_with_info(const spmats<_T,_Index>& A,
                                      const std::size_t k,
                                      const eig_options<_T>& options,
                                      const Preconditioner& M)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (spmatrix_is_complex<_T>::value)
        return complex_preconditioned_eigs_unsupported_<_T,_Index>(k, options);
    eig_options<_T> active = resolve_eigs_options_<_T,_Index>(A, options);
    validate_eigs_input_<_T,_Index>(A, k, "spmats::eigs_with_info(with preconditioner)", active.method);
    if (active.max_iter == 0 || active.tol <= scalar_real_type(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "spmats::eigs_with_info(with preconditioner): invalid iteration option");
    }
    const scalar_real_type sigma = active.shift;
    const bool is_shift_invert =
        (active.method == eig_solver_method::shift_invert_lanczos)
        || (active.method == eig_solver_method::shift_invert_arnoldi)
        || (active.use_shift && active.method == eig_solver_method::lanczos)
        || (active.use_shift && active.method == eig_solver_method::arnoldi);

    if (!is_shift_invert) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_preconditioner_for_method";
        result.failure_reason =
            "preconditioner overload is supported only for shift-invert methods"
            " (shift_invert_lanczos, shift_invert_arnoldi, or use_shift=true);"
            " for non-shift Lanczos/Arnoldi, pass no preconditioner";
        result.message = result.failure_reason;
        result.method = active.method;
        result.used_method = eig_method_to_string_<_T,_Index>(active.method);
        result.used_shift_invert = false;
        result.used_dense_fallback = false;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const bool use_lanczos =
        (active.method == eig_solver_method::shift_invert_lanczos)
        || (active.use_shift && active.method == eig_solver_method::lanczos);

    eig_result<_T> result;
    if (use_lanczos)
        result = shift_invert_lanczos_eigs_with_prec_<_T,_Index>(A, k, active, sigma, M);
    else
        result = shift_invert_arnoldi_eigs_with_prec_<_T,_Index>(A, k, active, sigma, M);
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// policy_generalized_eigs_with_info  (no preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
eig_result<_T> policy_generalized_eigs_with_info(const spmats<_T,_Index>& A,
                                                   const spmats<_T,_Index>& B,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (spmatrix_is_complex<_T>::value)
        return complex_generalized_eigs_unsupported_<_T,_Index>(k, options);
    validate_generalized_eig_input_<_T,_Index>(A, B, "spmats::eigs(A,B)");

    const std::size_t n = static_cast<std::size_t>(A.rowsize());

    if (k == 0) {
        eig_result<_T> result;
        result.requested_count = 0;
        result.converged = true;
        result.status = "converged";
        result.message = "converged";
        result.used_generalized_operator = true;
        result.used_dense_fallback = false;
        result.used_shift_invert = false;
        result.method = options.method;
        if (options.method == eig_solver_method::lanczos
            && options.structure == matrix_structure_hint::symmetric) {
            result.used_method = "b_inner_lanczos";
        } else {
            result.used_method = eig_method_to_string_<_T,_Index>(options.method);
        }
        set_result_counts_<_T,_Index>(result, 0);
        return result;
    }

    if (n == 0) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "not_converged";
        result.failure_reason = "matrix is empty (n == 0) but k > 0";
        result.message = result.failure_reason;
        result.used_generalized_operator = true;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const std::size_t k_eff = (k < n) ? k : n;

    if (options.tol <= scalar_real_type(0)) {
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs(A,B): tol must be positive");
    }

    if (is_diagonal_matrix_<_T,_Index>(A) && is_diagonal_matrix_<_T,_Index>(B)) {
        if (options.max_iter == 0) {
            vcp::throw_error<vcp::invalid_argument>("spmats::eigs(A,B): max_iter must be positive");
        }
        eig_result<_T> result = generalized_diagonal_eigs_<_T,_Index>(A, B, k_eff, options);
        result.requested_count = k;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const bool wants_shift_invert =
        (options.method == eig_solver_method::shift_invert_lanczos)
        || (options.method == eig_solver_method::shift_invert_arnoldi)
        || (options.use_shift && options.method == eig_solver_method::lanczos)
        || (options.use_shift && options.method == eig_solver_method::arnoldi);
    if (wants_shift_invert) {
        if (options.max_iter == 0) {
            vcp::throw_error<vcp::invalid_argument>("spmats::eigs(A,B): max_iter must be positive");
        }
        eig_result<_T> result = generalized_shift_invert_arnoldi_eigs_<_T,_Index>(A, B, k_eff, options);
        result.requested_count = k;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    if (options.method == eig_solver_method::lanczos
        && options.structure == matrix_structure_hint::symmetric) {
        eig_result<_T> result = b_inner_lanczos_generalized_eigs_<_T,_Index>(A, B, k_eff, options);
        result.requested_count = k;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    if (options.max_iter == 0) {
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs(A,B): max_iter must be positive");
    }
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported_generalized_non_diagonal";
    result.used_generalized_operator = true;
    result.used_dense_fallback = false;
    result.method = options.method;
    result.used_method = eig_method_to_string_<_T,_Index>(options.method);
    result.failure_reason = "non-diagonal generalized sparse eigs without shift-invert is not supported;"
        " use shift_invert_arnoldi, shift_invert_lanczos, or set use_shift=true with symmetric structure";
    result.message = result.failure_reason;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// policy_generalized_eigs_with_info  (with preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class Preconditioner>
eig_result<_T> policy_generalized_eigs_with_info(const spmats<_T,_Index>& A,
                                                   const spmats<_T,_Index>& B,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options,
                                                   const Preconditioner& M)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (spmatrix_is_complex<_T>::value)
        return complex_generalized_eigs_unsupported_<_T,_Index>(k, options);
    validate_generalized_eig_input_<_T,_Index>(A, B, "spmats::eigs_with_info(A,B,with preconditioner)");

    const std::size_t n = static_cast<std::size_t>(A.rowsize());

    if (k == 0) {
        eig_result<_T> result;
        result.requested_count = 0;
        result.converged = true;
        result.status = "converged";
        result.message = "converged";
        result.used_generalized_operator = true;
        result.used_dense_fallback = false;
        result.used_shift_invert = false;
        result.method = options.method;
        result.used_method = eig_method_to_string_<_T,_Index>(options.method);
        set_result_counts_<_T,_Index>(result, 0);
        return result;
    }

    if (n == 0) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "not_converged";
        result.failure_reason = "matrix is empty (n == 0) but k > 0";
        result.message = result.failure_reason;
        result.used_generalized_operator = true;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    const std::size_t k_eff = (k < n) ? k : n;

    if (options.tol <= scalar_real_type(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "spmats::eigs_with_info(A,B,with preconditioner): tol must be positive");
    }

    const bool wants_shift_invert =
        (options.method == eig_solver_method::shift_invert_lanczos)
        || (options.method == eig_solver_method::shift_invert_arnoldi)
        || (options.use_shift && options.method == eig_solver_method::lanczos)
        || (options.use_shift && options.method == eig_solver_method::arnoldi);

    if (!wants_shift_invert) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_preconditioner_for_method";
        result.failure_reason =
            "preconditioner overload for generalized eigs is supported only"
            " for shift-invert methods; B-inner Lanczos and non-shift paths"
            " do not support external preconditioners in Phase 6";
        result.message = result.failure_reason;
        result.method = options.method;
        result.used_method = eig_method_to_string_<_T,_Index>(options.method);
        result.used_generalized_operator = true;
        result.used_shift_invert = false;
        result.used_dense_fallback = false;
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }

    if (options.max_iter == 0) {
        vcp::throw_error<vcp::invalid_argument>(
            "spmats::eigs_with_info(A,B,with preconditioner): max_iter must be positive");
    }

    eig_result<_T> result =
        generalized_shift_invert_arnoldi_eigs_with_prec_<_T,_Index>(A, B, k_eff, options, M);
    result.requested_count = k;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ===========================================================================
// spmats<_T,_Index> member function definitions for policy_eigs_with_info
// and policy_generalized_eigs_with_info.
// These delegate to the free function implementations above.
// ===========================================================================

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_eigs_with_info(
    const spmats<_T,_Index>& A,
    std::size_t k,
    const eig_options<_T>& opt) const
{
    return vcp::policy_eigs_with_info<_T,_Index>(A, k, opt);
}

template <typename _T, typename _Index>
template <class Prec>
eig_result<_T> spmats<_T,_Index>::policy_eigs_with_info(
    const spmats<_T,_Index>& A,
    std::size_t k,
    const eig_options<_T>& opt,
    const Prec& M) const
{
    return vcp::policy_eigs_with_info<_T,_Index,Prec>(A, k, opt, M);
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eigs_with_info(
    const spmats<_T,_Index>& A,
    const spmats<_T,_Index>& B,
    std::size_t k,
    const eig_options<_T>& opt) const
{
    return vcp::policy_generalized_eigs_with_info<_T,_Index>(A, B, k, opt);
}

template <typename _T, typename _Index>
template <class Prec>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eigs_with_info(
    const spmats<_T,_Index>& A,
    const spmats<_T,_Index>& B,
    std::size_t k,
    const eig_options<_T>& opt,
    const Prec& M) const
{
    return vcp::policy_generalized_eigs_with_info<_T,_Index,Prec>(A, B, k, opt, M);
}

// ===========================================================================
// Strict policy methods: policy_eig, policy_eigs, policy_generalized_eig,
// policy_generalized_eigs (and preconditioner overloads).
// All convergence/count checking lives here, not in spmatrix.hpp.
// ===========================================================================

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_eig(
    const spmats<_T,_Index>& A,
    const eig_options<_T>& opt) const
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    // default policy: full dense eig is only for real scalar types
    if (spmatrix_is_complex<_T>::value)
        vcp::throw_error<vcp::domain_error>("spmats::policy_eig: complex scalar is not supported for full dense eig; use eigs_with_info for complex Hermitian matrices");
    if (A.rowsize() != A.columnsize())
        vcp::throw_error<vcp::dimension_error>("spmats::policy_eig: matrix must be square");
    if (opt.max_iter == 0 || opt.tol <= scalar_real_type(0))
        vcp::throw_error<vcp::invalid_argument>("spmats::policy_eig: invalid options (max_iter or tol)");
    if (opt.method != eig_solver_method::dense_fallback_explicit)
        vcp::throw_error<vcp::invalid_argument>("spmats::policy_eig: full dense eig requires dense_fallback_explicit method");
    eig_result<_T> result = policy_eigs_with_info(A, static_cast<std::size_t>(A.rowsize()), opt);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_eig: eigensolver did not converge");
    return result;
}

template <typename _T, typename _Index>
std::vector<_T> spmats<_T,_Index>::policy_eigs(
    const spmats<_T,_Index>& A, std::size_t k,
    const eig_options<_T>& opt) const
{
    eig_result<_T> result = policy_eigs_with_info(A, k, opt);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_eigs: eigensolver did not converge");
    const std::size_t n = static_cast<std::size_t>(A.rowsize());
    const std::size_t k_eff = (k < n) ? k : n;
    if (result.returned_real_count < k_eff)
        vcp::throw_error<vcp::state_error>("spmats::policy_eigs: insufficient eigenvalues returned");
    return result.eigenvalues;
}

template <typename _T, typename _Index>
template <class Prec>
std::vector<_T> spmats<_T,_Index>::policy_eigs(
    const spmats<_T,_Index>& A, std::size_t k,
    const eig_options<_T>& opt, const Prec& M) const
{
    eig_result<_T> result = policy_eigs_with_info(A, k, opt, M);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_eigs(with preconditioner): eigensolver did not converge");
    if (result.returned_real_count < k)
        vcp::throw_error<vcp::state_error>("spmats::policy_eigs(with preconditioner): insufficient eigenvalues returned");
    return result.eigenvalues;
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eig(
    const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt) const
{
    eig_result<_T> result = policy_generalized_eigs_with_info(A, B, k, opt);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_generalized_eig: eigensolver did not converge");
    return result;
}

template <typename _T, typename _Index>
std::vector<_T> spmats<_T,_Index>::policy_generalized_eigs(
    const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt) const
{
    eig_result<_T> result = policy_generalized_eigs_with_info(A, B, k, opt);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_generalized_eigs: eigensolver did not converge");
    const std::size_t n = static_cast<std::size_t>(A.rowsize());
    const std::size_t k_eff = (k < n) ? k : n;
    if (result.returned_real_count < k_eff)
        vcp::throw_error<vcp::state_error>("spmats::policy_generalized_eigs: insufficient eigenvalues returned");
    return result.eigenvalues;
}

template <typename _T, typename _Index>
template <class Prec>
std::vector<_T> spmats<_T,_Index>::policy_generalized_eigs(
    const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt, const Prec& M) const
{
    eig_result<_T> result = policy_generalized_eigs_with_info(A, B, k, opt, M);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_generalized_eigs(with preconditioner): eigensolver did not converge");
    const std::size_t n = static_cast<std::size_t>(A.rowsize());
    const std::size_t k_eff = (k < n) ? k : n;
    if (result.returned_real_count < k_eff)
        vcp::throw_error<vcp::state_error>("spmats::policy_generalized_eigs(with preconditioner): insufficient eigenvalues returned");
    return result.eigenvalues;
}

} // namespace vcp

#endif // VCP_SPMATS_EIGS_HPP
