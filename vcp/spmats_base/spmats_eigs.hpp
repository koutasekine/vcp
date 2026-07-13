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
#include <vcp/tsparse/tsparse_honest_termination.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>
#include <vcp/tsparse/tsparse_krylov_schur.hpp>   // EIG-3: arnoldi/si successor
#include <vcp/tsparse/tsparse_thick_restart_lanczos.hpp>   // EIG-4: TRL promotion (T-1/T-2)
#include <vcp/tsparse/tsparse_factorization.hpp>
#include <vcp/tsparse/tsparse_preconditioner.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_generalized_shift_invert.hpp>
#include <vcp/tsparse/tsparse_b_inner_lanczos.hpp>
#include <vcp/tsparse/tsparse_hermitian_lanczos.hpp>
#include <vcp/tsparse/tsparse_dense_fallback.hpp>
#include <vcp/tsparse/tsparse_dense_schur_driver.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_solvers.hpp>
// NOTE: spmats.hpp includes this file after spmats<_T,_Index> is defined.
// Do NOT #include <vcp/spmats.hpp> here to avoid circular dependency.

// EIG-4: forward declarations for the TRL driver symbols.  When a TU includes
// tsparse_thick_restart_lanczos.hpp FIRST, its include chain (tsparse_restart
// -> spmatrix -> spmats -> this header) re-enters before the driver namespace
// is populated (the include above is then a guard no-op).  The dispatch
// wrappers below only need these declarations to parse; instantiation happens
// after the full driver definition is available in the TU.
namespace vcp {
namespace tsparse_experimental {
template <typename T> struct trl_is_real_floating;
template <class Apply, class T>
vcp::eig_result<T> thick_restart_lanczos_eigs(
    const Apply& apply, std::size_t n, std::size_t k,
    const vcp::eig_options<T>& options);
} // namespace tsparse_experimental
} // namespace vcp

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
// 1b. is_certainly_symmetric_   (EIG-4 T-2; B-28 の唯一の実装)
//
// auto_select の対称判定: 決定的・O(nnz log nnz)・T の素の == のみを使う。
// GT1 P1 の比較規約により、kv::interval の == は certainly 等号(両端点一致の
// 退化区間同士でのみ真)なので、認証不能な対称性は自動的に false(= KS 側)に
// 落ちる。「たぶん対称」で TRL に送ることはできない。片側欠落エントリは
// 明示格納値 T(0) との certainly 等号で判定する(明示ゼロの片側パターンは
// 対称と認証できる)。tol ベースの is_symmetric_value_(明示 lanczos の
// 入口検査)とは役割が異なり、置き換えない(B-26)。
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static bool is_certainly_symmetric_(const spmats<_T,_Index>& A)
{
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
            // certainly ==(interval では退化区間の一致のみ真)。否定は
            // 「等しいと保証できない」= 非対称側(KS)へ倒す(B-28)。
            if (!(val[static_cast<std::size_t>(p)] == mirrored)) return false;
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
    case eig_solver_method::thick_restart_lanczos:    return "thick_restart_lanczos";
    case eig_solver_method::krylov_schur:             return "krylov_schur";
    case eig_solver_method::auto_select:              return "auto_select";
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
                > vcp::tsparse_scalar::epsilon<scalar_real_type>()) {
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
// 10. populate_real_complex_eigenvalues_  (declared before select_eigenpairs_aligned_ which calls it)
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
// 9. select_eigenpairs_aligned_   (EIG-4 T-2 導入、EIG-4.1 で唯一の選択実装)
//
// 値・ベクトル・残差を同一置換で選択する(構築した選択済み配列を無条件に
// 採用し、非対応の配列はクリアする)。旧 select_eigenpairs_ は末尾の swap
// ガードが「選択前サイズ == 選択後サイズ」を要求し、k < n の選択で
// eigenvectors / residuals が未選択のまま残る嘘クラス欠陥(D-16)を持って
// いたため、EIG-4.1 で全呼び出し面(実 dense・complex dense の 2 箇所)を
// 本整列版へ統一の上、旧実装を削除した。
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void select_eigenpairs_aligned_(eig_result<_T>& result,
                                       const std::size_t k,
                                       const eig_target target,
                                       const typename vcp::tsparse_scalar::real_type<_T>::type& shift)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::vector<std::size_t> order =
        vcp::tsparse_eigen_selection::select_real_eigenpairs(result.eigenvalues, k, target, shift);
    const bool have_vectors = (result.eigenvectors.size()       == result.eigenvalues.size());
    const bool have_resabs  = (result.residuals_absolute.size() == result.eigenvalues.size());
    const bool have_resrel  = (result.residuals_relative.size() == result.eigenvalues.size());
    std::vector<_T> values;
    std::vector<std::vector<_T> > vectors;
    std::vector<scalar_real_type> residuals_abs;
    std::vector<scalar_real_type> residuals_rel;
    for (std::size_t i = 0; i < order.size(); i++) {
        const std::size_t j = order[i];
        values.push_back(result.eigenvalues[j]);
        if (have_vectors) vectors.push_back(result.eigenvectors[j]);
        if (have_resabs)  residuals_abs.push_back(result.residuals_absolute[j]);
        if (have_resrel)  residuals_rel.push_back(result.residuals_relative[j]);
    }
    result.eigenvalues.swap(values);
    if (have_vectors) result.eigenvectors.swap(vectors);
    else              result.eigenvectors.clear();
    if (have_resabs) result.residuals_absolute.swap(residuals_abs);
    else             result.residuals_absolute.clear();
    if (have_resrel) result.residuals_relative.swap(residuals_rel);
    else             result.residuals_relative.clear();
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
// 13b. matrix_inf_norm_value_   (EIG-4 T-4: 改訂 C-1 scale 用の決定的 ‖A‖∞)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename vcp::tsparse_scalar::real_type<_T>::type
matrix_inf_norm_value_(const spmats<_T,_Index>& A)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> C = A.as_csr();
    const std::vector<_Index>& outer = C.outer_index();
    const std::vector<_T>& val = C.values();
    scalar_real_type best(0);
    for (_Index i = 0; i < C.rowsize(); i++) {
        scalar_real_type s(0);
        for (_Index p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++)
            s += vcp::tsparse_scalar::abs_value(val[static_cast<std::size_t>(p)]);
        if (s > best) best = s;
    }
    return best;
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
        + vcp::tsparse_scalar::epsilon<scalar_real_type>();
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
    // SLU-GT1 D5: defensive branch only -- every caller guards
    // !eigenvectors.empty() and derives `converged` independently of this
    // value, so scalar_real_type(0) is a placeholder, not a sentinel.
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
        return scalar_real_type(0);
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
// EIG-2 Phase 3 (D-13/D-15): the nonsymmetric REAL branch is the new
// balancing + Hessenberg + real Schur core driver.  The complex-scalar
// nonsymmetric branch keeps the legacy qr_eig_dense (non-Hermitian complex
// dense is outside the EIG-2 scope; the real Schur core is real-T only).
// The symmetric branch (jacobi_eig_dense) is untouched (B-21).
template <typename _T, typename _Index>
static eig_result<_T> dense_nonsymmetric_eig_(const std::vector<std::vector<_T> >& dense,
                                              const eig_options<_T>& options,
                                              std::true_type /* is_complex */)
{
    return convert_dense_result_<_T,_Index>(
        vcp::tsparse_dense_linalg::qr_eig_dense(dense, options.max_iter, options.tol),
        eig_solver_method::dense_fallback_explicit);
}

template <typename _T, typename _Index>
static eig_result<_T> dense_nonsymmetric_eig_(const std::vector<std::vector<_T> >& dense,
                                              const eig_options<_T>& options,
                                              std::false_type /* is_complex */)
{
    std::string reason;
    eig_result<_T> result = convert_dense_result_<_T,_Index>(
        vcp::tsparse_dense_schur::real_schur_eig_dense(dense, options.tol, reason),
        eig_solver_method::dense_fallback_explicit);
    if (!result.converged && !reason.empty()) {
        result.status         = "not_converged";
        result.failure_reason = reason;
        result.message        = reason;
    }
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> dense_eig_(std::vector<std::vector<_T> > dense, const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    eig_result<_T> result = vcp::tsparse_dense_linalg::is_dense_symmetric(dense, options.tol * scalar_real_type(10))
        ? convert_dense_result_<_T,_Index>(
              vcp::tsparse_dense_linalg::jacobi_eig_dense(dense, options.max_iter, options.tol),
              eig_solver_method::dense_fallback_explicit)
        : dense_nonsymmetric_eig_<_T,_Index>(dense, options,
              typename std::integral_constant<bool, spmatrix_is_complex<_T>::value>::type());
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
        + vcp::tsparse_scalar::epsilon<scalar_real_type>();
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
    // SLU-GT1 D5: defensive branch only -- see max_eigenpair_residual_value_.
    if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
        return scalar_real_type(0);
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
    const eig_solver_method method,
    const typename vcp::tsparse_scalar::real_type<_T>::type& tol)
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
    bool exact_residual_rejected = false;
    if (!result.eigenvectors.empty()) {
        spmats<_T,_Index> A = self.as_csr();
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        // EIG-1 F-4/F-5 (EIG-0 C-1): 直前で再計算した「元の A に対する λ 空間の
        // 厳密残差」を converged 判定へ反映する(降格のみ — false を true には
        // しない)。shift-invert lanczos では内側反復は μ 空間だが、受理は
        // この λ 空間残差で判定する(G-2.1 承認案 ①〜④)。
        // EIG-4 T-4 (D4-4 / R-1): 受理は「既存式 ∨ res_abs ≤ tol·scale」、
        // scale = max(1+|θ|, ‖A‖∞ 厳密値)(共有ヘルパ。緩和方向のみ — B-29)
        std::vector<scalar_real_type> theta_abs_c1;
        theta_abs_c1.reserve(result.eigenvalues.size());
        for (std::size_t i2 = 0; i2 < result.eigenvalues.size(); i2++)
            theta_abs_c1.push_back(vcp::tsparse_scalar::abs_value(
                vcp::tsparse_scalar::real_part(result.eigenvalues[i2])));
        if (result.converged &&
            !vcp::tsparse::residual_acceptance_check_scaled_(
                result.residuals_absolute, result.residuals_relative, tol,
                theta_abs_c1, matrix_inf_norm_value_<_T,_Index>(A))) {
            exact_residual_rejected = true;
            result.converged = false;
            result.failure_reason =
                "end-of-run exact residual (lambda space) failed acceptance (C-1)";
        }
    }
    set_eig_diagnostics_<_T,_Index>(result, method, result.returned_count, result.matrix_vector_products,
        result.breakdown_reason, result.failure_reason);
    if (exact_residual_rejected) result.status = "residual_check_failed";
    result.used_shift_invert = (method == eig_solver_method::shift_invert_lanczos
                             || method == eig_solver_method::shift_invert_arnoldi);
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// demote_arnoldi_core_converged_claim_   (EIG-1 F-5; G-2.1 承認 案 (a))
//
// EIG-0 C-2(返却集合の内側に未収束候補が certainly 存在しないこと)は
// converged=true の必要条件だが、arnoldi コアは active 候補集合を結果
// パッケージとして輸出しておらず(かつ B-13 により編集禁止・EIG-3 で置換予定)、
// arnoldi コア shift-invert 経路の back-half では C-2 を評価できない。
// したがって同経路の converged 主張は正直な not_converged へ降格する。
//
// 判定順序(a-2): まず λ 空間の厳密残差受理(C-1)を評価し、
//   C-1 不合格 → status="residual_check_failed"
//   C-1 合格かつ C-2 証拠なし → (a-1) の明示文言で降格
// 値と λ 空間厳密残差は診断として結果に残る(B-15)。enum 追加なし(文字列のみ)。
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void demote_arnoldi_core_converged_claim_(
    eig_result<_T>& result,
    const typename vcp::tsparse_scalar::real_type<_T>::type& tol)
{
    const bool would_claim = result.converged;
    result.converged = false;
    if (!would_claim) return;   // 元々 not_converged: 触らない
    if (!vcp::tsparse::residual_acceptance_check_(
            result.residuals_absolute, result.residuals_relative, tol)) {
        result.status = "residual_check_failed";
        result.failure_reason =
            "end-of-run exact residual (lambda space) failed acceptance (C-1)";
    } else {
        result.status = "not_converged";
        result.failure_reason =
            "C-2 evidence unavailable on arnoldi-core shift-invert path (pending EIG-3)";
    }
    result.message = result.failure_reason;
}

// ---------------------------------------------------------------------------
// EIG-3 T-3: shift-invert back-halves on the rebuilt Krylov-Schur core.
//
// The KS driver (tsparse_krylov_schur.hpp) enforces C-1 (exact mu-space
// residuals), C-2 (shared honest_termination_check_complex_), D3-2 (complex
// pairs never returned as converged) and freshness internally, and exports
// the C-2 evidence (D3-4 / EIG-1 a-5).  The back-half additionally
// RE-EVALUATES the shared C-2 check on the exported evidence (design D3-4:
// "back-half passes the evidence to the shared check"), then the lambda-space
// C-1 acceptance below remains the final outer gate.  This lifts the EIG-1
// arnoldi-core demotion (demote_arnoldi_core_converged_claim_, now unused).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static void lambda_c1_acceptance_gate_(eig_result<_T>& result,
    const typename vcp::tsparse_scalar::real_type<_T>::type& tol)
{
    if (!result.converged) return;
    if (!vcp::tsparse::residual_acceptance_check_(
            result.residuals_absolute, result.residuals_relative, tol)) {
        result.converged = false;
        result.status = "residual_check_failed";
        result.failure_reason =
            "end-of-run exact residual (lambda space) failed acceptance (C-1)";
        result.message = result.failure_reason;
    }
}

// EIG-4 T-4 (D4-4 / R-1): 標準問題用の改訂 scale 版(A x = λ x の si 経路)。
// 受理は共有ヘルパ residual_acceptance_check_scaled_(既存式 ∨ tol·scale、
// scale = max(1+|θ|, ‖A‖∞))。
// (旧注記「一般化経路は従来ゲート維持」は EIG-8 T-3/D8-3 で共有形へ置換 —
//  下の 4 引数 overload。)
template <typename _T, typename _Index>
static void lambda_c1_acceptance_gate_(eig_result<_T>& result,
    const typename vcp::tsparse_scalar::real_type<_T>::type& tol,
    const spmats<_T,_Index>& A_for_scale)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (!result.converged) return;
    std::vector<scalar_real_type> theta_abs_c1;
    theta_abs_c1.reserve(result.eigenvalues.size());
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++)
        theta_abs_c1.push_back(vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(result.eigenvalues[i])));
    if (!vcp::tsparse::residual_acceptance_check_scaled_(
            result.residuals_absolute, result.residuals_relative, tol,
            theta_abs_c1, matrix_inf_norm_value_<_T,_Index>(A_for_scale))) {
        result.converged = false;
        result.status = "residual_check_failed";
        result.failure_reason =
            "end-of-run exact residual (lambda space) failed acceptance (C-1)";
        result.message = result.failure_reason;
    }
}

// EIG-8 T-3 (D8-3、停止報告裁定 s-3/s-4 の再位置づけ込み): 一般化問題用の
// 共有 scale 版(A x = λ B x)。受理は D8-2 共有ヘルパ
// residual_acceptance_check_generalized_scaled_(既存式 ∨ tol·max(1+|θ|,
// ‖A‖∞ + |θ|·‖B‖∞))。既存式の rel 分岐は Frobenius 後退正規化
// (generalized_eigenpair_relative_residual_norm_value_)を既に持つため、
// 本置換の挙動遷移期待はゼロ — scale 定義の共有集約(契約整合の完成)である。
template <typename _T, typename _Index>
static void lambda_c1_acceptance_gate_(eig_result<_T>& result,
    const typename vcp::tsparse_scalar::real_type<_T>::type& tol,
    const spmats<_T,_Index>& A_for_scale,
    const spmats<_T,_Index>& B_for_scale)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (!result.converged) return;
    std::vector<scalar_real_type> theta_abs_c1;
    theta_abs_c1.reserve(result.eigenvalues.size());
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++)
        theta_abs_c1.push_back(vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(result.eigenvalues[i])));
    if (!vcp::tsparse::residual_acceptance_check_generalized_scaled_(
            result.residuals_absolute, result.residuals_relative, tol,
            theta_abs_c1, matrix_inf_norm_value_<_T,_Index>(A_for_scale),
            matrix_inf_norm_value_<_T,_Index>(B_for_scale))) {
        result.converged = false;
        result.status = "residual_check_failed";
        result.failure_reason =
            "end-of-run exact residual (lambda space) failed acceptance (C-1)";
        result.message = result.failure_reason;
    }
}

// ---------------------------------------------------------------------------
// EIG-6 F-1: si back-half の返却前磨き(条件付き・single-round)
//
// 発火条件(G-0.1 承認 §3.0): 各 back-half の「既存の λ 空間 C-1 最終ゲート」が
// 不合格(status == "residual_check_failed")を出した場合のみ。現行ゲートを
// 通過する経路は磨き経路に入らず、実行軌跡・mv・値はバイト同一に保たれる。
//
// 磨き(D-17a 機構への対処): 返却予定の全対 (θ, v) へ逆反復 1 回
//   w = op.apply(v)(標準: (A−σI)^{-1} v / 一般化: (A−σB)^{-1} B v)
//   w ← w/‖w‖、θ' = Rayleigh 商(標準: wᵀAw / 一般化: wᵀAw / wᵀBw)
// を適用し、磨き後の (θ', w) に対して**同一の**受理ゲートを再判定する。
// λ 空間残差の床 ~|λ|·‖A‖·eps → ~‖A‖·eps に引き下げる(4〜5 桁改善)。
//
// 規律(EIG-6 設計 §3 F-1 (i)〜(v) + G-0.1 付帯 5 点):
//  (i)   single-round: 磨きは 1 回限り。再不合格なら従来どおり
//        residual_check_failed(磨きループ禁止)。
//  (ii)  発火時は返却予定の全対を磨く(部分磨き不可)。磨き後に全対で
//        ゲート再判定 + target 順再確定。
//  (iii) 磨きの solve は mv / linear_solves に 1:1 計上(B-38)。予算超過に
//        なる場合(mv + 対数 > max_iter)は磨かない(B-1: 無料磨きの禁止)。
//  (iv)  決定的(乱数なし)。
//  (v)   単調性ガードは対ごと: 磨き後残差が厳密に改善(certified <)した対のみ
//        採用、悪化・不確定(interval の indeterminate を含む)は磨き前を採用。
//  (vi)  受理式・scale・tol は 1 ビットも変えない(B-40)。ゲートは各 back-half
//        の現行判定そのもの(標準 = 改訂 scale 版 / 一般化 = scale なし版)。
//  (vii) 複素対が返却面にある場合(complex_pair_count != 0 または
//        eigenvalues_imag に非零)は磨かない(実ベクトル逆反復の前提外。
//        正直な不合格がそのまま残る)。
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static bool si_polish_returned_pairs_all_real_(const eig_result<_T>& result)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (result.complex_pair_count != 0) return false;
    if (result.eigenvalues_imag.size() == result.eigenvalues.size()) {
        for (std::size_t i = 0; i < result.eigenvalues_imag.size(); i++) {
            const scalar_real_type im = vcp::tsparse_scalar::abs_value(
                vcp::tsparse_scalar::real_part(result.eigenvalues_imag[i]));
            if (im > scalar_real_type(0)) return false;
        }
    }
    return true;
}

// 磨き本体(標準/一般化共通)。B_ptr == nullptr なら標準問題。
// 戻り値: 磨きを実施したか(mv 計上済みか)。ゲート再判定は呼び出し側で行う。
template <typename _T, typename _Index, class ApplyOp>
static bool si_polish_pairs_once_(eig_result<_T>& result,
                                  const spmats<_T,_Index>& A_csr,
                                  const spmats<_T,_Index>* B_ptr,
                                  const eig_options<_T>& options,
                                  ApplyOp& apply_si)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t P = result.eigenvalues.size();
    if (P == 0) return false;
    if (result.eigenvectors.size() != P) return false;
    if (result.residuals_absolute.size() != P) return false;
    if (!si_polish_returned_pairs_all_real_<_T,_Index>(result)) return false;
    // (iii) 予算規律: 磨き分の solve が予算を超えるなら磨かない(正直不合格を維持)
    if (result.matrix_vector_products + P > options.max_iter) return false;

    for (std::size_t p = 0; p < P; p++) {
        std::vector<_T> w;
        apply_si(result.eigenvectors[p], w);          // 逆反復 1 回 = 1 solve
        result.matrix_vector_products += 1;           // (iii) B-38: 1:1 計上
        const scalar_real_type wn = norm_value_<_T,_Index>(w);
        if (!(wn > scalar_real_type(0))) continue;    // 退化(certified >0 のみ採用側へ)
        for (std::size_t i = 0; i < w.size(); i++) w[i] = w[i] / _T(wn);
        // Rayleigh 商(標準: wᵀAw、一般化: wᵀAw / wᵀBw)
        std::vector<_T> Aw = A_csr.mul_vec(w);
        _T num = _T(0);
        for (std::size_t i = 0; i < w.size(); i++) num += w[i] * Aw[i];
        _T theta_new = num;
        if (B_ptr != 0) {
            std::vector<_T> Bw = B_ptr->mul_vec(w);
            _T den = _T(0);
            for (std::size_t i = 0; i < w.size(); i++) den += w[i] * Bw[i];
            const scalar_real_type den_abs = vcp::tsparse_scalar::abs_value(
                vcp::tsparse_scalar::real_part(den));
            if (!(den_abs > scalar_real_type(0))) continue;   // B 内積退化: 磨き前を採用
            theta_new = num / den;
        }
        // (v) 対ごと単調性ガード: certified に改善した対のみ採用
        const scalar_real_type res_new = (B_ptr != 0)
            ? generalized_eigenpair_residual_norm_value_<_T,_Index>(A_csr, *B_ptr, theta_new, w)
            : eigenpair_residual_norm_value_<_T,_Index>(A_csr, theta_new, w);
        if (res_new < result.residuals_absolute[p]) {
            result.eigenvalues[p]  = theta_new;
            result.eigenvectors[p] = w;
        }
    }
    // (ii) 磨き後の全対で残差を再計算(採用/非採用が混在しても一様に)
    if (B_ptr != 0) {
        result.residuals_absolute = generalized_eigenpair_residuals_<_T,_Index>(A_csr, *B_ptr, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = generalized_eigenpair_relative_residuals_<_T,_Index>(A_csr, *B_ptr, result.eigenvalues, result.eigenvectors);
    } else {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A_csr, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A_csr, result.eigenvalues, result.eigenvectors);
    }
    if (!result.residuals_absolute.empty()) {
        result.residual_norm_absolute =
            *std::max_element(result.residuals_absolute.begin(), result.residuals_absolute.end());
    }
    // (ii) target 順の再確定(θ' の微小変化での順序逆転に備える)
    vcp::tsparse_eigen_selection::sort_eigenpairs_by_target(
        result.eigenvalues, result.eigenvectors,
        result.residuals_absolute, result.residuals_relative,
        options.target, options.shift);
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    return true;
}

// 標準問題用: 現行の改訂 scale 版 C-1 ゲート(lambda_c1_acceptance_gate_ の
// scale 版と同一式)で再判定する。合格なら converged=true へ昇格。
template <typename _T, typename _Index, class ApplyOp>
static void si_polish_rescue_standard_(eig_result<_T>& result,
                                       const spmats<_T,_Index>& self,
                                       const eig_options<_T>& options,
                                       ApplyOp& apply_si)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (result.converged) return;
    if (result.status != "residual_check_failed") return;   // 発火条件(既存ゲート不合格のみ)
    spmats<_T,_Index> A = self.as_csr();
    if (!si_polish_pairs_once_<_T,_Index>(result, A, static_cast<const spmats<_T,_Index>*>(0), options, apply_si))
        return;
    std::vector<scalar_real_type> theta_abs_c1;
    theta_abs_c1.reserve(result.eigenvalues.size());
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++)
        theta_abs_c1.push_back(vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(result.eigenvalues[i])));
    if (vcp::tsparse::residual_acceptance_check_scaled_(
            result.residuals_absolute, result.residuals_relative, options.tol,
            theta_abs_c1, matrix_inf_norm_value_<_T,_Index>(A))) {
        result.converged = true;
        result.status = "converged";
        result.failure_reason.clear();
        result.message = "converged (C-1 rescue: one-step inverse-iteration polish; EIG-6 F-1)";
    } else {
        // (i) single-round: 再不合格は従来どおり(磨き改善分は診断として残る — B-15)
        result.failure_reason =
            "end-of-run exact residual (lambda space) failed acceptance (C-1); one-step polish attempted (EIG-6 F-1)";
        result.message = result.failure_reason;
    }
}

// 一般化問題用: EIG-8 T-3(D8-3、G-0.1 承認事項 d)— 再判定も最終ゲートと
// 同一の D8-2 共有 scale 形に統一する(発火条件が共有形不合格である以上、
// 再判定だけ旧形では非対称が残るため。EIG-6 G-0.1 (v) の「統一しない」判断は
// D8-3 の統一タスク承認で上書き)。
template <typename _T, typename _Index, class ApplyOp>
static void si_polish_rescue_generalized_(eig_result<_T>& result,
                                          const spmats<_T,_Index>& self,
                                          const spmats<_T,_Index>& B,
                                          const eig_options<_T>& options,
                                          ApplyOp& apply_si)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (result.converged) return;
    if (result.status != "residual_check_failed") return;   // 発火条件(既存ゲート不合格のみ)
    spmats<_T,_Index> A = self.as_csr();
    spmats<_T,_Index> B_csr = B.as_csr();
    if (!si_polish_pairs_once_<_T,_Index>(result, A, &B_csr, options, apply_si))
        return;
    std::vector<scalar_real_type> theta_abs_g;
    theta_abs_g.reserve(result.eigenvalues.size());
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++)
        theta_abs_g.push_back(vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(result.eigenvalues[i])));
    if (vcp::tsparse::residual_acceptance_check_generalized_scaled_(
            result.residuals_absolute, result.residuals_relative, options.tol,
            theta_abs_g, matrix_inf_norm_value_<_T,_Index>(A),
            matrix_inf_norm_value_<_T,_Index>(B_csr))) {
        result.converged = true;
        result.status = "converged";
        result.failure_reason.clear();
        result.message = "converged (C-1 rescue: one-step inverse-iteration polish; EIG-6 F-1)";
    } else {
        result.failure_reason =
            "end-of-run exact residual (lambda space) failed acceptance (C-1); one-step polish attempted (EIG-6 F-1)";
        result.message = result.failure_reason;
    }
}

// ---------------------------------------------------------------------------
// EIG-6 F-2' θ シフト磨き(G-2B.1 承認事項・E-A1 si_lanczos front 限定)
//
// 発火条件: converged だが「scaled-abs 分岐(res_abs ≤ tol·max(1+|θ|,‖A‖∞))を
// 満たさない対がある」場合のみ。契約ゲート(無変更)は rel 分岐(Frobenius
// 相対)で受理し得るが、μ ロック品質の λ 空間残差床 ~|λ|·‖A‖·r_μ は近接
// クラスタで scaled-abs を外れ得る(grid_large 実測 3.5e-9 > 1.2e-9)。
// σ=0 の逆反復磨き(F-1)は近接固有値の誤差成分をほぼ減衰させない
// (減衰比 λ_i/λ_j ≈ 0.997)ため、対ごとの Rayleigh 商 θ_p をシフトに使う
// 古典的な shifted inverse iteration で 1 回で eps 床へ落とす:
//   w = (A − θ_p I)^{-1} v_p(E-A1 sparse LU を対ごとに新規分解)、
//   w ← w/‖w‖、θ' = wᵀAw、対ごと単調性ガード(certified に改善した対のみ採用)。
// 規律: single-round(P4 同様の防御)・全対走査・solve は mv/linear_solves に
// 1:1 計上(P1)・決定的・複素対除外(P3)・予算内のみ(mv + 対数 ≤ max_iter)。
// 分解は mv 通貨の対象外(E-A1 と同じ扱い)だが、対ごとの追加分解 ≤ k 回を
// 診断メッセージに明記する。θ_p が固有値に極近く (A−θI) が数値特異でも、
// 逆反復の誤差は目的固有方向に落ちる(古典理論)— 分解失敗対は磨き前を採用。
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static typename std::enable_if<std::is_signed<_Index>::value, void>::type
si_polish_shifted_rescue_(eig_result<_T>& result,
                          const spmats<_T,_Index>& self,
                          const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    typedef vcp::tsparse::lu_shift_invert_operator<spmats<_T,_Index> > LUOp;
    if (!result.converged) return;
    const std::size_t P = result.eigenvalues.size();
    if (P == 0) return;
    if (result.eigenvectors.size() != P) return;
    if (result.residuals_absolute.size() != P) return;
    if (!si_polish_returned_pairs_all_real_<_T,_Index>(result)) return;
    spmats<_T,_Index> A = self.as_csr();
    const scalar_real_type anorm = matrix_inf_norm_value_<_T,_Index>(A);
    // 発火判定: scaled-abs 全対クリーンなら何もしない(既存 OK 面はバイト不変)
    bool all_clean = true;
    for (std::size_t i = 0; i < P; i++) {
        const scalar_real_type th = vcp::tsparse_scalar::abs_value(
            vcp::tsparse_scalar::real_part(result.eigenvalues[i]));
        const scalar_real_type sc = vcp::tsparse::c1_revised_scale_(th, anorm);
        if (!(result.residuals_absolute[i] <= options.tol * sc)) { all_clean = false; break; }
    }
    if (all_clean) return;
    if (result.matrix_vector_products + P > options.max_iter) return;   // 予算規律

    std::size_t polish_factorizations = 0;
    for (std::size_t p = 0; p < P; p++) {
        LUOp op(self, result.eigenvalues[p], options.shift_invert_lu);
        polish_factorizations++;
        if (!op.factorization_ok()) continue;              // 磨き前を採用
        std::vector<_T> w;
        op.apply(result.eigenvectors[p], w);               // 1 solve
        result.matrix_vector_products += 1;                // P1: 1:1 計上
        result.linear_solves += 1;
        const scalar_real_type wn = norm_value_<_T,_Index>(w);
        if (!(wn > scalar_real_type(0))) continue;
        for (std::size_t i = 0; i < w.size(); i++) w[i] = w[i] / _T(wn);
        std::vector<_T> Aw = A.mul_vec(w);
        _T num = _T(0);
        for (std::size_t i = 0; i < w.size(); i++) num += w[i] * Aw[i];
        const _T theta_new = num;
        const scalar_real_type res_new =
            eigenpair_residual_norm_value_<_T,_Index>(A, theta_new, w);
        if (res_new < result.residuals_absolute[p]) {      // 対ごと単調性(certified)
            result.eigenvalues[p]  = theta_new;
            result.eigenvectors[p] = w;
        }
    }
    result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    if (!result.residuals_absolute.empty()) {
        result.residual_norm_absolute =
            *std::max_element(result.residuals_absolute.begin(), result.residuals_absolute.end());
    }
    vcp::tsparse_eigen_selection::sort_eigenpairs_by_target(
        result.eigenvalues, result.eigenvectors,
        result.residuals_absolute, result.residuals_relative,
        options.target, options.shift);
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    if (polish_factorizations > 0) {
        // (f-3) 追加分解回数の正式公開(構造化フィールド)+ 人間可読の注記
        result.polish_factorizations = polish_factorizations;
        result.message = "converged (theta-shifted inverse-iteration polish; EIG-6 F-2'; "
            "extra factorizations=" + std::to_string(polish_factorizations) + ")";
    }
}

template <typename _T, typename _Index>
static typename std::enable_if<!std::is_signed<_Index>::value, void>::type
si_polish_shifted_rescue_(eig_result<_T>&,
                          const spmats<_T,_Index>&,
                          const eig_options<_T>&)
{
    // unsigned Index は sparse LU 不対応(E-A1 と同じ制約)— 磨きなし(現状維持)
}

// Field-compatible package (mirrors the consumed subset of the old
// arnoldi_result_package) produced by the KS core in mu space.
template <typename _T>
struct ks_mu_pkg_ {
    typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
    std::vector<_T> eigenvalues;                     // mu space
    std::vector<std::vector<_T> > eigenvectors;
    std::vector<R> history_abs;
    std::vector<R> history_rel;
    std::vector<std::pair<R, R> > complex_eigenvalues;  // lambda space (re, +im)
    std::size_t iterations;
    std::size_t mv_count;
    std::size_t converged_count;
    std::size_t used_subspace_dim;
    std::string breakdown_reason;
    std::string failure_reason;
    // Always false: the KS core already folds the complex-window verdict into
    // `converged` + failure_reason (D3-2); the legacy has_complex gating of
    // the old back-halves must not re-fire on mere complex diagnostics.
    bool has_complex;
    bool converged;   // KS converged AND back-half C-2 evidence re-check
    ks_mu_pkg_() : iterations(0), mv_count(0), converged_count(0),
                   used_subspace_dim(0), has_complex(false), converged(false) {}
};

template <typename _T, class Apply, class LambdaLockGate>
static ks_mu_pkg_<_T> ks_mu_core_drive_impl_(const std::size_t /*n*/, const std::size_t /*k*/,
                                             const eig_options<_T>& /*options*/,
                                             const typename vcp::tsparse_scalar::real_type<_T>::type& /*sigma*/,
                                             Apply& /*apply_si*/,
                                             LambdaLockGate /*lambda_lock_gate*/,
                                             const bool /*use_lambda_lock_gate*/,
                                             std::size_t* /*lambda_gate_products*/,
                                             std::true_type /* is_complex */)
{
    // Fence for complex instantiations (runtime dispatch never reaches here).
    ks_mu_pkg_<_T> pkg;
    pkg.failure_reason = "shift-invert krylov_schur core: complex scalar unsupported";
    return pkg;
}

template <typename _T, class Apply, class LambdaLockGate>
static ks_mu_pkg_<_T> ks_mu_core_drive_impl_(const std::size_t n, const std::size_t k,
                                             const eig_options<_T>& options,
                                             const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                             Apply& apply_si,
                                             LambdaLockGate lambda_lock_gate,
                                             const bool use_lambda_lock_gate,
                                             std::size_t* lambda_gate_products,
                                             std::false_type /* is_complex */)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
    eig_options<_T> mu_opts = options;
    mu_opts.method = eig_solver_method::arnoldi;
    mu_opts.target = eig_target::largest_magnitude;   // inner target (mu space)
    mu_opts.shift = R(0);
    mu_opts.use_shift = false;
    // options.max_iter passes through as the TOTAL inner-apply budget
    // (D-6 unification; the old max_iter*(k+1) restart expansion is deleted).
    // EIG-8 T-1: μ 面 opt-in の λ 形式 lock ゲートを KS コアへ転送(既定は
    // ks_lambda_lock_gate_none + false = 従来挙動バイト同一)。
    vcp::tsparse_experimental::krylov_schur_result<_T> d =
        vcp::tsparse_experimental::krylov_schur_eigs_with_diagnostics<Apply, _T, LambdaLockGate>(
            apply_si, n, k, mu_opts,
            lambda_lock_gate, use_lambda_lock_gate, lambda_gate_products);

    ks_mu_pkg_<_T> pkg;
    pkg.eigenvalues = d.eigs.eigenvalues;
    pkg.eigenvectors = d.eigs.eigenvectors;
    pkg.history_abs = d.eigs.residual_history_absolute;
    pkg.history_rel = d.eigs.residual_history_relative;
    pkg.iterations = d.eigs.iterations;
    pkg.mv_count = d.eigs.matrix_vector_products;
    pkg.converged_count = d.eigs.converged_count;
    pkg.used_subspace_dim = d.eigs.used_subspace_dim;
    pkg.breakdown_reason = d.eigs.breakdown_reason;
    pkg.failure_reason = d.eigs.failure_reason;
    // complex diagnostics: transform mu -> lambda = 1/mu + sigma when the
    // modulus is certifiably positive (pairs are stored adjacent +s/-s in the
    // KS result; export one (re, +im) entry per pair, legacy convention)
    for (std::size_t i = 0; i + 1 < d.eigs.complex_eigenvalues.size(); i += 2) {
        const R a = d.eigs.complex_eigenvalues[i].real();
        const R b = d.eigs.complex_eigenvalues[i].imag();
        const R den = a * a + b * b;
        if (den > R(0)) {
            const R re = a / den + sigma;
            R im = -b / den;
            if (im < R(0)) im = -im;
            pkg.complex_eigenvalues.push_back(std::pair<R, R>(re, im));
        } else {
            R babs = b; if (babs < R(0)) babs = -babs;
            pkg.complex_eigenvalues.push_back(std::pair<R, R>(a, babs));
        }
    }
    // D3-4: re-evaluate the shared C-2 check on the exported evidence
    bool evid_ok = false;
    if (d.eigs.converged && d.c2_evidence.exported && d.c2_evidence.fresh) {
        std::vector<R> mu_locked;
        for (std::size_t i = 0; i < d.eigs.eigenvalues.size(); i++)
            mu_locked.push_back(vcp::tsparse_scalar::real_part(d.eigs.eigenvalues[i]));
        evid_ok = vcp::tsparse::honest_termination_check_complex_<R>(
            d.c2_evidence.candidate_real, d.c2_evidence.candidate_imag,
            d.c2_evidence.candidate_converged, mu_locked,
            mu_locked.size(), eig_target::largest_magnitude, R(0));
        if (!evid_ok && pkg.failure_reason.empty())
            pkg.failure_reason =
                "back-half C-2 re-check on exported evidence failed";
    }
    pkg.converged = d.eigs.converged && evid_ok;
    return pkg;
}

template <typename _T, class Apply,
          class LambdaLockGate = vcp::tsparse_experimental::ks_lambda_lock_gate_none>
static ks_mu_pkg_<_T> ks_mu_core_drive_(const std::size_t n, const std::size_t k,
                                        const eig_options<_T>& options,
                                        const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                        Apply& apply_si,
                                        LambdaLockGate lambda_lock_gate = LambdaLockGate(),
                                        const bool use_lambda_lock_gate = false,
                                        std::size_t* lambda_gate_products = 0)
{
    return ks_mu_core_drive_impl_<_T, Apply, LambdaLockGate>(n, k, options, sigma, apply_si,
        lambda_lock_gate, use_lambda_lock_gate, lambda_gate_products,
        typename std::integral_constant<bool, spmatrix_is_complex<_T>::value>::type());
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
    (void)k;
    validate_eig_input_<_T,_Index>(A, routine);
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
    const scalar_real_type pivot_tol = vcp::tsparse_scalar::epsilon<scalar_real_type>();
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
        + vcp::tsparse_scalar::epsilon<scalar_real_type>();
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
    // EIG-3 T-4 (D-6, Q1 ruling 2026-07-06): max_iter is the TOTAL mv budget on
    // every path.  Pure unit conversion to restarts with the worst per-restart
    // cost 2*sdim (expansion sdim + lock-scan exact residuals <= sdim;
    // tsparse_lanczos.hpp mv_count++ sites); the inflation term "+ k + 1" is
    // deleted.  The restart loop is INCLUSIVE (restart <= max_restarts), so a
    // floor of 0 already yields one restart; no max(1,.) floor (it would
    // double the minimum work and break the bound for max_iter < 2*sdim).
    // Machine-checked overshoot bound: mv <= max_iter + 2*sdim.
    const std::size_t max_restarts_l = (sdim > 0)
        ? options.max_iter / (2 * sdim)
        : options.max_iter;
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
                                                         const eig_options<_T>& options_in)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    // EIG-4 T-2 (R-4 と同趣旨): 複素 T の auto_select は旧既定(lanczos =
    // hermitian 経路)へ正規化する(挙動は旧既定と同一)。実 T 専用の
    // thick_restart_lanczos / krylov_schur は正直拒否。以下の本体は
    // options 名の付け替えのみで無変更(B-26)。
    eig_options<_T> options_norm = options_in;
    if (options_norm.method == eig_solver_method::auto_select)
        options_norm.method = eig_solver_method::lanczos;
    if (k > 0
     && (options_norm.method == eig_solver_method::thick_restart_lanczos
      || options_norm.method == eig_solver_method::krylov_schur)) {
        eig_result<_T> result;
        result.requested_count = k;
        result.converged = false;
        result.status = "unsupported_complex_non_hermitian";
        result.failure_reason =
            "thick_restart_lanczos / krylov_schur are real-scalar drivers;"
            " for complex Hermitian matrices use method=lanczos";
        result.message = result.failure_reason;
        result.method = options_norm.method;
        result.used_method = eig_method_to_string_<_T,_Index>(options_norm.method);
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    const eig_options<_T>& options = options_norm;
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
        // EIG-G1 entrance guard (TASK-3 ledger §E policy;
        // eig_dense_complex_legacy_issue.md): this branch was the one
        // runtime-reachable entrance by which complex-scalar matrices flowed
        // silently into the legacy complex qr_eig_dense dense path (D-13 /
        // D-15 unrepaired there).  Complex scalars are outside the mats /
        // spmats design scope, so reject explicitly at the user-facing API
        // (honest termination: misuse -> throw).  The legacy path body
        // (dense_eig_ complex branch / qr_eig_dense) is kept, merely
        // unreachable for complex scalars.
        vcp::throw_error<vcp::invalid_argument>(
            "spmats::eigs(complex, dense_fallback_explicit): complex scalars"
            " are unsupported; use a real coupled formulation"
            " (see VCP_task_list §E policy)");
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
// E-A1: shared shift-invert back halves and sparse_lu dispatch helpers.
//
// The *_drive_ helpers are byte-preserving mechanical extractions of the
// legacy back halves of shift_invert_lanczos_eigs_ / shift_invert_arnoldi_
// eigs_ / generalized_shift_invert_arnoldi_eigs_ (Krylov driver call ->
// lambda = sigma + 1/mu inversion -> result assembly -> residuals).  Both
// the legacy ILU(0)+GMRES branch and the E-A1 sparse_lu branch call them,
// so the post-processing stays single-source (B-10: the extraction does not
// change any of the inversion / residual / selection semantics).
//
// The *_sparse_lu_ helpers implement the E-A1 default path: build a
// vcp::tsparse::lu_shift_invert_operator (factorize (A - sigma*B) once,
// direct-solve per apply), check factorization_ok() immediately, and on
// failure return status="factorization_failed" with the sparse_lu
// diagnostics (D-4 / B-9: no silent fallback to ILU+GMRES).  They are
// SFINAE-split on the signedness of _Index because sparse LU requires a
// signed Index (same precedent as dispatch_sparse_lu_ in spmats_lss.hpp);
// for unsigned _Index a controlled factorization_failed result directs the
// caller to the ilu0_gmres compatibility solver.
// ---------------------------------------------------------------------------

// --- back half #1: standard shift-invert Lanczos --------------------------
// EIG-10 D-17d 委譲用の前方宣言(back half #2 は本ファイル後方で定義)
template <typename _T, typename _Index, class Apply>
static eig_result<_T> shift_invert_arnoldi_drive_(const spmats<_T,_Index>& self,
                                                    const std::size_t k,
                                                    const eig_options<_T>& options,
                                                    const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                                    Apply& apply_si,
                                                    bool& pkg_converged);

template <typename _T, typename _Index, class Apply>
static eig_result<_T> shift_invert_lanczos_drive_(const spmats<_T,_Index>& self,
                                                    const std::size_t k,
                                                    const eig_options<_T>& options,
                                                    const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                                    Apply& apply_si)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    const std::size_t sdim = (options.subspace_dim == 0)
        ? std::max(k + 5, std::min(n, std::size_t(30)))
        : options.subspace_dim;
    // EIG-3 T-4 (D-6, Q1 ruling 2026-07-06): max_iter is the TOTAL mv budget on
    // every path.  Pure unit conversion to restarts with the worst per-restart
    // cost 2*sdim (expansion sdim + lock-scan exact residuals <= sdim;
    // tsparse_lanczos.hpp mv_count++ sites); the inflation term "+ k + 1" is
    // deleted.  The restart loop is INCLUSIVE (restart <= max_restarts), so a
    // floor of 0 already yields one restart; no max(1,.) floor (it would
    // double the minimum work and break the bound for max_iter < 2*sdim).
    // Machine-checked overshoot bound: mv <= max_iter + 2*sdim.
    const std::size_t max_restarts_si = (sdim > 0)
        ? options.max_iter / (2 * sdim)
        : options.max_iter;
    // EIG-6 F-2' (G-2.1 承認): si_lanczos back-half は μ コアの opt-in
    // リスタート制御(b' warm start / d' 走査キャップ)+ λ 形式 lock 併記 (e)
    // を有効化する。λ ゲートは契約 C-1 の共有ヘルパ
    // (eigenpair_residual_norm_value_ = λ 空間厳密残差 / c1_revised_scale_ =
    // 改訂 scale の唯一定義)をそのまま呼ぶ(e-1。受理は scaled 分岐のみ =
    // 契約受理の厳密部分集合。最終ゲート lanczos_package_to_result_ は無変更の
    // まま全対を再検査する — e-3)。A·x は mv に 1:1 計上 + 別建てカウント(e-2)。
    spmats<_T,_Index> A_for_gate = self.as_csr();
    const scalar_real_type anorm_gate = matrix_inf_norm_value_<_T,_Index>(A_for_gate);
    struct SiLambdaGate {
        const spmats<_T,_Index>* A;
        scalar_real_type anorm;
        scalar_real_type tol;
        scalar_real_type sigma;
        bool operator()(const std::vector<_T>& x, const _T& mu, std::size_t& mv) const {
            const scalar_real_type mu_re = vcp::tsparse_scalar::real_part(mu);
            if (!(vcp::tsparse_scalar::abs_value(mu_re) > scalar_real_type(0))) return false;
            const scalar_real_type lam = scalar_real_type(1) / mu_re + sigma;
            const scalar_real_type r_abs =
                eigenpair_residual_norm_value_<_T,_Index>(*A, _T(lam), x);
            mv += 1;   // λ 判定の A·x 積(B-38: 1:1 計上)
            const scalar_real_type scale = vcp::tsparse::c1_revised_scale_(
                vcp::tsparse_scalar::abs_value(lam), anorm);
            return r_abs <= tol * scale;   // certified-≤(interval は不確定を失敗側へ)
        }
    } si_lambda_gate = { &A_for_gate, anorm_gate, options.tol, sigma };
    std::size_t lambda_gate_count = 0;
    auto pkg = vcp::tsparse_lanczos::lanczos_eigs_standard_si_<_T, Apply, SiLambdaGate>(
        n, k, sdim, max_restarts_si, options.tol,
        options.random_seed, options.random_start,
        eig_target::largest_magnitude, scalar_real_type(0), options.compute_residual_history,
        apply_si, si_lambda_gate, &lambda_gate_count);

    for (std::size_t i = 0; i < pkg.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(pkg.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0)) {
            pkg.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
        }
    }
    eig_result<_T> result = lanczos_package_to_result_<_T,_Index>(pkg, self, k, eig_solver_method::shift_invert_lanczos, options.tol);
    result.lambda_gate_products = lambda_gate_count;   // (e-2) 別建て計上
    // EIG-6 F-1: 条件付き磨き(既存 C-1 ゲート不合格時のみ発火)。sparse_lu /
    // legacy ilu0_gmres の両フロントを被覆(linear_solves はフロント側が
    // apply 内カウンタ / op から再同期するため磨き分も自動計上される)。
    si_polish_rescue_standard_<_T,_Index>(result, self, options, apply_si);
    // -----------------------------------------------------------------------
    // EIG-10 D-17d 委譲(G-0.1 承認 案 (b)・(c-1) 改訂・additive route)。
    //
    // 発火署名(すべて既存フィールドの読み取りのみ・新定数ゼロ):
    //   現行パイプライン(μ コア → C-1 λ 空間最終ゲート → F-1 磨き)が
    //   完全失敗(status == residual_check_failed)し、かつ残予算がある場合のみ。
    //   凍結スイートに本署名で終わる si_lanczos セルは存在しない(G-0.1 §8)
    //   ため、既存挙動はビット不変(休眠は構造的 — G-2.3 で全数機械証明)。
    //
    // 機構(D-17d): μ 空間絶対 lock 基準(r_μ ≤ tol)は |μ| ≪ 1 で分解能
    // 不正直(過受理)となり、core は honest termination で予算 ~1% の時点で
    // converged 終了 → C-1 が正直棄却 → 再入機構不在で残予算 99% が死蔵される。
    // 委譲先の KS μ コアは同条件の縮退多重度を予算内で完全発見できることが
    // 実証済み(EIG-10 Phase 0 §5)。
    //
    // 規律:
    //  - 1 回限り・決定的(乱数なし・ループ禁止)。
    //  - all-or-nothing((c-1)): 委譲結果が「全ゲート合格の converged ∧
    //    返却全対実 ∧ 複製方向なし(R4 と同一ヘルパ)」の場合のみ丸ごと採用。
    //    それ以外は元の正直結果をそのまま返す(対集合の混合・部分採用禁止)。
    //  - allow_complex_pairs=true の内部指定は「過渡複素窓の通過許可」であって
    //    「複素対の返却許可」ではない: FP 丸め非対称の対称宣言行列で KS が
    //    D3-2 の正直拒否により空振りするのを防ぐためのみ。複素対を含む委譲
    //    結果は全実ガードが不採用に倒すため、structure=symmetric の呼び出しに
    //    複素対が漏れる経路は構造的に存在しない。
    //  - mv は採否に依らず実消費を 1:1 計上(B-38)し、ks_rescue_products に
    //    別建て記録((c-2) の恒等式 mv_total == si 分 + 委譲分)。KS コアは
    //    mv ≤ 残予算を厳守するため mv_total ≤ max_iter(B-1)。
    //    linear_solves はフロントが op カウンタから再同期(磨きと同じ経路)。
    // -----------------------------------------------------------------------
    if (!result.converged
        && result.status == "residual_check_failed"
        && options.max_iter > result.matrix_vector_products) {
        eig_options<_T> ks_opts = options;
        ks_opts.max_iter = options.max_iter - result.matrix_vector_products;
        ks_opts.allow_complex_pairs = true;   // 過渡複素窓の通過許可(上記規律)
        bool ks_conv = false;
        eig_result<_T> alt = shift_invert_arnoldi_drive_<_T,_Index>(
            self, k, ks_opts, sigma, apply_si, ks_conv);
        alt.converged = ks_conv && alt.eigenvalues.size() >= k;
        // 帰結配線は E-A1 arnoldi front と同一(B-43 再利用): C-1 改訂 scale
        // 最終ゲート → F-1 条件付き磨き(不合格時のみ発火・same-helper。
        // 磨き solve は apply_si 経由 = ls はフロント再同期で自動整合、
        // mv は helper 内で 1:1 計上され下の rescue 合算に含まれる)。
        lambda_c1_acceptance_gate_<_T,_Index>(alt, options.tol, self);
        si_polish_rescue_standard_<_T,_Index>(alt, self, ks_opts, apply_si);
        result.ks_rescue_products = alt.matrix_vector_products;
        result.matrix_vector_products += alt.matrix_vector_products;
        result.lambda_gate_products += alt.lambda_gate_products;
        if (alt.converged && si_polish_returned_pairs_all_real_<_T,_Index>(alt)) {
            // 採用前の基底整形: 返却集合の MGS 直交正規化(2 パス)。縮退固有
            // 空間の基底の向きを整えるだけで、対集合の混合・部分採用ではない
            // ((c-1) の all-or-nothing は集合単位で維持)。整形後に残差を
            // 再計算し、最終ゲートと同一の受理式で全対を再判定してから採用する
            // (整形で品質が立たない場合・基底が退化する場合は丸ごと不採用)。
            typedef typename vcp::tsparse_scalar::real_type<_T>::type R__;
            std::vector<std::vector<_T> > ortho = alt.eigenvectors;
            bool ortho_ok = (ortho.size() == alt.eigenvalues.size()) && !ortho.empty();
            for (std::size_t i = 0; ortho_ok && i < ortho.size(); i++) {
                for (int pass = 0; pass < 2; pass++) {
                    for (std::size_t j = 0; j < i; j++) {
                        const R__ c = vcp::tsparse_scalar::real_dot_value(ortho[j], ortho[i]);
                        for (std::size_t t = 0; t < ortho[i].size(); t++)
                            ortho[i][t] -= _T(c) * ortho[j][t];
                    }
                }
                const R__ nv = vcp::tsparse_scalar::real_norm_value(ortho[i]);
                if (!(nv > R__(0))) { ortho_ok = false; break; }   // certified > 0 のみ採用側
                for (std::size_t t = 0; t < ortho[i].size(); t++)
                    ortho[i][t] = ortho[i][t] / _T(nv);
            }
            if (ortho_ok) {
                std::vector<R__> res_abs_o = eigenpair_residuals_<_T,_Index>(
                    A_for_gate, alt.eigenvalues, ortho);
                std::vector<R__> res_rel_o = eigenpair_relative_residuals_<_T,_Index>(
                    A_for_gate, alt.eigenvalues, ortho);
                std::vector<R__> theta_abs_o;
                theta_abs_o.reserve(alt.eigenvalues.size());
                for (std::size_t i = 0; i < alt.eigenvalues.size(); i++)
                    theta_abs_o.push_back(vcp::tsparse_scalar::abs_value(
                        vcp::tsparse_scalar::real_part(alt.eigenvalues[i])));
                if (vcp::tsparse::residual_acceptance_check_scaled_(
                        res_abs_o, res_rel_o, options.tol, theta_abs_o, anorm_gate)
                    && !vcp::tsparse::returned_pair_duplicate_direction_found_<_T>(ortho)) {
                    result.eigenvalues          = alt.eigenvalues;
                    result.eigenvalues_imag     = alt.eigenvalues_imag;
                    result.complex_pair_count   = alt.complex_pair_count;
                    result.complex_eigenvalues  = alt.complex_eigenvalues;
                    result.eigenvectors         = ortho;
                    result.residuals_absolute   = res_abs_o;
                    result.residuals_relative   = res_rel_o;
                    if (!res_abs_o.empty()) {
                        result.residual_norm_absolute = *std::max_element(
                            res_abs_o.begin(), res_abs_o.end());
                    }
                    result.converged        = true;
                    result.converged_count  = alt.converged_count;
                    result.status           = "converged";
                    result.failure_reason.clear();
                    result.breakdown_reason.clear();
                    result.message = "converged (D-17d rescue: ks_mu_core delegation; EIG-10)";
                    result.iterations += alt.iterations;
                    // method / used_method / used_subspace_dim は公開面
                    // (shift_invert_lanczos)のまま(G-0.1 §6.1: enum・文字列追加なし)。
                    populate_real_complex_eigenvalues_<_T,_Index>(result);
                    set_result_counts_<_T,_Index>(result, k);
                }
            }
        }
        // 不採用: 元の正直結果を維持(値・status 不変。消費計上のみ反映済み)。
    }
    return result;
}

// --- back half #2: standard shift-invert Arnoldi ---------------------------
template <typename _T, typename _Index, class Apply>
static eig_result<_T> shift_invert_arnoldi_drive_(const spmats<_T,_Index>& self,
                                                    const std::size_t k,
                                                    const eig_options<_T>& options,
                                                    const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                                    Apply& apply_si,
                                                    bool& pkg_converged)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    // EIG-8 T-1 (e): μ 面 opt-in の λ 形式 lock ゲート(標準問題)。
    // 契約 C-1 の改訂 scale 受理式そのもの(共有ヘルパ
    // eigenpair_residual_norm_value_ / c1_revised_scale_ — e-1)。
    // A·x は mv に 1:1 計上 + lambda_gate_products 別建て(e-2)。
    // certified-≤(interval は不確定を失敗側へ)。
    const scalar_real_type anorm_gate_ks = matrix_inf_norm_value_<_T,_Index>(A);
    struct KsLambdaGate {
        const spmats<_T,_Index>* A;
        scalar_real_type anorm;
        scalar_real_type tol;
        scalar_real_type sigma;
        bool operator()(const std::vector<_T>& x, const _T& mu, std::size_t& mv) const {
            const scalar_real_type mu_re = vcp::tsparse_scalar::real_part(mu);
            if (!(vcp::tsparse_scalar::abs_value(mu_re) > scalar_real_type(0))) return false;
            const scalar_real_type lam = scalar_real_type(1) / mu_re + sigma;
            const scalar_real_type r_abs =
                eigenpair_residual_norm_value_<_T,_Index>(*A, _T(lam), x);
            mv += 1;   // λ 判定の A·x 積(B-24: 1:1 計上)
            const scalar_real_type scale = vcp::tsparse::c1_revised_scale_(
                vcp::tsparse_scalar::abs_value(lam), anorm);
            return r_abs <= tol * scale;
        }
    } ks_lambda_gate = { &A, anorm_gate_ks, options.tol, sigma };
    std::size_t ks_lambda_gate_count = 0;
    // EIG-3 T-3: KS core in mu space (old arnoldi core + max_iter*(k+1)
    // expansion deleted; options.max_iter = total inner-apply budget)
    const ks_mu_pkg_<_T> pkg =
        ks_mu_core_drive_<_T>(n, k, options, sigma, apply_si,
                              ks_lambda_gate, true, &ks_lambda_gate_count);
    const std::size_t sdim = pkg.used_subspace_dim;
    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::shift_invert_arnoldi;
    result.lambda_gate_products = ks_lambda_gate_count;   // (e-2) 別建て計上
    result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_arnoldi);
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
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
    for (std::size_t i = 0; i < pkg.complex_eigenvalues.size(); i++) {
        result.complex_eigenvalues.push_back(
            typename eig_result<_T>::eigenvalue_type(
                pkg.complex_eigenvalues[i].first,
                pkg.complex_eigenvalues[i].second));
    }
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    // EIG-8 R4(裁定 s-2・必須併設): 返却対の複製方向監査(demote-only)。
    // λ ゲート発火ラン限定(tol regime は無評価 = 挙動不変)。検出時は正直
    // not_converged へ降格するのみ(間引き・自動修復は C-3 抵触で禁止)。
    if (pkg.converged && ks_lambda_gate_count > 0 &&
        vcp::tsparse::returned_pair_duplicate_direction_found_<_T>(result.eigenvectors)) {
        if (result.failure_reason.empty())
            result.failure_reason = "duplicate direction among returned pairs "
                "(ghost-copy audit, EIG-8 R4): honest not_converged";
        pkg_converged = false;
    } else {
        pkg_converged = pkg.converged;
    }
    return result;
}

// --- back half #3: generalized shift-invert Arnoldi ------------------------
// Op is any shift-invert operator exposing apply() and the shared diagnostic
// accessors (generalized_shift_invert_operator / lu_shift_invert_operator).
template <typename _T, typename _Index, class Op>
static eig_result<_T> generalized_shift_invert_drive_(const spmats<_T,_Index>& self,
                                                        const spmats<_T,_Index>& B,
                                                        const std::size_t k,
                                                        const eig_options<_T>& options,
                                                        const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                                        const eig_solver_method actual_method,
                                                        const bool promoted_from_lanczos,
                                                        Op& op,
                                                        bool& pkg_converged,
                                                        bool& has_complex)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());

    struct ApplyFn {
        Op* op;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
    } apply_fn = { &op };

    // EIG-8 T-1 (e): μ 面 opt-in の λ 形式 lock ゲート(一般化 A x = λ B x)。
    // 受理は D8-2 の共有ヘルパ(c1_generalized_scale_ — B-48)による
    // r = ‖Ax − λBx‖ ≤ tol·max(1+|λ|, ‖A‖∞ + |λ|·‖B‖∞)。
    // A·x と B·x はゲート内で mv に 1:1 計上(B-24)+ 別建て(e-2)。
    const scalar_real_type anorm_gate_ks = matrix_inf_norm_value_<_T,_Index>(self);
    const scalar_real_type bnorm_gate_ks = matrix_inf_norm_value_<_T,_Index>(B);
    struct KsGenLambdaGate {
        const spmats<_T,_Index>* A;
        const spmats<_T,_Index>* Bm;
        scalar_real_type anorm;
        scalar_real_type bnorm;
        scalar_real_type tol;
        scalar_real_type sigma;
        bool operator()(const std::vector<_T>& x, const _T& mu, std::size_t& mv) const {
            const scalar_real_type mu_re = vcp::tsparse_scalar::real_part(mu);
            if (!(vcp::tsparse_scalar::abs_value(mu_re) > scalar_real_type(0))) return false;
            const scalar_real_type lam = scalar_real_type(1) / mu_re + sigma;
            const scalar_real_type r_abs =
                generalized_eigenpair_residual_norm_value_<_T,_Index>(*A, *Bm, _T(lam), x);
            mv += 2;   // λ 判定の A·x と B·x(B-24: 1:1 計上)
            const scalar_real_type scale = vcp::tsparse::c1_generalized_scale_(
                vcp::tsparse_scalar::abs_value(lam), anorm, bnorm);
            return r_abs <= tol * scale;   // certified-≤(interval は失敗側へ)
        }
    } ks_lambda_gate = { &self, &B, anorm_gate_ks, bnorm_gate_ks, options.tol, sigma };
    std::size_t ks_lambda_gate_count = 0;
    // EIG-3 T-3: KS core in mu space (old arnoldi core + max_iter*(k+1)
    // expansion deleted; options.max_iter = total inner-apply budget)
    const ks_mu_pkg_<_T> pkg =
        ks_mu_core_drive_<_T>(n, k, options, sigma, apply_fn,
                              ks_lambda_gate, true, &ks_lambda_gate_count);
    const std::size_t sdim = pkg.used_subspace_dim;

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
    result.lambda_gate_products = ks_lambda_gate_count;   // (e-2) 別建て計上
    result.linear_solves = op.linear_solves();
    result.inner_iterations = op.inner_iterations();
    result.inner_failure_count = op.inner_failure_count();
    result.inner_residual_norm = op.inner_residual_norm();
    result.factorization_diagnostics = op.factorization_diagnostics();
    result.factorization_zero_pivots = op.factorization_zero_pivots();
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

    // EIG-8 R4(裁定 s-2・必須併設): 返却対の複製方向監査(demote-only、
    // λ ゲート発火ラン限定 — 標準面 back half #2 と同一規約)。
    if (pkg.converged && ks_lambda_gate_count > 0 &&
        vcp::tsparse::returned_pair_duplicate_direction_found_<_T>(result.eigenvectors)) {
        if (result.failure_reason.empty())
            result.failure_reason = "duplicate direction among returned pairs "
                "(ghost-copy audit, EIG-8 R4): honest not_converged";
        pkg_converged = false;
    } else {
        pkg_converged = pkg.converged;
    }
    has_complex = pkg.has_complex;
    return result;
}

// --- E-A1 sparse_lu path: standard shift-invert Lanczos --------------------
template <typename _T, typename _Index>
static typename std::enable_if<std::is_signed<_Index>::value, eig_result<_T> >::type
shift_invert_lanczos_sparse_lu_(const spmats<_T,_Index>& self,
                                  const std::size_t k,
                                  const eig_options<_T>& options,
                                  const typename vcp::tsparse_scalar::real_type<_T>::type& sigma)
{
    typedef vcp::tsparse::lu_shift_invert_operator<spmats<_T,_Index> > LUOp;
    LUOp op(self, _T(sigma), options.shift_invert_lu);
    if (!op.factorization_ok()) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = eig_solver_method::shift_invert_lanczos;
        result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_lanczos);
        result.used_shift_invert = true;
        result.used_dense_fallback = false;
        result.status = "factorization_failed";
        result.failure_reason = "sparse LU factorization of (A - sigma*I) failed: "
            + op.factorization_diagnostics();
        result.message = result.failure_reason;
        result.factorization_diagnostics = op.factorization_diagnostics();
        result.factorization_zero_pivots = op.factorization_zero_pivots();
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    struct ApplyLU {
        LUOp* op;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
    } apply_lu = { &op };
    eig_result<_T> result = shift_invert_lanczos_drive_<_T,_Index>(self, k, options, sigma, apply_lu);
    // E-A1 E4 diagnostics: IR totals; IR non-convergence is a diagnostic
    // count only (no "inner_solve_failed" status on the sparse_lu path).
    result.linear_solves = op.linear_solves();
    result.inner_iterations = op.inner_iterations();
    result.inner_failure_count = op.inner_failure_count();
    result.inner_residual_norm = op.inner_residual_norm();
    result.factorization_diagnostics = op.factorization_diagnostics();
    result.factorization_zero_pivots = op.factorization_zero_pivots();
    // EIG-6 F-2': θ シフト磨き(converged だが scaled-abs 未クリーンの対が
    // ある場合のみ発火 — ヘルパ冒頭の設計コメント参照。E-A1 front 限定)。
    si_polish_shifted_rescue_<_T,_Index>(result, self, options);
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static typename std::enable_if<!std::is_signed<_Index>::value, eig_result<_T> >::type
shift_invert_lanczos_sparse_lu_(const spmats<_T,_Index>& self,
                                  const std::size_t k,
                                  const eig_options<_T>& options,
                                  const typename vcp::tsparse_scalar::real_type<_T>::type& sigma)
{
    (void)self; (void)options; (void)sigma;
    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::shift_invert_lanczos;
    result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_lanczos);
    result.used_shift_invert = true;
    result.used_dense_fallback = false;
    result.status = "factorization_failed";
    result.failure_reason = "shift-invert sparse_lu solver requires a signed Index type;"
        " set eig_options::shift_invert_solver = eig_shift_invert_solver::ilu0_gmres";
    result.message = result.failure_reason;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// --- E-A1 sparse_lu path: standard shift-invert Arnoldi --------------------
template <typename _T, typename _Index>
static typename std::enable_if<std::is_signed<_Index>::value, eig_result<_T> >::type
shift_invert_arnoldi_sparse_lu_(const spmats<_T,_Index>& self,
                                  const std::size_t k,
                                  const eig_options<_T>& options,
                                  const typename vcp::tsparse_scalar::real_type<_T>::type& sigma)
{
    typedef vcp::tsparse::lu_shift_invert_operator<spmats<_T,_Index> > LUOp;
    LUOp op(self, _T(sigma), options.shift_invert_lu);
    if (!op.factorization_ok()) {
        eig_result<_T> result;
        result.requested_count = k;
        result.method = eig_solver_method::shift_invert_arnoldi;
        result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_arnoldi);
        result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
        result.used_shift_invert = true;
        result.used_dense_fallback = false;
        result.status = "factorization_failed";
        result.failure_reason = "sparse LU factorization of (A - sigma*I) failed: "
            + op.factorization_diagnostics();
        result.message = result.failure_reason;
        result.factorization_diagnostics = op.factorization_diagnostics();
        result.factorization_zero_pivots = op.factorization_zero_pivots();
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    struct ApplyLU {
        LUOp* op;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
    } apply_lu = { &op };
    bool pkg_converged = false;
    eig_result<_T> result = shift_invert_arnoldi_drive_<_T,_Index>(self, k, options, sigma, apply_lu, pkg_converged);
    result.linear_solves = op.linear_solves();
    result.inner_iterations = op.inner_iterations();
    result.inner_failure_count = op.inner_failure_count();
    result.inner_residual_norm = op.inner_residual_norm();
    result.factorization_diagnostics = op.factorization_diagnostics();
    result.factorization_zero_pivots = op.factorization_zero_pivots();
    // E-A1 E4: direct solve; IR non-convergence does not gate convergence.
    result.converged = pkg_converged && result.eigenvalues.size() >= k;
    set_eig_diagnostics_<_T,_Index>(result, eig_solver_method::shift_invert_arnoldi, result.used_subspace_dim,
        result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
    // EIG-1 F-5 (G-2.1 案 (a)): arnoldi コア経路の converged 主張は C-1 判定後、
    // EIG-3 T-3: KS core exports C-2 evidence (D3-4) and the back-half
    // re-checked it in mu space; the EIG-1 arnoldi-core demotion is lifted.
    // Lambda-space C-1 remains the final acceptance gate.
    lambda_c1_acceptance_gate_<_T,_Index>(result, options.tol, self);   // EIG-4 T-4: 標準 si は改訂 scale
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    si_polish_rescue_standard_<_T,_Index>(result, self, options, apply_lu);
    result.linear_solves = op.linear_solves();
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static typename std::enable_if<!std::is_signed<_Index>::value, eig_result<_T> >::type
shift_invert_arnoldi_sparse_lu_(const spmats<_T,_Index>& self,
                                  const std::size_t k,
                                  const eig_options<_T>& options,
                                  const typename vcp::tsparse_scalar::real_type<_T>::type& sigma)
{
    (void)self; (void)sigma;
    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::shift_invert_arnoldi;
    result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::shift_invert_arnoldi);
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.used_shift_invert = true;
    result.used_dense_fallback = false;
    result.status = "factorization_failed";
    result.failure_reason = "shift-invert sparse_lu solver requires a signed Index type;"
        " set eig_options::shift_invert_solver = eig_shift_invert_solver::ilu0_gmres";
    result.message = result.failure_reason;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// --- E-A1 sparse_lu path: generalized shift-invert -------------------------
template <typename _T, typename _Index>
static typename std::enable_if<std::is_signed<_Index>::value, eig_result<_T> >::type
generalized_shift_invert_sparse_lu_(const spmats<_T,_Index>& self,
                                      const spmats<_T,_Index>& B,
                                      const std::size_t k,
                                      const eig_options<_T>& options,
                                      const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                      const eig_solver_method actual_method,
                                      const bool promoted_from_lanczos)
{
    typedef vcp::tsparse::lu_shift_invert_operator<spmats<_T,_Index> > LUOp;
    LUOp op(self, B, _T(sigma), options.shift_invert_lu);
    if (!op.factorization_ok()) {
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
        result.failure_reason = "sparse LU factorization of (A - sigma*B) failed: "
            + op.factorization_diagnostics();
        result.message = result.failure_reason;
        result.factorization_diagnostics = op.factorization_diagnostics();
        result.factorization_zero_pivots = op.factorization_zero_pivots();
        set_result_counts_<_T,_Index>(result, k);
        return result;
    }
    bool pkg_converged = false;
    bool has_complex = false;
    eig_result<_T> result = generalized_shift_invert_drive_<_T,_Index>(
        self, B, k, options, sigma, actual_method, promoted_from_lanczos,
        op, pkg_converged, has_complex);
    // E-A1 E4: direct solve; IR non-convergence does not gate convergence.
    result.converged = pkg_converged && !has_complex && (result.eigenvalues.size() >= k);
    if (has_complex) {
        result.status = "complex_ritz_values";
        result.message = "complex Ritz values detected in generalized shift-invert";
        if (result.failure_reason.empty())
            result.failure_reason = "complex Ritz values in requested subset";
    } else {
        set_eig_diagnostics_<_T,_Index>(result, actual_method, result.used_subspace_dim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
        if (promoted_from_lanczos)
            result.used_method = "shift_invert_arnoldi(promoted_from_lanczos)";
    }
    // EIG-3 T-3: KS core exports C-2 evidence (D3-4) and the back-half
    // re-checked it in mu space; the EIG-1 arnoldi-core demotion is lifted.
    // Lambda-space C-1 remains the final acceptance gate.
    // EIG-8 T-3 (D8-3): 一般化は D8-2 共有 scale 版(挙動遷移期待ゼロ —
    // 既存 rel 分岐が Frobenius 後退正規化を内包。契約整合の完成)。
    lambda_c1_acceptance_gate_<_T,_Index>(result, options.tol, self, B);
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    {
        struct ApplyPolish {
            LUOp* op;
            void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
        } apply_polish = { &op };
        si_polish_rescue_generalized_<_T,_Index>(result, self, B, options, apply_polish);
        result.linear_solves = op.linear_solves();
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static typename std::enable_if<!std::is_signed<_Index>::value, eig_result<_T> >::type
generalized_shift_invert_sparse_lu_(const spmats<_T,_Index>& self,
                                      const spmats<_T,_Index>& B,
                                      const std::size_t k,
                                      const eig_options<_T>& options,
                                      const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                      const eig_solver_method actual_method,
                                      const bool promoted_from_lanczos)
{
    (void)self; (void)B; (void)options; (void)sigma;
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
    result.failure_reason = "shift-invert sparse_lu solver requires a signed Index type;"
        " set eig_options::shift_invert_solver = eig_shift_invert_solver::ilu0_gmres";
    result.message = result.failure_reason;
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
    // E-A1 D-1: solver dispatch.  Default = sparse_lu direct solve; the
    // legacy ILU(0)+GMRES implementation below is the ilu0_gmres opt-in
    // compatibility path and is unchanged (B-8).
    if (options.shift_invert_solver == eig_shift_invert_solver::sparse_lu) {
        return shift_invert_lanczos_sparse_lu_<_T,_Index>(self, k, options, sigma);
    }
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

    // E-A1: back half (driver -> lambda inversion -> result assembly)
    // mechanically extracted to shift_invert_lanczos_drive_; the legacy-path
    // computation and results are unchanged (shared with the sparse_lu path).
    eig_result<_T> result = shift_invert_lanczos_drive_<_T,_Index>(self, k, options, sigma, apply_si);
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

    // EIG-3 T-4 (D-6, Q1 ruling 2026-07-06): max_iter is the TOTAL mv budget on
    // every path.  Pure unit conversion to restarts with the worst per-restart
    // cost 2*sdim (expansion sdim + lock-scan exact residuals <= sdim;
    // tsparse_lanczos.hpp mv_count++ sites); the inflation term "+ k + 1" is
    // deleted.  The restart loop is INCLUSIVE (restart <= max_restarts), so a
    // floor of 0 already yields one restart; no max(1,.) floor (it would
    // double the minimum work and break the bound for max_iter < 2*sdim).
    // Machine-checked overshoot bound: mv <= max_iter + 2*sdim.
    const std::size_t max_restarts_l = (sdim > 0)
        ? options.max_iter / (2 * sdim)
        : options.max_iter;
    auto pkg = vcp::tsparse_lanczos::lanczos_eigs_standard<_T, apply_fn>(
        n, k, sdim, max_restarts_l, options.tol,
        options.random_seed, options.random_start,
        options.target, shift_val, options.compute_residual_history, apply);

    return lanczos_package_to_result_<_T,_Index>(pkg, self, k, eig_solver_method::lanczos, options.tol);
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
    // E-A1 D-1: solver dispatch.  Default = sparse_lu direct solve; the
    // legacy ILU(0)+GMRES implementation below is the ilu0_gmres opt-in
    // compatibility path and is unchanged (B-8).
    if (options.shift_invert_solver == eig_shift_invert_solver::sparse_lu) {
        return shift_invert_arnoldi_sparse_lu_<_T,_Index>(self, k, options, sigma);
    }
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
    // E-A1: back half (driver -> lambda inversion -> result assembly ->
    // residuals) mechanically extracted to shift_invert_arnoldi_drive_; the
    // legacy-path computation and results are unchanged (shared with the
    // sparse_lu path).  sdim is available as result.used_subspace_dim.
    bool pkg_converged = false;
    eig_result<_T> result = shift_invert_arnoldi_drive_<_T,_Index>(self, k, options, sigma, apply_si, pkg_converged);
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_failure_count = inner_failure_count;
    result.inner_residual_norm = inner_residual_norm;
    result.factorization_diagnostics = ilu.diagnostics;
    result.factorization_zero_pivots = ilu.zero_pivots;
    result.converged = pkg_converged && result.eigenvalues.size() >= k && inner_failure_count == 0;
    if (inner_failure_count != 0) {
        result.status = "inner_solve_failed";
        result.failure_reason = "inner GMRES solve failed";
        result.inner_failure_reason = result.failure_reason;
        result.message = result.failure_reason;
    } else {
        set_eig_diagnostics_<_T,_Index>(result, eig_solver_method::shift_invert_arnoldi, result.used_subspace_dim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
    }
    // EIG-1 F-5 (G-2.1 案 (a)): arnoldi コア経路の降格(legacy ilu0_gmres 分岐も
    // EIG-3 T-3: KS core exports C-2 evidence (D3-4) and the back-half
    // re-checked it in mu space; the EIG-1 arnoldi-core demotion is lifted.
    // Lambda-space C-1 remains the final acceptance gate.
    lambda_c1_acceptance_gate_<_T,_Index>(result, options.tol, self);   // EIG-4 T-4: 標準 si は改訂 scale
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    si_polish_rescue_standard_<_T,_Index>(result, self, options, apply_si);
    result.linear_solves = linear_solve_count;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// arnoldi_eigs_new_
//
// EIG-3 T-2 (D3-1): the public `arnoldi` enum now dispatches to the rebuilt
// Krylov-Schur driver (tsparse_krylov_schur.hpp, EIG-2 real Schur core).
// The old Ritz-restart arnoldi implementation (zero convergence record,
// D-3 value dedup, D-5 unbudgeted restarts) is replaced, not preserved.
// Budget: options.max_iter is passed through as the TOTAL matrix-vector
// product budget (D-6 unification; the old max_iter*(k+1) restart expansion
// is deleted).  Declared behavior change: timeouts/lies -> OK or honest
// not_converged.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> arnoldi_ks_drive_(const spmats<_T,_Index>& /*self*/,
                                        const std::size_t k,
                                        const eig_options<_T>& /*options*/,
                                        std::true_type /* is_complex */)
{
    // Fence for complex instantiations (runtime dispatch returns earlier via
    // complex_standard_eigs_with_info_; this body is never executed).
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    result.failure_reason = "arnoldi(krylov_schur): complex scalar unsupported";
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> arnoldi_ks_drive_(const spmats<_T,_Index>& self,
                                        const std::size_t k,
                                        const eig_options<_T>& options,
                                        std::false_type /* is_complex */)
{
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    };
    apply_fn apply_op = { &A };
    vcp::tsparse_experimental::krylov_schur_result<_T> d =
        vcp::tsparse_experimental::krylov_schur_eigs_with_diagnostics<apply_fn, _T>(
            apply_op, n, k, options);
    eig_result<_T> result = d.eigs;
    result.method = eig_solver_method::arnoldi;
    result.used_method = "arnoldi(krylov_schur)";
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> arnoldi_eigs_new_(const spmats<_T,_Index>& self,
                                          const std::size_t k,
                                          const eig_options<_T>& options)
{
    if (options.method == eig_solver_method::shift_invert_arnoldi) {
        return shift_invert_arnoldi_eigs_<_T,_Index>(self, k, options, options.shift);
    }
    return arnoldi_ks_drive_<_T,_Index>(self, k, options,
        typename std::integral_constant<bool, spmatrix_is_complex<_T>::value>::type());
}

// ---------------------------------------------------------------------------
// EIG-4 T-1: explicit thick_restart_lanczos / krylov_schur wrappers.
//
// Low-level explicit selection (D4-1): the TRL wrapper does NOT reject
// nonsymmetric input up front (same convention as the v1 `trl` cases) --
// misuse is demoted to an honest failure by the TRL end-of-run exact C-1
// check (D-12 line).  The driver-internal used_method strings
// ("thick_restart_lanczos_experimental" / "krylov_schur") are left unchanged
// for direct-driver callers; only the dispatch-level result is renamed.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> trl_drive_(const spmats<_T,_Index>& /*self*/,
                                 const std::size_t k,
                                 const eig_options<_T>& /*options*/,
                                 std::false_type /* trl_supported */)
{
    // Fence for scalar types the TRL driver cannot instantiate (complex and
    // non-builtin real scalars: kv::dd / kv::mpfr / kv::interval -- the
    // driver static_asserts trl_is_real_floating).  Honest explicit refusal
    // (C-5); the auto_select routing never reaches this (it sends such
    // scalars to the KS side), only explicit method=thick_restart_lanczos.
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    result.failure_reason =
        "thick_restart_lanczos: scalar type unsupported (builtin real"
        " floating-point only); use krylov_schur or lanczos";
    result.message = result.failure_reason;
    result.method = eig_solver_method::thick_restart_lanczos;
    result.used_method = "thick_restart_lanczos";
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> trl_drive_(const spmats<_T,_Index>& self,
                                 const std::size_t k,
                                 const eig_options<_T>& options,
                                 std::true_type /* trl_supported */)
{
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    };
    apply_fn apply_op = { &A };
    eig_result<_T> result =
        vcp::tsparse_experimental::thick_restart_lanczos_eigs<apply_fn, _T>(
            apply_op, n, k, options);
    result.method = eig_solver_method::thick_restart_lanczos;
    result.used_method = "thick_restart_lanczos";
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> trl_eigs_new_(const spmats<_T,_Index>& self,
                                    const std::size_t k,
                                    const eig_options<_T>& options)
{
    return trl_drive_<_T,_Index>(self, k, options,
        typename std::integral_constant<bool,
            vcp::tsparse_experimental::trl_is_real_floating<_T>::value>::type());
}

template <typename _T, typename _Index>
static eig_result<_T> ks_explicit_eigs_(const spmats<_T,_Index>& self,
                                        const std::size_t k,
                                        const eig_options<_T>& options)
{
    eig_result<_T> result = arnoldi_ks_drive_<_T,_Index>(self, k, options,
        typename std::integral_constant<bool, spmatrix_is_complex<_T>::value>::type());
    result.method = eig_solver_method::krylov_schur;
    result.used_method = "krylov_schur";
    return result;
}

// ---------------------------------------------------------------------------
// EIG-4 T-3: dense return path with complex-pair opt-in (real scalars only;
// reached only when options.allow_complex_pairs == true -- B-27).
// Selection is pair-aware slot filling (keep-together: a straddled pair
// returns k+1 values -- R-6).  converged carries the full-spectrum verdict of
// real_schur_eig_dense_pairs (legacy dense acceptance scale).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> dense_pairs_eig_impl_(const std::vector<std::vector<_T> >& /*dense*/,
                                            const std::size_t k,
                                            const eig_options<_T>& /*options*/,
                                            std::true_type /* is_complex */)
{
    // Fence: complex scalars never reach the pairs path (runtime guarded).
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    result.failure_reason = "allow_complex_pairs: complex scalar unsupported";
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> dense_pairs_eig_impl_(const std::vector<std::vector<_T> >& dense,
                                            const std::size_t k,
                                            const eig_options<_T>& options,
                                            std::false_type /* is_complex */)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    std::string reason;
    vcp::tsparse_dense_schur::real_schur_pairs_output<_T> po;
    vcp::tsparse_dense_linalg::dense_eigen_result<_T> base =
        vcp::tsparse_dense_schur::real_schur_eig_dense_pairs(dense, options.tol, reason, po);

    eig_result<_T> result;
    result.requested_count = k;
    result.iterations = base.iterations;
    result.complex_eigenvalues = base.complex_eigenvalues;
    result.residual_norm_absolute = base.residual_norm;

    const std::size_t nreal = base.eigenvalues.size();
    std::vector<scalar_real_type> re, im;
    for (std::size_t i = 0; i < nreal; i++) {
        re.push_back(vcp::tsparse_scalar::real_part(base.eigenvalues[i]));
        im.push_back(scalar_real_type(0));
    }
    for (std::size_t j = 0; j < po.pair_re.size(); j++) {
        re.push_back(vcp::tsparse_scalar::real_part(po.pair_re[j]));
        im.push_back(vcp::tsparse_scalar::real_part(po.pair_im[j]));
    }
    const std::vector<std::size_t> sel =
        vcp::tsparse_experimental::ks_detail::select_slot_indices_pairs<scalar_real_type>(
            re, im, k, options.target,
            vcp::tsparse_scalar::real_part(options.shift));
    std::size_t pair_count = 0;
    for (std::size_t s = 0; s < sel.size(); s++) {
        const std::size_t idx = sel[s];
        if (idx < nreal) {
            result.eigenvalues.push_back(base.eigenvalues[idx]);
            result.eigenvalues_imag.push_back(_T(0));
            if (idx < base.eigenvectors.size())
                result.eigenvectors.push_back(base.eigenvectors[idx]);
            if (idx < base.residuals.size()) {
                result.residuals_absolute.push_back(base.residuals[idx]);
                result.residuals_relative.push_back(base.residuals[idx]
                    / (scalar_real_type(1)
                       + vcp::tsparse_scalar::abs_value(base.eigenvalues[idx])));
            }
        } else {
            const std::size_t j = idx - nreal;
            result.eigenvalues.push_back(po.pair_re[j]);
            result.eigenvalues.push_back(po.pair_re[j]);
            result.eigenvalues_imag.push_back(po.pair_im[j]);
            result.eigenvalues_imag.push_back(-po.pair_im[j]);
            result.eigenvectors.push_back(po.pair_u[j]);
            result.eigenvectors.push_back(po.pair_v[j]);
            const scalar_real_type mag = vcp::tsparse_scalar::sqrt_value(
                vcp::tsparse_scalar::real_part(po.pair_re[j]) * vcp::tsparse_scalar::real_part(po.pair_re[j])
              + vcp::tsparse_scalar::real_part(po.pair_im[j]) * vcp::tsparse_scalar::real_part(po.pair_im[j]));
            result.residuals_absolute.push_back(po.pair_res[j]);
            result.residuals_absolute.push_back(po.pair_res[j]);
            result.residuals_relative.push_back(po.pair_res[j] / (scalar_real_type(1) + mag));
            result.residuals_relative.push_back(po.pair_res[j] / (scalar_real_type(1) + mag));
            pair_count++;
        }
    }
    if (pair_count == 0) result.eigenvalues_imag.clear();
    result.complex_pair_count = pair_count;

    result.converged = base.converged && !sel.empty();
    if (result.converged) {
        result.status = "converged";
        result.message = "converged";
    } else {
        result.status = "not_converged";
        result.failure_reason = reason.empty()
            ? "dense (complex-pair opt-in) did not satisfy the acceptance test"
            : reason;
        result.message = result.failure_reason;
    }
    result.method = eig_solver_method::dense_fallback_explicit;
    result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::dense_fallback_explicit);
    result.used_dense_fallback = true;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> dense_pairs_eig_(const std::vector<std::vector<_T> >& dense,
                                       const std::size_t k,
                                       const eig_options<_T>& options)
{
    return dense_pairs_eig_impl_<_T,_Index>(dense, k, options,
        typename std::integral_constant<bool, spmatrix_is_complex<_T>::value>::type());
}

// ---------------------------------------------------------------------------
// EIG-4 T-2: auto_select dense branch.
//
// The jacobi / real-Schur split uses the SAME certainly check as the auto
// symmetric routing (B-28; approved R-5) -- NOT the tol-fuzzy
// is_dense_symmetric of the explicit dense_fallback_explicit path (which is
// untouched, B-26).  The symmetric branch applies the jacobi rotation-budget
// floor max(max_iter, eig_auto_jacobi_budget_factor*n^2) (G-1.1 c-1): the
// dense route is mv-free (outside B-1), and the floor realizes the D4-2
// intent that the dense region answers deterministically; adversarial
// spectra beyond the floor end in an HONEST not_converged.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> auto_dense_eigs_(const spmats<_T,_Index>& A,
                                       const std::size_t k,
                                       const eig_options<_T>& options,
                                       const bool certainly_sym)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(A.rowsize());
    std::vector<std::vector<_T> > dense = to_dense_impl_<_T,_Index>(A);
    eig_result<_T> result;
    if (certainly_sym) {
        const std::size_t floor_budget = eig_auto_jacobi_budget_factor * n * n;
        const std::size_t budget = (options.max_iter > floor_budget)
            ? options.max_iter : floor_budget;
        result = convert_dense_result_<_T,_Index>(
            vcp::tsparse_dense_linalg::jacobi_eig_dense(dense, budget, options.tol),
            eig_solver_method::dense_fallback_explicit);
    } else {
        // EIG-4 T-3 (additive, B-27): opt-in complex pairs (real scalars)
        if (options.allow_complex_pairs && !spmatrix_is_complex<_T>::value) {
            return dense_pairs_eig_<_T,_Index>(dense, k, options);
        }
        result = dense_nonsymmetric_eig_<_T,_Index>(dense, options,
            typename std::integral_constant<bool, spmatrix_is_complex<_T>::value>::type());
    }
    result.method = eig_solver_method::dense_fallback_explicit;
    result.used_method = eig_method_to_string_<_T,_Index>(eig_solver_method::dense_fallback_explicit);
    result.used_dense_fallback = true;
    // aligned selection (9b): the auto default must return value/vector pairs
    // that actually correspond (the legacy helper's guard skips the vector
    // trim -- kept byte-identical on explicit paths per B-26, issue filed)
    select_eigenpairs_aligned_<_T,_Index>(result, k, options.target, options.shift);
    // EIG-4: exact end-of-run residuals of the RETURNED pairs against the
    // sparse A (C-1 spirit on the new default path; the jacobi/Schur internal
    // verdicts are additionally gated below -- demotion only, never promotion)
    if (!result.eigenvectors.empty()
        && result.eigenvectors.size() == result.eigenvalues.size()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        if (result.converged) {
            std::vector<scalar_real_type> theta_abs_c1;
            theta_abs_c1.reserve(result.eigenvalues.size());
            for (std::size_t i = 0; i < result.eigenvalues.size(); i++)
                theta_abs_c1.push_back(vcp::tsparse_scalar::abs_value(
                    vcp::tsparse_scalar::real_part(result.eigenvalues[i])));
            // dense の既存受理スケール(tol·10n)を含めて広義に判定する
            scalar_real_type anorm = matrix_inf_norm_value_<_T,_Index>(A);
            const scalar_real_type tol10n =
                options.tol * scalar_real_type(static_cast<int>(n) * 10);
            if (tol10n > options.tol * anorm) anorm = tol10n / options.tol;
            if (!vcp::tsparse::residual_acceptance_check_scaled_(
                    result.residuals_absolute, result.residuals_relative,
                    options.tol, theta_abs_c1, anorm)) {
                result.converged = false;
                result.status = "residual_check_failed";
                result.failure_reason =
                    "end-of-run exact residual failed acceptance (C-1, auto dense)";
                result.message = result.failure_reason;
            }
        }
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// EIG-4 T-2 (D4-2): auto_select routing.  ADDITIVE dispatch only (B-26):
// every branch re-enters an existing path or a new wrapper; no pre-existing
// branch is rewritten.  used_method wraps the routed diagnostic as
// "auto_select(<routed>)" for mechanical routing verification (v2 suite).
//
//   1. use_shift == true                  -> legacy default shift path
//      (method := shift_invert_lanczos, E-A1 LU; identical to the old
//      lanczos-default + use_shift behavior)
//   2. n <= eig_auto_dense_threshold and dense conversion permitted
//                                          -> dense (jacobi / real Schur by
//                                             the B-28 certainly check)
//   3. certainly symmetric (hint or check) -> thick_restart_lanczos
//   4. otherwise (incl. uncertifiable)     -> krylov_schur (B-28: KS is
//                                             correct for symmetric input too)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> auto_select_eigs_(const spmats<_T,_Index>& A,
                                        const std::size_t k,
                                        const eig_options<_T>& options)
{
    const std::size_t n = static_cast<std::size_t>(A.rowsize());

    if (options.use_shift) {
        eig_options<_T> active = options;
        active.method = eig_solver_method::shift_invert_lanczos;
        eig_result<_T> result = lanczos_eigs_new_<_T,_Index>(A, k, active);
        result.method = eig_solver_method::auto_select;
        result.used_method = "auto_select(" + result.used_method + ")";
        return result;
    }

    if (n <= eig_auto_dense_threshold
        && options.allow_dense_conversion
        && n * n <= options.max_dense_size) {
        const bool sym = (options.structure == matrix_structure_hint::symmetric)
            || (options.structure == matrix_structure_hint::auto_detect
                && is_certainly_symmetric_<_T,_Index>(A));
        eig_result<_T> result = auto_dense_eigs_<_T,_Index>(A, k, options, sym);
        result.method = eig_solver_method::auto_select;
        result.used_method = "auto_select(" + result.used_method + ")";
        return result;
    }

    bool sym;
    if (options.structure == matrix_structure_hint::symmetric)     sym = true;
    else if (options.structure == matrix_structure_hint::general)  sym = false;
    else sym = is_certainly_symmetric_<_T,_Index>(A);

    // B-28 と同方向の保守則: TRL が対応しないスカラー型(kv::dd / kv::mpfr /
    // kv::interval — 組込み浮動小数点以外)は対称でも KS 側に倒す(KS は
    // 対称入力でも正しく、汎用スカラーで動く。auto が型を理由に拒否結果へ
    // ルートすることはない)。
    const bool trl_ok =
        vcp::tsparse_experimental::trl_is_real_floating<_T>::value;

    eig_result<_T> result = (sym && trl_ok)
        ? trl_eigs_new_<_T,_Index>(A, k, options)
        : ks_explicit_eigs_<_T,_Index>(A, k, options);
    result.method = eig_solver_method::auto_select;
    result.used_method = "auto_select(" + result.used_method + ")";
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

    // E-A1 D-1: solver dispatch.  Default = sparse_lu direct solve; the
    // legacy ILU(0)+GMRES implementation below is the ilu0_gmres opt-in
    // compatibility path and is unchanged (B-8; the inner_max_iter /
    // inner_tol / inner_restart computations moved INTO the legacy branch
    // with values and expressions unchanged, because the sparse_lu path
    // does not read them -- D-2).
    if (options.shift_invert_solver == eig_shift_invert_solver::sparse_lu) {
        return generalized_shift_invert_sparse_lu_<_T,_Index>(
            self, B, k, options, sigma, actual_method, promoted_from_lanczos);
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

    // E-A1: back half (apply wrapper -> driver -> lambda inversion -> result
    // assembly -> residuals) mechanically extracted to
    // generalized_shift_invert_drive_; the legacy-path computation and
    // results are unchanged (shared with the sparse_lu path).
    bool pkg_converged = false;
    bool has_complex = false;
    eig_result<_T> result = generalized_shift_invert_drive_<_T,_Index>(
        self, B, k, options, sigma, actual_method, promoted_from_lanczos,
        gsi_op, pkg_converged, has_complex);

    const bool inner_ok = (gsi_op.inner_failure_count() == 0);
    result.converged = pkg_converged && !has_complex && (result.eigenvalues.size() >= k) && inner_ok;

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
        set_eig_diagnostics_<_T,_Index>(result, actual_method, result.used_subspace_dim,
            result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
        if (promoted_from_lanczos)
            result.used_method = "shift_invert_arnoldi(promoted_from_lanczos)";
    }
    // EIG-3 T-3: KS core exports C-2 evidence (D3-4) and the back-half
    // re-checked it in mu space; the EIG-1 arnoldi-core demotion is lifted.
    // Lambda-space C-1 remains the final acceptance gate.
    // EIG-8 T-3 (D8-3): 一般化は D8-2 共有 scale 版(挙動遷移期待ゼロ —
    // 既存 rel 分岐が Frobenius 後退正規化を内包。契約整合の完成)。
    lambda_c1_acceptance_gate_<_T,_Index>(result, options.tol, self, B);
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    {
        struct ApplyPolish {
            GSIOperator* op;
            void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
        } apply_polish = { &gsi_op };
        si_polish_rescue_generalized_<_T,_Index>(result, self, B, options, apply_polish);
        result.linear_solves = gsi_op.linear_solves();
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
    // EIG-3 T-4 (D-6, Q1 ruling 2026-07-06): max_iter is the TOTAL mv budget on
    // every path.  Pure unit conversion to restarts with the worst per-restart
    // cost 2*sdim (expansion sdim + lock-scan exact residuals <= sdim;
    // tsparse_lanczos.hpp mv_count++ sites); the inflation term "+ k + 1" is
    // deleted.  The restart loop is INCLUSIVE (restart <= max_restarts), so a
    // floor of 0 already yields one restart; no max(1,.) floor (it would
    // double the minimum work and break the bound for max_iter < 2*sdim).
    // Machine-checked overshoot bound: mv <= max_iter + 2*sdim.
    const std::size_t max_restarts_si = (sdim > 0)
        ? options.max_iter / (2 * sdim)
        : options.max_iter;

    // EIG-6 F-2': with_prec si_lanczos も同一の opt-in(drive_ と同じ λ ゲート
    // 構成 — 共有ヘルパ e-1 / 1:1 計上 e-2 / 最終ゲート再検査 e-3)。
    const scalar_real_type anorm_gate_p = matrix_inf_norm_value_<_T,_Index>(A);
    struct SiLambdaGateP {
        const spmats<_T,_Index>* Am;
        scalar_real_type anorm;
        scalar_real_type tol;
        scalar_real_type sigma;
        bool operator()(const std::vector<_T>& x, const _T& mu, std::size_t& mv) const {
            const scalar_real_type mu_re = vcp::tsparse_scalar::real_part(mu);
            if (!(vcp::tsparse_scalar::abs_value(mu_re) > scalar_real_type(0))) return false;
            const scalar_real_type lam = scalar_real_type(1) / mu_re + sigma;
            const scalar_real_type r_abs =
                eigenpair_residual_norm_value_<_T,_Index>(*Am, _T(lam), x);
            mv += 1;
            const scalar_real_type scale = vcp::tsparse::c1_revised_scale_(
                vcp::tsparse_scalar::abs_value(lam), anorm);
            return r_abs <= tol * scale;
        }
    } si_lambda_gate_p = { &A, anorm_gate_p, options.tol, sigma };
    std::size_t lambda_gate_count_p = 0;
    typedef vcp::tsparse_lanczos::lanczos_result_package<_T, SIApply> LPkg;
    LPkg pkg = vcp::tsparse_lanczos::lanczos_eigs_standard_si_<_T, SIApply, SiLambdaGateP>(
        n, k, sdim, max_restarts_si, options.tol,
        options.random_seed, options.random_start,
        eig_target::largest_magnitude, scalar_real_type(0),
        options.compute_residual_history, apply_si, si_lambda_gate_p, &lambda_gate_count_p);

    for (std::size_t i = 0; i < pkg.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(pkg.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            pkg.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }

    eig_result<_T> result = lanczos_package_to_result_<_T,_Index>(pkg, self, k, eig_solver_method::shift_invert_lanczos, options.tol);
    result.lambda_gate_products = lambda_gate_count_p;   // (e-2)
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
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    si_polish_rescue_standard_<_T,_Index>(result, self, options, apply_si);
    // -----------------------------------------------------------------------
    // EIG-10.1 D-17d-wp 委譲(G-0.1/G-1.1 承認・EIG-10 D-17d ブロックの
    // 逐語鏡映・additive route)。
    //
    // 発火署名は EIG-10 と同一(既存フィールドの読み取りのみ・新定数ゼロ):
    //   !converged ∧ status == residual_check_failed ∧ 残予算 > 0。
    //   inner_solve_failed / preconditioner_failed は署名不一致 = 構造休眠。
    //   凍結スイートに preconditioner オーバーロード呼び出しは存在しない
    //   (EIG-10.1 G-0.1 §4)ため既存挙動はビット不変。
    //
    // 委譲先は plain 面と同じ shift_invert_arnoldi_drive_(same-helper)に
    // **wp 既存の apply_si(内側 GMRES+前処理 closure)をそのまま渡す**。
    // 前処理オブジェクトの引き回しは closure 経由で暗黙に実現され、委譲先の
    // 意味論は不変(G-0.1 承認 2)。plain 面との正直な差異: 内側解が実質
    // 厳密な族(diag 族)のみ委譲が採用に至り、mock/rpp6c 級は内側 GMRES
    // 精度により C-1 正直棄却 → all-or-nothing 不採用(値・status 不変・
    // 消費のみ計上 = EIG-10 mixed/interval と同じ受理済みクラス)。
    //
    // 規律((c-1) 全継承):
    //  - 1 回限り・決定的(乱数なし・ループ禁止)。
    //  - all-or-nothing: 委譲結果が「全ゲート合格の converged ∧ 返却全対実 ∧
    //    複製方向なし(R4 と同一ヘルパ)」の場合のみ丸ごと採用。それ以外は
    //    元の正直結果をそのまま返す(対集合の混合・部分採用禁止)。
    //  - allow_complex_pairs=true の内部指定は「過渡複素窓の通過許可」であって
    //    「複素対の返却許可」ではない(全実ガードにより structure=symmetric
    //    呼び出しへ複素対が漏れる経路は構造的に不存在)。
    //  - mv は採否に依らず実消費を 1:1 計上(B-38)し ks_rescue_products に
    //    別建て記録。KS コアは mv ≤ 残予算を厳守するため mv_total ≤ max_iter
    //    (B-1)。
    //  - wp 会計 (i)(G-0.1 裁定・恒等式維持の必須要件): ls /
    //    inner_iterations / inner_residual_norm / inner_failure_count は
    //    ブロック末尾で共有カウンタから 1:1 再同期(委譲・磨きの内側消費を
    //    含む物理量)。
    //  - wp 会計 (ii)(G-0.1 裁定): 不採用時は status=rcf のまま消費のみ
    //    反映。診断組合せ「status==residual_check_failed ∧
    //    inner_failure_count > 0」は「委譲を試みたが内側精度で正直に不採用」
    //    の意味(嘘ではない — 採用は C-1 再ゲート合格のみが通す)。
    // -----------------------------------------------------------------------
    if (!result.converged
        && result.status == "residual_check_failed"
        && options.max_iter > result.matrix_vector_products) {
        eig_options<_T> ks_opts = options;
        ks_opts.max_iter = options.max_iter - result.matrix_vector_products;
        ks_opts.allow_complex_pairs = true;   // 過渡複素窓の通過許可(上記規律)
        bool ks_conv = false;
        eig_result<_T> alt = shift_invert_arnoldi_drive_<_T,_Index>(
            self, k, ks_opts, sigma, apply_si, ks_conv);
        alt.converged = ks_conv && alt.eigenvalues.size() >= k;
        // 帰結配線は EIG-10 ブロックと同一(B-43 再利用): C-1 改訂 scale
        // 最終ゲート → F-1 条件付き磨き(不合格時のみ発火・same-helper。
        // 磨き solve は apply_si 経由 = 内側カウンタへ自動計上、mv は helper
        // 内で 1:1 計上され下の rescue 合算に含まれる)。
        lambda_c1_acceptance_gate_<_T,_Index>(alt, options.tol, self);
        si_polish_rescue_standard_<_T,_Index>(alt, self, ks_opts, apply_si);
        result.ks_rescue_products = alt.matrix_vector_products;
        result.matrix_vector_products += alt.matrix_vector_products;
        result.lambda_gate_products += alt.lambda_gate_products;
        if (alt.converged && si_polish_returned_pairs_all_real_<_T,_Index>(alt)) {
            // 採用前の基底整形: 返却集合の MGS 直交正規化(2 パス)。縮退固有
            // 空間の基底の向きを整えるだけで、対集合の混合・部分採用ではない
            // ((c-1) の all-or-nothing は集合単位で維持)。整形後に残差を
            // 再計算し、最終ゲートと同一の受理式で全対を再判定してから採用する
            // (整形で品質が立たない場合・基底が退化する場合は丸ごと不採用)。
            typedef typename vcp::tsparse_scalar::real_type<_T>::type R__;
            std::vector<std::vector<_T> > ortho = alt.eigenvectors;
            bool ortho_ok = (ortho.size() == alt.eigenvalues.size()) && !ortho.empty();
            for (std::size_t i = 0; ortho_ok && i < ortho.size(); i++) {
                for (int pass = 0; pass < 2; pass++) {
                    for (std::size_t j = 0; j < i; j++) {
                        const R__ c = vcp::tsparse_scalar::real_dot_value(ortho[j], ortho[i]);
                        for (std::size_t t = 0; t < ortho[i].size(); t++)
                            ortho[i][t] -= _T(c) * ortho[j][t];
                    }
                }
                const R__ nv = vcp::tsparse_scalar::real_norm_value(ortho[i]);
                if (!(nv > R__(0))) { ortho_ok = false; break; }   // certified > 0 のみ採用側
                for (std::size_t t = 0; t < ortho[i].size(); t++)
                    ortho[i][t] = ortho[i][t] / _T(nv);
            }
            if (ortho_ok) {
                std::vector<R__> res_abs_o = eigenpair_residuals_<_T,_Index>(
                    A, alt.eigenvalues, ortho);
                std::vector<R__> res_rel_o = eigenpair_relative_residuals_<_T,_Index>(
                    A, alt.eigenvalues, ortho);
                std::vector<R__> theta_abs_o;
                theta_abs_o.reserve(alt.eigenvalues.size());
                for (std::size_t i = 0; i < alt.eigenvalues.size(); i++)
                    theta_abs_o.push_back(vcp::tsparse_scalar::abs_value(
                        vcp::tsparse_scalar::real_part(alt.eigenvalues[i])));
                if (vcp::tsparse::residual_acceptance_check_scaled_(
                        res_abs_o, res_rel_o, options.tol, theta_abs_o, anorm_gate_p)
                    && !vcp::tsparse::returned_pair_duplicate_direction_found_<_T>(ortho)) {
                    result.eigenvalues          = alt.eigenvalues;
                    result.eigenvalues_imag     = alt.eigenvalues_imag;
                    result.complex_pair_count   = alt.complex_pair_count;
                    result.complex_eigenvalues  = alt.complex_eigenvalues;
                    result.eigenvectors         = ortho;
                    result.residuals_absolute   = res_abs_o;
                    result.residuals_relative   = res_rel_o;
                    if (!res_abs_o.empty()) {
                        result.residual_norm_absolute = *std::max_element(
                            res_abs_o.begin(), res_abs_o.end());
                    }
                    result.converged        = true;
                    result.converged_count  = alt.converged_count;
                    result.status           = "converged";
                    result.failure_reason.clear();
                    result.breakdown_reason.clear();
                    result.message = "converged (D-17d-wp rescue: ks_mu_core delegation; EIG-10.1)";
                    result.iterations += alt.iterations;
                    // method / used_method / used_subspace_dim は公開面
                    // ("shift_invert_lanczos+user_preconditioner")のまま
                    // (EIG-10 G-0.1 §6.1 と同じ: enum・route 文字列追加なし)。
                    populate_real_complex_eigenvalues_<_T,_Index>(result);
                    set_result_counts_<_T,_Index>(result, k);
                }
            }
        }
        // 不採用: 元の正直結果を維持(値・status 不変。消費計上のみ反映済み)。
    }
    // wp 会計 (i): 内側消費の 1:1 再同期(委譲・磨き分を含む。恒等式
    // mv == ls + lgate の維持に必須 — G-0.1 裁定で凍結)。
    result.linear_solves = linear_solve_count;
    result.inner_iterations = inner_iteration_count;
    result.inner_residual_norm = inner_residual_norm;
    result.inner_failure_count = inner_failure_count;
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

    // EIG-8 T-1 (e): μ 面 opt-in の λ 形式 lock ゲート(標準問題。
    // shift_invert_arnoldi_drive_ と同一式 — 共有ヘルパ e-1、mv 1:1 e-2)。
    const scalar_real_type anorm_gate_ks = matrix_inf_norm_value_<_T,_Index>(A);
    struct KsLambdaGateP {
        const spmats<_T,_Index>* A;
        scalar_real_type anorm;
        scalar_real_type tol;
        scalar_real_type sigma;
        bool operator()(const std::vector<_T>& x, const _T& mu, std::size_t& mv) const {
            const scalar_real_type mu_re = vcp::tsparse_scalar::real_part(mu);
            if (!(vcp::tsparse_scalar::abs_value(mu_re) > scalar_real_type(0))) return false;
            const scalar_real_type lam = scalar_real_type(1) / mu_re + sigma;
            const scalar_real_type r_abs =
                eigenpair_residual_norm_value_<_T,_Index>(*A, _T(lam), x);
            mv += 1;
            const scalar_real_type scale = vcp::tsparse::c1_revised_scale_(
                vcp::tsparse_scalar::abs_value(lam), anorm);
            return r_abs <= tol * scale;
        }
    } ks_lambda_gate = { &A, anorm_gate_ks, options.tol, sigma };
    std::size_t ks_lambda_gate_count = 0;
    // EIG-3 T-3: KS core in mu space (old arnoldi core + max_iter*(k+1)
    // expansion deleted; options.max_iter = total inner-apply budget)
    const ks_mu_pkg_<_T> pkg =
        ks_mu_core_drive_<_T>(n, k, options, sigma, apply_si,
                              ks_lambda_gate, true, &ks_lambda_gate_count);
    const std::size_t sdim = pkg.used_subspace_dim;

    eig_result<_T> result;
    result.requested_count = k;
    result.method = eig_solver_method::shift_invert_arnoldi;
    result.used_method = "shift_invert_arnoldi+user_preconditioner";
    result.used_orthogonalization = orthogonalization_to_string_<_T,_Index>(options.orthogonalization);
    result.iterations = pkg.iterations;
    result.matrix_vector_products = pkg.mv_count;
    result.lambda_gate_products = ks_lambda_gate_count;   // (e-2) 別建て計上
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
    // EIG-8 R4(裁定 s-2・必須併設): 返却対の複製方向監査(demote-only、
    // λ ゲート発火ラン限定)。
    if (result.converged && ks_lambda_gate_count > 0 &&
        vcp::tsparse::returned_pair_duplicate_direction_found_<_T>(result.eigenvectors)) {
        result.converged = false;
        if (result.failure_reason.empty())
            result.failure_reason = "duplicate direction among returned pairs "
                "(ghost-copy audit, EIG-8 R4): honest not_converged";
    }
    if (result.converged) {
        result.status = "converged";
        result.message = "converged";
    } else {
        result.status = "not_converged";
        if (result.failure_reason.empty()) result.failure_reason = "Arnoldi shift-invert did not converge";
        result.message = result.failure_reason;
    }
    // EIG-3 T-3: KS core exports C-2 evidence (D3-4) and the back-half
    // re-checked it in mu space; the EIG-1 arnoldi-core demotion is lifted.
    // Lambda-space C-1 remains the final acceptance gate.
    lambda_c1_acceptance_gate_<_T,_Index>(result, options.tol, self);   // EIG-4 T-4: 標準 si は改訂 scale
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    si_polish_rescue_standard_<_T,_Index>(result, self, options, apply_si);
    result.linear_solves = linear_solve_count;
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

    // EIG-8 T-1 (e): μ 面 opt-in の λ 形式 lock ゲート(一般化。
    // generalized_shift_invert_drive_ と同一式 — D8-2 共有ヘルパ、B-48)。
    const scalar_real_type anorm_gate_ks = matrix_inf_norm_value_<_T,_Index>(self);
    const scalar_real_type bnorm_gate_ks = matrix_inf_norm_value_<_T,_Index>(B);
    struct KsGenLambdaGateP {
        const spmats<_T,_Index>* A;
        const spmats<_T,_Index>* Bm;
        scalar_real_type anorm;
        scalar_real_type bnorm;
        scalar_real_type tol;
        scalar_real_type sigma;
        bool operator()(const std::vector<_T>& x, const _T& mu, std::size_t& mv) const {
            const scalar_real_type mu_re = vcp::tsparse_scalar::real_part(mu);
            if (!(vcp::tsparse_scalar::abs_value(mu_re) > scalar_real_type(0))) return false;
            const scalar_real_type lam = scalar_real_type(1) / mu_re + sigma;
            const scalar_real_type r_abs =
                generalized_eigenpair_residual_norm_value_<_T,_Index>(*A, *Bm, _T(lam), x);
            mv += 2;   // λ 判定の A·x と B·x(B-24: 1:1 計上)
            const scalar_real_type scale = vcp::tsparse::c1_generalized_scale_(
                vcp::tsparse_scalar::abs_value(lam), anorm, bnorm);
            return r_abs <= tol * scale;
        }
    } ks_lambda_gate = { &self, &B, anorm_gate_ks, bnorm_gate_ks, options.tol, sigma };
    std::size_t ks_lambda_gate_count = 0;
    // EIG-3 T-3: KS core in mu space (old arnoldi core + max_iter*(k+1)
    // expansion deleted; options.max_iter = total inner-apply budget)
    const ks_mu_pkg_<_T> pkg =
        ks_mu_core_drive_<_T>(n, k, options, sigma, apply_gsi,
                              ks_lambda_gate, true, &ks_lambda_gate_count);
    const std::size_t sdim = pkg.used_subspace_dim;

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
    result.lambda_gate_products = ks_lambda_gate_count;   // (e-2) 別建て計上
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
    // EIG-8 R4(裁定 s-2・必須併設): 返却対の複製方向監査(demote-only、
    // λ ゲート発火ラン限定)。
    if (result.converged && ks_lambda_gate_count > 0 &&
        vcp::tsparse::returned_pair_duplicate_direction_found_<_T>(result.eigenvectors)) {
        result.converged = false;
        if (result.failure_reason.empty())
            result.failure_reason = "duplicate direction among returned pairs "
                "(ghost-copy audit, EIG-8 R4): honest not_converged";
    }

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
    // EIG-3 T-3: KS core exports C-2 evidence (D3-4) and the back-half
    // re-checked it in mu space; the EIG-1 arnoldi-core demotion is lifted.
    // Lambda-space C-1 remains the final acceptance gate.
    // EIG-8 T-3 (D8-3): 一般化は D8-2 共有 scale 版(挙動遷移期待ゼロ —
    // 既存 rel 分岐が Frobenius 後退正規化を内包。契約整合の完成)。
    lambda_c1_acceptance_gate_<_T,_Index>(result, options.tol, self, B);
    // EIG-6 F-1: 条件付き磨き(ゲート不合格時のみ)+ 磨き solve の再同期(B-38)
    si_polish_rescue_generalized_<_T,_Index>(result, self, B, options, apply_gsi);
    result.linear_solves = linear_solve_count;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// ---------------------------------------------------------------------------
// make_empty_eigs_success_ - k==0 empty request: no solver launched
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> make_empty_eigs_success_(const spmats<_T,_Index>&,
                                                const eig_options<_T>& options)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    eig_result<_T> result;
    result.converged                = true;
    result.requested_count          = 0;
    result.returned_count           = 0;
    result.returned_real_count      = 0;
    result.returned_complex_count   = 0;
    result.converged_count          = 0;
    result.iterations               = 0;
    result.matrix_vector_products   = 0;
    result.linear_solves            = 0;
    result.residual_norm_absolute   = scalar_real_type(0);
    result.residual_norm_relative   = scalar_real_type(0);
    result.status                   = "empty_success";
    result.message                  = "requested zero eigenvalues";
    result.used_dense_fallback      = false;
    result.used_shift_invert        = false;
    result.used_generalized_operator = false;
    result.method                   = options.method;
    result.used_method              = eig_method_to_string_<_T,_Index>(options.method);
    return result;
}

// ---------------------------------------------------------------------------
// make_internal_error_eigs_result_  (SLU-GT1 D6)
//
// Result for the last-resort exception net of the policy_*_with_info
// boundaries.  Reached only if a std::exception (kv domain error, bad_alloc,
// ...) escapes the algorithm body: with the certified gates (D1/D3) in place
// this should not happen, so firing indicates a gate leak (investigate).
// Numerical singularity is NOT reported here -- the D3 gates report it as
// numerical_singularity / zero_pivot BEFORE any throw can occur.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> make_internal_error_eigs_result_(const std::size_t k,
                                                        const char* what)
{
    eig_result<_T> result;
    result.converged = false;
    result.requested_count = k;
    result.status = "internal_error";
    result.failure_reason = (what != 0) ? what : "unknown exception";
    result.message = result.failure_reason;
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
    // SLU-GT1 D6 net: vcp::error (validation/misuse) is rethrown unchanged;
    // any other std::exception maps to status "internal_error" (see
    // make_internal_error_eigs_result_).  Body indentation kept as-is.
    try {
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (k == 0) return make_empty_eigs_success_<_T,_Index>(A, options);
    if (spmatrix_is_complex<_T>::value)
        return complex_standard_eigs_with_info_<_T,_Index>(A, k, options);
    eig_options<_T> active = resolve_eigs_options_<_T,_Index>(A, options);
    validate_eigs_input_<_T,_Index>(A, k, "spmats::eigs", active.method);
    if (active.max_iter == 0 || active.tol <= scalar_real_type(0)) {
        vcp::throw_error<vcp::invalid_argument>("spmats::eigs: invalid iteration option");
    }
    // EIG-4 T-1/T-2: additive dispatch branches (B-26 -- the pre-existing
    // branches below are untouched).  auto_select is the new default (D4-2).
    if (active.method == eig_solver_method::auto_select)
        return auto_select_eigs_<_T,_Index>(A, k, active);
    if (active.method == eig_solver_method::thick_restart_lanczos)
        return trl_eigs_new_<_T,_Index>(A, k, active);
    if (active.method == eig_solver_method::krylov_schur)
        return ks_explicit_eigs_<_T,_Index>(A, k, active);
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
    // EIG-4 T-3 (additive, B-27): opt-in complex pairs on the explicit dense
    // path (real scalars, nonsymmetric spectra only -- the legacy branch
    // below is untouched when the opt-in is off).
    if (active.allow_complex_pairs && !spmatrix_is_complex<_T>::value
        && !vcp::tsparse_dense_linalg::is_dense_symmetric(dense, active.tol * scalar_real_type(10))) {
        return dense_pairs_eig_<_T,_Index>(dense, k, active);
    }
    eig_result<_T> result = dense_eig_<_T,_Index>(dense, active);
    result.used_dense_fallback = true;
    select_eigenpairs_aligned_<_T,_Index>(result, k, active.target, active.shift);
    set_result_counts_<_T,_Index>(result, k);
    return result;
    } catch (const vcp::error&) {
        throw;   // misuse contract preserved
    } catch (const std::exception& e) {
        return make_internal_error_eigs_result_<_T,_Index>(k, e.what());
    }
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
    // SLU-GT1 D6 net (same convention as the no-preconditioner overload).
    try {
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    if (k == 0) { (void)M; return make_empty_eigs_success_<_T,_Index>(A, options); }
    if (spmatrix_is_complex<_T>::value)
        return complex_preconditioned_eigs_unsupported_<_T,_Index>(k, options);
    eig_options<_T> active = resolve_eigs_options_<_T,_Index>(A, options);
    validate_eigs_input_<_T,_Index>(A, k, "spmats::eigs_with_info(with preconditioner)", active.method);
    if (active.max_iter == 0 || active.tol <= scalar_real_type(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "spmats::eigs_with_info(with preconditioner): invalid iteration option");
    }
    const scalar_real_type sigma = active.shift;
    // EIG-4 T-2: auto_select + use_shift on the preconditioner overload
    // follows the legacy default (lanczos) shift-invert path (additive OR
    // clauses only; the non-shift auto_select falls into the existing
    // unsupported_preconditioner branch below).
    const bool is_shift_invert =
        (active.method == eig_solver_method::shift_invert_lanczos)
        || (active.method == eig_solver_method::shift_invert_arnoldi)
        || (active.use_shift && active.method == eig_solver_method::lanczos)
        || (active.use_shift && active.method == eig_solver_method::arnoldi)
        || (active.use_shift && active.method == eig_solver_method::auto_select);

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
        || (active.use_shift && active.method == eig_solver_method::lanczos)
        || (active.use_shift && active.method == eig_solver_method::auto_select);

    eig_result<_T> result;
    if (use_lanczos)
        result = shift_invert_lanczos_eigs_with_prec_<_T,_Index>(A, k, active, sigma, M);
    else
        result = shift_invert_arnoldi_eigs_with_prec_<_T,_Index>(A, k, active, sigma, M);
    // EIG-4 T-2: routing diagnostic for the auto branch (additive)
    if (options.method == eig_solver_method::auto_select) {
        result.method = eig_solver_method::auto_select;
        result.used_method = "auto_select(" + result.used_method + ")";
    }
    set_result_counts_<_T,_Index>(result, k);
    return result;
    } catch (const vcp::error&) {
        throw;   // misuse contract preserved
    } catch (const std::exception& e) {
        return make_internal_error_eigs_result_<_T,_Index>(k, e.what());
    }
}

// ---------------------------------------------------------------------------
// policy_generalized_eigs_with_info  (no preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
eig_result<_T> policy_generalized_eigs_with_info(const spmats<_T,_Index>& A,
                                                   const spmats<_T,_Index>& B,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options_in)
{
    // SLU-GT1 D6 net (same convention as policy_eigs_with_info).
    try {
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    // EIG-4 T-2 (R-4, 設計 §4「一般化 不変」): 一般化 eigs の auto_select は
    // 入口で旧既定(lanczos)へ正規化し、挙動・診断文字列とも旧既定と同一に
    // する。以下の本体は options 名の付け替えのみで無変更(B-26)。
    eig_options<_T> options_norm = options_in;
    if (options_norm.method == eig_solver_method::auto_select)
        options_norm.method = eig_solver_method::lanczos;
    const eig_options<_T>& options = options_norm;
    if (k == 0) return make_empty_eigs_success_<_T,_Index>(A, options);
    if (spmatrix_is_complex<_T>::value)
        return complex_generalized_eigs_unsupported_<_T,_Index>(k, options);
    validate_generalized_eig_input_<_T,_Index>(A, B, "spmats::eigs(A,B)");

    const std::size_t n = static_cast<std::size_t>(A.rowsize());

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
    } catch (const vcp::error&) {
        throw;   // misuse contract preserved
    } catch (const std::exception& e) {
        return make_internal_error_eigs_result_<_T,_Index>(k, e.what());
    }
}

// ---------------------------------------------------------------------------
// policy_generalized_eigs_with_info  (with preconditioner)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index, class Preconditioner>
eig_result<_T> policy_generalized_eigs_with_info(const spmats<_T,_Index>& A,
                                                   const spmats<_T,_Index>& B,
                                                   const std::size_t k,
                                                   const eig_options<_T>& options_in,
                                                   const Preconditioner& M)
{
    // SLU-GT1 D6 net (same convention as policy_eigs_with_info).
    try {
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    // EIG-4 T-2 (R-4): 一般化の auto_select = 旧既定(lanczos)正規化(上の
    // 無前処理 overload と同一規則。本体無変更)。
    eig_options<_T> options_norm = options_in;
    if (options_norm.method == eig_solver_method::auto_select)
        options_norm.method = eig_solver_method::lanczos;
    const eig_options<_T>& options = options_norm;
    if (k == 0) { (void)B; (void)M; return make_empty_eigs_success_<_T,_Index>(A, options); }
    if (spmatrix_is_complex<_T>::value)
        return complex_generalized_eigs_unsupported_<_T,_Index>(k, options);
    validate_generalized_eig_input_<_T,_Index>(A, B, "spmats::eigs_with_info(A,B,with preconditioner)");

    const std::size_t n = static_cast<std::size_t>(A.rowsize());

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
    } catch (const vcp::error&) {
        throw;   // misuse contract preserved
    } catch (const std::exception& e) {
        return make_internal_error_eigs_result_<_T,_Index>(k, e.what());
    }
}

// ===========================================================================
// spmats<_T,_Index> member function definitions for policy_eigs_with_info
// and policy_generalized_eigs_with_info.
//
// NVI pattern (see sandbox/docs/design/spmats_finalize_policy.md §3):
// the non-Preconditioner overloads are non-virtual outers that finalize
// their spmats arguments and then delegate to a virtual _impl. Must never
// be overridden themselves -- override the _impl instead. The
// Preconditioner-templated overloads cannot be virtual (C++ forbids virtual
// template member functions), so they get a plain auto-finalize check
// instead of an _impl split; the Preconditioner M itself is never a bare
// spmats/spmatrix (it must expose apply(r,z), which spmats does not), and
// every concrete preconditioner (jacobi_preconditioner, ilu0_preconditioner,
// identity_preconditioner, function_preconditioner) already consumes its
// source matrix at construction time via get() or as_csr(), both of which
// are finalize-state-agnostic -- so M needs no finalize handling here.
// ===========================================================================

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_eigs_with_info(
    std::size_t k,
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    return policy_eigs_with_info_impl(k, opt);
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_eigs_with_info_impl(
    std::size_t k,
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    return vcp::policy_eigs_with_info<_T,_Index>(A, k, opt);
}

template <typename _T, typename _Index>
template <class Prec>
eig_result<_T> spmats<_T,_Index>::policy_eigs_with_info(
    std::size_t k,
    const eig_options<_T>& opt,
    const Prec& M) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    return vcp::policy_eigs_with_info<_T,_Index,Prec>(A, k, opt, M);
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eigs_with_info(
    const spmats<_T,_Index>& B,
    std::size_t k,
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    if (!B.is_finalized()) B.finalize();
    return policy_generalized_eigs_with_info_impl(B, k, opt);
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eigs_with_info_impl(
    const spmats<_T,_Index>& B,
    std::size_t k,
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    return vcp::policy_generalized_eigs_with_info<_T,_Index>(A, B, k, opt);
}

template <typename _T, typename _Index>
template <class Prec>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eigs_with_info(
    const spmats<_T,_Index>& B,
    std::size_t k,
    const eig_options<_T>& opt,
    const Prec& M) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    if (!B.is_finalized()) B.finalize();
    return vcp::policy_generalized_eigs_with_info<_T,_Index,Prec>(A, B, k, opt, M);
}

// ===========================================================================
// Strict policy methods: policy_eig, policy_eigs, policy_generalized_eig,
// policy_generalized_eigs (and preconditioner overloads).
// All convergence/count checking lives here, not in spmatrix.hpp.
//
// policy_eig / policy_eigs / policy_generalized_eig / policy_generalized_eigs
// (no Preconditioner) are NVI outers too: they finalize their spmats
// arguments and delegate to a virtual _impl, even though their _impl bodies
// only call back into policy_eigs_with_info / policy_generalized_eigs_with_info
// (already finalize-safe on their own). This keeps all six eigs entry points
// independently safe and overridable, per spmats_finalize_policy.md §3.
// ===========================================================================

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_eig(
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    return policy_eig_impl(opt);
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_eig_impl(
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
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
    eig_result<_T> result = policy_eigs_with_info(static_cast<std::size_t>(A.rowsize()), opt);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_eig: eigensolver did not converge");
    return result;
}

template <typename _T, typename _Index>
std::vector<_T> spmats<_T,_Index>::policy_eigs(
    std::size_t k,
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    return policy_eigs_impl(k, opt);
}

template <typename _T, typename _Index>
std::vector<_T> spmats<_T,_Index>::policy_eigs_impl(
    std::size_t k,
    const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    eig_result<_T> result = policy_eigs_with_info(k, opt);
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
    std::size_t k,
    const eig_options<_T>& opt, const Prec& M) const
{
    // policy_eigs_with_info(k,opt,M) below already auto-finalizes *this.
    eig_result<_T> result = policy_eigs_with_info(k, opt, M);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_eigs(with preconditioner): eigensolver did not converge");
    if (result.returned_real_count < k)
        vcp::throw_error<vcp::state_error>("spmats::policy_eigs(with preconditioner): insufficient eigenvalues returned");
    return result.eigenvalues;
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eig(
    const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    if (!B.is_finalized()) B.finalize();
    return policy_generalized_eig_impl(B, k, opt);
}

template <typename _T, typename _Index>
eig_result<_T> spmats<_T,_Index>::policy_generalized_eig_impl(
    const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt) const
{
    eig_result<_T> result = policy_generalized_eigs_with_info(B, k, opt);
    if (!result.converged)
        vcp::throw_error<vcp::state_error>("spmats::policy_generalized_eig: eigensolver did not converge");
    return result;
}

template <typename _T, typename _Index>
std::vector<_T> spmats<_T,_Index>::policy_generalized_eigs(
    const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    if (!A.is_finalized()) A.finalize();
    if (!B.is_finalized()) B.finalize();
    return policy_generalized_eigs_impl(B, k, opt);
}

template <typename _T, typename _Index>
std::vector<_T> spmats<_T,_Index>::policy_generalized_eigs_impl(
    const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    eig_result<_T> result = policy_generalized_eigs_with_info(B, k, opt);
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
    const spmats<_T,_Index>& B,
    std::size_t k, const eig_options<_T>& opt, const Prec& M) const
{
    const spmats<_T,_Index>& A = *this;   // WFIX-2: subject is *this
    // policy_generalized_eigs_with_info(B,k,opt,M) below already
    // auto-finalizes *this and B.
    eig_result<_T> result = policy_generalized_eigs_with_info(B, k, opt, M);
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
