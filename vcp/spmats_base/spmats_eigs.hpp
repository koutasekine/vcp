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
#include <vcp/tsparse/tsparse_krylov_schur.hpp>   // EIG-3: arnoldi/si successor
#include <vcp/tsparse/tsparse_thick_restart_lanczos.hpp>   // EIG-4: TRL promotion (T-1/T-2)
#include <vcp/tsparse/tsparse_factorization.hpp>
#include <vcp/tsparse/tsparse_preconditioner.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_generalized_shift_invert.hpp>
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
template <class T> struct trl_is_supported_scalar_;   // EIG-T2 D-1
template <class Apply, class T>
vcp::eig_result<T> thick_restart_lanczos_eigs(
    const Apply& apply, std::size_t n, std::size_t k,
    const vcp::eig_options<T>& options);
// EIG-T1 T-2: the kernel without its end-of-run residual applications
template <class Apply, class T>
vcp::eig_result<T> thick_restart_lanczos_eigs_no_residual_(
    const Apply& apply, std::size_t n, std::size_t k,
    const vcp::eig_options<T>& options);
// EIG-T1 P5a: the kernel in the B inner product
template <class T> struct thick_restart_lanczos_result;
template <class Apply, class BMul, class T>
thick_restart_lanczos_result<T> thick_restart_lanczos_eigs_b_inner_(
    const Apply& apply, const BMul& bmul, std::size_t n, std::size_t k,
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
    const spmats<_T,_Index>* Cp; spmats<_T,_Index> Cc_storage;
    if (A.is_finalized() && A.format() == vcp::sparse_csr) { Cp = &A; }
    else { Cc_storage = A.as_csr(); Cp = &Cc_storage; }
    const spmats<_T,_Index>& C = *Cp;
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
// EIG-T1 P5b: the verdict of every route and every scalar type is the
// solver's own (T-9).  The lambda-space end-of-run gates
// (lambda_c1_acceptance_gate_) and the one-step polish behind them
// (si_polish_*) were removed; the residual fields stay as diagnostics.
// ---------------------------------------------------------------------------

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
// normalize_eig_result_counts_: the count fields describe the returned arrays.
//   returned_count         = eigenvalues.size()
//   returned_complex_count = entries of eigenvalues_imag that are not zero
//   returned_real_count    = returned_count - returned_complex_count
// ---------------------------------------------------------------------------
template <typename _T>
static eig_result<_T> normalize_eig_result_counts_(eig_result<_T> result)
{
    std::size_t complex_entries = 0;
    for (std::size_t i = 0; i < result.eigenvalues_imag.size(); i++) {
        if (!(result.eigenvalues_imag[i] == _T(0))) complex_entries++;
    }
    result.returned_count = result.eigenvalues.size();
    result.returned_complex_count = complex_entries;
    result.returned_real_count = result.returned_count - complex_entries;
    return result;
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
    // EIG-T1 (G-1 ruling A): the same input check for the two other
    // symmetric-only methods when they are named explicitly (after the
    // use_shift resolution above, so lanczos + use_shift is covered too).
    // With structure = symmetric the caller is trusted, as for lanczos.
    // A check of the precondition on the input, not of a residual.
    if (active.structure == matrix_structure_hint::auto_detect
     && (active.method == eig_solver_method::thick_restart_lanczos
      || active.method == eig_solver_method::shift_invert_lanczos)
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
    // EIG-T1 P6b: a request that is a shift-invert request on the general
    // routes (use_shift, or shift_invert_lanczos / shift_invert_arnoldi named)
    // returns the k values nearest sigma = options.shift whatever the target
    // is, as every other route does (generalized / standard shift-invert,
    // spumar).  Other requests select by options.target as before.
    const bool nearest_shift =
        options.use_shift
        || options.method == eig_solver_method::shift_invert_lanczos
        || options.method == eig_solver_method::shift_invert_arnoldi;
    std::vector<std::size_t> order;
    if (nearest_shift) {
        std::vector<scalar_real_type> dist(n);
        for (std::size_t i = 0; i < n; i++)
            dist[i] = vcp::tsparse_scalar::abs_value(
                scalar_real_type(vcp::tsparse_scalar::real_part(values[i]) - options.shift));
        order.resize(n);
        for (std::size_t i = 0; i < n; i++) order[i] = i;
        std::stable_sort(order.begin(), order.end(),
            [&dist](const std::size_t x, const std::size_t y) { return dist[x] < dist[y]; });
        if (order.size() > k) order.resize(k);
    } else {
        order = vcp::tsparse_eigen_selection::select_real_eigenpairs(values, k, options.target, options.shift);
    }
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
    return vcp::tsparse_scalar::hermitian_norm_value(r);
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
    const scalar_real_type vn = vcp::tsparse_scalar::hermitian_norm_value(v);
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
//
// EIG-T2 D-3: the complex Hermitian standard route (no shift) runs on the
// thick-restart kernel, in the shape of the real symmetric route T-2c
// (t2_lanczos_plain_).  Declared here for complex_standard_eigs_with_info_;
// defined after the T-7c helpers it uses.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> hermitian_lanczos_eigs_impl_(const spmats<_T,_Index>& self,
                                                     const std::size_t k,
                                                     const eig_options<_T>& options);

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
    // hermitian 経路)へ正規化する(挙動は旧既定と同一)。複素 T には提供
    // しない method = thick_restart_lanczos / krylov_schur の明示指定は
    // 正直拒否(EIG-T2 以降、thick-restart の核そのものは複素エルミートを
    // 受け付け、method = lanczos の経路が使う。krylov_schur の核は実 T
    // 専用のまま)。以下の本体は options 名の付け替えのみで無変更(B-26)。
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
            "method=thick_restart_lanczos / krylov_schur is not offered for"
            " complex scalars; for complex Hermitian matrices use method=lanczos";
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
// EIG-T1 T-2: the four Lanczos-type routes on the T-1 thick-restart kernel
// (they replaced the plain-restart cores lanczos_eigs_standard /
// lanczos_eigs_standard_si_ / b_inner_lanczos_eigs).
//
//   T-2a / T-2b  standard shift-invert (with or without a user
//                preconditioner): the operator (A - sigma I)^{-1} goes to the
//                kernel (largest |mu|), lambda = sigma + 1 / mu
//   T-2c         standard, method = lanczos, no shift: the operator A
//   T-2d         generalized, no shift: B = P^T L L^T P (sparse Cholesky),
//                the operator C = L^{-1} P A P^T L^{-T}, x = P^T L^{-T} y
//
// max_iter is the total number of operator applications.  The residual
// fields are true residuals of the original problem, computed after the run
// as diagnostics; they do not decide `converged`.
// Tag: t2_supported_ (the real point scalars the kernel accepts,
// trl_is_real_floating; since EIG-T2 the kernel also accepts their
// std::complex, which has its own route hermitian_trl_route_ and does not
// use these four); the false overloads are compile-time fences for the
// other scalar types (interval), which are refused before they are reached
// (t2_unsupported_scalar_).
// ---------------------------------------------------------------------------
template <typename _T>
struct t2_supported_ {
    static const bool value =
        vcp::tsparse_experimental::trl_is_real_floating<_T>::value;
};

// EIG-T1 (G-1 ruling E): scalar types the kernel does not accept (interval
// scalars) have no Lanczos-type route.  lanczos / shift_invert_lanczos (with
// or without preconditioner) and the generalized problem without shift are
// refused with status "unsupported" (the shape of the thick_restart_lanczos
// refusal in trl_drive_).  krylov_schur / arnoldi / dense stay available.
template <typename _T, typename _Index>
static eig_result<_T> t2_unsupported_scalar_(const std::size_t k,
                                             const eig_solver_method method,
                                             const std::string& used_method)
{
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    result.failure_reason = used_method +
        ": scalar type unsupported (real point scalars only);"
        " use krylov_schur / arnoldi (shift_invert_arnoldi with a shift)";
    result.message = result.failure_reason;
    result.method = method;
    result.used_method = used_method;
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// status / message of the lanczos-type routes ("converged" /
// "not_converged") from a kernel result
template <typename _T, typename _Index>
static void t2_finish_lanczos_type_(eig_result<_T>& result,
                                    const eig_solver_method method,
                                    const std::size_t k)
{
    // These routes never report `converged` with fewer pairs than requested
    // (the kernel clamps k to n and judges the clamped count).
    if (result.converged && result.eigenvalues.size() < k) {
        result.converged = false;
        if (result.failure_reason.empty())
            result.failure_reason = "fewer eigenpairs than requested (k exceeds the dimension)";
    }
    const std::string kernel_failure = result.failure_reason;
    const std::size_t sdim = result.used_subspace_dim;
    const std::size_t mv = result.matrix_vector_products;
    result.requested_count = k;
    set_eig_diagnostics_<_T,_Index>(result, method, sdim, mv,
        result.breakdown_reason, kernel_failure);
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    set_result_counts_<_T,_Index>(result, k);
}

// --- T-2a / T-2b -----------------------------------------------------------
template <typename _T, typename _Index, class Apply>
static eig_result<_T> t2_si_lanczos_core_(const spmats<_T,_Index>&,
                                          const std::size_t k,
                                          const eig_options<_T>&,
                                          const typename vcp::tsparse_scalar::real_type<_T>::type&,
                                          Apply&,
                                          std::false_type /* t2_supported */)
{
    eig_result<_T> result;   // fence (never executed)
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    return result;
}

template <typename _T, typename _Index, class Apply>
static eig_result<_T> t2_si_lanczos_core_(const spmats<_T,_Index>& self,
                                          const std::size_t k,
                                          const eig_options<_T>& options,
                                          const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                          Apply& apply_si,
                                          std::true_type /* t2_supported */)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    eig_options<_T> mu_opts = options;
    mu_opts.target = eig_target::largest_magnitude;   // nearest sigma in lambda
    mu_opts.shift = scalar_real_type(0);
    mu_opts.use_shift = false;
    eig_result<_T> result =
        vcp::tsparse_experimental::thick_restart_lanczos_eigs_no_residual_<Apply, _T>(
            apply_si, n, k, mu_opts);
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(result.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            result.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }
    if (!result.eigenvectors.empty()) {
        spmats<_T,_Index> A = self.as_csr();
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    t2_finish_lanczos_type_<_T,_Index>(result, eig_solver_method::shift_invert_lanczos, k);
    result.used_shift_invert = true;
    return result;
}

// --- T-2c --------------------------------------------------------------------
template <typename _T, typename _Index>
static eig_result<_T> t2_lanczos_plain_(const spmats<_T,_Index>&,
                                        const std::size_t k,
                                        const eig_options<_T>&,
                                        std::false_type /* t2_supported */)
{
    eig_result<_T> result;   // fence (never executed)
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> t2_lanczos_plain_(const spmats<_T,_Index>& self,
                                        const std::size_t k,
                                        const eig_options<_T>& options,
                                        std::true_type /* t2_supported */)
{
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    };
    apply_fn apply_op = { &A };
    eig_result<_T> result =
        vcp::tsparse_experimental::thick_restart_lanczos_eigs_no_residual_<apply_fn, _T>(
            apply_op, n, k, options);
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    t2_finish_lanczos_type_<_T,_Index>(result, eig_solver_method::lanczos, k);
    result.used_shift_invert = false;
    return result;
}

// --- T-2d --------------------------------------------------------------------
template <typename _T, typename _Index>
static bool t2_generalized_chol_(const spmats<_T,_Index>&,
                                 const spmats<_T,_Index>&,
                                 const std::size_t,
                                 const eig_options<_T>&,
                                 eig_result<_T>&,
                                 std::false_type /* t2_supported */)
{
    return false;   // fence: the caller refuses the scalar type
}

// Returns true when the route produced `result` (also for an honest
// factorization failure); false = not applicable.
template <typename _T, typename _Index>
static bool t2_generalized_chol_(const spmats<_T,_Index>& self,
                                 const spmats<_T,_Index>& B,
                                 const std::size_t k,
                                 const eig_options<_T>& options,
                                 eig_result<_T>& result,
                                 std::true_type /* t2_supported */)
{
    // The sparse Cholesky needs a signed index type; the factor is held with
    // std::ptrdiff_t indices whatever _Index is (index conversion only).
    typedef std::ptrdiff_t CI;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());

    vcp::sparse_chol_result<_T, CI> f;
    bool chol_ok = false;
    std::string chol_note;
    try {
        const spmats<_T,_Index> Bc = B.as_csc();
        const std::vector<_Index>& bo = Bc.outer_index();
        const std::vector<_Index>& bi = Bc.inner_index();
        std::vector<CI> cp(bo.size()), ri(bi.size());
        for (std::size_t i = 0; i < bo.size(); i++) cp[i] = static_cast<CI>(bo[i]);
        for (std::size_t i = 0; i < bi.size(); i++) ri[i] = static_cast<CI>(bi[i]);
        vcp::sparse_chol_options<_T> copt;
        copt.check_symmetry = false;   // symmetry of B was checked by the caller (tolerance 1e-10)
        f = vcp::sparse_chol_factorize_with_info<_T, CI>(
            static_cast<CI>(n), cp, ri, Bc.values(), copt);
        chol_ok = (f.status == vcp::sparse_chol_status::success);
        if (!chol_ok) {
            if (f.status == vcp::sparse_chol_status::not_positive_definite)
                chol_note = "B is not positive definite";
            else if (f.status == vcp::sparse_chol_status::inconclusive_pivot_test)
                chol_note = "positive definiteness of B could not be decided";
            else
                chol_note = "sparse Cholesky of B did not complete";
        }
    } catch (const std::exception& e) {
        chol_ok = false;
        chol_note = std::string("sparse Cholesky of B threw: ") + e.what();
    }
    if (!chol_ok) {
        result.converged = false;
        result.status = "factorization_failed";
        result.failure_reason = "generalized eigs without shift: " + chol_note;
        result.message = result.failure_reason;
        set_result_counts_<_T,_Index>(result, k);
        return true;
    }

    spmats<_T,_Index> A = self.as_csr();
    struct chol_apply {
        const spmats<_T,_Index>* A;
        const vcp::sparse_chol_result<_T, CI>* f;
        std::size_t n;
        // z <- L^{-T} z
        void lt_solve(std::vector<_T>& z) const {
            for (std::size_t jj = n; jj > 0; jj--) {
                const std::size_t j = jj - 1;
                const std::size_t p0 = static_cast<std::size_t>(f->L_col_ptr[j]);
                const std::size_t p1 = static_cast<std::size_t>(f->L_col_ptr[j + 1]);
                _T sum = z[j];
                for (std::size_t q = p0 + 1; q < p1; q++)
                    sum -= f->L_val[q] * z[static_cast<std::size_t>(f->L_row_ind[q])];
                z[j] = sum / f->L_val[p0];
            }
        }
        // z <- L^{-1} z
        void l_solve(std::vector<_T>& z) const {
            for (std::size_t j = 0; j < n; j++) {
                const std::size_t p0 = static_cast<std::size_t>(f->L_col_ptr[j]);
                const std::size_t p1 = static_cast<std::size_t>(f->L_col_ptr[j + 1]);
                z[j] = z[j] / f->L_val[p0];
                for (std::size_t q = p0 + 1; q < p1; q++)
                    z[static_cast<std::size_t>(f->L_row_ind[q])] -= f->L_val[q] * z[j];
            }
        }
        // x = P^T L^{-T} y
        std::vector<_T> back(const std::vector<_T>& y) const {
            std::vector<_T> z(y), x(n, _T(0));
            lt_solve(z);
            for (std::size_t i = 0; i < n; i++) x[static_cast<std::size_t>(f->perm[i])] = z[i];
            return x;
        }
        // y = L^{-1} P A P^T L^{-T} x
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
            const std::vector<_T> u = back(x);
            const std::vector<_T> t = A->mul_vec(u);
            y.resize(n);
            for (std::size_t i = 0; i < n; i++) y[i] = t[static_cast<std::size_t>(f->perm[i])];
            l_solve(y);
        }
    };
    chol_apply apply_c = { &A, &f, n };

    const eig_result<_T> kr =
        vcp::tsparse_experimental::thick_restart_lanczos_eigs_no_residual_<chol_apply, _T>(
            apply_c, n, k, options);

    result.iterations = kr.iterations;
    result.matrix_vector_products = kr.matrix_vector_products;
    result.used_subspace_dim = kr.used_subspace_dim;
    result.breakdown_reason = kr.breakdown_reason;
    result.failure_reason = kr.failure_reason;
    result.converged_count = kr.converged_count;
    result.eigenvalues = kr.eigenvalues;
    result.eigenvectors.clear();
    for (std::size_t i = 0; i < kr.eigenvectors.size(); i++)
        result.eigenvectors.push_back(apply_c.back(kr.eigenvectors[i]));
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = generalized_eigenpair_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = generalized_eigenpair_relative_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
    }
    result.converged = kr.converged;
    if (options.compute_residual_history) {
        result.residual_history_absolute = kr.residual_history_absolute;
        result.residual_history_relative = kr.residual_history_relative;
    }
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    if (result.converged) {
        result.status = "converged";
        result.message = "converged";
        result.failure_reason.clear();
    } else if (kr.status == "max_iter_exhausted") {
        result.status = "max_iter_exhausted";
        if (result.failure_reason.empty())
            result.failure_reason = "MV budget (max_iter) exhausted before convergence";
        result.message = result.failure_reason;
    } else {
        result.status = "not_converged";
        if (result.failure_reason.empty())
            result.failure_reason = "generalized thick-restart Lanczos did not converge";
        result.message = result.failure_reason;
    }
    set_result_counts_<_T,_Index>(result, k);
    return true;
}

// ---------------------------------------------------------------------------
// EIG-T1 P5a: generalized shift-invert of a symmetric pair on the T-1 kernel
//
// A x = lambda B x with A, B symmetric and a Lanczos-type request
// (shift_invert_lanczos, lanczos + use_shift, auto_select + use_shift; all
// three arrive here as promoted_from_lanczos).  The iteration has the shape
// of ARPACK mode 3: OP = (A - sigma B)^{-1} B is symmetric in the B inner
// product, the kernel runs on OP with that inner product (largest |mu|),
// lambda = sigma + 1 / mu.  B is not factorized and may be semidefinite.
// Convergence / kept pairs / termination are the kernel's (T-1, T-7, T-8,
// T-9); no end-of-run residual test, no polish, no purification.
// When a B inner product turns out negative or not finite the kernel stops
// and the request is served by the Krylov-Schur route as before.
// Requests outside the conditions (not symmetric, shift_invert_arnoldi named,
// scalar types the kernel does not accept) keep the Krylov-Schur route.
// used_method / status strings and the meaning of max_iter (total number of
// applications of OP, the one for the start vector included) are unchanged.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
static bool p5a_symmetric_pair_(const spmats<_T,_Index>& A,
                                const spmats<_T,_Index>& B,
                                const eig_options<_T>& options)
{
    // the symmetry decision of auto_select_eigs_ (hint first, else certain)
    if (options.structure == matrix_structure_hint::symmetric) return true;
    if (options.structure == matrix_structure_hint::general)   return false;
    return is_certainly_symmetric_<_T,_Index>(A) && is_certainly_symmetric_<_T,_Index>(B);
}

template <typename _T, typename _Index>
static bool p5a_route_(const spmats<_T,_Index>& A,
                       const spmats<_T,_Index>& B,
                       const eig_options<_T>& options,
                       const bool promoted_from_lanczos)
{
    if (!t2_supported_<_T>::value) return false;
    if (!promoted_from_lanczos) return false;
    return p5a_symmetric_pair_<_T,_Index>(A, B, options);
}

// Symmetric factorization of K = A - sigma B for the default (direct) entry:
//   (a) sparse Cholesky; (b) sparse LDL^T (supernodal, Bunch-Kaufman; after
//   pivot_out_of_panel once more with baseline_dynamic); kind == 0 when
//   neither succeeded (the caller then uses the sparse LU).
// solve(): K^{-1} b in place.  LDL^T: P^T K P = L D L^T (perm new->old) --
// permute, unit lower forward substitution, 1x1 / 2x2 block diagonal,
// back substitution, permute back.  No iterative refinement.
template <typename _T, typename _Index>
struct p5a_sym_factor_ {
    typedef std::ptrdiff_t CI;
    int kind;                    // 0 none, 1 Cholesky, 2 LDL^T
    std::string diagnostics;
    std::size_t n;
    vcp::sparse_chol_result<_T, CI> chol;
    std::vector<std::size_t> l_ptr, l_row, perm;   // strict lower part of L by columns
    std::vector<_T> l_val, d_diag, d_sub;          // d_sub[k] = D(k+1, k)
    std::vector<char> d_two;                       // 1: a 2x2 block starts at k

    p5a_sym_factor_() : kind(0), n(0) {}

    void factorize(const spmats<_T,_Index>& K)
    {
        kind = 0;
        n = static_cast<std::size_t>(K.rowsize());
        std::string note_a, note_b;
        // (a) Cholesky
        try {
            const spmats<_T,_Index> Kc = K.as_csc();
            const std::vector<_Index>& ko = Kc.outer_index();
            const std::vector<_Index>& ki = Kc.inner_index();
            std::vector<CI> cp(ko.size()), ri(ki.size());
            for (std::size_t i = 0; i < ko.size(); i++) cp[i] = static_cast<CI>(ko[i]);
            for (std::size_t i = 0; i < ki.size(); i++) ri[i] = static_cast<CI>(ki[i]);
            vcp::sparse_chol_options<_T> copt;
            copt.check_symmetry = false;   // symmetry was decided by the caller
            chol = vcp::sparse_chol_factorize_with_info<_T, CI>(
                static_cast<CI>(n), cp, ri, Kc.values(), copt);
            if (chol.status == vcp::sparse_chol_status::success) {
                kind = 1;
                diagnostics = "shift_invert_factorization=cholesky";
                return;
            }
            note_a = vcp::sparse_chol_status_to_string(chol.status);
        } catch (const std::exception& e) {
            note_a = std::string("exception: ") + e.what();
        }
        // (b) LDL^T
        for (int attempt = 0; attempt < 2; attempt++) {
            const bool sn = (attempt == 0);
            vcp::sparse_ldl_status st = vcp::sparse_ldl_status::internal_error;
            try {
                vcp::ldl_options<_T> lopt;
                lopt.method = sn ? vcp::sparse_ldl_method::supernodal
                                 : vcp::sparse_ldl_method::baseline_dynamic;
                lopt.pivoting = vcp::sparse_ldl_pivoting::bk;
                spmats<_T,_Index> L, D;
                std::vector<_Index> pm;
                const vcp::ldl_result<_T,_Index> lr = K.policy_ldl_with_info(L, D, pm, lopt);
                st = lr.status;
                if (st == vcp::sparse_ldl_status::success && pm.size() == n) {
                    load_ldl_(L, D, pm);
                    kind = 2;
                    diagnostics = std::string("shift_invert_factorization=ldlt(method=")
                        + (sn ? "supernodal" : "baseline_dynamic")
                        + ",pivoting=bk); cholesky: " + note_a
                        + (note_b.empty() ? std::string() : ("; " + note_b));
                    return;
                }
                note_b += std::string(note_b.empty() ? "" : ", ")
                    + "ldlt(" + (sn ? "supernodal" : "baseline_dynamic") + "): "
                    + vcp::sparse_ldl_status_to_string(st);
            } catch (const std::exception& e) {
                note_b += std::string(note_b.empty() ? "" : ", ")
                    + "ldlt(" + (sn ? "supernodal" : "baseline_dynamic") + "): exception: " + e.what();
                break;
            }
            if (!(sn && st == vcp::sparse_ldl_status::pivot_out_of_panel)) break;
        }
        diagnostics = "cholesky: " + note_a + "; " + note_b;
    }

    void load_ldl_(const spmats<_T,_Index>& L, const spmats<_T,_Index>& D,
                   const std::vector<_Index>& pm)
    {
        const spmats<_T,_Index> Lc = L.as_csc();
        const std::vector<_Index>& lo = Lc.outer_index();
        const std::vector<_Index>& li = Lc.inner_index();
        const std::vector<_T>& lv = Lc.values();
        l_ptr.assign(n + 1, 0);
        l_row.clear();
        l_val.clear();
        for (std::size_t j = 0; j < n; j++) {
            for (std::size_t q = static_cast<std::size_t>(lo[j]); q < static_cast<std::size_t>(lo[j + 1]); q++) {
                const std::size_t i = static_cast<std::size_t>(li[q]);
                if (i <= j) continue;   // the explicit unit diagonal
                l_row.push_back(i);
                l_val.push_back(lv[q]);
            }
            l_ptr[j + 1] = l_row.size();
        }
        d_diag.assign(n, _T(0));
        d_sub.assign(n, _T(0));
        d_two.assign(n, char(0));
        {
            // an entry (o, inner) is (row, col) or (col, row); D is symmetric
            const std::vector<_Index>& dout = D.outer_index();
            const std::vector<_Index>& din = D.inner_index();
            const std::vector<_T>& dv = D.values();
            for (std::size_t o = 0; o < n; o++) {
                for (std::size_t q = static_cast<std::size_t>(dout[o]); q < static_cast<std::size_t>(dout[o + 1]); q++) {
                    const std::size_t i = static_cast<std::size_t>(din[q]);
                    if (i == o) d_diag[o] = dv[q];
                    else if (i + 1 == o) { d_sub[i] = dv[q]; d_two[i] = char(1); }
                    else if (o + 1 == i) { d_sub[o] = dv[q]; d_two[o] = char(1); }
                }
            }
        }
        perm.resize(n);
        for (std::size_t i = 0; i < n; i++) perm[i] = static_cast<std::size_t>(pm[i]);
    }

    void solve(std::vector<_T>& b) const
    {
        std::vector<_T> z(n);
        if (kind == 1) {
            for (std::size_t i = 0; i < n; i++) z[i] = b[static_cast<std::size_t>(chol.perm[i])];
            for (std::size_t j = 0; j < n; j++) {
                const std::size_t p0 = static_cast<std::size_t>(chol.L_col_ptr[j]);
                const std::size_t p1 = static_cast<std::size_t>(chol.L_col_ptr[j + 1]);
                z[j] = z[j] / chol.L_val[p0];
                for (std::size_t q = p0 + 1; q < p1; q++)
                    z[static_cast<std::size_t>(chol.L_row_ind[q])] -= chol.L_val[q] * z[j];
            }
            for (std::size_t jj = n; jj > 0; jj--) {
                const std::size_t j = jj - 1;
                const std::size_t p0 = static_cast<std::size_t>(chol.L_col_ptr[j]);
                const std::size_t p1 = static_cast<std::size_t>(chol.L_col_ptr[j + 1]);
                _T sum = z[j];
                for (std::size_t q = p0 + 1; q < p1; q++)
                    sum -= chol.L_val[q] * z[static_cast<std::size_t>(chol.L_row_ind[q])];
                z[j] = sum / chol.L_val[p0];
            }
            for (std::size_t i = 0; i < n; i++) b[static_cast<std::size_t>(chol.perm[i])] = z[i];
            return;
        }
        for (std::size_t i = 0; i < n; i++) z[i] = b[perm[i]];
        for (std::size_t j = 0; j < n; j++) {
            const _T zj = z[j];
            for (std::size_t q = l_ptr[j]; q < l_ptr[j + 1]; q++) z[l_row[q]] -= l_val[q] * zj;
        }
        for (std::size_t kk = 0; kk < n; ) {
            if (d_two[kk] != char(0) && kk + 1 < n) {
                // [a e; e c]^{-1} with the scaling of LAPACK's xSYTRS
                const _T e = d_sub[kk];
                const _T am = d_diag[kk] / e;
                const _T cm = d_diag[kk + 1] / e;
                const _T den = am * cm - _T(1);
                const _T b0 = z[kk] / e;
                const _T b1 = z[kk + 1] / e;
                z[kk]     = (cm * b0 - b1) / den;
                z[kk + 1] = (am * b1 - b0) / den;
                kk += 2;
            } else {
                z[kk] = z[kk] / d_diag[kk];
                kk += 1;
            }
        }
        for (std::size_t jj = n; jj > 0; jj--) {
            const std::size_t j = jj - 1;
            _T sum = z[j];
            for (std::size_t q = l_ptr[j]; q < l_ptr[j + 1]; q++) sum -= l_val[q] * z[l_row[q]];
            z[j] = sum;
        }
        for (std::size_t i = 0; i < n; i++) b[perm[i]] = z[i];
    }
};

// The kernel run in the B inner product on apply_si: y = (A - sigma B)^{-1} B x.
// Returns false when the route did not produce a result (B inner product
// negative / not finite, or the scalar fence); the caller then continues
// with the Krylov-Schur route.  On true, `result` holds values, Euclidean
// unit vectors, true residuals of the original pair (diagnostics) and the
// kernel's verdict in `converged`; method / status are the caller's.
template <typename _T, typename _Index, class Apply>
static bool p5a_gsi_lanczos_core_(const spmats<_T,_Index>&,
                                  const spmats<_T,_Index>&,
                                  const std::size_t,
                                  const eig_options<_T>&,
                                  const typename vcp::tsparse_scalar::real_type<_T>::type&,
                                  Apply&,
                                  eig_result<_T>&,
                                  std::false_type /* t2_supported */)
{
    return false;   // fence (p5a_route_ is false for these scalar types)
}

template <typename _T, typename _Index, class Apply>
static bool p5a_gsi_lanczos_core_(const spmats<_T,_Index>& self,
                                  const spmats<_T,_Index>& B,
                                  const std::size_t k,
                                  const eig_options<_T>& options,
                                  const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                  Apply& apply_si,
                                  eig_result<_T>& result,
                                  std::true_type /* t2_supported */)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    const spmats<_T,_Index> Bc = B.as_csr();
    struct BMul {
        const spmats<_T,_Index>* Bm;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = Bm->mul_vec(x); }
    } bmul = { &Bc };
    eig_options<_T> mu_opts = options;
    mu_opts.target = eig_target::largest_magnitude;   // nearest sigma in lambda
    mu_opts.shift = scalar_real_type(0);
    mu_opts.use_shift = false;
    const vcp::tsparse_experimental::thick_restart_lanczos_result<_T> tr =
        vcp::tsparse_experimental::thick_restart_lanczos_eigs_b_inner_<Apply, BMul, _T>(
            apply_si, bmul, n, k, mu_opts);
    if (tr.b_inner_invalid) return false;
    result = tr.eigs;
    for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
        const scalar_real_type mu = vcp::tsparse_scalar::real_part(result.eigenvalues[i]);
        if (vcp::tsparse_scalar::abs_value(mu) > scalar_real_type(0))
            result.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
    }
    for (std::size_t i = 0; i < result.eigenvectors.size(); i++) {
        std::vector<_T>& v = result.eigenvectors[i];
        const scalar_real_type nv = vcp::tsparse_scalar::real_norm_value(v);
        if (nv > scalar_real_type(0))
            for (std::size_t j = 0; j < v.size(); j++) v[j] = v[j] / _T(nv);
    }
    result.residuals_absolute.clear();
    result.residuals_relative.clear();
    result.residual_norm_absolute = scalar_real_type(0);
    result.residual_norm_relative = scalar_real_type(0);
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = generalized_eigenpair_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = generalized_eigenpair_relative_residuals_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        result.residual_norm_absolute = max_generalized_eigenpair_residual_value_<_T,_Index>(self, B, result.eigenvalues, result.eigenvectors);
        for (std::size_t i = 0; i < result.residuals_relative.size(); i++)
            if (result.residuals_relative[i] > result.residual_norm_relative)
                result.residual_norm_relative = result.residuals_relative[i];
    }
    if (result.converged && result.eigenvalues.size() < k) {
        result.converged = false;
        if (result.failure_reason.empty())
            result.failure_reason = "fewer eigenpairs than requested (k exceeds the dimension)";
    }
    result.requested_count = k;
    result.used_shift_invert = true;
    result.used_generalized_operator = true;
    result.used_dense_fallback = false;
    populate_real_complex_eigenvalues_<_T,_Index>(result);
    return true;
}

// method / used_method / status of a P5a result, as the Krylov-Schur route
// writes them for the same request
template <typename _T, typename _Index>
static void p5a_finish_(eig_result<_T>& result,
                        const eig_solver_method actual_method,
                        const std::string& used_method)
{
    const std::string kernel_failure = result.failure_reason;
    set_eig_diagnostics_<_T,_Index>(result, actual_method, result.used_subspace_dim,
        result.matrix_vector_products, result.breakdown_reason, kernel_failure);
    result.used_method = used_method;
}

// Default (direct) entry.  Returns false when the caller must continue with
// its present code (Krylov-Schur on the sparse LU; also reports an LU that
// cannot be factorized).
template <typename _T, typename _Index>
static bool p5a_gsi_direct_(const spmats<_T,_Index>& self,
                            const spmats<_T,_Index>& B,
                            const std::size_t k,
                            const eig_options<_T>& options,
                            const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                            const eig_solver_method actual_method,
                            eig_result<_T>& result)
{
    typedef typename std::integral_constant<bool, t2_supported_<_T>::value>::type tag_t;
    const std::string used = "shift_invert_arnoldi(promoted_from_lanczos)";
    p5a_sym_factor_<_T,_Index> fac;
    const spmats<_T,_Index> Bc = B.as_csr();
    {
        const spmats<_T,_Index> K = vcp::tsparse::subtract_scaled_sparse(self, B, _T(sigma));
        fac.factorize(K);
    }
    if (fac.kind != 0) {
        std::size_t solves = 0;
        struct SymApply {
            const p5a_sym_factor_<_T,_Index>* f;
            const spmats<_T,_Index>* Bm;
            std::size_t* count;
            void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
                y = Bm->mul_vec(x);
                f->solve(y);
                (*count)++;
            }
        } apply_sym = { &fac, &Bc, &solves };
        if (!p5a_gsi_lanczos_core_<_T,_Index>(self, B, k, options, sigma, apply_sym, result, tag_t()))
            return false;
        result.linear_solves = solves;
        result.factorization_diagnostics = fac.diagnostics;
        p5a_finish_<_T,_Index>(result, actual_method, used);
        return true;
    }
    // (c) the sparse LU as the solver of the same iteration
    typedef vcp::tsparse::lu_shift_invert_operator<spmats<_T,_Index> > LUOp;
    LUOp op(self, B, _T(sigma), options.shift_invert_lu);
    if (!op.factorization_ok()) return false;
    struct LUApply {
        LUOp* op;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
    } apply_lu = { &op };
    if (!p5a_gsi_lanczos_core_<_T,_Index>(self, B, k, options, sigma, apply_lu, result, tag_t()))
        return false;
    result.linear_solves = op.linear_solves();
    result.inner_iterations = op.inner_iterations();
    result.inner_failure_count = op.inner_failure_count();
    result.inner_residual_norm = op.inner_residual_norm();
    result.factorization_diagnostics = "shift_invert_factorization=sparse_lu ("
        + fac.diagnostics + "): " + op.factorization_diagnostics();
    result.factorization_zero_pivots = op.factorization_zero_pivots();
    p5a_finish_<_T,_Index>(result, actual_method, used);
    return true;
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
// EIG-T1 T-2a: the T-1 kernel on (A - sigma I)^{-1}; no lambda-space gate,
// no end-of-run demotion, no polish, no Krylov-Schur delegation (T-9b).
template <typename _T, typename _Index, class Apply>
static eig_result<_T> shift_invert_lanczos_drive_(const spmats<_T,_Index>& self,
                                                    const std::size_t k,
                                                    const eig_options<_T>& options,
                                                    const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                                    Apply& apply_si)
{
    return t2_si_lanczos_core_<_T,_Index>(self, k, options, sigma, apply_si,
        typename std::integral_constant<bool, t2_supported_<_T>::value>::type());
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
static eig_result<_T>
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
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// --- E-A1 sparse_lu path: standard shift-invert Arnoldi --------------------
template <typename _T, typename _Index>
static eig_result<_T>
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
    // EIG-T1 P5b: no end-of-run residual gate and no polish (T-9, every scalar type)
    result.linear_solves = op.linear_solves();
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

// --- E-A1 sparse_lu path: generalized shift-invert -------------------------
template <typename _T, typename _Index>
static eig_result<_T>
generalized_shift_invert_sparse_lu_(const spmats<_T,_Index>& self,
                                      const spmats<_T,_Index>& B,
                                      const std::size_t k,
                                      const eig_options<_T>& options,
                                      const typename vcp::tsparse_scalar::real_type<_T>::type& sigma,
                                      const eig_solver_method actual_method,
                                      const bool promoted_from_lanczos)
{
    // EIG-T1 P5a: symmetric pair + Lanczos-type request -> the T-1 kernel in
    // the B inner product on Cholesky / LDL^T / sparse LU (in this order)
    if (p5a_route_<_T,_Index>(self, B, options, promoted_from_lanczos)) {
        eig_result<_T> r5;
        if (p5a_gsi_direct_<_T,_Index>(self, B, k, options, sigma, actual_method, r5)) {
            set_result_counts_<_T,_Index>(r5, k);
            return r5;
        }
    }
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
    // EIG-T1 P5b: no end-of-run residual gate and no polish (T-9, every scalar type)
    result.linear_solves = op.linear_solves();
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
                av_fn, pv_fn, x, (inner_max / 10 + 1) * std::min(nn, inner_max), inner_tol, std::min(nn, inner_max));
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
    const scalar_real_type shift_val = options.shift;

    // EIG-T1 (G-1 ruling E): no Lanczos-type route for interval scalars
    if (!t2_supported_<_T>::value) {
        const eig_solver_method m_e =
            (options.method == eig_solver_method::shift_invert_lanczos)
                ? eig_solver_method::shift_invert_lanczos : eig_solver_method::lanczos;
        eig_result<_T> r_e = t2_unsupported_scalar_<_T,_Index>(k, m_e, eig_method_to_string_<_T,_Index>(m_e));
        r_e.used_shift_invert = (m_e == eig_solver_method::shift_invert_lanczos);
        return r_e;
    }

    if (options.method == eig_solver_method::shift_invert_lanczos) {
        return shift_invert_lanczos_eigs_<_T,_Index>(self, k, options, shift_val);
    }
    // EIG-T1 T-2c: the T-1 kernel on A (method / used_method stay "lanczos")
    return t2_lanczos_plain_<_T,_Index>(self, k, options,
        typename std::integral_constant<bool, t2_supported_<_T>::value>::type());
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
                av_fn, pv_fn, x, (inner_max / 10 + 1) * std::min(nn, inner_max), inner_tol, std::min(nn, inner_max));
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
    // EIG-T1 P5b: no end-of-run residual gate and no polish (T-9, every scalar type)
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
    // Fence for the scalar types outside trl_is_real_floating, the trait this
    // real route is selected by (trl_eigs_new_).  The kernel's own
    // static_assert is trl_is_supported_scalar_ (EIG-T2: real point scalars
    // with a numeric_limits specialization and their std::complex).
    // What actually arrives here: scalars without a numeric_limits
    // specialization (kv::interval).  float / double / long double / kv::dd /
    // kv::mpfr take the true overload (EIG-T1 T-1f).  Complex scalars never
    // arrive: policy_eigs_with_info sends them to
    // complex_standard_eigs_with_info_ first, which refuses an explicit
    // method=thick_restart_lanczos itself.  Honest explicit refusal (C-5);
    // the auto_select routing never reaches this (it sends such scalars to
    // the KS side), only explicit method=thick_restart_lanczos.
    eig_result<_T> result;
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    result.failure_reason =
        "thick_restart_lanczos: scalar type unsupported (real point"
        " scalars only); use krylov_schur or lanczos";
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
    // sparse A (diagnostics; EIG-T1 P5b: they do not demote the verdict of
    // the dense solver, for any scalar type)
    if (!result.eigenvectors.empty()
        && result.eigenvectors.size() == result.eigenvalues.size()) {
        result.residuals_absolute = eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
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
        // EIG-T1 T-9b: without the end-of-run residual check and the
        // Krylov-Schur delegation behind it, the symmetric driver must only
        // receive symmetric input (as dsaupd / dnaupd are chosen in ARPACK).
        // The same symmetry decision as the no-shift routing below; a matrix
        // that is not certainly symmetric goes to shift_invert_arnoldi.
        // (G-1 ruling E) scalars without a Lanczos-type route also go to the
        // Krylov-Schur side: auto never routes to a refusal because of the type.
        {
            bool sym_si;
            if (options.structure == matrix_structure_hint::symmetric)     sym_si = true;
            else if (options.structure == matrix_structure_hint::general)  sym_si = false;
            else sym_si = is_certainly_symmetric_<_T,_Index>(A);
            if (!t2_supported_<_T>::value) sym_si = false;
            if (!sym_si) {
                active.method = eig_solver_method::shift_invert_arnoldi;
                eig_result<_T> result = arnoldi_eigs_new_<_T,_Index>(A, k, active);
                result.method = eig_solver_method::auto_select;
                result.used_method = "auto_select(" + result.used_method + ")";
                return result;
            }
        }
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

    // B-28 と同方向の保守則: TRL の実の経路が対応しないスカラー型
    // (trl_is_real_floating が false = numeric_limits の特殊化が無い型、
    // kv::interval)は対称でも KS 側に倒す(KS は対称入力でも正しく、
    // 汎用スカラーで動く。auto が型を理由に拒否結果へルートすることは
    // ない)。float / double / long double / kv::dd / kv::mpfr は
    // trl_is_real_floating が true で、対称なら核(TRL)に行く(EIG-T1 T-1f)。
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

    // EIG-T1 P5a: symmetric pair + Lanczos-type request -> the T-1 kernel in
    // the B inner product; the inner solver (ILU(0) + GMRES) is unchanged
    if (p5a_route_<_T,_Index>(self, B, options, promoted_from_lanczos)) {
        struct ApplyP5 {
            GSIOperator* op;
            void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { op->apply(x, y); }
        } apply_p5 = { &gsi_op };
        eig_result<_T> r5;
        if (p5a_gsi_lanczos_core_<_T,_Index>(self, B, k, options, sigma, apply_p5, r5,
                typename std::integral_constant<bool, t2_supported_<_T>::value>::type())) {
            r5.linear_solves = gsi_op.linear_solves();
            r5.inner_iterations = gsi_op.inner_iterations();
            r5.inner_failure_count = gsi_op.inner_failure_count();
            r5.inner_residual_norm = gsi_op.inner_residual_norm();
            r5.factorization_diagnostics = gsi_op.factorization_diagnostics();
            r5.factorization_zero_pivots = gsi_op.factorization_zero_pivots();
            const bool inner_ok5 = (gsi_op.inner_failure_count() == 0);
            r5.converged = r5.converged && inner_ok5;
            if (!inner_ok5) {
                r5.method = actual_method;
                r5.used_method = "shift_invert_arnoldi(promoted_from_lanczos)";
                r5.status = "inner_solve_failed";
                r5.failure_reason = "inner GMRES solve failed in generalized shift-invert";
                r5.inner_failure_reason = r5.failure_reason;
                r5.message = r5.failure_reason;
            } else {
                p5a_finish_<_T,_Index>(r5, actual_method, "shift_invert_arnoldi(promoted_from_lanczos)");
            }
            set_result_counts_<_T,_Index>(r5, k);
            return r5;
        }
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
    // EIG-T1 P5b: no end-of-run residual gate and no polish (T-9, every scalar type)
    result.linear_solves = gsi_op.linear_solves();
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

    // EIG-T1 T-2d: sparse Cholesky of B and the T-1 kernel on
    // L^{-1} P A P^T L^{-T} (used_method stays "b_inner_lanczos").
    if (t2_generalized_chol_<_T,_Index>(self, B, k, options, result,
            typename std::integral_constant<bool, t2_supported_<_T>::value>::type())) {
        return result;
    }
    // EIG-T1 (G-1 ruling E): no generalized no-shift route for interval scalars
    eig_result<_T> r_e = t2_unsupported_scalar_<_T,_Index>(k, eig_solver_method::lanczos, "b_inner_lanczos");
    r_e.used_generalized_operator = true;
    return r_e;
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
    // EIG-T1 (G-1 ruling E): no Lanczos-type route for interval scalars
    if (!t2_supported_<_T>::value) {
        eig_result<_T> r_e = t2_unsupported_scalar_<_T,_Index>(k,
            eig_solver_method::shift_invert_lanczos, "shift_invert_lanczos+user_preconditioner");
        r_e.used_shift_invert = true;
        return r_e;
    }
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
    const scalar_real_type inner_tol = options.tol / scalar_real_type(10);
    // 相対(初回前処理後残差 β₀ 基準)。係数 1/10 の根拠は SPU-M1 設計書の実測。
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
                av_fn, pv_fn, x, (imax / 10 + 1) * std::min(nn, imax), itol, std::min(nn, imax), true);
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

    // EIG-T1 T-2b: the T-1 kernel on the preconditioned inner solve (same
    // core as T-2a)
    eig_result<_T> result = t2_si_lanczos_core_<_T,_Index>(self, k, options, sigma, apply_si,
        typename std::integral_constant<bool, t2_supported_<_T>::value>::type());
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
    // 内側の消費(求解の回数・反復数・残差・失敗の数)を共有カウンタから同期する
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
    const scalar_real_type inner_tol = options.tol / scalar_real_type(10);
    // 相対(初回前処理後残差 β₀ 基準)。係数 1/10 の根拠は SPU-M1 設計書の実測。
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
                av_fn, pv_fn, x, (imax / 10 + 1) * std::min(nn, imax), itol, std::min(nn, imax), true);
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
    // EIG-T1 P5b: no end-of-run residual gate and no polish (T-9, every scalar type)
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
    const scalar_real_type inner_tol = options.tol / scalar_real_type(10);
    // 相対(初回前処理後残差 β₀ 基準)。係数 1/10 の根拠は SPU-M1 設計書の実測。
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
            GR gr = vcp::tsparse_factorization::gmres_solve<_T, AV, PV>(av_fn, pv_fn, rhs, act_max * std::min(act_restart, nn), itol, act_restart, true);
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

    // EIG-T1 P5a: symmetric pair + Lanczos-type request -> the T-1 kernel in
    // the B inner product; the inner solver (preconditioned GMRES) is unchanged
    if (p5a_route_<_T,_Index>(self, B, options, promoted_from_lanczos)) {
        eig_result<_T> r5;
        if (p5a_gsi_lanczos_core_<_T,_Index>(self, B, k, options, sigma, apply_gsi, r5,
                typename std::integral_constant<bool, t2_supported_<_T>::value>::type())) {
            r5.linear_solves = linear_solve_count;
            r5.inner_iterations = inner_iteration_count;
            r5.inner_failure_count = inner_failure_count;
            r5.inner_residual_norm = inner_residual_norm;
            r5.factorization_diagnostics = vcp::tsparse_prec_traits::get_diagnostics(M);
            r5.factorization_zero_pivots = vcp::tsparse_prec_traits::get_zero_pivots(M);
            if (inner_failure_count != 0) {
                r5.method = actual_method;
                r5.used_method = used_method_str;
                r5.converged = false;
                r5.eigenvalues.clear();
                r5.eigenvectors.clear();
                r5.status = "inner_solve_failed";
                r5.failure_reason = prec_failure_reason.empty()
                    ? "inner GMRES solve failed in generalized shift-invert with user preconditioner"
                    : prec_failure_reason;
                r5.inner_failure_reason = r5.failure_reason;
                r5.message = r5.failure_reason;
            } else {
                p5a_finish_<_T,_Index>(r5, actual_method, used_method_str);
            }
            set_result_counts_<_T,_Index>(r5, k);
            return r5;
        }
    }

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
    // EIG-T1 P5b: no end-of-run residual gate and no polish (T-9, every scalar type)
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

// ---------------------------------------------------------------------------
// EIG-T1 T-7c: entry normalization helpers.
//
// The point entries (standard / generalized, with or without a user
// preconditioner; real point scalars only -- the complex Hermitian route
// applies the same scaling inside hermitian_trl_route_, EIG-T2) solve A / s
// with s a power of two
// near ||A||_inf and map the result back, so that every path behind the
// entry (thick-restart, Krylov-Schur, dense, shift-invert) works on a matrix
// of norm about 1.  Powers of two do not round.  ||A||_inf = 0 and
// non-finite norms are not normalized.
// ---------------------------------------------------------------------------
template <typename _T>
struct t7c_enabled_ {
    static const bool value =
        !spmatrix_is_complex<_T>::value
        && !vcp::tsparse_scalar::numeric_limits_is_unspecialized<
               typename vcp::tsparse_scalar::real_type<_T>::type>::value;
};

// s = 2^e with s <= norm < 2 s (norm > 0); returns 1 (e = 0) otherwise.
template <typename R>
static R t7c_pow2_floor_(const R& norm, int& e)
{
    e = 0;
    R s = R(1);
    if (!(norm > R(0))) return s;
    int guard = 0;
    while (s * R(2) <= norm && guard < 5000) { s = s * R(2); e++; guard++; }
    while (s > norm && guard < 10000) { s = s / R(2); e--; guard++; }
    if (!(s <= norm && norm < s * R(2))) { e = 0; return R(1); }
    return s;
}

// 2^round(log2 norm): the power of two nearest in the log scale.
template <typename R>
static R t7c_pow2_scale_(const R& norm)
{
    int e = 0;
    R s = t7c_pow2_floor_(norm, e);
    if (!(norm > R(0))) return s;
    const R q = norm / s;
    if (q * q > R(2)) s = s * R(2);
    return s;
}

// A power of FOUR s with s / 2 <= norm < 2 s; root = sqrt(s) (a power of
// two), so that B-normalized eigenvectors map back without rounding.
template <typename R>
static R t7c_pow4_scale_(const R& norm, R& root)
{
    int e = 0;
    R s = t7c_pow2_floor_(norm, e);
    root = R(1);
    if (!(norm > R(0))) return s;
    if (e % 2 != 0) { s = s * R(2); e++; }
    for (int i = 0; i < e / 2; i++) root = root * R(2);
    for (int i = 0; i > e / 2; i--) root = root / R(2);
    return s;
}

template <typename _T, typename _Index>
static spmats<_T,_Index> t7c_scaled_copy_(const spmats<_T,_Index>& A,
    const typename vcp::tsparse_scalar::real_type<_T>::type& s)
{
    spmats<_T,_Index> C = A.as_csr();
    std::vector<_T>& v = C.values();
    for (std::size_t i = 0; i < v.size(); i++) v[i] = v[i] / _T(s);
    return C;
}

// eigenvalues x lam_factor, absolute residuals x res_factor, eigenvectors x
// vec_factor (relative residuals are dimensionless and stay as computed).
template <typename _T>
static void t7c_unscale_result_(eig_result<_T>& r,
    const typename vcp::tsparse_scalar::real_type<_T>::type& lam_factor,
    const typename vcp::tsparse_scalar::real_type<_T>::type& res_factor,
    const typename vcp::tsparse_scalar::real_type<_T>::type& vec_factor)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
    for (std::size_t i = 0; i < r.eigenvalues.size(); i++)
        r.eigenvalues[i] = r.eigenvalues[i] * _T(lam_factor);
    for (std::size_t i = 0; i < r.eigenvalues_imag.size(); i++)
        r.eigenvalues_imag[i] = r.eigenvalues_imag[i] * _T(lam_factor);
    for (std::size_t i = 0; i < r.complex_eigenvalues.size(); i++)
        r.complex_eigenvalues[i] = typename eig_result<_T>::eigenvalue_type(
            r.complex_eigenvalues[i].real() * lam_factor,
            r.complex_eigenvalues[i].imag() * lam_factor);
    for (std::size_t i = 0; i < r.residuals_absolute.size(); i++)
        r.residuals_absolute[i] = r.residuals_absolute[i] * res_factor;
    for (std::size_t i = 0; i < r.residual_history_absolute.size(); i++)
        r.residual_history_absolute[i] = r.residual_history_absolute[i] * res_factor;
    r.residual_norm_absolute = r.residual_norm_absolute * res_factor;
    if (vec_factor != R(1)) {
        for (std::size_t i = 0; i < r.eigenvectors.size(); i++)
            for (std::size_t j = 0; j < r.eigenvectors[i].size(); j++)
                r.eigenvectors[i][j] = r.eigenvectors[i][j] * _T(vec_factor);
    }
}

// ---------------------------------------------------------------------------
// EIG-T2 D-3: hermitian_lanczos_eigs_impl_ (complex Hermitian, standard, no
// shift) on the thick-restart kernel.  The steps are those of the real
// symmetric route (policy_eigs_with_info -> lanczos_eigs_new_ ->
// t2_lanczos_plain_), in the same order:
//   1. T-7c: solve A / s with s = 2^round(log2 ||A||_inf), shift / s
//   2. the kernel without its end-of-run applications, on the operator A / s,
//      with the caller's options (max_iter = total number of applications;
//      tol, subspace_dim, target, shift, seed, history)
//   3. true residuals of the returned pairs as diagnostics (Hermitian norm);
//      they do not decide `converged`
//   4. t2_finish_lanczos_type_ (status "converged" / "not_converged",
//      method = lanczos, counts); used_method stays "hermitian_lanczos"
//   5. map back: eigenvalues and absolute residuals x s
// The verdict is the kernel's own; no acceptance, demotion or rescue by the
// true residuals.
// Tag: complex scalars the kernel accepts; the false overload is a fence for
// every other scalar type (complex_standard_eigs_with_info_ is instantiated
// for real scalars too, but never executed for them).
// ---------------------------------------------------------------------------
template <typename _T>
struct hermitian_trl_supported_ {
    static const bool value =
        spmatrix_is_complex<_T>::value
        && vcp::tsparse_experimental::trl_is_supported_scalar_<_T>::value;
};

template <typename _T, typename _Index>
static eig_result<_T> hermitian_trl_route_(const spmats<_T,_Index>&,
                                           const std::size_t k,
                                           const eig_options<_T>& options,
                                           std::false_type /* hermitian_trl_supported */)
{
    eig_result<_T> result;   // fence
    result.requested_count = k;
    result.converged = false;
    result.status = "unsupported";
    result.failure_reason =
        "hermitian_lanczos: scalar type unsupported (std::complex of a real"
        " point scalar only)";
    result.message = result.failure_reason;
    result.method = options.method;
    result.used_method = "hermitian_lanczos";
    set_result_counts_<_T,_Index>(result, k);
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> hermitian_trl_route_(const spmats<_T,_Index>& self,
                                           const std::size_t k,
                                           const eig_options<_T>& options,
                                           std::true_type /* hermitian_trl_supported */)
{
    typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
    spmats<_T,_Index> A = self.as_csr();
    const std::size_t n = static_cast<std::size_t>(self.rowsize());
    // 1. T-7c (as the real entry; a power of two does not round)
    const R s7 = t7c_pow2_scale_<R>(matrix_inf_norm_value_<_T,_Index>(A));
    eig_options<_T> o7 = options;
    if (s7 != R(1)) {
        std::vector<_T>& v = A.values();
        for (std::size_t i = 0; i < v.size(); i++) v[i] = v[i] / s7;
        o7.shift = options.shift / s7;
    }
    // 2. the kernel (Hermitian inner product, EIG-T2 D-1)
    struct apply_fn {
        const spmats<_T,_Index>* mat;
        void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
    };
    apply_fn apply_op = { &A };
    eig_result<_T> result =
        vcp::tsparse_experimental::thick_restart_lanczos_eigs_no_residual_<apply_fn, _T>(
            apply_op, n, k, o7);
    // 3. diagnostics
    if (!result.eigenvectors.empty()) {
        result.residuals_absolute = hermitian_eigenpair_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
        result.residuals_relative = hermitian_eigenpair_relative_residuals_<_T,_Index>(A, result.eigenvalues, result.eigenvectors);
    }
    // 4. status / counts
    t2_finish_lanczos_type_<_T,_Index>(result, eig_solver_method::lanczos, k);
    result.used_method = "hermitian_lanczos";
    result.used_shift_invert = false;
    // 5. map back
    if (s7 != R(1)) t7c_unscale_result_<_T>(result, s7, s7, R(1));
    return result;
}

template <typename _T, typename _Index>
static eig_result<_T> hermitian_lanczos_eigs_impl_(const spmats<_T,_Index>& self,
                                                     const std::size_t k,
                                                     const eig_options<_T>& options)
{
    return hermitian_trl_route_<_T,_Index>(self, k, options,
        typename std::integral_constant<bool, hermitian_trl_supported_<_T>::value>::type());
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
    // EIG-T1 T-7c: solve A / s (s = 2^round(log2 ||A||_inf)), shift sigma / s,
    // and map back.  The call below re-enters this function with a matrix
    // whose scale is exactly 1, so it does not normalize again.
    if (t7c_enabled_<_T>::value) {
        const scalar_real_type s7 =
            t7c_pow2_scale_<scalar_real_type>(matrix_inf_norm_value_<_T,_Index>(A));
        if (s7 != scalar_real_type(1)) {
            const spmats<_T,_Index> A7 = t7c_scaled_copy_<_T,_Index>(A, s7);
            eig_options<_T> o7 = options;
            o7.shift = options.shift / s7;
            eig_result<_T> r7 = policy_eigs_with_info<_T,_Index>(A7, k, o7);
            t7c_unscale_result_<_T>(r7, s7, s7, scalar_real_type(1));
            return r7;
        }
    }
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
    // EIG-T1 T-7c: as in the overload without preconditioner.  M is passed
    // unchanged: the inner systems (A - sigma I) / s y = x have the same
    // residuals and the same Krylov spaces as the unscaled ones.
    if (t7c_enabled_<_T>::value) {
        const scalar_real_type s7 =
            t7c_pow2_scale_<scalar_real_type>(matrix_inf_norm_value_<_T,_Index>(A));
        if (s7 != scalar_real_type(1)) {
            const spmats<_T,_Index> A7 = t7c_scaled_copy_<_T,_Index>(A, s7);
            eig_options<_T> o7 = options;
            o7.shift = options.shift / s7;
            eig_result<_T> r7 = policy_eigs_with_info<_T,_Index,Preconditioner>(A7, k, o7, M);
            t7c_unscale_result_<_T>(r7, s7, s7, scalar_real_type(1));
            return r7;
        }
    }
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

    bool use_lanczos =
        (active.method == eig_solver_method::shift_invert_lanczos)
        || (active.use_shift && active.method == eig_solver_method::lanczos)
        || (active.use_shift && active.method == eig_solver_method::auto_select);
    // EIG-T1 T-9b: auto_select hands the symmetric driver symmetric input only
    // (same decision as auto_select_eigs_).
    if (active.use_shift && active.method == eig_solver_method::auto_select) {
        if (active.structure == matrix_structure_hint::general) use_lanczos = false;
        else if (active.structure != matrix_structure_hint::symmetric)
            use_lanczos = is_certainly_symmetric_<_T,_Index>(A);
    }
    // (G-1 ruling E) auto_select: scalars without a Lanczos-type route go to the Krylov-Schur side
    if (!t2_supported_<_T>::value
        && active.use_shift && active.method == eig_solver_method::auto_select) use_lanczos = false;

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
    // EIG-T1 T-7c: solve (A / sA, B / sB) with sA a power of two and sB a
    // power of four; eigenvalues x sA / sB, shift x sB / sA.  Vectors of the
    // B-normalized route (used_method "b_inner_lanczos") x 1 / sqrt(sB), so
    // that x^T B x = 1 still holds; absolute residuals follow.
    if (t7c_enabled_<_T>::value) {
        scalar_real_type rootB7(1);
        const scalar_real_type sA7 =
            t7c_pow2_scale_<scalar_real_type>(matrix_inf_norm_value_<_T,_Index>(A));
        const scalar_real_type sB7 =
            t7c_pow4_scale_<scalar_real_type>(matrix_inf_norm_value_<_T,_Index>(B), rootB7);
        if (sA7 != scalar_real_type(1) || sB7 != scalar_real_type(1)) {
            const spmats<_T,_Index> A7 = t7c_scaled_copy_<_T,_Index>(A, sA7);
            const spmats<_T,_Index> B7 = t7c_scaled_copy_<_T,_Index>(B, sB7);
            eig_options<_T> o7 = options_in;
            o7.shift = options_in.shift * sB7 / sA7;
            eig_result<_T> r7 = policy_generalized_eigs_with_info<_T,_Index>(A7, B7, k, o7);
            const bool bnorm7 = (r7.used_method == "b_inner_lanczos");
            const scalar_real_type vf7 = bnorm7 ? scalar_real_type(scalar_real_type(1) / rootB7) : scalar_real_type(1);
            t7c_unscale_result_<_T>(r7, sA7 / sB7, sA7 * vf7, vf7);
            return r7;
        }
    }

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
    // EIG-T1 T-7c: as in the overload without preconditioner (M unchanged:
    // (A - sigma B) / sA has the same inner residuals and Krylov spaces).
    if (t7c_enabled_<_T>::value) {
        scalar_real_type rootB7(1);
        const scalar_real_type sA7 =
            t7c_pow2_scale_<scalar_real_type>(matrix_inf_norm_value_<_T,_Index>(A));
        const scalar_real_type sB7 =
            t7c_pow4_scale_<scalar_real_type>(matrix_inf_norm_value_<_T,_Index>(B), rootB7);
        if (sA7 != scalar_real_type(1) || sB7 != scalar_real_type(1)) {
            const spmats<_T,_Index> A7 = t7c_scaled_copy_<_T,_Index>(A, sA7);
            const spmats<_T,_Index> B7 = t7c_scaled_copy_<_T,_Index>(B, sB7);
            eig_options<_T> o7 = options_in;
            o7.shift = options_in.shift * sB7 / sA7;
            eig_result<_T> r7 = policy_generalized_eigs_with_info<_T,_Index,Preconditioner>(A7, B7, k, o7, M);
            const bool bnorm7 = (r7.used_method == "b_inner_lanczos");
            const scalar_real_type vf7 = bnorm7 ? scalar_real_type(scalar_real_type(1) / rootB7) : scalar_real_type(1);
            t7c_unscale_result_<_T>(r7, sA7 / sB7, sA7 * vf7, vf7);
            return r7;
        }
    }

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
    return normalize_eig_result_counts_<_T>(policy_eigs_with_info_impl(k, opt));
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
    return normalize_eig_result_counts_<_T>(vcp::policy_eigs_with_info<_T,_Index,Prec>(A, k, opt, M));
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
    return normalize_eig_result_counts_<_T>(policy_generalized_eigs_with_info_impl(B, k, opt));
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
    return normalize_eig_result_counts_<_T>(vcp::policy_generalized_eigs_with_info<_T,_Index,Prec>(A, B, k, opt, M));
}

// ===========================================================================
// Strict policy methods: policy_eigs, policy_generalized_eigs (and
// preconditioner overloads).
// All convergence/count checking lives here, not in spmatrix.hpp.
//
// policy_eigs / policy_generalized_eigs (no Preconditioner) are NVI outers
// too: they finalize their spmats arguments and delegate to a virtual _impl,
// even though their _impl bodies only call back into policy_eigs_with_info /
// policy_generalized_eigs_with_info (already finalize-safe on their own).
// This keeps all four eigs entry points independently safe and overridable,
// per spmats_finalize_policy.md §3.
// ===========================================================================

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
