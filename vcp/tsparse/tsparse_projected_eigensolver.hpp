// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// vcp/tsparse/tsparse_projected_eigensolver.hpp
//
// Phase 8: internal helper for small dense projected eigenproblems.
//
// Provides:
//   - projected_eigenpair<T>          internal result type (NOT eig_result<T>)
//   - projected_eigensolver_result<T> internal result type
//   - solve_real_symmetric_projected  robust Jacobi solver for real symmetric H
//   - reconstruct_ritz_vector_safe    dimension-checked Ritz vector reconstruction
//   - select_projected_indices        deterministic target selection
//   - compute_projected_residuals     ||H y - theta y|| in type T (no double fallback)
//
// None of these types or functions are part of the public vcp::spmatrix API.
// Do NOT expose via eig_result<T>, eig_options<T>, or eig_method enum.

#pragma once

#ifndef VCP_TSPARSE_PROJECTED_EIGENSOLVER_HPP
#define VCP_TSPARSE_PROJECTED_EIGENSOLVER_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_projected {

// ---------------------------------------------------------------------------
// projected_eigenpair<T>
//
// Stores one eigenpair of a small projected dense problem.
// residual_estimate = ||H y - theta y|| computed in type real_type<T>.
// NOT related to eig_result::residuals_absolute.
// ---------------------------------------------------------------------------
template <class T>
struct projected_eigenpair {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
    T value;
    std::vector<T> vector;
    real_type residual_estimate;

    projected_eigenpair()
        : value(T(0)), residual_estimate(real_type(0)) {}
};

// ---------------------------------------------------------------------------
// projected_eigensolver_result<T>
//
// Stores all eigenpairs of a small projected dense problem.
// 'success' is false if the underlying solver did not converge.
// Partial results are still returned for diagnostic use.
// ---------------------------------------------------------------------------
template <class T>
struct projected_eigensolver_result {
    bool success;
    bool converged;
    std::string status;
    std::string message;
    std::vector<projected_eigenpair<T> > pairs;
    std::size_t iterations;

    projected_eigensolver_result()
        : success(false), converged(false), iterations(0) {}
};

// ---------------------------------------------------------------------------
// compute_projected_residuals
//
// For each (value, vector) pair: computes ||H y - theta y|| in real_type<T>.
// H and y must have matching sizes.
// Does NOT cast to double; uses T arithmetic throughout.
// ---------------------------------------------------------------------------
template <class T>
std::vector<typename vcp::tsparse_scalar::real_type<T>::type>
compute_projected_residuals(
    const std::vector<std::vector<T> >& H,
    const std::vector<T>& eigenvalues,
    const std::vector<std::vector<T> >& eigenvectors)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = H.size();
    const std::size_t m = eigenvalues.size();
    std::vector<R> residuals;
    residuals.reserve(m);
    for (std::size_t p = 0; p < m; p++) {
        if (p >= eigenvectors.size() || eigenvectors[p].size() != n) {
            residuals.push_back((std::numeric_limits<R>::infinity)());
            continue;
        }
        const std::vector<T>& y = eigenvectors[p];
        std::vector<T> hy(n, T(0));
        for (std::size_t row = 0; row < n; row++) {
            if (H[row].size() < n) {
                residuals.push_back((std::numeric_limits<R>::infinity)());
                goto next_pair;
            }
            for (std::size_t col = 0; col < n; col++) {
                hy[row] += H[row][col] * y[col];
            }
            hy[row] -= eigenvalues[p] * y[row];
        }
        residuals.push_back(vcp::tsparse_scalar::real_norm_value(hy));
        next_pair:;
    }
    return residuals;
}

// ---------------------------------------------------------------------------
// solve_real_symmetric_projected
//
// Solves the small real symmetric eigenproblem H y = theta y.
// H is n x n real symmetric; n is allowed to be 0 or 1.
//
// Guarantees:
//   - n == 0: empty success
//   - n == 1: single pair, residual == 0
//   - Repeated eigenvalues are preserved (no value-based dedup)
//   - Eigenvectors are normalized
//   - Projected residuals ||H y - theta y|| are computed in T / real_type<T>
//   - Double is NOT used; all arithmetic is in T / real_type<T>
//   - If Jacobi does not converge, success=false, status="jacobi_not_converged"
//     but partial pairs are still returned
//
// max_iter == 0 uses the default: max(100*n*n, 1000)
// ---------------------------------------------------------------------------
template <class T>
projected_eigensolver_result<T> solve_real_symmetric_projected(
    const std::vector<std::vector<T> >& H,
    const std::size_t max_iter,
    const typename vcp::tsparse_scalar::real_type<T>::type& tol)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    projected_eigensolver_result<T> result;
    const std::size_t n = H.size();

    // --- n == 0 ---
    if (n == 0) {
        result.success = true;
        result.converged = true;
        result.status = "success";
        result.message = "empty projected problem (n=0)";
        result.iterations = 0;
        return result;
    }

    // --- n == 1 ---
    if (n == 1) {
        projected_eigenpair<T> pair;
        pair.value = H[0][0];
        pair.vector.assign(1, T(1));
        pair.residual_estimate = R(0);
        result.pairs.push_back(pair);
        result.success = true;
        result.converged = true;
        result.status = "success";
        result.message = "trivial 1x1 projected problem";
        result.iterations = 0;
        return result;
    }

    // --- General case: Jacobi eigensolver ---
    const std::size_t actual_max_iter = (max_iter == 0)
        ? std::max(n * n * std::size_t(100), std::size_t(1000))
        : max_iter;

    vcp::tsparse_dense_linalg::dense_eigen_result<T> dense =
        vcp::tsparse_dense_linalg::jacobi_eig_dense(H, actual_max_iter, tol);

    result.converged = dense.converged;
    result.iterations = dense.iterations;

    if (!dense.converged) {
        result.success = false;
        result.status = "jacobi_not_converged";
        result.message = "Jacobi projected eigensolver did not converge within max_iter sweeps; "
                         "partial eigenpairs are returned for diagnostic use";
    } else {
        result.success = true;
        result.status = "success";
    }

    // Compute projected residuals ||H y - theta y|| in T / real_type<T>
    // Uses compute_projected_residuals to avoid code duplication
    std::vector<R> res_vec;
    if (dense.eigenvalues.size() == dense.eigenvectors.size() && !dense.eigenvalues.empty()) {
        res_vec = compute_projected_residuals(H, dense.eigenvalues, dense.eigenvectors);
    }

    const std::size_t m = dense.eigenvalues.size();
    result.pairs.resize(m);
    for (std::size_t i = 0; i < m; i++) {
        result.pairs[i].value = dense.eigenvalues[i];

        if (i < dense.eigenvectors.size()) {
            result.pairs[i].vector = dense.eigenvectors[i];
            // Ensure unit norm (Jacobi produces unit-norm vectors but be defensive)
            const R nrm = vcp::tsparse_scalar::real_norm_value(result.pairs[i].vector);
            if (nrm > R(0)) {
                const R inv_nrm = R(1) / nrm;
                if (inv_nrm != R(1)) {
                    for (std::size_t j = 0; j < result.pairs[i].vector.size(); j++) {
                        result.pairs[i].vector[j] *= T(inv_nrm);
                    }
                }
            }
        } else {
            result.pairs[i].vector.assign(n, T(0));
        }

        result.pairs[i].residual_estimate = (i < res_vec.size())
            ? res_vec[i]
            : (std::numeric_limits<R>::infinity)();
    }

    return result;
}

// ---------------------------------------------------------------------------
// reconstruct_ritz_vector_safe<T, Coeff>
//
// Computes v = sum_{j=0}^{m-1}  V[j] * T(y[j])
// where V contains m basis vectors of dimension n, and y is the coefficient
// vector of dimension m.
//
// Coeff is typically T (real Lanczos) or real_type<T> (Hermitian Lanczos).
//
// Throws vcp::dimension_error if:
//   - y.size() != V.size()
//   - any V[j].size() != n
//
// Throws vcp::numerical_error if the resulting vector has norm <= zero_tol.
//
// On success, returns the normalized Ritz vector.
// ---------------------------------------------------------------------------
template <class T, class Coeff>
std::vector<T> reconstruct_ritz_vector_safe(
    const std::vector<std::vector<T> >& V,
    const std::vector<Coeff>& y,
    const std::size_t n,
    const typename vcp::tsparse_scalar::real_type<T>::type& zero_tol)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    if (y.size() != V.size()) {
        vcp::throw_error<vcp::dimension_error>(
            "tsparse_projected::reconstruct_ritz_vector_safe: "
            "coefficient size does not match basis size");
    }

    std::vector<T> v(n, T(0));
    for (std::size_t j = 0; j < V.size(); j++) {
        if (V[j].size() != n) {
            vcp::throw_error<vcp::dimension_error>(
                "tsparse_projected::reconstruct_ritz_vector_safe: "
                "basis vector dimension does not match n");
        }
        const T coeff = T(y[j]);
        for (std::size_t i = 0; i < n; i++) {
            v[i] += V[j][i] * coeff;
        }
    }

    const R nrm = vcp::tsparse_scalar::real_norm_value(v);
    if (nrm <= zero_tol) {
        vcp::throw_error<vcp::numerical_error>(
            "tsparse_projected::reconstruct_ritz_vector_safe: "
            "Ritz vector has near-zero norm");
    }

    const R inv_nrm = R(1) / nrm;
    for (std::size_t i = 0; i < n; i++) v[i] *= T(inv_nrm);
    return v;
}

// Convenience: T = Coeff specialization (most common case)
template <class T>
std::vector<T> reconstruct_ritz_vector_safe(
    const std::vector<std::vector<T> >& V,
    const std::vector<T>& y,
    const std::size_t n,
    const typename vcp::tsparse_scalar::real_type<T>::type& zero_tol)
{
    return reconstruct_ritz_vector_safe<T, T>(V, y, n, zero_tol);
}

// ---------------------------------------------------------------------------
// select_projected_indices
//
// Returns indices into pairs[] sorted by target (best first).
// k == 0            → empty vector
// k > pairs.size()  → clamped to pairs.size() (all indices returned)
// Repeated eigenvalues are NOT deduped; tie-breaking is deterministic
// (stable by original index).
// ---------------------------------------------------------------------------
template <class T>
std::vector<std::size_t> select_projected_indices(
    const std::vector<projected_eigenpair<T> >& pairs,
    const std::size_t k,
    const eig_target target,
    const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    if (k == 0 || pairs.empty()) return std::vector<std::size_t>();

    std::vector<std::complex<R> > ceigs;
    ceigs.reserve(pairs.size());
    for (std::size_t i = 0; i < pairs.size(); i++) {
        ceigs.push_back(std::complex<R>(
            vcp::tsparse_scalar::real_part(pairs[i].value), R(0)));
    }

    const std::size_t k_clamped = std::min(k, pairs.size());
    return vcp::tsparse_eigen_selection::select_eigen_indices(ceigs, k_clamped, target, shift);
}

// Variant operating on raw real eigenvalue arrays (for use inside Krylov solvers)
template <class T>
std::vector<std::size_t> select_projected_real_indices(
    const std::vector<T>& eigenvalues,
    const std::size_t k,
    const eig_target target,
    const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    if (k == 0 || eigenvalues.empty()) return std::vector<std::size_t>();
    const std::size_t k_clamped = std::min(k, eigenvalues.size());

    std::vector<std::complex<R> > ceigs;
    ceigs.reserve(eigenvalues.size());
    for (std::size_t i = 0; i < eigenvalues.size(); i++) {
        ceigs.push_back(std::complex<R>(
            vcp::tsparse_scalar::real_part(eigenvalues[i]), R(0)));
    }
    return vcp::tsparse_eigen_selection::select_eigen_indices(ceigs, k_clamped, target, shift);
}

} // namespace tsparse_projected
} // namespace vcp

#endif // VCP_TSPARSE_PROJECTED_EIGENSOLVER_HPP
