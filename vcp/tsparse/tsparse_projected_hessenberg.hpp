// VCP Library
// http://verified.computation.jp
//
// vcp/tsparse/tsparse_projected_hessenberg.hpp
//
// Phase 3 helper: dense Hessenberg reduction for projected eigenproblem.
// Used by tsparse_krylov_schur.hpp when the restart compressed matrix
// is not already in upper Hessenberg form.
//
// Namespace: vcp::tsparse_proj_hess
//
// NOT part of the public vcp::spmatrix API.

#pragma once

#ifndef VCP_TSPARSE_PROJECTED_HESSENBERG_HPP
#define VCP_TSPARSE_PROJECTED_HESSENBERG_HPP

#include <cstddef>
#include <limits>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_proj_hess {

// Reduce a general n x n real matrix A to upper Hessenberg form in-place
// using Householder reflections.  The transformation is NOT accumulated.
// Only the upper Hessenberg portion of A is meaningful after this call.
// No-op for n <= 2 (already Hessenberg).
template <typename T>
void reduce_to_hessenberg_inplace(std::vector<std::vector<T> >& A)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = A.size();
    if (n <= 2) return;

    for (std::size_t k = 0; k + 2 <= n; k++) {
        // Build Householder vector u to zero out A[k+2 .. n-1][k].
        // The vector lives in rows k+1 .. n-1 (length len = n-k-1).
        const std::size_t len = n - k - 1;

        R sigma = R(0);
        for (std::size_t i = k + 1; i < n; i++) {
            const R xi = vcp::tsparse_scalar::real_part(A[i][k]);
            sigma += xi * xi;
        }
        sigma = vcp::tsparse_scalar::sqrt_value(sigma);

        const R eps10 = vcp::tsparse_scalar::epsilon<R>() * R(10);
        if (sigma <= eps10) continue;

        const R x0 = vcp::tsparse_scalar::real_part(A[k + 1][k]);
        const R sign_x0 = (x0 >= R(0)) ? R(1) : R(-1);
        const R u0 = x0 + sign_x0 * sigma;

        std::vector<R> u(len, R(0));
        u[0] = u0;
        for (std::size_t i = k + 2; i < n; i++) {
            u[i - k - 1] = vcp::tsparse_scalar::real_part(A[i][k]);
        }

        R u_norm2 = R(0);
        for (std::size_t i = 0; i < len; i++) u_norm2 += u[i] * u[i];
        if (u_norm2 <= eps10 * eps10) continue;
        const R inv_u2 = R(1) / u_norm2;

        // Apply P = I - 2 u u^T / ||u||^2 from the left (rows k+1 .. n-1)
        for (std::size_t j = k; j < n; j++) {
            R dot = R(0);
            for (std::size_t i = 0; i < len; i++) {
                dot += u[i] * vcp::tsparse_scalar::real_part(A[k + 1 + i][j]);
            }
            dot *= R(2) * inv_u2;
            for (std::size_t i = 0; i < len; i++) {
                A[k + 1 + i][j] -= T(dot * u[i]);
            }
        }

        // Apply P from the right (cols k+1 .. n-1)
        for (std::size_t i = 0; i < n; i++) {
            R dot = R(0);
            for (std::size_t jj = 0; jj < len; jj++) {
                dot += vcp::tsparse_scalar::real_part(A[i][k + 1 + jj]) * u[jj];
            }
            dot *= R(2) * inv_u2;
            for (std::size_t jj = 0; jj < len; jj++) {
                A[i][k + 1 + jj] -= T(dot * u[jj]);
            }
        }
    }
}

// Reduce a general n x n real matrix A to upper Hessenberg form in-place,
// accumulating the orthogonal transformation Q such that
//   Q^T * A_original * Q = A_hessenberg.
// Q is initialized to identity before the first reflection.
// After the call, A contains the upper Hessenberg form.
// No-op for n <= 2.
template <typename T>
void reduce_to_hessenberg_with_q(
    std::vector<std::vector<T> >& A,
    std::vector<std::vector<T> >& Q)
{
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t n = A.size();

    // Initialize Q = I_n
    Q.assign(n, std::vector<T>(n, T(0)));
    for (std::size_t i = 0; i < n; i++) Q[i][i] = T(1);

    if (n <= 2) return;

    for (std::size_t k = 0; k + 2 <= n; k++) {
        const std::size_t len = n - k - 1;

        R sigma = R(0);
        for (std::size_t i = k + 1; i < n; i++) {
            const R xi = vcp::tsparse_scalar::real_part(A[i][k]);
            sigma += xi * xi;
        }
        sigma = vcp::tsparse_scalar::sqrt_value(sigma);

        const R eps10 = vcp::tsparse_scalar::epsilon<R>() * R(10);
        if (sigma <= eps10) continue;

        const R x0 = vcp::tsparse_scalar::real_part(A[k + 1][k]);
        const R sign_x0 = (x0 >= R(0)) ? R(1) : R(-1);
        const R u0 = x0 + sign_x0 * sigma;

        std::vector<R> u(len, R(0));
        u[0] = u0;
        for (std::size_t i = k + 2; i < n; i++) {
            u[i - k - 1] = vcp::tsparse_scalar::real_part(A[i][k]);
        }

        R u_norm2 = R(0);
        for (std::size_t i = 0; i < len; i++) u_norm2 += u[i] * u[i];
        if (u_norm2 <= eps10 * eps10) continue;
        const R inv_u2 = R(1) / u_norm2;

        // Apply P from the left to A
        for (std::size_t j = k; j < n; j++) {
            R dot = R(0);
            for (std::size_t i = 0; i < len; i++) {
                dot += u[i] * vcp::tsparse_scalar::real_part(A[k + 1 + i][j]);
            }
            dot *= R(2) * inv_u2;
            for (std::size_t i = 0; i < len; i++) {
                A[k + 1 + i][j] -= T(dot * u[i]);
            }
        }

        // Apply P from the right to A
        for (std::size_t i = 0; i < n; i++) {
            R dot = R(0);
            for (std::size_t jj = 0; jj < len; jj++) {
                dot += vcp::tsparse_scalar::real_part(A[i][k + 1 + jj]) * u[jj];
            }
            dot *= R(2) * inv_u2;
            for (std::size_t jj = 0; jj < len; jj++) {
                A[i][k + 1 + jj] -= T(dot * u[jj]);
            }
        }

        // Accumulate Q: Q_new = Q * P (apply P from the right to Q)
        for (std::size_t i = 0; i < n; i++) {
            R dot = R(0);
            for (std::size_t jj = 0; jj < len; jj++) {
                dot += vcp::tsparse_scalar::real_part(Q[i][k + 1 + jj]) * u[jj];
            }
            dot *= R(2) * inv_u2;
            for (std::size_t jj = 0; jj < len; jj++) {
                Q[i][k + 1 + jj] -= T(dot * u[jj]);
            }
        }
    }
}

} // namespace tsparse_proj_hess
} // namespace vcp

#endif // VCP_TSPARSE_PROJECTED_HESSENBERG_HPP
