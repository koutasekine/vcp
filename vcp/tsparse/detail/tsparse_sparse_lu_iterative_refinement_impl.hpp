// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// ---------------------------------------------------------------------------
// O4.1a: iterative refinement (IR) -- solve-time post-process.
//
// Classical same-type iterative refinement on the ORIGINAL system A x = b,
// reusing an existing factorization:
//
//     x  = solve(b)                 (initial solve)
//     for k = 1..maxit:
//         r = b - A x               (residual on ORIGINAL A, in working type T)
//         d = solve(r)              (reuse the existing factorization)
//         x = x + d
//         stop when ||r|| / max(||b||,1) <= tol, or k == maxit,
//                or ||r|| stops decreasing (non-convergence guard).
//
// Invariants (directive O4.1a):
//   * Opt-in: runs only when opt.iterative_refinement == true.  When false this
//     reduces to a single fac.solve(b) -- byte-identical to the default path.
//   * Standalone: never modifies the factorization; pure post-process.  Works
//     with any pivoting mode (default threshold_partial, static_mc64, ...).
//   * Deterministic + bounded: always stops by maxit; the non-convergence guard
//     stops as soon as the residual fails to decrease (no divergence/loop).
//   * No-harm: tracks the BEST iterate (minimum residual) and returns it, so a
//     refinement step that does not help cannot degrade the solution.
//   * T-agnostic: residual and update are ordinary scalar arithmetic in T; no
//     mixed / higher-precision residual is introduced.
//
// This header has no `namespace vcp` wrapper; it is injected into namespace vcp
// by the include site in tsparse_sparse_lu.hpp (after the factorization class,
// the CSC convert helpers, and the solve impl are all available).
// ---------------------------------------------------------------------------
#ifndef VCP_TSPARSE_SPARSE_LU_ITERATIVE_REFINEMENT_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_ITERATIVE_REFINEMENT_IMPL_HPP

namespace sparse_lu_detail {

// Residual r = b - A*x, with A in CSC and all arithmetic in the working type T.
// A is the ORIGINAL (unpermuted) matrix; solve() already un-applies Dr/Dc/perm,
// so x and b live in the original coordinate system.
template <class T, class Index>
std::vector<T> sparse_lu_residual_b_minus_Ax(
    const csc_storage<T, Index>& A,
    std::size_t                  n,
    const std::vector<T>&        b,
    const std::vector<T>&        x)
{
    std::vector<T> r = b;  // r <- b, then subtract A*x column by column
    for (std::size_t j = 0; j < n; ++j) {
        const T xj = x[j];
        const Index p_end = A.col_ptr[j + 1u];
        for (Index p = A.col_ptr[j]; p < p_end; ++p) {
            r[static_cast<std::size_t>(A.row_ind[static_cast<std::size_t>(p)])]
                -= A.values[static_cast<std::size_t>(p)] * xj;
        }
    }
    return r;
}

}  // namespace sparse_lu_detail

// ---------------------------------------------------------------------------
// sparse_lu_solve_refined
//
// Opt-in iterative refinement solve.  Returns the (possibly refined) solution
// of A x = b and, optionally, IR diagnostics via ir_info.
//
//   * opt.iterative_refinement == false  -> returns fac.solve(b) unchanged
//                                           (byte-identical default path).
//   * opt.iterative_refinement == true   -> runs the IR loop described above,
//                                           returning the best iterate.
// ---------------------------------------------------------------------------
template <class Matrix>
std::vector<typename Matrix::value_type>
sparse_lu_solve_refined(
    const Matrix&                                                       A,
    const sparse_lu_factorization<typename Matrix::value_type,
                                  typename Matrix::index_type>&         fac,
    const std::vector<typename Matrix::value_type>&                     b,
    const sparse_lu_options<typename Matrix::value_type>&               opt,
    sparse_lu_refinement_info<typename Matrix::value_type>*             ir_info = 0)
{
    typedef typename Matrix::value_type                       T;
    typedef typename Matrix::index_type                       Index;
    typedef typename sparse_lu_scalar_policy<T>::real_type    Real;

    // Initial solve.  This is exactly the plain solve; on the default (opt-out)
    // path nothing else happens, so the result is byte-identical to fac.solve(b).
    std::vector<T> x = fac.solve(b);

    sparse_lu_refinement_info<T> info;  // default: performed=false, status=success

    if (!opt.iterative_refinement) {
        if (ir_info) *ir_info = info;
        return x;
    }

    // ---- Opt-in iterative refinement path ----
    info.performed = true;

    const std::size_t n = b.size();

    // Original A in CSC for the residual matvec (converted once, reused).
    csc_storage<T, Index> Acsc = sparse_lu_make_csc_storage(A);
    if (Acsc.col_ptr.size() != n + 1u) {
        vcp::throw_error<vcp::dimension_error>(
            "sparse_lu_solve_refined: matrix/RHS dimension mismatch");
    }

    const Real bnorm = vcp::tsparse_scalar::hermitian_norm_value(b);
    const Real denom = (bnorm > Real(0)) ? bnorm : Real(1);
    const Real tol   = opt.iterative_refinement_tolerance;
    const std::size_t maxit = opt.iterative_refinement_max_iterations;

    std::vector<T> r =
        sparse_lu_detail::sparse_lu_residual_b_minus_Ax(Acsc, n, b, x);
    Real rnorm = vcp::tsparse_scalar::hermitian_norm_value(r);

    info.initial_residual          = rnorm;
    info.initial_relative_residual = rnorm / denom;

    // Best-iterate tracking (no-harm guarantee).
    std::vector<T> x_best     = x;
    Real           best_rnorm = rnorm;
    Real           prev_rnorm = rnorm;

    bool converged = (rnorm <= tol * denom);
    std::size_t k = 0;

    while (k < maxit && !converged) {
        // Correction d solves A d = r using the existing factorization.
        std::vector<T> d = fac.solve(r);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += d[i];
        }
        ++k;

        r = sparse_lu_detail::sparse_lu_residual_b_minus_Ax(Acsc, n, b, x);
        rnorm = vcp::tsparse_scalar::hermitian_norm_value(r);

        if (rnorm < best_rnorm) {
            best_rnorm = rnorm;
            x_best     = x;
        }

        if (rnorm <= tol * denom) {
            converged = true;
            break;
        }
        // Non-convergence guard: stop when the residual stops decreasing.
        if (rnorm >= prev_rnorm) {
            break;
        }
        prev_rnorm = rnorm;
    }

    info.iterations             = k;
    info.converged              = converged;
    info.final_residual         = best_rnorm;
    info.final_relative_residual = best_rnorm / denom;
    info.status = converged ? sparse_lu_status::success
                            : sparse_lu_status::numerical_instability_suspected;

    if (ir_info) *ir_info = info;
    return x_best;
}

#endif  // VCP_TSPARSE_SPARSE_LU_ITERATIVE_REFINEMENT_IMPL_HPP
