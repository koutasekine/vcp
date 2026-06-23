// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_GENERALIZED_SHIFT_INVERT_HPP
#define VCP_TSPARSE_GENERALIZED_SHIFT_INVERT_HPP

#include <cstddef>
#include <sstream>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_factorization.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_spgemm.hpp>

namespace vcp {
namespace tsparse {

// ---------------------------------------------------------------------------
// subtract_scaled_sparse: compute C = A - sigma * B as a sparse matrix.
// Both A and B must be square with the same dimensions.
// Uses the csr_csr_linear_combination kernel.
// ---------------------------------------------------------------------------
template <class SparseMatrix>
SparseMatrix subtract_scaled_sparse(
    const SparseMatrix& A,
    const SparseMatrix& B,
    typename SparseMatrix::value_type sigma)
{
    typedef typename SparseMatrix::value_type value_type;
    typedef typename SparseMatrix::index_type index_type;

    if (A.rowsize() != A.columnsize() || B.rowsize() != B.columnsize()
     || A.rowsize() != B.rowsize()) {
        vcp::throw_error<vcp::dimension_error>(
            "tsparse::subtract_scaled_sparse: dimension mismatch");
    }

    SparseMatrix Ac = A.as_csr();
    SparseMatrix Bc = B.as_csr();
    SparseMatrix C;
    C.resize(Ac.rowsize(), Ac.columnsize());
    C.reserve(Ac.nnz() + Bc.nnz());

    // C = 1*A + (-sigma)*B  via  csr_csr_linear_combination
    vcp::tsparse_spgemm::csr_csr_linear_combination(
        Ac.rowsize(), Ac.columnsize(),
        Ac.outer_index(), Ac.inner_index(), Ac.values(),
        Bc.outer_index(), Bc.inner_index(), Bc.values(),
        value_type(1), -sigma,
        [&](const index_type i, const index_type j, const value_type& v) {
            C.add(i, j, v);
        });
    C.finalize();
    return C;
}

// ---------------------------------------------------------------------------
// generalized_shift_invert_operator
//
// Implements the operator  (A - sigma*B)^{-1} B
// for the Arnoldi shift-invert iteration on  A x = lambda B x.
//
// Solve strategy: ILU(0) preconditioned GMRES on (A - sigma*B).
// B x is a plain SpMV; no dense matrix is formed.
// ---------------------------------------------------------------------------
template <class SparseMatrix>
class generalized_shift_invert_operator {
public:
    typedef typename SparseMatrix::value_type   value_type;
    typedef typename SparseMatrix::index_type   index_type;
    typedef typename vcp::tsparse_scalar::real_type<value_type>::type real_type;

    // Construct and immediately factorize (A - sigma*B).
    // inner_opts: controls inner GMRES (max_iter, tol, restart, relative).
    generalized_shift_invert_operator(
        const SparseMatrix& A,
        const SparseMatrix& B,
        value_type sigma,
        const std::size_t inner_max_iter,
        const real_type   inner_tol,
        const std::size_t inner_restart)
        : B_csr_(B.as_csr())
        , n_(static_cast<std::size_t>(A.rowsize()))
        , sigma_(sigma)
        , inner_max_iter_(inner_max_iter)
        , inner_tol_(inner_tol)
        , inner_restart_(inner_restart)
        , linear_solves_(0)
        , inner_iterations_(0)
        , inner_failure_count_(0)
        , inner_residual_norm_(real_type(0))
        , factorization_ok_(false)
        , factorization_zero_pivots_(0)
    {
        // Build (A - sigma*B) as sparse matrix
        SparseMatrix shifted = subtract_scaled_sparse(A, B, sigma);
        SparseMatrix shiftedCSR = shifted.as_csr();

        // ILU(0) factorization
        ilu_ = vcp::tsparse_factorization::ilu0_factorize<value_type, index_type>(
            shiftedCSR.outer_index(),
            shiftedCSR.inner_index(),
            shiftedCSR.values(),
            n_,
            vcp::tsparse_scalar::decimal_power_negative<real_type>(14));

        factorization_zero_pivots_ = ilu_.zero_pivots;
        factorization_diagnostics_ = ilu_.diagnostics;
        factorization_ok_ = !ilu_.singular_or_unstable;

        // Store (A - sigma*B) CSR for the inner GMRES operator
        shifted_csr_ = shiftedCSR;
    }

    std::size_t rows() const { return n_; }
    std::size_t cols() const { return n_; }

    // apply: y = (A - sigma*B)^{-1} (B x)
    void apply(const std::vector<value_type>& x, std::vector<value_type>& y)
    {
        if (!factorization_ok_) {
            vcp::throw_error<vcp::state_error>(
                "tsparse::generalized_shift_invert_operator::apply: "
                "factorization failed; cannot apply operator");
        }

        // rhs = B * x  (sparse SpMV)
        std::vector<value_type> rhs = B_csr_.mul_vec(x);

        // Solve (A - sigma*B) y = rhs  via ILU(0)-preconditioned GMRES
        struct AV {
            const SparseMatrix* mat;
            void operator()(const std::vector<value_type>& u,
                            std::vector<value_type>& v) const
            { v = mat->mul_vec(u); }
        } av = { &shifted_csr_ };

        struct PV {
            const vcp::tsparse_factorization::ilu0_data<value_type, index_type>* p;
            void operator()(const std::vector<value_type>& r,
                            std::vector<value_type>& z) const
            { z = vcp::tsparse_factorization::ilu0_solve(*p, r); }
        } pv = { &ilu_ };

        const std::size_t actual_restart =
            (inner_restart_ == 0) ? std::min(n_, std::size_t(30)) : inner_restart_;
        const std::size_t actual_max = (inner_max_iter_ == 0) ? (actual_restart + 1) : inner_max_iter_;

        typedef vcp::tsparse_factorization::gmres_result<value_type, AV, PV> gr_t;
        gr_t gr = vcp::tsparse_factorization::gmres_solve<value_type, AV, PV>(
            av, pv, rhs, actual_max, inner_tol_, actual_restart);

        linear_solves_++;
        inner_iterations_ += gr.iterations;
        if (gr.residual_norm > inner_residual_norm_) inner_residual_norm_ = gr.residual_norm;
        if (!gr.converged) {
            inner_failure_count_++;
        }

        y = gr.x;
    }

    void operator()(const std::vector<value_type>& x, std::vector<value_type>& y)
    { apply(x, y); }

    std::size_t linear_solves()       const { return linear_solves_; }
    std::size_t inner_iterations()    const { return inner_iterations_; }
    std::size_t inner_failure_count() const { return inner_failure_count_; }
    real_type   inner_residual_norm() const { return inner_residual_norm_; }

    bool        factorization_ok()         const { return factorization_ok_; }
    std::size_t factorization_zero_pivots() const { return factorization_zero_pivots_; }

    std::string diagnostics() const {
        std::ostringstream os;
        os << "linear_solves=" << linear_solves_
           << "; inner_iterations=" << inner_iterations_
           << "; inner_failure_count=" << inner_failure_count_
           << "; inner_residual_norm=" << inner_residual_norm_
           << "; factorization_ok=" << (factorization_ok_ ? "true" : "false")
           << "; factorization_zero_pivots=" << factorization_zero_pivots_;
        return os.str();
    }

    std::string factorization_diagnostics() const {
        return factorization_diagnostics_;
    }

private:
    SparseMatrix  B_csr_;
    SparseMatrix  shifted_csr_;
    std::size_t   n_;
    value_type    sigma_;
    std::size_t   inner_max_iter_;
    real_type     inner_tol_;
    std::size_t   inner_restart_;

    vcp::tsparse_factorization::ilu0_data<value_type, index_type> ilu_;

    // Mutable counters (apply() increments these)
    mutable std::size_t linear_solves_;
    mutable std::size_t inner_iterations_;
    mutable std::size_t inner_failure_count_;
    mutable real_type   inner_residual_norm_;

    bool        factorization_ok_;
    std::size_t factorization_zero_pivots_;
    std::string factorization_diagnostics_;
};

} // namespace tsparse
} // namespace vcp

#endif // VCP_TSPARSE_GENERALIZED_SHIFT_INVERT_HPP
