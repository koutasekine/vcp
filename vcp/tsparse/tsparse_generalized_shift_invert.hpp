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
#include <vcp/tsparse/tsparse_sparse_lu.hpp>

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

// ---------------------------------------------------------------------------
// lu_shift_invert_operator  (E-A1)
//
// Implements the shift-invert operator with a supernodal sparse LU direct
// solver as the inner linear solver ("factorize once, solve repeatedly"):
//
//   standard    problem:  y = (A - sigma*I)^{-1} x
//   generalized problem:  y = (A - sigma*B)^{-1} (B x)
//
// Solve strategy: sparse_lu_factorize_with_info() at construction, then a
// const fac_.solve() per apply; when slu_opts.iterative_refinement is true,
// each apply uses sparse_lu_solve_refined() (solve-time IR) instead.
//
// Design notes (E-A1 design doc):
//  * finalize classification: this operator follows the finalize_policy §2 Q6
//    pattern (same as the preconditioners): all sparse inputs are consumed at
//    CONSTRUCTION time via as_csr()/finalize into owned internal copies; no
//    raw caller matrix is scanned during the iteration.
//  * shifted_ (= A - sigma*I or A - sigma*B) is RETAINED as a member because
//    solve-time iterative refinement needs the shifted matrix itself for the
//    residual r = b - (A - sigma*B) x.
//  * D-4: when the LU factorization fails, this operator does NOT fall back
//    to any iterative solver.  It records factorization_ok() == false and a
//    sparse_lu_status-based diagnostics string; the CALLER must check
//    factorization_ok() right after construction and report the failure
//    (status = "factorization_failed").  apply() on a failed operator is a
//    contract violation and throws vcp::state_error.
// ---------------------------------------------------------------------------
template <class SparseMatrix>
class lu_shift_invert_operator {
public:
    typedef typename SparseMatrix::value_type   value_type;
    typedef typename SparseMatrix::index_type   index_type;
    typedef typename vcp::tsparse_scalar::real_type<value_type>::type real_type;

    // Standard problem:  y = (A - sigma*I)^{-1} x
    lu_shift_invert_operator(
        const SparseMatrix& A,
        value_type sigma,
        const vcp::sparse_lu_options<value_type>& slu_opts)
        : has_B_(false)
        , slu_opts_(slu_opts)
        , n_(static_cast<std::size_t>(A.rowsize()))
        , sigma_(sigma)
        , linear_solves_(0)
        , inner_iterations_(0)
        , inner_failure_count_(0)
        , inner_residual_norm_(real_type(0))
        , factorization_ok_(false)
    {
        // Same diagonal-loop construction as the legacy standard shift-invert
        // paths (get -> set(cur - sigma) -> finalize; C-3 O(1) amortized set).
        SparseMatrix A_csr = A.as_csr();
        SparseMatrix shifted = A_csr;
        for (std::size_t i = 0; i < n_; i++) {
            const value_type cur = shifted.get(static_cast<index_type>(i),
                                               static_cast<index_type>(i));
            shifted.set(static_cast<index_type>(i), static_cast<index_type>(i),
                        cur - sigma_);
        }
        shifted.finalize();
        shifted_ = shifted.as_csr();
        factorize_();
    }

    // Generalized problem:  y = (A - sigma*B)^{-1} (B x)
    lu_shift_invert_operator(
        const SparseMatrix& A,
        const SparseMatrix& B,
        value_type sigma,
        const vcp::sparse_lu_options<value_type>& slu_opts)
        : has_B_(true)
        , slu_opts_(slu_opts)
        , n_(static_cast<std::size_t>(A.rowsize()))
        , sigma_(sigma)
        , linear_solves_(0)
        , inner_iterations_(0)
        , inner_failure_count_(0)
        , inner_residual_norm_(real_type(0))
        , factorization_ok_(false)
    {
        shifted_ = subtract_scaled_sparse(A, B, sigma_);
        B_csr_ = B.as_csr();
        factorize_();
    }

    std::size_t rows() const { return n_; }
    std::size_t cols() const { return n_; }

    // apply: y = (A - sigma*B)^{-1} (B x)   (standard: B = I, i.e. rhs = x)
    //
    // Precondition: factorization_ok() == true (the caller must have checked
    // right after construction; violating this is a caller bug).
    void apply(const std::vector<value_type>& x, std::vector<value_type>& y)
    {
        if (!factorization_ok_) {
            vcp::throw_error<vcp::state_error>(
                "tsparse::lu_shift_invert_operator::apply: "
                "factorization failed; cannot apply operator");
        }

        std::vector<value_type> t;
        if (has_B_) t = B_csr_.mul_vec(x);
        const std::vector<value_type>& rhs = has_B_ ? t : x;

        if (slu_opts_.iterative_refinement) {
            vcp::sparse_lu_refinement_info<value_type> ir;
            y = vcp::sparse_lu_solve_refined(shifted_, fac_, rhs, slu_opts_, &ir);
            inner_iterations_ += ir.iterations;
            if (!ir.converged) inner_failure_count_++;
            inner_residual_norm_ = ir.final_residual;
        } else {
            y = fac_.solve(rhs);
        }
        linear_solves_++;
    }

    void operator()(const std::vector<value_type>& x, std::vector<value_type>& y)
    { apply(x, y); }

    bool factorization_ok() const { return factorization_ok_; }

    // sparse_lu_status string + factorization summary (E-A1 E4).
    std::string factorization_diagnostics() const { return factorization_diagnostics_; }

    // G-1.2 [REVISED BY SLU-CLN1 C1, 2026-07-05]: previously sourced from
    // sparse_lu_info::within_panel_zero_pivot_count, a transitional §17.2(B)
    // prototype diagnostic REMOVED together with the prototype pass.  On the
    // production paths (supernodal true-numeric AND baseline GP alike) zero
    // pivots abort the factorization and are reported through info().status,
    // so a usable factorization always has zero such events; return 0
    // (no fabricated value), preserving the E-A1 E4 accessor contract.
    std::size_t factorization_zero_pivots() const {
        return 0u;
    }

    // Diagnostic accumulators (E-A1 E4 semantics for the sparse_lu path):
    //   linear_solves()       -- number of apply() solves performed
    //   inner_iterations()    -- TOTAL solve-time IR iterations (0 if IR off)
    //   inner_failure_count() -- number of solves whose IR did not reach the
    //                            IR tolerance (0 if IR off)
    //   inner_residual_norm() -- IR final_residual of the LAST solve (0 if IR off)
    std::size_t linear_solves()       const { return linear_solves_; }
    std::size_t inner_iterations()    const { return inner_iterations_; }
    std::size_t inner_failure_count() const { return inner_failure_count_; }
    real_type   inner_residual_norm() const { return inner_residual_norm_; }

private:
    // Factorize shifted_ and record status.  Never throws on numerical
    // failure: sparse_lu_factorize_with_info reports through info().status
    // (SLU-GT1 P3 non-throw contract) and the failure is exposed to the
    // caller via factorization_ok() / factorization_diagnostics() (D-4).
    void factorize_()
    {
        fac_ = vcp::sparse_lu_factorize_with_info(shifted_, slu_opts_);
        factorization_ok_ = fac_.info().success;
        std::ostringstream os;
        os << "sparse_lu status="
           << vcp::sparse_lu_status_to_string(fac_.info().status)
           << "; n=" << fac_.info().n
           << "; nnz_L=" << fac_.info().nnz_L
           << "; nnz_U=" << fac_.info().nnz_U
           << "; supernodes=" << fac_.info().number_of_supernodes
           << "; iterative_refinement="
           << (slu_opts_.iterative_refinement ? "on" : "off");
        factorization_diagnostics_ = os.str();
    }

    SparseMatrix  shifted_;   // A - sigma*I / A - sigma*B (kept for IR residual)
    SparseMatrix  B_csr_;     // generalized only (empty for standard)
    bool          has_B_;
    vcp::sparse_lu_factorization<value_type, index_type> fac_;
    vcp::sparse_lu_options<value_type> slu_opts_;
    std::size_t   n_;
    value_type    sigma_;

    std::size_t   linear_solves_;
    std::size_t   inner_iterations_;
    std::size_t   inner_failure_count_;
    real_type     inner_residual_norm_;

    bool          factorization_ok_;
    std::string   factorization_diagnostics_;
};

} // namespace tsparse
} // namespace vcp

#endif // VCP_TSPARSE_GENERALIZED_SHIFT_INVERT_HPP
