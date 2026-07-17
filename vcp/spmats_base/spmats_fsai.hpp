// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_fsai.hpp
// Policy-layer types for the static FSAI factored sparse approximate inverse
// (FSAI-1): sqrt-free factor triple R ~ A^{-1} = P U D^{-1} U^T P^T after
// [JFSG15] (C. Janna, M. Ferronato, F. Sartoretto, G. Gambolati, "FSAIPACK:
// A Software Package for High-Performance Factored Sparse Approximate
// Inverse Preconditioning", ACM TOMS 41(2), Art. 10, 2015).  Types:
// fsai_status (+to_string) / fsai_options<T> / fsai_result<T,Index> (design
// fsai_design_v0.md SS2.1).
// This file is included by spmats.hpp BEFORE the spmats class body (the
// types appear in policy method signatures).  The out-of-line policy method
// definitions live in spmats_base/spmats_fsai_impl.hpp (included after the
// class body, like spmats_ainv_impl.hpp).
//
// This header also hosts the sparse-layer OpenMP guard macro (design SS6.1;
// first OpenMP use in the sparse policy layer, bfem 7c generation mirrored
// from vcp/bfem/fe_space.hpp L45-56 with the VCP_SPMATS_* names) so it is
// available to any spmats_base implementation header.

#ifndef VCP_SPMATS_FSAI_HPP
#define VCP_SPMATS_FSAI_HPP

// --- sparse-layer OpenMP guard (design SS6.1, bfem 7c generation) ---------
// VCP_NOMP        : global kill switch (implies VCP_SPMATS_NOMP)
// VCP_SPMATS_NOMP : sparse-layer kill switch
// VCP_SPMATS_USE_OPENMP : 0/1, the ONLY symbol implementation code tests
#ifdef VCP_NOMP
#  ifndef VCP_SPMATS_NOMP
#    define VCP_SPMATS_NOMP
#  endif
#endif
#if defined(_OPENMP) && !defined(VCP_SPMATS_NOMP)
#  define VCP_SPMATS_USE_OPENMP 1
#  include <omp.h>
#else
#  define VCP_SPMATS_USE_OPENMP 0
#endif

#include <cstddef>

#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_sparse_chol.hpp>   // sparse_chol_ordering (reused enum, F-D8)

namespace vcp {

// fsai_status: non-throwing status of the FSAI construction / apply /
// residual-norm-estimate entry points (design SS1.2: numerical events --
// tiny pivots, non-SPD input, singular input -- are NEVER a failure status;
// pivot lifting is applied silently and reported through the result
// diagnostics, F-D3/F-D4).  Independent enum, NOT shared with ainv_status
// (module independence, design OPEN-F10; the three values are isomorphic).
enum class fsai_status {
    success,        // outputs are valid (pivot modifications included)
    invalid_input,  // non-square input / dimension mismatch / bad perm /
                    // non-diagonal or structurally-zero-diagonal D
    internal_error
};

inline const char* fsai_status_to_string(fsai_status s) {
    switch (s) {
    case fsai_status::success:        return "success";
    case fsai_status::invalid_input:  return "invalid_input";
    case fsai_status::internal_error: return "internal_error";
    }
    return "unknown";
}

// fsai_options<T>: static-FSAI controls (design SS2.1).  All tolerances are
// REAL-typed and T-dependent constants come from decimal_power_negative (no
// double literals).  The ordering enum is the chol/ldl one reused as-is
// (F-D8: same four choices natural / rcm / amd / nested_dissection; the
// perm convention is also chol's: new->old, P(p[k],k) = 1).
template <class T>
struct fsai_options {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    // --- pattern generation ([JFSG15] Algorithm 2) ---
    std::size_t pattern_power;        // k: S = Low((A')^k) pattern union diagonal
    real_type   prefilter_tolerance;  // tau_pre: drop |a_ij| < tau_pre*sqrt(|a_ii||a_jj|)
                                      // when building A' (implemented as the
                                      // squared comparison |a_ij|^2 <
                                      // tau_pre^2*|a_ii|*|a_jj| -- sqrt-free
                                      // equivalent form, design OPEN-F2).
                                      // 0 = no prefilter.
    real_type   max_density;          // mu_max: recursion stops once
                                      // nnz(B_i) >= mu_max*nnz(A).  0 = unlimited.

    // --- ordering (F-D8; FSAI is built for P^T A P) ---
    sparse_chol_ordering ordering;    // auto_select resolves to natural here
                                      // (conservative FSAI default; NOT chol's
                                      // auto->amd, design OPEN-F3 / OPEN-F13)

    // --- post-filtration ([JFSG15] Algorithm 5) ---
    real_type   postfilter_tolerance;      // tau_post: drop |g_ij| < tau_post*||g_i||_2
                                           // (off-diagonal only).  0 = off.
    std::size_t postfilter_max_nnz_per_row; // mmax, unit diagonal not counted.  0 = unlimited.

    // --- row-solve safeguarding (F-D3/F-D4, AINV SS4.3 thinking) ---
    real_type   pivot_small_threshold;   // |p| <= this -> sign-preserving lift
    real_type   pivot_lift_value;        // lifted pivot magnitude

    fsai_options()
        : pattern_power(1),                                                   // provisional (design v0 OPEN-F1)
          prefilter_tolerance(real_type(0)),                                  // provisional (design v0 OPEN-F1)
          max_density(real_type(5)),                                          // provisional (design v0 OPEN-F1)
          ordering(sparse_chol_ordering::natural),                            // provisional (design v0 OPEN-F3)
          postfilter_tolerance(real_type(0)),                                 // provisional (design v0 OPEN-F4)
          postfilter_max_nnz_per_row(0),                                      // provisional (design v0 OPEN-F4)
          pivot_small_threshold(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(14)),    // provisional (design v0 OPEN-F5)
          pivot_lift_value(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(3)) {}   // provisional (design v0 OPEN-F5)
};

// fsai_result<T,Index>: status + integer diagnostics only (U / D / perm are
// out parameters of policy_fsai_with_info).  All diagnostics are advisory
// (F-D4).  R = P U D^{-1} U^T P^T is never materialized.
template <class T, class Index>
struct fsai_result {
    fsai_status status;

    Index nnz_U;                    // factor fill (post-filtration applied)
    Index nnz_pattern;              // nnz of the pattern S (before filtration)
    Index n_pivot_modifications;    // row-solve lift count, advisory
    Index first_modified_row;       // first (lowest, permuted) row with a lift; -1 = none
    Index n_postfilter_dropped;     // components rejected by post-filtration
    Index n_postfilter_restore_skipped;  // rows whose D-scale restoration was
                                         // skipped because 1 + eps^T A eps was
                                         // not certified positive (design SS3.5;
                                         // field added over the SS2.1 draft --
                                         // OPEN-F13 in the completion report)

    fsai_result()
        : status(fsai_status::internal_error),
          nnz_U(Index(0)), nnz_pattern(Index(0)),
          n_pivot_modifications(Index(0)),
          first_modified_row(Index(-1)),
          n_postfilter_dropped(Index(0)),
          n_postfilter_restore_skipped(Index(0)) {}
};

} // namespace vcp

#endif // VCP_SPMATS_FSAI_HPP
