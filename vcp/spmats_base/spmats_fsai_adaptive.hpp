// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_fsai_adaptive.hpp
// Policy-layer types for the ADAPTIVE FSAI factored sparse approximate
// inverse (FSAI-2): [JFSG15] Algorithm 3 (adaptive pattern generation,
// C. Janna, M. Ferronato, F. Sartoretto, G. Gambolati, "FSAIPACK: ...",
// ACM TOMS 41(2), Art. 10, 2015) producing the SAME sqrt-free factor triple
// as static FSAI-1: U (unit upper triangular), Dhat (diagonal), perm
// (new->old), with R := P U Dhat^{-1} U^T P^T ~ A^{-1} never materialized
// (design fsai_adaptive_design_v0.md SS2.1, F2-D6).  Types:
// fsai_adaptive_options<T> / fsai_adaptive_result<T,Index>.
// fsai_status and the sparse-layer OpenMP guard are REUSED from
// spmats_base/spmats_fsai.hpp (no new enum, design SS2.1); this header is
// included by spmats.hpp BEFORE the spmats class body.  The out-of-line
// policy method definitions live in spmats_base/spmats_fsai_adaptive_impl.hpp
// (included after the class body).  The FSAI-1 files are NOT modified by
// this track (directive invariant 10).

#ifndef VCP_SPMATS_FSAI_ADAPTIVE_HPP
#define VCP_SPMATS_FSAI_ADAPTIVE_HPP

#include <cstddef>

#include <vcp/spmats_base/spmats_fsai.hpp>   // fsai_status + OpenMP guard + sparse_chol_ordering (reused)

namespace vcp {

// fsai_adaptive_options<T>: adaptive-FSAI controls (design SS2.1).  All
// tolerances are REAL-typed; T-dependent constants come from
// decimal_power_negative (no double literals).  The ordering enum is the
// chol/ldl one reused as-is; the initial-value overload IGNORES it and
// reuses perm0 (design SS3.6).
template <class T>
struct fsai_adaptive_options {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    // --- adaptive iteration ([JFSG15] Algorithm 3) ---
    std::size_t max_iterations;       // k_iter: pattern-adaptation rounds per row
    std::size_t entries_per_step;     // s: candidates added per round
    real_type   exit_tolerance;       // eps: row is cut off once
                                      // |psi_k| <= eps * |psi_0| ([JFSG15]
                                      // eq. (25); real_type abs comparison,
                                      // design SS3.3)
    real_type   inner_drop_tolerance; // tau of the in-iteration drop
                                      // ([JFSG15] eq. (26)) applied in the
                                      // SQUARED comparison form
                                      // |g_j|^2 <= tau^2 * ||g||_2^2
                                      // (sqrt-free equivalent, directive
                                      // invariant 11).  0 = no drop.

    // --- ordering / post-filtration / safeguarding (same meaning as FSAI-1) ---
    sparse_chol_ordering ordering;             // diagonal-init overload only;
                                               // auto_select resolves to
                                               // natural (FSAI-1 OPEN-F13
                                               // convention kept)
    real_type   postfilter_tolerance;          // terminal single application
    std::size_t postfilter_max_nnz_per_row;    // (design SS3.5); 0/0 = off
    real_type   pivot_small_threshold;         // dense-solve pivots AND psi
    real_type   pivot_lift_value;              // (design SS3.3 safeguarding)

    fsai_adaptive_options()
        : max_iterations(10),                                                 // provisional (design v0 OPEN-A1)
          entries_per_step(5),                                                // provisional (design v0 OPEN-A1)
          exit_tolerance(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(3)),     // provisional (design v0 OPEN-A1)
          inner_drop_tolerance(real_type(0)),                                 // provisional (design v0 OPEN-A2)
          ordering(sparse_chol_ordering::natural),                            // provisional (FSAI-1 OPEN-F3 inherited)
          postfilter_tolerance(real_type(0)),                                 // provisional (FSAI-1 OPEN-F4 inherited)
          postfilter_max_nnz_per_row(0),                                      // provisional (FSAI-1 OPEN-F4 inherited)
          pivot_small_threshold(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(14)),    // provisional (design v0 OPEN-A4 = OPEN-F5)
          pivot_lift_value(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(3)) {}   // provisional (design v0 OPEN-A4 = OPEN-F5)
};

// fsai_adaptive_result<T,Index>: status + integer diagnostics only (U / D /
// perm are out parameters).  All diagnostics are advisory (F-D4 inherited).
// n_postfilter_restore_skipped is added over the design SS2.1 draft for
// parity with FSAI-1's fsai_result (OPEN-F14 precedent; recorded as a new
// OPEN item in the completion report).
template <class T, class Index>
struct fsai_adaptive_result {
    fsai_status status;

    Index nnz_U;                    // final fill (post-filtration applied)
    Index n_pivot_modifications;    // dense-solve + psi lift count, advisory
    Index first_modified_row;       // first (lowest, permuted) row with a lift; -1 = none
    Index n_postfilter_dropped;     // components rejected by the terminal post-filtration
    Index n_postfilter_restore_skipped;  // rows whose Dhat-scale restoration
                                         // was skipped (sign test failed)
    Index n_rows_converged;         // rows cut off by the eq. (25) test (advisory)
    Index n_rows_capped;            // rows that ran the full k_iter budget (advisory)

    fsai_adaptive_result()
        : status(fsai_status::internal_error),
          nnz_U(Index(0)),
          n_pivot_modifications(Index(0)),
          first_modified_row(Index(-1)),
          n_postfilter_dropped(Index(0)),
          n_postfilter_restore_skipped(Index(0)),
          n_rows_converged(Index(0)),
          n_rows_capped(Index(0)) {}
};

} // namespace vcp

#endif // VCP_SPMATS_FSAI_ADAPTIVE_HPP
