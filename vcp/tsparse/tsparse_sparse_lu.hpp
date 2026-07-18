// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_SPARSE_LU_HPP
#define VCP_TSPARSE_SPARSE_LU_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <new>
#include <set>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_convert.hpp>
#include <vcp/tsparse/tsparse_format.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

// SLU-8R.1: tblas/tlapack required for sparse_lu_dense_kernel<T> implementations.
// These are pure template headers; they define things in namespace vcp.
// They MUST be included before namespace vcp opens to avoid double-namespace nesting.
#include <vcp/tblas/tblas.hpp>
#include <vcp/tlapack/tlapack.hpp>
// SLU-K1: blocked GEMM kernels for the opt-in supernode_panel path (options
// field panel_gemm_kernel below).  Pure template header in namespace vcp;
// same include discipline as tblas/tlapack above (before namespace vcp opens).
// Non-specialized scalar types forward verbatim to vcp::tgemm inside it.
#include <vcp/tblas/tblas_blocked.hpp>

namespace vcp {

// ===========================================================================
// Enumerations
// ===========================================================================

enum class sparse_lu_method {
    auto_select,
    baseline_gp,
    supernodal,
    // SLU-SP1 (appended LAST; existing enumerator values unchanged): opt-in
    // left-looking supernode-panel path with zero-padding GEMM updates
    // (design SLU-SP1 §2, Li 2005 §2.3 technique ported to the sequential
    // left-looking factorization).  NEVER selected by auto_select; the
    // default behavior of every existing path is byte-identical (D-3).
    supernode_panel
};

enum class sparse_lu_ordering {
    auto_select,
    natural,
    rcm,
    amd,
    colamd,
    // SLU-MF4: nested dissection (recursive graph bisection, self-contained).
    // Produces large separators -> wide dense fronts -> larger average front
    // width k for the multifrontal numeric source.  Pattern-only column
    // permutation Q (col_perm[new]=old), like rcm/amd/colamd.  Opt-in; the
    // existing orderings are unchanged.
    nested_dissection
};

enum class sparse_lu_pivoting {
    none,
    diagonal,
    threshold_partial,
    static_mc64
};

enum class sparse_lu_status {
    success,
    structural_singularity,
    numerical_singularity,
    zero_pivot,
    near_zero_pivot,
    inconclusive_pivot_test,
    pivot_rejected,
    numerical_instability_suspected,
    unsupported_scalar_type,
    unsupported_pivoting_mode,
    not_implemented,
    memory_allocation_failed,
    invalid_input,       // user/input-domain error (e.g. non-square matrix)
    internal_error       // implementation bug or unexpected state
};

enum class sparse_lu_storage_kind {
    baseline_csc,
    supernodal
};

enum class pivot_decision {
    acceptable,
    reject,
    inconclusive
};

// SLU-8R.4.1: Within-panel factorization event status.
// Reports the worst event seen during §17.2(B) within-panel factorization.
// Priority (highest = most severe): zero_pivot > rejected_pivot > inconclusive_pivot > near_zero_warning > success > not_run.
// Any status != success (and != not_run) means a notable event occurred; status must NOT be
// silently consumed as "success" even when the factor-level solve still uses CSC baseline.
enum class within_panel_factor_status {
    not_run,            // within-panel factorization was not executed (baseline path)
    success,            // all pivots accepted, no abnormal events
    near_zero_warning,  // at least one pivot was below near_zero_tolerance but > zero_tolerance
    inconclusive_pivot, // non-finite value encountered; fallback pivot used
    rejected_pivot,     // diagonal rejected by threshold; max-abs row was used as pivot
    zero_pivot,         // at least one column had all-zero (or effectively zero) pivot
    failed              // other failure (storage invalid, etc.)
};

// SLU-8R.5.1 / B2+ Option A: Storage-native solve dispatch status.
// Records WHY the supernodal solve used native vs fallback path.
// Set at factorization time (before any solve() call) via set_supernodal_solve_info_().
// Exposed in sparse_lu_info::supernodal_solve_status_value.
//
// Priority for fallback (most to least restrictive):
//   invalid_storage > true_numeric_source_false > fallback_csc_used
//
// Semantics:
//   not_attempted:             method != supernodal (baseline/auto path)
//   storage_native_success:    supernodal_solve_native == true; A_eff-origin
//                              true-numeric storage accepted; native §18.2 solve used
//   fallback_csc_used:         fallback for unknown reason (catch-all)
//   true_numeric_source_false: production transitional storage (key R.5 fallback reason)
//   invalid_storage:           storage.valid == false or source_of_truth_storage == false
//   unsupported_shape:         supernodes.empty()
//   zero_pivot:                retained for genuine A_eff-origin pivot failure (reflected
//                              in true_numeric_source==false); NOT set because
//                              transitional within_panel_status was zero_pivot after
//                              true_numeric_source has been accepted.
//   failed:                    retained for genuine solve or storage failure; NOT set
//                              because transitional within_panel_status was failed after
//                              true_numeric_source has been accepted.
enum class supernodal_solve_status {
    not_attempted,               // baseline/auto path; supernodal solve not entered
    storage_native_success,      // native §18.2 solve will be used (supernodal_solve_native==true)
    fallback_csc_used,           // fallback to CSC (catch-all reason)
    true_numeric_source_false,   // production transitional storage (key R.5 fallback reason)
    invalid_storage,             // storage not valid or not source_of_truth
    unsupported_shape,           // supernodes.empty()
    zero_pivot,                  // genuine A_eff-origin pivot failure (true_numeric_source==false)
    failed                       // genuine storage or solve failure
};

// SLU-8R.5.5: Supernodal numeric source origin kind.
// Tracks where the numeric values in panel_values / U_segments came from.
//   none:                     not yet assigned (pre-bootstrap or invalid storage)
//   csc_bootstrap_transitional: values copied from baseline CSC L/U (bootstrap only)
//   a_eff_true_numeric:       values computed from A_eff origin (SLU-8R.5.5)
//   injected_true_numeric:    values injected by test factory (oracle)
enum class supernodal_numeric_source_kind {
    none,
    csc_bootstrap_transitional,
    a_eff_true_numeric,
    injected_true_numeric
};

// SLU-8R.5.5: True numeric factorization status.
// Records outcome of the A_eff-origin supernodal factorization attempt.
//   not_attempted:          factorize_supernodal_from_a_eff not called
//   success:                A_eff-origin factorization succeeded + residual passed
//   unsupported_structure:  storage structure unsuitable for A-origin path
//   missing_a_entry:        A_eff value required for panel init but not found in A_csc
//   pivot_failure:          zero pivot during within-panel factorization
//   residual_not_checked:   sparse residual could not run (degenerate n / malformed
//                           CSC); NOT accepted as true-numeric (SLU-8R.8: no skip-accept)
//   residual_failed:        factorization ran but ||L*U - A_eff|| / scale too large
//   used_csc_numeric_source: CSC numeric values were copied (prohibited in true-numeric)
//   failed:                 other failure (storage invalid, etc.)
enum class supernodal_true_numeric_status {
    not_attempted,
    success,
    unsupported_structure,
    missing_a_entry,
    pivot_failure,
    residual_not_checked,
    residual_failed,
    used_csc_numeric_source,
    failed
};

// ===========================================================================
// sparse_lu_status helper
// ===========================================================================

inline const char* sparse_lu_status_to_string(sparse_lu_status s) {
    switch (s) {
    case sparse_lu_status::success:                          return "success";
    case sparse_lu_status::structural_singularity:           return "structural_singularity";
    case sparse_lu_status::numerical_singularity:            return "numerical_singularity";
    case sparse_lu_status::zero_pivot:                       return "zero_pivot";
    case sparse_lu_status::near_zero_pivot:                  return "near_zero_pivot";
    case sparse_lu_status::inconclusive_pivot_test:          return "inconclusive_pivot_test";
    case sparse_lu_status::pivot_rejected:                   return "pivot_rejected";
    case sparse_lu_status::numerical_instability_suspected:  return "numerical_instability_suspected";
    case sparse_lu_status::unsupported_scalar_type:          return "unsupported_scalar_type";
    case sparse_lu_status::unsupported_pivoting_mode:        return "unsupported_pivoting_mode";
    case sparse_lu_status::not_implemented:                  return "not_implemented";
    case sparse_lu_status::memory_allocation_failed:         return "memory_allocation_failed";
    case sparse_lu_status::invalid_input:                    return "invalid_input";
    case sparse_lu_status::internal_error:                   return "internal_error";
    }
    return "unknown";
}

// ===========================================================================
// SLU-K1: GEMM kernel selector for the opt-in method=supernode_panel path
// (design SLU-K1 D-2).  reference = the pre-K1 vcp::tgemm<T> call (byte-
// identical escape hatch, T-1 gate); blocked = vcp::tblas_blocked::gemm<T>
// (double is blocked; non-specialized scalars forward verbatim to the
// reference inside tblas_blocked, so the selector is byte-neutral for them).
// Read ONLY by the supernode_panel numeric; every other path ignores it.
// ===========================================================================

enum class sparse_lu_panel_gemm_kernel {
    reference,
    blocked
};

// ===========================================================================
// sparse_lu_options<T>
// ===========================================================================

template <class T>
struct sparse_lu_options {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    sparse_lu_method   method;
    sparse_lu_ordering ordering;
    sparse_lu_pivoting pivoting;

    real_type pivot_threshold;
    real_type absolute_pivot_tolerance;  // absolute threshold for pivot acceptance
    real_type zero_tolerance;
    real_type near_zero_tolerance;

    bool equilibration;
    bool iterative_refinement;

    std::size_t blas_min_block_size;
    std::size_t supernode_relaxation;
    std::size_t panel_size;

    real_type growth_threshold;
    real_type static_pivot_perturbation;

    bool compute_condition_estimate;

    // O4.1a: iterative refinement (IR) controls (opt-in; only read when
    // iterative_refinement == true).  These fields are additive and do NOT
    // change any default behavior: the factorization never reads them, and the
    // default solve path is unaffected.  IR is a solve-time post-process applied
    // via sparse_lu_solve_refined().
    //   iterative_refinement_max_iterations: hard upper bound on refinement
    //     iterations (bounded / deterministic; IR always stops by this count).
    //   iterative_refinement_tolerance: target RELATIVE residual
    //     ||b - A x|| / max(||b||, 1).  Default = working-type epsilon.
    std::size_t iterative_refinement_max_iterations;
    real_type   iterative_refinement_tolerance;

    // SLU-MF: select the numeric source for method=supernodal.
    //   false (default): multifrontal driver (production numeric source).
    //   true:            retained left-looking A_eff-origin driver (SLU-8R.5.5),
    //                    kept for diagnostic before/after comparison only.
    // Additive and non-destructive: read ONLY on the method=supernodal path;
    // baseline_gp / auto_select are unaffected.
    bool supernodal_numeric_diagnostic_leftlooking;

    // SLU-MF2: use the self-symbolic dense-front structure builder for the
    // multifrontal numeric source (method=supernodal) instead of the GP-exact
    // CSC bootstrap structure. The self-symbolic structure contains the dense
    // multifrontal (relaxed-amalgamation / AMD) fill, so wide fronts run native
    // (no CSC fallback). Additive and non-destructive: read ONLY on the
    // method=supernodal multifrontal path; the CSC bootstrap remains the
    // fallback structure when this is false or the builder returns invalid.
    bool supernodal_self_symbolic;

    // SLU-MF3: use the in-place frontal-update multifrontal driver, which
    // accumulates each front's Schur complement directly into its parent front's
    // persistent dense buffer, eliminating the separate contribution-block
    // copy-out (cb_emit) and forwarding-replay, and reusing front buffers instead
    // of allocating a fresh frontal matrix per supernode. Requires a self-symbolic
    // (source_of_truth) structure so routing is fully contained (no forwarding).
    // Additive and non-destructive: read ONLY on the method=supernodal
    // multifrontal path; falls back to the standard MF driver / CSC bootstrap when
    // false, when the structure is not self-symbolic, or on any scatter miss.
    bool supernodal_inplace_frontal;

    // SLU-MF7: opt-in of the in-place multifrontal post-factorization residual
    // acceptance gate (||L*U - A_eff|| / scale).  DEFAULT false = the residual
    // check does NOT run (a column-accumulator L*U reconstruction that, on narrow
    // fronts, can cost several times the numeric factorization itself; approximate
    // solvers do not require this self-verification).  Set true to enable the gate:
    // the factorization is accepted only when rel_res <= 1e-6 or abs_res <= 1e-6,
    // and falls back to the GP path when the check fails — byte-identical to the
    // original native path.  Read ONLY on the in-place self-symbolic native path.
    bool supernodal_native_check_residual;

    // SLU-SP1 (appended LAST; additive, default changes nothing): maximum
    // supernode width of the opt-in method=supernode_panel path.  Supernodes
    // may span panel boundaries up to this cap (Li 2005 maxsup).  Read ONLY
    // by the supernode_panel numeric; every other path ignores it.
    std::size_t supernode_panel_maxsup;

    // SLU-K1 (appended LAST): GEMM kernel of the opt-in supernode_panel
    // panel update (design D-2).  Default blocked; reference restores the
    // pre-K1 byte-exact behavior.  Read ONLY by the supernode_panel numeric.
    sparse_lu_panel_gemm_kernel panel_gemm_kernel;

    sparse_lu_options()
        : method(sparse_lu_method::auto_select),
          ordering(sparse_lu_ordering::auto_select),
          pivoting(sparse_lu_pivoting::threshold_partial),
          pivot_threshold(real_type(1e-1)),
          absolute_pivot_tolerance(real_type(0)),
          zero_tolerance(vcp::tsparse_scalar::epsilon<real_type>()),
          near_zero_tolerance(
              vcp::tsparse_scalar::sqrt_value(vcp::tsparse_scalar::epsilon<real_type>())),
          equilibration(false),
          iterative_refinement(true),
          blas_min_block_size(16),
          // SLU-14R.5: default 0 = fundamental-supernode partition (no relaxed
          // amalgamation).  Relaxation is opt-in via supernode_relaxation > 0.
          // The field was previously read nowhere (default 4 was inert); keeping
          // the default at no-merge preserves the SLU-14R.4 rollback baseline
          // exactly (directive SLU-14R.5 §3 L1 / §7).
          supernode_relaxation(0),
          panel_size(8),
          growth_threshold(real_type(1e4)),
          static_pivot_perturbation(
              vcp::tsparse_scalar::sqrt_value(vcp::tsparse_scalar::epsilon<real_type>())),
          compute_condition_estimate(false),
          iterative_refinement_max_iterations(20),
          iterative_refinement_tolerance(
              vcp::tsparse_scalar::epsilon<real_type>()),
          supernodal_numeric_diagnostic_leftlooking(false),
          supernodal_self_symbolic(false),
          supernodal_inplace_frontal(false),
          supernodal_native_check_residual(false),
          supernode_panel_maxsup(64),
          panel_gemm_kernel(sparse_lu_panel_gemm_kernel::blocked) {}
};

// ===========================================================================
// sparse_lu_info<T, Index>
// ===========================================================================

template <class T, class Index>
struct sparse_lu_info {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    bool              success;
    sparse_lu_status  status;
    sparse_lu_method  method_used;   // method actually used for factorization

    Index n;
    Index nnz_L;
    Index nnz_U;
    Index number_of_supernodes;

    real_type growth_factor;
    real_type estimated_condition;    // negative sentinel if not computed
    real_type perturbation_magnitude; // negative sentinel if unused

    std::size_t symbolic_time_ticks;
    std::size_t numeric_time_ticks;
    std::size_t solve_time_ticks;
    std::size_t dense_kernel_time_ticks;

    // Status flags for the supernodal factorization path (method=supernodal).
    //
    // TWO-TIER FLAG SYSTEM — read this note before interpreting individual flags.
    //
    //   Genuine-state flags (track actual operational status):
    //     supernodal_storage_is_numeric_source — true iff the A_eff-origin interleaved
    //       factorization (SLU-8R.5.5) succeeded and supernodal storage is the numeric
    //       source of truth.
    //     supernodal_solve_native              — true iff §18.2 storage-native solve
    //       is in use (requires supernodal_storage_is_numeric_source == true).
    //     supernodal_true_numeric_success      — true iff A_eff-origin factorization
    //       succeeded and the ||L*U - A_eff|| residual check passed.
    //     When all three are true, genuine supernodal numeric source + solve are active.
    //
    //   Historical / frozen conformance markers — ALWAYS false:
    //     true_supernodal_numeric
    //     supernode_panel_update_is_numeric_source
    //     within_panel_factorization_is_numeric_source
    //     A false frozen marker does NOT mean the feature is absent.
    //     These markers tracked prototype-phase sub-goals defined in SLU-8R.0–8R.4.
    //     Genuine numeric source was achieved via the A_eff-origin interleaved driver
    //     (SLU-8R.5.5), which bypasses the prototype conformance gates those markers
    //     were designed to flip.  They are hard-coded false at all assignment sites
    //     and must not be flipped; do not read them as "feature not implemented."
    //
    // uses_supernodal_prototype: true iff method=supernodal was requested and the
    //   factorization entered the supernodal path (SLU-10 + SLU-8R.1+ chain).
    //   The CSC baseline GP factorization runs first as numeric foundation; supernodal
    //   storage is then bootstrapped (SLU-8R.2) and an A_eff-origin interleaved
    //   factorization is attempted (SLU-8R.5.5).  False for baseline/auto paths.
    //
    // true_supernodal_numeric: ALWAYS false.  Historical SLU-8R.0 prototype conformance
    //   marker.  Was designed to flip true only when the two-stage prototype approach
    //   (CSC-bootstrap-first §17.2 pipeline) completed and became sole numeric source —
    //   a design gate superseded by the A_eff-origin interleaved driver (SLU-8R.5.5).
    //   Hard-coded false at all assignment sites; must not be flipped.  Genuine
    //   operational status is reported by supernodal_true_numeric_success and
    //   supernodal_storage_is_numeric_source.
    //
    // dense_kernel_connected: true iff the production path called
    //   sparse_lu_dense_kernel<T> on real data (SLU-8R.1+).
    //   True for method=supernodal when the dense kernel adapter was invoked during
    //   the transitional §17.2(A)/(B) runs or the A_eff-origin factorization.
    //   False for baseline/auto paths.
    //
    // SLU-8R.0.1:
    //   dense_kernel_connected == true means: the production Sparse LU path has
    //   called the dense kernel adapter and real dense kernel work can be measured.
    //   dense_kernel_connected == true does NOT imply:
    //     - true_supernodal_numeric (historical prototype conformance gate; always false)
    //     - dense kernel time is the main factorization component
    //     - A_eff-origin factorization succeeded (check supernodal_true_numeric_success)
    //     - §17.2(A)/(B) supernode-panel update / within-panel factorization done
    //     - §18.2 storage-native supernodal solve in use (check supernodal_solve_native)
    //     - SLU-8 full conformance (as defined by the prototype-era design gate)
    bool uses_supernodal_prototype;
    bool true_supernodal_numeric;
    bool dense_kernel_connected;

    // SLU-8R.2: supernodal storage source-of-truth conversion diagnostics.
    //
    // has_supernodal_storage: true iff the factor owns a supernodal_lu_storage
    //   object (storage_kind == supernodal and supernodal storage is valid).
    //   True for method=supernodal after SLU-8R.2; false for baseline/auto.
    //   Note: storage_kind == supernodal means supernodal storage exists.
    //         It does NOT mean §17.2 true supernodal numeric is complete.
    //
    // supernodal_storage_bootstrapped_from_csc: true iff supernodal storage was
    //   populated from baseline CSC L/U (SLU-8R.2 transition), not by §17.2.
    //
    // supernodal_storage_is_numeric_source: true iff A_eff-origin interleaved
    //   factorization (SLU-8R.5.5) succeeded; equivalent to true_numeric_source
    //   in supernodal_lu_storage.  Set by set_true_numeric_info_.
    //
    // SLU-8R.5 solve diagnostics (factorization-time, not per-call):
    //   supernodal_solve_native: true iff solve() will use storage-native §18.2 path.
    //     Requires: has_supernodal_storage,
    //               accepted source_of_truth supernodal storage,
    //               storage.true_numeric_source == true.
    //     info_.within_panel_status may describe the transitional CSC-bootstrapped
    //     §17.2(B) run and does not gate native solve once A_eff-origin
    //     true_numeric_source storage has been accepted (B2+ Option A).
    //     For production transitional storage (true_numeric_source == false):
    //       supernodal_solve_native == false (CSC fallback used).
    //     For accepted A_eff-origin storage (true_numeric_source == true):
    //       supernodal_solve_native == true.
    //   supernodal_solve_attempted: true iff solve() enters the supernodal dispatch path.
    //   supernodal_solve_fallback_csc: true iff solve() falls back to CSC baseline.
    //   supernodal_solve_used_csc_factor_data: true iff CSC L/U is the actual solve source.
    //   supernodal_solve_used_injected_storage: true for test factory factors.
    //
    // SLU-8R.5.1: supernodal_solve_status_value records WHY native/fallback was chosen.
    //   not_attempted:             baseline/auto method
    //   storage_native_success:    native solve will be used
    //   true_numeric_source_false: key R.5 reason — production transitional storage fallback
    //   invalid_storage:           storage validation failed
    //   zero_pivot / failed:       fatal within-panel status
    //   fallback_csc_used:         fallback, other reason
    bool has_supernodal_storage;
    bool supernodal_storage_bootstrapped_from_csc;
    bool supernodal_storage_is_numeric_source;
    bool supernodal_solve_native;
    bool supernodal_solve_attempted;
    bool supernodal_solve_fallback_csc;
    bool supernodal_solve_used_csc_factor_data;
    bool supernodal_solve_used_injected_storage;
    supernodal_solve_status supernodal_solve_status_value; // SLU-8R.5.1: fallback reason enum

    std::size_t supernodal_panel_value_count;  // panel_values.size() for method=supernodal
    std::size_t supernodal_u_segment_count;    // U_segments.row_ind.size()
    std::size_t supernodal_storage_bytes;      // approximate byte footprint of supernodal storage

    // [SLU-CLN1 C1, 2026-07-05] The transitional §17.2(A)/(B) prototype
    // diagnostics (supernode_panel_update_* / within_panel_* fields,
    // symmetric_pruning_* hooks) were REMOVED together with the prototype
    // pass.  Production §17.2 diagnostics are the supernodal_true_numeric_*
    // counters below.

    // SLU-8R.5.5: True numeric source switch diagnostics.
    //
    // supernodal_true_numeric_attempted: true iff factorize_supernodal_from_a_eff was called.
    // supernodal_true_numeric_success:   true iff A_eff-origin factorization succeeded + residual passed.
    // supernodal_values_initialized_from_A: true iff panel_values / U_segments initialized from A_eff.
    // supernodal_values_initialized_from_csc_numeric: true iff values came from CSC L/U (transitional).
    // supernodal_factorization_residual_checked: true iff ||L*U - A_eff|| was computed.
    // supernodal_factorization_residual_passed:  true iff residual is within tolerance.
    // supernodal_numeric_source: origin of current panel_values / U_segments numeric values.
    // supernodal_true_numeric_status_value: detailed outcome of the true-numeric attempt.
    // supernodal_true_numeric_supernodes_processed: supernodes completed in interleaved driver.
    // supernodal_true_numeric_panel_update_count: (k,j) panel update pairs in interleaved driver.
    // supernodal_true_numeric_within_panel_count: within-panel factorizations in interleaved driver.
    // supernodal_true_numeric_trsm_count: trsm calls in interleaved §17.2(A).
    // supernodal_true_numeric_gemm_count: gemm calls in interleaved §17.2(A).
    // supernodal_true_numeric_gemv_count: gemv calls in interleaved §17.2(A).
    // supernodal_true_numeric_ger_count:  ger calls in interleaved §17.2(B).
    // supernodal_factorization_residual_abs: absolute ||L*U - A_eff||_F residual.
    // supernodal_factorization_residual_rel: relative residual / ||A_eff||_F.
    bool        supernodal_true_numeric_attempted;
    bool        supernodal_true_numeric_success;
    bool        supernodal_values_initialized_from_A;
    bool        supernodal_values_initialized_from_csc_numeric;
    bool        supernodal_factorization_residual_checked;
    bool        supernodal_factorization_residual_passed;
    supernodal_numeric_source_kind    supernodal_numeric_source;
    supernodal_true_numeric_status    supernodal_true_numeric_status_value;
    std::size_t supernodal_true_numeric_supernodes_processed;
    std::size_t supernodal_true_numeric_panel_update_count;
    std::size_t supernodal_true_numeric_within_panel_count;
    std::size_t supernodal_true_numeric_trsm_count;
    std::size_t supernodal_true_numeric_gemm_count;
    std::size_t supernodal_true_numeric_gemv_count;
    std::size_t supernodal_true_numeric_ger_count;
    // SLU-GT1 D2: residual diagnostics carry real_type (requirement-set
    // arithmetic).  Valid ONLY when supernodal_factorization_residual_checked
    // is true; when the flag is false these fields are undefined and must not
    // be read.
    real_type   supernodal_factorization_residual_abs;
    real_type   supernodal_factorization_residual_rel;

    // SLU-8R.6.1: Gate 2 repair — true-numeric path dedicated timing (nanoseconds).
    // These fields are EXCLUSIVELY for the A_eff-origin true-numeric factorization path.
    // They are SEPARATE from transitional-phase timing (run on CSC-bootstrapped
    // storage before the A_eff-origin step; diagnostic only):
    //   dense_kernel_time_ticks        = SLU-8R.1 U-diagonal block kernel +
    //                                    transitional §17.2(A)/(B) kernel ticks
    //   supernode_panel_update_ticks   = transitional §17.2(A) kernel ticks
    //   within_panel_update_ticks      = transitional §17.2(B) kernel ticks
    // No forced nonzero: real measurements only.
    std::size_t supernodal_true_numeric_total_ticks;          // factorize_supernodal_from_a_eff wall-clock
    std::size_t supernodal_true_numeric_dense_kernel_ticks;   // §17.2(A)+(B) dense kernel sum in true-numeric
    std::size_t supernodal_true_numeric_panel_update_ticks;   // §17.2(A) dense kernel in true-numeric
    std::size_t supernodal_true_numeric_within_panel_ticks;   // §17.2(B) dense kernel in true-numeric

    // SLU-SN-OPT: numeric-time breakdown (nanoseconds, same clock). Partitions
    // supernodal_true_numeric_total_ticks into factorization phases.  axis-1 is
    // reported against factorization_ticks (= total - residual), since the
    // residual sanity check is a post-factorization verification phase, not part
    // of the LU factorization ("分解時間", design §17.2 line 894).
    std::size_t supernodal_true_numeric_init_scatter_ticks;
    std::size_t supernodal_true_numeric_symbolic_ticks;
    std::size_t supernodal_true_numeric_panel_nonkernel_ticks;
    std::size_t supernodal_true_numeric_within_nonkernel_ticks;
    std::size_t supernodal_true_numeric_residual_ticks;
    std::size_t supernodal_true_numeric_factorization_ticks;

    // SLU-MF3: multifrontal assembly breakdown (subset of panel_nonkernel_ticks).
    std::size_t supernodal_true_numeric_mf_aeff_ticks;
    std::size_t supernodal_true_numeric_mf_extend_add_ticks;
    std::size_t supernodal_true_numeric_mf_cb_emit_ticks;

    // SLU-PERF (design §17.2 line 894 accounting): dense-kernel FLOP location.
    //
    // These are FLOATING-POINT OPERATION COUNTS (not time), accumulated at every
    // dense-kernel adapter call site in the A_eff-origin true-numeric path. They
    // are independent of BLAS implementation quality (MKL vs reference): they
    // measure WHERE the arithmetic work is, not how fast it runs.
    //
    // FLOP conventions (per design §17.2 measurement spec):
    //   gemm(m,n,k) -> 2*m*n*k      gemv(m,n) -> 2*m*n
    //   trsm(m,n)   -> m*m*n        ger(m,n)  -> 2*m*n
    //
    // The §17.2(A) panel update emits trsm + gemm/gemv; the §17.2(B) within-panel
    // factorization emits ger (rank-1, BLAS-2). gemm is the only BLAS-3 kernel,
    // so its share of total dense-kernel FLOP is the primary "work is concentrated
    // in large gemm" discriminator for the line-894 verdict.
    double      supernodal_true_numeric_flop_trsm;
    double      supernodal_true_numeric_flop_gemm;
    double      supernodal_true_numeric_flop_gemv;
    double      supernodal_true_numeric_flop_ger;

    // gemm call-shape accounting (§17.2(A) only). m=n_off, n=w_j, k=n_inter.
    // Sums let the harness report average gemm shape; max reports the largest call;
    // gemm_dim_hist buckets calls by min(m,n,k) so a "thin supernode / tiny gemm"
    // regime is visible directly. Buckets:
    //   [0]=1 [1]=2 [2]=3..4 [3]=5..8 [4]=9..16 [5]=17..32 [6]=33..64 [7]=65+
    std::size_t supernodal_true_numeric_gemm_m_sum;
    std::size_t supernodal_true_numeric_gemm_n_sum;
    std::size_t supernodal_true_numeric_gemm_k_sum;
    std::size_t supernodal_true_numeric_gemm_max_m;
    std::size_t supernodal_true_numeric_gemm_max_n;
    std::size_t supernodal_true_numeric_gemm_max_k;
    std::size_t supernodal_true_numeric_gemm_dim_hist[8];

    sparse_lu_info()
        : success(false),
          status(sparse_lu_status::not_implemented),
          method_used(sparse_lu_method::auto_select),
          n(Index(0)), nnz_L(Index(0)), nnz_U(Index(0)), number_of_supernodes(Index(0)),
          growth_factor(real_type(0)),
          estimated_condition(real_type(-1)),
          perturbation_magnitude(real_type(-1)),
          symbolic_time_ticks(0), numeric_time_ticks(0),
          solve_time_ticks(0), dense_kernel_time_ticks(0),
          uses_supernodal_prototype(false),
          true_supernodal_numeric(false),
          dense_kernel_connected(false),
          has_supernodal_storage(false),
          supernodal_storage_bootstrapped_from_csc(false),
          supernodal_storage_is_numeric_source(false),
          supernodal_solve_native(false),
          supernodal_solve_attempted(false),
          supernodal_solve_fallback_csc(false),
          supernodal_solve_used_csc_factor_data(false),
          supernodal_solve_used_injected_storage(false),
          supernodal_solve_status_value(supernodal_solve_status::not_attempted),
          supernodal_panel_value_count(0),
          supernodal_u_segment_count(0),
          supernodal_storage_bytes(0),
          supernodal_true_numeric_attempted(false),
          supernodal_true_numeric_success(false),
          supernodal_values_initialized_from_A(false),
          supernodal_values_initialized_from_csc_numeric(false),
          supernodal_factorization_residual_checked(false),
          supernodal_factorization_residual_passed(false),
          supernodal_numeric_source(supernodal_numeric_source_kind::none),
          supernodal_true_numeric_status_value(supernodal_true_numeric_status::not_attempted),
          supernodal_true_numeric_supernodes_processed(0),
          supernodal_true_numeric_panel_update_count(0),
          supernodal_true_numeric_within_panel_count(0),
          supernodal_true_numeric_trsm_count(0),
          supernodal_true_numeric_gemm_count(0),
          supernodal_true_numeric_gemv_count(0),
          supernodal_true_numeric_ger_count(0),
          supernodal_factorization_residual_abs(real_type(0)),
          supernodal_factorization_residual_rel(real_type(0)),
          supernodal_true_numeric_total_ticks(0),
          supernodal_true_numeric_dense_kernel_ticks(0),
          supernodal_true_numeric_panel_update_ticks(0),
          supernodal_true_numeric_within_panel_ticks(0),
          supernodal_true_numeric_init_scatter_ticks(0),
          supernodal_true_numeric_symbolic_ticks(0),
          supernodal_true_numeric_panel_nonkernel_ticks(0),
          supernodal_true_numeric_within_nonkernel_ticks(0),
          supernodal_true_numeric_residual_ticks(0),
          supernodal_true_numeric_factorization_ticks(0),
          supernodal_true_numeric_mf_aeff_ticks(0),
          supernodal_true_numeric_mf_extend_add_ticks(0),
          supernodal_true_numeric_mf_cb_emit_ticks(0),
          supernodal_true_numeric_flop_trsm(0.0),
          supernodal_true_numeric_flop_gemm(0.0),
          supernodal_true_numeric_flop_gemv(0.0),
          supernodal_true_numeric_flop_ger(0.0),
          supernodal_true_numeric_gemm_m_sum(0),
          supernodal_true_numeric_gemm_n_sum(0),
          supernodal_true_numeric_gemm_k_sum(0),
          supernodal_true_numeric_gemm_max_m(0),
          supernodal_true_numeric_gemm_max_n(0),
          supernodal_true_numeric_gemm_max_k(0)
    {
        for (int i = 0; i < 8; ++i) supernodal_true_numeric_gemm_dim_hist[i] = 0;
    }
};

// ===========================================================================
// O4.1a: iterative refinement (IR) diagnostics
//
// Returned (optionally) by sparse_lu_solve_refined().  IR is a SOLVE-TIME
// post-process: it reuses an existing factorization and never modifies it, so
// its diagnostics live here rather than in sparse_lu_info (which is set at
// factorization time).
//
// Honesty contract:
//   - performed  : true iff opt.iterative_refinement was honored and the IR loop
//                  ran (false => the call reduced to the plain solve, byte-
//                  identical to fac.solve(b)).
//   - converged  : true iff the relative residual reached the tolerance within
//                  the iteration budget.  When false, the returned solution is
//                  the BEST iterate seen (minimum ||b - A x||); the caller is
//                  told non-convergence honestly (status != success), never a
//                  silent claim of success.
//   - residuals  : measured in the working type T (T-agnostic; no mixed/higher
//                  precision residual).
// ===========================================================================
template <class T>
struct sparse_lu_refinement_info {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    bool        performed;   // true iff IR actually ran (opt-in honored)
    bool        converged;   // true iff relative residual reached tolerance
    std::size_t iterations;  // refinement iterations actually performed

    real_type   initial_residual;           // ||b - A x0||           (in T)
    real_type   final_residual;             // ||b - A x_best||        (in T)
    real_type   initial_relative_residual;  // initial_residual / max(||b||, 1)
    real_type   final_relative_residual;    // final_residual  / max(||b||, 1)

    // success if converged within budget; numerical_instability_suspected if the
    // tolerance was not reached (best iterate returned, non-convergence reported).
    sparse_lu_status status;

    sparse_lu_refinement_info()
        : performed(false),
          converged(false),
          iterations(0),
          initial_residual(real_type(0)),
          final_residual(real_type(0)),
          initial_relative_residual(real_type(0)),
          final_relative_residual(real_type(0)),
          status(sparse_lu_status::success) {}
};

// ===========================================================================
// Owning storage types
// All Index-typed storage requires a signed Index (sentinel -1 used throughout).
// ===========================================================================

template <class T, class Index>
struct csc_storage {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    std::vector<Index> col_ptr;
    std::vector<Index> row_ind;
    std::vector<T>     values;
};

template <class T, class Index>
struct baseline_lu_storage {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    csc_storage<T, Index> L;   // unit diagonal (implicit)
    csc_storage<T, Index> U;   // explicit diagonal

    std::vector<Index> row_perm;      // row_perm[new_row]     = old_row
    std::vector<Index> inv_row_perm;  // inv_row_perm[old_row] = new_row
    std::vector<Index> col_perm;      // col_perm[new_col]     = old_col
    std::vector<Index> inv_col_perm;  // inv_col_perm[old_col] = new_col

    std::vector<T> Dr;   // row equilibration scaling
    std::vector<T> Dc;   // column equilibration scaling
};

template <class Index>
struct supernode_desc {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");

    // Column range: [first_col, first_col + num_cols)  (exclusive upper bound).
    // Invariant: num_cols > 0.
    Index first_col;
    Index num_cols;  // == last_col - first_col (exclusive convention)

    // row_indices[local_row] = global permuted row id (for gather/scatter).
    // Panel layout: rows 0..num_cols-1 = diagonal block rows,
    //               rows num_cols..row_count-1 = off-diagonal L rows.
    // Invariant: leading_dimension >= row_indices.size().
    std::vector<Index> row_indices;

    // Offset and leading dimension into supernodal_lu_storage::panel_values.
    // panel_values[values_offset + c * leading_dimension + r] for col c, row r.
    Index values_offset;
    Index leading_dimension;

    // Index range into supernodal_lu_storage::U_segments for this supernode.
    // U_segments.seg_ptr[u_segment_start] .. U_segments.seg_ptr[u_segment_start + u_segment_count]
    // gives off-diagonal U entries (rows < first_col) stored column-by-column.
    // Column association within segments is implicit (see supernodal_solve_impl §18.2).
    Index u_segment_start;
    Index u_segment_count;

    // SLU-8R.3.1: per-column entry count within the U_segment for this supernode.
    // u_seg_col_ptr[c] = cumulative entry count for local cols 0..c-1 (relative to segment start).
    // u_seg_col_ptr[num_cols] = u_segment_count.
    // Absolute U_segments index for local col c:
    //   U_segments[u_segment_start + u_seg_col_ptr[c] .. u_segment_start + u_seg_col_ptr[c+1])
    // Size: num_cols + 1.  Empty if num_cols == 0 or no U entries exist.
    std::vector<Index> u_seg_col_ptr;

    supernode_desc()
        : first_col(Index(0)), num_cols(Index(0)),
          values_offset(Index(0)), leading_dimension(Index(0)),
          u_segment_start(Index(0)), u_segment_count(Index(0)) {}
};

template <class T, class Index>
struct u_segment_storage {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    std::vector<Index> seg_ptr;
    std::vector<Index> row_ind;
    std::vector<T>     values;
};

template <class T, class Index>
struct supernodal_lu_storage {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");

    // Supernode descriptors: one per supernode.
    std::vector<supernode_desc<Index> > supernodes;

    // Combined dense panel values (column-major per supernode).
    // Layout per supernode s:
    //   panel_values[values_offset + c * ld + r] = (col c, panel row r) value
    // Upper num_cols x num_cols block: combined L\U (U on/above diag, L below diag).
    // Rows num_cols..row_count-1: off-diagonal L entries.
    std::vector<T> panel_values;

    // Sparse off-diagonal U contributions outside the diagonal block.
    // seg_ptr[s]..seg_ptr[s+1] indexes row_ind/values for supernode s.
    // Entries stored column-by-column (col_begin first) within each segment.
    // Column association is implicit (see supernodal_solve_impl §18.2).
    u_segment_storage<T, Index> U_segments;

    // Permutations and equilibration scaling (copied from baseline CSC L/U for bootstrap).
    // Convention: row_perm[new_row] = old_row, inv_row_perm[old_row] = new_row.
    //             col_perm[new_col] = old_col, inv_col_perm[old_col] = new_col.
    std::vector<Index> row_perm;
    std::vector<Index> inv_row_perm;
    std::vector<Index> col_perm;
    std::vector<Index> inv_col_perm;

    // Equilibration scaling vectors.
    // Dr[i] = row scaling factor for permuted row i.
    // Dc[j] = column scaling factor for permuted column j.
    // Size n or empty (empty = identity scaling).
    std::vector<T> Dr;
    std::vector<T> Dc;

    // SLU-8R.2 transition status flags.
    //
    // valid: true iff this storage object has been successfully populated.
    //
    // bootstrapped_from_csc: true iff values were derived from baseline CSC L/U
    //   (SLU-8R.2 transition). False once §17.2(A)/(B) fills this directly.
    //
    // source_of_truth_storage: true iff the factor object owns this as its
    //   supernodal storage representation.  Does NOT imply true_numeric_source.
    //
    // true_numeric_source: false until §17.2(A)/(B) supernodal numeric path
    //   writes directly to this storage.  MUST remain false in SLU-8R.2.
    //
    // SLU-8R.7.1: storage/info diagnostic contract cleanup.
    //   The §17.2(A) diagnostic-only flag `panel_update_applied` was REMOVED from
    //   storage.  After the SLU-8R.7 re-bootstrap, the final supernodal storage is
    //   the A_eff-origin true-numeric SOURCE storage, not the prototype diagnostic
    //   storage; a storage-resident §17.2(A) flag would be reset by re-bootstrap and
    //   diverge from the authoritative info-side flags.  §17.2(A)/(B) diagnostics
    //   are now authoritative in sparse_lu_info
    //   (info.supernode_panel_update_applied / info.supernode_panel_update_executed).
    bool valid;
    bool bootstrapped_from_csc;
    bool source_of_truth_storage;
    bool true_numeric_source;

    // SLU-8R.5.5: value origin tracking.
    // values_initialized_from_A: true iff panel_values/U_segments set from A_eff.
    // values_initialized_from_csc_numeric: true iff values came from CSC L/U bootstrap.
    // numeric_source_kind: detailed classification of value origin.
    bool values_initialized_from_A;
    bool values_initialized_from_csc_numeric;
    supernodal_numeric_source_kind numeric_source_kind;

    // SLU-MF6: per-front U-right column sets (Ur_s), exported by the self-symbolic
    // builder (build_supernodal_self_symbolic_storage) so the in-place
    // multifrontal driver does not rebuild them from the U_segments transpose.
    // The builder already computes this set per front; it was previously
    // discarded. Empty unless this storage is self-symbolic (source_of_truth).
    // When size() == supernodes.size(), mf_uright_of[s] is exactly the column set
    // the driver would otherwise recompute (sorted ascending, unique, cols > the
    // supernode's pivot range). Does not affect baseline_gp / solve (unused there).
    std::vector<std::vector<Index> > mf_uright_of;

    supernodal_lu_storage()
        : valid(false),
          bootstrapped_from_csc(false),
          source_of_truth_storage(false),
          true_numeric_source(false),
          values_initialized_from_A(false),
          values_initialized_from_csc_numeric(false),
          numeric_source_kind(supernodal_numeric_source_kind::none) {}
};

// ===========================================================================
// SLU-7: Symbolic supernode metadata descriptor
// Carries supernode partition, column map, supernode parent forest, and
// structural row patterns.  Built by sparse_lu_build_supernode_symbolic_csc.
// ===========================================================================

template <class Index>
struct sparse_lu_supernode_symbolic {
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");

    // supernode_ptr[s], supernode_ptr[s+1]) gives the column range of supernode s.
    std::vector<Index> supernode_ptr;        // size nsup + 1

    std::vector<Index> column_to_supernode;  // size n; column -> supernode index

    // parent[s] == -1 (root) or s < parent[s] < nsup (forward supernode parent).
    std::vector<Index> parent;               // size nsup

    // Structural row pattern per supernode (sorted unique).
    // Supernode s occupies row_ind[row_ptr[s] .. row_ptr[s+1]).
    std::vector<Index> row_ptr;              // size nsup + 1
    std::vector<Index> row_ind;              // sorted unique per supernode

    bool valid;

    sparse_lu_supernode_symbolic() : valid(false) {}
};

// ===========================================================================
// SLU-8: Symbolic supernode reach metadata descriptor
// Carries per-supernode ancestor reach in the supernode etree, child lists,
// and symbolic panel row patterns. Built by sparse_lu_build_supernode_reach_symbolic.
//
// SLU-8 symbolic panel row pattern.
// This metadata is a structural over-approximation used for future
// supernodal symbolic planning. It is not the final numeric row structure
// after threshold partial pivoting.
// ===========================================================================

template <class Index>
struct sparse_lu_supernode_reach_symbolic {
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");

    // Per-supernode ancestor reach over the supernode etree (Approach B).
    // reach_ind[reach_ptr[s]..reach_ptr[s+1]) = ancestor supernodes of s.
    // Entries are sorted unique and strictly > s (no self-reach).
    std::vector<Index> reach_ptr;      // size nsup + 1
    std::vector<Index> reach_ind;      // sorted unique, values in (s, nsup)

    // Children of each supernode in the supernode etree.
    // child_ind[child_ptr[s]..child_ptr[s+1]) = children of s (sorted unique).
    std::vector<Index> child_ptr;      // size nsup + 1
    std::vector<Index> child_ind;      // sorted unique, values in [0, nsup)

    // SLU-8 symbolic panel row pattern.
    // This metadata is a structural over-approximation used for future
    // supernodal symbolic planning. It is not the final numeric row structure
    // after threshold partial pivoting.
    std::vector<Index> panel_row_ptr;  // size nsup + 1
    std::vector<Index> panel_row_ind;  // sorted unique, values in [0, n)

    bool valid;

    sparse_lu_supernode_reach_symbolic() : valid(false) {}
};

// ===========================================================================
// SLU-10: Numeric supernode metadata descriptor
// Carries supernode numeric block descriptors derived from actual CSC L/U factors.
// Built by build_supernode_numeric_from_csc after baseline sparse GP factorization.
//
// IMPORTANT: This struct carries CSC-derived metadata only; it is NOT the genuine
// supernodal numeric source.  All row patterns / values derive from actual CSC L/U.
// Symbolic panel_row_ind is NOT used as the final numeric row structure.
// Genuine supernodal numeric source is supernodal_lu_storage (via SLU-8R.5.5).
// ===========================================================================

template <class T, class Index>
struct sparse_lu_supernode_numeric {
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_supernode_numeric: Index must be signed");

    bool valid;

    // Supernode partition mirrors SLU-7 partition.
    std::vector<Index> supernode_ptr;        // size nsup + 1
    std::vector<Index> column_to_supernode;  // size n

    // Column start indices in baseline CSC L/U for each supernode's first column.
    std::vector<Index> l_col_start;  // size nsup
    std::vector<Index> u_col_start;  // size nsup

    // Actual numeric row patterns per supernode, derived from actual CSC L/U.
    // NOT derived from symbolic panel_row_ind.
    std::vector<Index> l_row_ptr;  // size nsup + 1, monotone
    std::vector<Index> l_row_ind;  // sorted unique rows per supernode from actual L
    std::vector<Index> u_row_ptr;  // size nsup + 1, monotone
    std::vector<Index> u_row_ind;  // sorted unique rows per supernode from actual U

    // Diagonal block values from actual U (column-major layout).
    // Supernode s: diag_block_values[diag_block_ptr[s]..diag_block_ptr[s+1])
    // is the width x width panel (width = supernode_ptr[s+1] - supernode_ptr[s]).
    std::vector<Index> diag_block_ptr;    // size nsup + 1
    std::vector<T>     diag_block_values; // column-major

    sparse_lu_supernode_numeric() : valid(false) {}
};

// ===========================================================================
// SLU-11: Supernode numeric metadata diagnostics
// Returned by sparse_lu_make_supernode_numeric_diagnostics and
// sparse_lu_factorization::supernode_numeric_diagnostics().
// Provides a summary of the supernode numeric metadata state and its
// consistency with the actual CSC L/U factors.
// ===========================================================================

template <class Index>
struct sparse_lu_supernode_numeric_diagnostics {
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_supernode_numeric_diagnostics: Index must be signed");

    bool  valid;    // numeric metadata is structurally valid
    Index n;        // matrix size
    Index nsup;     // number of supernodes

    Index nnz_L_rows_total;     // total entries in l_row_ind
    Index nnz_U_rows_total;     // total entries in u_row_ind
    Index diag_block_count;     // number of diagonal blocks (== nsup)
    Index diag_block_value_count; // total entries in diag_block_values

    // csc_backed: true for metadata built from CSC L/U (SLU-10).
    // Always true for valid supernodal numeric metadata derived from CSC L/U.
    bool csc_backed;

    // actual_csc_consistent: l_row_ind/u_row_ind per supernode match the
    // sorted unique union of actual CSC L/U row indices for those columns.
    bool actual_csc_consistent;

    // diag_blocks_consistent: diag_block_values match the actual U CSC entries
    // for the diagonal block of each supernode (column-major layout).
    bool diag_blocks_consistent;

    // symbolic_panel_rows_used_as_numeric_rows: should ALWAYS be false.
    // panel_row_ind is symbolic-only; numeric rows derive from actual CSC L/U.
    bool symbolic_panel_rows_used_as_numeric_rows;

    sparse_lu_supernode_numeric_diagnostics()
        : valid(false),
          n(Index(0)), nsup(Index(0)),
          nnz_L_rows_total(Index(0)), nnz_U_rows_total(Index(0)),
          diag_block_count(Index(0)), diag_block_value_count(Index(0)),
          csc_backed(false),
          actual_csc_consistent(false),
          diag_blocks_consistent(false),
          symbolic_panel_rows_used_as_numeric_rows(false) {}
};

// ===========================================================================
// sparse_lu_symbolic_result<Index>
// Diagnostic Case A: carries success/status so two-stage API can propagate
// symbolic failure to the numeric stage without additional throw paths.
// ===========================================================================

template <class Index>
struct sparse_lu_symbolic_result {
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    Index n;

    bool             success; // false = symbolic failed or not yet run
    sparse_lu_status status;  // set by sparse_lu_symbolic; not_implemented in SLU-0

    std::vector<Index> col_perm;           // col_perm[new_col]     = old_col
    std::vector<Index> inv_col_perm;       // inv_col_perm[old_col] = new_col

    std::vector<Index> col_etree;          // elimination tree parent, -1 sentinel
    std::vector<Index> relaxed_supernodes; // symbolic candidate supernode boundaries

    sparse_lu_supernode_symbolic<Index>      supernode_info;       // SLU-7 symbolic metadata

    Index estimated_fill_upper_bound;
    bool  structural_singularity;

    sparse_lu_symbolic_result()
        : n(Index(0)),
          success(false),
          status(sparse_lu_status::not_implemented),
          estimated_fill_upper_bound(Index(0)),
          structural_singularity(false) {}
};

// ===========================================================================
// sparse_lu_scalar_policy<T>
// ===========================================================================

template <class T>
struct sparse_lu_scalar_policy {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    static real_type abs_value(const T& x) {
        return vcp::tsparse_scalar::abs_value(x);
    }

    // Exact zero: abs_value(x) == 0. Works for complex (magnitude == 0).
    static bool is_exact_zero(const T& x) {
        return abs_value(x) == real_type(0);
    }

    // Legacy: !(abs_value(x) > tol). Does not throw on negative tol.
    // Certified-only idiom (SLU-GT1 P1): for interval types this is
    // "cannot certify abs_value(x) > tol", which falls to the safe side.
    // For totally ordered types it is identical to abs_value(x) <= tol.
    static bool is_zero(const T& x, const real_type& tol) {
        return !(abs_value(x) > tol);
    }

    // Near-zero: !(abs_value(x) > tol). Throws if tol < 0 (invalid input).
    // Certified-only idiom (SLU-GT1 P1); same semantics as is_zero above.
    static bool is_near_zero(const T& x, const real_type& tol) {
        if (tol < real_type(0)) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu: is_near_zero: tolerance must be non-negative");
        }
        return !(abs_value(x) > tol);
    }

    // Magnitude comparisons: compare abs_value(a) vs abs_value(b).
    // Do not require T to have operator<; always compare via real_type.
    static bool abs_less(const T& a, const T& b) {
        return abs_value(a) < abs_value(b);
    }

    static bool abs_less_equal(const T& a, const T& b) {
        return abs_value(a) <= abs_value(b);
    }

    static bool abs_greater_equal(const T& a, const T& b) {
        return abs_value(a) >= abs_value(b);
    }

    // Compare magnitudes: is abs(a) >= abs(b)?
    // Returns inconclusive when either value is non-finite (NaN etc.)
    static pivot_decision compare_abs(const T& a, const T& b) {
        const real_type ra = abs_value(a);
        const real_type rb = abs_value(b);
        if (!vcp::tsparse_scalar::is_finite(ra) || !vcp::tsparse_scalar::is_finite(rb)) {
            return pivot_decision::inconclusive;
        }
        return (ra >= rb) ? pivot_decision::acceptable : pivot_decision::reject;
    }

    // Threshold partial pivot test: abs(pivot) >= threshold * column_max
    static pivot_decision acceptable_pivot(
        const T&         pivot,
        const real_type& column_max,
        const real_type& threshold)
    {
        const real_type ap = abs_value(pivot);
        if (!vcp::tsparse_scalar::is_finite(ap) ||
            !vcp::tsparse_scalar::is_finite(column_max)) {
            return pivot_decision::inconclusive;
        }
        return (ap >= threshold * column_max) ? pivot_decision::acceptable
                                              : pivot_decision::reject;
    }
};

// ===========================================================================
// SLU-3: Pivot parameter validation and acceptable-pivot helper
// ===========================================================================

// Throws invalid_argument if threshold or absolute_tol are out of valid range.
// threshold must be in [0, 1]; absolute_tol must be >= 0.
template <class Real>
void sparse_lu_validate_pivot_parameters(Real threshold, Real absolute_tol)
{
    if (threshold < Real(0) || threshold > Real(1)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu: pivot_threshold must be in [0, 1]");
    }
    if (absolute_tol < Real(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu: absolute_pivot_tolerance must be non-negative");
    }
}

// Returns true iff:
//   1. abs_value(pivot) > absolute_tol
//   2. abs_value(pivot) >= threshold * abs_value(column_max)
// Throws invalid_argument if threshold not in [0,1] or absolute_tol < 0.
//
// SLU-GT1 P1 note: condition 1 is written as !(abs_pivot > absolute_tol) ->
// reject.  The comparison is a certainly-> relation, so for interval scalars a
// pivot whose magnitude cannot be certified positive (e.g. an interval
// containing 0) is rejected here, BEFORE any division can occur.  This is the
// certified direction; do not rewrite it as abs_pivot <= absolute_tol.
template <class T>
bool sparse_lu_is_acceptable_pivot(
    const T& pivot,
    const T& column_max,
    typename sparse_lu_scalar_policy<T>::real_type threshold,
    typename sparse_lu_scalar_policy<T>::real_type absolute_tol)
{
    typedef typename sparse_lu_scalar_policy<T>::real_type Real;
    sparse_lu_validate_pivot_parameters(threshold, absolute_tol);
    const Real abs_pivot  = sparse_lu_scalar_policy<T>::abs_value(pivot);
    const Real abs_colmax = sparse_lu_scalar_policy<T>::abs_value(column_max);
    if (!(abs_pivot > absolute_tol)) return false;
    return abs_pivot >= threshold * abs_colmax;
}

// ===========================================================================
// sparse_lu_dense_kernel<T>  -- SLU-8R.1 dense kernel adapter
//
// Delegates to vcp::tblas (gemm/trsm/gemv/ger) and vcp::tlapack (getrf/getrs).
// tblas / tlapack are NOT modified.
//
// SLU-8R.1: connected to production factorization path (method=supernodal).
//   - info.dense_kernel_connected is true when production path called the adapter.
//   - info.dense_kernel_time_ticks is the measured time of adapter calls.
//   - info.true_supernodal_numeric is always false (historical conformance marker).
//
// dense_kernel_connected == true does NOT imply:
//   - true_supernodal_numeric (historical prototype conformance gate; always false)
//   - §17.2(A)/(B) supernode-panel update / within-panel factorization done
//   - §18.2 storage-native supernodal solve in use (check supernodal_solve_native)
//   - SLU-8 full conformance (as defined by the prototype-era design gate)
//
// ipiv convention: 0-based (matches tgetrf/tgetrs in tlapack).
// All methods: column-major storage.
//
// Implementations are in:
//   vcp/tsparse/detail/tsparse_sparse_lu_dense_kernel_impl.hpp
// (included below, inside namespace vcp)
// ===========================================================================

template <class T>
struct sparse_lu_dense_kernel {

    // General matrix multiply: C = alpha*A*B + beta*C
    // A: m x k column-major (lda >= m), B: k x n (ldb >= k), C: m x n (ldc >= m).
    static void gemm(
        std::size_t m, std::size_t n, std::size_t k,
        const T& alpha, const T* A, std::size_t lda,
        const T* B, std::size_t ldb,
        const T& beta, T* C, std::size_t ldc);

    // Triangular solve (legacy): unit-lower, left side, no transpose.
    // Solves L * X = B in-place.  L: rows x rows (lda >= rows), B: rows x cols (ldb >= rows).
    static void trsm(
        std::size_t rows, std::size_t cols,
        const T* L, std::size_t lda,
        T* B, std::size_t ldb);

    // Triangular solve (full LAPACK-style, SLU-8R.1+).
    // Solves op(A)*X = alpha*B (side='L') or X*op(A) = alpha*B (side='R') in-place.
    // A: triangular (uplo, diag).  op(A) = A ('N') or A^T ('T','C').
    // m: rows of B, n: cols of B.  A is m x m for side='L', n x n for side='R'.
    // For §17.2(A)/(B) and §18.2: ('L','L','N','U') = unit-lower forward,
    //                              ('L','U','N','N') = non-unit upper backward.
    static void trsm(
        char side, char uplo, char trans, char diag,
        std::size_t m, std::size_t n,
        const T& alpha, const T* A, std::size_t lda,
        T* B, std::size_t ldb);

    // Matrix-vector product: y = alpha*A*x + beta*y
    // A: m x n column-major (lda >= m), x: length n, y: length m.
    static void gemv(
        std::size_t m, std::size_t n,
        const T& alpha, const T* A, std::size_t lda,
        const T* x, const T& beta, T* y);

    // Rank-1 update: A = alpha*x*y^T + A
    // A: m x n column-major (lda >= m), x: length m, y: length n.
    static void ger(
        std::size_t m, std::size_t n,
        const T& alpha, const T* x, const T* y,
        T* A, std::size_t lda);

    // Dense LU factorization (LAPACK getrf equivalent).
    // Returns INFO: 0 = success, >0 = zero pivot at INFO-th step.
    // ipiv: 0-based pivot array (size min(rows,cols)).
    // SLU-8R.1 scope: small dense blocks only (NOT supernodal pivoting wholesale).
    // Active panel-height pivot search (§17.2(B)) is NOT delegated here.
    static int getrf(
        std::size_t rows, std::size_t cols,
        T* data, std::size_t lda, int* ipiv);

    // Dense triangular solve using getrf factorization result.
    // Returns INFO: 0 = success.
    // ipiv: 0-based pivot array from getrf.
    // b: n x nrhs column-major (ldb >= n). Overwritten with solution.
    static int getrs(
        std::size_t n, std::size_t nrhs,
        const T* data, std::size_t lda,
        const int* ipiv, T* b, std::size_t ldb);
};

// SLU-8R.1: include out-of-line implementations (uses tblas/tlapack in namespace vcp).
#include <vcp/tsparse/detail/tsparse_sparse_lu_dense_kernel_impl.hpp>

// ===========================================================================
// Forward declarations of API functions (for friend declarations)
// ===========================================================================

template <class T, class Index>
class sparse_lu_factorization;

// ---------------------------------------------------------------------------
// Forward declarations needed for solve dispatch (defined in detail include).
// Placed before the class so solve() can call them inline.
// ---------------------------------------------------------------------------
namespace sparse_lu_detail {

template <class T, class Index>
std::vector<T>
solve_baseline_storage(
    const baseline_lu_storage<T, Index>& lu,
    Index n,
    const std::vector<T>& b);

// SLU-12: supernode-aware solve path -- forward declaration for solve() dispatch
template <class T, class Index>
std::vector<T>
solve_baseline_storage_supernode_aware(
    const baseline_lu_storage<T, Index>& lu,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernode_info,
    const std::vector<T>& b);

// SLU-14: diag-block solve path -- forward declaration for solve() dispatch
template <class T, class Index>
std::vector<T>
solve_baseline_storage_supernode_aware_diag_block(
    const baseline_lu_storage<T, Index>& lu,
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernode_info,
    const sparse_lu_supernode_numeric<T, Index>& sn_numeric,
    const std::vector<T>& b);

// SLU-8R.3: §17.2(A) supernode-panel update statistics.
// Accumulates counts and timing for trsm + gemm/gemv adapter calls during
// panel updates.  Shared struct: filled by the production true-numeric driver
// and by test helpers.  (Its original consumer, the transitional driver
// run_supernode_panel_leftlooking_update + set_panel_update_info_, was
// removed by SLU-CLN1, 2026-07-05.)
// SLU-8R.3.1: added driver_called, update_applied, dense_kernel_called.
struct supernode_panel_update_stats {
    std::size_t panel_count;         // supernodes with non-empty update sets
    std::size_t update_count;        // total (k,j) pairs applied
    std::size_t trsm_count;          // trsm calls via dense kernel adapter
    std::size_t gemm_count;          // gemm calls via dense kernel adapter
    std::size_t gemv_count;          // gemv calls (w_j==1) via dense kernel adapter
    std::size_t scatter_count;       // scatter-back operations
    std::size_t dense_kernel_ticks;  // nanoseconds for trsm+gemm/gemv
    // SLU-PERF: dense-kernel FLOP accounting (BLAS-implementation-independent).
    double      flop_trsm;           // sum of m*m*n over trsm calls
    double      flop_gemm;           // sum of 2*m*n*k over gemm calls
    double      flop_gemv;           // sum of 2*m*n over gemv calls
    std::size_t gemm_m_sum;          // sum of m (=n_off)  over gemm calls
    std::size_t gemm_n_sum;          // sum of n (=w_j)    over gemm calls
    std::size_t gemm_k_sum;          // sum of k (=n_inter) over gemm calls
    std::size_t gemm_max_m;
    std::size_t gemm_max_n;
    std::size_t gemm_max_k;
    std::size_t gemm_dim_hist[8];    // bucket by min(m,n,k); see sparse_lu_info
    bool used_conservative_reach;    // false: exact U-segment criterion used
    bool symmetric_pruning_active;   // ALWAYS false in SLU-8R.3 (hook only)
    bool driver_called;              // [SLU-CLN1] historical; prototype driver removed (production driver does not set this)
    bool update_applied;             // update_count > 0 && scatter_count > 0
    bool dense_kernel_called;        // trsm > 0 && (gemm + gemv) > 0

    supernode_panel_update_stats()
        : panel_count(0u), update_count(0u), trsm_count(0u),
          gemm_count(0u), gemv_count(0u), scatter_count(0u),
          dense_kernel_ticks(0u),
          flop_trsm(0.0), flop_gemm(0.0), flop_gemv(0.0),
          gemm_m_sum(0u), gemm_n_sum(0u), gemm_k_sum(0u),
          gemm_max_m(0u), gemm_max_n(0u), gemm_max_k(0u),
          used_conservative_reach(false),
          symmetric_pruning_active(false),
          driver_called(false),
          update_applied(false),
          dense_kernel_called(false)
    {
        for (int i = 0; i < 8; ++i) gemm_dim_hist[i] = 0u;
    }
};

// SLU-PERF: gemm shape histogram bucket by smallest dimension min(m,n,k).
//   [0]=1 [1]=2 [2]=3..4 [3]=5..8 [4]=9..16 [5]=17..32 [6]=33..64 [7]=65+
inline int gemm_dim_bucket(std::size_t d) {
    if (d <= 1u)  return 0;
    if (d == 2u)  return 1;
    if (d <= 4u)  return 2;
    if (d <= 8u)  return 3;
    if (d <= 16u) return 4;
    if (d <= 32u) return 5;
    if (d <= 64u) return 6;
    return 7;
}

// SLU-8R.4: §17.2(B) within-panel factorization statistics.
// Defined here (not in the impl file) so that sparse_lu_factorization::set_within_panel_factor_info_
// can use it inline inside the class body.
// Accumulates counts and timing for §17.2(B) pivot search, row swap, scale, and ger calls.
struct within_panel_factor_stats {
    bool executed;           // [SLU-CLN1] set by factorize_within_panel_single (production §17.2(B))
    bool completed;          // all panels completed without zero-pivot abort
    bool is_numeric_source;  // always false (historical sub-phase marker; transitional stats)
    bool used_getrf;         // always false (§25 prohibition on getrf as pivot search)

    std::size_t panel_count;                // supernodes processed
    std::size_t pivot_search_count;         // explicit pivot searches over active rows
    std::size_t pivot_accept_count;         // diagonal accepted by threshold
    std::size_t pivot_reject_count;         // diagonal rejected; max-abs row used
    std::size_t row_swap_count;             // row swaps applied
    std::size_t scale_count;               // columns where L-multiplier scaling applied
    std::size_t ger_count;                 // rank-1 ger updates via dense_kernel adapter
    std::size_t gemm_count;               // batch gemm updates (unused in SLU-8R.4)
    std::size_t scalar_update_count;      // reserved
    std::size_t zero_pivot_count;          // pivot columns with abs(p) <= zero_tolerance
    std::size_t near_zero_pivot_count;     // pivot columns with abs(p) <= near_zero_tolerance
    std::size_t inconclusive_pivot_count;  // non-finite in active column
    std::size_t dense_kernel_ticks;        // nanoseconds for ger adapter calls
    double      flop_ger;                  // SLU-PERF: sum of 2*m*n over ger calls

    // SLU-8R.4.1: worst event seen (priority: zero_pivot > rejected > inconclusive > near_zero > success).
    // != success when any abnormal pivot event occurs.  != not_run when factorization ran.
    within_panel_factor_status worst_status;

    within_panel_factor_stats()
        : executed(false), completed(false), is_numeric_source(false), used_getrf(false),
          panel_count(0u), pivot_search_count(0u), pivot_accept_count(0u),
          pivot_reject_count(0u), row_swap_count(0u), scale_count(0u),
          ger_count(0u), gemm_count(0u), scalar_update_count(0u),
          zero_pivot_count(0u), near_zero_pivot_count(0u),
          inconclusive_pivot_count(0u), dense_kernel_ticks(0u),
          flop_ger(0.0),
          worst_status(within_panel_factor_status::not_run) {}
};

// SLU-8R.5.5: Supernodal true numeric factorization statistics.
// Returned by factorize_supernodal_from_a_eff().
// Defined here (in sparse_lu_detail) so set_true_numeric_info_ can use it inline.
// SLU-GT1 D2: templated on R = real_type of the module scalar so the residual
// report fields carry R (requirement-set arithmetic, no double conversion).
// Only factorization_residual_abs/rel are R; flop_* stay double (counters).
template <class R>
struct supernodal_true_numeric_stats {
    bool attempted;                     // factorize_supernodal_from_a_eff was called
    bool success;                       // A_eff-origin factorization + residual passed
    bool values_initialized_from_A;     // panel_values / U_segments set from A_eff (not CSC)
    bool values_initialized_from_csc_numeric; // false for A-origin path
    bool factorization_residual_checked;
    bool factorization_residual_passed;
    supernodal_true_numeric_status status;

    std::size_t supernodes_processed;
    std::size_t panel_update_count;     // (k,j) update pairs in interleaved §17.2(A)
    std::size_t within_panel_count;     // within-panel factorizations (nsup if success)
    std::size_t trsm_count;             // trsm calls during §17.2(A)
    std::size_t gemm_count;             // gemm calls during §17.2(A)
    std::size_t gemv_count;             // gemv calls during §17.2(A)
    std::size_t ger_count;              // ger calls during §17.2(B)
    std::size_t a_entry_scatter_count;  // entries written during A_eff initialization
    std::size_t a_entry_missing_count;  // structural zero positions in A_eff init
    R factorization_residual_abs;
    R factorization_residual_rel;

    // SLU-8R.6.1: Gate 2 repair — true-numeric path dedicated timing (nanoseconds).
    // Measured by std::chrono::steady_clock in factorize_supernodal_from_a_eff.
    // No forced nonzero: these reflect real call durations only.
    // bridge / prototype timing (dense_kernel_time_ticks, supernode_panel_update_ticks,
    // within_panel_update_ticks) is SEPARATE; these fields are true-numeric-path-only.
    std::size_t total_ticks;                      // wall-clock for entire A_eff factorization
    std::size_t dense_kernel_ticks;               // sum of dense kernel adapter timing in true-numeric
    std::size_t panel_update_dense_kernel_ticks;  // §17.2(A) portion of dense_kernel_ticks
    std::size_t within_panel_dense_kernel_ticks;  // §17.2(B) portion of dense_kernel_ticks

    // SLU-SN-OPT: numeric-time breakdown (nanoseconds, same clock as total_ticks).
    // These partition total_ticks into the factorization phases so axis-1
    // (dense_kernel_time / numeric_time) can be reported against the PURE
    // factorization time (factorization_ticks), excluding the post-factorization
    // verification step (residual_ticks).  Design §17.2 line 894 speaks of
    // "分解時間" (factorization time); the residual sanity check is a separate
    // acceptance/verification phase, not part of the LU factorization itself.
    //   total_ticks = init_scatter_ticks + symbolic_ticks
    //               + panel_nonkernel_ticks + within_nonkernel_ticks
    //               + dense_kernel_ticks + residual_ticks (+ small unattributed)
    //   factorization_ticks = total_ticks - residual_ticks
    std::size_t init_scatter_ticks;     // Step 1: A_eff -> panel/U scatter
    std::size_t symbolic_ticks;         // compute_panel_update_set over all supernodes
    std::size_t panel_nonkernel_ticks;  // §17.2(A) gather/scatter/Z-build (excl dense kernel)
    std::size_t within_nonkernel_ticks; // §17.2(B) gather/pivot/scale/scatter (excl dense kernel)
    std::size_t residual_ticks;         // Step 3: post-factorization residual check
    std::size_t factorization_ticks;    // total_ticks - residual_ticks (pure factorization)

    // SLU-MF3: fine-grained breakdown of panel_nonkernel_ticks (multifrontal
    // assembly). Lets us separate the A_eff scatter, the children extend-add, and
    // the Schur-complement contribution-block emit/copy-out, so the in-place
    // frontal update can be targeted and measured (before/after).
    std::size_t mf_aeff_ticks;          // (a) A_eff scatter into the frontal matrix
    std::size_t mf_extend_add_ticks;    // (b) children contribution-block extend-add
    std::size_t mf_cb_emit_ticks;       // (6) Schur-complement contribution-block emit

    // SLU-PERF: dense-kernel FLOP accounting (see sparse_lu_info for conventions).
    double      flop_trsm;
    double      flop_gemm;
    double      flop_gemv;
    double      flop_ger;
    std::size_t gemm_m_sum;
    std::size_t gemm_n_sum;
    std::size_t gemm_k_sum;
    std::size_t gemm_max_m;
    std::size_t gemm_max_n;
    std::size_t gemm_max_k;
    std::size_t gemm_dim_hist[8];

    supernodal_true_numeric_stats()
        : attempted(false), success(false),
          values_initialized_from_A(false), values_initialized_from_csc_numeric(false),
          factorization_residual_checked(false), factorization_residual_passed(false),
          status(supernodal_true_numeric_status::not_attempted),
          supernodes_processed(0), panel_update_count(0), within_panel_count(0),
          trsm_count(0), gemm_count(0), gemv_count(0), ger_count(0),
          a_entry_scatter_count(0), a_entry_missing_count(0),
          factorization_residual_abs(R(0)), factorization_residual_rel(R(0)),
          total_ticks(0), dense_kernel_ticks(0),
          panel_update_dense_kernel_ticks(0), within_panel_dense_kernel_ticks(0),
          init_scatter_ticks(0), symbolic_ticks(0),
          panel_nonkernel_ticks(0), within_nonkernel_ticks(0),
          residual_ticks(0), factorization_ticks(0),
          mf_aeff_ticks(0), mf_extend_add_ticks(0), mf_cb_emit_ticks(0),
          flop_trsm(0.0), flop_gemm(0.0), flop_gemv(0.0), flop_ger(0.0),
          gemm_m_sum(0u), gemm_n_sum(0u), gemm_k_sum(0u),
          gemm_max_m(0u), gemm_max_n(0u), gemm_max_k(0u)
    {
        for (int i = 0; i < 8; ++i) gemm_dim_hist[i] = 0u;
    }
};

// SLU-8R.5: §18.2 storage-native supernodal solve diagnostics.
// Defined here (in sparse_lu_detail, before class body) so solve() can
// declare a local stats object inline without a separate forward declaration.
// Accumulates counts of dense kernel adapter calls and timing for native solve.
struct supernodal_storage_solve_stats {
    bool        attempted;          // native solve was attempted
    bool        succeeded;          // native solve completed without error
    bool        fallback_csc;       // CSC baseline fallback was used
    std::string fallback_reason;    // reason for fallback (empty if native succeeded)

    std::size_t l_block_count;      // L diagonal blocks solved via trsm adapter
    std::size_t l_update_count;     // off-diagonal L updates (gemv adapter calls)
    std::size_t u_block_count;      // U diagonal blocks solved via trsm adapter
    std::size_t u_seg_apply_count;  // U_segment entries applied (sparse scatter)
    std::size_t trsm_count;         // total trsm adapter calls (L + U blocks)
    std::size_t gemv_count;         // gemv adapter calls (off-diagonal L, nrhs==1)
    std::size_t gemm_count;         // gemm adapter calls (off-diagonal L, nrhs>1)
    std::size_t ticks;              // nanoseconds for native solve
    std::size_t nrhs;               // number of RHS columns solved

    supernodal_storage_solve_stats()
        : attempted(false), succeeded(false), fallback_csc(false),
          l_block_count(0), l_update_count(0), u_block_count(0),
          u_seg_apply_count(0), trsm_count(0), gemv_count(0), gemm_count(0),
          ticks(0), nrhs(0) {}
};

// SLU-8R.5 / B2+ Option A: Forward declarations for §18.2 native solve helpers.
// Implementations in vcp/tsparse/detail/tsparse_sparse_lu_supernodal_solve_impl.hpp.
//
// Native supernodal solve eligibility is based on the accepted A_eff-origin
// storage contract (valid / source_of_truth_storage / true_numeric_source /
// nonempty supernodes).  info_.within_panel_status may describe the transitional
// CSC-bootstrapped §17.2(B) run; it must not gate native solve once
// storage.true_numeric_source has been accepted.
template <class T, class Index>
bool can_use_supernodal_storage_solve(
    const supernodal_lu_storage<T, Index>& storage,
    std::string*                           reason_out);

template <class T, class Index>
std::vector<T>
supernodal_storage_solve_single_rhs(
    const supernodal_lu_storage<T, Index>& storage,
    Index                                  n,
    const std::vector<T>&                  b,
    supernodal_storage_solve_stats&        stats);

} // namespace sparse_lu_detail

// SLU-12: supernode numeric invariant validator -- forward declaration for solve() preconditions
template <class T, class Index>
bool sparse_lu_is_valid_supernode_numeric(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& numeric);

// SLU-12: factor-level CSC consistency check -- forward declaration for solve() preconditions
template <class T, class Index>
bool sparse_lu_verify_supernode_numeric_against_csc(
    const sparse_lu_factorization<T, Index>& fac);

// Forward declaration for the SLU-2 test factory (defined in detail include).
template <class T, class Index>
sparse_lu_factorization<T, Index>
sparse_lu_make_baseline_factor_for_testing(
    Index n,
    const baseline_lu_storage<T, Index>& storage);

// SLU-8R.5: forward declarations for §18.2 solve free functions.
// Defined in vcp/tsparse/detail/tsparse_sparse_lu_supernodal_solve_impl.hpp (included below).
template <class T, class Index>
sparse_lu_factorization<T, Index>
sparse_lu_make_supernodal_factor_for_testing(
    Index n,
    const supernodal_lu_storage<T, Index>& storage);

template <class T, class Index>
std::vector<T>
sparse_lu_supernodal_storage_solve(
    const sparse_lu_factorization<T, Index>&          fac,
    const std::vector<T>&                             b,
    sparse_lu_detail::supernodal_storage_solve_stats* stats_out = 0);

// SLU-8R.5.1: multi-RHS native solve — solves for each column in B.
// Returns X[j] = A^{-1}*B[j] for j in [0, B.size()). stats.nrhs == B.size().
// This is the canonical multi-RHS entry point; not a wrapper over sequential calls.
template <class T, class Index>
std::vector<std::vector<T>>
sparse_lu_supernodal_storage_solve_all_rhs(
    const sparse_lu_factorization<T, Index>&                fac,
    const std::vector<std::vector<T>>&                      B,
    sparse_lu_detail::supernodal_storage_solve_stats*       stats_out = 0);

template <class Matrix>
sparse_lu_symbolic_result<typename Matrix::index_type>
sparse_lu_symbolic(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt);

template <class Matrix>
sparse_lu_factorization<typename Matrix::value_type, typename Matrix::index_type>
sparse_lu_numeric(
    const Matrix& A,
    const sparse_lu_symbolic_result<typename Matrix::index_type>& sym,
    const sparse_lu_options<typename Matrix::value_type>& opt);

template <class Matrix>
sparse_lu_factorization<typename Matrix::value_type, typename Matrix::index_type>
sparse_lu_factorize_with_info(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt);

template <class Matrix>
sparse_lu_factorization<typename Matrix::value_type, typename Matrix::index_type>
sparse_lu_factorize(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt);

// SLU-8R.5.5: forward declaration for A_eff-origin true numeric factorization.
// Defined in vcp/tsparse/detail/tsparse_sparse_lu_true_numeric_impl.hpp (included below).
// Takes A_csc (column-permuted A), csc_lu (for row_perm/Dr/Dc), supernodal storage,
// and opt.  Modifies storage in-place (may-modify contract, SLU-8R.5.5.2 Option B):
//   on success: sets true_numeric_source=true, numeric_source_kind=a_eff_true_numeric.
//   on failure: panel_values/U_segments may be modified (zeroed + A_eff scatter);
//               true_numeric_source remains false, native solve disabled.
template <class T, class Index>
sparse_lu_detail::supernodal_true_numeric_stats<
    typename vcp::tsparse_scalar::real_type<T>::type>
sparse_lu_factorize_supernodal_from_a_eff(
    const csc_storage<T, Index>&           A_csc,
    const baseline_lu_storage<T, Index>&   csc_lu,
    supernodal_lu_storage<T, Index>&       storage,
    const sparse_lu_options<T>&            opt);

// SLU-9: supernode-aware triangular solve helper (forward declaration)
template <class T, class Index>
std::vector<T>
sparse_lu_solve_supernode_aware(
    const sparse_lu_factorization<T, Index>& fac,
    const std::vector<T>& b);

// ===========================================================================
// sparse_lu_factorization<T, Index>
// ===========================================================================

template <class T, class Index>
class sparse_lu_factorization {
public:
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");

    typedef T     value_type;
    typedef Index index_type;

    sparse_lu_factorization()
        : storage_kind_(sparse_lu_storage_kind::baseline_csc),
          uses_supernodal_prototype_(false) {}

    bool valid() const { return info_.success; }

    sparse_lu_storage_kind storage_kind() const { return storage_kind_; }

    const sparse_lu_info<T, Index>& info() const { return info_; }

    // SLU-10: returns true iff this factor was produced via the explicit
    // supernodal prototype path (method=supernodal).  False for all baseline paths.
    bool uses_supernodal_prototype() const { return uses_supernodal_prototype_; }

    // SLU-8R.2: returns true iff factor owns a valid supernodal_lu_storage object.
    // True for method=supernodal after SLU-8R.2 bootstrap.
    // storage_kind()==supernodal iff has_supernodal_storage()==true.
    // Note: has_supernodal_storage==true does NOT imply true_supernodal_numeric.
    bool has_supernodal_storage() const {
        return supernodal_.valid && supernodal_.source_of_truth_storage;
    }

    // SLU-8R.2: returns const reference to the supernodal storage object.
    // Caller must check has_supernodal_storage() before use.
    const supernodal_lu_storage<T, Index>& supernodal_storage() const {
        return supernodal_;
    }

    // SLU-8R.2: returns true iff supernodal storage was bootstrapped from CSC L/U.
    // True in SLU-8R.2 transition; false once §17.2 fills storage directly.
    bool supernodal_storage_bootstrapped_from_csc() const {
        return supernodal_.bootstrapped_from_csc && supernodal_.valid;
    }

    // SLU-8R.5.5: returns true iff A_eff-origin interleaved factorization succeeded
    // and supernodal storage is the numeric source of truth (true_numeric_source == true).
    bool supernodal_storage_is_numeric_source() const {
        return supernodal_.true_numeric_source;
    }

    // SLU-8R.5 / B2+ Option A: returns true iff this factor supports storage-native
    // supernodal solve.  True when valid supernodal storage exists as source_of_truth
    // and storage.true_numeric_source == true.
    // Set by set_supernodal_solve_info_() at factorization time.
    // info_.within_panel_status does NOT affect this flag; that status may describe
    // the transitional CSC-bootstrapped §17.2(B) run and must not gate native solve
    // once A_eff-origin true_numeric_source storage has been accepted.
    bool supernodal_solve_native() const {
        return info_.supernodal_solve_native;
    }

    // SLU-10: returns the numeric supernode metadata built from actual CSC L/U.
    // valid() is false for non-supernodal factors.
    const sparse_lu_supernode_numeric<T, Index>& supernode_numeric_info() const {
        return supernode_numeric_info_;
    }

    // SLU-11: returns the symbolic supernode metadata stored in this factor.
    // Used by validator tests to call sparse_lu_is_valid_supernode_numeric.
    const sparse_lu_supernode_symbolic<Index>& supernode_info() const {
        return supernode_info_;
    }

    // [SLU-RQ1] reach retention removed; derive from supernode_info() on demand.

    // SLU-11: diagnostic-only accessor for the CSC-backed baseline LU storage.
    // Allows verification helpers to access the actual L/U factors.
    const baseline_lu_storage<T, Index>& baseline_storage() const {
        return baseline_;
    }

    // SLU-11: returns diagnostics for the supernode numeric metadata.
    // For baseline/auto factors (uses_supernodal_prototype==false), returns
    // a diagnostics struct with valid==false.
    // Defined out-of-line after sparse_lu_make_supernode_numeric_diagnostics.
    sparse_lu_supernode_numeric_diagnostics<Index> supernode_numeric_diagnostics() const;

    // SLU-12: returns true iff fac.solve(b) will use the supernode-aware CSC path.
    // True for valid explicit supernodal prototype factors with valid numeric metadata.
    // False for baseline/auto factors and failed factors.
    bool solve_uses_supernode_aware_path() const {
        return valid() && uses_supernodal_prototype_ && supernode_numeric_info_.valid;
    }

    // SLU-14: returns true iff fac.solve(b) will use the validated diag-block
    // U local solve prototype.  True exactly when solve_uses_supernode_aware_path()
    // is true: same preconditions apply (valid factor, explicit supernodal, valid
    // numeric metadata).  False for baseline/auto factors and failed factors.
    bool solve_uses_diag_block_path() const {
        return solve_uses_supernode_aware_path();
    }

    // Solve A*x = b using the stored factor.
    //
    // Dispatch priority (SLU-8R.5.5.1):
    //   1. storage_kind_==supernodal && true_numeric_source==true
    //        -> §18.2 storage-native supernodal solve (highest priority)
    //           Bypasses uses_supernodal_prototype_ so that accepted A_eff-origin
    //           storage always uses the native path, bypassing the transitional CSC-backed solve.
    //   2. uses_supernodal_prototype_ && true_numeric_source==false
    //        -> transitional CSC-backed supernode-aware solve (SLU-12/14)
    //   3. storage_kind_==baseline_csc
    //        -> baseline CSC solve (SLU-2)
    //
    // Throws if factor is not valid or no valid solve path exists.
    std::vector<T> solve(const std::vector<T>& b) const {
        if (!valid()) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_factorization::solve: factor is not valid (status=",
                sparse_lu_status_to_string(info_.status), ")");
        }
        // [SLU-CLN2 C1, 2026-07-05] Validate the RHS size once for ALL dispatch
        // priorities.  Priority 1 (supernodal native) previously called the
        // size-unchecked internal helper directly, silently returning an
        // indeterminate solution from an out-of-bounds read when
        // b.size() != n (issue_SLU_native_solve_rhs_size_gap.md).
        if (static_cast<Index>(b.size()) != info_.n) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_factorization::solve: b.size() != n");
        }
        // SLU-8R.5.5.1 DISPATCH PRIORITY 1: §18.2 storage-native supernodal solve.
        // Checked BEFORE uses_supernodal_prototype_ so that accepted A_eff-origin
        // storage (true_numeric_source==true) always uses the native path.
        // can_use_supernodal_storage_solve() rejects transitional storage
        // (true_numeric_source==false), so bootstrap storage falls through below.
        // B2+ Option A: info_.within_panel_status is NOT passed; the transitional
        // §17.2(B) status must not gate native solve once true_numeric_source
        // has been accepted.
        if (storage_kind_ == sparse_lu_storage_kind::supernodal) {
            std::string fallback_reason;
            if (sparse_lu_detail::can_use_supernodal_storage_solve(
                    supernodal_, &fallback_reason)) {
                sparse_lu_detail::supernodal_storage_solve_stats solve_stats;
                return sparse_lu_detail::supernodal_storage_solve_single_rhs(
                    supernodal_, info_.n, b, solve_stats);
            }
            // true_numeric_source==false or invalid storage:
            // fall through to transitional prototype or CSC fallback below.
        }
        // SLU-8R.5.5.1 DISPATCH PRIORITY 2: transitional supernodal prototype.
        // Reached only when true_numeric_source==false (native solve not available)
        // or when storage_kind_!=supernodal.
        // For accepted true-numeric storage (true_numeric_source==true), Priority 1
        // returns before this block is reached.
        if (uses_supernodal_prototype_) {
            // Precondition 1: numeric metadata must be valid
            if (!supernode_numeric_info_.valid) {
                vcp::throw_error<vcp::state_error>(
                    "sparse_lu_factorization::solve: "
                    "supernode numeric metadata is not valid (uses_supernodal_prototype=true)");
            }
            // Precondition 2: metadata must satisfy all structural invariants
            if (!sparse_lu_is_valid_supernode_numeric(
                    info_.n, supernode_info_, supernode_numeric_info_)) {
                vcp::throw_error<vcp::state_error>(
                    "sparse_lu_factorization::solve: "
                    "supernode numeric metadata fails invariant check");
            }
            // Precondition 3: metadata must be consistent with actual CSC L/U
            // (CSC is still the numeric source of truth in transitional path)
            if (!sparse_lu_verify_supernode_numeric_against_csc(*this)) {
                vcp::throw_error<vcp::state_error>(
                    "sparse_lu_factorization::solve: "
                    "supernode numeric metadata is inconsistent with CSC L/U factors");
            }
            // SLU-14: dispatch to diag-block solve (uses validated diag_block_values
            // for U diagonal blocks; CSC data for off-diagonal contributions).
            // Numerically equivalent to solve_baseline_storage_supernode_aware
            // when diag_block_values is CSC-consistent (verified above).
            return sparse_lu_detail::solve_baseline_storage_supernode_aware_diag_block(
                baseline_, info_.n, supernode_info_, supernode_numeric_info_, b);
        }
        // SLU-8R.5.5.1 DISPATCH PRIORITY 3: baseline CSC solve.
        // Also last-resort fallback for supernodal without native solve and no prototype.
        switch (storage_kind_) {
        case sparse_lu_storage_kind::baseline_csc:
            return sparse_lu_detail::solve_baseline_storage(baseline_, info_.n, b);
        case sparse_lu_storage_kind::supernodal:
            if (!baseline_.L.col_ptr.empty()) {
                return sparse_lu_detail::solve_baseline_storage(baseline_, info_.n, b);
            }
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_factorization::solve: supernodal native solve unavailable "
                "and no baseline CSC fallback");
        }
        return std::vector<T>(); // unreachable
    }

    // Solve A*X = B for multiple RHS columns.
    // Throws if factor is not valid or solve is not yet implemented.
    template <class DenseMatrix>
    DenseMatrix solve_multiple_rhs(const DenseMatrix& B) const {
        (void)B;
        if (!valid()) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_factorization::solve_multiple_rhs: factor is not valid (status=",
                sparse_lu_status_to_string(info_.status), ")");
        }
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_factorization::solve_multiple_rhs: not implemented (SLU-9)");
        return DenseMatrix(); // unreachable (throw_error is [[noreturn]])
    }

private:
    sparse_lu_storage_kind              storage_kind_;
    baseline_lu_storage<T, Index>       baseline_;
    supernodal_lu_storage<T, Index>     supernodal_;
    sparse_lu_info<T, Index>            info_;
    sparse_lu_supernode_symbolic<Index> supernode_info_; // SLU-9: stored for supernode-aware solve

    // SLU-10: supernodal prototype fields
    bool                                       uses_supernodal_prototype_;
    sparse_lu_supernode_numeric<T, Index>      supernode_numeric_info_;

    // Package-internal setter used by API free functions
    void set_info_(const sparse_lu_info<T, Index>& i) { info_ = i; }

    // SLU-2 testing hook: injects manually constructed baseline storage.
    // NOT a numeric factorization result; growth_factor is set to 0 sentinel.
    void set_baseline_storage_for_internal_use_(
        Index n,
        const baseline_lu_storage<T, Index>& storage)
    {
        typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
        storage_kind_ = sparse_lu_storage_kind::baseline_csc;
        baseline_     = storage;
        sparse_lu_info<T, Index> info;
        info.success      = true;
        info.status       = sparse_lu_status::success;
        info.n            = n;
        info.nnz_L        = static_cast<Index>(storage.L.row_ind.size());
        info.nnz_U        = static_cast<Index>(storage.U.row_ind.size());
        info.growth_factor = real_type(0);
        // SLU-8R.0: test-injection path -- not a supernodal prototype.
        info.uses_supernodal_prototype = false;
        info.true_supernodal_numeric   = false;
        info.dense_kernel_connected    = false;
        // SLU-8R.2: no supernodal storage for baseline test-injection path.
        info.has_supernodal_storage                   = false;
        info.supernodal_storage_bootstrapped_from_csc = false;
        info.supernodal_storage_is_numeric_source     = false;
        info.supernodal_solve_native                  = false;
        info.supernodal_solve_attempted               = false;
        info.supernodal_solve_fallback_csc            = false;
        info.supernodal_solve_used_csc_factor_data    = false;
        info.supernodal_solve_used_injected_storage   = false;
        info.supernodal_solve_status_value            = supernodal_solve_status::not_attempted;
        info.supernodal_panel_value_count             = 0;
        info.supernodal_u_segment_count               = 0;
        info.supernodal_storage_bytes                 = 0;
        // [SLU-CLN1 C1] transitional §17.2(A)/(B) prototype fields removed.
        info_ = info;
    }

    // SLU-9: store symbolic supernode metadata for supernode-aware solve.
    void set_supernode_info_(const sparse_lu_supernode_symbolic<Index>& si) {
        supernode_info_ = si;
    }

    // SLU-10: store supernodal prototype numeric metadata and mark factor accordingly.
    // Must be called AFTER set_baseline_storage_with_diagnostics_ so info_ is valid.
    //
    // SLU-8R.1: Added kernel_called and kernel_ticks parameters for production
    //   dense kernel connection tracking.
    //
    //   dense_kernel_connected = kernel_called:
    //     true iff production path called sparse_lu_dense_kernel<T> on real data.
    //     Does NOT imply true_supernodal_numeric (historical marker; always false).
    //     Does NOT imply genuine A_eff-origin factorization (check supernodal_true_numeric_success).
    //     Does NOT imply §18.2 storage-native solve in use (check supernodal_solve_native).
    //
    //   true_supernodal_numeric remains false: historical conformance marker; always false.
    //   uses_supernodal_prototype remains true: supernodal path was entered.
    void set_supernodal_prototype_info_(
        const sparse_lu_supernode_numeric<T, Index>& sn_num,
        bool kernel_called = false,
        std::size_t kernel_ticks = 0)
    {
        uses_supernodal_prototype_ = true;
        supernode_numeric_info_    = sn_num;
        info_.method_used             = sparse_lu_method::supernodal;
        // SLU-8R.1: prototype status flags.
        info_.uses_supernodal_prototype = true;
        info_.true_supernodal_numeric   = false; // historical conformance marker; always false
        // SLU-8R.1: dense_kernel_connected = true only if production adapter was called.
        // At this point, CSC baseline is numeric source; A_eff-origin step follows later.
        info_.dense_kernel_connected    = kernel_called;
        info_.dense_kernel_time_ticks   = kernel_ticks;
        // number_of_supernodes from actual symbolic partition (not a true supernodal count).
        if (sn_num.valid && sn_num.supernode_ptr.size() >= 1u) {
            info_.number_of_supernodes =
                static_cast<Index>(sn_num.supernode_ptr.size()) - Index(1);
        }
    }

    // SLU-8R.2: store bootstrapped supernodal storage and update diagnostics.
    // Must be called AFTER set_supernodal_prototype_info_ so info_ is already valid.
    //
    // Sets storage_kind_=supernodal if sn_storage.valid.
    // Updates has_supernodal_storage, supernodal_storage_bootstrapped_from_csc,
    // supernodal_storage_is_numeric_source, supernodal_solve_native, and size fields.
    //
    // SLU-8R.2 constraints:
    //   sn_storage.true_numeric_source is false at this bootstrap stage.
    //   supernodal_solve_native is finalized later by set_supernodal_solve_info_().
    void set_supernodal_storage_info_(
        const supernodal_lu_storage<T, Index>& sn_storage)
    {
        supernodal_ = sn_storage;
        if (sn_storage.valid && sn_storage.source_of_truth_storage) {
            storage_kind_ = sparse_lu_storage_kind::supernodal;
        }
        info_.has_supernodal_storage =
            sn_storage.valid && sn_storage.source_of_truth_storage;
        info_.supernodal_storage_bootstrapped_from_csc = sn_storage.bootstrapped_from_csc;
        // At this bootstrap stage, true_numeric_source is false; storage values come
        // from CSC L/U.  A_eff-origin step (SLU-8R.5.5) sets true_numeric_source.
        info_.supernodal_storage_is_numeric_source = sn_storage.true_numeric_source;
        // SLU-8R.5: supernodal_solve_native is finalized by set_supernodal_solve_info_().
        // For production transitional storage (true_numeric_source==false), it will be false.
        info_.supernodal_solve_native                 = false; // finalized by set_supernodal_solve_info_()
        info_.supernodal_solve_attempted              = false; // finalized by set_supernodal_solve_info_()
        info_.supernodal_solve_fallback_csc           = false; // finalized by set_supernodal_solve_info_()
        info_.supernodal_solve_used_csc_factor_data   = false; // finalized by set_supernodal_solve_info_()
        info_.supernodal_solve_used_injected_storage  = false; // production path
        info_.supernodal_solve_status_value = supernodal_solve_status::not_attempted; // finalized by set_supernodal_solve_info_()
        info_.supernodal_panel_value_count = sn_storage.panel_values.size();
        info_.supernodal_u_segment_count   = sn_storage.U_segments.row_ind.size();
        // Approximate byte footprint: panel_values + U_segments + descriptors.
        std::size_t bytes = sn_storage.panel_values.size() * sizeof(T);
        bytes += sn_storage.U_segments.seg_ptr.size() * sizeof(Index);
        bytes += sn_storage.U_segments.row_ind.size()  * sizeof(Index);
        bytes += sn_storage.U_segments.values.size()   * sizeof(T);
        bytes += sn_storage.supernodes.size() * sizeof(supernode_desc<Index>);
        bytes += sn_storage.row_perm.size()     * sizeof(Index);
        bytes += sn_storage.inv_row_perm.size() * sizeof(Index);
        bytes += sn_storage.col_perm.size()     * sizeof(Index);
        bytes += sn_storage.inv_col_perm.size() * sizeof(Index);
        bytes += sn_storage.Dr.size() * sizeof(T);
        bytes += sn_storage.Dc.size() * sizeof(T);
        info_.supernodal_storage_bytes = bytes;
    }

    // [SLU-CLN1 C1, 2026-07-05] set_panel_update_info_ /
    // set_within_panel_factor_info_ (transitional §17.2(A)/(B) prototype
    // diagnostics capture) REMOVED together with the prototype pass and the
    // prototype-only info fields.  Production §17.2 diagnostics are the
    // supernodal_true_numeric_* counters (set_true_numeric_info_ below).

    // SLU-8R.5.5: store true-numeric factorization diagnostics.
    // Called AFTER factorize_supernodal_from_a_eff() returns, passing the
    // modified sn_storage so that supernodal_ is refreshed with the
    // updated true_numeric_source / numeric_source_kind flags.
    // Must be called BEFORE set_supernodal_solve_info_() so that the
    // supernodal_.true_numeric_source is visible to the solve gate check.
    void set_true_numeric_info_(
        const sparse_lu_detail::supernodal_true_numeric_stats<
            typename vcp::tsparse_scalar::real_type<T>::type>& tns,
        const supernodal_lu_storage<T, Index>&                 updated_storage)
    {
        // Refresh supernodal_ to pick up true_numeric_source / numeric_source_kind.
        supernodal_ = updated_storage;
        // Re-derive storage_is_numeric_source from the refreshed storage.
        if (supernodal_.valid && supernodal_.source_of_truth_storage) {
            info_.supernodal_storage_is_numeric_source = supernodal_.true_numeric_source;
        }
        info_.supernodal_true_numeric_attempted               = tns.attempted;
        info_.supernodal_true_numeric_success                 = tns.success;
        info_.supernodal_values_initialized_from_A            = tns.values_initialized_from_A;
        info_.supernodal_values_initialized_from_csc_numeric  = tns.values_initialized_from_csc_numeric;
        info_.supernodal_factorization_residual_checked       = tns.factorization_residual_checked;
        info_.supernodal_factorization_residual_passed        = tns.factorization_residual_passed;
        info_.supernodal_numeric_source                       = supernodal_.numeric_source_kind;
        info_.supernodal_true_numeric_status_value            = tns.status;
        info_.supernodal_true_numeric_supernodes_processed    = tns.supernodes_processed;
        info_.supernodal_true_numeric_panel_update_count      = tns.panel_update_count;
        info_.supernodal_true_numeric_within_panel_count      = tns.within_panel_count;
        info_.supernodal_true_numeric_trsm_count              = tns.trsm_count;
        info_.supernodal_true_numeric_gemm_count              = tns.gemm_count;
        info_.supernodal_true_numeric_gemv_count              = tns.gemv_count;
        info_.supernodal_true_numeric_ger_count               = tns.ger_count;
        info_.supernodal_factorization_residual_abs           = tns.factorization_residual_abs;
        info_.supernodal_factorization_residual_rel           = tns.factorization_residual_rel;
        // SLU-8R.6.1: true-numeric path dedicated timing (separated from bridge/prototype).
        info_.supernodal_true_numeric_total_ticks             = tns.total_ticks;
        info_.supernodal_true_numeric_dense_kernel_ticks      = tns.dense_kernel_ticks;
        info_.supernodal_true_numeric_panel_update_ticks      = tns.panel_update_dense_kernel_ticks;
        info_.supernodal_true_numeric_within_panel_ticks      = tns.within_panel_dense_kernel_ticks;
        // SLU-SN-OPT: numeric-time phase breakdown.
        info_.supernodal_true_numeric_init_scatter_ticks      = tns.init_scatter_ticks;
        info_.supernodal_true_numeric_symbolic_ticks          = tns.symbolic_ticks;
        info_.supernodal_true_numeric_panel_nonkernel_ticks   = tns.panel_nonkernel_ticks;
        info_.supernodal_true_numeric_within_nonkernel_ticks  = tns.within_nonkernel_ticks;
        info_.supernodal_true_numeric_residual_ticks          = tns.residual_ticks;
        info_.supernodal_true_numeric_factorization_ticks     = tns.factorization_ticks;
        info_.supernodal_true_numeric_mf_aeff_ticks           = tns.mf_aeff_ticks;
        info_.supernodal_true_numeric_mf_extend_add_ticks     = tns.mf_extend_add_ticks;
        info_.supernodal_true_numeric_mf_cb_emit_ticks        = tns.mf_cb_emit_ticks;

        // SLU-PERF (design §17.2 line 894): fill the generic numeric_time_ticks
        // with the A_eff-origin numeric-factorization wall-clock so that the
        // time-ratio dense_kernel_time_ticks/numeric_time_ticks can be formed in
        // the SAME unit (nanoseconds) and SAME clock (steady_clock) used for the
        // dense-kernel timing. Both come from factorize_supernodal_from_a_eff.
        // numeric_time_ticks was previously left at 0 for the supernodal path.
        if (tns.attempted) {
            info_.numeric_time_ticks = tns.total_ticks;
        }
        // SLU-PERF: dense-kernel FLOP location + gemm call-shape distribution.
        info_.supernodal_true_numeric_flop_trsm  = tns.flop_trsm;
        info_.supernodal_true_numeric_flop_gemm  = tns.flop_gemm;
        info_.supernodal_true_numeric_flop_gemv  = tns.flop_gemv;
        info_.supernodal_true_numeric_flop_ger   = tns.flop_ger;
        info_.supernodal_true_numeric_gemm_m_sum = tns.gemm_m_sum;
        info_.supernodal_true_numeric_gemm_n_sum = tns.gemm_n_sum;
        info_.supernodal_true_numeric_gemm_k_sum = tns.gemm_k_sum;
        info_.supernodal_true_numeric_gemm_max_m = tns.gemm_max_m;
        info_.supernodal_true_numeric_gemm_max_n = tns.gemm_max_n;
        info_.supernodal_true_numeric_gemm_max_k = tns.gemm_max_k;
        for (int hi = 0; hi < 8; ++hi)
            info_.supernodal_true_numeric_gemm_dim_hist[hi] = tns.gemm_dim_hist[hi];
    }

    // SLU-8R.5 / B2+ Option A: finalize supernodal solve diagnostic flags after
    // all §17.2 info is set.  Must be called AFTER set_within_panel_factor_info_()
    // so info_.within_panel_status is up-to-date.
    //
    // SLU-8R.5 / B2+ Option A contract:
    //   supernodal_solve_native = true when:
    //     - valid supernodal storage exists as source_of_truth
    //     - storage.true_numeric_source == true (A_eff-origin factorization accepted)
    //
    //   info_.within_panel_status may describe the transitional CSC-bootstrapped
    //   §17.2(B) run and must NOT gate native solve once A_eff-origin
    //   true_numeric_source storage has been accepted.  A_eff-origin pivot failure
    //   is already reflected in storage.true_numeric_source == false.
    //
    // For production transitional storage (true_numeric_source == false):
    //   supernodal_solve_native = false  (CSC baseline fallback will be used)
    //   supernodal_solve_fallback_csc = true
    //
    // For A_eff-origin accepted storage (true_numeric_source == true):
    //   supernodal_solve_native = true
    //   supernodal_solve_fallback_csc = false
    //
    // Gate 6: tied to true_numeric_source, not merely to supernodal_solve_native.
    void set_supernodal_solve_info_() {
        const bool is_supernodal = (storage_kind_ == sparse_lu_storage_kind::supernodal);
        info_.supernodal_solve_attempted = is_supernodal;
        info_.supernodal_solve_used_injected_storage = false;

        if (!is_supernodal) {
            info_.supernodal_solve_native           = false;
            info_.supernodal_solve_fallback_csc     = false;
            info_.supernodal_solve_used_csc_factor_data = false;
            info_.supernodal_solve_status_value = supernodal_solve_status::not_attempted;
            return;
        }

        // Check storage validity
        if (!info_.has_supernodal_storage ||
            !supernodal_.valid ||
            !supernodal_.source_of_truth_storage) {
            info_.supernodal_solve_native           = false;
            info_.supernodal_solve_fallback_csc     = true;
            info_.supernodal_solve_used_csc_factor_data = true;
            info_.supernodal_solve_status_value = supernodal_solve_status::invalid_storage;
            return;
        }

        // SLU-8R.5 key gate: production transitional storage must fallback
        if (!supernodal_.true_numeric_source) {
            info_.supernodal_solve_native           = false;
            info_.supernodal_solve_fallback_csc     = true;
            info_.supernodal_solve_used_csc_factor_data = true;
            info_.supernodal_solve_status_value =
                supernodal_solve_status::true_numeric_source_false;
            return;
        }

        // NOTE: info_.within_panel_status reflects the TRANSITIONAL §17.2(B) run
        // (on CSC-bootstrapped storage, executed before the A_eff-origin run).
        // By this point, true_numeric_source == true, meaning the A_eff-origin
        // factorization succeeded without any fatal pivot failure.  The A_eff
        // §17.2(B) validates via residual; any zero_pivot seen in the transitional
        // run does not apply to the A_eff-origin storage and must not gate native
        // solve.  The transitional status checks are therefore omitted here.
        // (T-generic fix: kv::mpfr<N> high-precision arithmetic can produce exact
        // zero in the CSC-bootstrapped transitional path while the A_eff path
        // succeeds; blocking native solve on stale transitional data is incorrect.)

        // All conditions met: native solve will be used
        info_.supernodal_solve_native               = true;
        info_.supernodal_solve_fallback_csc         = false;
        info_.supernodal_solve_used_csc_factor_data = false;
        info_.supernodal_solve_status_value = supernodal_solve_status::storage_native_success;
    }

    // SLU-8R.5 testing hook: injects manually constructed supernodal storage.
    // NOT a numeric factorization result; use only in tests.
    //
    // Design intent: the injected storage is ASSUMED to contain correct L/U values.
    // Therefore, true_numeric_source is FORCED to true (regardless of what the
    // caller set in storage.true_numeric_source), allowing native solve.
    //
    // Sets baseline_ to empty (no CSC fallback).
    // Sets supernodal_solve_native = true for valid storage with non-empty supernodes.
    // Sets supernodal_solve_used_injected_storage = true.
    // Called by sparse_lu_make_supernodal_factor_for_testing (in impl file).
    void set_supernodal_factor_for_testing_(
        Index n,
        const supernodal_lu_storage<T, Index>& storage)
    {
        typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
        supernodal_    = storage;
        // Force true_numeric_source = true for test factory (storage is assumed correct).
        // This enables can_use_supernodal_storage_solve() to return true, allowing
        // the native §18.2 solve path to be tested.
        supernodal_.true_numeric_source = true;
        storage_kind_  = sparse_lu_storage_kind::supernodal;
        uses_supernodal_prototype_ = false;
        // baseline_ remains default (empty CSC -- no fallback for test factory)
        sparse_lu_info<T, Index> info;
        info.success       = true;
        info.status        = sparse_lu_status::success;
        info.n             = n;
        info.growth_factor = real_type(0);
        // SLU-8R.5 test factory: not a supernodal prototype.
        info.uses_supernodal_prototype = false;
        info.true_supernodal_numeric   = false; // historical conformance marker; always false
        info.dense_kernel_connected    = false;
        // Supernodal storage flags: use forced true_numeric_source
        info.has_supernodal_storage =
            supernodal_.valid && supernodal_.source_of_truth_storage;
        info.supernodal_storage_bootstrapped_from_csc = supernodal_.bootstrapped_from_csc;
        info.supernodal_storage_is_numeric_source     = supernodal_.true_numeric_source; // true (forced)
        info.supernodal_panel_value_count             = supernodal_.panel_values.size();
        info.supernodal_u_segment_count               = supernodal_.U_segments.row_ind.size();
        info.supernodal_storage_bytes                 = 0;
        // SLU-8R.5: native solve available because true_numeric_source is forced true
        const bool valid_for_native =
            supernodal_.valid && supernodal_.source_of_truth_storage &&
            !supernodal_.supernodes.empty() && supernodal_.true_numeric_source;
        info.supernodal_solve_native              = valid_for_native;
        info.supernodal_solve_attempted           = true; // test factory always supernodal
        info.supernodal_solve_fallback_csc        = !valid_for_native;
        info.supernodal_solve_used_csc_factor_data = !valid_for_native;
        info.supernodal_solve_used_injected_storage = true; // test factory marker
        // SLU-8R.5.1: fallback reason enum for injected storage
        info.supernodal_solve_status_value = valid_for_native
            ? supernodal_solve_status::storage_native_success
            : supernodal_solve_status::invalid_storage;
        // [SLU-CLN1 C1] transitional §17.2(A)/(B) prototype fields removed.
        info_ = info;
    }

    // SLU-4.1: setter for production factorization results.
    // Accepts precomputed growth_factor (max|U| / max|A_eff|).
    void set_baseline_storage_with_diagnostics_(
        Index n,
        const baseline_lu_storage<T, Index>& storage,
        typename vcp::tsparse_scalar::real_type<T>::type growth_factor_val)
    {
        storage_kind_ = sparse_lu_storage_kind::baseline_csc;
        baseline_     = storage;
        sparse_lu_info<T, Index> info;
        info.success       = true;
        info.status        = sparse_lu_status::success;
        // SLU-L1 L-3: the baseline GP numeric produced this factor; report it.
        // For method=supernodal the later set_supernodal_prototype_info_ call
        // overwrites this with sparse_lu_method::supernodal (existing behavior).
        info.method_used   = sparse_lu_method::baseline_gp;
        info.n             = n;
        info.nnz_L         = static_cast<Index>(storage.L.row_ind.size());
        info.nnz_U         = static_cast<Index>(storage.U.row_ind.size());
        info.growth_factor = growth_factor_val;
        // SLU-8R.0: baseline/auto path -- not a supernodal prototype.
        info.uses_supernodal_prototype = false;
        info.true_supernodal_numeric   = false; // historical conformance marker; always false
        info.dense_kernel_connected    = false; // baseline path does not call dense kernel
        // SLU-8R.2: no supernodal storage for baseline/auto path.
        info.has_supernodal_storage                   = false;
        info.supernodal_storage_bootstrapped_from_csc = false;
        info.supernodal_storage_is_numeric_source     = false;
        info.supernodal_solve_native                  = false;
        info.supernodal_solve_attempted               = false;
        info.supernodal_solve_fallback_csc            = false;
        info.supernodal_solve_used_csc_factor_data    = false;
        info.supernodal_solve_used_injected_storage   = false;
        info.supernodal_solve_status_value            = supernodal_solve_status::not_attempted;
        info.supernodal_panel_value_count             = 0;
        info.supernodal_u_segment_count               = 0;
        info.supernodal_storage_bytes                 = 0;
        // [SLU-CLN1 C1] transitional §17.2(A)/(B) prototype fields removed.
        info_ = info;
    }

    // Friend declarations for all API entry points
    template <class Matrix>
    friend sparse_lu_factorization<typename Matrix::value_type,
                                   typename Matrix::index_type>
    sparse_lu_numeric(
        const Matrix&,
        const sparse_lu_symbolic_result<typename Matrix::index_type>&,
        const sparse_lu_options<typename Matrix::value_type>&);

    template <class Matrix>
    friend sparse_lu_factorization<typename Matrix::value_type,
                                   typename Matrix::index_type>
    sparse_lu_factorize_with_info(
        const Matrix&,
        const sparse_lu_options<typename Matrix::value_type>&);

    // SLU-2 test factory (all instantiations are friends)
    template <class T2, class Idx2>
    friend sparse_lu_factorization<T2, Idx2>
    sparse_lu_make_baseline_factor_for_testing(
        Idx2 n,
        const baseline_lu_storage<T2, Idx2>& storage);

    // SLU-9: supernode-aware solve helper needs access to private storage
    template <class T2, class Idx2>
    friend std::vector<T2>
    sparse_lu_solve_supernode_aware(
        const sparse_lu_factorization<T2, Idx2>& fac,
        const std::vector<T2>& b);

    // SLU-8R.5: test factory and explicit native-solve free function
    template <class T2, class Idx2>
    friend sparse_lu_factorization<T2, Idx2>
    sparse_lu_make_supernodal_factor_for_testing(
        Idx2 n,
        const supernodal_lu_storage<T2, Idx2>& storage);

    template <class T2, class Idx2>
    friend std::vector<T2>
    sparse_lu_supernodal_storage_solve(
        const sparse_lu_factorization<T2, Idx2>&          fac,
        const std::vector<T2>&                             b,
        sparse_lu_detail::supernodal_storage_solve_stats* stats_out);
};

// ===========================================================================
// Internal: option validation helpers
// ===========================================================================

namespace sparse_lu_detail {

// validate_symbolic_options: checks pivot params, ordering, pivoting, method,
// and future feature flags.
// Returns invalid_input for invalid pivot_threshold / absolute_pivot_tolerance.
// Returns not_implemented for explicitly requested unimplemented features
// (equilibration, compute_condition_estimate).  iterative_refinement (O4.1a) is
// accepted: it is a solve-time post-process and does not affect factorization.
// rcm/amd/colamd orderings are all implemented
// column pre-permutations and are accepted.
// Explicit supernodal is allowed for all orderings (genuine numeric source via SLU-8R.5.5).
// SLU-SP1: method == supernode_panel is ACCEPTED here (the symbolic phase --
// column ordering / etree -- is shared with the other methods); the Phase 1
// skeleton answers not_implemented at the NUMERIC dispatch in
// sparse_lu_numeric / sparse_lu_factorize_with_info, replaced by the real
// driver in Phase 2/3.
// Returns success for natural/auto_select/baseline_gp + valid params + defaults.
template <class T>
inline sparse_lu_status validate_symbolic_options(
    const sparse_lu_options<T>& opt)
{
    typedef typename sparse_lu_options<T>::real_type Real;
    // Validate numeric pivot parameters first: invalid_input before not_implemented.
    if (opt.pivot_threshold < Real(0) || opt.pivot_threshold > Real(1)) {
        return sparse_lu_status::invalid_input;
    }
    if (opt.absolute_pivot_tolerance < Real(0)) {
        return sparse_lu_status::invalid_input;
    }
    // Method and ordering validation.
    // Ordering Track O1/O2/O3: rcm + amd + colamd are un-stubbed column
    // pre-permutations (§2A).  The earlier "colamd out of scope per roadmap §2.1"
    // note is corrected/withdrawn by the O3/O4 roadmap §0; colamd is now a
    // pattern-only A^T A column ordering.  All ordering modes are accepted for
    // both the supernodal and baseline branches; no ordering mode is rejected.
    // Pivoting Track O4: static_mc64 is an implemented, opt-in static pivoting
    // (MC64 maximum-weight matching + Dr/Dc scaling -> zero-free diagonal, then
    // the existing threshold partial pivoting refines on the matched matrix).
    // The default pivoting (threshold_partial) is untouched (S-C).  No pivoting
    // mode is rejected here.
    // O4.2a: equilibration is implemented for baseline_gp / natural|amd|colamd /
    // threshold_partial.  It is incompatible with static_mc64 because both paths
    // own Dr/Dc in storage; the composition rule is undefined until a later task.
    // Default (false) passes through; static_mc64 combination is rejected.
    if (opt.equilibration &&
        opt.pivoting == sparse_lu_pivoting::static_mc64) {
        return sparse_lu_status::invalid_input;
    }
    // O4.1a: iterative_refinement is un-stubbed.  It is a SOLVE-TIME post-process
    // (sparse_lu_solve_refined) that reuses the factorization and does NOT change
    // the factorization or the default solve path; the flag is therefore accepted
    // for every pivoting mode (incl. the default threshold_partial).
    // Only the tolerance is range-checked, and ONLY when IR is enabled: a negative
    // target residual is a meaningless option contract.  tol == 0 stays valid
    // (unreachable but legitimate); the false (default) path is never rejected,
    // preserving the byte-identical default behaviour.
    if (opt.iterative_refinement &&
        opt.iterative_refinement_tolerance < Real(0)) {
        return sparse_lu_status::invalid_input;
    }
    if (opt.compute_condition_estimate) {
        return sparse_lu_status::not_implemented;
    }
    return sparse_lu_status::success;
}

// validate_options: comprehensive check for both symbolic and numeric path.
// Delegates to validate_symbolic_options which covers all SLU-4 constraints.
template <class T>
inline sparse_lu_status validate_options(const sparse_lu_options<T>& opt) {
    return validate_symbolic_options(opt);
}

// ---------------------------------------------------------------------------
// sparse_lu_validate_csc_pattern_for_etree
//
// Validates CSC pattern before elimination tree computation.
// Checks: n >= 0, col_ptr size, col_ptr[0]==0, monotonicity, row index range.
// Throws vcp::invalid_argument on any malformed input; never causes UB.
// ---------------------------------------------------------------------------
template <class Index>
void sparse_lu_validate_csc_pattern_for_etree(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_validate_csc_pattern_for_etree: Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_elimination_tree_csc: n must be non-negative");
    }
    const std::size_t un = static_cast<std::size_t>(n);
    if (col_ptr.size() != un + 1u) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_elimination_tree_csc: col_ptr size must be n+1");
    }
    if (col_ptr[0] != Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_elimination_tree_csc: col_ptr[0] must be zero");
    }
    for (std::size_t j = 0u; j < un; ++j) {
        if (col_ptr[j] < Index(0) || col_ptr[j + 1u] < Index(0)) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_elimination_tree_csc: negative col_ptr entry");
        }
        if (col_ptr[j] > col_ptr[j + 1u]) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_elimination_tree_csc: non-monotone col_ptr");
        }
    }
    const std::size_t nnz = static_cast<std::size_t>(col_ptr[un]);
    if (nnz > row_ind.size()) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_elimination_tree_csc: col_ptr[n] exceeds row_ind.size()");
    }
    for (std::size_t p = 0u; p < nnz; ++p) {
        if (row_ind[p] < Index(0) || row_ind[p] >= n) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_elimination_tree_csc: row index out of [0,n) range");
        }
    }
}

// ---------------------------------------------------------------------------
// sparse_lu_validate_etree_reach_inputs
//
// Validates inputs for sparse_lu_column_reach_from_etree before any vector
// indexing.  Checks: n >= 0, 0 <= j <= n, parent.size() == n,
// parent[k] == -1 or k < parent[k] < n for all k, 0 <= seed < n for all seeds.
// Throws vcp::invalid_argument on malformed input; never causes UB.
// ---------------------------------------------------------------------------
template <class Index>
void sparse_lu_validate_etree_reach_inputs(
    Index n,
    Index j,
    const std::vector<Index>& parent,
    const std::vector<Index>& active_rows)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_validate_etree_reach_inputs: Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_column_reach_from_etree: n must be non-negative");
    }
    if (j < Index(0) || j > n) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_column_reach_from_etree: j must satisfy 0 <= j <= n");
    }
    const std::size_t un = static_cast<std::size_t>(n);
    if (parent.size() != un) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_column_reach_from_etree: parent size must equal n");
    }
    for (Index k = Index(0); k < n; ++k) {
        const std::size_t sk = static_cast<std::size_t>(k);
        const Index p = parent[sk];
        if (!(p == Index(-1) || (k < p && p < n))) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_column_reach_from_etree: invalid parent link"
                " (must be -1 or k < parent[k] < n)");
        }
    }
    for (std::size_t ai = 0u; ai < active_rows.size(); ++ai) {
        const Index seed = active_rows[ai];
        if (seed < Index(0) || seed >= n) {
            vcp::throw_error<vcp::invalid_argument>(
                "sparse_lu_column_reach_from_etree: active row out of [0,n) range");
        }
    }
}

} // namespace sparse_lu_detail

// ===========================================================================
// SLU-1: Input conversion and owning storage helpers
// Included here, inside namespace vcp, after all storage types are defined.
// The detail implementation header has no namespace vcp wrapper and is
// injected into this namespace by this include.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_convert_impl.hpp>

// ===========================================================================
// Ordering Track O1: RCM fill-reducing ordering (pattern-only, deterministic).
// Injected here after convert helpers (csc_storage / permutation helpers).
// Produces a column permutation Q for the symbolic col_perm slot (§2A).
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_ordering_impl.hpp>

// ===========================================================================
// SLU-2: Baseline triangular solve helpers and test factory
// Injected here after convert helpers are available.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_solve_impl.hpp>

// ===========================================================================
// SLU-4: Baseline numeric factorization helper
// Injected here after SLU-2 solve helpers and scalar policy are in scope.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_numeric_impl.hpp>

// ===========================================================================
// SLU-SP1: opt-in left-looking supernode-panel numeric (method =
// supernode_panel).  Injected here with the same prerequisite set as the
// GP numeric helper above (csc_storage, baseline_lu_storage, scalar policy,
// pivot acceptability, permutation helpers).
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_supernode_panel_impl.hpp>

// ===========================================================================
// O4: MC64 / static pivoting (zero-free diagonal; saddle-point solvability).
// Injected after the numeric helpers; uses csc_storage, sparse_lu_scalar_policy,
// and the identity/inverse permutation helpers.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_mc64_impl.hpp>

// ===========================================================================
// O4.2a: Two-sided max-norm equilibration helper.
// Injected after mc64; uses csc_storage, tsparse_scalar, convert helpers.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_equilibration_impl.hpp>

// ===========================================================================
// SLU-6: Symbolic elimination tree and column reachability helpers (public)
// These are standalone symbolic utilities; production GP uses the internal
// L-structure DFS (sparse_lu_column_reach_from_L in sparse_lu_detail).
// ===========================================================================

// ---------------------------------------------------------------------------
// sparse_lu_elimination_tree_csc
//
// Computes the column elimination tree of A^T A from the CSC sparsity pattern.
// parent[k] = -1  (root, column k has no parent)
//           = p   where k < p < n (column p is the parent of column k).
//
// Algorithm (correctness-first, C++11):
//   1. Validate CSC pattern; throw invalid_argument on malformed input.
//   2. Build row -> sorted-unique-columns incidence (rows_to_cols[r]).
//   3. For each column k in order 0..n-1:
//        For each row r in column k:
//          For each column j < k sharing row r (from rows_to_cols[r]):
//            Path-compress j up toward k in the ancestor forest.
//            When an unset ancestor is found, set parent[root]=k and
//            ancestor[root]=k (link this tree to k).
// This is the standard column-etree algorithm with Tarjan path compression,
// applied to the structural pattern of A^T A (no explicit A^T A is formed).
//
// Invariant: parent[k] == -1  OR  k < parent[k] < n  (natural-order forest).
// Result is deterministic for any given CSC pattern.
//
// Throws invalid_argument on malformed CSC input (via validate helper).
// Template parameter Index must be signed (sentinel -1 used for roots).
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index>
sparse_lu_elimination_tree_csc(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_elimination_tree_csc: Index must be signed");

    // Validate input before any access (Finding 2 fix)
    sparse_lu_detail::sparse_lu_validate_csc_pattern_for_etree(n, col_ptr, row_ind);

    const std::size_t un = static_cast<std::size_t>(n);

    if (n == Index(0)) {
        return std::vector<Index>();
    }

    // Step 1: Build row -> sorted unique columns incidence.
    std::vector<std::vector<Index> > rows_to_cols(un);
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index p = col_ptr[sj]; p < col_ptr[sj + 1u]; ++p) {
            const Index r = row_ind[static_cast<std::size_t>(p)];
            rows_to_cols[static_cast<std::size_t>(r)].push_back(j);
        }
    }
    for (std::size_t r = 0u; r < un; ++r) {
        std::sort(rows_to_cols[r].begin(), rows_to_cols[r].end());
        rows_to_cols[r].erase(
            std::unique(rows_to_cols[r].begin(), rows_to_cols[r].end()),
            rows_to_cols[r].end());
    }

    // Step 2: Compute A^T A column etree with Tarjan path compression.
    // ancestor[k] = current compressed root of tree rooted at or above k.
    std::vector<Index> parent(un, Index(-1));
    std::vector<Index> ancestor(un, Index(-1));

    for (Index k = Index(0); k < n; ++k) {
        const std::size_t sk = static_cast<std::size_t>(k);
        // For each row r in column k:
        for (Index p = col_ptr[sk]; p < col_ptr[sk + 1u]; ++p) {
            const Index r = row_ind[static_cast<std::size_t>(p)];
            const std::vector<Index>& cols =
                rows_to_cols[static_cast<std::size_t>(r)];
            // For each column j < k also containing row r,
            // link j's tree component to k in the etree forest.
            for (std::size_t b = 0u; b < cols.size() && cols[b] < k; ++b) {
                Index i = cols[b];
                // Walk ancestor chain from i toward the representative;
                // compress path by pointing intermediate nodes directly at k.
                while (ancestor[static_cast<std::size_t>(i)] != Index(-1) &&
                       ancestor[static_cast<std::size_t>(i)] != k) {
                    const Index next = ancestor[static_cast<std::size_t>(i)];
                    ancestor[static_cast<std::size_t>(i)] = k; // path compression
                    i = next;
                }
                if (ancestor[static_cast<std::size_t>(i)] == Index(-1)) {
                    // i is the root of its component; attach to k
                    parent[static_cast<std::size_t>(i)]   = k;
                    ancestor[static_cast<std::size_t>(i)] = k;
                }
                // else ancestor[i]==k already: component already linked to k
            }
        }
    }

    return parent;
}

// ---------------------------------------------------------------------------
// sparse_lu_column_reach_from_etree
//
// Computes column reachability for column j via elimination tree parent walks.
// For each seed in active_rows with seed < j, walks the parent chain:
//   seed → parent[seed] → parent[parent[seed]] → ...
// until the node reaches >= j or -1 (root).  All visited nodes are in reach.
//
// Result: sorted in increasing order, no duplicates.
//
// Note: this is the etree-based approach (Approach B in design).  It is used
// for standalone symbolic reach tests (Tests 3,4).  The production GP path
// uses the internal sparse_lu_column_reach_from_L (L-structure DFS) instead,
// which handles multi-off-diagonal columns correctly.
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index>
sparse_lu_column_reach_from_etree(
    Index n,
    Index j,
    const std::vector<Index>& parent,
    const std::vector<Index>& active_rows)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_column_reach_from_etree: Index must be signed");
    sparse_lu_detail::sparse_lu_validate_etree_reach_inputs(
        n, j, parent, active_rows);
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<bool>  visited(un, false);
    std::vector<Index> reach;

    for (std::size_t ai = 0u; ai < active_rows.size(); ++ai) {
        Index cur = active_rows[ai];
        // Walk parent chain until node >= j or root (-1)
        while (cur >= Index(0) && cur < j) {
            const std::size_t sc = static_cast<std::size_t>(cur);
            if (visited[sc]) break;  // already included; stop chain walk
            visited[sc] = true;
            reach.push_back(cur);
            cur = parent[sc];
        }
    }

    std::sort(reach.begin(), reach.end());
    return reach;
}

// ===========================================================================
// SLU-7: Symbolic supernode metadata helpers
// ===========================================================================

// ---------------------------------------------------------------------------
// sparse_lu_is_valid_supernode_symbolic
//
// Validates all invariants of a sparse_lu_supernode_symbolic<Index> object
// for a matrix of size n x n.  Returns false if any invariant is violated;
// true if the object is consistent and safe to use.
//
// Invariants checked:
//   sym.valid must be true
//   n >= 0
//   supernode_ptr: size nsup+1, front==0, back==n, strictly increasing
//   column_to_supernode: size n, each entry in [0,nsup), consistent with supernode_ptr
//   parent: size nsup, each entry -1 or (s < parent[s] < nsup)
//   row_ptr: size nsup+1, front==0, monotone, back==row_ind.size()
//   row_ind: per supernode sorted unique, all entries in [0,n)
// ---------------------------------------------------------------------------
template <class Index>
bool sparse_lu_is_valid_supernode_symbolic(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& sym)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_is_valid_supernode_symbolic: Index must be signed");

    if (!sym.valid) return false;
    if (n < Index(0)) return false;

    const std::size_t un = static_cast<std::size_t>(n);

    // supernode_ptr: must have at least 1 element
    if (sym.supernode_ptr.empty()) return false;
    if (sym.supernode_ptr.front() != Index(0)) return false;
    if (sym.supernode_ptr.back() != n) return false;

    const std::size_t nsup = sym.supernode_ptr.size() - 1u;

    // supernode_ptr must be strictly increasing
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (sym.supernode_ptr[s] >= sym.supernode_ptr[s + 1u]) return false;
    }

    // column_to_supernode: size n
    if (sym.column_to_supernode.size() != un) return false;
    for (std::size_t c = 0u; c < un; ++c) {
        const Index s = sym.column_to_supernode[c];
        if (s < Index(0) || static_cast<std::size_t>(s) >= nsup) return false;
        // Consistency: column c must lie within the column range of supernode s
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index b  = sym.supernode_ptr[ss];
        const Index e  = sym.supernode_ptr[ss + 1u];
        const Index ic = static_cast<Index>(c);
        if (ic < b || ic >= e) return false;
    }

    // parent: size nsup; each -1 or forward supernode index
    if (sym.parent.size() != nsup) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index p = sym.parent[s];
        if (p == Index(-1)) continue;
        if (static_cast<Index>(s) >= p) return false;         // must be strictly forward
        if (static_cast<std::size_t>(p) >= nsup) return false; // must be < nsup
    }

    // row_ptr: size nsup+1, front 0, monotone
    if (sym.row_ptr.size() != nsup + 1u) return false;
    if (sym.row_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (sym.row_ptr[s] > sym.row_ptr[s + 1u]) return false;
    }
    // back must match row_ind.size()
    if (sym.row_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(sym.row_ptr[nsup]) != sym.row_ind.size()) return false;

    // row_ind: per supernode sorted unique, all entries in [0, n)
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = sym.row_ptr[s];
        const Index e = sym.row_ptr[s + 1u];
        for (Index p = b; p < e; ++p) {
            const std::size_t sp = static_cast<std::size_t>(p);
            const Index r = sym.row_ind[sp];
            if (r < Index(0) || r >= n) return false;
            if (p + Index(1) < e) {
                if (sym.row_ind[sp] >= sym.row_ind[sp + 1u]) return false;
            }
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// sparse_lu_is_valid_supernode_reach_symbolic
//
// Validates all invariants of a sparse_lu_supernode_reach_symbolic<Index> object
// for a matrix of size n x n with the given supernode partition.
// Returns false if any invariant is violated; true if the object is consistent.
//
// Invariants checked:
//   supernodes.valid and sparse_lu_is_valid_supernode_symbolic pass
//   reach.valid must be true
//   reach_ptr: size nsup+1, front==0, monotone, back==reach_ind.size()
//   reach_ind: per supernode sorted unique, values in (s, nsup) (no self-reach)
//   child_ptr: size nsup+1, front==0, monotone, back==child_ind.size()
//   child_ind: per supernode sorted unique, values in [0, nsup)
//   child/parent consistency: child c in children[s] iff supernodes.parent[c]==s
//   panel_row_ptr: size nsup+1, front==0, monotone, back==panel_row_ind.size()
//   panel_row_ind: per supernode sorted unique, values in [0, n)
// ---------------------------------------------------------------------------
template <class Index>
bool sparse_lu_is_valid_supernode_reach_symbolic(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_reach_symbolic<Index>& reach)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_is_valid_supernode_reach_symbolic: Index must be signed");

    if (!supernodes.valid) return false;
    if (!sparse_lu_is_valid_supernode_symbolic(n, supernodes)) return false;
    if (!reach.valid) return false;
    if (n < Index(0)) return false;

    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    // reach_ptr: size nsup+1, front 0, monotone, back == reach_ind.size()
    if (reach.reach_ptr.size() != nsup + 1u) return false;
    if (reach.reach_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (reach.reach_ptr[s] > reach.reach_ptr[s + 1u]) return false;
    }
    if (reach.reach_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(reach.reach_ptr[nsup]) != reach.reach_ind.size())
        return false;

    // reach_ind: per supernode sorted unique, all entries in (s, nsup) (no self-reach)
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = reach.reach_ptr[s];
        const Index e = reach.reach_ptr[s + 1u];
        Index prev = Index(-1);
        for (Index p = b; p < e; ++p) {
            const Index r = reach.reach_ind[static_cast<std::size_t>(p)];
            if (r < Index(0) || static_cast<std::size_t>(r) >= nsup) return false;
            if (r <= static_cast<Index>(s)) return false; // must be strictly > s
            if (prev >= Index(0) && r <= prev) return false; // strictly increasing
            prev = r;
        }
    }

    // child_ptr: size nsup+1, front 0, monotone, back == child_ind.size()
    if (reach.child_ptr.size() != nsup + 1u) return false;
    if (reach.child_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (reach.child_ptr[s] > reach.child_ptr[s + 1u]) return false;
    }
    if (reach.child_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(reach.child_ptr[nsup]) != reach.child_ind.size())
        return false;

    // child_ind: per supernode sorted unique, values in [0, nsup), parent consistent
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = reach.child_ptr[s];
        const Index e = reach.child_ptr[s + 1u];
        Index prev = Index(-1);
        for (Index p = b; p < e; ++p) {
            const Index c = reach.child_ind[static_cast<std::size_t>(p)];
            if (c < Index(0) || static_cast<std::size_t>(c) >= nsup) return false;
            if (prev >= Index(0) && c <= prev) return false; // strictly increasing
            prev = c;
            // parent consistency: parent of c must equal s
            if (supernodes.parent[static_cast<std::size_t>(c)] != static_cast<Index>(s))
                return false;
        }
    }

    // child/parent completeness: for every c with parent[c]==s, c must be in children[s]
    for (std::size_t c = 0u; c < nsup; ++c) {
        const Index p = supernodes.parent[c];
        if (p == Index(-1)) continue;
        const std::size_t sp = static_cast<std::size_t>(p);
        const Index b = reach.child_ptr[sp];
        const Index e = reach.child_ptr[sp + 1u];
        bool found = false;
        for (Index q = b; q < e; ++q) {
            if (reach.child_ind[static_cast<std::size_t>(q)] == static_cast<Index>(c)) {
                found = true;
                break;
            }
        }
        if (!found) return false;
    }

    // panel_row_ptr: size nsup+1, front 0, monotone, back == panel_row_ind.size()
    if (reach.panel_row_ptr.size() != nsup + 1u) return false;
    if (reach.panel_row_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (reach.panel_row_ptr[s] > reach.panel_row_ptr[s + 1u]) return false;
    }
    if (reach.panel_row_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(reach.panel_row_ptr[nsup]) != reach.panel_row_ind.size())
        return false;

    // panel_row_ind: per supernode sorted unique, all entries in [0, n)
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = reach.panel_row_ptr[s];
        const Index e = reach.panel_row_ptr[s + 1u];
        Index prev = Index(-1);
        for (Index p = b; p < e; ++p) {
            const Index r = reach.panel_row_ind[static_cast<std::size_t>(p)];
            if (r < Index(0) || r >= n) return false;
            if (prev >= Index(0) && r <= prev) return false; // strictly increasing
            prev = r;
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// sparse_lu_build_supernode_symbolic_csc
//
// Constructs symbolic supernode metadata from the CSC sparsity pattern of A
// and a pre-computed A^T A column elimination tree.
//
// Supernode detection rule (conservative structural check on A pattern):
//   Columns j and j+1 may merge iff:
//     1. column_parent[j] == j + 1  (j and j+1 are adjacent in the etree)
//     2. {r in A_col[j] : r > j} \ {j+1} == {r in A_col[j+1] : r > j+1}
//        (the below-diagonal patterns match after accounting for the diagonal shift)
//   This is a conservative approximation; false negatives (singleton fallback) are
//   acceptable.  False positives (incorrect merge) are never produced.
//
// Supernode parent:
//   parent[s] = column_to_supernode[column_parent[last_col_of_s]],
//   or -1 if last_col_of_s has no column parent.
//   Always satisfies: parent[s] == -1  OR  s < parent[s] < nsup.
//
// Row pattern per supernode:
//   Sorted unique union of all structural rows from the CSC columns in the
//   supernode.  Sufficient for symbolic use; does not include fill-in rows.
//
// Returns invalid (valid==false) if:
//   - n < 0 or CSC structure is malformed
//   - column_parent.size() != n
//   - any column_parent entry violates the forward-parent convention
//
// Template parameter Index must be signed (sentinel -1 used for roots/invalid).
// ---------------------------------------------------------------------------
template <class Index>
sparse_lu_supernode_symbolic<Index>
sparse_lu_build_supernode_symbolic_csc(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<Index>& column_parent,
    std::size_t supernode_relaxation = 0u)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_build_supernode_symbolic_csc: Index must be signed");

    sparse_lu_supernode_symbolic<Index> result;

    // Validate CSC pattern; return invalid on any malformed input
    try {
        sparse_lu_detail::sparse_lu_validate_csc_pattern_for_etree(
            n, col_ptr, row_ind);
    } catch (...) {
        result.valid = false;
        return result;
    }

    const std::size_t un = static_cast<std::size_t>(n);

    // Validate column_parent size
    if (column_parent.size() != un) {
        result.valid = false;
        return result;
    }

    // Validate column_parent links: -1 or strictly forward (k < parent[k] < n)
    for (Index k = Index(0); k < n; ++k) {
        const std::size_t sk = static_cast<std::size_t>(k);
        const Index p = column_parent[sk];
        if (!(p == Index(-1) || (k < p && p < n))) {
            result.valid = false;
            return result;
        }
    }

    // n == 0: empty result
    if (n == Index(0)) {
        result.supernode_ptr.push_back(Index(0));
        result.row_ptr.push_back(Index(0));
        result.valid = true;
        return result;
    }

    // Build sorted-unique column patterns from CSC input
    std::vector<std::vector<Index> > col_patterns(un);
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        const Index pb = col_ptr[sj];
        const Index pe = col_ptr[sj + 1u];
        for (Index p = pb; p < pe; ++p) {
            col_patterns[sj].push_back(row_ind[static_cast<std::size_t>(p)]);
        }
        std::sort(col_patterns[sj].begin(), col_patterns[sj].end());
        col_patterns[sj].erase(
            std::unique(col_patterns[sj].begin(), col_patterns[sj].end()),
            col_patterns[sj].end());
    }

    // Supernode detection: same_sn[j] == true means columns j and j+1 belong
    // to the same supernode.
    const std::size_t n_pairs = static_cast<std::size_t>(n) - 1u;
    std::vector<bool> same_sn(n_pairs, false);

    for (Index j = Index(0); j + Index(1) < n; ++j) {
        const std::size_t sj  = static_cast<std::size_t>(j);
        const std::size_t sj1 = sj + 1u;

        // Condition 1: j+1 must be the etree parent of j
        if (column_parent[sj] != j + Index(1)) continue;

        // Condition 2: rows strictly below j in col j, excluding j+1,
        //              must equal rows strictly below j+1 in col j+1.
        std::vector<Index> below_j_trimmed;
        {
            const std::vector<Index>& cp = col_patterns[sj];
            for (std::size_t bi = 0u; bi < cp.size(); ++bi) {
                const Index r = cp[bi];
                if (r > j && r != j + Index(1)) below_j_trimmed.push_back(r);
            }
        }

        std::vector<Index> below_j1;
        {
            const std::vector<Index>& cp = col_patterns[sj1];
            for (std::size_t bi = 0u; bi < cp.size(); ++bi) {
                const Index r = cp[bi];
                if (r > j + Index(1)) below_j1.push_back(r);
            }
        }

        if (below_j_trimmed == below_j1) {
            same_sn[sj] = true;
        }
    }

    // Build fundamental supernode column boundaries first.
    std::vector<Index> fundamental_sn_ptr;
    fundamental_sn_ptr.push_back(Index(0));
    for (Index j = Index(0); j + Index(1) < n; ++j) {
        if (!same_sn[static_cast<std::size_t>(j)]) {
            fundamental_sn_ptr.push_back(j + Index(1));
        }
    }
    fundamental_sn_ptr.push_back(n);

    // SLU-14R.5 relaxed amalgamation:
    // Greedily merge adjacent fundamental supernodes when the merge adds no more
    // than supernode_relaxation panel slots per column relative to each
    // fundamental panel envelope.  supernode_relaxation == 0 preserves the
    // fundamental partition exactly.
    std::vector<Index> sn_ptr = fundamental_sn_ptr;
    if (supernode_relaxation > 0u && fundamental_sn_ptr.size() > 2u) {
        const std::size_t nf = fundamental_sn_ptr.size() - 1u;

        std::vector<std::vector<Index> > fund_rows(nf);
        for (std::size_t f = 0u; f < nf; ++f) {
            const Index b = fundamental_sn_ptr[f];
            const Index e = fundamental_sn_ptr[f + 1u];
            std::vector<Index> rows;
            for (Index c = b; c < e; ++c) {
                rows.push_back(c);
                const std::vector<Index>& cp =
                    col_patterns[static_cast<std::size_t>(c)];
                rows.insert(rows.end(), cp.begin(), cp.end());
            }
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
            fund_rows[f].swap(rows);
        }

        std::vector<Index> relaxed_sn_ptr;
        relaxed_sn_ptr.push_back(Index(0));

        std::size_t f = 0u;
        while (f < nf) {
            std::vector<Index> current_rows = fund_rows[f];
            std::size_t next = f + 1u;

            while (next < nf) {
                std::vector<Index> candidate_rows = current_rows;
                candidate_rows.insert(candidate_rows.end(),
                                      fund_rows[next].begin(),
                                      fund_rows[next].end());
                std::sort(candidate_rows.begin(), candidate_rows.end());
                candidate_rows.erase(
                    std::unique(candidate_rows.begin(), candidate_rows.end()),
                    candidate_rows.end());

                bool merge_ok = true;
                for (std::size_t g = f; g <= next && merge_ok; ++g) {
                    std::size_t added = 0u;
                    std::size_t ib = 0u;
                    for (std::size_t ia = 0u; ia < candidate_rows.size(); ++ia) {
                        const Index r = candidate_rows[ia];
                        while (ib < fund_rows[g].size() && fund_rows[g][ib] < r) {
                            ++ib;
                        }
                        if (ib >= fund_rows[g].size() || fund_rows[g][ib] != r) {
                            ++added;
                            if (added > supernode_relaxation) {
                                merge_ok = false;
                                break;
                            }
                        }
                    }
                }

                if (!merge_ok) break;
                current_rows.swap(candidate_rows);
                ++next;
            }

            relaxed_sn_ptr.push_back(fundamental_sn_ptr[next]);
            f = next;
        }

        sn_ptr.swap(relaxed_sn_ptr);
    }

    const Index nsup = static_cast<Index>(sn_ptr.size()) - Index(1);

    // Build column_to_supernode
    std::vector<Index> col_to_sn(un, Index(-1));
    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index b = sn_ptr[ss];
        const Index e = sn_ptr[ss + 1u];
        for (Index c = b; c < e; ++c) {
            col_to_sn[static_cast<std::size_t>(c)] = s;
        }
    }

    // Build supernode parent array
    // parent[s] = supernode of column_parent[last_col_of_s], or -1 if root.
    // The merge rule guarantees column_parent[last_col] is outside supernode s,
    // so ps != s always for valid inputs; the while-loop below is defensive.
    std::vector<Index> sn_parent(static_cast<std::size_t>(nsup), Index(-1));
    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index last_col = sn_ptr[ss + 1u] - Index(1);
        Index pcol = column_parent[static_cast<std::size_t>(last_col)];
        if (pcol == Index(-1)) {
            sn_parent[ss] = Index(-1);
        } else {
            Index ps = col_to_sn[static_cast<std::size_t>(pcol)];
            // Defensive: walk column parent until we reach a different supernode
            while (ps == s) {
                pcol = column_parent[static_cast<std::size_t>(pcol)];
                if (pcol == Index(-1)) { ps = Index(-1); break; }
                ps = col_to_sn[static_cast<std::size_t>(pcol)];
            }
            sn_parent[ss] = ps;
        }
    }

    // Build structural row patterns per supernode:
    // sorted unique union of all A-column rows within the supernode.
    std::vector<Index> sn_row_ptr;
    std::vector<Index> sn_row_ind;
    sn_row_ptr.push_back(Index(0));

    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index b = sn_ptr[ss];
        const Index e = sn_ptr[ss + 1u];

        std::vector<Index> rows;
        for (Index c = b; c < e; ++c) {
            const std::vector<Index>& cp =
                col_patterns[static_cast<std::size_t>(c)];
            rows.insert(rows.end(), cp.begin(), cp.end());
        }
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

        for (std::size_t ri = 0u; ri < rows.size(); ++ri) {
            sn_row_ind.push_back(rows[ri]);
        }
        sn_row_ptr.push_back(static_cast<Index>(sn_row_ind.size()));
    }

    result.supernode_ptr       = sn_ptr;
    result.column_to_supernode = col_to_sn;
    result.parent              = sn_parent;
    result.row_ptr             = sn_row_ptr;
    result.row_ind             = sn_row_ind;
    result.valid               = true;

    return result;
}

// ---------------------------------------------------------------------------
// sparse_lu_build_supernode_symbolic_from_l   [SLU-SNA1 P1-B, D-2/D-4]
//
// Constructs supernode symbolic metadata for the SUPERNODAL path from the
// ACTUAL GP L pattern (baseline CSC L, pivoted row space), replacing the
// A-pattern approximation of sparse_lu_build_supernode_symbolic_csc as the
// partition source of the supernodal pipeline (D-3: baseline paths keep the
// A-pattern builder).
//
// Motivation (SLU-SNA1 Phase 0): the A-pattern merge rule degenerates via
// trivially-satisfied empty-set comparisons (below sets both empty), producing
// width-n (SNQ1) or width-0.065n (SNQ3) supernodes whose w x (w + |L union|)
// panel envelopes explode.  The L-pattern fundamental rule merges j and j+1 iff
//   1. j+1 == min below(j)                (L parent adjacency; below(j) nonempty)
//   2. below(j) \ {j+1} == below(j+1)     (exact SET equality, not counts)
// where below(j) = { r : L(r, j) != 0 } (strictly lower, pivoted rows).
// Condition 1 requires a nonempty below(j), so the empty-set degeneracy is
// structurally impossible.  Condition 2 gives exact nesting within a supernode:
//   below(first) ⊇ below(first+1) ⊇ ... ⊇ below(last),
// hence the per-supernode L-row union consumed by the storage bootstrap
// (build_supernode_numeric_from_csc -> l_row_ind) equals below(first_col)
// exactly: the panel is the true fill envelope (D-4) with zero slack.
//
// Output mirrors sparse_lu_build_supernode_symbolic_csc:
//   supernode_ptr / column_to_supernode : partition (fundamental + relaxation)
//   parent[s]  = supernode of min below(last_col_of_s), or -1 (root)
//   row_ind[s] = sorted unique ORIGINAL row ids (row_perm[new] = old) of
//                { below(c) : c in supernode } ∪ diagonal block rows.
//     Mapping to original ids keeps the retained-supernode_info convention
//     consumed by gate5 accountability (which maps back by inv_row_perm);
//     coverage of the bootstrap panel rows holds with equality (design §7).
//   supernode_relaxation: same greedy amalgamation semantics as the A-pattern
//     builder (added panel slots per column <= relaxation), evaluated on
//     L-pattern rows.  Default 0 keeps the fundamental partition exactly.
//
// Returns invalid (valid==false) if n < 0, CSC L is malformed, or row_perm
// size mismatches; callers fall back to the A-pattern info in that case.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_supernode_symbolic<Index>
sparse_lu_build_supernode_symbolic_from_l(
    Index n,
    const baseline_lu_storage<T, Index>& csc_lu,
    std::size_t supernode_relaxation = 0u)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_build_supernode_symbolic_from_l: Index must be signed");

    sparse_lu_supernode_symbolic<Index> result;
    result.valid = false;

    if (n < Index(0)) return result;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.L, n, n)) return result;

    const std::size_t un = static_cast<std::size_t>(n);
    if (csc_lu.row_perm.size() != un) return result;

    if (n == Index(0)) {
        result.supernode_ptr.push_back(Index(0));
        result.row_ptr.push_back(Index(0));
        result.valid = true;
        return result;
    }

    // ---- Sorted-unique per-column below sets from actual CSC L ----
    // below(j) = { r : L(r, j) != 0 }; strictly lower (r > j) by the GP
    // emission invariant (kept defensively by the r > j filter).
    std::vector<Index> below_ptr(un + 1u, Index(0));
    std::vector<Index> below_rows;
    below_rows.reserve(csc_lu.L.row_ind.size());
    {
        std::vector<Index> tmp;
        for (Index j = Index(0); j < n; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            tmp.clear();
            for (Index k = csc_lu.L.col_ptr[sj];
                 k < csc_lu.L.col_ptr[sj + 1u]; ++k) {
                const Index r = csc_lu.L.row_ind[static_cast<std::size_t>(k)];
                if (r > j) tmp.push_back(r);
            }
            std::sort(tmp.begin(), tmp.end());
            tmp.erase(std::unique(tmp.begin(), tmp.end()), tmp.end());
            below_ptr[sj] = static_cast<Index>(below_rows.size());
            below_rows.insert(below_rows.end(), tmp.begin(), tmp.end());
        }
        below_ptr[un] = static_cast<Index>(below_rows.size());
    }

    // ---- Fundamental merge rule on actual L below sets ----
    const std::size_t n_pairs = un - 1u;
    std::vector<bool> same_sn(n_pairs, false);
    for (Index j = Index(0); j + Index(1) < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        const Index pb = below_ptr[sj];
        const Index pe = below_ptr[sj + 1u];
        // Condition 1: below(j) nonempty and its minimum is j+1.
        if (pb == pe) continue;
        if (below_rows[static_cast<std::size_t>(pb)] != j + Index(1)) continue;
        // Condition 2: below(j) \ {j+1} == below(j+1) as sorted sets.
        const Index qb = below_ptr[sj + 1u];
        const Index qe = below_ptr[sj + 2u];
        if ((pe - pb) - Index(1) != (qe - qb)) continue;
        bool equal = true;
        for (Index t = Index(0); t < qe - qb; ++t) {
            if (below_rows[static_cast<std::size_t>(pb + Index(1) + t)] !=
                below_rows[static_cast<std::size_t>(qb + t)]) {
                equal = false;
                break;
            }
        }
        if (equal) same_sn[sj] = true;
    }

    // ---- Fundamental supernode column boundaries ----
    std::vector<Index> fundamental_sn_ptr;
    fundamental_sn_ptr.push_back(Index(0));
    for (Index j = Index(0); j + Index(1) < n; ++j) {
        if (!same_sn[static_cast<std::size_t>(j)]) {
            fundamental_sn_ptr.push_back(j + Index(1));
        }
    }
    fundamental_sn_ptr.push_back(n);

    // ---- Relaxed amalgamation (same semantics as the A-pattern builder):
    // greedily merge adjacent fundamental supernodes when the merge adds no
    // more than supernode_relaxation panel slots per column relative to each
    // fundamental panel footprint ({c} ∪ below(c) per column, pivoted rows).
    std::vector<Index> sn_ptr = fundamental_sn_ptr;
    if (supernode_relaxation > 0u && fundamental_sn_ptr.size() > 2u) {
        const std::size_t nf = fundamental_sn_ptr.size() - 1u;

        std::vector<std::vector<Index> > fund_rows(nf);
        for (std::size_t f = 0u; f < nf; ++f) {
            const Index b = fundamental_sn_ptr[f];
            const Index e = fundamental_sn_ptr[f + 1u];
            std::vector<Index> rows;
            for (Index c = b; c < e; ++c) {
                rows.push_back(c);
                const std::size_t sc = static_cast<std::size_t>(c);
                for (Index k = below_ptr[sc]; k < below_ptr[sc + 1u]; ++k) {
                    rows.push_back(below_rows[static_cast<std::size_t>(k)]);
                }
            }
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
            fund_rows[f].swap(rows);
        }

        std::vector<Index> relaxed_sn_ptr;
        relaxed_sn_ptr.push_back(Index(0));

        std::size_t f = 0u;
        while (f < nf) {
            std::vector<Index> current_rows = fund_rows[f];
            std::size_t next = f + 1u;

            while (next < nf) {
                std::vector<Index> candidate_rows = current_rows;
                candidate_rows.insert(candidate_rows.end(),
                                      fund_rows[next].begin(),
                                      fund_rows[next].end());
                std::sort(candidate_rows.begin(), candidate_rows.end());
                candidate_rows.erase(
                    std::unique(candidate_rows.begin(), candidate_rows.end()),
                    candidate_rows.end());

                bool merge_ok = true;
                for (std::size_t g = f; g <= next && merge_ok; ++g) {
                    std::size_t added = 0u;
                    std::size_t ib = 0u;
                    for (std::size_t ia = 0u; ia < candidate_rows.size(); ++ia) {
                        const Index r = candidate_rows[ia];
                        while (ib < fund_rows[g].size() && fund_rows[g][ib] < r) {
                            ++ib;
                        }
                        if (ib >= fund_rows[g].size() || fund_rows[g][ib] != r) {
                            ++added;
                            if (added > supernode_relaxation) {
                                merge_ok = false;
                                break;
                            }
                        }
                    }
                }

                if (!merge_ok) break;
                current_rows.swap(candidate_rows);
                ++next;
            }

            relaxed_sn_ptr.push_back(fundamental_sn_ptr[next]);
            f = next;
        }

        sn_ptr.swap(relaxed_sn_ptr);
    }

    const Index nsup = static_cast<Index>(sn_ptr.size()) - Index(1);

    // ---- column_to_supernode ----
    std::vector<Index> col_to_sn(un, Index(-1));
    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index b = sn_ptr[ss];
        const Index e = sn_ptr[ss + 1u];
        for (Index c = b; c < e; ++c) {
            col_to_sn[static_cast<std::size_t>(c)] = s;
        }
    }

    // ---- Supernode parent: supernode of min below(last_col), or -1.
    // min below(last_col) > last_col >= col_end - 1, so the parent column is
    // always outside supernode s (forward-parent invariant holds); the walk
    // below is defensive, mirroring the A-pattern builder.
    std::vector<Index> sn_parent(static_cast<std::size_t>(nsup), Index(-1));
    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        Index pcol;
        {
            const Index last_col = sn_ptr[ss + 1u] - Index(1);
            const std::size_t sl = static_cast<std::size_t>(last_col);
            pcol = (below_ptr[sl] == below_ptr[sl + 1u])
                       ? Index(-1)
                       : below_rows[static_cast<std::size_t>(below_ptr[sl])];
        }
        if (pcol == Index(-1)) {
            sn_parent[ss] = Index(-1);
        } else {
            Index ps = col_to_sn[static_cast<std::size_t>(pcol)];
            while (ps == s) {
                const std::size_t sp = static_cast<std::size_t>(pcol);
                pcol = (below_ptr[sp] == below_ptr[sp + 1u])
                           ? Index(-1)
                           : below_rows[static_cast<std::size_t>(below_ptr[sp])];
                if (pcol == Index(-1)) { ps = Index(-1); break; }
                ps = col_to_sn[static_cast<std::size_t>(pcol)];
            }
            sn_parent[ss] = ps;
        }
    }

    // ---- Row pattern per supernode: union of below(c) over the supernode's
    // columns plus the diagonal block rows, mapped to ORIGINAL row ids via
    // row_perm[new] = old, sorted unique.  For the fundamental partition the
    // union equals below(first_col) by nesting; the union form also stays
    // correct under relaxation-merged supernodes.
    std::vector<Index> sn_row_ptr;
    std::vector<Index> sn_row_ind;
    sn_row_ptr.push_back(Index(0));

    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index b = sn_ptr[ss];
        const Index e = sn_ptr[ss + 1u];

        std::vector<Index> rows;
        for (Index c = b; c < e; ++c) {
            rows.push_back(csc_lu.row_perm[static_cast<std::size_t>(c)]);
            const std::size_t sc = static_cast<std::size_t>(c);
            for (Index k = below_ptr[sc]; k < below_ptr[sc + 1u]; ++k) {
                rows.push_back(csc_lu.row_perm[static_cast<std::size_t>(
                    below_rows[static_cast<std::size_t>(k)])]);
            }
        }
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

        for (std::size_t ri = 0u; ri < rows.size(); ++ri) {
            sn_row_ind.push_back(rows[ri]);
        }
        sn_row_ptr.push_back(static_cast<Index>(sn_row_ind.size()));
    }

    result.supernode_ptr       = sn_ptr;
    result.column_to_supernode = col_to_sn;
    result.parent              = sn_parent;
    result.row_ptr             = sn_row_ptr;
    result.row_ind             = sn_row_ind;
    result.valid               = true;

    return result;
}

// ===========================================================================
// SLU-8: Symbolic supernode reach metadata builder
// ===========================================================================

// ---------------------------------------------------------------------------
// sparse_lu_for_each_supernode_ancestor
//
// Visits the ancestors of supernode s in the supernode etree, in strictly
// increasing order (forward-parent invariant: parent[t] > t), by walking
// parent[] from parent[s] to the root.  O(depth(s)) time, O(1) space.
// This is the lazy replacement for the materialized ancestor reach
// (sparse_lu_supernode_reach_symbolic): the reach of s is exactly the
// visited sequence.  Visitor: void(Index ancestor).
//
// Preconditions (caller responsibility, matching the builder's validation):
//   parent.size() == nsup; parent[t] == -1 or t < parent[t] < nsup.
// ---------------------------------------------------------------------------
template <class Index, class Visitor>
void sparse_lu_for_each_supernode_ancestor(
    const std::vector<Index>& parent, const Index s, Visitor visit)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_for_each_supernode_ancestor: Index must be signed");
    Index t = parent[static_cast<std::size_t>(s)];
    while (t != Index(-1)) {
        visit(t);
        t = parent[static_cast<std::size_t>(t)];
    }
}

// ---------------------------------------------------------------------------
// sparse_lu_build_supernode_reach_symbolic
//
// Constructs symbolic supernode reach metadata from validated supernode partition.
//
// Reach definition (Approach B: supernode etree ancestor reach):
//   For each supernode s:
//     Walk the parent chain from parent[s] upward to the root.
//     All ancestors encountered are collected as reach[s].
//   This is conservative (over-approximation) but safe:
//     - no false negatives (all true structural dependencies are included)
//     - may include false positives (conservative ancestor reach)
//     - entries are strictly > s (no self-reach), sorted unique
//
// Child metadata:
//   For each supernode c with parent[c] == p (p != -1):
//     c is added to the children of p.
//   Children are sorted unique.
//
// Panel row pattern (SLU-8 symbolic panel row pattern):
//   panel_row_pattern[s] = SLU-7 structural row pattern for supernode s.
//   This is a structural over-approximation for future supernodal symbolic
//   planning. It is NOT the final numeric row structure under pivoting.
//
// Returns invalid (valid==false) if:
//   - supernodes.valid == false
//   - sparse_lu_is_valid_supernode_symbolic(n, supernodes) == false
//   - n < 0
//
// Template parameter Index must be signed (sentinel -1 used for roots).
//
// ON-DEMAND DIAGNOSTIC (SLU-RQ1): this builder materializes the full
// ancestor reach, which is O(sum of ancestor-chain lengths) — O(nsup^2)
// time and memory on chain etrees.  It is NOT called on any production
// path; production code derives reach lazily via
// sparse_lu_for_each_supernode_ancestor.  Intended for tests and
// small-fixture diagnostics only.
// ---------------------------------------------------------------------------
template <class Index>
sparse_lu_supernode_reach_symbolic<Index>
sparse_lu_build_supernode_reach_symbolic(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_build_supernode_reach_symbolic: Index must be signed");

    sparse_lu_supernode_reach_symbolic<Index> result;

    // Validate supernode metadata before any access
    if (!supernodes.valid) { result.valid = false; return result; }
    if (n < Index(0))      { result.valid = false; return result; }
    if (!sparse_lu_is_valid_supernode_symbolic(n, supernodes)) {
        result.valid = false;
        return result;
    }

    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    // Empty case (nsup==0, n==0)
    if (nsup == 0u) {
        result.reach_ptr.push_back(Index(0));
        result.child_ptr.push_back(Index(0));
        result.panel_row_ptr.push_back(Index(0));
        result.valid = true;
        return result;
    }

    // Build ancestor reach per supernode (Approach B: walk parent chain to root).
    // The forward-parent invariant guarantees parent[s] > s, so the chain is
    // strictly increasing -- no cycle is possible. Walk until -1 (root).
    result.reach_ptr.resize(nsup + 1u, Index(0));
    for (std::size_t s = 0u; s < nsup; ++s) {
        sparse_lu_for_each_supernode_ancestor<Index>(
            supernodes.parent, static_cast<Index>(s),
            [&result](const Index t) { result.reach_ind.push_back(t); });
        // Parent chain is already strictly increasing by forward-parent invariant.
        // Sort and unique for defensive determinism (handles any edge cases).
        const std::size_t b = static_cast<std::size_t>(result.reach_ptr[s]);
        const std::size_t e = result.reach_ind.size();
        std::sort(result.reach_ind.begin() + static_cast<std::ptrdiff_t>(b),
                  result.reach_ind.end());
        result.reach_ind.erase(
            std::unique(result.reach_ind.begin() + static_cast<std::ptrdiff_t>(b),
                        result.reach_ind.end()),
            result.reach_ind.end());
        (void)e; // suppress unused-variable warning
        result.reach_ptr[s + 1u] = static_cast<Index>(result.reach_ind.size());
    }

    // Build children per supernode: for each c with parent[c] != -1, add c to children[parent[c]].
    std::vector<std::vector<Index> > child_sets(nsup);
    for (std::size_t c = 0u; c < nsup; ++c) {
        const Index p = supernodes.parent[c];
        if (p != Index(-1)) {
            child_sets[static_cast<std::size_t>(p)].push_back(static_cast<Index>(c));
        }
    }
    result.child_ptr.resize(nsup + 1u, Index(0));
    for (std::size_t s = 0u; s < nsup; ++s) {
        // child indices are always < s (forward-parent invariant), already sorted ascending
        std::sort(child_sets[s].begin(), child_sets[s].end());
        for (std::size_t ci = 0u; ci < child_sets[s].size(); ++ci) {
            result.child_ind.push_back(child_sets[s][ci]);
        }
        result.child_ptr[s + 1u] = static_cast<Index>(result.child_ind.size());
    }

    // Build symbolic panel row patterns from SLU-7 structural row patterns.
    // SLU-8 symbolic panel row pattern.
    // This metadata is a structural over-approximation used for future
    // supernodal symbolic planning. It is not the final numeric row structure
    // after threshold partial pivoting.
    result.panel_row_ptr.resize(nsup + 1u, Index(0));
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = supernodes.row_ptr[s];
        const Index e = supernodes.row_ptr[s + 1u];
        for (Index p = b; p < e; ++p) {
            result.panel_row_ind.push_back(
                supernodes.row_ind[static_cast<std::size_t>(p)]);
        }
        result.panel_row_ptr[s + 1u] = static_cast<Index>(result.panel_row_ind.size());
        // SLU-7 row_ind per supernode is already sorted unique; no additional sort needed.
    }

    result.valid = true;
    return result;
}

// ===========================================================================
// SLU-9: Supernode-aware triangular solve helper
// Included here after SLU-7/8 declarations so sparse_lu_is_valid_supernode_symbolic
// and sparse_lu_supernode_symbolic are in scope.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_supernode_solve_impl.hpp>

// ===========================================================================
// SLU-10.1: Numeric supernode metadata validator
//
// sparse_lu_is_valid_supernode_numeric: validates all invariants of a
// sparse_lu_supernode_numeric<T, Index> object built for matrix size n and
// given symbolic supernode partition.
//
// Checks (numeric.valid must be true first):
//   n >= 0 and supernodes.valid and sparse_lu_is_valid_supernode_symbolic pass
//   supernode_ptr, column_to_supernode match symbolic partition
//   l_row_ptr/u_row_ptr: size nsup+1, front==0, monotone, back==*_row_ind.size()
//   l_row_ind/u_row_ind: per supernode sorted unique, entries in [0,n)
//   diag_block_ptr: size nsup+1, front==0, monotone, back==diag_block_values.size()
//   each diag block size == width^2 (width = supernode column count)
//
// Placed before the SLU-10 impl include so build_supernode_numeric_from_csc
// can use it as a postcondition check.
// ===========================================================================
template <class T, class Index>
bool sparse_lu_is_valid_supernode_numeric(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& numeric)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_is_valid_supernode_numeric: Index must be signed");

    if (!numeric.valid) return false;
    if (!supernodes.valid) return false;
    if (!sparse_lu_is_valid_supernode_symbolic(n, supernodes)) return false;
    if (n < Index(0)) return false;

    const std::size_t un   = static_cast<std::size_t>(n);
    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    // Partition must match symbolic
    if (numeric.supernode_ptr.size() != nsup + 1u) return false;
    for (std::size_t i = 0u; i <= nsup; ++i) {
        if (numeric.supernode_ptr[i] != supernodes.supernode_ptr[i]) return false;
    }
    if (numeric.column_to_supernode.size() != un) return false;
    for (std::size_t i = 0u; i < un; ++i) {
        if (numeric.column_to_supernode[i] != supernodes.column_to_supernode[i])
            return false;
    }

    // l_row_ptr: size nsup+1, front==0, monotone, back==l_row_ind.size()
    if (numeric.l_row_ptr.size() != nsup + 1u) return false;
    if (numeric.l_row_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (numeric.l_row_ptr[s] > numeric.l_row_ptr[s + 1u]) return false;
    }
    if (numeric.l_row_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(numeric.l_row_ptr[nsup]) != numeric.l_row_ind.size())
        return false;

    // l_row_ind: per supernode sorted unique, entries in [0,n)
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = numeric.l_row_ptr[s];
        const Index e = numeric.l_row_ptr[s + 1u];
        for (Index p = b; p < e; ++p) {
            const std::size_t sp = static_cast<std::size_t>(p);
            const Index r = numeric.l_row_ind[sp];
            if (r < Index(0) || r >= n) return false;
            if (p + Index(1) < e) {
                if (numeric.l_row_ind[sp] >= numeric.l_row_ind[sp + 1u]) return false;
            }
        }
    }

    // u_row_ptr: size nsup+1, front==0, monotone, back==u_row_ind.size()
    if (numeric.u_row_ptr.size() != nsup + 1u) return false;
    if (numeric.u_row_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (numeric.u_row_ptr[s] > numeric.u_row_ptr[s + 1u]) return false;
    }
    if (numeric.u_row_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(numeric.u_row_ptr[nsup]) != numeric.u_row_ind.size())
        return false;

    // u_row_ind: per supernode sorted unique, entries in [0,n)
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index b = numeric.u_row_ptr[s];
        const Index e = numeric.u_row_ptr[s + 1u];
        for (Index p = b; p < e; ++p) {
            const std::size_t sp = static_cast<std::size_t>(p);
            const Index r = numeric.u_row_ind[sp];
            if (r < Index(0) || r >= n) return false;
            if (p + Index(1) < e) {
                if (numeric.u_row_ind[sp] >= numeric.u_row_ind[sp + 1u]) return false;
            }
        }
    }

    // diag_block_ptr: size nsup+1, front==0, monotone, back==diag_block_values.size()
    if (numeric.diag_block_ptr.size() != nsup + 1u) return false;
    if (numeric.diag_block_ptr.front() != Index(0)) return false;
    for (std::size_t s = 0u; s < nsup; ++s) {
        if (numeric.diag_block_ptr[s] > numeric.diag_block_ptr[s + 1u]) return false;
    }
    if (numeric.diag_block_ptr[nsup] < Index(0)) return false;
    if (static_cast<std::size_t>(numeric.diag_block_ptr[nsup]) !=
            numeric.diag_block_values.size()) return false;

    // each diagonal block size must equal width^2
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index width = supernodes.supernode_ptr[s + 1u] - supernodes.supernode_ptr[s];
        const Index blk   = numeric.diag_block_ptr[s + 1u] - numeric.diag_block_ptr[s];
        if (blk != width * width) return false;
    }

    return true;
}

// ===========================================================================
// SLU-10: Supernodal prototype numeric metadata builder
// Included here after sparse_lu_is_valid_supernode_numeric is defined so the
// builder can use it as a postcondition check before setting valid = true.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_supernode_numeric_impl.hpp>

// ===========================================================================
// SLU-8R.2: Supernodal storage bootstrap from CSC L/U
// Included after tsparse_sparse_lu_supernode_numeric_impl.hpp so that
// sparse_lu_supernode_numeric<T, Index> is in scope.
// Provides bootstrap_supernodal_storage_from_csc in namespace sparse_lu_detail.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_supernodal_storage_bootstrap_impl.hpp>

// ===========================================================================
// SLU-8R.3: Supernode-panel left-looking update -- §17.2(A) implementation.
// Included after tsparse_sparse_lu_supernodal_storage_bootstrap_impl.hpp so
// that supernodal_lu_storage<T, Index> and sparse_lu_dense_kernel<T> are
// fully defined and implemented (dense_kernel_impl.hpp included at line ~788).
// Provides the shared panel-update machinery in namespace sparse_lu_detail:
// compute_panel_update_set (3-arg), panel workspace gather/apply/scatter,
// sparse_lu_build_col_to_supernode_, and the SN-OPT helpers.  These are called
// by the production true-numeric driver (tsparse_sparse_lu_true_numeric_impl.hpp).
// (The transitional §17.2(A) driver run_supernode_panel_leftlooking_update was
// removed by SLU-CLN1, 2026-07-05.)
//
// SCOPE BOUNDARY:
//   §17.2(A): IMPLEMENTED (supernode-panel update, trsm + gemm/gemv via adapter).
//   §17.2(B): implemented in tsparse_sparse_lu_within_panel_factor_impl.hpp (below).
//   true_supernodal_numeric: always false (historical conformance marker).
//   Gate 6: resolved by A_eff-origin factorization in SLU-8R.5.5.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_supernode_panel_update_impl.hpp>

// ===========================================================================
// SLU-8R.4: §17.2(B) Within-panel factorization -- explicit pivot search.
// Included after tsparse_sparse_lu_supernode_panel_update_impl.hpp so that
// supernodal_lu_storage<T, Index>, sparse_lu_dense_kernel<T>, and
// sparse_lu_scalar_policy<T> are all in scope.
// Provides factorize_within_panel_single and the shared within-panel parts
// (workspace, pivot search, row swap, L-multiplier scale, rank-1 update,
// scatter) in namespace sparse_lu_detail, called per-panel by the production
// true-numeric driver (tsparse_sparse_lu_true_numeric_impl.hpp).
// (The transitional §17.2(B) driver run_within_panel_factorization was
// removed by SLU-CLN1, 2026-07-05.)
//
// SCOPE BOUNDARY:
//   §17.2(B): IMPLEMENTED (pivot search over active panel height, threshold pivoting,
//             row swap, L-multiplier scale, ger rank-1 update via adapter).
//   §18.2: storage-native solve is implemented (tsparse_sparse_lu_supernodal_solve_impl.hpp).
//     Activated after A_eff-origin factorization (SLU-8R.5.5) sets true_numeric_source.
//   true_supernodal_numeric: always false (historical conformance marker; genuine status
//     tracked by supernodal_true_numeric_success).
//   within_panel_factorization_is_numeric_source: always false (historical sub-phase marker).
//   within_panel_used_getrf: always false (§25 prohibition on getrf as pivot search).
//   Gate 6: resolved by A_eff-origin factorization in SLU-8R.5.5.
//   issue_SLU8_contract_violation.md: OPEN.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_within_panel_factor_impl.hpp>
// SLU-8R.5: §18.2 storage-native supernodal solve implementation.
// Provides: can_use_supernodal_storage_solve, supernodal_l_solve_single_rhs,
// supernodal_u_solve_single_rhs, supernodal_storage_solve_single_rhs,
// sparse_lu_make_supernodal_factor_for_testing, sparse_lu_supernodal_storage_solve.
#include <vcp/tsparse/detail/tsparse_sparse_lu_supernodal_solve_impl.hpp>
// SLU-8R.5.5: A_eff-origin true numeric supernodal factorization.
// Provides: sparse_lu_factorize_supernodal_from_a_eff,
//           sparse_lu_initialize_supernodal_from_a_eff,
//           sparse_lu_compute_supernodal_factorization_residual.
#include <vcp/tsparse/detail/tsparse_sparse_lu_true_numeric_impl.hpp>
// SLU-MF: multifrontal numeric source (production) + numeric-source dispatcher.
// Included after the left-looking driver so it can reuse a_eff_lookup and
// compute_supernodal_lu_residual_sparse, and dispatch to the diagnostic driver.
// Provides: sparse_lu_factorize_supernodal_multifrontal,
//           sparse_lu_factorize_supernodal_numeric_source.
#include <vcp/tsparse/detail/tsparse_sparse_lu_multifrontal_impl.hpp>

// ---------------------------------------------------------------------------
// SLU-10.1: Public wrapper for sparse_lu_build_supernode_numeric_from_csc.
// Exposes the detail builder for direct testing of malformed CSC input.
// Parameter order: (n, sym, csc_lu) -- matches test expectations.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_supernode_numeric<T, Index>
sparse_lu_build_supernode_numeric_from_csc(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& sym,
    const baseline_lu_storage<T, Index>& csc_lu)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_build_supernode_numeric_from_csc: Index must be signed");
    return sparse_lu_detail::build_supernode_numeric_from_csc(n, csc_lu, sym);
}

// ---------------------------------------------------------------------------
// sparse_lu_solve_supernode_aware
//
// SLU-9/SLU-13: Free function to solve A*x = b using stored L/U CSC factors,
// traversing columns in supernode-grouped order.
//
// SLU-13 contract:
//   if fac.uses_supernodal_prototype() == true:
//     Validates supernode numeric metadata and CSC consistency (same precondition
//     checks as fac.solve(b)), then dispatches to supernode-aware CSC solve path.
//     Throws state_error if metadata is invalid or inconsistent with CSC L/U.
//   if fac.uses_supernodal_prototype() == false:
//     Uses supernode-aware CSC solve with singleton fallback (mathematically
//     equivalent to fac.solve(b) for baseline/auto factors).
//   In all cases:
//     Validates factor success and RHS size.
//     Does not perform unchecked metadata access.
//     Does not silently produce wrong results.
//
// Throws:
//   - vcp::state_error if factor is not valid or metadata invariant violated
//   - vcp::invalid_argument if b.size() != n or factor storage is malformed
// ---------------------------------------------------------------------------
template <class T, class Index>
std::vector<T>
sparse_lu_solve_supernode_aware(
    const sparse_lu_factorization<T, Index>& fac,
    const std::vector<T>& b)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_solve_supernode_aware: Index must be signed");

    if (!fac.valid()) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_solve_supernode_aware: factor is not valid (status=",
            sparse_lu_status_to_string(fac.info_.status), ")");
    }

    // SLU-13: explicit supernodal prototype path -- mirror fac.solve(b) preconditions
    // to avoid unchecked metadata access.
    if (fac.uses_supernodal_prototype_) {
        if (!fac.supernode_numeric_info_.valid) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_solve_supernode_aware: "
                "supernode numeric metadata is not valid (uses_supernodal_prototype=true)");
        }
        if (!sparse_lu_is_valid_supernode_numeric(
                fac.info_.n, fac.supernode_info_, fac.supernode_numeric_info_)) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_solve_supernode_aware: "
                "supernode numeric metadata fails invariant check");
        }
        if (!sparse_lu_verify_supernode_numeric_against_csc(fac)) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_solve_supernode_aware: "
                "supernode numeric metadata is inconsistent with CSC L/U factors");
        }
        return sparse_lu_detail::solve_baseline_storage_supernode_aware(
            fac.baseline_, fac.info_.n, fac.supernode_info_, b);
    }

    // SLU-13: non-supernodal (baseline/auto) path.
    // Supernode-aware CSC solve with valid-metadata check and singleton fallback.
    // Mathematically equivalent to fac.solve(b) for baseline/auto factors.
    if (fac.storage_kind_ != sparse_lu_storage_kind::baseline_csc) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_solve_supernode_aware: "
            "unsupported storage kind (only baseline_csc supported)");
    }
    return sparse_lu_detail::solve_baseline_storage_supernode_aware(
        fac.baseline_, fac.info_.n, fac.supernode_info_, b);
}

// ===========================================================================
// SLU-11: Supernode numeric metadata cross-check and diagnostics helpers
// ===========================================================================

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// slu11_verify_row_patterns
//
// Checks that numeric.l_row_ind / numeric.u_row_ind per supernode exactly
// match the sorted unique union of the actual CSC L/U row indices for those
// columns.  Returns false on any mismatch.
//
// Source of truth: actual emitted L.row_ind / U.row_ind.
// NOT derived from symbolic panel_row_ind.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool slu11_verify_row_patterns(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& numeric,
    const baseline_lu_storage<T, Index>& csc_lu)
{
    if (!numeric.valid) return false;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.L, n, n)) return false;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.U, n, n)) return false;

    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index col_begin = supernodes.supernode_ptr[s];
        const Index col_end   = supernodes.supernode_ptr[s + 1u];

        // Expected L rows: sorted unique union over actual L CSC columns
        std::vector<Index> exp_L;
        for (Index j = col_begin; j < col_end; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const Index kb = csc_lu.L.col_ptr[sj];
            const Index ke = csc_lu.L.col_ptr[sj + 1u];
            for (Index k = kb; k < ke; ++k) {
                exp_L.push_back(csc_lu.L.row_ind[static_cast<std::size_t>(k)]);
            }
        }
        std::sort(exp_L.begin(), exp_L.end());
        exp_L.erase(std::unique(exp_L.begin(), exp_L.end()), exp_L.end());

        const Index lb = numeric.l_row_ptr[s];
        const Index le = numeric.l_row_ptr[s + 1u];
        if (static_cast<std::size_t>(le - lb) != exp_L.size()) return false;
        for (std::size_t i = 0u; i < exp_L.size(); ++i) {
            if (numeric.l_row_ind[static_cast<std::size_t>(lb) + i] != exp_L[i])
                return false;
        }

        // Expected U rows: sorted unique union over actual U CSC columns
        std::vector<Index> exp_U;
        for (Index j = col_begin; j < col_end; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const Index kb = csc_lu.U.col_ptr[sj];
            const Index ke = csc_lu.U.col_ptr[sj + 1u];
            for (Index k = kb; k < ke; ++k) {
                exp_U.push_back(csc_lu.U.row_ind[static_cast<std::size_t>(k)]);
            }
        }
        std::sort(exp_U.begin(), exp_U.end());
        exp_U.erase(std::unique(exp_U.begin(), exp_U.end()), exp_U.end());

        const Index ub = numeric.u_row_ptr[s];
        const Index ue = numeric.u_row_ptr[s + 1u];
        if (static_cast<std::size_t>(ue - ub) != exp_U.size()) return false;
        for (std::size_t i = 0u; i < exp_U.size(); ++i) {
            if (numeric.u_row_ind[static_cast<std::size_t>(ub) + i] != exp_U[i])
                return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// slu11_verify_diag_blocks
//
// Checks that numeric.diag_block_values[diag_block_ptr[s] + lc*width + lr]
// equals the actual U CSC value at (col_begin+lr, col_begin+lc), or T(0)
// if that entry is structurally absent.
//
// Layout: column-major (lc outer, lr inner).
// width = supernode_ptr[s+1] - supernode_ptr[s].
// ---------------------------------------------------------------------------
template <class T, class Index>
bool slu11_verify_diag_blocks(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& numeric,
    const baseline_lu_storage<T, Index>& csc_lu)
{
    if (!numeric.valid) return false;
    if (!sparse_lu_is_valid_csc_storage(csc_lu.U, n, n)) return false;

    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index col_begin = supernodes.supernode_ptr[s];
        const Index col_end   = supernodes.supernode_ptr[s + 1u];
        const Index width     = col_end - col_begin;
        const Index blk_b     = numeric.diag_block_ptr[s];

        // column-major: column lc is outer, row lr is inner
        for (Index lc = 0; lc < width; ++lc) {
            const Index j = col_begin + lc;
            const std::size_t sj = static_cast<std::size_t>(j);
            const Index kb = csc_lu.U.col_ptr[sj];
            const Index ke = csc_lu.U.col_ptr[sj + 1u];

            for (Index lr = 0; lr < width; ++lr) {
                const Index i = col_begin + lr;
                bool found = false;
                T expected = T(0);
                for (Index k = kb; k < ke; ++k) {
                    const std::size_t sk = static_cast<std::size_t>(k);
                    if (csc_lu.U.row_ind[sk] == i) {
                        expected = csc_lu.U.values[sk];
                        found = true;
                        break;
                    }
                }
                // SLU-11.1: required diagonal entry U(j,j) must exist in U CSC.
                // i == global row (col_begin+lr), j == global col (col_begin+lc)
                if (i == j && !found) return false;
                const std::size_t idx =
                    static_cast<std::size_t>(blk_b + lc * width + lr);
                if (!(numeric.diag_block_values[idx] == expected)) return false;
            }
        }
    }
    return true;
}

} // namespace sparse_lu_detail

// ---------------------------------------------------------------------------
// sparse_lu_verify_supernode_numeric_against_csc
//
// SLU-11: Cross-checks that sparse_lu_supernode_numeric row patterns and
// diagonal block values are consistent with the actual CSC L/U factors.
//
// Specifically verifies:
//   For each supernode s with columns [col_begin, col_end):
//     numeric.l_row_ind slice == sorted unique union of L.row_ind in those columns
//     numeric.u_row_ind slice == sorted unique union of U.row_ind in those columns
//     diag_block_values (column-major) == actual U entries or zero
//
// Source of truth: actual emitted L.row_ind / U.row_ind.
// Symbolic panel_row_ind is NOT used.
//
// Low-level overload: takes explicit n, symbolic, numeric, and CSC storage.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool sparse_lu_verify_supernode_numeric_against_csc(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& numeric,
    const baseline_lu_storage<T, Index>& csc_lu)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_verify_supernode_numeric_against_csc: Index must be signed");
    if (!numeric.valid) return false;
    if (!supernodes.valid) return false;
    if (!sparse_lu_is_valid_supernode_numeric(n, supernodes, numeric)) return false;
    return sparse_lu_detail::slu11_verify_row_patterns(n, supernodes, numeric, csc_lu) &&
           sparse_lu_detail::slu11_verify_diag_blocks(n, supernodes, numeric, csc_lu);
}

// ---------------------------------------------------------------------------
// sparse_lu_verify_supernode_numeric_against_csc -- factor overload
//
// Convenience overload: extracts n, symbolic, numeric, and CSC storage from
// a sparse_lu_factorization object.  Returns false if the factor does not use
// the supernodal prototype path.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool sparse_lu_verify_supernode_numeric_against_csc(
    const sparse_lu_factorization<T, Index>& fac)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_verify_supernode_numeric_against_csc: Index must be signed");
    if (!fac.valid()) return false;
    if (!fac.uses_supernodal_prototype()) return false;
    return sparse_lu_verify_supernode_numeric_against_csc(
        fac.info().n,
        fac.supernode_info(),
        fac.supernode_numeric_info(),
        fac.baseline_storage());
}

// ---------------------------------------------------------------------------
// sparse_lu_make_supernode_numeric_diagnostics
//
// SLU-11: Builds a sparse_lu_supernode_numeric_diagnostics<Index> from
// explicit supernode numeric metadata and the actual CSC L/U storage.
//
// Populates:
//   valid, n, nsup: from numeric metadata and symbolic partition
//   nnz_L_rows_total, nnz_U_rows_total: from l_row_ind / u_row_ind sizes
//   diag_block_count, diag_block_value_count: from diag_block ptr/values
//   csc_backed: always true for valid SLU-10 numeric metadata
//   actual_csc_consistent: row patterns AND diag blocks consistent with CSC
//   diag_blocks_consistent: diag_block_values consistent with actual U CSC
//   symbolic_panel_rows_used_as_numeric_rows: always false
//
// Low-level overload: takes explicit n, symbolic, numeric, and CSC storage.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_supernode_numeric_diagnostics<Index>
sparse_lu_make_supernode_numeric_diagnostics(
    Index n,
    const sparse_lu_supernode_symbolic<Index>& supernodes,
    const sparse_lu_supernode_numeric<T, Index>& numeric,
    const baseline_lu_storage<T, Index>& csc_lu)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_make_supernode_numeric_diagnostics: Index must be signed");

    sparse_lu_supernode_numeric_diagnostics<Index> d;

    if (!numeric.valid) return d;
    if (!supernodes.valid) return d;
    if (!sparse_lu_is_valid_supernode_numeric(n, supernodes, numeric)) return d;

    const std::size_t nsup = supernodes.supernode_ptr.size() - 1u;

    d.valid   = true;
    d.n       = n;
    d.nsup    = static_cast<Index>(nsup);
    d.nnz_L_rows_total       = static_cast<Index>(numeric.l_row_ind.size());
    d.nnz_U_rows_total       = static_cast<Index>(numeric.u_row_ind.size());
    d.diag_block_count       = static_cast<Index>(nsup);
    d.diag_block_value_count = static_cast<Index>(numeric.diag_block_values.size());

    // CSC-backed: true for valid SLU-10 prototype metadata (always built from CSC)
    d.csc_backed = true;

    // Numeric rows derive from actual CSC L/U; panel_row_ind is symbolic-only
    d.symbolic_panel_rows_used_as_numeric_rows = false;

    // Cross-check row patterns against actual CSC L/U
    const bool row_ok =
        sparse_lu_detail::slu11_verify_row_patterns(n, supernodes, numeric, csc_lu);
    // Cross-check diag block values against actual U CSC
    const bool diag_ok =
        sparse_lu_detail::slu11_verify_diag_blocks(n, supernodes, numeric, csc_lu);

    d.diag_blocks_consistent  = diag_ok;
    d.actual_csc_consistent   = row_ok && diag_ok;

    return d;
}

// ---------------------------------------------------------------------------
// sparse_lu_factorization<T, Index>::supernode_numeric_diagnostics()
//
// SLU-11: Out-of-line definition; requires sparse_lu_make_supernode_numeric_diagnostics
// to be defined before this point.  For baseline/auto factors that do not use
// the supernodal prototype (uses_supernodal_prototype_==false), returns an
// invalid diagnostics struct with valid==false.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_supernode_numeric_diagnostics<Index>
sparse_lu_factorization<T, Index>::supernode_numeric_diagnostics() const
{
    if (!uses_supernodal_prototype_) {
        return sparse_lu_supernode_numeric_diagnostics<Index>();
    }
    return sparse_lu_make_supernode_numeric_diagnostics(
        info_.n, supernode_info_, supernode_numeric_info_, baseline_);
}

// ===========================================================================
// API function implementations
// ===========================================================================

// ---------------------------------------------------------------------------
// sparse_lu_symbolic -- SLU-1: builds the column permutation and col_etree
// skeleton.  SLU-L1 L-1: ordering=auto_select resolves to amd; an identity
// permutation is installed only for explicit ordering=natural.
// Diagnostic variant: always returns sym.success/status
// so that sparse_lu_numeric can propagate failure without extra throw paths.
// ---------------------------------------------------------------------------
template <class Matrix>
sparse_lu_symbolic_result<typename Matrix::index_type>
sparse_lu_symbolic(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    typedef typename Matrix::index_type Index;
    typedef typename Matrix::value_type T;
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");

    sparse_lu_symbolic_result<Index> sym;
    try {
        // Non-square matrices cannot be LU-factored
        if (A.rowsize() != A.columnsize()) {
            sym.success = false;
            sym.status  = sparse_lu_status::invalid_input;
            return sym;
        }
        const Index n = static_cast<Index>(A.rowsize());
        sym.n = n;

        // Option validation: reject explicitly unsupported ordering/pivoting
        const sparse_lu_status opt_st =
            sparse_lu_detail::validate_symbolic_options(opt);
        if (opt_st != sparse_lu_status::success) {
            sym.success = false;
            sym.status  = opt_st;
            return sym;
        }

        // Column pre-permutation slot (§2A).  SLU-L1 L-1 (D-1): auto_select now
        // RESOLVES TO amd (measured dominant on all subject matrices, design
        // §1; the escape hatch is an explicit ordering=natural, which keeps the
        // identity permutation byte-identically).  Ordering Track O1
        // (rcm), O2 (amd) on the A + A^T pattern,
        // and O3 (colamd) on the A^T A column-intersection
        // pattern install a deterministic, pattern-only fill-reducing permutation
        // Q here, leaving the downstream etree / supernode / numeric / solve /
        // storage pipeline unchanged.  Q enters ONLY this column slot; numeric
        // threshold partial pivoting and row_perm are untouched (S-4).
        const csc_storage<T, Index> A_csc_nat = sparse_lu_make_csc_storage(A);
        const sparse_lu_ordering effective_ordering =
            (opt.ordering == sparse_lu_ordering::auto_select)
                ? sparse_lu_ordering::amd : opt.ordering;
        const bool ordering_active =
            (effective_ordering == sparse_lu_ordering::rcm) ||
            (effective_ordering == sparse_lu_ordering::amd) ||
            (effective_ordering == sparse_lu_ordering::colamd) ||
            (effective_ordering == sparse_lu_ordering::nested_dissection);
        if (effective_ordering == sparse_lu_ordering::rcm) {
            // S-1: pattern-only (reads col_ptr/row_ind, never values).
            sym.col_perm = sparse_lu_rcm_ordering(
                n, A_csc_nat.col_ptr, A_csc_nat.row_ind);
        } else if (effective_ordering == sparse_lu_ordering::amd) {
            // S-1: pattern-only (reads col_ptr/row_ind, never values).
            sym.col_perm = sparse_lu_amd_ordering(
                n, A_csc_nat.col_ptr, A_csc_nat.row_ind);
        } else if (effective_ordering == sparse_lu_ordering::colamd) {
            // S-1: pattern-only (reads col_ptr/row_ind, never values).  A^T A is
            // NOT formed explicitly (design §634); rows act as quotient-graph
            // elements.
            sym.col_perm = sparse_lu_colamd_ordering(
                n, A_csc_nat.col_ptr, A_csc_nat.row_ind);
        } else if (effective_ordering == sparse_lu_ordering::nested_dissection) {
            // SLU-MF4: pattern-only (reads col_ptr/row_ind, never values).
            // Recursive graph bisection; separators numbered last -> wide fronts.
            sym.col_perm = sparse_lu_nested_dissection_ordering(
                n, A_csc_nat.col_ptr, A_csc_nat.row_ind);
        } else {
            // natural (explicit): identity permutation, byte-identical legacy.
            sym.col_perm = sparse_lu_identity_permutation(n);
        }
        // inv_col_perm[old] = new (validates Q is a proper permutation).
        sym.inv_col_perm = sparse_lu_inverse_permutation(sym.col_perm);

        // SLU-6: compute A^T A column elimination tree from sparsity pattern.
        // parent[k] = -1 (root) or k < parent[k] < n.
        // SLU-7: also build symbolic supernode metadata from A pattern and etree.
        // Production GP numeric path uses L-structure DFS (independent of col_etree).
        // The symbolic structure is built from the Q-permuted pattern so that it
        // matches the permuted CSC the numeric phase factorizes.
        {
            const csc_storage<T, Index> A_csc_sym =
                ordering_active
                    ? sparse_lu_apply_column_permutation_csc(
                          A_csc_nat, n, n, sym.col_perm)
                    : A_csc_nat;
            sym.col_etree = sparse_lu_elimination_tree_csc(
                n, A_csc_sym.col_ptr, A_csc_sym.row_ind);
            const std::size_t relaxation =
                (opt.method == sparse_lu_method::supernodal)
                    ? opt.supernode_relaxation
                    : 0u;
            sym.supernode_info = sparse_lu_build_supernode_symbolic_csc(
                n, A_csc_sym.col_ptr, A_csc_sym.row_ind, sym.col_etree,
                relaxation);
        }

        // Relaxed supernodes: empty (boundary array superseded by supernode_info)
        sym.relaxed_supernodes.clear();

        sym.estimated_fill_upper_bound = Index(0);
        sym.structural_singularity     = false;

        sym.success = true;
        sym.status  = sparse_lu_status::success;

    } catch (const std::bad_alloc&) {
        sym.success = false;
        sym.status  = sparse_lu_status::memory_allocation_failed;
    }
    return sym;
}

// ---------------------------------------------------------------------------
// sparse_lu_symbolic_with_info -- diagnostic alias for sparse_lu_symbolic.
// Returns the same symbolic result with success/status fields set.
// Provided for API symmetry with sparse_lu_factorize_with_info.
// ---------------------------------------------------------------------------
template <class Matrix>
sparse_lu_symbolic_result<typename Matrix::index_type>
sparse_lu_symbolic_with_info(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    return sparse_lu_symbolic(A, opt);
}

// ---------------------------------------------------------------------------
// O4: effective-matrix numeric driver.
//
// Builds the matrix the baseline numeric factorization runs on and returns the
// factorization result with row_perm/Dr/Dc finalized:
//   - default pivoting: A_eff = column-permuted A (byte-identical to the prior
//     inline path); row_perm = the dynamic pivoting permutation; Dr/Dc empty.
//   - pivoting == static_mc64 (O4): A_eff = P_static * Dr * A * Dc * Qc with a
//     zero-free diagonal; the numeric runs with threshold partial pivoting and
//     the returned row_perm is the static MC64 permutation COMPOSED with the
//     dynamic pivoting (final_row_perm[new] = p_static[dynamic_row_perm[new]]),
//     with Dr/Dc installed in original coordinates.  A structurally singular
//     system (no perfect matching) returns pre_fail with structural_singularity
//     (S-D safe failure); a numerically singular B fails inside the numeric.
// ---------------------------------------------------------------------------
namespace sparse_lu_detail {

template <class Matrix>
struct effective_numeric_outcome {
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    bool                                          pre_fail;
    sparse_lu_status                              pre_status;
    csc_storage<T, Index>                         A_eff;     // matrix that was factored
    baseline_reference_factorize_result<T, Index> num;
    real_type                                     max_abs_A;

    effective_numeric_outcome()
        : pre_fail(false),
          pre_status(sparse_lu_status::success),
          max_abs_A(real_type(0)) {}
};

template <class Matrix>
effective_numeric_outcome<Matrix>
run_effective_baseline_numeric(
    const Matrix& A,
    const sparse_lu_symbolic_result<typename Matrix::index_type>& sym,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;

    effective_numeric_outcome<Matrix> r;
    const Index n = sym.n;

    if (opt.pivoting == sparse_lu_pivoting::static_mc64) {
        // O4: MC64 maximum-weight matching + Dr/Dc scaling -> zero-free diagonal.
        const csc_storage<T, Index> A_nat = sparse_lu_make_csc_storage(A);
        sparse_lu_mc64_transform<T, Index> tf =
            sparse_lu_mc64_make_transform(A_nat, n, sym.col_perm);
        if (!tf.success) {
            // No perfect matching: structurally singular -> safe failure (S-D).
            r.pre_fail   = true;
            r.pre_status = tf.status;
            return r;
        }
        r.A_eff     = tf.B;
        r.max_abs_A = sparse_lu_max_abs_csc(r.A_eff);

        // Numeric runs on the matched matrix with the existing threshold partial
        // pivoting (roadmap §3.2): the static matching gives a zero-free diagonal;
        // dynamic pivoting refines for stability.
        sparse_lu_options<T> opt_internal = opt;
        opt_internal.pivoting = sparse_lu_pivoting::threshold_partial;
        r.num = baseline_sparse_gp_lu_factorize(
            r.A_eff, n, sym.col_perm, sym.inv_col_perm, opt_internal);

        if (r.num.success) {
            // Compose static + dynamic row permutations and install Dr/Dc so the
            // solve un-permutes/un-scales back to the ORIGINAL system (S-A/S-B).
            const std::vector<Index> composed =
                sparse_lu_mc64_compose_row_perm(tf.p_static, r.num.storage.row_perm);
            r.num.storage.row_perm     = composed;
            r.num.storage.inv_row_perm = sparse_lu_inverse_permutation(composed);
            r.num.storage.Dr           = tf.Dr;
            r.num.storage.Dc           = tf.Dc;
        }
        return r;
    }

    // O4.2a: equilibration path (opt-in; default is byte-identical).
    // Computes two-sided max-norm scaling Dr/Dc from original A, builds
    // Aeq = diag(Dr)*A*diag(Dc), then applies the column ordering permutation
    // and runs the existing baseline GP numeric.  Dr/Dc are stored in original
    // (unpermuted) coordinates so solve_baseline_storage's existing un-apply
    // (Step 1: rhs1=Dr*b, Step 6: x=Dc*x_orig) works without modification.
    // static_mc64 is already excluded by validate_symbolic_options.
    if (opt.equilibration) {
        const csc_storage<T, Index> A_nat = sparse_lu_make_csc_storage(A);
        std::vector<T> Dr, Dc;
        sparse_lu_detail::sparse_lu_compute_equilibration(A_nat, n, Dr, Dc);
        const csc_storage<T, Index> A_eq =
            sparse_lu_detail::sparse_lu_apply_equilibration(A_nat, Dr, Dc, n);
        r.A_eff     = sparse_lu_apply_column_permutation_csc(
                          A_eq, n, n, sym.col_perm);
        r.max_abs_A = sparse_lu_max_abs_csc(r.A_eff);
        r.num = baseline_sparse_gp_lu_factorize(
            r.A_eff, n, sym.col_perm, sym.inv_col_perm, opt);
        if (r.num.success) {
            r.num.storage.Dr = Dr;
            r.num.storage.Dc = Dc;
        }
        return r;
    }

    // Default path (pivoting untouched, S-C): byte-identical to the prior inline
    // build + numeric call.
    r.A_eff     = sparse_lu_make_permuted_csc_storage(A, sym);
    r.max_abs_A = sparse_lu_max_abs_csc(r.A_eff);
    r.num = baseline_sparse_gp_lu_factorize(
        r.A_eff, n, sym.col_perm, sym.inv_col_perm, opt);
    return r;
}

// ===========================================================================
// SLU-MF7: GP-less native MC64 supernodal path.
//
// Diagnosis (SLU investigation): method=supernodal currently runs the full
// baseline GP numeric factorization (~60% of wall) whose L/U values are then
// COMPLETELY UNUSED by the native multifrontal driver (it rebuilds the frontal
// matrices from A_eff).  The only essential thing GP supplies is the row
// permutation (pivot order).  The MC64 maximum-product matching already supplies
// a zero-free, large-magnitude diagonal ordering (O4, existing asset), so the GP
// numeric can be skipped entirely on the native path.
//
// Dr/Dc scaling is intentionally FORGONE on the native path.  The native §18.2
// solve indexes Dr/Dc in ORIGINAL coordinates (Step 1: rhs=Dr*b, Step 6:
// x=Dc*x_orig) while the multifrontal residual gate / a_eff_lookup index them in
// POST-PIVOT (new) coordinates; under a non-identity row/column permutation a
// single stored Dr/Dc vector cannot be consistent for both.  The matching alone
// (placing the largest-magnitude entries on the diagonal) carries the dominant
// stability benefit; numerically hard / badly-scaled systems that the matching
// cannot factor are caught by the residual acceptance gate and fall back to the
// full GP+MC64 path (Dr/Dc baked into B, un-applied by the CSC solve), which is
// retained unchanged as the safety net.
// ===========================================================================

// build_mc64_skeleton_storage -- partition-only supernodal storage.
//
// Carries ONLY the supernode partition (from symbolic), the MC64 static row
// permutation as row_perm, and the fill-reducing column ordering as col_perm.
// Dr/Dc are left empty (see header).  Consumed exclusively as the `templ`
// partition source by build_supernodal_self_symbolic_storage, which rebuilds
// row_indices / U_segments / panels from A_eff; panels here stay empty.
// ---------------------------------------------------------------------------
template <class T, class Index>
supernodal_lu_storage<T, Index>
build_mc64_skeleton_storage(
    Index                                       n,
    const sparse_lu_supernode_symbolic<Index>&  sym_info,
    const std::vector<Index>&                   p_static,
    const std::vector<Index>&                   col_perm,
    const std::vector<Index>&                   inv_col_perm)
{
    static_assert(std::is_signed<Index>::value,
                  "build_mc64_skeleton_storage: Index must be signed");
    supernodal_lu_storage<T, Index> skel;
    if (n <= Index(0)) return skel;
    if (!sym_info.valid || sym_info.supernode_ptr.size() < 2u) return skel;
    if (static_cast<Index>(p_static.size()) != n) return skel;

    const std::size_t nsup = sym_info.supernode_ptr.size() - 1u;
    skel.supernodes.resize(nsup);
    for (std::size_t s = 0u; s < nsup; ++s) {
        const Index cb = sym_info.supernode_ptr[s];
        const Index ce = sym_info.supernode_ptr[s + 1u];
        if (cb < Index(0) || ce <= cb ||
            static_cast<std::size_t>(ce) > static_cast<std::size_t>(n)) {
            return supernodal_lu_storage<T, Index>(); // invalid partition -> invalid skeleton
        }
        supernode_desc<Index>& d = skel.supernodes[s];
        d.first_col = cb;
        d.num_cols  = ce - cb;
    }
    skel.row_perm     = p_static;
    skel.inv_row_perm = sparse_lu_inverse_permutation(p_static);
    skel.col_perm     = col_perm;
    skel.inv_col_perm = inv_col_perm;
    // Dr / Dc intentionally empty (matching-only native path).
    skel.valid                   = true;
    skel.source_of_truth_storage = false; // partition skeleton, not a numeric source
    skel.bootstrapped_from_csc   = false;
    skel.true_numeric_source     = false;
    return skel;
}

// native_mc64_supernodal_outcome -- result of the GP-less native attempt.
//   ok:              a complete native supernodal numeric source was produced.
//   matching_failed: MC64 found no perfect matching (structurally singular).
//   storage / tns:   the (possibly failed) self-symbolic storage + numeric stats.
template <class T, class Index>
struct native_mc64_supernodal_outcome {
    bool                            ok;
    bool                            matching_failed;
    supernodal_lu_storage<T, Index> storage;
    supernodal_true_numeric_stats<
        typename vcp::tsparse_scalar::real_type<T>::type> tns;
    native_mc64_supernodal_outcome() : ok(false), matching_failed(false) {}
};

// try_build_native_mc64_supernodal -- MC64 matching + self-symbolic + in-place
// multifrontal numeric, with NO baseline GP numeric.  Returns ok=true only when
// the residual acceptance gate accepts the A_eff-origin storage
// (storage.true_numeric_source == true).  Any failure leaves ok=false so the
// caller can fall through to the existing GP path (safety net).
// ---------------------------------------------------------------------------
template <class Matrix>
native_mc64_supernodal_outcome<typename Matrix::value_type,
                               typename Matrix::index_type>
try_build_native_mc64_supernodal(
    const Matrix& A,
    const sparse_lu_symbolic_result<typename Matrix::index_type>& sym,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    native_mc64_supernodal_outcome<T, Index> out;

    const Index n = sym.n;
    if (n <= Index(0)) return out;            // trivial -> let GP path handle it
    if (!sym.supernode_info.valid) return out;

    // 1. MC64 maximum-product matching: supplies the static row permutation
    //    p_static (zero-free diagonal).  The native path is matching-only -- it
    //    never uses Dr/Dc or the assembled B (MF7 §1.3) -- so SLU-MF8 computes the
    //    matching with the sparse heap-SSP path (sparse_lu_mc64_match_native),
    //    which produces the SAME p_static without building any dense n*n cost
    //    matrix (O(n+nnz) memory).  The dense sparse_lu_mc64_make_transform is left
    //    untouched for the non-native GP+MC64 path that does consume B/Dr/Dc.
    const csc_storage<T, Index> A_nat = sparse_lu_make_csc_storage(A);
    sparse_lu_detail::sparse_lu_mc64_matching<Index> mt =
        sparse_lu_detail::sparse_lu_mc64_match_native(A_nat, n, sym.col_perm);
    if (!mt.success) { out.matching_failed = true; return out; }

    // 2. Native effective matrix = column-permuted A (original rows, no scaling).
    //    The matched diagonal is reached through row_perm = p_static; structurally
    //    identical to the default native path's A_eff.
    const csc_storage<T, Index> A_eff = sparse_lu_make_permuted_csc_storage(A, sym);

    // 3. Partition-only skeleton carrying p_static as the row permutation.
    supernodal_lu_storage<T, Index> skel =
        build_mc64_skeleton_storage<T, Index>(
            n, sym.supernode_info, mt.p_static, sym.col_perm, sym.inv_col_perm);
    if (!skel.valid) return out;

    // 4. Self-symbolic dense-front structure from A_eff (GP-independent).
    supernodal_lu_storage<T, Index> ss =
        build_supernodal_self_symbolic_storage(A_eff, skel, n);
    if (!ss.valid) return out;

    // 5. In-place multifrontal numeric directly from A_eff (no GP L/U).  The
    //    in-place driver does not read csc_lu ((void)csc_lu); pass an empty one.
    baseline_lu_storage<T, Index> empty_csc;
    out.tns = sparse_lu_factorize_supernodal_numeric_source(
                  A_eff, empty_csc, ss, opt);
    out.storage = ss;
    if (!out.tns.success || !ss.true_numeric_source) {
        return out;           // residual / structural failure -> GP fallback
    }
    out.ok = true;
    return out;
}

// ---------------------------------------------------------------------------
// SLU-L1 L-3 (D-3): structural nnz of a supernodal_lu_storage.
//
// Definition: number of STORED entries of the supernodal representation,
// excluding alignment padding (panel rows in [row_count, leading_dimension)
// are padding and are NOT counted).  Conventions match the baseline CSC
// definition (info.nnz_L/nnz_U = row_ind.size()):
//   - L: unit diagonal implicit (not counted); strict lower part of the
//        diagonal block + all off-diagonal panel rows (rows num_cols..
//        row_count-1, every column).
//   - U: explicit diagonal; upper triangle (incl. diagonal) of the diagonal
//        block + all U_segments entries (off-diagonal-block U rows).
// Note (D-3, unverified): the panel is a dense model, so exact-zero fill
// created during elimination is still a stored entry here; MATLAB nnz is
// value-based and could report fewer if such exact zeros exist.  Whether
// they actually occur has NOT been verified.
// ---------------------------------------------------------------------------
template <class T, class Index>
struct supernodal_structural_nnz_result {
    Index nnz_L;
    Index nnz_U;
    supernodal_structural_nnz_result() : nnz_L(Index(0)), nnz_U(Index(0)) {}
};

template <class T, class Index>
supernodal_structural_nnz_result<T, Index>
sparse_lu_supernodal_structural_nnz_(
    const supernodal_lu_storage<T, Index>& storage)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse_lu_supernodal_structural_nnz_: Index must be signed");
    supernodal_structural_nnz_result<T, Index> r;
    std::size_t nl = 0u, nu = 0u;
    for (std::size_t s = 0u; s < storage.supernodes.size(); ++s) {
        const supernode_desc<Index>& d = storage.supernodes[s];
        if (d.num_cols <= Index(0)) continue;
        const std::size_t w    = static_cast<std::size_t>(d.num_cols);
        const std::size_t rcnt = d.row_indices.size();
        // diagonal block: U upper incl. diag = w(w+1)/2; L strict lower = w(w-1)/2.
        nu += (w * (w + 1u)) / 2u;
        nl += (w * (w - 1u)) / 2u;
        // off-diagonal L rows (structural; padding rows beyond rcnt excluded).
        if (rcnt > w) nl += (rcnt - w) * w;
    }
    nu += storage.U_segments.row_ind.size();
    r.nnz_L = static_cast<Index>(nl);
    r.nnz_U = static_cast<Index>(nu);
    return r;
}

} // namespace sparse_lu_detail

// ---------------------------------------------------------------------------
// sparse_lu_numeric -- numeric factorization given symbolic result (SLU-5)
// Propagates symbolic failure status; validates options independently.
// Uses the same sparse GP baseline path as sparse_lu_factorize_with_info.
// ---------------------------------------------------------------------------
template <class Matrix>
sparse_lu_factorization<typename Matrix::value_type, typename Matrix::index_type>
sparse_lu_numeric(
    const Matrix& A,
    const sparse_lu_symbolic_result<typename Matrix::index_type>& sym,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");

    sparse_lu_factorization<T, Index> fac;
    sparse_lu_info<T, Index> info;
    try {
        // Propagate symbolic failure: preserve sym.status (e.g. structural_singularity,
        // invalid_input for non-square or inconsistent symbolic input).
        if (!sym.success) {
            info.status = sym.status;
            fac.set_info_(info);
            return fac;
        }

        // Dimension check: must match symbolic result
        if (A.rowsize() != A.columnsize()) {
            info.status = sparse_lu_status::invalid_input;
            fac.set_info_(info);
            return fac;
        }
        const Index n = static_cast<Index>(A.rowsize());
        if (n != sym.n) {
            info.status = sparse_lu_status::invalid_input;
            fac.set_info_(info);
            return fac;
        }
        info.n = n;

        // Option validation: method, ordering, pivoting, feature flags
        const sparse_lu_status opt_st =
            sparse_lu_detail::validate_symbolic_options(opt);
        if (opt_st != sparse_lu_status::success) {
            info.status = opt_st;
            fac.set_info_(info);
            return fac;
        }

        // SLU-SP1 Phase 2: opt-in supernode_panel numeric (left-looking
        // supernode-panel factorization; design §2).  Default-pivoting path
        // only: static_mc64 (D-6: the new path is MC64-independent) and
        // equilibration are answered with an honest not_implemented.  The
        // result is baseline CSC storage, so the existing solve / IR / LUX
        // consumers work unchanged; method_used reports supernode_panel.
        if (opt.method == sparse_lu_method::supernode_panel) {
            if (opt.pivoting == sparse_lu_pivoting::static_mc64 ||
                opt.equilibration) {
                info.success = false;
                info.status  = sparse_lu_status::not_implemented;
                info.method_used = sparse_lu_method::supernode_panel;
                fac.set_info_(info);
                return fac;
            }
            const csc_storage<T, Index> A_eff =
                sparse_lu_make_permuted_csc_storage(A, sym);
            const real_type max_abs_A = sparse_lu_max_abs_csc(A_eff);
            sparse_lu_detail::supernode_panel_factorize_result<T, Index> pres =
                sparse_lu_detail::supernode_panel_lu_factorize(
                    A_eff, n, sym.col_perm, sym.inv_col_perm, opt);
            if (!pres.success) {
                info.success = false;
                info.status  = pres.status;
                info.method_used = sparse_lu_method::supernode_panel;
                fac.set_info_(info);
                return fac;
            }
            const real_type max_abs_U =
                sparse_lu_max_abs_csc(pres.storage.U);
            const real_type gf = (max_abs_A > real_type(0))
                               ? max_abs_U / max_abs_A
                               : real_type(1);
            fac.set_baseline_storage_with_diagnostics_(n, pres.storage, gf);
            fac.set_supernode_info_(sym.supernode_info);
            {
                sparse_lu_info<T, Index> ninfo = fac.info();
                ninfo.method_used          = sparse_lu_method::supernode_panel;
                ninfo.number_of_supernodes = pres.number_of_supernodes;
                fac.set_info_(ninfo);
            }
            return fac;
        }

        // SLU-MF7: GP-less native MC64 supernodal path.  For method=supernodal
        // with static_mc64 matching + self-symbolic + in-place frontal, build the
        // factor directly from the MC64 matching + multifrontal numeric, skipping
        // the baseline GP numeric entirely.  On any failure (no matching,
        // structural inconsistency, residual rejection) fall through to the
        // existing GP path below (safety net; also re-handles MC64 safe-failure).
        if (opt.method == sparse_lu_method::supernodal &&
            opt.pivoting == sparse_lu_pivoting::static_mc64 &&
            opt.supernodal_self_symbolic &&
            opt.supernodal_inplace_frontal &&
            !opt.supernodal_numeric_diagnostic_leftlooking) {
            sparse_lu_detail::native_mc64_supernodal_outcome<T, Index> mf7 =
                sparse_lu_detail::try_build_native_mc64_supernodal(A, sym, opt);
            if (mf7.ok) {
                sparse_lu_info<T, Index> ninfo;
                ninfo.success = true;
                ninfo.status  = sparse_lu_status::success;
                ninfo.n       = n;
                // SLU-L1 L-3: MF7 native path diagnostics.  Structural nnz per
                // D-3 (padding-excluded stored entries; see the helper's note on
                // the value-based MATLAB nnz difference, unverified).
                // growth_factor is NOT computed on this path (D-4): it stays at
                // its default; computing it would need a max|U| scan (L-2+).
                ninfo.method_used = sparse_lu_method::supernodal;
                {
                    const sparse_lu_detail::supernodal_structural_nnz_result<T, Index>
                        snnz = sparse_lu_detail::sparse_lu_supernodal_structural_nnz_(
                            mf7.storage);
                    ninfo.nnz_L = snnz.nnz_L;
                    ninfo.nnz_U = snnz.nnz_U;
                }
                fac.set_info_(ninfo);
                fac.set_supernode_info_(sym.supernode_info);
                fac.set_supernodal_storage_info_(mf7.storage);
                fac.set_true_numeric_info_(mf7.tns, mf7.storage);
                fac.set_supernodal_solve_info_();
                return fac;
            }
            // else: fall through to the GP path below.
        }

        // O4: build the effective matrix and run the baseline numeric.  For
        // pivoting=static_mc64 this is B = P_static*Dr*A*Dc*Qc (zero-free
        // diagonal) with the row permutation composed and Dr/Dc installed;
        // otherwise it is the column-permuted A exactly as before.
        sparse_lu_detail::effective_numeric_outcome<Matrix> eff =
            sparse_lu_detail::run_effective_baseline_numeric(A, sym, opt);
        if (eff.pre_fail) {
            info.status = eff.pre_status;   // MC64 structural singularity (safe)
            fac.set_info_(info);
            return fac;
        }
        const csc_storage<T, Index>& A_csc = eff.A_eff;
        const real_type max_abs_A = eff.max_abs_A;
        sparse_lu_detail::baseline_reference_factorize_result<T, Index>& num_result =
            eff.num;

        if (num_result.success) {
            const real_type max_abs_U =
                sparse_lu_max_abs_csc(num_result.storage.U);
            const real_type gf = (max_abs_A > real_type(0))
                               ? max_abs_U / max_abs_A
                               : real_type(1);
            fac.set_baseline_storage_with_diagnostics_(n, num_result.storage, gf);
            // SLU-9: store supernode metadata for supernode-aware solve
            fac.set_supernode_info_(sym.supernode_info);
            // SLU-10: if explicit supernodal requested, build prototype numeric metadata
            // from actual CSC L/U (not from symbolic panel_row_ind).
            if (opt.method == sparse_lu_method::supernodal) {
                // SLU-SNA1 P1-B (D-2/D-3/D-4): the supernodal path derives its
                // partition from the ACTUAL GP L pattern instead of the A-pattern
                // approximation (whose empty-set merge degeneracy caused SNQ1/SNQ3).
                // The retained supernode_info is replaced so gate5 / supernode-aware
                // consumers see the partition the panels actually use (design §7).
                // Baseline paths keep the A-pattern info set above (D-3); on builder
                // failure the A-pattern info remains the safe fallback.
                sparse_lu_supernode_symbolic<Index> sn_sym_l =
                    sparse_lu_build_supernode_symbolic_from_l(
                        n, num_result.storage, opt.supernode_relaxation);
                const sparse_lu_supernode_symbolic<Index>& sn_partition =
                    sn_sym_l.valid ? sn_sym_l : sym.supernode_info;
                fac.set_supernode_info_(sn_partition);
                sparse_lu_supernode_numeric<T, Index> sn_num =
                    sparse_lu_detail::build_supernode_numeric_from_csc(
                        n, num_result.storage, sn_partition);
                // SLU-8R.1: Production dense kernel connection.
                // Apply adapter to actual U diagonal blocks from CSC-backed sn_num metadata.
                // At this step, CSC baseline is still numeric source (transitional phase);
                // A_eff-origin step below becomes the genuine numeric source.
                bool kernel_called = false;
                std::size_t kernel_ticks = 0;
                if (sn_num.valid) {
                    sparse_lu_detail::run_dense_kernel_on_sn_blocks(
                        sn_num, kernel_called, kernel_ticks);
                }
                fac.set_supernodal_prototype_info_(sn_num, kernel_called, kernel_ticks);
                // SLU-8R.2: Bootstrap supernodal storage source-of-truth from CSC L/U.
                // bootstrapped_from_csc = true (values from baseline GP, not §17.2).
                // source_of_truth_storage = true (factor owns this storage).
                // true_numeric_source = false at this bootstrap step; set true by A_eff step.
                // storage_kind changes to supernodal once storage is valid.
                const bool add_alignment_padding = (opt.supernode_relaxation > 0u);
                supernodal_lu_storage<T, Index> sn_storage =
                    sparse_lu_detail::bootstrap_supernodal_storage_from_csc(
                        n, num_result.storage, sn_num, add_alignment_padding);
                // [SLU-CLN1 C1, 2026-07-05] The transitional §17.2(A)/(B) prototype
                // pass (run_supernode_panel_leftlooking_update +
                // run_within_panel_factorization + re-bootstrap) that previously ran
                // here was REMOVED: its storage mutations were discarded by an
                // idempotent re-bootstrap (proven bit-exact,
                // sandbox/tmp/cln1_rebootstrap_idem.cpp), so the production
                // true-numeric input below is unchanged.  Production §17.2(A)/(B)
                // diagnostics live in the supernodal_true_numeric_* (tn) counters.
                fac.set_supernodal_storage_info_(sn_storage);
                // SLU-MF2: opt-in self-symbolic dense-front structure for the
                // multifrontal numeric source. Replaces the GP-exact bootstrap
                // structure with the dense multifrontal (relaxed/AMD) fill so wide
                // fronts run native. Non-destructive: only when requested and not
                // running the left-looking diagnostic; keeps the bootstrap as
                // fallback if the builder returns invalid.
                if (sn_storage.valid && opt.supernodal_self_symbolic &&
                    !opt.supernodal_numeric_diagnostic_leftlooking) {
                    supernodal_lu_storage<T, Index> ss_storage =
                        sparse_lu_detail::build_supernodal_self_symbolic_storage(
                            A_csc, sn_storage, n);
                    if (ss_storage.valid) sn_storage = std::move(ss_storage);
                }
                // SLU-8R.5.5: A_eff-origin true numeric source switch.
                // Initializes panel_values/U_segments from A_eff and runs interleaved
                // §17.2(A)/(B) driver. On success, sets sn_storage.true_numeric_source=true
                // and sn_storage.numeric_source_kind=a_eff_true_numeric.
                // On failure (residual too large, zero pivot, etc.): in-place may-modify
                // contract (SLU-8R.5.5.2 Option B): panel_values/U_segments may be
                // modified (zeroed + scatter-filled from A_eff); true_numeric_source
                // stays false, native solve stays disabled, CSC fallback continues.
                if (sn_storage.valid) {
                    // SLU-MF: production numeric source is the multifrontal driver
                    // (left-looking retained as opt-in diagnostic).
                    sparse_lu_detail::supernodal_true_numeric_stats<
                        typename vcp::tsparse_scalar::real_type<T>::type> tns =
                        sparse_lu_factorize_supernodal_numeric_source(
                            A_csc, num_result.storage, sn_storage, opt);
                    // Pass updated sn_storage so that supernodal_ is refreshed with
                    // the new true_numeric_source / numeric_source_kind flags.
                    fac.set_true_numeric_info_(tns, sn_storage);
                }
                // SLU-8R.5: finalize supernodal_solve_native based on storage + status.
                // true_numeric_source is now set by set_true_numeric_info_ if A-origin
                // factorization succeeded, enabling native solve.
                fac.set_supernodal_solve_info_();
            }
        } else {
            info.success = false;
            info.status  = num_result.status;
            fac.set_info_(info);
        }

    } catch (const std::bad_alloc&) {
        info.success = false;
        info.status  = sparse_lu_status::memory_allocation_failed;
        fac.set_info_(info);
    } catch (const vcp::error&) {
        info.success = false;
        info.status  = sparse_lu_status::invalid_input;
        fac.set_info_(info);
    } catch (const std::exception&) {
        // certified ゲート(SLU-GT1 D1/D3)が正しければ到達しない最終防護網。
        // 発火は「ゲートの取りこぼし」= 実装上の予期しない状態を意味する。
        // 数値的特異性は D3 ゲートが throw 前に numerical_singularity /
        // zero_pivot として報告する(例外網を数値失敗の正規経路にしない)。
        info.success = false;
        info.status  = sparse_lu_status::internal_error;
        fac.set_info_(info);
    }
    return fac;
}

// ---------------------------------------------------------------------------
// sparse_lu_factorize_with_info -- diagnostic API (SLU-5 production path)
// Performs: dimension check -> symbolic -> CSC conversion -> sparse GP numeric.
// Catches std::bad_alloc and maps it to memory_allocation_failed.
// ---------------------------------------------------------------------------
template <class Matrix>
sparse_lu_factorization<typename Matrix::value_type, typename Matrix::index_type>
sparse_lu_factorize_with_info(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
    static_assert(std::is_signed<Index>::value,
                  "sparse LU Index must be signed");

    sparse_lu_factorization<T, Index> fac;
    sparse_lu_info<T, Index> info;
    try {
        // Dimension check: non-square matrices are invalid input
        if (A.rowsize() != A.columnsize()) {
            info.status = sparse_lu_status::invalid_input;
            fac.set_info_(info);
            return fac;
        }
        info.n = static_cast<Index>(A.rowsize());

        // Symbolic phase: builds an identity col_perm skeleton for explicit
        // natural and a fill-reducing col_perm for rcm/amd/colamd (SLU-L1 L-1:
        // auto_select resolves to amd inside sparse_lu_symbolic). static_mc64 is
        // handled later in the numeric effective-matrix path.
        sparse_lu_symbolic_result<Index> sym = sparse_lu_symbolic(A, opt);
        if (!sym.success) {
            info.status = sym.status;
            fac.set_info_(info);
            return fac;
        }

        // SLU-SP1 Phase 2: opt-in supernode_panel numeric (same contract and
        // wiring as the matching branch in sparse_lu_numeric -- see there).
        if (opt.method == sparse_lu_method::supernode_panel) {
            if (opt.pivoting == sparse_lu_pivoting::static_mc64 ||
                opt.equilibration) {
                info.success = false;
                info.status  = sparse_lu_status::not_implemented;
                info.method_used = sparse_lu_method::supernode_panel;
                fac.set_info_(info);
                return fac;
            }
            const csc_storage<T, Index> A_eff =
                sparse_lu_make_permuted_csc_storage(A, sym);
            const real_type max_abs_A = sparse_lu_max_abs_csc(A_eff);
            sparse_lu_detail::supernode_panel_factorize_result<T, Index> pres =
                sparse_lu_detail::supernode_panel_lu_factorize(
                    A_eff, info.n, sym.col_perm, sym.inv_col_perm, opt);
            if (!pres.success) {
                info.success = false;
                info.status  = pres.status;
                info.method_used = sparse_lu_method::supernode_panel;
                fac.set_info_(info);
                return fac;
            }
            const real_type max_abs_U =
                sparse_lu_max_abs_csc(pres.storage.U);
            const real_type gf = (max_abs_A > real_type(0))
                               ? max_abs_U / max_abs_A
                               : real_type(1);
            fac.set_baseline_storage_with_diagnostics_(info.n, pres.storage, gf);
            fac.set_supernode_info_(sym.supernode_info);
            {
                sparse_lu_info<T, Index> ninfo = fac.info();
                ninfo.method_used          = sparse_lu_method::supernode_panel;
                ninfo.number_of_supernodes = pres.number_of_supernodes;
                fac.set_info_(ninfo);
            }
            return fac;
        }

        // SLU-MF7: GP-less native MC64 supernodal path (see sparse_lu_numeric and
        // try_build_native_mc64_supernodal).  Skips the baseline GP numeric when
        // method=supernodal + static_mc64 + self-symbolic + in-place frontal; on
        // any failure falls through to the GP path below (safety net).
        if (opt.method == sparse_lu_method::supernodal &&
            opt.pivoting == sparse_lu_pivoting::static_mc64 &&
            opt.supernodal_self_symbolic &&
            opt.supernodal_inplace_frontal &&
            !opt.supernodal_numeric_diagnostic_leftlooking) {
            sparse_lu_detail::native_mc64_supernodal_outcome<T, Index> mf7 =
                sparse_lu_detail::try_build_native_mc64_supernodal(A, sym, opt);
            if (mf7.ok) {
                sparse_lu_info<T, Index> ninfo;
                ninfo.success = true;
                ninfo.status  = sparse_lu_status::success;
                ninfo.n       = info.n;
                // SLU-L1 L-3: MF7 native path diagnostics (same as the matching
                // block in sparse_lu_numeric).  Structural nnz per D-3;
                // growth_factor NOT computed on this path (D-4).
                ninfo.method_used = sparse_lu_method::supernodal;
                {
                    const sparse_lu_detail::supernodal_structural_nnz_result<T, Index>
                        snnz = sparse_lu_detail::sparse_lu_supernodal_structural_nnz_(
                            mf7.storage);
                    ninfo.nnz_L = snnz.nnz_L;
                    ninfo.nnz_U = snnz.nnz_U;
                }
                fac.set_info_(ninfo);
                fac.set_supernode_info_(sym.supernode_info);
                fac.set_supernodal_storage_info_(mf7.storage);
                fac.set_true_numeric_info_(mf7.tns, mf7.storage);
                fac.set_supernodal_solve_info_();
                return fac;
            }
            // else: fall through to the GP path below.
        }

        // O4: build the effective matrix and run the baseline numeric.  Default
        // pivoting -> column-permuted A (byte-identical to before).  static_mc64
        // -> B = P_static*Dr*A*Dc*Qc with zero-free diagonal; row permutation
        // composed (static x dynamic) and Dr/Dc installed so the solve recovers
        // the ORIGINAL system.  No perfect matching -> safe failure (S-D).
        sparse_lu_detail::effective_numeric_outcome<Matrix> eff =
            sparse_lu_detail::run_effective_baseline_numeric(A, sym, opt);
        if (eff.pre_fail) {
            info.status = eff.pre_status;
            fac.set_info_(info);
            return fac;
        }
        const csc_storage<T, Index>& A_csc = eff.A_eff;

        // Track max_abs(A_effective) for growth_factor = max|U| / max|A_eff|.
        const real_type max_abs_A = eff.max_abs_A;

        // SLU-5: sparse GP baseline factorization (production path).
        sparse_lu_detail::baseline_reference_factorize_result<T, Index>& num_result =
            eff.num;

        if (num_result.success) {
            const real_type max_abs_U =
                sparse_lu_max_abs_csc(num_result.storage.U);
            const real_type gf = (max_abs_A > real_type(0))
                               ? max_abs_U / max_abs_A
                               : real_type(1);
            fac.set_baseline_storage_with_diagnostics_(
                info.n, num_result.storage, gf);
            // SLU-9: store supernode metadata for supernode-aware solve
            fac.set_supernode_info_(sym.supernode_info);
            // SLU-10: if explicit supernodal requested, build prototype numeric metadata
            // from actual CSC L/U (not from symbolic panel_row_ind).
            if (opt.method == sparse_lu_method::supernodal) {
                // SLU-SNA1 P1-B (D-2/D-3/D-4): L-pattern partition for the
                // supernodal path; same wiring as sparse_lu_numeric (see there).
                sparse_lu_supernode_symbolic<Index> sn_sym_l =
                    sparse_lu_build_supernode_symbolic_from_l(
                        info.n, num_result.storage, opt.supernode_relaxation);
                const sparse_lu_supernode_symbolic<Index>& sn_partition =
                    sn_sym_l.valid ? sn_sym_l : sym.supernode_info;
                fac.set_supernode_info_(sn_partition);
                sparse_lu_supernode_numeric<T, Index> sn_num =
                    sparse_lu_detail::build_supernode_numeric_from_csc(
                        info.n, num_result.storage, sn_partition);
                // SLU-8R.1: Production dense kernel connection.
                // Apply adapter to actual U diagonal blocks from CSC-backed sn_num metadata.
                // At this step, CSC baseline is still numeric source (transitional phase);
                // A_eff-origin step below becomes the genuine numeric source.
                bool kernel_called = false;
                std::size_t kernel_ticks = 0;
                if (sn_num.valid) {
                    sparse_lu_detail::run_dense_kernel_on_sn_blocks(
                        sn_num, kernel_called, kernel_ticks);
                }
                fac.set_supernodal_prototype_info_(sn_num, kernel_called, kernel_ticks);
                // SLU-8R.2: Bootstrap supernodal storage source-of-truth from CSC L/U.
                // bootstrapped_from_csc = true (values from baseline GP, not §17.2).
                // source_of_truth_storage = true (factor owns this storage).
                // true_numeric_source = false at this bootstrap step; set true by A_eff step.
                // storage_kind changes to supernodal once storage is valid.
                const bool add_alignment_padding = (opt.supernode_relaxation > 0u);
                supernodal_lu_storage<T, Index> sn_storage =
                    sparse_lu_detail::bootstrap_supernodal_storage_from_csc(
                        info.n, num_result.storage, sn_num, add_alignment_padding);
                // [SLU-CLN1 C1, 2026-07-05] Transitional §17.2(A)/(B) prototype pass
                // REMOVED (see the matching block in sparse_lu_numeric): its output
                // was discarded by an idempotent re-bootstrap, so the production
                // true-numeric input below is unchanged.
                fac.set_supernodal_storage_info_(sn_storage);
                // SLU-MF2: opt-in self-symbolic dense-front structure (see the
                // matching block in sparse_lu_factorize). Non-destructive.
                if (sn_storage.valid && opt.supernodal_self_symbolic &&
                    !opt.supernodal_numeric_diagnostic_leftlooking) {
                    supernodal_lu_storage<T, Index> ss_storage =
                        sparse_lu_detail::build_supernodal_self_symbolic_storage(
                            A_csc, sn_storage, info.n);
                    if (ss_storage.valid) sn_storage = std::move(ss_storage);
                }
                // SLU-8R.5.5: A_eff-origin true numeric source switch.
                // In-place may-modify contract (SLU-8R.5.5.2 Option B):
                //   panel_values/U_segments may be modified regardless of outcome.
                //   On failure: true_numeric_source stays false, native solve disabled.
                if (sn_storage.valid) {
                    // SLU-MF: production numeric source is the multifrontal driver
                    // (left-looking retained as opt-in diagnostic).
                    sparse_lu_detail::supernodal_true_numeric_stats<
                        typename vcp::tsparse_scalar::real_type<T>::type> tns =
                        sparse_lu_factorize_supernodal_numeric_source(
                            A_csc, num_result.storage, sn_storage, opt);
                    fac.set_true_numeric_info_(tns, sn_storage);
                }
                // SLU-8R.5: finalize supernodal_solve_native based on storage + status.
                fac.set_supernodal_solve_info_();
            }
        } else {
            info.success = false;
            info.status  = num_result.status;
            fac.set_info_(info);
        }

    } catch (const std::bad_alloc&) {
        info.success = false;
        info.status  = sparse_lu_status::memory_allocation_failed;
        fac.set_info_(info);
    } catch (const vcp::error&) {
        // Covers invalid_argument / index_error / state_error from conversion:
        //   - out-of-range COO entry      → invalid_input
        //   - inconsistent COO buffers    → invalid_input
        //   - duplicate permutation entry → invalid_input
        //   - unsupported finalized format → invalid_input
        info.success = false;
        info.status  = sparse_lu_status::invalid_input;
        fac.set_info_(info);
    } catch (const std::exception&) {
        // certified ゲート(SLU-GT1 D1/D3)が正しければ到達しない最終防護網。
        // 発火は「ゲートの取りこぼし」= 実装上の予期しない状態を意味する。
        // 数値的特異性は D3 ゲートが throw 前に numerical_singularity /
        // zero_pivot として報告する(例外網を数値失敗の正規経路にしない)。
        info.success = false;
        info.status  = sparse_lu_status::internal_error;
        fac.set_info_(info);
    }
    return fac;
}

// ---------------------------------------------------------------------------
// sparse_lu_factorize -- strict API (throws on failure)
// Delegates to factorize_with_info; throws state_error if not success.
// ---------------------------------------------------------------------------
template <class Matrix>
sparse_lu_factorization<typename Matrix::value_type, typename Matrix::index_type>
sparse_lu_factorize(
    const Matrix& A,
    const sparse_lu_options<typename Matrix::value_type>& opt)
{
    typedef typename Matrix::value_type T;
    typedef typename Matrix::index_type Index;

    sparse_lu_factorization<T, Index> fac = sparse_lu_factorize_with_info(A, opt);
    if (!fac.info().success) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_factorize: factorization failed: ",
            sparse_lu_status_to_string(fac.info().status));
    }
    return fac;
}

// ===========================================================================
// O4.1a: iterative refinement (IR) -- solve-time post-process.
// Injected here, after the factorization class, the CSC convert helpers, and
// the baseline solve impl are all available.
// ===========================================================================
#include <vcp/tsparse/detail/tsparse_sparse_lu_iterative_refinement_impl.hpp>

} // namespace vcp

#endif // VCP_TSPARSE_SPARSE_LU_HPP
