// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_EIGS_HPP
#define VCP_TSPARSE_EIGS_HPP

#include <cstddef>

namespace vcp {

// ---------------------------------------------------------------------------
// Solver method
// ---------------------------------------------------------------------------
enum class eig_solver_method {
	lanczos,
	// NOTE (EIG-3 D3-1): `arnoldi` is an alias of the rebuilt Krylov-Schur
	// driver (used_method = "arnoldi(krylov_schur)").  The old Ritz-restart
	// arnoldi core was deleted in EIG-3; no deprecation is planned in EIG-4.
	arnoldi,
	shift_invert_lanczos,
	shift_invert_arnoldi,
	dense_fallback_explicit,
	// EIG-4 T-1: explicit low-level selection of the promoted solvers.
	thick_restart_lanczos,
	krylov_schur,
	// EIG-4 T-2 (D4-2): automatic routing -- the eig_options<T> default.
	//   use_shift == true                  -> shift-invert (E-A1, unchanged)
	//   n <= eig_auto_dense_threshold      -> dense (jacobi / real Schur)
	//   certainly symmetric (B-28)         -> thick_restart_lanczos
	//   otherwise (incl. uncertifiable)    -> krylov_schur
	auto_select
};

// ---------------------------------------------------------------------------
// EIG-4 T-2 tuning constants (B-30: single named definition, no per-solver
// hardcoding; rationale = sandbox/docs/reports/EIG-4_G2.1_ndense_measurement.md)
// ---------------------------------------------------------------------------
// auto_select routes to the dense fallback when n <= this threshold.  There is
// no wall-clock crossover vs the iterative solvers (they are faster whenever
// they converge); the threshold bounds the region where the default path is
// CERTAIN to produce an answer at acceptable worst-case cost (~4.4s at n=300).
constexpr std::size_t eig_auto_dense_threshold = 300;
// Jacobi rotation-budget floor for the auto dense SYMMETRIC branch:
// budget = max(options.max_iter, eig_auto_jacobi_budget_factor * n * n).
// (~10 sweeps; measured need at n=300 is <= 2e5 rotations, floor = 4.5e5.)
// The dense route is mv-free (outside the B-1 matrix-vector budget); on
// adversarial dense spectra exceeding the floor the outcome is an HONEST
// not_converged.  The explicit dense_fallback_explicit path keeps the legacy
// max_iter interpretation unchanged (B-26).
constexpr std::size_t eig_auto_jacobi_budget_factor = 5;

// backward compat alias
typedef eig_solver_method eig_method;

// ---------------------------------------------------------------------------
// Eigenvalue target
// ---------------------------------------------------------------------------
enum class eig_target {
	largest_magnitude,
	smallest_magnitude,
	largest_algebraic,
	smallest_algebraic,
	target_magnitude,
	target_real
};

// backward compat alias
typedef eig_target eigs_target;

// ---------------------------------------------------------------------------
// Orthogonalization method (extended)
// ---------------------------------------------------------------------------
enum class orthogonalization_method {
	modified_gram_schmidt,
	classical_gram_schmidt_twice
};

// ---------------------------------------------------------------------------
// Matrix structure hint
// ---------------------------------------------------------------------------
enum class matrix_structure_hint {
	auto_detect,
	general,
	symmetric,
	hermitian
};

// ---------------------------------------------------------------------------
// Generalized eig method (kept for backward compat)
// ---------------------------------------------------------------------------
enum class generalized_eig_method {
	diagonal_or_operator
};

} // namespace vcp

#endif
