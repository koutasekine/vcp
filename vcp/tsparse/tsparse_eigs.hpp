// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_EIGS_HPP
#define VCP_TSPARSE_EIGS_HPP

namespace vcp {

// ---------------------------------------------------------------------------
// Solver method
// ---------------------------------------------------------------------------
enum class eig_solver_method {
	lanczos,
	arnoldi,
	shift_invert_lanczos,
	shift_invert_arnoldi,
	dense_fallback_explicit
};

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
