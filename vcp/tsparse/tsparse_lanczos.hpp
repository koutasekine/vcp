// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_LANCZOS_HPP
#define VCP_TSPARSE_LANCZOS_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_lanczos {

// ---------------------------------------------------------------------------
// Generate a deterministic starting vector (based on indices)
// ---------------------------------------------------------------------------
template <typename T>
std::vector<T> deterministic_start_vector(const std::size_t n, const unsigned int seed) {
	typedef typename tsparse_scalar::real_type<T>::type R;
	std::vector<T> v(n, T(0));
	// LCG-based pseudorandom (simple, deterministic, no <random> dependency issue)
	unsigned int state = seed + 1u;
	for (std::size_t i = 0; i < n; i++) {
		state = state * 1664525u + 1013904223u;
		// Map to [-1, 1]
		const R val = R(static_cast<int>(state >> 16)) / R(32768) - R(1);
		v[i] = T(val);
	}
	const R nv = tsparse_scalar::real_norm_value(v);
	if (nv > R(0)) {
		for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	}
	return v;
}

// ---------------------------------------------------------------------------
// Result package of a Lanczos-type core (kept for callers; the plain-restart
// cores lanczos_eigs / lanczos_eigs_standard / lanczos_eigs_standard_si_ that
// produced it were removed in EIG-T1 T-6 -- the thick-restart kernel
// tsparse_thick_restart_lanczos.hpp replaced them).
// ---------------------------------------------------------------------------
template <typename T, class ApplyA>
struct lanczos_result_package {
	std::vector<T> eigenvalues;
	std::vector<std::vector<T> > eigenvectors;
		std::vector<typename tsparse_scalar::real_type<T>::type> residuals_abs;
		std::vector<typename tsparse_scalar::real_type<T>::type> residuals_rel;
		std::vector<typename tsparse_scalar::real_type<T>::type> history_abs;
		std::vector<typename tsparse_scalar::real_type<T>::type> history_rel;
	bool converged;
	std::size_t converged_count;
	std::size_t returned_count;
	std::size_t iterations;
	std::size_t restarts;
	std::size_t mv_count;
	std::string breakdown_reason;
	std::string failure_reason;
	std::string used_method;
};

} // namespace tsparse_lanczos
} // namespace vcp

#endif
