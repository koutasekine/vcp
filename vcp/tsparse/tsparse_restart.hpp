// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_RESTART_HPP
#define VCP_TSPARSE_RESTART_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_honest_termination.hpp>
#include <vcp/spmatrix.hpp>

namespace vcp {
namespace tsparse {

// ---------------------------------------------------------------------------
// ritz_pair<T>
// Active Ritz pair produced during a Krylov iteration.
// ---------------------------------------------------------------------------
template <class T>
struct ritz_pair {
	typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

	T value;
	std::vector<T> vector;
	real_type residual_absolute;
	real_type residual_relative;
	bool converged;

	// SLU-GT1 D5: residual fields initialized to real_type(0); `converged`
	// (and the writer that computes the residual) is the validity witness.
	// Before a residual is written these fields are undefined -- do not read.
	ritz_pair()
		: value(T(0)),
		  vector(),
		  residual_absolute(real_type(0)),
		  residual_relative(real_type(0)),
		  converged(false) {}
};

// ---------------------------------------------------------------------------
// locked_pair<T>
// A Ritz pair that has been declared converged and moved to the locked set.
// ---------------------------------------------------------------------------
template <class T>
struct locked_pair {
	typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

	T value;
	std::vector<T> vector;
	real_type residual_absolute;
	real_type residual_relative;

	// SLU-GT1 D5: a locked pair is by construction converged; its residual
	// fields are written at lock time (see lock_converged_pairs).  The
	// real_type(0) initializer is a placeholder, not a sentinel -- a
	// default-constructed locked_pair must not have its residuals read.
	locked_pair()
		: value(T(0)),
		  vector(),
		  residual_absolute(real_type(0)),
		  residual_relative(real_type(0)) {}
};

// ---------------------------------------------------------------------------
// krylov_restart_state<T>
// Full state for a restart-capable Krylov solver.
// ---------------------------------------------------------------------------
template <class T>
struct krylov_restart_state {
	std::vector<std::vector<T> > basis;
	std::vector<ritz_pair<T> > active_ritz_pairs;
	std::vector<locked_pair<T> > locked_pairs;
	std::size_t restarts;
	std::size_t matrix_vector_products;

	krylov_restart_state() : basis(), active_ritz_pairs(), locked_pairs(),
		restarts(0), matrix_vector_products(0) {}

	void clear() {
		basis.clear();
		active_ritz_pairs.clear();
		locked_pairs.clear();
		restarts = 0;
		matrix_vector_products = 0;
	}
};

// ---------------------------------------------------------------------------
// vector_norm2
// Returns the 2-norm of x using the Hermitian inner product.
// Returns 0 if norm is non-positive (never NaN).
// ---------------------------------------------------------------------------
template <class T>
typename vcp::tsparse_scalar::real_type<T>::type
vector_norm2(const std::vector<T>& x) {
	return vcp::tsparse_scalar::hermitian_norm_value(x);
}

// ---------------------------------------------------------------------------
// scale_vector
// Multiplies every element of x by alpha (real scalar).
// ---------------------------------------------------------------------------
template <class T>
void scale_vector(std::vector<T>& x,
                  const typename vcp::tsparse_scalar::real_type<T>::type& alpha)
{
	for (std::size_t i = 0; i < x.size(); i++) {
		x[i] = T(alpha) * x[i];
	}
}

// ---------------------------------------------------------------------------
// orthogonalize_against_locked
//
// Projects v orthogonal to every locked vector via full (double-pass)
// Hermitian Gram-Schmidt using the correct projection formula
//   v <- v - (<q, v> / <q, q>) * q
// then optionally normalizes v.
//
// Returns the 2-norm of v after orthogonalization (before normalization).
// Returns 0 if v is near-zero after projection; never returns NaN.
//
// Throws vcp::dimension_error   if any locked vector size != v.size().
// Throws vcp::numerical_error   if any locked vector has zero or near-zero norm.
// ---------------------------------------------------------------------------
template <class T>
typename vcp::tsparse_scalar::real_type<T>::type
orthogonalize_against_locked(
	std::vector<T>& v,
	const std::vector<locked_pair<T> >& locked,
	bool normalize)
{
	typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
	const std::size_t n = v.size();
	const real_type eps = vcp::tsparse_scalar::epsilon<real_type>();

	// Pre-check dimensions and precompute <qk, qk> for each locked vector.
	std::vector<real_type> q_norms2(locked.size());
	for (std::size_t k = 0; k < locked.size(); k++) {
		if (locked[k].vector.size() != n) {
			vcp::throw_error<vcp::dimension_error>(
				"tsparse::orthogonalize_against_locked: dimension mismatch");
		}
		real_type q_norm2 = real_type(0);
		const std::vector<T>& qk = locked[k].vector;
		for (std::size_t i = 0; i < n; i++) {
			q_norm2 += vcp::tsparse_scalar::real_part(
				vcp::tsparse_scalar::conjugate_if_needed(qk[i]) * qk[i]);
		}
		if (!(q_norm2 > eps)) {
			vcp::throw_error<vcp::numerical_error>(
				"tsparse::orthogonalize_against_locked: locked vector has zero or near-zero norm");
		}
		q_norms2[k] = q_norm2;
	}

	// Two-pass orthogonalization: v <- v - (<qk, v> / <qk, qk>) * qk
	for (int pass = 0; pass < 2; pass++) {
		for (std::size_t k = 0; k < locked.size(); k++) {
			const std::vector<T>& qk = locked[k].vector;
			T dot = T(0);
			for (std::size_t i = 0; i < n; i++) {
				dot += vcp::tsparse_scalar::conjugate_if_needed(qk[i]) * v[i];
			}
			const T coeff = dot / T(q_norms2[k]);
			for (std::size_t i = 0; i < n; i++) {
				v[i] -= coeff * qk[i];
			}
		}
	}

	const real_type nrm = vector_norm2(v);
	if (normalize) {
		if (nrm > real_type(0)) {
			const real_type inv = real_type(1) / nrm;
			scale_vector(v, inv);
		}
	}
	return nrm;
}

// ---------------------------------------------------------------------------
// deflate_against_locked
//
// Orthogonalizes v against all locked vectors and normalizes.
// Thin wrapper around orthogonalize_against_locked(v, locked, true).
// ---------------------------------------------------------------------------
template <class T>
typename vcp::tsparse_scalar::real_type<T>::type
deflate_against_locked(
	std::vector<T>& v,
	const std::vector<locked_pair<T> >& locked)
{
	return orthogonalize_against_locked(v, locked, true);
}

// ---------------------------------------------------------------------------
// lock_converged_pairs
//
// Moves converged active Ritz pairs to the locked set.
// - Moves at most max_to_lock pairs in a single call.
// - Pairs are moved in their original order within active.
// - Repeated eigenvalues are preserved (multiplicity is never collapsed).
// - Returns the number of newly locked pairs.
// ---------------------------------------------------------------------------
template <class T>
std::size_t lock_converged_pairs(
	std::vector<ritz_pair<T> >& active,
	std::vector<locked_pair<T> >& locked,
	std::size_t max_to_lock)
{
	if (max_to_lock == 0) return 0;

	std::size_t locked_count = 0;
	std::vector<ritz_pair<T> > remaining;
	remaining.reserve(active.size());

	for (std::size_t i = 0; i < active.size(); i++) {
		if (active[i].converged && locked_count < max_to_lock) {
			locked_pair<T> lp;
			lp.value = active[i].value;
			lp.vector = active[i].vector;
			lp.residual_absolute = active[i].residual_absolute;
			lp.residual_relative = active[i].residual_relative;
			locked.push_back(lp);
			locked_count++;
		} else {
			remaining.push_back(active[i]);
		}
	}

	active.swap(remaining);
	return locked_count;
}

// ---------------------------------------------------------------------------
// take_locked_prefix
//
// Returns the first k locked pairs (or all of them if k >= locked.size()).
// - Input order is preserved.
// - Repeated eigenvalues are never collapsed.
// ---------------------------------------------------------------------------
template <class T>
std::vector<locked_pair<T> > take_locked_prefix(
	const std::vector<locked_pair<T> >& locked,
	std::size_t k)
{
	const std::size_t take = (k < locked.size()) ? k : locked.size();
	return std::vector<locked_pair<T> >(locked.begin(), locked.begin() + static_cast<std::ptrdiff_t>(take));
}

// ---------------------------------------------------------------------------
// honest_termination_check_   (EIG-1 F-3-1 / F-4; EIG-0 C-2)
//
// コア実装は tsparse_honest_termination.hpp の 1 箇所のみ(コピー禁止)。
// ここには locked_pair / ritz_pair を扱う thin overload だけを置く。
// ---------------------------------------------------------------------------
template <class T>
bool honest_termination_check_(
	const std::vector<T>& active_values,
	const std::vector<bool>& active_converged,
	const std::vector<locked_pair<T> >& locked,
	const std::size_t k,
	const eig_target target,
	const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
	std::vector<T> locked_values;
	locked_values.reserve(locked.size());
	for (std::size_t i = 0; i < locked.size(); i++) {
		locked_values.push_back(locked[i].value);
	}
	return honest_termination_check_(
		active_values, active_converged, locked_values, k, target, shift);
}

// ritz_pair 版の thin overload(TRL / lanczos 系の active 集合をそのまま渡す)
template <class T>
bool honest_termination_check_(
	const std::vector<ritz_pair<T> >& active,
	const std::vector<locked_pair<T> >& locked,
	const std::size_t k,
	const eig_target target,
	const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
	std::vector<T> values;
	std::vector<bool> conv;
	values.reserve(active.size());
	conv.reserve(active.size());
	for (std::size_t i = 0; i < active.size(); i++) {
		values.push_back(active[i].value);
		conv.push_back(active[i].converged);
	}
	return honest_termination_check_(values, conv, locked, k, target, shift);
}

// residual_acceptance_check_(EIG-0 C-1)のコア実装も
// tsparse_honest_termination.hpp にある(同 namespace vcp::tsparse)。

// ---------------------------------------------------------------------------
// append_locked_to_result
//
// Appends all locked pairs to the existing eig_result<T>.
// Updates eigenvalues, eigenvectors, residuals_absolute, residuals_relative,
// returned_count, returned_real_count, converged_count, and residual norms.
//
// Does NOT add new fields to eig_result<T>.
// ---------------------------------------------------------------------------
template <class T>
void append_locked_to_result(
	const std::vector<locked_pair<T> >& locked,
	vcp::eig_result<T>& result)
{
	typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

	for (std::size_t i = 0; i < locked.size(); i++) {
		result.eigenvalues.push_back(locked[i].value);
		result.eigenvectors.push_back(locked[i].vector);
		result.residuals_absolute.push_back(locked[i].residual_absolute);
		result.residuals_relative.push_back(locked[i].residual_relative);
		result.returned_count++;
		result.returned_real_count++;
		result.converged_count++;
	}

	// Recompute residual norms as max over all returned pairs.
	if (!result.residuals_absolute.empty()) {
		real_type max_abs = result.residuals_absolute[0];
		real_type max_rel = result.residuals_relative.empty()
			? real_type(0) : result.residuals_relative[0];
		for (std::size_t i = 1; i < result.residuals_absolute.size(); i++) {
			if (result.residuals_absolute[i] > max_abs)
				max_abs = result.residuals_absolute[i];
		}
		for (std::size_t i = 1; i < result.residuals_relative.size(); i++) {
			if (result.residuals_relative[i] > max_rel)
				max_rel = result.residuals_relative[i];
		}
		result.residual_norm_absolute = max_abs;
		result.residual_norm_relative = max_rel;
	}

	if (result.converged_count > 0 &&
	    result.converged_count >= result.requested_count) {
		result.converged = true;
	}
}

} // namespace tsparse
} // namespace vcp

#endif
