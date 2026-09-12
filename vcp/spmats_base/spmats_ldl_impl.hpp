// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_ldl_impl.hpp
// Out-of-line policy method definitions for the sparse LDL^T factorization
// (LDL-3; design v2 SS5).  This file is included inside spmats.hpp AFTER the
// closing brace of spmats<_T,_Index>, alongside spmats_lss.hpp.  The types
// live in spmats_base/spmats_ldl.hpp (included before the class body).

#ifndef VCP_SPMATS_LDL_IMPL_HPP
#define VCP_SPMATS_LDL_IMPL_HPP

#include <exception>
#include <type_traits>
#include <vector>

#include <vcp/spmats_base/spmats_ldl.hpp>

namespace vcp {

namespace spmats_ldl_detail {

	// -----------------------------------------------------------------------
	// dispatch_sparse_ldl_: SFINAE-guarded helper, same pattern as
	// spmats_lss_detail::dispatch_sparse_lu_ -- the signed-Index body is never
	// instantiated for unsigned Index, and the unsigned path reports exactly
	// like the SLU dispatch (throws vcp::state_error).
	// -----------------------------------------------------------------------

	// copy status + integer diagnostics from the tsparse result
	template <typename _T, typename _Index>
	inline ldl_result<_T, _Index> make_ldl_result_(
	    const sparse_ldl_result<_T, _Index>& r)
	{
		ldl_result<_T, _Index> out;
		out.status              = r.status;
		out.n_pivots_1x1        = r.n_pivots_1x1;
		out.n_pivots_2x2        = r.n_pivots_2x2;
		out.first_zero_pivot    = r.first_zero_pivot;
		out.inconclusive_at     = r.inconclusive_at;
		out.structural_empty_at = r.structural_empty_at;
		out.nnz_L               = r.nnz_L;
		out.ordering_used       = r.ordering_used;
		out.method_used         = r.method_used;
		out.dense_delegated     = r.dense_delegated;
		// SLDL-SP diagnostics conduit (design v1 SS4: wired from the start,
		// not retrofitted).  Every field is copied verbatim; the policy layer
		// adds no interpretation of its own.
		out.pivot_mode_used     = r.pivot_mode_used;
		out.diag_kernel_used    = r.diag_kernel_used;
		out.n_supernodes        = r.n_supernodes;
		out.max_supernode_width = r.max_supernode_width;
		out.n_boundary_splits   = r.n_boundary_splits;
		out.nnz_L_static        = r.nnz_L_static;
		out.n_zero_skips        = r.n_zero_skips;
		out.out_of_panel_at     = r.out_of_panel_at;
		out.gemm_call_count     = r.gemm_call_count;
		out.gemm_time_ns        = r.gemm_time_ns;
		out.growth_log2         = r.growth_log2;
		out.growth_valid        = r.growth_valid;
		return out;
	}

	// signed Index path: spmats -> CSC -> sparse_ldl_factorize_with_info ->
	// L / D spmats construction (design v2 SS5.2) + perm.
	template <typename _T, typename _Index>
	inline ldl_result<_T, _Index>
	dispatch_sparse_ldl_(
	    const spmats<_T, _Index>& A,
	    spmats<_T, _Index>& L,
	    spmats<_T, _Index>& D,
	    std::vector<_Index>& perm,
	    const ldl_options<_T>& opt)
	{
		const _Index n = A.rowsize();

		// spmats -> CSC via the SLU conversion helper (include-only reuse;
		// same pattern as the sparse_lu dispatch path).
		const csc_storage<_T, _Index> C = sparse_lu_make_csc_storage(A);

		const sparse_ldl_result<_T, _Index> r =
		    sparse_ldl_factorize_with_info<_T, _Index>(
		        n, C.col_ptr, C.row_ind, C.values, opt);

		ldl_result<_T, _Index> out = make_ldl_result_<_T, _Index>(r);

		// L / D / perm are valid outputs only for success / zero_pivot
		// (design v2 SS4.2 contract); otherwise the out parameters are
		// left empty.
		L.resize(_Index(0), _Index(0));
		D.resize(_Index(0), _Index(0));
		perm.clear();
		if (r.status != sparse_ldl_status::success &&
		    r.status != sparse_ldl_status::zero_pivot) {
			return out;
		}

		// L: kernel CSC is column-sorted, unit diagonal explicit, certified
		// zeros already dropped -> assign_csc accepts it directly.
		L.assign_csc(n, n, r.L_col_ptr, r.L_row_ind, r.L_val);

		// D: symmetric block diagonal from the internal triple; strictly
		// zero values are structural zeros and are NOT stored (design v2
		// SS5.2: subdiagonal stored on BOTH sides; certified-zero drop uses
		// the same x == T(0) rule as the kernels).
		{
			std::vector<_Index> d_ptr(static_cast<std::size_t>(n) + 1u, _Index(0));
			std::vector<_Index> d_ind;
			std::vector<_T>     d_val;
			const std::size_t un = static_cast<std::size_t>(n);
			for (std::size_t k = 0; k < un; ++k) {
				// column k entries, ascending row order
				if (k > 0 && r.D_block2[k - 1u] != char(0)) {
					// k is the SECOND column of a 2x2 block: (k-1,k) = e
					const _T& e = r.D_sub[k - 1u];
					if (e == _T(0)) { /* certified zero: not stored */ }
					else {
						d_ind.push_back(static_cast<_Index>(k - 1u));
						d_val.push_back(e);
					}
				}
				const _T& dk = r.D_diag[k];
				if (dk == _T(0)) { /* certified zero: structural zero */ }
				else {
					d_ind.push_back(static_cast<_Index>(k));
					d_val.push_back(dk);
				}
				if (r.D_block2[k] != char(0) && k + 1u < un) {
					// k leads a 2x2 block: (k+1,k) = e
					const _T& e = r.D_sub[k];
					if (e == _T(0)) { /* certified zero: not stored */ }
					else {
						d_ind.push_back(static_cast<_Index>(k + 1u));
						d_val.push_back(e);
					}
				}
				d_ptr[k + 1u] = static_cast<_Index>(d_ind.size());
			}
			D.assign_csc(n, n, d_ptr, d_ind, d_val);
		}

		perm = r.perm;
		return out;
	}


} // namespace spmats_ldl_detail

// ---------------------------------------------------------------------------
// policy_ldl_with_info: NVI outer (non-virtual).  Finalize guarantee (same
// auto-finalize as the existing NVI outers) + squareness entry check, then
// delegates to the virtual policy_ldl_with_info_impl.  Must never be
// overridden -- override policy_ldl_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ldl_result<_T, _Index> spmats<_T, _Index>::policy_ldl_with_info(
	spmats<_T, _Index>& L,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const ldl_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_ldl_with_info: matrix must be square");
	return policy_ldl_with_info_impl(L, D, perm, opt);
}

// ---------------------------------------------------------------------------
// policy_ldl_with_info_impl: virtual algorithm body (default: signed-Index
// guard -> CSC conversion -> tsparse factorization -> L/D/perm assembly).
// Runtime failure is a status; misuse (vcp::error) keeps its throwing
// contract; the final std::exception net maps to internal_error (P3).
// No scalar-type branching exists on this path (design v2 SS6).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ldl_result<_T, _Index> spmats<_T, _Index>::policy_ldl_with_info_impl(
	spmats<_T, _Index>& L,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const ldl_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	try {
		return spmats_ldl_detail::dispatch_sparse_ldl_<_T, _Index>(A, L, D, perm, opt);
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		// P3 final protection net: firing means a certified gate was missed.
		ldl_result<_T, _Index> out;
		out.status = sparse_ldl_status::internal_error;
		return out;
	}
}

// ===========================================================================
// inertia (LDL-4, design v2 SS7)
// ===========================================================================

// ---------------------------------------------------------------------------
// policy_inertia_from_block_diagonal: non-virtual shared D-consumer.
// Scans a general symmetric 1x1/2x2 block diagonal matrix (BK origin is NOT
// assumed).  SS7.2 rules:
//   1. any stored nonzero with |i-j| > 1        -> not_block_diagonal
//   2. subdiagonal e_k marks a 2x2 block (k,k+1): certified sign counting on
//      det = d_k d_{k+1} - e_k^2 (unstored diagonals read as 0) and trace
//   3. lone diagonal d_k: certified |d_k| <= tol -> n_zero; certified sign
//      -> n_pos / n_neg
//   4. structurally empty column                 -> 1x1 zero block (n_zero)
//   5. self-check n_pos + n_neg + n_zero == n    (violation -> internal_error)
// Counting (success declaration) happens ONLY on certified comparisons; an
// undecidable case returns inconclusive_sign with its position (P1).  For
// totally ordered scalars the inconclusive branch is unreachable.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
inertia_result<_Index> spmats<_T, _Index>::policy_inertia_from_block_diagonal(
	const scalar_real_type& tol) const
{
	const spmats<_T, _Index>& D = *this;   // WFIX-2 (W2-3): the scanned D is *this
	using std::abs;
	typedef scalar_real_type R;

	inertia_result<_Index> out;
	out.status = inertia_status::invalid_input;
	if (D.rowsize() != D.columnsize() || D.rowsize() < _Index(0)) return out;
	if (!D.is_finalized()) D.finalize();

	const std::size_t un = static_cast<std::size_t>(D.rowsize());

	// gather the band: d[k] = (k,k), e[k] = 2x2 marker/value for (k+1,k)
	// (either stored triangle marks the block; the lower value wins when both
	// are stored -- the D produced by policy_ldl_with_info is symmetric).
	std::vector<_T> dv(un, _T(0)), ev(un > 0 ? un - 1u : 0u, _T(0));
	std::vector<char> dhas(un, char(0)), ehas(un > 0 ? un - 1u : 0u, char(0));
	{
		const std::vector<_Index>& outer = D.outer_index();
		const std::vector<_Index>& inner = D.inner_index();
		const std::vector<_T>& value = D.values();
		// format-agnostic for the symmetric band scan: an entry (o, inner)
		// is (row, col) for CSR and (col, row) for CSC; both orientations
		// classify identically by |i - j|.
		for (std::size_t o = 0; o < un; ++o) {
			for (_Index q = outer[o]; q < outer[o + 1u]; ++q) {
				const std::size_t i = static_cast<std::size_t>(inner[static_cast<std::size_t>(q)]);
				const _T& v = value[static_cast<std::size_t>(q)];
				if (i == o) { dv[o] = v; dhas[o] = char(1); }
				else if (i + 1u == o) {   // (o, o-1): subdiagonal in one orientation
					if (ehas[i] == char(0)) { ev[i] = v; ehas[i] = char(1); }
				}
				else if (o + 1u == i) {   // (o, o+1): subdiagonal in the other
					if (ehas[o] == char(0)) { ev[o] = v; ehas[o] = char(1); }
				}
				else {
					out.status = inertia_status::not_block_diagonal;
					return out;
				}
			}
		}
	}

	_Index npos = _Index(0), nneg = _Index(0), nzero = _Index(0);
	std::size_t k = 0;
	while (k < un) {
		if (k + 1u < un && ehas[k] != char(0)) {
			// ---- 2x2 block (k, k+1); unstored diagonals read as 0
			const _T dk  = dv[k];
			const _T dk1 = dv[k + 1u];
			const _T e   = ev[k];
			const _T det = dk * dk1 - e * e;
			const _T tr  = dk + dk1;
			// sign decisions go through R = real_type<T> (identity for real
			// and interval scalars; keeps the complex instantiation of the
			// virtual _impl compilable -- LDL^H itself is out of scope)
			const R rdet = vcp::tsparse_scalar::real_part(det);
			const R rtr  = vcp::tsparse_scalar::real_part(tr);
			const R absdet = abs(det);
			if (absdet <= tol) {
				// certified zero det: one zero eigenvalue; the other has the
				// certified sign of the trace
				nzero = nzero + _Index(1);
				const R abstr = abs(tr);
				if (abstr <= tol) { nzero = nzero + _Index(1); }
				else if (rtr > R(0)) { npos = npos + _Index(1); }
				else if (rtr < R(0)) { nneg = nneg + _Index(1); }
				else {
					out.status = inertia_status::inconclusive_sign;
					out.inconclusive_at = static_cast<_Index>(k);
					return out;
				}
			} else if (rdet < R(0)) {
				npos = npos + _Index(1);
				nneg = nneg + _Index(1);
			} else if (rdet > R(0)) {
				if (rtr > R(0)) { npos = npos + _Index(2); }
				else if (rtr < R(0)) { nneg = nneg + _Index(2); }
				else {
					out.status = inertia_status::inconclusive_sign;
					out.inconclusive_at = static_cast<_Index>(k);
					return out;
				}
			} else {
				out.status = inertia_status::inconclusive_sign;
				out.inconclusive_at = static_cast<_Index>(k);
				return out;
			}
			k += 2;
		} else if (dhas[k] != char(0)) {
			// ---- lone 1x1 diagonal
			const _T dk = dv[k];
			const R rdk = vcp::tsparse_scalar::real_part(dk);
			const R absd = abs(dk);
			if (absd <= tol) { nzero = nzero + _Index(1); }
			else if (rdk > R(0)) { npos = npos + _Index(1); }
			else if (rdk < R(0)) { nneg = nneg + _Index(1); }
			else {
				out.status = inertia_status::inconclusive_sign;
				out.inconclusive_at = static_cast<_Index>(k);
				return out;
			}
			k += 1;
		} else {
			// ---- structurally empty column: 1x1 zero block (rule 4; no
			// certified question arises, the zero is structural)
			nzero = nzero + _Index(1);
			k += 1;
		}
	}

	if (npos + nneg + nzero != static_cast<_Index>(un)) {
		out.status = inertia_status::internal_error;   // rule 5 self-check
		return out;
	}
	out.n_pos = npos;
	out.n_neg = nneg;
	out.n_zero = nzero;
	out.status = inertia_status::success;
	return out;
}

// ---------------------------------------------------------------------------
// policy_inertia_with_info: NVI outer (finalize + squareness entry checks).
// Must never be overridden -- override policy_inertia_with_info_impl.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
inertia_result<_Index> spmats<_T, _Index>::policy_inertia_with_info(
	const inertia_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_inertia_with_info: matrix must be square");
	return policy_inertia_with_info_impl(opt);
}

// ---------------------------------------------------------------------------
// policy_inertia_with_info_impl: default = "call policy_ldl_with_info, feed
// the D it returns to policy_inertia_from_block_diagonal" (design decision
// 4).  Internal ldl status zero_pivot is NOT a failure (the zero eigenvalues
// appear as structural zeros of D and land in n_zero); every other
// non-success -- inconclusive_pivot_test included -- maps to
// factorization_failed with the internal status kept in the diagnostics.
// Virtual: derived policies may replace the computation.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
inertia_result<_Index> spmats<_T, _Index>::policy_inertia_with_info_impl(
	const inertia_options<_T>& opt) const
{
	try {
		spmats<_T, _Index> L, D;
		std::vector<_Index> perm;
		const ldl_result<_T, _Index> lr =
		    policy_ldl_with_info(L, D, perm, opt.ldl);
		if (lr.status == sparse_ldl_status::success ||
		    lr.status == sparse_ldl_status::zero_pivot) {
			inertia_result<_Index> out =
			    D.policy_inertia_from_block_diagonal(opt.zero_tol);
			out.ldl_status = lr.status;
			return out;
		}
		inertia_result<_Index> out;
		out.status = inertia_status::factorization_failed;
		out.ldl_status = lr.status;
		return out;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		inertia_result<_Index> out;
		out.status = inertia_status::internal_error;
		return out;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_LDL_IMPL_HPP
