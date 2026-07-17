// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_chol_impl.hpp
// Out-of-line policy method definitions for the sparse LL^T (Cholesky)
// factorization (CHOL-2; chol design v1 SS3/SS5).  This file is included
// inside spmats.hpp AFTER the closing brace of spmats<_T,_Index>, alongside
// spmats_ldl_impl.hpp.  The types live in spmats_base/spmats_chol.hpp
// (included before the class body).

#ifndef VCP_SPMATS_CHOL_IMPL_HPP
#define VCP_SPMATS_CHOL_IMPL_HPP

#include <exception>
#include <type_traits>
#include <vector>

#include <vcp/spmats_base/spmats_chol.hpp>

namespace vcp {

namespace spmats_chol_detail {

	// -----------------------------------------------------------------------
	// dispatch_sparse_chol_: SFINAE-guarded helper, same pattern as
	// spmats_ldl_detail::dispatch_sparse_ldl_ -- the signed-Index body is
	// never instantiated for unsigned Index, and the unsigned path reports
	// exactly like the LDL dispatch (throws vcp::state_error).
	// -----------------------------------------------------------------------

	// copy status + integer diagnostics from the tsparse result
	template <typename _T, typename _Index>
	inline chol_result<_T, _Index> make_chol_result_(
	    const sparse_chol_result<_T, _Index>& r)
	{
		chol_result<_T, _Index> out;
		out.status              = r.status;
		out.failure_at          = r.failure_at;
		out.inconclusive_at     = r.inconclusive_at;
		out.structural_empty_at = r.structural_empty_at;
		out.nnz_L               = r.nnz_L;
		out.ordering_used       = r.ordering_used;
		out.method_used         = r.method_used;
		return out;
	}

	// signed Index path: spmats -> CSC -> sparse_chol_factorize_with_info ->
	// L spmats construction + perm.
	template <typename _T, typename _Index>
	inline typename std::enable_if<std::is_signed<_Index>::value,
	                               chol_result<_T, _Index> >::type
	dispatch_sparse_chol_(
	    const spmats<_T, _Index>& A,
	    spmats<_T, _Index>& L,
	    std::vector<_Index>& perm,
	    const chol_options<_T>& opt)
	{
		const _Index n = A.rowsize();

		// spmats -> CSC via the SLU conversion helper (include-only reuse;
		// same pattern as the sparse_ldl dispatch path).
		const csc_storage<_T, _Index> C = sparse_lu_make_csc_storage(A);

		const sparse_chol_result<_T, _Index> r =
		    sparse_chol_factorize_with_info<_T, _Index>(
		        n, C.col_ptr, C.row_ind, C.values, opt);

		chol_result<_T, _Index> out = make_chol_result_<_T, _Index>(r);

		// L / perm are valid outputs only on success (design D-3: chol has
		// no LDL-style "completed failure"); otherwise the out parameters
		// are left empty.
		L.resize(_Index(0), _Index(0));
		perm.clear();
		if (r.status != sparse_chol_status::success) {
			return out;
		}

		// L: the kernel CSC keeps numerically cancelled zeros (stored count
		// == nnz_L, design SS4.3); the spmats invariant forbids explicit
		// zeros, so entries CERTIFIED equal to zero (e == T(0), the same
		// rule as the LDL boundary layer) are dropped here.  Column order
		// and ascending rows are preserved, so assign_csc accepts the
		// filtered triple directly.  Materialized stored count <= nnz_L.
		{
			const std::size_t un = static_cast<std::size_t>(n);
			std::vector<_Index> l_ptr(un + 1u, _Index(0));
			std::vector<_Index> l_ind;
			std::vector<_T>     l_val;
			l_ind.reserve(r.L_row_ind.size());
			l_val.reserve(r.L_val.size());
			for (std::size_t j = 0; j < un; ++j) {
				for (_Index k = r.L_col_ptr[j]; k < r.L_col_ptr[j + 1u]; ++k) {
					const _T& e = r.L_val[static_cast<std::size_t>(k)];
					if (e == _T(0)) { /* certified zero: not stored */ }
					else {
						l_ind.push_back(r.L_row_ind[static_cast<std::size_t>(k)]);
						l_val.push_back(e);
					}
				}
				l_ptr[j + 1u] = static_cast<_Index>(l_ind.size());
			}
			L.assign_csc(n, n, l_ptr, l_ind, l_val);
		}

		perm = r.perm;
		return out;
	}

	// unsigned Index path: sparse CHOL cannot be used (Index must be
	// signed); same reporting convention as dispatch_sparse_ldl_.
	template <typename _T, typename _Index>
	inline typename std::enable_if<!std::is_signed<_Index>::value,
	                               chol_result<_T, _Index> >::type
	dispatch_sparse_chol_(
	    const spmats<_T, _Index>& A,
	    spmats<_T, _Index>& L,
	    std::vector<_Index>& perm,
	    const chol_options<_T>& opt)
	{
		(void)A; (void)L; (void)perm; (void)opt;
		vcp::throw_error<vcp::state_error>(
		    "spmats::policy_chol_with_info: sparse CHOL requires a signed Index type");
		return chol_result<_T, _Index>();
	}

} // namespace spmats_chol_detail

// ---------------------------------------------------------------------------
// policy_chol_with_info: NVI outer (non-virtual).  Finalize guarantee (same
// auto-finalize as the existing NVI outers) + squareness entry check, then
// delegates to the virtual policy_chol_with_info_impl.  Must never be
// overridden -- override policy_chol_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
chol_result<_T, _Index> spmats<_T, _Index>::policy_chol_with_info(
	spmats<_T, _Index>& L,
	std::vector<_Index>& perm,
	const chol_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // subject is *this (WFIX-2 parity)
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_chol_with_info: matrix must be square");
	return policy_chol_with_info_impl(L, perm, opt);
}

// ---------------------------------------------------------------------------
// policy_chol_with_info_impl: virtual algorithm body (default: signed-Index
// guard -> CSC conversion -> tsparse factorization -> L/perm assembly).
// This is the designated replacement point for external-backend policies
// (spumar / CHOLMOD delegation, design SS7.6).  Runtime failure is a status;
// misuse (vcp::error) keeps its throwing contract; the final std::exception
// net maps to internal_error (P3).  No scalar-type branching exists on this
// path (design SS6).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
chol_result<_T, _Index> spmats<_T, _Index>::policy_chol_with_info_impl(
	spmats<_T, _Index>& L,
	std::vector<_Index>& perm,
	const chol_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // subject is *this (WFIX-2 parity)
	try {
		return spmats_chol_detail::dispatch_sparse_chol_<_T, _Index>(A, L, perm, opt);
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		// P3 final protection net: firing means a certified gate was missed.
		chol_result<_T, _Index> out;
		out.status = sparse_chol_status::internal_error;
		return out;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_CHOL_IMPL_HPP
