// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_lu_extract_impl.hpp
// Out-of-line policy method definitions for the LU factor extraction
// (LUX-1; design lux_design_v0 SS2) and the factor consuming layer
// (LUX-2; design SS2a: policy_lu_solve_with_info /
// policy_lu_inverse_row_with_info).  This file is included inside
// spmats.hpp AFTER the closing brace of spmats<_T,_Index>, alongside
// spmats_ldl_impl.hpp.  The types live in spmats_base/spmats_lu_extract.hpp
// (included before the class body).

#ifndef VCP_SPMATS_LU_EXTRACT_IMPL_HPP
#define VCP_SPMATS_LU_EXTRACT_IMPL_HPP

#include <exception>
#include <type_traits>
#include <vector>

#include <vcp/spmats_base/spmats_lu_extract.hpp>

namespace vcp {

namespace spmats_lu_extract_detail {

	// -----------------------------------------------------------------------
	// dispatch_sparse_lu_extract_: SFINAE-guarded helper, same pattern as
	// spmats_ldl_detail::dispatch_sparse_ldl_ -- the signed-Index body is
	// never instantiated for unsigned Index, and the unsigned path reports
	// exactly like the SLU dispatch (throws vcp::state_error).
	// -----------------------------------------------------------------------

	// signed Index path: P-5 equilibration entry check ->
	// sparse_lu_factorize_with_info -> sparse_lu_extract_factors ->
	// L / U spmats construction (assign_csc: extraction output is
	// column-sorted with exact zeros already dropped) + p / q.
	template <typename _T, typename _Index>
	inline typename std::enable_if<std::is_signed<_Index>::value,
	                               lu_extract_result<_T, _Index> >::type
	dispatch_sparse_lu_extract_(
	    const spmats<_T, _Index>& A,
	    spmats<_T, _Index>& L,
	    spmats<_T, _Index>& U,
	    std::vector<_Index>& p,
	    std::vector<_Index>& q,
	    const lu_extract_options<_T>& opt)
	{
		lu_extract_result<_T, _Index> out;
		out.method_used   = opt.slu.method;     // echo; overwritten below
		out.ordering_used = opt.slu.ordering;   // request record (see type doc)

		// out parameters are valid only on success (LDL contract style)
		L.resize(_Index(0), _Index(0));
		U.resize(_Index(0), _Index(0));
		p.clear();
		q.clear();

		// P-5: the SSC convention has no scaling; reject equilibration at
		// the entry WITHOUT running the factorization.
		if (opt.slu.equilibration) {
			out.status = sparse_lu_extract_status::unsupported_options;
			return out;
		}

		const _Index n = A.rowsize();

		const sparse_lu_factorization<_T, _Index> fac =
		    sparse_lu_factorize_with_info(A, opt.slu);
		out.method_used = fac.info().method_used;

		sparse_lu_extracted<_T, _Index> ex = sparse_lu_extract_factors(fac);
		out.status = ex.status;
		if (ex.status != sparse_lu_extract_status::success) {
			return out;
		}

		L.assign_csc(n, n, ex.L_col_ptr, ex.L_row_ind, ex.L_val);
		U.assign_csc(n, n, ex.U_col_ptr, ex.U_row_ind, ex.U_val);
		p.swap(ex.p);
		q.swap(ex.q);
		out.nnz_L = ex.nnz_L;
		out.nnz_U = ex.nnz_U;
		return out;
	}

	// unsigned Index path: sparse LU cannot be used (Index must be signed);
	// same reporting convention as dispatch_sparse_lu_ / dispatch_sparse_ldl_.
	template <typename _T, typename _Index>
	inline typename std::enable_if<!std::is_signed<_Index>::value,
	                               lu_extract_result<_T, _Index> >::type
	dispatch_sparse_lu_extract_(
	    const spmats<_T, _Index>& A,
	    spmats<_T, _Index>& L,
	    spmats<_T, _Index>& U,
	    std::vector<_Index>& p,
	    std::vector<_Index>& q,
	    const lu_extract_options<_T>& opt)
	{
		(void)A; (void)L; (void)U; (void)p; (void)q; (void)opt;
		vcp::throw_error<vcp::state_error>(
		    "spmats::policy_lu_with_info: sparse LU requires a signed Index type");
		return lu_extract_result<_T, _Index>();
	}

} // namespace spmats_lu_extract_detail

// ---------------------------------------------------------------------------
// policy_lu_with_info: NVI outer (non-virtual).  Finalize guarantee (same
// auto-finalize as the existing NVI outers) + squareness entry check, then
// delegates to the virtual policy_lu_with_info_impl.  Must never be
// overridden -- override policy_lu_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_extract_result<_T, _Index> spmats<_T, _Index>::policy_lu_with_info(
	spmats<_T, _Index>& L,
	spmats<_T, _Index>& U,
	std::vector<_Index>& p,
	std::vector<_Index>& q,
	const lu_extract_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_lu_with_info: matrix must be square");
	return policy_lu_with_info_impl(L, U, p, q, opt);
}

// ---------------------------------------------------------------------------
// policy_lu_with_info_impl: virtual algorithm body (default: signed-Index
// guard -> equilibration entry check -> SLU factorization -> factor
// extraction -> L/U/p/q assembly).  Runtime failure is a status; misuse
// (vcp::error) keeps its throwing contract; the final std::exception net
// maps to internal_error (same P3 pattern as the LDL policy).  This virtual
// is the designated replacement point for external-backend policies
// (spumar / UMFPACK delegation, design SS2).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_extract_result<_T, _Index> spmats<_T, _Index>::policy_lu_with_info_impl(
	spmats<_T, _Index>& L,
	spmats<_T, _Index>& U,
	std::vector<_Index>& p,
	std::vector<_Index>& q,
	const lu_extract_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	try {
		return spmats_lu_extract_detail::dispatch_sparse_lu_extract_<_T, _Index>(
		    A, L, U, p, q, opt);
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		lu_extract_result<_T, _Index> out;
		out.status = sparse_lu_extract_status::internal_error;
		return out;
	}
}

// ===========================================================================
// factor consuming layer (LUX-2, design lux_design_v0 SS2a)
// ===========================================================================

namespace spmats_lu_extract_detail {

	// certified three-branch division gate (P-9; LDL SS3.1 standard form
	// with tol = 0): certified nonzero -> success (divide) / certified
	// zero -> singular_factor / undecidable -> inconclusive_division.
	// For totally ordered scalars the third branch is unreachable.
	template <typename _T>
	inline lu_apply_status certified_division_gate_(const _T& d)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		const R absd = vcp::tsparse_scalar::abs_value(d);
		if (absd > R(0))  return lu_apply_status::success;
		if (absd <= R(0)) return lu_apply_status::singular_factor;
		return lu_apply_status::inconclusive_division;
	}

	// bijection check on [0, n) (index-only)
	template <typename _Index>
	inline bool is_permutation_(const _Index n, const std::vector<_Index>& p)
	{
		if (static_cast<_Index>(p.size()) != n) return false;
		std::vector<char> seen(static_cast<std::size_t>(n), char(0));
		for (std::size_t k = 0; k < p.size(); ++k) {
			if (!(p[k] >= _Index(0)) || !(p[k] < n)) return false;
			const std::size_t s = static_cast<std::size_t>(p[k]);
			if (seen[s]) return false;
			seen[s] = char(1);
		}
		return true;
	}

	// L structural contract (P-9): CSC, every column has an explicit unit
	// diagonal (value exactly T(1)) and no upper (row < col) entry.
	// Column-internal order is not assumed.
	template <typename _T, typename _Index>
	inline bool check_unit_lower_(const spmats<_T, _Index>& Lc, const _Index n)
	{
		const std::vector<_Index>& cp = Lc.outer_index();
		const std::vector<_Index>& ri = Lc.inner_index();
		const std::vector<_T>&     v  = Lc.values();
		for (_Index j = _Index(0); j < n; ++j) {
			const std::size_t sj = static_cast<std::size_t>(j);
			bool diag = false;
			for (_Index k = cp[sj]; k < cp[sj + 1u]; ++k) {
				const std::size_t sk = static_cast<std::size_t>(k);
				const _Index r = ri[sk];
				if (r < j) return false;              // upper entry in L
				if (r == j) {
					if (diag) return false;           // duplicate diagonal
					if (!(v[sk] == _T(1))) return false;   // non-unit diagonal
					diag = true;
				}
			}
			if (!diag) return false;                  // missing unit diagonal
		}
		return true;
	}

	// U structural contract: CSC, no lower (row > col) entry.  Diagonal
	// presence/zero-ness is a NUMERICAL property handled by the certified
	// division gate during the solves (missing diagonal == structural zero
	// -> singular_factor).
	template <typename _T, typename _Index>
	inline bool check_upper_(const spmats<_T, _Index>& Uc, const _Index n)
	{
		const std::vector<_Index>& cp = Uc.outer_index();
		const std::vector<_Index>& ri = Uc.inner_index();
		for (_Index j = _Index(0); j < n; ++j) {
			const std::size_t sj = static_cast<std::size_t>(j);
			for (_Index k = cp[sj]; k < cp[sj + 1u]; ++k) {
				if (ri[static_cast<std::size_t>(k)] > j) return false;
			}
		}
		return true;
	}

} // namespace spmats_lu_extract_detail

// ---------------------------------------------------------------------------
// policy_lu_solve_with_info: NVI outer (non-virtual).  Finalize guarantee on
// the factor arguments + ALL dimension checks (P-7: dimension checking is
// the outer's responsibility; reported as dimension_mismatch, non-throwing).
// Must never be overridden -- override policy_lu_solve_with_info_impl.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_apply_result spmats<_T, _Index>::policy_lu_solve_with_info(
	const spmats<_T, _Index>& L,
	const spmats<_T, _Index>& U,
	const std::vector<_Index>& p,
	const std::vector<_Index>& q,
	const std::vector<_T>& b,
	std::vector<_T>& x) const
{
	if (!L.is_finalized()) L.finalize();
	if (!U.is_finalized()) U.finalize();
	const _Index n = L.rowsize();
	if (L.columnsize() != n || U.rowsize() != n || U.columnsize() != n ||
	    static_cast<_Index>(p.size()) != n ||
	    static_cast<_Index>(q.size()) != n ||
	    static_cast<_Index>(b.size()) != n) {
		x.clear();
		lu_apply_result out;
		out.status = lu_apply_status::dimension_mismatch;
		return out;
	}
	return policy_lu_solve_with_info_impl(L, U, p, q, b, x);
}

// ---------------------------------------------------------------------------
// policy_lu_solve_with_info_impl: virtual algorithm body (default).
// P-8:  x = Q U^{-1} L^{-1} P b  under the SSC convention P A Q = L U:
//   y[k] = b[p[k]]  ->  L forward (explicit unit diagonal, skipped)
//   ->  U backward (certified division gate per column, P-9)
//   ->  x[q[k]] = z[k].
// Consumes ONLY the argument factors (P-7: no policy state is read), so any
// derived policy -- spumar included -- inherits this implementation
// unchanged; a derived policy may override it to change the computation.
// x is a valid output only when the returned status is success (empty
// otherwise).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_apply_result spmats<_T, _Index>::policy_lu_solve_with_info_impl(
	const spmats<_T, _Index>& L,
	const spmats<_T, _Index>& U,
	const std::vector<_Index>& p,
	const std::vector<_Index>& q,
	const std::vector<_T>& b,
	std::vector<_T>& x) const
{
	try {
		lu_apply_result out;
		x.clear();

		const _Index n = L.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		// structural contracts (P-9; index-level + the L unit-diagonal read)
		const spmats<_T, _Index> Lc = L.as_csc();
		const spmats<_T, _Index> Uc = U.as_csc();
		if (!spmats_lu_extract_detail::is_permutation_<_Index>(n, p) ||
		    !spmats_lu_extract_detail::is_permutation_<_Index>(n, q) ||
		    !spmats_lu_extract_detail::check_unit_lower_<_T, _Index>(Lc, n) ||
		    !spmats_lu_extract_detail::check_upper_<_T, _Index>(Uc, n)) {
			out.status = lu_apply_status::invalid_input;
			return out;
		}

		if (n == _Index(0)) {
			out.status = lu_apply_status::success;
			return out;
		}

		// entry: y = P b  ((Pb)[k] = b[p[k]])
		std::vector<_T> y(un);
		for (std::size_t k = 0; k < un; ++k) {
			y[k] = b[static_cast<std::size_t>(p[k])];
		}

		const std::vector<_Index>& lp = Lc.outer_index();
		const std::vector<_Index>& li = Lc.inner_index();
		const std::vector<_T>&     lv = Lc.values();

		// forward: L y' = y (unit diagonal contract: stored diagonal skipped)
		for (_Index j = _Index(0); j < n; ++j) {
			const std::size_t sj = static_cast<std::size_t>(j);
			for (_Index k = lp[sj]; k < lp[sj + 1u]; ++k) {
				const std::size_t sk = static_cast<std::size_t>(k);
				const _Index r = li[sk];
				if (r == j) continue;
				y[static_cast<std::size_t>(r)] -= lv[sk] * y[sj];
			}
		}

		const std::vector<_Index>& up = Uc.outer_index();
		const std::vector<_Index>& ui = Uc.inner_index();
		const std::vector<_T>&     uv = Uc.values();

		// backward: U z = y' (explicit diagonal, certified gate per column)
		for (_Index jj = n; jj-- > _Index(0); ) {
			const std::size_t sj = static_cast<std::size_t>(jj);
			bool found = false;
			_T d = _T(0);
			for (_Index k = up[sj]; k < up[sj + 1u]; ++k) {
				const std::size_t sk = static_cast<std::size_t>(k);
				if (ui[sk] == jj) { d = uv[sk]; found = true; }
			}
			if (!found) {
				out.status = lu_apply_status::singular_factor;
				return out;
			}
			const lu_apply_status g =
			    spmats_lu_extract_detail::certified_division_gate_<_T>(d);
			if (g != lu_apply_status::success) {
				out.status = g;
				return out;
			}
			y[sj] /= d;
			for (_Index k = up[sj]; k < up[sj + 1u]; ++k) {
				const std::size_t sk = static_cast<std::size_t>(k);
				const _Index r = ui[sk];
				if (r == jj) continue;
				y[static_cast<std::size_t>(r)] -= uv[sk] * y[sj];
			}
		}

		// exit: x[q[k]] = z[k]
		x.assign(un, _T(0));
		for (std::size_t k = 0; k < un; ++k) {
			x[static_cast<std::size_t>(q[k])] = y[k];
		}
		out.status = lu_apply_status::success;
		return out;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		x.clear();
		lu_apply_result out;
		out.status = lu_apply_status::internal_error;
		return out;
	}
}

// ---------------------------------------------------------------------------
// policy_lu_inverse_row_with_info: NVI outer (dimension + index-range
// checks; P-7).  Must never be overridden -- override the _impl.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_apply_result spmats<_T, _Index>::policy_lu_inverse_row_with_info(
	const spmats<_T, _Index>& L,
	const spmats<_T, _Index>& U,
	const std::vector<_Index>& p,
	const std::vector<_Index>& q,
	const _Index i,
	std::vector<_T>& row) const
{
	if (!L.is_finalized()) L.finalize();
	if (!U.is_finalized()) U.finalize();
	const _Index n = L.rowsize();
	if (L.columnsize() != n || U.rowsize() != n || U.columnsize() != n ||
	    static_cast<_Index>(p.size()) != n ||
	    static_cast<_Index>(q.size()) != n ||
	    !(i >= _Index(0)) || !(i < n)) {
		row.clear();
		lu_apply_result out;
		out.status = lu_apply_status::dimension_mismatch;
		return out;
	}
	return policy_lu_inverse_row_with_info_impl(L, U, p, q, i, row);
}

// ---------------------------------------------------------------------------
// policy_lu_inverse_row_with_info_impl: virtual algorithm body (default).
// P-8:  row_i(A^{-1})^T = (A^{-1})^T e_i = P^T L^{-T} U^{-T} Q^T e_i.
// The transposed triangular solves read the CSC arrays row-wise (column j
// of U is row j of U^T; same for L) -- no transpose is materialized:
//   w = Q^T e_i = e_{k0} with q[k0] = i
//   U^T y = w   (lower triangular, explicit diagonal; certified gate)
//   L^T z = y   (upper triangular, unit diagonal; no division)
//   row[p[k]] = z[k].
// Consumes only the argument factors (P-7).  row is a valid output only on
// success (empty otherwise).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
lu_apply_result spmats<_T, _Index>::policy_lu_inverse_row_with_info_impl(
	const spmats<_T, _Index>& L,
	const spmats<_T, _Index>& U,
	const std::vector<_Index>& p,
	const std::vector<_Index>& q,
	const _Index i,
	std::vector<_T>& row) const
{
	try {
		lu_apply_result out;
		row.clear();

		const _Index n = L.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		const spmats<_T, _Index> Lc = L.as_csc();
		const spmats<_T, _Index> Uc = U.as_csc();
		if (!spmats_lu_extract_detail::is_permutation_<_Index>(n, p) ||
		    !spmats_lu_extract_detail::is_permutation_<_Index>(n, q) ||
		    !spmats_lu_extract_detail::check_unit_lower_<_T, _Index>(Lc, n) ||
		    !spmats_lu_extract_detail::check_upper_<_T, _Index>(Uc, n)) {
			out.status = lu_apply_status::invalid_input;
			return out;
		}

		// seed: w = Q^T e_i = e_{k0},  q[k0] = i
		_Index k0 = _Index(0);
		for (_Index k = _Index(0); k < n; ++k) {
			if (q[static_cast<std::size_t>(k)] == i) { k0 = k; break; }
		}

		const std::vector<_Index>& up = Uc.outer_index();
		const std::vector<_Index>& ui = Uc.inner_index();
		const std::vector<_T>&     uv = Uc.values();

		// stage 1: U^T y = e_{k0}.  Row j of U^T = column j of U:
		//   U(j,j) y[j] + sum_{r<j} U(r,j) y[r] = delta_{j,k0}.
		// y[j] = 0 for j < k0, so the loop starts at k0; every processed
		// column divides through the certified gate.
		std::vector<_T> y(un, _T(0));
		for (_Index j = k0; j < n; ++j) {
			const std::size_t sj = static_cast<std::size_t>(j);
			_T s = (j == k0) ? _T(1) : _T(0);
			bool found = false;
			_T d = _T(0);
			for (_Index k = up[sj]; k < up[sj + 1u]; ++k) {
				const std::size_t sk = static_cast<std::size_t>(k);
				const _Index r = ui[sk];
				if (r == j) { d = uv[sk]; found = true; }
				else        { s -= uv[sk] * y[static_cast<std::size_t>(r)]; }
			}
			if (!found) {
				out.status = lu_apply_status::singular_factor;
				return out;
			}
			const lu_apply_status g =
			    spmats_lu_extract_detail::certified_division_gate_<_T>(d);
			if (g != lu_apply_status::success) {
				out.status = g;
				return out;
			}
			y[sj] = s / d;
		}

		const std::vector<_Index>& lp = Lc.outer_index();
		const std::vector<_Index>& li = Lc.inner_index();
		const std::vector<_T>&     lv = Lc.values();

		// stage 2: L^T z = y (in place).  Row j of L^T = column j of L;
		// unit diagonal, so no division:
		//   z[j] = y[j] - sum_{r>j} L(r,j) z[r],  j = n-1 .. 0.
		for (_Index jj = n; jj-- > _Index(0); ) {
			const std::size_t sj = static_cast<std::size_t>(jj);
			for (_Index k = lp[sj]; k < lp[sj + 1u]; ++k) {
				const std::size_t sk = static_cast<std::size_t>(k);
				const _Index r = li[sk];
				if (r == jj) continue;
				y[sj] -= lv[sk] * y[static_cast<std::size_t>(r)];
			}
		}

		// exit: row[p[k]] = z[k]
		row.assign(un, _T(0));
		for (std::size_t k = 0; k < un; ++k) {
			row[static_cast<std::size_t>(p[k])] = y[k];
		}
		out.status = lu_apply_status::success;
		return out;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		row.clear();
		lu_apply_result out;
		out.status = lu_apply_status::internal_error;
		return out;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_LU_EXTRACT_IMPL_HPP
