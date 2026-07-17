// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_ainv_impl.hpp
// Out-of-line policy method definitions for the AINV approximate inverse
// (AINV-1; design sandbox/docs/design/ainv_design_v1.2.md).  Incomplete
// biconjugation after [BT98] (M. Benzi, M. Tuma, SIAM J. Sci. Comput.
// 19(3), 1998, pp. 968-994), SDS left-looking formulation ([BT98] §7):
//   policy_ainv_with_info            construction (factors only, D-13)
//   policy_ainv_apply                z = Z (D^{-1} (W^T r)) (D-15)
//   policy_ainv_residual_norm_estimate  non-guaranteed ||I - R A||_inf (D-11)
// R = Z D^{-1} W^T is NEVER materialized, and neither is W^T (the W^T
// products are scatter scans over W's stored lines).  This file is included
// inside spmats.hpp AFTER the closing brace of spmats<_T,_Index>, alongside
// spmats_ldl_impl.hpp.  The types live in spmats_base/spmats_ainv.hpp
// (included before the class body).

#ifndef VCP_SPMATS_AINV_IMPL_HPP
#define VCP_SPMATS_AINV_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <vector>

#include <vcp/spmats_base/spmats_ainv.hpp>

namespace vcp {

namespace spmats_ainv_detail {

	// -----------------------------------------------------------------------
	// ainv_biconjugate_side_: one incomplete biconjugation sweep ([BT98] §7,
	// SDS left-looking) over the CSR of M (M = A for the Z side, M = A^T for
	// the W side; both sides run this same function -- directive A2 step 5).
	// Produces the unit upper triangular factor F (CSC arrays, columns in
	// ascending row order, exact zeros not stored) and the lifted per-column
	// pivots (p_i for the Z side, q_i for the W side).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void ainv_biconjugate_side_(
	    const spmats<_T, _Index>& Mcsr,       // finalized CSR, n x n
	    const ainv_options<_T>& opt,
	    std::vector<_Index>& f_col_ptr,
	    std::vector<_Index>& f_row_ind,
	    std::vector<_T>& f_val,
	    std::vector<_T>& pivots,
	    _Index& n_modifications,
	    _Index& first_modified)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		const _Index n = Mcsr.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);
		const std::vector<_Index>& mp = Mcsr.outer_index();
		const std::vector<_Index>& mi = Mcsr.inner_index();
		const std::vector<_T>&     mv = Mcsr.values();

		f_col_ptr.assign(un + 1u, _Index(0));
		f_row_ind.clear();
		f_val.clear();
		pivots.assign(un, _T(0));
		n_modifications = _Index(0);
		first_modified = _Index(-1);

		// dense work column + occupied index list (one n-length work vector,
		// only the occupied positions are reset per column -- design §3.2)
		std::vector<_T> work(un, _T(0));
		std::vector<char> inwork(un, char(0));
		std::vector<_Index> occ;    // insertion order; occ[0] is the diagonal
		std::vector<_Index> keep;

		for (_Index i = _Index(0); i < n; ++i) {
			const std::size_t si = static_cast<std::size_t>(i);
			occ.clear();
			work[si] = _T(1);       // explicit unit diagonal (untouchable, §3.3)
			inwork[si] = char(1);
			occ.push_back(i);

			// left-looking updates against the final columns j < i (naive
			// leading-column sweep, design OPEN-10; SDS re-computes the
			// inner products instead of keeping DDS update lists)
			for (_Index j = _Index(0); j < i; ++j) {
				const std::size_t sj = static_cast<std::size_t>(j);
				// p = m_j^T z_i: sparse dot of M row j with the work column
				_T p = _T(0);
				for (_Index q = mp[sj]; q < mp[sj + 1u]; ++q) {
					const std::size_t sq = static_cast<std::size_t>(q);
					const std::size_t k = static_cast<std::size_t>(mi[sq]);
					if (inwork[k] != char(0)) p += mv[sq] * work[k];
				}
				// A coefficient-zero update contributes only exact-zero
				// fill-in (dropped by |0| < tau, or never representable in
				// the output: explicit zeros are filtered at store time), so
				// skipping it leaves every output value unchanged.
				if (p == _T(0)) continue;
				const _T alpha = p / pivots[sj];   // pivot lifted: division defined
				// z_i <- z_i - alpha * f_j, drop applied per update
				// ([BT98] §7 / design §3.3, D-16): only NEW fill-in is
				// tested against tau; existing occupied positions (unit
				// diagonal included) accumulate unconditionally and are
				// never dropped, even when the result becomes small
				// (OPEN-23 reading of the [BT98] wording -- interpretation
				// of the paper's text, not stated as an algorithm there).
				for (_Index q = f_col_ptr[sj]; q < f_col_ptr[sj + 1u]; ++q) {
					const std::size_t sq = static_cast<std::size_t>(q);
					const std::size_t k = static_cast<std::size_t>(f_row_ind[sq]);
					if (inwork[k] != char(0)) {
						work[k] -= alpha * f_val[sq];
					} else {
						const _T cand = -(alpha * f_val[sq]);
						const R mag = vcp::tsparse_scalar::abs_value(cand);
						if (mag < opt.drop_tolerance) continue;   // new fill-in dropped
						inwork[k] = char(1);
						work[k] = cand;
						occ.push_back(static_cast<_Index>(k));
					}
				}
			}

			// column confirmation 1/3: max_nnz_per_column top-magnitude
			// selection, applied once here (design §3.3).  The unit diagonal
			// is not counted and always survives.  Equal-magnitude ties keep
			// the earlier-arrived candidate (stable sort on insertion
			// order).  The dual-threshold selection algorithm and this tie
			// rule are NOT in [BT98] (ILUT-style, [BT98] §8 lists dual
			// thresholds as future work) -- this implementation's decision.
			if (opt.max_nnz_per_column > 0 &&
			    occ.size() - 1u > opt.max_nnz_per_column) {
				keep.assign(occ.begin() + 1, occ.end());
				std::stable_sort(keep.begin(), keep.end(),
					[&work](const _Index a, const _Index b) {
						return vcp::tsparse_scalar::abs_value(
						           work[static_cast<std::size_t>(a)]) >
						       vcp::tsparse_scalar::abs_value(
						           work[static_cast<std::size_t>(b)]);
					});
				for (std::size_t t = opt.max_nnz_per_column; t < keep.size(); ++t) {
					const std::size_t k = static_cast<std::size_t>(keep[t]);
					inwork[k] = char(0);
					work[k] = _T(0);
				}
				occ.resize(1u);
				occ.insert(occ.end(), keep.begin(),
				           keep.begin() + static_cast<std::ptrdiff_t>(opt.max_nnz_per_column));
			}

			// column confirmation 2/3: pivot on the SELECTED column, then
			// the sign-preserving lift (design §4.3).  The lift target
			// magnitude follows [BT98] §7 (eps -> 1e-3, "no theoretical
			// justification" there either), but sign preservation and the
			// exact-zero -> +lift rule are NOT in [BT98] -- this
			// implementation's decision.  The sign is read through
			// real_part (identity for real scalars; keeps the complex
			// instantiation of this virtual-reachable body compilable,
			// mirroring the inertia precedent).
			_T piv = _T(0);
			for (_Index q = mp[si]; q < mp[si + 1u]; ++q) {
				const std::size_t sq = static_cast<std::size_t>(q);
				const std::size_t k = static_cast<std::size_t>(mi[sq]);
				if (inwork[k] != char(0)) piv += mv[sq] * work[k];
			}
			{
				const R mag = vcp::tsparse_scalar::abs_value(piv);
				if (mag <= opt.pivot_small_threshold) {
					const R rpiv = vcp::tsparse_scalar::real_part(piv);
					piv = (rpiv < R(0)) ? _T(-opt.pivot_lift_value)
					                    : _T(opt.pivot_lift_value);
					n_modifications = n_modifications + _Index(1);
					if (first_modified < _Index(0)) first_modified = i;
				}
			}
			pivots[si] = piv;

			// column confirmation 3/3: store the column in ascending row
			// order (assign_csc contract) and reset the work column.
			// A surviving position whose value is exactly zero
			// (cancellation) has no stored representation -- the spmats
			// invariant forbids explicit zeros, and an unstored position is
			// the same implicit zero (D-16's "never drop an existing
			// component" is the tau test, which this does not touch).
			std::sort(occ.begin(), occ.end());
			for (std::size_t t = 0; t < occ.size(); ++t) {
				const std::size_t k = static_cast<std::size_t>(occ[t]);
				if (!(work[k] == _T(0))) {
					f_row_ind.push_back(occ[t]);
					f_val.push_back(work[k]);
				}
				inwork[k] = char(0);
				work[k] = _T(0);
			}
			f_col_ptr[si + 1u] = static_cast<_Index>(f_val.size());
		}
	}

	// -----------------------------------------------------------------------
	// ainv_validate_factors_: shared input validation of the factor triple
	// (single implementation used by apply AND the estimate -- directive
	// A3-2 step 1).  Z / W / D all n x n; D strictly diagonal with a stored
	// (hence nonzero, by the spmats invariant) entry in every diagonal
	// position: an off-diagonal entry, or a structurally zero diagonal, is
	// invalid_input (design §2.3).  A finalized CSR and CSC store a diagonal
	// matrix identically, so outer line o must hold exactly the entry (o,o).
	// On success d receives the diagonal values.  Callers must have
	// finalized the arguments.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline ainv_status ainv_validate_factors_(
	    const spmats<_T, _Index>& Z, const spmats<_T, _Index>& W,
	    const spmats<_T, _Index>& D, std::vector<_T>& d)
	{
		const _Index n = Z.rowsize();
		if (Z.columnsize() != n || W.rowsize() != n || W.columnsize() != n ||
		    D.rowsize() != n || D.columnsize() != n) {
			return ainv_status::invalid_input;
		}
		const std::vector<_Index>& dp = D.outer_index();
		const std::vector<_Index>& di = D.inner_index();
		const std::vector<_T>& dv = D.values();
		const std::size_t un = static_cast<std::size_t>(n);
		d.assign(un, _T(0));
		for (std::size_t o = 0; o < un; ++o) {
			const _Index first = dp[o];
			const _Index last = dp[o + 1u];
			if (last - first != _Index(1)) {
				// 0 entries: structurally zero diagonal; >1: off-diagonal
				return ainv_status::invalid_input;
			}
			if (di[static_cast<std::size_t>(first)] != static_cast<_Index>(o)) {
				return ainv_status::invalid_input;   // off-diagonal entry
			}
			d[o] = dv[static_cast<std::size_t>(first)];
		}
		return ainv_status::success;
	}

	// -----------------------------------------------------------------------
	// ainv_scatter_outer_lines_: the ONE W^T-product / row-product scatter
	// kernel shared by apply and the estimate (directive A3-2 step 4: the
	// implementation is not written twice).  For each ACTIVE outer line o of
	// the given finalized storage, accumulate
	//     acc[inner[q]] += values[q] * coeff[o].
	//   apply, stage 1:  arrays = CSR of W, active = all rows k, coeff = r
	//                    -> acc = W^T r  (design §2.4 scatter; W^T never
	//                    materialized)
	//   estimate:        arrays = CSC of W (whose lines ARE the rows of
	//                    W^T), active = supp(v), coeff = v -> acc = W v;
	//                    then arrays = CSR of A, active = supp(r_i),
	//                    coeff = r_i -> acc = A^T r_i = (r_i^T A)^T.
	// Newly touched positions are appended to occ / marked in inocc when
	// track_occ is true (the estimate's O(n)-reset bookkeeping).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void ainv_scatter_outer_lines_(
	    const std::vector<_Index>& outer, const std::vector<_Index>& inner,
	    const std::vector<_T>& val,
	    const std::vector<_Index>& active, const std::vector<_T>& coeff,
	    std::vector<_T>& acc,
	    const bool track_occ, std::vector<char>& inocc, std::vector<_Index>& occ)
	{
		for (std::size_t t = 0; t < active.size(); ++t) {
			const std::size_t o = static_cast<std::size_t>(active[t]);
			const _T c = coeff[o];
			for (_Index q = outer[o]; q < outer[o + 1u]; ++q) {
				const std::size_t sq = static_cast<std::size_t>(q);
				const std::size_t k = static_cast<std::size_t>(inner[sq]);
				acc[k] += val[sq] * c;
				if (track_occ && inocc[k] == char(0)) {
					inocc[k] = char(1);
					occ.push_back(static_cast<_Index>(k));
				}
			}
		}
	}

	// -----------------------------------------------------------------------
	// ainv_scale_diag_inv_: the ONE D^{-1} componentwise scaling shared by
	// apply (active = all indices) and the estimate (active = supp(u)).
	// d holds validated, nonzero diagonal values.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void ainv_scale_diag_inv_(
	    const std::vector<_T>& d, const std::vector<_Index>& active,
	    std::vector<_T>& acc)
	{
		for (std::size_t t = 0; t < active.size(); ++t) {
			const std::size_t k = static_cast<std::size_t>(active[t]);
			acc[k] /= d[k];
		}
	}

} // namespace spmats_ainv_detail

// ---------------------------------------------------------------------------
// policy_ainv_with_info: NVI outer (non-virtual).  Finalize guarantee (same
// auto-finalize as the existing NVI outers) + squareness entry check --
// NON-throwing: a non-square input is invalid_input (design §1.3), unlike
// the throwing ldl / lu outers.  Must never be overridden -- override
// policy_ainv_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ainv_result<_T, _Index> spmats<_T, _Index>::policy_ainv_with_info(
	spmats<_T, _Index>& Z,
	spmats<_T, _Index>& W,
	spmats<_T, _Index>& D,
	const ainv_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize()) {
		Z.resize(_Index(0), _Index(0));
		W.resize(_Index(0), _Index(0));
		D.resize(_Index(0), _Index(0));
		ainv_result<_T, _Index> out;
		out.status = ainv_status::invalid_input;
		return out;
	}
	return policy_ainv_with_info_impl(Z, W, D, opt);
}

// ---------------------------------------------------------------------------
// policy_ainv_with_info_impl: virtual algorithm body (default: [BT98] §7
// SDS incomplete biconjugation, one shared sweep per side).  Numerical
// events (tiny pivots) are lifted silently (D-4) and never a failure
// status; a singular input completes with success and is exposed by the
// downstream gate (D-3).  Misuse (vcp::error) keeps its throwing contract;
// the final std::exception net maps to internal_error (P3 pattern).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ainv_result<_T, _Index> spmats<_T, _Index>::policy_ainv_with_info_impl(
	spmats<_T, _Index>& Z,
	spmats<_T, _Index>& W,
	spmats<_T, _Index>& D,
	const ainv_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	ainv_result<_T, _Index> out;
	try {
		const _Index n = A.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		// a_i^T = CSR row i of A; c_i^T = row i of A^T (= column i of A).
		// transpose() returns a finalized CSR copy (spmats.hpp), built once
		// per construction (design §3.2).
		const spmats<_T, _Index> Ac = A.as_csr();
		const spmats<_T, _Index> At = A.transpose();

		std::vector<_Index> zp, zi, wp, wi;
		std::vector<_T> zv, wv, pivz, pivw;
		_Index modz = _Index(0), modw = _Index(0);
		_Index firstz = _Index(-1), firstw = _Index(-1);

		spmats_ainv_detail::ainv_biconjugate_side_<_T, _Index>(
		    Ac, opt, zp, zi, zv, pivz, modz, firstz);
		spmats_ainv_detail::ainv_biconjugate_side_<_T, _Index>(
		    At, opt, wp, wi, wv, pivw, modw, firstw);

		// born-finalized outputs (ldl / lu W-1/W-2 precedent)
		Z.assign_csc(n, n, zp, zi, zv);
		W.assign_csc(n, n, wp, wi, wv);

		// D(i,i) = Z-side lifted pivot p_i (design §4.5: the W-side q_i is
		// used only inside the W sweep and is NOT written to D -- storage
		// rule not in [BT98], this implementation's decision).  An
		// exact-zero pivot (possible only with a user-set zero lift value)
		// has no stored representation (spmats invariant); such a D is then
		// rejected by the apply / estimate validation as a zero diagonal.
		{
			std::vector<_Index> dp(un + 1u, _Index(0)), di;
			std::vector<_T> dv;
			for (std::size_t k = 0; k < un; ++k) {
				if (!(pivz[k] == _T(0))) {
					di.push_back(static_cast<_Index>(k));
					dv.push_back(pivz[k]);
				}
				dp[k + 1u] = static_cast<_Index>(di.size());
			}
			D.assign_csc(n, n, dp, di, dv);
		}

		out.status = ainv_status::success;
		out.n_pivot_modifications = modz + modw;   // Z-side p + W-side q (§4.5)
		// first (lowest) column index at which any lift occurred, over both
		// sides (-1 = none); advisory only (D-4)
		if (firstz >= _Index(0) && firstw >= _Index(0)) {
			out.first_modified_pivot = (firstz < firstw) ? firstz : firstw;
		} else {
			out.first_modified_pivot = (firstz >= _Index(0)) ? firstz : firstw;
		}
		out.nnz_Z = static_cast<_Index>(zv.size());
		out.nnz_W = static_cast<_Index>(wv.size());
		return out;
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		out = ainv_result<_T, _Index>();
		out.status = ainv_status::internal_error;
		return out;
	}
}

// ---------------------------------------------------------------------------
// policy_ainv_apply: NVI outer (non-virtual, consumer type --
// policy_lu_solve_with_info precedent).  Finalize guarantee on the factor
// arguments + ALL dimension / structure checks (invalid_input,
// non-throwing).  Must never be overridden -- override
// policy_ainv_apply_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ainv_status spmats<_T, _Index>::policy_ainv_apply(
	const spmats<_T, _Index>& Z,
	const spmats<_T, _Index>& W,
	const spmats<_T, _Index>& D,
	const std::vector<_T>& r,
	std::vector<_T>& z) const
{
	if (!Z.is_finalized()) Z.finalize();
	if (!W.is_finalized()) W.finalize();
	if (!D.is_finalized()) D.finalize();
	std::vector<_T> d;
	if (spmats_ainv_detail::ainv_validate_factors_<_T, _Index>(Z, W, D, d)
	        != ainv_status::success ||
	    r.size() != static_cast<std::size_t>(Z.rowsize())) {
		z.clear();
		return ainv_status::invalid_input;
	}
	return policy_ainv_apply_impl(Z, W, D, r, z);
}

// ---------------------------------------------------------------------------
// policy_ainv_apply_impl: virtual algorithm body (default).  Three sparse
// stages, [BT98] eq. (6) structure x = Z D^{-1} W^T b:
//   t = W^T r    scatter over the CSR rows of W (W^T never materialized)
//   t <- D^{-1} t componentwise (shared helper)
//   z = Z t      ordinary CSR matvec
// One apply is O(nnz(Z) + nnz(W)) plus the O(n) scaling ([BT98] §1 design
// goal: an application costs about one matvec with A).  z is a valid output
// only on success (empty otherwise).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ainv_status spmats<_T, _Index>::policy_ainv_apply_impl(
	const spmats<_T, _Index>& Z,
	const spmats<_T, _Index>& W,
	const spmats<_T, _Index>& D,
	const std::vector<_T>& r,
	std::vector<_T>& z) const
{
	try {
		const _Index n = Z.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		// diagonal extraction through the shared validation helper (the
		// outer already accepted the triple; a direct _impl caller gets the
		// same gate)
		std::vector<_T> d;
		if (spmats_ainv_detail::ainv_validate_factors_<_T, _Index>(Z, W, D, d)
		        != ainv_status::success ||
		    r.size() != un) {
			z.clear();
			return ainv_status::invalid_input;
		}

		// stage 1: t = W^T r -- the shared scatter kernel over W's CSR rows
		const spmats<_T, _Index> Wr = W.as_csr();
		std::vector<_T> t(un, _T(0));
		std::vector<_Index> all(un);
		for (std::size_t k = 0; k < un; ++k) all[k] = static_cast<_Index>(k);
		std::vector<char> no_marks;
		std::vector<_Index> no_occ;
		spmats_ainv_detail::ainv_scatter_outer_lines_<_T, _Index>(
		    Wr.outer_index(), Wr.inner_index(), Wr.values(),
		    all, r, t, false, no_marks, no_occ);

		// stage 2: t <- D^{-1} t (shared componentwise scaling)
		spmats_ainv_detail::ainv_scale_diag_inv_<_T, _Index>(d, all, t);

		// stage 3: z = Z t (ordinary CSR matvec)
		const spmats<_T, _Index> Zr = Z.as_csr();
		z.assign(un, _T(0));
		if (n > _Index(0)) Zr.mul_vec(t.data(), z.data());
		return ainv_status::success;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		z.clear();
		return ainv_status::internal_error;
	}
}

// ---------------------------------------------------------------------------
// policy_ainv_residual_norm_estimate: NON-GUARANTEED estimate of
// ||I - Z D^{-1} W^T A||_inf, computed in plain T point arithmetic (D-11).
// No statement about the true norm is implied; accepting or rejecting the
// factors on the basis of this value is entirely the downstream consumer's
// responsibility (D-1).  Row by row and factor-based -- R is never
// materialized (design §2.3): for each row i,
//   u = Z^T e_i (CSR row i of Z), v = D^{-1} u (shared scaling),
//   r_i = W v (shared scatter over W's CSC lines = rows of W^T),
//   s_i^T = r_i^T A (shared scatter over A's CSR rows in supp(r_i)),
//   est_i = sum_j |s_ij - delta_ij|;  est = max_i est_i.
// All rows are computed, no early exit (design OPEN-18).  Working memory
// O(n) (dense work vectors + occupied lists, reset per row).  Single
// non-virtual method (design §2.3); consumer-type checks inline,
// non-throwing.  (a)-type: A is *this.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ainv_status spmats<_T, _Index>::policy_ainv_residual_norm_estimate(
	const spmats<_T, _Index>& Z,
	const spmats<_T, _Index>& W,
	const spmats<_T, _Index>& D,
	typename vcp::tsparse_scalar::real_type<_T>::type& est) const
{
	typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	try {
		if (!A.is_finalized()) A.finalize();
		if (!Z.is_finalized()) Z.finalize();
		if (!W.is_finalized()) W.finalize();
		if (!D.is_finalized()) D.finalize();
		std::vector<_T> d;
		if (spmats_ainv_detail::ainv_validate_factors_<_T, _Index>(Z, W, D, d)
		        != ainv_status::success ||
		    A.rowsize() != Z.rowsize() || A.columnsize() != Z.rowsize()) {
			return ainv_status::invalid_input;
		}

		const _Index n = Z.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		const spmats<_T, _Index> Zr = Z.as_csr();   // Z^T e_i = CSR row i of Z
		const spmats<_T, _Index> Wc = W.as_csc();   // CSC lines of W = rows of W^T
		const spmats<_T, _Index> Ac = A.as_csr();

		const std::vector<_Index>& zp = Zr.outer_index();
		const std::vector<_Index>& zi = Zr.inner_index();
		const std::vector<_T>&     zv = Zr.values();

		std::vector<_T> workv(un, _T(0)), workr(un, _T(0)), works(un, _T(0));
		std::vector<char> markv(un, char(0)), markr(un, char(0)), marks(un, char(0));
		std::vector<_Index> occv, occr, occs;

		R best = R(0);
		for (_Index i = _Index(0); i < n; ++i) {
			const std::size_t si = static_cast<std::size_t>(i);
			occv.clear(); occr.clear(); occs.clear();

			// v = D^{-1} (Z^T e_i): read CSR row i of Z (no copy), divide
			// componentwise through the shared scaling helper
			for (_Index q = zp[si]; q < zp[si + 1u]; ++q) {
				const std::size_t sq = static_cast<std::size_t>(q);
				const std::size_t j = static_cast<std::size_t>(zi[sq]);
				workv[j] = zv[sq];
				markv[j] = char(1);
				occv.push_back(zi[sq]);
			}
			spmats_ainv_detail::ainv_scale_diag_inv_<_T, _Index>(d, occv, workv);

			// r_i = W v: the SAME scatter kernel as apply stage 1, fed with
			// W's CSC lines (column j of W = row j of W^T)
			spmats_ainv_detail::ainv_scatter_outer_lines_<_T, _Index>(
			    Wc.outer_index(), Wc.inner_index(), Wc.values(),
			    occv, workv, workr, true, markr, occr);

			// s_i^T = r_i^T A: scatter over A's CSR rows k in supp(r_i)
			spmats_ainv_detail::ainv_scatter_outer_lines_<_T, _Index>(
			    Ac.outer_index(), Ac.inner_index(), Ac.values(),
			    occr, workr, works, true, marks, occs);

			// est_i = sum_j |s_ij - delta_ij|: untouched positions are
			// zero, so only j = i contributes (|0 - 1| = 1) when it was
			// never touched
			R esti = R(0);
			bool diag_seen = false;
			for (std::size_t t = 0; t < occs.size(); ++t) {
				const std::size_t j = static_cast<std::size_t>(occs[t]);
				if (j == si) {
					esti += vcp::tsparse_scalar::abs_value(works[j] - _T(1));
					diag_seen = true;
				} else {
					esti += vcp::tsparse_scalar::abs_value(works[j]);
				}
			}
			if (!diag_seen) esti += R(1);
			if (esti > best) best = esti;

			// O(n) reset through the occupied lists
			for (std::size_t t = 0; t < occv.size(); ++t) {
				const std::size_t k = static_cast<std::size_t>(occv[t]);
				workv[k] = _T(0); markv[k] = char(0);
			}
			for (std::size_t t = 0; t < occr.size(); ++t) {
				const std::size_t k = static_cast<std::size_t>(occr[t]);
				workr[k] = _T(0); markr[k] = char(0);
			}
			for (std::size_t t = 0; t < occs.size(); ++t) {
				const std::size_t k = static_cast<std::size_t>(occs[t]);
				works[k] = _T(0); marks[k] = char(0);
			}
		}
		est = best;
		return ainv_status::success;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		return ainv_status::internal_error;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_AINV_IMPL_HPP
