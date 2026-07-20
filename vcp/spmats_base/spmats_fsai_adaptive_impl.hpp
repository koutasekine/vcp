// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_fsai_adaptive_impl.hpp
// Out-of-line policy method definitions for the ADAPTIVE FSAI (FSAI-2;
// design sandbox/docs/design/fsai_adaptive_design_v0.md), [JFSG15]
// Algorithm 3: per-row pattern adaptation driven by the Kaporin-oriented
// quality psi_k,i = G_k[i,:] A G_k[i,:]^T (eq. (22)), gradient-based
// candidate selection (eq. (24)), exit test psi_k/psi_0 <= eps (eq. (25)),
// optional in-iteration drop (eq. (26)), terminal post-filtration, output
// as the SAME sqrt-free triple as FSAI-1:
//   U (unit upper triangular, = Ghat^T), Dhat (diagonal), perm (new->old)
// with Dhat_ii := psi_k,i (design SS3.3 -- THIS DESIGN'S derivation; the
// normalization line of the paper's Algorithm 3 is NOT the normative source,
// directive invariant 12) and R := P U Dhat^{-1} U^T P^T NEVER materialized.
//   policy_fsai_adaptive_with_info (diagonal init)     construction
//   policy_fsai_adaptive_with_info (U0, perm0 init)    combined strategy
// Estimate / apply: use the EXISTING policy_fsai_residual_norm_estimate /
// policy_fsai_apply on the returned triple (identical output contract,
// F2-D6 -- no new helpers).
// This file is included inside spmats.hpp AFTER the closing brace of
// spmats<_T,_Index>, alongside spmats_fsai_impl.hpp.  The FSAI-1 files are
// READ-ONLY for this track (directive invariant 10): the detail helpers
// fsai_perm_is_bijection_ / fsai_perm_is_identity_ /
// fsai_build_permuted_copy_ are reused by call; the dense solve and the
// post-filtration could NOT be reused as functions (the FSAI-1 solve
// hard-codes the e_m right-hand side; the FSAI-1 post-filtration is fused
// into fsai_process_row_), so minimal replicas live in the detail namespace
// below (completion report OPEN-A9 / OPEN-A10).
//
// OpenMP contract (design SS4, F2-D4): rows are fully independent; the
// dynamic per-row fill is handled by CAPPED PRE-ALLOCATED SLICES fixed in a
// serial pass 1 (cap_i = |P0_i| + k_iter*s + 1), so every row writes only
// its own slice and per-row diagnostic cells.  The parallel run is
// BIT-IDENTICAL to the serial run by construction: no shared accumulation,
// no atomic, no reduction, no thread-count-dependent branch; diagnostics
// are per-row arrays combined serially in row order after the loop.

#ifndef VCP_SPMATS_FSAI_ADAPTIVE_IMPL_HPP
#define VCP_SPMATS_FSAI_ADAPTIVE_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <type_traits>
#include <vector>

#include <vcp/spmats_base/spmats_fsai_adaptive.hpp>
#include <vcp/spmats_base/spmats_fsai_impl.hpp>   // spmats_fsai_detail helpers (read-only reuse)

namespace vcp {

namespace spmats_fsai_adaptive_detail {

	// -----------------------------------------------------------------------
	// fsai2_dense_solve_: dense in-place factorization + forward/back
	// substitution for A_p[Pbar,Pbar] g = b with a GENERAL right-hand side.
	// Minimal replica of spmats_fsai_detail::fsai_dense_ldl_solve_ (FSAI-1,
	// read-only -- completion report OPEN-A9): same unit-lower elimination
	// WITHOUT interchanges (Doolittle; sqrt-free LDL^T on a symmetric
	// block), same sign-preserving pivot lift (|p| <= threshold -> +-lift,
	// exact zero -> +lift, sign read through real_part).  The FSAI-1 kernel
	// hard-codes the e_m right-hand side (back substitution only), which the
	// adaptive row system A_p[Pbar,Pbar] g = -A_p[Pbar,i] cannot use, hence
	// this replica adds the unit-lower FORWARD substitution.
	// M is column-major m x m and is destroyed; x is the right-hand side on
	// entry and the solution on return.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai2_dense_solve_(
	    std::vector<_T>& M, const std::size_t m,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& threshold,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& lift,
	    std::vector<_T>& x,
	    _Index& nlifts)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		for (std::size_t k = 0; k < m; ++k) {
			_T p = M[k + k * m];
			if (vcp::tsparse_scalar::abs_value(p) <= threshold) {
				const R rp = vcp::tsparse_scalar::real_part(p);
				p = (rp < R(0)) ? _T(-lift) : _T(lift);
				M[k + k * m] = p;
				nlifts = nlifts + _Index(1);
			}
			for (std::size_t r = k + 1u; r < m; ++r) {
				const _T l = M[r + k * m] / p;
				M[r + k * m] = l;
				for (std::size_t c = k + 1u; c < m; ++c) {
					M[r + c * m] -= l * M[k + c * m];
				}
			}
		}
		// forward substitution (unit lower triangle of M)
		for (std::size_t r = 1u; r < m; ++r) {
			_T s = _T(0);
			for (std::size_t c = 0; c < r; ++c) {
				s += M[r + c * m] * x[c];
			}
			x[r] -= s;
		}
		// back substitution (upper triangle of M, pivots on the diagonal)
		x[m - 1u] = x[m - 1u] / M[(m - 1u) + (m - 1u) * m];
		for (std::size_t r = m - 1u; r-- > 0u; ) {
			_T s = _T(0);
			for (std::size_t c = r + 1u; c < m; ++c) {
				s += M[r + c * m] * x[c];
			}
			x[r] = (x[r] - s) / M[r + r * m];
		}
	}

	// -----------------------------------------------------------------------
	// fsai2_psi_lift_: sign-preserving lift of psi (design SS3.3: same
	// regime as an elimination pivot; keeps Dhat = psi nonzero).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai2_psi_lift_(
	    _T& psi,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& threshold,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& lift,
	    _Index& nlifts)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		if (vcp::tsparse_scalar::abs_value(psi) <= threshold) {
			const R rp = vcp::tsparse_scalar::real_part(psi);
			psi = (rp < R(0)) ? _T(-lift) : _T(lift);
			nlifts = nlifts + _Index(1);
		}
	}

	// -----------------------------------------------------------------------
	// fsai2_u0_is_unit_upper_: entry validation of the initial factor U0
	// (design SS3.6: square/dimension checked by the caller; here each CSC
	// column must be upper triangular with an EXPLICIT unit diagonal).
	// Finalized CSC columns are row-ascending and unique, so the diagonal,
	// when present, is the last stored entry of its column.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline bool fsai2_u0_is_unit_upper_(
	    const std::vector<_Index>& up,
	    const std::vector<_Index>& ui,
	    const std::vector<_T>&     uv,
	    const _Index n)
	{
		for (_Index j = _Index(0); j < n; ++j) {
			const std::size_t sj = static_cast<std::size_t>(j);
			const _Index qb = up[sj];
			const _Index qe = up[sj + 1u];
			if (qe <= qb) return false;               // empty column: no diagonal
			for (_Index q = qb; q < qe; ++q) {
				if (ui[static_cast<std::size_t>(q)] > j) return false;   // below diagonal
			}
			const std::size_t last = static_cast<std::size_t>(qe) - 1u;
			if (ui[last] != j) return false;          // diagonal not stored
			if (!(uv[last] == _T(1))) return false;   // diagonal not exactly one
		}
		return true;
	}

	// -----------------------------------------------------------------------
	// fsai2_candidate: gradient candidate (|w_j|, j) -- sorted by magnitude
	// descending with EXACT ties broken by the smaller index (design
	// OPEN-A3, not in [JFSG15]).  The comparator is a total order, so the
	// sorted result is independent of the collection order.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	struct fsai2_candidate {
		typename vcp::tsparse_scalar::real_type<_T>::type aw;
		_Index j;
	};

	template <typename _T, typename _Index>
	struct fsai2_candidate_before {
		bool operator()(const fsai2_candidate<_T, _Index>& a,
		                const fsai2_candidate<_T, _Index>& b) const {
			if (a.aw > b.aw) return true;
			if (b.aw > a.aw) return false;
			return a.j < b.j;
		}
	};

	// -----------------------------------------------------------------------
	// fsai2_row_workspace: per-thread scratch of the row loop (design SS4:
	// thread-local; in the parallel build one instance lives inside each
	// thread's parallel region).  The n-sized arrays (inpat / wmark / w)
	// obey a CLEAN-ON-EXIT invariant: every row resets exactly the cells it
	// touched, so no cross-row state survives and the produced numbers
	// cannot depend on which thread ran which row.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	struct fsai2_row_workspace {
		std::vector<_Index> pat;     // current off-diagonal pattern Pbar (ascending, < i)
		std::vector<_T>     gt;      // g values aligned with pat
		std::vector<_T>     ci;      // A_p(pat, i) values aligned with pat
		std::vector<char>   inpat;   // size n, pattern membership marker
		std::vector<_T>     w;       // size n, dense gradient work
		std::vector<char>   wmark;   // size n, touched marker for w
		std::vector<_Index> touched; // touched cells of w this iteration
		std::vector<fsai2_candidate<_T, _Index> > cand;
		std::vector<_T>     M;       // m x m gather / factorization buffer
		std::vector<_T>     b;       // rhs / solution buffer
		std::vector<_T>     u;       // final unit-diagonal row (post-filter input)
		std::vector<char>   keep;    // post-filtration keep flags
		std::vector<_Index> sel;     // mmax selection scratch
		std::vector<_Index> eloc;    // dropped locals, ascending
		std::vector<_Index> eglob;   // dropped globals, ascending

		void ensure_n(const std::size_t n) {
			if (inpat.size() != n) {
				inpat.assign(n, char(0));
				wmark.assign(n, char(0));
				w.assign(n, _T(0));
			}
		}
	};

	// -----------------------------------------------------------------------
	// fsai2_gather_col_i_: ws.ci[t] = A_p(pat[t], i) by a two-pointer merge
	// of CSC column i with the ascending pattern (absent entries are zero).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai2_gather_col_i_(
	    const std::size_t i,
	    const std::vector<_Index>& ap,
	    const std::vector<_Index>& ai,
	    const std::vector<_T>&     av,
	    const std::vector<_Index>& pat,
	    std::vector<_T>& ci)
	{
		const std::size_t m = pat.size();
		ci.assign(m, _T(0));
		_Index q = ap[i];
		const _Index qe = ap[i + 1u];
		std::size_t t = 0;
		while (q < qe && t < m) {
			const _Index r = ai[static_cast<std::size_t>(q)];
			if (r < pat[t])      { ++q; }
			else if (r > pat[t]) { ++t; }
			else {
				ci[t] = av[static_cast<std::size_t>(q)];
				++q; ++t;
			}
		}
	}

	// -----------------------------------------------------------------------
	// fsai2_process_row_: the complete per-row adaptive computation.  Reads
	// ONLY A_p (CSC + CSR), the optional U0 column and opt; writes ONLY the
	// row-i slice [slice_off[i], slice_off[i] + cap_i) of u_ind/u_val and
	// the row-i cells of the per-row diagnostic arrays -- the parallel
	// bit-identity contract rests on this (design SS4).
	//
	// Loop shape (equivalent reading of directive SSA2 step 2 fixed by the
	// cap formula cap_i = |P0_i| + k_iter*s + 1): the row STATE (pattern +
	// values + psi) starts from the initial pattern WITHOUT a solve
	// (diagonal init: empty Pbar, psi_0 = a_ii; initial-value path: g from
	// the U0 column, psi_0 recomputed by the SS3.2 inner product -- design
	// SS3.6, no D0 needed), then runs up to k_iter rounds of
	//   exit test (eq. 25) -> inner drop (eq. 26) -> gradient (eq. 24) ->
	//   add top-s candidates -> re-solve (SS3.1) -> psi (SS3.2),
	// so every candidate addition is followed by a solve and the total
	// number of additions is bounded by k_iter*s (the cap).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai2_process_row_(
	    const std::size_t i,
	    const std::size_t n,
	    const std::vector<_Index>& ap,    // CSC col_ptr of A_p
	    const std::vector<_Index>& ai,    // CSC row_ind of A_p
	    const std::vector<_T>&     av,    // CSC values of A_p
	    const std::vector<_Index>& rp,    // CSR row_ptr of A_p
	    const std::vector<_Index>& ri,    // CSR col_ind of A_p
	    const std::vector<_T>&     rv,    // CSR values of A_p
	    const bool has_u0,
	    const std::vector<_Index>* u0p,   // U0 CSC arrays (permuted frame), used
	    const std::vector<_Index>* u0i,   //   only when has_u0
	    const std::vector<_T>*     u0v,
	    const fsai_adaptive_options<_T>& opt,
	    fsai2_row_workspace<_T, _Index>& ws,
	    const std::vector<std::size_t>& slice_off,
	    std::vector<_Index>& u_ind,       // slice target (row-i slice only)
	    std::vector<_T>&     u_val,       // slice target (row-i slice only)
	    std::vector<_Index>& row_count,
	    std::vector<_T>&     d_hat,
	    std::vector<_Index>& row_lifts,
	    std::vector<_Index>& row_dropped,
	    std::vector<_Index>& row_restore_skipped,
	    std::vector<_Index>& row_converged,
	    std::vector<_Index>& row_capped)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		ws.ensure_n(n);
		_Index nl = _Index(0);

		// --- initial pattern P0 (design SS3.6): diagonal only (empty Pbar)
		// or the off-diagonal support/values of U0 column i (already in the
		// permuted frame -- the factor triple lives in permuted coordinates)
		ws.pat.clear();
		ws.gt.clear();
		if (has_u0) {
			for (_Index q = (*u0p)[i]; q < (*u0p)[i + 1u]; ++q) {
				const _Index r = (*u0i)[static_cast<std::size_t>(q)];
				if (r != static_cast<_Index>(i)) {
					ws.pat.push_back(r);
					ws.gt.push_back((*u0v)[static_cast<std::size_t>(q)]);
				}
			}
		}
		for (std::size_t t = 0; t < ws.pat.size(); ++t) {
			ws.inpat[static_cast<std::size_t>(ws.pat[t])] = char(1);
		}

		// a_ii from CSC column i (structurally missing diagonal reads as 0)
		_T aii = _T(0);
		for (_Index q = ap[i]; q < ap[i + 1u]; ++q) {
			if (ai[static_cast<std::size_t>(q)] == static_cast<_Index>(i)) {
				aii = av[static_cast<std::size_t>(q)];
				break;
			}
		}

		// psi_0 = a_ii + g^T A_p[Pbar, i] (design SS3.2 -- THIS DESIGN'S
		// two-line derivation from the normal equations, not stated in
		// [JFSG15]; single inner product, no extra solve)
		fsai2_gather_col_i_<_T, _Index>(i, ap, ai, av, ws.pat, ws.ci);
		_T psi = aii;
		for (std::size_t t = 0; t < ws.gt.size(); ++t) {
			psi += ws.gt[t] * ws.ci[t];
		}
		fsai2_psi_lift_<_T, _Index>(psi, opt.pivot_small_threshold,
		                            opt.pivot_lift_value, nl);
		const _T psi0 = psi;
		const R  apsi0 = vcp::tsparse_scalar::abs_value(psi0);

		bool done = false;
		bool converged = false;
		for (std::size_t it = 0; it < opt.max_iterations && !done; ++it) {
			// (c) exit test, eq. (25): |psi_k| <= eps * |psi_0| (real abs)
			if (vcp::tsparse_scalar::abs_value(psi) <=
			    opt.exit_tolerance * apsi0) {
				converged = true;
				done = true;
				break;
			}
			// (d) in-iteration drop, eq. (26) in the squared comparison
			// form |g_j|^2 <= tau^2 * ||g||_2^2 (sqrt-free equivalent for
			// nonnegative operands; ||g||_2 is the norm of the off-diagonal
			// solution vector).  Dropped positions may re-enter as later
			// candidates.
			if (opt.inner_drop_tolerance > R(0) && !ws.pat.empty()) {
				R nrm2 = R(0);
				for (std::size_t t = 0; t < ws.gt.size(); ++t) {
					const R a = vcp::tsparse_scalar::abs_value(ws.gt[t]);
					nrm2 += a * a;
				}
				const R tau2 = opt.inner_drop_tolerance * opt.inner_drop_tolerance;
				std::size_t kept = 0;
				for (std::size_t t = 0; t < ws.pat.size(); ++t) {
					const R a = vcp::tsparse_scalar::abs_value(ws.gt[t]);
					if (a * a <= tau2 * nrm2) {
						ws.inpat[static_cast<std::size_t>(ws.pat[t])] = char(0);
					} else {
						ws.pat[kept] = ws.pat[t];
						ws.gt[kept]  = ws.gt[t];
						ws.ci[kept]  = ws.ci[t];
						++kept;
					}
				}
				ws.pat.resize(kept);
				ws.gt.resize(kept);
				ws.ci.resize(kept);
			}
			// (e) gradient, eq. (24): the work vector accumulates
			// A_p^T g + A_p[:, i] (rows of A_p scattered with the g
			// coefficients + CSC column i; for symmetric A_p this is the
			// eq. (24) quantity, and NONSYMMETRIC input evaluates the same
			// expression as-is -- design OPEN-A7, quality policed
			// downstream).  The paper's factor 2 is rank-invariant and
			// omitted (directive SSA2 note).
			for (std::size_t t = 0; t < ws.pat.size(); ++t) {
				const std::size_t r = static_cast<std::size_t>(ws.pat[t]);
				const _T coeff = ws.gt[t];
				for (_Index q = rp[r]; q < rp[r + 1u]; ++q) {
					const std::size_t j = static_cast<std::size_t>(
					    ri[static_cast<std::size_t>(q)]);
					if (ws.wmark[j] == char(0)) {
						ws.wmark[j] = char(1);
						ws.touched.push_back(static_cast<_Index>(j));
					}
					ws.w[j] += rv[static_cast<std::size_t>(q)] * coeff;
				}
			}
			for (_Index q = ap[i]; q < ap[i + 1u]; ++q) {
				const std::size_t j = static_cast<std::size_t>(
				    ai[static_cast<std::size_t>(q)]);
				if (ws.wmark[j] == char(0)) {
					ws.wmark[j] = char(1);
					ws.touched.push_back(static_cast<_Index>(j));
				}
				ws.w[j] += av[static_cast<std::size_t>(q)];
			}
			// candidates: j < i, outside the current pattern, with a
			// NONZERO gradient magnitude (an exactly zero component carries
			// no first-order descent information and is skipped --
			// implementation decision, completion report OPEN-A11)
			ws.cand.clear();
			for (std::size_t t = 0; t < ws.touched.size(); ++t) {
				const _Index j = ws.touched[t];
				if (j < static_cast<_Index>(i) &&
				    ws.inpat[static_cast<std::size_t>(j)] == char(0)) {
					const R aw = vcp::tsparse_scalar::abs_value(
					    ws.w[static_cast<std::size_t>(j)]);
					if (aw > R(0)) {
						fsai2_candidate<_T, _Index> c;
						c.aw = aw;
						c.j = j;
						ws.cand.push_back(c);
					}
				}
			}
			// clean-on-exit for w / wmark (deterministic, row-local)
			for (std::size_t t = 0; t < ws.touched.size(); ++t) {
				const std::size_t j = static_cast<std::size_t>(ws.touched[t]);
				ws.wmark[j] = char(0);
				ws.w[j] = _T(0);
			}
			ws.touched.clear();
			if (ws.cand.empty()) {
				done = true;   // row fixed: no admissible candidate
				break;
			}
			std::sort(ws.cand.begin(), ws.cand.end(),
			          fsai2_candidate_before<_T, _Index>());
			const std::size_t nadd =
			    (ws.cand.size() < opt.entries_per_step) ? ws.cand.size()
			                                            : opt.entries_per_step;
			// (f) add the top-s candidates and re-solve (design SS3.1):
			// A_p[Pbar,Pbar] g = -A_p[Pbar,i], full-block gather as in
			// FSAI-1 (both triangles as stored, F-D3)
			for (std::size_t t = 0; t < nadd; ++t) {
				ws.pat.push_back(ws.cand[t].j);
				ws.inpat[static_cast<std::size_t>(ws.cand[t].j)] = char(1);
			}
			std::sort(ws.pat.begin(), ws.pat.end());
			const std::size_t m = ws.pat.size();
			ws.M.assign(m * m, _T(0));
			for (std::size_t c = 0; c < m; ++c) {
				const std::size_t j = static_cast<std::size_t>(ws.pat[c]);
				_Index q = ap[j];
				const _Index qe = ap[j + 1u];
				std::size_t t = 0;
				while (q < qe && t < m) {
					const _Index r = ai[static_cast<std::size_t>(q)];
					if (r < ws.pat[t])      { ++q; }
					else if (r > ws.pat[t]) { ++t; }
					else {
						ws.M[t + c * m] = av[static_cast<std::size_t>(q)];
						++q; ++t;
					}
				}
			}
			fsai2_gather_col_i_<_T, _Index>(i, ap, ai, av, ws.pat, ws.ci);
			ws.b.assign(m, _T(0));
			for (std::size_t t = 0; t < m; ++t) {
				ws.b[t] = _T(0) - ws.ci[t];
			}
			fsai2_dense_solve_<_T, _Index>(ws.M, m, opt.pivot_small_threshold,
			                               opt.pivot_lift_value, ws.b, nl);
			ws.gt.assign(ws.b.begin(), ws.b.end());
			// psi = a_ii + g^T A_p[Pbar, i] (design SS3.2, one inner product)
			psi = aii;
			for (std::size_t t = 0; t < m; ++t) {
				psi += ws.gt[t] * ws.ci[t];
			}
			fsai2_psi_lift_<_T, _Index>(psi, opt.pivot_small_threshold,
			                            opt.pivot_lift_value, nl);
		}
		if (converged) row_converged[i] = _Index(1);
		if (!done)     row_capped[i]    = _Index(1);   // full k_iter budget used

		// --- terminal post-filtration, applied ONCE (design SS3.5).
		// Minimal replica of the FSAI-1 dual threshold + Dhat-scale
		// restoration fused inside fsai_process_row_ (read-only, completion
		// report OPEN-A10).  The first threshold is evaluated in the SQUARED
		// comparison form |u_t|^2 < tau_post^2 * ||u||_2^2 -- mathematically
		// equivalent to the FSAI-1 |u_t| < tau_post * ||u||_2 for
		// nonnegative operands and sqrt-free (directive invariant 11; same
		// rewrite precedent as the FSAI-1 prefilter, design OPEN-F2).
		_T dh = psi;   // Dhat_ii := psi_k,i (design SS3.3, F2-D6)
		_Index dropped = _Index(0);
		_Index skipped = _Index(0);
		const std::size_t m = ws.pat.size() + 1u;   // + unit diagonal
		ws.u.assign(ws.gt.begin(), ws.gt.end());
		ws.u.push_back(_T(1));                      // explicit unit diagonal
		ws.keep.assign(m, char(1));
		if ((opt.postfilter_tolerance > R(0) ||
		     opt.postfilter_max_nnz_per_row > 0u) && m > 1u) {
			if (opt.postfilter_tolerance > R(0)) {
				R nrm2 = R(0);
				for (std::size_t t = 0; t < m; ++t) {
					const R a = vcp::tsparse_scalar::abs_value(ws.u[t]);
					nrm2 += a * a;
				}
				const R tau2 = opt.postfilter_tolerance * opt.postfilter_tolerance;
				for (std::size_t t = 0; t < m - 1u; ++t) {
					const R a = vcp::tsparse_scalar::abs_value(ws.u[t]);
					if (a * a < tau2 * nrm2) {
						ws.keep[t] = char(0);
					}
				}
			}
			// mmax top-magnitude selection over the survivors (diagonal not
			// counted, never dropped); equal magnitudes keep the earlier
			// (= lower column index) candidate -- stable sort on insertion
			// order, FSAI-1 / AINV tie rule
			if (opt.postfilter_max_nnz_per_row > 0u) {
				ws.sel.clear();
				for (std::size_t t = 0; t < m - 1u; ++t) {
					if (ws.keep[t] != char(0)) ws.sel.push_back(static_cast<_Index>(t));
				}
				if (ws.sel.size() > opt.postfilter_max_nnz_per_row) {
					fsai2_row_workspace<_T, _Index>& wsr = ws;
					std::stable_sort(ws.sel.begin(), ws.sel.end(),
						[&wsr](const _Index a, const _Index b) {
							return vcp::tsparse_scalar::abs_value(
							           wsr.u[static_cast<std::size_t>(a)]) >
							       vcp::tsparse_scalar::abs_value(
							           wsr.u[static_cast<std::size_t>(b)]);
						});
					for (std::size_t t = opt.postfilter_max_nnz_per_row;
					     t < ws.sel.size(); ++t) {
						ws.keep[static_cast<std::size_t>(ws.sel[t])] = char(0);
					}
				}
			}
			// dropped set eps_i (ascending) and the sqrt-free scale
			// restoration: Dhat_ii <- Dhat_ii * (1 + (u[E]^T A_p[E,E] u[E])
			// / psi) -- the FSAI-1 form 1 + d * quad with d = 1/Dhat
			// written for Dhat = psi (design SS3.5; psi is lifted nonzero)
			ws.eloc.clear();
			ws.eglob.clear();
			for (std::size_t t = 0; t < m - 1u; ++t) {
				if (ws.keep[t] == char(0)) {
					ws.eloc.push_back(static_cast<_Index>(t));
					ws.eglob.push_back(ws.pat[t]);
				}
			}
			dropped = static_cast<_Index>(ws.eloc.size());
			if (!ws.eloc.empty()) {
				_T squad = _T(0);
				const std::size_t ne = ws.eloc.size();
				for (std::size_t bb = 0; bb < ne; ++bb) {
					const std::size_t j = static_cast<std::size_t>(ws.eglob[bb]);
					const _T ub = ws.u[static_cast<std::size_t>(ws.eloc[bb])];
					_Index q = ap[j];
					const _Index qe = ap[j + 1u];
					std::size_t t = 0;
					while (q < qe && t < ne) {
						const _Index r = ai[static_cast<std::size_t>(q)];
						if (r < ws.eglob[t])      { ++q; }
						else if (r > ws.eglob[t]) { ++t; }
						else {
							squad += ws.u[static_cast<std::size_t>(ws.eloc[t])]
							         * av[static_cast<std::size_t>(q)] * ub;
							++q; ++t;
						}
					}
				}
				const _T cfac = _T(1) + squad / psi;
				// restoration only when the sign test holds (success-side
				// gate; NaN falls to the skip side) -- FSAI-1 regime
				if (vcp::tsparse_scalar::real_part(cfac) > R(0)) {
					dh = dh * cfac;
				} else {
					skipped = _Index(1);
				}
			}
		}

		// --- slice write: kept entries ascending (pattern < i, then the
		// diagonal); exact zeros have no stored representation (spmats
		// invariant, same rule as FSAI-1)
		const std::size_t base = slice_off[i];
		std::size_t cnt = 0;
		for (std::size_t t = 0; t < m; ++t) {
			if (ws.keep[t] == char(0)) continue;
			const _T v = ws.u[t];
			if (!(v == _T(0))) {
				u_ind[base + cnt] = (t + 1u < m) ? ws.pat[t] : static_cast<_Index>(i);
				u_val[base + cnt] = v;
				++cnt;
			}
		}
		row_count[i] = static_cast<_Index>(cnt);
		d_hat[i] = dh;
		row_lifts[i] = nl;
		row_dropped[i] = dropped;
		row_restore_skipped[i] = skipped;

		// clean-on-exit for inpat
		for (std::size_t t = 0; t < ws.pat.size(); ++t) {
			ws.inpat[static_cast<std::size_t>(ws.pat[t])] = char(0);
		}
	}

	// -----------------------------------------------------------------------
	// fsai2_core_: shared construction body of both overloads.  perm is
	// already validated (bijection); U0, when present, is finalized and
	// validated (square, dimension n, unit upper triangular) and lives in
	// the SAME permuted frame as the output factors.
	//   pass 1 (serial):   caps cap_i = |P0_i| + k_iter*s + 1 -> slice offsets
	//   pass 2 (parallel): per-row adaptive computation into own slices
	//   pass 3 (serial):   deterministic row-order compaction + assembly
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai2_core_(
	    const spmats<_T, _Index>& A,
	    const std::vector<_Index>& perm,
	    const spmats<_T, _Index>* U0,
	    const fsai_adaptive_options<_T>& opt,
	    spmats<_T, _Index>& U,
	    spmats<_T, _Index>& D,
	    fsai_adaptive_result<_T, _Index>& out)
	{
		const _Index n = A.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		// A_p = P^T A P, explicit (FSAI-1 helper reused; identity keeps A)
		spmats<_T, _Index> Ap_store;
		const bool ident = spmats_fsai_detail::fsai_perm_is_identity_<_Index>(perm);
		if (!ident) {
			spmats_fsai_detail::fsai_build_permuted_copy_<_T, _Index>(A, perm, Ap_store);
		}
		const spmats<_T, _Index>& Ap = ident ? A : Ap_store;
		const spmats<_T, _Index> Apc = Ap.as_csc();
		const spmats<_T, _Index> Apr = Ap.as_csr();

		const bool has_u0 = (U0 != 0);
		spmats<_T, _Index> U0c_store;
		if (has_u0) U0c_store = U0->as_csc();
		const std::vector<_Index>* u0p = has_u0 ? &U0c_store.outer_index() : 0;
		const std::vector<_Index>* u0i = has_u0 ? &U0c_store.inner_index() : 0;
		const std::vector<_T>*     u0v = has_u0 ? &U0c_store.values()      : 0;

		// --- pass 1 (serial): capped slice offsets (design SS4; F2-D4).
		// |P0_i| counts the stored entries of the initial column INCLUDING
		// its unit diagonal (1 for the diagonal init), so the cap has one
		// spare cell over the reachable maximum -- writes stay disjoint by
		// construction.
		const std::size_t addcap = opt.max_iterations * opt.entries_per_step;
		std::vector<std::size_t> slice_off(un + 1u, 0u);
		for (std::size_t i = 0; i < un; ++i) {
			const std::size_t p0 = has_u0
			    ? static_cast<std::size_t>((*u0p)[i + 1u]) -
			      static_cast<std::size_t>((*u0p)[i])
			    : 1u;
			slice_off[i + 1u] = slice_off[i] + p0 + addcap + 1u;
		}

		// --- pass 2: per-row adaptive solves.  Writes go only to the row's
		// own cap slice and per-row cells -- see the bit-identity contract
		// in the file header.
		std::vector<_Index> u_ind(slice_off[un]);
		std::vector<_T>     u_val(slice_off[un]);
		std::vector<_Index> row_count(un, _Index(0));
		std::vector<_T>     d_hat(un, _T(0));
		std::vector<_Index> row_lifts(un, _Index(0));
		std::vector<_Index> row_dropped(un, _Index(0));
		std::vector<_Index> row_restore_skipped(un, _Index(0));
		std::vector<_Index> row_converged(un, _Index(0));
		std::vector<_Index> row_capped(un, _Index(0));
		{
			const std::vector<_Index>& ap = Apc.outer_index();
			const std::vector<_Index>& ai = Apc.inner_index();
			const std::vector<_T>&     av = Apc.values();
			const std::vector<_Index>& rp = Apr.outer_index();
			const std::vector<_Index>& ri = Apr.inner_index();
			const std::vector<_T>&     rv = Apr.values();
			const long long n_ll = static_cast<long long>(un);
#if VCP_SPMATS_USE_OPENMP
			// design SS4: rows are fully independent; each thread owns one
			// workspace (declared inside the parallel region) and writes
			// only its rows' own cap slices, so the parallel run is
			// bit-identical to the serial one regardless of thread count
			// and schedule.  Small problems stay serial (n >= 1024
			// threshold inherited from FSAI-1, provisional OPEN-F7).
			#pragma omp parallel if(n_ll >= 1024)
			{
				fsai2_row_workspace<_T, _Index> ws;
				#pragma omp for schedule(dynamic)
				for (long long i = 0; i < n_ll; ++i) {
					fsai2_process_row_<_T, _Index>(
					    static_cast<std::size_t>(i), un, ap, ai, av, rp, ri, rv,
					    has_u0, u0p, u0i, u0v, opt, ws, slice_off,
					    u_ind, u_val, row_count, d_hat, row_lifts,
					    row_dropped, row_restore_skipped,
					    row_converged, row_capped);
				}
			}
#else
			fsai2_row_workspace<_T, _Index> ws;
			for (long long i = 0; i < n_ll; ++i) {
				fsai2_process_row_<_T, _Index>(
				    static_cast<std::size_t>(i), un, ap, ai, av, rp, ri, rv,
				    has_u0, u0p, u0i, u0v, opt, ws, slice_off,
				    u_ind, u_val, row_count, d_hat, row_lifts,
				    row_dropped, row_restore_skipped,
				    row_converged, row_capped);
			}
#endif
		}

		// --- pass 3 (serial): deterministic compaction in row order
		std::vector<_Index> ucp(un + 1u, _Index(0));
		for (std::size_t i = 0; i < un; ++i) {
			ucp[i + 1u] = ucp[i] + row_count[i];
		}
		const std::size_t nnz_u = static_cast<std::size_t>(ucp[un]);
		std::vector<_Index> uri(nnz_u);
		std::vector<_T>     uva(nnz_u);
		for (std::size_t i = 0; i < un; ++i) {
			const std::size_t src = slice_off[i];
			const std::size_t dst = static_cast<std::size_t>(ucp[i]);
			const std::size_t cnt = static_cast<std::size_t>(row_count[i]);
			for (std::size_t t = 0; t < cnt; ++t) {
				uri[dst + t] = u_ind[src + t];
				uva[dst + t] = u_val[src + t];
			}
		}
		U.assign_csc(n, n, ucp, uri, uva);

		// Dhat: diagonal store, exact zeros unstored (FSAI-1 / AINV
		// precedent -- possible only through overflow / underflow, since
		// psi itself is lifted nonzero)
		{
			std::vector<_Index> dp(un + 1u, _Index(0)), di;
			std::vector<_T> dv;
			for (std::size_t k = 0; k < un; ++k) {
				if (!(d_hat[k] == _T(0))) {
					di.push_back(static_cast<_Index>(k));
					dv.push_back(d_hat[k]);
				}
				dp[k + 1u] = static_cast<_Index>(di.size());
			}
			D.assign_csc(n, n, dp, di, dv);
		}

		// diagnostics: combined serially in row order (deterministic; the
		// parallel loop never touches these aggregates)
		out.status = fsai_status::success;
		out.nnz_U = static_cast<_Index>(nnz_u);
		out.n_pivot_modifications = _Index(0);
		out.first_modified_row = _Index(-1);
		out.n_postfilter_dropped = _Index(0);
		out.n_postfilter_restore_skipped = _Index(0);
		out.n_rows_converged = _Index(0);
		out.n_rows_capped = _Index(0);
		for (std::size_t i = 0; i < un; ++i) {
			out.n_pivot_modifications = out.n_pivot_modifications + row_lifts[i];
			if (row_lifts[i] > _Index(0) && out.first_modified_row < _Index(0)) {
				out.first_modified_row = static_cast<_Index>(i);
			}
			out.n_postfilter_dropped = out.n_postfilter_dropped + row_dropped[i];
			out.n_postfilter_restore_skipped =
			    out.n_postfilter_restore_skipped + row_restore_skipped[i];
			out.n_rows_converged = out.n_rows_converged + row_converged[i];
			out.n_rows_capped = out.n_rows_capped + row_capped[i];
		}
	}

	// -----------------------------------------------------------------------
	// fsai2_empty_success_: shared n = 0 result (both overloads)
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline fsai_adaptive_result<_T, _Index> fsai2_empty_success_(
	    spmats<_T, _Index>& U, spmats<_T, _Index>& D,
	    std::vector<_Index>& perm)
	{
		const std::vector<_Index> ep(1u, _Index(0));
		U.assign_csc(_Index(0), _Index(0), ep, std::vector<_Index>(), std::vector<_T>());
		D.assign_csc(_Index(0), _Index(0), ep, std::vector<_Index>(), std::vector<_T>());
		perm.clear();
		fsai_adaptive_result<_T, _Index> out;
		out.status = fsai_status::success;
		out.nnz_U = _Index(0);
		out.n_pivot_modifications = _Index(0);
		out.first_modified_row = _Index(-1);
		out.n_postfilter_dropped = _Index(0);
		out.n_postfilter_restore_skipped = _Index(0);
		out.n_rows_converged = _Index(0);
		out.n_rows_capped = _Index(0);
		return out;
	}

	// -----------------------------------------------------------------------
	// fsai2_reject_invalid_: shared invalid_input rejection (outputs
	// cleared, non-throwing -- FSAI-1 outer regime)
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline fsai_adaptive_result<_T, _Index> fsai2_reject_invalid_(
	    spmats<_T, _Index>& U, spmats<_T, _Index>& D,
	    std::vector<_Index>& perm)
	{
		U.resize(_Index(0), _Index(0));
		D.resize(_Index(0), _Index(0));
		perm.clear();
		fsai_adaptive_result<_T, _Index> out;
		out.status = fsai_status::invalid_input;
		return out;
	}

} // namespace spmats_fsai_adaptive_detail

// ---------------------------------------------------------------------------
// policy_fsai_adaptive_with_info (diagonal init): NVI outer (non-virtual).
// Finalize guarantee + squareness entry check -- NON-throwing (invalid_input,
// FSAI-1 outer regime).  Must never be overridden -- override
// policy_fsai_adaptive_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_adaptive_result<_T, _Index> spmats<_T, _Index>::policy_fsai_adaptive_with_info(
	spmats<_T, _Index>& U,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const fsai_adaptive_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize()) {
		return spmats_fsai_adaptive_detail::fsai2_reject_invalid_<_T, _Index>(U, D, perm);
	}
	return policy_fsai_adaptive_with_info_impl(U, D, perm, opt);
}

// ---------------------------------------------------------------------------
// policy_fsai_adaptive_with_info (initial-value overload, F2-D3): NVI outer.
// U0 is a unit upper triangular factor (FSAI-1 / FSAI-2 output) in the SAME
// permuted frame as perm0; opt.ordering is IGNORED and perm = perm0 is
// returned as-is (design SS3.6 -- perm consistency is the CALLER'S
// responsibility).  Entry validation (all invalid_input, non-throwing):
// A square, U0 square of the same dimension, perm0 a bijection, U0 upper
// triangular with an explicit unit diagonal.  No D0 is taken: psi_0 is
// recomputed from the U0 values (design SS3.3 consequence, OPEN-A5).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_adaptive_result<_T, _Index> spmats<_T, _Index>::policy_fsai_adaptive_with_info(
	spmats<_T, _Index>& U,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const spmats<_T, _Index>& U0,
	const std::vector<_Index>& perm0,
	const fsai_adaptive_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (!U0.is_finalized()) U0.finalize();
	const _Index n = A.rowsize();
	if (A.rowsize() != A.columnsize() ||
	    U0.rowsize() != U0.columnsize() || U0.rowsize() != n ||
	    !spmats_fsai_detail::fsai_perm_is_bijection_<_Index>(perm0, n)) {
		return spmats_fsai_adaptive_detail::fsai2_reject_invalid_<_T, _Index>(U, D, perm);
	}
	{
		const spmats<_T, _Index> U0c = U0.as_csc();
		if (!spmats_fsai_adaptive_detail::fsai2_u0_is_unit_upper_<_T, _Index>(
		        U0c.outer_index(), U0c.inner_index(), U0c.values(), n)) {
			return spmats_fsai_adaptive_detail::fsai2_reject_invalid_<_T, _Index>(U, D, perm);
		}
	}
	return policy_fsai_adaptive_with_info_impl(U, D, perm, U0, perm0, opt);
}

// ---------------------------------------------------------------------------
// policy_fsai_adaptive_with_info_impl (diagonal init): virtual algorithm
// body.  Ordering resolution and the bijection firewall follow the FSAI-1
// impl verbatim (auto_select -> natural, OPEN-F13 convention; unknown enum
// -> invalid_input; firewall failure -> internal_error).  Numerical events
// are lifted silently and never a failure status (F-D3/F-D4 inherited).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_adaptive_result<_T, _Index> spmats<_T, _Index>::policy_fsai_adaptive_with_info_impl(
	spmats<_T, _Index>& U,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const fsai_adaptive_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;
	fsai_adaptive_result<_T, _Index> out;
	try {
		// [SLU-HK1 STOP-4 ruling] FSAI x unsigned Index guard: same contract
		// as the static-FSAI entry guard (spmats_fsai_impl.hpp) -- the shared
		// SLU ordering layer is signed-Index only, so unsigned TUs must stay
		// buildable (see fsai_ordering_perm_) and get the honest runtime
		// failure through the existing rejection helper BEFORE any work.
		// The fsai_adaptive track should formally adopt this sign dispatch
		// in its own design revision.
		if (!std::is_signed<_Index>::value) {
			return spmats_fsai_adaptive_detail::fsai2_reject_invalid_<_T, _Index>(U, D, perm);
		}

		const _Index n = A.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);
		if (n == _Index(0)) {
			return spmats_fsai_adaptive_detail::fsai2_empty_success_<_T, _Index>(U, D, perm);
		}
		sparse_chol_ordering ord = opt.ordering;
		if (ord == sparse_chol_ordering::auto_select) {
			ord = sparse_chol_ordering::natural;
		}
		perm.resize(un);
		for (std::size_t k = 0; k < un; ++k) perm[k] = static_cast<_Index>(k);
		if (ord != sparse_chol_ordering::natural) {
			// [SLU-HK1 STOP-4 ruling] SFINAE-split ordering resolution shared
			// with the static FSAI (spmats_fsai_detail::fsai_ordering_perm_):
			// identical calls in identical order for signed Index; false =
			// unknown enum value -> invalid_input (unchanged).
			if (!spmats_fsai_detail::fsai_ordering_perm_<_T, _Index>(A, ord, n, perm)) {
				return spmats_fsai_adaptive_detail::fsai2_reject_invalid_<_T, _Index>(U, D, perm);
			}
		}
		if (!spmats_fsai_detail::fsai_perm_is_bijection_<_Index>(perm, n)) {
			U.resize(_Index(0), _Index(0));
			D.resize(_Index(0), _Index(0));
			perm.clear();
			out = fsai_adaptive_result<_T, _Index>();
			out.status = fsai_status::internal_error;
			return out;
		}
		spmats_fsai_adaptive_detail::fsai2_core_<_T, _Index>(
		    A, perm, static_cast<const spmats<_T, _Index>*>(0), opt, U, D, out);
		return out;
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		out = fsai_adaptive_result<_T, _Index>();
		out.status = fsai_status::internal_error;
		return out;
	}
}

// ---------------------------------------------------------------------------
// policy_fsai_adaptive_with_info_impl (initial-value overload): virtual
// algorithm body.  perm0 is reused verbatim (opt.ordering ignored, design
// SS3.6); all structural validation already ran in the outer.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_adaptive_result<_T, _Index> spmats<_T, _Index>::policy_fsai_adaptive_with_info_impl(
	spmats<_T, _Index>& U,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const spmats<_T, _Index>& U0,
	const std::vector<_Index>& perm0,
	const fsai_adaptive_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;
	fsai_adaptive_result<_T, _Index> out;
	try {
		// [SLU-HK1 STOP-4 ruling] FSAI x unsigned Index guard: same honest
		// rejection as the diagonal-init entry above -- keeps the unsigned
		// contract consistent across both construction entries (this path
		// does not reach the ordering layer, but an unsigned adaptive FSAI
		// has never been validated by any track).
		if (!std::is_signed<_Index>::value) {
			return spmats_fsai_adaptive_detail::fsai2_reject_invalid_<_T, _Index>(U, D, perm);
		}

		const _Index n = A.rowsize();
		if (n == _Index(0)) {
			return spmats_fsai_adaptive_detail::fsai2_empty_success_<_T, _Index>(U, D, perm);
		}
		perm = perm0;
		spmats_fsai_adaptive_detail::fsai2_core_<_T, _Index>(
		    A, perm, &U0, opt, U, D, out);
		return out;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		out = fsai_adaptive_result<_T, _Index>();
		out.status = fsai_status::internal_error;
		return out;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_FSAI_ADAPTIVE_IMPL_HPP
