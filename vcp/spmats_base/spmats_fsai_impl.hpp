// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_fsai_impl.hpp
// Out-of-line policy method definitions for the static FSAI factored sparse
// approximate inverse (FSAI-1; design sandbox/docs/design/fsai_design_v0.md).
// Static FSAI after [JFSG15] (Janna/Ferronato/Sartoretto/Gambolati, ACM TOMS
// 41(2), 2015): pattern generation ([JFSG15] Algorithm 2), per-row dense
// solves A_p[P_i,P_i] g_i = e_m ([JFSG15] Algorithm 1, eq. (11)),
// post-filtration ([JFSG15] Algorithm 5), output as the sqrt-free triple
//   U (unit upper triangular, = Ghat^T), Dhat (diagonal), perm (new->old)
// with R := P U Dhat^{-1} U^T P^T ~ A^{-1} NEVER materialized (design SS3.3;
// the sqrt-free conversion is THIS DESIGN'S derivation, not in [JFSG15]).
//   policy_fsai_with_info               construction (factors only)
//   policy_fsai_residual_norm_estimate  non-guaranteed ||I - R A||_inf,
//                                       delegated to the AINV helper with
//                                       Z = W = U (F-D10)
//   policy_fsai_apply                   z = R r, delegated likewise
// This file is included inside spmats.hpp AFTER the closing brace of
// spmats<_T,_Index>, alongside spmats_ainv_impl.hpp.  The types live in
// spmats_base/spmats_fsai.hpp (included before the class body), which also
// hosts the sparse-layer OpenMP guard VCP_SPMATS_USE_OPENMP (design SS6.1).
//
// OpenMP contract (design SS6.3): the parallel row loops below are
// BIT-IDENTICAL to the serial execution by construction -- every row reads
// only A_p and the pattern S, and writes only its own slice fixed by pass 1;
// there is no shared accumulation, no atomic, no reduction, and no
// thread-count-dependent branch.  Diagnostics are per-row arrays combined
// serially in row order after the loop.

#ifndef VCP_SPMATS_FSAI_IMPL_HPP
#define VCP_SPMATS_FSAI_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <exception>
#include <type_traits>
#include <vector>

#include <vcp/spmats_base/spmats_fsai.hpp>
#include <vcp/spmats_base/spmats_ainv_impl.hpp>   // spmats_ainv_detail helpers (read-only reuse)

namespace vcp {

namespace spmats_fsai_detail {

	// -----------------------------------------------------------------------
	// fsai_from_ainv_status_: status mapping for the AINV delegation (the two
	// enums are isomorphic by design, OPEN-F10).
	// -----------------------------------------------------------------------
	inline fsai_status fsai_from_ainv_status_(const ainv_status s) {
		switch (s) {
		case ainv_status::success:        return fsai_status::success;
		case ainv_status::invalid_input:  return fsai_status::invalid_input;
		case ainv_status::internal_error: return fsai_status::internal_error;
		}
		return fsai_status::internal_error;
	}

	// -----------------------------------------------------------------------
	// fsai_perm_is_bijection_: O(n) permutation validation (chol SS7.5-1
	// firewall precedent; also the consumer-side perm check of the estimate /
	// apply wrappers).
	// -----------------------------------------------------------------------
	template <typename _Index>
	inline bool fsai_perm_is_bijection_(const std::vector<_Index>& perm,
	                                    const _Index n)
	{
		const std::size_t un = static_cast<std::size_t>(n);
		if (perm.size() != un) return false;
		std::vector<char> seen(un, char(0));
		for (std::size_t k = 0; k < un; ++k) {
			const _Index p = perm[k];
			if (p < _Index(0) || p >= n) return false;
			const std::size_t sp = static_cast<std::size_t>(p);
			if (seen[sp] != char(0)) return false;
			seen[sp] = char(1);
		}
		return true;
	}

	template <typename _Index>
	inline bool fsai_perm_is_identity_(const std::vector<_Index>& perm) {
		for (std::size_t k = 0; k < perm.size(); ++k) {
			if (perm[k] != static_cast<_Index>(k)) return false;
		}
		return true;
	}

	// -----------------------------------------------------------------------
	// fsai_build_permuted_copy_: A_p = P^T A P built EXPLICITLY (design F2-1:
	// one reordering copy of nnz(A) size -- unlike R, this does not destroy
	// sparsity).  perm is new->old with P(perm[k],k) = 1 (chol convention),
	// so entry (r,c) of A lands at (pinv[r], pinv[c]).  Entries are unique,
	// therefore the COO finalize produces a canonical result independent of
	// insertion order (deterministic).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai_build_permuted_copy_(
	    const spmats<_T, _Index>& A,
	    const std::vector<_Index>& perm,
	    spmats<_T, _Index>& Ap)
	{
		const _Index n = A.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);
		std::vector<_Index> pinv(un);   // old -> new
		for (std::size_t k = 0; k < un; ++k) {
			pinv[static_cast<std::size_t>(perm[k])] = static_cast<_Index>(k);
		}
		const spmats<_T, _Index> Ac = A.as_csc();
		const std::vector<_Index>& ap = Ac.outer_index();
		const std::vector<_Index>& ai = Ac.inner_index();
		const std::vector<_T>&     av = Ac.values();
		Ap.resize(n, n);
		for (_Index j = _Index(0); j < n; ++j) {
			const std::size_t sj = static_cast<std::size_t>(j);
			for (_Index q = ap[sj]; q < ap[sj + 1u]; ++q) {
				const std::size_t sq = static_cast<std::size_t>(q);
				Ap.add(pinv[static_cast<std::size_t>(ai[sq])],
				       pinv[sj], av[sq]);
			}
		}
		Ap.finalize();
	}

	// -----------------------------------------------------------------------
	// fsai_prefilter_rows_: A' row patterns from the CSR of A_p ([JFSG15]
	// eq. (17) prefiltration).  The drop test is the SQUARED comparison
	//   |a_ij|^2 < tau_pre^2 * |a_ii| * |a_jj|
	// (mathematically equivalent to the eq. (17) sqrt form for nonnegative
	// operands -- sqrt completely avoided, design OPEN-F2).  The absolute
	// values on the diagonal are a safeguard for zero / negative diagonals
	// (NOT in [JFSG15], which assumes SPD -- this implementation's decision,
	// design SS3.1).  A structurally missing diagonal gives |a_ii| = 0 and
	// then no drop against row i (safe direction).  Diagonal entries are
	// never dropped.  tau_pre = 0 keeps every entry (A' = A).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai_prefilter_rows_(
	    const spmats<_T, _Index>& Apr,   // finalized CSR of A_p
	    const typename vcp::tsparse_scalar::real_type<_T>::type& tau_pre,
	    std::vector<_Index>& arp,        // out: A' row pointers (n+1)
	    std::vector<_Index>& ari)        // out: A' column indices (ascending per row)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		const _Index n = Apr.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);
		const std::vector<_Index>& rp = Apr.outer_index();
		const std::vector<_Index>& ri = Apr.inner_index();
		const std::vector<_T>&     rv = Apr.values();

		arp.assign(un + 1u, _Index(0));
		ari.clear();

		std::vector<R> adiag(un, R(0));
		for (std::size_t i = 0; i < un; ++i) {
			for (_Index q = rp[i]; q < rp[i + 1u]; ++q) {
				const std::size_t sq = static_cast<std::size_t>(q);
				if (ri[sq] == static_cast<_Index>(i)) {
					adiag[i] = vcp::tsparse_scalar::abs_value(rv[sq]);
					break;
				}
			}
		}
		const R tau2 = tau_pre * tau_pre;
		for (std::size_t i = 0; i < un; ++i) {
			for (_Index q = rp[i]; q < rp[i + 1u]; ++q) {
				const std::size_t sq = static_cast<std::size_t>(q);
				const _Index j = ri[sq];
				if (j != static_cast<_Index>(i)) {
					const R a = vcp::tsparse_scalar::abs_value(rv[sq]);
					if (a * a < tau2 * adiag[i] * adiag[static_cast<std::size_t>(j)]) {
						continue;   // prefiltered
					}
				}
				ari.push_back(j);
			}
			arp[i + 1u] = static_cast<_Index>(ari.size());
		}
	}

	// -----------------------------------------------------------------------
	// fsai_pattern_: symbolic power recursion B_0 = I, B_i = Low(B_{i-1} A')
	// ([JFSG15] eq. (16), Huckle secondary citation; pattern-only, no
	// values).  Each step is a marker-based row merge followed by an
	// ascending sort (deterministic).  The recursion stops early once
	// nnz(B_i) >= mu_max * nnz(A) (density guard; mu_max = 0 means
	// unlimited).  The final pattern is S = pattern(B_k) union the diagonal
	// (the diagonal is always included -- existence-guarantee premise,
	// design SS3.1), stored as ascending rows P_i with row pointers (this is
	// the pass-1 output fixing the per-row output slices, design SS6.2).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai_pattern_(
	    const _Index n,
	    const std::vector<_Index>& arp,   // A' row pointers
	    const std::vector<_Index>& ari,   // A' column indices
	    const std::size_t pattern_power,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& mu_max,
	    const _Index nnz_A,
	    std::vector<_Index>& s_ptr,       // out: S row pointers (n+1)
	    std::vector<_Index>& s_ind)       // out: S column indices (ascending, diag included)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		const std::size_t un = static_cast<std::size_t>(n);

		// B_0 = I
		std::vector<std::vector<_Index> > b(un);
		for (std::size_t i = 0; i < un; ++i) b[i].assign(1u, static_cast<_Index>(i));

		// nnz(A) as a real value for the density guard (counts are integers,
		// exactly representable through double well beyond any practical nnz;
		// no double LITERAL is involved)
		const R nnzA_r = R(static_cast<double>(nnz_A));

		std::vector<char> mark(un, char(0));
		std::vector<_Index> row;
		for (std::size_t step = 0; step < pattern_power; ++step) {
			std::vector<std::vector<_Index> > bnew(un);
			std::size_t nnz_b = 0;
			for (std::size_t i = 0; i < un; ++i) {
				row.clear();
				for (std::size_t t = 0; t < b[i].size(); ++t) {
					const std::size_t c = static_cast<std::size_t>(b[i][t]);
					for (_Index q = arp[c]; q < arp[c + 1u]; ++q) {
						const _Index j = ari[static_cast<std::size_t>(q)];
						if (j > static_cast<_Index>(i)) continue;   // Low(): lower triangle only
						const std::size_t sjj = static_cast<std::size_t>(j);
						if (mark[sjj] == char(0)) {
							mark[sjj] = char(1);
							row.push_back(j);
						}
					}
				}
				std::sort(row.begin(), row.end());
				for (std::size_t t = 0; t < row.size(); ++t) {
					mark[static_cast<std::size_t>(row[t])] = char(0);
				}
				bnew[i] = row;
				nnz_b += row.size();
			}
			b.swap(bnew);
			if (mu_max > R(0) &&
			    R(static_cast<double>(nnz_b)) >= mu_max * nnzA_r) {
				break;   // density guard: keep B_i, stop recursing
			}
		}

		// S = pattern(B_k) union diagonal
		s_ptr.assign(un + 1u, _Index(0));
		s_ind.clear();
		for (std::size_t i = 0; i < un; ++i) {
			bool has_diag = false;
			for (std::size_t t = 0; t < b[i].size(); ++t) {
				if (b[i][t] == static_cast<_Index>(i)) { has_diag = true; }
				s_ind.push_back(b[i][t]);
			}
			if (!has_diag) s_ind.push_back(static_cast<_Index>(i));
			// rows of Low() are <= i and sorted; appending the diagonal (the
			// maximum) preserves ascending order
			s_ptr[i + 1u] = static_cast<_Index>(s_ind.size());
		}
	}

	// -----------------------------------------------------------------------
	// fsai_dense_ldl_solve_: dense in-place sqrt-free factorization + back
	// substitution for A_p[P_i,P_i] g = e_m ([JFSG15] eq. (11)).  The
	// factorization is the unit-lower elimination WITHOUT interchanges
	// (Doolittle); for a symmetric block this IS the sqrt-free LDL^T
	// (U = D L^T), and running it on the full gathered block keeps the
	// nonsymmetric case (F-D3) factoring the block as-is -- full-matrix
	// elimination is this implementation's decision, [JFSG15] assumes SPD
	// (noted per directive invariant 8).  Every diagonal pivot |p| <=
	// threshold is lifted sign-preservingly to +-lift (exact zero -> +lift;
	// AINV SS4.3 regime, not in [JFSG15]); the sign is read through
	// real_part (identity for real scalars).  The right-hand side e_m passes
	// the unit-lower forward solve unchanged (all leading components are
	// zero), so only the back substitution is performed.
	// M is column-major m x m and is destroyed.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai_dense_ldl_solve_(
	    std::vector<_T>& M, const std::size_t m,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& threshold,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& lift,
	    std::vector<_T>& g,
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
		g.assign(m, _T(0));
		g[m - 1u] = _T(1) / M[(m - 1u) + (m - 1u) * m];
		for (std::size_t r = m - 1u; r-- > 0u; ) {
			_T s = _T(0);
			for (std::size_t c = r + 1u; c < m; ++c) {
				s += M[r + c * m] * g[c];
			}
			g[r] = (_T(0) - s) / M[r + r * m];
		}
	}

	// -----------------------------------------------------------------------
	// fsai_row_workspace: per-thread scratch of the row loop (design SS6.2:
	// gather buffer + solve work + filter marks are thread-local; in the
	// parallel build one instance lives inside each thread's parallel
	// region).  All members are reset per row with row-determined values, so
	// the produced numbers cannot depend on which thread ran which row.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	struct fsai_row_workspace {
		std::vector<_T> M;          // m x m gather / factorization buffer
		std::vector<_T> g;          // solve result g_i
		std::vector<_T> u;          // unit-scaled row (g / g_m), diagonal 1
		std::vector<char> keep;     // post-filtration keep flags (size m)
		std::vector<_Index> sel;    // kept off-diagonal locals (mmax selection)
		std::vector<_Index> eloc;   // dropped locals, ascending
		std::vector<_Index> eglob;  // dropped globals, ascending
	};

	// -----------------------------------------------------------------------
	// fsai_process_row_: the complete per-row computation (gather -> solve ->
	// sqrt-free transform -> post-filtration -> slice write).  Reads ONLY
	// A_p (CSC), the pattern S and opt; writes ONLY the row-i slice
	// u_ind/u_val[s_ptr[i] .. s_ptr[i+1])) and the row-i cells of the per-row
	// output arrays -- the parallel bit-identity contract rests on this
	// (design SS6.2/6.3).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void fsai_process_row_(
	    const std::size_t i,
	    const std::vector<_Index>& ap,    // CSC col_ptr of A_p
	    const std::vector<_Index>& ai,    // CSC row_ind of A_p
	    const std::vector<_T>&     av,    // CSC values of A_p
	    const std::vector<_Index>& s_ptr,
	    const std::vector<_Index>& s_ind,
	    const fsai_options<_T>& opt,
	    fsai_row_workspace<_T, _Index>& ws,
	    std::vector<_Index>& u_ind,       // slice target (row-i slice only)
	    std::vector<_T>&     u_val,       // slice target (row-i slice only)
	    std::vector<_Index>& row_count,
	    std::vector<_T>&     d_hat,
	    std::vector<_Index>& row_lifts,
	    std::vector<_Index>& row_dropped,
	    std::vector<_Index>& row_restore_skipped)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		const std::size_t base = static_cast<std::size_t>(s_ptr[i]);
		const std::size_t m = static_cast<std::size_t>(s_ptr[i + 1u]) - base;
		const _Index* P = &s_ind[base];   // ascending, P[m-1] == i

		// gather: dense M(r,c) = A_p(P[r], P[c]) from BOTH triangles as
		// stored (the block is well defined for nonsymmetric input, F-D3;
		// [JFSG15] gathers only the lower triangle under its SPD assumption
		// -- full gather is this implementation's decision, directive F3-1).
		// Column c is a two-pointer merge of CSC column P[c] with P.
		ws.M.assign(m * m, _T(0));
		for (std::size_t c = 0; c < m; ++c) {
			const std::size_t j = static_cast<std::size_t>(P[c]);
			_Index q = ap[j];
			const _Index qe = ap[j + 1u];
			std::size_t t = 0;
			while (q < qe && t < m) {
				const _Index r = ai[static_cast<std::size_t>(q)];
				if (r < P[t])      { ++q; }
				else if (r > P[t]) { ++t; }
				else {
					ws.M[t + c * m] = av[static_cast<std::size_t>(q)];
					++q; ++t;
				}
			}
		}

		// solve A_p[P,P] g = e_m
		_Index nl = _Index(0);
		fsai_dense_ldl_solve_<_T, _Index>(ws.M, m,
		    opt.pivot_small_threshold, opt.pivot_lift_value, ws.g, nl);

		// sqrt-free transform (design SS3.3, THIS DESIGN's derivation, not
		// in [JFSG15]): d = g(m), Dhat_ii = 1/d, U column i = unit-scaled
		// row g/d.  A tiny |d| gets the same sign-preserving lift as an
		// elimination pivot (design SS3.4 tail rule -- this design's
		// safeguard, unconfirmed by the paper) so Dhat never receives a
		// structural zero from this path.
		_T d = ws.g[m - 1u];
		if (vcp::tsparse_scalar::abs_value(d) <= opt.pivot_small_threshold) {
			const R rd = vcp::tsparse_scalar::real_part(d);
			d = (rd < R(0)) ? _T(-opt.pivot_lift_value) : _T(opt.pivot_lift_value);
			nl = nl + _Index(1);
		}
		ws.u.assign(m, _T(0));
		for (std::size_t t = 0; t < m - 1u; ++t) ws.u[t] = ws.g[t] / d;
		ws.u[m - 1u] = _T(1);   // explicit unit diagonal

		// post-filtration ([JFSG15] Algorithm 5 dual threshold).  The test
		// |g_ij| < tau_post * ||g_i||_2 is SCALE-INVARIANT within the row,
		// so it is evaluated on the unit-scaled row u (equivalent to the
		// paper's G row; equivalence note per directive invariant 8).  The
		// 2-norm is the only real_type sqrt in this module (design OPEN-F2).
		_T dh = _T(1) / d;
		_Index dropped = _Index(0);
		_Index skipped = _Index(0);
		ws.keep.assign(m, char(1));
		if ((opt.postfilter_tolerance > R(0) ||
		     opt.postfilter_max_nnz_per_row > 0u) && m > 1u) {
			if (opt.postfilter_tolerance > R(0)) {
				R nrm2 = R(0);
				for (std::size_t t = 0; t < m; ++t) {
					const R a = vcp::tsparse_scalar::abs_value(ws.u[t]);
					nrm2 += a * a;
				}
				const R nrm = vcp::tsparse_scalar::sqrt_value(nrm2);
				for (std::size_t t = 0; t < m - 1u; ++t) {
					if (vcp::tsparse_scalar::abs_value(ws.u[t]) <
					    opt.postfilter_tolerance * nrm) {
						ws.keep[t] = char(0);
					}
				}
			}
			// mmax top-magnitude selection over the survivors (diagonal not
			// counted, never dropped).  Equal magnitudes keep the
			// earlier-arrived (= lower column index) candidate: stable sort
			// on insertion order -- tie rule NOT in [JFSG15], same decision
			// as AINV.
			if (opt.postfilter_max_nnz_per_row > 0u) {
				ws.sel.clear();
				for (std::size_t t = 0; t < m - 1u; ++t) {
					if (ws.keep[t] != char(0)) ws.sel.push_back(static_cast<_Index>(t));
				}
				if (ws.sel.size() > opt.postfilter_max_nnz_per_row) {
					std::stable_sort(ws.sel.begin(), ws.sel.end(),
						[&ws](const _Index a, const _Index b) {
							return vcp::tsparse_scalar::abs_value(
							           ws.u[static_cast<std::size_t>(a)]) >
							       vcp::tsparse_scalar::abs_value(
							           ws.u[static_cast<std::size_t>(b)]);
						});
					for (std::size_t t = opt.postfilter_max_nnz_per_row;
					     t < ws.sel.size(); ++t) {
						ws.keep[static_cast<std::size_t>(ws.sel[t])] = char(0);
					}
				}
			}
			// dropped set eps_i (ascending) and the sqrt-free scale
			// restoration (design SS3.5, THIS DESIGN's equivalent of
			// [JFSG15] eq. (38): Dhat_ii <- Dhat_ii * (1 + eps^T A eps),
			// U column unchanged).  In the paper's G scaling eps =
			// sqrt(g_m) * u[E], hence eps^T A eps = d * (u[E]^T A_p[E,E]
			// u[E]) -- computed here on the unit-scaled values.
			ws.eloc.clear();
			ws.eglob.clear();
			for (std::size_t t = 0; t < m - 1u; ++t) {
				if (ws.keep[t] == char(0)) {
					ws.eloc.push_back(static_cast<_Index>(t));
					ws.eglob.push_back(P[t]);
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
				const _T cfac = _T(1) + d * squad;
				// pathological (non-SPD) restoration factor: apply only when
				// certified positive; otherwise skip and record (design
				// SS3.5 safeguard -- this design's decision).  The positive
				// test is the success-side gate (NaN falls to the skip side).
				if (vcp::tsparse_scalar::real_part(cfac) > R(0)) {
					dh = dh * cfac;
				} else {
					skipped = _Index(1);
				}
			}
		}

		// slice write: kept entries in ascending row order; exact zeros have
		// no stored representation (spmats invariant -- same rule as the
		// AINV store).  The unit diagonal is always kept and nonzero.
		std::size_t cnt = 0;
		for (std::size_t t = 0; t < m; ++t) {
			if (ws.keep[t] == char(0)) continue;
			const _T v = ws.u[t];
			if (!(v == _T(0))) {
				u_ind[base + cnt] = P[t];
				u_val[base + cnt] = v;
				++cnt;
			}
		}
		row_count[i] = static_cast<_Index>(cnt);
		d_hat[i] = dh;
		row_lifts[i] = nl;
		row_dropped[i] = dropped;
		row_restore_skipped[i] = skipped;
	}

	// -----------------------------------------------------------------------
	// fsai_ordering_perm_: SFINAE-split ordering resolution (SLU-HK1 D-B),
	// same pattern as spmats_lu_extract_detail::dispatch_sparse_lu_extract_.
	// The shared SLU ordering routines are signed-Index only (static_assert
	// in tsparse_sparse_lu_ordering_impl.hpp), so the switch below must never
	// be instantiated in an unsigned-Index TU (the whole TU would fail to
	// build).  The runtime entry guard in policy_fsai_with_info_impl reports
	// the honest failure (fsai_status::invalid_input) before this helper can
	// be reached, so the unsigned overload is compile-only plumbing.  Returns
	// false for an unknown enum value (-> invalid_input at the caller,
	// unchanged semantics for signed Index).
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline bool
	fsai_ordering_perm_(const spmats<_T, _Index>& A,
	                    const sparse_chol_ordering ord,
	                    const _Index n,
	                    std::vector<_Index>& perm)
	{
		const spmats<_T, _Index> Acsc = A.as_csc();
		switch (ord) {
		case sparse_chol_ordering::rcm:
			perm = sparse_lu_rcm_ordering(n, Acsc.outer_index(), Acsc.inner_index());
			return true;
		case sparse_chol_ordering::amd:
			perm = sparse_lu_amd_ordering(n, Acsc.outer_index(), Acsc.inner_index());
			return true;
		case sparse_chol_ordering::nested_dissection:
			perm = sparse_lu_nested_dissection_ordering(n, Acsc.outer_index(), Acsc.inner_index());
			return true;
		default:
			// unknown enum value (natural / auto_select handled by the caller)
			return false;
		}
	}

} // namespace spmats_fsai_detail

// ---------------------------------------------------------------------------
// policy_fsai_with_info: NVI outer (non-virtual).  Finalize guarantee (same
// auto-finalize as the existing NVI outers) + squareness entry check --
// NON-throwing: a non-square input is invalid_input (design SS1.3).  Must
// never be overridden -- override policy_fsai_with_info_impl instead.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_result<_T, _Index> spmats<_T, _Index>::policy_fsai_with_info(
	spmats<_T, _Index>& U,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const fsai_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize()) {
		U.resize(_Index(0), _Index(0));
		D.resize(_Index(0), _Index(0));
		perm.clear();
		fsai_result<_T, _Index> out;
		out.status = fsai_status::invalid_input;
		return out;
	}
	return policy_fsai_with_info_impl(U, D, perm, opt);
}

// ---------------------------------------------------------------------------
// policy_fsai_with_info_impl: virtual algorithm body (default: static FSAI,
// [JFSG15] Algorithms 1+2+5, sqrt-free triple output).  Numerical events
// (tiny pivots, non-SPD, singular input) are lifted / skipped silently
// (F-D3/F-D4) and never a failure status; correctness is policed 100% by
// the downstream gate (F-D1).  Misuse (vcp::error) keeps its throwing
// contract; the final std::exception net maps to internal_error (P3
// pattern).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_result<_T, _Index> spmats<_T, _Index>::policy_fsai_with_info_impl(
	spmats<_T, _Index>& U,
	spmats<_T, _Index>& D,
	std::vector<_Index>& perm,
	const fsai_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	fsai_result<_T, _Index> out;
	try {

		const _Index n = A.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);

		if (n == _Index(0)) {
			const std::vector<_Index> ep(1u, _Index(0));
			U.assign_csc(_Index(0), _Index(0), ep, std::vector<_Index>(), std::vector<_T>());
			D.assign_csc(_Index(0), _Index(0), ep, std::vector<_Index>(), std::vector<_T>());
			perm.clear();
			out.status = fsai_status::success;
			out.nnz_U = _Index(0);
			out.nnz_pattern = _Index(0);
			out.n_pivot_modifications = _Index(0);
			out.first_modified_row = _Index(-1);
			out.n_postfilter_dropped = _Index(0);
			out.n_postfilter_restore_skipped = _Index(0);
			return out;
		}

		// --- ordering (F2-1).  auto_select resolves to natural (FSAI's
		// conservative default, design OPEN-F3 -- NOT chol's auto->amd;
		// recorded as OPEN-F13).  The ordering functions are the shared SLU
		// pattern-only routines (chol SS5.5 reuse); they symmetrize every
		// off-diagonal edge internally, so the full CSC pattern of A is a
		// valid direct input (for symmetric A this is the same graph as
		// chol's lower-triangle feed).
		sparse_chol_ordering ord = opt.ordering;
		if (ord == sparse_chol_ordering::auto_select) {
			ord = sparse_chol_ordering::natural;
		}
		perm.resize(un);
		for (std::size_t k = 0; k < un; ++k) perm[k] = static_cast<_Index>(k);
		if (ord != sparse_chol_ordering::natural) {
			// SFINAE-split helper (SLU-HK1 D-B): identical calls in identical
			// order for signed Index; false = unknown enum value (natural /
			// auto_select handled above)
			if (!spmats_fsai_detail::fsai_ordering_perm_<_T, _Index>(A, ord, n, perm)) {
				U.resize(_Index(0), _Index(0));
				D.resize(_Index(0), _Index(0));
				perm.clear();
				out = fsai_result<_T, _Index>();
				out.status = fsai_status::invalid_input;
				return out;
			}
		}
		// bijection-verification firewall (chol SS7.5-1 precedent)
		if (!spmats_fsai_detail::fsai_perm_is_bijection_<_Index>(perm, n)) {
			U.resize(_Index(0), _Index(0));
			D.resize(_Index(0), _Index(0));
			perm.clear();
			out = fsai_result<_T, _Index>();
			out.status = fsai_status::internal_error;
			return out;
		}

		// --- A_p = P^T A P, explicit (F2-1).  natural keeps A itself.
		spmats<_T, _Index> Ap_store;
		const bool ident = spmats_fsai_detail::fsai_perm_is_identity_<_Index>(perm);
		if (!ident) {
			spmats_fsai_detail::fsai_build_permuted_copy_<_T, _Index>(A, perm, Ap_store);
		}
		const spmats<_T, _Index>& Ap = ident ? A : Ap_store;
		const spmats<_T, _Index> Apc = Ap.as_csc();
		const spmats<_T, _Index> Apr = Ap.as_csr();

		// --- pattern generation (F2-2/3/4): prefilter + power recursion
		std::vector<_Index> arp, ari;
		spmats_fsai_detail::fsai_prefilter_rows_<_T, _Index>(
		    Apr, opt.prefilter_tolerance, arp, ari);
		std::vector<_Index> s_ptr, s_ind;
		spmats_fsai_detail::fsai_pattern_<_T, _Index>(
		    n, arp, ari, opt.pattern_power, opt.max_density,
		    static_cast<_Index>(Apc.values().size()), s_ptr, s_ind);
		const std::size_t nnz_pattern = s_ind.size();

		// --- pass 2: per-row solves (F3/F4).  Writes go only to the row's
		// own slice fixed by s_ptr (pass 1) and to per-row cells -- see the
		// bit-identity contract in the file header.
		std::vector<_Index> u_ind(nnz_pattern);
		std::vector<_T>     u_val(nnz_pattern);
		std::vector<_Index> row_count(un, _Index(0));
		std::vector<_T>     d_hat(un, _T(0));
		std::vector<_Index> row_lifts(un, _Index(0));
		std::vector<_Index> row_dropped(un, _Index(0));
		std::vector<_Index> row_restore_skipped(un, _Index(0));
		{
			const std::vector<_Index>& ap = Apc.outer_index();
			const std::vector<_Index>& ai = Apc.inner_index();
			const std::vector<_T>&     av = Apc.values();
			const long long n_ll = static_cast<long long>(un);
#if VCP_SPMATS_USE_OPENMP
			// design SS6.2/6.3: rows are fully independent; each thread owns
			// one workspace (declared inside the parallel region, NOT
			// firstprivate) and writes only its rows' own slices, so the
			// parallel run is bit-identical to the serial one regardless of
			// thread count and schedule.  Small problems stay serial
			// (n >= 1024 threshold, provisional -- design v0 OPEN-F7).
			#pragma omp parallel if(n_ll >= 1024)
			{
				spmats_fsai_detail::fsai_row_workspace<_T, _Index> ws;
				#pragma omp for schedule(dynamic)
				for (long long i = 0; i < n_ll; ++i) {
					spmats_fsai_detail::fsai_process_row_<_T, _Index>(
					    static_cast<std::size_t>(i), ap, ai, av, s_ptr, s_ind,
					    opt, ws, u_ind, u_val, row_count, d_hat,
					    row_lifts, row_dropped, row_restore_skipped);
				}
			}
#else
			spmats_fsai_detail::fsai_row_workspace<_T, _Index> ws;
			for (long long i = 0; i < n_ll; ++i) {
				spmats_fsai_detail::fsai_process_row_<_T, _Index>(
				    static_cast<std::size_t>(i), ap, ai, av, s_ptr, s_ind,
				    opt, ws, u_ind, u_val, row_count, d_hat,
				    row_lifts, row_dropped, row_restore_skipped);
			}
#endif
		}

		// --- pass 3 (serial): deterministic assembly in row order
		std::vector<_Index> ucp(un + 1u, _Index(0));
		for (std::size_t i = 0; i < un; ++i) {
			ucp[i + 1u] = ucp[i] + row_count[i];
		}
		const std::size_t nnz_u = static_cast<std::size_t>(ucp[un]);
		std::vector<_Index> uri(nnz_u);
		std::vector<_T>     uva(nnz_u);
		for (std::size_t i = 0; i < un; ++i) {
			const std::size_t src = static_cast<std::size_t>(s_ptr[i]);
			const std::size_t dst = static_cast<std::size_t>(ucp[i]);
			const std::size_t cnt = static_cast<std::size_t>(row_count[i]);
			for (std::size_t t = 0; t < cnt; ++t) {
				uri[dst + t] = u_ind[src + t];
				uva[dst + t] = u_val[src + t];
			}
		}
		U.assign_csc(n, n, ucp, uri, uva);

		// Dhat: diagonal store, exact zeros unstored (AINV precedent: such a
		// D is then rejected by the apply / estimate validation as a
		// structurally zero diagonal -- possible only through overflow /
		// underflow of 1/d, since d itself is lifted nonzero)
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
		out.nnz_pattern = static_cast<_Index>(nnz_pattern);
		out.n_pivot_modifications = _Index(0);
		out.first_modified_row = _Index(-1);
		out.n_postfilter_dropped = _Index(0);
		out.n_postfilter_restore_skipped = _Index(0);
		for (std::size_t i = 0; i < un; ++i) {
			out.n_pivot_modifications = out.n_pivot_modifications + row_lifts[i];
			if (row_lifts[i] > _Index(0) && out.first_modified_row < _Index(0)) {
				out.first_modified_row = static_cast<_Index>(i);
			}
			out.n_postfilter_dropped = out.n_postfilter_dropped + row_dropped[i];
			out.n_postfilter_restore_skipped =
			    out.n_postfilter_restore_skipped + row_restore_skipped[i];
		}
		return out;
	} catch (const vcp::error&) {
		// misuse / state errors keep their throwing contract (unchanged)
		throw;
	} catch (const std::exception&) {
		out = fsai_result<_T, _Index>();
		out.status = fsai_status::internal_error;
		return out;
	}
}

// ---------------------------------------------------------------------------
// policy_fsai_residual_norm_estimate: NON-GUARANTEED estimate of
// ||I - R A||_inf with R = P U Dhat^{-1} U^T P^T, computed by DELEGATION to
// the AINV helper with Z = W = U (design SS2.3, F-D10).  Permutation
// invariance (2-line identity):
//   I - R A = I - P (U Dhat^{-1} U^T) P^T A = P (I - U Dhat^{-1} U^T A_p) P^T
//   with A_p = P^T A P, and ||P X P^T||_inf = ||X||_inf,
// so the estimate equals the AINV estimate of (U, U, Dhat) against A_p.
// natural (identity perm) delegates against A itself; otherwise a permuted
// copy A_p is built once (nnz(A)-sized, design OPEN-F9 resolution).
// Single non-virtual method; consumer-type checks inline, non-throwing.
// (a)-type: A is *this.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_status spmats<_T, _Index>::policy_fsai_residual_norm_estimate(
	const spmats<_T, _Index>& U,
	const spmats<_T, _Index>& D,
	const std::vector<_Index>& perm,
	typename vcp::tsparse_scalar::real_type<_T>::type& est) const
{
	const spmats<_T, _Index>& A = *this;   // (a)-type: subject is *this
	try {
		if (!A.is_finalized()) A.finalize();
		if (!U.is_finalized()) U.finalize();
		if (!D.is_finalized()) D.finalize();
		const _Index n = U.rowsize();
		if (A.rowsize() != A.columnsize() || A.rowsize() != n) {
			return fsai_status::invalid_input;
		}
		if (!spmats_fsai_detail::fsai_perm_is_bijection_<_Index>(perm, n)) {
			return fsai_status::invalid_input;
		}
		if (spmats_fsai_detail::fsai_perm_is_identity_<_Index>(perm)) {
			return spmats_fsai_detail::fsai_from_ainv_status_(
			    A.policy_ainv_residual_norm_estimate(U, U, D, est));
		}
		spmats<_T, _Index> Ap;
		spmats_fsai_detail::fsai_build_permuted_copy_<_T, _Index>(A, perm, Ap);
		return spmats_fsai_detail::fsai_from_ainv_status_(
		    Ap.policy_ainv_residual_norm_estimate(U, U, D, est));
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		return fsai_status::internal_error;
	}
}

// ---------------------------------------------------------------------------
// policy_fsai_apply: z = R r = P (U Dhat^{-1} U^T) P^T r by DELEGATION to
// policy_ainv_apply with Z = W = U (design SS2.3):
//   r_p = P^T r  (r_p[k] = r[perm[k]], since P(perm[k], k) = 1)
//   z_p = U Dhat^{-1} U^T r_p   (AINV three-stage apply)
//   z = P z_p    (z[perm[k]] = z_p[k])
// The permutation wrappers are O(n) copies (design OPEN-F9: no permuted A
// is needed here).  Single non-virtual method; consumer-type checks inline,
// non-throwing.  z is a valid output only on success (empty otherwise).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
fsai_status spmats<_T, _Index>::policy_fsai_apply(
	const spmats<_T, _Index>& U,
	const spmats<_T, _Index>& D,
	const std::vector<_Index>& perm,
	const std::vector<_T>& r, std::vector<_T>& z) const
{
	try {
		if (!U.is_finalized()) U.finalize();
		if (!D.is_finalized()) D.finalize();
		const _Index n = U.rowsize();
		const std::size_t un = static_cast<std::size_t>(n);
		if (!spmats_fsai_detail::fsai_perm_is_bijection_<_Index>(perm, n) ||
		    r.size() != un) {
			z.clear();
			return fsai_status::invalid_input;
		}
		if (spmats_fsai_detail::fsai_perm_is_identity_<_Index>(perm)) {
			return spmats_fsai_detail::fsai_from_ainv_status_(
			    policy_ainv_apply(U, U, D, r, z));
		}
		std::vector<_T> rp(un, _T(0)), zp;
		for (std::size_t k = 0; k < un; ++k) {
			rp[k] = r[static_cast<std::size_t>(perm[k])];
		}
		const ainv_status st = policy_ainv_apply(U, U, D, rp, zp);
		if (st != ainv_status::success) {
			z.clear();
			return spmats_fsai_detail::fsai_from_ainv_status_(st);
		}
		z.assign(un, _T(0));
		for (std::size_t k = 0; k < un; ++k) {
			z[static_cast<std::size_t>(perm[k])] = zp[k];
		}
		return fsai_status::success;
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		z.clear();
		return fsai_status::internal_error;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_FSAI_IMPL_HPP
