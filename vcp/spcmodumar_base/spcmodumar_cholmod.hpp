// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spcmodumar_base/spcmodumar_cholmod.hpp
// CHOLMOD delegation layer for the sparse LL^T (Cholesky) factorization
// (SPCM campaign; design: sandbox/docs/design/spcmodumar_design_v1.md v1.0).
//
// This is the ONLY header in the repository that includes the external
// <suitesparse/cholmod.h> (dependency isolation, design SS3 / SS7-3; the
// same discipline as spumar_base/spumar_umfpack.hpp for umfpack.h).  It may
// include vcp core headers, but no vcp core header may include it.
//
// Contents:
//   - RAII guards: cholmod_common start/finish, cholmod_factor and
//     cholmod_sparse lifetime (freed on every path, design SS2-8),
//   - the normative enum -> CHOLMOD-constant mapping tables of design SS4
//     (ordering / method / failure-status),
//   - the delegation body called by spcmodumar::policy_chol_with_info_impl
//     (SPCM-1).
//
// Environment facts recorded during SPCM-0 (2026-07-17, this machine):
//   - CHOLMOD 5.2.0 (libsuitesparse-dev 1:7.6.1, libcholmod.so.5).
//   - cholmod.h of this version has NO CHOLMOD_HAS_PARTITION macro; the
//     Partition-module detection required by design D-5 is RUNTIME detection:
//     "These routines still exist if CHOLMOD is compiled with -DNPARTITION,
//     but they return Common->status = CHOLMOD_NOT_INSTALLED in that case."
//     (cholmod.h, Partition module section).  Accordingly, a NESDIS request
//     on a Partition-less build is detected by CHOLMOD_NOT_INSTALLED from
//     the analyze stage and mapped to invalid_options (D-5's explicit,
//     non-#error refusal).

#pragma once

#ifndef VCP_SPCMODUMAR_CHOLMOD_HPP
#define VCP_SPCMODUMAR_CHOLMOD_HPP

#include <suitesparse/cholmod.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include <vcp/spmats.hpp>

namespace vcp {
namespace spcmodumar_detail {

	// =====================================================================
	// RAII guards (design SS2-8: cholmod_finish / free are guaranteed on
	// every path by scope, never by manual calls on each early return).
	// =====================================================================

	// cholmod_common start/finish guard.  int32 API (cholmod_start), D-13.
	class cholmod_env {
	public:
		cholmod_env() : started_(false) {
			started_ = (cholmod_start(&common_) == 1);
			common_.print = 0;   // SPCM-F1: suppress CHOLMOD console messages -- failures are
			                     // reported through the VCP status mapping (D-11), and the
			                     // base policy is silent on npd; keep observable behaviour
			                     // aligned across policies.
		}
		~cholmod_env() {
			if (started_) cholmod_finish(&common_);
		}
		cholmod_env(const cholmod_env&) = delete;
		cholmod_env& operator=(const cholmod_env&) = delete;

		bool ok() const { return started_; }
		cholmod_common* common() { return &common_; }

	private:
		cholmod_common common_;
		bool started_;
	};

	// cholmod_factor lifetime guard.
	class cholmod_factor_guard {
	public:
		explicit cholmod_factor_guard(cholmod_common* c)
			: factor_(0), common_(c) {}
		~cholmod_factor_guard() { reset(); }
		cholmod_factor_guard(const cholmod_factor_guard&) = delete;
		cholmod_factor_guard& operator=(const cholmod_factor_guard&) = delete;

		void adopt(cholmod_factor* f) { reset(); factor_ = f; }
		cholmod_factor* get() const { return factor_; }
		void reset() {
			if (factor_ != 0) {
				cholmod_free_factor(&factor_, common_);
				factor_ = 0;
			}
		}

	private:
		cholmod_factor* factor_;
		cholmod_common* common_;
	};

	// cholmod_sparse lifetime guard (for cholmod_factor_to_sparse output).
	class cholmod_sparse_guard {
	public:
		explicit cholmod_sparse_guard(cholmod_common* c)
			: sparse_(0), common_(c) {}
		~cholmod_sparse_guard() { reset(); }
		cholmod_sparse_guard(const cholmod_sparse_guard&) = delete;
		cholmod_sparse_guard& operator=(const cholmod_sparse_guard&) = delete;

		void adopt(cholmod_sparse* s) { reset(); sparse_ = s; }
		cholmod_sparse* get() const { return sparse_; }
		void reset() {
			if (sparse_ != 0) {
				cholmod_free_sparse(&sparse_, common_);
				sparse_ = 0;
			}
		}

	private:
		cholmod_sparse* sparse_;
		cholmod_common* common_;
	};

	// =====================================================================
	// Normative mapping tables (design SS4).  These are the ONLY place the
	// VCP<->CHOLMOD enum correspondence is defined (G4: single definition).
	// =====================================================================

	// VCP ordering -> Common->method[0].ordering constant.
	//   natural            -> CHOLMOD_NATURAL                      (D-6)
	//   amd / auto_select  -> CHOLMOD_AMD  (auto resolved to amd BEFORE this
	//                         map by the caller; D-3, pure AMD, nmethods=1)
	//   rcm                -> CHOLMOD_GIVEN (user perm = VCP rcm via
	//                         cholmod_analyze_p; D-4)
	//   nested_dissection  -> CHOLMOD_NESDIS                        (D-5)
	// Returns false for an unresolved / unknown enum (caller must have
	// resolved auto_select already; a false return is a caller bug ->
	// internal_error, same convention as the native kernel dispatch).
	inline bool chol_ordering_to_cholmod_(
		const sparse_chol_ordering ordering, int& cholmod_ordering)
	{
		switch (ordering) {
		case sparse_chol_ordering::natural:
			cholmod_ordering = CHOLMOD_NATURAL; return true;
		case sparse_chol_ordering::amd:
			cholmod_ordering = CHOLMOD_AMD; return true;
		case sparse_chol_ordering::rcm:
			cholmod_ordering = CHOLMOD_GIVEN; return true;
		case sparse_chol_ordering::nested_dissection:
			cholmod_ordering = CHOLMOD_NESDIS; return true;
		default:
			return false;   // auto_select unresolved / unknown
		}
	}

	// VCP method -> Common->supernodal constant (D-8).
	//   auto_select          -> CHOLMOD_AUTO
	//   simplicial_uplooking -> CHOLMOD_SIMPLICIAL (forced)
	// Every other input value returns false and is rejected as
	// invalid_options by the caller.  In particular the enum value
	// `supernodal` (the D-9 pure addition) is a method_used RECORDING value
	// only -- it is NOT accepted as an input method (design SS4 lists
	// exactly two input methods; unmapped input -> invalid_options, the
	// same convention as the native kernel's default branch).
	inline bool chol_method_to_cholmod_(
		const sparse_chol_method method, int& cholmod_supernodal)
	{
		switch (method) {
		case sparse_chol_method::auto_select:
			cholmod_supernodal = CHOLMOD_AUTO; return true;
		case sparse_chol_method::simplicial_uplooking:
			cholmod_supernodal = CHOLMOD_SIMPLICIAL; return true;
		default:
			return false;   // incl. supernodal-as-input (D-9)
		}
	}

	// CHOLMOD outcome -> VCP status (D-11 / design SS4 second table).
	// To be called AFTER cholmod_factorize with the final Common->status and
	// L->minor.  n / minor are size_t on the CHOLMOD side.
	//   status == CHOLMOD_OK, minor == n          -> success
	//   status == CHOLMOD_NOT_POSDEF, minor <  n  -> not_positive_definite,
	//                                                failure_at = minor
	//   anything else (negative status: memory, invalid, not-installed...)
	//                                             -> internal_error
	// The D-5 refinement (NESDIS on a Partition-less build ->
	// CHOLMOD_NOT_INSTALLED -> invalid_options) is applied by the caller at
	// the analyze stage, before this map is reached.
	// inconclusive_pivot_test is NEVER produced on the delegation path
	// (double-only, D-11; stated in the contract).
	inline sparse_chol_status map_cholmod_outcome_(
		const int common_status, const std::size_t minor, const std::size_t n,
		int& failure_at)
	{
		failure_at = -1;
		if (common_status == CHOLMOD_OK && minor == n) {
			return sparse_chol_status::success;
		}
		if (common_status == CHOLMOD_NOT_POSDEF && minor < n) {
			failure_at = static_cast<int>(minor);
			return sparse_chol_status::not_positive_definite;
		}
		return sparse_chol_status::internal_error;
	}

	// =====================================================================
	// Delegation body (SPCM-1; design SS2, FIXED pipeline order 1..8).
	//
	// Called by spcmodumar::policy_chol_with_info_impl; A arrives finalized
	// and square (the non-virtual outer's guarantee).  Runtime failure is a
	// status; this function itself throws only through the vcp core helpers
	// it calls (misuse contract) -- the caller wraps it in the P3 net.
	//
	// Self-contained correctness (design G4): the delegate is not trusted --
	//   * the VCP certified symmetry check runs BEFORE delegation (D-10),
	//   * the rcm permutation is bijection-verified before sending AND the
	//     L->Perm read-back is bijection-verified after (D-4),
	//   * the returned factor is structurally checked (lower-triangular,
	//     ascending rows) at the boundary copy.
	// =====================================================================
	inline chol_result<double, int> cholmod_chol_delegate_(
		const spmats<double, int>& A,
		spmats<double, int>& L,
		std::vector<int>& perm,
		const chol_options<double>& opt)
	{
		typedef int Index;
		chol_result<double, Index> out;   // status defaults to internal_error
		out.ordering_used = opt.ordering;
		out.method_used   = opt.method;

		// out parameters are valid only on success (D-3, base contract)
		L.resize(0, 0);
		perm.clear();

		const Index n = A.rowsize();

		// ---- 1. options check -> invalid_options (design SS2-1) ----------
		// pd_tol != 0: CHOLMOD has no corresponding knob; rejected explicitly
		// instead of silently ignored (D-7).  The default 0 is accepted.
		if (opt.pd_tol != 0.0) {
			out.status = sparse_chol_status::invalid_options;
			return out;
		}
		int cm_supernodal = CHOLMOD_AUTO;
		if (!chol_method_to_cholmod_(opt.method, cm_supernodal)) {
			// unmapped input method, incl. `supernodal` (recording-only
			// value, D-9) and unknown enum values.
			out.status = sparse_chol_status::invalid_options;
			return out;
		}
		sparse_chol_ordering resolved = opt.ordering;
		switch (opt.ordering) {
		case sparse_chol_ordering::auto_select:
			resolved = sparse_chol_ordering::amd;   // D-3: pure AMD, reported
			break;
		case sparse_chol_ordering::natural:
		case sparse_chol_ordering::rcm:
		case sparse_chol_ordering::amd:
		case sparse_chol_ordering::nested_dissection:
			break;
		default:
			out.status = sparse_chol_status::invalid_options;
			return out;
		}
		out.ordering_used = resolved;
		// D-5 (Partition absence -> invalid_options) is a RUNTIME detection
		// in CHOLMOD 5.x (no CHOLMOD_HAS_PARTITION macro): it is applied at
		// the analyze stage below via CHOLMOD_NOT_INSTALLED.

		// ---- 2. input validation, CSC form, int32 guard (SS2-2, D-13) ----
		// spmats -> full CSC via the same conversion helper as the native
		// dispatch (spmats_chol_impl.hpp).  NOTE (recorded): the
		// spumar_convert.hpp asset (umfpack_csc_view) is an UMFPACK-shaped
		// full-CSC view without the lower-triangle/stype=-1 form needed
		// here, so the reused asset is sparse_lu_make_csc_storage instead.
		const csc_storage<double, Index> C = sparse_lu_make_csc_storage(A);
		if (!sparse_chol_detail::sparse_chol_validate_csc_<double, Index>(
			n, C.col_ptr, C.row_ind, C.values)) {
			out.status = sparse_chol_status::invalid_input;
			return out;
		}
		if (C.values.size() > static_cast<std::size_t>(
			std::numeric_limits<int>::max())) {
			// int32 API guard (D-13): n fits by construction (Index = int,
			// checked >= 0 above); nnz is guarded here.
			out.status = sparse_chol_status::invalid_input;
			return out;
		}

		// ---- 3. VCP certified symmetry check (SS2-3, D-10) ---------------
		if (opt.check_symmetry) {
			if (!sparse_chol_detail::sparse_chol_symmetry_certified_<double, Index>(
				n, C.col_ptr, C.row_ind, C.values, opt.symmetry_tol)) {
				out.status = sparse_chol_status::not_symmetric;
				return out;
			}
		}

		// ---- 4. structural-empty scan (SS2-4; same status mapping as the
		//         native kernel: not_positive_definite + structural_empty_at,
		//         failure_at stays -1) ------------------------------------
		{
			const Index empty_at =
				sparse_chol_detail::sparse_chol_scan_structural_empty_<Index>(
					n, C.col_ptr, C.row_ind);
			if (empty_at >= Index(0)) {
				out.status = sparse_chol_status::not_positive_definite;
				out.structural_empty_at = empty_at;
				return out;
			}
		}

		// ---- trivial n == 0 short-circuit --------------------------------
		// CHOLMOD 5.2 rejects an empty view with NULL value array ("invalid
		// xtype or dtype", measured during SPCM-1), so the empty problem is
		// answered before delegation: success with empty outputs, method_used
		// = simplicial_uplooking (base-kernel parity: the native pipeline
		// reports the resolved simplicial method for n = 0).
		if (n == 0) {
			std::vector<Index> zp(1, Index(0));
			L.assign_csc(0, 0, zp, std::vector<Index>(), std::vector<double>());
			out.method_used = sparse_chol_method::simplicial_uplooking;
			out.nnz_L = 0;
			out.status = sparse_chol_status::success;
			return out;
		}

		// ---- 5. ordering preparation (SS2-5) -----------------------------
		int cm_ordering = CHOLMOD_NATURAL;
		if (!chol_ordering_to_cholmod_(resolved, cm_ordering)) {
			out.status = sparse_chol_status::internal_error;   // caller bug
			return out;
		}
		std::vector<Index> rcm_perm;
		if (resolved == sparse_chol_ordering::rcm) {
			// VCP rcm on the SAME full pattern the native kernel uses
			// (SS5.5: the ordering builder symmetrizes internally), then the
			// send-side bijection verification (D-4).
			rcm_perm = sparse_lu_rcm_ordering<Index>(n, C.col_ptr, C.row_ind);
			if (!sparse_chol_detail::sparse_chol_verify_permutation_<Index>(
				n, rcm_perm)) {
				out.status = sparse_chol_status::internal_error;
				return out;
			}
		}

		// ---- 6. CHOLMOD delegation (SS2-6) -------------------------------
		// Lower-triangle CSC (r >= c only; strictly-upper stored entries are
		// ignored, same read scope as the native numeric stage, SS1.2) for
		// the stype = -1 view (D-13).
		std::vector<int>    lo_ptr(static_cast<std::size_t>(n) + 1u, 0);
		std::vector<int>    lo_ind;
		std::vector<double> lo_val;
		{
			lo_ind.reserve(C.row_ind.size());
			lo_val.reserve(C.values.size());
			const std::size_t un = static_cast<std::size_t>(n);
			for (std::size_t c = 0; c < un; ++c) {
				for (Index k = C.col_ptr[c]; k < C.col_ptr[c + 1u]; ++k) {
					const Index r = C.row_ind[static_cast<std::size_t>(k)];
					if (r >= static_cast<Index>(c)) {
						lo_ind.push_back(r);
						lo_val.push_back(C.values[static_cast<std::size_t>(k)]);
					}
				}
				lo_ptr[c + 1u] = static_cast<int>(lo_ind.size());
			}
		}

		cholmod_env env;
		if (!env.ok()) {
			out.status = sparse_chol_status::internal_error;
			return out;
		}
		cholmod_common* cm = env.common();
		cm->final_ll   = 1;               // D-12: LL' forced (simplicial LDL'
		                                  // converted; supernodal is LL' natively)
		cm->supernodal = cm_supernodal;   // D-8
		cm->nmethods   = 1;               // single method, no fallback chain
		cm->method[0].ordering = cm_ordering;

		cholmod_sparse Aview;             // non-owning view, int32/real/double
		Aview.nrow  = static_cast<std::size_t>(n);
		Aview.ncol  = static_cast<std::size_t>(n);
		Aview.nzmax = lo_ind.size();
		Aview.p = static_cast<void*>(lo_ptr.data());
		Aview.i = static_cast<void*>(lo_ind.empty() ? static_cast<int*>(0) : lo_ind.data());
		Aview.nz = 0;
		Aview.x = static_cast<void*>(lo_val.empty() ? static_cast<double*>(0) : lo_val.data());
		Aview.z = 0;
		Aview.stype  = -1;                // symmetric, LOWER stored (D-13)
		Aview.itype  = CHOLMOD_INT;
		Aview.xtype  = CHOLMOD_REAL;
		Aview.dtype  = CHOLMOD_DOUBLE;
		Aview.sorted = 1;                 // ascending rows per column (CSC)
		Aview.packed = 1;

		cholmod_factor_guard Lf(cm);
		if (resolved == sparse_chol_ordering::rcm) {
			Lf.adopt(cholmod_analyze_p(&Aview, rcm_perm.data(),
				static_cast<int*>(0), 0, cm));
		} else {
			Lf.adopt(cholmod_analyze(&Aview, cm));
		}
		if (Lf.get() == 0 || cm->status < CHOLMOD_OK) {
			if (resolved == sparse_chol_ordering::nested_dissection &&
				cm->status == CHOLMOD_NOT_INSTALLED) {
				// D-5: NESDIS requested on a Partition-less CHOLMOD build --
				// explicit, notifying refusal (never #error).
				out.status = sparse_chol_status::invalid_options;
				return out;
			}
			out.status = sparse_chol_status::internal_error;
			return out;
		}

		const int fac_rc = cholmod_factorize(&Aview, Lf.get(), cm);

		// method_used: what ACTUALLY ran (D-8/D-9), from L->is_super.  Read
		// right after factorize -- the extraction below converts the factor
		// to simplicial and would destroy this bit.
		out.method_used = (Lf.get()->is_super != 0)
			? sparse_chol_method::supernodal
			: sparse_chol_method::simplicial_uplooking;

		// ---- 7. status / minor mapping (SS2-7, D-11) ---------------------
		int failure_at = -1;
		const sparse_chol_status mapped = map_cholmod_outcome_(
			cm->status, Lf.get()->minor, static_cast<std::size_t>(n), failure_at);
		if (mapped != sparse_chol_status::success) {
			out.status = mapped;
			out.failure_at = failure_at;
			return out;
		}
		if (fac_rc == 0) {
			// defensive: FALSE return with an OK-looking Common->status must
			// never be declared a success.
			out.status = sparse_chol_status::internal_error;
			return out;
		}

		// ---- success: perm read-back + verification (D-4) ----------------
		std::vector<Index> perm_out(static_cast<std::size_t>(n));
		{
			const int* p = static_cast<const int*>(Lf.get()->Perm);
			if (p == 0 && n > 0) {
				out.status = sparse_chol_status::internal_error;
				return out;
			}
			for (std::size_t i = 0; i < static_cast<std::size_t>(n); ++i) {
				perm_out[i] = p[i];
			}
		}
		if (!sparse_chol_detail::sparse_chol_verify_permutation_<Index>(n, perm_out)) {
			out.status = sparse_chol_status::internal_error;
			return out;
		}

		// ---- factor extraction (D-12): simplicial / packed / monotonic
		// LL', then CSC copy-out ------------------------------------------
		if (cholmod_change_factor(CHOLMOD_REAL, /*to_ll*/1, /*to_super*/0,
			/*to_packed*/1, /*to_monotonic*/1, Lf.get(), cm) == 0 ||
			cm->status < CHOLMOD_OK) {
			out.status = sparse_chol_status::internal_error;
			return out;
		}
		cholmod_sparse_guard Ls(cm);
		Ls.adopt(cholmod_factor_to_sparse(Lf.get(), cm));
		if (Ls.get() == 0 || cm->status < CHOLMOD_OK ||
			Ls.get()->nrow != static_cast<std::size_t>(n) ||
			Ls.get()->ncol != static_cast<std::size_t>(n) ||
			Ls.get()->itype != CHOLMOD_INT ||
			Ls.get()->xtype != CHOLMOD_REAL ||
			Ls.get()->dtype != CHOLMOD_DOUBLE ||
			Ls.get()->packed == 0) {
			out.status = sparse_chol_status::internal_error;
			return out;
		}

		// Boundary copy with the same rules as the native boundary layer:
		// lower-triangular check, ascending rows per column (sorted here if
		// the backend returned them unsorted), certified zero (e == 0.0)
		// NOT stored (spmats explicit-zero invariant; materialized stored
		// count <= backend stored count).
		{
			const int*    sp = static_cast<const int*>(Ls.get()->p);
			const int*    si = static_cast<const int*>(Ls.get()->i);
			const double* sx = static_cast<const double*>(Ls.get()->x);
			const std::size_t un = static_cast<std::size_t>(n);

			out.nnz_L = sp[un];   // backend stored count (pre zero-drop)

			std::vector<Index>  l_ptr(un + 1u, Index(0));
			std::vector<Index>  l_ind;
			std::vector<double> l_val;
			l_ind.reserve(static_cast<std::size_t>(sp[un]));
			l_val.reserve(static_cast<std::size_t>(sp[un]));
			std::vector<std::pair<int, double> > colbuf;
			for (std::size_t j = 0; j < un; ++j) {
				colbuf.clear();
				bool ascending = true;
				int prev = -1;
				for (int k = sp[j]; k < sp[j + 1u]; ++k) {
					const int r = si[static_cast<std::size_t>(k)];
					if (r < static_cast<int>(j) || r >= n) {
						// not lower triangular / out of range: never accept
						out.status = sparse_chol_status::internal_error;
						return out;
					}
					if (r <= prev) ascending = false;
					prev = r;
					colbuf.push_back(std::pair<int, double>(
						r, sx[static_cast<std::size_t>(k)]));
				}
				if (!ascending) {
					std::sort(colbuf.begin(), colbuf.end());
					for (std::size_t t = 1; t < colbuf.size(); ++t) {
						if (colbuf[t].first == colbuf[t - 1u].first) {
							out.status = sparse_chol_status::internal_error;
							return out;   // duplicate row index
						}
					}
				}
				for (std::size_t t = 0; t < colbuf.size(); ++t) {
					const double e = colbuf[t].second;
					if (e == 0.0) { /* certified zero: not stored */ }
					else {
						l_ind.push_back(colbuf[t].first);
						l_val.push_back(e);
					}
				}
				l_ptr[j + 1u] = static_cast<Index>(l_ind.size());
			}
			L.assign_csc(n, n, l_ptr, l_ind, l_val);
		}

		perm = perm_out;
		out.failure_at = -1;
		out.status = sparse_chol_status::success;
		return out;
	}

} // namespace spcmodumar_detail
} // namespace vcp

#endif // VCP_SPCMODUMAR_CHOLMOD_HPP
