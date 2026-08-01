// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_ldl_shift.hpp
// SLDL-SH (Phase B / B2): the A - sigma*B LDL^T shift-iteration handle
// (ldl_shift_handle) -- setup once, then inertia_at(sigma) [primary] /
// factor_at(sigma) [secondary] repeatedly (design SLDL-SH v1 SS2.2).
//
// TWO-PHASE HEADER (single file, two include points -- the one-new-file
// analogue of the spmats_ldl.hpp / spmats_ldl_impl.hpp pair):
//   phase 1 (included by spmats.hpp BEFORE the spmats class body): the
//     option helper, result / workspace types, the handle class and the
//     T-only detail helpers (the types appear in the NVI signatures);
//   phase 2 (included by spmats.hpp AFTER the class body): the setup
//     dispatch and the out-of-line policy method definitions, which need
//     the complete spmats<_T,_Index>.
// No "#pragma once" here by design: the second include must enter phase 2.
// Do NOT include this file directly (phase 2 needs the class body in
// between); include <vcp/spmats.hpp> or <vcp/spmatrix.hpp> instead.
//
// Concurrency contract (design SS2.4): after setup the handle is READ-ONLY
// (no mutable members, no caches).  inertia_at / factor_at take an explicit
// per-thread workspace, so different sigma values may run concurrently on
// const-shared handles; all per-call diagnostics live in the workspace and
// the return value.  A sigma-level failure (H-3) never invalidates the
// handle: every later call proceeds normally.
//
// Failure reporting is info-only (H-2): no strict variants exist in this
// module.  Misuse (non-square input, dimension mismatch, unsigned Index)
// keeps the throwing vcp::error contract of the existing policy entries.
//
// Allocation contract (STOP-1 ruling, recorded in
// sandbox/docs/issues/SLDL-SH_stop1_issue.md):
//   A-1  the B2 layer itself (plan / handle / assembled-value buffer /
//        inertia scan) performs ZERO dynamic allocation in inertia_at once
//        set up and given an explicit reused workspace;
//   A-2  the total allocation count/bytes of one inertia_at call are
//        CONSTANT over sigma and over repetition (no growth, no leak);
//   A-3  the call-local scratch inside the B1 numeric kernels is a known,
//        accepted allocation source and is out of scope for this track.

// ===========================================================================
// PHASE 1 -- types (include point: before the spmats class body)
// ===========================================================================
#ifndef VCP_SPMATS_LDL_SHIFT_HPP_TYPES
#define VCP_SPMATS_LDL_SHIFT_HPP_TYPES

#include <cstddef>
#include <exception>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_sparse_ldl.hpp>
#include <vcp/tsparse/detail/tsparse_sparse_ldl_shift_plan_impl.hpp>
#include <vcp/spmats_base/spmats_ldl.hpp>

namespace vcp {

template <typename _T, typename _Index> class spmats;

namespace spmats_ldl_shift_detail {
	// Sole writer of ldl_shift_handle private state (befriended below);
	// used only by the policy implementations of phase 2.
	struct shift_handle_access;
} // namespace spmats_ldl_shift_detail

// ---------------------------------------------------------------------------
// Handle default options (H-1): supernodal x none.  The one-shot LDL default
// (auto_select -> baseline_dynamic, pivoting bk) is NOT changed (H-6); this
// helper only seeds the handle entry points, and callers may override any
// field.  none is the only mode whose symbolic analysis is fully reusable
// across sigma (bk exchanges make the structure sigma-dependent), which is
// exactly the workload of the shift iteration.
// ---------------------------------------------------------------------------
template <typename _T>
inline ldl_options<_T> ldl_shift_default_options() {
	ldl_options<_T> opt;
	opt.method   = sparse_ldl_method::supernodal;
	opt.pivoting = sparse_ldl_pivoting::none;
	return opt;
}

// ---------------------------------------------------------------------------
// inertia_at result (implementation guide SS2-2): the existing
// inertia_result<Index> (counts + scan status + inconclusive position +
// the per-sigma numeric status in ldl_status) plus the growth diagnostic
// (D-5 conduit).  growth_log2 is valid iff growth_valid (P4).
// ---------------------------------------------------------------------------
template <typename _Index>
struct ldl_shift_inertia_result {
	inertia_result<_Index> inertia;
	int  growth_log2;
	bool growth_valid;

	ldl_shift_inertia_result()
		: inertia(), growth_log2(0), growth_valid(false) {}
};

// ---------------------------------------------------------------------------
// Per-thread workspace: the B1 numeric workspace, the assembled A - sigma*B
// value buffer and the reused tsparse result object (its vectors keep their
// capacity across calls; all per-call numeric diagnostics stay readable here
// -- design SS2.4 "diagnostics live on the workspace / return-value side").
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
struct ldl_shift_workspace {
	static_assert(std::is_signed<_Index>::value,
	              "ldl_shift_workspace requires a signed Index type");

	sparse_ldl_numeric_workspace<_T, _Index> numeric;
	std::vector<_T>                          val;   // assembled lower-CSC values
	sparse_ldl_result<_T, _Index>            res;   // reused numeric output
};

namespace spmats_ldl_shift_detail {

	// -----------------------------------------------------------------------
	// reset_result_scalars_: restore every scalar/enum diagnostic of a reused
	// sparse_ldl_result to its default-constructed value WITHOUT touching the
	// vectors (their capacity is the allocation-reuse contract).  The numeric
	// kernels assign every field they own, but the frozen baseline kernels
	// predate the workspace reuse pattern, so a stale field must never leak
	// from call k into call k+1.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline void reset_result_scalars_(sparse_ldl_result<_T, _Index>& r) {
		r.status              = sparse_ldl_status::internal_error;
		r.n_pivots_1x1        = _Index(0);
		r.n_pivots_2x2        = _Index(0);
		r.first_zero_pivot    = _Index(-1);
		r.inconclusive_at     = _Index(-1);
		r.structural_empty_at = _Index(-1);
		r.nnz_L               = _Index(0);
		r.ordering_used       = sparse_ldl_ordering::auto_select;
		r.method_used         = sparse_ldl_method::auto_select;
		r.dense_delegated     = false;
		r.pivot_mode_used     = sparse_ldl_pivoting::bk;
		r.diag_kernel_used    = sparse_ldl_kernel_used::not_applicable;
		r.n_supernodes        = _Index(0);
		r.max_supernode_width = _Index(0);
		r.n_boundary_splits   = _Index(0);
		r.nnz_L_static        = _Index(0);
		r.n_zero_skips        = _Index(0);
		r.out_of_panel_at     = _Index(-1);
		r.gemm_call_count     = 0;
		r.gemm_time_ns        = 0;
		r.growth_log2         = 0;
		r.growth_valid        = false;
	}

	// -----------------------------------------------------------------------
	// inertia_from_internal_: certified inertia scan over the INTERNAL block
	// diagonal triple (D_diag / D_sub / D_block2) of a sparse_ldl_result.
	//
	// STOP-1 secondary ruling: this is the FUNCTION-FORM verbatim mirror of
	// the counting branches of spmats::policy_inertia_from_block_diagonal
	// (spmats_base/spmats_ldl_impl.hpp, design v2 SS7.2 rules 2-5) -- kept
	// line-for-line where the representations allow, because extracting the
	// member's loop would touch a file outside this track's change scope.
	// Representation mapping (equivalence machine-checked by acceptance 1:
	// inertia_at == the existing member run on the materialized D):
	//   - 2x2 block head:  ehas[k] (stored subdiagonal)  ->  D_block2[k] != 0
	//     (the kernels set D_block2 exactly for the executed/zero-skipped 2x2
	//     pivots; a zero-skipped pair has d = e = 0 and lands on the same
	//     counts through the certified branches);
	//   - structurally empty column (member rule 4):  the internal triple
	//     always stores a value, an exact zero -- certified |0| <= tol is the
	//     same n_zero outcome (tol >= 0), so rule 4 is subsumed by rule 3;
	//   - every certified comparison, the R = real_type mapping and the
	//     inconclusive third branches mirror the member verbatim.
	// The ldl_status field of the returned inertia_result is left at its
	// default (success); the caller records the per-sigma numeric status.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline inertia_result<_Index> inertia_from_internal_(
	    const _Index n,
	    const std::vector<_T>&   D_diag,
	    const std::vector<_T>&   D_sub,
	    const std::vector<char>& D_block2,
	    const typename vcp::tsparse_scalar::real_type<_T>::type& tol)
	{
		using std::abs;
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;

		inertia_result<_Index> out;
		out.status = inertia_status::invalid_input;
		if (n < _Index(0)) return out;
		const std::size_t un = static_cast<std::size_t>(n);
		if (D_diag.size() != un || D_sub.size() != un || D_block2.size() != un) {
			return out;
		}

		_Index npos = _Index(0), nneg = _Index(0), nzero = _Index(0);
		std::size_t k = 0;
		while (k < un) {
			if (k + 1u < un && D_block2[k] != char(0)) {
				// ---- 2x2 block (k, k+1) -- member scan mirrored verbatim
				const _T dk  = D_diag[k];
				const _T dk1 = D_diag[k + 1u];
				const _T e   = D_sub[k];
				const _T det = dk * dk1 - e * e;
				const _T tr  = dk + dk1;
				const R rdet = vcp::tsparse_scalar::real_part(det);
				const R rtr  = vcp::tsparse_scalar::real_part(tr);
				const R absdet = abs(det);
				if (absdet <= tol) {
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
			} else {
				// ---- 1x1 diagonal (member rule 3; rule 4 subsumed, see above)
				const _T dk = D_diag[k];
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

} // namespace spmats_ldl_shift_detail

// ---------------------------------------------------------------------------
// ldl_shift_handle<_T,_Index> (design SS2.2)
//
// Produced by spmats::policy_ldl_shift_setup_with_info (NVI; the spmatrix
// transfer API is ldl_shift_setup_with_info).  Owns every sigma-INDEPENDENT
// stage of the shift iteration, computed once at setup:
//   spmats -> CSC conversion, certified symmetry checks (B: symmetry only,
//   H-4), the pattern-merge plan, the ordering and the full symbolic
//   analysis (const-shared B1 struct).
// inertia_at(sigma [, ws [, zero_tol]]) then performs only the
// sigma-DEPENDENT work: O(nnz) value assembly -> B1 numeric (workspace
// reuse) -> direct certified scan of the internal D triple.  No spmatrix /
// spmats L, D or perm object is ever constructed on this path.
//
// valid() / status(): setup outcome.  status() is only meaningful after a
// setup call (a default-constructed handle reports internal_error and
// valid() == false).  Setup failures reuse the one-shot status vocabulary
// (invalid_options / not_symmetric / structural_singularity / invalid_input
// / internal_error); symmetry_failed_on() tells which input failed the
// certified symmetry check (0 = none, 1 = A, 2 = B).
//
// unsigned _Index: the handle type stays INSTANTIABLE (it is the return
// type of the policy entry, whose unsigned path throws vcp::state_error at
// runtime -- same SFINAE discipline as lu_factor_handle); the signed-only
// bodies are never instantiated for unsigned _Index.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
class ldl_shift_handle {
public:
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	typedef typename std::conditional<std::is_signed<_Index>::value, _Index,
		typename std::make_signed<_Index>::type>::type factor_index_type;
	typedef ldl_shift_workspace<_T, factor_index_type>       workspace_type;
	typedef ldl_shift_inertia_result<factor_index_type>      inertia_result_type;

	ldl_shift_handle()
		: plan_(), sym_(), valA_(), valB_(), opt_(),
		  method_used_(sparse_ldl_method::auto_select),
		  ordering_used_(sparse_ldl_ordering::auto_select),
		  n_(_Index(0)), valid_(false),
		  setup_status_(sparse_ldl_status::internal_error),
		  symmetry_failed_on_(0) {}

	bool valid() const { return valid_; }
	sparse_ldl_status status() const { return setup_status_; }
	_Index n() const { return n_; }
	bool has_B() const { return plan_.has_B; }
	int symmetry_failed_on() const { return symmetry_failed_on_; }
	sparse_ldl_method   method_used() const { return method_used_; }
	sparse_ldl_ordering ordering_used() const { return ordering_used_; }
	// sigma-independent structure diagnostics (valid() only)
	factor_index_type nnz_merged() const { return plan_.nnz(); }
	const sparse_ldl_symbolic_result<factor_index_type>& symbolic() const { return sym_; }

	// ---- primary API (H-5): inertia only, no factor materialization.
	// Convenience form: call-local workspace (H-3: the handle itself never
	// holds mutable state; omitting ws trades the allocation-reuse contract
	// for convenience, implementation guide SS2-3).
	inertia_result_type inertia_at(const _T& sigma) const {
		workspace_type ws;
		return inertia_at(sigma, ws, real_type(0));
	}
	inertia_result_type inertia_at(const _T& sigma, workspace_type& ws) const {
		return inertia_at(sigma, ws, real_type(0));
	}
	// zero_tol gates the certified zero test of the D scan (same meaning as
	// inertia_options::zero_tol; default 0 = certified exact zeros only).
	inertia_result_type inertia_at(const _T& sigma, workspace_type& ws,
	                               const real_type& zero_tol) const {
		return inertia_at_impl_(sigma, ws, zero_tol);
	}

	// ---- secondary API (H-5): factor materialization (SH-2).
	// Output contract IDENTICAL to policy_ldl_with_info / ldl_with_info
	// (design v2 SS5.2/SS5.4): L unit lower (explicit unit diagonal), D
	// symmetric 1x1/2x2 block diagonal with exact zeros NOT stored, perm
	// new->old; L / D / p are valid only when the returned status is
	// success or zero_pivot, and are emptied otherwise.  Byte-identity
	// (acceptance 1): with the handle's own options the outputs are
	// byte-identical to the one-shot ldl_with_info of the merged matrix.
	// SpMat is any matrix type exposing resize/assign_csc with this
	// handle's factor_index_type (vcp::spmats and vcp::spmatrix both do).
	// Definitions live in phase 2 of this header (they reuse the one-shot
	// conduit make_ldl_result_ and need the complete spmats).
	template <class SpMat>
	ldl_result<_T, factor_index_type> factor_at(
	    const _T& sigma, SpMat& L, SpMat& D,
	    std::vector<factor_index_type>& p, workspace_type& ws) const;
	template <class SpMat>
	ldl_result<_T, factor_index_type> factor_at(
	    const _T& sigma, SpMat& L, SpMat& D,
	    std::vector<factor_index_type>& p) const {
		workspace_type ws;
		return factor_at(sigma, L, D, p, ws);
	}

	// Explicit retry (D-6 integration): re-factor the SAME handle data with
	// an explicitly chosen method x pivoting -- the documented escape route
	// after pivot_out_of_panel (typically baseline_dynamic x bk, which
	// handles exchanges at any distance).  This is the caller's EXPLICIT
	// decision; the library never switches silently.  The handle's
	// converted/merged data and ordering are reused; when the requested
	// method needs a deeper symbolic analysis than setup produced (e.g.
	// supernodal from a baseline handle), the full analysis is computed
	// locally for this call only (the handle stays unchanged).
	template <class SpMat>
	ldl_result<_T, factor_index_type> factor_at_with(
	    const _T& sigma, const sparse_ldl_method method,
	    const sparse_ldl_pivoting pivoting, SpMat& L, SpMat& D,
	    std::vector<factor_index_type>& p, workspace_type& ws) const;
	template <class SpMat>
	ldl_result<_T, factor_index_type> factor_at_with(
	    const _T& sigma, const sparse_ldl_method method,
	    const sparse_ldl_pivoting pivoting, SpMat& L, SpMat& D,
	    std::vector<factor_index_type>& p) const {
		workspace_type ws;
		return factor_at_with(sigma, method, pivoting, L, D, p, ws);
	}

private:
	// shared sigma-dependent factor body (declared here, defined in phase
	// 2): assembly -> numeric under opt2 -> L/D/p materialization.
	template <class SpMat>
	ldl_result<_T, factor_index_type> factor_core_(
	    const _T& sigma, const sparse_ldl_options<_T>& opt2,
	    const sparse_ldl_method method_used2, SpMat& L, SpMat& D,
	    std::vector<factor_index_type>& p, workspace_type& ws) const;

	// signed-Index body: assembly -> numeric (by delegation) -> direct scan.
	template <typename _I = _Index>
	typename std::enable_if<std::is_signed<_I>::value, inertia_result_type>::type
	inertia_at_impl_(const _T& sigma, workspace_type& ws,
	                 const real_type& zero_tol) const
	{
		inertia_result_type out;
		if (!valid_) {
			// info-only reporting (H-2): an unset / failed handle is a
			// controlled invalid_input, with the setup status attached.
			out.inertia.status = inertia_status::invalid_input;
			out.inertia.ldl_status = setup_status_;
			return out;
		}
		try {
			if (!plan_.build_values(sigma, valA_, valB_, ws.val)) {
				out.inertia.status = inertia_status::internal_error;
				return out;
			}
			spmats_ldl_shift_detail::reset_result_scalars_(ws.res);
			// the numeric entry expects pre-resolved method/ordering
			ws.res.method_used   = method_used_;
			ws.res.ordering_used = ordering_used_;
			sparse_ldl_factorize_numeric_with_info(
			    static_cast<factor_index_type>(n_),
			    plan_.col_ptr, plan_.row_ind, ws.val,
			    sym_, opt_, ws.numeric, ws.res);

			out.growth_log2  = ws.res.growth_log2;
			out.growth_valid = ws.res.growth_valid;
			if (ws.res.status == sparse_ldl_status::success ||
			    ws.res.status == sparse_ldl_status::zero_pivot) {
				out.inertia = spmats_ldl_shift_detail::inertia_from_internal_(
				    static_cast<factor_index_type>(n_),
				    ws.res.D_diag, ws.res.D_sub, ws.res.D_block2, zero_tol);
				out.inertia.ldl_status = ws.res.status;
			} else {
				// sigma-level failure (H-3): reported, handle stays alive.
				out.inertia.status = inertia_status::factorization_failed;
				out.inertia.ldl_status = ws.res.status;
			}
			return out;
		} catch (const std::exception&) {
			// P3 final protection net (firing = defect)
			out.inertia.status = inertia_status::internal_error;
			return out;
		}
	}

	// unsigned-Index body: unreachable through a valid handle (setup throws
	// before one is produced); split so the signed-only code above is never
	// instantiated for unsigned _Index.
	template <typename _I = _Index>
	typename std::enable_if<!std::is_signed<_I>::value, inertia_result_type>::type
	inertia_at_impl_(const _T&, workspace_type&, const real_type&) const {
		vcp::throw_error<vcp::state_error>(
		    "spmats::ldl_shift_handle::inertia_at: sparse LDL requires a signed Index type");
		return inertia_result_type();
	}

	sparse_ldl_shift_plan<_T, factor_index_type>  plan_;
	sparse_ldl_symbolic_result<factor_index_type> sym_;
	std::vector<_T>       valA_;   // CSC values of A (conversion order)
	std::vector<_T>       valB_;   // CSC values of B; empty when B omitted
	sparse_ldl_options<_T> opt_;   // frozen at setup
	sparse_ldl_method     method_used_;
	sparse_ldl_ordering   ordering_used_;
	_Index n_;
	bool valid_;
	sparse_ldl_status setup_status_;
	int symmetry_failed_on_;       // 0 = none, 1 = A, 2 = B

	friend struct spmats_ldl_shift_detail::shift_handle_access;
};

} // namespace vcp

// ===========================================================================
// PHASE 2 -- setup dispatch + policy method definitions (include point:
// after the spmats class body)
// ===========================================================================
#else
#ifndef VCP_SPMATS_LDL_SHIFT_HPP_IMPL
#define VCP_SPMATS_LDL_SHIFT_HPP_IMPL

namespace vcp {

namespace spmats_ldl_shift_detail {

	struct shift_handle_access {
		// success path: move the sigma-independent stages into the handle
		template <typename _T, typename _Index>
		static void assign(
		    ldl_shift_handle<_T, _Index>& h,
		    sparse_ldl_shift_plan<_T,
		        typename ldl_shift_handle<_T, _Index>::factor_index_type>&& plan,
		    sparse_ldl_symbolic_result<
		        typename ldl_shift_handle<_T, _Index>::factor_index_type>&& sym,
		    std::vector<_T>&& valA, std::vector<_T>&& valB,
		    const sparse_ldl_options<_T>& opt,
		    const sparse_ldl_method method_used,
		    const sparse_ldl_ordering ordering_used,
		    const _Index n)
		{
			h.plan_ = std::move(plan);
			h.sym_  = std::move(sym);
			h.valA_ = std::move(valA);
			h.valB_ = std::move(valB);
			h.opt_  = opt;
			h.method_used_   = method_used;
			h.ordering_used_ = ordering_used;
			h.n_ = n;
			h.valid_ = true;
			h.setup_status_ = sparse_ldl_status::success;
			h.symmetry_failed_on_ = 0;
		}

		// failure path: honest non-valid handle with the reason recorded
		template <typename _T, typename _Index>
		static void fail(ldl_shift_handle<_T, _Index>& h,
		                 const sparse_ldl_status why,
		                 const _Index n,
		                 const int symmetry_failed_on)
		{
			h.n_ = n;
			h.valid_ = false;
			h.setup_status_ = why;
			h.symmetry_failed_on_ = symmetry_failed_on;
		}
	};

	// -----------------------------------------------------------------------
	// dispatch_ldl_shift_setup_: signed-Index setup body.  Stage order
	// mirrors the one-shot entry (options -> conversion -> symmetry ->
	// merge -> structural singularity -> symbolic; tsparse_sparse_ldl.hpp
	// SS "processing order"), so every failure mode reuses the one-shot
	// status vocabulary.  B == 0 selects the A - sigma*I form.
	// -----------------------------------------------------------------------
	template <typename _T, typename _Index>
	inline typename std::enable_if<std::is_signed<_Index>::value,
	                               ldl_shift_handle<_T, _Index> >::type
	dispatch_ldl_shift_setup_(
	    const spmats<_T, _Index>& A,
	    const spmats<_T, _Index>* B,
	    const sparse_ldl_options<_T>& opt)
	{
		typedef typename vcp::tsparse_scalar::real_type<_T>::type R;
		ldl_shift_handle<_T, _Index> h;
		const _Index n = A.rowsize();

		// ---- options check (mirror of the one-shot entry step 2)
		if (opt.pivot_threshold == R(0)) { /* certified zero: accepted */ }
		else {
			shift_handle_access::fail(h, sparse_ldl_status::invalid_options, n, 0);
			return h;
		}
		sparse_ldl_method method_used;
		switch (opt.method) {
		case sparse_ldl_method::auto_select:
		case sparse_ldl_method::baseline_dynamic:
		case sparse_ldl_method::supernodal:
			method_used = sparse_ldl_detail::resolve_auto_method(opt.method);
			break;
		default:
			shift_handle_access::fail(h, sparse_ldl_status::invalid_options, n, 0);
			return h;
		}
		switch (opt.pivoting) {
		case sparse_ldl_pivoting::bk:
		case sparse_ldl_pivoting::none:
			break;
		default:
			shift_handle_access::fail(h, sparse_ldl_status::invalid_options, n, 0);
			return h;
		}
		switch (opt.diag_kernel) {
		case sparse_ldl_diag_kernel::auto_select:
		case sparse_ldl_diag_kernel::gemmtr:
		case sparse_ldl_diag_kernel::gemm:
			break;
		default:
			shift_handle_access::fail(h, sparse_ldl_status::invalid_options, n, 0);
			return h;
		}
		if (opt.ldl_min_block_size < 1) {
			shift_handle_access::fail(h, sparse_ldl_status::invalid_options, n, 0);
			return h;
		}
		switch (opt.ordering) {
		case sparse_ldl_ordering::auto_select:
		case sparse_ldl_ordering::natural:
		case sparse_ldl_ordering::rcm:
		case sparse_ldl_ordering::amd:
		case sparse_ldl_ordering::nested_dissection:
		case sparse_ldl_ordering::nested_dissection_ml:
			break;
		default:
			shift_handle_access::fail(h, sparse_ldl_status::invalid_options, n, 0);
			return h;
		}

		// ---- spmats -> CSC (full both-triangle CSC; the plan ignores the
		// strictly-upper entries through its -1 maps, so the lower merge is
		// canonical while the symmetry check below still sees the full input)
		csc_storage<_T, _Index> CA = sparse_lu_make_csc_storage(A);
		csc_storage<_T, _Index> CB;
		if (B) CB = sparse_lu_make_csc_storage(*B);

		// ---- certified symmetry checks (A always; B: symmetry ONLY, never
		// an SPD test -- H-4)
		if (opt.check_symmetry) {
			if (!sparse_ldl_detail::sparse_ldl_symmetry_certified_(
			        n, CA.col_ptr, CA.row_ind, CA.values, opt.symmetry_tol)) {
				shift_handle_access::fail(h, sparse_ldl_status::not_symmetric, n, 1);
				return h;
			}
			if (B &&
			    !sparse_ldl_detail::sparse_ldl_symmetry_certified_(
			        n, CB.col_ptr, CB.row_ind, CB.values, opt.symmetry_tol)) {
				shift_handle_access::fail(h, sparse_ldl_status::not_symmetric, n, 2);
				return h;
			}
		}

		// ---- pattern-merge plan (SH-0 layer)
		sparse_ldl_shift_plan<_T, _Index> plan =
		    B ? sparse_ldl_shift_plan<_T, _Index>(n, CA.col_ptr, CA.row_ind,
		                                          CB.col_ptr, CB.row_ind)
		      : sparse_ldl_shift_plan<_T, _Index>(n, CA.col_ptr, CA.row_ind);
		if (plan.status != sparse_ldl_shift_plan_status::success) {
			shift_handle_access::fail(h, sparse_ldl_status::invalid_input, n, 0);
			return h;
		}

		// ---- structural-singularity check on the merged lower pattern
		// (mirror of the one-shot entry step 4; sigma-independent, so it
		// belongs to setup)
		{
			const std::size_t un = static_cast<std::size_t>(n);
			std::vector<char> occupied(un, char(0));
			for (std::size_t c = 0; c < un; ++c) {
				for (_Index k = plan.col_ptr[c]; k < plan.col_ptr[c + 1u]; ++k) {
					occupied[static_cast<std::size_t>(
					    plan.row_ind[static_cast<std::size_t>(k)])] = char(1);
					occupied[c] = char(1);
				}
			}
			for (std::size_t i = 0; i < un; ++i) {
				if (occupied[i] == char(0)) {
					shift_handle_access::fail(
					    h, sparse_ldl_status::structural_singularity, n, 0);
					return h;
				}
			}
		}

		// ---- symbolic phase on the merged pattern (const-shared afterwards)
		sparse_ldl_symbolic_options sopt;
		sopt.ordering = opt.ordering;
		sopt.level = (method_used == sparse_ldl_method::supernodal)
		           ? sparse_ldl_symbolic_level::full
		           : sparse_ldl_symbolic_level::ordering_only;
		sparse_ldl_symbolic_result<_Index> sym =
		    sparse_ldl_symbolic_analyze(n, plan.col_ptr, plan.row_ind, sopt);
		if (sym.status != sparse_ldl_symbolic_status::success) {
			shift_handle_access::fail(
			    h,
			    (sym.status == sparse_ldl_symbolic_status::invalid_options)
			        ? sparse_ldl_status::invalid_options
			        : ((sym.status == sparse_ldl_symbolic_status::invalid_input)
			               ? sparse_ldl_status::invalid_input
			               : sparse_ldl_status::internal_error),
			    n, 0);
			return h;
		}
		const sparse_ldl_ordering ordering_used = sym.ordering_used;

		shift_handle_access::assign(
		    h, std::move(plan), std::move(sym),
		    std::move(CA.values),
		    B ? std::move(CB.values) : std::vector<_T>(),
		    opt, method_used, ordering_used, n);
		return h;
	}

	// unsigned-Index path: sparse LDL cannot be used (Index must be signed);
	// same reporting convention as the other sparse dispatches.
	template <typename _T, typename _Index>
	inline typename std::enable_if<!std::is_signed<_Index>::value,
	                               ldl_shift_handle<_T, _Index> >::type
	dispatch_ldl_shift_setup_(
	    const spmats<_T, _Index>& A,
	    const spmats<_T, _Index>* B,
	    const sparse_ldl_options<_T>& opt)
	{
		(void)A; (void)B; (void)opt;
		vcp::throw_error<vcp::state_error>(
		    "spmats::policy_ldl_shift_setup_with_info: sparse LDL requires a signed Index type");
		return ldl_shift_handle<_T, _Index>();
	}

} // namespace spmats_ldl_shift_detail

// ---------------------------------------------------------------------------
// factor_core_ (SH-2): shared sigma-dependent factor body -- assembly ->
// numeric under the given options -> L/D/perm materialization.  The
// materialization block is the verbatim mirror of the one-shot assembly in
// spmats_ldl_detail::dispatch_sparse_ldl_ (byte-identity contract of
// acceptance 1; machine-checked by the G-S2 tests), and the diagnostics
// conduit REUSES the one-shot's make_ldl_result_ (include-order reuse, no
// second copy).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
template <class SpMat>
ldl_result<_T, typename ldl_shift_handle<_T, _Index>::factor_index_type>
ldl_shift_handle<_T, _Index>::factor_core_(
	const _T& sigma, const sparse_ldl_options<_T>& opt2,
	const sparse_ldl_method method_used2, SpMat& L, SpMat& D,
	std::vector<factor_index_type>& p, workspace_type& ws) const
{
	typedef factor_index_type FI;
	ldl_result<_T, FI> out;                    // default status internal_error
	L.resize(FI(0), FI(0));
	D.resize(FI(0), FI(0));
	p.clear();
	if (!valid_) {
		// info-only (H-2): an unset / failed handle answers with its setup
		// status (internal_error for a default-constructed handle).
		out.status = setup_status_;
		return out;
	}
	try {
		if (!plan_.build_values(sigma, valA_, valB_, ws.val)) {
			out.status = sparse_ldl_status::internal_error;
			return out;
		}
		spmats_ldl_shift_detail::reset_result_scalars_(ws.res);
		ws.res.method_used   = method_used2;
		ws.res.ordering_used = ordering_used_;

		// symbolic: const-reuse of the handle's analysis; a retry that
		// RAISES the method to supernodal from an ordering_only handle
		// computes the full analysis locally (this call only -- the handle
		// itself never mutates, H-3 / concurrency contract).
		const sparse_ldl_symbolic_result<FI>* symp = &sym_;
		sparse_ldl_symbolic_result<FI> local_sym;
		if (method_used2 == sparse_ldl_method::supernodal &&
		    sym_.level != sparse_ldl_symbolic_level::full) {
			sparse_ldl_symbolic_options sopt;
			sopt.ordering = opt_.ordering;
			sopt.level = sparse_ldl_symbolic_level::full;
			local_sym = sparse_ldl_symbolic_analyze(
			    static_cast<FI>(n_), plan_.col_ptr, plan_.row_ind, sopt);
			if (local_sym.status != sparse_ldl_symbolic_status::success ||
			    local_sym.ordering_used != ordering_used_) {
				out.status = sparse_ldl_status::internal_error;
				return out;
			}
			symp = &local_sym;
		}

		sparse_ldl_factorize_numeric_with_info(
		    static_cast<FI>(n_), plan_.col_ptr, plan_.row_ind, ws.val,
		    *symp, opt2, ws.numeric, ws.res);

		out = spmats_ldl_detail::make_ldl_result_<_T, FI>(ws.res);
		if (ws.res.status != sparse_ldl_status::success &&
		    ws.res.status != sparse_ldl_status::zero_pivot) {
			return out;   // outputs stay empty (same validity rule as one-shot)
		}

		// ---- L: kernel CSC is column-sorted, unit diagonal explicit,
		// certified zeros already dropped -> assign_csc accepts it directly.
		const FI n = static_cast<FI>(n_);
		L.assign_csc(n, n, ws.res.L_col_ptr, ws.res.L_row_ind, ws.res.L_val);

		// ---- D: symmetric block diagonal from the internal triple;
		// strictly zero values are structural zeros and are NOT stored
		// (design v2 SS5.2; same x == T(0) rule as the kernels).
		{
			const sparse_ldl_result<_T, FI>& r = ws.res;
			std::vector<FI> d_ptr(static_cast<std::size_t>(n) + 1u, FI(0));
			std::vector<FI> d_ind;
			std::vector<_T> d_val;
			const std::size_t un = static_cast<std::size_t>(n);
			for (std::size_t k = 0; k < un; ++k) {
				if (k > 0 && r.D_block2[k - 1u] != char(0)) {
					const _T& e = r.D_sub[k - 1u];
					if (e == _T(0)) { /* certified zero: not stored */ }
					else {
						d_ind.push_back(static_cast<FI>(k - 1u));
						d_val.push_back(e);
					}
				}
				const _T& dk = r.D_diag[k];
				if (dk == _T(0)) { /* certified zero: structural zero */ }
				else {
					d_ind.push_back(static_cast<FI>(k));
					d_val.push_back(dk);
				}
				if (r.D_block2[k] != char(0) && k + 1u < un) {
					const _T& e = r.D_sub[k];
					if (e == _T(0)) { /* certified zero: not stored */ }
					else {
						d_ind.push_back(static_cast<FI>(k + 1u));
						d_val.push_back(e);
					}
				}
				d_ptr[k + 1u] = static_cast<FI>(d_ind.size());
			}
			D.assign_csc(n, n, d_ptr, d_ind, d_val);
		}

		p = ws.res.perm;
		return out;
	} catch (const vcp::error&) {
		throw;   // misuse / state errors keep their throwing contract
	} catch (const std::exception&) {
		// P3 final protection net (firing = defect)
		ldl_result<_T, FI> bad;
		bad.status = sparse_ldl_status::internal_error;
		L.resize(FI(0), FI(0));
		D.resize(FI(0), FI(0));
		p.clear();
		return bad;
	}
}

// factor_at: the handle's own (setup-frozen) options.
template <typename _T, typename _Index>
template <class SpMat>
ldl_result<_T, typename ldl_shift_handle<_T, _Index>::factor_index_type>
ldl_shift_handle<_T, _Index>::factor_at(
	const _T& sigma, SpMat& L, SpMat& D,
	std::vector<factor_index_type>& p, workspace_type& ws) const
{
	return factor_core_(sigma, opt_, method_used_, L, D, p, ws);
}

// factor_at_with: explicit method x pivoting override (D-6).  Only the two
// enums are overridable -- everything else stays setup-frozen, and no
// automatic switching happens anywhere on this path.
template <typename _T, typename _Index>
template <class SpMat>
ldl_result<_T, typename ldl_shift_handle<_T, _Index>::factor_index_type>
ldl_shift_handle<_T, _Index>::factor_at_with(
	const _T& sigma, const sparse_ldl_method method,
	const sparse_ldl_pivoting pivoting, SpMat& L, SpMat& D,
	std::vector<factor_index_type>& p, workspace_type& ws) const
{
	typedef factor_index_type FI;
	switch (method) {
	case sparse_ldl_method::auto_select:
	case sparse_ldl_method::baseline_dynamic:
	case sparse_ldl_method::supernodal:
		break;
	default: {
		ldl_result<_T, FI> out;
		out.status = sparse_ldl_status::invalid_options;
		L.resize(FI(0), FI(0)); D.resize(FI(0), FI(0)); p.clear();
		return out;
	}
	}
	switch (pivoting) {
	case sparse_ldl_pivoting::bk:
	case sparse_ldl_pivoting::none:
		break;
	default: {
		ldl_result<_T, FI> out;
		out.status = sparse_ldl_status::invalid_options;
		L.resize(FI(0), FI(0)); D.resize(FI(0), FI(0)); p.clear();
		return out;
	}
	}
	sparse_ldl_options<_T> opt2 = opt_;
	opt2.method   = method;
	opt2.pivoting = pivoting;
	return factor_core_(
	    sigma, opt2, sparse_ldl_detail::resolve_auto_method(method),
	    L, D, p, ws);
}

// ---------------------------------------------------------------------------
// policy_ldl_shift_setup_with_info: NVI outers (non-virtual).  Finalize
// guarantee + squareness / dimension entry checks (throwing, misuse), then
// delegation to the virtual _impl.  Must never be overridden -- override
// the _impl (derived policies inherit the default unchanged; design SS2.2).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ldl_shift_handle<_T, _Index> spmats<_T, _Index>::policy_ldl_shift_setup_with_info(
	const ldl_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (A.rowsize() != A.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_ldl_shift_setup_with_info: matrix must be square");
	return policy_ldl_shift_setup_with_info_impl(opt);
}

template <typename _T, typename _Index>
ldl_shift_handle<_T, _Index> spmats<_T, _Index>::policy_ldl_shift_setup_with_info(
	const spmats<_T, _Index>& B, const ldl_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	if (!A.is_finalized()) A.finalize();
	if (!B.is_finalized()) B.finalize();
	if (A.rowsize() != A.columnsize() || B.rowsize() != B.columnsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_ldl_shift_setup_with_info: matrices must be square");
	if (B.rowsize() != A.rowsize())
		vcp::throw_error<vcp::dimension_error>(
		    "spmats::policy_ldl_shift_setup_with_info: A and B dimensions differ");
	return policy_ldl_shift_setup_with_info_impl(B, opt);
}

// ---------------------------------------------------------------------------
// policy_ldl_shift_setup_with_info_impl: virtual algorithm bodies (default:
// signed-Index guard -> conversion / symmetry / merge / symbolic -> handle
// assembly).  Runtime failure is a non-valid handle; misuse (vcp::error)
// keeps its throwing contract; the final std::exception net returns a
// non-valid internal_error handle (P3, firing = defect).
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
ldl_shift_handle<_T, _Index> spmats<_T, _Index>::policy_ldl_shift_setup_with_info_impl(
	const ldl_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	try {
		return spmats_ldl_shift_detail::dispatch_ldl_shift_setup_<_T, _Index>(
		    A, static_cast<const spmats<_T, _Index>*>(0), opt);
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		ldl_shift_handle<_T, _Index> h;
		spmats_ldl_shift_detail::shift_handle_access::fail(
		    h, sparse_ldl_status::internal_error, A.rowsize(), 0);
		return h;
	}
}

template <typename _T, typename _Index>
ldl_shift_handle<_T, _Index> spmats<_T, _Index>::policy_ldl_shift_setup_with_info_impl(
	const spmats<_T, _Index>& B, const ldl_options<_T>& opt) const
{
	const spmats<_T, _Index>& A = *this;   // WFIX-2: subject is *this
	try {
		return spmats_ldl_shift_detail::dispatch_ldl_shift_setup_<_T, _Index>(
		    A, &B, opt);
	} catch (const vcp::error&) {
		throw;
	} catch (const std::exception&) {
		ldl_shift_handle<_T, _Index> h;
		spmats_ldl_shift_detail::shift_handle_access::fail(
		    h, sparse_ldl_status::internal_error, A.rowsize(), 0);
		return h;
	}
}

} // namespace vcp

#endif // VCP_SPMATS_LDL_SHIFT_HPP_IMPL
#endif // VCP_SPMATS_LDL_SHIFT_HPP_TYPES
