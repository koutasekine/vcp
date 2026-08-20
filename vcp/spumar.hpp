// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// ---------------------------------------------------------------------------
// License note (3-line summary; full text: sandbox/docs/usage/spumar_usage.md §7):
//  1. This header itself is BSD; it delegates to external backends whose
//     licenses are: UMFPACK+AMD = LGPL-2.1+ (NOT GPL), ARPACK-NG = BSD.
//  2. spumar を include もリンクもしないコードに LGPL の影響は一切ない
//     (code that neither includes nor links spumar is BSD-only, entirely
//     unaffected by the LGPL).
//  3. Dynamic linking against the distro shared libraries (apt/brew) imposes
//     no source-disclosure obligation on your code; check the LGPL text
//     yourself for static linking / modified redistribution.
// ---------------------------------------------------------------------------
//
// spumar.hpp — external-delegation sparse policy (UMFPACK + ARPACK),
// design: sandbox/docs/design/spumar_design_v1.1.md.
//
//   class spumar : public spmats<double>
//   usage: vcp::spmatrix<double, vcp::spumar>   (T = double ONLY, B-3)
//
// Overrides exactly the four virtual _impl replacement points of design
// §1.1a (NVI discipline; the non-virtual outers are never overridden):
//   1. policy_lss_with_info_impl                (UMFPACK, this file; SPU-0)
//   2. policy_lu_with_info_impl                 (UMFPACK get_numeric; SPU-1)
//   3. policy_eigs_with_info_impl               (ARPACK; SPU-2)
//   4. policy_generalized_eigs_with_info_impl   (ARPACK mode 3; SPU-3)
// Everything else (LDL, inertia, LU consumers, mul, strict shells) runs the
// base spmats<double> implementation unchanged — in particular LDL / inertia
// have NO external oracle in this configuration and use the (slower) own
// code even through spumar.
//
// Dependency isolation (G4 / B-2): this header includes the external
// <suitesparse/umfpack.h> (via spumar_base/spumar_umfpack.hpp) and requires
// linking -lumfpack.  No vcp core header includes this file; code that does
// not opt in builds and runs without UMFPACK/ARPACK installed.
//
// Thread safety: set_default_options() mutates a process-global default and
// last_delegate_info is a mutable diagnostic cache written from const
// methods — NEITHER is thread-safe.  Do not call spumar policy methods on
// the same object from multiple threads, and do not race
// set_default_options() against spumar construction.

#pragma once

#ifndef VCP_SPUMAR_HPP
#define VCP_SPUMAR_HPP

#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <vcp/spmats.hpp>
#include <vcp/spumar_base/spumar_convert.hpp>
#include <vcp/spumar_base/spumar_umfpack.hpp>
#include <vcp/spumar_base/spumar_arpack.hpp>

namespace vcp {
namespace spumar_detail {

	// ---- exact-operator evaluation helpers (B-4 / P-6) -------------------
	// These mirror the base dispatch-layer formulas (spmats_eigs.hpp helpers
	// 13/13b/14/16) so that the C-1 verdict below is decided under the SAME
	// conditions as the base solvers.  The C-1 acceptance CHECK itself is
	// the shared vcp::tsparse::residual_acceptance_check_scaled_ (B-29: the
	// scale definition is never re-implemented locally).

	inline double vec_norm2_(const std::vector<double>& v) {
		double s = 0.0;
		for (std::size_t i = 0; i < v.size(); i++) s += v[i] * v[i];
		return std::sqrt(s);
	}

	// exact ||A||_inf (max abs row sum) — base helper 13b mirror
	inline double mat_inf_norm_(const vcp::spmats<double, int>& A) {
		vcp::spmats<double, int> C = A.as_csr();
		const std::vector<int>& outer = C.outer_index();
		const std::vector<double>& val = C.values();
		double best = 0.0;
		for (int i = 0; i < C.rowsize(); i++) {
			double s = 0.0;
			for (int p = outer[static_cast<std::size_t>(i)];
			     p < outer[static_cast<std::size_t>(i) + 1]; p++) {
				s += std::fabs(val[static_cast<std::size_t>(p)]);
			}
			if (s > best) best = s;
		}
		return best;
	}

	// ||A||_F — base helper 13 mirror
	inline double mat_frob_norm_(const vcp::spmats<double, int>& A) {
		vcp::spmats<double, int> C = A.as_csr();
		const std::vector<double>& val = C.values();
		double s = 0.0;
		for (std::size_t i = 0; i < val.size(); i++) s += val[i] * val[i];
		return std::sqrt(s);
	}

	// real eigenpair: ||A v - lambda v||_2 — base helper 14 mirror
	inline double pair_residual_abs_(const vcp::spmats<double, int>& A,
	                                 const double lambda, const std::vector<double>& v) {
		std::vector<double> r = A.mul_vec(v);
		for (std::size_t i = 0; i < r.size(); i++) r[i] -= lambda * v[i];
		return vec_norm2_(r);
	}

	// base helper 16 mirror: abs / (||A||_F ||v|| + |lambda| ||v|| + eps)
	inline double pair_residual_rel_(const double abs_res, const double frobA,
	                                 const double lambda_abs, const double vnorm) {
		const double denom = frobA * vnorm + lambda_abs * vnorm +
			std::numeric_limits<double>::epsilon();
		return abs_res / denom;
	}

	// generalized pair: ||A v - lambda B v||_2 — base helper 24 mirror
	inline double gen_pair_residual_abs_(const vcp::spmats<double, int>& A,
	                                     const vcp::spmats<double, int>& B,
	                                     const double lambda, const std::vector<double>& v) {
		std::vector<double> r = A.mul_vec(v);
		const std::vector<double> bv = B.mul_vec(v);
		for (std::size_t i = 0; i < r.size(); i++) r[i] -= lambda * bv[i];
		return vec_norm2_(r);
	}

	// base helper 25 mirror: abs / (||A||_F ||v|| + |lambda| ||B||_F ||v|| + eps)
	inline double gen_pair_residual_rel_(const double abs_res, const double frobA,
	                                     const double frobB, const double lambda_abs,
	                                     const double vnorm) {
		const double denom = frobA * vnorm + lambda_abs * frobB * vnorm +
			std::numeric_limits<double>::epsilon();
		return abs_res / denom;
	}

	// complex conjugate pair (u, v; lambda = re +- i*im), eig_result
	// convention: || A [u v] - [u v] B ||_F with B = [[re, im], [-im, re]].
	inline double pair_residual_abs_complex_(const vcp::spmats<double, int>& A,
	                                         const double re, const double im,
	                                         const std::vector<double>& u,
	                                         const std::vector<double>& v) {
		const std::vector<double> Au = A.mul_vec(u);
		const std::vector<double> Av = A.mul_vec(v);
		double s = 0.0;
		for (std::size_t i = 0; i < u.size(); i++) {
			const double r1 = Au[i] - (u[i] * re + v[i] * (-im));
			const double r2 = Av[i] - (u[i] * im + v[i] * re);
			s += r1 * r1 + r2 * r2;
		}
		return std::sqrt(s);
	}

} // namespace spumar_detail
} // namespace vcp

namespace vcp {

	// -----------------------------------------------------------------------
	// spumar_options (design §1.4, S-2): thin wrapper over the ARPACK ncv
	// knob and the main UMFPACK Control knobs.  Sentinel values mean "leave
	// the backend default untouched".
	// -----------------------------------------------------------------------
	struct spumar_options {
		// ARPACK subspace dimension ncv; 0 = automatic (design §3 P-1).
		std::size_t ncv;
		// UMFPACK Control[UMFPACK_IRSTEP]; negative = UMFPACK default.
		// Applies to the solve delegation path only (the factor-extraction
		// path is SCALE_NONE and IR-free by the §C convention, S-3).
		long long umf_irstep;
		// UMFPACK Control[UMFPACK_PIVOT_TOLERANCE]; negative = default.
		double umf_pivot_tol;

		spumar_options() : ncv(0), umf_irstep(-1), umf_pivot_tol(-1.0) {}
	};

	// -----------------------------------------------------------------------
	// spumar_delegate_info (S-5): raw diagnostic cache of the LAST delegated
	// call.  Written from const policy methods through a mutable member;
	// NOT thread-safe.  Raw backend codes are kept in full so that honest
	// mapping (G5) never loses information.
	// -----------------------------------------------------------------------
	struct spumar_delegate_info {
		std::string backend;              // "umfpack" / "arpack" / "" (none)
		std::string note;                 // e.g. SLU-knob non-application record
		long long umfpack_symbolic_status;
		long long umfpack_numeric_status;
		long long umfpack_solve_status;
		std::vector<double> umfpack_info;     // full Info[UMFPACK_INFO] snapshot
		std::vector<double> umfpack_control;  // full Control[UMFPACK_CONTROL] used
		long long arpack_info;            // dsaupd/dnaupd info (SPU-2+)
		long long arpack_ierr;            // dseupd/dneupd ierr (SPU-2+)
		std::vector<int> arpack_iparam;   // full iparam snapshot (SPU-2+)

		spumar_delegate_info()
			: umfpack_symbolic_status(0), umfpack_numeric_status(0),
			  umfpack_solve_status(0), arpack_info(0), arpack_ierr(0) {}

		void clear() {
			backend.clear();
			note.clear();
			umfpack_symbolic_status = 0;
			umfpack_numeric_status = 0;
			umfpack_solve_status = 0;
			umfpack_info.clear();
			umfpack_control.clear();
			arpack_info = 0;
			arpack_ierr = 0;
			arpack_iparam.clear();
		}
	};

	// -----------------------------------------------------------------------
	// spumar
	// -----------------------------------------------------------------------
	class spumar : public spmats<double> {
	public:
		typedef spmats<double> base_type;

		// Per-object options (design §1.4).  When used through
		// spmatrix<double, vcp::spumar>, the process-wide static default is
		// copied into muar_opt at construction time (below).
		spumar_options muar_opt;

		// Raw diagnostics of the last delegated call (S-5; mutable cache,
		// written from const _impl overrides; not thread-safe).
		mutable spumar_delegate_info last_delegate_info;

		spumar() : base_type(), muar_opt(default_options_ref()) {}

		// SPUM-FIX1 (design §2.1): spmatrix arithmetic and other
		// result-returning APIs assign a base spmats<double> value through
		// `static_cast<_P&>(C) = policy_xxx(...)` [spmatrix.hpp L452 etc.,
		// 28 lines].  A derived policy needs this base-to-derived assignment
		// hook, defined as assignment to the base subobject only
		// (spimats.hpp L238 と同一規律).  muar_opt / last_delegate_info are
		// NOT touched (design R-1 / R-2, SPI-1 precedent); the implicit copy
		// assignment operator=(const spumar&) remains implicitly defined.
		spumar& operator=(const base_type& rhs) {
			base_type::operator=(rhs);
			return *this;
		}

		// Process-wide default options (design §1.4; NOT thread-safe).
		static void set_default_options(const spumar_options& o) {
			default_options_ref() = o;
		}
		static spumar_options get_default_options() {
			return default_options_ref();
		}

		// -------------------------------------------------------------------
		// T = double enforcement (B-3).  spmatrix<T, P> does not check
		// T == P::value_type statically, so the mismatch is caught at the
		// first element insertion: spmatrix<T, spumar>::add/set forward
		// their own T here, and any non-double floating-point T fires the
		// static_assert below with an explicit message.  Integral arguments
		// (e.g. A.add(i, j, 1) from double-typed code) convert to double
		// exactly and are allowed through.
		// -------------------------------------------------------------------
		template <typename S>
		void add(const index_type i, const index_type j, const S& a) {
			static_assert(std::is_same<S, double>::value || std::is_integral<S>::value,
				"vcp::spumar is double-only: spmatrix<T, vcp::spumar> requires T == double (use spmatrix<T> / spmats<T> for other scalar types)");
			base_type::add(i, j, static_cast<double>(a));
		}
		template <typename S>
		void set(const index_type i, const index_type j, const S& a) {
			static_assert(std::is_same<S, double>::value || std::is_integral<S>::value,
				"vcp::spumar is double-only: spmatrix<T, vcp::spumar> requires T == double (use spmatrix<T> / spmats<T> for other scalar types)");
			base_type::set(i, j, static_cast<double>(a));
		}

		// -------------------------------------------------------------------
		// Override #1 (design §1.1a): linear system solve via UMFPACK.
		// Called through the non-virtual outer policy_lss_with_info, which
		// owns the finalize guarantee; *this arrives finalized.
		//
		// Contract (design §2, G5):
		//   - the requested opt.method is NOT dispatched: every request is
		//     delegated to the UMFPACK direct solve.  result.method keeps
		//     the REQUESTED value (S-3 spirit: existing result types are
		//     not extended); the actual backend is recorded in
		//     last_delegate_info.
		//   - UMFPACK defaults are allowed on this path (scaling + iterative
		//     refinement).  converged is decided by the UMFPACK
		//     factorization / solve STATUS, exactly as the base
		//     spmats dispatch_sparse_lu_ decides it from
		//     fac.info().success (D-7 parity).  opt.tol is the stopping
		//     criterion of the ITERATIVE solvers and is not an acceptance
		//     criterion for a direct solve, so it is not used here.
		//     Honesty is preserved by re-evaluating the residual with the
		//     exact operator (*this) and REPORTING it (report-only): the
		//     caller judges solution quality from result.residual_norm.
		//   - status mapping: UMFPACK_OK -> re-evaluate residual;
		//     WARNING_singular_matrix -> converged = false (x = 0, solve not
		//     attempted); any other status -> converged = false (x = 0).
		//     Raw statuses / Info / Control are kept in last_delegate_info.
		// -------------------------------------------------------------------
		virtual linear_solve_result<double> policy_lss_with_info_impl(
			const std::vector<double>& b,
			const linear_solve_options<double>& opt) const override
		{
			// Input validation: mirror of the base policy_lss_with_info_impl
			// entry checks (same exception types).
			if (this->rowsize() != this->columnsize()) {
				vcp::throw_error<vcp::dimension_error>("spumar::policy_lss_with_info: matrix must be square");
			}
			if (b.size() != static_cast<std::size_t>(this->rowsize())) {
				vcp::throw_error<vcp::dimension_error>("spumar::policy_lss_with_info: rhs dimension mismatch");
			}
			if (opt.max_iter == 0) {
				vcp::throw_error<vcp::invalid_argument>("spumar::policy_lss_with_info: max_iter must be positive");
			}
			if (!(opt.tol > 0.0)) {
				vcp::throw_error<vcp::invalid_argument>("spumar::policy_lss_with_info: tol must be positive");
			}

			last_delegate_info.clear();
			last_delegate_info.backend = "umfpack";
			last_delegate_info.note =
				"backend=umfpack (direct LU solve); requested linear_solver_method delegated as-is; "
				"converged decided by factorization/solve status (D-7 parity with spmats); "
				"residual re-evaluated with the exact operator and reported (report-only)";

			linear_solve_result<double> result;
			result.method = opt.method; // requested value preserved (S-3 spirit)
			result.iterations = 0;

			const std::size_t n = static_cast<std::size_t>(this->rowsize());
			if (n == 0) {
				// Trivial empty system: nothing to delegate.
				result.x.clear();
				result.solution.clear();
				result.converged = true;
				return result;
			}

			spumar_detail::umfpack_lu lu;
			if (muar_opt.umf_irstep >= 0) {
				lu.control[UMFPACK_IRSTEP] = static_cast<double>(muar_opt.umf_irstep);
			}
			if (muar_opt.umf_pivot_tol >= 0.0) {
				lu.control[UMFPACK_PIVOT_TOLERANCE] = muar_opt.umf_pivot_tol;
			}

			const int fac_status = lu.factorize(*this);
			last_delegate_info.umfpack_symbolic_status = lu.symbolic_status();
			last_delegate_info.umfpack_numeric_status = lu.numeric_status();
			last_delegate_info.umfpack_control.assign(lu.control, lu.control + UMFPACK_CONTROL);
			last_delegate_info.umfpack_info.assign(lu.info, lu.info + UMFPACK_INFO);

			if (fac_status != UMFPACK_OK) {
				// WARNING_singular_matrix and every hard error map to a
				// non-converged result (design §2); solve is not attempted
				// on a singular factorization (its solution would contain
				// Inf/NaN by the UMFPACK contract).
				result.x.assign(n, 0.0);
				result.solution = result.x;
				result.converged = false;
				return result;
			}

			std::vector<double> x;
			const int solve_status = lu.solve(b, x);
			last_delegate_info.umfpack_solve_status = solve_status;
			last_delegate_info.umfpack_info.assign(lu.info, lu.info + UMFPACK_INFO);

			if (solve_status != UMFPACK_OK) {
				result.x.assign(n, 0.0);
				result.solution = result.x;
				result.converged = false;
				return result;
			}

			// Self re-evaluation (B-4 / G5): residual with the exact operator
			// *this, reported but NOT used as the verdict.
			result.x = x;
			result.solution = x;
			const double ir_taken = lu.info[UMFPACK_IR_TAKEN];
			result.iterations = ir_taken > 0.0 ? static_cast<std::size_t>(ir_taken) : 0;
			// SPUM-CONV: 直接法の合否は分解・求解のステータスで決まる
			// (基底 spmats の dispatch_sparse_lu_ と同じ D-7 規約)。
			// opt.tol は反復解法の停止条件であって直接法の受理基準ではない
			// ため、ここでは使わない。残差フィールドは
			// set_linear_residual_fields が report-only で埋める
			// (呼び出し側が result.residual_norm を見る責任を負う —
			//  spmats が既にその規約)。
			// is_finite は防御的な最終ゲートであり、現状の実装では発火
			// しない(SPUM-CONV 設計 §3.3 の実測を参照)。
			vcp::tsparse_solvers::set_linear_residual_fields(result, *this, b);
			result.converged = vcp::tsparse_scalar::is_finite(result.residual_norm);
			return result;
		}

		// -------------------------------------------------------------------
		// Override #2 (design §1.1a, SPU-1): LU factor extraction via
		// UMFPACK get_numeric, mapped to the §C (SSC) convention
		// P A Q = L U, p/q new->old, L unit lower (explicit unit diagonal),
		// U upper, exact zeros not stored, columns ascending.
		// Called through the non-virtual outer policy_lu_with_info (finalize
		// guarantee + squareness check); the matrix-form A.lu(L,U,P,Q)
		// reaches this override through the base default materialization
		// (delegation chain, design (b)).
		//
		// Contract (S-3 / design §2):
		//   - Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE is FORCED (the §C
		//     form has no scaling); UMFPACK's permutation convention
		//     (PAQ)[k][l] = A[P[k]][Q[l]] coincides with SSC new->old.
		//   - opt.slu.equilibration == true -> unsupported_options WITHOUT
		//     running the factorization (base-parity entry rejection).
		//   - SLU-specific knobs (method / ordering) have no UMFPACK
		//     counterpart: they are ignored but NOT silently — method_used /
		//     ordering_used echo the REQUEST and last_delegate_info records
		//     the non-application (existing result types unchanged).
		//   - status mapping: UMFPACK_OK -> success (after factor checks) /
		//     WARNING_singular_matrix (or a structurally deficient U
		//     diagonal) -> invalid_factorization / other -> internal_error;
		//     raw codes kept in last_delegate_info.
		//   - L arrives from UMFPACK in row (CSR) form and is transposed
		//     into CSC by a counting pass (ascending rows by construction);
		//     exact zeros are dropped on both factors (§C storage rule).
		// -------------------------------------------------------------------
		virtual lu_extract_result<double, int> policy_lu_with_info_impl(
			spmats<double, int>& L, spmats<double, int>& U,
			std::vector<int>& p, std::vector<int>& q,
			const lu_extract_options<double>& opt) const override
		{
			lu_extract_result<double, int> out;
			out.method_used = opt.slu.method;     // request echoed (S-3)
			out.ordering_used = opt.slu.ordering; // request record (S-3)

			// out parameters are valid only on success (base contract)
			L.resize(0, 0);
			U.resize(0, 0);
			p.clear();
			q.clear();

			last_delegate_info.clear();
			last_delegate_info.backend = "umfpack";
			last_delegate_info.note =
				"backend=umfpack; SLU knobs (method/ordering) have no UMFPACK counterpart and were NOT applied "
				"(method_used/ordering_used echo the request); SCALE_NONE forced (§C)";

			// P-5 parity: reject equilibration at the entry, factorization
			// not attempted.
			if (opt.slu.equilibration) {
				out.status = sparse_lu_extract_status::unsupported_options;
				return out;
			}

			const int n = this->rowsize();
			if (n == 0) {
				std::vector<int> ptr(1, 0);
				L.assign_csc(0, 0, ptr, std::vector<int>(), std::vector<double>());
				U.assign_csc(0, 0, ptr, std::vector<int>(), std::vector<double>());
				out.nnz_L = 0;
				out.nnz_U = 0;
				out.status = sparse_lu_extract_status::success;
				return out;
			}

			try {
				spumar_detail::umfpack_lu lu;
				lu.control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE; // §C forced
				if (muar_opt.umf_pivot_tol >= 0.0) {
					lu.control[UMFPACK_PIVOT_TOLERANCE] = muar_opt.umf_pivot_tol;
				}

				const int fac_status = lu.factorize(*this);
				last_delegate_info.umfpack_symbolic_status = lu.symbolic_status();
				last_delegate_info.umfpack_numeric_status = lu.numeric_status();
				last_delegate_info.umfpack_control.assign(lu.control, lu.control + UMFPACK_CONTROL);
				last_delegate_info.umfpack_info.assign(lu.info, lu.info + UMFPACK_INFO);

				if (fac_status == UMFPACK_WARNING_singular_matrix) {
					out.status = sparse_lu_extract_status::invalid_factorization;
					return out;
				}
				if (fac_status != UMFPACK_OK) {
					out.status = sparse_lu_extract_status::internal_error;
					return out;
				}

				int lnz = 0, unz = 0, n_row = 0, n_col = 0, nz_udiag = 0;
				int status = umfpack_di_get_lunz(&lnz, &unz, &n_row, &n_col,
					&nz_udiag, lu.numeric_handle());
				if (status != UMFPACK_OK || n_row != n || n_col != n) {
					last_delegate_info.umfpack_solve_status = status;
					out.status = sparse_lu_extract_status::internal_error;
					return out;
				}
				if (nz_udiag < n) {
					// structurally deficient U diagonal = singular factor
					out.status = sparse_lu_extract_status::invalid_factorization;
					return out;
				}

				std::vector<int> Lp(static_cast<std::size_t>(n) + 1),
					Lj(static_cast<std::size_t>(lnz)),
					Up(static_cast<std::size_t>(n) + 1),
					Ui(static_cast<std::size_t>(unz)),
					Pv(static_cast<std::size_t>(n)), Qv(static_cast<std::size_t>(n));
				std::vector<double> Lx(static_cast<std::size_t>(lnz)),
					Ux(static_cast<std::size_t>(unz));
				status = umfpack_di_get_numeric(Lp.data(), Lj.data(), Lx.data(),
					Up.data(), Ui.data(), Ux.data(), Pv.data(), Qv.data(),
					static_cast<double*>(0), static_cast<int*>(0),
					static_cast<double*>(0), lu.numeric_handle());
				if (status != UMFPACK_OK) {
					last_delegate_info.umfpack_solve_status = status;
					out.status = sparse_lu_extract_status::internal_error;
					return out;
				}

				// ---- L: CSR (row form) -> CSC by counting transpose;
				//      exact zeros dropped (the unit diagonal is 1.0 and
				//      survives by value). Ascending rows per column follow
				//      from the ascending row sweep.
				std::vector<int> Lc_ptr(static_cast<std::size_t>(n) + 1, 0);
				for (int i = 0; i < n; i++) {
					for (int k = Lp[static_cast<std::size_t>(i)];
					     k < Lp[static_cast<std::size_t>(i) + 1]; k++) {
						if (Lx[static_cast<std::size_t>(k)] == 0.0) continue;
						Lc_ptr[static_cast<std::size_t>(Lj[static_cast<std::size_t>(k)]) + 1]++;
					}
				}
				for (int j = 0; j < n; j++) {
					Lc_ptr[static_cast<std::size_t>(j) + 1] += Lc_ptr[static_cast<std::size_t>(j)];
				}
				const std::size_t l_stored = static_cast<std::size_t>(Lc_ptr[static_cast<std::size_t>(n)]);
				std::vector<int> Lc_row(l_stored);
				std::vector<double> Lc_val(l_stored);
				{
					std::vector<int> next(Lc_ptr.begin(), Lc_ptr.end() - 1);
					for (int i = 0; i < n; i++) {
						for (int k = Lp[static_cast<std::size_t>(i)];
						     k < Lp[static_cast<std::size_t>(i) + 1]; k++) {
							const double v = Lx[static_cast<std::size_t>(k)];
							if (v == 0.0) continue;
							const int j = Lj[static_cast<std::size_t>(k)];
							const std::size_t dst = static_cast<std::size_t>(next[static_cast<std::size_t>(j)]++);
							Lc_row[dst] = i;
							Lc_val[dst] = v;
						}
					}
				}

				// ---- U: already CSC (column form, ascending rows); drop
				//      exact zeros.
				std::vector<int> Uc_ptr(static_cast<std::size_t>(n) + 1, 0);
				std::vector<int> Uc_row;
				std::vector<double> Uc_val;
				Uc_row.reserve(static_cast<std::size_t>(unz));
				Uc_val.reserve(static_cast<std::size_t>(unz));
				for (int j = 0; j < n; j++) {
					for (int k = Up[static_cast<std::size_t>(j)];
					     k < Up[static_cast<std::size_t>(j) + 1]; k++) {
						const double v = Ux[static_cast<std::size_t>(k)];
						if (v == 0.0) continue;
						Uc_row.push_back(Ui[static_cast<std::size_t>(k)]);
						Uc_val.push_back(v);
					}
					Uc_ptr[static_cast<std::size_t>(j) + 1] = static_cast<int>(Uc_val.size());
				}

				// A numerically zero U diagonal entry would have been
				// dropped above: re-verify the explicit diagonal survived
				// (defense in depth on top of nz_udiag).
				for (int j = 0; j < n; j++) {
					const int b = Uc_ptr[static_cast<std::size_t>(j)];
					const int e = Uc_ptr[static_cast<std::size_t>(j) + 1];
					if (b >= e || Uc_row[static_cast<std::size_t>(e) - 1] != j) {
						out.status = sparse_lu_extract_status::invalid_factorization;
						return out;
					}
				}

				L.assign_csc(n, n, Lc_ptr, Lc_row, Lc_val);
				U.assign_csc(n, n, Uc_ptr, Uc_row, Uc_val);
				p.assign(Pv.begin(), Pv.end());
				q.assign(Qv.begin(), Qv.end());
				out.nnz_L = static_cast<int>(l_stored);
				out.nnz_U = static_cast<int>(Uc_val.size());
				out.status = sparse_lu_extract_status::success;
				return out;
			} catch (const vcp::error&) {
				throw; // misuse keeps its throwing contract (base parity)
			} catch (const std::exception&) {
				lu_extract_result<double, int> fresh;
				fresh.method_used = opt.slu.method;
				fresh.ordering_used = opt.slu.ordering;
				fresh.status = sparse_lu_extract_status::internal_error;
				return fresh;
			}
		}

		// -------------------------------------------------------------------
		// Override #3 (design §1.1a, SPU-2): standard eigenproblem via ARPACK.
		// Called through the non-virtual outer policy_eigs_with_info (*this
		// arrives finalized).  Routing (S-4 / design §3):
		//   P-2 which mapping:      target                 sym      nonsym
		//                           largest_magnitude      LM       LM
		//                           largest_algebraic      LA       LR
		//                           smallest_algebraic     SA       SR
		//                           smallest_magnitude     (SI only, else fallback)
		//                           target_magnitude/real  (SI only, else fallback)
		//   P-3: opt.use_shift == true -> mode 3, OP = (A - sigma I)^{-1}
		//        (UMFPACK, factorize once before the loop), which = LM.
		//   Fallback (explicit qualified base call
		//   spmats<double>::policy_eigs_with_info_impl): unmappable target /
		//   SM or target_* without shift / dense_fallback_explicit /
		//   k > n-2 / k == 0 / n == 0.
		//   P-1 ncv: default min(n, max(2k+1, 20)); with muar_opt.ncv != 0:
		//        min(n, max(ncv, k+2)).
		//   P-5: nonsymmetric -> dnaupd/dneupd; complex pairs follow the
		//        base convention (adjacent (re,re)/(+im,-im), u/v vectors;
		//        refused honestly when allow_complex_pairs == false, D3-2).
		//   P-6: converged is decided ONLY by re-evaluated exact residuals
		//        under the shared C-1 scaled acceptance
		//        (residual_acceptance_check_scaled_, anorm = exact ||A||_inf);
		//        C-2 remains backend(ARPACK)-dependent and is disclosed in
		//        result.message.
		//   info = 1 (max_iter exhausted): honest termination — converged =
		//        false, converged_count = nconv, converged portion returned.
		// -------------------------------------------------------------------
		virtual eig_result<double> policy_eigs_with_info_impl(
			std::size_t k,
			const eig_options<double>& opt) const override
		{
			const std::size_t n = static_cast<std::size_t>(this->rowsize());

			// ---- routing (fallback reasons per P-2) ----
			std::string fb;
			if (opt.method == eig_solver_method::dense_fallback_explicit) {
				fb = "dense_fallback_explicit requested";
			} else if (n == 0 || k == 0 || this->rowsize() != this->columnsize()) {
				fb = "degenerate/invalid sizes are the base's contract";
			} else if (k > n - 2) {
				fb = "k > n-2 (outside ARPACK nev bounds)";
			} else if (!opt.use_shift &&
			           (opt.target == eig_target::smallest_magnitude ||
			            opt.target == eig_target::target_magnitude ||
			            opt.target == eig_target::target_real)) {
				fb = "target unmappable without shift-invert (SM/target_* need use_shift)";
			}
			if (!fb.empty()) {
				last_delegate_info.clear();
				last_delegate_info.backend = "spmats-base-fallback";
				last_delegate_info.note = "eigs fallback to own implementation: " + fb;
				return spmats<double>::policy_eigs_with_info_impl(k, opt);
			}

			last_delegate_info.clear();
			last_delegate_info.backend = "arpack";

			// ---- structure decision ----
			bool symmetric;
			switch (opt.structure) {
			case matrix_structure_hint::symmetric:
			case matrix_structure_hint::hermitian:
				symmetric = true; break;
			case matrix_structure_hint::general:
				symmetric = false; break;
			case matrix_structure_hint::auto_detect:
			default:
				symmetric = this->is_symmetric(); break;
			}

			// ---- which / mode (P-2, P-3) ----
			const int mode = opt.use_shift ? 3 : 1;
			const double sigma = opt.use_shift ? opt.shift : 0.0;
			const char* which = "LM";
			if (!opt.use_shift) {
				switch (opt.target) {
				case eig_target::largest_magnitude:  which = "LM"; break;
				case eig_target::largest_algebraic:  which = symmetric ? "LA" : "LR"; break;
				case eig_target::smallest_algebraic: which = symmetric ? "SA" : "SR"; break;
				default:                             which = "LM"; break; // unreachable (routed)
				}
			}

			// ---- ncv (P-1) ----
			std::size_t ncv_s = (muar_opt.ncv != 0)
				? (muar_opt.ncv > k + 2 ? muar_opt.ncv : k + 2)
				: (2 * k + 1 > 20 ? 2 * k + 1 : 20);
			if (ncv_s > n) ncv_s = n;
			const int ncv = static_cast<int>(ncv_s);
			const int max_iter = opt.max_iter > (static_cast<std::size_t>(1) << 30)
				? (1 << 30) : static_cast<int>(opt.max_iter);

			eig_result<double> result;
			result.requested_count = k;
			result.method = opt.method;
			result.used_method = symmetric ? "spumar/arpack(dsaupd)" : "spumar/arpack(dnaupd)";
			result.used_shift_invert = (mode == 3);
			result.used_subspace_dim = ncv_s;

			// ---- OP construction ----
			std::size_t op_count = 0;
			spumar_detail::umfpack_lu si_lu;
			if (mode == 3) {
				// A - sigma I, factorized once before the loop (P-3)
				spmats<double, int> S(*this);
				for (int i = 0; i < static_cast<int>(n); i++) S.add(i, i, -sigma);
				S.finalize();
				const int st = si_lu.factorize(S);
				last_delegate_info.umfpack_symbolic_status = si_lu.symbolic_status();
				last_delegate_info.umfpack_numeric_status = si_lu.numeric_status();
				last_delegate_info.umfpack_control.assign(si_lu.control, si_lu.control + UMFPACK_CONTROL);
				last_delegate_info.umfpack_info.assign(si_lu.info, si_lu.info + UMFPACK_INFO);
				if (st != UMFPACK_OK) {
					result.converged = false;
					result.status = "factorization_failed";
					result.failure_reason =
						"shift-invert factorization failed (A - sigma*I; raw UMFPACK status in last_delegate_info)";
					result.message = result.failure_reason;
					return result;
				}
			}
			const spmats<double, int>* self = this;
			spumar_detail::arpack_op_fn opx;
			if (mode == 1) {
				opx = [self, &op_count](const double* x, double* y, const double*) {
					self->mul_vec(x, y);
					op_count++;
				};
			} else {
				const std::size_t nn = n;
				opx = [&si_lu, &op_count, nn](const double* x, double* y, const double*) {
					std::vector<double> b(x, x + nn), sol;
					si_lu.solve(b, sol);
					for (std::size_t i = 0; i < nn; i++) y[i] = sol[i];
					op_count++;
				};
			}
			spumar_detail::arpack_b_fn bx =
				[](const double*, double*) {}; // bmat='I': never called

			// ---- drive ----
			spumar_detail::arpack_outcome oc = symmetric
				? spumar_detail::arpack_drive_symmetric(static_cast<int>(n),
					static_cast<int>(k), ncv, which, 'I', mode, sigma,
					opt.tol, max_iter, opx, bx)
				: spumar_detail::arpack_drive_nonsymmetric(static_cast<int>(n),
					static_cast<int>(k), ncv, which, 'I', mode, sigma, 0.0,
					opt.tol, max_iter, opx, bx);

			last_delegate_info.arpack_info = oc.aupd_info;
			last_delegate_info.arpack_ierr = oc.eupd_ierr;
			last_delegate_info.arpack_iparam = oc.iparam;

			return assemble_arpack_result_(result, oc, *this, k, opt,
				symmetric, mode, op_count);
		}

		// -------------------------------------------------------------------
		// Override #4 (design §1.1a, SPU-3): generalized eigenproblem via
		// ARPACK, ALWAYS mode 3 (P-4: bmat='G', OP = (A - sigma B)^{-1} B,
		// default sigma = 0; mode 2 is never used).  B is the operand
		// argument (WFIX-2 (a) form).
		//
		// Delegable window (honesty first — mode 3 with which='LM' finds the
		// eigenvalues NEAREST sigma, so only sigma-expressible requests are
		// delegated):
		//   - opt.use_shift == true                  -> sigma = opt.shift
		//   - opt.target == smallest_magnitude       -> sigma = 0
		// Everything else (largest_*, smallest_algebraic, target_* without
		// shift, dense_fallback_explicit, k > n-2, degenerate sizes) falls
		// back to spmats<double>::policy_generalized_eigs_with_info_impl
		// (explicit qualified call, reason recorded in last_delegate_info).
		//
		// A - sigma B is factorized ONCE by UMFPACK before the loop; a
		// singular pencil shift maps to an honest factorization_failed
		// result (G-S3.3).  converged is decided ONLY by re-evaluated exact
		// generalized residuals ||A x - lambda B x|| under the shared
		// generalized C-1 scaled acceptance (B-48: the scale is
		// vcp::tsparse::residual_acceptance_check_generalized_scaled_, never
		// re-implemented).  A complex pair inside the returned window is
		// refused honestly on this path (generalized complex-pair return is
		// out of spumar's scope; documented).
		// -------------------------------------------------------------------
		virtual eig_result<double> policy_generalized_eigs_with_info_impl(
			const spmats<double, int>& B,
			std::size_t k, const eig_options<double>& opt) const override
		{
			const std::size_t n = static_cast<std::size_t>(this->rowsize());

			std::string fb;
			if (opt.method == eig_solver_method::dense_fallback_explicit) {
				fb = "dense_fallback_explicit requested";
			} else if (n == 0 || k == 0 || this->rowsize() != this->columnsize() ||
			           B.rowsize() != this->rowsize() || B.rowsize() != B.columnsize()) {
				fb = "degenerate/invalid sizes are the base's contract";
			} else if (k > n - 2) {
				fb = "k > n-2 (outside ARPACK nev bounds)";
			} else if (!opt.use_shift && opt.target != eig_target::smallest_magnitude) {
				fb = "target not sigma-expressible in mode 3 (delegable: use_shift or smallest_magnitude)";
			}
			if (!fb.empty()) {
				last_delegate_info.clear();
				last_delegate_info.backend = "spmats-base-fallback";
				last_delegate_info.note = "generalized eigs fallback to own implementation: " + fb;
				return spmats<double>::policy_generalized_eigs_with_info_impl(B, k, opt);
			}

			last_delegate_info.clear();
			last_delegate_info.backend = "arpack";
			last_delegate_info.note = "generalized delegation: ARPACK mode 3 (bmat='G'), OP=(A-sigma*B)^{-1}B";

			bool symmetric;
			switch (opt.structure) {
			case matrix_structure_hint::symmetric:
			case matrix_structure_hint::hermitian:
				symmetric = true; break;
			case matrix_structure_hint::general:
				symmetric = false; break;
			case matrix_structure_hint::auto_detect:
			default:
				symmetric = this->is_symmetric(); break;
			}
			symmetric = symmetric && B.is_symmetric();

			const double sigma = opt.use_shift ? opt.shift : 0.0;

			std::size_t ncv_s = (muar_opt.ncv != 0)
				? (muar_opt.ncv > k + 2 ? muar_opt.ncv : k + 2)
				: (2 * k + 1 > 20 ? 2 * k + 1 : 20);
			if (ncv_s > n) ncv_s = n;
			const int ncv = static_cast<int>(ncv_s);
			const int max_iter = opt.max_iter > (static_cast<std::size_t>(1) << 30)
				? (1 << 30) : static_cast<int>(opt.max_iter);

			eig_result<double> result;
			result.requested_count = k;
			result.method = opt.method;
			result.used_method = symmetric ? "spumar/arpack(dsaupd,gen)" : "spumar/arpack(dnaupd,gen)";
			result.used_shift_invert = true;
			result.used_generalized_operator = true;
			result.used_subspace_dim = ncv_s;

			// A - sigma B, factorized once (P-4).
			spumar_detail::umfpack_lu si_lu;
			{
				spmats<double, int> S(*this);
				if (sigma != 0.0) {
					spmats<double, int> Bc = B.as_csr();
					const std::vector<int>& outer = Bc.outer_index();
					const std::vector<int>& inner = Bc.inner_index();
					const std::vector<double>& val = Bc.values();
					for (int i = 0; i < Bc.rowsize(); i++) {
						for (int p = outer[static_cast<std::size_t>(i)];
						     p < outer[static_cast<std::size_t>(i) + 1]; p++) {
							S.add(i, inner[static_cast<std::size_t>(p)],
								-sigma * val[static_cast<std::size_t>(p)]);
						}
					}
				}
				S.finalize();
				const int st = si_lu.factorize(S);
				last_delegate_info.umfpack_symbolic_status = si_lu.symbolic_status();
				last_delegate_info.umfpack_numeric_status = si_lu.numeric_status();
				last_delegate_info.umfpack_control.assign(si_lu.control, si_lu.control + UMFPACK_CONTROL);
				last_delegate_info.umfpack_info.assign(si_lu.info, si_lu.info + UMFPACK_INFO);
				if (st != UMFPACK_OK) {
					result.converged = false;
					result.status = "factorization_failed";
					result.failure_reason =
						"generalized shift-invert factorization failed (A - sigma*B singular?; "
						"raw UMFPACK status in last_delegate_info)";
					result.message = result.failure_reason;
					return result;
				}
			}

			std::size_t solve_count = 0, bmul_count = 0;
			const spmats<double, int>* Bp = &B;
			const std::size_t nn = n;
			spumar_detail::arpack_op_fn opx =
				[&si_lu, Bp, &solve_count, &bmul_count, nn](const double* x, double* y, const double* bx_hint) {
					std::vector<double> rhs;
					if (bx_hint != 0) {
						rhs.assign(bx_hint, bx_hint + nn);
					} else {
						rhs.assign(nn, 0.0);
						Bp->mul_vec(x, rhs.data());
						bmul_count++;
					}
					std::vector<double> sol;
					si_lu.solve(rhs, sol);
					solve_count++;
					for (std::size_t i = 0; i < nn; i++) y[i] = sol[i];
				};
			spumar_detail::arpack_b_fn bx =
				[Bp, &bmul_count](const double* x, double* y) {
					Bp->mul_vec(x, y);
					bmul_count++;
				};

			spumar_detail::arpack_outcome oc = symmetric
				? spumar_detail::arpack_drive_symmetric(static_cast<int>(n),
					static_cast<int>(k), ncv, "LM", 'G', 3, sigma,
					opt.tol, max_iter, opx, bx)
				: spumar_detail::arpack_drive_nonsymmetric(static_cast<int>(n),
					static_cast<int>(k), ncv, "LM", 'G', 3, sigma, 0.0,
					opt.tol, max_iter, opx, bx);

			last_delegate_info.arpack_info = oc.aupd_info;
			last_delegate_info.arpack_ierr = oc.eupd_ierr;
			last_delegate_info.arpack_iparam = oc.iparam;

			result.iterations = static_cast<std::size_t>(oc.niter > 0 ? oc.niter : 0);
			result.linear_solves = solve_count;
			result.matrix_vector_products = bmul_count;
			result.converged_count = static_cast<std::size_t>(oc.nconv > 0 ? oc.nconv : 0);

			if ((oc.aupd_info != 0 && oc.aupd_info != 1) || oc.eupd_ierr != 0) {
				result.converged = false;
				result.status = "backend_error";
				result.failure_reason =
					"ARPACK error on the generalized path (raw codes in last_delegate_info)";
				result.message = result.failure_reason;
				return result;
			}

			const std::size_t m = oc.dr.size();
			bool has_complex = false;
			for (std::size_t i = 0; i < m; i++) {
				if (oc.di[i] != 0.0) { has_complex = true; break; }
			}
			if (has_complex) {
				result.converged = false;
				result.status = "complex_pair_refused";
				result.failure_reason =
					"complex conjugate pair inside the generalized returned window; "
					"refused honestly (generalized complex-pair return is outside spumar's scope)";
				result.message = result.failure_reason;
				return result;
			}

			const double frobA = spumar_detail::mat_frob_norm_(*this);
			const double frobB = spumar_detail::mat_frob_norm_(B);
			const double anorm = spumar_detail::mat_inf_norm_(*this);
			const double bnorm = spumar_detail::mat_inf_norm_(B);
			result.eigenvalues.assign(oc.dr.begin(), oc.dr.end());
			result.eigenvectors = oc.z;
			result.residuals_absolute.assign(m, 0.0);
			result.residuals_relative.assign(m, 0.0);
			for (std::size_t i = 0; i < m; i++) {
				const double lam = result.eigenvalues[i];
				const double abs_res = spumar_detail::gen_pair_residual_abs_(*this, B, lam, result.eigenvectors[i]);
				const double vn = spumar_detail::vec_norm2_(result.eigenvectors[i]);
				result.residuals_absolute[i] = abs_res;
				result.residuals_relative[i] =
					spumar_detail::gen_pair_residual_rel_(abs_res, frobA, frobB, std::fabs(lam), vn);
			}
			result.matrix_vector_products += 2 * m; // exact-operator re-evaluation cost
			vcp::tsparse_eigen_selection::sort_eigenpairs_by_target(
				result.eigenvalues, result.eigenvectors,
				result.residuals_absolute, result.residuals_relative,
				opt.target, opt.shift);
			std::vector<double> theta_abs(m, 0.0);
			for (std::size_t i = 0; i < m; i++) theta_abs[i] = std::fabs(result.eigenvalues[i]);
			result.complex_eigenvalues.clear();
			for (std::size_t i = 0; i < m; i++) {
				result.complex_eigenvalues.push_back(
					eig_result<double>::eigenvalue_type(result.eigenvalues[i]));
			}
			result.returned_count = m;
			result.returned_real_count = m;
			if (!result.residuals_absolute.empty()) {
				result.residual_norm_absolute = result.residuals_absolute[0];
				result.residual_norm_relative = result.residuals_relative[0];
				for (std::size_t i = 1; i < m; i++) {
					if (result.residuals_absolute[i] > result.residual_norm_absolute)
						result.residual_norm_absolute = result.residuals_absolute[i];
					if (result.residuals_relative[i] > result.residual_norm_relative)
						result.residual_norm_relative = result.residuals_relative[i];
				}
			}

			const bool c1_ok = vcp::tsparse::residual_acceptance_check_generalized_scaled_(
				result.residuals_absolute, result.residuals_relative, opt.tol,
				theta_abs, anorm, bnorm);
			const bool count_ok = (m >= k);

			if (oc.honest_nonconv()) {
				result.converged = false;
				result.status = "not_converged";
				result.failure_reason =
					"ARPACK reached max_iter (info=1); converged portion returned (honest termination)";
				result.message = result.failure_reason;
				return result;
			}
			if (c1_ok && count_ok) {
				result.converged = true;
				result.status = "converged";
				result.message = "converged (generalized C-1 re-evaluated with the exact operators; "
					"C-2 rests on ARPACK's internal deflation and is backend-dependent — "
					"disclosed per design §3 P-6)";
			} else {
				result.converged = false;
				result.status = count_ok ? "residual_check_failed" : "insufficient_count";
				result.failure_reason = count_ok
					? "end-of-run exact generalized residual failed the shared C-1 scaled acceptance"
					: "fewer pairs returned than requested";
				result.message = result.failure_reason;
			}
			return result;
		}

	private:
		// Shared result assembly for the standard eigs delegation (also the
		// error/honest-termination mapping).  Kept private and non-virtual.
		eig_result<double> assemble_arpack_result_(
			eig_result<double>& result,
			const spumar_detail::arpack_outcome& oc,
			const spmats<double, int>& A,
			const std::size_t k,
			const eig_options<double>& opt,
			const bool symmetric,
			const int mode,
			const std::size_t op_count) const
		{
			result.iterations = static_cast<std::size_t>(oc.niter > 0 ? oc.niter : 0);
			if (mode == 3) result.linear_solves = op_count;
			else result.matrix_vector_products = op_count;
			result.converged_count = static_cast<std::size_t>(oc.nconv > 0 ? oc.nconv : 0);

			// hard backend errors: no eigenpairs to report
			if ((oc.aupd_info != 0 && oc.aupd_info != 1) || oc.eupd_ierr != 0) {
				result.converged = false;
				result.status = "backend_error";
				result.failure_reason = std::string("ARPACK error (") +
					(oc.eupd_ierr != 0 ? "dseupd/dneupd ierr" : "dsaupd/dnaupd info") +
					" != 0; raw codes in last_delegate_info)";
				result.message = result.failure_reason;
				return result;
			}

			const std::size_t m = oc.dr.size();

			// ---- complex-window handling (P-5, base D3-2 convention) ----
			std::size_t pair_count = 0;
			bool has_complex = false;
			for (std::size_t i = 0; i < m; i++) {
				if (oc.di[i] != 0.0) { has_complex = true; break; }
			}
			if (has_complex && !opt.allow_complex_pairs) {
				// honest refusal: return the REAL entries as diagnostics only
				std::vector<double> vals;
				std::vector<std::vector<double> > vecs;
				for (std::size_t i = 0; i < m; i++) {
					if (oc.di[i] == 0.0) { vals.push_back(oc.dr[i]); vecs.push_back(oc.z[i]); }
				}
				fill_real_pairs_diagnostics_(result, A, vals, vecs, opt);
				result.converged = false;
				result.status = "complex_pair_refused";
				result.failure_reason =
					"complex conjugate pair inside the returned window; refused honestly "
					"(allow_complex_pairs == false; base D3-2 convention)";
				result.message = result.failure_reason;
				return result;
			}

			// ---- assemble eigenpairs ----
			const double frobA = spumar_detail::mat_frob_norm_(A);
			const double anorm = spumar_detail::mat_inf_norm_(A);
			result.eigenvalues.assign(oc.dr.begin(), oc.dr.end());
			result.eigenvectors = oc.z;
			result.residuals_absolute.assign(m, 0.0);
			result.residuals_relative.assign(m, 0.0);
			std::vector<double> theta_abs(m, 0.0);
			std::size_t extra_mv = 0;

			if (!has_complex) {
				for (std::size_t i = 0; i < m; i++) {
					const double lam = oc.dr[i];
					const double abs_res = spumar_detail::pair_residual_abs_(A, lam, oc.z[i]);
					const double vn = spumar_detail::vec_norm2_(oc.z[i]);
					result.residuals_absolute[i] = abs_res;
					result.residuals_relative[i] =
						spumar_detail::pair_residual_rel_(abs_res, frobA, std::fabs(lam), vn);
					theta_abs[i] = std::fabs(lam);
					extra_mv += 1;
				}
				// target-order re-fix (base parity); adjacency not a concern
				vcp::tsparse_eigen_selection::sort_eigenpairs_by_target(
					result.eigenvalues, result.eigenvectors,
					result.residuals_absolute, result.residuals_relative,
					opt.target, opt.shift);
				for (std::size_t i = 0; i < m; i++) theta_abs[i] = std::fabs(result.eigenvalues[i]);
				result.complex_eigenvalues.clear();
				for (std::size_t i = 0; i < m; i++) {
					result.complex_eigenvalues.push_back(
						eig_result<double>::eigenvalue_type(result.eigenvalues[i]));
				}
				result.returned_real_count = m;
				result.returned_complex_count = 0;
			} else {
				// opt-in complex pairs: keep dneupd adjacency (re,re)/(+im,-im)
				result.eigenvalues_imag.assign(oc.di.begin(), oc.di.end());
				result.complex_eigenvalues.clear();
				std::size_t complex_entries = 0;
				for (std::size_t i = 0; i < m; i++) {
					const double re = oc.dr[i];
					const double im = oc.di[i];
					result.complex_eigenvalues.push_back(
						eig_result<double>::eigenvalue_type(re, im));
					theta_abs[i] = std::sqrt(re * re + im * im);
					if (im > 0.0 && i + 1 < m && oc.di[i + 1] == -im) {
						// pair (i, i+1): u = z[i], v = z[i+1]
						const double abs_res = spumar_detail::pair_residual_abs_complex_(
							A, re, im, oc.z[i], oc.z[i + 1]);
						double uv = 0.0;
						{
							const double nu = spumar_detail::vec_norm2_(oc.z[i]);
							const double nv = spumar_detail::vec_norm2_(oc.z[i + 1]);
							uv = std::sqrt(nu * nu + nv * nv); // ||[u v]||_F
						}
						result.residuals_absolute[i] = abs_res;
						result.residuals_absolute[i + 1] = abs_res;
						result.residuals_relative[i] =
							spumar_detail::pair_residual_rel_(abs_res, frobA, theta_abs[i], uv);
						result.residuals_relative[i + 1] = result.residuals_relative[i];
						pair_count++;
						complex_entries += 2;
						extra_mv += 2;
					} else if (im == 0.0) {
						const double abs_res = spumar_detail::pair_residual_abs_(A, re, oc.z[i]);
						const double vn = spumar_detail::vec_norm2_(oc.z[i]);
						result.residuals_absolute[i] = abs_res;
						result.residuals_relative[i] =
							spumar_detail::pair_residual_rel_(abs_res, frobA, std::fabs(re), vn);
						extra_mv += 1;
					}
				}
				result.complex_pair_count = pair_count;
				result.returned_complex_count = complex_entries;
				result.returned_real_count = m - complex_entries;
				result.message =
					"complex pairs kept in dneupd adjacency order (no target re-sort); ";
			}
			result.matrix_vector_products += extra_mv; // exact-operator re-evaluation cost
			result.returned_count = m;
			if (!result.residuals_absolute.empty()) {
				result.residual_norm_absolute = result.residuals_absolute[0];
				result.residual_norm_relative = result.residuals_relative[0];
				for (std::size_t i = 1; i < m; i++) {
					if (result.residuals_absolute[i] > result.residual_norm_absolute)
						result.residual_norm_absolute = result.residuals_absolute[i];
					if (result.residuals_relative[i] > result.residual_norm_relative)
						result.residual_norm_relative = result.residuals_relative[i];
				}
			}

			// ---- verdict (P-6: own exact residuals under the shared C-1
			//      scaled acceptance; backend self-report never decides) ----
			const bool c1_ok = vcp::tsparse::residual_acceptance_check_scaled_(
				result.residuals_absolute, result.residuals_relative, opt.tol,
				theta_abs, anorm);
			const bool count_ok = (m >= k);

			if (oc.honest_nonconv()) {
				result.converged = false;
				result.status = "not_converged";
				result.failure_reason =
					"ARPACK reached max_iter (info=1); converged portion returned "
					"(honest termination)";
				result.message += result.failure_reason;
				return result;
			}
			if (c1_ok && count_ok) {
				result.converged = true;
				result.status = "converged";
				result.message += "converged (C-1 re-evaluated with the exact operator; "
					"C-2 rests on ARPACK's internal deflation and is backend-dependent — "
					"disclosed per design §3 P-6)";
			} else {
				result.converged = false;
				result.status = count_ok ? "residual_check_failed" : "insufficient_count";
				result.failure_reason = count_ok
					? "end-of-run exact residual failed the shared C-1 scaled acceptance"
					: "fewer pairs returned than requested";
				result.message += result.failure_reason;
			}
			return result;
		}

		// residual diagnostics for the real subset (used by the honest
		// complex-window refusal path)
		void fill_real_pairs_diagnostics_(
			eig_result<double>& result,
			const spmats<double, int>& A,
			const std::vector<double>& vals,
			const std::vector<std::vector<double> >& vecs,
			const eig_options<double>& opt) const
		{
			(void)opt;
			const double frobA = spumar_detail::mat_frob_norm_(A);
			result.eigenvalues = vals;
			result.eigenvectors = vecs;
			const std::size_t m = vals.size();
			result.residuals_absolute.assign(m, 0.0);
			result.residuals_relative.assign(m, 0.0);
			for (std::size_t i = 0; i < m; i++) {
				const double abs_res = spumar_detail::pair_residual_abs_(A, vals[i], vecs[i]);
				const double vn = spumar_detail::vec_norm2_(vecs[i]);
				result.residuals_absolute[i] = abs_res;
				result.residuals_relative[i] =
					spumar_detail::pair_residual_rel_(abs_res, frobA, std::fabs(vals[i]), vn);
			}
			result.returned_count = m;
			result.returned_real_count = m;
			result.matrix_vector_products += m;
		}

		static spumar_options& default_options_ref() {
			static spumar_options o;
			return o;
		}
	};

} // namespace vcp

#endif // VCP_SPUMAR_HPP
