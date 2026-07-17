// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spmats_lu_factor.hpp
// LSS-1 P-4: reusable sparse-LU factorization handle (lu_factor_handle).
// TYPE header: included by spmats.hpp BEFORE the spmats<_T,_Index> class body
// (same layering as spmats_lu_extract.hpp).  The factorizing policy methods
// (policy_lu_factorize_with_info / _impl) are DECLARED in spmats.hpp and
// DEFINED at the end of spmats_base/spmats_lss.hpp, next to the sparse-LU
// solve dispatch whose result semantics the handle mirrors.

#ifndef VCP_SPMATS_LU_FACTOR_HPP
#define VCP_SPMATS_LU_FACTOR_HPP

#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_solvers.hpp>
#include <vcp/tsparse/tsparse_sparse_lu.hpp>
#include <vcp/spmats_base/spmats_eigs_types.hpp>

namespace vcp {

template <typename _T, typename _Index> class spmats;

namespace spmats_lu_factor_detail {
	// Sole writer of lu_factor_handle private state (befriended below).
	// Used only by the policy implementations in spmats_base/spmats_lss.hpp.
	struct handle_access;
} // namespace spmats_lu_factor_detail

// ---------------------------------------------------------------------------
// lu_factor_handle<_T,_Index> (design LSS-1 §3.1)
//
// Owns one sparse-LU factorization plus everything needed to solve
// repeatedly: the factorization-time sparse_lu_options (IR on/off, tol,
// maxit are FROZEN at factorization; there is no solve-time option API in
// v1) and a finalized CSR copy of A (required by sparse_lu_solve_refined
// for the IR residual; nnz-order memory, smaller than the factors).
//
// solve_with_info result semantics are IDENTICAL to the policy_lss
// sparse_lu dispatch (spmats_lss_detail::dispatch_sparse_lu_, post-P-6):
//   IR on : residual fields from sparse_lu_refinement_info,
//           iterations = IR iterations.
//   IR off: honest computed residual, initial = absolute (0-iteration
//           reading), iterations = 0.
//   converged = factorization success (D-7).
// Contracts:
//   valid()==false + solve_with_info : NO throw; converged=false, x = 0
//                                      (D5 contract).  solve() throws.
//   rhs size != n                    : vcp::dimension_error (policy_lss
//                                      contract), valid or not.
// Byte-identity contract (§3.5, T-10): for the same A and the same
// sparse_lu_options, solve_with_info(b).x is byte-identical to
// policy_lss_with_info(b, {method=sparse_lu, same sparse_lu}).x -- both
// paths run the same factorize + solve function chain.
//
// unsigned _Index: the handle type must remain INSTANTIABLE (it is the
// return type of policy_lu_factorize_with_info, whose unsigned path throws
// vcp::state_error at RUNTIME -- same SFINAE discipline as
// dispatch_sparse_lu_).  The sparse-LU storage types static_assert on
// signed Index, so the internal factorization is stored with the
// signed-mapped factor_index_type, and the signed-only solve body is
// SFINAE-split (never instantiated for unsigned _Index).  For signed
// _Index -- the only usable case -- factor_index_type == _Index, so all
// public signatures match the design verbatim.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
class lu_factor_handle {
public:
	typedef typename vcp::tsparse_scalar::real_type<_T>::type real_type;
	typedef typename std::conditional<std::is_signed<_Index>::value, _Index,
		typename std::make_signed<_Index>::type>::type factor_index_type;

	lu_factor_handle() : fac_(), A_(), opt_(), n_(_Index(0)), valid_(false) {}

	// factorization success (usable for solves)
	bool valid() const { return valid_; }

	// factorization status (sparse_lu_status::success iff valid())
	sparse_lu_status status() const { return fac_.info().status; }

	const sparse_lu_info<_T, factor_index_type>& info() const { return fac_.info(); }

	linear_solve_result<_T> solve_with_info(const std::vector<_T>& b) const {
		if (b.size() != static_cast<std::size_t>(n_)) {
			vcp::throw_error<vcp::dimension_error>(
				"spmats::lu_factor_handle::solve_with_info: rhs dimension mismatch");
		}
		linear_solve_result<_T> result;
		result.method = linear_solver_method::sparse_lu;
		result.iterations = 0;
		if (!valid_) {
			// D5 contract: no throw; converged==false marks the residual
			// fields as undefined (do not read).
			result.converged = false;
			result.x.assign(static_cast<std::size_t>(n_), _T(0));
			result.solution = result.x;
			return result;
		}
		run_solve_(b, result);
		result.solution  = result.x;
		result.converged = true;
		return result;
	}

	// throwing wrapper (same contract as policy_lss)
	std::vector<_T> solve(const std::vector<_T>& b) const {
		linear_solve_result<_T> result = solve_with_info(b);
		if (!result.converged) {
			vcp::throw_error<vcp::state_error>(
				"spmats::lu_factor_handle::solve: factorization is not valid");
		}
		return result.x;
	}

private:
	// signed Index path: the actual solve (mirrors dispatch_sparse_lu_).
	template <typename _I = _Index>
	typename std::enable_if<std::is_signed<_I>::value, void>::type
	run_solve_(const std::vector<_T>& b, linear_solve_result<_T>& result) const {
		if (opt_.iterative_refinement) {
			// IR path: identical field mapping to dispatch_sparse_lu_.
			vcp::sparse_lu_refinement_info<_T> ir_info;
			result.x = vcp::sparse_lu_solve_refined(A_, fac_, b, opt_, &ir_info);
			result.iterations             = ir_info.iterations;
			result.initial_residual_norm  = ir_info.initial_residual;
			result.residual_norm          = ir_info.final_residual;
			result.absolute_residual_norm = ir_info.final_residual;
			result.relative_residual_norm = ir_info.final_relative_residual;
		} else {
			// Plain path: honest residual (P-6 semantics shared).
			result.x = fac_.solve(b);
			vcp::tsparse_solvers::set_linear_residual_fields(result, A_, b);
			result.initial_residual_norm = result.absolute_residual_norm;
		}
	}

	// unsigned Index path: unreachable at runtime (valid_ is always false --
	// the factorize dispatch throws before a handle is produced), split so
	// the signed-only code above is never instantiated for unsigned _Index.
	template <typename _I = _Index>
	typename std::enable_if<!std::is_signed<_I>::value, void>::type
	run_solve_(const std::vector<_T>& b, linear_solve_result<_T>& result) const {
		(void)b; (void)result;
		vcp::throw_error<vcp::state_error>(
			"spmats::lu_factor_handle::solve_with_info: sparse_lu requires a signed Index type");
	}

	sparse_lu_factorization<_T, factor_index_type> fac_;
	spmats<_T, _Index> A_;          // finalized copy for the IR residual (§3.3)
	sparse_lu_options<_T> opt_;     // options frozen at factorization time
	_Index n_;
	bool valid_;

	friend struct spmats_lu_factor_detail::handle_access;
};

namespace spmats_lu_factor_detail {

	struct handle_access {
		template <typename _T, typename _Index>
		static void assign(lu_factor_handle<_T, _Index>& handle,
		                   sparse_lu_factorization<_T,
		                       typename lu_factor_handle<_T, _Index>::factor_index_type>&& fac,
		                   spmats<_T, _Index>&& A_finalized,
		                   const sparse_lu_options<_T>& opt,
		                   const _Index n,
		                   const bool valid)
		{
			handle.fac_   = std::move(fac);
			handle.A_     = std::move(A_finalized);
			handle.opt_   = opt;
			handle.n_     = n;
			handle.valid_ = valid;
		}
	};

} // namespace spmats_lu_factor_detail

} // namespace vcp

#endif // VCP_SPMATS_LU_FACTOR_HPP
