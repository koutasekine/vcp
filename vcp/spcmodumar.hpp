// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// ---------------------------------------------------------------------------
// License note (summary; full text: docs/spmatrix.md, "外部ライブラリのライセンスに関する注意"):
//  1. This header itself is BSD.  It delegates the sparse Cholesky
//     factorization to CHOLMOD, whose modules carry different licenses
//     (Check / Cholesky / Core(Utility) / Partition = LGPL-2.1-or-later,
//     Supernodal / MatrixOps / Modify = GPL-2.0-or-later), and, through
//     spumar, uses UMFPACK (GPL-2.0-or-later), AMD (BSD-3-Clause) and
//     ARPACK-NG (BSD).
//  2. A combined work that uses spcmodumar links UMFPACK and is therefore
//     subject to the GPL, even with a CHOLMOD build without the GPL modules.
//  3. Code that neither includes nor links spcmodumar or spumar is not
//     affected by the GPL or the LGPL.  Check the license texts before
//     redistribution.  This note is informational, not legal advice.
// ---------------------------------------------------------------------------
//
// spcmodumar.hpp — external-delegation sparse policy adding a CHOLMOD-backed
// sparse Cholesky on top of spumar (UMFPACK + ARPACK),
// design: sandbox/docs/design/spcmodumar_design_v1.md v1.0.
//
//   class spcmodumar : public spumar
//   usage: vcp::spmatrix<double, vcp::spcmodumar>   (T = double ONLY, D-14;
//          enforcement inherited from spumar's add/set static_assert)
//
// Overrides exactly ONE virtual _impl replacement point (design D-1; NVI
// discipline -- the non-virtual outers are never overridden):
//   policy_chol_with_info_impl   (CHOLMOD analyze/factorize; SPCM-1)
// Everything else -- LU solve / LU factors / eigs / generalized eigs
// (UMFPACK + ARPACK via spumar) and LDL / inertia / LU consumers / mul
// (base spmats<double>) -- is INHERITED UNCHANGED (design G1).
//
// Contract transparency (design G3): the spmatrix/spmats chol contract
// (4 members, chol_options / chol_result, P^T A P = L L^T with perm
// new->old) is unchanged; all 4 orderings + auto_select are accepted
// (rcm included, via CHOLMOD_GIVEN with the VCP rcm permutation, D-4).
// Differences allowed between policies: the permutation DETAILS (CHOLMOD
// may compose a postorder -- the returned perm is always the L->Perm
// read-back), rounding of L, performance, and method_used (which may
// report `supernodal`, the D-9 enum value reserved for delegation).
//
// Dependency isolation (design SS7-3): this header includes the external
// <suitesparse/cholmod.h> (via spcmodumar_base/spcmodumar_cholmod.hpp) and
// requires linking -lcholmod.  No vcp core header includes this file; code
// that does not opt in builds and runs without CHOLMOD installed.

#pragma once

#ifndef VCP_SPCMODUMAR_HPP
#define VCP_SPCMODUMAR_HPP

#include <vector>

#include <vcp/spumar.hpp>
#include <vcp/spcmodumar_base/spcmodumar_cholmod.hpp>

namespace vcp {

	// -----------------------------------------------------------------------
	// spcmodumar_basic<_Index>: only _Index = std::int64_t is supported
	// (cholmod_l_*, CHOLMOD_LONG; spcmodumar_detail::cholmod_api<_Index>);
	// any other index type stops at the static_assert below.
	// vcp::spcmodumar is an alias (below).
	// -----------------------------------------------------------------------
	template <typename _Index>
	class spcmodumar_basic : public spumar_basic<_Index> {
		static_assert(std::is_same<_Index, std::int64_t>::value,
			"vcp::spcmodumar supports only the index type std::int64_t (it calls "
			"the 64-bit CHOLMOD interface cholmod_l_*); use vcp::spcmodumar");
	public:
		typedef spumar_basic<_Index> base_type;

		spcmodumar_basic() : base_type() {}

		// SPUM-FIX1 (design §2.2 案 A): the implicitly generated copy
		// assignment operator=(const spcmodumar&) hides the base's
		// operator=(const spmats<double>&) by name, so fixing spumar alone
		// does not fix this class.  Re-expose every base assignment operator;
		// the implicit copy assignment still wins for spcmodumar arguments.
		using base_type::operator=;

		// -------------------------------------------------------------------
		// Override (design D-1): sparse LL^T Cholesky via CHOLMOD.  Called
		// through the non-virtual outer policy_chol_with_info, which owns the
		// finalize guarantee and the squareness entry check; *this arrives
		// finalized and square.  Signature transcribed from the base
		// declaration (spmats.hpp; T = double, Index = _Index).
		//
		// The delegation pipeline (design SS2, fixed order) is implemented in
		// SPCM-1; runtime failure is a status, never an exception (P3).
		// -------------------------------------------------------------------
		virtual chol_result<double, _Index> policy_chol_with_info_impl(
			spmats<double, _Index>& L,
			std::vector<_Index>& perm,
			const chol_options<double>& opt) const override;
	};

	// vcp::spcmodumar: the std::int64_t index version (cholmod_l_*).
	typedef spcmodumar_basic<std::int64_t> spcmodumar;

	// ---------------------------------------------------------------------
	// SPCM-1 body: the fixed delegation pipeline of design SS2 lives in
	// spcmodumar_detail::cholmod_chol_delegate_ (spcmodumar_base).  This
	// wrapper adds the exception discipline of the base _impl
	// (spmats_chol_impl.hpp, mirrored verbatim): misuse (vcp::error) keeps
	// its throwing contract; any other escaped exception is the P3 final
	// net and maps to internal_error (runtime failure is a status, never a
	// throw; SS2-8).
	// ---------------------------------------------------------------------
	template <typename _Index>
	inline chol_result<double, _Index> spcmodumar_basic<_Index>::policy_chol_with_info_impl(
		spmats<double, _Index>& L,
		std::vector<_Index>& perm,
		const chol_options<double>& opt) const
	{
		const spmats<double, _Index>& A = *this;   // subject is *this (WFIX-2 parity)
		try {
			return spcmodumar_detail::cholmod_chol_delegate_(A, L, perm, opt);
		} catch (const vcp::error&) {
			// misuse / state errors keep their throwing contract (base parity)
			throw;
		} catch (const std::exception&) {
			// P3 final protection net
			L.resize(0, 0);
			perm.clear();
			chol_result<double, _Index> out;
			out.status = sparse_chol_status::internal_error;
			return out;
		}
	}

} // namespace vcp

#endif // VCP_SPCMODUMAR_HPP
