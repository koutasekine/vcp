// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// ---------------------------------------------------------------------------
// License note (3-line summary; full text: sandbox/docs/usage/spcmodumar_usage.md §7):
//  1. This header itself is BSD; it delegates the sparse Cholesky
//     factorization to CHOLMOD, whose modules carry DIFFERENT licenses:
//     Check / Cholesky / Core(Utility) / Partition = LGPL,
//     Supernodal / MatrixOps / Modify = GPL.  A combined work linking the
//     default full CHOLMOD build is subject to the GPL conditions; an NGPL
//     build (GPL modules excluded) is LGPL but simplicial-only.  Which to
//     choose is the USER's decision; this library does not interfere (no
//     #error / #warning enforcement -- notification only, design G2).
//  2. spcmodumar を include もリンクもしないコード(spumar 利用を含む)に
//     GPL/LGPL の影響は一切ない(code that neither includes nor links
//     spcmodumar -- including plain spumar users -- is entirely unaffected;
//     include isolation, machine-checked by acceptance criterion (f)).
//  3. Dynamic linking against the distro shared library (apt/brew) vs static
//     linking / redistribution have different obligations under GPL/LGPL;
//     check the license texts yourself.  This note is informational, not
//     legal advice.
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
	// spcmodumar
	// -----------------------------------------------------------------------
	class spcmodumar : public spumar {
	public:
		typedef spumar base_type;

		spcmodumar() : base_type() {}

		// -------------------------------------------------------------------
		// Override (design D-1): sparse LL^T Cholesky via CHOLMOD.  Called
		// through the non-virtual outer policy_chol_with_info, which owns the
		// finalize guarantee and the squareness entry check; *this arrives
		// finalized and square.  Signature transcribed from the base
		// declaration (spmats.hpp; T = double, Index = int).
		//
		// The delegation pipeline (design SS2, fixed order) is implemented in
		// SPCM-1; runtime failure is a status, never an exception (P3).
		// -------------------------------------------------------------------
		virtual chol_result<double, int> policy_chol_with_info_impl(
			spmats<double, int>& L,
			std::vector<int>& perm,
			const chol_options<double>& opt) const override;
	};

	// ---------------------------------------------------------------------
	// SPCM-1 body: the fixed delegation pipeline of design SS2 lives in
	// spcmodumar_detail::cholmod_chol_delegate_ (spcmodumar_base).  This
	// wrapper adds the exception discipline of the base _impl
	// (spmats_chol_impl.hpp, mirrored verbatim): misuse (vcp::error) keeps
	// its throwing contract; any other escaped exception is the P3 final
	// net and maps to internal_error (runtime failure is a status, never a
	// throw; SS2-8).
	// ---------------------------------------------------------------------
	inline chol_result<double, int> spcmodumar::policy_chol_with_info_impl(
		spmats<double, int>& L,
		std::vector<int>& perm,
		const chol_options<double>& opt) const
	{
		const spmats<double, int>& A = *this;   // subject is *this (WFIX-2 parity)
		try {
			return spcmodumar_detail::cholmod_chol_delegate_(A, L, perm, opt);
		} catch (const vcp::error&) {
			// misuse / state errors keep their throwing contract (base parity)
			throw;
		} catch (const std::exception&) {
			// P3 final protection net
			L.resize(0, 0);
			perm.clear();
			chol_result<double, int> out;
			out.status = sparse_chol_status::internal_error;
			return out;
		}
	}

} // namespace vcp

#endif // VCP_SPCMODUMAR_HPP
