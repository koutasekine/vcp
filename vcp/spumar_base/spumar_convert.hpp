// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spumar_base/spumar_convert.hpp
// spmats (0-origin CSC) <-> UMFPACK (0-origin CSC, 64-bit indices)
// conversion.
// Both sides share the same storage convention (0-origin CSC, column pointers
// + ascending unique row indices after finalize), so "conversion" is a
// finalized CSC copy whose vectors are handed to UMFPACK as raw pointers.
//
// This header is part of the spumar delegation layer: it may include vcp
// core headers, but no vcp core header may include it (dependency isolation,
// spumar design v1.1 G4 / B-2).

#pragma once

#ifndef VCP_SPUMAR_CONVERT_HPP
#define VCP_SPUMAR_CONVERT_HPP

#include <cstdint>
#include <type_traits>
#include <vector>

#include <vcp/spmats.hpp>

namespace vcp {
namespace spumar_detail {

	// suitesparse_index_cast: the ONLY place where an index array of vcp
	// (std::int64_t) is handed to the 64-bit interface of SuiteSparse
	// (umfpack_dl_* / cholmod_l_*) without a copy.  The 64-bit index type
	// declared by the SuiteSparse headers is `long` (SuiteSparse 5.x:
	// SuiteSparse_long) or `int64_t` (SuiteSparse 7.x), which may be a type
	// distinct from std::int64_t (macOS: long vs long long).  Both are
	// checked to be signed 8-byte integers.
	template <typename _To, typename _From>
	inline _To* suitesparse_index_cast(_From* p) {
		static_assert(std::is_integral<_To>::value && std::is_signed<_To>::value
			&& sizeof(_To) == 8,
			"suitesparse_index_cast: the SuiteSparse index type must be a signed 8-byte integer");
		static_assert(std::is_integral<_From>::value && std::is_signed<_From>::value
			&& sizeof(_From) == 8,
			"suitesparse_index_cast: the vcp index type must be a signed 8-byte integer");
		static_assert(std::is_const<_To>::value || !std::is_const<_From>::value,
			"suitesparse_index_cast: must not cast away const");
		return reinterpret_cast<_To*>(p);
	}

	// umfpack_csc_view_basic<_Index>: owns a finalized CSC copy of A and exposes the raw
	// arrays in the exact form umfpack_dl_* (_Index = std::int64_t, through
	// suitesparse_index_cast) expects.  No raw pointers are
	// stored (they are recomputed from the owned vectors on each call), so
	// the object is safely copyable/movable.  The view must outlive every
	// UMFPACK call that receives its pointers (umfpack_dl_solve re-reads
	// Ap/Ai/Ax for iterative refinement).
	template <typename _Index>
	class umfpack_csc_view_basic {
	public:
		explicit umfpack_csc_view_basic(const vcp::spmats<double, _Index>& A)
			: csc_(A.as_csc()) {}

		_Index n_rows() const { return csc_.rowsize(); }
		_Index n_cols() const { return csc_.columnsize(); }
		const _Index* Ap() const { return csc_.outer_index().data(); }
		const _Index* Ai() const { return csc_.inner_index().data(); }
		const double* Ax() const { return csc_.values().data(); }

	private:
		vcp::spmats<double, _Index> csc_;
	};

} // namespace spumar_detail
} // namespace vcp

#endif // VCP_SPUMAR_CONVERT_HPP
