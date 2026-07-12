// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spumar_base/spumar_convert.hpp
// spmats (0-origin CSC) <-> UMFPACK (0-origin CSC, int32 indices) conversion.
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

#include <vector>

#include <vcp/spmats.hpp>

namespace vcp {
namespace spumar_detail {

	// umfpack_csc_view: owns a finalized CSC copy of A and exposes the raw
	// arrays in the exact form umfpack_di_* expects.  No raw pointers are
	// stored (they are recomputed from the owned vectors on each call), so
	// the object is safely copyable/movable.  The view must outlive every
	// UMFPACK call that receives its pointers (umfpack_di_solve re-reads
	// Ap/Ai/Ax for iterative refinement).
	class umfpack_csc_view {
	public:
		explicit umfpack_csc_view(const vcp::spmats<double, int>& A)
			: csc_(A.as_csc()) {}

		int n_rows() const { return csc_.rowsize(); }
		int n_cols() const { return csc_.columnsize(); }
		const int* Ap() const { return csc_.outer_index().data(); }
		const int* Ai() const { return csc_.inner_index().data(); }
		const double* Ax() const { return csc_.values().data(); }

	private:
		vcp::spmats<double, int> csc_;
	};

} // namespace spumar_detail
} // namespace vcp

#endif // VCP_SPUMAR_CONVERT_HPP
