// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spumar_base/spumar_umfpack.hpp
// RAII wrappers for the UMFPACK di (int32 / double) interface: two-stage
// symbolic / numeric handles plus a bundled factorization object designed
// for "factorize once, solve many" (the solve delegation and the future
// shift-invert reverse-communication loop both reuse the numeric handle).
//
// This header is part of the spumar delegation layer: it includes the
// external <suitesparse/umfpack.h>.  No vcp core header may include it
// (dependency isolation, spumar design v1.1 G4 / B-2).

#pragma once

#ifndef VCP_SPUMAR_UMFPACK_HPP
#define VCP_SPUMAR_UMFPACK_HPP

#include <suitesparse/umfpack.h>

#include <cstddef>
#include <memory>
#include <vector>

#include <vcp/spmats.hpp>
#include <vcp/spumar_base/spumar_convert.hpp>

namespace vcp {
namespace spumar_detail {

	// Stage-1 RAII handle: umfpack_di_symbolic result.
	class umfpack_symbolic_handle {
	public:
		umfpack_symbolic_handle() : ptr_(0) {}
		~umfpack_symbolic_handle() { reset(); }
		umfpack_symbolic_handle(const umfpack_symbolic_handle&) = delete;
		umfpack_symbolic_handle& operator=(const umfpack_symbolic_handle&) = delete;

		void** slot() { reset(); return &ptr_; }
		void* get() const { return ptr_; }
		void reset() {
			if (ptr_ != 0) { umfpack_di_free_symbolic(&ptr_); ptr_ = 0; }
		}

	private:
		void* ptr_;
	};

	// Stage-2 RAII handle: umfpack_di_numeric result.
	class umfpack_numeric_handle {
	public:
		umfpack_numeric_handle() : ptr_(0) {}
		~umfpack_numeric_handle() { reset(); }
		umfpack_numeric_handle(const umfpack_numeric_handle&) = delete;
		umfpack_numeric_handle& operator=(const umfpack_numeric_handle&) = delete;

		void** slot() { reset(); return &ptr_; }
		void* get() const { return ptr_; }
		void reset() {
			if (ptr_ != 0) { umfpack_di_free_numeric(&ptr_); ptr_ = 0; }
		}

	private:
		void* ptr_;
	};

	// umfpack_lu: one factorization = one object.
	//   umfpack_lu lu;
	//   lu.control[UMFPACK_...] = ...;   // adjust knobs BEFORE factorize()
	//   if (lu.factorize(A) == UMFPACK_OK) lu.solve(b, x);  // repeatable
	// The CSC copy of A is kept alive inside (umfpack_di_solve re-reads
	// Ap/Ai/Ax for iterative refinement).  info[] reflects the most recent
	// umfpack_di_* call; the per-stage statuses are kept separately.
	class umfpack_lu {
	public:
		double control[UMFPACK_CONTROL];
		double info[UMFPACK_INFO];

		umfpack_lu()
			: symbolic_status_(UMFPACK_ERROR_invalid_Symbolic_object),
			  numeric_status_(UMFPACK_ERROR_invalid_Numeric_object) {
			umfpack_di_defaults(control);
			for (std::size_t i = 0; i < UMFPACK_INFO; i++) info[i] = 0.0;
		}

		umfpack_lu(const umfpack_lu&) = delete;
		umfpack_lu& operator=(const umfpack_lu&) = delete;

		// symbolic + numeric.  Returns UMFPACK_OK on full success; on a
		// symbolic failure the numeric stage is not attempted and the
		// symbolic status is returned.  UMFPACK_WARNING_singular_matrix
		// from the numeric stage is returned as-is (caller maps it).
		int factorize(const vcp::spmats<double, int>& A) {
			view_.reset(new umfpack_csc_view(A));
			const int n_row = view_->n_rows();
			const int n_col = view_->n_cols();
			symbolic_status_ = umfpack_di_symbolic(n_row, n_col,
				view_->Ap(), view_->Ai(), view_->Ax(),
				symbolic_.slot(), control, info);
			if (symbolic_status_ != UMFPACK_OK) return symbolic_status_;
			numeric_status_ = umfpack_di_numeric(
				view_->Ap(), view_->Ai(), view_->Ax(),
				symbolic_.get(), numeric_.slot(), control, info);
			return numeric_status_;
		}

		// x = A^{-1} b (UMFPACK_A).  Requires a prior successful factorize().
		int solve(const std::vector<double>& b, std::vector<double>& x) {
			x.assign(b.size(), 0.0);
			return umfpack_di_solve(UMFPACK_A,
				view_->Ap(), view_->Ai(), view_->Ax(),
				x.data(), b.data(), numeric_.get(), control, info);
		}

		int symbolic_status() const { return symbolic_status_; }
		int numeric_status() const { return numeric_status_; }
		void* numeric_handle() const { return numeric_.get(); }

	private:
		std::unique_ptr<umfpack_csc_view> view_;
		umfpack_symbolic_handle symbolic_;
		umfpack_numeric_handle numeric_;
		int symbolic_status_;
		int numeric_status_;
	};

} // namespace spumar_detail
} // namespace vcp

#endif // VCP_SPUMAR_UMFPACK_HPP
