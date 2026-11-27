// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spumar_base/spumar_umfpack.hpp
// RAII wrappers for the UMFPACK dl (64-bit index / double) interface: two-stage
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
#include <cstdint>
#include <memory>
#include <vector>

#include <vcp/spmats.hpp>
#include <vcp/spumar_base/spumar_convert.hpp>

namespace vcp {
namespace spumar_detail {

	// umfpack_api<_Index>: the UMFPACK entry points for the index type
	// _Index.  Only _Index = std::int64_t (umfpack_dl_*) exists; any other
	// index type is a compile error.  Every status is returned as int (the
	// UMFPACK_* codes).
	template <typename _Index> struct umfpack_api;

	// The 64-bit index type declared by <suitesparse/umfpack.h> is read off
	// the signature of umfpack_dl_symbolic (SuiteSparse 5.x: SuiteSparse_long
	// = long, 7.x: int64_t).  Declaration only; never called.
	template <typename _R, typename _I, typename... _Rest>
	_I umfpack_dl_index_of_(_R (*)(_I, _Rest...));

	template <>
	struct umfpack_api<std::int64_t> {
		typedef decltype(umfpack_dl_index_of_(&umfpack_dl_symbolic)) ss_index;

		static const ss_index* in(const std::int64_t* p) {
			return suitesparse_index_cast<const ss_index>(p);
		}
		static ss_index* out(std::int64_t* p) {
			return suitesparse_index_cast<ss_index>(p);
		}

		static void defaults(double* control) { umfpack_dl_defaults(control); }
		static int symbolic(const std::int64_t n_row, const std::int64_t n_col,
			const std::int64_t* Ap, const std::int64_t* Ai, const double* Ax,
			void** symbolic, const double* control, double* info) {
			return static_cast<int>(umfpack_dl_symbolic(
				static_cast<ss_index>(n_row), static_cast<ss_index>(n_col),
				in(Ap), in(Ai), Ax, symbolic, control, info));
		}
		static int numeric(const std::int64_t* Ap, const std::int64_t* Ai, const double* Ax,
			void* symbolic, void** numeric, const double* control, double* info) {
			return static_cast<int>(umfpack_dl_numeric(
				in(Ap), in(Ai), Ax, symbolic, numeric, control, info));
		}
		static int solve(const std::int64_t* Ap, const std::int64_t* Ai, const double* Ax,
			double* x, const double* b, void* numeric, const double* control, double* info) {
			return static_cast<int>(umfpack_dl_solve(UMFPACK_A,
				in(Ap), in(Ai), Ax, x, b, numeric, control, info));
		}
		static void free_symbolic(void** symbolic) { umfpack_dl_free_symbolic(symbolic); }
		static void free_numeric(void** numeric) { umfpack_dl_free_numeric(numeric); }
		static int get_lunz(std::int64_t* lnz, std::int64_t* unz, std::int64_t* n_row,
			std::int64_t* n_col, std::int64_t* nz_udiag, void* numeric) {
			return static_cast<int>(umfpack_dl_get_lunz(out(lnz), out(unz), out(n_row),
				out(n_col), out(nz_udiag), numeric));
		}
		// L (CSR), U (CSC), P, Q; no D, no recip flag, no Rs.
		static int get_numeric(std::int64_t* Lp, std::int64_t* Lj, double* Lx,
			std::int64_t* Up, std::int64_t* Ui, double* Ux,
			std::int64_t* P, std::int64_t* Q, void* numeric) {
			return static_cast<int>(umfpack_dl_get_numeric(out(Lp), out(Lj), Lx,
				out(Up), out(Ui), Ux, out(P), out(Q),
				static_cast<double*>(0), static_cast<ss_index*>(0),
				static_cast<double*>(0), numeric));
		}
	};

	// Stage-1 RAII handle: umfpack_{di,dl}_symbolic result.
	template <typename _Index>
	class umfpack_symbolic_handle_basic {
	public:
		umfpack_symbolic_handle_basic() : ptr_(0) {}
		~umfpack_symbolic_handle_basic() { reset(); }
		umfpack_symbolic_handle_basic(const umfpack_symbolic_handle_basic&) = delete;
		umfpack_symbolic_handle_basic& operator=(const umfpack_symbolic_handle_basic&) = delete;

		void** slot() { reset(); return &ptr_; }
		void* get() const { return ptr_; }
		void reset() {
			if (ptr_ != 0) { umfpack_api<_Index>::free_symbolic(&ptr_); ptr_ = 0; }
		}

	private:
		void* ptr_;
	};

	// Stage-2 RAII handle: umfpack_{di,dl}_numeric result.
	template <typename _Index>
	class umfpack_numeric_handle_basic {
	public:
		umfpack_numeric_handle_basic() : ptr_(0) {}
		~umfpack_numeric_handle_basic() { reset(); }
		umfpack_numeric_handle_basic(const umfpack_numeric_handle_basic&) = delete;
		umfpack_numeric_handle_basic& operator=(const umfpack_numeric_handle_basic&) = delete;

		void** slot() { reset(); return &ptr_; }
		void* get() const { return ptr_; }
		void reset() {
			if (ptr_ != 0) { umfpack_api<_Index>::free_numeric(&ptr_); ptr_ = 0; }
		}

	private:
		void* ptr_;
	};

	// umfpack_lu_basic<_Index>: one factorization = one object.
	//   umfpack_lu_basic<_Index> lu;
	//   lu.control[UMFPACK_...] = ...;   // adjust knobs BEFORE factorize()
	//   if (lu.factorize(A) == UMFPACK_OK) lu.solve(b, x);  // repeatable
	// The CSC copy of A is kept alive inside (umfpack_{di,dl}_solve re-reads
	// Ap/Ai/Ax for iterative refinement).  info[] reflects the most recent
	// umfpack_{di,dl}_* call; the per-stage statuses are kept separately.
	template <typename _Index>
	class umfpack_lu_basic {
	public:
		double control[UMFPACK_CONTROL];
		double info[UMFPACK_INFO];

		umfpack_lu_basic()
			: symbolic_status_(UMFPACK_ERROR_invalid_Symbolic_object),
			  numeric_status_(UMFPACK_ERROR_invalid_Numeric_object) {
			umfpack_api<_Index>::defaults(control);
			for (std::size_t i = 0; i < UMFPACK_INFO; i++) info[i] = 0.0;
		}

		umfpack_lu_basic(const umfpack_lu_basic&) = delete;
		umfpack_lu_basic& operator=(const umfpack_lu_basic&) = delete;

		// symbolic + numeric.  Returns UMFPACK_OK on full success; on a
		// symbolic failure the numeric stage is not attempted and the
		// symbolic status is returned.  UMFPACK_WARNING_singular_matrix
		// from the numeric stage is returned as-is (caller maps it).
		int factorize(const vcp::spmats<double, _Index>& A) {
			view_.reset(new umfpack_csc_view_basic<_Index>(A));
			const _Index n_row = view_->n_rows();
			const _Index n_col = view_->n_cols();
			symbolic_status_ = umfpack_api<_Index>::symbolic(n_row, n_col,
				view_->Ap(), view_->Ai(), view_->Ax(),
				symbolic_.slot(), control, info);
			if (symbolic_status_ != UMFPACK_OK) return symbolic_status_;
			numeric_status_ = umfpack_api<_Index>::numeric(
				view_->Ap(), view_->Ai(), view_->Ax(),
				symbolic_.get(), numeric_.slot(), control, info);
			return numeric_status_;
		}

		// x = A^{-1} b (UMFPACK_A).  Requires a prior successful factorize().
		int solve(const std::vector<double>& b, std::vector<double>& x) {
			x.assign(b.size(), 0.0);
			return umfpack_api<_Index>::solve(
				view_->Ap(), view_->Ai(), view_->Ax(),
				x.data(), b.data(), numeric_.get(), control, info);
		}

		int symbolic_status() const { return symbolic_status_; }
		int numeric_status() const { return numeric_status_; }
		void* numeric_handle() const { return numeric_.get(); }

	private:
		std::unique_ptr<umfpack_csc_view_basic<_Index> > view_;
		umfpack_symbolic_handle_basic<_Index> symbolic_;
		umfpack_numeric_handle_basic<_Index> numeric_;
		int symbolic_status_;
		int numeric_status_;
	};

} // namespace spumar_detail
} // namespace vcp

#endif // VCP_SPUMAR_UMFPACK_HPP
