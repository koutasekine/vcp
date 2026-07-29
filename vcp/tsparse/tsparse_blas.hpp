// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_BLAS_HPP
#define VCP_TSPARSE_BLAS_HPP

#include <atomic>
#include <cstddef>
#include <exception>

#include <vcp/error.hpp>

// SPOMP-1 (design §2, D-1): two-stage OpenMP guard (replicated per file so
// each header stands alone; same block as spmats.hpp / tsparse_spgemm.hpp).
#ifdef VCP_NOMP
#  ifndef VCP_SPARSE_NOMP
#    define VCP_SPARSE_NOMP
#  endif
#endif
#if defined(_OPENMP) && !defined(VCP_SPARSE_NOMP)
#  define VCP_SPARSE_USE_OPENMP 1
#  include <omp.h>
#else
#  define VCP_SPARSE_USE_OPENMP 0
#endif

namespace vcp {

	template <typename _Index, typename _T>
	void tcsrmv(const char trans,
	            _Index rows, _Index cols,
	            const _T& alpha,
	            const _Index* row_ptr,
	            const _Index* col_ind,
	            const _T* val,
	            const _T* x,
	            const _T& beta,
	            _T* y) {
		if (rows < 0 || cols < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcsrmv: negative size");
		}
		if (trans == 'N' || trans == 'n') {
			for (_Index i = 0; i < rows; i++) y[i] *= beta;
			for (_Index i = 0; i < rows; i++) {
				_T sum = _T(0);
				for (_Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
					sum += val[p] * x[col_ind[p]];
				}
				y[i] += alpha * sum;
			}
		}
		else if (trans == 'T' || trans == 't') {
			for (_Index j = 0; j < cols; j++) y[j] *= beta;
			for (_Index i = 0; i < rows; i++) {
				const _T xi = x[i];
				for (_Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
					y[col_ind[p]] += alpha * val[p] * xi;
				}
			}
		}
		else {
			vcp::throw_error<vcp::invalid_argument>("tcsrmv: trans must be 'N' or 'T'");
		}
	}

	template <typename _Index, typename _T>
	void tcscmv(const char trans,
	            _Index rows, _Index cols,
	            const _T& alpha,
	            const _Index* col_ptr,
	            const _Index* row_ind,
	            const _T* val,
	            const _T* x,
	            const _T& beta,
	            _T* y) {
		if (rows < 0 || cols < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcscmv: negative size");
		}
		if (trans == 'N' || trans == 'n') {
			for (_Index i = 0; i < rows; i++) y[i] *= beta;
			for (_Index j = 0; j < cols; j++) {
				const _T xj = x[j];
				for (_Index p = col_ptr[j]; p < col_ptr[j + 1]; p++) {
					y[row_ind[p]] += alpha * val[p] * xj;
				}
			}
		}
		else if (trans == 'T' || trans == 't') {
			for (_Index j = 0; j < cols; j++) y[j] *= beta;
			for (_Index j = 0; j < cols; j++) {
				_T sum = _T(0);
				for (_Index p = col_ptr[j]; p < col_ptr[j + 1]; p++) {
					sum += val[p] * x[row_ind[p]];
				}
				y[j] += alpha * sum;
			}
		}
		else {
			vcp::throw_error<vcp::invalid_argument>("tcscmv: trans must be 'N' or 'T'");
		}
	}

	// -----------------------------------------------------------------------
	// SPOMP-1 (design §4.3, D-4/D-9): gather-side parallel SpMV kernels.
	// PURE ADDITIONS -- the sequential tcsrmv / tcscmv above are untouched
	// and keep serving every scatter-side call ('T' on CSR / 'N' on CSC) and
	// every non-spmats caller (e.g. tsparse_generalized_shift_invert.hpp).
	// tcsrmv_gather_par mirrors the 'N' branch of tcsrmv; tcscmv_gather_par
	// mirrors the 'T' branch of tcscmv.  Each y[i] (resp. y[j]) is produced
	// by the identical per-row (per-column) operation sequence of the
	// sequential branch, so the result is bit-identical for any thread
	// count.  Exceptions from the loop bodies (kv arithmetic) are captured
	// once (atomic election, design §6-1) and rethrown outside the parallel
	// region.  Loop variables are signed (design §6-2).  The work-threshold
	// decision is made by the spmats caller (design §4.3), not here.
	// When VCP_SPARSE_USE_OPENMP == 0 these compile to plain sequential
	// loops with the same operation sequence.
	// -----------------------------------------------------------------------
	template <typename _Index, typename _T>
	void tcsrmv_gather_par(_Index rows, _Index cols,
	                       const _T& alpha,
	                       const _Index* row_ptr,
	                       const _Index* col_ind,
	                       const _T* val,
	                       const _T* x,
	                       const _T& beta,
	                       _T* y) {
		if (rows < 0 || cols < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcsrmv_gather_par: negative size");
		}
		std::atomic<bool> caught(false);
		std::exception_ptr eptr;
		const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(rows);
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (std::ptrdiff_t i = 0; i < n; i++) {
			try {
				y[i] *= beta;
			}
			catch (...) {
				bool expected = false;
				if (caught.compare_exchange_strong(expected, true)) {
					eptr = std::current_exception();
				}
			}
		}
		if (caught.load()) std::rethrow_exception(eptr);
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (std::ptrdiff_t i = 0; i < n; i++) {
			try {
				_T sum = _T(0);
				for (_Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
					sum += val[p] * x[col_ind[p]];
				}
				y[i] += alpha * sum;
			}
			catch (...) {
				bool expected = false;
				if (caught.compare_exchange_strong(expected, true)) {
					eptr = std::current_exception();
				}
			}
		}
		if (caught.load()) std::rethrow_exception(eptr);
	}

	template <typename _Index, typename _T>
	void tcscmv_gather_par(_Index rows, _Index cols,
	                       const _T& alpha,
	                       const _Index* col_ptr,
	                       const _Index* row_ind,
	                       const _T* val,
	                       const _T* x,
	                       const _T& beta,
	                       _T* y) {
		if (rows < 0 || cols < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcscmv_gather_par: negative size");
		}
		std::atomic<bool> caught(false);
		std::exception_ptr eptr;
		const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(cols);
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (std::ptrdiff_t j = 0; j < n; j++) {
			try {
				y[j] *= beta;
			}
			catch (...) {
				bool expected = false;
				if (caught.compare_exchange_strong(expected, true)) {
					eptr = std::current_exception();
				}
			}
		}
		if (caught.load()) std::rethrow_exception(eptr);
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (std::ptrdiff_t j = 0; j < n; j++) {
			try {
				_T sum = _T(0);
				for (_Index p = col_ptr[j]; p < col_ptr[j + 1]; p++) {
					sum += val[p] * x[row_ind[p]];
				}
				y[j] += alpha * sum;
			}
			catch (...) {
				bool expected = false;
				if (caught.compare_exchange_strong(expected, true)) {
					eptr = std::current_exception();
				}
			}
		}
		if (caught.load()) std::rethrow_exception(eptr);
	}
}

#endif
