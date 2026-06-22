// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_BLAS_HPP
#define VCP_TSPARSE_BLAS_HPP

#include <vcp/error.hpp>

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
}

#endif
