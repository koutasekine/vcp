// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and / or other materials provided with the distribution.
// * Neither the name of the Kouta Sekine nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL KOUTA SEKINE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//
// ----------------------------------------------------------------------------
// This file contains code derived from the reference BLAS contained in the
// reference LAPACK 3.12.1 distribution (https://www.netlib.org/lapack/),
// Copyright (c) 1992-2023 The University of Tennessee and The University
//                         of Tennessee Research Foundation,
// Copyright (c) 2000-2023 The University of California Berkeley,
// Copyright (c) 2006-2023 The University of Colorado Denver.
// All rights reserved.
// Distributed under the modified BSD license; see vlapack/LICENSE_LAPACK.txt
// for the full license text, the list of conditions and the disclaimer.
// Modified by Kouta Sekine (2026): translated to C++, rounding-mode control
// added, O(N^3) kernels replaced by rmatmul-based blocked implementations.
// ----------------------------------------------------------------------------

#pragma once

#ifndef VBLAS_RDBLAS_LEVEL2_HPP
#define VBLAS_RDBLAS_LEVEL2_HPP

// Level 2 BLAS (丸めモード指定付き)．行列は column-major，index は 0-based．
// 引数は reference BLAS と同順 + 末尾に rounding_mode (1 上向き / -1 下向き / 0 最近点)．
// gemv / ger / symv / syr / syr2 は OpenMP 並列 (各 worker thread が丸めモードを
// 保存・変更・復元)，banded / packed / 三角系は逐次 (O(N^2) 以下で律速にならない)．

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <omp.h>

#include "rdblas_common.hpp"

#pragma STDC FENV_ACCESS ON

namespace vblas_rdblas_detail {

// y(leny, incy 付き) := beta*y (beta==0 は exact 0 埋め)
inline void scale_vector(const int leny, const double beta, double* y, const int incy, const int fe_mode) {
	if (leny <= 0 || beta == 1.0) {
		return;
	}
	int iy = vec_start(leny, incy);
	if (beta == 0.0) {
		for (int i = 0; i < leny; i++) {
			y[iy] = 0.0;
			iy += incy;
		}
		return;
	}
	RoundingGuard guard(fe_mode);
	for (int i = 0; i < leny; i++) {
		y[iy] = beta * y[iy];
		iy += incy;
	}
}

// packed 格納の列先頭 offset (upper: 列 j は AP[j(j+1)/2] から，
// lower: 列 j は AP[jn - j(j-1)/2] から)
inline std::size_t packed_offset_upper(const int j) {
	return static_cast<std::size_t>(j) * (j + 1) / 2;
}

inline std::size_t packed_offset_lower(const int n, const int j) {
	return static_cast<std::size_t>(j) * n - static_cast<std::size_t>(j) * (j - 1) / 2;
}

} // namespace vblas_rdblas_detail

namespace vcp {

// y := alpha*op(A)*x + beta*y, op(A) = A ('N') or A^T ('T','C'), A: m x n
inline void rdgemv(
	const char trans, const int m, const int n,
	const double alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double beta, double* y, const int incy,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool ntrans = det::option_is(trans, 'N');
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::rdblas_error("rdgemv: invalid trans");
	}
	if (m < 0 || n < 0 || lda < std::max(1, m) || incx == 0 || incy == 0) {
		det::rdblas_error("rdgemv: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	const int lenx = ntrans ? n : m;
	const int leny = ntrans ? m : n;
	det::scale_vector(leny, beta, y, incy, fe);
	if (alpha == 0.0) {
		return;
	}
	const int kx = det::vec_start(lenx, incx);
	const int ky = det::vec_start(leny, incy);
	const int threads = det::threads_for_flops(2.0 * m * n);

	if (ntrans) {
		// y_i (行方向) で並列化し，各 thread が全列 j を走査する
#pragma omp parallel num_threads(threads)
		{
			det::RoundingGuard guard(fe);
			const int nth = omp_get_num_threads();
			const int tid = omp_get_thread_num();
			const int chunk = (m + nth - 1) / nth;
			const int i0 = std::min(m, tid * chunk);
			const int i1 = std::min(m, i0 + chunk);
			for (int j = 0; j < n; j++) {
				const double temp = alpha * x[kx + j * incx];
				if (temp == 0.0) {
					continue;
				}
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (incy == 1) {
					for (int i = i0; i < i1; i++) {
						y[i] = std::fma(a[i], temp, y[i]);
					}
				}
				else {
					for (int i = i0; i < i1; i++) {
						y[ky + i * incy] = std::fma(a[i], temp, y[ky + i * incy]);
					}
				}
			}
		}
	}
	else {
#pragma omp parallel num_threads(threads)
		{
			det::RoundingGuard guard(fe);
#pragma omp for schedule(static)
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = 0.0;
				if (incx == 1) {
					for (int i = 0; i < m; i++) {
						temp = std::fma(a[i], x[i], temp);
					}
				}
				else {
					for (int i = 0; i < m; i++) {
						temp = std::fma(a[i], x[kx + i * incx], temp);
					}
				}
				y[ky + j * incy] = std::fma(alpha, temp, y[ky + j * incy]);
			}
		}
	}
}

// y := alpha*op(A)*x + beta*y, A: band 行列 (kl 下帯域, ku 上帯域, band 格納)
inline void rdgbmv(
	const char trans, const int m, const int n, const int kl, const int ku,
	const double alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double beta, double* y, const int incy,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool ntrans = det::option_is(trans, 'N');
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::rdblas_error("rdgbmv: invalid trans");
	}
	if (m < 0 || n < 0 || kl < 0 || ku < 0 || lda < kl + ku + 1 || incx == 0 || incy == 0) {
		det::rdblas_error("rdgbmv: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	const int lenx = ntrans ? n : m;
	const int leny = ntrans ? m : n;
	det::scale_vector(leny, beta, y, incy, fe);
	if (alpha == 0.0) {
		return;
	}
	const int kx = det::vec_start(lenx, incx);
	const int ky = det::vec_start(leny, incy);
	det::RoundingGuard guard(fe);
	for (int j = 0; j < n; j++) {
		const int i0 = std::max(0, j - ku);
		const int i1 = std::min(m - 1, j + kl);
		const double* a = A + static_cast<std::size_t>(lda) * j;
		if (ntrans) {
			const double temp = alpha * x[kx + j * incx];
			if (temp == 0.0) {
				continue;
			}
			for (int i = i0; i <= i1; i++) {
				y[ky + i * incy] = std::fma(a[ku + i - j], temp, y[ky + i * incy]);
			}
		}
		else {
			double temp = 0.0;
			for (int i = i0; i <= i1; i++) {
				temp = std::fma(a[ku + i - j], x[kx + i * incx], temp);
			}
			y[ky + j * incy] = std::fma(alpha, temp, y[ky + j * incy]);
		}
	}
}

// y := alpha*A*x + beta*y, A: n x n 対称 (uplo の三角のみ参照)
inline void rdsymv(
	const char uplo, const int n,
	const double alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double beta, double* y, const int incy,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsymv: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0 || incy == 0) {
		det::rdblas_error("rdsymv: invalid argument");
	}
	if (n == 0 || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	if (alpha == 0.0) {
		det::scale_vector(n, beta, y, incy, fe);
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	const int threads = det::threads_for_flops(2.0 * n * n);

	// 出力行 i ごとの dot 形式 (行ごとに独立なので並列化できる)
#pragma omp parallel num_threads(threads)
	{
		det::RoundingGuard guard(fe);
#pragma omp for schedule(static)
		for (int i = 0; i < n; i++) {
			double temp = 0.0;
			if (upper) {
				const double* ai = A + static_cast<std::size_t>(lda) * i;
				for (int j = 0; j < i; j++) {
					temp = std::fma(ai[j], x[kx + j * incx], temp);
				}
				for (int j = i; j < n; j++) {
					temp = std::fma(A[i + static_cast<std::size_t>(lda) * j], x[kx + j * incx], temp);
				}
			}
			else {
				for (int j = 0; j <= i; j++) {
					temp = std::fma(A[i + static_cast<std::size_t>(lda) * j], x[kx + j * incx], temp);
				}
				const double* ai = A + static_cast<std::size_t>(lda) * i;
				for (int j = i + 1; j < n; j++) {
					temp = std::fma(ai[j], x[kx + j * incx], temp);
				}
			}
			double* yi = y + (ky + i * incy);
			if (beta == 0.0) {
				*yi = alpha * temp;
			}
			else {
				*yi = std::fma(beta, *yi, alpha * temp);
			}
		}
	}
}

// y := alpha*A*x + beta*y, A: n x n 対称 band (帯域 k, band 格納)
inline void rdsbmv(
	const char uplo, const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double beta, double* y, const int incy,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsbmv: invalid uplo");
	}
	if (n < 0 || k < 0 || lda < k + 1 || incx == 0 || incy == 0) {
		det::rdblas_error("rdsbmv: invalid argument");
	}
	if (n == 0 || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	det::scale_vector(n, beta, y, incy, fe);
	if (alpha == 0.0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	det::RoundingGuard guard(fe);
	for (int j = 0; j < n; j++) {
		const double temp1 = alpha * x[kx + j * incx];
		double temp2 = 0.0;
		const double* a = A + static_cast<std::size_t>(lda) * j;
		if (upper) {
			for (int i = std::max(0, j - k); i < j; i++) {
				y[ky + i * incy] = std::fma(a[k + i - j], temp1, y[ky + i * incy]);
				temp2 = std::fma(a[k + i - j], x[kx + i * incx], temp2);
			}
			y[ky + j * incy] = std::fma(temp1, a[k], y[ky + j * incy]);
			y[ky + j * incy] = std::fma(alpha, temp2, y[ky + j * incy]);
		}
		else {
			y[ky + j * incy] = std::fma(temp1, a[0], y[ky + j * incy]);
			for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
				y[ky + i * incy] = std::fma(a[i - j], temp1, y[ky + i * incy]);
				temp2 = std::fma(a[i - j], x[kx + i * incx], temp2);
			}
			y[ky + j * incy] = std::fma(alpha, temp2, y[ky + j * incy]);
		}
	}
}

// y := alpha*A*x + beta*y, A: n x n 対称 packed 格納
inline void rdspmv(
	const char uplo, const int n,
	const double alpha, const double* AP,
	const double* x, const int incx,
	const double beta, double* y, const int incy,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdspmv: invalid uplo");
	}
	if (n < 0 || incx == 0 || incy == 0) {
		det::rdblas_error("rdspmv: invalid argument");
	}
	if (n == 0 || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	det::scale_vector(n, beta, y, incy, fe);
	if (alpha == 0.0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	det::RoundingGuard guard(fe);
	for (int j = 0; j < n; j++) {
		const double temp1 = alpha * x[kx + j * incx];
		double temp2 = 0.0;
		if (upper) {
			const double* ap = AP + det::packed_offset_upper(j);
			for (int i = 0; i < j; i++) {
				y[ky + i * incy] = std::fma(ap[i], temp1, y[ky + i * incy]);
				temp2 = std::fma(ap[i], x[kx + i * incx], temp2);
			}
			y[ky + j * incy] = std::fma(temp1, ap[j], y[ky + j * incy]);
			y[ky + j * incy] = std::fma(alpha, temp2, y[ky + j * incy]);
		}
		else {
			const double* ap = AP + det::packed_offset_lower(n, j);
			y[ky + j * incy] = std::fma(temp1, ap[0], y[ky + j * incy]);
			for (int i = j + 1; i < n; i++) {
				y[ky + i * incy] = std::fma(ap[i - j], temp1, y[ky + i * incy]);
				temp2 = std::fma(ap[i - j], x[kx + i * incx], temp2);
			}
			y[ky + j * incy] = std::fma(alpha, temp2, y[ky + j * incy]);
		}
	}
}

// x := op(A)*x, A: n x n 三角行列
inline void rdtrmv(
	const char uplo, const char trans, const char diag, const int n,
	const double* A, const int lda, double* x, const int incx,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdtrmv: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::rdblas_error("rdtrmv: invalid trans");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::rdblas_error("rdtrmv: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0) {
		det::rdblas_error("rdtrmv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	if (ntrans) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const double temp = x[kx + j * incx];
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != 0.0) {
					for (int i = 0; i < j; i++) {
						x[kx + i * incx] = std::fma(temp, a[i], x[kx + i * incx]);
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * a[j];
				}
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const double temp = x[kx + j * incx];
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != 0.0) {
					for (int i = n - 1; i > j; i--) {
						x[kx + i * incx] = std::fma(temp, a[i], x[kx + i * incx]);
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * a[j];
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[j];
				}
				for (int i = j - 1; i >= 0; i--) {
					temp = std::fma(a[i], x[kx + i * incx], temp);
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[j];
				}
				for (int i = j + 1; i < n; i++) {
					temp = std::fma(a[i], x[kx + i * incx], temp);
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// op(A)*x = b を解いて x := A^{-1} 形の解で上書き，A: n x n 三角行列
inline void rdtrsv(
	const char uplo, const char trans, const char diag, const int n,
	const double* A, const int lda, double* x, const int incx,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdtrsv: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::rdblas_error("rdtrsv: invalid trans");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::rdblas_error("rdtrsv: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0) {
		det::rdblas_error("rdtrsv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	if (ntrans) {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != 0.0) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[j];
					}
					const double temp = x[kx + j * incx];
					for (int i = j - 1; i >= 0; i--) {
						x[kx + i * incx] = std::fma(-temp, a[i], x[kx + i * incx]);
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != 0.0) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[j];
					}
					const double temp = x[kx + j * incx];
					for (int i = j + 1; i < n; i++) {
						x[kx + i * incx] = std::fma(-temp, a[i], x[kx + i * incx]);
					}
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				for (int i = 0; i < j; i++) {
					temp = std::fma(-a[i], x[kx + i * incx], temp);
				}
				if (nounit) {
					temp = temp / a[j];
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				for (int i = n - 1; i > j; i--) {
					temp = std::fma(-a[i], x[kx + i * incx], temp);
				}
				if (nounit) {
					temp = temp / a[j];
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// x := op(A)*x, A: n x n 三角 band 行列 (帯域 k, band 格納)
inline void rdtbmv(
	const char uplo, const char trans, const char diag, const int n, const int k,
	const double* A, const int lda, double* x, const int incx,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || k < 0 || lda < k + 1 || incx == 0) {
		det::rdblas_error("rdtbmv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	if (ntrans) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const double temp = x[kx + j * incx];
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != 0.0) {
					for (int i = std::max(0, j - k); i < j; i++) {
						x[kx + i * incx] = std::fma(temp, a[k + i - j], x[kx + i * incx]);
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * a[k];
				}
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const double temp = x[kx + j * incx];
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != 0.0) {
					for (int i = std::min(n - 1, j + k); i > j; i--) {
						x[kx + i * incx] = std::fma(temp, a[i - j], x[kx + i * incx]);
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * a[0];
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[k];
				}
				for (int i = j - 1; i >= std::max(0, j - k); i--) {
					temp = std::fma(a[k + i - j], x[kx + i * incx], temp);
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[0];
				}
				for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
					temp = std::fma(a[i - j], x[kx + i * incx], temp);
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// op(A)*x = b を解く，A: n x n 三角 band 行列 (帯域 k, band 格納)
inline void rdtbsv(
	const char uplo, const char trans, const char diag, const int n, const int k,
	const double* A, const int lda, double* x, const int incx,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || k < 0 || lda < k + 1 || incx == 0) {
		det::rdblas_error("rdtbsv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	if (ntrans) {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != 0.0) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[k];
					}
					const double temp = x[kx + j * incx];
					for (int i = j - 1; i >= std::max(0, j - k); i--) {
						x[kx + i * incx] = std::fma(-temp, a[k + i - j], x[kx + i * incx]);
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != 0.0) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[0];
					}
					const double temp = x[kx + j * incx];
					for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
						x[kx + i * incx] = std::fma(-temp, a[i - j], x[kx + i * incx]);
					}
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				for (int i = std::max(0, j - k); i < j; i++) {
					temp = std::fma(-a[k + i - j], x[kx + i * incx], temp);
				}
				if (nounit) {
					temp = temp / a[k];
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const double* a = A + static_cast<std::size_t>(lda) * j;
				double temp = x[kx + j * incx];
				for (int i = std::min(n - 1, j + k); i > j; i--) {
					temp = std::fma(-a[i - j], x[kx + i * incx], temp);
				}
				if (nounit) {
					temp = temp / a[0];
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// x := op(A)*x, A: n x n 三角行列 (packed 格納)
inline void rdtpmv(
	const char uplo, const char trans, const char diag, const int n,
	const double* AP, double* x, const int incx,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || incx == 0) {
		det::rdblas_error("rdtpmv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	if (ntrans) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const double temp = x[kx + j * incx];
				const double* ap = AP + det::packed_offset_upper(j);
				if (temp != 0.0) {
					for (int i = 0; i < j; i++) {
						x[kx + i * incx] = std::fma(temp, ap[i], x[kx + i * incx]);
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * ap[j];
				}
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const double temp = x[kx + j * incx];
				const double* ap = AP + det::packed_offset_lower(n, j);
				if (temp != 0.0) {
					for (int i = n - 1; i > j; i--) {
						x[kx + i * incx] = std::fma(temp, ap[i - j], x[kx + i * incx]);
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * ap[0];
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const double* ap = AP + det::packed_offset_upper(j);
				double temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * ap[j];
				}
				for (int i = j - 1; i >= 0; i--) {
					temp = std::fma(ap[i], x[kx + i * incx], temp);
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const double* ap = AP + det::packed_offset_lower(n, j);
				double temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * ap[0];
				}
				for (int i = j + 1; i < n; i++) {
					temp = std::fma(ap[i - j], x[kx + i * incx], temp);
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// op(A)*x = b を解く，A: n x n 三角行列 (packed 格納)
inline void rdtpsv(
	const char uplo, const char trans, const char diag, const int n,
	const double* AP, double* x, const int incx,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || incx == 0) {
		det::rdblas_error("rdtpsv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	if (ntrans) {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const double* ap = AP + det::packed_offset_upper(j);
				if (x[kx + j * incx] != 0.0) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / ap[j];
					}
					const double temp = x[kx + j * incx];
					for (int i = j - 1; i >= 0; i--) {
						x[kx + i * incx] = std::fma(-temp, ap[i], x[kx + i * incx]);
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const double* ap = AP + det::packed_offset_lower(n, j);
				if (x[kx + j * incx] != 0.0) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / ap[0];
					}
					const double temp = x[kx + j * incx];
					for (int i = j + 1; i < n; i++) {
						x[kx + i * incx] = std::fma(-temp, ap[i - j], x[kx + i * incx]);
					}
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const double* ap = AP + det::packed_offset_upper(j);
				double temp = x[kx + j * incx];
				for (int i = 0; i < j; i++) {
					temp = std::fma(-ap[i], x[kx + i * incx], temp);
				}
				if (nounit) {
					temp = temp / ap[j];
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const double* ap = AP + det::packed_offset_lower(n, j);
				double temp = x[kx + j * incx];
				for (int i = n - 1; i > j; i--) {
					temp = std::fma(-ap[i - j], x[kx + i * incx], temp);
				}
				if (nounit) {
					temp = temp / ap[0];
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// A := alpha*x*y^T + A, A: m x n
inline void rdger(
	const int m, const int n,
	const double alpha, const double* x, const int incx,
	const double* y, const int incy, double* A, const int lda,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m) || incx == 0 || incy == 0) {
		det::rdblas_error("rdger: invalid argument");
	}
	if (m == 0 || n == 0 || alpha == 0.0) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	const int kx = det::vec_start(m, incx);
	const int ky = det::vec_start(n, incy);
	const int threads = det::threads_for_flops(2.0 * m * n);
#pragma omp parallel num_threads(threads)
	{
		det::RoundingGuard guard(fe);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			const double temp = alpha * y[ky + j * incy];
			if (temp == 0.0) {
				continue;
			}
			double* a = A + static_cast<std::size_t>(lda) * j;
			if (incx == 1) {
				for (int i = 0; i < m; i++) {
					a[i] = std::fma(x[i], temp, a[i]);
				}
			}
			else {
				for (int i = 0; i < m; i++) {
					a[i] = std::fma(x[kx + i * incx], temp, a[i]);
				}
			}
		}
	}
}

// A := alpha*x*x^T + A, A: n x n 対称 (uplo の三角のみ更新)
inline void rdsyr(
	const char uplo, const int n,
	const double alpha, const double* x, const int incx,
	double* A, const int lda,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsyr: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0) {
		det::rdblas_error("rdsyr: invalid argument");
	}
	if (n == 0 || alpha == 0.0) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	const int kx = det::vec_start(n, incx);
	const int threads = det::threads_for_flops(1.0 * n * n);
#pragma omp parallel num_threads(threads)
	{
		det::RoundingGuard guard(fe);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			const double temp = alpha * x[kx + j * incx];
			if (temp == 0.0) {
				continue;
			}
			double* a = A + static_cast<std::size_t>(lda) * j;
			if (upper) {
				for (int i = 0; i <= j; i++) {
					a[i] = std::fma(x[kx + i * incx], temp, a[i]);
				}
			}
			else {
				for (int i = j; i < n; i++) {
					a[i] = std::fma(x[kx + i * incx], temp, a[i]);
				}
			}
		}
	}
}

// A := alpha*x*x^T + A, A: n x n 対称 packed 格納
inline void rdspr(
	const char uplo, const int n,
	const double alpha, const double* x, const int incx,
	double* AP,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdspr: invalid uplo");
	}
	if (n < 0 || incx == 0) {
		det::rdblas_error("rdspr: invalid argument");
	}
	if (n == 0 || alpha == 0.0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	for (int j = 0; j < n; j++) {
		const double temp = alpha * x[kx + j * incx];
		if (temp == 0.0) {
			continue;
		}
		if (upper) {
			double* ap = AP + det::packed_offset_upper(j);
			for (int i = 0; i <= j; i++) {
				ap[i] = std::fma(x[kx + i * incx], temp, ap[i]);
			}
		}
		else {
			double* ap = AP + det::packed_offset_lower(n, j);
			for (int i = j; i < n; i++) {
				ap[i - j] = std::fma(x[kx + i * incx], temp, ap[i - j]);
			}
		}
	}
}

// A := alpha*x*y^T + alpha*y*x^T + A, A: n x n 対称 (uplo の三角のみ更新)
inline void rdsyr2(
	const char uplo, const int n,
	const double alpha, const double* x, const int incx,
	const double* y, const int incy, double* A, const int lda,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsyr2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0 || incy == 0) {
		det::rdblas_error("rdsyr2: invalid argument");
	}
	if (n == 0 || alpha == 0.0) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	const int threads = det::threads_for_flops(2.0 * n * n);
#pragma omp parallel num_threads(threads)
	{
		det::RoundingGuard guard(fe);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			const double temp1 = alpha * y[ky + j * incy];
			const double temp2 = alpha * x[kx + j * incx];
			if (temp1 == 0.0 && temp2 == 0.0) {
				continue;
			}
			double* a = A + static_cast<std::size_t>(lda) * j;
			if (upper) {
				for (int i = 0; i <= j; i++) {
					a[i] = std::fma(x[kx + i * incx], temp1, a[i]);
					a[i] = std::fma(y[ky + i * incy], temp2, a[i]);
				}
			}
			else {
				for (int i = j; i < n; i++) {
					a[i] = std::fma(x[kx + i * incx], temp1, a[i]);
					a[i] = std::fma(y[ky + i * incy], temp2, a[i]);
				}
			}
		}
	}
}

// A := alpha*x*y^T + alpha*y*x^T + A, A: n x n 対称 packed 格納
inline void rdspr2(
	const char uplo, const int n,
	const double alpha, const double* x, const int incx,
	const double* y, const int incy, double* AP,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdspr2: invalid uplo");
	}
	if (n < 0 || incx == 0 || incy == 0) {
		det::rdblas_error("rdspr2: invalid argument");
	}
	if (n == 0 || alpha == 0.0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	det::RoundingGuard guard(det::fe_rounding(rounding_mode));
	for (int j = 0; j < n; j++) {
		const double temp1 = alpha * y[ky + j * incy];
		const double temp2 = alpha * x[kx + j * incx];
		if (temp1 == 0.0 && temp2 == 0.0) {
			continue;
		}
		if (upper) {
			double* ap = AP + det::packed_offset_upper(j);
			for (int i = 0; i <= j; i++) {
				ap[i] = std::fma(x[kx + i * incx], temp1, ap[i]);
				ap[i] = std::fma(y[ky + i * incy], temp2, ap[i]);
			}
		}
		else {
			double* ap = AP + det::packed_offset_lower(n, j);
			for (int i = j; i < n; i++) {
				ap[i - j] = std::fma(x[kx + i * incx], temp1, ap[i - j]);
				ap[i - j] = std::fma(y[ky + i * incy], temp2, ap[i - j]);
			}
		}
	}
}

} // namespace vcp

#endif // VBLAS_RDBLAS_LEVEL2_HPP
