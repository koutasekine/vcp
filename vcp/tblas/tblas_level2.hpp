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
// Modified by Kouta Sekine (2026): translated to C++ templates over a generic
// element type T (no rounding-mode control, no fma).
// ----------------------------------------------------------------------------

#pragma once

#ifndef TBLAS_TBLAS_LEVEL2_HPP
#define TBLAS_TBLAS_LEVEL2_HPP

// Level 2 BLAS (template 版)．行列は column-major，index は 0-based．
// 引数は reference BLAS と同順 (rounding_mode なし)，scalar は const T& 渡し．
// tgemv / tsymv / tger / tsyr / tsyr2 は演算量が閾値を超えると OpenMP 並列，
// banded / packed / 三角系は逐次．

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "tblas_common.hpp"

namespace vcp {

// y := alpha*op(A)*x + beta*y, op(A) = A ('N') or A^T ('T','C'), A: m x n
template <typename T>
inline void tgemv(
	const char trans, const int m, const int n,
	const T& alpha, const T* A, const int lda,
	const T* x, const int incx,
	const T& beta, T* y, const int incy
) {
	namespace det = tblas_detail;
	const bool ntrans = det::option_is(trans, 'N');
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("tgemv: invalid trans");
	}
	if (m < 0 || n < 0 || lda < std::max(1, m) || incx == 0 || incy == 0) {
		det::tblas_error("tgemv: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == T(0) && beta == T(1))) {
		return;
	}
	const int lenx = ntrans ? n : m;
	const int leny = ntrans ? m : n;
	det::scale_vector(leny, beta, y, incy);
	if (alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(lenx, incx);
	const int ky = det::vec_start(leny, incy);
	const bool par = det::use_parallel(2.0 * m * n);

	if (ntrans) {
		// y_i (行方向) で並列化し，各 thread が全列 j を走査する
#ifdef _OPENMP
#pragma omp parallel if (par)
#endif
		{
			const int nth = det::region_threads();
			const int tid = det::region_thread_id();
			const int chunk = (m + nth - 1) / nth;
			const int i0 = std::min(m, tid * chunk);
			const int i1 = std::min(m, i0 + chunk);
			for (int j = 0; j < n; j++) {
				const T temp = alpha * x[kx + j * incx];
				if (temp == T(0)) {
					continue;
				}
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (incy == 1) {
					for (int i = i0; i < i1; i++) {
						y[i] = a[i] * temp + y[i];
					}
				}
				else {
					for (int i = i0; i < i1; i++) {
						y[ky + i * incy] = a[i] * temp + y[ky + i * incy];
					}
				}
			}
		}
	}
	else {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (par)
#endif
		for (int j = 0; j < n; j++) {
			const T* a = A + static_cast<std::size_t>(lda) * j;
			T temp = T(0);
			if (incx == 1) {
				for (int i = 0; i < m; i++) {
					temp = a[i] * x[i] + temp;
				}
			}
			else {
				for (int i = 0; i < m; i++) {
					temp = a[i] * x[kx + i * incx] + temp;
				}
			}
			y[ky + j * incy] = alpha * temp + y[ky + j * incy];
		}
	}
}

// y := alpha*op(A)*x + beta*y, A: band 行列 (kl 下帯域, ku 上帯域, band 格納)
template <typename T>
inline void tgbmv(
	const char trans, const int m, const int n, const int kl, const int ku,
	const T& alpha, const T* A, const int lda,
	const T* x, const int incx,
	const T& beta, T* y, const int incy
) {
	namespace det = tblas_detail;
	const bool ntrans = det::option_is(trans, 'N');
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("tgbmv: invalid trans");
	}
	if (m < 0 || n < 0 || kl < 0 || ku < 0 || lda < kl + ku + 1 || incx == 0 || incy == 0) {
		det::tblas_error("tgbmv: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == T(0) && beta == T(1))) {
		return;
	}
	const int lenx = ntrans ? n : m;
	const int leny = ntrans ? m : n;
	det::scale_vector(leny, beta, y, incy);
	if (alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(lenx, incx);
	const int ky = det::vec_start(leny, incy);
	for (int j = 0; j < n; j++) {
		const int i0 = std::max(0, j - ku);
		const int i1 = std::min(m - 1, j + kl);
		const T* a = A + static_cast<std::size_t>(lda) * j;
		if (ntrans) {
			const T temp = alpha * x[kx + j * incx];
			if (temp == T(0)) {
				continue;
			}
			for (int i = i0; i <= i1; i++) {
				y[ky + i * incy] = a[ku + i - j] * temp + y[ky + i * incy];
			}
		}
		else {
			T temp = T(0);
			for (int i = i0; i <= i1; i++) {
				temp = a[ku + i - j] * x[kx + i * incx] + temp;
			}
			y[ky + j * incy] = alpha * temp + y[ky + j * incy];
		}
	}
}

// y := alpha*A*x + beta*y, A: n x n 対称 (uplo の三角のみ参照)
template <typename T>
inline void tsymv(
	const char uplo, const int n,
	const T& alpha, const T* A, const int lda,
	const T* x, const int incx,
	const T& beta, T* y, const int incy
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsymv: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0 || incy == 0) {
		det::tblas_error("tsymv: invalid argument");
	}
	if (n == 0 || (alpha == T(0) && beta == T(1))) {
		return;
	}
	if (alpha == T(0)) {
		det::scale_vector(n, beta, y, incy);
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	const bool zero_beta = (beta == T(0));

	// 出力行 i ごとの dot 形式 (行ごとに独立なので並列化できる)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(2.0 * n * n))
#endif
	for (int i = 0; i < n; i++) {
		T temp = T(0);
		if (upper) {
			const T* ai = A + static_cast<std::size_t>(lda) * i;
			for (int j = 0; j < i; j++) {
				temp = ai[j] * x[kx + j * incx] + temp;
			}
			for (int j = i; j < n; j++) {
				temp = A[i + static_cast<std::size_t>(lda) * j] * x[kx + j * incx] + temp;
			}
		}
		else {
			for (int j = 0; j <= i; j++) {
				temp = A[i + static_cast<std::size_t>(lda) * j] * x[kx + j * incx] + temp;
			}
			const T* ai = A + static_cast<std::size_t>(lda) * i;
			for (int j = i + 1; j < n; j++) {
				temp = ai[j] * x[kx + j * incx] + temp;
			}
		}
		T* yi = y + (ky + i * incy);
		if (zero_beta) {
			*yi = alpha * temp;
		}
		else {
			*yi = beta * *yi + alpha * temp;
		}
	}
}

// y := alpha*A*x + beta*y, A: n x n 対称 band (帯域 k, band 格納)
template <typename T>
inline void tsbmv(
	const char uplo, const int n, const int k,
	const T& alpha, const T* A, const int lda,
	const T* x, const int incx,
	const T& beta, T* y, const int incy
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsbmv: invalid uplo");
	}
	if (n < 0 || k < 0 || lda < k + 1 || incx == 0 || incy == 0) {
		det::tblas_error("tsbmv: invalid argument");
	}
	if (n == 0 || (alpha == T(0) && beta == T(1))) {
		return;
	}
	det::scale_vector(n, beta, y, incy);
	if (alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	for (int j = 0; j < n; j++) {
		const T temp1 = alpha * x[kx + j * incx];
		T temp2 = T(0);
		const T* a = A + static_cast<std::size_t>(lda) * j;
		if (upper) {
			for (int i = std::max(0, j - k); i < j; i++) {
				y[ky + i * incy] = a[k + i - j] * temp1 + y[ky + i * incy];
				temp2 = a[k + i - j] * x[kx + i * incx] + temp2;
			}
			y[ky + j * incy] = temp1 * a[k] + y[ky + j * incy];
			y[ky + j * incy] = alpha * temp2 + y[ky + j * incy];
		}
		else {
			y[ky + j * incy] = temp1 * a[0] + y[ky + j * incy];
			for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
				y[ky + i * incy] = a[i - j] * temp1 + y[ky + i * incy];
				temp2 = a[i - j] * x[kx + i * incx] + temp2;
			}
			y[ky + j * incy] = alpha * temp2 + y[ky + j * incy];
		}
	}
}

// y := alpha*A*x + beta*y, A: n x n 対称 packed 格納
template <typename T>
inline void tspmv(
	const char uplo, const int n,
	const T& alpha, const T* AP,
	const T* x, const int incx,
	const T& beta, T* y, const int incy
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tspmv: invalid uplo");
	}
	if (n < 0 || incx == 0 || incy == 0) {
		det::tblas_error("tspmv: invalid argument");
	}
	if (n == 0 || (alpha == T(0) && beta == T(1))) {
		return;
	}
	det::scale_vector(n, beta, y, incy);
	if (alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	for (int j = 0; j < n; j++) {
		const T temp1 = alpha * x[kx + j * incx];
		T temp2 = T(0);
		if (upper) {
			const T* ap = AP + det::packed_offset_upper(j);
			for (int i = 0; i < j; i++) {
				y[ky + i * incy] = ap[i] * temp1 + y[ky + i * incy];
				temp2 = ap[i] * x[kx + i * incx] + temp2;
			}
			y[ky + j * incy] = temp1 * ap[j] + y[ky + j * incy];
			y[ky + j * incy] = alpha * temp2 + y[ky + j * incy];
		}
		else {
			const T* ap = AP + det::packed_offset_lower(n, j);
			y[ky + j * incy] = temp1 * ap[0] + y[ky + j * incy];
			for (int i = j + 1; i < n; i++) {
				y[ky + i * incy] = ap[i - j] * temp1 + y[ky + i * incy];
				temp2 = ap[i - j] * x[kx + i * incx] + temp2;
			}
			y[ky + j * incy] = alpha * temp2 + y[ky + j * incy];
		}
	}
}

// x := op(A)*x, A: n x n 三角行列
template <typename T>
inline void ttrmv(
	const char uplo, const char trans, const char diag, const int n,
	const T* A, const int lda, T* x, const int incx
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("ttrmv: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("ttrmv: invalid trans");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::tblas_error("ttrmv: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0) {
		det::tblas_error("ttrmv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	if (ntrans) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const T temp = x[kx + j * incx];
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != T(0)) {
					for (int i = 0; i < j; i++) {
						x[kx + i * incx] = temp * a[i] + x[kx + i * incx];
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * a[j];
				}
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const T temp = x[kx + j * incx];
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != T(0)) {
					for (int i = n - 1; i > j; i--) {
						x[kx + i * incx] = temp * a[i] + x[kx + i * incx];
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
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[j];
				}
				for (int i = j - 1; i >= 0; i--) {
					temp = a[i] * x[kx + i * incx] + temp;
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[j];
				}
				for (int i = j + 1; i < n; i++) {
					temp = a[i] * x[kx + i * incx] + temp;
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// op(A)*x = b を解いて x を解で上書き，A: n x n 三角行列
template <typename T>
inline void ttrsv(
	const char uplo, const char trans, const char diag, const int n,
	const T* A, const int lda, T* x, const int incx
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("ttrsv: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("ttrsv: invalid trans");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::tblas_error("ttrsv: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0) {
		det::tblas_error("ttrsv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	if (ntrans) {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != T(0)) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[j];
					}
					const T temp = x[kx + j * incx];
					for (int i = j - 1; i >= 0; i--) {
						x[kx + i * incx] = x[kx + i * incx] - temp * a[i];
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != T(0)) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[j];
					}
					const T temp = x[kx + j * incx];
					for (int i = j + 1; i < n; i++) {
						x[kx + i * incx] = x[kx + i * incx] - temp * a[i];
					}
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				for (int i = 0; i < j; i++) {
					temp = temp - a[i] * x[kx + i * incx];
				}
				if (nounit) {
					temp = temp / a[j];
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				for (int i = n - 1; i > j; i--) {
					temp = temp - a[i] * x[kx + i * incx];
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
template <typename T>
inline void ttbmv(
	const char uplo, const char trans, const char diag, const int n, const int k,
	const T* A, const int lda, T* x, const int incx
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || k < 0 || lda < k + 1 || incx == 0) {
		det::tblas_error("ttbmv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	if (ntrans) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const T temp = x[kx + j * incx];
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != T(0)) {
					for (int i = std::max(0, j - k); i < j; i++) {
						x[kx + i * incx] = temp * a[k + i - j] + x[kx + i * incx];
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * a[k];
				}
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const T temp = x[kx + j * incx];
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (temp != T(0)) {
					for (int i = std::min(n - 1, j + k); i > j; i--) {
						x[kx + i * incx] = temp * a[i - j] + x[kx + i * incx];
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
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[k];
				}
				for (int i = j - 1; i >= std::max(0, j - k); i--) {
					temp = a[k + i - j] * x[kx + i * incx] + temp;
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * a[0];
				}
				for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
					temp = a[i - j] * x[kx + i * incx] + temp;
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// op(A)*x = b を解く，A: n x n 三角 band 行列 (帯域 k, band 格納)
template <typename T>
inline void ttbsv(
	const char uplo, const char trans, const char diag, const int n, const int k,
	const T* A, const int lda, T* x, const int incx
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || k < 0 || lda < k + 1 || incx == 0) {
		det::tblas_error("ttbsv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	if (ntrans) {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != T(0)) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[k];
					}
					const T temp = x[kx + j * incx];
					for (int i = j - 1; i >= std::max(0, j - k); i--) {
						x[kx + i * incx] = x[kx + i * incx] - temp * a[k + i - j];
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				if (x[kx + j * incx] != T(0)) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / a[0];
					}
					const T temp = x[kx + j * incx];
					for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
						x[kx + i * incx] = x[kx + i * incx] - temp * a[i - j];
					}
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				for (int i = std::max(0, j - k); i < j; i++) {
					temp = temp - a[k + i - j] * x[kx + i * incx];
				}
				if (nounit) {
					temp = temp / a[k];
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const T* a = A + static_cast<std::size_t>(lda) * j;
				T temp = x[kx + j * incx];
				for (int i = std::min(n - 1, j + k); i > j; i--) {
					temp = temp - a[i - j] * x[kx + i * incx];
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
template <typename T>
inline void ttpmv(
	const char uplo, const char trans, const char diag, const int n,
	const T* AP, T* x, const int incx
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || incx == 0) {
		det::tblas_error("ttpmv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	if (ntrans) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const T temp = x[kx + j * incx];
				const T* ap = AP + det::packed_offset_upper(j);
				if (temp != T(0)) {
					for (int i = 0; i < j; i++) {
						x[kx + i * incx] = temp * ap[i] + x[kx + i * incx];
					}
				}
				if (nounit) {
					x[kx + j * incx] = temp * ap[j];
				}
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const T temp = x[kx + j * incx];
				const T* ap = AP + det::packed_offset_lower(n, j);
				if (temp != T(0)) {
					for (int i = n - 1; i > j; i--) {
						x[kx + i * incx] = temp * ap[i - j] + x[kx + i * incx];
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
				const T* ap = AP + det::packed_offset_upper(j);
				T temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * ap[j];
				}
				for (int i = j - 1; i >= 0; i--) {
					temp = ap[i] * x[kx + i * incx] + temp;
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const T* ap = AP + det::packed_offset_lower(n, j);
				T temp = x[kx + j * incx];
				if (nounit) {
					temp = temp * ap[0];
				}
				for (int i = j + 1; i < n; i++) {
					temp = ap[i - j] * x[kx + i * incx] + temp;
				}
				x[kx + j * incx] = temp;
			}
		}
	}
}

// op(A)*x = b を解く，A: n x n 三角行列 (packed 格納)
template <typename T>
inline void ttpsv(
	const char uplo, const char trans, const char diag, const int n,
	const T* AP, T* x, const int incx
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if ((!upper && !det::option_is(uplo, 'L')) ||
	    (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) ||
	    (!nounit && !det::option_is(diag, 'U')) ||
	    n < 0 || incx == 0) {
		det::tblas_error("ttpsv: invalid argument");
	}
	if (n == 0) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	if (ntrans) {
		if (upper) {
			for (int j = n - 1; j >= 0; j--) {
				const T* ap = AP + det::packed_offset_upper(j);
				if (x[kx + j * incx] != T(0)) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / ap[j];
					}
					const T temp = x[kx + j * incx];
					for (int i = j - 1; i >= 0; i--) {
						x[kx + i * incx] = x[kx + i * incx] - temp * ap[i];
					}
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				const T* ap = AP + det::packed_offset_lower(n, j);
				if (x[kx + j * incx] != T(0)) {
					if (nounit) {
						x[kx + j * incx] = x[kx + j * incx] / ap[0];
					}
					const T temp = x[kx + j * incx];
					for (int i = j + 1; i < n; i++) {
						x[kx + i * incx] = x[kx + i * incx] - temp * ap[i - j];
					}
				}
			}
		}
	}
	else {
		if (upper) {
			for (int j = 0; j < n; j++) {
				const T* ap = AP + det::packed_offset_upper(j);
				T temp = x[kx + j * incx];
				for (int i = 0; i < j; i++) {
					temp = temp - ap[i] * x[kx + i * incx];
				}
				if (nounit) {
					temp = temp / ap[j];
				}
				x[kx + j * incx] = temp;
			}
		}
		else {
			for (int j = n - 1; j >= 0; j--) {
				const T* ap = AP + det::packed_offset_lower(n, j);
				T temp = x[kx + j * incx];
				for (int i = n - 1; i > j; i--) {
					temp = temp - ap[i - j] * x[kx + i * incx];
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
template <typename T>
inline void tger(
	const int m, const int n,
	const T& alpha, const T* x, const int incx,
	const T* y, const int incy, T* A, const int lda
) {
	namespace det = tblas_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m) || incx == 0 || incy == 0) {
		det::tblas_error("tger: invalid argument");
	}
	if (m == 0 || n == 0 || alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(m, incx);
	const int ky = det::vec_start(n, incy);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(2.0 * m * n))
#endif
	for (int j = 0; j < n; j++) {
		const T temp = alpha * y[ky + j * incy];
		if (temp == T(0)) {
			continue;
		}
		T* a = A + static_cast<std::size_t>(lda) * j;
		if (incx == 1) {
			for (int i = 0; i < m; i++) {
				a[i] = x[i] * temp + a[i];
			}
		}
		else {
			for (int i = 0; i < m; i++) {
				a[i] = x[kx + i * incx] * temp + a[i];
			}
		}
	}
}

// A := alpha*x*x^T + A, A: n x n 対称 (uplo の三角のみ更新)
template <typename T>
inline void tsyr(
	const char uplo, const int n,
	const T& alpha, const T* x, const int incx,
	T* A, const int lda
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsyr: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0) {
		det::tblas_error("tsyr: invalid argument");
	}
	if (n == 0 || alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(n, incx);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(1.0 * n * n))
#endif
	for (int j = 0; j < n; j++) {
		const T temp = alpha * x[kx + j * incx];
		if (temp == T(0)) {
			continue;
		}
		T* a = A + static_cast<std::size_t>(lda) * j;
		if (upper) {
			for (int i = 0; i <= j; i++) {
				a[i] = x[kx + i * incx] * temp + a[i];
			}
		}
		else {
			for (int i = j; i < n; i++) {
				a[i] = x[kx + i * incx] * temp + a[i];
			}
		}
	}
}

// A := alpha*x*x^T + A, A: n x n 対称 packed 格納
template <typename T>
inline void tspr(
	const char uplo, const int n,
	const T& alpha, const T* x, const int incx,
	T* AP
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tspr: invalid uplo");
	}
	if (n < 0 || incx == 0) {
		det::tblas_error("tspr: invalid argument");
	}
	if (n == 0 || alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	for (int j = 0; j < n; j++) {
		const T temp = alpha * x[kx + j * incx];
		if (temp == T(0)) {
			continue;
		}
		if (upper) {
			T* ap = AP + det::packed_offset_upper(j);
			for (int i = 0; i <= j; i++) {
				ap[i] = x[kx + i * incx] * temp + ap[i];
			}
		}
		else {
			T* ap = AP + det::packed_offset_lower(n, j);
			for (int i = j; i < n; i++) {
				ap[i - j] = x[kx + i * incx] * temp + ap[i - j];
			}
		}
	}
}

// A := alpha*x*y^T + alpha*y*x^T + A, A: n x n 対称 (uplo の三角のみ更新)
template <typename T>
inline void tsyr2(
	const char uplo, const int n,
	const T& alpha, const T* x, const int incx,
	const T* y, const int incy, T* A, const int lda
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsyr2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || incx == 0 || incy == 0) {
		det::tblas_error("tsyr2: invalid argument");
	}
	if (n == 0 || alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(2.0 * n * n))
#endif
	for (int j = 0; j < n; j++) {
		const T temp1 = alpha * y[ky + j * incy];
		const T temp2 = alpha * x[kx + j * incx];
		if (temp1 == T(0) && temp2 == T(0)) {
			continue;
		}
		T* a = A + static_cast<std::size_t>(lda) * j;
		if (upper) {
			for (int i = 0; i <= j; i++) {
				a[i] = x[kx + i * incx] * temp1 + a[i];
				a[i] = y[ky + i * incy] * temp2 + a[i];
			}
		}
		else {
			for (int i = j; i < n; i++) {
				a[i] = x[kx + i * incx] * temp1 + a[i];
				a[i] = y[ky + i * incy] * temp2 + a[i];
			}
		}
	}
}

// A := alpha*x*y^T + alpha*y*x^T + A, A: n x n 対称 packed 格納
template <typename T>
inline void tspr2(
	const char uplo, const int n,
	const T& alpha, const T* x, const int incx,
	const T* y, const int incy, T* AP
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tspr2: invalid uplo");
	}
	if (n < 0 || incx == 0 || incy == 0) {
		det::tblas_error("tspr2: invalid argument");
	}
	if (n == 0 || alpha == T(0)) {
		return;
	}
	const int kx = det::vec_start(n, incx);
	const int ky = det::vec_start(n, incy);
	for (int j = 0; j < n; j++) {
		const T temp1 = alpha * y[ky + j * incy];
		const T temp2 = alpha * x[kx + j * incx];
		if (temp1 == T(0) && temp2 == T(0)) {
			continue;
		}
		if (upper) {
			T* ap = AP + det::packed_offset_upper(j);
			for (int i = 0; i <= j; i++) {
				ap[i] = x[kx + i * incx] * temp1 + ap[i];
				ap[i] = y[ky + i * incy] * temp2 + ap[i];
			}
		}
		else {
			T* ap = AP + det::packed_offset_lower(n, j);
			for (int i = j; i < n; i++) {
				ap[i - j] = x[kx + i * incx] * temp1 + ap[i - j];
				ap[i - j] = y[ky + i * incy] * temp2 + ap[i - j];
			}
		}
	}
}

} // namespace vcp

#endif // TBLAS_TBLAS_LEVEL2_HPP
