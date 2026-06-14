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

#pragma once

#ifndef TBLAS_TBLAS_COMMON_HPP
#define TBLAS_TBLAS_COMMON_HPP

// tblas (template BLAS) の共通基盤．
// - 要素型 T の要件 (R1: 値セマンティクス / T(int) 構築 / + - * 単項- == != ,
//   R2: 除算, R3: 順序比較, R4: ADL の abs / sqrt) は tblas/README_tblas.md 参照
// - 丸めモードは一切扱わない (vblas/rdblas との最大の違い)
// - 行列は column-major，index は 0-based
// - OpenMP が有効な場合，演算量が閾値を超える独立 loop のみ並列化する
//   (T の演算が thread-safe な値計算であることを前提とする)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

// OpenMP 並列化を発動するスカラー演算回数の閾値 (T の 1 演算が重い型を想定して低め)
#ifndef TBLAS_OMP_THRESHOLD
#define TBLAS_OMP_THRESHOLD 20000.0
#endif

namespace tblas_detail {

inline void tblas_error(const char* message) {
#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
	throw std::invalid_argument(message);
#else
	(void)message;
	std::abort();
#endif
}

// BLAS の option 文字比較 (lsame 相当，target は大文字で渡す)
inline bool option_is(const char option, const char target) {
	return option == target || option == target + ('a' - 'A');
}

// 増分 incx の vector の先頭 index (BLAS の負増分規約)
inline int vec_start(const int n, const int inc) {
	return inc > 0 ? 0 : (n - 1) * (-inc);
}

// packed 格納の列先頭 offset (upper: 列 j は AP[j(j+1)/2] から，
// lower: 列 j は AP[jn - j(j-1)/2] から)
inline std::size_t packed_offset_upper(const int j) {
	return static_cast<std::size_t>(j) * (j + 1) / 2;
}

inline std::size_t packed_offset_lower(const int n, const int j) {
	return static_cast<std::size_t>(j) * n - static_cast<std::size_t>(j) * (j - 1) / 2;
}

// 演算量 ops が閾値を超えたら OpenMP 並列化する (if 句用)
inline bool use_parallel(const double ops) {
#ifdef _OPENMP
	return ops >= TBLAS_OMP_THRESHOLD;
#else
	(void)ops;
	return false;
#endif
}

// OpenMP 無効時にも compile できるようにする shim
inline int region_threads() {
#ifdef _OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

inline int region_thread_id() {
#ifdef _OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

// y(leny, incy 付き) := beta*y (beta==0 は exact 0 埋め: BLAS 規約で y は読まない)
template <typename T>
inline void scale_vector(const int leny, const T& beta, T* y, const int incy) {
	if (leny <= 0 || beta == T(1)) {
		return;
	}
	int iy = vec_start(leny, incy);
	if (beta == T(0)) {
		for (int i = 0; i < leny; i++) {
			y[iy] = T(0);
			iy += incy;
		}
		return;
	}
	for (int i = 0; i < leny; i++) {
		y[iy] = beta * y[iy];
		iy += incy;
	}
}

// C(m x n, ldc 付き) := beta*C (beta==0 は exact 0 埋め，beta==1 は何もしない)
template <typename T>
inline void scale_matrix(const int m, const int n, const T& beta, T* C, const int ldc) {
	if (m <= 0 || n <= 0 || beta == T(1)) {
		return;
	}
	const bool zero = (beta == T(0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_parallel(static_cast<double>(m) * n))
#endif
	for (int j = 0; j < n; j++) {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		if (zero) {
			for (int i = 0; i < m; i++) {
				c[i] = T(0);
			}
		}
		else {
			for (int i = 0; i < m; i++) {
				c[i] = beta * c[i];
			}
		}
	}
}

// 三角部分のみ C := beta*C (beta==0 は exact 0 埋め)
template <typename T>
inline void scale_triangle(const bool upper, const int n, const T& beta, T* C, const int ldc) {
	if (n <= 0 || beta == T(1)) {
		return;
	}
	const bool zero = (beta == T(0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (use_parallel(0.5 * n * n))
#endif
	for (int j = 0; j < n; j++) {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		const int i0 = upper ? 0 : j;
		const int i1 = upper ? j + 1 : n;
		if (zero) {
			for (int i = i0; i < i1; i++) {
				c[i] = T(0);
			}
		}
		else {
			for (int i = i0; i < i1; i++) {
				c[i] = beta * c[i];
			}
		}
	}
}

// 対称行列 (uplo の三角のみ格納) を full dense (n x n 連続) に展開する (copy のみ)
template <typename T>
inline void expand_symmetric(const bool upper, const int n, const T* A, const int lda, T* dst) {
	for (int j = 0; j < n; j++) {
		T* d = dst + static_cast<std::size_t>(n) * j;
		if (upper) {
			const T* a = A + static_cast<std::size_t>(lda) * j;
			for (int i = 0; i <= j; i++) {
				d[i] = a[i];
			}
			for (int i = j + 1; i < n; i++) {
				d[i] = A[j + static_cast<std::size_t>(lda) * i];
			}
		}
		else {
			for (int i = 0; i < j; i++) {
				d[i] = A[j + static_cast<std::size_t>(lda) * i];
			}
			const T* a = A + static_cast<std::size_t>(lda) * j;
			for (int i = j; i < n; i++) {
				d[i] = a[i];
			}
		}
	}
}

// 三角行列を 0 詰め (unit diag は 1) の full dense (n x n 連続) に展開する (copy のみ)
template <typename T>
inline void expand_triangular(const bool upper, const bool nounit, const int n, const T* A, const int lda, T* dst) {
	for (int j = 0; j < n; j++) {
		const T* a = A + static_cast<std::size_t>(lda) * j;
		T* d = dst + static_cast<std::size_t>(n) * j;
		if (upper) {
			for (int i = 0; i < j; i++) {
				d[i] = a[i];
			}
			d[j] = nounit ? a[j] : T(1);
			for (int i = j + 1; i < n; i++) {
				d[i] = T(0);
			}
		}
		else {
			for (int i = 0; i < j; i++) {
				d[i] = T(0);
			}
			d[j] = nounit ? a[j] : T(1);
			for (int i = j + 1; i < n; i++) {
				d[i] = a[i];
			}
		}
	}
}

} // namespace tblas_detail

#endif // TBLAS_TBLAS_COMMON_HPP
