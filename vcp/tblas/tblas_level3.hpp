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
// element type T (no rounding-mode control, no fma, no SIMD kernels).
// ----------------------------------------------------------------------------

#pragma once

#ifndef TBLAS_TBLAS_LEVEL3_HPP
#define TBLAS_TBLAS_LEVEL3_HPP

// Level 3 BLAS (template 版)．行列は column-major．
// 引数は reference BLAS と同順 (rounding_mode なし)，scalar は const T& 渡し．
//
// 実装方針:
// - tgemm は reference dgemm と同じ列方向 algorithm (列 j ごとに独立なので並列化)
// - tsymm / ttrmm は対称・三角格納を full dense に展開して (copy のみ，
//   unit diag は 1, 三角の外は exact 0) tgemm に帰着する (rdblas と同方式)
// - tsyrk / tsyr2k / tgemmtr は三角部分のみを dot 形式で直接計算する
// - ttrsm は vblas/rdblas の leaf solver (列 sweep) の非再帰移植
// - O(N^3) loop は列方向に独立なものだけ OpenMP 並列化する

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "tblas_common.hpp"

namespace vcp {

// C := alpha*op(A)*op(B) + beta*C
// op(A): m x k, op(B): k x n, C: m x n
template <typename T>
inline void tgemm(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const T& alpha, const T* A, const int lda,
	const T* B, const int ldb,
	const T& beta, T* C, const int ldc
) {
	namespace det = tblas_detail;
	const bool ta = !det::option_is(transa, 'N');
	const bool tb = !det::option_is(transb, 'N');
	if (ta && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("tgemm: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::tblas_error("tgemm: invalid transb");
	}
	const int nrowa = ta ? k : m;
	const int nrowb = tb ? n : k;
	if (m < 0 || n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, m)) {
		det::tblas_error("tgemm: invalid argument");
	}
	if (m == 0 || n == 0 || (k == 0 && beta == T(1)) || (alpha == T(0) && beta == T(1))) {
		return;
	}
	if (alpha == T(0) || k == 0) {
		det::scale_matrix(m, n, beta, C, ldc);
		return;
	}
	const bool zero_beta = (beta == T(0));
	const bool one_beta = (beta == T(1));
	const bool par = det::use_parallel(2.0 * m * n * k);

	if (!ta) {
		// C(:,j) := beta*C(:,j) + sum_l (alpha*opB(l,j)) * A(:,l)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (par)
#endif
		for (int j = 0; j < n; j++) {
			T* c = C + static_cast<std::size_t>(ldc) * j;
			if (zero_beta) {
				for (int i = 0; i < m; i++) {
					c[i] = T(0);
				}
			}
			else if (!one_beta) {
				for (int i = 0; i < m; i++) {
					c[i] = beta * c[i];
				}
			}
			for (int l = 0; l < k; l++) {
				const T blj = tb ? B[j + static_cast<std::size_t>(ldb) * l] : B[l + static_cast<std::size_t>(ldb) * j];
				if (blj == T(0)) {
					continue;
				}
				const T temp = alpha * blj;
				const T* a = A + static_cast<std::size_t>(lda) * l;
				for (int i = 0; i < m; i++) {
					c[i] = temp * a[i] + c[i];
				}
			}
		}
	}
	else {
		// C(i,j) := alpha * (A(:,i)^T . opB(:,j)) + beta*C(i,j)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (par)
#endif
		for (int j = 0; j < n; j++) {
			T* c = C + static_cast<std::size_t>(ldc) * j;
			for (int i = 0; i < m; i++) {
				const T* ai = A + static_cast<std::size_t>(lda) * i;
				T temp = T(0);
				if (!tb) {
					const T* b = B + static_cast<std::size_t>(ldb) * j;
					for (int l = 0; l < k; l++) {
						temp = ai[l] * b[l] + temp;
					}
				}
				else {
					for (int l = 0; l < k; l++) {
						temp = ai[l] * B[j + static_cast<std::size_t>(ldb) * l] + temp;
					}
				}
				if (zero_beta) {
					c[i] = alpha * temp;
				}
				else {
					c[i] = beta * c[i] + alpha * temp;
				}
			}
		}
	}
}

// C := alpha*A*B + beta*C (side 'L') または C := alpha*B*A + beta*C (side 'R')
// A は対称行列 (uplo の三角のみ参照)，C: m x n
template <typename T>
inline void tsymm(
	const char side, const char uplo,
	const int m, const int n,
	const T& alpha, const T* A, const int lda,
	const T* B, const int ldb,
	const T& beta, T* C, const int ldc
) {
	namespace det = tblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	if (!lside && !det::option_is(side, 'R')) {
		det::tblas_error("tsymm: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsymm: invalid uplo");
	}
	const int ka = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, ka) || ldb < std::max(1, m) || ldc < std::max(1, m)) {
		det::tblas_error("tsymm: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == T(0) && beta == T(1))) {
		return;
	}
	if (alpha == T(0)) {
		det::scale_matrix(m, n, beta, C, ldc);
		return;
	}

	// 対称三角格納を full dense に展開 (copy のみ) して tgemm に帰着する
	std::vector<T> Afull(static_cast<std::size_t>(ka) * ka, T(0));
	det::expand_symmetric(upper, ka, A, lda, Afull.data());
	if (lside) {
		tgemm('N', 'N', m, n, m, alpha, Afull.data(), m, B, ldb, beta, C, ldc);
	}
	else {
		tgemm('N', 'N', m, n, n, alpha, B, ldb, Afull.data(), n, beta, C, ldc);
	}
}

// C := alpha*op(A)*op(A)^T + beta*C (uplo の三角のみ更新)
// trans 'N': op(A) = A (n x k), trans 'T'/'C': op(A) = A^T (A は k x n)
template <typename T>
inline void tsyrk(
	const char uplo, const char trans,
	const int n, const int k,
	const T& alpha, const T* A, const int lda,
	const T& beta, T* C, const int ldc
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsyrk: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("tsyrk: invalid trans");
	}
	const int nrowa = ntrans ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldc < std::max(1, n)) {
		det::tblas_error("tsyrk: invalid argument");
	}
	if (n == 0 || (alpha == T(0) && beta == T(1)) || (k == 0 && beta == T(1))) {
		return;
	}
	if (alpha == T(0) || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc);
		return;
	}
	const bool zero_beta = (beta == T(0));

	// 三角部分のみ dot 形式: C(i,j) = alpha * sum_l opA(i,l)*opA(j,l) + beta*C(i,j)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(1.0 * n * n * k))
#endif
	for (int j = 0; j < n; j++) {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		const int i0 = upper ? 0 : j;
		const int i1 = upper ? j + 1 : n;
		for (int i = i0; i < i1; i++) {
			T temp = T(0);
			if (ntrans) {
				for (int l = 0; l < k; l++) {
					temp = A[i + static_cast<std::size_t>(lda) * l] * A[j + static_cast<std::size_t>(lda) * l] + temp;
				}
			}
			else {
				const T* ai = A + static_cast<std::size_t>(lda) * i;
				const T* aj = A + static_cast<std::size_t>(lda) * j;
				for (int l = 0; l < k; l++) {
					temp = ai[l] * aj[l] + temp;
				}
			}
			if (zero_beta) {
				c[i] = alpha * temp;
			}
			else {
				c[i] = beta * c[i] + alpha * temp;
			}
		}
	}
}

// C := alpha*op(A)*op(B)^T + alpha*op(B)*op(A)^T + beta*C (uplo の三角のみ更新)
template <typename T>
inline void tsyr2k(
	const char uplo, const char trans,
	const int n, const int k,
	const T& alpha, const T* A, const int lda,
	const T* B, const int ldb,
	const T& beta, T* C, const int ldc
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsyr2k: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("tsyr2k: invalid trans");
	}
	const int nrowa = ntrans ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowa) || ldc < std::max(1, n)) {
		det::tblas_error("tsyr2k: invalid argument");
	}
	if (n == 0 || (alpha == T(0) && beta == T(1)) || (k == 0 && beta == T(1))) {
		return;
	}
	if (alpha == T(0) || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc);
		return;
	}
	const bool zero_beta = (beta == T(0));

	// 三角部分のみ dot 形式:
	// C(i,j) = alpha*sum_l opA(i,l)*opB(j,l) + alpha*sum_l opB(i,l)*opA(j,l) + beta*C(i,j)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(2.0 * n * n * k))
#endif
	for (int j = 0; j < n; j++) {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		const int i0 = upper ? 0 : j;
		const int i1 = upper ? j + 1 : n;
		for (int i = i0; i < i1; i++) {
			T temp1 = T(0);
			T temp2 = T(0);
			if (ntrans) {
				for (int l = 0; l < k; l++) {
					const T ail = A[i + static_cast<std::size_t>(lda) * l];
					const T ajl = A[j + static_cast<std::size_t>(lda) * l];
					const T bil = B[i + static_cast<std::size_t>(ldb) * l];
					const T bjl = B[j + static_cast<std::size_t>(ldb) * l];
					temp1 = ail * bjl + temp1;
					temp2 = bil * ajl + temp2;
				}
			}
			else {
				const T* ai = A + static_cast<std::size_t>(lda) * i;
				const T* aj = A + static_cast<std::size_t>(lda) * j;
				const T* bi = B + static_cast<std::size_t>(ldb) * i;
				const T* bj = B + static_cast<std::size_t>(ldb) * j;
				for (int l = 0; l < k; l++) {
					temp1 = ai[l] * bj[l] + temp1;
					temp2 = bi[l] * aj[l] + temp2;
				}
			}
			if (zero_beta) {
				c[i] = alpha * temp1 + alpha * temp2;
			}
			else {
				c[i] = beta * c[i] + (alpha * temp1 + alpha * temp2);
			}
		}
	}
}

// B := alpha*op(A)*B (side 'L') または B := alpha*B*op(A) (side 'R')
// A は三角行列，B: m x n
template <typename T>
inline void ttrmm(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const T& alpha, const T* A, const int lda,
	T* B, const int ldb
) {
	namespace det = tblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(transa, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!lside && !det::option_is(side, 'R')) {
		det::tblas_error("ttrmm: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("ttrmm: invalid uplo");
	}
	if (!ntrans && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("ttrmm: invalid transa");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::tblas_error("ttrmm: invalid diag");
	}
	const int nrowa = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, m)) {
		det::tblas_error("ttrmm: invalid argument");
	}
	if (m == 0 || n == 0) {
		return;
	}
	if (alpha == T(0)) {
		det::scale_matrix(m, n, T(0), B, ldb);
		return;
	}

	// 三角行列を 0 詰め dense に展開し (copy のみ)，temp := op(A)*B (または B*op(A))
	// を tgemm で計算して B := alpha*temp で書き戻す．
	// 三角の外は exact な 0 との積和なので結果に影響しない (区間型でも幅は増えない)
	const int nb = lside ? m : n;
	std::vector<T> Adense(static_cast<std::size_t>(nb) * nb, T(0));
	det::expand_triangular(upper, nounit, nb, A, lda, Adense.data());
	std::vector<T> temp(static_cast<std::size_t>(m) * n, T(0));
	if (lside) {
		tgemm(ntrans ? 'N' : 'T', 'N', m, n, m, T(1), Adense.data(), m, B, ldb, T(0), temp.data(), m);
	}
	else {
		tgemm('N', ntrans ? 'N' : 'T', m, n, n, T(1), B, ldb, Adense.data(), n, T(0), temp.data(), m);
	}
	const bool one_alpha = (alpha == T(1));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(static_cast<double>(m) * n))
#endif
	for (int j = 0; j < n; j++) {
		const T* t = temp.data() + static_cast<std::size_t>(m) * j;
		T* b = B + static_cast<std::size_t>(ldb) * j;
		if (one_alpha) {
			for (int i = 0; i < m; i++) {
				b[i] = t[i];
			}
		}
		else {
			for (int i = 0; i < m; i++) {
				b[i] = alpha * t[i];
			}
		}
	}
}

// op(A)*X = alpha*B (side 'L') または X*op(A) = alpha*B (side 'R') を解き，
// B を解 X で上書きする．A は三角行列 (non-unit diag は正則を仮定)，B: m x n
template <typename T>
inline void ttrsm(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const T& alpha, const T* A, const int lda,
	T* B, const int ldb
) {
	namespace det = tblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(transa, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!lside && !det::option_is(side, 'R')) {
		det::tblas_error("ttrsm: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("ttrsm: invalid uplo");
	}
	if (!ntrans && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("ttrsm: invalid transa");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::tblas_error("ttrsm: invalid diag");
	}
	const int nrowa = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, m)) {
		det::tblas_error("ttrsm: invalid argument");
	}
	if (m == 0 || n == 0) {
		return;
	}
	if (alpha == T(0)) {
		det::scale_matrix(m, n, T(0), B, ldb);
		return;
	}
	const bool one_alpha = (alpha == T(1));

	if (lside) {
		// B の各列を独立に三角 solve する (列ごとに独立なので並列化できる)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(static_cast<double>(m) * m * n))
#endif
		for (int j = 0; j < n; j++) {
			T* b = B + static_cast<std::size_t>(ldb) * j;
			if (!one_alpha) {
				for (int i = 0; i < m; i++) {
					b[i] = alpha * b[i];
				}
			}
			if (ntrans) {
				if (upper) {
					for (int jj = m - 1; jj >= 0; jj--) {
						const T* a = A + static_cast<std::size_t>(lda) * jj;
						if (b[jj] != T(0)) {
							if (nounit) {
								b[jj] = b[jj] / a[jj];
							}
							const T temp = b[jj];
							for (int i = 0; i < jj; i++) {
								b[i] = b[i] - temp * a[i];
							}
						}
					}
				}
				else {
					for (int jj = 0; jj < m; jj++) {
						const T* a = A + static_cast<std::size_t>(lda) * jj;
						if (b[jj] != T(0)) {
							if (nounit) {
								b[jj] = b[jj] / a[jj];
							}
							const T temp = b[jj];
							for (int i = jj + 1; i < m; i++) {
								b[i] = b[i] - temp * a[i];
							}
						}
					}
				}
			}
			else {
				if (upper) {
					// A^T は lower: 前進代入 (dot 形式)
					for (int jj = 0; jj < m; jj++) {
						const T* a = A + static_cast<std::size_t>(lda) * jj;
						T temp = b[jj];
						for (int i = 0; i < jj; i++) {
							temp = temp - a[i] * b[i];
						}
						if (nounit) {
							temp = temp / a[jj];
						}
						b[jj] = temp;
					}
				}
				else {
					// A^T は upper: 後退代入 (dot 形式)
					for (int jj = m - 1; jj >= 0; jj--) {
						const T* a = A + static_cast<std::size_t>(lda) * jj;
						T temp = b[jj];
						for (int i = jj + 1; i < m; i++) {
							temp = temp - a[i] * b[i];
						}
						if (nounit) {
							temp = temp / a[jj];
						}
						b[jj] = temp;
					}
				}
			}
		}
	}
	else {
		// 列 sweep solve: X_j = (alpha*B_j - sum_{先行 i} X_i * coef(i,j)) / diag(j)
		// 列の処理順は先行列のみに依存するように選ぶ (列間に依存があるため逐次)
		const bool ascending = (upper == ntrans);
		for (int jj = 0; jj < n; jj++) {
			const int j = ascending ? jj : n - 1 - jj;
			T* bj = B + static_cast<std::size_t>(ldb) * j;
			if (!one_alpha) {
				for (int r = 0; r < m; r++) {
					bj[r] = alpha * bj[r];
				}
			}
			const int i0 = ascending ? 0 : j + 1;
			const int i1 = ascending ? j : n;
			for (int i = i0; i < i1; i++) {
				// coef(i,j): op(A)(i,j)，ntrans なら A(i,j)，trans なら A(j,i)
				const T aij = ntrans ? A[i + static_cast<std::size_t>(lda) * j] : A[j + static_cast<std::size_t>(lda) * i];
				if (aij == T(0)) {
					continue;
				}
				const T* bi = B + static_cast<std::size_t>(ldb) * i;
				for (int r = 0; r < m; r++) {
					bj[r] = bj[r] - aij * bi[r];
				}
			}
			if (nounit) {
				const T ajj = A[j + static_cast<std::size_t>(lda) * j];
				for (int r = 0; r < m; r++) {
					bj[r] = bj[r] / ajj;
				}
			}
		}
	}
}

// C の uplo 三角部分のみ := alpha*op(A)*op(B) + beta*C (C: n x n, op(A): n x k,
// op(B): k x n)．reference BLAS の GEMMTR (LAPACK 3.12.1 で追加) に対応する．
template <typename T>
inline void tgemmtr(
	const char uplo, const char transa, const char transb,
	const int n, const int k,
	const T& alpha, const T* A, const int lda,
	const T* B, const int ldb,
	const T& beta, T* C, const int ldc
) {
	namespace det = tblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tgemmtr: invalid uplo");
	}
	const bool ta = !det::option_is(transa, 'N');
	const bool tb = !det::option_is(transb, 'N');
	if (ta && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("tgemmtr: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::tblas_error("tgemmtr: invalid transb");
	}
	const int nrowa = ta ? k : n;
	const int nrowb = tb ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, n)) {
		det::tblas_error("tgemmtr: invalid argument");
	}
	if (n == 0 || ((alpha == T(0) || k == 0) && beta == T(1))) {
		return;
	}
	if (alpha == T(0) || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc);
		return;
	}
	const bool zero_beta = (beta == T(0));

	// 三角部分のみ dot 形式: C(i,j) = alpha * sum_l opA(i,l)*opB(l,j) + beta*C(i,j)
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(1.0 * n * n * k))
#endif
	for (int j = 0; j < n; j++) {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		const int i0 = upper ? 0 : j;
		const int i1 = upper ? j + 1 : n;
		for (int i = i0; i < i1; i++) {
			T temp = T(0);
			for (int l = 0; l < k; l++) {
				const T ail = ta ? A[l + static_cast<std::size_t>(lda) * i] : A[i + static_cast<std::size_t>(lda) * l];
				const T blj = tb ? B[j + static_cast<std::size_t>(ldb) * l] : B[l + static_cast<std::size_t>(ldb) * j];
				temp = ail * blj + temp;
			}
			if (zero_beta) {
				c[i] = alpha * temp;
			}
			else {
				c[i] = beta * c[i] + alpha * temp;
			}
		}
	}
}

} // namespace vcp

#endif // TBLAS_TBLAS_LEVEL3_HPP
