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
// This file contains code derived from reference LAPACK 3.12.1
// (https://www.netlib.org/lapack/),
// Copyright (c) 1992-2023 The University of Tennessee and The University
//                         of Tennessee Research Foundation,
// Copyright (c) 2000-2023 The University of California Berkeley,
// Copyright (c) 2006-2023 The University of Colorado Denver.
// All rights reserved.
// Distributed under the modified BSD license; see vlapack/LICENSE_LAPACK.txt
// for the full license text, the list of conditions and the disclaimer.
// Modified by Kouta Sekine (2026): translated to C++, BLAS calls replaced by
// tblas, no rounding-mode control, workspace arguments removed.
// ----------------------------------------------------------------------------

#pragma once

#ifndef TLAPACK_TLAPACK_SY_HPP
#define TLAPACK_TLAPACK_SY_HPP

// 対称不定値 (Bunch-Kaufman) と一般化対称固有値問題 (template 版)．
// reference LAPACK 3.12.1 の dsytf2/dlasyf/dsytrf/dsytrs/dsysv/
// dsygs2/dsygst/dsygv を移植．
//
// - tsytf2/tsytrf の ipiv は LAPACK と同じ 1-based (2x2 block は負値) (注意)．
// - tlasyf の GEMMTR (三角部分のみの GEMM, LAPACK 3.12 で導入) は
//   LAPACK 3.11 の dlasyf と同じ blocked 更新 (gemv + gemm) で代替している．
// - tsysv は dsytrs2 ではなく常に tsytrs を使う．

#include "tlapack_geev.hpp"

namespace vcp {

// Bunch-Kaufman 分解 (unblocked)．ipiv は 1-based 符号付き
template <typename T>
inline int tsytf2(const char uplo, const int n, T* A, const int lda, int* ipiv) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsytf2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsytf2: invalid n/lda");
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const T alpha = (T(1) + tlapack_detail::tsqrt(T(17))) / T(8);
	int info = 0;
	if (upper) {
		int k = n;
		while (k >= 1) {
			int kstep = 1;
			int kp;
			const T absakk = tlapack_detail::tabs(A[(k - 1) + L * (k - 1)]);
			int imax = 0;
			T colmax = T(0);
			if (k > 1) {
				imax = itamax(k - 1, A + L * (k - 1), 1) + 1;
				colmax = tlapack_detail::tabs(A[(imax - 1) + L * (k - 1)]);
			}
			if (std::max(absakk, colmax) == T(0) || tlapack_detail::tisnan(absakk)) {
				if (info == 0) {
					info = k;
				}
				kp = k;
			}
			else {
				if (absakk >= alpha * colmax) {
					kp = k;
				}
				else {
					int jmax = imax + itamax(k - imax, A + (imax - 1) + L * imax, lda) + 1;
					T rowmax = tlapack_detail::tabs(A[(imax - 1) + L * (jmax - 1)]);
					if (imax > 1) {
						jmax = itamax(imax - 1, A + L * (imax - 1), 1) + 1;
						rowmax = std::max(rowmax, tlapack_detail::tabs(A[(jmax - 1) + L * (imax - 1)]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (tlapack_detail::tabs(A[(imax - 1) + L * (imax - 1)]) >= alpha * rowmax) {
						kp = imax;
					}
					else {
						kp = imax;
						kstep = 2;
					}
				}
				const int kk = k - kstep + 1;
				if (kp != kk) {
					tswap(kp - 1, A + L * (kk - 1), 1, A + L * (kp - 1), 1);
					tswap(kk - kp - 1, A + kp + L * (kk - 1), 1, A + (kp - 1) + L * kp, lda);
					T t = A[(kk - 1) + L * (kk - 1)];
					A[(kk - 1) + L * (kk - 1)] = A[(kp - 1) + L * (kp - 1)];
					A[(kp - 1) + L * (kp - 1)] = t;
					if (kstep == 2) {
						t = A[(k - 2) + L * (k - 1)];
						A[(k - 2) + L * (k - 1)] = A[(kp - 1) + L * (k - 1)];
						A[(kp - 1) + L * (k - 1)] = t;
					}
				}
				if (kstep == 1) {
					const T r1 = T(1) / A[(k - 1) + L * (k - 1)];
					tsyr('U', k - 1, -r1, A + L * (k - 1), 1, A, lda);
					tscal(k - 1, r1, A + L * (k - 1), 1);
				}
				else {
					if (k > 2) {
						T d12 = A[(k - 2) + L * (k - 1)];
						const T d22 = A[(k - 2) + L * (k - 2)] / d12;
						const T d11 = A[(k - 1) + L * (k - 1)] / d12;
						const T t = T(1) / (d11 * d22 - T(1));
						d12 = t / d12;
						for (int j = k - 2; j >= 1; j--) {
							const T wkm1 = d12 * (d11 * A[(j - 1) + L * (k - 2)] - A[(j - 1) + L * (k - 1)]);
							const T wk = d12 * (d22 * A[(j - 1) + L * (k - 1)] - A[(j - 1) + L * (k - 2)]);
							for (int i = j; i >= 1; i--) {
								A[(i - 1) + L * (j - 1)] = A[(i - 1) + L * (j - 1)] -
									A[(i - 1) + L * (k - 1)] * wk - A[(i - 1) + L * (k - 2)] * wkm1;
							}
							A[(j - 1) + L * (k - 1)] = wk;
							A[(j - 1) + L * (k - 2)] = wkm1;
						}
					}
				}
			}
			if (kstep == 1) {
				ipiv[k - 1] = kp;
			}
			else {
				ipiv[k - 1] = -kp;
				ipiv[k - 2] = -kp;
			}
			k -= kstep;
		}
	}
	else {
		int k = 1;
		while (k <= n) {
			int kstep = 1;
			int kp;
			const T absakk = tlapack_detail::tabs(A[(k - 1) + L * (k - 1)]);
			int imax = 0;
			T colmax = T(0);
			if (k < n) {
				imax = k + itamax(n - k, A + k + L * (k - 1), 1) + 1;
				colmax = tlapack_detail::tabs(A[(imax - 1) + L * (k - 1)]);
			}
			if (std::max(absakk, colmax) == T(0) || tlapack_detail::tisnan(absakk)) {
				if (info == 0) {
					info = k;
				}
				kp = k;
			}
			else {
				if (absakk >= alpha * colmax) {
					kp = k;
				}
				else {
					int jmax = k - 1 + itamax(imax - k, A + (imax - 1) + L * (k - 1), lda) + 1;
					T rowmax = tlapack_detail::tabs(A[(imax - 1) + L * (jmax - 1)]);
					if (imax < n) {
						jmax = imax + itamax(n - imax, A + imax + L * (imax - 1), 1) + 1;
						rowmax = std::max(rowmax, tlapack_detail::tabs(A[(jmax - 1) + L * (imax - 1)]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (tlapack_detail::tabs(A[(imax - 1) + L * (imax - 1)]) >= alpha * rowmax) {
						kp = imax;
					}
					else {
						kp = imax;
						kstep = 2;
					}
				}
				const int kk = k + kstep - 1;
				if (kp != kk) {
					if (kp < n) {
						tswap(n - kp, A + kp + L * (kk - 1), 1, A + kp + L * (kp - 1), 1);
					}
					tswap(kp - kk - 1, A + kk + L * (kk - 1), 1, A + (kp - 1) + L * kk, lda);
					T t = A[(kk - 1) + L * (kk - 1)];
					A[(kk - 1) + L * (kk - 1)] = A[(kp - 1) + L * (kp - 1)];
					A[(kp - 1) + L * (kp - 1)] = t;
					if (kstep == 2) {
						t = A[k + L * (k - 1)];
						A[k + L * (k - 1)] = A[(kp - 1) + L * (k - 1)];
						A[(kp - 1) + L * (k - 1)] = t;
					}
				}
				if (kstep == 1) {
					if (k < n) {
						const T d11 = T(1) / A[(k - 1) + L * (k - 1)];
						tsyr('L', n - k, -d11, A + k + L * (k - 1), 1, A + k + L * k, lda);
						tscal(n - k, d11, A + k + L * (k - 1), 1);
					}
				}
				else {
					if (k < n - 1) {
						T d21 = A[k + L * (k - 1)];
						const T d11 = A[k + L * k] / d21;
						const T d22 = A[(k - 1) + L * (k - 1)] / d21;
						const T t = T(1) / (d11 * d22 - T(1));
						d21 = t / d21;
						for (int j = k + 2; j <= n; j++) {
							const T wk = d21 * (d11 * A[(j - 1) + L * (k - 1)] - A[(j - 1) + L * k]);
							const T wkp1 = d21 * (d22 * A[(j - 1) + L * k] - A[(j - 1) + L * (k - 1)]);
							for (int i = j; i <= n; i++) {
								A[(i - 1) + L * (j - 1)] = A[(i - 1) + L * (j - 1)] -
									A[(i - 1) + L * (k - 1)] * wk - A[(i - 1) + L * k] * wkp1;
							}
							A[(j - 1) + L * (k - 1)] = wk;
							A[(j - 1) + L * k] = wkp1;
						}
					}
				}
			}
			if (kstep == 1) {
				ipiv[k - 1] = kp;
			}
			else {
				ipiv[k - 1] = -kp;
				ipiv[k] = -kp;
			}
			k += kstep;
		}
	}
	return info;
}

// Bunch-Kaufman 分解の panel (blocked 用)．kb に処理済み列数を返す
template <typename T>
inline int tlasyf(const char uplo, const int n, const int nb, int& kb, T* A, const int lda,
	int* ipiv, T* W, const int ldw) {
	namespace detail = tlapack_detail;
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LW = static_cast<std::size_t>(ldw);
	const T alpha = (T(1) + tlapack_detail::tsqrt(T(17))) / T(8);
	int info = 0;
	if (detail::option_is(uplo, 'U')) {
		int k = n;
		int kw = 0;
		for (;;) {
			kw = nb + k - n;
			if ((k <= n - nb + 1 && nb < n) || k < 1) {
				break;
			}
			tcopy(k, A + L * (k - 1), 1, W + LW * (kw - 1), 1);
			if (k < n) {
				tgemv('N', k, n - k, -T(1), A + L * k, lda, W + (k - 1) + LW * kw, ldw,
					T(1), W + LW * (kw - 1), 1);
			}
			int kstep = 1;
			int kp;
			const T absakk = tlapack_detail::tabs(W[(k - 1) + LW * (kw - 1)]);
			int imax = 0;
			T colmax = T(0);
			if (k > 1) {
				imax = itamax(k - 1, W + LW * (kw - 1), 1) + 1;
				colmax = tlapack_detail::tabs(W[(imax - 1) + LW * (kw - 1)]);
			}
			if (std::max(absakk, colmax) == T(0)) {
				if (info == 0) {
					info = k;
				}
				kp = k;
			}
			else {
				if (absakk >= alpha * colmax) {
					kp = k;
				}
				else {
					tcopy(imax, A + L * (imax - 1), 1, W + LW * (kw - 2), 1);
					tcopy(k - imax, A + (imax - 1) + L * imax, lda, W + imax + LW * (kw - 2), 1);
					if (k < n) {
						tgemv('N', k, n - k, -T(1), A + L * k, lda, W + (imax - 1) + LW * kw, ldw,
							T(1), W + LW * (kw - 2), 1);
					}
					int jmax = imax + itamax(k - imax, W + imax + LW * (kw - 2), 1) + 1;
					T rowmax = tlapack_detail::tabs(W[(jmax - 1) + LW * (kw - 2)]);
					if (imax > 1) {
						jmax = itamax(imax - 1, W + LW * (kw - 2), 1) + 1;
						rowmax = std::max(rowmax, tlapack_detail::tabs(W[(jmax - 1) + LW * (kw - 2)]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (tlapack_detail::tabs(W[(imax - 1) + LW * (kw - 2)]) >= alpha * rowmax) {
						kp = imax;
						tcopy(k, W + LW * (kw - 2), 1, W + LW * (kw - 1), 1);
					}
					else {
						kp = imax;
						kstep = 2;
					}
				}
				const int kk = k - kstep + 1;
				const int kkw = nb + kk - n;
				if (kp != kk) {
					A[(kp - 1) + L * (kp - 1)] = A[(kk - 1) + L * (kk - 1)];
					tcopy(kk - 1 - kp, A + kp + L * (kk - 1), 1, A + (kp - 1) + L * kp, lda);
					if (kp > 1) {
						tcopy(kp - 1, A + L * (kk - 1), 1, A + L * (kp - 1), 1);
					}
					if (k < n) {
						tswap(n - k, A + (kk - 1) + L * k, lda, A + (kp - 1) + L * k, lda);
					}
					tswap(n - kk + 1, W + (kk - 1) + LW * (kkw - 1), ldw, W + (kp - 1) + LW * (kkw - 1), ldw);
				}
				if (kstep == 1) {
					tcopy(k, W + LW * (kw - 1), 1, A + L * (k - 1), 1);
					const T r1 = T(1) / A[(k - 1) + L * (k - 1)];
					tscal(k - 1, r1, A + L * (k - 1), 1);
				}
				else {
					if (k > 2) {
						T d21 = W[(k - 2) + LW * (kw - 1)];
						const T d11 = W[(k - 1) + LW * (kw - 1)] / d21;
						const T d22 = W[(k - 2) + LW * (kw - 2)] / d21;
						const T t = T(1) / (d11 * d22 - T(1));
						d21 = t / d21;
						for (int j = 1; j <= k - 2; j++) {
							A[(j - 1) + L * (k - 2)] = d21 * (d11 * W[(j - 1) + LW * (kw - 2)] -
								W[(j - 1) + LW * (kw - 1)]);
							A[(j - 1) + L * (k - 1)] = d21 * (d22 * W[(j - 1) + LW * (kw - 1)] -
								W[(j - 1) + LW * (kw - 2)]);
						}
					}
					A[(k - 2) + L * (k - 2)] = W[(k - 2) + LW * (kw - 2)];
					A[(k - 2) + L * (k - 1)] = W[(k - 2) + LW * (kw - 1)];
					A[(k - 1) + L * (k - 1)] = W[(k - 1) + LW * (kw - 1)];
				}
			}
			if (kstep == 1) {
				ipiv[k - 1] = kp;
			}
			else {
				ipiv[k - 1] = -kp;
				ipiv[k - 2] = -kp;
			}
			k -= kstep;
		}
		// A(1:k, 1:k) の上三角を A11 := A11 - U12*W^T で更新 (GEMMTR 相当)
		for (int j = ((k - 1) / nb) * nb + 1; j >= 1; j -= nb) {
			const int jb = std::min(nb, k - j + 1);
			for (int jj = j; jj <= j + jb - 1; jj++) {
				tgemv('N', jj - j + 1, n - k, -T(1), A + (j - 1) + L * k, lda,
					W + (jj - 1) + LW * kw, ldw, T(1), A + (j - 1) + L * (jj - 1), 1);
			}
			if (j >= 2) {
				tgemm('N', 'T', j - 1, jb, n - k, -T(1), A + L * k, lda,
					W + (j - 1) + LW * kw, ldw, T(1), A + L * (j - 1), lda);
			}
		}
		// 部分的に行った交換を戻す
		{
			int j = k + 1;
			do {
				const int jj = j;
				int jp = ipiv[j - 1];
				if (jp < 0) {
					jp = -jp;
					j = j + 1;
				}
				j = j + 1;
				if (jp != jj && j <= n) {
					tswap(n - j + 1, A + (jp - 1) + L * (j - 1), lda, A + (jj - 1) + L * (j - 1), lda);
				}
			} while (j < n);
		}
		kb = n - k;
	}
	else {
		int k = 1;
		for (;;) {
			if ((k >= nb && nb < n) || k > n) {
				break;
			}
			tcopy(n - k + 1, A + (k - 1) + L * (k - 1), 1, W + (k - 1) + LW * (k - 1), 1);
			tgemv('N', n - k + 1, k - 1, -T(1), A + (k - 1), lda, W + (k - 1), ldw,
				T(1), W + (k - 1) + LW * (k - 1), 1);
			int kstep = 1;
			int kp;
			const T absakk = tlapack_detail::tabs(W[(k - 1) + LW * (k - 1)]);
			int imax = 0;
			T colmax = T(0);
			if (k < n) {
				imax = k + itamax(n - k, W + k + LW * (k - 1), 1) + 1;
				colmax = tlapack_detail::tabs(W[(imax - 1) + LW * (k - 1)]);
			}
			if (std::max(absakk, colmax) == T(0)) {
				if (info == 0) {
					info = k;
				}
				kp = k;
			}
			else {
				if (absakk >= alpha * colmax) {
					kp = k;
				}
				else {
					tcopy(imax - k, A + (imax - 1) + L * (k - 1), lda, W + (k - 1) + LW * k, 1);
					tcopy(n - imax + 1, A + (imax - 1) + L * (imax - 1), 1, W + (imax - 1) + LW * k, 1);
					tgemv('N', n - k + 1, k - 1, -T(1), A + (k - 1), lda, W + (imax - 1), ldw,
						T(1), W + (k - 1) + LW * k, 1);
					int jmax = k - 1 + itamax(imax - k, W + (k - 1) + LW * k, 1) + 1;
					T rowmax = tlapack_detail::tabs(W[(jmax - 1) + LW * k]);
					if (imax < n) {
						jmax = imax + itamax(n - imax, W + imax + LW * k, 1) + 1;
						rowmax = std::max(rowmax, tlapack_detail::tabs(W[(jmax - 1) + LW * k]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (tlapack_detail::tabs(W[(imax - 1) + LW * k]) >= alpha * rowmax) {
						kp = imax;
						tcopy(n - k + 1, W + (k - 1) + LW * k, 1, W + (k - 1) + LW * (k - 1), 1);
					}
					else {
						kp = imax;
						kstep = 2;
					}
				}
				const int kk = k + kstep - 1;
				if (kp != kk) {
					A[(kp - 1) + L * (kp - 1)] = A[(kk - 1) + L * (kk - 1)];
					tcopy(kp - kk - 1, A + kk + L * (kk - 1), 1, A + (kp - 1) + L * kk, lda);
					if (kp < n) {
						tcopy(n - kp, A + kp + L * (kk - 1), 1, A + kp + L * (kp - 1), 1);
					}
					if (k > 1) {
						tswap(k - 1, A + (kk - 1), lda, A + (kp - 1), lda);
					}
					tswap(kk, W + (kk - 1), ldw, W + (kp - 1), ldw);
				}
				if (kstep == 1) {
					tcopy(n - k + 1, W + (k - 1) + LW * (k - 1), 1, A + (k - 1) + L * (k - 1), 1);
					if (k < n) {
						const T r1 = T(1) / A[(k - 1) + L * (k - 1)];
						tscal(n - k, r1, A + k + L * (k - 1), 1);
					}
				}
				else {
					if (k < n - 1) {
						T d21 = W[k + LW * (k - 1)];
						const T d11 = W[k + LW * k] / d21;
						const T d22 = W[(k - 1) + LW * (k - 1)] / d21;
						const T t = T(1) / (d11 * d22 - T(1));
						d21 = t / d21;
						for (int j = k + 2; j <= n; j++) {
							A[(j - 1) + L * (k - 1)] = d21 * (d11 * W[(j - 1) + LW * (k - 1)] -
								W[(j - 1) + LW * k]);
							A[(j - 1) + L * k] = d21 * (d22 * W[(j - 1) + LW * k] -
								W[(j - 1) + LW * (k - 1)]);
						}
					}
					A[(k - 1) + L * (k - 1)] = W[(k - 1) + LW * (k - 1)];
					A[k + L * (k - 1)] = W[k + LW * (k - 1)];
					A[k + L * k] = W[k + LW * k];
				}
			}
			if (kstep == 1) {
				ipiv[k - 1] = kp;
			}
			else {
				ipiv[k - 1] = -kp;
				ipiv[k] = -kp;
			}
			k += kstep;
		}
		// A(k:n, k:n) の下三角を A22 := A22 - L21*W^T で更新 (GEMMTR 相当)
		for (int j = k; j <= n; j += nb) {
			const int jb = std::min(nb, n - j + 1);
			for (int jj = j; jj <= j + jb - 1; jj++) {
				tgemv('N', j + jb - jj, k - 1, -T(1), A + (jj - 1), lda, W + (jj - 1), ldw,
					T(1), A + (jj - 1) + L * (jj - 1), 1);
			}
			if (j + jb <= n) {
				tgemm('N', 'T', n - j - jb + 1, jb, k - 1, -T(1), A + (j + jb - 1), lda,
					W + (j - 1), ldw, T(1), A + (j + jb - 1) + L * (j - 1), lda);
			}
		}
		// 部分的に行った交換を戻す
		{
			int j = k - 1;
			while (j >= 1) {
				const int jj = j;
				int jp = ipiv[j - 1];
				if (jp < 0) {
					jp = -jp;
					j = j - 1;
				}
				j = j - 1;
				if (jp != jj && j >= 1) {
					tswap(j, A + (jp - 1), lda, A + (jj - 1), lda);
				}
				if (j <= 1) {
					break;
				}
			}
		}
		kb = k - 1;
	}
	return info;
}

// Bunch-Kaufman 分解 (blocked)．ipiv は 1-based 符号付き
template <typename T>
inline int tsytrf(const char uplo, const int n, T* A, const int lda, int* ipiv) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsytrf: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsytrf: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	int nb = detail::ilaenv(1, "SYTRF", n, -1, -1, -1);
	const int nbmin = detail::ilaenv(2, "SYTRF", n, -1, -1, -1);
	if (nb < nbmin || nb >= n) {
		nb = n;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int ldwork = n;
	std::vector<T> work(nb < n ? static_cast<std::size_t>(ldwork) * nb : 1, T(0));
	int info = 0;
	if (upper) {
		int k = n;
		while (k >= 1) {
			int kb;
			int iinfo;
			if (k > nb) {
				iinfo = tlasyf(uplo, k, nb, kb, A, lda, ipiv, work.data(), ldwork);
			}
			else {
				iinfo = tsytf2(uplo, k, A, lda, ipiv);
				kb = k;
			}
			if (info == 0 && iinfo > 0) {
				info = iinfo;
			}
			k -= kb;
		}
	}
	else {
		int k = 1;
		while (k <= n) {
			int kb;
			int iinfo;
			if (k <= n - nb) {
				iinfo = tlasyf(uplo, n - k + 1, nb, kb, A + (k - 1) + L * (k - 1), lda,
					ipiv + (k - 1), work.data(), ldwork);
			}
			else {
				iinfo = tsytf2(uplo, n - k + 1, A + (k - 1) + L * (k - 1), lda, ipiv + (k - 1));
				kb = n - k + 1;
			}
			if (info == 0 && iinfo > 0) {
				info = iinfo + k - 1;
			}
			for (int j = k; j <= k + kb - 1; j++) {
				if (ipiv[j - 1] > 0) {
					ipiv[j - 1] = ipiv[j - 1] + k - 1;
				}
				else {
					ipiv[j - 1] = ipiv[j - 1] - k + 1;
				}
			}
			k += kb;
		}
	}
	return info;
}

// Bunch-Kaufman 分解を使った連立一次方程式の求解
template <typename T>
inline int tsytrs(const char uplo, const int n, const int nrhs, const T* A, const int lda,
	const int* ipiv, T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsytrs: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsytrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LB = static_cast<std::size_t>(ldb);
	if (upper) {
		// B := U^-Tf D^-1 U^-1 B (後退方向 → 前進方向)
		int k = n;
		while (k >= 1) {
			if (ipiv[k - 1] > 0) {
				const int kp = ipiv[k - 1];
				if (kp != k) {
					tswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				tger(k - 1, nrhs, -T(1), A + L * (k - 1), 1, B + (k - 1), ldb, B, ldb);
				tscal(nrhs, T(1) / A[(k - 1) + L * (k - 1)], B + (k - 1), ldb);
				k = k - 1;
			}
			else {
				const int kp = -ipiv[k - 1];
				if (kp != k - 1) {
					tswap(nrhs, B + (k - 2), ldb, B + (kp - 1), ldb);
				}
				tger(k - 2, nrhs, -T(1), A + L * (k - 1), 1, B + (k - 1), ldb, B, ldb);
				tger(k - 2, nrhs, -T(1), A + L * (k - 2), 1, B + (k - 2), ldb, B, ldb);
				const T akm1k = A[(k - 2) + L * (k - 1)];
				const T akm1 = A[(k - 2) + L * (k - 2)] / akm1k;
				const T ak = A[(k - 1) + L * (k - 1)] / akm1k;
				const T denom = akm1 * ak - T(1);
				for (int j = 1; j <= nrhs; j++) {
					const T bkm1 = B[(k - 2) + LB * (j - 1)] / akm1k;
					const T bk = B[(k - 1) + LB * (j - 1)] / akm1k;
					B[(k - 2) + LB * (j - 1)] = (ak * bkm1 - bk) / denom;
					B[(k - 1) + LB * (j - 1)] = (akm1 * bk - bkm1) / denom;
				}
				k = k - 2;
			}
		}
		k = 1;
		while (k <= n) {
			if (ipiv[k - 1] > 0) {
				tgemv('T', k - 1, nrhs, -T(1), B, ldb, A + L * (k - 1), 1, T(1), B + (k - 1), ldb);
				const int kp = ipiv[k - 1];
				if (kp != k) {
					tswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k + 1;
			}
			else {
				tgemv('T', k - 1, nrhs, -T(1), B, ldb, A + L * (k - 1), 1, T(1), B + (k - 1), ldb);
				tgemv('T', k - 1, nrhs, -T(1), B, ldb, A + L * k, 1, T(1), B + k, ldb);
				const int kp = -ipiv[k - 1];
				if (kp != k) {
					tswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k + 2;
			}
		}
	}
	else {
		int k = 1;
		while (k <= n) {
			if (ipiv[k - 1] > 0) {
				const int kp = ipiv[k - 1];
				if (kp != k) {
					tswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				if (k < n) {
					tger(n - k, nrhs, -T(1), A + k + L * (k - 1), 1, B + (k - 1), ldb, B + k, ldb);
				}
				tscal(nrhs, T(1) / A[(k - 1) + L * (k - 1)], B + (k - 1), ldb);
				k = k + 1;
			}
			else {
				const int kp = -ipiv[k - 1];
				if (kp != k + 1) {
					tswap(nrhs, B + k, ldb, B + (kp - 1), ldb);
				}
				if (k < n - 1) {
					tger(n - k - 1, nrhs, -T(1), A + (k + 1) + L * (k - 1), 1, B + (k - 1), ldb,
						B + (k + 1), ldb);
					tger(n - k - 1, nrhs, -T(1), A + (k + 1) + L * k, 1, B + k, ldb,
						B + (k + 1), ldb);
				}
				const T akm1k = A[k + L * (k - 1)];
				const T akm1 = A[(k - 1) + L * (k - 1)] / akm1k;
				const T ak = A[k + L * k] / akm1k;
				const T denom = akm1 * ak - T(1);
				for (int j = 1; j <= nrhs; j++) {
					const T bkm1 = B[(k - 1) + LB * (j - 1)] / akm1k;
					const T bk = B[k + LB * (j - 1)] / akm1k;
					B[(k - 1) + LB * (j - 1)] = (ak * bkm1 - bk) / denom;
					B[k + LB * (j - 1)] = (akm1 * bk - bkm1) / denom;
				}
				k = k + 2;
			}
		}
		k = n;
		while (k >= 1) {
			if (ipiv[k - 1] > 0) {
				if (k < n) {
					tgemv('T', n - k, nrhs, -T(1), B + k, ldb, A + k + L * (k - 1), 1,
						T(1), B + (k - 1), ldb);
				}
				const int kp = ipiv[k - 1];
				if (kp != k) {
					tswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k - 1;
			}
			else {
				if (k < n) {
					tgemv('T', n - k, nrhs, -T(1), B + k, ldb, A + k + L * (k - 1), 1,
						T(1), B + (k - 1), ldb);
					tgemv('T', n - k, nrhs, -T(1), B + k, ldb, A + k + L * (k - 2), 1,
						T(1), B + (k - 2), ldb);
				}
				const int kp = -ipiv[k - 1];
				if (kp != k) {
					tswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k - 2;
			}
		}
	}
	return 0;
}

// 対称不定値の連立一次方程式 A*X = B の求解 (driver)
template <typename T>
inline int tsysv(const char uplo, const int n, const int nrhs, T* A, const int lda,
	int* ipiv, T* B, const int ldb) {
	namespace detail = tlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsysv: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsysv: invalid argument");
	}
	const int info = tsytrf(uplo, n, A, lda, ipiv);
	if (info == 0) {
		tsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb);
	}
	return info;
}

// 一般化対称固有値問題の標準形への変換 (unblocked)
// itype = 1: A := inv(U^T)*A*inv(U) など, 2/3: A := U*A*U^T など (B は tpotrf 済み)
template <typename T>
inline int tsygs2(const int itype, const char uplo, const int n, T* A, const int lda,
	const T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::tlapack_error("tsygs2: invalid itype");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsygs2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsygs2: invalid argument");
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LBB = static_cast<std::size_t>(ldb);
	if (itype == 1) {
		if (upper) {
			for (int k = 1; k <= n; k++) {
				T akk = A[(k - 1) + L * (k - 1)];
				const T bkk = B[(k - 1) + LBB * (k - 1)];
				akk = akk / (bkk * bkk);
				A[(k - 1) + L * (k - 1)] = akk;
				if (k < n) {
					tscal(n - k, T(1) / bkk, A + (k - 1) + L * k, lda);
					const T ct = -(T(1) / T(2)) * akk;
					taxpy(n - k, ct, B + (k - 1) + LBB * k, ldb, A + (k - 1) + L * k, lda);
					tsyr2('U', n - k, -T(1), A + (k - 1) + L * k, lda, B + (k - 1) + LBB * k, ldb,
						A + k + L * k, lda);
					taxpy(n - k, ct, B + (k - 1) + LBB * k, ldb, A + (k - 1) + L * k, lda);
					ttrsv('U', 'T', 'N', n - k, B + k + LBB * k, ldb, A + (k - 1) + L * k, lda);
				}
			}
		}
		else {
			for (int k = 1; k <= n; k++) {
				T akk = A[(k - 1) + L * (k - 1)];
				const T bkk = B[(k - 1) + LBB * (k - 1)];
				akk = akk / (bkk * bkk);
				A[(k - 1) + L * (k - 1)] = akk;
				if (k < n) {
					tscal(n - k, T(1) / bkk, A + k + L * (k - 1), 1);
					const T ct = -(T(1) / T(2)) * akk;
					taxpy(n - k, ct, B + k + LBB * (k - 1), 1, A + k + L * (k - 1), 1);
					tsyr2('L', n - k, -T(1), A + k + L * (k - 1), 1, B + k + LBB * (k - 1), 1,
						A + k + L * k, lda);
					taxpy(n - k, ct, B + k + LBB * (k - 1), 1, A + k + L * (k - 1), 1);
					ttrsv('L', 'N', 'N', n - k, B + k + LBB * k, ldb, A + k + L * (k - 1), 1);
				}
			}
		}
	}
	else {
		if (upper) {
			for (int k = 1; k <= n; k++) {
				const T akk = A[(k - 1) + L * (k - 1)];
				const T bkk = B[(k - 1) + LBB * (k - 1)];
				ttrmv('U', 'N', 'N', k - 1, B, ldb, A + L * (k - 1), 1);
				const T ct = (T(1) / T(2)) * akk;
				taxpy(k - 1, ct, B + LBB * (k - 1), 1, A + L * (k - 1), 1);
				tsyr2('U', k - 1, T(1), A + L * (k - 1), 1, B + LBB * (k - 1), 1, A, lda);
				taxpy(k - 1, ct, B + LBB * (k - 1), 1, A + L * (k - 1), 1);
				tscal(k - 1, bkk, A + L * (k - 1), 1);
				A[(k - 1) + L * (k - 1)] = akk * (bkk * bkk);
			}
		}
		else {
			for (int k = 1; k <= n; k++) {
				const T akk = A[(k - 1) + L * (k - 1)];
				const T bkk = B[(k - 1) + LBB * (k - 1)];
				ttrmv('L', 'T', 'N', k - 1, B, ldb, A + (k - 1), lda);
				const T ct = (T(1) / T(2)) * akk;
				taxpy(k - 1, ct, B + (k - 1), ldb, A + (k - 1), lda);
				tsyr2('L', k - 1, T(1), A + (k - 1), lda, B + (k - 1), ldb, A, lda);
				taxpy(k - 1, ct, B + (k - 1), ldb, A + (k - 1), lda);
				tscal(k - 1, bkk, A + (k - 1), lda);
				A[(k - 1) + L * (k - 1)] = akk * (bkk * bkk);
			}
		}
	}
	return 0;
}

// 一般化対称固有値問題の標準形への変換 (blocked)
template <typename T>
inline int tsygst(const int itype, const char uplo, const int n, T* A, const int lda,
	const T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::tlapack_error("tsygst: invalid itype");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsygst: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsygst: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LBB = static_cast<std::size_t>(ldb);
	const int nb = detail::ilaenv(1, "SYGST", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		return tsygs2(itype, uplo, n, A, lda, B, ldb);
	}
	if (itype == 1) {
		if (upper) {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				tsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb);
				if (k + kb <= n) {
					ttrsm('L', 'U', 'T', 'N', kb, n - k - kb + 1, T(1), B + (k - 1) + LBB * (k - 1), ldb,
						A + (k - 1) + L * (k + kb - 1), lda);
					tsymm('L', 'U', kb, n - k - kb + 1, -(T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
						B + (k - 1) + LBB * (k + kb - 1), ldb, T(1), A + (k - 1) + L * (k + kb - 1), lda);
					tsyr2k('U', 'T', n - k - kb + 1, kb, -T(1), A + (k - 1) + L * (k + kb - 1), lda,
						B + (k - 1) + LBB * (k + kb - 1), ldb, T(1), A + (k + kb - 1) + L * (k + kb - 1), lda);
					tsymm('L', 'U', kb, n - k - kb + 1, -(T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
						B + (k - 1) + LBB * (k + kb - 1), ldb, T(1), A + (k - 1) + L * (k + kb - 1), lda);
					ttrsm('R', 'U', 'N', 'N', kb, n - k - kb + 1, T(1),
						B + (k + kb - 1) + LBB * (k + kb - 1), ldb, A + (k - 1) + L * (k + kb - 1), lda);
				}
			}
		}
		else {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				tsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb);
				if (k + kb <= n) {
					ttrsm('R', 'L', 'T', 'N', n - k - kb + 1, kb, T(1), B + (k - 1) + LBB * (k - 1), ldb,
						A + (k + kb - 1) + L * (k - 1), lda);
					tsymm('R', 'L', n - k - kb + 1, kb, -(T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
						B + (k + kb - 1) + LBB * (k - 1), ldb, T(1), A + (k + kb - 1) + L * (k - 1), lda);
					tsyr2k('L', 'N', n - k - kb + 1, kb, -T(1), A + (k + kb - 1) + L * (k - 1), lda,
						B + (k + kb - 1) + LBB * (k - 1), ldb, T(1), A + (k + kb - 1) + L * (k + kb - 1), lda);
					tsymm('R', 'L', n - k - kb + 1, kb, -(T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
						B + (k + kb - 1) + LBB * (k - 1), ldb, T(1), A + (k + kb - 1) + L * (k - 1), lda);
					ttrsm('L', 'L', 'N', 'N', n - k - kb + 1, kb, T(1),
						B + (k + kb - 1) + LBB * (k + kb - 1), ldb, A + (k + kb - 1) + L * (k - 1), lda);
				}
			}
		}
	}
	else {
		if (upper) {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				ttrmm('L', 'U', 'N', 'N', k - 1, kb, T(1), B, ldb, A + L * (k - 1), lda);
				tsymm('R', 'U', k - 1, kb, (T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
					B + LBB * (k - 1), ldb, T(1), A + L * (k - 1), lda);
				tsyr2k('U', 'N', k - 1, kb, T(1), A + L * (k - 1), lda, B + LBB * (k - 1), ldb,
					T(1), A, lda);
				tsymm('R', 'U', k - 1, kb, (T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
					B + LBB * (k - 1), ldb, T(1), A + L * (k - 1), lda);
				ttrmm('R', 'U', 'T', 'N', k - 1, kb, T(1), B + (k - 1) + LBB * (k - 1), ldb,
					A + L * (k - 1), lda);
				tsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb);
			}
		}
		else {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				ttrmm('R', 'L', 'N', 'N', kb, k - 1, T(1), B, ldb, A + (k - 1), lda);
				tsymm('L', 'L', kb, k - 1, (T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
					B + (k - 1), ldb, T(1), A + (k - 1), lda);
				tsyr2k('L', 'T', k - 1, kb, T(1), A + (k - 1), lda, B + (k - 1), ldb, T(1), A, lda);
				tsymm('L', 'L', kb, k - 1, (T(1) / T(2)), A + (k - 1) + L * (k - 1), lda,
					B + (k - 1), ldb, T(1), A + (k - 1), lda);
				ttrmm('L', 'L', 'T', 'N', kb, k - 1, T(1), B + (k - 1) + LBB * (k - 1), ldb,
					A + (k - 1), lda);
				tsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb);
			}
		}
	}
	return 0;
}

// 一般化対称固有値問題 A*x = lambda*B*x (itype=1), A*B*x = lambda*x (2),
// B*A*x = lambda*x (3) の driver (B は対称正定値，上書きされる)
template <typename T>
inline int tsygv(const int itype, const char jobz, const char uplo, const int n,
	T* A, const int lda, T* B, const int ldb, T* w) {
	namespace detail = tlapack_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::tlapack_error("tsygv: invalid itype");
	}
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::tlapack_error("tsygv: invalid jobz");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsygv: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsygv: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	int info = tpotrf(uplo, n, B, ldb);
	if (info != 0) {
		return n + info;
	}
	tsygst(itype, uplo, n, A, lda, B, ldb);
	info = tsyev(jobz, uplo, n, A, lda, w);
	if (wantz) {
		int neig = n;
		if (info > 0) {
			neig = info - 1;
		}
		if (itype == 1 || itype == 2) {
			const char trans = upper ? 'N' : 'T';
			ttrsm('L', uplo, trans, 'N', n, neig, T(1), B, ldb, A, lda);
		}
		else {
			const char trans = upper ? 'T' : 'N';
			ttrmm('L', uplo, trans, 'N', n, neig, T(1), B, ldb, A, lda);
		}
	}
	return info;
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_SY_HPP
