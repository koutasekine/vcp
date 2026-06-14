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
// rdblas, rounding-mode argument added, workspace arguments removed.
// ----------------------------------------------------------------------------

#pragma once

#ifndef VBLAS_RDLAPACK_SY_HPP
#define VBLAS_RDLAPACK_SY_HPP

// 対称不定値 (Bunch-Kaufman) と一般化対称固有値問題 (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dsytf2/dlasyf/dsytrf/dsytrs/dsysv/
// dsygs2/dsygst/dsygv を移植．
//
// - rdsytf2/rdsytrf の ipiv は LAPACK と同じ 1-based (2x2 block は負値) (注意)．
// - rdlasyf の GEMMTR (三角部分のみの GEMM, LAPACK 3.12 で導入) は
//   LAPACK 3.11 の dlasyf と同じ blocked 更新 (gemv + gemm) で代替している．
// - rdsysv は dsytrs2 ではなく常に rdsytrs を使う．

#include "rdlapack_geev.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// Bunch-Kaufman 分解 (unblocked)．ipiv は 1-based 符号付き
inline int rdsytf2(const char uplo, const int n, double* A, const int lda, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsytf2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdsytf2: invalid n/lda");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	const double alpha = (1.0 + std::sqrt(17.0)) / 8.0;
	int info = 0;
	if (upper) {
		int k = n;
		while (k >= 1) {
			int kstep = 1;
			int kp;
			const double absakk = std::fabs(A[(k - 1) + L * (k - 1)]);
			int imax = 0;
			double colmax = 0.0;
			if (k > 1) {
				imax = idamax(k - 1, A + L * (k - 1), 1) + 1;
				colmax = std::fabs(A[(imax - 1) + L * (k - 1)]);
			}
			if (std::max(absakk, colmax) == 0.0 || std::isnan(absakk)) {
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
					int jmax = imax + idamax(k - imax, A + (imax - 1) + L * imax, lda) + 1;
					double rowmax = std::fabs(A[(imax - 1) + L * (jmax - 1)]);
					if (imax > 1) {
						jmax = idamax(imax - 1, A + L * (imax - 1), 1) + 1;
						rowmax = std::max(rowmax, std::fabs(A[(jmax - 1) + L * (imax - 1)]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (std::fabs(A[(imax - 1) + L * (imax - 1)]) >= alpha * rowmax) {
						kp = imax;
					}
					else {
						kp = imax;
						kstep = 2;
					}
				}
				const int kk = k - kstep + 1;
				if (kp != kk) {
					dswap(kp - 1, A + L * (kk - 1), 1, A + L * (kp - 1), 1);
					dswap(kk - kp - 1, A + kp + L * (kk - 1), 1, A + (kp - 1) + L * kp, lda);
					double t = A[(kk - 1) + L * (kk - 1)];
					A[(kk - 1) + L * (kk - 1)] = A[(kp - 1) + L * (kp - 1)];
					A[(kp - 1) + L * (kp - 1)] = t;
					if (kstep == 2) {
						t = A[(k - 2) + L * (k - 1)];
						A[(k - 2) + L * (k - 1)] = A[(kp - 1) + L * (k - 1)];
						A[(kp - 1) + L * (k - 1)] = t;
					}
				}
				if (kstep == 1) {
					const double r1 = 1.0 / A[(k - 1) + L * (k - 1)];
					rdsyr('U', k - 1, -r1, A + L * (k - 1), 1, A, lda, rm);
					rdscal(k - 1, r1, A + L * (k - 1), 1, rm);
				}
				else {
					if (k > 2) {
						double d12 = A[(k - 2) + L * (k - 1)];
						const double d22 = A[(k - 2) + L * (k - 2)] / d12;
						const double d11 = A[(k - 1) + L * (k - 1)] / d12;
						const double t = 1.0 / (d11 * d22 - 1.0);
						d12 = t / d12;
						for (int j = k - 2; j >= 1; j--) {
							const double wkm1 = d12 * (d11 * A[(j - 1) + L * (k - 2)] - A[(j - 1) + L * (k - 1)]);
							const double wk = d12 * (d22 * A[(j - 1) + L * (k - 1)] - A[(j - 1) + L * (k - 2)]);
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
			const double absakk = std::fabs(A[(k - 1) + L * (k - 1)]);
			int imax = 0;
			double colmax = 0.0;
			if (k < n) {
				imax = k + idamax(n - k, A + k + L * (k - 1), 1) + 1;
				colmax = std::fabs(A[(imax - 1) + L * (k - 1)]);
			}
			if (std::max(absakk, colmax) == 0.0 || std::isnan(absakk)) {
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
					int jmax = k - 1 + idamax(imax - k, A + (imax - 1) + L * (k - 1), lda) + 1;
					double rowmax = std::fabs(A[(imax - 1) + L * (jmax - 1)]);
					if (imax < n) {
						jmax = imax + idamax(n - imax, A + imax + L * (imax - 1), 1) + 1;
						rowmax = std::max(rowmax, std::fabs(A[(jmax - 1) + L * (imax - 1)]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (std::fabs(A[(imax - 1) + L * (imax - 1)]) >= alpha * rowmax) {
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
						dswap(n - kp, A + kp + L * (kk - 1), 1, A + kp + L * (kp - 1), 1);
					}
					dswap(kp - kk - 1, A + kk + L * (kk - 1), 1, A + (kp - 1) + L * kk, lda);
					double t = A[(kk - 1) + L * (kk - 1)];
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
						const double d11 = 1.0 / A[(k - 1) + L * (k - 1)];
						rdsyr('L', n - k, -d11, A + k + L * (k - 1), 1, A + k + L * k, lda, rm);
						rdscal(n - k, d11, A + k + L * (k - 1), 1, rm);
					}
				}
				else {
					if (k < n - 1) {
						double d21 = A[k + L * (k - 1)];
						const double d11 = A[k + L * k] / d21;
						const double d22 = A[(k - 1) + L * (k - 1)] / d21;
						const double t = 1.0 / (d11 * d22 - 1.0);
						d21 = t / d21;
						for (int j = k + 2; j <= n; j++) {
							const double wk = d21 * (d11 * A[(j - 1) + L * (k - 1)] - A[(j - 1) + L * k]);
							const double wkp1 = d21 * (d22 * A[(j - 1) + L * k] - A[(j - 1) + L * (k - 1)]);
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
inline int rdlasyf(const char uplo, const int n, const int nb, int& kb, double* A, const int lda,
	int* ipiv, double* W, const int ldw, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LW = static_cast<std::size_t>(ldw);
	const double alpha = (1.0 + std::sqrt(17.0)) / 8.0;
	int info = 0;
	if (detail::option_is(uplo, 'U')) {
		int k = n;
		int kw = 0;
		for (;;) {
			kw = nb + k - n;
			if ((k <= n - nb + 1 && nb < n) || k < 1) {
				break;
			}
			dcopy(k, A + L * (k - 1), 1, W + LW * (kw - 1), 1);
			if (k < n) {
				rdgemv('N', k, n - k, -1.0, A + L * k, lda, W + (k - 1) + LW * kw, ldw,
					1.0, W + LW * (kw - 1), 1, rm);
			}
			int kstep = 1;
			int kp;
			const double absakk = std::fabs(W[(k - 1) + LW * (kw - 1)]);
			int imax = 0;
			double colmax = 0.0;
			if (k > 1) {
				imax = idamax(k - 1, W + LW * (kw - 1), 1) + 1;
				colmax = std::fabs(W[(imax - 1) + LW * (kw - 1)]);
			}
			if (std::max(absakk, colmax) == 0.0) {
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
					dcopy(imax, A + L * (imax - 1), 1, W + LW * (kw - 2), 1);
					dcopy(k - imax, A + (imax - 1) + L * imax, lda, W + imax + LW * (kw - 2), 1);
					if (k < n) {
						rdgemv('N', k, n - k, -1.0, A + L * k, lda, W + (imax - 1) + LW * kw, ldw,
							1.0, W + LW * (kw - 2), 1, rm);
					}
					int jmax = imax + idamax(k - imax, W + imax + LW * (kw - 2), 1) + 1;
					double rowmax = std::fabs(W[(jmax - 1) + LW * (kw - 2)]);
					if (imax > 1) {
						jmax = idamax(imax - 1, W + LW * (kw - 2), 1) + 1;
						rowmax = std::max(rowmax, std::fabs(W[(jmax - 1) + LW * (kw - 2)]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (std::fabs(W[(imax - 1) + LW * (kw - 2)]) >= alpha * rowmax) {
						kp = imax;
						dcopy(k, W + LW * (kw - 2), 1, W + LW * (kw - 1), 1);
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
					dcopy(kk - 1 - kp, A + kp + L * (kk - 1), 1, A + (kp - 1) + L * kp, lda);
					if (kp > 1) {
						dcopy(kp - 1, A + L * (kk - 1), 1, A + L * (kp - 1), 1);
					}
					if (k < n) {
						dswap(n - k, A + (kk - 1) + L * k, lda, A + (kp - 1) + L * k, lda);
					}
					dswap(n - kk + 1, W + (kk - 1) + LW * (kkw - 1), ldw, W + (kp - 1) + LW * (kkw - 1), ldw);
				}
				if (kstep == 1) {
					dcopy(k, W + LW * (kw - 1), 1, A + L * (k - 1), 1);
					const double r1 = 1.0 / A[(k - 1) + L * (k - 1)];
					rdscal(k - 1, r1, A + L * (k - 1), 1, rm);
				}
				else {
					if (k > 2) {
						double d21 = W[(k - 2) + LW * (kw - 1)];
						const double d11 = W[(k - 1) + LW * (kw - 1)] / d21;
						const double d22 = W[(k - 2) + LW * (kw - 2)] / d21;
						const double t = 1.0 / (d11 * d22 - 1.0);
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
				rdgemv('N', jj - j + 1, n - k, -1.0, A + (j - 1) + L * k, lda,
					W + (jj - 1) + LW * kw, ldw, 1.0, A + (j - 1) + L * (jj - 1), 1, rm);
			}
			if (j >= 2) {
				rdgemm('N', 'T', j - 1, jb, n - k, -1.0, A + L * k, lda,
					W + (j - 1) + LW * kw, ldw, 1.0, A + L * (j - 1), lda, rm);
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
					dswap(n - j + 1, A + (jp - 1) + L * (j - 1), lda, A + (jj - 1) + L * (j - 1), lda);
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
			dcopy(n - k + 1, A + (k - 1) + L * (k - 1), 1, W + (k - 1) + LW * (k - 1), 1);
			rdgemv('N', n - k + 1, k - 1, -1.0, A + (k - 1), lda, W + (k - 1), ldw,
				1.0, W + (k - 1) + LW * (k - 1), 1, rm);
			int kstep = 1;
			int kp;
			const double absakk = std::fabs(W[(k - 1) + LW * (k - 1)]);
			int imax = 0;
			double colmax = 0.0;
			if (k < n) {
				imax = k + idamax(n - k, W + k + LW * (k - 1), 1) + 1;
				colmax = std::fabs(W[(imax - 1) + LW * (k - 1)]);
			}
			if (std::max(absakk, colmax) == 0.0) {
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
					dcopy(imax - k, A + (imax - 1) + L * (k - 1), lda, W + (k - 1) + LW * k, 1);
					dcopy(n - imax + 1, A + (imax - 1) + L * (imax - 1), 1, W + (imax - 1) + LW * k, 1);
					rdgemv('N', n - k + 1, k - 1, -1.0, A + (k - 1), lda, W + (imax - 1), ldw,
						1.0, W + (k - 1) + LW * k, 1, rm);
					int jmax = k - 1 + idamax(imax - k, W + (k - 1) + LW * k, 1) + 1;
					double rowmax = std::fabs(W[(jmax - 1) + LW * k]);
					if (imax < n) {
						jmax = imax + idamax(n - imax, W + imax + LW * k, 1) + 1;
						rowmax = std::max(rowmax, std::fabs(W[(jmax - 1) + LW * k]));
					}
					if (absakk >= alpha * colmax * (colmax / rowmax)) {
						kp = k;
					}
					else if (std::fabs(W[(imax - 1) + LW * k]) >= alpha * rowmax) {
						kp = imax;
						dcopy(n - k + 1, W + (k - 1) + LW * k, 1, W + (k - 1) + LW * (k - 1), 1);
					}
					else {
						kp = imax;
						kstep = 2;
					}
				}
				const int kk = k + kstep - 1;
				if (kp != kk) {
					A[(kp - 1) + L * (kp - 1)] = A[(kk - 1) + L * (kk - 1)];
					dcopy(kp - kk - 1, A + kk + L * (kk - 1), 1, A + (kp - 1) + L * kk, lda);
					if (kp < n) {
						dcopy(n - kp, A + kp + L * (kk - 1), 1, A + kp + L * (kp - 1), 1);
					}
					if (k > 1) {
						dswap(k - 1, A + (kk - 1), lda, A + (kp - 1), lda);
					}
					dswap(kk, W + (kk - 1), ldw, W + (kp - 1), ldw);
				}
				if (kstep == 1) {
					dcopy(n - k + 1, W + (k - 1) + LW * (k - 1), 1, A + (k - 1) + L * (k - 1), 1);
					if (k < n) {
						const double r1 = 1.0 / A[(k - 1) + L * (k - 1)];
						rdscal(n - k, r1, A + k + L * (k - 1), 1, rm);
					}
				}
				else {
					if (k < n - 1) {
						double d21 = W[k + LW * (k - 1)];
						const double d11 = W[k + LW * k] / d21;
						const double d22 = W[(k - 1) + LW * (k - 1)] / d21;
						const double t = 1.0 / (d11 * d22 - 1.0);
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
				rdgemv('N', j + jb - jj, k - 1, -1.0, A + (jj - 1), lda, W + (jj - 1), ldw,
					1.0, A + (jj - 1) + L * (jj - 1), 1, rm);
			}
			if (j + jb <= n) {
				rdgemm('N', 'T', n - j - jb + 1, jb, k - 1, -1.0, A + (j + jb - 1), lda,
					W + (j - 1), ldw, 1.0, A + (j + jb - 1) + L * (j - 1), lda, rm);
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
					dswap(j, A + (jp - 1), lda, A + (jj - 1), lda);
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
inline int rdsytrf(const char uplo, const int n, double* A, const int lda, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsytrf: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdsytrf: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	int nb = detail::ilaenv(1, "SYTRF", n, -1, -1, -1);
	const int nbmin = detail::ilaenv(2, "SYTRF", n, -1, -1, -1);
	if (nb < nbmin || nb >= n) {
		nb = n;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int ldwork = n;
	std::vector<double> work(nb < n ? static_cast<std::size_t>(ldwork) * nb : 1);
	int info = 0;
	if (upper) {
		int k = n;
		while (k >= 1) {
			int kb;
			int iinfo;
			if (k > nb) {
				iinfo = rdlasyf(uplo, k, nb, kb, A, lda, ipiv, work.data(), ldwork, rounding_mode);
			}
			else {
				iinfo = rdsytf2(uplo, k, A, lda, ipiv, rounding_mode);
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
				iinfo = rdlasyf(uplo, n - k + 1, nb, kb, A + (k - 1) + L * (k - 1), lda,
					ipiv + (k - 1), work.data(), ldwork, rounding_mode);
			}
			else {
				iinfo = rdsytf2(uplo, n - k + 1, A + (k - 1) + L * (k - 1), lda, ipiv + (k - 1), rounding_mode);
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
inline int rdsytrs(const char uplo, const int n, const int nrhs, const double* A, const int lda,
	const int* ipiv, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsytrs: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdsytrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LB = static_cast<std::size_t>(ldb);
	if (upper) {
		// B := U^-T D^-1 U^-1 B (後退方向 → 前進方向)
		int k = n;
		while (k >= 1) {
			if (ipiv[k - 1] > 0) {
				const int kp = ipiv[k - 1];
				if (kp != k) {
					dswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				rdger(k - 1, nrhs, -1.0, A + L * (k - 1), 1, B + (k - 1), ldb, B, ldb, rm);
				rdscal(nrhs, 1.0 / A[(k - 1) + L * (k - 1)], B + (k - 1), ldb, rm);
				k = k - 1;
			}
			else {
				const int kp = -ipiv[k - 1];
				if (kp != k - 1) {
					dswap(nrhs, B + (k - 2), ldb, B + (kp - 1), ldb);
				}
				rdger(k - 2, nrhs, -1.0, A + L * (k - 1), 1, B + (k - 1), ldb, B, ldb, rm);
				rdger(k - 2, nrhs, -1.0, A + L * (k - 2), 1, B + (k - 2), ldb, B, ldb, rm);
				const double akm1k = A[(k - 2) + L * (k - 1)];
				const double akm1 = A[(k - 2) + L * (k - 2)] / akm1k;
				const double ak = A[(k - 1) + L * (k - 1)] / akm1k;
				const double denom = akm1 * ak - 1.0;
				for (int j = 1; j <= nrhs; j++) {
					const double bkm1 = B[(k - 2) + LB * (j - 1)] / akm1k;
					const double bk = B[(k - 1) + LB * (j - 1)] / akm1k;
					B[(k - 2) + LB * (j - 1)] = (ak * bkm1 - bk) / denom;
					B[(k - 1) + LB * (j - 1)] = (akm1 * bk - bkm1) / denom;
				}
				k = k - 2;
			}
		}
		k = 1;
		while (k <= n) {
			if (ipiv[k - 1] > 0) {
				rdgemv('T', k - 1, nrhs, -1.0, B, ldb, A + L * (k - 1), 1, 1.0, B + (k - 1), ldb, rm);
				const int kp = ipiv[k - 1];
				if (kp != k) {
					dswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k + 1;
			}
			else {
				rdgemv('T', k - 1, nrhs, -1.0, B, ldb, A + L * (k - 1), 1, 1.0, B + (k - 1), ldb, rm);
				rdgemv('T', k - 1, nrhs, -1.0, B, ldb, A + L * k, 1, 1.0, B + k, ldb, rm);
				const int kp = -ipiv[k - 1];
				if (kp != k) {
					dswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
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
					dswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				if (k < n) {
					rdger(n - k, nrhs, -1.0, A + k + L * (k - 1), 1, B + (k - 1), ldb, B + k, ldb, rm);
				}
				rdscal(nrhs, 1.0 / A[(k - 1) + L * (k - 1)], B + (k - 1), ldb, rm);
				k = k + 1;
			}
			else {
				const int kp = -ipiv[k - 1];
				if (kp != k + 1) {
					dswap(nrhs, B + k, ldb, B + (kp - 1), ldb);
				}
				if (k < n - 1) {
					rdger(n - k - 1, nrhs, -1.0, A + (k + 1) + L * (k - 1), 1, B + (k - 1), ldb,
						B + (k + 1), ldb, rm);
					rdger(n - k - 1, nrhs, -1.0, A + (k + 1) + L * k, 1, B + k, ldb,
						B + (k + 1), ldb, rm);
				}
				const double akm1k = A[k + L * (k - 1)];
				const double akm1 = A[(k - 1) + L * (k - 1)] / akm1k;
				const double ak = A[k + L * k] / akm1k;
				const double denom = akm1 * ak - 1.0;
				for (int j = 1; j <= nrhs; j++) {
					const double bkm1 = B[(k - 1) + LB * (j - 1)] / akm1k;
					const double bk = B[k + LB * (j - 1)] / akm1k;
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
					rdgemv('T', n - k, nrhs, -1.0, B + k, ldb, A + k + L * (k - 1), 1,
						1.0, B + (k - 1), ldb, rm);
				}
				const int kp = ipiv[k - 1];
				if (kp != k) {
					dswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k - 1;
			}
			else {
				if (k < n) {
					rdgemv('T', n - k, nrhs, -1.0, B + k, ldb, A + k + L * (k - 1), 1,
						1.0, B + (k - 1), ldb, rm);
					rdgemv('T', n - k, nrhs, -1.0, B + k, ldb, A + k + L * (k - 2), 1,
						1.0, B + (k - 2), ldb, rm);
				}
				const int kp = -ipiv[k - 1];
				if (kp != k) {
					dswap(nrhs, B + (k - 1), ldb, B + (kp - 1), ldb);
				}
				k = k - 2;
			}
		}
	}
	return 0;
}

// 対称不定値の連立一次方程式 A*X = B の求解 (driver)
inline int rdsysv(const char uplo, const int n, const int nrhs, double* A, const int lda,
	int* ipiv, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsysv: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdsysv: invalid argument");
	}
	const int info = rdsytrf(uplo, n, A, lda, ipiv, rounding_mode);
	if (info == 0) {
		rdsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, rounding_mode);
	}
	return info;
}

// 一般化対称固有値問題の標準形への変換 (unblocked)
// itype = 1: A := inv(U^T)*A*inv(U) など, 2/3: A := U*A*U^T など (B は rdpotrf 済み)
inline int rdsygs2(const int itype, const char uplo, const int n, double* A, const int lda,
	const double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::rdlapack_error("rdsygs2: invalid itype");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsygs2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdsygs2: invalid argument");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LBB = static_cast<std::size_t>(ldb);
	if (itype == 1) {
		if (upper) {
			for (int k = 1; k <= n; k++) {
				double akk = A[(k - 1) + L * (k - 1)];
				const double bkk = B[(k - 1) + LBB * (k - 1)];
				akk = akk / (bkk * bkk);
				A[(k - 1) + L * (k - 1)] = akk;
				if (k < n) {
					rdscal(n - k, 1.0 / bkk, A + (k - 1) + L * k, lda, rm);
					const double ct = -0.5 * akk;
					rdaxpy(n - k, ct, B + (k - 1) + LBB * k, ldb, A + (k - 1) + L * k, lda, rm);
					rdsyr2('U', n - k, -1.0, A + (k - 1) + L * k, lda, B + (k - 1) + LBB * k, ldb,
						A + k + L * k, lda, rm);
					rdaxpy(n - k, ct, B + (k - 1) + LBB * k, ldb, A + (k - 1) + L * k, lda, rm);
					rdtrsv('U', 'T', 'N', n - k, B + k + LBB * k, ldb, A + (k - 1) + L * k, lda, rm);
				}
			}
		}
		else {
			for (int k = 1; k <= n; k++) {
				double akk = A[(k - 1) + L * (k - 1)];
				const double bkk = B[(k - 1) + LBB * (k - 1)];
				akk = akk / (bkk * bkk);
				A[(k - 1) + L * (k - 1)] = akk;
				if (k < n) {
					rdscal(n - k, 1.0 / bkk, A + k + L * (k - 1), 1, rm);
					const double ct = -0.5 * akk;
					rdaxpy(n - k, ct, B + k + LBB * (k - 1), 1, A + k + L * (k - 1), 1, rm);
					rdsyr2('L', n - k, -1.0, A + k + L * (k - 1), 1, B + k + LBB * (k - 1), 1,
						A + k + L * k, lda, rm);
					rdaxpy(n - k, ct, B + k + LBB * (k - 1), 1, A + k + L * (k - 1), 1, rm);
					rdtrsv('L', 'N', 'N', n - k, B + k + LBB * k, ldb, A + k + L * (k - 1), 1, rm);
				}
			}
		}
	}
	else {
		if (upper) {
			for (int k = 1; k <= n; k++) {
				const double akk = A[(k - 1) + L * (k - 1)];
				const double bkk = B[(k - 1) + LBB * (k - 1)];
				rdtrmv('U', 'N', 'N', k - 1, B, ldb, A + L * (k - 1), 1, rm);
				const double ct = 0.5 * akk;
				rdaxpy(k - 1, ct, B + LBB * (k - 1), 1, A + L * (k - 1), 1, rm);
				rdsyr2('U', k - 1, 1.0, A + L * (k - 1), 1, B + LBB * (k - 1), 1, A, lda, rm);
				rdaxpy(k - 1, ct, B + LBB * (k - 1), 1, A + L * (k - 1), 1, rm);
				rdscal(k - 1, bkk, A + L * (k - 1), 1, rm);
				A[(k - 1) + L * (k - 1)] = akk * (bkk * bkk);
			}
		}
		else {
			for (int k = 1; k <= n; k++) {
				const double akk = A[(k - 1) + L * (k - 1)];
				const double bkk = B[(k - 1) + LBB * (k - 1)];
				rdtrmv('L', 'T', 'N', k - 1, B, ldb, A + (k - 1), lda, rm);
				const double ct = 0.5 * akk;
				rdaxpy(k - 1, ct, B + (k - 1), ldb, A + (k - 1), lda, rm);
				rdsyr2('L', k - 1, 1.0, A + (k - 1), lda, B + (k - 1), ldb, A, lda, rm);
				rdaxpy(k - 1, ct, B + (k - 1), ldb, A + (k - 1), lda, rm);
				rdscal(k - 1, bkk, A + (k - 1), lda, rm);
				A[(k - 1) + L * (k - 1)] = akk * (bkk * bkk);
			}
		}
	}
	return 0;
}

// 一般化対称固有値問題の標準形への変換 (blocked)
inline int rdsygst(const int itype, const char uplo, const int n, double* A, const int lda,
	const double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::rdlapack_error("rdsygst: invalid itype");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsygst: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdsygst: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LBB = static_cast<std::size_t>(ldb);
	const int nb = detail::ilaenv(1, "SYGST", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		return rdsygs2(itype, uplo, n, A, lda, B, ldb, rounding_mode);
	}
	if (itype == 1) {
		if (upper) {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				rdsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb, rm);
				if (k + kb <= n) {
					rdtrsm('L', 'U', 'T', 'N', kb, n - k - kb + 1, 1.0, B + (k - 1) + LBB * (k - 1), ldb,
						A + (k - 1) + L * (k + kb - 1), lda, rm);
					rdsymm('L', 'U', kb, n - k - kb + 1, -0.5, A + (k - 1) + L * (k - 1), lda,
						B + (k - 1) + LBB * (k + kb - 1), ldb, 1.0, A + (k - 1) + L * (k + kb - 1), lda, rm);
					rdsyr2k('U', 'T', n - k - kb + 1, kb, -1.0, A + (k - 1) + L * (k + kb - 1), lda,
						B + (k - 1) + LBB * (k + kb - 1), ldb, 1.0, A + (k + kb - 1) + L * (k + kb - 1), lda, rm);
					rdsymm('L', 'U', kb, n - k - kb + 1, -0.5, A + (k - 1) + L * (k - 1), lda,
						B + (k - 1) + LBB * (k + kb - 1), ldb, 1.0, A + (k - 1) + L * (k + kb - 1), lda, rm);
					rdtrsm('R', 'U', 'N', 'N', kb, n - k - kb + 1, 1.0,
						B + (k + kb - 1) + LBB * (k + kb - 1), ldb, A + (k - 1) + L * (k + kb - 1), lda, rm);
				}
			}
		}
		else {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				rdsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb, rm);
				if (k + kb <= n) {
					rdtrsm('R', 'L', 'T', 'N', n - k - kb + 1, kb, 1.0, B + (k - 1) + LBB * (k - 1), ldb,
						A + (k + kb - 1) + L * (k - 1), lda, rm);
					rdsymm('R', 'L', n - k - kb + 1, kb, -0.5, A + (k - 1) + L * (k - 1), lda,
						B + (k + kb - 1) + LBB * (k - 1), ldb, 1.0, A + (k + kb - 1) + L * (k - 1), lda, rm);
					rdsyr2k('L', 'N', n - k - kb + 1, kb, -1.0, A + (k + kb - 1) + L * (k - 1), lda,
						B + (k + kb - 1) + LBB * (k - 1), ldb, 1.0, A + (k + kb - 1) + L * (k + kb - 1), lda, rm);
					rdsymm('R', 'L', n - k - kb + 1, kb, -0.5, A + (k - 1) + L * (k - 1), lda,
						B + (k + kb - 1) + LBB * (k - 1), ldb, 1.0, A + (k + kb - 1) + L * (k - 1), lda, rm);
					rdtrsm('L', 'L', 'N', 'N', n - k - kb + 1, kb, 1.0,
						B + (k + kb - 1) + LBB * (k + kb - 1), ldb, A + (k + kb - 1) + L * (k - 1), lda, rm);
				}
			}
		}
	}
	else {
		if (upper) {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				rdtrmm('L', 'U', 'N', 'N', k - 1, kb, 1.0, B, ldb, A + L * (k - 1), lda, rm);
				rdsymm('R', 'U', k - 1, kb, 0.5, A + (k - 1) + L * (k - 1), lda,
					B + LBB * (k - 1), ldb, 1.0, A + L * (k - 1), lda, rm);
				rdsyr2k('U', 'N', k - 1, kb, 1.0, A + L * (k - 1), lda, B + LBB * (k - 1), ldb,
					1.0, A, lda, rm);
				rdsymm('R', 'U', k - 1, kb, 0.5, A + (k - 1) + L * (k - 1), lda,
					B + LBB * (k - 1), ldb, 1.0, A + L * (k - 1), lda, rm);
				rdtrmm('R', 'U', 'T', 'N', k - 1, kb, 1.0, B + (k - 1) + LBB * (k - 1), ldb,
					A + L * (k - 1), lda, rm);
				rdsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb, rm);
			}
		}
		else {
			for (int k = 1; k <= n; k += nb) {
				const int kb = std::min(n - k + 1, nb);
				rdtrmm('R', 'L', 'N', 'N', kb, k - 1, 1.0, B, ldb, A + (k - 1), lda, rm);
				rdsymm('L', 'L', kb, k - 1, 0.5, A + (k - 1) + L * (k - 1), lda,
					B + (k - 1), ldb, 1.0, A + (k - 1), lda, rm);
				rdsyr2k('L', 'T', k - 1, kb, 1.0, A + (k - 1), lda, B + (k - 1), ldb, 1.0, A, lda, rm);
				rdsymm('L', 'L', kb, k - 1, 0.5, A + (k - 1) + L * (k - 1), lda,
					B + (k - 1), ldb, 1.0, A + (k - 1), lda, rm);
				rdtrmm('L', 'L', 'T', 'N', kb, k - 1, 1.0, B + (k - 1) + LBB * (k - 1), ldb,
					A + (k - 1), lda, rm);
				rdsygs2(itype, uplo, kb, A + (k - 1) + L * (k - 1), lda, B + (k - 1) + LBB * (k - 1), ldb, rm);
			}
		}
	}
	return 0;
}

// 一般化対称固有値問題 A*x = lambda*B*x (itype=1), A*B*x = lambda*x (2),
// B*A*x = lambda*x (3) の driver (B は対称正定値，上書きされる)
inline int rdsygv(const int itype, const char jobz, const char uplo, const int n,
	double* A, const int lda, double* B, const int ldb, double* w, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::rdlapack_error("rdsygv: invalid itype");
	}
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::rdlapack_error("rdsygv: invalid jobz");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsygv: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdsygv: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	int info = rdpotrf(uplo, n, B, ldb, rounding_mode);
	if (info != 0) {
		return n + info;
	}
	rdsygst(itype, uplo, n, A, lda, B, ldb, rounding_mode);
	info = rdsyev(jobz, uplo, n, A, lda, w, rounding_mode);
	if (wantz) {
		int neig = n;
		if (info > 0) {
			neig = info - 1;
		}
		if (itype == 1 || itype == 2) {
			const char trans = upper ? 'N' : 'T';
			rdtrsm('L', uplo, trans, 'N', n, neig, 1.0, B, ldb, A, lda, rounding_mode);
		}
		else {
			const char trans = upper ? 'T' : 'N';
			rdtrmm('L', uplo, trans, 'N', n, neig, 1.0, B, ldb, A, lda, rounding_mode);
		}
	}
	return info;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_SY_HPP
