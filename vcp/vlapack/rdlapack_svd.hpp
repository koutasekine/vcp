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

#ifndef VBLAS_RDLAPACK_SVD_HPP
#define VBLAS_RDLAPACK_SVD_HPP

// 特異値分解 (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dgebd2/dlabrd/dgebrd/dorgbr/dormbr/dbdsqr を移植．
//
// reference からの変更点 (ドキュメントにも記載):
// - rdbdsqr は特異値のみ (ncvt = nru = ncc = 0) の場合も dlasq1 (dqds) ではなく
//   QR 反復で計算する (dlasq 系は移植していない．結果の意味は同じ)．
// - rdgesvd は基本経路 (A を直接二重対角化) のみで，m >> n / n >> m 用の
//   QR/LQ 前処理経路 (mnthr 分岐) は持たない．結果は同等で演算量のみ異なる．

#include "rdlapack_eig.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// 二重対角化 (unblocked): A = Q * B * P^T
inline int rdgebd2(const int m, const int n, double* A, const int lda,
	double* d, double* e, double* tauq, double* taup, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgebd2: invalid argument");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (m >= n) {
		for (int i = 0; i < n; i++) {
			rdlarfg(m - i, A[i + L * i], A + std::min(i + 1, m - 1) + L * i, 1, tauq[i], rounding_mode);
			d[i] = A[i + L * i];
			if (i < n - 1) {
				rdlarf1f('L', m - i, n - i - 1, A + i + L * i, 1, tauq[i], A + i + L * (i + 1), lda, rounding_mode);
				rdlarfg(n - i - 1, A[i + L * (i + 1)], A + i + L * std::min(i + 2, n - 1), lda, taup[i], rounding_mode);
				e[i] = A[i + L * (i + 1)];
				rdlarf1f('R', m - i - 1, n - i - 1, A + i + L * (i + 1), lda, taup[i],
					A + (i + 1) + L * (i + 1), lda, rounding_mode);
			}
			else {
				taup[i] = 0.0;
			}
		}
	}
	else {
		for (int i = 0; i < m; i++) {
			rdlarfg(n - i, A[i + L * i], A + i + L * std::min(i + 1, n - 1), lda, taup[i], rounding_mode);
			d[i] = A[i + L * i];
			if (i < m - 1) {
				rdlarf1f('R', m - i - 1, n - i, A + i + L * i, lda, taup[i], A + (i + 1) + L * i, lda, rounding_mode);
				rdlarfg(m - i - 1, A[(i + 1) + L * i], A + std::min(i + 2, m - 1) + L * i, 1, tauq[i], rounding_mode);
				e[i] = A[(i + 1) + L * i];
				rdlarf1f('L', m - i - 1, n - i - 1, A + (i + 1) + L * i, 1, tauq[i],
					A + (i + 1) + L * (i + 1), lda, rounding_mode);
			}
			else {
				tauq[i] = 0.0;
			}
		}
	}
	return 0;
}

// 二重対角化の先頭 nb 段 (blocked 用 panel, X/Y は更新行列)
inline void rdlabrd(const int m, const int n, const int nb, double* A, const int lda,
	double* d, double* e, double* tauq, double* taup,
	double* X, const int ldx, double* Y, const int ldy, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m <= 0 || n <= 0) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LX = static_cast<std::size_t>(ldx);
	const std::size_t LY = static_cast<std::size_t>(ldy);
	const int rm = rounding_mode;
	if (m >= n) {
		for (int i = 0; i < nb; i++) {
			rdgemv('N', m - i, i, -1.0, A + i, lda, Y + i, ldy, 1.0, A + i + L * i, 1, rm);
			rdgemv('N', m - i, i, -1.0, X + i, ldx, A + L * i, 1, 1.0, A + i + L * i, 1, rm);
			rdlarfg(m - i, A[i + L * i], A + std::min(i + 1, m - 1) + L * i, 1, tauq[i], rm);
			d[i] = A[i + L * i];
			if (i < n - 1) {
				A[i + L * i] = 1.0;
				rdgemv('T', m - i, n - i - 1, 1.0, A + i + L * (i + 1), lda, A + i + L * i, 1,
					0.0, Y + (i + 1) + LY * i, 1, rm);
				rdgemv('T', m - i, i, 1.0, A + i, lda, A + i + L * i, 1, 0.0, Y + LY * i, 1, rm);
				rdgemv('N', n - i - 1, i, -1.0, Y + (i + 1), ldy, Y + LY * i, 1,
					1.0, Y + (i + 1) + LY * i, 1, rm);
				rdgemv('T', m - i, i, 1.0, X + i, ldx, A + i + L * i, 1, 0.0, Y + LY * i, 1, rm);
				rdgemv('T', i, n - i - 1, -1.0, A + L * (i + 1), lda, Y + LY * i, 1,
					1.0, Y + (i + 1) + LY * i, 1, rm);
				rdscal(n - i - 1, tauq[i], Y + (i + 1) + LY * i, 1, rm);
				rdgemv('N', n - i - 1, i + 1, -1.0, Y + (i + 1), ldy, A + i, lda,
					1.0, A + i + L * (i + 1), lda, rm);
				rdgemv('T', i, n - i - 1, -1.0, A + L * (i + 1), lda, X + i, ldx,
					1.0, A + i + L * (i + 1), lda, rm);
				rdlarfg(n - i - 1, A[i + L * (i + 1)], A + i + L * std::min(i + 2, n - 1), lda, taup[i], rm);
				e[i] = A[i + L * (i + 1)];
				A[i + L * (i + 1)] = 1.0;
				rdgemv('N', m - i - 1, n - i - 1, 1.0, A + (i + 1) + L * (i + 1), lda,
					A + i + L * (i + 1), lda, 0.0, X + (i + 1) + LX * i, 1, rm);
				rdgemv('T', n - i - 1, i + 1, 1.0, Y + (i + 1), ldy, A + i + L * (i + 1), lda,
					0.0, X + LX * i, 1, rm);
				rdgemv('N', m - i - 1, i + 1, -1.0, A + (i + 1), lda, X + LX * i, 1,
					1.0, X + (i + 1) + LX * i, 1, rm);
				rdgemv('N', i, n - i - 1, 1.0, A + L * (i + 1), lda, A + i + L * (i + 1), lda,
					0.0, X + LX * i, 1, rm);
				rdgemv('N', m - i - 1, i, -1.0, X + (i + 1), ldx, X + LX * i, 1,
					1.0, X + (i + 1) + LX * i, 1, rm);
				rdscal(m - i - 1, taup[i], X + (i + 1) + LX * i, 1, rm);
			}
		}
	}
	else {
		for (int i = 0; i < nb; i++) {
			rdgemv('N', n - i, i, -1.0, Y + i, ldy, A + i, lda, 1.0, A + i + L * i, lda, rm);
			rdgemv('T', i, n - i, -1.0, A + L * i, lda, X + i, ldx, 1.0, A + i + L * i, lda, rm);
			rdlarfg(n - i, A[i + L * i], A + i + L * std::min(i + 1, n - 1), lda, taup[i], rm);
			d[i] = A[i + L * i];
			if (i < m - 1) {
				A[i + L * i] = 1.0;
				rdgemv('N', m - i - 1, n - i, 1.0, A + (i + 1) + L * i, lda, A + i + L * i, lda,
					0.0, X + (i + 1) + LX * i, 1, rm);
				rdgemv('T', n - i, i, 1.0, Y + i, ldy, A + i + L * i, lda, 0.0, X + LX * i, 1, rm);
				rdgemv('N', m - i - 1, i, -1.0, A + (i + 1), lda, X + LX * i, 1,
					1.0, X + (i + 1) + LX * i, 1, rm);
				rdgemv('N', i, n - i, 1.0, A + L * i, lda, A + i + L * i, lda, 0.0, X + LX * i, 1, rm);
				rdgemv('N', m - i - 1, i, -1.0, X + (i + 1), ldx, X + LX * i, 1,
					1.0, X + (i + 1) + LX * i, 1, rm);
				rdscal(m - i - 1, taup[i], X + (i + 1) + LX * i, 1, rm);
				rdgemv('N', m - i - 1, i, -1.0, A + (i + 1), lda, Y + i, ldy,
					1.0, A + (i + 1) + L * i, 1, rm);
				rdgemv('N', m - i - 1, i + 1, -1.0, X + (i + 1), ldx, A + L * i, 1,
					1.0, A + (i + 1) + L * i, 1, rm);
				rdlarfg(m - i - 1, A[(i + 1) + L * i], A + std::min(i + 2, m - 1) + L * i, 1, tauq[i], rm);
				e[i] = A[(i + 1) + L * i];
				A[(i + 1) + L * i] = 1.0;
				rdgemv('T', m - i - 1, n - i - 1, 1.0, A + (i + 1) + L * (i + 1), lda,
					A + (i + 1) + L * i, 1, 0.0, Y + (i + 1) + LY * i, 1, rm);
				rdgemv('T', m - i - 1, i, 1.0, A + (i + 1), lda, A + (i + 1) + L * i, 1,
					0.0, Y + LY * i, 1, rm);
				rdgemv('N', n - i - 1, i, -1.0, Y + (i + 1), ldy, Y + LY * i, 1,
					1.0, Y + (i + 1) + LY * i, 1, rm);
				rdgemv('T', m - i - 1, i + 1, 1.0, X + (i + 1), ldx, A + (i + 1) + L * i, 1,
					0.0, Y + LY * i, 1, rm);
				rdgemv('T', i + 1, n - i - 1, -1.0, A + L * (i + 1), lda, Y + LY * i, 1,
					1.0, Y + (i + 1) + LY * i, 1, rm);
				rdscal(n - i - 1, tauq[i], Y + (i + 1) + LY * i, 1, rm);
			}
		}
	}
}

// 二重対角化 (blocked)
inline int rdgebrd(const int m, const int n, double* A, const int lda,
	double* d, double* e, double* tauq, double* taup, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgebrd: invalid argument");
	}
	const int minmn = std::min(m, n);
	if (minmn == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	int nb = std::max(1, detail::ilaenv(1, "GEBRD", m, n, -1, -1));
	int nx = minmn;
	if (nb > 1 && nb < minmn) {
		nx = std::max(nb, detail::ilaenv(3, "GEBRD", m, n, -1, -1));
	}
	int i = 0;
	if (nx < minmn) {
		const int ldwrkx = m;
		const int ldwrky = n;
		std::vector<double> X(static_cast<std::size_t>(ldwrkx) * nb);
		std::vector<double> Y(static_cast<std::size_t>(ldwrky) * nb);
		for (; i < minmn - nx; i += nb) {
			rdlabrd(m - i, n - i, nb, A + i + L * i, lda, d + i, e + i, tauq + i, taup + i,
				X.data(), ldwrkx, Y.data(), ldwrky, rounding_mode);
			rdgemm('N', 'T', m - i - nb, n - i - nb, nb, -1.0, A + (i + nb) + L * i, lda,
				Y.data() + nb, ldwrky, 1.0, A + (i + nb) + L * (i + nb), lda, rounding_mode);
			rdgemm('N', 'N', m - i - nb, n - i - nb, nb, -1.0, X.data() + nb, ldwrkx,
				A + i + L * (i + nb), lda, 1.0, A + (i + nb) + L * (i + nb), lda, rounding_mode);
			if (m >= n) {
				for (int j = i; j < i + nb; j++) {
					A[j + L * j] = d[j];
					A[j + L * (j + 1)] = e[j];
				}
			}
			else {
				for (int j = i; j < i + nb; j++) {
					A[j + L * j] = d[j];
					A[(j + 1) + L * j] = e[j];
				}
			}
		}
	}
	rdgebd2(m - i, n - i, A + i + L * i, lda, d + i, e + i, tauq + i, taup + i, rounding_mode);
	return 0;
}

// rdgebrd の Q ('Q') または P^T ('P') の生成
inline int rdorgbr(const char vect, const int m, const int n, const int k,
	double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool wantq = detail::option_is(vect, 'Q');
	if (!wantq && !detail::option_is(vect, 'P')) {
		detail::rdlapack_error("rdorgbr: invalid vect");
	}
	if (m < 0 || k < 0 || lda < std::max(1, m) ||
		(wantq && (n > m || n < std::min(m, k))) ||
		(!wantq && (m > n || m < std::min(n, k)))) {
		detail::rdlapack_error("rdorgbr: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (wantq) {
		if (m >= k) {
			rdorgqr(m, n, k, A, lda, tau, rounding_mode);
		}
		else {
			for (int j = m - 1; j >= 1; j--) {
				A[L * j] = 0.0;
				for (int i = j + 1; i < m; i++) {
					A[i + L * j] = A[i + L * (j - 1)];
				}
			}
			A[0] = 1.0;
			for (int i = 1; i < m; i++) {
				A[i] = 0.0;
			}
			if (m > 1) {
				rdorgqr(m - 1, m - 1, m - 1, A + 1 + L, lda, tau, rounding_mode);
			}
		}
	}
	else {
		if (k < n) {
			rdorglq(m, n, k, A, lda, tau, rounding_mode);
		}
		else {
			A[0] = 1.0;
			for (int i = 1; i < n; i++) {
				A[i] = 0.0;
			}
			for (int j = 1; j < n; j++) {
				for (int i = j - 1; i >= 1; i--) {
					A[i + L * j] = A[(i - 1) + L * j];
				}
				A[L * j] = 0.0;
			}
			if (n > 1) {
				rdorglq(n - 1, n - 1, n - 1, A + 1 + L, lda, tau, rounding_mode);
			}
		}
	}
	return 0;
}

// rdgebrd の Q ('Q') または P^T ('P') の適用
inline int rdormbr(const char vect, const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool applyq = detail::option_is(vect, 'Q');
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!applyq && !detail::option_is(vect, 'P')) || (!left && !detail::option_is(side, 'R')) ||
		(!notran && !detail::option_is(trans, 'T')) || m < 0 || n < 0 || k < 0 ||
		(applyq && lda < std::max(1, nq)) || (!applyq && lda < std::max(1, std::min(nq, k))) ||
		ldc < std::max(1, m)) {
		detail::rdlapack_error("rdormbr: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	if (applyq) {
		if (nq >= k) {
			rdormqr(side, trans, m, n, k, A, lda, tau, C, ldc, rounding_mode);
		}
		else if (nq > 1) {
			const int mi = left ? m - 1 : m;
			const int ni = left ? n : n - 1;
			double* C1 = C + (left ? 1 : LC);
			rdormqr(side, trans, mi, ni, nq - 1, A + 1, lda, tau, C1, ldc, rounding_mode);
		}
	}
	else {
		const char transt = notran ? 'T' : 'N';
		if (nq > k) {
			rdormlq(side, transt, m, n, k, A, lda, tau, C, ldc, rounding_mode);
		}
		else if (nq > 1) {
			const int mi = left ? m - 1 : m;
			const int ni = left ? n : n - 1;
			double* C1 = C + (left ? 1 : LC);
			rdormlq(side, transt, mi, ni, nq - 1, A + L, lda, tau, C1, ldc, rounding_mode);
		}
	}
	return 0;
}

// 二重対角行列の SVD (implicit QR / 0-shift QR)．
// d, e は 1-based の Fortran 流 index で記述 (access は [idx-1])．
// 注: 特異値のみの場合も dlasq1 ではなく QR 反復を使う (reference との相違)．
inline int rdbdsqr(const char uplo, const int n, const int ncvt, const int nru, const int ncc,
	double* d, double* e, double* VT, const int ldvt, double* U, const int ldu,
	double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool lower = detail::option_is(uplo, 'L');
	if (!detail::option_is(uplo, 'U') && !lower) {
		detail::rdlapack_error("rdbdsqr: invalid uplo");
	}
	if (n < 0 || ncvt < 0 || nru < 0 || ncc < 0 ||
		(ncvt == 0 && ldvt < 1) || (ncvt > 0 && ldvt < std::max(1, n)) ||
		ldu < std::max(1, nru) ||
		(ncc == 0 && ldc < 1) || (ncc > 0 && ldc < std::max(1, n))) {
		detail::rdlapack_error("rdbdsqr: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t LVT = static_cast<std::size_t>(ldvt);
	const std::size_t LU = static_cast<std::size_t>(ldu);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	int info = 0;
	const int maxitr = 6;

	if (n > 1) {
		const int nm1 = n - 1;
		const int nm12 = nm1 + nm1;
		const int nm13 = nm12 + nm1;
		std::vector<double> work(static_cast<std::size_t>(4) * n);
		double* wk = work.data(); // wk[X-1] が Fortran の WORK(X)
		int idir = 0;
		const double eps = DBL_EPSILON;
		const double unfl = DBL_MIN;
		if (lower) {
			// 下二重対角を回転で上二重対角へ変換
			for (int i = 1; i <= n - 1; i++) {
				double cs, sn, r;
				rdlartg(d[i - 1], e[i - 1], cs, sn, r, rm);
				d[i - 1] = r;
				e[i - 1] = sn * d[i];
				d[i] = cs * d[i];
				wk[i - 1] = cs;
				wk[nm1 + i - 1] = sn;
			}
			if (nru > 0) {
				rdlasr('R', 'V', 'F', nru, n, wk, wk + (n - 1), U, ldu, rm);
			}
			if (ncc > 0) {
				rdlasr('L', 'V', 'F', n, ncc, wk, wk + (n - 1), C, ldc, rm);
			}
		}
		const double tolmul = std::max(10.0, std::min(100.0, std::pow(eps, -0.125)));
		const double tol = tolmul * eps;
		double smax = 0.0;
		for (int i = 1; i <= n; i++) {
			smax = std::max(smax, std::fabs(d[i - 1]));
		}
		for (int i = 1; i <= n - 1; i++) {
			smax = std::max(smax, std::fabs(e[i - 1]));
		}
		double smin = 0.0;
		double thresh;
		{
			double sminoa = std::fabs(d[0]);
			if (sminoa != 0.0) {
				double mu = sminoa;
				for (int i = 2; i <= n; i++) {
					mu = std::fabs(d[i - 1]) * (mu / (mu + std::fabs(e[i - 2])));
					sminoa = std::min(sminoa, mu);
					if (sminoa == 0.0) {
						break;
					}
				}
			}
			sminoa = sminoa / std::sqrt(static_cast<double>(n));
			thresh = std::max(tol * sminoa, maxitr * (n * (n * unfl)));
		}
		const int maxitdivn = maxitr * n;
		int iterdivn = 0;
		int iter = -1;
		int oldll = -1;
		int oldm = -1;
		int m = n;
		int ll = 0;

	L60:
		if (m <= 1) {
			goto L160;
		}
		if (iter >= n) {
			iter -= n;
			iterdivn++;
			if (iterdivn >= maxitdivn) {
				// 収束失敗
				for (int i = 1; i <= n - 1; i++) {
					if (e[i - 1] != 0.0) {
						info++;
					}
				}
				return info;
			}
		}
		smax = std::fabs(d[m - 1]);
		{
			bool found = false;
			for (int lll = 1; lll <= m - 1; lll++) {
				ll = m - lll;
				const double abss = std::fabs(d[ll - 1]);
				const double abse = std::fabs(e[ll - 1]);
				if (abse <= thresh) {
					found = true;
					break;
				}
				smax = std::max(std::max(smax, abss), abse);
			}
			if (!found) {
				ll = 0;
			}
			else {
				e[ll - 1] = 0.0;
				if (ll == m - 1) {
					m = m - 1;
					goto L60;
				}
			}
		}
		ll = ll + 1;
		if (ll == m - 1) {
			// 2x2 block は直接 SVD
			double sigmn, sigmx, sinr, cosr, sinl, cosl;
			rdlasv2(d[m - 2], e[m - 2], d[m - 1], sigmn, sigmx, sinr, cosr, sinl, cosl, rm);
			d[m - 2] = sigmx;
			e[m - 2] = 0.0;
			d[m - 1] = sigmn;
			if (ncvt > 0) {
				rdrot(ncvt, VT + (m - 2), ldvt, VT + (m - 1), ldvt, cosr, sinr, rm);
			}
			if (nru > 0) {
				rdrot(nru, U + LU * (m - 2), 1, U + LU * (m - 1), 1, cosl, sinl, rm);
			}
			if (ncc > 0) {
				rdrot(ncc, C + (m - 2), ldc, C + (m - 1), ldc, cosl, sinl, rm);
			}
			m = m - 2;
			goto L60;
		}
		if (ll > oldm || m < oldll) {
			idir = (std::fabs(d[ll - 1]) >= std::fabs(d[m - 1])) ? 1 : 2;
		}
		if (idir == 1) {
			if (std::fabs(e[m - 2]) <= std::fabs(tol) * std::fabs(d[m - 1])) {
				e[m - 2] = 0.0;
				goto L60;
			}
			{
				double mu = std::fabs(d[ll - 1]);
				smin = mu;
				for (int lll = ll; lll <= m - 1; lll++) {
					if (std::fabs(e[lll - 1]) <= tol * mu) {
						e[lll - 1] = 0.0;
						goto L60;
					}
					mu = std::fabs(d[lll]) * (mu / (mu + std::fabs(e[lll - 1])));
					smin = std::min(smin, mu);
				}
			}
		}
		else {
			if (std::fabs(e[ll - 1]) <= std::fabs(tol) * std::fabs(d[ll - 1])) {
				e[ll - 1] = 0.0;
				goto L60;
			}
			{
				double mu = std::fabs(d[m - 1]);
				smin = mu;
				for (int lll = m - 1; lll >= ll; lll--) {
					if (std::fabs(e[lll - 1]) <= tol * mu) {
						e[lll - 1] = 0.0;
						goto L60;
					}
					mu = std::fabs(d[lll - 1]) * (mu / (mu + std::fabs(e[lll - 1])));
					smin = std::min(smin, mu);
				}
			}
		}
		oldll = ll;
		oldm = m;
		{
			double shift = 0.0;
			double r;
			if (!(n * tol * (smin / smax) <= std::max(eps, 0.01 * tol))) {
				double sll;
				if (idir == 1) {
					sll = std::fabs(d[ll - 1]);
					rdlas2(d[m - 2], e[m - 2], d[m - 1], shift, r, rm);
				}
				else {
					sll = std::fabs(d[m - 1]);
					rdlas2(d[ll - 1], e[ll - 1], d[ll], shift, r, rm);
				}
				if (sll > 0.0) {
					if ((shift / sll) * (shift / sll) < eps) {
						shift = 0.0;
					}
				}
			}
			iter += m - ll;
			if (shift == 0.0) {
				// 0-shift QR (高精度)
				if (idir == 1) {
					double cs = 1.0;
					double oldcs = 1.0;
					double sn = 0.0, oldsn = 0.0;
					for (int i = ll; i <= m - 1; i++) {
						rdlartg(d[i - 1] * cs, e[i - 1], cs, sn, r, rm);
						if (i > ll) {
							e[i - 2] = oldsn * r;
						}
						rdlartg(oldcs * r, d[i] * sn, oldcs, oldsn, d[i - 1], rm);
						wk[i - ll] = cs;
						wk[i - ll + nm1] = sn;
						wk[i - ll + nm12] = oldcs;
						wk[i - ll + nm13] = oldsn;
					}
					const double h = d[m - 1] * cs;
					d[m - 1] = h * oldcs;
					e[m - 2] = h * oldsn;
					if (ncvt > 0) {
						rdlasr('L', 'V', 'F', m - ll + 1, ncvt, wk, wk + (n - 1), VT + (ll - 1), ldvt, rm);
					}
					if (nru > 0) {
						rdlasr('R', 'V', 'F', nru, m - ll + 1, wk + nm12, wk + nm13, U + LU * (ll - 1), ldu, rm);
					}
					if (ncc > 0) {
						rdlasr('L', 'V', 'F', m - ll + 1, ncc, wk + nm12, wk + nm13, C + (ll - 1), ldc, rm);
					}
					if (std::fabs(e[m - 2]) <= thresh) {
						e[m - 2] = 0.0;
					}
				}
				else {
					double cs = 1.0;
					double oldcs = 1.0;
					double sn = 0.0, oldsn = 0.0;
					for (int i = m; i >= ll + 1; i--) {
						rdlartg(d[i - 1] * cs, e[i - 2], cs, sn, r, rm);
						if (i < m) {
							e[i - 1] = oldsn * r;
						}
						rdlartg(oldcs * r, d[i - 2] * sn, oldcs, oldsn, d[i - 1], rm);
						wk[i - ll - 1] = cs;
						wk[i - ll - 1 + nm1] = -sn;
						wk[i - ll - 1 + nm12] = oldcs;
						wk[i - ll - 1 + nm13] = -oldsn;
					}
					const double h = d[ll - 1] * cs;
					d[ll - 1] = h * oldcs;
					e[ll - 1] = h * oldsn;
					if (ncvt > 0) {
						rdlasr('L', 'V', 'B', m - ll + 1, ncvt, wk + nm12, wk + nm13, VT + (ll - 1), ldvt, rm);
					}
					if (nru > 0) {
						rdlasr('R', 'V', 'B', nru, m - ll + 1, wk, wk + (n - 1), U + LU * (ll - 1), ldu, rm);
					}
					if (ncc > 0) {
						rdlasr('L', 'V', 'B', m - ll + 1, ncc, wk, wk + (n - 1), C + (ll - 1), ldc, rm);
					}
					if (std::fabs(e[ll - 1]) <= thresh) {
						e[ll - 1] = 0.0;
					}
				}
			}
			else {
				// shift 付き implicit QR
				if (idir == 1) {
					double f = (std::fabs(d[ll - 1]) - shift) *
						(detail::f_sign(1.0, d[ll - 1]) + shift / d[ll - 1]);
					double g = e[ll - 1];
					double cosr = 0.0, sinr = 0.0, cosl = 0.0, sinl = 0.0;
					for (int i = ll; i <= m - 1; i++) {
						rdlartg(f, g, cosr, sinr, r, rm);
						if (i > ll) {
							e[i - 2] = r;
						}
						f = cosr * d[i - 1] + sinr * e[i - 1];
						e[i - 1] = cosr * e[i - 1] - sinr * d[i - 1];
						g = sinr * d[i];
						d[i] = cosr * d[i];
						rdlartg(f, g, cosl, sinl, r, rm);
						d[i - 1] = r;
						f = cosl * e[i - 1] + sinl * d[i];
						d[i] = cosl * d[i] - sinl * e[i - 1];
						if (i < m - 1) {
							g = sinl * e[i];
							e[i] = cosl * e[i];
						}
						wk[i - ll] = cosr;
						wk[i - ll + nm1] = sinr;
						wk[i - ll + nm12] = cosl;
						wk[i - ll + nm13] = sinl;
					}
					e[m - 2] = f;
					if (ncvt > 0) {
						rdlasr('L', 'V', 'F', m - ll + 1, ncvt, wk, wk + (n - 1), VT + (ll - 1), ldvt, rm);
					}
					if (nru > 0) {
						rdlasr('R', 'V', 'F', nru, m - ll + 1, wk + nm12, wk + nm13, U + LU * (ll - 1), ldu, rm);
					}
					if (ncc > 0) {
						rdlasr('L', 'V', 'F', m - ll + 1, ncc, wk + nm12, wk + nm13, C + (ll - 1), ldc, rm);
					}
					if (std::fabs(e[m - 2]) <= thresh) {
						e[m - 2] = 0.0;
					}
				}
				else {
					double f = (std::fabs(d[m - 1]) - shift) *
						(detail::f_sign(1.0, d[m - 1]) + shift / d[m - 1]);
					double g = e[m - 2];
					double cosr = 0.0, sinr = 0.0, cosl = 0.0, sinl = 0.0;
					for (int i = m; i >= ll + 1; i--) {
						rdlartg(f, g, cosr, sinr, r, rm);
						if (i < m) {
							e[i - 1] = r;
						}
						f = cosr * d[i - 1] + sinr * e[i - 2];
						e[i - 2] = cosr * e[i - 2] - sinr * d[i - 1];
						g = sinr * d[i - 2];
						d[i - 2] = cosr * d[i - 2];
						rdlartg(f, g, cosl, sinl, r, rm);
						d[i - 1] = r;
						f = cosl * e[i - 2] + sinl * d[i - 2];
						d[i - 2] = cosl * d[i - 2] - sinl * e[i - 2];
						if (i > ll + 1) {
							g = sinl * e[i - 3];
							e[i - 3] = cosl * e[i - 3];
						}
						wk[i - ll - 1] = cosr;
						wk[i - ll - 1 + nm1] = -sinr;
						wk[i - ll - 1 + nm12] = cosl;
						wk[i - ll - 1 + nm13] = -sinl;
					}
					e[ll - 1] = f;
					if (std::fabs(e[ll - 1]) <= thresh) {
						e[ll - 1] = 0.0;
					}
					if (ncvt > 0) {
						rdlasr('L', 'V', 'B', m - ll + 1, ncvt, wk + nm12, wk + nm13, VT + (ll - 1), ldvt, rm);
					}
					if (nru > 0) {
						rdlasr('R', 'V', 'B', nru, m - ll + 1, wk, wk + (n - 1), U + LU * (ll - 1), ldu, rm);
					}
					if (ncc > 0) {
						rdlasr('L', 'V', 'B', m - ll + 1, ncc, wk, wk + (n - 1), C + (ll - 1), ldc, rm);
					}
				}
			}
		}
		goto L60;
	}

L160:
	// 特異値を正にして降順に sort
	for (int i = 1; i <= n; i++) {
		if (d[i - 1] < 0.0) {
			d[i - 1] = -d[i - 1];
			if (ncvt > 0) {
				rdscal(ncvt, -1.0, VT + (i - 1), ldvt, rm);
			}
		}
	}
	for (int i = 1; i <= n - 1; i++) {
		int isub = 1;
		double smin2 = d[0];
		for (int j = 2; j <= n + 1 - i; j++) {
			if (d[j - 1] <= smin2) {
				isub = j;
				smin2 = d[j - 1];
			}
		}
		if (isub != n + 1 - i) {
			d[isub - 1] = d[n - i];
			d[n - i] = smin2;
			if (ncvt > 0) {
				dswap(ncvt, VT + (isub - 1), ldvt, VT + (n - i), ldvt);
			}
			if (nru > 0) {
				dswap(nru, U + LU * (isub - 1), 1, U + LU * (n - i), 1);
			}
			if (ncc > 0) {
				dswap(ncc, C + (isub - 1), ldc, C + (n - i), ldc);
			}
		}
	}
	return info;
}

// 特異値分解 A = U * S * V^T (driver)
// jobu/jobvt: 'A' 全列, 'S' min(m,n) 列, 'O' A に上書き, 'N' 計算しない
// (jobu と jobvt を同時に 'O' にはできない)
// 注: reference の m >> n 用 QR 前処理経路は持たない基本経路のみの実装．
inline int rdgesvd(const char jobu, const char jobvt, const int m, const int n,
	double* A, const int lda, double* s, double* U, const int ldu,
	double* VT, const int ldvt, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const int minmn = std::min(m, n);
	const bool wntua = detail::option_is(jobu, 'A');
	const bool wntus = detail::option_is(jobu, 'S');
	const bool wntuo = detail::option_is(jobu, 'O');
	const bool wntun = detail::option_is(jobu, 'N');
	const bool wntva = detail::option_is(jobvt, 'A');
	const bool wntvs = detail::option_is(jobvt, 'S');
	const bool wntvo = detail::option_is(jobvt, 'O');
	const bool wntvn = detail::option_is(jobvt, 'N');
	if (!(wntua || wntus || wntuo || wntun) || !(wntva || wntvs || wntvo || wntvn) ||
		(wntuo && wntvo)) {
		detail::rdlapack_error("rdgesvd: invalid jobu/jobvt");
	}
	if (m < 0 || n < 0 || lda < std::max(1, m) ||
		((wntua || wntus) && ldu < std::max(1, m)) ||
		(wntva && ldvt < std::max(1, n)) || (wntvs && ldvt < std::max(1, minmn))) {
		detail::rdlapack_error("rdgesvd: invalid argument");
	}
	if (minmn == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	// scaling (dgesvd と同じ)
	const double eps = DBL_EPSILON * 0.5;
	const double smlnum = std::sqrt(DBL_MIN) / eps;
	const double bignum = 1.0 / smlnum;
	const double anrm = rdlange('M', m, n, A, lda, rm);
	int iscale = 0;
	if (anrm > 0.0 && anrm < smlnum) {
		iscale = 1;
		rdlascl('G', 0, 0, anrm, smlnum, m, n, A, lda, rm);
	}
	else if (anrm > bignum) {
		iscale = 2;
		rdlascl('G', 0, 0, anrm, bignum, m, n, A, lda, rm);
	}
	// 二重対角化
	std::vector<double> d(static_cast<std::size_t>(minmn));
	std::vector<double> e(static_cast<std::size_t>(std::max(1, minmn)));
	std::vector<double> tauq(static_cast<std::size_t>(minmn));
	std::vector<double> taup(static_cast<std::size_t>(minmn));
	rdgebrd(m, n, A, lda, d.data(), e.data(), tauq.data(), taup.data(), rm);
	// U の構築
	const bool wantu = wntua || wntus || wntuo;
	const bool wantvt = wntva || wntvs || wntvo;
	std::vector<double> ubuf, vbuf;
	double* Umat = U;
	int ldumat = ldu;
	int ncu = 0;
	if (wantu) {
		ncu = wntua ? m : minmn;
		if (wntuo) {
			ubuf.assign(static_cast<std::size_t>(m) * minmn, 0.0);
			Umat = ubuf.data();
			ldumat = m;
		}
		dlacpy('L', m, minmn, A, lda, Umat, ldumat);
		rdorgbr('Q', m, ncu, n, Umat, ldumat, tauq.data(), rm);
	}
	double* Vmat = VT;
	int ldvmat = ldvt;
	int nrvt = 0;
	if (wantvt) {
		nrvt = wntva ? n : minmn;
		if (wntvo) {
			vbuf.assign(static_cast<std::size_t>(minmn) * n, 0.0);
			Vmat = vbuf.data();
			ldvmat = minmn;
		}
		dlacpy('U', minmn, n, A, lda, Vmat, ldvmat);
		rdorgbr('P', nrvt, n, m, Vmat, ldvmat, taup.data(), rm);
	}
	// 二重対角 SVD
	const char bduplo = (m >= n) ? 'U' : 'L';
	double dummy = 0.0;
	const int info = rdbdsqr(bduplo, minmn, wantvt ? n : 0, wantu ? m : 0, 0,
		d.data(), e.data(), wantvt ? Vmat : &dummy, wantvt ? ldvmat : 1,
		wantu ? Umat : &dummy, wantu ? ldumat : 1, &dummy, 1, rm);
	for (int i = 0; i < minmn; i++) {
		s[i] = d[i];
	}
	// 'O' の書き戻し
	if (wntuo) {
		dlacpy('A', m, minmn, Umat, ldumat, A, lda);
	}
	if (wntvo) {
		dlacpy('A', minmn, n, Vmat, ldvmat, A, lda);
	}
	// scaling を戻す
	if (iscale == 1) {
		rdlascl('G', 0, 0, smlnum, anrm, minmn, 1, s, minmn, rm);
		if (info != 0) {
			rdlascl('G', 0, 0, smlnum, anrm, minmn - 1, 1, e.data(), minmn, rm);
		}
	}
	else if (iscale == 2) {
		rdlascl('G', 0, 0, bignum, anrm, minmn, 1, s, minmn, rm);
		if (info != 0) {
			rdlascl('G', 0, 0, bignum, anrm, minmn - 1, 1, e.data(), minmn, rm);
		}
	}
	return info;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_SVD_HPP
