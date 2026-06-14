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

#ifndef VBLAS_RDLAPACK_EIG_HPP
#define VBLAS_RDLAPACK_EIG_HPP

// 対称固有値問題 (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dsytd2/dlatrd/dsytrd/dorgtr/dormtr/dsterf/
// dsteqr/dsyev を移植．
// dsterf/dsteqr は Fortran の goto 構造を保つため 1-based の局所 index で
// 記述している (配列 access は [idx-1])．

#include "rdlapack_qr.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// 対称行列の三重対角化 (unblocked)
inline int rdsytd2(const char uplo, const int n, double* A, const int lda,
	double* d, double* e, double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsytd2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdsytd2: invalid n/lda");
	}
	if (n <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int i = n - 2; i >= 0; i--) {
			double taui;
			rdlarfg(i + 1, A[i + L * (i + 1)], A + L * (i + 1), 1, taui, rounding_mode);
			e[i] = A[i + L * (i + 1)];
			if (taui != 0.0) {
				A[i + L * (i + 1)] = 1.0;
				rdsymv('U', i + 1, taui, A, lda, A + L * (i + 1), 1, 0.0, tau, 1, rounding_mode);
				const double alpha = -0.5 * taui * rddot(i + 1, tau, 1, A + L * (i + 1), 1, rounding_mode);
				rdaxpy(i + 1, alpha, A + L * (i + 1), 1, tau, 1, rounding_mode);
				rdsyr2('U', i + 1, -1.0, A + L * (i + 1), 1, tau, 1, A, lda, rounding_mode);
				A[i + L * (i + 1)] = e[i];
			}
			d[i + 1] = A[(i + 1) + L * (i + 1)];
			tau[i] = taui;
		}
		d[0] = A[0];
	}
	else {
		for (int i = 0; i < n - 1; i++) {
			double taui;
			rdlarfg(n - i - 1, A[(i + 1) + L * i], A + std::min(i + 2, n - 1) + L * i, 1, taui, rounding_mode);
			e[i] = A[(i + 1) + L * i];
			if (taui != 0.0) {
				A[(i + 1) + L * i] = 1.0;
				rdsymv('L', n - i - 1, taui, A + (i + 1) + L * (i + 1), lda, A + (i + 1) + L * i, 1,
					0.0, tau + i, 1, rounding_mode);
				const double alpha = -0.5 * taui * rddot(n - i - 1, tau + i, 1, A + (i + 1) + L * i, 1, rounding_mode);
				rdaxpy(n - i - 1, alpha, A + (i + 1) + L * i, 1, tau + i, 1, rounding_mode);
				rdsyr2('L', n - i - 1, -1.0, A + (i + 1) + L * i, 1, tau + i, 1,
					A + (i + 1) + L * (i + 1), lda, rounding_mode);
				A[(i + 1) + L * i] = e[i];
			}
			d[i] = A[i + L * i];
			tau[i] = taui;
		}
		d[n - 1] = A[(n - 1) + L * (n - 1)];
	}
	return 0;
}

// 三重対角化の先頭/末尾 nb 段 (blocked 用の panel 計算)
inline void rdlatrd(const char uplo, const int n, const int nb, double* A, const int lda,
	double* e, double* tau, double* W, const int ldw, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n <= 0) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LW = static_cast<std::size_t>(ldw);
	if (detail::option_is(uplo, 'U')) {
		for (int i = n - 1; i >= n - nb; i--) {
			const int iw = i - n + nb;
			if (i < n - 1) {
				rdgemv('N', i + 1, n - 1 - i, -1.0, A + L * (i + 1), lda, W + i + LW * (iw + 1), ldw,
					1.0, A + L * i, 1, rounding_mode);
				rdgemv('N', i + 1, n - 1 - i, -1.0, W + LW * (iw + 1), ldw, A + i + L * (i + 1), lda,
					1.0, A + L * i, 1, rounding_mode);
			}
			if (i > 0) {
				rdlarfg(i, A[(i - 1) + L * i], A + L * i, 1, tau[i - 1], rounding_mode);
				e[i - 1] = A[(i - 1) + L * i];
				A[(i - 1) + L * i] = 1.0;
				rdsymv('U', i, 1.0, A, lda, A + L * i, 1, 0.0, W + LW * iw, 1, rounding_mode);
				if (i < n - 1) {
					rdgemv('T', i, n - 1 - i, 1.0, W + LW * (iw + 1), ldw, A + L * i, 1,
						0.0, W + (i + 1) + LW * iw, 1, rounding_mode);
					rdgemv('N', i, n - 1 - i, -1.0, A + L * (i + 1), lda, W + (i + 1) + LW * iw, 1,
						1.0, W + LW * iw, 1, rounding_mode);
					rdgemv('T', i, n - 1 - i, 1.0, A + L * (i + 1), lda, A + L * i, 1,
						0.0, W + (i + 1) + LW * iw, 1, rounding_mode);
					rdgemv('N', i, n - 1 - i, -1.0, W + LW * (iw + 1), ldw, W + (i + 1) + LW * iw, 1,
						1.0, W + LW * iw, 1, rounding_mode);
				}
				rdscal(i, tau[i - 1], W + LW * iw, 1, rounding_mode);
				const double alpha = -0.5 * tau[i - 1] * rddot(i, W + LW * iw, 1, A + L * i, 1, rounding_mode);
				rdaxpy(i, alpha, A + L * i, 1, W + LW * iw, 1, rounding_mode);
			}
		}
	}
	else {
		for (int i = 0; i < nb; i++) {
			rdgemv('N', n - i, i, -1.0, A + i, lda, W + i, ldw, 1.0, A + i + L * i, 1, rounding_mode);
			rdgemv('N', n - i, i, -1.0, W + i, ldw, A + i, lda, 1.0, A + i + L * i, 1, rounding_mode);
			if (i < n - 1) {
				rdlarfg(n - i - 1, A[(i + 1) + L * i], A + std::min(i + 2, n - 1) + L * i, 1, tau[i], rounding_mode);
				e[i] = A[(i + 1) + L * i];
				A[(i + 1) + L * i] = 1.0;
				rdsymv('L', n - i - 1, 1.0, A + (i + 1) + L * (i + 1), lda, A + (i + 1) + L * i, 1,
					0.0, W + (i + 1) + LW * i, 1, rounding_mode);
				rdgemv('T', n - i - 1, i, 1.0, W + (i + 1), ldw, A + (i + 1) + L * i, 1,
					0.0, W + LW * i, 1, rounding_mode);
				rdgemv('N', n - i - 1, i, -1.0, A + (i + 1), lda, W + LW * i, 1,
					1.0, W + (i + 1) + LW * i, 1, rounding_mode);
				rdgemv('T', n - i - 1, i, 1.0, A + (i + 1), lda, A + (i + 1) + L * i, 1,
					0.0, W + LW * i, 1, rounding_mode);
				rdgemv('N', n - i - 1, i, -1.0, W + (i + 1), ldw, W + LW * i, 1,
					1.0, W + (i + 1) + LW * i, 1, rounding_mode);
				rdscal(n - i - 1, tau[i], W + (i + 1) + LW * i, 1, rounding_mode);
				const double alpha = -0.5 * tau[i] * rddot(n - i - 1, W + (i + 1) + LW * i, 1,
					A + (i + 1) + L * i, 1, rounding_mode);
				rdaxpy(n - i - 1, alpha, A + (i + 1) + L * i, 1, W + (i + 1) + LW * i, 1, rounding_mode);
			}
		}
	}
}

// 対称行列の三重対角化 (blocked)
inline int rdsytrd(const char uplo, const int n, double* A, const int lda,
	double* d, double* e, double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdsytrd: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdsytrd: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	int nb = detail::ilaenv(1, "SYTRD", n, -1, -1, -1);
	int nx = n;
	if (nb > 1 && nb < n) {
		nx = std::max(nb, detail::ilaenv(3, "SYTRD", n, -1, -1, -1));
		if (nx >= n) {
			nb = 1;
		}
	}
	else {
		nb = 1;
	}
	if (nb == 1 || nx >= n) {
		return rdsytd2(uplo, n, A, lda, d, e, tau, rounding_mode);
	}
	const int ldwork = n;
	std::vector<double> work(static_cast<std::size_t>(ldwork) * nb);
	if (upper) {
		const int kk = n - ((n - nx + nb - 1) / nb) * nb;
		for (int i = n - nb; i >= kk; i -= nb) {
			rdlatrd('U', i + nb, nb, A, lda, e, tau, work.data(), ldwork, rounding_mode);
			rdsyr2k('U', 'N', i, nb, -1.0, A + L * i, lda, work.data(), ldwork, 1.0, A, lda, rounding_mode);
			for (int j = i; j < i + nb; j++) {
				A[(j - 1) + L * j] = e[j - 1];
				d[j] = A[j + L * j];
			}
		}
		rdsytd2('U', kk, A, lda, d, e, tau, rounding_mode);
	}
	else {
		int i = 0;
		for (; i < n - nx; i += nb) {
			rdlatrd('L', n - i, nb, A + i + L * i, lda, e + i, tau + i, work.data(), ldwork, rounding_mode);
			rdsyr2k('L', 'N', n - i - nb, nb, -1.0, A + (i + nb) + L * i, lda,
				work.data() + nb, ldwork, 1.0, A + (i + nb) + L * (i + nb), lda, rounding_mode);
			for (int j = i; j < i + nb; j++) {
				A[(j + 1) + L * j] = e[j];
				d[j] = A[j + L * j];
			}
		}
		rdsytd2('L', n - i, A + i + L * i, lda, d + i, e + i, tau + i, rounding_mode);
	}
	return 0;
}

// rdsytrd の Q の生成
inline int rdorgtr(const char uplo, const int n, double* A, const int lda,
	const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdorgtr: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdorgtr: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int j = 0; j < n - 1; j++) {
			for (int i = 0; i < j; i++) {
				A[i + L * j] = A[i + L * (j + 1)];
			}
			A[(n - 1) + L * j] = 0.0;
		}
		for (int i = 0; i < n - 1; i++) {
			A[i + L * (n - 1)] = 0.0;
		}
		A[(n - 1) + L * (n - 1)] = 1.0;
		rdorgql(n - 1, n - 1, n - 1, A, lda, tau, rounding_mode);
	}
	else {
		for (int j = n - 1; j >= 1; j--) {
			A[L * j] = 0.0;
			for (int i = j + 1; i < n; i++) {
				A[i + L * j] = A[i + L * (j - 1)];
			}
		}
		A[0] = 1.0;
		for (int i = 1; i < n; i++) {
			A[i] = 0.0;
		}
		if (n > 1) {
			rdorgqr(n - 1, n - 1, n - 1, A + 1 + L, lda, tau, rounding_mode);
		}
	}
	return 0;
}

// rdsytrd の Q の適用
inline int rdormtr(const char side, const char uplo, const char trans, const int m, const int n,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool upper = detail::option_is(uplo, 'U');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!upper && !detail::option_is(uplo, 'L')) ||
		(!detail::option_is(trans, 'N') && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdormtr: invalid argument");
	}
	if (m == 0 || n == 0 || nq == 1) {
		return 0;
	}
	const int mi = left ? m - 1 : m;
	const int ni = left ? n : n - 1;
	if (upper) {
		rdormql(side, trans, mi, ni, nq - 1, A + static_cast<std::size_t>(lda), lda, tau, C, ldc, rounding_mode);
	}
	else {
		double* C1 = C + (left ? 1 : static_cast<std::size_t>(ldc));
		rdormqr(side, trans, mi, ni, nq - 1, A + 1, lda, tau, C1, ldc, rounding_mode);
	}
	return 0;
}

// 対称三重対角行列の全固有値 (Pal-Walker-Kahan QL/QR, 固有 vector なし)
// d, e は 1-based の Fortran 流 index で記述 (access は [idx-1])
inline int rdsterf(const int n, double* d, double* e, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0) {
		detail::rdlapack_error("rdsterf: invalid n");
	}
	if (n <= 1) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	int info = 0;
	const double eps = DBL_EPSILON * 0.5;
	const double eps2 = eps * eps;
	const double safmin = DBL_MIN;
	const double safmax = 1.0 / safmin;
	const double ssfmax = std::sqrt(safmax) / 3.0;
	const double ssfmin = std::sqrt(safmin) / eps2;
	const int nmaxit = n * 30;
	int jtot = 0;
	int l1 = 1;
	int l = 0, lsv = 0, lend = 0, lendsv = 0, iscale = 0;
	double anorm = 0.0;

L10:
	if (l1 > n) {
		dlasrt('I', n, d);
		return info;
	}
	if (l1 > 1) {
		e[l1 - 2] = 0.0;
	}
	{
		int m = n;
		for (int mm = l1; mm <= n - 1; mm++) {
			if (std::fabs(e[mm - 1]) <= (std::sqrt(std::fabs(d[mm - 1])) * std::sqrt(std::fabs(d[mm]))) * eps) {
				e[mm - 1] = 0.0;
				m = mm;
				break;
			}
		}
		l = l1;
		lsv = l;
		lend = m;
		lendsv = lend;
		l1 = m + 1;
		if (lend == l) {
			goto L10;
		}
	}
	anorm = rdlanst('M', lend - l + 1, d + (l - 1), e + (l - 1), rm);
	iscale = 0;
	if (anorm == 0.0) {
		goto L10;
	}
	if (anorm > ssfmax) {
		iscale = 1;
		rdlascl('G', 0, 0, anorm, ssfmax, lend - l + 1, 1, d + (l - 1), n, rm);
		rdlascl('G', 0, 0, anorm, ssfmax, lend - l, 1, e + (l - 1), n, rm);
	}
	else if (anorm < ssfmin) {
		iscale = 2;
		rdlascl('G', 0, 0, anorm, ssfmin, lend - l + 1, 1, d + (l - 1), n, rm);
		rdlascl('G', 0, 0, anorm, ssfmin, lend - l, 1, e + (l - 1), n, rm);
	}
	for (int i = l; i <= lend - 1; i++) {
		e[i - 1] = e[i - 1] * e[i - 1];
	}
	if (std::fabs(d[lend - 1]) < std::fabs(d[l - 1])) {
		lend = lsv;
		l = lendsv;
	}
	if (lend >= l) {
		// QL 反復 (下から上)
	L50:
		{
			int m = lend;
			if (l != lend) {
				for (int mm = l; mm <= lend - 1; mm++) {
					if (std::fabs(e[mm - 1]) <= eps2 * std::fabs(d[mm - 1] * d[mm])) {
						m = mm;
						break;
					}
				}
			}
			if (m < lend) {
				e[m - 1] = 0.0;
			}
			double p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l + 1;
				if (l <= lend) {
					goto L50;
				}
				goto L150;
			}
			if (m == l + 1) {
				double rt1, rt2;
				const double rte = std::sqrt(e[l - 1]);
				rdlae2(d[l - 1], rte, d[l], rt1, rt2, rm);
				d[l - 1] = rt1;
				d[l] = rt2;
				e[l - 1] = 0.0;
				l = l + 2;
				if (l <= lend) {
					goto L50;
				}
				goto L150;
			}
			if (jtot == nmaxit) {
				goto L150;
			}
			jtot++;
			const double rte = std::sqrt(e[l - 1]);
			double sigma = (d[l] - p) / (2.0 * rte);
			double r = rdlapy2(sigma, 1.0, rm);
			sigma = p - (rte / (sigma + detail::f_sign(r, sigma)));
			double c = 1.0;
			double s = 0.0;
			double gamma = d[m - 1] - sigma;
			p = gamma * gamma;
			for (int i = m - 1; i >= l; i--) {
				const double bb = e[i - 1];
				r = p + bb;
				if (i != m - 1) {
					e[i] = s * r;
				}
				const double oldc = c;
				c = p / r;
				s = bb / r;
				const double oldgam = gamma;
				const double alpha = d[i - 1];
				gamma = c * (alpha - sigma) - s * oldgam;
				d[i] = oldgam + (alpha - gamma);
				if (c != 0.0) {
					p = (gamma * gamma) / c;
				}
				else {
					p = oldc * bb;
				}
			}
			e[l - 1] = s * p;
			d[l - 1] = sigma + gamma;
			goto L50;
		}
	}
	else {
		// QR 反復 (上から下)
	L100:
		{
			int m = lend;
			for (int mm = l; mm >= lend + 1; mm--) {
				if (std::fabs(e[mm - 2]) <= eps2 * std::fabs(d[mm - 1] * d[mm - 2])) {
					m = mm;
					break;
				}
			}
			if (m > lend) {
				e[m - 2] = 0.0;
			}
			double p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l - 1;
				if (l >= lend) {
					goto L100;
				}
				goto L150;
			}
			if (m == l - 1) {
				double rt1, rt2;
				const double rte = std::sqrt(e[l - 2]);
				rdlae2(d[l - 1], rte, d[l - 2], rt1, rt2, rm);
				d[l - 1] = rt1;
				d[l - 2] = rt2;
				e[l - 2] = 0.0;
				l = l - 2;
				if (l >= lend) {
					goto L100;
				}
				goto L150;
			}
			if (jtot == nmaxit) {
				goto L150;
			}
			jtot++;
			const double rte = std::sqrt(e[l - 2]);
			double sigma = (d[l - 2] - p) / (2.0 * rte);
			double r = rdlapy2(sigma, 1.0, rm);
			sigma = p - (rte / (sigma + detail::f_sign(r, sigma)));
			double c = 1.0;
			double s = 0.0;
			double gamma = d[m - 1] - sigma;
			p = gamma * gamma;
			for (int i = m; i <= l - 1; i++) {
				const double bb = e[i - 1];
				r = p + bb;
				if (i != m) {
					e[i - 2] = s * r;
				}
				const double oldc = c;
				c = p / r;
				s = bb / r;
				const double oldgam = gamma;
				const double alpha = d[i];
				gamma = c * (alpha - sigma) - s * oldgam;
				d[i - 1] = oldgam + (alpha - gamma);
				if (c != 0.0) {
					p = (gamma * gamma) / c;
				}
				else {
					p = oldc * bb;
				}
			}
			e[l - 2] = s * p;
			d[l - 1] = sigma + gamma;
			goto L100;
		}
	}

L150:
	if (iscale == 1) {
		rdlascl('G', 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n, rm);
	}
	if (iscale == 2) {
		rdlascl('G', 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n, rm);
	}
	if (jtot < nmaxit) {
		goto L10;
	}
	for (int i = 1; i <= n - 1; i++) {
		if (e[i - 1] != 0.0) {
			info++;
		}
	}
	return info;
}

// 対称三重対角行列の固有値・固有 vector (implicit QL/QR)
// compz: 'N' 固有値のみ, 'V' Z に三重対角化の Q を渡して元行列の固有 vector,
//        'I' Z を単位行列で初期化して三重対角行列の固有 vector
inline int rdsteqr(const char compz, const int n, double* d, double* e,
	double* Z, const int ldz, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	int icompz;
	if (detail::option_is(compz, 'N')) {
		icompz = 0;
	}
	else if (detail::option_is(compz, 'V')) {
		icompz = 1;
	}
	else if (detail::option_is(compz, 'I')) {
		icompz = 2;
	}
	else {
		detail::rdlapack_error("rdsteqr: invalid compz");
		return -1;
	}
	if (n < 0 || ldz < 1 || (icompz > 0 && ldz < std::max(1, n))) {
		detail::rdlapack_error("rdsteqr: invalid n/ldz");
	}
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		if (icompz == 2) {
			Z[0] = 1.0;
		}
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t LZ = static_cast<std::size_t>(ldz);
	int info = 0;
	const double eps = DBL_EPSILON * 0.5;
	const double eps2 = eps * eps;
	const double safmin = DBL_MIN;
	const double safmax = 1.0 / safmin;
	const double ssfmax = std::sqrt(safmax) / 3.0;
	const double ssfmin = std::sqrt(safmin) / eps2;
	if (icompz == 2) {
		dlaset('F', n, n, 0.0, 1.0, Z, ldz);
	}
	std::vector<double> work(icompz > 0 ? static_cast<std::size_t>(2 * n - 2) : 1);
	const int nmaxit = n * 30;
	int jtot = 0;
	int l1 = 1;
	int l = 0, lsv = 0, lend = 0, lendsv = 0, iscale = 0;
	double anorm = 0.0;

L10:
	if (l1 > n) {
		goto L160;
	}
	if (l1 > 1) {
		e[l1 - 2] = 0.0;
	}
	{
		int m = n;
		if (l1 <= n - 1) {
			for (int mm = l1; mm <= n - 1; mm++) {
				const double tst = std::fabs(e[mm - 1]);
				if (tst == 0.0) {
					m = mm;
					break;
				}
				if (tst <= (std::sqrt(std::fabs(d[mm - 1])) * std::sqrt(std::fabs(d[mm]))) * eps) {
					e[mm - 1] = 0.0;
					m = mm;
					break;
				}
			}
		}
		l = l1;
		lsv = l;
		lend = m;
		lendsv = lend;
		l1 = m + 1;
		if (lend == l) {
			goto L10;
		}
	}
	anorm = rdlanst('M', lend - l + 1, d + (l - 1), e + (l - 1), rm);
	iscale = 0;
	if (anorm == 0.0) {
		goto L10;
	}
	if (anorm > ssfmax) {
		iscale = 1;
		rdlascl('G', 0, 0, anorm, ssfmax, lend - l + 1, 1, d + (l - 1), n, rm);
		rdlascl('G', 0, 0, anorm, ssfmax, lend - l, 1, e + (l - 1), n, rm);
	}
	else if (anorm < ssfmin) {
		iscale = 2;
		rdlascl('G', 0, 0, anorm, ssfmin, lend - l + 1, 1, d + (l - 1), n, rm);
		rdlascl('G', 0, 0, anorm, ssfmin, lend - l, 1, e + (l - 1), n, rm);
	}
	if (std::fabs(d[lend - 1]) < std::fabs(d[l - 1])) {
		lend = lsv;
		l = lendsv;
	}
	if (lend > l) {
		// QL 反復 (下から上)
	L40:
		{
			int m = lend;
			if (l != lend) {
				for (int mm = l; mm <= lend - 1; mm++) {
					const double tst = std::fabs(e[mm - 1]) * std::fabs(e[mm - 1]);
					if (tst <= (eps2 * std::fabs(d[mm - 1])) * std::fabs(d[mm]) + safmin) {
						m = mm;
						break;
					}
				}
			}
			if (m < lend) {
				e[m - 1] = 0.0;
			}
			double p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l + 1;
				if (l <= lend) {
					goto L40;
				}
				goto L140;
			}
			if (m == l + 1) {
				double rt1, rt2, c, s;
				if (icompz > 0) {
					rdlaev2(d[l - 1], e[l - 1], d[l], rt1, rt2, c, s, rm);
					work[l - 1] = c;
					work[n - 1 + l - 1] = s;
					rdlasr('R', 'V', 'B', n, 2, work.data() + (l - 1), work.data() + (n - 1 + l - 1),
						Z + LZ * (l - 1), ldz, rm);
				}
				else {
					rdlae2(d[l - 1], e[l - 1], d[l], rt1, rt2, rm);
				}
				d[l - 1] = rt1;
				d[l] = rt2;
				e[l - 1] = 0.0;
				l = l + 2;
				if (l <= lend) {
					goto L40;
				}
				goto L140;
			}
			if (jtot == nmaxit) {
				goto L140;
			}
			jtot++;
			double g = (d[l] - p) / (2.0 * e[l - 1]);
			double r = rdlapy2(g, 1.0, rm);
			g = d[m - 1] - p + (e[l - 1] / (g + detail::f_sign(r, g)));
			double s = 1.0;
			double c = 1.0;
			p = 0.0;
			for (int i = m - 1; i >= l; i--) {
				double f = s * e[i - 1];
				const double b = c * e[i - 1];
				rdlartg(g, f, c, s, r, rm);
				if (i != m - 1) {
					e[i] = r;
				}
				g = d[i] - p;
				r = (d[i - 1] - g) * s + 2.0 * c * b;
				p = s * r;
				d[i] = g + p;
				g = c * r - b;
				if (icompz > 0) {
					work[i - 1] = c;
					work[n - 1 + i - 1] = -s;
				}
			}
			if (icompz > 0) {
				const int mm = m - l + 1;
				rdlasr('R', 'V', 'B', n, mm, work.data() + (l - 1), work.data() + (n - 1 + l - 1),
					Z + LZ * (l - 1), ldz, rm);
			}
			d[l - 1] = d[l - 1] - p;
			e[l - 1] = g;
			goto L40;
		}
	}
	else {
		// QR 反復 (上から下)
	L90:
		{
			int m = lend;
			if (l != lend) {
				for (int mm = l; mm >= lend + 1; mm--) {
					const double tst = std::fabs(e[mm - 2]) * std::fabs(e[mm - 2]);
					if (tst <= (eps2 * std::fabs(d[mm - 1])) * std::fabs(d[mm - 2]) + safmin) {
						m = mm;
						break;
					}
				}
			}
			if (m > lend) {
				e[m - 2] = 0.0;
			}
			double p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l - 1;
				if (l >= lend) {
					goto L90;
				}
				goto L140;
			}
			if (m == l - 1) {
				double rt1, rt2, c, s;
				if (icompz > 0) {
					rdlaev2(d[l - 2], e[l - 2], d[l - 1], rt1, rt2, c, s, rm);
					work[m - 1] = c;
					work[n - 1 + m - 1] = s;
					rdlasr('R', 'V', 'F', n, 2, work.data() + (m - 1), work.data() + (n - 1 + m - 1),
						Z + LZ * (l - 2), ldz, rm);
				}
				else {
					rdlae2(d[l - 2], e[l - 2], d[l - 1], rt1, rt2, rm);
				}
				d[l - 2] = rt1;
				d[l - 1] = rt2;
				e[l - 2] = 0.0;
				l = l - 2;
				if (l >= lend) {
					goto L90;
				}
				goto L140;
			}
			if (jtot == nmaxit) {
				goto L140;
			}
			jtot++;
			double g = (d[l - 2] - p) / (2.0 * e[l - 2]);
			double r = rdlapy2(g, 1.0, rm);
			g = d[m - 1] - p + (e[l - 2] / (g + detail::f_sign(r, g)));
			double s = 1.0;
			double c = 1.0;
			p = 0.0;
			for (int i = m; i <= l - 1; i++) {
				double f = s * e[i - 1];
				const double b = c * e[i - 1];
				rdlartg(g, f, c, s, r, rm);
				if (i != m) {
					e[i - 2] = r;
				}
				g = d[i - 1] - p;
				r = (d[i] - g) * s + 2.0 * c * b;
				p = s * r;
				d[i - 1] = g + p;
				g = c * r - b;
				if (icompz > 0) {
					work[i - 1] = c;
					work[n - 1 + i - 1] = s;
				}
			}
			if (icompz > 0) {
				const int mm = l - m + 1;
				rdlasr('R', 'V', 'F', n, mm, work.data() + (m - 1), work.data() + (n - 1 + m - 1),
					Z + LZ * (m - 1), ldz, rm);
			}
			d[l - 1] = d[l - 1] - p;
			e[l - 2] = g;
			goto L90;
		}
	}

L140:
	if (iscale == 1) {
		rdlascl('G', 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n, rm);
		rdlascl('G', 0, 0, ssfmax, anorm, lendsv - lsv, 1, e + (lsv - 1), n, rm);
	}
	else if (iscale == 2) {
		rdlascl('G', 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n, rm);
		rdlascl('G', 0, 0, ssfmin, anorm, lendsv - lsv, 1, e + (lsv - 1), n, rm);
	}
	if (jtot < nmaxit) {
		goto L10;
	}
	for (int i = 1; i <= n - 1; i++) {
		if (e[i - 1] != 0.0) {
			info++;
		}
	}
	return info;

L160:
	// 固有値を昇順に sort (vector 付きは選択 sort で Z の列も入れ替え)
	if (icompz == 0) {
		dlasrt('I', n, d);
	}
	else {
		for (int ii = 2; ii <= n; ii++) {
			const int i = ii - 1;
			int k = i;
			double p = d[i - 1];
			for (int j = ii; j <= n; j++) {
				if (d[j - 1] < p) {
					k = j;
					p = d[j - 1];
				}
			}
			if (k != i) {
				d[k - 1] = d[i - 1];
				d[i - 1] = p;
				dswap(n, Z + LZ * (i - 1), 1, Z + LZ * (k - 1), 1);
			}
		}
	}
	return info;
}

// 対称行列の固有値・固有 vector (driver)
// jobz: 'N' 固有値のみ, 'V' 固有 vector も計算 (A に上書き)
// w: 固有値 (昇順, n 個)
inline int rdsyev(const char jobz, const char uplo, const int n, double* A, const int lda,
	double* w, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	const bool lower = detail::option_is(uplo, 'L');
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::rdlapack_error("rdsyev: invalid jobz");
	}
	if (!lower && !detail::option_is(uplo, 'U')) {
		detail::rdlapack_error("rdsyev: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdsyev: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		w[0] = A[0];
		if (wantz) {
			A[0] = 1.0;
		}
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const double safmin = DBL_MIN;
	const double eps = DBL_EPSILON;
	const double smlnum = safmin / eps;
	const double bignum = 1.0 / smlnum;
	const double rmin = std::sqrt(smlnum);
	const double rmax = std::sqrt(bignum);
	const double anrm = rdlansy('M', uplo, n, A, lda, rm);
	int iscale = 0;
	double sigma = 1.0;
	if (anrm > 0.0 && anrm < rmin) {
		iscale = 1;
		sigma = rmin / anrm;
	}
	else if (anrm > rmax) {
		iscale = 1;
		sigma = rmax / anrm;
	}
	if (iscale == 1) {
		rdlascl(uplo, 0, 0, 1.0, sigma, n, n, A, lda, rm);
	}
	std::vector<double> e(static_cast<std::size_t>(n));
	std::vector<double> tau(static_cast<std::size_t>(n));
	rdsytrd(uplo, n, A, lda, w, e.data(), tau.data(), rm);
	int info;
	if (!wantz) {
		info = rdsterf(n, w, e.data(), rm);
	}
	else {
		rdorgtr(uplo, n, A, lda, tau.data(), rm);
		info = rdsteqr(jobz, n, w, e.data(), A, lda, rm);
	}
	if (iscale == 1) {
		const int imax = (info == 0) ? n : info - 1;
		rdscal(imax, 1.0 / sigma, w, 1, rm);
	}
	return info;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_EIG_HPP
