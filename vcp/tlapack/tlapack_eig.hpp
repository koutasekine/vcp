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

#ifndef TLAPACK_TLAPACK_EIG_HPP
#define TLAPACK_TLAPACK_EIG_HPP

// 対称固有値問題 (template 版)．
// reference LAPACK 3.12.1 の dsytd2/dlatrd/dsytrd/dorgtr/dormtr/dsterf/
// dsteqr/dsyev を移植．
// dsterf/dsteqr は Fortran の goto 構造を保つため 1-based の局所 index で
// 記述している (配列 access は [idx-1])．

#include "tlapack_qr.hpp"


// 対称行列の三重対角化 (unblocked)
template <typename T>
inline int tsytd2(const char uplo, const int n, T* A, const int lda,
	T* d, T* e, T* tau) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsytd2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsytd2: invalid n/lda");
	}
	if (n <= 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int i = n - 2; i >= 0; i--) {
			T taui;
			tlarfg(i + 1, A[i + L * (i + 1)], A + L * (i + 1), 1, taui);
			e[i] = A[i + L * (i + 1)];
			if (taui != T(0)) {
				A[i + L * (i + 1)] = T(1);
				tsymv('U', i + 1, taui, A, lda, A + L * (i + 1), 1, T(0), tau, 1);
				const T alpha = -(T(1) / T(2)) * taui * tdot(i + 1, tau, 1, A + L * (i + 1), 1);
				taxpy(i + 1, alpha, A + L * (i + 1), 1, tau, 1);
				tsyr2('U', i + 1, -T(1), A + L * (i + 1), 1, tau, 1, A, lda);
				A[i + L * (i + 1)] = e[i];
			}
			d[i + 1] = A[(i + 1) + L * (i + 1)];
			tau[i] = taui;
		}
		d[0] = A[0];
	}
	else {
		for (int i = 0; i < n - 1; i++) {
			T taui;
			tlarfg(n - i - 1, A[(i + 1) + L * i], A + std::min(i + 2, n - 1) + L * i, 1, taui);
			e[i] = A[(i + 1) + L * i];
			if (taui != T(0)) {
				A[(i + 1) + L * i] = T(1);
				tsymv('L', n - i - 1, taui, A + (i + 1) + L * (i + 1), lda, A + (i + 1) + L * i, 1,
					T(0), tau + i, 1);
				const T alpha = -(T(1) / T(2)) * taui * tdot(n - i - 1, tau + i, 1, A + (i + 1) + L * i, 1);
				taxpy(n - i - 1, alpha, A + (i + 1) + L * i, 1, tau + i, 1);
				tsyr2('L', n - i - 1, -T(1), A + (i + 1) + L * i, 1, tau + i, 1,
					A + (i + 1) + L * (i + 1), lda);
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
template <typename T>
inline void tlatrd(const char uplo, const int n, const int nb, T* A, const int lda,
	T* e, T* tau, T* W, const int ldw) {
	namespace detail = tlapack_detail;
	if (n <= 0) {
		return;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LW = static_cast<std::size_t>(ldw);
	if (detail::option_is(uplo, 'U')) {
		for (int i = n - 1; i >= n - nb; i--) {
			const int iw = i - n + nb;
			if (i < n - 1) {
				tgemv('N', i + 1, n - 1 - i, -T(1), A + L * (i + 1), lda, W + i + LW * (iw + 1), ldw,
					T(1), A + L * i, 1);
				tgemv('N', i + 1, n - 1 - i, -T(1), W + LW * (iw + 1), ldw, A + i + L * (i + 1), lda,
					T(1), A + L * i, 1);
			}
			if (i > 0) {
				tlarfg(i, A[(i - 1) + L * i], A + L * i, 1, tau[i - 1]);
				e[i - 1] = A[(i - 1) + L * i];
				A[(i - 1) + L * i] = T(1);
				tsymv('U', i, T(1), A, lda, A + L * i, 1, T(0), W + LW * iw, 1);
				if (i < n - 1) {
					tgemv('T', i, n - 1 - i, T(1), W + LW * (iw + 1), ldw, A + L * i, 1,
						T(0), W + (i + 1) + LW * iw, 1);
					tgemv('N', i, n - 1 - i, -T(1), A + L * (i + 1), lda, W + (i + 1) + LW * iw, 1,
						T(1), W + LW * iw, 1);
					tgemv('T', i, n - 1 - i, T(1), A + L * (i + 1), lda, A + L * i, 1,
						T(0), W + (i + 1) + LW * iw, 1);
					tgemv('N', i, n - 1 - i, -T(1), W + LW * (iw + 1), ldw, W + (i + 1) + LW * iw, 1,
						T(1), W + LW * iw, 1);
				}
				tscal(i, tau[i - 1], W + LW * iw, 1);
				const T alpha = -(T(1) / T(2)) * tau[i - 1] * tdot(i, W + LW * iw, 1, A + L * i, 1);
				taxpy(i, alpha, A + L * i, 1, W + LW * iw, 1);
			}
		}
	}
	else {
		for (int i = 0; i < nb; i++) {
			tgemv('N', n - i, i, -T(1), A + i, lda, W + i, ldw, T(1), A + i + L * i, 1);
			tgemv('N', n - i, i, -T(1), W + i, ldw, A + i, lda, T(1), A + i + L * i, 1);
			if (i < n - 1) {
				tlarfg(n - i - 1, A[(i + 1) + L * i], A + std::min(i + 2, n - 1) + L * i, 1, tau[i]);
				e[i] = A[(i + 1) + L * i];
				A[(i + 1) + L * i] = T(1);
				tsymv('L', n - i - 1, T(1), A + (i + 1) + L * (i + 1), lda, A + (i + 1) + L * i, 1,
					T(0), W + (i + 1) + LW * i, 1);
				tgemv('T', n - i - 1, i, T(1), W + (i + 1), ldw, A + (i + 1) + L * i, 1,
					T(0), W + LW * i, 1);
				tgemv('N', n - i - 1, i, -T(1), A + (i + 1), lda, W + LW * i, 1,
					T(1), W + (i + 1) + LW * i, 1);
				tgemv('T', n - i - 1, i, T(1), A + (i + 1), lda, A + (i + 1) + L * i, 1,
					T(0), W + LW * i, 1);
				tgemv('N', n - i - 1, i, -T(1), W + (i + 1), ldw, W + LW * i, 1,
					T(1), W + (i + 1) + LW * i, 1);
				tscal(n - i - 1, tau[i], W + (i + 1) + LW * i, 1);
				const T alpha = -(T(1) / T(2)) * tau[i] * tdot(n - i - 1, W + (i + 1) + LW * i, 1,
					A + (i + 1) + L * i, 1);
				taxpy(n - i - 1, alpha, A + (i + 1) + L * i, 1, W + (i + 1) + LW * i, 1);
			}
		}
	}
}

// 対称行列の三重対角化 (blocked)
template <typename T>
inline int tsytrd(const char uplo, const int n, T* A, const int lda,
	T* d, T* e, T* tau) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsytrd: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsytrd: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
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
		return tsytd2(uplo, n, A, lda, d, e, tau);
	}
	const int ldwork = n;
	std::vector<T> work(static_cast<std::size_t>(ldwork) * nb, T(0));
	if (upper) {
		const int kk = n - ((n - nx + nb - 1) / nb) * nb;
		for (int i = n - nb; i >= kk; i -= nb) {
			tlatrd('U', i + nb, nb, A, lda, e, tau, work.data(), ldwork);
			tsyr2k('U', 'N', i, nb, -T(1), A + L * i, lda, work.data(), ldwork, T(1), A, lda);
			for (int j = i; j < i + nb; j++) {
				A[(j - 1) + L * j] = e[j - 1];
				d[j] = A[j + L * j];
			}
		}
		tsytd2('U', kk, A, lda, d, e, tau);
	}
	else {
		int i = 0;
		for (; i < n - nx; i += nb) {
			tlatrd('L', n - i, nb, A + i + L * i, lda, e + i, tau + i, work.data(), ldwork);
			tsyr2k('L', 'N', n - i - nb, nb, -T(1), A + (i + nb) + L * i, lda,
				work.data() + nb, ldwork, T(1), A + (i + nb) + L * (i + nb), lda);
			for (int j = i; j < i + nb; j++) {
				A[(j + 1) + L * j] = e[j];
				d[j] = A[j + L * j];
			}
		}
		tsytd2('L', n - i, A + i + L * i, lda, d + i, e + i, tau + i);
	}
	return 0;
}

// tsytrd の Q の生成
template <typename T>
inline int torgtr(const char uplo, const int n, T* A, const int lda,
	const T* tau) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("torgtr: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("torgtr: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int j = 0; j < n - 1; j++) {
			for (int i = 0; i < j; i++) {
				A[i + L * j] = A[i + L * (j + 1)];
			}
			A[(n - 1) + L * j] = T(0);
		}
		for (int i = 0; i < n - 1; i++) {
			A[i + L * (n - 1)] = T(0);
		}
		A[(n - 1) + L * (n - 1)] = T(1);
		torgql(n - 1, n - 1, n - 1, A, lda, tau);
	}
	else {
		for (int j = n - 1; j >= 1; j--) {
			A[L * j] = T(0);
			for (int i = j + 1; i < n; i++) {
				A[i + L * j] = A[i + L * (j - 1)];
			}
		}
		A[0] = T(1);
		for (int i = 1; i < n; i++) {
			A[i] = T(0);
		}
		if (n > 1) {
			torgqr(n - 1, n - 1, n - 1, A + 1 + L, lda, tau);
		}
	}
	return 0;
}

// tsytrd の Q の適用
template <typename T>
inline int tormtr(const char side, const char uplo, const char trans, const int m, const int n,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool upper = detail::option_is(uplo, 'U');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!upper && !detail::option_is(uplo, 'L')) ||
		(!detail::option_is(trans, 'N') && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::tlapack_error("tormtr: invalid argument");
	}
	if (m == 0 || n == 0 || nq == 1) {
		return 0;
	}
	const int mi = left ? m - 1 : m;
	const int ni = left ? n : n - 1;
	if (upper) {
		tormql(side, trans, mi, ni, nq - 1, A + static_cast<std::size_t>(lda), lda, tau, C, ldc);
	}
	else {
		T* C1 = C + (left ? 1 : static_cast<std::size_t>(ldc));
		tormqr(side, trans, mi, ni, nq - 1, A + 1, lda, tau, C1, ldc);
	}
	return 0;
}

// 対称三重対角行列の全固有値 (Pal-Walker-Kahan QL/QR, 固有 vector なし)
// d, e は 1-based の Fortran 流 index で記述 (access は [idx-1])
template <typename T>
inline int tsterf(const int n, T* d, T* e) {
	namespace detail = tlapack_detail;
	if (n < 0) {
		detail::tlapack_error("tsterf: invalid n");
	}
	if (n <= 1) {
		return 0;
	}
	int info = 0;
	const T eps = tlamch<T>('E');
	const T eps2 = eps * eps;
	const T safmin = tlamch<T>('S');
	const T safmax = T(1) / safmin;
	const T ssfmax = tlapack_detail::tsqrt(safmax) / T(3);
	const T ssfmin = tlapack_detail::tsqrt(safmin) / eps2;
	const int nmaxit = n * 30;
	int jtot = 0;
	int l1 = 1;
	int l = 0, lsv = 0, lend = 0, lendsv = 0, iscale = 0;
	T anorm = T(0);

L10:
	if (l1 > n) {
		tlasrt('I', n, d);
		return info;
	}
	if (l1 > 1) {
		e[l1 - 2] = T(0);
	}
	{
		int m = n;
		for (int mm = l1; mm <= n - 1; mm++) {
			if (tlapack_detail::tabs(e[mm - 1]) <= (tlapack_detail::tsqrt(tlapack_detail::tabs(d[mm - 1])) * tlapack_detail::tsqrt(tlapack_detail::tabs(d[mm]))) * eps) {
				e[mm - 1] = T(0);
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
	anorm = tlanst('M', lend - l + 1, d + (l - 1), e + (l - 1));
	iscale = 0;
	if (anorm == T(0)) {
		goto L10;
	}
	if (anorm > ssfmax) {
		iscale = 1;
		tlascl('G', 0, 0, anorm, ssfmax, lend - l + 1, 1, d + (l - 1), n);
		tlascl('G', 0, 0, anorm, ssfmax, lend - l, 1, e + (l - 1), n);
	}
	else if (anorm < ssfmin) {
		iscale = 2;
		tlascl('G', 0, 0, anorm, ssfmin, lend - l + 1, 1, d + (l - 1), n);
		tlascl('G', 0, 0, anorm, ssfmin, lend - l, 1, e + (l - 1), n);
	}
	for (int i = l; i <= lend - 1; i++) {
		e[i - 1] = e[i - 1] * e[i - 1];
	}
	if (tlapack_detail::tabs(d[lend - 1]) < tlapack_detail::tabs(d[l - 1])) {
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
					if (tlapack_detail::tabs(e[mm - 1]) <= eps2 * tlapack_detail::tabs(d[mm - 1] * d[mm])) {
						m = mm;
						break;
					}
				}
			}
			if (m < lend) {
				e[m - 1] = T(0);
			}
			T p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l + 1;
				if (l <= lend) {
					goto L50;
				}
				goto L150;
			}
			if (m == l + 1) {
				T rt1, rt2;
				const T rte = tlapack_detail::tsqrt(e[l - 1]);
				tlae2(d[l - 1], rte, d[l], rt1, rt2);
				d[l - 1] = rt1;
				d[l] = rt2;
				e[l - 1] = T(0);
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
			const T rte = tlapack_detail::tsqrt(e[l - 1]);
			T sigma = (d[l] - p) / (T(2) * rte);
			T r = tlapy2(sigma, T(1));
			sigma = p - (rte / (sigma + detail::f_sign(r, sigma)));
			T c = T(1);
			T s = T(0);
			T gamma = d[m - 1] - sigma;
			p = gamma * gamma;
			for (int i = m - 1; i >= l; i--) {
				const T bb = e[i - 1];
				r = p + bb;
				if (i != m - 1) {
					e[i] = s * r;
				}
				const T oldc = c;
				c = p / r;
				s = bb / r;
				const T oldgam = gamma;
				const T alpha = d[i - 1];
				gamma = c * (alpha - sigma) - s * oldgam;
				d[i] = oldgam + (alpha - gamma);
				if (c != T(0)) {
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
				if (tlapack_detail::tabs(e[mm - 2]) <= eps2 * tlapack_detail::tabs(d[mm - 1] * d[mm - 2])) {
					m = mm;
					break;
				}
			}
			if (m > lend) {
				e[m - 2] = T(0);
			}
			T p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l - 1;
				if (l >= lend) {
					goto L100;
				}
				goto L150;
			}
			if (m == l - 1) {
				T rt1, rt2;
				const T rte = tlapack_detail::tsqrt(e[l - 2]);
				tlae2(d[l - 1], rte, d[l - 2], rt1, rt2);
				d[l - 1] = rt1;
				d[l - 2] = rt2;
				e[l - 2] = T(0);
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
			const T rte = tlapack_detail::tsqrt(e[l - 2]);
			T sigma = (d[l - 2] - p) / (T(2) * rte);
			T r = tlapy2(sigma, T(1));
			sigma = p - (rte / (sigma + detail::f_sign(r, sigma)));
			T c = T(1);
			T s = T(0);
			T gamma = d[m - 1] - sigma;
			p = gamma * gamma;
			for (int i = m; i <= l - 1; i++) {
				const T bb = e[i - 1];
				r = p + bb;
				if (i != m) {
					e[i - 2] = s * r;
				}
				const T oldc = c;
				c = p / r;
				s = bb / r;
				const T oldgam = gamma;
				const T alpha = d[i];
				gamma = c * (alpha - sigma) - s * oldgam;
				d[i - 1] = oldgam + (alpha - gamma);
				if (c != T(0)) {
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
		tlascl('G', 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n);
	}
	if (iscale == 2) {
		tlascl('G', 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n);
	}
	if (jtot < nmaxit) {
		goto L10;
	}
	for (int i = 1; i <= n - 1; i++) {
		if (e[i - 1] != T(0)) {
			info++;
		}
	}
	return info;
}

// 対称三重対角行列の固有値・固有 vector (implicit QL/QR)
// compz: 'N' 固有値のみ, 'V' Z に三重対角化の Q を渡して元行列の固有 vector,
//        'I' Z を単位行列で初期化して三重対角行列の固有 vector
template <typename T>
inline int tsteqr(const char compz, const int n, T* d, T* e,
	T* Z, const int ldz) {
	namespace detail = tlapack_detail;
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
		detail::tlapack_error("tsteqr: invalid compz");
		return -1;
	}
	if (n < 0 || ldz < 1 || (icompz > 0 && ldz < std::max(1, n))) {
		detail::tlapack_error("tsteqr: invalid n/ldz");
	}
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		if (icompz == 2) {
			Z[0] = T(1);
		}
		return 0;
	}
	const std::size_t LZ = static_cast<std::size_t>(ldz);
	int info = 0;
	const T eps = tlamch<T>('E');
	const T eps2 = eps * eps;
	const T safmin = tlamch<T>('S');
	const T safmax = T(1) / safmin;
	const T ssfmax = tlapack_detail::tsqrt(safmax) / T(3);
	const T ssfmin = tlapack_detail::tsqrt(safmin) / eps2;
	if (icompz == 2) {
		tlaset('F', n, n, T(0), T(1), Z, ldz);
	}
	std::vector<T> work(icompz > 0 ? static_cast<std::size_t>(2 * n - 2) : 1, T(0));
	const int nmaxit = n * 30;
	int jtot = 0;
	int l1 = 1;
	int l = 0, lsv = 0, lend = 0, lendsv = 0, iscale = 0;
	T anorm = T(0);

L10:
	if (l1 > n) {
		goto L160;
	}
	if (l1 > 1) {
		e[l1 - 2] = T(0);
	}
	{
		int m = n;
		if (l1 <= n - 1) {
			for (int mm = l1; mm <= n - 1; mm++) {
				const T tst = tlapack_detail::tabs(e[mm - 1]);
				if (tst == T(0)) {
					m = mm;
					break;
				}
				if (tst <= (tlapack_detail::tsqrt(tlapack_detail::tabs(d[mm - 1])) * tlapack_detail::tsqrt(tlapack_detail::tabs(d[mm]))) * eps) {
					e[mm - 1] = T(0);
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
	anorm = tlanst('M', lend - l + 1, d + (l - 1), e + (l - 1));
	iscale = 0;
	if (anorm == T(0)) {
		goto L10;
	}
	if (anorm > ssfmax) {
		iscale = 1;
		tlascl('G', 0, 0, anorm, ssfmax, lend - l + 1, 1, d + (l - 1), n);
		tlascl('G', 0, 0, anorm, ssfmax, lend - l, 1, e + (l - 1), n);
	}
	else if (anorm < ssfmin) {
		iscale = 2;
		tlascl('G', 0, 0, anorm, ssfmin, lend - l + 1, 1, d + (l - 1), n);
		tlascl('G', 0, 0, anorm, ssfmin, lend - l, 1, e + (l - 1), n);
	}
	if (tlapack_detail::tabs(d[lend - 1]) < tlapack_detail::tabs(d[l - 1])) {
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
					const T tst = tlapack_detail::tabs(e[mm - 1]) * tlapack_detail::tabs(e[mm - 1]);
					if (tst <= (eps2 * tlapack_detail::tabs(d[mm - 1])) * tlapack_detail::tabs(d[mm]) + safmin) {
						m = mm;
						break;
					}
				}
			}
			if (m < lend) {
				e[m - 1] = T(0);
			}
			T p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l + 1;
				if (l <= lend) {
					goto L40;
				}
				goto L140;
			}
			if (m == l + 1) {
				T rt1, rt2, c, s;
				if (icompz > 0) {
					tlaev2(d[l - 1], e[l - 1], d[l], rt1, rt2, c, s);
					work[l - 1] = c;
					work[n - 1 + l - 1] = s;
					tlasr('R', 'V', 'B', n, 2, work.data() + (l - 1), work.data() + (n - 1 + l - 1),
						Z + LZ * (l - 1), ldz);
				}
				else {
					tlae2(d[l - 1], e[l - 1], d[l], rt1, rt2);
				}
				d[l - 1] = rt1;
				d[l] = rt2;
				e[l - 1] = T(0);
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
			T g = (d[l] - p) / (T(2) * e[l - 1]);
			T r = tlapy2(g, T(1));
			g = d[m - 1] - p + (e[l - 1] / (g + detail::f_sign(r, g)));
			T s = T(1);
			T c = T(1);
			p = T(0);
			for (int i = m - 1; i >= l; i--) {
				T f = s * e[i - 1];
				const T b = c * e[i - 1];
				tlartg(g, f, c, s, r);
				if (i != m - 1) {
					e[i] = r;
				}
				g = d[i] - p;
				r = (d[i - 1] - g) * s + T(2) * c * b;
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
				tlasr('R', 'V', 'B', n, mm, work.data() + (l - 1), work.data() + (n - 1 + l - 1),
					Z + LZ * (l - 1), ldz);
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
					const T tst = tlapack_detail::tabs(e[mm - 2]) * tlapack_detail::tabs(e[mm - 2]);
					if (tst <= (eps2 * tlapack_detail::tabs(d[mm - 1])) * tlapack_detail::tabs(d[mm - 2]) + safmin) {
						m = mm;
						break;
					}
				}
			}
			if (m > lend) {
				e[m - 2] = T(0);
			}
			T p = d[l - 1];
			if (m == l) {
				d[l - 1] = p;
				l = l - 1;
				if (l >= lend) {
					goto L90;
				}
				goto L140;
			}
			if (m == l - 1) {
				T rt1, rt2, c, s;
				if (icompz > 0) {
					tlaev2(d[l - 2], e[l - 2], d[l - 1], rt1, rt2, c, s);
					work[m - 1] = c;
					work[n - 1 + m - 1] = s;
					tlasr('R', 'V', 'F', n, 2, work.data() + (m - 1), work.data() + (n - 1 + m - 1),
						Z + LZ * (l - 2), ldz);
				}
				else {
					tlae2(d[l - 2], e[l - 2], d[l - 1], rt1, rt2);
				}
				d[l - 2] = rt1;
				d[l - 1] = rt2;
				e[l - 2] = T(0);
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
			T g = (d[l - 2] - p) / (T(2) * e[l - 2]);
			T r = tlapy2(g, T(1));
			g = d[m - 1] - p + (e[l - 2] / (g + detail::f_sign(r, g)));
			T s = T(1);
			T c = T(1);
			p = T(0);
			for (int i = m; i <= l - 1; i++) {
				T f = s * e[i - 1];
				const T b = c * e[i - 1];
				tlartg(g, f, c, s, r);
				if (i != m) {
					e[i - 2] = r;
				}
				g = d[i - 1] - p;
				r = (d[i] - g) * s + T(2) * c * b;
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
				tlasr('R', 'V', 'F', n, mm, work.data() + (m - 1), work.data() + (n - 1 + m - 1),
					Z + LZ * (m - 1), ldz);
			}
			d[l - 1] = d[l - 1] - p;
			e[l - 2] = g;
			goto L90;
		}
	}

L140:
	if (iscale == 1) {
		tlascl('G', 0, 0, ssfmax, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n);
		tlascl('G', 0, 0, ssfmax, anorm, lendsv - lsv, 1, e + (lsv - 1), n);
	}
	else if (iscale == 2) {
		tlascl('G', 0, 0, ssfmin, anorm, lendsv - lsv + 1, 1, d + (lsv - 1), n);
		tlascl('G', 0, 0, ssfmin, anorm, lendsv - lsv, 1, e + (lsv - 1), n);
	}
	if (jtot < nmaxit) {
		goto L10;
	}
	for (int i = 1; i <= n - 1; i++) {
		if (e[i - 1] != T(0)) {
			info++;
		}
	}
	return info;

L160:
	// 固有値を昇順に sort (vector 付きは選択 sort で Z の列も入れ替え)
	if (icompz == 0) {
		tlasrt('I', n, d);
	}
	else {
		for (int ii = 2; ii <= n; ii++) {
			const int i = ii - 1;
			int k = i;
			T p = d[i - 1];
			for (int j = ii; j <= n; j++) {
				if (d[j - 1] < p) {
					k = j;
					p = d[j - 1];
				}
			}
			if (k != i) {
				d[k - 1] = d[i - 1];
				d[i - 1] = p;
				tswap(n, Z + LZ * (i - 1), 1, Z + LZ * (k - 1), 1);
			}
		}
	}
	return info;
}

// 対称行列の固有値・固有 vector (driver)
// jobz: 'N' 固有値のみ, 'V' 固有 vector も計算 (A に上書き)
// w: 固有値 (昇順, n 個)
template <typename T>
inline int tsyev(const char jobz, const char uplo, const int n, T* A, const int lda,
	T* w) {
	namespace detail = tlapack_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	const bool lower = detail::option_is(uplo, 'L');
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::tlapack_error("tsyev: invalid jobz");
	}
	if (!lower && !detail::option_is(uplo, 'U')) {
		detail::tlapack_error("tsyev: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsyev: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		w[0] = A[0];
		if (wantz) {
			A[0] = T(1);
		}
		return 0;
	}
	const T safmin = tlamch<T>('S');
	const T eps = tlamch<T>('P');
	const T smlnum = safmin / eps;
	const T bignum = T(1) / smlnum;
	const T rmin = tlapack_detail::tsqrt(smlnum);
	const T rmax = tlapack_detail::tsqrt(bignum);
	const T anrm = tlansy('M', uplo, n, A, lda);
	int iscale = 0;
	T sigma = T(1);
	if (anrm > T(0) && anrm < rmin) {
		iscale = 1;
		sigma = rmin / anrm;
	}
	else if (anrm > rmax) {
		iscale = 1;
		sigma = rmax / anrm;
	}
	if (iscale == 1) {
		tlascl(uplo, 0, 0, T(1), sigma, n, n, A, lda);
	}
	std::vector<T> e(static_cast<std::size_t>(n), T(0));
	std::vector<T> tau(static_cast<std::size_t>(n), T(0));
	tsytrd(uplo, n, A, lda, w, e.data(), tau.data());
	int info;
	if (!wantz) {
		info = tsterf(n, w, e.data());
	}
	else {
		torgtr(uplo, n, A, lda, tau.data());
		info = tsteqr(jobz, n, w, e.data(), A, lda);
	}
	if (iscale == 1) {
		const int imax = (info == 0) ? n : info - 1;
		tscal(imax, T(1) / sigma, w, 1);
	}
	return info;
}

#endif // TLAPACK_TLAPACK_EIG_HPP
