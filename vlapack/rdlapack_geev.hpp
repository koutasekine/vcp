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

#ifndef VBLAS_RDLAPACK_GEEV_HPP
#define VBLAS_RDLAPACK_GEEV_HPP

// 非対称固有値問題 (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dgebal/dgebak/dgehd2/dlahr2/dgehrd/dorghr/
// dlahqr/dlaln2/dtrevc/dgeev を移植．
//
// - ilo/ihi と scale の swap 先 index は LAPACK と同じ 1-based (注意)．
// - rdhseqr は dlahqr (double-shift QR) を全 size に使う簡略版
//   (reference の dlaqr0 系 multishift は移植していない．結果の意味は同じ)．
// - rdtrevc は howmny = 'A' / 'B' のみ (SELECT 指定は未対応)．

#include "rdlapack_svd.hpp"

#pragma STDC FENV_ACCESS ON

// 行列の balancing (permute + scale)．ilo/ihi は 1-based で返る
inline int rdgebal(const char job, const int n, double* A, const int lda,
	int& ilo, int& ihi, double* scale, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (!detail::option_is(job, 'N') && !detail::option_is(job, 'P') &&
		!detail::option_is(job, 'S') && !detail::option_is(job, 'B')) {
		detail::rdlapack_error("rdgebal: invalid job");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdgebal: invalid n/lda");
	}
	if (n == 0) {
		ilo = 1;
		ihi = 0;
		return 0;
	}
	if (detail::option_is(job, 'N')) {
		for (int i = 0; i < n; i++) {
			scale[i] = 1.0;
		}
		ilo = 1;
		ihi = n;
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	int k = 1;
	int l = n;
	if (!detail::option_is(job, 'S')) {
		// permutation: 行/列の zero 構造で対角 block へ分離
		bool noconv = true;
		while (noconv) {
			noconv = false;
			for (int i = l; i >= 1; i--) {
				bool canswap = true;
				for (int j = 1; j <= l; j++) {
					if (i != j && A[(i - 1) + L * (j - 1)] != 0.0) {
						canswap = false;
						break;
					}
				}
				if (canswap) {
					scale[l - 1] = static_cast<double>(i);
					if (i != l) {
						dswap(l, A + L * (i - 1), 1, A + L * (l - 1), 1);
						dswap(n - k + 1, A + (i - 1) + L * (k - 1), lda, A + (l - 1) + L * (k - 1), lda);
					}
					noconv = true;
					if (l == 1) {
						ilo = 1;
						ihi = 1;
						return 0;
					}
					l = l - 1;
				}
			}
		}
		noconv = true;
		while (noconv) {
			noconv = false;
			for (int j = k; j <= l; j++) {
				bool canswap = true;
				for (int i = k; i <= l; i++) {
					if (i != j && A[(i - 1) + L * (j - 1)] != 0.0) {
						canswap = false;
						break;
					}
				}
				if (canswap) {
					scale[k - 1] = static_cast<double>(j);
					if (j != k) {
						dswap(l, A + L * (j - 1), 1, A + L * (k - 1), 1);
						dswap(n - k + 1, A + (j - 1) + L * (k - 1), lda, A + (k - 1) + L * (k - 1), lda);
					}
					noconv = true;
					k = k + 1;
				}
			}
		}
	}
	for (int i = k; i <= l; i++) {
		scale[i - 1] = 1.0;
	}
	if (detail::option_is(job, 'P')) {
		ilo = k;
		ihi = l;
		return 0;
	}
	// scaling (基数 2 の冪で行/列 norm を揃える)
	const double sclfac = 2.0;
	const double factor = 0.95;
	const double sfmin1 = DBL_MIN / DBL_EPSILON;
	const double sfmax1 = 1.0 / sfmin1;
	const double sfmin2 = sfmin1 * sclfac;
	const double sfmax2 = 1.0 / sfmin2;
	bool noconv = true;
	while (noconv) {
		noconv = false;
		for (int i = k; i <= l; i++) {
			double c = rdnrm2(l - k + 1, A + (k - 1) + L * (i - 1), 1, rm);
			double r = rdnrm2(l - k + 1, A + (i - 1) + L * (k - 1), lda, rm);
			const int ica = idamax(l, A + L * (i - 1), 1) + 1;
			double ca = std::fabs(A[(ica - 1) + L * (i - 1)]);
			const int ira = idamax(n - k + 1, A + (i - 1) + L * (k - 1), lda) + 1;
			double ra = std::fabs(A[(i - 1) + L * (ira + k - 2)]);
			if (c == 0.0 || r == 0.0) {
				continue;
			}
			if (std::isnan(c + ca + r + ra)) {
				detail::rdlapack_error("rdgebal: NaN in matrix");
			}
			double g = r / sclfac;
			double f = 1.0;
			const double s = c + r;
			while (c < g && std::max(std::max(f, c), ca) < sfmax2 &&
				std::min(std::min(r, g), ra) > sfmin2) {
				f = f * sclfac;
				c = c * sclfac;
				ca = ca * sclfac;
				r = r / sclfac;
				g = g / sclfac;
				ra = ra / sclfac;
			}
			g = c / sclfac;
			while (g >= r && std::max(r, ra) < sfmax2 &&
				std::min(std::min(std::min(f, c), g), ca) > sfmin2) {
				f = f / sclfac;
				c = c / sclfac;
				g = g / sclfac;
				ca = ca / sclfac;
				r = r * sclfac;
				ra = ra * sclfac;
			}
			if ((c + r) >= factor * s) {
				continue;
			}
			if (f < 1.0 && scale[i - 1] < 1.0) {
				if (f * scale[i - 1] <= sfmin1) {
					continue;
				}
			}
			if (f > 1.0 && scale[i - 1] > 1.0) {
				if (scale[i - 1] >= sfmax1 / f) {
					continue;
				}
			}
			g = 1.0 / f;
			scale[i - 1] = scale[i - 1] * f;
			noconv = true;
			rdscal(n - k + 1, g, A + (i - 1) + L * (k - 1), lda, rm);
			rdscal(l, f, A + L * (i - 1), 1, rm);
		}
	}
	ilo = k;
	ihi = l;
	return 0;
}

// balancing の逆変換を固有 vector へ適用 (ilo/ihi は 1-based)
inline int rdgebak(const char job, const char side, const int n, const int ilo, const int ihi,
	const double* scale, const int m, double* V, const int ldv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool rightv = detail::option_is(side, 'R');
	const bool leftv = detail::option_is(side, 'L');
	if ((!detail::option_is(job, 'N') && !detail::option_is(job, 'P') &&
		!detail::option_is(job, 'S') && !detail::option_is(job, 'B')) || (!rightv && !leftv)) {
		detail::rdlapack_error("rdgebak: invalid job/side");
	}
	if (n < 0 || m < 0 || ldv < std::max(1, n)) {
		detail::rdlapack_error("rdgebak: invalid argument");
	}
	if (n == 0 || m == 0 || detail::option_is(job, 'N')) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t LV = static_cast<std::size_t>(ldv);
	if (ilo != ihi && (detail::option_is(job, 'S') || detail::option_is(job, 'B'))) {
		if (rightv) {
			for (int i = ilo; i <= ihi; i++) {
				rdscal(m, scale[i - 1], V + (i - 1), ldv, rounding_mode);
			}
		}
		if (leftv) {
			for (int i = ilo; i <= ihi; i++) {
				rdscal(m, 1.0 / scale[i - 1], V + (i - 1), ldv, rounding_mode);
			}
		}
	}
	if (detail::option_is(job, 'P') || detail::option_is(job, 'B')) {
		for (int pass = 0; pass < 2; pass++) {
			if ((pass == 0 && !rightv) || (pass == 1 && !leftv)) {
				continue;
			}
			for (int ii = 1; ii <= n; ii++) {
				int i = ii;
				if (i >= ilo && i <= ihi) {
					continue;
				}
				if (i < ilo) {
					i = ilo - ii;
				}
				const int kk = static_cast<int>(scale[i - 1]);
				if (kk == i) {
					continue;
				}
				dswap(m, V + (i - 1), ldv, V + (kk - 1), ldv);
			}
		}
	}
	return 0;
}

// Hessenberg 化 (unblocked, ilo/ihi は 1-based)
inline int rdgehd2(const int n, const int ilo, const int ihi, double* A, const int lda,
	double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		lda < std::max(1, n)) {
		detail::rdlapack_error("rdgehd2: invalid argument");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	for (int i = ilo; i <= ihi - 1; i++) {
		rdlarfg(ihi - i, A[i + L * (i - 1)], A + std::min(i + 1, n - 1) + L * (i - 1), 1,
			tau[i - 1], rounding_mode);
		rdlarf1f('R', ihi, ihi - i, A + i + L * (i - 1), 1, tau[i - 1], A + L * i, lda, rounding_mode);
		rdlarf1f('L', ihi - i, n - i, A + i + L * (i - 1), 1, tau[i - 1], A + i + L * i, lda, rounding_mode);
	}
	return 0;
}

// Hessenberg 化の先頭 nb 段 (blocked 用 panel, k は 1-based)
inline void rdlahr2(const int n, const int k, const int nb, double* A, const int lda,
	double* tau, double* T, const int ldt, double* Y, const int ldy, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n <= 1) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LT = static_cast<std::size_t>(ldt);
	const std::size_t LY = static_cast<std::size_t>(ldy);
	double ei = 0.0;
	for (int i = 1; i <= nb; i++) {
		if (i > 1) {
			rdgemv('N', n - k, i - 1, -1.0, Y + k, ldy, A + (k + i - 2), lda,
				1.0, A + k + L * (i - 1), 1, rm);
			dcopy(i - 1, A + k + L * (i - 1), 1, T + LT * (nb - 1), 1);
			rdtrmv('L', 'T', 'U', i - 1, A + k, lda, T + LT * (nb - 1), 1, rm);
			rdgemv('T', n - k - i + 1, i - 1, 1.0, A + (k + i - 1), lda,
				A + (k + i - 1) + L * (i - 1), 1, 1.0, T + LT * (nb - 1), 1, rm);
			rdtrmv('U', 'T', 'N', i - 1, T, ldt, T + LT * (nb - 1), 1, rm);
			rdgemv('N', n - k - i + 1, i - 1, -1.0, A + (k + i - 1), lda, T + LT * (nb - 1), 1,
				1.0, A + (k + i - 1) + L * (i - 1), 1, rm);
			rdtrmv('L', 'N', 'U', i - 1, A + k, lda, T + LT * (nb - 1), 1, rm);
			rdaxpy(i - 1, -1.0, T + LT * (nb - 1), 1, A + k + L * (i - 1), 1, rm);
			A[(k + i - 2) + L * (i - 2)] = ei;
		}
		rdlarfg(n - k - i + 1, A[(k + i - 1) + L * (i - 1)],
			A + std::min(k + i, n - 1) + L * (i - 1), 1, tau[i - 1], rm);
		ei = A[(k + i - 1) + L * (i - 1)];
		A[(k + i - 1) + L * (i - 1)] = 1.0;
		rdgemv('N', n - k, n - k - i + 1, 1.0, A + k + L * i, lda,
			A + (k + i - 1) + L * (i - 1), 1, 0.0, Y + k + LY * (i - 1), 1, rm);
		rdgemv('T', n - k - i + 1, i - 1, 1.0, A + (k + i - 1), lda,
			A + (k + i - 1) + L * (i - 1), 1, 0.0, T + LT * (i - 1), 1, rm);
		rdgemv('N', n - k, i - 1, -1.0, Y + k, ldy, T + LT * (i - 1), 1,
			1.0, Y + k + LY * (i - 1), 1, rm);
		rdscal(n - k, tau[i - 1], Y + k + LY * (i - 1), 1, rm);
		rdscal(i - 1, -tau[i - 1], T + LT * (i - 1), 1, rm);
		rdtrmv('U', 'N', 'N', i - 1, T, ldt, T + LT * (i - 1), 1, rm);
		T[(i - 1) + LT * (i - 1)] = tau[i - 1];
	}
	A[(k + nb - 1) + L * (nb - 1)] = ei;
	dlacpy('A', k, nb, A + L, lda, Y, ldy);
	rdtrmm('R', 'L', 'N', 'U', k, nb, 1.0, A + k, lda, Y, ldy, rm);
	if (n > k + nb) {
		rdgemm('N', 'N', k, nb, n - k - nb, 1.0, A + L * (nb + 1), lda, A + (k + nb), lda,
			1.0, Y, ldy, rm);
	}
	rdtrmm('R', 'U', 'N', 'N', k, nb, 1.0, T, ldt, Y, ldy, rm);
}

// Hessenberg 化 (blocked, ilo/ihi は 1-based)
inline int rdgehrd(const int n, const int ilo, const int ihi, double* A, const int lda,
	double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		lda < std::max(1, n)) {
		detail::rdlapack_error("rdgehrd: invalid argument");
	}
	for (int i = 1; i <= ilo - 1; i++) {
		tau[i - 1] = 0.0;
	}
	for (int i = std::max(1, ihi); i <= n - 1; i++) {
		tau[i - 1] = 0.0;
	}
	const int nh = ihi - ilo + 1;
	if (nh <= 1) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = std::min(64, detail::ilaenv(1, "GEHRD", n, ilo, ihi, -1));
	int i = ilo;
	if (nb >= 2 && nb < nh) {
		const int nx = std::max(nb, detail::ilaenv(3, "GEHRD", n, ilo, ihi, -1));
		if (nx < nh) {
			const int ldwork = n;
			std::vector<double> T(static_cast<std::size_t>(nb) * nb);
			std::vector<double> work(static_cast<std::size_t>(ldwork) * nb);
			for (; i <= ihi - 1 - nx; i += nb) {
				const int ib = std::min(nb, ihi - i);
				rdlahr2(ihi, i, ib, A + L * (i - 1), lda, tau + (i - 1), T.data(), nb,
					work.data(), ldwork, rounding_mode);
				const double ei = A[(i + ib - 1) + L * (i + ib - 2)];
				A[(i + ib - 1) + L * (i + ib - 2)] = 1.0;
				rdgemm('N', 'T', ihi, ihi - i - ib + 1, ib, -1.0, work.data(), ldwork,
					A + (i + ib - 1) + L * (i - 1), lda, 1.0, A + L * (i + ib - 1), lda, rounding_mode);
				A[(i + ib - 1) + L * (i + ib - 2)] = ei;
				rdtrmm('R', 'L', 'T', 'U', i, ib - 1, 1.0, A + i + L * (i - 1), lda,
					work.data(), ldwork, rounding_mode);
				for (int j = 0; j <= ib - 2; j++) {
					rdaxpy(i, -1.0, work.data() + static_cast<std::size_t>(ldwork) * j, 1,
						A + L * (i + j), 1, rounding_mode);
				}
				rdlarfb('L', 'T', 'F', 'C', ihi - i, n - i - ib + 1, ib, A + i + L * (i - 1), lda,
					T.data(), nb, A + i + L * (i + ib - 1), lda, rounding_mode);
			}
		}
	}
	rdgehd2(n, i, ihi, A, lda, tau, rounding_mode);
	return 0;
}

// rdgehrd の Q の生成 (ilo/ihi は 1-based)
inline int rdorghr(const int n, const int ilo, const int ihi, double* A, const int lda,
	const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		lda < std::max(1, n)) {
		detail::rdlapack_error("rdorghr: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nh = ihi - ilo;
	// 列を 1 つ右へ shift して reflector を配置し，外側は単位行列にする
	for (int j = ihi; j >= ilo + 1; j--) {
		for (int i = 1; i <= j - 1; i++) {
			A[(i - 1) + L * (j - 1)] = 0.0;
		}
		for (int i = j + 1; i <= ihi; i++) {
			A[(i - 1) + L * (j - 1)] = A[(i - 1) + L * (j - 2)];
		}
		for (int i = ihi + 1; i <= n; i++) {
			A[(i - 1) + L * (j - 1)] = 0.0;
		}
	}
	for (int j = 1; j <= ilo; j++) {
		for (int i = 1; i <= n; i++) {
			A[(i - 1) + L * (j - 1)] = 0.0;
		}
		A[(j - 1) + L * (j - 1)] = 1.0;
	}
	for (int j = ihi + 1; j <= n; j++) {
		for (int i = 1; i <= n; i++) {
			A[(i - 1) + L * (j - 1)] = 0.0;
		}
		A[(j - 1) + L * (j - 1)] = 1.0;
	}
	if (nh > 0) {
		rdorgqr(nh, nh, nh, A + ilo + L * ilo, lda, tau + (ilo - 1), rounding_mode);
	}
	return 0;
}

// Hessenberg 行列の double-shift QR (小規模用 algorithm を全 size に使用)
// ilo/ihi/iloz/ihiz は 1-based．返り値 info > 0 は非収束 (LAPACK と同じ意味)
inline int rdlahqr(const bool wantt, const bool wantz, const int n, const int ilo, const int ihi,
	double* H, const int ldh, double* wr, double* wi, const int iloz, const int ihiz,
	double* Z, const int ldz, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t LH = static_cast<std::size_t>(ldh);
	const std::size_t LZ = static_cast<std::size_t>(ldz);
	if (ilo == ihi) {
		wr[ilo - 1] = H[(ilo - 1) + LH * (ilo - 1)];
		wi[ilo - 1] = 0.0;
		return 0;
	}
	for (int j = ilo; j <= ihi - 3; j++) {
		H[(j + 1) + LH * (j - 1)] = 0.0;
		H[(j + 2) + LH * (j - 1)] = 0.0;
	}
	if (ilo <= ihi - 2) {
		H[(ihi - 1) + LH * (ihi - 3)] = 0.0;
	}
	const int nh = ihi - ilo + 1;
	const int nz = ihiz - iloz + 1;
	const double safmin = DBL_MIN;
	const double ulp = DBL_EPSILON;
	const double smlnum = safmin * (static_cast<double>(nh) / ulp);
	int i1 = 0, i2 = 0;
	if (wantt) {
		i1 = 1;
		i2 = n;
	}
	const int itmax = 30 * std::max(10, nh);
	int kdefl = 0;
	int i = ihi;
	double v[3];

L20:
	{
		int l = ilo;
		if (i < ilo) {
			return 0;
		}
		int its;
		int k = 0;
		int m = 0;
		for (its = 0; its <= itmax; its++) {
			// 副対角の無視できる要素を探す
			bool deflate = false;
			for (k = i; k >= l + 1; k--) {
				if (std::fabs(H[(k - 1) + LH * (k - 2)]) <= smlnum) {
					deflate = true;
					break;
				}
				double tst = std::fabs(H[(k - 2) + LH * (k - 2)]) + std::fabs(H[(k - 1) + LH * (k - 1)]);
				if (tst == 0.0) {
					if (k - 2 >= ilo) {
						tst = tst + std::fabs(H[(k - 2) + LH * (k - 3)]);
					}
					if (k + 1 <= ihi) {
						tst = tst + std::fabs(H[k + LH * (k - 1)]);
					}
				}
				if (std::fabs(H[(k - 1) + LH * (k - 2)]) <= ulp * tst) {
					const double ab = std::max(std::fabs(H[(k - 1) + LH * (k - 2)]),
						std::fabs(H[(k - 2) + LH * (k - 1)]));
					const double ba = std::min(std::fabs(H[(k - 1) + LH * (k - 2)]),
						std::fabs(H[(k - 2) + LH * (k - 1)]));
					const double aa = std::max(std::fabs(H[(k - 1) + LH * (k - 1)]),
						std::fabs(H[(k - 2) + LH * (k - 2)] - H[(k - 1) + LH * (k - 1)]));
					const double bb = std::min(std::fabs(H[(k - 1) + LH * (k - 1)]),
						std::fabs(H[(k - 2) + LH * (k - 2)] - H[(k - 1) + LH * (k - 1)]));
					const double s = aa + ab;
					if (ba * (ab / s) <= std::max(smlnum, ulp * (bb * (aa / s)))) {
						deflate = true;
						break;
					}
				}
			}
			if (!deflate) {
				k = l;
			}
			l = k;
			if (l > ilo) {
				H[(l - 1) + LH * (l - 2)] = 0.0;
			}
			if (l >= i - 1) {
				goto L150;
			}
			kdefl++;
			if (!wantt) {
				i1 = l;
				i2 = i;
			}
			// shift の決定
			double h11, h12, h21, h22;
			if (kdefl % 20 == 0) {
				// exceptional shift
				const double s = std::fabs(H[(i - 1) + LH * (i - 2)]) + std::fabs(H[(i - 2) + LH * (i - 3)]);
				h11 = 0.75 * s + H[(i - 1) + LH * (i - 1)];
				h12 = -0.4375 * s;
				h21 = s;
				h22 = h11;
			}
			else if (kdefl % 10 == 0) {
				const double s = std::fabs(H[l + LH * (l - 1)]) + std::fabs(H[(l + 1) + LH * l]);
				h11 = 0.75 * s + H[(l - 1) + LH * (l - 1)];
				h12 = -0.4375 * s;
				h21 = s;
				h22 = h11;
			}
			else {
				h11 = H[(i - 2) + LH * (i - 2)];
				h21 = H[(i - 1) + LH * (i - 2)];
				h12 = H[(i - 2) + LH * (i - 1)];
				h22 = H[(i - 1) + LH * (i - 1)];
			}
			double rt1r, rt1i, rt2r, rt2i;
			{
				const double s = std::fabs(h11) + std::fabs(h12) + std::fabs(h21) + std::fabs(h22);
				if (s == 0.0) {
					rt1r = rt1i = rt2r = rt2i = 0.0;
				}
				else {
					h11 = h11 / s;
					h21 = h21 / s;
					h12 = h12 / s;
					h22 = h22 / s;
					const double tr = (h11 + h22) / 2.0;
					const double det = (h11 - tr) * (h22 - tr) - h12 * h21;
					const double rtdisc = std::sqrt(std::fabs(det));
					if (det >= 0.0) {
						rt1r = tr * s;
						rt2r = rt1r;
						rt1i = rtdisc * s;
						rt2i = -rt1i;
					}
					else {
						rt1r = tr + rtdisc;
						rt2r = tr - rtdisc;
						if (std::fabs(rt1r - h22) <= std::fabs(rt2r - h22)) {
							rt1r = rt1r * s;
							rt2r = rt1r;
						}
						else {
							rt2r = rt2r * s;
							rt1r = rt2r;
						}
						rt1i = 0.0;
						rt2i = 0.0;
					}
				}
			}
			// 最初の列を決める位置 m を探す
			for (m = i - 2; m >= l; m--) {
				double h21s = H[m + LH * (m - 1)];
				double s = std::fabs(H[(m - 1) + LH * (m - 1)] - rt2r) + std::fabs(rt2i) + std::fabs(h21s);
				h21s = H[m + LH * (m - 1)] / s;
				v[0] = h21s * H[(m - 1) + LH * m] + (H[(m - 1) + LH * (m - 1)] - rt1r) *
					((H[(m - 1) + LH * (m - 1)] - rt2r) / s) - rt1i * (rt2i / s);
				v[1] = h21s * (H[(m - 1) + LH * (m - 1)] + H[m + LH * m] - rt1r - rt2r);
				v[2] = h21s * H[(m + 1) + LH * m];
				s = std::fabs(v[0]) + std::fabs(v[1]) + std::fabs(v[2]);
				v[0] = v[0] / s;
				v[1] = v[1] / s;
				v[2] = v[2] / s;
				if (m == l) {
					break;
				}
				if (std::fabs(H[(m - 1) + LH * (m - 2)]) * (std::fabs(v[1]) + std::fabs(v[2])) <=
					ulp * std::fabs(v[0]) * (std::fabs(H[(m - 2) + LH * (m - 2)]) +
						std::fabs(H[(m - 1) + LH * (m - 1)]) + std::fabs(H[m + LH * m]))) {
					break;
				}
			}
			// double-shift QR sweep
			for (k = m; k <= i - 1; k++) {
				const int nr = std::min(3, i - k + 1);
				double t1;
				if (k > m) {
					dcopy(nr, H + (k - 1) + LH * (k - 2), 1, v, 1);
				}
				rdlarfg(nr, v[0], v + 1, 1, t1, rm);
				if (k > m) {
					H[(k - 1) + LH * (k - 2)] = v[0];
					H[k + LH * (k - 2)] = 0.0;
					if (k < i - 1) {
						H[(k + 1) + LH * (k - 2)] = 0.0;
					}
				}
				else if (m > l) {
					H[(k - 1) + LH * (k - 2)] = H[(k - 1) + LH * (k - 2)] * (1.0 - t1);
				}
				const double v2 = v[1];
				const double t2 = t1 * v2;
				if (nr == 3) {
					const double v3 = v[2];
					const double t3 = t1 * v3;
					for (int j = k; j <= i2; j++) {
						const double sum = H[(k - 1) + LH * (j - 1)] + v2 * H[k + LH * (j - 1)] +
							v3 * H[(k + 1) + LH * (j - 1)];
						H[(k - 1) + LH * (j - 1)] -= sum * t1;
						H[k + LH * (j - 1)] -= sum * t2;
						H[(k + 1) + LH * (j - 1)] -= sum * t3;
					}
					for (int j = i1; j <= std::min(k + 3, i); j++) {
						const double sum = H[(j - 1) + LH * (k - 1)] + v2 * H[(j - 1) + LH * k] +
							v3 * H[(j - 1) + LH * (k + 1)];
						H[(j - 1) + LH * (k - 1)] -= sum * t1;
						H[(j - 1) + LH * k] -= sum * t2;
						H[(j - 1) + LH * (k + 1)] -= sum * t3;
					}
					if (wantz) {
						for (int j = iloz; j <= ihiz; j++) {
							const double sum = Z[(j - 1) + LZ * (k - 1)] + v2 * Z[(j - 1) + LZ * k] +
								v3 * Z[(j - 1) + LZ * (k + 1)];
							Z[(j - 1) + LZ * (k - 1)] -= sum * t1;
							Z[(j - 1) + LZ * k] -= sum * t2;
							Z[(j - 1) + LZ * (k + 1)] -= sum * t3;
						}
					}
				}
				else if (nr == 2) {
					for (int j = k; j <= i2; j++) {
						const double sum = H[(k - 1) + LH * (j - 1)] + v2 * H[k + LH * (j - 1)];
						H[(k - 1) + LH * (j - 1)] -= sum * t1;
						H[k + LH * (j - 1)] -= sum * t2;
					}
					for (int j = i1; j <= i; j++) {
						const double sum = H[(j - 1) + LH * (k - 1)] + v2 * H[(j - 1) + LH * k];
						H[(j - 1) + LH * (k - 1)] -= sum * t1;
						H[(j - 1) + LH * k] -= sum * t2;
					}
					if (wantz) {
						for (int j = iloz; j <= ihiz; j++) {
							const double sum = Z[(j - 1) + LZ * (k - 1)] + v2 * Z[(j - 1) + LZ * k];
							Z[(j - 1) + LZ * (k - 1)] -= sum * t1;
							Z[(j - 1) + LZ * k] -= sum * t2;
						}
					}
				}
			}
		}
		return i; // 非収束

	L150:
		if (l == i) {
			wr[i - 1] = H[(i - 1) + LH * (i - 1)];
			wi[i - 1] = 0.0;
		}
		else if (l == i - 1) {
			double cs, sn;
			rdlanv2(H[(i - 2) + LH * (i - 2)], H[(i - 2) + LH * (i - 1)],
				H[(i - 1) + LH * (i - 2)], H[(i - 1) + LH * (i - 1)],
				wr[i - 2], wi[i - 2], wr[i - 1], wi[i - 1], cs, sn, rm);
			if (wantt) {
				if (i2 > i) {
					rdrot(i2 - i, H + (i - 2) + LH * i, ldh, H + (i - 1) + LH * i, ldh, cs, sn, rm);
				}
				rdrot(i - i1 - 1, H + (i1 - 1) + LH * (i - 2), 1, H + (i1 - 1) + LH * (i - 1), 1, cs, sn, rm);
			}
			if (wantz) {
				rdrot(nz, Z + (iloz - 1) + LZ * (i - 2), 1, Z + (iloz - 1) + LZ * (i - 1), 1, cs, sn, rm);
			}
		}
		kdefl = 0;
		i = l - 1;
		goto L20;
	}
}

// Hessenberg 行列の固有値 (と Schur 形) — rdlahqr を全 size に使う簡略版
// job: 'E' 固有値のみ / 'S' Schur 形も, compz: 'N'/'I'/'V', ilo/ihi は 1-based
inline int rdhseqr(const char job, const char compz, const int n, const int ilo, const int ihi,
	double* H, const int ldh, double* wr, double* wi, double* Z, const int ldz, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool wantt = detail::option_is(job, 'S');
	const bool initz = detail::option_is(compz, 'I');
	const bool wantz = initz || detail::option_is(compz, 'V');
	if (!wantt && !detail::option_is(job, 'E')) {
		detail::rdlapack_error("rdhseqr: invalid job");
	}
	if (!wantz && !detail::option_is(compz, 'N')) {
		detail::rdlapack_error("rdhseqr: invalid compz");
	}
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		ldh < std::max(1, n) || (wantz && ldz < std::max(1, n)) || (!wantz && ldz < 1)) {
		detail::rdlapack_error("rdhseqr: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	if (initz) {
		dlaset('F', n, n, 0.0, 1.0, Z, ldz);
	}
	// ilo..ihi の外側の固有値は対角要素
	for (int i = 1; i <= ilo - 1; i++) {
		wr[i - 1] = H[(i - 1) + static_cast<std::size_t>(ldh) * (i - 1)];
		wi[i - 1] = 0.0;
	}
	for (int i = ihi + 1; i <= n; i++) {
		wr[i - 1] = H[(i - 1) + static_cast<std::size_t>(ldh) * (i - 1)];
		wi[i - 1] = 0.0;
	}
	return rdlahqr(wantt, wantz, n, ilo, ihi, H, ldh, wr, wi, 1, n, Z, ldz, rounding_mode);
}

// 1x1 / 2x2 の (ca*A - w*D) x = s*b を解く (固有 vector 計算用)
inline int rdlaln2(const bool ltrans, const int na, const int nw, const double smin, const double ca,
	const double* A, const int lda, const double d1, const double d2, const double* B, const int ldb,
	const double wr, const double wi, double* X, const int ldx, double& scale, double& xnorm,
	const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t LA = static_cast<std::size_t>(lda);
	const std::size_t LB = static_cast<std::size_t>(ldb);
	const std::size_t LX = static_cast<std::size_t>(ldx);
	static const bool zswap[4] = { false, false, true, true };
	static const bool rswap[4] = { false, true, false, true };
	static const int ipivot[4][4] = { { 1, 2, 3, 4 }, { 2, 1, 4, 3 }, { 3, 4, 1, 2 }, { 4, 3, 2, 1 } };
	const double smlnum = 2.0 * DBL_MIN;
	const double bignum = 1.0 / smlnum;
	const double smini = std::max(smin, smlnum);
	int info = 0;
	scale = 1.0;
	if (na == 1) {
		if (nw == 1) {
			double csr = ca * A[0] - wr * d1;
			double cnorm = std::fabs(csr);
			if (cnorm < smini) {
				csr = smini;
				cnorm = smini;
				info = 1;
			}
			const double bnorm = std::fabs(B[0]);
			if (cnorm < 1.0 && bnorm > 1.0) {
				if (bnorm > bignum * cnorm) {
					scale = 1.0 / bnorm;
				}
			}
			X[0] = (B[0] * scale) / csr;
			xnorm = std::fabs(X[0]);
		}
		else {
			double csr = ca * A[0] - wr * d1;
			double csi = -wi * d1;
			double cnorm = std::fabs(csr) + std::fabs(csi);
			if (cnorm < smini) {
				csr = smini;
				csi = 0.0;
				cnorm = smini;
				info = 1;
			}
			const double bnorm = std::fabs(B[0]) + std::fabs(B[LB]);
			if (cnorm < 1.0 && bnorm > 1.0) {
				if (bnorm > bignum * cnorm) {
					scale = 1.0 / bnorm;
				}
			}
			rdladiv(scale * B[0], scale * B[LB], csr, csi, X[0], X[LX], rounding_mode);
			xnorm = std::fabs(X[0]) + std::fabs(X[LX]);
		}
	}
	else {
		// 2x2: crv = (CR(1,1), CR(2,1), CR(1,2), CR(2,2))
		double crv[4];
		crv[0] = ca * A[0] - wr * d1;
		crv[3] = ca * A[1 + LA] - wr * d2;
		if (ltrans) {
			crv[2] = ca * A[1];
			crv[1] = ca * A[LA];
		}
		else {
			crv[1] = ca * A[1];
			crv[2] = ca * A[LA];
		}
		if (nw == 1) {
			double cmax = 0.0;
			int icmax = 0;
			for (int j = 1; j <= 4; j++) {
				if (std::fabs(crv[j - 1]) > cmax) {
					cmax = std::fabs(crv[j - 1]);
					icmax = j;
				}
			}
			if (cmax < smini) {
				const double bnorm = std::max(std::fabs(B[0]), std::fabs(B[1]));
				if (smini < 1.0 && bnorm > 1.0) {
					if (bnorm > bignum * smini) {
						scale = 1.0 / bnorm;
					}
				}
				const double temp = scale / smini;
				X[0] = temp * B[0];
				X[1] = temp * B[1];
				xnorm = temp * bnorm;
				return 1;
			}
			const double ur11 = crv[icmax - 1];
			const double cr21 = crv[ipivot[icmax - 1][1] - 1];
			const double ur12 = crv[ipivot[icmax - 1][2] - 1];
			const double cr22 = crv[ipivot[icmax - 1][3] - 1];
			const double ur11r = 1.0 / ur11;
			const double lr21 = ur11r * cr21;
			double ur22 = cr22 - ur12 * lr21;
			if (std::fabs(ur22) < smini) {
				ur22 = smini;
				info = 1;
			}
			double br1, br2;
			if (rswap[icmax - 1]) {
				br1 = B[1];
				br2 = B[0];
			}
			else {
				br1 = B[0];
				br2 = B[1];
			}
			br2 = br2 - lr21 * br1;
			const double bbnd = std::max(std::fabs(br1 * (ur22 * ur11r)), std::fabs(br2));
			if (bbnd > 1.0 && std::fabs(ur22) < 1.0) {
				if (bbnd >= bignum * std::fabs(ur22)) {
					scale = 1.0 / bbnd;
				}
			}
			const double xr2 = (br2 * scale) / ur22;
			const double xr1 = (scale * br1) * ur11r - xr2 * (ur11r * ur12);
			if (zswap[icmax - 1]) {
				X[0] = xr2;
				X[1] = xr1;
			}
			else {
				X[0] = xr1;
				X[1] = xr2;
			}
			xnorm = std::max(std::fabs(xr1), std::fabs(xr2));
			if (xnorm > 1.0 && cmax > 1.0) {
				if (xnorm > bignum / cmax) {
					const double temp = cmax / bignum;
					X[0] = temp * X[0];
					X[1] = temp * X[1];
					xnorm = temp * xnorm;
					scale = temp * scale;
				}
			}
		}
		else {
			double civ[4];
			civ[0] = -wi * d1;
			civ[1] = 0.0;
			civ[2] = 0.0;
			civ[3] = -wi * d2;
			double cmax = 0.0;
			int icmax = 0;
			for (int j = 1; j <= 4; j++) {
				if (std::fabs(crv[j - 1]) + std::fabs(civ[j - 1]) > cmax) {
					cmax = std::fabs(crv[j - 1]) + std::fabs(civ[j - 1]);
					icmax = j;
				}
			}
			if (cmax < smini) {
				const double bnorm = std::max(std::fabs(B[0]) + std::fabs(B[LB]),
					std::fabs(B[1]) + std::fabs(B[1 + LB]));
				if (smini < 1.0 && bnorm > 1.0) {
					if (bnorm > bignum * smini) {
						scale = 1.0 / bnorm;
					}
				}
				const double temp = scale / smini;
				X[0] = temp * B[0];
				X[1] = temp * B[1];
				X[LX] = temp * B[LB];
				X[1 + LX] = temp * B[1 + LB];
				xnorm = temp * bnorm;
				return 1;
			}
			const double ur11 = crv[icmax - 1];
			const double ui11 = civ[icmax - 1];
			const double cr21 = crv[ipivot[icmax - 1][1] - 1];
			const double ci21 = civ[ipivot[icmax - 1][1] - 1];
			const double ur12 = crv[ipivot[icmax - 1][2] - 1];
			const double ui12 = civ[ipivot[icmax - 1][2] - 1];
			const double cr22 = crv[ipivot[icmax - 1][3] - 1];
			const double ci22 = civ[ipivot[icmax - 1][3] - 1];
			double ur11r, ui11r, lr21, li21, ur12s, ui12s, ur22, ui22;
			if (icmax == 1 || icmax == 4) {
				if (std::fabs(ur11) > std::fabs(ui11)) {
					const double temp = ui11 / ur11;
					ur11r = 1.0 / (ur11 * (1.0 + temp * temp));
					ui11r = -temp * ur11r;
				}
				else {
					const double temp = ur11 / ui11;
					ui11r = -1.0 / (ui11 * (1.0 + temp * temp));
					ur11r = -temp * ui11r;
				}
				lr21 = cr21 * ur11r;
				li21 = cr21 * ui11r;
				ur12s = ur12 * ur11r;
				ui12s = ur12 * ui11r;
				ur22 = cr22 - ur12 * lr21;
				ui22 = ci22 - ur12 * li21;
			}
			else {
				ur11r = 1.0 / ur11;
				ui11r = 0.0;
				lr21 = cr21 * ur11r;
				li21 = ci21 * ur11r;
				ur12s = ur12 * ur11r;
				ui12s = ui12 * ur11r;
				ur22 = cr22 - ur12 * lr21 + ui12 * li21;
				ui22 = -ur12 * li21 - ui12 * lr21;
			}
			const double u22abs = std::fabs(ur22) + std::fabs(ui22);
			if (u22abs < smini) {
				ur22 = smini;
				ui22 = 0.0;
				info = 1;
			}
			double br1, br2, bi1, bi2;
			if (rswap[icmax - 1]) {
				br2 = B[0];
				br1 = B[1];
				bi2 = B[LB];
				bi1 = B[1 + LB];
			}
			else {
				br1 = B[0];
				br2 = B[1];
				bi1 = B[LB];
				bi2 = B[1 + LB];
			}
			br2 = br2 - lr21 * br1 + li21 * bi1;
			bi2 = bi2 - li21 * br1 - lr21 * bi1;
			const double bbnd = std::max((std::fabs(br1) + std::fabs(bi1)) *
				(u22abs * (std::fabs(ur11r) + std::fabs(ui11r))),
				std::fabs(br2) + std::fabs(bi2));
			if (bbnd > 1.0 && u22abs < 1.0) {
				if (bbnd >= bignum * u22abs) {
					scale = 1.0 / bbnd;
					br1 = scale * br1;
					bi1 = scale * bi1;
					br2 = scale * br2;
					bi2 = scale * bi2;
				}
			}
			double xr2, xi2;
			rdladiv(br2, bi2, ur22, ui22, xr2, xi2, rounding_mode);
			const double xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
			const double xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
			if (zswap[icmax - 1]) {
				X[0] = xr2;
				X[1] = xr1;
				X[LX] = xi2;
				X[1 + LX] = xi1;
			}
			else {
				X[0] = xr1;
				X[1] = xr2;
				X[LX] = xi1;
				X[1 + LX] = xi2;
			}
			xnorm = std::max(std::fabs(xr1) + std::fabs(xi1), std::fabs(xr2) + std::fabs(xi2));
			if (xnorm > 1.0 && cmax > 1.0) {
				if (xnorm > bignum / cmax) {
					const double temp = cmax / bignum;
					X[0] = temp * X[0];
					X[1] = temp * X[1];
					X[LX] = temp * X[LX];
					X[1 + LX] = temp * X[1 + LX];
					xnorm = temp * xnorm;
					scale = temp * scale;
				}
			}
		}
	}
	return info;
}

// 実 Schur 形 T の固有 vector (howmny: 'A' 全て / 'B' 全てを back-transform)
// side: 'R' 右 / 'L' 左 / 'B' 両方．SELECT 指定 ('S') は未対応．
inline int rdtrevc(const char side, const char howmny, const int n, const double* T, const int ldt,
	double* VL, const int ldvl, double* VR, const int ldvr, const int mm, int& m, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool bothv = detail::option_is(side, 'B');
	const bool rightv = detail::option_is(side, 'R') || bothv;
	const bool leftv = detail::option_is(side, 'L') || bothv;
	const bool over = detail::option_is(howmny, 'B');
	if ((!rightv && !leftv) || (!detail::option_is(howmny, 'A') && !over)) {
		detail::rdlapack_error("rdtrevc: invalid side/howmny (SELECT は未対応)");
	}
	m = n;
	if (n < 0 || ldt < std::max(1, n) || ldvl < 1 || (leftv && ldvl < n) ||
		ldvr < 1 || (rightv && ldvr < n) || mm < m) {
		detail::rdlapack_error("rdtrevc: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t LT = static_cast<std::size_t>(ldt);
	const std::size_t LVL = static_cast<std::size_t>(ldvl);
	const std::size_t LVR = static_cast<std::size_t>(ldvr);
	const double unfl = DBL_MIN;
	const double ulp = DBL_EPSILON;
	const double smlnum = unfl * (n / ulp);
	const double bignum = (1.0 - ulp) / smlnum;
	std::vector<double> work(static_cast<std::size_t>(3) * n + 1);
	double* wk = work.data(); // wk[X-1] が Fortran の WORK(X)
	double x[4]; // X(2,2), ldx=2
	wk[0] = 0.0;
	for (int j = 2; j <= n; j++) {
		wk[j - 1] = 0.0;
		for (int i = 1; i <= j - 1; i++) {
			wk[j - 1] += std::fabs(T[(i - 1) + LT * (j - 1)]);
		}
	}
	const int n2 = 2 * n;
	double scale2, xnorm;
	if (rightv) {
		int ip = 0;
		int is = m;
		for (int ki = n; ki >= 1; ki--) {
			bool skip = false;
			if (ip == 1) {
				skip = true;
			}
			else if (ki != 1 && T[(ki - 1) + LT * (ki - 2)] != 0.0) {
				ip = -1;
			}
			if (!skip) {
				const double wr = T[(ki - 1) + LT * (ki - 1)];
				double wi = 0.0;
				if (ip != 0) {
					wi = std::sqrt(std::fabs(T[(ki - 2) + LT * (ki - 1)])) *
						std::sqrt(std::fabs(T[(ki - 1) + LT * (ki - 2)]));
				}
				const double smin = std::max(ulp * (std::fabs(wr) + std::fabs(wi)), smlnum);
				if (ip == 0) {
					// 実固有値
					wk[ki + n - 1] = 1.0;
					for (int k = 1; k <= ki - 1; k++) {
						wk[k + n - 1] = -T[(k - 1) + LT * (ki - 1)];
					}
					int jnxt = ki - 1;
					for (int j = ki - 1; j >= 1; j--) {
						if (j > jnxt) {
							continue;
						}
						int j1 = j;
						jnxt = j - 1;
						if (j > 1 && T[(j - 1) + LT * (j - 2)] != 0.0) {
							j1 = j - 1;
							jnxt = j - 2;
						}
						if (j1 == j) {
							rdlaln2(false, 1, 1, smin, 1.0, T + (j - 1) + LT * (j - 1), ldt, 1.0, 1.0,
								wk + (j + n - 1), n, wr, 0.0, x, 2, scale2, xnorm, rm);
							if (xnorm > 1.0 && wk[j - 1] > bignum / xnorm) {
								x[0] = x[0] / xnorm;
								scale2 = scale2 / xnorm;
							}
							if (scale2 != 1.0) {
								rdscal(ki, scale2, wk + n, 1, rm);
							}
							wk[j + n - 1] = x[0];
							rdaxpy(j - 1, -x[0], T + LT * (j - 1), 1, wk + n, 1, rm);
						}
						else {
							rdlaln2(false, 2, 1, smin, 1.0, T + (j - 2) + LT * (j - 2), ldt, 1.0, 1.0,
								wk + (j + n - 2), n, wr, 0.0, x, 2, scale2, xnorm, rm);
							if (xnorm > 1.0) {
								const double beta = std::max(wk[j - 2], wk[j - 1]);
								if (beta > bignum / xnorm) {
									x[0] = x[0] / xnorm;
									x[1] = x[1] / xnorm;
									scale2 = scale2 / xnorm;
								}
							}
							if (scale2 != 1.0) {
								rdscal(ki, scale2, wk + n, 1, rm);
							}
							wk[j + n - 2] = x[0];
							wk[j + n - 1] = x[1];
							rdaxpy(j - 2, -x[0], T + LT * (j - 2), 1, wk + n, 1, rm);
							rdaxpy(j - 2, -x[1], T + LT * (j - 1), 1, wk + n, 1, rm);
						}
					}
					if (!over) {
						dcopy(ki, wk + n, 1, VR + LVR * (is - 1), 1);
						const int ii = idamax(ki, VR + LVR * (is - 1), 1) + 1;
						const double remax = 1.0 / std::fabs(VR[(ii - 1) + LVR * (is - 1)]);
						rdscal(ki, remax, VR + LVR * (is - 1), 1, rm);
						for (int k = ki + 1; k <= n; k++) {
							VR[(k - 1) + LVR * (is - 1)] = 0.0;
						}
					}
					else {
						if (ki > 1) {
							rdgemv('N', n, ki - 1, 1.0, VR, ldvr, wk + n, 1, wk[ki + n - 1],
								VR + LVR * (ki - 1), 1, rm);
						}
						const int ii = idamax(n, VR + LVR * (ki - 1), 1) + 1;
						const double remax = 1.0 / std::fabs(VR[(ii - 1) + LVR * (ki - 1)]);
						rdscal(n, remax, VR + LVR * (ki - 1), 1, rm);
					}
				}
				else {
					// 複素共役対
					if (std::fabs(T[(ki - 2) + LT * (ki - 1)]) >= std::fabs(T[(ki - 1) + LT * (ki - 2)])) {
						wk[ki + n - 2] = 1.0;
						wk[ki + n2 - 1] = wi / T[(ki - 2) + LT * (ki - 1)];
					}
					else {
						wk[ki + n - 2] = -wi / T[(ki - 1) + LT * (ki - 2)];
						wk[ki + n2 - 1] = 1.0;
					}
					wk[ki + n - 1] = 0.0;
					wk[ki + n2 - 2] = 0.0;
					for (int k = 1; k <= ki - 2; k++) {
						wk[k + n - 1] = -wk[ki + n - 2] * T[(k - 1) + LT * (ki - 2)];
						wk[k + n2 - 1] = -wk[ki + n2 - 1] * T[(k - 1) + LT * (ki - 1)];
					}
					int jnxt = ki - 2;
					for (int j = ki - 2; j >= 1; j--) {
						if (j > jnxt) {
							continue;
						}
						int j1 = j;
						jnxt = j - 1;
						if (j > 1 && T[(j - 1) + LT * (j - 2)] != 0.0) {
							j1 = j - 1;
							jnxt = j - 2;
						}
						if (j1 == j) {
							rdlaln2(false, 1, 2, smin, 1.0, T + (j - 1) + LT * (j - 1), ldt, 1.0, 1.0,
								wk + (j + n - 1), n, wr, wi, x, 2, scale2, xnorm, rm);
							if (xnorm > 1.0 && wk[j - 1] > bignum / xnorm) {
								x[0] = x[0] / xnorm;
								x[2] = x[2] / xnorm;
								scale2 = scale2 / xnorm;
							}
							if (scale2 != 1.0) {
								rdscal(ki, scale2, wk + n, 1, rm);
								rdscal(ki, scale2, wk + n2, 1, rm);
							}
							wk[j + n - 1] = x[0];
							wk[j + n2 - 1] = x[2];
							rdaxpy(j - 1, -x[0], T + LT * (j - 1), 1, wk + n, 1, rm);
							rdaxpy(j - 1, -x[2], T + LT * (j - 1), 1, wk + n2, 1, rm);
						}
						else {
							rdlaln2(false, 2, 2, smin, 1.0, T + (j - 2) + LT * (j - 2), ldt, 1.0, 1.0,
								wk + (j + n - 2), n, wr, wi, x, 2, scale2, xnorm, rm);
							if (xnorm > 1.0) {
								const double beta = std::max(wk[j - 2], wk[j - 1]);
								if (beta > bignum / xnorm) {
									const double rec = 1.0 / xnorm;
									x[0] *= rec;
									x[2] *= rec;
									x[1] *= rec;
									x[3] *= rec;
									scale2 *= rec;
								}
							}
							if (scale2 != 1.0) {
								rdscal(ki, scale2, wk + n, 1, rm);
								rdscal(ki, scale2, wk + n2, 1, rm);
							}
							wk[j + n - 2] = x[0];
							wk[j + n - 1] = x[1];
							wk[j + n2 - 2] = x[2];
							wk[j + n2 - 1] = x[3];
							rdaxpy(j - 2, -x[0], T + LT * (j - 2), 1, wk + n, 1, rm);
							rdaxpy(j - 2, -x[1], T + LT * (j - 1), 1, wk + n, 1, rm);
							rdaxpy(j - 2, -x[2], T + LT * (j - 2), 1, wk + n2, 1, rm);
							rdaxpy(j - 2, -x[3], T + LT * (j - 1), 1, wk + n2, 1, rm);
						}
					}
					if (!over) {
						dcopy(ki, wk + n, 1, VR + LVR * (is - 2), 1);
						dcopy(ki, wk + n2, 1, VR + LVR * (is - 1), 1);
						double emax = 0.0;
						for (int k = 1; k <= ki; k++) {
							emax = std::max(emax, std::fabs(VR[(k - 1) + LVR * (is - 2)]) +
								std::fabs(VR[(k - 1) + LVR * (is - 1)]));
						}
						const double remax = 1.0 / emax;
						rdscal(ki, remax, VR + LVR * (is - 2), 1, rm);
						rdscal(ki, remax, VR + LVR * (is - 1), 1, rm);
						for (int k = ki + 1; k <= n; k++) {
							VR[(k - 1) + LVR * (is - 2)] = 0.0;
							VR[(k - 1) + LVR * (is - 1)] = 0.0;
						}
					}
					else {
						if (ki > 2) {
							rdgemv('N', n, ki - 2, 1.0, VR, ldvr, wk + n, 1, wk[ki + n - 2],
								VR + LVR * (ki - 2), 1, rm);
							rdgemv('N', n, ki - 2, 1.0, VR, ldvr, wk + n2, 1, wk[ki + n2 - 1],
								VR + LVR * (ki - 1), 1, rm);
						}
						else {
							rdscal(n, wk[ki + n - 2], VR + LVR * (ki - 2), 1, rm);
							rdscal(n, wk[ki + n2 - 1], VR + LVR * (ki - 1), 1, rm);
						}
						double emax = 0.0;
						for (int k = 1; k <= n; k++) {
							emax = std::max(emax, std::fabs(VR[(k - 1) + LVR * (ki - 2)]) +
								std::fabs(VR[(k - 1) + LVR * (ki - 1)]));
						}
						const double remax = 1.0 / emax;
						rdscal(n, remax, VR + LVR * (ki - 2), 1, rm);
						rdscal(n, remax, VR + LVR * (ki - 1), 1, rm);
					}
				}
				is = is - 1;
				if (ip != 0) {
					is = is - 1;
				}
			}
			if (ip == 1) {
				ip = 0;
			}
			else if (ip == -1) {
				ip = 1;
			}
		}
	}
	if (leftv) {
		int ip = 0;
		int is = 1;
		for (int ki = 1; ki <= n; ki++) {
			bool skip = false;
			if (ip == -1) {
				skip = true;
			}
			else if (ki != n && T[ki + LT * (ki - 1)] != 0.0) {
				ip = 1;
			}
			if (!skip) {
				const double wr = T[(ki - 1) + LT * (ki - 1)];
				double wi = 0.0;
				if (ip != 0) {
					wi = std::sqrt(std::fabs(T[(ki - 1) + LT * ki])) *
						std::sqrt(std::fabs(T[ki + LT * (ki - 1)]));
				}
				const double smin = std::max(ulp * (std::fabs(wr) + std::fabs(wi)), smlnum);
				if (ip == 0) {
					// 実固有値
					wk[ki + n - 1] = 1.0;
					for (int k = ki + 1; k <= n; k++) {
						wk[k + n - 1] = -T[(ki - 1) + LT * (k - 1)];
					}
					double vmax = 1.0;
					double vcrit = bignum;
					int jnxt = ki + 1;
					for (int j = ki + 1; j <= n; j++) {
						if (j < jnxt) {
							continue;
						}
						int j2 = j;
						jnxt = j + 1;
						if (j < n && T[j + LT * (j - 1)] != 0.0) {
							j2 = j + 1;
							jnxt = j + 2;
						}
						if (j2 == j) {
							if (wk[j - 1] > vcrit) {
								const double rec = 1.0 / vmax;
								rdscal(n - ki + 1, rec, wk + (ki + n - 1), 1, rm);
								vmax = 1.0;
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								rddot(j - ki - 1, T + ki + LT * (j - 1), 1, wk + (ki + n), 1, rm);
							rdlaln2(true, 1, 1, smin, 1.0, T + (j - 1) + LT * (j - 1), ldt, 1.0, 1.0,
								wk + (j + n - 1), n, wr, 0.0, x, 2, scale2, xnorm, rm);
							if (scale2 != 1.0) {
								rdscal(n - ki + 1, scale2, wk + (ki + n - 1), 1, rm);
							}
							wk[j + n - 1] = x[0];
							vmax = std::max(std::fabs(wk[j + n - 1]), vmax);
							vcrit = bignum / vmax;
						}
						else {
							const double beta = std::max(wk[j - 1], wk[j]);
							if (beta > vcrit) {
								const double rec = 1.0 / vmax;
								rdscal(n - ki + 1, rec, wk + (ki + n - 1), 1, rm);
								vmax = 1.0;
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								rddot(j - ki - 1, T + ki + LT * (j - 1), 1, wk + (ki + n), 1, rm);
							wk[j + n] = wk[j + n] -
								rddot(j - ki - 1, T + ki + LT * j, 1, wk + (ki + n), 1, rm);
							rdlaln2(true, 2, 1, smin, 1.0, T + (j - 1) + LT * (j - 1), ldt, 1.0, 1.0,
								wk + (j + n - 1), n, wr, 0.0, x, 2, scale2, xnorm, rm);
							if (scale2 != 1.0) {
								rdscal(n - ki + 1, scale2, wk + (ki + n - 1), 1, rm);
							}
							wk[j + n - 1] = x[0];
							wk[j + n] = x[1];
							vmax = std::max(std::max(std::fabs(wk[j + n - 1]), std::fabs(wk[j + n])), vmax);
							vcrit = bignum / vmax;
						}
					}
					if (!over) {
						dcopy(n - ki + 1, wk + (ki + n - 1), 1, VL + (ki - 1) + LVL * (is - 1), 1);
						const int ii = idamax(n - ki + 1, VL + (ki - 1) + LVL * (is - 1), 1) + ki;
						const double remax = 1.0 / std::fabs(VL[(ii - 1) + LVL * (is - 1)]);
						rdscal(n - ki + 1, remax, VL + (ki - 1) + LVL * (is - 1), 1, rm);
						for (int k = 1; k <= ki - 1; k++) {
							VL[(k - 1) + LVL * (is - 1)] = 0.0;
						}
					}
					else {
						if (ki < n) {
							rdgemv('N', n, n - ki, 1.0, VL + LVL * ki, ldvl, wk + (ki + n), 1,
								wk[ki + n - 1], VL + LVL * (ki - 1), 1, rm);
						}
						const int ii = idamax(n, VL + LVL * (ki - 1), 1) + 1;
						const double remax = 1.0 / std::fabs(VL[(ii - 1) + LVL * (ki - 1)]);
						rdscal(n, remax, VL + LVL * (ki - 1), 1, rm);
					}
				}
				else {
					// 複素共役対
					if (std::fabs(T[(ki - 1) + LT * ki]) >= std::fabs(T[ki + LT * (ki - 1)])) {
						wk[ki + n - 1] = wi / T[(ki - 1) + LT * ki];
						wk[ki + n2] = 1.0;
					}
					else {
						wk[ki + n - 1] = 1.0;
						wk[ki + n2] = -wi / T[ki + LT * (ki - 1)];
					}
					wk[ki + n] = 0.0;
					wk[ki + n2 - 1] = 0.0;
					for (int k = ki + 2; k <= n; k++) {
						wk[k + n - 1] = -wk[ki + n - 1] * T[(ki - 1) + LT * (k - 1)];
						wk[k + n2 - 1] = -wk[ki + n2] * T[ki + LT * (k - 1)];
					}
					double vmax = 1.0;
					double vcrit = bignum;
					int jnxt = ki + 2;
					for (int j = ki + 2; j <= n; j++) {
						if (j < jnxt) {
							continue;
						}
						int j2 = j;
						jnxt = j + 1;
						if (j < n && T[j + LT * (j - 1)] != 0.0) {
							j2 = j + 1;
							jnxt = j + 2;
						}
						if (j2 == j) {
							if (wk[j - 1] > vcrit) {
								const double rec = 1.0 / vmax;
								rdscal(n - ki + 1, rec, wk + (ki + n - 1), 1, rm);
								rdscal(n - ki + 1, rec, wk + (ki + n2 - 1), 1, rm);
								vmax = 1.0;
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								rddot(j - ki - 2, T + (ki + 1) + LT * (j - 1), 1, wk + (ki + n + 1), 1, rm);
							wk[j + n2 - 1] = wk[j + n2 - 1] -
								rddot(j - ki - 2, T + (ki + 1) + LT * (j - 1), 1, wk + (ki + n2 + 1), 1, rm);
							rdlaln2(false, 1, 2, smin, 1.0, T + (j - 1) + LT * (j - 1), ldt, 1.0, 1.0,
								wk + (j + n - 1), n, wr, -wi, x, 2, scale2, xnorm, rm);
							if (scale2 != 1.0) {
								rdscal(n - ki + 1, scale2, wk + (ki + n - 1), 1, rm);
								rdscal(n - ki + 1, scale2, wk + (ki + n2 - 1), 1, rm);
							}
							wk[j + n - 1] = x[0];
							wk[j + n2 - 1] = x[2];
							vmax = std::max(std::max(std::fabs(wk[j + n - 1]), std::fabs(wk[j + n2 - 1])), vmax);
							vcrit = bignum / vmax;
						}
						else {
							const double beta = std::max(wk[j - 1], wk[j]);
							if (beta > vcrit) {
								const double rec = 1.0 / vmax;
								rdscal(n - ki + 1, rec, wk + (ki + n - 1), 1, rm);
								rdscal(n - ki + 1, rec, wk + (ki + n2 - 1), 1, rm);
								vmax = 1.0;
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								rddot(j - ki - 2, T + (ki + 1) + LT * (j - 1), 1, wk + (ki + n + 1), 1, rm);
							wk[j + n2 - 1] = wk[j + n2 - 1] -
								rddot(j - ki - 2, T + (ki + 1) + LT * (j - 1), 1, wk + (ki + n2 + 1), 1, rm);
							wk[j + n] = wk[j + n] -
								rddot(j - ki - 2, T + (ki + 1) + LT * j, 1, wk + (ki + n + 1), 1, rm);
							wk[j + n2] = wk[j + n2] -
								rddot(j - ki - 2, T + (ki + 1) + LT * j, 1, wk + (ki + n2 + 1), 1, rm);
							rdlaln2(true, 2, 2, smin, 1.0, T + (j - 1) + LT * (j - 1), ldt, 1.0, 1.0,
								wk + (j + n - 1), n, wr, -wi, x, 2, scale2, xnorm, rm);
							if (scale2 != 1.0) {
								rdscal(n - ki + 1, scale2, wk + (ki + n - 1), 1, rm);
								rdscal(n - ki + 1, scale2, wk + (ki + n2 - 1), 1, rm);
							}
							wk[j + n - 1] = x[0];
							wk[j + n2 - 1] = x[2];
							wk[j + n] = x[1];
							wk[j + n2] = x[3];
							vmax = std::max(std::max(std::fabs(x[0]), std::fabs(x[2])),
								std::max(std::max(std::fabs(x[1]), std::fabs(x[3])), vmax));
							vcrit = bignum / vmax;
						}
					}
					if (!over) {
						dcopy(n - ki + 1, wk + (ki + n - 1), 1, VL + (ki - 1) + LVL * (is - 1), 1);
						dcopy(n - ki + 1, wk + (ki + n2 - 1), 1, VL + (ki - 1) + LVL * is, 1);
						double emax = 0.0;
						for (int k = ki; k <= n; k++) {
							emax = std::max(emax, std::fabs(VL[(k - 1) + LVL * (is - 1)]) +
								std::fabs(VL[(k - 1) + LVL * is]));
						}
						const double remax = 1.0 / emax;
						rdscal(n - ki + 1, remax, VL + (ki - 1) + LVL * (is - 1), 1, rm);
						rdscal(n - ki + 1, remax, VL + (ki - 1) + LVL * is, 1, rm);
						for (int k = 1; k <= ki - 1; k++) {
							VL[(k - 1) + LVL * (is - 1)] = 0.0;
							VL[(k - 1) + LVL * is] = 0.0;
						}
					}
					else {
						if (ki < n - 1) {
							rdgemv('N', n, n - ki - 1, 1.0, VL + LVL * (ki + 1), ldvl,
								wk + (ki + n + 1), 1, wk[ki + n - 1], VL + LVL * (ki - 1), 1, rm);
							rdgemv('N', n, n - ki - 1, 1.0, VL + LVL * (ki + 1), ldvl,
								wk + (ki + n2 + 1), 1, wk[ki + n2], VL + LVL * ki, 1, rm);
						}
						else {
							rdscal(n, wk[ki + n - 1], VL + LVL * (ki - 1), 1, rm);
							rdscal(n, wk[ki + n2], VL + LVL * ki, 1, rm);
						}
						double emax = 0.0;
						for (int k = 1; k <= n; k++) {
							emax = std::max(emax, std::fabs(VL[(k - 1) + LVL * (ki - 1)]) +
								std::fabs(VL[(k - 1) + LVL * ki]));
						}
						const double remax = 1.0 / emax;
						rdscal(n, remax, VL + LVL * (ki - 1), 1, rm);
						rdscal(n, remax, VL + LVL * ki, 1, rm);
					}
				}
				is = is + 1;
				if (ip != 0) {
					is = is + 1;
				}
			}
			if (ip == -1) {
				ip = 0;
			}
			else if (ip == 1) {
				ip = -1;
			}
		}
	}
	return 0;
}

// 非対称行列の固有値・固有 vector (driver)
// jobvl/jobvr: 'V' 計算する / 'N' しない．固有値は wr + i*wi (複素対は連続して
// wi[k] > 0, wi[k+1] < 0)．複素対の固有 vector は VL/VR の連続 2 列に
// (実部, 虚部) で格納される (LAPACK と同じ)．
inline int rdgeev(const char jobvl, const char jobvr, const int n, double* A, const int lda,
	double* wr, double* wi, double* VL, const int ldvl, double* VR, const int ldvr,
	const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool wantvl = detail::option_is(jobvl, 'V');
	const bool wantvr = detail::option_is(jobvr, 'V');
	if ((!wantvl && !detail::option_is(jobvl, 'N')) || (!wantvr && !detail::option_is(jobvr, 'N'))) {
		detail::rdlapack_error("rdgeev: invalid jobvl/jobvr");
	}
	if (n < 0 || lda < std::max(1, n) || ldvl < 1 || (wantvl && ldvl < n) ||
		ldvr < 1 || (wantvr && ldvr < n)) {
		detail::rdlapack_error("rdgeev: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t LVL = static_cast<std::size_t>(ldvl);
	const std::size_t LVR = static_cast<std::size_t>(ldvr);
	const double eps = DBL_EPSILON;
	double smlnum = DBL_MIN;
	smlnum = std::sqrt(smlnum) / eps;
	const double bignum = 1.0 / smlnum;
	const double anrm = rdlange('M', n, n, A, lda, rm);
	bool scalea = false;
	double cscale = 1.0;
	if (anrm > 0.0 && anrm < smlnum) {
		scalea = true;
		cscale = smlnum;
	}
	else if (anrm > bignum) {
		scalea = true;
		cscale = bignum;
	}
	if (scalea) {
		rdlascl('G', 0, 0, anrm, cscale, n, n, A, lda, rm);
	}
	std::vector<double> balscale(static_cast<std::size_t>(n));
	std::vector<double> tau(static_cast<std::size_t>(std::max(1, n - 1)));
	int ilo, ihi;
	rdgebal('B', n, A, lda, ilo, ihi, balscale.data(), rm);
	rdgehrd(n, ilo, ihi, A, lda, tau.data(), rm);
	int info;
	char side = 'N';
	if (wantvl) {
		side = 'L';
		dlacpy('L', n, n, A, lda, VL, ldvl);
		rdorghr(n, ilo, ihi, VL, ldvl, tau.data(), rm);
		info = rdhseqr('S', 'V', n, ilo, ihi, A, lda, wr, wi, VL, ldvl, rm);
		if (wantvr) {
			side = 'B';
			dlacpy('F', n, n, VL, ldvl, VR, ldvr);
		}
	}
	else if (wantvr) {
		side = 'R';
		dlacpy('L', n, n, A, lda, VR, ldvr);
		rdorghr(n, ilo, ihi, VR, ldvr, tau.data(), rm);
		info = rdhseqr('S', 'V', n, ilo, ihi, A, lda, wr, wi, VR, ldvr, rm);
	}
	else {
		info = rdhseqr('E', 'N', n, ilo, ihi, A, lda, wr, wi, VR, ldvr, rm);
	}
	if (info == 0) {
		if (wantvl || wantvr) {
			int nout;
			rdtrevc(side, 'B', n, A, lda, VL, ldvl, VR, ldvr, n, nout, rm);
		}
		std::vector<double> tmp(static_cast<std::size_t>(n));
		for (int pass = 0; pass < 2; pass++) {
			double* V = (pass == 0) ? VL : VR;
			const std::size_t LV = (pass == 0) ? LVL : LVR;
			const int ldv = (pass == 0) ? ldvl : ldvr;
			if ((pass == 0 && !wantvl) || (pass == 1 && !wantvr)) {
				continue;
			}
			rdgebak('B', (pass == 0) ? 'L' : 'R', n, ilo, ihi, balscale.data(), n, V, ldv, rm);
			for (int i = 1; i <= n; i++) {
				if (wi[i - 1] == 0.0) {
					const double scl = 1.0 / rdnrm2(n, V + LV * (i - 1), 1, rm);
					rdscal(n, scl, V + LV * (i - 1), 1, rm);
				}
				else if (wi[i - 1] > 0.0) {
					const double scl = 1.0 / rdlapy2(rdnrm2(n, V + LV * (i - 1), 1, rm),
						rdnrm2(n, V + LV * i, 1, rm), rm);
					rdscal(n, scl, V + LV * (i - 1), 1, rm);
					rdscal(n, scl, V + LV * i, 1, rm);
					for (int k = 1; k <= n; k++) {
						tmp[k - 1] = V[(k - 1) + LV * (i - 1)] * V[(k - 1) + LV * (i - 1)] +
							V[(k - 1) + LV * i] * V[(k - 1) + LV * i];
					}
					const int k = idamax(n, tmp.data(), 1) + 1;
					double cs, sn, r;
					rdlartg(V[(k - 1) + LV * (i - 1)], V[(k - 1) + LV * i], cs, sn, r, rm);
					rdrot(n, V + LV * (i - 1), 1, V + LV * i, 1, cs, sn, rm);
					V[(k - 1) + LV * i] = 0.0;
				}
			}
		}
	}
	if (scalea) {
		rdlascl('G', 0, 0, cscale, anrm, n - info, 1, wr + info, std::max(n - info, 1), rm);
		rdlascl('G', 0, 0, cscale, anrm, n - info, 1, wi + info, std::max(n - info, 1), rm);
		if (info > 0) {
			rdlascl('G', 0, 0, cscale, anrm, ilo - 1, 1, wr, n, rm);
			rdlascl('G', 0, 0, cscale, anrm, ilo - 1, 1, wi, n, rm);
		}
	}
	return info;
}

#endif // VBLAS_RDLAPACK_GEEV_HPP
