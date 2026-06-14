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

#ifndef TLAPACK_TLAPACK_GEEV_HPP
#define TLAPACK_TLAPACK_GEEV_HPP

// 非対称固有値問題 (template 版)．
// reference LAPACK 3.12.1 の dgebal/dgebak/dgehd2/dlahr2/dgehrd/dorghr/
// dlahqr/dlaln2/dtrevc/dgeev を移植．
//
// - ilo/ihi と scale の swap 先 index は LAPACK と同じ 1-based (注意)．
// - thseqr は dlahqr (T-shift QR) を全 size に使う簡略版
//   (reference の dlaqr0 系 multishift は移植していない．結果の意味は同じ)．
// - ttrevc は howmny = 'A' / 'B' のみ (SELECT 指定は未対応)．

#include "tlapack_svd.hpp"


// 行列の balancing (permute + scale)．ilo/ihi は 1-based で返る
template <typename T>
inline int tgebal(const char job, const int n, T* A, const int lda,
	int& ilo, int& ihi, T* scale) {
	namespace detail = tlapack_detail;
	if (!detail::option_is(job, 'N') && !detail::option_is(job, 'P') &&
		!detail::option_is(job, 'S') && !detail::option_is(job, 'B')) {
		detail::tlapack_error("tgebal: invalid job");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tgebal: invalid n/lda");
	}
	if (n == 0) {
		ilo = 1;
		ihi = 0;
		return 0;
	}
	if (detail::option_is(job, 'N')) {
		for (int i = 0; i < n; i++) {
			scale[i] = T(1);
		}
		ilo = 1;
		ihi = n;
		return 0;
	}
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
					if (i != j && A[(i - 1) + L * (j - 1)] != T(0)) {
						canswap = false;
						break;
					}
				}
				if (canswap) {
					scale[l - 1] = static_cast<T>(i);
					if (i != l) {
						tswap(l, A + L * (i - 1), 1, A + L * (l - 1), 1);
						tswap(n - k + 1, A + (i - 1) + L * (k - 1), lda, A + (l - 1) + L * (k - 1), lda);
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
					if (i != j && A[(i - 1) + L * (j - 1)] != T(0)) {
						canswap = false;
						break;
					}
				}
				if (canswap) {
					scale[k - 1] = static_cast<T>(j);
					if (j != k) {
						tswap(l, A + L * (j - 1), 1, A + L * (k - 1), 1);
						tswap(n - k + 1, A + (j - 1) + L * (k - 1), lda, A + (k - 1) + L * (k - 1), lda);
					}
					noconv = true;
					k = k + 1;
				}
			}
		}
	}
	for (int i = k; i <= l; i++) {
		scale[i - 1] = T(1);
	}
	if (detail::option_is(job, 'P')) {
		ilo = k;
		ihi = l;
		return 0;
	}
	// scaling (基数 2 の冪で行/列 norm を揃える)
	const T sclfac = T(2);
	const T factor = (T(95) / T(100));
	const T sfmin1 = tlamch<T>('S') / tlamch<T>('P');
	const T sfmax1 = T(1) / sfmin1;
	const T sfmin2 = sfmin1 * sclfac;
	const T sfmax2 = T(1) / sfmin2;
	bool noconv = true;
	while (noconv) {
		noconv = false;
		for (int i = k; i <= l; i++) {
			T c = tnrm2(l - k + 1, A + (k - 1) + L * (i - 1), 1);
			T r = tnrm2(l - k + 1, A + (i - 1) + L * (k - 1), lda);
			const int ica = itamax(l, A + L * (i - 1), 1) + 1;
			T ca = tlapack_detail::tabs(A[(ica - 1) + L * (i - 1)]);
			const int ira = itamax(n - k + 1, A + (i - 1) + L * (k - 1), lda) + 1;
			T ra = tlapack_detail::tabs(A[(i - 1) + L * (ira + k - 2)]);
			if (c == T(0) || r == T(0)) {
				continue;
			}
			if (tlapack_detail::tisnan(c + ca + r + ra)) {
				detail::tlapack_error("tgebal: NaN in matrix");
			}
			T g = r / sclfac;
			T f = T(1);
			const T s = c + r;
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
			if (f < T(1) && scale[i - 1] < T(1)) {
				if (f * scale[i - 1] <= sfmin1) {
					continue;
				}
			}
			if (f > T(1) && scale[i - 1] > T(1)) {
				if (scale[i - 1] >= sfmax1 / f) {
					continue;
				}
			}
			g = T(1) / f;
			scale[i - 1] = scale[i - 1] * f;
			noconv = true;
			tscal(n - k + 1, g, A + (i - 1) + L * (k - 1), lda);
			tscal(l, f, A + L * (i - 1), 1);
		}
	}
	ilo = k;
	ihi = l;
	return 0;
}

// balancing の逆変換を固有 vector へ適用 (ilo/ihi は 1-based)
template <typename T>
inline int tgebak(const char job, const char side, const int n, const int ilo, const int ihi,
	const T* scale, const int m, T* V, const int ldv) {
	namespace detail = tlapack_detail;
	const bool rightv = detail::option_is(side, 'R');
	const bool leftv = detail::option_is(side, 'L');
	if ((!detail::option_is(job, 'N') && !detail::option_is(job, 'P') &&
		!detail::option_is(job, 'S') && !detail::option_is(job, 'B')) || (!rightv && !leftv)) {
		detail::tlapack_error("tgebak: invalid job/side");
	}
	if (n < 0 || m < 0 || ldv < std::max(1, n)) {
		detail::tlapack_error("tgebak: invalid argument");
	}
	if (n == 0 || m == 0 || detail::option_is(job, 'N')) {
		return 0;
	}
	const std::size_t LV = static_cast<std::size_t>(ldv);
	if (ilo != ihi && (detail::option_is(job, 'S') || detail::option_is(job, 'B'))) {
		if (rightv) {
			for (int i = ilo; i <= ihi; i++) {
				tscal(m, scale[i - 1], V + (i - 1), ldv);
			}
		}
		if (leftv) {
			for (int i = ilo; i <= ihi; i++) {
				tscal(m, T(1) / scale[i - 1], V + (i - 1), ldv);
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
				tswap(m, V + (i - 1), ldv, V + (kk - 1), ldv);
			}
		}
	}
	return 0;
}

// Hessenberg 化 (unblocked, ilo/ihi は 1-based)
template <typename T>
inline int tgehd2(const int n, const int ilo, const int ihi, T* A, const int lda,
	T* tau) {
	namespace detail = tlapack_detail;
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		lda < std::max(1, n)) {
		detail::tlapack_error("tgehd2: invalid argument");
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	for (int i = ilo; i <= ihi - 1; i++) {
		tlarfg(ihi - i, A[i + L * (i - 1)], A + std::min(i + 1, n - 1) + L * (i - 1), 1,
			tau[i - 1]);
		tlarf1f('R', ihi, ihi - i, A + i + L * (i - 1), 1, tau[i - 1], A + L * i, lda);
		tlarf1f('L', ihi - i, n - i, A + i + L * (i - 1), 1, tau[i - 1], A + i + L * i, lda);
	}
	return 0;
}

// Hessenberg 化の先頭 nb 段 (blocked 用 panel, k は 1-based)
template <typename T>
inline void tlahr2(const int n, const int k, const int nb, T* A, const int lda,
	T* tau, T* Tf, const int ldt, T* Y, const int ldy) {
	namespace detail = tlapack_detail;
	if (n <= 1) {
		return;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LT = static_cast<std::size_t>(ldt);
	const std::size_t LY = static_cast<std::size_t>(ldy);
	T ei = T(0);
	for (int i = 1; i <= nb; i++) {
		if (i > 1) {
			tgemv('N', n - k, i - 1, -T(1), Y + k, ldy, A + (k + i - 2), lda,
				T(1), A + k + L * (i - 1), 1);
			tcopy(i - 1, A + k + L * (i - 1), 1, Tf + LT * (nb - 1), 1);
			ttrmv('L', 'T', 'U', i - 1, A + k, lda, Tf + LT * (nb - 1), 1);
			tgemv('T', n - k - i + 1, i - 1, T(1), A + (k + i - 1), lda,
				A + (k + i - 1) + L * (i - 1), 1, T(1), Tf + LT * (nb - 1), 1);
			ttrmv('U', 'T', 'N', i - 1, Tf, ldt, Tf + LT * (nb - 1), 1);
			tgemv('N', n - k - i + 1, i - 1, -T(1), A + (k + i - 1), lda, Tf + LT * (nb - 1), 1,
				T(1), A + (k + i - 1) + L * (i - 1), 1);
			ttrmv('L', 'N', 'U', i - 1, A + k, lda, Tf + LT * (nb - 1), 1);
			taxpy(i - 1, -T(1), Tf + LT * (nb - 1), 1, A + k + L * (i - 1), 1);
			A[(k + i - 2) + L * (i - 2)] = ei;
		}
		tlarfg(n - k - i + 1, A[(k + i - 1) + L * (i - 1)],
			A + std::min(k + i, n - 1) + L * (i - 1), 1, tau[i - 1]);
		ei = A[(k + i - 1) + L * (i - 1)];
		A[(k + i - 1) + L * (i - 1)] = T(1);
		tgemv('N', n - k, n - k - i + 1, T(1), A + k + L * i, lda,
			A + (k + i - 1) + L * (i - 1), 1, T(0), Y + k + LY * (i - 1), 1);
		tgemv('T', n - k - i + 1, i - 1, T(1), A + (k + i - 1), lda,
			A + (k + i - 1) + L * (i - 1), 1, T(0), Tf + LT * (i - 1), 1);
		tgemv('N', n - k, i - 1, -T(1), Y + k, ldy, Tf + LT * (i - 1), 1,
			T(1), Y + k + LY * (i - 1), 1);
		tscal(n - k, tau[i - 1], Y + k + LY * (i - 1), 1);
		tscal(i - 1, -tau[i - 1], Tf + LT * (i - 1), 1);
		ttrmv('U', 'N', 'N', i - 1, Tf, ldt, Tf + LT * (i - 1), 1);
		Tf[(i - 1) + LT * (i - 1)] = tau[i - 1];
	}
	A[(k + nb - 1) + L * (nb - 1)] = ei;
	tlacpy('A', k, nb, A + L, lda, Y, ldy);
	ttrmm('R', 'L', 'N', 'U', k, nb, T(1), A + k, lda, Y, ldy);
	if (n > k + nb) {
		tgemm('N', 'N', k, nb, n - k - nb, T(1), A + L * (nb + 1), lda, A + (k + nb), lda,
			T(1), Y, ldy);
	}
	ttrmm('R', 'U', 'N', 'N', k, nb, T(1), Tf, ldt, Y, ldy);
}

// Hessenberg 化 (blocked, ilo/ihi は 1-based)
template <typename T>
inline int tgehrd(const int n, const int ilo, const int ihi, T* A, const int lda,
	T* tau) {
	namespace detail = tlapack_detail;
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		lda < std::max(1, n)) {
		detail::tlapack_error("tgehrd: invalid argument");
	}
	for (int i = 1; i <= ilo - 1; i++) {
		tau[i - 1] = T(0);
	}
	for (int i = std::max(1, ihi); i <= n - 1; i++) {
		tau[i - 1] = T(0);
	}
	const int nh = ihi - ilo + 1;
	if (nh <= 1) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = std::min(64, detail::ilaenv(1, "GEHRD", n, ilo, ihi, -1));
	int i = ilo;
	if (nb >= 2 && nb < nh) {
		const int nx = std::max(nb, detail::ilaenv(3, "GEHRD", n, ilo, ihi, -1));
		if (nx < nh) {
			const int ldwork = n;
			std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
			std::vector<T> work(static_cast<std::size_t>(ldwork) * nb, T(0));
			for (; i <= ihi - 1 - nx; i += nb) {
				const int ib = std::min(nb, ihi - i);
				tlahr2(ihi, i, ib, A + L * (i - 1), lda, tau + (i - 1), Tf.data(), nb,
					work.data(), ldwork);
				const T ei = A[(i + ib - 1) + L * (i + ib - 2)];
				A[(i + ib - 1) + L * (i + ib - 2)] = T(1);
				tgemm('N', 'T', ihi, ihi - i - ib + 1, ib, -T(1), work.data(), ldwork,
					A + (i + ib - 1) + L * (i - 1), lda, T(1), A + L * (i + ib - 1), lda);
				A[(i + ib - 1) + L * (i + ib - 2)] = ei;
				ttrmm('R', 'L', 'T', 'U', i, ib - 1, T(1), A + i + L * (i - 1), lda,
					work.data(), ldwork);
				for (int j = 0; j <= ib - 2; j++) {
					taxpy(i, -T(1), work.data() + static_cast<std::size_t>(ldwork) * j, 1,
						A + L * (i + j), 1);
				}
				tlarfb('L', 'T', 'F', 'C', ihi - i, n - i - ib + 1, ib, A + i + L * (i - 1), lda,
					Tf.data(), nb, A + i + L * (i + ib - 1), lda);
			}
		}
	}
	tgehd2(n, i, ihi, A, lda, tau);
	return 0;
}

// tgehrd の Q の生成 (ilo/ihi は 1-based)
template <typename T>
inline int torghr(const int n, const int ilo, const int ihi, T* A, const int lda,
	const T* tau) {
	namespace detail = tlapack_detail;
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		lda < std::max(1, n)) {
		detail::tlapack_error("torghr: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nh = ihi - ilo;
	// 列を 1 つ右へ shift して reflector を配置し，外側は単位行列にする
	for (int j = ihi; j >= ilo + 1; j--) {
		for (int i = 1; i <= j - 1; i++) {
			A[(i - 1) + L * (j - 1)] = T(0);
		}
		for (int i = j + 1; i <= ihi; i++) {
			A[(i - 1) + L * (j - 1)] = A[(i - 1) + L * (j - 2)];
		}
		for (int i = ihi + 1; i <= n; i++) {
			A[(i - 1) + L * (j - 1)] = T(0);
		}
	}
	for (int j = 1; j <= ilo; j++) {
		for (int i = 1; i <= n; i++) {
			A[(i - 1) + L * (j - 1)] = T(0);
		}
		A[(j - 1) + L * (j - 1)] = T(1);
	}
	for (int j = ihi + 1; j <= n; j++) {
		for (int i = 1; i <= n; i++) {
			A[(i - 1) + L * (j - 1)] = T(0);
		}
		A[(j - 1) + L * (j - 1)] = T(1);
	}
	if (nh > 0) {
		torgqr(nh, nh, nh, A + ilo + L * ilo, lda, tau + (ilo - 1));
	}
	return 0;
}

// Hessenberg 行列の T-shift QR (小規模用 algorithm を全 size に使用)
// ilo/ihi/iloz/ihiz は 1-based．返り値 info > 0 は非収束 (LAPACK と同じ意味)
template <typename T>
inline int tlahqr(const bool wantt, const bool wantz, const int n, const int ilo, const int ihi,
	T* H, const int ldh, T* wr, T* wi, const int iloz, const int ihiz,
	T* Z, const int ldz) {
	namespace detail = tlapack_detail;
	if (n == 0) {
		return 0;
	}
	const std::size_t LH = static_cast<std::size_t>(ldh);
	const std::size_t LZ = static_cast<std::size_t>(ldz);
	if (ilo == ihi) {
		wr[ilo - 1] = H[(ilo - 1) + LH * (ilo - 1)];
		wi[ilo - 1] = T(0);
		return 0;
	}
	for (int j = ilo; j <= ihi - 3; j++) {
		H[(j + 1) + LH * (j - 1)] = T(0);
		H[(j + 2) + LH * (j - 1)] = T(0);
	}
	if (ilo <= ihi - 2) {
		H[(ihi - 1) + LH * (ihi - 3)] = T(0);
	}
	const int nh = ihi - ilo + 1;
	const int nz = ihiz - iloz + 1;
	const T safmin = tlamch<T>('S');
	const T ulp = tlamch<T>('P');
	const T smlnum = safmin * (static_cast<T>(nh) / ulp);
	int i1 = 0, i2 = 0;
	if (wantt) {
		i1 = 1;
		i2 = n;
	}
	const int itmax = 30 * std::max(10, nh);
	int kdefl = 0;
	int i = ihi;
	T v[3];

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
				if (tlapack_detail::tabs(H[(k - 1) + LH * (k - 2)]) <= smlnum) {
					deflate = true;
					break;
				}
				T tst = tlapack_detail::tabs(H[(k - 2) + LH * (k - 2)]) + tlapack_detail::tabs(H[(k - 1) + LH * (k - 1)]);
				if (tst == T(0)) {
					if (k - 2 >= ilo) {
						tst = tst + tlapack_detail::tabs(H[(k - 2) + LH * (k - 3)]);
					}
					if (k + 1 <= ihi) {
						tst = tst + tlapack_detail::tabs(H[k + LH * (k - 1)]);
					}
				}
				if (tlapack_detail::tabs(H[(k - 1) + LH * (k - 2)]) <= ulp * tst) {
					const T ab = std::max(tlapack_detail::tabs(H[(k - 1) + LH * (k - 2)]),
						tlapack_detail::tabs(H[(k - 2) + LH * (k - 1)]));
					const T ba = std::min(tlapack_detail::tabs(H[(k - 1) + LH * (k - 2)]),
						tlapack_detail::tabs(H[(k - 2) + LH * (k - 1)]));
					const T aa = std::max(tlapack_detail::tabs(H[(k - 1) + LH * (k - 1)]),
						tlapack_detail::tabs(H[(k - 2) + LH * (k - 2)] - H[(k - 1) + LH * (k - 1)]));
					const T bb = std::min(tlapack_detail::tabs(H[(k - 1) + LH * (k - 1)]),
						tlapack_detail::tabs(H[(k - 2) + LH * (k - 2)] - H[(k - 1) + LH * (k - 1)]));
					const T s = aa + ab;
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
				H[(l - 1) + LH * (l - 2)] = T(0);
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
			T h11, h12, h21, h22;
			if (kdefl % 20 == 0) {
				// exceptional shift
				const T s = tlapack_detail::tabs(H[(i - 1) + LH * (i - 2)]) + tlapack_detail::tabs(H[(i - 2) + LH * (i - 3)]);
				h11 = (T(3) / T(4)) * s + H[(i - 1) + LH * (i - 1)];
				h12 = -(T(7) / T(16)) * s;
				h21 = s;
				h22 = h11;
			}
			else if (kdefl % 10 == 0) {
				const T s = tlapack_detail::tabs(H[l + LH * (l - 1)]) + tlapack_detail::tabs(H[(l + 1) + LH * l]);
				h11 = (T(3) / T(4)) * s + H[(l - 1) + LH * (l - 1)];
				h12 = -(T(7) / T(16)) * s;
				h21 = s;
				h22 = h11;
			}
			else {
				h11 = H[(i - 2) + LH * (i - 2)];
				h21 = H[(i - 1) + LH * (i - 2)];
				h12 = H[(i - 2) + LH * (i - 1)];
				h22 = H[(i - 1) + LH * (i - 1)];
			}
			T rt1r, rt1i, rt2r, rt2i;
			{
				const T s = tlapack_detail::tabs(h11) + tlapack_detail::tabs(h12) + tlapack_detail::tabs(h21) + tlapack_detail::tabs(h22);
				if (s == T(0)) {
					rt1r = rt1i = rt2r = rt2i = T(0);
				}
				else {
					h11 = h11 / s;
					h21 = h21 / s;
					h12 = h12 / s;
					h22 = h22 / s;
					const T tr = (h11 + h22) / T(2);
					const T det = (h11 - tr) * (h22 - tr) - h12 * h21;
					const T rtdisc = tlapack_detail::tsqrt(tlapack_detail::tabs(det));
					if (det >= T(0)) {
						rt1r = tr * s;
						rt2r = rt1r;
						rt1i = rtdisc * s;
						rt2i = -rt1i;
					}
					else {
						rt1r = tr + rtdisc;
						rt2r = tr - rtdisc;
						if (tlapack_detail::tabs(rt1r - h22) <= tlapack_detail::tabs(rt2r - h22)) {
							rt1r = rt1r * s;
							rt2r = rt1r;
						}
						else {
							rt2r = rt2r * s;
							rt1r = rt2r;
						}
						rt1i = T(0);
						rt2i = T(0);
					}
				}
			}
			// 最初の列を決める位置 m を探す
			for (m = i - 2; m >= l; m--) {
				T h21s = H[m + LH * (m - 1)];
				T s = tlapack_detail::tabs(H[(m - 1) + LH * (m - 1)] - rt2r) + tlapack_detail::tabs(rt2i) + tlapack_detail::tabs(h21s);
				h21s = H[m + LH * (m - 1)] / s;
				v[0] = h21s * H[(m - 1) + LH * m] + (H[(m - 1) + LH * (m - 1)] - rt1r) *
					((H[(m - 1) + LH * (m - 1)] - rt2r) / s) - rt1i * (rt2i / s);
				v[1] = h21s * (H[(m - 1) + LH * (m - 1)] + H[m + LH * m] - rt1r - rt2r);
				v[2] = h21s * H[(m + 1) + LH * m];
				s = tlapack_detail::tabs(v[0]) + tlapack_detail::tabs(v[1]) + tlapack_detail::tabs(v[2]);
				v[0] = v[0] / s;
				v[1] = v[1] / s;
				v[2] = v[2] / s;
				if (m == l) {
					break;
				}
				if (tlapack_detail::tabs(H[(m - 1) + LH * (m - 2)]) * (tlapack_detail::tabs(v[1]) + tlapack_detail::tabs(v[2])) <=
					ulp * tlapack_detail::tabs(v[0]) * (tlapack_detail::tabs(H[(m - 2) + LH * (m - 2)]) +
						tlapack_detail::tabs(H[(m - 1) + LH * (m - 1)]) + tlapack_detail::tabs(H[m + LH * m]))) {
					break;
				}
			}
			// T-shift QR sweep
			for (k = m; k <= i - 1; k++) {
				const int nr = std::min(3, i - k + 1);
				T t1;
				if (k > m) {
					tcopy(nr, H + (k - 1) + LH * (k - 2), 1, v, 1);
				}
				tlarfg(nr, v[0], v + 1, 1, t1);
				if (k > m) {
					H[(k - 1) + LH * (k - 2)] = v[0];
					H[k + LH * (k - 2)] = T(0);
					if (k < i - 1) {
						H[(k + 1) + LH * (k - 2)] = T(0);
					}
				}
				else if (m > l) {
					H[(k - 1) + LH * (k - 2)] = H[(k - 1) + LH * (k - 2)] * (T(1) - t1);
				}
				const T v2 = v[1];
				const T t2 = t1 * v2;
				if (nr == 3) {
					const T v3 = v[2];
					const T t3 = t1 * v3;
					for (int j = k; j <= i2; j++) {
						const T sum = H[(k - 1) + LH * (j - 1)] + v2 * H[k + LH * (j - 1)] +
							v3 * H[(k + 1) + LH * (j - 1)];
						H[(k - 1) + LH * (j - 1)] -= sum * t1;
						H[k + LH * (j - 1)] -= sum * t2;
						H[(k + 1) + LH * (j - 1)] -= sum * t3;
					}
					for (int j = i1; j <= std::min(k + 3, i); j++) {
						const T sum = H[(j - 1) + LH * (k - 1)] + v2 * H[(j - 1) + LH * k] +
							v3 * H[(j - 1) + LH * (k + 1)];
						H[(j - 1) + LH * (k - 1)] -= sum * t1;
						H[(j - 1) + LH * k] -= sum * t2;
						H[(j - 1) + LH * (k + 1)] -= sum * t3;
					}
					if (wantz) {
						for (int j = iloz; j <= ihiz; j++) {
							const T sum = Z[(j - 1) + LZ * (k - 1)] + v2 * Z[(j - 1) + LZ * k] +
								v3 * Z[(j - 1) + LZ * (k + 1)];
							Z[(j - 1) + LZ * (k - 1)] -= sum * t1;
							Z[(j - 1) + LZ * k] -= sum * t2;
							Z[(j - 1) + LZ * (k + 1)] -= sum * t3;
						}
					}
				}
				else if (nr == 2) {
					for (int j = k; j <= i2; j++) {
						const T sum = H[(k - 1) + LH * (j - 1)] + v2 * H[k + LH * (j - 1)];
						H[(k - 1) + LH * (j - 1)] -= sum * t1;
						H[k + LH * (j - 1)] -= sum * t2;
					}
					for (int j = i1; j <= i; j++) {
						const T sum = H[(j - 1) + LH * (k - 1)] + v2 * H[(j - 1) + LH * k];
						H[(j - 1) + LH * (k - 1)] -= sum * t1;
						H[(j - 1) + LH * k] -= sum * t2;
					}
					if (wantz) {
						for (int j = iloz; j <= ihiz; j++) {
							const T sum = Z[(j - 1) + LZ * (k - 1)] + v2 * Z[(j - 1) + LZ * k];
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
			wi[i - 1] = T(0);
		}
		else if (l == i - 1) {
			T cs, sn;
			tlanv2(H[(i - 2) + LH * (i - 2)], H[(i - 2) + LH * (i - 1)],
				H[(i - 1) + LH * (i - 2)], H[(i - 1) + LH * (i - 1)],
				wr[i - 2], wi[i - 2], wr[i - 1], wi[i - 1], cs, sn);
			if (wantt) {
				if (i2 > i) {
					trot(i2 - i, H + (i - 2) + LH * i, ldh, H + (i - 1) + LH * i, ldh, cs, sn);
				}
				trot(i - i1 - 1, H + (i1 - 1) + LH * (i - 2), 1, H + (i1 - 1) + LH * (i - 1), 1, cs, sn);
			}
			if (wantz) {
				trot(nz, Z + (iloz - 1) + LZ * (i - 2), 1, Z + (iloz - 1) + LZ * (i - 1), 1, cs, sn);
			}
		}
		kdefl = 0;
		i = l - 1;
		goto L20;
	}
}

// Hessenberg 行列の固有値 (と Schur 形) — tlahqr を全 size に使う簡略版
// job: 'E' 固有値のみ / 'S' Schur 形も, compz: 'N'/'I'/'V', ilo/ihi は 1-based
template <typename T>
inline int thseqr(const char job, const char compz, const int n, const int ilo, const int ihi,
	T* H, const int ldh, T* wr, T* wi, T* Z, const int ldz) {
	namespace detail = tlapack_detail;
	const bool wantt = detail::option_is(job, 'S');
	const bool initz = detail::option_is(compz, 'I');
	const bool wantz = initz || detail::option_is(compz, 'V');
	if (!wantt && !detail::option_is(job, 'E')) {
		detail::tlapack_error("thseqr: invalid job");
	}
	if (!wantz && !detail::option_is(compz, 'N')) {
		detail::tlapack_error("thseqr: invalid compz");
	}
	if (n < 0 || ilo < 1 || ilo > std::max(1, n) || ihi < std::min(ilo, n) || ihi > n ||
		ldh < std::max(1, n) || (wantz && ldz < std::max(1, n)) || (!wantz && ldz < 1)) {
		detail::tlapack_error("thseqr: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	if (initz) {
		tlaset('F', n, n, T(0), T(1), Z, ldz);
	}
	// ilo..ihi の外側の固有値は対角要素
	for (int i = 1; i <= ilo - 1; i++) {
		wr[i - 1] = H[(i - 1) + static_cast<std::size_t>(ldh) * (i - 1)];
		wi[i - 1] = T(0);
	}
	for (int i = ihi + 1; i <= n; i++) {
		wr[i - 1] = H[(i - 1) + static_cast<std::size_t>(ldh) * (i - 1)];
		wi[i - 1] = T(0);
	}
	return tlahqr(wantt, wantz, n, ilo, ihi, H, ldh, wr, wi, 1, n, Z, ldz);
}

// 1x1 / 2x2 の (ca*A - w*D) x = s*b を解く (固有 vector 計算用)
template <typename T>
inline int tlaln2(const bool ltrans, const int na, const int nw, const T& smin, const T& ca,
	const T* A, const int lda, const T& d1, const T& d2, const T* B, const int ldb,
	const T& wr, const T& wi, T* X, const int ldx, T& scale, T& xnorm) {
	namespace detail = tlapack_detail;
	const std::size_t LA = static_cast<std::size_t>(lda);
	const std::size_t LB = static_cast<std::size_t>(ldb);
	const std::size_t LX = static_cast<std::size_t>(ldx);
	static const bool zswap[4] = { false, false, true, true };
	static const bool rswap[4] = { false, true, false, true };
	static const int ipivot[4][4] = { { 1, 2, 3, 4 }, { 2, 1, 4, 3 }, { 3, 4, 1, 2 }, { 4, 3, 2, 1 } };
	const T smlnum = T(2) * tlamch<T>('S');
	const T bignum = T(1) / smlnum;
	const T smini = std::max(smin, smlnum);
	int info = 0;
	scale = T(1);
	if (na == 1) {
		if (nw == 1) {
			T csr = ca * A[0] - wr * d1;
			T cnorm = tlapack_detail::tabs(csr);
			if (cnorm < smini) {
				csr = smini;
				cnorm = smini;
				info = 1;
			}
			const T bnorm = tlapack_detail::tabs(B[0]);
			if (cnorm < T(1) && bnorm > T(1)) {
				if (bnorm > bignum * cnorm) {
					scale = T(1) / bnorm;
				}
			}
			X[0] = (B[0] * scale) / csr;
			xnorm = tlapack_detail::tabs(X[0]);
		}
		else {
			T csr = ca * A[0] - wr * d1;
			T csi = -wi * d1;
			T cnorm = tlapack_detail::tabs(csr) + tlapack_detail::tabs(csi);
			if (cnorm < smini) {
				csr = smini;
				csi = T(0);
				cnorm = smini;
				info = 1;
			}
			const T bnorm = tlapack_detail::tabs(B[0]) + tlapack_detail::tabs(B[LB]);
			if (cnorm < T(1) && bnorm > T(1)) {
				if (bnorm > bignum * cnorm) {
					scale = T(1) / bnorm;
				}
			}
			tladiv(scale * B[0], scale * B[LB], csr, csi, X[0], X[LX]);
			xnorm = tlapack_detail::tabs(X[0]) + tlapack_detail::tabs(X[LX]);
		}
	}
	else {
		// 2x2: crv = (CR(1,1), CR(2,1), CR(1,2), CR(2,2))
		T crv[4];
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
			T cmax = T(0);
			int icmax = 0;
			for (int j = 1; j <= 4; j++) {
				if (tlapack_detail::tabs(crv[j - 1]) > cmax) {
					cmax = tlapack_detail::tabs(crv[j - 1]);
					icmax = j;
				}
			}
			if (cmax < smini) {
				const T bnorm = std::max(tlapack_detail::tabs(B[0]), tlapack_detail::tabs(B[1]));
				if (smini < T(1) && bnorm > T(1)) {
					if (bnorm > bignum * smini) {
						scale = T(1) / bnorm;
					}
				}
				const T temp = scale / smini;
				X[0] = temp * B[0];
				X[1] = temp * B[1];
				xnorm = temp * bnorm;
				return 1;
			}
			const T ur11 = crv[icmax - 1];
			const T cr21 = crv[ipivot[icmax - 1][1] - 1];
			const T ur12 = crv[ipivot[icmax - 1][2] - 1];
			const T cr22 = crv[ipivot[icmax - 1][3] - 1];
			const T ur11r = T(1) / ur11;
			const T lr21 = ur11r * cr21;
			T ur22 = cr22 - ur12 * lr21;
			if (tlapack_detail::tabs(ur22) < smini) {
				ur22 = smini;
				info = 1;
			}
			T br1, br2;
			if (rswap[icmax - 1]) {
				br1 = B[1];
				br2 = B[0];
			}
			else {
				br1 = B[0];
				br2 = B[1];
			}
			br2 = br2 - lr21 * br1;
			const T bbnd = std::max(tlapack_detail::tabs(br1 * (ur22 * ur11r)), tlapack_detail::tabs(br2));
			if (bbnd > T(1) && tlapack_detail::tabs(ur22) < T(1)) {
				if (bbnd >= bignum * tlapack_detail::tabs(ur22)) {
					scale = T(1) / bbnd;
				}
			}
			const T xr2 = (br2 * scale) / ur22;
			const T xr1 = (scale * br1) * ur11r - xr2 * (ur11r * ur12);
			if (zswap[icmax - 1]) {
				X[0] = xr2;
				X[1] = xr1;
			}
			else {
				X[0] = xr1;
				X[1] = xr2;
			}
			xnorm = std::max(tlapack_detail::tabs(xr1), tlapack_detail::tabs(xr2));
			if (xnorm > T(1) && cmax > T(1)) {
				if (xnorm > bignum / cmax) {
					const T temp = cmax / bignum;
					X[0] = temp * X[0];
					X[1] = temp * X[1];
					xnorm = temp * xnorm;
					scale = temp * scale;
				}
			}
		}
		else {
			T civ[4];
			civ[0] = -wi * d1;
			civ[1] = T(0);
			civ[2] = T(0);
			civ[3] = -wi * d2;
			T cmax = T(0);
			int icmax = 0;
			for (int j = 1; j <= 4; j++) {
				if (tlapack_detail::tabs(crv[j - 1]) + tlapack_detail::tabs(civ[j - 1]) > cmax) {
					cmax = tlapack_detail::tabs(crv[j - 1]) + tlapack_detail::tabs(civ[j - 1]);
					icmax = j;
				}
			}
			if (cmax < smini) {
				const T bnorm = std::max(tlapack_detail::tabs(B[0]) + tlapack_detail::tabs(B[LB]),
					tlapack_detail::tabs(B[1]) + tlapack_detail::tabs(B[1 + LB]));
				if (smini < T(1) && bnorm > T(1)) {
					if (bnorm > bignum * smini) {
						scale = T(1) / bnorm;
					}
				}
				const T temp = scale / smini;
				X[0] = temp * B[0];
				X[1] = temp * B[1];
				X[LX] = temp * B[LB];
				X[1 + LX] = temp * B[1 + LB];
				xnorm = temp * bnorm;
				return 1;
			}
			const T ur11 = crv[icmax - 1];
			const T ui11 = civ[icmax - 1];
			const T cr21 = crv[ipivot[icmax - 1][1] - 1];
			const T ci21 = civ[ipivot[icmax - 1][1] - 1];
			const T ur12 = crv[ipivot[icmax - 1][2] - 1];
			const T ui12 = civ[ipivot[icmax - 1][2] - 1];
			const T cr22 = crv[ipivot[icmax - 1][3] - 1];
			const T ci22 = civ[ipivot[icmax - 1][3] - 1];
			T ur11r, ui11r, lr21, li21, ur12s, ui12s, ur22, ui22;
			if (icmax == 1 || icmax == 4) {
				if (tlapack_detail::tabs(ur11) > tlapack_detail::tabs(ui11)) {
					const T temp = ui11 / ur11;
					ur11r = T(1) / (ur11 * (T(1) + temp * temp));
					ui11r = -temp * ur11r;
				}
				else {
					const T temp = ur11 / ui11;
					ui11r = -T(1) / (ui11 * (T(1) + temp * temp));
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
				ur11r = T(1) / ur11;
				ui11r = T(0);
				lr21 = cr21 * ur11r;
				li21 = ci21 * ur11r;
				ur12s = ur12 * ur11r;
				ui12s = ui12 * ur11r;
				ur22 = cr22 - ur12 * lr21 + ui12 * li21;
				ui22 = -ur12 * li21 - ui12 * lr21;
			}
			const T u22abs = tlapack_detail::tabs(ur22) + tlapack_detail::tabs(ui22);
			if (u22abs < smini) {
				ur22 = smini;
				ui22 = T(0);
				info = 1;
			}
			T br1, br2, bi1, bi2;
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
			const T bbnd = std::max((tlapack_detail::tabs(br1) + tlapack_detail::tabs(bi1)) *
				(u22abs * (tlapack_detail::tabs(ur11r) + tlapack_detail::tabs(ui11r))),
				tlapack_detail::tabs(br2) + tlapack_detail::tabs(bi2));
			if (bbnd > T(1) && u22abs < T(1)) {
				if (bbnd >= bignum * u22abs) {
					scale = T(1) / bbnd;
					br1 = scale * br1;
					bi1 = scale * bi1;
					br2 = scale * br2;
					bi2 = scale * bi2;
				}
			}
			T xr2, xi2;
			tladiv(br2, bi2, ur22, ui22, xr2, xi2);
			const T xr1 = ur11r * br1 - ui11r * bi1 - ur12s * xr2 + ui12s * xi2;
			const T xi1 = ui11r * br1 + ur11r * bi1 - ui12s * xr2 - ur12s * xi2;
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
			xnorm = std::max(tlapack_detail::tabs(xr1) + tlapack_detail::tabs(xi1), tlapack_detail::tabs(xr2) + tlapack_detail::tabs(xi2));
			if (xnorm > T(1) && cmax > T(1)) {
				if (xnorm > bignum / cmax) {
					const T temp = cmax / bignum;
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

// 実 Schur 形 Tf の固有 vector (howmny: 'A' 全て / 'B' 全てを back-transform)
// side: 'R' 右 / 'L' 左 / 'B' 両方．SELECT 指定 ('S') は未対応．
template <typename T>
inline int ttrevc(const char side, const char howmny, const int n, const T* Tf, const int ldt,
	T* VL, const int ldvl, T* VR, const int ldvr, const int mm, int& m) {
	namespace detail = tlapack_detail;
	const bool bothv = detail::option_is(side, 'B');
	const bool rightv = detail::option_is(side, 'R') || bothv;
	const bool leftv = detail::option_is(side, 'L') || bothv;
	const bool over = detail::option_is(howmny, 'B');
	if ((!rightv && !leftv) || (!detail::option_is(howmny, 'A') && !over)) {
		detail::tlapack_error("ttrevc: invalid side/howmny (SELECT は未対応)");
	}
	m = n;
	if (n < 0 || ldt < std::max(1, n) || ldvl < 1 || (leftv && ldvl < n) ||
		ldvr < 1 || (rightv && ldvr < n) || mm < m) {
		detail::tlapack_error("ttrevc: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t LT = static_cast<std::size_t>(ldt);
	const std::size_t LVL = static_cast<std::size_t>(ldvl);
	const std::size_t LVR = static_cast<std::size_t>(ldvr);
	const T unfl = tlamch<T>('S');
	const T ulp = tlamch<T>('P');
	const T smlnum = unfl * (n / ulp);
	const T bignum = (T(1) - ulp) / smlnum;
	std::vector<T> work(static_cast<std::size_t>(3) * n + 1, T(0));
	T* wk = work.data(); // wk[X-1] が Fortran の WORK(X)
	T x[4]; // X(2,2), ldx=2
	wk[0] = T(0);
	for (int j = 2; j <= n; j++) {
		wk[j - 1] = T(0);
		for (int i = 1; i <= j - 1; i++) {
			wk[j - 1] += tlapack_detail::tabs(Tf[(i - 1) + LT * (j - 1)]);
		}
	}
	const int n2 = 2 * n;
	T scale2, xnorm;
	if (rightv) {
		int ip = 0;
		int is = m;
		for (int ki = n; ki >= 1; ki--) {
			bool skip = false;
			if (ip == 1) {
				skip = true;
			}
			else if (ki != 1 && Tf[(ki - 1) + LT * (ki - 2)] != T(0)) {
				ip = -1;
			}
			if (!skip) {
				const T wr = Tf[(ki - 1) + LT * (ki - 1)];
				T wi = T(0);
				if (ip != 0) {
					wi = tlapack_detail::tsqrt(tlapack_detail::tabs(Tf[(ki - 2) + LT * (ki - 1)])) *
						tlapack_detail::tsqrt(tlapack_detail::tabs(Tf[(ki - 1) + LT * (ki - 2)]));
				}
				const T smin = std::max(ulp * (tlapack_detail::tabs(wr) + tlapack_detail::tabs(wi)), smlnum);
				if (ip == 0) {
					// 実固有値
					wk[ki + n - 1] = T(1);
					for (int k = 1; k <= ki - 1; k++) {
						wk[k + n - 1] = -Tf[(k - 1) + LT * (ki - 1)];
					}
					int jnxt = ki - 1;
					for (int j = ki - 1; j >= 1; j--) {
						if (j > jnxt) {
							continue;
						}
						int j1 = j;
						jnxt = j - 1;
						if (j > 1 && Tf[(j - 1) + LT * (j - 2)] != T(0)) {
							j1 = j - 1;
							jnxt = j - 2;
						}
						if (j1 == j) {
							tlaln2(false, 1, 1, smin, T(1), Tf + (j - 1) + LT * (j - 1), ldt, T(1), T(1),
								wk + (j + n - 1), n, wr, T(0), x, 2, scale2, xnorm);
							if (xnorm > T(1) && wk[j - 1] > bignum / xnorm) {
								x[0] = x[0] / xnorm;
								scale2 = scale2 / xnorm;
							}
							if (scale2 != T(1)) {
								tscal(ki, scale2, wk + n, 1);
							}
							wk[j + n - 1] = x[0];
							taxpy(j - 1, -x[0], Tf + LT * (j - 1), 1, wk + n, 1);
						}
						else {
							tlaln2(false, 2, 1, smin, T(1), Tf + (j - 2) + LT * (j - 2), ldt, T(1), T(1),
								wk + (j + n - 2), n, wr, T(0), x, 2, scale2, xnorm);
							if (xnorm > T(1)) {
								const T beta = std::max(wk[j - 2], wk[j - 1]);
								if (beta > bignum / xnorm) {
									x[0] = x[0] / xnorm;
									x[1] = x[1] / xnorm;
									scale2 = scale2 / xnorm;
								}
							}
							if (scale2 != T(1)) {
								tscal(ki, scale2, wk + n, 1);
							}
							wk[j + n - 2] = x[0];
							wk[j + n - 1] = x[1];
							taxpy(j - 2, -x[0], Tf + LT * (j - 2), 1, wk + n, 1);
							taxpy(j - 2, -x[1], Tf + LT * (j - 1), 1, wk + n, 1);
						}
					}
					if (!over) {
						tcopy(ki, wk + n, 1, VR + LVR * (is - 1), 1);
						const int ii = itamax(ki, VR + LVR * (is - 1), 1) + 1;
						const T remax = T(1) / tlapack_detail::tabs(VR[(ii - 1) + LVR * (is - 1)]);
						tscal(ki, remax, VR + LVR * (is - 1), 1);
						for (int k = ki + 1; k <= n; k++) {
							VR[(k - 1) + LVR * (is - 1)] = T(0);
						}
					}
					else {
						if (ki > 1) {
							tgemv('N', n, ki - 1, T(1), VR, ldvr, wk + n, 1, wk[ki + n - 1],
								VR + LVR * (ki - 1), 1);
						}
						const int ii = itamax(n, VR + LVR * (ki - 1), 1) + 1;
						const T remax = T(1) / tlapack_detail::tabs(VR[(ii - 1) + LVR * (ki - 1)]);
						tscal(n, remax, VR + LVR * (ki - 1), 1);
					}
				}
				else {
					// 複素共役対
					if (tlapack_detail::tabs(Tf[(ki - 2) + LT * (ki - 1)]) >= tlapack_detail::tabs(Tf[(ki - 1) + LT * (ki - 2)])) {
						wk[ki + n - 2] = T(1);
						wk[ki + n2 - 1] = wi / Tf[(ki - 2) + LT * (ki - 1)];
					}
					else {
						wk[ki + n - 2] = -wi / Tf[(ki - 1) + LT * (ki - 2)];
						wk[ki + n2 - 1] = T(1);
					}
					wk[ki + n - 1] = T(0);
					wk[ki + n2 - 2] = T(0);
					for (int k = 1; k <= ki - 2; k++) {
						wk[k + n - 1] = -wk[ki + n - 2] * Tf[(k - 1) + LT * (ki - 2)];
						wk[k + n2 - 1] = -wk[ki + n2 - 1] * Tf[(k - 1) + LT * (ki - 1)];
					}
					int jnxt = ki - 2;
					for (int j = ki - 2; j >= 1; j--) {
						if (j > jnxt) {
							continue;
						}
						int j1 = j;
						jnxt = j - 1;
						if (j > 1 && Tf[(j - 1) + LT * (j - 2)] != T(0)) {
							j1 = j - 1;
							jnxt = j - 2;
						}
						if (j1 == j) {
							tlaln2(false, 1, 2, smin, T(1), Tf + (j - 1) + LT * (j - 1), ldt, T(1), T(1),
								wk + (j + n - 1), n, wr, wi, x, 2, scale2, xnorm);
							if (xnorm > T(1) && wk[j - 1] > bignum / xnorm) {
								x[0] = x[0] / xnorm;
								x[2] = x[2] / xnorm;
								scale2 = scale2 / xnorm;
							}
							if (scale2 != T(1)) {
								tscal(ki, scale2, wk + n, 1);
								tscal(ki, scale2, wk + n2, 1);
							}
							wk[j + n - 1] = x[0];
							wk[j + n2 - 1] = x[2];
							taxpy(j - 1, -x[0], Tf + LT * (j - 1), 1, wk + n, 1);
							taxpy(j - 1, -x[2], Tf + LT * (j - 1), 1, wk + n2, 1);
						}
						else {
							tlaln2(false, 2, 2, smin, T(1), Tf + (j - 2) + LT * (j - 2), ldt, T(1), T(1),
								wk + (j + n - 2), n, wr, wi, x, 2, scale2, xnorm);
							if (xnorm > T(1)) {
								const T beta = std::max(wk[j - 2], wk[j - 1]);
								if (beta > bignum / xnorm) {
									const T rec = T(1) / xnorm;
									x[0] *= rec;
									x[2] *= rec;
									x[1] *= rec;
									x[3] *= rec;
									scale2 *= rec;
								}
							}
							if (scale2 != T(1)) {
								tscal(ki, scale2, wk + n, 1);
								tscal(ki, scale2, wk + n2, 1);
							}
							wk[j + n - 2] = x[0];
							wk[j + n - 1] = x[1];
							wk[j + n2 - 2] = x[2];
							wk[j + n2 - 1] = x[3];
							taxpy(j - 2, -x[0], Tf + LT * (j - 2), 1, wk + n, 1);
							taxpy(j - 2, -x[1], Tf + LT * (j - 1), 1, wk + n, 1);
							taxpy(j - 2, -x[2], Tf + LT * (j - 2), 1, wk + n2, 1);
							taxpy(j - 2, -x[3], Tf + LT * (j - 1), 1, wk + n2, 1);
						}
					}
					if (!over) {
						tcopy(ki, wk + n, 1, VR + LVR * (is - 2), 1);
						tcopy(ki, wk + n2, 1, VR + LVR * (is - 1), 1);
						T emax = T(0);
						for (int k = 1; k <= ki; k++) {
							emax = std::max(emax, tlapack_detail::tabs(VR[(k - 1) + LVR * (is - 2)]) +
								tlapack_detail::tabs(VR[(k - 1) + LVR * (is - 1)]));
						}
						const T remax = T(1) / emax;
						tscal(ki, remax, VR + LVR * (is - 2), 1);
						tscal(ki, remax, VR + LVR * (is - 1), 1);
						for (int k = ki + 1; k <= n; k++) {
							VR[(k - 1) + LVR * (is - 2)] = T(0);
							VR[(k - 1) + LVR * (is - 1)] = T(0);
						}
					}
					else {
						if (ki > 2) {
							tgemv('N', n, ki - 2, T(1), VR, ldvr, wk + n, 1, wk[ki + n - 2],
								VR + LVR * (ki - 2), 1);
							tgemv('N', n, ki - 2, T(1), VR, ldvr, wk + n2, 1, wk[ki + n2 - 1],
								VR + LVR * (ki - 1), 1);
						}
						else {
							tscal(n, wk[ki + n - 2], VR + LVR * (ki - 2), 1);
							tscal(n, wk[ki + n2 - 1], VR + LVR * (ki - 1), 1);
						}
						T emax = T(0);
						for (int k = 1; k <= n; k++) {
							emax = std::max(emax, tlapack_detail::tabs(VR[(k - 1) + LVR * (ki - 2)]) +
								tlapack_detail::tabs(VR[(k - 1) + LVR * (ki - 1)]));
						}
						const T remax = T(1) / emax;
						tscal(n, remax, VR + LVR * (ki - 2), 1);
						tscal(n, remax, VR + LVR * (ki - 1), 1);
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
			else if (ki != n && Tf[ki + LT * (ki - 1)] != T(0)) {
				ip = 1;
			}
			if (!skip) {
				const T wr = Tf[(ki - 1) + LT * (ki - 1)];
				T wi = T(0);
				if (ip != 0) {
					wi = tlapack_detail::tsqrt(tlapack_detail::tabs(Tf[(ki - 1) + LT * ki])) *
						tlapack_detail::tsqrt(tlapack_detail::tabs(Tf[ki + LT * (ki - 1)]));
				}
				const T smin = std::max(ulp * (tlapack_detail::tabs(wr) + tlapack_detail::tabs(wi)), smlnum);
				if (ip == 0) {
					// 実固有値
					wk[ki + n - 1] = T(1);
					for (int k = ki + 1; k <= n; k++) {
						wk[k + n - 1] = -Tf[(ki - 1) + LT * (k - 1)];
					}
					T vmax = T(1);
					T vcrit = bignum;
					int jnxt = ki + 1;
					for (int j = ki + 1; j <= n; j++) {
						if (j < jnxt) {
							continue;
						}
						int j2 = j;
						jnxt = j + 1;
						if (j < n && Tf[j + LT * (j - 1)] != T(0)) {
							j2 = j + 1;
							jnxt = j + 2;
						}
						if (j2 == j) {
							if (wk[j - 1] > vcrit) {
								const T rec = T(1) / vmax;
								tscal(n - ki + 1, rec, wk + (ki + n - 1), 1);
								vmax = T(1);
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								tdot(j - ki - 1, Tf + ki + LT * (j - 1), 1, wk + (ki + n), 1);
							tlaln2(true, 1, 1, smin, T(1), Tf + (j - 1) + LT * (j - 1), ldt, T(1), T(1),
								wk + (j + n - 1), n, wr, T(0), x, 2, scale2, xnorm);
							if (scale2 != T(1)) {
								tscal(n - ki + 1, scale2, wk + (ki + n - 1), 1);
							}
							wk[j + n - 1] = x[0];
							vmax = std::max(tlapack_detail::tabs(wk[j + n - 1]), vmax);
							vcrit = bignum / vmax;
						}
						else {
							const T beta = std::max(wk[j - 1], wk[j]);
							if (beta > vcrit) {
								const T rec = T(1) / vmax;
								tscal(n - ki + 1, rec, wk + (ki + n - 1), 1);
								vmax = T(1);
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								tdot(j - ki - 1, Tf + ki + LT * (j - 1), 1, wk + (ki + n), 1);
							wk[j + n] = wk[j + n] -
								tdot(j - ki - 1, Tf + ki + LT * j, 1, wk + (ki + n), 1);
							tlaln2(true, 2, 1, smin, T(1), Tf + (j - 1) + LT * (j - 1), ldt, T(1), T(1),
								wk + (j + n - 1), n, wr, T(0), x, 2, scale2, xnorm);
							if (scale2 != T(1)) {
								tscal(n - ki + 1, scale2, wk + (ki + n - 1), 1);
							}
							wk[j + n - 1] = x[0];
							wk[j + n] = x[1];
							vmax = std::max(std::max(tlapack_detail::tabs(wk[j + n - 1]), tlapack_detail::tabs(wk[j + n])), vmax);
							vcrit = bignum / vmax;
						}
					}
					if (!over) {
						tcopy(n - ki + 1, wk + (ki + n - 1), 1, VL + (ki - 1) + LVL * (is - 1), 1);
						const int ii = itamax(n - ki + 1, VL + (ki - 1) + LVL * (is - 1), 1) + ki;
						const T remax = T(1) / tlapack_detail::tabs(VL[(ii - 1) + LVL * (is - 1)]);
						tscal(n - ki + 1, remax, VL + (ki - 1) + LVL * (is - 1), 1);
						for (int k = 1; k <= ki - 1; k++) {
							VL[(k - 1) + LVL * (is - 1)] = T(0);
						}
					}
					else {
						if (ki < n) {
							tgemv('N', n, n - ki, T(1), VL + LVL * ki, ldvl, wk + (ki + n), 1,
								wk[ki + n - 1], VL + LVL * (ki - 1), 1);
						}
						const int ii = itamax(n, VL + LVL * (ki - 1), 1) + 1;
						const T remax = T(1) / tlapack_detail::tabs(VL[(ii - 1) + LVL * (ki - 1)]);
						tscal(n, remax, VL + LVL * (ki - 1), 1);
					}
				}
				else {
					// 複素共役対
					if (tlapack_detail::tabs(Tf[(ki - 1) + LT * ki]) >= tlapack_detail::tabs(Tf[ki + LT * (ki - 1)])) {
						wk[ki + n - 1] = wi / Tf[(ki - 1) + LT * ki];
						wk[ki + n2] = T(1);
					}
					else {
						wk[ki + n - 1] = T(1);
						wk[ki + n2] = -wi / Tf[ki + LT * (ki - 1)];
					}
					wk[ki + n] = T(0);
					wk[ki + n2 - 1] = T(0);
					for (int k = ki + 2; k <= n; k++) {
						wk[k + n - 1] = -wk[ki + n - 1] * Tf[(ki - 1) + LT * (k - 1)];
						wk[k + n2 - 1] = -wk[ki + n2] * Tf[ki + LT * (k - 1)];
					}
					T vmax = T(1);
					T vcrit = bignum;
					int jnxt = ki + 2;
					for (int j = ki + 2; j <= n; j++) {
						if (j < jnxt) {
							continue;
						}
						int j2 = j;
						jnxt = j + 1;
						if (j < n && Tf[j + LT * (j - 1)] != T(0)) {
							j2 = j + 1;
							jnxt = j + 2;
						}
						if (j2 == j) {
							if (wk[j - 1] > vcrit) {
								const T rec = T(1) / vmax;
								tscal(n - ki + 1, rec, wk + (ki + n - 1), 1);
								tscal(n - ki + 1, rec, wk + (ki + n2 - 1), 1);
								vmax = T(1);
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								tdot(j - ki - 2, Tf + (ki + 1) + LT * (j - 1), 1, wk + (ki + n + 1), 1);
							wk[j + n2 - 1] = wk[j + n2 - 1] -
								tdot(j - ki - 2, Tf + (ki + 1) + LT * (j - 1), 1, wk + (ki + n2 + 1), 1);
							tlaln2(false, 1, 2, smin, T(1), Tf + (j - 1) + LT * (j - 1), ldt, T(1), T(1),
								wk + (j + n - 1), n, wr, -wi, x, 2, scale2, xnorm);
							if (scale2 != T(1)) {
								tscal(n - ki + 1, scale2, wk + (ki + n - 1), 1);
								tscal(n - ki + 1, scale2, wk + (ki + n2 - 1), 1);
							}
							wk[j + n - 1] = x[0];
							wk[j + n2 - 1] = x[2];
							vmax = std::max(std::max(tlapack_detail::tabs(wk[j + n - 1]), tlapack_detail::tabs(wk[j + n2 - 1])), vmax);
							vcrit = bignum / vmax;
						}
						else {
							const T beta = std::max(wk[j - 1], wk[j]);
							if (beta > vcrit) {
								const T rec = T(1) / vmax;
								tscal(n - ki + 1, rec, wk + (ki + n - 1), 1);
								tscal(n - ki + 1, rec, wk + (ki + n2 - 1), 1);
								vmax = T(1);
								vcrit = bignum;
							}
							wk[j + n - 1] = wk[j + n - 1] -
								tdot(j - ki - 2, Tf + (ki + 1) + LT * (j - 1), 1, wk + (ki + n + 1), 1);
							wk[j + n2 - 1] = wk[j + n2 - 1] -
								tdot(j - ki - 2, Tf + (ki + 1) + LT * (j - 1), 1, wk + (ki + n2 + 1), 1);
							wk[j + n] = wk[j + n] -
								tdot(j - ki - 2, Tf + (ki + 1) + LT * j, 1, wk + (ki + n + 1), 1);
							wk[j + n2] = wk[j + n2] -
								tdot(j - ki - 2, Tf + (ki + 1) + LT * j, 1, wk + (ki + n2 + 1), 1);
							tlaln2(true, 2, 2, smin, T(1), Tf + (j - 1) + LT * (j - 1), ldt, T(1), T(1),
								wk + (j + n - 1), n, wr, -wi, x, 2, scale2, xnorm);
							if (scale2 != T(1)) {
								tscal(n - ki + 1, scale2, wk + (ki + n - 1), 1);
								tscal(n - ki + 1, scale2, wk + (ki + n2 - 1), 1);
							}
							wk[j + n - 1] = x[0];
							wk[j + n2 - 1] = x[2];
							wk[j + n] = x[1];
							wk[j + n2] = x[3];
							vmax = std::max(std::max(tlapack_detail::tabs(x[0]), tlapack_detail::tabs(x[2])),
								std::max(std::max(tlapack_detail::tabs(x[1]), tlapack_detail::tabs(x[3])), vmax));
							vcrit = bignum / vmax;
						}
					}
					if (!over) {
						tcopy(n - ki + 1, wk + (ki + n - 1), 1, VL + (ki - 1) + LVL * (is - 1), 1);
						tcopy(n - ki + 1, wk + (ki + n2 - 1), 1, VL + (ki - 1) + LVL * is, 1);
						T emax = T(0);
						for (int k = ki; k <= n; k++) {
							emax = std::max(emax, tlapack_detail::tabs(VL[(k - 1) + LVL * (is - 1)]) +
								tlapack_detail::tabs(VL[(k - 1) + LVL * is]));
						}
						const T remax = T(1) / emax;
						tscal(n - ki + 1, remax, VL + (ki - 1) + LVL * (is - 1), 1);
						tscal(n - ki + 1, remax, VL + (ki - 1) + LVL * is, 1);
						for (int k = 1; k <= ki - 1; k++) {
							VL[(k - 1) + LVL * (is - 1)] = T(0);
							VL[(k - 1) + LVL * is] = T(0);
						}
					}
					else {
						if (ki < n - 1) {
							tgemv('N', n, n - ki - 1, T(1), VL + LVL * (ki + 1), ldvl,
								wk + (ki + n + 1), 1, wk[ki + n - 1], VL + LVL * (ki - 1), 1);
							tgemv('N', n, n - ki - 1, T(1), VL + LVL * (ki + 1), ldvl,
								wk + (ki + n2 + 1), 1, wk[ki + n2], VL + LVL * ki, 1);
						}
						else {
							tscal(n, wk[ki + n - 1], VL + LVL * (ki - 1), 1);
							tscal(n, wk[ki + n2], VL + LVL * ki, 1);
						}
						T emax = T(0);
						for (int k = 1; k <= n; k++) {
							emax = std::max(emax, tlapack_detail::tabs(VL[(k - 1) + LVL * (ki - 1)]) +
								tlapack_detail::tabs(VL[(k - 1) + LVL * ki]));
						}
						const T remax = T(1) / emax;
						tscal(n, remax, VL + LVL * (ki - 1), 1);
						tscal(n, remax, VL + LVL * ki, 1);
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
template <typename T>
inline int tgeev(const char jobvl, const char jobvr, const int n, T* A, const int lda,
	T* wr, T* wi, T* VL, const int ldvl, T* VR, const int ldvr) {
	namespace detail = tlapack_detail;
	const bool wantvl = detail::option_is(jobvl, 'V');
	const bool wantvr = detail::option_is(jobvr, 'V');
	if ((!wantvl && !detail::option_is(jobvl, 'N')) || (!wantvr && !detail::option_is(jobvr, 'N'))) {
		detail::tlapack_error("tgeev: invalid jobvl/jobvr");
	}
	if (n < 0 || lda < std::max(1, n) || ldvl < 1 || (wantvl && ldvl < n) ||
		ldvr < 1 || (wantvr && ldvr < n)) {
		detail::tlapack_error("tgeev: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t LVL = static_cast<std::size_t>(ldvl);
	const std::size_t LVR = static_cast<std::size_t>(ldvr);
	const T eps = tlamch<T>('P');
	T smlnum = tlamch<T>('S');
	smlnum = tlapack_detail::tsqrt(smlnum) / eps;
	const T bignum = T(1) / smlnum;
	const T anrm = tlange('M', n, n, A, lda);
	bool scalea = false;
	T cscale = T(1);
	if (anrm > T(0) && anrm < smlnum) {
		scalea = true;
		cscale = smlnum;
	}
	else if (anrm > bignum) {
		scalea = true;
		cscale = bignum;
	}
	if (scalea) {
		tlascl('G', 0, 0, anrm, cscale, n, n, A, lda);
	}
	std::vector<T> balscale(static_cast<std::size_t>(n), T(0));
	std::vector<T> tau(static_cast<std::size_t>(std::max(1, n - 1)), T(0));
	int ilo, ihi;
	tgebal('B', n, A, lda, ilo, ihi, balscale.data());
	tgehrd(n, ilo, ihi, A, lda, tau.data());
	int info;
	char side = 'N';
	if (wantvl) {
		side = 'L';
		tlacpy('L', n, n, A, lda, VL, ldvl);
		torghr(n, ilo, ihi, VL, ldvl, tau.data());
		info = thseqr('S', 'V', n, ilo, ihi, A, lda, wr, wi, VL, ldvl);
		if (wantvr) {
			side = 'B';
			tlacpy('F', n, n, VL, ldvl, VR, ldvr);
		}
	}
	else if (wantvr) {
		side = 'R';
		tlacpy('L', n, n, A, lda, VR, ldvr);
		torghr(n, ilo, ihi, VR, ldvr, tau.data());
		info = thseqr('S', 'V', n, ilo, ihi, A, lda, wr, wi, VR, ldvr);
	}
	else {
		info = thseqr('E', 'N', n, ilo, ihi, A, lda, wr, wi, VR, ldvr);
	}
	if (info == 0) {
		if (wantvl || wantvr) {
			int nout;
			ttrevc(side, 'B', n, A, lda, VL, ldvl, VR, ldvr, n, nout);
		}
		std::vector<T> tmp(static_cast<std::size_t>(n), T(0));
		for (int pass = 0; pass < 2; pass++) {
			T* V = (pass == 0) ? VL : VR;
			const std::size_t LV = (pass == 0) ? LVL : LVR;
			const int ldv = (pass == 0) ? ldvl : ldvr;
			if ((pass == 0 && !wantvl) || (pass == 1 && !wantvr)) {
				continue;
			}
			tgebak('B', (pass == 0) ? 'L' : 'R', n, ilo, ihi, balscale.data(), n, V, ldv);
			for (int i = 1; i <= n; i++) {
				if (wi[i - 1] == T(0)) {
					const T scl = T(1) / tnrm2(n, V + LV * (i - 1), 1);
					tscal(n, scl, V + LV * (i - 1), 1);
				}
				else if (wi[i - 1] > T(0)) {
					const T scl = T(1) / tlapy2(tnrm2(n, V + LV * (i - 1), 1),
						tnrm2(n, V + LV * i, 1));
					tscal(n, scl, V + LV * (i - 1), 1);
					tscal(n, scl, V + LV * i, 1);
					for (int k = 1; k <= n; k++) {
						tmp[k - 1] = V[(k - 1) + LV * (i - 1)] * V[(k - 1) + LV * (i - 1)] +
							V[(k - 1) + LV * i] * V[(k - 1) + LV * i];
					}
					const int k = itamax(n, tmp.data(), 1) + 1;
					T cs, sn, r;
					tlartg(V[(k - 1) + LV * (i - 1)], V[(k - 1) + LV * i], cs, sn, r);
					trot(n, V + LV * (i - 1), 1, V + LV * i, 1, cs, sn);
					V[(k - 1) + LV * i] = T(0);
				}
			}
		}
	}
	if (scalea) {
		tlascl('G', 0, 0, cscale, anrm, n - info, 1, wr + info, std::max(n - info, 1));
		tlascl('G', 0, 0, cscale, anrm, n - info, 1, wi + info, std::max(n - info, 1));
		if (info > 0) {
			tlascl('G', 0, 0, cscale, anrm, ilo - 1, 1, wr, n);
			tlascl('G', 0, 0, cscale, anrm, ilo - 1, 1, wi, n);
		}
	}
	return info;
}

#endif // TLAPACK_TLAPACK_GEEV_HPP
