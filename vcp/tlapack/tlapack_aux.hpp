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

#ifndef TLAPACK_TLAPACK_AUX_HPP
#define TLAPACK_TLAPACK_AUX_HPP

// LAPACK 補助 routine (template 版)．reference LAPACK 3.12.1 からの移植．
// 関数名は LAPACK の d を t に読み替えたもの (丸めモードは扱わない)．
// FP 演算を含まない routine (tlaset, tlacpy, tlaswp, tlasrt, ilatlr, ilatlc) は
// LAPACK の d を t に読み替えた同名 routine．index は 0-based．

#include <algorithm>
#include <functional>

#include "tlapack_common.hpp"

namespace vcp {

// ---- 演算なし ----

// 末尾の zero 行を除いた行数 (A の非 zero 要素を含む最終行 + 1)
template <typename T>
inline int ilatlr(const int m, const int n, const T* A, const int lda) {
	if (m == 0) {
		return 0;
	}
	if (A[(m - 1)] != T(0) || A[(m - 1) + static_cast<std::size_t>(lda) * (n - 1)] != T(0)) {
		return m;
	}
	int result = 0;
	for (int j = 0; j < n; j++) {
		int i = m;
		while (i >= 1 && A[(i - 1) + static_cast<std::size_t>(lda) * j] == T(0)) {
			i--;
		}
		result = std::max(result, i);
	}
	return result;
}

// 末尾の zero 列を除いた列数
template <typename T>
inline int ilatlc(const int m, const int n, const T* A, const int lda) {
	if (n == 0) {
		return 0;
	}
	if (A[static_cast<std::size_t>(lda) * (n - 1)] != T(0) ||
		A[(m - 1) + static_cast<std::size_t>(lda) * (n - 1)] != T(0)) {
		return n;
	}
	for (int j = n; j >= 1; j--) {
		for (int i = 0; i < m; i++) {
			if (A[i + static_cast<std::size_t>(lda) * (j - 1)] != T(0)) {
				return j;
			}
		}
	}
	return 0;
}

// A の off-diagonal を alpha，diagonal を beta に設定 (uplo: 'U' 上三角のみ /
// 'L' 下三角のみ / その他 全体)
template <typename T>
inline void tlaset(const char uplo, const int m, const int n, const T& alpha, const T& beta, T* A, const int lda) {
	namespace detail = tlapack_detail;
	if (detail::option_is(uplo, 'U')) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < std::min(j, m); i++) {
				A[i + static_cast<std::size_t>(lda) * j] = alpha;
			}
		}
	}
	else if (detail::option_is(uplo, 'L')) {
		for (int j = 0; j < std::min(m, n); j++) {
			for (int i = j + 1; i < m; i++) {
				A[i + static_cast<std::size_t>(lda) * j] = alpha;
			}
		}
	}
	else {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				A[i + static_cast<std::size_t>(lda) * j] = alpha;
			}
		}
	}
	for (int i = 0; i < std::min(m, n); i++) {
		A[i + static_cast<std::size_t>(lda) * i] = beta;
	}
}

// B := A (uplo: 'U' 上三角 / 'L' 下三角 / その他 全体)
template <typename T>
inline void tlacpy(const char uplo, const int m, const int n, const T* A, const int lda, T* B, const int ldb) {
	namespace detail = tlapack_detail;
	if (detail::option_is(uplo, 'U')) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i <= std::min(j, m - 1); i++) {
				B[i + static_cast<std::size_t>(ldb) * j] = A[i + static_cast<std::size_t>(lda) * j];
			}
		}
	}
	else if (detail::option_is(uplo, 'L')) {
		for (int j = 0; j < n; j++) {
			for (int i = j; i < m; i++) {
				B[i + static_cast<std::size_t>(ldb) * j] = A[i + static_cast<std::size_t>(lda) * j];
			}
		}
	}
	else {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				B[i + static_cast<std::size_t>(ldb) * j] = A[i + static_cast<std::size_t>(lda) * j];
			}
		}
	}
}

// 行交換: i = k1, ..., k2-1 (0-based 半開区間) について A の行 i と行 ipiv[i] を
// 交換する (incx = 1: 昇順, incx = -1: 降順)．ipiv は 0-based．
template <typename T>
inline void tlaswp(const int n, T* A, const int lda, const int k1, const int k2, const int* ipiv, const int incx) {
	if (n <= 0) {
		return;
	}
	if (incx > 0) {
		for (int i = k1; i < k2; i++) {
			const int ip = ipiv[i];
			if (ip != i) {
				for (int j = 0; j < n; j++) {
					const std::size_t c = static_cast<std::size_t>(lda) * j;
					const T tmp = A[i + c];
					A[i + c] = A[ip + c];
					A[ip + c] = tmp;
				}
			}
		}
	}
	else if (incx < 0) {
		for (int i = k2 - 1; i >= k1; i--) {
			const int ip = ipiv[i];
			if (ip != i) {
				for (int j = 0; j < n; j++) {
					const std::size_t c = static_cast<std::size_t>(lda) * j;
					const T tmp = A[i + c];
					A[i + c] = A[ip + c];
					A[ip + c] = tmp;
				}
			}
		}
	}
}

// d を昇順 ('I') または降順 ('D') に sort (比較のみで丸め誤差なし)
template <typename T>
inline void tlasrt(const char id, const int n, T* d) {
	namespace detail = tlapack_detail;
	if (n < 0) {
		detail::tlapack_error("tlasrt: n < 0");
	}
	if (detail::option_is(id, 'I')) {
		std::sort(d, d + n);
	}
	else if (detail::option_is(id, 'D')) {
		std::sort(d, d + n, std::greater<T>());
	}
	else {
		detail::tlapack_error("tlasrt: invalid id");
	}
}

// ---- 演算あり ----

// scaled sum of squares: scale^2*sumsq += sum x_i^2 (overflow/underflow 回避付き)
// reference LAPACK 3.12 の Blue's algorithm は 2^k 定数 (2^-511 など) が型に
// 依存するため使わず，古典的な scale/ssq 更新 (LAPACK 3.8 以前の DLASSQ と
// 同方式) で実装する．scale^2*sumsq の不変量と呼び出し規約は同じ．
template <typename T>
inline void tlassq(const int n, const T* x, const int incx, T& scale, T& sumsq) {
	if (tlapack_detail::tisnan(scale) || tlapack_detail::tisnan(sumsq)) {
		return;
	}
	if (sumsq == T(0)) {
		scale = T(1);
	}
	if (scale == T(0)) {
		scale = T(1);
		sumsq = T(0);
	}
	if (n <= 0) {
		return;
	}
	int ix = 0;
	if (incx < 0) {
		ix = -(n - 1) * incx;
	}
	for (int i = 0; i < n; i++) {
		const T ax = tlapack_detail::tabs(x[ix]);
		if (ax != T(0)) { // NaN も != が true になりこの枝で伝播する
			if (scale < ax) {
				const T r = scale / ax;
				sumsq = T(1) + sumsq * (r * r);
				scale = ax;
			}
			else {
				const T r = ax / scale;
				sumsq = sumsq + r * r;
			}
		}
		ix += incx;
	}
}

// sqrt(x^2 + y^2) (overflow 回避付き)
template <typename T>
inline T tlapy2(const T& x, const T& y) {
	namespace detail = tlapack_detail;
	if (tlapack_detail::tisnan(x)) {
		return x;
	}
	if (tlapack_detail::tisnan(y)) {
		return y;
	}
	const T xabs = tlapack_detail::tabs(x);
	const T yabs = tlapack_detail::tabs(y);
	const T w = std::max(xabs, yabs);
	const T z = std::min(xabs, yabs);
	if (z == T(0) || w > tlamch<T>('O')) {
		return w;
	}
	return w * tlapack_detail::tsqrt(T(1) + (z / w) * (z / w));
}

// sqrt(x^2 + y^2 + z^2) (overflow 回避付き)
template <typename T>
inline T tlapy3(const T& x, const T& y, const T& z) {
	namespace detail = tlapack_detail;
	const T xabs = tlapack_detail::tabs(x);
	const T yabs = tlapack_detail::tabs(y);
	const T zabs = tlapack_detail::tabs(z);
	const T w = std::max(std::max(xabs, yabs), zabs);
	if (w == T(0) || w > tlamch<T>('O')) {
		return xabs + yabs + zabs;
	}
	return w * tlapack_detail::tsqrt((xabs / w) * (xabs / w) + (yabs / w) * (yabs / w) + (zabs / w) * (zabs / w));
}

} // namespace vcp

namespace tlapack_detail {

template <typename T>
inline T tladiv2(const T& a, const T& b, const T& c, const T& d, const T& r, const T& t) {
	if (r != T(0)) {
		const T br = b * r;
		if (br != T(0)) {
			return (a + br) * t;
		}
		return a * t + (b * t) * r;
	}
	return (a + d * (b / c)) * t;
}

template <typename T>
inline void tladiv1(const T& a, const T& b, const T& c, const T& d, T& p, T& q) {
	const T r = d / c;
	const T t = T(1) / (c + d * r);
	p = tladiv2(a, b, c, d, r, t);
	q = tladiv2(b, -a, c, d, r, t);
}

} // namespace tlapack_detail

namespace vcp {

// 複素数除算 p + iq := (a + ib) / (c + id) (overflow 回避付き, Baudin-Smith)
template <typename T>
inline void tladiv(const T& a, const T& b, const T& c, const T& d, T& p, T& q) {
	namespace detail = tlapack_detail;
	T aa = a;
	T bb = b;
	T cc = c;
	T dd = d;
	const T ab = std::max(tlapack_detail::tabs(a), tlapack_detail::tabs(b));
	const T cd = std::max(tlapack_detail::tabs(c), tlapack_detail::tabs(d));
	T s = T(1);
	const T ov = tlamch<T>('O');
	const T un = tlamch<T>('S');
	const T eps = tlamch<T>('E');
	const T be = T(2) / (eps * eps);
	if (ab >= (T(1) / T(2)) * ov) {
		aa = (T(1) / T(2)) * aa;
		bb = (T(1) / T(2)) * bb;
		s = T(2) * s;
	}
	if (cd >= (T(1) / T(2)) * ov) {
		cc = (T(1) / T(2)) * cc;
		dd = (T(1) / T(2)) * dd;
		s = (T(1) / T(2)) * s;
	}
	if (ab <= un * T(2) / eps) {
		aa = aa * be;
		bb = bb * be;
		s = s / be;
	}
	if (cd <= un * T(2) / eps) {
		cc = cc * be;
		dd = dd * be;
		s = s * be;
	}
	if (tlapack_detail::tabs(d) <= tlapack_detail::tabs(c)) {
		detail::tladiv1(aa, bb, cc, dd, p, q);
	}
	else {
		detail::tladiv1(bb, aa, dd, cc, p, q);
		q = -q;
	}
	p = p * s;
	q = q * s;
}

// Givens 回転の生成: [c s; -s c]^T [f; g] = [r; 0]
template <typename T>
inline void tlartg(const T& f, const T& g, T& cs, T& sn, T& r) {
	namespace detail = tlapack_detail;
	const T safmin = tlamch<T>('S');
	const T safmax = T(1) / tlamch<T>('S');
	const T rtmin = tlapack_detail::tsqrt(safmin);
	const T rtmax = tlapack_detail::tsqrt(safmax / T(2));
	const T f1 = tlapack_detail::tabs(f);
	const T g1 = tlapack_detail::tabs(g);
	if (g == T(0)) {
		cs = T(1);
		sn = T(0);
		r = f;
	}
	else if (f == T(0)) {
		cs = T(0);
		sn = detail::f_sign(T(1), g);
		r = g1;
	}
	else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
		const T d = tlapack_detail::tsqrt(f * f + g * g);
		cs = f1 / d;
		r = detail::f_sign(d, f);
		sn = g / r;
	}
	else {
		const T u = std::min(safmax, std::max(std::max(safmin, f1), g1));
		const T fs = f / u;
		const T gs = g / u;
		const T d = tlapack_detail::tsqrt(fs * fs + gs * gs);
		cs = tlapack_detail::tabs(fs) / d;
		r = detail::f_sign(d, f);
		sn = gs / r;
		r = r * u;
	}
}

// 対称 2x2 行列 [a b; b c] の固有値 (rt1 >= rt2)
template <typename T>
inline void tlae2(const T& a, const T& b, const T& c, T& rt1, T& rt2) {
	namespace detail = tlapack_detail;
	const T sm = a + c;
	const T df = a - c;
	const T adf = tlapack_detail::tabs(df);
	const T tb = b + b;
	const T ab = tlapack_detail::tabs(tb);
	T acmx, acmn;
	if (tlapack_detail::tabs(a) > tlapack_detail::tabs(c)) {
		acmx = a;
		acmn = c;
	}
	else {
		acmx = c;
		acmn = a;
	}
	T rt;
	if (adf > ab) {
		rt = adf * tlapack_detail::tsqrt(T(1) + (ab / adf) * (ab / adf));
	}
	else if (adf < ab) {
		rt = ab * tlapack_detail::tsqrt(T(1) + (adf / ab) * (adf / ab));
	}
	else {
		rt = ab * tlapack_detail::tsqrt(T(2));
	}
	if (sm < T(0)) {
		rt1 = (T(1) / T(2)) * (sm - rt);
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else if (sm > T(0)) {
		rt1 = (T(1) / T(2)) * (sm + rt);
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else {
		rt1 = (T(1) / T(2)) * rt;
		rt2 = -(T(1) / T(2)) * rt;
	}
}

// 対称 2x2 行列 [a b; b c] の固有値と固有 vector (回転 cs1, sn1)
template <typename T>
inline void tlaev2(const T& a, const T& b, const T& c, T& rt1, T& rt2, T& cs1, T& sn1) {
	namespace detail = tlapack_detail;
	const T sm = a + c;
	const T df = a - c;
	const T adf = tlapack_detail::tabs(df);
	const T tb = b + b;
	const T ab = tlapack_detail::tabs(tb);
	T acmx, acmn;
	if (tlapack_detail::tabs(a) > tlapack_detail::tabs(c)) {
		acmx = a;
		acmn = c;
	}
	else {
		acmx = c;
		acmn = a;
	}
	T rt;
	if (adf > ab) {
		rt = adf * tlapack_detail::tsqrt(T(1) + (ab / adf) * (ab / adf));
	}
	else if (adf < ab) {
		rt = ab * tlapack_detail::tsqrt(T(1) + (adf / ab) * (adf / ab));
	}
	else {
		rt = ab * tlapack_detail::tsqrt(T(2));
	}
	int sgn1;
	if (sm < T(0)) {
		rt1 = (T(1) / T(2)) * (sm - rt);
		sgn1 = -1;
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else if (sm > T(0)) {
		rt1 = (T(1) / T(2)) * (sm + rt);
		sgn1 = 1;
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else {
		rt1 = (T(1) / T(2)) * rt;
		rt2 = -(T(1) / T(2)) * rt;
		sgn1 = 1;
	}
	T cs;
	int sgn2;
	if (df >= T(0)) {
		cs = df + rt;
		sgn2 = 1;
	}
	else {
		cs = df - rt;
		sgn2 = -1;
	}
	const T acs = tlapack_detail::tabs(cs);
	if (acs > ab) {
		const T ct = -tb / cs;
		sn1 = T(1) / tlapack_detail::tsqrt(T(1) + ct * ct);
		cs1 = ct * sn1;
	}
	else {
		if (ab == T(0)) {
			cs1 = T(1);
			sn1 = T(0);
		}
		else {
			const T tn = -cs / tb;
			cs1 = T(1) / tlapack_detail::tsqrt(T(1) + tn * tn);
			sn1 = tn * cs1;
		}
	}
	if (sgn1 == sgn2) {
		const T tn = cs1;
		cs1 = -sn1;
		sn1 = tn;
	}
}

// 2x2 上三角行列 [f g; 0 h] の特異値 (ssmin <= ssmax)
template <typename T>
inline void tlas2(const T& f, const T& g, const T& h, T& ssmin, T& ssmax) {
	namespace detail = tlapack_detail;
	const T fa = tlapack_detail::tabs(f);
	const T ga = tlapack_detail::tabs(g);
	const T ha = tlapack_detail::tabs(h);
	const T fhmn = std::min(fa, ha);
	const T fhmx = std::max(fa, ha);
	if (fhmn == T(0)) {
		ssmin = T(0);
		if (fhmx == T(0)) {
			ssmax = ga;
		}
		else {
			const T mn = std::min(fhmx, ga);
			const T mx = std::max(fhmx, ga);
			ssmax = mx * tlapack_detail::tsqrt(T(1) + (mn / mx) * (mn / mx));
		}
	}
	else {
		if (ga < fhmx) {
			const T as = T(1) + fhmn / fhmx;
			const T at = (fhmx - fhmn) / fhmx;
			const T au = (ga / fhmx) * (ga / fhmx);
			const T c = T(2) / (tlapack_detail::tsqrt(as * as + au) + tlapack_detail::tsqrt(at * at + au));
			ssmin = fhmn * c;
			ssmax = fhmx / c;
		}
		else {
			const T au = fhmx / ga;
			if (au == T(0)) {
				ssmin = (fhmn * fhmx) / ga;
				ssmax = ga;
			}
			else {
				const T as = T(1) + fhmn / fhmx;
				const T at = (fhmx - fhmn) / fhmx;
				const T c = T(1) / (tlapack_detail::tsqrt(T(1) + (as * au) * (as * au)) +
					tlapack_detail::tsqrt(T(1) + (at * au) * (at * au)));
				ssmin = (fhmn * c) * au;
				ssmin = ssmin + ssmin;
				ssmax = ga / (c + c);
			}
		}
	}
}

// 2x2 上三角行列 [f g; 0 h] の SVD (左右の回転付き)
template <typename T>
inline void tlasv2(const T& f, const T& g, const T& h,
	T& ssmin, T& ssmax, T& snr, T& csr, T& snl, T& csl) {
	namespace detail = tlapack_detail;
	T ft = f;
	T fa = tlapack_detail::tabs(ft);
	T ht = h;
	T ha = tlapack_detail::tabs(h);
	int pmax = 1;
	const bool swap = (ha > fa);
	if (swap) {
		pmax = 3;
		T temp = ft;
		ft = ht;
		ht = temp;
		temp = fa;
		fa = ha;
		ha = temp;
	}
	const T gt = g;
	const T ga = tlapack_detail::tabs(gt);
	T clt = T(0), crt = T(0), slt = T(0), srt = T(0);
	if (ga == T(0)) {
		ssmin = ha;
		ssmax = fa;
		clt = T(1);
		crt = T(1);
		slt = T(0);
		srt = T(0);
	}
	else {
		bool gasmal = true;
		if (ga > fa) {
			pmax = 2;
			if ((fa / ga) < (tlamch<T>('E'))) {
				gasmal = false;
				ssmax = ga;
				if (ha > T(1)) {
					ssmin = fa / (ga / ha);
				}
				else {
					ssmin = (fa / ga) * ha;
				}
				clt = T(1);
				slt = ht / gt;
				srt = T(1);
				crt = ft / gt;
			}
		}
		if (gasmal) {
			const T d = fa - ha;
			T l;
			if (d == fa) {
				l = T(1);
			}
			else {
				l = d / fa;
			}
			const T m = gt / ft;
			T t = T(2) - l;
			const T mm = m * m;
			const T tt = t * t;
			const T s = tlapack_detail::tsqrt(tt + mm);
			T r;
			if (l == T(0)) {
				r = tlapack_detail::tabs(m);
			}
			else {
				r = tlapack_detail::tsqrt(l * l + mm);
			}
			const T a = (T(1) / T(2)) * (s + r);
			ssmin = ha / a;
			ssmax = fa * a;
			if (mm == T(0)) {
				if (l == T(0)) {
					t = detail::f_sign(T(2), ft) * detail::f_sign(T(1), gt);
				}
				else {
					t = gt / detail::f_sign(d, ft) + m / t;
				}
			}
			else {
				t = (m / (s + t) + m / (r + l)) * (T(1) + a);
			}
			const T ll = tlapack_detail::tsqrt(t * t + T(4));
			crt = T(2) / ll;
			srt = t / ll;
			clt = (crt + srt * m) / a;
			slt = (ht / ft) * srt / a;
		}
	}
	if (swap) {
		csl = srt;
		snl = crt;
		csr = slt;
		snr = clt;
	}
	else {
		csl = clt;
		snl = slt;
		csr = crt;
		snr = srt;
	}
	T tsign = T(0);
	if (pmax == 1) {
		tsign = detail::f_sign(T(1), csr) * detail::f_sign(T(1), csl) * detail::f_sign(T(1), f);
	}
	if (pmax == 2) {
		tsign = detail::f_sign(T(1), snr) * detail::f_sign(T(1), csl) * detail::f_sign(T(1), g);
	}
	if (pmax == 3) {
		tsign = detail::f_sign(T(1), snr) * detail::f_sign(T(1), snl) * detail::f_sign(T(1), h);
	}
	ssmax = detail::f_sign(ssmax, tsign);
	ssmin = detail::f_sign(ssmin, tsign * detail::f_sign(T(1), f) * detail::f_sign(T(1), h));
}

// 2x2 行列 [a b; c d] を実 Schur 標準形へ変換する回転 (cs, sn) と固有値
template <typename T>
inline void tlanv2(T& a, T& b, T& c, T& d,
	T& rt1r, T& rt1i, T& rt2r, T& rt2i, T& cs, T& sn) {
	namespace detail = tlapack_detail;
	const T eps = tlamch<T>('P');
	// reference LAPACK の safmn2 = 2^int(log(safmin/eps)/log(2)/2) (T では 2^-485)
	// に対応する scaling 閾値．2^k への丸めは型に依存するため，値そのもの
	// sqrt(safmin/eps) を使う (閾値・scaling 係数としての役割は同じで，
	// 極端な大きさの入力での scaling 時に丸めが 1 回増えるだけ)
	const T safmn2 = tlapack_detail::tsqrt(tlamch<T>('S') / eps);
	const T safmx2 = T(1) / safmn2;
	if (c == T(0)) {
		cs = T(1);
		sn = T(0);
	}
	else if (b == T(0)) {
		cs = T(0);
		sn = T(1);
		const T temp = d;
		d = a;
		a = temp;
		b = -c;
		c = T(0);
	}
	else if ((a - d) == T(0) && detail::f_sign(T(1), b) != detail::f_sign(T(1), c)) {
		cs = T(1);
		sn = T(0);
	}
	else {
		T temp = a - d;
		T p = (T(1) / T(2)) * temp;
		const T bcmax = std::max(tlapack_detail::tabs(b), tlapack_detail::tabs(c));
		const T bcmis = std::min(tlapack_detail::tabs(b), tlapack_detail::tabs(c)) * detail::f_sign(T(1), b) * detail::f_sign(T(1), c);
		T scale = std::max(tlapack_detail::tabs(p), bcmax);
		T z = (p / scale) * p + (bcmax / scale) * bcmis;
		if (z >= T(4) * eps) {
			z = p + detail::f_sign(tlapack_detail::tsqrt(scale) * tlapack_detail::tsqrt(z), p);
			a = d + z;
			d = d - (bcmax / z) * bcmis;
			const T tau = tlapy2(c, z);
			cs = z / tau;
			sn = c / tau;
			b = b - c;
			c = T(0);
		}
		else {
			int count = 0;
			T sigma = b + c;
			for (;;) {
				count++;
				scale = std::max(tlapack_detail::tabs(temp), tlapack_detail::tabs(sigma));
				if (scale >= safmx2) {
					sigma = sigma * safmn2;
					temp = temp * safmn2;
					if (count <= 20) {
						continue;
					}
				}
				if (scale <= safmn2) {
					sigma = sigma * safmx2;
					temp = temp * safmx2;
					if (count <= 20) {
						continue;
					}
				}
				break;
			}
			p = (T(1) / T(2)) * temp;
			const T tau = tlapy2(sigma, temp);
			cs = tlapack_detail::tsqrt((T(1) / T(2)) * (T(1) + tlapack_detail::tabs(sigma) / tau));
			sn = -(p / (tau * cs)) * detail::f_sign(T(1), sigma);
			const T aa = a * cs + b * sn;
			const T bb = -a * sn + b * cs;
			const T cc = c * cs + d * sn;
			const T dd = -c * sn + d * cs;
			a = aa * cs + cc * sn;
			b = (bb * cs) + (dd * sn);
			c = -(aa * sn) + (cc * cs);
			d = -bb * sn + dd * cs;
			temp = (T(1) / T(2)) * (a + d);
			a = temp;
			d = temp;
			if (c != T(0)) {
				if (b != T(0)) {
					if (detail::f_sign(T(1), b) == detail::f_sign(T(1), c)) {
						const T sab = tlapack_detail::tsqrt(tlapack_detail::tabs(b));
						const T sac = tlapack_detail::tsqrt(tlapack_detail::tabs(c));
						p = detail::f_sign(sab * sac, c);
						const T tau2 = T(1) / tlapack_detail::tsqrt(tlapack_detail::tabs(b + c));
						a = temp + p;
						d = temp - p;
						b = b - c;
						c = T(0);
						const T cs1 = sab * tau2;
						const T sn1 = sac * tau2;
						temp = cs * cs1 - sn * sn1;
						sn = cs * sn1 + sn * cs1;
						cs = temp;
					}
				}
				else {
					b = -c;
					c = T(0);
					temp = cs;
					cs = -sn;
					sn = temp;
				}
			}
		}
	}
	rt1r = a;
	rt2r = d;
	if (c == T(0)) {
		rt1i = T(0);
		rt2i = T(0);
	}
	else {
		rt1i = tlapack_detail::tsqrt(tlapack_detail::tabs(b)) * tlapack_detail::tsqrt(tlapack_detail::tabs(c));
		rt2i = -rt1i;
	}
}

// A := A * (cto/cfrom) (overflow/underflow を避けて段階的に乗算)
// type: 'G' 全体, 'L' 下三角, 'U' 上三角, 'H' Hessenberg,
//       'B' 対称帯下三角格納, 'Q' 対称帯上三角格納, 'Z' 帯 (LU 用)
template <typename T>
inline void tlascl(const char type, const int kl, const int ku, const T& cfrom, const T& cto,
	const int m, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	int itype;
	if (detail::option_is(type, 'G')) {
		itype = 0;
	}
	else if (detail::option_is(type, 'L')) {
		itype = 1;
	}
	else if (detail::option_is(type, 'U')) {
		itype = 2;
	}
	else if (detail::option_is(type, 'H')) {
		itype = 3;
	}
	else if (detail::option_is(type, 'B')) {
		itype = 4;
	}
	else if (detail::option_is(type, 'Q')) {
		itype = 5;
	}
	else if (detail::option_is(type, 'Z')) {
		itype = 6;
	}
	else {
		detail::tlapack_error("tlascl: invalid type");
		return;
	}
	if (cfrom == T(0) || tlapack_detail::tisnan(cfrom) || tlapack_detail::tisnan(cto)) {
		detail::tlapack_error("tlascl: invalid cfrom/cto");
	}
	if (m < 0 || n < 0) {
		detail::tlapack_error("tlascl: invalid m/n");
	}
	if (n == 0 || m == 0) {
		return;
	}
	const T smlnum = tlamch<T>('S');
	const T bignum = T(1) / smlnum;
	T cfromc = cfrom;
	T ctoc = cto;
	bool done = false;
	while (!done) {
		T mul;
		const T cfrom1 = cfromc * smlnum;
		if (cfrom1 == cfromc) {
			// cfromc は inf: mul は正規の値
			mul = ctoc / cfromc;
			done = true;
		}
		else {
			const T cto1 = ctoc / bignum;
			if (cto1 == ctoc) {
				// ctoc は 0 または inf
				mul = ctoc;
				done = true;
				cfromc = T(1);
			}
			else if (tlapack_detail::tabs(cfrom1) > tlapack_detail::tabs(ctoc) && ctoc != T(0)) {
				mul = smlnum;
				done = false;
				cfromc = cfrom1;
			}
			else if (tlapack_detail::tabs(cto1) > tlapack_detail::tabs(cfromc)) {
				mul = bignum;
				done = false;
				ctoc = cto1;
			}
			else {
				mul = ctoc / cfromc;
				done = true;
				if (mul == T(1)) {
					return;
				}
			}
		}
		if (itype == 0) {
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i < m; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 1) {
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = j; i < m; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 2) {
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i <= std::min(j, m - 1); i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 3) {
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i <= std::min(j + 1, m - 1); i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 4) {
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i <= std::min(kl, n - 1 - j); i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 5) {
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = std::max(ku - j, 0); i <= ku; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else { // itype == 6
			for (int j = 0; j < n; j++) {
				T* a = A + static_cast<std::size_t>(lda) * j;
				const int lo = std::max(kl + ku - j, kl);
				const int hi = std::min(2 * kl + ku, kl + ku + m - 1 - j);
				for (int i = lo; i <= hi; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
	}
}

// 平面回転列の適用: side='L' で A := P*A, side='R' で A := A*P^T
// pivot: 'V' 隣接 (i,i+1), 'T' 先頭 (1,i), 'B' 末尾 (i,last)
// direct: 'F' 前進, 'B' 後退
template <typename T>
inline void tlasr(const char side, const char pivot, const char direct, const int m, const int n,
	const T* c, const T* s, T* A, const int lda) {
	namespace detail = tlapack_detail;
	if (!(detail::option_is(side, 'L') || detail::option_is(side, 'R'))) {
		detail::tlapack_error("tlasr: invalid side");
	}
	if (!(detail::option_is(pivot, 'V') || detail::option_is(pivot, 'T') || detail::option_is(pivot, 'B'))) {
		detail::tlapack_error("tlasr: invalid pivot");
	}
	if (!(detail::option_is(direct, 'F') || detail::option_is(direct, 'B'))) {
		detail::tlapack_error("tlasr: invalid direct");
	}
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tlasr: invalid m/n/lda");
	}
	if (m == 0 || n == 0) {
		return;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (detail::option_is(side, 'L')) {
		if (detail::option_is(pivot, 'V')) {
			const int j0 = detail::option_is(direct, 'F') ? 0 : m - 2;
			const int j1 = detail::option_is(direct, 'F') ? m - 1 : -1;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const T ctemp = c[j];
				const T stemp = s[j];
				if (ctemp != T(1) || stemp != T(0)) {
					for (int i = 0; i < n; i++) {
						const T temp = A[(j + 1) + L * i];
						A[(j + 1) + L * i] = ctemp * temp - stemp * A[j + L * i];
						A[j + L * i] = stemp * temp + ctemp * A[j + L * i];
					}
				}
			}
		}
		else if (detail::option_is(pivot, 'T')) {
			const int j0 = detail::option_is(direct, 'F') ? 1 : m - 1;
			const int j1 = detail::option_is(direct, 'F') ? m : 0;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const T ctemp = c[j - 1];
				const T stemp = s[j - 1];
				if (ctemp != T(1) || stemp != T(0)) {
					for (int i = 0; i < n; i++) {
						const T temp = A[j + L * i];
						A[j + L * i] = ctemp * temp - stemp * A[0 + L * i];
						A[0 + L * i] = stemp * temp + ctemp * A[0 + L * i];
					}
				}
			}
		}
		else { // pivot 'B'
			const int j0 = detail::option_is(direct, 'F') ? 0 : m - 2;
			const int j1 = detail::option_is(direct, 'F') ? m - 1 : -1;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const T ctemp = c[j];
				const T stemp = s[j];
				if (ctemp != T(1) || stemp != T(0)) {
					for (int i = 0; i < n; i++) {
						const T temp = A[j + L * i];
						A[j + L * i] = stemp * A[(m - 1) + L * i] + ctemp * temp;
						A[(m - 1) + L * i] = ctemp * A[(m - 1) + L * i] - stemp * temp;
					}
				}
			}
		}
	}
	else {
		if (detail::option_is(pivot, 'V')) {
			const int j0 = detail::option_is(direct, 'F') ? 0 : n - 2;
			const int j1 = detail::option_is(direct, 'F') ? n - 1 : -1;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const T ctemp = c[j];
				const T stemp = s[j];
				if (ctemp != T(1) || stemp != T(0)) {
					for (int i = 0; i < m; i++) {
						const T temp = A[i + L * (j + 1)];
						A[i + L * (j + 1)] = ctemp * temp - stemp * A[i + L * j];
						A[i + L * j] = stemp * temp + ctemp * A[i + L * j];
					}
				}
			}
		}
		else if (detail::option_is(pivot, 'T')) {
			const int j0 = detail::option_is(direct, 'F') ? 1 : n - 1;
			const int j1 = detail::option_is(direct, 'F') ? n : 0;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const T ctemp = c[j - 1];
				const T stemp = s[j - 1];
				if (ctemp != T(1) || stemp != T(0)) {
					for (int i = 0; i < m; i++) {
						const T temp = A[i + L * j];
						A[i + L * j] = ctemp * temp - stemp * A[i + L * 0];
						A[i + L * 0] = stemp * temp + ctemp * A[i + L * 0];
					}
				}
			}
		}
		else { // pivot 'B'
			const int j0 = detail::option_is(direct, 'F') ? 0 : n - 2;
			const int j1 = detail::option_is(direct, 'F') ? n - 1 : -1;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const T ctemp = c[j];
				const T stemp = s[j];
				if (ctemp != T(1) || stemp != T(0)) {
					for (int i = 0; i < m; i++) {
						const T temp = A[i + L * j];
						A[i + L * j] = stemp * A[i + L * (n - 1)] + ctemp * temp;
						A[i + L * (n - 1)] = ctemp * A[i + L * (n - 1)] - stemp * temp;
					}
				}
			}
		}
	}
}

} // namespace vcp

namespace tlapack_detail {

// LAPACK の norm 系で使う max 更新 (NaN を伝播させる)
template <typename T>
inline void norm_max(T& value, const T& sum) {
	if (value < sum || tlapack_detail::tisnan(sum)) {
		value = sum;
	}
}

} // namespace tlapack_detail

namespace vcp {

// 一般行列の norm ('M' max|a_ij|, '1'/'O' 1-norm, 'I' inf-norm, 'F'/'E' Frobenius)
template <typename T>
inline T tlange(const char norm, const int m, const int n, const T* A, const int lda) {
	namespace detail = tlapack_detail;
	if (std::min(m, n) == 0) {
		return T(0);
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	T value = T(0);
	if (detail::option_is(norm, 'M')) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				detail::norm_max(value, tlapack_detail::tabs(A[i + L * j]));
			}
		}
	}
	else if (detail::option_is(norm, 'O') || norm == '1') {
		for (int j = 0; j < n; j++) {
			T sum = T(0);
			for (int i = 0; i < m; i++) {
				sum += tlapack_detail::tabs(A[i + L * j]);
			}
			detail::norm_max(value, sum);
		}
	}
	else if (detail::option_is(norm, 'I')) {
		std::vector<T> work(static_cast<std::size_t>(m), T(0));
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				work[i] += tlapack_detail::tabs(A[i + L * j]);
			}
		}
		for (int i = 0; i < m; i++) {
			detail::norm_max(value, work[i]);
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		T scale = T(0);
		T sum = T(1);
		for (int j = 0; j < n; j++) {
			tlassq(m, A + L * j, 1, scale, sum);
		}
		value = scale * tlapack_detail::tsqrt(sum);
	}
	else {
		detail::tlapack_error("tlange: invalid norm");
	}
	return value;
}

// 対称行列の norm
template <typename T>
inline T tlansy(const char norm, const char uplo, const int n, const T* A, const int lda) {
	namespace detail = tlapack_detail;
	if (n == 0) {
		return T(0);
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const bool upper = detail::option_is(uplo, 'U');
	T value = T(0);
	if (detail::option_is(norm, 'M')) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= j; i++) {
					detail::norm_max(value, tlapack_detail::tabs(A[i + L * j]));
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				for (int i = j; i < n; i++) {
					detail::norm_max(value, tlapack_detail::tabs(A[i + L * j]));
				}
			}
		}
	}
	else if (detail::option_is(norm, 'I') || detail::option_is(norm, 'O') || norm == '1') {
		std::vector<T> work(static_cast<std::size_t>(n), T(0));
		if (upper) {
			for (int j = 0; j < n; j++) {
				T sum = T(0);
				for (int i = 0; i < j; i++) {
					const T absa = tlapack_detail::tabs(A[i + L * j]);
					sum += absa;
					work[i] += absa;
				}
				work[j] = sum + tlapack_detail::tabs(A[j + L * j]);
			}
			for (int i = 0; i < n; i++) {
				detail::norm_max(value, work[i]);
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				T sum = work[j] + tlapack_detail::tabs(A[j + L * j]);
				for (int i = j + 1; i < n; i++) {
					const T absa = tlapack_detail::tabs(A[i + L * j]);
					sum += absa;
					work[i] += absa;
				}
				detail::norm_max(value, sum);
			}
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		T scale = T(0);
		T sum = T(1);
		if (upper) {
			for (int j = 1; j < n; j++) {
				tlassq(j, A + L * j, 1, scale, sum);
			}
		}
		else {
			for (int j = 0; j < n - 1; j++) {
				tlassq(n - j - 1, A + (j + 1) + L * j, 1, scale, sum);
			}
		}
		sum = T(2) * sum;
		tlassq(n, A, lda + 1, scale, sum);
		value = scale * tlapack_detail::tsqrt(sum);
	}
	else {
		detail::tlapack_error("tlansy: invalid norm");
	}
	return value;
}

// 対称三重対角行列の norm (d: 対角 n 個, e: 副対角 n-1 個)
template <typename T>
inline T tlanst(const char norm, const int n, const T* d, const T* e) {
	namespace detail = tlapack_detail;
	if (n <= 0) {
		return T(0);
	}
	T value = T(0);
	if (detail::option_is(norm, 'M')) {
		value = tlapack_detail::tabs(d[n - 1]);
		for (int i = 0; i < n - 1; i++) {
			detail::norm_max(value, tlapack_detail::tabs(d[i]));
			detail::norm_max(value, tlapack_detail::tabs(e[i]));
		}
	}
	else if (detail::option_is(norm, 'O') || norm == '1' || detail::option_is(norm, 'I')) {
		if (n == 1) {
			value = tlapack_detail::tabs(d[0]);
		}
		else {
			value = tlapack_detail::tabs(d[0]) + tlapack_detail::tabs(e[0]);
			detail::norm_max(value, tlapack_detail::tabs(e[n - 2]) + tlapack_detail::tabs(d[n - 1]));
			for (int i = 1; i < n - 1; i++) {
				detail::norm_max(value, tlapack_detail::tabs(d[i]) + tlapack_detail::tabs(e[i]) + tlapack_detail::tabs(e[i - 1]));
			}
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		T scale = T(0);
		T sum = T(1);
		if (n > 1) {
			tlassq(n - 1, e, 1, scale, sum);
			sum = T(2) * sum;
		}
		tlassq(n, d, 1, scale, sum);
		value = scale * tlapack_detail::tsqrt(sum);
	}
	else {
		detail::tlapack_error("tlanst: invalid norm");
	}
	return value;
}

// 一般帯行列 (n x n, 帯幅 kl/ku, 帯格納) の norm
template <typename T>
inline T tlangb(const char norm, const int n, const int kl, const int ku, const T* AB, const int ldab) {
	namespace detail = tlapack_detail;
	if (n == 0) {
		return T(0);
	}
	const std::size_t L = static_cast<std::size_t>(ldab);
	T value = T(0);
	if (detail::option_is(norm, 'M')) {
		for (int j = 0; j < n; j++) {
			for (int i = std::max(ku - j, 0); i <= std::min(n + ku - 1 - j, kl + ku); i++) {
				detail::norm_max(value, tlapack_detail::tabs(AB[i + L * j]));
			}
		}
	}
	else if (detail::option_is(norm, 'O') || norm == '1') {
		for (int j = 0; j < n; j++) {
			T sum = T(0);
			for (int i = std::max(ku - j, 0); i <= std::min(n + ku - 1 - j, kl + ku); i++) {
				sum += tlapack_detail::tabs(AB[i + L * j]);
			}
			detail::norm_max(value, sum);
		}
	}
	else if (detail::option_is(norm, 'I')) {
		std::vector<T> work(static_cast<std::size_t>(n), T(0));
		for (int j = 0; j < n; j++) {
			for (int i = std::max(0, j - ku); i <= std::min(n - 1, j + kl); i++) {
				work[i] += tlapack_detail::tabs(AB[(ku + i - j) + L * j]);
			}
		}
		for (int i = 0; i < n; i++) {
			detail::norm_max(value, work[i]);
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		T scale = T(0);
		T sum = T(1);
		for (int j = 0; j < n; j++) {
			const int l = std::max(0, j - ku);
			const int k = ku - j + l;
			tlassq(std::min(n - 1, j + kl) - l + 1, AB + k + L * j, 1, scale, sum);
		}
		value = scale * tlapack_detail::tsqrt(sum);
	}
	else {
		detail::tlapack_error("tlangb: invalid norm");
	}
	return value;
}

// 対称帯行列 (帯幅 k, 帯格納) の norm
template <typename T>
inline T tlansb(const char norm, const char uplo, const int n, const int k, const T* AB, const int ldab) {
	namespace detail = tlapack_detail;
	if (n == 0) {
		return T(0);
	}
	const std::size_t L = static_cast<std::size_t>(ldab);
	const bool upper = detail::option_is(uplo, 'U');
	T value = T(0);
	if (detail::option_is(norm, 'M')) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				for (int i = std::max(k - j, 0); i <= k; i++) {
					detail::norm_max(value, tlapack_detail::tabs(AB[i + L * j]));
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= std::min(n - 1 - j, k); i++) {
					detail::norm_max(value, tlapack_detail::tabs(AB[i + L * j]));
				}
			}
		}
	}
	else if (detail::option_is(norm, 'I') || detail::option_is(norm, 'O') || norm == '1') {
		std::vector<T> work(static_cast<std::size_t>(n), T(0));
		if (upper) {
			for (int j = 0; j < n; j++) {
				T sum = T(0);
				for (int i = std::max(0, j - k); i < j; i++) {
					const T absa = tlapack_detail::tabs(AB[(k + i - j) + L * j]);
					sum += absa;
					work[i] += absa;
				}
				work[j] = sum + tlapack_detail::tabs(AB[k + L * j]);
			}
			for (int i = 0; i < n; i++) {
				detail::norm_max(value, work[i]);
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				T sum = work[j] + tlapack_detail::tabs(AB[0 + L * j]);
				for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
					const T absa = tlapack_detail::tabs(AB[(i - j) + L * j]);
					sum += absa;
					work[i] += absa;
				}
				detail::norm_max(value, sum);
			}
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		T scale = T(0);
		T sum = T(1);
		if (k > 0) {
			if (upper) {
				for (int j = 1; j < n; j++) {
					tlassq(std::min(j, k), AB + std::max(k - j, 0) + L * j, 1, scale, sum);
				}
			}
			else {
				for (int j = 0; j < n - 1; j++) {
					tlassq(std::min(n - 1 - j, k), AB + 1 + L * j, 1, scale, sum);
				}
			}
			sum = T(2) * sum;
		}
		tlassq(n, AB + (upper ? k : 0), ldab, scale, sum);
		value = scale * tlapack_detail::tsqrt(sum);
	}
	else {
		detail::tlapack_error("tlansb: invalid norm");
	}
	return value;
}

// Householder reflector の生成: H * [alpha; x] = [beta; 0]
// alpha は beta に上書き，x は v(2:n) に上書きされる
template <typename T>
inline void tlarfg(const int n, T& alpha, T* x, const int incx, T& tau) {
	namespace detail = tlapack_detail;
	if (n <= 1) {
		tau = T(0);
		return;
	}
	T xnorm = tnrm2(n - 1, x, incx);
	if (xnorm == T(0)) {
		tau = T(0);
	}
	else {
		T beta = -detail::f_sign(tlapy2(alpha, xnorm), alpha);
		const T safmin = tlamch<T>('S') / (tlamch<T>('E'));
		int knt = 0;
		if (tlapack_detail::tabs(beta) < safmin) {
			const T rsafmn = T(1) / safmin;
			do {
				knt++;
				tscal(n - 1, rsafmn, x, incx);
				beta = beta * rsafmn;
				alpha = alpha * rsafmn;
			} while (tlapack_detail::tabs(beta) < safmin && knt < 20);
			xnorm = tnrm2(n - 1, x, incx);
			beta = -detail::f_sign(tlapy2(alpha, xnorm), alpha);
		}
		tau = (beta - alpha) / beta;
		tscal(n - 1, T(1) / (alpha - beta), x, incx);
		for (int j = 0; j < knt; j++) {
			beta = beta * safmin;
		}
		alpha = beta;
	}
}

// Householder reflector の適用: C := H*C (side='L') または C*H (side='R'),
// H = I - tau*v*v^T
template <typename T>
inline void tlarf(const char side, const int m, const int n, const T* v, const int incv,
	const T& tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool applyleft = detail::option_is(side, 'L');
	int lastv = 0;
	int lastc = 0;
	if (tau != T(0)) {
		lastv = applyleft ? m : n;
		int i = (incv > 0) ? (lastv - 1) * incv : 0;
		while (lastv > 0 && v[i] == T(0)) {
			lastv--;
			i -= incv;
		}
		if (applyleft) {
			lastc = ilatlc(lastv, n, C, ldc);
		}
		else {
			lastc = ilatlr(m, lastv, C, ldc);
		}
	}
	if (lastv <= 0) {
		return;
	}
	std::vector<T> work(static_cast<std::size_t>(lastc) + 1, T(0));
	if (applyleft) {
		tgemv('T', lastv, lastc, T(1), C, ldc, v, incv, T(0), work.data(), 1);
		tger(lastv, lastc, -tau, v, incv, work.data(), 1, C, ldc);
	}
	else {
		tgemv('N', lastc, lastv, T(1), C, ldc, v, incv, T(0), work.data(), 1);
		tger(lastc, lastv, -tau, work.data(), 1, v, incv, C, ldc);
	}
}

// tlarf の v(1) = 1 を暗黙とする版 (v の先頭要素は参照されない)
template <typename T>
inline void tlarf1f(const char side, const int m, const int n, const T* v, const int incv,
	const T& tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool applyleft = detail::option_is(side, 'L');
	int lastv = 1;
	int lastc = 0;
	if (tau != T(0)) {
		lastv = applyleft ? m : n;
		int i = (incv > 0) ? (lastv - 1) * incv : 0;
		while (lastv > 1 && v[i] == T(0)) {
			lastv--;
			i -= incv;
		}
		if (applyleft) {
			lastc = ilatlc(lastv, n, C, ldc);
		}
		else {
			lastc = ilatlr(m, lastv, C, ldc);
		}
	}
	if (lastc == 0) {
		return;
	}
	if (applyleft) {
		if (lastv == 1) {
			tscal(lastc, T(1) - tau, C, ldc);
		}
		else {
			std::vector<T> work(static_cast<std::size_t>(lastc), T(0));
			tgemv('T', lastv - 1, lastc, T(1), C + 1, ldc, v + incv, incv, T(0), work.data(), 1);
			taxpy(lastc, T(1), C, ldc, work.data(), 1);
			taxpy(lastc, -tau, work.data(), 1, C, ldc);
			tger(lastv - 1, lastc, -tau, v + incv, incv, work.data(), 1, C + 1, ldc);
		}
	}
	else {
		if (lastv == 1) {
			tscal(lastc, T(1) - tau, C, 1);
		}
		else {
			std::vector<T> work(static_cast<std::size_t>(lastc), T(0));
			tgemv('N', lastc, lastv - 1, T(1), C + static_cast<std::size_t>(ldc), ldc, v + incv, incv, T(0), work.data(), 1);
			taxpy(lastc, T(1), C, 1, work.data(), 1);
			taxpy(lastc, -tau, work.data(), 1, C, 1);
			tger(lastc, lastv - 1, -tau, work.data(), 1, v + incv, incv, C + static_cast<std::size_t>(ldc), ldc);
		}
	}
}

// tlarf の v(last) = 1 を暗黙とする版 (v の末尾要素は参照されない)
template <typename T>
inline void tlarf1l(const char side, const int m, const int n, const T* v, const int incv,
	const T& tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool applyleft = detail::option_is(side, 'L');
	int firstv = 1;
	int lastv = 0;
	int lastc = 0;
	int i = 0;
	if (tau != T(0)) {
		lastv = applyleft ? m : n;
		while (lastv > firstv && v[i] == T(0)) {
			firstv++;
			i += incv;
		}
		if (applyleft) {
			lastc = ilatlc(lastv, n, C, ldc);
		}
		else {
			lastc = ilatlr(m, lastv, C, ldc);
		}
	}
	if (lastc == 0 || lastv <= 0) {
		return;
	}
	const std::size_t L = static_cast<std::size_t>(ldc);
	if (applyleft) {
		if (lastv == firstv) {
			tscal(lastc, T(1) - tau, C + (firstv - 1), ldc);
		}
		else {
			std::vector<T> work(static_cast<std::size_t>(lastc), T(0));
			tgemv('T', lastv - firstv, lastc, T(1), C + (firstv - 1), ldc, v + i, incv, T(0), work.data(), 1);
			taxpy(lastc, T(1), C + (lastv - 1), ldc, work.data(), 1);
			taxpy(lastc, -tau, work.data(), 1, C + (lastv - 1), ldc);
			tger(lastv - firstv, lastc, -tau, v + i, incv, work.data(), 1, C + (firstv - 1), ldc);
		}
	}
	else {
		if (lastv == firstv) {
			// H = I - tau*e*e^T (e は第 lastv 単位 vector): 列 lastv のみ scale
			// (reference 3.12.1 は C(1,1) を scale するが，正しくは列 lastv)
			tscal(lastc, T(1) - tau, C + L * (lastv - 1), 1);
		}
		else {
			std::vector<T> work(static_cast<std::size_t>(lastc), T(0));
			tgemv('N', lastc, lastv - firstv, T(1), C + L * (firstv - 1), ldc, v + i, incv, T(0), work.data(), 1);
			taxpy(lastc, T(1), C + L * (lastv - 1), 1, work.data(), 1);
			taxpy(lastc, -tau, work.data(), 1, C + L * (lastv - 1), 1);
			tger(lastc, lastv - firstv, -tau, work.data(), 1, v + i, incv, C + L * (firstv - 1), ldc);
		}
	}
}

// block reflector の三角因子 Tf の生成 (recursive, reference LAPACK 3.12.1)
// direct: 'F'/'B', storev: 'C'/'R'
template <typename T>
inline void tlarft(const char direct, const char storev, const int n, const int k,
	const T* V, const int ldv, const T* tau, T* Tf, const int ldt) {
	namespace detail = tlapack_detail;
	if (n == 0 || k == 0) {
		return;
	}
	if (n == 1 || k == 1) {
		Tf[0] = tau[0];
		return;
	}
	const std::size_t LV = static_cast<std::size_t>(ldv);
	const std::size_t LT = static_cast<std::size_t>(ldt);
	const int l = k / 2;
	const bool dirf = detail::option_is(direct, 'F');
	const bool colv = detail::option_is(storev, 'C');
	if (dirf && colv) { // QR
		tlarft(direct, storev, n, l, V, ldv, tau, Tf, ldt);
		tlarft(direct, storev, n - l, k - l, V + l + LV * l, ldv, tau + l, Tf + l + LT * l, ldt);
		for (int j = 0; j < l; j++) {
			for (int i = 0; i < k - l; i++) {
				Tf[j + LT * (l + i)] = V[(l + i) + LV * j];
			}
		}
		ttrmm('R', 'L', 'N', 'U', l, k - l, T(1), V + l + LV * l, ldv, Tf + LT * l, ldt);
		if (n > k) {
			tgemm('T', 'N', l, k - l, n - k, T(1), V + k, ldv, V + k + LV * l, ldv, T(1), Tf + LT * l, ldt);
		}
		ttrmm('L', 'U', 'N', 'N', l, k - l, -T(1), Tf, ldt, Tf + LT * l, ldt);
		ttrmm('R', 'U', 'N', 'N', l, k - l, T(1), Tf + l + LT * l, ldt, Tf + LT * l, ldt);
	}
	else if (dirf && !colv) { // LQ
		tlarft(direct, storev, n, l, V, ldv, tau, Tf, ldt);
		tlarft(direct, storev, n - l, k - l, V + l + LV * l, ldv, tau + l, Tf + l + LT * l, ldt);
		tlacpy('A', l, k - l, V + LV * l, ldv, Tf + LT * l, ldt);
		ttrmm('R', 'U', 'T', 'U', l, k - l, T(1), V + l + LV * l, ldv, Tf + LT * l, ldt);
		if (n > k) {
			tgemm('N', 'T', l, k - l, n - k, T(1), V + LV * k, ldv, V + l + LV * k, ldv, T(1), Tf + LT * l, ldt);
		}
		ttrmm('L', 'U', 'N', 'N', l, k - l, -T(1), Tf, ldt, Tf + LT * l, ldt);
		ttrmm('R', 'U', 'N', 'N', l, k - l, T(1), Tf + l + LT * l, ldt, Tf + LT * l, ldt);
	}
	else if (!dirf && colv) { // QL
		tlarft(direct, storev, n - l, k - l, V, ldv, tau, Tf, ldt);
		tlarft(direct, storev, n, l, V + LV * (k - l), ldv, tau + (k - l), Tf + (k - l) + LT * (k - l), ldt);
		for (int j = 0; j < k - l; j++) {
			for (int i = 0; i < l; i++) {
				Tf[(k - l + i) + LT * j] = V[(n - k + j) + LV * (k - l + i)];
			}
		}
		ttrmm('R', 'U', 'N', 'U', l, k - l, T(1), V + (n - k), ldv, Tf + (k - l), ldt);
		if (n > k) {
			tgemm('T', 'N', l, k - l, n - k, T(1), V + LV * (k - l), ldv, V, ldv, T(1), Tf + (k - l), ldt);
		}
		ttrmm('L', 'L', 'N', 'N', l, k - l, -T(1), Tf + (k - l) + LT * (k - l), ldt, Tf + (k - l), ldt);
		ttrmm('R', 'L', 'N', 'N', l, k - l, T(1), Tf, ldt, Tf + (k - l), ldt);
	}
	else { // RQ (direct='B', storev='R')
		tlarft(direct, storev, n - l, k - l, V, ldv, tau, Tf, ldt);
		tlarft(direct, storev, n, l, V + (k - l), ldv, tau + (k - l), Tf + (k - l) + LT * (k - l), ldt);
		tlacpy('A', l, k - l, V + (k - l) + LV * (n - k), ldv, Tf + (k - l), ldt);
		ttrmm('R', 'L', 'T', 'U', l, k - l, T(1), V + LV * (n - k), ldv, Tf + (k - l), ldt);
		if (n > k) {
			tgemm('N', 'T', l, k - l, n - k, T(1), V + (k - l), ldv, V, ldv, T(1), Tf + (k - l), ldt);
		}
		ttrmm('L', 'L', 'N', 'N', l, k - l, -T(1), Tf + (k - l) + LT * (k - l), ldt, Tf + (k - l), ldt);
		ttrmm('R', 'L', 'N', 'N', l, k - l, T(1), Tf, ldt, Tf + (k - l), ldt);
	}
}

// block reflector の適用: C := H*C, H^T*C, C*H, C*H^T
template <typename T>
inline void tlarfb(const char side, const char trans, const char direct, const char storev,
	const int m, const int n, const int k, const T* V, const int ldv,
	const T* Tf, const int ldt, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	if (m <= 0 || n <= 0) {
		return;
	}
	const char transt = detail::option_is(trans, 'N') ? 'T' : 'N';
	const std::size_t LV = static_cast<std::size_t>(ldv);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const int ldwork = std::max(1, detail::option_is(side, 'L') ? n : m);
	const std::size_t LW = static_cast<std::size_t>(ldwork);
	std::vector<T> work(LW * static_cast<std::size_t>(k), T(0));
	T* W = work.data();
	if (detail::option_is(storev, 'C')) {
		if (detail::option_is(direct, 'F')) {
			if (detail::option_is(side, 'L')) {
				// W := C1^T * V1, C2 -= V2*W^T, ...
				for (int j = 0; j < k; j++) {
					tcopy(n, C + j, ldc, W + LW * j, 1);
				}
				ttrmm('R', 'L', 'N', 'U', n, k, T(1), V, ldv, W, ldwork);
				if (m > k) {
					tgemm('T', 'N', n, k, m - k, T(1), C + k, ldc, V + k, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'U', transt, 'N', n, k, T(1), Tf, ldt, W, ldwork);
				if (m > k) {
					tgemm('N', 'T', m - k, n, k, -T(1), V + k, ldv, W, ldwork, T(1), C + k, ldc);
				}
				ttrmm('R', 'L', 'T', 'U', n, k, T(1), V, ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[j + LC * i] = C[j + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					tcopy(m, C + LC * j, 1, W + LW * j, 1);
				}
				ttrmm('R', 'L', 'N', 'U', m, k, T(1), V, ldv, W, ldwork);
				if (n > k) {
					tgemm('N', 'N', m, k, n - k, T(1), C + LC * k, ldc, V + k, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'U', trans, 'N', m, k, T(1), Tf, ldt, W, ldwork);
				if (n > k) {
					tgemm('N', 'T', m, n - k, k, -T(1), W, ldwork, V + k, ldv, T(1), C + LC * k, ldc);
				}
				ttrmm('R', 'L', 'T', 'U', m, k, T(1), V, ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < m; i++) {
						C[i + LC * j] = C[i + LC * j] - W[i + LW * j];
					}
				}
			}
		}
		else {
			if (detail::option_is(side, 'L')) {
				for (int j = 0; j < k; j++) {
					tcopy(n, C + (m - k + j), ldc, W + LW * j, 1);
				}
				ttrmm('R', 'U', 'N', 'U', n, k, T(1), V + (m - k), ldv, W, ldwork);
				if (m > k) {
					tgemm('T', 'N', n, k, m - k, T(1), C, ldc, V, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'L', transt, 'N', n, k, T(1), Tf, ldt, W, ldwork);
				if (m > k) {
					tgemm('N', 'T', m - k, n, k, -T(1), V, ldv, W, ldwork, T(1), C, ldc);
				}
				ttrmm('R', 'U', 'T', 'U', n, k, T(1), V + (m - k), ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[(m - k + j) + LC * i] = C[(m - k + j) + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					tcopy(m, C + LC * (n - k + j), 1, W + LW * j, 1);
				}
				ttrmm('R', 'U', 'N', 'U', m, k, T(1), V + (n - k), ldv, W, ldwork);
				if (n > k) {
					tgemm('N', 'N', m, k, n - k, T(1), C, ldc, V, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'L', trans, 'N', m, k, T(1), Tf, ldt, W, ldwork);
				if (n > k) {
					tgemm('N', 'T', m, n - k, k, -T(1), W, ldwork, V, ldv, T(1), C, ldc);
				}
				ttrmm('R', 'U', 'T', 'U', m, k, T(1), V + (n - k), ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < m; i++) {
						C[i + LC * (n - k + j)] = C[i + LC * (n - k + j)] - W[i + LW * j];
					}
				}
			}
		}
	}
	else {
		if (detail::option_is(direct, 'F')) {
			if (detail::option_is(side, 'L')) {
				for (int j = 0; j < k; j++) {
					tcopy(n, C + j, ldc, W + LW * j, 1);
				}
				ttrmm('R', 'U', 'T', 'U', n, k, T(1), V, ldv, W, ldwork);
				if (m > k) {
					tgemm('T', 'T', n, k, m - k, T(1), C + k, ldc, V + LV * k, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'U', transt, 'N', n, k, T(1), Tf, ldt, W, ldwork);
				if (m > k) {
					tgemm('T', 'T', m - k, n, k, -T(1), V + LV * k, ldv, W, ldwork, T(1), C + k, ldc);
				}
				ttrmm('R', 'U', 'N', 'U', n, k, T(1), V, ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[j + LC * i] = C[j + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					tcopy(m, C + LC * j, 1, W + LW * j, 1);
				}
				ttrmm('R', 'U', 'T', 'U', m, k, T(1), V, ldv, W, ldwork);
				if (n > k) {
					tgemm('N', 'T', m, k, n - k, T(1), C + LC * k, ldc, V + LV * k, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'U', trans, 'N', m, k, T(1), Tf, ldt, W, ldwork);
				if (n > k) {
					tgemm('N', 'N', m, n - k, k, -T(1), W, ldwork, V + LV * k, ldv, T(1), C + LC * k, ldc);
				}
				ttrmm('R', 'U', 'N', 'U', m, k, T(1), V, ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < m; i++) {
						C[i + LC * j] = C[i + LC * j] - W[i + LW * j];
					}
				}
			}
		}
		else {
			if (detail::option_is(side, 'L')) {
				for (int j = 0; j < k; j++) {
					tcopy(n, C + (m - k + j), ldc, W + LW * j, 1);
				}
				ttrmm('R', 'L', 'T', 'U', n, k, T(1), V + LV * (m - k), ldv, W, ldwork);
				if (m > k) {
					tgemm('T', 'T', n, k, m - k, T(1), C, ldc, V, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'L', transt, 'N', n, k, T(1), Tf, ldt, W, ldwork);
				if (m > k) {
					tgemm('T', 'T', m - k, n, k, -T(1), V, ldv, W, ldwork, T(1), C, ldc);
				}
				ttrmm('R', 'L', 'N', 'U', n, k, T(1), V + LV * (m - k), ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[(m - k + j) + LC * i] = C[(m - k + j) + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					tcopy(m, C + LC * (n - k + j), 1, W + LW * j, 1);
				}
				ttrmm('R', 'L', 'T', 'U', m, k, T(1), V + LV * (n - k), ldv, W, ldwork);
				if (n > k) {
					tgemm('N', 'T', m, k, n - k, T(1), C, ldc, V, ldv, T(1), W, ldwork);
				}
				ttrmm('R', 'L', trans, 'N', m, k, T(1), Tf, ldt, W, ldwork);
				if (n > k) {
					tgemm('N', 'N', m, n - k, k, -T(1), W, ldwork, V, ldv, T(1), C, ldc);
				}
				ttrmm('R', 'L', 'N', 'U', m, k, T(1), V + LV * (n - k), ldv, W, ldwork);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < m; i++) {
						C[i + LC * (n - k + j)] = C[i + LC * (n - k + j)] - W[i + LW * j];
					}
				}
			}
		}
	}
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_AUX_HPP
