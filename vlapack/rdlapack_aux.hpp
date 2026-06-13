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

#ifndef VBLAS_RDLAPACK_AUX_HPP
#define VBLAS_RDLAPACK_AUX_HPP

// LAPACK 補助 routine (丸めモード指定付き)．reference LAPACK 3.12.1 からの移植．
// FP 演算を含む routine は rd 接頭辞 + 末尾引数 rounding_mode．
// FP 演算を含まない routine (dlaset, dlacpy, dlaswp, dlasrt, iladlr, iladlc) は
// LAPACK と同名 (rounding_mode なし)．index は 0-based．

#include <algorithm>
#include <functional>

#include "rdlapack_common.hpp"

#pragma STDC FENV_ACCESS ON

// ---- 演算なし ----

// 末尾の zero 行を除いた行数 (A の非 zero 要素を含む最終行 + 1)
inline int iladlr(const int m, const int n, const double* A, const int lda) {
	if (m == 0) {
		return 0;
	}
	if (A[(m - 1)] != 0.0 || A[(m - 1) + static_cast<std::size_t>(lda) * (n - 1)] != 0.0) {
		return m;
	}
	int result = 0;
	for (int j = 0; j < n; j++) {
		int i = m;
		while (i >= 1 && A[(i - 1) + static_cast<std::size_t>(lda) * j] == 0.0) {
			i--;
		}
		result = std::max(result, i);
	}
	return result;
}

// 末尾の zero 列を除いた列数
inline int iladlc(const int m, const int n, const double* A, const int lda) {
	if (n == 0) {
		return 0;
	}
	if (A[static_cast<std::size_t>(lda) * (n - 1)] != 0.0 ||
		A[(m - 1) + static_cast<std::size_t>(lda) * (n - 1)] != 0.0) {
		return n;
	}
	for (int j = n; j >= 1; j--) {
		for (int i = 0; i < m; i++) {
			if (A[i + static_cast<std::size_t>(lda) * (j - 1)] != 0.0) {
				return j;
			}
		}
	}
	return 0;
}

// A の off-diagonal を alpha，diagonal を beta に設定 (uplo: 'U' 上三角のみ /
// 'L' 下三角のみ / その他 全体)
inline void dlaset(const char uplo, const int m, const int n, const double alpha, const double beta, double* A, const int lda) {
	namespace detail = vblas_rdlapack_detail;
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
inline void dlacpy(const char uplo, const int m, const int n, const double* A, const int lda, double* B, const int ldb) {
	namespace detail = vblas_rdlapack_detail;
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
inline void dlaswp(const int n, double* A, const int lda, const int k1, const int k2, const int* ipiv, const int incx) {
	if (n <= 0) {
		return;
	}
	if (incx > 0) {
		for (int i = k1; i < k2; i++) {
			const int ip = ipiv[i];
			if (ip != i) {
				for (int j = 0; j < n; j++) {
					const std::size_t c = static_cast<std::size_t>(lda) * j;
					const double tmp = A[i + c];
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
					const double tmp = A[i + c];
					A[i + c] = A[ip + c];
					A[ip + c] = tmp;
				}
			}
		}
	}
}

// d を昇順 ('I') または降順 ('D') に sort (比較のみで丸め誤差なし)
inline void dlasrt(const char id, const int n, double* d) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0) {
		detail::rdlapack_error("dlasrt: n < 0");
	}
	if (detail::option_is(id, 'I')) {
		std::sort(d, d + n);
	}
	else if (detail::option_is(id, 'D')) {
		std::sort(d, d + n, std::greater<double>());
	}
	else {
		detail::rdlapack_error("dlasrt: invalid id");
	}
}

// ---- 演算あり (rd 接頭辞 + rounding_mode) ----

// scaled sum of squares: scale^2*sumsq += sum x_i^2 (overflow/underflow 回避付き)
inline void rdlassq(const int n, const double* x, const int incx, double& scale, double& sumsq, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (std::isnan(scale) || std::isnan(sumsq)) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	// LA_CONSTANTS (double): tsml = 2^-511, tbig = 2^486, ssml = 2^537, sbig = 2^-538
	const double tsml = std::ldexp(1.0, -511);
	const double tbig = std::ldexp(1.0, 486);
	const double ssml = std::ldexp(1.0, 537);
	const double sbig = std::ldexp(1.0, -538);
	if (sumsq == 0.0) {
		scale = 1.0;
	}
	if (scale == 0.0) {
		scale = 1.0;
		sumsq = 0.0;
	}
	if (n <= 0) {
		return;
	}
	bool notbig = true;
	double asml = 0.0;
	double amed = 0.0;
	double abig = 0.0;
	int ix = 0;
	if (incx < 0) {
		ix = -(n - 1) * incx;
	}
	for (int i = 0; i < n; i++) {
		const double ax = std::fabs(x[ix]);
		if (ax > tbig) {
			abig = abig + (ax * sbig) * (ax * sbig);
			notbig = false;
		}
		else if (ax < tsml) {
			if (notbig) {
				asml = asml + (ax * ssml) * (ax * ssml);
			}
		}
		else {
			amed = amed + ax * ax;
		}
		ix += incx;
	}
	if (sumsq > 0.0) {
		const double ax = scale * std::sqrt(sumsq);
		if (ax > tbig) {
			if (scale > 1.0) {
				scale = scale * sbig;
				abig = abig + scale * (scale * sumsq);
			}
			else {
				abig = abig + scale * (scale * (sbig * (sbig * sumsq)));
			}
		}
		else if (ax < tsml) {
			if (notbig) {
				if (scale < 1.0) {
					scale = scale * ssml;
					asml = asml + scale * (scale * sumsq);
				}
				else {
					asml = asml + scale * (scale * (ssml * (ssml * sumsq)));
				}
			}
		}
		else {
			amed = amed + scale * (scale * sumsq);
		}
	}
	if (abig > 0.0) {
		if (amed > 0.0 || std::isnan(amed)) {
			abig = abig + (amed * sbig) * sbig;
		}
		scale = 1.0 / sbig;
		sumsq = abig;
	}
	else if (asml > 0.0) {
		if (amed > 0.0 || std::isnan(amed)) {
			amed = std::sqrt(amed);
			asml = std::sqrt(asml) / ssml;
			double ymin, ymax;
			if (asml > amed) {
				ymin = amed;
				ymax = asml;
			}
			else {
				ymin = asml;
				ymax = amed;
			}
			scale = 1.0;
			sumsq = ymax * ymax * (1.0 + (ymin / ymax) * (ymin / ymax));
		}
		else {
			scale = 1.0 / ssml;
			sumsq = asml;
		}
	}
	else {
		scale = 1.0;
		sumsq = amed;
	}
}

// sqrt(x^2 + y^2) (overflow 回避付き)
inline double rdlapy2(const double x, const double y, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (std::isnan(x)) {
		return x;
	}
	if (std::isnan(y)) {
		return y;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double xabs = std::fabs(x);
	const double yabs = std::fabs(y);
	const double w = std::max(xabs, yabs);
	const double z = std::min(xabs, yabs);
	if (z == 0.0 || w > DBL_MAX) {
		return w;
	}
	return w * std::sqrt(1.0 + (z / w) * (z / w));
}

// sqrt(x^2 + y^2 + z^2) (overflow 回避付き)
inline double rdlapy3(const double x, const double y, const double z, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double xabs = std::fabs(x);
	const double yabs = std::fabs(y);
	const double zabs = std::fabs(z);
	const double w = std::max(std::max(xabs, yabs), zabs);
	if (w == 0.0 || w > DBL_MAX) {
		return xabs + yabs + zabs;
	}
	return w * std::sqrt((xabs / w) * (xabs / w) + (yabs / w) * (yabs / w) + (zabs / w) * (zabs / w));
}

namespace vblas_rdlapack_detail {

inline double dladiv2(const double a, const double b, const double c, const double d, const double r, const double t) {
	if (r != 0.0) {
		const double br = b * r;
		if (br != 0.0) {
			return (a + br) * t;
		}
		return a * t + (b * t) * r;
	}
	return (a + d * (b / c)) * t;
}

inline void dladiv1(const double a, const double b, const double c, const double d, double& p, double& q) {
	const double r = d / c;
	const double t = 1.0 / (c + d * r);
	p = dladiv2(a, b, c, d, r, t);
	q = dladiv2(b, -a, c, d, r, t);
}

} // namespace vblas_rdlapack_detail

// 複素数除算 p + iq := (a + ib) / (c + id) (overflow 回避付き, Baudin-Smith)
inline void rdladiv(const double a, const double b, const double c, const double d, double& p, double& q, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	double aa = a;
	double bb = b;
	double cc = c;
	double dd = d;
	const double ab = std::max(std::fabs(a), std::fabs(b));
	const double cd = std::max(std::fabs(c), std::fabs(d));
	double s = 1.0;
	const double ov = DBL_MAX;
	const double un = DBL_MIN;
	const double eps = DBL_EPSILON * 0.5;
	const double be = 2.0 / (eps * eps);
	if (ab >= 0.5 * ov) {
		aa = 0.5 * aa;
		bb = 0.5 * bb;
		s = 2.0 * s;
	}
	if (cd >= 0.5 * ov) {
		cc = 0.5 * cc;
		dd = 0.5 * dd;
		s = 0.5 * s;
	}
	if (ab <= un * 2.0 / eps) {
		aa = aa * be;
		bb = bb * be;
		s = s / be;
	}
	if (cd <= un * 2.0 / eps) {
		cc = cc * be;
		dd = dd * be;
		s = s * be;
	}
	if (std::fabs(d) <= std::fabs(c)) {
		detail::dladiv1(aa, bb, cc, dd, p, q);
	}
	else {
		detail::dladiv1(bb, aa, dd, cc, p, q);
		q = -q;
	}
	p = p * s;
	q = q * s;
}

// Givens 回転の生成: [c s; -s c]^T [f; g] = [r; 0]
inline void rdlartg(const double f, const double g, double& cs, double& sn, double& r, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double safmin = DBL_MIN;
	const double safmax = 1.0 / DBL_MIN;
	const double rtmin = std::sqrt(safmin);
	const double rtmax = std::sqrt(safmax / 2.0);
	const double f1 = std::fabs(f);
	const double g1 = std::fabs(g);
	if (g == 0.0) {
		cs = 1.0;
		sn = 0.0;
		r = f;
	}
	else if (f == 0.0) {
		cs = 0.0;
		sn = detail::f_sign(1.0, g);
		r = g1;
	}
	else if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
		const double d = std::sqrt(f * f + g * g);
		cs = f1 / d;
		r = detail::f_sign(d, f);
		sn = g / r;
	}
	else {
		const double u = std::min(safmax, std::max(std::max(safmin, f1), g1));
		const double fs = f / u;
		const double gs = g / u;
		const double d = std::sqrt(fs * fs + gs * gs);
		cs = std::fabs(fs) / d;
		r = detail::f_sign(d, f);
		sn = gs / r;
		r = r * u;
	}
}

// 対称 2x2 行列 [a b; b c] の固有値 (rt1 >= rt2)
inline void rdlae2(const double a, const double b, const double c, double& rt1, double& rt2, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double sm = a + c;
	const double df = a - c;
	const double adf = std::fabs(df);
	const double tb = b + b;
	const double ab = std::fabs(tb);
	double acmx, acmn;
	if (std::fabs(a) > std::fabs(c)) {
		acmx = a;
		acmn = c;
	}
	else {
		acmx = c;
		acmn = a;
	}
	double rt;
	if (adf > ab) {
		rt = adf * std::sqrt(1.0 + (ab / adf) * (ab / adf));
	}
	else if (adf < ab) {
		rt = ab * std::sqrt(1.0 + (adf / ab) * (adf / ab));
	}
	else {
		rt = ab * std::sqrt(2.0);
	}
	if (sm < 0.0) {
		rt1 = 0.5 * (sm - rt);
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else if (sm > 0.0) {
		rt1 = 0.5 * (sm + rt);
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else {
		rt1 = 0.5 * rt;
		rt2 = -0.5 * rt;
	}
}

// 対称 2x2 行列 [a b; b c] の固有値と固有 vector (回転 cs1, sn1)
inline void rdlaev2(const double a, const double b, const double c, double& rt1, double& rt2, double& cs1, double& sn1, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double sm = a + c;
	const double df = a - c;
	const double adf = std::fabs(df);
	const double tb = b + b;
	const double ab = std::fabs(tb);
	double acmx, acmn;
	if (std::fabs(a) > std::fabs(c)) {
		acmx = a;
		acmn = c;
	}
	else {
		acmx = c;
		acmn = a;
	}
	double rt;
	if (adf > ab) {
		rt = adf * std::sqrt(1.0 + (ab / adf) * (ab / adf));
	}
	else if (adf < ab) {
		rt = ab * std::sqrt(1.0 + (adf / ab) * (adf / ab));
	}
	else {
		rt = ab * std::sqrt(2.0);
	}
	int sgn1;
	if (sm < 0.0) {
		rt1 = 0.5 * (sm - rt);
		sgn1 = -1;
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else if (sm > 0.0) {
		rt1 = 0.5 * (sm + rt);
		sgn1 = 1;
		rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
	}
	else {
		rt1 = 0.5 * rt;
		rt2 = -0.5 * rt;
		sgn1 = 1;
	}
	double cs;
	int sgn2;
	if (df >= 0.0) {
		cs = df + rt;
		sgn2 = 1;
	}
	else {
		cs = df - rt;
		sgn2 = -1;
	}
	const double acs = std::fabs(cs);
	if (acs > ab) {
		const double ct = -tb / cs;
		sn1 = 1.0 / std::sqrt(1.0 + ct * ct);
		cs1 = ct * sn1;
	}
	else {
		if (ab == 0.0) {
			cs1 = 1.0;
			sn1 = 0.0;
		}
		else {
			const double tn = -cs / tb;
			cs1 = 1.0 / std::sqrt(1.0 + tn * tn);
			sn1 = tn * cs1;
		}
	}
	if (sgn1 == sgn2) {
		const double tn = cs1;
		cs1 = -sn1;
		sn1 = tn;
	}
}

// 2x2 上三角行列 [f g; 0 h] の特異値 (ssmin <= ssmax)
inline void rdlas2(const double f, const double g, const double h, double& ssmin, double& ssmax, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double fa = std::fabs(f);
	const double ga = std::fabs(g);
	const double ha = std::fabs(h);
	const double fhmn = std::min(fa, ha);
	const double fhmx = std::max(fa, ha);
	if (fhmn == 0.0) {
		ssmin = 0.0;
		if (fhmx == 0.0) {
			ssmax = ga;
		}
		else {
			const double mn = std::min(fhmx, ga);
			const double mx = std::max(fhmx, ga);
			ssmax = mx * std::sqrt(1.0 + (mn / mx) * (mn / mx));
		}
	}
	else {
		if (ga < fhmx) {
			const double as = 1.0 + fhmn / fhmx;
			const double at = (fhmx - fhmn) / fhmx;
			const double au = (ga / fhmx) * (ga / fhmx);
			const double c = 2.0 / (std::sqrt(as * as + au) + std::sqrt(at * at + au));
			ssmin = fhmn * c;
			ssmax = fhmx / c;
		}
		else {
			const double au = fhmx / ga;
			if (au == 0.0) {
				ssmin = (fhmn * fhmx) / ga;
				ssmax = ga;
			}
			else {
				const double as = 1.0 + fhmn / fhmx;
				const double at = (fhmx - fhmn) / fhmx;
				const double c = 1.0 / (std::sqrt(1.0 + (as * au) * (as * au)) +
					std::sqrt(1.0 + (at * au) * (at * au)));
				ssmin = (fhmn * c) * au;
				ssmin = ssmin + ssmin;
				ssmax = ga / (c + c);
			}
		}
	}
}

// 2x2 上三角行列 [f g; 0 h] の SVD (左右の回転付き)
inline void rdlasv2(const double f, const double g, const double h,
	double& ssmin, double& ssmax, double& snr, double& csr, double& snl, double& csl,
	const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	double ft = f;
	double fa = std::fabs(ft);
	double ht = h;
	double ha = std::fabs(h);
	int pmax = 1;
	const bool swap = (ha > fa);
	if (swap) {
		pmax = 3;
		double temp = ft;
		ft = ht;
		ht = temp;
		temp = fa;
		fa = ha;
		ha = temp;
	}
	const double gt = g;
	const double ga = std::fabs(gt);
	double clt = 0.0, crt = 0.0, slt = 0.0, srt = 0.0;
	if (ga == 0.0) {
		ssmin = ha;
		ssmax = fa;
		clt = 1.0;
		crt = 1.0;
		slt = 0.0;
		srt = 0.0;
	}
	else {
		bool gasmal = true;
		if (ga > fa) {
			pmax = 2;
			if ((fa / ga) < (DBL_EPSILON * 0.5)) {
				gasmal = false;
				ssmax = ga;
				if (ha > 1.0) {
					ssmin = fa / (ga / ha);
				}
				else {
					ssmin = (fa / ga) * ha;
				}
				clt = 1.0;
				slt = ht / gt;
				srt = 1.0;
				crt = ft / gt;
			}
		}
		if (gasmal) {
			const double d = fa - ha;
			double l;
			if (d == fa) {
				l = 1.0;
			}
			else {
				l = d / fa;
			}
			const double m = gt / ft;
			double t = 2.0 - l;
			const double mm = m * m;
			const double tt = t * t;
			const double s = std::sqrt(tt + mm);
			double r;
			if (l == 0.0) {
				r = std::fabs(m);
			}
			else {
				r = std::sqrt(l * l + mm);
			}
			const double a = 0.5 * (s + r);
			ssmin = ha / a;
			ssmax = fa * a;
			if (mm == 0.0) {
				if (l == 0.0) {
					t = detail::f_sign(2.0, ft) * detail::f_sign(1.0, gt);
				}
				else {
					t = gt / detail::f_sign(d, ft) + m / t;
				}
			}
			else {
				t = (m / (s + t) + m / (r + l)) * (1.0 + a);
			}
			const double ll = std::sqrt(t * t + 4.0);
			crt = 2.0 / ll;
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
	double tsign = 0.0;
	if (pmax == 1) {
		tsign = detail::f_sign(1.0, csr) * detail::f_sign(1.0, csl) * detail::f_sign(1.0, f);
	}
	if (pmax == 2) {
		tsign = detail::f_sign(1.0, snr) * detail::f_sign(1.0, csl) * detail::f_sign(1.0, g);
	}
	if (pmax == 3) {
		tsign = detail::f_sign(1.0, snr) * detail::f_sign(1.0, snl) * detail::f_sign(1.0, h);
	}
	ssmax = detail::f_sign(ssmax, tsign);
	ssmin = detail::f_sign(ssmin, tsign * detail::f_sign(1.0, f) * detail::f_sign(1.0, h));
}

// 2x2 行列 [a b; c d] を実 Schur 標準形へ変換する回転 (cs, sn) と固有値
inline void rdlanv2(double& a, double& b, double& c, double& d,
	double& rt1r, double& rt1i, double& rt2r, double& rt2i, double& cs, double& sn,
	const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double eps = DBL_EPSILON;
	// safmn2 = 2^int(log(safmin/eps)/log(2)/2) = 2^-485 (reference LAPACK と同値)
	const double safmn2 = std::ldexp(1.0, -485);
	const double safmx2 = 1.0 / safmn2;
	if (c == 0.0) {
		cs = 1.0;
		sn = 0.0;
	}
	else if (b == 0.0) {
		cs = 0.0;
		sn = 1.0;
		const double temp = d;
		d = a;
		a = temp;
		b = -c;
		c = 0.0;
	}
	else if ((a - d) == 0.0 && detail::f_sign(1.0, b) != detail::f_sign(1.0, c)) {
		cs = 1.0;
		sn = 0.0;
	}
	else {
		double temp = a - d;
		double p = 0.5 * temp;
		const double bcmax = std::max(std::fabs(b), std::fabs(c));
		const double bcmis = std::min(std::fabs(b), std::fabs(c)) * detail::f_sign(1.0, b) * detail::f_sign(1.0, c);
		double scale = std::max(std::fabs(p), bcmax);
		double z = (p / scale) * p + (bcmax / scale) * bcmis;
		if (z >= 4.0 * eps) {
			z = p + detail::f_sign(std::sqrt(scale) * std::sqrt(z), p);
			a = d + z;
			d = d - (bcmax / z) * bcmis;
			const double tau = rdlapy2(c, z, rounding_mode);
			cs = z / tau;
			sn = c / tau;
			b = b - c;
			c = 0.0;
		}
		else {
			int count = 0;
			double sigma = b + c;
			for (;;) {
				count++;
				scale = std::max(std::fabs(temp), std::fabs(sigma));
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
			p = 0.5 * temp;
			const double tau = rdlapy2(sigma, temp, rounding_mode);
			cs = std::sqrt(0.5 * (1.0 + std::fabs(sigma) / tau));
			sn = -(p / (tau * cs)) * detail::f_sign(1.0, sigma);
			const double aa = a * cs + b * sn;
			const double bb = -a * sn + b * cs;
			const double cc = c * cs + d * sn;
			const double dd = -c * sn + d * cs;
			a = aa * cs + cc * sn;
			b = (bb * cs) + (dd * sn);
			c = -(aa * sn) + (cc * cs);
			d = -bb * sn + dd * cs;
			temp = 0.5 * (a + d);
			a = temp;
			d = temp;
			if (c != 0.0) {
				if (b != 0.0) {
					if (detail::f_sign(1.0, b) == detail::f_sign(1.0, c)) {
						const double sab = std::sqrt(std::fabs(b));
						const double sac = std::sqrt(std::fabs(c));
						p = detail::f_sign(sab * sac, c);
						const double tau2 = 1.0 / std::sqrt(std::fabs(b + c));
						a = temp + p;
						d = temp - p;
						b = b - c;
						c = 0.0;
						const double cs1 = sab * tau2;
						const double sn1 = sac * tau2;
						temp = cs * cs1 - sn * sn1;
						sn = cs * sn1 + sn * cs1;
						cs = temp;
					}
				}
				else {
					b = -c;
					c = 0.0;
					temp = cs;
					cs = -sn;
					sn = temp;
				}
			}
		}
	}
	rt1r = a;
	rt2r = d;
	if (c == 0.0) {
		rt1i = 0.0;
		rt2i = 0.0;
	}
	else {
		rt1i = std::sqrt(std::fabs(b)) * std::sqrt(std::fabs(c));
		rt2i = -rt1i;
	}
}

// A := A * (cto/cfrom) (overflow/underflow を避けて段階的に乗算)
// type: 'G' 全体, 'L' 下三角, 'U' 上三角, 'H' Hessenberg,
//       'B' 対称帯下三角格納, 'Q' 対称帯上三角格納, 'Z' 帯 (LU 用)
inline void rdlascl(const char type, const int kl, const int ku, const double cfrom, const double cto,
	const int m, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
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
		detail::rdlapack_error("rdlascl: invalid type");
		return;
	}
	if (cfrom == 0.0 || std::isnan(cfrom) || std::isnan(cto)) {
		detail::rdlapack_error("rdlascl: invalid cfrom/cto");
	}
	if (m < 0 || n < 0) {
		detail::rdlapack_error("rdlascl: invalid m/n");
	}
	if (n == 0 || m == 0) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const double smlnum = DBL_MIN;
	const double bignum = 1.0 / smlnum;
	double cfromc = cfrom;
	double ctoc = cto;
	bool done = false;
	while (!done) {
		double mul;
		const double cfrom1 = cfromc * smlnum;
		if (cfrom1 == cfromc) {
			// cfromc は inf: mul は正規の値
			mul = ctoc / cfromc;
			done = true;
		}
		else {
			const double cto1 = ctoc / bignum;
			if (cto1 == ctoc) {
				// ctoc は 0 または inf
				mul = ctoc;
				done = true;
				cfromc = 1.0;
			}
			else if (std::fabs(cfrom1) > std::fabs(ctoc) && ctoc != 0.0) {
				mul = smlnum;
				done = false;
				cfromc = cfrom1;
			}
			else if (std::fabs(cto1) > std::fabs(cfromc)) {
				mul = bignum;
				done = false;
				ctoc = cto1;
			}
			else {
				mul = ctoc / cfromc;
				done = true;
				if (mul == 1.0) {
					return;
				}
			}
		}
		if (itype == 0) {
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i < m; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 1) {
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = j; i < m; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 2) {
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i <= std::min(j, m - 1); i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 3) {
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i <= std::min(j + 1, m - 1); i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 4) {
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = 0; i <= std::min(kl, n - 1 - j); i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else if (itype == 5) {
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
				for (int i = std::max(ku - j, 0); i <= ku; i++) {
					a[i] = a[i] * mul;
				}
			}
		}
		else { // itype == 6
			for (int j = 0; j < n; j++) {
				double* a = A + static_cast<std::size_t>(lda) * j;
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
inline void rdlasr(const char side, const char pivot, const char direct, const int m, const int n,
	const double* c, const double* s, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (!(detail::option_is(side, 'L') || detail::option_is(side, 'R'))) {
		detail::rdlapack_error("rdlasr: invalid side");
	}
	if (!(detail::option_is(pivot, 'V') || detail::option_is(pivot, 'T') || detail::option_is(pivot, 'B'))) {
		detail::rdlapack_error("rdlasr: invalid pivot");
	}
	if (!(detail::option_is(direct, 'F') || detail::option_is(direct, 'B'))) {
		detail::rdlapack_error("rdlasr: invalid direct");
	}
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdlasr: invalid m/n/lda");
	}
	if (m == 0 || n == 0) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (detail::option_is(side, 'L')) {
		if (detail::option_is(pivot, 'V')) {
			const int j0 = detail::option_is(direct, 'F') ? 0 : m - 2;
			const int j1 = detail::option_is(direct, 'F') ? m - 1 : -1;
			const int dj = detail::option_is(direct, 'F') ? 1 : -1;
			for (int j = j0; j != j1; j += dj) {
				const double ctemp = c[j];
				const double stemp = s[j];
				if (ctemp != 1.0 || stemp != 0.0) {
					for (int i = 0; i < n; i++) {
						const double temp = A[(j + 1) + L * i];
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
				const double ctemp = c[j - 1];
				const double stemp = s[j - 1];
				if (ctemp != 1.0 || stemp != 0.0) {
					for (int i = 0; i < n; i++) {
						const double temp = A[j + L * i];
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
				const double ctemp = c[j];
				const double stemp = s[j];
				if (ctemp != 1.0 || stemp != 0.0) {
					for (int i = 0; i < n; i++) {
						const double temp = A[j + L * i];
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
				const double ctemp = c[j];
				const double stemp = s[j];
				if (ctemp != 1.0 || stemp != 0.0) {
					for (int i = 0; i < m; i++) {
						const double temp = A[i + L * (j + 1)];
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
				const double ctemp = c[j - 1];
				const double stemp = s[j - 1];
				if (ctemp != 1.0 || stemp != 0.0) {
					for (int i = 0; i < m; i++) {
						const double temp = A[i + L * j];
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
				const double ctemp = c[j];
				const double stemp = s[j];
				if (ctemp != 1.0 || stemp != 0.0) {
					for (int i = 0; i < m; i++) {
						const double temp = A[i + L * j];
						A[i + L * j] = stemp * A[i + L * (n - 1)] + ctemp * temp;
						A[i + L * (n - 1)] = ctemp * A[i + L * (n - 1)] - stemp * temp;
					}
				}
			}
		}
	}
}

namespace vblas_rdlapack_detail {

// LAPACK の norm 系で使う max 更新 (NaN を伝播させる)
inline void norm_max(double& value, const double sum) {
	if (value < sum || std::isnan(sum)) {
		value = sum;
	}
}

} // namespace vblas_rdlapack_detail

// 一般行列の norm ('M' max|a_ij|, '1'/'O' 1-norm, 'I' inf-norm, 'F'/'E' Frobenius)
inline double rdlange(const char norm, const int m, const int n, const double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (std::min(m, n) == 0) {
		return 0.0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	double value = 0.0;
	if (detail::option_is(norm, 'M')) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				detail::norm_max(value, std::fabs(A[i + L * j]));
			}
		}
	}
	else if (detail::option_is(norm, 'O') || norm == '1') {
		for (int j = 0; j < n; j++) {
			double sum = 0.0;
			for (int i = 0; i < m; i++) {
				sum += std::fabs(A[i + L * j]);
			}
			detail::norm_max(value, sum);
		}
	}
	else if (detail::option_is(norm, 'I')) {
		std::vector<double> work(static_cast<std::size_t>(m), 0.0);
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				work[i] += std::fabs(A[i + L * j]);
			}
		}
		for (int i = 0; i < m; i++) {
			detail::norm_max(value, work[i]);
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		double scale = 0.0;
		double sum = 1.0;
		for (int j = 0; j < n; j++) {
			rdlassq(m, A + L * j, 1, scale, sum, rounding_mode);
		}
		value = scale * std::sqrt(sum);
	}
	else {
		detail::rdlapack_error("rdlange: invalid norm");
	}
	return value;
}

// 対称行列の norm
inline double rdlansy(const char norm, const char uplo, const int n, const double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n == 0) {
		return 0.0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const bool upper = detail::option_is(uplo, 'U');
	double value = 0.0;
	if (detail::option_is(norm, 'M')) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= j; i++) {
					detail::norm_max(value, std::fabs(A[i + L * j]));
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				for (int i = j; i < n; i++) {
					detail::norm_max(value, std::fabs(A[i + L * j]));
				}
			}
		}
	}
	else if (detail::option_is(norm, 'I') || detail::option_is(norm, 'O') || norm == '1') {
		std::vector<double> work(static_cast<std::size_t>(n), 0.0);
		if (upper) {
			for (int j = 0; j < n; j++) {
				double sum = 0.0;
				for (int i = 0; i < j; i++) {
					const double absa = std::fabs(A[i + L * j]);
					sum += absa;
					work[i] += absa;
				}
				work[j] = sum + std::fabs(A[j + L * j]);
			}
			for (int i = 0; i < n; i++) {
				detail::norm_max(value, work[i]);
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				double sum = work[j] + std::fabs(A[j + L * j]);
				for (int i = j + 1; i < n; i++) {
					const double absa = std::fabs(A[i + L * j]);
					sum += absa;
					work[i] += absa;
				}
				detail::norm_max(value, sum);
			}
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		double scale = 0.0;
		double sum = 1.0;
		if (upper) {
			for (int j = 1; j < n; j++) {
				rdlassq(j, A + L * j, 1, scale, sum, rounding_mode);
			}
		}
		else {
			for (int j = 0; j < n - 1; j++) {
				rdlassq(n - j - 1, A + (j + 1) + L * j, 1, scale, sum, rounding_mode);
			}
		}
		sum = 2.0 * sum;
		rdlassq(n, A, lda + 1, scale, sum, rounding_mode);
		value = scale * std::sqrt(sum);
	}
	else {
		detail::rdlapack_error("rdlansy: invalid norm");
	}
	return value;
}

// 対称三重対角行列の norm (d: 対角 n 個, e: 副対角 n-1 個)
inline double rdlanst(const char norm, const int n, const double* d, const double* e, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n <= 0) {
		return 0.0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	double value = 0.0;
	if (detail::option_is(norm, 'M')) {
		value = std::fabs(d[n - 1]);
		for (int i = 0; i < n - 1; i++) {
			detail::norm_max(value, std::fabs(d[i]));
			detail::norm_max(value, std::fabs(e[i]));
		}
	}
	else if (detail::option_is(norm, 'O') || norm == '1' || detail::option_is(norm, 'I')) {
		if (n == 1) {
			value = std::fabs(d[0]);
		}
		else {
			value = std::fabs(d[0]) + std::fabs(e[0]);
			detail::norm_max(value, std::fabs(e[n - 2]) + std::fabs(d[n - 1]));
			for (int i = 1; i < n - 1; i++) {
				detail::norm_max(value, std::fabs(d[i]) + std::fabs(e[i]) + std::fabs(e[i - 1]));
			}
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		double scale = 0.0;
		double sum = 1.0;
		if (n > 1) {
			rdlassq(n - 1, e, 1, scale, sum, rounding_mode);
			sum = 2.0 * sum;
		}
		rdlassq(n, d, 1, scale, sum, rounding_mode);
		value = scale * std::sqrt(sum);
	}
	else {
		detail::rdlapack_error("rdlanst: invalid norm");
	}
	return value;
}

// 一般帯行列 (n x n, 帯幅 kl/ku, 帯格納) の norm
inline double rdlangb(const char norm, const int n, const int kl, const int ku, const double* AB, const int ldab, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n == 0) {
		return 0.0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(ldab);
	double value = 0.0;
	if (detail::option_is(norm, 'M')) {
		for (int j = 0; j < n; j++) {
			for (int i = std::max(ku - j, 0); i <= std::min(n + ku - 1 - j, kl + ku); i++) {
				detail::norm_max(value, std::fabs(AB[i + L * j]));
			}
		}
	}
	else if (detail::option_is(norm, 'O') || norm == '1') {
		for (int j = 0; j < n; j++) {
			double sum = 0.0;
			for (int i = std::max(ku - j, 0); i <= std::min(n + ku - 1 - j, kl + ku); i++) {
				sum += std::fabs(AB[i + L * j]);
			}
			detail::norm_max(value, sum);
		}
	}
	else if (detail::option_is(norm, 'I')) {
		std::vector<double> work(static_cast<std::size_t>(n), 0.0);
		for (int j = 0; j < n; j++) {
			for (int i = std::max(0, j - ku); i <= std::min(n - 1, j + kl); i++) {
				work[i] += std::fabs(AB[(ku + i - j) + L * j]);
			}
		}
		for (int i = 0; i < n; i++) {
			detail::norm_max(value, work[i]);
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		double scale = 0.0;
		double sum = 1.0;
		for (int j = 0; j < n; j++) {
			const int l = std::max(0, j - ku);
			const int k = ku - j + l;
			rdlassq(std::min(n - 1, j + kl) - l + 1, AB + k + L * j, 1, scale, sum, rounding_mode);
		}
		value = scale * std::sqrt(sum);
	}
	else {
		detail::rdlapack_error("rdlangb: invalid norm");
	}
	return value;
}

// 対称帯行列 (帯幅 k, 帯格納) の norm
inline double rdlansb(const char norm, const char uplo, const int n, const int k, const double* AB, const int ldab, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n == 0) {
		return 0.0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(ldab);
	const bool upper = detail::option_is(uplo, 'U');
	double value = 0.0;
	if (detail::option_is(norm, 'M')) {
		if (upper) {
			for (int j = 0; j < n; j++) {
				for (int i = std::max(k - j, 0); i <= k; i++) {
					detail::norm_max(value, std::fabs(AB[i + L * j]));
				}
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= std::min(n - 1 - j, k); i++) {
					detail::norm_max(value, std::fabs(AB[i + L * j]));
				}
			}
		}
	}
	else if (detail::option_is(norm, 'I') || detail::option_is(norm, 'O') || norm == '1') {
		std::vector<double> work(static_cast<std::size_t>(n), 0.0);
		if (upper) {
			for (int j = 0; j < n; j++) {
				double sum = 0.0;
				for (int i = std::max(0, j - k); i < j; i++) {
					const double absa = std::fabs(AB[(k + i - j) + L * j]);
					sum += absa;
					work[i] += absa;
				}
				work[j] = sum + std::fabs(AB[k + L * j]);
			}
			for (int i = 0; i < n; i++) {
				detail::norm_max(value, work[i]);
			}
		}
		else {
			for (int j = 0; j < n; j++) {
				double sum = work[j] + std::fabs(AB[0 + L * j]);
				for (int i = j + 1; i <= std::min(n - 1, j + k); i++) {
					const double absa = std::fabs(AB[(i - j) + L * j]);
					sum += absa;
					work[i] += absa;
				}
				detail::norm_max(value, sum);
			}
		}
	}
	else if (detail::option_is(norm, 'F') || detail::option_is(norm, 'E')) {
		double scale = 0.0;
		double sum = 1.0;
		if (k > 0) {
			if (upper) {
				for (int j = 1; j < n; j++) {
					rdlassq(std::min(j, k), AB + std::max(k - j, 0) + L * j, 1, scale, sum, rounding_mode);
				}
			}
			else {
				for (int j = 0; j < n - 1; j++) {
					rdlassq(std::min(n - 1 - j, k), AB + 1 + L * j, 1, scale, sum, rounding_mode);
				}
			}
			sum = 2.0 * sum;
		}
		rdlassq(n, AB + (upper ? k : 0), ldab, scale, sum, rounding_mode);
		value = scale * std::sqrt(sum);
	}
	else {
		detail::rdlapack_error("rdlansb: invalid norm");
	}
	return value;
}

// Householder reflector の生成: H * [alpha; x] = [beta; 0]
// alpha は beta に上書き，x は v(2:n) に上書きされる
inline void rdlarfg(const int n, double& alpha, double* x, const int incx, double& tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n <= 1) {
		tau = 0.0;
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	double xnorm = rdnrm2(n - 1, x, incx, rounding_mode);
	if (xnorm == 0.0) {
		tau = 0.0;
	}
	else {
		double beta = -detail::f_sign(rdlapy2(alpha, xnorm, rounding_mode), alpha);
		const double safmin = DBL_MIN / (DBL_EPSILON * 0.5);
		int knt = 0;
		if (std::fabs(beta) < safmin) {
			const double rsafmn = 1.0 / safmin;
			do {
				knt++;
				rdscal(n - 1, rsafmn, x, incx, rounding_mode);
				beta = beta * rsafmn;
				alpha = alpha * rsafmn;
			} while (std::fabs(beta) < safmin && knt < 20);
			xnorm = rdnrm2(n - 1, x, incx, rounding_mode);
			beta = -detail::f_sign(rdlapy2(alpha, xnorm, rounding_mode), alpha);
		}
		tau = (beta - alpha) / beta;
		rdscal(n - 1, 1.0 / (alpha - beta), x, incx, rounding_mode);
		for (int j = 0; j < knt; j++) {
			beta = beta * safmin;
		}
		alpha = beta;
	}
}

// Householder reflector の適用: C := H*C (side='L') または C*H (side='R'),
// H = I - tau*v*v^T
inline void rdlarf(const char side, const int m, const int n, const double* v, const int incv,
	const double tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const bool applyleft = detail::option_is(side, 'L');
	int lastv = 0;
	int lastc = 0;
	if (tau != 0.0) {
		lastv = applyleft ? m : n;
		int i = (incv > 0) ? (lastv - 1) * incv : 0;
		while (lastv > 0 && v[i] == 0.0) {
			lastv--;
			i -= incv;
		}
		if (applyleft) {
			lastc = iladlc(lastv, n, C, ldc);
		}
		else {
			lastc = iladlr(m, lastv, C, ldc);
		}
	}
	if (lastv <= 0) {
		return;
	}
	std::vector<double> work(static_cast<std::size_t>(lastc) + 1);
	if (applyleft) {
		rdgemv('T', lastv, lastc, 1.0, C, ldc, v, incv, 0.0, work.data(), 1, rounding_mode);
		rdger(lastv, lastc, -tau, v, incv, work.data(), 1, C, ldc, rounding_mode);
	}
	else {
		rdgemv('N', lastc, lastv, 1.0, C, ldc, v, incv, 0.0, work.data(), 1, rounding_mode);
		rdger(lastc, lastv, -tau, work.data(), 1, v, incv, C, ldc, rounding_mode);
	}
}

// rdlarf の v(1) = 1 を暗黙とする版 (v の先頭要素は参照されない)
inline void rdlarf1f(const char side, const int m, const int n, const double* v, const int incv,
	const double tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const bool applyleft = detail::option_is(side, 'L');
	int lastv = 1;
	int lastc = 0;
	if (tau != 0.0) {
		lastv = applyleft ? m : n;
		int i = (incv > 0) ? (lastv - 1) * incv : 0;
		while (lastv > 1 && v[i] == 0.0) {
			lastv--;
			i -= incv;
		}
		if (applyleft) {
			lastc = iladlc(lastv, n, C, ldc);
		}
		else {
			lastc = iladlr(m, lastv, C, ldc);
		}
	}
	if (lastc == 0) {
		return;
	}
	if (applyleft) {
		if (lastv == 1) {
			rdscal(lastc, 1.0 - tau, C, ldc, rounding_mode);
		}
		else {
			std::vector<double> work(static_cast<std::size_t>(lastc));
			rdgemv('T', lastv - 1, lastc, 1.0, C + 1, ldc, v + incv, incv, 0.0, work.data(), 1, rounding_mode);
			rdaxpy(lastc, 1.0, C, ldc, work.data(), 1, rounding_mode);
			rdaxpy(lastc, -tau, work.data(), 1, C, ldc, rounding_mode);
			rdger(lastv - 1, lastc, -tau, v + incv, incv, work.data(), 1, C + 1, ldc, rounding_mode);
		}
	}
	else {
		if (lastv == 1) {
			rdscal(lastc, 1.0 - tau, C, 1, rounding_mode);
		}
		else {
			std::vector<double> work(static_cast<std::size_t>(lastc));
			rdgemv('N', lastc, lastv - 1, 1.0, C + static_cast<std::size_t>(ldc), ldc, v + incv, incv, 0.0, work.data(), 1, rounding_mode);
			rdaxpy(lastc, 1.0, C, 1, work.data(), 1, rounding_mode);
			rdaxpy(lastc, -tau, work.data(), 1, C, 1, rounding_mode);
			rdger(lastc, lastv - 1, -tau, work.data(), 1, v + incv, incv, C + static_cast<std::size_t>(ldc), ldc, rounding_mode);
		}
	}
}

// rdlarf の v(last) = 1 を暗黙とする版 (v の末尾要素は参照されない)
inline void rdlarf1l(const char side, const int m, const int n, const double* v, const int incv,
	const double tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const bool applyleft = detail::option_is(side, 'L');
	int firstv = 1;
	int lastv = 0;
	int lastc = 0;
	int i = 0;
	if (tau != 0.0) {
		lastv = applyleft ? m : n;
		while (lastv > firstv && v[i] == 0.0) {
			firstv++;
			i += incv;
		}
		if (applyleft) {
			lastc = iladlc(lastv, n, C, ldc);
		}
		else {
			lastc = iladlr(m, lastv, C, ldc);
		}
	}
	if (lastc == 0 || lastv <= 0) {
		return;
	}
	const std::size_t L = static_cast<std::size_t>(ldc);
	if (applyleft) {
		if (lastv == firstv) {
			rdscal(lastc, 1.0 - tau, C + (firstv - 1), ldc, rounding_mode);
		}
		else {
			std::vector<double> work(static_cast<std::size_t>(lastc));
			rdgemv('T', lastv - firstv, lastc, 1.0, C + (firstv - 1), ldc, v + i, incv, 0.0, work.data(), 1, rounding_mode);
			rdaxpy(lastc, 1.0, C + (lastv - 1), ldc, work.data(), 1, rounding_mode);
			rdaxpy(lastc, -tau, work.data(), 1, C + (lastv - 1), ldc, rounding_mode);
			rdger(lastv - firstv, lastc, -tau, v + i, incv, work.data(), 1, C + (firstv - 1), ldc, rounding_mode);
		}
	}
	else {
		if (lastv == firstv) {
			// H = I - tau*e*e^T (e は第 lastv 単位 vector): 列 lastv のみ scale
			// (reference 3.12.1 は C(1,1) を scale するが，正しくは列 lastv)
			rdscal(lastc, 1.0 - tau, C + L * (lastv - 1), 1, rounding_mode);
		}
		else {
			std::vector<double> work(static_cast<std::size_t>(lastc));
			rdgemv('N', lastc, lastv - firstv, 1.0, C + L * (firstv - 1), ldc, v + i, incv, 0.0, work.data(), 1, rounding_mode);
			rdaxpy(lastc, 1.0, C + L * (lastv - 1), 1, work.data(), 1, rounding_mode);
			rdaxpy(lastc, -tau, work.data(), 1, C + L * (lastv - 1), 1, rounding_mode);
			rdger(lastc, lastv - firstv, -tau, work.data(), 1, v + i, incv, C + L * (firstv - 1), ldc, rounding_mode);
		}
	}
}

// block reflector の三角因子 T の生成 (recursive, reference LAPACK 3.12.1)
// direct: 'F'/'B', storev: 'C'/'R'
inline void rdlarft(const char direct, const char storev, const int n, const int k,
	const double* V, const int ldv, const double* tau, double* T, const int ldt, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n == 0 || k == 0) {
		return;
	}
	if (n == 1 || k == 1) {
		T[0] = tau[0];
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t LV = static_cast<std::size_t>(ldv);
	const std::size_t LT = static_cast<std::size_t>(ldt);
	const int l = k / 2;
	const bool dirf = detail::option_is(direct, 'F');
	const bool colv = detail::option_is(storev, 'C');
	if (dirf && colv) { // QR
		rdlarft(direct, storev, n, l, V, ldv, tau, T, ldt, rounding_mode);
		rdlarft(direct, storev, n - l, k - l, V + l + LV * l, ldv, tau + l, T + l + LT * l, ldt, rounding_mode);
		for (int j = 0; j < l; j++) {
			for (int i = 0; i < k - l; i++) {
				T[j + LT * (l + i)] = V[(l + i) + LV * j];
			}
		}
		rdtrmm('R', 'L', 'N', 'U', l, k - l, 1.0, V + l + LV * l, ldv, T + LT * l, ldt, rounding_mode);
		if (n > k) {
			rdgemm('T', 'N', l, k - l, n - k, 1.0, V + k, ldv, V + k + LV * l, ldv, 1.0, T + LT * l, ldt, rounding_mode);
		}
		rdtrmm('L', 'U', 'N', 'N', l, k - l, -1.0, T, ldt, T + LT * l, ldt, rounding_mode);
		rdtrmm('R', 'U', 'N', 'N', l, k - l, 1.0, T + l + LT * l, ldt, T + LT * l, ldt, rounding_mode);
	}
	else if (dirf && !colv) { // LQ
		rdlarft(direct, storev, n, l, V, ldv, tau, T, ldt, rounding_mode);
		rdlarft(direct, storev, n - l, k - l, V + l + LV * l, ldv, tau + l, T + l + LT * l, ldt, rounding_mode);
		dlacpy('A', l, k - l, V + LV * l, ldv, T + LT * l, ldt);
		rdtrmm('R', 'U', 'T', 'U', l, k - l, 1.0, V + l + LV * l, ldv, T + LT * l, ldt, rounding_mode);
		if (n > k) {
			rdgemm('N', 'T', l, k - l, n - k, 1.0, V + LV * k, ldv, V + l + LV * k, ldv, 1.0, T + LT * l, ldt, rounding_mode);
		}
		rdtrmm('L', 'U', 'N', 'N', l, k - l, -1.0, T, ldt, T + LT * l, ldt, rounding_mode);
		rdtrmm('R', 'U', 'N', 'N', l, k - l, 1.0, T + l + LT * l, ldt, T + LT * l, ldt, rounding_mode);
	}
	else if (!dirf && colv) { // QL
		rdlarft(direct, storev, n - l, k - l, V, ldv, tau, T, ldt, rounding_mode);
		rdlarft(direct, storev, n, l, V + LV * (k - l), ldv, tau + (k - l), T + (k - l) + LT * (k - l), ldt, rounding_mode);
		for (int j = 0; j < k - l; j++) {
			for (int i = 0; i < l; i++) {
				T[(k - l + i) + LT * j] = V[(n - k + j) + LV * (k - l + i)];
			}
		}
		rdtrmm('R', 'U', 'N', 'U', l, k - l, 1.0, V + (n - k), ldv, T + (k - l), ldt, rounding_mode);
		if (n > k) {
			rdgemm('T', 'N', l, k - l, n - k, 1.0, V + LV * (k - l), ldv, V, ldv, 1.0, T + (k - l), ldt, rounding_mode);
		}
		rdtrmm('L', 'L', 'N', 'N', l, k - l, -1.0, T + (k - l) + LT * (k - l), ldt, T + (k - l), ldt, rounding_mode);
		rdtrmm('R', 'L', 'N', 'N', l, k - l, 1.0, T, ldt, T + (k - l), ldt, rounding_mode);
	}
	else { // RQ (direct='B', storev='R')
		rdlarft(direct, storev, n - l, k - l, V, ldv, tau, T, ldt, rounding_mode);
		rdlarft(direct, storev, n, l, V + (k - l), ldv, tau + (k - l), T + (k - l) + LT * (k - l), ldt, rounding_mode);
		dlacpy('A', l, k - l, V + (k - l) + LV * (n - k), ldv, T + (k - l), ldt);
		rdtrmm('R', 'L', 'T', 'U', l, k - l, 1.0, V + LV * (n - k), ldv, T + (k - l), ldt, rounding_mode);
		if (n > k) {
			rdgemm('N', 'T', l, k - l, n - k, 1.0, V + (k - l), ldv, V, ldv, 1.0, T + (k - l), ldt, rounding_mode);
		}
		rdtrmm('L', 'L', 'N', 'N', l, k - l, -1.0, T + (k - l) + LT * (k - l), ldt, T + (k - l), ldt, rounding_mode);
		rdtrmm('R', 'L', 'N', 'N', l, k - l, 1.0, T, ldt, T + (k - l), ldt, rounding_mode);
	}
}

// block reflector の適用: C := H*C, H^T*C, C*H, C*H^T
inline void rdlarfb(const char side, const char trans, const char direct, const char storev,
	const int m, const int n, const int k, const double* V, const int ldv,
	const double* T, const int ldt, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m <= 0 || n <= 0) {
		return;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const char transt = detail::option_is(trans, 'N') ? 'T' : 'N';
	const std::size_t LV = static_cast<std::size_t>(ldv);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const int ldwork = std::max(1, detail::option_is(side, 'L') ? n : m);
	const std::size_t LW = static_cast<std::size_t>(ldwork);
	std::vector<double> work(LW * static_cast<std::size_t>(k));
	double* W = work.data();
	const int rm = rounding_mode;
	if (detail::option_is(storev, 'C')) {
		if (detail::option_is(direct, 'F')) {
			if (detail::option_is(side, 'L')) {
				// W := C1^T * V1, C2 -= V2*W^T, ...
				for (int j = 0; j < k; j++) {
					dcopy(n, C + j, ldc, W + LW * j, 1);
				}
				rdtrmm('R', 'L', 'N', 'U', n, k, 1.0, V, ldv, W, ldwork, rm);
				if (m > k) {
					rdgemm('T', 'N', n, k, m - k, 1.0, C + k, ldc, V + k, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'U', transt, 'N', n, k, 1.0, T, ldt, W, ldwork, rm);
				if (m > k) {
					rdgemm('N', 'T', m - k, n, k, -1.0, V + k, ldv, W, ldwork, 1.0, C + k, ldc, rm);
				}
				rdtrmm('R', 'L', 'T', 'U', n, k, 1.0, V, ldv, W, ldwork, rm);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[j + LC * i] = C[j + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					dcopy(m, C + LC * j, 1, W + LW * j, 1);
				}
				rdtrmm('R', 'L', 'N', 'U', m, k, 1.0, V, ldv, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'N', m, k, n - k, 1.0, C + LC * k, ldc, V + k, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'U', trans, 'N', m, k, 1.0, T, ldt, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'T', m, n - k, k, -1.0, W, ldwork, V + k, ldv, 1.0, C + LC * k, ldc, rm);
				}
				rdtrmm('R', 'L', 'T', 'U', m, k, 1.0, V, ldv, W, ldwork, rm);
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
					dcopy(n, C + (m - k + j), ldc, W + LW * j, 1);
				}
				rdtrmm('R', 'U', 'N', 'U', n, k, 1.0, V + (m - k), ldv, W, ldwork, rm);
				if (m > k) {
					rdgemm('T', 'N', n, k, m - k, 1.0, C, ldc, V, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'L', transt, 'N', n, k, 1.0, T, ldt, W, ldwork, rm);
				if (m > k) {
					rdgemm('N', 'T', m - k, n, k, -1.0, V, ldv, W, ldwork, 1.0, C, ldc, rm);
				}
				rdtrmm('R', 'U', 'T', 'U', n, k, 1.0, V + (m - k), ldv, W, ldwork, rm);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[(m - k + j) + LC * i] = C[(m - k + j) + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					dcopy(m, C + LC * (n - k + j), 1, W + LW * j, 1);
				}
				rdtrmm('R', 'U', 'N', 'U', m, k, 1.0, V + (n - k), ldv, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'N', m, k, n - k, 1.0, C, ldc, V, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'L', trans, 'N', m, k, 1.0, T, ldt, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'T', m, n - k, k, -1.0, W, ldwork, V, ldv, 1.0, C, ldc, rm);
				}
				rdtrmm('R', 'U', 'T', 'U', m, k, 1.0, V + (n - k), ldv, W, ldwork, rm);
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
					dcopy(n, C + j, ldc, W + LW * j, 1);
				}
				rdtrmm('R', 'U', 'T', 'U', n, k, 1.0, V, ldv, W, ldwork, rm);
				if (m > k) {
					rdgemm('T', 'T', n, k, m - k, 1.0, C + k, ldc, V + LV * k, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'U', transt, 'N', n, k, 1.0, T, ldt, W, ldwork, rm);
				if (m > k) {
					rdgemm('T', 'T', m - k, n, k, -1.0, V + LV * k, ldv, W, ldwork, 1.0, C + k, ldc, rm);
				}
				rdtrmm('R', 'U', 'N', 'U', n, k, 1.0, V, ldv, W, ldwork, rm);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[j + LC * i] = C[j + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					dcopy(m, C + LC * j, 1, W + LW * j, 1);
				}
				rdtrmm('R', 'U', 'T', 'U', m, k, 1.0, V, ldv, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'T', m, k, n - k, 1.0, C + LC * k, ldc, V + LV * k, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'U', trans, 'N', m, k, 1.0, T, ldt, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'N', m, n - k, k, -1.0, W, ldwork, V + LV * k, ldv, 1.0, C + LC * k, ldc, rm);
				}
				rdtrmm('R', 'U', 'N', 'U', m, k, 1.0, V, ldv, W, ldwork, rm);
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
					dcopy(n, C + (m - k + j), ldc, W + LW * j, 1);
				}
				rdtrmm('R', 'L', 'T', 'U', n, k, 1.0, V + LV * (m - k), ldv, W, ldwork, rm);
				if (m > k) {
					rdgemm('T', 'T', n, k, m - k, 1.0, C, ldc, V, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'L', transt, 'N', n, k, 1.0, T, ldt, W, ldwork, rm);
				if (m > k) {
					rdgemm('T', 'T', m - k, n, k, -1.0, V, ldv, W, ldwork, 1.0, C, ldc, rm);
				}
				rdtrmm('R', 'L', 'N', 'U', n, k, 1.0, V + LV * (m - k), ldv, W, ldwork, rm);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < n; i++) {
						C[(m - k + j) + LC * i] = C[(m - k + j) + LC * i] - W[i + LW * j];
					}
				}
			}
			else {
				for (int j = 0; j < k; j++) {
					dcopy(m, C + LC * (n - k + j), 1, W + LW * j, 1);
				}
				rdtrmm('R', 'L', 'T', 'U', m, k, 1.0, V + LV * (n - k), ldv, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'T', m, k, n - k, 1.0, C, ldc, V, ldv, 1.0, W, ldwork, rm);
				}
				rdtrmm('R', 'L', trans, 'N', m, k, 1.0, T, ldt, W, ldwork, rm);
				if (n > k) {
					rdgemm('N', 'N', m, n - k, k, -1.0, W, ldwork, V, ldv, 1.0, C, ldc, rm);
				}
				rdtrmm('R', 'L', 'N', 'U', m, k, 1.0, V + LV * (n - k), ldv, W, ldwork, rm);
				for (int j = 0; j < k; j++) {
					for (int i = 0; i < m; i++) {
						C[i + LC * (n - k + j)] = C[i + LC * (n - k + j)] - W[i + LW * j];
					}
				}
			}
		}
	}
}

#endif // VBLAS_RDLAPACK_AUX_HPP
