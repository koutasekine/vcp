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
// element type T (no rounding-mode control, no fma).
// ----------------------------------------------------------------------------

#pragma once

#ifndef TBLAS_TBLAS_LEVEL1_HPP
#define TBLAS_TBLAS_LEVEL1_HPP

// Level 1 BLAS (template 版)．BLAS 名の d を t に読み替えた命名 (idamax -> itamax)．
// - itamax の返り値は C 流の 0-based index (Fortran BLAS の 1-based とは異なる)
// - abs / sqrt は `using std::abs;` + 非修飾呼び出しで ADL 解決する
//   (kv::dd / kv::mpfr / kv::interval は namespace kv の abs / sqrt が使われる)
// - itamax / tnrm2 / trotg / trotmg は順序比較で分岐するため，区間型では非推奨
//   (tblas/README_tblas.md 参照)

#include <cmath>
#include <cstddef>

#include "tblas_common.hpp"

// ---- 演算なし ----

// y := x
template <typename T>
inline void tcopy(const int n, const T* x, const int incx, T* y, const int incy) {
	if (n <= 0) {
		return;
	}
	int ix = tblas_detail::vec_start(n, incx);
	int iy = tblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		y[iy] = x[ix];
		ix += incx;
		iy += incy;
	}
}

// x <-> y
template <typename T>
inline void tswap(const int n, T* x, const int incx, T* y, const int incy) {
	if (n <= 0) {
		return;
	}
	int ix = tblas_detail::vec_start(n, incx);
	int iy = tblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		const T tmp = x[ix];
		x[ix] = y[iy];
		y[iy] = tmp;
		ix += incx;
		iy += incy;
	}
}

// |x_i| が最大となる 0-based index を返す (n<=0 または incx<=0 のとき -1)
template <typename T>
inline int itamax(const int n, const T* x, const int incx) {
	using std::abs;
	if (n <= 0 || incx <= 0) {
		return -1;
	}
	int result = 0;
	T tmax = abs(x[0]);
	int ix = incx;
	for (int i = 1; i < n; i++) {
		const T v = abs(x[ix]);
		if (v > tmax) {
			tmax = v;
			result = i;
		}
		ix += incx;
	}
	return result;
}

// ---- 演算あり ----

// x := alpha*x
template <typename T>
inline void tscal(const int n, const T& alpha, T* x, const int incx) {
	if (n <= 0 || incx <= 0 || alpha == T(1)) {
		return;
	}
	int ix = 0;
	for (int i = 0; i < n; i++) {
		x[ix] = alpha * x[ix];
		ix += incx;
	}
}

// y := alpha*x + y
template <typename T>
inline void taxpy(const int n, const T& alpha, const T* x, const int incx, T* y, const int incy) {
	if (n <= 0 || alpha == T(0)) {
		return;
	}
	int ix = tblas_detail::vec_start(n, incx);
	int iy = tblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		y[iy] = alpha * x[ix] + y[iy];
		ix += incx;
		iy += incy;
	}
}

// return x^T y (逐次和)
template <typename T>
inline T tdot(const int n, const T* x, const int incx, const T* y, const int incy) {
	if (n <= 0) {
		return T(0);
	}
	int ix = tblas_detail::vec_start(n, incx);
	int iy = tblas_detail::vec_start(n, incy);
	T s = T(0);
	for (int i = 0; i < n; i++) {
		s = s + x[ix] * y[iy];
		ix += incx;
		iy += incy;
	}
	return s;
}

// return sum |x_i|
template <typename T>
inline T tasum(const int n, const T* x, const int incx) {
	using std::abs;
	if (n <= 0 || incx <= 0) {
		return T(0);
	}
	T s = T(0);
	int ix = 0;
	for (int i = 0; i < n; i++) {
		s = s + abs(x[ix]);
		ix += incx;
	}
	return s;
}

// return ||x||_2 (overflow/underflow 回避の scale 付き，reference dnrm2 と同方式)
template <typename T>
inline T tnrm2(const int n, const T* x, const int incx) {
	using std::abs;
	using std::sqrt;
	if (n <= 0 || incx <= 0) {
		return T(0);
	}
	if (n == 1) {
		return abs(x[0]);
	}
	T scale = T(0);
	T ssq = T(1);
	int ix = 0;
	for (int i = 0; i < n; i++) {
		const T v = abs(x[ix]);
		if (v != T(0)) {
			if (scale < v) {
				const T r = scale / v;
				ssq = ssq * (r * r) + T(1);
				scale = v;
			}
			else {
				const T r = v / scale;
				ssq = r * r + ssq;
			}
		}
		ix += incx;
	}
	return scale * sqrt(ssq);
}

// 平面回転の適用: (x_i, y_i) := (c*x_i + s*y_i, c*y_i - s*x_i)
template <typename T>
inline void trot(const int n, T* x, const int incx, T* y, const int incy, const T& c, const T& s) {
	if (n <= 0) {
		return;
	}
	int ix = tblas_detail::vec_start(n, incx);
	int iy = tblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		const T tmp = c * x[ix] + s * y[iy];
		y[iy] = c * y[iy] - s * x[ix];
		x[ix] = tmp;
		ix += incx;
		iy += incy;
	}
}

// Givens 回転の生成 (reference drotg と同方式)
template <typename T>
inline void trotg(T* a, T* b, T* c, T* s) {
	using std::abs;
	using std::sqrt;
	T roe = *b;
	if (abs(*a) > abs(*b)) {
		roe = *a;
	}
	const T scale = abs(*a) + abs(*b);
	T r, z;
	if (scale == T(0)) {
		*c = T(1);
		*s = T(0);
		r = T(0);
		z = T(0);
	}
	else {
		const T as = *a / scale;
		const T bs = *b / scale;
		r = scale * sqrt(as * as + bs * bs);
		if (roe < T(0)) {
			r = -r;
		}
		*c = *a / r;
		*s = *b / r;
		z = T(1);
		if (abs(*a) > abs(*b)) {
			z = *s;
		}
		if (abs(*b) >= abs(*a) && *c != T(0)) {
			z = T(1) / *c;
		}
	}
	*a = r;
	*b = z;
}

// modified Givens 回転の適用 (reference drotm と同方式)
// param[0] (flag) は trotmg が生成する -2, -1, 0, 1 のいずれかを仮定する
template <typename T>
inline void trotm(const int n, T* x, const int incx, T* y, const int incy, const T* param) {
	const T flag = param[0];
	if (n <= 0 || flag == T(-2)) {
		return;
	}
	int ix = tblas_detail::vec_start(n, incx);
	int iy = tblas_detail::vec_start(n, incy);
	if (flag == T(0)) {
		const T h12 = param[3];
		const T h21 = param[2];
		for (int i = 0; i < n; i++) {
			const T w = x[ix];
			const T z = y[iy];
			x[ix] = w + h12 * z;
			y[iy] = h21 * w + z;
			ix += incx;
			iy += incy;
		}
	}
	else if (flag == T(1)) {
		const T h11 = param[1];
		const T h22 = param[4];
		for (int i = 0; i < n; i++) {
			const T w = x[ix];
			const T z = y[iy];
			x[ix] = h11 * w + z;
			y[iy] = h22 * z - w;
			ix += incx;
			iy += incy;
		}
	}
	else {
		const T h11 = param[1];
		const T h12 = param[3];
		const T h21 = param[2];
		const T h22 = param[4];
		for (int i = 0; i < n; i++) {
			const T w = x[ix];
			const T z = y[iy];
			x[ix] = h11 * w + h12 * z;
			y[iy] = h21 * w + h22 * z;
			ix += incx;
			iy += incy;
		}
	}
}

// modified Givens 回転の生成 (reference drotmg と同方式)
// rescaling 閾値は gam = 4096 から gamsq = gam^2, rgamsq = 1/gam^2 を T の演算で
// 構成する (reference の double 定数 5.9604645e-8 とは末位が僅かに異なる)
template <typename T>
inline void trotmg(T* d1, T* d2, T* x1, const T& y1, T* param) {
	using std::abs;
	const T gam = T(4096);
	const T gamsq = gam * gam;
	const T rgamsq = T(1) / gamsq;
	T flag, h11 = T(0), h12 = T(0), h21 = T(0), h22 = T(0);

	if (*d1 < T(0)) {
		flag = T(-1);
		*d1 = T(0);
		*d2 = T(0);
		*x1 = T(0);
	}
	else {
		const T p2 = *d2 * y1;
		if (p2 == T(0)) {
			param[0] = T(-2);
			return;
		}
		const T p1 = *d1 * *x1;
		const T q2 = p2 * y1;
		const T q1 = p1 * *x1;
		if (abs(q1) > abs(q2)) {
			h21 = -y1 / *x1;
			h12 = p2 / p1;
			const T u = T(1) - h12 * h21;
			if (u > T(0)) {
				flag = T(0);
				*d1 = *d1 / u;
				*d2 = *d2 / u;
				*x1 = *x1 * u;
			}
			else {
				flag = T(-1);
				h11 = T(0);
				h12 = T(0);
				h21 = T(0);
				h22 = T(0);
				*d1 = T(0);
				*d2 = T(0);
				*x1 = T(0);
			}
		}
		else {
			if (q2 < T(0)) {
				flag = T(-1);
				h11 = T(0);
				h12 = T(0);
				h21 = T(0);
				h22 = T(0);
				*d1 = T(0);
				*d2 = T(0);
				*x1 = T(0);
			}
			else {
				flag = T(1);
				h11 = p1 / p2;
				h22 = *x1 / y1;
				const T u = T(1) + h11 * h22;
				const T tmp = *d2 / u;
				*d2 = *d1 / u;
				*d1 = tmp;
				*x1 = y1 * u;
			}
		}
		if (*d1 != T(0)) {
			while (*d1 <= rgamsq || *d1 >= gamsq) {
				if (flag == T(0)) {
					h11 = T(1);
					h22 = T(1);
					flag = T(-1);
				}
				else {
					h21 = T(-1);
					h12 = T(1);
					flag = T(-1);
				}
				if (*d1 <= rgamsq) {
					*d1 = *d1 * (gam * gam);
					*x1 = *x1 / gam;
					h11 = h11 / gam;
					h12 = h12 / gam;
				}
				else {
					*d1 = *d1 / (gam * gam);
					*x1 = *x1 * gam;
					h11 = h11 * gam;
					h12 = h12 * gam;
				}
			}
		}
		if (*d2 != T(0)) {
			while (abs(*d2) <= rgamsq || abs(*d2) >= gamsq) {
				if (flag == T(0)) {
					h11 = T(1);
					h22 = T(1);
					flag = T(-1);
				}
				else {
					h21 = T(-1);
					h12 = T(1);
					flag = T(-1);
				}
				if (abs(*d2) <= rgamsq) {
					*d2 = *d2 * (gam * gam);
					h21 = h21 / gam;
					h22 = h22 / gam;
				}
				else {
					*d2 = *d2 / (gam * gam);
					h21 = h21 * gam;
					h22 = h22 * gam;
				}
			}
		}
	}
	if (flag < T(0)) {
		param[1] = h11;
		param[2] = h21;
		param[3] = h12;
		param[4] = h22;
	}
	else if (flag == T(0)) {
		param[2] = h21;
		param[3] = h12;
	}
	else {
		param[1] = h11;
		param[4] = h22;
	}
	param[0] = flag;
}

#endif // TBLAS_TBLAS_LEVEL1_HPP
