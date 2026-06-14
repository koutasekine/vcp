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
// Modified by Kouta Sekine (2026): translated to C++, rounding-mode control
// added, O(N^3) kernels replaced by rmatmul-based blocked implementations.
// ----------------------------------------------------------------------------

#pragma once

#ifndef VBLAS_RDBLAS_LEVEL1_HPP
#define VBLAS_RDBLAS_LEVEL1_HPP

// Level 1 BLAS (丸めモード指定付き)
// FP 演算を含む関数は r 接頭辞 + 末尾引数 rounding_mode (1 上向き / -1 下向き / 0 最近点)．
// FP 演算を含まない dcopy / dswap / idamax は BLAS と同名 (rounding_mode なし)．
// idamax の返り値は C 流の 0-based index (Fortran BLAS の 1-based とは異なる)．

#include <cmath>
#include <cstddef>

#include "rdblas_common.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// ---- 演算なし (丸めモード引数なし) ----

// y := x
inline void dcopy(const int n, const double* x, const int incx, double* y, const int incy) {
	if (n <= 0) {
		return;
	}
	int ix = vblas_rdblas_detail::vec_start(n, incx);
	int iy = vblas_rdblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		y[iy] = x[ix];
		ix += incx;
		iy += incy;
	}
}

// x <-> y
inline void dswap(const int n, double* x, const int incx, double* y, const int incy) {
	if (n <= 0) {
		return;
	}
	int ix = vblas_rdblas_detail::vec_start(n, incx);
	int iy = vblas_rdblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		const double tmp = x[ix];
		x[ix] = y[iy];
		y[iy] = tmp;
		ix += incx;
		iy += incy;
	}
}

// |x_i| が最大となる 0-based index を返す (n<=0 または incx<=0 のとき -1)．
// |.| と比較のみで丸め誤差は生じない．
inline int idamax(const int n, const double* x, const int incx) {
	if (n <= 0 || incx <= 0) {
		return -1;
	}
	int result = 0;
	double dmax = std::fabs(x[0]);
	int ix = incx;
	for (int i = 1; i < n; i++) {
		const double v = std::fabs(x[ix]);
		if (v > dmax) {
			dmax = v;
			result = i;
		}
		ix += incx;
	}
	return result;
}

// ---- 演算あり (r 接頭辞 + rounding_mode) ----

// x := alpha*x
inline void rdscal(const int n, const double alpha, double* x, const int incx, const int rounding_mode) {
	if (n <= 0 || incx <= 0 || alpha == 1.0) {
		return;
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	int ix = 0;
	for (int i = 0; i < n; i++) {
		x[ix] = alpha * x[ix];
		ix += incx;
	}
}

// y := alpha*x + y (要素ごとに fma 1 回)
inline void rdaxpy(const int n, const double alpha, const double* x, const int incx, double* y, const int incy, const int rounding_mode) {
	if (n <= 0 || alpha == 0.0) {
		return;
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	int ix = vblas_rdblas_detail::vec_start(n, incx);
	int iy = vblas_rdblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		y[iy] = std::fma(alpha, x[ix], y[iy]);
		ix += incx;
		iy += incy;
	}
}

// return x^T y (fma による逐次和)
inline double rddot(const int n, const double* x, const int incx, const double* y, const int incy, const int rounding_mode) {
	if (n <= 0) {
		return 0.0;
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	int ix = vblas_rdblas_detail::vec_start(n, incx);
	int iy = vblas_rdblas_detail::vec_start(n, incy);
	double s = 0.0;
	for (int i = 0; i < n; i++) {
		s = std::fma(x[ix], y[iy], s);
		ix += incx;
		iy += incy;
	}
	return s;
}

// return sum |x_i|
inline double rdasum(const int n, const double* x, const int incx, const int rounding_mode) {
	if (n <= 0 || incx <= 0) {
		return 0.0;
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	double s = 0.0;
	int ix = 0;
	for (int i = 0; i < n; i++) {
		s = s + std::fabs(x[ix]);
		ix += incx;
	}
	return s;
}

// return ||x||_2 (overflow/underflow 回避の scale 付き，reference dnrm2 と同方式)
inline double rdnrm2(const int n, const double* x, const int incx, const int rounding_mode) {
	if (n <= 0 || incx <= 0) {
		return 0.0;
	}
	if (n == 1) {
		return std::fabs(x[0]);
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	double scale = 0.0;
	double ssq = 1.0;
	int ix = 0;
	for (int i = 0; i < n; i++) {
		const double v = std::fabs(x[ix]);
		if (v != 0.0) {
			if (scale < v) {
				const double r = scale / v;
				ssq = std::fma(ssq, r * r, 1.0);
				scale = v;
			}
			else {
				const double r = v / scale;
				ssq = std::fma(r, r, ssq);
			}
		}
		ix += incx;
	}
	return scale * std::sqrt(ssq);
}

// 平面回転の適用: (x_i, y_i) := (c*x_i + s*y_i, c*y_i - s*x_i)
inline void rdrot(const int n, double* x, const int incx, double* y, const int incy, const double c, const double s, const int rounding_mode) {
	if (n <= 0) {
		return;
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	int ix = vblas_rdblas_detail::vec_start(n, incx);
	int iy = vblas_rdblas_detail::vec_start(n, incy);
	for (int i = 0; i < n; i++) {
		const double tmp = std::fma(s, y[iy], c * x[ix]);
		y[iy] = std::fma(c, y[iy], -(s * x[ix]));
		x[ix] = tmp;
		ix += incx;
		iy += incy;
	}
}

// Givens 回転の生成 (reference drotg と同方式)
inline void rdrotg(double* a, double* b, double* c, double* s, const int rounding_mode) {
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	double roe = *b;
	if (std::fabs(*a) > std::fabs(*b)) {
		roe = *a;
	}
	const double scale = std::fabs(*a) + std::fabs(*b);
	double r, z;
	if (scale == 0.0) {
		*c = 1.0;
		*s = 0.0;
		r = 0.0;
		z = 0.0;
	}
	else {
		const double as = *a / scale;
		const double bs = *b / scale;
		r = scale * std::sqrt(std::fma(as, as, bs * bs));
		if (roe < 0.0) {
			r = -r;
		}
		*c = *a / r;
		*s = *b / r;
		z = 1.0;
		if (std::fabs(*a) > std::fabs(*b)) {
			z = *s;
		}
		if (std::fabs(*b) >= std::fabs(*a) && *c != 0.0) {
			z = 1.0 / *c;
		}
	}
	*a = r;
	*b = z;
}

// modified Givens 回転の適用 (reference drotm と同方式)
inline void rdrotm(const int n, double* x, const int incx, double* y, const int incy, const double* param, const int rounding_mode) {
	const double flag = param[0];
	if (n <= 0 || flag == -2.0) {
		return;
	}
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	int ix = vblas_rdblas_detail::vec_start(n, incx);
	int iy = vblas_rdblas_detail::vec_start(n, incy);
	if (flag == 0.0) {
		const double h12 = param[3];
		const double h21 = param[2];
		for (int i = 0; i < n; i++) {
			const double w = x[ix];
			const double z = y[iy];
			x[ix] = std::fma(h12, z, w);
			y[iy] = std::fma(h21, w, z);
			ix += incx;
			iy += incy;
		}
	}
	else if (flag == 1.0) {
		const double h11 = param[1];
		const double h22 = param[4];
		for (int i = 0; i < n; i++) {
			const double w = x[ix];
			const double z = y[iy];
			x[ix] = std::fma(h11, w, z);
			y[iy] = std::fma(h22, z, -w);
			ix += incx;
			iy += incy;
		}
	}
	else {
		const double h11 = param[1];
		const double h12 = param[3];
		const double h21 = param[2];
		const double h22 = param[4];
		for (int i = 0; i < n; i++) {
			const double w = x[ix];
			const double z = y[iy];
			x[ix] = std::fma(h11, w, h12 * z);
			y[iy] = std::fma(h21, w, h22 * z);
			ix += incx;
			iy += incy;
		}
	}
}

// modified Givens 回転の生成 (reference drotmg と同方式)
inline void rdrotmg(double* d1, double* d2, double* x1, const double y1, double* param, const int rounding_mode) {
	vblas_rdblas_detail::RoundingGuard guard(vblas_rdblas_detail::fe_rounding(rounding_mode));
	const double gam = 4096.0;
	const double gamsq = 16777216.0;
	const double rgamsq = 5.9604645e-8;
	double flag, h11 = 0.0, h12 = 0.0, h21 = 0.0, h22 = 0.0;

	if (*d1 < 0.0) {
		flag = -1.0;
		*d1 = 0.0;
		*d2 = 0.0;
		*x1 = 0.0;
	}
	else {
		const double p2 = *d2 * y1;
		if (p2 == 0.0) {
			param[0] = -2.0;
			return;
		}
		const double p1 = *d1 * *x1;
		const double q2 = p2 * y1;
		const double q1 = p1 * *x1;
		if (std::fabs(q1) > std::fabs(q2)) {
			h21 = -y1 / *x1;
			h12 = p2 / p1;
			const double u = 1.0 - h12 * h21;
			if (u > 0.0) {
				flag = 0.0;
				*d1 = *d1 / u;
				*d2 = *d2 / u;
				*x1 = *x1 * u;
			}
			else {
				flag = -1.0;
				h11 = 0.0;
				h12 = 0.0;
				h21 = 0.0;
				h22 = 0.0;
				*d1 = 0.0;
				*d2 = 0.0;
				*x1 = 0.0;
			}
		}
		else {
			if (q2 < 0.0) {
				flag = -1.0;
				h11 = 0.0;
				h12 = 0.0;
				h21 = 0.0;
				h22 = 0.0;
				*d1 = 0.0;
				*d2 = 0.0;
				*x1 = 0.0;
			}
			else {
				flag = 1.0;
				h11 = p1 / p2;
				h22 = *x1 / y1;
				const double u = 1.0 + h11 * h22;
				const double tmp = *d2 / u;
				*d2 = *d1 / u;
				*d1 = tmp;
				*x1 = y1 * u;
			}
		}
		if (*d1 != 0.0) {
			while (*d1 <= rgamsq || *d1 >= gamsq) {
				if (flag == 0.0) {
					h11 = 1.0;
					h22 = 1.0;
					flag = -1.0;
				}
				else {
					h21 = -1.0;
					h12 = 1.0;
					flag = -1.0;
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
		if (*d2 != 0.0) {
			while (std::fabs(*d2) <= rgamsq || std::fabs(*d2) >= gamsq) {
				if (flag == 0.0) {
					h11 = 1.0;
					h22 = 1.0;
					flag = -1.0;
				}
				else {
					h21 = -1.0;
					h12 = 1.0;
					flag = -1.0;
				}
				if (std::fabs(*d2) <= rgamsq) {
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
	if (flag < 0.0) {
		param[1] = h11;
		param[2] = h21;
		param[3] = h12;
		param[4] = h22;
	}
	else if (flag == 0.0) {
		param[2] = h21;
		param[3] = h12;
	}
	else {
		param[1] = h11;
		param[4] = h22;
	}
	param[0] = flag;
}

} // namespace vcp

#endif // VBLAS_RDBLAS_LEVEL1_HPP
