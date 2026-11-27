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

#ifndef TLAPACK_TLAPACK_EIGDC_HPP
#define TLAPACK_TLAPACK_EIGDC_HPP

// 対称固有値問題の分割統治法 (template 版)．
// reference LAPACK 3.12.1 の dlamrg/dlaed5/dlaed6/dlaed4/dlaed2/dlaed3/dlaed1/
// dlaeda/dlaed8/dlaed9/dlaed7/dlaed0/dstedc/dsyevd/dsygvd を移植．
//
// - 元の goto 構造と作業領域内の位置 (IQ, IZ など) を保つため 1-based の
//   局所 index で記述している (配列 access は [idx-1])．
// - index (tlamrg の出力，indxq, perm, givcol など) は 1-based．
// - WORK/IWORK 引数は持たず，各 routine の内部で確保する．
// - ILAENV(9) (分割統治から QR 法へ切り替える大きさ) は 25．
// - tlaed6 のしきい値 SMALL1 = BASE^INT(LOG(SAFMIN)/LOG(BASE)/3) は，対数を使わず，
//   T の乗算と比較だけで求める (eigdc_ilog2 と pow2i．基数は 2 とする)．

#include <vcp/tlapack/tlapack_sy.hpp>

namespace tlapack_detail {

// LAPACK の ILAENV( 9, ... ) の既定値 (reference LAPACK 3.12.1)
const int eigdc_smlsiz = 25;

// 2^e <= x となる最大の整数 e (0 < x < 1)．T の乗算と比較だけで求める
template <typename T>
inline long long eigdc_ilog2(const T& x) {
	long long lo = -1;
	while (!(pow2i<T>(lo) <= x)) {
		lo = lo * 2;
	}
	long long hi = lo / 2;
	while (hi - lo > 1) {
		const long long mid = lo + (hi - lo) / 2;
		if (pow2i<T>(mid) <= x) {
			lo = mid;
		}
		else {
			hi = mid;
		}
	}
	return lo;
}

} // namespace tlapack_detail

namespace vcp {

// 昇順または降順に整列済みの 2 つの列 a(1:n1), a(n1+1:n1+n2) を昇順に併合する置換を作る．
// dtrd1, dtrd2: 1 なら昇順，-1 なら降順に格納されている．index は 1-based
template <typename T>
inline void tlamrg(const int n1, const int n2, const T* a, const int dtrd1, const int dtrd2, int* index) {
	int n1sv = n1;
	int n2sv = n2;
	int ind1, ind2;
	if (dtrd1 > 0) {
		ind1 = 1;
	}
	else {
		ind1 = n1;
	}
	if (dtrd2 > 0) {
		ind2 = 1 + n1;
	}
	else {
		ind2 = n1 + n2;
	}
	int i = 1;
L10:
	if (n1sv > 0 && n2sv > 0) {
		if (a[ind1 - 1] <= a[ind2 - 1]) {
			index[i - 1] = ind1;
			i = i + 1;
			ind1 = ind1 + dtrd1;
			n1sv = n1sv - 1;
		}
		else {
			index[i - 1] = ind2;
			i = i + 1;
			ind2 = ind2 + dtrd2;
			n2sv = n2sv - 1;
		}
		goto L10;
	}
	if (n1sv == 0) {
		for (n1sv = 1; n1sv <= n2sv; n1sv++) {
			index[i - 1] = ind2;
			i = i + 1;
			ind2 = ind2 + dtrd2;
		}
	}
	else {
		for (n2sv = 1; n2sv <= n1sv; n2sv++) {
			index[i - 1] = ind1;
			i = i + 1;
			ind1 = ind1 + dtrd1;
		}
	}
}

// 2x2 の secular 方程式の i 番目の根 (diag(d) + rho*z*z^T, d(1) < d(2), rho > 0)
template <typename T>
inline void tlaed5(const int i, const T* d, const T* z, T* delta,
	const T rho, T& dlam) {
	namespace detail = tlapack_detail;
	T b, c, tau, temp, w;
	const T del = d[1] - d[0];
	if (i == 1) {
		w = T(1) + T(2) * rho * (z[1] * z[1] - z[0] * z[0]) / del;
		if (w > T(0)) {
			b = del + rho * (z[0] * z[0] + z[1] * z[1]);
			c = rho * z[0] * z[0] * del;
			tau = T(2) * c / (b + tlapack_detail::tsqrt(tlapack_detail::tabs(b * b - T(4) * c)));
			dlam = d[0] + tau;
			delta[0] = -(z[0] / tau);
			delta[1] = z[1] / (del - tau);
		}
		else {
			b = -del + rho * (z[0] * z[0] + z[1] * z[1]);
			c = rho * z[1] * z[1] * del;
			if (b > T(0)) {
				tau = -(T(2) * c / (b + tlapack_detail::tsqrt(b * b + T(4) * c)));
			}
			else {
				tau = (b - tlapack_detail::tsqrt(b * b + T(4) * c)) / T(2);
			}
			dlam = d[1] + tau;
			delta[0] = -(z[0] / (del + tau));
			delta[1] = -(z[1] / tau);
		}
		temp = tlapack_detail::tsqrt(delta[0] * delta[0] + delta[1] * delta[1]);
		delta[0] = delta[0] / temp;
		delta[1] = delta[1] / temp;
	}
	else {
		b = -del + rho * (z[0] * z[0] + z[1] * z[1]);
		c = rho * z[1] * z[1] * del;
		if (b > T(0)) {
			tau = (b + tlapack_detail::tsqrt(b * b + T(4) * c)) / T(2);
		}
		else {
			tau = T(2) * c / (-b + tlapack_detail::tsqrt(b * b + T(4) * c));
		}
		dlam = d[1] + tau;
		delta[0] = -(z[0] / (del + tau));
		delta[1] = -(z[1] / tau);
		temp = tlapack_detail::tsqrt(delta[0] * delta[0] + delta[1] * delta[1]);
		delta[0] = delta[0] / temp;
		delta[1] = delta[1] / temp;
	}
}

// tlaed4 の中で使う 3 項の secular 方程式
//   f(x) = rho + z(1)/(d(1)-x) + z(2)/(d(2)-x) + z(3)/(d(3)-x)
// の根．返り値 1 は収束しなかったことを表す
template <typename T>
inline int tlaed6(const int kniter, const bool orgati, const T rho, const T* d,
	const T* z, const T finit, T& tau) {
	namespace detail = tlapack_detail;
	const int maxit = 40;
	T dscale[3], zscale[3];
	bool scale;
	int i, iter, niter;
	T a, b, c, ddf, df, eps, erretm, eta, f, fc, sclfac, sclinv, small1, small2, sminv1,
		sminv2, temp, temp1, temp2, temp3, temp4, lbd, ubd;
	int info = 0;
	sclfac = T(0);
	sclinv = T(0);

	if (orgati) {
		lbd = d[1];
		ubd = d[2];
	}
	else {
		lbd = d[0];
		ubd = d[1];
	}
	if (finit < T(0)) {
		lbd = T(0);
	}
	else {
		ubd = T(0);
	}

	niter = 1;
	tau = T(0);
	if (kniter == 2) {
		if (orgati) {
			temp = (d[2] - d[1]) / T(2);
			c = rho + z[0] / ((d[0] - d[1]) - temp);
			a = c * (d[1] + d[2]) + z[1] + z[2];
			b = c * d[1] * d[2] + z[1] * d[2] + z[2] * d[1];
		}
		else {
			temp = (d[0] - d[1]) / T(2);
			c = rho + z[2] / ((d[2] - d[1]) - temp);
			a = c * (d[0] + d[1]) + z[0] + z[1];
			b = c * d[0] * d[1] + z[0] * d[1] + z[1] * d[0];
		}
		temp = std::max(std::max(tlapack_detail::tabs(a), tlapack_detail::tabs(b)), tlapack_detail::tabs(c));
		a = a / temp;
		b = b / temp;
		c = c / temp;
		if (c == T(0)) {
			tau = b / a;
		}
		else if (a <= T(0)) {
			tau = (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
		}
		else {
			tau = T(2) * b / (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
		}
		if (tau < lbd || tau > ubd) {
			tau = (lbd + ubd) / T(2);
		}
		if (d[0] == tau || d[1] == tau || d[2] == tau) {
			tau = T(0);
		}
		else {
			temp = finit + tau * z[0] / (d[0] * (d[0] - tau)) +
				tau * z[1] / (d[1] * (d[1] - tau)) +
				tau * z[2] / (d[2] * (d[2] - tau));
			if (temp <= T(0)) {
				lbd = tau;
			}
			else {
				ubd = tau;
			}
			if (tlapack_detail::tabs(finit) <= tlapack_detail::tabs(temp)) {
				tau = T(0);
			}
		}
	}

	// 値が小さすぎるときに拡大するための，しきい値と倍率
	eps = tlamch<T>('E');
	small1 = detail::pow2i<T>(detail::eigdc_ilog2(tlamch<T>('S')) / 3);
	sminv1 = T(1) / small1;
	small2 = small1 * small1;
	sminv2 = sminv1 * sminv1;

	if (orgati) {
		temp = std::min(tlapack_detail::tabs(d[1] - tau), tlapack_detail::tabs(d[2] - tau));
	}
	else {
		temp = std::min(tlapack_detail::tabs(d[0] - tau), tlapack_detail::tabs(d[1] - tau));
	}
	scale = false;
	if (temp <= small1) {
		scale = true;
		if (temp <= small2) {
			sclfac = sminv2;
			sclinv = small2;
		}
		else {
			sclfac = sminv1;
			sclinv = small1;
		}
		for (i = 1; i <= 3; i++) {
			dscale[i - 1] = d[i - 1] * sclfac;
			zscale[i - 1] = z[i - 1] * sclfac;
		}
		tau = tau * sclfac;
		lbd = lbd * sclfac;
		ubd = ubd * sclfac;
	}
	else {
		for (i = 1; i <= 3; i++) {
			dscale[i - 1] = d[i - 1];
			zscale[i - 1] = z[i - 1];
		}
	}

	fc = T(0);
	df = T(0);
	ddf = T(0);
	for (i = 1; i <= 3; i++) {
		temp = T(1) / (dscale[i - 1] - tau);
		temp1 = zscale[i - 1] * temp;
		temp2 = temp1 * temp;
		temp3 = temp2 * temp;
		fc = fc + temp1 / dscale[i - 1];
		df = df + temp2;
		ddf = ddf + temp3;
	}
	f = finit + tau * fc;

	if (tlapack_detail::tabs(f) <= T(0)) {
		goto L60;
	}
	if (f <= T(0)) {
		lbd = tau;
	}
	else {
		ubd = tau;
	}

	iter = niter + 1;

	for (niter = iter; niter <= maxit; niter++) {
		if (orgati) {
			temp1 = dscale[1] - tau;
			temp2 = dscale[2] - tau;
		}
		else {
			temp1 = dscale[0] - tau;
			temp2 = dscale[1] - tau;
		}
		a = (temp1 + temp2) * f - temp1 * temp2 * df;
		b = temp1 * temp2 * f;
		c = f - (temp1 + temp2) * df + temp1 * temp2 * ddf;
		temp = std::max(std::max(tlapack_detail::tabs(a), tlapack_detail::tabs(b)), tlapack_detail::tabs(c));
		a = a / temp;
		b = b / temp;
		c = c / temp;
		if (c == T(0)) {
			eta = b / a;
		}
		else if (a <= T(0)) {
			eta = (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
		}
		else {
			eta = T(2) * b / (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
		}
		if (f * eta >= T(0)) {
			eta = -(f / df);
		}

		tau = tau + eta;
		if (tau < lbd || tau > ubd) {
			tau = (lbd + ubd) / T(2);
		}

		fc = T(0);
		erretm = T(0);
		df = T(0);
		ddf = T(0);
		for (i = 1; i <= 3; i++) {
			if ((dscale[i - 1] - tau) != T(0)) {
				temp = T(1) / (dscale[i - 1] - tau);
				temp1 = zscale[i - 1] * temp;
				temp2 = temp1 * temp;
				temp3 = temp2 * temp;
				temp4 = temp1 / dscale[i - 1];
				fc = fc + temp4;
				erretm = erretm + tlapack_detail::tabs(temp4);
				df = df + temp2;
				ddf = ddf + temp3;
			}
			else {
				goto L60;
			}
		}
		f = finit + tau * fc;
		erretm = T(8) * (tlapack_detail::tabs(finit) + tlapack_detail::tabs(tau) * erretm) +
			tlapack_detail::tabs(tau) * df;
		if ((tlapack_detail::tabs(f) <= T(4) * eps * erretm) ||
			((ubd - lbd) <= T(4) * eps * tlapack_detail::tabs(tau))) {
			goto L60;
		}
		if (f <= T(0)) {
			lbd = tau;
		}
		else {
			ubd = tau;
		}
	}
	info = 1;
L60:

	// 拡大していたら戻す
	if (scale) {
		tau = tau * sclinv;
	}
	return info;
}

// secular 方程式の i 番目の根 (diag(d) + rho*z*z^T の i 番目の固有値, d は昇順, rho > 0)．
// delta(j) には d(j) - dlam が入る．返り値 1 は収束しなかったことを表す
template <typename T>
inline int tlaed4(const int n, const int i, const T* d, const T* z, T* delta,
	const T rho, T& dlam) {
	namespace detail = tlapack_detail;
	const int maxit = 30;
	bool orgati, swtch, swtch3;
	int ii, iim1, iip1, ip1, iter, j, niter;
	T a, b, c, del, dltlb, dltub, dphi, dpsi, dw, eps, erretm, eta, midpt, phi, prew, psi,
		rhoinv, tau, temp, temp1, w;
	T zz[3];
	int info = 0;

	if (n == 1) {
		// n = 1 のとき
		dlam = d[0] + rho * z[0] * z[0];
		delta[0] = T(1);
		return info;
	}
	if (n == 2) {
		tlaed5(i, d, z, delta, rho, dlam);
		return info;
	}

	eps = tlamch<T>('E');
	rhoinv = T(1) / rho;

	if (i == n) {
		// 最後の根 (i = n)
		ii = n - 1;
		niter = 1;

		midpt = rho / T(2);
		for (j = 1; j <= n; j++) {
			delta[j - 1] = (d[j - 1] - d[i - 1]) - midpt;
		}

		psi = T(0);
		for (j = 1; j <= n - 2; j++) {
			psi = psi + z[j - 1] * z[j - 1] / delta[j - 1];
		}

		c = rhoinv + psi;
		w = c + z[ii - 1] * z[ii - 1] / delta[ii - 1] +
			z[n - 1] * z[n - 1] / delta[n - 1];

		if (w <= T(0)) {
			temp = z[n - 2] * z[n - 2] / (d[n - 1] - d[n - 2] + rho) +
				z[n - 1] * z[n - 1] / rho;
			if (c <= temp) {
				tau = rho;
			}
			else {
				del = d[n - 1] - d[n - 2];
				a = -(c * del) + z[n - 2] * z[n - 2] + z[n - 1] * z[n - 1];
				b = z[n - 1] * z[n - 1] * del;
				if (a < T(0)) {
					tau = T(2) * b / (tlapack_detail::tsqrt(a * a + T(4) * b * c) - a);
				}
				else {
					tau = (a + tlapack_detail::tsqrt(a * a + T(4) * b * c)) / (T(2) * c);
				}
			}
			dltlb = midpt;
			dltub = rho;
		}
		else {
			del = d[n - 1] - d[n - 2];
			a = -(c * del) + z[n - 2] * z[n - 2] + z[n - 1] * z[n - 1];
			b = z[n - 1] * z[n - 1] * del;
			if (a < T(0)) {
				tau = T(2) * b / (tlapack_detail::tsqrt(a * a + T(4) * b * c) - a);
			}
			else {
				tau = (a + tlapack_detail::tsqrt(a * a + T(4) * b * c)) / (T(2) * c);
			}
			dltlb = T(0);
			dltub = midpt;
		}

		for (j = 1; j <= n; j++) {
			delta[j - 1] = (d[j - 1] - d[i - 1]) - tau;
		}

		// secular 方程式の値 w と誤差の見積もり
		dpsi = T(0);
		psi = T(0);
		erretm = T(0);
		for (j = 1; j <= ii; j++) {
			temp = z[j - 1] / delta[j - 1];
			psi = psi + z[j - 1] * temp;
			dpsi = dpsi + temp * temp;
			erretm = erretm + psi;
		}
		erretm = tlapack_detail::tabs(erretm);

		temp = z[n - 1] / delta[n - 1];
		phi = z[n - 1] * temp;
		dphi = temp * temp;
		erretm = T(8) * (-phi - psi) + erretm - phi + rhoinv +
			tlapack_detail::tabs(tau) * (dpsi + dphi);

		w = rhoinv + phi + psi;

		if (tlapack_detail::tabs(w) <= eps * erretm) {
			dlam = d[i - 1] + tau;
			goto L250;
		}

		if (w <= T(0)) {
			dltlb = std::max(dltlb, tau);
		}
		else {
			dltub = std::min(dltub, tau);
		}

		niter = niter + 1;
		c = w - delta[n - 2] * dpsi - delta[n - 1] * dphi;
		a = (delta[n - 2] + delta[n - 1]) * w -
			delta[n - 2] * delta[n - 1] * (dpsi + dphi);
		b = delta[n - 2] * delta[n - 1] * w;
		if (c < T(0)) {
			c = tlapack_detail::tabs(c);
		}
		if (c == T(0)) {
			eta = -(w / (dpsi + dphi));
		}
		else if (a >= T(0)) {
			eta = (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
		}
		else {
			eta = T(2) * b / (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
		}

		if (w * eta > T(0)) {
			eta = -(w / (dpsi + dphi));
		}
		temp = tau + eta;
		if (temp > dltub || temp < dltlb) {
			if (w < T(0)) {
				eta = (dltub - tau) / T(2);
			}
			else {
				eta = (dltlb - tau) / T(2);
			}
		}
		for (j = 1; j <= n; j++) {
			delta[j - 1] = delta[j - 1] - eta;
		}

		tau = tau + eta;

		dpsi = T(0);
		psi = T(0);
		erretm = T(0);
		for (j = 1; j <= ii; j++) {
			temp = z[j - 1] / delta[j - 1];
			psi = psi + z[j - 1] * temp;
			dpsi = dpsi + temp * temp;
			erretm = erretm + psi;
		}
		erretm = tlapack_detail::tabs(erretm);

		temp = z[n - 1] / delta[n - 1];
		phi = z[n - 1] * temp;
		dphi = temp * temp;
		erretm = T(8) * (-phi - psi) + erretm - phi + rhoinv +
			tlapack_detail::tabs(tau) * (dpsi + dphi);

		w = rhoinv + phi + psi;

		iter = niter + 1;

		for (niter = iter; niter <= maxit; niter++) {
			if (tlapack_detail::tabs(w) <= eps * erretm) {
				dlam = d[i - 1] + tau;
				goto L250;
			}

			if (w <= T(0)) {
				dltlb = std::max(dltlb, tau);
			}
			else {
				dltub = std::min(dltub, tau);
			}

			c = w - delta[n - 2] * dpsi - delta[n - 1] * dphi;
			a = (delta[n - 2] + delta[n - 1]) * w -
				delta[n - 2] * delta[n - 1] * (dpsi + dphi);
			b = delta[n - 2] * delta[n - 1] * w;
			if (a >= T(0)) {
				eta = (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
			}
			else {
				eta = T(2) * b / (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
			}

			if (w * eta > T(0)) {
				eta = -(w / (dpsi + dphi));
			}
			temp = tau + eta;
			if (temp > dltub || temp < dltlb) {
				if (w < T(0)) {
					eta = (dltub - tau) / T(2);
				}
				else {
					eta = (dltlb - tau) / T(2);
				}
			}
			for (j = 1; j <= n; j++) {
				delta[j - 1] = delta[j - 1] - eta;
			}

			tau = tau + eta;

			dpsi = T(0);
			psi = T(0);
			erretm = T(0);
			for (j = 1; j <= ii; j++) {
				temp = z[j - 1] / delta[j - 1];
				psi = psi + z[j - 1] * temp;
				dpsi = dpsi + temp * temp;
				erretm = erretm + psi;
			}
			erretm = tlapack_detail::tabs(erretm);

			temp = z[n - 1] / delta[n - 1];
			phi = z[n - 1] * temp;
			dphi = temp * temp;
			erretm = T(8) * (-phi - psi) + erretm - phi + rhoinv +
				tlapack_detail::tabs(tau) * (dpsi + dphi);

			w = rhoinv + phi + psi;
		}

		// 収束しなかった
		info = 1;
		dlam = d[i - 1] + tau;
		goto L250;
	}
	else {
		// i < n のとき
		niter = 1;
		ip1 = i + 1;

		del = d[ip1 - 1] - d[i - 1];
		midpt = del / T(2);
		for (j = 1; j <= n; j++) {
			delta[j - 1] = (d[j - 1] - d[i - 1]) - midpt;
		}

		psi = T(0);
		for (j = 1; j <= i - 1; j++) {
			psi = psi + z[j - 1] * z[j - 1] / delta[j - 1];
		}

		phi = T(0);
		for (j = n; j >= i + 2; j--) {
			phi = phi + z[j - 1] * z[j - 1] / delta[j - 1];
		}
		c = rhoinv + psi + phi;
		w = c + z[i - 1] * z[i - 1] / delta[i - 1] +
			z[ip1 - 1] * z[ip1 - 1] / delta[ip1 - 1];

		if (w > T(0)) {
			// d(i) < 根 < (d(i)+d(i+1))/2．原点を d(i) に取る
			orgati = true;
			a = c * del + z[i - 1] * z[i - 1] + z[ip1 - 1] * z[ip1 - 1];
			b = z[i - 1] * z[i - 1] * del;
			if (a > T(0)) {
				tau = T(2) * b / (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
			}
			else {
				tau = (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
			}
			dltlb = T(0);
			dltub = midpt;
		}
		else {
			// (d(i)+d(i+1))/2 <= 根 < d(i+1)．原点を d(i+1) に取る
			orgati = false;
			a = c * del - z[i - 1] * z[i - 1] - z[ip1 - 1] * z[ip1 - 1];
			b = z[ip1 - 1] * z[ip1 - 1] * del;
			if (a < T(0)) {
				tau = T(2) * b / (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a + T(4) * b * c)));
			}
			else {
				tau = -((a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a + T(4) * b * c))) / (T(2) * c));
			}
			dltlb = -midpt;
			dltub = T(0);
		}

		if (orgati) {
			for (j = 1; j <= n; j++) {
				delta[j - 1] = (d[j - 1] - d[i - 1]) - tau;
			}
		}
		else {
			for (j = 1; j <= n; j++) {
				delta[j - 1] = (d[j - 1] - d[ip1 - 1]) - tau;
			}
		}
		if (orgati) {
			ii = i;
		}
		else {
			ii = i + 1;
		}
		iim1 = ii - 1;
		iip1 = ii + 1;

		// secular 方程式の値 w と誤差の見積もり
		dpsi = T(0);
		psi = T(0);
		erretm = T(0);
		for (j = 1; j <= iim1; j++) {
			temp = z[j - 1] / delta[j - 1];
			psi = psi + z[j - 1] * temp;
			dpsi = dpsi + temp * temp;
			erretm = erretm + psi;
		}
		erretm = tlapack_detail::tabs(erretm);

		dphi = T(0);
		phi = T(0);
		for (j = n; j >= iip1; j--) {
			temp = z[j - 1] / delta[j - 1];
			phi = phi + z[j - 1] * temp;
			dphi = dphi + temp * temp;
			erretm = erretm + phi;
		}

		w = rhoinv + phi + psi;

		// w の符号が見込みと逆なら，3 項の補間 (tlaed6) に切り替える
		swtch3 = false;
		if (orgati) {
			if (w < T(0)) {
				swtch3 = true;
			}
		}
		else {
			if (w > T(0)) {
				swtch3 = true;
			}
		}
		if (ii == 1 || ii == n) {
			swtch3 = false;
		}

		temp = z[ii - 1] / delta[ii - 1];
		dw = dpsi + dphi + temp * temp;
		temp = z[ii - 1] * temp;
		w = w + temp;
		erretm = T(8) * (phi - psi) + erretm + T(2) * rhoinv +
			T(3) * tlapack_detail::tabs(temp) + tlapack_detail::tabs(tau) * dw;

		if (tlapack_detail::tabs(w) <= eps * erretm) {
			if (orgati) {
				dlam = d[i - 1] + tau;
			}
			else {
				dlam = d[ip1 - 1] + tau;
			}
			goto L250;
		}

		if (w <= T(0)) {
			dltlb = std::max(dltlb, tau);
		}
		else {
			dltub = std::min(dltub, tau);
		}

		niter = niter + 1;
		if (!swtch3) {
			if (orgati) {
				temp1 = z[i - 1] / delta[i - 1];
				c = w - delta[ip1 - 1] * dw - (d[i - 1] - d[ip1 - 1]) *
					(temp1 * temp1);
			}
			else {
				temp1 = z[ip1 - 1] / delta[ip1 - 1];
				c = w - delta[i - 1] * dw - (d[ip1 - 1] - d[i - 1]) *
					(temp1 * temp1);
			}
			a = (delta[i - 1] + delta[ip1 - 1]) * w -
				delta[i - 1] * delta[ip1 - 1] * dw;
			b = delta[i - 1] * delta[ip1 - 1] * w;
			if (c == T(0)) {
				if (a == T(0)) {
					if (orgati) {
						a = z[i - 1] * z[i - 1] + delta[ip1 - 1] * delta[ip1 - 1] *
							(dpsi + dphi);
					}
					else {
						a = z[ip1 - 1] * z[ip1 - 1] + delta[i - 1] * delta[i - 1] *
							(dpsi + dphi);
					}
				}
				eta = b / a;
			}
			else if (a <= T(0)) {
				eta = (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
			}
			else {
				eta = T(2) * b / (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
			}
		}
		else {
			// 3 項の補間
			temp = rhoinv + psi + phi;
			if (orgati) {
				temp1 = z[iim1 - 1] / delta[iim1 - 1];
				temp1 = temp1 * temp1;
				c = temp - delta[iip1 - 1] * (dpsi + dphi) -
					(d[iim1 - 1] - d[iip1 - 1]) * temp1;
				zz[0] = z[iim1 - 1] * z[iim1 - 1];
				zz[2] = delta[iip1 - 1] * delta[iip1 - 1] *
					((dpsi - temp1) + dphi);
			}
			else {
				temp1 = z[iip1 - 1] / delta[iip1 - 1];
				temp1 = temp1 * temp1;
				c = temp - delta[iim1 - 1] * (dpsi + dphi) -
					(d[iip1 - 1] - d[iim1 - 1]) * temp1;
				zz[0] = delta[iim1 - 1] * delta[iim1 - 1] *
					(dpsi + (dphi - temp1));
				zz[2] = z[iip1 - 1] * z[iip1 - 1];
			}
			zz[1] = z[ii - 1] * z[ii - 1];
			info = tlaed6(niter, orgati, c, delta + (iim1 - 1), zz, w, eta);
			if (info != 0) {
				goto L250;
			}
		}

		if (w * eta >= T(0)) {
			eta = -(w / dw);
		}
		temp = tau + eta;
		if (temp > dltub || temp < dltlb) {
			if (w < T(0)) {
				eta = (dltub - tau) / T(2);
			}
			else {
				eta = (dltlb - tau) / T(2);
			}
		}

		prew = w;

		for (j = 1; j <= n; j++) {
			delta[j - 1] = delta[j - 1] - eta;
		}

		dpsi = T(0);
		psi = T(0);
		erretm = T(0);
		for (j = 1; j <= iim1; j++) {
			temp = z[j - 1] / delta[j - 1];
			psi = psi + z[j - 1] * temp;
			dpsi = dpsi + temp * temp;
			erretm = erretm + psi;
		}
		erretm = tlapack_detail::tabs(erretm);

		dphi = T(0);
		phi = T(0);
		for (j = n; j >= iip1; j--) {
			temp = z[j - 1] / delta[j - 1];
			phi = phi + z[j - 1] * temp;
			dphi = dphi + temp * temp;
			erretm = erretm + phi;
		}

		temp = z[ii - 1] / delta[ii - 1];
		dw = dpsi + dphi + temp * temp;
		temp = z[ii - 1] * temp;
		w = rhoinv + phi + psi + temp;
		erretm = T(8) * (phi - psi) + erretm + T(2) * rhoinv +
			T(3) * tlapack_detail::tabs(temp) + tlapack_detail::tabs(tau + eta) * dw;

		swtch = false;
		if (orgati) {
			if (-w > tlapack_detail::tabs(prew) / T(10)) {
				swtch = true;
			}
		}
		else {
			if (w > tlapack_detail::tabs(prew) / T(10)) {
				swtch = true;
			}
		}

		tau = tau + eta;

		iter = niter + 1;

		for (niter = iter; niter <= maxit; niter++) {
			if (tlapack_detail::tabs(w) <= eps * erretm) {
				if (orgati) {
					dlam = d[i - 1] + tau;
				}
				else {
					dlam = d[ip1 - 1] + tau;
				}
				goto L250;
			}

			if (w <= T(0)) {
				dltlb = std::max(dltlb, tau);
			}
			else {
				dltub = std::min(dltub, tau);
			}

			if (!swtch3) {
				if (!swtch) {
					if (orgati) {
						temp1 = z[i - 1] / delta[i - 1];
						c = w - delta[ip1 - 1] * dw -
							(d[i - 1] - d[ip1 - 1]) * (temp1 * temp1);
					}
					else {
						temp1 = z[ip1 - 1] / delta[ip1 - 1];
						c = w - delta[i - 1] * dw - (d[ip1 - 1] - d[i - 1]) *
							(temp1 * temp1);
					}
				}
				else {
					temp = z[ii - 1] / delta[ii - 1];
					if (orgati) {
						dpsi = dpsi + temp * temp;
					}
					else {
						dphi = dphi + temp * temp;
					}
					c = w - delta[i - 1] * dpsi - delta[ip1 - 1] * dphi;
				}
				a = (delta[i - 1] + delta[ip1 - 1]) * w -
					delta[i - 1] * delta[ip1 - 1] * dw;
				b = delta[i - 1] * delta[ip1 - 1] * w;
				if (c == T(0)) {
					if (a == T(0)) {
						if (!swtch) {
							if (orgati) {
								a = z[i - 1] * z[i - 1] + delta[ip1 - 1] *
									delta[ip1 - 1] * (dpsi + dphi);
							}
							else {
								a = z[ip1 - 1] * z[ip1 - 1] +
									delta[i - 1] * delta[i - 1] * (dpsi + dphi);
							}
						}
						else {
							a = delta[i - 1] * delta[i - 1] * dpsi +
								delta[ip1 - 1] * delta[ip1 - 1] * dphi;
						}
					}
					eta = b / a;
				}
				else if (a <= T(0)) {
					eta = (a - tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c))) / (T(2) * c);
				}
				else {
					eta = T(2) * b / (a + tlapack_detail::tsqrt(tlapack_detail::tabs(a * a - T(4) * b * c)));
				}
			}
			else {
				// 3 項の補間
				temp = rhoinv + psi + phi;
				if (swtch) {
					c = temp - delta[iim1 - 1] * dpsi - delta[iip1 - 1] * dphi;
					zz[0] = delta[iim1 - 1] * delta[iim1 - 1] * dpsi;
					zz[2] = delta[iip1 - 1] * delta[iip1 - 1] * dphi;
				}
				else {
					if (orgati) {
						temp1 = z[iim1 - 1] / delta[iim1 - 1];
						temp1 = temp1 * temp1;
						c = temp - delta[iip1 - 1] * (dpsi + dphi) -
							(d[iim1 - 1] - d[iip1 - 1]) * temp1;
						zz[0] = z[iim1 - 1] * z[iim1 - 1];
						zz[2] = delta[iip1 - 1] * delta[iip1 - 1] *
							((dpsi - temp1) + dphi);
					}
					else {
						temp1 = z[iip1 - 1] / delta[iip1 - 1];
						temp1 = temp1 * temp1;
						c = temp - delta[iim1 - 1] * (dpsi + dphi) -
							(d[iip1 - 1] - d[iim1 - 1]) * temp1;
						zz[0] = delta[iim1 - 1] * delta[iim1 - 1] *
							(dpsi + (dphi - temp1));
						zz[2] = z[iip1 - 1] * z[iip1 - 1];
					}
				}
				info = tlaed6(niter, orgati, c, delta + (iim1 - 1), zz, w, eta);
				if (info != 0) {
					goto L250;
				}
			}

			if (w * eta >= T(0)) {
				eta = -(w / dw);
			}
			temp = tau + eta;
			if (temp > dltub || temp < dltlb) {
				if (w < T(0)) {
					eta = (dltub - tau) / T(2);
				}
				else {
					eta = (dltlb - tau) / T(2);
				}
			}

			for (j = 1; j <= n; j++) {
				delta[j - 1] = delta[j - 1] - eta;
			}

			tau = tau + eta;
			prew = w;

			dpsi = T(0);
			psi = T(0);
			erretm = T(0);
			for (j = 1; j <= iim1; j++) {
				temp = z[j - 1] / delta[j - 1];
				psi = psi + z[j - 1] * temp;
				dpsi = dpsi + temp * temp;
				erretm = erretm + psi;
			}
			erretm = tlapack_detail::tabs(erretm);

			dphi = T(0);
			phi = T(0);
			for (j = n; j >= iip1; j--) {
				temp = z[j - 1] / delta[j - 1];
				phi = phi + z[j - 1] * temp;
				dphi = dphi + temp * temp;
				erretm = erretm + phi;
			}

			temp = z[ii - 1] / delta[ii - 1];
			dw = dpsi + dphi + temp * temp;
			temp = z[ii - 1] * temp;
			w = rhoinv + phi + psi + temp;
			erretm = T(8) * (phi - psi) + erretm + T(2) * rhoinv +
				T(3) * tlapack_detail::tabs(temp) + tlapack_detail::tabs(tau) * dw;
			if (w * prew > T(0) && tlapack_detail::tabs(w) > tlapack_detail::tabs(prew) / T(10)) {
				swtch = !swtch;
			}
		}

		// 収束しなかった
		info = 1;
		if (orgati) {
			dlam = d[i - 1] + tau;
		}
		else {
			dlam = d[ip1 - 1] + tau;
		}
	}

L250:
	return info;
}

// 併合の前処理 (deflation)．2 つの部分問題の固有値を整列し，z の成分が小さいもの，
// 固有値が近いものを取り除いて secular 方程式の大きさ k を決める (compz = I の経路)．
// 終了時の coltyp(1:4) は 4 種類の列の本数
template <typename T>
inline int tlaed2(int& k, const int n, const int n1, T* d, T* Q, const int ldq,
	int* indxq, T& rho, T* z, T* dlambda, T* w, T* q2,
	int* indx, int* indxc, int* indxp, int* coltyp) {
	namespace detail = tlapack_detail;
	if (n < 0) {
		detail::tlapack_error("tlaed2: invalid n");
	}
	else if (ldq < std::max(1, n)) {
		detail::tlapack_error("tlaed2: invalid ldq");
	}
	else if (std::min(1, (n / 2)) > n1 || (n / 2) < n1) {
		detail::tlapack_error("tlaed2: invalid n1");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t LQ = static_cast<std::size_t>(ldq);
	int ctot[4], psm[4];
	int ct, i, imax, iq1, iq2, j, jmax, js, k2, n1p1, n2, nj, pj;
	T c, eps, s, t, tau, tol;
	pj = 0;

	n2 = n - n1;
	n1p1 = n1 + 1;

	if (rho < T(0)) {
		tscal(n2, -T(1), z + (n1p1 - 1), 1);
	}

	// z を正規化する (2 つの部分問題の z はそれぞれノルム 1 なので，全体は sqrt(2))
	t = T(1) / tlapack_detail::tsqrt(T(2));
	tscal(n, t, z, 1);

	// rho = |2 * rho|
	rho = tlapack_detail::tabs(T(2) * rho);

	// 2 つの部分問題の固有値を昇順に並べる
	for (i = n1p1; i <= n; i++) {
		indxq[i - 1] = indxq[i - 1] + n1;
	}

	for (i = 1; i <= n; i++) {
		dlambda[i - 1] = d[indxq[i - 1] - 1];
	}
	tlamrg(n1, n2, dlambda, 1, 1, indxc);
	for (i = 1; i <= n; i++) {
		indx[i - 1] = indxq[indxc[i - 1] - 1];
	}

	// deflation の許容値
	imax = itamax(n, z, 1) + 1;
	jmax = itamax(n, d, 1) + 1;
	eps = tlamch<T>('E');
	tol = T(8) * eps * std::max(tlapack_detail::tabs(d[jmax - 1]), tlapack_detail::tabs(z[imax - 1]));

	// rho が小さければ，並べ替えるだけで終わる
	if (rho * tlapack_detail::tabs(z[imax - 1]) <= tol) {
		k = 0;
		iq2 = 1;
		for (j = 1; j <= n; j++) {
			i = indx[j - 1];
			tcopy(n, Q + LQ * (i - 1), 1, q2 + (iq2 - 1), 1);
			dlambda[j - 1] = d[i - 1];
			iq2 = iq2 + n;
		}
		tlacpy('A', n, n, q2, n, Q, ldq);
		tcopy(n, dlambda, 1, d, 1);
		goto L190;
	}

	// 列の種類: 1 上半分だけ非零, 2 密, 3 下半分だけ非零, 4 deflation 済み
	for (i = 1; i <= n1; i++) {
		coltyp[i - 1] = 1;
	}
	for (i = n1p1; i <= n; i++) {
		coltyp[i - 1] = 3;
	}

	k = 0;
	k2 = n + 1;
	for (j = 1; j <= n; j++) {
		nj = indx[j - 1];
		if (rho * tlapack_detail::tabs(z[nj - 1]) <= tol) {
			// z の成分が小さいので deflation する
			k2 = k2 - 1;
			coltyp[nj - 1] = 4;
			indxp[k2 - 1] = nj;
			if (j == n) {
				goto L100;
			}
		}
		else {
			pj = nj;
			goto L80;
		}
	}
L80:
	j = j + 1;
	if (j > n) {
		goto L100;
	}
	nj = indx[j - 1];
	if (rho * tlapack_detail::tabs(z[nj - 1]) <= tol) {
		// z の成分が小さいので deflation する
		k2 = k2 - 1;
		coltyp[nj - 1] = 4;
		indxp[k2 - 1] = nj;
	}
	else {
		// 固有値が近ければ，Givens 回転で z の成分を 1 つ消して deflation する
		s = z[pj - 1];
		c = z[nj - 1];
		tau = tlapy2(c, s);
		t = d[nj - 1] - d[pj - 1];
		c = c / tau;
		s = -(s / tau);
		if (tlapack_detail::tabs(t * c * s) <= tol) {
			z[nj - 1] = tau;
			z[pj - 1] = T(0);
			if (coltyp[nj - 1] != coltyp[pj - 1]) {
				coltyp[nj - 1] = 2;
			}
			coltyp[pj - 1] = 4;
			trot(n, Q + LQ * (pj - 1), 1, Q + LQ * (nj - 1), 1, c, s);
			t = d[pj - 1] * (c * c) + d[nj - 1] * (s * s);
			d[nj - 1] = d[pj - 1] * (s * s) + d[nj - 1] * (c * c);
			d[pj - 1] = t;
			k2 = k2 - 1;
			i = 1;
L90:
			if (k2 + i <= n) {
				if (d[pj - 1] < d[indxp[k2 + i - 1] - 1]) {
					indxp[k2 + i - 2] = indxp[k2 + i - 1];
					indxp[k2 + i - 1] = pj;
					i = i + 1;
					goto L90;
				}
				else {
					indxp[k2 + i - 2] = pj;
				}
			}
			else {
				indxp[k2 + i - 2] = pj;
			}
			pj = nj;
		}
		else {
			k = k + 1;
			dlambda[k - 1] = d[pj - 1];
			w[k - 1] = z[pj - 1];
			indxp[k - 1] = pj;
			pj = nj;
		}
	}
	goto L80;
L100:

	// 最後の固有値を記録する
	k = k + 1;
	dlambda[k - 1] = d[pj - 1];
	w[k - 1] = z[pj - 1];
	indxp[k - 1] = pj;

	// 列の種類ごとの本数と，種類ごとにまとめる置換
	for (j = 1; j <= 4; j++) {
		ctot[j - 1] = 0;
	}
	for (j = 1; j <= n; j++) {
		ct = coltyp[j - 1];
		ctot[ct - 1] = ctot[ct - 1] + 1;
	}

	psm[0] = 1;
	psm[1] = 1 + ctot[0];
	psm[2] = psm[1] + ctot[1];
	psm[3] = psm[2] + ctot[2];
	k = n - ctot[3];

	for (j = 1; j <= n; j++) {
		js = indxp[j - 1];
		ct = coltyp[js - 1];
		indx[psm[ct - 1] - 1] = js;
		indxc[psm[ct - 1] - 1] = j;
		psm[ct - 1] = psm[ct - 1] + 1;
	}

	// 固有 vector を種類ごとに q2 に詰め，固有値を z に並べる
	i = 1;
	iq1 = 1;
	iq2 = 1 + (ctot[0] + ctot[1]) * n1;
	for (j = 1; j <= ctot[0]; j++) {
		js = indx[i - 1];
		tcopy(n1, Q + LQ * (js - 1), 1, q2 + (iq1 - 1), 1);
		z[i - 1] = d[js - 1];
		i = i + 1;
		iq1 = iq1 + n1;
	}

	for (j = 1; j <= ctot[1]; j++) {
		js = indx[i - 1];
		tcopy(n1, Q + LQ * (js - 1), 1, q2 + (iq1 - 1), 1);
		tcopy(n2, Q + n1 + LQ * (js - 1), 1, q2 + (iq2 - 1), 1);
		z[i - 1] = d[js - 1];
		i = i + 1;
		iq1 = iq1 + n1;
		iq2 = iq2 + n2;
	}

	for (j = 1; j <= ctot[2]; j++) {
		js = indx[i - 1];
		tcopy(n2, Q + n1 + LQ * (js - 1), 1, q2 + (iq2 - 1), 1);
		z[i - 1] = d[js - 1];
		i = i + 1;
		iq2 = iq2 + n2;
	}

	iq1 = iq2;
	for (j = 1; j <= ctot[3]; j++) {
		js = indx[i - 1];
		tcopy(n, Q + LQ * (js - 1), 1, q2 + (iq2 - 1), 1);
		iq2 = iq2 + n;
		z[i - 1] = d[js - 1];
		i = i + 1;
	}

	// deflation した固有値と固有 vector を，d と Q の後ろに戻す
	if (k < n) {
		tlacpy('A', n, ctot[3], q2 + (iq1 - 1), n, Q + LQ * k, ldq);
		tcopy(n - k, z + k, 1, d + k, 1);
	}

	// tlaed3 に渡すため，種類ごとの本数を coltyp に入れる
	for (j = 1; j <= 4; j++) {
		coltyp[j - 1] = ctot[j - 1];
	}

L190:
	return 0;
}

// secular 方程式の根 (tlaed4) から，併合した問題の固有値と固有 vector を作る
// (compz = I の経路)．返り値 1 は固有値が収束しなかったことを表す
template <typename T>
inline int tlaed3(const int k, const int n, const int n1, T* d, T* Q, const int ldq,
	const T rho, const T* dlambda, const T* q2, const int* indx, const int* ctot,
	T* w, T* s) {
	namespace detail = tlapack_detail;
	if (k < 0) {
		detail::tlapack_error("tlaed3: invalid k");
	}
	else if (n < k) {
		detail::tlapack_error("tlaed3: invalid n");
	}
	else if (ldq < std::max(1, n)) {
		detail::tlapack_error("tlaed3: invalid ldq");
	}
	if (k == 0) {
		return 0;
	}
	const std::size_t LQ = static_cast<std::size_t>(ldq);
	int i, ii, iq2, j, n12, n2, n23;
	T temp;
	int info = 0;

	for (j = 1; j <= k; j++) {
		info = tlaed4(k, j, dlambda, w, Q + LQ * (j - 1), rho, d[j - 1]);
		// 収束しなければ戻る
		if (info != 0) {
			goto L120;
		}
	}

	if (k == 1) {
		goto L110;
	}
	if (k == 2) {
		for (j = 1; j <= k; j++) {
			w[0] = Q[LQ * (j - 1)];
			w[1] = Q[1 + LQ * (j - 1)];
			ii = indx[0];
			Q[LQ * (j - 1)] = w[ii - 1];
			ii = indx[1];
			Q[1 + LQ * (j - 1)] = w[ii - 1];
		}
		goto L110;
	}

	// z を作り直す (Gu-Eisenstat)
	tcopy(k, w, 1, s, 1);

	// w を Q の対角で初期化する
	tcopy(k, Q, ldq + 1, w, 1);
	for (j = 1; j <= k; j++) {
		for (i = 1; i <= j - 1; i++) {
			w[i - 1] = w[i - 1] * (Q[(i - 1) + LQ * (j - 1)] / (dlambda[i - 1] - dlambda[j - 1]));
		}
		for (i = j + 1; i <= k; i++) {
			w[i - 1] = w[i - 1] * (Q[(i - 1) + LQ * (j - 1)] / (dlambda[i - 1] - dlambda[j - 1]));
		}
	}
	for (i = 1; i <= k; i++) {
		w[i - 1] = detail::f_sign(tlapack_detail::tsqrt(-w[i - 1]), s[i - 1]);
	}

	// 併合した問題の固有 vector (secular 方程式の側)
	for (j = 1; j <= k; j++) {
		for (i = 1; i <= k; i++) {
			s[i - 1] = w[i - 1] / Q[(i - 1) + LQ * (j - 1)];
		}
		temp = tnrm2(k, s, 1);
		for (i = 1; i <= k; i++) {
			ii = indx[i - 1];
			Q[(i - 1) + LQ * (j - 1)] = s[ii - 1] / temp;
		}
	}

	// 元の問題の固有 vector に戻す
L110:

	n2 = n - n1;
	n12 = ctot[0] + ctot[1];
	n23 = ctot[1] + ctot[2];

	tlacpy('A', n23, k, Q + ctot[0], ldq, s, n23);
	iq2 = n1 * n12 + 1;
	if (n23 != 0) {
		tgemm('N', 'N', n2, k, n23, T(1), q2 + (iq2 - 1), n2, s, n23, T(0), Q + n1, ldq);
	}
	else {
		tlaset('A', n2, k, T(0), T(0), Q + n1, ldq);
	}

	tlacpy('A', n12, k, Q, ldq, s, n12);
	if (n12 != 0) {
		tgemm('N', 'N', n1, k, n12, T(1), q2, n1, s, n12, T(0), Q, ldq);
	}
	else {
		tlaset('A', n1, k, T(0), T(0), Q, ldq);
	}

L120:
	return info;
}

// 併合 1 回 (compz = I の経路)．隣り合う 2 つの部分問題の固有分解から，
// 階数 1 の修正 (rho) を加えた問題の固有分解を作る．
// indxq は入力で 2 つの部分問題を，出力で併合した問題を昇順に並べる置換 (1-based)
template <typename T>
inline int tlaed1(const int n, T* d, T* Q, const int ldq, int* indxq, T& rho,
	const int cutpnt) {
	namespace detail = tlapack_detail;
	if (n < 0) {
		detail::tlapack_error("tlaed1: invalid n");
	}
	else if (ldq < std::max(1, n)) {
		detail::tlapack_error("tlaed1: invalid ldq");
	}
	else if (std::min(1, n / 2) > cutpnt || (n / 2) < cutpnt) {
		detail::tlapack_error("tlaed1: invalid cutpnt");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t LQ = static_cast<std::size_t>(ldq);
	std::vector<T> work(4 * static_cast<std::size_t>(n) + static_cast<std::size_t>(n) * static_cast<std::size_t>(n), T(0));
	std::vector<int> iwork(4 * static_cast<std::size_t>(n));
	int coltyp, i, idlmda, indx, indxc, indxp, iq2, is, iw, iz, k, n1, n2, zpp1;
	int info = 0;

	// 作業領域の中の位置
	iz = 1;
	idlmda = iz + n;
	iw = idlmda + n;
	iq2 = iw + n;

	indx = 1;
	indxc = indx + n;
	coltyp = indxc + n;
	indxp = coltyp + n;

	// z vector: 上の部分問題の固有 vector の最後の行と，下の部分問題の最初の行
	tcopy(cutpnt, Q + (cutpnt - 1), ldq, work.data() + (iz - 1), 1);
	zpp1 = cutpnt + 1;
	tcopy(n - cutpnt, Q + (zpp1 - 1) + LQ * (zpp1 - 1), ldq, work.data() + (iz + cutpnt - 1), 1);

	// deflation
	info = tlaed2(k, n, cutpnt, d, Q, ldq, indxq, rho, work.data() + (iz - 1),
		work.data() + (idlmda - 1), work.data() + (iw - 1), work.data() + (iq2 - 1),
		iwork.data() + (indx - 1), iwork.data() + (indxc - 1), iwork.data() + (indxp - 1),
		iwork.data() + (coltyp - 1));

	if (info != 0) {
		goto L20;
	}

	// secular 方程式を解く
	if (k != 0) {
		is = (iwork[coltyp - 1] + iwork[coltyp]) * cutpnt +
			(iwork[coltyp] + iwork[coltyp + 1]) * (n - cutpnt) + iq2;
		info = tlaed3(k, n, cutpnt, d, Q, ldq, rho, work.data() + (idlmda - 1),
			work.data() + (iq2 - 1), iwork.data() + (indxc - 1), iwork.data() + (coltyp - 1),
			work.data() + (iw - 1), work.data() + (is - 1));
		if (info != 0) {
			goto L20;
		}

		// 併合した問題を昇順に並べる置換
		n1 = k;
		n2 = n - k;
		tlamrg(n1, n2, d, 1, -1, indxq);
	}
	else {
		for (i = 1; i <= n; i++) {
			indxq[i - 1] = i;
		}
	}

L20:
	return info;
}

// 併合に使う z vector を，下の階で記録した変換 (Givens 回転，置換，固有 vector) から作る
// (tlaed7 の経路)．givcol, givnum は 2 行の配列 (column-major)
template <typename T>
inline int tlaeda(const int n, const int tlvls, const int curlvl, const int curpbm,
	const int* prmptr, const int* perm, const int* givptr, const int* givcol, const T* givnum,
	const T* q, const int* qptr, T* z, T* ztemp) {
	namespace detail = tlapack_detail;
	if (n < 0) {
		detail::tlapack_error("tlaeda: invalid n");
	}
	if (n == 0) {
		return 0;
	}
	int bsiz1, bsiz2, curr, i, k, mid, psiz1, psiz2, ptr, zptr1;

	// 中央の位置
	mid = n / 2 + 1;

	// 現在の問題の，一番下の階での位置
	ptr = 1;
	curr = ptr + curpbm * (1 << curlvl) + (1 << (curlvl - 1)) - 1;

	// 一番下の階の 2 つの固有 vector の行列から，上の最後の行と下の最初の行を取り出す．
	// 行列の大きさは，記録してある「大きさの 2 乗」から戻す
	bsiz1 = static_cast<int>(0.5 + std::sqrt(static_cast<double>(qptr[curr] - qptr[curr - 1])));
	bsiz2 = static_cast<int>(0.5 + std::sqrt(static_cast<double>(qptr[curr + 1] - qptr[curr])));
	for (k = 1; k <= mid - bsiz1 - 1; k++) {
		z[k - 1] = T(0);
	}
	tcopy(bsiz1, q + (qptr[curr - 1] + bsiz1 - 2), bsiz1, z + (mid - bsiz1 - 1), 1);
	tcopy(bsiz2, q + (qptr[curr] - 1), bsiz2, z + (mid - 1), 1);
	for (k = mid + bsiz2; k <= n; k++) {
		z[k - 1] = T(0);
	}

	// 上の階へ順に，回転・置換を掛け，その階の固有 vector の行列を掛ける
	ptr = (1 << tlvls) + 1;
	for (k = 1; k <= curlvl - 1; k++) {
		curr = ptr + curpbm * (1 << (curlvl - k)) + (1 << (curlvl - k - 1)) - 1;
		psiz1 = prmptr[curr] - prmptr[curr - 1];
		psiz2 = prmptr[curr + 1] - prmptr[curr];
		zptr1 = mid - psiz1;

		// Givens 回転を掛ける
		for (i = givptr[curr - 1]; i <= givptr[curr] - 1; i++) {
			trot(1, z + (zptr1 + givcol[2 * (i - 1)] - 2), 1,
				z + (zptr1 + givcol[1 + 2 * (i - 1)] - 2), 1, givnum[2 * (i - 1)],
				givnum[1 + 2 * (i - 1)]);
		}
		for (i = givptr[curr]; i <= givptr[curr + 1] - 1; i++) {
			trot(1, z + (mid - 2 + givcol[2 * (i - 1)]), 1,
				z + (mid - 2 + givcol[1 + 2 * (i - 1)]), 1, givnum[2 * (i - 1)],
				givnum[1 + 2 * (i - 1)]);
		}
		psiz1 = prmptr[curr] - prmptr[curr - 1];
		psiz2 = prmptr[curr + 1] - prmptr[curr];
		for (i = 0; i <= psiz1 - 1; i++) {
			ztemp[i] = z[zptr1 + perm[prmptr[curr - 1] + i - 1] - 2];
		}
		for (i = 0; i <= psiz2 - 1; i++) {
			ztemp[psiz1 + i] = z[mid + perm[prmptr[curr] + i - 1] - 2];
		}

		// 固有 vector の行列 (deflation していない部分) を掛ける
		bsiz1 = static_cast<int>(0.5 + std::sqrt(static_cast<double>(qptr[curr] - qptr[curr - 1])));
		bsiz2 = static_cast<int>(0.5 + std::sqrt(static_cast<double>(qptr[curr + 1] - qptr[curr])));
		if (bsiz1 > 0) {
			tgemv('T', bsiz1, bsiz1, T(1), q + (qptr[curr - 1] - 1), bsiz1, ztemp, 1, T(0),
				z + (zptr1 - 1), 1);
		}
		tcopy(psiz1 - bsiz1, ztemp + bsiz1, 1, z + (zptr1 + bsiz1 - 1), 1);
		if (bsiz2 > 0) {
			tgemv('T', bsiz2, bsiz2, T(1), q + (qptr[curr] - 1), bsiz2, ztemp + psiz1, 1, T(0),
				z + (mid - 1), 1);
		}
		tcopy(psiz2 - bsiz2, ztemp + (psiz1 + bsiz2), 1, z + (mid + bsiz2 - 1), 1);

		ptr = ptr + (1 << (tlvls - k));
	}

	return 0;
}

// 併合の前処理 (deflation)．tlaed7 の経路で，置換 perm と Givens 回転
// (givptr, givcol, givnum) を記録する．icompq: 0 固有値だけ, 1 固有 vector も
template <typename T>
inline int tlaed8(const int icompq, int& k, const int n, const int qsiz, T* d, T* Q,
	const int ldq, int* indxq, T& rho, const int cutpnt, T* z, T* dlambda,
	T* q2, const int ldq2, T* w, int* perm, int& givptr, int* givcol, T* givnum,
	int* indxp, int* indx) {
	namespace detail = tlapack_detail;
	if (icompq < 0 || icompq > 1) {
		detail::tlapack_error("tlaed8: invalid icompq");
	}
	else if (n < 0) {
		detail::tlapack_error("tlaed8: invalid n");
	}
	else if (icompq == 1 && qsiz < n) {
		detail::tlapack_error("tlaed8: invalid qsiz");
	}
	else if (ldq < std::max(1, n)) {
		detail::tlapack_error("tlaed8: invalid ldq");
	}
	else if (cutpnt < std::min(1, n) || cutpnt > n) {
		detail::tlapack_error("tlaed8: invalid cutpnt");
	}
	else if (ldq2 < std::max(1, n)) {
		detail::tlapack_error("tlaed8: invalid ldq2");
	}

	// n = 0 のときにも値を持たせる
	givptr = 0;

	if (n == 0) {
		return 0;
	}
	const std::size_t LQ = static_cast<std::size_t>(ldq);
	const std::size_t LQ2 = static_cast<std::size_t>(ldq2);
	int i, imax, j, jlam, jmax, jp, k2, n1, n1p1, n2;
	T c, eps, s, t, tau, tol;
	jlam = 0;

	n1 = cutpnt;
	n2 = n - n1;
	n1p1 = n1 + 1;

	if (rho < T(0)) {
		tscal(n2, -T(1), z + (n1p1 - 1), 1);
	}

	// z を正規化する
	t = T(1) / tlapack_detail::tsqrt(T(2));
	for (j = 1; j <= n; j++) {
		indx[j - 1] = j;
	}
	tscal(n, t, z, 1);
	rho = tlapack_detail::tabs(T(2) * rho);

	// 2 つの部分問題の固有値を昇順に並べる
	for (i = cutpnt + 1; i <= n; i++) {
		indxq[i - 1] = indxq[i - 1] + cutpnt;
	}
	for (i = 1; i <= n; i++) {
		dlambda[i - 1] = d[indxq[i - 1] - 1];
		w[i - 1] = z[indxq[i - 1] - 1];
	}
	i = 1;
	j = cutpnt + 1;
	tlamrg(n1, n2, dlambda, 1, 1, indx);
	for (i = 1; i <= n; i++) {
		d[i - 1] = dlambda[indx[i - 1] - 1];
		z[i - 1] = w[indx[i - 1] - 1];
	}

	// deflation の許容値
	imax = itamax(n, z, 1) + 1;
	jmax = itamax(n, d, 1) + 1;
	eps = tlamch<T>('E');
	tol = T(8) * eps * tlapack_detail::tabs(d[jmax - 1]);

	// rho が小さければ，並べ替えるだけで終わる
	if (rho * tlapack_detail::tabs(z[imax - 1]) <= tol) {
		k = 0;
		if (icompq == 0) {
			for (j = 1; j <= n; j++) {
				perm[j - 1] = indxq[indx[j - 1] - 1];
			}
		}
		else {
			for (j = 1; j <= n; j++) {
				perm[j - 1] = indxq[indx[j - 1] - 1];
				tcopy(qsiz, Q + LQ * (perm[j - 1] - 1), 1, q2 + LQ2 * (j - 1), 1);
			}
			tlacpy('A', qsiz, n, q2, ldq2, Q, ldq);
		}
		return 0;
	}

	k = 0;
	k2 = n + 1;
	for (j = 1; j <= n; j++) {
		if (rho * tlapack_detail::tabs(z[j - 1]) <= tol) {
			// z の成分が小さいので deflation する
			k2 = k2 - 1;
			indxp[k2 - 1] = j;
			if (j == n) {
				goto L110;
			}
		}
		else {
			jlam = j;
			goto L80;
		}
	}
L80:
	j = j + 1;
	if (j > n) {
		goto L100;
	}
	if (rho * tlapack_detail::tabs(z[j - 1]) <= tol) {
		// z の成分が小さいので deflation する
		k2 = k2 - 1;
		indxp[k2 - 1] = j;
	}
	else {
		// 固有値が近ければ，Givens 回転で z の成分を 1 つ消して deflation する
		s = z[jlam - 1];
		c = z[j - 1];
		tau = tlapy2(c, s);
		t = d[j - 1] - d[jlam - 1];
		c = c / tau;
		s = -(s / tau);
		if (tlapack_detail::tabs(t * c * s) <= tol) {
			z[j - 1] = tau;
			z[jlam - 1] = T(0);

			// 回転を記録する
			givptr = givptr + 1;
			givcol[2 * (givptr - 1)] = indxq[indx[jlam - 1] - 1];
			givcol[1 + 2 * (givptr - 1)] = indxq[indx[j - 1] - 1];
			givnum[2 * (givptr - 1)] = c;
			givnum[1 + 2 * (givptr - 1)] = s;
			if (icompq == 1) {
				trot(qsiz, Q + LQ * (indxq[indx[jlam - 1] - 1] - 1), 1,
					Q + LQ * (indxq[indx[j - 1] - 1] - 1), 1, c, s);
			}
			t = d[jlam - 1] * c * c + d[j - 1] * s * s;
			d[j - 1] = d[jlam - 1] * s * s + d[j - 1] * c * c;
			d[jlam - 1] = t;
			k2 = k2 - 1;
			i = 1;
L90:
			if (k2 + i <= n) {
				if (d[jlam - 1] < d[indxp[k2 + i - 1] - 1]) {
					indxp[k2 + i - 2] = indxp[k2 + i - 1];
					indxp[k2 + i - 1] = jlam;
					i = i + 1;
					goto L90;
				}
				else {
					indxp[k2 + i - 2] = jlam;
				}
			}
			else {
				indxp[k2 + i - 2] = jlam;
			}
			jlam = j;
		}
		else {
			k = k + 1;
			w[k - 1] = z[jlam - 1];
			dlambda[k - 1] = d[jlam - 1];
			indxp[k - 1] = jlam;
			jlam = j;
		}
	}
	goto L80;
L100:

	// 最後の固有値を記録する
	k = k + 1;
	w[k - 1] = z[jlam - 1];
	dlambda[k - 1] = d[jlam - 1];
	indxp[k - 1] = jlam;

L110:

	// 固有値と固有 vector を，deflation していないものが前に来る順に並べる
	if (icompq == 0) {
		for (j = 1; j <= n; j++) {
			jp = indxp[j - 1];
			dlambda[j - 1] = d[jp - 1];
			perm[j - 1] = indxq[indx[jp - 1] - 1];
		}
	}
	else {
		for (j = 1; j <= n; j++) {
			jp = indxp[j - 1];
			dlambda[j - 1] = d[jp - 1];
			perm[j - 1] = indxq[indx[jp - 1] - 1];
			tcopy(qsiz, Q + LQ * (perm[j - 1] - 1), 1, q2 + LQ2 * (j - 1), 1);
		}
	}

	// deflation した固有値と固有 vector を，d と Q の後ろに戻す
	if (k < n) {
		if (icompq == 0) {
			tcopy(n - k, dlambda + k, 1, d + k, 1);
		}
		else {
			tcopy(n - k, dlambda + k, 1, d + k, 1);
			tlacpy('A', qsiz, n - k, q2 + LQ2 * k, ldq2, Q + LQ * k, ldq);
		}
	}

	return 0;
}

// secular 方程式の根 (tlaed4) から，併合した問題の固有値と，secular 方程式の側の
// 固有 vector (S) を作る (tlaed7 の経路)．返り値 1 は固有値が収束しなかったことを表す
template <typename T>
inline int tlaed9(const int k, const int kstart, const int kstop, const int n, T* d,
	T* Q, const int ldq, const T rho, const T* dlambda, T* w,
	T* S, const int lds) {
	namespace detail = tlapack_detail;
	if (k < 0) {
		detail::tlapack_error("tlaed9: invalid k");
	}
	else if (kstart < 1 || kstart > std::max(1, k)) {
		detail::tlapack_error("tlaed9: invalid kstart");
	}
	else if (std::max(1, kstop) < kstart || kstop > std::max(1, k)) {
		detail::tlapack_error("tlaed9: invalid kstop");
	}
	else if (n < k) {
		detail::tlapack_error("tlaed9: invalid n");
	}
	else if (ldq < std::max(1, k)) {
		detail::tlapack_error("tlaed9: invalid ldq");
	}
	else if (lds < std::max(1, k)) {
		detail::tlapack_error("tlaed9: invalid lds");
	}
	if (k == 0) {
		return 0;
	}
	const std::size_t LQ = static_cast<std::size_t>(ldq);
	const std::size_t LS = static_cast<std::size_t>(lds);
	int i, j;
	T temp;
	int info = 0;

	for (j = kstart; j <= kstop; j++) {
		info = tlaed4(k, j, dlambda, w, Q + LQ * (j - 1), rho, d[j - 1]);
		// 収束しなければ戻る
		if (info != 0) {
			goto L120;
		}
	}

	if (k == 1 || k == 2) {
		for (i = 1; i <= k; i++) {
			for (j = 1; j <= k; j++) {
				S[(j - 1) + LS * (i - 1)] = Q[(j - 1) + LQ * (i - 1)];
			}
		}
		goto L120;
	}

	// z を作り直す (Gu-Eisenstat)
	tcopy(k, w, 1, S, 1);

	// w を Q の対角で初期化する
	tcopy(k, Q, ldq + 1, w, 1);
	for (j = 1; j <= k; j++) {
		for (i = 1; i <= j - 1; i++) {
			w[i - 1] = w[i - 1] * (Q[(i - 1) + LQ * (j - 1)] / (dlambda[i - 1] - dlambda[j - 1]));
		}
		for (i = j + 1; i <= k; i++) {
			w[i - 1] = w[i - 1] * (Q[(i - 1) + LQ * (j - 1)] / (dlambda[i - 1] - dlambda[j - 1]));
		}
	}
	for (i = 1; i <= k; i++) {
		w[i - 1] = detail::f_sign(tlapack_detail::tsqrt(-w[i - 1]), S[i - 1]);
	}

	// 併合した問題の固有 vector (secular 方程式の側)
	for (j = 1; j <= k; j++) {
		for (i = 1; i <= k; i++) {
			Q[(i - 1) + LQ * (j - 1)] = w[i - 1] / Q[(i - 1) + LQ * (j - 1)];
		}
		temp = tnrm2(k, Q + LQ * (j - 1), 1);
		for (i = 1; i <= k; i++) {
			S[(i - 1) + LS * (j - 1)] = Q[(i - 1) + LQ * (j - 1)] / temp;
		}
	}

L120:
	return info;
}

// 併合 1 回 (compz = V の経路と，固有値だけの経路)．変換を qstore, qptr, prmptr, perm,
// givptr, givcol, givnum に記録し，上の階の併合 (tlaeda) で使う．
// icompq: 0 固有値だけ, 1 元の密行列の固有 vector も
template <typename T>
inline int tlaed7(const int icompq, const int n, const int qsiz, const int tlvls, const int curlvl,
	const int curpbm, T* d, T* Q, const int ldq, int* indxq, T& rho, const int cutpnt,
	T* qstore, int* qptr, int* prmptr, int* perm, int* givptr, int* givcol, T* givnum) {
	namespace detail = tlapack_detail;
	if (icompq < 0 || icompq > 1) {
		detail::tlapack_error("tlaed7: invalid icompq");
	}
	else if (n < 0) {
		detail::tlapack_error("tlaed7: invalid n");
	}
	else if (icompq == 1 && qsiz < n) {
		detail::tlapack_error("tlaed7: invalid qsiz");
	}
	else if (ldq < std::max(1, n)) {
		detail::tlapack_error("tlaed7: invalid ldq");
	}
	else if (std::min(1, n) > cutpnt || n < cutpnt) {
		detail::tlapack_error("tlaed7: invalid cutpnt");
	}
	if (n == 0) {
		return 0;
	}
	int coltyp, curr, i, idlmda, indx, indxc, indxp, iq2, is, iw, iz, k, ldq2, n1, n2, ptr;
	int info = 0;

	// 作業領域の中の位置
	if (icompq == 1) {
		ldq2 = qsiz;
	}
	else {
		ldq2 = n;
	}

	iz = 1;
	idlmda = iz + n;
	iw = idlmda + n;
	iq2 = iw + n;
	is = iq2 + n * ldq2;

	indx = 1;
	indxc = indx + n;
	coltyp = indxc + n;
	indxp = coltyp + n;
	(void)indxc;
	(void)coltyp;

	std::vector<T> work(3 * static_cast<std::size_t>(n) +
		static_cast<std::size_t>(n) * static_cast<std::size_t>(ldq2) +
		static_cast<std::size_t>(n) * static_cast<std::size_t>(n), T(0));
	std::vector<int> iwork(4 * static_cast<std::size_t>(n));

	// z vector を作る．現在の問題の記録の位置は curr
	ptr = 1 + (1 << tlvls);
	for (i = 1; i <= curlvl - 1; i++) {
		ptr = ptr + (1 << (tlvls - i));
	}
	curr = ptr + curpbm;
	tlaeda(n, tlvls, curlvl, curpbm, prmptr, perm, givptr, givcol, givnum, qstore, qptr,
		work.data() + (iz - 1), work.data() + (iz + n - 1));

	// 一番上の階では，下の階の記録はもう要らないので，先頭から書く
	if (curlvl == tlvls) {
		qptr[curr - 1] = 1;
		prmptr[curr - 1] = 1;
		givptr[curr - 1] = 1;
	}

	// deflation
	tlaed8(icompq, k, n, qsiz, d, Q, ldq, indxq, rho, cutpnt, work.data() + (iz - 1),
		work.data() + (idlmda - 1), work.data() + (iq2 - 1), ldq2, work.data() + (iw - 1),
		perm + (prmptr[curr - 1] - 1), givptr[curr], givcol + 2 * (givptr[curr - 1] - 1),
		givnum + 2 * (givptr[curr - 1] - 1), iwork.data() + (indxp - 1), iwork.data() + (indx - 1));
	prmptr[curr] = prmptr[curr - 1] + n;
	givptr[curr] = givptr[curr] + givptr[curr - 1];

	// secular 方程式を解く
	if (k != 0) {
		info = tlaed9(k, 1, k, n, d, work.data() + (is - 1), k, rho, work.data() + (idlmda - 1),
			work.data() + (iw - 1), qstore + (qptr[curr - 1] - 1), k);
		if (info != 0) {
			goto L30;
		}
		if (icompq == 1) {
			tgemm('N', 'N', qsiz, k, k, T(1), work.data() + (iq2 - 1), ldq2,
				qstore + (qptr[curr - 1] - 1), k, T(0), Q, ldq);
		}
		qptr[curr] = qptr[curr - 1] + k * k;

		// 併合した問題を昇順に並べる置換
		n1 = k;
		n2 = n - k;
		tlamrg(n1, n2, d, 1, -1, indxq);
	}
	else {
		qptr[curr] = qptr[curr - 1];
		for (i = 1; i <= n; i++) {
			indxq[i - 1] = i;
		}
	}

L30:
	return info;
}

// 分割統治の本体．三重対角行列 (d, e) を eigdc_smlsiz 以下の部分問題に分け，
// 下の階から順に併合する．
// icompq: 0 固有値だけ, 1 Q (qsiz 行) に三重対角化の直交行列を渡して元の密行列の
//         固有 vector, 2 三重対角行列の固有 vector (Q は n x n)．
// qstore は icompq = 1 のときだけ使う (qsiz 行 n 列)．
// 返り値 > 0 は，部分問題の固有値が求まらなかったことを表す
// (行と列が info/(n+1) から info%(n+1) の部分行列)
template <typename T>
inline int tlaed0(const int icompq, const int qsiz, const int n, T* d, T* e,
	T* Q, const int ldq, T* qstore, const int ldqs) {
	namespace detail = tlapack_detail;
	if (icompq < 0 || icompq > 2) {
		detail::tlapack_error("tlaed0: invalid icompq");
	}
	else if ((icompq == 1) && (qsiz < std::max(0, n))) {
		detail::tlapack_error("tlaed0: invalid qsiz");
	}
	else if (n < 0) {
		detail::tlapack_error("tlaed0: invalid n");
	}
	else if (ldq < std::max(1, n)) {
		detail::tlapack_error("tlaed0: invalid ldq");
	}
	else if (ldqs < std::max(1, n)) {
		detail::tlapack_error("tlaed0: invalid ldqs");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t LQ = static_cast<std::size_t>(ldq);
	const std::size_t LQS = static_cast<std::size_t>(ldqs);
	const std::size_t NZ = static_cast<std::size_t>(n);
	int curlvl, curprb, curr, i, igivcl, igivnm, igivpt, indxq, iperm, iprmpt, iq, iqptr, j, k,
		lgn, matsiz, msd2, smlsiz, smm1, spm1, spm2, submat, subpbs, tlvls;
	int info = 0;
	curprb = 0;
	igivcl = 0;
	igivnm = 0;
	igivpt = 0;
	iperm = 0;
	iprmpt = 0;
	iq = 0;
	iqptr = 0;
	submat = 0;
	matsiz = 0;
	std::vector<T> work;
	std::vector<int> iwork(3 + 5 * NZ);

	smlsiz = detail::eigdc_smlsiz;

	// 部分問題の大きさを決める (iwork(1:subpbs) に，はじめ大きさ，後で末尾の位置)
	iwork[0] = n;
	subpbs = 1;
	tlvls = 0;
L10:
	if (iwork[subpbs - 1] > smlsiz) {
		for (j = subpbs; j >= 1; j--) {
			iwork[2 * j - 1] = (iwork[j - 1] + 1) / 2;
			iwork[2 * j - 2] = iwork[j - 1] / 2;
		}
		tlvls = tlvls + 1;
		subpbs = 2 * subpbs;
		goto L10;
	}
	for (j = 2; j <= subpbs; j++) {
		iwork[j - 1] = iwork[j - 1] + iwork[j - 2];
	}

	// 階数 1 の修正を引いて，部分問題に分ける
	spm1 = subpbs - 1;
	for (i = 1; i <= spm1; i++) {
		submat = iwork[i - 1] + 1;
		smm1 = submat - 1;
		d[smm1 - 1] = d[smm1 - 1] - tlapack_detail::tabs(e[smm1 - 1]);
		d[submat - 1] = d[submat - 1] - tlapack_detail::tabs(e[smm1 - 1]);
	}

	indxq = 4 * n + 3;
	if (icompq != 2) {
		// 作業領域の中の位置．lgn は 2^lgn >= n となる最小の整数
		lgn = 0;
		while ((1 << lgn) < n) {
			lgn = lgn + 1;
		}
		iprmpt = indxq + n + 1;
		iperm = iprmpt + n * lgn;
		iqptr = iperm + n * lgn;
		igivpt = iqptr + n + 2;
		igivcl = igivpt + n * lgn;

		igivnm = 1;
		iq = igivnm + 2 * n * lgn;
		iwork.resize(6 + 6 * NZ + 5 * NZ * static_cast<std::size_t>(lgn), T(0));
		work.resize(2 * NZ * static_cast<std::size_t>(lgn) + NZ * NZ + 1, T(0));

		// 記録の先頭の位置を初期化する
		for (i = 0; i <= subpbs; i++) {
			iwork[iprmpt + i - 1] = 1;
			iwork[igivpt + i - 1] = 1;
		}
		iwork[iqptr - 1] = 1;
	}
	else {
		work.resize(NZ + NZ * NZ, T(0));
	}

	// 一番下の階の部分問題を QR 法で解く
	curr = 0;
	for (i = 0; i <= spm1; i++) {
		if (i == 0) {
			submat = 1;
			matsiz = iwork[0];
		}
		else {
			submat = iwork[i - 1] + 1;
			matsiz = iwork[i] - iwork[i - 1];
		}
		if (icompq == 2) {
			info = tsteqr('I', matsiz, d + (submat - 1), e + (submat - 1),
				Q + (submat - 1) + LQ * (submat - 1), ldq);
			if (info != 0) {
				goto L130;
			}
		}
		else {
			info = tsteqr('I', matsiz, d + (submat - 1), e + (submat - 1),
				work.data() + (iq - 1 + iwork[iqptr + curr - 1] - 1), matsiz);
			if (info != 0) {
				goto L130;
			}
			if (icompq == 1) {
				tgemm('N', 'N', qsiz, matsiz, matsiz, T(1), Q + LQ * (submat - 1), ldq,
					work.data() + (iq - 1 + iwork[iqptr + curr - 1] - 1), matsiz, T(0),
					qstore + LQS * (submat - 1), ldqs);
			}
			iwork[iqptr + curr] = iwork[iqptr + curr - 1] + matsiz * matsiz;
			curr = curr + 1;
		}
		k = 1;
		for (j = submat; j <= iwork[i]; j++) {
			iwork[indxq + j - 1] = k;
			k = k + 1;
		}
	}

	// 隣り合う部分問題を 2 つずつ併合することを，1 つになるまで繰り返す
	curlvl = 1;
L80:
	if (subpbs > 1) {
		spm2 = subpbs - 2;
		for (i = 0; i <= spm2; i += 2) {
			if (i == 0) {
				submat = 1;
				matsiz = iwork[1];
				msd2 = iwork[0];
				curprb = 0;
			}
			else {
				submat = iwork[i - 1] + 1;
				matsiz = iwork[i + 1] - iwork[i - 1];
				msd2 = matsiz / 2;
				curprb = curprb + 1;
			}

			// 併合する．icompq = 2 は tlaed1，それ以外は変換を記録する tlaed7
			if (icompq == 2) {
				info = tlaed1(matsiz, d + (submat - 1), Q + (submat - 1) + LQ * (submat - 1),
					ldq, iwork.data() + (indxq + submat - 1), e[submat + msd2 - 2], msd2);
			}
			else {
				info = tlaed7(icompq, matsiz, qsiz, tlvls, curlvl, curprb,
					d + (submat - 1), qstore + LQS * (submat - 1), ldqs,
					iwork.data() + (indxq + submat - 1), e[submat + msd2 - 2],
					msd2, work.data() + (iq - 1), iwork.data() + (iqptr - 1),
					iwork.data() + (iprmpt - 1), iwork.data() + (iperm - 1),
					iwork.data() + (igivpt - 1), iwork.data() + (igivcl - 1),
					work.data() + (igivnm - 1));
			}
			if (info != 0) {
				goto L130;
			}
			iwork[i / 2] = iwork[i + 1];
		}
		subpbs = subpbs / 2;
		curlvl = curlvl + 1;
		goto L80;
	}

	// 固有値を昇順に並べ，固有 vector も同じ順に並べる
	if (icompq == 1) {
		for (i = 1; i <= n; i++) {
			j = iwork[indxq + i - 1];
			work[i - 1] = d[j - 1];
			tcopy(qsiz, qstore + LQS * (j - 1), 1, Q + LQ * (i - 1), 1);
		}
		tcopy(n, work.data(), 1, d, 1);
	}
	else if (icompq == 2) {
		for (i = 1; i <= n; i++) {
			j = iwork[indxq + i - 1];
			work[i - 1] = d[j - 1];
			tcopy(n, Q + LQ * (j - 1), 1, work.data() + NZ * i, 1);
		}
		tcopy(n, work.data(), 1, d, 1);
		tlacpy('A', n, n, work.data() + NZ, n, Q, ldq);
	}
	else {
		for (i = 1; i <= n; i++) {
			j = iwork[indxq + i - 1];
			work[i - 1] = d[j - 1];
		}
		tcopy(n, work.data(), 1, d, 1);
	}
	goto L140;

L130:
	info = submat * (n + 1) + submat + matsiz - 1;

L140:
	return info;
}

// 対称三重対角行列の固有値・固有 vector (分割統治法)
// compz: 'N' 固有値のみ, 'V' Z に三重対角化の Q を渡して元行列の固有 vector,
//        'I' Z を単位行列で初期化して三重対角行列の固有 vector
// 返り値 > 0 は，行と列が info/(n+1) から info%(n+1) の部分行列の固有値が
// 求まらなかったことを表す
template <typename T>
inline int tstedc(const char compz, const int n, T* d, T* e,
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
		detail::tlapack_error("tstedc: invalid compz");
		return -1;
	}
	if (n < 0 || ldz < 1 || (icompz > 0 && ldz < std::max(1, n))) {
		detail::tlapack_error("tstedc: invalid n/ldz");
	}
	if (n == 0) {
		return 0;
	}
	if (n == 1) {
		if (icompz != 0) {
			Z[0] = T(1);
		}
		return 0;
	}
	const std::size_t LZ = static_cast<std::size_t>(ldz);
	const std::size_t NZ = static_cast<std::size_t>(n);
	const int smlsiz = detail::eigdc_smlsiz;
	int finish, i, ii, j, k, m, start, storez, strtrw;
	T eps, orgnrm, p, tiny;
	int info = 0;
	std::vector<T> work;

	// 固有値だけなら dsterf を使う
	if (icompz == 0) {
		info = tsterf(n, d, e);
		goto L50;
	}

	// 小さい問題は QR 法 (dsteqr) で解く
	if (n <= smlsiz) {
		info = tsteqr(compz, n, d, e, Z, ldz);
	}
	else {
		// compz = 'V' のとき，作業領域の前半 (n x n) に固有 vector を置く
		if (icompz == 1) {
			storez = 1 + n * n;
			work.resize(NZ * NZ + NZ * static_cast<std::size_t>(smlsiz), T(0));
		}
		else {
			storez = 1;
		}

		if (icompz == 2) {
			tlaset('F', n, n, T(0), T(1), Z, ldz);
		}

		// 零行列なら終わり
		orgnrm = tlanst('M', n, d, e);
		if (orgnrm == T(0)) {
			goto L50;
		}

		eps = tlamch<T>('E');

		start = 1;

L10:
		if (start <= n) {
			// 副対角が小さいところで，独立な部分行列に分ける
			finish = start;
L20:
			if (finish < n) {
				tiny = eps * tlapack_detail::tsqrt(tlapack_detail::tabs(d[finish - 1])) *
					tlapack_detail::tsqrt(tlapack_detail::tabs(d[finish]));
				if (tlapack_detail::tabs(e[finish - 1]) > tiny) {
					finish = finish + 1;
					goto L20;
				}
			}

			// 部分行列 (start から finish) を解く
			m = finish - start + 1;
			if (m == 1) {
				start = finish + 1;
				goto L10;
			}
			if (m > smlsiz) {
				// 最大絶対値が 1 になるように scaling して，分割統治で解く
				orgnrm = tlanst('M', m, d + (start - 1), e + (start - 1));
				tlascl('G', 0, 0, orgnrm, T(1), m, 1, d + (start - 1), m);
				tlascl('G', 0, 0, orgnrm, T(1), m - 1, 1, e + (start - 1), m - 1);

				if (icompz == 1) {
					strtrw = 1;
				}
				else {
					strtrw = start;
				}
				info = tlaed0(icompz, n, m, d + (start - 1), e + (start - 1),
					Z + (strtrw - 1) + LZ * (start - 1), ldz, work.data(), n);
				if (info != 0) {
					info = (info / (m + 1) + start - 1) * (n + 1) +
						(info % (m + 1)) + start - 1;
					goto L50;
				}

				// scaling を戻す
				tlascl('G', 0, 0, T(1), orgnrm, m, 1, d + (start - 1), m);
			}
			else {
				if (icompz == 1) {
					// 小さい部分行列は QR 法で解き，Z の対応する列に固有 vector を掛ける
					info = tsteqr('I', m, d + (start - 1), e + (start - 1), work.data(), m);
					tlacpy('A', n, m, Z + LZ * (start - 1), ldz, work.data() + (storez - 1), n);
					tgemm('N', 'N', n, m, m, T(1), work.data() + (storez - 1), n, work.data(), m,
						T(0), Z + LZ * (start - 1), ldz);
				}
				else if (icompz == 2) {
					info = tsteqr('I', m, d + (start - 1), e + (start - 1),
						Z + (start - 1) + LZ * (start - 1), ldz);
				}
				else {
					info = tsterf(m, d + (start - 1), e + (start - 1));
				}
				if (info != 0) {
					info = start * (n + 1) + finish;
					goto L50;
				}
			}

			start = finish + 1;
			goto L10;
		}

		if (icompz == 0) {
			// 固有値を昇順に並べる
			tlasrt('I', n, d);
		}
		else {
			// 固有値を昇順に並べ，固有 vector も入れ替える (選択法)
			for (ii = 2; ii <= n; ii++) {
				i = ii - 1;
				k = i;
				p = d[i - 1];
				for (j = ii; j <= n; j++) {
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
	}

L50:
	return info;
}

// 対称行列の固有値・固有 vector (driver，固有 vector は分割統治法)
// jobz: 'N' 固有値のみ, 'V' 固有 vector も計算 (A に上書き)
// w: 固有値 (昇順, n 個)
template <typename T>
inline int tsyevd(const char jobz, const char uplo, const int n, T* A, const int lda,
	T* w) {
	if (tlapack_detail::double_route<T>::on) {
		return tlapack_detail::double_route<T>::syevd(jobz, uplo, n, A, lda, w);
	}
	namespace detail = tlapack_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	const bool lower = detail::option_is(uplo, 'L');
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::tlapack_error("tsyevd: invalid jobz");
	}
	if (!lower && !detail::option_is(uplo, 'U')) {
		detail::tlapack_error("tsyevd: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsyevd: invalid n/lda");
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
		// 三重対角行列の固有 vector を分割統治で求め，三重対角化の Q を掛ける
		std::vector<T> work(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), T(0));
		info = tstedc('I', n, w, e.data(), work.data(), n);
		tormtr('L', uplo, 'N', n, n, A, lda, tau.data(), work.data(), n);
		tlacpy('A', n, n, work.data(), n, A, lda);
	}
	if (iscale == 1) {
		tscal(n, T(1) / sigma, w, 1);
	}
	return info;
}

// 一般化対称固有値問題 A*x = lambda*B*x (itype=1), A*B*x = lambda*x (2),
// B*A*x = lambda*x (3) の driver (固有 vector は分割統治法．B は対称正定値，上書きされる)
template <typename T>
inline int tsygvd(const int itype, const char jobz, const char uplo, const int n,
	T* A, const int lda, T* B, const int ldb, T* w) {
	if (tlapack_detail::double_route<T>::on) {
		return tlapack_detail::double_route<T>::sygvd(itype, jobz, uplo, n, A, lda, B, ldb, w);
	}
	namespace detail = tlapack_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	const bool upper = detail::option_is(uplo, 'U');
	if (itype < 1 || itype > 3) {
		detail::tlapack_error("tsygvd: invalid itype");
	}
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::tlapack_error("tsygvd: invalid jobz");
	}
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsygvd: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsygvd: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	int info = tpotrf(uplo, n, B, ldb);
	if (info != 0) {
		return n + info;
	}
	tsygst(itype, uplo, n, A, lda, B, ldb);
	info = tsyevd(jobz, uplo, n, A, lda, w);
	if (wantz && info == 0) {
		if (itype == 1 || itype == 2) {
			const char trans = upper ? 'N' : 'T';
			ttrsm('L', uplo, trans, 'N', n, n, T(1), B, ldb, A, lda);
		}
		else if (itype == 3) {
			const char trans = upper ? 'T' : 'N';
			ttrmm('L', uplo, trans, 'N', n, n, T(1), B, ldb, A, lda);
		}
	}
	return info;
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_EIGDC_HPP
