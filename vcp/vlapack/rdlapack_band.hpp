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

#ifndef VBLAS_RDLAPACK_BAND_HPP
#define VBLAS_RDLAPACK_BAND_HPP

// 帯行列の連立一次方程式 (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dgbtf2/dgbtrf/dgbtrs/dgbsv/
// dpbtf2/dpbtrf/dpbtrs/dpbsv を移植．
// AB は LAPACK の帯格納 (dgb 系: ldab >= 2*kl+ku+1 の LU 用格納)．
// ipiv は 0-based．

#include "rdlapack_sy.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// 帯行列の LU 分解 (unblocked)．ipiv は 0-based
inline int rdgbtf2(const int m, const int n, const int kl, const int ku,
	double* AB, const int ldab, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const int kv = ku + kl;
	if (m < 0 || n < 0 || kl < 0 || ku < 0 || ldab < kl + kv + 1) {
		detail::rdlapack_error("rdgbtf2: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(ldab);
	int info = 0;
	// fill-in 部分を 0 に初期化
	for (int j = ku + 2; j <= std::min(kv, n); j++) {
		for (int i = kv - j + 2; i <= kl; i++) {
			AB[(i - 1) + L * (j - 1)] = 0.0;
		}
	}
	int ju = 1;
	for (int j = 1; j <= std::min(m, n); j++) {
		if (j + kv <= n) {
			for (int i = 1; i <= kl; i++) {
				AB[(i - 1) + L * (j + kv - 1)] = 0.0;
			}
		}
		const int km = std::min(kl, m - j);
		const int jp = idamax(km + 1, AB + kv + L * (j - 1), 1) + 1;
		ipiv[j - 1] = (jp + j - 1) - 1; // 0-based
		if (AB[(kv + jp - 1) + L * (j - 1)] != 0.0) {
			ju = std::max(ju, std::min(j + ku + jp - 1, n));
			if (jp != 1) {
				dswap(ju - j + 1, AB + (kv + jp - 1) + L * (j - 1), ldab - 1,
					AB + kv + L * (j - 1), ldab - 1);
			}
			if (km > 0) {
				rdscal(km, 1.0 / AB[kv + L * (j - 1)], AB + (kv + 1) + L * (j - 1), 1, rm);
				if (ju > j) {
					rdger(km, ju - j, -1.0, AB + (kv + 1) + L * (j - 1), 1,
						AB + (kv - 1) + L * j, ldab - 1, AB + kv + L * j, ldab - 1, rm);
				}
			}
		}
		else {
			if (info == 0) {
				info = j;
			}
		}
	}
	return info;
}

// 帯行列の LU 分解 (blocked)．ipiv は 0-based
inline int rdgbtrf(const int m, const int n, const int kl, const int ku,
	double* AB, const int ldab, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const int kv = ku + kl;
	if (m < 0 || n < 0 || kl < 0 || ku < 0 || ldab < kl + kv + 1) {
		detail::rdlapack_error("rdgbtrf: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(ldab);
	const int nbmax = 64;
	const int ldwork = nbmax + 1;
	int nb = std::min(detail::ilaenv(1, "GBTRF", m, n, kl, ku), nbmax);
	if (nb <= 1 || nb > kl) {
		return rdgbtf2(m, n, kl, ku, AB, ldab, ipiv, rounding_mode);
	}
	std::vector<double> work13(static_cast<std::size_t>(ldwork) * nbmax, 0.0);
	std::vector<double> work31(static_cast<std::size_t>(ldwork) * nbmax, 0.0);
	double* W13 = work13.data();
	double* W31 = work31.data();
	const std::size_t LW = static_cast<std::size_t>(ldwork);
	int info = 0;
	for (int j = ku + 2; j <= std::min(kv, n); j++) {
		for (int i = kv - j + 2; i <= kl; i++) {
			AB[(i - 1) + L * (j - 1)] = 0.0;
		}
	}
	int ju = 1;
	for (int j = 1; j <= std::min(m, n); j += nb) {
		const int jb = std::min(nb, std::min(m, n) - j + 1);
		const int i2 = std::min(kl - jb, m - j - jb + 1);
		const int i3 = std::min(jb, m - j - kl + 1);
		// panel の分解 (pivot は band 内と work31 に分けて適用)
		for (int jj = j; jj <= j + jb - 1; jj++) {
			if (jj + kv <= n) {
				for (int i = 1; i <= kl; i++) {
					AB[(i - 1) + L * (jj + kv - 1)] = 0.0;
				}
			}
			const int km = std::min(kl, m - jj);
			const int jp = idamax(km + 1, AB + kv + L * (jj - 1), 1) + 1;
			ipiv[jj - 1] = jp + jj - j; // 一時的に block 内 1-based
			if (AB[(kv + jp - 1) + L * (jj - 1)] != 0.0) {
				ju = std::max(ju, std::min(jj + ku + jp - 1, n));
				if (jp != 1) {
					if (jp + jj - 1 < j + kl) {
						dswap(jb, AB + (kv + jj - j) + L * (j - 1), ldab - 1,
							AB + (kv + jp + jj - j - 1) + L * (j - 1), ldab - 1);
					}
					else {
						dswap(jj - j, AB + (kv + jj - j) + L * (j - 1), ldab - 1,
							W31 + (jp + jj - j - kl - 1), ldwork);
						dswap(j + jb - jj, AB + kv + L * (jj - 1), ldab - 1,
							AB + (kv + jp - 1) + L * (jj - 1), ldab - 1);
					}
				}
				rdscal(km, 1.0 / AB[kv + L * (jj - 1)], AB + (kv + 1) + L * (jj - 1), 1, rm);
				const int jm = std::min(ju, j + jb - 1);
				if (jm > jj) {
					rdger(km, jm - jj, -1.0, AB + (kv + 1) + L * (jj - 1), 1,
						AB + (kv - 1) + L * jj, ldab - 1, AB + kv + L * jj, ldab - 1, rm);
				}
			}
			else {
				if (info == 0) {
					info = jj;
				}
			}
			const int nw = std::min(jj - j + 1, i3);
			if (nw > 0) {
				dcopy(nw, AB + (kv + kl - jj + j) + L * (jj - 1), 1, W31 + LW * (jj - j), 1);
			}
		}
		if (j + jb <= n) {
			const int j2 = std::min(ju - j + 1, kv) - jb;
			const int j3 = std::max(0, ju - j - kv + 1);
			// block 内 pivot を A12 (band 内) に適用 (ipiv は 0-based block 相対へ)
			{
				std::vector<int> piv0(static_cast<std::size_t>(jb));
				for (int i = 0; i < jb; i++) {
					piv0[i] = ipiv[j - 1 + i] - 1;
				}
				dlaswp(j2, AB + (kv - jb) + L * (j + jb - 1), ldab - 1, 0, jb, piv0.data(), 1);
			}
			for (int i = j; i <= j + jb - 1; i++) {
				ipiv[i - 1] = (ipiv[i - 1] + j - 1) - 1; // 0-based global へ
			}
			const int k2 = j - 1 + jb + j2;
			for (int i = 1; i <= j3; i++) {
				const int jj = k2 + i;
				for (int ii = j + i - 1; ii <= j + jb - 1; ii++) {
					const int ip = ipiv[ii - 1] + 1; // 1-based
					if (ip != ii) {
						const double temp = AB[(kv + ii - jj) + L * (jj - 1)];
						AB[(kv + ii - jj) + L * (jj - 1)] = AB[(kv + ip - jj) + L * (jj - 1)];
						AB[(kv + ip - jj) + L * (jj - 1)] = temp;
					}
				}
			}
			if (j2 > 0) {
				rdtrsm('L', 'L', 'N', 'U', jb, j2, 1.0, AB + kv + L * (j - 1), ldab - 1,
					AB + (kv - jb) + L * (j + jb - 1), ldab - 1, rm);
				if (i2 > 0) {
					rdgemm('N', 'N', i2, j2, jb, -1.0, AB + (kv + jb) + L * (j - 1), ldab - 1,
						AB + (kv - jb) + L * (j + jb - 1), ldab - 1, 1.0,
						AB + kv + L * (j + jb - 1), ldab - 1, rm);
				}
				if (i3 > 0) {
					rdgemm('N', 'N', i3, j2, jb, -1.0, W31, ldwork,
						AB + (kv - jb) + L * (j + jb - 1), ldab - 1, 1.0,
						AB + (kv + kl - jb) + L * (j + jb - 1), ldab - 1, rm);
				}
			}
			if (j3 > 0) {
				for (int jj = 1; jj <= j3; jj++) {
					for (int ii = jj; ii <= jb; ii++) {
						W13[(ii - 1) + LW * (jj - 1)] = AB[(ii - jj) + L * (jj + j + kv - 2)];
					}
				}
				rdtrsm('L', 'L', 'N', 'U', jb, j3, 1.0, AB + kv + L * (j - 1), ldab - 1,
					W13, ldwork, rm);
				if (i2 > 0) {
					rdgemm('N', 'N', i2, j3, jb, -1.0, AB + (kv + jb) + L * (j - 1), ldab - 1,
						W13, ldwork, 1.0, AB + jb + L * (j + kv - 1), ldab - 1, rm);
				}
				if (i3 > 0) {
					rdgemm('N', 'N', i3, j3, jb, -1.0, W31, ldwork, W13, ldwork, 1.0,
						AB + kl + L * (j + kv - 1), ldab - 1, rm);
				}
				for (int jj = 1; jj <= j3; jj++) {
					for (int ii = jj; ii <= jb; ii++) {
						AB[(ii - jj) + L * (jj + j + kv - 2)] = W13[(ii - 1) + LW * (jj - 1)];
					}
				}
			}
		}
		else {
			for (int i = j; i <= j + jb - 1; i++) {
				ipiv[i - 1] = (ipiv[i - 1] + j - 1) - 1; // 0-based global へ
			}
		}
		// panel 内で work31 へ退避した部分の swap を band 構造に合うよう戻す
		for (int jj = j + jb - 1; jj >= j; jj--) {
			const int jp = (ipiv[jj - 1] + 1) - jj + 1; // 1-based 相対
			if (jp != 1) {
				if (jp + jj - 1 < j + kl) {
					dswap(jj - j, AB + (kv + jj - j) + L * (j - 1), ldab - 1,
						AB + (kv + jp + jj - j - 1) + L * (j - 1), ldab - 1);
				}
				else {
					dswap(jj - j, AB + (kv + jj - j) + L * (j - 1), ldab - 1,
						W31 + (jp + jj - j - kl - 1), ldwork);
				}
			}
			const int nw = std::min(i3, jj - j + 1);
			if (nw > 0) {
				dcopy(nw, W31 + LW * (jj - j), 1, AB + (kv + kl - jj + j) + L * (jj - 1), 1);
			}
		}
	}
	return info;
}

// 帯行列の LU 分解を使った連立一次方程式の求解 (ipiv は 0-based)
inline int rdgbtrs(const char trans, const int n, const int kl, const int ku, const int nrhs,
	const double* AB, const int ldab, const int* ipiv, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool notran = detail::option_is(trans, 'N');
	if (!notran && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::rdlapack_error("rdgbtrs: invalid trans");
	}
	if (n < 0 || kl < 0 || ku < 0 || nrhs < 0 || ldab < 2 * kl + ku + 1 || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdgbtrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(ldab);
	const std::size_t LB = static_cast<std::size_t>(ldb);
	const int kd = ku + kl + 1;
	const bool lnoti = kl > 0;
	if (notran) {
		if (lnoti) {
			for (int j = 1; j <= n - 1; j++) {
				const int lm = std::min(kl, n - j);
				const int l = ipiv[j - 1] + 1; // 1-based
				if (l != j) {
					dswap(nrhs, B + (l - 1), ldb, B + (j - 1), ldb);
				}
				rdger(lm, nrhs, -1.0, AB + kd + L * (j - 1), 1, B + (j - 1), ldb, B + j, ldb, rm);
			}
		}
		for (int i = 1; i <= nrhs; i++) {
			rdtbsv('U', 'N', 'N', n, kl + ku, AB, ldab, B + LB * (i - 1), 1, rm);
		}
	}
	else {
		for (int i = 1; i <= nrhs; i++) {
			rdtbsv('U', 'T', 'N', n, kl + ku, AB, ldab, B + LB * (i - 1), 1, rm);
		}
		if (lnoti) {
			for (int j = n - 1; j >= 1; j--) {
				const int lm = std::min(kl, n - j);
				rdgemv('T', lm, nrhs, -1.0, B + j, ldb, AB + kd + L * (j - 1), 1,
					1.0, B + (j - 1), ldb, rm);
				const int l = ipiv[j - 1] + 1;
				if (l != j) {
					dswap(nrhs, B + (l - 1), ldb, B + (j - 1), ldb);
				}
			}
		}
	}
	return 0;
}

// 帯行列の連立一次方程式 A*X = B の求解 (driver)．ipiv は 0-based
inline int rdgbsv(const int n, const int kl, const int ku, const int nrhs,
	double* AB, const int ldab, int* ipiv, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0 || kl < 0 || ku < 0 || nrhs < 0 || ldab < 2 * kl + ku + 1 || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdgbsv: invalid argument");
	}
	const int info = rdgbtrf(n, n, kl, ku, AB, ldab, ipiv, rounding_mode);
	if (info == 0) {
		rdgbtrs('N', n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, rounding_mode);
	}
	return info;
}

// 対称正定値帯行列の Cholesky 分解 (unblocked)
inline int rdpbtf2(const char uplo, const int n, const int kd, double* AB, const int ldab, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpbtf2: invalid uplo");
	}
	if (n < 0 || kd < 0 || ldab < kd + 1) {
		detail::rdlapack_error("rdpbtf2: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(ldab);
	const int kld = std::max(1, ldab - 1);
	if (upper) {
		for (int j = 1; j <= n; j++) {
			double ajj = AB[kd + L * (j - 1)];
			if (ajj <= 0.0) {
				return j;
			}
			ajj = std::sqrt(ajj);
			AB[kd + L * (j - 1)] = ajj;
			const int kn = std::min(kd, n - j);
			if (kn > 0) {
				rdscal(kn, 1.0 / ajj, AB + (kd - 1) + L * j, kld, rm);
				rdsyr('U', kn, -1.0, AB + (kd - 1) + L * j, kld, AB + kd + L * j, kld, rm);
			}
		}
	}
	else {
		for (int j = 1; j <= n; j++) {
			double ajj = AB[L * (j - 1)];
			if (ajj <= 0.0) {
				return j;
			}
			ajj = std::sqrt(ajj);
			AB[L * (j - 1)] = ajj;
			const int kn = std::min(kd, n - j);
			if (kn > 0) {
				rdscal(kn, 1.0 / ajj, AB + 1 + L * (j - 1), 1, rm);
				rdsyr('L', kn, -1.0, AB + 1 + L * (j - 1), 1, AB + L * j, kld, rm);
			}
		}
	}
	return 0;
}

// 対称正定値帯行列の Cholesky 分解 (blocked)
inline int rdpbtrf(const char uplo, const int n, const int kd, double* AB, const int ldab, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpbtrf: invalid uplo");
	}
	if (n < 0 || kd < 0 || ldab < kd + 1) {
		detail::rdlapack_error("rdpbtrf: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int rm = rounding_mode;
	const std::size_t L = static_cast<std::size_t>(ldab);
	const int nbmax = 32;
	const int ldwork = nbmax + 1;
	int nb = std::min(detail::ilaenv(1, "PBTRF", n, kd, -1, -1), nbmax);
	if (nb <= 1 || nb > kd) {
		return rdpbtf2(uplo, n, kd, AB, ldab, rounding_mode);
	}
	std::vector<double> work(static_cast<std::size_t>(ldwork) * nbmax, 0.0);
	double* W = work.data();
	const std::size_t LW = static_cast<std::size_t>(ldwork);
	if (upper) {
		for (int i = 1; i <= n; i += nb) {
			const int ib = std::min(nb, n - i + 1);
			const int ii = rdpotf2('U', ib, AB + kd + L * (i - 1), ldab - 1, rm);
			if (ii != 0) {
				return i + ii - 1;
			}
			if (i + ib <= n) {
				const int i2 = std::min(kd - ib, n - i - ib + 1);
				const int i3 = std::min(ib, n - i - kd + 1);
				if (i2 > 0) {
					rdtrsm('L', 'U', 'T', 'N', ib, i2, 1.0, AB + kd + L * (i - 1), ldab - 1,
						AB + (kd - ib) + L * (i + ib - 1), ldab - 1, rm);
					rdsyrk('U', 'T', i2, ib, -1.0, AB + (kd - ib) + L * (i + ib - 1), ldab - 1,
						1.0, AB + kd + L * (i + ib - 1), ldab - 1, rm);
				}
				if (i3 > 0) {
					for (int jj = 1; jj <= i3; jj++) {
						for (int iii = jj; iii <= ib; iii++) {
							W[(iii - 1) + LW * (jj - 1)] = AB[(iii - jj) + L * (jj + i + kd - 2)];
						}
					}
					rdtrsm('L', 'U', 'T', 'N', ib, i3, 1.0, AB + kd + L * (i - 1), ldab - 1,
						W, ldwork, rm);
					if (i2 > 0) {
						rdgemm('T', 'N', i2, i3, ib, -1.0, AB + (kd - ib) + L * (i + ib - 1), ldab - 1,
							W, ldwork, 1.0, AB + ib + L * (i + kd - 1), ldab - 1, rm);
					}
					rdsyrk('U', 'T', i3, ib, -1.0, W, ldwork, 1.0, AB + kd + L * (i + kd - 1), ldab - 1, rm);
					for (int jj = 1; jj <= i3; jj++) {
						for (int iii = jj; iii <= ib; iii++) {
							AB[(iii - jj) + L * (jj + i + kd - 2)] = W[(iii - 1) + LW * (jj - 1)];
						}
					}
				}
			}
		}
	}
	else {
		for (int i = 1; i <= n; i += nb) {
			const int ib = std::min(nb, n - i + 1);
			const int ii = rdpotf2('L', ib, AB + L * (i - 1), ldab - 1, rm);
			if (ii != 0) {
				return i + ii - 1;
			}
			if (i + ib <= n) {
				const int i2 = std::min(kd - ib, n - i - ib + 1);
				const int i3 = std::min(ib, n - i - kd + 1);
				if (i2 > 0) {
					rdtrsm('R', 'L', 'T', 'N', i2, ib, 1.0, AB + L * (i - 1), ldab - 1,
						AB + ib + L * (i - 1), ldab - 1, rm);
					rdsyrk('L', 'N', i2, ib, -1.0, AB + ib + L * (i - 1), ldab - 1,
						1.0, AB + L * (i + ib - 1), ldab - 1, rm);
				}
				if (i3 > 0) {
					for (int jj = 1; jj <= ib; jj++) {
						for (int iii = 1; iii <= std::min(jj, i3); iii++) {
							W[(iii - 1) + LW * (jj - 1)] = AB[(kd - jj + iii) + L * (jj + i - 2)];
						}
					}
					rdtrsm('R', 'L', 'T', 'N', i3, ib, 1.0, AB + L * (i - 1), ldab - 1,
						W, ldwork, rm);
					if (i2 > 0) {
						rdgemm('N', 'T', i3, i2, ib, -1.0, W, ldwork, AB + ib + L * (i - 1), ldab - 1,
							1.0, AB + (kd - ib) + L * (i + ib - 1), ldab - 1, rm);
					}
					rdsyrk('L', 'N', i3, ib, -1.0, W, ldwork, 1.0, AB + L * (i + kd - 1), ldab - 1, rm);
					for (int jj = 1; jj <= ib; jj++) {
						for (int iii = 1; iii <= std::min(jj, i3); iii++) {
							AB[(kd - jj + iii) + L * (jj + i - 2)] = W[(iii - 1) + LW * (jj - 1)];
						}
					}
				}
			}
		}
	}
	return 0;
}

// 対称正定値帯行列の Cholesky 分解を使った求解
inline int rdpbtrs(const char uplo, const int n, const int kd, const int nrhs,
	const double* AB, const int ldab, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpbtrs: invalid uplo");
	}
	if (n < 0 || kd < 0 || nrhs < 0 || ldab < kd + 1 || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdpbtrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	const std::size_t LB = static_cast<std::size_t>(ldb);
	for (int j = 0; j < nrhs; j++) {
		if (upper) {
			rdtbsv('U', 'T', 'N', n, kd, AB, ldab, B + LB * j, 1, rounding_mode);
			rdtbsv('U', 'N', 'N', n, kd, AB, ldab, B + LB * j, 1, rounding_mode);
		}
		else {
			rdtbsv('L', 'N', 'N', n, kd, AB, ldab, B + LB * j, 1, rounding_mode);
			rdtbsv('L', 'T', 'N', n, kd, AB, ldab, B + LB * j, 1, rounding_mode);
		}
	}
	return 0;
}

// 対称正定値帯行列の連立一次方程式 A*X = B の求解 (driver)
inline int rdpbsv(const char uplo, const int n, const int kd, const int nrhs,
	double* AB, const int ldab, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpbsv: invalid uplo");
	}
	if (n < 0 || kd < 0 || nrhs < 0 || ldab < kd + 1 || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdpbsv: invalid argument");
	}
	const int info = rdpbtrf(uplo, n, kd, AB, ldab, rounding_mode);
	if (info == 0) {
		rdpbtrs(uplo, n, kd, nrhs, AB, ldab, B, ldb, rounding_mode);
	}
	return info;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_BAND_HPP
