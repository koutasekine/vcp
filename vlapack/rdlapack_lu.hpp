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

#ifndef VBLAS_RDLAPACK_LU_HPP
#define VBLAS_RDLAPACK_LU_HPP

// LU 分解と三角行列 routine (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dgetf2/dgetrf2/dgetrf/dgetrs/dgesv/dgetri/
// dtrti2/dtrtri/dtrtrs/dlauu2/dlauum を移植．
// ipiv は 0-based．返り値は LAPACK の INFO (0: 成功, >0: LAPACK と同じ意味の
// 1-based 位置)．引数 error は std::invalid_argument を投げる．

#include "rdlapack_aux.hpp"

#pragma STDC FENV_ACCESS ON

// LU 分解 (部分 pivot 付き, unblocked)
inline int rdgetf2(const int m, const int n, double* A, const int lda, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgetf2: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const double sfmin = DBL_MIN;
	int info = 0;
	for (int j = 0; j < std::min(m, n); j++) {
		const int jp = j + idamax(m - j, A + j + L * j, 1);
		ipiv[j] = jp;
		if (A[jp + L * j] != 0.0) {
			if (jp != j) {
				dswap(n, A + j, lda, A + jp, lda);
			}
			if (j < m - 1) {
				if (std::fabs(A[j + L * j]) >= sfmin) {
					rdscal(m - j - 1, 1.0 / A[j + L * j], A + (j + 1) + L * j, 1, rounding_mode);
				}
				else {
					for (int i = 1; i <= m - j - 1; i++) {
						A[(j + i) + L * j] = A[(j + i) + L * j] / A[j + L * j];
					}
				}
			}
		}
		else if (info == 0) {
			info = j + 1;
		}
		if (j < std::min(m, n) - 1) {
			rdger(m - j - 1, n - j - 1, -1.0, A + (j + 1) + L * j, 1, A + j + L * (j + 1), lda,
				A + (j + 1) + L * (j + 1), lda, rounding_mode);
		}
	}
	return info;
}

// LU 分解 (部分 pivot 付き, recursive)
inline int rdgetrf2(const int m, const int n, double* A, const int lda, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgetrf2: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	int info = 0;
	if (m == 1) {
		ipiv[0] = 0;
		if (A[0] == 0.0) {
			info = 1;
		}
	}
	else if (n == 1) {
		const double sfmin = DBL_MIN;
		const int i = idamax(m, A, 1);
		ipiv[0] = i;
		if (A[i] != 0.0) {
			if (i != 0) {
				const double temp = A[0];
				A[0] = A[i];
				A[i] = temp;
			}
			if (std::fabs(A[0]) >= sfmin) {
				rdscal(m - 1, 1.0 / A[0], A + 1, 1, rounding_mode);
			}
			else {
				for (int k = 1; k < m; k++) {
					A[k] = A[k] / A[0];
				}
			}
		}
		else {
			info = 1;
		}
	}
	else {
		const int n1 = std::min(m, n) / 2;
		const int n2 = n - n1;
		int iinfo = rdgetrf2(m, n1, A, lda, ipiv, rounding_mode);
		if (info == 0 && iinfo > 0) {
			info = iinfo;
		}
		dlaswp(n2, A + L * n1, lda, 0, n1, ipiv, 1);
		rdtrsm('L', 'L', 'N', 'U', n1, n2, 1.0, A, lda, A + L * n1, lda, rounding_mode);
		rdgemm('N', 'N', m - n1, n2, n1, -1.0, A + n1, lda, A + L * n1, lda, 1.0,
			A + n1 + L * n1, lda, rounding_mode);
		iinfo = rdgetrf2(m - n1, n2, A + n1 + L * n1, lda, ipiv + n1, rounding_mode);
		if (info == 0 && iinfo > 0) {
			info = iinfo + n1;
		}
		for (int i = n1; i < std::min(m, n); i++) {
			ipiv[i] = ipiv[i] + n1;
		}
		dlaswp(n1, A, lda, n1, std::min(m, n), ipiv, 1);
	}
	return info;
}

// LU 分解 (部分 pivot 付き, blocked)
inline int rdgetrf(const int m, const int n, double* A, const int lda, int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgetrf: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GETRF", m, n, -1, -1);
	int info = 0;
	if (nb <= 1 || nb >= std::min(m, n)) {
		return rdgetrf2(m, n, A, lda, ipiv, rounding_mode);
	}
	for (int j = 0; j < std::min(m, n); j += nb) {
		const int jb = std::min(std::min(m, n) - j, nb);
		const int iinfo = rdgetrf2(m - j, jb, A + j + L * j, lda, ipiv + j, rounding_mode);
		if (info == 0 && iinfo > 0) {
			info = iinfo + j;
		}
		for (int i = j; i < std::min(m, j + jb); i++) {
			ipiv[i] = j + ipiv[i];
		}
		dlaswp(j, A, lda, j, j + jb, ipiv, 1);
		if (j + jb < n) {
			dlaswp(n - j - jb, A + L * (j + jb), lda, j, j + jb, ipiv, 1);
			rdtrsm('L', 'L', 'N', 'U', jb, n - j - jb, 1.0, A + j + L * j, lda,
				A + j + L * (j + jb), lda, rounding_mode);
			if (j + jb < m) {
				rdgemm('N', 'N', m - j - jb, n - j - jb, jb, -1.0, A + (j + jb) + L * j, lda,
					A + j + L * (j + jb), lda, 1.0, A + (j + jb) + L * (j + jb), lda, rounding_mode);
			}
		}
	}
	return info;
}

// LU 分解を使った連立一次方程式の求解 (A は rdgetrf の出力)
inline int rdgetrs(const char trans, const int n, const int nrhs, const double* A, const int lda,
	const int* ipiv, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool notran = detail::option_is(trans, 'N');
	if (!notran && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::rdlapack_error("rdgetrs: invalid trans");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdgetrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	if (notran) {
		dlaswp(nrhs, B, ldb, 0, n, ipiv, 1);
		rdtrsm('L', 'L', 'N', 'U', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
		rdtrsm('L', 'U', 'N', 'N', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
	}
	else {
		rdtrsm('L', 'U', 'T', 'N', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
		rdtrsm('L', 'L', 'T', 'U', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
		dlaswp(nrhs, B, ldb, 0, n, ipiv, -1);
	}
	return 0;
}

// 連立一次方程式 A*X = B の求解 (driver)
inline int rdgesv(const int n, const int nrhs, double* A, const int lda, int* ipiv,
	double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdgesv: invalid argument");
	}
	const int info = rdgetrf(n, n, A, lda, ipiv, rounding_mode);
	if (info == 0) {
		rdgetrs('N', n, nrhs, A, lda, ipiv, B, ldb, rounding_mode);
	}
	return info;
}

// 三角行列の逆行列 (unblocked)
inline int rdtrti2(const char uplo, const char diag, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	const bool nounit = detail::option_is(diag, 'N');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdtrti2: invalid uplo");
	}
	if (!nounit && !detail::option_is(diag, 'U')) {
		detail::rdlapack_error("rdtrti2: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdtrti2: invalid n/lda");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int j = 0; j < n; j++) {
			double ajj;
			if (nounit) {
				A[j + L * j] = 1.0 / A[j + L * j];
				ajj = -A[j + L * j];
			}
			else {
				ajj = -1.0;
			}
			rdtrmv('U', 'N', diag, j, A, lda, A + L * j, 1, rounding_mode);
			rdscal(j, ajj, A + L * j, 1, rounding_mode);
		}
	}
	else {
		for (int j = n - 1; j >= 0; j--) {
			double ajj;
			if (nounit) {
				A[j + L * j] = 1.0 / A[j + L * j];
				ajj = -A[j + L * j];
			}
			else {
				ajj = -1.0;
			}
			if (j < n - 1) {
				rdtrmv('L', 'N', diag, n - j - 1, A + (j + 1) + L * (j + 1), lda,
					A + (j + 1) + L * j, 1, rounding_mode);
				rdscal(n - j - 1, ajj, A + (j + 1) + L * j, 1, rounding_mode);
			}
		}
	}
	return 0;
}

// 三角行列の逆行列 (blocked)
inline int rdtrtri(const char uplo, const char diag, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	const bool nounit = detail::option_is(diag, 'N');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdtrtri: invalid uplo");
	}
	if (!nounit && !detail::option_is(diag, 'U')) {
		detail::rdlapack_error("rdtrtri: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdtrtri: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (nounit) {
		for (int i = 0; i < n; i++) {
			if (A[i + L * i] == 0.0) {
				return i + 1;
			}
		}
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const int nb = detail::ilaenv(1, "TRTRI", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		rdtrti2(uplo, diag, n, A, lda, rounding_mode);
	}
	else if (upper) {
		for (int j = 0; j < n; j += nb) {
			const int jb = std::min(nb, n - j);
			rdtrmm('L', 'U', 'N', diag, j, jb, 1.0, A, lda, A + L * j, lda, rounding_mode);
			rdtrsm('R', 'U', 'N', diag, j, jb, -1.0, A + j + L * j, lda, A + L * j, lda, rounding_mode);
			rdtrti2('U', diag, jb, A + j + L * j, lda, rounding_mode);
		}
	}
	else {
		const int nn = ((n - 1) / nb) * nb;
		for (int j = nn; j >= 0; j -= nb) {
			const int jb = std::min(nb, n - j);
			if (j + jb < n) {
				rdtrmm('L', 'L', 'N', diag, n - j - jb, jb, 1.0, A + (j + jb) + L * (j + jb), lda,
					A + (j + jb) + L * j, lda, rounding_mode);
				rdtrsm('R', 'L', 'N', diag, n - j - jb, jb, -1.0, A + j + L * j, lda,
					A + (j + jb) + L * j, lda, rounding_mode);
			}
			rdtrti2('L', diag, jb, A + j + L * j, lda, rounding_mode);
		}
	}
	return 0;
}

// 三角行列の連立一次方程式 (特異性 check 付き)
inline int rdtrtrs(const char uplo, const char trans, const char diag, const int n, const int nrhs,
	const double* A, const int lda, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool nounit = detail::option_is(diag, 'N');
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdtrtrs: invalid uplo");
	}
	if (!detail::option_is(trans, 'N') && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::rdlapack_error("rdtrtrs: invalid trans");
	}
	if (!nounit && !detail::option_is(diag, 'U')) {
		detail::rdlapack_error("rdtrtrs: invalid diag");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdtrtrs: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	if (nounit) {
		for (int i = 0; i < n; i++) {
			if (A[i + static_cast<std::size_t>(lda) * i] == 0.0) {
				return i + 1;
			}
		}
	}
	rdtrsm('L', uplo, trans, diag, n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
	return 0;
}

// U*U^T (uplo='U') または L^T*L (uplo='L') の三角部分計算 (unblocked)
inline int rdlauu2(const char uplo, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdlauu2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdlauu2: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int i = 0; i < n; i++) {
			const double aii = A[i + L * i];
			if (i < n - 1) {
				A[i + L * i] = rddot(n - i, A + i + L * i, lda, A + i + L * i, lda, rounding_mode);
				rdgemv('N', i, n - i - 1, 1.0, A + L * (i + 1), lda, A + i + L * (i + 1), lda,
					aii, A + L * i, 1, rounding_mode);
			}
			else {
				rdscal(i + 1, aii, A + L * i, 1, rounding_mode);
			}
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			const double aii = A[i + L * i];
			if (i < n - 1) {
				A[i + L * i] = rddot(n - i, A + i + L * i, 1, A + i + L * i, 1, rounding_mode);
				rdgemv('T', n - i - 1, i, 1.0, A + (i + 1), lda, A + (i + 1) + L * i, 1,
					aii, A + i, lda, rounding_mode);
			}
			else {
				rdscal(i + 1, aii, A + i, lda, rounding_mode);
			}
		}
	}
	return 0;
}

// U*U^T または L^T*L の三角部分計算 (blocked)
inline int rdlauum(const char uplo, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdlauum: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdlauum: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "LAUUM", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		rdlauu2(uplo, n, A, lda, rounding_mode);
	}
	else if (upper) {
		for (int i = 0; i < n; i += nb) {
			const int ib = std::min(nb, n - i);
			rdtrmm('R', 'U', 'T', 'N', i, ib, 1.0, A + i + L * i, lda, A + L * i, lda, rounding_mode);
			rdlauu2('U', ib, A + i + L * i, lda, rounding_mode);
			if (i + ib < n) {
				rdgemm('N', 'T', i, ib, n - i - ib, 1.0, A + L * (i + ib), lda,
					A + i + L * (i + ib), lda, 1.0, A + L * i, lda, rounding_mode);
				rdsyrk('U', 'N', ib, n - i - ib, 1.0, A + i + L * (i + ib), lda, 1.0,
					A + i + L * i, lda, rounding_mode);
			}
		}
	}
	else {
		for (int i = 0; i < n; i += nb) {
			const int ib = std::min(nb, n - i);
			rdtrmm('L', 'L', 'T', 'N', ib, i, 1.0, A + i + L * i, lda, A + i, lda, rounding_mode);
			rdlauu2('L', ib, A + i + L * i, lda, rounding_mode);
			if (i + ib < n) {
				rdgemm('T', 'N', ib, i, n - i - ib, 1.0, A + (i + ib) + L * i, lda,
					A + (i + ib), lda, 1.0, A + i, lda, rounding_mode);
				rdsyrk('L', 'T', ib, n - i - ib, 1.0, A + (i + ib) + L * i, lda, 1.0,
					A + i + L * i, lda, rounding_mode);
			}
		}
	}
	return 0;
}

// 逆行列 (A は rdgetrf の出力, ipiv 0-based)
inline int rdgetri(const int n, double* A, const int lda, const int* ipiv, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdgetri: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const int info = rdtrtri('U', 'N', n, A, lda, rounding_mode);
	if (info > 0) {
		return info;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GETRI", n, -1, -1, -1);
	const int ldwork = n;
	if (nb <= 1 || nb >= n) {
		// unblocked: L^-1 を右から解く
		std::vector<double> work(static_cast<std::size_t>(n));
		for (int j = n - 1; j >= 0; j--) {
			for (int i = j + 1; i < n; i++) {
				work[i] = A[i + L * j];
				A[i + L * j] = 0.0;
			}
			if (j < n - 1) {
				rdgemv('N', n, n - j - 1, -1.0, A + L * (j + 1), lda, work.data() + (j + 1), 1,
					1.0, A + L * j, 1, rounding_mode);
			}
		}
	}
	else {
		std::vector<double> work(static_cast<std::size_t>(ldwork) * nb);
		const int nn = ((n - 1) / nb) * nb;
		for (int j = nn; j >= 0; j -= nb) {
			const int jb = std::min(nb, n - j);
			for (int jj = j; jj < j + jb; jj++) {
				for (int i = jj + 1; i < n; i++) {
					work[i + static_cast<std::size_t>(ldwork) * (jj - j)] = A[i + L * jj];
					A[i + L * jj] = 0.0;
				}
			}
			if (j + jb < n) {
				rdgemm('N', 'N', n, jb, n - j - jb, -1.0, A + L * (j + jb), lda,
					work.data() + (j + jb), ldwork, 1.0, A + L * j, lda, rounding_mode);
			}
			rdtrsm('R', 'L', 'N', 'U', n, jb, 1.0, work.data() + j, ldwork, A + L * j, lda, rounding_mode);
		}
	}
	// 列の pivot を戻す
	for (int j = n - 2; j >= 0; j--) {
		const int jp = ipiv[j];
		if (jp != j) {
			dswap(n, A + L * j, 1, A + L * jp, 1);
		}
	}
	return 0;
}

#endif // VBLAS_RDLAPACK_LU_HPP
