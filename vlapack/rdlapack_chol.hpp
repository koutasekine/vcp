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

#ifndef VBLAS_RDLAPACK_CHOL_HPP
#define VBLAS_RDLAPACK_CHOL_HPP

// Cholesky 分解 routine (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dpotf2/dpotrf2/dpotrf/dpotrs/dposv/dpotri を移植．

#include "rdlapack_lu.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// Cholesky 分解 (unblocked)
inline int rdpotf2(const char uplo, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpotf2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdpotf2: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int j = 0; j < n; j++) {
			double ajj = A[j + L * j] - rddot(j, A + L * j, 1, A + L * j, 1, rounding_mode);
			if (ajj <= 0.0 || std::isnan(ajj)) {
				A[j + L * j] = ajj;
				return j + 1;
			}
			ajj = std::sqrt(ajj);
			A[j + L * j] = ajj;
			if (j < n - 1) {
				rdgemv('T', j, n - j - 1, -1.0, A + L * (j + 1), lda, A + L * j, 1,
					1.0, A + j + L * (j + 1), lda, rounding_mode);
				rdscal(n - j - 1, 1.0 / ajj, A + j + L * (j + 1), lda, rounding_mode);
			}
		}
	}
	else {
		for (int j = 0; j < n; j++) {
			double ajj = A[j + L * j] - rddot(j, A + j, lda, A + j, lda, rounding_mode);
			if (ajj <= 0.0 || std::isnan(ajj)) {
				A[j + L * j] = ajj;
				return j + 1;
			}
			ajj = std::sqrt(ajj);
			A[j + L * j] = ajj;
			if (j < n - 1) {
				rdgemv('N', n - j - 1, j, -1.0, A + (j + 1), lda, A + j, lda,
					1.0, A + (j + 1) + L * j, 1, rounding_mode);
				rdscal(n - j - 1, 1.0 / ajj, A + (j + 1) + L * j, 1, rounding_mode);
			}
		}
	}
	return 0;
}

// Cholesky 分解 (recursive)
inline int rdpotrf2(const char uplo, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpotrf2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdpotrf2: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (n == 1) {
		if (A[0] <= 0.0 || std::isnan(A[0])) {
			return 1;
		}
		A[0] = std::sqrt(A[0]);
		return 0;
	}
	const int n1 = n / 2;
	const int n2 = n - n1;
	int iinfo = rdpotrf2(uplo, n1, A, lda, rounding_mode);
	if (iinfo != 0) {
		return iinfo;
	}
	if (upper) {
		rdtrsm('L', 'U', 'T', 'N', n1, n2, 1.0, A, lda, A + L * n1, lda, rounding_mode);
		rdsyrk('U', 'T', n2, n1, -1.0, A + L * n1, lda, 1.0, A + n1 + L * n1, lda, rounding_mode);
	}
	else {
		rdtrsm('R', 'L', 'T', 'N', n2, n1, 1.0, A, lda, A + n1, lda, rounding_mode);
		rdsyrk('L', 'N', n2, n1, -1.0, A + n1, lda, 1.0, A + n1 + L * n1, lda, rounding_mode);
	}
	iinfo = rdpotrf2(uplo, n2, A + n1 + L * n1, lda, rounding_mode);
	if (iinfo != 0) {
		return iinfo + n1;
	}
	return 0;
}

// Cholesky 分解 (blocked)
inline int rdpotrf(const char uplo, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpotrf: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdpotrf: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "POTRF", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		return rdpotrf2(uplo, n, A, lda, rounding_mode);
	}
	if (upper) {
		for (int j = 0; j < n; j += nb) {
			const int jb = std::min(nb, n - j);
			rdsyrk('U', 'T', jb, j, -1.0, A + L * j, lda, 1.0, A + j + L * j, lda, rounding_mode);
			const int iinfo = rdpotrf2('U', jb, A + j + L * j, lda, rounding_mode);
			if (iinfo != 0) {
				return iinfo + j;
			}
			if (j + jb < n) {
				rdgemm('T', 'N', jb, n - j - jb, j, -1.0, A + L * j, lda, A + L * (j + jb), lda,
					1.0, A + j + L * (j + jb), lda, rounding_mode);
				rdtrsm('L', 'U', 'T', 'N', jb, n - j - jb, 1.0, A + j + L * j, lda,
					A + j + L * (j + jb), lda, rounding_mode);
			}
		}
	}
	else {
		for (int j = 0; j < n; j += nb) {
			const int jb = std::min(nb, n - j);
			rdsyrk('L', 'N', jb, j, -1.0, A + j, lda, 1.0, A + j + L * j, lda, rounding_mode);
			const int iinfo = rdpotrf2('L', jb, A + j + L * j, lda, rounding_mode);
			if (iinfo != 0) {
				return iinfo + j;
			}
			if (j + jb < n) {
				rdgemm('N', 'T', n - j - jb, jb, j, -1.0, A + (j + jb), lda, A + j, lda,
					1.0, A + (j + jb) + L * j, lda, rounding_mode);
				rdtrsm('R', 'L', 'T', 'N', n - j - jb, jb, 1.0, A + j + L * j, lda,
					A + (j + jb) + L * j, lda, rounding_mode);
			}
		}
	}
	return 0;
}

// Cholesky 分解を使った連立一次方程式の求解 (A は rdpotrf の出力)
inline int rdpotrs(const char uplo, const int n, const int nrhs, const double* A, const int lda,
	double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpotrs: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdpotrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	if (upper) {
		rdtrsm('L', 'U', 'T', 'N', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
		rdtrsm('L', 'U', 'N', 'N', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
	}
	else {
		rdtrsm('L', 'L', 'N', 'N', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
		rdtrsm('L', 'L', 'T', 'N', n, nrhs, 1.0, A, lda, B, ldb, rounding_mode);
	}
	return 0;
}

// 対称正定値の連立一次方程式 A*X = B の求解 (driver)
inline int rdposv(const char uplo, const int n, const int nrhs, double* A, const int lda,
	double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdposv: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::rdlapack_error("rdposv: invalid argument");
	}
	const int info = rdpotrf(uplo, n, A, lda, rounding_mode);
	if (info == 0) {
		rdpotrs(uplo, n, nrhs, A, lda, B, ldb, rounding_mode);
	}
	return info;
}

// 対称正定値行列の逆行列 (A は rdpotrf の出力)
inline int rdpotri(const char uplo, const int n, double* A, const int lda, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::rdlapack_error("rdpotri: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::rdlapack_error("rdpotri: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const int info = rdtrtri(uplo, 'N', n, A, lda, rounding_mode);
	if (info > 0) {
		return info;
	}
	rdlauum(uplo, n, A, lda, rounding_mode);
	return 0;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_CHOL_HPP
