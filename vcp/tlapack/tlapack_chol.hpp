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

#ifndef TLAPACK_TLAPACK_CHOL_HPP
#define TLAPACK_TLAPACK_CHOL_HPP

// Cholesky 分解 routine (template 版)．
// reference LAPACK 3.12.1 の dpotf2/dpotrf2/dpotrf/dpotrs/dposv/dpotri を移植．

#include "tlapack_lu.hpp"


// Cholesky 分解 (unblocked)
template <typename T>
inline int tpotf2(const char uplo, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotf2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tpotf2: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int j = 0; j < n; j++) {
			T ajj = A[j + L * j] - tdot(j, A + L * j, 1, A + L * j, 1);
			if (ajj <= T(0) || tlapack_detail::tisnan(ajj)) {
				A[j + L * j] = ajj;
				return j + 1;
			}
			ajj = tlapack_detail::tsqrt(ajj);
			A[j + L * j] = ajj;
			if (j < n - 1) {
				tgemv('T', j, n - j - 1, -T(1), A + L * (j + 1), lda, A + L * j, 1,
					T(1), A + j + L * (j + 1), lda);
				tscal(n - j - 1, T(1) / ajj, A + j + L * (j + 1), lda);
			}
		}
	}
	else {
		for (int j = 0; j < n; j++) {
			T ajj = A[j + L * j] - tdot(j, A + j, lda, A + j, lda);
			if (ajj <= T(0) || tlapack_detail::tisnan(ajj)) {
				A[j + L * j] = ajj;
				return j + 1;
			}
			ajj = tlapack_detail::tsqrt(ajj);
			A[j + L * j] = ajj;
			if (j < n - 1) {
				tgemv('N', n - j - 1, j, -T(1), A + (j + 1), lda, A + j, lda,
					T(1), A + (j + 1) + L * j, 1);
				tscal(n - j - 1, T(1) / ajj, A + (j + 1) + L * j, 1);
			}
		}
	}
	return 0;
}

// Cholesky 分解 (recursive)
template <typename T>
inline int tpotrf2(const char uplo, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotrf2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tpotrf2: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (n == 1) {
		if (A[0] <= T(0) || tlapack_detail::tisnan(A[0])) {
			return 1;
		}
		A[0] = tlapack_detail::tsqrt(A[0]);
		return 0;
	}
	const int n1 = n / 2;
	const int n2 = n - n1;
	int iinfo = tpotrf2(uplo, n1, A, lda);
	if (iinfo != 0) {
		return iinfo;
	}
	if (upper) {
		ttrsm('L', 'U', 'T', 'N', n1, n2, T(1), A, lda, A + L * n1, lda);
		tsyrk('U', 'T', n2, n1, -T(1), A + L * n1, lda, T(1), A + n1 + L * n1, lda);
	}
	else {
		ttrsm('R', 'L', 'T', 'N', n2, n1, T(1), A, lda, A + n1, lda);
		tsyrk('L', 'N', n2, n1, -T(1), A + n1, lda, T(1), A + n1 + L * n1, lda);
	}
	iinfo = tpotrf2(uplo, n2, A + n1 + L * n1, lda);
	if (iinfo != 0) {
		return iinfo + n1;
	}
	return 0;
}

// Cholesky 分解 (blocked)
template <typename T>
inline int tpotrf(const char uplo, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotrf: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tpotrf: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "POTRF", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		return tpotrf2(uplo, n, A, lda);
	}
	if (upper) {
		for (int j = 0; j < n; j += nb) {
			const int jb = std::min(nb, n - j);
			tsyrk('U', 'T', jb, j, -T(1), A + L * j, lda, T(1), A + j + L * j, lda);
			const int iinfo = tpotrf2('U', jb, A + j + L * j, lda);
			if (iinfo != 0) {
				return iinfo + j;
			}
			if (j + jb < n) {
				tgemm('T', 'N', jb, n - j - jb, j, -T(1), A + L * j, lda, A + L * (j + jb), lda,
					T(1), A + j + L * (j + jb), lda);
				ttrsm('L', 'U', 'T', 'N', jb, n - j - jb, T(1), A + j + L * j, lda,
					A + j + L * (j + jb), lda);
			}
		}
	}
	else {
		for (int j = 0; j < n; j += nb) {
			const int jb = std::min(nb, n - j);
			tsyrk('L', 'N', jb, j, -T(1), A + j, lda, T(1), A + j + L * j, lda);
			const int iinfo = tpotrf2('L', jb, A + j + L * j, lda);
			if (iinfo != 0) {
				return iinfo + j;
			}
			if (j + jb < n) {
				tgemm('N', 'T', n - j - jb, jb, j, -T(1), A + (j + jb), lda, A + j, lda,
					T(1), A + (j + jb) + L * j, lda);
				ttrsm('R', 'L', 'T', 'N', n - j - jb, jb, T(1), A + j + L * j, lda,
					A + (j + jb) + L * j, lda);
			}
		}
	}
	return 0;
}

// Cholesky 分解を使った連立一次方程式の求解 (A は tpotrf の出力)
template <typename T>
inline int tpotrs(const char uplo, const int n, const int nrhs, const T* A, const int lda,
	T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotrs: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tpotrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	if (upper) {
		ttrsm('L', 'U', 'T', 'N', n, nrhs, T(1), A, lda, B, ldb);
		ttrsm('L', 'U', 'N', 'N', n, nrhs, T(1), A, lda, B, ldb);
	}
	else {
		ttrsm('L', 'L', 'N', 'N', n, nrhs, T(1), A, lda, B, ldb);
		ttrsm('L', 'L', 'T', 'N', n, nrhs, T(1), A, lda, B, ldb);
	}
	return 0;
}

// 対称正定値の連立一次方程式 A*X = B の求解 (driver)
template <typename T>
inline int tposv(const char uplo, const int n, const int nrhs, T* A, const int lda,
	T* B, const int ldb) {
	namespace detail = tlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tposv: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tposv: invalid argument");
	}
	const int info = tpotrf(uplo, n, A, lda);
	if (info == 0) {
		tpotrs(uplo, n, nrhs, A, lda, B, ldb);
	}
	return info;
}

// 対称正定値行列の逆行列 (A は tpotrf の出力)
template <typename T>
inline int tpotri(const char uplo, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotri: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tpotri: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const int info = ttrtri(uplo, 'N', n, A, lda);
	if (info > 0) {
		return info;
	}
	tlauum(uplo, n, A, lda);
	return 0;
}

#endif // TLAPACK_TLAPACK_CHOL_HPP
