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

#ifndef TLAPACK_TLAPACK_LU_HPP
#define TLAPACK_TLAPACK_LU_HPP

// LU 分解と三角行列 routine (template 版)．
// reference LAPACK 3.12.1 の dgetf2/dgetrf2/dgetrf/dgetrs/dgesv/dgetri/
// dtrti2/dtrtri/dtrtrs/dlauu2/dlauum を移植．
// ipiv は 0-based．返り値は LAPACK の INFO (0: 成功, >0: LAPACK と同じ意味の
// 1-based 位置)．引数 error は std::invalid_argument を投げる．

#include "tlapack_aux.hpp"


// LU 分解 (部分 pivot 付き, unblocked)
template <typename T>
inline int tgetf2(const int m, const int n, T* A, const int lda, int* ipiv) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgetf2: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const T sfmin = tlamch<T>('S');
	int info = 0;
	for (int j = 0; j < std::min(m, n); j++) {
		const int jp = j + itamax(m - j, A + j + L * j, 1);
		ipiv[j] = jp;
		if (A[jp + L * j] != T(0)) {
			if (jp != j) {
				tswap(n, A + j, lda, A + jp, lda);
			}
			if (j < m - 1) {
				if (tlapack_detail::tabs(A[j + L * j]) >= sfmin) {
					tscal(m - j - 1, T(1) / A[j + L * j], A + (j + 1) + L * j, 1);
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
			tger(m - j - 1, n - j - 1, -T(1), A + (j + 1) + L * j, 1, A + j + L * (j + 1), lda,
				A + (j + 1) + L * (j + 1), lda);
		}
	}
	return info;
}

// LU 分解 (部分 pivot 付き, recursive)
template <typename T>
inline int tgetrf2(const int m, const int n, T* A, const int lda, int* ipiv) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgetrf2: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	int info = 0;
	if (m == 1) {
		ipiv[0] = 0;
		if (A[0] == T(0)) {
			info = 1;
		}
	}
	else if (n == 1) {
		const T sfmin = tlamch<T>('S');
		const int i = itamax(m, A, 1);
		ipiv[0] = i;
		if (A[i] != T(0)) {
			if (i != 0) {
				const T temp = A[0];
				A[0] = A[i];
				A[i] = temp;
			}
			if (tlapack_detail::tabs(A[0]) >= sfmin) {
				tscal(m - 1, T(1) / A[0], A + 1, 1);
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
		int iinfo = tgetrf2(m, n1, A, lda, ipiv);
		if (info == 0 && iinfo > 0) {
			info = iinfo;
		}
		tlaswp(n2, A + L * n1, lda, 0, n1, ipiv, 1);
		ttrsm('L', 'L', 'N', 'U', n1, n2, T(1), A, lda, A + L * n1, lda);
		tgemm('N', 'N', m - n1, n2, n1, -T(1), A + n1, lda, A + L * n1, lda, T(1),
			A + n1 + L * n1, lda);
		iinfo = tgetrf2(m - n1, n2, A + n1 + L * n1, lda, ipiv + n1);
		if (info == 0 && iinfo > 0) {
			info = iinfo + n1;
		}
		for (int i = n1; i < std::min(m, n); i++) {
			ipiv[i] = ipiv[i] + n1;
		}
		tlaswp(n1, A, lda, n1, std::min(m, n), ipiv, 1);
	}
	return info;
}

// LU 分解 (部分 pivot 付き, blocked)
template <typename T>
inline int tgetrf(const int m, const int n, T* A, const int lda, int* ipiv) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgetrf: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GETRF", m, n, -1, -1);
	int info = 0;
	if (nb <= 1 || nb >= std::min(m, n)) {
		return tgetrf2(m, n, A, lda, ipiv);
	}
	for (int j = 0; j < std::min(m, n); j += nb) {
		const int jb = std::min(std::min(m, n) - j, nb);
		const int iinfo = tgetrf2(m - j, jb, A + j + L * j, lda, ipiv + j);
		if (info == 0 && iinfo > 0) {
			info = iinfo + j;
		}
		for (int i = j; i < std::min(m, j + jb); i++) {
			ipiv[i] = j + ipiv[i];
		}
		tlaswp(j, A, lda, j, j + jb, ipiv, 1);
		if (j + jb < n) {
			tlaswp(n - j - jb, A + L * (j + jb), lda, j, j + jb, ipiv, 1);
			ttrsm('L', 'L', 'N', 'U', jb, n - j - jb, T(1), A + j + L * j, lda,
				A + j + L * (j + jb), lda);
			if (j + jb < m) {
				tgemm('N', 'N', m - j - jb, n - j - jb, jb, -T(1), A + (j + jb) + L * j, lda,
					A + j + L * (j + jb), lda, T(1), A + (j + jb) + L * (j + jb), lda);
			}
		}
	}
	return info;
}

// LU 分解を使った連立一次方程式の求解 (A は tgetrf の出力)
template <typename T>
inline int tgetrs(const char trans, const int n, const int nrhs, const T* A, const int lda,
	const int* ipiv, T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool notran = detail::option_is(trans, 'N');
	if (!notran && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::tlapack_error("tgetrs: invalid trans");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tgetrs: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	if (notran) {
		tlaswp(nrhs, B, ldb, 0, n, ipiv, 1);
		ttrsm('L', 'L', 'N', 'U', n, nrhs, T(1), A, lda, B, ldb);
		ttrsm('L', 'U', 'N', 'N', n, nrhs, T(1), A, lda, B, ldb);
	}
	else {
		ttrsm('L', 'U', 'T', 'N', n, nrhs, T(1), A, lda, B, ldb);
		ttrsm('L', 'L', 'T', 'U', n, nrhs, T(1), A, lda, B, ldb);
		tlaswp(nrhs, B, ldb, 0, n, ipiv, -1);
	}
	return 0;
}

// 連立一次方程式 A*X = B の求解 (driver)
template <typename T>
inline int tgesv(const int n, const int nrhs, T* A, const int lda, int* ipiv,
	T* B, const int ldb) {
	namespace detail = tlapack_detail;
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tgesv: invalid argument");
	}
	const int info = tgetrf(n, n, A, lda, ipiv);
	if (info == 0) {
		tgetrs('N', n, nrhs, A, lda, ipiv, B, ldb);
	}
	return info;
}

// 三角行列の逆行列 (unblocked)
template <typename T>
inline int ttrti2(const char uplo, const char diag, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	const bool nounit = detail::option_is(diag, 'N');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("ttrti2: invalid uplo");
	}
	if (!nounit && !detail::option_is(diag, 'U')) {
		detail::tlapack_error("ttrti2: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("ttrti2: invalid n/lda");
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int j = 0; j < n; j++) {
			T ajj;
			if (nounit) {
				A[j + L * j] = T(1) / A[j + L * j];
				ajj = -A[j + L * j];
			}
			else {
				ajj = -T(1);
			}
			ttrmv('U', 'N', diag, j, A, lda, A + L * j, 1);
			tscal(j, ajj, A + L * j, 1);
		}
	}
	else {
		for (int j = n - 1; j >= 0; j--) {
			T ajj;
			if (nounit) {
				A[j + L * j] = T(1) / A[j + L * j];
				ajj = -A[j + L * j];
			}
			else {
				ajj = -T(1);
			}
			if (j < n - 1) {
				ttrmv('L', 'N', diag, n - j - 1, A + (j + 1) + L * (j + 1), lda,
					A + (j + 1) + L * j, 1);
				tscal(n - j - 1, ajj, A + (j + 1) + L * j, 1);
			}
		}
	}
	return 0;
}

// 三角行列の逆行列 (blocked)
template <typename T>
inline int ttrtri(const char uplo, const char diag, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	const bool nounit = detail::option_is(diag, 'N');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("ttrtri: invalid uplo");
	}
	if (!nounit && !detail::option_is(diag, 'U')) {
		detail::tlapack_error("ttrtri: invalid diag");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("ttrtri: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (nounit) {
		for (int i = 0; i < n; i++) {
			if (A[i + L * i] == T(0)) {
				return i + 1;
			}
		}
	}
	const int nb = detail::ilaenv(1, "TRTRI", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		ttrti2(uplo, diag, n, A, lda);
	}
	else if (upper) {
		for (int j = 0; j < n; j += nb) {
			const int jb = std::min(nb, n - j);
			ttrmm('L', 'U', 'N', diag, j, jb, T(1), A, lda, A + L * j, lda);
			ttrsm('R', 'U', 'N', diag, j, jb, -T(1), A + j + L * j, lda, A + L * j, lda);
			ttrti2('U', diag, jb, A + j + L * j, lda);
		}
	}
	else {
		const int nn = ((n - 1) / nb) * nb;
		for (int j = nn; j >= 0; j -= nb) {
			const int jb = std::min(nb, n - j);
			if (j + jb < n) {
				ttrmm('L', 'L', 'N', diag, n - j - jb, jb, T(1), A + (j + jb) + L * (j + jb), lda,
					A + (j + jb) + L * j, lda);
				ttrsm('R', 'L', 'N', diag, n - j - jb, jb, -T(1), A + j + L * j, lda,
					A + (j + jb) + L * j, lda);
			}
			ttrti2('L', diag, jb, A + j + L * j, lda);
		}
	}
	return 0;
}

// 三角行列の連立一次方程式 (特異性 check 付き)
template <typename T>
inline int ttrtrs(const char uplo, const char trans, const char diag, const int n, const int nrhs,
	const T* A, const int lda, T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool nounit = detail::option_is(diag, 'N');
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("ttrtrs: invalid uplo");
	}
	if (!detail::option_is(trans, 'N') && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::tlapack_error("ttrtrs: invalid trans");
	}
	if (!nounit && !detail::option_is(diag, 'U')) {
		detail::tlapack_error("ttrtrs: invalid diag");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("ttrtrs: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	if (nounit) {
		for (int i = 0; i < n; i++) {
			if (A[i + static_cast<std::size_t>(lda) * i] == T(0)) {
				return i + 1;
			}
		}
	}
	ttrsm('L', uplo, trans, diag, n, nrhs, T(1), A, lda, B, ldb);
	return 0;
}

// U*U^T (uplo='U') または L^T*L (uplo='L') の三角部分計算 (unblocked)
template <typename T>
inline int tlauu2(const char uplo, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tlauu2: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tlauu2: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (upper) {
		for (int i = 0; i < n; i++) {
			const T aii = A[i + L * i];
			if (i < n - 1) {
				A[i + L * i] = tdot(n - i, A + i + L * i, lda, A + i + L * i, lda);
				tgemv('N', i, n - i - 1, T(1), A + L * (i + 1), lda, A + i + L * (i + 1), lda,
					aii, A + L * i, 1);
			}
			else {
				tscal(i + 1, aii, A + L * i, 1);
			}
		}
	}
	else {
		for (int i = 0; i < n; i++) {
			const T aii = A[i + L * i];
			if (i < n - 1) {
				A[i + L * i] = tdot(n - i, A + i + L * i, 1, A + i + L * i, 1);
				tgemv('T', n - i - 1, i, T(1), A + (i + 1), lda, A + (i + 1) + L * i, 1,
					aii, A + i, lda);
			}
			else {
				tscal(i + 1, aii, A + i, lda);
			}
		}
	}
	return 0;
}

// U*U^T または L^T*L の三角部分計算 (blocked)
template <typename T>
inline int tlauum(const char uplo, const int n, T* A, const int lda) {
	namespace detail = tlapack_detail;
	const bool upper = detail::option_is(uplo, 'U');
	if (!upper && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tlauum: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tlauum: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "LAUUM", n, -1, -1, -1);
	if (nb <= 1 || nb >= n) {
		tlauu2(uplo, n, A, lda);
	}
	else if (upper) {
		for (int i = 0; i < n; i += nb) {
			const int ib = std::min(nb, n - i);
			ttrmm('R', 'U', 'T', 'N', i, ib, T(1), A + i + L * i, lda, A + L * i, lda);
			tlauu2('U', ib, A + i + L * i, lda);
			if (i + ib < n) {
				tgemm('N', 'T', i, ib, n - i - ib, T(1), A + L * (i + ib), lda,
					A + i + L * (i + ib), lda, T(1), A + L * i, lda);
				tsyrk('U', 'N', ib, n - i - ib, T(1), A + i + L * (i + ib), lda, T(1),
					A + i + L * i, lda);
			}
		}
	}
	else {
		for (int i = 0; i < n; i += nb) {
			const int ib = std::min(nb, n - i);
			ttrmm('L', 'L', 'T', 'N', ib, i, T(1), A + i + L * i, lda, A + i, lda);
			tlauu2('L', ib, A + i + L * i, lda);
			if (i + ib < n) {
				tgemm('T', 'N', ib, i, n - i - ib, T(1), A + (i + ib) + L * i, lda,
					A + (i + ib), lda, T(1), A + i, lda);
				tsyrk('L', 'T', ib, n - i - ib, T(1), A + (i + ib) + L * i, lda, T(1),
					A + i + L * i, lda);
			}
		}
	}
	return 0;
}

// 逆行列 (A は tgetrf の出力, ipiv 0-based)
template <typename T>
inline int tgetri(const int n, T* A, const int lda, const int* ipiv) {
	namespace detail = tlapack_detail;
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tgetri: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const int info = ttrtri('U', 'N', n, A, lda);
	if (info > 0) {
		return info;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GETRI", n, -1, -1, -1);
	const int ldwork = n;
	if (nb <= 1 || nb >= n) {
		// unblocked: L^-1 を右から解く
		std::vector<T> work(static_cast<std::size_t>(n), T(0));
		for (int j = n - 1; j >= 0; j--) {
			for (int i = j + 1; i < n; i++) {
				work[i] = A[i + L * j];
				A[i + L * j] = T(0);
			}
			if (j < n - 1) {
				tgemv('N', n, n - j - 1, -T(1), A + L * (j + 1), lda, work.data() + (j + 1), 1,
					T(1), A + L * j, 1);
			}
		}
	}
	else {
		std::vector<T> work(static_cast<std::size_t>(ldwork) * nb, T(0));
		const int nn = ((n - 1) / nb) * nb;
		for (int j = nn; j >= 0; j -= nb) {
			const int jb = std::min(nb, n - j);
			for (int jj = j; jj < j + jb; jj++) {
				for (int i = jj + 1; i < n; i++) {
					work[i + static_cast<std::size_t>(ldwork) * (jj - j)] = A[i + L * jj];
					A[i + L * jj] = T(0);
				}
			}
			if (j + jb < n) {
				tgemm('N', 'N', n, jb, n - j - jb, -T(1), A + L * (j + jb), lda,
					work.data() + (j + jb), ldwork, T(1), A + L * j, lda);
			}
			ttrsm('R', 'L', 'N', 'U', n, jb, T(1), work.data() + j, ldwork, A + L * j, lda);
		}
	}
	// 列の pivot を戻す
	for (int j = n - 2; j >= 0; j--) {
		const int jp = ipiv[j];
		if (jp != j) {
			tswap(n, A + L * j, 1, A + L * jp, 1);
		}
	}
	return 0;
}

#endif // TLAPACK_TLAPACK_LU_HPP
