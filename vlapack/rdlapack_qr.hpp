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

#ifndef VBLAS_RDLAPACK_QR_HPP
#define VBLAS_RDLAPACK_QR_HPP

// QR / LQ / QL 分解と直交変換の生成・適用，最小二乗 driver (丸めモード指定付き)．
// reference LAPACK 3.12.1 の dgeqr2/dgeqrf/dorg2r/dorgqr/dorm2r/dormqr/
// dgelq2/dgelqf/dorgl2/dorglq/dorml2/dormlq/dorg2l/dorgql/dorm2l/dormql/dgels
// を移植．WORK/LWORK は持たず内部で最適 size を確保する．

#include "rdlapack_lu.hpp"

#pragma STDC FENV_ACCESS ON

namespace vcp {

// QR 分解 (unblocked)
inline int rdgeqr2(const int m, const int n, double* A, const int lda, double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgeqr2: invalid argument");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int k = std::min(m, n);
	for (int i = 0; i < k; i++) {
		rdlarfg(m - i, A[i + L * i], A + std::min(i + 1, m - 1) + L * i, 1, tau[i], rounding_mode);
		if (i < n - 1) {
			rdlarf1f('L', m - i, n - i - 1, A + i + L * i, 1, tau[i], A + i + L * (i + 1), lda, rounding_mode);
		}
	}
	return 0;
}

// QR 分解 (blocked)
inline int rdgeqrf(const int m, const int n, double* A, const int lda, double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgeqrf: invalid argument");
	}
	const int k = std::min(m, n);
	if (k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GEQRF", m, n, -1, -1);
	const int nx = std::max(0, detail::ilaenv(3, "GEQRF", m, n, -1, -1));
	int i = 0;
	if (nb >= 2 && nb < k && nx < k) {
		std::vector<double> T(static_cast<std::size_t>(nb) * nb);
		for (; i < k - nx; i += nb) {
			const int ib = std::min(k - i, nb);
			rdgeqr2(m - i, ib, A + i + L * i, lda, tau + i, rounding_mode);
			if (i + ib < n) {
				rdlarft('F', 'C', m - i, ib, A + i + L * i, lda, tau + i, T.data(), nb, rounding_mode);
				rdlarfb('L', 'T', 'F', 'C', m - i, n - i - ib, ib, A + i + L * i, lda,
					T.data(), nb, A + i + L * (i + ib), lda, rounding_mode);
			}
		}
	}
	if (i < k) {
		rdgeqr2(m - i, n - i, A + i + L * i, lda, tau + i, rounding_mode);
	}
	return 0;
}

// QR の Q の生成 (unblocked)
inline int rdorg2r(const int m, const int n, const int k, double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::rdlapack_error("rdorg2r: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	for (int j = k; j < n; j++) {
		for (int l = 0; l < m; l++) {
			A[l + L * j] = 0.0;
		}
		A[j + L * j] = 1.0;
	}
	for (int i = k - 1; i >= 0; i--) {
		if (i < n - 1) {
			rdlarf1f('L', m - i, n - i - 1, A + i + L * i, 1, tau[i], A + i + L * (i + 1), lda, rounding_mode);
		}
		if (i < m - 1) {
			rdscal(m - i - 1, -tau[i], A + (i + 1) + L * i, 1, rounding_mode);
		}
		A[i + L * i] = 1.0 - tau[i];
		for (int l = 0; l < i; l++) {
			A[l + L * i] = 0.0;
		}
	}
	return 0;
}

// QR の Q の生成 (blocked)
inline int rdorgqr(const int m, const int n, const int k, double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::rdlapack_error("rdorgqr: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "ORGQR", m, n, k, -1);
	const int nx = std::max(0, detail::ilaenv(3, "ORGQR", m, n, k, -1));
	int kk = 0;
	int ki = 0;
	if (nb >= 2 && nb < k && nx < k) {
		ki = ((k - nx - 1) / nb) * nb;
		kk = std::min(k, ki + nb);
		for (int j = kk; j < n; j++) {
			for (int i = 0; i < kk; i++) {
				A[i + L * j] = 0.0;
			}
		}
	}
	if (kk < n) {
		rdorg2r(m - kk, n - kk, k - kk, A + kk + L * kk, lda, tau + kk, rounding_mode);
	}
	if (kk > 0) {
		std::vector<double> T(static_cast<std::size_t>(nb) * nb);
		for (int i = ki; i >= 0; i -= nb) {
			const int ib = std::min(nb, k - i);
			if (i + ib < n) {
				rdlarft('F', 'C', m - i, ib, A + i + L * i, lda, tau + i, T.data(), nb, rounding_mode);
				rdlarfb('L', 'N', 'F', 'C', m - i, n - i - ib, ib, A + i + L * i, lda,
					T.data(), nb, A + i + L * (i + ib), lda, rounding_mode);
			}
			rdorg2r(m - i, ib, ib, A + i + L * i, lda, tau + i, rounding_mode);
			for (int j = i; j < i + ib; j++) {
				for (int l = 0; l < i; l++) {
					A[l + L * j] = 0.0;
				}
			}
		}
	}
	return 0;
}

// QR の Q の適用 (unblocked): C := Q*C, Q^T*C, C*Q, C*Q^T
inline int rdorm2r(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdorm2r: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const bool forward = (left && !notran) || (!left && notran);
	const int i1 = forward ? 0 : k - 1;
	const int i2 = forward ? k : -1;
	const int i3 = forward ? 1 : -1;
	for (int i = i1; i != i2; i += i3) {
		const int mi = left ? m - i : m;
		const int ni = left ? n : n - i;
		const int ic = left ? i : 0;
		const int jc = left ? 0 : i;
		rdlarf1f(side, mi, ni, A + i + L * i, 1, tau[i], C + ic + LC * jc, ldc, rounding_mode);
	}
	return 0;
}

// QR の Q の適用 (blocked)
inline int rdormqr(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdormqr: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const int nb = std::min(64, detail::ilaenv(1, "ORMQR", m, n, k, -1));
	if (nb < 2 || nb >= k) {
		return rdorm2r(side, trans, m, n, k, A, lda, tau, C, ldc, rounding_mode);
	}
	std::vector<double> T(static_cast<std::size_t>(nb) * nb);
	const bool forward = (left && !notran) || (!left && notran);
	const int i1 = forward ? 0 : ((k - 1) / nb) * nb;
	const int i3 = forward ? nb : -nb;
	for (int i = i1; forward ? (i < k) : (i >= 0); i += i3) {
		const int ib = std::min(nb, k - i);
		rdlarft('F', 'C', nq - i, ib, A + i + L * i, lda, tau + i, T.data(), nb, rounding_mode);
		const int mi = left ? m - i : m;
		const int ni = left ? n : n - i;
		const int ic = left ? i : 0;
		const int jc = left ? 0 : i;
		rdlarfb(side, trans, 'F', 'C', mi, ni, ib, A + i + L * i, lda, T.data(), nb,
			C + ic + LC * jc, ldc, rounding_mode);
	}
	return 0;
}

// LQ 分解 (unblocked)
inline int rdgelq2(const int m, const int n, double* A, const int lda, double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgelq2: invalid argument");
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int k = std::min(m, n);
	for (int i = 0; i < k; i++) {
		rdlarfg(n - i, A[i + L * i], A + i + L * std::min(i + 1, n - 1), lda, tau[i], rounding_mode);
		if (i < m - 1) {
			rdlarf1f('R', m - i - 1, n - i, A + i + L * i, lda, tau[i], A + (i + 1) + L * i, lda, rounding_mode);
		}
	}
	return 0;
}

// LQ 分解 (blocked)
inline int rdgelqf(const int m, const int n, double* A, const int lda, double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::rdlapack_error("rdgelqf: invalid argument");
	}
	const int k = std::min(m, n);
	if (k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GELQF", m, n, -1, -1);
	const int nx = std::max(0, detail::ilaenv(3, "GELQF", m, n, -1, -1));
	int i = 0;
	if (nb >= 2 && nb < k && nx < k) {
		std::vector<double> T(static_cast<std::size_t>(nb) * nb);
		for (; i < k - nx; i += nb) {
			const int ib = std::min(k - i, nb);
			rdgelq2(ib, n - i, A + i + L * i, lda, tau + i, rounding_mode);
			if (i + ib < m) {
				rdlarft('F', 'R', n - i, ib, A + i + L * i, lda, tau + i, T.data(), nb, rounding_mode);
				rdlarfb('R', 'N', 'F', 'R', m - i - ib, n - i, ib, A + i + L * i, lda,
					T.data(), nb, A + (i + ib) + L * i, lda, rounding_mode);
			}
		}
	}
	if (i < k) {
		rdgelq2(m - i, n - i, A + i + L * i, lda, tau + i, rounding_mode);
	}
	return 0;
}

// LQ の Q の生成 (unblocked)
inline int rdorgl2(const int m, const int n, const int k, double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < m || k < 0 || k > m || lda < std::max(1, m)) {
		detail::rdlapack_error("rdorgl2: invalid argument");
	}
	if (m <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	if (k < m) {
		for (int j = 0; j < n; j++) {
			for (int l = k; l < m; l++) {
				A[l + L * j] = 0.0;
			}
			if (j >= k && j < m) {
				A[j + L * j] = 1.0;
			}
		}
	}
	for (int i = k - 1; i >= 0; i--) {
		if (i < n - 1) {
			if (i < m - 1) {
				rdlarf1f('R', m - i - 1, n - i, A + i + L * i, lda, tau[i], A + (i + 1) + L * i, lda, rounding_mode);
			}
			rdscal(n - i - 1, -tau[i], A + i + L * (i + 1), lda, rounding_mode);
		}
		A[i + L * i] = 1.0 - tau[i];
		for (int l = 0; l < i; l++) {
			A[i + L * l] = 0.0;
		}
	}
	return 0;
}

// LQ の Q の生成 (blocked)
inline int rdorglq(const int m, const int n, const int k, double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < m || k < 0 || k > m || lda < std::max(1, m)) {
		detail::rdlapack_error("rdorglq: invalid argument");
	}
	if (m <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "ORGLQ", m, n, k, -1);
	const int nx = std::max(0, detail::ilaenv(3, "ORGLQ", m, n, k, -1));
	int kk = 0;
	int ki = 0;
	if (nb >= 2 && nb < k && nx < k) {
		ki = ((k - nx - 1) / nb) * nb;
		kk = std::min(k, ki + nb);
		for (int j = 0; j < kk; j++) {
			for (int i = kk; i < m; i++) {
				A[i + L * j] = 0.0;
			}
		}
	}
	if (kk < m) {
		rdorgl2(m - kk, n - kk, k - kk, A + kk + L * kk, lda, tau + kk, rounding_mode);
	}
	if (kk > 0) {
		std::vector<double> T(static_cast<std::size_t>(nb) * nb);
		for (int i = ki; i >= 0; i -= nb) {
			const int ib = std::min(nb, k - i);
			if (i + ib < m) {
				rdlarft('F', 'R', n - i, ib, A + i + L * i, lda, tau + i, T.data(), nb, rounding_mode);
				rdlarfb('R', 'T', 'F', 'R', m - i - ib, n - i, ib, A + i + L * i, lda,
					T.data(), nb, A + (i + ib) + L * i, lda, rounding_mode);
			}
			rdorgl2(ib, n - i, ib, A + i + L * i, lda, tau + i, rounding_mode);
			for (int j = 0; j < i; j++) {
				for (int l = i; l < i + ib; l++) {
					A[l + L * j] = 0.0;
				}
			}
		}
	}
	return 0;
}

// LQ の Q の適用 (unblocked)
inline int rdorml2(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, k) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdorml2: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : k - 1;
	const int i2 = forward ? k : -1;
	const int i3 = forward ? 1 : -1;
	for (int i = i1; i != i2; i += i3) {
		const int mi = left ? m - i : m;
		const int ni = left ? n : n - i;
		const int ic = left ? i : 0;
		const int jc = left ? 0 : i;
		rdlarf1f(side, mi, ni, A + i + L * i, lda, tau[i], C + ic + LC * jc, ldc, rounding_mode);
	}
	return 0;
}

// LQ の Q の適用 (blocked)
inline int rdormlq(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, k) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdormlq: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const int nb = std::min(64, detail::ilaenv(1, "ORMLQ", m, n, k, -1));
	if (nb < 2 || nb >= k) {
		return rdorml2(side, trans, m, n, k, A, lda, tau, C, ldc, rounding_mode);
	}
	std::vector<double> T(static_cast<std::size_t>(nb) * nb);
	const char transt = notran ? 'T' : 'N';
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : ((k - 1) / nb) * nb;
	const int i3 = forward ? nb : -nb;
	for (int i = i1; forward ? (i < k) : (i >= 0); i += i3) {
		const int ib = std::min(nb, k - i);
		rdlarft('F', 'R', nq - i, ib, A + i + L * i, lda, tau + i, T.data(), nb, rounding_mode);
		const int mi = left ? m - i : m;
		const int ni = left ? n : n - i;
		const int ic = left ? i : 0;
		const int jc = left ? 0 : i;
		rdlarfb(side, transt, 'F', 'R', mi, ni, ib, A + i + L * i, lda, T.data(), nb,
			C + ic + LC * jc, ldc, rounding_mode);
	}
	return 0;
}

// QL の Q の生成 (unblocked)
inline int rdorg2l(const int m, const int n, const int k, double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::rdlapack_error("rdorg2l: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	for (int j = 0; j < n - k; j++) {
		for (int l = 0; l < m; l++) {
			A[l + L * j] = 0.0;
		}
		A[(m - n + j) + L * j] = 1.0;
	}
	for (int i = 0; i < k; i++) {
		const int ii = n - k + i; // 0-based column
		rdlarf1l('L', m - n + ii + 1, ii, A + L * ii, 1, tau[i], A, lda, rounding_mode);
		rdscal(m - n + ii, -tau[i], A + L * ii, 1, rounding_mode);
		A[(m - n + ii) + L * ii] = 1.0 - tau[i];
		for (int l = m - n + ii + 1; l < m; l++) {
			A[l + L * ii] = 0.0;
		}
	}
	return 0;
}

// QL の Q の生成 (blocked)
inline int rdorgql(const int m, const int n, const int k, double* A, const int lda, const double* tau, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::rdlapack_error("rdorgql: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "ORGQL", m, n, k, -1);
	const int nx = std::max(0, detail::ilaenv(3, "ORGQL", m, n, k, -1));
	int kk = 0;
	if (nb >= 2 && nb < k && nx < k) {
		kk = std::min(k, ((k - nx + nb - 1) / nb) * nb);
		for (int j = 0; j < n - kk; j++) {
			for (int i = m - kk; i < m; i++) {
				A[i + L * j] = 0.0;
			}
		}
	}
	rdorg2l(m - kk, n - kk, k - kk, A, lda, tau, rounding_mode);
	if (kk > 0) {
		std::vector<double> T(static_cast<std::size_t>(nb) * nb);
		for (int i = k - kk; i < k; i += nb) {
			const int ib = std::min(nb, k - i);
			const int col = n - k + i; // block の先頭列 (0-based)
			if (col > 0) {
				rdlarft('B', 'C', m - k + i + ib, ib, A + L * col, lda, tau + i, T.data(), nb, rounding_mode);
				rdlarfb('L', 'N', 'B', 'C', m - k + i + ib, col, ib, A + L * col, lda,
					T.data(), nb, A, lda, rounding_mode);
			}
			rdorg2l(m - k + i + ib, ib, ib, A + L * col, lda, tau + i, rounding_mode);
			for (int j = col; j < col + ib; j++) {
				for (int l = m - k + i + ib; l < m; l++) {
					A[l + L * j] = 0.0;
				}
			}
		}
	}
	return 0;
}

// QL の Q の適用 (unblocked)
inline int rdorm2l(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdorm2l: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : k - 1;
	const int i2 = forward ? k : -1;
	const int i3 = forward ? 1 : -1;
	for (int i = i1; i != i2; i += i3) {
		const int mi = left ? m - k + i + 1 : m;
		const int ni = left ? n : n - k + i + 1;
		rdlarf1l(side, mi, ni, A + L * i, 1, tau[i], C, ldc, rounding_mode);
	}
	return 0;
}

// QL の Q の適用 (blocked)
inline int rdormql(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::rdlapack_error("rdormql: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = std::min(64, detail::ilaenv(1, "ORMQL", m, n, k, -1));
	if (nb < 2 || nb >= k) {
		return rdorm2l(side, trans, m, n, k, A, lda, tau, C, ldc, rounding_mode);
	}
	std::vector<double> T(static_cast<std::size_t>(nb) * nb);
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : ((k - 1) / nb) * nb;
	const int i3 = forward ? nb : -nb;
	for (int i = i1; forward ? (i < k) : (i >= 0); i += i3) {
		const int ib = std::min(nb, k - i);
		rdlarft('B', 'C', nq - k + i + ib, ib, A + L * i, lda, tau + i, T.data(), nb, rounding_mode);
		const int mi = left ? m - k + i + ib : m;
		const int ni = left ? n : n - k + i + ib;
		rdlarfb(side, trans, 'B', 'C', mi, ni, ib, A + L * i, lda, T.data(), nb, C, ldc, rounding_mode);
	}
	return 0;
}

// 最小二乗問題 min ||B - A*X|| または min ||B - A^T*X|| の求解 (driver)
// B は max(m,n) x nrhs (解は先頭 n 行または m 行に上書き)
inline int rdgels(const char trans, const int m, const int n, const int nrhs,
	double* A, const int lda, double* B, const int ldb, const int rounding_mode) {
	namespace detail = vblas_rdlapack_detail;
	const bool tpsd = !detail::option_is(trans, 'N');
	if (tpsd && !detail::option_is(trans, 'T')) {
		detail::rdlapack_error("rdgels: invalid trans");
	}
	if (m < 0 || n < 0 || nrhs < 0 || lda < std::max(1, m) || ldb < std::max(1, std::max(m, n))) {
		detail::rdlapack_error("rdgels: invalid argument");
	}
	if (std::min(std::min(m, n), nrhs) == 0) {
		dlaset('F', std::max(m, n), nrhs, 0.0, 0.0, B, ldb);
		return 0;
	}
	detail::RoundingGuard guard(detail::fe_rounding(rounding_mode));
	const std::size_t LB = static_cast<std::size_t>(ldb);
	const int mn = std::min(m, n);
	std::vector<double> tau(static_cast<std::size_t>(mn));
	const double smlnum = DBL_MIN / DBL_EPSILON;
	const double bignum = 1.0 / smlnum;
	const double anrm = rdlange('M', m, n, A, lda, rounding_mode);
	int iascl = 0;
	if (anrm > 0.0 && anrm < smlnum) {
		rdlascl('G', 0, 0, anrm, smlnum, m, n, A, lda, rounding_mode);
		iascl = 1;
	}
	else if (anrm > bignum) {
		rdlascl('G', 0, 0, anrm, bignum, m, n, A, lda, rounding_mode);
		iascl = 2;
	}
	else if (anrm == 0.0) {
		// A = 0: 最小二乗解は X = 0
		dlaset('F', std::max(m, n), nrhs, 0.0, 0.0, B, ldb);
		return 0;
	}
	const int brow = tpsd ? n : m;
	const double bnrm = rdlange('M', brow, nrhs, B, ldb, rounding_mode);
	int ibscl = 0;
	if (bnrm > 0.0 && bnrm < smlnum) {
		rdlascl('G', 0, 0, bnrm, smlnum, brow, nrhs, B, ldb, rounding_mode);
		ibscl = 1;
	}
	else if (bnrm > bignum) {
		rdlascl('G', 0, 0, bnrm, bignum, brow, nrhs, B, ldb, rounding_mode);
		ibscl = 2;
	}
	int scllen = 0;
	if (m >= n) {
		rdgeqrf(m, n, A, lda, tau.data(), rounding_mode);
		if (!tpsd) {
			// 最小二乗解: B := inv(R) * Q^T * B
			rdormqr('L', 'T', m, nrhs, n, A, lda, tau.data(), B, ldb, rounding_mode);
			const int info = rdtrtrs('U', 'N', 'N', n, nrhs, A, lda, B, ldb, rounding_mode);
			if (info > 0) {
				return info;
			}
			scllen = n;
		}
		else {
			// underdetermined A^T * X = B の最小 norm 解
			const int info = rdtrtrs('U', 'T', 'N', n, nrhs, A, lda, B, ldb, rounding_mode);
			if (info > 0) {
				return info;
			}
			for (int j = 0; j < nrhs; j++) {
				for (int i = n; i < m; i++) {
					B[i + LB * j] = 0.0;
				}
			}
			rdormqr('L', 'N', m, nrhs, n, A, lda, tau.data(), B, ldb, rounding_mode);
			scllen = m;
		}
	}
	else {
		rdgelqf(m, n, A, lda, tau.data(), rounding_mode);
		if (!tpsd) {
			// underdetermined A * X = B の最小 norm 解
			const int info = rdtrtrs('L', 'N', 'N', m, nrhs, A, lda, B, ldb, rounding_mode);
			if (info > 0) {
				return info;
			}
			for (int j = 0; j < nrhs; j++) {
				for (int i = m; i < n; i++) {
					B[i + LB * j] = 0.0;
				}
			}
			rdormlq('L', 'T', n, nrhs, m, A, lda, tau.data(), B, ldb, rounding_mode);
			scllen = n;
		}
		else {
			// overdetermined A^T * X = B の最小二乗解
			rdormlq('L', 'N', n, nrhs, m, A, lda, tau.data(), B, ldb, rounding_mode);
			const int info = rdtrtrs('L', 'T', 'N', m, nrhs, A, lda, B, ldb, rounding_mode);
			if (info > 0) {
				return info;
			}
			scllen = m;
		}
	}
	if (iascl == 1) {
		rdlascl('G', 0, 0, anrm, smlnum, scllen, nrhs, B, ldb, rounding_mode);
	}
	else if (iascl == 2) {
		rdlascl('G', 0, 0, anrm, bignum, scllen, nrhs, B, ldb, rounding_mode);
	}
	if (ibscl == 1) {
		rdlascl('G', 0, 0, smlnum, bnrm, scllen, nrhs, B, ldb, rounding_mode);
	}
	else if (ibscl == 2) {
		rdlascl('G', 0, 0, bignum, bnrm, scllen, nrhs, B, ldb, rounding_mode);
	}
	return 0;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_QR_HPP
