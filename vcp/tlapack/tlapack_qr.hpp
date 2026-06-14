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

#ifndef TLAPACK_TLAPACK_QR_HPP
#define TLAPACK_TLAPACK_QR_HPP

// QR / LQ / QL 分解と直交変換の生成・適用，最小二乗 driver (template 版)．
// reference LAPACK 3.12.1 の dgeqr2/dgeqrf/dorg2r/dorgqr/dorm2r/dormqr/
// dgelq2/dgelqf/dorgl2/dorglq/dorml2/dormlq/dorg2l/dorgql/dorm2l/dormql/dgels
// を移植．WORK/LWORK は持たず内部で最適 size を確保する．

#include "tlapack_lu.hpp"

namespace vcp {

// QR 分解 (unblocked)
template <typename T>
inline int tgeqr2(const int m, const int n, T* A, const int lda, T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgeqr2: invalid argument");
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int k = std::min(m, n);
	for (int i = 0; i < k; i++) {
		tlarfg(m - i, A[i + L * i], A + std::min(i + 1, m - 1) + L * i, 1, tau[i]);
		if (i < n - 1) {
			tlarf1f('L', m - i, n - i - 1, A + i + L * i, 1, tau[i], A + i + L * (i + 1), lda);
		}
	}
	return 0;
}

// QR 分解 (blocked)
template <typename T>
inline int tgeqrf(const int m, const int n, T* A, const int lda, T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgeqrf: invalid argument");
	}
	const int k = std::min(m, n);
	if (k == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GEQRF", m, n, -1, -1);
	const int nx = std::max(0, detail::ilaenv(3, "GEQRF", m, n, -1, -1));
	int i = 0;
	if (nb >= 2 && nb < k && nx < k) {
		std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
		for (; i < k - nx; i += nb) {
			const int ib = std::min(k - i, nb);
			tgeqr2(m - i, ib, A + i + L * i, lda, tau + i);
			if (i + ib < n) {
				tlarft('F', 'C', m - i, ib, A + i + L * i, lda, tau + i, Tf.data(), nb);
				tlarfb('L', 'T', 'F', 'C', m - i, n - i - ib, ib, A + i + L * i, lda,
					Tf.data(), nb, A + i + L * (i + ib), lda);
			}
		}
	}
	if (i < k) {
		tgeqr2(m - i, n - i, A + i + L * i, lda, tau + i);
	}
	return 0;
}

// QR の Q の生成 (unblocked)
template <typename T>
inline int torg2r(const int m, const int n, const int k, T* A, const int lda, const T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::tlapack_error("torg2r: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	for (int j = k; j < n; j++) {
		for (int l = 0; l < m; l++) {
			A[l + L * j] = T(0);
		}
		A[j + L * j] = T(1);
	}
	for (int i = k - 1; i >= 0; i--) {
		if (i < n - 1) {
			tlarf1f('L', m - i, n - i - 1, A + i + L * i, 1, tau[i], A + i + L * (i + 1), lda);
		}
		if (i < m - 1) {
			tscal(m - i - 1, -tau[i], A + (i + 1) + L * i, 1);
		}
		A[i + L * i] = T(1) - tau[i];
		for (int l = 0; l < i; l++) {
			A[l + L * i] = T(0);
		}
	}
	return 0;
}

// QR の Q の生成 (blocked)
template <typename T>
inline int torgqr(const int m, const int n, const int k, T* A, const int lda, const T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::tlapack_error("torgqr: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
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
				A[i + L * j] = T(0);
			}
		}
	}
	if (kk < n) {
		torg2r(m - kk, n - kk, k - kk, A + kk + L * kk, lda, tau + kk);
	}
	if (kk > 0) {
		std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
		for (int i = ki; i >= 0; i -= nb) {
			const int ib = std::min(nb, k - i);
			if (i + ib < n) {
				tlarft('F', 'C', m - i, ib, A + i + L * i, lda, tau + i, Tf.data(), nb);
				tlarfb('L', 'N', 'F', 'C', m - i, n - i - ib, ib, A + i + L * i, lda,
					Tf.data(), nb, A + i + L * (i + ib), lda);
			}
			torg2r(m - i, ib, ib, A + i + L * i, lda, tau + i);
			for (int j = i; j < i + ib; j++) {
				for (int l = 0; l < i; l++) {
					A[l + L * j] = T(0);
				}
			}
		}
	}
	return 0;
}

// QR の Q の適用 (unblocked): C := Q*C, Q^T*C, C*Q, C*Q^T
template <typename T>
inline int torm2r(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::tlapack_error("torm2r: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
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
		tlarf1f(side, mi, ni, A + i + L * i, 1, tau[i], C + ic + LC * jc, ldc);
	}
	return 0;
}

// QR の Q の適用 (blocked)
template <typename T>
inline int tormqr(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::tlapack_error("tormqr: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const int nb = std::min(64, detail::ilaenv(1, "ORMQR", m, n, k, -1));
	if (nb < 2 || nb >= k) {
		return torm2r(side, trans, m, n, k, A, lda, tau, C, ldc);
	}
	std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
	const bool forward = (left && !notran) || (!left && notran);
	const int i1 = forward ? 0 : ((k - 1) / nb) * nb;
	const int i3 = forward ? nb : -nb;
	for (int i = i1; forward ? (i < k) : (i >= 0); i += i3) {
		const int ib = std::min(nb, k - i);
		tlarft('F', 'C', nq - i, ib, A + i + L * i, lda, tau + i, Tf.data(), nb);
		const int mi = left ? m - i : m;
		const int ni = left ? n : n - i;
		const int ic = left ? i : 0;
		const int jc = left ? 0 : i;
		tlarfb(side, trans, 'F', 'C', mi, ni, ib, A + i + L * i, lda, Tf.data(), nb,
			C + ic + LC * jc, ldc);
	}
	return 0;
}

// LQ 分解 (unblocked)
template <typename T>
inline int tgelq2(const int m, const int n, T* A, const int lda, T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgelq2: invalid argument");
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int k = std::min(m, n);
	for (int i = 0; i < k; i++) {
		tlarfg(n - i, A[i + L * i], A + i + L * std::min(i + 1, n - 1), lda, tau[i]);
		if (i < m - 1) {
			tlarf1f('R', m - i - 1, n - i, A + i + L * i, lda, tau[i], A + (i + 1) + L * i, lda);
		}
	}
	return 0;
}

// LQ 分解 (blocked)
template <typename T>
inline int tgelqf(const int m, const int n, T* A, const int lda, T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgelqf: invalid argument");
	}
	const int k = std::min(m, n);
	if (k == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "GELQF", m, n, -1, -1);
	const int nx = std::max(0, detail::ilaenv(3, "GELQF", m, n, -1, -1));
	int i = 0;
	if (nb >= 2 && nb < k && nx < k) {
		std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
		for (; i < k - nx; i += nb) {
			const int ib = std::min(k - i, nb);
			tgelq2(ib, n - i, A + i + L * i, lda, tau + i);
			if (i + ib < m) {
				tlarft('F', 'R', n - i, ib, A + i + L * i, lda, tau + i, Tf.data(), nb);
				tlarfb('R', 'N', 'F', 'R', m - i - ib, n - i, ib, A + i + L * i, lda,
					Tf.data(), nb, A + (i + ib) + L * i, lda);
			}
		}
	}
	if (i < k) {
		tgelq2(m - i, n - i, A + i + L * i, lda, tau + i);
	}
	return 0;
}

// LQ の Q の生成 (unblocked)
template <typename T>
inline int torgl2(const int m, const int n, const int k, T* A, const int lda, const T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < m || k < 0 || k > m || lda < std::max(1, m)) {
		detail::tlapack_error("torgl2: invalid argument");
	}
	if (m <= 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	if (k < m) {
		for (int j = 0; j < n; j++) {
			for (int l = k; l < m; l++) {
				A[l + L * j] = T(0);
			}
			if (j >= k && j < m) {
				A[j + L * j] = T(1);
			}
		}
	}
	for (int i = k - 1; i >= 0; i--) {
		if (i < n - 1) {
			if (i < m - 1) {
				tlarf1f('R', m - i - 1, n - i, A + i + L * i, lda, tau[i], A + (i + 1) + L * i, lda);
			}
			tscal(n - i - 1, -tau[i], A + i + L * (i + 1), lda);
		}
		A[i + L * i] = T(1) - tau[i];
		for (int l = 0; l < i; l++) {
			A[i + L * l] = T(0);
		}
	}
	return 0;
}

// LQ の Q の生成 (blocked)
template <typename T>
inline int torglq(const int m, const int n, const int k, T* A, const int lda, const T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < m || k < 0 || k > m || lda < std::max(1, m)) {
		detail::tlapack_error("torglq: invalid argument");
	}
	if (m <= 0) {
		return 0;
	}
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
				A[i + L * j] = T(0);
			}
		}
	}
	if (kk < m) {
		torgl2(m - kk, n - kk, k - kk, A + kk + L * kk, lda, tau + kk);
	}
	if (kk > 0) {
		std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
		for (int i = ki; i >= 0; i -= nb) {
			const int ib = std::min(nb, k - i);
			if (i + ib < m) {
				tlarft('F', 'R', n - i, ib, A + i + L * i, lda, tau + i, Tf.data(), nb);
				tlarfb('R', 'T', 'F', 'R', m - i - ib, n - i, ib, A + i + L * i, lda,
					Tf.data(), nb, A + (i + ib) + L * i, lda);
			}
			torgl2(ib, n - i, ib, A + i + L * i, lda, tau + i);
			for (int j = 0; j < i; j++) {
				for (int l = i; l < i + ib; l++) {
					A[l + L * j] = T(0);
				}
			}
		}
	}
	return 0;
}

// LQ の Q の適用 (unblocked)
template <typename T>
inline int torml2(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, k) || ldc < std::max(1, m)) {
		detail::tlapack_error("torml2: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
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
		tlarf1f(side, mi, ni, A + i + L * i, lda, tau[i], C + ic + LC * jc, ldc);
	}
	return 0;
}

// LQ の Q の適用 (blocked)
template <typename T>
inline int tormlq(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, k) || ldc < std::max(1, m)) {
		detail::tlapack_error("tormlq: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const std::size_t LC = static_cast<std::size_t>(ldc);
	const int nb = std::min(64, detail::ilaenv(1, "ORMLQ", m, n, k, -1));
	if (nb < 2 || nb >= k) {
		return torml2(side, trans, m, n, k, A, lda, tau, C, ldc);
	}
	std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
	const char transt = notran ? 'T' : 'N';
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : ((k - 1) / nb) * nb;
	const int i3 = forward ? nb : -nb;
	for (int i = i1; forward ? (i < k) : (i >= 0); i += i3) {
		const int ib = std::min(nb, k - i);
		tlarft('F', 'R', nq - i, ib, A + i + L * i, lda, tau + i, Tf.data(), nb);
		const int mi = left ? m - i : m;
		const int ni = left ? n : n - i;
		const int ic = left ? i : 0;
		const int jc = left ? 0 : i;
		tlarfb(side, transt, 'F', 'R', mi, ni, ib, A + i + L * i, lda, Tf.data(), nb,
			C + ic + LC * jc, ldc);
	}
	return 0;
}

// QL の Q の生成 (unblocked)
template <typename T>
inline int torg2l(const int m, const int n, const int k, T* A, const int lda, const T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::tlapack_error("torg2l: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	for (int j = 0; j < n - k; j++) {
		for (int l = 0; l < m; l++) {
			A[l + L * j] = T(0);
		}
		A[(m - n + j) + L * j] = T(1);
	}
	for (int i = 0; i < k; i++) {
		const int ii = n - k + i; // 0-based column
		tlarf1l('L', m - n + ii + 1, ii, A + L * ii, 1, tau[i], A, lda);
		tscal(m - n + ii, -tau[i], A + L * ii, 1);
		A[(m - n + ii) + L * ii] = T(1) - tau[i];
		for (int l = m - n + ii + 1; l < m; l++) {
			A[l + L * ii] = T(0);
		}
	}
	return 0;
}

// QL の Q の生成 (blocked)
template <typename T>
inline int torgql(const int m, const int n, const int k, T* A, const int lda, const T* tau) {
	namespace detail = tlapack_detail;
	if (m < 0 || n < 0 || n > m || k < 0 || k > n || lda < std::max(1, m)) {
		detail::tlapack_error("torgql: invalid argument");
	}
	if (n <= 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = detail::ilaenv(1, "ORGQL", m, n, k, -1);
	const int nx = std::max(0, detail::ilaenv(3, "ORGQL", m, n, k, -1));
	int kk = 0;
	if (nb >= 2 && nb < k && nx < k) {
		kk = std::min(k, ((k - nx + nb - 1) / nb) * nb);
		for (int j = 0; j < n - kk; j++) {
			for (int i = m - kk; i < m; i++) {
				A[i + L * j] = T(0);
			}
		}
	}
	torg2l(m - kk, n - kk, k - kk, A, lda, tau);
	if (kk > 0) {
		std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
		for (int i = k - kk; i < k; i += nb) {
			const int ib = std::min(nb, k - i);
			const int col = n - k + i; // block の先頭列 (0-based)
			if (col > 0) {
				tlarft('B', 'C', m - k + i + ib, ib, A + L * col, lda, tau + i, Tf.data(), nb);
				tlarfb('L', 'N', 'B', 'C', m - k + i + ib, col, ib, A + L * col, lda,
					Tf.data(), nb, A, lda);
			}
			torg2l(m - k + i + ib, ib, ib, A + L * col, lda, tau + i);
			for (int j = col; j < col + ib; j++) {
				for (int l = m - k + i + ib; l < m; l++) {
					A[l + L * j] = T(0);
				}
			}
		}
	}
	return 0;
}

// QL の Q の適用 (unblocked)
template <typename T>
inline int torm2l(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::tlapack_error("torm2l: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : k - 1;
	const int i2 = forward ? k : -1;
	const int i3 = forward ? 1 : -1;
	for (int i = i1; i != i2; i += i3) {
		const int mi = left ? m - k + i + 1 : m;
		const int ni = left ? n : n - k + i + 1;
		tlarf1l(side, mi, ni, A + L * i, 1, tau[i], C, ldc);
	}
	return 0;
}

// QL の Q の適用 (blocked)
template <typename T>
inline int tormql(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc) {
	namespace detail = tlapack_detail;
	const bool left = detail::option_is(side, 'L');
	const bool notran = detail::option_is(trans, 'N');
	const int nq = left ? m : n;
	if ((!left && !detail::option_is(side, 'R')) || (!notran && !detail::option_is(trans, 'T')) ||
		m < 0 || n < 0 || k < 0 || k > nq || lda < std::max(1, nq) || ldc < std::max(1, m)) {
		detail::tlapack_error("tormql: invalid argument");
	}
	if (m == 0 || n == 0 || k == 0) {
		return 0;
	}
	const std::size_t L = static_cast<std::size_t>(lda);
	const int nb = std::min(64, detail::ilaenv(1, "ORMQL", m, n, k, -1));
	if (nb < 2 || nb >= k) {
		return torm2l(side, trans, m, n, k, A, lda, tau, C, ldc);
	}
	std::vector<T> Tf(static_cast<std::size_t>(nb) * nb, T(0));
	const bool forward = (left && notran) || (!left && !notran);
	const int i1 = forward ? 0 : ((k - 1) / nb) * nb;
	const int i3 = forward ? nb : -nb;
	for (int i = i1; forward ? (i < k) : (i >= 0); i += i3) {
		const int ib = std::min(nb, k - i);
		tlarft('B', 'C', nq - k + i + ib, ib, A + L * i, lda, tau + i, Tf.data(), nb);
		const int mi = left ? m - k + i + ib : m;
		const int ni = left ? n : n - k + i + ib;
		tlarfb(side, trans, 'B', 'C', mi, ni, ib, A + L * i, lda, Tf.data(), nb, C, ldc);
	}
	return 0;
}

// 最小二乗問題 min ||B - A*X|| または min ||B - A^T*X|| の求解 (driver)
// B は max(m,n) x nrhs (解は先頭 n 行または m 行に上書き)
template <typename T>
inline int tgels(const char trans, const int m, const int n, const int nrhs,
	T* A, const int lda, T* B, const int ldb) {
	namespace detail = tlapack_detail;
	const bool tpsd = !detail::option_is(trans, 'N');
	if (tpsd && !detail::option_is(trans, 'T')) {
		detail::tlapack_error("tgels: invalid trans");
	}
	if (m < 0 || n < 0 || nrhs < 0 || lda < std::max(1, m) || ldb < std::max(1, std::max(m, n))) {
		detail::tlapack_error("tgels: invalid argument");
	}
	if (std::min(std::min(m, n), nrhs) == 0) {
		tlaset('F', std::max(m, n), nrhs, T(0), T(0), B, ldb);
		return 0;
	}
	const std::size_t LB = static_cast<std::size_t>(ldb);
	const int mn = std::min(m, n);
	std::vector<T> tau(static_cast<std::size_t>(mn), T(0));
	const T smlnum = tlamch<T>('S') / tlamch<T>('P');
	const T bignum = T(1) / smlnum;
	const T anrm = tlange('M', m, n, A, lda);
	int iascl = 0;
	if (anrm > T(0) && anrm < smlnum) {
		tlascl('G', 0, 0, anrm, smlnum, m, n, A, lda);
		iascl = 1;
	}
	else if (anrm > bignum) {
		tlascl('G', 0, 0, anrm, bignum, m, n, A, lda);
		iascl = 2;
	}
	else if (anrm == T(0)) {
		// A = 0: 最小二乗解は X = 0
		tlaset('F', std::max(m, n), nrhs, T(0), T(0), B, ldb);
		return 0;
	}
	const int brow = tpsd ? n : m;
	const T bnrm = tlange('M', brow, nrhs, B, ldb);
	int ibscl = 0;
	if (bnrm > T(0) && bnrm < smlnum) {
		tlascl('G', 0, 0, bnrm, smlnum, brow, nrhs, B, ldb);
		ibscl = 1;
	}
	else if (bnrm > bignum) {
		tlascl('G', 0, 0, bnrm, bignum, brow, nrhs, B, ldb);
		ibscl = 2;
	}
	int scllen = 0;
	if (m >= n) {
		tgeqrf(m, n, A, lda, tau.data());
		if (!tpsd) {
			// 最小二乗解: B := inv(R) * Q^T * B
			tormqr('L', 'T', m, nrhs, n, A, lda, tau.data(), B, ldb);
			const int info = ttrtrs('U', 'N', 'N', n, nrhs, A, lda, B, ldb);
			if (info > 0) {
				return info;
			}
			scllen = n;
		}
		else {
			// underdetermined A^T * X = B の最小 norm 解
			const int info = ttrtrs('U', 'T', 'N', n, nrhs, A, lda, B, ldb);
			if (info > 0) {
				return info;
			}
			for (int j = 0; j < nrhs; j++) {
				for (int i = n; i < m; i++) {
					B[i + LB * j] = T(0);
				}
			}
			tormqr('L', 'N', m, nrhs, n, A, lda, tau.data(), B, ldb);
			scllen = m;
		}
	}
	else {
		tgelqf(m, n, A, lda, tau.data());
		if (!tpsd) {
			// underdetermined A * X = B の最小 norm 解
			const int info = ttrtrs('L', 'N', 'N', m, nrhs, A, lda, B, ldb);
			if (info > 0) {
				return info;
			}
			for (int j = 0; j < nrhs; j++) {
				for (int i = m; i < n; i++) {
					B[i + LB * j] = T(0);
				}
			}
			tormlq('L', 'T', n, nrhs, m, A, lda, tau.data(), B, ldb);
			scllen = n;
		}
		else {
			// overdetermined A^T * X = B の最小二乗解
			tormlq('L', 'N', n, nrhs, m, A, lda, tau.data(), B, ldb);
			const int info = ttrtrs('L', 'T', 'N', m, nrhs, A, lda, B, ldb);
			if (info > 0) {
				return info;
			}
			scllen = m;
		}
	}
	if (iascl == 1) {
		tlascl('G', 0, 0, anrm, smlnum, scllen, nrhs, B, ldb);
	}
	else if (iascl == 2) {
		tlascl('G', 0, 0, anrm, bignum, scllen, nrhs, B, ldb);
	}
	if (ibscl == 1) {
		tlascl('G', 0, 0, smlnum, bnrm, scllen, nrhs, B, ldb);
	}
	else if (ibscl == 2) {
		tlascl('G', 0, 0, bignum, bnrm, scllen, nrhs, B, ldb);
	}
	return 0;
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_QR_HPP
