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

#ifndef VBLAS_RDBLAS_LEVEL3_HPP
#define VBLAS_RDBLAS_LEVEL3_HPP

// Level 3 BLAS (丸めモード指定付き)．行列は column-major．
// 引数は reference BLAS と同順 + 末尾に rounding_mode (1 上向き / -1 下向き / 0 最近点)．
//
// O(N^3) の核は rmatmul (AVX-512 / AVX2+FMA / NEON / no-SIMD を compile 時に
// 自動選択) に集約する．rdsyrk / rdsyr2k / rdtrmm / rdtrsm は再帰 2 分割で
// off-diagonal block を vcp::rdgemm に帰着させ，小さな対角 leaf のみ個別に処理する
// (余分な演算は O(N^2 * LEAF) で全体の数 % 以下)．
//
// 丸めモードの扱い:
// - rmatmul 内部は backend 毎に embedded rounding (AVX-512) または各 worker
//   thread での fesetround により指定丸めで計算し，終了後に各 thread の
//   丸めモードを復元する (vblas/README.md 参照)．
// - rmatmul 外の epilogue (alpha/beta 更新・三角更新・leaf solve) も全て
//   OpenMP parallel region 内で RoundingGuard により thread 毎に
//   保存 -> 変更 -> 復元する．

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <omp.h>

#include "rdblas_common.hpp"
#include "rmatmul.hpp"

#pragma STDC FENV_ACCESS ON

namespace vblas_rdblas_detail {

// O(N^3) 系の再帰分割 leaf size (off-diagonal は vcp::rdgemm 経由で rmatmul に乗る)
#ifndef VBLAS_RDBLAS_L3_LEAF
#define VBLAS_RDBLAS_L3_LEAF 128
#endif
// trsm の leaf は scalar solve になるため小さめにして scalar 部の比率を抑える
#ifndef VBLAS_RDBLAS_TRSM_LEAF
#define VBLAS_RDBLAS_TRSM_LEAF 64
#endif
const int RDBLAS_L3_LEAF = VBLAS_RDBLAS_L3_LEAF;
const int RDBLAS_TRSM_LEAF = VBLAS_RDBLAS_TRSM_LEAF;

inline int split_half(const int n) {
	const int half = ((n / 2 + 63) / 64) * 64;
	return std::min(half, n - 1);
}

// T(m x n 連続) := op(A)*op(B)．op が必要な場合や leading dimension が
// 詰まっていない場合は FP 演算なしの copy で連続 buffer に詰め，
// O(N^3) 本体は rmatmul (指定丸め) で計算する．
inline void matmul_into(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const double* A, const int lda,
	const double* B, const int ldb,
	double* T, const int rounding_mode
) {
	const bool ta = !option_is(transa, 'N');
	const bool tb = !option_is(transb, 'N');

	std::vector<double> buf_a;
	std::vector<double> buf_b;
	const double* a = A;
	const double* b = B;

	if (ta) {
		buf_a.resize(static_cast<std::size_t>(m) * k);
		pack_trans(m, k, A, lda, buf_a.data());
		a = buf_a.data();
	}
	else if (lda != m) {
		buf_a.resize(static_cast<std::size_t>(m) * k);
		pack_no_trans(m, k, A, lda, buf_a.data());
		a = buf_a.data();
	}

	if (tb) {
		buf_b.resize(static_cast<std::size_t>(k) * n);
		pack_trans(k, n, B, ldb, buf_b.data());
		b = buf_b.data();
	}
	else if (ldb != k) {
		buf_b.resize(static_cast<std::size_t>(k) * n);
		pack_no_trans(k, n, B, ldb, buf_b.data());
		b = buf_b.data();
	}

	// rmatmul の引数規約: rmatmul(M, N, K) で C(M x K) = A(M x N) * B(N x K)
	rmatmul(m, k, n, a, b, T, rounding_mode);
}

// 三角部分のみ C := beta*C (beta==0 は exact 0 埋め)
inline void scale_triangle(const bool upper, const int n, const double beta, double* C, const int ldc, const int fe_mode) {
	if (n <= 0 || beta == 1.0) {
		return;
	}
	const int threads = threads_for_flops(0.5 * n * n);
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe_mode);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			double* c = C + static_cast<std::size_t>(ldc) * j;
			const int i0 = upper ? 0 : j;
			const int i1 = upper ? j + 1 : n;
			if (beta == 0.0) {
				for (int i = i0; i < i1; i++) {
					c[i] = 0.0;
				}
			}
			else {
				for (int i = i0; i < i1; i++) {
					c[i] = beta * c[i];
				}
			}
		}
	}
}

// syrk/syr2k の band 用 epilogue:
// T1/T2 (ib x ncols 連続, ncols = i0+ib) は global 行 i0..i0+ib-1, 列 0..ncols-1 の
// 丸め付き積和の値を持つ．lower は C[i0+r, j] (j <= i0+r) をそのまま更新し，
// upper は「rmatmul が各要素を k 方向昇順の同一順序で積和する」ことから
// 対称位置の値が bit 単位で一致するため，同じ T の値で C[j, i0+r] を更新する．
inline void update_band_triangle(
	const bool upper, const int i0, const int ib, const int ncols,
	const double alpha, const double* T1, const double* T2,
	const double beta, double* C, const int ldc,
	const int fe_mode
) {
	const int threads = threads_for_flops(2.0 * ib * ncols);
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe_mode);
		if (!upper) {
#pragma omp for schedule(static)
			for (int j = 0; j < ncols; j++) {
				const int r0 = j <= i0 ? 0 : j - i0;
				double* col = C + i0 + static_cast<std::size_t>(ldc) * j;
				const double* t1 = T1 + static_cast<std::size_t>(ib) * j;
				if (T2 == NULL) {
					if (beta == 0.0) {
						for (int r = r0; r < ib; r++) {
							col[r] = alpha * t1[r];
						}
					}
					else {
						for (int r = r0; r < ib; r++) {
							col[r] = std::fma(beta, col[r], alpha * t1[r]);
						}
					}
				}
				else {
					const double* t2 = T2 + static_cast<std::size_t>(ib) * j;
					if (beta == 0.0) {
						for (int r = r0; r < ib; r++) {
							col[r] = std::fma(alpha, t2[r], alpha * t1[r]);
						}
					}
					else {
						for (int r = r0; r < ib; r++) {
							col[r] = std::fma(beta, col[r], std::fma(alpha, t2[r], alpha * t1[r]));
						}
					}
				}
			}
		}
		else {
#pragma omp for schedule(static)
			for (int r = 0; r < ib; r++) {
				double* col = C + static_cast<std::size_t>(ldc) * (i0 + r);
				const double* t1 = T1 + r;
				const int jend = i0 + r;
				if (T2 == NULL) {
					if (beta == 0.0) {
						for (int j = 0; j <= jend; j++) {
							col[j] = alpha * t1[static_cast<std::size_t>(ib) * j];
						}
					}
					else {
						for (int j = 0; j <= jend; j++) {
							col[j] = std::fma(beta, col[j], alpha * t1[static_cast<std::size_t>(ib) * j]);
						}
					}
				}
				else {
					const double* t2 = T2 + r;
					if (beta == 0.0) {
						for (int j = 0; j <= jend; j++) {
							col[j] = std::fma(alpha, t2[static_cast<std::size_t>(ib) * j], alpha * t1[static_cast<std::size_t>(ib) * j]);
						}
					}
					else {
						for (int j = 0; j <= jend; j++) {
							col[j] = std::fma(beta, col[j], std::fma(alpha, t2[static_cast<std::size_t>(ib) * j], alpha * t1[static_cast<std::size_t>(ib) * j]));
						}
					}
				}
			}
		}
	}
}

// 対称行列 (uplo の三角のみ格納) を full dense (n x n 連続) に展開する (copy のみ)
inline void expand_symmetric(const bool upper, const int n, const double* A, const int lda, double* dst) {
	const int threads = threads_for_flops(static_cast<double>(n) * n);
#pragma omp parallel for schedule(static) num_threads(threads)
	for (int j = 0; j < n; j++) {
		double* d = dst + static_cast<std::size_t>(n) * j;
		if (upper) {
			const double* a = A + static_cast<std::size_t>(lda) * j;
			for (int i = 0; i <= j; i++) {
				d[i] = a[i];
			}
			for (int i = j + 1; i < n; i++) {
				d[i] = A[j + static_cast<std::size_t>(lda) * i];
			}
		}
		else {
			for (int i = 0; i < j; i++) {
				d[i] = A[j + static_cast<std::size_t>(lda) * i];
			}
			const double* a = A + static_cast<std::size_t>(lda) * j;
			for (int i = j; i < n; i++) {
				d[i] = a[i];
			}
		}
	}
}

// 三角行列を 0 詰め (unit diag は 1.0) の full dense (n x n 連続) に展開する (copy のみ)
inline void expand_triangular(const bool upper, const bool nounit, const int n, const double* A, const int lda, double* dst) {
	for (int j = 0; j < n; j++) {
		const double* a = A + static_cast<std::size_t>(lda) * j;
		double* d = dst + static_cast<std::size_t>(n) * j;
		if (upper) {
			for (int i = 0; i < j; i++) {
				d[i] = a[i];
			}
			d[j] = nounit ? a[j] : 1.0;
			for (int i = j + 1; i < n; i++) {
				d[i] = 0.0;
			}
		}
		else {
			for (int i = 0; i < j; i++) {
				d[i] = 0.0;
			}
			d[j] = nounit ? a[j] : 1.0;
			for (int i = j + 1; i < n; i++) {
				d[i] = a[i];
			}
		}
	}
}

} // namespace vblas_rdblas_detail

namespace vcp {

// C := alpha*op(A)*op(B) + beta*C
// op(A): m x k, op(B): k x n, C: m x n
inline void rdgemm(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool ta = !det::option_is(transa, 'N');
	const bool tb = !det::option_is(transb, 'N');
	if (ta && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::rdblas_error("vcp::rdgemm: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::rdblas_error("vcp::rdgemm: invalid transb");
	}
	const int nrowa = ta ? k : m;
	const int nrowb = tb ? n : k;
	if (m < 0 || n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, m)) {
		det::rdblas_error("vcp::rdgemm: invalid argument");
	}
	if (m == 0 || n == 0 || (k == 0 && beta == 1.0) || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	if (alpha == 0.0 || k == 0) {
		det::scale_matrix(m, n, beta, C, ldc, fe);
		return;
	}

	if (alpha == 1.0 && beta == 0.0 && ldc == m) {
		det::matmul_into(transa, transb, m, n, k, A, lda, B, ldb, C, rounding_mode);
		return;
	}

	std::vector<double> T(static_cast<std::size_t>(m) * n);
	det::matmul_into(transa, transb, m, n, k, A, lda, B, ldb, T.data(), rounding_mode);
	det::apply_alpha_beta(m, n, alpha, T.data(), beta, C, ldc, fe);
}

// C := alpha*A*B + beta*C (side 'L') または C := alpha*B*A + beta*C (side 'R')
// A は対称行列 (uplo の三角のみ参照)，C: m x n
inline void rdsymm(
	const char side, const char uplo,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	if (!lside && !det::option_is(side, 'R')) {
		det::rdblas_error("rdsymm: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsymm: invalid uplo");
	}
	const int ka = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, ka) || ldb < std::max(1, m) || ldc < std::max(1, m)) {
		det::rdblas_error("rdsymm: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	if (alpha == 0.0) {
		det::scale_matrix(m, n, beta, C, ldc, fe);
		return;
	}

	// 対称三角格納を full dense に展開 (copy のみ) して gemm に帰着する
	std::vector<double> Afull(static_cast<std::size_t>(ka) * ka);
	det::expand_symmetric(upper, ka, A, lda, Afull.data());
	if (lside) {
		vcp::rdgemm('N', 'N', m, n, m, alpha, Afull.data(), m, B, ldb, beta, C, ldc, rounding_mode);
	}
	else {
		vcp::rdgemm('N', 'N', m, n, n, alpha, B, ldb, Afull.data(), n, beta, C, ldc, rounding_mode);
	}
}

} // namespace vcp

namespace vblas_rdblas_detail {

// syrk/syr2k 用 band 幅: 大きいほど rmatmul 呼び出しが太くなり，
// 一方で band 対角 tile の使わない三角部分の余分な flops (~ nb/(2n)) が増える．
// 既定 (0) は n に応じた自動調整 (n/4 を 64 の倍数に丸めて [128, 512] に clamp)
#ifndef VBLAS_RDBLAS_SYRK_BAND
#define VBLAS_RDBLAS_SYRK_BAND 0
#endif

inline int syrk_band_width(const int n) {
	if (VBLAS_RDBLAS_SYRK_BAND > 0) {
		return VBLAS_RDBLAS_SYRK_BAND;
	}
	const int rounded = ((n / 4 + 63) / 64) * 64;
	return std::max(128, std::min(512, rounded));
}

// op(A)^T (k x n) を一度だけ連続 buffer に作る (copy のみ)．
// 以後の band gemm の右オペランドはこの buffer の連続部分列なので再 pack 不要．
inline const double* prepare_op_trans(
	const bool ntrans, const int n, const int k,
	const double* A, const int lda, std::vector<double>& buf
) {
	if (ntrans) {
		buf.resize(static_cast<std::size_t>(k) * n);
		pack_trans(k, n, A, lda, buf.data());
		return buf.data();
	}
	if (lda == k) {
		return A;
	}
	buf.resize(static_cast<std::size_t>(k) * n);
	pack_no_trans(k, n, A, lda, buf.data());
	return buf.data();
}

// op(A) の行 i0..i0+ib-1 (ib x k) を連続 buffer に詰める (copy のみ)
inline void pack_op_rows(
	const bool ntrans, const int i0, const int ib, const int k,
	const double* A, const int lda, double* dst
) {
	if (ntrans) {
		pack_no_trans(ib, k, A + i0, lda, dst);
	}
	else {
		pack_trans(ib, k, A + static_cast<std::size_t>(lda) * i0, lda, dst);
	}
}

// 行 band ごとに T = op(A)[band] * op(A)^T[:, 0:i0+ib] を 1 回の rmatmul で計算し，
// 三角部分のみを alpha/beta 更新する．op(A)^T は事前に 1 回だけ pack するので
// band ごとの再 pack は左オペランド (ib x k) のみで済む．
inline void syrk_panel(
	const bool upper, const bool ntrans,
	const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	const int fe = fe_rounding(rounding_mode);
	std::vector<double> Rbuf;
	const double* R = prepare_op_trans(ntrans, n, k, A, lda, Rbuf);

	const int band = syrk_band_width(n);
	std::vector<double> L(static_cast<std::size_t>(std::min(band, n)) * k);
	std::vector<double> T;
	for (int i0 = 0; i0 < n; i0 += band) {
		const int ib = std::min(band, n - i0);
		const int ncols = i0 + ib;
		pack_op_rows(ntrans, i0, ib, k, A, lda, L.data());
		T.resize(static_cast<std::size_t>(ib) * ncols);
		rmatmul(ib, k, ncols, L.data(), R, T.data(), rounding_mode);
		update_band_triangle(upper, i0, ib, ncols, alpha, T.data(), NULL, beta, C, ldc, fe);
	}
}

inline void syr2k_panel(
	const bool upper, const bool ntrans,
	const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	const int fe = fe_rounding(rounding_mode);
	std::vector<double> RAbuf, RBbuf;
	const double* RA = prepare_op_trans(ntrans, n, k, A, lda, RAbuf);
	const double* RB = prepare_op_trans(ntrans, n, k, B, ldb, RBbuf);

	const int band = syrk_band_width(n);
	std::vector<double> LA(static_cast<std::size_t>(std::min(band, n)) * k);
	std::vector<double> LB(static_cast<std::size_t>(std::min(band, n)) * k);
	std::vector<double> T1, T2;
	for (int i0 = 0; i0 < n; i0 += band) {
		const int ib = std::min(band, n - i0);
		const int ncols = i0 + ib;
		pack_op_rows(ntrans, i0, ib, k, A, lda, LA.data());
		pack_op_rows(ntrans, i0, ib, k, B, ldb, LB.data());
		T1.resize(static_cast<std::size_t>(ib) * ncols);
		T2.resize(static_cast<std::size_t>(ib) * ncols);
		rmatmul(ib, k, ncols, LA.data(), RB, T1.data(), rounding_mode);
		rmatmul(ib, k, ncols, LB.data(), RA, T2.data(), rounding_mode);
		update_band_triangle(upper, i0, ib, ncols, alpha, T1.data(), T2.data(), beta, C, ldc, fe);
	}
}

// trmm leaf: 三角 block を 0 詰め dense に展開し (copy のみ)，
// T = op(A)*B (または B*op(A)) を rmatmul で計算して B := alpha*T で書き戻す．
// 0 との fma は結果を変えないため，余分な演算は結果に影響しない．
inline void trmm_leaf(
	const bool lside, const bool upper, const bool ntrans, const bool nounit,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	const int nb = lside ? m : n;
	if (nb <= 0 || m <= 0 || n <= 0) {
		return;
	}
	const int fe = fe_rounding(rounding_mode);
	std::vector<double> Adense(static_cast<std::size_t>(nb) * nb);
	expand_triangular(upper, nounit, nb, A, lda, Adense.data());
	std::vector<double> T(static_cast<std::size_t>(m) * n);
	if (lside) {
		matmul_into(ntrans ? 'N' : 'T', 'N', m, n, m, Adense.data(), m, B, ldb, T.data(), rounding_mode);
	}
	else {
		matmul_into('N', ntrans ? 'N' : 'T', m, n, n, B, ldb, Adense.data(), n, T.data(), rounding_mode);
	}
	apply_alpha_beta(m, n, alpha, T.data(), 0.0, B, ldb, fe);
}

inline void trmm_recursive(
	const bool lside, const bool upper, const bool ntrans, const bool nounit,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	if (m <= 0 || n <= 0) {
		return;
	}
	const int na = lside ? m : n;
	if (na <= RDBLAS_L3_LEAF) {
		trmm_leaf(lside, upper, ntrans, nounit, m, n, alpha, A, lda, B, ldb, rounding_mode);
		return;
	}
	const int n1 = split_half(na);
	const int n2 = na - n1;
	const double* A11 = A;
	const double* A21 = A + n1;
	const double* A12 = A + static_cast<std::size_t>(lda) * n1;
	const double* A22 = A + n1 + static_cast<std::size_t>(lda) * n1;

	if (lside) {
		double* B1 = B;
		double* B2 = B + n1;
		if (upper == ntrans) {
			// upper&N: B1 := alpha*A11*B1 + alpha*A12*B2 / lower&T: B1 := alpha*A11^T*B1 + alpha*A21^T*B2
			trmm_recursive(lside, upper, ntrans, nounit, n1, n, alpha, A11, lda, B1, ldb, rounding_mode);
			vcp::rdgemm(ntrans ? 'N' : 'T', 'N', n1, n, n2, alpha, ntrans ? A12 : A21, lda, B2, ldb, 1.0, B1, ldb, rounding_mode);
			trmm_recursive(lside, upper, ntrans, nounit, n2, n, alpha, A22, lda, B2, ldb, rounding_mode);
		}
		else {
			// lower&N: B2 := alpha*A21*B1 + alpha*A22*B2 / upper&T: B2 := alpha*A12^T*B1 + alpha*A22^T*B2
			trmm_recursive(lside, upper, ntrans, nounit, n2, n, alpha, A22, lda, B2, ldb, rounding_mode);
			vcp::rdgemm(ntrans ? 'N' : 'T', 'N', n2, n, n1, alpha, ntrans ? A21 : A12, lda, B1, ldb, 1.0, B2, ldb, rounding_mode);
			trmm_recursive(lside, upper, ntrans, nounit, n1, n, alpha, A11, lda, B1, ldb, rounding_mode);
		}
	}
	else {
		double* B1 = B;
		double* B2 = B + static_cast<std::size_t>(ldb) * n1;
		if (upper == ntrans) {
			// upper&N: B2 := alpha*B1*A12 + alpha*B2*A22 / lower&T: B2 := alpha*B1*A21^T + alpha*B2*A22^T
			trmm_recursive(lside, upper, ntrans, nounit, m, n2, alpha, A22, lda, B2, ldb, rounding_mode);
			vcp::rdgemm('N', ntrans ? 'N' : 'T', m, n2, n1, alpha, B1, ldb, ntrans ? A12 : A21, lda, 1.0, B2, ldb, rounding_mode);
			trmm_recursive(lside, upper, ntrans, nounit, m, n1, alpha, A11, lda, B1, ldb, rounding_mode);
		}
		else {
			// lower&N: B1 := alpha*B1*A11 + alpha*B2*A21 / upper&T: B1 := alpha*B1*A11^T + alpha*B2*A12^T
			trmm_recursive(lside, upper, ntrans, nounit, m, n1, alpha, A11, lda, B1, ldb, rounding_mode);
			vcp::rdgemm('N', ntrans ? 'N' : 'T', m, n1, n2, alpha, B2, ldb, ntrans ? A21 : A12, lda, 1.0, B1, ldb, rounding_mode);
			trmm_recursive(lside, upper, ntrans, nounit, m, n2, alpha, A22, lda, B2, ldb, rounding_mode);
		}
	}
}

// trsm の left 側 leaf: B の各列を独立に三角 solve する (列方向に OpenMP 並列)
inline void trsm_left_leaf(
	const bool upper, const bool ntrans, const bool nounit,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	const int fe = fe_rounding(rounding_mode);
	const int threads = threads_for_flops(static_cast<double>(m) * m * n);
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			double* b = B + static_cast<std::size_t>(ldb) * j;
			if (alpha != 1.0) {
				for (int i = 0; i < m; i++) {
					b[i] = alpha * b[i];
				}
			}
			if (ntrans) {
				if (upper) {
					for (int jj = m - 1; jj >= 0; jj--) {
						const double* a = A + static_cast<std::size_t>(lda) * jj;
						if (b[jj] != 0.0) {
							if (nounit) {
								b[jj] = b[jj] / a[jj];
							}
							const double temp = b[jj];
							for (int i = 0; i < jj; i++) {
								b[i] = std::fma(-temp, a[i], b[i]);
							}
						}
					}
				}
				else {
					for (int jj = 0; jj < m; jj++) {
						const double* a = A + static_cast<std::size_t>(lda) * jj;
						if (b[jj] != 0.0) {
							if (nounit) {
								b[jj] = b[jj] / a[jj];
							}
							const double temp = b[jj];
							for (int i = jj + 1; i < m; i++) {
								b[i] = std::fma(-temp, a[i], b[i]);
							}
						}
					}
				}
			}
			else {
				if (upper) {
					// A^T は lower: 前進代入 (dot 形式)
					for (int jj = 0; jj < m; jj++) {
						const double* a = A + static_cast<std::size_t>(lda) * jj;
						double temp = b[jj];
						for (int i = 0; i < jj; i++) {
							temp = std::fma(-a[i], b[i], temp);
						}
						if (nounit) {
							temp = temp / a[jj];
						}
						b[jj] = temp;
					}
				}
				else {
					// A^T は upper: 後退代入 (dot 形式)
					for (int jj = m - 1; jj >= 0; jj--) {
						const double* a = A + static_cast<std::size_t>(lda) * jj;
						double temp = b[jj];
						for (int i = jj + 1; i < m; i++) {
							temp = std::fma(-a[i], b[i], temp);
						}
						if (nounit) {
							temp = temp / a[jj];
						}
						b[jj] = temp;
					}
				}
			}
		}
	}
}

// trsm の right 側 leaf: B の行 block ごとに独立な列 sweep solve
// X_j = (alpha*B_j - sum_{先行 i} X_i * coef(i,j)) / diag(j)
inline void trsm_right_leaf(
	const bool upper, const bool ntrans, const bool nounit,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	const int fe = fe_rounding(rounding_mode);
	const int threads = threads_for_flops(static_cast<double>(m) * n * n);
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe);
		const int nth = omp_get_num_threads();
		const int tid = omp_get_thread_num();
		const int chunk = (m + nth - 1) / nth;
		const int r0 = std::min(m, tid * chunk);
		const int r1 = std::min(m, r0 + chunk);
		// 列の処理順: 先行列のみに依存するように選ぶ
		const bool ascending = (upper == ntrans);
		for (int jj = 0; jj < n; jj++) {
			const int j = ascending ? jj : n - 1 - jj;
			double* bj = B + static_cast<std::size_t>(ldb) * j;
			if (alpha != 1.0) {
				for (int r = r0; r < r1; r++) {
					bj[r] = alpha * bj[r];
				}
			}
			const int i0 = ascending ? 0 : j + 1;
			const int i1 = ascending ? j : n;
			for (int i = i0; i < i1; i++) {
				// coef(i,j): op(A)(i,j)，ntrans なら A(i,j)，trans なら A(j,i)
				const double aij = ntrans ? A[i + static_cast<std::size_t>(lda) * j] : A[j + static_cast<std::size_t>(lda) * i];
				if (aij == 0.0) {
					continue;
				}
				const double* bi = B + static_cast<std::size_t>(ldb) * i;
				for (int r = r0; r < r1; r++) {
					bj[r] = std::fma(-aij, bi[r], bj[r]);
				}
			}
			if (nounit) {
				const double ajj = A[j + static_cast<std::size_t>(lda) * j];
				for (int r = r0; r < r1; r++) {
					bj[r] = bj[r] / ajj;
				}
			}
		}
	}
}

inline void trsm_recursive(
	const bool lside, const bool upper, const bool ntrans, const bool nounit,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	if (m <= 0 || n <= 0) {
		return;
	}
	const int na = lside ? m : n;
	if (na <= RDBLAS_TRSM_LEAF) {
		if (lside) {
			trsm_left_leaf(upper, ntrans, nounit, m, n, alpha, A, lda, B, ldb, rounding_mode);
		}
		else {
			trsm_right_leaf(upper, ntrans, nounit, m, n, alpha, A, lda, B, ldb, rounding_mode);
		}
		return;
	}
	const int n1 = split_half(na);
	const int n2 = na - n1;
	const double* A11 = A;
	const double* A21 = A + n1;
	const double* A12 = A + static_cast<std::size_t>(lda) * n1;
	const double* A22 = A + n1 + static_cast<std::size_t>(lda) * n1;

	if (lside) {
		double* B1 = B;
		double* B2 = B + n1;
		if (upper == ntrans) {
			// upper&N / lower&T: 後退型 (B2 を先に solve)
			trsm_recursive(lside, upper, ntrans, nounit, n2, n, alpha, A22, lda, B2, ldb, rounding_mode);
			vcp::rdgemm(ntrans ? 'N' : 'T', 'N', n1, n, n2, -1.0, ntrans ? A12 : A21, lda, B2, ldb, alpha, B1, ldb, rounding_mode);
			trsm_recursive(lside, upper, ntrans, nounit, n1, n, 1.0, A11, lda, B1, ldb, rounding_mode);
		}
		else {
			// lower&N / upper&T: 前進型 (B1 を先に solve)
			trsm_recursive(lside, upper, ntrans, nounit, n1, n, alpha, A11, lda, B1, ldb, rounding_mode);
			vcp::rdgemm(ntrans ? 'N' : 'T', 'N', n2, n, n1, -1.0, ntrans ? A21 : A12, lda, B1, ldb, alpha, B2, ldb, rounding_mode);
			trsm_recursive(lside, upper, ntrans, nounit, n2, n, 1.0, A22, lda, B2, ldb, rounding_mode);
		}
	}
	else {
		double* B1 = B;
		double* B2 = B + static_cast<std::size_t>(ldb) * n1;
		if (upper == ntrans) {
			// upper&N / lower&T: B1 を先に solve
			trsm_recursive(lside, upper, ntrans, nounit, m, n1, alpha, A11, lda, B1, ldb, rounding_mode);
			vcp::rdgemm('N', ntrans ? 'N' : 'T', m, n2, n1, -1.0, B1, ldb, ntrans ? A12 : A21, lda, alpha, B2, ldb, rounding_mode);
			trsm_recursive(lside, upper, ntrans, nounit, m, n2, 1.0, A22, lda, B2, ldb, rounding_mode);
		}
		else {
			// lower&N / upper&T: B2 を先に solve
			trsm_recursive(lside, upper, ntrans, nounit, m, n2, alpha, A22, lda, B2, ldb, rounding_mode);
			vcp::rdgemm('N', ntrans ? 'N' : 'T', m, n1, n2, -1.0, B2, ldb, ntrans ? A21 : A12, lda, alpha, B1, ldb, rounding_mode);
			trsm_recursive(lside, upper, ntrans, nounit, m, n1, 1.0, A11, lda, B1, ldb, rounding_mode);
		}
	}
}

} // namespace vblas_rdblas_detail

namespace vcp {

// C := alpha*op(A)*op(A)^T + beta*C (uplo の三角のみ更新)
// trans 'N': op(A) = A (n x k), trans 'T'/'C': op(A) = A^T (A は k x n)
inline void rdsyrk(
	const char uplo, const char trans,
	const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsyrk: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::rdblas_error("rdsyrk: invalid trans");
	}
	const int nrowa = ntrans ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldc < std::max(1, n)) {
		det::rdblas_error("rdsyrk: invalid argument");
	}
	if (n == 0 || (alpha == 0.0 && beta == 1.0) || (k == 0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	if (alpha == 0.0 || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc, fe);
		return;
	}
	det::syrk_panel(upper, ntrans, n, k, alpha, A, lda, beta, C, ldc, rounding_mode);
}

// C := alpha*op(A)*op(B)^T + alpha*op(B)*op(A)^T + beta*C (uplo の三角のみ更新)
inline void rdsyr2k(
	const char uplo, const char trans,
	const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdsyr2k: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::rdblas_error("rdsyr2k: invalid trans");
	}
	const int nrowa = ntrans ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowa) || ldc < std::max(1, n)) {
		det::rdblas_error("rdsyr2k: invalid argument");
	}
	if (n == 0 || (alpha == 0.0 && beta == 1.0) || (k == 0 && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	if (alpha == 0.0 || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc, fe);
		return;
	}
	det::syr2k_panel(upper, ntrans, n, k, alpha, A, lda, B, ldb, beta, C, ldc, rounding_mode);
}

// B := alpha*op(A)*B (side 'L') または B := alpha*B*op(A) (side 'R')
// A は三角行列，B: m x n
inline void rdtrmm(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(transa, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!lside && !det::option_is(side, 'R')) {
		det::rdblas_error("rdtrmm: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdtrmm: invalid uplo");
	}
	if (!ntrans && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::rdblas_error("rdtrmm: invalid transa");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::rdblas_error("rdtrmm: invalid diag");
	}
	const int nrowa = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, m)) {
		det::rdblas_error("rdtrmm: invalid argument");
	}
	if (m == 0 || n == 0) {
		return;
	}
	if (alpha == 0.0) {
		det::scale_matrix(m, n, 0.0, B, ldb, det::fe_rounding(rounding_mode));
		return;
	}
	det::trmm_recursive(lside, upper, ntrans, nounit, m, n, alpha, A, lda, B, ldb, rounding_mode);
}

// op(A)*X = alpha*B (side 'L') または X*op(A) = alpha*B (side 'R') を解き，
// B を解 X で上書きする．A は三角行列，B: m x n
inline void rdtrsm(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const double alpha, const double* A, const int lda,
	double* B, const int ldb,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(transa, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!lside && !det::option_is(side, 'R')) {
		det::rdblas_error("rdtrsm: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdtrsm: invalid uplo");
	}
	if (!ntrans && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::rdblas_error("rdtrsm: invalid transa");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::rdblas_error("rdtrsm: invalid diag");
	}
	const int nrowa = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, m)) {
		det::rdblas_error("rdtrsm: invalid argument");
	}
	if (m == 0 || n == 0) {
		return;
	}
	if (alpha == 0.0) {
		det::scale_matrix(m, n, 0.0, B, ldb, det::fe_rounding(rounding_mode));
		return;
	}
	det::trsm_recursive(lside, upper, ntrans, nounit, m, n, alpha, A, lda, B, ldb, rounding_mode);
}

} // namespace vcp

namespace vblas_rdblas_detail {

// 対角 block (nb x nb) の三角部分のみ C := alpha*T + beta*C (gemmtr 用 epilogue)
// 演算は apply_alpha_beta と同じ形 (beta==0 のとき C は読まない)
inline void gemmtr_diag_update(
	const bool upper, const int nb,
	const double alpha, const double* T,
	const double beta, double* C, const int ldc,
	const int fe_mode
) {
	const int threads = threads_for_flops(static_cast<double>(nb) * nb);
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe_mode);
#pragma omp for schedule(static)
		for (int j = 0; j < nb; j++) {
			const double* t = T + static_cast<std::size_t>(nb) * j;
			double* c = C + static_cast<std::size_t>(ldc) * j;
			const int i0 = upper ? 0 : j;
			const int i1 = upper ? j : nb - 1;
			if (beta == 0.0) {
				if (alpha == 1.0) {
					for (int i = i0; i <= i1; i++) {
						c[i] = t[i];
					}
				}
				else {
					for (int i = i0; i <= i1; i++) {
						c[i] = alpha * t[i];
					}
				}
			}
			else if (alpha == 1.0) {
				for (int i = i0; i <= i1; i++) {
					c[i] = std::fma(beta, c[i], t[i]);
				}
			}
			else {
				for (int i = i0; i <= i1; i++) {
					c[i] = std::fma(beta, c[i], alpha * t[i]);
				}
			}
		}
	}
}

} // namespace vblas_rdblas_detail

namespace vcp {

// C の uplo 三角部分のみ := alpha*op(A)*op(B) + beta*C (C: n x n, op(A): n x k,
// op(B): k x n)．reference BLAS の GEMMTR (LAPACK 3.12.1 で追加) に対応する．
// 列 block ごとに長方形部分を vcp::rdgemm，対角 block は temp に積を作って
// 三角部分のみ指定丸めの epilogue で更新する．
inline void rdgemmtr(
	const char uplo, const char transa, const char transb,
	const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc,
	const int rounding_mode
) {
	namespace det = vblas_rdblas_detail;
	const bool upper = det::option_is(uplo, 'U');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::rdblas_error("rdgemmtr: invalid uplo");
	}
	const bool ta = !det::option_is(transa, 'N');
	const bool tb = !det::option_is(transb, 'N');
	if (ta && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::rdblas_error("rdgemmtr: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::rdblas_error("rdgemmtr: invalid transb");
	}
	const int nrowa = ta ? k : n;
	const int nrowb = tb ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, n)) {
		det::rdblas_error("rdgemmtr: invalid argument");
	}
	if (n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0)) {
		return;
	}
	const int fe = det::fe_rounding(rounding_mode);
	if (alpha == 0.0 || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc, fe);
		return;
	}
	const int nb = det::RDBLAS_L3_LEAF;
	std::vector<double> T(static_cast<std::size_t>(std::min(nb, n)) * std::min(nb, n));
	for (int j0 = 0; j0 < n; j0 += nb) {
		const int jb = std::min(nb, n - j0);
		const double* opb = tb ? B + j0 : B + static_cast<std::size_t>(ldb) * j0;
		if (upper) {
			// 長方形部分 C(0:j0, j0:j0+jb)
			if (j0 > 0) {
				vcp::rdgemm(transa, transb, j0, jb, k, alpha, A, lda, opb, ldb,
					beta, C + static_cast<std::size_t>(ldc) * j0, ldc, rounding_mode);
			}
		}
		else {
			// 長方形部分 C(j0+jb:n, j0:j0+jb)
			const int m2 = n - j0 - jb;
			if (m2 > 0) {
				const double* opa2 = ta ? A + static_cast<std::size_t>(lda) * (j0 + jb) : A + (j0 + jb);
				vcp::rdgemm(transa, transb, m2, jb, k, alpha, opa2, lda, opb, ldb,
					beta, C + (j0 + jb) + static_cast<std::size_t>(ldc) * j0, ldc, rounding_mode);
			}
		}
		// 対角 block C(j0:j0+jb, j0:j0+jb) の三角部分
		const double* opa = ta ? A + static_cast<std::size_t>(lda) * j0 : A + j0;
		det::matmul_into(transa, transb, jb, jb, k, opa, lda, opb, ldb, T.data(), rounding_mode);
		det::gemmtr_diag_update(upper, jb, alpha, T.data(), beta,
			C + j0 + static_cast<std::size_t>(ldc) * j0, ldc, fe);
	}
}

} // namespace vcp

#endif // VBLAS_RDBLAS_LEVEL3_HPP
