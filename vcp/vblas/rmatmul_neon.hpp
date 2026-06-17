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

#pragma once

#ifndef VBLAS_RMATMUL_NEON_HPP
#define VBLAS_RMATMUL_NEON_HPP

#include <algorithm>
#include <cfenv>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#if !(defined(__aarch64__) || defined(__arm64__)) || !(defined(__ARM_NEON) || defined(__ARM_NEON__))
#error "rmatmul_neon.hpp requires AArch64 NEON compiler support."
#endif

#include <arm_neon.h>
#include <omp.h>

#include "rmatmul_common.hpp"

#pragma STDC FENV_ACCESS ON

namespace vblas_rmatmul_neon_detail {

constexpr int VEC = 2;
constexpr int MR = 8;
constexpr int NR = 4;
constexpr int MR_CHUNKS = MR / VEC;

namespace common = vblas_rmatmul_common_detail;

using common::CpuInfo;
using common::Blocking;
using common::ParallelMode;
using common::PARALLEL_2D;
using common::PARALLEL_M_ONLY;
using common::ceil_div;
using common::round_up_to;
using common::round_down_to;
using common::cpu_info;
using common::limited_threads;
using common::set_rounding_mode;
using common::get_rounding_mode;
using common::reset_rounding_mode;

// blocking 決定は rmatmul_common.hpp の共通実装に委譲する
inline Blocking select_blocking(const int m, const int n, const int k) {
	return common::select_blocking<MR, NR>(m, n, k);
}

inline float64x2_t load_partial_f64(const double* src, const int rows) {
	double tmp[VEC] = {0.0, 0.0};
	if (rows >= 1) {
		tmp[0] = src[0];
	}
	if (rows >= 2) {
		tmp[1] = src[1];
	}
	return vld1q_f64(tmp);
}

inline void store_partial_f64(double* dst, const int rows, const float64x2_t value) {
	double tmp[VEC];
	vst1q_f64(tmp, value);
	if (rows >= 1) {
		dst[0] = tmp[0];
	}
	if (rows >= 2) {
		dst[1] = tmp[1];
	}
}

template<int ROUNDING>
inline void kernel_mrx4(
	const int pb,
	const double* A,
	const int lda,
	const double* B,
	const int ldb,
	double* CU,
	const int ldc,
	const bool init_zero
) {
	const float64x2_t zero = vdupq_n_f64(0.0);
	float64x2_t c00 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 0);
	float64x2_t c01 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 2);
	float64x2_t c02 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 4);
	float64x2_t c03 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 6);
	float64x2_t c10 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 0);
	float64x2_t c11 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 2);
	float64x2_t c12 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 4);
	float64x2_t c13 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 6);
	float64x2_t c20 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 0);
	float64x2_t c21 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 2);
	float64x2_t c22 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 4);
	float64x2_t c23 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 6);
	float64x2_t c30 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 0);
	float64x2_t c31 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 2);
	float64x2_t c32 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 4);
	float64x2_t c33 = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 6);

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		const float64x2_t a0 = vld1q_f64(ap + 0);
		const float64x2_t a1 = vld1q_f64(ap + 2);
		const float64x2_t a2 = vld1q_f64(ap + 4);
		const float64x2_t a3 = vld1q_f64(ap + 6);

		const float64x2_t b0 = vdupq_n_f64(B[p + static_cast<std::size_t>(ldb) * 0]);
		c00 = vfmaq_f64(c00, a0, b0);
		c01 = vfmaq_f64(c01, a1, b0);
		c02 = vfmaq_f64(c02, a2, b0);
		c03 = vfmaq_f64(c03, a3, b0);

		const float64x2_t b1 = vdupq_n_f64(B[p + static_cast<std::size_t>(ldb) * 1]);
		c10 = vfmaq_f64(c10, a0, b1);
		c11 = vfmaq_f64(c11, a1, b1);
		c12 = vfmaq_f64(c12, a2, b1);
		c13 = vfmaq_f64(c13, a3, b1);

		const float64x2_t b2 = vdupq_n_f64(B[p + static_cast<std::size_t>(ldb) * 2]);
		c20 = vfmaq_f64(c20, a0, b2);
		c21 = vfmaq_f64(c21, a1, b2);
		c22 = vfmaq_f64(c22, a2, b2);
		c23 = vfmaq_f64(c23, a3, b2);

		const float64x2_t b3 = vdupq_n_f64(B[p + static_cast<std::size_t>(ldb) * 3]);
		c30 = vfmaq_f64(c30, a0, b3);
		c31 = vfmaq_f64(c31, a1, b3);
		c32 = vfmaq_f64(c32, a2, b3);
		c33 = vfmaq_f64(c33, a3, b3);
	}

	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 0, c00);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 2, c01);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 4, c02);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 0 + 6, c03);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 0, c10);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 2, c11);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 4, c12);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 1 + 6, c13);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 0, c20);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 2, c21);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 4, c22);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 2 + 6, c23);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 0, c30);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 2, c31);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 4, c32);
	vst1q_f64(CU + static_cast<std::size_t>(ldc) * 3 + 6, c33);
}

template<int ROUNDING>
inline void kernel_mrx4_masked(
	const int pb,
	const int mr,
	const int nr,
	const double* A,
	const int lda,
	const double* B,
	const int ldb,
	double* CU,
	const int ldc,
	const bool init_zero
) {
	const float64x2_t zero = vdupq_n_f64(0.0);
	float64x2_t cu[NR][MR_CHUNKS];

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			const int rows = std::min(VEC, std::max(0, mr - VEC * r));
			cu[j][r] = init_zero ? zero : load_partial_f64(CU + static_cast<std::size_t>(ldc) * j + VEC * r, rows);
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		float64x2_t a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			const int rows = std::min(VEC, std::max(0, mr - VEC * r));
			a[r] = load_partial_f64(ap + VEC * r, rows);
		}
		for (int j = 0; j < nr; j++) {
			const float64x2_t b = vdupq_n_f64(B[p + static_cast<std::size_t>(ldb) * j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = vfmaq_f64(cu[j][r], a[r], b);
			}
		}
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			const int rows = std::min(VEC, std::max(0, mr - VEC * r));
			store_partial_f64(CU + static_cast<std::size_t>(ldc) * j + VEC * r, rows, cu[j][r]);
		}
	}
}

// packed micro-panel 用 kernel: A は stride MR, B は stride NR の完全連続アクセス
template<int ROUNDING>
inline void kernel_packed_mrx4(
	const int pb,
	const double* A,
	const double* B,
	double* CU,
	const int ldc,
	const bool init_zero
) {
	const float64x2_t zero = vdupq_n_f64(0.0);
	float64x2_t cu[NR][MR_CHUNKS];

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = init_zero ? zero : vld1q_f64(CU + static_cast<std::size_t>(ldc) * j + VEC * r);
		}
	}

	for (int p = 0; p < pb; p++) {
		float64x2_t a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = vld1q_f64(A + VEC * r);
		}
		for (int j = 0; j < NR; j++) {
			const float64x2_t b = vdupq_n_f64(B[j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = vfmaq_f64(cu[j][r], a[r], b);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			vst1q_f64(CU + static_cast<std::size_t>(ldc) * j + VEC * r, cu[j][r]);
		}
	}
}

// 端 tile 用: A の端数行と B の端数列は pack 時に 0 詰め済みなので
// load は常に full 幅で行い，CU の load/store だけ部分アクセスにする．
// fma の 0 オペランドは結果を変えない (正確に恒等) ため計算結果は不変．
template<int ROUNDING>
inline void kernel_packed_mrx4_masked(
	const int pb,
	const int mr,
	const int nr,
	const double* A,
	const double* B,
	double* CU,
	const int ldc,
	const bool init_zero
) {
	const float64x2_t zero = vdupq_n_f64(0.0);
	float64x2_t cu[NR][MR_CHUNKS];

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			if (init_zero || j >= nr) {
				cu[j][r] = zero;
			}
			else {
				const int rows = std::min(VEC, std::max(0, mr - VEC * r));
				cu[j][r] = load_partial_f64(CU + static_cast<std::size_t>(ldc) * j + VEC * r, rows);
			}
		}
	}

	for (int p = 0; p < pb; p++) {
		float64x2_t a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = vld1q_f64(A + VEC * r);
		}
		for (int j = 0; j < NR; j++) {
			const float64x2_t b = vdupq_n_f64(B[j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = vfmaq_f64(cu[j][r], a[r], b);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			const int rows = std::min(VEC, std::max(0, mr - VEC * r));
			store_partial_f64(CU + static_cast<std::size_t>(ldc) * j + VEC * r, rows, cu[j][r]);
		}
	}
}

// micro-panel packing は rmatmul_common.hpp の共通実装に委譲する
inline void pack_a_micro(const int m, const int mb, const int pb, const double* A, double* packed) {
	common::pack_a_micro<MR>(m, mb, pb, A, packed);
}

inline void pack_b_micro_parallel(const int n, const int pb, const int nb, const double* B, double* packed) {
	common::pack_b_micro_parallel<NR>(n, pb, nb, B, packed);
}

template<int ROUNDING>
inline void compute_packed_ic_block(
	const int m,
	const int pb,
	const int ng,
	const int ic,
	const int pc,
	const int jcg,
	const Blocking& blocking,
	const double* A,
	const double* packed_b,
	double* packed_a,
	double* CU,
	const bool init_zero
) {
	const int mb = std::min(blocking.mc, m - ic);
	pack_a_micro(m, mb, pb, A + ic + static_cast<std::size_t>(m) * pc, packed_a);

	for (int jc_offset = 0; jc_offset < ng; jc_offset += blocking.column_tile) {
		const int nb = std::min(blocking.column_tile, ng - jc_offset);
		for (int jr = 0; jr < nb; jr += NR) {
			const int nr = std::min(NR, nb - jr);
			const double* bp = packed_b + static_cast<std::size_t>((jc_offset + jr) / NR) * pb * NR;
			for (int ir = 0; ir < mb; ir += MR) {
				const int mr = std::min(MR, mb - ir);
				const double* ap = packed_a + static_cast<std::size_t>(ir / MR) * MR * pb;
				double* cup = CU + (ic + ir) + static_cast<std::size_t>(m) * (jcg + jc_offset + jr);
				if (mr == MR && nr == NR) {
					kernel_packed_mrx4<ROUNDING>(pb, ap, bp, cup, m, init_zero);
				}
				else {
					kernel_packed_mrx4_masked<ROUNDING>(pb, mr, nr, ap, bp, cup, m, init_zero);
				}
			}
		}
	}
}

template<int ROUNDING>
inline void compute_packed_ic_column_tile(
	const int m,
	const int pb,
	const int nb,
	const int ic,
	const int pc,
	const int j_base,
	const Blocking& blocking,
	const double* A,
	const double* packed_b_tile,
	double* packed_a,
	double* CU,
	const bool init_zero
) {
	const int mb = std::min(blocking.mc, m - ic);
	pack_a_micro(m, mb, pb, A + ic + static_cast<std::size_t>(m) * pc, packed_a);

	for (int jr = 0; jr < nb; jr += NR) {
		const int nr = std::min(NR, nb - jr);
		const double* bp = packed_b_tile + static_cast<std::size_t>(jr / NR) * pb * NR;
		for (int ir = 0; ir < mb; ir += MR) {
			const int mr = std::min(MR, mb - ir);
			const double* ap = packed_a + static_cast<std::size_t>(ir / MR) * MR * pb;
			double* cup = CU + (ic + ir) + static_cast<std::size_t>(m) * (j_base + jr);
			if (mr == MR && nr == NR) {
				kernel_packed_mrx4<ROUNDING>(pb, ap, bp, cup, m, init_zero);
			}
			else {
				kernel_packed_mrx4_masked<ROUNDING>(pb, mr, nr, ap, bp, cup, m, init_zero);
			}
		}
	}
}

inline void zero_output(const int m, const int k, double* CU, const int threads) {
	const long long size = static_cast<long long>(m) * k;
#pragma omp parallel for schedule(static) num_threads(threads)
	for (long long i = 0; i < size; i++) {
		CU[i] = 0.0;
	}
}

template<int ROUNDING>
inline void rmatmul_direct_neon(
	const int m,
	const int n,
	const int k,
	const double* A,
	const double* B,
	double* CU,
	const int threads
) {
#pragma omp parallel num_threads(threads)
	{
		const int previous_rounding = get_rounding_mode();
		set_rounding_mode(ROUNDING);
#pragma omp for collapse(2) schedule(static)
		for (int jr = 0; jr < k; jr += NR) {
			for (int ir = 0; ir < m; ir += MR) {
				const int nr = std::min(NR, k - jr);
				const int mr = std::min(MR, m - ir);
				const double* ap = A + ir;
				const double* bp = B + static_cast<std::size_t>(n) * jr;
				double* cup = CU + ir + static_cast<std::size_t>(m) * jr;
				if (mr == MR && nr == NR) {
					kernel_mrx4<ROUNDING>(n, ap, m, bp, n, cup, m, true);
				}
				else {
					kernel_mrx4_masked<ROUNDING>(n, mr, nr, ap, m, bp, n, cup, m, true);
				}
			}
		}
		reset_rounding_mode(previous_rounding);
	}
}

template<int ROUNDING>
inline void rmatmul_blocked_m_only_neon(
	const int m,
	const int n,
	const int k,
	const double* A,
	const double* B,
	double* CU,
	const Blocking& blocking
) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * round_up_to(std::min(blocking.nc_group, k), NR));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);
		const int previous_rounding = get_rounding_mode();
		set_rounding_mode(ROUNDING);

		for (int pc = 0; pc < n; pc += blocking.kc) {
			const int pb = std::min(blocking.kc, n - pc);
			const bool init_zero = pc == 0;

			for (int jcg = 0; jcg < k; jcg += blocking.nc_group) {
				const int ng = std::min(blocking.nc_group, k - jcg);

#pragma omp single
				{
					packed_b.resize(static_cast<std::size_t>(pb) * round_up_to(ng, NR));
				}
				pack_b_micro_parallel(n, pb, ng, B + pc + static_cast<std::size_t>(n) * jcg, packed_b.data());

#pragma omp for schedule(static)
				for (int ic = 0; ic < m; ic += blocking.mc) {
					compute_packed_ic_block<ROUNDING>(m, pb, ng, ic, pc, jcg, blocking, A, packed_b.data(), packed_a.data(), CU, init_zero);
				}
			}
		}
		reset_rounding_mode(previous_rounding);
	}
}

template<int ROUNDING>
inline void rmatmul_blocked_2d_neon(
	const int m,
	const int n,
	const int k,
	const double* A,
	const double* B,
	double* CU,
	const Blocking& blocking
) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * round_up_to(std::min(blocking.nc_group, k), NR));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);
		const int previous_rounding = get_rounding_mode();
		set_rounding_mode(ROUNDING);

		for (int pc = 0; pc < n; pc += blocking.kc) {
			const int pb = std::min(blocking.kc, n - pc);
			const bool init_zero = pc == 0;

			for (int jcg = 0; jcg < k; jcg += blocking.nc_group) {
				const int ng = std::min(blocking.nc_group, k - jcg);

#pragma omp single
				{
					packed_b.resize(static_cast<std::size_t>(pb) * round_up_to(ng, NR));
				}
				pack_b_micro_parallel(n, pb, ng, B + pc + static_cast<std::size_t>(n) * jcg, packed_b.data());

#pragma omp for collapse(2) schedule(static)
				for (int jc_offset = 0; jc_offset < ng; jc_offset += blocking.column_tile) {
					for (int ic = 0; ic < m; ic += blocking.mc) {
						const int nb = std::min(blocking.column_tile, ng - jc_offset);
						compute_packed_ic_column_tile<ROUNDING>(
							m, pb, nb, ic, pc, jcg + jc_offset,
							blocking, A, packed_b.data() + static_cast<std::size_t>(jc_offset / NR) * pb * NR, packed_a.data(), CU, init_zero);
					}
				}
			}
		}
		reset_rounding_mode(previous_rounding);
	}
}

} // namespace vblas_rmatmul_neon_detail

template<int ROUNDING>
inline void rmatmul_impl_neon(int m, int n, int k, const double* A, const double* B, double* CU) {
	using namespace vblas_rmatmul_neon_detail;

	if (m <= 0 || k <= 0) {
		return;
	}

	const Blocking blocking = select_blocking(m, n, k);
	if (n <= 0) {
		zero_output(m, k, CU, blocking.threads);
		return;
	}

	if (blocking.direct) {
		rmatmul_direct_neon<ROUNDING>(m, n, k, A, B, CU, blocking.threads);
	}
	else if (blocking.mode == PARALLEL_M_ONLY) {
		rmatmul_blocked_m_only_neon<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else {
		rmatmul_blocked_2d_neon<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
}

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU is overwritten with A*B using rounding_mode: 1 upward, 0 nearest, -1 downward.
inline void rmatmul_neon(int m, int n, int k, const double* A, const double* B, double* CU, const int rounding_mode) {
	const int previous_rounding = vblas_rmatmul_neon_detail::get_rounding_mode();
	if (rounding_mode == 1) {
		rmatmul_impl_neon<FE_UPWARD>(m, n, k, A, B, CU);
	}
	else if (rounding_mode == -1) {
		rmatmul_impl_neon<FE_DOWNWARD>(m, n, k, A, B, CU);
	}
	else {
		rmatmul_impl_neon<FE_TONEAREST>(m, n, k, A, B, CU);
	}
	vblas_rmatmul_neon_detail::reset_rounding_mode(previous_rounding);
}

#endif // VBLAS_RMATMUL_NEON_HPP
