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

#ifndef VBLAS_RMATMUL_AVX512_HPP
#define VBLAS_RMATMUL_AVX512_HPP

#include <algorithm>
#include <cfenv>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <immintrin.h>
#include <omp.h>

#include "rmatmul_common.hpp"

namespace vblas_rmatmul_avx512_detail {

constexpr int VEC = 8;
constexpr int MR = 16;
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

inline __mmask8 row_mask(const int rows) {
	if (rows <= 0) {
		return static_cast<__mmask8>(0);
	}
	if (rows >= VEC) {
		return static_cast<__mmask8>(0xff);
	}
	return static_cast<__mmask8>((1u << rows) - 1u);
}

// micro-panel packing は rmatmul_common.hpp の共通実装に委譲する
inline void pack_a_micro(const int m, const int mb, const int pb, const double* A, double* packed) {
	common::pack_a_micro<MR>(m, mb, pb, A, packed);
}

inline void pack_b_micro_parallel(const int n, const int pb, const int nb, const double* B, double* packed) {
	common::pack_b_micro_parallel<NR>(n, pb, nb, B, packed);
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
	const __m512d zero = _mm512_setzero_pd();
	__m512d cu0l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 0);
	__m512d cu0h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 0 + VEC);
	__m512d cu1l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 1);
	__m512d cu1h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 1 + VEC);
	__m512d cu2l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 2);
	__m512d cu2h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 2 + VEC);
	__m512d cu3l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 3);
	__m512d cu3h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 3 + VEC);

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		const __m512d a0 = _mm512_loadu_pd(ap);
		const __m512d a1 = _mm512_loadu_pd(ap + VEC);
		const __m512d b0 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 0]);
		const __m512d b1 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 1]);
		const __m512d b2 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 2]);
		const __m512d b3 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 3]);

		cu0l = _mm512_fmadd_round_pd(a0, b0, cu0l, ROUNDING);
		cu0h = _mm512_fmadd_round_pd(a1, b0, cu0h, ROUNDING);
		cu1l = _mm512_fmadd_round_pd(a0, b1, cu1l, ROUNDING);
		cu1h = _mm512_fmadd_round_pd(a1, b1, cu1h, ROUNDING);
		cu2l = _mm512_fmadd_round_pd(a0, b2, cu2l, ROUNDING);
		cu2h = _mm512_fmadd_round_pd(a1, b2, cu2h, ROUNDING);
		cu3l = _mm512_fmadd_round_pd(a0, b3, cu3l, ROUNDING);
		cu3h = _mm512_fmadd_round_pd(a1, b3, cu3h, ROUNDING);
	}

	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 0, cu0l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 0 + VEC, cu0h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 1, cu1l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 1 + VEC, cu1h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 2, cu2l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 2 + VEC, cu2h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 3, cu3l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 3 + VEC, cu3h);
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
	const __m512d zero = _mm512_setzero_pd();
	__m512d cu[NR][MR_CHUNKS];
	__mmask8 mask[MR_CHUNKS];

	for (int r = 0; r < MR_CHUNKS; r++) {
		mask[r] = row_mask(std::min(VEC, std::max(0, mr - VEC * r)));
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = init_zero ? zero : _mm512_maskz_loadu_pd(mask[r], CU + static_cast<std::size_t>(ldc) * j + VEC * r);
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		__m512d a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = _mm512_maskz_loadu_pd(mask[r], ap + VEC * r);
		}
		for (int j = 0; j < nr; j++) {
			const __m512d b = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = _mm512_fmadd_round_pd(a[r], b, cu[j][r], ROUNDING);
			}
		}
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			_mm512_mask_storeu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, mask[r], cu[j][r]);
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
	const __m512d zero = _mm512_setzero_pd();
	__m512d cu0l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 0);
	__m512d cu0h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 0 + VEC);
	__m512d cu1l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 1);
	__m512d cu1h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 1 + VEC);
	__m512d cu2l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 2);
	__m512d cu2h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 2 + VEC);
	__m512d cu3l = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 3);
	__m512d cu3h = init_zero ? zero : _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 3 + VEC);

	for (int p = 0; p < pb; p++) {
		const __m512d a0 = _mm512_loadu_pd(A);
		const __m512d a1 = _mm512_loadu_pd(A + VEC);
		const __m512d b0 = _mm512_set1_pd(B[0]);
		const __m512d b1 = _mm512_set1_pd(B[1]);
		const __m512d b2 = _mm512_set1_pd(B[2]);
		const __m512d b3 = _mm512_set1_pd(B[3]);

		cu0l = _mm512_fmadd_round_pd(a0, b0, cu0l, ROUNDING);
		cu0h = _mm512_fmadd_round_pd(a1, b0, cu0h, ROUNDING);
		cu1l = _mm512_fmadd_round_pd(a0, b1, cu1l, ROUNDING);
		cu1h = _mm512_fmadd_round_pd(a1, b1, cu1h, ROUNDING);
		cu2l = _mm512_fmadd_round_pd(a0, b2, cu2l, ROUNDING);
		cu2h = _mm512_fmadd_round_pd(a1, b2, cu2h, ROUNDING);
		cu3l = _mm512_fmadd_round_pd(a0, b3, cu3l, ROUNDING);
		cu3h = _mm512_fmadd_round_pd(a1, b3, cu3h, ROUNDING);

		A += MR;
		B += NR;
	}

	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 0, cu0l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 0 + VEC, cu0h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 1, cu1l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 1 + VEC, cu1h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 2, cu2l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 2 + VEC, cu2h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 3, cu3l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 3 + VEC, cu3h);
}

// 端 tile 用: A の端数行と B の端数列は pack 時に 0 詰め済みなので
// load は常に full 幅で行い，CU の load/store だけ mask する．
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
	const __m512d zero = _mm512_setzero_pd();
	__m512d cu[NR][MR_CHUNKS];
	__mmask8 mask[MR_CHUNKS];

	for (int r = 0; r < MR_CHUNKS; r++) {
		mask[r] = row_mask(std::min(VEC, std::max(0, mr - VEC * r)));
	}

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = (init_zero || j >= nr) ? zero : _mm512_maskz_loadu_pd(mask[r], CU + static_cast<std::size_t>(ldc) * j + VEC * r);
		}
	}

	for (int p = 0; p < pb; p++) {
		__m512d a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = _mm512_loadu_pd(A + VEC * r);
		}
		for (int j = 0; j < NR; j++) {
			const __m512d b = _mm512_set1_pd(B[j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = _mm512_fmadd_round_pd(a[r], b, cu[j][r], ROUNDING);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			_mm512_mask_storeu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, mask[r], cu[j][r]);
		}
	}
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
inline void rmatmul_direct_avx512(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
#pragma omp parallel num_threads(blocking.threads)
	{
		// embedded rounding のため丸めモードは変更しないが，
		// 各 thread で実行前の丸めモードを保存し実行後に復元することを保証する
		const int previous_rounding = std::fegetround();
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
		std::fesetround(previous_rounding);
	}
}

template<int ROUNDING>
inline void rmatmul_blocked_m_only_avx512(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * round_up_to(std::min(blocking.nc_group, k), NR));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);
		const int previous_rounding = std::fegetround();

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
		std::fesetround(previous_rounding);
	}
}

template<int ROUNDING>
inline void rmatmul_blocked_2d_avx512(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * round_up_to(std::min(blocking.nc_group, k), NR));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);
		const int previous_rounding = std::fegetround();

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
		std::fesetround(previous_rounding);
	}
}

} // namespace vblas_rmatmul_avx512_detail

template<int ROUNDING>
inline void rmatmul_impl_avx512(int m, int n, int k, const double* A, const double* B, double* CU) {
	using namespace vblas_rmatmul_avx512_detail;

	if (m <= 0 || k <= 0) {
		return;
	}

	const CpuInfo& info = cpu_info();
	if (n <= 0) {
		zero_output(m, k, CU, limited_threads(std::min(info.logical_threads, info.omp_max_threads)));
		return;
	}

	const Blocking blocking = select_blocking(m, n, k);
	if (blocking.direct) {
		rmatmul_direct_avx512<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else if (blocking.mode == PARALLEL_M_ONLY) {
		rmatmul_blocked_m_only_avx512<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else {
		rmatmul_blocked_2d_avx512<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
}

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU is overwritten with A*B using rounding_mode: 1 upward, 0 nearest, -1 downward.
inline void rmatmul_avx512(int m, int n, int k, const double* A, const double* B, double* CU, const int rounding_mode) {
	// 呼び出し元 thread の丸めモードを保存し，終了後に復元する
	const int previous_rounding = std::fegetround();
	if (rounding_mode == 1) {
		rmatmul_impl_avx512<_MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC>(m, n, k, A, B, CU);
	}
	else if (rounding_mode == -1) {
		rmatmul_impl_avx512<_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC>(m, n, k, A, B, CU);
	}
	else {
		rmatmul_impl_avx512<_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC>(m, n, k, A, B, CU);
	}
	std::fesetround(previous_rounding);
}

#endif // VBLAS_RMATMUL_AVX512_HPP
