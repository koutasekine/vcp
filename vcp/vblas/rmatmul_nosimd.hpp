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

#ifndef VBLAS_RMATMUL_NOSIMD_HPP
#define VBLAS_RMATMUL_NOSIMD_HPP

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <omp.h>

#include "rmatmul_common.hpp"

#pragma STDC FENV_ACCESS ON

namespace vblas_rmatmul_nosimd_detail {

constexpr int MR = 8;
constexpr int NR = 4;

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

template<int ROUNDING>
inline double rounded_madd(const double acc, const double a, const double b) {
	return std::fma(a, b, acc);
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
	double cu[NR][MR];

	for (int j = 0; j < NR; j++) {
		for (int i = 0; i < MR; i++) {
			cu[j][i] = init_zero ? 0.0 : CU[static_cast<std::size_t>(ldc) * j + i];
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		for (int j = 0; j < NR; j++) {
			const double b = B[p + static_cast<std::size_t>(ldb) * j];
			for (int i = 0; i < MR; i++) {
				cu[j][i] = rounded_madd<ROUNDING>(cu[j][i], ap[i], b);
			}
		}
	}

	for (int j = 0; j < NR; j++) {
		for (int i = 0; i < MR; i++) {
			CU[static_cast<std::size_t>(ldc) * j + i] = cu[j][i];
		}
	}
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
	double cu[NR][MR];

	for (int j = 0; j < nr; j++) {
		for (int i = 0; i < mr; i++) {
			cu[j][i] = init_zero ? 0.0 : CU[static_cast<std::size_t>(ldc) * j + i];
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		for (int j = 0; j < nr; j++) {
			const double b = B[p + static_cast<std::size_t>(ldb) * j];
			for (int i = 0; i < mr; i++) {
				cu[j][i] = rounded_madd<ROUNDING>(cu[j][i], ap[i], b);
			}
		}
	}

	for (int j = 0; j < nr; j++) {
		for (int i = 0; i < mr; i++) {
			CU[static_cast<std::size_t>(ldc) * j + i] = cu[j][i];
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
	double cu[NR][MR];

	for (int j = 0; j < NR; j++) {
		for (int i = 0; i < MR; i++) {
			cu[j][i] = init_zero ? 0.0 : CU[static_cast<std::size_t>(ldc) * j + i];
		}
	}

	for (int p = 0; p < pb; p++) {
		for (int j = 0; j < NR; j++) {
			const double b = B[j];
			for (int i = 0; i < MR; i++) {
				cu[j][i] = rounded_madd<ROUNDING>(cu[j][i], A[i], b);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < NR; j++) {
		for (int i = 0; i < MR; i++) {
			CU[static_cast<std::size_t>(ldc) * j + i] = cu[j][i];
		}
	}
}

// 端 tile 用: packed A/B を stride MR/NR で読み，有効な mr 行 nr 列のみ計算する
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
	double cu[NR][MR];

	for (int j = 0; j < nr; j++) {
		for (int i = 0; i < mr; i++) {
			cu[j][i] = init_zero ? 0.0 : CU[static_cast<std::size_t>(ldc) * j + i];
		}
	}

	for (int p = 0; p < pb; p++) {
		for (int j = 0; j < nr; j++) {
			const double b = B[j];
			for (int i = 0; i < mr; i++) {
				cu[j][i] = rounded_madd<ROUNDING>(cu[j][i], A[i], b);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < nr; j++) {
		for (int i = 0; i < mr; i++) {
			CU[static_cast<std::size_t>(ldc) * j + i] = cu[j][i];
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
inline void rmatmul_direct_nosimd(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
#pragma omp parallel num_threads(blocking.threads)
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
inline void rmatmul_blocked_m_only_nosimd(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
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
inline void rmatmul_blocked_2d_nosimd(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
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

} // namespace vblas_rmatmul_nosimd_detail

template<int ROUNDING>
inline void rmatmul_impl_nosimd(int m, int n, int k, const double* A, const double* B, double* CU) {
	using namespace vblas_rmatmul_nosimd_detail;

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
		rmatmul_direct_nosimd<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else if (blocking.mode == PARALLEL_M_ONLY) {
		rmatmul_blocked_m_only_nosimd<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else {
		rmatmul_blocked_2d_nosimd<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
}

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU is overwritten with A*B using rounding_mode: 1 upward, 0 nearest, -1 downward.
inline void rmatmul_nosimd(int m, int n, int k, const double* A, const double* B, double* CU, const int rounding_mode) {
	const int previous_rounding = vblas_rmatmul_nosimd_detail::get_rounding_mode();
	if (rounding_mode == 1) {
		rmatmul_impl_nosimd<FE_UPWARD>(m, n, k, A, B, CU);
	}
	else if (rounding_mode == -1) {
		rmatmul_impl_nosimd<FE_DOWNWARD>(m, n, k, A, B, CU);
	}
	else {
		rmatmul_impl_nosimd<FE_TONEAREST>(m, n, k, A, B, CU);
	}
	vblas_rmatmul_nosimd_detail::reset_rounding_mode(previous_rounding);
}

#endif // VBLAS_RMATMUL_NOSIMD_HPP
