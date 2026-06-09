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
#include <vector>

#if !(defined(__aarch64__) || defined(__arm64__)) || !(defined(__ARM_NEON) || defined(__ARM_NEON__))
#error "rmatmul_neon.hpp requires AArch64 NEON compiler support."
#endif

#include <arm_neon.h>
#include <omp.h>

#pragma STDC FENV_ACCESS ON

namespace vblas_rmatmul_neon_detail {

constexpr int VEC = 2;
constexpr int MR = 8;
constexpr int NR = 4;
constexpr int MR_CHUNKS = MR / VEC;

enum ParallelMode {
	PARALLEL_2D,
	PARALLEL_M_ONLY
};

struct Blocking {
	int mc;
	int kc;
	int nc_group;
	int column_tile;
	int threads;
	ParallelMode mode;
	bool direct;
};

inline int ceil_div(const int n, const int block) {
	return (n + block - 1) / block;
}

inline int round_up_to(const int n, const int block) {
	return ((n + block - 1) / block) * block;
}

inline int round_down_to(const int n, const int block) {
	return (n / block) * block;
}

inline int limited_threads(const int requested) {
	return std::max(1, std::min(requested, std::max(1, omp_get_max_threads())));
}

inline int default_thread_count(const int m, const int n, const int k) {
	const int max_threads = limited_threads(omp_get_num_procs());
	const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
	if (work < 1.0e6) {
		return 1;
	}
	if (work >= 4.0e6 * static_cast<double>(max_threads)) {
		return max_threads;
	}
	return std::max(1, std::min(max_threads, static_cast<int>(work / 4.0e6) + 1));
}

inline Blocking sanitize_blocking(const Blocking& b) {
	Blocking r = b;
	r.mc = std::max(MR, round_up_to(r.mc, MR));
	r.kc = std::max(NR, round_up_to(r.kc, NR));
	r.nc_group = std::max(NR, round_up_to(r.nc_group, NR));
	r.column_tile = std::max(NR, round_up_to(r.column_tile, NR));
	r.column_tile = std::min(r.column_tile, r.nc_group);
	r.threads = limited_threads(r.threads);
	return r;
}

inline Blocking cap_threads_for_work(Blocking b, const int m, const int k) {
	if (b.direct) {
		const int tasks = ceil_div(m, MR) * ceil_div(k, NR);
		b.threads = std::max(1, std::min(b.threads, std::max(1, tasks)));
		return b;
	}

	const int m_blocks = ceil_div(m, b.mc);
	const int column_blocks = ceil_div(std::min(k, b.nc_group), b.column_tile);
	const int tasks = b.mode == PARALLEL_2D ? m_blocks * std::max(1, column_blocks) : m_blocks;
	b.threads = std::max(1, std::min(b.threads, std::max(1, tasks)));
	return b;
}

inline Blocking adjust_mc_for_m_only(Blocking b, const int m) {
	if (b.mode != PARALLEL_M_ONLY || b.threads <= 1 || m <= 0) {
		return b;
	}

	const int initial_blocks = ceil_div(m, b.mc);
	if (initial_blocks < b.threads) {
		const int raw = m / b.threads;
		const int mc_up = std::max(MR, round_up_to(raw == 0 ? MR : raw, MR));
		if (ceil_div(m, mc_up) >= b.threads) {
			b.mc = mc_up;
		}
		else {
			const int rounded_down = (raw / MR) * MR;
			b.mc = std::max(MR, rounded_down > 0 ? rounded_down : MR);
		}
	}

	const int blocks = ceil_div(m, b.mc);
	const int max_load = ceil_div(blocks, b.threads);
	const double current_eff = static_cast<double>(blocks) / (static_cast<double>(max_load) * b.threads);
	if (current_eff >= 0.90) {
		return b;
	}

	const int base_mc = b.mc;
	const int min_mc = std::max(MR, round_up_to(base_mc * 3 / 4, MR));
	int best_mc = base_mc;
	double best_eff = current_eff;

	for (int mc = base_mc; mc >= min_mc; mc -= MR) {
		const int ic = ceil_div(m, mc);
		if (ic < b.threads) {
			break;
		}
		const int load = ceil_div(ic, b.threads);
		const double eff = static_cast<double>(ic) / (static_cast<double>(load) * b.threads);
		if (eff > best_eff) {
			best_eff = eff;
			best_mc = mc;
			if (eff >= 0.99) {
				break;
			}
		}
	}

	b.mc = best_mc;
	return b;
}

inline Blocking finish_blocking(Blocking b, const int m, const int k) {
	b = sanitize_blocking(b);
	b = adjust_mc_for_m_only(b, m);
	b = cap_threads_for_work(b, m, k);
	return sanitize_blocking(b);
}

inline bool use_direct_path(const int m, const int n, const int k) {
	const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
	if (work <= 5.0e6) {
		return true;
	}
	if (k <= 64 && m <= 8192 && n <= 8192) {
		return true;
	}
	return static_cast<double>(m) * static_cast<double>(n) <= 65536.0 && n <= 1024;
}

inline Blocking select_blocking(const int m, const int n, const int k) {
	const int threads = default_thread_count(m, n, k);

	if (use_direct_path(m, n, k)) {
		const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
		const int direct_threads = work >= 1.0e7 ? threads : std::min(threads, 8);
		return finish_blocking(Blocking{MR, std::max(NR, n), std::max(NR, k), NR, direct_threads, PARALLEL_2D, true}, m, k);
	}

	int kc = 256;
	int column_tile = 128;
	int mc = 128;
	int nc_group = 1024;
	ParallelMode mode = ceil_div(m, mc) >= std::max(1, threads / 2) ? PARALLEL_M_ONLY : PARALLEL_2D;

	const int max_nk = std::max(n, k);
	const int max_mn = std::max(m, n);
	const int max_dim = std::max(m, max_nk);

	if (m >= 4096 && m >= 4 * max_nk) {
		mc = 192;
		column_tile = 192;
		mode = PARALLEL_M_ONLY;
	}
	else if (k >= 4096 && k >= 4 * max_mn) {
		mc = 128;
		column_tile = 96;
		nc_group = 768;
		mode = PARALLEL_2D;
	}
	else if (n >= 4096 && n >= 4 * std::max(m, k)) {
		kc = 192;
		mc = 128;
		column_tile = 96;
		nc_group = 768;
		mode = PARALLEL_2D;
	}
	else if (max_dim <= 1500) {
		kc = 192;
		mc = 96;
		column_tile = 96;
		nc_group = 768;
		mode = ceil_div(m, mc) >= std::max(1, threads / 2) ? PARALLEL_M_ONLY : PARALLEL_2D;
	}

	nc_group = std::max(column_tile, nc_group);
	return finish_blocking(Blocking{mc, kc, nc_group, column_tile, threads, mode, false}, m, k);
}

inline void set_rounding_mode(const int rounding) {
	std::fesetround(rounding);
}

inline void reset_rounding_mode() {
	std::fesetround(FE_TONEAREST);
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

inline void pack_a_panel(
	const int m,
	const int mb,
	const int pb,
	const int mb_padded,
	const double* A,
	double* packed
) {
	for (int p = 0; p < pb; p++) {
		double* dst = packed + static_cast<std::size_t>(p) * mb_padded;
		const double* src = A + static_cast<std::size_t>(m) * p;
		for (int i = 0; i < mb; i++) {
			dst[i] = src[i];
		}
		for (int i = mb; i < mb_padded; i++) {
			dst[i] = 0.0;
		}
	}
}

inline void pack_b_panel_parallel(
	const int n,
	const int pb,
	const int nb,
	const double* B,
	double* packed
) {
#pragma omp for schedule(static)
	for (int j = 0; j < nb; j++) {
		double* dst = packed + static_cast<std::size_t>(pb) * j;
		const double* src = B + static_cast<std::size_t>(n) * j;
		for (int p = 0; p < pb; p++) {
			dst[p] = src[p];
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
	const int mb_padded = round_up_to(mb, MR);
	pack_a_panel(m, mb, pb, mb_padded, A + ic + static_cast<std::size_t>(m) * pc, packed_a);

	for (int jc_offset = 0; jc_offset < ng; jc_offset += blocking.column_tile) {
		const int nb = std::min(blocking.column_tile, ng - jc_offset);
		for (int jr = 0; jr < nb; jr += NR) {
			const int nr = std::min(NR, nb - jr);
			const double* bp = packed_b + static_cast<std::size_t>(pb) * (jc_offset + jr);
			for (int ir = 0; ir < mb; ir += MR) {
				const int mr = std::min(MR, mb - ir);
				const double* ap = packed_a + ir;
				double* cup = CU + (ic + ir) + static_cast<std::size_t>(m) * (jcg + jc_offset + jr);
				if (mr == MR && nr == NR) {
					kernel_mrx4<ROUNDING>(pb, ap, mb_padded, bp, pb, cup, m, init_zero);
				}
				else {
					kernel_mrx4_masked<ROUNDING>(pb, mr, nr, ap, mb_padded, bp, pb, cup, m, init_zero);
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
	const int mb_padded = round_up_to(mb, MR);
	pack_a_panel(m, mb, pb, mb_padded, A + ic + static_cast<std::size_t>(m) * pc, packed_a);

	for (int jr = 0; jr < nb; jr += NR) {
		const int nr = std::min(NR, nb - jr);
		const double* bp = packed_b_tile + static_cast<std::size_t>(pb) * jr;
		for (int ir = 0; ir < mb; ir += MR) {
			const int mr = std::min(MR, mb - ir);
			const double* ap = packed_a + ir;
			double* cup = CU + (ic + ir) + static_cast<std::size_t>(m) * (j_base + jr);
			if (mr == MR && nr == NR) {
				kernel_mrx4<ROUNDING>(pb, ap, mb_padded, bp, pb, cup, m, init_zero);
			}
			else {
				kernel_mrx4_masked<ROUNDING>(pb, mr, nr, ap, mb_padded, bp, pb, cup, m, init_zero);
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
		reset_rounding_mode();
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
	packed_b.reserve(static_cast<std::size_t>(max_pb) * std::min(blocking.nc_group, k));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);
		set_rounding_mode(ROUNDING);

		for (int pc = 0; pc < n; pc += blocking.kc) {
			const int pb = std::min(blocking.kc, n - pc);
			const bool init_zero = pc == 0;

			for (int jcg = 0; jcg < k; jcg += blocking.nc_group) {
				const int ng = std::min(blocking.nc_group, k - jcg);

#pragma omp single
				{
					packed_b.resize(static_cast<std::size_t>(pb) * ng);
				}
				pack_b_panel_parallel(n, pb, ng, B + pc + static_cast<std::size_t>(n) * jcg, packed_b.data());

#pragma omp for schedule(static)
				for (int ic = 0; ic < m; ic += blocking.mc) {
					compute_packed_ic_block<ROUNDING>(m, pb, ng, ic, pc, jcg, blocking, A, packed_b.data(), packed_a.data(), CU, init_zero);
				}
			}
		}
		reset_rounding_mode();
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
	packed_b.reserve(static_cast<std::size_t>(max_pb) * std::min(blocking.nc_group, k));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);
		set_rounding_mode(ROUNDING);

		for (int pc = 0; pc < n; pc += blocking.kc) {
			const int pb = std::min(blocking.kc, n - pc);
			const bool init_zero = pc == 0;

			for (int jcg = 0; jcg < k; jcg += blocking.nc_group) {
				const int ng = std::min(blocking.nc_group, k - jcg);

#pragma omp single
				{
					packed_b.resize(static_cast<std::size_t>(pb) * ng);
				}
				pack_b_panel_parallel(n, pb, ng, B + pc + static_cast<std::size_t>(n) * jcg, packed_b.data());

#pragma omp for collapse(2) schedule(static)
				for (int jc_offset = 0; jc_offset < ng; jc_offset += blocking.column_tile) {
					for (int ic = 0; ic < m; ic += blocking.mc) {
						const int nb = std::min(blocking.column_tile, ng - jc_offset);
						compute_packed_ic_column_tile<ROUNDING>(
							m, pb, nb, ic, pc, jcg + jc_offset,
							blocking, A, packed_b.data() + static_cast<std::size_t>(pb) * jc_offset, packed_a.data(), CU, init_zero);
					}
				}
			}
		}
		reset_rounding_mode();
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
	if (rounding_mode == 1) {
		rmatmul_impl_neon<FE_UPWARD>(m, n, k, A, B, CU);
	}
	else if (rounding_mode == -1) {
		rmatmul_impl_neon<FE_DOWNWARD>(m, n, k, A, B, CU);
	}
	else {
		rmatmul_impl_neon<FE_TONEAREST>(m, n, k, A, B, CU);
	}
	vblas_rmatmul_neon_detail::reset_rounding_mode();
}

#endif // VBLAS_RMATMUL_NEON_HPP
