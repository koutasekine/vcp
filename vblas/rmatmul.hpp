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

#ifndef VBLAS_RMATMUL_HPP
#define VBLAS_RMATMUL_HPP

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <immintrin.h>
#include <omp.h>

namespace vblas_rmatmul_detail {

constexpr int VEC = 8;
constexpr int MR = 16;
constexpr int NR = 4;
constexpr int MR_CHUNKS = MR / VEC;
constexpr int MATRIX_BYTES = 8;

enum ParallelMode {
	PARALLEL_2D,
	PARALLEL_M_ONLY
};

struct CpuInfo {
	long long l1d_bytes;
	long long l2_bytes;
	long long l3_bytes;
	int logical_threads;
	int omp_max_threads;
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

inline int clamp_rounded_block(const int value, const int min_value, const int max_value, const int block) {
	const int rounded_min = std::max(block, round_up_to(min_value, block));
	const int rounded_max = std::max(rounded_min, round_down_to(max_value, block));
	int rounded = round_down_to(value, block);
	if (rounded < rounded_min) {
		rounded = rounded_min;
	}
	if (rounded > rounded_max) {
		rounded = rounded_max;
	}
	return rounded;
}

inline std::string trim_copy(const std::string& value) {
	std::string::size_type first = 0;
	while (first < value.size() && (value[first] == ' ' || value[first] == '\t' || value[first] == '\n' || value[first] == '\r')) {
		first++;
	}
	std::string::size_type last = value.size();
	while (last > first && (value[last - 1] == ' ' || value[last - 1] == '\t' || value[last - 1] == '\n' || value[last - 1] == '\r')) {
		last--;
	}
	return value.substr(first, last - first);
}

inline bool read_first_line(const std::string& path, std::string& line) {
	std::ifstream file(path.c_str());
	if (!file) {
		return false;
	}
	std::getline(file, line);
	line = trim_copy(line);
	return true;
}

inline long long parse_cache_size(const std::string& value) {
	const std::string text = trim_copy(value);
	char* end = NULL;
	const long long number = std::strtoll(text.c_str(), &end, 10);
	if (end == text.c_str()) {
		return 0;
	}
	while (*end == ' ' || *end == '\t') {
		end++;
	}
	if (*end == 'K' || *end == 'k') {
		return number * 1024LL;
	}
	if (*end == 'M' || *end == 'm') {
		return number * 1024LL * 1024LL;
	}
	return number;
}

inline CpuInfo detect_cpu_info() {
	CpuInfo info = {0, 0, 0, std::max(1, omp_get_num_procs()), std::max(1, omp_get_max_threads())};

	for (int index = 0; index < 16; index++) {
		const std::string base = "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(index) + "/";
		std::string level_text;
		std::string type_text;
		std::string size_text;
		if (!read_first_line(base + "level", level_text) ||
		    !read_first_line(base + "type", type_text) ||
		    !read_first_line(base + "size", size_text)) {
			continue;
		}

		const int level = std::atoi(level_text.c_str());
		const long long bytes = parse_cache_size(size_text);
		if (level == 1 && type_text == "Data") {
			info.l1d_bytes = bytes;
		}
		else if (level == 2) {
			info.l2_bytes = bytes;
		}
		else if (level == 3) {
			info.l3_bytes = std::max(info.l3_bytes, bytes);
		}
	}

	return info;
}

inline const CpuInfo& cpu_info() {
	static const CpuInfo info = detect_cpu_info();
	return info;
}

inline int limited_threads(const int requested) {
	return std::max(1, std::min(requested, std::max(1, omp_get_max_threads())));
}

inline int default_thread_count(const CpuInfo& info, const int m, const int n, const int k) {
	const int max_threads = limited_threads(std::min(info.logical_threads, info.omp_max_threads));
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

inline int cache_based_kc(const CpuInfo& info, const int fallback) {
	if (info.l1d_bytes <= 0) {
		return fallback;
	}

	const long long budget = info.l1d_bytes * 3LL / 4LL;
	const long long fixed_c = static_cast<long long>(MR) * NR * MATRIX_BYTES;
	const long long per_k = static_cast<long long>(MR + NR) * MATRIX_BYTES;
	const long long raw = budget > fixed_c ? (budget - fixed_c) / per_k : fallback;
	return clamp_rounded_block(static_cast<int>(raw), 64, 512, NR);
}

inline int cache_based_column_tile(const CpuInfo& info, const int kc, const int fallback) {
	if (info.l2_bytes <= 0 || kc <= 0) {
		return fallback;
	}

	const int seed_mc = 128;
	const long long budget_doubles = (info.l2_bytes * 3LL / 4LL) / MATRIX_BYTES;
	const long long fixed_a = static_cast<long long>(seed_mc) * kc;
	const long long raw = budget_doubles > fixed_a ? (budget_doubles - fixed_a) / (kc + seed_mc) : fallback;
	return clamp_rounded_block(static_cast<int>(raw), 32, 512, NR);
}

inline int cache_based_mc(const CpuInfo& info, const int kc, const int column_tile, const int fallback) {
	if (info.l2_bytes <= 0 || kc <= 0 || column_tile <= 0) {
		return fallback;
	}

	const long long budget_doubles = (info.l2_bytes * 3LL / 4LL) / MATRIX_BYTES;
	const long long fixed_b = static_cast<long long>(kc) * column_tile;
	const long long raw = budget_doubles > fixed_b ? (budget_doubles - fixed_b) / (kc + column_tile) : fallback;
	return clamp_rounded_block(static_cast<int>(raw), 64, 512, MR);
}

inline int cache_based_nc_group(const CpuInfo& info, const int kc, const int mc, const int threads, const int fallback) {
	if (info.l3_bytes <= 0 || kc <= 0 || mc <= 0) {
		return fallback;
	}

	const int active_threads = std::max(1, threads);
	const long long budget_doubles = (info.l3_bytes * 2LL / 3LL) / MATRIX_BYTES;
	const long long active_a = static_cast<long long>(active_threads) * mc * kc;
	const long long per_column = kc;
	const long long raw = budget_doubles > active_a ? (budget_doubles - active_a) / per_column : fallback;
	return clamp_rounded_block(static_cast<int>(raw), 128, 4096, NR);
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
	const CpuInfo& info = cpu_info();
	const int threads = default_thread_count(info, m, n, k);

	if (use_direct_path(m, n, k)) {
		const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
		const int direct_threads = work >= 1.0e7 ? threads : std::min(threads, 8);
		return finish_blocking(Blocking{MR, std::max(NR, n), std::max(NR, k), NR, direct_threads, PARALLEL_2D, true}, m, k);
	}

	const int kc = cache_based_kc(info, 160);
	int column_tile = cache_based_column_tile(info, kc, 128);
	int mc = cache_based_mc(info, kc, column_tile, 128);
	int nc_group = cache_based_nc_group(info, kc, mc, threads, 1024);
	ParallelMode mode = ceil_div(m, std::max(MR, mc)) >= std::max(1, threads / 2) ? PARALLEL_M_ONLY : PARALLEL_2D;

	const int max_nk = std::max(n, k);
	const int max_mn = std::max(m, n);
	const int max_dim = std::max(m, max_nk);

	if (m >= 4096 && m >= 4 * max_nk) {
		mc = std::max(mc, 128);
		column_tile = std::max(NR * 8, std::min(column_tile, 256));
		mode = PARALLEL_M_ONLY;
	}
	else if (k >= 4096 && k >= 4 * max_mn) {
		mc = std::min(mc, 160);
		column_tile = std::min(column_tile, 128);
		mode = PARALLEL_2D;
	}
	else if (n >= 4096 && n >= 4 * std::max(m, k)) {
		mc = std::min(mc, 160);
		column_tile = std::min(column_tile, 128);
		mode = PARALLEL_2D;
	}
	else if (max_dim <= 1500) {
		mc = std::min(mc, 128);
		column_tile = std::min(column_tile, 128);
		mode = ceil_div(m, std::max(MR, mc)) >= std::max(1, threads / 2) ? PARALLEL_M_ONLY : PARALLEL_2D;
	}

	nc_group = std::max(column_tile, nc_group);
	return finish_blocking(Blocking{mc, kc, nc_group, column_tile, threads, mode, false}, m, k);
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

template<int ROUNDING>
inline void compute_block(
	const int m,
	const int n,
	const int pb,
	const int ic,
	const int mb,
	const int pc,
	const int jc,
	const int nb,
	const double* A,
	const double* B,
	double* CU,
	const bool init_zero
) {
	for (int jr = 0; jr < nb; jr += NR) {
		const int nr = std::min(NR, nb - jr);
		for (int ir = 0; ir < mb; ir += MR) {
			const int mr = std::min(MR, mb - ir);
			const double* ap = A + (ic + ir) + static_cast<std::size_t>(m) * pc;
			const double* bp = B + pc + static_cast<std::size_t>(n) * (jc + jr);
			double* cup = CU + (ic + ir) + static_cast<std::size_t>(m) * (jc + jr);
			if (mr == MR && nr == NR) {
				kernel_mrx4<ROUNDING>(pb, ap, m, bp, n, cup, m, init_zero);
			}
			else {
				kernel_mrx4_masked<ROUNDING>(pb, mr, nr, ap, m, bp, n, cup, m, init_zero);
			}
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
inline void rmatmul_direct(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
#pragma omp parallel num_threads(blocking.threads)
	{
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
	}
}

template<int ROUNDING>
inline void rmatmul_blocked_m_only(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * std::min(blocking.nc_group, k));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);

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
	}
}

template<int ROUNDING>
inline void rmatmul_blocked_2d(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * std::min(blocking.nc_group, k));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);

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
	}
}

} // namespace vblas_rmatmul_detail

template<int ROUNDING>
inline void rmatmul_impl(int m, int n, int k, const double* A, const double* B, double* CU) {
	using namespace vblas_rmatmul_detail;

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
		rmatmul_direct<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else if (blocking.mode == PARALLEL_M_ONLY) {
		rmatmul_blocked_m_only<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else {
		rmatmul_blocked_2d<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
}

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU is overwritten with A*B using rounding_mode: 1 upward, 0 nearest, -1 downward.
inline void rmatmul(int m, int n, int k, const double* A, const double* B, double* CU, const int rounding_mode) {
	if (rounding_mode == 1) {
		rmatmul_impl<_MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC>(m, n, k, A, B, CU);
	}
	else if (rounding_mode == -1) {
		rmatmul_impl<_MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC>(m, n, k, A, B, CU);
	}
	else {
		rmatmul_impl<_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC>(m, n, k, A, B, CU);
	}
}

#endif // VBLAS_RMATMUL_HPP
