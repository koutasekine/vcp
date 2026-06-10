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

#ifndef VBLAS_RMATMUL_AVX2_HPP
#define VBLAS_RMATMUL_AVX2_HPP

#include <algorithm>
#include <cfenv>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <immintrin.h>
#include <omp.h>

#pragma STDC FENV_ACCESS ON

#if !defined(__AVX2__) || !defined(__FMA__)
#error "rmatmul_avx2.hpp requires AVX2 and FMA compiler support."
#endif

namespace vblas_rmatmul_avx2_detail {

constexpr int VEC = 4;
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
	int threads_per_core;
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

// "0-1" や "0,8" 形式の cpu list から CPU 数を数える
inline int parse_cpu_list_count(const std::string& text) {
	int count = 0;
	std::string::size_type pos = 0;
	while (pos < text.size()) {
		const std::string::size_type comma = text.find(',', pos);
		const std::string token = text.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
		const std::string::size_type dash = token.find('-');
		if (dash == std::string::npos) {
			if (!token.empty()) {
				count++;
			}
		}
		else {
			const int first = std::atoi(token.substr(0, dash).c_str());
			const int last = std::atoi(token.substr(dash + 1).c_str());
			count += std::max(0, last - first + 1);
		}
		if (comma == std::string::npos) {
			break;
		}
		pos = comma + 1;
	}
	return count > 0 ? count : 1;
}

inline CpuInfo detect_cpu_info() {
	CpuInfo info = {0, 0, 0, std::max(1, omp_get_num_procs()), std::max(1, omp_get_max_threads()), 1};

	std::string siblings;
	if (read_first_line("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list", siblings)) {
		info.threads_per_core = parse_cpu_list_count(siblings);
	}

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

// B micro-panel (kc*NR) が L1d の約 1/2 に収まるように kc を決める．
// kc が大きいほど C への蓄積パス回数 (n/kc) が減り C のメモリ往復が減る．
inline int cache_based_kc(const CpuInfo& info) {
	const long long l1 = info.l1d_bytes > 0 ? info.l1d_bytes : 32LL * 1024LL;
	const long long budget = l1 / 2LL;
	const long long per_k = static_cast<long long>(NR) * MATRIX_BYTES;
	return clamp_rounded_block(static_cast<int>(budget / per_k), 64, 512, NR);
}

// packed A block (mc*kc) が L2 の約 1/2 に収まるように mc を決める
// (残りは B column tile と C の通過分に充てる)
inline int cache_based_mc(const CpuInfo& info, const int kc) {
	const long long l2 = info.l2_bytes > 0 ? info.l2_bytes : 256LL * 1024LL;
	// L1/L2 は同一 core 上の hyper-thread で共有されるため，core あたりの
	// thread 数で budget を分割する
	const long long budget_doubles = l2 / 2LL / MATRIX_BYTES / std::max(1, info.threads_per_core);
	return clamp_rounded_block(static_cast<int>(budget_doubles / kc), MR, 1024, MR);
}

// B column tile (kc*column_tile) が L2 の約 1/4 に収まるように決める
inline int cache_based_column_tile(const CpuInfo& info, const int kc) {
	const long long l2 = info.l2_bytes > 0 ? info.l2_bytes : 256LL * 1024LL;
	const long long budget_doubles = l2 / 4LL / MATRIX_BYTES;
	return clamp_rounded_block(static_cast<int>(budget_doubles / kc), NR, 512, NR);
}

// packed B group (kc*nc_group, 全 thread で共有) が L3 の約 1/2 に収まるように決める
inline int cache_based_nc_group(const CpuInfo& info, const int kc, const int column_tile) {
	const long long l3 = info.l3_bytes > 0 ? info.l3_bytes : 8LL * 1024LL * 1024LL;
	const long long budget_doubles = l3 / 2LL / MATRIX_BYTES;
	const int raw = static_cast<int>(budget_doubles / kc);
	return clamp_rounded_block(std::max(column_tile, raw), 128, 8192, NR);
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

	// cache size と thread 数のみから block size を導出する (環境特化の分岐は置かない)
	const int kc = cache_based_kc(info);
	const int mc = cache_based_mc(info, kc);
	const int column_tile = cache_based_column_tile(info, kc);
	const int nc_group = cache_based_nc_group(info, kc, column_tile);

	// m 方向の block 数が thread 数以上なら m 方向のみの並列で負荷分散できる
	const ParallelMode mode = ceil_div(m, std::max(MR, mc)) >= std::max(1, threads) ? PARALLEL_M_ONLY : PARALLEL_2D;
	return finish_blocking(Blocking{mc, kc, nc_group, column_tile, threads, mode, false}, m, k);
}

inline __m256i row_mask(const int rows) {
	const long long on = -1LL;
	const long long off = 0LL;
	return _mm256_set_epi64x(
		rows >= 4 ? on : off,
		rows >= 3 ? on : off,
		rows >= 2 ? on : off,
		rows >= 1 ? on : off
	);
}

inline void set_rounding_mode(const int rounding) {
	std::fesetround(rounding);
}

inline int get_rounding_mode() {
	return std::fegetround();
}

inline void reset_rounding_mode(const int rounding) {
	std::fesetround(rounding);
}

// A block (mb x pb) を MR 行の micro-panel 単位で連続配置に pack する．
// layout: packed[(ir/MR)*(MR*pb) + p*MR + i]
// kernel が A を stride MR の完全連続 load で読めるようにし，
// 大きな leading dimension による L1 set 競合を避ける．端数行は 0 詰め．
inline void pack_a_micro(
	const int m,
	const int mb,
	const int pb,
	const double* A,
	double* packed
) {
	const int full_panels = mb / MR;
	for (int q = 0; q < full_panels; q++) {
		double* dst = packed + static_cast<std::size_t>(q) * MR * pb;
		const double* src = A + static_cast<std::size_t>(q) * MR;
		for (int p = 0; p < pb; p++) {
			const double* s = src + static_cast<std::size_t>(m) * p;
			double* d = dst + static_cast<std::size_t>(p) * MR;
			for (int i = 0; i < MR; i++) {
				d[i] = s[i];
			}
		}
	}

	const int ir = full_panels * MR;
	const int rows = mb - ir;
	if (rows > 0) {
		double* dst = packed + static_cast<std::size_t>(full_panels) * MR * pb;
		const double* src = A + ir;
		for (int p = 0; p < pb; p++) {
			const double* s = src + static_cast<std::size_t>(m) * p;
			double* d = dst + static_cast<std::size_t>(p) * MR;
			for (int i = 0; i < rows; i++) {
				d[i] = s[i];
			}
			for (int i = rows; i < MR; i++) {
				d[i] = 0.0;
			}
		}
	}
}

// B panel (pb x nb) を NR 列の micro-panel 単位で連続配置に pack する．
// layout: packed[(jr/NR)*(pb*NR) + p*NR + j]
// kernel が B を stride NR の単一連続 stream で読めるようにする．端数列は 0 詰め．
inline void pack_b_micro_parallel(
	const int n,
	const int pb,
	const int nb,
	const double* B,
	double* packed
) {
	const int panels = ceil_div(nb, NR);
#pragma omp for schedule(static)
	for (int jp = 0; jp < panels; jp++) {
		const int jr = jp * NR;
		const int cols = std::min(NR, nb - jr);
		double* dst = packed + static_cast<std::size_t>(jp) * pb * NR;
		const double* src = B + static_cast<std::size_t>(n) * jr;
		if (cols == NR) {
			for (int p = 0; p < pb; p++) {
				double* d = dst + static_cast<std::size_t>(p) * NR;
				for (int j = 0; j < NR; j++) {
					d[j] = src[p + static_cast<std::size_t>(n) * j];
				}
			}
		}
		else {
			for (int p = 0; p < pb; p++) {
				double* d = dst + static_cast<std::size_t>(p) * NR;
				for (int j = 0; j < cols; j++) {
					d[j] = src[p + static_cast<std::size_t>(n) * j];
				}
				for (int j = cols; j < NR; j++) {
					d[j] = 0.0;
				}
			}
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
	const __m256d zero = _mm256_setzero_pd();
	__m256d cu[NR][MR_CHUNKS];

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = init_zero ? zero : _mm256_loadu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r);
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		__m256d a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = _mm256_loadu_pd(ap + VEC * r);
		}
		for (int j = 0; j < NR; j++) {
			const __m256d b = _mm256_set1_pd(B[p + static_cast<std::size_t>(ldb) * j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = _mm256_fmadd_pd(a[r], b, cu[j][r]);
			}
		}
	}

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			_mm256_storeu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, cu[j][r]);
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
	const __m256d zero = _mm256_setzero_pd();
	__m256d cu[NR][MR_CHUNKS];
	__m256i mask[MR_CHUNKS];

	for (int r = 0; r < MR_CHUNKS; r++) {
		mask[r] = row_mask(std::min(VEC, std::max(0, mr - VEC * r)));
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = init_zero ? zero : _mm256_maskload_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, mask[r]);
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		__m256d a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = _mm256_maskload_pd(ap + VEC * r, mask[r]);
		}
		for (int j = 0; j < nr; j++) {
			const __m256d b = _mm256_set1_pd(B[p + static_cast<std::size_t>(ldb) * j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = _mm256_fmadd_pd(a[r], b, cu[j][r]);
			}
		}
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			_mm256_maskstore_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, mask[r], cu[j][r]);
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
	const __m256d zero = _mm256_setzero_pd();
	__m256d cu[NR][MR_CHUNKS];

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = init_zero ? zero : _mm256_loadu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r);
		}
	}

	for (int p = 0; p < pb; p++) {
		__m256d a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = _mm256_loadu_pd(A + VEC * r);
		}
		for (int j = 0; j < NR; j++) {
			const __m256d b = _mm256_set1_pd(B[j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = _mm256_fmadd_pd(a[r], b, cu[j][r]);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			_mm256_storeu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, cu[j][r]);
		}
	}
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
	const __m256d zero = _mm256_setzero_pd();
	__m256d cu[NR][MR_CHUNKS];
	__m256i mask[MR_CHUNKS];

	for (int r = 0; r < MR_CHUNKS; r++) {
		mask[r] = row_mask(std::min(VEC, std::max(0, mr - VEC * r)));
	}

	for (int j = 0; j < NR; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			cu[j][r] = (init_zero || j >= nr) ? zero : _mm256_maskload_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, mask[r]);
		}
	}

	for (int p = 0; p < pb; p++) {
		__m256d a[MR_CHUNKS];
		for (int r = 0; r < MR_CHUNKS; r++) {
			a[r] = _mm256_loadu_pd(A + VEC * r);
		}
		for (int j = 0; j < NR; j++) {
			const __m256d b = _mm256_set1_pd(B[j]);
			for (int r = 0; r < MR_CHUNKS; r++) {
				cu[j][r] = _mm256_fmadd_pd(a[r], b, cu[j][r]);
			}
		}
		A += MR;
		B += NR;
	}

	for (int j = 0; j < nr; j++) {
		for (int r = 0; r < MR_CHUNKS; r++) {
			_mm256_maskstore_pd(CU + static_cast<std::size_t>(ldc) * j + VEC * r, mask[r], cu[j][r]);
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
inline void rmatmul_direct_avx2(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
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
inline void rmatmul_blocked_m_only_avx2(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
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
inline void rmatmul_blocked_2d_avx2(const int m, const int n, const int k, const double* A, const double* B, double* CU, const Blocking& blocking) {
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

} // namespace vblas_rmatmul_avx2_detail

template<int ROUNDING>
inline void rmatmul_impl_avx2(int m, int n, int k, const double* A, const double* B, double* CU) {
	using namespace vblas_rmatmul_avx2_detail;

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
		rmatmul_direct_avx2<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else if (blocking.mode == PARALLEL_M_ONLY) {
		rmatmul_blocked_m_only_avx2<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
	else {
		rmatmul_blocked_2d_avx2<ROUNDING>(m, n, k, A, B, CU, blocking);
	}
}

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU is overwritten with A*B using rounding_mode: 1 upward, 0 nearest, -1 downward.
inline void rmatmul_avx2(int m, int n, int k, const double* A, const double* B, double* CU, const int rounding_mode) {
	const int previous_rounding = vblas_rmatmul_avx2_detail::get_rounding_mode();
	if (rounding_mode == 1) {
		rmatmul_impl_avx2<FE_UPWARD>(m, n, k, A, B, CU);
	}
	else if (rounding_mode == -1) {
		rmatmul_impl_avx2<FE_DOWNWARD>(m, n, k, A, B, CU);
	}
	else {
		rmatmul_impl_avx2<FE_TONEAREST>(m, n, k, A, B, CU);
	}
	vblas_rmatmul_avx2_detail::reset_rounding_mode(previous_rounding);
}

#endif // VBLAS_RMATMUL_AVX2_HPP
