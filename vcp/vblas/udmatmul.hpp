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

#ifndef VBLAS_UDMATMUL_HPP
#define VBLAS_UDMATMUL_HPP

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

#include <immintrin.h>
#include <omp.h>


namespace vblas_udmatmul_detail {

constexpr int VEC = 8;
constexpr int MR = 16;
constexpr int NR = 4;

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
	std::string model_name;
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

inline std::string read_model_name() {
	std::ifstream file("/proc/cpuinfo");
	std::string line;
	while (std::getline(file, line)) {
		const std::string key = "model name";
		if (line.compare(0, key.size(), key) == 0) {
			const std::string::size_type colon = line.find(':');
			if (colon != std::string::npos) {
				return trim_copy(line.substr(colon + 1));
			}
		}
	}
	return std::string();
}

inline CpuInfo detect_cpu_info() {
	CpuInfo info = {0, 0, 0, std::max(1, omp_get_num_procs()), std::max(1, omp_get_max_threads()), read_model_name()};

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

inline int default_thread_cap(const CpuInfo& info) {
	if (info.logical_threads <= 8) {
		return std::max(1, info.logical_threads);
	}
	return 16;
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
	const int kc_raw = static_cast<int>(info.l1d_bytes / (MR * 8LL));
	const int kc = (std::max(NR, kc_raw) / NR) * NR;
	return std::max(64, std::min(256, kc));
}

inline int cache_based_column_tile(const CpuInfo& info, const int kc, const int fallback) {
	if (info.l2_bytes <= 0 || kc <= 0) {
		return fallback;
	}
	const int ct_raw = static_cast<int>(info.l2_bytes / (4LL * kc * 8LL));
	const int ct = (std::max(NR, ct_raw) / NR) * NR;
	return std::max(NR, std::min(256, ct));
}

inline bool use_direct_path(const int m, const int n, const int k) {
	const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
	if (work <= 2.5e6) {
		return true;
	}
	if (k <= 64 && m <= 8192 && n <= 8192) {
		return true;
	}
	return static_cast<double>(m) * static_cast<double>(n) <= 32768.0 && n <= 512;
}

inline bool prefer_backend(const int m, const int n, const int k) {
	const int max_dim = std::max(m, std::max(n, k));
	return max_dim >= 8000 && !use_direct_path(m, n, k);
}

inline Blocking select_blocking(const int m, const int n, const int k) {
	const CpuInfo& info = cpu_info();
	const int t = limited_threads(default_thread_cap(info));
	const int max_dim = std::max(m, std::max(n, k));
	const int max_nk = std::max(n, k);
	const int max_mn = std::max(m, n);

	if (use_direct_path(m, n, k)) {
		const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
		const int direct_threads = work >= 1.0e7 ? t : std::min(t, 8);
		return finish_blocking({MR, n, 0, NR, direct_threads, PARALLEL_2D, true}, m, k);
	}

	const int kc = cache_based_kc(info, 192);
	const int ct = cache_based_column_tile(info, kc, 96);
	const int nc_group = info.l3_bytes >= 8LL * 1024LL * 1024LL ? 2048 : 1024;

	if (m >= 4096 && m >= 4 * max_nk) {
		return finish_blocking({128, kc, nc_group, std::max(64, ct), t, PARALLEL_M_ONLY, false}, m, k);
	}
	if (k >= 4096 && k >= 4 * max_mn) {
		return finish_blocking({64, kc, std::min(nc_group, 1024), std::min(ct, 96), t, PARALLEL_2D, false}, m, k);
	}
	if (n >= 4096 && n >= 4 * std::max(m, k)) {
		return finish_blocking({64, kc, nc_group, std::min(ct, 96), t, PARALLEL_2D, false}, m, k);
	}
	if (max_dim <= 640) {
		return finish_blocking({64, std::min(kc, 192), nc_group, std::min(ct, 96), t, PARALLEL_2D, false}, m, k);
	}
	if (max_dim <= 1500) {
		return finish_blocking({64, std::min(kc, 192), nc_group, 128, t, PARALLEL_2D, false}, m, k);
	}

	const ParallelMode mode = ceil_div(m, 128) >= std::max(1, t / 2) ? PARALLEL_M_ONLY : PARALLEL_2D;
	return finish_blocking({128, kc, nc_group, std::max(64, std::min(ct, 128)), t, mode, false}, m, k);
}

inline const char* mode_name(const ParallelMode mode) {
	return mode == PARALLEL_M_ONLY ? "m_only" : "2d";
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

inline void pack_b_panel(
	const int n,
	const int pb,
	const int nb,
	const double* B,
	double* packed
) {
	for (int j = 0; j < nb; j++) {
		double* dst = packed + static_cast<std::size_t>(pb) * j;
		const double* src = B + static_cast<std::size_t>(n) * j;
		for (int p = 0; p < pb; p++) {
			dst[p] = src[p];
		}
	}
}

inline void kernel_16x4(
	const int pb,
	const double* A,
	const int lda,
	const double* B,
	const int ldb,
	double* CU,
	double* CD,
	const int ldc
) {
	__m512d cu0l = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 0);
	__m512d cu0h = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 0 + VEC);
	__m512d cu1l = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 1);
	__m512d cu1h = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 1 + VEC);
	__m512d cu2l = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 2);
	__m512d cu2h = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 2 + VEC);
	__m512d cu3l = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 3);
	__m512d cu3h = _mm512_loadu_pd(CU + static_cast<std::size_t>(ldc) * 3 + VEC);

	__m512d cd0l = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 0);
	__m512d cd0h = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 0 + VEC);
	__m512d cd1l = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 1);
	__m512d cd1h = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 1 + VEC);
	__m512d cd2l = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 2);
	__m512d cd2h = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 2 + VEC);
	__m512d cd3l = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 3);
	__m512d cd3h = _mm512_loadu_pd(CD + static_cast<std::size_t>(ldc) * 3 + VEC);

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		const __m512d a0 = _mm512_loadu_pd(ap);
		const __m512d a1 = _mm512_loadu_pd(ap + VEC);
		const __m512d b0 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 0]);
		const __m512d b1 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 1]);
		const __m512d b2 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 2]);
		const __m512d b3 = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * 3]);

		cu0l = _mm512_fmadd_round_pd(a0, b0, cu0l, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu0h = _mm512_fmadd_round_pd(a1, b0, cu0h, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu1l = _mm512_fmadd_round_pd(a0, b1, cu1l, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu1h = _mm512_fmadd_round_pd(a1, b1, cu1h, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu2l = _mm512_fmadd_round_pd(a0, b2, cu2l, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu2h = _mm512_fmadd_round_pd(a1, b2, cu2h, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu3l = _mm512_fmadd_round_pd(a0, b3, cu3l, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
		cu3h = _mm512_fmadd_round_pd(a1, b3, cu3h, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);

		cd0l = _mm512_fmadd_round_pd(a0, b0, cd0l, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd0h = _mm512_fmadd_round_pd(a1, b0, cd0h, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd1l = _mm512_fmadd_round_pd(a0, b1, cd1l, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd1h = _mm512_fmadd_round_pd(a1, b1, cd1h, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd2l = _mm512_fmadd_round_pd(a0, b2, cd2l, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd2h = _mm512_fmadd_round_pd(a1, b2, cd2h, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd3l = _mm512_fmadd_round_pd(a0, b3, cd3l, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
		cd3h = _mm512_fmadd_round_pd(a1, b3, cd3h, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
	}

	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 0, cu0l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 0 + VEC, cu0h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 1, cu1l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 1 + VEC, cu1h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 2, cu2l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 2 + VEC, cu2h);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 3, cu3l);
	_mm512_storeu_pd(CU + static_cast<std::size_t>(ldc) * 3 + VEC, cu3h);

	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 0, cd0l);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 0 + VEC, cd0h);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 1, cd1l);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 1 + VEC, cd1h);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 2, cd2l);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 2 + VEC, cd2h);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 3, cd3l);
	_mm512_storeu_pd(CD + static_cast<std::size_t>(ldc) * 3 + VEC, cd3h);
}

inline void kernel_16x4_packed_masked(
	const int pb,
	const int mr,
	const int nr,
	const double* A,
	const int lda,
	const double* B,
	const int ldb,
	double* CU,
	double* CD,
	const int ldc
) {
	const int rows_lo = std::min(mr, VEC);
	const int rows_hi = std::max(0, mr - VEC);
	const __mmask8 mask_lo = row_mask(rows_lo);
	const __mmask8 mask_hi = row_mask(rows_hi);
	__m512d cu_lo[NR];
	__m512d cu_hi[NR];
	__m512d cd_lo[NR];
	__m512d cd_hi[NR];

	for (int j = 0; j < nr; j++) {
		cu_lo[j] = _mm512_maskz_loadu_pd(mask_lo, CU + static_cast<std::size_t>(ldc) * j);
		cd_lo[j] = _mm512_maskz_loadu_pd(mask_lo, CD + static_cast<std::size_t>(ldc) * j);
		if (rows_hi > 0) {
			cu_hi[j] = _mm512_maskz_loadu_pd(mask_hi, CU + static_cast<std::size_t>(ldc) * j + VEC);
			cd_hi[j] = _mm512_maskz_loadu_pd(mask_hi, CD + static_cast<std::size_t>(ldc) * j + VEC);
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		const __m512d a0 = _mm512_loadu_pd(ap);
		const __m512d a1 = rows_hi > 0 ? _mm512_loadu_pd(ap + VEC) : _mm512_setzero_pd();
		for (int j = 0; j < nr; j++) {
			const __m512d b = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * j]);
			cu_lo[j] = _mm512_fmadd_round_pd(a0, b, cu_lo[j], _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
			cd_lo[j] = _mm512_fmadd_round_pd(a0, b, cd_lo[j], _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			if (rows_hi > 0) {
				cu_hi[j] = _mm512_fmadd_round_pd(a1, b, cu_hi[j], _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				cd_hi[j] = _mm512_fmadd_round_pd(a1, b, cd_hi[j], _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			}
		}
	}

	for (int j = 0; j < nr; j++) {
		_mm512_mask_storeu_pd(CU + static_cast<std::size_t>(ldc) * j, mask_lo, cu_lo[j]);
		_mm512_mask_storeu_pd(CD + static_cast<std::size_t>(ldc) * j, mask_lo, cd_lo[j]);
		if (rows_hi > 0) {
			_mm512_mask_storeu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC, mask_hi, cu_hi[j]);
			_mm512_mask_storeu_pd(CD + static_cast<std::size_t>(ldc) * j + VEC, mask_hi, cd_hi[j]);
		}
	}
}

inline void kernel_16x4_direct_masked(
	const int pb,
	const int mr,
	const int nr,
	const double* A,
	const int lda,
	const double* B,
	const int ldb,
	double* CU,
	double* CD,
	const int ldc
) {
	const int rows_lo = std::min(mr, VEC);
	const int rows_hi = std::max(0, mr - VEC);
	const __mmask8 mask_lo = row_mask(rows_lo);
	const __mmask8 mask_hi = row_mask(rows_hi);
	__m512d cu_lo[NR];
	__m512d cu_hi[NR];
	__m512d cd_lo[NR];
	__m512d cd_hi[NR];

	for (int j = 0; j < nr; j++) {
		cu_lo[j] = _mm512_maskz_loadu_pd(mask_lo, CU + static_cast<std::size_t>(ldc) * j);
		cd_lo[j] = _mm512_maskz_loadu_pd(mask_lo, CD + static_cast<std::size_t>(ldc) * j);
		if (rows_hi > 0) {
			cu_hi[j] = _mm512_maskz_loadu_pd(mask_hi, CU + static_cast<std::size_t>(ldc) * j + VEC);
			cd_hi[j] = _mm512_maskz_loadu_pd(mask_hi, CD + static_cast<std::size_t>(ldc) * j + VEC);
		}
	}

	for (int p = 0; p < pb; p++) {
		const double* ap = A + static_cast<std::size_t>(p) * lda;
		const __m512d a0 = _mm512_maskz_loadu_pd(mask_lo, ap);
		const __m512d a1 = rows_hi > 0 ? _mm512_maskz_loadu_pd(mask_hi, ap + VEC) : _mm512_setzero_pd();
		for (int j = 0; j < nr; j++) {
			const __m512d b = _mm512_set1_pd(B[p + static_cast<std::size_t>(ldb) * j]);
			cu_lo[j] = _mm512_fmadd_round_pd(a0, b, cu_lo[j], _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
			cd_lo[j] = _mm512_fmadd_round_pd(a0, b, cd_lo[j], _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			if (rows_hi > 0) {
				cu_hi[j] = _mm512_fmadd_round_pd(a1, b, cu_hi[j], _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
				cd_hi[j] = _mm512_fmadd_round_pd(a1, b, cd_hi[j], _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
			}
		}
	}

	for (int j = 0; j < nr; j++) {
		_mm512_mask_storeu_pd(CU + static_cast<std::size_t>(ldc) * j, mask_lo, cu_lo[j]);
		_mm512_mask_storeu_pd(CD + static_cast<std::size_t>(ldc) * j, mask_lo, cd_lo[j]);
		if (rows_hi > 0) {
			_mm512_mask_storeu_pd(CU + static_cast<std::size_t>(ldc) * j + VEC, mask_hi, cu_hi[j]);
			_mm512_mask_storeu_pd(CD + static_cast<std::size_t>(ldc) * j + VEC, mask_hi, cd_hi[j]);
		}
	}
}

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
	double* CD
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
				double* cdp = CD + (ic + ir) + static_cast<std::size_t>(m) * (jcg + jc_offset + jr);
				if (mr == MR && nr == NR) {
					kernel_16x4(pb, ap, mb_padded, bp, pb, cup, cdp, m);
				}
				else {
					kernel_16x4_packed_masked(pb, mr, nr, ap, mb_padded, bp, pb, cup, cdp, m);
				}
			}
		}
	}
}

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
	double* CD
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
			double* cdp = CD + (ic + ir) + static_cast<std::size_t>(m) * (j_base + jr);
			if (mr == MR && nr == NR) {
				kernel_16x4(pb, ap, mb_padded, bp, pb, cup, cdp, m);
			}
			else {
				kernel_16x4_packed_masked(pb, mr, nr, ap, mb_padded, bp, pb, cup, cdp, m);
			}
		}
	}
}

inline void udmatmul_direct(const int m, const int n, const int k, double* A, double* B, double* CU, double* CD, const Blocking& blocking) {
#pragma omp parallel for collapse(2) schedule(static) num_threads(blocking.threads)
	for (int jr = 0; jr < k; jr += NR) {
		for (int ir = 0; ir < m; ir += MR) {
			const int nr = std::min(NR, k - jr);
			const int mr = std::min(MR, m - ir);
			double* cup = CU + ir + static_cast<std::size_t>(m) * jr;
			double* cdp = CD + ir + static_cast<std::size_t>(m) * jr;
			const double* ap = A + ir;
			const double* bp = B + static_cast<std::size_t>(n) * jr;
			if (mr == MR && nr == NR) {
				kernel_16x4(n, ap, m, bp, n, cup, cdp, m);
			}
			else {
				kernel_16x4_direct_masked(n, mr, nr, ap, m, bp, n, cup, cdp, m);
			}
		}
	}
}

inline void udmatmul_packed_m_only(const int m, const int n, const int k, double* A, double* B, double* CU, double* CD, const Blocking& blocking) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * std::min(blocking.nc_group, k));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);

		for (int pc = 0; pc < n; pc += blocking.kc) {
			const int pb = std::min(blocking.kc, n - pc);

			for (int jcg = 0; jcg < k; jcg += blocking.nc_group) {
				const int ng = std::min(blocking.nc_group, k - jcg);

#pragma omp single
				{
					packed_b.resize(static_cast<std::size_t>(pb) * ng);
					pack_b_panel(n, pb, ng, B + pc + static_cast<std::size_t>(n) * jcg, packed_b.data());
				}

#pragma omp for schedule(static)
				for (int ic = 0; ic < m; ic += blocking.mc) {
					compute_packed_ic_block(m, pb, ng, ic, pc, jcg, blocking, A, packed_b.data(), packed_a.data(), CU, CD);
				}
			}
		}
	}
}

inline void udmatmul_packed_2d(const int m, const int n, const int k, double* A, double* B, double* CU, double* CD, const Blocking& blocking) {
	const int max_mb_padded = round_up_to(std::min(blocking.mc, m), MR);
	const int max_pb = std::min(blocking.kc, n);
	std::vector<double> packed_b;
	packed_b.reserve(static_cast<std::size_t>(max_pb) * std::min(blocking.nc_group, k));

#pragma omp parallel num_threads(blocking.threads)
	{
		std::vector<double> packed_a(static_cast<std::size_t>(max_mb_padded) * max_pb);

		for (int pc = 0; pc < n; pc += blocking.kc) {
			const int pb = std::min(blocking.kc, n - pc);

			for (int jcg = 0; jcg < k; jcg += blocking.nc_group) {
				const int ng = std::min(blocking.nc_group, k - jcg);

#pragma omp single
				{
					packed_b.resize(static_cast<std::size_t>(pb) * ng);
					pack_b_panel(n, pb, ng, B + pc + static_cast<std::size_t>(n) * jcg, packed_b.data());
				}

#pragma omp for collapse(2) schedule(static)
				for (int jc_offset = 0; jc_offset < ng; jc_offset += blocking.column_tile) {
					for (int ic = 0; ic < m; ic += blocking.mc) {
						const int nb = std::min(blocking.column_tile, ng - jc_offset);
						compute_packed_ic_column_tile(
							m, pb, nb, ic, pc, jcg + jc_offset,
							blocking, A, packed_b.data() + static_cast<std::size_t>(pb) * jc_offset, packed_a.data(), CU, CD);
					}
				}
			}
		}
	}
}

} // namespace vblas_udmatmul_detail

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU and CD are updated in-place with upward/downward rounded fused FMA products.
inline void udmatmul_fused(int m, int n, int k, double* A, double* B, double* CU, double* CD) {
	using namespace vblas_udmatmul_detail;

	if (m <= 0 || n <= 0 || k <= 0) {
		return;
	}

	const Blocking blocking = select_blocking(m, n, k);
	if (blocking.direct) {
		udmatmul_direct(m, n, k, A, B, CU, CD, blocking);
	}
	else if (blocking.mode == PARALLEL_M_ONLY) {
		udmatmul_packed_m_only(m, n, k, A, B, CU, CD, blocking);
	}
	else {
		udmatmul_packed_2d(m, n, k, A, B, CU, CD, blocking);
	}
}

inline void udmatmul(int m, int n, int k, double* A, double* B, double* CU, double* CD) {
	udmatmul_fused(m, n, k, A, B, CU, CD);
}

#endif // VBLAS_UDMATMUL_HPP
