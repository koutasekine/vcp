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

#ifndef VBLAS_RMATMUL_COMMON_HPP
#define VBLAS_RMATMUL_COMMON_HPP

// rmatmul_avx512 / rmatmul_avx2 / rmatmul_neon / rmatmul_nosimd の
// instruction set に依存しない共通部分 (CPU/cache 検出，blocking 決定，
// micro-panel packing)．SIMD kernel 本体は各 backend header に置く．

#include <algorithm>
#include <cfenv>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <string>

#include <omp.h>

#if defined(__APPLE__)
#include <sys/sysctl.h>
// <sys/sysctl.h> が取り込む <sys/param.h> の関数形式マクロが
// ユーザコード (kv::hwround::roundup() など) と衝突するため解除する
#if defined(roundup)
#undef roundup
#endif
#if defined(howmany)
#undef howmany
#endif
#if defined(powerof2)
#undef powerof2
#endif
#endif

namespace vblas_rmatmul_common_detail {

constexpr int MATRIX_BYTES = 8;

// cache 検出に失敗した場合に用いる固定 block size
constexpr int FALLBACK_KC = 160;
constexpr int FALLBACK_MC = 128;
constexpr int FALLBACK_COLUMN_TILE = 128;
constexpr int FALLBACK_NC_GROUP = 1024;

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
	// L2 を共有する CPU 数．x86 では同一 core 上の hyper-thread 数と一致し，
	// Apple Silicon では cluster 内で L2 を共有する core 数になる．
	int l2_shared_cpus;
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

inline long long sysctl_value(const char* name) {
#if defined(__APPLE__)
	long long value = 0;
	std::size_t size = sizeof(value);
	if (sysctlbyname(name, &value, &size, NULL, 0) == 0 && value > 0) {
		return value;
	}
#else
	(void)name;
#endif
	return 0;
}

inline CpuInfo detect_cpu_info() {
	CpuInfo info = {0, 0, 0, std::max(1, omp_get_num_procs()), std::max(1, omp_get_max_threads()), 1};

	// Linux: sysfs から cache 構成を読む
	bool l2_shared_found = false;
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
			std::string shared_text;
			if (read_first_line(base + "shared_cpu_list", shared_text)) {
				info.l2_shared_cpus = parse_cpu_list_count(shared_text);
				l2_shared_found = true;
			}
		}
		else if (level == 3) {
			info.l3_bytes = std::max(info.l3_bytes, bytes);
		}
	}

	// L2 の shared_cpu_list が読めない場合は hyper-thread sibling 数で代用する
	// (per-core L2 の x86 では両者は一致する)
	if (!l2_shared_found) {
		std::string siblings;
		if (read_first_line("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list", siblings)) {
			info.l2_shared_cpus = parse_cpu_list_count(siblings);
		}
	}

	// macOS (Apple Silicon): sysctl から cache 構成を読む
	if (info.l1d_bytes <= 0) {
		long long l1 = sysctl_value("hw.perflevel0.l1dcachesize");
		if (l1 <= 0) {
			l1 = sysctl_value("hw.l1dcachesize");
		}
		info.l1d_bytes = l1;
	}
	if (info.l2_bytes <= 0) {
		long long l2 = sysctl_value("hw.perflevel0.l2cachesize");
		if (l2 <= 0) {
			l2 = sysctl_value("hw.l2cachesize");
		}
		info.l2_bytes = l2;
		// Apple Silicon の L2 は cluster 内の複数 core で共有される
		const long long cpus_per_l2 = sysctl_value("hw.perflevel0.cpusperl2");
		if (cpus_per_l2 > 0) {
			info.l2_shared_cpus = static_cast<int>(cpus_per_l2);
		}
	}
	if (info.l3_bytes <= 0) {
		info.l3_bytes = sysctl_value("hw.l3cachesize");
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

template<int MR, int NR>
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

template<int MR, int NR>
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

template<int MR>
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

template<int MR, int NR>
inline Blocking finish_blocking(Blocking b, const int m, const int k) {
	b = sanitize_blocking<MR, NR>(b);
	b = adjust_mc_for_m_only<MR>(b, m);
	b = cap_threads_for_work<MR, NR>(b, m, k);
	return sanitize_blocking<MR, NR>(b);
}

// B micro-panel (kc*NR) が L1d の約 1/2 に収まるように kc を決める．
// kc が大きいほど C への蓄積パス回数 (n/kc) が減り C のメモリ往復が減る．
// L1d が検出できない場合は固定値 FALLBACK_KC に退避する．
inline int cache_based_kc(const CpuInfo& info, const int nr) {
	if (info.l1d_bytes <= 0) {
		return clamp_rounded_block(FALLBACK_KC, 64, 512, nr);
	}
	const long long budget = info.l1d_bytes / 2LL;
	const long long per_k = static_cast<long long>(nr) * MATRIX_BYTES;
	return clamp_rounded_block(static_cast<int>(budget / per_k), 64, 512, nr);
}

// packed A block (mc*kc) が L2 の約 1/2 に収まるように mc を決める
// (残りは B column tile と C の通過分に充てる)．
// L2 が検出できない場合は固定値 FALLBACK_MC に退避する．
inline int cache_based_mc(const CpuInfo& info, const int kc, const int mr) {
	if (info.l2_bytes <= 0) {
		return clamp_rounded_block(FALLBACK_MC, mr, 1024, mr);
	}
	// L2 を複数 CPU で共有する構成 (hyper-thread や Apple Silicon の
	// cluster 共有 L2) では，共有 CPU 数で budget を分割する
	const long long budget_doubles = info.l2_bytes / 2LL / MATRIX_BYTES / std::max(1, info.l2_shared_cpus);
	return clamp_rounded_block(static_cast<int>(budget_doubles / kc), mr, 1024, mr);
}

// B column tile (kc*column_tile) が L2 の約 1/4 に収まるように決める．
// 縦長行列で B の通過回数が増えすぎないよう下限は NR*8 とする．
// L2 が検出できない場合は固定値 FALLBACK_COLUMN_TILE に退避する．
inline int cache_based_column_tile(const CpuInfo& info, const int kc, const int nr) {
	if (info.l2_bytes <= 0) {
		return clamp_rounded_block(FALLBACK_COLUMN_TILE, nr * 8, 512, nr);
	}
	const long long budget_doubles = info.l2_bytes / 4LL / MATRIX_BYTES;
	return clamp_rounded_block(static_cast<int>(budget_doubles / kc), nr * 8, 512, nr);
}

// packed B group (kc*nc_group, 全 thread で共有) が最外殻 cache の約 1/2 に
// 収まるように決める．全 thread の packed A block (threads*mc*kc) が同じ
// cache を占有するため，その分を budget から除く．
// L3 が無い構成 (Apple Silicon など) は L2 合計で代用し，
// どちらも検出できない場合は固定値 FALLBACK_NC_GROUP に退避する．
inline int cache_based_nc_group(const CpuInfo& info, const int kc, const int mc, const int column_tile, const int threads, const int nr) {
	long long outer = info.l3_bytes;
	if (outer <= 0) {
		outer = info.l2_bytes;
	}
	if (outer <= 0) {
		return clamp_rounded_block(std::max(column_tile, FALLBACK_NC_GROUP), 128, 8192, nr);
	}
	const long long budget_doubles = outer / 2LL / MATRIX_BYTES;
	const long long active_a = static_cast<long long>(std::max(1, threads)) * mc * kc;
	const long long raw = budget_doubles > active_a ? (budget_doubles - active_a) / kc : FALLBACK_NC_GROUP;
	const long long lifted = std::max(static_cast<long long>(column_tile), raw);
	return clamp_rounded_block(static_cast<int>(std::min(lifted, 8192LL)), 128, 8192, nr);
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

template<int MR, int NR>
inline Blocking select_blocking(const int m, const int n, const int k) {
	const CpuInfo& info = cpu_info();
	const int threads = default_thread_count(info, m, n, k);

	if (use_direct_path(m, n, k)) {
		const double work = static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
		const int direct_threads = work >= 1.0e7 ? threads : std::min(threads, 8);
		return finish_blocking<MR, NR>(Blocking{MR, std::max(NR, n), std::max(NR, k), NR, direct_threads, PARALLEL_2D, true}, m, k);
	}

	// cache size と thread 数のみから block size を導出する (環境特化の分岐は置かない)
	const int kc = cache_based_kc(info, NR);
	const int mc = cache_based_mc(info, kc, MR);
	const int column_tile = cache_based_column_tile(info, kc, NR);
	const int nc_group = cache_based_nc_group(info, kc, mc, column_tile, threads, NR);

	// m 方向の block 数が thread 数以上なら m 方向のみの並列で負荷分散できる
	const ParallelMode mode = ceil_div(m, std::max(MR, mc)) >= std::max(1, threads) ? PARALLEL_M_ONLY : PARALLEL_2D;
	return finish_blocking<MR, NR>(Blocking{mc, kc, nc_group, column_tile, threads, mode, false}, m, k);
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
template<int MR>
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
template<int NR>
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

} // namespace vblas_rmatmul_common_detail

#endif // VBLAS_RMATMUL_COMMON_HPP
