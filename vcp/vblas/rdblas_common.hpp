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

#ifndef VBLAS_RDBLAS_COMMON_HPP
#define VBLAS_RDBLAS_COMMON_HPP

// rdblas (丸めモード指定付き double BLAS) の共通基盤．
// - rounding_mode の規約は rmatmul と同じ: 1 上向き，-1 下向き，それ以外は最近点
// - 丸めモードの変更は「呼び出した thread と OpenMP の各 worker thread」の
//   両方で行い，計算終了後はそれぞれの thread が計算前のモードへ復元する
//   (RoundingGuard による RAII)

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>

#include <omp.h>

#include "rmatmul_common.hpp"

#pragma STDC FENV_ACCESS ON

namespace vblas_rdblas_detail {

inline void rdblas_error(const char* message) {
#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
	throw std::invalid_argument(message);
#else
	(void)message;
	std::abort();
#endif
}

// rounding_mode: 1 上向き，-1 下向き，それ以外は最近点 (rmatmul と同じ規約)
inline int fe_rounding(const int rounding_mode) {
	if (rounding_mode == 1) {
		return FE_UPWARD;
	}
	if (rounding_mode == -1) {
		return FE_DOWNWARD;
	}
	return FE_TONEAREST;
}

// 現在の thread の丸めモードを構築時に保存して fe_mode に変更し，
// 破棄時に保存しておいたモードへ復元する．
// OpenMP parallel region の内側で構築すれば，丸めモードは thread 毎の状態
// なので，全 worker thread について保存・変更・復元が行われる．
class RoundingGuard {
public:
	explicit RoundingGuard(const int fe_mode) : previous_(vblas_rmatmul_common_detail::get_rounding_mode()) {
		vblas_rmatmul_common_detail::set_rounding_mode(fe_mode);
	}
	~RoundingGuard() {
		vblas_rmatmul_common_detail::reset_rounding_mode(previous_);
	}

private:
	int previous_;
	RoundingGuard(const RoundingGuard&);
	RoundingGuard& operator=(const RoundingGuard&);
};

// BLAS の option 文字比較 (lsame 相当，target は大文字で渡す)
inline bool option_is(const char option, const char target) {
	return option == target || option == target + ('a' - 'A');
}

// 増分 incx の vector の先頭 index (BLAS の負増分規約)
inline int vec_start(const int n, const int inc) {
	return inc > 0 ? 0 : (n - 1) * (-inc);
}

// O(N) ~ O(N^2) 演算用の thread 数 (演算量 flops に応じて制限)
inline int threads_for_flops(const double flops) {
	namespace common = vblas_rmatmul_common_detail;
	const common::CpuInfo& info = common::cpu_info();
	const int max_threads = common::limited_threads(std::min(info.logical_threads, info.omp_max_threads));
	if (flops < 4.0e5) {
		return 1;
	}
	if (flops >= 4.0e5 * static_cast<double>(max_threads)) {
		return max_threads;
	}
	return std::max(1, std::min(max_threads, static_cast<int>(flops / 4.0e5) + 1));
}

// C(m x n, ldc 付き) := beta*C．beta==0 は exact な 0 埋め (BLAS 規約: C は読まない)，
// beta==1 は何もしない．乗算は fe_mode の丸めで実行する．
inline void scale_matrix(const int m, const int n, const double beta, double* C, const int ldc, const int fe_mode) {
	if (m <= 0 || n <= 0 || beta == 1.0) {
		return;
	}
	const int threads = threads_for_flops(static_cast<double>(m) * n);
	if (beta == 0.0) {
#pragma omp parallel for schedule(static) num_threads(threads)
		for (int j = 0; j < n; j++) {
			double* c = C + static_cast<std::size_t>(ldc) * j;
			for (int i = 0; i < m; i++) {
				c[i] = 0.0;
			}
		}
		return;
	}
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe_mode);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			double* c = C + static_cast<std::size_t>(ldc) * j;
			for (int i = 0; i < m; i++) {
				c[i] = beta * c[i];
			}
		}
	}
}

// C(m x n, ldc 付き) := alpha*T + beta*C．T は m x n の連続格納．
// 各要素は t1 = alpha*T(i,j) (1 丸め)，C(i,j) = fma(beta, C(i,j), t1) (1 丸め) で
// 全演算を fe_mode の丸めで実行する．beta==0 のとき C は読まない (BLAS 規約)．
inline void apply_alpha_beta(
	const int m,
	const int n,
	const double alpha,
	const double* T,
	const double beta,
	double* C,
	const int ldc,
	const int fe_mode
) {
	if (m <= 0 || n <= 0) {
		return;
	}
	const int threads = threads_for_flops(2.0 * m * n);
#pragma omp parallel num_threads(threads)
	{
		RoundingGuard guard(fe_mode);
#pragma omp for schedule(static)
		for (int j = 0; j < n; j++) {
			const double* t = T + static_cast<std::size_t>(m) * j;
			double* c = C + static_cast<std::size_t>(ldc) * j;
			if (beta == 0.0) {
				if (alpha == 1.0) {
					for (int i = 0; i < m; i++) {
						c[i] = t[i];
					}
				}
				else {
					for (int i = 0; i < m; i++) {
						c[i] = alpha * t[i];
					}
				}
			}
			else if (alpha == 1.0) {
				for (int i = 0; i < m; i++) {
					c[i] = std::fma(beta, c[i], t[i]);
				}
			}
			else {
				for (int i = 0; i < m; i++) {
					c[i] = std::fma(beta, c[i], alpha * t[i]);
				}
			}
		}
	}
}

// dst(m x n 連続) := src(m x n, ld 付き)．FP 演算なしの copy．
inline void pack_no_trans(const int m, const int n, const double* src, const int ld, double* dst) {
	const int threads = threads_for_flops(static_cast<double>(m) * n);
#pragma omp parallel for schedule(static) num_threads(threads)
	for (int j = 0; j < n; j++) {
		const double* s = src + static_cast<std::size_t>(ld) * j;
		double* d = dst + static_cast<std::size_t>(m) * j;
		for (int i = 0; i < m; i++) {
			d[i] = s[i];
		}
	}
}

// dst(m x n 連続) := src(n x m, ld 付き)^T．FP 演算なしの blocked transpose copy．
inline void pack_trans(const int m, const int n, const double* src, const int ld, double* dst) {
	const int TB = 64;
	const int jblocks = (n + TB - 1) / TB;
	const int iblocks = (m + TB - 1) / TB;
	const int threads = threads_for_flops(static_cast<double>(m) * n);
#pragma omp parallel for collapse(2) schedule(static) num_threads(threads)
	for (int jb = 0; jb < jblocks; jb++) {
		for (int ib = 0; ib < iblocks; ib++) {
			const int j0 = jb * TB;
			const int i0 = ib * TB;
			const int j1 = std::min(n, j0 + TB);
			const int i1 = std::min(m, i0 + TB);
			for (int j = j0; j < j1; j++) {
				const double* s = src + j;
				double* d = dst + static_cast<std::size_t>(m) * j;
				for (int i = i0; i < i1; i++) {
					d[i] = s[static_cast<std::size_t>(ld) * i];
				}
			}
		}
	}
}

} // namespace vblas_rdblas_detail

#endif // VBLAS_RDBLAS_COMMON_HPP
