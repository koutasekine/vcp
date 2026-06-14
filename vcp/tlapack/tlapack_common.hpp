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
// This file contains code derived from reference LAPACK 3.12.1
// (https://www.netlib.org/lapack/),
// Copyright (c) 1992-2023 The University of Tennessee and The University
//                         of Tennessee Research Foundation,
// Copyright (c) 2000-2023 The University of California Berkeley,
// Copyright (c) 2006-2023 The University of Colorado Denver.
// All rights reserved.
// Distributed under the modified BSD license; see vlapack/LICENSE_LAPACK.txt
// for the full license text, the list of conditions and the disclaimer.
// Modified by Kouta Sekine (2026): translated to C++ templates over a generic
// element type T, BLAS calls replaced by tblas (no rounding-mode control),
// workspace arguments removed.
// ----------------------------------------------------------------------------

#pragma once

#ifndef TLAPACK_TLAPACK_COMMON_HPP
#define TLAPACK_TLAPACK_COMMON_HPP

// tlapack (template LAPACK) の共通基盤．vlapack/rdlapack (reference LAPACK 3.12.1
// 由来) を，BLAS 呼び出しを tblas に置き換えながら template 化したもの．
//
// - 丸めモードは扱わない (丸めモード指定が必要な場合は vlapack/rdlapack を使う)
// - 要素型 T の要件は tlapack/README_tlapack.md (tblas の R1-R4 に加えて
//   R5: std::numeric_limits<T> の特殊化，R6: NaN との == が false)
// - 行列は column-major，index (ipiv 等) の 0-based / 1-based 規約は rdlapack と同一
// - LAPACK の INFO は返り値 (int)．INFO < 0 に相当する引数 error は
//   std::invalid_argument を投げる (例外無効時は abort)
// - LAPACK の WORK/LWORK 引数は持たず，作業領域は内部で std::vector<T> により確保する

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include <vcp/tblas/tblas.hpp>

namespace tlapack_detail {

using tblas_detail::option_is;

inline void tlapack_error(const char* message) {
	tblas_detail::tblas_error(message);
}

// abs / sqrt を ADL で解決する (組み込み型は std::，kv 型は namespace kv の関数)
template <typename T>
inline T tabs(const T& x) {
	using std::abs;
	return abs(x);
}

template <typename T>
inline T tsqrt(const T& x) {
	using std::sqrt;
	return sqrt(x);
}

// NaN 判定 (要件 R6)．組み込み型は std::isnan，それ以外は既定で !(x == x)
// (NaN との == が false を返す型で機能する．kv::dd / kv::mpfr は満たす)．
// 個別の型でより適切な判定が必要なら，この overload より優先される
// 非 template の overload を tlapack_detail に追加すればよい．
inline bool tisnan(const float x) { return std::isnan(x); }
inline bool tisnan(const double x) { return std::isnan(x); }
inline bool tisnan(const long double x) { return std::isnan(x); }

template <typename T>
inline bool tisnan(const T& x) {
	return !(x == x);
}

// Fortran の SIGN(a, b): b の符号を |a| に付ける．
// rdlapack の f_sign (std::signbit) と異なり b = -0.0 を正と扱うが (要件 R7)，
// 使用箇所 (Householder の符号選択等) では結果の正しさに影響しない．
template <typename T>
inline T f_sign(const T& a, const T& b) {
	return (b < T(0)) ? -tabs(a) : tabs(a);
}

// 2^e を 2 進冪乗で構成する (O(log|e|) 回の乗算，T の演算のみ)．
// kv::mpfr のように指数 range が極端に広い型でも使えるよう long long を取る．
template <typename T>
inline T pow2i(const long long e) {
	T r = T(1);
	T b = (e >= 0) ? T(2) : (T(1) / T(2));
	unsigned long long k = (e >= 0) ? static_cast<unsigned long long>(e)
	                                : static_cast<unsigned long long>(-(e + 1)) + 1ULL;
	while (k != 0ULL) {
		if (k & 1ULL) {
			r = r * b;
		}
		if (k > 1ULL) {
			b = b * b;
		}
		k >>= 1;
	}
	return r;
}

// LAPACK の ILAENV 相当 (実数 routine のみ，reference LAPACK 3.12.1 の既定値)．
// name は精度文字を除いた大文字名 ("GETRF" など)．
// ispec = 1: block size NB, 2: 最小 block size NBMIN, 3: crossover point NX
inline int ilaenv(const int ispec, const char* name, const int n1, const int n2, const int n3, const int n4) {
	(void)n1;
	(void)n3;
	if (ispec == 1) {
		if (std::strcmp(name, "GETRF") == 0 || std::strcmp(name, "GETRI") == 0) {
			return 64;
		}
		if (std::strcmp(name, "GEQRF") == 0 || std::strcmp(name, "GERQF") == 0 ||
			std::strcmp(name, "GELQF") == 0 || std::strcmp(name, "GEQLF") == 0) {
			return 32;
		}
		if (std::strcmp(name, "GEHRD") == 0 || std::strcmp(name, "GEBRD") == 0) {
			return 32;
		}
		if (std::strcmp(name, "POTRF") == 0 || std::strcmp(name, "SYTRF") == 0 ||
			std::strcmp(name, "SYGST") == 0) {
			return 64;
		}
		if (std::strcmp(name, "SYTRD") == 0) {
			return 32;
		}
		if (std::strncmp(name, "ORG", 3) == 0 || std::strncmp(name, "ORM", 3) == 0) {
			return 32;
		}
		if (std::strcmp(name, "GBTRF") == 0) {
			return n4 <= 64 ? 1 : 32;
		}
		if (std::strcmp(name, "PBTRF") == 0) {
			return n2 <= 64 ? 1 : 32;
		}
		if (std::strcmp(name, "TRTRI") == 0 || std::strcmp(name, "LAUUM") == 0 ||
			std::strcmp(name, "TREVC") == 0) {
			return 64;
		}
		return 1;
	}
	if (ispec == 2) {
		if (std::strcmp(name, "SYTRF") == 0) {
			return 8;
		}
		return 2;
	}
	if (ispec == 3) {
		if (std::strcmp(name, "GEQRF") == 0 || std::strcmp(name, "GERQF") == 0 ||
			std::strcmp(name, "GELQF") == 0 || std::strcmp(name, "GEQLF") == 0 ||
			std::strcmp(name, "GEHRD") == 0 || std::strcmp(name, "GEBRD") == 0) {
			return 128;
		}
		if (std::strcmp(name, "SYTRD") == 0) {
			return 32;
		}
		if (std::strncmp(name, "ORG", 3) == 0) {
			return 128;
		}
		return 0;
	}
	return -1;
}

} // namespace tlapack_detail

// LAPACK の DLAMCH の template 版: machine parameter を std::numeric_limits<T> から
// 構成して返す (要件 R5)．
//   'E' eps (丸め単位 = epsilon()/2)，'S' safe minimum (1/sfmin が overflow しない
//   最小の正数)，'B' base，'P' eps*base = epsilon()，'N' 仮数部 bit 数，
//   'R' rounding (1)，'M' emin，'U' underflow threshold，'L' emax，
//   'O' overflow threshold
template <typename T>
inline T tlamch(const char cmach) {
	typedef std::numeric_limits<T> lim;
	using tlapack_detail::option_is;
	if (option_is(cmach, 'E')) {
		return lim::epsilon() / T(2);
	}
	if (option_is(cmach, 'S')) {
		// reference DLAMCH と同じ: 1/max() >= min() なら (1 + eps)/max() を返す
		T sfmin = lim::min();
		const T small_ = T(1) / lim::max();
		if (small_ >= sfmin) {
			sfmin = small_ * (T(1) + lim::epsilon());
		}
		return sfmin;
	}
	if (option_is(cmach, 'B')) {
		return T(lim::radix);
	}
	if (option_is(cmach, 'P')) {
		return lim::epsilon();
	}
	if (option_is(cmach, 'N')) {
		return T(lim::digits);
	}
	if (option_is(cmach, 'R')) {
		return T(1);
	}
	if (option_is(cmach, 'M')) {
		return T(lim::min_exponent);
	}
	if (option_is(cmach, 'U')) {
		return lim::min();
	}
	if (option_is(cmach, 'L')) {
		return T(lim::max_exponent);
	}
	if (option_is(cmach, 'O')) {
		return lim::max();
	}
	tlapack_detail::tlapack_error("tlamch: invalid cmach");
	return T(0);
}

#endif // TLAPACK_TLAPACK_COMMON_HPP
