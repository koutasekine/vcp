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
// Modified by Kouta Sekine (2026): translated to C++, BLAS calls replaced by
// rdblas, rounding-mode argument added, workspace arguments removed.
// ----------------------------------------------------------------------------

#pragma once

#ifndef VBLAS_RDLAPACK_COMMON_HPP
#define VBLAS_RDLAPACK_COMMON_HPP

// rdlapack (丸めモード指定付き double LAPACK) の共通基盤．
// reference LAPACK 3.12.1 (BSD license, http://www.netlib.org/lapack/) を
// BLAS 呼び出しを rdblas に置き換えながら C++ へ移植したもの．
//
// - rounding_mode の規約は rdblas と同じ: 1 上向き，-1 下向き，それ以外は最近点
// - 各 routine は入口で呼び出し thread の丸めモードを保存して指定丸めへ変更し，
//   終了時に復元する (RoundingGuard)．O(N^3) 部分は rdblas が実行し，
//   rdblas は OpenMP の各 worker thread についても同じ保存・変更・復元を行う．
// - 行列は column-major，index (ipiv 等) は 0-based．
// - LAPACK の INFO は返り値 (int) で返す．INFO < 0 に相当する引数 error は
//   rdblas と同様に std::invalid_argument を投げる (xerbla は呼ばない)．
// - LAPACK の WORK/LWORK 引数は持たず，作業領域は内部で std::vector により
//   最適 size を確保する．

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <vector>

#include <vcp/vblas/rdblas.hpp>

#pragma STDC FENV_ACCESS ON

namespace vblas_rdlapack_detail {

using vblas_rdblas_detail::RoundingGuard;
using vblas_rdblas_detail::fe_rounding;
using vblas_rdblas_detail::option_is;

inline void rdlapack_error(const char* message) {
	vblas_rdblas_detail::rdblas_error(message);
}

// Fortran の SIGN(a, b): b の符号 (IEEE, -0 は負) を |a| に付ける
inline double f_sign(const double a, const double b) {
	return std::signbit(b) ? -std::fabs(a) : std::fabs(a);
}

// LAPACK の ILAENV 相当 (double 実数 routine のみ，reference LAPACK 3.12.1 の既定値)．
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

} // namespace vblas_rdlapack_detail

namespace vcp {

// LAPACK の DLAMCH: machine parameter を返す (定数のみで丸めモードに依存しない)
//   'E' eps (相対 machine epsilon, 2^-53), 'S' safe minimum, 'B' base,
//   'P' eps*base, 'N' 仮数部 bit 数, 'R' rounding (1.0), 'M' emin,
//   'U' underflow threshold, 'L' emax, 'O' overflow threshold
inline double dlamch(const char cmach) {
	using vblas_rdlapack_detail::option_is;
	if (option_is(cmach, 'E')) {
		return DBL_EPSILON * 0.5;
	}
	if (option_is(cmach, 'S')) {
		return DBL_MIN;
	}
	if (option_is(cmach, 'B')) {
		return 2.0;
	}
	if (option_is(cmach, 'P')) {
		return DBL_EPSILON;
	}
	if (option_is(cmach, 'N')) {
		return 53.0;
	}
	if (option_is(cmach, 'R')) {
		return 1.0;
	}
	if (option_is(cmach, 'M')) {
		return -1021.0;
	}
	if (option_is(cmach, 'U')) {
		return DBL_MIN;
	}
	if (option_is(cmach, 'L')) {
		return 1024.0;
	}
	if (option_is(cmach, 'O')) {
		return DBL_MAX;
	}
	vblas_rdlapack_detail::rdlapack_error("dlamch: invalid cmach");
	return 0.0;
}

} // namespace vcp

#endif // VBLAS_RDLAPACK_COMMON_HPP
