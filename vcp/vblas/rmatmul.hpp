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

#include <cstdlib>
#include <stdexcept>

#if defined(__AVX512F__)
#define VBLAS_RMATMUL_HAS_AVX512 1
#include "rmatmul_avx512.hpp"
#else
#define VBLAS_RMATMUL_HAS_AVX512 0
#endif

#if defined(__AVX2__) && defined(__FMA__)
#define VBLAS_RMATMUL_HAS_AVX2 1
#include "rmatmul_avx2.hpp"
#else
#define VBLAS_RMATMUL_HAS_AVX2 0
#endif

#if (defined(__aarch64__) || defined(__arm64__)) && (defined(__ARM_NEON) || defined(__ARM_NEON__))
#define VBLAS_RMATMUL_HAS_NEON 1
#include "rmatmul_neon.hpp"
#else
#define VBLAS_RMATMUL_HAS_NEON 0
#endif

#if !VBLAS_RMATMUL_HAS_AVX512 && !VBLAS_RMATMUL_HAS_AVX2 && !VBLAS_RMATMUL_HAS_NEON
#define VBLAS_RMATMUL_HAS_NOSIMD 1
#include "rmatmul_nosimd.hpp"
#else
#define VBLAS_RMATMUL_HAS_NOSIMD 0
#endif

namespace vblas_rmatmul_detail {

enum Backend {
	BACKEND_NONE,
	BACKEND_AVX512,
	BACKEND_AVX2,
	BACKEND_NEON,
	BACKEND_NOSIMD
};

inline Backend selected_backend() {
#if VBLAS_RMATMUL_HAS_AVX512
	return BACKEND_AVX512;
#elif VBLAS_RMATMUL_HAS_AVX2
	return BACKEND_AVX2;
#elif VBLAS_RMATMUL_HAS_NEON
	return BACKEND_NEON;
#elif VBLAS_RMATMUL_HAS_NOSIMD
	return BACKEND_NOSIMD;
#else
	return BACKEND_NONE;
#endif
}

inline void unavailable_backend() {
#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
	throw std::runtime_error("vblas rmatmul: no available implementation for this target");
#else
	std::abort();
#endif
}

} // namespace vblas_rmatmul_detail

// A: m*n matrix, B: n*k matrix, all matrices are column-major.
// CU is overwritten with A*B using rounding_mode: 1 upward, 0 nearest, -1 downward.
inline void rmatmul(int m, int n, int k, const double* A, const double* B, double* CU, const int rounding_mode) {
	switch (vblas_rmatmul_detail::selected_backend()) {
#if VBLAS_RMATMUL_HAS_AVX512
	case vblas_rmatmul_detail::BACKEND_AVX512:
		rmatmul_avx512(m, n, k, A, B, CU, rounding_mode);
		return;
#endif
#if VBLAS_RMATMUL_HAS_AVX2
	case vblas_rmatmul_detail::BACKEND_AVX2:
		rmatmul_avx2(m, n, k, A, B, CU, rounding_mode);
		return;
#endif
#if VBLAS_RMATMUL_HAS_NEON
	case vblas_rmatmul_detail::BACKEND_NEON:
		rmatmul_neon(m, n, k, A, B, CU, rounding_mode);
		return;
#endif
#if VBLAS_RMATMUL_HAS_NOSIMD
	case vblas_rmatmul_detail::BACKEND_NOSIMD:
		rmatmul_nosimd(m, n, k, A, B, CU, rounding_mode);
		return;
#endif
	default:
		vblas_rmatmul_detail::unavailable_backend();
	}
}

#endif // VBLAS_RMATMUL_HPP
