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

#ifndef VBLAS_RDBLAS_OMP_HPP
#define VBLAS_RDBLAS_OMP_HPP

// --- rdblas-layer OpenMP guard (IOP-1 B, same form as vcp/tblas/tblas_common.hpp) ---
// VCP_NOMP              : global kill switch (implies VCP_RDBLAS_NOMP)
// VCP_TBLAS_NOMP        : tblas-layer kill switch (implies VCP_RDBLAS_NOMP: tblas / tlapack
//                         for double call rdblas)
// VCP_RDBLAS_NOMP       : kill switch of vcp/vblas/ (rdblas, rmatmul, udmatmul)
// VCP_RDBLAS_USE_OPENMP : 0/1, the ONLY symbol implementation code tests
// Without -fopenmp (_OPENMP is not defined) or with a kill switch, <omp.h> is not included
// and every routine runs on the calling thread only.
#ifdef VCP_NOMP
#  ifndef VCP_RDBLAS_NOMP
#    define VCP_RDBLAS_NOMP
#  endif
#endif
#ifdef VCP_TBLAS_NOMP
#  ifndef VCP_RDBLAS_NOMP
#    define VCP_RDBLAS_NOMP
#  endif
#endif
#if defined(_OPENMP) && !defined(VCP_RDBLAS_NOMP)
#  define VCP_RDBLAS_USE_OPENMP 1
#  include <omp.h>
#else
#  define VCP_RDBLAS_USE_OPENMP 0
#endif

// Stand-ins for the omp_ functions so that the code compiles without OpenMP.
namespace vblas_omp_detail {

inline int max_threads() {
#if VCP_RDBLAS_USE_OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}

inline int num_procs() {
#if VCP_RDBLAS_USE_OPENMP
	return omp_get_num_procs();
#else
	return 1;
#endif
}

inline int region_threads() {
#if VCP_RDBLAS_USE_OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

inline int region_thread_id() {
#if VCP_RDBLAS_USE_OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

} // namespace vblas_omp_detail

#endif // VBLAS_RDBLAS_OMP_HPP
