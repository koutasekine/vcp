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

#ifndef VBLAS_RDBLAS_HPP
#define VBLAS_RDBLAS_HPP

// rdblas: 丸めモード指定付き double BLAS (LAPACK 構築に必要な全 routine)
//
// - FP 演算を含む routine は BLAS 名に r を付け，引数の末尾に rounding_mode を
//   追加する (1: 上向き丸め, -1: 下向き丸め, それ以外: 最近点丸め)．
//   例: dgemm -> rdgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb,
//                       beta, C, ldc, rounding_mode)
// - FP 演算を含まない routine は BLAS と同名のまま: dcopy, dswap, idamax
//   (idamax の返り値は 0-based index)
// - 行列は column-major．全 routine について，呼び出し thread と OpenMP の
//   各 worker thread の丸めモードを計算前に保存し，指定丸めへ変更して計算，
//   終了後に各 thread とも元の丸めモードへ復元する．
// - O(N^3) routine (rdgemm / rdgemmtr / rdsymm / rdsyrk / rdsyr2k / rdtrmm /
//   rdtrsm) の核は rmatmul (AVX-512 / AVX2+FMA / NEON / no-SIMD を compile 時に
//   自動選択) を使う．
//
// Level 1: rdscal rdaxpy rddot rdasum rdnrm2 rdrot rdrotg rdrotm rdrotmg
//          dcopy dswap idamax
// Level 2: rdgemv rdgbmv rdsymv rdsbmv rdspmv rdtrmv rdtbmv rdtpmv
//          rdtrsv rdtbsv rdtpsv rdger rdsyr rdspr rdsyr2 rdspr2
// Level 3: rdgemm rdgemmtr rdsymm rdsyrk rdsyr2k rdtrmm rdtrsm
//          (rdgemmtr は LAPACK 3.12.1 で追加された GEMMTR に対応)

#include "rdblas_common.hpp"
#include "rdblas_level1.hpp"
#include "rdblas_level2.hpp"
#include "rdblas_level3.hpp"

#endif // VBLAS_RDBLAS_HPP
