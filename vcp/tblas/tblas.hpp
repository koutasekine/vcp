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

#ifndef TBLAS_TBLAS_HPP
#define TBLAS_TBLAS_HPP

// tblas: 要素型 T を template parameter とする BLAS (umbrella header)
//
// - 関数名は double BLAS の d を t に読み替えたもの (dgemm -> tgemm<T> など，
//   idamax -> itamax)．引数は reference BLAS と同順 (rounding_mode はない)
// - 行列は column-major，index は 0-based，scalar は const T& 渡し
// - 丸めモードは扱わない．丸めモード指定が必要な場合は vblas/rdblas を使うこと
// - 要素型 T の要件 (R1-R4) と各 routine の対応は tblas/README_tblas.md 参照
//
// 提供 routine:
//   Level 1: tcopy tswap itamax tscal taxpy tdot tasum tnrm2
//            trot trotg trotm trotmg
//   Level 2: tgemv tgbmv tsymv tsbmv tspmv ttrmv ttbmv ttpmv ttrsv ttbsv ttpsv
//            tger tsyr tspr tsyr2 tspr2
//   Level 3: tgemm tgemmtr tsymm tsyrk tsyr2k ttrmm ttrsm

#include "tblas_common.hpp"
#include "tblas_level1.hpp"
#include "tblas_level2.hpp"
#include "tblas_level3.hpp"

#endif // TBLAS_TBLAS_HPP
