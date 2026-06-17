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

#ifndef TLAPACK_TLAPACK_HPP
#define TLAPACK_TLAPACK_HPP

// tlapack: 要素型 T を template parameter とする LAPACK
// (reference LAPACK 3.12.1 由来の vlapack/rdlapack の template 化)
//
// - 関数名は double LAPACK の d を t に読み替えたもの (dgesv -> tgesv<T> など)．
//   引数は rdlapack から rounding_mode を除いたもの (= reference LAPACK と同順)．
// - 丸めモードは扱わない．丸めモード指定が必要な場合は vlapack/rdlapack を使う．
// - 要素型 T の要件は tlapack/README_tlapack.md:
//   tblas の R1-R4 (四則・比較・ADL の abs/sqrt) に加えて
//   R5: std::numeric_limits<T> の特殊化 (epsilon/min/max/radix/digits)，
//   R6: NaN との == が false (既定の tisnan = !(x == x) が機能すること)．
//   double / float / long double / kv::dd / kv::mpfr<N> は追加作業なしで使える．
//   区間型は対象外 (分岐が多く数学的に不適切)．
// - 行列は column-major，index は 0-based (例外: tsytrf 系の ipiv は
//   LAPACK と同じ 1-based 符号付き，tgebal/tgehrd 系の ilo/ihi は 1-based)．
// - LAPACK の INFO は返り値 (int)．引数 error は std::invalid_argument を投げる．
// - WORK/LWORK 引数は持たず，作業領域は内部で確保する．
//
// 提供 routine (駆動 routine):
//   tgesv  tposv  tsysv  tgbsv  tpbsv      (連立一次方程式)
//   tgels                                  (最小二乗)
//   tgetri tpotri ttrtri                   (逆行列)
//   tsyev  tsygv                           (対称固有値)
//   tgeev                                  (非対称固有値)
//   tgesvd                                 (特異値分解)
// 各分解・補助 routine は各 header を参照．

#include "tlapack_common.hpp"
#include "tlapack_aux.hpp"
#include "tlapack_lu.hpp"
#include "tlapack_chol.hpp"
#include "tlapack_qr.hpp"
#include "tlapack_eig.hpp"
#include "tlapack_svd.hpp"
#include "tlapack_geev.hpp"
#include "tlapack_sy.hpp"
#include "tlapack_band.hpp"

#endif // TLAPACK_TLAPACK_HPP
