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

#ifndef VBLAS_RDLAPACK_HPP
#define VBLAS_RDLAPACK_HPP

// rdlapack: 丸めモード指定付き double LAPACK (reference LAPACK 3.12.1 の移植)
//
// - FP 演算を含む routine は LAPACK 名の先頭に r を付け，引数の末尾に
//   rounding_mode を追加する (1: 上向き, -1: 下向き, それ以外: 最近点)．
//   例: dgesv -> rdgesv(n, nrhs, A, lda, ipiv, B, ldb, rounding_mode)
// - 全ての浮動小数点演算 (四則演算と平方根, BLAS 部分も scalar 部分も) が
//   指定した同一の丸めモードで実行され，終了時に呼び出し前の丸めモードへ
//   復元される (BLAS 部分は rdblas が OpenMP の全 worker thread についても保証)．
// - これは「全演算を同一丸めモードで計算する」だけであり，
//   精度保証付き数値計算 (区間演算による厳密な誤差評価) ではないことに注意．
// - 行列は column-major，index は 0-based (例外: rdsytrf 系の ipiv は
//   LAPACK と同じ 1-based 符号付き，rdgebal/rdgehrd 系の ilo/ihi は 1-based)．
// - LAPACK の INFO は返り値 (int)．引数 error は std::invalid_argument を投げる．
// - WORK/LWORK 引数は持たず，作業領域は内部で確保する．
//
// 提供 routine (駆動 routine):
//   rdgesv  rdposv  rdsysv  rdgbsv  rdpbsv      (連立一次方程式)
//   rdgels                                      (最小二乗)
//   rdgetri rdpotri rdtrtri                     (逆行列)
//   rdsyev  rdsygv                              (対称固有値)
//   rdgeev                                      (非対称固有値)
//   rdgesvd                                     (特異値分解)
// 各分解・補助 routine は各 header を参照．

#include "rdlapack_common.hpp"
#include "rdlapack_aux.hpp"
#include "rdlapack_lu.hpp"
#include "rdlapack_chol.hpp"
#include "rdlapack_qr.hpp"
#include "rdlapack_eig.hpp"
#include "rdlapack_svd.hpp"
#include "rdlapack_geev.hpp"
#include "rdlapack_sy.hpp"
#include "rdlapack_band.hpp"

#endif // VBLAS_RDLAPACK_HPP
