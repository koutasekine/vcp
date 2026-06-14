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

#ifndef VBLAS_DBLAS_HPP
#define VBLAS_DBLAS_HPP

// dblas: rdblas を呼び出して構築する「通常の」double BLAS (Fortran 互換 wrapper)
//
// - 関数名は Fortran BLAS に合わせて末尾 underscore (dgemm_ など)，
//   引数は全て pointer 渡し (Fortran 互換 ABI，INTEGER は 32bit int / LP64)．
// - 各関数は呼び出し時点の丸めモードを std::fegetround() で取得し，
//     FE_DOWNWARD -> -1, FE_UPWARD -> 1, それ以外 -> 0 (最近点)
//   を rdblas の末尾引数 rounding_mode として渡す．
//   つまり「現在の丸めモードを尊重して計算する BLAS」になる
//   (外部 BLAS と異なり，全 OpenMP worker thread でも同じ丸めが保証され，
//   計算後は各 thread の丸めモードが呼び出し前の状態へ復元される)．
// - idamax_ の返り値は Fortran BLAS と同じ 1-based index (n<=0 等では 0)．
// - 引数 error は rdblas に従い std::invalid_argument を投げる (xerbla は呼ばない)．
// - lsame_ / xerbla_ は提供しない (LAPACK や他の BLAS 側の実装を使うこと)．
//
// header-only で extern "C" inline 定義のため，Fortran object と link する
// library を作る場合は，1 つの翻訳単位で VBLAS_DBLAS_EMIT_SYMBOLS を定義して
// 本 header を include し，全関数の symbol を強制的に発行させる:
//
//   // dblas_symbols.cpp
//   #define VBLAS_DBLAS_EMIT_SYMBOLS
//   #include "vblas/dblas.hpp"
//
// このとき MKL 等の他の BLAS と同名 symbol が衝突しないよう，link する
// library 構成に注意すること．

#include <cfenv>

#include "rdblas.hpp"

#pragma STDC FENV_ACCESS ON

namespace vblas_dblas_detail {

// 現在の丸めモードを rdblas の rounding_mode へ変換する
// (FE_DOWNWARD -> -1, FE_UPWARD -> 1, それ以外 -> 0)
inline int current_rounding_mode() {
	const int mode = std::fegetround();
	if (mode == FE_DOWNWARD) {
		return -1;
	}
	if (mode == FE_UPWARD) {
		return 1;
	}
	return 0;
}

} // namespace vblas_dblas_detail

extern "C" {

// ---- Level 1 ----

inline void dcopy_(const int* n, const double* x, const int* incx, double* y, const int* incy) {
	vcp::dcopy(*n, x, *incx, y, *incy);
}

inline void dswap_(const int* n, double* x, const int* incx, double* y, const int* incy) {
	vcp::dswap(*n, x, *incx, y, *incy);
}

// Fortran BLAS と同じ 1-based index を返す (n<=0 や incx<=0 では 0)
inline int idamax_(const int* n, const double* x, const int* incx) {
	return vcp::idamax(*n, x, *incx) + 1;
}

inline void dscal_(const int* n, const double* alpha, double* x, const int* incx) {
	vcp::rdscal(*n, *alpha, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void daxpy_(const int* n, const double* alpha, const double* x, const int* incx, double* y, const int* incy) {
	vcp::rdaxpy(*n, *alpha, x, *incx, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline double ddot_(const int* n, const double* x, const int* incx, const double* y, const int* incy) {
	return vcp::rddot(*n, x, *incx, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline double dasum_(const int* n, const double* x, const int* incx) {
	return vcp::rdasum(*n, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline double dnrm2_(const int* n, const double* x, const int* incx) {
	return vcp::rdnrm2(*n, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void drot_(const int* n, double* x, const int* incx, double* y, const int* incy, const double* c, const double* s) {
	vcp::rdrot(*n, x, *incx, y, *incy, *c, *s, vblas_dblas_detail::current_rounding_mode());
}

inline void drotg_(double* a, double* b, double* c, double* s) {
	vcp::rdrotg(a, b, c, s, vblas_dblas_detail::current_rounding_mode());
}

inline void drotm_(const int* n, double* x, const int* incx, double* y, const int* incy, const double* param) {
	vcp::rdrotm(*n, x, *incx, y, *incy, param, vblas_dblas_detail::current_rounding_mode());
}

inline void drotmg_(double* d1, double* d2, double* x1, const double* y1, double* param) {
	vcp::rdrotmg(d1, d2, x1, *y1, param, vblas_dblas_detail::current_rounding_mode());
}

// ---- Level 2 ----

inline void dgemv_(const char* trans, const int* m, const int* n,
	const double* alpha, const double* a, const int* lda,
	const double* x, const int* incx,
	const double* beta, double* y, const int* incy) {
	vcp::rdgemv(*trans, *m, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline void dgbmv_(const char* trans, const int* m, const int* n, const int* kl, const int* ku,
	const double* alpha, const double* a, const int* lda,
	const double* x, const int* incx,
	const double* beta, double* y, const int* incy) {
	vcp::rdgbmv(*trans, *m, *n, *kl, *ku, *alpha, a, *lda, x, *incx, *beta, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline void dsymv_(const char* uplo, const int* n,
	const double* alpha, const double* a, const int* lda,
	const double* x, const int* incx,
	const double* beta, double* y, const int* incy) {
	vcp::rdsymv(*uplo, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline void dsbmv_(const char* uplo, const int* n, const int* k,
	const double* alpha, const double* a, const int* lda,
	const double* x, const int* incx,
	const double* beta, double* y, const int* incy) {
	vcp::rdsbmv(*uplo, *n, *k, *alpha, a, *lda, x, *incx, *beta, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline void dspmv_(const char* uplo, const int* n,
	const double* alpha, const double* ap,
	const double* x, const int* incx,
	const double* beta, double* y, const int* incy) {
	vcp::rdspmv(*uplo, *n, *alpha, ap, x, *incx, *beta, y, *incy, vblas_dblas_detail::current_rounding_mode());
}

inline void dtrmv_(const char* uplo, const char* trans, const char* diag, const int* n,
	const double* a, const int* lda, double* x, const int* incx) {
	vcp::rdtrmv(*uplo, *trans, *diag, *n, a, *lda, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void dtbmv_(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
	const double* a, const int* lda, double* x, const int* incx) {
	vcp::rdtbmv(*uplo, *trans, *diag, *n, *k, a, *lda, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void dtpmv_(const char* uplo, const char* trans, const char* diag, const int* n,
	const double* ap, double* x, const int* incx) {
	vcp::rdtpmv(*uplo, *trans, *diag, *n, ap, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void dtrsv_(const char* uplo, const char* trans, const char* diag, const int* n,
	const double* a, const int* lda, double* x, const int* incx) {
	vcp::rdtrsv(*uplo, *trans, *diag, *n, a, *lda, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void dtbsv_(const char* uplo, const char* trans, const char* diag, const int* n, const int* k,
	const double* a, const int* lda, double* x, const int* incx) {
	vcp::rdtbsv(*uplo, *trans, *diag, *n, *k, a, *lda, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void dtpsv_(const char* uplo, const char* trans, const char* diag, const int* n,
	const double* ap, double* x, const int* incx) {
	vcp::rdtpsv(*uplo, *trans, *diag, *n, ap, x, *incx, vblas_dblas_detail::current_rounding_mode());
}

inline void dger_(const int* m, const int* n,
	const double* alpha, const double* x, const int* incx,
	const double* y, const int* incy, double* a, const int* lda) {
	vcp::rdger(*m, *n, *alpha, x, *incx, y, *incy, a, *lda, vblas_dblas_detail::current_rounding_mode());
}

inline void dsyr_(const char* uplo, const int* n,
	const double* alpha, const double* x, const int* incx,
	double* a, const int* lda) {
	vcp::rdsyr(*uplo, *n, *alpha, x, *incx, a, *lda, vblas_dblas_detail::current_rounding_mode());
}

inline void dspr_(const char* uplo, const int* n,
	const double* alpha, const double* x, const int* incx, double* ap) {
	vcp::rdspr(*uplo, *n, *alpha, x, *incx, ap, vblas_dblas_detail::current_rounding_mode());
}

inline void dsyr2_(const char* uplo, const int* n,
	const double* alpha, const double* x, const int* incx,
	const double* y, const int* incy, double* a, const int* lda) {
	vcp::rdsyr2(*uplo, *n, *alpha, x, *incx, y, *incy, a, *lda, vblas_dblas_detail::current_rounding_mode());
}

inline void dspr2_(const char* uplo, const int* n,
	const double* alpha, const double* x, const int* incx,
	const double* y, const int* incy, double* ap) {
	vcp::rdspr2(*uplo, *n, *alpha, x, *incx, y, *incy, ap, vblas_dblas_detail::current_rounding_mode());
}

// ---- Level 3 ----

inline void dgemm_(const char* transa, const char* transb,
	const int* m, const int* n, const int* k,
	const double* alpha, const double* a, const int* lda,
	const double* b, const int* ldb,
	const double* beta, double* c, const int* ldc) {
	vcp::rdgemm(*transa, *transb, *m, *n, *k, *alpha, a, *lda, b, *ldb, *beta, c, *ldc, vblas_dblas_detail::current_rounding_mode());
}

// LAPACK 3.12.1 で追加された GEMMTR (C の uplo 三角部分のみ更新する GEMM)
inline void dgemmtr_(const char* uplo, const char* transa, const char* transb,
	const int* n, const int* k,
	const double* alpha, const double* a, const int* lda,
	const double* b, const int* ldb,
	const double* beta, double* c, const int* ldc) {
	vcp::rdgemmtr(*uplo, *transa, *transb, *n, *k, *alpha, a, *lda, b, *ldb, *beta, c, *ldc, vblas_dblas_detail::current_rounding_mode());
}

inline void dsymm_(const char* side, const char* uplo,
	const int* m, const int* n,
	const double* alpha, const double* a, const int* lda,
	const double* b, const int* ldb,
	const double* beta, double* c, const int* ldc) {
	vcp::rdsymm(*side, *uplo, *m, *n, *alpha, a, *lda, b, *ldb, *beta, c, *ldc, vblas_dblas_detail::current_rounding_mode());
}

inline void dsyrk_(const char* uplo, const char* trans,
	const int* n, const int* k,
	const double* alpha, const double* a, const int* lda,
	const double* beta, double* c, const int* ldc) {
	vcp::rdsyrk(*uplo, *trans, *n, *k, *alpha, a, *lda, *beta, c, *ldc, vblas_dblas_detail::current_rounding_mode());
}

inline void dsyr2k_(const char* uplo, const char* trans,
	const int* n, const int* k,
	const double* alpha, const double* a, const int* lda,
	const double* b, const int* ldb,
	const double* beta, double* c, const int* ldc) {
	vcp::rdsyr2k(*uplo, *trans, *n, *k, *alpha, a, *lda, b, *ldb, *beta, c, *ldc, vblas_dblas_detail::current_rounding_mode());
}

inline void dtrmm_(const char* side, const char* uplo, const char* transa, const char* diag,
	const int* m, const int* n,
	const double* alpha, const double* a, const int* lda,
	double* b, const int* ldb) {
	vcp::rdtrmm(*side, *uplo, *transa, *diag, *m, *n, *alpha, a, *lda, b, *ldb, vblas_dblas_detail::current_rounding_mode());
}

inline void dtrsm_(const char* side, const char* uplo, const char* transa, const char* diag,
	const int* m, const int* n,
	const double* alpha, const double* a, const int* lda,
	double* b, const int* ldb) {
	vcp::rdtrsm(*side, *uplo, *transa, *diag, *m, *n, *alpha, a, *lda, b, *ldb, vblas_dblas_detail::current_rounding_mode());
}

} // extern "C"

// Fortran object と link する library を作る場合に，1 つの翻訳単位で
// VBLAS_DBLAS_EMIT_SYMBOLS を定義して全 symbol を強制発行させる
#ifdef VBLAS_DBLAS_EMIT_SYMBOLS
namespace vblas_dblas_detail {

typedef void (*emit_fn)();

// __attribute__((used)) がないと未使用 static 配列として最適化で削除され，
// inline 関数の symbol が発行されない
__attribute__((used)) static const emit_fn emit_table[] = {
	reinterpret_cast<emit_fn>(dcopy_),
	reinterpret_cast<emit_fn>(dswap_),
	reinterpret_cast<emit_fn>(idamax_),
	reinterpret_cast<emit_fn>(dscal_),
	reinterpret_cast<emit_fn>(daxpy_),
	reinterpret_cast<emit_fn>(ddot_),
	reinterpret_cast<emit_fn>(dasum_),
	reinterpret_cast<emit_fn>(dnrm2_),
	reinterpret_cast<emit_fn>(drot_),
	reinterpret_cast<emit_fn>(drotg_),
	reinterpret_cast<emit_fn>(drotm_),
	reinterpret_cast<emit_fn>(drotmg_),
	reinterpret_cast<emit_fn>(dgemv_),
	reinterpret_cast<emit_fn>(dgbmv_),
	reinterpret_cast<emit_fn>(dsymv_),
	reinterpret_cast<emit_fn>(dsbmv_),
	reinterpret_cast<emit_fn>(dspmv_),
	reinterpret_cast<emit_fn>(dtrmv_),
	reinterpret_cast<emit_fn>(dtbmv_),
	reinterpret_cast<emit_fn>(dtpmv_),
	reinterpret_cast<emit_fn>(dtrsv_),
	reinterpret_cast<emit_fn>(dtbsv_),
	reinterpret_cast<emit_fn>(dtpsv_),
	reinterpret_cast<emit_fn>(dger_),
	reinterpret_cast<emit_fn>(dsyr_),
	reinterpret_cast<emit_fn>(dspr_),
	reinterpret_cast<emit_fn>(dsyr2_),
	reinterpret_cast<emit_fn>(dspr2_),
	reinterpret_cast<emit_fn>(dgemm_),
	reinterpret_cast<emit_fn>(dgemmtr_),
	reinterpret_cast<emit_fn>(dsymm_),
	reinterpret_cast<emit_fn>(dsyrk_),
	reinterpret_cast<emit_fn>(dsyr2k_),
	reinterpret_cast<emit_fn>(dtrmm_),
	reinterpret_cast<emit_fn>(dtrsm_)
};

} // namespace vblas_dblas_detail
#endif // VBLAS_DBLAS_EMIT_SYMBOLS

#endif // VBLAS_DBLAS_HPP
