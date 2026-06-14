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

#ifndef VLAPACK_DLAPACK_HPP
#define VLAPACK_DLAPACK_HPP

// dlapack: rdlapack を呼び出して構築する「通常の」double LAPACK (Fortran 互換 wrapper)
//
// - 関数名は Fortran LAPACK に合わせて末尾 underscore (dgesv_ など)，
//   引数は全て pointer 渡し (Fortran 互換 ABI，INTEGER は 32bit int / LP64)．
// - 各関数は呼び出し時点の丸めモードを取得し，
//     FE_DOWNWARD -> -1, FE_UPWARD -> 1, それ以外 -> 0 (最近点)
//   を rdlapack の末尾引数 rounding_mode として渡す (dblas と同じ規約)．
//   つまり「現在の丸めモードを尊重して計算する LAPACK」になる．
// - 内部の BLAS は rdlapack が rdblas を C++ 直接呼び出しで使う (固定)．
//   dgemm_ などの BLAS シンボルを参照も提供もしないため，リンクされた既存
//   BLAS (MKL 等) とは干渉しないが，その速度の恩恵も受けない
//   (詳細は vlapack/README_dlapack.md)．
//
// LAPACK との互換性の差異 (詳細は README_dlapack.md):
// - ipiv は LAPACK と同じ 1-based に変換して入出力する
//   (dsytrf 系は rdlapack 自体が LAPACK 互換のためそのまま)．
// - WORK は使用しない (作業領域は内部確保)．LWORK = -1 の workspace query
//   には推奨 size を WORK(1) に返す．それ以外の LWORK の値は検査せず無視する．
// - 引数 error (LAPACK の INFO < 0) は rdlapack が投げる例外を catch して
//   INFO = -1 に変換する (どの引数かは特定しない)．
// - dtrevc_ は HOWMNY = 'A'/'B' のみ (SELECT 指定 'S' は INFO = -2)．
// - 提供しない routine (複素数/単精度，*con/*rfs 等) は本 header にはない．
//   liblapack と同時リンクした場合，それらはそちらに解決される (混在に注意)．
//
// header-only で extern "C" inline 定義のため，Fortran object と link する
// library を作る場合は，1 つの翻訳単位で VLAPACK_DLAPACK_EMIT_SYMBOLS を
// 定義して本 header を include し，全関数の symbol を強制的に発行させる．
// このとき liblapack や MKL と同名 symbol の衝突に注意すること．

#include <cfenv>
#include <stdexcept>

#include "rdlapack.hpp"

#pragma STDC FENV_ACCESS ON

namespace vlapack_dlapack_detail {

// 現在の丸めモードを rdlapack の rounding_mode へ変換する
// (FE_DOWNWARD -> -1, FE_UPWARD -> 1, それ以外 -> 0)
inline int current_rounding_mode() {
	const int mode = vblas_rmatmul_common_detail::get_rounding_mode();
	if (mode == FE_DOWNWARD) {
		return -1;
	}
	if (mode == FE_UPWARD) {
		return 1;
	}
	return 0;
}

// ipiv: rdlapack (0-based) -> LAPACK (1-based)
inline void ipiv_to_fortran(const int n, int* ipiv) {
	for (int i = 0; i < n; i++) {
		ipiv[i] += 1;
	}
}

// ipiv: LAPACK (1-based) -> rdlapack (0-based) の一時 copy
inline std::vector<int> ipiv_to_c(const int n, const int* ipiv) {
	std::vector<int> ip0(static_cast<std::size_t>(n > 0 ? n : 0));
	for (int i = 0; i < n; i++) {
		ip0[i] = ipiv[i] - 1;
	}
	return ip0;
}

} // namespace vlapack_dlapack_detail

// 例外 (引数 error) を LAPACK の INFO = -1 に変換する共通 macro
#define VLAPACK_DLAPACK_TRY(stmt) \
	try { stmt } catch (const std::invalid_argument&) { *info = -1; return; }

extern "C" {

// ================= machine parameter / 補助 (FP 演算の有無は rdlapack に従う) =================

inline double dlamch_(const char* cmach) {
	try {
		return vcp::dlamch(*cmach);
	}
	catch (const std::invalid_argument&) {
		return 0.0;
	}
}

inline double dlapy2_(const double* x, const double* y) {
	return vcp::rdlapy2(*x, *y, vlapack_dlapack_detail::current_rounding_mode());
}

inline double dlapy3_(const double* x, const double* y, const double* z) {
	return vcp::rdlapy3(*x, *y, *z, vlapack_dlapack_detail::current_rounding_mode());
}

inline void dlassq_(const int* n, const double* x, const int* incx, double* scale, double* sumsq) {
	vcp::rdlassq(*n, x, *incx, *scale, *sumsq, vlapack_dlapack_detail::current_rounding_mode());
}

inline void dlartg_(const double* f, const double* g, double* cs, double* sn, double* r) {
	vcp::rdlartg(*f, *g, *cs, *sn, *r, vlapack_dlapack_detail::current_rounding_mode());
}

inline void dlarfg_(const int* n, double* alpha, double* x, const int* incx, double* tau) {
	vcp::rdlarfg(*n, *alpha, x, *incx, *tau, vlapack_dlapack_detail::current_rounding_mode());
}

inline void dlascl_(const char* type, const int* kl, const int* ku, const double* cfrom, const double* cto,
	const int* m, const int* n, double* a, const int* lda, int* info) {
	*info = 0;
	VLAPACK_DLAPACK_TRY(
		vcp::rdlascl(*type, *kl, *ku, *cfrom, *cto, *m, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dlaset_(const char* uplo, const int* m, const int* n, const double* alpha, const double* beta,
	double* a, const int* lda) {
	vcp::dlaset(*uplo, *m, *n, *alpha, *beta, a, *lda);
}

inline void dlacpy_(const char* uplo, const int* m, const int* n, const double* a, const int* lda,
	double* b, const int* ldb) {
	vcp::dlacpy(*uplo, *m, *n, a, *lda, b, *ldb);
}

inline void dlasrt_(const char* id, const int* n, double* d, int* info) {
	*info = 0;
	VLAPACK_DLAPACK_TRY(vcp::dlasrt(*id, *n, d);)
}

// k1/k2/ipiv は LAPACK と同じ 1-based (incx = 1 / -1 のみ対応)
inline void dlaswp_(const int* n, double* a, const int* lda, const int* k1, const int* k2,
	const int* ipiv, const int* incx) {
	const std::vector<int> ip0 = vlapack_dlapack_detail::ipiv_to_c(*k2, ipiv);
	vcp::dlaswp(*n, a, *lda, *k1 - 1, *k2, ip0.data(), *incx);
}

inline double dlange_(const char* norm, const int* m, const int* n, const double* a, const int* lda,
	double* work) {
	(void)work;
	try {
		return vcp::rdlange(*norm, *m, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	}
	catch (const std::invalid_argument&) {
		return 0.0;
	}
}

inline double dlangb_(const char* norm, const int* n, const int* kl, const int* ku,
	const double* ab, const int* ldab, double* work) {
	(void)work;
	try {
		return vcp::rdlangb(*norm, *n, *kl, *ku, ab, *ldab, vlapack_dlapack_detail::current_rounding_mode());
	}
	catch (const std::invalid_argument&) {
		return 0.0;
	}
}

inline double dlansy_(const char* norm, const char* uplo, const int* n, const double* a, const int* lda,
	double* work) {
	(void)work;
	try {
		return vcp::rdlansy(*norm, *uplo, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	}
	catch (const std::invalid_argument&) {
		return 0.0;
	}
}

inline double dlansb_(const char* norm, const char* uplo, const int* n, const int* k,
	const double* ab, const int* ldab, double* work) {
	(void)work;
	try {
		return vcp::rdlansb(*norm, *uplo, *n, *k, ab, *ldab, vlapack_dlapack_detail::current_rounding_mode());
	}
	catch (const std::invalid_argument&) {
		return 0.0;
	}
}

inline double dlanst_(const char* norm, const int* n, const double* d, const double* e) {
	try {
		return vcp::rdlanst(*norm, *n, d, e, vlapack_dlapack_detail::current_rounding_mode());
	}
	catch (const std::invalid_argument&) {
		return 0.0;
	}
}

// ================= LU / 三角行列系 (ipiv は 1-based に変換) =================

inline void dgetf2_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgetf2(*m, *n, a, *lda, ipiv, vlapack_dlapack_detail::current_rounding_mode());
		vlapack_dlapack_detail::ipiv_to_fortran(*m < *n ? *m : *n, ipiv);
	)
}

inline void dgetrf_(const int* m, const int* n, double* a, const int* lda, int* ipiv, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgetrf(*m, *n, a, *lda, ipiv, vlapack_dlapack_detail::current_rounding_mode());
		vlapack_dlapack_detail::ipiv_to_fortran(*m < *n ? *m : *n, ipiv);
	)
}

inline void dgetrs_(const char* trans, const int* n, const int* nrhs, const double* a, const int* lda,
	const int* ipiv, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		const std::vector<int> ip0 = vlapack_dlapack_detail::ipiv_to_c(*n, ipiv);
		*info = vcp::rdgetrs(*trans, *n, *nrhs, a, *lda, ip0.data(), b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv,
	double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgesv(*n, *nrhs, a, *lda, ipiv, b, *ldb, vlapack_dlapack_detail::current_rounding_mode());
		vlapack_dlapack_detail::ipiv_to_fortran(*n, ipiv);
	)
}

inline void dgetri_(const int* n, double* a, const int* lda, const int* ipiv,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 64 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		const std::vector<int> ip0 = vlapack_dlapack_detail::ipiv_to_c(*n, ipiv);
		*info = vcp::rdgetri(*n, a, *lda, ip0.data(), vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dtrti2_(const char* uplo, const char* diag, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdtrti2(*uplo, *diag, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dtrtri_(const char* uplo, const char* diag, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdtrtri(*uplo, *diag, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dtrtrs_(const char* uplo, const char* trans, const char* diag, const int* n, const int* nrhs,
	const double* a, const int* lda, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdtrtrs(*uplo, *trans, *diag, *n, *nrhs, a, *lda, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dlauu2_(const char* uplo, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdlauu2(*uplo, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dlauum_(const char* uplo, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdlauum(*uplo, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= Cholesky 系 =================

inline void dpotf2_(const char* uplo, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpotf2(*uplo, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dpotrf_(const char* uplo, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpotrf(*uplo, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dpotrs_(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda,
	double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpotrs(*uplo, *n, *nrhs, a, *lda, b, *ldb, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dposv_(const char* uplo, const int* n, const int* nrhs, double* a, const int* lda,
	double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdposv(*uplo, *n, *nrhs, a, *lda, b, *ldb, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dpotri_(const char* uplo, const int* n, double* a, const int* lda, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpotri(*uplo, *n, a, *lda, vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= 対称不定値 (Bunch-Kaufman) 系 (ipiv は rdlapack 自体が LAPACK 互換) =================

inline void dsytf2_(const char* uplo, const int* n, double* a, const int* lda, int* ipiv, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsytf2(*uplo, *n, a, *lda, ipiv, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsytrf_(const char* uplo, const int* n, double* a, const int* lda, int* ipiv,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 64 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsytrf(*uplo, *n, a, *lda, ipiv, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsytrs_(const char* uplo, const int* n, const int* nrhs, const double* a, const int* lda,
	const int* ipiv, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsytrs(*uplo, *n, *nrhs, a, *lda, ipiv, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a, const int* lda,
	int* ipiv, double* b, const int* ldb, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 64 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsysv(*uplo, *n, *nrhs, a, *lda, ipiv, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= 帯行列系 (ipiv は 1-based に変換) =================

inline void dgbtf2_(const int* m, const int* n, const int* kl, const int* ku,
	double* ab, const int* ldab, int* ipiv, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgbtf2(*m, *n, *kl, *ku, ab, *ldab, ipiv, vlapack_dlapack_detail::current_rounding_mode());
		vlapack_dlapack_detail::ipiv_to_fortran(*m < *n ? *m : *n, ipiv);
	)
}

inline void dgbtrf_(const int* m, const int* n, const int* kl, const int* ku,
	double* ab, const int* ldab, int* ipiv, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgbtrf(*m, *n, *kl, *ku, ab, *ldab, ipiv, vlapack_dlapack_detail::current_rounding_mode());
		vlapack_dlapack_detail::ipiv_to_fortran(*m < *n ? *m : *n, ipiv);
	)
}

inline void dgbtrs_(const char* trans, const int* n, const int* kl, const int* ku, const int* nrhs,
	const double* ab, const int* ldab, const int* ipiv, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		const std::vector<int> ip0 = vlapack_dlapack_detail::ipiv_to_c(*n, ipiv);
		*info = vcp::rdgbtrs(*trans, *n, *kl, *ku, *nrhs, ab, *ldab, ip0.data(), b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgbsv_(const int* n, const int* kl, const int* ku, const int* nrhs,
	double* ab, const int* ldab, int* ipiv, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgbsv(*n, *kl, *ku, *nrhs, ab, *ldab, ipiv, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
		vlapack_dlapack_detail::ipiv_to_fortran(*n, ipiv);
	)
}

inline void dpbtf2_(const char* uplo, const int* n, const int* kd, double* ab, const int* ldab, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpbtf2(*uplo, *n, *kd, ab, *ldab, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dpbtrf_(const char* uplo, const int* n, const int* kd, double* ab, const int* ldab, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpbtrf(*uplo, *n, *kd, ab, *ldab, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dpbtrs_(const char* uplo, const int* n, const int* kd, const int* nrhs,
	const double* ab, const int* ldab, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpbtrs(*uplo, *n, *kd, *nrhs, ab, *ldab, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dpbsv_(const char* uplo, const int* n, const int* kd, const int* nrhs,
	double* ab, const int* ldab, double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdpbsv(*uplo, *n, *kd, *nrhs, ab, *ldab, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= QR / LQ / QL / 最小二乗系 =================

inline void dgeqr2_(const int* m, const int* n, double* a, const int* lda, double* tau,
	double* work, int* info) {
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgeqr2(*m, *n, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgeqrf_(const int* m, const int* n, double* a, const int* lda, double* tau,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgeqrf(*m, *n, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgelq2_(const int* m, const int* n, double* a, const int* lda, double* tau,
	double* work, int* info) {
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgelq2(*m, *n, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgelqf_(const int* m, const int* n, double* a, const int* lda, double* tau,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*m > 0 ? *m * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgelqf(*m, *n, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dorgqr_(const int* m, const int* n, const int* k, double* a, const int* lda,
	const double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdorgqr(*m, *n, *k, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dorglq_(const int* m, const int* n, const int* k, double* a, const int* lda,
	const double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*m > 0 ? *m * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdorglq(*m, *n, *k, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dorgql_(const int* m, const int* n, const int* k, double* a, const int* lda,
	const double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdorgql(*m, *n, *k, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dormqr_(const char* side, const char* trans, const int* m, const int* n, const int* k,
	const double* a, const int* lda, const double* tau, double* c, const int* ldc,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m > *n ? *m : *n) > 0 ? (*m > *n ? *m : *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdormqr(*side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dormlq_(const char* side, const char* trans, const int* m, const int* n, const int* k,
	const double* a, const int* lda, const double* tau, double* c, const int* ldc,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m > *n ? *m : *n) > 0 ? (*m > *n ? *m : *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdormlq(*side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dormql_(const char* side, const char* trans, const int* m, const int* n, const int* k,
	const double* a, const int* lda, const double* tau, double* c, const int* ldc,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m > *n ? *m : *n) > 0 ? (*m > *n ? *m : *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdormql(*side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgels_(const char* trans, const int* m, const int* n, const int* nrhs,
	double* a, const int* lda, double* b, const int* ldb,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		const int mn = *m < *n ? *m : *n;
		const int mx = *m > *n ? *m : *n;
		const int nb = mx > *nrhs ? mx : *nrhs;
		work[0] = static_cast<double>(mn + nb * 32 > 0 ? mn + nb * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgels(*trans, *m, *n, *nrhs, a, *lda, b, *ldb,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= 対称固有値系 =================

inline void dsytd2_(const char* uplo, const int* n, double* a, const int* lda,
	double* d, double* e, double* tau, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsytd2(*uplo, *n, a, *lda, d, e, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsytrd_(const char* uplo, const int* n, double* a, const int* lda,
	double* d, double* e, double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsytrd(*uplo, *n, a, *lda, d, e, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dorgtr_(const char* uplo, const int* n, double* a, const int* lda, const double* tau,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 1 ? (*n - 1) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdorgtr(*uplo, *n, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dormtr_(const char* side, const char* uplo, const char* trans, const int* m, const int* n,
	const double* a, const int* lda, const double* tau, double* c, const int* ldc,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m > *n ? *m : *n) > 0 ? (*m > *n ? *m : *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdormtr(*side, *uplo, *trans, *m, *n, a, *lda, tau, c, *ldc,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsterf_(const int* n, double* d, double* e, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsterf(*n, d, e, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsteqr_(const char* compz, const int* n, double* d, double* e, double* z, const int* ldz,
	double* work, int* info) {
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsteqr(*compz, *n, d, e, z, *ldz, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsyev_(const char* jobz, const char* uplo, const int* n, double* a, const int* lda,
	double* w, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? 34 * *n : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsyev(*jobz, *uplo, *n, a, *lda, w, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsygs2_(const int* itype, const char* uplo, const int* n, double* a, const int* lda,
	const double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsygs2(*itype, *uplo, *n, a, *lda, b, *ldb, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsygst_(const int* itype, const char* uplo, const int* n, double* a, const int* lda,
	const double* b, const int* ldb, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsygst(*itype, *uplo, *n, a, *lda, b, *ldb, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dsygv_(const int* itype, const char* jobz, const char* uplo, const int* n,
	double* a, const int* lda, double* b, const int* ldb, double* w,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? 34 * *n : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdsygv(*itype, *jobz, *uplo, *n, a, *lda, b, *ldb, w,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= SVD 系 =================

inline void dgebd2_(const int* m, const int* n, double* a, const int* lda,
	double* d, double* e, double* tauq, double* taup, double* work, int* info) {
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgebd2(*m, *n, a, *lda, d, e, tauq, taup, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgebrd_(const int* m, const int* n, double* a, const int* lda,
	double* d, double* e, double* tauq, double* taup, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m + *n) > 0 ? (*m + *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgebrd(*m, *n, a, *lda, d, e, tauq, taup, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dorgbr_(const char* vect, const int* m, const int* n, const int* k, double* a, const int* lda,
	const double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m < *n ? *m : *n) > 0 ? (*m < *n ? *m : *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdorgbr(*vect, *m, *n, *k, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dormbr_(const char* vect, const char* side, const char* trans, const int* m, const int* n,
	const int* k, const double* a, const int* lda, const double* tau, double* c, const int* ldc,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>((*m > *n ? *m : *n) > 0 ? (*m > *n ? *m : *n) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdormbr(*vect, *side, *trans, *m, *n, *k, a, *lda, tau, c, *ldc,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dbdsqr_(const char* uplo, const int* n, const int* ncvt, const int* nru, const int* ncc,
	double* d, double* e, double* vt, const int* ldvt, double* u, const int* ldu,
	double* c, const int* ldc, double* work, int* info) {
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdbdsqr(*uplo, *n, *ncvt, *nru, *ncc, d, e, vt, *ldvt, u, *ldu, c, *ldc,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
	double* a, const int* lda, double* s, double* u, const int* ldu, double* vt, const int* ldvt,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		const int mn = *m < *n ? *m : *n;
		const int mx = *m > *n ? *m : *n;
		const int w1 = 3 * mn + mx;
		const int w2 = 5 * mn;
		work[0] = static_cast<double>((w1 > w2 ? w1 : w2) > 0 ? (w1 > w2 ? w1 : w2) : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgesvd(*jobu, *jobvt, *m, *n, a, *lda, s, u, *ldu, vt, *ldvt,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

// ================= 非対称固有値系 (ilo/ihi は LAPACK と同じ 1-based) =================

inline void dgebal_(const char* job, const int* n, double* a, const int* lda,
	int* ilo, int* ihi, double* scale, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgebal(*job, *n, a, *lda, *ilo, *ihi, scale, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgebak_(const char* job, const char* side, const int* n, const int* ilo, const int* ihi,
	const double* scale, const int* m, double* v, const int* ldv, int* info) {
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgebak(*job, *side, *n, *ilo, *ihi, scale, *m, v, *ldv,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgehd2_(const int* n, const int* ilo, const int* ihi, double* a, const int* lda,
	double* tau, double* work, int* info) {
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgehd2(*n, *ilo, *ihi, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgehrd_(const int* n, const int* ilo, const int* ihi, double* a, const int* lda,
	double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgehrd(*n, *ilo, *ihi, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dorghr_(const int* n, const int* ilo, const int* ihi, double* a, const int* lda,
	const double* tau, double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*ihi - *ilo > 0 ? (*ihi - *ilo) * 32 : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdorghr(*n, *ilo, *ihi, a, *lda, tau, vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dhseqr_(const char* job, const char* compz, const int* n, const int* ilo, const int* ihi,
	double* h, const int* ldh, double* wr, double* wi, double* z, const int* ldz,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? *n : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdhseqr(*job, *compz, *n, *ilo, *ihi, h, *ldh, wr, wi, z, *ldz,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

// howmny = 'A' / 'B' のみ対応 ('S' は info = -2)．select は参照されない
inline void dtrevc_(const char* side, const char* howmny, const int* select, const int* n,
	const double* t, const int* ldt, double* vl, const int* ldvl, double* vr, const int* ldvr,
	const int* mm, int* m, double* work, int* info) {
	(void)select;
	(void)work;
	if (*howmny == 'S' || *howmny == 's') {
		*info = -2;
		return;
	}
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdtrevc(*side, *howmny, *n, t, *ldt, vl, *ldvl, vr, *ldvr, *mm, *m,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

inline void dgeev_(const char* jobvl, const char* jobvr, const int* n, double* a, const int* lda,
	double* wr, double* wi, double* vl, const int* ldvl, double* vr, const int* ldvr,
	double* work, const int* lwork, int* info) {
	*info = 0;
	if (*lwork == -1) {
		work[0] = static_cast<double>(*n > 0 ? 4 * *n : 1);
		return;
	}
	(void)work;
	VLAPACK_DLAPACK_TRY(
		*info = vcp::rdgeev(*jobvl, *jobvr, *n, a, *lda, wr, wi, vl, *ldvl, vr, *ldvr,
			vlapack_dlapack_detail::current_rounding_mode());
	)
}

} // extern "C"

#undef VLAPACK_DLAPACK_TRY

// Fortran object と link する library を作る場合に，1 つの翻訳単位で
// VLAPACK_DLAPACK_EMIT_SYMBOLS を定義して全 symbol を強制発行させる
#ifdef VLAPACK_DLAPACK_EMIT_SYMBOLS
namespace vlapack_dlapack_detail {

typedef void (*emit_fn)();

// __attribute__((used)) がないと未使用 static 配列として最適化で削除され，
// inline 関数の symbol が発行されない
__attribute__((used)) static const emit_fn emit_table[] = {
	reinterpret_cast<emit_fn>(dlamch_),
	reinterpret_cast<emit_fn>(dlapy2_),
	reinterpret_cast<emit_fn>(dlapy3_),
	reinterpret_cast<emit_fn>(dlassq_),
	reinterpret_cast<emit_fn>(dlartg_),
	reinterpret_cast<emit_fn>(dlarfg_),
	reinterpret_cast<emit_fn>(dlascl_),
	reinterpret_cast<emit_fn>(dlaset_),
	reinterpret_cast<emit_fn>(dlacpy_),
	reinterpret_cast<emit_fn>(dlasrt_),
	reinterpret_cast<emit_fn>(dlaswp_),
	reinterpret_cast<emit_fn>(dlange_),
	reinterpret_cast<emit_fn>(dlangb_),
	reinterpret_cast<emit_fn>(dlansy_),
	reinterpret_cast<emit_fn>(dlansb_),
	reinterpret_cast<emit_fn>(dlanst_),
	reinterpret_cast<emit_fn>(dgetf2_),
	reinterpret_cast<emit_fn>(dgetrf_),
	reinterpret_cast<emit_fn>(dgetrs_),
	reinterpret_cast<emit_fn>(dgesv_),
	reinterpret_cast<emit_fn>(dgetri_),
	reinterpret_cast<emit_fn>(dtrti2_),
	reinterpret_cast<emit_fn>(dtrtri_),
	reinterpret_cast<emit_fn>(dtrtrs_),
	reinterpret_cast<emit_fn>(dlauu2_),
	reinterpret_cast<emit_fn>(dlauum_),
	reinterpret_cast<emit_fn>(dpotf2_),
	reinterpret_cast<emit_fn>(dpotrf_),
	reinterpret_cast<emit_fn>(dpotrs_),
	reinterpret_cast<emit_fn>(dposv_),
	reinterpret_cast<emit_fn>(dpotri_),
	reinterpret_cast<emit_fn>(dsytf2_),
	reinterpret_cast<emit_fn>(dsytrf_),
	reinterpret_cast<emit_fn>(dsytrs_),
	reinterpret_cast<emit_fn>(dsysv_),
	reinterpret_cast<emit_fn>(dgbtf2_),
	reinterpret_cast<emit_fn>(dgbtrf_),
	reinterpret_cast<emit_fn>(dgbtrs_),
	reinterpret_cast<emit_fn>(dgbsv_),
	reinterpret_cast<emit_fn>(dpbtf2_),
	reinterpret_cast<emit_fn>(dpbtrf_),
	reinterpret_cast<emit_fn>(dpbtrs_),
	reinterpret_cast<emit_fn>(dpbsv_),
	reinterpret_cast<emit_fn>(dgeqr2_),
	reinterpret_cast<emit_fn>(dgeqrf_),
	reinterpret_cast<emit_fn>(dgelq2_),
	reinterpret_cast<emit_fn>(dgelqf_),
	reinterpret_cast<emit_fn>(dorgqr_),
	reinterpret_cast<emit_fn>(dorglq_),
	reinterpret_cast<emit_fn>(dorgql_),
	reinterpret_cast<emit_fn>(dormqr_),
	reinterpret_cast<emit_fn>(dormlq_),
	reinterpret_cast<emit_fn>(dormql_),
	reinterpret_cast<emit_fn>(dgels_),
	reinterpret_cast<emit_fn>(dsytd2_),
	reinterpret_cast<emit_fn>(dsytrd_),
	reinterpret_cast<emit_fn>(dorgtr_),
	reinterpret_cast<emit_fn>(dormtr_),
	reinterpret_cast<emit_fn>(dsterf_),
	reinterpret_cast<emit_fn>(dsteqr_),
	reinterpret_cast<emit_fn>(dsyev_),
	reinterpret_cast<emit_fn>(dsygs2_),
	reinterpret_cast<emit_fn>(dsygst_),
	reinterpret_cast<emit_fn>(dsygv_),
	reinterpret_cast<emit_fn>(dgebd2_),
	reinterpret_cast<emit_fn>(dgebrd_),
	reinterpret_cast<emit_fn>(dorgbr_),
	reinterpret_cast<emit_fn>(dormbr_),
	reinterpret_cast<emit_fn>(dbdsqr_),
	reinterpret_cast<emit_fn>(dgesvd_),
	reinterpret_cast<emit_fn>(dgebal_),
	reinterpret_cast<emit_fn>(dgebak_),
	reinterpret_cast<emit_fn>(dgehd2_),
	reinterpret_cast<emit_fn>(dgehrd_),
	reinterpret_cast<emit_fn>(dorghr_),
	reinterpret_cast<emit_fn>(dhseqr_),
	reinterpret_cast<emit_fn>(dtrevc_),
	reinterpret_cast<emit_fn>(dgeev_)
};

} // namespace vlapack_dlapack_detail
#endif // VLAPACK_DLAPACK_EMIT_SYMBOLS

#endif // VLAPACK_DLAPACK_HPP
