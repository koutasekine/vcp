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
// - 要素型 T の要件は docs/kernels/tlapack.md (tblas の R1-R4 に加えて
//   R5: std::numeric_limits<T> の特殊化，R6: NaN との == が false)
// - 行列は column-major，index (ipiv 等) の 0-based / 1-based 規約は rdlapack と同一
// - LAPACK の INFO は返り値 (int)．INFO < 0 に相当する引数 error は
//   vcp::invalid_argument を投げる (例外無効時は abort)
// - LAPACK の WORK/LWORK 引数は持たず，作業領域は内部で std::vector<T> により確保する

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

#include <vcp/tblas/tblas.hpp>
#include <vcp/vlapack/rdlapack.hpp>

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

// double without tblas_double.hpp / tlapack_double.hpp: the generic functions below hand a double call over to
// VCP's own double routines (the same ones dblas.hpp / dlapack.hpp call), with the current rounding mode.
// The explicit specializations of tblas_double.hpp / tlapack_double.hpp replace the whole function, so a
// translation unit that includes them does not come here.
namespace tlapack_detail {

// T other than double: no route (the members are never called)
template <typename T>
struct double_route {
	enum { on = 0 };
	template <typename... A> static T lamch(const A&...) { return T(0); }
	template <typename... A> static void laset(const A&...) {}
	template <typename... A> static void lacpy(const A&...) {}
	template <typename... A> static void laswp(const A&...) {}
	template <typename... A> static void lasrt(const A&...) {}
	template <typename... A> static void lassq(const A&...) {}
	template <typename... A> static T lapy2(const A&...) { return T(0); }
	template <typename... A> static T lapy3(const A&...) { return T(0); }
	template <typename... A> static void lartg(const A&...) {}
	template <typename... A> static void lascl(const A&...) {}
	template <typename... A> static T lange(const A&...) { return T(0); }
	template <typename... A> static T lansy(const A&...) { return T(0); }
	template <typename... A> static T lanst(const A&...) { return T(0); }
	template <typename... A> static T langb(const A&...) { return T(0); }
	template <typename... A> static T lansb(const A&...) { return T(0); }
	template <typename... A> static void larfg(const A&...) {}
	template <typename... A> static int getf2(const A&...) { return 0; }
	template <typename... A> static int getrf(const A&...) { return 0; }
	template <typename... A> static int getrs(const A&...) { return 0; }
	template <typename... A> static int gesv(const A&...) { return 0; }
	template <typename... A> static int getri(const A&...) { return 0; }
	template <typename... A> static int trti2(const A&...) { return 0; }
	template <typename... A> static int trtri(const A&...) { return 0; }
	template <typename... A> static int trtrs(const A&...) { return 0; }
	template <typename... A> static int lauu2(const A&...) { return 0; }
	template <typename... A> static int lauum(const A&...) { return 0; }
	template <typename... A> static int potf2(const A&...) { return 0; }
	template <typename... A> static int potrf(const A&...) { return 0; }
	template <typename... A> static int potrs(const A&...) { return 0; }
	template <typename... A> static int posv(const A&...) { return 0; }
	template <typename... A> static int potri(const A&...) { return 0; }
	template <typename... A> static int sytf2(const A&...) { return 0; }
	template <typename... A> static int sytrf(const A&...) { return 0; }
	template <typename... A> static int sytrs(const A&...) { return 0; }
	template <typename... A> static int sysv(const A&...) { return 0; }
	template <typename... A> static int sygs2(const A&...) { return 0; }
	template <typename... A> static int sygst(const A&...) { return 0; }
	template <typename... A> static int sygv(const A&...) { return 0; }
	template <typename... A> static int sygvd(const A&...) { return 0; }
	template <typename... A> static int gbtf2(const A&...) { return 0; }
	template <typename... A> static int gbtrf(const A&...) { return 0; }
	template <typename... A> static int gbtrs(const A&...) { return 0; }
	template <typename... A> static int gbsv(const A&...) { return 0; }
	template <typename... A> static int pbtf2(const A&...) { return 0; }
	template <typename... A> static int pbtrf(const A&...) { return 0; }
	template <typename... A> static int pbtrs(const A&...) { return 0; }
	template <typename... A> static int pbsv(const A&...) { return 0; }
	template <typename... A> static int geqr2(const A&...) { return 0; }
	template <typename... A> static int geqrf(const A&...) { return 0; }
	template <typename... A> static int gelq2(const A&...) { return 0; }
	template <typename... A> static int gelqf(const A&...) { return 0; }
	template <typename... A> static int orgqr(const A&...) { return 0; }
	template <typename... A> static int orglq(const A&...) { return 0; }
	template <typename... A> static int orgql(const A&...) { return 0; }
	template <typename... A> static int ormqr(const A&...) { return 0; }
	template <typename... A> static int ormlq(const A&...) { return 0; }
	template <typename... A> static int ormql(const A&...) { return 0; }
	template <typename... A> static int gels(const A&...) { return 0; }
	template <typename... A> static int sytd2(const A&...) { return 0; }
	template <typename... A> static int sytrd(const A&...) { return 0; }
	template <typename... A> static int orgtr(const A&...) { return 0; }
	template <typename... A> static int ormtr(const A&...) { return 0; }
	template <typename... A> static int sterf(const A&...) { return 0; }
	template <typename... A> static int steqr(const A&...) { return 0; }
	template <typename... A> static int syev(const A&...) { return 0; }
	template <typename... A> static int syevd(const A&...) { return 0; }
	template <typename... A> static int gebd2(const A&...) { return 0; }
	template <typename... A> static int gebrd(const A&...) { return 0; }
	template <typename... A> static int orgbr(const A&...) { return 0; }
	template <typename... A> static int ormbr(const A&...) { return 0; }
	template <typename... A> static int bdsqr(const A&...) { return 0; }
	template <typename... A> static int gesvd(const A&...) { return 0; }
	template <typename... A> static int gebal(const A&...) { return 0; }
	template <typename... A> static int gebak(const A&...) { return 0; }
	template <typename... A> static int gehd2(const A&...) { return 0; }
	template <typename... A> static int gehrd(const A&...) { return 0; }
	template <typename... A> static int orghr(const A&...) { return 0; }
	template <typename... A> static int hseqr(const A&...) { return 0; }
	template <typename... A> static int trevc(const A&...) { return 0; }
	template <typename... A> static int geev(const A&...) { return 0; }
};

template <>
struct double_route<double> {
	enum { on = 1 };
	static double lamch(const char cmach) {
		return vcp::dlamch(cmach);
	}
	static void laset(const char uplo, const int m, const int n, const double& alpha, const double& beta, double* A, const int lda) {
		vcp::dlaset(uplo, m, n, alpha, beta, A, lda);
	}
	static void lacpy(const char uplo, const int m, const int n, const double* A, const int lda, double* B, const int ldb) {
		vcp::dlacpy(uplo, m, n, A, lda, B, ldb);
	}
	static void laswp(const int n, double* A, const int lda, const int k1, const int k2, const int* ipiv, const int incx) {
		vcp::dlaswp(n, A, lda, k1, k2, ipiv, incx);
	}
	static void lasrt(const char id, const int n, double* d) {
		vcp::dlasrt(id, n, d);
	}
	static void lassq(const int n, const double* x, const int incx, double& scale, double& sumsq) {
		vcp::rdlassq(n, x, incx, scale, sumsq, vblas_rdblas_detail::current_rounding_mode());
	}
	static double lapy2(const double& x, const double& y) {
		return vcp::rdlapy2(x, y, vblas_rdblas_detail::current_rounding_mode());
	}
	static double lapy3(const double& x, const double& y, const double& z) {
		return vcp::rdlapy3(x, y, z, vblas_rdblas_detail::current_rounding_mode());
	}
	static void lartg(const double& f, const double& g, double& cs, double& sn, double& r) {
		vcp::rdlartg(f, g, cs, sn, r, vblas_rdblas_detail::current_rounding_mode());
	}
	static void lascl(const char type, const int kl, const int ku, const double& cfrom, const double& cto, const int m, const int n, double* A, const int lda) {
		vcp::rdlascl(type, kl, ku, cfrom, cto, m, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static double lange(const char norm, const int m, const int n, const double* A, const int lda) {
		return vcp::rdlange(norm, m, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static double lansy(const char norm, const char uplo, const int n, const double* A, const int lda) {
		return vcp::rdlansy(norm, uplo, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static double lanst(const char norm, const int n, const double* d, const double* e) {
		return vcp::rdlanst(norm, n, d, e, vblas_rdblas_detail::current_rounding_mode());
	}
	static double langb(const char norm, const int n, const int kl, const int ku, const double* AB, const int ldab) {
		return vcp::rdlangb(norm, n, kl, ku, AB, ldab, vblas_rdblas_detail::current_rounding_mode());
	}
	static double lansb(const char norm, const char uplo, const int n, const int k, const double* AB, const int ldab) {
		return vcp::rdlansb(norm, uplo, n, k, AB, ldab, vblas_rdblas_detail::current_rounding_mode());
	}
	static void larfg(const int n, double& alpha, double* x, const int incx, double& tau) {
		vcp::rdlarfg(n, alpha, x, incx, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int getf2(const int m, const int n, double* A, const int lda, int* ipiv) {
		return vcp::rdgetf2(m, n, A, lda, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int getrf(const int m, const int n, double* A, const int lda, int* ipiv) {
		return vcp::rdgetrf(m, n, A, lda, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int getrs(const char trans, const int n, const int nrhs, const double* A, const int lda, const int* ipiv, double* B, const int ldb) {
		return vcp::rdgetrs(trans, n, nrhs, A, lda, ipiv, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gesv(const int n, const int nrhs, double* A, const int lda, int* ipiv, double* B, const int ldb) {
		return vcp::rdgesv(n, nrhs, A, lda, ipiv, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int getri(const int n, double* A, const int lda, const int* ipiv) {
		return vcp::rdgetri(n, A, lda, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int trti2(const char uplo, const char diag, const int n, double* A, const int lda) {
		return vcp::rdtrti2(uplo, diag, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int trtri(const char uplo, const char diag, const int n, double* A, const int lda) {
		return vcp::rdtrtri(uplo, diag, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int trtrs(const char uplo, const char trans, const char diag, const int n, const int nrhs, const double* A, const int lda, double* B, const int ldb) {
		return vcp::rdtrtrs(uplo, trans, diag, n, nrhs, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int lauu2(const char uplo, const int n, double* A, const int lda) {
		return vcp::rdlauu2(uplo, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int lauum(const char uplo, const int n, double* A, const int lda) {
		return vcp::rdlauum(uplo, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int potf2(const char uplo, const int n, double* A, const int lda) {
		return vcp::rdpotf2(uplo, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int potrf(const char uplo, const int n, double* A, const int lda) {
		return vcp::rdpotrf(uplo, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int potrs(const char uplo, const int n, const int nrhs, const double* A, const int lda, double* B, const int ldb) {
		return vcp::rdpotrs(uplo, n, nrhs, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int posv(const char uplo, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb) {
		return vcp::rdposv(uplo, n, nrhs, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int potri(const char uplo, const int n, double* A, const int lda) {
		return vcp::rdpotri(uplo, n, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sytf2(const char uplo, const int n, double* A, const int lda, int* ipiv) {
		return vcp::rdsytf2(uplo, n, A, lda, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sytrf(const char uplo, const int n, double* A, const int lda, int* ipiv) {
		return vcp::rdsytrf(uplo, n, A, lda, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sytrs(const char uplo, const int n, const int nrhs, const double* A, const int lda, const int* ipiv, double* B, const int ldb) {
		return vcp::rdsytrs(uplo, n, nrhs, A, lda, ipiv, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sysv(const char uplo, const int n, const int nrhs, double* A, const int lda, int* ipiv, double* B, const int ldb) {
		return vcp::rdsysv(uplo, n, nrhs, A, lda, ipiv, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sygs2(const int itype, const char uplo, const int n, double* A, const int lda, const double* B, const int ldb) {
		return vcp::rdsygs2(itype, uplo, n, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sygst(const int itype, const char uplo, const int n, double* A, const int lda, const double* B, const int ldb) {
		return vcp::rdsygst(itype, uplo, n, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sygv(const int itype, const char jobz, const char uplo, const int n, double* A, const int lda, double* B, const int ldb, double* w) {
		// same choice as tlapack_double.hpp: divide and conquer with eigenvectors when its workspace length fits in int
		const long long lwmin = 1LL + 6LL * n + 2LL * n * n;
		if ((jobz == 'V' || jobz == 'v') && lwmin <= static_cast<long long>(std::numeric_limits<int>::max())) {
			return vcp::rdsygvd(itype, jobz, uplo, n, A, lda, B, ldb, w, vblas_rdblas_detail::current_rounding_mode());
		}
		return vcp::rdsygv(itype, jobz, uplo, n, A, lda, B, ldb, w, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sygvd(const int itype, const char jobz, const char uplo, const int n, double* A, const int lda, double* B, const int ldb, double* w) {
		// same choice as tlapack_double.hpp: divide and conquer with eigenvectors when its workspace length fits in int
		const long long lwmin = 1LL + 6LL * n + 2LL * n * n;
		if ((jobz == 'V' || jobz == 'v') && lwmin <= static_cast<long long>(std::numeric_limits<int>::max())) {
			return vcp::rdsygvd(itype, jobz, uplo, n, A, lda, B, ldb, w, vblas_rdblas_detail::current_rounding_mode());
		}
		return vcp::rdsygv(itype, jobz, uplo, n, A, lda, B, ldb, w, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gbtf2(const int m, const int n, const int kl, const int ku, double* AB, const int ldab, int* ipiv) {
		return vcp::rdgbtf2(m, n, kl, ku, AB, ldab, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gbtrf(const int m, const int n, const int kl, const int ku, double* AB, const int ldab, int* ipiv) {
		return vcp::rdgbtrf(m, n, kl, ku, AB, ldab, ipiv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gbtrs(const char trans, const int n, const int kl, const int ku, const int nrhs, const double* AB, const int ldab, const int* ipiv, double* B, const int ldb) {
		return vcp::rdgbtrs(trans, n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gbsv(const int n, const int kl, const int ku, const int nrhs, double* AB, const int ldab, int* ipiv, double* B, const int ldb) {
		return vcp::rdgbsv(n, kl, ku, nrhs, AB, ldab, ipiv, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int pbtf2(const char uplo, const int n, const int kd, double* AB, const int ldab) {
		return vcp::rdpbtf2(uplo, n, kd, AB, ldab, vblas_rdblas_detail::current_rounding_mode());
	}
	static int pbtrf(const char uplo, const int n, const int kd, double* AB, const int ldab) {
		return vcp::rdpbtrf(uplo, n, kd, AB, ldab, vblas_rdblas_detail::current_rounding_mode());
	}
	static int pbtrs(const char uplo, const int n, const int kd, const int nrhs, const double* AB, const int ldab, double* B, const int ldb) {
		return vcp::rdpbtrs(uplo, n, kd, nrhs, AB, ldab, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int pbsv(const char uplo, const int n, const int kd, const int nrhs, double* AB, const int ldab, double* B, const int ldb) {
		return vcp::rdpbsv(uplo, n, kd, nrhs, AB, ldab, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int geqr2(const int m, const int n, double* A, const int lda, double* tau) {
		return vcp::rdgeqr2(m, n, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int geqrf(const int m, const int n, double* A, const int lda, double* tau) {
		return vcp::rdgeqrf(m, n, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gelq2(const int m, const int n, double* A, const int lda, double* tau) {
		return vcp::rdgelq2(m, n, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gelqf(const int m, const int n, double* A, const int lda, double* tau) {
		return vcp::rdgelqf(m, n, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int orgqr(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
		return vcp::rdorgqr(m, n, k, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int orglq(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
		return vcp::rdorglq(m, n, k, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int orgql(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
		return vcp::rdorgql(m, n, k, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int ormqr(const char side, const char trans, const int m, const int n, const int k, const double* A, const int lda, const double* tau, double* C, const int ldc) {
		return vcp::rdormqr(side, trans, m, n, k, A, lda, tau, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static int ormlq(const char side, const char trans, const int m, const int n, const int k, const double* A, const int lda, const double* tau, double* C, const int ldc) {
		return vcp::rdormlq(side, trans, m, n, k, A, lda, tau, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static int ormql(const char side, const char trans, const int m, const int n, const int k, const double* A, const int lda, const double* tau, double* C, const int ldc) {
		return vcp::rdormql(side, trans, m, n, k, A, lda, tau, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gels(const char trans, const int m, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb) {
		return vcp::rdgels(trans, m, n, nrhs, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sytd2(const char uplo, const int n, double* A, const int lda, double* d, double* e, double* tau) {
		return vcp::rdsytd2(uplo, n, A, lda, d, e, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sytrd(const char uplo, const int n, double* A, const int lda, double* d, double* e, double* tau) {
		return vcp::rdsytrd(uplo, n, A, lda, d, e, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int orgtr(const char uplo, const int n, double* A, const int lda, const double* tau) {
		return vcp::rdorgtr(uplo, n, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int ormtr(const char side, const char uplo, const char trans, const int m, const int n, const double* A, const int lda, const double* tau, double* C, const int ldc) {
		return vcp::rdormtr(side, uplo, trans, m, n, A, lda, tau, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static int sterf(const int n, double* d, double* e) {
		return vcp::rdsterf(n, d, e, vblas_rdblas_detail::current_rounding_mode());
	}
	static int steqr(const char compz, const int n, double* d, double* e, double* Z, const int ldz) {
		return vcp::rdsteqr(compz, n, d, e, Z, ldz, vblas_rdblas_detail::current_rounding_mode());
	}
	static int syev(const char jobz, const char uplo, const int n, double* A, const int lda, double* w) {
		// same choice as tlapack_double.hpp: divide and conquer with eigenvectors when its workspace length fits in int
		const long long lwmin = 1LL + 6LL * n + 2LL * n * n;
		if ((jobz == 'V' || jobz == 'v') && lwmin <= static_cast<long long>(std::numeric_limits<int>::max())) {
			return vcp::rdsyevd(jobz, uplo, n, A, lda, w, vblas_rdblas_detail::current_rounding_mode());
		}
		return vcp::rdsyev(jobz, uplo, n, A, lda, w, vblas_rdblas_detail::current_rounding_mode());
	}
	static int syevd(const char jobz, const char uplo, const int n, double* A, const int lda, double* w) {
		// same choice as tlapack_double.hpp: divide and conquer with eigenvectors when its workspace length fits in int
		const long long lwmin = 1LL + 6LL * n + 2LL * n * n;
		if ((jobz == 'V' || jobz == 'v') && lwmin <= static_cast<long long>(std::numeric_limits<int>::max())) {
			return vcp::rdsyevd(jobz, uplo, n, A, lda, w, vblas_rdblas_detail::current_rounding_mode());
		}
		return vcp::rdsyev(jobz, uplo, n, A, lda, w, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gebd2(const int m, const int n, double* A, const int lda, double* d, double* e, double* tauq, double* taup) {
		return vcp::rdgebd2(m, n, A, lda, d, e, tauq, taup, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gebrd(const int m, const int n, double* A, const int lda, double* d, double* e, double* tauq, double* taup) {
		return vcp::rdgebrd(m, n, A, lda, d, e, tauq, taup, vblas_rdblas_detail::current_rounding_mode());
	}
	static int orgbr(const char vect, const int m, const int n, const int k, double* A, const int lda, const double* tau) {
		return vcp::rdorgbr(vect, m, n, k, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int ormbr(const char vect, const char side, const char trans, const int m, const int n, const int k, const double* A, const int lda, const double* tau, double* C, const int ldc) {
		return vcp::rdormbr(vect, side, trans, m, n, k, A, lda, tau, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static int bdsqr(const char uplo, const int n, const int ncvt, const int nru, const int ncc, double* d, double* e, double* VT, const int ldvt, double* U, const int ldu, double* C, const int ldc) {
		return vcp::rdbdsqr(uplo, n, ncvt, nru, ncc, d, e, VT, ldvt, U, ldu, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gesvd(const char jobu, const char jobvt, const int m, const int n, double* A, const int lda, double* s, double* U, const int ldu, double* VT, const int ldvt) {
		return vcp::rdgesvd(jobu, jobvt, m, n, A, lda, s, U, ldu, VT, ldvt, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gebal(const char job, const int n, double* A, const int lda, int& ilo, int& ihi, double* scale) {
		return vcp::rdgebal(job, n, A, lda, ilo, ihi, scale, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gebak(const char job, const char side, const int n, const int ilo, const int ihi, const double* scale, const int m, double* V, const int ldv) {
		return vcp::rdgebak(job, side, n, ilo, ihi, scale, m, V, ldv, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gehd2(const int n, const int ilo, const int ihi, double* A, const int lda, double* tau) {
		return vcp::rdgehd2(n, ilo, ihi, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int gehrd(const int n, const int ilo, const int ihi, double* A, const int lda, double* tau) {
		return vcp::rdgehrd(n, ilo, ihi, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int orghr(const int n, const int ilo, const int ihi, double* A, const int lda, const double* tau) {
		return vcp::rdorghr(n, ilo, ihi, A, lda, tau, vblas_rdblas_detail::current_rounding_mode());
	}
	static int hseqr(const char job, const char compz, const int n, const int ilo, const int ihi, double* H, const int ldh, double* wr, double* wi, double* Z, const int ldz) {
		return vcp::rdhseqr(job, compz, n, ilo, ihi, H, ldh, wr, wi, Z, ldz, vblas_rdblas_detail::current_rounding_mode());
	}
	static int trevc(const char side, const char howmny, const int n, const double* Tf, const int ldt, double* VL, const int ldvl, double* VR, const int ldvr, const int mm, int& m) {
		if (howmny == 'S' || howmny == 's') {
			return -2;   // as dtrevc_ of dlapack.hpp: selected eigenvectors are not supported
		}
		return vcp::rdtrevc(side, howmny, n, Tf, ldt, VL, ldvl, VR, ldvr, mm, m, vblas_rdblas_detail::current_rounding_mode());
	}
	static int geev(const char jobvl, const char jobvr, const int n, double* A, const int lda, double* wr, double* wi, double* VL, const int ldvl, double* VR, const int ldvr) {
		return vcp::rdgeev(jobvl, jobvr, n, A, lda, wr, wi, VL, ldvl, VR, ldvr, vblas_rdblas_detail::current_rounding_mode());
	}
};

} // namespace tlapack_detail

namespace vcp {

// LAPACK の DLAMCH の template 版: machine parameter を std::numeric_limits<T> から
// 構成して返す (要件 R5)．
//   'E' eps (丸め単位 = epsilon()/2)，'S' safe minimum (1/sfmin が overflow しない
//   最小の正数)，'B' base，'P' eps*base = epsilon()，'N' 仮数部 bit 数，
//   'R' rounding (1)，'M' emin，'U' underflow threshold，'L' emax，
//   'O' overflow threshold
template <typename T>
inline T tlamch(const char cmach) {
	if (tlapack_detail::double_route<T>::on) {
		return tlapack_detail::double_route<T>::lamch(cmach);
	}
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
	// VCP-C1 C1-1: the unary plus turns the static constexpr integer members into
	// prvalues, so a T constructor taking a const reference (kv::dd / kv::mpfr<N>)
	// does not odr-use them (no out-of-class definition exists in kv; -O0 link error).
	if (option_is(cmach, 'B')) {
		return T(+lim::radix);
	}
	if (option_is(cmach, 'P')) {
		return lim::epsilon();
	}
	if (option_is(cmach, 'N')) {
		return T(+lim::digits);
	}
	if (option_is(cmach, 'R')) {
		return T(1);
	}
	if (option_is(cmach, 'M')) {
		return T(+lim::min_exponent);
	}
	if (option_is(cmach, 'U')) {
		return lim::min();
	}
	if (option_is(cmach, 'L')) {
		return T(+lim::max_exponent);
	}
	if (option_is(cmach, 'O')) {
		return lim::max();
	}
	tlapack_detail::tlapack_error("tlamch: invalid cmach");
	return T(0);
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_COMMON_HPP
