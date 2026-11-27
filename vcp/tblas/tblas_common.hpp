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

#ifndef TBLAS_TBLAS_COMMON_HPP
#define TBLAS_TBLAS_COMMON_HPP

// tblas (template BLAS) の共通基盤．
// - 要素型 T の要件 (R1: 値セマンティクス / T(int) 構築 / + - * 単項- == != ,
//   R2: 除算, R3: 順序比較, R4: ADL の abs / sqrt) は docs/kernels/tblas.md 参照
// - 丸めモードは一切扱わない (vblas/rdblas との最大の違い)
// - 行列は column-major，index は 0-based
// - OpenMP が有効な場合，演算量が閾値を超える独立 loop のみ並列化する
//   (T の演算が thread-safe な値計算であることを前提とする)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <atomic>
#include <exception>
#include <vcp/error.hpp>
#include <vcp/vblas/rdblas.hpp>

// --- tblas-layer OpenMP guard (SLU-K1F, mirrors spmats SS6.1 / ldbase) ---
// VCP_NOMP             : global kill switch (implies VCP_TBLAS_NOMP)
// VCP_TBLAS_NOMP       : tblas-layer kill switch
// VCP_TBLAS_USE_OPENMP : 0/1, the ONLY symbol implementation code tests
#ifdef VCP_NOMP
#  ifndef VCP_TBLAS_NOMP
#    define VCP_TBLAS_NOMP
#  endif
#endif
#if defined(_OPENMP) && !defined(VCP_TBLAS_NOMP)
#  define VCP_TBLAS_USE_OPENMP 1
#  include <omp.h>
#else
#  define VCP_TBLAS_USE_OPENMP 0
#endif

// OpenMP 並列化を発動するスカラー演算回数の閾値 (T の 1 演算が重い型を想定して低め)
#ifndef TBLAS_OMP_THRESHOLD
#define TBLAS_OMP_THRESHOLD 20000.0
#endif

namespace tblas_detail {

inline void tblas_error(const char* message) {
#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
	vcp::throw_error<vcp::invalid_argument>(message);
#else
	(void)message;
	std::abort();
#endif
}

// BLAS の option 文字比較 (lsame 相当，target は大文字で渡す)
inline bool option_is(const char option, const char target) {
	return option == target || option == target + ('a' - 'A');
}

// 増分 incx の vector の先頭 index (BLAS の負増分規約)
inline int vec_start(const int n, const int inc) {
	return inc > 0 ? 0 : (n - 1) * (-inc);
}

// packed 格納の列先頭 offset (upper: 列 j は AP[j(j+1)/2] から，
// lower: 列 j は AP[jn - j(j-1)/2] から)
inline std::size_t packed_offset_upper(const int j) {
	return static_cast<std::size_t>(j) * (j + 1) / 2;
}

inline std::size_t packed_offset_lower(const int n, const int j) {
	return static_cast<std::size_t>(j) * n - static_cast<std::size_t>(j) * (j - 1) / 2;
}

// 演算量 ops が閾値を超えたら OpenMP 並列化する (if 句用)
inline bool use_parallel(const double ops) {
#if VCP_TBLAS_USE_OPENMP
	return ops >= TBLAS_OMP_THRESHOLD;
#else
	(void)ops;
	return false;
#endif
}

// OpenMP 無効時にも compile できるようにする shim
inline int region_threads() {
#if VCP_TBLAS_USE_OPENMP
	return omp_get_num_threads();
#else
	return 1;
#endif
}

inline int region_thread_id() {
#if VCP_TBLAS_USE_OPENMP
	return omp_get_thread_num();
#else
	return 0;
#endif
}

// y(leny, incy 付き) := beta*y (beta==0 は exact 0 埋め: BLAS 規約で y は読まない)
template <typename T>
inline void scale_vector(const int leny, const T& beta, T* y, const int incy) {
	if (leny <= 0 || beta == T(1)) {
		return;
	}
	int iy = vec_start(leny, incy);
	if (beta == T(0)) {
		for (int i = 0; i < leny; i++) {
			y[iy] = T(0);
			iy += incy;
		}
		return;
	}
	for (int i = 0; i < leny; i++) {
		y[iy] = beta * y[iy];
		iy += incy;
	}
}

// C(m x n, ldc 付き) := beta*C (beta==0 は exact 0 埋め，beta==1 は何もしない)
template <typename T>
inline void scale_matrix(const int m, const int n, const T& beta, T* C, const int ldc) {
	if (m <= 0 || n <= 0 || beta == T(1)) {
		return;
	}
	const bool zero = (beta == T(0));
	{ // OMP-E1: capture exceptions thrown inside the OpenMP region
	std::atomic<bool> omp_caught(false);
	std::exception_ptr omp_eptr;
#if VCP_TBLAS_USE_OPENMP
#pragma omp parallel for schedule(static) if (use_parallel(static_cast<double>(m) * n))
#endif
	for (int j = 0; j < n; j++) {
		try {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		if (zero) {
			for (int i = 0; i < m; i++) {
				c[i] = T(0);
			}
		}
		else {
			for (int i = 0; i < m; i++) {
				c[i] = beta * c[i];
			}
		}
		}
		catch (...) {
			bool omp_expected = false;
			if (omp_caught.compare_exchange_strong(omp_expected, true)) omp_eptr = std::current_exception();
		}
	}
	if (omp_caught.load()) std::rethrow_exception(omp_eptr);
	} // OMP-E1
}

// 三角部分のみ C := beta*C (beta==0 は exact 0 埋め)
template <typename T>
inline void scale_triangle(const bool upper, const int n, const T& beta, T* C, const int ldc) {
	if (n <= 0 || beta == T(1)) {
		return;
	}
	const bool zero = (beta == T(0));
	{ // OMP-E1: capture exceptions thrown inside the OpenMP region
	std::atomic<bool> omp_caught(false);
	std::exception_ptr omp_eptr;
#if VCP_TBLAS_USE_OPENMP
#pragma omp parallel for schedule(static) if (use_parallel(0.5 * n * n))
#endif
	for (int j = 0; j < n; j++) {
		try {
		T* c = C + static_cast<std::size_t>(ldc) * j;
		const int i0 = upper ? 0 : j;
		const int i1 = upper ? j + 1 : n;
		if (zero) {
			for (int i = i0; i < i1; i++) {
				c[i] = T(0);
			}
		}
		else {
			for (int i = i0; i < i1; i++) {
				c[i] = beta * c[i];
			}
		}
		}
		catch (...) {
			bool omp_expected = false;
			if (omp_caught.compare_exchange_strong(omp_expected, true)) omp_eptr = std::current_exception();
		}
	}
	if (omp_caught.load()) std::rethrow_exception(omp_eptr);
	} // OMP-E1
}

// 対称行列 (uplo の三角のみ格納) を full dense (n x n 連続) に展開する (copy のみ)
template <typename T>
inline void expand_symmetric(const bool upper, const int n, const T* A, const int lda, T* dst) {
	for (int j = 0; j < n; j++) {
		T* d = dst + static_cast<std::size_t>(n) * j;
		if (upper) {
			const T* a = A + static_cast<std::size_t>(lda) * j;
			for (int i = 0; i <= j; i++) {
				d[i] = a[i];
			}
			for (int i = j + 1; i < n; i++) {
				d[i] = A[j + static_cast<std::size_t>(lda) * i];
			}
		}
		else {
			for (int i = 0; i < j; i++) {
				d[i] = A[j + static_cast<std::size_t>(lda) * i];
			}
			const T* a = A + static_cast<std::size_t>(lda) * j;
			for (int i = j; i < n; i++) {
				d[i] = a[i];
			}
		}
	}
}

// 三角行列を 0 詰め (unit diag は 1) の full dense (n x n 連続) に展開する (copy のみ)
template <typename T>
inline void expand_triangular(const bool upper, const bool nounit, const int n, const T* A, const int lda, T* dst) {
	for (int j = 0; j < n; j++) {
		const T* a = A + static_cast<std::size_t>(lda) * j;
		T* d = dst + static_cast<std::size_t>(n) * j;
		if (upper) {
			for (int i = 0; i < j; i++) {
				d[i] = a[i];
			}
			d[j] = nounit ? a[j] : T(1);
			for (int i = j + 1; i < n; i++) {
				d[i] = T(0);
			}
		}
		else {
			for (int i = 0; i < j; i++) {
				d[i] = T(0);
			}
			d[j] = nounit ? a[j] : T(1);
			for (int i = j + 1; i < n; i++) {
				d[i] = a[i];
			}
		}
	}
}

} // namespace tblas_detail

// double without tblas_double.hpp / tlapack_double.hpp: the generic functions below hand a double call over to
// VCP's own double routines (the same ones dblas.hpp / dlapack.hpp call), with the current rounding mode.
// The explicit specializations of tblas_double.hpp / tlapack_double.hpp replace the whole function, so a
// translation unit that includes them does not come here.
namespace tblas_detail {

// T other than double: no route (the members are never called)
template <typename T>
struct double_route {
	enum { on = 0 };
	template <typename... A> static void copy(const A&...) {}
	template <typename... A> static void swap(const A&...) {}
	template <typename... A> static int iamax(const A&...) { return 0; }
	template <typename... A> static void scal(const A&...) {}
	template <typename... A> static void axpy(const A&...) {}
	template <typename... A> static T dot(const A&...) { return T(0); }
	template <typename... A> static T asum(const A&...) { return T(0); }
	template <typename... A> static T nrm2(const A&...) { return T(0); }
	template <typename... A> static void rot(const A&...) {}
	template <typename... A> static void rotg(const A&...) {}
	template <typename... A> static void rotm(const A&...) {}
	template <typename... A> static void rotmg(const A&...) {}
	template <typename... A> static void gemv(const A&...) {}
	template <typename... A> static void gbmv(const A&...) {}
	template <typename... A> static void symv(const A&...) {}
	template <typename... A> static void sbmv(const A&...) {}
	template <typename... A> static void spmv(const A&...) {}
	template <typename... A> static void trmv(const A&...) {}
	template <typename... A> static void trsv(const A&...) {}
	template <typename... A> static void tbmv(const A&...) {}
	template <typename... A> static void tbsv(const A&...) {}
	template <typename... A> static void tpmv(const A&...) {}
	template <typename... A> static void tpsv(const A&...) {}
	template <typename... A> static void ger(const A&...) {}
	template <typename... A> static void syr(const A&...) {}
	template <typename... A> static void spr(const A&...) {}
	template <typename... A> static void syr2(const A&...) {}
	template <typename... A> static void spr2(const A&...) {}
	template <typename... A> static void gemm(const A&...) {}
	template <typename... A> static void symm(const A&...) {}
	template <typename... A> static void syrk(const A&...) {}
	template <typename... A> static void syr2k(const A&...) {}
	template <typename... A> static void trmm(const A&...) {}
	template <typename... A> static void trsm(const A&...) {}
	template <typename... A> static void gemmtr(const A&...) {}
};

template <>
struct double_route<double> {
	enum { on = 1 };
	static void copy(const int n, const double* x, const int incx, double* y, const int incy) {
		vcp::dcopy(n, x, incx, y, incy);
	}
	static void swap(const int n, double* x, const int incx, double* y, const int incy) {
		vcp::dswap(n, x, incx, y, incy);
	}
	static int iamax(const int n, const double* x, const int incx) {
		return vcp::idamax(n, x, incx);
	}
	static void scal(const int n, const double& alpha, double* x, const int incx) {
		vcp::rdscal(n, alpha, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void axpy(const int n, const double& alpha, const double* x, const int incx, double* y, const int incy) {
		vcp::rdaxpy(n, alpha, x, incx, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static double dot(const int n, const double* x, const int incx, const double* y, const int incy) {
		return vcp::rddot(n, x, incx, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static double asum(const int n, const double* x, const int incx) {
		return vcp::rdasum(n, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static double nrm2(const int n, const double* x, const int incx) {
		return vcp::rdnrm2(n, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void rot(const int n, double* x, const int incx, double* y, const int incy, const double& c, const double& s) {
		vcp::rdrot(n, x, incx, y, incy, c, s, vblas_rdblas_detail::current_rounding_mode());
	}
	static void rotg(double* a, double* b, double* c, double* s) {
		vcp::rdrotg(a, b, c, s, vblas_rdblas_detail::current_rounding_mode());
	}
	static void rotm(const int n, double* x, const int incx, double* y, const int incy, const double* param) {
		vcp::rdrotm(n, x, incx, y, incy, param, vblas_rdblas_detail::current_rounding_mode());
	}
	static void rotmg(double* d1, double* d2, double* x1, const double& y1, double* param) {
		vcp::rdrotmg(d1, d2, x1, y1, param, vblas_rdblas_detail::current_rounding_mode());
	}
	static void gemv(const char trans, const int m, const int n, const double& alpha, const double* A, const int lda, const double* x, const int incx, const double& beta, double* y, const int incy) {
		vcp::rdgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static void gbmv(const char trans, const int m, const int n, const int kl, const int ku, const double& alpha, const double* A, const int lda, const double* x, const int incx, const double& beta, double* y, const int incy) {
		vcp::rdgbmv(trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static void symv(const char uplo, const int n, const double& alpha, const double* A, const int lda, const double* x, const int incx, const double& beta, double* y, const int incy) {
		vcp::rdsymv(uplo, n, alpha, A, lda, x, incx, beta, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static void sbmv(const char uplo, const int n, const int k, const double& alpha, const double* A, const int lda, const double* x, const int incx, const double& beta, double* y, const int incy) {
		vcp::rdsbmv(uplo, n, k, alpha, A, lda, x, incx, beta, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static void spmv(const char uplo, const int n, const double& alpha, const double* AP, const double* x, const int incx, const double& beta, double* y, const int incy) {
		vcp::rdspmv(uplo, n, alpha, AP, x, incx, beta, y, incy, vblas_rdblas_detail::current_rounding_mode());
	}
	static void trmv(const char uplo, const char trans, const char diag, const int n, const double* A, const int lda, double* x, const int incx) {
		vcp::rdtrmv(uplo, trans, diag, n, A, lda, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void trsv(const char uplo, const char trans, const char diag, const int n, const double* A, const int lda, double* x, const int incx) {
		vcp::rdtrsv(uplo, trans, diag, n, A, lda, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void tbmv(const char uplo, const char trans, const char diag, const int n, const int k, const double* A, const int lda, double* x, const int incx) {
		vcp::rdtbmv(uplo, trans, diag, n, k, A, lda, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void tbsv(const char uplo, const char trans, const char diag, const int n, const int k, const double* A, const int lda, double* x, const int incx) {
		vcp::rdtbsv(uplo, trans, diag, n, k, A, lda, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void tpmv(const char uplo, const char trans, const char diag, const int n, const double* AP, double* x, const int incx) {
		vcp::rdtpmv(uplo, trans, diag, n, AP, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void tpsv(const char uplo, const char trans, const char diag, const int n, const double* AP, double* x, const int incx) {
		vcp::rdtpsv(uplo, trans, diag, n, AP, x, incx, vblas_rdblas_detail::current_rounding_mode());
	}
	static void ger(const int m, const int n, const double& alpha, const double* x, const int incx, const double* y, const int incy, double* A, const int lda) {
		vcp::rdger(m, n, alpha, x, incx, y, incy, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static void syr(const char uplo, const int n, const double& alpha, const double* x, const int incx, double* A, const int lda) {
		vcp::rdsyr(uplo, n, alpha, x, incx, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static void spr(const char uplo, const int n, const double& alpha, const double* x, const int incx, double* AP) {
		vcp::rdspr(uplo, n, alpha, x, incx, AP, vblas_rdblas_detail::current_rounding_mode());
	}
	static void syr2(const char uplo, const int n, const double& alpha, const double* x, const int incx, const double* y, const int incy, double* A, const int lda) {
		vcp::rdsyr2(uplo, n, alpha, x, incx, y, incy, A, lda, vblas_rdblas_detail::current_rounding_mode());
	}
	static void spr2(const char uplo, const int n, const double& alpha, const double* x, const int incx, const double* y, const int incy, double* AP) {
		vcp::rdspr2(uplo, n, alpha, x, incx, y, incy, AP, vblas_rdblas_detail::current_rounding_mode());
	}
	static void gemm(const char transa, const char transb, const int m, const int n, const int k, const double& alpha, const double* A, const int lda, const double* B, const int ldb, const double& beta, double* C, const int ldc) {
		vcp::rdgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static void symm(const char side, const char uplo, const int m, const int n, const double& alpha, const double* A, const int lda, const double* B, const int ldb, const double& beta, double* C, const int ldc) {
		vcp::rdsymm(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static void syrk(const char uplo, const char trans, const int n, const int k, const double& alpha, const double* A, const int lda, const double& beta, double* C, const int ldc) {
		vcp::rdsyrk(uplo, trans, n, k, alpha, A, lda, beta, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static void syr2k(const char uplo, const char trans, const int n, const int k, const double& alpha, const double* A, const int lda, const double* B, const int ldb, const double& beta, double* C, const int ldc) {
		vcp::rdsyr2k(uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
	static void trmm(const char side, const char uplo, const char transa, const char diag, const int m, const int n, const double& alpha, const double* A, const int lda, double* B, const int ldb) {
		vcp::rdtrmm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static void trsm(const char side, const char uplo, const char transa, const char diag, const int m, const int n, const double& alpha, const double* A, const int lda, double* B, const int ldb) {
		vcp::rdtrsm(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb, vblas_rdblas_detail::current_rounding_mode());
	}
	static void gemmtr(const char uplo, const char transa, const char transb, const int n, const int k, const double& alpha, const double* A, const int lda, const double* B, const int ldb, const double& beta, double* C, const int ldc) {
		vcp::rdgemmtr(uplo, transa, transb, n, k, alpha, A, lda, B, ldb, beta, C, ldc, vblas_rdblas_detail::current_rounding_mode());
	}
};

} // namespace tblas_detail

#endif // TBLAS_TBLAS_COMMON_HPP
