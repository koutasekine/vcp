// VCP Library
// double specializations for tblas templates.

#pragma once

#ifndef TBLAS_TBLAS_DOUBLE_HPP
#define TBLAS_TBLAS_DOUBLE_HPP

#include <vcp/dblas_dlapack.hpp>
#include <vcp/tblas/tblas.hpp>

namespace vcp {

template <>
inline void tcopy<double>(const int n, const double* x, const int incx, double* y, const int incy) {
	dcopy_(&n, x, &incx, y, &incy);
}

template <>
inline void tswap<double>(const int n, double* x, const int incx, double* y, const int incy) {
	dswap_(&n, x, &incx, y, &incy);
}

template <>
inline int itamax<double>(const int n, const double* x, const int incx) {
	return idamax_(&n, x, &incx) - 1;
}

template <>
inline void tscal<double>(const int n, const double& alpha, double* x, const int incx) {
	dscal_(&n, &alpha, x, &incx);
}

template <>
inline void taxpy<double>(const int n, const double& alpha, const double* x, const int incx, double* y, const int incy) {
	daxpy_(&n, &alpha, x, &incx, y, &incy);
}

template <>
inline double tdot<double>(const int n, const double* x, const int incx, const double* y, const int incy) {
	return ddot_(&n, x, &incx, y, &incy);
}

template <>
inline double tasum<double>(const int n, const double* x, const int incx) {
	return dasum_(&n, x, &incx);
}

template <>
inline double tnrm2<double>(const int n, const double* x, const int incx) {
	return dnrm2_(&n, x, &incx);
}

template <>
inline void trot<double>(const int n, double* x, const int incx, double* y, const int incy, const double& c, const double& s) {
	drot_(&n, x, &incx, y, &incy, &c, &s);
}

template <>
inline void trotg<double>(double* a, double* b, double* c, double* s) {
	drotg_(a, b, c, s);
}

template <>
inline void trotm<double>(const int n, double* x, const int incx, double* y, const int incy, const double* param) {
	drotm_(&n, x, &incx, y, &incy, param);
}

template <>
inline void trotmg<double>(double* d1, double* d2, double* x1, const double& y1, double* param) {
	drotmg_(d1, d2, x1, &y1, param);
}

template <>
inline void tgemv<double>(
	const char trans, const int m, const int n,
	const double& alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double& beta, double* y, const int incy
) {
	dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

template <>
inline void tgbmv<double>(
	const char trans, const int m, const int n, const int kl, const int ku,
	const double& alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double& beta, double* y, const int incy
) {
	dgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

template <>
inline void tsymv<double>(
	const char uplo, const int n,
	const double& alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double& beta, double* y, const int incy
) {
	dsymv_(&uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

template <>
inline void tsbmv<double>(
	const char uplo, const int n, const int k,
	const double& alpha, const double* A, const int lda,
	const double* x, const int incx,
	const double& beta, double* y, const int incy
) {
	dsbmv_(&uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

template <>
inline void tspmv<double>(
	const char uplo, const int n,
	const double& alpha, const double* AP,
	const double* x, const int incx,
	const double& beta, double* y, const int incy
) {
	dspmv_(&uplo, &n, &alpha, AP, x, &incx, &beta, y, &incy);
}

template <>
inline void ttrmv<double>(
	const char uplo, const char trans, const char diag, const int n,
	const double* A, const int lda, double* x, const int incx
) {
	dtrmv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx);
}

template <>
inline void ttrsv<double>(
	const char uplo, const char trans, const char diag, const int n,
	const double* A, const int lda, double* x, const int incx
) {
	dtrsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx);
}

template <>
inline void ttbmv<double>(
	const char uplo, const char trans, const char diag, const int n, const int k,
	const double* A, const int lda, double* x, const int incx
) {
	dtbmv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx);
}

template <>
inline void ttbsv<double>(
	const char uplo, const char trans, const char diag, const int n, const int k,
	const double* A, const int lda, double* x, const int incx
) {
	dtbsv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx);
}

template <>
inline void ttpmv<double>(
	const char uplo, const char trans, const char diag, const int n,
	const double* AP, double* x, const int incx
) {
	dtpmv_(&uplo, &trans, &diag, &n, AP, x, &incx);
}

template <>
inline void ttpsv<double>(
	const char uplo, const char trans, const char diag, const int n,
	const double* AP, double* x, const int incx
) {
	dtpsv_(&uplo, &trans, &diag, &n, AP, x, &incx);
}

template <>
inline void tger<double>(
	const int m, const int n,
	const double& alpha, const double* x, const int incx,
	const double* y, const int incy, double* A, const int lda
) {
	dger_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

template <>
inline void tsyr<double>(
	const char uplo, const int n,
	const double& alpha, const double* x, const int incx,
	double* A, const int lda
) {
	dsyr_(&uplo, &n, &alpha, x, &incx, A, &lda);
}

template <>
inline void tspr<double>(
	const char uplo, const int n,
	const double& alpha, const double* x, const int incx,
	double* AP
) {
	dspr_(&uplo, &n, &alpha, x, &incx, AP);
}

template <>
inline void tsyr2<double>(
	const char uplo, const int n,
	const double& alpha, const double* x, const int incx,
	const double* y, const int incy, double* A, const int lda
) {
	dsyr2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

template <>
inline void tspr2<double>(
	const char uplo, const int n,
	const double& alpha, const double* x, const int incx,
	const double* y, const int incy, double* AP
) {
	dspr2_(&uplo, &n, &alpha, x, &incx, y, &incy, AP);
}

template <>
inline void tgemm<double>(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const double& alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double& beta, double* C, const int ldc
) {
	dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

template <>
inline void tsymm<double>(
	const char side, const char uplo,
	const int m, const int n,
	const double& alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double& beta, double* C, const int ldc
) {
	dsymm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

template <>
inline void tsyrk<double>(
	const char uplo, const char trans,
	const int n, const int k,
	const double& alpha, const double* A, const int lda,
	const double& beta, double* C, const int ldc
) {
	dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc);
}

template <>
inline void tsyr2k<double>(
	const char uplo, const char trans,
	const int n, const int k,
	const double& alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double& beta, double* C, const int ldc
) {
	dsyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

template <>
inline void ttrmm<double>(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const double& alpha, const double* A, const int lda,
	double* B, const int ldb
) {
	dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
}

template <>
inline void ttrsm<double>(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const double& alpha, const double* A, const int lda,
	double* B, const int ldb
) {
	dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
}

template <>
inline void tgemmtr<double>(
	const char uplo, const char transa, const char transb,
	const int n, const int k,
	const double& alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double& beta, double* C, const int ldc
) {
	dgemmtr_(&uplo, &transa, &transb, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

} // namespace vcp

#endif // TBLAS_TBLAS_DOUBLE_HPP
