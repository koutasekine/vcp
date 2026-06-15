// VCP Library
// kv::dd specializations for tlapack templates using double LAPACK as a
// preconditioner and kv::dd residual refinement.

#pragma once

#ifndef TLAPACK_TLAPACK_DD_HPP
#define TLAPACK_TLAPACK_DD_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <kv/dd.hpp>
#include <vcp/tblas/tblas_double.hpp>
#include <vcp/tblas/tblas_dd.hpp>
#include <vcp/tlapack/tlapack_double.hpp>
#include <vcp/tlapack/tlapack.hpp>

namespace vcp {
namespace tlapack_dd_detail {

inline double to_double(const kv::dd& x) {
	return x.a1;
}

inline double abs_double(const kv::dd& x) {
	return std::fabs(x.a1) + std::fabs(x.a2);
}

inline std::vector<double> to_double_matrix(const int m, const int n, const kv::dd* A, const int lda) {
	std::vector<double> B(static_cast<std::size_t>(m) * n, 0.0);
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			B[i + static_cast<std::size_t>(m) * j] = to_double(A[i + static_cast<std::size_t>(lda) * j]);
		}
	}
	return B;
}

inline std::vector<kv::dd> to_dd_matrix(const int m, const int n, const double* A, const int lda) {
	std::vector<kv::dd> B(static_cast<std::size_t>(m) * n, kv::dd(0.0));
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			B[i + static_cast<std::size_t>(m) * j] = kv::dd(A[i + static_cast<std::size_t>(lda) * j]);
		}
	}
	return B;
}

inline double fro_norm(const int m, const int n, const kv::dd* A, const int lda) {
	long double s = 0.0L;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			const kv::dd& a = A[i + static_cast<std::size_t>(lda) * j];
			const long double v = static_cast<long double>(a.a1) + static_cast<long double>(a.a2);
			s += v * v;
		}
	}
	return static_cast<double>(std::sqrt(s));
}

inline double max_abs(const int m, const int n, const kv::dd* A, const int lda) {
	double r = 0.0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			r = std::max(r, abs_double(A[i + static_cast<std::size_t>(lda) * j]));
		}
	}
	return r;
}

inline void copy_dd_matrix(const int m, const int n, const kv::dd* A, const int lda, kv::dd* B, const int ldb) {
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			B[i + static_cast<std::size_t>(ldb) * j] = A[i + static_cast<std::size_t>(lda) * j];
		}
	}
}

inline void convert_double_to_dd_matrix(const int m, const int n, const double* A, const int lda, kv::dd* B, const int ldb) {
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			B[i + static_cast<std::size_t>(ldb) * j] = kv::dd(A[i + static_cast<std::size_t>(lda) * j]);
		}
	}
}

inline std::vector<kv::dd> reconstruct_lu_matrix(const int n, const kv::dd* LU, const int lda, const int* ipiv) {
	std::vector<kv::dd> L(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> U(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> A(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			const kv::dd a = LU[i + static_cast<std::size_t>(lda) * j];
			if (i > j) {
				L[i + static_cast<std::size_t>(n) * j] = a;
			}
			else {
				U[i + static_cast<std::size_t>(n) * j] = a;
				if (i == j) {
					L[i + static_cast<std::size_t>(n) * j] = kv::dd(1.0);
				}
			}
		}
	}
	tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), L.data(), n, U.data(), n, kv::dd(0.0), A.data(), n);
	for (int j = n - 1; j >= 0; j--) {
		const int jp = ipiv[j];
		if (jp != j) {
			for (int k = 0; k < n; k++) {
				std::swap(A[j + static_cast<std::size_t>(n) * k], A[jp + static_cast<std::size_t>(n) * k]);
			}
		}
	}
	return A;
}

inline void residual_general(
	const char trans, const int n, const int nrhs,
	const kv::dd* A, const int lda,
	const kv::dd* X, const int ldx,
	const kv::dd* B0, const int ldb0,
	kv::dd* R, const int ldr
) {
	copy_dd_matrix(n, nrhs, B0, ldb0, R, ldr);
	tgemm<kv::dd>(trans, 'N', n, nrhs, n, kv::dd(-1.0), A, lda, X, ldx, kv::dd(1.0), R, ldr);
}

inline void correction_by_inverse(
	const char trans, const int n, const int nrhs,
	const std::vector<kv::dd>& Ainv,
	const kv::dd* R, const int ldr,
	kv::dd* X, const int ldx
) {
	const char itrans = tlapack_detail::option_is(trans, 'N') ? 'N' : 'T';
	std::vector<kv::dd> DX(static_cast<std::size_t>(n) * nrhs, kv::dd(0.0));
	tgemm<kv::dd>(itrans, 'N', n, nrhs, n, kv::dd(1.0), Ainv.data(), n, R, ldr, kv::dd(0.0), DX.data(), n);
	for (int j = 0; j < nrhs; j++) {
		for (int i = 0; i < n; i++) {
			X[i + static_cast<std::size_t>(ldx) * j] += DX[i + static_cast<std::size_t>(n) * j];
		}
	}
}

inline int refine_general_with_inverse(
	const char trans, const int n, const int nrhs,
	const kv::dd* Aorig, const int lda,
	const kv::dd* B0, const int ldb0,
	const std::vector<kv::dd>& Ainv,
	kv::dd* X, const int ldx
) {
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	std::vector<kv::dd> R(static_cast<std::size_t>(n) * nrhs, kv::dd(0.0));
	double prev = 0.0;
	const int max_iter = 4;
	for (int iter = 0; iter < max_iter; iter++) {
		residual_general(trans, n, nrhs, Aorig, lda, X, ldx, B0, ldb0, R.data(), n);
		const double rn = fro_norm(n, nrhs, R.data(), n);
		if (rn == 0.0) {
			break;
		}
		if (iter > 0 && rn >= prev * 0.95) {
			break;
		}
		prev = rn;
		correction_by_inverse(trans, n, nrhs, Ainv, R.data(), n, X, ldx);
	}
	return 0;
}

inline void triangular_residual(
	const char uplo, const char trans, const char diag,
	const int n, const int nrhs,
	const kv::dd* A, const int lda,
	const kv::dd* X, const int ldx,
	const kv::dd* B0, const int ldb0,
	kv::dd* R, const int ldr
) {
	const bool upper = tlapack_detail::option_is(uplo, 'U');
	const bool notran = tlapack_detail::option_is(trans, 'N');
	const bool unit = tlapack_detail::option_is(diag, 'U');
	copy_dd_matrix(n, nrhs, B0, ldb0, R, ldr);
	for (int j = 0; j < nrhs; j++) {
		for (int i = 0; i < n; i++) {
			kv::dd s(0.0);
			for (int k = 0; k < n; k++) {
				bool in_tri;
				kv::dd a;
				if (notran) {
					in_tri = upper ? (k >= i) : (k <= i);
					a = (unit && i == k) ? kv::dd(1.0) : A[i + static_cast<std::size_t>(lda) * k];
				}
				else {
					in_tri = upper ? (i >= k) : (i <= k);
					a = (unit && i == k) ? kv::dd(1.0) : A[k + static_cast<std::size_t>(lda) * i];
				}
				if (in_tri) {
					s += a * X[k + static_cast<std::size_t>(ldx) * j];
				}
			}
			R[i + static_cast<std::size_t>(ldr) * j] -= s;
		}
	}
}

inline std::vector<kv::dd> reconstruct_cholesky_matrix(const char uplo, const int n, const kv::dd* C, const int ldc) {
	const bool upper = tlapack_detail::option_is(uplo, 'U');
	std::vector<kv::dd> T(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> A(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			if ((upper && i <= j) || (!upper && i >= j)) {
				T[i + static_cast<std::size_t>(n) * j] = C[i + static_cast<std::size_t>(ldc) * j];
			}
		}
	}
	if (upper) {
		tgemm<kv::dd>('T', 'N', n, n, n, kv::dd(1.0), T.data(), n, T.data(), n, kv::dd(0.0), A.data(), n);
	}
	else {
		tgemm<kv::dd>('N', 'T', n, n, n, kv::dd(1.0), T.data(), n, T.data(), n, kv::dd(0.0), A.data(), n);
	}
	return A;
}

inline int refine_spd_with_cholesky(
	const char uplo, const int n, const int nrhs,
	const kv::dd* Aorig, const int lda,
	const kv::dd* B0, const int ldb0,
	const std::vector<double>& Cd,
	kv::dd* X, const int ldx
) {
	std::vector<kv::dd> R(static_cast<std::size_t>(n) * nrhs, kv::dd(0.0));
	double prev = 0.0;
	for (int iter = 0; iter < 4; iter++) {
		residual_general('N', n, nrhs, Aorig, lda, X, ldx, B0, ldb0, R.data(), n);
		const double rn = fro_norm(n, nrhs, R.data(), n);
		if (rn == 0.0 || (iter > 0 && rn >= prev * 0.95)) {
			break;
		}
		prev = rn;
		std::vector<double> Rd = to_double_matrix(n, nrhs, R.data(), n);
		const int info = tpotrs<double>(uplo, n, nrhs, Cd.data(), n, Rd.data(), n);
		if (info != 0) {
			break;
		}
		for (int j = 0; j < nrhs; j++) {
			for (int i = 0; i < n; i++) {
				X[i + static_cast<std::size_t>(ldx) * j] += kv::dd(Rd[i + static_cast<std::size_t>(n) * j]);
			}
		}
	}
	return 0;
}

inline int refine_symmetric_eigenvectors(
	const char uplo, const int n,
	const std::vector<kv::dd>& Aorig,
	kv::dd* X, const int ldx,
	kv::dd* w
) {
	(void)uplo;
	if (n == 0) {
		return 0;
	}
	const double anorm = std::max(1.0, fro_norm(n, n, Aorig.data(), n));
	std::vector<kv::dd> AX(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> R(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> S(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> E(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> XE(static_cast<std::size_t>(n) * n, kv::dd(0.0));

	double prev = 0.0;
	const int max_iter = 3;
	for (int iter = 0; iter < max_iter; iter++) {
		tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), Aorig.data(), n, X, ldx, kv::dd(0.0), AX.data(), n);
		tgemm<kv::dd>('T', 'N', n, n, n, kv::dd(1.0), X, ldx, AX.data(), n, kv::dd(0.0), S.data(), n);
		tgemm<kv::dd>('T', 'N', n, n, n, kv::dd(-1.0), X, ldx, X, ldx, kv::dd(0.0), R.data(), n);
		for (int i = 0; i < n; i++) {
			R[i + static_cast<std::size_t>(n) * i] += kv::dd(1.0);
		}
		for (int i = 0; i < n; i++) {
			const kv::dd den = kv::dd(1.0) - R[i + static_cast<std::size_t>(n) * i];
			w[i] = S[i + static_cast<std::size_t>(n) * i] / den;
		}
		std::vector<kv::dd> SD = S;
		for (int i = 0; i < n; i++) {
			SD[i + static_cast<std::size_t>(n) * i] -= w[i];
		}
		const double delta = 2.0 * (fro_norm(n, n, SD.data(), n) + anorm * fro_norm(n, n, R.data(), n));
		const double measure = fro_norm(n, n, SD.data(), n) + fro_norm(n, n, R.data(), n);
		if (measure == 0.0) {
			break;
		}
		if (iter > 0 && measure >= prev * 0.95) {
			break;
		}
		prev = measure;

		std::fill(E.begin(), E.end(), kv::dd(0.0));
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				const kv::dd rij = R[i + static_cast<std::size_t>(n) * j];
				if (i == j) {
					E[i + static_cast<std::size_t>(n) * j] = rij / kv::dd(2.0);
				}
				else {
					const kv::dd gap = w[j] - w[i];
					if (abs_double(gap) <= delta) {
						E[i + static_cast<std::size_t>(n) * j] = rij / kv::dd(2.0);
					}
					else {
						const kv::dd sij = S[i + static_cast<std::size_t>(n) * j];
						E[i + static_cast<std::size_t>(n) * j] = (sij + w[j] * rij) / gap;
					}
				}
			}
		}
		tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), X, ldx, E.data(), n, kv::dd(0.0), XE.data(), n);
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				X[i + static_cast<std::size_t>(ldx) * j] += XE[i + static_cast<std::size_t>(n) * j];
			}
		}
	}
	return 0;
}

} // namespace tlapack_dd_detail

template <>
inline int tgetrs<kv::dd>(const char trans, const int n, const int nrhs, const kv::dd* A, const int lda,
	const int* ipiv, kv::dd* B, const int ldb) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	const bool notran = detail::option_is(trans, 'N');
	if (!notran && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::tlapack_error("tgetrs<kv::dd>: invalid trans");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tgetrs<kv::dd>: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	const std::vector<kv::dd> B0(B, B + static_cast<std::size_t>(ldb) * nrhs);
	const std::vector<kv::dd> Aorig = ddd::reconstruct_lu_matrix(n, A, lda, ipiv);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> Bd = ddd::to_double_matrix(n, nrhs, B, ldb);
	const int info = tgetrs<double>(trans, n, nrhs, Ad.data(), n, ipiv, Bd.data(), n);
	if (info != 0) {
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, nrhs, Bd.data(), n, B, ldb);
	const int info_inv = tgetri<double>(n, Ad.data(), n, ipiv);
	if (info_inv != 0) {
		return 0;
	}
	const std::vector<kv::dd> Ainv = ddd::to_dd_matrix(n, n, Ad.data(), n);
	return ddd::refine_general_with_inverse(trans, n, nrhs, Aorig.data(), n, B0.data(), ldb, Ainv, B, ldb);
}

template <>
inline int tgesv<kv::dd>(const int n, const int nrhs, kv::dd* A, const int lda, int* ipiv,
	kv::dd* B, const int ldb) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tgesv<kv::dd>: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	std::vector<kv::dd> Aorig(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	ddd::copy_dd_matrix(n, n, A, lda, Aorig.data(), n);
	const std::vector<kv::dd> B0(B, B + static_cast<std::size_t>(ldb) * nrhs);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> Bd = ddd::to_double_matrix(n, nrhs, B, ldb);
	const int info = tgesv<double>(n, nrhs, Ad.data(), n, ipiv, Bd.data(), n);
	if (info != 0) {
		ddd::convert_double_to_dd_matrix(n, n, Ad.data(), n, A, lda);
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, nrhs, Bd.data(), n, B, ldb);
	std::vector<double> Ainvd = Ad;
	const int info_inv = tgetri<double>(n, Ainvd.data(), n, ipiv);
	if (info_inv == 0) {
		const std::vector<kv::dd> Ainv = ddd::to_dd_matrix(n, n, Ainvd.data(), n);
		ddd::refine_general_with_inverse('N', n, nrhs, Aorig.data(), n, B0.data(), ldb, Ainv, B, ldb);
	}
	ddd::convert_double_to_dd_matrix(n, n, Ad.data(), n, A, lda);
	return 0;
}

template <>
inline int ttrtrs<kv::dd>(const char uplo, const char trans, const char diag, const int n, const int nrhs,
	const kv::dd* A, const int lda, kv::dd* B, const int ldb) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("ttrtrs<kv::dd>: invalid uplo");
	}
	if (!detail::option_is(trans, 'N') && !detail::option_is(trans, 'T') && !detail::option_is(trans, 'C')) {
		detail::tlapack_error("ttrtrs<kv::dd>: invalid trans");
	}
	if (!detail::option_is(diag, 'N') && !detail::option_is(diag, 'U')) {
		detail::tlapack_error("ttrtrs<kv::dd>: invalid diag");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("ttrtrs<kv::dd>: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	const std::vector<kv::dd> B0(B, B + static_cast<std::size_t>(ldb) * nrhs);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> Bd = ddd::to_double_matrix(n, nrhs, B, ldb);
	const int info = ttrtrs<double>(uplo, trans, diag, n, nrhs, Ad.data(), n, Bd.data(), n);
	if (info != 0) {
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, nrhs, Bd.data(), n, B, ldb);
	std::vector<kv::dd> R(static_cast<std::size_t>(n) * nrhs, kv::dd(0.0));
	double prev = 0.0;
	for (int iter = 0; iter < 4; iter++) {
		ddd::triangular_residual(uplo, trans, diag, n, nrhs, A, lda, B, ldb, B0.data(), ldb, R.data(), n);
		const double rn = ddd::fro_norm(n, nrhs, R.data(), n);
		if (rn == 0.0 || (iter > 0 && rn >= prev * 0.95)) {
			break;
		}
		prev = rn;
		std::vector<double> Rd = ddd::to_double_matrix(n, nrhs, R.data(), n);
		const int cinfo = ttrtrs<double>(uplo, trans, diag, n, nrhs, Ad.data(), n, Rd.data(), n);
		if (cinfo != 0) {
			break;
		}
		for (int j = 0; j < nrhs; j++) {
			for (int i = 0; i < n; i++) {
				B[i + static_cast<std::size_t>(ldb) * j] += kv::dd(Rd[i + static_cast<std::size_t>(n) * j]);
			}
		}
	}
	return 0;
}

template <>
inline int tpotrs<kv::dd>(const char uplo, const int n, const int nrhs, const kv::dd* A, const int lda,
	kv::dd* B, const int ldb) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotrs<kv::dd>: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tpotrs<kv::dd>: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	const std::vector<kv::dd> B0(B, B + static_cast<std::size_t>(ldb) * nrhs);
	const std::vector<kv::dd> Aorig = ddd::reconstruct_cholesky_matrix(uplo, n, A, lda);
	std::vector<double> Cd = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> Bd = ddd::to_double_matrix(n, nrhs, B, ldb);
	const int info = tpotrs<double>(uplo, n, nrhs, Cd.data(), n, Bd.data(), n);
	if (info != 0) {
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, nrhs, Bd.data(), n, B, ldb);
	return ddd::refine_spd_with_cholesky(uplo, n, nrhs, Aorig.data(), n, B0.data(), ldb, Cd, B, ldb);
}

template <>
inline int tposv<kv::dd>(const char uplo, const int n, const int nrhs, kv::dd* A, const int lda,
	kv::dd* B, const int ldb) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tposv<kv::dd>: invalid uplo");
	}
	if (n < 0 || nrhs < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tposv<kv::dd>: invalid argument");
	}
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	std::vector<kv::dd> Aorig(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	ddd::copy_dd_matrix(n, n, A, lda, Aorig.data(), n);
	const std::vector<kv::dd> B0(B, B + static_cast<std::size_t>(ldb) * nrhs);
	std::vector<double> Cd = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> Bd = ddd::to_double_matrix(n, nrhs, B, ldb);
	const int info = tposv<double>(uplo, n, nrhs, Cd.data(), n, Bd.data(), n);
	ddd::convert_double_to_dd_matrix(n, n, Cd.data(), n, A, lda);
	if (info != 0) {
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, nrhs, Bd.data(), n, B, ldb);
	return ddd::refine_spd_with_cholesky(uplo, n, nrhs, Aorig.data(), n, B0.data(), ldb, Cd, B, ldb);
}

template <>
inline int tsyev<kv::dd>(const char jobz, const char uplo, const int n, kv::dd* A, const int lda, kv::dd* w) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::tlapack_error("tsyev<kv::dd>: invalid jobz");
	}
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsyev<kv::dd>: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tsyev<kv::dd>: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	std::vector<kv::dd> Aorig(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	ddd::copy_dd_matrix(n, n, A, lda, Aorig.data(), n);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> wd(static_cast<std::size_t>(n), 0.0);
	const int info = tsyev<double>(jobz, uplo, n, Ad.data(), n, wd.data());
	if (info != 0) {
		return info;
	}
	for (int i = 0; i < n; i++) {
		w[i] = kv::dd(wd[i]);
	}
	if (!wantz) {
		return 0;
	}
	ddd::convert_double_to_dd_matrix(n, n, Ad.data(), n, A, lda);
	return ddd::refine_symmetric_eigenvectors(uplo, n, Aorig, A, lda, w);
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_DD_HPP
