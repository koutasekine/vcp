// VCP Library
// kv::dd specializations for tlapack templates using double LAPACK as a
// preconditioner and kv::dd residual refinement.

#pragma once

#ifndef TLAPACK_TLAPACK_DD_HPP
#define TLAPACK_TLAPACK_DD_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
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

// Frobenius norm.  The entries are divided by the largest absolute value
// before they are squared, so that the sum of squares neither overflows nor
// underflows.
inline double fro_norm(const int m, const int n, const kv::dd* A, const int lda) {
	double amax = 0.0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			const kv::dd& a = A[i + static_cast<std::size_t>(lda) * j];
			amax = std::max(amax, std::fabs(a.a1 + a.a2));
		}
	}
	if (amax == 0.0) {
		return 0.0;
	}
	double s = 0.0;
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			const kv::dd& a = A[i + static_cast<std::size_t>(lda) * j];
			const double v = (a.a1 + a.a2) / amax;
			s += v * v;
		}
	}
	return amax * std::sqrt(s);
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

inline std::vector<kv::dd> symmetric_full_matrix(const char uplo, const int n, const kv::dd* A, const int lda) {
	const bool upper = tlapack_detail::option_is(uplo, 'U');
	std::vector<kv::dd> B(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			if ((upper && i <= j) || (!upper && i >= j)) {
				B[i + static_cast<std::size_t>(n) * j] = A[i + static_cast<std::size_t>(lda) * j];
			}
			else {
				B[i + static_cast<std::size_t>(n) * j] = A[j + static_cast<std::size_t>(lda) * i];
			}
		}
	}
	return B;
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

inline std::map<const kv::dd*, std::vector<kv::dd> >& lu_original_cache() {
	static std::map<const kv::dd*, std::vector<kv::dd> > cache;
	return cache;
}

inline void save_lu_original(const kv::dd* key, const int m, const int n, const kv::dd* A, const int lda) {
	std::vector<kv::dd> Aorig(static_cast<std::size_t>(m) * n, kv::dd(0.0));
	copy_dd_matrix(m, n, A, lda, Aorig.data(), m);
	lu_original_cache()[key] = Aorig;
}

inline bool load_lu_original(const kv::dd* key, const int n, std::vector<kv::dd>& Aorig) {
	std::map<const kv::dd*, std::vector<kv::dd> >& cache = lu_original_cache();
	std::map<const kv::dd*, std::vector<kv::dd> >::const_iterator it = cache.find(key);
	if (it == cache.end() || static_cast<int>(it->second.size()) != n * n) {
		return false;
	}
	Aorig = it->second;
	return true;
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

inline int refine_general_with_double_lu(
	const char trans, const int n, const int nrhs,
	const kv::dd* Aorig, const int lda,
	const kv::dd* B0, const int ldb0,
	const std::vector<double>& LUd,
	const int* ipiv,
	kv::dd* X, const int ldx
) {
	if (n == 0 || nrhs == 0) {
		return 0;
	}
	std::vector<kv::dd> R(static_cast<std::size_t>(n) * nrhs, kv::dd(0.0));
	std::vector<double> Rd(static_cast<std::size_t>(n) * nrhs, 0.0);
	double prev = std::numeric_limits<double>::infinity();
	const int max_iter = 10;
	for (int iter = 0; iter < max_iter; iter++) {
		residual_general(trans, n, nrhs, Aorig, lda, X, ldx, B0, ldb0, R.data(), n);
		double normr = 0.0;
		double normx = 0.0;
		for (int j = 0; j < nrhs; j++) {
			for (int i = 0; i < n; i++) {
				const std::size_t ridx = i + static_cast<std::size_t>(n) * j;
				const std::size_t xidx = i + static_cast<std::size_t>(ldx) * j;
				Rd[ridx] = R[ridx].a1;
				normr = std::max(normr, std::fabs(Rd[ridx]));
				normx = std::max(normx, abs_double(X[xidx]));
			}
		}
		if (normr == 0.0 || normr <= 0.25e-32 * std::max(1.0, normx) || normr >= 0.5 * prev) {
			break;
		}
		prev = normr;
		std::vector<double> Cd = LUd;
		const int info = tgetrs<double>(trans, n, nrhs, Cd.data(), n, ipiv, Rd.data(), n);
		if (info != 0) {
			return info;
		}
		for (int j = 0; j < nrhs; j++) {
			for (int i = 0; i < n; i++) {
				X[i + static_cast<std::size_t>(ldx) * j] += kv::dd(Rd[i + static_cast<std::size_t>(n) * j]);
			}
		}
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

inline int refine_cholesky_factor(
	const char uplo, const int n,
	const std::vector<kv::dd>& Aorig,
	kv::dd* C, const int ldc
) {
	const bool upper = tlapack_detail::option_is(uplo, 'U');
	std::vector<kv::dd> T(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> P(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> E(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> W(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<kv::dd> F(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	std::vector<double> Cd(static_cast<std::size_t>(n) * n, 0.0);
	std::vector<double> G(static_cast<std::size_t>(n) * n, 0.0);

	const int max_iter = 2;
	for (int iter = 0; iter < max_iter; iter++) {
		std::fill(T.begin(), T.end(), kv::dd(0.0));
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				if ((upper && i <= j) || (!upper && i >= j)) {
					T[i + static_cast<std::size_t>(n) * j] = C[i + static_cast<std::size_t>(ldc) * j];
				}
			}
		}
		if (upper) {
			tgemm<kv::dd>('T', 'N', n, n, n, kv::dd(1.0), T.data(), n, T.data(), n, kv::dd(0.0), P.data(), n);
		}
		else {
			tgemm<kv::dd>('N', 'T', n, n, n, kv::dd(1.0), T.data(), n, T.data(), n, kv::dd(0.0), P.data(), n);
		}
		for (int idx = 0; idx < n * n; idx++) {
			E[static_cast<std::size_t>(idx)] = Aorig[static_cast<std::size_t>(idx)] - P[static_cast<std::size_t>(idx)];
			Cd[static_cast<std::size_t>(idx)] = T[static_cast<std::size_t>(idx)].a1;
			G[static_cast<std::size_t>(idx)] = E[static_cast<std::size_t>(idx)].a1;
		}
		if (max_abs(n, n, E.data(), n) == 0.0) {
			break;
		}

		if (upper) {
			ttrsm<double>('L', 'U', 'T', 'N', n, n, 1.0, Cd.data(), n, G.data(), n);
			ttrsm<double>('R', 'U', 'N', 'N', n, n, 1.0, Cd.data(), n, G.data(), n);
			std::fill(W.begin(), W.end(), kv::dd(0.0));
			for (int j = 0; j < n; j++) {
				for (int i = 0; i <= j; i++) {
					W[i + static_cast<std::size_t>(n) * j] = kv::dd(G[i + static_cast<std::size_t>(n) * j] * (i == j ? 0.5 : 1.0));
				}
			}
			tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), W.data(), n, T.data(), n, kv::dd(0.0), F.data(), n);
		}
		else {
			ttrsm<double>('L', 'L', 'N', 'N', n, n, 1.0, Cd.data(), n, G.data(), n);
			ttrsm<double>('R', 'L', 'T', 'N', n, n, 1.0, Cd.data(), n, G.data(), n);
			std::fill(W.begin(), W.end(), kv::dd(0.0));
			for (int j = 0; j < n; j++) {
				for (int i = j; i < n; i++) {
					W[i + static_cast<std::size_t>(n) * j] = kv::dd(G[i + static_cast<std::size_t>(n) * j] * (i == j ? 0.5 : 1.0));
				}
			}
			tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), T.data(), n, W.data(), n, kv::dd(0.0), F.data(), n);
		}
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				if ((upper && i <= j) || (!upper && i >= j)) {
					C[i + static_cast<std::size_t>(ldc) * j] += F[i + static_cast<std::size_t>(n) * j];
				}
			}
		}
	}
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < n; i++) {
			if ((upper && i > j) || (!upper && i < j)) {
				C[i + static_cast<std::size_t>(ldc) * j] = kv::dd(0.0);
			}
		}
	}
	return 0;
}

// Indices 0..l-1 ordered so that w is ascending.
inline std::vector<int> ascending_index(const int l, const kv::dd* w) {
	std::vector<int> idx(static_cast<std::size_t>(l), 0);
	for (int i = 0; i < l; i++) {
		idx[i] = i;
	}
	std::stable_sort(idx.begin(), idx.end(), [w](const int a, const int b) { return w[a] < w[b]; });
	return idx;
}

// Algorithm 1 (RefSyEv) of T. Ogita, K. Aishima, "Iterative refinement for
// symmetric eigenvalue decomposition II: clustered eigenvalues", Japan J.
// Indust. Appl. Math. 36 (2019), 435-459.
// A is n x n, X is n x l (l <= n) and is replaced by X' = X + X E.  B = 0:
// R = I - X^T X.  B != 0 (Remark 3): R = I - X^T B X.  lam receives the
// approximate eigenvalues, enorm = ||E||_F, omega the criterion of line 6
// (Frobenius norms).
inline void refsyev(
	const int n, const int l,
	const kv::dd* A, const kv::dd* B,
	kv::dd* X, const int ldx,
	std::vector<kv::dd>& lam, double& enorm, double& omega
) {
	const std::size_t N = static_cast<std::size_t>(n);
	const std::size_t L = static_cast<std::size_t>(l);
	std::vector<kv::dd> AX(N * L, kv::dd(0.0));
	std::vector<kv::dd> R(L * L, kv::dd(0.0));
	std::vector<kv::dd> S(L * L, kv::dd(0.0));
	std::vector<kv::dd> E(L * L, kv::dd(0.0));
	std::vector<kv::dd> XE(N * L, kv::dd(0.0));

	// line 2: R = I - X^T X (or I - X^T B X), line 3: S = X^T A X
	if (B != 0) {
		tgemm<kv::dd>('N', 'N', n, l, n, kv::dd(1.0), B, n, X, ldx, kv::dd(0.0), AX.data(), n);
		tgemm<kv::dd>('T', 'N', l, l, n, kv::dd(-1.0), X, ldx, AX.data(), n, kv::dd(0.0), R.data(), l);
	}
	else {
		tgemm<kv::dd>('T', 'N', l, l, n, kv::dd(-1.0), X, ldx, X, ldx, kv::dd(0.0), R.data(), l);
	}
	for (int i = 0; i < l; i++) {
		R[i + L * i] += kv::dd(1.0);
	}
	tgemm<kv::dd>('N', 'N', n, l, n, kv::dd(1.0), A, n, X, ldx, kv::dd(0.0), AX.data(), n);
	tgemm<kv::dd>('T', 'N', l, l, n, kv::dd(1.0), X, ldx, AX.data(), n, kv::dd(0.0), S.data(), l);

	// line 4: lambda_i = s_ii / (1 - r_ii)
	lam.assign(L, kv::dd(0.0));
	for (int i = 0; i < l; i++) {
		lam[i] = S[i + L * i] / (kv::dd(1.0) - R[i + L * i]);
	}

	// line 6: omega = 2 (||S - D|| + ||A|| ||R||)
	std::vector<kv::dd> SD = S;
	for (int i = 0; i < l; i++) {
		SD[i + L * i] -= lam[i];
	}
	omega = 2.0 * (fro_norm(l, l, SD.data(), l) + fro_norm(n, n, A, n) * fro_norm(l, l, R.data(), l));

	// line 7
	for (int j = 0; j < l; j++) {
		for (int i = 0; i < l; i++) {
			const kv::dd rij = R[i + L * j];
			const kv::dd gap = lam[j] - lam[i];
			if (abs_double(gap) > omega) {
				E[i + L * j] = (S[i + L * j] + lam[j] * rij) / gap;
			}
			else {
				E[i + L * j] = rij / kv::dd(2.0);
			}
		}
	}
	enorm = fro_norm(l, l, E.data(), l);

	// line 8: X' = X + X E
	tgemm<kv::dd>('N', 'N', n, l, l, kv::dd(1.0), X, ldx, E.data(), l, kv::dd(0.0), XE.data(), n);
	for (int j = 0; j < l; j++) {
		for (int i = 0; i < n; i++) {
			X[i + static_cast<std::size_t>(ldx) * j] += XE[i + N * j];
		}
	}
}

// Algorithm 2 (RefSyEvCL) of the same paper, one call.  X (n x n) is replaced
// by X'.  enorm receives ||E||_F of line 2.  Returns false when the algorithm
// returned at line 3 or line 15 (improvement cannot be expected) or when the
// double eigensolver of line 11 failed; X then holds the X' of that moment.
inline bool refsyevcl(
	const int n,
	const kv::dd* A, const kv::dd* B,
	kv::dd* X, const int ldx,
	std::vector<kv::dd>& lam, double& enorm
) {
	const std::size_t N = static_cast<std::size_t>(n);
	double omega = 0.0;
	refsyev(n, n, A, B, X, ldx, lam, enorm, omega);   // line 2
	if (!(enorm < 1.0)) {                             // line 3
		return false;
	}
	// line 4: the index sets J_k of (27) are the runs of the ascending
	// eigenvalues whose neighbouring differences are at most omega
	const std::vector<int> idx = ascending_index(n, lam.data());
	int p = 0;
	while (p < n) {
		int q = p;
		while (q + 1 < n && abs_double(lam[idx[q + 1]] - lam[idx[q]]) <= omega) {
			q++;
		}
		const int nk = q - p + 1;
		if (nk >= 2) {
			const std::size_t NK = static_cast<std::size_t>(nk);
			// line 6: V_k = X'(:, J_k)
			std::vector<kv::dd> V(N * NK, kv::dd(0.0));
			for (int c = 0; c < nk; c++) {
				for (int i = 0; i < n; i++) {
					V[i + N * c] = X[i + static_cast<std::size_t>(ldx) * idx[p + c]];
				}
			}
			// line 7: mu_k, line 8: A_k = A - mu_k I (Remark 3: A - mu_k B)
			const kv::dd mu = (lam[idx[p]] + lam[idx[q]]) / kv::dd(2.0);
			std::vector<kv::dd> Ak(A, A + N * N);
			if (B != 0) {
				for (std::size_t x = 0; x < N * N; x++) {
					Ak[x] -= mu * B[x];
				}
			}
			else {
				for (int i = 0; i < n; i++) {
					Ak[i + N * i] -= mu;
				}
			}
			// line 9: C_k = V_k^T A_k V_k, line 10: T_k = fl(C_k), line 11: W_k = eig(T_k)
			std::vector<kv::dd> AV(N * NK, kv::dd(0.0));
			std::vector<kv::dd> C(NK * NK, kv::dd(0.0));
			tgemm<kv::dd>('N', 'N', n, nk, n, kv::dd(1.0), Ak.data(), n, V.data(), n, kv::dd(0.0), AV.data(), n);
			tgemm<kv::dd>('T', 'N', nk, nk, n, kv::dd(1.0), V.data(), n, AV.data(), n, kv::dd(0.0), C.data(), nk);
			std::vector<double> Tk = to_double_matrix(nk, nk, C.data(), nk);
			std::vector<double> wk(NK, 0.0);
			if (tsyev<double>('V', 'U', nk, Tk.data(), nk, wk.data()) != 0) {
				return false;
			}
			// line 12: V_k = V_k W_k
			const std::vector<kv::dd> W = to_dd_matrix(nk, nk, Tk.data(), nk);
			tgemm<kv::dd>('N', 'N', n, nk, nk, kv::dd(1.0), V.data(), n, W.data(), nk, kv::dd(0.0), AV.data(), n);
			V.swap(AV);
			// lines 13-17.  Unlike the wording of Remark 3, B is also used here (R = I - V_k^T B V_k).
			// In addition to the condition of line 17 the loop ends when ||E_k|| no
			// longer decreases, and after 10 passes.
			std::vector<kv::dd> lamk;
			double ek = 0.0;
			double ekprev = 0.0;
			double omegak = 0.0;
			for (int pass = 0; pass < 10; pass++) {
				refsyev(n, nk, Ak.data(), B, V.data(), n, lamk, ek, omegak);   // line 14
				if (!(ek < 1.0)) {                                           // line 15
					return false;
				}
				if (ek <= enorm) {                                           // line 17
					break;
				}
				if (pass > 0 && !(ek < ekprev)) {
					break;
				}
				ekprev = ek;
			}
			// line 18
			for (int c = 0; c < nk; c++) {
				for (int i = 0; i < n; i++) {
					X[i + static_cast<std::size_t>(ldx) * idx[p + c]] = V[i + N * c];
				}
			}
		}
		p = q + 1;
	}
	return true;
}

// Refinement of the eigenvectors X (n x n) of A x = lambda x (B = 0) or
// A x = lambda B x: Algorithm 2 at most 6 times (it ends earlier when
// ||E|| < 2^-52 after an update, or when Algorithm 2 returned without
// improvement), orthonormalization steps, Rayleigh quotients, ascending order.
inline int refine_eigenvectors_clustered(
	const int n,
	const kv::dd* A, const kv::dd* B,
	kv::dd* X, const int ldx,
	kv::dd* w
) {
	if (n == 0) {
		return 0;
	}
	const std::size_t N = static_cast<std::size_t>(n);
	std::vector<kv::dd> lam;
	std::vector<kv::dd> Xkeep(N * N, kv::dd(0.0));
	const int max_iter = 6;
	for (int iter = 0; iter < max_iter; iter++) {
		double enorm = 0.0;
		const bool go_on = refsyevcl(n, A, B, X, ldx, lam, enorm);
		if (!go_on) {
			break;
		}
		if (enorm < std::ldexp(1.0, -52)) {
			// the next correction would be below the unit of kv::dd
			break;
		}
	}

	// Steps X = X (I + R / 2), R = I - X^T B X: formula (25) of the paper used for all pairs.
	// Repeated (at most 4 times) until the R that was used satisfies ||R||_F <= 2^-45.
	std::vector<kv::dd> AX(N * N, kv::dd(0.0));
	std::vector<kv::dd> R(N * N, kv::dd(0.0));
	for (int step = 0; step < 4; step++) {
		if (B != 0) {
			tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), B, n, X, ldx, kv::dd(0.0), AX.data(), n);
			tgemm<kv::dd>('T', 'N', n, n, n, kv::dd(-1.0), X, ldx, AX.data(), n, kv::dd(0.0), R.data(), n);
		}
		else {
			tgemm<kv::dd>('T', 'N', n, n, n, kv::dd(-1.0), X, ldx, X, ldx, kv::dd(0.0), R.data(), n);
		}
		for (int i = 0; i < n; i++) {
			R[i + N * i] += kv::dd(1.0);
		}
		const double rnorm = fro_norm(n, n, R.data(), n);
		for (std::size_t x = 0; x < N * N; x++) {
			R[x] = R[x] / kv::dd(2.0);
		}
		tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), X, ldx, R.data(), n, kv::dd(0.0), AX.data(), n);
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				X[i + static_cast<std::size_t>(ldx) * j] += AX[i + N * j];
			}
		}
		if (rnorm <= std::ldexp(1.0, -45)) {
			break;
		}
	}

	// Rayleigh quotients of the returned X
	std::vector<kv::dd> BX;
	tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), A, n, X, ldx, kv::dd(0.0), AX.data(), n);
	if (B != 0) {
		BX.assign(N * N, kv::dd(0.0));
		tgemm<kv::dd>('N', 'N', n, n, n, kv::dd(1.0), B, n, X, ldx, kv::dd(0.0), BX.data(), n);
	}
	lam.assign(N, kv::dd(0.0));
	for (int j = 0; j < n; j++) {
		kv::dd num(0.0);
		kv::dd den(0.0);
		for (int i = 0; i < n; i++) {
			const kv::dd x = X[i + static_cast<std::size_t>(ldx) * j];
			num += x * AX[i + N * j];
			den += x * (B != 0 ? BX[i + N * j] : x);
		}
		lam[j] = num / den;
	}

	// ascending order
	const std::vector<int> idx = ascending_index(n, lam.data());
	copy_dd_matrix(n, n, X, ldx, Xkeep.data(), n);
	for (int c = 0; c < n; c++) {
		w[c] = lam[idx[c]];
		for (int i = 0; i < n; i++) {
			X[i + static_cast<std::size_t>(ldx) * c] = Xkeep[i + N * idx[c]];
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
	return refine_eigenvectors_clustered(n, Aorig.data(), static_cast<const kv::dd*>(0), X, ldx, w);
}

inline int refine_generalized_symmetric_eigenvectors(
	const int n,
	const std::vector<kv::dd>& Aorig,
	const std::vector<kv::dd>& Borig,
	kv::dd* X, const int ldx,
	kv::dd* w
) {
	return refine_eigenvectors_clustered(n, Aorig.data(), Borig.data(), X, ldx, w);
}

} // namespace tlapack_dd_detail

template <>
inline int tgetrf<kv::dd>(const int m, const int n, kv::dd* A, const int lda, int* ipiv) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (m < 0 || n < 0 || lda < std::max(1, m)) {
		detail::tlapack_error("tgetrf<kv::dd>: invalid argument");
	}
	if (m == 0 || n == 0) {
		return 0;
	}
	ddd::save_lu_original(A, m, n, A, lda);
	std::vector<double> Ad = ddd::to_double_matrix(m, n, A, lda);
	const int info = tgetrf<double>(m, n, Ad.data(), m, ipiv);
	ddd::convert_double_to_dd_matrix(m, n, Ad.data(), m, A, lda);
	return info;
}

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
	std::vector<kv::dd> Aorig;
	if (!ddd::load_lu_original(A, n, Aorig)) {
		Aorig = ddd::reconstruct_lu_matrix(n, A, lda, ipiv);
	}
	std::vector<double> Ad = ddd::to_double_matrix(n, n, A, lda);
	std::vector<double> Bd = ddd::to_double_matrix(n, nrhs, B, ldb);
	const int info = tgetrs<double>(trans, n, nrhs, Ad.data(), n, ipiv, Bd.data(), n);
	if (info != 0) {
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, nrhs, Bd.data(), n, B, ldb);
	return ddd::refine_general_with_double_lu(trans, n, nrhs, Aorig.data(), n, B0.data(), ldb, Ad, ipiv, B, ldb);
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
	ddd::refine_general_with_double_lu('N', n, nrhs, Aorig.data(), n, B0.data(), ldb, Ad, ipiv, B, ldb);
	ddd::convert_double_to_dd_matrix(n, n, Ad.data(), n, A, lda);
	return 0;
}

template <>
inline int tgetri<kv::dd>(const int n, kv::dd* A, const int lda, const int* ipiv) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tgetri<kv::dd>: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	std::vector<kv::dd> Aorig;
	if (!ddd::load_lu_original(A, n, Aorig)) {
		Aorig = ddd::reconstruct_lu_matrix(n, A, lda, ipiv);
	}
	std::vector<kv::dd> X(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	for (int i = 0; i < n; i++) {
		X[i + static_cast<std::size_t>(n) * i] = kv::dd(1.0);
	}
	std::vector<int> ipiv_work(static_cast<std::size_t>(n), 0);
	const int info = tgesv<kv::dd>(n, n, Aorig.data(), n, ipiv_work.data(), X.data(), n);
	if (info != 0) {
		return info;
	}
	ddd::copy_dd_matrix(n, n, X.data(), n, A, lda);
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
inline int tpotrf<kv::dd>(const char uplo, const int n, kv::dd* A, const int lda) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tpotrf<kv::dd>: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n)) {
		detail::tlapack_error("tpotrf<kv::dd>: invalid n/lda");
	}
	if (n == 0) {
		return 0;
	}
	const std::vector<kv::dd> Aorig = ddd::symmetric_full_matrix(uplo, n, A, lda);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, Aorig.data(), n);
	const int info = tpotrf<double>(uplo, n, Ad.data(), n);
	if (info != 0) {
		ddd::convert_double_to_dd_matrix(n, n, Ad.data(), n, A, lda);
		return info;
	}
	ddd::convert_double_to_dd_matrix(n, n, Ad.data(), n, A, lda);
	return ddd::refine_cholesky_factor(uplo, n, Aorig, A, lda);
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
	const std::vector<kv::dd> Aorig = ddd::symmetric_full_matrix(uplo, n, A, lda);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, Aorig.data(), n);
	std::vector<double> wd(static_cast<std::size_t>(n), 0.0);
	const int info = tsyev<double>('V', uplo, n, Ad.data(), n, wd.data());
	if (info != 0) {
		return info;
	}
	for (int i = 0; i < n; i++) {
		w[i] = kv::dd(wd[i]);
	}
	std::vector<kv::dd> X = ddd::to_dd_matrix(n, n, Ad.data(), n);
	const int rinfo = ddd::refine_symmetric_eigenvectors(uplo, n, Aorig, X.data(), n, w);
	if (wantz) {
		ddd::copy_dd_matrix(n, n, X.data(), n, A, lda);
	}
	return rinfo;
}

template <>
inline int tsygv<kv::dd>(const int itype, const char jobz, const char uplo, const int n,
	kv::dd* A, const int lda, kv::dd* B, const int ldb, kv::dd* w) {
	namespace detail = tlapack_detail;
	namespace ddd = tlapack_dd_detail;
	const bool wantz = detail::option_is(jobz, 'V');
	if (itype != 1) {
		detail::tlapack_error("tsygv<kv::dd>: only itype=1 is supported");
	}
	if (!wantz && !detail::option_is(jobz, 'N')) {
		detail::tlapack_error("tsygv<kv::dd>: invalid jobz");
	}
	if (!detail::option_is(uplo, 'U') && !detail::option_is(uplo, 'L')) {
		detail::tlapack_error("tsygv<kv::dd>: invalid uplo");
	}
	if (n < 0 || lda < std::max(1, n) || ldb < std::max(1, n)) {
		detail::tlapack_error("tsygv<kv::dd>: invalid argument");
	}
	if (n == 0) {
		return 0;
	}
	const std::vector<kv::dd> Aorig = ddd::symmetric_full_matrix(uplo, n, A, lda);
	const std::vector<kv::dd> Borig = ddd::symmetric_full_matrix(uplo, n, B, ldb);
	std::vector<double> Ad = ddd::to_double_matrix(n, n, Aorig.data(), n);
	std::vector<double> Bd = ddd::to_double_matrix(n, n, Borig.data(), n);
	std::vector<double> wd(static_cast<std::size_t>(n), 0.0);
	const int info = tsygv<double>(1, 'V', uplo, n, Ad.data(), n, Bd.data(), n, wd.data());
	if (info != 0) {
		ddd::convert_double_to_dd_matrix(n, n, Bd.data(), n, B, ldb);
		return info;
	}
	for (int i = 0; i < n; i++) {
		w[i] = kv::dd(wd[i]);
	}
	std::vector<kv::dd> X = ddd::to_dd_matrix(n, n, Ad.data(), n);
	const int rinfo = ddd::refine_generalized_symmetric_eigenvectors(n, Aorig, Borig, X.data(), n, w);
	if (wantz) {
		ddd::copy_dd_matrix(n, n, X.data(), n, A, lda);
	}
	ddd::copy_dd_matrix(n, n, Borig.data(), n, B, ldb);
	const int cinfo = tpotrf<kv::dd>(uplo, n, B, ldb);
	return rinfo != 0 ? rinfo : cinfo;
}

// Same computation as tsyev<kv::dd> / tsygv<kv::dd>: the double approximation
// (divide and conquer with eigenvectors) followed by the refinement in kv::dd.
template <>
inline int tsyevd<kv::dd>(const char jobz, const char uplo, const int n, kv::dd* A, const int lda, kv::dd* w) {
	return tsyev<kv::dd>(jobz, uplo, n, A, lda, w);
}

template <>
inline int tsygvd<kv::dd>(const int itype, const char jobz, const char uplo, const int n,
	kv::dd* A, const int lda, kv::dd* B, const int ldb, kv::dd* w) {
	return tsygv<kv::dd>(itype, jobz, uplo, n, A, lda, B, ldb, w);
}

} // namespace vcp

#endif // TLAPACK_TLAPACK_DD_HPP
