// VCP Library
// kv::dd specializations for tblas Level 3 routines using the Ozaki scheme.

#pragma once

#ifndef TBLAS_TBLAS_DD_HPP
#define TBLAS_TBLAS_DD_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include <kv/dd.hpp>
#include <vcp/tblas/tblas.hpp>

namespace vcp {
namespace tblas_dd_detail {

inline int ceil_log2(const int k) {
	int r = 0;
	int t = 1;
	while (t < k) {
		t *= 2;
		r++;
	}
	return r;
}

inline int split_shift(const int k) {
	return (51 + ceil_log2(k) + 1) / 2;
}

inline int split_bits(const int k) {
	int b = 52 - split_shift(k);
	if (b < 1) {
		b = 1;
	}
	return b;
}

inline int max_slices(const int k) {
	return 110 / split_bits(k) + 2;
}

inline void pack_op_a(
	const bool transa, const int m, const int k,
	const kv::dd* A, const int lda,
	std::vector< kv::dd >& opA
) {
	opA.assign(static_cast<std::size_t>(m) * k, kv::dd(0.0));
	for (int j = 0; j < k; j++) {
		for (int i = 0; i < m; i++) {
			opA[i + static_cast<std::size_t>(m) * j] =
				transa ? A[j + static_cast<std::size_t>(lda) * i] : A[i + static_cast<std::size_t>(lda) * j];
		}
	}
}

inline void pack_op_b(
	const bool transb, const int k, const int n,
	const kv::dd* B, const int ldb,
	std::vector< kv::dd >& opB
) {
	opB.assign(static_cast<std::size_t>(k) * n, kv::dd(0.0));
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < k; i++) {
			opB[i + static_cast<std::size_t>(k) * j] =
				transb ? B[j + static_cast<std::size_t>(ldb) * i] : B[i + static_cast<std::size_t>(ldb) * j];
		}
	}
}

inline void split_row(
	const std::vector< kv::dd >& A, const int m, const int k,
	std::vector< std::vector< double > >& slices
) {
	using std::fabs;
	using std::frexp;
	using std::ldexp;
	const int s = split_shift(k);
	const int smax = max_slices(k);
	std::vector< kv::dd > rem = A;
	slices.clear();
	for (int p = 0; p < smax; p++) {
		std::vector< double > S(static_cast<std::size_t>(m) * k, 0.0);
		bool nonzero = false;
		for (int i = 0; i < m; i++) {
			double mu = 0.0;
			for (int j = 0; j < k; j++) {
				const double a = fabs(rem[i + static_cast<std::size_t>(m) * j].a1);
				if (a > mu) {
					mu = a;
				}
			}
			if (mu == 0.0) {
				continue;
			}
			int e;
			frexp(mu, &e);
			const double sigma = ldexp(1.0, e + s);
			for (int j = 0; j < k; j++) {
				const std::size_t idx = i + static_cast<std::size_t>(m) * j;
				const double q = (rem[idx].a1 + sigma) - sigma;
				if (q != 0.0) {
					S[idx] = q;
					rem[idx] -= kv::dd(q);
					nonzero = true;
				}
			}
		}
		if (!nonzero) {
			break;
		}
		slices.push_back(S);
	}
}

inline void split_col(
	const std::vector< kv::dd >& A, const int k, const int n,
	std::vector< std::vector< double > >& slices
) {
	using std::fabs;
	using std::frexp;
	using std::ldexp;
	const int s = split_shift(k);
	const int smax = max_slices(k);
	std::vector< kv::dd > rem = A;
	slices.clear();
	for (int p = 0; p < smax; p++) {
		std::vector< double > S(static_cast<std::size_t>(k) * n, 0.0);
		bool nonzero = false;
		for (int j = 0; j < n; j++) {
			double mu = 0.0;
			for (int i = 0; i < k; i++) {
				const double a = fabs(rem[i + static_cast<std::size_t>(k) * j].a1);
				if (a > mu) {
					mu = a;
				}
			}
			if (mu == 0.0) {
				continue;
			}
			int e;
			frexp(mu, &e);
			const double sigma = ldexp(1.0, e + s);
			for (int i = 0; i < k; i++) {
				const std::size_t idx = i + static_cast<std::size_t>(k) * j;
				const double q = (rem[idx].a1 + sigma) - sigma;
				if (q != 0.0) {
					S[idx] = q;
					rem[idx] -= kv::dd(q);
					nonzero = true;
				}
			}
		}
		if (!nonzero) {
			break;
		}
		slices.push_back(S);
	}
}

inline void ozaki_product(
	const std::vector< kv::dd >& opA, const std::vector< kv::dd >& opB,
	const int m, const int n, const int k,
	std::vector< kv::dd >& C
) {
	std::vector< std::vector< double > > DA;
	std::vector< std::vector< double > > DB;
	split_row(opA, m, k, DA);
	split_col(opB, k, n, DB);

	C.assign(static_cast<std::size_t>(m) * n, kv::dd(0.0));
	if (DA.empty() || DB.empty()) {
		return;
	}

	const double one = 1.0;
	const double zero = 0.0;
	const int limit = max_slices(k);
	std::vector< double > P(static_cast<std::size_t>(m) * n, 0.0);
	for (int t = 0; t < limit; t++) {
		for (int p = 0; p <= t; p++) {
			const int q = t - p;
			if (p >= static_cast<int>(DA.size()) || q >= static_cast<int>(DB.size())) {
				continue;
			}
			tgemm<double>('N', 'N', m, n, k, one, DA[p].data(), m, DB[q].data(), k, zero, P.data(), m);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (tblas_detail::use_parallel(static_cast<double>(m) * n))
#endif
			for (int idx = 0; idx < m * n; idx++) {
				if (P[static_cast<std::size_t>(idx)] != 0.0) {
					C[static_cast<std::size_t>(idx)] += kv::dd(P[static_cast<std::size_t>(idx)]);
				}
			}
		}
	}
}

inline void merge_full_to_triangle(
	const bool upper, const int n,
	const kv::dd& beta, const std::vector< kv::dd >& T,
	kv::dd* C, const int ldc
) {
	const bool zero_beta = (beta == kv::dd(0.0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (tblas_detail::use_parallel(0.5 * n * n))
#endif
	for (int j = 0; j < n; j++) {
		const int i0 = upper ? 0 : j;
		const int i1 = upper ? j + 1 : n;
		for (int i = i0; i < i1; i++) {
			const std::size_t tidx = i + static_cast<std::size_t>(n) * j;
			const std::size_t cidx = i + static_cast<std::size_t>(ldc) * j;
			C[cidx] = zero_beta ? T[tidx] : beta * C[cidx] + T[tidx];
		}
	}
}

} // namespace tblas_dd_detail

template <>
inline void tgemm<kv::dd>(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const kv::dd& alpha, const kv::dd* A, const int lda,
	const kv::dd* B, const int ldb,
	const kv::dd& beta, kv::dd* C, const int ldc
) {
	namespace det = tblas_detail;
	namespace ddd = tblas_dd_detail;
	const bool ta = !det::option_is(transa, 'N');
	const bool tb = !det::option_is(transb, 'N');
	if (ta && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("tgemm<kv::dd>: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::tblas_error("tgemm<kv::dd>: invalid transb");
	}
	const int nrowa = ta ? k : m;
	const int nrowb = tb ? n : k;
	if (m < 0 || n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, m)) {
		det::tblas_error("tgemm<kv::dd>: invalid argument");
	}
	if (m == 0 || n == 0 || (k == 0 && beta == kv::dd(1.0)) || (alpha == kv::dd(0.0) && beta == kv::dd(1.0))) {
		return;
	}
	if (alpha == kv::dd(0.0) || k == 0) {
		det::scale_matrix(m, n, beta, C, ldc);
		return;
	}

	std::vector< kv::dd > opA;
	std::vector< kv::dd > opB;
	std::vector< kv::dd > prod;
	ddd::pack_op_a(ta, m, k, A, lda, opA);
	ddd::pack_op_b(tb, k, n, B, ldb, opB);
	ddd::ozaki_product(opA, opB, m, n, k, prod);

	const bool zero_beta = (beta == kv::dd(0.0));
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(static_cast<double>(m) * n))
#endif
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			const std::size_t pidx = i + static_cast<std::size_t>(m) * j;
			const std::size_t cidx = i + static_cast<std::size_t>(ldc) * j;
			C[cidx] = zero_beta ? alpha * prod[pidx] : beta * C[cidx] + alpha * prod[pidx];
		}
	}
}

template <>
inline void tsymm<kv::dd>(
	const char side, const char uplo,
	const int m, const int n,
	const kv::dd& alpha, const kv::dd* A, const int lda,
	const kv::dd* B, const int ldb,
	const kv::dd& beta, kv::dd* C, const int ldc
) {
	namespace det = tblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	if (!lside && !det::option_is(side, 'R')) {
		det::tblas_error("tsymm<kv::dd>: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsymm<kv::dd>: invalid uplo");
	}
	const int ka = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, ka) || ldb < std::max(1, m) || ldc < std::max(1, m)) {
		det::tblas_error("tsymm<kv::dd>: invalid argument");
	}
	if (m == 0 || n == 0 || (alpha == kv::dd(0.0) && beta == kv::dd(1.0))) {
		return;
	}
	if (alpha == kv::dd(0.0)) {
		det::scale_matrix(m, n, beta, C, ldc);
		return;
	}

	std::vector< kv::dd > Afull(static_cast<std::size_t>(ka) * ka, kv::dd(0.0));
	det::expand_symmetric(upper, ka, A, lda, Afull.data());
	if (lside) {
		tgemm<kv::dd>('N', 'N', m, n, m, alpha, Afull.data(), m, B, ldb, beta, C, ldc);
	}
	else {
		tgemm<kv::dd>('N', 'N', m, n, n, alpha, B, ldb, Afull.data(), n, beta, C, ldc);
	}
}

template <>
inline void tsyrk<kv::dd>(
	const char uplo, const char trans,
	const int n, const int k,
	const kv::dd& alpha, const kv::dd* A, const int lda,
	const kv::dd& beta, kv::dd* C, const int ldc
) {
	namespace det = tblas_detail;
	namespace ddd = tblas_dd_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsyrk<kv::dd>: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("tsyrk<kv::dd>: invalid trans");
	}
	const int nrowa = ntrans ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldc < std::max(1, n)) {
		det::tblas_error("tsyrk<kv::dd>: invalid argument");
	}
	if (n == 0 || (alpha == kv::dd(0.0) && beta == kv::dd(1.0)) || (k == 0 && beta == kv::dd(1.0))) {
		return;
	}
	if (alpha == kv::dd(0.0) || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc);
		return;
	}

	std::vector< kv::dd > T(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	if (ntrans) {
		tgemm<kv::dd>('N', 'T', n, n, k, alpha, A, lda, A, lda, kv::dd(0.0), T.data(), n);
	}
	else {
		tgemm<kv::dd>('T', 'N', n, n, k, alpha, A, lda, A, lda, kv::dd(0.0), T.data(), n);
	}
	ddd::merge_full_to_triangle(upper, n, beta, T, C, ldc);
}

template <>
inline void tsyr2k<kv::dd>(
	const char uplo, const char trans,
	const int n, const int k,
	const kv::dd& alpha, const kv::dd* A, const int lda,
	const kv::dd* B, const int ldb,
	const kv::dd& beta, kv::dd* C, const int ldc
) {
	namespace det = tblas_detail;
	namespace ddd = tblas_dd_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(trans, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tsyr2k<kv::dd>: invalid uplo");
	}
	if (!ntrans && !det::option_is(trans, 'T') && !det::option_is(trans, 'C')) {
		det::tblas_error("tsyr2k<kv::dd>: invalid trans");
	}
	const int nrowa = ntrans ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowa) || ldc < std::max(1, n)) {
		det::tblas_error("tsyr2k<kv::dd>: invalid argument");
	}
	if (n == 0 || (alpha == kv::dd(0.0) && beta == kv::dd(1.0)) || (k == 0 && beta == kv::dd(1.0))) {
		return;
	}
	if (alpha == kv::dd(0.0) || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc);
		return;
	}

	std::vector< kv::dd > T(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	if (ntrans) {
		tgemm<kv::dd>('N', 'T', n, n, k, alpha, A, lda, B, ldb, kv::dd(0.0), T.data(), n);
		tgemm<kv::dd>('N', 'T', n, n, k, alpha, B, ldb, A, lda, kv::dd(1.0), T.data(), n);
	}
	else {
		tgemm<kv::dd>('T', 'N', n, n, k, alpha, A, lda, B, ldb, kv::dd(0.0), T.data(), n);
		tgemm<kv::dd>('T', 'N', n, n, k, alpha, B, ldb, A, lda, kv::dd(1.0), T.data(), n);
	}
	ddd::merge_full_to_triangle(upper, n, beta, T, C, ldc);
}

template <>
inline void ttrmm<kv::dd>(
	const char side, const char uplo, const char transa, const char diag,
	const int m, const int n,
	const kv::dd& alpha, const kv::dd* A, const int lda,
	kv::dd* B, const int ldb
) {
	namespace det = tblas_detail;
	const bool lside = det::option_is(side, 'L');
	const bool upper = det::option_is(uplo, 'U');
	const bool ntrans = det::option_is(transa, 'N');
	const bool nounit = det::option_is(diag, 'N');
	if (!lside && !det::option_is(side, 'R')) {
		det::tblas_error("ttrmm<kv::dd>: invalid side");
	}
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("ttrmm<kv::dd>: invalid uplo");
	}
	if (!ntrans && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("ttrmm<kv::dd>: invalid transa");
	}
	if (!nounit && !det::option_is(diag, 'U')) {
		det::tblas_error("ttrmm<kv::dd>: invalid diag");
	}
	const int nrowa = lside ? m : n;
	if (m < 0 || n < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, m)) {
		det::tblas_error("ttrmm<kv::dd>: invalid argument");
	}
	if (m == 0 || n == 0) {
		return;
	}
	if (alpha == kv::dd(0.0)) {
		det::scale_matrix(m, n, kv::dd(0.0), B, ldb);
		return;
	}

	const int nb = lside ? m : n;
	std::vector< kv::dd > Adense(static_cast<std::size_t>(nb) * nb, kv::dd(0.0));
	std::vector< kv::dd > temp(static_cast<std::size_t>(m) * n, kv::dd(0.0));
	det::expand_triangular(upper, nounit, nb, A, lda, Adense.data());
	if (lside) {
		tgemm<kv::dd>(ntrans ? 'N' : 'T', 'N', m, n, m, alpha, Adense.data(), m, B, ldb, kv::dd(0.0), temp.data(), m);
	}
	else {
		tgemm<kv::dd>('N', ntrans ? 'N' : 'T', m, n, n, alpha, B, ldb, Adense.data(), n, kv::dd(0.0), temp.data(), m);
	}
#ifdef _OPENMP
#pragma omp parallel for schedule(static) if (det::use_parallel(static_cast<double>(m) * n))
#endif
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			B[i + static_cast<std::size_t>(ldb) * j] = temp[i + static_cast<std::size_t>(m) * j];
		}
	}
}

template <>
inline void tgemmtr<kv::dd>(
	const char uplo, const char transa, const char transb,
	const int n, const int k,
	const kv::dd& alpha, const kv::dd* A, const int lda,
	const kv::dd* B, const int ldb,
	const kv::dd& beta, kv::dd* C, const int ldc
) {
	namespace det = tblas_detail;
	namespace ddd = tblas_dd_detail;
	const bool upper = det::option_is(uplo, 'U');
	const bool ta = !det::option_is(transa, 'N');
	const bool tb = !det::option_is(transb, 'N');
	if (!upper && !det::option_is(uplo, 'L')) {
		det::tblas_error("tgemmtr<kv::dd>: invalid uplo");
	}
	if (ta && !det::option_is(transa, 'T') && !det::option_is(transa, 'C')) {
		det::tblas_error("tgemmtr<kv::dd>: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::tblas_error("tgemmtr<kv::dd>: invalid transb");
	}
	const int nrowa = ta ? k : n;
	const int nrowb = tb ? n : k;
	if (n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, n)) {
		det::tblas_error("tgemmtr<kv::dd>: invalid argument");
	}
	if (n == 0 || ((alpha == kv::dd(0.0) || k == 0) && beta == kv::dd(1.0))) {
		return;
	}
	if (alpha == kv::dd(0.0) || k == 0) {
		det::scale_triangle(upper, n, beta, C, ldc);
		return;
	}

	std::vector< kv::dd > T(static_cast<std::size_t>(n) * n, kv::dd(0.0));
	tgemm<kv::dd>(transa, transb, n, n, k, alpha, A, lda, B, ldb, kv::dd(0.0), T.data(), n);
	ddd::merge_full_to_triangle(upper, n, beta, T, C, ldc);
}

} // namespace vcp

#endif // TBLAS_TBLAS_DD_HPP
