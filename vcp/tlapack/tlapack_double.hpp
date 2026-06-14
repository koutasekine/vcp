// VCP Library
// double specializations for tlapack templates.

#pragma once

#ifndef TLAPACK_TLAPACK_DOUBLE_HPP
#define TLAPACK_TLAPACK_DOUBLE_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <vcp/dblas_dlapack.hpp>
#include <vcp/tblas/tblas_double.hpp>
#include <vcp/tlapack/tlapack.hpp>

#ifndef USE_VCP_LAPACK
extern "C" {
	void dgetrf2_(const int*, const int*, double*, const int*, int*, int*);
	void dpotrf2_(const char*, const int*, double*, const int*, int*);
	void dorg2r_(const int*, const int*, const int*, double*, const int*, const double*, double*, int*);
	void dorm2r_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, int*);
	void dorgl2_(const int*, const int*, const int*, double*, const int*, const double*, double*, int*);
	void dorml2_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, int*);
	void dorg2l_(const int*, const int*, const int*, double*, const int*, const double*, double*, int*);
	void dorm2l_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, int*);
	void dlae2_(const double*, const double*, const double*, double*, double*);
	void dlaev2_(const double*, const double*, const double*, double*, double*, double*, double*);
	void dlas2_(const double*, const double*, const double*, double*, double*);
	void dlasv2_(const double*, const double*, const double*, double*, double*, double*, double*, double*, double*);
	void dlanv2_(double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);
	void dlasr_(const char*, const char*, const char*, const int*, const int*,
		const double*, const double*, double*, const int*);
	void dlarf_(const char*, const int*, const int*, const double*, const int*,
		const double*, double*, const int*, double*);
	void dlarft_(const char*, const char*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*);
	void dlarfb_(const char*, const char*, const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, const int*, double*, const int*, double*, const int*);
	void dlasyf_(const char*, const int*, const int*, int*, double*, const int*, int*, double*, const int*, int*);
	void dlatrd_(const char*, const int*, const int*, double*, const int*, double*, double*, double*, const int*);
	void dlahr2_(const int*, const int*, const int*, double*, const int*, double*, double*, const int*, double*, const int*);
	void dlabrd_(const int*, const int*, const int*, double*, const int*, double*, double*, double*, double*,
		double*, const int*, double*, const int*);
	void dlahqr_(const int*, const int*, const int*, const int*, const int*, double*, const int*,
		double*, double*, const int*, const int*, double*, const int*, int*);
	void dlaln2_(const int*, const int*, const int*, const double*, const double*, const double*, const int*,
		const double*, const double*, const double*, const int*, const double*, const double*, double*, const int*,
		double*, double*, int*);
}
#endif

namespace tlapack_double_detail {

inline int work_size(const double query, const int fallback) {
	if (query > static_cast<double>(fallback)) {
		return static_cast<int>(query);
	}
	return fallback;
}

inline void check_info(const int info, const char* name) {
	if (info < 0) {
		tlapack_detail::tlapack_error(name);
	}
}

inline std::vector<int> piv_to_fortran(const int n, const int* ipiv) {
	std::vector<int> piv(static_cast<std::size_t>(std::max(0, n)));
	for (int i = 0; i < n; i++) {
		piv[static_cast<std::size_t>(i)] = ipiv[i] + 1;
	}
	return piv;
}

inline void piv_to_c(const int n, int* ipiv) {
	for (int i = 0; i < n; i++) {
		ipiv[i] -= 1;
	}
}

inline int minmn(const int m, const int n) {
	return m < n ? m : n;
}

} // namespace tlapack_double_detail

namespace vcp {

template <>
inline double tlamch<double>(const char cmach) {
	return dlamch_(&cmach);
}

template <>
inline void tlaset<double>(const char uplo, const int m, const int n, const double& alpha, const double& beta, double* A, const int lda) {
	dlaset_(&uplo, &m, &n, &alpha, &beta, A, &lda);
}

template <>
inline void tlacpy<double>(const char uplo, const int m, const int n, const double* A, const int lda, double* B, const int ldb) {
	dlacpy_(&uplo, &m, &n, A, &lda, B, &ldb);
}

template <>
inline void tlaswp<double>(const int n, double* A, const int lda, const int k1, const int k2, const int* ipiv, const int incx) {
	const std::vector<int> piv = tlapack_double_detail::piv_to_fortran(k2, ipiv);
	const int fk1 = k1 + 1;
	dlaswp_(&n, A, &lda, &fk1, &k2, piv.data(), &incx);
}

template <>
inline void tlasrt<double>(const char id, const int n, double* d) {
	int info = 0;
	dlasrt_(&id, &n, d, &info);
	tlapack_double_detail::check_info(info, "tlasrt<double>");
}

template <>
inline void tlassq<double>(const int n, const double* x, const int incx, double& scale, double& sumsq) {
	dlassq_(&n, x, &incx, &scale, &sumsq);
}

template <>
inline double tlapy2<double>(const double& x, const double& y) {
	return dlapy2_(&x, &y);
}

template <>
inline double tlapy3<double>(const double& x, const double& y, const double& z) {
	return dlapy3_(&x, &y, &z);
}

template <>
inline void tlartg<double>(const double& f, const double& g, double& cs, double& sn, double& r) {
	dlartg_(&f, &g, &cs, &sn, &r);
}

template <>
inline void tlascl<double>(const char type, const int kl, const int ku, const double& cfrom, const double& cto,
	const int m, const int n, double* A, const int lda) {
	int info = 0;
	dlascl_(&type, &kl, &ku, &cfrom, &cto, &m, &n, A, &lda, &info);
	tlapack_double_detail::check_info(info, "tlascl<double>");
}

template <>
inline double tlange<double>(const char norm, const int m, const int n, const double* A, const int lda) {
	std::vector<double> work(static_cast<std::size_t>(std::max(1, m)));
	return dlange_(&norm, &m, &n, A, &lda, work.data());
}

template <>
inline double tlansy<double>(const char norm, const char uplo, const int n, const double* A, const int lda) {
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	return dlansy_(&norm, &uplo, &n, A, &lda, work.data());
}

template <>
inline double tlanst<double>(const char norm, const int n, const double* d, const double* e) {
	return dlanst_(&norm, &n, d, e);
}

template <>
inline double tlangb<double>(const char norm, const int n, const int kl, const int ku, const double* AB, const int ldab) {
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	return dlangb_(&norm, &n, &kl, &ku, AB, &ldab, work.data());
}

template <>
inline double tlansb<double>(const char norm, const char uplo, const int n, const int k, const double* AB, const int ldab) {
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	return dlansb_(&norm, &uplo, &n, &k, AB, &ldab, work.data());
}

template <>
inline void tlarfg<double>(const int n, double& alpha, double* x, const int incx, double& tau) {
	dlarfg_(&n, &alpha, x, &incx, &tau);
}

template <>
inline int tgetf2<double>(const int m, const int n, double* A, const int lda, int* ipiv) {
	int info = 0;
	dgetf2_(&m, &n, A, &lda, ipiv, &info);
	tlapack_double_detail::piv_to_c(tlapack_double_detail::minmn(m, n), ipiv);
	return info;
}

template <>
inline int tgetrf<double>(const int m, const int n, double* A, const int lda, int* ipiv) {
	int info = 0;
	dgetrf_(&m, &n, A, &lda, ipiv, &info);
	tlapack_double_detail::piv_to_c(tlapack_double_detail::minmn(m, n), ipiv);
	return info;
}

template <>
inline int tgetrs<double>(const char trans, const int n, const int nrhs, const double* A, const int lda,
	const int* ipiv, double* B, const int ldb) {
	std::vector<int> piv = tlapack_double_detail::piv_to_fortran(n, ipiv);
	int info = 0;
	dgetrs_(&trans, &n, &nrhs, A, &lda, piv.data(), B, &ldb, &info);
	return info;
}

template <>
inline int tgesv<double>(const int n, const int nrhs, double* A, const int lda, int* ipiv, double* B, const int ldb) {
	int info = 0;
	dgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
	tlapack_double_detail::piv_to_c(n, ipiv);
	return info;
}

template <>
inline int tgetri<double>(const int n, double* A, const int lda, const int* ipiv) {
	std::vector<int> piv = tlapack_double_detail::piv_to_fortran(n, ipiv);
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgetri_(&n, A, &lda, piv.data(), &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgetri_(&n, A, &lda, piv.data(), work.data(), &lwork, &info);
	return info;
}

template <>
inline int ttrti2<double>(const char uplo, const char diag, const int n, double* A, const int lda) {
	int info = 0;
	dtrti2_(&uplo, &diag, &n, A, &lda, &info);
	return info;
}

template <>
inline int ttrtri<double>(const char uplo, const char diag, const int n, double* A, const int lda) {
	int info = 0;
	dtrtri_(&uplo, &diag, &n, A, &lda, &info);
	return info;
}

template <>
inline int ttrtrs<double>(const char uplo, const char trans, const char diag, const int n, const int nrhs,
	const double* A, const int lda, double* B, const int ldb) {
	int info = 0;
	dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, A, &lda, B, &ldb, &info);
	return info;
}

template <>
inline int tlauu2<double>(const char uplo, const int n, double* A, const int lda) {
	int info = 0;
	dlauu2_(&uplo, &n, A, &lda, &info);
	return info;
}

template <>
inline int tlauum<double>(const char uplo, const int n, double* A, const int lda) {
	int info = 0;
	dlauum_(&uplo, &n, A, &lda, &info);
	return info;
}

template <>
inline int tpotf2<double>(const char uplo, const int n, double* A, const int lda) {
	int info = 0;
	dpotf2_(&uplo, &n, A, &lda, &info);
	return info;
}

template <>
inline int tpotrf<double>(const char uplo, const int n, double* A, const int lda) {
	int info = 0;
	dpotrf_(&uplo, &n, A, &lda, &info);
	return info;
}

template <>
inline int tpotrs<double>(const char uplo, const int n, const int nrhs, const double* A, const int lda, double* B, const int ldb) {
	int info = 0;
	dpotrs_(&uplo, &n, &nrhs, A, &lda, B, &ldb, &info);
	return info;
}

template <>
inline int tposv<double>(const char uplo, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb) {
	int info = 0;
	dposv_(&uplo, &n, &nrhs, A, &lda, B, &ldb, &info);
	return info;
}

template <>
inline int tpotri<double>(const char uplo, const int n, double* A, const int lda) {
	int info = 0;
	dpotri_(&uplo, &n, A, &lda, &info);
	return info;
}

template <>
inline int tsytf2<double>(const char uplo, const int n, double* A, const int lda, int* ipiv) {
	int info = 0;
	dsytf2_(&uplo, &n, A, &lda, ipiv, &info);
	return info;
}

template <>
inline int tsytrf<double>(const char uplo, const int n, double* A, const int lda, int* ipiv) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dsytrf_(&uplo, &n, A, &lda, ipiv, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dsytrf_(&uplo, &n, A, &lda, ipiv, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tsytrs<double>(const char uplo, const int n, const int nrhs, const double* A, const int lda,
	const int* ipiv, double* B, const int ldb) {
	int info = 0;
	dsytrs_(&uplo, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
	return info;
}

template <>
inline int tsysv<double>(const char uplo, const int n, const int nrhs, double* A, const int lda,
	int* ipiv, double* B, const int ldb) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dsysv_(&uplo, &n, &nrhs, A, &lda, ipiv, B, &ldb, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dsysv_(&uplo, &n, &nrhs, A, &lda, ipiv, B, &ldb, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tsygs2<double>(const int itype, const char uplo, const int n, double* A, const int lda,
	const double* B, const int ldb) {
	int info = 0;
	dsygs2_(&itype, &uplo, &n, A, &lda, B, &ldb, &info);
	return info;
}

template <>
inline int tsygst<double>(const int itype, const char uplo, const int n, double* A, const int lda,
	const double* B, const int ldb) {
	int info = 0;
	dsygst_(&itype, &uplo, &n, A, &lda, B, &ldb, &info);
	return info;
}

template <>
inline int tsygv<double>(const int itype, const char jobz, const char uplo, const int n,
	double* A, const int lda, double* B, const int ldb, double* w) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dsygv_(&itype, &jobz, &uplo, &n, A, &lda, B, &ldb, w, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, 3 * n - 1));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dsygv_(&itype, &jobz, &uplo, &n, A, &lda, B, &ldb, w, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tgbtf2<double>(const int m, const int n, const int kl, const int ku,
	double* AB, const int ldab, int* ipiv) {
	int info = 0;
	dgbtf2_(&m, &n, &kl, &ku, AB, &ldab, ipiv, &info);
	tlapack_double_detail::piv_to_c(tlapack_double_detail::minmn(m, n), ipiv);
	return info;
}

template <>
inline int tgbtrf<double>(const int m, const int n, const int kl, const int ku,
	double* AB, const int ldab, int* ipiv) {
	int info = 0;
	dgbtrf_(&m, &n, &kl, &ku, AB, &ldab, ipiv, &info);
	tlapack_double_detail::piv_to_c(tlapack_double_detail::minmn(m, n), ipiv);
	return info;
}

template <>
inline int tgbtrs<double>(const char trans, const int n, const int kl, const int ku, const int nrhs,
	const double* AB, const int ldab, const int* ipiv, double* B, const int ldb) {
	std::vector<int> piv = tlapack_double_detail::piv_to_fortran(n, ipiv);
	int info = 0;
	dgbtrs_(&trans, &n, &kl, &ku, &nrhs, AB, &ldab, piv.data(), B, &ldb, &info);
	return info;
}

template <>
inline int tgbsv<double>(const int n, const int kl, const int ku, const int nrhs,
	double* AB, const int ldab, int* ipiv, double* B, const int ldb) {
	int info = 0;
	dgbsv_(&n, &kl, &ku, &nrhs, AB, &ldab, ipiv, B, &ldb, &info);
	tlapack_double_detail::piv_to_c(n, ipiv);
	return info;
}

template <>
inline int tpbtf2<double>(const char uplo, const int n, const int kd, double* AB, const int ldab) {
	int info = 0;
	dpbtf2_(&uplo, &n, &kd, AB, &ldab, &info);
	return info;
}

template <>
inline int tpbtrf<double>(const char uplo, const int n, const int kd, double* AB, const int ldab) {
	int info = 0;
	dpbtrf_(&uplo, &n, &kd, AB, &ldab, &info);
	return info;
}

template <>
inline int tpbtrs<double>(const char uplo, const int n, const int kd, const int nrhs,
	const double* AB, const int ldab, double* B, const int ldb) {
	int info = 0;
	dpbtrs_(&uplo, &n, &kd, &nrhs, AB, &ldab, B, &ldb, &info);
	return info;
}

template <>
inline int tpbsv<double>(const char uplo, const int n, const int kd, const int nrhs,
	double* AB, const int ldab, double* B, const int ldb) {
	int info = 0;
	dpbsv_(&uplo, &n, &kd, &nrhs, AB, &ldab, B, &ldb, &info);
	return info;
}

template <>
inline int tgeqr2<double>(const int m, const int n, double* A, const int lda, double* tau) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	dgeqr2_(&m, &n, A, &lda, tau, work.data(), &info);
	return info;
}

template <>
inline int tgeqrf<double>(const int m, const int n, double* A, const int lda, double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgeqrf_(&m, &n, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgeqrf_(&m, &n, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tgelq2<double>(const int m, const int n, double* A, const int lda, double* tau) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, m)));
	dgelq2_(&m, &n, A, &lda, tau, work.data(), &info);
	return info;
}

template <>
inline int tgelqf<double>(const int m, const int n, double* A, const int lda, double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgelqf_(&m, &n, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, m));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgelqf_(&m, &n, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int torgqr<double>(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dorgqr_(&m, &n, &k, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorgqr_(&m, &n, &k, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int torglq<double>(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dorglq_(&m, &n, &k, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, m));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorglq_(&m, &n, &k, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int torgql<double>(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dorgql_(&m, &n, &k, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorgql_(&m, &n, &k, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tormqr<double>(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dormqr_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dormqr_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tormlq<double>(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dormlq_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, m));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dormlq_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tormql<double>(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dormql_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dormql_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tgels<double>(const char trans, const int m, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgels_(&trans, &m, &n, &nrhs, A, &lda, B, &ldb, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, tlapack_double_detail::minmn(m, n) + std::max(tlapack_double_detail::minmn(m, n), nrhs)));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgels_(&trans, &m, &n, &nrhs, A, &lda, B, &ldb, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tsytd2<double>(const char uplo, const int n, double* A, const int lda, double* d, double* e, double* tau) {
	int info = 0;
	dsytd2_(&uplo, &n, A, &lda, d, e, tau, &info);
	return info;
}

template <>
inline int tsytrd<double>(const char uplo, const int n, double* A, const int lda, double* d, double* e, double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dsytrd_(&uplo, &n, A, &lda, d, e, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dsytrd_(&uplo, &n, A, &lda, d, e, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int torgtr<double>(const char uplo, const int n, double* A, const int lda, const double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dorgtr_(&uplo, &n, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorgtr_(&uplo, &n, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tormtr<double>(const char side, const char uplo, const char trans, const int m, const int n,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dormtr_(&side, &uplo, &trans, &m, &n, A, &lda, tau, C, &ldc, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dormtr_(&side, &uplo, &trans, &m, &n, A, &lda, tau, C, &ldc, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tsterf<double>(const int n, double* d, double* e) {
	int info = 0;
	dsterf_(&n, d, e, &info);
	return info;
}

template <>
inline int tsteqr<double>(const char compz, const int n, double* d, double* e, double* Z, const int ldz) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, 2 * n - 2)));
	dsteqr_(&compz, &n, d, e, Z, &ldz, work.data(), &info);
	return info;
}

template <>
inline int tsyev<double>(const char jobz, const char uplo, const int n, double* A, const int lda, double* w) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dsyev_(&jobz, &uplo, &n, A, &lda, w, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, 3 * n - 1));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dsyev_(&jobz, &uplo, &n, A, &lda, w, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tgebd2<double>(const int m, const int n, double* A, const int lda,
	double* d, double* e, double* tauq, double* taup) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, std::max(m, n))));
	dgebd2_(&m, &n, A, &lda, d, e, tauq, taup, work.data(), &info);
	return info;
}

template <>
inline int tgebrd<double>(const int m, const int n, double* A, const int lda,
	double* d, double* e, double* tauq, double* taup) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgebrd_(&m, &n, A, &lda, d, e, tauq, taup, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, std::max(m, n)));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgebrd_(&m, &n, A, &lda, d, e, tauq, taup, work.data(), &lwork, &info);
	return info;
}

template <>
inline int torgbr<double>(const char vect, const int m, const int n, const int k,
	double* A, const int lda, const double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dorgbr_(&vect, &m, &n, &k, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, m));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorgbr_(&vect, &m, &n, &k, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tormbr<double>(const char vect, const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dormbr_(&vect, &side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dormbr_(&vect, &side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tbdsqr<double>(const char uplo, const int n, const int ncvt, const int nru, const int ncc,
	double* d, double* e, double* VT, const int ldvt, double* U, const int ldu, double* C, const int ldc) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, 4 * n)));
	dbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, VT, &ldvt, U, &ldu, C, &ldc, work.data(), &info);
	return info;
}

template <>
inline int tgesvd<double>(const char jobu, const char jobvt, const int m, const int n,
	double* A, const int lda, double* s, double* U, const int ldu, double* VT, const int ldvt) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VT, &ldvt, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, 5 * tlapack_double_detail::minmn(m, n)));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgesvd_(&jobu, &jobvt, &m, &n, A, &lda, s, U, &ldu, VT, &ldvt, work.data(), &lwork, &info);
	return info;
}

template <>
inline int tgebal<double>(const char job, const int n, double* A, const int lda, int& ilo, int& ihi, double* scale) {
	int info = 0;
	dgebal_(&job, &n, A, &lda, &ilo, &ihi, scale, &info);
	return info;
}

template <>
inline int tgebak<double>(const char job, const char side, const int n, const int ilo, const int ihi,
	const double* scale, const int m, double* V, const int ldv) {
	int info = 0;
	dgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, V, &ldv, &info);
	return info;
}

template <>
inline int tgehd2<double>(const int n, const int ilo, const int ihi, double* A, const int lda, double* tau) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	dgehd2_(&n, &ilo, &ihi, A, &lda, tau, work.data(), &info);
	return info;
}

template <>
inline int tgehrd<double>(const int n, const int ilo, const int ihi, double* A, const int lda, double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgehrd_(&n, &ilo, &ihi, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgehrd_(&n, &ilo, &ihi, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int torghr<double>(const int n, const int ilo, const int ihi, double* A, const int lda, const double* tau) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dorghr_(&n, &ilo, &ihi, A, &lda, tau, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorghr_(&n, &ilo, &ihi, A, &lda, tau, work.data(), &lwork, &info);
	return info;
}

template <>
inline int thseqr<double>(const char job, const char compz, const int n, const int ilo, const int ihi,
	double* H, const int ldh, double* wr, double* wi, double* Z, const int ldz) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dhseqr_(&job, &compz, &n, &ilo, &ihi, H, &ldh, wr, wi, Z, &ldz, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dhseqr_(&job, &compz, &n, &ilo, &ihi, H, &ldh, wr, wi, Z, &ldz, work.data(), &lwork, &info);
	return info;
}

template <>
inline int ttrevc<double>(const char side, const char howmny, const int n, const double* Tf, const int ldt,
	double* VL, const int ldvl, double* VR, const int ldvr, const int mm, int& m) {
	int info = 0;
	std::vector<int> select(static_cast<std::size_t>(std::max(1, n)), 0);
	std::vector<double> work(static_cast<std::size_t>(std::max(1, 3 * n)));
	dtrevc_(&side, &howmny, select.data(), &n, Tf, &ldt, VL, &ldvl, VR, &ldvr, &mm, &m, work.data(), &info);
	return info;
}

template <>
inline int tgeev<double>(const char jobvl, const char jobvr, const int n, double* A, const int lda,
	double* wr, double* wi, double* VL, const int ldvl, double* VR, const int ldvr) {
	int info = 0;
	int lwork = -1;
	double query = 0.0;
	dgeev_(&jobvl, &jobvr, &n, A, &lda, wr, wi, VL, &ldvl, VR, &ldvr, &query, &lwork, &info);
	if (info != 0) {
		return info;
	}
	lwork = tlapack_double_detail::work_size(query, std::max(1, 4 * n));
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dgeev_(&jobvl, &jobvr, &n, A, &lda, wr, wi, VL, &ldvl, VR, &ldvr, work.data(), &lwork, &info);
	return info;
}

#ifndef USE_VCP_LAPACK
template <>
inline int tgetrf2<double>(const int m, const int n, double* A, const int lda, int* ipiv) {
	int info = 0;
	dgetrf2_(&m, &n, A, &lda, ipiv, &info);
	tlapack_double_detail::piv_to_c(tlapack_double_detail::minmn(m, n), ipiv);
	return info;
}

template <>
inline int tpotrf2<double>(const char uplo, const int n, double* A, const int lda) {
	int info = 0;
	dpotrf2_(&uplo, &n, A, &lda, &info);
	return info;
}

template <>
inline void tlae2<double>(const double& a, const double& b, const double& c, double& rt1, double& rt2) {
	dlae2_(&a, &b, &c, &rt1, &rt2);
}

template <>
inline void tlaev2<double>(const double& a, const double& b, const double& c, double& rt1, double& rt2, double& cs1, double& sn1) {
	dlaev2_(&a, &b, &c, &rt1, &rt2, &cs1, &sn1);
}

template <>
inline void tlas2<double>(const double& f, const double& g, const double& h, double& ssmin, double& ssmax) {
	dlas2_(&f, &g, &h, &ssmin, &ssmax);
}

template <>
inline void tlasv2<double>(const double& f, const double& g, const double& h,
	double& ssmin, double& ssmax, double& snr, double& csr, double& snl, double& csl) {
	dlasv2_(&f, &g, &h, &ssmin, &ssmax, &snr, &csr, &snl, &csl);
}

template <>
inline void tlanv2<double>(double& a, double& b, double& c, double& d,
	double& rt1r, double& rt1i, double& rt2r, double& rt2i, double& cs, double& sn) {
	dlanv2_(&a, &b, &c, &d, &rt1r, &rt1i, &rt2r, &rt2i, &cs, &sn);
}

template <>
inline void tlasr<double>(const char side, const char pivot, const char direct, const int m, const int n,
	const double* c, const double* s, double* A, const int lda) {
	dlasr_(&side, &pivot, &direct, &m, &n, c, s, A, &lda);
}

template <>
inline void tlarf<double>(const char side, const int m, const int n, const double* v, const int incv,
	const double& tau, double* C, const int ldc) {
	const int lwork = tblas_detail::option_is(side, 'L') ? std::max(1, n) : std::max(1, m);
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dlarf_(&side, &m, &n, v, &incv, &tau, C, &ldc, work.data());
}

template <>
inline void tlarf1f<double>(const char side, const int m, const int n, const double* v, const int incv,
	const double& tau, double* C, const int ldc) {
	tlarf<double>(side, m, n, v, incv, tau, C, ldc);
}

template <>
inline void tlarf1l<double>(const char side, const int m, const int n, const double* v, const int incv,
	const double& tau, double* C, const int ldc) {
	tlarf<double>(side, m, n, v, incv, tau, C, ldc);
}

template <>
inline void tlarft<double>(const char direct, const char storev, const int n, const int k,
	const double* V, const int ldv, const double* tau, double* Tf, const int ldt) {
	dlarft_(&direct, &storev, &n, &k, V, &ldv, tau, Tf, &ldt);
}

template <>
inline void tlarfb<double>(const char side, const char trans, const char direct, const char storev,
	const int m, const int n, const int k, const double* V, const int ldv,
	const double* Tf, const int ldt, double* C, const int ldc) {
	const int ldwork = tblas_detail::option_is(side, 'L') ? std::max(1, n) : std::max(1, m);
	std::vector<double> work(static_cast<std::size_t>(ldwork) * std::max(1, k));
	dlarfb_(&side, &trans, &direct, &storev, &m, &n, &k, V, &ldv, Tf, &ldt, C, &ldc, work.data(), &ldwork);
}

template <>
inline int torg2r<double>(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	dorg2r_(&m, &n, &k, A, &lda, tau, work.data(), &info);
	return info;
}

template <>
inline int torm2r<double>(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	const int lwork = tblas_detail::option_is(side, 'L') ? std::max(1, n) : std::max(1, m);
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorm2r_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &info);
	return info;
}

template <>
inline int torgl2<double>(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, m)));
	dorgl2_(&m, &n, &k, A, &lda, tau, work.data(), &info);
	return info;
}

template <>
inline int torml2<double>(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	const int lwork = tblas_detail::option_is(side, 'L') ? std::max(1, n) : std::max(1, m);
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorml2_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &info);
	return info;
}

template <>
inline int torg2l<double>(const int m, const int n, const int k, double* A, const int lda, const double* tau) {
	int info = 0;
	std::vector<double> work(static_cast<std::size_t>(std::max(1, n)));
	dorg2l_(&m, &n, &k, A, &lda, tau, work.data(), &info);
	return info;
}

template <>
inline int torm2l<double>(const char side, const char trans, const int m, const int n, const int k,
	const double* A, const int lda, const double* tau, double* C, const int ldc) {
	int info = 0;
	const int lwork = tblas_detail::option_is(side, 'L') ? std::max(1, n) : std::max(1, m);
	std::vector<double> work(static_cast<std::size_t>(lwork));
	dorm2l_(&side, &trans, &m, &n, &k, A, &lda, tau, C, &ldc, work.data(), &info);
	return info;
}

template <>
inline int tlasyf<double>(const char uplo, const int n, const int nb, int& kb, double* A, const int lda,
	int* ipiv, double* W, const int ldw) {
	int info = 0;
	dlasyf_(&uplo, &n, &nb, &kb, A, &lda, ipiv, W, &ldw, &info);
	return info;
}

template <>
inline void tlatrd<double>(const char uplo, const int n, const int nb, double* A, const int lda,
	double* e, double* tau, double* W, const int ldw) {
	dlatrd_(&uplo, &n, &nb, A, &lda, e, tau, W, &ldw);
}

template <>
inline void tlahr2<double>(const int n, const int k, const int nb, double* A, const int lda,
	double* tau, double* Tf, const int ldt, double* Y, const int ldy) {
	dlahr2_(&n, &k, &nb, A, &lda, tau, Tf, &ldt, Y, &ldy);
}

template <>
inline void tlabrd<double>(const int m, const int n, const int nb, double* A, const int lda,
	double* d, double* e, double* tauq, double* taup, double* X, const int ldx, double* Y, const int ldy) {
	dlabrd_(&m, &n, &nb, A, &lda, d, e, tauq, taup, X, &ldx, Y, &ldy);
}

template <>
inline int tlahqr<double>(const bool wantt, const bool wantz, const int n, const int ilo, const int ihi,
	double* H, const int ldh, double* wr, double* wi, const int iloz, const int ihiz, double* Z, const int ldz) {
	int info = 0;
	const int iwantt = wantt ? 1 : 0;
	const int iwantz = wantz ? 1 : 0;
	dlahqr_(&iwantt, &iwantz, &n, &ilo, &ihi, H, &ldh, wr, wi, &iloz, &ihiz, Z, &ldz, &info);
	return info;
}

template <>
inline int tlaln2<double>(const bool ltrans, const int na, const int nw, const double& smin, const double& ca,
	const double* A, const int lda, const double& d1, const double& d2, const double* B, const int ldb,
	const double& wr, const double& wi, double* X, const int ldx, double& scale, double& xnorm) {
	int info = 0;
	const int iltrans = ltrans ? 1 : 0;
	dlaln2_(&iltrans, &na, &nw, &smin, &ca, A, &lda, &d1, &d2, B, &ldb, &wr, &wi, X, &ldx, &scale, &xnorm, &info);
	return info;
}
#endif // USE_VCP_LAPACK

} // namespace vcp

#endif // TLAPACK_TLAPACK_DOUBLE_HPP
