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

#ifndef VCP_DBLAS_DLAPACK_HPP
#define VCP_DBLAS_DLAPACK_HPP

#ifdef USE_VCP_BLAS
#  include <vcp/vblas/dblas.hpp>
#endif

#ifdef USE_VCP_LAPACK
#  include <vcp/vlapack/dlapack.hpp>
#endif

#if !defined(USE_VCP_BLAS) || !defined(USE_VCP_LAPACK)
extern "C" {
#  ifndef USE_VCP_BLAS
	void dcopy_(const int*, const double*, const int*, double*, const int*);
	void dswap_(const int*, double*, const int*, double*, const int*);
	int idamax_(const int*, const double*, const int*);
	void dscal_(const int*, const double*, double*, const int*);
	void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
	double ddot_(const int*, const double*, const int*, const double*, const int*);
	double dasum_(const int*, const double*, const int*);
	double dnrm2_(const int*, const double*, const int*);
	void drot_(const int*, double*, const int*, double*, const int*, const double*, const double*);
	void drotg_(double*, double*, double*, double*);
	void drotm_(const int*, double*, const int*, double*, const int*, const double*);
	void drotmg_(double*, double*, double*, const double*, double*);

	void dgemv_(const char*, const int*, const int*, const double*, const double*, const int*,
		const double*, const int*, const double*, double*, const int*);
	void dgbmv_(const char*, const int*, const int*, const int*, const int*, const double*,
		const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void dsymv_(const char*, const int*, const double*, const double*, const int*,
		const double*, const int*, const double*, double*, const int*);
	void dsbmv_(const char*, const int*, const int*, const double*, const double*, const int*,
		const double*, const int*, const double*, double*, const int*);
	void dspmv_(const char*, const int*, const double*, const double*, const double*,
		const int*, const double*, double*, const int*);
	void dtrmv_(const char*, const char*, const char*, const int*, const double*,
		const int*, double*, const int*);
	void dtbmv_(const char*, const char*, const char*, const int*, const int*, const double*,
		const int*, double*, const int*);
	void dtpmv_(const char*, const char*, const char*, const int*, const double*, double*, const int*);
	void dtrsv_(const char*, const char*, const char*, const int*, const double*,
		const int*, double*, const int*);
	void dtbsv_(const char*, const char*, const char*, const int*, const int*, const double*,
		const int*, double*, const int*);
	void dtpsv_(const char*, const char*, const char*, const int*, const double*, double*, const int*);
	void dger_(const int*, const int*, const double*, const double*, const int*,
		const double*, const int*, double*, const int*);
	void dsyr_(const char*, const int*, const double*, const double*, const int*, double*, const int*);
	void dspr_(const char*, const int*, const double*, const double*, const int*, double*);
	void dsyr2_(const char*, const int*, const double*, const double*, const int*,
		const double*, const int*, double*, const int*);
	void dspr2_(const char*, const int*, const double*, const double*, const int*,
		const double*, const int*, double*);

	void dgemm_(const char*, const char*, const int*, const int*, const int*, const double*,
		const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void dgemmtr_(const char*, const char*, const char*, const int*, const int*, const double*,
		const double*, const int*, const double*, const int*, const double*, double*, const int*);
	void dsymm_(const char*, const char*, const int*, const int*, const double*, const double*,
		const int*, const double*, const int*, const double*, double*, const int*);
	void dsyrk_(const char*, const char*, const int*, const int*, const double*, const double*,
		const int*, const double*, double*, const int*);
	void dsyr2k_(const char*, const char*, const int*, const int*, const double*, const double*,
		const int*, const double*, const int*, const double*, double*, const int*);
	void dtrmm_(const char*, const char*, const char*, const char*, const int*, const int*,
		const double*, const double*, const int*, double*, const int*);
	void dtrsm_(const char*, const char*, const char*, const char*, const int*, const int*,
		const double*, const double*, const int*, double*, const int*);
#  endif

#  ifndef USE_VCP_LAPACK
	double dlamch_(const char*);
	double dlapy2_(const double*, const double*);
	double dlapy3_(const double*, const double*, const double*);
	void dlassq_(const int*, const double*, const int*, double*, double*);
	void dlartg_(const double*, const double*, double*, double*, double*);
	void dlarfg_(const int*, double*, double*, const int*, double*);
	void dlascl_(const char*, const int*, const int*, const double*, const double*,
		const int*, const int*, double*, const int*, int*);
	void dlaset_(const char*, const int*, const int*, const double*, const double*, double*, const int*);
	void dlacpy_(const char*, const int*, const int*, const double*, const int*, double*, const int*);
	void dlasrt_(const char*, const int*, double*, int*);
	void dlaswp_(const int*, double*, const int*, const int*, const int*, const int*, const int*);
	double dlange_(const char*, const int*, const int*, const double*, const int*, double*);
	double dlangb_(const char*, const int*, const int*, const int*, const double*, const int*, double*);
	double dlansy_(const char*, const char*, const int*, const double*, const int*, double*);
	double dlansb_(const char*, const char*, const int*, const int*, const double*, const int*, double*);
	double dlanst_(const char*, const int*, const double*, const double*);

	void dgetf2_(const int*, const int*, double*, const int*, int*, int*);
	void dgetrf_(const int*, const int*, double*, const int*, int*, int*);
	void dgetrs_(const char*, const int*, const int*, const double*, const int*,
		const int*, double*, const int*, int*);
	void dgesv_(const int*, const int*, double*, const int*, int*, double*, const int*, int*);
	void dgetri_(const int*, double*, const int*, const int*, double*, const int*, int*);
	void dtrti2_(const char*, const char*, const int*, double*, const int*, int*);
	void dtrtri_(const char*, const char*, const int*, double*, const int*, int*);
	void dtrtrs_(const char*, const char*, const char*, const int*, const int*,
		const double*, const int*, double*, const int*, int*);
	void dlauu2_(const char*, const int*, double*, const int*, int*);
	void dlauum_(const char*, const int*, double*, const int*, int*);

	void dpotf2_(const char*, const int*, double*, const int*, int*);
	void dpotrf_(const char*, const int*, double*, const int*, int*);
	void dpotrs_(const char*, const int*, const int*, const double*, const int*, double*, const int*, int*);
	void dposv_(const char*, const int*, const int*, double*, const int*, double*, const int*, int*);
	void dpotri_(const char*, const int*, double*, const int*, int*);

	void dsytf2_(const char*, const int*, double*, const int*, int*, int*);
	void dsytrf_(const char*, const int*, double*, const int*, int*, double*, const int*, int*);
	void dsytrs_(const char*, const int*, const int*, const double*, const int*,
		const int*, double*, const int*, int*);
	void dsysv_(const char*, const int*, const int*, double*, const int*, int*,
		double*, const int*, double*, const int*, int*);

	void dgbtf2_(const int*, const int*, const int*, const int*, double*, const int*, int*, int*);
	void dgbtrf_(const int*, const int*, const int*, const int*, double*, const int*, int*, int*);
	void dgbtrs_(const char*, const int*, const int*, const int*, const int*,
		const double*, const int*, const int*, double*, const int*, int*);
	void dgbsv_(const int*, const int*, const int*, const int*, double*, const int*,
		int*, double*, const int*, int*);
	void dpbtf2_(const char*, const int*, const int*, double*, const int*, int*);
	void dpbtrf_(const char*, const int*, const int*, double*, const int*, int*);
	void dpbtrs_(const char*, const int*, const int*, const int*, const double*,
		const int*, double*, const int*, int*);
	void dpbsv_(const char*, const int*, const int*, const int*, double*, const int*,
		double*, const int*, int*);

	void dgeqr2_(const int*, const int*, double*, const int*, double*, double*, int*);
	void dgeqrf_(const int*, const int*, double*, const int*, double*, double*, const int*, int*);
	void dgelq2_(const int*, const int*, double*, const int*, double*, double*, int*);
	void dgelqf_(const int*, const int*, double*, const int*, double*, double*, const int*, int*);
	void dorgqr_(const int*, const int*, const int*, double*, const int*,
		const double*, double*, const int*, int*);
	void dorglq_(const int*, const int*, const int*, double*, const int*,
		const double*, double*, const int*, int*);
	void dorgql_(const int*, const int*, const int*, double*, const int*,
		const double*, double*, const int*, int*);
	void dormqr_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, const int*, int*);
	void dormlq_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, const int*, int*);
	void dormql_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, const int*, int*);
	void dgels_(const char*, const int*, const int*, const int*, double*, const int*,
		double*, const int*, double*, const int*, int*);

	void dsytd2_(const char*, const int*, double*, const int*, double*, double*, double*, int*);
	void dsytrd_(const char*, const int*, double*, const int*, double*, double*, double*,
		double*, const int*, int*);
	void dorgtr_(const char*, const int*, double*, const int*, const double*, double*, const int*, int*);
	void dormtr_(const char*, const char*, const char*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, const int*, int*);
	void dsterf_(const int*, double*, double*, int*);
	void dsteqr_(const char*, const int*, double*, double*, double*, const int*, double*, int*);
	void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);
	void dsygs2_(const int*, const char*, const int*, double*, const int*, const double*, const int*, int*);
	void dsygst_(const int*, const char*, const int*, double*, const int*, const double*, const int*, int*);
	void dsygv_(const int*, const char*, const char*, const int*, double*, const int*,
		double*, const int*, double*, double*, const int*, int*);

	void dgebd2_(const int*, const int*, double*, const int*, double*, double*, double*, double*, double*, int*);
	void dgebrd_(const int*, const int*, double*, const int*, double*, double*, double*,
		double*, double*, const int*, int*);
	void dorgbr_(const char*, const int*, const int*, const int*, double*, const int*,
		const double*, double*, const int*, int*);
	void dormbr_(const char*, const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, const double*, double*, const int*, double*, const int*, int*);
	void dbdsqr_(const char*, const int*, const int*, const int*, const int*,
		double*, double*, double*, const int*, double*, const int*, double*, const int*, double*, int*);
	void dgesvd_(const char*, const char*, const int*, const int*, double*, const int*,
		double*, double*, const int*, double*, const int*, double*, const int*, int*);

	void dgebal_(const char*, const int*, double*, const int*, int*, int*, double*, int*);
	void dgebak_(const char*, const char*, const int*, const int*, const int*,
		const double*, const int*, double*, const int*, int*);
	void dgehd2_(const int*, const int*, const int*, double*, const int*, double*, double*, int*);
	void dgehrd_(const int*, const int*, const int*, double*, const int*, double*, double*, const int*, int*);
	void dorghr_(const int*, const int*, const int*, double*, const int*,
		const double*, double*, const int*, int*);
	void dhseqr_(const char*, const char*, const int*, const int*, const int*,
		double*, const int*, double*, double*, double*, const int*, double*, const int*, int*);
	void dtrevc_(const char*, const char*, const int*, const int*, const double*, const int*,
		double*, const int*, double*, const int*, const int*, int*, double*, int*);
	void dgeev_(const char*, const char*, const int*, double*, const int*, double*, double*,
		double*, const int*, double*, const int*, double*, const int*, int*);
#  endif
}
#endif

#endif // VCP_DBLAS_DLAPACK_HPP
