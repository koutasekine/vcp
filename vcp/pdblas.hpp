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

#ifndef VCP_PDBLAS_HPP
#define VCP_PDBLAS_HPP

#include <vcp/mats.hpp>

extern "C" {
	double ddot_(int*, const double*, int*, const double*, int*);
	void dgemm_(char*, char*, int*, int*, int*, double*, const double*, int*, const double*, int*, double*, double*, int*);
	void dsymm_(char*, char*, int*, int*, double*, const double*, int*, const double*, int*, double*, double*, int*);
	void dgemv_(char*, int*, int*, double*, const double*, int*, const double*, int*, double*, double*, int*);
	void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);

	void dgetrf_(int*, int*, double*, int*, int*, int*);
	void dgetri_(int*, double*, int*, int*, double*, int*, int*);
	void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
	void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
	void dsygv_(int*, char*, char*, int*, double*, int*, double*, int*, double*, double*, int*, int*);
	void dpotrf_(char*, int*, double*, int*, int*);
}

namespace vcp {
	class pdblas : public mats< double > {
	public:
		void mulmm(const pdblas& B, pdblas& c) const{
			c.row =(*this).row;
			c.column = B.column;
			c.n = c.row * c.column;
			c.v.resize(c.n);
			int inca = 1, incb = 1, incc = 1;
			char trans[] = "N";
			int Arow =(*this).row;
			int Acolumn = (*this).column;
			int NA =(*this).row * (*this).column;
			int Brow = B.row;
			int Bcolumn = B.column;
			double alpha = 1.0, beta = 0.0;
			if ((*this).type == 'S' && B.type == 'S') {
				c.type = 'S';
				c.v[0] = (*this).v[0] * B.v[0];
				return ;
			}
			else if ((*this).type == 'R' && B.type == 'C' && (*this).column == B.row) {
				c.type = 'S';
				c.v[0] = ddot_(&NA, &(*this).v.front(), &inca, &B.v.front(), &incb);
				return ;
			}
			else if ((*this).type == 'C' && B.type == 'R') {
				c.type = 'M';
				dgemm_(trans, trans, &Arow, &Bcolumn, &Acolumn, &alpha, &(*this).v.front(), &Arow, &B.v.front(), &Brow, &beta, &c.v.front(), &Arow);
				return ;
			}
			else if ((*this).type == 'M' && B.type == 'C' && (*this).column == B.row) {
				c.type = 'C';
				dgemv_(trans, &Arow, &Acolumn, &alpha, &(*this).v.front(), &Arow, &B.v.front(), &incb, &beta, &c.v.front(), &incc);
				return ;
			}
			else if ((*this).type == 'R' && B.type == 'M' && (*this).column == B.row) {
				c.type = 'R';
				dgemm_(trans, trans, &Arow, &Bcolumn, &Acolumn, &alpha, &(*this).v.front(), &Arow, &B.v.front(), &Brow, &beta, &c.v.front(), &Arow);
				return ;
			}
			else if ((*this).type == 'M' && B.type == 'M' && (*this).column == B.row) {
				c.type = 'M';
				dgemm_(trans, trans, &Arow, &Bcolumn, &Acolumn, &alpha, &(*this).v.front(), &Arow, &B.v.front(), &Brow, &beta, &c.v.front(), &Arow);
				return ;
			}
			else {
				std::cout << "*:error " << (*this).type << " , " << B.type << " , " << (*this).column << " , " << B.row << std::endl;
				exit(1);
			}
		}
		// C = transpose(A)*A : multiplication left side transpose
		void mulltmm(pdblas& c)const {
			c.row = column;
			c.column = column;
			c.n = c.row * c.column;
			c.v.resize(c.n);
			if ((*this).type == 'S') {
				using std::pow;
				c.type = 'S';
				c.v[0] = (*this).v[0] * (*this).v[0];
				return;
			}
			else if ((*this).type == 'C') {
				using std::pow;
				c.type = 'S';
				int inca = 1, incb = 1;
				int NA = (*this).row * (*this).column;
				c.v[0] = ddot_(&NA, &(*this).v.front(), &inca, &(*this).v.front(), &incb);
			}
			else if ((*this).type == 'M' || (*this).type == 'R') {
				c.type = 'M';
				char transN[] = "N";
				char transT[] = "T";
				int Arow = (*this).column;
				int Acolumn = (*this).row;
				int Brow = (*this).row;
				int Bcolumn = (*this).column;
				double alpha = 1.0, beta = 0.0;
				dgemm_(transT, transN, &Arow, &Bcolumn, &Acolumn, &alpha, &(*this).v.front(), &Acolumn, &(*this).v.front(), &Brow, &beta, &c.v.front(), &Arow);

				for (int i = 0; i < (*this).column; i++) {
					for (int j = i + 1; j < (*this).column; j++) {
						c.v[j + (*this).column*i] = c.v[i + (*this).column*j];
					}
				}
			}
			else {
				std::cout << "*:error " << (*this).type << std::endl;
				exit(1);
			}
		}

		void linearsolve(pdblas& b, pdblas& x) {
			if ((*this).row != (*this).column || b.row != (*this).row) {
				std::cout << "error : linearsolve row != column || b.row != row" << std::endl;
				exit(1);
			}

			int Asize = (*this).row;
			int Bsize = b.column;
			int lda = (*this).row, ldb = b.row;
			int* ipiv = new int[(*this).row];
			int info;

			x = b;
			dgesv_(&Asize, &Bsize, &(*this).v.front(), &lda, ipiv, &x.v.front(), &ldb, &info);
			delete[] ipiv;

			if (info != 0) {
				std::cout << "error : linearsolve :: info = " << info << std::endl;
				exit(1);
			}
		}

		void inv() override {
			if ((*this).row != (*this).column) {
				std::cout << "error : inv row != column" << std::endl;
				exit(1);
			}
			int* ipiv = new int[(*this).n];
			double* work = new double[(*this).n];
			int lwork = (*this).n;
			int lda = (*this).row;
			int info;
			
			dgetrf_(&(*this).row, &(*this).column, &(*this).v.front(), &lda, ipiv, &info);
			dgetri_(&(*this).row, &(*this).v.front(), &lda, ipiv, work, &lwork, &info);
			delete[] ipiv;
			delete[] work;

			if (info != 0) {
				std::cout << "error : inv :: info = " << info << std::endl;
				exit(1);
			}
		}

		void Cholesky() override {
			if (!(*this).is_symmetric()) {
				std::cout << "error : Cholesky : no symmetric matrix" << std::endl;
				exit(1);
			}
			char uplo[] = "U";
			int N = (*this).row;
			int lda = (*this).row;
			int info;
			dpotrf_(uplo, &N, &(*this).v.front(), &lda, &info);
			if (info != 0) {
				std::cout << "error : Cholesky :: info = " << info << std::endl;
				exit(1);
			}
			for (int i = 0; i < (*this).row; i++) {
				for (int j = 0; j < i; j++) {
					(*this).v[i + row * j] = 0.0;
				}
			}
		}

		void eigsym(int itep = 4) override {

			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix" << std::endl;
				exit(1);
			}

			char N[] = "N";
			char U[] = "U";

			double lw;
			int m1 = -1;
			int an = (*this).row;
			int lwork, info;

			pdblas E;
			E.zeros(an,1);
			dsyev_(N, U, &an, &(*this).v.front(), &an, &E.v.front(), &lw, &m1, &info);
			lwork = int(lw);
			double* work = new double[lwork];
			dsyev_(N, U, &an, &(*this).v.front(), &an, &E.v.front(), work, &lwork, &info);

			if (info != 0) {
				std::cout << "error : : eigsym :: info = " << info << std::endl;
			}
			(*this).zeros(an, an);

			for (int i = 0; i < an; i++) {
				(*this).v[i + (*this).row * i] = E.v[i];
			}
		}
		void eigsym(pdblas& V, int itep = 4) {
			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix" << std::endl;
				exit(1);
			}

			char VV[] = "V";
			char U[] = "U";
			double lw;
			int m1 = -1;
			int an = (*this).row;
			int lwork, info;

			pdblas E;
			E.zeros(an, 1);
			dsyev_(VV, U, &an, &(*this).v.front(), &an, &E.v.front(), &lw, &m1, &info);
			lwork = int(lw);
			double* work = new double[lwork];
			dsyev_(VV, U, &an, &(*this).v.front(), &an, &E.v.front(), work, &lwork, &info);
			if (info != 0) {
				std::cout << "error : : eigsym :: info = " << info << std::endl;
			}
			V = (*this);
			(*this).zeros(an, an);
			for (int i = 0; i < an; i++) {
				(*this).v[i + (*this).row * i] = E.v[i];
			}
			
		}

		void eigsymge(pdblas& B, int itep = 1) {
			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix A" << std::endl;
				exit(1);
			}
			if (!B.is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix B" << std::endl;
				exit(1);
			}
			if ((*this).row != B.row || (*this).column != B.column) {
				std::cout << "error : eigsym : no muching matrix size matrix A, B" << std::endl;
				exit(1);
			}
			
			int itype = 1;
			char jbobz[] = "N";
			char uplo[] = "U";
			int n = (*this).row;
			
			int lda = (*this).row;
			int ldb = B.row;
			double lw;
			int m1 = -1;			
			int lwork, info;

			pdblas W;
			W.zeros((*this).row, 1);

			dsygv_(&itype, jbobz, uplo, &n, &(*this).v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), &lw, &m1, &info);
			lwork = int(lw);
			double* work = new double[lwork];
			dsygv_(&itype, jbobz, uplo, &n, &(*this).v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), work, &lwork, &info);
			if (info != 0) {
				std::cout << "error : : eigsymge :: info = " << info << std::endl;
			}
			(*this).zeros((*this).row, (*this).row);

			for (int i = 0; i < (*this).row; i++) {
				(*this).v[i + (*this).row * i] = W.v[i];
			}
		}
		void eigsymge(pdblas& B, pdblas& V, int itep = 1) {
				if (!(*this).is_symmetric()) {
					std::cout << "error : eigsym : no symmetric matrix A" << std::endl;
					exit(1);
				}
				if (!B.is_symmetric()) {
					std::cout << "error : eigsym : no symmetric matrix B" << std::endl;
					exit(1);
				}
				if ((*this).row != B.row || (*this).column != B.column) {
					std::cout << "error : eigsym : no muching matrix size matrix A, B" << std::endl;
					exit(1);
				}

				int itype = 1;
				char jbobz[] = "V";
				char uplo[] = "U";
				int n = (*this).row;

				int lda = (*this).row;
				int ldb = B.row;
				double lw;
				int m1 = -1;
				int lwork, info;

				pdblas W;
				W.zeros((*this).row, 1);

				dsygv_(&itype, jbobz, uplo, &n, &(*this).v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), &lw, &m1, &info);
				lwork = int(lw);
				double* work = new double[lwork];
				dsygv_(&itype, jbobz, uplo, &n, &(*this).v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), work, &lwork, &info);
				if (info != 0) {
					std::cout << "error : : eigsymge :: info = " << info << std::endl;
				}
				V = (*this);
				(*this).zeros((*this).row, (*this).row);
				for (int i = 0; i < (*this).row; i++) {
					(*this).v[i + (*this).row * i] = W.v[i];
				}
		}
	};
}

#endif // VCP_PDBLAS_HPP
