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
#include <vcp/dblas_dlapack.hpp>

namespace vcp {
	class pdblas : public mats< double > {
	public:
		void mulmm(const pdblas& B, pdblas& c) const{
			if (this->type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c = B;
				for (int i = 0; i < B.n; i++) {
					c.v[i] *= this->v[0];
				}
				return;
			}
			else if ((this->type == 'C' || this->type == 'R' || this->type == 'M') && B.type == 'S') {
				c = *this;
				for (int i = 0; i < this->n; i++) {
					c.v[i] *= B.v[0];
				}
				return;
			}
			c.row = this->row;
			c.column = B.column;
			c.n = c.row * c.column;
			c.v.resize(c.n);
			int inca = 1, incb = 1, incc = 1;
			char trans[] = "N";
			int Arow = this->row;
			int Acolumn = this->column;
			int NA = this->row * this->column;
			int Brow = B.row;
			int Bcolumn = B.column;
			double alpha = 1.0, beta = 0.0;
			if (this->type == 'S' && B.type == 'S') {
				c.type = 'S';
				c.v[0] = this->v[0] * B.v[0];
				return ;
			}
			else if (this->type == 'R' && B.type == 'C' && this->column == B.row) {
				c.type = 'S';
				c.v[0] = ddot_(&NA, &this->v.front(), &inca, &B.v.front(), &incb);
				return ;
			}
			else if (this->type == 'C' && B.type == 'R') {
				c.type = 'M';
				dgemm_(trans, trans, &Arow, &Bcolumn, &Acolumn, &alpha, &this->v.front(), &Arow, &B.v.front(), &Brow, &beta, &c.v.front(), &Arow);
				return ;
			}
			else if (this->type == 'M' && B.type == 'C' && this->column == B.row) {
				c.type = 'C';
				dgemv_(trans, &Arow, &Acolumn, &alpha, &this->v.front(), &Arow, &B.v.front(), &incb, &beta, &c.v.front(), &incc);
				return ;
			}
			else if (this->type == 'R' && B.type == 'M' && this->column == B.row) {
				c.type = 'R';
				dgemm_(trans, trans, &Arow, &Bcolumn, &Acolumn, &alpha, &this->v.front(), &Arow, &B.v.front(), &Brow, &beta, &c.v.front(), &Arow);
				return ;
			}
			else if (this->type == 'M' && B.type == 'M' && this->column == B.row) {
				c.type = 'M';
				dgemm_(trans, trans, &Arow, &Bcolumn, &Acolumn, &alpha, &this->v.front(), &Arow, &B.v.front(), &Brow, &beta, &c.v.front(), &Arow);
				return ;
			}
			else {
				vcp::throw_error<vcp::dimension_error>(
					"mulmm: unsupported matrix product: type ", this->type, " * ", B.type,
					", size (", this->row, ", ", this->column, ") * (", B.row, ", ", B.column, ")");
			}
		}
		// C = transpose(A)*A : multiplication left side transpose
		void mulltmm(pdblas& c)const {
			c.row = column;
			c.column = column;
			c.n = c.row * c.column;
			c.v.resize(c.n);
			if (this->type == 'S') {
				using std::pow;
				c.type = 'S';
				c.v[0] = this->v[0] * this->v[0];
				return;
			}
			else if (this->type == 'C') {
				using std::pow;
				c.type = 'S';
				int inca = 1, incb = 1;
				int NA = this->row * this->column;
				c.v[0] = ddot_(&NA, &this->v.front(), &inca, &this->v.front(), &incb);
			}
			else if (this->type == 'M' || this->type == 'R') {
				c.type = 'M';
				char uplo[] = "U";
				char transT[] = "T";
				int n = this->column;
				int k = this->row;
				int lda = this->row;
				int ldc = this->column;
				double alpha = 1.0, beta = 0.0;
				dsyrk_(uplo, transT, &n, &k, &alpha, &this->v.front(), &lda, &beta, &c.v.front(), &ldc);

				for (int i = 0; i < this->column; i++) {
					for (int j = i + 1; j < this->column; j++) {
						c.v[j + this->column*i] = c.v[i + this->column*j];
					}
				}
			}
			else {
				vcp::throw_error<vcp::state_error>("mulltmm: unsupported matrix type: ", this->type);
			}
		}

		void linearsolve(pdblas& b, pdblas& x) {
			if (this->row != this->column || b.row != this->row) {
				vcp::throw_error<vcp::dimension_error>(
					"linearsolve: invalid dimensions: A=(", this->row, ", ", this->column,
					"), b=(", b.row, ", ", b.column, ")");
			}

			int Asize = this->row;
			int Bsize = b.column;
			int lda = this->row, ldb = b.row;
			std::vector<int> ipiv(this->row);
			int info;

			x = b;
			dgesv_(&Asize, &Bsize, &this->v.front(), &lda, ipiv.data(), &x.v.front(), &ldb, &info);

			if (info != 0) {
				throw vcp::lapack_error("dgesv", info, "linearsolve: dgesv failed");
			}
		}

		void inv() override {
			if (this->row != this->column) {
				vcp::throw_error<vcp::dimension_error>(
					"inv: matrix must be square: ", this->row, " != ", this->column);
			}
			std::vector<int> ipiv(this->row);
			double lw;
			int lwork = -1;
			int lda = this->row;
			int info;
			
			dgetrf_(&this->row, &this->column, &this->v.front(), &lda, ipiv.data(), &info);
			if (info != 0) {
				throw vcp::lapack_error("dgetrf", info, "inv: dgetrf failed");
			}
			dgetri_(&this->row, &this->v.front(), &lda, ipiv.data(), &lw, &lwork, &info);
			if (info != 0) {
				throw vcp::lapack_error("dgetri", info, "inv: dgetri workspace query failed");
			}
			lwork = int(lw);
			std::vector<double> work(lwork);
			dgetri_(&this->row, &this->v.front(), &lda, ipiv.data(), work.data(), &lwork, &info);
			if (info != 0) {
				throw vcp::lapack_error("dgetri", info, "inv: dgetri failed");
			}
		}

		void Cholesky() override {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("Cholesky: matrix must be symmetric");
			}
			char uplo[] = "U";
			int N = this->row;
			int lda = this->row;
			int info;
			dpotrf_(uplo, &N, &this->v.front(), &lda, &info);
			if (info != 0) {
				throw vcp::lapack_error("dpotrf", info, "Cholesky: dpotrf failed");
			}
			for (int i = 0; i < this->row; i++) {
				for (int j = 0; j < i; j++) {
					this->v[i + row * j] = 0.0;
				}
			}
		}

		void eigsym(int itep = 4) override {

			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("eigsym: matrix must be symmetric");
			}

			char N[] = "N";
			char U[] = "U";

			double lw;
			int m1 = -1;
			int an = this->row;
			int lwork, info;

			pdblas E;
			E.zeros(an,1);
			dsyev_(N, U, &an, &this->v.front(), &an, &E.v.front(), &lw, &m1, &info);
			if (info != 0) {
				throw vcp::lapack_error("dsyev", info, "eigsym: dsyev workspace query failed");
			}
			lwork = int(lw);
			std::vector<double> work(lwork);
			dsyev_(N, U, &an, &this->v.front(), &an, &E.v.front(), work.data(), &lwork, &info);

			if (info != 0) {
				throw vcp::lapack_error("dsyev", info, "eigsym: dsyev failed");
			}
			this->zeros(an, an);

			for (int i = 0; i < an; i++) {
				this->v[i + this->row * i] = E.v[i];
			}
		}
		void eigsym(pdblas& V, int itep = 4) {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("eigsym: matrix must be symmetric");
			}

			char VV[] = "V";
			char U[] = "U";
			double lw;
			int m1 = -1;
			int an = this->row;
			int lwork, info;

			pdblas E;
			E.zeros(an, 1);
			dsyev_(VV, U, &an, &this->v.front(), &an, &E.v.front(), &lw, &m1, &info);
			if (info != 0) {
				throw vcp::lapack_error("dsyev", info, "eigsym: dsyev workspace query failed");
			}
			lwork = int(lw);
			std::vector<double> work(lwork);
			dsyev_(VV, U, &an, &this->v.front(), &an, &E.v.front(), work.data(), &lwork, &info);
			if (info != 0) {
				throw vcp::lapack_error("dsyev", info, "eigsym: dsyev failed");
			}
			V = (*this);
			this->zeros(an, an);
			for (int i = 0; i < an; i++) {
				this->v[i + this->row * i] = E.v[i];
			}
			
		}

		void eigsymge(pdblas& B, int itep = 1) {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("eigsymge: matrix A must be symmetric");
			}
			if (!B.is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("eigsymge: matrix B must be symmetric");
			}
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>(
					"eigsymge: dimension mismatch: A=(", this->row, ", ", this->column,
					"), B=(", B.row, ", ", B.column, ")");
			}
			
			int itype = 1;
			char jbobz[] = "N";
			char uplo[] = "U";
			int n = this->row;
			
			int lda = this->row;
			int ldb = B.row;
			double lw;
			int m1 = -1;			
			int lwork, info;

			pdblas W;
			W.zeros(this->row, 1);

			dsygv_(&itype, jbobz, uplo, &n, &this->v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), &lw, &m1, &info);
			if (info != 0) {
				throw vcp::lapack_error("dsygv", info, "eigsymge: dsygv workspace query failed");
			}
			lwork = int(lw);
			std::vector<double> work(lwork);
			dsygv_(&itype, jbobz, uplo, &n, &this->v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), work.data(), &lwork, &info);
			if (info != 0) {
				throw vcp::lapack_error("dsygv", info, "eigsymge: dsygv failed");
			}
			this->zeros(this->row, this->row);

			for (int i = 0; i < this->row; i++) {
				this->v[i + this->row * i] = W.v[i];
			}
		}
		void eigsymge(pdblas& B, pdblas& V, int itep = 1) {
				if (!this->is_symmetric()) {
					vcp::throw_error<vcp::domain_error>("eigsymge: matrix A must be symmetric");
				}
				if (!B.is_symmetric()) {
					vcp::throw_error<vcp::domain_error>("eigsymge: matrix B must be symmetric");
				}
				if (this->row != B.row || this->column != B.column) {
					vcp::throw_error<vcp::dimension_error>(
						"eigsymge: dimension mismatch: A=(", this->row, ", ", this->column,
						"), B=(", B.row, ", ", B.column, ")");
				}

				int itype = 1;
				char jbobz[] = "V";
				char uplo[] = "U";
				int n = this->row;

				int lda = this->row;
				int ldb = B.row;
				double lw;
				int m1 = -1;
				int lwork, info;

				pdblas W;
				W.zeros(this->row, 1);

				dsygv_(&itype, jbobz, uplo, &n, &this->v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), &lw, &m1, &info);
				if (info != 0) {
					throw vcp::lapack_error("dsygv", info, "eigsymge: dsygv workspace query failed");
				}
				lwork = int(lw);
				std::vector<double> work(lwork);
				dsygv_(&itype, jbobz, uplo, &n, &this->v.front(), &lda, &B.v.front(), &ldb, &W.v.front(), work.data(), &lwork, &info);
				if (info != 0) {
					throw vcp::lapack_error("dsygv", info, "eigsymge: dsygv failed");
				}
				V = (*this);
				this->zeros(this->row, this->row);
				for (int i = 0; i < this->row; i++) {
					this->v[i + this->row * i] = W.v[i];
				}
		}
	};
}

#endif // VCP_PDBLAS_HPP
