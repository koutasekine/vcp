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

#ifndef VCP_PDDBLAS_HPP
#define VCP_PDDBLAS_HPP

#include <cmath>
#include <vector>
#include <utility>

#include <kv/dd.hpp>

#include <vcp/pdblas.hpp>

extern "C" {
	void dgetrs_(char*, int*, int*, double*, int*, int*, double*, int*, int*);
	void dtrsm_(char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);
}

namespace vcp {
	// Helper functions for the Ozaki scheme:
	// A dd matrix is split into a sum of double matrices ("slices") such that
	// every dgemm of two slices is error-free. The exact dgemm results are
	// accumulated in dd arithmetic.
	namespace pddblas_assist {
		inline int ceil_log2(const int k) {
			int r = 0;
			int t = 1;
			while (t < k) {
				t *= 2;
				r++;
			}
			return r;
		}
		// Shift s of the splitting: with rho = ceil(log2(k)) the slice products
		// summed over k terms remain exact iff 104 + rho - 2*s <= 53.
		inline int split_shift(const int k) {
			return (51 + ceil_log2(k) + 1) / 2;
		}
		// Significant bits carried by one slice
		inline int split_bits(const int k) {
			int b = 52 - split_shift(k);
			if (b < 1) b = 1;
			return b;
		}
		// Number of slices so that the dropped part is below the dd precision
		// (106 bits) with some margin
		inline int max_slices(const int k) {
			return 110 / split_bits(k) + 2;
		}

		// Row-wise split (one sigma per row) of an (m x k) dd array (column-major).
		// Used for the left factor of op(A)*B with op(A) = A.
		inline void split_row(const std::vector< kv::dd >& A, const int m, const int k, std::vector< vcp::pdblas >& D) {
			using std::fabs;
			using std::frexp;
			using std::ldexp;
			const int s = split_shift(k);
			const int smax = max_slices(k);
			std::vector< kv::dd > T = A;
			D.clear();
			for (int p = 0; p < smax; p++) {
				vcp::pdblas S;
				S.zeros(m, k);
				bool nonzero = false;
				for (int i = 0; i < m; i++) {
					double mu = 0.0;
					for (int j = 0; j < k; j++) {
						const double a = fabs(T[i + m * j].a1);
						if (a > mu) mu = a;
					}
					if (mu == 0.0) continue;
					int e;
					frexp(mu, &e); // 2^(e-1) <= mu < 2^e
					const double sigma = ldexp(1.0, e + s);
					for (int j = 0; j < k; j++) {
						const double q = (T[i + m * j].a1 + sigma) - sigma;
						if (q != 0.0) {
							S.v[i + m * j] = q;
							T[i + m * j] -= kv::dd(q);
							nonzero = true;
						}
					}
				}
				if (!nonzero) break;
				D.push_back(std::move(S));
			}
		}
		// Column-wise split (one sigma per column) of a (k x n) dd array.
		// Used for the right factor B, and (transposed) for the left factor of
		// transpose(A)*B since a column of A is a row of transpose(A).
		inline void split_col(const std::vector< kv::dd >& A, const int k, const int n, std::vector< vcp::pdblas >& D) {
			using std::fabs;
			using std::frexp;
			using std::ldexp;
			const int s = split_shift(k);
			const int smax = max_slices(k);
			std::vector< kv::dd > T = A;
			D.clear();
			for (int p = 0; p < smax; p++) {
				vcp::pdblas S;
				S.zeros(k, n);
				bool nonzero = false;
				for (int j = 0; j < n; j++) {
					double mu = 0.0;
					for (int i = 0; i < k; i++) {
						const double a = fabs(T[i + k * j].a1);
						if (a > mu) mu = a;
					}
					if (mu == 0.0) continue;
					int e;
					frexp(mu, &e);
					const double sigma = ldexp(1.0, e + s);
					for (int i = 0; i < k; i++) {
						const double q = (T[i + k * j].a1 + sigma) - sigma;
						if (q != 0.0) {
							S.v[i + k * j] = q;
							T[i + k * j] -= kv::dd(q);
							nonzero = true;
						}
					}
				}
				if (!nonzero) break;
				D.push_back(std::move(S));
			}
		}

		// C (dd, m x n) += sum over slice pairs of op(DA[p]) * DB[q].
		// op(DA[p]) is m x k (stored k x m if transA), DB[q] is k x n.
		// Each dgemm is exact by construction of the splitting; the exact double
		// results are summed in dd arithmetic, larger pairs (p+q small) first.
		inline void ozaki_acc(const std::vector< vcp::pdblas >& DA, const bool transA,
			const std::vector< vcp::pdblas >& DB,
			const int m, const int n, const int k,
			std::vector< kv::dd >& C) {
			if (DA.empty() || DB.empty()) return;
			char TA[] = "N";
			if (transA) TA[0] = 'T';
			char TB[] = "N";
			double alpha = 1.0, beta = 0.0;
			int M = m, N = n, K = k;
			int lda = transA ? k : m;
			int ldb = k, ldc = m;
			const int limit = max_slices(k);
			const int mn = m * n;
			vcp::pdblas P;
			P.zeros(m, n);
			for (int t = 0; t < limit; t++) {
				for (int p = 0; p <= t; p++) {
					const int q = t - p;
					if (p >= static_cast<int>(DA.size()) || q >= static_cast<int>(DB.size())) continue;
					dgemm_(TA, TB, &M, &N, &K, &alpha, &DA[p].v.front(), &lda, &DB[q].v.front(), &ldb, &beta, &P.v.front(), &ldc);
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
					#pragma omp parallel for
#endif
#endif
					for (int idx = 0; idx < mn; idx++) {
						if (P.v[idx] != 0.0) {
							C[idx] += kv::dd(P.v[idx]);
						}
					}
				}
			}
		}
	}

	// Policy for fast approximate kv::dd computations based on double BLAS/LAPACK:
	//  - matrix product       : Ozaki scheme (error-free dgemm slices)
	//  - linear system        : LU in double (dgetrf/dgetrs) + residual iteration
	//  - eigenvalue problems  : dsyev/dsygv in double + Ogita-Aishima refinement
	//    (T. Ogita, K. Aishima, Japan J. Indust. Appl. Math. (2018),
	//     https://link.springer.com/article/10.1007/s13160-018-0310-3)
	// All double computations are managed through the vcp::pdblas policy or the
	// BLAS/LAPACK routines already declared in pdblas.hpp.
	class pddblas : public mats< kv::dd > {
	protected:
		static void dd_to_d(const vcp::mats< kv::dd >& A, vcp::pdblas& B) {
			B.row = A.row;
			B.column = A.column;
			B.n = A.n;
			B.type = A.type;
			B.v.resize(A.n);
			for (int i = 0; i < A.n; i++) {
				B.v[i] = A.v[i].a1;
			}
		}

		// One step of the Ogita-Aishima refinement for the symmetric (standard or
		// generalized) eigenvalue problem. DA / DB are row-wise Ozaki splits of A
		// and B (DB = NULL for the standard problem B = I). X (n x n, columns are
		// approximate eigenvectors with X^T B X ~ I) and lam are updated.
		static void oa_step(const std::vector< vcp::pdblas >& DA,
			const std::vector< vcp::pdblas >* DB,
			std::vector< kv::dd >& X, std::vector< kv::dd >& lam, const int n) {
			using std::fabs;
			const kv::dd dd0 = kv::dd(0.0);
			const int nn = n * n;

			std::vector< vcp::pdblas > DX;
			pddblas_assist::split_col(X, n, n, DX);

			// AX = A*X,  S = X^T * AX
			std::vector< kv::dd > AX(nn, dd0);
			pddblas_assist::ozaki_acc(DA, false, DX, n, n, n, AX);
			std::vector< vcp::pdblas > DAX;
			pddblas_assist::split_col(AX, n, n, DAX);
			std::vector< kv::dd > S(nn, dd0);
			pddblas_assist::ozaki_acc(DX, true, DAX, n, n, n, S);

			// G = X^T * B * X
			std::vector< kv::dd > G(nn, dd0);
			if (DB != NULL) {
				std::vector< kv::dd > BX(nn, dd0);
				pddblas_assist::ozaki_acc(*DB, false, DX, n, n, n, BX);
				std::vector< vcp::pdblas > DBX;
				pddblas_assist::split_col(BX, n, n, DBX);
				pddblas_assist::ozaki_acc(DX, true, DBX, n, n, n, G);
			}
			else {
				pddblas_assist::ozaki_acc(DX, true, DX, n, n, n, G);
			}

			// lambda_i = S_ii / (X^T B X)_ii,   R = I - G
			for (int i = 0; i < n; i++) {
				lam[i] = S[i + n * i] / G[i + n * i];
			}
			std::vector< kv::dd >& R = G;
			for (int idx = 0; idx < nn; idx++) {
				R[idx] = -R[idx];
			}
			for (int i = 0; i < n; i++) {
				R[i + n * i] += kv::dd(1.0);
			}

			// Threshold delta for nearly multiple eigenvalues
			double normR = 0.0, offS = 0.0, lmax = 0.0;
			for (int i = 0; i < n; i++) {
				double s1 = 0.0, s2 = 0.0;
				for (int j = 0; j < n; j++) {
					s1 += fabs(R[i + n * j].a1);
					if (i != j) s2 += fabs(S[i + n * j].a1);
				}
				if (s1 > normR) normR = s1;
				if (s2 > offS) offS = s2;
				const double al = fabs(lam[i].a1);
				if (al > lmax) lmax = al;
			}
			const double delta = 2.0 * (offS + normR * lmax);

			// E_ij = (S_ij + lambda_j R_ij)/(lambda_j - lambda_i)  (separated case)
			// E_ij = R_ij/2                                        (otherwise)
			std::vector< kv::dd > E(nn, dd0);
			const kv::dd half = kv::dd(0.5);
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < n; i++) {
					if (i == j) {
						E[i + n * i] = R[i + n * i] * half;
					}
					else {
						kv::dd d = lam[j] - lam[i];
						if (fabs(d.a1) > delta) {
							E[i + n * j] = (S[i + n * j] + lam[j] * R[i + n * j]) / d;
						}
						else {
							E[i + n * j] = R[i + n * j] * half;
						}
					}
				}
			}

			// X = X + X*E
			std::vector< vcp::pdblas > DE, DXr;
			pddblas_assist::split_col(E, n, n, DE);
			pddblas_assist::split_row(X, n, n, DXr);
			std::vector< kv::dd > XE(nn, dd0);
			pddblas_assist::ozaki_acc(DXr, false, DE, n, n, n, XE);
			for (int idx = 0; idx < nn; idx++) {
				X[idx] += XE[idx];
			}
		}

	public:
		// C = A*B by the Ozaki scheme
		void mulmm(const pddblas& B, pddblas& c) const {
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
			else if (this->type == 'S' && B.type == 'S') {
				c.zeros(1, 1);
				c.v[0] = this->v[0] * B.v[0];
				return;
			}

			char ctype;
			if (this->type == 'R' && B.type == 'C' && this->column == B.row) {
				ctype = 'S';
			}
			else if (this->type == 'C' && B.type == 'R') {
				ctype = 'M';
			}
			else if (this->type == 'M' && B.type == 'C' && this->column == B.row) {
				ctype = 'C';
			}
			else if (this->type == 'R' && B.type == 'M' && this->column == B.row) {
				ctype = 'R';
			}
			else if (this->type == 'M' && B.type == 'M' && this->column == B.row) {
				ctype = 'M';
			}
			else {
				vcp::throw_error<vcp::dimension_error>(
					"mulmm: unsupported matrix product: type ", this->type, " * ", B.type,
					", size (", this->row, ", ", this->column, ") * (", B.row, ", ", B.column, ")");
				return;
			}

			const int m = this->row;
			const int k = this->column;
			const int nc = B.column;
			std::vector< vcp::pdblas > DA, DB;
			pddblas_assist::split_row(this->v, m, k, DA);
			pddblas_assist::split_col(B.v, k, nc, DB);
			std::vector< kv::dd > Cv(m * nc, kv::dd(0.0));
			pddblas_assist::ozaki_acc(DA, false, DB, m, nc, k, Cv);

			c.row = m;
			c.column = nc;
			c.n = m * nc;
			c.type = ctype;
			c.v = std::move(Cv);
		}

		// C = transpose(A)*A by the Ozaki scheme
		void mulltmm(pddblas& c) const {
			if (this->type == 'S') {
				c.zeros(1, 1);
				c.v[0] = this->v[0] * this->v[0];
				return;
			}
			else if (this->type == 'C' || this->type == 'R' || this->type == 'M') {
				const int m = this->row;
				const int nc = this->column;
				std::vector< vcp::pdblas > DA;
				pddblas_assist::split_col(this->v, m, nc, DA);
				std::vector< kv::dd > Cv(nc * nc, kv::dd(0.0));
				pddblas_assist::ozaki_acc(DA, true, DA, nc, nc, m, Cv);

				c.row = nc;
				c.column = nc;
				c.n = nc * nc;
				c.type = (nc == 1) ? 'S' : 'M';
				c.v = std::move(Cv);

				// enforce exact symmetry as in pdblas
				for (int i = 0; i < nc; i++) {
					for (int j = i + 1; j < nc; j++) {
						c.v[j + nc * i] = c.v[i + nc * j];
					}
				}
			}
			else {
				vcp::throw_error<vcp::state_error>("mulltmm: unsupported matrix type: ", this->type);
			}
		}

		// Solve A*x = b : LU in double precision, then residual iteration
		// (Newton-like) with accurate residuals by the Ozaki scheme
		void linearsolve(pddblas& b, pddblas& x) {
			if (this->row != this->column || b.row != this->row) {
				vcp::throw_error<vcp::dimension_error>(
					"linearsolve: invalid dimensions: A=(", this->row, ", ", this->column,
					"), b=(", b.row, ", ", b.column, ")");
			}
			using std::fabs;
			const int an = this->row;
			const int nrhs = b.column;
			const int nb = an * nrhs;

			// LU factorization of double(A), reused in every iteration
			vcp::pdblas Ad;
			dd_to_d(*this, Ad);
			std::vector<int> ipiv(an);
			int N = an, NRHS = nrhs, lda = an, ldb = an, info;
			dgetrf_(&N, &N, &Ad.v.front(), &lda, ipiv.data(), &info);
			if (info != 0) {
				throw vcp::lapack_error("dgetrf", info, "linearsolve: dgetrf failed");
			}

			// Ozaki split of A, reused for all residuals
			std::vector< vcp::pdblas > DA;
			pddblas_assist::split_row(this->v, an, an, DA);

			x.zeros(an, nrhs);
			std::vector< kv::dd > r = b.v; // r = b - A*x with x = 0
			vcp::pdblas rd;
			rd.zeros(an, nrhs);
			char TN[] = "N";
			const int itmax = 10;
			double prevnormd = std::numeric_limits<double>::infinity();
			for (int it = 0; it < itmax; it++) {
				for (int i = 0; i < nb; i++) {
					rd.v[i] = r[i].a1;
				}
				dgetrs_(TN, &N, &NRHS, &Ad.v.front(), &lda, ipiv.data(), &rd.v.front(), &ldb, &info);
				if (info != 0) {
					throw vcp::lapack_error("dgetrs", info, "linearsolve: dgetrs failed");
				}
				double normd = 0.0, normx = 0.0;
				for (int i = 0; i < nb; i++) {
					x.v[i] += kv::dd(rd.v[i]);
					const double ad = fabs(rd.v[i]);
					if (ad > normd) normd = ad;
					const double ax = fabs(x.v[i].a1);
					if (ax > normx) normx = ax;
				}
				// converged to dd accuracy or stagnated
				if (normd == 0.0 || normd <= 0.25e-32 * normx || normd >= 0.5 * prevnormd) break;
				prevnormd = normd;

				// r = b - A*x with accurate product
				std::vector< vcp::pdblas > DX;
				pddblas_assist::split_col(x.v, an, nrhs, DX);
				std::vector< kv::dd > Ax(nb, kv::dd(0.0));
				pddblas_assist::ozaki_acc(DA, false, DX, an, nrhs, an, Ax);
				for (int i = 0; i < nb; i++) {
					r[i] = b.v[i] - Ax[i];
				}
			}
		}

		// A = inv(A) : solve A*X = I with residual iteration
		void inv() override {
			if (this->row != this->column) {
				vcp::throw_error<vcp::dimension_error>(
					"inv: matrix must be square: ", this->row, " != ", this->column);
			}
			pddblas b, x;
			b.eye(this->row);
			this->linearsolve(b, x);
			this->v = std::move(x.v);
		}

		// Cholesky factor (upper triangular) : dpotrf in double precision via the
		// pdblas policy, then refinement steps with accurate residuals
		void Cholesky() override {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("Cholesky: matrix must be symmetric");
			}
			const int an = this->row;
			const int nn = an * an;

			vcp::pdblas Ud;
			dd_to_d(*this, Ud);
			Ud.Cholesky(); // dpotrf, throws lapack_error unless positive definite

			std::vector< kv::dd > U(nn);
			for (int i = 0; i < nn; i++) {
				U[i] = kv::dd(Ud.v[i]);
			}

			char cL[] = "L", cR[] = "R", cU[] = "U", cT[] = "T", cN[] = "N";
			double one = 1.0;
			int N = an;
			const int itep = 2;
			for (int it = 0; it < itep; it++) {
				// E = A - transpose(U)*U  (accurate)
				std::vector< vcp::pdblas > DU;
				pddblas_assist::split_col(U, an, an, DU);
				std::vector< kv::dd > E(nn, kv::dd(0.0));
				pddblas_assist::ozaki_acc(DU, true, DU, an, an, an, E);
				for (int i = 0; i < nn; i++) {
					E[i] = this->v[i] - E[i];
				}

				// G = Ud^{-T} * E * Ud^{-1} in double precision (Ud = double(U))
				for (int i = 0; i < nn; i++) {
					Ud.v[i] = U[i].a1;
				}
				vcp::pdblas G;
				G.zeros(an, an);
				for (int i = 0; i < nn; i++) {
					G.v[i] = E[i].a1;
				}
				dtrsm_(cL, cU, cT, cN, &N, &N, &one, &Ud.v.front(), &N, &G.v.front(), &N);
				dtrsm_(cR, cU, cN, cN, &N, &N, &one, &Ud.v.front(), &N, &G.v.front(), &N);

				// W = triu(G,1) + diag(G)/2, then U = U + W*U
				for (int j = 0; j < an; j++) {
					for (int i = 0; i < an; i++) {
						if (i > j) {
							G.v[i + an * j] = 0.0;
						}
						else if (i == j) {
							G.v[i + an * j] *= 0.5;
						}
					}
				}
				std::vector< kv::dd > W(nn);
				for (int i = 0; i < nn; i++) {
					W[i] = kv::dd(G.v[i]);
				}
				std::vector< vcp::pdblas > DW;
				pddblas_assist::split_row(W, an, an, DW);
				std::vector< kv::dd > F(nn, kv::dd(0.0));
				pddblas_assist::ozaki_acc(DW, false, DU, an, an, an, F);
				for (int i = 0; i < nn; i++) {
					U[i] += F[i];
				}
			}

			// exact zeros below the diagonal
			for (int j = 0; j < an; j++) {
				for (int i = j + 1; i < an; i++) {
					U[i + an * j] = kv::dd(0.0);
				}
			}
			this->v = std::move(U);
		}

		// Symmetric eigenvalue problem : dsyev in double precision via the pdblas
		// policy, then Ogita-Aishima iterative refinement
		void eigsym(int itep = 4) override {
			pddblas V;
			this->eigsym(V, itep);
		}
		void eigsym(pddblas& V, int itep = 4) {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("eigsym: matrix must be symmetric");
			}
			const int an = this->row;
			const int nn = an * an;

			// initial guess by LAPACK dsyev (double)
			vcp::pdblas Ad, Vd;
			dd_to_d(*this, Ad);
			Ad.eigsym(Vd);

			std::vector< kv::dd > X(nn), lam(an);
			for (int i = 0; i < nn; i++) {
				X[i] = kv::dd(Vd.v[i]);
			}
			for (int i = 0; i < an; i++) {
				lam[i] = kv::dd(Ad.v[i + an * i]);
			}

			std::vector< vcp::pdblas > DA;
			pddblas_assist::split_row(this->v, an, an, DA);

			// at least 2 refinement steps to reach full dd accuracy
			const int kmax = (itep < 2) ? 2 : itep;
			for (int it = 0; it < kmax; it++) {
				oa_step(DA, NULL, X, lam, an);
			}

			V.zeros(an, an);
			V.v = std::move(X);
			this->zeros(an, an);
			for (int i = 0; i < an; i++) {
				this->v[i + an * i] = lam[i];
			}
		}

		// Generalized symmetric eigenvalue problem A*x = lambda*B*x : dsygv in
		// double precision via the pdblas policy, then Ogita-Aishima refinement
		// with the B-inner product
		void eigsymge(pddblas& B, int itep = 1) {
			pddblas V;
			this->eigsymge(B, V, itep);
		}
		void eigsymge(pddblas& B, pddblas& V, int itep = 1) {
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
			const int an = this->row;
			const int nn = an * an;

			// initial guess by LAPACK dsygv (double), eigenvectors with V^T B V = I
			vcp::pdblas Ad, Bd, Vd;
			dd_to_d(*this, Ad);
			dd_to_d(B, Bd);
			Ad.eigsymge(Bd, Vd);

			std::vector< kv::dd > X(nn), lam(an);
			for (int i = 0; i < nn; i++) {
				X[i] = kv::dd(Vd.v[i]);
			}
			for (int i = 0; i < an; i++) {
				lam[i] = kv::dd(Ad.v[i + an * i]);
			}

			std::vector< vcp::pdblas > DA, DB;
			pddblas_assist::split_row(this->v, an, an, DA);
			pddblas_assist::split_row(B.v, an, an, DB);

			const int kmax = (itep < 2) ? 2 : itep;
			for (int it = 0; it < kmax; it++) {
				oa_step(DA, &DB, X, lam, an);
			}

			V.zeros(an, an);
			V.v = std::move(X);
			this->zeros(an, an);
			for (int i = 0; i < an; i++) {
				this->v[i + an * i] = lam[i];
			}
		}
	};
}

#endif // VCP_PDDBLAS_HPP
