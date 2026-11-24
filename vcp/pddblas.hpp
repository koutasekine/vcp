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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>
#include <utility>

#include <kv/dd.hpp>

#include <vcp/pdblas.hpp>

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
		// Shift s of the splitting: a slice entry q = fl((a + sigma) - sigma),
		// sigma = 2^(e+s), |a| < 2^e, is a multiple of 2^(e+s-53) (for a < 0 the
		// sum a + sigma lies in [sigma/2, sigma)) and |q| <= 2^e. With
		// rho = ceil(log2(k)) the slice products summed over k terms remain
		// exact iff 106 + rho - 2*s <= 53.
		inline int split_shift(const int k) {
			return (53 + ceil_log2(k) + 1) / 2;
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
					for (vcp::index_t j = 0; j < k; j++) {
						const double a = fabs(T[i + m * j].a1);
						if (a > mu) mu = a;
					}
					if (mu == 0.0) continue;
					int e;
					frexp(mu, &e); // 2^(e-1) <= mu < 2^e
					const double sigma = ldexp(1.0, e + s);
					for (vcp::index_t j = 0; j < k; j++) {
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
				for (vcp::index_t j = 0; j < n; j++) {
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
			const vcp::index_t mn = static_cast<vcp::index_t>(m) * n;
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
					for (vcp::index_t idx = 0; idx < mn; idx++) {
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
			for (vcp::index_t i = 0; i < A.n; i++) {
				B.v[i] = A.v[i].a1;
			}
		}

		// Frobenius norm (leading parts) of a dd array. The entries are divided
		// by the largest absolute value before they are squared, so that the
		// sum of squares neither overflows nor underflows.
		static double fnorm(const std::vector< kv::dd >& A) {
			using std::fabs;
			double amax = 0.0;
			for (std::size_t i = 0; i < A.size(); i++) {
				const double a = fabs(A[i].a1);
				if (a > amax) amax = a;
			}
			if (amax == 0.0) return 0.0;
			double s = 0.0;
			for (std::size_t i = 0; i < A.size(); i++) {
				const double t = A[i].a1 / amax;
				s += t * t;
			}
			return amax * std::sqrt(s);
		}

		// Indices 0..l-1 ordered so that lam is ascending
		static std::vector< int > ascending_index(const std::vector< kv::dd >& lam) {
			std::vector< int > idx(lam.size());
			for (std::size_t i = 0; i < lam.size(); i++) {
				idx[i] = static_cast<int>(i);
			}
			const kv::dd* w = lam.data();
			std::stable_sort(idx.begin(), idx.end(), [w](const int a, const int b) { return w[a] < w[b]; });
			return idx;
		}

		// Algorithm 1 (RefSyEv) of T. Ogita, K. Aishima, "Iterative refinement
		// for symmetric eigenvalue decomposition II: clustered eigenvalues",
		// Japan J. Indust. Appl. Math. 36 (2019), 435-459, for the standard or
		// generalized problem. DA / DB are row-wise Ozaki splits of A and B
		// (DB = NULL for the standard problem B = I), anorm = ||A||_F. X (n x l,
		// columns are approximate eigenvectors with X^T B X ~ I) is replaced by
		// X + X*E. lam receives the approximate eigenvalues, enorm = ||E||_F,
		// omega the criterion of line 6 (Frobenius norms).
		static void oa_step(const std::vector< vcp::pdblas >& DA,
			const std::vector< vcp::pdblas >* DB, const double anorm,
			std::vector< kv::dd >& X, std::vector< kv::dd >& lam, const int n, const int l,
			double& enorm, double& omega) {
			using std::fabs;
			const kv::dd dd0 = kv::dd(0.0);
			const int nl = n * l;
			const int ll = l * l;

			std::vector< vcp::pdblas > DX;
			pddblas_assist::split_col(X, n, l, DX);

			// AX = A*X,  S = X^T * AX
			std::vector< kv::dd > AX(nl, dd0);
			pddblas_assist::ozaki_acc(DA, false, DX, n, l, n, AX);
			std::vector< vcp::pdblas > DAX;
			pddblas_assist::split_col(AX, n, l, DAX);
			std::vector< kv::dd > S(ll, dd0);
			pddblas_assist::ozaki_acc(DX, true, DAX, l, l, n, S);

			// G = X^T * B * X
			std::vector< kv::dd > G(ll, dd0);
			if (DB != NULL) {
				std::vector< kv::dd > BX(nl, dd0);
				pddblas_assist::ozaki_acc(*DB, false, DX, n, l, n, BX);
				std::vector< vcp::pdblas > DBX;
				pddblas_assist::split_col(BX, n, l, DBX);
				pddblas_assist::ozaki_acc(DX, true, DBX, l, l, n, G);
			}
			else {
				pddblas_assist::ozaki_acc(DX, true, DX, l, l, n, G);
			}

			// lambda_i = S_ii / (X^T B X)_ii,   R = I - G
			lam.assign(l, dd0);
			for (int i = 0; i < l; i++) {
				lam[i] = S[i + l * i] / G[i + l * i];
			}
			std::vector< kv::dd >& R = G;
			for (int idx = 0; idx < ll; idx++) {
				R[idx] = -R[idx];
			}
			for (int i = 0; i < l; i++) {
				R[i + l * i] += kv::dd(1.0);
			}

			// omega = 2 (||S - D||_F + ||A||_F ||R||_F)
			std::vector< kv::dd > SD = S;
			for (int i = 0; i < l; i++) {
				SD[i + l * i] -= lam[i];
			}
			omega = 2.0 * (fnorm(SD) + anorm * fnorm(R));

			// E_ij = (S_ij + lambda_j R_ij)/(lambda_j - lambda_i)  (separated case)
			// E_ij = R_ij/2                                        (otherwise)
			std::vector< kv::dd > E(ll, dd0);
			const kv::dd half = kv::dd(0.5);
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < l; j++) {
				for (int i = 0; i < l; i++) {
					if (i == j) {
						E[i + l * i] = R[i + l * i] * half;
					}
					else {
						kv::dd d = lam[j] - lam[i];
						if (fabs(d.a1) + fabs(d.a2) > omega) {
							E[i + l * j] = (S[i + l * j] + lam[j] * R[i + l * j]) / d;
						}
						else {
							E[i + l * j] = R[i + l * j] * half;
						}
					}
				}
			}
			enorm = fnorm(E);

			// X = X + X*E
			std::vector< vcp::pdblas > DE, DXr;
			pddblas_assist::split_col(E, l, l, DE);
			pddblas_assist::split_row(X, n, l, DXr);
			std::vector< kv::dd > XE(nl, dd0);
			pddblas_assist::ozaki_acc(DXr, false, DE, n, l, l, XE);
			for (int idx = 0; idx < nl; idx++) {
				X[idx] += XE[idx];
			}
		}

		// Algorithm 2 (RefSyEvCL) of the same paper, one call. Av / Bv are the
		// matrices (Bv = NULL for the standard problem), DA / DB their row-wise
		// splits. X (n x n) is replaced by X'. enorm receives ||E||_F of line 2.
		// Returns false when the algorithm returned at line 3 or line 15
		// (improvement cannot be expected); X then holds the X' of that moment.
		static bool oa_cluster_step(const std::vector< kv::dd >& Av, const std::vector< vcp::pdblas >& DA,
			const std::vector< kv::dd >* Bv, const std::vector< vcp::pdblas >* DB, const double anorm,
			std::vector< kv::dd >& X, std::vector< kv::dd >& lam, const int n, double& enorm) {
			using std::fabs;
			const kv::dd dd0 = kv::dd(0.0);
			double omega = 0.0;
			oa_step(DA, DB, anorm, X, lam, n, n, enorm, omega);   // line 2
			if (!(enorm < 1.0)) return false;                     // line 3
			// line 4: the index sets J_k of (27) are the runs of the ascending
			// eigenvalues whose neighbouring differences are at most omega
			const std::vector< int > idx = ascending_index(lam);
			int p = 0;
			while (p < n) {
				int q = p;
				while (q + 1 < n) {
					const kv::dd d = lam[idx[q + 1]] - lam[idx[q]];
					if (!(fabs(d.a1) + fabs(d.a2) <= omega)) break;
					q++;
				}
				const int nk = q - p + 1;
				if (nk >= 2) {
					// line 6: V_k = X'(:, J_k)
					std::vector< kv::dd > V(n * nk, dd0);
					for (int c = 0; c < nk; c++) {
						for (int i = 0; i < n; i++) {
							V[i + n * c] = X[i + n * idx[p + c]];
						}
					}
					// line 7: mu_k, line 8: A_k = A - mu_k I (Remark 3: A - mu_k B)
					const kv::dd mu = (lam[idx[p]] + lam[idx[q]]) * kv::dd(0.5);
					std::vector< kv::dd > Ak = Av;
					if (Bv != NULL) {
						for (int x = 0; x < n * n; x++) {
							Ak[x] -= mu * (*Bv)[x];
						}
					}
					else {
						for (int i = 0; i < n; i++) {
							Ak[i + n * i] -= mu;
						}
					}
					std::vector< vcp::pdblas > DAk;
					pddblas_assist::split_row(Ak, n, n, DAk);
					const double anormk = fnorm(Ak);
					// line 9: C_k = V_k^T A_k V_k, line 10: T_k = fl(C_k), line 11: W_k = eig(T_k)
					std::vector< vcp::pdblas > DV;
					pddblas_assist::split_col(V, n, nk, DV);
					std::vector< kv::dd > AV(n * nk, dd0);
					pddblas_assist::ozaki_acc(DAk, false, DV, n, nk, n, AV);
					std::vector< vcp::pdblas > DAV;
					pddblas_assist::split_col(AV, n, nk, DAV);
					std::vector< kv::dd > C(nk * nk, dd0);
					pddblas_assist::ozaki_acc(DV, true, DAV, nk, nk, n, C);
					vcp::pdblas Td, Wd;
					Td.zeros(nk, nk);
					for (int j = 0; j < nk; j++) {
						for (int i = 0; i <= j; i++) {
							Td.v[i + nk * j] = C[i + nk * j].a1;
							Td.v[j + nk * i] = C[i + nk * j].a1;
						}
					}
					Td.eigsym(Wd);
					// line 12: V_k = V_k W_k
					std::vector< kv::dd > W(nk * nk);
					for (int i = 0; i < nk * nk; i++) {
						W[i] = kv::dd(Wd.v[i]);
					}
					std::vector< vcp::pdblas > DVr, DW;
					pddblas_assist::split_row(V, n, nk, DVr);
					pddblas_assist::split_col(W, nk, nk, DW);
					std::vector< kv::dd > VW(n * nk, dd0);
					pddblas_assist::ozaki_acc(DVr, false, DW, n, nk, nk, VW);
					V.swap(VW);
					// lines 13-17. Unlike the wording of Remark 3, B is also used here (R = I - V_k^T B V_k).
					// In addition to the condition of line 17 the loop ends when
					// ||E_k|| no longer decreases, and after 10 passes.
					std::vector< kv::dd > lamk;
					double ek = 0.0, ekprev = 0.0, omegak = 0.0;
					for (int pass = 0; pass < 10; pass++) {
						oa_step(DAk, DB, anormk, V, lamk, n, nk, ek, omegak);   // line 14
						if (!(ek < 1.0)) return false;                        // line 15
						if (ek <= enorm) break;                               // line 17
						if (pass > 0 && !(ek < ekprev)) break;
						ekprev = ek;
					}
					// line 18
					for (int c = 0; c < nk; c++) {
						for (int i = 0; i < n; i++) {
							X[i + n * idx[p + c]] = V[i + n * c];
						}
					}
				}
				p = q + 1;
			}
			return true;
		}

		// Refinement of the eigenvectors X (n x n) of A x = lambda x (Bv = NULL)
		// or A x = lambda B x: Algorithm 2 at most max(6, itep) times (it ends
		// earlier when ||E|| < 2^-52 after an update, or when Algorithm 2
		// returned without improvement), orthonormalization steps, Rayleigh
		// quotients, ascending order.
		static void oa_refine(const std::vector< kv::dd >& Av, const std::vector< kv::dd >* Bv,
			std::vector< kv::dd >& X, std::vector< kv::dd >& lam, const int n, const int itep) {
			const kv::dd dd0 = kv::dd(0.0);
			const int nn = n * n;
			std::vector< vcp::pdblas > DA, DBs;
			pddblas_assist::split_row(Av, n, n, DA);
			const std::vector< vcp::pdblas >* DB = NULL;
			if (Bv != NULL) {
				pddblas_assist::split_row(*Bv, n, n, DBs);
				DB = &DBs;
			}
			const double anorm = fnorm(Av);

			const int kmax = (itep < 6) ? 6 : itep;
			std::vector< kv::dd > Xkeep;
			for (int it = 0; it < kmax; it++) {
				double enorm = 0.0;
				const bool go_on = oa_cluster_step(Av, DA, Bv, DB, anorm, X, lam, n, enorm);
				if (!go_on) break;
				// the next correction would be below the unit of kv::dd
				if (enorm < std::ldexp(1.0, -52)) break;
			}

			// Steps X = X (I + R / 2), R = I - X^T B X: formula (25) of the paper used for all pairs.
			// Repeated (at most 4 times) until the R that was used satisfies ||R||_F <= 2^-45.
			for (int step = 0; step < 4; step++) {
				std::vector< vcp::pdblas > DX;
				pddblas_assist::split_col(X, n, n, DX);
				std::vector< kv::dd > R(nn, dd0);
				if (DB != NULL) {
					std::vector< kv::dd > BX(nn, dd0);
					pddblas_assist::ozaki_acc(*DB, false, DX, n, n, n, BX);
					std::vector< vcp::pdblas > DBX;
					pddblas_assist::split_col(BX, n, n, DBX);
					pddblas_assist::ozaki_acc(DX, true, DBX, n, n, n, R);
				}
				else {
					pddblas_assist::ozaki_acc(DX, true, DX, n, n, n, R);
				}
				const kv::dd half = kv::dd(0.5);
				for (int idx = 0; idx < nn; idx++) {
					R[idx] = -R[idx];
				}
				for (int i = 0; i < n; i++) {
					R[i + n * i] += kv::dd(1.0);
				}
				const double rnorm = fnorm(R);
				for (int idx = 0; idx < nn; idx++) {
					R[idx] = R[idx] * half;
				}
				std::vector< vcp::pdblas > DR, DXr;
				pddblas_assist::split_col(R, n, n, DR);
				pddblas_assist::split_row(X, n, n, DXr);
				std::vector< kv::dd > XR(nn, dd0);
				pddblas_assist::ozaki_acc(DXr, false, DR, n, n, n, XR);
				for (int idx = 0; idx < nn; idx++) {
					X[idx] += XR[idx];
				}
				if (rnorm <= std::ldexp(1.0, -45)) break;
			}

			// Rayleigh quotients of the returned X
			{
				std::vector< vcp::pdblas > DX;
				pddblas_assist::split_col(X, n, n, DX);
				std::vector< kv::dd > AX(nn, dd0), BX;
				pddblas_assist::ozaki_acc(DA, false, DX, n, n, n, AX);
				if (DB != NULL) {
					BX.assign(nn, dd0);
					pddblas_assist::ozaki_acc(*DB, false, DX, n, n, n, BX);
				}
				lam.assign(n, dd0);
				for (int j = 0; j < n; j++) {
					kv::dd num = dd0, den = dd0;
					for (int i = 0; i < n; i++) {
						const kv::dd x = X[i + n * j];
						num += x * AX[i + n * j];
						den += x * ((DB != NULL) ? BX[i + n * j] : x);
					}
					lam[j] = num / den;
				}
			}

			// ascending order
			const std::vector< int > idx = ascending_index(lam);
			Xkeep = X;
			const std::vector< kv::dd > lkeep = lam;
			for (int c = 0; c < n; c++) {
				lam[c] = lkeep[idx[c]];
				for (int i = 0; i < n; i++) {
					X[i + n * c] = Xkeep[i + n * idx[c]];
				}
			}
		}

	public:
		// C = A*B by the Ozaki scheme
		void mulmm(const pddblas& B, pddblas& c) const {
			if (this->type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c = B;
				for (vcp::index_t i = 0; i < B.n; i++) {
					c.v[i] *= this->v[0];
				}
				return;
			}
			else if ((this->type == 'C' || this->type == 'R' || this->type == 'M') && B.type == 'S') {
				c = *this;
				for (vcp::index_t i = 0; i < this->n; i++) {
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

			const int m = static_cast<int>(this->row);
			const int k = static_cast<int>(this->column);
			const int nc = static_cast<int>(B.column);
			std::vector< vcp::pdblas > DA, DB;
			pddblas_assist::split_row(this->v, m, k, DA);
			pddblas_assist::split_col(B.v, k, nc, DB);
			std::vector< kv::dd > Cv(static_cast<std::size_t>(m) * nc, kv::dd(0.0));
			pddblas_assist::ozaki_acc(DA, false, DB, m, nc, k, Cv);

			c.row = m;
			c.column = nc;
			c.n = static_cast<vcp::index_t>(m) * nc;
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
				const int m = static_cast<int>(this->row);
				const int nc = static_cast<int>(this->column);
				std::vector< vcp::pdblas > DA;
				pddblas_assist::split_col(this->v, m, nc, DA);
				std::vector< kv::dd > Cv(static_cast<std::size_t>(nc) * nc, kv::dd(0.0));
				pddblas_assist::ozaki_acc(DA, true, DA, nc, nc, m, Cv);

				c.row = nc;
				c.column = nc;
				c.n = static_cast<vcp::index_t>(nc) * nc;
				c.type = (nc == 1) ? 'S' : 'M';
				c.v = std::move(Cv);

				// enforce exact symmetry as in pdblas
				for (vcp::index_t i = 0; i < nc; i++) {
					for (vcp::index_t j = i + 1; j < nc; j++) {
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
			const int an = static_cast<int>(this->row);
			const int nrhs = static_cast<int>(b.column);
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
			b.eye(static_cast<int>(this->row));
			this->linearsolve(b, x);
			this->v = std::move(x.v);
		}

		// Cholesky factor (upper triangular) : dpotrf in double precision via the
		// pdblas policy, then refinement steps with accurate residuals
		void Cholesky() override {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("Cholesky: matrix must be symmetric");
			}
			const int an = static_cast<int>(this->row);
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
		// policy, then Ogita-Aishima iterative refinement (with the treatment of
		// clustered eigenvalues); itep is an upper bound of the number of
		// refinement steps (at least 6 are allowed)
		void eigsym(int itep = 4) override {
			pddblas V;
			this->eigsym(V, itep);
		}
		void eigsym(pddblas& V, int itep = 4) {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("eigsym: matrix must be symmetric");
			}
			const int an = static_cast<int>(this->row);
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

			oa_refine(this->v, NULL, X, lam, an, itep);

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
			const int an = static_cast<int>(this->row);
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

			oa_refine(this->v, &B.v, X, lam, an, itep);

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
