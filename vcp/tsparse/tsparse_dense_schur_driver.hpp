// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_DENSE_SCHUR_DRIVER_HPP
#define VCP_TSPARSE_DENSE_SCHUR_DRIVER_HPP

// ---------------------------------------------------------------------------
// EIG-2 Phase 3: dense nonsymmetric eigensolver driver (D-13 / D-15 kill)
//
// Pipeline (EIG-2 design §5, incl. the 2026-07-06 G-1.1 revision (3)):
//   radix-2 diagonal balancing (dgebal 'S' style; exact FP similarity)
//     -> Householder Hessenberg reduction (orthogonal accumulation Q0)
//     -> real Schur core (tsparse_real_schur.hpp; 30*m cap, honest failure)
//     -> real eigenvectors by quasi-triangular back-substitution (dtrevc-like)
//        back-transformed v = D * Q0 * Z * y and normalized.
//
// Result contract (existing dense_eigen_result<T>):
//   * eigenvalues / eigenvectors carry the REAL eigenvalues (1x1 Schur blocks)
//     with their vectors; residuals are exact lambda-space residuals against
//     the ORIGINAL input matrix (EIG-0 C-1 evidence).
//   * complex conjugate pairs (certified 2x2 blocks) cannot be represented as
//     converged pairs by the real-only result API.  They are reported through
//     complex_eigenvalues as diagnostics and force converged=false with an
//     explicit reason (honest non-convergence; a complex eigenpair API is
//     EIG-3/EIG-4 scope).  This is what keeps the strongly non-normal cases
//     (e.g. toep n=200, eigenvalue condition ~3^100) honest: any backward
//     stable method scatters such spectra into complex pairs, and pretending
//     the real subset converged would be a set-selection lie.
//   * converged=true requires: Schur success AND no complex pairs AND every
//     returned pair passing the legacy dense acceptance
//     max residual <= max(tol, tol*n*10) (the pre-existing dense scale,
//     documented per C-1).
//
// Balancing note (approval condition, G-1.1 (3)): scaling factors are exact
// powers of two (radix-2), so D^-1 A D is an EXACT floating-point similarity
// for double/dd/mpfr and an exact point scaling for kv::interval.  The
// permutation part of gebal is NOT adopted: it only short-circuits already
// isolated eigenvalues (a performance nicety), while the correctness need —
// killing scale-gradient-induced ill conditioning — is covered by scaling
// alone; omitting it keeps the back-transform a pure diagonal.  (Reported as
// the implementer decision requested by the approval record.)
//
// T-generic under SLU-GT1 P1-P6; every success-side gate is certified, and
// non-certifiable situations degrade to honest failure.
// ---------------------------------------------------------------------------

#include <cstddef>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_real_schur.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
	namespace tsparse_dense_schur {

		namespace dense_schur_detail {

			template <typename T>
			bool is_exact_zero_(const T& x) {
				typedef typename tsparse_scalar::real_type<T>::type R;
				return tsparse_scalar::abs_value(x) <= R(0);   // certainly |x| <= 0
			}

			// dgebal 'S' style radix-2 balancing: A <- D^-1 A D, d[i] = 2^k exactly.
			template <typename T>
			void balance_radix2_(std::vector<std::vector<T> >& A, std::vector<T>& d) {
				typedef typename tsparse_scalar::real_type<T>::type R;
				using tsparse_scalar::abs_value;
				const std::size_t m = A.size();
				d.assign(m, T(1));
				const R two = R(2);
				const std::size_t max_pass = 64 + 4 * m;   // safety bound; dgebal has none
				for (std::size_t pass = 0; pass < max_pass; pass++) {
					bool changed = false;
					for (std::size_t i = 0; i < m; i++) {
						R c(0), r(0);
						for (std::size_t j = 0; j < m; j++) {
							if (j == i) continue;
							c += abs_value(A[j][i]);
							r += abs_value(A[i][j]);
						}
						if (!(c > R(0)) || !(r > R(0))) continue;
						R f(1);
						int nup = 0, ndown = 0;
						// grow f while (c*f)*2 < r/f  <=>  2*c*f^2 < r
						while (nup < 512) {
							const R lhs = two * ((c * f) * f);
							if (lhs < r) { f = f * two; nup++; } else break;
						}
						// shrink f while (r/f)*2 < c*f  <=>  2*r < c*f^2
						while (ndown < 512) {
							const R rhs = (c * f) * f;
							if (two * r < rhs) { f = f / two; ndown++; } else break;
						}
						if (nup == 0 && ndown == 0) continue;
						// apply only on a certified strict gain (dgebal FACTOR=0.95)
						const R gained = c * f + r / f;
						const R factor = R(19) / R(20);
						if (!(gained < factor * (c + r))) continue;
						const T tf = T(f);
						for (std::size_t j = 0; j < m; j++) {
							A[i][j] = A[i][j] / tf;   // row i scaled by 1/f
							A[j][i] = A[j][i] * tf;   // column i scaled by f
						}
						d[i] = d[i] * tf;
						changed = true;
					}
					if (!changed) break;
				}
			}

			// Householder Hessenberg reduction with orthogonal accumulation:
			// A <- Q0^T A Q0 (upper Hessenberg, exact zeros below the subdiagonal).
			// Returns false when a nonzero column cannot be certifiably eliminated
			// (interval-only situation -> the caller fails honestly).
			template <typename T>
			bool hessenberg_reduce_(std::vector<std::vector<T> >& A,
			                        std::vector<std::vector<T> >& Q0) {
				typedef typename tsparse_scalar::real_type<T>::type R;
				using tsparse_scalar::abs_value;
				using tsparse_scalar::sqrt_value;
				const std::size_t m = A.size();
				Q0.assign(m, std::vector<T>(m, T(0)));
				for (std::size_t i = 0; i < m; i++) Q0[i][i] = T(1);
				if (m < 3) return true;
				std::vector<T> u(m, T(0));
				for (std::size_t k = 0; k + 2 < m; k++) {
					R scale(0);
					for (std::size_t i = k + 1; i < m; i++) {
						const R a = abs_value(A[i][k]);
						if (a > scale) scale = a;
					}
					if (!(scale > R(0))) {
						// column not certifiably nonzero: acceptable only if it is
						// exactly zero (already reduced); otherwise the similarity
						// cannot be certified.
						bool all_zero = true;
						for (std::size_t i = k + 2; i < m; i++)
							if (!is_exact_zero_(A[i][k])) all_zero = false;
						if (!all_zero) return false;
						continue;
					}
					R nrm2(0);
					for (std::size_t i = 0; i < m; i++) u[i] = T(0);
					for (std::size_t i = k + 1; i < m; i++) {
						u[i] = A[i][k] / T(scale);
						const R au = abs_value(u[i]);
						nrm2 += au * au;
					}
					if (!(nrm2 > R(0))) return false;
					const T nrm = T(sqrt_value(nrm2));
					const T sgn = (u[k + 1] >= T(0)) ? T(1) : T(-1);
					const T alpha = -sgn * nrm;              // scaled target subdiagonal
					u[k + 1] = u[k + 1] - alpha;
					R uu_r(0);
					for (std::size_t i = k + 1; i < m; i++) {
						const R au = abs_value(u[i]);
						uu_r += au * au;
					}
					if (!(uu_r > R(0))) return false;
					const T uu = T(uu_r);
					const T two = T(2);
					// A <- P A (rows k+1..m-1, all columns)
					for (std::size_t j = 0; j < m; j++) {
						T s(0);
						for (std::size_t i = k + 1; i < m; i++) s += u[i] * A[i][j];
						const T fct = two * s / uu;
						for (std::size_t i = k + 1; i < m; i++) A[i][j] -= fct * u[i];
					}
					// A <- A P (columns k+1..m-1, all rows)
					for (std::size_t i = 0; i < m; i++) {
						T s(0);
						for (std::size_t j = k + 1; j < m; j++) s += A[i][j] * u[j];
						const T fct = two * s / uu;
						for (std::size_t j = k + 1; j < m; j++) A[i][j] -= fct * u[j];
					}
					// Q0 <- Q0 P
					for (std::size_t i = 0; i < m; i++) {
						T s(0);
						for (std::size_t j = k + 1; j < m; j++) s += Q0[i][j] * u[j];
						const T fct = two * s / uu;
						for (std::size_t j = k + 1; j < m; j++) Q0[i][j] -= fct * u[j];
					}
					// eliminated column: analytic values assigned exactly
					A[k + 1][k] = T(scale) * alpha;
					for (std::size_t i = k + 2; i < m; i++) A[i][k] = T(0);
				}
				return true;
			}

			// Real eigenvector of the quasi-triangular Schur form for the real
			// eigenvalue at 1x1 block position j (dtrevc-like back substitution).
			template <typename T>
			bool schur_real_eigenvector_(const std::vector<std::vector<T> >& S,
			                             const std::size_t j,
			                             const T& lambda,
			                             std::vector<T>& y) {
				typedef typename tsparse_scalar::real_type<T>::type R;
				using tsparse_scalar::abs_value;
				const std::size_t m = S.size();
				const R eps = vcp::tsparse_scalar::epsilon<R>();
				y.assign(m, T(0));
				y[j] = T(1);
				if (j == 0) return true;
				std::size_t i = j;   // process rows j-1 .. 0 (1-based walk via i-1)
				while (i > 0) {
					const std::size_t q = i - 1;   // current row to determine
					const bool coupled = (q > 0) && !is_exact_zero_(S[q][q - 1]);
					if (coupled) {
						const std::size_t p = q - 1;
						T rp(0), rq(0);
						for (std::size_t t = q + 1; t <= j; t++) {
							rp += S[p][t] * y[t];
							rq += S[q][t] * y[t];
						}
						const T a11 = S[p][p] - lambda;
						const T a12 = S[p][q];
						const T a21 = S[q][p];
						const T a22 = S[q][q] - lambda;
						T det = a11 * a22 - a12 * a21;
						const R adet = abs_value(det);
						const R flo = eps * (abs_value(a11) + abs_value(a12)
						                   + abs_value(a21) + abs_value(a22));
						if (!(adet > flo)) {
							if (flo > R(0)) det = T(flo);
							else if (is_exact_zero_(rp) && is_exact_zero_(rq)) {
								i = p;   // rows solved by zero
								continue;
							}
							else return false;
						}
						y[p] = (a12 * rq - a22 * rp) / det;
						y[q] = (a21 * rp - a11 * rq) / det;
						i = p;
					} else {
						T r(0);
						for (std::size_t t = q + 1; t <= j; t++) r += S[q][t] * y[t];
						T den = S[q][q] - lambda;
						const R aden = abs_value(den);
						const R flo = eps * (abs_value(S[q][q]) + abs_value(lambda));
						if (!(aden > flo)) {
							if (flo > R(0)) den = T(flo);
							else if (is_exact_zero_(r)) { y[q] = T(0); i = q; continue; }
							else return false;
						}
						y[q] = -r / den;
						i = q;
					}
					// growth control: renormalize if the freshly computed entries
					// exceed 1/eps (double-style guard; skipped when eps == 0)
					if (eps > R(0)) {
						const R big = R(1) / eps;
						R cur = abs_value(y[i]);
						if (i + 1 <= j) {
							const R c2 = abs_value(y[i + 1]);
							if (c2 > cur) cur = c2;
						}
						if (cur > big) {
							const T inv = T(1) / T(cur);
							for (std::size_t t = i; t <= j; t++) y[t] = y[t] * inv;
						}
					}
				}
				return true;
			}

		} // namespace dense_schur_detail

		// Full-spectrum dense nonsymmetric eigensolver.  `reason` receives a
		// human-readable explanation when converged stays false for a reason
		// beyond the generic residual failure.
		template <typename T>
		vcp::tsparse_dense_linalg::dense_eigen_result<T>
		real_schur_eig_dense(const std::vector<std::vector<T> >& A_in,
		                     const typename tsparse_scalar::real_type<T>::type& tol,
		                     std::string& reason) {
			typedef typename tsparse_scalar::real_type<T>::type R;
			using tsparse_scalar::abs_value;
			namespace det = dense_schur_detail;

			vcp::tsparse_dense_linalg::dense_eigen_result<T> result;
			reason.clear();
			const std::size_t m = A_in.size();
			if (m == 0) { result.converged = true; return result; }

			// 1. balancing (exact radix-2 diagonal similarity)
			std::vector<std::vector<T> > B = A_in;
			std::vector<T> d;
			det::balance_radix2_(B, d);

			// 2. Hessenberg reduction with accumulation
			std::vector<std::vector<T> > Q0;
			if (!det::hessenberg_reduce_(B, Q0)) {
				result.converged = false;
				reason = "hessenberg reduction not certifiable for this scalar type";
				return result;
			}

			// 3. real Schur core (30*m cap, honest failure)
			vcp::tsparse_real_schur::real_schur_result<T> schur =
				vcp::tsparse_real_schur::real_schur_decompose<T>(B, true);
			result.iterations = schur.iterations;
			if (!schur.success) {
				result.converged = false;
				reason = "real Schur core: " + schur.failure_reason;
				return result;
			}

			// 4. classify blocks; complex pairs are diagnostics only (see header)
			std::size_t complex_pairs = 0;
			std::vector<std::size_t> real_positions;
			for (std::size_t p = 0; p < m; p++) {
				if (det::is_exact_zero_(schur.eig_imag[p])) real_positions.push_back(p);
				else if (p + 1 < m) { complex_pairs++; p++; }
			}
			result.complex_eigenvalues.reserve(m);
			for (std::size_t p = 0; p < m; p++) {
				result.complex_eigenvalues.push_back(
					typename vcp::tsparse_dense_linalg::dense_eigen_result<T>::eigenvalue_type(
						tsparse_scalar::real_part(schur.eig_real[p]),
						tsparse_scalar::real_part(schur.eig_imag[p])));
			}

			// 5. real eigenvectors: back substitution + v = D Q0 Z y, normalized
			bool vectors_ok = true;
			for (std::size_t idx = 0; idx < real_positions.size(); idx++) {
				const std::size_t p = real_positions[idx];
				const T lambda = schur.eig_real[p];
				std::vector<T> y;
				std::vector<T> v(m, T(0));
				bool ok = det::schur_real_eigenvector_(schur.schur_form, p, lambda, y);
				if (ok) {
					std::vector<T> w(m, T(0));
					for (std::size_t r0 = 0; r0 < m; r0++) {   // w = Z y (y is 0 above p)
						T s(0);
						for (std::size_t t = 0; t <= p; t++) s += schur.schur_vectors[r0][t] * y[t];
						w[r0] = s;
					}
					for (std::size_t r0 = 0; r0 < m; r0++) {   // v = Q0 w
						T s(0);
						for (std::size_t t = 0; t < m; t++) s += Q0[r0][t] * w[t];
						v[r0] = s;
					}
					for (std::size_t r0 = 0; r0 < m; r0++) v[r0] = d[r0] * v[r0];
					const R nv = tsparse_scalar::real_norm_value(v);
					if (nv > R(0)) {
						for (std::size_t r0 = 0; r0 < m; r0++) v[r0] = v[r0] / T(nv);
					} else {
						ok = false;
					}
				}
				if (!ok) {
					vectors_ok = false;
					v.assign(m, T(0));
				}
				result.eigenvalues.push_back(lambda);
				result.eigenvectors.push_back(v);
			}

			// 6. C-1 residuals against the ORIGINAL matrix
			R max_res(0);
			result.residuals.reserve(result.eigenvalues.size());
			for (std::size_t idx = 0; idx < result.eigenvalues.size(); idx++) {
				const std::vector<T>& v = result.eigenvectors[idx];
				const T lambda = result.eigenvalues[idx];
				std::vector<T> rvec(m, T(0));
				for (std::size_t r0 = 0; r0 < m; r0++) {
					T s(0);
					for (std::size_t t = 0; t < m; t++) s += A_in[r0][t] * v[t];
					rvec[r0] = s - lambda * v[r0];
				}
				const R rn = tsparse_scalar::real_norm_value(rvec);
				result.residuals.push_back(rn);
				if (rn > max_res) max_res = rn;
			}
			result.residual_norm = max_res;

			// 7. converged verdict (C-1; legacy dense acceptance scale
			//    max(tol, tol * n * 10) — the pre-existing qr_eig_dense convention)
			R residual_limit = tol;
			{
				const R tol_n = tol * R(static_cast<int>(m) * 10);
				if (tol_n > residual_limit) residual_limit = tol_n;
			}

			const bool residuals_pass = !(max_res > residual_limit) && vectors_ok;
			if (complex_pairs > 0) {
				result.converged = false;
				reason = "complex conjugate pairs present (" + std::to_string(complex_pairs)
				       + "): real-only dense eig_result cannot return them as converged pairs"
				         " (complex eigenpair API = EIG-3/EIG-4); honest not_converged";
			} else if (!residuals_pass) {
				result.converged = false;
				reason = vectors_ok ? "" : "eigenvector back-substitution not certifiable";
			} else {
				result.converged = true;
			}

			vcp::tsparse_dense_linalg::sort_eigenpairs(result);
			return result;
		}

	} // namespace tsparse_dense_schur
} // namespace vcp

#endif
