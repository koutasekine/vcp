// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_REAL_SCHUR_HPP
#define VCP_TSPARSE_REAL_SCHUR_HPP

// ---------------------------------------------------------------------------
// EIG-2: real Schur QR core (LAPACK dhseqr/dlahqr style, written from scratch)
//
// Root-cause replacement for the D-QR defect family (frozen no-op steps,
// bottom-2-only deflation, cap burning, denormal persistence — see
// sandbox/docs/issues/eig_solver_defect_audit.md) and the D-13 unresponsive
// dense nonsymmetric path.  Spec: sandbox/docs/design/EIG-2_design.md §4
// (including the 2026-07-06 G-1.1 revisions).  This file deliberately shares
// no code with the legacy core (tsparse_eigensolvers.hpp) — B-17.
//
// Contract summary:
//  * Input: upper Hessenberg H (m x m, copied by value; entries below the
//    first subdiagonal are treated as zero).  Optional Schur vector
//    accumulation Z with H0 = Z T Z^T.
//  * Deflation (§4.1): every sweep scans ALL subdiagonals of the active
//    region with the scale-adaptive LAPACK criterion
//        certainly( |h[k][k-1]| <= eps * (|h[k-1][k-1]| + |h[k][k]|) )
//    (neighbor-augmented when the diagonal scale is exactly zero), then
//    assigns an EXACT T(0) — the algebraic kill of denormal persistence
//    (D-8/K-1) and of the tol^2 small-scale false deflation (D-15).
//    No fixed absolute threshold exists anywhere.
//  * Certified-only gates (SLU-GT1 P1, GT1.1 addendum): deflation and the
//    2x2 real/complex discriminant are success-side gates and use certified
//    comparisons with an explicit third branch (uncertifiable -> honest
//    failure).  For double the third branch is unreachable.
//  * Shifts (§4.3): Francis implicit double shift (Wilkinson, trailing 2x2);
//    ad-hoc exceptional shift (LAPACK dat1/dat2 style, built
//    deterministically from subdiagonal magnitudes) after 10 and 20 sweeps
//    without progress on the same bottom eigenvalue; afterwards plain shifts
//    until the iteration cap.  A negligible initial shift column escalates
//    to the exceptional shift instead of silently no-op-ing (§4.4 — the
//    direct prohibition of the D-QR freeze mechanism).
//  * Standardization (§4.2, G-1.1 revision (5)): every 2x2 block of the
//    final Schur form is a certified complex pair; certified real pairs are
//    split into two 1x1 blocks by a deterministic rotation (dlanv2-like).
//  * Budget (§4.6, B-20): total sweep cap = 30*m, never re-initialized or
//    bypassed.  Cap reached with work remaining = honest failure carrying
//    partial results: trailing eigenvalues [unconverged_count, m) are valid,
//    the leading block [0, unconverged_count) is undetermined (values R(0),
//    flag semantics — GT1 P4).
//  * Determinism (§4.7): no randomness; identical input -> identical output.
//  * T-generic (GT1 P1-P6): module-scalar operations only (ADL abs/sqrt via
//    tsparse_scalar, no numeric_limits direct calls, no double casts, no
//    SWO-dependent std algorithms on T).
// ---------------------------------------------------------------------------

#include <cstddef>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
	namespace tsparse_real_schur {

		template <typename T>
		struct real_schur_result {
			std::vector<std::vector<T> > schur_form;     // quasi-upper-triangular T (m x m)
			std::vector<std::vector<T> > schur_vectors;  // Z (m x m) if requested, else empty
			std::vector<T> eig_real;                     // length m; valid on [unconverged_count, m)
			std::vector<T> eig_imag;                     // complex pairs adjacent: +s at i, -s at i+1
			bool success;                                // all m eigenvalues certified
			std::size_t iterations;                      // consumed QR sweeps (<= 30*m)
			std::size_t unconverged_count;               // leading rows [0, uc) undetermined
			std::string failure_reason;                  // non-empty iff !success

			real_schur_result() : success(false), iterations(0), unconverged_count(0) {}
		};

		namespace real_schur_detail {

			// Apply the Householder similarity for reflector u (rows/cols k..k+len-1)
			// H := P H P with P = I - (2/uu) u u^T, restricted to the structurally
			// nonzero ranges (full trailing columns / leading rows: Schur form and Z
			// are maintained for the whole matrix, not only the active block).
			template <typename T>
			void apply_householder_full(std::vector<std::vector<T> >& H,
			                            std::vector<std::vector<T> >* Z,
			                            const std::size_t m,
			                            const std::size_t k,
			                            const std::size_t len,
			                            const std::size_t col_start,
			                            const std::size_t row_end,      // inclusive
			                            const T& u0, const T& u1, const T& u2,
			                            const T& uu) {
				const T two = T(2);
				// left: rows k..k+len-1, columns col_start..m-1
				for (std::size_t j = col_start; j < m; j++) {
					T d = u0 * H[k][j] + u1 * H[k + 1][j];
					if (len == 3) d += u2 * H[k + 2][j];
					const T f = two * d / uu;
					H[k][j] -= f * u0;
					H[k + 1][j] -= f * u1;
					if (len == 3) H[k + 2][j] -= f * u2;
				}
				// right: columns k..k+len-1, rows 0..row_end
				for (std::size_t r = 0; r <= row_end; r++) {
					T d = H[r][k] * u0 + H[r][k + 1] * u1;
					if (len == 3) d += H[r][k + 2] * u2;
					const T f = two * d / uu;
					H[r][k] -= f * u0;
					H[r][k + 1] -= f * u1;
					if (len == 3) H[r][k + 2] -= f * u2;
				}
				if (Z != 0) {
					std::vector<std::vector<T> >& Zr = *Z;
					for (std::size_t r = 0; r < m; r++) {
						T d = Zr[r][k] * u0 + Zr[r][k + 1] * u1;
						if (len == 3) d += Zr[r][k + 2] * u2;
						const T f = two * d / uu;
						Zr[r][k] -= f * u0;
						Zr[r][k + 1] -= f * u1;
						if (len == 3) Zr[r][k + 2] -= f * u2;
					}
				}
			}

			// Apply the Givens similarity G = [[c, -s], [s, c]] on index pair
			// (p, p+1): H := G^T H G, Z := Z G.
			template <typename T>
			void apply_givens_full(std::vector<std::vector<T> >& H,
			                       std::vector<std::vector<T> >* Z,
			                       const std::size_t m,
			                       const std::size_t p,
			                       const T& c, const T& s) {
				// left (G^T): rows p, p+1, columns p..m-1 (columns < p are exact zeros)
				for (std::size_t j = p; j < m; j++) {
					const T hp = H[p][j];
					const T hq = H[p + 1][j];
					H[p][j] = c * hp + s * hq;
					H[p + 1][j] = c * hq - s * hp;
				}
				// right (G): columns p, p+1, rows 0..p+1 (rows > p+1 are exact zeros)
				for (std::size_t r = 0; r <= p + 1; r++) {
					const T hp = H[r][p];
					const T hq = H[r][p + 1];
					H[r][p] = c * hp + s * hq;
					H[r][p + 1] = c * hq - s * hp;
				}
				if (Z != 0) {
					std::vector<std::vector<T> >& Zr = *Z;
					for (std::size_t r = 0; r < m; r++) {
						const T zp = Zr[r][p];
						const T zq = Zr[r][p + 1];
						Zr[r][p] = c * zp + s * zq;
						Zr[r][p + 1] = c * zq - s * zp;
					}
				}
			}

		} // namespace real_schur_detail

		template <typename T>
		real_schur_result<T> real_schur_decompose(const std::vector<std::vector<T> >& H_in,
		                                          const bool accumulate_schur_vectors) {
			typedef typename tsparse_scalar::real_type<T>::type R;
			using vcp::tsparse_scalar::abs_value;
			using vcp::tsparse_scalar::sqrt_value;

			real_schur_result<T> res;
			const std::size_t m = H_in.size();
			for (std::size_t i = 0; i < m; i++) {
				if (H_in[i].size() != m) {
					vcp::throw_error<vcp::dimension_error>(
						"tsparse_real_schur::real_schur_decompose: input matrix is not square");
				}
			}

			res.eig_real.assign(m, T(0));
			res.eig_imag.assign(m, T(0));
			if (m == 0) {
				res.success = true;
				return res;
			}

			// working copy: keep the Hessenberg envelope, force exact zeros below
			std::vector<std::vector<T> >& H = res.schur_form;
			H.assign(m, std::vector<T>(m, T(0)));
			for (std::size_t i = 0; i < m; i++)
				for (std::size_t j = (i > 0 ? i - 1 : 0); j < m; j++) H[i][j] = H_in[i][j];

			std::vector<std::vector<T> >* Zp = 0;
			if (accumulate_schur_vectors) {
				res.schur_vectors.assign(m, std::vector<T>(m, T(0)));
				for (std::size_t i = 0; i < m; i++) res.schur_vectors[i][i] = T(1);
				Zp = &res.schur_vectors;
			}

			// subdiag_zero[k] == 1  <=>  H[k][k-1] holds an exact structural zero.
			// Block navigation uses only this bookkeeping (never a straddling
			// comparison), so an interval subdiagonal that merely contains zero is
			// treated as nonzero and can only leave through certified deflation.
			std::vector<char> subdiag_zero(m, 0);
			for (std::size_t k = 1; k < m; k++) {
				if (abs_value(H[k][k - 1]) <= R(0)) {   // certainly |h| <= 0  <=>  h == 0 exactly
					H[k][k - 1] = T(0);
					subdiag_zero[k] = 1;
				}
			}

			const R eps = vcp::tsparse_scalar::epsilon<R>();
			const T dat1 = T(3) / T(4);      // LAPACK dlahqr exceptional shift constants
			const T dat2 = T(-7) / T(16);    // (exact in binary floating point)
			const std::size_t cap = 30 * m;  // §4.6 / B-20: the only wall-clock levee

			std::size_t bottom = m;          // rows [bottom, m) are finalized
			std::size_t its = 0;             // sweeps without progress on current bottom

			// finalize helpers operate on the block ending at row i (inclusive)
			while (bottom > 0) {
				const std::size_t i = bottom - 1;

				// ---- deflation pass (§4.1): all subdiagonals of the active region ----
				for (std::size_t k = 1; k <= i; k++) {
					if (subdiag_zero[k]) continue;
					R tst = abs_value(H[k - 1][k - 1]) + abs_value(H[k][k]);
					if (!(tst > R(0))) {
						// diagonal scale not certifiably positive: augment with the
						// neighboring pattern (LAPACK dlahqr protection floor built
						// from the matrix itself; no absolute constant)
						if (k >= 2 && !subdiag_zero[k - 1]) tst += abs_value(H[k - 1][k - 2]);
						if (k + 1 <= i && !subdiag_zero[k + 1]) tst += abs_value(H[k + 1][k]);
					}
					if (abs_value(H[k][k - 1]) <= eps * tst) {   // certified deflation only
						H[k][k - 1] = T(0);                       // exact zero assignment
						subdiag_zero[k] = 1;
					}
				}

				// ---- bottom active block [l..i] ----
				std::size_t l = i;
				while (l > 0 && !subdiag_zero[l]) l--;

				if (l == i) {
					// 1x1 block: eigenvalue confirmed
					res.eig_real[i] = H[i][i];
					res.eig_imag[i] = T(0);
					bottom = i;
					its = 0;
					continue;
				}

				if (l + 1 == i) {
					// 2x2 block: certified three-branch (GT1.1 form) + standardization.
					// All products are formed on entries normalized by the block scale
					// so that squares neither underflow (~1e-300 inputs) nor overflow;
					// the discriminant sign and the rotation are scale-invariant.
					R bs = abs_value(H[l][l]);
					{
						const R t1 = abs_value(H[l][i]);
						const R t2 = abs_value(H[i][l]);
						const R t3 = abs_value(H[i][i]);
						if (t1 > bs) bs = t1;
						if (t2 > bs) bs = t2;
						if (t3 > bs) bs = t3;
					}
					if (!(bs > R(0))) {
						// cannot certify a positive block scale: for double this means
						// the block is exactly zero except the (nonzero) subdiagonal —
						// impossible since bs >= |H[i][l]|; for interval it means the
						// certification fails -> honest failure.
						res.success = false;
						res.unconverged_count = i + 1;
						res.failure_reason =
							"2x2 block: scale not certifiably positive";
						return res;
					}
					const T a = H[l][l] / T(bs);
					const T b = H[l][i] / T(bs);
					const T c = H[i][l] / T(bs);
					const T d = H[i][i] / T(bs);
					const T tr = a + d;
					const T amd = a - d;
					const T disc = amd * amd + T(4) * (b * c);
					if (disc >= T(0)) {
						// certified real pair -> split into 1x1 x 2 by a deterministic
						// rotation (§4.2 standardization; dlanv2-like via eigenvector)
						const T sd = sqrt_value(disc);
						const T lam1 = (amd >= T(0)) ? (tr + sd) / T(2) : (tr - sd) / T(2);
						// eigenvector of scaled [[a,b],[c,d]] for lam1: v = (lam1-d, c);
						// the rotation is invariant under the block scaling
						const T v0 = lam1 - d;
						const R av0 = abs_value(v0);
						const R avc = abs_value(c);
						R vs = av0;
						if (avc > vs) vs = avc;
						if (!(vs > R(0))) {
							res.success = false;
							res.unconverged_count = i + 1;
							res.failure_reason =
								"real 2x2 block: rotation not certifiable (zero eigenvector norm)";
							return res;
						}
						const T v0s = v0 / T(vs);
						const T cs = c / T(vs);
						const R av0s = abs_value(v0s);
						const R acs = abs_value(cs);
						const R vv = av0s * av0s + acs * acs;
						if (!(vv > R(0))) {
							res.success = false;
							res.unconverged_count = i + 1;
							res.failure_reason =
								"real 2x2 block: rotation not certifiable (zero eigenvector norm)";
							return res;
						}
						const T vn = T(sqrt_value(vv));
						const T gc = v0s / vn;
						const T gs = cs / vn;
						real_schur_detail::apply_givens_full(H, Zp, m, l, gc, gs);
						H[i][l] = T(0);                  // annihilated analytically -> exact zero
						subdiag_zero[i] = 1;
						res.eig_real[l] = H[l][l];
						res.eig_imag[l] = T(0);
						res.eig_real[i] = H[i][i];
						res.eig_imag[i] = T(0);
					} else if (disc < T(0)) {
						// certified complex pair: block stays 2x2
						const T s = sqrt_value(-disc);
						res.eig_real[l] = (H[l][l] + H[i][i]) / T(2);
						res.eig_real[i] = res.eig_real[l];
						res.eig_imag[l] = T(bs) * s / T(2);
						res.eig_imag[i] = -res.eig_imag[l];
					} else {
						// 0-straddling discriminant: certification impossible (interval);
						// honest failure — unreachable for double (total order).
						res.success = false;
						res.unconverged_count = i + 1;
						res.failure_reason =
							"2x2 block: discriminant sign not certifiable (0-straddle)";
						return res;
					}
					bottom = l;
					its = 0;
					continue;
				}

				// ---- QR sweep needed: budget check first (§4.6, B-20) ----
				if (res.iterations >= cap) {
					res.success = false;
					res.unconverged_count = i + 1;
					res.failure_reason = "iteration cap 30*m reached (honest failure)";
					return res;
				}

				// ---- shift selection (§4.3) ----
				// The shift column is formed on entries normalized by a local scale
				// (LAPACK dlahqr style) so that the products neither underflow at
				// ~1e-300 input scale nor overflow at huge scale.  Only the DIRECTION
				// of (x, y, z) matters for the implicit-shift reflector.
				R lsc = abs_value(H[l][l]);
				{
					const R t1 = abs_value(H[l][l + 1]);
					const R t2 = abs_value(H[l + 1][l]);
					const R t3 = abs_value(H[l + 1][l + 1]);
					const R t4 = abs_value(H[i - 1][i - 1]);
					const R t5 = abs_value(H[i][i]);
					const R t6 = abs_value(H[i - 1][i]);
					const R t7 = abs_value(H[i][i - 1]);
					const R t8 = (i >= l + 2) ? abs_value(H[l + 2][l + 1]) : R(0);
					if (t1 > lsc) lsc = t1;
					if (t2 > lsc) lsc = t2;
					if (t3 > lsc) lsc = t3;
					if (t4 > lsc) lsc = t4;
					if (t5 > lsc) lsc = t5;
					if (t6 > lsc) lsc = t6;
					if (t7 > lsc) lsc = t7;
					if (t8 > lsc) lsc = t8;
				}
				bool exceptional = (its == 10 || its == 20);
				for (int attempt = 0; attempt < 2; attempt++) {
					T x0 = T(0), y0 = T(0), z0 = T(0);
					R amax = R(0);
					if (lsc > R(0)) {
						// scaled entries (all O(1) in the certified-positive-scale case)
						const T al = H[l][l] / T(lsc);
						const T bl = H[l][l + 1] / T(lsc);
						const T cl = H[l + 1][l] / T(lsc);
						const T dl = H[l + 1][l + 1] / T(lsc);
						const T el = (i >= l + 2) ? H[l + 2][l + 1] / T(lsc) : T(0);
						T s1s, s2s;   // s1/lsc and s2/lsc^2
						if (exceptional) {
							const T sms = (abs_value(H[i][i - 1])
							            + ((i >= l + 2) ? abs_value(H[i - 1][i - 2]) : T(0)))
							            / T(lsc);
							const T aas = dat1 * sms + H[i][i] / T(lsc);
							s1s = aas + aas;
							s2s = aas * aas - (dat2 * sms) * sms;
						} else {
							const T p1 = H[i - 1][i - 1] / T(lsc);
							const T p2 = H[i][i] / T(lsc);
							const T q1 = H[i - 1][i] / T(lsc);
							const T q2 = H[i][i - 1] / T(lsc);
							s1s = p1 + p2;
							s2s = p1 * p2 - q1 * q2;
						}
						x0 = al * al + bl * cl - s1s * al + s2s;
						y0 = cl * (al + dl - s1s);
						z0 = cl * el;
						amax = abs_value(x0);
						const R ay = abs_value(y0);
						const R az = abs_value(z0);
						if (ay > amax) amax = ay;
						if (az > amax) amax = az;
					}
					if (!(amax > R(0))) {
						// negligible shift vector: never a silent no-op (§4.4).
						// Escalate to the exceptional shift; if that is negligible
						// too, the sweep below degenerates but still consumes budget
						// (honest progress accounting -> the cap terminates).
						if (!exceptional) { exceptional = true; continue; }
						x0 = T(0); y0 = T(0); z0 = T(0);
					}
					// perform the bulge chase with this (scaled) shift column
					T x = x0, y = y0, z = z0;
					for (std::size_t k = l; k < i; k++) {
						const std::size_t len = (k + 2 <= i) ? 3 : 2;
						R scale = abs_value(x);
						{
							const R ay = abs_value(y);
							if (ay > scale) scale = ay;
							if (len == 3) {
								const R az = abs_value(z);
								if (az > scale) scale = az;
							}
						}
						if (!(scale > R(0))) {
							// bulge vanished exactly (double) or is not certifiable
							// (interval; the sweep then cannot make certified progress
							// and the cap provides honest termination)
							x = H[k + 1][k];
							y = (k + 2 <= i) ? H[k + 2][k] : T(0);
							z = (k + 3 <= i) ? H[k + 3][k] : T(0);
							continue;
						}
						const T vx = x / T(scale);
						const T vy = y / T(scale);
						const T vz = (len == 3) ? z / T(scale) : T(0);
						const R avx = abs_value(vx);
						const R avy = abs_value(vy);
						const R avz = abs_value(vz);
						const R nrm2 = (len == 3) ? (avx * avx + avy * avy + avz * avz)
						                          : (avx * avx + avy * avy);
						if (!(nrm2 > R(0))) {
							x = H[k + 1][k];
							y = (k + 2 <= i) ? H[k + 2][k] : T(0);
							z = (k + 3 <= i) ? H[k + 3][k] : T(0);
							continue;
						}
						const T nrm = T(sqrt_value(nrm2));
						const T sgn = (vx >= T(0)) ? T(1) : T(-1);
						const T u0 = vx + sgn * nrm;
						const T u1 = vy;
						const T u2 = (len == 3) ? vz : T(0);
						const R au0 = abs_value(u0);
						const R au1 = abs_value(u1);
						const R au2 = abs_value(u2);
						const R uu_r = (len == 3) ? (au0 * au0 + au1 * au1 + au2 * au2)
						                          : (au0 * au0 + au1 * au1);
						if (!(uu_r > R(0))) {
							x = H[k + 1][k];
							y = (k + 2 <= i) ? H[k + 2][k] : T(0);
							z = (k + 3 <= i) ? H[k + 3][k] : T(0);
							continue;
						}
						const T uu = T(uu_r);
						const std::size_t col_start = (k > l) ? k - 1 : l;
						const std::size_t row_end = (k + len < i) ? k + len : i;
						real_schur_detail::apply_householder_full(H, Zp, m, k, len,
						                                          col_start, row_end,
						                                          u0, u1, u2, uu);
						if (k > l) {
							// annihilated bulge column entries are analytically zero
							H[k + 1][k - 1] = T(0);
							if (len == 3) H[k + 2][k - 1] = T(0);
						}
						x = H[k + 1][k];
						y = (k + 2 <= i) ? H[k + 2][k] : T(0);
						z = (k + 3 <= i) ? H[k + 3][k] : T(0);
					}
					break;   // shifts were applied (or honestly skipped); leave attempt loop
				}

				res.iterations++;
				its++;
			}

			res.success = true;
			res.unconverged_count = 0;
			res.failure_reason.clear();
			return res;
		}

	} // namespace tsparse_real_schur
} // namespace vcp

#endif
