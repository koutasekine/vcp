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

#ifndef TBLAS_TBLAS_BLOCKED_HPP
#define TBLAS_TBLAS_BLOCKED_HPP

// SLU-K1: blocked GEMM kernels (design sandbox/docs/design/SLU-K1_design.md).
//
//   vcp::tblas_blocked::gemm<T>(transa, transb, m, n, k,
//                               alpha, A, lda, B, ldb, beta, C, ldc)
//
// - column-major, same argument order and edge-case semantics as vcp::tgemm.
// - ONLY the ('N','N') double path is blocked in this header (design §2 v1):
//     * skinny shapes (k <= 16 or n <= 16, the supernode-panel shapes) go to
//       a packing-free row-block kernel (design: "パッキング省略、tall-loop 直行");
//     * all other ('N','N') shapes go to a reduced Goto/BLIS 3-level blocking
//       (k_c, m_c, n_r x m_r register tile) with A/B packing;
//     * any transposed case forwards VERBATIM to the reference vcp::tgemm<T>.
// - every non-specialized T forwards VERBATIM to vcp::tgemm<T> (design D-7),
//   so the type matrix {double, dd, interval<double>, interval<dd>} compiles
//   and runs through this entry unconditionally.  The kv::dd specialization
//   (blocked-Ozaki, design §3) is added by SLU-K1 Phase 3 in this header.
// - numerical contract (design D-3): deterministic (fixed blocking, disjoint
//   output tiles, order independent of thread count), componentwise error
//   band |C_blocked - C_ref| <= c*k*eps*(|A|*|B|) -- NOT bit-identical to the
//   reference (summation order differs).
// - portable C++ only: no intrinsics, no inline assembly, no architecture
//   macros (K2's domain).  __restrict__ ("restrict 相当", design §2) plus
//   simple innermost loops carry the compiler auto-vectorization.
//
// Blocking parameters (design §2: one compile-time set, chosen from measured
// candidates on the implementation machine; per-machine tuning is K2).
// Measured 2026-07-18, i7-11700 (WSL2, g++ 13.3.0), single thread, min of 5,
// sluk1_02_blocked_bench sweep mode (deterministic LCG data), GFLOP/s:
//
//   square @256, candidates (MR x NR, KC, MC),  -O2  /  -O3 -march=native,
//   measured with the kernels compiled as ISOLATED functions (see
//   VCP_TBLAS_BLOCKED_KERNEL_FN below -- when the kernels were still
//   inlined into the entry, -O2 codegen oscillated 7.4 <-> 10.4 GFLOP/s
//   with unrelated header edits; isolation makes the numbers build-stable):
//     4 x 4  256 128 :  10.5 / 28.7
//     4 x 8  256 128 :  11.4 / 31.2   <- chosen
//     8 x 4  256 128 :   7.5 / 45.9   (native best, but -O2 spills the
//                                      32-double tile deterministically
//                                      to 1.40x ref < the 1.5x band)
//     8 x 8  256 128 :  ~10  /  5.7   (native COLLAPSE: g++ 13 mis-
//                                      vectorizes the 64-acc tile)
//     reference tgemm :  5.3 / 14.8
//
//   chosen MR=4, NR=8, KC=256, MC=128: the only leading candidate meeting
//   BOTH D-6 (1) bands deterministically -- -O2 2.14x ref (>= 1.5x MET),
//   native 2.11x ref (>= 1.0x MET).  8x4's 45.9 native is the K2
//   (arch-specific) upside, forfeited here for the portable-single-set
//   contract.  KC/MC insensitive in the measured range (@128/@512 equal).
//
//   skinny tile (MRS x NRS), m=4096 n=16, blocked/ref per k,  -O2 | native:
//     8 x 4 (chosen): k=2  3.3/5.8 | 11.2/9.4   k=4  4.4/6.2 | 21.4/8.8
//                     k=8  5.1/4.5 | 29.6/9.5   k=16 5.5/4.4 | 27.7/9.6
//     8 x 2         : k=2  4.5/5.8 |  7.5/9.4   k=8  7.3/4.5 | 12.3/9.5
//     16 x 4        : k=2  5.3/5.8 |  4.2/9.4   k=8  6.3/4.5 |  4.4/9.5
//   chosen MRS=8, NRS=4: native-dominant (2.4-3.1x ref at k >= 4; the
//   supernode-panel k distribution). At plain -O2 (SSE2 baseline) the
//   32-double accumulator tile spills (16 xmm registers), so small k
//   (<= 4) runs 0.6-0.7x ref there -- recorded as the D-6 (2) -O2 yellow;
//   no measured candidate reaches 1.3x at k <= 4 under -O2 because the
//   reference runs those shapes entirely from L2 (A block <= 128 KB).

#include <algorithm>
#include <cstddef>
#include <vector>

#include <vcp/tblas/tblas.hpp>

#if defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
#define VCP_TBLAS_BLOCKED_RESTRICT __restrict__
#else
#define VCP_TBLAS_BLOCKED_RESTRICT
#endif

// Compiler-portable inline barrier (same discipline as __restrict__ above:
// a compiler feature, not an architecture switch).  The packed-block and
// skinny kernels are compiled as isolated functions so that the size of the
// dispatching gemm<double> entry cannot perturb their code generation
// (measured at -O2: fully inlined into the entry, the 8x4 packed kernel
// dropped from 10.3-10.5 to 7.0-7.4 GFLOP/s when the entry grew by the
// density-dispatch code; the barrier restores the isolated-kernel numbers).
#if defined(__GNUC__) || defined(__clang__)
#define VCP_TBLAS_BLOCKED_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define VCP_TBLAS_BLOCKED_NOINLINE __declspec(noinline)
#else
#define VCP_TBLAS_BLOCKED_NOINLINE
#endif

namespace vcp {
namespace tblas_blocked {
namespace detail {

// -------------------------------------------------------------------------
// micro kernel: acc(MR x NR, column-major) := sum_{l<kc} ap(:,l) * bp(l,:)
// ap: MR x kc packed column panels (column l contiguous MR doubles),
// bp: kc x NR packed row panels (row l contiguous NR doubles).
// Fully accumulated in a local array so the compiler can keep the tile in
// registers; the caller applies alpha/beta at write-back.
// -------------------------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)
#define VCP_TBLAS_BLOCKED_KERNEL_FN __attribute__((noinline, aligned(64)))
#elif defined(_MSC_VER)
#define VCP_TBLAS_BLOCKED_KERNEL_FN __declspec(noinline)
#else
#define VCP_TBLAS_BLOCKED_KERNEL_FN
#endif

template <int MR, int NR>
VCP_TBLAS_BLOCKED_KERNEL_FN void micro_kernel_nn(
	const int kc,
	const double* VCP_TBLAS_BLOCKED_RESTRICT ap,
	const double* VCP_TBLAS_BLOCKED_RESTRICT bp,
	double* VCP_TBLAS_BLOCKED_RESTRICT acc
) {
	double c[MR * NR];
	for (int x = 0; x < MR * NR; x++) {
		c[x] = 0.0;
	}
	for (int l = 0; l < kc; l++) {
		const double* VCP_TBLAS_BLOCKED_RESTRICT a = ap + static_cast<std::size_t>(l) * MR;
		const double* VCP_TBLAS_BLOCKED_RESTRICT b = bp + static_cast<std::size_t>(l) * NR;
		for (int j = 0; j < NR; j++) {
			const double blj = b[j];
			for (int i = 0; i < MR; i++) {
				c[i + j * MR] += a[i] * blj;
			}
		}
	}
	for (int x = 0; x < MR * NR; x++) {
		acc[x] = c[x];
	}
}

// write back one register tile into C (only the effective mr_eff x nr_eff
// part; zero-padded lanes are dropped here, so padding never contaminates C)
template <int MR, int NR>
inline void write_tile_nn(
	const int mr_eff, const int nr_eff,
	const double* VCP_TBLAS_BLOCKED_RESTRICT acc,
	const double alpha, const double beta, const bool first_k_block,
	double* C, const int ldc
) {
	for (int j = 0; j < nr_eff; j++) {
		double* c = C + static_cast<std::size_t>(ldc) * j;
		const double* a = acc + j * MR;
		if (!first_k_block) {
			for (int i = 0; i < mr_eff; i++) {
				c[i] += alpha * a[i];
			}
		}
		else if (beta == 0.0) {
			// BLAS convention: beta == 0 never reads C
			for (int i = 0; i < mr_eff; i++) {
				c[i] = alpha * a[i];
			}
		}
		else if (beta == 1.0) {
			for (int i = 0; i < mr_eff; i++) {
				c[i] = c[i] + alpha * a[i];
			}
		}
		else {
			for (int i = 0; i < mr_eff; i++) {
				c[i] = beta * c[i] + alpha * a[i];
			}
		}
	}
}

// pack A(ic:ic+mc, pc:pc+kc) into MR-row panels: panel ir holds
// ap[ir/MR][l*MR + i] = A(ic+ir+i, pc+l), zero-padded past mc.
inline void pack_a_nn(
	const int mc, const int kc, const int MR,
	const double* A, const int lda,
	double* VCP_TBLAS_BLOCKED_RESTRICT ap
) {
	const int mpanels = (mc + MR - 1) / MR;
	for (int p = 0; p < mpanels; p++) {
		const int i0 = p * MR;
		const int ilim = std::min(MR, mc - i0);
		double* dst = ap + static_cast<std::size_t>(p) * MR * kc;
		for (int l = 0; l < kc; l++) {
			const double* a = A + i0 + static_cast<std::size_t>(lda) * l;
			double* d = dst + static_cast<std::size_t>(l) * MR;
			for (int i = 0; i < ilim; i++) {
				d[i] = a[i];
			}
			for (int i = ilim; i < MR; i++) {
				d[i] = 0.0;
			}
		}
	}
}

// pack B(pc:pc+kc, 0:n) into NR-column panels: panel jr holds
// bp[jr/NR][l*NR + j] = B(pc+l, jr+j), zero-padded past n.
inline void pack_b_nn(
	const int kc, const int n, const int NR,
	const double* B, const int ldb,
	double* VCP_TBLAS_BLOCKED_RESTRICT bp
) {
	const int npanels = (n + NR - 1) / NR;
	for (int p = 0; p < npanels; p++) {
		const int j0 = p * NR;
		const int jlim = std::min(NR, n - j0);
		double* dst = bp + static_cast<std::size_t>(p) * NR * kc;
		for (int l = 0; l < kc; l++) {
			double* d = dst + static_cast<std::size_t>(l) * NR;
			for (int j = 0; j < jlim; j++) {
				d[j] = B[l + static_cast<std::size_t>(ldb) * (j0 + j)];
			}
			for (int j = jlim; j < NR; j++) {
				d[j] = 0.0;
			}
		}
	}
}

// -------------------------------------------------------------------------
// main blocked ('N','N') kernel: C := alpha*A*B + beta*C.
// 3-level blocking: KC (k), MC (m), MR x NR register tiles; B packed once
// per k-block over the full n, A packed per (k-block, m-block).
// KC/MC are runtime parameters so the parameter-selection benchmark can
// sweep candidates; production entry uses the chosen constants.
// Deterministic: output tiles are disjoint and each tile's accumulation
// order is fixed, so the result is independent of the thread count.
// -------------------------------------------------------------------------
template <int MR, int NR>
VCP_TBLAS_BLOCKED_NOINLINE void gemm_nn_blocked(
	const int m, const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc,
	const int KC, const int MC
) {
	static thread_local std::vector<double> ap_buf;
	static thread_local std::vector<double> bp_buf;

	const int npanels = (n + NR - 1) / NR;
	bp_buf.resize(static_cast<std::size_t>(npanels) * NR * KC);

	for (int pc = 0; pc < k; pc += KC) {
		const int kc = std::min(KC, k - pc);
		const bool first_k_block = (pc == 0);
		pack_b_nn(kc, n, NR, B + pc, ldb, bp_buf.data());

		for (int ic = 0; ic < m; ic += MC) {
			const int mc = std::min(MC, m - ic);
			const int mpanels = (mc + MR - 1) / MR;
			ap_buf.resize(static_cast<std::size_t>(mpanels) * MR * kc);
			pack_a_nn(mc, kc, MR,
			          A + ic + static_cast<std::size_t>(lda) * pc, lda,
			          ap_buf.data());

			const double* apd = ap_buf.data();
			const double* bpd = bp_buf.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
	if (tblas_detail::use_parallel(2.0 * mc * n * kc))
#endif
			for (int jp = 0; jp < npanels; jp++) {
				const int j0 = jp * NR;
				const int nr_eff = std::min(NR, n - j0);
				const double* bpanel =
					bpd + static_cast<std::size_t>(jp) * NR * kc;
				double acc[MR * NR];
				for (int ip = 0; ip < mpanels; ip++) {
					const int i0 = ip * MR;
					const int mr_eff = std::min(MR, mc - i0);
					micro_kernel_nn<MR, NR>(
						kc,
						apd + static_cast<std::size_t>(ip) * MR * kc,
						bpanel, acc);
					write_tile_nn<MR, NR>(
						mr_eff, nr_eff, acc, alpha, beta, first_k_block,
						C + (ic + i0) +
							static_cast<std::size_t>(ldc) * j0, ldc);
				}
			}
		}
	}
}

// -------------------------------------------------------------------------
// skinny ('N','N') kernel (k <= 16 or n <= 16 -- the supernode-panel
// shapes): packing-free register-tile kernel ("tall-loop 直行").
// Column groups of NRS are processed with an MRS x NRS register tile that
// runs the WHOLE k loop before touching C, so
//   * each A row segment is loaded once per column group (the reference
//     re-streams all of A per single column, and re-reads/rewrites the C
//     column per k step -- at small k that C traffic dominates);
//   * C is written exactly once (and never read when beta == 0).
// A column group whose NRS B entries are all zero at some l skips that l
// (the supernode-panel fragments are zero-padded per column, design
// SLU-SP1 §2.2; the reference kernel skips those columns the same way).
// -------------------------------------------------------------------------
template <int MRS, int NRS>
VCP_TBLAS_BLOCKED_NOINLINE void skinny_col_group(
	const int m, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,   // B(0:k, j0:j0+NRS): first column base
	const double beta, double* C, const int ldc   // C(:, j0) base
) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
	if (tblas_detail::use_parallel(2.0 * m * NRS * k))
#endif
	for (int ib = 0; ib < m; ib += MRS) {
		const int mb = std::min(MRS, m - ib);
		double acc[MRS * NRS];
		for (int x = 0; x < MRS * NRS; x++) {
			acc[x] = 0.0;
		}
		if (mb == MRS) {
			for (int l = 0; l < k; l++) {
				double bl[NRS];
				bool all_zero = true;
				for (int j = 0; j < NRS; j++) {
					bl[j] = B[l + static_cast<std::size_t>(ldb) * j];
					if (bl[j] != 0.0) {
						all_zero = false;
					}
				}
				if (all_zero) {
					continue;
				}
				const double* VCP_TBLAS_BLOCKED_RESTRICT a =
					A + ib + static_cast<std::size_t>(lda) * l;
				for (int j = 0; j < NRS; j++) {
					const double blj = bl[j];
					for (int i = 0; i < MRS; i++) {
						acc[i + j * MRS] += a[i] * blj;
					}
				}
			}
		}
		else {
			for (int l = 0; l < k; l++) {
				const double* a =
					A + ib + static_cast<std::size_t>(lda) * l;
				for (int j = 0; j < NRS; j++) {
					const double blj =
						B[l + static_cast<std::size_t>(ldb) * j];
					for (int i = 0; i < mb; i++) {
						acc[i + j * MRS] += a[i] * blj;
					}
				}
			}
		}
		for (int j = 0; j < NRS; j++) {
			double* c = C + ib + static_cast<std::size_t>(ldc) * j;
			const double* a = acc + j * MRS;
			if (beta == 0.0) {
				for (int i = 0; i < mb; i++) {
					c[i] = alpha * a[i];
				}
			}
			else if (beta == 1.0) {
				for (int i = 0; i < mb; i++) {
					c[i] = c[i] + alpha * a[i];
				}
			}
			else {
				for (int i = 0; i < mb; i++) {
					c[i] = beta * c[i] + alpha * a[i];
				}
			}
		}
	}
}

// runtime-width column tail (nr < NRS): plain loops, same traffic pattern
template <int MRS, int NRS>
VCP_TBLAS_BLOCKED_NOINLINE void skinny_col_tail(
	const int m, const int nr, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc
) {
	for (int ib = 0; ib < m; ib += MRS) {
		const int mb = std::min(MRS, m - ib);
		double acc[MRS * NRS];
		for (int x = 0; x < MRS * NRS; x++) {
			acc[x] = 0.0;
		}
		for (int l = 0; l < k; l++) {
			const double* a = A + ib + static_cast<std::size_t>(lda) * l;
			for (int j = 0; j < nr; j++) {
				const double blj = B[l + static_cast<std::size_t>(ldb) * j];
				if (blj == 0.0) {
					continue;
				}
				for (int i = 0; i < mb; i++) {
					acc[i + j * MRS] += a[i] * blj;
				}
			}
		}
		for (int j = 0; j < nr; j++) {
			double* c = C + ib + static_cast<std::size_t>(ldc) * j;
			const double* a = acc + j * MRS;
			if (beta == 0.0) {
				for (int i = 0; i < mb; i++) {
					c[i] = alpha * a[i];
				}
			}
			else if (beta == 1.0) {
				for (int i = 0; i < mb; i++) {
					c[i] = c[i] + alpha * a[i];
				}
			}
			else {
				for (int i = 0; i < mb; i++) {
					c[i] = beta * c[i] + alpha * a[i];
				}
			}
		}
	}
}

template <int MRS, int NRS>
inline void gemm_nn_skinny(
	const int m, const int n, const int k,
	const double alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double beta, double* C, const int ldc
) {
	int j0 = 0;
	for (; j0 + NRS <= n; j0 += NRS) {
		skinny_col_group<MRS, NRS>(
			m, k, alpha, A, lda,
			B + static_cast<std::size_t>(ldb) * j0, ldb,
			beta, C + static_cast<std::size_t>(ldc) * j0, ldc);
	}
	if (j0 < n) {
		skinny_col_tail<MRS, NRS>(
			m, n - j0, k, alpha, A, lda,
			B + static_cast<std::size_t>(ldb) * j0, ldb,
			beta, C + static_cast<std::size_t>(ldc) * j0, ldc);
	}
}

// chosen production parameters (see the header comment measurement table)
const int gemm_mr = 4;
const int gemm_nr = 8;
const int gemm_kc = 256;
const int gemm_mc = 128;
const int skinny_mrs = 8;
const int skinny_nrs = 4;
const int skinny_max_k = 16;
const int skinny_max_n = 16;

// Shape/density dispatch inside the skinny path (SP1-realistic shapes,
// measured 2026-07-18 on the implementation machine, n = 8,
// m in {32,128,512,2048}, ref-vs-tile ratio, -O2 | native):
//   dense B      : k=2 tile 0.46-0.71x (LOSES both tiers)
//                  k=4 tile 0.78-1.09x (parity)
//                  k>=8 tile 1.10-1.43x | 1.30-2.53x (wins)
//   50%-zero B   : reference wins almost everywhere (its per-entry zero
//                  skip halves its work; the tile pays full flops) --
//                  k=8: 0.45-0.67x | 0.65-0.79x, parity only from k>=16.
// Rules (structure/value-determined only -> deterministic, D-3):
//   k <= skinny_ref_max_k(4)          -> reference (vcp::tgemm)
//   4*nnz(B) < 3*k*n (density < 3/4)  -> reference (zero-skip advantage)
//   otherwise                          -> MRS x NRS register tile
// The B scan is O(k*n), negligible against 2*m*n*k kernel work.
const int skinny_ref_max_k = 4;

} // namespace detail

// -------------------------------------------------------------------------
// gemm<T>: blocked GEMM entry point (design §2).
// Generic template: forwards VERBATIM to the reference vcp::tgemm<T>
// (design D-7 -- every T compiles and runs through this entry; only the
// double ('N','N') path below, and the kv::dd path added in Phase 3, are
// actually blocked).
// -------------------------------------------------------------------------
template <typename T>
inline void gemm(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const T& alpha, const T* A, const int lda,
	const T* B, const int ldb,
	const T& beta, T* C, const int ldc
) {
	vcp::tgemm<T>(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

// double specialization: ('N','N') only is blocked; transposed cases
// forward VERBATIM to the reference (v1 scope -- the supernode-panel call
// site only issues ('N','N'); the forward keeps the full tgemm contract).
template <>
inline void gemm<double>(
	const char transa, const char transb,
	const int m, const int n, const int k,
	const double& alpha, const double* A, const int lda,
	const double* B, const int ldb,
	const double& beta, double* C, const int ldc
) {
	namespace det = tblas_detail;
	if (!det::option_is(transa, 'N') || !det::option_is(transb, 'N')) {
		vcp::tgemm<double>(transa, transb, m, n, k, alpha, A, lda,
		                   B, ldb, beta, C, ldc);
		return;
	}
	// same argument validation and edge-case semantics as the reference
	if (m < 0 || n < 0 || k < 0 || lda < std::max(1, m) ||
	    ldb < std::max(1, k) || ldc < std::max(1, m)) {
		det::tblas_error("tblas_blocked::gemm<double>: invalid argument");
	}
	if (m == 0 || n == 0 || (k == 0 && beta == 1.0) ||
	    (alpha == 0.0 && beta == 1.0)) {
		return;
	}
	if (alpha == 0.0 || k == 0) {
		det::scale_matrix(m, n, beta, C, ldc);
		return;
	}
	if (k <= detail::skinny_max_k || n <= detail::skinny_max_n) {
		bool use_reference = (k <= detail::skinny_ref_max_k);
		if (!use_reference) {
			// density scan (see the dispatch comment in detail above)
			std::size_t nnz = 0u;
			for (int j = 0; j < n; j++) {
				const double* b = B + static_cast<std::size_t>(ldb) * j;
				for (int l = 0; l < k; l++) {
					if (b[l] != 0.0) {
						++nnz;
					}
				}
			}
			use_reference =
				(4u * nnz < 3u * static_cast<std::size_t>(k) *
				                static_cast<std::size_t>(n));
		}
		if (use_reference) {
			vcp::tgemm<double>('N', 'N', m, n, k, alpha, A, lda,
			                   B, ldb, beta, C, ldc);
		}
		else {
			detail::gemm_nn_skinny<detail::skinny_mrs, detail::skinny_nrs>(
				m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
		}
		return;
	}
	detail::gemm_nn_blocked<detail::gemm_mr, detail::gemm_nr>(
		m, n, k, alpha, A, lda, B, ldb, beta, C, ldc,
		detail::gemm_kc, detail::gemm_mc);
}

// -------------------------------------------------------------------------
// has_blocked_kernel<T>: true when gemm<T> above is an actually-blocked
// kernel (not the verbatim reference forward).  Read by the supernode_panel
// type-dependent GEMM gate (SLU-K1 Phase 3) so the gate can distinguish
// "fast blocked path exists" from "reference forward" WITHOUT tsparse
// naming any kv type.  The kv::dd specialization lives in the opt-in
// guarded section below, so its visibility follows the same include-order
// discipline as the kernel itself.
// -------------------------------------------------------------------------
template <typename T>
struct has_blocked_kernel {
	static const bool value = false;
};
template <>
struct has_blocked_kernel<double> {
	static const bool value = true;
};

} // namespace tblas_blocked
} // namespace vcp

// ===========================================================================
// SLU-K1 Phase 3 (v1.3 ruling, plan 1): blocked-Ozaki kv::dd gemm.
//
// OPT-IN section, same discipline as tlapack_dd.hpp: it is compiled ONLY
// when <vcp/tblas/tblas_dd.hpp> has been included BEFORE this header (e.g.
// before <vcp/spmatrix.hpp> in a factorization TU).  tsparse headers stay
// kv-agnostic; a TU without tblas_dd.hpp gets the verbatim reference
// forward for kv::dd through the generic gemm<T> above.
//
// Measured basis (T-7 before, sluk1_07, i7-11700, full-dd random @256):
// the double-gemm bundle is 86-93% of the current Ozaki tgemm<kv::dd> time
// (split 5-12%, merge 2-4%) -- the old D-4 split/merge hypothesis was
// REFUTED (STOP-6, owner ruling 2026-07-18: proceed with plan 1).  This
// section therefore keeps the Ozaki mathematics COMPLETELY UNCHANGED --
// pack / split / slice count / truncation / pair order / per-pair merge are
// the verbatim tblas_dd_detail helpers and loops -- and swaps ONLY the
// internal double gemm to the blocked kernel above.
// ===========================================================================
#if defined(TBLAS_TBLAS_DD_HPP)

namespace vcp {
namespace tblas_blocked {
namespace detail {

// ---------------------------------------------------------------------
// Array-oriented split / merge (design §3 "上積み"; adopted because the
// T-7 after-swap share of split+merge measured 25-34% at native > the
// 15% skip line).  The MATHEMATICS is unchanged: same split_shift /
// max_slices / per-slice extraction expression / break condition / merge
// accumulation chain as tblas_dd_detail.  Value-identity argument: the
// branch-free Knuth twosum below is an EXACT transform (x = fl(a+b),
// x + y == a + b), so for all finite data it produces the same (x, y)
// pair as kv::dd's magnitude-branch twosum, and every rem / accumulator
// update below replicates the corresponding kv::dd operator chain
// value-for-value (the kv infinity short-cuts are not replicated; near-
// overflow data is outside the Ozaki contract in either implementation).
// The dd arrays are kept as separate hi/lo double arrays so the hot
// loops are unit-stride, branch-poor and auto-vectorizable.
// ---------------------------------------------------------------------
inline void ozaki_twosum(const double a, const double b, double& x, double& y) {
	x = a + b;
	const double bb = x - a;
	y = (a - (x - bb)) + (b - bb);
}

// row-wise slice extraction of an m x k dd array (column-major),
// value-identical to tblas_dd_detail::split_row
inline void ozaki_split_rows(
	const std::vector< kv::dd >& A, const int m, const int k,
	std::vector< std::vector< double > >& slices
) {
	using std::fabs;
	using std::frexp;
	using std::ldexp;
	namespace ddd = ::vcp::tblas_dd_detail;
	const int s = ddd::split_shift(k);
	const int smax = ddd::max_slices(k);
	const std::size_t mk = static_cast<std::size_t>(m) * k;
	std::vector<double> rem_hi(mk), rem_lo(mk), mu(static_cast<std::size_t>(m)),
	    sigma(static_cast<std::size_t>(m));
	for (std::size_t x = 0; x < mk; x++) {
		rem_hi[x] = A[x].a1;
		rem_lo[x] = A[x].a2;
	}
	slices.clear();
	for (int p = 0; p < smax; p++) {
		for (int i = 0; i < m; i++) {
			mu[static_cast<std::size_t>(i)] = 0.0;
		}
		for (int j = 0; j < k; j++) {
			const double* VCP_TBLAS_BLOCKED_RESTRICT rh =
				&rem_hi[static_cast<std::size_t>(m) * j];
			double* VCP_TBLAS_BLOCKED_RESTRICT mv = &mu[0];
			for (int i = 0; i < m; i++) {
				const double d = fabs(rh[i]);
				if (d > mv[i]) {
					mv[i] = d;
				}
			}
		}
		for (int i = 0; i < m; i++) {
			int e;
			frexp(mu[static_cast<std::size_t>(i)], &e);
			sigma[static_cast<std::size_t>(i)] = ldexp(1.0, e + s);
		}
		std::vector<double> S(mk, 0.0);
		double nz = 0.0;
		for (int j = 0; j < k; j++) {
			const std::size_t base = static_cast<std::size_t>(m) * j;
			double* VCP_TBLAS_BLOCKED_RESTRICT rh = &rem_hi[base];
			double* VCP_TBLAS_BLOCKED_RESTRICT rl = &rem_lo[base];
			double* VCP_TBLAS_BLOCKED_RESTRICT sl = &S[base];
			const double* VCP_TBLAS_BLOCKED_RESTRICT sg = &sigma[0];
			for (int i = 0; i < m; i++) {
				const double q = (rh[i] + sg[i]) - sg[i];
				sl[i] = q;
				const double t = rh[i] - q;   // exact (q = leading bits)
				double z3, z4;
				ozaki_twosum(t, rl[i], z3, z4);
				const bool upd = (q != 0.0);
				rh[i] = upd ? z3 : rh[i];
				rl[i] = upd ? z4 : rl[i];
				nz += upd ? 1.0 : 0.0;
			}
		}
		if (nz == 0.0) {
			break;
		}
		slices.push_back(std::vector<double>());
		slices.back().swap(S);
	}
}

// column-wise slice extraction of a k x n dd array (column-major),
// value-identical to tblas_dd_detail::split_col
inline void ozaki_split_cols(
	const std::vector< kv::dd >& A, const int k, const int n,
	std::vector< std::vector< double > >& slices
) {
	using std::fabs;
	using std::frexp;
	using std::ldexp;
	namespace ddd = ::vcp::tblas_dd_detail;
	const int s = ddd::split_shift(k);
	const int smax = ddd::max_slices(k);
	const std::size_t kn = static_cast<std::size_t>(k) * n;
	std::vector<double> rem_hi(kn), rem_lo(kn);
	for (std::size_t x = 0; x < kn; x++) {
		rem_hi[x] = A[x].a1;
		rem_lo[x] = A[x].a2;
	}
	slices.clear();
	for (int p = 0; p < smax; p++) {
		std::vector<double> S(kn, 0.0);
		double nz = 0.0;
		for (int j = 0; j < n; j++) {
			const std::size_t base = static_cast<std::size_t>(k) * j;
			double* VCP_TBLAS_BLOCKED_RESTRICT rh = &rem_hi[base];
			double* VCP_TBLAS_BLOCKED_RESTRICT rl = &rem_lo[base];
			double* VCP_TBLAS_BLOCKED_RESTRICT sl = &S[base];
			double muj = 0.0;
			for (int i = 0; i < k; i++) {
				const double d = fabs(rh[i]);
				if (d > muj) {
					muj = d;
				}
			}
			if (muj == 0.0) {
				continue;   // same skip as split_col (column untouched)
			}
			int e;
			frexp(muj, &e);
			const double sg = ldexp(1.0, e + s);
			for (int i = 0; i < k; i++) {
				const double q = (rh[i] + sg) - sg;
				sl[i] = q;
				const double t = rh[i] - q;   // exact
				double z3, z4;
				ozaki_twosum(t, rl[i], z3, z4);
				const bool upd = (q != 0.0);
				rh[i] = upd ? z3 : rh[i];
				rl[i] = upd ? z4 : rl[i];
				nz += upd ? 1.0 : 0.0;
			}
		}
		if (nz == 0.0) {
			break;
		}
		slices.push_back(std::vector<double>());
		slices.back().swap(S);
	}
}

// C(hi,lo) += P elementwise with the exact kv::dd(+= double) chain and
// the same exact-zero skip as tblas_dd_detail::ozaki_product's merge
inline void ozaki_merge_accumulate(
	std::vector<double>& C_hi, std::vector<double>& C_lo,
	const std::vector<double>& P, const int mn
) {
	double* VCP_TBLAS_BLOCKED_RESTRICT ch = &C_hi[0];
	double* VCP_TBLAS_BLOCKED_RESTRICT cl = &C_lo[0];
	const double* VCP_TBLAS_BLOCKED_RESTRICT pp = &P[0];
#ifdef _OPENMP
#pragma omp parallel for schedule(static) \
	if (tblas_detail::use_parallel(static_cast<double>(mn)))
#endif
	for (int idx = 0; idx < mn; idx++) {
		const double p = pp[idx];
		double z1, z2, z3, z4;
		ozaki_twosum(ch[idx], p, z1, z2);
		z2 += cl[idx];
		ozaki_twosum(z1, z2, z3, z4);
		const bool upd = (p != 0.0);
		ch[idx] = upd ? z3 : ch[idx];
		cl[idx] = upd ? z4 : cl[idx];
	}
}

// tblas_dd_detail::ozaki_product with (a) the internal double gemm bundle
// swapped to the blocked kernel (the D-4 v1.3 main work: the bundle is
// 83-93% of the Ozaki time) and (b) split/merge in the array-oriented
// form above (the "上積み"); slice count / truncation / pair order are
// the verbatim original policy.
inline void ozaki_product_blocked(
	const std::vector< kv::dd >& opA, const std::vector< kv::dd >& opB,
	const int m, const int n, const int k,
	std::vector< kv::dd >& C
) {
	namespace ddd = ::vcp::tblas_dd_detail;
	std::vector< std::vector< double > > DA;
	std::vector< std::vector< double > > DB;
	ozaki_split_rows(opA, m, k, DA);
	ozaki_split_cols(opB, k, n, DB);

	C.assign(static_cast<std::size_t>(m) * n, kv::dd(0.0));
	if (DA.empty() || DB.empty()) {
		return;
	}

	const double one = 1.0;
	const double zero = 0.0;
	const int limit = ddd::max_slices(k);
	const std::size_t smn = static_cast<std::size_t>(m) * n;
	std::vector< double > P(smn, 0.0);
	std::vector< double > C_hi(smn, 0.0);
	std::vector< double > C_lo(smn, 0.0);
	for (int t = 0; t < limit; t++) {
		for (int p = 0; p <= t; p++) {
			const int q = t - p;
			if (p >= static_cast<int>(DA.size()) || q >= static_cast<int>(DB.size())) {
				continue;
			}
			::vcp::tblas_blocked::gemm<double>(
				'N', 'N', m, n, k, one, DA[p].data(), m,
				DB[q].data(), k, zero, P.data(), m);
			ozaki_merge_accumulate(C_hi, C_lo, P, m * n);
		}
	}
	for (std::size_t x = 0; x < smn; x++) {
		C[x] = kv::dd(C_hi[x], C_lo[x]);
	}
}

} // namespace detail

// blocked-Ozaki kv::dd gemm: verbatim vcp::tgemm<kv::dd> entry semantics
// (validation, edge cases, transposes via pack_op, alpha/beta application)
// with detail::ozaki_product_blocked as the product core.
template <>
inline void gemm<kv::dd>(
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
		det::tblas_error("tblas_blocked::gemm<kv::dd>: invalid transa");
	}
	if (tb && !det::option_is(transb, 'T') && !det::option_is(transb, 'C')) {
		det::tblas_error("tblas_blocked::gemm<kv::dd>: invalid transb");
	}
	const int nrowa = ta ? k : m;
	const int nrowb = tb ? n : k;
	if (m < 0 || n < 0 || k < 0 || lda < std::max(1, nrowa) || ldb < std::max(1, nrowb) || ldc < std::max(1, m)) {
		det::tblas_error("tblas_blocked::gemm<kv::dd>: invalid argument");
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
	detail::ozaki_product_blocked(opA, opB, m, n, k, prod);

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
struct has_blocked_kernel<kv::dd> {
	static const bool value = true;
};

} // namespace tblas_blocked
} // namespace vcp

#endif // TBLAS_TBLAS_DD_HPP (blocked-Ozaki dd section)

#endif // TBLAS_TBLAS_BLOCKED_HPP
