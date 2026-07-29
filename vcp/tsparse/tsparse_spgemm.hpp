// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_SPGEMM_HPP
#define VCP_TSPARSE_SPGEMM_HPP

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <vector>

#include <vcp/error.hpp>

// SPOMP-1 (design §2, D-1): two-stage OpenMP guard (replicated per file so
// each header stands alone; same block as spmats.hpp / tsparse_blas.hpp).
#ifdef VCP_NOMP
#  ifndef VCP_SPARSE_NOMP
#    define VCP_SPARSE_NOMP
#  endif
#endif
#if defined(_OPENMP) && !defined(VCP_SPARSE_NOMP)
#  define VCP_SPARSE_USE_OPENMP 1
#  include <omp.h>
#else
#  define VCP_SPARSE_USE_OPENMP 0
#endif

namespace vcp {
	namespace tsparse_spgemm {
		template <typename Index, typename T, class Emit>
		void csr_csr_linear_combination(Index rows, Index cols,
		                                const std::vector<Index>& ao, const std::vector<Index>& ai, const std::vector<T>& av,
		                                const std::vector<Index>& bo, const std::vector<Index>& bi, const std::vector<T>& bv,
		                                const T& alpha, const T& beta, Emit emit) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("csr_csr_linear_combination: negative size");
			}
			std::vector<char> used(static_cast<std::size_t>(cols), 0);
			std::vector<T> workspace(static_cast<std::size_t>(cols), T(0));
			std::vector<Index> touched;
			for (Index i = 0; i < rows; i++) {
				touched.clear();
				for (Index ap = ao[static_cast<std::size_t>(i)]; ap < ao[static_cast<std::size_t>(i + 1)]; ap++) {
					const Index j = ai[static_cast<std::size_t>(ap)];
					const std::size_t js = static_cast<std::size_t>(j);
					if (!used[js]) {
						used[js] = 1;
						touched.push_back(j);
					}
					workspace[js] += alpha * av[static_cast<std::size_t>(ap)];
				}
				for (Index bp = bo[static_cast<std::size_t>(i)]; bp < bo[static_cast<std::size_t>(i + 1)]; bp++) {
					const Index j = bi[static_cast<std::size_t>(bp)];
					const std::size_t js = static_cast<std::size_t>(j);
					if (!used[js]) {
						used[js] = 1;
						touched.push_back(j);
					}
					workspace[js] += beta * bv[static_cast<std::size_t>(bp)];
				}
				std::sort(touched.begin(), touched.end());
				for (std::size_t p = 0; p < touched.size(); p++) {
					const Index j = touched[p];
					const std::size_t js = static_cast<std::size_t>(j);
					if (!(workspace[js] == T(0))) emit(i, j, workspace[js]);
					workspace[js] = T(0);
					used[js] = 0;
				}
			}
		}

		template <typename Index, typename T, class Emit>
		void csr_csr_multiply(Index rows, Index inner_dim, Index cols,
		                      const std::vector<Index>& ao, const std::vector<Index>& ai, const std::vector<T>& av,
		                      const std::vector<Index>& bo, const std::vector<Index>& bi, const std::vector<T>& bv,
		                      Emit emit) {
			(void)inner_dim;
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("csr_csr_multiply: negative size");
			}
			std::vector<char> used(static_cast<std::size_t>(cols), 0);
			std::vector<T> workspace(static_cast<std::size_t>(cols), T(0));
			std::vector<Index> touched;
			for (Index i = 0; i < rows; i++) {
				touched.clear();
				for (Index ap = ao[static_cast<std::size_t>(i)]; ap < ao[static_cast<std::size_t>(i + 1)]; ap++) {
					const Index k = ai[static_cast<std::size_t>(ap)];
					const T aik = av[static_cast<std::size_t>(ap)];
					for (Index bp = bo[static_cast<std::size_t>(k)]; bp < bo[static_cast<std::size_t>(k + 1)]; bp++) {
						const Index j = bi[static_cast<std::size_t>(bp)];
						const std::size_t js = static_cast<std::size_t>(j);
						if (!used[js]) {
							used[js] = 1;
							touched.push_back(j);
						}
						workspace[js] += aik * bv[static_cast<std::size_t>(bp)];
					}
				}
				std::sort(touched.begin(), touched.end());
				for (std::size_t p = 0; p < touched.size(); p++) {
					const Index j = touched[p];
					const std::size_t js = static_cast<std::size_t>(j);
					if (!(workspace[js] == T(0))) emit(i, j, workspace[js]);
					workspace[js] = T(0);
					used[js] = 0;
				}
			}
		}

		// -------------------------------------------------------------------
		// SPOMP-1 (design §4.1, D-9): two-pass parallel CSR builders.
		// PURE ADDITIONS -- the emit-style kernels above are untouched and
		// keep serving every existing caller.  Output nnz is value-dependent
		// (exact zeros are dropped), so a one-pass parallel emit cannot fix
		// output positions; instead:
		//   pass 1 (parallel): rows are split into CONTIGUOUS blocks (one
		//     parallel-for iteration per block = explicit static partition);
		//     each block runs the IDENTICAL per-row algorithm of the
		//     sequential kernel (thread-local workspace/used/touched
		//     allocated inside the parallel region), pushing (inner, value)
		//     into a block-local buffer and recording per-row output counts.
		//   sequential prefix sum over the row counts fixes the output outer.
		//   pass 2 (parallel): each block copies its local buffer to the
		//     final inner/value at its fixed offset (no arithmetic).
		// Row order inside a block and the per-row operation sequence equal
		// the sequential kernel, so the built CSR is bit-identical for any
		// thread count and satisfies the same output contract as the emit
		// kernels (ascending rows, ascending columns within a row,
		// duplicate-merged, exact-zero-free).  Exceptions are captured once
		// (atomic election) and rethrown outside the region (design §6-1);
		// loop variables are signed (design §6-2).  When
		// VCP_SPARSE_USE_OPENMP == 0 this runs as a single block,
		// sequentially, with the same operation sequence.
		// -------------------------------------------------------------------
		template <typename Index, typename T>
		void csr_csr_linear_combination_par(Index rows, Index cols,
		                                    const std::vector<Index>& ao, const std::vector<Index>& ai, const std::vector<T>& av,
		                                    const std::vector<Index>& bo, const std::vector<Index>& bi, const std::vector<T>& bv,
		                                    const T& alpha, const T& beta,
		                                    std::vector<Index>& co, std::vector<Index>& ci, std::vector<T>& cv) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("csr_csr_linear_combination_par: negative size");
			}
			const std::ptrdiff_t nrows = static_cast<std::ptrdiff_t>(rows);
			std::ptrdiff_t nblocks = 1;
#if VCP_SPARSE_USE_OPENMP
			nblocks = static_cast<std::ptrdiff_t>(omp_get_max_threads());
			if (nblocks < 1) nblocks = 1;
			if (nblocks > nrows) nblocks = (nrows > 0) ? nrows : 1;
#endif
			std::vector<std::vector<Index> > block_inner(static_cast<std::size_t>(nblocks));
			std::vector<std::vector<T> > block_value(static_cast<std::size_t>(nblocks));
			std::vector<Index> row_cnt(static_cast<std::size_t>(nrows), Index(0));
			std::atomic<bool> caught(false);
			std::exception_ptr eptr;
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (std::ptrdiff_t b = 0; b < nblocks; b++) {
				try {
					const std::ptrdiff_t r0 = b * nrows / nblocks;
					const std::ptrdiff_t r1 = (b + 1) * nrows / nblocks;
					std::vector<char> used(static_cast<std::size_t>(cols), 0);
					std::vector<T> workspace(static_cast<std::size_t>(cols), T(0));
					std::vector<Index> touched;
					std::vector<Index>& li = block_inner[static_cast<std::size_t>(b)];
					std::vector<T>& lv = block_value[static_cast<std::size_t>(b)];
					for (std::ptrdiff_t i = r0; i < r1; i++) {
						touched.clear();
						const std::size_t before = li.size();
						for (Index ap = ao[static_cast<std::size_t>(i)]; ap < ao[static_cast<std::size_t>(i + 1)]; ap++) {
							const Index j = ai[static_cast<std::size_t>(ap)];
							const std::size_t js = static_cast<std::size_t>(j);
							if (!used[js]) {
								used[js] = 1;
								touched.push_back(j);
							}
							workspace[js] += alpha * av[static_cast<std::size_t>(ap)];
						}
						for (Index bp = bo[static_cast<std::size_t>(i)]; bp < bo[static_cast<std::size_t>(i + 1)]; bp++) {
							const Index j = bi[static_cast<std::size_t>(bp)];
							const std::size_t js = static_cast<std::size_t>(j);
							if (!used[js]) {
								used[js] = 1;
								touched.push_back(j);
							}
							workspace[js] += beta * bv[static_cast<std::size_t>(bp)];
						}
						std::sort(touched.begin(), touched.end());
						for (std::size_t p = 0; p < touched.size(); p++) {
							const Index j = touched[p];
							const std::size_t js = static_cast<std::size_t>(j);
							if (!(workspace[js] == T(0))) {
								li.push_back(j);
								lv.push_back(workspace[js]);
							}
							workspace[js] = T(0);
							used[js] = 0;
						}
						row_cnt[static_cast<std::size_t>(i)] = static_cast<Index>(li.size() - before);
					}
				}
				catch (...) {
					bool expected = false;
					if (caught.compare_exchange_strong(expected, true)) {
						eptr = std::current_exception();
					}
				}
			}
			if (caught.load()) std::rethrow_exception(eptr);
			co.assign(static_cast<std::size_t>(nrows) + 1, Index(0));
			for (std::ptrdiff_t i = 0; i < nrows; i++) {
				co[static_cast<std::size_t>(i) + 1] = co[static_cast<std::size_t>(i)] + row_cnt[static_cast<std::size_t>(i)];
			}
			ci.resize(static_cast<std::size_t>(co[static_cast<std::size_t>(nrows)]));
			cv.resize(static_cast<std::size_t>(co[static_cast<std::size_t>(nrows)]), T(0));
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (std::ptrdiff_t b = 0; b < nblocks; b++) {
				try {
					const std::ptrdiff_t r0 = b * nrows / nblocks;
					const std::vector<Index>& li = block_inner[static_cast<std::size_t>(b)];
					const std::vector<T>& lv = block_value[static_cast<std::size_t>(b)];
					const std::size_t dst = static_cast<std::size_t>(co[static_cast<std::size_t>(r0)]);
					for (std::size_t k = 0; k < li.size(); k++) {
						ci[dst + k] = li[k];
						cv[dst + k] = lv[k];
					}
				}
				catch (...) {
					bool expected = false;
					if (caught.compare_exchange_strong(expected, true)) {
						eptr = std::current_exception();
					}
				}
			}
			if (caught.load()) std::rethrow_exception(eptr);
		}

		template <typename Index, typename T>
		void csr_csr_multiply_par(Index rows, Index inner_dim, Index cols,
		                          const std::vector<Index>& ao, const std::vector<Index>& ai, const std::vector<T>& av,
		                          const std::vector<Index>& bo, const std::vector<Index>& bi, const std::vector<T>& bv,
		                          std::vector<Index>& co, std::vector<Index>& ci, std::vector<T>& cv) {
			(void)inner_dim;
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("csr_csr_multiply_par: negative size");
			}
			const std::ptrdiff_t nrows = static_cast<std::ptrdiff_t>(rows);
			std::ptrdiff_t nblocks = 1;
#if VCP_SPARSE_USE_OPENMP
			nblocks = static_cast<std::ptrdiff_t>(omp_get_max_threads());
			if (nblocks < 1) nblocks = 1;
			if (nblocks > nrows) nblocks = (nrows > 0) ? nrows : 1;
#endif
			std::vector<std::vector<Index> > block_inner(static_cast<std::size_t>(nblocks));
			std::vector<std::vector<T> > block_value(static_cast<std::size_t>(nblocks));
			std::vector<Index> row_cnt(static_cast<std::size_t>(nrows), Index(0));
			std::atomic<bool> caught(false);
			std::exception_ptr eptr;
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (std::ptrdiff_t b = 0; b < nblocks; b++) {
				try {
					const std::ptrdiff_t r0 = b * nrows / nblocks;
					const std::ptrdiff_t r1 = (b + 1) * nrows / nblocks;
					std::vector<char> used(static_cast<std::size_t>(cols), 0);
					std::vector<T> workspace(static_cast<std::size_t>(cols), T(0));
					std::vector<Index> touched;
					std::vector<Index>& li = block_inner[static_cast<std::size_t>(b)];
					std::vector<T>& lv = block_value[static_cast<std::size_t>(b)];
					for (std::ptrdiff_t i = r0; i < r1; i++) {
						touched.clear();
						const std::size_t before = li.size();
						for (Index ap = ao[static_cast<std::size_t>(i)]; ap < ao[static_cast<std::size_t>(i + 1)]; ap++) {
							const Index k = ai[static_cast<std::size_t>(ap)];
							const T aik = av[static_cast<std::size_t>(ap)];
							for (Index bp = bo[static_cast<std::size_t>(k)]; bp < bo[static_cast<std::size_t>(k + 1)]; bp++) {
								const Index j = bi[static_cast<std::size_t>(bp)];
								const std::size_t js = static_cast<std::size_t>(j);
								if (!used[js]) {
									used[js] = 1;
									touched.push_back(j);
								}
								workspace[js] += aik * bv[static_cast<std::size_t>(bp)];
							}
						}
						std::sort(touched.begin(), touched.end());
						for (std::size_t p = 0; p < touched.size(); p++) {
							const Index j = touched[p];
							const std::size_t js = static_cast<std::size_t>(j);
							if (!(workspace[js] == T(0))) {
								li.push_back(j);
								lv.push_back(workspace[js]);
							}
							workspace[js] = T(0);
							used[js] = 0;
						}
						row_cnt[static_cast<std::size_t>(i)] = static_cast<Index>(li.size() - before);
					}
				}
				catch (...) {
					bool expected = false;
					if (caught.compare_exchange_strong(expected, true)) {
						eptr = std::current_exception();
					}
				}
			}
			if (caught.load()) std::rethrow_exception(eptr);
			co.assign(static_cast<std::size_t>(nrows) + 1, Index(0));
			for (std::ptrdiff_t i = 0; i < nrows; i++) {
				co[static_cast<std::size_t>(i) + 1] = co[static_cast<std::size_t>(i)] + row_cnt[static_cast<std::size_t>(i)];
			}
			ci.resize(static_cast<std::size_t>(co[static_cast<std::size_t>(nrows)]));
			cv.resize(static_cast<std::size_t>(co[static_cast<std::size_t>(nrows)]), T(0));
#if VCP_SPARSE_USE_OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (std::ptrdiff_t b = 0; b < nblocks; b++) {
				try {
					const std::ptrdiff_t r0 = b * nrows / nblocks;
					const std::vector<Index>& li = block_inner[static_cast<std::size_t>(b)];
					const std::vector<T>& lv = block_value[static_cast<std::size_t>(b)];
					const std::size_t dst = static_cast<std::size_t>(co[static_cast<std::size_t>(r0)]);
					for (std::size_t k = 0; k < li.size(); k++) {
						ci[dst + k] = li[k];
						cv[dst + k] = lv[k];
					}
				}
				catch (...) {
					bool expected = false;
					if (caught.compare_exchange_strong(expected, true)) {
						eptr = std::current_exception();
					}
				}
			}
			if (caught.load()) std::rethrow_exception(eptr);
		}
	}
}

#endif
