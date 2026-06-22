// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_SPGEMM_HPP
#define VCP_TSPARSE_SPGEMM_HPP

#include <algorithm>
#include <vector>

#include <vcp/error.hpp>

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
	}
}

#endif
