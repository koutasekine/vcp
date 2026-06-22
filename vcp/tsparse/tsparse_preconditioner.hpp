// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_PRECONDITIONER_HPP
#define VCP_TSPARSE_PRECONDITIONER_HPP

#include <vector>

#include <vcp/error.hpp>

namespace vcp {

	template <typename _T> struct identity_preconditioner {
		typedef _T value_type;

		void apply(const std::vector<value_type>& r, std::vector<value_type>& z) const {
			z = r;
		}
	};

	template <class SparseMatrix> class jacobi_preconditioner {
	public:
		typedef typename SparseMatrix::value_type value_type;
		typedef typename SparseMatrix::index_type index_type;

		explicit jacobi_preconditioner(const SparseMatrix& A) {
			const index_type n = A.rowsize();
			inv_diag_.assign(static_cast<std::size_t>(n), value_type(0));
			for (index_type i = 0; i < n; i++) {
				const value_type diag = A.get(i, i);
				if (diag == value_type(0)) {
					vcp::throw_error<vcp::numerical_error>("jacobi_preconditioner: zero diagonal");
				}
				inv_diag_[static_cast<std::size_t>(i)] = value_type(1) / diag;
			}
		}

		void apply(const std::vector<value_type>& r, std::vector<value_type>& z) const {
			if (r.size() != inv_diag_.size()) {
				vcp::throw_error<vcp::dimension_error>("jacobi_preconditioner::apply: dimension mismatch");
			}
			z.assign(r.size(), value_type(0));
			for (std::size_t i = 0; i < r.size(); i++) z[i] = inv_diag_[i] * r[i];
		}

	private:
		std::vector<value_type> inv_diag_;
	};
}

#endif
