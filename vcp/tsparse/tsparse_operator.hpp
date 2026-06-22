// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_OPERATOR_HPP
#define VCP_TSPARSE_OPERATOR_HPP

#include <functional>
#include <vector>

#include <vcp/error.hpp>

namespace vcp {

	template <typename _T, typename _Index = int> class linear_operator {
	public:
		typedef _T value_type;
		typedef _Index index_type;

		virtual ~linear_operator() {}
		virtual index_type rowsize() const = 0;
		virtual index_type columnsize() const = 0;
		virtual void apply(const std::vector<value_type>& x, std::vector<value_type>& y) const = 0;
	};

	template <class SparseMatrix> class sparse_matrix_operator {
	public:
		typedef typename SparseMatrix::value_type value_type;
		typedef typename SparseMatrix::index_type index_type;

		explicit sparse_matrix_operator(const SparseMatrix& A) : A_(&A) {}

		index_type rowsize() const { return A_->rowsize(); }
		index_type columnsize() const { return A_->columnsize(); }

		void apply(const std::vector<value_type>& x, std::vector<value_type>& y) const {
			y = A_->mul_vec(x);
		}

	private:
		const SparseMatrix* A_;
	};

	template <typename _T, typename _Index = int> class function_linear_operator : public linear_operator<_T, _Index> {
	public:
		typedef _T value_type;
		typedef _Index index_type;
		typedef std::function<void(const std::vector<value_type>&, std::vector<value_type>&)> apply_type;

		function_linear_operator(index_type rows, index_type cols, const apply_type& apply)
			: rows_(rows), cols_(cols), apply_(apply) {}

		index_type rowsize() const { return rows_; }
		index_type columnsize() const { return cols_; }

		void apply(const std::vector<value_type>& x, std::vector<value_type>& y) const {
			apply_(x, y);
		}

	private:
		index_type rows_;
		index_type cols_;
		apply_type apply_;
	};
}

#endif
