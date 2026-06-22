// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_FORMAT_HPP
#define VCP_TSPARSE_FORMAT_HPP

namespace vcp {

	enum sparse_format {
		sparse_coo,
		sparse_csr,
		sparse_csc
	};

	template <typename _Index, typename _T> struct coo_view {
		typedef _Index index_type;
		typedef _T value_type;

		_Index rows;
		_Index cols;
		_Index nnz;
		const _Index* row;
		const _Index* col;
		const _T* value;
	};

	template <typename _Index, typename _T> struct csr_view {
		typedef _Index index_type;
		typedef _T value_type;

		_Index rows;
		_Index cols;
		_Index nnz;
		const _Index* row_ptr;
		const _Index* col_ind;
		const _T* value;
	};

	template <typename _Index, typename _T> struct csc_view {
		typedef _Index index_type;
		typedef _T value_type;

		_Index rows;
		_Index cols;
		_Index nnz;
		const _Index* col_ptr;
		const _Index* row_ind;
		const _T* value;
	};
}

#endif
