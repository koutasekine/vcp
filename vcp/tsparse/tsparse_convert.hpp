// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_CONVERT_HPP
#define VCP_TSPARSE_CONVERT_HPP

#include <algorithm>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_format.hpp>

namespace vcp {
	namespace tsparse_detail {
		template <typename _Index, typename _T> struct coo_entry {
			_Index row;
			_Index col;
			_T value;
		};
	}

	template <typename _Index, typename _T>
	void tcoo_sort(_Index nnz, _Index* row, _Index* col, _T* value) {
		if (nnz < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcoo_sort: negative nnz");
		}
		std::vector<tsparse_detail::coo_entry<_Index, _T> > entries(static_cast<std::size_t>(nnz));
		for (_Index k = 0; k < nnz; k++) {
			entries[static_cast<std::size_t>(k)].row = row[k];
			entries[static_cast<std::size_t>(k)].col = col[k];
			entries[static_cast<std::size_t>(k)].value = value[k];
		}
		std::sort(entries.begin(), entries.end(),
			[](const tsparse_detail::coo_entry<_Index, _T>& a, const tsparse_detail::coo_entry<_Index, _T>& b) {
				return (a.row < b.row) || (a.row == b.row && a.col < b.col);
			});
		for (_Index k = 0; k < nnz; k++) {
			row[k] = entries[static_cast<std::size_t>(k)].row;
			col[k] = entries[static_cast<std::size_t>(k)].col;
			value[k] = entries[static_cast<std::size_t>(k)].value;
		}
	}

	template <typename _Index, typename _T>
	_Index tcoo_sum_duplicates(_Index nnz, _Index* row, _Index* col, _T* value) {
		if (nnz < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcoo_sum_duplicates: negative nnz");
		}
		if (nnz == 0) return 0;
		_Index out = 0;
		for (_Index k = 0; k < nnz; k++) {
			if (out > 0 && row[k] == row[out - 1] && col[k] == col[out - 1]) {
				value[out - 1] += value[k];
			}
			else {
				row[out] = row[k];
				col[out] = col[k];
				value[out] = value[k];
				out++;
			}
		}
		return out;
	}

	template <typename _Index, typename _T>
	_Index tcoo_remove_zeros(_Index nnz, _Index* row, _Index* col, _T* value) {
		if (nnz < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcoo_remove_zeros: negative nnz");
		}
		_Index out = 0;
		for (_Index k = 0; k < nnz; k++) {
			if (!(value[k] == _T(0))) {
				row[out] = row[k];
				col[out] = col[k];
				value[out] = value[k];
				out++;
			}
		}
		return out;
	}

	template <typename _Index, typename _T>
	void tcoo_to_csr(_Index rows, _Index cols, _Index nnz,
	                 const _Index* coo_row, const _Index* coo_col, const _T* coo_val,
	                 _Index* csr_row_ptr, _Index* csr_col_ind, _T* csr_val) {
		if (rows < 0 || cols < 0 || nnz < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcoo_to_csr: negative size");
		}
		for (_Index i = 0; i <= rows; i++) csr_row_ptr[i] = 0;
		for (_Index k = 0; k < nnz; k++) {
			if (coo_row[k] < 0 || coo_row[k] >= rows || coo_col[k] < 0 || coo_col[k] >= cols) {
				vcp::throw_error<vcp::index_error>("tcoo_to_csr: COO index out of range");
			}
			csr_row_ptr[coo_row[k] + 1]++;
		}
		for (_Index i = 0; i < rows; i++) csr_row_ptr[i + 1] += csr_row_ptr[i];
		std::vector<_Index> next(csr_row_ptr, csr_row_ptr + rows);
		for (_Index k = 0; k < nnz; k++) {
			_Index p = next[coo_row[k]]++;
			csr_col_ind[p] = coo_col[k];
			csr_val[p] = coo_val[k];
		}
	}

	template <typename _Index, typename _T>
	void tcoo_to_csc(_Index rows, _Index cols, _Index nnz,
	                 const _Index* coo_row, const _Index* coo_col, const _T* coo_val,
	                 _Index* csc_col_ptr, _Index* csc_row_ind, _T* csc_val) {
		if (rows < 0 || cols < 0 || nnz < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcoo_to_csc: negative size");
		}
		for (_Index j = 0; j <= cols; j++) csc_col_ptr[j] = 0;
		for (_Index k = 0; k < nnz; k++) {
			if (coo_row[k] < 0 || coo_row[k] >= rows || coo_col[k] < 0 || coo_col[k] >= cols) {
				vcp::throw_error<vcp::index_error>("tcoo_to_csc: COO index out of range");
			}
			csc_col_ptr[coo_col[k] + 1]++;
		}
		for (_Index j = 0; j < cols; j++) csc_col_ptr[j + 1] += csc_col_ptr[j];
		std::vector<_Index> next(csc_col_ptr, csc_col_ptr + cols);
		for (_Index k = 0; k < nnz; k++) {
			_Index p = next[coo_col[k]]++;
			csc_row_ind[p] = coo_row[k];
			csc_val[p] = coo_val[k];
		}
	}

	template <typename _Index, typename _T>
	void tcsr_to_coo(_Index rows, _Index cols, _Index nnz,
	                 const _Index* csr_row_ptr, const _Index* csr_col_ind, const _T* csr_val,
	                 _Index* coo_row, _Index* coo_col, _T* coo_val) {
		(void)cols;
		if (rows < 0 || nnz < 0) {
			vcp::throw_error<vcp::invalid_argument>("tcsr_to_coo: negative size");
		}
		for (_Index i = 0; i < rows; i++) {
			for (_Index p = csr_row_ptr[i]; p < csr_row_ptr[i + 1]; p++) {
				coo_row[p] = i;
				coo_col[p] = csr_col_ind[p];
				coo_val[p] = csr_val[p];
			}
		}
	}

	template <typename _Index, typename _T>
	void tcsr_to_csc(_Index rows, _Index cols, _Index nnz,
	                 const _Index* csr_row_ptr, const _Index* csr_col_ind, const _T* csr_val,
	                 _Index* csc_col_ptr, _Index* csc_row_ind, _T* csc_val) {
		std::vector<_Index> coo_row(static_cast<std::size_t>(nnz));
		std::vector<_Index> coo_col(static_cast<std::size_t>(nnz));
		std::vector<_T> coo_val(static_cast<std::size_t>(nnz));
		tcsr_to_coo(rows, cols, nnz, csr_row_ptr, csr_col_ind, csr_val, coo_row.data(), coo_col.data(), coo_val.data());
		tcoo_to_csc(rows, cols, nnz, coo_row.data(), coo_col.data(), coo_val.data(), csc_col_ptr, csc_row_ind, csc_val);
	}
}

#endif
