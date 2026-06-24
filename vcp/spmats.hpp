// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATS_HPP
#define VCP_SPMATS_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse.hpp>
#include <vcp/spmats_eigs_types.hpp>
#include <vcp/spmats_policy_traits.hpp>

namespace vcp {

	template <typename _T, typename _Index = int> class spmats {
	public:
		typedef _Index index_type;
		typedef _T value_type;
		typedef vcp::sparse_format format_type;

		spmats()
			: row(0), column(0), fmt(vcp::sparse_coo), finalized(false), sorted(true), unique(true) {}

		spmats(const spmats&) = default;
		spmats(spmats&&) = default;
		spmats& operator=(const spmats&) = default;
		spmats& operator=(spmats&&) = default;
		virtual ~spmats() = default;

		index_type rowsize() const { return row; }
		index_type columnsize() const { return column; }
		index_type stored_nnz() const {
			return finalized ? size_to_index(value.size(), "spmats::stored_nnz") : size_to_index(coo_value.size(), "spmats::stored_nnz");
		}
		index_type nnz() const {
			if (finalized) return size_to_index(value.size(), "spmats::nnz");
			index_type count = 0;
			for (std::size_t k = 0; k < coo_value.size(); k++) {
				if (!(coo_value[k] == _T(0))) count++;
			}
			return count;
		}
		bool is_finalized() const { return finalized; }
		bool is_sorted() const { return sorted; }
		bool is_unique() const { return unique; }
		format_type format() const { return fmt; }

		void resize(const index_type rows, const index_type cols) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::resize: negative size");
			}
			row = rows;
			column = cols;
			clear_storage();
			fmt = vcp::sparse_coo;
			finalized = false;
			sorted = true;
			unique = true;
		}

		void clear() {
			row = 0;
			column = 0;
			clear_storage();
			fmt = vcp::sparse_coo;
			finalized = false;
			sorted = true;
			unique = true;
		}

		void reserve(const index_type n) {
			if (n < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::reserve: negative size");
			}
			const std::size_t m = index_to_size(n, "spmats::reserve");
			coo_row.reserve(m);
			coo_col.reserve(m);
			coo_value.reserve(m);
			inner.reserve(m);
			value.reserve(m);
		}

		void add(const index_type i, const index_type j, const _T& a) {
			check_index(i, j, "spmats::add");
			ensure_coo_buffer();
			coo_row.push_back(i);
			coo_col.push_back(j);
			coo_value.push_back(a);
			finalized = false;
			fmt = vcp::sparse_coo;
			sorted = false;
			unique = false;
		}

		void set(const index_type i, const index_type j, const _T& a) {
			check_index(i, j, "spmats::set");
			ensure_coo_buffer();
			std::vector<index_type> nr;
			std::vector<index_type> nc;
			std::vector<_T> nv;
			nr.reserve(coo_value.size() + 1);
			nc.reserve(coo_value.size() + 1);
			nv.reserve(coo_value.size() + 1);
			for (std::size_t k = 0; k < coo_value.size(); k++) {
				if (!(coo_row[k] == i && coo_col[k] == j)) {
					nr.push_back(coo_row[k]);
					nc.push_back(coo_col[k]);
					nv.push_back(coo_value[k]);
				}
			}
			if (!(a == _T(0))) {
				nr.push_back(i);
				nc.push_back(j);
				nv.push_back(a);
			}
			coo_row.swap(nr);
			coo_col.swap(nc);
			coo_value.swap(nv);
			finalized = false;
			fmt = vcp::sparse_coo;
			sorted = false;
			unique = true;
		}

		_T get(const index_type i, const index_type j) const {
			check_index(i, j, "spmats::get");
			if (!finalized) {
				_T sum = _T(0);
				for (std::size_t k = 0; k < coo_value.size(); k++) {
					if (coo_row[k] == i && coo_col[k] == j) sum += coo_value[k];
				}
				return sum;
			}
			if (fmt == vcp::sparse_csr) {
				const index_type first = outer[i];
				const index_type last = outer[i + 1];
				const typename std::vector<index_type>::const_iterator begin = inner.begin() + first;
				const typename std::vector<index_type>::const_iterator end = inner.begin() + last;
				typename std::vector<index_type>::const_iterator it = std::lower_bound(begin, end, j);
				if (it != end && *it == j) {
					return value[static_cast<std::size_t>(it - inner.begin())];
				}
				return _T(0);
			}
			if (fmt == vcp::sparse_csc) {
				const index_type first = outer[j];
				const index_type last = outer[j + 1];
				const typename std::vector<index_type>::const_iterator begin = inner.begin() + first;
				const typename std::vector<index_type>::const_iterator end = inner.begin() + last;
				typename std::vector<index_type>::const_iterator it = std::lower_bound(begin, end, i);
				if (it != end && *it == i) {
					return value[static_cast<std::size_t>(it - inner.begin())];
				}
				return _T(0);
			}
			return _T(0);
		}

		void finalize() { to_csr(); }

		void sort_coo() {
			ensure_coo_buffer();
			const index_type n = size_to_index(coo_value.size(), "spmats::sort_coo");
			if (n > 0) {
				vcp::tcoo_sort(n, coo_row.data(), coo_col.data(), coo_value.data());
			}
			sorted = true;
		}

		void normalize_coo() {
			ensure_coo_buffer();
			sort_coo();
			index_type n = size_to_index(coo_value.size(), "spmats::normalize_coo");
			if (n > 0) {
				n = vcp::tcoo_sum_duplicates(n, coo_row.data(), coo_col.data(), coo_value.data());
				n = vcp::tcoo_remove_zeros(n, coo_row.data(), coo_col.data(), coo_value.data());
			}
			coo_row.resize(static_cast<std::size_t>(n));
			coo_col.resize(static_cast<std::size_t>(n));
			coo_value.resize(static_cast<std::size_t>(n));
			sorted = true;
			unique = true;
		}

		void to_csr() {
			ensure_coo_buffer();
			normalize_coo();
			const index_type n = size_to_index(coo_value.size(), "spmats::to_csr");
			outer.assign(index_to_size(checked_plus_one(row, "spmats::to_csr"), "spmats::to_csr"), 0);
			inner.assign(index_to_size(n, "spmats::to_csr"), 0);
			value.assign(index_to_size(n, "spmats::to_csr"), _T(0));
			if (n > 0 || row >= 0) {
				vcp::tcoo_to_csr(row, column, n, coo_row.data(), coo_col.data(), coo_value.data(), outer.data(), inner.data(), value.data());
			}
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			fmt = vcp::sparse_csr;
			finalized = true;
			sorted = true;
			unique = true;
		}

		void to_csc() {
			ensure_coo_buffer();
			normalize_coo();
			const index_type n = size_to_index(coo_value.size(), "spmats::to_csc");
			outer.assign(index_to_size(checked_plus_one(column, "spmats::to_csc"), "spmats::to_csc"), 0);
			inner.assign(index_to_size(n, "spmats::to_csc"), 0);
			value.assign(index_to_size(n, "spmats::to_csc"), _T(0));
			if (n > 0 || column >= 0) {
				vcp::tcoo_to_csc(row, column, n, coo_row.data(), coo_col.data(), coo_value.data(), outer.data(), inner.data(), value.data());
			}
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			fmt = vcp::sparse_csc;
			finalized = true;
			sorted = true;
			unique = true;
		}

		spmats as_csr() const {
			spmats tmp(*this);
			tmp.to_csr();
			return tmp;
		}

		spmats as_csc() const {
			spmats tmp(*this);
			tmp.to_csc();
			return tmp;
		}

		std::vector<_T> mul_vec(const std::vector<_T>& x) const {
			if (column != size_to_index(x.size(), "spmats::mul_vec")) {
				vcp::throw_error<vcp::dimension_error>("spmats::mul_vec: dimension mismatch");
			}
			std::vector<_T> y(index_to_size(row, "spmats::mul_vec"), _T(0));
			mul_vec(x.data(), y.data());
			return y;
		}

		void mul_vec(const _T* x, _T* y) const {
			require_finalized("spmats::mul_vec");
			if (fmt == vcp::sparse_csr) {
				vcp::tcsrmv('N', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else if (fmt == vcp::sparse_csc) {
				vcp::tcscmv('N', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else {
				vcp::throw_error<vcp::state_error>("spmats::mul_vec: unsupported sparse format");
			}
		}

		std::vector<_T> trans_mul_vec(const std::vector<_T>& x) const {
			if (row != size_to_index(x.size(), "spmats::trans_mul_vec")) {
				vcp::throw_error<vcp::dimension_error>("spmats::trans_mul_vec: dimension mismatch");
			}
			std::vector<_T> y(index_to_size(column, "spmats::trans_mul_vec"), _T(0));
			trans_mul_vec(x.data(), y.data());
			return y;
		}

		void trans_mul_vec(const _T* x, _T* y) const {
			require_finalized("spmats::trans_mul_vec");
			if (fmt == vcp::sparse_csr) {
				vcp::tcsrmv('T', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else if (fmt == vcp::sparse_csc) {
				vcp::tcscmv('T', row, column, _T(1), outer.data(), inner.data(), value.data(), x, _T(0), y);
			}
			else {
				vcp::throw_error<vcp::state_error>("spmats::trans_mul_vec: unsupported sparse format");
			}
		}

		spmats transpose() const {
			spmats tmp(*this);
			tmp.ensure_coo_buffer();
			spmats out;
			out.resize(column, row);
			out.reserve(size_to_index(tmp.coo_value.size(), "spmats::transpose"));
			for (std::size_t k = 0; k < tmp.coo_value.size(); k++) {
				out.coo_row.push_back(tmp.coo_col[k]);
				out.coo_col.push_back(tmp.coo_row[k]);
				out.coo_value.push_back(tmp.coo_value[k]);
			}
			out.finalize();
			return out;
		}

		void assign_csr(const index_type rows, const index_type cols,
		                const std::vector<index_type>& row_ptr,
		                const std::vector<index_type>& col_ind,
		                const std::vector<_T>& val) {
			validate_csr(rows, cols, row_ptr, col_ind, val);
			row = rows;
			column = cols;
			outer = row_ptr;
			inner = col_ind;
			value = val;
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			fmt = vcp::sparse_csr;
			finalized = true;
			sorted = true;
			unique = true;
		}

		void assign_csc(const index_type rows, const index_type cols,
		                const std::vector<index_type>& col_ptr,
		                const std::vector<index_type>& row_ind,
		                const std::vector<_T>& val) {
			validate_csc(rows, cols, col_ptr, row_ind, val);
			row = rows;
			column = cols;
			outer = col_ptr;
			inner = row_ind;
			value = val;
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
			fmt = vcp::sparse_csc;
			finalized = true;
			sorted = true;
			unique = true;
		}

		const std::vector<index_type>& outer_index() const { return outer; }
		const std::vector<index_type>& inner_index() const { return inner; }
		const std::vector<_T>& values() const { return value; }
		const std::vector<index_type>& coo_rows() const { return coo_row; }
		const std::vector<index_type>& coo_columns() const { return coo_col; }
		const std::vector<_T>& coo_values() const { return coo_value; }

	protected:
		index_type row;
		index_type column;
		format_type fmt;
		bool finalized;
		bool sorted;
		bool unique;

		std::vector<index_type> outer;
		std::vector<index_type> inner;
		std::vector<_T> value;

		std::vector<index_type> coo_row;
		std::vector<index_type> coo_col;
		std::vector<_T> coo_value;

		void clear_storage() {
			outer.clear();
			inner.clear();
			value.clear();
			coo_row.clear();
			coo_col.clear();
			coo_value.clear();
		}

		void check_index(const index_type i, const index_type j, const char* routine) const {
			if (i < 0 || i >= row || j < 0 || j >= column) {
				vcp::throw_error<vcp::index_error>(routine, ": index out of range");
			}
		}

		void require_finalized(const char* routine) const {
			if (!finalized) {
				vcp::throw_error<vcp::state_error>(routine, ": matrix is not finalized");
			}
		}

		static std::size_t index_to_size(const index_type n, const char* routine) {
			if (n < 0) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": negative size");
			}
			return static_cast<std::size_t>(n);
		}

		static index_type size_to_index(const std::size_t n, const char* routine) {
			if (n > static_cast<std::size_t>((std::numeric_limits<index_type>::max)())) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": size exceeds index_type range");
			}
			return static_cast<index_type>(n);
		}

		static index_type checked_plus_one(const index_type n, const char* routine) {
			if (n < 0) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": negative size");
			}
			if (n == (std::numeric_limits<index_type>::max)()) {
				vcp::throw_error<vcp::invalid_argument>(routine, ": index_type overflow");
			}
			return static_cast<index_type>(n + 1);
		}

		void ensure_coo_buffer() {
			if (!finalized && fmt == vcp::sparse_coo) return;

			std::vector<index_type> nr;
			std::vector<index_type> nc;
			std::vector<_T> nv;
			const index_type n = size_to_index(value.size(), "spmats::ensure_coo_buffer");
			nr.resize(index_to_size(n, "spmats::ensure_coo_buffer"));
			nc.resize(index_to_size(n, "spmats::ensure_coo_buffer"));
			nv.resize(index_to_size(n, "spmats::ensure_coo_buffer"));

			if (fmt == vcp::sparse_csr) {
				vcp::tcsr_to_coo(row, column, n, outer.data(), inner.data(), value.data(), nr.data(), nc.data(), nv.data());
			}
			else if (fmt == vcp::sparse_csc) {
				index_type p = 0;
				for (index_type j = 0; j < column; j++) {
					for (index_type k = outer[j]; k < outer[j + 1]; k++) {
						nr[static_cast<std::size_t>(p)] = inner[k];
						nc[static_cast<std::size_t>(p)] = j;
						nv[static_cast<std::size_t>(p)] = value[k];
						p++;
					}
				}
			}
			else {
				vcp::throw_error<vcp::state_error>("spmats::ensure_coo_buffer: unsupported sparse format");
			}
			coo_row.swap(nr);
			coo_col.swap(nc);
			coo_value.swap(nv);
			outer.clear();
			inner.clear();
			value.clear();
			fmt = vcp::sparse_coo;
			finalized = false;
			sorted = true;
			unique = true;
		}

		static void validate_csr(const index_type rows, const index_type cols,
		                         const std::vector<index_type>& row_ptr,
		                         const std::vector<index_type>& col_ind,
		                         const std::vector<_T>& val) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: negative size");
			}
			if (row_ptr.size() != index_to_size(checked_plus_one(rows, "spmats::assign_csr"), "spmats::assign_csr") || col_ind.size() != val.size()) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: invalid array size");
			}
			if (row_ptr.empty() || row_ptr[0] != 0 || row_ptr.back() != size_to_index(col_ind.size(), "spmats::assign_csr")) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: invalid row_ptr boundary");
			}
			for (index_type i = 0; i < rows; i++) {
				if (row_ptr[i] > row_ptr[i + 1]) {
					vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: row_ptr must be monotone");
				}
				index_type previous = -1;
				for (index_type p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
					if (col_ind[p] < 0 || col_ind[p] >= cols) {
						vcp::throw_error<vcp::index_error>("spmats::assign_csr: col_ind out of range");
					}
					if (col_ind[p] <= previous) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: columns must be sorted and unique in each row");
					}
					if (val[p] == _T(0)) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csr: explicit zero is not allowed");
					}
					previous = col_ind[p];
				}
			}
		}

		static void validate_csc(const index_type rows, const index_type cols,
		                         const std::vector<index_type>& col_ptr,
		                         const std::vector<index_type>& row_ind,
		                         const std::vector<_T>& val) {
			if (rows < 0 || cols < 0) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: negative size");
			}
			if (col_ptr.size() != index_to_size(checked_plus_one(cols, "spmats::assign_csc"), "spmats::assign_csc") || row_ind.size() != val.size()) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: invalid array size");
			}
			if (col_ptr.empty() || col_ptr[0] != 0 || col_ptr.back() != size_to_index(row_ind.size(), "spmats::assign_csc")) {
				vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: invalid col_ptr boundary");
			}
			for (index_type j = 0; j < cols; j++) {
				if (col_ptr[j] > col_ptr[j + 1]) {
					vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: col_ptr must be monotone");
				}
				index_type previous = -1;
				for (index_type p = col_ptr[j]; p < col_ptr[j + 1]; p++) {
					if (row_ind[p] < 0 || row_ind[p] >= rows) {
						vcp::throw_error<vcp::index_error>("spmats::assign_csc: row_ind out of range");
					}
					if (row_ind[p] <= previous) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: rows must be sorted and unique in each column");
					}
					if (val[p] == _T(0)) {
						vcp::throw_error<vcp::invalid_argument>("spmats::assign_csc: explicit zero is not allowed");
					}
					previous = row_ind[p];
				}
			}
		}
		// ------------------------------------------------------------------
		// Symmetry check (used by eigs dispatch and lss)
		// ------------------------------------------------------------------
	public:
		typedef typename vcp::tsparse_scalar::real_type<_T>::type scalar_real_type;

		bool is_symmetric() const {
			return is_symmetric(vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(12));
		}

		bool is_symmetric(const scalar_real_type& tol) const {
			if (tol <= scalar_real_type(0))
				vcp::throw_error<vcp::invalid_argument>("spmats::is_symmetric: tol must be positive");
			if (row != column) return false;
			spmats A = this->as_csr();
			const std::vector<_Index>& outer = A.outer;
			const std::vector<_Index>& inner = A.inner;
			const std::vector<_T>& val = A.value;
			for (_Index i = 0; i < A.row; i++) {
				for (_Index p = outer[static_cast<std::size_t>(i)];
				     p < outer[static_cast<std::size_t>(i + 1)]; p++) {
					const _Index j = inner[static_cast<std::size_t>(p)];
					if (i == j) continue;
					const _Index first = outer[static_cast<std::size_t>(j)];
					const _Index last = outer[static_cast<std::size_t>(j + 1)];
					const typename std::vector<_Index>::const_iterator begin = inner.begin() + first;
					const typename std::vector<_Index>::const_iterator end = inner.begin() + last;
					typename std::vector<_Index>::const_iterator it = std::lower_bound(begin, end, i);
					_T mirrored = _T(0);
					if (it != end && *it == i)
						mirrored = val[static_cast<std::size_t>(it - inner.begin())];
					const scalar_real_type diff = vcp::tsparse_scalar::abs_value(
						val[static_cast<std::size_t>(p)] - mirrored);
					if (diff > tol) return false;
				}
			}
			return true;
		}

		// ------------------------------------------------------------------
		// Phase 7.7: policy methods for to_dense / is_symmetric forwarding
		// spmatrix<T,P> delegates A.to_dense() / A.is_symmetric() here.
		// Custom policy P can override these to change dense conversion or
		// symmetry-check semantics (e.g. verified / interval enclosure).
		// ------------------------------------------------------------------
		typedef std::vector<std::vector<_T> > dense_matrix_type;

		// policy_to_dense: CSR traversal → row-major dense 2-D vector.
		// Missing entries become T(0).  Behaviour matches the pre-7.7
		// spmatrix::to_dense() implementation.
		template <class Matrix>
		dense_matrix_type policy_to_dense(const Matrix& A) const {
			dense_matrix_type dense(
				static_cast<std::size_t>(A.rowsize()),
				std::vector<_T>(static_cast<std::size_t>(A.columnsize()), _T(0)));
			Matrix Acsr = A.as_csr();
			const std::vector<_Index>& outerv = Acsr.outer_index();
			const std::vector<_Index>& innerv = Acsr.inner_index();
			const std::vector<_T>& val = Acsr.values();
			for (_Index i = 0; i < Acsr.rowsize(); i++) {
				for (_Index p = outerv[static_cast<std::size_t>(i)];
				     p < outerv[static_cast<std::size_t>(i + 1)]; p++) {
					dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(innerv[static_cast<std::size_t>(p)])] =
						val[static_cast<std::size_t>(p)];
				}
			}
			return dense;
		}

		// policy_is_symmetric (no-tol): delegates to tolerance overload.
		template <class Matrix>
		bool policy_is_symmetric(const Matrix& A) const {
			return policy_is_symmetric(A, vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(12));
		}

		// policy_is_symmetric (with tol): complex-symmetric check, NOT Hermitian.
		// For complex T, compares A(i,j) with A(j,i), not conj(A(j,i)).
		// Behaviour matches the pre-7.7 spmatrix::is_symmetric(tol) implementation.
		template <class Matrix>
		bool policy_is_symmetric(const Matrix& A, const scalar_real_type& tol) const {
			if (tol <= scalar_real_type(0))
				vcp::throw_error<vcp::invalid_argument>("spmats::policy_is_symmetric: tol must be positive");
			if (A.rowsize() != A.columnsize()) return false;
			Matrix Acsr = A.as_csr();
			const std::vector<_Index>& outerv = Acsr.outer_index();
			const std::vector<_Index>& innerv = Acsr.inner_index();
			const std::vector<_T>& val = Acsr.values();
			for (_Index i = 0; i < Acsr.rowsize(); i++) {
				for (_Index p = outerv[static_cast<std::size_t>(i)];
				     p < outerv[static_cast<std::size_t>(i + 1)]; p++) {
					const _Index j = innerv[static_cast<std::size_t>(p)];
					if (i == j) continue;
					const _Index first = outerv[static_cast<std::size_t>(j)];
					const _Index last  = outerv[static_cast<std::size_t>(j + 1)];
					const typename std::vector<_Index>::const_iterator begin = innerv.begin() + first;
					const typename std::vector<_Index>::const_iterator end   = innerv.begin() + last;
					typename std::vector<_Index>::const_iterator it = std::lower_bound(begin, end, i);
					_T mirrored = _T(0);
					if (it != end && *it == i)
						mirrored = val[static_cast<std::size_t>(it - innerv.begin())];
					const scalar_real_type diff = vcp::tsparse_scalar::abs_value(
						val[static_cast<std::size_t>(p)] - mirrored);
					if (diff > tol) return false;
				}
			}
			return true;
		}

		// ------------------------------------------------------------------
		// Private helpers for policy methods
		// ------------------------------------------------------------------
	private:
		// solve_* private helpers (called from policy_lss_with_info)
		linear_solve_result<_T> policy_solve_jacobi_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol, bool use_relative) const;
		linear_solve_result<_T> policy_solve_gauss_seidel_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol, bool use_relative) const;
		linear_solve_result<_T> policy_solve_cg_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol,
			bool check_symmetric, preconditioner_type prec, bool use_relative) const;
		linear_solve_result<_T> policy_solve_bicgstab_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol,
			bool use_relative, preconditioner_type prec) const;
		linear_solve_result<_T> policy_solve_gmres_with_info_(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			std::size_t max_iter, const scalar_real_type& tol,
			std::size_t restart, bool use_relative, preconditioner_type prec) const;

	public:
		// ------------------------------------------------------------------
		// Policy methods: arithmetic
		// ------------------------------------------------------------------
		spmats<_T,_Index> policy_add(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;
		spmats<_T,_Index> policy_sub(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;
		spmats<_T,_Index> policy_mul(const spmats<_T,_Index>& A, const spmats<_T,_Index>& B) const;
		std::vector<_T> policy_mul_vec(const spmats<_T,_Index>& A, const std::vector<_T>& x) const;
		std::vector<_T> policy_left_mul_vec(const std::vector<_T>& x, const spmats<_T,_Index>& A) const;
		spmats<_T,_Index> policy_scalar_mul(const spmats<_T,_Index>& A, const _T& alpha) const;
		spmats<_T,_Index> policy_scalar_div(const spmats<_T,_Index>& A, const _T& alpha) const;
		spmats<_T,_Index> policy_neg(const spmats<_T,_Index>& A) const;

		// ------------------------------------------------------------------
		// Policy methods: linear system solve
		// ------------------------------------------------------------------
		linear_solve_result<_T> policy_lss_with_info(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			const linear_solve_options<_T>& opt) const;
		std::vector<_T> policy_lss(
			const spmats<_T,_Index>& A, const std::vector<_T>& b,
			const linear_solve_options<_T>& opt) const;

		// ------------------------------------------------------------------
		// Policy methods: eigenvalue _with_info (non-throwing, definitions in spmats_eigs.hpp)
		// ------------------------------------------------------------------
		eig_result<_T> policy_eigs_with_info(
			const spmats<_T,_Index>& A, std::size_t k,
			const eig_options<_T>& opt) const;

		template <class Prec>
		eig_result<_T> policy_eigs_with_info(
			const spmats<_T,_Index>& A, std::size_t k,
			const eig_options<_T>& opt, const Prec& M) const;

		eig_result<_T> policy_generalized_eigs_with_info(
			const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;

		template <class Prec>
		eig_result<_T> policy_generalized_eigs_with_info(
			const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt, const Prec& M) const;

		// ------------------------------------------------------------------
		// Policy methods: eigenvalue strict (throwing on failure)
		// These own all convergence/count checking; spmatrix.hpp does none.
		// ------------------------------------------------------------------
		eig_result<_T> policy_eig(
			const spmats<_T,_Index>& A,
			const eig_options<_T>& opt) const;

		std::vector<_T> policy_eigs(
			const spmats<_T,_Index>& A, std::size_t k,
			const eig_options<_T>& opt) const;

		template <class Prec>
		std::vector<_T> policy_eigs(
			const spmats<_T,_Index>& A, std::size_t k,
			const eig_options<_T>& opt, const Prec& M) const;

		eig_result<_T> policy_generalized_eig(
			const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;

		std::vector<_T> policy_generalized_eigs(
			const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt) const;

		template <class Prec>
		std::vector<_T> policy_generalized_eigs(
			const spmats<_T,_Index>& A, const spmats<_T,_Index>& B,
			std::size_t k, const eig_options<_T>& opt, const Prec& M) const;
	};
}

// Include policy method implementations (out-of-line definitions)
#include <vcp/spmats_product.hpp>
#include <vcp/spmats_lss.hpp>
#include <vcp/spmats_eigs.hpp>

#endif
