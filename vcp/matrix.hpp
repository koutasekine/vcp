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

#ifndef VCP_MATRIX_HPP
#define VCP_MATRIX_HPP

#include <initializer_list>
#include <type_traits>
#include <utility>

#include <vcp/mats.hpp>
#include <vcp/vcp_converter.hpp>

namespace vcp {
	// SPC-2: forward declaration only. matrix.hpp must NOT include
	// spmatrix.hpp (design SPC-2_design.md S2-2): every use of spmatrix in
	// the sparse-to-dense operator= below is a dependent expression whose
	// names resolve at instantiation time, and any translation unit that
	// instantiates it necessarily includes spmatrix.hpp itself.  No default
	// argument here -- spmatrix.hpp declares it.
	template <typename _T, class _P> class spmatrix;

	// SUB-1: block-write proxy (defined after matrix below).
	template <typename _T, class _P> class matrix_block;

	template <typename _T, class _P = mats< _T >> class matrix : protected _P {
	public:
		matrix() {
			this->row = 0;
			this->column = 0;
			this->n = 0;
			this->type = 'N';
		}
		~matrix() = default;
		matrix(const matrix&) = default;
		matrix(matrix&&) = default;
		matrix& operator=(const matrix& A) = default;
		matrix& operator=(matrix&& A) = default;

		// (1) Policy conversion constructor: T same, P different
		template <class _P2,
		          typename std::enable_if<!std::is_same<_P, _P2>::value, int>::type = 0>
		matrix(const matrix<_T, _P2>& A) {
			this->row = 0; this->column = 0; this->n = 0; this->type = 'N';
			this->zeros(A.rowsize(), A.columnsize());
			for (int i = 0; i < A.rowsize(); i++)
				for (int j = 0; j < A.columnsize(); j++)
					(*this)(i, j) = A(i, j);
		}

		// (2) Type conversion constructor: T different, P arbitrary
		template <typename _T2, class _P2,
		          typename std::enable_if<!std::is_same<_T, _T2>::value, int>::type = 0>
		matrix(const matrix<_T2, _P2>& A) {
			this->row = 0; this->column = 0; this->n = 0; this->type = 'N';
			this->zeros(A.rowsize(), A.columnsize());
			for (int i = 0; i < A.rowsize(); i++)
				for (int j = 0; j < A.columnsize(); j++)
					vcp::convert(A(i, j), (*this)(i, j));
		}

		// (3) Policy conversion assignment: T same, P different
		template <class _P2,
		          typename std::enable_if<!std::is_same<_P, _P2>::value, int>::type = 0>
		matrix<_T, _P>& operator=(const matrix<_T, _P2>& A) {
			this->zeros(A.rowsize(), A.columnsize());
			for (int i = 0; i < A.rowsize(); i++)
				for (int j = 0; j < A.columnsize(); j++)
					(*this)(i, j) = A(i, j);
			return *this;
		}

		// (4) Type conversion assignment: T different, P arbitrary
		template <typename _T2, class _P2,
		          typename std::enable_if<!std::is_same<_T, _T2>::value, int>::type = 0>
		matrix<_T, _P>& operator=(const matrix<_T2, _P2>& A) {
			this->zeros(A.rowsize(), A.columnsize());
			for (int i = 0; i < A.rowsize(); i++)
				for (int j = 0; j < A.columnsize(); j++)
					vcp::convert(A(i, j), (*this)(i, j));
			return *this;
		}

	private:
		// SPC-2: value-conversion helper for the sparse-to-dense operator=
		// below (directive §4).  Same scalar type: plain assignment (does
		// not rely on an identity vcp::convert overload).  Different
		// scalar type: delegate to the vcp::convert scalar layer
		// (point->point widening exact / narrowing nearest, ->interval
		// inclusion-preserving, interval->point mid).
		template <typename _S>
		static void convert_value_(const _S& x, _S& y) { y = x; }
		template <typename _S2, typename _S,
		          typename std::enable_if<!std::is_same<_S, _S2>::value, int>::type = 0>
		static void convert_value_(const _S2& x, _S& y) { vcp::convert(x, y); }

		// SUB-1: selector validation/normalization for the write proxy
		// (matrix_block).  Conditions, order and message text mirror
		// mats::submat exactly (the read path keeps its own identical
		// checks inside mats::submat itself, which is deliberately left
		// untouched).  Validation is completed here, at proxy construction,
		// before any write can begin.
		static void submat_check_sizes_(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) {
			if (list1.size() > 3 || list2.size() > 3) {
				vcp::throw_error<vcp::index_error>(
					"submat: invalid selector size: ", list1.size(), ", ", list2.size());
			}
		}
		static void submat_normalize_(const std::initializer_list<int>& list, const int dim, const bool is_row,
		                              int& start, int& stride, int& count, bool& full) {
			const std::vector<int> l = list;
			full = false;
			if (l.size() == 0) {
				start = 0;
				stride = 1;
				count = dim;
				full = true;
			}
			else if (l.size() == 1) {
				if (l[0] < 0 || l[0] >= dim) {
					if (is_row) {
						vcp::throw_error<vcp::index_error>("submat: row index out of range: ", l[0]);
					}
					vcp::throw_error<vcp::index_error>("submat: column index out of range: ", l[0]);
				}
				start = l[0];
				stride = 1;
				count = 1;
			}
			else if (l.size() == 2) {
				if (l[0] > l[1] || l[0] < 0 || l[1] >= dim) {
					if (is_row) {
						vcp::throw_error<vcp::index_error>(
							"submat: invalid row range: ", l[0], ":", l[1]);
					}
					vcp::throw_error<vcp::index_error>(
						"submat: invalid column range: ", l[0], ":", l[1]);
				}
				start = l[0];
				stride = 1;
				count = l[1] - l[0] + 1;
			}
			else {
				if (l[0] > l[2] || l[0] < 0 || l[1] < 1 || l[2] >= dim) {
					if (is_row) {
						vcp::throw_error<vcp::index_error>(
							"submat: invalid row range: ", l[0], ":", l[1], ":", l[2]);
					}
					vcp::throw_error<vcp::index_error>(
						"submat: invalid column range: ", l[0], ":", l[1], ":", l[2]);
				}
				int k = 0;
				for (int i = l[0]; i <= l[2]; i += l[1]) {
					k++;
				}
				start = l[0];
				stride = l[1];
				count = k;
			}
		}

	public:
		// SPC-2: sparse-to-dense conversion assignment
		// (design SPC-2_design.md v1.0 S2-1..S2-3).
		// Densifies A: O(rowsize*columnsize) memory.  Only fires on an
		// explicit assignment statement; does not participate in implicit
		// conversions (no converting constructor is provided).
		// A is read through an as_csr()/as_csc() copy only (finalize
		// policy category 4): the source is never modified, a finalized
		// CSC source is read as CSC to avoid the O(nnz log nnz) re-sort,
		// and an unfinalized (COO) source is normalized on the copy.
		// Complexity O(rows*cols + nnz).
		template <typename _T2, class _P2>
		matrix<_T, _P>& operator=(const spmatrix<_T2, _P2>& A) {
			typedef typename spmatrix<_T2, _P2>::index_type src_index_type;
			typedef typename spmatrix<_T2, _P2>::format_type src_format_type;
			// C++11 [dcl.enum]/11: an enumerator is reachable through the
			// scope of its enumeration type, so src_format_type::sparse_csc
			// is the same enumerator as vcp::sparse_csc but spelled as a
			// dependent name -- matrix.hpp stays compilable on its own
			// (vcp::sparse_csc itself is not declared here).
			const bool use_csc = A.is_finalized()
			                     && A.format() == src_format_type::sparse_csc;
			const spmatrix<_T2, _P2> src = use_csc ? A.as_csc() : A.as_csr();
			const int rows = static_cast<int>(src.rowsize());
			const int cols = static_cast<int>(src.columnsize());
			const std::vector<src_index_type>& src_outer = src.outer_index();
			const std::vector<src_index_type>& src_inner = src.inner_index();
			const std::vector<_T2>& src_value = src.values();
			this->zeros(rows, cols);
			const int nouter = use_csc ? cols : rows;
			for (int i = 0; i < nouter; i++) {
				const std::size_t ui = static_cast<std::size_t>(i);
				for (src_index_type k = src_outer[ui]; k < src_outer[ui + 1]; k++) {
					const std::size_t kk = static_cast<std::size_t>(k);
					const int q = static_cast<int>(src_inner[kk]);
					if (use_csc) {
						convert_value_(src_value[kk], (*this)(q, i));
					}
					else {
						convert_value_(src_value[kk], (*this)(i, q));
					}
				}
			}
			return *this;
		}

		_T& operator () (const int i) {
			return this->v[i];
		}
		_T operator () (const int i) const {
			return this->v[i];
		}
		_T& operator () (const int i, const int j) {
			if (this->type == 'R') {
				return this->v[j];
			}
			else {
				return this->v[i + this->row*j];
			}
		}
		_T operator () (const int i, const int j)const {
			if (this->type == 'R') {
				return this->v[j];
			}
			else {
				return this->v[i + this->row*j];
			}
		}
	/*
		_T& operator [] (const int i) {
			return this->v[i];
		}
		_T operator [] (const int i) const {
			return this->v[i];
		}
		*/
		
		matrix< _T, _P > submatrix(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
			matrix< _T, _P > A;
			this->submat(A, list1, list2);
			return A;
		}

		// SUB-1: const operator() sugar -- pure forwarding to submatrix
		// (single implementation).  Selector grammar is mats::submat's
		// (0-based): {} whole axis, {i} single index, {a,b} CLOSED range
		// a..b (both ends included -- NOT half-open), {a,s,b} stride
		// a, a+s, ... while <= b; invalid selectors throw vcp::index_error.
		// Return-type rule (S11): bare ints only = element (_T, A(i,j));
		// ANY braced selector = submatrix: A({1},{2}), A(1,{2}), A({1},2)
		// are 1x1 matrices, A(1,2) is the element.
		matrix< _T, _P > operator () (const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
			return this->submatrix(list1, list2);
		}
		// mixed forms (S11): the int is wrapped as the single selector {i}
		matrix< _T, _P > operator () (const int i, const std::initializer_list<int>& list2) const {
			return this->submatrix({i}, list2);
		}
		matrix< _T, _P > operator () (const std::initializer_list<int>& list1, const int j) const {
			return this->submatrix(list1, {j});
		}

		// SUB-1 Phase W: the non-const forms return the block-write proxy
		// (matrix_block, defined after this class): A({..},{..}) = B
		// (same-type matrix, block assignment) or = scalar (fill).
		// Selector validation happens HERE, at proxy construction.
		// See matrix_block for the semantics and the lifetime warning
		// (`auto x = A({..},{..});` captures the proxy, not a matrix).
		matrix_block< _T, _P > operator () (const std::initializer_list<int>& list1, const std::initializer_list<int>& list2);
		matrix_block< _T, _P > operator () (const int i, const std::initializer_list<int>& list2);
		matrix_block< _T, _P > operator () (const std::initializer_list<int>& list1, const int j);

		int elementsize()const {
			if (this->n > static_cast<vcp::index_t>(2147483647)) {
				vcp::throw_error<vcp::dimension_error>(
					"elementsize(): n exceeds INT_MAX; use within-INT_MAX matrices or await 64-bit accessor (n = ", this->n, ")");
			}
			return static_cast<int>(this->n);
		}
		int columnsize()const { return static_cast<int>(this->column); }
		int rowsize()const { return static_cast<int>(this->row); }
		char matstype()const {
			return this->type;
		}
		const std::vector< _T >& vecpointer()const {
			return this->v;
		}

		_T* data() {
			return this->v.data();
		}
		const _T* data() const {
			return this->v.data();
		}

		void eye(const int r) { _P::eye(r); }
		void ones(const int i) { _P::ones(i); }
		void ones(const int r, const int c) { _P::ones(r, c); }
		void zeros(const int i) { _P::zeros(i); }
		void zeros(const int r, const int c) { _P::zeros(r, c); }
		void rand(const int i) { _P::rand(i); }
		void rand(const int r, const int c) { _P::rand(r, c); }
		void resize(const int i, const int j) { _P::resize(i, j); }
		void clear() { _P::clear(); }

		//***************** Operator Overload *****************//
		friend matrix< _T, _P > operator+(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			C = A;
			C.addmm(B);
			return C;
		}
		friend matrix< _T, _P > operator+(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			A.addmm(B);
			return std::move(A);
		}
		friend matrix< _T, _P > operator+(const matrix< _T, _P >& A, matrix< _T, _P >&& B) {
			B.addmm(A);
			return std::move(B);
		}
		friend matrix< _T, _P > operator+(matrix< _T, _P >&& A, matrix< _T, _P >&& B) {
			A.addmm(B);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.addsm(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.addsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.addms(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator+(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.addms(Ta);
			return std::move(B);
		}
		friend matrix< _T, _P > operator+(const matrix< _T, _P >& A) {
			//		A.plusm();
			return A;
		}
		friend matrix< _T, _P > operator+(matrix< _T, _P >&& A) {
			//		A.plusm();
			return std::move(A);
		}

		friend matrix< _T, _P >& operator+=(matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			A.addmm(B);
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator+=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.addms(Ta);
			return A;
		}

		friend matrix< _T, _P > operator-(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			C = A;
			C.subsmmA(B);
			return C;
		}
		friend matrix< _T, _P > operator-(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		friend matrix< _T, _P > operator-(const matrix< _T, _P >& A, matrix< _T, _P >&& B) {
			B.subsmmB(A);
			return std::move(B);
		}
		friend matrix< _T, _P > operator-(matrix< _T, _P >&& A, matrix< _T, _P >&& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.subssm(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.subssm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.subsms(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator-(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.subsms(Ta);
			return std::move(B);
		}
		friend matrix< _T, _P > operator-(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.minusm();
			return C;
		}
		friend matrix< _T, _P > operator-(matrix< _T, _P >&& A) {
			A.minusm();
			return std::move(A);
		}

		friend matrix< _T, _P >& operator-=(matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			A.subsmmA(B);
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator-=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.subsms(Ta);
			return A;
		}

		friend matrix< _T, _P > operator*(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			return C;
		}
		friend matrix< _T, _P > operator*(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			A = std::move(C);
			return std::move(A);
		}
		friend matrix< _T, _P > operator*(const matrix< _T, _P >& A, matrix< _T, _P >&& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			B = std::move(C);
			return std::move(B);
		}
		friend matrix< _T, _P > operator*(matrix< _T, _P >&& A, matrix< _T, _P >&& B) {
			matrix< _T, _P > C;
			A.mulmm(B, C);
			A = std::move(C);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.mulsm(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.mulsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.mulms(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator*(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.mulms(Ta);
			return std::move(B);
		}

		friend matrix< _T, _P >& operator*=(matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			A = A * B;
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator*=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.mulms(Ta);
			return A;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.divsm(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.divsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.divms(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type operator/(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.divms(Ta);
			return std::move(B);
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P >& >::type operator/=(matrix< _T, _P >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.divms(Ta);
			return A;
		}

		friend mbool operator>(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			mbool C;
			A.gt(B, C);
			return C;
		}
		friend mbool operator>=(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			mbool C;
			A.ge(B, C);
			return C;
		}
		friend mbool operator<(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			mbool C;
			A.lt(B, C);
			return C;
		}
		friend mbool operator<=(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			mbool C;
			A.le(B, C);
			return C;
		}
		friend mbool operator==(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			mbool C;
			A.eq(B, C);
			return C;
		}
		friend mbool operator!=(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			mbool C;
			A.neq(B, C);
			return C;
		}

		friend matrix< _T, _P > pow(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			C = A;
			C.powmmA(B);
			return C;
		}
		friend matrix< _T, _P > pow(matrix< _T, _P >&& A, const matrix< _T, _P >& B) {
			A.powmmA(B);
			return std::move(A);
		}
		friend matrix< _T, _P > pow(const matrix< _T, _P > A, matrix< _T, _P >&& B) {
			B.powmmB(A);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(const matrix< _T, _P >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.powms(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(matrix< _T, _P >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.powms(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(const _Tm a, const matrix< _T, _P >& B) {
			_T Ta = _T(a);
			matrix< _T, _P > C;
			C = B;
			C.powsm(Ta);
			return C;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, _P > >::type pow(const _Tm a, matrix< _T, _P >&& B) {
			_T Ta = _T(a);
			B.powsm(Ta);
			return std::move(B);
		}

		//***************** Special Matrix Arithemetic *****************//
		friend matrix< _T, _P > ltransmul(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			A.mulltmm(C);
			return C;
		}
		friend matrix< _T, _P > ltransmul(matrix< _T, _P >&& A) {
			matrix< _T, _P > C;
			A.mulltmm(C);
			A = std::move(C);
			return std::move(A);
		}

		//***************** Math functions *****************//
		friend matrix< _T, _P > abs(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.abs();
			return C;
		}
		friend matrix< _T, _P > abs(matrix< _T, _P >&& A) {
			A.abs();
			return std::move(A);
		}
		friend matrix< _T, _P > sqrt(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.sqrt();
			return C;
		}
		friend matrix< _T, _P > sqrt(matrix< _T, _P >&& A) {
			A.sqrt();
			return std::move(A);
		}
		friend matrix< _T, _P > sin(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.sin();
			return C;
		}
		friend matrix< _T, _P > sin(matrix< _T, _P >&& A) {
			A.sin();
			return std::move(A);
		}
		friend matrix< _T, _P > cos(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.cos();
			return C;
		}
		friend matrix< _T, _P > cos(matrix< _T, _P >&& A) {
			A.cos();
			return std::move(A);
		}
		friend matrix< _T, _P > exp(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.exp();
			return C;
		}
		friend matrix< _T, _P > exp(matrix< _T, _P >&& A) {
			A.exp();
			return std::move(A);
		}
		friend matrix< _T, _P > log(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			C = A;
			C.log();
			return C;
		}
		friend matrix< _T, _P > log(matrix< _T, _P >&& A) {
			A.log();
			return std::move(A);
		}

		//************* matlab like functions *************//	
		friend matrix< _T, _P > sum(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.sum(c);
			return c;
		}
		friend matrix< _T, _P > diag(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			A.diag(C);
			return C;
		}
		friend matrix< _T, _P > transpose(const matrix< _T, _P >& A) {
			matrix< _T, _P > C;
			A.transpose(C);
			return C;
		}
		friend matrix< _T, _P > max(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.max(c);
			return c;
		}
		friend matrix< _T, _P > min(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.min(c);
			return c;
		}
		friend matrix< _T, _P > normone(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.normone(c);
			return c;
		}
		friend matrix< _T, _P > normtwo(const matrix< _T, _P >& A) {
			matrix< _T, _P > c = A;
			c.normtwo();
			return c;
		}
		friend matrix< _T, _P > normtwo(matrix< _T, _P >&& A) {
			A.normtwo();
			return std::move(A);
		}
		friend matrix< _T, _P > norminf(const matrix< _T, _P >& A) {
			matrix< _T, _P > c;
			A.norminf(c);
			return c;
		}

		friend matrix< _T, _P > tril(const matrix< _T, _P >& A) {
			matrix< _T, _P > c = A;
			c.tril();
			return c;
		}
		friend matrix< _T, _P > tril(matrix< _T, _P >&& A) {
			A.tril();
			return std::move(A);
		}
		friend matrix< _T, _P > triu(const matrix< _T, _P >& A) {
			matrix< _T, _P > c = A;
			c.triu();
			return c;
		}
		friend matrix< _T, _P > triu(matrix< _T, _P >&& A) {
			A.triu();
			return std::move(A);
		}

		friend matrix< _T, _P > SymTridiagonalization(const matrix< _T, _P >& A) {
			matrix< _T, _P > AA;
			AA = A;
			AA.SymTridiagonalization();
			return AA;
		}
		friend matrix< _T, _P > SymTridiagonalization(matrix< _T, _P >&& A) {
			A.SymTridiagonalization();
			return A;
		}

		friend matrix< _T, _P > HessenbergTrans(const matrix< _T, _P >& A) {
			matrix< _T, _P > AA;
			AA = A;
			AA.HessenbergTrans();
			return AA;
		}
		friend matrix< _T, _P > HessenbergTrans(matrix< _T, _P >&& A) {
			A.HessenbergTrans();
			return A;
		}

		//  L = tril(ipiv*A) + eye
		//  U = triu(ipiv*A) + eye
		friend void lu(matrix< _T, _P >& A, matrix< int >& ipiv) {
			A.ludecomposition(ipiv);
		}
		// [A, ipiv] = lu(A)
		//  L = tril(ipiv*LU) + eye
		//  U = triu(ipiv*LU) + eye
		friend void lu(const matrix< _T, _P >& A, matrix< _T, _P >& LU, matrix< int >& ipiv) {
			LU = A;
			LU.ludecomposition(ipiv);
		}
		// ipiv * A = L * U
		friend void lu(const matrix< _T, _P >& A, matrix< _T, _P >& L, matrix< _T, _P >& U, matrix< int >& ipiv) {
			U = A;
			U.ludecomposition(ipiv);
			U.LUtoLandU(L, ipiv);
		}

		friend void Trilu(const matrix< _T, _P >& A, matrix< _T, _P >& LU, matrix< int >& ipiv) {
			LU = A;
			LU.TriLudecomposition(ipiv);
		}

		// [Q, R] = qr(A), A=Q*R
		friend void qr(const matrix< _T, _P >& A, matrix< _T, _P >& Q, matrix< _T, _P >& R) {
			R = A;
			R.Householder_qrdecomposition(Q);
		}
		// [Q, R] = qr(A), A=Q*R Hessenberg
		friend void qr_Hessenberg(const matrix< _T, _P >& A, matrix< _T, _P >& Q, matrix< _T, _P >& R) {
			R = A;
			R.Householder_qrdecomposition_Hessenberg(Q);
		}
		// [Q, R] = qr(A), A=Q*R Tridiagonal
		friend void qr_Tridiag(const matrix< _T, _P >& A, matrix< _T, _P >& Q, matrix< _T, _P >& R) {
			R = A;
			R.Householder_qrdecomposition_Tridiag(Q);
		}

		// x = A\b
		friend void lss(const matrix< _T, _P >& A, const matrix< _T, _P >& b, matrix< _T, _P >& x) {
			matrix< _T, _P > AA = A;
			matrix< _T, _P > bb = b;
			AA.linearsolve(bb, x);
		}
		friend matrix< _T, _P > lss(const matrix< _T, _P >& A, const matrix< _T, _P >& b) {
			matrix< _T, _P > AA = A;
			matrix< _T, _P > bb = b;
			matrix< _T, _P > x;
			AA.linearsolve(bb, x);
			return x;
		}
		friend matrix< _T, _P > lss(matrix< _T, _P >&& A, const matrix< _T, _P >& b) {
			matrix< _T, _P > bb = b;
			matrix< _T, _P > x;
			A.linearsolve(bb, x);
			return x;
		}
		friend matrix< _T, _P > lss(const matrix< _T, _P >& A, matrix< _T, _P >&& b) {
			matrix< _T, _P > AA = A;
			matrix< _T, _P > x;
			AA.linearsolve(b, x);
			return x;
		}
		friend matrix< _T, _P > lss(matrix< _T, _P >&& A, matrix< _T, _P >&& b) {
			matrix< _T, _P > x;
			A.linearsolve(b, x);
			return x;
		}
		// R = inv(A)
		friend matrix< _T, _P > inv(const matrix< _T, _P >& A) {
			matrix< _T, _P > R = A;
			R.inv();
			return R;
		}
		friend matrix< _T, _P > inv(matrix< _T, _P >&& A) {
			A.inv();
			return std::move(A);
		}
		friend matrix< _T, _P > Cholesky(const matrix< _T, _P >& A) {
			matrix< _T, _P > C = A;
			C.Cholesky();
			return C;
		}
		friend matrix< _T, _P > Cholesky(matrix< _T, _P >&& A) {
			A.Cholesky();
			return std::move(A);
		}
		friend void eigsym(const matrix< _T, _P >& A, matrix< _T, _P >& E, int itep = 1) {
			E = A;
			E.eigsym(itep);
		}
		friend void eigsym(const matrix< _T, _P >& A, matrix< _T, _P >& E, matrix< _T, _P >& V, int itep = 1) {
			E = A;
			E.eigsym(V, itep);
		}
		friend void eigsymge(const matrix< _T, _P >& A, const matrix< _T, _P >& B, matrix< _T, _P >& E, int itep = 1) {
			E = A;
			matrix< _T, _P > BD = B;
			E.eigsymge(BD, itep);
		}
		friend void eigsymge(const matrix< _T, _P >& A, const matrix< _T, _P >& B, matrix< _T, _P >& E, matrix< _T, _P >& V, int itep = 1) {
			E = A;
			matrix< _T, _P > BD = B;
			E.eigsymge(BD, V, itep);
		}

		//Matlab C = [A,B];
		friend matrix< _T, _P > horzcat(const matrix< _T, _P >& A) {
			return A;
		}
		friend matrix< _T, _P > horzcat(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.horzcat(B, C);
			return C;
		}
		template<typename... Args> friend matrix< _T, _P > horzcat(const matrix< _T, _P >& A, const matrix< _T, _P >& B, const Args&... args) {
			matrix< _T, _P > C;
			A.horzcat(B, C);
			return horzcat(C, args...);
		}
		
		//Matlab C = [A;B];
		friend matrix< _T, _P > vercat(const matrix< _T, _P >& A) {
			return A;
		}
		friend matrix< _T, _P > vercat(const matrix< _T, _P >& A, const matrix< _T, _P >& B) {
			matrix< _T, _P > C;
			A.vercat(B, C);
			return C;
		}
		template<typename... Args> friend matrix< _T, _P > vercat(const matrix< _T, _P >& A, const matrix< _T, _P >& B, const Args&... args) {
			matrix< _T, _P > C;
			A.vercat(B, C);
			return vercat(C, args...);
		}

		friend int length(matrix< _T, _P >& A) {
			return A.length();
		}

		//**************** display function ***************//
		friend std::ostream& operator<<(std::ostream& os, const matrix< _T, _P >& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, matrix< _T, _P >&& A) {
			return A.display(os);
		}
	};

	// SUB-1: block-write proxy returned by the non-const
	// matrix::operator()({..},{..}) and the mixed int/list forms.
	// Holds a reference to the parent and the NORMALIZED selectors only;
	// selector validation is completed at proxy construction
	// (vcp::index_error), and operator= validates the right-hand-side
	// dimensions BEFORE any element is written, so a throwing assignment
	// leaves the parent untouched.  Supported operations (SUB-1_design.md
	// v1.1 §3):
	//   proxy = matrix<_T,_P>  -- block assignment.  The right-hand side is
	//                             materialized into a temporary FIRST, so
	//                             aliased/overlapping assignments such as
	//                             A({0,2},{0,2}) = A({1,3},{1,3}) follow
	//                             the MATLAB semantics.
	//   proxy = scalar _T      -- fill the whole block with the value.
	//   matrix<_T,_P>(proxy)   -- read; delegates to matrix::submatrix
	//                             (single implementation of extraction).
	// LIFETIME WARNING: `auto x = A({..},{..});` captures the PROXY, not a
	// matrix -- it must not outlive the parent A.  Use an explicit
	// matrix<_T,_P> variable to take a copy of the block.
	template <typename _T, class _P> class matrix_block {
	public:
		matrix_block(matrix< _T, _P >& A,
		             const int r0, const int rs, const int rn, const bool rfull,
		             const int c0, const int cs, const int cn, const bool cfull)
			: A_(A), r0_(r0), rs_(rs), rn_(rn), rfull_(rfull),
			  c0_(c0), cs_(cs), cn_(cn), cfull_(cfull) {}

		matrix_block(const matrix_block&) = default;

		// block assignment (right-hand side materialized first, S6)
		matrix_block& operator=(const matrix< _T, _P >& B) {
			if (B.rowsize() != rn_ || B.columnsize() != cn_) {
				vcp::throw_error<vcp::dimension_error>(
					"matrix_block: block assignment size mismatch: ",
					B.rowsize(), "x", B.columnsize(), " != ", rn_, "x", cn_);
			}
			const matrix< _T, _P > tmp = B;
			for (int i = 0; i < rn_; i++) {
				for (int j = 0; j < cn_; j++) {
					A_(r0_ + i * rs_, c0_ + j * cs_) = tmp(i, j);
				}
			}
			return *this;
		}

		// scalar fill
		matrix_block& operator=(const _T& s) {
			for (int i = 0; i < rn_; i++) {
				for (int j = 0; j < cn_; j++) {
					A_(r0_ + i * rs_, c0_ + j * cs_) = s;
				}
			}
			return *this;
		}

		// proxy = proxy (e.g. A({0,2},{0,2}) = A({1,3},{1,3})): materialize
		// the right-hand block first, then block-assign.  Without this
		// overload the implicitly-deleted copy assignment would win the
		// overload resolution (same reason as spmats_element_proxy).
		matrix_block& operator=(const matrix_block& other) {
			return (*this) = static_cast<matrix< _T, _P > >(other);
		}

		// read conversion -- delegates to matrix::submatrix, rebuilding the
		// normalized selectors as {} / {start, stride, last} (identical
		// selections by construction).
		operator matrix< _T, _P >() const {
			if (rfull_ && cfull_) {
				return A_.submatrix({}, {});
			}
			if (rfull_) {
				return A_.submatrix({}, {c0_, cs_, c0_ + (cn_ - 1) * cs_});
			}
			if (cfull_) {
				return A_.submatrix({r0_, rs_, r0_ + (rn_ - 1) * rs_}, {});
			}
			return A_.submatrix({r0_, rs_, r0_ + (rn_ - 1) * rs_},
			                    {c0_, cs_, c0_ + (cn_ - 1) * cs_});
		}

	private:
		matrix< _T, _P >& A_;
		int r0_, rs_, rn_;
		bool rfull_;
		int c0_, cs_, cn_;
		bool cfull_;
	};

	// SUB-1: out-of-line definitions of the proxy-returning operator()
	// overloads (declared inside matrix; matrix_block must be complete
	// here).  The mixed forms wrap the int as the single selector {i}.
	template <typename _T, class _P>
	matrix_block< _T, _P > matrix< _T, _P >::operator () (const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) {
		submat_check_sizes_(list1, list2);
		int r0, rs, rn, c0, cs, cn;
		bool rfull, cfull;
		submat_normalize_(list1, this->row, true, r0, rs, rn, rfull);
		submat_normalize_(list2, this->column, false, c0, cs, cn, cfull);
		return matrix_block< _T, _P >(*this, r0, rs, rn, rfull, c0, cs, cn, cfull);
	}
	template <typename _T, class _P>
	matrix_block< _T, _P > matrix< _T, _P >::operator () (const int i, const std::initializer_list<int>& list2) {
		return (*this)({i}, list2);
	}
	template <typename _T, class _P>
	matrix_block< _T, _P > matrix< _T, _P >::operator () (const std::initializer_list<int>& list1, const int j) {
		return (*this)(list1, {j});
	}

	template <> class matrix< bool >{
	protected:
		vcp::index_t row;
		vcp::index_t column;
		vcp::index_t n;
		char type;      //'N':NULL  'S':Scala  'R' Row Vector 'C':Column Vector 'M':Matrix
		std::vector< bool > v;

		// A = and(A,B)
		void mats_and(const matrix< bool >& B) {
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>("&&: dimension mismatch");
			}
			if (this->type == 'S') {
				this->v[0] = this->v[0] && B.v[0];
				return;
			}
			else {
				for (vcp::index_t i = 0; i < n; i++) {
					this->v[i] = this->v[i] && B.v[i];
				}
				return;
			}
		}
		// A = or(A,B)
		void mats_or(const matrix< bool >& B) {
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>("||: dimension mismatch");
			}
			if (this->type == 'S') {
				this->v[0] = this->v[0] || B.v[0];
				return;
			}
			else {
				for (vcp::index_t i = 0; i < n; i++) {
					this->v[i] = this->v[i] || B.v[i];
				}
				return;
			}
		}
		int length()const {
			using std::max;
			return static_cast<int>(max(column, row));
		}
		std::ostream& display(std::ostream& os)const {
			if (type == 'S') {
				os << v[0] << "\n";
			}
			else if (type == 'C') {
				for (int i = 0; i <= row - 1; i++) {
					os << v[i] << "\n";
				}
			}
			else if (type == 'R') {
				for (int i = 0; i <= column - 1; i++) {
					os << v[i] << " ";
				}
				os << "\n";
			}
			else if (type == 'M') {
				for (int j = 0; j <= row - 1; j++) {
					for (vcp::index_t i = j; i <= row*(column - 1) + j; i = i + row) {
						os << v[i] << "  ";
					}
					os << "\n";
				}
			}
			else {
				os << "display error";
			}
			return os;
		}
	
		bool all() const {
			for (vcp::index_t i = 0; i < this->n; i++) {
				if (!this->v[i]) {
					return false;
				}
			}
			return true;
		}
		bool any() const {
			for (vcp::index_t i = 0; i < this->n; i++) {
				if (this->v[i]) {
					return true;
				}
			}
			return false;
		}
		bool none() const {
			for (vcp::index_t i = 0; i < this->n; i++) {
				if (this->v[i]) {
					return false;
				}
			}
			return true;
		}


	public:
		matrix() {
			this->row = 0;
			this->column = 0;
			this->n = 0;
			this->type = 'N';
		}
		matrix(const mbool& A) {
			this->row = A.row;
			this->column = A.column;
			this->n = A.n;
			this->type = A.type;
			this->v = A.v;
		}
		matrix(mbool&& A) {
			this->row = A.row;
			this->column = A.column;
			this->n = A.n;
			this->type = A.type;
			this->v = std::move(A.v);
		}
		~matrix() = default;
		matrix(const matrix&) = default;
		matrix(matrix&&) = default;
		matrix& operator=(const matrix& A) = default;
		matrix& operator=(matrix&& A) = default;
	
		std::vector< bool >::reference operator () (const int i) {
			return this->v[i];
		}
		std::vector< bool >::const_reference operator () (const int i) const {
			return this->v[i];
		}
		std::vector< bool >::reference operator () (const int i, const int j) {
			if (this->type == 'R') {
				return this->v[j];
			}
			else {
				return this->v[i + this->row*j];
			}
		}
		std::vector< bool >::const_reference operator () (const int i, const int j)const {
			if (this->type == 'R') {
				return this->v[j];
			}
			else {
				return this->v[i + this->row*j];
			}
		}

		int elementsize()const {
			if (this->n > static_cast<vcp::index_t>(2147483647)) {
				vcp::throw_error<vcp::dimension_error>(
					"elementsize(): n exceeds INT_MAX; use within-INT_MAX matrices or await 64-bit accessor (n = ", this->n, ")");
			}
			return static_cast<int>(this->n);
		}
		int columnsize()const { return static_cast<int>(this->column); }
		int rowsize()const { return static_cast<int>(this->row); }
		char matstype()const {
			return this->type;
		}
		const std::vector< bool >& vecpointer()const {
			return this->v;
		}
		friend int length(matrix < bool > & A) {
			return A.length();
		}
		
		void alltrue(const int i) {
			row = i;
			column = i;
			n = static_cast<vcp::index_t>(i) * i;
			if (i == 1) {
				type = 'S';
			}
			else {
				type = 'M';
			}
			v.resize(n);
			for (vcp::index_t j = 0; j < n; j++) {
				v[j] = true;
			}
		}
		void alltrue(const int r, const int c) {
			row = r;
			column = c;
			n = row*column;
			if (row == 1 && column == 1) {
				type = 'S';
			}
			else if (column == 1) {
				type = 'C';
			}
			else if (row == 1) {
				type = 'R';
			}
			else if (row > 1 && column > 1) {
				type = 'M';
			}
			v.resize(n);
			for (int j = 0; j < column; j++) {
				for (int i = 0; i < row; i++) {
					v[i + row*j] = true;
				}
			}
		}
		void set(const int i) {
			this->alltrue(i);
		}
		void set(const int r, const int c) {
			this->alltrue(r, c);
		}

		void allfalse(const int i) {
			row = i;
			column = i;
			n = static_cast<vcp::index_t>(i) * i;
			if (i == 1) {
				type = 'S';
			}
			else {
				type = 'M';
			}
			v.resize(n);
			for (vcp::index_t j = 0; j < n; j++) {
				v[j] = false;
			}
		}
		void allfalse(const int r, const int c) {
			row = r;
			column = c;
			n = row*column;
			if (row == 1 && column == 1) {
				type = 'S';
			}
			else if (column == 1) {
				type = 'C';
			}
			else if (row == 1) {
				type = 'R';
			}
			else if (row > 1 && column > 1) {
				type = 'M';
			}
			v.resize(n);
			for (int j = 0; j < column; j++) {
				for (int i = 0; i < row; i++) {
					v[i + row*j] = false;
				}
			}
		}
		void reset(const int i) {
			this->allfalse(i);
		}
		void reset(const int r, const int c) {
			this->allfalse(r, c);
		}

		void flip() {
			this->v.flip();
		}

		void resize(const int i, const int j) {
			vcp::index_t nn = static_cast<vcp::index_t>(i) * j;
			int orow, ocolumn, on;
			orow = row;
			ocolumn = column;
			on = n;

			if (row > i || column > j) {
				vcp::throw_error<vcp::dimension_error>(
					"resize: new size is smaller than current size: (",
					row, ", ", column, ") > (", i, ", ", j, ")");
			}
			n = nn;
			row = i;
			column = j;

			if (row == 1 && column == 1) {
				type = 'S';
			}
			else if (column == 1) {
				type = 'C';
			}
			else if (row == 1) {
				type = 'R';
			}
			else if (row > 1 && column > 1) {
				type = 'M';
			}
			v.resize(nn);
			if (row > orow) {
				for (int jj = ocolumn - 1; jj >= 1; jj--) {
					for (int ii = orow - 1; ii >= 0; ii--) {
						v[ii + row*jj] = v[ii + orow*jj];
						v[ii + orow*jj] = false;
					}
				}
			}
		}
		void clear() {
			this->v.clear();
			this->row = 0;
			this->column = 0;
			this->n = 0;
			this->type = 'N';
		}

		
		friend bool all(const matrix< bool >& A) {
			return A.all();
		}
		friend bool any(const matrix< bool >& A) {
			return A.any();
		}
		friend bool none(const matrix< bool >& A) {
			return A.none();
		}

		friend matrix< bool > operator&&(const matrix< bool >& A, const matrix< bool >& B) {
			matrix< bool > C;
			C = A;
			C.mats_and(B);
			return C;
		}
		friend matrix< bool > operator&&(matrix< bool >&& A, const matrix< bool >& B) {
			A.mats_and(B);
			return std::move(A);
		}
		friend matrix< bool > operator&&(const matrix< bool >& A, matrix< bool >&& B) {
			B.mats_and(A);
			return std::move(B);
		}
		friend matrix< bool > operator&&(matrix< bool >&& A, matrix< bool >&& B) {
			A.mats_and(B);
			return std::move(A);
		}
		
		friend matrix< bool > operator||(const matrix< bool >& A, const matrix< bool >& B) {
			matrix< bool > C;
			C = A;
			C.mats_or(B);
			return C;
		}
		friend matrix< bool > operator||(matrix< bool >&& A, const matrix< bool >& B) {
			A.mats_or(B);
			return std::move(A);
		}
		friend matrix< bool > operator||(const matrix< bool >& A, matrix< bool >&& B) {
			B.mats_or(A);
			return std::move(B);
		}
		friend matrix< bool > operator||(matrix< bool >&& A, matrix< bool >&& B) {
			A.mats_or(B);
			return std::move(A);
		}

		friend matrix< bool > operator!(const matrix< bool >& A) {
			matrix< bool > C;
			C = A;
			C.flip();
			return C;
		}
		friend matrix< bool > operator!(matrix< bool >&& A) {
			A.flip();
			return std::move(A);
		}
		
		friend std::ostream& operator<<(std::ostream& os, const matrix< bool >& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, matrix< bool >&& A) {
			return A.display(os);
		}
	};
}
#endif // VCP_MATRIX_HPP
