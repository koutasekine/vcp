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

#ifndef VCP_MBOOL_HPP
#define VCP_MBOOL_HPP

#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>

#include <vcp/error.hpp>

namespace vcp {
	class mbool {
	public:
		int row;
		int column;
		int n;
		char type;      //'N':NULL  'S':Scala  'R' Row Vector 'C':Column Vector 'M':Matrix
		std::vector< bool > v;

		mbool() {
			this->row = 0;
			this->column = 0;
			this->n = 0;
			this->type = 'N';
		}

		~mbool() = default;
		mbool(const mbool&) = default;
		mbool(mbool&&) = default;
		mbool& operator=(const mbool& A) = default;
		mbool& operator=(mbool&& A) = default;

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

		int elementsize()const { return this->n; }
		int columnsize()const { return this->column; }
		int rowsize()const { return this->row; }
		char matstype()const {
			return this->type;
		}
		const std::vector< bool >& vecpointer()const {
			return this->v;
		}
		friend int length(mbool& A) {
			return A.length();
		}

		void alltrue(const int i) {
			row = i;
			column = i;
			n = i*i;
			if (i == 1) {
				type = 'S';
			}
			else {
				type = 'M';
			}
			v.resize(n);
			for (int j = 0; j < n; j++) {
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
			this->alltrue(r,c);
		}

		void allfalse(const int i) {
			row = i;
			column = i;
			n = i*i;
			if (i == 1) {
				type = 'S';
			}
			else {
				type = 'M';
			}
			v.resize(n);
			for (int j = 0; j < n; j++) {
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

		bool all() const {
			for (int i = 0; i < this->n; i++) {
				if (!this->v[i]) {
					return false;
				}
			}
			return true;
		}
		bool any() const {
			for (int i = 0; i < this->n; i++) {
				if (this->v[i]) {
					return true;
				}
			}
			return false;
		}
		bool none() const {
			for (int i = 0; i < this->n; i++) {
				if (this->v[i]) {
					return false;
				}
			}
			return true;
		}

		friend bool all(const mbool& A) {
			return A.all();
		}
		friend bool any(const mbool& A) {
			return A.any();
		}
		friend bool none(const mbool& A) {
			return A.none();
		}

		// A = and(A,B)
		void mats_and(const mbool& B) {
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>("&&: dimension mismatch");
			}
			if (this->type == 'S') {
				this->v[0] = this->v[0] && B.v[0];
				return;
			}
			else {
				for (int i = 0; i < n; i++) {
					this->v[i] = this->v[i] && B.v[i];
				}
				return;
			}
		}
		// A = or(A,B)
		void mats_or(const mbool& B) {
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>("||: dimension mismatch");
			}
			if (this->type == 'S') {
				this->v[0] = this->v[0] || B.v[0];
				return;
			}
			else {
				for (int i = 0; i < n; i++) {
					this->v[i] = this->v[i] || B.v[i];
				}
				return;
			}
		}

		friend mbool operator&&(const mbool& A, const mbool& B) {
			mbool C;
			C = A;
			C.mats_and(B);
			return C;
		}
		friend mbool operator&&(mbool&& A, const mbool& B) {
			A.mats_and(B);
			return std::move(A);
		}
		friend mbool operator&&(const mbool& A, mbool&& B) {
			B.mats_and(A);
			return std::move(B);
		}
		friend mbool operator&&(mbool&& A, mbool&& B) {
			A.mats_and(B);
			return std::move(A);
		}

		friend mbool operator||(const mbool& A, const mbool& B) {
			mbool C;
			C = A;
			C.mats_or(B);
			return C;
		}
		friend mbool operator||(mbool&& A, const mbool& B) {
			A.mats_or(B);
			return std::move(A);
		}
		friend mbool operator||(const mbool& A, mbool&& B) {
			B.mats_or(A);
			return std::move(B);
		}
		friend mbool operator||(mbool&& A, mbool&& B) {
			A.mats_or(B);
			return std::move(A);
		}

		friend mbool operator!(const mbool& A) {
			mbool C;
			C = A;
			C.flip();
			return C;
		}
		friend mbool operator!(mbool&& A) {
			A.flip();
			return std::move(A);
		}

		int length()const {
			using std::max;
			return max(column, row);
		}

		void resize(const int i, const int j) {
			int nn = i*j;
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
					for (int i = j; i <= row*(column - 1) + j; i = i + row) {
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
		friend std::ostream& operator<<(std::ostream& os, const mbool& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, mbool&& A) {
			return A.display(os);
		}
	};

} // namespace vcp

#endif // VCP_MBOOL_HPP
