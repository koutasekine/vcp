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

#ifndef VCP_MATS_HPP
#define VCP_MATS_HPP

#ifdef VCP_NOMP
#ifndef VCP_MATS_NOMP
#define VCP_MATS_NOMP
#endif
#endif

//If you define Micro of rand_debug, then the function rand display a 'warning'.
//#define rand_debug

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <initializer_list>

namespace vcp{
	class mbool {
	public:
		int row;
		int column;
		int n;
		char type;      //'N':NULL  'S':Scala  'R' Row Vector 'C':Column Vector 'M':Matrix
		std::vector< bool > v;

		mbool() {
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
		}

		virtual ~mbool() = default;
		mbool(const mbool&) = default;
		mbool(mbool&&) = default;
		mbool& operator=(const mbool& A) = default;
		mbool& operator=(mbool&& A) = default;

		std::vector< bool >::reference operator () (const int i) {
			return (*this).v[i];
		}
		std::vector< bool >::const_reference operator () (const int i) const {
			return (*this).v[i];
		}
		std::vector< bool >::reference operator () (const int i, const int j) {
			if ((*this).type == 'R') {
				return (*this).v[j];
			}
			else {
				return (*this).v[i + (*this).row*j];
			}
		}
		std::vector< bool >::const_reference operator () (const int i, const int j)const {
			if ((*this).type == 'R') {
				return (*this).v[j];
			}
			else {
				return (*this).v[i + (*this).row*j];
			}
		}

		int elementsize()const { return (*this).n; }
		int columnsize()const { return (*this).column; }
		int rowsize()const { return (*this).row; }
		char matstype()const {
			return (*this).type;
		}
		std::vector< bool > vecpointer()const {
			return (*this).v;
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
			v.shrink_to_fit();
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
			v.shrink_to_fit();
		}
		void set(const int i) {
			(*this).alltrue(i);
		}
		void set(const int r, const int c) {
			(*this).alltrue(r,c);
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
			v.shrink_to_fit();
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
			v.shrink_to_fit();
		}
		void reset(const int i) {
			(*this).allfalse(i);
		}
		void reset(const int r, const int c) {
			(*this).allfalse(r, c);
		}

		void flip() {
			(*this).v.flip();
		}

		bool all() const {
			for (int i = 0; i < (*this).n; i++) {
				if (!(*this).v[i]) {
					return false;
				}
			}
			return true;
		}
		bool any() const {
			for (int i = 0; i < (*this).n; i++) {
				if ((*this).v[i]) {
					return true;
				}
			}
			return false;
		}
		bool none() const {
			for (int i = 0; i < (*this).n; i++) {
				if ((*this).v[i]) {
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
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "&&:error " << std::endl;
				exit(1);
			}
			if ((*this).type == 'S') {
				(*this).v[0] = (*this).v[0] && B.v[0];
				return;
			}
			else {
				for (int i = 0; i < n; i++) {
					(*this).v[i] = (*this).v[i] && B.v[i];
				}
				return;
			}
		}
		// A = or(A,B)
		void mats_or(const mbool& B) {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "&&:error " << std::endl;
				exit(1);
			}
			if ((*this).type == 'S') {
				(*this).v[0] = (*this).v[0] || B.v[0];
				return;
			}
			else {
				for (int i = 0; i < n; i++) {
					(*this).v[i] = (*this).v[i] || B.v[i];
				}
				return;
			}
		}

		friend mbool operator&&(const mbool& A, const mbool& B) {
			mbool C;
			C = A;
			C.mats_and(B);
			return std::move(C);
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
			return std::move(C);
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
			return std::move(C);
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
				std::cout << "error : resize : row > i || column > j" << std::endl;
				exit(1);
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
			v.shrink_to_fit();
		}
		void clear() {
			(*this).v.clear();
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
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
		friend std::ostream& operator<<(std::ostream& os, mbool& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, mbool&& A) {
			return A.display(os);
		}
	};
	
	template <typename _T> class mats {
	public:
		int row;
		int column;
		int n;
		char type;      //'N':NULL  'S':Scala  'R' Row Vector 'C':Column Vector 'M':Matrix
		std::vector< _T > v;

		mats< _T >() {
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
		}
		virtual ~mats< _T >() = default;
		mats< _T >(const mats< _T >&) = default;
		mats< _T >(mats< _T >&&) = default;
		mats< _T >& operator=(const mats< _T >& A) = default;
		mats< _T >& operator=(mats< _T >&& A) = default;

		// A = A+B
		void addmm(const mats< _T >& B) {
			if (row != B.row && column != B.column) {
				std::cout << "+:error " << std::endl;
				exit(1);
			}
			if (type == 'S') {
				v[0] += B.v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif	
				for (int i = 0; i < n; i++) {
					v[i] += B.v[i];
				}
				return;
			}
		}
		// A = a+A
		void addsm(const _T a) {
			if (type == 'S') {
				v[0] = a + v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif	
				for (int i = 0; i < n; i++) {
					v[i] = a + v[i];
				}
				return;
			}
		}
		// A = A+a
		void addms(const _T a) {
			if (type == 'S') {
				v[0] = v[0] + a;
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = v[i] + a;
				}
				return;
			}
		}
		// A = +A
		void plusm() {
			return *this;
		}

		// A = A-B
		void subsmmA(const mats< _T >& B) {
			if (row != B.row && column != B.column) {
				std::cout << "-:error " << std::endl;
				exit(1);
			}
			if (type == 'S') {
				v[0] = v[0] - B.v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = v[i] - B.v[i];
				}
				return;
			}
		}
		// A = B-A
		void subsmmB(const mats< _T >& B) {
			if (row != B.row && column != B.column) {
				std::cout << "-:error " << std::endl;
				exit(1);
			}
			if (type == 'S') {
				v[0] = B.v[0] - v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = B.v[i] - v[i];
				}
				return;
			}
		}
		// A = a-A
		void subssm(const _T a) {
			if (type == 'S') {
				v[0] = a - v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = a - v[i];
				}
				return;
			}
		}
		// A = A-a
		void subsms(const _T a) {
			if (type == 'S') {
				v[0] = v[0] - a;
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = v[i] - a;
				}
				return;
			}
		}
		// A = -A
		void minusm() {
			if (type == 'S') {
				v[0] = -v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = -v[i];
				}
				return;
			}
		}

		// C = A*B
		virtual void mulmm(const mats< _T >& B, mats< _T >& c)const {
			if (type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c = B;
				for (int i = 0; i < B.n; i++) {
					c.v[i] *= v[0];
				}
				return;
			}
			else if ((type == 'C' || type == 'R' || type == 'M') && B.type == 'S') {
				c = *this;
				for (int i = 0; i < n; i++) {
					c.v[i] *= B.v[0];
				}
				return;
			}
			c.row = row;
			c.column = B.column;
			c.n = c.row * c.column;
			c.v.resize(c.n);
			if (type == 'S' && B.type == 'S') {
				c.type = 'S';
				c.v[0] = v[0] * B.v[0];
				return;
			}
			else if (type == 'R' && B.type == 'C' && column == B.row) {
				c.type = 'S';
				c.v[0] = v[0] * B.v[0];
				for (int i = 1; i < n; i++) {
					c.v[0] = v[i] * B.v[i] + c.v[0];
				}
				return;
			}
			else if (type == 'C' && B.type == 'R') {
				c.type = 'M';
				int k = 0;
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < B.n; j++) {
						c.v[k] = v[i] * B.v[j];
						k++;
					}
				}
				return;
			}
			else if (type == 'M' && B.type == 'C' && column == B.row) {
				c.type = 'C';
				for (int i = 0; i < row; i++) {
					c.v[i] = v[i] * B.v[0];
					for (int j = 1; j < column; j++) {
						c.v[i] += v[i + row*j] * B.v[j];
					}
				}
				return;
			}
			else if (type == 'R' && B.type == 'M' && column == B.row) {
				c.type = 'R';
				for (int i = 0; i < B.column; i++) {
					c.v[i] = v[0] * B.v[B.row*i];
					for (int j = 1; j < B.row; j++) {
						c.v[i] += v[j] * B.v[j + B.row*i];
					}
				}
				return;
			}
			else if (type == 'M' && B.type == 'M' && column == B.row) {
				c.type = 'M';

				for (int i = 0; i < c.n; i++) {
					c.v[i] = _T(0);
				}
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int k = 0; k < B.column; k++) {
					for (int j = 0; j < B.row; j++) {
						for (int i = 0; i < row; i++) {
							c.v[i + c.row*k] += v[i + row*j] * B.v[j + B.row*k];
						}
					}
				}
				return;
			}
			else {
				std::cout << "*:error " << (*this).type << " , " << B.type << std::endl;
				exit(1);
			}
		}
		// C = transpose(A)*A : multiplication left side transpose
		virtual void mulltmm(mats< _T >& c)const {
			c.zeros(column);
			if ((*this).type == 'S') {
				using std::pow;
				c.type = 'S';
				c.v[0] = pow((*this).v[0],2);
				return;
			}
			else if ((*this).type == 'C') {
				using std::pow;
				c.type = 'S';
				c.v[0] = _T(0);
				for (int i = 0; i < (*this).n; i++) {
					c.v[0] += pow((*this).v[i], 2);
				}
			}
			else if ((*this).type == 'R') {
				c.type = 'M';
				for (int i = 0; i < (*this).column; i++) {
					for (int j = i; j < (*this).column; j++) {
						c.v[i + (*this).column*j] = (*this).v[i] * (*this).v[j];
					}
				}
				for (int i = 0; i < (*this).column; i++) {
					for (int j = i+1; j < (*this).column; j++) {
						c.v[j + (*this).column*i] = c.v[i + (*this).column*j];
					}
				}
			}
			else if ((*this).type == 'M') {
				c.type = 'M';
				for (int i = 0; i < (*this).column; i++) {
					for (int j = i; j < (*this).column; j++) {
						for (int k = 0; k < (*this).row; k++) {
							c.v[i + (*this).column*j] += (*this).v[k + (*this).row*i] * (*this).v[k + (*this).row*j];
						}
					}
				}
				for (int i = 0; i < (*this).column; i++) {
					for (int j = i+1; j < (*this).column; j++) {
						c.v[j + (*this).column*i] = c.v[i + (*this).column*j];
					}
				}
			}
			else {
				std::cout << "*:error " << (*this).type << std::endl;
				exit(1);
			}
		}
		// A = a*A
		void mulsm(const _T a) {
			if (type == 'S') {
				v[0] = a * v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = a * v[i];
				}
				return;
			}
		}
		// A = A*a
		void mulms(const _T a) {
			if (type == 'S') {
				v[0] = v[0] * a;
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = v[i] * a;
				}
				return;
			}
		}

		// A = a/A
		void divsm(const _T a) {
			if (type == 'S') {
				v[0] = a / v[0];
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = a / v[i];
				}
				return;
			}
		}
		// A = A/a
		void divms(const _T a) {
			if (type == 'S') {
				v[0] = v[0] / a;
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = v[i] / a;
				}
				return;
			}
		}

		// C = A > B
		void gt(const mats< _T >& B, mbool& C) const {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "gt:error " << std::endl;
				exit(1);
			}
			C.allfalse((*this).row, (*this).column);
			for (int i = 0; i < (*this).n; i++) {
				C.v[i] = (*this).v[i] > B.v[i];
			}
		}
		// C = A >= B
		void ge(const mats< _T >& B, mbool& C) const {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "gt:error " << std::endl;
				exit(1);
			}
			C.allfalse((*this).row, (*this).column);
			for (int i = 0; i < (*this).n; i++) {
				C.v[i] = (*this).v[i] >= B.v[i];
			}
		}
		// C = A < B
		void lt(const mats< _T >& B, mbool& C) const {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "gt:error " << std::endl;
				exit(1);
			}
			C.allfalse((*this).row, (*this).column);
			for (int i = 0; i < (*this).n; i++) {
				C.v[i] = (*this).v[i] < B.v[i];
			}
		}
		// C = A <= B
		void le(const mats< _T >& B, mbool& C) const {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "gt:error " << std::endl;
				exit(1);
			}
			C.allfalse((*this).row, (*this).column);
			for (int i = 0; i < (*this).n; i++) {
				C.v[i] = (*this).v[i] <= B.v[i];
			}
		}
		// C = A == B
		void eq(const mats< _T >& B, mbool& C) const {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "gt:error " << std::endl;
				exit(1);
			}
			C.allfalse((*this).row, (*this).column);
			for (int i = 0; i < (*this).n; i++) {
				C.v[i] = (*this).v[i] == B.v[i];
			}
		}
		// C = A != B
		void neq(const mats< _T >& B, mbool& C) const {
			if ((*this).row != B.row && (*this).column != B.column) {
				std::cout << "gt:error " << std::endl;
				exit(1);
			}
			C.allfalse((*this).row, (*this).column);
			for (int i = 0; i < (*this).n; i++) {
				C.v[i] = (*this).v[i] != B.v[i];
			}
		}

		// A = pow(A,B)
		void powmmA(const mats< _T >& B) {
			if (row != B.row && column != B.column) {
				std::cout << "+:error " << std::endl;
				exit(1);
			}
			if (type == 'S') {
				using std::pow;
				v[0] = pow(v[0], B.v[0]);
				return;
			}
			else {
				using std::pow;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = pow(v[i], B.v[i]);
				}
				return;
			}
		}
		// A = pow(B,A)
		void powmmB(const mats< _T >& B) {
			if (row != B.row && column != B.column) {
				std::cout << "+:error " << std::endl;
				exit(1);
			}
			if (type == 'S') {
				using std::pow;
				v[0] = pow(B.v[0], v[0]);
				return;
			}
			else {
				using std::pow;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					v[i] = pow(B.v[i], v[i]);
				}
				return;
			}
		}
		// A = pow(A,a)
		void powms(const _T a) {
			if (type == 'S') {
				using std::pow;
				v[0] = pow(v[0], a);
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					using std::pow;
					v[i] = pow(v[i], a);
				}
				return;
			}
		}
		// A = pow(a,A)
		void powsm(const _T a) {
			if (type == 'S') {
				using std::pow;
				v[0] = pow(a, v[0]);
				return;
			}
			else {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < n; i++) {
					using std::pow;
					v[i] = pow(a, v[i]);
				}
				return;
			}
		}

		// A = abs(A)
		void abs() {
			using std::abs;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < n; i++) {
				v[i] = abs(v[i]);
			}
		}
		// A = sqrt(A)
		void sqrt() {
			using std::sqrt;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < n; i++) {
				v[i] = sqrt(v[i]);
			}
		}
		// A = sin(A)
		void sin() {
			using std::sin;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < n; i++) {
				v[i] = sin(v[i]);
			}
		}
		// A = cos(A)
		void cos() {
			using std::cos;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < n; i++) {
				v[i] = cos(v[i]);
			}
		}
		// A = exp(A)
		void exp() {
			using std::exp;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < n; i++) {
				v[i] = exp(v[i]);
			}
		}
		// A = log(A)
		void log() {
			using std::log;
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < n; i++) {
				v[i] = log(v[i]);
			}
		}

		void diag(mats< _T >& B)const {
			if (type == 'S') {
				B = *this;
				return;
			}
			else if (type == 'R' || type == 'C') {
				B.zeros(n);
				for (int i = 0; i < n; i++) {
					B.v[i + B.row*i] = v[i];
				}
				return;
			}
			else if (type == 'M') {
				using std::min;
				int nn = min(column, row);
				B.zeros(nn, 1);
				for (int i = 0; i < nn; i++) {
					B.v[i] = v[i + row*i];
				}
				return;
			}
			else {
				std::cout << "error : diag" << std::endl;
				exit(1);
			}
		}
		void transpose(mats< _T >& B)const {
			if (type == 'S') {
				B = *this;
				return;
			}
			else if (type == 'R') {
				B = *this;
				B.row = (*this).column;
				B.column = (*this).row;
				B.type = 'C';
			}
			else if (type == 'C') {
				B = *this;
				B.row = (*this).column;
				B.column = (*this).row;
				B.type = 'R';
			}
			else if (type == 'M') {
				B.zeros(column, row);
				for (int i = 0; i < row; i++) {
					for (int j = 0; j < column; j++) {
						B.v[j + B.row*i] = v[i + row*j];
					}
				}
			}
			else {
				std::cout << "error : transpose" << std::endl;
				exit(1);
			}
		}

		// matlab C = [A,B]
		void horzcat(const mats< _T >& B, mats< _T >& C)const {
			int An, Am, Bn, Bm;
			An = row;
			Am = column;
			Bn = B.row;
			Bm = B.column;

			if (An != Bn) {
				std::cout << "Error : horzcat : An != Bn" << std::endl;
				exit(1);
			}
			C.zeros(An, Am);
			for (int i = 0; i < n; i++) {
				C.v[i] = v[i];
			}
			C.resize(An, Am + Bm);
			for (int i = 0; i < An; i++) {
				for (int j = 0; j < Bm; j++) {
					C.v[i + C.row*(j + Am)] = B.v[i + B.row*j];
				}
			}
			return;
		}
		// matlab [A;B]
		void vercat(const mats< _T >& B, mats< _T >& C)const {
			int An, Am, Bn, Bm;
			An = row;
			Am = column;
			Bn = B.row;
			Bm = B.column;

			if (Am != Bm) {
				std::cout << "Error : vercat : Am != Bm" << std::endl;
				exit(1);
			}
			C.zeros(An, Am);
			for (int i = 0; i < n; i++) {
				C.v[i] = v[i];
			}
			C.resize(An + Bn, Am);
			for (int i = 0; i < Bn; i++) {
				for (int j = 0; j < Am; j++) {
					C.v[(i + An) + C.row*j] = B.v[i + B.row*j];
				}
			}
			return;
		}

		void submat(mats< _T >& B, const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
			std::vector<int> l1 = list1;
			std::vector<int> l2 = list2;

			if (list1.size() > 3 || list2.size() > 3 || list1.size() < 0 || list2.size() < 0) {
				std::cout << "ERROR : submat : size error" << std::endl;
				exit(1);
			}
			if (list1.size() == 0 ) {
				if (list2.size() == 0) {
					B = (*this);
					return;
				}
				else if (list2.size() == 1) {
					if (l2[0] < 0 || l2[0] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << std::endl;
						exit(1);
					}
					B.zeros((*this).row, 1);
					for (int i = 0; i < (*this).row; i++) {
						B.v[i] = (*this).v[i + (*this).row*l2[0]];
					}
					return;
				}
				else if (list2.size() == 2) {
					if (l2[0] > l2[1] || l2[0] < 0 || l2[1] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[1] << std::endl;
						exit(1);
					}
					int l2i = l2[1] - l2[0] + 1;
					B.zeros((*this).row, l2i);
					int k = 0;
					for (int i = 0; i < (*this).row; i++) {
						for (int j = l2[0]; j <= l2[1]; j++) {
							B.v[i + B.row*k] = (*this).v[i + (*this).row*j];
							k++;
						}
						k = 0;
					}
					return;
				}
				else if (list2.size() == 3) {
					if (l2[0] > l2[2] || l2[0] < 0 || l2[1] < 1 || l2[2] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[2] << " || " << l2[1] << " < 1" << std::endl;
						exit(1);
					}
					int k = 0;
					for (int i = l2[0]; i <= l2[2]; i += l2[1]) {
						k++;
					}

					B.zeros((*this).row, k);
					k = 0;
					for (int i = 0; i < (*this).row; i++) {
						for (int j = l2[0]; j <= l2[2]; j += l2[1]) {
							B.v[i + B.row*k] = (*this).v[i + (*this).row*j];
							k++;
						}
						k = 0;
					}
					return;
				}
			}
			else if (list1.size() == 1) {
				if (l1[0] < 0 || l1[0] >= (*this).row) {
					std::cout << "ERROR : submat :" << l1[0] << std::endl;
					exit(1);
				}
				if (list2.size() == 0) {
					B.zeros(1, (*this).column);
					for (int i = 0; i < (*this).column; i++) {
						B.v[i] = (*this).v[l1[0] + (*this).row*i];
					}
					return;
				}
				else if (list2.size() == 1 || l2[0] >= (*this).column) {
					if (l2[0] < 0) {
						std::cout << "ERROR : submat :" << l2[0] << std::endl;
						exit(1);
					}
					B.zeros(1, 1);
					B.v[0] = (*this).v[l1[0] + (*this).row*l2[0]];
					return;
				}
				else if (list2.size() == 2) {
					if (l2[0] > l2[1] || l2[0] < 0 || l2[1] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[1] << std::endl;
						exit(1);
					}
					int l2i = l2[1] - l2[0] + 1;
					B.zeros(1, l2i);
					int k = 0;
					for (int j = l2[0]; j <= l2[1]; j++) {
						B.v[k] = (*this).v[l1[0] + (*this).row*j];
						k++;
					}
					return;
				}
				else if (list2.size() == 3) {
					if (l2[0] > l2[2] || l2[0] < 0 || l2[1] < 1 || l2[2] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[2] << " || " << l2[1] << " < 1" << std::endl;
						exit(1);
					}
					int k = 0;
					for (int i = l2[0]; i <= l2[2]; i += l2[1]) {
						k++;
					}

					B.zeros(1, k);
					k = 0;
					for (int j = l2[0]; j <= l2[2]; j += l2[1]) {
						B.v[k] = (*this).v[l1[0] + (*this).row*j];
						k++;
					}
					return;
				}
			}
			else if (list1.size() == 2) {
				if (l1[0] > l1[1] || l1[0] < 0 || l1[1] >= (*this).row) {
					std::cout << "ERROR : submat :" << l1[0] << " >= " << l1[1] << std::endl;
					exit(1);
				}
				int l1i = l1[1] - l1[0] + 1;
				if (list2.size() == 0) {
					B.zeros(l1i, (*this).column);
					int k = 0;
					for (int i = l1[0]; i <= l1[1]; i++) {
						for (int j = 0; j < (*this).column; j++) {
							B.v[k + B.row*j] = (*this).v[i + (*this).row*j];	
						}
						k++;
					}
					return;
				}
				else if (list2.size() == 1) {
					if (l2[0] < 0 || l2[0] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << std::endl;
						exit(1);
					}
					B.zeros(l1i, 1);
					int k = 0;
					for (int i = l1[0]; i <= l1[1]; i++) {
						B.v[k] = (*this).v[i + (*this).row*l2[0]];
						k++;
					}
					return;
				}
				else if (list2.size() == 2) {
					if (l2[0] > l2[1] || l2[0] < 0 || l2[1] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[1] << std::endl;
						exit(1);
					}
					int l2i = l2[1] - l2[0] + 1;
					B.zeros(l1i, l2i);
					int ki = 0;
					int kj = 0;
					for (int i = l1[0]; i <= l1[1]; i++) {
						for (int j = l2[0]; j <= l2[1]; j++) {
							B.v[ki + B.row*kj] = (*this).v[i + (*this).row*j];
							kj++;
						}
						kj = 0;
						ki++;
					}
					return;
				}
				else if (list2.size() == 3) {
					if (l2[0] > l2[2] || l2[0] < 0 || l2[1] < 1 || l2[2] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[2] << " || " << l2[1] << " < 1" << std::endl;
						exit(1);
					}
					int ki = 0;
					for (int i = l2[0]; i <= l2[2]; i += l2[1]) {
						ki++;
					}

					B.zeros(l1i, ki);
					ki = 0;
					int kj = 0;
					for (int i = l1[0]; i <= l1[1]; i++) {
						for (int j = l2[0]; j <= l2[2]; j += l2[1]) {
							B.v[ki + B.row*kj] = (*this).v[i + (*this).row*j];
							kj++;
						}
						kj = 0;
						ki++;
					}
					return;
				}
			}
			else if (list1.size() == 3) {
				if (l1[0] > l1[2] || l1[0] < 0 || l1[1] < 1 || l1[2] >= (*this).row) {
					std::cout << "ERROR : submat :" << l1[0] << " >= " << l1[2] << " || " << l1[1] << " < 1" << std::endl;
					exit(1);
				}
				int l1i = 0;
				for (int i = l1[0]; i <= l1[2]; i += l1[1]) {
					l1i++;
				}
				if (list2.size() == 0) {
					B.zeros(l1i, (*this).column);
					int k = 0;
					for (int i = l1[0]; i <= l1[2]; i += l1[1]) {
						for (int j = 0; j < (*this).column; j++) {
							B.v[k + B.row*j] = (*this).v[i + (*this).row*j];
						}
						k++;
					}
					return;
				}
				else if (list2.size() == 1) {
					if (l2[0] < 0 || l2[0] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << std::endl;
						exit(1);
					}
					B.zeros(l1i, 1);
					int k = 0;
					for (int i = l1[0]; i <= l1[2]; i += l1[1]) {
						B.v[k] = (*this).v[i + (*this).row*l2[0]];
						k++;
					}
					return;
				}
				else if (list2.size() == 2) {
					if (l2[0] > l2[1] || l2[0] < 0 || l2[1] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[1] << std::endl;
						exit(1);
					}
					int l2i = l2[1] - l2[0] + 1;
					B.zeros(l1i, l2i);
					int ki = 0;
					int kj = 0;
					for (int i = l1[0]; i <= l1[2]; i += l1[1]) {
						for (int j = l2[0]; j <= l2[1]; j++) {
							B.v[ki + B.row*kj] = (*this).v[i + (*this).row*j];
							kj++;
						}
						kj = 0;
						ki++;
					}
					return;
				}
				else if (list2.size() == 3) {
					if (l2[0] > l2[2] || l2[0] < 0 || l2[1] < 1 || l2[2] >= (*this).column) {
						std::cout << "ERROR : submat :" << l2[0] << " >= " << l2[2] << " || " << l2[1] << " < 1" << std::endl;
						exit(1);
					}
					int ki = 0;
					for (int i = l2[0]; i <= l2[2]; i += l2[1]) {
						ki++;
					}

					B.zeros(l1i, ki);
					ki = 0;
					int kj = 0;
					for (int i = l1[0]; i <= l1[2]; i += l1[1]) {
						for (int j = l2[0]; j <= l2[2]; j += l2[1]) {
							B.v[ki + B.row*kj] = (*this).v[i + (*this).row*j];
							kj++;
						}
						kj = 0;
						ki++;
					}
					return;
				}
			}
		}

		void sum(mats< _T >& B)const {
			if (type == 'S') {
				B.zeros(1);
				B.v[0] = v[0];
				return;
			}
			else if (type == 'C' || type == 'R') {
				B.zeros(1);
				B.v[0] = v[0];
				for (int i = 1; i < n; i++) {
					B.v[0] += v[i];
				}
				return;
			}
			else if (type == 'M') {
				B.zeros(1, column);
				for (int i = 0; i < column; i++) {
					int ii = row*i;
					B.v[i] = v[ii];
					for (int j = 1; j < row; j++) {
						B.v[i] = B.v[i] + v[j + ii];
					}
				}
				return;
			}
			else {
				std::cout << "error : sum" << std::endl;
				exit(1);
			}
		}
		void max(mats< _T >& B)const {
			if (type == 'S') {
				B.zeros(1);
				B.v[0] = v[0];
				return;
			}
			else if (type == 'C' || type == 'R') {
				B.zeros(1);
				B.v[0] = v[0];
				using std::max;
				for (int i = 1; i < n; i++) {
					B.v[0] = max(v[i], B.v[0]);
				}
				return;
			}
			else if (type == 'M') {
				B.zeros(1, column);
				using std::max;
				for (int i = 0; i < column; i++) {
					int ii = row*i;
					B.v[i] = v[ii];
					for (int j = 1; j < row; j++) {
						B.v[i] = max(B.v[i], v[j + ii]);
					}
				}
				return;
			}
			else {
				std::cout << "error : max" << std::endl;
				exit(1);
			}
		}
		void min(mats< _T >& B)const {
			if (type == 'S') {
				B.zeros(1);
				B.v[0] = v[0];
				return;
			}
			else if (type == 'C' || type == 'R') {
				B.zeros(1);
				B.v[0] = v[0];
				using std::min;
				for (int i = 1; i < n; i++) {
					B.v[0] = min(v[i], B.v[0]);
				}
				return;
			}
			else if (type == 'M') {
				B.zeros(1, column);
				using std::min;
				for (int i = 0; i < column; i++) {
					int ii = row*i;
					B.v[i] = v[ii];
					for (int j = 1; j < row; j++) {
						B.v[i] = min(B.v[i], v[j + ii]);
					}
				}
				return;
			}
			else {
				std::cout << "error : max" << std::endl;
				exit(1);
			}
		}
		void normone(mats< _T >& B)const {
			if (type == 'S') {
				B.zeros(1);
				using std::abs;
				B.v[0] = abs(v[0]);
				return;
			}
			else if (type == 'C' || type == 'R') {
				B.zeros(1);
				using std::max;
				using std::abs;
				B.v[0] = abs(v[0]);
				for (int i = 1; i < n; i++) {
					B.v[0] += abs(v[i]);
				}
				return;
			}
			else if (type == 'M') {
				B.zeros(1, column);
				using std::max;
				using std::abs;
				for (int i = 0; i < column; i++) {
					int ii = row*i;
					B.v[i] = abs(v[ii]);
					for (int j = 1; j < row; j++) {
						B.v[i] += abs(v[j + ii]);
					}
				}
				mats< _T > C;
				C.ones(1, 1);
				C.v[0] = B.v[0];
				for (int i = 1; i < B.n; i++) {
					C.v[0] = max(B.v[i], C.v[0]);
				}
				B = C;
				return;
			}
			else {
				std::cout << "error : norminf" << std::endl;
				exit(1);
			}
		}
		void norminf(mats< _T >& B)const {
			if (type == 'S') {
				B.zeros(1);
				using std::abs;
				B.v[0] = abs(v[0]);
				return;
			}
			else if (type == 'C' || type == 'R') {
				B.zeros(1);
				using std::max;
				using std::abs;
				B.v[0] = abs(v[0]);
				for (int i = 1; i < n; i++) {
					B.v[0] = max(abs(v[i]), B.v[0]);
				}
				return;
			}
			else if (type == 'M') {
				B.zeros(row, 1);
				using std::max;
				using std::abs;
				for (int i = 0; i < row; i++) {
					B.v[i] = abs(v[i]);
					for (int j = 1; j < column; j++) {
						B.v[i] += abs(v[i + row*j]);
					}
				}
				mats< _T > C;
				C.ones(1, 1);
				C.v[0] = B.v[0];
				for (int i = 1; i < B.n; i++) {
					C.v[0] = max(B.v[i], C.v[0]);
				}
				B = C;
				return;
			}
			else {
				std::cout << "error : norminf" << std::endl;
				exit(1);
			}
		}
		int length()const {
			using std::max;
			return max(column, row);
		}

		void tril() {
			_T T0 = _T(0);
			if (row <= column) {
				for (int i = 0; i < row; i++) {
					for (int j = i + 1; j < column; j++) {
						v[i + row*j] = T0;
					}
				}
			}
			else {
				for (int i = 0; i < column; i++) {
					for (int j = i + 1; j < column; j++) {
						v[i + row*j] = T0;
					}
				}


			}
		}
		void triu() {
			_T T0 = _T(0);
			if (row <= column) {
				for (int i = 1; i < row; i++) {
					for (int j = 0; j < i; j++) {
						v[i + row*j] = T0;
					}
				}
			}
			else {
				for (int i = 1; i < column; i++) {
					for (int j = 0; j < i; j++) {
						v[i + row*j] = T0;
					}
				}
				for (int i = column; i < row; i++) {
					for (int j = 0; j < column; j++) {
						v[i + row*j] = T0;
					}
				}

			}
		}
		virtual void LUtoLandU(mats< _T >& L, mats< int >& ipiv) {
			if (column != row) {
				std::cout << "error : LUtoLandU row != column" << std::endl;
				exit(1);
			}
			L = *this;
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < row; j++) {
					L.v[i + row*j] = v[ipiv.v[i] + row*j];
				}
			}
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < row; j++) {
					v[i + row*j] = L.v[i + row*j];
				}
			}
			L.tril();
			_T T1 = _T(1);
			for (int i = 0; i < row; i++) {
				L.v[i + row*i] = T1;
			}
			(*this).triu();
		}
		void ludecomposition(mats< int >& ipiv) {
			if (column != row) {
				std::cout << "error : ludeconposition row != column" << std::endl;
				exit(1);
			}
			int nn = row;
			ipiv.zeros(nn, 1);
			for (int i = 0; i < nn; i++) {
				ipiv.v[i] = i;
			}

			_T amax, dd;
			int j, imax;
			using std::abs;
			for (int k = 0; k < nn - 1; k++) {
				amax = abs(v[ipiv.v[k] + row*k]);
				imax = k;
				for (int i = k + 1; i < nn; i++) {
					if (amax < abs(v[ipiv.v[i] + row*k])) {
						amax = abs(v[ipiv.v[i] + row*k]);
						imax = i;
					}
				}
				if (imax != k) {
					j = ipiv.v[k];
					ipiv.v[k] = ipiv.v[imax];
					ipiv.v[imax] = j;
				}
				for (int i = k + 1; i < nn; i++) {
					dd = v[ipiv.v[i] + row*k] / v[ipiv.v[k] + row*k];
					for (j = k + 1; j < nn; j++) {
						v[ipiv.v[i] + row*j] = v[ipiv.v[i] + row*j] - dd * v[ipiv.v[k] + row*j];
					}
					v[ipiv.v[i] + row*k] = dd;
				}
			}
		}
		void luonedsolve(mats< _T >& b, mats< _T >& x, mats< int >& ipiv) const {
			if (row != column || b.column != 1 || b.row != row) {
				std::cout << "error : luonedsolve row != column || b.column != 1 || b.row != row" << std::endl;
				exit(1);
			}
			x.zeros(row, 1);
			for (int j = 0; j < row - 1; j++) {
				for (int i = j + 1; i < row; i++) {
					b.v[ipiv.v[i]] = b.v[ipiv.v[i]] - v[ipiv.v[i] + row*j] * b.v[ipiv.v[j]];
				}
			}
			_T dd;
			x.v[row - 1] = b.v[ipiv.v[row - 1]] / v[ipiv.v[row - 1] + row*(row - 1)];
			for (int i = row - 2; i >= 0; i--) {
				dd = b.v[ipiv.v[i]];
				for (int j = i + 1; j < row; j++) {
					dd = dd - v[ipiv.v[i] + row*j] * x.v[j];
				}
				x.v[i] = dd / v[ipiv.v[i] + row*i];
			}
		}
		void linearsolve(mats< _T >& b, mats< _T >& x) {
			if (row != column || b.row != row) {
				std::cout << "error : linearsolve row != column || b.row != row" << std::endl;
				exit(1);
			}
			mats< int > ipiv;
			(*this).ludecomposition(ipiv);
			mats< _T > bb, xx;
			x.zeros(row, b.column);
			bb.zeros(row, 1);
			xx.ones(row, 1);
			for (int j = 0; j < b.column; j++) {
				for (int i = 0; i < b.row; i++) {
					bb.v[i] = b.v[i + row*j];
				}
				(*this).luonedsolve(bb, xx, ipiv);
				for (int i = 0; i < b.row; i++) {
					x.v[i + row*j] = xx.v[i];
				}
			}
		}
		virtual void inv() {
			if (row != column) {
				std::cout << "error : inv row != column" << std::endl;
				exit(1);
			}
			mats< int > ipiv;
			(*this).ludecomposition(ipiv);
			mats< _T > b;
			b.eye(row);
			mats< _T > bb, xx;
			bb.zeros(row, 1);
			xx.ones(row, 1);
			for (int j = 0; j < b.column; j++) {
				for (int i = 0; i < b.row; i++) {
					bb.v[i] = b.v[i + row*j];
				}
				(*this).luonedsolve(bb, xx, ipiv);
				for (int i = 0; i < b.row; i++) {
					b.v[i + row*j] = xx.v[i];
				}
			}
			for (int j = 0; j < b.column; j++) {
				for (int i = 0; i < b.row; i++) {
					v[i + row*j] = b.v[i + row*j];
				}
			}
		}
		// Cholesky decomposition (Cholesky factor : upper triangular)
		virtual void Cholesky() {
			if (!(*this).is_symmetric()) {
				std::cout << "error : Cholesky : no symmetric matrix" << std::endl;
				exit(1);
			}
			using std::sqrt;

			for (int i = 0; i < row; i++) {
				v[i + row * i] = v[i + row * i];
				for (int k = 0; k < i; k++)
					v[i + row * i] -= v[k + row * i] * v[k + row * i];
				if (v[i + row * i] < _T(0)) {
					std::cout << "error : negative value in sqrt" << std::endl;
				}
				v[i + row * i] = sqrt(v[i + row * i]);
				for (int j = i + 1; j < row; j++) {
					v[i + row * j] = v[i + row * j];
					for (int k = 0; k < i; k++)
						v[i + row * j] -= v[k + row * i] * v[k + row * j];
					v[i + row * j] /= v[i + row * i];
				}
			}

			for (int i = 0; i < row; i++) {
				for (int j = 0; j < i; j++) {
					v[i + row * j] = _T(0);
				}
			}

		}
		void TriLudecomposition(mats< int >& ipiv) {
			if (column != row) {
				std::cout << "error : TriLudecomposition row != column" << std::endl;
				exit(1);
			}
			int nn = row;
			ipiv.zeros(nn, 1);
			for (int i = 0; i < nn; i++) {
				ipiv.v[i] = i;
			}

			_T amax, dd;
			int j, imax;
			using std::abs;
			using std::min;
			for (int k = 0; k < nn - 1; k++) {
				amax = abs(v[ipiv.v[k] + row*k]);
				imax = k;
				for (int i = k + 1; i < min(k + 2, nn); i++) {
					if (amax < abs(v[ipiv.v[i] + row*k])) {
						amax = abs(v[ipiv.v[i] + row*k]);
						imax = i;
					}
				}
				if (imax != k) {
					j = ipiv.v[k];
					ipiv.v[k] = ipiv.v[imax];
					ipiv.v[imax] = j;
				}
				for (int i = k + 1; i < min(k + 3, nn); i++) {
					dd = v[ipiv.v[i] + row*k] / v[ipiv.v[k] + row*k];
					for (j = k + 1; j < min(k + 3, nn); j++) {
						v[ipiv.v[i] + row*j] = v[ipiv.v[i] + row*j] - dd * v[ipiv.v[k] + row*j];
					}
					v[ipiv.v[i] + row*k] = dd;
				}
			}
		}

		// Householder tansformation for a Symmetric Matrix
		void SymTridiagonalization() {
			if (row != column) {
				std::cout << "error : SymTridiagonalization row != column" << std::endl;
				exit(1);
			}
			_T s = _T(0);
			_T c, tmp;
			mats< _T > w, q;
			q.zeros(row, 1);
			using std::sqrt;
			for (int i = 0; i < row - 1; i++) {
				// make s and s2
				_T s2 = _T(0);
				for (int j = i + 1; j < row; j++) {
					s2 += v[j + row*i] * v[j + row*i];
				}
				s = sqrt(s2);
				if (v[(i + 1) + row*i] < 0) {
					s = -s;
				}
				// make c
				c = _T(1) / (s2 + v[(i + 1) + row*i] * s);
				// make w
				w.zeros(row, 1);
				w.v[i + 1] = v[(i + 1) + row*i] + s;
				for (int k = i + 2; k < row; k++) {
					w.v[k] = v[k + row*i];
				}
				// make q = c*A*w
				for (int ii = 0; ii < row; ii++) {
					q.v[ii] = v[ii] * w.v[0];
					for (int jj = 1; jj < column; jj++) {
						q.v[ii] += v[ii + row*jj] * w.v[jj];
					}
				}
				for (int ii = 0; ii < row; ii++) {
					q.v[ii] *= c;
				}
				// make tmp = q^T * w
				tmp = q.v[0] * w.v[0];
				for (int ii = 1; ii < row; ii++) {
					tmp += q.v[ii] * w.v[ii];
				}
				tmp *= c / _T(2);
				for (int ii = 0; ii < row; ii++) {
					q.v[ii] = q.v[ii] - tmp*w.v[ii];
				}
				for (int jj = 0; jj < row; jj++) {
					for (int ii = 0; ii < row; ii++) {
						v[ii + row*jj] = v[ii + row*jj] - (q.v[ii] * w.v[jj] + w.v[ii] * q.v[jj]);
					}
				}
			}
			s = _T(0);
			for (int ii = 2; ii < row; ii++) {
				for (int jj = 0; jj < ii - 1; jj++) {
					v[ii + row*jj] = s;
					v[jj + row*ii] = s;
				}
			}
		}
		// Householder tansformation for a Symmetric Matrix with Orthogonal Matrix
		//TriA = P*A*transpose(P)
		void SymTridiagonalization(mats< _T >& P) {
			if (row != column) {
				std::cout << "error : SymTridiagonalization row != column" << std::endl;
				exit(1);
			}
			_T s = _T(0);
			_T c, tmp;
			mats< _T > w, q, PP, Ptmp;
			q.zeros(row, 1);
			P.eye(row);
			using std::sqrt;
			for (int i = 0; i < row - 1; i++) {
				// make s and s2
				_T s2 = _T(0);
				for (int j = i + 1; j < row; j++) {
					s2 += v[j + row*i] * v[j + row*i];
				}
				s = sqrt(s2);
				if (v[(i + 1) + row*i] < 0) {
					s = -s;
				}
				// make c
				c = _T(1) / (s2 + v[(i + 1) + row*i] * s);
				// make w
				w.zeros(row, 1);
				w.v[i + 1] = v[(i + 1) + row*i] + s;
				for (int k = i + 2; k < row; k++) {
					w.v[k] = v[k + row*i];
				}
				// make P = I-cw*w^T
				PP.eye(row - i);
				for (int jj = 0; jj < row - i; jj++) {
					for (int ii = 0; ii < row - i; ii++) {
						PP.v[ii + (row - i)*jj] -= c*w.v[ii + i] * w.v[jj + i];
					}
				}
				Ptmp.zeros(row - i, row);
				for (int kk = 0; kk < row; kk++) {
					for (int jj = i; jj < row; jj++) {
						for (int ii = i; ii < row; ii++) {
							Ptmp.v[ii - i + (row - i)*kk] += PP.v[ii - i + (row - i)*(jj - i)] * P.v[jj + row*kk];
						}
					}
				}
				for (int kk = 0; kk < row; kk++) {
					for (int ii = i; ii < row; ii++) {
						P.v[ii + row*kk] = Ptmp.v[ii - i + (row - i)*kk];
					}
				}

				// make q = c*A*w
				for (int ii = 0; ii < row; ii++) {
					q.v[ii] = v[ii] * w.v[0];
					for (int jj = 1; jj < column; jj++) {
						q.v[ii] += v[ii + row*jj] * w.v[jj];
					}
				}
				for (int ii = 0; ii < row; ii++) {
					q.v[ii] *= c;
				}
				// make tmp = q^T * w
				tmp = q.v[0] * w.v[0];
				for (int ii = 1; ii < row; ii++) {
					tmp += q.v[ii] * w.v[ii];
				}
				tmp *= c / _T(2);
				for (int ii = 0; ii < row; ii++) {
					q.v[ii] = q.v[ii] - tmp*w.v[ii];
				}
				for (int jj = 0; jj < row; jj++) {
					for (int ii = 0; ii < row; ii++) {
						v[ii + row*jj] = v[ii + row*jj] - (q.v[ii] * w.v[jj] + w.v[ii] * q.v[jj]);
					}
				}
			}
			s = _T(0);
			for (int ii = 2; ii < row; ii++) {
				for (int jj = 0; jj < ii - 1; jj++) {
					v[ii + row*jj] = s;
					v[jj + row*ii] = s;
				}
			}
		}
		// Householder tansformation for a General Matrix
		void HessenbergTrans() {
			if (row != column) {
				std::cout << "error : HessenbergTrans row != column" << std::endl;
				exit(1);
			}
			_T s = _T(0);
			_T c;
			mats< _T > w, P;
			using std::sqrt;
			for (int i = 0; i < row - 1; i++) {
				// make s and s2
				_T s2 = _T(0);
				for (int j = i + 1; j < row; j++) {
					s2 += v[j + row*i] * v[j + row*i];
				}
				s = sqrt(s2);
				if (v[(i + 1) + row*i] < 0) {
					s = -s;
				}
				// make c
				c = _T(1) / (s2 + v[(i + 1) + row*i] * s);
				// make w
				w.zeros(row, 1);
				w.v[i + 1] = v[(i + 1) + row*i] + s;
				for (int k = i + 2; k < row; k++) {
					w.v[k] = v[k + row*i];
				}
				// make P = I-cw*w^T
				P.eye(row);
				for (int ii = 0; ii < row; ii++) {
					for (int jj = 0; jj < row; jj++) {
						P.v[ii + row*jj] -= c*w.v[ii] * w.v[jj];
					}
				}
				//Please correct
				P.mulmm((*this), w);
				w.mulmm(P, (*this));
//				(*this) = P*(*this)*P;
			}
			s = _T(0);
			for (int ii = 2; ii < row; ii++) {
				for (int jj = 0; jj < ii - 1; jj++) {
					v[ii + row*jj] = s;
				}
			}
		}

		// QR deconposition for a General Matrix using Householder tansformation
		void Householder_qrdecomposition(mats< _T >& Q) {
			if (row != column) {
				std::cout << "error : Householder_qrdecomposition row != column" << std::endl;
				exit(1);
			}
			using std::pow;
			using std::sqrt;
			mats< _T > V, H, QRs;
			mats< _T > u;
			Q.eye(row);
			_T xynorm, xynorm2;
			_T zero = _T(0);

			for (int j = 0; j < row - 1; j++) {
				xynorm = zero;
				xynorm2 = zero;
				for (int i = j + 1; i < row; i++) {
					xynorm2 += pow(v[i + row*j], 2);
				}
				xynorm = xynorm2 + pow(v[j + row*j], 2);
				u.zeros(row - j, 1);

				if (v[j + row*j] >= 0) {
					u.v[0] = v[j + row*j] + sqrt(xynorm);
				}
				else if (v[j + row*j] < 0) {
					u.v[0] = v[j + row*j] - sqrt(xynorm);
				}
				xynorm = sqrt(xynorm2 + pow(u.v[0], 2));
				u.v[0] /= xynorm;
				for (int i = 1; i < row - j; i++) {
					u.v[i] = v[j + i + row*j] / xynorm;
				}
				H.eye(row - j);
				for (int i = 0; i < row - j; i++) {
					for (int k = 0; k < row - j; k++) {
						H.v[i + (row - j)*k] -= 2 * u.v[i] * u.v[k];
					}
				}
				QRs.zeros(row - j, row);
				for (int k = 0; k < row; k++) {
					for (int i = 0; i < row - j; i++) {
						QRs.v[i + (row - j) * k] = v[i + j + row*k];
						v[i + j + row*k] = zero;
					}
				}
				for (int ii = 0; ii < row - j; ii++) {
					for (int jj = 0; jj < row - j; jj++) {
						for (int kk = j; kk < row; kk++) {
							v[ii + j + row*kk] += H.v[ii + (row - j)*jj] * QRs.v[jj + (row - j)*kk];
						}
					}
				}

				for (int k = 0; k < row; k++) {
					for (int i = 0; i < row - j; i++) {
						QRs.v[i + (row - j) * k] = Q.v[i + j + row*k];
						Q.v[i + j + row*k] = zero;
					}
				}
				for (int kk = 0; kk < row; kk++) {
					for (int jj = 0; jj < row - j; jj++) {
						for (int ii = 0; ii < row - j; ii++) {
							Q.v[ii + j + row*kk] += H.v[ii + (row - j)*jj] * QRs.v[jj + (row - j)*kk];
						}
					}
				}

			}
			_T tmp;
			for (int ii = 0; ii < row; ii++) {
				for (int jj = ii + 1; jj < row; jj++) {
					tmp = Q.v[ii + row* jj];
					Q.v[ii + row* jj] = Q.v[jj + row* ii];
					Q.v[jj + row* ii] = tmp;
				}
			}
			v[row - 1 + row*(row - 2)] = zero;
		}
		// QR deconposition for a Hessenberg Matrix using Householder tansformation
		void Householder_qrdecomposition_Hessenberg(mats< _T >& Q) {
			if (row != column) {
				std::cout << "error : Householder_qrdecomposition_Hessenberg row != column" << std::endl;
				exit(1);
			}
			using std::pow;
			using std::sqrt;
			mats< _T > H, QRs;
			mats< _T > u;
			Q.eye(row);
			H.zeros(2, 2);
			QRs.zeros(2, row);
			u.zeros(2, 1);
			_T xynorm, xynorm2;
			_T zero = _T(0);
			_T one = _T(1);
			_T two = _T(2);

			for (int j = 0; j < row - 1; j++) {
				xynorm2 = pow(v[j + 1 + row*j], 2);
				xynorm = xynorm2 + pow(v[j + row*j], 2);

				if (v[j + row*j] >= 0) {
					u.v[0] = v[j + row*j] + sqrt(xynorm);
				}
				else if (v[j + row*j] < 0) {
					u.v[0] = v[j + row*j] - sqrt(xynorm);
				}
				xynorm = sqrt(xynorm2 + pow(u.v[0], 2));
				u.v[0] /= xynorm;
				u.v[1] = v[j + 1 + row*j] / xynorm;

				H.v[0] = one - two * u.v[0] * u.v[0];
				H.v[1] = -two * u.v[1] * u.v[0];
				H.v[2] = -two * u.v[0] * u.v[1];
				H.v[3] = one - two * u.v[1] * u.v[1];

				for (int k = 0; k < row; k++) {
					for (int i = 0; i < 2; i++) {
						QRs.v[i + 2 * k] = v[i + j + row*k];
						v[i + j + row*k] = zero;
					}
				}
				for (int ii = 0; ii < 2; ii++) {
					for (int jj = 0; jj < 2; jj++) {
						for (int kk = j; kk < row; kk++) {
							v[ii + j + row*kk] += H.v[ii + 2 * jj] * QRs.v[jj + 2 * kk];
						}
					}
				}

				for (int k = 0; k < row; k++) {
					for (int i = 0; i < 2; i++) {
						QRs.v[i + 2 * k] = Q.v[i + j + row*k];
						Q.v[i + j + row*k] = zero;
					}
				}

				for (int kk = 0; kk < row; kk++) {
					for (int jj = 0; jj < 2; jj++) {
						for (int ii = 0; ii < 2; ii++) {
							Q.v[ii + j + row*kk] += H.v[ii + 2 * jj] * QRs.v[jj + 2 * kk];
						}
					}
				}

			}
			_T tmp;
			for (int ii = 0; ii < row; ii++) {
				for (int jj = ii + 1; jj < row; jj++) {
					tmp = Q.v[ii + row* jj];
					Q.v[ii + row* jj] = Q.v[jj + row* ii];
					Q.v[jj + row* ii] = tmp;
				}
			}
			v[row - 1 + row*(row - 2)] = zero;
		}
		// QR deconposition for a Tridiagonal Matrix using Householder tansformation
		void Householder_qrdecomposition_Tridiag(mats< _T >& Q) {
			if (row != column) {
				std::cout << "error : Householder_qrdecomposition_Tridiag row != column" << std::endl;
				exit(1);
			}
			using std::pow;
			using std::sqrt;
			using std::min;
			mats< _T > H, QRs;
			mats< _T > u;
			Q.eye(row);
			H.zeros(2, 2);
			QRs.zeros(2, row);
			u.zeros(2, 1);
			_T xynorm, xynorm2;
			_T zero = _T(0);
			_T one = _T(1);
			_T two = _T(2);

			for (int j = 0; j < row - 1; j++) {
				xynorm2 = pow(v[j + 1 + row*j], 2);
				xynorm = xynorm2 + pow(v[j + row*j], 2);

				if (v[j + row*j] >= 0) {
					u.v[0] = v[j + row*j] + sqrt(xynorm);
				}
				else if (v[j + row*j] < 0) {
					u.v[0] = v[j + row*j] - sqrt(xynorm);
				}
				xynorm = sqrt(xynorm2 + pow(u.v[0], 2));
				u.v[0] /= xynorm;
				u.v[1] = v[j + 1 + row*j] / xynorm;

				H.v[0] = one - two * u.v[0] * u.v[0];
				H.v[1] = -two * u.v[1] * u.v[0];
				H.v[2] = -two * u.v[0] * u.v[1];
				H.v[3] = one - two * u.v[1] * u.v[1];

				for (int k = 0; k < row; k++) {
					for (int i = 0; i < 2; i++) {
						QRs.v[i + 2 * k] = v[i + j + row*k];
						v[i + j + row*k] = zero;
					}
				}

				for (int ii = 0; ii < 2; ii++) {
					for (int jj = 0; jj < 2; jj++) {
						for (int kk = j; kk < min(j + 3, row); kk++) {
							v[ii + j + row*kk] += H.v[ii + 2 * jj] * QRs.v[jj + 2 * kk];
						}
					}
				}
				for (int k = 0; k < row; k++) {
					for (int i = 0; i < 2; i++) {
						QRs.v[i + 2 * k] = Q.v[i + j + row*k];
						Q.v[i + j + row*k] = zero;
					}
				}
				for (int kk = 0; kk < row; kk++) {
					for (int jj = 0; jj < 2; jj++) {
						for (int ii = 0; ii < 2; ii++) {
							Q.v[ii + j + row*kk] += H.v[ii + 2 * jj] * QRs.v[jj + 2 * kk];
						}
					}
				}

			}
			_T tmp;
			for (int ii = 0; ii < row; ii++) {
				for (int jj = ii + 1; jj < row; jj++) {
					tmp = Q.v[ii + row* jj];
					Q.v[ii + row* jj] = Q.v[jj + row* ii];
					Q.v[jj + row* ii] = tmp;
				}
			}
			v[row - 1 + row*(row - 2)] = zero;
		}

		void eigtrisym(int itep = 1) {
			using std::min;
			using std::max;
			using std::abs;
			using std::pow;
			using std::sqrt;
			mats< _T > E, Q, C;
			mats< _T > evo;
			mats< _T > I;
			mats< _T > Tmp;
			_T zero = _T(0);
			_T two = _T(2);
			_T four = _T(4);
			_T epsilon = std::numeric_limits< _T >::epsilon();
			epsilon *= _T(n);
			_T tr, det, e1, e2, shift;
			E = (*this);
			int N = row; // N is size of E      
			int ite = 0; // iteration times


			I.eye(N);
			evo.zeros(row, 1);


			while (N > 1) {

				ite++;

//				std::cout << evo << std::endl;
//				std::cout << "\nstep : " << ite << std::endl;

				if (abs(E.v[N - 2 + E.row*(N - 1)]) < epsilon) {
//				if (abs(E(N - 2, N - 1)) < epsilon) {

					evo.v[N - 1] = E.v[(N - 1) + N * (N - 1)];
					N--;
//					std::cout << "\ndeflation n = " << N + 1 << " to " << N << std::endl;

					// iteration finish if N=1
					if (N == 1) {
						evo.v[0] = E.v[0];
//						std::cout << "\n" << std::endl;
						break;
					}

					// resize E (N+1 to N)
					Tmp = E; // temporal copy
					E.zeros(N);
					for (int ii = 0; ii < N; ii++) {
						for (int jj = 0; jj < N; jj++) {
							E.v[ii + N * jj] = Tmp.v[ii + (N + 1) * jj];
						}
					}

					I.eye(N);
				}


				// Wilkinson shift

				tr = E.v[(N - 2) + N * (N - 2)] + E.v[(N - 1) + N * (N - 1)];
				det = E.v[(N - 2) + N * (N - 2)] * E.v[(N - 1) + N * (N - 1)] - E.v[(N - 2) + N * (N - 1)] * E.v[(N - 1) + N * (N - 2)];

				// eigenvalues of 2-by-2 matrix (solve quadratic equation)
				if (tr > 0) {
					e1 = (tr + sqrt(tr * tr - four * det)) / two;
					e2 = det / e1;
				}
				else {
					e1 = (tr - sqrt(tr * tr - four * det)) / two;
					e2 = det / e1;
				}


				if (abs(e1 - E.v[(N - 1) + N * (N - 1)]) < abs(e2 - E.v[(N - 1) + N * (N - 1)])) {
					shift = e1;
				}
				else {
					shift = e2;
				}

//				std::cout << "\nshift = " << shift << std::endl;

				// QR method

				// E = E - shift * I (E -= shift * I)
				for (int ii = 0; ii < N; ii++) {
					E.v[ii + N * ii] -= shift;
				}

				// QR decomposition of E (E = Q*E)
				E.Householder_qrdecomposition(Q);

				// E = E * Q + shift * I;
				C = E;
				C.mulmm(Q,E);
//				E = E * Q;
				for (int ii = 0; ii < N; ii++) {
					E.v[ii + N * ii] += shift;
				}

			}

			// A's diagonal elements are eigenvalues    
			for (int ii = 0; ii < row; ii++) {
				for (int jj = 0; jj < row; jj++) {
					if (ii != jj) {
						v[ii + row * jj] = zero;
					}
					else v[ii + row * jj] = evo.v[ii];
				}
			}


		}
		virtual void eigsym(int itep = 1) {
			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix" << std::endl;
				exit(1);
			}
			(*this).SymTridiagonalization();
			(*this).eigtrisym(itep);
		}
		void eigsym(mats< _T >& V, int itep = 1) {
			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix" << std::endl;
				exit(1);
			}
			using std::abs;
			using std::pow;
			using std::sqrt;
			mats< _T > TriA, Eigval, Tritmp, Vtmp, tmpeig1, tmpeig2;
			mats< int > ipiv;
			_T ynorm;
			_T zero = _T(0);
			_T one = _T(1);
			_T epsilon = std::numeric_limits< _T >::epsilon();
			epsilon *= _T(n);

			_T rate = _T(0.6); // rate for setting mu

//			std::cout << "Check1" << std::endl;
			(*this).SymTridiagonalization(V);

			// transpose V
			Vtmp = V;
			for (int i = 0; i < row; i++) {
				for (int j = 0; j < row; j++) {
					V.v[i + V.row*j] = Vtmp.v[j + Vtmp.row*i];
//					V(i, j) = Vtmp(j, i);
				}
			}
			Vtmp = V;

//			std::cout << "Check2" << std::endl;
			TriA = (*this);
			(*this).eigtrisym(itep);

//			std::cout << "Check3" << std::endl;
			(*this).diag(Eigval);
			tmpeig1 = Eigval; // eigenvalues by QR method
//			std::cout << "eigenvalues by QR method" << std::endl;
//			std::cout << Eigval << std::endl;


//			std::cout << "sign of eigenvalues" << std::endl;
			mats< _T > sign;
			sign.ones(row, 1); // one : positive eigenvalue
			for (int i = 0; i < row; i++) {
				if (Eigval.v[i] < 0) {
					sign.v[i] = -one; // negative eigenvalue
				}
			}
//			std::cout << sign << std::endl;


			Eigval.abs();

			tmpeig2 = Eigval; // absolute eigenvalues         
			std::sort(Eigval.v.begin(), Eigval.v.end()); // sort (ascending order)
//			std::cout << "absolute eigenvalues" << std::endl;
//			std::cout << Eigval << std::endl;


//			std::cout << "absolute mu" << std::endl;
			mats< _T > mu;
			mu.zeros(row, 1);
			for (int i = 0; i < row; i++) {
				if (i == 0) {
					mu.v[i] = rate * Eigval.v[i];
				}
				else {
					mu.v[i] = Eigval.v[i - 1] + rate * (Eigval.v[i] - Eigval.v[i - 1]);
				}
			}
//			std::cout << mu << std::endl;


			_T d;
			mats< int > iarray;
			iarray.zeros(row, 1);
			int flag;
			for (int j = 0; j < row; j++) {
				d = abs(mu.v[j] - tmpeig2.v[0]);
				flag = 0;
				for (int i = 1; i < row; i++) {
					if (abs(mu.v[j] - tmpeig2.v[i]) < d) {
						d = abs(mu.v[j] - tmpeig2.v[i]);
						flag = i;
					}
				}
				iarray.v[j] = flag;
			}
			// iarray should be replacement of 1,2,...,row
//			std::cout << iarray << std::endl;


			// set sign of mu
			for (int i = 0; i < row; i++) {
				if (tmpeig1.v[iarray.v[i]] < 0) {
					mu.v[i] *= -one;
				}
			}
//			std::cout << "mu" << std::endl;
//			std::cout << mu << std::endl;


			// inverse iterative method
			mats < _T > A, I, Tmp, vv, y, eig;
			eig.zeros(row, 1);


			for (int l = 0; l < row; l++) {

				I.eye(row);
				A = TriA;

				// A -= mu.v[l] * I;
				for (int i = 0; i < row; i++) {
					A.v[i + row * i] -= mu.v[l];
				}

				Tmp = A;
				vv.rand(row, 1); // random vector (as initial vector)
				_T d, dold;
				d = zero;

				int k = 0; // counter of iteration

				// calculate eigenvectors (of transformed matrix)
				do {

					dold = d;
					k++;

					// Solve Ax = v
					A = Tmp;
					A.TriLudecomposition(ipiv);
					A.luonedsolve(vv, y, ipiv);

					// normalization of y
					ynorm = zero;
					for (int i = 0; i < row; i++) {
						ynorm += pow(y.v[i], 2);
					}
					ynorm = sqrt(ynorm);
					
					//y /= ynorm;
					y.divms(ynorm);

					d = abs(y.v[0]);
					for (int j = 1; j < row; j++) {
						if (abs(y.v[j])>d) {
							d = abs(y.v[j]);
						}
					}

					// renew v
					vv = y;

				} while (abs((dold - d) / d) > epsilon);

//				std::cout << "iteration : " << k << std::endl;
				_T a1, a2, lambda;
				a1 = zero;
				a2 = zero;
				mats < _T > x1, x2;
				TriA.mulmm(vv, x1);
//				x1 = TriA * v;
				x2 = vv;

				// calculate norm
				for (int i = 0; i < row; i++) {
					a1 += pow(x1.v[i], 2);
					a2 += pow(x2.v[i], 2);
				}

				// eigenvalue
				lambda = sqrt(a1 / a2);
				if (mu.v[l] < 0) {
					lambda *= -_T(1);
				}

//				std::cout << "lambda" << std::endl;
//				std::cout << lambda << "\n" << std::endl;

				eig.v[l] = lambda;

				for (int i = 0; i < row; i++) {
					V.v[i + row * l] = vv.v[i];
				}

			}

			vv = V;
			// eigenvectors of original matrix
//			V = Vtmp * V;
			Vtmp.mulmm(vv, V);
			// eigenvalues (diagonal elements)
			(*this).zeros(row);
			for (int i = 0; i < row; i++) {
				(*this).v[i + row * i] = eig.v[i];
			}

		}
		void eigsymge(mats< _T >& B, int itep = 1) {
			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix A" << std::endl;
				exit(1);
			}
			if (!B.is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix B" << std::endl;
				exit(1);
			}
			if ((*this).row != B.row || (*this).column != B.column) {
				std::cout << "error : eigsym : no muching matrix size matrix A, B" << std::endl;
				exit(1);
			}
			B.Cholesky();
			B.inv();
			mats< _T > C;
			(*this).mulmm(B, C);
			(*this) = B;
			(*this).transpose(B);
			B.mulmm(C, (*this));
			
			for (int i = 0; i < (*this).row; i++) {
				for (int j = i + 1; j < (*this).column; j++) {
					(*this).v[i + row*j] = (*this).v[j + row*i];
				}
			}

			(*this).eigsym();
		}
		void eigsymge(mats< _T >& B, mats< _T >& V, int itep = 1) {
			if (!(*this).is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix A" << std::endl;
				exit(1);
			}
			if (!B.is_symmetric()) {
				std::cout << "error : eigsym : no symmetric matrix B" << std::endl;
				exit(1);
			}
			if ((*this).row != B.row || (*this).column != B.column) {
				std::cout << "error : eigsym : no muching matrix size matrix A, B" << std::endl;
				exit(1);
			}
			B.Cholesky();
			B.inv();
			mats< _T > C;
			(*this).mulmm(B, C);
			(*this) = B;
			(*this).transpose(B);
			B.mulmm(C, (*this));

			for (int i = 0; i < (*this).row; i++) {
				for (int j = i + 1; j < (*this).column; j++) {
					(*this).v[i + row*j] = (*this).v[j + row*i];
				}
			}

			(*this).eigsym(V);
			B.transpose(C);
			B = V;
			C.mulmm(B, V);
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

		std::vector< _T > vecpointer()const {
			return v;
		}

		void eye(const int r) {
			row = r;
			column = r;
			n = r*r;
			if (r == 1) {
				type = 'S';
			}
			else {
				type = 'M';
			}
			v.resize(n);
			_T a1 = _T(1);
			_T a0 = _T(0);
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < column; j++) {
				for (int i = 0; i < row; i++) {
					if (i == j) {
						v[i + row*j] = a1;
					}
					else {
						v[i + row*j] = a0;
					}
				}
			}
			v.shrink_to_fit();
		}
		void ones(const int i) {
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
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < n; j++) {
				v[j] = _T(1);
			}
			v.shrink_to_fit();
		}
		void ones(const int r, const int c) {
			_T a0 = _T(1);
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
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < column; j++) {
				for (int i = 0; i < row; i++) {
					v[i + row*j] = a0;
				}
			}
			v.shrink_to_fit();
		}
		void zeros(const int i) {
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
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < n; j++) {
				v[j] = _T(0);
			}
			v.shrink_to_fit();
		}
		void zeros(const int r, const int c) {
			_T a0 = _T(0);
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
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < column; j++) {
				for (int i = 0; i < row; i++) {
					v[i + row*j] = a0;
				}
			}
			v.shrink_to_fit();
		}
		void rand(const int i) {
#ifdef rand_debug
			std::cout << "warning : A data type of the function rand is double." << std::endl;
#endif
			std::random_device seed_gen;
			std::mt19937 engine(seed_gen());
			std::normal_distribution< double > dist(0.0, 1.0);
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
				v[j] = _T(dist(engine));
			}
			v.shrink_to_fit();
		}
		void rand(const int r, const int c) {
#ifdef rand_debug
			std::cout << "warning : A data type of the function rand is double." << std::endl;
#endif
			std::random_device seed_gen;
			std::mt19937 engine(seed_gen());
			std::normal_distribution< double > dist(0.0, 1.0);
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
					v[i + row*j] = _T(dist(engine));
				}
			}
			v.shrink_to_fit();
		}

		void resize(const int i, const int j) {
			int nn = i*j;
			int orow, ocolumn, on;
			_T a0 = _T(0);
			orow = row;
			ocolumn = column;
			on = n;

			if (row > i || column > j) {
				std::cout << "error : resize : row > i || column > j" << std::endl;
				exit(1);
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
						v[ii + orow*jj] = a0;
					}
				}
			}
			v.shrink_to_fit();
		}
		
		void clear() {
			(*this).v.clear();
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
		}

		bool is_symmetric() {
			if ((*this).row != (*this).column){
				std::cout << "error : is_symmetric : row != column" << std::endl;
			}

			for (int i = 0; i < (*this).row; i++) {
				for (int j = i + 1; j < (*this).column; j++) {
					if ((*this).v[i + (*this).row*j] != (*this).v[j + (*this).row*i]) {
						return false;
					}
				}
			}
			return true;
		}

		friend std::ostream& operator<<(std::ostream& os, mats< _T >& A) {
			return A.display(os);
		}
	};
}
#endif // VCP_MATS_HPP