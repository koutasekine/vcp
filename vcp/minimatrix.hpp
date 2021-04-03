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

#ifndef VCP_MINIMATRIX_HPP
#define VCP_MINIMATRIX_HPP

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
#include <vcp/matrix.hpp>

namespace vcp{
	template <typename _T> class minimats {
	public:
		int row;
		int column;
		int n;
		char type;      //'N':NULL  'S':Scala  'R' Row Vector 'C':Column Vector 'M':Matrix
		std::vector< _T > v;

		minimats< _T >() {
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
		}
		virtual ~minimats< _T >() = default;
		minimats< _T >(const minimats< _T >&) = default;
		minimats< _T >(minimats< _T >&&) = default;
		minimats< _T >& operator=(const minimats< _T >& A) = default;
		minimats< _T >& operator=(minimats< _T >&& A) = default;

		// A = A+B
		void addmm(const minimats< _T >& B) {
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
		}

		// A = A-B
		void subsmmA(const minimats< _T >& B) {
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
		void subsmmB(const minimats< _T >& B) {
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
		virtual void mulmm(const minimats< _T >& B, minimats< _T >& c)const {
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

		void diag(minimats< _T >& B)const {
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
		void transpose(minimats< _T >& B)const {
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
		void horzcat(const minimats< _T >& B, minimats< _T >& C)const {
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
		void vercat(const minimats< _T >& B, minimats< _T >& C)const {
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

		void submat(minimats< _T >& B, const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
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
		
		int length()const {
			using std::max;
			return max(column, row);
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

		friend std::ostream& operator<<(std::ostream& os, const minimats< _T >& A) {
			return A.display(os);
		}
	};

	template <typename _T> class matrix< _T, vcp::minimats< _T > > : protected vcp::minimats< _T > {
	public:
		matrix< _T, vcp::minimats< _T > >() {
			(*this).row = 0;
			(*this).column = 0;
			(*this).n = 0;
			(*this).type = 'N';
		}
		~matrix< _T, vcp::minimats< _T > >() = default;
		matrix< _T, vcp::minimats< _T > >(const matrix< _T, vcp::minimats< _T > >&) = default;
		matrix< _T, vcp::minimats< _T > >(matrix< _T, vcp::minimats< _T > >&&) = default;
		matrix< _T, vcp::minimats< _T > >& operator=(const matrix< _T, vcp::minimats< _T > >& A) = default;
		matrix< _T, vcp::minimats< _T > >& operator=(matrix< _T, vcp::minimats< _T > >&& A) = default;

		_T& operator () (const int i) {
			return (*this).v[i];
		}
		_T operator () (const int i) const {
			return (*this).v[i];
		}
		_T& operator () (const int i, const int j) {
			if ((*this).type == 'R') {
				return (*this).v[j];
			}
			else {
				return (*this).v[i + (*this).row*j];
			}
		}
		_T operator () (const int i, const int j)const {
			if ((*this).type == 'R') {
				return (*this).v[j];
			}
			else {
				return (*this).v[i + (*this).row*j];
			}
		}
	/*
		_T& operator [] (const int i) {
			return (*this).v[i];
		}
		_T operator [] (const int i) const {
			return (*this).v[i];
		}
		*/
		
//***************** Operator Overload *****************//
		friend matrix< _T, vcp::minimats< _T > > operator+(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			matrix< _T, vcp::minimats< _T > > C;
			C = A;
			C.addmm(B);
			return std::move(C);
		}
		friend matrix< _T, vcp::minimats< _T > > operator+(matrix< _T, vcp::minimats< _T > >&& A, const matrix< _T, vcp::minimats< _T > >& B) {
			A.addmm(B);
			return std::move(A);
		}
		friend matrix< _T, vcp::minimats< _T > > operator+(const matrix< _T, vcp::minimats< _T > >& A, matrix< _T, vcp::minimats< _T > >&& B) {
			B.addmm(A);
			return std::move(B);
		}
		friend matrix< _T, vcp::minimats< _T > > operator+(matrix< _T, vcp::minimats< _T > >&& A, matrix< _T, vcp::minimats< _T > >&& B) {
			A.addmm(B);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator+(const _Tm a, const matrix< _T, vcp::minimats< _T > >& B) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.addsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator+(const _Tm a, matrix< _T, vcp::minimats< _T > >&& B) {
			_T Ta = _T(a);
			B.addsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator+(const matrix< _T, vcp::minimats< _T > >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.addms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator+(matrix< _T, vcp::minimats< _T > >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.addms(Ta);
			return std::move(B);
		}
		friend matrix< _T, vcp::minimats< _T > > operator+(const matrix< _T, vcp::minimats< _T > >& A) {
			//		A.plusm();
			return A;
		}
		friend matrix< _T, vcp::minimats< _T > > operator+(matrix< _T, vcp::minimats< _T > >&& A) {
			//		A.plusm();
			return std::move(A);
		}

		friend matrix< _T, vcp::minimats< _T > >& operator+=(matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			A.addmm(B);
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > >& >::type operator+=(matrix< _T, vcp::minimats< _T > >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.addms(Ta);
			return A;
		}

		friend matrix< _T, vcp::minimats< _T > > operator-(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			matrix< _T, vcp::minimats< _T > > C;
			C = A;
			C.subsmmA(B);
			return std::move(C);
		}
		friend matrix< _T, vcp::minimats< _T > > operator-(matrix< _T, vcp::minimats< _T > >&& A, const matrix< _T, vcp::minimats< _T > >& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		friend matrix< _T, vcp::minimats< _T > > operator-(const matrix< _T, vcp::minimats< _T > >& A, matrix< _T, vcp::minimats< _T > >&& B) {
			B.subsmmB(A);
			return std::move(B);
		}
		friend matrix< _T, vcp::minimats< _T > > operator-(matrix< _T, vcp::minimats< _T > >&& A, matrix< _T, vcp::minimats< _T > >&& B) {
			A.subsmmA(B);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator-(const _Tm a, const matrix< _T, vcp::minimats< _T > >& B) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.subssm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator-(const _Tm a, matrix< _T, vcp::minimats< _T > >&& B) {
			_T Ta = _T(a);
			B.subssm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator-(const matrix< _T, vcp::minimats< _T > >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.subsms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator-(matrix< _T, vcp::minimats< _T > >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.subsms(Ta);
			return std::move(B);
		}
		friend matrix< _T, vcp::minimats< _T > > operator-(const matrix< _T, vcp::minimats< _T > >& A) {
			matrix< _T, vcp::minimats< _T > > C;
			C = A;
			C.minusm();
			return std::move(C);
		}
		friend matrix< _T, vcp::minimats< _T > > operator-(matrix< _T, vcp::minimats< _T > >&& A) {
			A.minusm();
			return std::move(A);
		}

		friend matrix< _T, vcp::minimats< _T > >& operator-=(matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			A.subsmmA(B);
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > >& >::type operator-=(matrix< _T, vcp::minimats< _T > >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.subsms(Ta);
			return A;
		}

		friend matrix< _T, vcp::minimats< _T > > operator*(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			matrix< _T, vcp::minimats< _T > > C;
			A.mulmm(B, C);
			return std::move(C);
		}
		friend matrix< _T, vcp::minimats< _T > > operator*(matrix< _T, vcp::minimats< _T > >&& A, const matrix< _T, vcp::minimats< _T > >& B) {
			matrix< _T, vcp::minimats< _T > > C;
			A.mulmm(B, C);
			A = std::move(C);
			return std::move(A);
		}
		friend matrix< _T, vcp::minimats< _T > > operator*(const matrix< _T, vcp::minimats< _T > >& A, matrix< _T, vcp::minimats< _T > >&& B) {
			matrix< _T, vcp::minimats< _T > > C;
			A.mulmm(B, C);
			B = std::move(C);
			return std::move(B);
		}
		friend matrix< _T, vcp::minimats< _T > > operator*(matrix< _T, vcp::minimats< _T > >&& A, matrix< _T, vcp::minimats< _T > >&& B) {
			matrix< _T, vcp::minimats< _T > > C;
			A.mulmm(B, C);
			A = std::move(C);
			return std::move(A);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator*(const _Tm a, const matrix< _T, vcp::minimats< _T > >& B) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.mulsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator*(const _Tm a, matrix< _T, vcp::minimats< _T > >&& B) {
			_T Ta = _T(a);
			B.mulsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator*(const matrix< _T, vcp::minimats< _T > >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.mulms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator*(matrix< _T, vcp::minimats< _T > >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.mulms(Ta);
			return std::move(B);
		}

		friend matrix< _T, vcp::minimats< _T > >& operator*=(matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			A = A * B;
			return A;
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > >& >::type operator*=(matrix< _T, vcp::minimats< _T > >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.mulms(Ta);
			return A;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator/(const _Tm a, const matrix< _T, vcp::minimats< _T > >& B) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.divsm(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator/(const _Tm a, matrix< _T, vcp::minimats< _T > >&& B) {
			_T Ta = _T(a);
			B.divsm(Ta);
			return std::move(B);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator/(const matrix< _T, vcp::minimats< _T > >& B, const _Tm a) {
			_T Ta = _T(a);
			matrix< _T, vcp::minimats< _T > > C;
			C = B;
			C.divms(Ta);
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > > >::type operator/(matrix< _T, vcp::minimats< _T > >&& B, const _Tm a) {
			_T Ta = _T(a);
			B.divms(Ta);
			return std::move(B);
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, matrix< _T, vcp::minimats< _T > >& >::type operator/=(matrix< _T, vcp::minimats< _T > >& A, const _Tm& a) {
			_T Ta = _T(a);
			A.divms(Ta);
			return A;
		}

		matrix< _T, vcp::minimats< _T > > submatrix(const std::initializer_list<int>& list1, const std::initializer_list<int>& list2) const {
			matrix< _T, vcp::minimats< _T > > A;
			(*this).submat(A, list1, list2);
			return std::move(A);
		}

		int elementsize()const { return (*this).n; }
		int columnsize()const { return (*this).column; }
		int rowsize()const { return (*this).row; }
		char matstype()const {
			return (*this).type;
		}
		std::vector< _T > vecpointer()const {
			return (*this).v;
		}

		void zeros(const int i) { minimats< _T >::zeros(i); }
		void zeros(const int r, const int c) { minimats< _T >::zeros(r, c); }
		void ones(const int i) { minimats< _T >::ones(i); }
		void ones(const int r, const int c) { minimats< _T >::ones(r, c); }
		void resize(const int i, const int j) { minimats< _T >::resize(i, j); }
		void clear() { minimats< _T >::clear(); }

		friend matrix< _T, vcp::minimats< _T > > transpose(const matrix< _T, vcp::minimats< _T > >& A) {
			matrix< _T, vcp::minimats< _T > > C;
			A.transpose(C);
			return std::move(C);
		}

		//Matlab C = [A,B];
		friend matrix< _T, vcp::minimats< _T > > horzcat(const matrix< _T, vcp::minimats< _T > >& A) {
			return A;
		}
		friend matrix< _T, vcp::minimats< _T > > horzcat(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			matrix< _T, vcp::minimats< _T > > C;
			A.horzcat(B, C);
			return std::move(C);
		}
		template<typename... Args> friend matrix< _T, vcp::minimats< _T > > horzcat(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B, const Args&... args) {
			matrix< _T, vcp::minimats< _T > > C;
			A.horzcat(B, C);
			return std::move(horzcat(C, args...));
		}
		
		//Matlab C = [A;B];
		friend matrix< _T, vcp::minimats< _T > > vercat(const matrix< _T, vcp::minimats< _T > >& A) {
			return A;
		}
		friend matrix< _T, vcp::minimats< _T > > vercat(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B) {
			matrix< _T, vcp::minimats< _T > > C;
			A.vercat(B, C);
			return std::move(C);
		}
		template<typename... Args> friend matrix< _T, vcp::minimats< _T > > vercat(const matrix< _T, vcp::minimats< _T > >& A, const matrix< _T, vcp::minimats< _T > >& B, const Args&... args) {
			matrix< _T, vcp::minimats< _T > > C;
			A.vercat(B, C);
			return std::move(vercat(C, args...));
		}

		friend int length(matrix< _T, vcp::minimats< _T > >& A) {
			return A.length();
		}


		//**************** display function ***************//
		friend std::ostream& operator<<(std::ostream& os, matrix< _T, const vcp::minimats< _T > >& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, matrix< _T, vcp::minimats< _T > >&& A) {
			return A.display(os);
		}
	};
}
#endif // VCP_MATS_HPP