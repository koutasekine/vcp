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

#ifndef VCP_FOURIER_SERIES_HPP
#define VCP_FOURIER_SERIES_HPP

#ifdef VCP_NOMP
#ifndef VCP_FOURIER_NOMP
#define VCP_FOURIER_NOMP
#endif
#endif

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <initializer_list>

#include <vcp/sincos_class.hpp>

// Basis {1, sinx, cosx, sin2x, cos2x, ... , sinnx, cosnx}
// phi0 = 1
// phi1 = sinx
// phi2 = cosx
// ...

namespace vcp {
    template <typename _T>
    class fourier_series{
protected:
		int n;
		_T a0;
		std::vector< vcp::sincos< _T > > c;

		// A = A + B
		void AplusB(const fourier_series< _T >& B){
			this->a0 += B.a0;			
			if (this->n >= B.n){
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < B.n; i++){
					this->c[i] += B.c[i];
				}
			}
			else{

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < this->n; i++){
					this->c[i] += B.c[i];
				}

				for (int i = this->n; i < B.n; i++){
					this->c.push_back( B.c[i] );
				}
				this->n = B.n;
                this->c.shrink_to_fit();
			}
		}

		// A = A - B
		void AminusB(const fourier_series< _T >& B){
			this->a0 -= B.a0;			
			if (this->n >= B.n){
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < B.n; i++){
					this->c[i] -= B.c[i];
				}
			}
			else{

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < this->n; i++){
					this->c[i] -= B.c[i];
				}

				for (int i = this->n; i < B.n; i++){
					this->c.push_back( -B.c[i] );
				}
				this->n = B.n;
                this->c.shrink_to_fit();
			}
		}

		void minusAplusB(const fourier_series< _T >& B){
			this->a0 = - this->a0 + B.a0;			
			if (this->n >= B.n){
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < B.n; i++){
					this->c[i] = - this->c[i] + B.c[i];
				}
			}
			else{

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
				for (int i = 0; i < this->n; i++){
					this->c[i] = - this->c[i] + B.c[i];
				}

				for (int i = this->n; i < B.n; i++){
					this->c.push_back( B.c[i] );
				}
				this->n = B.n;
                this->c.shrink_to_fit();
			}
		}

		void mul(const fourier_series< _T >& B, fourier_series< _T >& C) const {
			C.zeros((*this).n + B.n);
			C.a0 = (*this).a0*B.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for ( int i = 0; i < (*this).n; i++){
				C.c[i] = (*this).c[i]*B.a0;
			}
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for ( int j = 0; j < B.n; j++){
				C.c[j] += (*this).a0*B.c[j];
			}
			
			// sin(a) * sin(b) =-(cos(a+b) - cos(a-b))/2
			// sin(a) * cos(b) = (sin(a+b) + sin(a-b))/2
			// cos(a) * cos(b) = (cos(a+b) + cos(a-b))/2
			for ( int i = 1; i <= (*this).n; i++){
				for ( int j = 1; j <= B.n; j++){
					int ii = i - 1;
					int jj = j - 1;
					int ipj = i + j - 1;
					int imj = i - j - 1;
					if ( i > j ){
						// ai*sin(i)*aj*sin(j) = -ai*aj*( cos(i+j) - cos(i-j) )/2
						C.c[ipj].b -= ((*this).c[ii].a*B.c[jj].a)/_T(2);
						C.c[imj].b += ((*this).c[ii].a*B.c[jj].a)/_T(2);


						// ai*sin(i)*bj*cos(j) = ai*bj*(sin(i+j) + sin(i-j))/2
						C.c[ipj].a += ((*this).c[ii].a*B.c[jj].b)/_T(2);
						C.c[imj].a += ((*this).c[ii].a*B.c[jj].b)/_T(2);


						// bi*cos(i)*aj*sin(j) = bi*aj*(sin(i+j) - sin(i-j))/2
						C.c[ipj].a += ((*this).c[ii].b*B.c[jj].a)/_T(2);
						C.c[imj].a -= ((*this).c[ii].b*B.c[jj].a)/_T(2);


						// bi*cos(i)*bj*cos(j) = bi*bj*(cos(i+j) + cos(i-j))/2
						C.c[ipj].b += ((*this).c[ii].b*B.c[jj].b)/_T(2);
						C.c[imj].b += ((*this).c[ii].b*B.c[jj].b)/_T(2);

					}
					else if( i == j ){
						// sin(0)=0, cos(0)=1
						// ai*sin(i) * aj*sin(j) = -ai*aj*( cos(i+j) - 1 )/2 = ai*aj/2 - ai*aj*cos(i+j)/2
						C.c[ipj].b -= ((*this).c[ii].a*B.c[jj].a)/_T(2);
						C.a0       += ((*this).c[ii].a*B.c[jj].a)/_T(2);

						// ai*sin(i)*bj*cos(j) = ai*bj*(sin(i+j) + sin(0))/2 = ai*bj*sin(i+j)/2
						C.c[ipj].a += ((*this).c[ii].a*B.c[jj].b)/_T(2);

						// bi*cos(i)*aj*sin(j) = bi*aj*(sin(i+j) - sin(0))/2 = bi*aj*sin(i+j)/2
						C.c[ipj].a += ((*this).c[ii].b*B.c[jj].a)/_T(2);

						// bi*cos(i)*bj*cos(j) = bi*bj*(cos(i+j) + 1)/2 = bi*bj/2 + bi*bj*cos(i+j)/2
						C.c[ipj].b += ((*this).c[ii].b*B.c[jj].b)/_T(2);
						C.a0       += ((*this).c[ii].b*B.c[jj].b)/_T(2);

					}
					else if( i < j ){
						// sin(-a)=-sin(a), cos(-a)=cos(a)
						// ai*sin(i)*aj*sin(j) = -ai*aj*( cos(i+j) - cos(i-j) )/2 = -ai*aj*( cos(i+j) - cos(j-i) )/2
						C.c[ipj].b   -= ((*this).c[ii].a*B.c[jj].a)/_T(2);
						C.c[j-i-1].b += ((*this).c[ii].a*B.c[jj].a)/_T(2);


						// ai*sin(i)*bj*cos(j) = ai*bj*(sin(i+j) + sin(i-j))/2 = ai*bj*(sin(i+j) - sin(j-i))/2
						C.c[ipj].a   += ((*this).c[ii].a*B.c[jj].b)/_T(2);
						C.c[j-i-1].a -= ((*this).c[ii].a*B.c[jj].b)/_T(2);


						// bi*cos(i)*aj*sin(j) = bi*aj*(sin(i+j) + sin(j-i))/2
						C.c[ipj].a += ((*this).c[ii].b*B.c[jj].a)/_T(2);
						C.c[j-i-1].a += ((*this).c[ii].b*B.c[jj].a)/_T(2);


						// bi*cos(i)*bj*cos(j) = bi*bj*(cos(i+j) + cos(i-j))/2 = bi*bj*(cos(i+j) + cos(j-i))/2
						C.c[ipj].b   += ((*this).c[ii].b*B.c[jj].b)/_T(2);
						C.c[j-i-1].b += ((*this).c[ii].b*B.c[jj].b)/_T(2);
					}
				}
			}			
		}


		void mulsin( const int& N, fourier_series< _T >& C ) const {
			C.zeros((*this).n + N);

			if ( N <= 0 ){
				std::cout << "Error: fourier_series: mulsin: N < 0..." << std::endl;
				exit(0);
			}
			else if( N > 0 ){
				// a0*sin(N)
				C.c[N-1].a = (*this).a0;

				// N = 1 : phi1 = sin(1t)
				// N = 2 : phi3 = sin(2t)
				// N = 3 : phi5 = sin(3t)
				// N = 4 : phi7 = sin(4t)
				for ( int i = 1; i <= (*this).n; i++){
					int ii = i - 1;
					int NN = N - 1;
					int ipN = i + N - 1;
					int imN = i - N - 1;
					_T tmp;

					if ( i > N ){
						// ai*sin(i) * sin(N) =-(ai*cos(i+N) - ai*cos(i-N))/2
						tmp = (*this).c[ii].a/_T(2);
						C.c[ipN].b -= tmp;
						C.c[imN].b += tmp;

						// sin(N) * bi*cos(i) = (bi*sin(i+N) + bi*sin(N-i))/2 = (bi*sin(i+N) - bi*sin(i-N)))/2
						tmp = (*this).c[ii].b/_T(2);
						C.c[ipN].a += tmp;
						C.c[imN].a -= tmp;
					}
					else if(i == N){
						// ai*sin(i) * sin(N) =-(ai*cos(i+N) - ai*cos(i-N))/2 = -(ai*cos(i+N) - ai)/2
						tmp = (*this).c[ii].a/_T(2);
						C.c[ipN].b -= tmp;
						C.a0 += tmp;

						// sin(N) * bi*cos(i) = (bi*sin(i+N) + bi*sin(N-i))/2 = bi*sin(i+N) /2
						C.c[ipN].a += (*this).c[ii].b/_T(2);

					}
					else if( i < N){
						// ai*sin(i) * sin(N) =-(ai*cos(i+N) - ai*cos(i-N))/2 = -(ai*cos(i+N) - ai*cos(N-i))/2
						tmp = (*this).c[ii].a/_T(2);
						C.c[ipN].b -= tmp;
						C.c[N-i-1].b += tmp;

						// sin(N) * bi*cos(i) = (bi*sin(i+N) + bi*sin(N-i))/2
						tmp = (*this).c[ii].b/_T(2);
						C.c[ipN].a += tmp;
						C.c[N-i-1].a += tmp;
					}
				}
			}	
		}

		
		void mulcos( const int& N, fourier_series< _T >& C ) const {
			C.zeros((*this).n + N);


			if ( N <= 0 ){
				std::cout << "Error: fourier_series: mulcos: N < 0..." << std::endl;
				exit(0);
			}
			else if( N > 0 ){
				// a0*sin(N)
				C.c[N-1].b = (*this).a0;

				// N = 1 : phi2 = cos(1t)
				// N = 2 : phi4 = cos(2t)
				// N = 3 : phi6 = cos(3t)
				// N = 4 : phi8 = cos(4t)
				// sin(a) * cos(b) = (sin(a+b) + sin(a-b))/2
				// cos(a) * cos(b) = (cos(a+b) + cos(a-b))/2
				for ( int i = 1; i <= (*this).n; i++){
					int ii = i - 1;
					int NN = N - 1;
					int ipN = i + N - 1;
					int imN = i - N - 1;
					_T tmp;

					if ( i > N ){
						// ai*sin(i) * cos(N) = ai*(sin(i+N) + sin(i-N))/2
						tmp = (*this).c[ii].a/_T(2);
						C.c[ipN].a += tmp;
						C.c[imN].a += tmp;

						// bi*cos(i) * cos(N) = bi*(cos(i+N) + cos(i-N))/2
						tmp = (*this).c[ii].b/_T(2);
						C.c[ipN].b += tmp;
						C.c[imN].b += tmp;
					}
					else if(i == N){
						// ai*sin(i) * cos(N) = ai*(sin(i+N) + sin(i-N))/2 = ai*sin(i+N)/2
						C.c[ipN].a += (*this).c[ii].a/_T(2);

						// bi*cos(i) * cos(N) = bi*(cos(i+N) + cos(i-N))/2 = bi/2*(cos(i+N) + 1)
						tmp = (*this).c[ii].b/_T(2);
						C.c[ipN].b += tmp;
						C.a0 += tmp;

					}
					else if( i < N){
						// ai*sin(i) * cos(N) = ai*(sin(i+N) + sin(i-N))/2 = ai*(sin(i+N) - sin(N-i))/2
						tmp = (*this).c[ii].a/_T(2);
						C.c[ipN].a += tmp;
						C.c[N-i-1].a -= tmp;

						// bi*cos(i) * cos(N) = bi*(cos(i+N) + cos(i-N))/2 = bi*(cos(i+N) + cos(N-i))/2
						tmp = (*this).c[ii].b/_T(2);
						C.c[ipN].b += tmp;
						C.c[N-i-1].b += tmp;
					}
				}
			}	
			
		}


		template <typename _Tm>
        typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type time_delay( const _Tm& tau){
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < n; j++) {
				(*this).c[j].time_delay((j+1)*tau);
			}
        }


public:
        fourier_series< _T >() {
			(*this).n = 0;
		}
		~fourier_series< _T >() = default;
		fourier_series< _T >(const fourier_series< _T >&) = default;
		fourier_series< _T >(fourier_series< _T >&&) = default;
		fourier_series< _T >& operator=(const fourier_series< _T >& A) = default;
		fourier_series< _T >& operator=(fourier_series< _T >&& A) = default;

		template <typename _Tm> typename std::enable_if<std::is_constructible< _T, _Tm >::value, void >::type set_a0( _Tm a ){
			(*this).a0 = _T(a);
		}

		template <typename _Tm> typename std::enable_if<std::is_constructible< _T, _Tm >::value, void >::type set_sinm( _Tm a, int m ){
			(*this).c[m-1].set_a( _T(a) );
		}

		template <typename _Tm> typename std::enable_if<std::is_constructible< _T, _Tm >::value, void >::type set_cosm( _Tm a, int m ){
			(*this).c[m-1].set_b( _T(a) );
		}

		_T get_sin(const int& NN) const {
			if ( NN <= this->n ){
				return (*this).c[NN - 1].a;
			}
			else{
				return _T(0);
			}
		}

		_T get_cos(const int& NN) const {
			if ( NN <= this->n ){
				return (*this).c[NN - 1].b;
			}
			else{
				return _T(0);
			}
		}

		_T get_a0(void) const {
			return (*this).a0;
		}

		int get_n( void ) const {
			return (*this).n;
		}

        void zeros(const int NN) {
			(*this).n = NN;
			(*this).a0 = _T(0);
			(*this).c.resize(NN);
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < (*this).n; j++) {
				(*this).c[j].set_ab(_T(0), _T(0));
			}
			(*this).c.shrink_to_fit();
		}

		void ones(const int NN) {
			(*this).n = NN;
			(*this).a0 = _T(1);
			(*this).c.resize(n);
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int j = 0; j < n; j++) {
				(*this).c[j].set_ab(_T(1), _T(1));
			}
			(*this).c.shrink_to_fit();
		}

        friend fourier_series< _T > operator+(const fourier_series< _T >& A) {
			return A;
		}
        
		friend fourier_series< _T > operator+(fourier_series< _T >&& A) {
			return std::move(A);
		}

        friend fourier_series< _T > operator+(const fourier_series< _T >& A, const fourier_series< _T >& B) {
			fourier_series< _T > C;
			C = A;
			C.AplusB(B);
			return std::move(C);
		}

		friend fourier_series< _T > operator+(fourier_series< _T >&& A, const fourier_series< _T >& B) {
			A.AplusB(B);
			return std::move(A);
		}

		friend fourier_series< _T > operator+(const fourier_series< _T >& A, fourier_series< _T >&& B) {
			B.AplusB(A);
			return std::move(B);
		}

		friend fourier_series< _T > operator+(fourier_series< _T >&& A, fourier_series< _T >&& B) {
			A.AplusB(B);
			return std::move(A);
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator+( const fourier_series<_T>& A, const _Tm& a ){
			fourier_series<_T> C;
			_T Ta = _T(a);
			C = A;
			C.a0 += Ta;
			return C;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator+( fourier_series<_T>&& A, const _Tm& a ){
			_T Ta = _T(a);
			A.a0 += Ta;
			return std::move(A);
		}


		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator+( const _Tm& a, const fourier_series<_T>& A ){
			fourier_series<_T> C;
			_T Ta = _T(a);
			C = +A;
			C.a0 += Ta;
			return C;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator+( const _Tm& a, fourier_series<_T>&& A ){
			_T Ta = _T(a);
			A.a0 += Ta;
			return std::move(A);
		}


        friend fourier_series< _T >& operator+=(fourier_series< _T >& A, const fourier_series< _T >& B) {
			A.AplusB(B);
			return A;
		}

        friend fourier_series< _T > operator-(const fourier_series< _T >& A) {
            fourier_series< _T > B;
			B.zeros(A.n);
			B.a0 = -A.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < A.n; i++){
				B.c[i] = -A.c[i];
			}
			return B;
		}
		
		friend fourier_series< _T > operator-(fourier_series< _T >&& A) {
			A.a0 = -A.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < A.n; i++){
				A.c[i] = -A.c[i];
			}
			return std::move(A);
		}

        friend fourier_series< _T > operator-(const fourier_series< _T >& A, const fourier_series< _T >& B) {
			fourier_series< _T > C;
			C = A;
			C.AminusB(B);
			return std::move(C);
		}

		friend fourier_series< _T > operator-(fourier_series< _T >&& A, const fourier_series< _T >& B) {
			A.AminusB(B);
			return std::move(A);
		}

		friend fourier_series< _T > operator-(const fourier_series< _T >& A, fourier_series< _T >&& B) {
			B.minusAplusB(A);
			return std::move(B);
		}

		friend fourier_series< _T > operator-(fourier_series< _T >&& A, fourier_series< _T >&& B) {
			A.AminusB(B);
			return std::move(A);
		}

		
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator-( const fourier_series<_T>& A, const _Tm& a ){
			fourier_series<_T> C;
			_T Ta = _T(a);
			C = A;
			C.a0 -= Ta;
			return C;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator-( fourier_series<_T>&& A, const _Tm& a ){
			_T Ta = _T(a);
			A.a0 -= Ta;
			return std::move(A);
		}


		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator-( const _Tm& a, const fourier_series<_T>& A ){
			fourier_series<_T> C;
			_T Ta = _T(a);
			C = -A;
			C.a0 += Ta;
			return C;
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator-( const _Tm& a, fourier_series<_T>&& A ){
			_T Ta = _T(a);
			A.a0 -= Ta;
			return std::move(-A);
		}

		friend fourier_series< _T > operator*(const fourier_series< _T >& A, const fourier_series< _T >& B) {
			fourier_series< _T > C;
			A.mul(B, C);
			return std::move(C);
		}

		friend fourier_series< _T > operator*(const fourier_series< _T >& A, fourier_series< _T >&& B) {
			fourier_series< _T > C;
			A.mul(B, C);
			B = std::move(C);
			return std::move(B);
		}

		friend fourier_series< _T > operator*(fourier_series< _T >&& A, const fourier_series< _T >& B) {
			fourier_series< _T > C;
			A.mul(B, C);
			A = std::move(C);
			return std::move(A);
		}

		friend fourier_series< _T > operator*(fourier_series< _T >&& A, fourier_series< _T >&& B) {
			fourier_series< _T > C;
			A.mul(B, C);
			A = std::move(C);
			return std::move(A);
		}

        template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator*(const fourier_series< _T >& B, const _Tm& a) {
			_T Ta = _T(a);
			fourier_series< _T > C;
			C = B;
			C.a0 *= Ta;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < C.n; i++){
				C.c[i] = C.c[i]*Ta;
			}
			return std::move(C);
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator*(fourier_series< _T >&& B, const _Tm& a) {
			_T Ta = _T(a);
			B.a0 *= Ta;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < B.n; i++){
				B.c[i] = B.c[i]*Ta;
			}
			return std::move(B);
		}

        template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator*(const _Tm& a, const fourier_series< _T >& B) {
			_T Ta = _T(a);
			fourier_series< _T > C;
			C = B;
			C.a0 *= Ta;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < C.n; i++){
				C.c[i] = Ta*C.c[i];
			}
			return std::move(C);
		}

		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator*(const _Tm& a, fourier_series< _T >&& B) {
			_T Ta = _T(a);
			B.a0 *= Ta;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < B.n; i++){
				B.c[i] = Ta*B.c[i];
			}
			return std::move(B);
		}

        template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator/(const fourier_series< _T >& B, const _Tm& a) {
			_T Ta = _T(a);
			fourier_series< _T > C;
			C = B;
			C.a0 /= Ta;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < C.n; i++){
				C.c[i] = C.c[i]/Ta;
			}
			return std::move(C);
		}
		
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type operator/(fourier_series< _T >&& B, const _Tm& a) {
			_T Ta = _T(a);
			B.a0 /= Ta;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < B.n; i++){
				B.c[i] = B.c[i]/Ta;
			}
			return std::move(B);
		}

		template <typename _Tm> typename std::enable_if<std::is_constructible< _T, _Tm >::value, fourier_series< _T > >::type delay( const _Tm& tau )const{
			fourier_series< _T > B;
			B = (*this);
			B.time_delay(tau);
			return B;
		}

		fourier_series< _T > mul_sin( const int& N )const{
			fourier_series< _T > B;
			(*this).mulsin(N, B);
			return B;
		}

		fourier_series< _T > mul_cos( const int& N )const{
			fourier_series< _T > B;
			(*this).mulcos(N, B);
			return B;
		}

		fourier_series< _T > diff( void ){
			fourier_series< _T > B;

			B = (*this);
			B.a0 = _T(0);
			
			for (int i = 1; i <= B.n; i++){
				_T tmp = B.c[i-1].a;
				B.c[i-1].a = -_T(i)*B.c[i-1].b;
				B.c[i-1].b = _T(i)*tmp;
			}

			return B;
		}

		fourier_series< _T > Pm(const int& m)const{
			fourier_series< _T > B;
			B.zeros(m);
			B.a0 = (*this).a0;
			for (int i = 0; i < m; i++){
				B.c[i] = (*this).c[i];
			}
			return B;
		}

		fourier_series< _T > I_minus_Pm(const int& m)const{
			fourier_series< _T > B;
			B = (*this);
			B.a0 = _T(0);
			for (int i = 0; i < m; i++){
				B.c[i].a = _T(0);
				B.c[i].b = _T(0);
			}
			return B;
		}
		
		void display(){
			std::cout << (*this).a0;
			for (int i = 0; i < (*this).n; i++){
				std::cout << " + " << (*this).c[i].a << "sin(" << i+1 << "t)";
				std::cout << " + " << (*this).c[i].b << "cos(" << i+1 << "t)";
			}
			std::cout << "\n" << std::endl;
		}

		_T L2norm(void) const {
			using std::pow;
			using std::sqrt;
			_T s = pow(_T(2)*a0, 2)/_T(2);
			for (int i = 0; i < this->n; i++){
				s += pow(this->c[i].a, 2);
				s += pow(this->c[i].b, 2);
			}
			return sqrt(s);
		}

		template <typename _Tm>
		typename std::enable_if<std::is_constructible< _Tm, _T >::value, _Tm >::type value( const _Tm& t) const {
			int n = this->get_n();
			_Tm s = this->get_a0();
			for (int j = 1; j <= n; j++){
				s += this->get_sin(j)*(sin(t*j));
				s += this->get_cos(j)*(cos(t*j));
			}
			return s;
		}

		_T Linfnorm( int N ){
			using std::pow;
			using std::sqrt;
			_T pi = kv::constants< _T >::pi();
			_T h = _T(2)*pi/_T(N);
			_T h0 = _T(0);

			for (int i = 0; i <= N; i++ ){
				h0 = max(h0, abs(this->value( i*h )));
			}

			_T c0 = _T(0);
			int n = this->n;
			for (int j = 1; j <= n; j++){
				c0 += abs(this->get_sin(j))*pow(_T(j),2);
				c0 += abs(this->get_cos(j))*pow(_T(j),2);
			}
			return h0 + c0/_T(2)*pow(h, 2);
		}

		std::ostream& display(std::ostream& os)const
		{
			os << (*this).a0;
			for (int i = 0; i < (*this).n; i++){
				if ((*this).c[i].b != _T(0)){
					if ((*this).c[i].b > _T(0)){
						os << " + " << abs((*this).c[i].b) << "*cos(" << i+1 << "t)";
					}
					else{
						os << " - " << abs((*this).c[i].b) << "*cos(" << i+1 << "t)";
					}
				}
				if ((*this).c[i].a != _T(0)){
					if ((*this).c[i].a > _T(0)){
						os << " + " << abs((*this).c[i].a) << "*sin(" << i+1 << "t)";
					}
					else{
						os << " - " << abs((*this).c[i].a) << "*sin(" << i+1 << "t)";
					}
				}
			}
			os << "\n" << std::endl;
			return os;
		}

		friend std::ostream& operator<<(std::ostream& os, const fourier_series< _T >& A) {
			return A.display(os);
		}
		friend std::ostream& operator<<(std::ostream& os, fourier_series< _T >&& A) {
			return A.display(os);
		}

		void clear(void){
			(*this).n = 0;
			(*this).a0 = 0;
			(*this).c.clear();
		}

    };
}

#endif // VCP_SINCOS_CLASS_HPP