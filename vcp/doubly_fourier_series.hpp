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

#ifndef VCP_DOUBLY_FOURIER_SERIES_HPP
#define VCP_DOUBLY_FOURIER_SERIES_HPP

#ifdef VCP_NOMP
#ifndef VCP_DOUBLY_FOURIER_NOMP
#define VCP_DOUBLY_FOURIER_NOMP
#endif
#endif

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <initializer_list>
#include <functional>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/matrix.hpp>

namespace vcp {
    template < typename _T >
    class doubly_fourier_series{
    public:
        int m;
        int elementsize;

		_T a0;
        _T omega1;
        _T omega2;

        _T T1;
        _T T2;

		std::vector< _T > a_cos;
        std::vector< _T > b_sin;

        std::vector< std::vector< int > > plist;

        void plist_sort(){
            sort( this->plist.begin(), this->plist.end(), [](const std::vector< int >& x, const std::vector< int >& y) { return (x[0] < y[0]) || (x[0] == y[0] && (x[1] < y[1])) ;});
            this->plist.shrink_to_fit();
        }

        void setting_list(){
            if ( this->m > 0 ){
                for (int r = 1; r <= this->m; r++){
                    for (int p1 = r; p1 >= 0; p1--){
                        for (int p2 = r; p2 >= 0; p2--){
                            if ( p1 != 0 ){
                                if ( r == p1 + p2 ){
                                    std::vector< int > tmp;
                                    tmp.resize(2);
                                    
                                    tmp[0] = p1; tmp[1] = p2;
                                    this->plist.push_back(tmp);

                                    if (p2 != 0){
                                        tmp[0] = p1; tmp[1] = -p2;
                                        this->plist.push_back(tmp);
                                    }
                                    else {
                                        tmp[0] = p2; tmp[1] = p1;
                                        this->plist.push_back(tmp);                                        
                                    }
/*
                                    if(0){
                                        if (p1 != 0 && p2 != 0 ){
                                            tmp[0] = -p1; tmp[1] = -p2;
                                            this->plist.push_back(tmp);
                                        }

                                        if (p1 != 0){
                                            tmp[0] = -p1; tmp[1] = p2;
                                            this->plist.push_back(tmp);
                                        }
                                    }
                                    */
                                }
                            }
                        }
                    }
                }

                // this->plist_sort();
                                
                // for (const auto& e : plist) std::cout << e[0] << ", " << e[1] << std::endl;                
                this->elementsize = plist.size();
                // std::cout << this->elementsize << std::endl;
            }
        }

        void init_element(){
            if (this->elementsize > 0){
                this->a_cos.resize(this->elementsize);
                this->a_cos.shrink_to_fit();

                this->b_sin.resize(this->elementsize);
                this->b_sin.shrink_to_fit();
                this->a0 = _T(0);
            }
        }

        void AplusB(const doubly_fourier_series< _T >& B){
			if ( this->m != B.m && this->elementsize != B.elementsize){
                std::cout << "Error: AplusB: Not match size or order" << std::endl;
                exit(1);
            }
            this->a0 += B.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
            for (int i = 0; i < this->elementsize; i++){
                this->a_cos[i] += B.a_cos[i];
                this->b_sin[i] += B.b_sin[i];
            }
        }

        void AminusB(const doubly_fourier_series< _T >& B){
			if ( this->m != B.m && this->elementsize != B.elementsize){
                std::cout << "Error: AminusB: Not match size or order" << std::endl;
                exit(1);
            }
            this->a0 -= B.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
            for (int i = 0; i < this->elementsize; i++){
                this->a_cos[i] -= B.a_cos[i];
                this->b_sin[i] -= B.b_sin[i];
            }
        }

        void minusAplusB(const doubly_fourier_series< _T >& B){
            if ( this->m != B.m && this->elementsize != B.elementsize){
                std::cout << "Error: minusAplusB: Not match size or order" << std::endl;
                exit(1);
            }
            this->a0 = - this->a0 + B.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
            for (int i = 0; i < this->elementsize; i++){
                this->a_cos[i] = -this->a_cos[i] + B.a_cos[i];
                this->b_sin[i] = -this->b_sin[i] + B.b_sin[i];
            }
        }


        void scalar_mul(const _T& c){
            this->a0 *= c;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
            for (int i = 0; i < this->elementsize; i++){
                this->a_cos[i] *= c;
                this->b_sin[i] *= c;
            }
        }

        _T triplet_value(const int& p, const _T& u1, const _T& u2){
            int p1 = this->plist[p][0];
            int p2 = this->plist[p][1];
            return p1*this->omega1*u1 + p2*this->omega2*u2;
        }

        _T sin_triplet(const int& i, const _T& u1, const _T& u2){
            using std::sin;
            _T trip_value = this->triplet_value( i, u1, u2);
            return sin(trip_value);
        }

        _T cos_triplet(const int& i, const _T& u1, const _T& u2){
            using std::cos;
            _T trip_value = this->triplet_value( i, u1, u2);
            return cos(trip_value);
        }

    public:
        doubly_fourier_series(){
            this->m = -1;
            this->elementsize = -1;
        }

        void setting_order( const int& mm ){
            this->m = mm;
            this->setting_list();
            this->init_element();
        }

        void setting_omega( const _T& omg1, const _T& omg2 ){
            _T pi = kv::constants< _T>::pi();
            this->omega1 = omg1;
            this->omega2 = omg2;
            this->T1 = 2*pi/this->omega1;
            this->T2 = 2*pi/this->omega2;
        }

        std::vector< std::vector< int > > output_list( void ) const {
            return this->plist;
        }


        doubly_fourier_series< _T > diff_u1( void ) const {
            doubly_fourier_series< _T > tmp;
            tmp.setting_order( this->m );
            tmp.setting_omega( this->omega1, this->omega2 );
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif            
            for (int i = 0; i < this->elementsize; i++ ){
                int p1 = this->plist[i][0];
                tmp.a_cos[ i ] =  p1*this->omega1*this->b_sin[i];
                tmp.b_sin[ i ] = -p1*this->omega1*this->a_cos[i];
            }
            return tmp;
        }

        doubly_fourier_series< _T > diff_u2( void ) const {
            doubly_fourier_series< _T > tmp;
            tmp.setting_order( this->m );
            tmp.setting_omega( this->omega1, this->omega2 );
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif            
            for (int i = 0; i < this->elementsize; i++ ){
                int p2 = this->plist[i][1];
                tmp.a_cos[ i ] =  p2*this->omega2*this->b_sin[i];
                tmp.b_sin[ i ] = -p2*this->omega2*this->a_cos[i];
            }
            return tmp;
        }

        // pu = p1*omega1*u1 + p2*omega2*u2
        // d  = p1*omega1*tau + p2*omega2*tau
		// d > 0
		// b*sin( p1*omega1*(u1+tau) + p2*omega2*(u2+tau) ) 
        //       = b*sin( pu + d )
        //       = b*( sin(pu)cos(d) + cos(pu)sin(d) )
        //       = b*cos(d)*sin(pu) + b*sin(d)*cos(pu)
        // a*cos( p1*omega1*(u1+tau) + p2*omega2*(u2+tau) )
        //       =  a*cos( pu + d )
        //       =  a*( cos(pu)cos(d) - sin(pu)sin(d) )
        //       = -a*sin(d)*sin(pu) + a*cos(d)*cos(pu)
        //
        // d < 0
		// b*sin( p1*omega1*(u1+tau) + p2*omega2*(u2+tau) ) 
        //       = b*sin( pu - abs(d) )
        //       = b*( sin(pu)cos(abs(d)) - cos(pu)sin(abs(d)) )
        //       = b*cos(abs(d))*sin(pu) - b*sin(abs(d))*cos(pu)
        // a*cos( p1*omega1*(u1+tau) + p2*omega2*(u2+tau) )
        //       =  a*cos( pu - abs(d) )
        //       =  a*( cos(pu)cos(abs(d)) + sin(pu)sin(abs(d)) )
        //       =  a*sin(abs(d))*sin(pu) + a*cos(abs(d))*cos(pu)
        template <typename _Tm>
        typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type delay( const _Tm& taum ) const {
            using std::sin;
            using std::cos;
            using std::abs;            
            _T tau = _T(taum);
            doubly_fourier_series< _T > tmp;
            tmp.setting_order( this->m );
            tmp.setting_omega( this->omega1, this->omega2 );
            tmp.a0 = this->a0;

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
            for (int i = 0; i < this->elementsize; i++){
                int p1 = this->plist[i][0];
                int p2 = this->plist[i][1];
                _T d  = p1*this->omega1*tau + p2*this->omega2*tau;
                _T abs_d = abs(d);
                if (d > 0){
                    tmp.a_cos[i] = this->a_cos[i]*sin(d) + this->b_sin[i]*cos(d);
                    tmp.b_sin[i] = this->a_cos[i]*cos(d) - this->b_sin[i]*sin(d);
                }
                else {
                    tmp.a_cos[i] = -this->a_cos[i]*sin(abs_d) + this->b_sin[i]*cos(abs_d);
                    tmp.b_sin[i] =  this->a_cos[i]*cos(abs_d) - this->b_sin[i]*sin(abs_d);
                }
            }
            return tmp;        
        }

		template <class _Pm>
        void setting_uh( const vcp::matrix< _T, _Pm >& uh ){
            if (uh.rowsize() > 2*plist.size()+1 ){
                std::cout << "setting_uh: Error...: No match size" << std::endl;
                exit(0);
            }
            this->a0 = uh(0);
            int j = 1;
            for (int i = 0; i < plist.size(); i++){
                this->a_cos[i] = uh(j);
                j++;
                this->b_sin[i] = uh(j);
                j++;
            }
        }
        template <class _Pm>
        vcp::matrix< _T, _Pm > output_matrix( void ) const {
            vcp::matrix< _T, _Pm > tmp;
            tmp.zeros( 2*this->plist.size()+1, 1 );
            tmp(0) = this->a0;
            int j = 1;
            for (int i = 0; i < this->plist.size(); i++){
                tmp(j) = this->a_cos[i];
                j++;
                tmp(j) = this->b_sin[i];
                j++;
            }
            return tmp;
        }

        _T uh_value(const _T& u1, const _T& u2) const {
            _T s = this->a0;
            for (int i = 0; i < this->plist.size(); i++ ){
                using std::sin;
                using std::cos;
                _T trip_value = this->triplet_value( i, u1, u2);
                s += a_cos[i]*cos(trip_value) + a_cos[i]*sin(trip_value);
            }
            return s;
        }

        _T uh_delay_value(const _T& u1, const _T& u2, const _T& tau) const {
            _T s = this->a0;
            for (int i = 0; i < this->plist.size(); i++ ){
                using std::sin;
                using std::cos;
                _T trip_value = this->triplet_value( i, u1 - tau, u2 - tau);
                s += a_cos[i]*cos(trip_value) + a_cos[i]*sin(trip_value);
            }
            return s;
        }

/*
        // int_0^2pi (u phi) dx dy
        vcp::matrix< _T, _P > integrate_fuphi(void) const {
            vcp::matrix< _T, _P > int_fu_phi;
            int_fu_phi.zeros(2*plist.size() + 1);
            int_fu_phi(0) = this->a0;
            int j = 1;
            for (int i = 0; i < plist.size(); i++) {
                int_fu_phi(j) = this->a_cos[i];
                j++;
                int_fu_phi(j) = this->b_sin[i];
                j++;
            }
            return int_fu_phi;
        }

        // int_0^2pi (u(tau) phi) dx dy
        vcp::matrix< _T, _P > integrate_fuphi(const _T& tau) const {
            vcp::matrix< _T, _P > int_fu_phi;

            _T pi = kv::constants< _T >::pi();
            _T zero = _T(0);
            int N = 100;

            _T h1 = this->T1/_T(N);
            _T h2 = this->T2/_T(N);

            int_fu_phi.zeros(2*plist.size() + 1);
            for (int i = 0; i < N; i++ ){
                _T u1 = i*h1;
                _T u1p1 = (i+1)*h1;
                for (int j = 0; j < N; j++){
                    _T yj = j*h2;
                    _T yjp1 = (j+1)*h2;
                    _T tmp = this->uh_value(u1, yj, tau);
                    int_fu_phi(0) += tmp;
                    for (int k = 1; k < 2*plist.size(); k += 2){
                        int_fu_phi(k) += tmp*cos_triplet(k, u1, yj);
                        int_fu_phi(k+1) += tmp*sin_triplet(k+1, u1, yj);
                    }
                }
            }

            int_fu_phi = h1*h2*int_fu_phi;
            return int_fu_phi;
        }

        // int_0^2pi (f(u) phi) dx dy
        vcp::matrix< _T, _P > integrate_fuphi(std::function<_T(_T)> fn) const {
            vcp::matrix< _T, _P > int_fu_phi;

            _T pi = kv::constants< _T>::pi();
            _T zero = _T(0);
            int N = 100;

            _T h1 = this->T1/_T(N);
            _T h2 = this->T2/_T(N);

            int_fu_phi.zeros(2*plist.size() + 1);
            for (int i = 0; i < N; i++ ){
                _T u1 = i*h1;
                _T u1p1 = (i+1)*h1;
                for (int j = 0; j < N; j++){
                    _T yj = j*h2;
                    _T yjp1 = (j+1)*h2;
                    _T tmp = fn(this->uh_value(u1, yj));
                    int_fu_phi(0) += tmp;
                    for (int k = 1; k < 2*plist.size(); k += 2){
                        int_fu_phi(k) += tmp*cos_triplet(k, u1, yj);
                        int_fu_phi(k+1) += tmp*sin_triplet(k+1, u1, yj);
                    }
                }
            }

            int_fu_phi = h1*h2*int_fu_phi;
            return int_fu_phi;
        }

        vcp::matrix< _T, _P > integrate_du1( void ) const {
            vcp::matrix< _T, _P > tmp;
            tmp.zeros( 2*(this->plist.size())+1 , 1);
            tmp(0) = _T(0);
            int j = 1;
            for (int i = 0; i < this->plist.size(); i++ ) {
                int p1 = plist[i][0];            
                tmp(j)   = p1*this->omega1*b_sin[i];
                j++;
                tmp(j+1) = -p1*this->omega1*a_cos[i];
                j++;
            }
            return tmp;
        }

        vcp::matrix< _T, _P > integrate_du2( void ) const {
            vcp::matrix< _T, _P > tmp;
            tmp.zeros( 2*(this->plist.size())+1 , 1);
            tmp(0) = _T(0);
            int j = 1;
            for (int i = 0; i < this->plist.size(); i++ ) {
                int p2 = plist[i][1];            
                tmp(j)   = p2*this->omega2*b_sin[i];
                j++;
                tmp(j+1) = -p2*this->omega2*a_cos[i];
                j++;
            }
            return tmp;
        }
*/

        friend doubly_fourier_series< _T > 
        operator+(const doubly_fourier_series< _T >& A) {
			return A;
		}
        
		friend doubly_fourier_series< _T > 
        operator+(doubly_fourier_series< _T >&& A) {
			return std::move(A);
		}

        friend doubly_fourier_series< _T > 
        operator+(const doubly_fourier_series< _T >& A, const doubly_fourier_series< _T >& B) {
			doubly_fourier_series< _T > C;
			C = A;
			C.AplusB(B);
			return std::move(C);
		}

		friend doubly_fourier_series< _T > 
        operator+(doubly_fourier_series< _T >&& A, const doubly_fourier_series< _T >& B) {
			A.AplusB(B);
			return std::move(A);
		}

		friend doubly_fourier_series< _T > 
        operator+(const doubly_fourier_series< _T >& A, doubly_fourier_series< _T >&& B) {
			B.AplusB(A);
			return std::move(B);
		}

		friend doubly_fourier_series< _T > 
        operator+(doubly_fourier_series< _T >&& A, doubly_fourier_series< _T >&& B) {
			A.AplusB(B);
			return std::move(A);
		}

		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator+( const doubly_fourier_series<_T>& A, const _Tm& a ){
			doubly_fourier_series< _T > C;
			_T Ta = _T(a);
			C = A;
			C.a0 += Ta;
			return C;
		}

		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator+( doubly_fourier_series<_T>&& A, const _Tm& a ){
			_T Ta = _T(a);
			A.a0 += Ta;
			return std::move(A);
		}


		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator+( const _Tm& a, const doubly_fourier_series<_T>& A ){
			doubly_fourier_series<_T> C;
			_T Ta = _T(a);
			C = A;
			C.a0 += Ta;
			return C;
		}

		template <typename _Tm>
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator+( const _Tm& a, doubly_fourier_series<_T>&& A ){
			_T Ta = _T(a);
			A.a0 += Ta;
			return std::move(A);
		}


        friend doubly_fourier_series< _T >& operator+=( doubly_fourier_series< _T >& A, const doubly_fourier_series< _T >& B) {
			A.AplusB(B);
			return A;
		}

        friend doubly_fourier_series< _T > operator-(const doubly_fourier_series< _T >& A) {
            doubly_fourier_series< _T > B;
			B.setting_order( A.m );
            B.setting_omega( A.omega1, A.omega2 );
			B.a0 = -A.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < A.elementsize; i++){
				B.a_cos[i] = -A.a_cos[i];
                B.b_sin[i] = -A.b_sin[i];
			}
			return B;
		}

        friend doubly_fourier_series< _T > operator-(doubly_fourier_series< _T >&& A) {
            A.a0 = -A.a0;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < A.elementsize; i++){
				A.a_cos[i] = -A.a_cos[i];
                A.b_sin[i] = -A.b_sin[i];
			}
			return std::move(A);
		}

        friend doubly_fourier_series< _T > operator-(const doubly_fourier_series< _T >& A, const doubly_fourier_series< _T >& B) {
			doubly_fourier_series< _T > C;
			C = A;
			C.AminusB(B);
			return std::move(C);
		}

		friend doubly_fourier_series< _T > operator-(doubly_fourier_series< _T >&& A, const doubly_fourier_series< _T >& B) {
			A.AminusB(B);
			return std::move(A);
		}

		friend doubly_fourier_series< _T > operator-(const doubly_fourier_series< _T >& A, doubly_fourier_series< _T >&& B) {
			B.minusAplusB(A);
			return std::move(B);
		}

		friend doubly_fourier_series< _T > operator-(doubly_fourier_series< _T >&& A, doubly_fourier_series< _T >&& B) {
			A.AminusB(B);
			return std::move(A);
		}

		
		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator-( const doubly_fourier_series<_T>& A, const _Tm& a ){
			doubly_fourier_series<_T> C;
			_T Ta = _T(a);
			C = A;
			C.a0 -= Ta;
			return C;
		}

		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator-( doubly_fourier_series<_T>&& A, const _Tm& a ){
			_T Ta = _T(a);
			A.a0 -= Ta;
			return std::move(A);
		}


		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator-( const _Tm& a, const doubly_fourier_series<_T>& A ){
			doubly_fourier_series<_T> C;
			_T Ta = _T(a);
			C = -A;
			C.a0 += Ta;
			return C;
		}

		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator-( const _Tm& a, doubly_fourier_series<_T>&& A ){
			_T Ta = _T(a);
			A.a0 -= Ta;
			return std::move(-A);
		}

        template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator*(const doubly_fourier_series< _T >& B, const _Tm& a) {
			_T Ta = _T(a);
			doubly_fourier_series< _T > C;
			C = B;
			C.scalar_mul(Ta);
			return std::move(C);
		}

		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator*(doubly_fourier_series< _T >&& B, const _Tm& a) {
			_T Ta = _T(a);
            B.scalar_mul(Ta);
			return std::move(B);
		}

        template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator*(const _Tm& a, const doubly_fourier_series< _T >& B) {
			_T Ta = _T(a);
			doubly_fourier_series< _T > C;
			C = B;
			C.scalar_mul(Ta);
			return std::move(C);
		}

		template <typename _Tm> 
        friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, doubly_fourier_series< _T > >::type operator*(const _Tm& a, doubly_fourier_series< _T >&& B) {
			_T Ta = _T(a);
            B.scalar_mul(Ta);
			return std::move(B);
		}

        std::ostream& display(std::ostream& os) const {
            os << "ω1 = " << this->omega1 << ", ω2 = " << this->omega2 << "\n";
			os << this->a0;
			for (int i = 0; i < this->elementsize; i++){
                int p1 = this->plist[i][0];
                int p2 = this->plist[i][1];                
				if ( this->a_cos[i] != _T(0)){
					if (this->a_cos[i] > _T(0)){
						os << " + " << abs(this->a_cos[i]) << "*cos(" << p1 << "ω1u1 + " << p2 << "ω2u2)";
					}
					else{
						os << " - " << abs(this->a_cos[i]) << "*cos(" << p1 << "ω1u1 + " << p2 << "ω2u2)";
					}
				}
				if (this->b_sin[i] != _T(0)){
					if (this->b_sin[i] > _T(0)){
						os << " + " << abs(this->b_sin[i]) << "*sin(" << p1 << "ω1u1 + " << p2 << "ω2u2)\n";
					}
					else{
						os << " - " << abs(this->b_sin[i]) << "*sin(" << p1 << "ω1u1 + " << p2 << "ω2u2)\n";
					}
				}
			}
			os << std::endl;
			return os;
		}

        friend std::ostream& operator<<(std::ostream& os, const doubly_fourier_series< _T >& A) {
			return A.display(os);
		}

		friend std::ostream& operator<<(std::ostream& os, doubly_fourier_series< _T >&& A) {
			return A.display(os);
		}

    };

// {1/2, cos{}, }
    template <typename _T, class _P>
	class doubly_fourier_basis{
		vcp::matrix< _T, _P > jvec;
		vcp::matrix< _T, _P > jmat;
        vcp::matrix< _T, _P > jvec_eota;
    protected:		



    };
}
#endif
