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
    template <typename _T, class _P>
    class doubly_fourier_series{
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
                    for (int p1 = 0; p1 <= r; p1++){
                        for (int p2 = 0; p2 <= r; p2++){
                            if (p1 != 0 && p2 != 1){
                                if ( r == p1 + p2 ){
                                    std::vector< int > tmp;
                                    tmp.resize(2);
                                    
                                        tmp[0] = p1; tmp[1] = p2;
                                        this->plist.push_back(tmp);

                                    if (p2 != 0){
                                        tmp[0] = p1; tmp[1] = -p2;
                                        this->plist.push_back(tmp);
                                    }

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
                                }
                            }
                        }
                    }
                }

                this->plist_sort();
                
                
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
            }
        }

        _T triplet_value(const int& p, const _T& xi, const _T& yi){
            int p1 = plist[p][0];
            int p2 = plist[p][1];
            return p1*this->omega1*xi + p2*this->omega2*yi;
        }

        _T sin_triplet(const int& i, const _T& xi, const _T& yi){
            using std::sin;
            _T trip_value = this->triplet_value( i, xi, yi);
            return sin(trip_value);
        }

        _T cos_triplet(const int& i, const _T& xi, const _T& yi){
            using std::cos;
            _T trip_value = this->triplet_value( i, xi, yi);
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
            this->omega1 = omg2;
            this->T1 = 2*pi/this->omega1;
            this->T2 = 2*pi/this->omega2;
        }

        void setting_uh( vcp::matrix< _T, _P > uh ){
            if (uh.rowsize() > 2*plist.size()+1 ){
                std::cout << "setting_uh: Error...: No match size" << std::endl;
                exit(0);
            }
            this->a0 = uh(0);
            int j = 1;
            for (int i = 0; i < uh.rowsize(); i++){
                this->a_cos[i] = uh(j);
                j++;
                this->b_sin[i] = uh(j);
                j++;
            }
        }

        _T uh_value(const _T& xi, const _T& yi) const {
            _T s = this->a0;
            for (int i = 0; i < plist.size(); i++ ){
                using std::sin;
                using std::cos;
                _T trip_value = this->triplet_value( i, xi, yi);
                s += a_cos[i]*cos(trip_value) + a_cos[i]*sin(trip_value);
            }
            return s;
        }

        _T uh_value(const _T& xi, const _T& yi, const _T& tau) const {
            _T s = this->a0;
            for (int i = 0; i < plist.size(); i++ ){
                using std::sin;
                using std::cos;
                _T trip_value = this->triplet_value( i, xi - tau, yi - tau);
                s += a_cos[i]*cos(trip_value) + a_cos[i]*sin(trip_value);
            }
            return s;
        }

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

            _T pi = kv::constants< _T>::pi();
            _T zero = _T(0);
            int N = 100;

            _T h1 = this->T1/_T(N);
            _T h2 = this->T2/_T(N);

            int_fu_phi.zeros(2*plist.size() + 1);
            for (int i = 0; i < N; i++ ){
                _T xi = i*h1;
                _T xip1 = (i+1)*h1;
                for (int j = 0; j < N; j++){
                    _T yj = j*h2;
                    _T yjp1 = (j+1)*h2;
                    _T tmp = this->uh_value(xi, yj, tau);
                    int_fu_phi(0) += tmp;
                    for (int k = 1; k < 2*plist.size(); k += 2){
                        int_fu_phi(k) += tmp*cos_triplet(k, xi, yj);
                        int_fu_phi(k+1) += tmp*sin_triplet(k+1, xi, yj);
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
                _T xi = i*h1;
                _T xip1 = (i+1)*h1;
                for (int j = 0; j < N; j++){
                    _T yj = j*h2;
                    _T yjp1 = (j+1)*h2;
                    _T tmp = fn(this->uh_value(xi, yj));
                    int_fu_phi(0) += tmp;
                    for (int k = 1; k < 2*plist.size(); k += 2){
                        int_fu_phi(k) += tmp*cos_triplet(k, xi, yj);
                        int_fu_phi(k+1) += tmp*sin_triplet(k+1, xi, yj);
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


    };
}
#endif

