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

#ifndef VCP_FOURIER_COS_BASIS_HPP
#define VCP_FOURIER_COS_BASIS_HPP

#ifdef VCP_NOMP
#ifndef VCP_FOURIER_NOMP
#define VCP_FOURIER_NOMP
#endif
#endif

#include <iostream>
#include <fstream>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/matrix.hpp>

#include <vcp/fourier_series.hpp>

// Basis {1/2, cosx, cos2x, ... , cosnx}
// phi0 = 1/2
// phi1 = cosx
// phi2 = cos2x
// ...

// 0 : a0
// 1 : cos(1t)
// 2 : cos(2t)
// 3 : cos(3t)

// (u,v)_L2 = 1/pi int_0^2pi u v dx 

namespace vcp {
    template <typename _T, class _P>
	class fourier_cos_basis {
protected:
        int m;

public:
        void setting_m(const int& M){
			this->m = M;
		}

        // (dphi, dphi)_L2
        // dphi = cos(ix)' = -i*sin(ix)
        // i^2 int sin(ix)*sin(ix) dx
        vcp::matrix< _T, _P > dphidphi(void){
            vcp::matrix< _T, _P > D;
            D.zeros( this->m + 1,  this->m + 1);

            for (int i = 1; i < this->m+1; i++) {
                D(i,i) = i*i;
            }
            return std::move(D);
        }

        vcp::matrix< _T, _P > phiphi(void){
            vcp::matrix< _T, _P > L;
            L.zeros( this->m + 1,  this->m + 1);

            L(0,0) = 1/_T(2);

            for (int i = 1; i < this->m+1; i++) {
                L(i,i) = _T(1);
            }

            return std::move(L);
        }

        vcp::matrix< _T, _P > fphi(const vcp::fourier_series< _T >& series){
            vcp::matrix< _T, _P > f;
            f.zeros( this->m + 1, 1);

            // (f, 1/2)_L2 = a0/2
            f(0) = series.get_a0();
            for (int i = 1; i < this->m+1; i++) {
                f(i) = series.get_cos(i);
            }

            return std::move(f);
        }

        vcp::matrix< _T, _P > Dfphiphi(const vcp::fourier_series< _T >& series){
            vcp::matrix< _T, _P > Df;
            Df.zeros( this->m + 1, this->m + 1);

            vcp::matrix< _T, _P > tmp = this->fphi(series/_T(2));
            for (int j = 0; j < this->m + 1; j++){
                Df( j , 0) = tmp(j);
            }

            for (int i = 1; i < this->m+1; i++) {
                tmp = this->fphi(series.mul_cos( i ));
                for (int j = 0; j < this->m + 1; j++){
                    Df( j , i) = tmp(j);
                }
            }
 
            return std::move(Df);
        }

        vcp::fourier_series< _T > vec_to_fourier_series(const vcp::matrix< _T, _P >& cof){
            vcp::fourier_series< _T > series;
			
			series.zeros(this->m);
            series.set_a0( cof(0)/_T(2) );
            for (int i = 1; i < this->m+1; i++){
                series.set_cosm( cof(i) , i);
            }
            return series;
        }

        vcp::matrix< _T, _P > graphics(const vcp::fourier_series< _T >& series, const int NN = 100){
            vcp::matrix< _T, _P > graf;
            graf.zeros(NN, 2);
            _T PI = kv::constants<_T>::pi();
            for (int i = 0; i < NN; i++){
                graf(i, 0) = i/_T(NN)*2*PI;
                graf(i, 1) = series.value( graf(i, 0) );
            }
            return graf;
        }
    };
}

#endif // VCP_FOURIER_COS_BASIS_HPP