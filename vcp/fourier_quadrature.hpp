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

#ifndef VCP_FOURIER_QUADRATURE_HPP
#define VCP_FOURIER_QUADRATURE_HPP

#ifdef VCP_NOMP
#ifndef VCP_FOURIER_NOMP
#define VCP_FOURIER_NOMP
#endif
#endif

#include <iostream>
#include <fstream>

#include <kv/autodif.hpp>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/matrix.hpp>

#include <functional>

#include <vcp/fourier_series.hpp>

// Basis {1/2, sinx, cosx, sin2x, cos2x, ... , sinnx, cosnx}
// phi0 = 1/2
// phi1 = sinx
// phi2 = cosx
// ...

// 0 : a0
// 1 : sin(1t)
// 2 : cos(1t)
// 3 : sin(2t)
// 4 : cos(2t)
// 5 : sin(3t) 
// 6 : cos(3t) 

namespace vcp {
	template <typename _T, class _P>
    class fourier_quadrature {
protected:
		int bks;
		vcp::matrix< _T, _P > fx;
		std::function< _T(_T) > func;
		std::function< kv::autodif< _T >(kv::autodif< _T >) > diff_func;
		_T pi;

public:
		fourier_quadrature() {
			this->pi = kv::constants< _T >::pi();
			this->bks = 12000;
			fx.zeros(this->bks, 1);
		}

		void set_bks( const int& bkss ){
			this->bks = bkss;
			fx.clear();
			fx.zeros(this->bks, 1);
		}

		void set_func( const std::function< _T(_T) >& ff ){
			this->func = ff;
#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif
			for( int i = 0; i < bks; i++ ){
				_T t=(_T(2)*this->pi/_T(bks))*i;
				fx(i) = ff(t);
			}
		}

		void set_autodif_func( const std::function< kv::autodif< _T >(kv::autodif< _T >) >& ff ){
			this->diff_func = ff;
		}

		/*
			数値積分: zh[m] = 1/pi integral fx * φm dt
			0-2piの積分
		*/
		vcp::matrix< _T, _P > output_vector( const int& order )const{
			vcp::matrix< _T, _P > zh;
			zh.zeros(2*order+1, 1);
			for( int i = 0; i < this->bks; i++ ) zh(0) += this->fx(i);
			zh(0) = zh(0)/_T(this->bks);

			using std::sin;
			using std::cos;


#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif	
			for(int m = 1; m <= order; m++){
				_T sum_sin = _T(0);
				_T sum_cos = _T(0);
				for(int i = 0; i < bks; i++){
					_T t=(_T(2)*pi/_T(bks))*i;
					sum_sin += this->fx(i)*sin(m*t);
					sum_cos += this->fx(i)*cos(m*t);
				}
				zh(2*m-1) = _T(2)*sum_sin/_T(this->bks);
				zh(2*m  ) = _T(2)*sum_cos/_T(this->bks);
			}

			return zh;
		}

		/*
			数値積分: Jmat(i, j) = 1/pi integral f'[x] * φi * φj dt
			0-2piの積分
		*/
		vcp::matrix< _T, _P > output_symmetric_Jacobi( const int& order )const{
			vcp::matrix< _T, _P > Jmat;
			Jmat.zeros(2*order+1, 2*order+1);
			for( int i = 0; i < this->bks; i++ ) Jmat(0, 0) += this->fx(i);
			Jmat(0, 0) = Jmat(0, 0)/_T(this->bks)/_T(2);

			using std::sin;
			using std::cos;

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif	
			for(int m = 1; m <= order; m++){
				_T sum_sin = _T(0);
				_T sum_cos = _T(0);
				for(int i = 0; i < bks; i++){
					_T t=(_T(2)*pi/_T(bks))*i;
					sum_sin += this->fx(i)*sin(m*t);
					sum_cos += this->fx(i)*cos(m*t);
				}
				Jmat( 0, 2*m-1) = sum_sin/_T(this->bks);
				Jmat( 0, 2*m  ) = sum_cos/_T(this->bks);
			}

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif	
			for(int m = 1; m <= order; m++){
				for(int k = m; k <= order; k++){
					_T sum_sinm_sink = _T(0);
					_T sum_sinm_cosk = _T(0);
					_T sum_cosm_sink = _T(0);
					_T sum_cosm_cosk = _T(0);
					for(int i = 0; i < bks; i++){
						_T t=(_T(2)*pi/_T(bks))*i;
						sum_sinm_sink += this->fx(i)*sin(m*t)*sin(k*t);
						sum_sinm_cosk += this->fx(i)*sin(m*t)*cos(k*t);
						sum_cosm_sink += this->fx(i)*cos(m*t)*sin(k*t);
						sum_cosm_cosk += this->fx(i)*cos(m*t)*cos(k*t);
					}
					Jmat(2*m-1, 2*k-1 ) = _T(2)*sum_sinm_sink/_T(this->bks);
					Jmat(2*m-1, 2*k   ) = _T(2)*sum_sinm_cosk/_T(this->bks);
					Jmat(2*m  , 2*k-1 ) = _T(2)*sum_cosm_sink/_T(this->bks);
					Jmat(2*m  , 2*k   ) = _T(2)*sum_cosm_cosk/_T(this->bks);
				}
			}

			for(int i = 0; i < 2*order+1; i++){
				for(int j = i + 1; j < 2*order+1; j++){
					Jmat( j, i ) = Jmat( i, j );
				}
			}

			return Jmat;
		}

		/*
			数値積分: Jmat(i, j) = 1/pi integral f'[x] * φi * φj dt
			0-2piの積分
		*/
		vcp::matrix< _T, _P > output_delay_Jacobi( const int& order, const _T dly )const{
			vcp::matrix< _T, _P > Jmat;
			Jmat.zeros(2*order+1, 2*order+1);
			for( int i = 0; i < this->bks; i++ ) Jmat(0, 0) += this->fx(i);
			Jmat(0, 0) = Jmat(0, 0)/_T(this->bks)/_T(2);

			using std::sin;
			using std::cos;

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif	
			for(int m = 1; m <= order; m++){
				_T sum_sin = _T(0);
				_T sum_cos = _T(0);
				_T sum_sin_dly = _T(0);
				_T sum_cos_dly = _T(0);
				for(int i = 0; i < bks; i++){
					_T t=(_T(2)*pi/_T(bks))*i;
					sum_sin     += this->fx(i)*sin(m*t);
					sum_cos     += this->fx(i)*cos(m*t);					
					sum_sin_dly += this->fx(i)*sin(m*(t + dly));
					sum_cos_dly += this->fx(i)*cos(m*(t + dly));
				}
				Jmat( 2*m-1, 0 ) = sum_sin/_T(this->bks);
				Jmat( 2*m  , 0 ) = sum_cos/_T(this->bks);
				Jmat( 0, 2*m-1 ) = sum_sin_dly/_T(this->bks);
				Jmat( 0, 2*m   ) = sum_cos_dly/_T(this->bks);
			}

#ifdef _OPENMP
#ifndef VCP_FOURIER_NOMP
			#pragma omp parallel for
#endif
#endif	
			for(int m = 1; m <= order; m++){
				for(int k = 1; k <= order; k++){
					_T sum_sinm_sink = _T(0);
					_T sum_sinm_cosk = _T(0);
					_T sum_cosm_sink = _T(0);
					_T sum_cosm_cosk = _T(0);
					for(int i = 0; i < bks; i++){
						_T t=(_T(2)*pi/_T(bks))*i;
						sum_sinm_sink += this->fx(i)*sin(m*t)*sin(k*(t + dly));
						sum_sinm_cosk += this->fx(i)*sin(m*t)*cos(k*(t + dly));
						sum_cosm_sink += this->fx(i)*cos(m*t)*sin(k*(t + dly));
						sum_cosm_cosk += this->fx(i)*cos(m*t)*cos(k*(t + dly));
					}
					Jmat(2*m-1, 2*k-1 ) = _T(2)*sum_sinm_sink/_T(this->bks);
					Jmat(2*m-1, 2*k   ) = _T(2)*sum_sinm_cosk/_T(this->bks);
					Jmat(2*m  , 2*k-1 ) = _T(2)*sum_cosm_sink/_T(this->bks);
					Jmat(2*m  , 2*k   ) = _T(2)*sum_cosm_cosk/_T(this->bks);
				}
			}

			return Jmat;
		}

    };
}



#endif // VCP_FOURIER_QUADRATURE_HPP