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

#ifndef VCP_FOURIER_BASIS_HPP
#define VCP_FOURIER_BASIS_HPP

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
    class fourier_basis_list {

public:
		int m;
		std::vector< int > list;
		std::vector< int > omit_list;
		char type;      //'N':NULL 'G':General  'F':Full  'E':Even 'O':Odd

public:
		fourier_basis_list() {
			(*this).m = 0;
			(*this).type = 'N';
		}

		void setting_m(const int& M){
			(*this).m = M;
		}

		void create_full_list( void ){
			(*this).type = 'F';
			for (int i = 0; i < 2*(*this).m+1; i++ ){
				(*this).list.push_back(i);
			}
		}

		void create_odd_list( void ){
			(*this).type = 'O';
			// sum(sin(i*t))
			int j = 0;
			for (int i = 1; i < 2*(*this).m+1; i+=2 ){
				(*this).list.push_back(i);
				j++;
			}
		}

		void create_even_list( void ){
			(*this).type = 'E';
			// a0 + sum(cos(i*t))
			int j = 0;
			for (int i = 0; i < 2*(*this).m+1; i+=2 ){
				(*this).list.push_back(i);
				j++;
			}
		}

		void add_omit_a0( void ){
			(*this).omit_list.push_back( 0 );
		}

		void add_omit_sin( const int& nn ){
			(*this).omit_list.push_back( nn*2 - 1 );
		}

		void add_omit_cos( const int& nn ){
			(*this).omit_list.push_back( nn*2 );
		}
		
		int list_size(void) const{
			return list.size();
		}
    };

	template <typename _T, class _P>
	class fourier_basis : public fourier_basis_list {
protected:		
		vcp::matrix< _T, _P > jmat; // Not omit
		vcp::matrix< _T, _P > jvec;
		vcp::matrix< _T, _P > fxvec;
		bool is_initial;

		void initial_matrix( void ){
			int nn = (*this).list.size();
			int mm = (*this).omit_list.size();
			(*this).jmat.zeros( nn, nn);
			(*this).jvec.zeros( nn, mm);
			(*this).fxvec.zeros( nn, 1);
			(*this).is_initial = false;
		}

		vcp::matrix< _T, _P > to_vector( const vcp::fourier_series<_T>& series, const _T& tt = _T(1)){
			vcp::matrix< _T, _P > V;
			V.zeros((*this).list.size(), 1);

			for (int i = 0; i < (*this).list.size(); i++ ){
				if ((*this).list.at(i) == 0){
					V(i) = tt*series.get_a0();
				}
				else if ((*this).list.at(i)%2 == 1){
					V(i) = series.get_sin((*this).list.at(i)/2+1);
				}
				else if ((*this).list.at(i)%2 == 0){
					V(i) = series.get_cos((*this).list.at(i)/2);
				}
			}
			return V;			
		}

public:
		fourier_basis< _T, _P >() {
			(*this).is_initial = true;
		}

		void clear_data( void ){
			(*this).jmat.clear();
			(*this).jvec.clear();
			(*this).fxvec.clear();
			(*this).is_initial = true;
		}


		void add_fourier_fx( const vcp::fourier_series< _T >& series ){
			if ((*this).is_initial) (*this).initial_matrix();

			vcp::matrix< _T, _P > tmp; 
			tmp = (*this).to_vector( series );
			for (int i = 0; i < (*this).jvec.rowsize(); i ++ ){
				(*this).fxvec(i) += tmp(i);
			}
		}

		template <typename _Tm> typename std::enable_if<std::is_constructible< _T, _Tm >::value, void >::type add_scalar_pt( const _Tm& ssc ){
			if ((*this).is_initial) (*this).initial_matrix();

			_T sc = _T(ssc);
			int nn = (*this).list.size();

			for (int i = 0; i < nn; i++ ){
				int ii = (*this).list.at(i);
				if (ii == 0) {
					(*this).jmat(i, i) += sc/(_T(2));
				}
				else {
					(*this).jmat(i, i) += sc;
				}
			}
		}

		void add_scalar_pt_delay( const _T& sc, const _T& dly ){
			if ((*this).type != 'F') {
				std::cout << "fourier_basis : add_scalar_pt_delay : Only use Full List of fourier_series " << std::endl;
				exit(1);
			}

			if ((*this).is_initial) (*this).initial_matrix();


			for (int i = 0; i < (*this).list.size(); i++ ){
				int ii = (*this).list.at(i);				

				if (ii == 0){
					(*this).jmat(i, i) += sc/(_T(2));
				}
				else if (ii%2 == 1){
					// sin(n(t + delay)) => a*sin(nt) + b*cos(nt)
					// sinterm.a, sinterm.b
					vcp::sincos< _T > sinterm;
					sinterm.set_a(1);
					sinterm.set_b(0);
					sinterm.time_delay( (ii/2+1)*dly );
					sinterm = sc*sinterm;
					(*this).jmat(i, i) += sinterm.a;
					(*this).jmat(i+1, i) += sinterm.b;
				}
				else if (ii%2 == 0){
					// cos(n(t + delay)) => a*sin(nt) + b*cos(nt)
					// costerm.a, costerm.b
					vcp::sincos< _T > costerm;
					costerm.set_a(0);
					costerm.set_b(1);
					costerm.time_delay( (ii/2)* dly );
					costerm = sc*costerm;

					(*this).jmat(i-1, i) += costerm.a;
					(*this).jmat(i, i) += costerm.b;					
				}
			}
		}

		void add_fourier_pt_delay( const vcp::fourier_series< _T >& series, const _T& dly ){
			if ((*this).type != 'F') {
				std::cout << "fourier_basis : add_fourier_pt_delay : Only use Full List of fourier_series " << std::endl;
				exit(1);
			}

			if ((*this).is_initial) (*this).initial_matrix();


			for ( int i = 0; i < (*this).list.size(); i++ ){
				int ii = (*this).list.at(i);
				vcp::fourier_series<_T> tfs;
				vcp::matrix< _T, _P > tmp;
				
				if (ii == 0){
				//	std::cout << "ii/2 = " << ii/2 << std::endl;
					tmp = (*this).to_vector( series/_T(2) );
				}
				else if (ii%2 == 1){
				//	std::cout << "ii/2+1 = " << ii/2 + 1 << std::endl;
					vcp::fourier_series<_T> ftmp;
					ftmp.zeros(ii/2 + 1);
					ftmp.set_sinm( 1, ii/2 + 1);
					ftmp = ftmp.delay(dly);
					tfs = series*ftmp;
					tmp = (*this).to_vector( tfs );
				}
				else if (ii%2 == 0){
				//	std::cout << "ii/2 = " << ii/2 << std::endl;
					vcp::fourier_series<_T> ftmp;
					ftmp.zeros(ii/2);
					ftmp.set_cosm( 1, ii/2);
					ftmp = ftmp.delay(dly);
					tfs = series*ftmp;
					
					tmp = (*this).to_vector( tfs );
				}
				for (int j = 0; j < (*this).list.size(); j++ ){
					(*this).jmat(j, i) += tmp(j);
				}
			}
		}


		void add_fourier_pt( const vcp::fourier_series< _T >& series ){
			if ((*this).is_initial) (*this).initial_matrix();
			int zeroind = -1;

			for ( int i = 0; i < (*this).list.size(); i++ ){
				int ii = (*this).list.at(i);
				vcp::fourier_series<_T> tfs;
				vcp::matrix< _T, _P > tmp;
				
				if (ii == 0){
				//	std::cout << "ii/2 = " << ii/2 << std::endl;
					tmp = (*this).to_vector( series/_T(2) );
				}
				else if (ii%2 == 1){
				//	std::cout << "ii/2+1 = " << ii/2 + 1 << std::endl;
					tfs = series.mul_sin(ii/2 + 1);
					tmp = (*this).to_vector( tfs );
				}
				else if (ii%2 == 0){
				//	std::cout << "ii/2 = " << ii/2 << std::endl;
					tfs = series.mul_cos(ii/2);
					tmp = (*this).to_vector( tfs );
				}
				for (int j = 0; j < (*this).list.size(); j++ ){
					(*this).jmat(j, i) += tmp(j);
				}
			}
		}


		void add_fourier_iota( const vcp::fourier_series< _T >& series, const int& nn = 1 ){
			if ((*this).is_initial) (*this).initial_matrix();
			vcp::matrix< _T, _P > tmp; 
			tmp = (*this).to_vector( series );
			for (int i = 0; i < (*this).jvec.rowsize(); i ++ ){
				(*this).jvec(i, nn-1) += tmp(i);
			}
		}

		template <typename _Tm> typename 
		std::enable_if<std::is_constructible< _T, _Tm >::value, void >::type add_dpt( const _Tm& ssc ){
			if ((*this).type != 'F') {
				std::cout << "fourier_basis : add_dpdx_Jacobi : Only use Full List of fourier_series " << std::endl;
				exit(1);
			}
			if ((*this).is_initial) (*this).initial_matrix();
			_T sc = _T(ssc);

			// 0 : a0
			// 1 : sin(1t) : 1%2 = 1, 1/2+1 = 1
			// 2 : cos(1t) : 2%2 = 0, 2/2   = 1
			// 3 : sin(2t) : 3%2 = 1, 3/2+1 = 2
			// 4 : cos(2t) : 4%2 = 0, 4/2   = 2
			
			for (int i = 0; i < (*this).list.size(); i++){
				int ii = (*this).list.at(i);
				if (ii%2 == 1){
					// (sin(nt))' => n*cos(nt)
					(*this).jmat(i+1, i) += sc*_T(ii/2+1);
				}
				else if ((ii%2 == 0) && (ii != 0)){
					// (cos(nt))' => -n*sin(nt)
					(*this).jmat(i-1, i) -= sc*_T(ii/2);
				}
			}
		}

		void add_dpt(void){
			this->add_dpt( _T(1) );
		}

		template <typename _Tm> typename 
		std::enable_if<std::is_constructible< _T, _Tm >::value, void >::type add_ddpt( const _Tm& ssc ){
			if ((*this).type != 'F') {
				std::cout << "fourier_basis : add_ddpdx_Jacobi : Only use Full List of fourier_series " << std::endl;
				exit(1);
			}
			if ((*this).is_initial) (*this).initial_matrix();
			_T sc = _T(ssc);
			
			for (int i = 0; i < (*this).list.size(); i++){
				int ii = (*this).list.at(i);
				if (ii%2 == 1){
					// (sin(nt))'' => -n*n*sin(nt)
					(*this).jmat(i, i) -= sc*_T(ii/2+1)*_T(ii/2+1);
				}
				else if (ii%2 == 0){
					// (cos(nt))'' => -n*n*cos(nt)
					(*this).jmat(i, i) -= sc*_T(ii/2)*_T(ii/2);
				}
			}
		}

		void add_ddpt(void){
			this->add_ddpt( _T(1) );
		}

		vcp::matrix< _T, _P > output_fx( void )const{
			return (*this).fxvec;
		}

		vcp::matrix< _T, _P > output_Jacobi( void )const{
			vcp::matrix< _T, _P > Jac;
			int nn = (*this).list.size();
			int mm = (*this).omit_list.size();
			bool flags = true;
			int i = 0;
			int jvecnum = 0;

			Jac.zeros( nn, nn );
			for ( int ii : (*this).list ){
				for (int kk : (*this).omit_list){
					if (ii == kk){ 
						flags = false;
						break;
					}
				}
				if (flags) {
					for (int j = 0; j < nn; j++){
						Jac(j, i) = (*this).jmat(j, ii);
					}
					i++;
				}
				flags = true;
			}
			for (int k = 0; k < mm; k++){
				for (int j = 0; j < nn; j++){
					Jac(j, i) = (*this).jvec(j, k);
				}
				i++;
			}

			return Jac;
		}

		vcp::fourier_series< _T > omit_vec_to_fourier_series( const vcp::matrix< _T, _P >& cof ){
			if ( (*this).list.size() - (*this).omit_list.size() != cof.rowsize() ){
				std::cout << "fourier_basis : omit_vec_to_fourier_series : " <<  (*this).list.size() << ", " <<  (*this).omit_list.size() << ", " << cof.rowsize() << std::endl;
				exit(1);
			}

			vcp::fourier_series< _T > series;
			bool flags = true;
			
			series.zeros((*this).m);
			int i = 0;
			for ( int ii : (*this).list ){

				for ( int kk : (*this).omit_list){
					if (ii == kk){ 
						flags = false;
						break;
					}
				}
				if (flags){
					if (ii == 0){
						series.set_a0( cof(i)/_T(2) );
						//series.set_a0( cof(i)/_T(2) );
					}
					else if (ii%2 == 1){
						// sin(nt) : 
						series.set_sinm( cof(i) , ii/2+1);
					}
					else if (ii%2 == 0){
						// (cos(nt))' => -n*sin(nt)
						series.set_cosm( cof(i) , ii/2);
					}
					i++;
				}
				flags = true;
			}
			return series;
		}

		void graphics( const vcp::fourier_series< _T >& x, const _T& t0 = _T(0), const _T& t1 = _T(100) , const _T& tau=_T(1)/_T(100) ){
			int n = x.get_n();
			for (_T i = t0; i <= t1; i += tau ){
				_T s = x.get_a0();
				for (int j = 1; j <= n; j++){
					s += x.get_sin(j)*(sin(j*i));
					s += x.get_cos(j)*(cos(j*i));
				}
				std::cout << i << ", " << s << std::endl;
			}

		}

		void save_graphics( const vcp::fourier_series< _T >& x, const _T& t0 = _T(0), const _T& t1 = _T(100) , const _T& tau=_T(1)/_T(100) ){
			std::string filename; 
	
			std::cout << "Please input filename for graphics" << std::endl;
			std::cin >> filename;
			std::ofstream writing_file;
			writing_file.open(filename, std::ios::out);
			std::cout << "writing " << filename << "..." << std::endl;

			int n = x.get_n();
			for (_T i = t0; i <= t1; i += tau ){
				_T s = x.get_a0();
				for (int j = 1; j <= n; j++){
					s += x.get_sin(j)*(sin(j*i));
					s += x.get_cos(j)*(cos(j*i));
				}
				writing_file << i << ", " << s << std::endl;
			}

		}
	};

}



#endif // VCP_FOURIER_BASIS_HPP