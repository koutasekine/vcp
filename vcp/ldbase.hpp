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

#ifndef VCP_LDBASE_HPP
#define VCP_LDBASE_HPP

#ifndef INTERVAL_HPP
#error Please include interval.hpp
#endif

#include <iostream>

#include <limits>
#include <stack>
#include <vector>

#include <kv/psa.hpp>
#include <kv/gamma.hpp>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/vcp_psa_assist.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/error.hpp>
#include <vcp/matrix.hpp>

#include <vcp/legendre_integral.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef VCP_DEBUG
#ifndef VCP_LEGENDRE_DEBUG
	#define VCP_LEGENDRE_DEBUG
#endif
#endif

#ifdef VCP_NOMP
#ifndef VCP_LEGENDRE_NOMP
	#define VCP_LEGENDRE_NOMP
#endif
#endif

#ifdef VCP_USE_DEATH_FUNC
#ifndef VCP_LEGENDRE_USE_DEATH_FUNC
	#define VCP_LEGENDRE_USE_DEATH_FUNC
#endif
#endif

namespace vcp {
	namespace LegendreBaseFunctions {
		inline int checked_int_mul(const int a, const int b, const char* context) {
			if (a < 0 || b < 0) {
				vcp::throw_error<vcp::invalid_argument>(context, ": negative integer size");
			}
			if (a != 0 && b > std::numeric_limits<int>::max() / a) {
				vcp::throw_error<vcp::dimension_error>(context, ": integer size overflow");
			}
			return a * b;
		}

		inline int checked_int_pow(const int base, const int exp, const char* context) {
			if (base < 0 || exp < 0) {
				vcp::throw_error<vcp::invalid_argument>(context, ": negative integer power argument");
			}
			int result = 1;
			for (int i = 0; i < exp; i++) {
				result = checked_int_mul(result, base, context);
			}
			return result;
		}

		inline std::vector<int> checked_int_power_table(const int base, const int max_exp, const char* context) {
			if (base < 0 || max_exp < 0) {
				vcp::throw_error<vcp::invalid_argument>(context, ": negative integer power argument");
			}
			std::vector<int> powers(max_exp + 1);
			powers[0] = 1;
			for (int i = 1; i <= max_exp; i++) {
				powers[i] = checked_int_mul(powers[i - 1], base, context);
			}
			return powers;
		}

		template <typename _T> void psaTodpsa(const std::vector< kv::psa< _T > >& inp, std::vector< kv::psa< _T > >& outp) {
			int n = inp.size();

			outp.resize(n);
			for (int i = 0; i < n; i++) {
				outp[i].v.resize(n);
				for (int j = 0; j < n; j++) {
					outp[i].v[j] = _T(0);
				}
			}
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n - 1; j++) {
					outp[i].v[j] = (j + 1)*inp[i].v[j + 1];
				}
			}
		}
		template <typename _T> void NmalLegendrePol(const int nn, std::vector< kv::psa< _T > >& outp) {
			int i, j, n;
			kv::psa< _T > x;
			n = nn + 1;

			outp.resize(n);
			for (i = 0; i < n; i++) {
				outp[i].v.resize(n);
				for (j = 0; j < n; j++) {
					outp[i].v[j] = _T(0);
				}
			}

			x.v.resize(n);
			for (i = 0; i < n; i++) {
				x.v[i] = _T(0);
			}
			x.v[1] = _T(1);
			outp[0].v[0] = _T(1);
			outp[1].v[0] = _T(-1);
			outp[1].v[1] = _T(2);
			for (i = 2; i < n; i++) {
				outp[i] = ((2 * i - 1)*(2 * x - 1)*outp[i - 1] - (i - 1)*outp[i - 2]) / i;
			}
		}
		template <typename _T> void makeLegendreBase(const int nn, std::vector< kv::psa< _T > >& phi) {
			std::vector< kv::psa< _T > > Npol;
			std::vector< kv::psa< _T > > DNpol;
			int i, j, n;
			kv::psa< _T > x;
			n = nn + 1;

			NmalLegendrePol(nn, Npol);
			psaTodpsa(Npol, DNpol);

			phi.resize(n);
			for (i = 0; i<n; i++) {
				phi[i].v.resize(n);
				for (j = 0; j < n; j++) {
					phi[i].v[j] = _T(0);
				}
			}

			x.v.resize(n);
			for (i = 0; i < n; i++) {
				x.v[i] = _T(0);
			}
			x.v[1] = _T(1);

			for (i = 2; i<n; i++) {
				phi[i] = (1 - x)*x*DNpol[i - 1] / (i*(i - 1));
			}
		}
		template <typename _T> void LegendrePointFunc(const std::vector< kv::psa< _T > >& phi, const std::vector< _T >& Point, std::vector< std::vector< _T > >& out) {

			out.resize(phi.size());
			for (int i = 0; i < phi.size(); i++) {
				out[i].resize(Point.size());
			}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < phi.size(); i++) {
				for (int j = 0; j < Point.size(); j++) {
					psa_value_Horner(phi[i], Point[j], out[i][j]);
				}
			}

		}

	}

	template <typename _T> class Legendre_base{
	protected:
		std::vector< kv::psa< _T > > phi;
		std::vector< kv::psa< _T > > dphi;
		std::vector< kv::psa< _T > > ddphi;
		// Order of Legendre Base
		int m;
		
	public:
		Legendre_base< _T >() {
		}
		void set(const int n) {
			m = n;
			int mm = m + 2;
			LegendreBaseFunctions::makeLegendreBase(mm, phi);
			LegendreBaseFunctions::psaTodpsa(phi, dphi);
			LegendreBaseFunctions::psaTodpsa(dphi, ddphi);
		}
	};
	
	template <typename _T, typename _TM, class _PM = mats< _TM >> class Legendre_Bases_Generator : protected interval_ld_weightpoint< _T > {
	protected:
//	public:
		int Order_of_Base;   // Order of Legendre Base
		int elementsize;     // integer power of this->phi.size() and this->dimension
		int variablesize;    // Number of Variables
		int variablesize_cx; // Number of Variables for cx
		int dimension;       // Dimension of Domein
		int mode;			 // 1: verification mode (Inverse operator), 2: verification mode (Residual), k > 2: approximation mode k*10
		int p;				 // -Delta u = u^p
		int Order_uh;		 // Order_uh
		int Order_cx;		 // Order_cx
		bool flag_order_uh;
		bool flag_order_cx;

		std::vector< kv::psa< _T > > phi;
		std::vector< kv::psa< _T > > dphi;
		std::vector< kv::psa< _T > > ddphi;

		matrix< int > list;
		matrix< int > Pointlist;
		matrix< _TM, _PM > Legendre_with_GLpoint;
		matrix< _TM, _PM > DifDifLegendre_with_GLpoint;
		matrix< _TM, _PM > uh;
		matrix< _TM, _PM > cx;
		matrix< _TM, _PM > phi_point;
		matrix< _TM, _PM > uhphi_point;
		matrix< _TM, _PM > cxphi_point;		
		matrix< _TM, _PM > DDuhphi_point;
		matrix< _TM, _PM > weight_point;

		void list_set() {
			int mmax = this->phi.size();
			int dim = this->dimension;
			int di;
			const std::vector<int> powers = LegendreBaseFunctions::checked_int_power_table(mmax, dim, "list_set");

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: List of bases------------------------------------------" << std::endl;
			std::cout << ">> Size of list is (" << powers[dim] << "," << dim << ") \n" << std::endl;
#endif
			this->list.zeros(powers[dim], dim);
			for (int d = 0; d < dim; d++) {
				di = dim - d;
				const int lower_block = powers[di - 1];
				const int upper_block = powers[di];
				for (int j = 0; j < powers[d]; j++) {
					for (int i = 0; i < mmax; i++) {
						for (int k = 0; k < lower_block; k++) {
							this->list(i*lower_block + j * upper_block + k, d) = i;
						}
					}
				}
			}
		}
		void even_list_set() {
			int mmax = this->phi.size();
			int mnumber = (mmax + 1) / 2;
			int dim = this->dimension;
			int di;
			const std::vector<int> powers = LegendreBaseFunctions::checked_int_power_table(mnumber, dim, "even_list_set");

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Even List of bases-------------------------------------" << std::endl;
			std::cout << ">> Size of list is (" << powers[dim] << "," << dim << ") \n" << std::endl;
#endif
			try {
				this->list.zeros(powers[dim], dim);
			}
			catch (const std::length_error& le){
				std::cout << "even_list_set:" << std::endl;
				std::cout << "Error: Matrix size: (" << powers[dim] << ", " << dim << " )" << std::endl;
				std::cerr << "Length error: " << le.what() << '\n';
			}
			int ii = 0;
			for (int d = 0; d < dim; d++) {
				di = dim - d;
				const int lower_block = powers[di - 1];
				const int upper_block = powers[di];
				for (int j = 0; j < powers[d]; j++) {
					for (int i = 0; i < mmax; i = i + 2) {
						for (int k = 0; k < lower_block; k++) {
							this->list(ii*lower_block + j * upper_block + k, d) = i;
						}
						ii++;
					}
					ii = 0;
				}
			}
		}
		void Pointlist_set() {
			int mmax = this->Legendre_with_GLpoint.columnsize();
			int dim = this->dimension;
			int di;
			const std::vector<int> powers = LegendreBaseFunctions::checked_int_power_table(mmax, dim, "Pointlist_set");
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: List of root of Legendre polynomial--------------------" << std::endl;
			std::cout << ">> Size of list is (" << powers[dim] << "," << dim << ") \n" << std::endl;
#endif
			this->Pointlist.zeros(powers[dim], dim);
			for (int d = 0; d < dim; d++) {
				di = dim - d;
				const int lower_block = powers[di - 1];
				const int upper_block = powers[di];
				for (int j = 0; j < powers[d]; j++) {
					for (int i = 0; i < mmax; i++) {
						for (int k = 0; k < lower_block; k++) {
							this->Pointlist(i*lower_block + j * upper_block + k, d) = i;
						}
					}
				}
			}
		}
		void weight_point_set() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Matrix for Weight_Point--------------------------------" << std::endl;
			std::cout << ">> Matrix size of Legendre Base is (" << 1 << "," << this->Pointlist.rowsize() << ") \n" << std::endl;
#endif
			this->weight_point.ones(1, this->Pointlist.rowsize());
			for (int i = 0; i < this->Pointlist.rowsize(); i++) {
				_T a = this->Weight[this->Pointlist(i, 0)];
				_TM b;
				convert(a, b);
				this->weight_point(0, i) = b;
				for (int j = 1; j < this->Pointlist.columnsize(); j++) {
					a = this->Weight[this->Pointlist(i, j)];
					convert(a, b);
					this->weight_point(0, i) *= b;
				}
			}
		}
		void phi_point_set() {
			int column_Pointlist = this->Pointlist.rowsize();
			int row_Point = this->list.rowsize();
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Matrix for Legendre base | Row:List | Column:PoinList--" << std::endl;
			std::cout << ">> Matrix size of Legendre Base is (" << row_Point << "," << column_Pointlist << ")" << std::endl;
#endif
			this->phi_point.ones(row_Point, column_Pointlist);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < row_Point; i++) {
				for (int j = 0; j < column_Pointlist; j++) {
					for (int d = 0; d < this->list.columnsize(); d++) {
						this->phi_point(i, j) *= this->Legendre_with_GLpoint(this->list(i, d), this->Pointlist(j, d));
					}
				}
			}
		}
		void uhphi_point_set() {
			int column_Pointlist = this->Pointlist.rowsize();
			int row_Point = this->list.rowsize();
			this->Order_uh = this->Order_of_Base;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: uh_i * phi_i with PointList----------------------------" << std::endl;
			std::cout << ">> uh size is (" << this->elementsize << "," << this->variablesize << ")" << std::endl;
			std::cout << ">> uh_i * phi_i size is (" << this->variablesize << "," << column_Pointlist << ")\n" << std::endl;
#endif
			this->uhphi_point.zeros(this->variablesize, column_Pointlist);
			this->DDuhphi_point.zeros(this->variablesize, column_Pointlist);
			if ( ( this->mode == 2 ) || ( this->mode == 0 ) ) {
				_TM Phi_Time;
				for (int v = 0; v < this->variablesize; v++) {
					for (int j = 0; j < this->Pointlist.rowsize(); j++) {
						for (int i = 0; i < this->list.rowsize(); i++) {
							Phi_Time = this->Legendre_with_GLpoint(this->list(i, 0), this->Pointlist(j, 0));
							for (int d = 1; d < this->list.columnsize(); d++) {
								Phi_Time *= this->Legendre_with_GLpoint(this->list(i, d), this->Pointlist(j, d));
							}
							this->uhphi_point(v, j) += uh(i, v) * Phi_Time;
						}
					}
				}

				_TM DDPhi_Time;
				for (int dd = 0; dd < this->list.columnsize(); dd++) {
					for (int v = 0; v < this->variablesize; v++) {
						for (int j = 0; j < this->Pointlist.rowsize(); j++) {
							for (int i = 0; i < this->list.rowsize(); i++) {
								DDPhi_Time = _TM(1);
								for (int d = 0; d < this->list.columnsize(); d++) {
									if (d == dd) {
										DDPhi_Time *= this->DifDifLegendre_with_GLpoint(this->list(i, d), this->Pointlist(j, d));
									}
									else {
										DDPhi_Time *= this->Legendre_with_GLpoint(this->list(i, d), this->Pointlist(j, d));
									}
								}
								this->DDuhphi_point(v, j) += uh(i, v) * DDPhi_Time;
							}
						}
					}
				}


			}
			else {
				for (int v = 0; v < this->variablesize; v++) {
					for (int i = 0; i < column_Pointlist; i++) {
						this->uhphi_point(v, i) = this->uh(0, v) * this->phi_point(0, i);
						for (int j = 1; j < row_Point; j++) {
							this->uhphi_point(v, i) += this->uh(j, v) * this->phi_point(j, i);
						}
					}
				}
			}
		}
		void uhphi_point_set(const matrix< int >& list_uh, int& order_uh) {
			this->Order_uh = order_uh;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base------------------------------------------" << std::endl;
			std::cout << ">> Order of Legendre base for uh is " << this->Order_uh << std::endl;
			std::cout << ">> Order of polynomial is " << this->Order_uh + 2 << " (= Order_uh + 2)\n" << std::endl;
#endif
			std::vector< kv::psa< _T > > phi_uh;
			LegendreBaseFunctions::makeLegendreBase(this->Order_uh + 2, phi_uh);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base <== Point of Gauss Legendre--------------" << std::endl;
			std::cout << ">> Size of Legendre_with_GLpoint is (" << this->Order_uh << "," << this->Point.size() << ")" << std::endl;
			std::cout << ">> Size of DifDifLegendre_with_GLpoint is (" << this->Order_uh << "," << this->Point.size() << ")\n" << std::endl;
#endif
			std::vector< std::vector< _T > > LBP2;
			LegendreBaseFunctions::LegendrePointFunc(phi_uh, this->Point, LBP2);

			matrix< _TM, _PM > Legendre_uh_with_GLpoint;
			Legendre_uh_with_GLpoint.zeros(this->Order_uh, this->Point.size());

			for (int i = 0; i < this->Order_uh; i++) {
				for (int j = 0; j < this->Point.size(); j++) {
					convert(LBP2[i + 2][j], Legendre_uh_with_GLpoint(i, j));
//					Legendre_uh_with_GLpoint(i, j) = LBP2[i + 2][j];
				}
			}
			for (int j = 0; j < phi_uh.size() - 2; j++) {
				phi_uh[j] = phi_uh[j + 2];
			}
			phi_uh.resize(this->Order_uh);

			int column_Pointlist = this->Pointlist.rowsize();
			this->uhphi_point.zeros(this->variablesize, column_Pointlist);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: uh_i * phi_i with PointList----------------------------" << std::endl;
			std::cout << ">> uh size is (" << list_uh.rowsize() << "," << this->variablesize << ")" << std::endl;
			std::cout << ">> uh_i * phi_i size is (" << this->variablesize << "," << column_Pointlist << ")\n" << std::endl;
#endif

			for (int v = 0; v < this->variablesize; v++) {
				for (int j = 0; j < this->Pointlist.rowsize(); j++) {
					for (int i = 0; i < list_uh.rowsize(); i++) {
						_TM Phi_Time = Legendre_uh_with_GLpoint(list_uh(i, 0), this->Pointlist(j, 0));
						for (int d = 1; d < list_uh.columnsize(); d++) {
							Phi_Time *= Legendre_uh_with_GLpoint(list_uh(i, d), this->Pointlist(j, d));
						}
						this->uhphi_point(v, j) += uh(i, v) * Phi_Time;
					}
				}
			}

		}

		void cxphi_point_set(const matrix< int >& list_cx, int& order_cx) {
			this->Order_cx = order_cx;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base------------------------------------------" << std::endl;
			std::cout << ">> Order of Legendre base for cx is " << this->Order_cx << std::endl;
			std::cout << ">> Order of polynomial is " << this->Order_cx + 2 << " (= Order_cx + 2)\n" << std::endl;
#endif
			std::vector< kv::psa< _T > > phi_cx;
			LegendreBaseFunctions::makeLegendreBase(this->Order_cx + 2, phi_cx);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base <== Point of Gauss Legendre--------------" << std::endl;
			std::cout << ">> Size of Legendre_with_GLpoint is (" << this->Order_cx << "," << this->Point.size() << ")" << std::endl;
			std::cout << ">> Size of DifDifLegendre_with_GLpoint is (" << this->Order_cx << "," << this->Point.size() << ")\n" << std::endl;
#endif
			std::vector< std::vector< _T > > LBP2;
			LegendreBaseFunctions::LegendrePointFunc(phi_cx, this->Point, LBP2);

			matrix< _TM, _PM > Legendre_cx_with_GLpoint;
			Legendre_cx_with_GLpoint.zeros(this->Order_cx, this->Point.size());

			for (int i = 0; i < this->Order_cx; i++) {
				for (int j = 0; j < this->Point.size(); j++) {
					convert(LBP2[i + 2][j], Legendre_cx_with_GLpoint(i, j));
//					Legendre_cx_with_GLpoint(i, j) = LBP2[i + 2][j];
				}
			}
			for (int j = 0; j < phi_cx.size() - 2; j++) {
				phi_cx[j] = phi_cx[j + 2];
			}
			phi_cx.resize(this->Order_cx);

			int column_Pointlist = this->Pointlist.rowsize();
			this->cxphi_point.zeros(this->variablesize, column_Pointlist);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: cx_i * phi_i with PointList----------------------------" << std::endl;
			std::cout << ">> cx size is (" << list_cx.rowsize() << "," << this->variablesize << ")" << std::endl;
			std::cout << ">> cx_i * phi_i size is (" << this->variablesize << "," << column_Pointlist << ")\n" << std::endl;
#endif

			for (int v = 0; v < this->variablesize_cx; v++) {
				for (int j = 0; j < this->Pointlist.rowsize(); j++) {
					for (int i = 0; i < list_cx.rowsize(); i++) {
						_TM Phi_Time = Legendre_cx_with_GLpoint(list_cx(i, 0), this->Pointlist(j, 0));
						for (int d = 1; d < list_cx.columnsize(); d++) {
							Phi_Time *= Legendre_cx_with_GLpoint(list_cx(i, d), this->Pointlist(j, d));
						}
						this->cxphi_point(v, j) += cx(i, v) * Phi_Time;
					}
				}
			}

		}

		template<typename... Args> void make_uh_p(matrix< _TM, _PM >& uh_p_phi_point) {
		}
		template<typename... Args> void make_uh_p(matrix< _TM, _PM >& uh_p_phi_point, const int& p, const Args&... args) {
			using std::pow;
			constexpr std::size_t vari_size = sizeof...(args);
			if (p == 0) {
				this->make_uh_p(uh_p_phi_point, args...);
			}
			else if (p < 0) {
				vcp::throw_error<vcp::invalid_argument>(">> ERROR : make_uh_p : input arguments is negative ( p = ", p, " )");
			}
			else {
				for (int i = 0; i < this->Pointlist.rowsize(); i++) {
					uh_p_phi_point(0, i) *= pow(this->uhphi_point(this->variablesize - (vari_size + 1), i), p);
				}
				this->make_uh_p(uh_p_phi_point, args...);
			}
		}

		template<typename... Args> void make_cx_p(matrix< _TM, _PM >& cx_p_phi_point) {
		}
		template<typename... Args> void make_cx_p(matrix< _TM, _PM >& cx_p_phi_point, const int& p, const Args&... args) {
			using std::pow;
			constexpr std::size_t vari_size = sizeof...(args);
			if (p == 0) {
				this->make_cx_p(cx_p_phi_point, args...);
			}
			else if (p < 0) {
				vcp::throw_error<vcp::invalid_argument>(">> ERROR : make_cx_p : input arguments is negative ( p = ", p, " )");
			}
			else {
				for (int i = 0; i < this->Pointlist.rowsize(); i++) {
					cx_p_phi_point(0, i) *= pow(this->cxphi_point(this->variablesize_cx - (vari_size + 1), i), p);
				}
				this->make_cx_p(cx_p_phi_point, args...);
			}
		}

		template<typename... Args> void disp_args(int& vs) {
		}
		template<typename... Args> void disp_args(int& vs, const int& p, const Args&... args) {
			vs += p;
#ifdef VCP_LEGENDRE_DEBUG
			constexpr std::size_t vari_size = sizeof...(args);
			int kk = this->variablesize - vari_size;
			std::cout << "u" << kk << "^" << p << "*";
#endif
			this->disp_args(vs, args...);
		}

		template <class _TC, class _PTCM = mats< _TC >> matrix< _TC, _PTCM> func(std::vector< _TC > x) {
			int mmax = this->phi.size();
			int dim = this->dimension;
			int var = this->variablesize;

			if (x.size() != dim) {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : func : dimension size != vector size");
			}

			vcp::matrix< _T > p_phi_T;
			vcp::matrix< _TC > p_phi_TC;

			p_phi_T.zeros(mmax, dim);
			for (int j = 0; j < dim; j++) {
				for (int i = 0; i < mmax; i++) {
					_T temp;
					vcp::convert(x[j], temp);
					psa_value(this->phi[i], temp, p_phi_T(i, j));
				}
			}
			convert(p_phi_T, p_phi_TC);
			matrix< _TC, _PTCM> output;
			output.zeros(1, var + dim);

			for (int i = 0; i < dim; i++) {
				output(0, i) = x[i];
			}
			for (int v = 0; v < var; v++) {
				_TC tmp;
				for (int k = 0; k < this->list.rowsize(); k++) {
					tmp = uh(k, v);
					for (int di = 0; di < dim; di++) {
						tmp *= p_phi_TC(this->list(k, di), di);
					}
					output(0, v + dim) += tmp;
				}
			}
			return output;
		}

		template < class _TC > _TC func(std::vector< _TC > x, int v) {
			int mmax = this->phi.size();
			int dim = this->dimension;
			int var = this->variablesize;

			if (x.size() != dim) {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : func : dimension size != vector size");
			}

			vcp::matrix< _T > p_phi_T;
			vcp::matrix< _TC > p_phi_TC;

			p_phi_T.zeros(mmax, dim);
			for (int j = 0; j < dim; j++) {
				for (int i = 0; i < mmax; i++) {
					_T temp;
					vcp::convert(x[j], temp);
					psa_value(this->phi[i], temp, p_phi_T(i, j));
				}
			}
			convert(p_phi_T, p_phi_TC);
			_TC output;
			output = _TC(0);


			_TC tmp;
			for (int k = 0; k < this->list.rowsize(); k++) {
				//tmp = uh(k, v);
				convert(uh(k,v), tmp);
				for (int di = 0; di < dim; di++) {
					tmp *= p_phi_TC(this->list(k, di), di);
				}
				output += tmp;
			}
			return output;
		}

	public:
		Legendre_Bases_Generator< _T, _TM, _PM >() {
			this->flag_order_uh = false;
		}

		void setting(const int mm, const int pp, const int dim, const int vs = 1, const int k = 2, const int uh_o = -1) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////////// Legendre_Bases_Generator :: setting(" << mm << "," << pp << "," << dim << "," << vs << "," << k << "," << uh_o << ") //////////////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			if (mm < 2 || pp < 1 || dim < 1 || vs < 1 || k < 0 || uh_o < -1 || uh_o == 0) {
				vcp::throw_error<vcp::invalid_argument>(">> ERROR : setting : Not Appropriate Input Arguments");
			}
			if (uh_o > 1 && (k != 1 && k != 0 ) ) {
				vcp::throw_error<vcp::state_error>(">> ERROR : Order of uh > 1 && Mode isn't With out Residual Mode  ");
			}
			if (uh_o > 1) {
				if (k == 0){
					this->flag_order_cx = true;
					this->Order_cx = uh_o;
				}
				else{
					this->flag_order_uh = true;
					this->Order_uh = uh_o;
				}

			}
			else {
				this->flag_order_uh = false;
				this->Order_uh = mm;
			}

			int m = mm + 2;

			this->Order_of_Base = mm;
			this->p = pp;
			this->dimension = dim;
			this->mode = k;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Gauss-Legendre Numerical Integration-------------------" << std::endl;
#endif
			int n;
			if (k == 1) {
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> Warning : You can only use (uh^p, phi)_{L^2} and (uh^(p-1) phi, phi)_{L^2}" << std::endl;
				std::cout << ">> Don't calculate the residual (uh^p, uh^p)_{L^2}" << std::endl;
#endif
				if (this->flag_order_uh) {
					//  2*nn1-1 > (p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2)   -->  n := ((p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2) + 1)/2 + 1
					int nn1 = ((p - 1)*(this->Order_uh + 2) + 2 * (this->Order_of_Base + 2) + 2) / 2 + 2;
					//  2*nn2-1 > p*(Order_uh + 2) + (Order_of_Base + 2)   -->  n := ((p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2) + 1)/2 + 1
					int nn2 = (p*(this->Order_uh + 2) + (this->Order_of_Base + 2) + 2) / 2 + 2;
					n = std::max(nn1, nn2);
					if (n % 2 != 0) {
						n = n + 1;
					}
#ifdef VCP_LEGENDRE_DEBUG
					std::cout << ">> Order of Gauss Legendre Integration is " << n << " = std::max(" << nn1 << "," << nn2 << ")\n" << std::endl;
#endif
				}
				else {
					//  2*n-1 > k * (Order_of_Base + 2) *(p + 1)  -->  n := (k * (Order_of_Base + 2) * p + 1)/2 + 1
					n = (k*(this->Order_of_Base + 2) * (p + 1) + 1) / 2 + 2;
					if (n % 2 != 0) {
						n = n + 1;
					}
#ifdef VCP_LEGENDRE_DEBUG
					std::cout << ">> Order of Gauss Legendre Integration is " << n << " = (" << k << "*(Order_of_Base + 2) * (p + 1) + 1)/2 + 2 \n" << std::endl;
#endif
				}
			}
			else if (k == 2) {
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> You can use the residual (uh^p, uh^p)_{L^2}" << std::endl;
#endif
				//  2*n-1 > k * (Order_of_Base + 2) * p  -->  n := (k * (Order_of_Base + 2) * p + 1)/2 + 1
				n = (k*(this->Order_of_Base + 2) * p + 1) / 2 + 2;
				if (n % 2 != 0) {
					n = n + 1;
				}
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> Order of Gauss Legendre Integration is " << n << " = (" << k << "*(Order_of_Base + 2) * p + 1)/2 + 2 \n" << std::endl;
#endif
			}
			else if (k == 0){
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> You can use the Goerisch theorem (L uhi, L uhj)_{L^2}" << std::endl;
#endif
				if (this->flag_order_cx) {
					//  2*nn1-1 > (p - 1)*(Order_cx + 2) + 2*(Order_of_Base + 2)   -->  n := ((p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2) + 1)/2 + 1
					int nn1 = (2*(p - 1)*(this->Order_cx + 2) + 2 * (this->Order_of_Base + 2) + 2) / 2 + 2;
					//  2*nn2-1 > p*(Order_uh + 2) + (Order_of_Base + 2)   -->  n := ((p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2) + 1)/2 + 1
					int nn2 = (2*p*(this->Order_cx + 2) + (this->Order_of_Base + 2) + 2) / 2 + 2;
					n = std::max(nn1, nn2);
				}
				else{
					//  2*n-1 > 2*(Order_of_Base + 2)  -->  n := ( 2*(Order_of_Base + 2) + 1 ) /2 + 2;
					n = (2*(this->Order_of_Base + 2) + 1) / 2 + 2;
				}
				if (n % 2 != 0) {
					n = n + 1;
				}
			}
			else {
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> Approximation mode = " << k << std::endl;
#endif
				n = k * 10;
				if (n % 2 != 0) {
					n = n + 1;
				}
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> Order of Gauss Legendre Integration is " << n << " = " << k << "*10 \n" << std::endl;
#endif
			}
			this->interval_ld_weightpoint< _T >::set(n);

			for (int i = 0; i < n; i++) {
				this->Point[i] = (this->Point[i] + 1) / 2;
				this->Weight[i] = this->Weight[i] / 2;
			}
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base------------------------------------------" << std::endl;
			std::cout << ">> Order of Legendre base m is " << mm << std::endl;
			std::cout << ">> Order of polynomial is " << mm + 2 << " (= m + 2)\n" << std::endl;
#endif
			LegendreBaseFunctions::makeLegendreBase(m, this->phi);
			LegendreBaseFunctions::psaTodpsa(this->phi, this->dphi);
			LegendreBaseFunctions::psaTodpsa(this->dphi, this->ddphi);


#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base <== Point of Gauss Legendre--------------" << std::endl;
			std::cout << ">> Size of Legendre_with_GLpoint is (" << mm << "," << this->Point.size() << ")" << std::endl;
			std::cout << ">> Size of DifDifLegendre_with_GLpoint is (" << mm << "," << this->Point.size() << ")\n" << std::endl;
#endif
			std::vector< std::vector< _T > > LBP2;
			std::vector< std::vector< _T > > LBDDP2;
			LegendreBaseFunctions::LegendrePointFunc(this->phi, this->Point, LBP2);
			LegendreBaseFunctions::LegendrePointFunc(this->ddphi, this->Point, LBDDP2);

			this->Legendre_with_GLpoint.zeros(mm, this->Point.size());
			this->DifDifLegendre_with_GLpoint.zeros(mm, this->Point.size());
			for (int i = 0; i < mm; i++) {
				for (int j = 0; j < this->Point.size(); j++) {
					_T a = LBP2[i + 2][j];
					_TM b;
					convert(a,b);
					this->Legendre_with_GLpoint(i, j) = b;

					a = LBDDP2[i + 2][j];
					convert(a, b);
					this->DifDifLegendre_with_GLpoint(i, j) = b;
				}
			}
			for (int j = 0; j < phi.size() - 2; j++) {
				phi[j] = phi[j + 2];
			}
			phi.resize(mm);

			this->Pointlist_set();
			this->weight_point_set();
			this->variablesize = vs;
		}

		void setting_list() {
			this->list_set();
			this->elementsize = this->list.rowsize();
			if ((this->mode != 2) && (this->mode != 0) ) {
				this->phi_point_set();
			}
		}

		void setting_evenlist() {
			this->even_list_set();
			this->elementsize = this->list.rowsize();
			if ((this->mode != 2) && (this->mode != 0) ) {
				this->phi_point_set();
			}
		}

		void setting_uh() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "///////////////////// Legendre_Bases_Generator :: setting_uh() /////////////////////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			if (this->flag_order_uh) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : setting_uh : Flag for order of uh is True...\n",
					">> Please check the Last Argument of the function .setting()");
			}

			this->uh.ones(this->elementsize, this->variablesize);
		//	if (this->mode != 2) {
				this->uhphi_point_set();
		//	}
		}
		void setting_uh(const matrix< _TM, _PM >& uuh) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////// Legendre_Bases_Generator :: setting_uh(vcp::matrix uh) //////////////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			if (this->flag_order_uh) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : setting_uh : Flag for order of uh is True...\n",
					">> Please check the Last Argument of the function .setting()");
			}

			if (uuh.rowsize() != this->elementsize || uuh.columnsize() != this->variablesize) {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : setting_uh : no much the matrix size : uh");
			}
			this->uh = uuh;
		//	if (this->mode != 2) {
				this->uhphi_point_set();
		//	}
		}
		void setting_uh(const matrix< _TM, _PM >& uuh, const matrix< int >& list_uh, int ind) {
			// ind: full list => 1 , even list => 2
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//// Legendre_Bases_Generator :: setting_uh(matrix<_T> uh, matrix<int> list_uh) ////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			matrix< int > uh_o = max(max(list_uh));
			int uh_order = uh_o(0) + ind;

			if (!this->flag_order_uh) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : setting_uh : Flag for order of uh is False...\n",
					">> Please check the Last Argument of the function .setting()");
			}
			if (this->Order_uh != uh_order) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : setting_uh : no much the uh_order : ", this->Order_uh, "!=", uh_order, "\n",
					">> Please check the Last Argument of the function .setting() or Order of the Argument uh");
			}
			
			const int expected_uh_rows = LegendreBaseFunctions::checked_int_pow(uh_order / ind, this->dimension, "setting_uh");
			if (uuh.rowsize() != expected_uh_rows || uuh.columnsize() != this->variablesize) {
				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : setting_uh : no much the matrix size : uh : ",
					uuh.rowsize(), "!=", expected_uh_rows, "||",
					uuh.columnsize(), "!=", this->variablesize);
			}
		
			this->uh = uuh;
			if (this->mode != 2) {
				this->uhphi_point_set(list_uh, uh_order);
			}
		}

		void setting_cx(const matrix< _TM, _PM >& ccx, const matrix< int >& list_cx, int ind) {
			// ind: full list => 1 , even list => 2
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//// Legendre_Bases_Generator :: setting_cx(matrix<_T> cx, matrix<int> list_cx) ////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			matrix< int > cx_o = max(max(list_cx));
			int cx_order = cx_o(0) + ind;

			if (!this->flag_order_cx) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : setting_cx : Flag for order of cx is False...\n",
					">> Please check the Last Argument of the function .setting()");
			}
			if (this->Order_cx != cx_order) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : setting_cx : no much the cx_order : ", this->Order_cx, "!=", cx_order, "\n",
					">> Please check the Last Argument of the function .setting() or Order of the Argument cx");
			}
			this->variablesize_cx = ccx.columnsize();
			const int expected_cx_rows = LegendreBaseFunctions::checked_int_pow(cx_order / ind, this->dimension, "setting_cx");
			if (ccx.rowsize() != expected_cx_rows) {
				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : setting_cx : no much the matrix size : cx : ",
					ccx.rowsize(), "!=", expected_cx_rows);
			}
		
			this->cx = ccx;
			if (this->mode != 2) {
				this->cxphi_point_set(list_cx, cx_order);
			}
		}

		template<typename... Args> matrix< _TM, _PM > uhphiphi(const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////// Legendre_Bases_Generator :: uhphiphi(Args... )//////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif
			constexpr std::size_t vari_size = sizeof...(args);
			if ( (this->mode == 2) || (this->mode == 0) ) {
				vcp::throw_error<vcp::state_error>(">> ERROR : uhphiphi : This function can not use the Residual mode (Mode = 2)");
			}

			if (vari_size != this->variablesize) {
				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : uhphiphi : size of input arguments != variablesize ( ",
					vari_size, " != ", this->variablesize, " )");
			}

			matrix< _TM, _PM > uh_p_phi_point;
			uh_p_phi_point.ones(1, this->Pointlist.rowsize());
			this->make_uh_p(uh_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the Matrix { (";
#endif
			this->disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "phi_i, phi_j )_{L^2} }_{i,j} \n" << std::endl;
#endif
			if (this->p < vs) {
				vcp::throw_error<vcp::invalid_argument>(">> ERROR : uhphiphi : p < vs( p=", this->p, ", vs=", vs, ")");
			}

			matrix< _TM, _PM > A;
			A.zeros(this->elementsize, this->elementsize);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < this->elementsize; i++) {
				for (int j = i; j < this->elementsize; j++) {
					for (int k = 0; k < this->Pointlist.rowsize(); k++) {
						A(i, j) += this->weight_point(0, k) * uh_p_phi_point(0, k) * this->phi_point(i, k) * this->phi_point(j, k);
					}
				}
			}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 1; i < this->elementsize; i++) {
				for (int j = 0; j < i; j++) {
					A(i, j) = A(j, i);
				}
			}
			return A;
		}
		
		template<typename... Args> matrix< _TM, _PM > uhphi(const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////////// Legendre_Bases_Generator :: uhphi(Args... )////////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif
			constexpr std::size_t vari_size = sizeof...(args);
			if ( (this->mode == 2) || (this->mode == 0) ) {
				vcp::throw_error<vcp::state_error>(">> ERROR : uhphi : This function can not use the Residual mode (Mode = 2)");
			}

			if (vari_size != this->variablesize) {

				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : uhphi : size of input arguments != variablesize ( ",
					vari_size, " != ", this->variablesize, " )");
			}

			matrix< _TM, _PM > uh_p_phi_point;
			uh_p_phi_point.ones(1, this->Pointlist.rowsize());
			this->make_uh_p(uh_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the Vector { (";
#endif
			this->disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ", phi_j )_{L^2} }_{i,j} \n" << std::endl;
#endif
			if (this->p < vs) {
				vcp::throw_error<vcp::invalid_argument>(">> ERROR : uhphi : p < vs( p=", this->p, ", vs=", vs, ")");
			}

			matrix< _TM, _PM > A;
			A.zeros(this->elementsize, 1);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < this->elementsize; i++) {
				for (int k = 0; k < this->Pointlist.rowsize(); k++) {
					A(i, 0) += this->weight_point(0, k) * uh_p_phi_point(0, k) * this->phi_point(i, k);
				}
			}
			return A;
		}

		matrix< _TM, _PM > phiphi() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////////////// Legendre_Bases_Generator :: phiphi()///////////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			if ( (this->mode == 2) || (this->mode == 0) ) {
				vcp::throw_error<vcp::state_error>(">> ERROR : phiphi : This function can not use the Residual mode (Mode = 2)");
			}
			matrix< _TM, _PM > L;
			L.zeros(this->elementsize, this->elementsize);

			matrix< int > lmax = max(max(this->list));
			int listmax = lmax(0) + 1;
			matrix< _TM, _PM > OneDimL;
			OneDimL.zeros(listmax, listmax);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < listmax; i++){
				for (int j = 0; j < listmax; j ++){
					if (i == j){
						OneDimL(i, j) = 1 / (_TM(2)*(2 * _TM(i) + 1)*(2 * _TM(i) + 5)*(2 * _TM(i) + 3));
					}
					else if (i - j == 2) {
						OneDimL(i, j) = -1 / (4 * (2 * _TM(j) + 5)*(2 * _TM(j) + 7)*(2 * _TM(j) + 3));
					}
					else if (j - i == 2) {
						OneDimL(i, j) = -1 / (4 * (2 * _TM(i) + 5)*(2 * _TM(i) + 7)*(2 * _TM(i) + 3));
					}
				}
			}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < this->elementsize; i++) {
				for (int j = i; j < this->elementsize; j++) {
					_TM Phi_Time = _TM(1);
					for (int d = 0; d < this->list.columnsize(); d++) {
						Phi_Time *= OneDimL(list(i, d), list(j, d));
					}
					L(i, j) = Phi_Time;
				}
			}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 1; i < this->elementsize; i++) {
				for (int j = 0; j < i; j++) {
					L(i, j) = L(j, i);
				}
			}
			return L;
		}

		matrix< _TM, _PM > fphi() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////////// Legendre_Bases_Generator :: fphi()////////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif
			if ( (this->mode == 2) || (this->mode == 0) ) {
				vcp::throw_error<vcp::state_error>(">> ERROR : uhphi : This function can not use the Residual mode (Mode = 2)");
			}
			matrix< _TM, _PM > A;
			A.zeros(this->elementsize, 1);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < this->elementsize; i++) {
				for (int k = 0; k < this->Pointlist.rowsize(); k++) {
					A(i, 0) += this->weight_point(0, k) * this->phi_point(i, k);
				}
			}
			return A;

		}

		matrix< _TM, _PM > dphidphi() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////// Legendre_Bases_Generator :: dphidhi()///////////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			if ( (this->mode == 2) || (this->mode == 0) ) {
				vcp::throw_error<vcp::state_error>(">> ERROR : dphidphi : This function can not use the Residual mode (Mode = 2)");
			}
			matrix< _TM, _PM > DL;
			DL.zeros(this->elementsize, this->elementsize);
			matrix< int > lmax = max(max(this->list));
			int listmax = lmax(0) + 1;
			matrix< _TM, _PM > OneDimL;
			OneDimL.zeros(listmax, listmax);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < listmax; i++) {
				for (int j = 0; j < listmax; j++) {
					if (i == j) {
						OneDimL(i, j) = 1 / (_TM(2)*(2 * _TM(i) + 1)*(2 * _TM(i) + 5)*(2 * _TM(i) + 3));
					}
					else if (i - j == 2) {
						OneDimL(i, j) = -1 / (4 * (2 * _TM(j) + 5)*(2 * _TM(j) + 7)*(2 * _TM(j) + 3));
					}
					else if (j - i == 2) {
						OneDimL(i, j) = -1 / (4 * (2 * _TM(i) + 5)*(2 * _TM(i) + 7)*(2 * _TM(i) + 3));
					}
				}
			}
			matrix< _TM, _PM > OneDimDL;
			OneDimDL.zeros(listmax, listmax);
			for (int i = 0; i < listmax; i++) {
				OneDimDL(i,i) = 1 / (2 * _TM(i) + 3);
			}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif			
			for (int i = 0; i < this->elementsize; i++) {
				for (int j = i; j < this->elementsize; j++) {
					for (int dc = 0; dc < this->list.columnsize(); dc++) {
						_TM Phi_Time = _TM(1);
						for (int dr = 0; dr < this->list.columnsize(); dr++) {
							if (dc == dr) {
								Phi_Time *= OneDimDL(list(i, dr), list(j, dr));
							}
							else {
								Phi_Time *= OneDimL(list(i, dr), list(j, dr));
							}
						}
						DL(i, j) += Phi_Time;
					}
				}
			}

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 1; i < this->elementsize; i++) {
				for (int j = 0; j < i; j++) {
					DL(i, j) = DL(j, i);
				}
			}
			return DL;
		}

		matrix< _TM, _PM > output_uh_for_graphics(const int Div_number = 100) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////////// Legendre_Bases_Generator :: output_uh_for_graphics()/////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif	

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Div_number =" << Div_number << std::endl;
#endif	
			if (this->mode <= 2) {
				vcp::throw_error<vcp::state_error>(">> ERROR : output_uh__for_graphics : This function only use the Approximation mode (Mode > 2)");
			}

			matrix< _TM, _PM > Output_Data;
			int dim = this->dimension;
			const std::vector<int> output_powers = LegendreBaseFunctions::checked_int_power_table(Div_number + 1, dim, "output_uh_for_graphics");
			int Output_data_row_size = output_powers[dim];
			matrix< int > list_of_Output_data;
			
			list_of_Output_data.zeros(Output_data_row_size, dim);
			
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int d = 0; d < dim; d++) {
				int di = dim - d;
				const int lower_block = output_powers[di - 1];
				const int upper_block = output_powers[di];
				for (int j = 0; j < output_powers[d]; j++) {
					for (int i = 0; i < Div_number + 1; i++) {
						for (int k = 0; k < lower_block; k++) {
							list_of_Output_data(i*lower_block + j * upper_block + k, d) = i;
						}
					}
				}
			}
			
			_T mesh_size_T = 1 / _T(Div_number);
			_TM mesh_size = 1 / _TM(Div_number);
			int mmax = this->phi.size();
			vcp::matrix< _T > p_phi_T;
			vcp::matrix< _TM, _PM > p_phi;
			p_phi_T.zeros(mmax, Div_number + 1);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < mmax; i++) {
				for (int j = 0; j < Div_number + 1; j++){
					psa_value_Horner(phi[i],j*mesh_size_T,p_phi_T(i, j));
				}
			}

			convert(p_phi_T, p_phi);
			Output_Data.zeros(Output_data_row_size,this->variablesize + dim);

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < Output_data_row_size; i++) {
				for (int di = 0; di < dim; di++){
					Output_Data(i, di) = list_of_Output_data(i, di)*mesh_size;
				}
				for (int v = 0; v < this->variablesize; v++) {
					_TM tmp;
					for (int k = 0; k < this->list.rowsize(); k++) {
						tmp = uh(k, v);
						for (int di = 0; di < dim; di++) {
							tmp *= p_phi(this->list(k, di), list_of_Output_data(i, di));
						}
						Output_Data(i, dim + v) += tmp;
					}
				}
			}

			return Output_Data;
		}

		template <class _TC > std::vector< _TC > global_min(const std::vector< kv::interval< _TC > > x, double minimal_mesh_size = std::pow(2.0,-11)){
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////////////// Legendre_Bases_Generator :: global_min() //////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif	

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> minimal_mesh_size = " << minimal_mesh_size << std::endl;
#endif

			int dim = this->dimension;
			int var = this->variablesize;

			if (this->mode == 2) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : global_min : This function can use the approximation mode or Inverse Mode (Mode = 1) (with verification)");
			}

			if (x.size() != dim) {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : global_min : x.size() != dim (", x.size(), " != ", dim, ")");
			}

			std::vector< _TC > gmin;
			gmin.resize(var);
			for (int i = 0; i < var; i++) {
				gmin[i] = _TC(INFINITY);
			}

			std::vector< kv::interval< _TC > > local_vector;


			matrix< int > div_list;
			const std::vector<int> div_powers = LegendreBaseFunctions::checked_int_power_table(2, dim, "global_min");
			int div_list_size = div_powers[dim];
			div_list.zeros(div_list_size, dim);
			for (int i = 0; i < 2; i++) {
				for (int d = 0; d < dim; d++) {
					int di = dim - d;
					const int lower_block = div_powers[di - 1];
					const int upper_block = div_powers[di];
					for (int j = 0; j < div_powers[d]; j++) {
						for (int k = 0; k < lower_block; k++) {
							div_list(i*lower_block + j * upper_block + k, d) = i;
						}
					}
				}
			}
			for (int v = 0; v < var; v++) {
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << "\n-------------- Variable : " <<  v + 1 << " Start -------------------" << std::endl;
#endif
				std::stack< std::vector< kv::interval< _TC > > > list_stack;
				local_vector.resize(dim);
				for (int i = 0; i < dim; i++) {
					local_vector[i] = x[i];
				}
				list_stack.push(local_vector);
				for (int i = 0; i < dim; i++) {
					local_vector[i] = kv::interval< _TC >(mid(local_vector[i]));
				}
				_TC delta = (this->template func<  kv::interval< _TC > >(local_vector, v)).lower();
				int active_workers = 0;
				bool done = false;

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel
			{
#endif
#endif	
				while (true) {
					std::vector< kv::interval< _TC > > local_point;
					bool have_work = false;
					bool should_break = false;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
				#pragma omp critical (global_min_work)
				{
#endif
#endif
					if (!list_stack.empty()) {
						local_point = list_stack.top();
						list_stack.pop();
						active_workers++;
						have_work = true;
					}
					else if (active_workers == 0) {
						done = true;
						should_break = true;
					}
					else if (done) {
						should_break = true;
					}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
				}
#endif
#endif
					if (should_break) {
						break;
					}
					if (!have_work) {
						continue;
					}
					
					kv::interval< _TC > tmp = this->template func<  kv::interval< _TC > >(local_point, v);
					_TC local_delta;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
					#pragma omp critical (global_min_delta)
					{
#endif
#endif
						local_delta = delta;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
					}
#endif
#endif
					
					if (local_delta < tmp.lower()) {
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						#pragma omp critical (global_min_work)
						{
#endif
#endif
							active_workers--;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						}
#endif
#endif
						continue;
					}
					if (width(local_point[0]) < minimal_mesh_size) {

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
					#pragma omp critical (global_min_result)
					{
#endif
#endif
						if (tmp.lower() < gmin[v]) {
#ifdef VCP_LEGENDRE_DEBUG
							std::cout << ">> Find candidate date of global_min : " << tmp << std::endl;
#endif
							gmin[v] = tmp.lower();
						}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
					}
#endif
#endif

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						#pragma omp critical (global_min_work)
						{
#endif
#endif
							active_workers--;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						}
#endif
#endif
						continue;
					}

					std::vector< kv::interval< _TC > > local_point_center = local_point;
					for (int i = 0; i < dim; i++) {
						local_point_center[i] = kv::interval< _TC >(mid(local_point_center[i]));
					}
					tmp = (this->template func<  kv::interval< _TC > >(local_point_center, v));
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
					#pragma omp critical (global_min_delta)
					{
#endif
#endif
						if (tmp.lower() < delta) {
							delta = tmp.lower();
						}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
					}
#endif
#endif

					std::vector< std::vector< kv::interval< _TC > > > div_interval_list;
					div_interval_list.resize(2);
					div_interval_list[0].resize(dim);
					div_interval_list[1].resize(dim);
					for (int d = 0; d < dim; d++) {
						div_interval_list[0][d] = kv::interval< _TC >(local_point[d].lower(), mid(local_point[d]));
						div_interval_list[1][d] = kv::interval< _TC >(mid(local_point[d]), local_point[d].upper());
					}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
				#pragma omp critical (global_min_work)
				{
#endif
#endif
					for (int i = 0; i < div_list_size; i++) {
						for (int d = 0; d < dim; d++) {
							local_point[d] = div_interval_list[div_list(i, d)][d];
						}
						list_stack.push(local_point);
					}					
					active_workers--;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
				}
#endif
#endif
				}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			}
#endif
#endif	

			}
			return gmin;
		}

		template <class _TC > std::vector< _TC > global_max(const std::vector< kv::interval< _TC > > x, double minimal_mesh_size = std::pow(2.0, -11)) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////////////// Legendre_Bases_Generator :: global_max() //////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif	

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> minimal_mesh_size = " << minimal_mesh_size << std::endl;
#endif
			
			int dim = this->dimension;
			int var = this->variablesize;

			if (this->mode == 2) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : global_max : This function can only use the approximation mode or Inverse Mode (Mode = 1) (with verification)");
			}

			if (x.size() != dim) {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : global_max : x.size() != dim (", x.size(), " != ", dim, ")");
			}

			std::vector< _TC > gmax;
			gmax.resize(var);
			for (int i = 0; i < var; i++) {
				gmax[i] = _TC(-INFINITY);
			}

			std::vector< kv::interval< _TC > > local_vector;


			matrix< int > div_list;
			const std::vector<int> div_powers = LegendreBaseFunctions::checked_int_power_table(2, dim, "global_max");
			int div_list_size = div_powers[dim];
			div_list.zeros(div_list_size, dim);
			for (int i = 0; i < 2; i++) {
				for (int d = 0; d < dim; d++) {
					int di = dim - d;
					const int lower_block = div_powers[di - 1];
					const int upper_block = div_powers[di];
					for (int j = 0; j < div_powers[d]; j++) {
						for (int k = 0; k < lower_block; k++) {
							div_list(i*lower_block + j * upper_block + k, d) = i;
						}
					}
				}
			}
			for (int v = 0; v < var; v++) {
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << "\n-------------- Variable : " << v + 1 << " Start -------------------" << std::endl;
#endif
				std::stack< std::vector< kv::interval< _TC > > > list_stack;
				local_vector.resize(dim);
				for (int i = 0; i < dim; i++) {
					local_vector[i] = x[i];
				}
				list_stack.push(local_vector);
				for (int i = 0; i < dim; i++) {
					local_vector[i] = kv::interval< _TC >(mid(local_vector[i]));
				}
				_TC delta = (this->template func<  kv::interval< _TC > >(local_vector, v)).upper();
				int active_workers = 0;
				bool done = false;

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel
				{
#endif
#endif	
					while (true) {
						std::vector< kv::interval< _TC > > local_point;
						bool have_work = false;
						bool should_break = false;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						#pragma omp critical (global_max_work)
						{
#endif
#endif
							if (!list_stack.empty()) {
								local_point = list_stack.top();
								list_stack.pop();
								active_workers++;
								have_work = true;
							}
							else if (active_workers == 0) {
								done = true;
								should_break = true;
							}
							else if (done) {
								should_break = true;
							}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						}
#endif
#endif
						if (should_break) {
							break;
						}
						if (!have_work) {
							continue;
						}

						kv::interval< _TC > tmp = this->template func<  kv::interval< _TC > >(local_point, v);
						_TC local_delta;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						#pragma omp critical (global_max_delta)
						{
#endif
#endif
							local_delta = delta;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						}
#endif
#endif

						if (local_delta > tmp.upper()) {
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
							#pragma omp critical (global_max_work)
							{
#endif
#endif
								active_workers--;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
							}
#endif
#endif
							continue;
					}
						if ( width(local_point[0]) < minimal_mesh_size ) {
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
							#pragma omp critical (global_max_result)
							{
#endif
#endif
								if (tmp.upper() > gmax[v]) {
#ifdef VCP_LEGENDRE_DEBUG
									std::cout << ">> Find candidate date of global_max : " << tmp << std::endl;
#endif
									gmax[v] = tmp.upper();
								}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
							}
#endif
#endif

#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
							#pragma omp critical (global_max_work)
							{
#endif
#endif
								active_workers--;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
							}
#endif
#endif
							continue;
						}

						std::vector< kv::interval< _TC > > local_point_center = local_point;
						for (int i = 0; i < dim; i++) {
							local_point_center[i] = kv::interval< _TC >(mid(local_point_center[i]));
						}
						tmp = (this->template func<  kv::interval< _TC > >(local_point_center, v));
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						#pragma omp critical (global_max_delta)
						{
#endif
#endif
							if (tmp.upper() > delta) {
								delta = tmp.upper();
							}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						}
#endif
#endif

						std::vector< std::vector< kv::interval< _TC > > > div_interval_list;
						div_interval_list.resize(2);
						div_interval_list[0].resize(dim);
						div_interval_list[1].resize(dim);
						for (int d = 0; d < dim; d++) {
							div_interval_list[0][d] = kv::interval< _TC >(local_point[d].lower(), mid(local_point[d]));
							div_interval_list[1][d] = kv::interval< _TC >(mid(local_point[d]), local_point[d].upper());
						}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						#pragma omp critical (global_max_work)
						{
#endif
#endif
							for (int i = 0; i < div_list_size; i++) {
								for (int d = 0; d < dim; d++) {
									local_point[d] = div_interval_list[div_list(i, d)][d];
								}
								list_stack.push(local_point);
							}
							active_workers--;
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
						}
#endif
#endif
				}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			}
#endif
#endif	

			}
			return gmax;
		}

		template<typename... Args> _TM integral_uh(const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "///////////////// Legendre_Bases_Generator :: integral_uh(Args... )/////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			constexpr std::size_t vari_size = sizeof...(args);
			if (this->mode != 2) {
				vcp::throw_error<vcp::state_error>(">> ERROR : integral_uh : This function only use the Residual mode (Mode = 2)");
			}

			if (vari_size != this->variablesize) {

				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : integral_uh : size of input arguments != variablesize ( ",
					vari_size, " != ", this->variablesize, " )");
			}

			matrix< _TM, _PM > uh_p_phi_point;
			uh_p_phi_point.ones(1, this->Pointlist.rowsize());
			this->make_uh_p(uh_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the value int(";
#endif
			this->disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ") dx \n" << std::endl;
#endif
			if (2 * this->p < vs) {
				vcp::throw_error<vcp::invalid_argument>(">> ERROR : integral_uh : 2*p < vs( 2*p=", 2 * this->p, ", vs=", vs, ")");
			}

			_TM res = _TM(0);

			for (int k = 0; k < this->Pointlist.rowsize(); k++) {
				res += this->weight_point(0, k) * uh_p_phi_point(0, k);
			}
			return res;
		}

		_TM integral_LuhLuh(int select_v) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "///////////// Legendre_Bases_Generator :: integral_LuhLuh(int select_v)/////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif
			if (this->mode != 2) {
				vcp::throw_error<vcp::state_error>(">> ERROR : integral_LuhLuh : This function only use the Residual mode (Mode = 2)");
			}

			if ( select_v + 1 >  this->variablesize ) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_LuhLuh : select valiable size != variablesize ( ",
					select_v + 1, " > ", this->variablesize, " )");
			}

			_TM res = _TM(0);
			using std::pow;

			for (int k = 0; k < this->Pointlist.rowsize(); k++) {
				res += this->weight_point(0, k) * pow( DDuhphi_point(select_v, k), 2);
			}
			return res;
		}

		_TM integral_LuhLuh(int select_v1, int select_v2) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////// Legendre_Bases_Generator :: integral_LuhLuh(v1, v2)   /////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif
			if ( ( this->mode != 2) && ( this->mode != 0) ) {
				vcp::throw_error<vcp::state_error>(
					">> ERROR : integral_LuhLuh : This function only use the Residual mode (Mode = 0 or 2) : ", this->mode);
			}

			if ( select_v1 + 1 >  this->variablesize ) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_LuhLuh : select valiable size != variablesize ( ",
					select_v1 + 1, " > ", this->variablesize, " )");
			}
			if ( select_v2 + 1 >  this->variablesize ) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_LuhLuh : select valiable size != variablesize ( ",
					select_v2 + 1, " > ", this->variablesize, " )");
			}

			_TM res = _TM(0);
			using std::pow;

			for (int k = 0; k < this->Pointlist.rowsize(); k++) {
				res += this->weight_point(0, k) * DDuhphi_point(select_v1, k)*DDuhphi_point(select_v2, k);
			}
			return res;
		}

		template<typename... Args> _TM integral_Luhuh(int select_v, const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////// Legendre_Bases_Generator :: integral_Luhuh(Args... )///////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			constexpr std::size_t vari_size = sizeof...(args);
			if (this->mode != 2) {
				vcp::throw_error<vcp::state_error>(">> ERROR : integral_Luhuh : This function only use the Residual mode (Mode = 2)");
			}

			if (vari_size != this->variablesize) {
				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : integral_Luhuh : size of input arguments != variablesize ( ",
					vari_size, " != ", this->variablesize, " )");
			}

			if (select_v + 1 >  this->variablesize) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_Luhuh : select valiable size != variablesize ( ",
					select_v + 1, " > ", this->variablesize, " )");
			}

			matrix< _TM, _PM > uh_p_phi_point;
			uh_p_phi_point.ones(1, this->Pointlist.rowsize());
			this->make_uh_p(uh_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the value int(";
#endif
			this->disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << " Laplace(uh)) dx \n" << std::endl;
#endif
			if (2*this->p < vs + 1) {
				vcp::throw_error<vcp::invalid_argument>(
					">> ERROR : integral_Luhuh : 2*p < vs + 1( 2*p=",
					2 * this->p, ", vs + 1=", vs + 1, ")");
			}

			_TM res = _TM(0);

			for (int k = 0; k < this->Pointlist.rowsize(); k++) {
				res += this->weight_point(0, k) * uh_p_phi_point(0, k) * DDuhphi_point(select_v, k);
			}
			return res;
		}

		template<typename... Args> _TM integral_Luhcxuh(int select_v_Luh, int select_v_uh, const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////// Legendre_Bases_Generator :: integral_Luhcxuh(Args... ) //////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			constexpr std::size_t vari_size = sizeof...(args);
			if (this->mode != 0) {
				vcp::throw_error<vcp::state_error>(">> ERROR : integral_Luhcxuh : This function only use the Goerisch mode (Mode = 0)");
			}

			if (vari_size != this->variablesize_cx) {
				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : integral_Luhcxuh : size of input arguments != variablesize_cx ( ",
					vari_size, " != ", this->variablesize_cx, " )");
			}

			if (select_v_Luh + 1 >  this->variablesize) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_Luhcxuh : select valiable size != variablesize ( ",
					select_v_Luh + 1, " > ", this->variablesize, " )");
			}

			if (select_v_uh + 1 >  this->variablesize) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_Luhcxuh : select valiable size != variablesize ( ",
					select_v_uh + 1, " > ", this->variablesize, " )");
			}

			matrix< _TM, _PM > cx_p_phi_point;
			cx_p_phi_point.ones(1, this->Pointlist.rowsize());
			this->make_cx_p(cx_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the value int(";
#endif
			this->disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << " Laplace(uh) uh) dx \n" << std::endl;
#endif
			if ( (this->p + 1) < vs + 1) {
				vcp::throw_error<vcp::invalid_argument>(
					">> ERROR : integral_Luhcxuh : p + 1 < vs + 1( p + 1=",
					this->p + 1, ", vs + 1=", vs + 1, ")");
			}

			_TM res = _TM(0);

			for (int k = 0; k < this->Pointlist.rowsize(); k++) {
				res += this->weight_point(0, k) * cx_p_phi_point(0, k) * DDuhphi_point(select_v_Luh, k) * uhphi_point(select_v_uh, k);
			}
			return res;
		}

		template<typename... Args> _TM integral_cxuhuh(int select_v_uh1, int select_v_uh2, const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////// Legendre_Bases_Generator :: integral_cxuhuh(Args... ) ///////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			constexpr std::size_t vari_size = sizeof...(args);
			if (this->mode != 0) {
				vcp::throw_error<vcp::state_error>(">> ERROR : integral_cxuhuh : This function only use the Goerisch mode (Mode = 0)");
			}

			if (vari_size != this->variablesize_cx) {
				vcp::throw_error<vcp::dimension_error>(
					">> ERROR : integral_cxuhuh : size of input arguments != variablesize_cx ( ",
					vari_size, " != ", this->variablesize_cx, " )");
			}

			if (select_v_uh1 + 1 >  this->variablesize) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_cxuhuh : select valiable size != variablesize ( ",
					select_v_uh1 + 1, " > ", this->variablesize, " )");
			}

			if (select_v_uh2 + 1 >  this->variablesize) {
				vcp::throw_error<vcp::index_error>(
					">> ERROR : integral_cxuhuh : select valiable size != variablesize ( ",
					select_v_uh2 + 1, " > ", this->variablesize, " )");
			}

			matrix< _TM, _PM > cx_p_phi_point;
			cx_p_phi_point.ones(1, this->Pointlist.rowsize());
			this->make_cx_p(cx_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the value int(";
#endif
			this->disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << " cx*ui*uj dx \n" << std::endl;
#endif
			if ( 2*(this->p -1) + 2 < vs + 1) {
				vcp::throw_error<vcp::invalid_argument>(
					">> ERROR : integral_cxuhuh : 2*(p-1)+2 < vs + 1( 2*(p-1)+2=",
					2*(this->p -1) + 2, ", vs + 1=", vs + 1, ")");
			}

			_TM res = _TM(0);

			for (int k = 0; k < this->Pointlist.rowsize(); k++) {
				res += this->weight_point(0, k) * cx_p_phi_point(0, k) * uhphi_point(select_v_uh1, k) * uhphi_point(select_v_uh2, k);
			}
			return res;
		}

		template <class _TC > _TC Ritz_projection_error(const int NN = -1) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//////////////// Legendre_Bases_Generator :: Ritz_projection_error()////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << ">> || Nabla( u - Rh u ) ||_L2 <= CN || Laplace u ||_L2" << std::endl;
#endif
			int N;

			if (NN <= 0) {
				N = this->Order_of_Base;
			}
			else {
				N = NN;
			}
			using std::sqrt;
			_TC sqrt2N_3 = sqrt(_TC(2 * N + 3));
			_TC sqrt2N_7 = sqrt(_TC(2 * N + 7));
			_TC sqrt2N_11 = sqrt(_TC(2 * N + 11));

			_TC L1 = 1 / _TC(2 * (2 * N + 1)*(2 * N + 5)) + 1 / ((4 * (2 * N + 5))*sqrt2N_3*sqrt2N_7);
			_TC L2 = 1 / (4 * (2 * N + 5)*sqrt2N_3*sqrt2N_7) + 1 / _TC(2 * (2 * N + 5)*(2 * N + 9)) + 1 / (4 * (2 * N + 9)*sqrt2N_7*sqrt2N_11);

			using std::max;

			return sqrt(_TC(max(L1.upper(), L2.upper())));

		}

		template <typename _TC, class _TMM > typename std::enable_if<std::is_constructible< _TC, _TMM >::value, _TC >::type weighted_Ritz_projection_error(_TMM w, const int NN = -1) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////// Legendre_Bases_Generator :: weighted_Ritz_projection_error()/////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << ">> sqrt(|| Nabla( u - Rh u ) ||_L2^2 + || u ||_L2^2 ) <= CNw || Laplace u + w u ||_L2" << std::endl;
#endif

			if (w < 0) {
				vcp::throw_error<vcp::invalid_argument>(
					">> ERROR : weighted_Ritz_projection_error : Please set weight >= 0 : weight < ", w);
			}

			_TC CN = Ritz_projection_error< _TC >(NN);
			using std::sqrt;

			return CN*sqrt(1 + _TC(w)*pow(CN, 2));
		}

		template <class _TC > _TC Poincare_constant() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////// Legendre_Bases_Generator :: Poincare_constant()/////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << ">> || u ||_L2 <= Cp || Nabla u ||_L2" << std::endl;
#endif
			using std::sqrt;
			return 1/(sqrt(_TC(this->dimension))*kv::constants< _TC >::pi());
		}

		template <typename _TC, class _TMM > typename std::enable_if<std::is_constructible< _TC, _TMM >::value, _TC >::type weighted_Poincare_constant( _TMM w ) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////// Legendre_Bases_Generator :: weighted_Poincare_constant()/////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << ">> || u ||_L2 <= Cwp sqrt( || Nabla u ||_L2^2 + sigma || u ||_L2^2 ) : weight = " << w << std::endl;
#endif
			if (w < 0) {
				vcp::throw_error<vcp::invalid_argument>(
					">> ERROR : weighted_Poincare_constant : Please set weight >= 0 : weight < ", w);
			}
			using std::pow;
			using std::sqrt;
			_TC lambda = this->dimension * pow(kv::constants< _TC >::pi(), 2) + _TC(w);

			return 1 / sqrt(lambda);
		}
		
		template <typename _TC, class _TMM > typename std::enable_if<std::is_constructible< _TC, _TMM >::value &&  vcp::is_interval< _TC >::value, _TC >::type Sobolev_constant( _TMM p ) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////// Legendre_Bases_Generator :: Sobolev_constant()//////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << ">> || u ||_Lp <= Csp || Nabla u ||_L2 : p = " << p << std::endl;
#endif
			using std::pow;

			// p = nq/(n - q), and n = 2 => q = np/(n + p)
			_TC n = _TC(this->dimension);

			if (this->dimension == 1) {
				if (p <= 0) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a Sobolev's constant : Please set p > 0  : p <= ", p);
				}
				return pow(2 / (_TC(p) + 2), 1/_TC(p));
			}
			else if (this->dimension == 2) {
				if (p <= 2) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a Sobolev's constant : Please set p > 2 : p <= ", p);
				}
			}
			else if (this->dimension > 2) {
				_TC	cc = n / (n - 1);
				if (p <= cc.upper()) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a Sobolev's constant : Please set p > 2 : p <= ", p);
				}

				cc = (2 * n) / (n - 2);
				if (p > cc.lower()) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a Sobolev's constant : Please set p > 2 : p <= ", p);
				}
			}
			else {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : Sobolev_constant :  Dimension Error: ", this->dimension);
			}
			_TC q = this->dimension*_TC(p) / (this->dimension + _TC(p));
			_TC Tp = pow(kv::constants< _TC >::pi(), -0.5)*pow(n, -1 / q)*pow((q - 1) / (n - q), 1 - 1 / q) * pow((kv::gamma(1 + n / 2) * kv::gamma(n)) / (kv::gamma(n / q) * kv::gamma(1 + n - n / q)), 1 / n);

			return Tp;
		}

		template <typename _TC, class _TMM, class _TMM2 > typename std::enable_if<std::is_constructible< _TC, _TMM >::value &&  vcp::is_interval< _TC >::value && std::is_constructible< _TC, _TMM2 >::value, _TC >::type weighted_Sobolev_constant(_TMM p, _TMM2 w) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////// Legendre_Bases_Generator :: weighted_Sobolev_constant()/////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << ">> || u ||_Lp <= Csp sqrt( || Nabla u ||_L2^2 + sigma || u ||_L2^2 ) : p = " << p << ", weight = " << w << std::endl;
#endif

			if (w < 0) {
				vcp::throw_error<vcp::invalid_argument>(
					">> ERROR : weighted_Sobolev_constant : Please set weight >= 0 : weight < ", w);
			}

			using std::pow;

			_TC n = _TC(this->dimension);
			_TC pp = _TC(p);

			if (this->dimension == 1) {
				if (p <= 0) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a Sobolev's constant : Please set p > 0  : p <= ", p);
				}
				return pow(2 / (_TC(p) + 2), 1 / _TC(p));
			}
			else if (this->dimension == 2) {
				if (p <= 2) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a Sobolev's constant : Please set p > 2 : p <= ", p);
				}
				_TC pp_harf = pp / 2;
				if (int(pp_harf.upper()) != int(pp_harf.lower())) {
					vcp::throw_error<vcp::verification_error>(
						">> ERROR : Not decide nu : nu.lower() = ",
						int(pp_harf.lower()), ", nu.upper() = ", int(pp_harf.upper()));
				}
				int nu = pp_harf.lower();
				_TC s;
				if (nu == 1) {
					s = _TC(1);
				}
				else {
					s = pp_harf;
					for (int i = 1; i <= nu - 2; i++) {
						s *= (pp_harf - i);
					}
					s = pow(s, 2 / pp);
				}

				_TC lambda = this->dimension * pow(kv::constants< _TC >::pi(), 2);
				return pow(_TC(0.5), 0.5 + (2 * nu - 3) / pp) * s / (pow(lambda + pp / 2 * _TC(w), 1 / pp));
			}
			else if (this->dimension > 2) {				
				if (p <= 2) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not calculate a weighted Sobolev's constant : Please set p > 2 : p <= ", p);
				}

				_TC	cc = (2 * n) / (n - 2);
				if (p > cc.lower()) {
					vcp::throw_error<vcp::domain_error>(
						">> ERROR : Not exist a weighted Sobolev's constant : Please set p > 2 : p <= ", p);
				}
				using std::sqrt;
				_TC lambda = this->dimension * pow(kv::constants< _TC >::pi(), 2);
				
				_TC s = n * (1/pp - 0.5 + 1 / n);
				_TC ss = _TC(0);
				if (s.lower() < ss.lower()) {	
					s.lower() = ss.lower();
				}

				return pow((n - 1) / (sqrt(n) * (n - 2)), 1 - s) * pow(s / ( s * lambda + w ), s / 2);
			}
			else {
				vcp::throw_error<vcp::dimension_error>(">> ERROR : Sobolev_constant :  Dimension Error: ", this->dimension);
			}
		}


		matrix< int > output_list() {
			return this->list;
		}

		void clear() {
			this->Order_of_Base = 0;
			this->elementsize = 0;
			this->variablesize = 0;
			this->dimension = 0; 
			this->mode = 0;
			this->p = 0;
			this->Order_uh = 0;
			this->flag_order_uh = false;

			this->phi.clear();
			this->dphi.clear();
			this->ddphi.clear();

			this->list.clear();
			this->Pointlist.clear();
			this->Legendre_with_GLpoint.clear();
			this->DifDifLegendre_with_GLpoint.clear();
			this->uh.clear();
			this->phi_point.clear();
			this->uhphi_point.clear();
			this->weight_point.clear();
		}

	};
}
#endif // VCP_LDBASE_HPP
