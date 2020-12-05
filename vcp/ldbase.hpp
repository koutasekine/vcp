#pragma once

#ifndef VCP_LDBASE_HPP
#define VCP_LDBASE_HPP

#include <iostream>

#include <kv/psa.hpp>
#include <vcp/vcp_converter.hpp>
#include <vector>
#include <vcp/matrix.hpp>

#include<omp.h>

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

namespace vcp {
	namespace LegendreBaseFunctions {
		template <typename _T> void psaTodpsa(const std::vector< kv::psa< _T > >& inp, std::vector< kv::psa< _T > >& outp) {
			int n = inp.size();
			int i, j;

			outp.resize(n);
			for (i = 0; i < n; i++) {
				outp[i].v.resize(n);
				for (j = 0; j < n; j++) {
					outp[i].v[j] = _T(0);
				}
			}
			for (i = 0; i < n; i++) {
				for (j = 0; j < n - 1; j++) {
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
		template <typename _T> void LegendreFunc(const kv::psa< _T >& pol, const _T x, _T& y) {
			int i;
			int n = pol.v.size();

			y = pol.v(n - 1);
			for (i = n - 2; i >= 0; i--) {
				y = y*x + pol.v(i);
			}
		}
		template <typename _T> void LegendrePointFunc(const std::vector< kv::psa< _T > >& phi, const std::vector< _T >& Point, std::vector< std::vector< _T > >& out) {
			int i, j;
			int n = Point.size();
			int m = phi.size();

			out.resize(m);
			for (i = 0; i < m; i++) {
				out[i].resize(n);
				for (j = 0; j < n; j++) {
					LegendreFunc(phi[i], Point[j], out[i][j]);
				}
			}

		}
	}

	namespace GaussLegendreIntegral {
		template <typename _T> void LegendrePol(const int nn, kv::psa< _T >& outp, kv::psa< _T >& outpm1) {
			kv::psa< _T > xp, oim1, oim2, ZZZ;
			int i, n;
			_T iT;
			n = nn + 1;
			xp.v.resize(n);
			oim1.v.resize(n);
			oim2.v.resize(n);
			outp.v.resize(n);
			outpm1.v.resize(n);
			ZZZ.v.resize(n);
			for (i = 0; i < n; i++) {
				xp.v(0) = _T(0);
				oim1.v(0) = _T(0);
				oim2.v(0) = _T(0);
				outp.v(0) = _T(0);
				outpm1.v(0) = _T(0);
				ZZZ.v(0) = _T(0);
			}
			xp.v(1) = _T(1);
			oim1.v(1) = _T(1);
			oim2.v(0) = _T(1);
			outp = 1 / _T(2) * (3 * xp*oim1 - oim2);

			for (i = 3; i < n; i++) {
				iT = _T(i);
				oim2 = oim1;
				oim1 = outp;
				outp = 1 / iT * ((2 * (iT - 1) + 1)*xp*oim1 - (iT - 1)*oim2);
				if (i == n - 2) {
					outpm1 = outp;
				}
			}
		}
		template <typename _T> void LegendreFunc(const kv::psa< _T >& pol, const _T x, _T& y) {
			int i;
			int n = pol.v.size();

			y = pol.v(n - 1);
			for (i = n - 2; i >= 0; i--) {
				y = y*x + pol.v(i);
			}
		}
		template <typename _T> void LegendreDifFunc(const kv::psa< _T >& pol, const _T x, _T& y) {
			int i;
			int n = pol.v.size();
			y = (n - 1)*pol.v(n - 1);
			for (i = n - 2; i>0; i--) {
				y = y*x + i*pol.v(i);
			}
		}
		template <typename _T> void psa_intvalTtoT(const kv::psa< kv::interval< _T > >& pol, kv::psa< _T >& polout) {
			int i;
			int n = pol.v.size();
			polout.v.resize(n);
			for (i = 0; i<n; i++) {
				polout.v(i) = mid(pol.v(i));
			}
		}
		template <typename _T> void AppNewton(const kv::psa< _T >& pol, const _T x, _T& y) {
			_T res, dif;
			y = x;
			for (int i = 0; i < 10; i++) {
				LegendreFunc(pol, y, res);
				LegendreDifFunc(pol, y, dif);
				y = y - res / dif;
			}
		}
		template <typename _T> void Krawczyk1d(const kv::psa< kv::interval< _T > >& pol, const _T x, kv::interval< _T >& y) {
			kv::interval< _T > res, Ky, dif;
			_T difmid, alpha;

			y = x;
			LegendreFunc(pol, y, res);
			LegendreDifFunc(pol, y, dif);
			difmid = mid(dif);
			alpha = mag(res / difmid);
			y = y + kv::interval< _T >(-2 * alpha, 2 * alpha);
			int i = 0;
			while (1) {
				LegendreFunc(pol, kv::interval< _T >(mid(y)), res);
				LegendreDifFunc(pol, y, dif);
				difmid = mid(dif);
				Ky = mid(y) - res / difmid + (1 - dif / difmid)*(y - mid(y));
				if (i == 1) {
					Ky = intersect(Ky, y);
					if (rad(Ky) > 0.9*rad(y)) break;
				}
				else if (proper_subset(Ky, y)) {
					i = 1;
				}
				else {
					std::cout << "Ahhhhhhhhhhh" << std::endl;
					std::cout << "y= " << y << std::endl;
					std::cout << "Ky= " << Ky << std::endl;
				}
				y = Ky;
			}
		}
		template <typename _T> void LPointWeight(const int n, std::vector< kv::interval< _T > >& Point, std::vector< kv::interval< _T > >& Weight) {
			_T x, y, K1, K2;
			_T PI = kv::constants< _T >::pi();
			kv::interval< _T > res, dif;
			kv::psa< _T > psapp;
			kv::psa< kv::interval< _T > > ps, psm1;

			Point.resize(n);
			Weight.resize(n);

			LegendrePol(n, ps, psm1);
			psa_intvalTtoT(ps, psapp);

			K1 = (n - 1) / (8 * n*n*_T(n));
			K2 = _T(4 * n + 2);

			std::cout.precision(32);
			for (int i = 1; i <= n; i++) {
				x = (1 - K1)*cos((4 * i - 1) / K2*PI);
				AppNewton(psapp, x, y);
				Krawczyk1d(ps, y, Point[i - 1]);
			}

			for (int i = 0; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {
					if (overlap(Point[i], Point[j])) {
						std::cout << "Ahhhhhhhhhhhhhhh" << std::endl;
					}
				}
			}

			for (int i = 1; i <= n; i++) {
				LegendreFunc(psm1, Point[i - 1], res);
				LegendreDifFunc(ps, Point[i - 1], dif);
				Weight[i - 1] = 2 / (n * res * dif);
			}
		}
	}

	template <typename _T> class interval_ld_weightpoint {
	public:
		std::vector< _T > Weight;
		std::vector< _T > Point;
		// Order of Legendre integral
		int n;

		interval_ld_weightpoint< _T >() {
		}

		void set(const int nn) {
			n = nn;
			Weight.resize(n);
			Point.reserve(n);
			GaussLegendreIntegral::LPointWeight(n, Point, Weight);
		}

	};

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
		int elementsize;     // (*this).elementsize = std::pow((*this).phi.size(), (*this).dimension);
		int variablesize;    // Number of Variables
		int dimension;       // Dimension of Domein
		int mode;			 // 1: verification mode without residual, 2: verification mode, k > 2: approximation mode k*10
		int p;				 // -Delta u = u^p
		int Order_uh;		 // Order_uh
		bool flag_order_uh;

		std::vector< kv::psa< _T > > phi;
		std::vector< kv::psa< _T > > dphi;
		std::vector< kv::psa< _T > > ddphi;

		matrix< int > list;
		matrix< int > Pointlist;
		matrix< _TM, _PM > Legendre_with_GLpoint;
		matrix< _TM, _PM > DifDifLegendre_with_GLpoint;
		matrix< _TM, _PM > uh;
		matrix< _TM, _PM > phi_point;
		matrix< _TM, _PM > uhphi_point;
		matrix< _TM, _PM > weight_point;

		void list_set() {
			int mmax = (*this).phi.size();
			int dim = (*this).dimension;
			int di;

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: List of bases------------------------------------------" << std::endl;
			std::cout << ">> Size of list is (" << std::pow(mmax, dim) << "," << dim << ") \n" << std::endl;
#endif
			(*this).list.zeros(std::pow(mmax, dim), dim);
			for (int d = 0; d < dim; d++) {
				di = dim - d;
				for (int j = 0; j < std::pow(mmax, d); j++) {
					for (int i = 0; i < mmax; i++) {
						for (int k = 0; k < std::pow(mmax, di - 1); k++) {
							(*this).list(i*std::pow(mmax, di - 1) + j * std::pow(mmax, di) + k, d) = i;
						}
					}
				}
			}
		}
		void even_list_set() {
			int mmax = (*this).phi.size();
			int mnumber = (mmax + 1) / 2;
			int dim = (*this).dimension;
			int di;

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Even List of bases-------------------------------------" << std::endl;
			std::cout << ">> Size of list is (" << std::pow(mnumber, dim) << "," << dim << ") \n" << std::endl;
#endif
			(*this).list.zeros(std::pow(mnumber, dim), dim);
			int ii = 0;
			for (int d = 0; d < dim; d++) {
				di = dim - d;
				for (int j = 0; j < std::pow(mnumber, d); j++) {
					for (int i = 0; i < mmax; i = i + 2) {
						for (int k = 0; k < std::pow(mnumber, di - 1); k++) {
							(*this).list(ii*std::pow(mnumber, di - 1) + j * std::pow(mnumber, di) + k, d) = i;
						}
						ii++;
					}
					ii = 0;
				}
			}
		}
		void Pointlist_set() {
			int mmax = (*this).Legendre_with_GLpoint.columnsize();
			int dim = (*this).dimension;
			int di;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: List of root of Legendre polynomial--------------------" << std::endl;
			std::cout << ">> Size of list is (" << std::pow(mmax, dim) << "," << dim << ") \n" << std::endl;
#endif
			(*this).Pointlist.zeros(std::pow(mmax, dim), dim);
			for (int d = 0; d < dim; d++) {
				di = dim - d;
				for (int j = 0; j < std::pow(mmax, d); j++) {
					for (int i = 0; i < mmax; i++) {
						for (int k = 0; k < std::pow(mmax, di - 1); k++) {
							(*this).Pointlist(i*std::pow(mmax, di - 1) + j * std::pow(mmax, di) + k, d) = i;
						}
					}
				}
			}
		}
		void weight_point_set() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Matrix for Weight_Point--------------------------------" << std::endl;
			std::cout << ">> Matrix size of Legendre Base is (" << 1 << "," << (*this).Pointlist.rowsize() << ") \n" << std::endl;
#endif
			(*this).weight_point.ones(1, (*this).Pointlist.rowsize());
			for (int i = 0; i < (*this).Pointlist.rowsize(); i++) {
				_T a = (*this).Weight[(*this).Pointlist(i, 0)];
				_TM b;
				convert(a, b);
				(*this).weight_point(0, i) = b;
				for (int j = 1; j < (*this).Pointlist.columnsize(); j++) {
					a = (*this).Weight[(*this).Pointlist(i, j)];
					convert(a, b);
					(*this).weight_point(0, i) *= b;
				}
			}
		}
		void phi_point_set() {
			int column_Pointlist = (*this).Pointlist.rowsize();
			int row_Point = (*this).list.rowsize();
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Matrix for Legendre base | Row:List | Column:PoinList--" << std::endl;
			std::cout << ">> Matrix size of Legendre Base is (" << row_Point << "," << column_Pointlist << ")" << std::endl;
#endif
			(*this).phi_point.ones(row_Point, column_Pointlist);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < row_Point; i++) {
				for (int j = 0; j < column_Pointlist; j++) {
					for (int d = 0; d < (*this).list.columnsize(); d++) {
						(*this).phi_point(i, j) *= (*this).Legendre_with_GLpoint((*this).list(i, d), (*this).Pointlist(j, d));
					}
				}
			}
		}
		void uhphi_point_set() {
			int column_Pointlist = (*this).Pointlist.rowsize();
			int row_Point = (*this).list.rowsize();
			(*this).Order_uh = (*this).Order_of_Base;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: uh_i * phi_i with PointList----------------------------" << std::endl;
			std::cout << ">> uh size is (" << (*this).elementsize << "," << (*this).variablesize << ")" << std::endl;
			std::cout << ">> uh_i * phi_i size is (" << (*this).variablesize << "," << column_Pointlist << ")\n" << std::endl;
#endif
			(*this).uhphi_point.zeros((*this).variablesize, column_Pointlist);
			if ((*this).mode == 2) {
				_TM Phi_Time;
				for (int v = 0; v < (*this).variablesize; v++) {
					for (int j = 0; j < (*this).Pointlist.rowsize(); j++) {
						for (int i = 0; i < (*this).list.rowsize(); i++) {
							Phi_Time = (*this).Legendre_with_GLpoint((*this).list(i, 0), (*this).Pointlist(j, 0));
							for (int d = 1; d < (*this).list.columnsize(); d++) {
								Phi_Time *= (*this).Legendre_with_GLpoint((*this).list(i, d), (*this).Pointlist(j, d));
							}
							(*this).uhphi_point(v, j) += uh(i, v) * Phi_Time;
						}
					}
				}
			}
			else {
				for (int v = 0; v < (*this).variablesize; v++) {
					for (int i = 0; i < column_Pointlist; i++) {
						(*this).uhphi_point(v, i) = (*this).uh(0, v) * (*this).phi_point(0, i);
						for (int j = 1; j < row_Point; j++) {
							(*this).uhphi_point(v, i) += (*this).uh(j, v) * (*this).phi_point(j, i);
						}
					}
				}
			}
		}
		void uhphi_point_set(const matrix< int >& list_uh, int& order_uh) {
			(*this).Order_uh = order_uh;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base------------------------------------------" << std::endl;
			std::cout << ">> Order of Legendre base for uh is " << (*this).Order_uh << std::endl;
			std::cout << ">> Order of polynomial is " << (*this).Order_uh + 2 << " (= Order_uh + 2)\n" << std::endl;
#endif
			std::vector< kv::psa< _T > > phi_uh;
			LegendreBaseFunctions::makeLegendreBase((*this).Order_uh + 2, phi_uh);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base <== Point of Gauss Legendre--------------" << std::endl;
			std::cout << ">> Size of Legendre_with_GLpoint is (" << (*this).Order_uh << "," << (*this).Point.size() << ")" << std::endl;
			std::cout << ">> Size of DifDifLegendre_with_GLpoint is (" << (*this).Order_uh << "," << (*this).Point.size() << ")\n" << std::endl;
#endif
			std::vector< std::vector< _T > > LBP2;
			LegendreBaseFunctions::LegendrePointFunc(phi_uh, (*this).Point, LBP2);

			matrix< _TM, _PM > Legendre_uh_with_GLpoint;
			Legendre_uh_with_GLpoint.zeros((*this).Order_uh, (*this).Point.size());

			for (int i = 0; i < (*this).Order_uh; i++) {
				for (int j = 0; j < (*this).Point.size(); j++) {
					convert(LBP2[i + 2][j], Legendre_uh_with_GLpoint(i, j));
//					Legendre_uh_with_GLpoint(i, j) = LBP2[i + 2][j];
				}
			}
			for (int j = 0; j < phi_uh.size() - 2; j++) {
				phi_uh[j] = phi_uh[j + 2];
			}
			phi_uh.resize((*this).Order_uh);

			int column_Pointlist = (*this).Pointlist.rowsize();
			(*this).uhphi_point.zeros((*this).variablesize, column_Pointlist);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: uh_i * phi_i with PointList----------------------------" << std::endl;
			std::cout << ">> uh size is (" << list_uh.rowsize() << "," << (*this).variablesize << ")" << std::endl;
			std::cout << ">> uh_i * phi_i size is (" << (*this).variablesize << "," << column_Pointlist << ")\n" << std::endl;
#endif

			for (int v = 0; v < (*this).variablesize; v++) {
				for (int j = 0; j < (*this).Pointlist.rowsize(); j++) {
					for (int i = 0; i < list_uh.rowsize(); i++) {
						_TM Phi_Time = Legendre_uh_with_GLpoint(list_uh(i, 0), (*this).Pointlist(j, 0));
						for (int d = 1; d < list_uh.columnsize(); d++) {
							Phi_Time *= Legendre_uh_with_GLpoint(list_uh(i, d), (*this).Pointlist(j, d));
						}
						(*this).uhphi_point(v, j) += uh(i, v) * Phi_Time;
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
				(*this).make_uh_p(uh_p_phi_point, args...);
			}
			else if (p < 0) {
				std::cout << ">> ERROR : make_uh_p : input arguments is negative ( p = " << p << " )" << std::endl;
				exit(1);
			}
			else {
				for (int i = 0; i < (*this).Pointlist.rowsize(); i++) {
					uh_p_phi_point(0, i) *= pow((*this).uhphi_point((*this).variablesize - (vari_size + 1), i), p);
				}
				(*this).make_uh_p(uh_p_phi_point, args...);
			}
		}

		template<typename... Args> void disp_args(int& vs) {
		}
		template<typename... Args> void disp_args(int& vs, const int& p, const Args&... args) {
			vs += p;
#ifdef VCP_LEGENDRE_DEBUG
			constexpr std::size_t vari_size = sizeof...(args);
			int kk = (*this).variablesize - vari_size;
			std::cout << "u" << kk << "^" << p << "*";
#endif
			(*this).disp_args(vs, args...);
		}

	public:
		Legendre_Bases_Generator< _T, _TM, _PM >() {
			(*this).flag_order_uh = false;
		}

		void setting(const int mm, const int pp, const int dim, const int vs = 1, const int k = 2, const int uh_o = -1) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////////// Legendre_Bases_Generator :: setting(" << mm << "," << pp << "," << dim << "," << vs << "," << k << "," << uh_o << ") //////////////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			if (mm < 2 || pp < 1 || dim < 1 || vs < 1 || k < 1 || uh_o < -1 || uh_o == 0) {
				std::cout << ">> ERROR : setting : Not Appropriate Input Arguments" << std::endl;
				exit(1);
			}
			if (uh_o > 1 && k != 1) {
				std::cout << ">> ERROR : Order of uh > 1 && Mode isn't With out Residual Mode  " << std::endl;
				exit(1);
			}
			if (uh_o > 1) {
				(*this).flag_order_uh = true;
				(*this).Order_uh = uh_o;
			}
			else {
				(*this).flag_order_uh = false;
				(*this).Order_uh = mm;
			}

			int m = mm + 2;

			(*this).Order_of_Base = mm;
			(*this).p = pp;
			(*this).dimension = dim;
			(*this).mode = k;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Gauss-Legendre Numerical Integration-------------------" << std::endl;
#endif
			int n;
			if (k == 1) {
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> Warning : You can only use (uh^p, phi)_{L^2} and (uh^(p-1) phi, phi)_{L^2}" << std::endl;
				std::cout << ">> Don't calculate the residual (uh^p, uh^p)_{L^2}" << std::endl;
#endif
				if ((*this).flag_order_uh) {
					//  2*nn1-1 > (p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2)   -->  n := ((p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2) + 1)/2 + 1
					int nn1 = ((p - 1)*((*this).Order_uh + 2) + 2 * ((*this).Order_of_Base + 2) + 2) / 2 + 2;
					//  2*nn2-1 > p*(Order_uh + 2) + (Order_of_Base + 2)   -->  n := ((p - 1)*(Order_uh + 2) + 2*(Order_of_Base + 2) + 1)/2 + 1
					int nn2 = (p*((*this).Order_uh + 2) + ((*this).Order_of_Base + 2) + 2) / 2 + 2;
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
					n = (k*((*this).Order_of_Base + 2) * (p + 1) + 1) / 2 + 2;
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
				n = (k*((*this).Order_of_Base + 2) * p + 1) / 2 + 2;
				if (n % 2 != 0) {
					n = n + 1;
				}
#ifdef VCP_LEGENDRE_DEBUG
				std::cout << ">> Order of Gauss Legendre Integration is " << n << " = (" << k << "*(Order_of_Base + 2) * p + 1)/2 + 2 \n" << std::endl;
#endif
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
			(*this).interval_ld_weightpoint< _T >::set(n);

			for (int i = 0; i < n; i++) {
				(*this).Point[i] = ((*this).Point[i] + 1) / 2;
				(*this).Weight[i] = (*this).Weight[i] / 2;
			}
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base------------------------------------------" << std::endl;
			std::cout << ">> Order of Legendre base m is " << mm << std::endl;
			std::cout << ">> Order of polynomial is " << mm + 2 << " (= m + 2)\n" << std::endl;
#endif
			LegendreBaseFunctions::makeLegendreBase(m, (*this).phi);
			LegendreBaseFunctions::psaTodpsa((*this).phi, (*this).dphi);
			LegendreBaseFunctions::psaTodpsa((*this).dphi, (*this).ddphi);

#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "--------------Generating :: Legendre Base <== Point of Gauss Legendre--------------" << std::endl;
			std::cout << ">> Size of Legendre_with_GLpoint is (" << mm << "," << (*this).Point.size() << ")" << std::endl;
			std::cout << ">> Size of DifDifLegendre_with_GLpoint is (" << mm << "," << (*this).Point.size() << ")\n" << std::endl;
#endif
			std::vector< std::vector< _T > > LBP2;
			std::vector< std::vector< _T > > LBDDP2;
			LegendreBaseFunctions::LegendrePointFunc((*this).phi, (*this).Point, LBP2);
			LegendreBaseFunctions::LegendrePointFunc((*this).ddphi, (*this).Point, LBDDP2);

			(*this).Legendre_with_GLpoint.zeros(mm, (*this).Point.size());
			(*this).DifDifLegendre_with_GLpoint.zeros(mm, (*this).Point.size());
			for (int i = 0; i < mm; i++) {
				for (int j = 0; j < (*this).Point.size(); j++) {
					_T a = LBP2[i + 2][j];
					_TM b;
					convert(a,b);
					(*this).Legendre_with_GLpoint(i, j) = b;

					a = LBDDP2[i + 2][j];
					convert(a, b);
					(*this).DifDifLegendre_with_GLpoint(i, j) = b;
				}
			}
			for (int j = 0; j < phi.size() - 2; j++) {
				phi[j] = phi[j + 2];
			}
			phi.resize(mm);

			(*this).Pointlist_set();
			(*this).weight_point_set();
			(*this).variablesize = vs;
		}

		void setting_list() {
			(*this).list_set();
			(*this).elementsize = (*this).list.rowsize();
			if ((*this).mode != 2) {
				(*this).phi_point_set();
			}
		}

		void setting_evenlist() {
			(*this).even_list_set();
			(*this).elementsize = (*this).list.rowsize();
			if ((*this).mode != 2) {
				(*this).phi_point_set();
			}
		}

		void setting_uh() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "///////////////////// Legendre_Bases_Generator :: setting_uh() /////////////////////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			if ((*this).flag_order_uh) {
				std::cout << ">> ERROR : setting_uh : Flag for order of uh is True..." << std::endl;
				std::cout << ">> Please check the Last Argument of the function .setting()" << std::endl;
				exit(1);
			}

			(*this).uh.ones((*this).elementsize, (*this).variablesize);
			(*this).uhphi_point_set();
		}
		void setting_uh(const matrix< _TM, _PM >& uuh) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "////////////// Legendre_Bases_Generator :: setting_uh(vcp::matrix uh) //////////////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			if ((*this).flag_order_uh) {
				std::cout << ">> ERROR : setting_uh : Flag for order of uh is True..." << std::endl;
				std::cout << ">> Please check the Last Argument of the function .setting()" << std::endl;
				exit(1);
			}

			if (uuh.rowsize() != (*this).elementsize || uuh.columnsize() != (*this).variablesize) {
				std::cout << ">> ERROR : setting_uh : no much the matrix size : uh" << std::endl;
				exit(1);
			}
			(*this).uh = uuh;
			(*this).uhphi_point_set();
		}
		void setting_uh(const matrix< _TM, _PM >& uuh, const matrix< int >& list_uh, int ind) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "//// Legendre_Bases_Generator :: setting_uh(matrix<_T> uh, matrix<int> list_uh) ////" << std::endl;
			std::cout << "//////////////////////////////////////////////////////////////////////////////////// \n" << std::endl;
#endif
			matrix< int > uh_o = max(max(list_uh));
			int uh_order = uh_o(0) + ind;

			if (!(*this).flag_order_uh) {
				std::cout << ">> ERROR : setting_uh : Flag for order of uh is False..." << std::endl;
				std::cout << ">> Please check the Last Argument of the function .setting()" << std::endl;
				exit(1);
			}
			if ((*this).Order_uh != uh_order) {
				std::cout << ">> ERROR : setting_uh : no much the uh_order : " << (*this).Order_uh << "!=" << uh_order << std::endl;
				std::cout << ">> Please check the Last Argument of the function .setting() or Order of the Argument uh" << std::endl;
				exit(1);
			}
			
			if (uuh.rowsize() != (std::pow(uh_order/ind, (*this).dimension)) || uuh.columnsize() != (*this).variablesize) {
				std::cout << ">> ERROR : setting_uh : no much the matrix size : uh : " << uuh.rowsize() << "!=" << std::pow(uh_order, (*this).dimension) << "||" << uuh.columnsize() << "!=" << (*this).variablesize << std::endl;
				exit(1);
			}
		
			(*this).uh = uuh;
			(*this).uhphi_point_set(list_uh, uh_order);
		}

		template<typename... Args> matrix< _TM, _PM > uhphiphi(const Args&... args) {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////// Legendre_Bases_Generator :: uhphiphi(Args... )//////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif
			constexpr std::size_t vari_size = sizeof...(args);
			if ((*this).mode == 2) {
				std::cout << ">> ERROR : uhphiphi : This function can not use the Residual mode (Mode = 2)" << std::endl;
				exit(1);
			}

			if (vari_size != (*this).variablesize) {
				std::cout << ">> ERROR : uhphiphi : size of input arguments != variablesize ( " << vari_size << " != " << (*this).variablesize << " )" << std::endl;
				exit(1);
			}

			matrix< _TM, _PM > uh_p_phi_point;
			uh_p_phi_point.ones(1, (*this).Pointlist.rowsize());
			(*this).make_uh_p(uh_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the Matrix { (";
#endif
			(*this).disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "phi_i, phi_j )_{L^2} }_{i,j} \n" << std::endl;
#endif
			if ((*this).p < vs) {
				std::cout << ">> ERROR : uhphiphi : p < vs" << "( p=" << (*this).p << ", vs=" << vs << ")" << std::endl;
				exit(1);
			}

			matrix< _TM, _PM > A;
			A.zeros((*this).elementsize, (*this).elementsize);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif
			for (int i = 0; i < (*this).elementsize; i++) {
				for (int j = i; j < (*this).elementsize; j++) {
					for (int k = 0; k < (*this).Pointlist.rowsize(); k++) {
						A(i, j) += (*this).weight_point(0, k) * uh_p_phi_point(0, k) * (*this).phi_point(i, k) * (*this).phi_point(j, k);
					}
				}
			}
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 1; i < (*this).elementsize; i++) {
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
			if ((*this).mode == 2) {
				std::cout << ">> ERROR : uhphi : This function can not use the Residual mode (Mode = 2)" << std::endl;
				exit(1);
			}

			if (vari_size != (*this).variablesize) {

				std::cout << ">> ERROR : uhphi : size of input arguments != variablesize ( " << vari_size << " != " << (*this).variablesize << " )" << std::endl;
				exit(1);
			}

			matrix< _TM, _PM > uh_p_phi_point;
			uh_p_phi_point.ones(1, (*this).Pointlist.rowsize());
			(*this).make_uh_p(uh_p_phi_point, args...);

			int vs = 0;
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ">> Generate the Vector { (";
#endif
			(*this).disp_args(vs, args...);
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << ", phi_j )_{L^2} }_{i,j} \n" << std::endl;
#endif
			if ((*this).p < vs) {
				std::cout << ">> ERROR : uhphi : p < vs" << "( p=" << (*this).p << ", vs=" << vs << ")" << std::endl;
				exit(1);
			}

			matrix< _TM, _PM > A;
			A.zeros((*this).elementsize, 1);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
			#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < (*this).elementsize; i++) {
				for (int k = 0; k < (*this).Pointlist.rowsize(); k++) {
					A(i, 0) += (*this).weight_point(0, k) * uh_p_phi_point(0, k) * (*this).phi_point(i, k);
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

			if ((*this).mode == 2) {
				std::cout << ">> ERROR : phiphi : This function can not use the Residual mode (Mode = 2)" << std::endl;
				exit(1);
			}
			matrix< _TM, _PM > L;
			L.zeros((*this).elementsize, (*this).elementsize);

			matrix< int > lmax = max(max((*this).list));
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
			for (int i = 0; i < (*this).elementsize; i++) {
				for (int j = i; j < (*this).elementsize; j++) {
					_TM Phi_Time = _TM(1);
					for (int d = 0; d < (*this).list.columnsize(); d++) {
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
			for (int i = 1; i < (*this).elementsize; i++) {
				for (int j = 0; j < i; j++) {
					L(i, j) = L(j, i);
				}
			}
			return L;
		}

		matrix< _TM, _PM > dphidphi() {
#ifdef VCP_LEGENDRE_DEBUG
			std::cout << "\n\n////////////////////////////////////////////////////////////////////////////////////" << std::endl;
			std::cout << "/////////////////////// Legendre_Bases_Generator :: dphidhi()///////////////////////" << std::endl;
			std::cout << "////////////////////////////////////////////////////////////////////////////////////" << std::endl;
#endif

			if ((*this).mode == 2) {
				std::cout << ">> ERROR : dphidphi : This function can not use the Residual mode (Mode = 2)" << std::endl;
				exit(1);
			}
			matrix< _TM, _PM > DL;
			DL.zeros((*this).elementsize, (*this).elementsize);
			matrix< int > lmax = max(max((*this).list));
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
			for (int i = 0; i < (*this).elementsize; i++) {
				for (int j = i; j < (*this).elementsize; j++) {
					for (int dc = 0; dc < (*this).list.columnsize(); dc++) {
						_TM Phi_Time = _TM(1);
						for (int dr = 0; dr < (*this).list.columnsize(); dr++) {
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
			for (int i = 1; i < (*this).elementsize; i++) {
				for (int j = 0; j < i; j++) {
					DL(i, j) = DL(j, i);
				}
			}
			return DL;
		}

		matrix< int > output_list() {
			return (*this).list;
		}

		void clear() {
			(*this).Order_of_Base = 0;
			(*this).elementsize = 0;
			(*this).variablesize = 0;
			(*this).dimension = 0; 
			(*this).mode = 0;
			(*this).p = 0;
			(*this).Order_uh = 0;
			(*this).flag_order_uh = false;

			(*this).phi.clear();
			(*this).dphi.clear();
			(*this).ddphi.clear();

			(*this).list.clear();
			(*this).Pointlist.clear();
			(*this).Legendre_with_GLpoint.clear();
			(*this).DifDifLegendre_with_GLpoint.clear();
			(*this).uh.clear();
			(*this).phi_point.clear();
			(*this).uhphi_point.clear();
			(*this).weight_point.clear();
		}

	};
}
#endif // VCP_LDBASE_HPP