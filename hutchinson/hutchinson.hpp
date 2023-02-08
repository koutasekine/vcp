#pragma once

#ifndef VCP_LDBASE_FAST_FIXED_MINMAX_HPP
#define VCP_LDBASE_FAST_FIXED_MINMAX_HPP

#include <iostream>
#include <fstream>
#include <cmath>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>


#include <vcp/fourier_basis.hpp>
#include <vcp/vcp_timer.hpp>

#include <vcp/newton.hpp>
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

typedef kv::dd ResData;
typedef vcp::mats< kv::dd > ResPOLICY;

typedef kv::interval< kv::dd > ResVData;
typedef vcp::imats< kv::dd > ResVPOLICY;

namespace vcp {
	template < typename _T, typename _P >
	struct hutchinson : public vcp::Newton< _T, _P > {
		_T alpha;
		_T tau;
		_T beta;
		_T K;
		_T omega;
		_T omega_n;

		int order, vec_size;

		vcp::fourier_series< _T > x;
		vcp::fourier_series< _T > Bsin_nt;
		vcp::fourier_basis< _T, _P > Generator;

		void setting_newton( vcp::matrix< _T, _P >& zh ) override {
			//x = Generator.omit_vec_to_fourier_series( zh.submatrix( {0, vec_size - 1}, {0} ) );
			this->x = Generator.omit_vec_to_fourier_series( zh );
			Generator.clear_data();
		}
		
		vcp::fourier_series< _T > func(){
			return x.diff() - alpha * (1. + Bsin_nt) * x * (1. - x.delay(-omega_n*tau) / K) / omega_n;
		}

		
		vcp::matrix< _T, _P > f() override {

			Generator.add_fourier_fx( this->func() );
			return Generator.output_fx();
		}

		vcp::matrix< _T, _P > Df() override {
			Generator.clear_data();

			Generator.add_dpt();
			// Generator.add_scalar_pt( -alpha /omega_n);
			Generator.add_fourier_pt( - alpha * (1. + Bsin_nt) * (1. - x.delay(-omega_n*tau) / K) / omega_n );
			Generator.add_fourier_pt_delay( (alpha * (1 + Bsin_nt) * x)/K/omega_n, -omega_n*tau );
			
			return Generator.output_Jacobi();
		}

		template < typename _TM >
		void set_parameter( _TM Alpha, _TM Tau, _TM Beta, _TM Ka, _TM Omega, int N){
			alpha = _T( Alpha );
			tau   = _T( Tau );
			beta  = _T( Beta );
			K = _T( Ka );
			omega = _T( Omega );
			omega_n = _T( Omega )/_T(N);

			Bsin_nt.zeros(N+1);
			Bsin_nt.set_sinm( _T( Beta ), N);

		}

		void set_order( const int& Order ){
			order = Order;
			this->x.zeros(order);
			this->Generator.setting_m(order);
			this->Generator.create_full_list();
			vec_size = Generator.list_size();
		}

		vcp::matrix< AppData, AppPolicy > FourierCoefficient(vcp::matrix< AppData, AppPolicy > u, vcp::matrix< AppData, AppPolicy > t, double tau, double T, int n){

		  using std::cos;
		  using std::sin;

		  double h = 1./4096;
		  int N = 10*1024;
		  int L = tau*4096;
		  int l = 1000*N;
		  int ll = l/2;
		  double pi = kv::constants<double>::pi();
		  double w1 = 2*pi*n/T;
		  int nn = w1/h+1;

		  vcp::matrix< AppData, AppPolicy > x, s;
		  x.zeros(nn+1, 1);
		  s.zeros(nn+1, 1);
			double a = 0;
		
		  for(int i=0; i<nn+1; i++){
		    x(i) = u(i+ll);
		    s(i) = 2.*pi*t(i+ll)/w1;
				a += x(i);
			}
		  a /= nn;

		  vcp::matrix< AppData, AppPolicy > f;
		  f.zeros(vec_size, 1);
		  f(0) = 2.*a;
		  double a1 = 0;
		  double b1 = 0;
		  for(int i=1; i<=(vec_size-1)/2; i++){
		    int j = 2*i;
		    a1 = 0;
		    b1 = 0;
		    for(int k=0; k<nn+1; k++){
		        a1 += x(k)*sin(i*s(k));
		        b1 += x(k)*cos(i*s(k));
		    }
		    a1 *= 2.;
		    a1 /= (double)nn;
		    b1 *= 2.;
		    b1 /= (double)nn;
		    f(j-1) = a1;
		    f(j) = b1;
		  }
		  return f;
		}
	
		//オイラー法
		vcp::matrix< AppData, AppPolicy > euler_vdp(double alpha, double tau, double beta, double K, double omega, int vec_size, int n){
			
			vcp::matrix< AppData, AppPolicy > v;

			std::ofstream result_file("euler.csv");
		
		  //乱数生成
		  std::random_device rnd;
		  std::mt19937_64 mt(rnd());
		  std::uniform_real_distribution<> rand100(0.0, 1.0);

		  using std::sin;

		  double T = omega;//T=1
		  double h = 1./4096;
		  int N = 10*1024;
		  int L = tau*4096;
		
		  v.zeros(1000*N, 1);

		  for(int i = 0;i <= 1000*N;i++){
			  v(i) =  4.*rand100(mt);
		  }
		 
		  double t,z,u;

		  u = 4.*rand100(mt);
			//
		    for(int i=1;i<= (998*N-1); i++){
		      t=i*h;
				if(i-L<=0){
					z = 4.*rand100(mt);
				}else{
					z = v(i-1-L);
				}
				
				v(i) = u + (alpha * (1. + beta * sin(t*T)) * u * (1. - z/K))*h;
 				
				u = v(i);
			}

			vcp::matrix< AppData, AppPolicy > t_;
			t_.zeros(1000*N, 1);
				for(int i=0; i<1000*N; i++){
				t_(i) = (i+1)*h;
			}

			int l = 1000*N;
			int ll = l/2;

			for(int i=1; i<=102000*6; i++){
				double uu = (v(ll+i)-v(ll+i-1))/h;
				result_file << v(ll+i-1) << "," << uu << std::endl;
			}

			return FourierCoefficient(v, t_, tau, T, n);
		}

		vcp::matrix< AppData, AppPolicy > initial_fc(int n){
			vcp::matrix< AppData, AppPolicy > F;
			F = euler_vdp(alpha, tau, beta, K, omega, vec_size, n);
			
			//vcp::matrix< AppData, AppPolicy > F;
		    // F = euler_vdp(alpha1, aplha2, gamma, tau, beta, omega*n, vec_size, n);

			// F = runge_kutta(alpha, gamma, tau, beta, omega*n, vec_size, n);
			// F = hein(alpha, gamma, tau, beta, omega*n, vec_size, n);
			return F;
		}

		//void disp_continue() override {}
	};

	//漸近優対角理論
	template < typename _T, typename _P >
	_T bdab( const vcp::matrix< _T, _P >& G, const int m_app, const int m_verify, const int var_num = 1 ){
		int sm = (m_app*4 + 1)*var_num;  //バンド幅
		int mm = (m_verify*2 + 1)*var_num;
		vcp::matrix< _T, _P > A, B, CDf, D;
		A = G.submatrix( {0, sm-1},  {0, sm-1}  );
		B = G.submatrix( {0, sm-1},  {sm, mm-1} );
		CDf = G.submatrix( {sm, mm-1}, {0, sm-1}  );
		D = G.submatrix( {sm, mm-1}, {sm, mm-1} );

		vcp::matrix< _T, _P > Dd, Df;
		Df = D;		
		Dd.zeros( mm-sm, mm-sm );
		Dd(0, 0) = D(0, 0);
		Dd(0, 1) = D(0, 1);

		Df(0, 0) = _T(0);
		Df(0, 1) = _T(0);

		for (int i = 1; i < mm-sm-1; i++){
			Dd(i, i-1) = D(i, i-1);
			Dd(i, i)   = D(i, i);
			Dd(i, i+1) = D(i, i+1);

			Df(i, i-1) = _T(0);
			Df(i, i) = _T(0);
			Df(i, i+1) = _T(0);     
		}

		Dd(mm-sm-1, mm-sm-2) = D(mm-sm-1, mm-sm-2);
		Dd(mm-sm-1, mm-sm-1) = D(mm-sm-1, mm-sm-1);
		
		Df(mm-sm-1, mm-sm-2) = _T(0);
		Df(mm-sm-1, mm-sm-1) = _T(0);
//		Df = D - Dd;

		vcp::matrix< _T, _P > invAB = abs(lss(A, B));
		_T kk1 = norminf( invAB )(0);
		std::cout << "(kk1): || A^{-1}*B ||_{inf} <= " << kk1 << std::endl;

		CDf = horzcat( CDf, Df );

		vcp::matrix< _T, _P > invDdCDf = lss( Dd, CDf );
		_T kk2 = norminf( invDdCDf )(0);
		std::cout << "(kk2): || Dd^{-1}*[C, Df] ||_{inf} <= " << kk2 << std::endl;

		_T kk = max(kk1,kk2);
		std::cout << "(kk): max( || A^{-1}*B ||_{inf}, || Dd^{-1}*[C, Df] ||_{inf}) <= " << kk << std::endl;

		if ( kk.upper() >= 1 ){
			std::cout << "kk >= 1..." << std::endl;
			exit(0);
		}
		else{
			std::cout << "kk < 1 !!" << std::endl;
		}
		_T norm_invA, norm_invDd;
		{
			vcp::matrix< _T, _P > invTMP, I;
			I.eye(A.rowsize());
			invTMP = lss(A, I);
			norm_invA = norminf(invTMP)(0);

			I.clear();
			I.eye(D.rowsize());
			invTMP = lss(Dd, I);
			norm_invDd = norminf(invTMP)(0);
		}
		_T ninvM2 = norm_invA/(1 - kk);
		std::cout << "ninvM2 <= " << ninvM2 << std::endl;

		_T ninvM3 = norm_invDd/(1 - kk);
		std::cout << "ninvM3 <= " << ninvM3 << std::endl;

		_T Mn = max(ninvM2, ninvM3);
		return Mn;
	}

	std::vector< int > change_list( const int m_verify, const int var_num ){
		int mm1 = m_verify*2 + 1; //変数が2つあるから2倍している
		int mm = mm1*var_num;

		std::vector< int > list;

		list.push_back( 0 );
		list.push_back( mm1 );
		for (int k = 1; k < mm1; k += 2){
			list.push_back( k );
			list.push_back( k+1 );
			list.push_back( mm1 + k );
			list.push_back( mm1 + k + 1 );
		}
		return list;
	}

	template < typename _T, typename _P >
	vcp::matrix< _T, _P > change_matrix( const vcp::matrix< _T, _P >& A, const std::vector< int >& clist ){
		vcp::matrix< _T, _P > B;
		B.zeros( A.rowsize(), A.columnsize() );

		std::cout << "A size: " << A.rowsize() << ", " << A.columnsize() << std::endl;
		std::cout << "clist size: " << clist.size() << std::endl;
		
		for(int i = 0; i < clist.size(); i++ ){
			int ii = clist[i];
			for(int j = 0; j < clist.size(); j++ ){
				int jj = clist[j];
				B(i, j) = A(ii, jj);
			}
		}
		
		return B;
	}
}

#endif