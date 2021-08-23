#pragma once

#ifndef VCP_LDBASE_FAST_FIXED_MINMAX_HPP
#define VCP_LDBASE_FAST_FIXED_MINMAX_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

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
#include <vcp/fourier_quadrature.hpp>

#include <vcp/vcp_timer.hpp>

#include <vcp/newton.hpp>
//いろんなものをinclude
/////////////////////////////////////////
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

typedef kv::dd ResData;
typedef vcp::mats< kv::dd > ResPOLICY;

typedef kv::interval< kv::dd > ResVData;
typedef vcp::imats< kv::dd > ResVPOLICY;

//型の宣言
/////////////////////////////////////////////
namespace vcp {
	template < typename _T, typename _P >
	struct MG : public vcp::Newton< _T, _P > {
		_T b;
		_T c;
		_T n;
		_T tau;
		_T omega;
		_T omega_n;
		_T beta;
		int order, vec_size;
		int Eular_N;
		
		vcp::fourier_series< _T > x, y;
		vcp::matrix< _T, _P > xhvec;
		vcp::fourier_series< _T > Bcost;
		vcp::fourier_basis< _T, _P > Generator1, Generator2;
		vcp::fourier_quadrature< _T, _P > Integration;


	//変数の宣言
/////////////////////////////////////////////////////////////
		void setting_newton( vcp::matrix< _T, _P >& zh ) override {
			this->x = Generator1.omit_vec_to_fourier_series( zh );
			this->xhvec = zh;
			this->Generator1.clear_data();
		}

////////////////////////////////////////////////////////////////////////////////
		vcp::matrix< _T, _P > f() override {
			Generator1.add_fourier_fx( x.diff() + c*x/omega_n - Bcost/omega_n ); // 
			vcp::matrix< _T, _P > zh = Generator1.output_fx();

			using std::pow;
			this->Integration.set_bks( this->order*this->n/omega_n*64 );
			vcp::fourier_series< _T > x_delay = x.delay(-omega_n*tau);
			this->Integration.set_func( [&]( const auto& t ){
				_T z = x_delay.value(t);
				return -((this->b*z)/( _T(1) + pow(z, this->n)))/this->omega_n;
			} );
			zh += this->Integration.output_vector( this->order );

			return zh;
		}
		
////////////////////////////////////////////////////////////////////////
		vcp::matrix< _T, _P > Df() override {
			Generator1.add_scalar_pt( c/omega_n );
			Generator1.add_dpt();
			vcp::matrix< _T, _P > Jmat = Generator1.output_Jacobi();

			using std::pow;
			this->Integration.set_bks( this->order*this->n/omega_n*16 );
			vcp::fourier_series< _T > x_delay = x.delay(-omega_n*tau);
			this->Integration.set_func( [&]( const auto& t ){
				//_T z = this->x.value(t-omega_n*tau);
				_T z = x_delay.value(t);
				z = pow(z, this->n);
				return -this->b/this->omega_n*(( _T(1) - (this->n - _T(1))*z )/pow( _T(1) + z, 2));
			} );

			Jmat += this->Integration.output_delay_Jacobi( this->order, -omega_n*tau );

			return Jmat;			
		}
//////////////////////////////////////////////////////////////////////////////////////////////////
		template < typename _TM >
		void set_parameter( _TM bb, _TM cc, _TM Tau, _TM Beta,  _TM nmm, _TM Omega, int N){
			this->b = _T( bb );
			this->c = _T( cc );
			this->tau   = _T( Tau );
			this->n  = _T( nmm );
			this->omega = _T( Omega );
			this->omega_n = _T( Omega )/_T(N);
			this->beta=_T(Beta);
			
			this->Bcost.zeros(N+1);
			this->Bcost.set_cosm( _T( Beta ), N);
			// this->Eular_N = 4096;
			this->Eular_N = std::pow(2,12);
		}

		void set_order( const int& Order ){
			order = Order;
			this->x.zeros(order);
		
			this->Generator1.setting_m(order);
			this->Generator1.create_full_list();

			this->Generator2.setting_m(order);
			this->Generator2.create_full_list();
			
			this->vec_size = Generator1.list_size();
		}

		vcp::matrix< AppData, AppPolicy > FourierCoefficient(vcp::matrix< AppData, AppPolicy > u, vcp::matrix< AppData, AppPolicy > t, double tau, double T, int nf){
			using std::cos;
			using std::sin;

			double h = 1./this->Eular_N;
			int N = 10*(this->Eular_N/4);
			int L = tau*this->Eular_N;
			int l = 1000*N;
			int ll = l/2;
			double pi = kv::constants<double>::pi();
			double w1 = 2*pi*nf/T;
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
			f.zeros(this->vec_size, 1);
			f(0) = 2*a;
			double a1 = 0;
			double b1 = 0;
			for(int i=1; i<=(this->vec_size-1)/2; i++){
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
//////////////////////////////////////////////////////////////////////////////
		double ff(double u,double z,double t,double T){
			return -( - this->b*z/( _T(1)+std::pow(z,this->n)) + this->c*u - this->beta*std::cos(t*T));
			//return -(k*(y*y-1)*u+mu*y+gamma*y*y*y-alpha*z-beta*cos(t*T));
		}

		///////////////////////////////////////////////////////////////
		vcp::matrix< AppData, AppPolicy > euler_vdp( int nn ){

			vcp::matrix< AppData, AppPolicy > x;

			std::ofstream result_file("euler3.csv");

			std::random_device rnd;
			std::mt19937_64 mt(rnd());
			std::uniform_real_distribution<> rand100(0.0, 1.0);

			using std::cos;

			double T = omega;//T=1
			double h = 1./this->Eular_N;
			int N = 10*(this->Eular_N/4);
			int L = tau*this->Eular_N;
			x.zeros(1000*N, 1);
			double u = rand100(mt);
			double t,z,y;

			for(int i=1; i<=(998*N-1); i++){
				t=i*h;
				if(i-L<=0){
					z = rand100(mt);
				}
				else{
					z = x(i-1-L);
				}
				u = u + ff(u,z,t,T)*h;
				x(i) = u;
			}

			vcp::matrix< AppData, AppPolicy > t_;
			t_.zeros(1000*N, 1);
			for(int i=0; i<1000*N; i++){
				t_(i) = (i+1)*h;
			}

			int l = 1000*N;
			int ll = l/2;

			for(int i=1; i<=102000*6; i++){
				double uu = (x(ll+i)-x(ll+i-1))/h;
				result_file << x(ll+i-1) << "," << uu << std::endl;
			}

			return FourierCoefficient(x, t_, tau, T, nn);
		}

/////////////////////////////////////////////////////////////////////////////////
		vcp::matrix< AppData, AppPolicy > initial_fc(int nn){
			vcp::matrix< AppData, AppPolicy > F;
		     F = euler_vdp(nn);//nn=分数調波
			//F = runge_kutta(b, c, tau, beta,n, omega*nn, vec_size, nn);
			// F = hein(b, c, tau, beta,n, omega*nn, vec_size, nn);
			return F;
		}

		//void disp_continue() override {}
	};

	template < typename _T, typename _P >
	_T bdab( const vcp::matrix< _T, _P >& G, const int m_app, const int m_verify, const int var_num = 1 ){
		int sm = (m_app*4 + 1)*var_num;
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
            
            //もしかしてこっちが正しい?
			//Dd(i, i) = D(i, i);
			//Dd(i, i+1) = D(i, i+1);
			//Dd(i+1, i) = D(i+1, i);
			//Dd(i+1, i+1)   = D(i+1, i+1);

			//Df(i, i) = _T(0);
			//Df(i, i+1) = _T(0);
			//Df(i+1, i) = _T(0);
			//Df(i+1, i+1) = _T(0);
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
		int mm1 = m_verify*2 + 1;
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
