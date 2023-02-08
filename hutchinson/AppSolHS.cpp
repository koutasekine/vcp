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

#include "hutchinson.hpp"

/* Approximate data type (Newton method) */
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

typedef kv::dd ResData;
typedef vcp::mats< kv::dd > ResPOLICY;

typedef kv::interval< kv::dd > ResVData;
typedef vcp::imats< kv::dd > ResVPOLICY;


int main(void){

  std::cout.precision(17);

  int m_app = 150;
  vcp::hutchinson< AppData, AppPolicy > HS;
  vcp::matrix< AppData, AppPolicy > xh;
  vcp::fourier_series< AppData > x;

  std::string alpha = "0.5";
  std::string beta = "0.6";
  std::string K = "2.0";
  std::string omega = "1.0";
  std::string tau = "3.0";

  int n = 2;

  double pi = kv::constants< double >::pi();

  std::cout << "Approximate Fourier order m_app = " << m_app << std::endl;
  std::cout << "alpha = " << alpha << std::endl;
  std::cout << "tau = " << tau << std::endl;
  std::cout << "beta = " << beta << std::endl;
  std::cout << "K = " << K << std::endl;
  std::cout << "omega = " << omega << std::endl;
  std::cout << "tilde{omega} = " << std::stod(omega)/n << std::endl;
  std::cout << "n = " << n << std::endl;

  HS.set_order(m_app);
  HS.set_parameter(
    std::stod(alpha),
    std::stod(tau),
    std::stod(beta),
    std::stod(K),
    std::stod(omega),
    n
  );

  xh.zeros(HS.vec_size, 1);
  
  xh = HS.initial_fc(n);

  HS.setting_newton( xh );

  std::ofstream beforensolv("beforensolv.csv");
  for(double w=0.; w<=n*pi; w+=0.001){
    beforensolv << n*w << "," << HS.x.value(w) << "," << (HS.x.diff()/n).value(w) << "\n";
  }
  beforensolv.close();

  HS.setting_newton_tol( 100 );
  xh = HS.solve_nls( xh );
  
  HS.setting_newton( xh );
  std::ofstream appsol("appsol.csv"); //IBD1.cpp,IBD2.cppのために係数を記録しておく
  //　係数はsin⇨cosの順番 fourier_series.hpp 要参照
  for(int i=0;i<xh.rowsize();i++){
      appsol << xh(i) << std::endl;
  }
  appsol << std::stod(tau) << std::endl; //最後にtauも保存しておく
  appsol.close();
  
  x = HS.x;

  std::ofstream afternsolv("afternsolv.csv");
  for(double w=0.; w<=n*pi; w+=0.001){
    afternsolv << n*w << "," << HS.x.value(w) << "," << (HS.x.diff()/n).value(w) << "\n";
  }
  afternsolv.close();

  /*
  std::cout << "\n========================================================" << std::endl;
	 std::cout << "Approximation of x: " << std::endl;
	 std::cout << "||x||_L2 = " << std::endl;
	 std::cout << x.L2norm() << std::endl;

	 std::cout << "x = " << std::endl;
	 std::cout << x << std::endl;
  */
}
