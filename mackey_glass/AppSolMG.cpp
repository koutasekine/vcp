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

#include "MG.hpp"
//色々include
//////////////////////////////////////

/* Approximate data type (Newton method) */
typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

typedef kv::dd ResData;
typedef vcp::mats< kv::dd > ResPOLICY;

typedef kv::interval< kv::dd > ResVData;
typedef vcp::imats< kv::dd > ResVPOLICY;
//型の定義
////////////////////////////////////////////
int main(void){

  std::cout.precision(17);

  int m_app = 100;//フーリエ級数の次数
 ////////////////////////////////////////////
  vcp::MG< AppData, AppPolicy > MaG;
  vcp::matrix< AppData, AppPolicy > xh, yh, zh;
//////////////////////////////////////////////////
  std::string b = "0.25";
  //std::string b = "2";
  //std::string c = "0.5";
  std::string c = "0.1";
  std::string nm = "10";
  std::string tau = "5.5";
  //std::string tau = "1";
  //std::string tau = "7.8125";
  std::string omega = "1";
  std::string beta = "0.1";
  int n = 3;

  double pi = kv::constants< double >::pi();

  std::cout<< "Approximate Fourier order m_app = " << m_app << std::endl;
  std::cout<< "b = " << b << std::endl;
  std::cout<< "c = " << c << std::endl;
  std::cout << "tau = " << tau << std::endl;
 std::cout << "nm = " << nm << std::endl;
 std::cout << "omega = " << omega << std::endl;
 std::cout << "beta = " << beta << std::endl;
  
  std::cout << "tilde{omega} = " << std::stod(omega)/n << std::endl;
  std::cout << "n = " << n << std::endl;

  //表示
////////////////////////////////////////////////////
  MaG.set_order(m_app);
  MaG.set_parameter(
    std::stod(b),
    std::stod(c),
    std::stod(tau),
    std::stod(beta),
    std::stod(nm),
    std::stod(omega),
    n
  );
  /*std::ofstream param("parameter.csv");
  param << alpha << std::endl;
  param << gamma << std::endl;
  param << tau << std::endl;
  param << beta << std::endl;
  param << mu << std::endl;
  param << k << std::endl;
  param << omega << std::endl;
  param << n << std::endl;
  param.close();
  */

  xh.zeros(MaG.vec_size, 1);
  
  xh = MaG.initial_fc(n);

  for ( int i = 0; i < xh.rowsize(); i++){
    if (abs(xh(i)) <= pow(10.0, -5) ) xh(i) = 0;
  }

  std::cout << "xh = \n" << xh << std::endl;

  //std::cout << "yh = \n" << yh << std::endl;

//////////////////////////書き込み
  MaG.setting_newton( xh );
    //std::ofstream beforensolv("beforensolv.csv");
     std::ofstream beforensolv("beforensolv5_5_3.csv");
  for(double w=0.; w<=2*n*pi; w+=0.001){
    //beforensolv << MaG.x.value(w) << "," << (MaG.x.diff()/n).value(w) << "\n";
    beforensolv << MaG.x.value(w) << "," << (MaG.x.diff()/n).value(w) << "\n";
  }
  beforensolv.close();

  MaG.setting_newton_tol( 16 );
  xh = MaG.solve_nls( xh );
  
  MaG.setting_newton( xh );
  std::ofstream appsol("appsol.csv"); //IBD1.cpp,IBD2.cppのために係数を記録しておく
  //　係数はsin⇨cosの順番 fourier_series.hpp 要参照
  for(int i=0;i<xh.rowsize();i++){
      appsol << xh(i) << std::endl;
  }
  xh(0) *= 2.0;
  appsol << std::stod(tau) << std::endl; //最後にtauも保存しておく
  appsol.close();
  

  vcp::fourier_series< AppData > x = MaG.x;
 // vcp::fourier_series< AppData > y = VDP.y;

 //std::ofstream afternsolv("afternsolv.csv");
 std::ofstream afternsolv("afternsolv5_5_3.csv");
  for(double w=0.; w<=2*n*pi; w+=0.001){
    //std::cout << MaG.x.value(w) << "," << (MaG.x.diff()/n).value(w) << "\n";
     //afternsolv << MaG.x.value(w) << "," << (MaG.x.diff()/n).value(w) << "\n";
      afternsolv << MaG.x.value(w) << "," << (MaG.x.diff()/n).value(w) << "\n";
  }
  afternsolv.close();

  std::cout << "\n========================================================" << std::endl;
	 std::cout << "Approximation: " << std::endl;
	 std::cout << "||x||_L2 = " << std::endl;
	 std::cout << x.L2norm() << std::endl;

	 std::cout << "x = " << std::endl;
	 std::cout << x << std::endl;
}
