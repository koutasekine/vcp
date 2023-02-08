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


#include <iostream>
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

#include <vcp/ldbase.hpp>

#include <vcp/vcp_timer.hpp>

//typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 500 > > DataType;

typedef kv::dd ResData;
typedef kv::interval< ResData > VResData;
typedef vcp::imats< ResData > VResPOLICY;

//typedef vcp::pdblas POLICY;
typedef vcp::mats< AppData > POLICY;

typedef vcp::pidblas VPOLICY;

#define LAMBDA 250

int main(void){

	std::cout.precision(17);
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/***********                                                          *************/" << std::endl;
	std::cout << "/*********** Verified computation for solution to Nagumo equation *************/" << std::endl;
	std::cout << "/***********                                                          *************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/\n" << std::endl;


	vcp::matrix< AppData, POLICY > uh;
	vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;

	int Order_legendre = 20;
	int uh_Order_legendre = 20;
	int p = 3;
	int Dimension = 2;
	int Number_of_variables = 1;
	double a = 0.01;//追加
	

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "Inverse Norm for Legendre Bases Order = " << Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

	std::cout << "Lambda = " << LAMBDA << std::endl;
	std::cout << "a = " << a << std::endl;

	// Setting of Approximate_Generator
	std::cout << "\nSetting the Generator by Approximate mode " << std::endl;
	vcp::time.tic();
	Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 5);
	// Setting the list of Approximate_Generator
	//Approximate_Generator.setting_list();
	Approximate_Generator.setting_evenlist();

	// output the list => list_uh
	vcp::matrix< int > list_uh = Approximate_Generator.output_list();

	// setting initialization value of uh
	//uh.ones(list_uh.rowsize(), Number_of_variables);
	//uh = 5*uh;
	//uh(0) = 10;
	//uh(1) = 5;
	vcp::load(uh, "Data_Nagumo/uh_Lambda200_Base20");

	std::cout << "Newton Method Start " << std::endl;
	{
		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

		vcp::matrix< AppData, POLICY > uhphi;
		vcp::matrix< AppData, POLICY > uh2phi;//追加
		vcp::matrix< AppData, POLICY > uh3phi;
		vcp::matrix< AppData, POLICY > uhphiphi;//追加
		vcp::matrix< AppData, POLICY > uh2phiphi;
		vcp::matrix< AppData, POLICY > DF;
		vcp::matrix< AppData, POLICY > F;
		vcp::matrix< AppData, POLICY > syuusei;
		vcp::matrix< AppData, POLICY > check;

		{
			AppData cc;
			while(1){
				Approximate_Generator.setting_uh(uh);
				uhphi = Approximate_Generator.uhphi(1);
				uh2phi = Approximate_Generator.uhphi(2);//追加
				uh3phi = Approximate_Generator.uhphi(3);
				uhphiphi = Approximate_Generator.uhphiphi(1);//追加
				uh2phiphi = Approximate_Generator.uhphiphi(2);

				using std::pow;
				DF = DL - LAMBDA * ( -a*L + 2.0*(1.0+a)*uhphiphi  - 3.0*uh2phiphi);
				F = DL * uh - LAMBDA *( -a*uhphi + (1.0+a)*uh2phi -uh3phi);
				syuusei = lss(DF, F);
				uh = uh - syuusei;
				check = max(abs(syuusei));
				cc = check(0);
				std::cout << cc << std::endl;
				if (cc < pow(2.0,-60)) {
					Approximate_Generator.setting_uh(uh);
					std::cout << "Convergence \n" << std::endl;
					break;
				}
			}
		}
	}

	vcp::time.toc();

	vcp::save(list_uh, "Data_Nagumo/list_Base20");
	vcp::save(uh, "Data_Nagumo/uh_Lambda250_Base20");

	// uh data for Grafics
	vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
	std::cout << Grafics << std::endl;//実験する時にコメント外す
//	std::cout << uh << std::endl;//実験する時にコメント外す

	// minimal and maximum value of approximate solution uh
	vcp::time.tic();
	std::cout << "\nCalculate the maximum and minimum value" << std::endl;
	std::vector< kv::interval< double > > x;
	
	x.resize(Dimension);
	for (int d = 0; d < Dimension; d++) {
		x[d] = kv::interval< double >(0, 0.5);
	}
	std::vector< double > uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
	std::vector< double > uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

	for (int i = 0; i < Number_of_variables; i++) {
		std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
	}

	Approximate_Generator.clear();

	vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/********************************** Eigenvalue of tilde{F} *************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData K = VData(0);
	VData CN, CNw, Cp, Cpw, Cs, Cs3, Cs4, Csw4, CpFtilde;
	VData uh_infsup;
	VData uh_Ftilde_norm;
	uh_infsup.lower() = uh_min[0];
	uh_infsup.upper() = uh_max[0];
	VData VLAMBDA = VData( LAMBDA );
	VData fdtildeuh_Linf_norm = abs(VLAMBDA*( -VData(a) - 3*pow(uh_infsup,2) ));
{
	std::cout << "\nEigenvalue of tilde{F}" << std::endl;
	vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
	Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
	//Verification_Generator.setting_list();
	Verification_Generator.setting_evenlist();
	vcp::matrix< VData, VPOLICY > uhi;
	vcp::convert(uh, uhi);
	// uh setting : Last Argument is list divide : full list => 1 , even list => 2 
	Verification_Generator.setting_uh(uhi, list_uh, 2);

	// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
	vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();
	vcp::matrix< VData, VPOLICY > uhphiphi = Verification_Generator.uhphiphi(1);
	vcp::matrix< VData, VPOLICY > uh2phiphi = Verification_Generator.uhphiphi(2);
	vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();
//	Verification_Generator.clear();

	
	vcp::matrix< VData, VPOLICY > G = DL - VLAMBDA * ( -VData(a)*L - VData(3)*uh2phiphi);
	DL.clear();
	uh_Ftilde_norm = (transpose(uhi)*G*uhi)(0);
	uh2phiphi.clear();

	vcp::matrix< VData, VPOLICY > E;
	eigsymge(L, G, E);
	L.clear();

	E = 1/diag(E);
	std::cout << "Minmum Eigenvalue of tilde{F}:" << std::endl;
	std::cout << min(E) << std::endl;
	
	CN = Verification_Generator.Ritz_projection_error< VData >();
	std::cout << "CN = " << CN << std::endl;

	Cs3 = Verification_Generator.Sobolev_constant< VData >(3);
	std::cout << "Cs3 = " << Cs3 << ", p =" << "3" << std::endl;

	Cs4 = Verification_Generator.Sobolev_constant< VData >(4);
	std::cout << "Cs4 = " << Cs4 << ", p =" << "4" << std::endl;

	VData CFtilde = CN*(1 + 1/min(E)(0)*fdtildeuh_Linf_norm);
	std::cout << "CFtilde = " << CFtilde << std::endl;

	VData lambda_Ftilde = min(E)(0)/(1+pow(CFtilde,2)*min(E)(0));
	std::cout << "lambda_Ftilde = " << lambda_Ftilde << std::endl;

	CpFtilde = 1/sqrt(lambda_Ftilde);
	CpFtilde.lower() = CpFtilde.upper();
	std::cout << "CpFtilde = " << CpFtilde << std::endl;

	eigsymge(uhphiphi, G, E);
	G.clear();
	uhphiphi.clear();

	E = 1/diag(E);

	VData CFtilde2 = sqrt(abs(uh_infsup))*CFtilde;
	CFtilde2.lower() = CFtilde2.upper();
	std::cout << "CFtilde2 = " << CFtilde2 << std::endl;
	for (int i = 0; i < E.rowsize(); i++){
		E(i).lower() = (E(i)/(1+pow(CFtilde2,2)*E(i))).lower();
	}
	
	for (int i = 0; i < E.rowsize(); i++){
		if ( E(i).upper() < (2*(VLAMBDA*(VData(1.0)+a))).lower() ){
			E(i).lower() = E(i).upper();
		}
		else if ( E(i).lower() > (2*(VLAMBDA*(VData(1.0)+a))).upper() ){
			E(i).upper() = E(i).lower();
		}
		else {
			std::cout << "Verification failed..." << std::endl;
		}
	}

	for (int i = 0; i < E.rowsize(); i++ ){
		E(i) = E(i)/(E(i)-2*(VLAMBDA*(VData(1.0)+a) ));
	}
	std::cout << "E = " << E << std::endl;	
	E = min(abs( E ));
	std::cout << "E = " << E << std::endl;	

	K = E(0).upper();
	std::cout << "K = " << K << std::endl;
}

	vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/******************* Calculate Residual Norm || Laplace(uh) - f(uh) ||_L2 **********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData Res = VData(0);
	{
		std::cout << "\nCalculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
		vcp::Legendre_Bases_Generator< DataType, VResData, VResPOLICY > Verification_Generator;
		Verification_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 2);
		// Setting the list of Verification_Generator 
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();
		vcp::matrix< VResData, VResPOLICY > uhi;
		vcp::interval(uh, uhi);
		//vcp::convert(uh, uhi);

		Verification_Generator.setting_uh(uhi);

		// || Laplace(uh) - f(uh) ||_L2 = sqrt( | (Laplace(uh), Laplace(uh))_L2 + 2(-Laplace(uh), f(uh))_L2 + (f(uh), f(uh))_L2 | )
		VResData uh2 = Verification_Generator.integral_uh(2);
		VResData uh3 = Verification_Generator.integral_uh(3);		
		VResData uh4 = Verification_Generator.integral_uh(4);
		VResData uh5 = Verification_Generator.integral_uh(5);
		VResData uh6 = Verification_Generator.integral_uh(6);
		std::cout << "uh2 = " << uh2 << std::endl;
		std::cout << "uh3 = " << uh3 << std::endl;
		std::cout << "uh4 = " << uh4 << std::endl;
		std::cout << "uh5 = " << uh5 << std::endl;
		std::cout << "uh6 = " << uh6 << std::endl;								

		VResData LuhLuh = Verification_Generator.integral_LuhLuh(0);
		std::cout << "LuhLuh = " << LuhLuh << std::endl;

		VResData Luh_uh1 = Verification_Generator.integral_Luhuh(0, 1);
		VResData Luh_uh2 = Verification_Generator.integral_Luhuh(0, 2);
		VResData Luh_uh3 = Verification_Generator.integral_Luhuh(0, 3);		
		std::cout << "Luh_uh1 = " << Luh_uh1 << std::endl;
		std::cout << "Luh_uh2 = " << Luh_uh2 << std::endl;
		std::cout << "Luh_uh3 = " << Luh_uh3 << std::endl;				

		VResData first = LuhLuh;
		VResData second = -VResData(2)*LAMBDA*( -a*Luh_uh1 + (VResData(1) + a)*Luh_uh2 - Luh_uh3 );
		VResData third = LAMBDA*LAMBDA*( VResData(a)*a*uh2 -VResData(2)*a*(VResData(1) + a)*uh3 + ((VResData(1)+a)*(VResData(1)+a) + VResData(2)*a)*uh4 - 2*(VResData(1)+a)*uh5 + uh6);
		std::cout << "first : " << first << std::endl;
		std::cout << "second : " << second << std::endl;		
		std::cout << "third : " << third << std::endl;
		{
			using std::sqrt;
			using std::abs;
			vcp::convert(sqrt(abs(first - second + third)), Res);
			std::cout << "Residual Norm : || Laplace(uh) - f(uh) ||_L2 <= " << Res << std::endl;
		}
		Res = CpFtilde * Res;
		std::cout << "Residual Norm : || F(uh) ||_(H-1) <= " << Res << std::endl;
	}
	vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/******************* G **********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData G = VData(0);
	{
		G = 2*VData(LAMBDA)*(1+VData(a))*pow(Cs3,3) + 3*VData(LAMBDA)*pow(Cs4,4)*(2*uh_Ftilde_norm + 4*K*Res);
		std::cout << "G = " << G << std::endl;
	}


/////////////////////////////////////////////////////////////////////////////////////////////////
/******************* Kantorovich **********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	VData Check = pow(K,2)*Res*G;
	std::cout << "K^2*delta*G = " << Check << std::endl;
	if (Check.upper() <= 0.5 ){
		std::cout << "Verification Succeed!" << std::endl;
		VData rho = (1 - sqrt(1 - 2*Check))/(K*G);
		std::cout << "|| u*-uh || <= " << rho << std::endl;
	}
	else {
		std::cout << "Verofocation failed..." << std::endl;
	}

	return 0;
}
