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

int main(void){

	std::cout.precision(17);
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/*************                                                       **************/" << std::endl;
	std::cout << "/************* Verified computation for solution to Emden's equation **************/" << std::endl;
	std::cout << "/*************                                                       **************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/\n" << std::endl;


	vcp::matrix< AppData, POLICY > uh;
	vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;

	int Order_legendre = 30;
	int uh_Order_legendre = 20;
	int p = 2;
	int Dimension = 2;
	int Number_of_variables = 1;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "Inverse Norm for Legendre Bases Order = " << Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

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
	uh.ones(list_uh.rowsize(), Number_of_variables);
	//	uh(0) = 30;
	uh = 200 * uh;
	uh(0) = 600;

	std::cout << "Newton Method Start " << std::endl;
	{
		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

		vcp::matrix< AppData, POLICY > uh2phi;
		vcp::matrix< AppData, POLICY > uhphiphi;
		vcp::matrix< AppData, POLICY > DF;
		vcp::matrix< AppData, POLICY > F;
		vcp::matrix< AppData, POLICY > syuusei;
		vcp::matrix< AppData, POLICY > check;

		{
			AppData cc;
			while(1){
				Approximate_Generator.setting_uh(uh);

				uh2phi = Approximate_Generator.uhphi(2);
				uhphiphi = Approximate_Generator.uhphiphi(1);

				DF = DL - 2*uhphiphi;
				F = DL * uh - uh2phi;
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

	// uh data for Grafics
	vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
//	std::cout << Grafics << std::endl;

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
/******************* Calculate Inverse Norm || F'[uh]^-1 ||_(H-1,H10) <= K *********************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData CN, CNw, Cp, Cpw, Cs, Cs3, Csw;
	VData fduh_Linf_norm = 2 * VData(uh_max[0]);
	VData K = VData(0);
	{
		std::cout << "\nCalculate Inverse Norm || F'[uh]^-1 ||_(H-1,H10) <= K" << std::endl;
		vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
		Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();
		vcp::matrix< VData, VPOLICY > uhi;
		vcp::convert(uh, uhi);
		// uh setting : Last Argument is list divide : full list => 1 , even list => 2 
		Verification_Generator.setting_uh(uhi, list_uh, 2);
			
		vcp::matrix< VData, VPOLICY > uhphiphi = Verification_Generator.uhphiphi(1);

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();

		// How to calculate some constants 
		Cp = Verification_Generator.Poincare_constant< VData >();
		std::cout << "Cp = " << Cp << std::endl;

		Cs3 = Verification_Generator.Sobolev_constant< VData >(3);
		std::cout << "Cs3 = " << Cs3 << ", p =" << "3" << std::endl;

		CN = Verification_Generator.Ritz_projection_error< VData >();
		std::cout << "CN = " << CN << std::endl;
		/* 
		double weight_value = 1234.567;
		CNw = Verification_Generator.weighted_Ritz_projection_error< VData >(weight_value);
		std::cout << "CNw = " << CNw << ", weight_value = " << weight_value << std::endl;

		Cpw = Verification_Generator.weighted_Poincare_constant< VData >(weight_value);
		std::cout << "Cpw = " << Cpw << ", weight_value = " << weight_value << std::endl;

		
		Cs = Verification_Generator.Sobolev_constant< VData >(4);
		std::cout << "Cs4 = " << Cs << ", p =" << "4" << std::endl;
		Cs = Verification_Generator.Sobolev_constant< VData >(VData("5.2"));
		std::cout << "Cs5.2 = " << Cs << ", p =" << "5.2" << std::endl;
		Cs = Verification_Generator.Sobolev_constant< VData >(6);
		std::cout << "Cs6 = " << Cs << ", p =" << "6" << std::endl;

		Csw = Verification_Generator.weighted_Sobolev_constant< VData >(4, weight_value);
		std::cout << "Csw4 = " << Csw << ", p =" << "4" << ", weight_value = " << weight_value << std::endl;
		Csw = Verification_Generator.weighted_Sobolev_constant< VData >(VData("5.2"), weight_value);
		std::cout << "Csw5.2 = " << Csw << ", p =" << "5.2" << ", weight_value = " << weight_value << std::endl;
		Csw = Verification_Generator.weighted_Sobolev_constant< VData >(6, weight_value);
		std::cout << "Csw6 = " << Csw << ", p =" << "6" << ", weight_value = " << weight_value << std::endl;
		*/
		Verification_Generator.clear();

		vcp::matrix< VData, VPOLICY >  E;
		eigsymge(2*uhphiphi, DL, E);

		E = diag(E);
		
		for (int i = 0; i < length(E); i++){
			using std::pow;
			using std::abs;
			using std::max;
			VData cc;
			cc = E(i);
			E(i) = E(i)/(1 + pow(CN,2)*fduh_Linf_norm*E(i));
			E(i).upper() = cc.upper();
			E(i) = abs(E(i)/( 1 - E(i) ));
			if (K.upper() < E(i).upper()){
				K = E(i);
			}
		}
		if (K.upper() < 1) {
			K = VData(1);
		}
		std::cout << "Inverse Norm || F'[uh]^-1 ||_(H-1,H10) <= " << K << std::endl;
		
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
		VResData uh4 = Verification_Generator.integral_uh(4);
		std::cout << "uh4 = " << uh4 << std::endl;

		VResData LuhLuh = Verification_Generator.integral_LuhLuh(0);
		std::cout << "LuhLuh = " << LuhLuh << std::endl;

		VResData Luh_uh2 = Verification_Generator.integral_Luhuh(0, 2);
		std::cout << "Luh_uh2 = " << Luh_uh2 << std::endl;

		{
			using std::sqrt;
			using std::abs;
			vcp::convert(sqrt(abs(LuhLuh + 2 * Luh_uh2 + uh4)), Res);
			std::cout << "Residual Norm : || Laplace(uh) - f(uh) ||_L2 <= " << Res << std::endl;
		}
		Res = Cp * Res;
		std::cout << "Residual Norm : || F(uh) ||_(H-1) <= " << Res << std::endl;
	}
	vcp::time.toc();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/***** Calculate Lipschitz constant || F'[w1] - F'[w2] ||_(H-1,H10) <= G || w1 - w2 ||_(H10) *****/
	/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData G;
	{
		std::cout << "\nCalculate Lipschitz constant || F'[w1] - F'[w2] ||_(H-1,H10) <= G || w1 - w2 ||_(H10)" << std::endl;
		using std::pow;
		G = 2 * pow(Cs3, 3);
		std::cout << "G = " << G << std::endl;

		/* 
		vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
		std::cout << "Check Start" << std::endl;
		Verification_Generator.setting(uh_Order_legendre, 1, Dimension, Number_of_variables, 1, uh_Order_legendre);
		std::cout << "Check End" << std::endl;
		// Setting the list of Verification_Generator 
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();

		vcp::matrix< VData, VPOLICY > uhi;
		vcp::convert(uh, uhi);

		vcp::matrix< VData, VPOLICY > uh_H10_norm = transpose(uhi)*(DL*uhi);
		vcp::matrix< VData, VPOLICY > uh_L2_norm = transpose(uhi)*(L*uhi);

		std::cout << "|| uh ||_(H10) = " << uh_H10_norm << std::endl;
		std::cout << "|| uh ||_(L2) = " << uh_L2_norm << std::endl;
		*/
	}
	vcp::time.toc();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/*************************** Calculate Kantorovich's Theorem :  ********************************/
	/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData Kantorovich_constant;
	VData error_bound;
	{
		std::cout << "\nCalculate Kantorovich's Theorem" << std::endl;

		VData alpha = K*Res;
		VData omega = K*G;
		Kantorovich_constant = alpha*omega;

		std::cout << "Check :: K^2 delta G <= 1/2" << std::endl;
		std::cout << "K^2 delta G = " << Kantorovich_constant << std::endl;

		if (Kantorovich_constant.upper() <= 0.5) {
			std::cout << "Verification is success!!" << std::endl;
		}
		else {
			std::cout << "Verification is fail..." << std::endl;
			exit(1);
		}
		using std::sqrt;
		error_bound = (1 - sqrt(1 - 2 * alpha*omega)) / omega;
		std::cout << "|| u^* - uh ||_(H10) <= " << error_bound << std::endl;
	}
	vcp::time.toc();

	return 0;
}
