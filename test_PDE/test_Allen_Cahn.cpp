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

#define EPS 0.08

int main(void){

	std::cout.precision(17);
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/***********                                                          *************/" << std::endl;
	std::cout << "/*********** Verified computation for solution to Allen-Cahn equation *************/" << std::endl;
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
	uh.rand(list_uh.rowsize(), Number_of_variables);
	//	uh(0) = 30;
	uh = 4 * uh;
	uh(0) = 20;

	std::cout << "Newton Method Start " << std::endl;
	{
		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

		vcp::matrix< AppData, POLICY > uhphi;
		vcp::matrix< AppData, POLICY > uh3phi;
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
				uh3phi = Approximate_Generator.uhphi(3);
				uh2phiphi = Approximate_Generator.uhphiphi(2);

				using std::pow;
				DF = DL - ( 1.0/(  pow(EPS, 2) )) * ( L - 3.0*uh2phiphi);
				F = DL * uh - ( 1.0/(  pow(EPS, 2) )) *( uhphi - uh3phi )  ;
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
	VData CN, CNw, Cp, Cpw, Cs, Cs3, Cs4, Csw4;
	VData uh_infsup;
	uh_infsup.lower() = uh_min[0];
	uh_infsup.upper() = uh_max[0];
	VData LA = 1.0/ pow(VData(EPS),2);
	VData fduh_Linf_norm = abs(LA*( VData(1) - 3*pow(uh_infsup,2) ));
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
			
		vcp::matrix< VData, VPOLICY > uh2phiphi = Verification_Generator.uhphiphi(2);
		vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();

		// How to calculate some constants 
		Cp = Verification_Generator.Poincare_constant< VData >();
		std::cout << "Cp = " << Cp << std::endl;

		Cs3 = Verification_Generator.Sobolev_constant< VData >(3);
		std::cout << "Cs3 = " << Cs3 << ", p =" << "3" << std::endl;

		CN = Verification_Generator.Ritz_projection_error< VData >();
		std::cout << "CN = " << CN << std::endl;
		
		double weight_value = fduh_Linf_norm.upper();
		std::cout << "Weight value(sigma) = " << weight_value << std::endl;
		
		CNw = Verification_Generator.weighted_Ritz_projection_error< VData >(weight_value);
		std::cout << "CNw = " << CNw << ", weight_value = " << weight_value << std::endl;

		Cpw = Verification_Generator.weighted_Poincare_constant< VData >(weight_value);
		std::cout << "Cpw = " << Cpw << ", weight_value = " << weight_value << std::endl;

		
		Cs4 = Verification_Generator.Sobolev_constant< VData >(4);
		std::cout << "Cs4 = " << Cs4 << ", p =" << "4" << std::endl;

		Csw4 = Verification_Generator.weighted_Sobolev_constant< VData >(4, weight_value);
		std::cout << "Csw4 = " << Csw4 << ", p =" << "4" << ", weight_value = " << weight_value << std::endl;		

	
		Verification_Generator.clear();

		vcp::matrix< VData, VPOLICY >  E;
		eigsymge(( 1/(  pow(VData(EPS), 2) )) * ( L - 3.0*uh2phiphi) + weight_value*L, DL + weight_value*L, E);

		E = 1/diag(E);
		std::cout << E(0) << std::endl;

		VData K2 = abs( weight_value + LA*( VData(1) - 3* pow(uh_infsup,2) ));

		for (int i = 0; i < length(E); i++){
			using std::pow;
			using std::abs;
			using std::max;
			VData cc;
			cc = E(i);
			E(i) = E(i)/(1 + pow(CNw,2)*VData(K2.upper())*E(i));
			E(i).upper() = cc.upper();
			E(i) = abs(E(i)/( E(i) - 1 ));
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

	return 0;
}
