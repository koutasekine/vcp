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

/* Approximate data type (Newton method) */
//typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

/* Approximate data type (Krawczyk method) */
typedef kv::interval< kv::dd > Kraw_Data;
typedef kv::interval< kv::mpfr< 1500 > > Kraw_DataType;
typedef vcp::imats< kv::dd > Kraw_POLICY;

typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 1500 > > DataType;

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

	vcp::matrix< int > list_uh;
	vcp::matrix< AppData, POLICY > uh;
	vcp::matrix< Kraw_Data, Kraw_POLICY > uhi;
	std::vector< double > uh_min;
	std::vector< double > uh_max;

	int uh_Order = 40;
	int VOrder = 40;
	int p = 2;
	int Dimension = 2;
	int Number_of_variables = 1;
	int mode;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order << std::endl;
	std::cout << "Inverse Norm for Legendre Bases Order = " << VOrder << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;



/////////////////////////////////////////////////////////////////////////////////////////////////
/***************************************** Create uh *******************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	{
		// Setting of Generator
		vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Generator;
		std::cout << "\nSetting the Generator by Approximate mode " << std::endl;
		std::cout << "Newton Method Start " << std::endl;
		vcp::time.tic();
		mode = 5;
		Generator.setting(uh_Order, p, Dimension, Number_of_variables, mode);

		// Setting the list of Generator
		//Generator.setting_list();
		Generator.setting_evenlist();

		// output the list => list_uh
		list_uh = Generator.output_list();
		// std::cout << "list_uh = " << std::endl;
		// std::cout <<  list_uh << std::endl;


		// setting initialization value of uh
		uh.ones(list_uh.rowsize(), Number_of_variables);
		uh(0) = 30;
		uh = 20 * uh;	
	
		vcp::matrix< AppData, POLICY > DF;
		vcp::matrix< Kraw_Data, Kraw_POLICY > kraw_init;
		{
			// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
			vcp::matrix< AppData, POLICY > D = Generator.dphidphi();
			// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
			vcp::matrix< AppData, POLICY > L = Generator.phiphi();

			vcp::matrix< AppData, POLICY > uh2phi;
			vcp::matrix< AppData, POLICY > uhphiphi;
			vcp::matrix< AppData, POLICY > F;
			vcp::matrix< AppData, POLICY > syuusei;
			vcp::matrix< AppData, POLICY > check;

			{
				AppData cc;
				while(1){
					Generator.setting_uh(uh);

					uh2phi = Generator.uhphi(2);
					uhphiphi = Generator.uhphiphi(1);

					DF = D - 2*uhphiphi;
					F = D * uh - uh2phi;
					syuusei = lss(DF, F);
					uh = uh - syuusei;
					check = max(abs(syuusei));
					cc = check(0);
					std::cout << cc << std::endl;
					if (cc < pow(2.0,-60)) {
						DF = inv(DF);
						syuusei = 2*abs(syuusei);
						kraw_init.zeros(syuusei.rowsize(), 1);
						for (int i = 0; i < kraw_init.rowsize(); i++){
							kraw_init(i).upper() = syuusei(i);
							kraw_init(i).lower() = -syuusei(i); 
						}
						kraw_init = 2*kraw_init;
						Generator.setting_uh(uh);
						std::cout << "Convergence \n" << std::endl;
						break;
					}
				}
			}
		}
	
		vcp::time.toc();
		vcp::matrix< AppData, POLICY > Grafics = Generator.output_uh_for_graphics(100);
		std::cout << "Grafics = " << std::endl;
		std::cout << Grafics << std::endl;

		Generator.clear();
		vcp::time.tic();

		std::cout << "Krawczyk Method Start" << std::endl;
		vcp::Legendre_Bases_Generator< Kraw_DataType, Kraw_Data, Kraw_POLICY > Generator2;
		mode = 1;
		Generator2.setting(uh_Order, p, Dimension, Number_of_variables, mode, uh_Order);

		
		vcp::interval(uh, uhi);
		Generator2.setting_evenlist();
		Generator2.setting_uh(uhi, list_uh, 2);
		vcp::matrix< Kraw_Data, Kraw_POLICY > D = Generator2.dphidphi();
		vcp::matrix< Kraw_Data, Kraw_POLICY > uh2phi = Generator2.uhphi(2);
		vcp::matrix< Kraw_Data, Kraw_POLICY > F = D * uhi - uh2phi;
		uh2phi.clear();
		
		vcp::matrix< Kraw_Data, Kraw_POLICY > R;
		vcp::interval(DF, R);
		DF.clear();

		vcp::matrix< Kraw_Data, Kraw_POLICY > Z = -R*F;

		int ite = 0;
		while(true){
			Generator2.setting_uh(uhi + kraw_init, list_uh, 2);
			vcp::matrix< Kraw_Data, Kraw_POLICY > uhphiphi = Generator2.uhphiphi(1);

			vcp::matrix< Kraw_Data, Kraw_POLICY > DFI = D - 2*uhphiphi;
		//	D.clear();
			uhphiphi.clear();

			vcp::matrix< Kraw_Data, Kraw_POLICY > I;
			I.eye(DFI.rowsize());
			Z = Z + (I - R*DFI)*kraw_init;
			bool Flag = true;
			for (int i = 0; i < kraw_init.rowsize(); i++){
				if (!subset( Z(i), kraw_init(i) )){
					// std::cout << "i = " << i << "," << Z(i) << "," << kraw_init(i) << std::endl;
					Flag = false;
				}
			}
			if (Flag) {
				std::cout << "Verification Succeeded!!!!" << std::endl;
				uhi = uhi + kraw_init;
				std::cout << "uhi = " << std::endl;
				std::cout << uhi << std::endl;
				break;
			}
			kraw_init = 4*kraw_init;
			std::cout << ite << std::endl;
			ite++;
		}
		vcp::time.toc();

		// minimal and maximum value of approximate solution uh
		vcp::time.tic();
		std::cout << "\nCalculate the maximum and minimum value" << std::endl;
		Generator2.setting_uh(uhi, list_uh, 2);

		std::vector< kv::interval< double > > x;
		
		x.resize(Dimension);
		for (int d = 0; d < Dimension; d++) {
			x[d] = kv::interval< double >(0, 0.5);
		}

		
		uh_min = Generator2.global_min(x, std::pow(2.0, -9));
		uh_max = Generator2.global_max(x, std::pow(2.0, -9));

		for (int i = 0; i < Number_of_variables; i++) {
			std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
		}
		Generator2.clear();
		vcp::time.toc();
	}

/////////////////////////////////////////////////////////////////////////////////////////////////
/************************** Compute C1, C2 and C3 in Lemma 1 and K *****************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData Ch, Cs2, Cs4;
	VData CfdX = 2 * VData(uh_max[0]);
	VData invS_norm = VData(0);
	VData tau = VData(0);
	VData C1, C2, C3, K;
	VData kappa = VData(0);
	vcp::matrix< VData, VPOLICY > phi_i_norm;
	vcp::matrix< VData, VPOLICY > D;
	vcp::matrix< VData, VPOLICY > DF;
	vcp::matrix< int > list_verification;
	{
		std::cout << "\nCompute C1, C2 and C3 in Lemma 1" << std::endl;
		vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Generator;
		Generator.setting(VOrder, p, Dimension, Number_of_variables, 1, uh_Order);
		//Generator.setting_list();
		Generator.setting_evenlist();
		vcp::matrix< VData, VPOLICY > uhi_rough;
		vcp::convert(uhi, uhi_rough);
		// uh setting : Last Argument is list divide : full list => 1 , even list => 2 
		Generator.setting_uh(uhi_rough, list_uh, 2);

		list_verification = Generator.output_list();
			
		vcp::matrix< VData, VPOLICY > uhphiphi = Generator.uhphiphi(1);

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		D = Generator.dphidphi();
		vcp::matrix< VData, VPOLICY > L = Generator.phiphi();

		// Compute some constants 
		Cs2 = Generator.Poincare_constant< VData >();
		std::cout << "Cs2 = " << Cs2 << std::endl;

		Cs4 = Generator.Sobolev_constant< VData >(4);
		std::cout << "Cs4 = " << Cs4 << ", p =" << "4" << std::endl;

		Ch = Generator.Ritz_projection_error< VData >();
		std::cout << "Ch = " << Ch << std::endl;
		Generator.clear();

		DF = D - 2*(uhphiphi + transpose(uhphiphi)) + uhphiphi*lss(D,transpose(uhphiphi));
		uhphiphi.clear();
		std::cout << "T is invertible!" << std::endl;

		vcp::matrix< VData, VPOLICY >  E;
		compsym(DF);
		compsym(L);
		eigsymge(DF, L, E);
		L.clear();

		E = abs(diag(E));
		tau = sqrt((1/min(E))(0));
		std::cout << "tau = " << tau << std::endl;

		C1 = Ch*tau*CfdX;
		C2 = Ch*Cs2*CfdX;
		C3 = Ch*Ch*CfdX;
		std::cout << "CfdX = " << CfdX << std::endl;
		std::cout << "C1 = " << C1 << std::endl;
		std::cout << "C2 = " << C2 << std::endl;
		std::cout << "C3 = " << C3 << std::endl;

		std::cout << "\nCompute K" << std::endl;
		K = 2*pow(Cs4,2)*sqrt(pow(tau,2) + pow(Ch, 2));
		std::cout << "K = " << K << std::endl;
	}
	vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/************************************** Compute delta ******************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData Res = VData(0);
	VData delta = VData(0);
	{
		std::cout << "\nCalculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
		vcp::Legendre_Bases_Generator< DataType, VResData, VResPOLICY > Generator;
		Generator.setting(uh_Order, p, Dimension, Number_of_variables, 2);
		// Setting the list of Generator 
		//Generator.setting_list();
		Generator.setting_evenlist();
		vcp::matrix< VResData, VResPOLICY > uhi;
		vcp::interval(uh, uhi);
		//vcp::convert(uh, uhi);

		Generator.setting_uh(uhi);

		// || Laplace(uh) - f(uh) ||_L2 = sqrt( | (Laplace(uh), Laplace(uh))_L2 + 2(-Laplace(uh), f(uh))_L2 + (f(uh), f(uh))_L2 | )
		VResData uh4 = Generator.integral_uh(4);
		std::cout << "uh4 = " << uh4 << std::endl;

		VResData LuhLuh = Generator.integral_LuhLuh(0);
		std::cout << "LuhLuh = " << LuhLuh << std::endl;

		VResData Luh_uh2 = Generator.integral_Luhuh(0, 2);
		std::cout << "Luh_uh2 = " << Luh_uh2 << std::endl;

		{
			using std::sqrt;
			using std::abs;
			vcp::convert(sqrt(abs(LuhLuh + 2 * Luh_uh2 + uh4)), Res);
			std::cout << "Residual Norm : || Laplace(uh) - f(uh) ||_L2 <= " << Res << std::endl;
		}
		delta = Ch * Res;
		std::cout << "delta <= " << delta << std::endl;
	}
	vcp::time.toc();


/////////////////////////////////////////////////////////////////////////////////////////////////
/********************************** Method 1 (Corollary 1) *************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	{
		std::cout << "\nMethod 1 (Corollary 1)" << std::endl;
		VData eta = delta;
		VData m = sqrt((pow(C1,2) + pow(C2,2) + pow(C3,2) + sqrt( pow(pow(C1,2) - pow(C2,2) + pow(C3,2), 2) + 4*pow(C2,2)*pow(C3,2) ))/2);
		
		std::cout << "eta = " << eta << std::endl;
		std::cout << "m = " << m << std::endl;
		std::cout << "K = " << K << std::endl;

		VData tmp = (2*K*eta)/pow(1-m,2);
		std::cout << "(2*K*eta)/(1-m)^2 = " << tmp << std::endl;

		if (m.upper() < 1 && tmp.upper() < 1){
			std::cout << "Verification succeed using Method 1!!" << std::endl;
			std::cout << "Error bound : " << (1 - m - sqrt(pow(1-m, 2) - 2*K*eta))/K << std::endl;
		}
	}
	vcp::time.toc();

/////////////////////////////////////////////////////////////////////////////////////////////////
/********************************** Method 2 (Corollary 2) *************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	{
		std::cout << "\nMethod 2 (Corollary 2)" << std::endl;

		VData kappa = C3+C1*C2;
		std::cout << "kappa = " << kappa << std::endl;

		VData eta = sqrt(1+pow(C1,2))/(1-kappa)*delta;
		VData m = VData(0);
		K = K/(1-kappa) * sqrt((1 + pow(C1,2) + pow(C2,2) + pow(1-C3,2) + sqrt( pow(pow(C2,2)+pow(1-C3,2)-pow(1-C1,2),2) + 4*pow(C1*(1-C3)+C2,2) ))/2);

		std::cout << "eta = " << eta << std::endl;
		std::cout << "m = " << m << std::endl;
		std::cout << "K = " << K << std::endl;

		VData tmp = 2*K*eta;
		std::cout << "2*K*eta = " << tmp << std::endl;

		if (kappa.upper() < 1 && tmp.upper() < 1){
			std::cout << "Verification succeed using Method 2!!" << std::endl;
			std::cout << "Error bound : " << (1 - sqrt(1 - 2*K*eta))/K << std::endl;
		}
	}
	vcp::time.toc();


	
	return 0;
}
