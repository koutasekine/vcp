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

/* Approximate data (Newton method) */
//typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

/* Approximate data (Krawczyk method) */
typedef kv::interval< kv::dd > Kraw_Data;
typedef kv::interval< kv::mpfr< 500 > > Kraw_DataType;
typedef vcp::imats< kv::dd > Kraw_POLICY;

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

	vcp::matrix< int > list_uh;
	vcp::matrix< AppData, POLICY > uh;
	vcp::matrix< Kraw_Data, Kraw_POLICY > uhi;
	std::vector< double > uh_min;
	std::vector< double > uh_max;

	int Order_legendre = 20;
	int uh_Order_legendre = 20;
	int p = 2;
	int Dimension = 1;
	int Number_of_variables = 1;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "Inverse Norm for Legendre Bases Order = " << Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

	{
		// Setting of Approximate_Generator
		vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
		std::cout << "\nSetting the Generator by Approximate mode " << std::endl;
		std::cout << "Newton Method Start " << std::endl;
		vcp::time.tic();
		Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 5);

		// Setting the list of Approximate_Generator
		//Approximate_Generator.setting_list();
		Approximate_Generator.setting_evenlist();

		// output the list => list_uh
		list_uh = Approximate_Generator.output_list();

		// setting initialization value of uh
		uh.ones(list_uh.rowsize(), Number_of_variables);
		//	uh(0) = 30;
		uh = 200 * uh;
		uh(0) = 600;
	
		vcp::matrix< AppData, POLICY > DF;
		vcp::matrix< Kraw_Data, Kraw_POLICY > kraw_init;
		{
			// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
			vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
			// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
			vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

			vcp::matrix< AppData, POLICY > uh2phi;
			vcp::matrix< AppData, POLICY > uhphiphi;
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
						DF = inv(DF);
						syuusei = 2*abs(syuusei);
						kraw_init.zeros(syuusei.rowsize(), 1);
						for (int i = 0; i < kraw_init.rowsize(); i++){
							kraw_init(i).upper() = syuusei(i);
							kraw_init(i).lower() = -syuusei(i); 
						}
						kraw_init = 2*kraw_init;
						Approximate_Generator.setting_uh(uh);
						std::cout << "Convergence \n" << std::endl;
						break;
					}
				}
			}
		}
	
		vcp::time.toc();
		Approximate_Generator.clear();

		vcp::time.tic();
	
		std::cout << "Krawczyk Method Start" << std::endl;
		vcp::Legendre_Bases_Generator< Kraw_DataType, Kraw_Data, Kraw_POLICY > Verification_Generator;		
		Verification_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);

		
		vcp::interval(uh, uhi);
		Verification_Generator.setting_evenlist();
		Verification_Generator.setting_uh(uhi, list_uh, 2);
		vcp::matrix< Kraw_Data, Kraw_POLICY > DL = Verification_Generator.dphidphi();
		vcp::matrix< Kraw_Data, Kraw_POLICY > uh2phi = Verification_Generator.uhphi(2);
		vcp::matrix< Kraw_Data, Kraw_POLICY > F = DL * uhi - uh2phi;
		uh2phi.clear();
		
		vcp::matrix< Kraw_Data, Kraw_POLICY > R;
		vcp::interval(DF, R);
		DF.clear();

		vcp::matrix< Kraw_Data, Kraw_POLICY > Z = -R*F;

		int ite = 0;
		while(true){
			Verification_Generator.setting_uh(uhi + kraw_init, list_uh, 2);
			vcp::matrix< Kraw_Data, Kraw_POLICY > uhphiphi = Verification_Generator.uhphiphi(1);

			vcp::matrix< Kraw_Data, Kraw_POLICY > DFI = DL - 2*uhphiphi;
		//	DL.clear();
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
		Verification_Generator.setting_uh(uhi, list_uh, 2);

		std::vector< kv::interval< double > > x;
		
		x.resize(Dimension);
		for (int d = 0; d < Dimension; d++) {
			x[d] = kv::interval< double >(0, 0.5);
		}

		
		uh_min = Verification_Generator.global_min(x, std::pow(2.0, -9));
		uh_max = Verification_Generator.global_max(x, std::pow(2.0, -9));

		for (int i = 0; i < Number_of_variables; i++) {
			std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
		}
		Verification_Generator.clear();
		vcp::time.toc();
	}

/////////////////////////////////////////////////////////////////////////////////////////////////
/************************** Calculate Theorem 1 : inv(T) and inv(S) ****************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData CN, CNw, Cp, Cpw, Cs, Cs4, Csw;
	VData fduh_Linf_norm = 2 * VData(uh_max[0]);
	VData invS_norm = VData(0);
	VData rhoh = VData(0);
	VData kappa = VData(0);
	vcp::matrix< VData, VPOLICY > phi_i_norm;
	vcp::matrix< VData, VPOLICY > DL;
	vcp::matrix< VData, VPOLICY > DF;
	vcp::matrix< int > list_verification;


	{
		std::cout << "\nCalculate Theorem 1 : inv(T) and inv(S) <= K" << std::endl;
		vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
		Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();
		vcp::matrix< VData, VPOLICY > uhi_rough;
		vcp::convert(uhi, uhi_rough);
		// uh setting : Last Argument is list divide : full list => 1 , even list => 2 
		Verification_Generator.setting_uh(uhi_rough, list_uh, 2);

		list_verification = Verification_Generator.output_list();
			
		vcp::matrix< VData, VPOLICY > uhphiphi = Verification_Generator.uhphiphi(1);

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		DL = Verification_Generator.dphidphi();
		vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();

		// How to calculate some constants 
		Cp = Verification_Generator.Poincare_constant< VData >();
		std::cout << "Cp = " << Cp << std::endl;

		Cs4 = Verification_Generator.Sobolev_constant< VData >(4);
		std::cout << "Cs4 = " << Cs4 << ", p =" << "4" << std::endl;

		CN = Verification_Generator.Ritz_projection_error< VData >();
		std::cout << "CN = " << CN << std::endl;
		Verification_Generator.clear();

		DF = DL - 2 * uhphiphi;
		uhphiphi.clear();

		inv(DF);
		std::cout << "T is invertible!" << std::endl;

		vcp::matrix< VData, VPOLICY >  E;
		compsym(DF);
		compsym(L);
		eigsymge(DF, L, E);
		phi_i_norm = sqrt(diag(L));
		L.clear();

		E = abs(diag(E));
		rhoh = (1/min(E))(0);
		std::cout << "rhoh = " << 1/min(E) << std::endl;

		kappa = pow(CN,2) * ( 1 + rhoh*fduh_Linf_norm )*fduh_Linf_norm;
		std::cout << "kappa = " << kappa << std::endl;
		if (kappa.upper() >= 1 ){
			std::cout << "Verification failed..." << std::endl;
		}

		std::cout << "S is invertible!" << std::endl;
		invS_norm = 1/(1 - kappa);
		std::cout << "\ninvS_norm = " << invS_norm << std::endl;
		vcp::time.toc();
	}

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
/************************************ Fixed Point Theorem **************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
vcp::time.tic();

{
		std::cout << "\nFixed Point Theorem Start" << std::endl;
		vcp::matrix< VData, VPOLICY > f2, g1, g2;
		f2 = czero_ball( ((CN * fduh_Linf_norm) * invS_norm * (CN * Res)) *  phi_i_norm );
		
		vcp::matrix< VData, VPOLICY > w_h_vec = 2 * lss(DF, f2);
		VData w_h_normH10 = sqrt(abs(transpose(w_h_vec)*DL*w_h_vec))(0);
		VData w_bot = 2 * invS_norm * CN * Res;

		std::cout << "Finite dimension || wh ||_H10 <= " << w_h_normH10 << std::endl;
		std::cout << "Inite dimension || w_bot ||_H10 <= " << w_bot << std::endl;

		vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
		Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, Order_legendre);
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();



		while (true) {
			w_h_normH10 = sqrt(abs(transpose(w_h_vec)*DL*w_h_vec))(0);
			Verification_Generator.setting_uh(w_h_vec, list_verification, 2);
			vcp::matrix< VData, VPOLICY > wh2phi = Verification_Generator.uhphi(2);
			vcp::matrix< VData, VPOLICY > wh2phiphi_diag = diag( Verification_Generator.uhphiphi(2));

//			std::cout << "wh2phi = " << std::endl;
//			std::cout << wh2phi << std::endl;

//			std::cout << "wh2phiphi_diag = " << std::endl;		
//			std::cout << wh2phiphi_diag << std::endl;

			g1 = wh2phi + czero_ball((2*CN*w_bot)*wh2phiphi_diag + pow(Cs4*w_bot,2)*phi_i_norm);
//			std::cout << "g1 = " << std::endl;		
//			std::cout << g1 << std::endl;

			g2 = czero_ball((pow(CN*Cs4,2)*fduh_Linf_norm*(fduh_Linf_norm*rhoh + 1)*invS_norm*(pow(w_h_normH10,2) + pow(w_bot,2)))*phi_i_norm);
//			std::cout << "g2 = " << std::endl;		
//			std::cout << g2 << std::endl;

			vcp::matrix< VData, VPOLICY > Kw_h = lss(DF, -f2 + g1 + g2);
			std::cout << "DFinv*f2=" << std::endl;
			std::cout << lss(DF,f2) << std::endl;
			std::cout << "End: DFinv*f2" << std::endl;
//			std::cout << "Kw_h = " << std::endl;		
//			std::cout << Kw_h << std::endl;

//			std::cout << "w_h_vec = " << std::endl;		
//			std::cout << w_h_vec << std::endl;

			bool FiniteFlag = true;
			for (int i = 0; i < w_h_vec.rowsize(); i++){
				if (!subset( Kw_h(i), w_h_vec(i) )){
					std::cout << "i = " << i << "," << Kw_h(i) << "," << w_h_vec(i) << std::endl;
					FiniteFlag = false;
				}
			}

			if (FiniteFlag) {
				std::cout << "Finite OK!" << std::endl;
				std::cout << " w_h in " << std::endl;				
				std::cout << w_h_vec << std::endl;
				std::cout << "End" << std::endl;				
			}

			bool InfiniteFlag = false;
			VData Kw_bot = invS_norm*( CN*Res + CN*pow(Cs4, 2)*( fduh_Linf_norm*rhoh + 1 )*(pow(w_h_normH10, 2) + pow(w_bot,2)) );
			std::cout << Kw_bot << std::endl;
			if (Kw_bot.upper() < w_bot.lower() ){
				std::cout << "Infinite OK!" << std::endl;
				InfiniteFlag = true;
			}

			if (FiniteFlag && InfiniteFlag){
				std::cout << "Error bound = " << std::endl;
				std::cout << sqrt(pow(w_h_normH10, 2) + pow(w_bot,2)) << std::endl;
				w_h_vec = Kw_h;
				w_bot = 1.01*Kw_bot;
			}
			else {
				break;
			}
		}
}
	
	return 0;
}
