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
#include <list>
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
#include "solution_list.hpp"
#include "constants.hpp"

#include <vcp/vcp_timer.hpp>

typedef double AppData;
//typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 4000 > > DataType;


typedef vcp::pdblas POLICY;
//typedef vcp::mats< AppData > POLICY;

typedef vcp::pidblas VPOLICY;

#define LAMBDA 100 

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

	int Order_legendre = 10;
	int uh_Order_legendre = Order_legendre;
	int p = 3;
	int Dimension = 1;
	int Number_of_variables = 1;

	vcp::constants< VData > constant;
	constant.alpha = 2;
	vcp::solution_list< VData, VPOLICY > List;
	vcp::matrix< VData, VPOLICY > D, L;
	vcp::matrix< VData, VPOLICY > norm_phii_L2;
	vcp::matrix< int > list_uh;


	std::cout << "lambda = " << LAMBDA << std::endl;
	std::cout << "alpha = " << constant.alpha << std::endl;
	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "Inverse Norm for Legendre Bases Order = " << Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;



	std::cout << "/*********** Initialization *************/" << std::endl;
	vcp::time.tic();
	{
		vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
		Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
		//Verification_Generator.setting_list();
		Verification_Generator.setting_list();
		list_uh = Verification_Generator.output_list();

		// Make the matrix (\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		D = Verification_Generator.dphidphi();

		// Make the matrix (\phi_i, \phi_j)_{L^2})_{i,j}
		L = Verification_Generator.phiphi();
		norm_phii_L2 = sqrt(diag(L));
		
		constant.CN = Verification_Generator.Ritz_projection_error< VData >();
		constant.CN.lower() = constant.CN.upper();

		constant.Cs2 = Verification_Generator.Poincare_constant< VData >();
		constant.Cs2.lower() = constant.Cs2.upper();

		constant.Cs6 = Verification_Generator.Sobolev_constant< VData >(3);
		constant.Cs6.lower() = constant.Cs6.upper();

		constant.g_alpha = sqrt(pow(constant.Cs2*constant.alpha,2) + pow(constant.Cs6*constant.alpha,6) )*LAMBDA;
		constant.g_alpha.lower() = constant.g_alpha.upper();


		std::cout << "g(alpha) =" << constant.g_alpha << std::endl;
		std::cout << "CN = " << constant.CN << std::endl;

		vcp::matrix< VData, VPOLICY > phi_vec = sqrt(diag(L));
	//	vcp::matrix< VData, VPOLICY > phi_i_L2_vec = phi_vec;
		phi_vec = abs(constant.g_alpha*phi_vec);


		for (int i = 0; i < phi_vec.rowsize(); i++){
			phi_vec(i).lower() = -phi_vec(i).upper();
		}

		vcp::matrix< VData, VPOLICY > Ustar = lss(D, phi_vec);
		List.push_back_unknown( Ustar );

		std::cout << "Ustar =" << std::endl;
		std::cout << Ustar << std::endl;

		constant.rho_o = constant.CN*constant.g_alpha;
		constant.rho_o.lower() = 0;
		std::cout << "rho_o =" << constant.rho_o << std::endl;

/*
		vcp::matrix< AppData, POLICY > tmpc;
		mid(D, tmpc);
		tmpc = Cholesky(tmpc);
		interval( tmpc, DmidChol );
		constant.DmidChol_Res = norminf(transpose( DmidChol )*DmidChol - D)(0);		
*/		
	}
	vcp::time.toc();


	std::cout << "/*********** Start Finding All solutions *************/" << std::endl;
	vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
	Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
	Verification_Generator.setting_list();	
	while ( List.unknown_size() > 0 ){
		std::cout << "remaining amount:" << List.unknown_size() << ", solutions:" << List.solutions_size() << ", Non existence:" << List.nonexistence_size() << ", Out of range:" << List.out_of_range_count_size() << std::endl;
		vcp::matrix< VData, VPOLICY > Uh = List.pop_back_unknown();
	//	std::cout << "Uh = " << std::endl;
	//	std::cout << Uh << std::endl;

		constant.rho_h = (transpose(Uh)*D*Uh)(0);

		if ( constant.rho_h.lower() < 0 ){
			constant.rho_h.lower() = 0;
		}

		constant.rho_h = sqrt(constant.rho_h);
		constant.rho = sqrt(pow(constant.rho_h, 2) + pow( constant.rho_o, 2));

		std::cout << "target: rho_h = " << constant.rho_h << ", rho_o = " << constant.rho_o << ", rho = " << constant.rho << std::endl;
		if ( constant.rho.lower() > constant.alpha.upper() ){
			List.out_of_range_counter();
			std::cout << Uh << std::endl;
			continue;
		}

		vcp::matrix< VData, VPOLICY > uhat;
		imid( Uh, uhat );
		Verification_Generator.setting_uh(uhat, list_uh, 1);
		auto uhphi = Verification_Generator.uhphi(1);
		auto uh3phi = Verification_Generator.uhphi(3);
		auto Fuh = D*uhat - LAMBDA*( uhphi - uh3phi );


//		std::cout << " ------------- Check: Nonexistence Part I ----------------" << std::endl;
		bool nonexistence_flag = false;
		{
			constant.Cfdw = LAMBDA*( constant.Cs2 + 3*pow(constant.Cs6,3)*pow(constant.rho,2) );
	//		std::cout << "Cfdw = " << constant.Cfdw << std::endl;

			// || uN-uh ||_H10^2
			auto tmp = (transpose(Uh - uhat)*D*(Uh - uhat))(0);
			tmp.lower() = tmp.upper();

			// sqrt(|| uN-uh ||_H10^2 + rho_o^2)
			tmp = sqrt(tmp + pow(constant.rho_o,2));
			auto tmpm = constant.Cfdw*tmp*norm_phii_L2;
			for (int i = 0; i < tmpm.rowsize(); i++ ) tmpm(i).lower() = -tmpm(i).upper();

			auto QN = Fuh + D*(Uh - uhat) + tmpm;
			for (int i = 0; i < QN.rowsize(); i++){
				if ( !(QN(i).upper() >= 0 && QN(i).lower() <= 0) ){
					nonexistence_flag = true;
				//	std::cout << "===============================================================================" << std::endl;
				//	std::cout << "=====================================Yatta=====================================" << std::endl;
				//	std::cout << "===============================================================================" << std::endl;
				//	std::cout << Uh << std::endl;
					List.nonexistence_counter();					
					break;
				}
			}
			std::cout << QN << std::endl;
		}
		if ( nonexistence_flag ) continue;


		auto uh2phiphi = Verification_Generator.uhphiphi(2);

		using std::pow;
		auto T = D - LAMBDA*( L - 3.0*uh2phiphi);
		auto TinvFuh = lss(T, Fuh);
		auto vhat = uhat - TinvFuh;

		/*
		{
			std::vector< kv::interval< double > > x;
			x.resize(Dimension);
			for (int d = 0; d < Dimension; d++) {
				x[d] = kv::interval< double >(0, 1);
			}

			auto uh_min = Verification_Generator.global_min(x, std::pow(2.0, -3));
			auto uh_max = Verification_Generator.global_max(x, std::pow(2.0, -3));
			
			for (int i = 0; i < Number_of_variables; i++) {
				std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
			}
			
			constant.uh_minmax.lower() = uh_min[0];
			constant.uh_minmax.upper() = uh_max[0];
		}
		*/
		constant.uh_minmax = constant.rho_h/2;
		constant.uh_minmax.lower() = -constant.uh_minmax.upper();

		constant.CfdX = abs( LAMBDA*(1 - 3*pow(constant.uh_minmax,2)) );
		constant.Cfd_N = constant.Cs2*constant.CfdX;
		constant.Cfd_O = constant.CN*constant.CfdX;		
		constant.C2 = constant.CN*constant.Cfd_N;
		constant.C3 = constant.CN*constant.Cfd_O;

		std::cout << "C2 = " << constant.C2 << std::endl;
		std::cout << "C3 = " << constant.C3 << std::endl;
		
		if ( constant.C3.upper() >= 1 ) {
			List.dvided_push_back_unknown( Uh );
			continue;	
		}
/*
		///////////////////// Compute delta ////////////////////////
		{
			VData Res1, Res2;
			VData uh6;
			{
			//	std::cout << "\nCalculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
				vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Residual_Generator;

				Residual_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 2);
				Residual_Generator.setting_list();	
				Residual_Generator.setting_uh(uhat);

				// -DL * uh - (a*uh + b*uh2 + c*uh3)        
				VData LuhLuh = Residual_Generator.integral_LuhLuh(0); // (DL uh, DL uh)
				VData Luh_uh = Residual_Generator.integral_Luhuh(0, 1); // (DL uh, uh)
				VData Luh_uh3 = Residual_Generator.integral_Luhuh(0, 3); // (DL uh, uh3)                
				VData uh2 = Residual_Generator.integral_uh(2); // (uh, uh)        
				VData uh4 = Residual_Generator.integral_uh(4); // (uh, uh3)
				uh6 = Residual_Generator.integral_uh(6); // (uh3, uh3)
				
				{
					using std::sqrt;
					using std::abs;
					Res1 = LuhLuh + VData(2)*( LAMBDA*Luh_uh + LAMBDA*Luh_uh3 ) + ( LAMBDA*LAMBDA*uh2 + VData(2)*LAMBDA*LAMBDA*uh4 + LAMBDA*LAMBDA*uh6 );
				}
			}
			{
				vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Residual_Generator;

				Residual_Generator.setting(Order_legendre, p, Dimension, Number_of_variables + 1, 2);
				Residual_Generator.setting_list();	
				Residual_Generator.setting_uh( horzcat(lss(D, Fuh), uhat) );

				// (Delta psi, -Delta uh - LAMBDA*uh - LAMBDA*uh3  + Delta psi )
				//    =   (Delta psi, Delta psi )
				//      - 2*(Delta psi, Delta uh)
				//      - 2*LAMBDA*(Delta psi, uh )
				//      - 2*LAMBDA*(Delta psi, uh3 )
				VData LpsiLpsi = Residual_Generator.integral_LuhLuh(0);
				VData LuhLpsi = Residual_Generator.integral_LuhLuh(0, 1);		
				VData Lphi_uh = Residual_Generator.integral_Luhuh(0, 0, 1);
				VData Lphi_uh3 = Residual_Generator.integral_Luhuh(0, 0, 3);

				Res2 = LpsiLpsi - VData(2)*( LuhLpsi + LAMBDA*Lphi_uh + LAMBDA*Lphi_uh3 );
			}
			constant.delta = constant.CN*sqrt(abs(Res1 + Res2));
		//	std::cout << "delta = " << constant.delta << std::endl;

			constant.uh_L6norm = pow(uh6, 1/VData(6));
		//	std::cout << "uh_L6norm = " << constant.uh_L6norm << std::endl;
		}
		
*/
		List.dvided_push_back_unknown( Uh );
	}

	return 0;
}
