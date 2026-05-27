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

#include <vcp/allsol/candidate_set.hpp>
#include <vcp/allsol/solution_list.hpp>
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

template <typename _T, typename _TM, class _PM = vcp::mats< _TM >>
struct Problem : public vcp::candidate_set< _T, _TM, _PM > {
	vcp::matrix< _TM, _PM > uh;

	void setting_uh( const vcp::matrix< _TM, _PM >& uhh ){
		
		this->Generator.setting_uh(uhh, this->constant.list_uh, 1);
		this->uh = uhh;
	}

	_TM g(const _TM& alpha )override{  
		return sqrt(pow(this->constant.Cs[2]*this->constant.alpha,2) + pow(this->constant.Cs[6]*this->constant.alpha,6) )*LAMBDA;
	}

	vcp::matrix< _TM, _PM > f() {
        vcp::matrix< _TM, _PM > uh3phi = this->Generator.uhphi(3);
        return this->D*this->uh - LAMBDA*(this->L*this->uh - uh3phi);
    }

    vcp::matrix< _TM, _PM > Df() {
        vcp::matrix< _TM, _PM > uh2phiphi = this->Generator.uhphiphi(2);
        return this->D - LAMBDA*( this->L - _TM(3)*uh2phiphi);
    }
};


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

	int N = 20;
	int p = 3;
	int Dimension = 1;
	int Number_of_variables = 1;

	VData alpha = VData(2);
	vcp::solution_list< VData, VPOLICY > Candidate_Sets_List;
	Problem< DataType, VData, VPOLICY > Target_Set;


	std::cout << "lambda = " << LAMBDA << std::endl;
	std::cout << "alpha = " << alpha << std::endl;
	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order N = " << N << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

	Target_Set.setting< VData >(N, p, Dimension, Number_of_variables);
	Target_Set.setting_alpha(alpha);
	Target_Set.calc_rho_o();
	Candidate_Sets_List.push_back_unknown( Target_Set.calc_UN() );

	
	std::cout << "/*********** Start Finding All solutions *************/" << std::endl;
	vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
	Verification_Generator.setting(N, p, Dimension, Number_of_variables, 1, N);
	Verification_Generator.setting_list();

	while ( Candidate_Sets_List.unknown_size() > 0 ){
		std::cout << "============================= Status ===================================" << std::endl;
		std::cout << "remaining amount:" << Candidate_Sets_List.unknown_size() << ", solutions:" << Candidate_Sets_List.solutions_size() << ", Non existence:" << Candidate_Sets_List.nonexistence_size() << ", Out of range:" << Candidate_Sets_List.out_of_range_count_size() << std::endl;
		std::cout << "========================================================================\n" << std::endl;
		std::cout << "========================================================================" << std::endl;
		std::cout << "========================= New Target Start =============================" << std::endl;
		std::cout << "========================================================================" << std::endl;

		vcp::matrix< VData, VPOLICY > Uh = Candidate_Sets_List.pop_back_unknown();
		std::cout << "Target: Uh = " << std::endl;
		std::cout << Uh << std::endl;
		std::cout << "Target: Uo = " << std::endl;
		std::cout << "{ u in H10 | || u || <= " <<Target_Set.constant.rho_o << " }" << std::endl;

		
		std::cout << "===== Check Target Range: alpha < || Uh ||^2 + rho_0 || u_o ||^2? ======" << std::endl;
		Target_Set.constant.rho_h = (transpose(Uh)*Target_Set.D*Uh)(0);

		if (Target_Set.constant.rho_h.lower() < 0 ){
			Target_Set.constant.rho_h.lower() = 0;
		}

		Target_Set.constant.rho_h = sqrt(Target_Set.constant.rho_h);
		Target_Set.constant.rho = sqrt(pow(Target_Set.constant.rho_h, 2) + pow(Target_Set.constant.rho_o, 2));

		std::cout << "target: rho_h = " <<Target_Set.constant.rho_h << ", rho_o = " <<Target_Set.constant.rho_o << ", rho = " <<Target_Set.constant.rho << std::endl;
		if (Target_Set.constant.rho.lower() >Target_Set.constant.alpha.upper() ){
			Candidate_Sets_List.out_of_range_counter();
			std::cout << "Out of range the Target Set!" << std::endl;
			continue;
		}

		vcp::matrix< VData, VPOLICY > uhat;
		imid( Uh, uhat );
		Target_Set.setting_uh(uhat);
		auto Fuh = Target_Set.f();


		std::cout << "=========================== Nonexistence Part I ========================" << std::endl;
		bool nonexistence_flag = false;
		{
			Target_Set.constant.Cfdw = LAMBDA*(Target_Set.constant.Cs[2] + 3*pow(Target_Set.constant.Cs[6],3)*pow(Target_Set.constant.rho,2) );
			std::cout << "Cfdw = " <<Target_Set.constant.Cfdw << std::endl;

			// || uN-uh ||_H10^2
			auto tmp = (transpose(Uh - uhat)*Target_Set.D*(Uh - uhat))(0);
			Target_Set.constant.W_N_H10_2 = tmp;
			tmp.lower() = tmp.upper();

			// sqrt(|| uN-uh ||_H10^2 + rho_o^2)
			tmp = sqrt(tmp + pow(Target_Set.constant.rho_o,2));
			auto tmpm =Target_Set.constant.Cfdw*tmp*Target_Set.norm_phii_L2;
			for (int i = 0; i < tmpm.rowsize(); i++ ) tmpm(i).lower() = -tmpm(i).upper();

			auto QN = Fuh + Target_Set.D*(Uh - uhat) + tmpm;
			for (int i = 0; i < QN.rowsize(); i++){
				if ( !(QN(i).upper() >= 0 && QN(i).lower() <= 0) ){
					nonexistence_flag = true;
					std::cout << "Non existence a solution in the Target Set!" << std::endl;
					Candidate_Sets_List.nonexistence_counter();					
					break;
				}
			}
		}
		if ( nonexistence_flag ) continue;

		using std::pow;
		auto T =Target_Set.Df();
		auto TinvFuh = lss(T, Fuh);

		std::cout << "=========================== Nonexistence Part II =======================" << std::endl;
		{
			std::vector< kv::interval< double > > x;
			x.resize(Dimension);
			for (int d = 0; d < Dimension; d++) {
				x[d] = kv::interval< double >(0, 1);
			}

			auto uh_min = Target_Set.Generator.global_min(x, std::pow(2.0, -3));
			auto uh_max = Target_Set.Generator.global_max(x, std::pow(2.0, -3));
			
			for (int i = 0; i < Number_of_variables; i++) {
				std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
			}
			
			Target_Set.constant.uh_minmax.lower() = uh_min[0];
			Target_Set.constant.uh_minmax.upper() = uh_max[0];
		}
		{
			Target_Set.constant.CfdX = abs( VData( LAMBDA )*(1 - 3*pow(Target_Set.constant.uh_minmax,2)) );
			Target_Set.constant.Cfd_N =Target_Set.constant.Cs[2]*Target_Set.constant.CfdX;
			Target_Set.constant.Cfd_O =Target_Set.constant.CN*Target_Set.constant.CfdX;		
			Target_Set.constant.C2 =Target_Set.constant.CN*Target_Set.constant.Cfd_N;
			Target_Set.constant.C3 =Target_Set.constant.CN*Target_Set.constant.Cfd_O;


			auto tmp = sqrt(abs((transpose(uhat)*Target_Set.D*uhat)(0)));
			auto W_rho = sqrt(abs(Target_Set.constant.W_N_H10_2 + pow(Target_Set.constant.rho_o,2)));
			Target_Set.constant.G =abs(3*VData(LAMBDA)*pow(Target_Set.constant.Cs[3],3)*(2*tmp + W_rho)*W_rho);

			auto f1 =  Target_Set.constant.rho*Target_Set.constant.Cfd_O*Target_Set.norm_phii_L2;
			auto f2 = Target_Set.constant.G*Target_Set.constant.rho*Target_Set.norm_phii_L2;

			auto KU = lss(T, Uh + Fuh + f1 + f2);

			std::cout << "Check condition Nonexistence II" << std::endl;
			for ( int i = 0; i < KU.rowsize(); i++ ) {

				if ( !overlap(Uh(i), KU(i)) ) {
					nonexistence_flag = true;
					std::cout << Uh(i) << ", " <<  KU(i) << std::endl;					
					std::cout << "Non existence a solution in the Target Set!" << std::endl;
					Candidate_Sets_List.nonexistence_counter();					
					break;
				}
			}
		}
		if ( nonexistence_flag ) continue;

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
			constant.delta =Target_Set.constant.CN*sqrt(abs(Res1 + Res2));
		//	std::cout << "delta = " <<Target_Set.constant.delta << std::endl;

			constant.uh_L6norm = pow(uh6, 1/VData(6));
		//	std::cout << "uh_L6norm = " <<Target_Set.constant.uh_L6norm << std::endl;
		}
		
*/
		Candidate_Sets_List.dvided_push_back_unknown( Uh );
	}

	return 0;
}
