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

#include <cmath>

//typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

typedef kv::interval<kv::mpfr<2000>> DataType;

//typedef vcp::pdblas POLICY;
typedef vcp::mats<AppData> POLICY;

typedef kv::interval<AppData> VAccData;
typedef vcp::imats<AppData> VAccPOLICY;

typedef double Data;
typedef kv::interval<Data> VData;
typedef vcp::pidblas VPOLICY;

typedef vcp::pdblas PDBLAS;

typedef kv::dd ResData;
typedef kv::interval<ResData> VResData;
typedef vcp::imats<ResData> VResPOLICY;

int main(void)
{
	std::cout.precision(17);
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/***********                                                          *************/" << std::endl;
	std::cout << "/*********** Verified computation for solution to Gray-Scott equation *************/" << std::endl;
	std::cout << "/***********                                                          *************/" << std::endl;
	std::cout << "/**********************************************************************************/" << std::endl;
	std::cout << "/**********************************************************************************/\n"
			  << std::endl;

	vcp::matrix<AppData, POLICY> uh, uh1, uh2, zh;
	vcp::Legendre_Bases_Generator<DataType, AppData, POLICY> Approximate_Generator;

	int epsilon = 100;
	std::cout << "Epsilon = " << epsilon << std::endl;

	int Order_legendre = 30;
	int uh_Order_legendre = 40;
	int p = 3;
	int Dimension = 2;
	int Number_of_variables = 2;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "Inverse Norm for Legendre Bases Order = " << Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

	// Setting of Approximate_Generator
	std::cout << "Setting the Generator by Approximate mode " << std::endl;
	Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 15);

	// Setting the list of Approximate_Generator
	//Approximate_Generator.setting_list();
	Approximate_Generator.setting_evenlist();

	// output the list => list_uh
	vcp::matrix<int> list_uh = Approximate_Generator.output_list();

	// setting initialization value of uh
	//	uh.rand(list_uh.rowsize(), Number_of_variables);
	//	uh = 30 * uh;
	//	uh(0, 0) = 40;
	//	uh(1, 0) = 50;
	//	uh(0, 1) = 40;
	//	uh(1, 1) = 50;

	vcp::load(uh, "Data_Gray_Scott/uh_solution1_eps55");
	//vcp::load(uh, "Data_Gray_Scott/uh_solution2_eps55");
	uh.resize(list_uh.rowsize(), Number_of_variables);

	std::cout << "Newton Method Start " << std::endl;
	{
		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix<AppData, POLICY> DL = Approximate_Generator.dphidphi();
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix<AppData, POLICY> L = Approximate_Generator.phiphi();

		vcp::matrix<AppData, POLICY> uh2vhphi;
		vcp::matrix<AppData, POLICY> uhphi;
		vcp::matrix<AppData, POLICY> vhphi;
		vcp::matrix<AppData, POLICY> fphi;
		vcp::matrix<AppData, POLICY> uhvhphiphi;
		vcp::matrix<AppData, POLICY> uh2phiphi;
		vcp::matrix<AppData, POLICY> DF;
		vcp::matrix<AppData, POLICY> F;
		vcp::matrix<AppData, POLICY> syuusei;
		vcp::matrix<AppData, POLICY> check;

		uh1 = uh.submatrix({0, uh.rowsize() - 1}, {0});
		uh2 = uh.submatrix({0, uh.rowsize() - 1}, {1});
		zh = vercat(uh1, uh2); // zh = [u; v]

		{
			AppData cc;
			while (1)
			{
				uh = horzcat(uh1, uh2); // uh = [u, v]

				Approximate_Generator.setting_uh(uh);
				uh2vhphi = Approximate_Generator.uhphi(2, 1);	  // vector (u^2 v, phi)
				uhphi = Approximate_Generator.uhphi(1, 0);		   // vector (u, phi)
				vhphi = Approximate_Generator.uhphi(0, 1);		   // vector (v, phi)
				fphi = Approximate_Generator.fphi();			   // vector (f, phi)
				uh2phiphi = Approximate_Generator.uhphiphi(2, 0);  // matrix (u^2 phi, phi)
				uhvhphiphi = Approximate_Generator.uhphiphi(1, 1); // matrix (u v phi, phi)

				DF = vercat(horzcat(DL - epsilon * (2 * uhvhphiphi - L), -epsilon * uh2phiphi), horzcat(2 * uhvhphiphi, DL - (-uh2phiphi - 10 * L)));
				F = vercat(DL * uh1 - epsilon * (uh2vhphi - uhphi), DL * uh2 - (-uh2vhphi + 10 * (fphi - vhphi)));
				syuusei = lss(DF, F);
				zh = zh - syuusei;
				check = max(abs(syuusei));
				cc = check(0);
				std::cout << cc << std::endl;
				if (cc < pow(2.0, -90))
				{
					std::cout << "Convergence \n"
							  << std::endl;
					uh1 = zh.submatrix({0, uh.rowsize() - 1}, {0});
					uh2 = zh.submatrix({uh.rowsize(), zh.rowsize() - 1}, {0});
					uh = horzcat(uh1, uh2);
					Approximate_Generator.setting_uh(uh);
					break;
				}
				uh1 = zh.submatrix({0, uh.rowsize() - 1}, {0});
				uh2 = zh.submatrix({uh.rowsize(), zh.rowsize() - 1}, {0});
			}
		}
	}
	// uh data for Grafics
	//	std::cout << "list_uh = " << std::endl;
	//	std::cout << list_uh << std::endl;
	//	vcp::save(list_uh, "Data_Gray_Scott/list_uh_solution3");

	//	std::cout << "uh = " << std::endl;
	//	std::cout << uh << std::endl;
	//	vcp::save(uh, "Data_Gray_Scott/uh_solution3");

	vcp::matrix<AppData, POLICY> Grafics = Approximate_Generator.output_uh_for_graphics(100);
	std::cout << "Grafics = " << std::endl;
	std::cout << Grafics << std::endl;
	//	vcp::save(Grafics, "Data_Gray_Scott/uh_Graf_solution3");

	// minimal and maximum value of approximate solution uh
	std::vector<kv::interval<double>> x;

	x.resize(Dimension);
	for (int d = 0; d < Dimension; d++)
	{
		x[d] = kv::interval<double>(0, 1);
	}
	std::vector<double> uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
	std::vector<double> uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

	for (int i = 0; i < Number_of_variables; i++)
	{
		std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/******************* Calculate Inverse Norm || F'[uh]^-1 ||_(H-1,H10) <= K *********************/
	/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData CN, Cs6, C2;
	VData K = VData(1);
	vcp::matrix<VData, VPOLICY> uh_max_min;
	uh_max_min.zeros(Number_of_variables, 1);
	for (int i = 0; i < Number_of_variables; i++)
	{
		uh_max_min(i) = kv::interval<double>(uh_min[i], uh_max[i]);
	}
	
	{
		std::cout << "\n>>> Calculate Inverse Norm || F'[uh]^-1 ||_(H-1,H10) <= K" << std::endl;
		vcp::Legendre_Bases_Generator<DataType, VAccData, VAccPOLICY> Verification_Generator;
		Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
		Verification_Generator.setting_list();
		//Verification_Generator.setting_evenlist();
		vcp::matrix<VAccData, VAccPOLICY> uhi;
		vcp::interval(uh, uhi);
		// uh setting : Last Argument is list divide : full list => 1 , even list => 2
		Verification_Generator.setting_uh(uhi, list_uh, 2);

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix<VAccData, VAccPOLICY> DLAcc = Verification_Generator.dphidphi();
		vcp::matrix<VAccData, VAccPOLICY> LAcc = Verification_Generator.phiphi();
		vcp::matrix<VAccData, VAccPOLICY> uh2phiphiAcc = Verification_Generator.uhphiphi(2, 0);  // matrix (u^2 phi, phi)
		vcp::matrix<VAccData, VAccPOLICY> uhvhphiphiAcc = Verification_Generator.uhphiphi(1, 1); // matrix (u v phi, phi)

		vcp::matrix<VAccData, VAccPOLICY> Zero;
		Zero.zeros(DLAcc.rowsize(), DLAcc.columnsize());

		DLAcc = vercat(horzcat(DLAcc, Zero), horzcat(Zero, DLAcc));

		vcp::matrix<VAccData, VAccPOLICY> NAcc = vercat(horzcat(epsilon * (2 * uhvhphiphiAcc - LAcc), epsilon * uh2phiphiAcc), horzcat(-2 * uhvhphiphiAcc, -uh2phiphiAcc - 10 * LAcc));
		vcp::matrix<VAccData, VAccPOLICY> NTAcc = vercat(horzcat(epsilon * (2 * uhvhphiphiAcc - LAcc), -2 * uhvhphiphiAcc), horzcat(epsilon * uh2phiphiAcc, -uh2phiphiAcc - 10 * LAcc));

		uh2phiphiAcc.clear();
		uhvhphiphiAcc.clear();

		LAcc = vercat(horzcat(LAcc, Zero), horzcat(Zero, LAcc));
		Zero.clear();

		vcp::matrix<VAccData, VAccPOLICY> DFAcc = DLAcc - (NAcc + NTAcc);
		vcp::matrix<VData, VPOLICY> DF, DL, L, N, NT;

		vcp::convert(DFAcc, DF);
		DFAcc.clear();
		vcp::convert(DLAcc, DL);
		DLAcc.clear();
		vcp::convert(LAcc, L);
		LAcc.clear();
		vcp::convert(NAcc, N);
		NAcc.clear();
		vcp::convert(NTAcc, NT);
		NTAcc.clear();

		DF = DF + N * lss(DL, NT);

		DL.clear();
		N.clear();
		NT.clear();

		vcp::matrix<VData, VPOLICY> E;
		compsym(DF); // Forced symmetrization
		compsym(L);
		eigsymge(DF, L, E);

		E = diag(E);
		std::cout << min(E) << std::endl;

		// Calculate some constants
		C2 = Verification_Generator.Poincare_constant<VData>();
		std::cout << "C2 = " << C2 << std::endl;

		Cs6 = Verification_Generator.Sobolev_constant<VData>(6);
		std::cout << "Cs6 = " << Cs6 << ", p ="
				  << "6" << std::endl;

		CN = Verification_Generator.Ritz_projection_error<VData>();
		std::cout << "CN = " << CN << std::endl;
		Verification_Generator.clear();

		vcp::matrix<VData, VPOLICY> Q_max_min;
		Q_max_min.zeros(Number_of_variables, Number_of_variables);
		Q_max_min(0, 0) = epsilon * (2 * uh_max_min(0) * uh_max_min(1) - 1);
		Q_max_min(0, 1) = epsilon * pow(uh_max_min(0), 2);
		Q_max_min(1, 0) = -2 * uh_max_min(0) * uh_max_min(1);
		Q_max_min(1, 1) = -pow(uh_max_min(0), 2) - 10;

		vcp::matrix<VData, VPOLICY> QQ_norm;
		Q_max_min = intervalmag(Q_max_min + transpose(Q_max_min));
		QQ_norm = normtwo(Q_max_min);

		std::cout << "QQ_norm = " << std::endl;
		std::cout << QQ_norm << std::endl;

		VData sigma = QQ_norm(0) + 0.01;
		std::cout << "sigma = " << std::endl;
		std::cout << sigma << std::endl;

		Q_max_min.zeros(Number_of_variables, Number_of_variables);
		Q_max_min(0, 0) = sigma - epsilon * (2 * uh_max_min(0) * uh_max_min(1) - 1);
		Q_max_min(0, 1) = -epsilon * pow(uh_max_min(0), 2);
		Q_max_min(1, 0) = 2 * uh_max_min(0) * uh_max_min(1);
		Q_max_min(1, 1) = sigma + pow(uh_max_min(0), 2) + 10;

		vcp::matrix<VData, VPOLICY> Qsigma_norm;
		Q_max_min = intervalmag(Q_max_min + transpose(Q_max_min));
		Qsigma_norm = normtwo(Q_max_min);

		std::cout << "Qsigma_norm = " << std::endl;
		std::cout << Qsigma_norm << std::endl;

		Q_max_min.zeros(Number_of_variables, Number_of_variables);
		Q_max_min(0, 0) = epsilon * (2 * uh_max_min(0) * uh_max_min(1) - 1);
		Q_max_min(0, 1) = epsilon * pow(uh_max_min(0), 2);
		Q_max_min(1, 0) = -2 * uh_max_min(0) * uh_max_min(1);
		Q_max_min(1, 1) = -pow(uh_max_min(0), 2) - 10;

		vcp::matrix<VData, VPOLICY> Q_norm;
		Q_max_min = intervalmag(Q_max_min);
		Q_norm = normtwo(Q_max_min);

		std::cout << "Q_norm = " << std::endl;
		std::cout << Q_norm << std::endl;

		VData mu_sh1 = E(0) + sigma;

		std::cout << "mu_sh1 = " << std::endl;
		std::cout << mu_sh1 << std::endl;

		VData C_Msigma = CN * sqrt(1 + pow(CN, 2) * (Qsigma_norm(0) + pow(C2 * Q_norm(0), 2)));

		std::cout << "C_Msigma = " << std::endl;
		std::cout << C_Msigma << std::endl;

		vcp::matrix< VData, VPOLICY > mu_sh;
		vcp::matrix< VData, VPOLICY > Lambda;
		Lambda = E;
        mu_sh = E + sigma;
		VData mu_low;
        for (int i = 0; i < mu_sh.rowsize(); i++){
            mu_low = mu_sh(i)/(1 + pow(C_Msigma,2)*mu_sh(i) );
            mu_sh(i).lower() = mu_low.lower();
			Lambda(i).lower() = (mu_sh(i) - sigma).lower();
        }

		VData lambda = Lambda(0);
		std::cout << "lambda = " << std::endl;
		std::cout << lambda << std::endl;

		if (lambda.lower() > 0){
            using std::sqrt;
            K = 1/sqrt(lambda);

            std::cout << "Pre K = " << std::endl;
            std::cout << K << std::endl;
            sigma = 0;
        }
        else {
            std::cout << "Can not estimate Pre K" << std::endl;
            sigma = VData(abs(lambda.lower())) + 1;
        }



///////////////////////////////////////////////////////////////
        std::cout << "Goerisch method Start!" << std::endl;
        std::cout << "M(u, v) = (Du, Dv) - ((Q + QT)u, v) + (Rh Ainv QT u, QT) + sigma(u, v) + r(u, v)" << std::endl;
        std::cout << "N(u, v) = (u, v)" << std::endl;
        std::cout << "M(u, v) = lambda N(u, v)" << std::endl;
        std::cout << "B = Q*Rh*Ainv*QT - ( Q+QT) + sigma*I " << std::endl;        

        std::cout << "sigma = " << sigma << std::endl;
        std::cout << ">>> Find rho_mush" << std::endl;
        mu_sh = Lambda + sigma;

        VData rho_mush;
        for (int i = mu_sh.rowsize() - 1; i > 0 ; i--){
            if (!overlap(mu_sh(i), mu_sh(i-1))){
                int j = i - 1; // from 0 to  N-1 
                rho_mush = VData( mu_sh(i).lower() );
                std::cout << "j+1 = " << j+1 << ", rho_mush = " << rho_mush << std::endl;
                break;
            }
        }
        std::cout << "mu_sh size =" << mu_sh.rowsize() << std::endl;
		
		std::cout << ">>> Create approximate eigenvalue and approximate eigenfunction" << std::endl;

        vcp::matrix< VData, VPOLICY > DFsigma, z;
        vcp::matrix< Data, PDBLAS > DFh, Lh, Eh, Uh;

        DFsigma = DF + sigma*L;
        compsym(DFsigma);

        mid(DFsigma, DFh);
        mid(L, Lh);

		eigsymge(DFh, Lh, Eh, Uh);
        DFh.clear();
        Lh.clear();
        Eh = diag(Eh);
        Eh = Eh.submatrix({ 0, N - 1 }, { 0 });
        interval(Eh, E);
        std::cout << "E = " << E << std::endl;

        std::cout << ">>> Create matrixes A0, A1, and A2" << std::endl;
        vcp::matrix< VData, VPOLICY > U, A0, A1, A2, Er, Er2, Ert, Mu, ones, A, B;


        VData r = VData(1);
        Uh = Uh.submatrix({}, { 0, N - 1 });
        interval(Uh, U);
        z = lss(D, transpose(Qmat)*U);

        A0 = transpose(U)*(DFsigma + r*L)*U;
        compsym(A0);

        A1 = transpose(U)*L*U;
        compsym(A1);

        vcp::matrix< VAccData, VAccPOLICY > u_eigenvector, z_acc;
        convert(Uh, u_eigenvector);
        convert(z, z_acc);
		

	}
	vcp::time.toc();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/******************* Calculate Residual Norm || Laplace(uh) - f(uh) ||_L2 **********************/
	/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData Res = VData(0);
	VData Res1 = VData(0);
	VData Res2 = VData(0);
	{
		std::cout << "\n>>> Calculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
		vcp::Legendre_Bases_Generator<DataType, VResData, VResPOLICY> Verification_Generator;
		Verification_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 2);
		// Setting the list of Verification_Generator
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();
		vcp::matrix<VResData, VResPOLICY> uhi;
		vcp::interval(uh, uhi);
		//vcp::convert(uh, uhi);

		Verification_Generator.setting_uh(uhi);

		// 1. First form
		// F1(uh, vh) = - DL * uh - epsilon*(uh2vhphi - uhphi)
		VResData LuhLuh = Verification_Generator.integral_LuhLuh(0);
		VResData Luh_uh2vh = Verification_Generator.integral_Luhuh(0, 2, 1);
		VResData Luh_uh = Verification_Generator.integral_Luhuh(0, 1, 0);
		VResData uh4vh2 = Verification_Generator.integral_uh(4, 2);
		VResData uh3vh = Verification_Generator.integral_uh(3, 1);
		VResData uh2 = Verification_Generator.integral_uh(2, 0);
		{
			using std::abs;
			using std::sqrt;
			vcp::convert(sqrt(abs(LuhLuh + 2 * epsilon * (Luh_uh2vh - Luh_uh) + epsilon * epsilon * (uh4vh2 - 2 * uh3vh + uh2))), Res1);
			std::cout << "First Residual Norm : || Laplace(uh) - f1(uh, vh) ||_L2 <= " << Res1 << std::endl;
		}

		// 2. Second form
		// F2(uh, vh) = - DL * vh - ( -uh2v + 10( 1 - vh ) )
		VResData LvhLvh = Verification_Generator.integral_LuhLuh(1);
		VResData Lvh_uh2vh = Verification_Generator.integral_Luhuh(1, 2, 1);
		VResData Lvh_f = Verification_Generator.integral_Luhuh(1, 0, 0);
		VResData Lvh_vh = Verification_Generator.integral_Luhuh(1, 0, 1);
		VResData uh2vh = Verification_Generator.integral_uh(2, 1);
		VResData uh2vh2 = Verification_Generator.integral_uh(2, 2);
		VResData vh = Verification_Generator.integral_uh(0, 1);
		VResData vh2 = Verification_Generator.integral_uh(0, 2);

		{
			using std::abs;
			using std::sqrt;
			vcp::convert(sqrt(abs(LvhLvh + 2 * (-Lvh_uh2vh + 10 * (Lvh_f - Lvh_vh)) + uh4vh2 + 20 * (uh2vh2 - uh2vh) + 100 * (1 + vh2) - 200 * vh)), Res2);
			std::cout << "Second Residual Norm : || Laplace(vh) - f2(uh, vh) ||_L2 <= " << Res2 << std::endl;
		}
		using std::sqrt;
		Res = sqrt(Res1 * Res1 + Res2 * Res2);
		std::cout << "Residual Norm : || F(uh, vh) ||_X <= " << Res << std::endl;
		//		Res = Cp * Res;
		//		std::cout << "Residual Norm : || F(uh) ||_(H-1) <= " << Res << std::endl;
	}
	vcp::time.toc();

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/**** Calculate Lipschitz constant || F'[w1] - F'[w2] ||_(H-1,H10) <= G || w1 - w2 ||_(H10) ****/
	/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData G = VData(0);
	{
		std::cout << "\n>>> Calculate Lipschitz constant || F'[w1] - F'[w2] ||_(H-1,H10) <= G || w1 - w2 ||_(H10)" << std::endl;
		using std::pow;

		vcp::Legendre_Bases_Generator<DataType, VData, VPOLICY> Verification_Generator;
		Verification_Generator.setting(uh_Order_legendre, 1, Dimension, Number_of_variables, 1, uh_Order_legendre);
		// Setting the list of Verification_Generator
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix<VData, VPOLICY> DL = Verification_Generator.dphidphi();
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix<VData, VPOLICY> L = Verification_Generator.phiphi();

		vcp::matrix<VData, VPOLICY> uhi;
		vcp::convert(uh, uhi);

		vcp::matrix<VData, VPOLICY> uh_H10_norm = sqrt(abs(transpose(uhi) * (DL * uhi)));

		std::cout << "|| uh ||_(H10) = " << uh_H10_norm << std::endl;

		vcp::matrix<VData, VPOLICY> Gm;
		vcp::matrix<VData, VPOLICY> E;
		Gm.zeros(Number_of_variables, Number_of_variables);

		std::cout << "Cs6 = " << Cs6 << ", p ="
				  << "6" << std::endl;

		Gm(0, 0) = 2 * epsilon * pow(Cs6, 3) * (2 * uh_H10_norm(0, 0) + uh_H10_norm(1, 1) + 6 * K * Res);
		Gm(0, 1) = epsilon * pow(Cs6, 3) * (uh_H10_norm(0, 0) + 2 * K * Res);
		Gm(1, 0) = 2 * pow(Cs6, 3) * (2 * uh_H10_norm(0, 0) + uh_H10_norm(1, 1) + 6 * K * Res);
		Gm(1, 1) = pow(Cs6, 3) * (uh_H10_norm(0, 0) + 2 * K * Res);

		std::cout << "Gm = " << std::endl;
		std::cout << Gm << std::endl;

		Gm = ltransmul(Gm);

		std::cout << "Gm^T * Gm = " << std::endl;
		std::cout << Gm << std::endl;

		eigsym(Gm, E);
		E = diag(E);

		G = sqrt(abs(E(1, 0)));

		std::cout << "G = " << std::endl;
		std::cout << G << std::endl;
	}
	vcp::time.toc();

	return 0;
}
