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
#include <vcp/newton.hpp>

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


int a = 110;
int b = 7;
int c = 100;
int d = -8;


template < typename _T, typename _P, typename _TT = DataType >
struct Lotka_Volterra2 : public vcp::Newton< _T, _P > {
    int VOrder;
	int uh_Order;
	int p;
	int Dimension;
	int Number_of_variables;
    int list_size;
	int list_type;

    vcp::Legendre_Bases_Generator< _TT, _T, _P > Generator;
    vcp::matrix< _T, _P > DL;
    vcp::matrix< _T, _P > L;
    vcp::matrix< _T, _P > uh1, uh2, uh, z;

	void first_execute(){
        (*this).newton_tol = _T(256) * std::numeric_limits< _T >::epsilon();
        uh_Order = 20;
        VOrder = uh_Order;	    
	    p = 2;
	    Dimension = 2;
	    Number_of_variables = 2;
    }

    void setting_newton( vcp::matrix< _T, _P >& zh ) override {
        uh1 = zh.submatrix({ 0, list_size - 1 }, { 0 });
		uh2 = zh.submatrix({ list_size, zh.rowsize() - 1 }, { 0 });
        uh = horzcat(uh1, uh2);
        Generator.setting_uh(uh);
    }
    vcp::matrix< _T, _P > f() override {
        vcp::matrix< _T, _P > uhvhphi = Generator.uhphi(1, 1); // vector (u v, phi)
		vcp::matrix< _T, _P > uhphi = Generator.uhphi(1, 0); // vector (u, phi)
		vcp::matrix< _T, _P > uh2phi = Generator.uhphi(2, 0); // vector (u^2, phi)
		vcp::matrix< _T, _P > vhphi = Generator.uhphi(0, 1); // vector (v, phi)
		vcp::matrix< _T, _P > vh2phi = Generator.uhphi(0, 2); // vector (v^2, phi)

        return vercat(DL * uh1 - a*uhphi + uh2phi + c*uhvhphi, DL * uh2 - b*vhphi + d*uhvhphi + vh2phi);
    }
    vcp::matrix< _T, _P > Df() override {
        return (*this).make_D() - (*this).make_N();
    }

    // (D u, D v)
    vcp::matrix< _T, _P > make_D(){
        vcp::matrix< _T, _P > zero;
        zero.zeros(DL.rowsize(), DL.columnsize());
        return vercat(horzcat(DL, zero), horzcat(zero, DL));
    }

    // ( u,  v)
    vcp::matrix< _T, _P > make_L(){
        vcp::matrix< _T, _P > zero;
        zero.zeros(L.rowsize(), L.columnsize());
        return vercat(horzcat(L, zero), horzcat(zero, L));
    }

    // ( f'[uh] u,  v)
    vcp::matrix< _T, _P > make_N(){
        vcp::matrix< _T, _P > uhphiphi = Generator.uhphiphi(1, 0);
		vcp::matrix< _T, _P > vhphiphi = Generator.uhphiphi(0, 1);        
        return vercat(horzcat( a * L - 2 * uhphiphi - c * vhphiphi, -c * uhphiphi), horzcat(-d * vhphiphi, b * L - d * uhphiphi - 2 * vhphiphi));
    }

    // ( f'[uh]^* u,  v)
    vcp::matrix< _T, _P > make_NT(){
        vcp::matrix< _T, _P > uhphiphi = Generator.uhphiphi(1, 0);
		vcp::matrix< _T, _P > vhphiphi = Generator.uhphiphi(0, 1);
        return vercat(horzcat( a * L - 2 * uhphiphi - c * vhphiphi, -d * vhphiphi), horzcat(-c * uhphiphi, b * L - d * uhphiphi - 2 * vhphiphi));
    }

	void setting_uh( vcp::matrix< _T, _P > uhh){
		uh = uhh;
		uh1 = uhh.submatrix({ }, { 0 });
		uh2 = uhh.submatrix({ }, { 1 });
	}

    void set_VOrder(const int n){
        VOrder = n;
    }

    void Create_Generator(int mode, int ltype){
		if (mode == 1){
        	std::cout << "\n>>> Start Verification mode " << std::endl;
        	Generator.setting(VOrder, p, Dimension, Number_of_variables, 1, uh_Order);;
		}
		if (mode == 2){
        	std::cout << "\n>>> Start Residual mode " << std::endl;
        	Generator.setting(uh_Order, p, Dimension, Number_of_variables, mode);
		}
		if (mode >= 3){
        	std::cout << "\n>>> Start Approximate mode " << std::endl;
        	Generator.setting(uh_Order, p, Dimension, Number_of_variables, mode);
		}
		list_type = ltype;
		if (list_type == 1){
	        Generator.setting_list();
		}
		else if (list_type == 2){
			Generator.setting_evenlist();
		}

		if (mode != 2 ){
			DL = Generator.dphidphi();
			L = Generator.phiphi();

			list_size = DL.rowsize();
		}
    }

	vcp::matrix< _T, _P > decompose(vcp::matrix< _T, _P > zh){
		vcp::matrix< _T, _P > zh1 = zh.submatrix({ 0, list_size - 1 }, { 0 });
		vcp::matrix< _T, _P > zh2 = zh.submatrix({ list_size, zh.rowsize() - 1 }, { 0 });
        return horzcat(zh1, zh2);
	}

	void display(){
        std::cout.precision(17);
        std::cout << "/**********************************************************************************/" << std::endl;
        std::cout << "/**********************************************************************************/" << std::endl;
        std::cout << "/****                                                                          ****/" << std::endl;
        std::cout << "/**** Verified computation for solution to Lotla-Volterra competition equation ****/" << std::endl;
        std::cout << "/****                                                                          ****/" << std::endl;
        std::cout << "/**********************************************************************************/" << std::endl;
        std::cout << "/**********************************************************************************/\n" << std::endl;
        
        std::cout << "Dimension = " << Dimension << std::endl;
        std::cout << "uh of Legendre Bases Order = " << uh_Order << std::endl;
        std::cout << "Inverse Norm for Legendre Bases Order = " << VOrder << std::endl;
        std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
        std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;
    }

	void clear(){
        (*this).Generator.clear();
    }
};

template < typename _T, typename _P> vcp::matrix< _T, _P > make_Q(const vcp::matrix< _T, _P > &u){
    vcp::matrix< _T, _P > Q;
    Q.zeros(u.rowsize());
    Q(0, 0) = a - 2 * u(0) - c * u(1);
    Q(0, 1) = -c * u(0);
    Q(1, 0) = -d * u(1);
    Q(1, 1) = b - d * u(0) - 2 * u(1);
    return Q;
}

int main(void){

	std::cout.precision(17);
	vcp::matrix< int > list_uh;
	vcp::matrix< AppData, POLICY > uh;
	vcp::matrix< Kraw_Data, Kraw_POLICY > uhi;
	std::vector< double > uh_min;
	std::vector< double > uh_max;
	vcp::matrix< VData, VPOLICY > uh_max_min;

	int mode;
	int list_type = 2;


/////////////////////////////////////////////////////////////////////////////////////////////////
/***************************************** Create uh *******************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	{
		Lotka_Volterra2< AppData, POLICY > LV;
		LV.first_execute();
        LV.display();
		mode = 15;
        LV.Create_Generator(mode, list_type);

        vcp::matrix< AppData, POLICY > zh;
        list_uh = LV.Generator.output_list();
        uh.ones(list_uh.rowsize(), LV.Number_of_variables);
        uh(0, 0) = 50;
        uh(1, 0) = 20;
        uh(0, 1) = 50;
        uh(1, 1) = -5;

        zh = vercat(uh.submatrix({ 0, uh.rowsize() - 1 }, { 0 }), uh.submatrix({ 0, uh.rowsize() - 1 }, { 1 }));

        zh = LV.solve_nls(zh); // Compute an approximate solution
        LV.setting_newton(zh); // For Graphics and minimal/maximum value
        uh = LV.uh; // Approximate solution uh
	}
	{
		std::cout << "Krawczyk Method Start" << std::endl;
		Lotka_Volterra2< Kraw_Data, Kraw_POLICY > LV;
		LV.first_execute();
		mode = 1;
		LV.Create_Generator(mode, list_type);

		vcp::interval(uh, uhi);
		LV.Generator.setting_uh(uhi, list_uh, 2);
		LV.setting_uh(uhi);


		vcp::matrix< Kraw_Data, Kraw_POLICY > F = LV.f();
		vcp::matrix< Kraw_Data, Kraw_POLICY > R;
		vcp::matrix< Kraw_Data, Kraw_POLICY > DFxh = LV.Df();
		{			
			vcp::matrix< AppData, POLICY > DFmid;
//			vcp::matrix< AppData, POLICY > I;
//			I.eye(DFmid.rowsize());			
			vcp::convert(DFxh, DFmid);
			DFmid = inv(DFmid);

			R.zeros(DFmid.rowsize(), DFmid.columnsize());
			for (int i = 0; i < DFmid.rowsize(); i++){
				for (int j = 0; j < DFmid.columnsize(); j++){
					R(i, j) = Kraw_Data(DFmid(i, j));
				}
			}
		}

		vcp::matrix< Kraw_Data, Kraw_POLICY > Z = -R*F;
		vcp::matrix< Kraw_Data, Kraw_POLICY > kraw_init = 2*czero_ball(Z);
		vcp::matrix< Kraw_Data, Kraw_POLICY > I;
		std::cout << "check = " << kraw_init << std::endl;

		I.eye(Z.rowsize());

		int ite = 0;
		while(true){
			LV.Generator.setting_uh(uhi + LV.decompose(kraw_init), list_uh, 2);
			vcp::matrix< Kraw_Data, Kraw_POLICY > DFI = LV.Df();

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
				uhi = uhi + LV.decompose(kraw_init);
				std::cout << "uhi = " << std::endl;
				std::cout << uhi << std::endl;
				break;
			}
			kraw_init = Z + Kraw_Data(-0.1, 0.1)*Z;
			std::cout << "kraw_init(0) = " << kraw_init(0) << std::endl;
			std::cout << ite << std::endl;
			ite++;
		}
		vcp::time.toc();

		// minimal and maximum value of approximate solution uh
		vcp::time.tic();
		std::cout << "\nCalculate the maximum and minimum value" << std::endl;
		LV.Generator.setting_uh(uhi, list_uh, list_type);

		std::vector< kv::interval< double > > x;
		
		x.resize(LV.Dimension);
		for (int d = 0; d < LV.Dimension; d++) {
			x[d] = kv::interval< double >(0, 0.5);
		}

		
		uh_min = LV.Generator.global_min(x, std::pow(2.0, -9));
		uh_max = LV.Generator.global_max(x, std::pow(2.0, -9));
		uh_max_min.zeros(LV.Number_of_variables, 1);
		for (int i = 0; i < LV.Number_of_variables; i++) {
			std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
			uh_max_min(i) = VData(uh_min[i], uh_max[i]);
		}
		LV.clear();
		vcp::time.toc();
	}

/////////////////////////////////////////////////////////////////////////////////////////////////
/************************** Compute C1, C2 and C3 in Lemma 1 and K *****************************/
/////////////////////////////////////////////////////////////////////////////////////////////////

	vcp::time.tic();
	VData Ch, Cs2, Cs4;
	VData CfdX = VData(0);
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
		Lotka_Volterra2< VData, VPOLICY > LV;
		LV.first_execute();
		mode = 1;
		LV.Create_Generator(mode, list_type);
		vcp::matrix< VData, VPOLICY > uhi_rough;
		vcp::convert(uhi, uhi_rough);
		LV.Generator.setting_uh(uhi_rough, list_uh, 2);
		LV.setting_uh(uhi_rough);


		vcp::matrix< VData, VPOLICY > DF, L, E, Q, QT;
		{
			vcp::matrix< VData, VPOLICY >  Q, QT;	
			D = LV.make_D();
			Q = LV.make_N();
			QT = transpose(Q);
			DF = D - (Q + QT);
			DF += Q*lss(D,QT);
		}
		L = LV.make_L();

		// Compute some constants 
		Cs2 = LV.Generator.Poincare_constant< VData >();
		std::cout << "Cs2 = " << Cs2 << std::endl;

		Cs4 = LV.Generator.Sobolev_constant< VData >(4);
		std::cout << "Cs4 = " << Cs4 << ", p =" << "4" << std::endl;

		Ch = LV.Generator.Ritz_projection_error< VData >();
		std::cout << "Ch = " << Ch << std::endl;
		LV.clear();

		std::cout << "T is invertible!" << std::endl;

		compsym(DF);
		compsym(L);
		eigsymge(DF, L, E);

		DF.clear();
		L.clear();

		E = abs(diag(E));
		tau = sqrt((1/min(E))(0));
		std::cout << "tau = " << tau << std::endl;

		Q = intervalmag(make_Q(uh_max_min));
        Q = normtwo(Q);
		CfdX = Q(0);
        std::cout << "CfdX = " << std::endl;
        std::cout << CfdX << std::endl;

		C1 = Ch*tau*CfdX;
		C2 = Ch*Cs2*CfdX;
		C3 = Ch*Ch*CfdX;
		std::cout << "CfdX = " << CfdX << std::endl;
		std::cout << "C1 = " << C1 << std::endl;
		std::cout << "C2 = " << C2 << std::endl;
		std::cout << "C3 = " << C3 << std::endl;

		std::cout << "\nCompute K" << std::endl;
		Q.clear();
		Q.zeros(2,2);
		Q(0, 0) = 2;
		Q(0, 1) = abs(c);
		Q(1, 0) = abs(d);
		Q(1, 1) = 2;
		Q = normtwo(Q);
		K = pow(Cs4,2)*sqrt(pow(tau,2) + pow(Ch, 2))*Q(0);
		std::cout << "K = " << K << std::endl;
	}
	vcp::time.toc();


/////////////////////////////////////////////////////////////////////////////////////////////////
/************************************** Compute delta ******************************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	vcp::time.tic();
	VData Res1, Res2;
	VData Res = VData(0);
	VData delta = VData(0);
	{
		std::cout << "\nCalculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
		Lotka_Volterra2< VResData, VResPOLICY > LV;
		LV.first_execute();
		mode = 2;
		LV.Create_Generator(mode, list_type);
		LV.Generator.setting_uh(uhi);

// 1. First
        // -DL * uh - a*uhphi + uh2phi + c*uhvhphi        
		VResData LuhLuh = LV.Generator.integral_LuhLuh(0); // (DL uh, DL uh)
		VResData Luh_uh = LV.Generator.integral_Luhuh(0, 1, 0); // (DL uh, uhphi)
		VResData Luh_uh2 = LV.Generator.integral_Luhuh(0, 2, 0); // (DL uh, uh2phi)
		VResData Luh_uhvh = LV.Generator.integral_Luhuh(0, 1, 1); // (DL uh, uhvhphi)                
		VResData uh2 = LV.Generator.integral_uh(2, 0); // (uhphi, uhphi)        
		VResData uh3 = LV.Generator.integral_uh(3, 0); // (uhphi, uh2phi)
		VResData uh2vh = LV.Generator.integral_uh(2, 1); // (uhphi, uhvhphi)
		VResData uh4 = LV.Generator.integral_uh(4, 0); // (uh2phi, uh2phi)
		VResData uh3vh1 = LV.Generator.integral_uh(3, 1); // (uh2phi, uhvhphi)    
        VResData uh2vh2 = LV.Generator.integral_uh(2, 2); // (uhvhphi, uhvhphi)                
        {
			using std::sqrt;
			using std::abs;
			vcp::convert(abs(LuhLuh - 2 * ( -a*Luh_uh + Luh_uh2 + c*Luh_uhvh + a*uh3 + a*(c*uh2vh) - c*uh3vh1) + a*(a * uh2) + uh4 + c*(c*uh2vh2)), Res1);
			std::cout << "First Residual Norm : || Laplace(uh) - f1(uh, vh) ||_L2^2 <= " << Res1 << std::endl;
		}

        // 1. Second
        // -DL * vh - b*vhphi + d*uhvhphi + vh2phi
        VResData LvhLvh = LV.Generator.integral_LuhLuh(1); // (DL vh, DL vh)
    	VResData Lvh_vh = LV.Generator.integral_Luhuh(1, 0, 1); // (DL vh, vhphi)
        VResData Lvh_uhvh = LV.Generator.integral_Luhuh(1, 1, 1); // (DL vh, uhvhphi)
        VResData Lvh_vh2 = LV.Generator.integral_Luhuh(1, 0, 2); // (DL vh, vh2phi)
        VResData vh2 = LV.Generator.integral_uh(0, 2); // (vhphi, vhphi)
        VResData uhvh2 = LV.Generator.integral_uh(1, 2); // (vhphi, uhvhphi)
        VResData vh3 = LV.Generator.integral_uh(0, 3); // (vhphi, vh2phi)
        uh2vh2 = uh2vh2;
        VResData uhvh3 = LV.Generator.integral_uh(1, 3); // (uhvhphi, vh2phi)
        VResData vh4 = LV.Generator.integral_uh(0, 4); // (vh2phi, vh2phi)
        {
			using std::sqrt;
			using std::abs;
			vcp::convert(abs(LvhLvh + 2 * ( b*Lvh_vh - d*Lvh_uhvh - Lvh_vh2 - b*(d*uhvh2) - b*vh3 + d*uhvh3 ) + b*(b * vh2) + d*(d*uh2vh2) + vh4 ), Res2);
			std::cout << "Second Residual Norm : || Laplace(vh) - f2(uh, vh) ||_L2^2 <= " << Res2 << std::endl;
		}
        {
            using std::sqrt;
		    Res = sqrt(Res1 + Res2);
		    std::cout << "Residual Norm : || F(uh, vh) ||_X <= " << Res << std::endl;
        }
		delta = Ch*Res;
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
