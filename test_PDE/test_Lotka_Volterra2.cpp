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

#include <vcp/newton.hpp>

/*--- Data type for Approximate solution ---*/
//typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

//typedef vcp::pdblas AppPOLICY;
typedef vcp::mats< AppData > AppPOLICY;

typedef kv::interval< kv::mpfr< 800 > > AppDataType;

typedef vcp::mats< AppData > POLICY;

/*--- Data type for Inverse Norm ---*/
typedef kv::interval< AppData > VAccData;
typedef vcp::imats< AppData > VAccPOLICY;

typedef double Data;
typedef kv::interval< Data > VData;
typedef vcp::pidblas VPOLICY;
typedef vcp::pdblas PDBLAS;
//typedef vcp::imats< double, vcp::pdblas > VPOLICY;

/*--- Data type for Residual Norm ---*/
typedef kv::interval< AppData > VResData;
typedef vcp::imats< AppData > VResPOLICY;


int a = 110;
int b = 7;
int c = 100;
int d = -8;

template < typename _T, typename _P, typename _TT = AppDataType >
struct Lotka_Volterra2 : public vcp::Newton< _T, _P > {
    int Inv_Order_legendre;
	int uh_Order_legendre;
	int p;
	int Dimension;
	int Number_of_variables;
    int list_size;
    int N_Goerisch;
    vcp::Legendre_Bases_Generator< _TT, _T, _P > Generator;
    vcp::matrix< _T, _P > DL;
    vcp::matrix< _T, _P > L;
    vcp::matrix< _T, _P > uh1, uh2, uh, z;

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

    void first_execute(){
        (*this).newton_tol = _T(256) * std::numeric_limits< _T >::epsilon();

        uh_Order_legendre = 40;
        Inv_Order_legendre = 20;	    
	    p = 2;
	    Dimension = 2;
	    Number_of_variables = 2;
    }

    void set_Inv_Order_legendre(const int n){
        Inv_Order_legendre = n;
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
        std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
        std::cout << "Inverse Norm for Legendre Bases Order = " << Inv_Order_legendre << std::endl;
        std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
        std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;
    }

    void Approximate_mode(int mode){
        std::cout << "\n>>> Start Approximate mode " << std::endl;
        Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, mode);
        Generator.setting_evenlist();
        DL = Generator.dphidphi();
        L = Generator.phiphi();

        list_size = DL.rowsize();
    }

    void Inverse_mode( const vcp::matrix< VAccData, VAccPOLICY > &uhi, const vcp::matrix< int > &list_uh ){
        std::cout << "\n>>> Start Inverse Norm mode " << std::endl;
        Generator.setting(Inv_Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
        //Generator.setting_evenlist();
        Generator.setting_list();
        
        // uh setting : Last Argument is list divide : full list => 1 , even list => 2 
		Generator.setting_uh(uhi, list_uh, 2);

        DL = Generator.dphidphi();
        L = Generator.phiphi();

        list_size = DL.rowsize();
    }

    void Goerisch_A2_mode( int n, const vcp::matrix< VAccData, VAccPOLICY > &eigenvector, const vcp::matrix< VAccData, VAccPOLICY > &z, const vcp::matrix< VAccData, VAccPOLICY > &cxi, const vcp::matrix< int > &list_cx ){
        std::cout << "\n>>> Goerisch's Method : Make A2" << std::endl;
        // p : cx^(p-1) => f'^(p-1) 
		Generator.setting(Inv_Order_legendre, p, Dimension, 4*n, 0, uh_Order_legendre);
		// Setting the list of Verification_Generator 
		Generator.setting_list();
		//Generator.setting_evenlist();

        list_size = ((*this).Generator.output_list()).rowsize();

        std::cout << "Generator: [ eigvec_u, eigvec_v, z_u, z_v ], z = Rh*Ainv*QT*u" << std::endl;
        // uh1 = [ eigvec_u, eigvec_v, z_u, z_v ]
        uh1 = horzcat(eigenvector.submatrix({ 0, list_size - 1 }, { }), eigenvector.submatrix({ list_size, eigenvector.rowsize() - 1 }, { }), z.submatrix({ 0, list_size - 1 }, { }), z.submatrix({ list_size, eigenvector.rowsize() - 1 }, { }));
        Generator.setting_cx(cxi, list_cx, 2);
        Generator.setting_uh(uh1);

        (*this).N_Goerisch = n;
    }


    vcp::matrix< _T, _P > make_LuhLuh(){
        std::cout << "\n>>> Start: make_LuhLuh" << std::endl;        
        vcp::matrix< _T, _P > LuhLuh1, LuhLuh2;

        LuhLuh1.zeros((*this).N_Goerisch, (*this).N_Goerisch);

        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                LuhLuh1(i, j) = (*this).Generator.integral_LuhLuh(i, j);
//                LuhLuh1(j, i) = LuhLuh1(i, j);
            }
        }

        LuhLuh2.zeros((*this).N_Goerisch, (*this).N_Goerisch);
        for (int i = (*this).N_Goerisch; i < 2*(*this).N_Goerisch; i++){
            for (int j = (*this).N_Goerisch; j < 2 * (*this).N_Goerisch; j++){
                int ii = i - (*this).N_Goerisch;
                int jj = j - (*this).N_Goerisch;
                LuhLuh2(ii, jj) = (*this).Generator.integral_LuhLuh(i, j);
 //               LuhLuh2(jj, ii) = LuhLuh2(ii, jj);
            }
        }

        LuhLuh1 += LuhLuh2;
        compsym(LuhLuh1);
        return LuhLuh1;
    }

    vcp::matrix< _T, _P > make_Luhcxuh(){
        std::cout << "\n>>> Start: make_Luhcxuh" << std::endl;                
        std::cout << "( L*ui, Q*zj ), (L*ui, (Q+QT)*uj), sigma(D ui, D uj)" << std::endl;        
        vcp::matrix< _T, _P > LuhQz, LuhQuh, LuhQTuh;

        // Generator.integral_Luhcxuh(i, j, k, l, ...) => (L ui, cx(k,l,...) uj) 
        std::cout << "1. ( L*ui, Q*zj ) = " << std::endl;
        std::cout << "( L*ui1, (a - 2uh - c*vh)*zj1 + (-c*uh)*zj2 )" << std::endl;
        std::cout << "( L*ui2, (-d*vh)*zj1 + (b - d*uh - 2*vh)*zj2 )" << std::endl;

        LuhQz.zeros((*this).N_Goerisch, (*this).N_Goerisch);

        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                int ui1 = i;
                int ui2 = i + (*this).N_Goerisch;
                int zj1 = j + 2*(*this).N_Goerisch;
                int zj2 = j + 3*(*this).N_Goerisch;
                //  u1 - z1
                LuhQz(i, j) += a*(*this).Generator.integral_Luhcxuh(ui1, zj1, 0, 0); // (L u1, a*z1)
                LuhQz(i, j) += - 2*(*this).Generator.integral_Luhcxuh(ui1, zj1, 1, 0); // (L u1, -2*uh*z1)
                LuhQz(i, j) += - c*(*this).Generator.integral_Luhcxuh(ui1, zj1, 0, 1); // (L u1, -2c*vh*z1)

                LuhQz(i, j) += -c*(*this).Generator.integral_Luhcxuh(ui1, zj2, 1, 0); // (L u1, -c*uh*z2)

                LuhQz(i, j) += -d*(*this).Generator.integral_Luhcxuh(ui2, zj1, 0, 1); // (L u2, -d*vh*z1)

                LuhQz(i, j) += b*(*this).Generator.integral_Luhcxuh(ui2, zj2, 0, 0); // (L u2, b*z2)
                LuhQz(i, j) += -d*(*this).Generator.integral_Luhcxuh(ui2, zj2, 1, 0); // (L u2, -d*uh*z2)                
                LuhQz(i, j) += -2*(*this).Generator.integral_Luhcxuh(ui2, zj2, 0, 1); // (L u2, -2*vh*z2)

//                LuhQz(j, i) = LuhQz(i, j);
            }
        }

        std::cout << "2. ( L*ui, Q*uj ) = " << std::endl;
        std::cout << "( L*ui1, (a - 2uh - c*vh)*uj1 + (-c*uh)*uj2 )" << std::endl;
        std::cout << "( L*ui2, (-d*vh)*uj1 + (b - d*uh - 2*vh)*uj2 )" << std::endl;

        LuhQuh.zeros((*this).N_Goerisch, (*this).N_Goerisch);

        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                int ui1 = i;
                int ui2 = i + (*this).N_Goerisch;
                int uj1 = j;
                int uj2 = j + (*this).N_Goerisch;
                //  u1 - z1
                LuhQuh(i, j) += a*(*this).Generator.integral_Luhcxuh(ui1, uj1, 0, 0); // (L u1, a*u1)
                LuhQuh(i, j) += - 2*(*this).Generator.integral_Luhcxuh(ui1, uj1, 1, 0); // (L u1, -2*uh*u1)
                LuhQuh(i, j) += - c*(*this).Generator.integral_Luhcxuh(ui1, uj1, 0, 1); // (L u1, -2c*vh*u1)

                LuhQuh(i, j) += -c*(*this).Generator.integral_Luhcxuh(ui1, uj2, 1, 0); // (L u1, -c*uh*u2)

                LuhQuh(i, j) += -d*(*this).Generator.integral_Luhcxuh(ui2, uj1, 0, 1); // (L u2, -d*vh*u1)

                LuhQuh(i, j) += b*(*this).Generator.integral_Luhcxuh(ui2, uj2, 0, 0); // (L u2, b*u2)
                LuhQuh(i, j) += -d*(*this).Generator.integral_Luhcxuh(ui2, uj2, 1, 0); // (L u2, -d*uh*u2)                
                LuhQuh(i, j) += -2*(*this).Generator.integral_Luhcxuh(ui2, uj2, 0, 1); // (L u2, -2*vh*u2)

//                LuhQuh(j, i) = LuhQuh(i, j);
            }
        }

        std::cout << "3. ( L*ui, QT*uj ) = " << std::endl;
        std::cout << "( L*ui1, (a - 2uh - c*vh)*uj1 + (-d*vh)*uj2 )" << std::endl;
        std::cout << "( L*ui2, (-c*uh)*uj1 + (b - d*uh - 2*vh)*uj2 )" << std::endl;

        LuhQTuh.zeros((*this).N_Goerisch, (*this).N_Goerisch);

        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                int ui1 = i;
                int ui2 = i + (*this).N_Goerisch;
                int uj1 = j;
                int uj2 = j + (*this).N_Goerisch;
                //  u1 - z1
                LuhQTuh(i, j) += a*(*this).Generator.integral_Luhcxuh(ui1, uj1, 0, 0); // (L u1, a*u1)
                LuhQTuh(i, j) += - 2*(*this).Generator.integral_Luhcxuh(ui1, uj1, 1, 0); // (L u1, -2*uh*u1)
                LuhQTuh(i, j) += - c*(*this).Generator.integral_Luhcxuh(ui1, uj1, 0, 1); // (L u1, -2c*vh*u1)

                LuhQTuh(i, j) += -d*(*this).Generator.integral_Luhcxuh(ui1, uj2, 0, 1); // (L u1, -d*vh*u2)

                LuhQTuh(i, j) += -c*(*this).Generator.integral_Luhcxuh(ui2, uj1, 1, 0); // (L u2, -c*uh*u1)

                LuhQTuh(i, j) += b*(*this).Generator.integral_Luhcxuh(ui2, uj2, 0, 0); // (L u2, b*u2)
                LuhQTuh(i, j) += -d*(*this).Generator.integral_Luhcxuh(ui2, uj2, 1, 0); // (L u2, -d*uh*u2)                
                LuhQTuh(i, j) += -2*(*this).Generator.integral_Luhcxuh(ui2, uj2, 0, 1); // (L u2, -2*vh*u2)

//                LuhQTuh(j, i) = LuhQTuh(i, j);
            }
        }

        return LuhQz - (LuhQuh + LuhQTuh);
    }

    vcp::matrix< _T, _P > make_cxuhuh(){
        std::cout << "\n>>> Start: make_cxuhuh" << std::endl;
        vcp::matrix< _T, _P > QTQzizj, QpQT_T_Q, QpQTT_QpQT;
        

        // Generator.integral_Luhcxuh(i, j, k, l, ...) => (L ui, cx(k,l,...) uj) 
        std::cout << "1. ( QTQ*zi, zj ) = " << std::endl;
        std::cout << "( (a^2 - 4*a*uh - 2*a*c*vh  + 4*uh^2 + (c^2 + d^2)*vh^2 + 4*c*uh*vh )*zi1 + ( - a*c*uh - b*d*vh + 2*c*uh^2 + 2*d*vh^2 + (c^2 + d^2)*uh*vh )*zi2 , zj1)" << std::endl;
        std::cout << "( (- a*c*uh - b*d*vh + 2*c*uh^2 + 2*d*vh^2 + (c^2 + d^2)*uh*vh)*zi1 + (b^2 - 2*b*d*uh - 4*b*vh + (c^2 + d^2)*uh^2 + 4*vh^2 + 4*d*uh*vh )*zi2, zj2 )" << std::endl;

        QTQzizj.zeros((*this).N_Goerisch, (*this).N_Goerisch);
        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                int zi1 = i + 2*(*this).N_Goerisch;
                int zi2 = i + 3*(*this).N_Goerisch;                
                int zj1 = j + 2*(*this).N_Goerisch;
                int zj2 = j + 3*(*this).N_Goerisch;

                //  zi1 - zj1 : (a^2 - 4*a*uh - 2*a*c*vh  + 4*uh^2 + (c^2 + d^2)*vh^2 + 4*c*uh*vh )*zi1, zi2
                QTQzizj(i, j) += pow(a,2)*(*this).Generator.integral_cxuhuh(zi1, zj1, 0, 0); // (a^2 zi1, zj1)
                QTQzizj(i, j) += -4*a*(*this).Generator.integral_cxuhuh(zi1, zj1, 1, 0); // (-4*a*uh zi1, zj1)
                QTQzizj(i, j) += -2*a*c*(*this).Generator.integral_cxuhuh(zi1, zj1, 0, 1); // (-2*a*c*vh zi1, zj1)                
                QTQzizj(i, j) += 4*(*this).Generator.integral_cxuhuh(zi1, zj1, 2, 0); // (4*uh^2 zi1, zj1)
                QTQzizj(i, j) += (pow(c,2) + pow(d,2))*(*this).Generator.integral_cxuhuh(zi1, zj1, 0, 2); // ((c^2+d^2)*vh^2 zi1, zj1)
                QTQzizj(i, j) += 4*c*(*this).Generator.integral_cxuhuh(zi1, zj1, 1, 1); // (4*c*uh*vh zi1, zj1)

                //  zi2 - zj1 : (- a*c*uh - b*d*vh + 2*c*uh^2 + 2*d*vh^2 + (c^2 + d^2)*uh*vh )*zi2, zj1
                QTQzizj(i, j) += -a*c*(*this).Generator.integral_cxuhuh(zi2, zj1, 1, 0); // (-a*c*uh zi2, zj1)
                QTQzizj(i, j) += -b*d*(*this).Generator.integral_cxuhuh(zi2, zj1, 0, 1); // (-b*d*vh zi2, zj1)                
                QTQzizj(i, j) += 2*c*(*this).Generator.integral_cxuhuh(zi2, zj1, 2, 0); // (2*c*uh^2 zi2, zj1)
                QTQzizj(i, j) += 2*d*(*this).Generator.integral_cxuhuh(zi2, zj1, 0, 2); // (2*d*vh^2 zi2, zj1)
                QTQzizj(i, j) += (pow(c,2) + pow(d,2))*(*this).Generator.integral_cxuhuh(zi2, zj1, 1, 1); // ((c^2+d^2)*uh*vh zi2, zj1)

                // zi1 - zj2 : (- a*c*uh - b*d*vh + 2*c*uh^2 + 2*d*vh^2 + (c^2 + d^2)*uh*vh)*zi1, zj2
                QTQzizj(i, j) += -a*c*(*this).Generator.integral_cxuhuh(zi1, zj2, 1, 0); // (-a*c*uh zi1, zj2)
                QTQzizj(i, j) += -b*d*(*this).Generator.integral_cxuhuh(zi1, zj2, 0, 1); // (-b*d*vh zi1, zj2)                
                QTQzizj(i, j) += 2*c*(*this).Generator.integral_cxuhuh(zi1, zj2, 2, 0); // (2*c*uh^2 zi1, zj2)
                QTQzizj(i, j) += 2*d*(*this).Generator.integral_cxuhuh(zi1, zj2, 0, 2); // (2*d*vh^2 zi1, zj2)
                QTQzizj(i, j) += (pow(c,2) + pow(d,2))*(*this).Generator.integral_cxuhuh(zi1, zj2, 1, 1); // ((c^2+d^2)*uh*vh zi1, zj2)

                // zi2 - zj2 : (b^2 - 2*b*d*uh - 4*b*vh + (c^2 + d^2)*uh^2 + 4*vh^2 + 4*d*uh*vh )*zi2, zj2
                QTQzizj(i, j) += pow(b,2)*(*this).Generator.integral_cxuhuh(zi2, zj2, 0, 0); // (b^2 zi2, zj2)
                QTQzizj(i, j) += -2*b*d*(*this).Generator.integral_cxuhuh(zi2, zj2, 1, 0); // (-2*b*d*uh zi2, zj2)
                QTQzizj(i, j) += -4*b*(*this).Generator.integral_cxuhuh(zi2, zj2, 0, 1); // (-4*b*vh zi2, zj2)
                QTQzizj(i, j) += (pow(c,2) + pow(d,2))*(*this).Generator.integral_cxuhuh(zi2, zj2, 2, 0); // ((c^2+d^2)*uh^2 zi2, zj2)
                QTQzizj(i, j) += 4*(*this).Generator.integral_cxuhuh(zi2, zj2, 0, 2); // (4*vh^2 zi2, zj2)
                QTQzizj(i, j) += 4*d*(*this).Generator.integral_cxuhuh(zi2, zj2, 1, 1); // (4*d*uh*vh zi2, zj2)
            }
        }

        std::cout << "2. ( (Q+QT)T*Qzi, uj ) = " << std::endl;
        std::cout << "( ( 2*a^2 - 4*a*c*vh - 8*a*uh + 2*c^2*vh^2 + (c*d+8*c)*uh*vh + d^2*vh^2 + 8*uh^2 )*zi1 + ( -(2*a*c+b*c)*uh - b*d*vh + (c*d+4*c)*uh^2 + (2*c^2+2*c+d^2)*uh*vh + 2*d*v^2 )*zi2, uj1)" << std::endl;
        std::cout << "( ( -a*c*uh - a*d*vh - 2*b*d*vh + 2*c*uh^2 + (c*d+4*d)*vh^2 + (2*d^2+c^2+2*d)*uh*vh )*zi1 + ( 2*b^2 - 4*b*d*uh - 8*b*vh + (c^2+2*d^2)*uh^2 + (c*d+8*d)*uh*vh + 8*vh^2 )*zi2,  uj2 )" << std::endl;

        QpQT_T_Q.zeros((*this).N_Goerisch, (*this).N_Goerisch);
        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                int zi1 = i + 2*(*this).N_Goerisch;
                int zi2 = i + 3*(*this).N_Goerisch;                
                int uj1 = j;
                int uj2 = j + (*this).N_Goerisch;

                // zi1 - uj1 : ( 2*a^2 - 8*a*uh - 4*a*c*vh + 8*uh^2 + (2*c^2+d^2)*vh^2 + (c*d+8*c)*uh*vh )*zi1, uj1
                QpQT_T_Q(i, j) += 2*pow(a,2)*(*this).Generator.integral_cxuhuh(zi1, uj1, 0, 0); // (2*a^2 zi1, uj1)
                QpQT_T_Q(i, j) += -8*a*(*this).Generator.integral_cxuhuh(zi1, uj1, 1, 0); // (-8*a*uh zi1, uj1)                
                QpQT_T_Q(i, j) += -4*a*c*(*this).Generator.integral_cxuhuh(zi1, uj1, 0, 1); // (-4*a*c*vh zi1, uj1)
                QpQT_T_Q(i, j) += 8*(*this).Generator.integral_cxuhuh(zi1, uj1, 2, 0); // (8*uh^2 zi1, uj1)
                QpQT_T_Q(i, j) += (2*pow(c,2)+pow(d,2))*(*this).Generator.integral_cxuhuh(zi1, uj1, 0, 2); // ((2*pow(c,2)+pow(d,2))*vh^2 zi1, uj1)
                QpQT_T_Q(i, j) += (c*d+8*c)*(*this).Generator.integral_cxuhuh(zi1, uj1, 1, 1); // ((c*d+8*c)*vh^2 zi1, uj1)

                // zi2 - uj1 : ( -(2*a*c+b*c)*uh + (c*d+4*c)*uh^2 - b*d*vh + 2*d*v^2 + (2*c^2+2*c+d^2)*uh*vh )*zi2, uj1
                QpQT_T_Q(i, j) += -(2*a*c + b*c)*(*this).Generator.integral_cxuhuh(zi2, uj1, 1, 0); // (-(2*a*c+b*c)*uh zi2, uj1)
                QpQT_T_Q(i, j) += -b*d*(*this).Generator.integral_cxuhuh(zi2, uj1, 0, 1); // (-b*d*vh zi2, uj1)                
                QpQT_T_Q(i, j) += (c*d+4*c)*(*this).Generator.integral_cxuhuh(zi2, uj1, 2, 0); // ((c*d+4*c)*uh^2 zi2, uj1)
                QpQT_T_Q(i, j) += 2*d*(*this).Generator.integral_cxuhuh(zi2, uj1, 0, 2); // (2*d*vh^2 zi2, uj1)
                QpQT_T_Q(i, j) += (2*pow(c,2) + 2*c + pow(d,2))*(*this).Generator.integral_cxuhuh(zi2, uj1, 1, 1); // ((2*c^2+2*c+d^2)*uh*vh zi2, uj1)

                // zi1 - uj2 : ( -a*c*uh - (a*d+2*b*d)*vh + 2*c*uh^2 + (c*d+4*d)*vh^2 + (2*d^2+c^2+2*d)*uh*vh )*zi1, uj2
                QpQT_T_Q(i, j) += -a*c*(*this).Generator.integral_cxuhuh(zi1, uj2, 1, 0); // (-a*c*uh zi1, uj2)
                QpQT_T_Q(i, j) += -(a*d+2*b*d)*(*this).Generator.integral_cxuhuh(zi1, uj2, 0, 1); // ((a*d+2*b*d)*uh zi1, uj2)
                QpQT_T_Q(i, j) += 2*c*(*this).Generator.integral_cxuhuh(zi1, uj2, 2, 0); // (2*c*uh^2 zi1, uj2)
                QpQT_T_Q(i, j) += (c*d+4*d)*(*this).Generator.integral_cxuhuh(zi1, uj2, 0, 2); // ((c*d+4*d)*vh^2 zi1, uj2)
                QpQT_T_Q(i, j) += (2*pow(d,2) + pow(c,2) + 2*d)*(*this).Generator.integral_cxuhuh(zi1, uj2, 1, 1); // ((2*d^2+c^2+2*d)*uh*vh zi1, uj2)

                // zi2 - uj2 : ( 2*b^2 - 4*b*d*uh - 8*b*vh + (c^2+2*d^2)*uh^2 + 8*vh^2 + (c*d+8*d)*uh*vh )*zi2, uj2
                QpQT_T_Q(i, j) += 2*pow(b,2)*(*this).Generator.integral_cxuhuh(zi2, uj2, 0, 0); // (2*b^2 zi2, uj2)
                QpQT_T_Q(i, j) += - 4*b*d*(*this).Generator.integral_cxuhuh(zi2, uj2, 1, 0); // (-4*b*d*uh zi2, uj2)
                QpQT_T_Q(i, j) += - 8*b*(*this).Generator.integral_cxuhuh(zi2, uj2, 0, 1); // (-8*b*vh zi2, uj2)
                QpQT_T_Q(i, j) += (pow(c,2) + 2*pow(d,2))*(*this).Generator.integral_cxuhuh(zi2, uj2, 2, 0); // ((c^2+2*d^2)*uh^2 zi2, uj2)
                QpQT_T_Q(i, j) += 8*(*this).Generator.integral_cxuhuh(zi2, uj2, 0, 2); // (8*vh^2 zi2, uj2)
                QpQT_T_Q(i, j) += (c*d + 8*d)*(*this).Generator.integral_cxuhuh(zi2, uj2, 1, 1); // ((c*d+8*d)*uh*vh zi2, uj2)                
            }
        }

        std::cout << "3. ( (Q+QT)T*(Q+QT)ui, uj ) = " << std::endl;
        std::cout << "( ( 4*a^2 - 16*a*uh - 8*a*c*vh + (c^2+16)*uh^2 + (4*c^2+d^2)*vh^2 + (2*c*d+16*c)*uh*vh )*ui1 + ( -(2*a*c+2*b*c)*uh - (2*a*d+2*b*d)*vh + (4*c+2*c*d)*uh^2 + (4*d+2*c*d)*vh^2 + (2*c^2 + 2*d^2 + 4*c + 4*d)*uh*vh )*ui2, uj1)" << std::endl;
        std::cout << "( ( -(2*a*c+2*b*c)*uh - (2*a*d+2*b*d)*vh + (4*c+2*c*d)*uh^2 + (4*d+2*c*d)*vh^2 + (2*c^2 + 2*d^2 + 4*c + 4*d)*uh*vh )*ui1 + ( 4*b^2 - 8*b*d*uh - 16*b*vh + (c^2+4*d^2)*uh^2 + (d^2+16)*vh^2 + (2*c*d+16*d)*uh*vh )*ui2,  uj2 )" << std::endl;
        QpQTT_QpQT.zeros((*this).N_Goerisch, (*this).N_Goerisch);
        for (int i = 0; i < (*this).N_Goerisch; i++){
            for (int j = 0; j < (*this).N_Goerisch; j++){
                int ui1 = i;
                int ui2 = i + (*this).N_Goerisch;         
                int uj1 = j;
                int uj2 = j + (*this).N_Goerisch;

                // ui1 - uj1 : ( 4*a^2 - 16*a*uh - 8*a*c*vh + (c^2+16)*uh^2 + (4*c^2+d^2)*vh^2 + (2*c*d+16*c)*uh*vh )*ui1, uj1
                QpQTT_QpQT(i, j) += 4*pow(a,2)*(*this).Generator.integral_cxuhuh(ui1, uj1, 0, 0); // (4*a^2 ui1, uj1)
                QpQTT_QpQT(i, j) += -16*a*(*this).Generator.integral_cxuhuh(ui1, uj1, 1, 0); // (16*a*uh ui1, uj1)
                QpQTT_QpQT(i, j) += -8*a*c*(*this).Generator.integral_cxuhuh(ui1, uj1, 0, 1); // (-8*a*c*vh ui1, uj1)
                QpQTT_QpQT(i, j) += (pow(c,2) + 16)*(*this).Generator.integral_cxuhuh(ui1, uj1, 2, 0); // ((c^2+16)*uh^2 ui1, uj1)
                QpQTT_QpQT(i, j) += (4*pow(c,2) + pow(d,2))*(*this).Generator.integral_cxuhuh(ui1, uj1, 0, 2); // ((4*c^2+d^2)*vh^2 ui1, uj1)
                QpQTT_QpQT(i, j) += (2*c*d + 16*c)*(*this).Generator.integral_cxuhuh(ui1, uj1, 1, 1); // ((2*c*d+16*c)*uh*vh ui1, uj1)

                // ui2 - uj1 : ( -(2*a*c+2*b*c)*uh - (2*a*d+2*b*d)*vh + (4*c+2*c*d)*uh^2 + (4*d+2*c*d)*vh^2 + (2*c^2 + 2*d^2 + 4*c + 4*d)*uh*vh )*ui2, uj1
                QpQTT_QpQT(i, j) += -(2*a*c + 2*b*c)*(*this).Generator.integral_cxuhuh(ui2, uj1, 1, 0); // (-(2*a*c+2*b*c)*uh ui2, uj1)
                QpQTT_QpQT(i, j) += -(2*a*d + 2*b*d)*(*this).Generator.integral_cxuhuh(ui2, uj1, 0, 1); // (-(2*a*d+2*b*d)*vh ui2, uj1)
                QpQTT_QpQT(i, j) += (4*c + 2*c*d)*(*this).Generator.integral_cxuhuh(ui2, uj1, 2, 0); // ((4*c+2*c*d)*uh^2 ui2, uj1)
                QpQTT_QpQT(i, j) += (4*d + 2*c*d)*(*this).Generator.integral_cxuhuh(ui2, uj1, 0, 2); // ((4*d+2*c*d)*uh^2 ui2, uj1)
                QpQTT_QpQT(i, j) += (2*pow(c,2) + 2*pow(d,2) + 4*c + 4*d) *(*this).Generator.integral_cxuhuh(ui2, uj1, 1, 1); // ((2*c^2 + 2*d^2 + 4*c + 4*d)*uh^2 ui2, uj1)

                // ui1 - uj2 : ( -(2*a*c+2*b*c)*uh - (2*a*d+2*b*d)*vh + (4*c+2*c*d)*uh^2 + (4*d+2*c*d)*vh^2 + (2*c^2 + 2*d^2 + 4*c + 4*d)*uh*vh )*ui1, uj2
                QpQTT_QpQT(i, j) += -(2*a*c + 2*b*c)*(*this).Generator.integral_cxuhuh(ui1, uj2, 1, 0); // (-(2*a*c+2*b*c)*uh ui1, uj2)
                QpQTT_QpQT(i, j) += -(2*a*d + 2*b*d)*(*this).Generator.integral_cxuhuh(ui1, uj2, 0, 1); // (-(2*a*d+2*b*d)*vh ui1, uj2)
                QpQTT_QpQT(i, j) += (4*c + 2*c*d)*(*this).Generator.integral_cxuhuh(ui1, uj2, 2, 0); // ((4*c+2*c*d)*uh^2 ui1, uj2)
                QpQTT_QpQT(i, j) += (4*d + 2*c*d)*(*this).Generator.integral_cxuhuh(ui1, uj2, 0, 2); // ((4*d+2*c*d)*uh^2 ui1, uj2)
                QpQTT_QpQT(i, j) += (2*pow(c,2) + 2*pow(d,2) + 4*c + 4*d)*(*this).Generator.integral_cxuhuh(ui1, uj2, 1, 1); // ((2*c^2 + 2*d^2 + 4*c + 4*d)*uh^2 ui1, uj2)

                // ui2 - uj2 : ( 4*b^2 - 8*b*d*uh - 16*b*vh + (c^2+4*d^2)*uh^2 + (d^2+16)*vh^2 + (2*c*d+16*d)*uh*vh )*ui2, uj2
                QpQTT_QpQT(i, j) += 4*pow(b,2)*(*this).Generator.integral_cxuhuh(ui2, uj2, 0, 0); // (4*b^2 ui2, uj2)
                QpQTT_QpQT(i, j) += -8*b*d*(*this).Generator.integral_cxuhuh(ui2, uj2, 1, 0); // (-8*b*d*uh ui2, uj2)
                QpQTT_QpQT(i, j) += -16*b*(*this).Generator.integral_cxuhuh(ui2, uj2, 0, 1); // (-16*b*vh ui2, uj2)
                QpQTT_QpQT(i, j) += (pow(c,2)+4*pow(d,2))*(*this).Generator.integral_cxuhuh(ui2, uj2, 2, 0); // ((c^2+4*d^2)*uh^2 ui2, uj2)
                QpQTT_QpQT(i, j) += (pow(d,2)+16)*(*this).Generator.integral_cxuhuh(ui2, uj2, 0, 2); // ((d^2+16)*vh^2 ui2, uj2)
                QpQTT_QpQT(i, j) += (2*c*d + 16*d)*(*this).Generator.integral_cxuhuh(ui2, uj2, 1, 1); // ((2*c*d+16*d)*uh*vh ui2, uj2)
            }
        }              

        return QTQzizj - (QpQT_T_Q + transpose(QpQT_T_Q)) + QpQTT_QpQT;
    }

    void Residual_mode( const vcp::matrix< VAccData, VAccPOLICY > &uhi ){
        std::cout << "\n>>> Calculate Residual Norm || Laplace(uh) - f(uh) ||_L2" << std::endl;
		Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 2);

		// Setting the list of Verification_Generator 
		//Generator.setting_list();
		Generator.setting_evenlist();
        Generator.setting_uh(uhi);
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

template < typename _T, typename _P> vcp::matrix< _T, _P > dot_mul(const vcp::matrix< _T, _P > &A, const vcp::matrix< _T, _P > &B){
    if ( (A.rowsize() != B.rowsize() ) || (A.columnsize() != B.columnsize() ) ) {
        std::cout << "Error: dot_mul : size error..." << std::endl;
    } 
    vcp::matrix< _T, _P > C;
    C.zeros(A.rowsize(), A.columnsize());

    for (int i = 0; i < A.rowsize(); i++){
        for (int j = 0; j < A.columnsize(); j++){
            C(i,j) = A(i, j) * B(i, j);
        } 
    }

    return C;
}

int main(void){
    vcp::matrix< AppData, AppPOLICY > uh;
    vcp::matrix< VData, VPOLICY > uh_max_min;
    vcp::matrix< int > list_uh;

    {
        // Compute an approximate solution uh
        Lotka_Volterra2< AppData, AppPOLICY > LV;
        LV.first_execute();
        LV.display();
        LV.Approximate_mode(15);

        vcp::matrix< AppData, AppPOLICY > zh;
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

        vcp::matrix< AppData, AppPOLICY > Graphics = LV.Generator.output_uh_for_graphics(100);
        std::cout << "Graphics = " << std::endl;
        std::cout << Graphics << std::endl;

        // minimal and maximum value of approximate solution uh
        
        vcp::time.tic();
        std::vector< VData > x;
        x.resize(LV.Dimension);
        for (int d = 0; d < LV.Dimension; d++) {
            x[d] = VData(0, 1);
        }
        std::vector< Data > uh_min = LV.Generator.global_min(x, std::pow(2.0, -8));
        std::vector< Data > uh_max = LV.Generator.global_max(x, std::pow(2.0, -8));

        uh_max_min.zeros(LV.Number_of_variables, 1);
        for (int i = 0; i < LV.Number_of_variables; i++) {
            std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
            uh_max_min(i) = VData(uh_min[i], uh_max[i]);
        }
        vcp::time.toc();
        
    }

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/******************* Calculate Inverse Norm || F'[uh]^-1 ||_(H-1,H10) <= K *********************/
	/////////////////////////////////////////////////////////////////////////////////////////////////
    vcp::time.tic();
    VData C2, Cs6, CN;
    VData G;
    {
        vcp::matrix< VAccData, VAccPOLICY > uhi;     
        vcp::interval(uh, uhi);

        Lotka_Volterra2< VAccData, VAccPOLICY > LV;
        LV.first_execute();
        LV.Inverse_mode(uhi, list_uh);

        vcp::matrix< VData, VPOLICY > DF, L, E, mu_sh, VG, Linv, D, Qmat;
        {
            vcp::matrix< VData, VPOLICY >  N, NT;
            {
                vcp::matrix< VAccData, VAccPOLICY > AccD = LV.make_D();
                vcp::matrix< VAccData, VAccPOLICY > AccN = LV.make_N();
                vcp::matrix< VAccData, VAccPOLICY > AccNT = transpose(AccN);
                vcp::matrix< VAccData, VAccPOLICY > AccDF = AccD - (AccN + AccNT);
                vcp::convert(AccDF, DF);
                AccDF.clear();
                vcp::convert(AccD, D);
                AccD.clear();
                vcp::convert(AccN, N);
                AccN.clear();
                vcp::convert(AccNT, NT);
                AccNT.clear();
            }
            compsym(D);
            DF = DF + N * lss(D, NT);
            Qmat = N;
        }

        // Calculate some constants 
		C2 = LV.Generator.Poincare_constant< VData >();
		std::cout << "C2 = " << C2 <<  std::endl;

		Cs6 = LV.Generator.Sobolev_constant< VData >(6);
		std::cout << "Cs6 = " << Cs6 << ", p =" << "6" << std::endl;

		CN = LV.Generator.Ritz_projection_error< VData >();
		std::cout << "CN = " << CN << std::endl;
		LV.Generator.clear();

        vcp::convert(LV.make_L(), L);
		compsym(DF); // Forced symmetrization
		compsym(L);

        std::cout << "DF size = " << DF.rowsize() << ", " << DF.columnsize() << std::endl;
		eigsymge(DF, L, E);
        std::cout << "DF size = " << DF.rowsize() << ", " << DF.columnsize() << std::endl;

        E = diag(E);
        std::cout << "approximate eigenvalue = " << std::endl;
		std::cout << min(E) << std::endl;


        vcp::matrix< VData, VPOLICY > Q;
        Q = make_Q(uh_max_min);
        Q = intervalmag(Q + transpose(Q));
		vcp::matrix< VData, VPOLICY > QQ_norm = normtwo(Q);
        std::cout << QQ_norm << std::endl;

        VData sigma = QQ_norm(0) + 0.0001;
        Q = make_Q(uh_max_min);
        Q = -( Q + transpose(Q) );
        for (int i = 0; i < LV.Number_of_variables; i++){
            Q(i,i) += sigma;
        }
        vcp::matrix< VData, VPOLICY > Qsigma_norm = normtwo(Q);
        std::cout << Qsigma_norm << std::endl;

        Q = intervalmag(make_Q(uh_max_min));
        vcp::matrix< VData, VPOLICY > Q_norm = normtwo(Q);
        std::cout << "Qnorm = " << std::endl;
        std::cout << Q_norm << std::endl;

        VData mu_sh1 = E(0) + sigma;

		std::cout << "mu_sh1 = " << std::endl;
		std::cout << mu_sh1 << std::endl;

		VData C_Msigma = CN*sqrt(1 + pow(CN, 2)*( Qsigma_norm(0) + pow(C2*Q_norm(0), 2) ));

		std::cout << "C_Msigma = " << std::endl;
		std::cout << C_Msigma << std::endl;

	//	VData mu_low = (mu_sh1 - pow(CN*Q_norm(0),2))/(1 + pow(C_Msigma,2)*(mu_sh1 + pow(CN*Q_norm(0),2)));
    	VData mu_low = mu_sh1/(1 + pow(C_Msigma,2)*mu_sh1 );
		std::cout << "mu_low = " << std::endl;
		std::cout << mu_low << std::endl;

        VData mu_up = mu_sh1 + pow(CN*Q_norm(0),2);
        std::cout << "mu_up = " << std::endl;
        std::cout << mu_up << std::endl;

        VData lambda_low = mu_low - sigma;
        VData lambda_up = mu_up - sigma;
        VData lambda = lambda_low;
        lambda.upper() = lambda_up.upper();

        std::cout << "lambda* = " << std::endl;
        std::cout << lambda << std::endl;


        vcp::matrix< VData, VPOLICY > Lambda;
		Lambda = E;
        mu_sh = E + sigma;
        for (int i = 0; i < mu_sh.rowsize(); i++){
            mu_low = mu_sh(i)/(1 + pow(C_Msigma,2)*mu_sh(i) );
            mu_sh(i).lower() = mu_low.lower();
			Lambda(i).lower() = (mu_sh(i) - sigma).lower();
        }

        lambda = Lambda(0);
		std::cout << "lambdah* = " << std::endl;
		std::cout << lambda << std::endl;

        if (lambda.lower() > 0){
            using std::sqrt;
            G = 1/sqrt(lambda);

            std::cout << "Pre G = " << std::endl;
            std::cout << G << std::endl;
            sigma = 0;
        }
        else {
            std::cout << "Can not estimate Pre G" << std::endl;
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

        int N = 0;
        VData rho_mush;
        for (int i = mu_sh.rowsize() - 1; i > 0 ; i--){
            if (!overlap(mu_sh(i), mu_sh(i-1))){
                N = i-1; // from 0 to  N-1 
                rho_mush = VData( mu_sh(i).lower() );
                std::cout << "N+1 = " << N+1 << ", rho_mush = " << rho_mush << std::endl;
                break;
            }
        }
        if ( N == 0){
            std::cout << "Could not find rho..." << mu_sh.rowsize() << std::endl;            
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

        Lotka_Volterra2< VAccData, VAccPOLICY > LV2;
        LV2.first_execute();
        LV2.Goerisch_A2_mode(N, u_eigenvector, z_acc, uhi, list_uh);

        ones.ones(E.rowsize(), 1);
        Er = 1/(E + r);
        Ert = (Er*transpose(ones) + ones*transpose(Er))/pow(r,2);
        Er2 = Er*transpose(Er);
        compsym(Er2);
        compsym(Ert);
        
        std::cout << "Er = " << Er << std::endl;
        std::cout << "Ert = " << Ert << std::endl;
        std::cout << "Er2 = " << Er2 << std::endl;        

        std::cout << ">>> Create matrixes  (Laplace ui, Laplace uj)" << std::endl;
        vcp::matrix< VData, VPOLICY > LuLu;
        {
            vcp::matrix< VAccData, VAccPOLICY > LuLu_acc;          
            // (Δui, Δuj)
            LuLu_acc = LV2.make_LuhLuh();
            convert(LuLu_acc, LuLu);
            std::cout << "LuLu = " << LuLu << std::endl;
        }

        std::cout << ">>> Create matrixes  (Laplace*ui, B*uj)" << std::endl;
        std::cout << "(Laplace ui, B uj) = (Laplace*ui, Q*zj) - (Laplace ui, (Q + QT)*uj) + sigma(ui, uj)" << std::endl;
        vcp::matrix< VData, VPOLICY > LuBu;
        {
            vcp::matrix< VAccData, VAccPOLICY > LuBu_acc;          
            // (Δui, Buj)
            LuBu_acc = LV2.make_Luhcxuh();
            convert(LuBu_acc, LuBu);
            LuBu -= sigma*(transpose(U)*D*U);
            std::cout << "LuBu = " << LuBu << std::endl;
        }

        std::cout << ">>> Create matrixes  (B*ui, B*uj)" << std::endl;
        vcp::matrix< VData, VPOLICY > BuBu;
        {
            vcp::matrix< VData, VPOLICY > tmp;
            vcp::matrix< VAccData, VAccPOLICY > BuBu_acc;          
            // (Δui, Buj)
            BuBu_acc = LV2.make_cxuhuh();
            convert(BuBu_acc, BuBu);
            compsym(BuBu);
            tmp = sigma*transpose(U)*Qmat*z;
            tmp += transpose(tmp);
            BuBu += tmp;


            tmp = sigma*transpose(U)*Qmat*U;
            tmp += transpose(tmp);
            tmp = 2*tmp;
            BuBu = BuBu - tmp;
            BuBu += pow(sigma,2)*A1;
            compsym(BuBu);

            vcp::matrix< double, vcp::pdblas > QQQ, EEE;
            convert(BuBu, QQQ);
            eigsym(QQQ, EEE);
            std::cout << "eig BuBu " << diag(EEE) << std::endl;
        }

        A2 = dot_mul(LuLu - (LuBu+transpose(LuBu)) + BuBu, Er2/pow(r,2));
        compsym(A2);
        std::cout << "A2 = " << A2 << std::endl;
        vcp::matrix< double, vcp::pdblas > QQQ, EEE;
        convert(A2, QQQ);
        eigsym(QQQ, EEE);
        std::cout << "eig A2 " << diag(EEE) << std::endl;

        A2 = dot_mul(transpose(U)*DFsigma*U, Er2 - Ert) + A1/pow(r,2) + dot_mul(LuLu - (LuBu + transpose(LuBu)) + BuBu, Er2/pow(r,2));
        compsym(A2);
        std::cout << "A2 = " << A2 << std::endl;

        convert(A2, QQQ);
        eigsym(QQQ, EEE);
        std::cout << "eig A2 " << diag(EEE) << std::endl;

        A = A0 - rho_mush*A1;
        compsym(A);

        convert(A, QQQ);
        eigsym(QQQ, EEE);
        std::cout << "eig A " << diag(EEE) << std::endl;

        B = A0 - 2*rho_mush*A1 + pow(rho_mush, 2)*A2;
        compsym(B);
        convert(B, QQQ);
        eigsym(QQQ, EEE);
        std::cout << "eig B " << diag(EEE) << std::endl;

        eigsymge(A, B, Mu);
        Mu = diag(Mu);

        std::cout << min(rho_mush - rho_mush/(1 - Mu) ) << ", " << min((rho_mush - rho_mush/(1 - Mu)) - sigma - r ) << std::endl;
        auto eiglow = min((rho_mush - rho_mush/(1 - Mu)) - sigma - r );

        if (eiglow(0).lower() > 0){
            using std::sqrt;
            G = 1/sqrt(eiglow(0));

            std::cout << "G = " << std::endl;
            std::cout << G << std::endl;

        }
        else {
            std::cout << "Can not estimate G..." << std::endl;
        }


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
        vcp::matrix< VResData, VResPOLICY > uhi;     
        vcp::interval(uh, uhi);

        Lotka_Volterra2< VResData, VResPOLICY > LV;
        LV.first_execute();
        LV.Residual_mode(uhi);

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
    }
    vcp::time.toc();

    return 0;
}
