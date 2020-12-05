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

template < typename _T, typename _P, typename _TT = kv::interval< kv::mpfr< 500 > > >
struct HOGE : public vcp::Newton< _T, _P > {
    int Order_legendre;
	int uh_Order_legendre;
	int p;
	int Dimension;
	int Number_of_variables;
    vcp::Legendre_Bases_Generator< _TT, _T, _P > Generator;
    vcp::matrix< _T, _P > DL;
    vcp::matrix< _T, _P > L;
    vcp::matrix< _T, _P > u;

    void setting_newton( vcp::matrix< _T, _P >& uh ) override {
        Generator.setting_uh(uh);
        u = uh;
    }
    vcp::matrix< _T, _P > f() override {
        vcp::matrix< _T, _P > uh2phi = Generator.uhphi(2);
        return DL*u - uh2phi;
    }
    vcp::matrix< _T, _P > Df() override {
        vcp::matrix< _T, _P > uhphiphi = Generator.uhphiphi(1);
        return DL - 2*uhphiphi;
    }

    void first_execute(int mode){
        Order_legendre = 30;
	    uh_Order_legendre = 20;
	    p = 2;
	    Dimension = 2;
	    Number_of_variables = 1;
        Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, mode);
        Generator.setting_evenlist();
        DL = Generator.dphidphi();
        L = Generator.phiphi();
    }
};

int main(void){
//    HOGE< double, vcp::pdblas > N;
//    vcp::matrix< double, vcp::pdblas > uh;
    HOGE< kv::dd, vcp::mats< kv::dd > > N;
    N.first_execute(5);
    vcp::matrix< kv::dd, vcp::mats< kv::dd > > uh;
    

    vcp::matrix< int > list_uh = N.Generator.output_list();
	uh.ones(list_uh.rowsize(), N.Number_of_variables);
	uh = 200 * uh;
	uh(0) = 600;
    uh = N.solve_nls(uh);
    std::cout << uh << std::endl;

    return 0;
}