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

#include <vcp/implicit_rk.hpp>
#include <vcp/vcp_timer.hpp>

typedef kv::dd T;
typedef vcp::mats< T > P;
typedef kv::interval< kv::mpfr< 500 > > TT;

// Please declare a structure with two template arguments : _T = Datatype, _P = Matrix Policy
template < typename _T, class _P >
struct Parabolic_ODE_FUNC {
	int uh_Order_legendre;
	int p;
	int Dimension;
	int Number_of_variables;
    vcp::Legendre_Bases_Generator< TT, _T, _P > Generator;
    vcp::matrix< _T, _P > MDL;
    vcp::matrix< _T, _P > L;

    // Please define two member functions 'func' and 'diffunc'
    vcp::matrix< _T, _P > func( const _T& t, const vcp::matrix< _T, _P >& uh ) {
        Generator.setting_uh(uh);
        vcp::matrix< _T, _P > uh2phi = Generator.uhphi(2);
        return MDL*uh + uh2phi;
    }
    
    vcp::matrix< _T, _P > diffunc( const _T& t, const  vcp::matrix< _T, _P >& uh ) {
        Generator.setting_uh(uh);
        vcp::matrix< _T, _P > uhphiphi = Generator.uhphiphi(1);
        return MDL + 2*uhphiphi;
    }

    void first_execute(const int mode){
	    uh_Order_legendre = 20;
	    p = 2;
	    Dimension = 2;
	    Number_of_variables = 1;

        Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, mode);
        Generator.setting_evenlist();
        MDL = - Generator.dphidphi();
        L = Generator.phiphi();
    }
};

int main(void){
	std::cout.precision(17);

    vcp::Implicit_RK< Parabolic_ODE_FUNC, T, P > IRK_Generator; // Inheritance : Parabolic_ODE_FUNC< T, P >
    IRK_Generator.first_execute(5);

    IRK_Generator.setting(20); // (Only even value) Runge-Kutte Legendre n th -> 2n order

    vcp::matrix< T, P > uh_initial;
    vcp::matrix< int > list_uh = IRK_Generator.Generator.output_list();
	uh_initial.zeros(list_uh.rowsize(), IRK_Generator.Number_of_variables);
	uh_initial(0) = 1;


    IRK_Generator.initial_value(uh_initial); // Setting initial value
//    IRK_Generator.setting_logger(false);
//    IRK_Generator.setting_global_time(T(0), T(3)); // argument: (1) start_time (defalt:0), (2) end_time (defalt:1000000)
    IRK_Generator.setting_local_time(1/T(2)); // argument: (1) step_size (defalt:1/8)

// Local Solver:
    for (int i = 0; i < 5; i++){
        std::cout << "Step : " << i << std::endl;
        IRK_Generator.local_solve_ode();

        auto cc = IRK_Generator.local_output(); // Type : std::tuple< std::vector< _T >, std::vector< vcp::matrix< _T, _P > > >
        std::cout << "Time = " << std::endl;
        for (auto &tmp : std::get<0>(cc)){
            std::cout << tmp << std::endl;
        }
        std::cout << "Value = " << std::endl;
        for (auto &tmp : std::get<1>(cc)){
            std::cout << tmp << std::endl;
        }

        IRK_Generator.next_step();
    }

/*
    auto data_logger = IRK_Generator.output(); // Type : std::vector< std::tuple< _T , vcp::matrix< _T, _P > > >
    for (auto &one_step : data_logger){
        std::cout << "t = " << std::get<0>(one_step) << std::endl;
        std::cout << "value = " << std::endl;
        std::cout << std::get<1>(one_step) << std::endl;
    }
*/
    return 0;
}