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

typedef vcp::pidblas VDP;
typedef kv::interval< double > VD;
typedef vcp::imats< kv::dd > VDDP;
typedef kv::interval< kv::dd > VDD;
typedef kv::interval< kv::mpfr< 500 > > VMPFR;

int main(void){
	std::cout.precision(17);

	int Order_legendre = 30;
	int uh_Order_legendre = Order_legendre;
	int p = 1;
	int Dimension = 1;
	int Number_of_variables = 1;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "Legendre Bases Order = " << Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;


/////////////////////////////////////////////////////////////////////////////////////////////////
/*************************** Find Eigenvalue such that Lu = lambda u ***************************/
/////////////////////////////////////////////////////////////////////////////////////////////////
	{
		vcp::Legendre_Bases_Generator< VMPFR, VDD, VDDP > Verification_Generator;
		Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
		//Verification_Generator.setting_list();
		Verification_Generator.setting_evenlist();

		// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
		vcp::matrix< VD, VDP > DL;
		convert(Verification_Generator.dphidphi(), DL);
		// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
		vcp::matrix< VD, VDP > L;
		convert(Verification_Generator.phiphi(), L);

		VD CN = Verification_Generator.Ritz_projection_error< VD >();
		std::cout << "CN = " << CN << std::endl;

		Verification_Generator.clear();

		vcp::matrix< VD, VDP >  E;
		eigsymge(L, DL, E);

		std::cout << "Approximate eigenvalues" << std::endl;
		std::cout << 1/diag(E) << std::endl;

		DL.clear();
		L.clear();

		vcp::matrix< VD, VDP > Emin;
		Emin = 1 / max(diag(E));
		double Emin_up = Emin(0).upper();


		std::cout << "Approximate minimal eigenvalue for L u = lambda u" << std::endl;
		std::cout << Emin << std::endl;

		// Compute a lower bound of minimal eigenvalue based on Liu-Oishi's Theorem
		using std::pow;
		Emin(0) = Emin(0) / (1 + pow(CN, 2)*Emin(0));

		// Compute a upper bound of minimal eigenvalue based on Rayleigh-Ritz bound
		Emin(0).upper() = Emin_up;

		std::cout << "Verified minimal eigenvalue for L u = lambda u" << std::endl;
		std::cout << Emin << std::endl;
	}

	return 0;
}
