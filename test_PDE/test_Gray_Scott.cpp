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



#include <cmath>

//typedef double AppData;
typedef kv::dd AppData;
//typedef kv::mpfr<110> AppData;

typedef kv::interval< double > VData;
typedef kv::interval< kv::mpfr< 1000 > > DataType;

//typedef vcp::pdblas POLICY;
typedef vcp::mats< AppData > POLICY;

typedef vcp::pidblas VPOLICY;

int main(void){
	std::cout.precision(17);

	vcp::matrix< AppData, POLICY > uh, uh1, uh2, zh;
	vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;

	int Order_legendre = 20;
	int uh_Order_legendre = 100;
	int p = 3;
	int Dimension = 1;
	int Number_of_variables = 2;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

	// Setting of Approximate_Generator
	std::cout << "Setting the Generator by Approximate mode " << std::endl;
	Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 15);

	// Setting the list of Approximate_Generator
	//Approximate_Generator.setting_list();
	Approximate_Generator.setting_evenlist();

	// output the list => list_uh
	vcp::matrix< int > list_uh = Approximate_Generator.output_list();

	// setting initialization value of uh
//	uh.rand(list_uh.rowsize(), Number_of_variables);
//	uh = 8 * uh;
//	uh(0, 0) = 22;
//	uh(1, 0) = 50;
//	uh(0, 1) = 12;
//	uh(1, 1) = 37;
	
	// vcp::load(uh, "Data_Gray_Scott_/uh_solution1_eps55");
	vcp::load(uh, "Data_Gray_Scott/uh_solution2_eps55");
	uh.resize(list_uh.rowsize(), Number_of_variables);

	std::cout << "Newton Method Start " << std::endl;
	{
	// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
	vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
	// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
	vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

	vcp::matrix< AppData, POLICY > uh2vhphi;
	vcp::matrix< AppData, POLICY > uhphi;
	vcp::matrix< AppData, POLICY > vhphi;
	vcp::matrix< AppData, POLICY > fphi;
	vcp::matrix< AppData, POLICY > uhvhphiphi;
	vcp::matrix< AppData, POLICY > uh2phiphi;
	vcp::matrix< AppData, POLICY > DF;
	vcp::matrix< AppData, POLICY > F;
	vcp::matrix< AppData, POLICY > syuusei;
	vcp::matrix< AppData, POLICY > check;

	uh1 = uh.submatrix({ 0, uh.rowsize() - 1 }, { 0 });
	uh2 = uh.submatrix({ 0, uh.rowsize() - 1 }, { 1 });
	zh = vercat(uh1, uh2); // zh = [u; v]
		
	{
		int epsilon = 55;
		
		AppData cc;
		while(1){
			uh = horzcat(uh1, uh2); // uh = [u, v]

			Approximate_Generator.setting_uh(uh);
			uh2vhphi = Approximate_Generator.uhphi(2,1); // vector (u^2 v, phi)
			uhphi = Approximate_Generator.uhphi(1, 0); // vector (u, phi)
			vhphi = Approximate_Generator.uhphi(0, 1); // vector (v, phi)
			fphi = Approximate_Generator.fphi(); // vector (f, phi)
			uh2phiphi = Approximate_Generator.uhphiphi(2,0); // matrix (u^2 phi, phi)
			uhvhphiphi = Approximate_Generator.uhphiphi(1, 1); // matrix (u v phi, phi)

			DF = vercat(horzcat(DL - epsilon*(2 * uhvhphiphi - L), -epsilon*uh2phiphi), horzcat(2 * uhvhphiphi, DL - (-uh2phiphi - 10 * L)));
			F = vercat(DL * uh1 - epsilon*(uh2vhphi - uhphi), DL * uh2 - (-uh2vhphi + 10*(fphi - vhphi)));
			syuusei = lss(DF, F);
			zh = zh - syuusei;
			check = max(abs(syuusei));
			cc = check(0);
			std::cout << cc << std::endl;
			if (cc < pow(2.0,-90)) {
				std::cout << "Convergence \n" << std::endl;
				uh1 = zh.submatrix({ 0,uh.rowsize() - 1}, { 0 });
				uh2 = zh.submatrix({ uh.rowsize(), zh.rowsize() - 1 }, { 0 });
				uh = horzcat(uh1, uh2);
				Approximate_Generator.setting_uh(uh);
				break;
			}
			uh1 = zh.submatrix({ 0, uh.rowsize() - 1 }, { 0 });
			uh2 = zh.submatrix({ uh.rowsize(), zh.rowsize() -1 }, { 0 });
		}
	}
	}
	// uh data for Grafics
	std::cout << "list_uh = " << std::endl;
	std::cout << list_uh << std::endl;
//	vcp::save(list_uh, "Data_Gray_Scott/list_uh_solution3");

	std::cout << "uh = " << std::endl;
	std::cout << uh << std::endl;
//	vcp::save(uh, "Data_Gray_Scott/uh_solution3");

	vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
	std::cout << "Grafics = " << std::endl;
	std::cout << Grafics << std::endl;
//	vcp::save(Grafics, "Data_Gray_Scott/uh_Graf_solution3");


// minimal and maximum value of approximate solution uh
	std::vector< kv::interval< double > > x;

	x.resize(Dimension);
	for (int d = 0; d < Dimension; d++) {
		x[d] = kv::interval< double >(0, 1);
	}
	std::vector< double > uh_min = Approximate_Generator.global_min(x, std::pow(2.0, -9));
	std::vector< double > uh_max = Approximate_Generator.global_max(x, std::pow(2.0, -9));

	for (int i = 0; i < Number_of_variables; i++) {
		std::cout << "uh in [" << uh_min[i] << ", " << uh_max[i] << "]" << std::endl;
	}

	return 0;
}
