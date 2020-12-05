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
typedef kv::interval< kv::mpfr< 300 > > DataType;

//typedef vcp::pdblas POLICY;
typedef vcp::mats< AppData > POLICY;

typedef vcp::pidblas VPOLICY;

int main(void){
	std::cout.precision(4);

	vcp::matrix< AppData, POLICY > uh;
	vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;

	int Order_legendre = 20;
	int uh_Order_legendre = 10;
	int p = 2;
	int Dimension = 2;
	int Number_of_variables = 1;

	std::cout << "Dimension = " << Dimension << std::endl;
	std::cout << "uh of Legendre Bases Order = " << uh_Order_legendre << std::endl;
	std::cout << "p (which is maximum order of uh^p)  = " << p << std::endl;
	std::cout << "Number of Variables (e.g. 2 => u and v) = " << Number_of_variables << std::endl;

	// Setting of Approximate_Generator
	std::cout << "Setting the Generator by Approximate mode " << std::endl;
	Approximate_Generator.setting(uh_Order_legendre, p, Dimension, Number_of_variables, 5);

	// Setting the list of Approximate_Generator
	//Approximate_Generator.setting_list();
	Approximate_Generator.setting_evenlist();

	// output the list => list_uh
	vcp::matrix< int > list_uh = Approximate_Generator.output_list();

	// setting initialization value of uh
	uh.ones(list_uh.rowsize(), Number_of_variables);
	//	uh(0) = 30;
	uh = 200 * uh;
	uh(0) = 600;

	std::cout << "Newton Method Staart " << std::endl;

	{
	// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
	vcp::matrix< AppData, POLICY > DL = Approximate_Generator.dphidphi();
	// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
	vcp::matrix< AppData, POLICY > L = Approximate_Generator.phiphi();

	vcp::matrix< AppData, POLICY > uh2phi;
	vcp::matrix< AppData, POLICY > uhphiphi;
	vcp::matrix< AppData, POLICY > DF;
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
			if (cc < pow(2.0,-30)) {
				std::cout << "Convergence \n" << std::endl;
				break;
			}
		}
	}
	}
	// uh data for Grafics
	vcp::matrix< AppData, POLICY > Grafics = Approximate_Generator.output_uh_for_graphics(100);
	std::cout << Grafics << std::endl;

	Approximate_Generator.clear();

	vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
	Verification_Generator.setting(Order_legendre, p, Dimension, Number_of_variables, 1, uh_Order_legendre);
	Verification_Generator.setting_list();

	vcp::matrix< VData, VPOLICY > uhi;
	vcp::convert(uh, uhi);
	// uh setting : Last Argument is list divide : full list => 1 , even list => 2 
	Verification_Generator.setting_uh(uhi, list_uh, 2);

	vcp::matrix< VData, VPOLICY > uh2phi;
	vcp::matrix< VData, VPOLICY > uhphiphi;
	vcp::matrix< VData, VPOLICY > DF;
	vcp::matrix< VData, VPOLICY > F;

	// Make the matrix ((\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
	vcp::matrix< VData, VPOLICY > DL = Verification_Generator.dphidphi();
	// Make the matrix ((phi_i, \phi_j)_{L^2})_{i,j}
	vcp::matrix< VData, VPOLICY > L = Verification_Generator.phiphi();
	uhphiphi = Verification_Generator.uhphiphi(1);
	Verification_Generator.clear();

	DF = DL - 2 * uhphiphi;
	DL.clear();
	uhphiphi.clear();

	vcp::matrix< VData, VPOLICY >  E;
	eigsymge(DF, L, E);

	std::cout << "Approximate minimal eigenvalue for F'[uh] u = lambda u" << std::endl;
	std::cout << min(diag(E)) << std::endl;

	return 0;
}
