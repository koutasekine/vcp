#include <iostream>

#include <omp.h>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/imats.hpp>
#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>


int main(void) {
	 vcp::matrix< bool > A, B, C;
	 //vcp::mbool A, B, C;

	A.allfalse(1, 10);
	A(0, 1) = true;
	B.allfalse(1, 10);
	B(0, 1) = true;
	B(0, 0) = true;

	std::cout << A << std::endl;
	std::cout << B << std::endl;

	C = A && B;
	std::cout << C << std::endl;

	C = !((A || B) || B);
	std::cout << C << std::endl;
	std::cout << !C(0, 1) << std::endl;

	C.resize(5, 11);
	std::cout << C << std::endl;

	C.clear();

	vcp::matrix< kv::interval< double >, vcp::pidblas > AA, BB, CC;
	AA.rand(10);
	BB.rand(10);
	CC.rand(10);
	std::cout << ((AA != BB) || (AA != CC)) << std::endl;;

	C(0, 0).flip();
	C = AA != BB;

	std::cout << all((C || (AA != BB)) && (AA == CC)) << std::endl;
	std::cout << any((C || (AA != BB)) && (AA == CC)) << std::endl;
	std::cout << none((C || (AA != BB)) && (AA == CC)) << std::endl;

	if (none((C || (AA != BB)) && (AA == CC))) {
		std::cout << "nya---" << std::endl;
	}
}