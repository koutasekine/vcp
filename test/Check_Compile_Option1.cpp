#include <iostream>

#include <omp.h>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>

#include <vcp/imats.hpp>
//#include <vcp/pdblas.hpp>
//#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main(void) {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
	std::cout << "Using Opem MP for Mats Policy" << std::endl;
#endif
#endif
	int n = 100;
	vcp::matrix< kv::interval< double >, vcp::imats< double > > A, B, C0;
//	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > A, B, C0;
//  vcp::matrix< kv::interval< double >, vcp::pidblas > A, B, C0;

	A.rand(n);
	B.rand(n);
	save(A, "A");
	save(B, "B");
	C0 = A*B;
	save(C0, "C0");
	return 0;
}