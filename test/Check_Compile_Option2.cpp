#include <iostream>

#include <omp.h>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>

#include <vcp/imats.hpp>
#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main(void) {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
	std::cout << "Using Opem MP for Mats Policy" << std::endl;
#endif
#endif

	bool flag = true;
//	vcp::matrix< kv::interval< double >, vcp::imats< double > > A, B, C0, C3;
//	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > A, B, C0, C3;
    vcp::matrix< kv::interval< double >, vcp::pidblas > A, B, C0, C3;

	load(A, "A");
	load(B, "B");
	load(C0, "C0");
	C3 = A*B;
	for (int i = 0; i < C3.rowsize(); i++) {
		for (int j = 0; j < C3.columnsize(); j++) {
			if (C0(i, j) != C3(i, j)) {
				flag = false;
			}
		}
	}
	if (flag) {
		std::cout << "OK!!" << std::endl;
	}
	else {
		std::cout << "NG..." << std::endl;
	}
	return 0;
}