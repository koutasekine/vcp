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
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
	std::cout << "Using Opem MP for Mats Policy" << std::endl;
#endif
#endif
/*---  Matrix size  ---*/
	int n = 100;
/*---  Select Data type  ---*/
/*---  1.Approximate Data type(double precision)  ---*/
	vcp::matrix< double > A, B, E, D, X, G, I;

/*---  2.Approximate Data type(dd or mpfr precision using KV)  ---*/
//	vcp::matrix< kv::dd > A, B, E, D, X, G, I;
//	vcp::matrix< kv::mpfr< 500 > > A, B, E, D, X, G, I;

/*---  3.Approximate Data type with BLAS and Lapack  ---*/
//	vcp::matrix< double, vcp::pdblas > A, B, E, D, X, G, I;

/*---  4.Verification Data type with KV  ---*/
//	vcp::matrix< kv::interval< double >, vcp::imats< double > > A, B, E, D, X, G, I;
//	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > A, B, E, D, X, G, I;
//	vcp::matrix< kv::interval< kv::mpfr< 300 > >, vcp::imats< kv::mpfr< 300 > > > A, B, E, D, X, G, I;

/*---  5.Verification Data type with KV (Approximate term is used pdblas which is not necessary to chenge the rounding mode on BLAS)  ---*/
//	vcp::matrix< kv::interval< double >, vcp::imats< double, vcp::pdblas > > A, B, E, D, X, G, I;

/*---  6.Verification Data type with KV, BLAS and Lapack (Please check to chenge the rounding mode on BLAS)  ---*/
//	vcp::matrix< kv::interval< double >, vcp::pidblas > A, B, E, D, X, G, I;

/*---  Make vectors and matrices ---*/
	// Zero column vector A = (0,0, ... , 0)^T
	A.zeros(n, 1);
	// Zero row vector A = (0,0, ... , 0)
	A.zeros(1, n);
	// Zero n*5 matrix vector A = (0, ... , 0; 0, ... ,0; ...)
	A.zeros(n, 5);
	// Zero n*n matrix vector A = (0, ... , 0; 0, ... ,0; ...)
	A.zeros(n);

	// One column vector A = (1,1, ... , 1)^T
	A.ones(n, 1);
	// One row vector A = (1,1, ... , 1)
	A.ones(1, n);
	// One n*5 matrix vector A = (1, ... , 1; 1, ... ,1; ...)
	A.ones(n, 5);
	// One n*n matrix vector A = (1, ... , 1; 1, ... ,1; ...)
	A.ones(n);

	// Randam :: double precison Mersenne twister
	// Randam column vector A = (?,?, ... , ?)^T
	A.rand(n, 1);
	// Randam row vector A = (?,?, ... , ?)
	A.rand(1, n);
	// Randam n*5 matrix vector A = (?, ... , ?; ?, ... ,?; ...)
	A.rand(n, 5);
	// Randam n*n matrix vector A = (?, ... , ?; ?, ... ,?; ...)
	A.rand(n);

	// n*n Identity matrix
	A.eye(n);

/*
	// file save and load 
	A.rand(10, 30);
	save(A, "test");
	load(B, "test");
	// file io check
	for (int i = 0; i < A.rowsize(); i++) {
		std::cout << i << std::endl;
		for (int j = 0; j < A.columnsize(); j++) {
			if (A(i, j) != B(i, j)) {
				std::cout << "Ahhhhhhhhhhhhh" << std::endl;
			}
		}
	}
*/

	std::cout << "---- Solve to System of Linear equation ----" << std::endl;
	A.rand(n);
	B.ones(n, 3);
	B = A*B;

	X = lss(A, B);
	G = norminf(A*X - B);
	std::cout << G << std::endl;
	std::cout << X(0,0) << std::endl;


	std::cout << "---- Solve to Symmetric Eigenvalue problem ---" << std::endl;
	A.rand(n);
	// ltransmul : A <= transpose(A)*A
	A = ltransmul(A);
	eigsym(A, E);
	std::cout << min(diag(E)) << std::endl;

	A.rand(n);
	A = transpose(A) + A;
	eigsym(A, E);
	std::cout << min(diag(E)) << std::endl;

	std::cout << "---- Solve to Generalized Symmetric Eigenvalue problem ---" << std::endl;
	A.rand(n);
	A = transpose(A) + A;
	B.rand(n);
	// ltransmul : B <= transpose(B)*B
	B = ltransmul(B);
	I.eye(n);
	// B <= B + I
	B += I;
	I.clear();
	eigsymge(A, B, E);
	std::cout << min(diag(E)) << std::endl;

	return 0;
}
