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
	//	vcp::matrix< double > A, B, E, D;

	/*---  2.Approximate Data type(dd or mpfr precision using KV)  ---*/
	//	vcp::matrix< kv::dd > A, B, E, D, X, G, I;
	//	vcp::matrix< kv::mpfr< 500 > > A, B, E, D;

	/*---  3.Approximate Data type with BLAS and Lapack  ---*/
	//	vcp::matrix< double, vcp::pdblas > A, B, E, D;

	/*---  4.Verification Data type with KV  ---*/
	//	vcp::matrix< kv::interval< double >, vcp::imats< double > > A, B, E, D;
	//	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > A, B, E, D;
	//	vcp::matrix< kv::interval< kv::mpfr< 300 > >, vcp::imats< kv::mpfr< 300 > > > A, B, E, D;

	/*---  5.Verification Data type with KV (Approximate term is used pdblas which is not necessary to chenge the rounding mode on BLAS)  ---*/
	//	vcp::matrix< kv::interval< double >, vcp::imats< double, vcp::pdblas > > A, B, E, D;

	/*---  6.Verification Data type with KV, BLAS and Lapack (Please check to chenge the rounding mode on BLAS)  ---*/
	vcp::matrix< kv::interval< double >, vcp::pidblas > A, B, E, D;

	bool flag = false;

	// OpenMP Check addmm
	A.rand(n, n);
	B.rand(n, n);
	E = A + B;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			B(i, j) = A(i, j) + B(i, j);
			if (E(i, j) != B(i, j)) {
				flag = true;
			}
		}
	}
	if (flag) {
		std::cout << "addmm : error...." << std::endl;
	}
	else {
		std::cout << "addmm : OK!!" << std::endl;
	}

	// OpenMP Check mulmm M
	flag = false;
	A.rand(n, n);
	B.rand(n, n);
	D.zeros(n, n);
	E = A * B;
	for (int k = 0; k < n; k++) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < n; i++) {
				D(i, k) += A(i, j) * B(j, k);
			}
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (E(i, j) != D(i, j)) {
				flag = true;
			}
		}
	}
	if (flag) {
		std::cout << "mulmm M : error...." << std::endl;
	}
	else {
		std::cout << "mulmm M : OK!!" << std::endl;
	}

	// OpenMP Check subsmmA
	flag = false;
	A.rand(n, n);
	B.rand(n, n);
	E = A - B;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			B(i, j) = A(i, j) - B(i, j);
			if (E(i, j) != B(i, j)) {
				flag = true;
			}
		}
	}
	if (flag) {
		std::cout << "subsmmA : error...." << std::endl;
	}
	else {
		std::cout << "subsmmA : OK!!" << std::endl;
	}

	// OpenMP Check math functions
	flag = false;
	A.rand(n, n);
	B.rand(n, n);
	E = exp(pow(abs(log(sqrt(abs(cos(abs(sin(A))))))), B));
	{
		using std::pow;
		using std::abs;
		using std::sqrt;
		using std::sin;
		using std::cos;
		using std::exp;
		using std::log;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				B(i, j) = exp(pow(abs(log(sqrt(abs(cos(abs(sin(A(i, j)))))))), B(i, j)));
				if (E(i, j) != B(i, j)) {
					flag = true;
				}
			}
		}
	}
	if (flag) {
		std::cout << "math functions : error...." << std::endl;
	}
	else {
		std::cout << "math functions : OK!!" << std::endl;
	}

	// OpenMP Check eye mat
	flag = false;
	A.eye(n);
	B.zeros(n);
	for (int j = 0; j < A.columnsize(); j++) {
		for (int i = 0; i < A.rowsize(); i++) {
			if (i == j) {
				B(i, j) = 1;
			}
			else {
				B(i, j) = 0;
			}
			if (A(i, j) != B(i, j)) {
				flag = true;
			}
		}
	}
	if (flag) {
		std::cout << "eye : error...." << std::endl;
	}
	else {
		std::cout << "eye : OK!!" << std::endl;
	}
}