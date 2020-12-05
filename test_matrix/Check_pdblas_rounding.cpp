#include <iostream>
#include <kv/hwround.hpp>
#include<vcp/pdblas.hpp>
#include <vcp/matrix.hpp>

int main(void) {

	kv::hwround::roundnear();

	int n = 300;
	vcp::matrix< double, vcp::pdblas> A, B, CU, CD;
//	vcp::matrix< double > A, B, CU, CD;


	for (int i = 0; i < n; i++) {
		std::cout << i << std::endl;
		kv::hwround::roundnear();
		A.zeros(n, n);
		B.zeros(n, n);
		for (int j = 0; j < n; j++) {
			A(j, i) = 1.0/3.0;
			B(i, j) = 1.0/3.0;
		}
		kv::hwround::roundup();
		CU = A*B;
		kv::hwround::rounddown();
		CD = A*B;
		
		for (int k1 = 0; k1 < n; k1++) {
			for (int k2 = 0; k2 < n; k2++) {
				if (CU(k1, k2) == CD(k1, k2)) {
					std::cout << "No Rounding..." << std::endl;
					kv::hwround::roundnear();
					exit(1);
				}
			}
		}
	}

	kv::hwround::roundnear();

	for (int i = 2; i < n; i++) {
		std::cout << i << std::endl;
		kv::hwround::roundnear();
		A.zeros(n, n);
		B.zeros(n, n);

		for (int j = 0; j < n; j++) {
			A(j, 1) = 1.0;
			B(1, j) = 1.0;
		}
		for (int j = 0; j < n; j++) {
			A(j, i) = std::pow(2.0,-27);
			B(i, j) = std::pow(2.0, -27);
		}

		kv::hwround::roundup();
		CU = A*B;
		kv::hwround::rounddown();
		CD = A*B;

		for (int k1 = 0; k1 < n; k1++) {
			for (int k2 = 0; k2 < n; k2++) {
				if (CU(k1, k2) == CD(k1, k2)) {
					std::cout << "No Rounding..." << std::endl;
					kv::hwround::roundnear();
					exit(1);
				}
			}
		}


	}

	std::cout << "This BLAS can be changed rounding mode!" << std::endl;
	kv::hwround::roundnear();

}
