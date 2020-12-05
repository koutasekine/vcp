#include <iostream>
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>

int main(void) {

	kv::hwround::roundnear();

	int n = 100;
	vcp::pdblas A, B, CU, CD;

	for (int i = 0; i < n; i++) {
		std::cout << i << std::endl;
		kv::hwround::roundnear();
		A.zeros(n, n);
		B.zeros(n, n);
		for (int j = 0; j < n; j++) {
			A.v[j+A.row*i] = 1.0/3.0;
			B.v[i+A.row*j] = 1.0/3.0;
		}
		kv::hwround::roundup();
		A.mulmm(B, CU);
		kv::hwround::rounddown();
		A.mulmm(B, CD);
		
		for (int k1 = 0; k1 < n; k1++) {
			for (int k2 = 0; k2 < n; k2++) {
				if (CU.v[k1 + A.row*k2] <= CD.v[k1 + A.row*k2]) {
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
			A.v[j + A.row*1] = 1.0;
			B.v[1 + A.row*j] = 1.0;
		}
		for (int j = 0; j < n; j++) {
			A.v[j + A.row*i] = std::pow(2.0,-27);
			B.v[i + A.row*j] = std::pow(2.0, -27);
		}

		kv::hwround::roundup();
		A.mulmm(B, CU);
		kv::hwround::rounddown();
		A.mulmm(B, CD);

		for (int k1 = 0; k1 < n; k1++) {
			for (int k2 = 0; k2 < n; k2++) {
				if (CU.v[k1 + A.row*k2] <= CD.v[k1 + A.row*k2]) {
					std::cout << "No Rounding..." << std::endl;
					kv::hwround::roundnear();
					exit(1);
				}
			}
		}
	}

	kv::hwround::roundnear();

}
