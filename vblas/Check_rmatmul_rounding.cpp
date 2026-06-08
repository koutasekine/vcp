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
#include <kv/hwround.hpp>
#include <vcp/matrix.hpp>
#include "rmatmul.hpp"

int main(void) {

	kv::hwround::roundnear();

	int n = 1137;
	vcp::matrix< double > A, B, CU, CD;


	for (int i = 0; i < n; i++) {
//		std::cout << i << std::endl;
		kv::hwround::roundnear();
		A.zeros(n, n);
		B.zeros(n, n);
		for (int j = 0; j < n; j++) {
			A(j, i) = 1.0/3.0;
			B(i, j) = 1.0/3.0;
		}
		CU.zeros(n, n);
		CD.zeros(n, n);
		rmatmul( n, n, n, A.data(), B.data(), CU.data(), 1 );
		rmatmul( n, n, n, A.data(), B.data(), CD.data(), -1 );
		
		for (int k1 = 0; k1 < n; k1++) {
			for (int k2 = 0; k2 < n; k2++) {
			//	std::cout << "check:" << CU(k1, k2) <<  CD(k1, k2) <<std::endl;
				if (CU(k1, k2) == CD(k1, k2)) {
					std::cout << "VBLAS: rmatmul : Cannot change rounding mode..." << std::endl;
					kv::hwround::roundnear();
					exit(1);
				}
			}
		}
	}

	kv::hwround::roundnear();

	for (int i = 2; i < n; i++) {
//		std::cout << i << std::endl;
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
		CU.zeros(n, n);
		CD.zeros(n, n);
		rmatmul( n, n, n, A.data(), B.data(), CU.data(), 1 );
		rmatmul( n, n, n, A.data(), B.data(), CD.data(), -1 );

		for (int k1 = 0; k1 < n; k1++) {
			for (int k2 = 0; k2 < n; k2++) {
			//	std::cout << "check:" << CU(k1, k2) <<  CD(k1, k2) <<std::endl;
				if (CU(k1, k2) == CD(k1, k2)) {
					std::cout << "VBLAS: rmatmul : Cannot change rounding mode..." << std::endl;
					kv::hwround::roundnear();
					exit(1);
				}
			}
		}
	}

	std::cout << "rmatmul can be changed rounding mode!" << std::endl;
	kv::hwround::roundnear();

}
