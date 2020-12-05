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

#include <omp.h>

#include <vcp/pdblas.hpp>
#include <vcp/matrix.hpp>

#include <vcp/vcp_timer.hpp>

/*
#define CPU_NAME "Intel Xeon Gold 6142 * 2"
#define CLOCK 2.1 
#define CORE 44
#define SIMD_BIT 512
#define FMA_FLAG true
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING false
*/

/*
#define CPU_NAME "Intel Core i7-4790"
#define CLOCK 3.6 
#define CORE 4
#define SIMD_BIT 256
#define FMA_FLAG true
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING false
*/

/*
#define CPU_NAME "Intel Core i7-7700"
#define CLOCK 3.6 
#define CORE 4
#define SIMD_BIT 256
#define FMA_FLAG true
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING true
*/

/*
#define CPU_NAME "Intel Core i7-6950X"
#define CLOCK 3.5
#define CORE 10
#define SIMD_BIT 256
#define FMA_FLAG true
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING true
*/

/*
#define CPU_NAME "Intel Xeon Phi Knights Landing 7250"
#define CLOCK 1.5
#define CORE 72
#define SIMD_BIT 512
#define FMA_FLAG true
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING true
*/


/*
#define CPU_NAME "Intel Xeon E5-2687W * 2"
#define CLOCK 3.1
#define CORE 16
#define SIMD_BIT 256
#define FMA_FLAG false
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING false
*/


#define CPU_NAME "Intel Xeon E7-4830 v2 * 4"
#define CLOCK 2.2
#define CORE 40
#define SIMD_BIT 256
#define FMA_FLAG false
#define SIMD_UNIT_NUM 2
#define HYPER_THREADING true

int main(void) {
/*---  Matrix size  ---*/
	int MatrixSize = 20000;

	vcp::matrix< double, vcp::pdblas > A, B, C;

/*---  Make vectors and matrices ---*/
	A.rand(MatrixSize);
	B.rand(MatrixSize);

/*---  Matrix Multiplication ---*/
	vcp::time.tic();
	C = A*B;
	vcp::time.toc();

/*---  Benchmark ---*/
	std::cout << "*******************************" << std::endl;
	std::cout << "Benchmark Start" << std::endl;
	std::cout << "CPU: " <<  CPU_NAME << std::endl;
	std::cout << "Clock: " <<  CLOCK << " [GHz]" << std::endl;
	std::cout << "Core: " <<  CORE << std::endl;
	std::cout << "SIMD: " <<  SIMD_BIT << " [bit]" << std::endl;
	if (FMA_FLAG) {
		std::cout << "FMA: True" << std::endl;
	}
	else{
		std::cout << "FMA: False" << std::endl;
	}
	std::cout << "SIMD Culculator number: " << SIMD_UNIT_NUM << std::endl;

	if (HYPER_THREADING) {
		std::cout << "Hyper threading: ON" << std::endl;
	}
	else{
		std::cout << "Hyper threading: OFF" << std::endl;
	}

	std::cout << "Matrix size: " << MatrixSize << std::endl;


	double calctime = (vcp::time.output())/1000.0; // [sec]
	std::cout << "Calculate Time: "  << calctime << " [sec]" << std::endl;

	if (FMA_FLAG){
		std::cout << "Theoretical value: " << 2* CLOCK * CORE * (SIMD_BIT / 64) * SIMD_UNIT_NUM << " [GFlops]" << std::endl;;
	}
	else {
		std::cout << "Theoretical value: " << CLOCK * CORE * (SIMD_BIT / 64) * SIMD_UNIT_NUM << " [GFlops]" << std::endl;;
	}

	int M = MatrixSize/1000;
	std::cout << "Observed value: " << (2*M*M*M)/calctime << " [GFlops]" << std::endl;
	std::cout << "*******************************" << std::endl;

	return 0;
}
