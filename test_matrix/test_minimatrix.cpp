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

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>
#include <kv/psa.hpp>

#include <vcp/matrix.hpp>
#include <vcp/minimatrix.hpp>

int main(void) {
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
	std::cout << "Using Open MP for Mats Policy" << std::endl;
#endif
#endif
/*---  Matrix size  ---*/
	int n = 3;
	// vcp::matrix< double, vcp::minimats<double> > A, B, E, D, X, G, I;
	vcp::matrix< kv::psa< kv::interval< double > > , vcp::minimats< kv::psa< kv::interval< double > > > > A, B, C;

/*---  Make vectors and matrices ---*/
	A.zeros(n,2);
	
	A(0,0).v.resize(2);
	A(0,0).v[0] = 1;
	A(1,0).v.resize(2);
	A(1,0).v[1] = 1;
	A=2*A;
	A = A*transpose(A);

	B.zeros(n,4);
	B(0,0).v.resize(2);
	B(0,0).v[1] = 1;
	B(0,1).v.resize(2);
	B(0,1).v[0] = 1;
	B = B*transpose(B);

	C = A+B;
	std::cout << A + B << std::endl;
	return 0;
}
