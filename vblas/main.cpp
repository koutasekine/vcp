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

#include <kv/hwround.hpp>

#include <vcp/pdblas.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

#include "udmatmul.hpp"
#include "rmatmul.hpp"

#include <vcp/vcp_timer.hpp>

int main(void){
    int m = 10059;
    int n = 10010;
    int k = 9123;
    vcp::matrix< double, vcp::pdblas > A, B, CU, CD, DU, DD;

    A.rand(m, n);
    B.rand(n, k);
    CU.zeros(m, k);
    CD.zeros(m, k);

    vcp::time.tic();
    udmatmul( m, n, k, A.data(), B.data(), CU.data(), CD.data() );
    vcp::time.toc();

    vcp::time.tic();
    rmatmul( m, n, k, A.data(), B.data(), CU.data(), 1);
    rmatmul( m, n, k, A.data(), B.data(), CD.data(), -1);
    vcp::time.toc();

    //std::cout << A << std::endl;
    //std::cout << B << std::endl;
 


    vcp::time.tic();
    kv::hwround::roundup();
    DU = A*B;
    kv::hwround::rounddown();
    DD = A*B;
    vcp::time.toc();

    std::cout << CU.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;
    std::cout << CD.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;    
    std::cout << DU.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;

/*
    std::cout << CU.submatrix({0,10}, {0,10}) << std::endl;
    std::cout << CD.submatrix({0,10}, {0,10}) << std::endl;
    std::cout << D.submatrix({0,10}, {0,10}) << std::endl;
*/
    //std::cout << D - C << std::endl;

    std::cout << "|| DU - CU || <=" << norminf(DU - CU) << std::endl;
    std::cout << "|| DD - CD || <=" << norminf(DD - CD) << std::endl;
    std::cout << "|| CU - CD || <=" << norminf(CU - CD) << std::endl;
    std::cout << "|| DU - DD || <=" << norminf(DU - DD) << std::endl;
}
