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

#include <vcp/pidblas.hpp>
#include <vcp/pidblas_fma.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

#include <vcp/vcp_timer.hpp>

int main(void){
    int n = 2000;
    vcp::matrix< kv::interval< double >, vcp::pidblas > ABLAS, BBLAS, CBLAS;
    vcp::matrix< kv::interval< double >, vcp::pidblas_fma > AFMA, BFMA, CFMA;
    vcp::matrix< kv::interval< double > > AFOR, BFOR, CFOR;

    ABLAS.rand(n, n);
    BBLAS.rand(n, n);

    vcp::time.tic();
    CBLAS = ABLAS*BBLAS;
    vcp::time.toc();

    vcp::convert(ABLAS, AFMA);
    vcp::convert(BBLAS, BFMA);

    vcp::time.tic();
    CFMA = AFMA*BFMA;
    vcp::time.toc();

    vcp::convert(AFMA, AFOR);
    vcp::convert(BFMA, BFOR);

    vcp::time.tic();
    CFOR = AFOR*BFOR;
    vcp::time.toc();    
    
}