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
#include <cmath>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/imats.hpp>
#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/matrix.hpp>

#include <vcp/doubly_fourier_series.hpp>

int main(void){
    vcp::doubly_fourier_series< double > xs, du1_xs;
    vcp::matrix< double, vcp::pdblas > xvec;
    xs.setting_order(5);
    xs.setting_omega( 1.2, 2.2 );
    xs.b_sin[5] = 2.0;

    std::cout << "---- List ---- " << std::endl;
    std::cout << "p1, p2" << std::endl;
    for (const auto& e : xs.plist) std::cout << e[0] << ", " << e[1] << std::endl;
    std::cout << "-------------- " << std::endl;

    xvec = xs.output_matrix< vcp::pdblas >();
    std::cout << xvec << std::endl;

    du1_xs = xs.diff_u1();

    
    xvec = du1_xs.output_matrix< vcp::pdblas >();
    std::cout << xvec << std::endl;

    std::cout << xs << std::endl;
    std::cout << du1_xs << std::endl;
    std::cout << xs + du1_xs.delay( 2.2 ) << std::endl;

    xs += du1_xs.delay( 2.2 );
    std::cout << -(3.0 - xs) - xs << std::endl;

    return 0;
}