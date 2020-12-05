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

#include <vcp/imats.hpp>
#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main(void) {
	int n = 100;
	vcp::matrix< double > Ad;
	vcp::matrix< kv::dd > Add;
	vcp::matrix< kv::mpfr< 300 > > Am300;
	vcp::matrix< kv::mpfr< 500 > > Am500;

	vcp::matrix< kv::interval< double >, vcp::pidblas > Aid;
	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > Aidd;
	vcp::matrix< kv::interval< kv::mpfr< 300 > >, vcp::imats< kv::mpfr< 300 > > > Aim300;
	vcp::matrix< kv::interval< kv::mpfr< 500 > >, vcp::imats< kv::mpfr< 500 > > > Aim500;

	Ad.rand(n);
	vcp::convert(Ad, Add);
	vcp::convert(Add, Am300);
	std::cout << Am300(0,0) << std::endl;
	vcp::convert(Am300, Am500);
	vcp::convert(Am500, Ad);
	std::cout << Ad(0, 0) << std::endl;

	vcp::convert(Am500, Aim300);
	std::cout << Aim300(0, 0) << std::endl;
	vcp::convert(Aim300, Add);
	std::cout << Ad(0, 0) << std::endl;
	vcp::convert(Aim300, Aidd);
	std::cout << Aidd(0, 0) << std::endl;
	vcp::convert(Aidd, Am500);
	std::cout << Am500(0, 0) << std::endl;

	return 0;
}