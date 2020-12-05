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
//#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/vcp_metafunction.hpp>

template< class _T > typename std::enable_if<  vcp::is_round_control< _T >::value, void >::type test_round_control_function(_T a) {
	std::cout << "is_round_control : A" << std::endl;
}

template< class _T > typename std::enable_if< !vcp::is_round_control< _T >::value, void >::type test_round_control_function(_T a) {
	std::cout << "is_round_control : B" << std::endl;
}

template< class _T > typename std::enable_if<  vcp::is_interval< _T >::value, void >::type test_interval_function(_T a) {
	std::cout << "is_interval : A" << std::endl;
}

template< class _T > typename std::enable_if< !vcp::is_interval< _T >::value, void >::type test_interval_function(_T a) {
	std::cout << "is_interval : B" << std::endl;
}

template< class _T > typename std::enable_if<  vcp::is_round_interval< _T, _T >::value, void >::type test_round_interval_function(_T a) {
	std::cout << "is_round_interval : A" << std::endl;
}

template< class _T > typename std::enable_if< !vcp::is_round_interval< _T, _T >::value, void >::type test_round_interval_function(_T a) {
	std::cout << "is_round_interval : B" << std::endl;
}


int main(void) {
	std::cout << "Check : is_round_control" << std::endl;
	std::cout << vcp::is_round_control< int >::value << std::endl;
	std::cout << vcp::is_round_control< double >::value << std::endl;
	std::cout << vcp::is_round_control< kv::dd >::value << std::endl;
	std::cout << vcp::is_round_control< kv::mpfr< 200 > >::value << std::endl;

	std::cout << "Check : is_interval" << std::endl;
	std::cout << vcp::is_interval< kv::mpfr< 200 > >::value << std::endl;
	std::cout << vcp::is_interval< kv::interval< double > > ::value << std::endl;
	std::cout << vcp::is_interval< kv::interval<  kv::mpfr< 200 > > > ::value << std::endl;

	std::cout << "Check : is_round_interval" << std::endl;
	std::cout << vcp::is_round_interval< kv::interval<  double > >::value << std::endl;
	std::cout << vcp::is_round_interval< kv::interval< double >, kv::interval< double > > ::value << std::endl;
	std::cout << vcp::is_round_interval< kv::interval<  kv::mpfr< 200 > > > ::value << std::endl;
	std::cout << vcp::is_round_interval< kv::interval<  kv::mpfr< 200 > >, kv::interval<  kv::mpfr< 200 > > > ::value << std::endl;
	

	std::cout << "Check : test_round_control_function" << std::endl;
	int a = 1;
	double b = 1;
	kv::dd c = kv::dd(1);
	kv::mpfr< 300 > d = kv::mpfr< 300 >(1);

	test_round_control_function(a);
	test_round_control_function(b);
	test_round_control_function(c);

	std::cout << "Check : test_interval_function" << std::endl;
	kv::interval< int > ia;
	kv::interval< double > ib;
	kv::interval< kv::dd > ic;
	kv::interval< kv::mpfr< 300 > > id;

	test_interval_function(ia);
	test_interval_function(ib);
	test_interval_function(ic);
	test_interval_function(id);

	std::cout << "Check : test_round_interval_function" << std::endl;
	test_round_interval_function(ia);
	test_round_interval_function(ib);
	test_round_interval_function(ic);
	test_round_interval_function(id);

	return 0;
}