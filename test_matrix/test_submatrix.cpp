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


#include <kv/interval.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>

int main(void) {
	vcp::matrix< kv::interval< kv::dd >, vcp::imats< kv::dd > > A, B;
	A.rand(15);
	std::cout << A << std::endl;
	B = A.submatrix({}, {});
	//std::cout << B << std::endl;
	B = A.submatrix({}, {4});
	//std::cout << B << std::endl;
	B = A.submatrix({}, {4,14});
	//std::cout << B << std::endl;
	B = A.submatrix({}, {4,4,14});
	//std::cout << B << std::endl;
	B = A.submatrix({14}, {});
	//std::cout << B << std::endl;
	B = A.submatrix({2}, {5});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2 }, { 5,6 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2 }, { 5,3,11 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 4 }, {});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 14 }, {14});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 4 }, { 5,7 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 8 }, { });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 10 }, {3});
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 10 }, { 3,5 });
	//std::cout << B << std::endl;
	B = A.submatrix({ 2, 3, 10 }, { 3,5,10 });
	//std::cout << B << std::endl;

	std::cout << A.submatrix({1}, {})*A.submatrix({}, {1}) << std::endl;
	std::cout << B[1, 2] << std::endl;
}