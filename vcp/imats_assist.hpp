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

#pragma once

#ifndef VCP_IMATS_ASSIST_HPP
#define VCP_IMATS_ASSIST_HPP



namespace vcp {
#ifdef INTERVAL_HPP
	template<typename _T, class _P1, class _P2> void mid(const vcp::matrix< kv::interval< _T >, _P1 >& A, vcp::matrix< _T, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				B(i, j) = mid(A(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> void rad(const vcp::matrix< kv::interval< _T >, _P1 >& A, vcp::matrix< _T, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				B(i, j) = rad(A(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> void midrad(const vcp::matrix< kv::interval< _T >, _P1 >& A, vcp::matrix< _T, _P2 >& B, vcp::matrix< _T, _P2 >& C) {
		B.zeros(A.rowsize(), A.columnsize());
		C.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				midrad(A(i, j), B(i, j), C(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> void interval(const vcp::matrix< _T, _P1 >& A, vcp::matrix< kv::interval< _T >, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				B(i, j) = kv::interval< _T >(A(i, j));
			}
		}
	}
	
	template<typename _T, class _P> typename std::enable_if< vcp::is_interval< _T >::value, void >::type compsym(vcp::matrix< _T, _P >& A) {
		if (A.columnsize() != A.rowsize()) {
			std::cout << "compsym: error " << A.columnsize() << " != " << A.rowsize() << std::endl;
			exit(1);
		}
		
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = i + 1; j < A.columnsize(); j++) {
				if (overlap(A(i, j), A(j, i))) {
					A(i, j) = intersect(A(i, j), A(j, i));
					A(j, i) = A(i, j);
				}
				else {
					std::cout << "compsym: error : not overlap" << std::endl;
					exit(1);
				}
			}
		}
	}

	template<typename _T, class _P1, class _P2> vcp::matrix< kv::interval< _T >, _P1 > operator*(const vcp::matrix< kv::interval< _T >, _P1 >& A, const  vcp::matrix< _T, _P2 >& B) {
		if (A.columnsize() != B.rowsize()) {
			std::cout << "vmatmul:error " << A.columnsize() << " != " << B.rowsize() << std::endl;
			exit(1);
		}

		vcp::matrix< kv::interval< _T >, _P1 > C;
		C.zeros(A.rowsize(), B.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				for (int k = 0; k < B.columnsize(); k++) {
					C(i, k) = A(i, j) * B(j, k);
				}
			}
		}

		return std::move(C);
	}
	template<typename _T, class _P1, class _P2> vcp::matrix< kv::interval< _T >, _P2 > operator*(const vcp::matrix< _T, _P1 >& A, const vcp::matrix< kv::interval< _T >, _P2 >& B) {
		if (A.columnsize() != B.rowsize()) {
			std::cout << "vmatmul:error " << A.columnsize() << " != " << B.rowsize() << std::endl;
			exit(1);
		}
		vcp::matrix< kv::interval< _T >, _P2 > C;
		C.zeros(A.rowsize(), B.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				for (int k = 0; k < B.columnsize(); k++) {
					C(i, k) = A(i, j) * B(j, k);
				}
			}
		}

		return std::move(C);
	}
#endif
}
#endif // VCP_IMATS_ASSIST_HPP