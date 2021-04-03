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

#pragma once

#ifndef VCP_HYBRID_APP_LSS_HPP
#define VCP_HYBRID_APP_LSS_HPP

namespace vcp {
	template<typename _T, class _P> vcp::matrix< _T, _P > hybrid_app_lss(const vcp::matrix< _T, _P >& A, vcp::matrix< _T, _P >& b) {
		vcp::matrix< _T, _P > R, x;
		
		{
			vcp::matrix< double, vcp::pdblas > Rd, xd;
			{
				vcp::matrix< double, vcp::pdblas > Ad;
				vcp::matrix< double, vcp::pdblas > I;
				vcp::convert(A, Ad);
				I.eye(Ad.rowsize());
				Rd = lss(Ad, I);
			}
			{
				vcp::matrix< double, vcp::pdblas > bd;
				vcp::convert(b, bd);
				xd = Rd*bd;
			}
			vcp::convert(Rd, R);
			vcp::convert(xd, x);
		}


		_T eps = std::numeric_limits< _T >::epsilon();
		
		_T tmp2;
		int n = 1000;
		int j = 0;

		for(int i = 1; i <= n; i++){
			vcp::matrix< _T, _P> s;
			vcp::matrix< _T, _P > res = A*x-b;
			res = R*res;
			
			s = max(abs(res));
            _T Correction_term = s(0);
            s = max(abs(x));
            Correction_term /= s(0);
			x = x - res;
			std::cout << "i = " << i << " : Convergence condition: " << Correction_term << " < " << _T(256)*std::numeric_limits< _T >::epsilon() << std::endl;
			if (j <= 1){
				if ( j == 0 ){
					tmp2 = Correction_term;
				}
				if ( j == 1){
					std::cout << "Compute convergence rate: " << std::endl;
					_T tmp = pow(_T(2), -52);
					tmp2 = Correction_term/tmp2;
					std::cout << "Convergence rate " <<  tmp2 << std::endl;
					n = 0;
					while (1){
						if (tmp < eps) break;
						tmp *= tmp2;
						n++;
					}
					n++;
					std::cout << "n = " <<  n << std::endl;
				}
				j++;
			}
			
			if( Correction_term <= _T(256)*std::numeric_limits< _T >::epsilon()){
				std::cout << "Convergence: hybrid_app_lss" << std::endl;
				return x;
			}
		}
		return x;
	}
}
#endif