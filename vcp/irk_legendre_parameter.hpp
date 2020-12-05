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

#ifndef VCP_IRK_LEGENDRE_PARAMETER_HPP
#define VCP_IRK_LEGENDRE_PARAMETER_HPP

#include <iostream>

#include <stack>
#include <vector>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <kv/psa.hpp>
#include <kv/gamma.hpp>

#include<vcp/vcp_metafunction.hpp>
#include <vcp/vcp_psa_assist.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/matrix.hpp>

#include <vcp/legendre_integral.hpp>

#include<omp.h>

namespace vcp {
    template < typename _T, class _P, typename _TA = kv::interval< kv::mpfr< 500 > > > class IRK_Legendre_Parameter : protected interval_ld_weightpoint< _TA > {
    	protected:
        int order;
        vcp::matrix< _T, _P > alpha, bweight, cpoint;
        
        public:
        void setting_parameter(int n){
            if (n < 3){
                std::cout << "Error: setting_parameter : Please set n > 4";
                exit(1);
            }

            if (n % 2 != 0) {
				n = n + 1;
			}
            order = n;
            std::vector< _TA > PPoint;
            PPoint.resize(n);

            (*this).alpha.zeros(n,n);
            (*this).bweight.zeros(n,1);
            (*this).cpoint.zeros(n,1);

            (*this).interval_ld_weightpoint< _TA >::set(n);
			for (int i = 0; i < n; i++) {
                (*this).Point[i] = ((*this).Point[i] + 1) / 2;
                PPoint[n - ( i + 1)] = (*this).Point[i];
                (*this).Weight[i] = (*this).Weight[i] / 2;
                convert( (*this).Point[i], (*this).cpoint(n - ( i + 1) ));
                convert( (*this).Weight[i], (*this).bweight(i));
			}

            std::vector< kv::psa< _TA > > lagrange;

            lagrange.resize(n);
            for (int i = 0; i < n; i++){
                lagrange[i].v.resize(n);
                lagrange[i].v[0] = _TA(1);
            }
            kv::psa< _TA > tmp;
            _TA stmp;
            tmp.v.resize(n);
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    if (i != j){
                        stmp = PPoint[i] - PPoint[j];
                        tmp.v[0] = - PPoint[j] / stmp;
                        tmp.v[1] = _TA(1) / stmp;
                        lagrange[i] *= tmp;
                    }
                }
            }
            for (int i = 0; i < n; i++){
                lagrange[i]  = integrate(lagrange[i]);
            }
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    convert(eval(lagrange[j], PPoint[i]), (*this).alpha(i, j));
                }
            }
        }

    };
}
#endif // VCP_IRK_LEGENDRE_PARAMETER_HPP