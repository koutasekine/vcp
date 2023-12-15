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

#ifndef VCP_IFMA_HPP
#define VCP_IFMA_HPP



#include <cmath>
#include <cfenv>
#ifdef VCP_USE_AVX512
    #include <immintrin.h>
#endif

namespace vcp {
    #ifdef VCP_USE_AVX512
        constexpr __mmask8 mmask15 = 15;  // 11110000
        constexpr __mmask8 mmask240 = 240;  // 00001111

        inline void ifma( const double& au, const double& ad, const double& bu, const double& bd, const double& cu, const double& cd, double& resup, double& resdown ){
            volatile __m512d c, d, g, res;
            c = _mm512_set_pd(ad, au, ad, au, -ad, -au, -ad, -au);
            d = _mm512_set_pd(bd, bd, bu, bu, bd, bd, bu, bu);
            g = _mm512_set_pd(cu, cu, cu, cu, -cd, -cd, -cd, -cd);
            res = _mm512_fmadd_round_pd(c, d, g, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
            
            resup = _mm512_mask_reduce_max_pd(mmask240, res);
            resdown = -_mm512_mask_reduce_max_pd(mmask15, res);
        }
    #else 
        inline void ifma( const double& au, const double& ad, const double& bu, const double& bd, const double& cu, const double& cd, double& resup, double& resdown ){
            volatile double auv = au;
            volatile double adv = ad;
            volatile double buv = bu;
            volatile double bdv = bd;
            volatile double cuv = cu;
            volatile double cdv = cd;
            volatile double resupv, resdownv;
            volatile double tmp;

            if ( adv >= 0 ){
                if ( bdv >= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(auv, buv, cuv);
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma( adv, bdv, cdv );
                }
                else if ( bdv <= 0 && buv >= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(auv, buv, cuv);            
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma( auv, bdv, cdv );
                }
                else if ( buv <= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(adv, buv, cuv);         
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma( auv, bdv, cdv );     
                }
            }
            else if ( adv <= 0 && auv >= 0 ){
                if ( bdv >= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(auv, buv, cuv);
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma( adv, buv, cdv );            
                }
                else if ( bdv <= 0 && buv >= 0 ){
                    fesetround(FE_UPWARD);
                    tmp = std::fma(auv, buv, cuv);
                    resupv = std::fma(adv, bdv, cuv);
                    if ( resupv < tmp ) resupv = tmp;

                    fesetround(FE_DOWNWARD);
                    tmp = std::fma(auv, bdv, cdv);
                    resdownv = std::fma(adv, buv, cdv);
                    if ( resdownv > tmp ) resdownv = tmp;
                }
                else if ( buv <= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(adv, bdv, cuv);
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma(auv, bdv, cdv);
                }
            }
            else if ( auv <= 0 ){
                if ( bdv >= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(auv, bdv, cuv);
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma(adv, buv, cdv);
                }
                else if ( bdv <= 0 && buv >= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(adv, bdv, cuv);            
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma(adv, buv, cdv);
                }
                else if ( buv <= 0 ){
                    fesetround(FE_UPWARD);
                    resupv = std::fma(adv, bdv, cuv);
                    fesetround(FE_DOWNWARD);
                    resdownv = std::fma(auv, buv, cdv);
                }        
            }

            fesetround(FE_TONEAREST);
            
            resup = resupv;
            resdown = resdownv;
        } 
    #endif
}
#endif // VCP_IFMA_HPP
