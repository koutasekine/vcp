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

#ifndef VCP_CANDIDATE_SET_HPP
#define VCP_CANDIDATE_SET_HPP

#include <vcp/matrix.hpp>
#include <vcp/ldbase.hpp>
#include <vcp/allsol/constants.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif
namespace vcp {
    template <typename _T, typename _TM, class _PM>
    struct candidate_set {
        vcp::Legendre_Bases_Generator< _T, _TM, _PM > Generator;
        vcp::matrix< _TM, _PM > D;
        vcp::matrix< _TM, _PM > L;
        vcp::matrix< _TM, _PM > norm_phii_L2;        
        vcp::constants< _TM > constant;


        virtual _TM g(const _TM& alpha ){  return _TM(0);   }
        
        template <typename _TT >
        void setting(const int& NN, const int& p, const int& Dimension, const int& Number_of_var){
            this->constant.N = NN;
            this->constant.maximam_p = p;
            this->constant.dim = Dimension;
            this->constant.Number_of_variables = Dimension;

            this->Generator.setting(this->constant.N, this->constant.maximam_p, this->constant.dim, this->constant.Number_of_variables, 1, this->constant.N);
            this->Generator.setting_list();
            this->constant.list_uh = this->Generator.output_list();

            this->constant.CN = this->Generator.template Ritz_projection_error< _TM >();
		    this->constant.CN.lower() = constant.CN.upper();

            this->constant.Cs.resize(this->constant.maximam_p*3+1);
            this->constant.Cs[2] = this->Generator.template Poincare_constant< _TM >();
            this->constant.Cs[2].lower() = this->constant.Cs[2].upper();
            for (int i = 3; i <= this->constant.maximam_p*3; i++ )  {
                this->constant.Cs[i] = this->Generator.template Sobolev_constant< _TM >(i);
                this->constant.Cs[i].lower() = this->constant.Cs[i].upper();
            }
        }

        void setting_alpha(const _TM& alphaa){
            this->constant.alpha = alphaa;
        }

        void calc_g_alpha(void){
            this->constant.g_alpha = this->g( this->constant.alpha );
        }

        void calc_rho_o(void){
            this->calc_g_alpha();
            this->constant.rho_o = this->constant.CN*this->g( this->constant.alpha );
            this->constant.rho_o.lower() = _TM(0).lower();
        }

        vcp::matrix< _TM, _PM > calc_UN(void){
            // Make the matrix (\nabla \phi_i, \nabla \phi_j)_{L^2})_{i,j}
    		this->D = this->Generator.dphidphi();
		    // Make the matrix (\phi_i, \phi_j)_{L^2})_{i,j}
		    this->L = this->Generator.phiphi();
            norm_phii_L2 = sqrt(diag(L));
            this->calc_g_alpha();
            
            norm_phii_L2 = abs(this->constant.g_alpha*norm_phii_L2);          
            for (int i = 0; i < norm_phii_L2.rowsize(); i++) norm_phii_L2(i).lower() = -norm_phii_L2(i).upper();

            return lss(D, norm_phii_L2);
        }
        
    };
}
#endif
