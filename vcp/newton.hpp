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

#ifndef VCP_NEWTON_HPP
#define VCP_NEWTON_HPP

#include <limits>
#include <vcp/vcp_metafunction.hpp>

namespace vcp {
    template < typename _T, typename _PM > class Newton {
        protected:
        bool mode; // mode: False => f and Df, True => f and auto diff
        bool flag_Convergence;
        int newton_max_iteration;
        int iteration_newton;
        _T newton_tol;
        _T Correction_term;
        

        public:
        Newton< _T, _PM >() {
			(*this).newton_max_iteration = 50;
            (*this).flag_Convergence = false;
            (*this).newton_tol = _T(4) * std::numeric_limits< _T >::epsilon();
		}
        virtual void setting_newton( vcp::matrix< _T, _PM >& xx ){
        }

        virtual vcp::matrix< _T, _PM > f() {
            vcp::matrix< _T, _PM> x;
            x.zeros(1);
            return x;
        }
        virtual vcp::matrix< _T, _PM > Df() {
            vcp::matrix< _T, _PM> x;
            x.zeros(1);
            return x;
        }

        vcp::matrix< _T, _PM> solve_nls(const vcp::matrix< _T, _PM>& uh){
            vcp::matrix< _T, _PM> x_old, x_new, fdi_fx;
            x_new = uh;
            (*this).flag_Convergence = false;
            using std::abs;
            for (iteration_newton = 0; iteration_newton < (*this).newton_max_iteration; iteration_newton++){
                x_old = x_new;

                (*this).setting_newton(x_new);
                vcp::matrix< _T, _PM > fuh = (*this).f();

                fdi_fx = lss((*this).Df(), fuh);
                x_new = x_old - fdi_fx; 
                vcp::matrix< _T, _PM> s = max(abs(fdi_fx));
                Correction_term = s(0);
                s = max(abs(x_old));
                Correction_term /= s(0);
                if( Correction_term <= (*this).newton_tol){
                    (*this).flag_Convergence = true;
                    (*this).disp_convergence();
                    return x_new;            
                }
                (*this).disp_continue();
            }
            (*this).flag_Convergence = false;
            std::cout << "Not Convergence : i = " << (*this).newton_max_iteration << ", " << Correction_term << " <= " << (*this).newton_tol << std::endl;
            return x_new;
        }

        virtual void disp_convergence(){
            // std::cout << "Convergence : i = " << iteration_newton << ", " << Correction_term << " <= " << (*this).newton_tol << std::endl;
        }

        virtual void disp_continue(){
            std::cout << "i = " << iteration_newton << ", " << Correction_term << " <= " << (*this).newton_tol << std::endl;
        }

    };
}

#endif // VCP_NEWTON_HPP