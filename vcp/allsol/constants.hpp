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

#ifndef VCP_CONSTANTS_HPP
#define VCP_CONSTANTS_HPP
namespace vcp {
    template < typename _T >
    struct constants{
        int N;
        int maximam_p; // for ldbase
        int dim; // for ldbase
        int Number_of_variables; // for ldbase
        _T alpha;
        _T g_alpha;
        _T CN;
        _T rho_o;
        _T rho_h;
        _T rho;
        _T uh_minmax;
        _T C1;
        _T C2;
        _T C3;
        _T Cfdw; 
        _T CfdX; 
        _T Cfd_N;
        _T Cfd_O;
        _T delta;
        _T G;
        _T W_N_H10_2; // || uN-uh ||_H10^2
        std::vector< _T > Cs;
        std::vector< _T > uh_Lpnorm;
        vcp::matrix< int > list_uh;

        constants() {
			this->N = -1;
            this->maximam_p = -1;
            this->dim = -1;
            this->alpha = _T(0);
            this->g_alpha = _T(0);
            this->CN = _T(0); 
            this->rho_o = _T(0);
            this->rho_h = _T(0);
            this->uh_minmax = _T(0);
            this->C1 = _T(0);
            this->C2 = _T(0);
            this->C3 = _T(0);
            this->Cfdw = _T(0);
            this->CfdX = _T(0);
            this->Cfd_N = _T(0);
            this->Cfd_O = _T(0);
            this->delta = _T(0);
		}
		~constants() = default;
		constants(const constants&) = default;
		constants(constants&&) = default;
		constants& operator=(const constants& A) = default;
		constants& operator=(constants&& A) = default;
    };
}
#endif
