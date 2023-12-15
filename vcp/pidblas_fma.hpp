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

#ifndef VCP_PIDBLAS_FMA_HPP
#define VCP_PIDBLAS_FMA_HPP

#include <kv/rdouble.hpp>
#include <kv/interval.hpp>

#include <vcp/pidblas.hpp>
#include <vcp/ifma.hpp>


namespace vcp {
    class pidblas_fma : public pidblas {
        public:
        //C = IA * IB
		void mulmm(const vcp::mats< kv::interval< double > >& B, vcp::mats< kv::interval< double > >& c)const override {
			if (type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c = B;
				for (int i = 0; i < B.n; i++) {
					c.v[i] *= this->v[0];
				}
				return;
			}
			else if ((type == 'C' || type == 'R' || type == 'M') && B.type == 'S') {
				c = *this;
				for (int i = 0; i < this->n; i++) {
					c.v[i] *= B.v[0];
				}
				return;
			}

            c.zeros(this->row, B.column);
#ifdef _OPENMP
#ifndef VCP_MATS_NOMP
				#pragma omp parallel for
#endif
#endif            
			for (int k = 0; k < B.column; k++){
				for (int j = 0; j < this->column; j++){
					for (int i = 0; i < this->row; i++){
							ifma(this->v[i+this->row*j].upper(), this->v[i+this->row*j].lower(), B.v[j+B.row*k].upper(), B.v[j+B.row*k].lower(), c.v[i+this->row*k].upper(), c.v[i+this->row*k].lower(), c.v[i+this->row*k].upper(), c.v[i+this->row*k].lower());
						}
					}
				}
			}


    };
}

#endif