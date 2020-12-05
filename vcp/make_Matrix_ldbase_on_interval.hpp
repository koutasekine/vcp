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

#ifndef VCP_MAKE_MARRIX_LDBASE_ON_INTERVAL_HPP
#define VCP_MAKE_MARRIX_LDBASE_ON_INTERVAL_HPP

namespace vcp {
#if (defined(RDOUBLE_HPP) || (defined(DD_HPP) && defined(RDD_HPP))) && defined(RMPFR_HPP) && defined(MPFR_HPP)
	template <int _N, typename _T = kv::dd >
	class make_Matrix_ldbase_on_interval : public Legendre_Bases_Generator< kv::interval< kv::mpfr< _N > >, double, vcp::mats< double > > {
		int Div_number;
		vcp::matrix< kv::interval< _T >, vcp::imats< _T > > p_phi;

	public:	
		void make(const int mesh_order = 9) {
			if ((*this).mode <= 2) {
				std::cout << ">> ERROR : make : This function only use the Approximation mode (Mode > 2)" << std::endl;
				exit(1);
			}
			kv::mpfr< _N > mesh_size = pow(kv::mpfr< _N >(2), -mesh_order);
			(*this).Div_number = std::pow(2, mesh_order);

			int mmax = (*this).phi.size();
			vcp::matrix< kv::interval< kv::mpfr< _N >>, vcp::imats< kv::mpfr< _N > > > p_phi_T;
			
			p_phi_T.zeros(mmax, (*this).Div_number);
#ifdef _OPENMP
#ifndef VCP_LEGENDRE_NOMP
#pragma omp parallel for
#endif
#endif	
			for (int i = 0; i < mmax; i++) {
				kv::interval< kv::mpfr< _N >> interval_j;
				for (int j = 0; j < (*this).Div_number; j++) {
					interval_j.lower() = j * mesh_size;
					interval_j.upper() = (j + 1) * mesh_size;
					LegendreBaseFunctions::LegendreFunc((*this).phi[i], interval_j, p_phi_T(i, j));
				}
			}
			convert(p_phi_T, (*this).p_phi);
		}
		void save(const char* name = "") {
			std::string filename = name;
			filename += "Matrix_ldbase_on_interval_LegenderOrder";
			filename += std::to_string((*this).Order_of_Base);
			filename += "_Div";
			filename += std::to_string((*this).Div_number);
			vcp::save((*this).p_phi, filename);
		}
	};

	template < int _N >
	void make_interval_value_of_ldbase_demo(int Legendre_Order, int mesh_order,  const char* name = "") {
		vcp::make_Matrix_ldbase_on_interval< _N > A;
		A.setting(Legendre_Order, 1, 1, 1, 3);
		A.make(mesh_order);
		A.save(name);
	}
#endif
}
#endif // VCP_MAKE_MARRIX_LDBASE_ON_INTERVAL_HPP