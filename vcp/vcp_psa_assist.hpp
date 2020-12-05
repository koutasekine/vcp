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

#ifndef VCP_PSA_ASSIST_HPP
#define VCP_PSA_ASSIST_HPP

#ifndef PSA_HPP
#error Please include psa.hpp
#endif

#include<vcp/vcp_metafunction.hpp>

namespace vcp {
	template <typename _T> void psa_value_Horner(const kv::psa< _T >& pol, const _T x, _T& y) {
		int i;
		int n = pol.v.size();

		y = pol.v(n - 1);
		for (i = n - 2; i >= 0; i--) {
			y = y*x + pol.v(i);
		}
	}
	template <typename _T> void poltaylor(const kv::psa < kv::interval< _T > > & phii, kv::psa< kv::interval< _T > >& phit, const kv::interval< _T >& ival) {
		int n = phii.v.size() - 1;
		_T imid = -mid(ival);
		kv::interval< _T > imidi = kv::interval< _T >(imid);

		phit = phii;
		for (int i = 1; i <= n; i++) {
			for (int j = n; j >= i; j--) {
				phit.v[j - 1] = phit.v[j - 1] - imidi*phit.v[j];
			}
		}
	}
	template <typename _T> typename std::enable_if<  vcp::is_interval< _T >::value, void >::type psa_value_Taylor(const kv::psa< _T >& pol, const _T x, _T& y) {
		kv::psa< _T > polt;
		poltaylor(pol, polt, x);
		psa_value_Horner(polt, x - mid(x), y);
	}

	template <typename _T> typename std::enable_if<  vcp::is_interval< _T >::value, void >::type psa_value(const kv::psa< _T >& pol, const _T x, _T& y) {
		psa_value_Taylor(pol, x, y);
	}
	template <typename _T> typename std::enable_if< !vcp::is_interval< _T >::value, void >::type psa_value(const kv::psa< _T >& pol, const _T x, _T& y) {
		psa_value_Horner(pol, x, y);
	}
	template <typename _T> kv::psa< _T > psa_diff(const kv::psa< _T >& fhx){
		int n = fhx.v.size();
		kv::psa< _T > dfhx;
		dfhx.v.resize(n);
		for (int i = n - 1; i > 0; i--) {
			dfhx.v(i - 1) = i*fhx.v(i);
		}
		return dfhx;
	}
	template <typename _T> kv::psa< _T > psa_integral(const kv::psa< _T >& fhx){
		int n = fhx.v.size();
		kv::psa< _T > ifhx;
		ifhx.v.resize(n+1);
		for (int i = 0; i < n ; i++) {
			ifhx.v(i + 1) = fhx.v(i)/(_T(i + 1));
		}
		return ifhx;
	}
}
#endif // VCP_PSA_ASSIST_HPP