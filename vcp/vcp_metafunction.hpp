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

#ifndef VCP_METAFUNCTION_HPP
#define VCP_METAFUNCTION_HPP


namespace vcp {
	template < typename T > struct is_round_control{
			static constexpr bool value = false;
	};

	#if defined(RDOUBLE_HPP)
	template < > struct  is_round_control< double > {
		static constexpr bool value = true;
	};
	#endif

	#if defined(DD_HPP) && defined(RDD_HPP)
	template < > struct  is_round_control< kv::dd > {
		static constexpr bool value = true;
	};
	#endif

	#if defined(RMPFR_HPP) && defined(MPFR_HPP)
	template < int _N > struct  is_round_control< kv::mpfr< _N > > {
		static constexpr bool value = true;
	};
	#endif



	template < typename T > struct is_interval {
		static constexpr bool value = false;
	};

	#if defined(INTERVAL_HPP)
	template < typename T > struct is_interval< kv::interval< T > > {
		static constexpr bool value = true;
	};
	#endif



	template < typename T, typename = void > struct is_round_interval {
		static constexpr bool value = false;
	};

	#if defined(INTERVAL_HPP)
	template < typename T > struct is_round_interval< kv::interval< T >, kv::interval< typename std::enable_if< vcp::is_round_control< T >::value, T >::type > > {
		static constexpr bool value = true;
	};
	#endif

	template < typename T > struct is_psa {
		static constexpr bool value = false;
	};

	#if defined(PSA_HPP)
	template < typename T > struct is_psa< kv::psa< T > > {
		static constexpr bool value = true;
	};
	#endif

	template < typename T > struct is_hermite_series {
		static constexpr bool value = false;
	};

}
#endif //VCP_METAFUNCTION_HPP