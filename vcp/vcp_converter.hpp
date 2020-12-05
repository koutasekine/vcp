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
//
//
//
// The vcp_convert.hpp is based on interval_conv.hpp in kv library :
// Copyright(c) 2015 Masahide Kashiwagi(kashi@waseda.jp)
// Released under the MIT Licenses.
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#ifndef VCP_CONVERTER_HPP
#define VCP_CONVERTER_HPP

namespace vcp {
#if defined(RDOUBLE_HPP) && defined(RMPFR_HPP) && defined(MPFR_HPP)
	template <int N> void convert(const kv::mpfr<N>& x, double& y, int rnd = 0)
	{
		mp_rnd_t mode;

		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		y = mpfr_get_d(x.a, mode);
	}
	template <int N> void convert(const double& x, kv::mpfr<N>& y) {
		y = x;
	}
#endif

#if defined(DD_HPP) && defined(RDD_HPP) && defined(RMPFR_HPP) && defined(MPFR_HPP)
	template <int N> void convert(const kv::dd& x, kv::mpfr<N>& y, int rnd = 0)
	{
		mp_rnd_t mode;
		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		// y = x.a;
		mpfr_set_d(y.a, x.a1, mode);
		mpfr_add_d(y.a, y.a, x.a2, mode);
	}
	template <int N> void convert(const kv::mpfr<N>& x, kv::dd& y, int rnd = 0)
	{
		kv::mpfr<N> mtmp;
		double dtmp1, dtmp2;

		mp_rnd_t mode;

		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		convert(x, dtmp1, 0);

		// mtmp = x - dtmp1;
		// theoretically no rounding error.
		// use rounded subtraction just to be safe.
		mpfr_sub_d(mtmp.a, x.a, dtmp1, mode);

		dtmp2 = mpfr_get_d(mtmp.a, mode);

		kv::dd::twosum(dtmp1, dtmp2, y.a1, y.a2);
	}
#endif

#if defined(MPFR_HPP) && defined(RMPFR_HPP)
	template <int N, int M> void convert(const kv::mpfr<N>& x, kv::mpfr<M>& y, int rnd = 0)
	{
		mp_rnd_t mode;

		if (rnd == 1) {
			mode = MPFR_RNDU;
		}
		else if (rnd == -1) {
			mode = MPFR_RNDD;
		}
		else {
			mode = MPFR_RNDN;
		}

		mpfr_set(y.a, x.a, mode);
	}
#endif

#if defined(DD_HPP) && defined(RDD_HPP) && defined(RDOUBLE_HPP)
	void convert(const kv::dd& x, double& y, int rnd = 0)
	{
		if (rnd == 1) {
			kv::rop<double>::begin();
			y = kv::rop<double>::add_up(x.a1, x.a2);
			kv::rop<double>::end();
		}
		else if (rnd == -1) {
			kv::rop<double>::begin();
			y = kv::rop<double>::add_down(x.a1, x.a2);
			kv::rop<double>::end();
		}
		else {
			y = x.a1 + x.a2;
		}
	}
	void convert(const double& x, kv::dd& y) {
		y = x;
	}
#endif

#if defined(INTERVAL_HPP) && defined(DD_HPP) && defined(RDD_HPP) && defined(RDOUBLE_HPP)
	void convert(const kv::interval< kv::dd >& x, kv::interval<double>& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
	void convert(const kv::interval<double>& x, kv::interval< kv::dd >& y)
	{
		y.lower() = x.lower();
		y.upper() = x.upper();
	}
#endif

#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP) && defined(RDOUBLE_HPP)
	template <int N> void convert(const kv::interval< kv::mpfr<N> >& x, kv::interval<double>& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
	template <int N> void convert(const kv::interval<double>& x, kv::interval< kv::mpfr<N> >& y)
	{
		y.lower() = x.lower();
		y.upper() = x.upper();
	}
#endif

#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP) && defined(DD_HPP) && defined(RDD_HPP)
	template <int N> void convert(const kv::interval< kv::mpfr<N> >& x, kv::interval< kv::dd >& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
	template <int N> void convert(const kv::interval< kv::dd >& x, kv::interval< kv::mpfr<N> >& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
#endif

#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP)
	template <int N, int M> void convert(const kv::interval< kv::mpfr<N> >& x, kv::interval< kv::mpfr<M> >& y)
	{
		convert(x.lower(), y.lower(), -1);
		convert(x.upper(), y.upper(), 1);
	}
#endif

	void convert(const int& x, double& y) {
		y = x;
	}

#if defined(INTERVAL_HPP) && defined(RDOUBLE_HPP)
	template<typename _T> void convert(const kv::interval< _T >& x, double& y) {
		convert(mid(x), y);
	}
	template<typename _T> void convert(const double& x, kv::interval< _T >& y) {
		kv::interval< double > yy = kv::interval< double >(x);
		convert(yy, y);
	}
	void convert(const int& x, kv::interval< double >& y) {
		y = x;
	}
#endif
#if defined(INTERVAL_HPP) && defined(DD_HPP) && defined(RDD_HPP)
	template<typename _T> void convert(const kv::interval< _T >& x, kv::dd& y) {
		convert(mid(x), y);
	}
	template<typename _T> void convert(const kv::dd& x, kv::interval< _T >& y) {
		kv::interval< kv::dd > yy = kv::interval< kv::dd >(x);
		convert(yy, y);
	}
	void convert(const int& x, kv::dd& y) {
		y = x;
	}
	void convert(const int& x, kv::interval< kv::dd >& y) {
		y = x;
	}
#endif
#if defined(INTERVAL_HPP) && defined(MPFR_HPP) && defined(RMPFR_HPP)
	template<typename _T, int N> void convert(const kv::interval< _T >& x, kv::mpfr< N >& y) {
		convert(mid(x), y);
	}
	template<typename _T, int N> void convert(const kv::mpfr< N >& x, kv::interval< _T >& y) {
		kv::interval< kv::mpfr< N > > yy = kv::interval< kv::mpfr< N > >(x);
		convert(yy, y);
	}
	template<int N> void convert(const int& x, kv::mpfr< N >& y) {
		y = x;
	}
	template<int N> void convert(const int& x, kv::interval< kv::mpfr< N > >& y) {
		y = x;
	}
#endif

	template<typename _T> void convert(const _T& x, _T& y) {
		y = x;
	}
}
#endif //VCP_CONVERTER_HPP