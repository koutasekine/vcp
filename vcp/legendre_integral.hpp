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

#ifndef VCP_LEGENDRE_INTEGRAL_HPP
#define VCP_LEGENDRE_INTEGRAL_HPP

#ifndef INTERVAL_HPP
#error Please include interval.hpp
#endif

#include <stack>
#include <vector>

#include <kv/psa.hpp>
#include <kv/gamma.hpp>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/vcp_psa_assist.hpp>
#include <vcp/vcp_converter.hpp>
#include <vcp/matrix.hpp>

#include<omp.h>

namespace vcp {
    namespace GaussLegendreIntegral {
		template <typename _T> void LegendrePol(const int nn, kv::psa< _T >& outp, kv::psa< _T >& outpm1) {
			kv::psa< _T > xp, oim1, oim2, ZZZ;
			int i, n;
			_T iT;
			n = nn + 1;
			xp.v.resize(n);
			oim1.v.resize(n);
			oim2.v.resize(n);
			outp.v.resize(n);
			outpm1.v.resize(n);
			ZZZ.v.resize(n);
			for (i = 0; i < n; i++) {
				xp.v(0) = _T(0);
				oim1.v(0) = _T(0);
				oim2.v(0) = _T(0);
				outp.v(0) = _T(0);
				outpm1.v(0) = _T(0);
				ZZZ.v(0) = _T(0);
			}
			xp.v(1) = _T(1);
			oim1.v(1) = _T(1);
			oim2.v(0) = _T(1);
			outp = 1 / _T(2) * (3 * xp*oim1 - oim2);

			for (i = 3; i < n; i++) {
				iT = _T(i);
				oim2 = oim1;
				oim1 = outp;
				outp = 1 / iT * ((2 * (iT - 1) + 1)*xp*oim1 - (iT - 1)*oim2);
				if (i == n - 2) {
					outpm1 = outp;
				}
			}
		}
		template <typename _T> void LegendreDifFunc(const kv::psa< _T >& pol, const _T x, _T& y) {
			int i;
			int n = pol.v.size();
			y = (n - 1)*pol.v(n - 1);
			for (i = n - 2; i>0; i--) {
				y = y*x + i*pol.v(i);
			}
		}
		template <typename _T> void psa_intvalTtoT(const kv::psa< kv::interval< _T > >& pol, kv::psa< _T >& polout) {
			int i;
			int n = pol.v.size();
			polout.v.resize(n);
			for (i = 0; i<n; i++) {
				polout.v(i) = mid(pol.v(i));
			}
		}
		template <typename _T> void AppNewton(const kv::psa< _T >& pol, const _T x, _T& y) {
			_T res, dif;
			y = x;
			for (int i = 0; i < 10; i++) {
				psa_value_Horner(pol, y, res);
				LegendreDifFunc(pol, y, dif);
				y = y - res / dif;
			}
		}
		template <typename _T> void Krawczyk1d(const kv::psa< kv::interval< _T > >& pol, const _T x, kv::interval< _T >& y) {
			kv::interval< _T > res, Ky, dif;
			_T difmid, alpha;

			y = x;
			psa_value_Horner(pol, y, res);
			LegendreDifFunc(pol, y, dif);
			difmid = mid(dif);
			alpha = mag(res / difmid);
			y = y + kv::interval< _T >(-2 * alpha, 2 * alpha);
			int i = 0;
			while (1) {
				psa_value_Horner(pol, kv::interval< _T >(mid(y)), res);
				LegendreDifFunc(pol, y, dif);
				difmid = mid(dif);
				Ky = mid(y) - res / difmid + (1 - dif / difmid)*(y - mid(y));
				if (i == 1) {
					Ky = intersect(Ky, y);
					if (rad(Ky) > 0.9*rad(y)) break;
				}
				else if (proper_subset(Ky, y)) {
					i = 1;
				}
				else {
					std::cout << "Ahhhhhhhhhhh" << std::endl;
					std::cout << "y= " << y << std::endl;
					std::cout << "Ky= " << Ky << std::endl;
				}
				y = Ky;
			}
		}
		template <typename _T> void LPointWeight(const int n, std::vector< kv::interval< _T > >& Point, std::vector< kv::interval< _T > >& Weight) {
			_T x, y, K1, K2;
			_T PI = kv::constants< _T >::pi();
			kv::interval< _T > res, dif;
			kv::psa< _T > psapp;
			kv::psa< kv::interval< _T > > ps, psm1;

			Point.resize(n);
			Weight.resize(n);

			LegendrePol(n, ps, psm1);
			psa_intvalTtoT(ps, psapp);

			K1 = (n - 1) / (8 * n*n*_T(n));
			K2 = _T(4 * n + 2);

			//std::cout.precision(32);
			for (int i = 1; i <= n; i++) {
				x = (1 - K1)*cos((4 * i - 1) / K2*PI);
				AppNewton(psapp, x, y);
				Krawczyk1d(ps, y, Point[i - 1]);
			}

			for (int i = 0; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {
					if (overlap(Point[i], Point[j])) {
						std::cout << "Ahhhhhhhhhhhhhhh" << std::endl;
					}
				}
			}

			for (int i = 1; i <= n; i++) {
				psa_value_Horner(psm1, Point[i - 1], res);
				LegendreDifFunc(ps, Point[i - 1], dif);
				Weight[i - 1] = 2 / (n * res * dif);
			}
		}
	}

	template <typename _T> class interval_ld_weightpoint {
	public:
		std::vector< _T > Weight;
		std::vector< _T > Point;
		// Order of Legendre integral
		int n;

		interval_ld_weightpoint< _T >() {
		}

		void set(const int nn) {
			n = nn;
			Weight.resize(n);
			Point.reserve(n);
			GaussLegendreIntegral::LPointWeight(n, Point, Weight);
		}

	};
}
#endif // VCP_LEGENDRE_INTEGRAL_HPP

