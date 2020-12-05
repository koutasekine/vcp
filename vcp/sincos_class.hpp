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

#ifndef VCP_SINCOS_CLASS_HPP
#define VCP_SINCOS_CLASS_HPP

#include <iostream>

namespace vcp {
    template <typename _T>
    class sincos{
        public:
        _T a; // a*sin
        _T b; // b*cos

        sincos< _T >() {
			(*this).a = _T(0);
			(*this).b = _T(0);
		}
		~sincos< _T >() = default;
		sincos< _T >(const sincos< _T >&) = default;
		sincos< _T >(sincos< _T >&&) = default;
		sincos< _T >& operator=(const sincos< _T >& A) = default;
		sincos< _T >& operator=(sincos< _T >&& A) = default;

        template <typename _Tm> 
        typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type set_a( const _Tm& aa){
            (*this).a = _T(aa);
        }

        template <typename _Tm>
        typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type set_b( const _Tm& bb){
            (*this).b = _T(bb);
        }

        template <typename _Tm>
        typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type set_ab( const _Tm& aa, const _Tm& bb){
            (*this).a = _T(aa);
            (*this).b = _T(bb);
        }
		// tau > 0
		// a*sin(nt + tau) =  a*cos(tau)sin(nt) + a*sin(tau)cos(nt)
		// b*cos(nt + tau) = -b*sin(tau)sin(nt) + b*cos(tau)cos(nt)
		// c1*sin(nt): c1 = a*cos(tau) - b*sin(tau)
		// c2*cos(nt): c2 = a*sin(tau) + b*cos(tau)

		// tau < 0
		// a*sin(nt - abstau) =  a*cos(abstau)sin(nt) - a*sin(abstau)cos(nt)
		// b*cos(nt - abstau) =  b*sin(abstau)sin(nt) + b*cos(abstau)cos(nt)
        template <typename _Tm>
        typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type time_delay( const _Tm& tau){
            using std::sin;
            using std::cos;
            using std::abs;
            _T abstau = abs(_T(tau));
            _T aa = (*this).a;
            _T bb = (*this).b;

            if (tau > 0){
                (*this).a = aa*cos(abstau) - bb*sin(abstau);
                (*this).b = bb*cos(abstau) + aa*sin(abstau);
            }
            else{
                (*this).a = aa*cos(abstau) + bb*sin(abstau);
                (*this).b = bb*cos(abstau) - aa*sin(abstau);
            }
        }

        friend sincos< _T > operator+(const sincos< _T >& A) {
			return A;
		}
        
		friend sincos< _T > operator+(sincos< _T >&& A) {
			return std::move(A);
		}

        friend sincos< _T > operator+(const sincos< _T >& A, const sincos< _T >& B) {
			sincos< _T > C;
			C = A;
			C.a += B.a;
            C.b += B.b;
			return std::move(C);
		}

		friend sincos< _T > operator+(sincos< _T >&& A, const sincos< _T >& B) {
			A.a += B.a;
            A.b += B.b;
			return std::move(A);
		}

		friend sincos< _T > operator+(const sincos< _T >& A, sincos< _T >&& B) {
			B.a += A.a;
            B.b += A.b;
			return std::move(B);
		}

		friend sincos< _T > operator+(sincos< _T >&& A, sincos< _T >&& B) {
			A.a += B.a;
            A.b += B.b;
			return std::move(A);
		}

        friend sincos< _T >& operator+=(sincos< _T >& A, const sincos< _T >& B) {
			A.a += B.a;
            A.b += B.b;
			return A;
		}

        friend sincos< _T > operator-(const sincos< _T >& A) {
            sincos< _T > B;
            B.a = -A.a;
            B.b = -A.b;
			return B;
		}
		friend sincos< _T > operator-(sincos< _T >&& A) {
            A.a = -A.a;
            A.b = -A.b;
			return std::move(A);
		}

        friend sincos< _T > operator-(const sincos< _T >& A, const sincos< _T >& B) {
			sincos< _T > C;
			C = A;
			C.a -= B.a;
            C.b -= B.b;
			return std::move(C);
		}

		friend sincos< _T > operator-(sincos< _T >&& A, const sincos< _T >& B) {
			A.a -= B.a;
            A.b -= B.b;
			return std::move(A);
		}

		friend sincos< _T > operator-(const sincos< _T >& A, sincos< _T >&& B) {
			B.a -= A.a;
            B.b -= A.b;
			return std::move(B);
		}

		friend sincos< _T > operator-(sincos< _T >&& A, sincos< _T >&& B) {
			A.a -= B.a;
            A.b -= B.b;
			return std::move(A);
		}

		friend sincos< _T >& operator-=(sincos< _T >& A, const sincos< _T >& B) {
			A.a -= B.a;
            A.b -= B.b;
			return A;
		}


        template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, sincos< _T > >::type operator*(const sincos< _T >& B, const _Tm& a) {
			_T Ta = _T(a);
			sincos< _T > C;
			C = B;
			C.a *= Ta;
            C.b *= Ta;
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, sincos< _T > >::type operator*(sincos< _T >&& B, const _Tm& a) {
			_T Ta = _T(a);
			B.a *= Ta;
            B.b *= Ta;
			return std::move(B);
		}

        template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, sincos< _T > >::type operator*(const _Tm& a, const sincos< _T >& B) {
			_T Ta = _T(a);
			sincos< _T > C;
			C = B;
			C.a *= Ta;
            C.b *= Ta;
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, sincos< _T > >::type operator*(const _Tm& a, sincos< _T >&& B) {
			_T Ta = _T(a);
			B.a *= Ta;
            B.b *= Ta;
			return std::move(B);
		}

        template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, sincos< _T > >::type operator/(const sincos< _T >& B, const _Tm& a) {
			_T Ta = _T(a);
			sincos< _T > C;
			C = B;
			C.a /= Ta;
            C.b /= Ta;
			return std::move(C);
		}
		template <typename _Tm> friend typename std::enable_if<std::is_constructible< _T, _Tm >::value, sincos< _T > >::type operator/(sincos< _T >&& B, const _Tm& a) {
			_T Ta = _T(a);
			B.a /= Ta;
            B.b /= Ta;
			return std::move(B);
		}

    };


}


#endif // VCP_SINCOS_CLASS_HPP