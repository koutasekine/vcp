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

#ifndef VCP_INIT_VECTOR_HPP
#define VCP_INIT_VECTOR_HPP

#include <initializer_list>

namespace vcp{
 
    template <typename _T, class _P, typename _Tm> 
    typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type init( matrix< _T, _P >& xh, const std::initializer_list< _Tm >& list1){
        if ( xh.matstype() == 'M' ){
            std::cout << "ERROR : init_vector :  xh is Matrix..." << std::endl;
            exit(1);
        } 
        if ( list1.size() > xh.rowsize() && xh.matstype() == 'C' ) {
            xh.clear();
        }
        if ( list1.size() > xh.columnsize() && xh.matstype() == 'R' ) {
            xh.clear();
        }
        if (xh.matstype() == 'N'){
            xh.zeros(list1.size(), 1);
        }
        std::vector<_Tm> l1 = list1;

        for (int i = 0; i < list1.size(); i++){
            xh(i) = _T(l1[i]);
        }
    }

    template <typename _T, class _P, typename _Tm> 
    typename std::enable_if< std::is_constructible< _T, _Tm >::value, void >::type init_vector( matrix< _T, _P >& xh, const std::initializer_list< _Tm >& list1){
        if ( xh.matstype() == 'M' ){
            std::cout << "ERROR : init_vector :  xh is Matrix..." << std::endl;
            exit(1);
        } 
        if ( list1.size() > xh.rowsize() && xh.matstype() == 'C' ) {
            xh.clear();
        }
        if ( list1.size() > xh.columnsize() && xh.matstype() == 'R' ) {
            xh.clear();
        }
        if (xh.matstype() == 'N'){
            xh.zeros(list1.size(), 1);
        }
        std::vector<_Tm> l1 = list1;

        xh(0) = 2*_T(l1[0]);
        for (int i = 1; i < list1.size(); i++){
            xh(i) = _T(l1[i]);
        }
    }
}
#endif