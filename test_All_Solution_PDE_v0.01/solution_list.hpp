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

#ifndef VCP_SOLUTION_LIST_HPP
#define VCP_SOLUTION_LIST_HPP

#include <iostream>
#include <list>
#include <cmath>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/pdblas.hpp>
#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

#include <vcp/vcp_timer.hpp>


namespace vcp {
template <typename _T, class _P>
    class solution_list{
    protected:
        std::list< vcp::matrix< _T, _P >  > unknowns;
        std::list< vcp::matrix< _T, _P >  > solutions;
        int nonexistence_count;
        int out_of_range_count;

        int max_length_index( const vcp::matrix< _T, _P >& tmp ){
            int imax = 0;
            auto tmax = rad(tmp(0));
            for (int i = 1; i < tmp.rowsize(); i++ ){
                auto t = rad(tmp(i));
                if ( tmax < t ){
                    tmax = t;
                    imax = i;
                }
            }
            return imax;
        }

    public:
        solution_list(){
            this->nonexistence_count = 0;
            this->out_of_range_count = 0;
        }

        void push_back_unknown( const vcp::matrix< _T, _P >& tmp ){
            this->unknowns.push_back( tmp );
        }

        vcp::matrix< _T, _P > pop_back_unknown( void ){
            vcp::matrix< _T, _P > tmp = this->unknowns.back();
            this->unknowns.pop_back();
            return tmp;
        }

        std::size_t unknown_size(){
            return this->unknowns.size();
        }

        std::size_t solutions_size(){
            return this->solutions.size();
        }

        std::size_t nonexistence_size(){
            return this->nonexistence_count;
        }

        void nonexistence_counter(){
            this->nonexistence_count++;
        }

        std::size_t out_of_range_count_size(){
            return this->out_of_range_count;
        }

        void out_of_range_counter(){
            this->out_of_range_count++;
        }

        void dvided_push_back_unknown( const vcp::matrix< _T, _P >& l ){
            int index = this->max_length_index( l );
            vcp::matrix< _T, _P > tmp = l;
            auto c = mid(l(index));
            if (c == 0){
                c = (l(index).upper() - c)/1.5;
            }
            tmp(index).upper() = c;
            this->push_back_unknown( tmp );
                
            tmp(index) = l(index);
            tmp(index).lower() = c;
            this->push_back_unknown( tmp );
        }

    };
}
#endif