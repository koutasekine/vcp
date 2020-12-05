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

#ifndef VCP_IMPLICIT_RK_HPP
#define VCP_IMPLICIT_RK_HPP

#include <vector>
#include <tuple>
#include <vcp/irk_legendre_parameter.hpp>
#include <vcp/newton.hpp>

namespace vcp {
    template < template < typename _TT, class _PP > class _FUNC, typename _T, class _P, class _PIRK = IRK_Legendre_Parameter< _T, _P, kv::interval< kv::mpfr< 500 > > > > 
    class Implicit_RK : public _FUNC< _T, _P >, public _PIRK, public Newton< _T, _P > {
        protected:
        vcp::matrix< _T, _P > x_t0;
        vcp::matrix< _T, _P > x_t1;
        std::vector< vcp::matrix< _T, _P > > x_vec;
        vcp::matrix< _T, _P > x;
        
        _T tau;
        _T t0;
        _T t1;
        _T t_start;
        _T t_end;

        int dim;

        bool logger_flag;
        bool tau_refinement_flag;
        int number_local_step;
        std::vector< std::tuple< _T , vcp::matrix< _T, _P > > > ode_logger;

        void x_to_x_vec(){ 
            for ( int i = 0; i < (*this).order; i++) {
                for ( int j = 0; j < (*this).dim; j++ ){
                    ((*this).x_vec[i])(j, 0) = x((*this).dim * i + j, 0);
                }
            }
        }

        void save_data_logger_ode(){
            if ((*this).logger_flag){
                if ((*this).number_local_step == 1){
                    ode_logger.push_back(std::make_tuple((*this).t0, x_t0));
                }
                for (int i = 0; i < (*this).order; i++){
                    ode_logger.push_back(std::make_tuple((*this).t0 + (*this).cpoint(i)*(*this).tau, (*this).x_vec[i]));
                }
                ode_logger.push_back(std::make_tuple((*this).t1, x_t1));
            }
        }

        public:
        Implicit_RK< _FUNC, _T, _P, _PIRK >() {
            (*this).t_start = _T(0);
            (*this).t_end = _T(1000000);
            (*this).t0 = _T(0);
            (*this).tau = 1/_T(8);
            (*this).t1 = (*this).t0 + (*this).tau;
            (*this).logger_flag = true;
            (*this).tau_refinement_flag = true;
            (*this).number_local_step = 0;
		}

        void setting(){
            (*this).setting_parameter();
            (*this).x_vec.resize((*this).order);
        }

        void setting( int n ){
            (*this).setting_parameter(n);
            (*this).x_vec.resize((*this).order);
        }

        void setting_local_time( _T hh ){
            (*this).tau = hh;
            (*this).t0 = (*this).t_start;
            (*this).t1 = (*this).t0 + (*this).tau;
        }

        void setting_local_time( _T hh, _T t ){
            (*this).tau = hh;
            (*this).t0 = t;
            (*this).t1 = (*this).t0 + (*this).tau;
        }

        void change_local_time( _T hh ){
            (*this).tau = hh;
            (*this).t1 = (*this).t0 + (*this).tau;
        }
        
        void setting_global_time(_T ts, _T te){
            (*this).t_start = ts;
            (*this).t_end = te;
        }

        void setting_order( int n = 4 ){
            (*this).setting_parameter(n);
            (*this).x_vec.resize((*this).order);
        }

        void setting_logger(bool lflag){
            (*this).logger_flag = lflag;
        }

        void setting_auto_refinement(bool lflag){
            (*this).tau_refinement_flag = lflag;
        }

        void initial_value(const vcp::matrix< _T, _P >& xh){
            (*this).x_t0 = xh;
            (*this).dim = xh.rowsize();
            (*this).x.zeros((*this).order * (*this).dim, 1);
            for ( int i = 0; i < (*this).order; i++) {
                 (*this).x_vec[i] = (*this).x_t0;
                 for (int j = 0; j < (*this).dim; j++ ){
                    x(i*(*this).dim + j, 0) = (*this).x_t0(j, 0);
                 }
            }
        }

        void setting_newton( vcp::matrix< _T, _P >& xh ) override {
            (*this).x = xh;
        }

        vcp::matrix< _T, _P > f() override {
            (*this).x_to_x_vec();
            std::vector< vcp::matrix< _T, _P > > fx_vec;
            vcp::matrix< _T, _P > xtmp;
            xtmp.zeros((*this).order * (*this).dim, 1);

            vcp::matrix< _T, _P > tmp;
            fx_vec.resize((*this).order);

            for (int i = 0; i < (*this).order; i++){
                fx_vec[i] = (*this).x_vec[i] - (*this).x_t0;
            }
            for (int j = 0; j < (*this).order; j++){
                tmp = (*this).tau * (*this).func( t0 + (*this).cpoint(j) *(*this).tau, (*this).x_vec[j] );
                for (int i = 0; i < (*this).order; i++ ){
                    fx_vec[i] -= (*this).alpha(i, j) * tmp;
                }
            }
            
            for ( int i = 0; i < (*this).order; i++) {
                for ( int j = 0; j < (*this).dim; j++ ){
                    xtmp((*this).dim * i + j, 0) = (fx_vec[i])(j, 0);
                }
            }
            return xtmp;
        }

        vcp::matrix< _T, _P > Df() override {
            (*this).x_to_x_vec();
            vcp::matrix< _T, _P > tmp1;
            vcp::matrix< _T, _P > tmp2;
            vcp::matrix< _T, _P > dfx;
            dfx.eye( (*this).order * (*this).dim);

            for (int j = 0; j < (*this).order; j++){
                tmp1 = (*this).tau * (*this).diffunc(t0 + (*this).cpoint(j) *(*this).tau, (*this).x_vec[j] );
                for (int i = 0; i < (*this).order; i++ ){
                    tmp2 = (*this).alpha(i, j) * tmp1;
                    for (int ii = 0; ii < (*this).dim; ii++ ){
                        for (int jj = 0; jj < (*this).dim; jj++){
                            dfx(i * (*this).dim + ii, j * (*this).dim + jj) -= tmp2(ii, jj);
                        }
                    }
                }
            }
            return dfx;
        }

        void disp_continue() override {
        //    std::cout << "i = " << iteration_newton << ", " << Correction_term << " <= " << (*this).newton_tol << std::endl;
        }

        bool local_solve_ode(){
            while (true){
                (*this).x = (*this).solve_nls((*this).x);
                
                if ( (*this).convergence_condition() ) break;
                if ((*this).tau_refinement_flag){
                    (*this).not_convergence();
                    (*this).initial_value((*this).x_t0);
                }
                else {
                    return false;
                }
            }
            (*this).x_t1 = (*this).x_t0;
            (*this).x_to_x_vec();
            for (int i = 0; i < (*this).order; i++){
                (*this).x_t1 += (*this).tau * (*this).bweight(i) * (*this).func( t0 + (*this).cpoint(i) *(*this).tau, (*this).x_vec[i] );
            }

            ++(*this).number_local_step;
            (*this).save_data_logger_ode();
            return true;
        }

        std::tuple< std::vector< _T >, std::vector< vcp::matrix< _T, _P > > > local_output(){
            std::vector< _T > t_local;
            t_local.resize((*this).order + 2);
            t_local[0] = (*this).t0;
            for (int i = 0; i < (*this).order; i++ ){
                t_local[i+1] = (*this).t0 + (*this).cpoint(i)*(*this).tau;
            }
            t_local[(*this).order + 1] = (*this).t1;
            
            std::vector< vcp::matrix< _T, _P > > x_local_output;
            x_local_output.resize((*this).order + 2);
            x_local_output[0] = (*this).x_t0;
            for (int i = 0; i < (*this).order; i++ ){
                x_local_output[i+1] = (*this).x_vec[i];
            }
            x_local_output[(*this).order + 1] = (*this).x_t1;

            return std::make_tuple(t_local, x_local_output);
        }

        void next_step(){
            (*this).setting_local_time( (*this).tau, (*this).t1 );
            (*this).initial_value((*this).x_t1);
        }

        std::vector< std::tuple< _T , vcp::matrix< _T, _P > > > output(){
            if (!(*this).logger_flag){
                std::cout << "Please logger flag = true (setting_logger(true))" << std::endl;
            }
            return (*this).ode_logger;
        }

        void solve_ode(){
            int step_num = 0; 
            while(true){
                (*this).local_solve_ode();
                ++step_num;
                std::cout << "Step : " << step_num << ", t in [" << (*this).t0 << ", " << (*this).t1 << "]"  << std::endl;
                if ((*this).t_end < (*this).t1) break;
                (*this).next_step();
            }
        }
        /*

        virtual vcp::matrix< _T, _P > func( vcp::matrix< _T, _P >& xh ){
            return xh;
        }
        virtual vcp::matrix< _T, _P > diffunc( vcp::matrix< _T, _P >& xh ){
            return xh;
        }
        */

        virtual void not_convergence(){ 
            (*this). setting_local_time( (*this).tau / _T(2), (*this).t0 );
            std::cout << "Refining tau : " << (*this).tau << std::endl;
        }

        virtual bool convergence_condition(){
            return (*this).flag_Convergence;
        }
    };
}
#endif // VCP_IMPLICIT_RK_HPP