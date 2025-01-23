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



#include <iostream>

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

#include <vcp/fourier_series.hpp>
#include <vcp/fourier_cos_basis.hpp>
#include <vcp/newton.hpp>

#include <vcp/vcp_timer.hpp>

#include <cmath>

/*  -Δu = -uv^2 + c(1-u)  */
/*  -Δv =  uv^2 - (c+k)v  */
template < typename _T, typename _P >
struct Gray_Scott : public vcp::Newton< _T, _P > {
    int m;
    _T c, k, iDu, iDv;
    vcp::fourier_series< _T > u, v;
    vcp::fourier_cos_basis< _T, _P > Generator;

    vcp::matrix< _T, _P > z;
    vcp::matrix< _T, _P > D, L;

    void setting_newton( vcp::matrix< _T, _P >& zh ) override {
        this->z = zh;
        this->u = Generator.vec_to_fourier_series( zh.submatrix( {0, zh.rowsize()/2 - 1}, {0}) );
        this->v = Generator.vec_to_fourier_series( zh.submatrix( {zh.rowsize()/2, zh.rowsize() - 1}, {0}) );
    }

    vcp::fourier_series< _T > func1(void){
        return this->iDu*(-this->u*this->v*this->v + this->c*(_T(1) - this->u));
    }

    vcp::fourier_series< _T > func2(void){
        return  this->iDv*(this->u*this->v*this->v - (this->c+this->k)*this->v);
    }

    vcp::matrix< _T, _P > f() override {			
        vcp::matrix< _T, _P > fzh = vercat( Generator.fphi( this->func1() ), Generator.fphi( this->func2() ) );
        return this->D*this->z - fzh;
    }

    /*  -v^2 - c, -2uv  */
    vcp::matrix< _T, _P > Q1(void) {
            return horzcat( Generator.Dfphiphi( this->iDu*(- this->v*this->v  - this->c)), Generator.Dfphiphi( this->iDu*(-_T(2)*this->u*this->v) ) );
    }

    /*  v^2 - (c+k), 2uv - (c+k)  */
    vcp::matrix< _T, _P > Q2(void) {
            return horzcat( Generator.Dfphiphi( this->iDv*(this->v*this->v) ), Generator.Dfphiphi( this->iDv*(_T(2)*this->u*this->v - (this->c + this->k) )) );
    }

    vcp::matrix< _T, _P > Df() override {
        return this->D - vercat( this->Q1(), this->Q2() );
    }


    void setting_m(const int& mm){
        this->m = mm;
        this->Generator.setting_m( this->m );

        vcp::matrix< _T, _P > tmp = this->Generator.dphidphi();
        
        vcp::matrix< _T, _P > zero;
        zero.zeros(tmp.rowsize(), tmp.columnsize());
        this->D = vercat( horzcat(tmp, zero), horzcat(zero, tmp) );
        
        tmp = this->Generator.phiphi();
        this->L = vercat( horzcat(tmp, zero), horzcat(zero, tmp) );
    }

    void setting_parameter(const _T& cc, const _T& kk, const _T& Du, const _T& Dv){
        this->c = cc;
        this->k = kk;
        this->iDu = _T(1)/Du;
        this->iDv = _T(1)/Dv;
    }
};

typedef double AppData;
typedef vcp::pdblas AppPolicy;

typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

int main(void){
	std::cout.precision(17);

	std::string c = "0.06";
	std::string k = "0.065";
    std::string Du = "0.02";
    std::string Dv = "0.01";

    int m_app = 100;

    vcp::matrix< VData, VPOLICY > izh;
    {
        vcp::matrix< AppData, AppPolicy > uh, vh, zh;
        Gray_Scott< AppData, AppPolicy > gray_scott;
        gray_scott.setting_parameter( std::stod(c), std::stod(k), std::stod(Du), std::stod(Dv) );
        gray_scott.setting_m( m_app );
        gray_scott.setting_newton_tol(32);
        
        do {
            uh.rand(m_app+1, 1);
            vh.rand(m_app+1, 1);
            uh(1,0) = 100;
            vh(1,0) = 100;

            zh = vercat( uh, vh );
            zh = gray_scott.solve_nls( zh );
        } while(!gray_scott.is_convergence());
    //    std::cout << zh << std::endl;
            
        auto graf = gray_scott.Generator.graphics(gray_scott.u, 500);
        std::cout << "u = " << std::endl;
        std::cout << graf << std::endl;

        graf = gray_scott.Generator.graphics(gray_scott.v, 500);
        std::cout << "v = " << std::endl;
        std::cout << graf << std::endl;

        vcp::interval( zh, izh );
    }

    {
        Gray_Scott< VData, VPOLICY > gray_scott;
        gray_scott.setting_parameter( VData(c), VData(k), VData(Du), VData(Dv) );
        gray_scott.setting_m( m_app );
        gray_scott.setting_newton( izh );

        vcp::matrix< VData, VPOLICY > G = gray_scott.Df();
        vcp::matrix< VData, VPOLICY > L = gray_scott.L;
        vcp::matrix< VData, VPOLICY > E;

        
        
    }


}