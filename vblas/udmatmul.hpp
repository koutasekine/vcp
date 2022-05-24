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

#ifndef VBLAS_UDMATMUL_HPP
#define VBLAS_UDMATMUL_HPP

#define NB_L1 16
#define NB_L2 32
#define NB_L3 64

#include <immintrin.h>

// A: m*n matrix
// B: n*k matrix

void udmatmul(int m, int n, int k, double *A, double *B, double *CU, double *CD){

    int NL1, NL2, NL0;
    int m_L3, n_L3, k_L3;

    m_L3 = m / NB_L3 * NB_L3;
    n_L3 = n / NB_L3 * NB_L3;
    k_L3 = k / NB_L3 * NB_L3;    

    NL2 = NB_L3 / NB_L2 * NB_L2;
    NL1 = NB_L2 / NB_L1 * NB_L1;
    NL0 = m/8 * 8;

    #pragma omp parallel for
    for (int iL3 = 0; iL3 < m_L3; iL3 += NB_L3){
    for (int kL3 = 0; kL3 < k_L3; kL3 += NB_L3){
    for (int jL3 = 0; jL3 < n_L3; jL3 += NB_L3){
        for (int jL2 = jL3; jL2 < jL3 + NB_L3; jL2+=NB_L2){
        for (int kL2 = kL3; kL2 < kL3 + NB_L3; kL2+=NB_L2){
        for (int iL2 = iL3; iL2 < iL3 + NB_L3; iL2+=NB_L2){
            for (int jL1 = jL2; jL1 < jL2 + NB_L2; jL1+=NB_L1){
            for (int kL1 = kL2; kL1 < kL2 + NB_L2; kL1+=NB_L1){
            for (int iL1 = iL2; iL1 < iL2 + NB_L2; iL1+=NB_L1){
                
                alignas(64) double AA[NB_L1*NB_L1];
                alignas(64) double CC[NB_L1*NB_L1];
                for (int ii = iL1; ii < iL1 + NB_L1; ii++){
                for (int jj = jL1; jj < jL1 + NB_L1; jj++){
                    AA[(ii - iL1) + NB_L1*(jj - jL1)] = A[ii + m*jj];
                }
                }
                

                for (int kk = kL1; kk < kL1 + NB_L1; kk++){
                for (int jj = jL1; jj < jL1 + NB_L1; jj++){
                    __m512d b = _mm512_set1_pd( B[jj + n*kk] );
                for (int ii = iL1; ii < iL1 + NB_L1; ii+=8){
                    /*
                        _MM_FROUND_TO_NEAREST_INT - rounds to nearest even
                        _MM_FROUND_TO_NEG_INF - rounds to negative infinity
                        _MM_FROUND_TO_POS_INF - rounds to positive infinity
                        _MM_FROUND_TO_ZERO - rounds to zero
                        _MM_FROUND_CUR_DIRECTION - rounds using default from MXCSR register
//                      __m512d a = _mm512_setr_pd( A[ii + m*jj], A[ii + 1 + m*jj], A[ii + 2 + m*jj], A[ii + 3 + m*jj], A[ii + 4 + m*jj], A[ii + 5 + m*jj], A[ii + 6 + m*jj], A[ii + 7 + m*jj] );
//                      c = _mm512_fmadd_pd(a, b, c);
                    */

                    __m512d a = _mm512_load_pd( AA + (ii - iL1) + NB_L1*(jj - jL1) );
                    __m512d c = _mm512_loadu_pd( CU + ii + m*kk );
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_storeu_pd(CU + ii + m*kk, c);

                    c = _mm512_loadu_pd( CD + ii + m*kk );
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_storeu_pd(CD + ii + m*kk, c);
/*
                    a = _mm512_loadu_pd( A + ii + 8 + m*jj );
                    c = _mm512_loadu_pd( C + ii + 8 + m*kk );
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_storeu_pd(C + ii + 8 + m*kk, c);
*/
                }}}
            }}}
        }}}
    }
        for (int jj = n_L3; jj < n; jj++){
            for (int kL2 = kL3; kL2 < kL3 + NB_L3; kL2+=NB_L2){
            for (int iL2 = iL3; iL2 < iL3 + NB_L3; iL2+=NB_L2){
                for (int kL1 = kL2; kL1 < kL2 + NB_L2; kL1+=NB_L1){
                for (int iL1 = iL2; iL1 < iL2 + NB_L2; iL1+=NB_L1){ 
                    for (int kk = kL1; kk < kL1 + NB_L1; kk++){
                        __m512d b = _mm512_set1_pd( B[jj + n*kk] );
                    for (int ii = iL1; ii < iL1 + NB_L1; ii+=8){                   
                        __m512d a = _mm512_loadu_pd( A + ii + m*jj );
                        __m512d c = _mm512_loadu_pd( CU + ii + m*kk );
                        c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                        _mm512_storeu_pd(CU + ii + m*kk, c);

                        c = _mm512_loadu_pd( CD + ii + m*kk );
                        c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                        _mm512_storeu_pd(CD + ii + m*kk, c);
                    }}
                }}
        }}}
    }
        for (int kk = k_L3;    kk < k; kk++){
        for (int jj = 0;    jj < n; jj++){
            __m512d b = _mm512_set1_pd( B[jj + n*kk] );
        for (int iL2 = iL3; iL2 < iL3 + NB_L3; iL2+=NB_L2){
            for (int iL1 = iL2; iL1 < iL2 + NB_L2; iL1+=NB_L1){
                for (int ii = iL1; ii < iL1 + NB_L1; ii+=8){
                    __m512d a = _mm512_loadu_pd( A + ii + m*jj );
                    __m512d c = _mm512_loadu_pd( CU + ii + m*kk );
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_storeu_pd(CU + ii + m*kk, c);

                    c = _mm512_loadu_pd( CD + ii + m*kk );
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_storeu_pd(CD + ii + m*kk, c);
                }
            }
        }}}
    }

    int tmp = 0;
    for (int i = 0; i < m - NL0; i++) tmp += pow(2,i);
    __mmask8 mmask = tmp;
    
    #pragma omp parallel for
    for (int kk = 0;    kk < k; kk++){
    for (int jj = 0;    jj < n; jj++){ 
        __m512d a, c;
        __m512d b = _mm512_set1_pd( B[jj + n*kk] );       
        for (int ii = m_L3; ii < NL0; ii+=8){
            a = _mm512_loadu_pd( A + ii + m*jj );
            c = _mm512_loadu_pd( CU + ii + m*kk );
            c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
            _mm512_storeu_pd(CU + ii + m*kk, c);

            c = _mm512_loadu_pd( CD + ii + m*kk );
            c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            _mm512_storeu_pd(CD + ii + m*kk, c);
        }
        
        a = _mm512_maskz_loadu_pd( mmask, A + NL0 + m*jj );
        c = _mm512_maskz_loadu_pd( mmask, CU + NL0 + m*kk );
        c = _mm512_maskz_fmadd_round_pd( mmask, a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
        _mm512_mask_storeu_pd(CU + NL0 + m*kk, mmask, c);

        c = _mm512_maskz_loadu_pd( mmask, CD + NL0 + m*kk );
        c = _mm512_maskz_fmadd_round_pd( mmask, a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
        _mm512_mask_storeu_pd(CD + NL0 + m*kk, mmask, c);
    }}
}

#endif // VBLAS_UDMATMUL_HPP