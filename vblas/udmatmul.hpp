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
// Created in collaboration with Mr. Obata, a 2022 graduating student in Sekine Laboratory
void udmatmul(int m, int n, int k, double *A, double *B, double *CU, double *CD){
    int mres, nres,kres,NL0;
    mres = m / NB_L2 * NB_L2;
    nres = n / NB_L2 * NB_L2;
    kres = k / NB_L2 * NB_L2;
    NL0 = m/8 * 8;
    int m8just=((m-mres)/8+1)*8;
    int tmp =0;
    for(int i = 0; i < m - NL0; i++)tmp += pow(2,i);
    __mmask8 mmask = tmp;

  #pragma omp parallel for
    for (int iii = 0; iii < mres; iii+=NB_L2){
    for (int jjj = 0; jjj < nres; jjj+=NB_L2){
    for (int zzz = 0; zzz < kres; zzz+=NB_L2){
        for (int ii = iii; ii < iii + NB_L2; ii+=NB_L1){
        for (int jj = jjj; jj < jjj + NB_L2; jj+=NB_L1){
        for (int zz = zzz; zz < zzz + NB_L2; zz+=NB_L1){
            alignas(64) double L1A[NB_L1*NB_L1];
            alignas(64) double L1CU[NB_L1*NB_L1];
            alignas(64) double L1CD[NB_L1*NB_L1];
            for (int j = jj; j < jj + NB_L1; j++){
                for(int i = ii; i < ii + NB_L1; i++){
                    L1A[i-ii + NB_L1*(j-jj)] = A[i + m*j];
                }}
            for(int z = zz; z < zz + NB_L1; z++){
                for(int i = ii; i < ii + NB_L1; i++){
                    L1CU[i-ii + NB_L1*(z-zz)] = CU[i + m*z];
                    L1CD[i-ii + NB_L1*(z-zz)] = CD[i + m*z];
                }}
            for(int z = 0; z <  NB_L1; z++){
            for(int j = 0; j <  NB_L1; j++){
                    __m512d b = _mm512_set1_pd(B[j + jj + n*(z + zz)]);
            for(int i = 0; i <  NB_L1; i+=8){
                    __m512d a = _mm512_load_pd(L1A + i + NB_L1*j);
                    __m512d c = _mm512_load_pd(L1CU + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC); 
                    _mm512_store_pd(L1CU + i + NB_L1*z, c);
                    c = _mm512_load_pd(L1CD + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + NB_L1*z, c);
            }}}
            for(int z = zz; z < zz + NB_L1; z++){
                for(int i = ii; i < ii + NB_L1; i++){
                    CU[i + m*z] = L1CU[i-ii + NB_L1*(z-zz)];
                    CD[i + m*z] = L1CD[i-ii + NB_L1*(z-zz)];
                }}
        }}}
    }
    for (int zres = kres; zres < k; zres+=(k-kres)){
        for (int ii = iii; ii < iii + NB_L2; ii+=NB_L1){
        for (int jj = jjj; jj < jjj + NB_L2; jj+=NB_L1){
            alignas(64) double L1A[NB_L1*NB_L1];
            alignas(64) double L1CU[NB_L1*(k-kres)];
            alignas(64) double L1CD[NB_L1*(k-kres)];
            for(int j = jj; j < jj + NB_L1; j++){
            for(int i = ii; i < ii + NB_L1; i++){
                    L1A[i-ii + NB_L1*(j-jj)] = A[i + m*j];
                }}
            for(int z = zres; z < k; z++){
            for(int i = ii; i < ii + NB_L1; i++){
                    L1CU[i-ii + NB_L1*(z-zres)] = CU[i + m*z];
                    L1CD[i-ii + NB_L1*(z-zres)] = CD[i + m*z];
                }}
            for(int z = 0; z <  k-kres; z++){
            for(int j = 0; j <  NB_L1; j++){
                    __m512d b = _mm512_set1_pd(B[j + jj + n*(z + zres)]);
            for(int i = 0; i <  NB_L1; i+=8){
                    __m512d a = _mm512_load_pd(L1A + i + NB_L1*j);
                    __m512d c = _mm512_load_pd(L1CU + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + NB_L1*z, c);
                    c = _mm512_load_pd(L1CD + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + NB_L1*z, c);
            }}}
            for(int z = zres; z < k; z++){
                for(int i = ii; i < ii + NB_L1; i++){
                    CU[i + m*z] = L1CU[i-ii + NB_L1*(z-zres)];
                    CD[i + m*z] = L1CD[i-ii + NB_L1*(z-zres)];
                }}
        }}
    }}
    for (int jres = nres; jres < n; jres+=(n-nres)){
    for (int zzz = 0; zzz < kres; zzz+=NB_L2){
        for (int ii = iii; ii < iii + NB_L2; ii+=NB_L1){
        for (int zz = zzz; zz < zzz + NB_L2; zz+=NB_L1){

            alignas(64) double L1A[NB_L1*(n-nres)];
            alignas(64) double L1CU[NB_L1*NB_L1];
            alignas(64) double L1CD[NB_L1*NB_L1];

            for(int j = jres; j < n; j++){
            for(int i = ii; i < ii + NB_L1; i++){
                    L1A[i-ii + NB_L1*(j-jres)] = A[i + m*j];
                }}
            for(int z = zz; z < zz + NB_L1; z++){
            for(int i = ii; i < ii + NB_L1; i++){
                    L1CU[i-ii + NB_L1*(z-zz)] = CU[i + m*z];
                    L1CD[i-ii + NB_L1*(z-zz)] = CD[i + m*z];
                }}
            for(int z = 0; z <  NB_L1; z++){
            for(int j = 0; j < n-nres; j++){
                    __m512d b = _mm512_set1_pd(B[j + jres + n*(z + zz)]);
            for(int i = 0; i <  NB_L1; i+=8){
                    __m512d a = _mm512_load_pd(L1A + i + NB_L1*j);
                    __m512d c = _mm512_load_pd(L1CU + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + NB_L1*z, c);
                    c = _mm512_load_pd(L1CD + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + NB_L1*z, c);
            }}}
            for(int z = zz; z < zz + NB_L1; z++){
                for(int i = ii; i < ii + NB_L1; i++){
                    CU[i + m*z] = L1CU[i-ii + NB_L1*(z-zz)];
                    CD[i + m*z] = L1CD[i-ii + NB_L1*(z-zz)];
                }}
        }}
    }
    for (int zres = kres; zres < k; zres+=(k-kres)){
        for (int ii = iii; ii < iii + NB_L2; ii+=NB_L1){
            alignas(64) double L1A[NB_L1*(n-nres)];
            alignas(64) double L1CU[NB_L1*(k-kres)];
            alignas(64) double L1CD[NB_L1*(k-kres)];
            for(int j = jres; j < n; j++){
            for(int i = ii; i < ii + NB_L1; i++){
                    L1A[i-ii + NB_L1*(j-jres)] = A[i + m*j];
                }}
            for(int z = zres; z < k; z++){
            for(int i = ii; i < ii + NB_L1; i++){
                    L1CU[i-ii + NB_L1*(z-zres)] = CU[i + m*z];
                    L1CD[i-ii + NB_L1*(z-zres)] = CD[i + m*z];
                }}
            for(int z = 0; z <  k-kres; z++){
            for(int j = 0; j <  n-nres; j++){
                    __m512d b = _mm512_set1_pd(B[j + jres + n*(z + zres)]);
            for(int i = 0; i <  NB_L1; i+=8){
                    __m512d a = _mm512_load_pd(L1A + i + NB_L1*j);
                    __m512d c = _mm512_load_pd(L1CU + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + NB_L1*z, c);
                    c = _mm512_load_pd(L1CD + i + NB_L1*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + NB_L1*z, c);
            }}}
            for(int z = zres; z < k; z++){
                for(int i = ii; i < ii + NB_L1; i++){
                    CU[i + m*z] = L1CU[i-ii + NB_L1*(z-zres)];
                    CD[i + m*z] = L1CD[i-ii + NB_L1*(z-zres)];
                }}
        }
    }}}
    for (int iii = mres; iii < m; iii+=NB_L2){
    for (int jjj = 0; jjj < nres; jjj+=NB_L2){ 
    #pragma omp parallel for
    for (int zzz = 0; zzz < kres; zzz+=NB_L2){           
        for (int jj = jjj; jj < jjj + NB_L2; jj+=NB_L1){
        for (int zz = zzz; zz < zzz + NB_L2; zz+=NB_L1){
            alignas(64) double L1A[NB_L1*m8just]; 
            alignas(64) double L1CU[NB_L1*m8just];
            alignas(64) double L1CD[NB_L1*m8just];

            for (int j = jj; j < jj + NB_L1; j++){
                for(int i = mres; i < m; i++){
                    L1A[i-mres + m8just*(j-jj)] = A[i + m*j];
            }}
            for(int z = zz; z < zz + NB_L1; z++){
                for(int i = mres; i < m; i++){
                    L1CU[i-mres + m8just*(z-zz)] = CU[i + m*z];
                    L1CD[i-mres + m8just*(z-zz)] = CD[i + m*z];
            }}
            
            for(int z = 0; z <  NB_L1; z++){
            for(int j = 0; j <  NB_L1; j++){
                    __m512d b = _mm512_set1_pd(B[j + jj + n*(z + zz)]);
                    __m512d a,c;
            for(int i = 0; i < NL0-mres; i+=8){  
                    a = _mm512_load_pd(L1A + i + m8just*j);
                    c = _mm512_load_pd(L1CU + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + m8just*z, c);
                    c = _mm512_load_pd(L1CD + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + m8just*z, c);
            }
                a = _mm512_maskz_load_pd( mmask,L1A + NL0-mres + m8just*j);
                c = _mm512_maskz_load_pd( mmask,L1CU + NL0-mres + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CU + NL0-mres + m8just*z, mmask, c);
                c = _mm512_maskz_load_pd( mmask,L1CD + NL0-mres + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CD + NL0-mres + m8just*z, mmask, c);
            }}
            for(int z = zz; z < zz + NB_L1; z++){
                for(int i = mres; i < m; i++){
                    CU[i + m*z] = L1CU[i-mres +m8just*(z-zz)];
                    CD[i + m*z] = L1CD[i-mres +m8just*(z-zz)];
                }}
        }}
    }
    for (int zres = kres; zres < k; zres+=(k-kres)){
        for (int jj = jjj; jj < jjj + NB_L2; jj+=NB_L1){
            alignas(64) double L1A[NB_L1*m8just];
            alignas(64) double L1CU[m8just*(k-kres)];
            alignas(64) double L1CD[m8just*(k-kres)];

            for(int j = jj; j < jj + NB_L1; j++){
            for(int i = mres; i < m; i++){
                    L1A[i-mres + m8just*(j-jj)] = A[i + m*j];
                }}

            for(int z = zres; z < k; z++){
            for(int i = mres; i < m; i++){
                    L1CU[i-mres + m8just*(z-zres)] = CU[i + m*z];
                    L1CD[i-mres + m8just*(z-zres)] = CD[i + m*z];
                }}
            for(int z = 0; z <  k-kres; z++){
            for(int j = 0; j <  NB_L1; j++){
                    __m512d b = _mm512_set1_pd(B[j + jj + n*(z + zres)]);
                    __m512d a,c;
            for(int i = 0; i <  NL0 - mres; i+=8){
                    a = _mm512_load_pd(L1A + i + m8just*j);
                    c = _mm512_load_pd(L1CU + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + m8just*z, c);
                    c = _mm512_load_pd(L1CD + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + m8just*z, c);
            }
                a = _mm512_maskz_load_pd( mmask,L1A + NL0-mres + m8just*j);
                c = _mm512_maskz_load_pd( mmask,L1CU + NL0-mres + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CU + NL0-mres + m8just*z, mmask, c);
                c = _mm512_maskz_load_pd( mmask,L1CD + NL0-mres + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CD + NL0-mres + m8just*z, mmask, c);
            }}
            for(int z = zres; z < k; z++){
                for(int i = mres; i < m; i++){
                    CU[i + m*z] = L1CU[i-mres + m8just*(z-zres)];
                    CD[i + m*z] = L1CD[i-mres + m8just*(z-zres)];
                }}
        }}
    }
    for (int jres = nres; jres < n; jres+=(n-nres)){
           #pragma omp parallel for
    for (int zzz = 0; zzz < kres; zzz+=NB_L2){
        for (int zz = zzz; zz < zzz + NB_L2; zz+=NB_L1){

            alignas(64) double L1A[m8just*(n-nres)];
            alignas(64) double L1CU[m8just*NB_L1];
            alignas(64) double L1CD[m8just*NB_L1];
            for(int j = jres; j < n; j++){
            for(int i = mres; i < m; i++){
                    L1A[i-mres + m8just*(j-jres)] = A[i + m*j];
                }}
            for(int z = zz; z < zz + NB_L1; z++){
            for(int i = mres; i < m; i++){
                    L1CU[i-mres + m8just*(z-zz)] = CU[i + m*z];
                    L1CD[i-mres + m8just*(z-zz)] = CD[i + m*z];
                }}
            for(int z = 0; z < NB_L1; z++){
            for(int j = 0; j < n-nres; j++){
                    __m512d b = _mm512_set1_pd(B[j + jres + n*(z + zz)]);
                    __m512d a,c;
                   int i = 0;
            for(; i <  m8just-8; i+=8){
                    a = _mm512_load_pd(L1A + i + m8just*j);
                    c = _mm512_load_pd(L1CU + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + m8just*z, c);
                    c = _mm512_load_pd(L1CD + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_storeu_pd(L1CD + i + m8just*z, c);
            }
                a = _mm512_maskz_load_pd( mmask,L1A + i + m8just*j);
                
                c = _mm512_maskz_load_pd( mmask,L1CU + i + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CU + i + m8just*z, mmask, c);
                c = _mm512_maskz_load_pd( mmask,L1CD + i + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CD + i + m8just*z, mmask, c);
            }
            }
            for(int z = zz; z < zz + NB_L1; z++){
                for(int i = mres; i < m; i++){
                    CU[i + m*z] = L1CU[i-mres + m8just*(z-zz)];
                    CD[i + m*z] = L1CD[i-mres + m8just*(z-zz)];
                }}
        }}
    for (int zres = kres; zres < k; zres+=(k-kres)){

            alignas(64) double L1A[m8just*(n-nres)];
            alignas(64) double L1CU[m8just*(k-kres)];
            alignas(64) double L1CD[m8just*(k-kres)];
            for(int j = jres; j < n; j++){
            for(int i = mres; i < m; i++){
                    L1A[i-mres + m8just*(j-jres)] = A[i + m*j];
                }}
            for(int z = zres; z < k; z++){
            for(int i = mres; i < m; i++){
                    L1CU[i-mres + m8just*(z-zres)] = CU[i + m*z];
                    L1CD[i-mres + m8just*(z-zres)] = CD[i + m*z];
                }}
            for(int z = 0; z <  k-kres; z++){
            for(int j = 0; j <  n-nres; j++){
                    __m512d b = _mm512_set1_pd(B[j + jres + n*(z + zres)]);
                    __m512d a,c;
            for(int i = 0; i <  NL0 - mres; i+=8){
                    a = _mm512_load_pd(L1A + i + m8just*j);
                    c = _mm512_load_pd(L1CU + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CU + i + m8just*z, c);
                    c = _mm512_load_pd(L1CD + i + m8just*z);
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(L1CD + i + m8just*z, c);
            }
                a = _mm512_maskz_load_pd( mmask,L1A + NL0-mres + m8just*j);
                c = _mm512_maskz_load_pd( mmask,L1CU + NL0-mres + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CU + NL0-mres + m8just*z, mmask, c);
                c = _mm512_maskz_load_pd( mmask,L1CD + NL0-mres + m8just*z);
                c = _mm512_maskz_fmadd_round_pd( mmask,a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                _mm512_mask_store_pd(L1CD + NL0-mres + m8just*z, mmask, c);            
            }}
            for(int z = zres; z < k; z++){
                for(int i = mres; i < m; i++){
                    CU[i + m*z] = L1CU[i-mres + m8just*(z-zres)];
                    CD[i + m*z] = L1CD[i-mres + m8just*(z-zres)];
                }}
        }
    }
}}


void oldudmatmul(int m, int n, int k, double *A, double *B, double *CU, double *CD){
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

            alignas(64) double AA[NB_L2*NB_L2];
            alignas(64) double CCU[NB_L2*NB_L2];
            alignas(64) double CCD[NB_L2*NB_L2];

            for (int ii = iL2; ii < iL2 + NB_L2; ii++){
            for (int kk = kL2; kk < kL2 + NB_L2; kk++){
                CCU[(ii - iL2) + NB_L2*(kk - kL2)] = CU[ii + m*kk];
                CCD[(ii - iL2) + NB_L2*(kk - kL2)] = CD[ii + m*kk];
            }}

            for (int ii = iL2; ii < iL2 + NB_L2; ii++){
            for (int jj = jL2; jj < jL2 + NB_L2; jj++){
                AA[(ii - iL2) + NB_L2*(jj - jL2)] = A[ii + m*jj];
            }}

            for (int jL1 = jL2; jL1 < jL2 + NB_L2; jL1+=NB_L1){
            for (int kL1 = kL2; kL1 < kL2 + NB_L2; kL1+=NB_L1){
            for (int iL1 = iL2; iL1 < iL2 + NB_L2; iL1+=NB_L1){                

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

                    __m512d a = _mm512_load_pd( AA + (ii - iL2) + NB_L2*(jj - jL2) );
                    __m512d c = _mm512_load_pd(CCU + (ii - iL2 ) + NB_L2*(kk - kL2));
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(CCU + (ii - iL2 ) + NB_L2*(kk -kL2), c);

                    c = _mm512_load_pd( CCD + (ii - iL2 ) + NB_L2*(kk - kL2) );
                    c = _mm512_fmadd_round_pd(a, b, c, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
                    _mm512_store_pd(CCD + (ii - iL2 ) + NB_L2*(kk - kL2), c);

                }}}
            }}}
            for (int ii = iL2; ii < iL2 + NB_L2; ii++){
            for (int kk = kL2; kk < kL2 + NB_L2; kk++){
                CU[ii + m*kk] = CCU[(ii - iL2) + NB_L2*(kk - kL2)];
                CD[ii + m*kk] = CCD[(ii - iL2) + NB_L2*(kk - kL2)];
            }}
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