#include <iostream>
#include <omp.h>

#include <kv/hwround.hpp>

#include <vcp/pdblas.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

#include "udmatmul.hpp"

#include <vcp/vcp_timer.hpp>

int main(void){
    int m = 1800;
    int n = 1500;
    int k = 1900;
    vcp::matrix< double, vcp::pdblas > A, B, CU, CD, DU, DD;

    A.rand(m, n);
    B.rand(n, k);
    CU.zeros(m, k);
    CD.zeros(m, k);

    vcp::time.tic();
    udmatmul( m, n, k, A.data(), B.data(), CU.data(), CD.data() );
    vcp::time.toc();

    //std::cout << A << std::endl;
    //std::cout << B << std::endl;
 


    vcp::time.tic();
    kv::hwround::roundup();
    DU = A*B;
    kv::hwround::rounddown();
    DD = A*B;
    vcp::time.toc();

    std::cout << CU.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;
    std::cout << CD.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;    
    std::cout << DU.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;

/*
    std::cout << CU.submatrix({0,10}, {0,10}) << std::endl;
    std::cout << CD.submatrix({0,10}, {0,10}) << std::endl;
    std::cout << D.submatrix({0,10}, {0,10}) << std::endl;
*/
    //std::cout << D - C << std::endl;

    std::cout << "|| DU - CU || <=" << norminf(DU - CU) << std::endl;
    std::cout << "|| DD - CD || <=" << norminf(DD - CD) << std::endl;
    std::cout << "|| CU - CD || <=" << norminf(CU - CD) << std::endl;
    std::cout << "|| DU - DD || <=" << norminf(DU - DD) << std::endl;
}