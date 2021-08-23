#include <iostream>
#include <omp.h>

#include <vcp/pdblas.hpp>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

#include "udmatmul.hpp"

#include <vcp/vcp_timer.hpp>

int main(void){
    int m = 2000;
    int n = 4000;
    int k = 4567;
    vcp::matrix< double, vcp::pdblas > A, B, CU, CD, D;

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
    D = A*B;
    D = A*B;
    vcp::time.toc();

    std::cout << CU.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;
    std::cout << CD.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;    
    std::cout << D.submatrix({m-11,m-1}, {k-11,k-1}) << std::endl;

/*
    std::cout << CU.submatrix({0,10}, {0,10}) << std::endl;
    std::cout << CD.submatrix({0,10}, {0,10}) << std::endl;
    std::cout << D.submatrix({0,10}, {0,10}) << std::endl;
*/
    //std::cout << D - C << std::endl;

    std::cout << max(max(D - CU)) << std::endl;
    std::cout << max(max(CU - CD)) << std::endl;
}