#include <iostream>
#include <vector>
#include <vcp/vcp_timer.hpp>
#include <vcp/ifma.hpp>
#include <random>

double randn(void){
    std::random_device seed_gen;
    std::mt19937 engine(seed_gen());
    std::normal_distribution< double > dist(0.0, 1.0);    
    return std::abs(dist(engine));
}

int main(void){
    double resup, resdown;

    int N = 64;
    int M = 100000000;

    std::vector< double > au, bu, cu, ad, bd, cd;
    au.resize(N);
    bu.resize(N);
    cu.resize(N);

    ad.resize(N);
    bd.resize(N);
    cd.resize(N);


// 2000000
    for (int i = 0; i < N; i++){
        au[i] = randn();
        bu[i] = randn();
        cu[i] = randn();

//        ad[i] = -vcp::randn();
//        bd[i] = -vcp::randn();
//        cd[i] = -vcp::randn();
        ad[i] = au[i] - au[i]/10;
        bd[i] = bu[i] - bu[i]/10;
        cd[i] = cu[i] - cu[i]/10;
    }

    std::cout << "N = " << N << std::endl;
    std::cout << "M = " << M << std::endl;

    std::cout << "a = [ " << ad[N-1] << ", " << au[N-1] << " ]" << std::endl;
    std::cout << "b = [ " << bd[N-1] << ", " << bu[N-1] << " ]" << std::endl;
    std::cout << "c = [ " << cd[N-1] << ", " << cu[N-1] << " ]" << std::endl;    
    

    std::cout << "===============================================" << std::endl;
    
    #ifdef VCP_USE_AVX512
        std::cout << "Use ifma with AVX512" << std::endl;
    #else
        std::cout << "Use ifma with std::fma" << std::endl;
    #endif

    vcp::time.tic();
    for (int j = 0; j < M; j++){
    for (int i = 0; i < N; i++){        
        vcp::ifma( au[i], ad[i], bu[i], bd[i], cu[i], cd[i], resup, resdown ); 
    }}
    vcp::time.toc();
    
    std::cout << "[" << resdown << ", " << resup << "]" << std::endl;
}