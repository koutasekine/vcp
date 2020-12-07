# VCP Library
[VCP Library](https://verified.computation.jp/)

The realization of computer-assisted existence proof of solutions cannot be completed without a computer environment.
However, various kinds of techniques are required to estimate the errors that occur in all computations.
Additionally, the computer-assisted proof for existence of solutions of PDEs requires numerical accuracy with a reasonable time of computation.
For this purpose, the Verified Computation for PDEs (VCP) library is introduced as a software library for the computer-assisted existence proof of solutions to PDEs.
The VCP library is a software library developed by the first author in the C++ programming language.
A feature of the matrix class of the VCP library is that it can be integrated with policy-based design, for example,
  * high-speed approximate computation by Intel(R) MKL with double-data type
  * high precision approximate computation using MPFR
  * numerical linear algebra with guaranteed accuracy using the above data type combined with the kv library

Additionally, because the VCP library has extensibility, which is one of the features of policy-based design, it is designed to withstand the computer-assisted proof of PDEs.

We require using the following libraries in the VCP library:
  * [BLAS and Lapack](http://www.netlib.org/lapack/)
  * [MPFR](https://www.mpfr.org/): multiple-precision floating-point number library
  * [kv library](http://verifiedby.me/kv/index-e.html): numerical verification library with guaranteed accuracy written in the C++ programming language

E-mail:
k.sekine@computation.jp

Kouta Sekine
