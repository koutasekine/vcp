# VCP Library

## About
The realization of computer-assisted existence proof to solutions of partial differential equations (PDEs) cannot be completed without a computer environment.
However, various kinds of techniques are required to estimate the errors that occur in all computations.
Additionally, the computer-assisted proof for existence of solutions of PDEs requires numerical accuracy with a reasonable time of computation.
For this purpose, the Verified Computation for PDEs (VCP) library is introduced as a software library for the computer-assisted existence proof of solutions to PDEs.
The VCP library is a software library developed by the first author in the C++ programming language.
A feature of the matrix class of the VCP library is that it can be integrated with policy-based design, for example,
  * high-speed approximate computation by Intel(R) MKL with double-data type
  * high precision approximate computation using MPFR
  * numerical linear algebra with guaranteed accuracy using the above data type combined with the kv library

Additionally, because the VCP library has extensibility, which is one of the features of policy-based design, it is designed to withstand the computer-assisted proof of PDEs.

VCP library comprises of
  * Matrix class
  * Newton method class
  * Legendre basis class (for elliptic PDEs)
  * Fourier series class (for delay ODEs)
  * Gauss-Legendre method: implicit runge kutta method (for ODEs and parabolic PDE, **Not verification**) 


We require using the following libraries in the VCP library:
  * [BLAS and Lapack](http://www.netlib.org/lapack/)
  * [MPFR](https://www.mpfr.org/): multiple-precision floating-point number library
  * [kv library](http://verifiedby.me/kv/index-e.html): numerical verification library with guaranteed accuracy written in the C++ programming language

Homepage (Japanese) [VCP Library](https://verified.computation.jp/)<br>
To cite VCP Library, please use <br>
*Kouta Sekine, Mitsuhiro T. Nakao, and Shin'ichi Oishi: “Numerical verification methods for a system of elliptic PDEs, and their software library”, NOLTA, IEICE, Vol.12, No.1, Jan., 2021.*

E-mail:
k.sekine@computation.jp

Kouta Sekine

## How to install
### Directory configuration
```
 ~/any/
  　├ vcp/
  　├ test_matrix/
  　├ test_PDE/
  　└ <your_folder>/
            └ <your_file.cpp>
```

### Directory configuration with kv library
```
 ~/any/
  　├ kv/
  　├ test/
　  ├ example/
  　├ vcp/
  　├ test_matrix/
  　├ test_PDE/
  　└ <your_folder>/
            └ <your_file.cpp>
```


### Ubuntu 20.04 (also compatible with WSL2 or Docker for Mac)

##### Recommended packages
```
sudo apt update
yes | sudo apt upgrade
yes | sudo apt install build-essential
yes | sudo apt install libboost-all-dev
yes | sudo apt install libgmp-dev
yes | sudo apt install libmpfr-dev
yes | sudo apt install liblapack-dev
yes | sudo apt install intel-mkl
```

##### Compile options
Minimum Compilation Options:<br>
`g++ -I.. <filename.cpp>`

Recommended Compile options with BLAS, Lapack, mpfr, OpenMP and kv library:<br>
`g++ -I.. -DNDEBUG -DKV_FASTROUND -O3 <filename.cpp> -llapack -lblas -lmpfr -fopenmp`

Note that compile options `-DNDEBUG` and `-DKV_FASTROUND` are optional for kv library.

##### Check if the rounding mode of BLAS's dgemm is changeable using the kv library
The rounding mode of BLAS is important for rigorous computing.
You have to make sure of the rounding mode of BLAS after installation.

```bash check_round.sh
  echo "######################################################"
  echo "Check for BLAS's rounding mode changes: Please wait..."
  cd "<folder_path>/test_matrix/"
  g++ -I.. -DNDEBUG -DKV_FASTROUND -O3 Check_pdblas_rounding.cpp -llapack -lblas -lmpfr -fopenmp && ./a.out
  g++ -I.. -DNDEBUG -DKV_FASTROUND -O3 Check_pidblas_rounding.cpp -llapack -lblas -lmpfr -fopenmp && ./a.out

  echo
  echo "######################################################"
  echo "Check for OpenMP's rounding mode changes: Please wait..."
  g++ -I.. -DNDEBUG -DKV_FASTROUND -O3 Check_OpenMP.cpp -llapack -lblas -lmpfr -fopenmp && ./a.out
```

If the rounding cannot be changed
 * Do not use the `pidblas` policy
 * Please check:
      - `sudo update-alternatives --config libblas.so-x86_64-linux-gnu`
      - `sudo update-alternatives --config liblapack.so-x86_64-linux-gnu`


### Ubuntu 16.04, 18.04, Centos 6, Centor 7

##### Recommended packages
```
sudo apt update
yes | sudo apt upgrade
yes | sudo apt install build-essential
yes | sudo apt install libboost-all-dev
yes | sudo apt install libgmp-dev
yes | sudo apt install libmpfr-dev
yes | sudo apt install liblapack-dev
```
with [MKL liblary](https://software.intel.com/content/www/us/en/develop/tools/performance-libraries.html).

How to install MKL library
  1) Download the MKL
  2) Extract the compressed folder of the MKL.
  3) Open the file of MKL, silent.cfg, with a text editor.
  4) Change `ACCEPT_EULA=decline` to `ACCEPT_EULA=accept` and close the editor.
  5) Run `sudo . /install.sh -s silent.cfg`
  6) Run `source /opt/intel/bin/compilervars.sh intel64`
  7) Run `source /opt/intel/mkl/bin/mklvars.sh intel64 lp64`


##### Compile options
Minimum Compilation Options:<br>
`g++ -I.. <filename.cpp>`

Recommended Compile options with MKL, mpfr, OpenMP and kv library:<br>
`g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 <filename.cpp> -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp`


## How to use VCP's matrix class
The matrix class of the VCP library is based on a policy using a template.
Therefore, when declaring the matrix, it is necessary to determine **data type for elements** and **algorithm policy** as template arguments.
The VCP library provides the following four Algorithm policies:

 * `vcp::mats< _T >`: performs approximate computations on general data type `_T`.
 * `vcp::pdblas`: performs fast approximate computations on `double` data type using BLAS and Lapack.
 * `vcp::imats< _T >`: performs verified computations on general data type `_T`.
 * `vcp::pidblas`: performs fast verified computations on `double` data type using BLAS and Lapack.


For example, we can use the VCP library as shown in Source code 1.
Line 6 in Source code \ref{basic_matrix1} declares matrix `A, b, x` with `double`-type and `vcp::mats<double>`-algorithm.
Random matrix `A` and random vector `b` are created on lines 7 and 8, respectively.
On line 9, the result of `A` times `b` is substituted for `b`.
Line 10 uses the function `lss` to determine the solution `x` of simultaneous linear equation `Ax = b`.
Finally, line 11 displays the solution `x`.


```cpp Source_code_1.cpp {.line-number .copy}
#include <iostream>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main(void){
   vcp::matrix< double, vcp::mats<double> > A, b, x;
   A.rand(10); // Create a 10*10 random matrix
   b.rand(10,1); // Create a 10*1 random vector
   b = A*b;  // Compute A times b
   x = lss(A, b); // Solve Ax = b
   std::cout << x << std::endl; // Display x
}
```

The feature of the VCP's matrix class is that it can be changed to verify a numerical computation by changing the `kv::interval< _T >`-data type and **algorithm policy**.
For example, when verifying a computation with high precision, the VCP's matrix class can be written as Source code 2.
Line 5 in Source code 2 declares matrix `A, b, x` with `kv::interval< kv::mpfr<300> >`-type and `vcp::imats< kv::mpfr<300> >`-algorithm, where `kv::mpfr<300>` is an MPFR type with a mantissa part of 300 bits, and `kv::interval< _T>` is kv's interval arithmetic class.
Because `vcp::imats< _T >` is selected as the algorithm policy in line 12, the matrix computation is an algorithm with guaranteed accuracy.
Therefore, lines 13 to 17 in Source code 2 are the same as lines 7 to 11 in Source code 1, but the results of the matrix-vector product on line 15 and the solution `x` on line 16 in Source code 2 contain exact solutions.

```cpp Source_code_2.cpp {.line-number .copy}
#include <iostream>

#include <kv/interval.hpp> // kv library
#include <kv/mpfr.hpp> // kv library
#include <kv/rmpfr.hpp> // kv library

#include <vcp/imats.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main(void){
   vcp::matrix< kv::interval< kv::mpfr<300> >, vcp::imats<kv::mpfr<300> > > A, b, x;
   A.rand(10); // Create a 10*10 random matrix
   b.rand(10,1); // Create a 10*1 random vector
   b = A*b;  // Compute A times b
   x = lss(A, b); // Solve Ax = b
   std::cout << x << std::endl; // Display x
}
```

Additionally, if we declare a matrix with `vcp::matrix<double, vcp::pdblas>`, it can be changed to an algorithm using BLAS and Lapack.
When using high-speed numerical computation with guaranteed accuracy using BLAS and Lapack, we declare a matrix similar to `vcp::matrix< kv::interval<double>, vcp::pidblas>`.
For other functions of VCP's matrix class, see:

```cpp
  //(1) Matrix initialization
  int n = 10;
  int m = 5;
  A.zeros(n); // Create an n*n zero matrix
  A.zeros(n, m); // Create an n*m zero matrix 
  A.ones(n); // Create an n*n matrix with all elements 1
  A.ones(n, m); // Create an n*m matrix with all elements 1 
  A.rand(n); // Create an n*n random matrix 
  A.rand(n, m); // Create an n*m random matrix 
  A.eye(n); // Create an n*n identity matrix 
  B.eye(n);

  //(2) Obtain the matrix size
  int row = A.rowsize();
  int column = A.columnsize();

  //(3) Access the elements
  A(0,0) = 10;
  A(5,3) = A(0,0); 

  //(4) Arithmetic operations
  C = 1 + A; // For all elements
  C = A + 2; // For all elements
  C = A + B; // Matrix addition
  C += A; // C = C + A

  C = 1 - A; // For all elements
  C = A - 2; // For all elements
  C = A - B; // Matrix subtraction
  C -= A; // C = C - A

  C = 1 * A; // For all elements
  C = A * 2; // For all elements
  C = A * B; // Matrix multiplication
  C *= A; // C = C*A
  C = ltransmul(A); // C = transpose(A)*A 

  A = A + 1; 
  C = 1 / A; // For all elements
  C = A / 2; // For all elements

  //(5) Elementary functions, etc. (Element wise)
  C = abs(A);
  C = sqrt(A);
  C = sin(A);
  C = cos(A);
  C = exp(A);
  C = log(A);

  //(6) MATLAB-like functions
  C = sum(A);
  C = diag(A);
  C = transpose(A);
  C = max(A);
  C = min(A);
  C = normone(A);
  C = norminf(A);

  //(7) Solve linear system Ax = b
  A.rand(n);
  b.ones(n,m);
  b = A*b;
  x = lss(A,b); // Find x s.t. Ax=b

  //(8) Solve symmetric eigenvalue problem Ax = lambda x
  A = transpose(A) + A; 
  eigsym(A, C); // Eigenvalue for the diagonal part of matrix C

  //(9) Solve generalized symmetric eigenvalue problem Ax = lambda B x
  B.rand(n); 
  A = transpose(A) + A; 
  B = ltransmul(B); 
  eigsymge(A, B, C); // Eigenvalue for the diagonal of matrix C

  //(10) Display matrix
  std::cout << "Matrix C = \n" << C << std::endl; //

  //(11) Concatenate the matrix
  D = vercat(A,B,C); // Vertically (similar to MATLAB's [A;B;C])
  D = horzcat(A,B,C); // Horizontally (similar to MATLAB's [A,B,C])
  C.rand(n);
  A = vercat(D, horzcat(A,B,C));

  //(12) Release the matrix 
  A.clear();
```