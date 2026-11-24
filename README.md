# VCP Library

VCP Library は、偏微分方程式に対する精度保証付き数値計算を支援する
C++11 以降で利用できるヘッダ中心のライブラリです。policy-based design の
密行列クラス `vcp::matrix` と疎行列クラス `vcp::spmatrix` を中心に、通常の
浮動小数点計算、BLAS/LAPACK を使う高速計算、別途導入した kv ライブラリによる
区間演算を同じインターフェースから扱えます。有限要素法(`vcp/bfem`、`vcp/btfem`)、
Legendre 基底(`vcp/ldbase.hpp`)、Fourier 級数、Newton 法も含みます。

- GitHub: <https://github.com/koutasekine/vcp>
- 旧 Web サイト: <https://verified.computation.jp/>
  現在の更新は GitHub に移行しており、この Web サイトは legacy information として
  参照してください。
- 論文: Kouta Sekine, Mitsuhiro T. Nakao, and Shin'ichi Oishi,
  "Numerical verification methods for a system of elliptic PDEs, and their
  software library", NOLTA, IEICE, Vol. 12, No. 1, pp. 41-74, 2021.
  DOI: <https://doi.org/10.1587/nolta.12.41>

## 最初に読む場所

インストールは [docs/install.md](docs/install.md) を読んでください。

## 各機能の詳細

| 文書 | 内容 |
| --- | --- |
| [docs/matrix.md](docs/matrix.md) | 密行列 `vcp::matrix` と policy |
| [docs/spmatrix.md](docs/spmatrix.md) | 疎行列 `vcp::spmatrix` と policy |
| [docs/file_io.md](docs/file_io.md) | 行列のファイル保存・読み込み |
| [docs/error_handling.md](docs/error_handling.md) | 例外と status の扱い |
| [docs/newton.md](docs/newton.md) | Newton 法 |
| [docs/ldbase.md](docs/ldbase.md) | Legendre 基底による単位立方体領域上の楕円型 PDE |
| [docs/fourier_series.md](docs/fourier_series.md) | Fourier 級数 |
| [docs/bfem.md](docs/bfem.md) | Bernstein 多項式による有限要素法 |
| [docs/btfem.md](docs/btfem.md) | 放物型方程式の全離散 Petrov–Galerkin 法 |
| [docs/kernels/](docs/kernels/) | 利用者が直接呼ばない下請けの層(tblas・tlapack・vblas・tsparse など)と、内部のアルゴリズム |

PDE の例は `test_PDE/` にあります。`test_PDE/test_Emden.cpp` が入口で、
`test_PDE/test_Emden_NOLTA2021.cpp` のように論文と紐づいた例もあります。
行列クラスの小さな例は `test_matrix/test_matrix.cpp` です。

## ディレクトリ構成

| ディレクトリ | 内容 |
| --- | --- |
| `vcp/` | VCP Library 本体のヘッダ |
| `vcp/tblas/`, `vcp/tlapack/` | テンプレート版 BLAS・LAPACK(`vcp::mats2<T>` の下請け) |
| `vcp/vblas/`, `vcp/vlapack/` | 丸めモード対応の double BLAS・LAPACK と、その Fortran 互換ラッパー |
| `vcp/tsparse/` | テンプレート版の疎行列カーネル(疎 LU・LDLᵀ・Cholesky・固有値解法・反復法) |
| `vcp/spmats_base/`, `vcp/spimats_base/`, `vcp/spumar_base/`, `vcp/spcmodumar_base/` | 疎行列 policy の下請け |
| `vcp/bfem/` | Bernstein 多項式による有限要素法 |
| `vcp/btfem/` | 放物型方程式の全離散 Petrov–Galerkin 法 |
| `test_matrix/` | 行列クラスの利用例と確認プログラム |
| `test_PDE/` | PDE の精度保証付き数値計算例 |
| `tools/` | ダウンロードや展開を補助する小さなツール |

## 依存ライブラリ

最小の `vcp::matrix<double>` と `vcp::spmatrix<double>` は、C++11 コンパイラだけで
利用できます。C++11 から C++20 までのコンパイルで、gcc と clang の両方で
動作を確認しています。Apple Silicon Mac や AMD 環境でも検証しており、
Intel 専用ではありません。選択する policy やスカラー型によって、次の依存が
必要になります。

| 用途 | 主な依存 |
| --- | --- |
| `vcp::mats2<T>`、`vcp::spmats<T>`、`vcp::spimats<T>` | 追加依存なし |
| `vcp::pdblas`, `vcp::pidblas` | BLAS/LAPACK または Intel MKL(`-DUSE_VCP_BLAS -DUSE_VCP_LAPACK` 指定時は不要) |
| `vcp::pddblas` | 上に加えて kv ライブラリ(`kv::dd`) |
| `vcp::mats2<double>`・`vcp::mats2<kv::dd>` の特殊化ヘッダ | BLAS/LAPACK または Intel MKL |
| `vcp::spumar` | UMFPACK(SuiteSparse)、ARPACK-NG |
| `vcp::spcmodumar` | UMFPACK・CHOLMOD(SuiteSparse)、ARPACK-NG |
| 区間演算、`kv::dd` | 別途取得した kv ライブラリ |
| `kv::mpfr` | MPFR |
| OpenMP による内部並列化 | OpenMP 対応コンパイラ |

`vcp::spumar` と `vcp::spcmodumar` を使う場合は、外部ライブラリのライセンスに
注意してください([docs/spmatrix.md](docs/spmatrix.md) の
「外部ライブラリのライセンスに関する注意」)。

## コンパイル例

以下の例では `-std=c++20` を指定しています。VCP Library 自体は C++11 以降の
どの規格でも使えるので、`-std=c++17` などに変えてもかまいません。
Ubuntu 22.04 の Boost 1.74 など古い Boost は、uBLAS が C++20 に対応していないため、
kv を使うコードを `-std=c++20` でコンパイルできないことがあります。その場合は
`-std=c++17` を指定するか、新しい Boost を使ってください。

`vcp/` そのものではなく、`vcp/` を含む親ディレクトリを include path に
指定します。kv ライブラリを使う場合は、`vcp/` と `kv/` を同じ親ディレクトリに
置くのが推奨の配置です。

```bash
g++ -std=c++20 -I/path/to/libs example.cpp
```

OpenMP・BLAS/LAPACK・MPFR も使う典型的な例:

```bash
g++ -std=c++20 -I/path/to/libs -O3 -DNDEBUG -DKV_FASTROUND example.cpp -llapack -lblas -lmpfr -fopenmp
```

外部 BLAS/LAPACK がない環境では、`-DUSE_VCP_BLAS -DUSE_VCP_LAPACK` を指定すると
VCP 内部の実装を代わりに使います。MKL、OpenBLAS、macOS などの具体例は
[docs/blas_build.md](docs/blas_build.md) を参照してください。

## 最小例

密行列(既定の policy `vcp::mats2<double>`):

```cpp
#include <iostream>

#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main() {
    vcp::matrix<double> A, b, x;

    A.rand(10);       // 10 x 10 random matrix
    b.rand(10, 1);    // 10 x 1 random vector

    b = A * b;
    x = lss(A, b);

    std::cout << x << std::endl;
}
```

疎行列(既定の policy `vcp::spmats<double>`):

```cpp
#include <iostream>
#include <vector>

#include <vcp/spmatrix.hpp>

int main() {
    const int n = 100;
    vcp::spmatrix<double> A;
    A.resize(n, n);
    for (int i = 0; i < n; i++) {
        A.add(i, i, 4.0);
        if (i + 1 < n) { A.add(i, i + 1, -1.0); A.add(i + 1, i, -1.0); }
    }
    A.finalize();

    std::vector<double> b(n, 1.0);
    std::vector<double> x = A.solve(b);   // 既定は疎 LU 分解

    std::cout << x[0] << std::endl;
}
```

## policy の選び方

| 目的 | 密行列 | 疎行列 |
| --- | --- | --- |
| まず普通の計算をしたい | `vcp::matrix<double>` | `vcp::spmatrix<double>` |
| 区間演算で包含を得たい | `vcp::matrix<kv::interval<double>, vcp::imats<double> >` | `vcp::spmatrix<kv::interval<double>, vcp::spimats<double> >` |
| 外部ライブラリで高速に計算したい | `vcp::matrix<double, vcp::pdblas>` | `vcp::spmatrix<double, vcp::spumar>` |
| `kv::interval<double>` を高速に扱いたい | `vcp::matrix<kv::interval<double>, vcp::pidblas>` | — |

密行列の `vcp::mats<T>` は非推奨です。`vcp::mats2<T>` を使ってください。
policy ごとの include、向いている用途、注意点は
[docs/matrix.md](docs/matrix.md) と [docs/spmatrix.md](docs/spmatrix.md) を
参照してください。

## エラー処理

VCP Library のエラーは例外で知らせます。VCP 由来の例外は `vcp::error` から
派生します。

```cpp
#include <vcp/error.hpp>

try {
    x = lss(A, b);
} catch (const vcp::error& e) {
    std::cerr << "VCP error: " << e.what() << std::endl;
}
```

例外の一覧、status を返す関数、破壊的な操作の注意は
[docs/error_handling.md](docs/error_handling.md) を参照してください。

## OpenMP

コンパイラで OpenMP が有効な場合、VCP 内部の一部ループが並列化されます。
外側のループを利用者側で並列化したい場合など、VCP 内部の OpenMP を
止めたいときは次を定義してください。

```bash
-DVCP_NOMP
```

部分ごとに止めることもできます。

| マクロ | 止める部分 |
| --- | --- |
| `VCP_MATS_NOMP` | 密行列の policy |
| `VCP_TBLAS_NOMP` | tblas |
| `VCP_SPMATS_NOMP` | 疎行列の policy(spmats・spimats とその下請け) |
| `VCP_SPARSE_NOMP` | 疎行列のカーネル(tsparse) |
| `VCP_FOURIER_NOMP` | Fourier 級数 |
| `VCP_LEGENDRE_NOMP` | Legendre 基底 |
| `VCP_BFEM_NOMP` | bfem・btfem |

## ライセンス

VCP Library は BSD 3-Clause License で配布されています。詳細は
[LICENSE](LICENSE) を参照してください。
