# VCP Library

VCP Library は、偏微分方程式に対する精度保証付き数値計算を支援する
C++11 以降で利用できるヘッダ中心のライブラリです。中心になるのは policy-based design の
密行列クラス `vcp::matrix` で、通常の浮動小数点計算、BLAS/LAPACK を使う
高速計算、別途導入した kv ライブラリによる区間演算を同じ行列インターフェースから扱えます。

- GitHub: <https://github.com/koutasekine/vcp>
- 旧 Web サイト: <https://verified.computation.jp/>
  現在の更新は GitHub に移行しており、この Web サイトは legacy information として
  参照してください。
- 論文: Kouta Sekine, Mitsuhiro T. Nakao, and Shin'ichi Oishi,
  "Numerical verification methods for a system of elliptic PDEs, and their
  software library", NOLTA, IEICE, Vol. 12, No. 1, pp. 41-74, 2021.
  DOI: <https://doi.org/10.1587/nolta.12.41>

## 最初に読む場所

インストールは [docs/install.md](docs/install.md) を読んでください。VCP Library は
ヘッダ中心のライブラリなので、基本的には `vcp/` ディレクトリを配置し、
`vcp/` を含むディレクトリを include path に指定します。

まず `vcp::matrix` を使う場合は、この README の「最小例」と
「Matrix policy の選び方」を読んでください。policy の詳しい違いは
[docs/matrix.md](docs/matrix.md) にまとめています。
行列のファイル保存・読み込みは [docs/file_io.md](docs/file_io.md) に
まとめています。既存形式の `save/load` と、環境間移行向けの
`save_portable/load_portable` を分けて説明しています。

コンパイラ、BLAS/LAPACK、OpenMP の設定例は
[docs/blas_build.md](docs/blas_build.md) にまとめています。Apple Silicon Mac で
丸めモード変更可能な OpenBLAS を使う設定例や、AMD/Intel の Linux 環境で
OpenBLAS を source build する設定例もそちらにあります。
BLAS の丸めモード変更を使う環境では、インストール後に
`test_matrix/Check_pdblas_rounding.cpp` と
`test_matrix/Check_pidblas_rounding.cpp` を実行することを推奨します。
Intel MKL を使う場合は、コンパイル前に oneAPI/MKL の環境設定も必要です。

PDE の基本的な使い方は、`test_PDE/test_Emden.cpp` が入口になります。
PDE の例は手法や論文ごとにファイルが分かれており、
`test_PDE/test_Emden_NOLTA2021.cpp` や
`test_PDE/test_Emden_Numerische_Mathematik2020.cpp` のように論文と紐づいた
構成の例もあります。
これらの PDE の例の中心になるのが `vcp/ldbase.hpp` です。Legendre 多項式を
変形して構成した H^1_0 の基底（Legendre 基底）により、単位立方体領域上の
連立半線形楕円型 PDE について、近似解の計算（Newton 法用の行列生成）から
誤差評価に必要な行列・積分値・各種定数の精度保証付き計算までを
1 つのクラス `vcp::Legendre_Bases_Generator` で扱えます。詳しい使い方は
[docs/ldbase.md](docs/ldbase.md) を参照してください。
Newton 法を独自の非線形問題に使う場合は、[docs/newton.md](docs/newton.md) を
参照してください。最小例は `test_PDE/test_newton1.cpp` です。
行列クラスの基本的な使い方は、`test_matrix/test_matrix.cpp` が小さく読みやすい例です。
Fourier 級数を係数列として扱う `vcp::fourier_series` については
[docs/fourier_series.md](docs/fourier_series.md) を参照してください。

## ディレクトリ構成

| ディレクトリ | 内容 |
| --- | --- |
| `vcp/` | VCP Library 本体のヘッダ |
| `vcp/vblas/` | 丸めモード対応 double BLAS（`rdblas`）と Fortran 互換ラッパー（`dblas`）。`-DUSE_VCP_BLAS` 使用時に `vcp/pdblas.hpp` 等から参照される |
| `vcp/vlapack/` | 丸めモード対応 double LAPACK（`rdlapack`）と Fortran 互換ラッパー（`dlapack`）。`-DUSE_VCP_LAPACK` 使用時に `vcp/pdblas.hpp` 等から参照される |
| `test_matrix/` | 行列クラスの利用例と確認プログラム |
| `test_PDE/` | PDE の精度保証付き数値計算例 |
| `tools/` | ダウンロードや展開を補助する小さなツール |

利用側のコンパイルでは、`vcp/` そのものではなく、`vcp/` を含む
ディレクトリを include path に追加してください。詳しくは
[docs/install.md](docs/install.md) を参照してください。区間演算や
`kv::mpfr` を使う場合は、別途取得した kv ライブラリの include path も
追加してください。

## 依存ライブラリ

最小の `vcp::matrix<double>` は C++11 コンパイラで利用できます。
現在、C++11 から C++20 までのコンパイルで動作を確認しています。
gcc と clang の両方で検証しています。Apple Silicon Mac や AMD 環境でも
検証しています。VCP Library 本体は Intel 専用ではありません。
選択する policy やスカラー型によって、次の依存が必要になります。

| 用途 | 主な依存 |
| --- | --- |
| `vcp::pdblas`, `vcp::pidblas` | BLAS/LAPACK または Intel MKL（`-DUSE_VCP_BLAS -DUSE_VCP_LAPACK` 指定時は不要） |
| `vcp::pddblas` | BLAS/LAPACK または Intel MKL（`-DUSE_VCP_BLAS -DUSE_VCP_LAPACK` 指定時は不要）、別途取得した kv ライブラリ（`kv::dd`） |
| OpenMP による内部並列化 | OpenMP 対応コンパイラ |
| `kv::mpfr` | MPFR |
| 区間演算 | 別途取得した kv ライブラリ |

開発環境、Boost、GMP、MPFR などの基本パッケージの導入例は
[docs/install.md](docs/install.md) にまとめています。

## コンパイル例

以下の例では、最小対応規格として `-std=c++11` を指定しています。
検証済みの環境では、`-std=c++14`、`-std=c++17`、`-std=c++20` も利用できます。

`vcp/` ディレクトリを含む親ディレクトリを include path に追加する例:

```bash
g++ -std=c++11 -I/path/to/libs example.cpp
```

`-I/path/to/libs/vcp` ではなく、`-I/path/to/libs` を指定します。

区間演算を使う場合は、kv ライブラリの include path も追加します。
推奨配置は、`vcp/` と `kv/` を同じ親ディレクトリに置く形です。

```bash
g++ -std=c++11 -I/path/to/libs example.cpp
```

BLAS/LAPACK を使う policy の場合:

```bash
g++ -std=c++11 -I/path/to/libs example.cpp -llapack -lblas
```

OpenMP と MPFR も使う典型的な例:

```bash
g++ -std=c++11 -I/path/to/libs -O3 -DNDEBUG -DKV_FASTROUND example.cpp -llapack -lblas -lmpfr -fopenmp
```

Intel MKL を使う場合は、`-lblas` や `-llapack` ではなく MKL の各ライブラリを
明示的にリンクする構成を推奨します。MKL のインストールは Intel oneAPI の
公式 Installation Guide を参照してください。
macOS、Apple Silicon、AMD/Intel Linux、OpenBLAS、clang、Homebrew GCC の具体例は
[docs/blas_build.md](docs/blas_build.md) を参照してください。

外部 BLAS/LAPACK がない環境では、`-DUSE_VCP_BLAS` と `-DUSE_VCP_LAPACK` を
指定することで、純 C++ 実装の BLAS/LAPACK（`vcp/vblas/dblas.hpp` と
`vcp/vlapack/dlapack.hpp`）を代わりに使用できます。

```bash
g++ -std=c++11 -I/path/to/libs -O3 -DNDEBUG -DKV_FASTROUND \
    -DUSE_VCP_BLAS -DUSE_VCP_LAPACK \
    example.cpp -lmpfr -fopenmp
```

`-DUSE_VCP_BLAS` のみで BLAS だけ、`-DUSE_VCP_LAPACK` のみで LAPACK だけを
置き換えることもできます。両フラグを指定した状態で外部ライブラリもリンクすると
シンボルが重複してリンカエラーになります。詳細は
[docs/matrix.md](docs/matrix.md) の「外部 BLAS/LAPACK を使わない場合」を
参照してください。

## 最小例

これは標準 policy `vcp::mats<double>` を使う例です。

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

## Matrix policy の選び方

`vcp::matrix` は次の形で使います。

```cpp
vcp::matrix<ElementType, Policy>
```

`Policy` を省略すると `vcp::mats<ElementType>` が使われます。

```cpp
vcp::matrix<double> A;
vcp::matrix<double, vcp::mats<double> > B;
```

よく使う選択は次の通りです。

| 目的 | 型 |
| --- | --- |
| まず普通の行列計算をしたい | `vcp::matrix<double>` |
| BLAS/LAPACK で高速に計算したい | `vcp::matrix<double, vcp::pdblas>` |
| `kv::dd`（double-double、約 32 桁）を高速に扱いたい | `vcp::matrix<kv::dd, vcp::pddblas>` |
| 汎用の区間行列を使いたい | `vcp::matrix<kv::interval<T>, vcp::imats<T> >` |
| `kv::interval<double>` を高速に扱いたい | `vcp::matrix<kv::interval<double>, vcp::pidblas>` |
| 依存を抑えた軽量 policy を使いたい | `vcp::matrix<T, vcp::minimats<T> >` |

policy ごとの include、向いている用途、注意点は
[docs/matrix.md](docs/matrix.md) を参照してください。

## 区間行列の例

```cpp
#include <iostream>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/imats.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

int main() {
    typedef kv::interval<double> I;
    typedef vcp::imats<double> Policy;

    vcp::matrix<I, Policy> A, b, x;

    A.rand(10);
    b.rand(10, 1);

    b = A * b;
    x = lss(A, b);

    std::cout << x << std::endl;
}
```

高精度計算を使う場合は、`double` の代わりに `kv::mpfr<N>` などの
kv の高精度型を使います。

## エラー処理

VCP Library のエラーは例外で通知されます。VCP 由来の例外は
`vcp::error` から派生します。

```cpp
#include <vcp/error.hpp>

try {
    x = lss(A, b);
} catch (const vcp::error& e) {
    std::cerr << "VCP error: " << e.what() << std::endl;
}
```

代表的な例外は次の通りです。

| 例外 | 意味 |
| --- | --- |
| `vcp::dimension_error` | 行列サイズの不一致 |
| `vcp::index_error` | 添字や部分行列範囲の誤り |
| `vcp::state_error` | 対応していない行列状態や policy 状態 |
| `vcp::domain_error` | 対称性など数学的仮定の破れ |
| `vcp::verification_error` | 精度保証の検証条件の失敗 |
| `vcp::lapack_error` | LAPACK routine が非ゼロの `info` を返した |
| `vcp::io_error` | ファイル入出力やバイナリ形式のエラー |

BLAS/LAPACK を使う一部の操作は、大きな行列で余分なメモリを使わないため
破壊的に実行されます。そのため、例外が投げられた時点で入力行列の内容が
backend routine によって上書き済みの場合があります。

## OpenMP

コンパイラで OpenMP が有効な場合、VCP 内部の一部ループが並列化されます。
外側のループを利用者側で並列化したい場合など、VCP 内部の OpenMP を
止めたいときは次を定義してください。

```bash
-DVCP_NOMP
```

コンポーネント単位では、次のような局所的な無効化も使えます。

```bash
-DVCP_MATS_NOMP
-DVCP_FOURIER_NOMP
-DVCP_LEGENDRE_NOMP
```

## ライセンス

VCP Library は BSD 3-Clause License で配布されています。詳細は
[LICENSE](LICENSE) を参照してください。
