# Installation

VCP Library はヘッダ中心のライブラリです。基本的なインストールは、
`vcp/` ディレクトリを利用したい場所に置き、コンパイル時の include path に
`vcp/` を含むディレクトリを指定するだけです。

区間演算や高精度計算で kv ライブラリも使う場合は、`vcp/` と `kv/` を
同じ親ディレクトリに置く配置を推奨します。この配置では include path を
1 つ指定するだけで、`#include <vcp/...>` と `#include <kv/...>` の両方を
解決できます。

`vcp/` ディレクトリそのものを include path に指定しない点に注意してください。

## 基本の考え方

VCP のヘッダは次のように include します。

```cpp
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`matrix.hpp` は `vcp::matrix` 本体のヘッダです。`matrix_assist.hpp` は
行列の利用でよく使う補助機能をまとめて読み込むヘッダで、型変換や
ファイル入出力などを使う場合にも利用します。

そのため、コンパイラには `vcp/` の中ではなく、`vcp/` が置かれている
親ディレクトリを渡します。kv を使う場合も、同じ親ディレクトリに `kv/` を
置いておくと同じ include path で利用できます。

```mermaid
flowchart LR
    A["compiler include path (-I)"] --> B["directory containing vcp/"]
    B --> C["vcp/"]
    C --> D["matrix.hpp"]
    C --> E["matrix_assist.hpp"]
    B --> F["kv/"]
    F --> G["interval.hpp"]
```

## 正しい配置例

例えば次のように配置したとします。

```text
/home/user/libs/
├── vcp/
│   ├── matrix.hpp
│   ├── matrix_assist.hpp
│   ├── mats.hpp
│   └── ...
└── kv/
    ├── interval.hpp
    ├── rdouble.hpp
    └── ...
```

この場合、include path に指定するのは
`/home/user/libs` です。

```bash
g++ -std=c++11 -I/home/user/libs example.cpp
```

`#include <vcp/matrix.hpp>` に対して、コンパイラは次のファイルを探します。

```text
/home/user/libs/vcp/matrix.hpp
```

`#include <kv/interval.hpp>` に対しては、次を探します。

```text
/home/user/libs/kv/interval.hpp
```

## 間違いやすい例

次のように `vcp/` ディレクトリそのものを include path に指定しないでください。

```bash
g++ -std=c++11 -I/home/user/libs/vcp example.cpp
```

この指定だと、`#include <vcp/matrix.hpp>` に対してコンパイラは次を探します。

```text
/home/user/libs/vcp/vcp/matrix.hpp
```

つまり `vcp/` が二重になり、ヘッダを見つけられません。

```mermaid
flowchart TD
    OK["OK: -I/home/user/libs"] --> OKFile["/home/user/libs/vcp/matrix.hpp"]
    NG["NG: -I/home/user/libs/vcp"] --> NGFile["/home/user/libs/vcp/vcp/matrix.hpp"]
```

## kv ライブラリを使う場合

区間演算や高精度計算を使う場合は、VCP Library とは別に kv ライブラリを
取得してください。kv のヘッダは次のように include します。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
```

そのため、kv についても `kv/` ディレクトリそのものではなく、
`kv/` が置かれているディレクトリを include path に指定します。
推奨配置では `vcp/` と `kv/` が同じ親ディレクトリにあるため、
include path は 1 つで済みます。

配置例:

```text
/home/user/libs/
├── vcp/
│   ├── matrix.hpp
│   └── ...
└── kv/
    ├── interval.hpp
    ├── rdouble.hpp
    └── ...
```

この場合のコンパイル例:

```bash
g++ -std=c++11 -I/home/user/libs example.cpp
```

コンパイラから見える必要があるファイルは次の 2 つです。

```text
/home/user/libs/vcp/matrix.hpp
/home/user/libs/kv/interval.hpp
```

## リポジトリ内のサンプルをコンパイルする場合

このリポジトリの `test_matrix/` などでサンプルをコンパイルする場合、
1 つ上のディレクトリが `vcp/` を含んでいます。kv を同じ親ディレクトリに
置いている場合は、`-I..` だけで `vcp/` と `kv/` の両方を参照できます。

```text
project-root/
├── vcp/
├── kv/
└── test_matrix/
    └── test_matrix.cpp
```

そのため、`test_matrix/` の中からコンパイルする例は次のようになります。

```bash
cd test_matrix
g++ -std=c++11 -I.. test_matrix.cpp
```

ここで `..` は `vcp/` と `kv/` を含む親ディレクトリを指します。

## 外部ライブラリが必要な場合

通常の `vcp::matrix<double>` だけなら、VCP のヘッダを配置して include path を
指定すれば利用できます。

一方、選択する policy やスカラー型によっては、追加の外部ライブラリが必要です。

| 用途 | 追加で必要なもの |
| --- | --- |
| 区間演算 | kv ライブラリ |
| `kv::mpfr` | MPFR |
| `vcp::pdblas`, `vcp::pidblas` | BLAS/LAPACK または Intel MKL |
| OpenMP を使うビルド | OpenMP 対応コンパイラとリンクオプション |

## 推奨する基本パッケージ

VCP Library 本体はヘッダを配置するだけで利用できますが、実際の数値計算では
開発環境、Boost、GMP、MPFR を入れておくと便利です。

Ubuntu の例:

```bash
sudo apt install build-essential libboost-all-dev libgmp-dev libmpfr-dev
```

Rocky Linux などの Red Hat 系 Linux の例:

```bash
sudo dnf groupinstall "Development Tools" -y
sudo dnf install boost-devel -y
sudo dnf install gmp-devel -y
sudo dnf install mpfr-devel -y
```

BLAS/LAPACK、Intel MKL、OpenBLAS、OpenMP の設定例は
[blas_build.md](blas_build.md) を参照してください。
