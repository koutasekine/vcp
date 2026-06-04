# Matrix / Policy Guide

この文書では、VCP Library の `vcp::matrix` と policy の関係を説明します。

基本形は次の通りです。

```cpp
vcp::matrix<ElementType, Policy>
```

`ElementType` は行列要素の型です。`double`、`kv::interval<double>`、
`kv::mpfr<N>` などを指定します。`Policy` は、行列積、線形方程式、逆行列、
固有値計算などをどの実装で行うかを指定します。

## 必要なヘッダ

区間演算を使う例では、VCP Library とは別に取得した kv ライブラリが
include path から見える必要があります。

`matrix.hpp` は `vcp::matrix` 本体と基本演算を定義するヘッダです。
`matrix_assist.hpp` は、行列を使う上でよく利用する補助機能をまとめて
読み込むヘッダです。現在は初期化補助、型変換、ファイル入出力用の
`vcp_fio.hpp` と `vcp_portable_fio.hpp` などを読み込みます。
通常の利用では、`matrix.hpp` の後に `matrix_assist.hpp` を include する
ことを推奨します。

通常の行列計算:

```cpp
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

BLAS/LAPACK を使う近似計算:

```cpp
#include <vcp/pdblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

汎用の区間行列:

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/imats.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

BLAS/LAPACK を利用する `kv::interval<double>` の区間行列:

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/pidblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`vcp::pidblas_fma` を使う場合は、追加で次を include します。

```cpp
#include <vcp/pidblas_fma.hpp>
```

## policy 一覧

| Policy | 典型的な行列型 | 主な用途 |
| --- | --- | --- |
| `vcp::mats<T>` | `vcp::matrix<T, vcp::mats<T> >` | 汎用の密行列計算 |
| `vcp::minimats<T>` | `vcp::matrix<T, vcp::minimats<T> >` | 機能を絞った軽量な密行列計算 |
| `vcp::pdblas` | `vcp::matrix<double, vcp::pdblas>` | BLAS/LAPACK を使う高速な近似計算 |
| `vcp::imats<T>` | `vcp::matrix<kv::interval<T>, vcp::imats<T> >` | 汎用の区間行列計算 |
| `vcp::imats<T, P>` | `vcp::matrix<kv::interval<T>, vcp::imats<T, P> >` | 中点計算などに policy `P` を使う区間行列計算 |
| `vcp::pidblas` | `vcp::matrix<kv::interval<double>, vcp::pidblas>` | BLAS/LAPACK と有向丸めを使う高速な区間 double 計算 |
| `vcp::pidblas_fma` | `vcp::matrix<kv::interval<double>, vcp::pidblas_fma>` | FMA を利用する実験的な区間行列積 |

`vcp::matrix<T>` の default policy は `vcp::mats<T>` です。

```cpp
vcp::matrix<double> A;
vcp::matrix<double, vcp::mats<double> > B;
```

この 2 つは同じ policy を使います。

## policy の選び方

まず動かす、または小さな検証例を書く場合は `vcp::matrix<double>` から
始めるのが簡単です。

```cpp
vcp::matrix<double> A;
```

大きな `double` 行列で速度が必要な場合は `vcp::pdblas` を使います。

```cpp
vcp::matrix<double, vcp::pdblas> A;
```

区間演算による包含が必要な場合は `vcp::imats<T>` を使います。

```cpp
typedef kv::interval<double> I;
vcp::matrix<I, vcp::imats<double> > A;
```

`kv::interval<double>` で速度も必要な場合は `vcp::pidblas` を使います。

```cpp
typedef kv::interval<double> I;
vcp::matrix<I, vcp::pidblas> A;
```

`pidblas` は一部の演算で kv の `hwround` interface を通して
ハードウェア丸めモードを直接変更します。VCP 側ではその範囲を
`hwround_guard` で保護し、演算後に round-to-nearest へ戻します。

## 行列の基本

行列要素は 0 始まりの添字でアクセスします。

```cpp
vcp::matrix<double> A, B, C;

A.zeros(3, 3);
B.eye(3);

A(0, 0) = 2.0;
A(1, 1) = 3.0;
A(2, 2) = 4.0;

C = A + B;
C = A * B;
C = transpose(A);
```

よく使う初期化:

```cpp
A.zeros(n);       // n x n zero matrix
A.zeros(n, m);    // n x m zero matrix
A.ones(n, m);     // n x m matrix filled with 1
A.eye(n);         // n x n identity matrix
A.rand(n, m);     // n x m random matrix
A.clear();        // matrix contents を解放
```

サイズとデータ:

```cpp
int r = A.rowsize();
int c = A.columnsize();
int n = A.elementsize();

double* p = A.data();
const std::vector<double>& v = A.vecpointer();
```

内部の格納は列優先です。`A(i, j)` は通常 `v[i + rowsize() * j]` に
対応します。

高速化のため、`operator()` による直接要素アクセスは境界チェックを
行いません。添字範囲を確認したい場合は、呼び出し側で明示的に確認して
ください。

## 線形方程式

`lss(A, b)` は `A x = b` を解きます。

```cpp
vcp::matrix<double, vcp::pdblas> A, b, x;

A.rand(100);
b.rand(100, 1);

x = lss(A, b);
```

区間行列でも同じ名前で使えます。

```cpp
typedef kv::interval<double> I;

vcp::matrix<I, vcp::pidblas> A, b, x;

A.rand(100);
b.rand(100, 1);

x = lss(A, b);
```

BLAS/LAPACK 系 policy では、大きな行列で余分なメモリを使わないため、
一部の内部 routine が対象行列を破壊的に上書きします。

## 固有値問題

対称行列の固有値:

```cpp
vcp::matrix<double, vcp::pdblas> A, E;

A.rand(100);
A = transpose(A) + A;

eigsym(A, E);
```

固有ベクトルも必要な場合:

```cpp
vcp::matrix<double, vcp::pdblas> A, E, V;

A.rand(100);
A = transpose(A) + A;

eigsym(A, E, V);
```

一般化対称固有値問題:

```cpp
vcp::matrix<double, vcp::pdblas> A, B, E;

A.rand(100);
B.rand(100);

A = transpose(A) + A;
B = ltransmul(B);

eigsymge(A, B, E);
```

## 区間演算での注意

`x` が区間になり得る場所では、`pow(x, 2)` を機械的に `x * x` に
置き換えないでください。

区間演算では、`x * x` は依存性による過大評価を受けることがあります。
一方で、区間の累乗関数は偶数乗の構造を利用できる場合があります。
そのため VCP では、区間値が入り得る数値計算上の `pow(x, 2)` を
意図的に残しています。

一方、配列サイズ、ループ回数、添字計算などの整数計算では `std::pow` を
使わず、整数演算で計算します。

## BLAS/LAPACK policy の注意

`vcp::pdblas` と `vcp::pidblas` は、選択された演算で BLAS/LAPACK routine を
呼びます。

重要な点:

- `pdblas` の要素型は `double` です。
- `pidblas` の要素型は `kv::interval<double>` です。
- 利用プログラムは BLAS/LAPACK または MKL とリンクする必要があります。
- LAPACK routine が非ゼロの `info` を返した場合、`vcp::lapack_error` を投げます。
- 一部の演算は破壊的で、例外発生時に入力行列がすでに上書きされている場合があります。

大きな行列を扱う用途を優先し、VCP は失敗時に必ず元の行列を保存するための
丸ごとのコピーを標準では行いません。

MKL の推奨リンクオプションや、Apple Silicon Mac で丸めモード変更可能な
OpenBLAS を使う設定例は [blas_build.md](blas_build.md) を参照してください。

Newton 法で非線形方程式を解く場合は、`matrix` で残差ベクトルと
ヤコビ行列を作り、`vcp::Newton` に渡す形になります。詳しくは
[newton.md](newton.md) を参照してください。

## ファイル入出力

行列の保存と読み込みには、従来形式の `save/load` と、環境間移行向けの
`save_portable/load_portable` があります。

```cpp
vcp::save(A, "A");
vcp::load(B, "A");

vcp::save_portable(A, "A");
vcp::load_portable(B, "A");
```

従来形式は既存の `.matrix_*` ファイルとの互換性を保つ形式です。
portable 形式は `.vcpmat` に保存し、型情報、payload byte size、CRC64 を
持つ binary 形式です。詳細は [file_io.md](file_io.md) を参照してください。

## エラー処理

Matrix 関連のエラーは `vcp::error` 派生の例外で通知されます。

```cpp
#include <vcp/error.hpp>

try {
    x = lss(A, b);
} catch (const vcp::dimension_error& e) {
    std::cerr << "dimension error: " << e.what() << std::endl;
} catch (const vcp::error& e) {
    std::cerr << "VCP error: " << e.what() << std::endl;
}
```

| 例外 | 意味 |
| --- | --- |
| `vcp::dimension_error` | 行列サイズの不一致 |
| `vcp::index_error` | 添字や部分行列範囲の誤り |
| `vcp::state_error` | 対応していない行列状態や policy 状態 |
| `vcp::domain_error` | 対称性など数学的仮定の破れ |
| `vcp::verification_error` | 精度保証の検証条件の失敗 |
| `vcp::backend_error` | 外部 backend や下位 routine のエラー |
| `vcp::blas_error` | BLAS 関連のエラー |
| `vcp::lapack_error` | LAPACK routine が非ゼロの `info` を返した |
| `vcp::io_error` | ファイル入出力やバイナリ形式のエラー |

## OpenMP

コンパイラが OpenMP を有効にしている場合、VCP 内部の一部ループは
OpenMP で並列化されます。

全体で無効化する場合:

```bash
-DVCP_NOMP
```

Matrix、Fourier、Legendre など、部分的に無効化する場合:

```bash
-DVCP_MATS_NOMP
-DVCP_FOURIER_NOMP
-DVCP_LEGENDRE_NOMP
```

外側の利用者コードで OpenMP 並列化を行う場合は、内側の VCP ループを
無効化することで nested OpenMP を避けられます。

## 関連するサンプル

| ファイル | 内容 |
| --- | --- |
| `test_matrix/test_matrix.cpp` | 基本的な行列操作 |
| `test_matrix/test_interval_matrix.cpp` | 区間行列の例 |
| `test_matrix/Check_OpenMP.cpp` | OpenMP 関連の確認 |
| `test_matrix/Check_pdblas_rounding.cpp` | `pdblas` と丸めの確認 |
| `test_matrix/Check_pidblas_rounding.cpp` | `pidblas` と丸めの確認 |
| `test_PDE/test_newton1.cpp` | `vcp::Newton` の最小例 |
| `test_PDE/test_newton2_Emden.cpp` | Legendre 基底を使う Newton 法の例 |
| `test_PDE/test_Emden.cpp` | PDE の基本的な使い方の例 |
| `test_PDE/test_Emden_NOLTA2021.cpp` | NOLTA 2021 論文に紐づく Emden 問題の例 |
| `test_PDE/test_Emden_Numerische_Mathematik2020.cpp` | Numerische Mathematik 2020 論文に紐づく Emden 問題の例 |
| `test_PDE/test_Allen_Cahn.cpp` | 非線形 PDE の例 |
