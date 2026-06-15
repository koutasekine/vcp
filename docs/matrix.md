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

BLAS/LAPACK を使う `kv::dd` の高速近似計算:

```cpp
#include <kv/dd.hpp>

#include <vcp/pddblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`vcp::mats2<T>`（tblas/tlapack 汎用実装。追加依存なし）:

```cpp
#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`vcp::mats2<double>`（double BLAS/LAPACK 特殊化を有効にした場合）:

```cpp
// tlapack_double.hpp は内部で tblas_double.hpp も include する
#include <vcp/tlapack/tlapack_double.hpp>

#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`vcp::mats2<kv::dd>`（kv::dd 特殊化を有効にした場合）:

```cpp
#include <kv/dd.hpp>

// tlapack_dd.hpp は内部で tblas_double.hpp, tblas_dd.hpp,
// tlapack_double.hpp も include するため、これだけで十分
#include <vcp/tlapack/tlapack_dd.hpp>

#include <vcp/mats2.hpp>
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
| `vcp::pddblas` | `vcp::matrix<kv::dd, vcp::pddblas>` | double の BLAS/LAPACK を使う `kv::dd` の高速な近似計算 |
| `vcp::mats2<T>` | `vcp::matrix<T, vcp::mats2<T> >` | tblas/tlapack を下請けとして使う近似計算。インクルード構成で特性が変わる |
| `vcp::imats<T>` | `vcp::matrix<kv::interval<T>, vcp::imats<T> >` | 汎用の区間行列計算 |
| `vcp::imats<T, P>` | `vcp::matrix<kv::interval<T>, vcp::imats<T, P> >` | 中点計算などに policy `P` を使う区間行列計算 |
| `vcp::pidblas` | `vcp::matrix<kv::interval<double>, vcp::pidblas>` | BLAS/LAPACK と有向丸めを使う高速な区間 double 計算 |
| `vcp::pidblas_fma` | `vcp::matrix<kv::interval<double>, vcp::pidblas_fma>` | FMA を利用する実験的な区間行列積 |

`vcp::matrix<T>` の default policy は `vcp::mats<T>` です。

`vcp::mats2<T>` は tblas/tlapack を下請けとして使う近似計算 policy で、
`vcp::mats<T>` を継承し、行列積・線形方程式・逆行列・Cholesky 分解・
対称固有値問題・一般化対称固有値問題をオーバーライドしています。
インクルードする特殊化ヘッダによって `mats<T>` 相当から `pdblas`・`pddblas`
相当まで動作が切り替わります。詳細は後述の「`mats2<T>` の詳細」節を参照してください。

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

`kv::dd`（double-double、約 32 桁）で速度が必要な場合は `vcp::pddblas` を
使います。アルゴリズムと制限は後述の
「`vcp::pddblas` の詳細」を参照してください。

```cpp
vcp::matrix<kv::dd, vcp::pddblas> A;
```

区間演算による包含が必要な場合は `vcp::imats<T>` を使います。

```cpp
typedef kv::interval<double> I;
vcp::matrix<I, vcp::imats<double> > A;
```

`vcp::imats<T, P>` の第 2 テンプレート引数 `P` には、精度保証付き
数値計算の内部で使う近似計算用の policy を指定します（省略時は
`vcp::mats<T>`）。`P` が使われるのは、`lss` の近似逆行列や
`eigsym`/`eigsymge` の近似固有対といった近似計算だけで、包含の計算は
すべて `kv::interval<T>` の区間演算で行われます。そのため `P` の選択は
包含の成否・包含の幅・計算速度にだけ影響し、精度保証の正しさには
影響しません。

たとえば `kv::interval<kv::dd>` の精度保証付き数値計算で、近似計算に
`vcp::pddblas` を使う場合は次のようにします。

```cpp
typedef kv::interval<kv::dd> I;
vcp::matrix<I, vcp::imats<kv::dd, vcp::pddblas> > A;
```

特に固有値問題では、近似固有対の精度が包含の幅に直結するため、
既定の `vcp::imats<kv::dd>` より高速になるだけでなく、包含幅も
小さくなります。詳細と制限は後述の
「`vcp::pddblas` の詳細」を参照してください。

`kv::interval<double>` で速度も必要な場合は `vcp::pidblas` を使います。

```cpp
typedef kv::interval<double> I;
vcp::matrix<I, vcp::pidblas> A;
```

`pidblas` は一部の演算で kv の `hwround` interface を通して
ハードウェア丸めモードを直接変更します。VCP 側ではその範囲を
`hwround_guard` で保護し、演算後に round-to-nearest へ戻します。

tblas/tlapack ベースの `vcp::mats2<T>` は、インクルードする特殊化ヘッダに
よって挙動が変わる点が特徴です。特殊化ヘッダなしでは追加の依存ライブラリ
なしに `mats<T>` 相当の計算ができ、特殊化ヘッダを追加することで
`pdblas`・`pddblas` 相当の高速計算に切り替わります。
`kv::mpfr<N>` など `double`・`kv::dd` 以外の型でも同じ policy インターフェース
を使える点も `pdblas`・`pddblas` との違いです。詳細は後述の
「`mats2<T>` の詳細」を参照してください。

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

`vcp::pdblas`、`vcp::pddblas`、`vcp::pidblas` は、選択された演算で
BLAS/LAPACK routine を呼びます。

重要な点:

- `pdblas` の要素型は `double` です。
- `pddblas` の要素型は `kv::dd` です（次節参照）。
- `pidblas` の要素型は `kv::interval<double>` です。
- 利用プログラムは BLAS/LAPACK または MKL とリンクする必要があります。外部ライブラリを使わない場合は後述の「外部 BLAS/LAPACK を使わない場合」を参照してください。
- LAPACK routine が非ゼロの `info` を返した場合、`vcp::lapack_error` を投げます。
- 一部の演算は破壊的で、例外発生時に入力行列がすでに上書きされている場合があります。

大きな行列を扱う用途を優先し、VCP は失敗時に必ず元の行列を保存するための
丸ごとのコピーを標準では行いません。

MKL の推奨リンクオプションや、Apple Silicon Mac で丸めモード変更可能な
OpenBLAS を使う設定例は [blas_build.md](blas_build.md) を参照してください。

Newton 法で非線形方程式を解く場合は、`matrix` で残差ベクトルと
ヤコビ行列を作り、`vcp::Newton` に渡す形になります。詳しくは
[newton.md](newton.md) を参照してください。

## 外部 BLAS/LAPACK を使わない場合 (`-DUSE_VCP_BLAS` / `-DUSE_VCP_LAPACK`)

`-DUSE_VCP_BLAS` または `-DUSE_VCP_LAPACK` をコンパイル時に指定すると、
`vcp::pdblas`・`vcp::pidblas`・`vcp::pddblas` が外部 BLAS/LAPACK ライブラリを
必要とせず、`vcp/vblas/dblas.hpp` と `vcp/vlapack/dlapack.hpp` の純 C++ 実装
（丸めモード対応 BLAS/LAPACK）を内部で使用します。これらの policy は
共通ヘッダー `vcp/dblas_dlapack.hpp` を経由して、BLAS/LAPACK 宣言と
backend 選択を共有しています。

| フラグ | 効果 |
| --- | --- |
| `-DUSE_VCP_BLAS` | BLAS 関数（`dgemm_` 等）を `vcp/vblas/dblas.hpp` で提供。外部 BLAS リンク不要 |
| `-DUSE_VCP_LAPACK` | LAPACK 関数（`dgesv_`、`dsyev_` 等）を `vcp/vlapack/dlapack.hpp` で提供。外部 LAPACK リンク不要 |
| 両方 | 外部 BLAS/LAPACK のどちらも不要 |

MKL や OpenBLAS がない環境でも次のようにコンパイルできます。

```bash
g++ -std=c++11 -I/path/to/libs -O3 -DNDEBUG -DKV_FASTROUND \
    -DUSE_VCP_BLAS -DUSE_VCP_LAPACK \
    example.cpp -lmpfr -fopenmp
```

`-DUSE_VCP_BLAS` のみで BLAS だけ、`-DUSE_VCP_LAPACK` のみで LAPACK だけを
置き換えることもできます。

注意点:

- 両フラグを指定した状態で外部 BLAS/LAPACK もリンクすると、シンボルが重複して
  リンカエラーになります。外部ライブラリとの混在はできません。
- `vcp/vblas/dblas.hpp` と `vcp/vlapack/dlapack.hpp` は `vcp/` 内の `vblas/`、`vlapack/`
  サブディレクトリに配置されています。`-I` で親ディレクトリを指定すれば
  `<vcp/vblas/dblas.hpp>` としてインクルードできます。
- `vcp/dblas_dlapack.hpp` は `pdblas` / `pidblas` / `pddblas` の内部実装用の
  共通ヘッダーです。通常の利用では `vcp/pdblas.hpp`、`vcp/pidblas.hpp`、
  `vcp/pddblas.hpp` をインクルードしてください。
- VCP Library の精度保証付き計算（`vcp::pidblas`）において、有向丸めで呼ばれる
  BLAS 関数は `dgemm_` のみです。`vcp/vblas/dblas.hpp` の `dgemm_` は呼び出し時点の
  丸めモードを自動的に読み取るため、`kv::hwround::roundup()` / `rounddown()` で
  丸めモードを設定してから呼び出す既存のコードと正しく連携します。

## `vcp::pddblas` の詳細

`vcp::pddblas` は、要素型 `kv::dd`（double-double、仮数部約 106 bit、
約 32 桁）の行列計算を、double の BLAS/LAPACK routine を使って
高速かつ高精度に行う近似計算 policy です。`vcp::mats<kv::dd>` を継承し、
`vcp::pdblas` と同じメソッド（行列積、線形方程式、逆行列、Cholesky 分解、
固有値問題、一般化固有値問題）を override しています。

```cpp
#include <kv/dd.hpp>

#include <vcp/pddblas.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

vcp::matrix<kv::dd, vcp::pddblas> A, b, x;
```

### 使用しているアルゴリズム

| 演算 | アルゴリズム |
| --- | --- |
| 行列積 `A * B`, `ltransmul(A)` | 尾崎スキーム（error-free transformation による行列分割） |
| 線形方程式 `lss(A, b)` | double の LU 分解（`dgetrf`/`dgetrs`）+ Newton 法型の残差反復 |
| 逆行列 `inv(A)` | 単位行列を右辺とする上記の残差反復 |
| Cholesky 分解 `Cholesky(A)` | double の `dpotrf` + 残差反復による Cholesky factor の改良 |
| 対称固有値問題 `eigsym` | double の `dsyev` + 荻田-相島の反復改良法 |
| 一般化対称固有値問題 `eigsymge` | double の `dsygv` + 荻田-相島の反復改良法（B 内積版） |

行列積は尾崎スキーム
（K. Ozaki, T. Ogita, S. Oishi, S. M. Rump,
"Error-free transformations of matrix multiplication by using fast routines
of matrix multiplication and its applications",
Numerical Algorithms, Vol. 59, pp. 95-118, 2012）を `kv::dd` に
適用したものです。dd 行列を行ごと・列ごとのスケーリングで複数の double
行列（スライス）の和に分割し、スライス同士の積を `dgemm` で計算します。
分割幅はスライス積の `dgemm` が丸め誤差なしで計算されるように選んであり、
誤差ゼロの `dgemm` 結果を dd 演算で総和することで、ほぼ dd 精度の
行列積が `dgemm` 数十回分のコストで得られます。

線形方程式は、まず double に丸めた係数行列を `dgetrf` で LU 分解して
double 精度の近似解を求め、その後 `b - A x` の残差（尾崎スキームで
高精度に計算）を `dgetrs` で解いて解を修正する残差反復を、dd 精度に
収束するまで繰り返します。LU 分解は最初の 1 回だけで、反復中は前進後退
代入と残差計算のみを行います。

固有値問題は、まず `dsyev`（一般化問題は `dsygv`）で double 精度の
近似固有対を求め、荻田-相島の反復改良法
（T. Ogita, K. Aishima, "Iterative refinement for symmetric eigenvalue
decomposition", Japan Journal of Industrial and Applied Mathematics,
Vol. 35, pp. 1007-1035, 2018,
DOI: <https://doi.org/10.1007/s13160-018-0310-3>）で固有値と固有ベクトルを
同時に高精度化します。1 回の反復で正しい桁数がほぼ倍増する
（二次収束する）ため、標準では 2 回の反復で dd 精度に到達します。
反復内の行列積はすべて尾崎スキームで計算します。一般化問題では
B 内積（`transpose(X) * B * X = I` の正規化）を使う同じ形の反復を
行います。`eigsym`/`eigsymge` の `itep` 引数は反復回数で、dd 精度を
確保するため内部では最低 2 回反復します。

### 精度保証付き数値計算での利用

`vcp::pddblas` 自体は近似計算 policy ですが、`vcp::imats<T, P>` の
第 2 テンプレート引数（近似計算用 policy）に指定することで、
`kv::interval<kv::dd>` の精度保証付き数値計算に利用できます。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>

#include <vcp/pddblas.hpp>
#include <vcp/imats.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

typedef kv::interval<kv::dd> I;
vcp::matrix<I, vcp::imats<kv::dd, vcp::pddblas> > A, b, x;
```

`vcp::imats` で `P` が使われるのは近似計算（`lss` の近似逆行列、
`eigsym`/`eigsymge` の近似固有対など）だけで、包含の計算はすべて
`kv::interval<kv::dd>` の区間演算で行われるため、精度保証の正しさは
保たれます。また、`vcp::imats` は hardware の丸めモードを変更しない
ため、round-to-nearest を前提とする `kv::dd` 演算や尾崎スキームとも
干渉しません。

特に固有値問題では、荻田-相島反復により近似固有対が dd 精度まで
改良されるため、既定の `vcp::imats<kv::dd>` と比べて検証が高速になる
だけでなく、残差 `||A V - V Lambda||` が小さくなり固有値の包含幅も
小さくなります。

精度保証付き数値計算で使う場合の留意点:

- 初期近似を double で計算する制限（次節参照）はそのまま残ります。
  条件数が約 `1e16` を超える問題では近似逆行列が破綻し、
  `||R A - I|| < 1` などの検証条件が成立せず
  `vcp::verification_error` になることがあります。その場合は既定の
  `vcp::imats<kv::dd>` や `kv::mpfr<N>` 系の利用を検討してください。
- 包含の計算自体は `kv::interval<kv::dd>` の区間演算のままなので、
  行列サイズが大きい場合の全体の実行時間は区間演算部が支配的に
  なります。`pddblas` が高速化するのは近似計算の部分です。

### 注意点と制限

- `pddblas` は近似計算 policy です。結果は高精度ですが、
  精度保証（包含）は行いません。包含が必要な場合は区間演算の policy を
  使ってください。
- **初期近似を double で計算するため、double で解けない問題は解けません。**
  たとえば線形方程式では、係数行列の条件数が約 `1e16` を超えると
  double の LU 分解が数値的に破綻し、残差反復が収束せず、dd 精度どころか
  解そのものが得られない（または `vcp::lapack_error` が投げられる）ことが
  あります。残差反復で改善できるのは、おおよそ
  `条件数 × 1e-16 < 1` を満たす問題に限られます。条件数がそれを超える
  問題では、`vcp::mats<kv::dd>`（dd 精度のままの LU 分解）や
  `kv::mpfr<N>` などの高精度型の利用を検討してください。
- 同様に、固有値問題でも固有値の分離が double 精度で判別できない
  （極端にクラスター化した）場合、荻田-相島反復は該当する固有ベクトルを
  個別には分離できず、その部分の改良は制限されます。
- 得られた解の精度は、`b - A x` などの残差を `pddblas` の高精度な
  行列積で計算して確認することを推奨します。
- 行列要素の絶対値が極端に大きい（`1e300` 程度を超える）場合、
  尾崎スキームの分割で指数部が overflow するため正しく動作しません。
- `pdblas` と同様、LAPACK routine が非ゼロの `info` を返した場合は
  `vcp::lapack_error` を投げます。
- 利用プログラムは BLAS/LAPACK または MKL とリンクする必要があり、
  `kv::dd` を使うため kv ライブラリの include path も必要です。

速度の目安として、`vcp::mats<kv::dd>` との比較で、線形方程式は
n = 1000 で約 30 倍、対称固有値問題は n = 100 で約 1000 倍以上高速です。
行列積は `dgemm` 数十回分のコストのため `vcp::pdblas`（double）より
30〜40 倍程度遅くなりますが、`vcp::mats<kv::dd>` よりは高速で、
かつ丸め誤差は大幅に小さくなります。

## mats2<T> の詳細

`vcp::mats2<T>` は、tblas（テンプレート版 BLAS、`vcp/tblas/`）と
tlapack（テンプレート版 LAPACK、`vcp/tlapack/`）を下請けとして使う
近似計算 policy です。`vcp::mats<T>` を継承し、行列積・線形方程式・逆行列・
Cholesky 分解・対称固有値問題・一般化対称固有値問題をオーバーライドしています。
その他の演算（要素ごとの演算、生成、整形、集約など）は `mats<T>` から継承します。

```cpp
#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

vcp::matrix<double, vcp::mats2<double> > A, b, x;

A.rand(100);
b.rand(100, 1);
x = lss(A, b);
```

### インクルード構成と動作モード

`mats2<T>` の最大の特徴は、**どの特殊化ヘッダをインクルードするかによって
計算の実装が変わる**点です。`mats2.hpp` 自体は特殊化ヘッダをインクルード
しないため、利用者が明示的に選択します。

#### モード 1: 汎用実装（追加依存なし）

```cpp
#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

外部 BLAS/LAPACK は不要です。`vcp/tblas/tblas.hpp` と
`vcp/tlapack/tlapack.hpp` の汎用 C++ 実装が使われ、精度・速度とも
`vcp::mats<T>` と同等です。`double` だけでなく `kv::dd`、`kv::mpfr<N>` など
`tlapack` が対象とする任意の型 `T` で利用できます。

#### モード 2: double BLAS/LAPACK 特殊化（`T=double` のとき `pdblas` 相当）

```cpp
// tlapack_double.hpp は内部で tblas_double.hpp も include する
#include <vcp/tlapack/tlapack_double.hpp>

#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`T=double` のとき、`vcp::pdblas` と同等の動作になります。
外部 BLAS/LAPACK または Intel MKL とのリンクが必要です
（`-DUSE_VCP_BLAS -DUSE_VCP_LAPACK` 指定時は不要）。

#### モード 3: kv::dd 特殊化（`T=kv::dd` のとき `pddblas` 相当）

```cpp
#include <kv/dd.hpp>

// tlapack_dd.hpp は内部で tblas_double.hpp, tblas_dd.hpp,
// tlapack_double.hpp も include するため、これだけで十分
#include <vcp/tlapack/tlapack_dd.hpp>

#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

モード 2 の内容を含んだ上で、さらに `T=kv::dd` のとき
`vcp::pddblas` と同等の動作になります。外部 BLAS/LAPACK または Intel MKL、
および kv ライブラリが必要です。

`tlapack_dd.hpp` のインクルード順は `mats2.hpp` より前であれば任意ですが、
明示的特殊化は対象 template の暗黙実体化より前に参照される必要があるため、
**`mats2.hpp` より先にインクルードしてください**。

### 使用するアルゴリズム

`mats2<T>` が tblas/tlapack を通じて呼び出す演算は以下の通りです。

| 演算 | 汎用実装（モード 1） | double 特殊化（モード 2） | kv::dd 特殊化（モード 3） |
| --- | --- | --- | --- |
| 行列積 `A * B` | tblas 汎用 `tgemm<T>` | BLAS `dgemm_` | 尾崎スキーム（`tgemm<kv::dd>`）|
| `transpose(A)*A` | tblas 汎用 `tsyrk<T>` | BLAS `dsyrk_` | 尾崎スキーム（`tsyrk<kv::dd>`）|
| 線形方程式 `lss(A, b)` | tlapack 汎用 `tgesv<T>` | LAPACK `dgesv_` | double 初期解 + kv::dd 残差反復 |
| 逆行列 `inv(A)` | tlapack 汎用 `tgetrf/tgetri<T>` | LAPACK `dgetrf_/dgetri_` | double LU + kv::dd 残差反復 |
| Cholesky 分解 | tlapack 汎用 `tpotrf<T>` | LAPACK `dpotrf_` | double Cholesky + kv::dd 残差補正 |
| 対称固有値 `eigsym` | tlapack 汎用 `tsyev<T>` | LAPACK `dsyev_` | double 初期固有対 + kv::dd 残差反復 |
| 一般化対称固有値 `eigsymge` | tlapack 汎用 `tsygv<T>` | LAPACK `dsygv_` | double 初期固有対 + kv::dd 残差反復 |

### 精度保証付き数値計算での利用

`vcp::mats2<T>` 自体は近似計算 policy ですが、`vcp::imats<T, P>` の
第 2 テンプレート引数に指定することで、区間演算による精度保証付き数値計算に
利用できます。

モード 2（double 特殊化有効）の `mats2<double>` は `vcp::pidblas` と同等に
なります。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/tlapack/tlapack_double.hpp>

#include <vcp/imats.hpp>
#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

typedef kv::interval<double> I;
vcp::matrix<I, vcp::imats<double, vcp::mats2<double> > > A, b, x;
```

モード 3（kv::dd 特殊化有効）の `mats2<kv::dd>` は
`vcp::imats<kv::dd, vcp::pddblas>` と同等になります。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>

#include <vcp/tlapack/tlapack_dd.hpp>

#include <vcp/imats.hpp>
#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>

typedef kv::interval<kv::dd> I;
vcp::matrix<I, vcp::imats<kv::dd, vcp::mats2<kv::dd> > > A, b, x;
```

### pdblas・pddblas との比較

| | `vcp::pdblas` / `vcp::pddblas` | `vcp::mats2<T>` |
| --- | --- | --- |
| 内部実装 | `vcp/vblas/`・`vcp/vlapack/` | `vcp/tblas/`・`vcp/tlapack/` |
| 高速化の切り替え方法 | コンパイルフラグ（`-DUSE_VCP_BLAS` 等）または外部 BLAS/LAPACK リンク | インクルードする特殊化ヘッダで切り替え |
| 外部依存（デフォルト） | BLAS/LAPACK またはその代替 | なし（汎用モード） |
| 対象型 | `double` 固定（pdblas）/ `kv::dd` 固定（pddblas） | 任意の `T`（特殊化は `double` と `kv::dd`） |
| 動作等価性 | — | モード 2 で `pdblas`、モード 3 で `pddblas` と同等 |

### 注意点と制限

- `mats2<bool>` は利用できません（`static_assert` でコンパイルエラーになります）。
- 区間型（`kv::interval<T>`）を `T` として直接使うことはできません。
  区間演算には `vcp::imats<T, vcp::mats2<T> >` を使ってください。
- `tlapack` は `kv::interval<T>` を対象としないため、`T` は `double`、
  `kv::dd`、`kv::mpfr<N>` などの点の型に限定されます。
- モード 3（kv::dd 特殊化）の制限は `vcp::pddblas` と同様です。
  初期近似を double で計算するため、条件数が約 `1e16` を超える問題では
  残差反復が停滞することがあります。

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
