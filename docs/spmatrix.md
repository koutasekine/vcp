# Sparse Matrix / Policy Guide

この文書では、VCP Library の疎行列クラス `vcp::spmatrix` と policy の関係を
説明します。密行列 `vcp::matrix` については [matrix.md](matrix.md) を参照して
ください。

## 概要

### この文書で分かること

- `vcp::spmatrix` の公開されている関数の全部の呼び方と、戻り値の型、失敗の
  ときの振る舞い(例外を投げるか、status を返すか)
- オプションの構造体の全部の欄と既定値
- policy(`spmats`・`spumar`・`spcmodumar`・`spimats`)ごとに、同じ呼び方が
  何を計算するか、どのオプションが効くか

基本形は次の通りです。

```cpp
vcp::spmatrix<ElementType, Policy>
```

`ElementType` は行列要素の型です。`double`、`kv::dd`、`kv::mpfr<N>`、
`kv::interval<double>` などを指定します。`Policy` は、格納、線形方程式、
分解、固有値計算などをどの実装で行うかを指定します。`Policy` を省略すると
`vcp::spmats<ElementType>` が使われます。

### 最小の例

疎行列を組み立て、連立一次方程式を解き、小さい側の固有値を 2 個求める
プログラムです。

```cpp
#include <iostream>
#include <vector>

#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

int main() {
    const int n = 5;

    // 三重対角行列(対角 2、副対角 -1)を組み立てる
    vcp::spmatrix<double> A;
    A.resize(n, n);
    for (int i = 0; i < n; i++) {
        A.add(i, i, 2.0);
        if (i + 1 < n) {
            A.add(i, i + 1, -1.0);
            A.add(i + 1, i, -1.0);
        }
    }
    A.finalize();

    // A x = b を解く(失敗したときは例外)
    std::vector<double> b(n, 1.0);
    std::vector<double> x = A.solve(b);
    for (int i = 0; i < n; i++) std::cout << "x[" << i << "] = " << x[i] << "\n";

    // 小さい側の固有値を 2 個求める(収束したかどうかを converged で確かめる)
    vcp::eig_result<double> r = A.eigs_with_info(2);
    if (!r.converged) return 1;
    std::cout << "lambda = " << r.eigenvalues[0] << ", " << r.eigenvalues[1] << "\n";

    // 同じ計算を、密行列で受ける書き方で行う(E は 2 x 1)
    vcp::matrix<double, vcp::mats2<double> > E;
    eigs(A, E, 2);
    std::cout << "E = " << E(0, 0) << ", " << E(1, 0) << "\n";
    return 0;
}
```

コンパイルの例(`/path/to/libs` は `vcp/` を置いたディレクトリ):

```bash
g++ -std=c++20 -I/path/to/libs example.cpp
```

### 必要なヘッダとコンパイル

区間演算や `kv::dd`・`kv::mpfr<N>` を使う例では、VCP Library とは別に
取得した kv ライブラリが include path から見える必要があります。

通常の疎行列計算(`vcp::spmats<T>`。追加依存なし):

```cpp
#include <vcp/spmatrix.hpp>
```

密行列との相互変換や、密行列を右辺に取る `lss`、密行列で受ける `eigs` を
使う場合は、`matrix.hpp` も include します。

```cpp
#include <vcp/mats2.hpp>
#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>
```

区間行列(`vcp::spimats<T>`):

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/spimats.hpp>
#include <vcp/spmatrix.hpp>
```

kv の丸めヘッダ(`kv/rdouble.hpp`、`kv::dd` なら `kv/rdd.hpp`、
`kv::mpfr<N>` なら `kv/rmpfr.hpp`)は `spimats.hpp` より前に include して
ください。

UMFPACK・ARPACK に委譲する policy(`vcp::spumar`):

```cpp
#include <vcp/spumar.hpp>
#include <vcp/spmatrix.hpp>
```

UMFPACK・ARPACK・CHOLMOD に委譲する policy(`vcp::spcmodumar`):

```cpp
#include <vcp/spcmodumar.hpp>
#include <vcp/spmatrix.hpp>
```

外部ライブラリとのリンクが必要になるのは、`spumar.hpp` または
`spcmodumar.hpp` を include した翻訳単位だけです。VCP の他のヘッダは
これらを include しません。リンクのオプションは次の通りです。

| policy | リンクのオプション |
| --- | --- |
| `vcp::spmats<T>`・`vcp::spimats<T>` | なし |
| `vcp::spumar` | `-lumfpack -larpack` |
| `vcp::spcmodumar` | `-lumfpack -larpack -lcholmod` |

### policy の一覧と選び方

| Policy | 典型的な行列型 | 主な用途 |
| --- | --- | --- |
| `vcp::spmats<T>` | `vcp::spmatrix<T, vcp::spmats<T> >` | 汎用の疎行列計算(近似計算)。追加依存なし |
| `vcp::spimats<T>` | `vcp::spmatrix<kv::interval<T>, vcp::spimats<T> >` | 区間疎行列の精度保証付き計算 |
| `vcp::spumar` | `vcp::spmatrix<double, vcp::spumar>` | 疎 LU を UMFPACK、疎固有値を ARPACK に委譲する近似計算 |
| `vcp::spcmodumar` | `vcp::spmatrix<double, vcp::spcmodumar>` | `spumar` に加えて、疎 Cholesky 分解を CHOLMOD に委譲する近似計算 |

`vcp::spmatrix<T>` の default policy は `vcp::spmats<T>` です。

```cpp
vcp::spmatrix<double> A;
vcp::spmatrix<double, vcp::spmats<double> > B;
```

この 2 つは同じ policy を使います。

`vcp::spumar` と `vcp::spcmodumar` は `double` 専用です。他の要素型を
指定すると、最初の `add` / `set` でコンパイルエラーになります。

まず疎行列計算をする場合は `vcp::spmats<T>` を使います。追加依存がなく、
`double` だけでなく `kv::dd`、`kv::mpfr<N>` などでも同じインターフェースを
使えます。区間演算による包含が必要な場合は `vcp::spimats<T>` を使います。

```cpp
typedef kv::interval<double> I;
vcp::spmatrix<I, vcp::spimats<double> > A;
```

`double` の大規模な問題で、fill-in の多い直接法や固有値計算が遅い場合は
`vcp::spumar` を使います。対称正定値行列の Cholesky 分解も外部に委譲したい
場合は `vcp::spcmodumar` を使います。どちらも外部ライブラリのライセンスに
注意してください(「外部ライブラリのライセンスに関する注意」の節)。

```cpp
vcp::spmatrix<double, vcp::spumar> A;
vcp::spmatrix<double, vcp::spcmodumar> B;
```

policy ごとの違いは「早見表」の節にまとめ、機能ごとの節の「policy ごとの違い」に
詳しく書いています。各 policy の内部で使われる下請けは [kernels/sparse_layers.md](kernels/sparse_layers.md) を参照してください。

### 目次

- [早見表](#早見表)
- [行列の基本](#行列の基本)
- [演算](#演算)
- [密行列・他の型との変換](#密行列他の型との変換)
- [線形方程式](#線形方程式)
- [LU 分解](#lu-分解)
- [LDLᵀ 分解と慣性](#ldlᵀ-分解と慣性)
- [Cholesky 分解](#cholesky-分解)
- [近似逆行列(AINV・FSAI)](#近似逆行列ainvfsai)
- [固有値問題](#固有値問題)
- [`spimats` の詳細](#spimats-の詳細)
- [`spumar` の詳細](#spumar-の詳細)
- [`spcmodumar` の詳細](#spcmodumar-の詳細)
- [外部ライブラリのライセンスに関する注意](#外部ライブラリのライセンスに関する注意)

「行列の基本」から「固有値問題」までの節は、どれも次の 4 つの小節から
なります。

- 版の一覧: 全部の呼び方、戻り値の型、失敗のときの振る舞い
- オプション: オプションの構造体の全部の欄と既定値
- 戻り値と status(または例外): 結果の読み方
- policy ごとの違い: 4 つの policy で何が変わるか

## 早見表

同じ呼び方が、policy によって何を計算するかの一覧です。セルの言葉の意味は
次の通りです。

- VCP の実装: `spmats` の実装がそのまま動きます(近似計算)。
- UMFPACK に委譲 / ARPACK に委譲 / CHOLMOD に委譲: 外部ライブラリを呼びます
  (近似計算)。
- 保証つき: 区間演算による包含、または証明を返します。
- 区間演算で実行(保証の主張なし): 区間を要素にした `spmats` の経路で、
  非推奨です。計算は区間演算で進みますが、結果が何かを保証するという主張は
  ありません。
- 使えない: コンパイルエラーになります。

| 機能 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| 行列の組み立て・格納・部分行列 | VCP の実装 | VCP の実装(`double` だけ) | VCP の実装(`double` だけ) | VCP の実装(要素は区間) |
| 四則・ベクトルとの積・構造・要素ごとの関数・集約 | VCP の実装 | VCP の実装 | VCP の実装 | 保証つき(区間演算による包含) |
| `normone` / `norminf` | VCP の実装 | VCP の実装 | VCP の実装 | 保証つき(上端が上界) |
| `normtwo` | VCP の実装 | VCP の実装 | VCP の実装 | 近似計算(保証の主張なし) |
| `is_symmetric()` | VCP の実装 | VCP の実装 | VCP の実装 | VCP の実装(区間の端点の一致による厳密な判定) |
| `is_symmetric(tol)` | VCP の実装 | VCP の実装 | VCP の実装 | 使えない |
| 密行列・他の型との変換 | VCP の実装 | VCP の実装 | VCP の実装 | VCP の実装(点から区間への変換は包含を保つ) |
| `solve` / `solve_with_info` / `lss` / `lss_with_info` | VCP の実装(解法は `method`。既定は疎 LU) | UMFPACK に委譲 | UMFPACK に委譲 | 保証つき(解の包含) |
| `solve_cg` などの短縮形 | VCP の実装(その解法) | UMFPACK に委譲 | UMFPACK に委譲 | 保証つき(解の包含) |
| `lu_factorize_with_info`(ハンドル) | VCP の実装 | VCP の実装 | VCP の実装 | 区間演算で実行(保証の主張なし) |
| `lu` / `lu_with_info` | VCP の実装 | UMFPACK に委譲 | UMFPACK に委譲 | 区間演算で実行(保証の主張なし) |
| `lu_solve` / `lu_inverse_row` | VCP の実装 | VCP の実装 | VCP の実装 | 区間演算で実行(保証の主張なし) |
| `ldl` / `ldl_with_info` / `inertia` / `ldl_shift_setup_with_info` | VCP の実装 | VCP の実装 | VCP の実装 | 区間演算で実行(保証の主張なし) |
| `chol` / `chol_with_info` / `Cholesky` | VCP の実装 | VCP の実装 | CHOLMOD に委譲 | 区間演算で実行(保証の主張なし) |
| `is_spd` | VCP の実装(近似の判定) | VCP の実装(近似の判定) | CHOLMOD に委譲(近似の判定) | 保証つき(`true` は証明) |
| `ainv_*` / `fsai_*` / `fsai_adaptive_*` | VCP の実装 | VCP の実装 | VCP の実装 | 区間演算で実行(保証の主張なし) |
| `eigs` / `eigs_with_info`(標準) | VCP の実装(解法は `method`) | ARPACK に委譲(一部は VCP の実装) | ARPACK に委譲(一部は VCP の実装) | 保証つき(対称行列の固有値の包含) |
| `eigs` / `eigs_with_info`(一般化) | VCP の実装(既定のオプションでは対角行列どうしだけ) | ARPACK に委譲(一部は VCP の実装) | ARPACK に委譲(一部は VCP の実装) | 保証つき(対称行列の固有値の包含) |

`spumar`・`spcmodumar` で ARPACK に委譲されない指定は「固有値問題」の節の
「policy ごとの違い」に書いています。

## 行列の基本

行列要素は 0 始まりの添字で指定します。要素は `add` / `set` で追加し、
組み立てが終わったら `finalize()` を呼びます。

```cpp
vcp::spmatrix<double> A;

A.resize(3, 3);
A.add(0, 0, 2.0);
A.add(1, 1, 3.0);
A.add(1, 1, 1.0);   // add は加算: A(1,1) = 4
A.set(2, 2, 5.0);   // set は上書き
A.finalize();
```

- `add(i, j, v)` は同じ位置の値に加算し、`set(i, j, v)` はそれまでの値を
  上書きします。有限要素法の組み立てでは `add` を使います。
- 組み立て中は COO 形式で保持し、`finalize()` で CSR 形式に変換します。
  `to_csr()` / `to_csc()` で格納形式を指定することもできます。
- `finalize()` 後に `add` / `set` を呼ぶと、COO 形式に戻ります。
- 値が厳密に零の要素は格納しません。足した結果が零になった要素や、零で
  上書きした要素も、`finalize()` のときに取り除かれます。

### 版の一覧

`T` は要素型、`index_type` は添字の型です。「失敗のとき」の欄に例外の型を
書いていないものは、失敗しません。

構成・大きさ・状態:

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `vcp::spmatrix<T, P> A;` | --- | 0×0 の行列 | --- |
| `vcp::spmatrix<T, P> A(rows, cols);` | --- | `resize(rows, cols)` と同じ | `resize` と同じ |
| `A.resize(rows, cols)` | `void` | 大きさを決め、格納要素を捨てる(COO 形式、未 finalize) | 負の大きさは `vcp::invalid_argument` |
| `A.clear()` | `void` | 0×0 に戻す | --- |
| `A.reserve(n)` | `void` | 格納数の予約 | 負の数は `vcp::invalid_argument` |
| `A.eye(n)` | `void` | 中身を n×n の単位行列に置き換える(finalize 済み) | --- |
| `A.zeros(r, c)` | `void` | r×c で格納要素なし(finalize 済み) | --- |
| `A.rowsize()` | `index_type` | 行数 | --- |
| `A.columnsize()` | `index_type` | 列数 | --- |
| `A.nnz()` | `index_type` | 零でない要素の数(未 finalize では、同じ位置に重ねて入れた分を別々に数える) | --- |
| `A.stored_nnz()` | `index_type` | 内部に持っている要素の数(未 finalize では、重ねて入れた分と、値が零のものも数える) | --- |
| `A.is_finalized()` | `bool` | CSR または CSC に変換済みか | --- |
| `A.is_sorted()` | `bool` | 格納要素が添字の順に並んでいるか | --- |
| `A.is_unique()` | `bool` | 同じ位置の要素が重複していないか | --- |
| `A.format()` | `format_type` | `vcp::sparse_coo` / `vcp::sparse_csr` / `vcp::sparse_csc` | --- |

要素の読み書き:

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `A.add(i, j, v)` | `void` | (i, j) に v を足す | 添字の範囲外は `vcp::index_error` |
| `A.set(i, j, v)` | `void` | (i, j) を v で上書きする | 同上 |
| `A.get(i, j)` | `T` | (i, j) の値。格納していない位置は 0 | 同上 |
| `T c = A(i, j)` | `T` | `get(i, j)` と同じ | 同上 |
| `A(i, j) = v` | --- | `set(i, j, v)` と同じ | 同上 |
| `A(i, j) += v` | --- | `add(i, j, v)` と同じ | 同上 |
| `A(i, j) -= v` | --- | `add(i, j, -v)` と同じ | 同上 |
| `A(i, j) *= v` | --- | 今の値に v を掛けて `set` | 同上 |
| `A(i, j) /= v` | --- | 今の値を v で割って `set` | 同上 |

`get` と `A(i, j)` の読み出しは、未 finalize でも使えます(finalize は
行いません)。`A(i, j)` への書き込みは `add` / `set` と同じなので、finalize
済みの行列は COO 形式に戻ります。書き込みのあとは `finalize()` を呼んで
ください。

格納形式:

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `A.finalize()` | `void` | CSR 形式に変換する | --- |
| `A.to_csr()` | `void` | CSR 形式に変換する | --- |
| `A.to_csc()` | `void` | CSC 形式に変換する | --- |
| `A.as_csr()` | `spmatrix` | CSR 形式の複製を返す(元の行列は変えない) | --- |
| `A.as_csc()` | `spmatrix` | CSC 形式の複製を返す(元の行列は変えない) | --- |
| `A.sort_coo()` | `void` | COO 形式のまま、要素を添字の順に並べる | --- |
| `A.normalize_coo()` | `void` | COO 形式のまま、重複をまとめ、零を取り除く | --- |
| `A.assign_csr(rows, cols, row_ptr, col_ind, value)` | `void` | CSR の配列から構成する(finalize 済み) | 並びや重複・零の違反は `vcp::invalid_argument`、添字の範囲外は `vcp::index_error` |
| `A.assign_csc(rows, cols, col_ptr, row_ind, value)` | `void` | CSC の配列から構成する(finalize 済み) | 同上 |
| `A.outer_index()` | `const std::vector<index_type>&` | CSR の行ポインタ(CSC では列ポインタ) | --- |
| `A.inner_index()` | `const std::vector<index_type>&` | CSR の列添字(CSC では行添字) | --- |
| `A.values()` | `const std::vector<T>&` | 格納要素の値 | --- |
| `A.coo_rows()` | `const std::vector<index_type>&` | COO の行添字(finalize 済みでは空) | --- |
| `A.coo_columns()` | `const std::vector<index_type>&` | COO の列添字(同上) | --- |
| `A.coo_values()` | `const std::vector<T>&` | COO の値(同上) | --- |

`sort_coo()` と `normalize_coo()` は `finalize()` の中でも行われるので、
ふつうは呼ぶ必要がありません。finalize 済みの行列に呼ぶと、COO 形式(未
finalize)に戻ります。`assign_csr` / `assign_csc` に渡す配列は、
各行(列)の添字が昇順・重複なしで、厳密零の値を含まない必要があります。
`outer_index()`・`inner_index()`・`values()` は finalize 済みの行列で
読んでください。

部分行列:

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `A.submatrix({...}, {...})` | `spmatrix` | 行と列を選んで取り出す | 選び方が不正なら `vcp::index_error` |
| `A({...}, {...})` | `spmatrix` | `submatrix` と同じ | 同上 |
| `A(i, {...})` | `spmatrix` | 行を 1 つ、列を選んで取り出す(1×列数) | 同上 |
| `A({...}, j)` | `spmatrix` | 行を選び、列を 1 つ取り出す(行数×1) | 同上 |
| `A({...}, {...}) = B` | --- | ブロックを疎行列 B で置き換える | 選び方が不正なら `vcp::index_error`、大きさが合わなければ `vcp::dimension_error` |
| `A({...}, {...}) = s` | --- | ブロックの全部の位置をスカラー s にする | 選び方が不正なら `vcp::index_error` |
| `A(i, {...}) = B` | --- | 行 i の選んだ列を B で置き換える | `A({...}, {...}) = B` と同じ |
| `A({...}, j) = B` | --- | 列 j の選んだ行を B で置き換える | 同上 |

行・列の選び方(`{...}` の中)は次の 4 通りです。添字は 0 始まりです。

| 書き方 | 意味 |
| --- | --- |
| `{}` | その軸の全部 |
| `{i}` | 添字 i だけ |
| `{a, b}` | a から b まで(b を含む) |
| `{a, s, b}` | a から s 刻みで、b を超えないところまで |

```cpp
vcp::spmatrix<double> S = A({0, 1}, {1, 2});   // 行 0〜1、列 1〜2 の 2×2
A({0, 1}, {1, 2}) = S;
```

- 取り出した結果は finalize 済みの新しい行列で、元の行列は変わりません。
  格納形式は元と同じです(元が CSC なら CSC、ほかは CSR)。
- 波括弧の付いた指定が 1 つでもあれば、結果は行列です。`A({1}, {2})` は
  1×1 の行列、`A(1, 2)` は要素です。
- `A({...}, {...}) = B` は置き換えです。B が零の位置にあった要素は消えます。
  書き込みのあと、行列は finalize 済みに戻ります。
- `A({...}, {...}) = s` に零でない s を渡すと、ブロックの全部の位置に要素が
  入ります(ブロックが密になります)。s が零なら、ブロックの要素が消えます。
- 非 const の行列で `auto x = A({0, 1}, {1, 2});` と書くと、x は部分行列の
  複製ではなく、A を指す書き込み用の代理になります。複製がほしいときは、
  上の例のように `vcp::spmatrix<double>` の変数で受けてください。

そのほか:

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `A.is_symmetric()` | `bool` | 対称か(許容誤差は 1e-12)。正方でなければ `false` | --- |
| `A.is_symmetric(tol)` | `bool` | 対称か(許容誤差を指定) | `tol` が正でなければ `vcp::invalid_argument` |
| `std::cout << A` | `std::ostream&` | 行列を出力する | --- |

`std::cout << A` は、1 行目に `spmatrix 行数 x 列数, nnz = 個数` を出力し、
2 行目からは `(i,j) 値` の形で格納要素を 1 行に 1 つずつ出力します。

```text
spmatrix 2 x 3, nnz = 2
(0,0) -1
(1,2) 2.5
```

`vcp::spmatrix_builder` は、組み立てと確定を分けて書きたいときの補助クラス
です。policy は `vcp::spmats<T, Index>` に固定されます。

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `vcp::spmatrix_builder<T, Index> b(rows, cols);` | --- | rows×cols の組み立て器(`Index` を省略すると `int`) | --- |
| `b.reserve(n)` | `void` | 格納数の予約 | 添字の型で表せない数は `vcp::invalid_argument` |
| `b.add(i, j, v)` | `void` | (i, j) に v を足す | `A.add` と同じ |
| `b.set(i, j, v)` | `void` | (i, j) を v で上書きする | `A.set` と同じ |
| `b.finalize_csr()` | `vcp::spmatrix<T, vcp::spmats<T, Index> >` | CSR 形式の行列を返す | --- |
| `b.finalize_csc()` | `vcp::spmatrix<T, vcp::spmats<T, Index> >` | CSC 形式の行列を返す | --- |

`vcp::spmatrix` は、型の別名を持っています。添字の型は `index_type` です
(既定は `int`)。変更する場合は `vcp::spmats<T, long long>` のように policy の
第 2 テンプレート引数で指定します。添字型は符号付き整数である必要があります。

| 別名 | 内容 |
| --- | --- |
| `value_type` / `policy_type` | 要素型 / policy |
| `index_type` | 添字の型 |
| `format_type` | 格納形式の型 |
| `scalar_real_type` | 許容誤差などに使う実数の型 |
| `dense_matrix_type` | `to_dense()` の戻り値の型 |
| `linear_solve_options_type` / `sparse_lu_options_type` / `lu_extract_options_type` / `ldl_options_type` / `inertia_options_type` / `chol_options_type` / `spd_options_type` / `ainv_options_type` / `fsai_options_type` / `fsai_adaptive_options_type` / `eig_options_type` | 各オプションの型 |
| `lu_extract_result_type` / `ldl_result_type` / `inertia_result_type` / `chol_result_type` / `ainv_result_type` / `fsai_result_type` / `fsai_adaptive_result_type` | 各結果の型 |
| `lu_factor_handle_type` / `ldl_shift_handle_type` / `ldl_shift_workspace_type` / `ldl_shift_inertia_result_type` | ハンドルとその作業領域・結果の型 |
| `element_proxy_type` | `A(i, j)` が返す書き込み用の代理の型 |

### オプション

オプションの構造体はありません。

### 戻り値と例外

- 添字の範囲外は `vcp::index_error`、不正な引数は `vcp::invalid_argument`、
  大きさの不一致は `vcp::dimension_error` です。例外の型の一覧は
  [error_handling.md](error_handling.md) を参照してください。
- `add` / `set` / `A(i, j)` への書き込みのあと、行列は未 finalize です。
  この文書の計算の関数(演算、`solve`、分解、`eigs` など)は、未 finalize の
  行列を受け取ると自分で finalize してから計算します。ただし、
  `mul_vec` / `trans_mul_vec` は finalize 済みの行列でないと
  `vcp::state_error` を投げます。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| 要素型 | 任意 | `double` だけ(ほかの型は `add` / `set` でコンパイルエラー) | `spumar` と同じ | `kv::interval<T>` |
| 組み立て・格納・部分行列 | VCP の実装 | 同じ | 同じ | 同じ(零の判定は、両端が 0 の区間だけを零とみなす) |
| `is_symmetric()` | 許容誤差 1e-12 で判定 | 同じ | 同じ | 厳密に判定(すべての (i, j) で、A(i, j) と A(j, i) の区間の両端が一致するとき `true`) |
| `is_symmetric(tol)` | 許容誤差 `tol` で判定 | 同じ | 同じ | 使えない(コンパイルエラー) |

`spimats` の `is_symmetric()` は、`is_spd` や `eigs` の入口の検査と同じ判定です
(入口の検査は、対称でなければ `vcp::domain_error` を投げます)。許容誤差つきの
`is_symmetric(tol)` は、区間の行列では意味を持たないので使えません。

内部の格納形式と演算のアルゴリズムは [kernels/tsparse.md](kernels/tsparse.md) を参照してください。

## 演算

### 版の一覧

`A`・`B`・`C` は同じ要素型・同じ policy の疎行列、`D` は密行列
`vcp::matrix`、`x` は `std::vector<T>`、`s` はスカラーです。

| 種類 | 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- | --- |
| 四則 | `A + B` | `spmatrix` | 大きさが合わなければ `vcp::dimension_error` |
| 四則 | `A - B` | `spmatrix` | 同上 |
| 四則 | `A * B` | `spmatrix` | 同上 |
| 四則 | `-A` | `spmatrix` | --- |
| 四則 | `s * A` | `spmatrix` | --- |
| 四則 | `A * s` | `spmatrix` | --- |
| 四則 | `A / s` | `spmatrix` | --- |
| 複合代入 | `A += B` | `spmatrix&` | 大きさが合わなければ `vcp::dimension_error` |
| 複合代入 | `A -= B` | `spmatrix&` | 同上 |
| 複合代入 | `A *= B` | `spmatrix&` | 同上 |
| 複合代入 | `A *= s` | `spmatrix&` | --- |
| 複合代入 | `A /= s` | `spmatrix&` | --- |
| 疎行列とベクトル | `A * x` | `std::vector<T>` | 長さが合わなければ `vcp::dimension_error` |
| 疎行列とベクトル | `x * A` | `std::vector<T>` | 同上 |
| 疎行列とベクトル | `A.mul_vec(x)` | `std::vector<T>` | 同上。未 finalize なら `vcp::state_error` |
| 疎行列とベクトル | `A.mul_vec(xp, yp)` | `void` | 未 finalize なら `vcp::state_error`(長さは検査しない) |
| 疎行列とベクトル | `A.trans_mul_vec(x)` | `std::vector<T>` | `mul_vec(x)` と同じ |
| 疎行列とベクトル | `A.trans_mul_vec(xp, yp)` | `void` | `mul_vec(xp, yp)` と同じ |
| 疎行列と密行列 | `A * D` | `vcp::matrix` | 大きさが合わなければ `vcp::dimension_error` |
| 疎行列と密行列 | `D * A` | `vcp::matrix` | 同上 |
| 構造 | `transpose(A)` | `spmatrix` | --- |
| 構造 | `A.transpose()` | `spmatrix` | --- |
| 構造 | `diag(A)` | `spmatrix` | --- |
| 構造 | `vcp::diag(v)` | `spmatrix` | --- |
| 構造 | `tril(A)` | `spmatrix` | --- |
| 構造 | `triu(A)` | `spmatrix` | --- |
| 構造 | `A.horzcat(B, C)` | `void` | 行数が合わなければ `vcp::dimension_error` |
| 構造 | `horzcat(A, B, ...)` | `spmatrix` | 同上 |
| 構造 | `A.vercat(B, C)` | `void` | 列数が合わなければ `vcp::dimension_error` |
| 構造 | `vercat(A, B, ...)` | `spmatrix` | 同上 |
| 要素ごとの関数 | `abs(A)` | `spmatrix` | --- |
| 要素ごとの関数 | `sqrt(A)` | `spmatrix` | 負の要素があれば `vcp::invalid_argument` |
| 要素ごとの関数 | `sin(A)` | `spmatrix` | --- |
| 要素ごとの関数 | `cos(A)` | `spmatrix` | --- |
| 要素ごとの関数 | `exp(A)` | `spmatrix` | --- |
| 要素ごとの関数 | `log(A)` | `spmatrix` | 負の要素があれば `vcp::invalid_argument` |
| 要素ごとの関数 | `pow(A, B)` | `spmatrix` | 大きさが合わなければ `vcp::dimension_error` |
| 要素ごとの関数 | `pow(A, n)` | `spmatrix` | --- |
| 要素ごとの関数 | `pow(A, s)` | `spmatrix` | --- |
| 要素ごとの関数 | `pow(s, A)` | `spmatrix` | --- |
| 集約 | `max(A)` | `spmatrix` | --- |
| 集約 | `min(A)` | `spmatrix` | --- |
| 集約 | `sum(A)` | `spmatrix` | --- |
| 集約 | `length(A)` | `int` | --- |
| ノルム | `normone(A)` | `spmatrix`(1×1) | --- |
| ノルム | `norminf(A)` | `spmatrix`(1×1) | --- |
| ノルム | `normtwo(A)` | `spmatrix`(1×1) | --- |
| その他 | `ltransmul(A)` | `spmatrix` | --- |

- `A.mul_vec(x)` は `A x`、`A.trans_mul_vec(x)` は `Aᵀ x` です。`xp`・`yp` は
  配列へのポインタ(`const T*`・`T*`)で、結果を `yp` の指す配列に書きます。
- `A.horzcat(B, C)` は `C = [A, B]`、`A.vercat(B, C)` は `C = [A; B]` です。
  自由関数の `horzcat(A, B, ...)` / `vercat(A, B, ...)` は、1 個以上の行列を
  いくつでも並べて、結果を戻り値で返します。
- `diag(A)` は、A の形で意味が変わります。1×1 ならそのまま、n×1 または
  1×n なら n×n の対角行列、それ以外の r×c なら対角成分を並べた
  min(r, c)×1 の行列を返します。
- `vcp::diag(v)` は `std::vector<T>` から n×n の対角行列を作ります。policy を
  選ぶときは `vcp::diag<T, P>(v)` と書きます。
- `pow(A, n)` の n は `int` です。`pow(A, B)` は要素ごとの冪です。
- `length(A)` は行数と列数の大きいほうです。
- `ltransmul(A)` は `Aᵀ A` です。
- `A + B` などの演算子は、右辺値(一時オブジェクト)を渡すと、その記憶領域を
  再利用します。結果は同じです。
- 呼び出し方と結果の形は `vcp::matrix` の同名の関数に合わせています。

### オプション

オプションの構造体はありません。

### 戻り値と例外

- 演算の結果は finalize 済みの行列です。結果が厳密に零になった要素は格納
  されません(`A - A` は格納要素なし)。
- `max`・`min`・`sum` は、行列には列ごとの値を 1×列数 の行列で返し、
  ベクトル(n×1 または 1×n)には 1×1 の行列を返します。格納していない
  位置は 0 として数えます。
- `normone`・`norminf`・`normtwo` は 1×1 の行列を返します。値は
  `normone(A).get(0, 0)` のように取り出します。行列には作用素ノルム(最大列和・
  最大行和・最大特異値)、ベクトルにはベクトルのノルムを返します。
- `cos`、`exp`、`log`、`pow(s, A)` のように `f(0) ≠ 0` となる関数は、
  格納していない位置にも適用されるので、結果は構造的にほぼ密になります。
  `log(A)` では、格納していない位置が `-inf` になります。`pow(A, 0)` も
  全部の位置が 1 になります。
- `A / s` に零を渡しても例外は投げず、要素が `inf` になります。

`std::cout << A` の出力の形は「行列の基本」の節に書いています。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| 四則・ベクトルとの積・密行列との積・構造・要素ごとの関数・集約 | VCP の実装 | 同じ(結果の policy も `spumar`) | 同じ(結果の policy も `spcmodumar`) | 区間演算で計算し、結果は真の値を含む |
| `normone` / `norminf` | VCP の実装 | 同じ | 同じ | 区間で返り、上端が保証付きの上界 |
| `normtwo` | VCP の実装(近似) | 同じ | 同じ | 近似計算。幅のない区間が返り、保証には使えない |

内部の格納形式と演算のアルゴリズムは [kernels/tsparse.md](kernels/tsparse.md) を参照してください。

## 密行列・他の型との変換

### 版の一覧

| 呼び方 | 戻り値の型 | 内容 | 失敗のとき |
| --- | --- | --- | --- |
| `A.to_dense()` | `dense_matrix_type` | 疎行列を密な 2 次元配列にする | --- |
| `S = D` | `spmatrix&` | 密行列 `vcp::matrix` を疎行列に代入する | --- |
| `vcp::spmatrix<T2, P2> B = A;` | --- | 要素型または policy の違う疎行列から初期化する | --- |
| `B = A` | `spmatrix&` | 要素型または policy の違う疎行列を代入する | --- |

疎行列から密行列へは `to_dense()` で変換します。戻り値の型は
`spmatrix<T>::dense_matrix_type`(4 つの policy とも
`std::vector<std::vector<T> >`)です。

密行列から疎行列へは代入で変換します。厳密に零の成分は格納されません。
密行列から疎行列を初期化する書き方(`vcp::spmatrix<double> S = D;`)は
ありません。

```cpp
vcp::matrix<double, vcp::mats2<double> > D;
vcp::spmatrix<double> S;
S = D;
```

疎行列の間では、要素型や policy が異なる場合も代入・初期化で変換できます。
要素の変換規則は `vcp::matrix` と同じです([matrix.md](matrix.md))。

```cpp
vcp::spmatrix<double> A;
vcp::spmatrix<kv::dd> B = A;                            // 要素型の変換
vcp::spmatrix<double, vcp::spumar> C = A;               // policy の変換
vcp::spmatrix<kv::interval<double>, vcp::spimats<double> > IA = A;   // 点区間化
```

### オプション

オプションの構造体はありません。

### 戻り値と例外

- 変換の結果は finalize 済みの行列です。元の行列は変わりません(未 finalize の
  行列から変換しても、元は未 finalize のままです)。
- 疎行列の間の変換は、格納形式を保ちます(元が CSC なら CSC、ほかは CSR)。
  密行列からの代入の結果は CSR です。
- 変換した結果が厳密に零になった要素は格納されません。
- 点から区間への変換は、元の値を含む区間になります(`double` から
  `kv::interval<double>` へは幅のない区間、`kv::dd` から
  `kv::interval<double>` へは元の値を含む幅のある区間)。
- 区間から点への変換は、区間の中点を取ります。中点が零になった要素
  (たとえば `[-1, 1]`)は格納されません。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `to_dense()` の戻り値の型 | `std::vector<std::vector<T> >` | `std::vector<std::vector<double> >` | `spumar` と同じ | `std::vector<std::vector<kv::interval<T> > >` |
| 変換できる要素型 | 任意 | `double` だけ | `double` だけ | `kv::interval<T>` |
| policy の固有の設定 | --- | 変換の結果は、既定の設定を持つ(「`spumar` の詳細」の節) | `spumar` と同じ | 変換の結果は、既定の設定を持つ(「`spimats` の詳細」の節) |

## 線形方程式

```cpp
std::vector<double> x = A.solve(b);
auto r = A.solve_with_info(b);   // r.x, r.converged, r.residual_norm など

lss(A, b, x);                    // x = A \ b
vcp::matrix<double, vcp::mats2<double> > X = lss(A, B);   // 密な右辺 B
```

`solve` は失敗時に例外を投げます。`solve_with_info` は、解けなかったことを
例外ではなく結果の `converged` で返します(引数の誤りは `solve_with_info` でも
例外です。`spimats` では、保証に失敗したときも例外です)。

### 版の一覧

`b`・`x` は `std::vector<T>`、`B`・`X` は密行列 `vcp::matrix`、`Bs` は疎行列、
`opt` は `vcp::linear_solve_options<T>` です。

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.solve(b)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve(b, opt)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve_with_info(b)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `A.solve_with_info(b, opt)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `x = lss(A, b)` | `std::vector<T>` | `A.solve(b)` と同じ |
| `x = lss(A, b, opt)` | `std::vector<T>` | `A.solve(b, opt)` と同じ |
| `lss(A, b, x)` | `void` | `A.solve(b)` と同じ |
| `X = lss(A, B)` | `vcp::matrix` | 例外(`vcp::state_error`) |
| `lss(A, B, X)` | `void` | 例外(`vcp::state_error`) |
| `lss(A, Bs, X)` | `void` | 例外(`vcp::state_error`) |
| `lss_with_info(A, b)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `lss_with_info(A, b, opt)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `lss_with_info(A, Bs)` | `std::vector<vcp::linear_solve_result<T> >` | 列ごとの `converged` |
| `lss_with_info(A, Bs, opt)` | `std::vector<vcp::linear_solve_result<T> >` | 列ごとの `converged` |
| `lss_with_info(A, B, res)` | `void` | 列ごとの `converged` |
| `lss_with_info(A, B, res, opt)` | `void` | 列ごとの `converged` |

- 密な右辺・疎な右辺を取る `lss`(戻り値または `X` で受ける 3 つの版)は、
  オプションを取りません。既定のオプションで解きます。解は密行列で返ります。
  どれかの列が解けなければ、全部の列を計算したあとで例外を投げます。
- 複数の右辺を取る `lss_with_info` は、列ごとの結果を返します。`res` は
  `std::vector<vcp::linear_solve_result<T> >` です。

反復法の短縮形(`max_iter` は反復回数の上限、`tol` は相対残差の閾値):

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.solve_jacobi(b, max_iter, tol)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve_gauss_seidel(b, max_iter, tol)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve_cg(b, max_iter, tol)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve_bicgstab(b, max_iter, tol)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve_gmres(b, max_iter, tol)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.solve_jacobi_with_info(b, max_iter, tol, use_relative_residual)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `A.solve_gauss_seidel_with_info(b, max_iter, tol, use_relative_residual)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `A.solve_cg_with_info(b, max_iter, tol, check_symmetric, preconditioner, use_relative_residual)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `A.solve_bicgstab_with_info(b, max_iter, tol, use_relative_residual, preconditioner)` | `vcp::linear_solve_result<T>` | `converged = false` |
| `A.solve_gmres_with_info(b, max_iter, tol, restart, use_relative_residual, preconditioner)` | `vcp::linear_solve_result<T>` | `converged = false` |

- 短縮形は、`linear_solve_options` の対応する欄を埋めて `solve` /
  `solve_with_info` を呼ぶのと同じです。`solve_gmres` の `restart` は 30 です。
- `_with_info` の短縮形の `use_relative_residual`(既定 `true`)、
  `check_symmetric`(既定 `true`)、`preconditioner`(既定
  `vcp::preconditioner_type::none`)は省略できます。`solve_gmres_with_info` の
  `restart` は省略できません。

同じ行列で何度も解く場合は、分解を保持するハンドルを使います。

```cpp
auto h = A.lu_factorize_with_info();
if (h.valid()) {
    auto r1 = h.solve_with_info(b1);
    auto r2 = h.solve_with_info(b2);
}
```

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.lu_factorize_with_info()` | `lu_factor_handle_type` | ハンドルの `valid()` が `false` |
| `A.lu_factorize_with_info(opt)` | `lu_factor_handle_type` | 同上(`opt` は `vcp::sparse_lu_options<T>`) |
| `h.valid()` | `bool` | --- |
| `h.status()` | `vcp::sparse_lu_status` | --- |
| `h.info()` | `const vcp::sparse_lu_info<T, index_type>&` | --- |
| `h.solve_with_info(b)` | `vcp::linear_solve_result<T>` | `converged = false`(解は零ベクトル) |
| `h.solve(b)` | `std::vector<T>` | 例外(`vcp::state_error`) |

ハンドルは、分解のときのオプション(反復改良の有無など)を覚えていて、
解くときには変えられません。`h.info()` は分解の記録で、主な欄は
`success`・`status`・`method_used`(実際に使った分解法)・`n`・`nnz_L`・
`nnz_U` です。

### オプション

`vcp::linear_solve_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 解法(下の表) |
| `max_iter` | 1000 | 反復法の反復回数の上限。0 は例外(直接法でも検査されます) |
| `tol` | 1e-12 | 反復法の残差の閾値。0 以下は例外(直接法でも検査されます) |
| `check_symmetric` | `true` | 共役勾配法で、行列が対称かを検査する |
| `use_relative_residual` | `true` | 反復法の判定に相対残差を使う(`false` なら絶対残差) |
| `restart` | 30 | GMRES の再開周期。0 は例外 |
| `preconditioner` | `none` | 前処理。`none` / `jacobi`。共役勾配法・BiCGSTAB 法・GMRES 法で使う |
| `sparse_lu` | 既定の構築 | 疎 LU で解くときのオプション(`sparse_lu_options<T>`。下の表) |

解法は `method`(`vcp::linear_solver_method`)で選びます。既定の `auto_select`
は疎 LU 直接法(`sparse_lu`)に解決されます。

| `method` | 解法 |
| --- | --- |
| `auto_select`(既定) | `sparse_lu` |
| `sparse_lu` | 疎 LU 分解による直接法 |
| `conjugate_gradient` | 共役勾配法(対称正定値) |
| `bicgstab` | BiCGSTAB 法 |
| `gmres` | GMRES 法(`restart` で再開周期を指定) |
| `jacobi` / `gauss_seidel` | 定常反復法 |

```cpp
vcp::linear_solve_options<double> opt;
opt.method = vcp::linear_solver_method::conjugate_gradient;
opt.max_iter = 1000;
opt.tol = 1e-12;
auto r = A.solve_with_info(b, opt);
```

反復法(`conjugate_gradient`・`bicgstab`・`gmres`・`jacobi`・`gauss_seidel`)と
前処理は、まだ実験的な実装です。

`vcp::sparse_lu_options<T>`(疎 LU のオプション)は、次の 4 か所から渡ります。

- `linear_solve_options<T>` の `sparse_lu`(`solve` など)
- `lu_factorize_with_info(opt)` の引数
- `lu_extract_options<T>` の `slu`(「LU 分解」の節)
- `eig_options<T>` の `shift_invert_lu`(「固有値問題」の節)

ふつう使う欄:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 分解法。`auto_select` / `baseline_gp` / `supernodal` / `supernode_panel`。`auto_select` は `supernode_panel` に解決されます(`equilibration = true` のときは `baseline_gp`) |
| `ordering` | `auto_select` | 順序付け。`auto_select` / `natural` / `rcm` / `amd` / `colamd` / `nested_dissection` / `nested_dissection_ml`。`auto_select` は `nested_dissection_ml` に解決されます。将来の版で選択が変わり得るため、再現性が必要な場合は明示してください |
| `pivoting` | `threshold_partial` | ピボット選択。`none` / `diagonal` / `threshold_partial` |
| `equilibration` | `false` | 行と列のスケーリング |
| `iterative_refinement` | `true` | 解いたあとの反復改良 |
| `iterative_refinement_max_iterations` | 20 | 反復改良の回数の上限 |
| `iterative_refinement_tolerance` | 機械イプシロン | 反復改良の目標(相対残差) |

調整用の欄(ふつうは既定値のまま使います):

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `absolute_pivot_tolerance` | 0 | ピボットを受け入れる絶対値の閾値 |
| `zero_tolerance` | 機械イプシロン | 零とみなす閾値 |
| `near_zero_tolerance` | 機械イプシロンの平方根 | 零に近いとみなす閾値 |
| `blas_min_block_size` | 16 | 密ブロックの演算に切り替える最小の幅 |
| `supernode_relaxation` | 0 | スーパーノードの緩和(0 は併合しない) |
| `panel_size` | 8 | パネルの幅 |
| `supernode_panel_maxsup` | 64 | `method = supernode_panel` のスーパーノードの幅の上限 |
| `panel_gemm_kernel` | `blocked` | `method = supernode_panel` の行列積の実装。`reference` / `blocked` |
| `supernodal_numeric_diagnostic_leftlooking` | `false` | `method = supernodal` の数値分解を、診断用の実装に切り替える |
| `supernodal_self_symbolic` | `false` | `method = supernodal` の構造の作り方を切り替える |
| `supernodal_inplace_frontal` | `false` | `method = supernodal` のフロントの更新の方法を切り替える |
| `supernodal_native_check_residual` | `false` | `method = supernodal` で、分解のあとに残差を検査する |
| `ndml_coarsen_stop` | −1 | `ordering = nested_dissection_ml` の調整(負の値はライブラリの既定) |
| `ndml_fm_passes` | −1 | 同上 |
| `ndml_balance_pct` | −1 | 同上 |
| `ndml_leaf_size` | −1 | 同上 |

### 戻り値と status

`vcp::linear_solve_result<T>`:

| 欄 | 意味 |
| --- | --- |
| `x` | 解 |
| `solution` | `x` と同じ内容 |
| `converged` | 直接法では分解と求解が成功したか、反復法では残差の判定を満たしたか |
| `iterations` | 反復法の反復回数。疎 LU では反復改良の回数 |
| `residual_norm` | 残差のノルム(`absolute_residual_norm` と同じ) |
| `absolute_residual_norm` | 残差 ‖b − A x‖ |
| `relative_residual_norm` | 相対残差 |
| `initial_residual_norm` | 最初の残差 |
| `method` | 使った解法 |

- `converged` が `false` のとき、残差の欄の値は意味を持ちません。
  読まないでください。
- 直接法が失敗したとき(特異な行列など)、`x` は零ベクトルです。反復法が
  収束しなかったときは、`x` に途中の近似が入ります。
- 直接法の `converged` は、分解と求解が成功したかどうかで決まります。
  `residual_norm` は参考値であり、残差が小さいことは解が正確であることを
  意味しません。精度が必要な場合は `vcp::spimats<T>` を使ってください。

例外になる場合(`solve_with_info` や `lss_with_info` でも投げます):

| 状況 | 例外 |
| --- | --- |
| 行列が正方でない、右辺の長さが合わない | `vcp::dimension_error` |
| `max_iter = 0`、`tol ≤ 0`、GMRES の `restart = 0` | `vcp::invalid_argument` |
| 共役勾配法に非対称の行列(`check_symmetric = true`) | `vcp::domain_error` |
| Jacobi 法・Gauss–Seidel 法・Jacobi 前処理で対角に零がある | `vcp::numerical_error` |

`solve` などの例外を投げる版は、これに加えて、解けなかったとき(直接法の
分解の失敗、反復法の不収束)に `vcp::state_error` を投げます。

ハンドルの `h.status()` は `vcp::sparse_lu_status` です。

| status | 意味 |
| --- | --- |
| `success` | 分解成功 |
| `structural_singularity` / `numerical_singularity` | 構造的に、または数値的に特異 |
| `zero_pivot` / `near_zero_pivot` / `pivot_rejected` | ピボットが零、零に近い、または受け入れられない |
| `inconclusive_pivot_test` | ピボットの判定を確定できなかった(区間型など) |
| `numerical_instability_suspected` | 数値的に不安定な疑い |
| `unsupported_scalar_type` / `unsupported_pivoting_mode` / `not_implemented` | 対応していない要素型・ピボット選択・オプションの組み合わせ |
| `memory_allocation_failed` | メモリの確保に失敗 |
| `invalid_input` / `internal_error` | 不正な入力 / 内部エラー |

無効なハンドル(`valid()` が `false`)に `solve_with_info` を呼んでも例外は
投げず、`converged = false` と零ベクトルを返します。右辺の長さが合わない
ときは、ハンドルが有効でも無効でも `vcp::dimension_error` です。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `solve` / `solve_with_info` / `lss` / `lss_with_info` | VCP の実装。解法は `method` | UMFPACK に委譲(直接法) | `spumar` と同じ | 保証つき。解の包含を返す |
| 短縮形(`solve_cg` など) | その解法で解く | UMFPACK に委譲(名前の解法は使わない) | `spumar` と同じ | 保証つき(名前の解法は使わない) |
| 複数の右辺 | 疎 LU なら分解 1 回を全部の列で使う | UMFPACK の分解 1 回を全部の列で使う | `spumar` と同じ | 分解 1 回を全部の列で使う |
| `lu_factorize_with_info` | VCP の実装 | VCP の実装(UMFPACK に委譲しない) | `spumar` と同じ | 区間演算で実行(保証の主張なし) |
| 結果の `method` | 使った解法(`auto_select` は `sparse_lu` になる) | 指定した値のまま(既定なら `auto_select`) | `spumar` と同じ | 指定した値のまま |
| 結果の `iterations` | 反復回数、または反復改良の回数 | UMFPACK が行った反復改良の回数 | `spumar` と同じ | 0 |
| 解けなかったとき | `converged = false` | `converged = false`(解は零ベクトル) | `spumar` と同じ | 例外(`_with_info` でも) |

オプションの欄の効き方:

| `linear_solve_options` の欄 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `method` | 効く | 無視される | 無視される | 無視される |
| `max_iter` | 反復法で効く | 0 の検査だけ | 0 の検査だけ | 0 の検査だけ |
| `tol` | 反復法で効く | 正かどうかの検査だけ | 正かどうかの検査だけ | 正かどうかの検査だけ |
| `check_symmetric` | 共役勾配法で効く | 無視される | 無視される | 無視される |
| `use_relative_residual` | 反復法で効く | 無視される | 無視される | 無視される |
| `restart` | GMRES 法で効く | 無視される | 無視される | 無視される |
| `preconditioner` | 共役勾配法・BiCGSTAB 法・GMRES 法で効く | 無視される | 無視される | 無視される |
| `sparse_lu` の全部の欄 | 疎 LU で効く | 無視される(UMFPACK の設定は「`spumar` の詳細」の節) | 無視される | 無視される |

`lu_factorize_with_info(opt)` の `sparse_lu_options` は、4 つの policy とも
VCP の疎 LU に渡ります。

`spimats` の `solve` の詳しい振る舞い(例外の種類、包含の意味)は
「`spimats` の詳細」の節を参照してください。

直接法の内部のアルゴリズムは [kernels/tsparse_lu.md](kernels/tsparse_lu.md) を参照してください。
反復法の内部のアルゴリズムは [kernels/tsparse_iterative.md](kernels/tsparse_iterative.md) を参照してください。

## LU 分解

正方行列 A に対し、行置換 P、列置換 Q、単位下三角行列 L(対角の 1 を明示
格納)、上三角行列 U を、次を満たすように構成します。MATLAB の
`[L,U,P,Q] = lu(A)` と同じ向きです。

```
P A Q = L U
```

- 置換ベクトル p, q は new→old の規約です。MATLAB 表記で `A(p,q) = L*U`。
- 置換行列は `P(k, p[k]) = 1`、`Q(q[k], k) = 1` です。LDLᵀ 分解とは行側の
  向きが転置の関係になります。

```cpp
vcp::spmatrix<double> L, U, P, Q;
std::vector<vcp::spmatrix<double>::index_type> p, q;

A.lu(L, U, P, Q);                  // 置換行列版
A.lu(L, U, p, q);                  // 置換ベクトル版
auto r = A.lu_with_info(L, U, p, q /*, opt*/);
```

### 版の一覧

`L`・`U`・`P`・`Q` は疎行列(出力)、`p`・`q` は
`std::vector<index_type>`(出力)、`opt` は `vcp::lu_extract_options<T>` です。

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.lu(L, U, P, Q)` | `void` | 例外(`vcp::numerical_error`) |
| `A.lu(L, U, P, Q, opt)` | `void` | 例外(`vcp::numerical_error`) |
| `A.lu(L, U, p, q)` | `void` | 例外(`vcp::numerical_error`) |
| `A.lu(L, U, p, q, opt)` | `void` | 例外(`vcp::numerical_error`) |
| `A.lu_with_info(L, U, p, q)` | `lu_extract_result_type` | status |
| `A.lu_with_info(L, U, p, q, opt)` | `lu_extract_result_type` | status |
| `A.lu_with_info(L, U, P, Q)` | `lu_extract_result_type` | status |
| `A.lu_with_info(L, U, P, Q, opt)` | `lu_extract_result_type` | status |
| `lu(A, L, U, p, q)` | `void` | 例外(`vcp::numerical_error`) |

`lu` は status が `success` 以外のとき例外を投げます。`lu_with_info` は
status(`vcp::sparse_lu_extract_status`)を返します。L / U / p / q(P / Q)が
有効なのは `success` のときだけです(失敗のときは空になります)。
自由関数 `lu(A, L, U, p, q)` は `A.lu(L, U, p, q)` と同じです。

取り出した因子を使って解く関数もあります(A⁻¹ = Q U⁻¹ L⁻¹ P)。

```cpp
std::vector<double> x   = A.lu_solve(L, U, p, q, b);        // A x = b
std::vector<double> row = A.lu_inverse_row(L, U, p, q, i);  // A^{-1} の第 i 行

std::vector<double> x2;
vcp::lu_apply_result r1 = A.lu_solve_with_info(L, U, p, q, b, x2);
```

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.lu_solve(L, U, p, q, b)` | `std::vector<T>` | 例外(大きさの不一致は `vcp::dimension_error`、そのほかは `vcp::numerical_error`) |
| `A.lu_solve_with_info(L, U, p, q, b, x)` | `vcp::lu_apply_result` | status |
| `A.lu_inverse_row(L, U, p, q, i)` | `std::vector<T>` | 例外(大きさの不一致と行番号の範囲外は `vcp::dimension_error`、そのほかは `vcp::numerical_error`) |
| `A.lu_inverse_row_with_info(L, U, p, q, i, row)` | `vcp::lu_apply_result` | status |

これらは引数の因子だけを使って計算します(呼び出しに使う行列 A の値は
読みません)。置換ベクトルの版だけがあります。

### オプション

オプション `vcp::lu_extract_options<T>` は、疎 LU のオプション
`sparse_lu_options<T>` をメンバ `slu` に持ちます。

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `slu` | 既定の構築 | 疎 LU のオプション。全部の欄は「線形方程式」の節の「オプション」の表 |

```cpp
vcp::spmatrix<double>::lu_extract_options_type opt;
opt.slu.ordering = vcp::sparse_lu_ordering::colamd;
```

`P A Q = L U` はスケーリングを含まないため、`opt.slu.equilibration = true`
は `unsupported_options` になります。

### 戻り値と status

`lu_extract_result_type`(`vcp::lu_extract_result<T, index_type>`):

| 欄 | 意味 |
| --- | --- |
| `status` | 下の表 |
| `nnz_L` / `nnz_U` | L・U の格納要素の数(`success` のときだけ意味を持つ) |
| `method_used` | 実際に使った分解法(`auto_select` を解決した結果) |
| `ordering_used` | 指定した順序付け(解決した結果ではなく、指定そのまま) |

| status | 意味 |
| --- | --- |
| `success` | 分解成功 |
| `invalid_factorization` | 分解自体が失敗(特異行列など) |
| `unsupported_options` | スケーリングを伴う指定(`equilibration = true` など) |
| `unsupported_storage` | 未対応の格納形式 |
| `internal_error` | 内部エラー |

`vcp::lu_apply_result` は `status`(`vcp::lu_apply_status`)だけを持ちます。
`x`・`row` が有効なのは `success` のときだけです。

| status | 意味 |
| --- | --- |
| `success` | 成功 |
| `dimension_mismatch` | 因子・置換・右辺の大きさ、または行の番号 i が合わない |
| `singular_factor` | U の対角が零 |
| `inconclusive_division` | U の対角が零でないことを確定できない(区間型など) |
| `invalid_input` | 因子の形が規約に合わない(L が単位下三角でない、p・q が置換でない、など) |
| `internal_error` | 内部エラー |

`lu_solve` と `lu_inverse_row` は、status が `success` 以外なら(大きさの
不一致でも)`vcp::numerical_error` を投げます。行列が正方でないときは、
`lu_with_info` でも `vcp::dimension_error` を投げます。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `lu` / `lu_with_info` | VCP の実装 | UMFPACK に委譲 | `spumar` と同じ | 区間演算で実行(保証の主張なし) |
| `lu_solve` / `lu_inverse_row` と `_with_info` | VCP の実装 | VCP の実装(委譲しない) | `spumar` と同じ | 区間演算で実行(保証の主張なし) |
| `slu.method` / `slu.ordering` / `slu.pivoting` | 効く | 無視される(因子は変わらない) | `spumar` と同じ | 効く |
| `slu.equilibration = true` | `unsupported_options` | `unsupported_options` | `unsupported_options` | `unsupported_options` |
| `slu` のほかの欄 | 効く | 無視される | 無視される | 効く |
| 結果の `method_used` | 実際に使った分解法(既定なら `supernode_panel`) | 指定した値のまま(既定なら `auto_select`) | `spumar` と同じ | 実際に使った分解法 |

`spumar`・`spcmodumar` では、因子は UMFPACK が決めた順序付けとピボット選択で
求まります。分解の規約(`P A Q = L U`、p・q は new→old)は同じなので、
返った因子で恒等式は成り立ちます。UMFPACK のピボットの閾値は
`spumar_options` の `umf_pivot_tol` で変えられます(「`spumar` の詳細」の節)。

LU 分解の内部のアルゴリズムは [kernels/tsparse_lu.md](kernels/tsparse_lu.md) を参照してください。

## LDLᵀ 分解と慣性

対称行列 A に対し、置換行列 P、単位下三角行列 L(対角の 1 を明示格納)、
1×1/2×2 ブロック対角行列 D を、次を満たすように構成します
(Bunch–Kaufman 型の対称ピボット)。

```
Pᵀ A P = L D Lᵀ
```

- 数値分解が読むのは下三角部(対角を含む)だけです。上三角は対称性の検査に
  だけ使います。
- 置換ベクトル p は new→old の規約です。MATLAB 表記で `A(p,p) = L*D*L'`。
  置換行列は `P(p[k], k) = 1` です。
- D の 2×2 ブロックは副対角を両側に格納します。零ピボットの対角は格納
  しません。

```cpp
vcp::spmatrix<double> L, D, P;
std::vector<vcp::spmatrix<double>::index_type> p;

A.ldl(L, D, P);
A.ldl(L, D, p);
auto r = A.ldl_with_info(L, D, p /*, opt*/);
```

### 版の一覧

`L`・`D`・`P` は疎行列(出力)、`p` は `std::vector<index_type>`(出力)です。

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.ldl(L, D, P)` | `void` | 例外(`vcp::numerical_error`) |
| `A.ldl(L, D, P, opt)` | `void` | 例外(`vcp::numerical_error`) |
| `A.ldl(L, D, p)` | `void` | 例外(`vcp::numerical_error`) |
| `A.ldl(L, D, p, opt)` | `void` | 例外(`vcp::numerical_error`) |
| `A.ldl_with_info(L, D, p)` | `ldl_result_type` | status |
| `A.ldl_with_info(L, D, p, opt)` | `ldl_result_type` | status |
| `A.ldl_with_info(L, D, P)` | `ldl_result_type` | status |
| `A.ldl_with_info(L, D, P, opt)` | `ldl_result_type` | status |

`opt` は `vcp::ldl_options<T>` です。`ldl` は status が `success` 以外のとき
(`zero_pivot` のときも)例外を投げます。`ldl_with_info` は
status(`vcp::sparse_ldl_status`)を返します。L / D / p(P)が有効なのは
`success` または `zero_pivot` のときだけです。

慣性(正・負・零の固有値の個数 n₊, n₋, n₀)は LDLᵀ 分解から求めます。

```cpp
auto ir  = A.inertia();             // 失敗時に例外
auto ir2 = A.inertia_with_info();   // status を返す
```

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.inertia()` | `inertia_result_type` | 例外(`vcp::numerical_error`) |
| `A.inertia(opt)` | `inertia_result_type` | 例外(`vcp::numerical_error`) |
| `A.inertia_with_info()` | `inertia_result_type` | status |
| `A.inertia_with_info(opt)` | `inertia_result_type` | status |

ここでの `opt` は `vcp::inertia_options<T>` です。

A − σB(B を省略すると A − σI)の慣性を、σ を変えて繰り返し求める場合は、
σ に依存しない処理を 1 回で済ませるハンドルを使います。

```cpp
auto h = A.ldl_shift_setup_with_info(B);
if (h.valid()) {
    auto i1 = h.inertia_at(sigma1);
    auto i2 = h.inertia_at(sigma2);
}
```

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.ldl_shift_setup_with_info()` | `ldl_shift_handle_type` | ハンドルの `valid()` が `false` |
| `A.ldl_shift_setup_with_info(opt)` | `ldl_shift_handle_type` | 同上 |
| `A.ldl_shift_setup_with_info(B)` | `ldl_shift_handle_type` | 同上 |
| `A.ldl_shift_setup_with_info(B, opt)` | `ldl_shift_handle_type` | 同上 |
| `vcp::ldl_shift_default_options<T>()` | `vcp::ldl_options<T>` | --- |
| `h.valid()` | `bool` | --- |
| `h.status()` | `vcp::sparse_ldl_status` | --- |
| `h.n()` | `index_type` | --- |
| `h.has_B()` | `bool` | --- |
| `h.symmetry_failed_on()` | `int` | --- |
| `h.method_used()` | `vcp::sparse_ldl_method` | --- |
| `h.ordering_used()` | `vcp::sparse_ldl_ordering` | --- |
| `h.nnz_merged()` | `index_type` | --- |
| `h.inertia_at(sigma)` | `ldl_shift_inertia_result_type` | status(`inertia.status`) |
| `h.inertia_at(sigma, ws)` | `ldl_shift_inertia_result_type` | 同上 |
| `h.inertia_at(sigma, ws, zero_tol)` | `ldl_shift_inertia_result_type` | 同上 |
| `h.factor_at(sigma, L, D, p)` | `vcp::ldl_result<T, index_type>` | status |
| `h.factor_at(sigma, L, D, p, ws)` | `vcp::ldl_result<T, index_type>` | status |
| `h.factor_at_with(sigma, method, pivoting, L, D, p)` | `vcp::ldl_result<T, index_type>` | status |
| `h.factor_at_with(sigma, method, pivoting, L, D, p, ws)` | `vcp::ldl_result<T, index_type>` | status |

- `ldl_shift_setup_with_info` の `opt` は `vcp::ldl_options<T>` です。省略
  したときの既定は `ldl` の既定と違い、`vcp::ldl_shift_default_options<T>()`
  (`method = supernodal`、`pivoting = none`)です。
- `ws` は作業領域(`ldl_shift_workspace_type`)です。同じ `ws` を渡し続けると、
  σ を変えるたびの記憶領域の確保を省けます。`zero_tol` は D の要素を零と
  みなす閾値(既定 0)です。
- `h.symmetry_failed_on()` は、対称性の検査で失敗した行列を返します(0: なし、
  1: A、2: B)。`h.nnz_merged()` は A と B を重ねたパターンの要素数です。
- `h.factor_at` は、A − σB の因子 L・D・p を `ldl_with_info` と同じ規約で
  取り出します。`h.factor_at_with` は、分解法(`vcp::sparse_ldl_method`)と
  ピボット選択(`vcp::sparse_ldl_pivoting`)を指定してやり直します。
- 既定の `pivoting = none` は対角の入れ替えをしないので、σ が固有値に一致
  するか、ごく近いと、`inertia_at` は失敗します(`inertia.status` が
  `factorization_failed`、`inertia.ldl_status` が `static_pivot_breakdown`
  など)。そのときは σ を少しずらすか、`h.factor_at_with` で
  `baseline_dynamic` と `bk` を指定してください。

### オプション

`vcp::ldl_options<T>`。ふつう使う欄:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 分解法。`auto_select` / `baseline_dynamic` / `supernodal`。`auto_select` は `baseline_dynamic` に解決されます |
| `pivoting` | `bk` | ピボット選択。`bk`(Bunch–Kaufman。1×1/2×2 と対称な入れ替え)/ `none`(入れ替えなし) |
| `ordering` | `auto_select` | 順序付け。`auto_select` / `natural` / `rcm` / `amd` / `nested_dissection` / `nested_dissection_ml`。`auto_select` は `nested_dissection_ml` に解決されます。選ばれた順序付けは `ordering_used` に記録されます。将来の版で選択が変わり得るため、再現性が必要な場合は明示してください |
| `check_symmetry` | `true` | 対称性の検査。下三角だけを与える場合は `false` |
| `symmetry_tol` | 1e-12 | 対称性の許容誤差 |
| `zero_pivot_tol` | 0 | 零ピボットとみなす閾値 |

調整用の欄(ふつうは既定値のまま使います):

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `ldl_min_block_size` | 32 | `method = supernodal` で密ブロックの演算を使う最小の幅 |
| `diag_kernel` | `auto_select` | `method = supernodal` の対角ブロックの更新の実装。`auto_select` / `gemmtr` / `gemm` |
| `ldl_amalgamation` | `off` | `method = supernodal` のスーパーノードの併合。`off` / `relaxed` |
| `ldl_amalg_n0` | 4 | 併合の規則の幅(1 段目) |
| `ldl_amalg_n1` | 16 | 併合の規則の幅(2 段目) |
| `ldl_amalg_n2` | 48 | 併合の規則の幅(3 段目) |
| `ldl_amalg_z0` | 80 | 併合の規則の零の割合(%、1 段目) |
| `ldl_amalg_z1` | 10 | 併合の規則の零の割合(%、2 段目) |
| `ldl_amalg_z2` | 5 | 併合の規則の零の割合(%、3 段目) |
| `pivot_threshold` | 0 | 予約。0 以外を指定すると `invalid_options` |
| `ndml_coarsen_stop` | −1 | `ordering = nested_dissection_ml` の調整(負の値はライブラリの既定) |
| `ndml_fm_passes` | −1 | 同上 |
| `ndml_balance_pct` | −1 | 同上 |
| `ndml_leaf_size` | −1 | 同上 |

`vcp::inertia_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `zero_tol` | 0 | D の要素を零とみなす閾値(`ldl.zero_pivot_tol` とは別の値) |
| `ldl` | 既定の構築 | 内部の LDLᵀ 分解のオプション(`ldl_options<T>`) |

### 戻り値と status

`ldl_result_type`(`vcp::ldl_result<T, index_type>`)の主な欄:

| 欄 | 意味 |
| --- | --- |
| `status` | 下の表 |
| `n_pivots_1x1` / `n_pivots_2x2` | 1×1・2×2 のピボットの数 |
| `first_zero_pivot` | 最初の零ピボットの位置(なければ −1) |
| `inconclusive_at` | 判定を確定できなかった列(なければ −1) |
| `structural_empty_at` | 構造的に空の行・列(なければ −1) |
| `breakdown_at` | `static_pivot_breakdown` で止まった列(なければ −1) |
| `nnz_L` | L の格納要素の数 |
| `ordering_used` / `method_used` | 実際に使った順序付け・分解法 |
| `pivot_mode_used` | 実際に使ったピボット選択 |

診断の欄:

| 欄 | 意味 |
| --- | --- |
| `dense_delegated` | 小さい行列で、密行列の実装に任せたか |
| `diag_kernel_used` | `supernodal` で使った対角ブロックの実装 |
| `n_supernodes` / `max_supernode_width` / `n_boundary_splits` | スーパーノードの数・最大の幅・2×2 ピボットで分かれた数 |
| `nnz_L_static` | 入れ替えのない分解での L の要素数 |
| `n_zero_skips` | 飛ばした零ピボットの数 |
| `out_of_panel_at` | パネルの外との入れ替えが必要になった列(なければ −1) |
| `gemm_call_count` / `gemm_time_ns` | 密ブロックの演算の回数と時間 |
| `growth_log2` / `growth_valid` | 要素の成長の目安(`growth_valid` が `true` のときだけ有効) |

| status | 意味 |
| --- | --- |
| `success` | 分解完了 |
| `structural_singularity` | 構造的に空の行・列がある(`structural_empty_at`) |
| `zero_pivot` | 零ピボットを飛ばして分解を最後まで続けた(`first_zero_pivot`) |
| `inconclusive_pivot_test` | ピボットの判定を確定できなかった(区間型で 0 を含む要素など。`inconclusive_at`) |
| `pivot_out_of_panel` | `method = supernodal` で、パネルの外の列との入れ替えが必要になり、分解を止めた(`out_of_panel_at`) |
| `not_symmetric` | 対称であることを確認できない |
| `invalid_options` / `invalid_input` | 不正なオプション・入力 |
| `internal_error` | 内部エラー |
| `static_pivot_breakdown` | `pivoting = none` で、零ピボットの列の対角より下に零でない要素があり、分解を止めた(`breakdown_at`) |

`method = supernodal`・`pivoting = bk` を指定して `ldl` / `ldl_with_info` を
呼ぶと、パネルの外との入れ替えが必要になったときに `baseline_dynamic` で
やり直します。そのとき `method_used` は `baseline_dynamic` になります。
`pivot_out_of_panel` が利用者に返るのは、ハンドルの `factor_at` と
`inertia_at`(`inertia.ldl_status`)です。

`inertia_result_type`(`vcp::inertia_result<index_type>`):

| 欄 | 意味 |
| --- | --- |
| `n_pos` / `n_neg` / `n_zero` | 正・負・零の固有値の個数 |
| `status` | 下の表(`vcp::inertia_status`) |
| `inconclusive_at` | 符号を確定できなかった位置(なければ −1) |
| `ldl_status` | 内部の LDLᵀ 分解の status |

| status | 意味 |
| --- | --- |
| `success` | 成功 |
| `factorization_failed` | 内部の LDLᵀ 分解が `success` でも `zero_pivot` でもない(理由は `ldl_status`) |
| `inconclusive_sign` | 符号を確定できなかった |
| `not_block_diagonal` | D がブロック対角でない |
| `invalid_input` / `internal_error` | 不正な入力 / 内部エラー |

分解が `zero_pivot` の場合は失敗ではなく、零固有値として n₀ に数えます
(このとき `inertia()` も例外を投げません)。構造的に空の行・列がある行列
(`structural_singularity`)と、対称でない行列(`not_symmetric`)では、慣性は
求まりません(`factorization_failed`)。

`ldl_shift_inertia_result_type` は、慣性の結果 `inertia`(上の
`inertia_result`)と、`growth_log2`・`growth_valid` を持ちます。

ハンドルの `h.status()` は設定の結果です。対称でない行列を渡すと
`not_symmetric` になり、`h.valid()` は `false` です。無効なハンドルに
`inertia_at` を呼んでも例外は投げず、`inertia.status` が `invalid_input` に
なります。行列が正方でない、B の大きさが合わないなどの誤りは、
`_with_info` の関数でも `vcp::dimension_error` を投げます。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `ldl` / `ldl_with_info` | VCP の実装 | VCP の実装(委譲しない) | VCP の実装(委譲しない) | 区間演算で実行(保証の主張なし) |
| `inertia` / `inertia_with_info` | VCP の実装 | VCP の実装 | VCP の実装 | 区間演算で実行(保証の主張なし) |
| `ldl_shift_setup_with_info` とハンドル | VCP の実装 | VCP の実装 | VCP の実装 | 区間演算で実行(保証の主張なし) |
| `ldl_options` / `inertia_options` の全部の欄 | 効く | 効く | 効く | 効く |
| 既定の `ordering_used` / `method_used` | `nested_dissection_ml` / `baseline_dynamic` | 同じ | 同じ | 同じ |

`spimats` では、幅のない区間(点)だけからなる行列なら `spmats` と同じように
分解が進みますが、幅のある区間を含む行列では、対称性の検査が通らず
`not_symmetric` になります(慣性は `factorization_failed`)。

LDLᵀ 分解と慣性の内部のアルゴリズムは [kernels/tsparse_ldl.md](kernels/tsparse_ldl.md) を参照してください。

## Cholesky 分解

対称正定値行列 A に対し、非単位下三角行列 L(対角は正)を、次を満たすように
構成します。MATLAB の `chol(A,'lower','vector')` と同じ恒等式です。

```
Pᵀ A P = L Lᵀ
```

置換ベクトル p は new→old の規約(`A(p,p) = L*L'`)、置換行列は
`P(p[k], k) = 1` です。

```cpp
vcp::spmatrix<double> L, P;
std::vector<vcp::spmatrix<double>::index_type> p;

A.chol(L, p);
auto r = A.chol_with_info(L, p /*, opt*/);
if (r.status == vcp::sparse_chol_status::success) { /* L, p が有効 */ }

vcp::spmatrix<double> L1 = Cholesky(A);      // A = L1 L1ᵀ(順序付けなし)
vcp::spmatrix<double> L2 = Cholesky(A, p);   // 既定の順序付け
```

### 版の一覧

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.chol(L, p)` | `void` | 例外(`vcp::numerical_error`) |
| `A.chol(L, p, opt)` | `void` | 例外(`vcp::numerical_error`) |
| `A.chol(L, P)` | `void` | 例外(`vcp::numerical_error`) |
| `A.chol(L, P, opt)` | `void` | 例外(`vcp::numerical_error`) |
| `A.chol_with_info(L, p)` | `chol_result_type` | status |
| `A.chol_with_info(L, p, opt)` | `chol_result_type` | status |
| `A.chol_with_info(L, P)` | `chol_result_type` | status |
| `A.chol_with_info(L, P, opt)` | `chol_result_type` | status |
| `Cholesky(A)` | `spmatrix` | 例外(`vcp::numerical_error`) |
| `Cholesky(A, p)` | `spmatrix` | 例外(`vcp::numerical_error`) |
| `A.is_spd()` | `bool` | `false` |
| `A.is_spd(opt)` | `bool` | `false` |

- `chol` の `opt` は `vcp::chol_options<T>`、`is_spd` の `opt` は
  `vcp::spd_options<T>` です。
- `chol` は status が `success` 以外のとき例外を投げます。`chol_with_info` は
  status(`vcp::sparse_chol_status`)を返します。L / p(P)が有効なのは
  `success` のときだけです(失敗のときは空になります)。
- `Cholesky(A)` は順序付けなし(`ordering = natural`)で分解し、
  `L Lᵀ = A` となる L を返します。順序付けをしないので、L の要素数は
  既定の順序付けより多くなることがあります。`Cholesky(A, p)` は既定の
  オプションの `A.chol(L, p)` と同じ L を返します。
- `A.is_spd()` は、行列が対称正定値かどうかを返します。`spmats`・`spumar`・
  `spcmodumar` では「Cholesky 分解が最後まで成功したか」という近似の判定で、
  `spimats` では証明です(「policy ごとの違い」)。

### オプション

`vcp::chol_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 分解法。指定できるのは `auto_select` と `simplicial_uplooking` です。`supernodal` は `method_used` に記録されるだけの値で、指定すると `invalid_options` になります |
| `ordering` | `auto_select` | 順序付け。`auto_select` / `natural` / `rcm` / `amd` / `nested_dissection` / `nested_dissection_ml`。`auto_select` の解決先は policy で違います(「policy ごとの違い」)。選ばれた順序付けは `ordering_used` に記録されます。`natural` を指定すると置換なしで分解します |
| `check_symmetry` | `true` | 対称性の検査。下三角だけを与える場合は `false` |
| `symmetry_tol` | 1e-12 | 対称性の許容誤差 |
| `pd_tol` | 0 | ピボットがこの値より大きいときに正とみなす |
| `ndml_coarsen_stop` | −1 | `ordering = nested_dissection_ml` の調整(負の値はライブラリの既定) |
| `ndml_fm_passes` | −1 | 同上 |
| `ndml_balance_pct` | −1 | 同上 |
| `ndml_leaf_size` | −1 | 同上 |

`vcp::spd_options<T>`(`is_spd` のオプション):

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `nosymcheck` | `false` | 入口の対称性の検査を省く(対称であることは呼び出す側が責任を持つ) |
| `ordering` | `auto_select` | 内部の Cholesky 分解の順序付け |

### 戻り値と status

`chol_result_type`(`vcp::chol_result<T, index_type>`):

| 欄 | 意味 |
| --- | --- |
| `status` | 下の表 |
| `failure_at` | 正定値でないと分かった位置(なければ −1。MATLAB の flag は `failure_at + 1`) |
| `inconclusive_at` | ピボットの符号を確定できなかった位置(なければ −1) |
| `structural_empty_at` | 構造的に空の行・列(なければ −1) |
| `nnz_L` | L の要素数 |
| `ordering_used` / `method_used` | 実際に使った順序付け・分解法 |

| status | 意味 |
| --- | --- |
| `success` | 分解成功 |
| `not_positive_definite` | 正定値でない(`failure_at`) |
| `inconclusive_pivot_test` | ピボットの符号を確定できなかった(`inconclusive_at`) |
| `not_symmetric` | 対称であることを確認できない |
| `invalid_input` / `invalid_options` | 不正な入力・オプション |
| `internal_error` | 内部エラー |

`is_spd` は `bool` だけを返します。

- 行列が正方でなければ `vcp::dimension_error`、対称でなければ
  `vcp::domain_error` を投げます(`nosymcheck = true` なら対称性は検査
  しません)。0×0 の行列は `true` です。
- `spmats`・`spumar`・`spcmodumar` の `false` は、「正定値でない」場合と
  「正定値だが浮動小数点の分解が途中で破綻した」場合を区別しません。
  `true` も証明ではありません。

行列が正方でないときは、`chol_with_info` でも `vcp::dimension_error` を
投げます。

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `chol` / `chol_with_info` / `Cholesky` | VCP の実装 | VCP の実装(委譲しない) | CHOLMOD に委譲 | 区間演算で実行(保証の主張なし) |
| `is_spd` | VCP の実装の Cholesky 分解が成功したか(近似の判定) | `spmats` と同じ | CHOLMOD の Cholesky 分解が成功したか(近似の判定) | 保証つき。`true` は、区間内のすべての対称行列が正定値であることの証明 |
| `ordering = auto_select` の解決先 | `nested_dissection_ml` | `nested_dissection_ml` | `amd` | `nested_dissection_ml` |
| `ordering = nested_dissection_ml` | 効く | 効く | 効く(VCP の順序を CHOLMOD に与える) | 効く |
| `ordering` のほかの値 | 効く | 効く | CHOLMOD の対応する設定に渡る(「`spcmodumar` の詳細」の節) | 効く |
| `method = simplicial_uplooking` | 効く | 効く | CHOLMOD の simplicial の分解を指定する | 効く |
| `check_symmetry` / `symmetry_tol` | 効く | 効く | 効く | 効く |
| `pd_tol` | 効く | 効く | 0 以外は `invalid_options` | 効く |
| `ndml_*` | `nested_dissection_ml` で効く | 同じ | 同じ | `spmats` と同じ |
| `method_used` | `simplicial_uplooking` | `simplicial_uplooking` | `simplicial_uplooking` または `supernodal`(CHOLMOD が選ぶ) | `simplicial_uplooking` |
| `spd_options::ordering` | 内部の Cholesky 分解に渡る | 同じ | CHOLMOD への委譲に渡る | 判定の前に行列へ適用する |

`spimats` の `is_spd` の `false` は「正定値でない」という主張ではなく、
十分条件が成り立たなかったことを意味します(「`spimats` の詳細」の節)。
`spimats` の `chol` は、幅のある区間を含む行列では対称性の検査が通らず
`not_symmetric` になります。

Cholesky 分解の内部のアルゴリズムは [kernels/tsparse_ldl.md](kernels/tsparse_ldl.md) を参照してください。

## 近似逆行列(AINV・FSAI)

反復法の前処理などに使う近似逆行列を構成します。近似逆行列 R そのものは
構成せず、ベクトルへの作用 `apply` と、‖I − R A‖∞ の推定値(保証なし)を
返す関数を提供します。

| 手法 | 構成 | 因子 | 作用 |
| --- | --- | --- | --- |
| AINV | `ainv_with_info` | Z, W は単位上三角、D は対角 | `ainv_apply`: z = Z D⁻¹ Wᵀ r |
| FSAI(対称正定値) | `fsai_with_info` | U は単位上三角、D は対角 | `fsai_apply`: z = P U D⁻¹ Uᵀ Pᵀ r |
| 適応 FSAI | `fsai_adaptive_with_info` | FSAI と同じ | `fsai_apply` を使う |

### 版の一覧

`Z`・`W`・`D`・`U` は疎行列(構成の関数では出力)、`perm` は
`std::vector<index_type>`、`r`・`z` は `std::vector<T>`、`est` は実数の変数
(出力)です。

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.ainv_with_info(Z, W, D)` | `ainv_result_type` | status |
| `A.ainv_with_info(Z, W, D, opt)` | `ainv_result_type` | status |
| `A.ainv_residual_norm_estimate(Z, W, D, est)` | `vcp::ainv_status` | status |
| `A.ainv_apply(Z, W, D, r, z)` | `vcp::ainv_status` | status |
| `A.fsai_with_info(U, D, perm)` | `fsai_result_type` | status |
| `A.fsai_with_info(U, D, perm, opt)` | `fsai_result_type` | status |
| `A.fsai_residual_norm_estimate(U, D, perm, est)` | `vcp::fsai_status` | status |
| `A.fsai_apply(U, D, perm, r, z)` | `vcp::fsai_status` | status |
| `A.fsai_adaptive_with_info(U, D, perm)` | `fsai_adaptive_result_type` | status |
| `A.fsai_adaptive_with_info(U, D, perm, opt)` | `fsai_adaptive_result_type` | status |
| `A.fsai_adaptive_with_info(U, D, perm, U0, perm0)` | `fsai_adaptive_result_type` | status |
| `A.fsai_adaptive_with_info(U, D, perm, U0, perm0, opt)` | `fsai_adaptive_result_type` | status |

- いずれも status を返し、例外は投げません(行列が正方でないときも
  `invalid_input` を返します)。例外を投げる版はありません。
- `perm` は new→old の置換ベクトルで、置換行列は `P(perm[k], k) = 1` です。
- `U0`・`perm0` を取る版は、すでにある因子(たとえば `fsai_with_info` の
  結果)を初期値にして適応 FSAI を行います。`perm0` をそのまま使い、
  `opt.ordering` は無視されます。
- 推定値 `est` は ‖I − R A‖∞ の近似で、保証はありません。

### オプション

`vcp::ainv_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `drop_tolerance` | 1e-2 | 絶対値がこれより小さい fill-in を捨てる(0 なら捨てない) |
| `max_nnz_per_column` | 0 | 因子の 1 列の要素数の上限(0 は無制限。対角の 1 は数えない) |
| `pivot_small_threshold` | 1e-14 | 絶対値がこれ以下のピボットを持ち上げる |
| `pivot_lift_value` | 1e-3 | 持ち上げたあとのピボットの大きさ |

`vcp::fsai_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `pattern_power` | 1 | 因子のパターンを A の何乗のパターンから作るか |
| `prefilter_tolerance` | 0 | パターンを作る前に、相対的に小さい要素を捨てる閾値(0 は捨てない) |
| `max_density` | 5 | パターンの要素数の上限(A の要素数の何倍まで。0 は無制限) |
| `ordering` | `natural` | 順序付け(`vcp::sparse_chol_ordering`)。`auto_select` は `natural` に解決されます |
| `postfilter_tolerance` | 0 | 求めた因子から、相対的に小さい要素を捨てる閾値(0 は捨てない) |
| `postfilter_max_nnz_per_row` | 0 | 因子の 1 行の要素数の上限(0 は無制限) |
| `pivot_small_threshold` | 1e-14 | 絶対値がこれ以下のピボットを持ち上げる |
| `pivot_lift_value` | 1e-3 | 持ち上げたあとのピボットの大きさ |

`vcp::fsai_adaptive_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `max_iterations` | 10 | 1 行あたりの、パターンを広げる回数の上限 |
| `entries_per_step` | 5 | 1 回に足す要素の候補の数 |
| `exit_tolerance` | 1e-3 | 行の更新を打ち切る閾値 |
| `inner_drop_tolerance` | 0 | 更新の途中で、相対的に小さい要素を捨てる閾値(0 は捨てない) |
| `ordering` | `natural` | 順序付け。`U0`・`perm0` を取る版では無視されます |
| `postfilter_tolerance` | 0 | 最後に、相対的に小さい要素を捨てる閾値(0 は捨てない) |
| `postfilter_max_nnz_per_row` | 0 | 因子の 1 行の要素数の上限(0 は無制限) |
| `pivot_small_threshold` | 1e-14 | 絶対値がこれ以下のピボットを持ち上げる |
| `pivot_lift_value` | 1e-3 | 持ち上げたあとのピボットの大きさ |

### 戻り値と status

`vcp::ainv_status` と `vcp::fsai_status` は、どちらも次の 3 つの値です。

| status | 意味 |
| --- | --- |
| `success` | 出力は有効 |
| `invalid_input` | 行列が正方でない、因子やベクトルの大きさが合わない、D が対角でない、`perm` が置換でない、など |
| `internal_error` | 内部エラー |

小さいピボットや特異な行列は失敗にはなりません。ピボットを持ち上げて
構成を続け、持ち上げた回数を結果に記録します。

`ainv_result_type`(`vcp::ainv_result<T, index_type>`):

| 欄 | 意味 |
| --- | --- |
| `status` | 上の表 |
| `n_pivot_modifications` | ピボットを持ち上げた回数 |
| `first_modified_pivot` | 最初に持ち上げた列(なければ −1) |
| `nnz_Z` / `nnz_W` | Z・W の要素数 |

`fsai_result_type`(`vcp::fsai_result<T, index_type>`):

| 欄 | 意味 |
| --- | --- |
| `status` | 上の表 |
| `nnz_U` | U の要素数 |
| `nnz_pattern` | 小さい要素を捨てる前のパターンの要素数 |
| `n_pivot_modifications` | ピボットを持ち上げた回数 |
| `first_modified_row` | 最初に持ち上げた行(なければ −1) |
| `n_postfilter_dropped` | あとから捨てた要素の数 |
| `n_postfilter_restore_skipped` | 捨てたあとの D の補正を見送った行の数 |

`fsai_adaptive_result_type`(`vcp::fsai_adaptive_result<T, index_type>`):

| 欄 | 意味 |
| --- | --- |
| `status` | 上の表 |
| `nnz_U` | U の要素数 |
| `n_pivot_modifications` | ピボットを持ち上げた回数 |
| `first_modified_row` | 最初に持ち上げた行(なければ −1) |
| `n_postfilter_dropped` | あとから捨てた要素の数 |
| `n_postfilter_restore_skipped` | 捨てたあとの D の補正を見送った行の数 |
| `n_rows_converged` | `exit_tolerance` で打ち切った行の数 |
| `n_rows_capped` | `max_iterations` まで回した行の数 |

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `ainv_*` / `fsai_*` / `fsai_adaptive_*` | VCP の実装 | VCP の実装(委譲しない) | VCP の実装(委譲しない) | 区間演算で実行(保証の主張なし) |
| オプションの全部の欄 | 効く | 効く | 効く | 効く |

近似逆行列の内部のアルゴリズムは [kernels/sparse_approx_inverse.md](kernels/sparse_approx_inverse.md) を参照してください。

## 固有値問題

疎行列の固有値問題では、求める個数 k を指定します。結果を密行列で受ける
自由関数の `eigs` が、ふつうの書き方です。

```cpp
vcp::matrix<double, vcp::mats2<double> > E, V;

eigs(A, E, k);              // 固有値。E は k x 1
eigs(A, E, V, k, opt);      // 固有値と固有ベクトル。V は n x k
eigs(A, B, E, V, k, opt);   // 一般化固有値問題 A x = λ B x
```

収束したかどうかを確かめたいとき、収束の情報がほしいとき、前処理を渡したい
ときは、メンバ関数を使います。

```cpp
std::vector<double> lam = A.eigs(k);     // 収束しなければ例外
auto r = A.eigs_with_info(k, opt);       // r.eigenvalues, r.eigenvectors, r.converged

auto rg = A.eigs_with_info(B, k, opt);   // 一般化固有値問題 A x = λ B x
```

自由関数の `eigs` とメンバの `A.eigs` は、収束しなかったとき、または返った
固有値が min(k, n) 個に足りないときに、例外(`vcp::state_error`)を投げます。
例外のとき、E と V は呼び出しの前のままです。収束しなかった理由や途中の値を
見たいときは、`A.eigs_with_info` を使ってください。

### 版の一覧

`E`・`V` は密行列 `vcp::matrix`(出力)、`B` は疎行列、`opt` は
`vcp::eig_options<T>`、`M` は前処理です。

自由関数(結果を密行列で受ける):

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `eigs(A, E, k)` | `void` | 例外(`vcp::state_error`) |
| `eigs(A, E, k, opt)` | `void` | 同上 |
| `eigs(A, E, V, k)` | `void` | 同上 |
| `eigs(A, E, V, k, opt)` | `void` | 同上 |
| `eigs(A, B, E, k)` | `void` | 同上 |
| `eigs(A, B, E, k, opt)` | `void` | 同上 |
| `eigs(A, B, E, V, k)` | `void` | 同上 |
| `eigs(A, B, E, V, k, opt)` | `void` | 同上 |

E は(返った個数)×1、V は n×(返った固有ベクトルの本数)で、V の第 j 列が
E の第 j 行の固有値の固有ベクトルです。E は min(k, n)×1、V は n×min(k, n) です
(`spimats` は固有ベクトルを返さないので、V は n×0)。
自由関数は、対応するメンバ `eigs_with_info` を呼び、メンバの `A.eigs` と同じ
判定(収束していて、min(k, n) 個が返っていること)を通ったときだけ、その
`eigenvalues` と `eigenvectors` を E と V に詰めます。判定を通らなければ
例外を投げ、E と V は書き換えません。

メンバ関数:

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.eigs(k)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.eigs(k, opt)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.eigs_with_info(k)` | `vcp::eig_result<T>` | `converged = false` と `status` |
| `A.eigs_with_info(k, opt)` | `vcp::eig_result<T>` | `converged = false` と `status` |
| `A.eigs(B, k)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.eigs(B, k, opt)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.eigs_with_info(B, k)` | `vcp::eig_result<T>` | `converged = false` と `status` |
| `A.eigs_with_info(B, k, opt)` | `vcp::eig_result<T>` | `converged = false` と `status` |

`A.eigs` は、収束しなかったとき、または返った固有値が min(k, n) 個に
足りないときに例外を投げます。

全部の固有値がほしいときは、k に n を指定します(`A.eigs_with_info(n)`)。
n が 300 以下なら、既定のオプションで密行列の解法が選ばれます。

前処理 M を渡す版(shift-invert の内側の反復に M を使う):

| 呼び方 | 戻り値の型 | 失敗のとき |
| --- | --- | --- |
| `A.eigs(k, opt, M)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.eigs_with_info(k, opt, M)` | `vcp::eig_result<T>` | `converged = false` と `status` |
| `A.eigs(B, k, opt, M)` | `std::vector<T>` | 例外(`vcp::state_error`) |
| `A.eigs_with_info(B, k, opt, M)` | `vcp::eig_result<T>` | `converged = false` と `status` |

M は、`apply(r, z)`(`r` は `const std::vector<T>&`、`z` は
`std::vector<T>&`。`z` に前処理をかけた結果を書く)を持つ型なら何でも
渡せます。既製のものは次の 4 つです。

| クラス | 作り方 | 内容 |
| --- | --- | --- |
| `vcp::identity_preconditioner<T>` | `vcp::identity_preconditioner<double> M;` | 何もしない |
| `vcp::jacobi_preconditioner<SpMat>` | `vcp::jacobi_preconditioner<vcp::spmatrix<double> > M(K);` | K の対角の逆数を掛ける |
| `vcp::ilu0_preconditioner<SpMat>` | `vcp::ilu0_preconditioner<vcp::spmatrix<double> > M(K);` | K の ILU(0) 分解 |
| `vcp::function_preconditioner<T, F>` | `vcp::function_preconditioner<double, F> M(f);` | 関数オブジェクト `f(r, z)` を呼ぶ |

M が使われるのは shift-invert のとき(`use_shift = true`、または `method` が
`shift_invert_lanczos` / `shift_invert_arnoldi`)だけです。それ以外の指定で
M を渡すと、計算せずに `status = "unsupported_preconditioner_for_method"` を
返します。M は A − σI(一般化では A − σB)の連立一次方程式を反復法で解く
ための前処理なので、その行列に合ったものを渡してください。弱い前処理
(恒等、Jacobi)では内側の反復が収束せず、`status = "inner_solve_failed"` に
なることがあります。

既定のオプションを返す関数:

| 呼び方 | 戻り値の型 |
| --- | --- |
| `A.default_eigs_options()` | `eig_options_type` |

`A.default_eigs_options()` は既定のオプション(`eig_options_type()` と同じ
内容)を返します。

### オプション

`vcp::eig_options<T>`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 解法(下の表) |
| `structure` | `auto_detect` | 構造の指定。`auto_detect` / `general` / `symmetric` / `hermitian` |
| `target` | `smallest_algebraic` | どの固有値を求めるか(下の表) |
| `tol` | 1e-12 | 収束判定。0 以下は例外 |
| `shift` | 0 | シフト σ(shift-invert、および `target_magnitude` / `target_real` の目標) |
| `use_shift` | `false` | shift-invert を使う |
| `max_iter` | 1000 | 反復の上限(行列とベクトルの積の回数。密行列の解法では回転の回数)。0 は例外 |
| `subspace_dim` | 0 | Krylov 部分空間の次元(0 は自動) |
| `allow_dense_conversion` | `true` | 密行列に変換して解くことを許す |
| `max_dense_size` | 1000000 | 密行列に変換してよい n² の上限 |
| `orthogonalization` | `modified_gram_schmidt` | 直交化。`modified_gram_schmidt` / `classical_gram_schmidt_twice` |
| `random_seed` | 0 | 開始ベクトルの乱数の種 |
| `random_start` | `false` | 開始ベクトルを乱数にする |
| `compute_residual_history` | `false` | 残差の履歴を記録する |
| `shift_invert_solver` | `sparse_lu` | shift-invert の内側の解法。`sparse_lu`(直接法)/ `ilu0_gmres`(ILU(0) 前処理つき GMRES) |
| `shift_invert_lu` | 既定の構築 | shift-invert の疎 LU のオプション(`sparse_lu_options<T>`。「線形方程式」の節) |
| `allow_complex_pairs` | `false` | 実行列の複素共役固有値の対を返すか |

`target`(`vcp::eig_target`):

| `target` | 求める固有値 | 返る順 |
| --- | --- | --- |
| `largest_magnitude` | 絶対値の大きいもの | 大きい順 |
| `smallest_magnitude` | 絶対値の小さいもの | 小さい順 |
| `largest_algebraic` | 値(実部)の大きいもの | 大きい順 |
| `smallest_algebraic`(既定) | 値(実部)の小さいもの | 小さい順 |
| `target_magnitude` | 絶対値が `shift` の絶対値に近いもの | 近い順 |
| `target_real` | 実部が `shift` に近いもの | 近い順 |

`use_shift = true` のときは、`target` によらず、`shift` に近い固有値を
求めます。「返る順」は `spmats` の場合です。ほかの policy では、同じ固有値が
違う順で返ることがあります。

`method`(`vcp::eig_solver_method`)の `auto_select` は、次の順に解法を
選びます。

| 条件 | 解法 |
| --- | --- |
| `use_shift == true` | shift-invert(対称なら Lanczos 法、ほかは Arnoldi 法) |
| n ≤ 300 で密行列への変換が許可されている(`allow_dense_conversion` が `true` で、n² ≤ `max_dense_size`) | 密行列の解法(対称: Jacobi 法、非対称: 実 Schur 分解) |
| 確実に対称と判定できる | thick-restart Lanczos 法 |
| それ以外 | Krylov–Schur 法 |

`method` に `lanczos`、`arnoldi`、`thick_restart_lanczos`、`krylov_schur`、
`shift_invert_lanczos`、`shift_invert_arnoldi`、`dense_fallback_explicit` を
明示することもできます(`arnoldi` は Krylov–Schur 法の別名です)。対称か
どうかは `structure` で指定でき、`auto_detect` のときは行列を調べて決めます。

一般化固有値問題では、`method = auto_select` は `lanczos` と同じ扱いです。
解ける指定は次の 3 通りで、それ以外は計算せずに空の結果
(`status = "unsupported_generalized_non_diagonal"`)を返します。

| 指定 | 解法 |
| --- | --- |
| A と B がどちらも対角行列 | 対角成分の比 |
| `structure = symmetric`(`method` は既定または `lanczos`) | B 内積の Lanczos 法 |
| `use_shift = true`、または `method` が `shift_invert_lanczos` / `shift_invert_arnoldi` | shift-invert |

つまり、既定のオプションでは、B が対角でない(または A が対角でない)場合に
空の結果を返します。`structure = symmetric` を指定するか、shift-invert を
使ってください。

### 戻り値と status

`vcp::eig_result<T>` の主な欄:

| 欄 | 意味 |
| --- | --- |
| `eigenvalues` | 固有値(複素共役の対を返すときは実部) |
| `eigenvalues_imag` | 固有値の虚部。空なら、全部が実数 |
| `eigenvectors` | 固有ベクトル。`eigenvectors[j]` が `eigenvalues[j]` の固有ベクトル |
| `converged` | 返した固有対がすべて残差の判定を満たしたか |
| `status` | 結果の種類を表す文字列(下の表) |
| `message` / `failure_reason` | 失敗の理由の説明 |
| `requested_count` | 求めた個数 k |
| `returned_count` | 返した個数。どの status でも `eigenvalues.size()` と同じです |
| `converged_count` | 解法が収束と判定した個数(密行列の解法は全部の固有値を解くので、`returned_count` より大きいことがあります) |
| `used_method` | 実際に使った解法の名前(`auto_select` で選ばれた解法は `auto_select(...)` の形) |
| `residuals_absolute` / `residuals_relative` | 固有対ごとの残差 |
| `residual_norm_absolute` / `residual_norm_relative` | 残差の最大値 |
| `iterations` | 反復の回数(密行列の解法では回転の回数) |
| `matrix_vector_products` | 行列とベクトルの積の回数 |
| `linear_solves` | 連立一次方程式を解いた回数(shift-invert) |

診断の欄:

| 欄 | 意味 |
| --- | --- |
| `complex_pair_count` | 返した複素共役の対の数 |
| `complex_eigenvalues` | 複素数の固有値(複素数の要素型のとき) |
| `returned_real_count` / `returned_complex_count` | 返した個数の内訳。`returned_complex_count` は、`eigenvalues_imag` が零でない要素の数(複素共役の対 1 つにつき 2)、`returned_real_count` は残りです。2 つの和は `returned_count` です |
| `residual_history_absolute` / `residual_history_relative` | 残差の履歴(`compute_residual_history = true` のとき) |
| `method` | 指定した解法 |
| `breakdown_reason` | 反復が破綻した理由 |
| `used_orthogonalization` | 使った直交化の名前 |
| `used_dense_fallback` / `used_shift_invert` / `used_generalized_operator` | 密行列の解法・shift-invert・一般化の経路を使ったか |
| `used_subspace_dim` | 使った部分空間の次元 |
| `inner_iterations` / `inner_failure_count` / `inner_residual_norm` / `inner_failure_reason` | shift-invert の内側の解法の記録 |
| `factorization_diagnostics` / `factorization_zero_pivots` | shift-invert の分解の記録 |
| `lambda_gate_products` / `polish_factorizations` / `confirm_restarts` / `ks_rescue_products` | 解法の内部の作業量の内訳 |
| `index_lo` / `index_hi` | 固有値の番号の範囲(`spimats` だけ。ほかの policy では空) |
| `rank_certified` | 求めた順位のものであることの検査の結果(`spimats` だけ。ほかの policy では 0) |

`status` の主な値:

| `status` | 意味 |
| --- | --- |
| `"converged"` | 収束した |
| `"empty_success"` | k = 0 を指定した(空の結果、`converged = true`) |
| `"max_iter_exhausted"` | `max_iter` に達した |
| `"not_converged"` | 収束しなかった |
| `"factorization_failed"` | shift-invert の分解に失敗した(σ が固有値に一致する場合など) |
| `"inner_solve_failed"` | shift-invert の内側の反復が収束しなかった |
| `"unsupported_generalized_non_diagonal"` | 一般化固有値問題で、解けない指定だった |
| `"unsupported_preconditioner_for_method"` | shift-invert でない指定に前処理 M を渡した |
| `"complex_pair_refused"` | 複素共役の対が現れたが、`allow_complex_pairs` が `false`(`spumar`) |
| `"internal_error"` | 内部エラー |

収束しなかったときの注意:

- `converged == false` のとき、`eigenvalues` は空のことも、収束していない
  近似値が入っていることもあります。解法によって違います(thick-restart
  Lanczos 法、Lanczos 法、密行列の解法は近似値を残し、Krylov–Schur 法と
  `spumar` の ARPACK は空にします)。`converged == false` の値は使わないで
  ください。
- 自由関数の `eigs` とメンバの `A.eigs` は、`converged == false` のとき
  例外を投げます。
- `converged == true` は、返した各固有対が残差の判定を満たすことを意味します。
  Krylov 部分空間法は開始ベクトルに直交する固有空間を見られないため、
  「求める側の k 個をすべて見つけた」ことは保証しません。

k と n:

- k = 0 は空の結果を返します(`converged = true`)。
- k が n より大きいときは、n 個を返します(shift-invert では、このとき
  `converged = false` になります)。

複素共役の対(実行列・非対称):

- `allow_complex_pairs = false`(既定)では、複素共役の対は返しません。求める
  範囲に対があると `converged = false` になり、実数の固有値だけが
  `eigenvalues` に残ります。
- `allow_complex_pairs = true` にすると、複素共役の対を
  `eigenvalues`(実部)と `eigenvalues_imag`(虚部)の隣り合う 2 つとして
  返します。固有ベクトル x = u + i v の実部 u と虚部 v は
  `eigenvectors` の隣り合う 2 本に入ります。
- 固有値だけを返す `A.eigs` と自由関数の `eigs` は、虚部を返せません。
  `allow_complex_pairs = true` を指定してこれらを呼ぶと、
  `vcp::invalid_argument` を投げます。複素共役の対を扱うときは
  `A.eigs_with_info` を使ってください。

例外になる場合(`eigs_with_info` でも投げます):

| 状況 | 例外 |
| --- | --- |
| 行列が正方でない、B の大きさが合わない | `vcp::dimension_error` |
| `max_iter = 0`、`tol ≤ 0` | `vcp::invalid_argument` |
| `A.eigs`・自由関数の `eigs` に `allow_complex_pairs = true` を指定した | `vcp::invalid_argument` |

### policy ごとの違い

| 項目 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| `eigs` / `eigs_with_info`(標準) | VCP の実装 | ARPACK に委譲(下の表) | `spumar` と同じ | 保証つき。対称行列の固有値の包含 |
| `eigs` / `eigs_with_info`(一般化) | VCP の実装(解ける指定は「オプション」の表) | ARPACK に委譲(下の表) | `spumar` と同じ | 保証つき。対称行列の固有値の包含 |
| 前処理 M つき | shift-invert の内側の GMRES に M を使う | shift-invert のとき、ARPACK と M つき GMRES で解く | `spumar` と同じ | 保証つき(M は近似の段にだけ渡る) |
| 固有ベクトル | 返す | 返す | 返す | 返さない(`eigenvectors` は空、V は n×0) |
| 非対称の行列 | 解ける | 解ける | 解ける | `vcp::domain_error` |
| 収束・保証に失敗したとき | `converged = false` | `converged = false`(`eigenvalues` は空) | `spumar` と同じ | 例外(`_with_info` でも) |
| 成功したときの `status` | `"converged"` | `"converged"` | `"converged"` | 空の文字列 |
| `used_method` | 解法の名前 | `spumar/arpack(dsaupd)` など | `spumar` と同じ | `spimats_veigs` |

`spumar`・`spcmodumar` の委譲先(`method` が `dense_fallback_explicit` でなく、
k ≤ n − 2 のとき):

| 問題 | 指定 | 委譲先 |
| --- | --- | --- |
| 標準 | `use_shift = true` | ARPACK(shift-invert。A − σI を UMFPACK で分解) |
| 標準 | `largest_magnitude` / `largest_algebraic` / `smallest_algebraic` | ARPACK(A を掛けるだけ) |
| 標準 | `smallest_magnitude`(シフトなし) | ARPACK(σ = 0 の shift-invert) |
| 標準 | `target_real`(シフトなし) | ARPACK(σ = `shift` の shift-invert) |
| 標準 | `target_magnitude`(シフトなし) | 委譲しない(VCP の実装) |
| 一般化 | `use_shift = true` | ARPACK(shift-invert。A − σB を UMFPACK で分解) |
| 一般化 | `smallest_magnitude`(シフトなし) | ARPACK(σ = 0 の shift-invert) |
| 一般化 | `largest_magnitude` / `largest_algebraic` / `smallest_algebraic`(シフトなし) | ARPACK(B を UMFPACK で分解し、B⁻¹A を掛ける) |
| 一般化 | `target_magnitude` / `target_real`(シフトなし) | 委譲しない(VCP の実装) |

- `method = dense_fallback_explicit` のとき、k > n − 2 のとき、k = 0 のときは、
  委譲せずに VCP の実装で計算します。
- 一般化固有値問題は、既定のオプション(`smallest_algebraic`)でも ARPACK に
  委譲されるので、`spmats` と違って B が対角でなくても解けます。
- 複素共役の対を返せるのは標準の固有値問題だけです。一般化では返しません。

`eig_options` の欄の効き方:

| 欄 | `spmats` | `spumar`・`spcmodumar` | `spimats` |
| --- | --- | --- | --- |
| `method` | 効く(解法の選択) | `dense_fallback_explicit` だけ効く(VCP の実装に切り替わる)。ほかの値は無視される | 近似の段に渡る |
| `structure` | 効く | 効く。`symmetric` / `hermitian` なら対称用、`general` なら非対称用の ARPACK のルーチンを使う。`auto_detect` は `is_symmetric()` で決める | 近似の段に渡る。行列は厳密に対称でなければならない |
| `target` | 効く | 効く(上の表) | 効く |
| `tol` | 効く | ARPACK に渡る | 近似の段に渡る |
| `shift` / `use_shift` | 効く | 効く | 近似の段に渡る(下の注意) |
| `max_iter` | 効く | 効く(ARPACK の再開の回数に換算する。少し超えることがある) | 近似の段に渡る |
| `subspace_dim` | 効く | 無視される(ARPACK の次元は `spumar_options` の `ncv`) | 近似の段に渡る |
| `allow_dense_conversion` / `max_dense_size` | 効く | 無視される(n が小さくても ARPACK) | 近似の段に渡る |
| `orthogonalization` / `random_seed` / `random_start` / `compute_residual_history` | 効く | 無視される | 近似の段に渡る |
| `shift_invert_solver` / `shift_invert_lu` | 効く | 無視される(常に UMFPACK) | 近似の段に渡る |
| `allow_complex_pairs` | 効く | 標準の固有値問題で効く | 意味を持たない(対称行列だけを扱う) |

近似計算の policy(`spmats`・`spumar`・`spcmodumar`)に共通の注意:

- `structure = symmetric` の指定は、行列が対称かどうかを検査せずに使われます。
  非対称の行列に指定すると、エラーにならずに、誤った値が `converged = true` で
  返ることがあります。対称かどうかが分からないときは `auto_detect`(既定)の
  ままにしてください。

`spumar`・`spcmodumar` の注意:

- `smallest_algebraic`(既定)や `largest_algebraic` は、固有値が密集している
  側では ARPACK の収束が遅く、既定の `max_iter` では収束しないことがあります。
  零に近い固有値がほしいときは、`smallest_magnitude` を指定するか、
  `use_shift = true` にしてください(shift-invert になり、速く収束します)。

`spimats` の注意(詳しくは「`spimats` の詳細」の節):

- 返る `eigenvalues` は区間で、それぞれが固有値を含みます。
- 返る個数が k より少ないことがあります(固有値を分離できなかった場合)。
  `A.eigs` と自由関数の `eigs` はそのとき例外を投げ、`A.eigs_with_info` は
  少ない個数を返します。
- シフトに近い固有値がほしいときは、`target = target_real`(または
  `target_magnitude`)と `shift` を指定してください。`spimats` では、どの
  固有値を保証するかは `target` で決まり、`use_shift` と `shift` は近似の段への
  手がかりにだけ使われます(近似の段の shift-invert は、`target` によらず
  `shift` に近い固有値を求めます)。そのため、`use_shift = true` だけを指定して
  `target` を既定(最小側)のままにすると、近似の段が最小の固有値を求めず、
  要求より少ない個数の、幅の広い区間が返ることがあります(返る区間は正しい
  包含です)。
- 全部の固有値の包含がほしいときは、k に n を指定します。

固有値問題の内部のアルゴリズムは [kernels/tsparse_eigs.md](kernels/tsparse_eigs.md) を参照してください。

## `spimats` の詳細

`vcp::spimats<T>` は区間疎行列の精度保証付き計算の policy です。次の
プログラムは、解の包含と固有値の包含を求めます。

```cpp
#include <iostream>
#include <vector>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/spimats.hpp>
#include <vcp/spmatrix.hpp>

int main() {
    typedef kv::interval<double> I;
    const int n = 5;

    vcp::spmatrix<I, vcp::spimats<double> > A;
    A.resize(n, n);
    for (int i = 0; i < n; i++) {
        A.add(i, i, I(2.0));
        if (i + 1 < n) {
            A.add(i, i + 1, I(-1.0));
            A.add(i + 1, i, I(-1.0));
        }
    }
    A.finalize();

    std::vector<I> b(n, I(1.0));
    std::vector<I> x = A.solve(b);    // 解の包含(保証に失敗したら例外)
    std::cout << "x[0] in " << x[0] << "\n";

    std::vector<I> lam = A.eigs(2);   // 小さい側の固有値 2 個の包含
    std::cout << "lambda[0] in " << lam[0] << "\n";

    std::cout << "is_spd = " << A.is_spd() << "\n";
    return 0;
}
```

```bash
g++ -std=c++20 -I/path/to/libs example.cpp
```

(`/path/to/libs` は `vcp/` と `kv/` を置いたディレクトリです。)

テンプレート引数は `vcp::spimats<T, AP = vcp::spmats<T>, DP = vcp::imats<T> >`
です。

- `T`: 点の要素型(`double`、`kv::dd`、`kv::mpfr<N>` など)。
- `AP`: 内部の近似計算に使う疎行列 policy。`AP = vcp::spumar` は
  `T = double` のときだけ使えます。
- `DP`: 内部で使う区間密行列 policy。`DP` だけを変える場合は
  `vcp::spimats_dp<double, vcp::pidblas>` のように書けます(`vcp/pidblas.hpp`
  の include が必要です)。

`spimats` で精度保証が付くのは次の計算です。

| 計算 | 返す値 |
| --- | --- |
| `solve` / `solve_with_info` / `lss` / `lss_with_info`(短縮形を含む) | 解の包含 |
| `eigs` / `eigs_with_info`(一般化を含む) | 対称行列の固有値の包含(固有ベクトルは返しません) |
| `is_spd()` | 区間内のすべての対称行列が正定値であることを証明できたら `true` |
| 四則・ベクトルとの積などの演算、`normone` / `norminf` | 真の値を含む区間 |

`is_spd()` の `false` は「正定値でない」という主張ではなく、十分条件が
成り立たなかったことを意味します。

`AP` と `DP` の選択は、包含の成否・包含の幅・計算速度にだけ影響し、
精度保証の正しさには影響しません。

区間の幅を持つ行列(`I(lo, hi)` を要素に入れた行列)では、区間に含まれる
すべての行列についての包含を返します。

### `spimats` の例外の扱い

例外の扱いは他の policy と異なります。保証に失敗したときは、
`solve_with_info` などの `_with_info` 系の関数でも例外を投げます。

| 状況 | 例外 |
| --- | --- |
| 次元の不一致などの誤用 | `vcp::dimension_error`、`vcp::invalid_argument` など |
| 固有値の計算で、行列が区間行列として厳密に対称でない | `vcp::domain_error` |
| 近似計算の失敗(中点の行列が特異で LU 分解できない、など) | `vcp::state_error` |
| 精度保証の十分条件が成り立たない(区間が特異な行列を含む、条件が悪い、など) | `vcp::verification_error` |

成功したときは `converged` は常に `true` です。結果のうち保証が付いているのは、
`solve` 系の `x` と固有値問題の `eigenvalues` だけです。

### `spimats` の線形方程式

- `linear_solve_options` の `method` などは使われません。どの指定でも、
  同じ保証つきの手順で解きます(`max_iter` と `tol` は、正の値かどうかだけを
  検査します)。`solve_cg` などの短縮形も同じです。
- 結果の `x` は区間で、真の解を含みます。`iterations` は 0、`method` は
  指定した値のままです。残差の欄は使いません。
- 複数の右辺を取る `lss` / `lss_with_info` は、分解を 1 回だけ行います。
  どれかの列で保証に失敗すると、例外を投げます。

### `spimats` の固有値問題

- 行列は、区間行列として厳密に対称(対になる要素の区間の両端が一致)で
  なければなりません。一般化固有値問題 A x = λ B x では、B が正定値である
  ことの証明もできなければなりません(できないときは
  `vcp::verification_error`)。
- 返る `eigenvalues` は区間で、求めた順位の固有値を含みます。固有ベクトルは
  返しません。
- `index_lo[j]`〜`index_hi[j]` は、`eigenvalues[j]` が含む固有値の番号
  (全部の固有値を小さい順に並べたときの番号。1 始まり)です。
- `rank_certified` は、返した固有値が求めた順位(最小、最大、`shift` に最も
  近い、など)のものであることを証明できたかどうかです(1: 証明できた、
  −1: 証明できなかった、0: 検査していない)。
- 固有値を分離できないとき(区間の幅が広い行列、近接した固有値)は、返る
  個数が k より少なくなることがあります。このとき `converged` は `true` の
  ままです。`A.eigs` は個数が足りないと `vcp::state_error` を投げます。
- k = 0 は `vcp::invalid_argument` です。k が n より大きいときは n 個を
  返します。
- 成功したときの `status` は空の文字列、`used_method` は `spimats_veigs` です。
- `eig_options` の欄は、内部の近似計算(`AP`)に渡ります。前処理 M も近似
  計算にだけ使われ、保証には関わりません。

### `spimats` の `is_spd`

`A.is_spd()` は `spmatrix` のメンバとして呼べます。

| 戻り値・例外 | 意味 |
| --- | --- |
| `true` | 区間内のすべての対称行列が正定値であることを証明できた |
| `false` | 証明できなかった(正定値でないとは限らない) |
| `vcp::dimension_error` | 行列が正方でない |
| `vcp::domain_error` | 区間行列として厳密に対称でない(`spd_options` の `nosymcheck = true` で検査を省ける) |

### 保証の主張がない計算

`lu`・`lu_factorize_with_info`・`lu_solve`・`ldl`・`inertia`・
`ldl_shift_setup_with_info`・`chol`・`ainv_*`・`fsai_*` は、`spimats` では
差し替えられていません。区間を要素にした `spmats` の実装がそのまま動きます
(区間演算で実行)。これらの結果には保証の主張がなく、この使い方は非推奨
です。

- 幅のない区間だけからなる行列では、`spmats` と同じように進み、status は
  `success` になります。
- 幅のある区間を含む行列では、`ldl`・`chol` は対称性の検査が通らず
  `not_symmetric` になり、`inertia` は `factorization_failed` になります。
- 連立一次方程式の解の包含がほしいときは `solve` を、正定値性の証明が
  ほしいときは `is_spd` を、固有値の包含がほしいときは `eigs` を使って
  ください。

`is_symmetric()` は、区間の両端が一致するかどうかで厳密に判定します。
`is_symmetric(tol)` は、`spimats` ではコンパイルエラーになります。

ノルムの上界には `norminf` / `normone` を使います。区間型では
`B.get(0, 0).upper()` が保証付きの上界です(`B` は `norminf(A)` などの
結果)。`normtwo` は近似計算であり、保証には使えません。

### policy のオブジェクトを直接使う

保証層のオプション `spi_opt`、計算の記録 `last_verify_info` は
`spimats` 固有のメンバです。`spmatrix` は policy を protected 継承している
ため、これらは policy のオブジェクトを直接使って呼び出します。

```cpp
vcp::spimats<double> A;
A.resize(n, n);
A.add(i, j, I(lo, hi));
A.finalize();
bool ok = A.is_spd();
int stage = A.last_verify_info.stage;
```

| メンバ | 内容 |
| --- | --- |
| `A.spi_opt` | 保証層のオプション(`vcp::spimats_options`。下の表) |
| `A.last_verify_info` | 直前の計算の記録(`vcp::spimats_verify_info<T>`) |
| `A.is_spd(opt)` | 行列自身が対象(`spmatrix` の `is_spd` と同じ判定) |
| `P.is_spd(Ap, opt)` | 点の疎行列 `Ap`(`vcp::spmats<T>`)が対象。`P` は `vcp::spimats<T>` のオブジェクト |

`vcp::spimats_options`:

| 欄 | 既定値 | 意味 |
| --- | --- | --- |
| `keep_ap_eigenvectors` | `false` | 近似固有ベクトルを `last_verify_info` に残す |
| `target_rank_check` | `record` | 返した固有値が求めた順位のものかを検査するか。`off`(検査しない)/ `record`(検査して結果を記録する)/ `require`(証明できなければ `vcp::verification_error`) |

`last_verify_info` には近似計算の結果などが点型のまま記録されます。これは
計算の記録であり、対象の行列についての主張ではありません。主な欄は、
進行の段階 `stage`(4 が成功)、近似解 `ap_lss_result`、近似固有値
`ap_eigenvalues` です。

精度保証の手順と混合演算は [kernels/spimats_kernel.md](kernels/spimats_kernel.md) を参照してください。

## `spumar` の詳細

`vcp::spumar` は `spmats<double>` を継承し、疎 LU を UMFPACK に、疎固有値を
ARPACK に委譲する policy です。

```cpp
#include <iostream>
#include <vector>

#include <vcp/spumar.hpp>
#include <vcp/spmatrix.hpp>

int main() {
    const int n = 50;

    vcp::spmatrix<double, vcp::spumar> A;
    A.resize(n, n);
    for (int i = 0; i < n; i++) {
        A.add(i, i, 2.0);
        if (i + 1 < n) {
            A.add(i, i + 1, -1.0);
            A.add(i + 1, i, -1.0);
        }
    }
    A.finalize();

    std::vector<double> b(n, 1.0);
    auto r = A.solve_with_info(b);          // UMFPACK
    if (!r.converged) return 1;
    std::cout << "x[0] = " << r.x[0] << "\n";

    vcp::eig_options<double> eopt;
    eopt.target = vcp::eig_target::smallest_magnitude;
    auto er = A.eigs_with_info(4, eopt);    // ARPACK
    if (!er.converged) return 1;
    std::cout << "lambda[0] = " << er.eigenvalues[0] << "\n";
    return 0;
}
```

インストールとリンクの例(Ubuntu):

```bash
sudo apt install libsuitesparse-dev libarpack2-dev
g++ -std=c++20 -I/path/to/libs example.cpp -lumfpack -larpack
```

| 操作 | 委譲先 |
| --- | --- |
| `solve` / `solve_with_info` / `lss` / `lss_with_info`(短縮形を含む) | UMFPACK(`method` の指定によらない) |
| `lu` / `lu_with_info` | UMFPACK |
| `eigs` / `eigs_with_info`(標準・一般化) | ARPACK。委譲されない指定は「固有値問題」の節の「policy ごとの違い」の表 |
| `lu_factorize_with_info`・LU 因子の消費(`lu_solve` など)・LDLᵀ 分解・慣性・Cholesky 分解・`is_spd`・近似逆行列・行列演算 | 委譲しない(`spmats` の実装) |

### `spumar` の線形方程式と LU 分解

- `solve` 系は、`linear_solve_options` の `method` が何であっても UMFPACK の
  直接法で解きます。`converged` は UMFPACK の分解と求解が成功したかどうかで
  決まり、`tol` は合否に使いません。結果の `method` は指定した値のまま
  (既定なら `auto_select`)で、`iterations` は UMFPACK が行った反復改良の
  回数です。
- 特異な行列では、`solve_with_info` は `converged = false` と零ベクトルを
  返し、`solve` は `vcp::state_error` を投げます。
- `linear_solve_options` の `sparse_lu` と、`lu_extract_options` の `slu` の
  `method` / `ordering` / `pivoting` などは、UMFPACK には対応する設定がない
  ため適用されません(`slu.equilibration = true` だけは `spmats` と同じく
  `unsupported_options` になります)。`lu_with_info` の結果の `method_used` /
  `ordering_used` には、指定した値がそのまま入ります。
- `lu` / `lu_with_info` は、UMFPACK のスケーリングを切って分解するので、
  返った因子で `P A Q = L U` が成り立ちます。

### `spumar` のオプション

オプションは `vcp::spumar_options` で、policy の公開メンバ `muar_opt` に
あります。

| フィールド | 既定 | 意味 |
| --- | --- | --- |
| `ncv` | 0(自動) | ARPACK の部分空間の次元。0 のときは min(n, max(2k + 1, 20))、指定したときは min(n, max(`ncv`, k + 2)) |
| `umf_irstep` | −1(UMFPACK の既定) | UMFPACK の反復改良回数。`solve` 系にだけ効く |
| `umf_pivot_tol` | −1.0(UMFPACK の既定) | UMFPACK のピボット閾値 |

`spmatrix` 経由で使う場合は、`vcp::spumar::set_default_options(o)` で
プロセス全体の既定値を設定します。この関数はスレッド安全ではありません。

| 呼び方 | 内容 |
| --- | --- |
| `vcp::spumar::set_default_options(o)` | これから作る行列の既定値を `o` にする |
| `vcp::spumar::get_default_options()` | 今の既定値を返す |
| `A.muar_opt` | policy のオブジェクト(`vcp::spumar A;`)のオプション。直接書き換えられる |
| `A.last_delegate_info` | policy のオブジェクトの、直前の委譲の記録(`backend`・`note` と、UMFPACK・ARPACK が返した状態) |

```cpp
vcp::spumar_options o;
o.ncv = 30;
vcp::spumar::set_default_options(o);

vcp::spmatrix<double, vcp::spumar> A;   // この行列から ncv = 30 が使われる
```

既定値は、行列を作るときに写し取られます。`set_default_options` を呼ぶ前に
作った行列には効きません。行列を複製すると、オプションも複製されます。

`last_delegate_info.backend` は、直前の呼び出しがどこへ行ったかを表します
(`"umfpack"`、`"arpack"`、VCP の実装に戻ったときは
`"spmats-base-fallback"`)。`muar_opt` と `last_delegate_info` は `spmatrix`
からは見えないので、使うときは policy のオブジェクト `vcp::spumar` を直接
使ってください。

UMFPACK・ARPACK への委譲の仕組みは [kernels/sparse_external.md](kernels/sparse_external.md) を参照してください。

## `spcmodumar` の詳細

`vcp::spcmodumar` は `vcp::spumar` を継承し、さらに疎 Cholesky 分解
(`chol` / `chol_with_info`)を CHOLMOD に委譲する policy です。`spumar` の
委譲(UMFPACK・ARPACK)はそのまま引き継ぎます。

```cpp
#include <iostream>
#include <vector>

#include <vcp/spcmodumar.hpp>
#include <vcp/spmatrix.hpp>

int main() {
    const int n = 50;

    vcp::spmatrix<double, vcp::spcmodumar> A, L;
    A.resize(n, n);
    for (int i = 0; i < n; i++) {
        A.add(i, i, 2.0);
        if (i + 1 < n) {
            A.add(i, i + 1, -1.0);
            A.add(i + 1, i, -1.0);
        }
    }
    A.finalize();

    std::vector<int> p;
    A.chol(L, p);                       // CHOLMOD(失敗したら例外)
    std::cout << "nnz(L) = " << L.nnz() << "\n";

    std::vector<double> b(n, 1.0);
    auto s = A.solve_with_info(b);      // UMFPACK
    if (!s.converged) return 1;
    std::cout << "x[0] = " << s.x[0] << "\n";
    return 0;
}
```

```bash
sudo apt install libsuitesparse-dev libarpack2-dev
g++ -std=c++20 -I/path/to/libs example.cpp -lumfpack -larpack -lcholmod
```

リンクには `-lumfpack -larpack -lcholmod` が必要です(`spumar` を継承する
ため、`-lcholmod` だけでは足りません)。

| 操作 | 委譲先 |
| --- | --- |
| `chol` / `chol_with_info` / `Cholesky` / `is_spd` | CHOLMOD |
| そのほか | `spumar` と同じ |

分解の規約(`Pᵀ A P = L Lᵀ`、p は new→old)は `spmats` と同じです。
ただし、CHOLMOD は順序付けに後処理の並べ替えを合成するため、返る p は
`spmats` と一致するとは限りません。返った p で恒等式は常に成り立ちます。

| `ordering` | CHOLMOD の設定 |
| --- | --- |
| `natural` | 自然順序 |
| `auto_select` / `amd` | AMD |
| `rcm` | VCP の RCM を与える |
| `nested_dissection` | NESDIS |
| `nested_dissection_ml` | VCP の多段 nested dissection の順序を与える |

- `auto_select` の解決先は `amd` です(`spmats` では
  `nested_dissection_ml`)。`ordering_used` に記録されます。
- `pd_tol` に 0 以外を指定すると `invalid_options` になります。
- `method` に指定できるのは `auto_select`(CHOLMOD が選ぶ)と
  `simplicial_uplooking` です。`supernodal` を指定すると `invalid_options` に
  なります。
- `method_used` に記録される `supernodal` は `spcmodumar` だけの値です。
- supernodal 経路の `nnz_L` は明示的に格納した零を含むため、simplicial
  経路の値と直接比べられません。
- `check_symmetry` と `symmetry_tol` は `spmats` と同じように効きます。
  `ndml_*` の欄は、`ordering = nested_dissection_ml` のときに `spmats` と同じ
  ように効きます。
- `is_spd` は CHOLMOD の分解が成功したかどうかを返します。

CHOLMOD への委譲の仕組みは [kernels/sparse_external.md](kernels/sparse_external.md) を参照してください。

## 外部ライブラリのライセンスに関する注意

VCP Library 本体と、`spumar` / `spcmodumar` のヘッダ(`vcp/spumar.hpp`、
`vcp/spumar_base/`、`vcp/spcmodumar.hpp`、`vcp/spcmodumar_base/`)は
BSD 3-Clause License です。UMFPACK・ARPACK・CHOLMOD のコードは同梱せず、
API を呼び出すだけです。

| ライブラリ | ライセンス |
| --- | --- |
| UMFPACK + AMD(SuiteSparse) | LGPL-2.1+ |
| ARPACK-NG | BSD |
| CHOLMOD(SuiteSparse) | モジュールごとに異なる。Check / Cholesky / Core(Utility) / Partition は LGPL、Supernodal / MatrixOps / Modify は GPL |

- `spumar.hpp` も `spcmodumar.hpp` も include・リンクしないコードは、
  LGPL・GPL の影響を受けません。`spumar` だけを使うコードは CHOLMOD の
  GPL の影響を受けません。
- 通常のフルビルドの CHOLMOD をリンクした結合著作物は GPL の条件に従います。
  GPL のモジュールを除いたビルドなら LGPL ですが、simplicial 経路だけに
  なります。どちらを選ぶかは利用者の判断です。
- 共有ライブラリへの動的リンクか、静的リンク・改変して同梱するかで義務が
  異なります。配布する場合は各ライブラリのライセンス原文を確認して
  ください。
- この注意は情報提供であり、法的助言ではありません。
