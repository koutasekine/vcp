# tsparse — テンプレート版の疎行列カーネル

`tsparse` は、疎行列 policy `vcp::spmats<T>` が下請けとして使う
header-only の疎行列カーネル群です(`vcp/tsparse/`)。要素型 `T` の
template として書かれており、`double`、`kv::dd`、`kv::mpfr<N>`、
`kv::interval<U>` などで使えます。利用者が直接呼ぶことは想定していません。
使い方は [spmatrix.md](../spmatrix.md)、policy との関係は
[sparse_layers.md](sparse_layers.md) を参照してください。

この文書では、格納形式、形式の変換、疎行列とベクトルの積、疎行列どうしの
和と積、要素型の補助を説明します。分解・固有値解法・反復法は、それぞれ
次の文書で説明します。

| 文書 | 内容 |
| --- | --- |
| [tsparse_lu.md](tsparse_lu.md) | 疎 LU 分解 |
| [tsparse_ldl.md](tsparse_ldl.md) | 疎 LDLᵀ 分解・Cholesky 分解・慣性 |
| [tsparse_eigs.md](tsparse_eigs.md) | 疎行列の固有値解法 |
| [tsparse_iterative.md](tsparse_iterative.md) | 作用素・前処理・反復法 |

## file 構成

`tsparse.hpp` は基本の header をまとめて include します。LDLᵀ 分解、
Cholesky 分解、LU 因子の抽出、逆行列の行・列の計算、thick-restart Lanczos
法、Krylov–Schur 法は `tsparse.hpp` に含まれないため、個別に include します。
`detail/` 以下の header は内部実装です。

この文書で扱う header は次の通りです。

| file | 内容 |
| --- | --- |
| `tsparse.hpp` | 基本 header の umbrella header |
| `tsparse_format.hpp` | 格納形式の列挙型 `sparse_format` と、生配列を参照する `coo_view` / `csr_view` / `csc_view` |
| `tsparse_convert.hpp` | COO の整列・重複の加算・零の除去、COO / CSR / CSC の相互変換 |
| `tsparse_blas.hpp` | 疎行列とベクトルの積 `tcsrmv` / `tcscmv` と、その並列版 |
| `tsparse_spgemm.hpp` | 疎行列どうしの和・積(名前空間 `tsparse_spgemm`) |
| `tsparse_scalar.hpp` | 要素型の補助(名前空間 `tsparse_scalar`) |

## 格納形式

すべて 0 始まりの添字です。

| 形式 | 配列 |
| --- | --- |
| COO | `row[k]`、`col[k]`、`value[k]`(k = 0, …, nnz−1) |
| CSR | `row_ptr`(長さ rows+1)、`col_ind`、`value`。第 i 行は `row_ptr[i]` から `row_ptr[i+1]` の手前まで |
| CSC | `col_ptr`(長さ cols+1)、`row_ind`、`value`。第 j 列は `col_ptr[j]` から `col_ptr[j+1]` の手前まで |

`coo_view` / `csr_view` / `csc_view` は、行数・列数・格納数と、これらの
配列への const pointer をまとめた構造体です。配列を所有しません。

### `spmats` の格納

`spmats` は、組み立て中は COO 形式で保持し、`finalize()`(`to_csr()`)で
CSR 形式に、`to_csc()` で CSC 形式に変換します。変換では次の順に COO を
正規化します。

1. `tcoo_sort` で (行, 列) の順に整列する。`set` で入れた要素を区別する
   印を付けた版を使い、同じ位置の要素は入れた順を保つ。
2. `tcoo_sum_duplicates_tagged` で同じ位置の要素をまとめる。`add` の値は
   加算し、`set` の値はそれまでの値を捨てて置き換える。
3. `tcoo_remove_zeros` で、値が厳密に零の要素を除く。
4. `tcoo_to_csr` / `tcoo_to_csc` で CSR / CSC に変換する。

その結果、変換後の各行(列)の添字は昇順で重複がなく、厳密零の値を含み
ません。

## 形式の変換

```cpp
#include <vcp/tsparse/tsparse_convert.hpp>

vcp::tcoo_sort(nnz, row, col, value);
nnz = vcp::tcoo_sum_duplicates(nnz, row, col, value);
nnz = vcp::tcoo_remove_zeros(nnz, row, col, value);
vcp::tcoo_to_csr(rows, cols, nnz, row, col, value, row_ptr, col_ind, val);
vcp::tcoo_to_csc(rows, cols, nnz, row, col, value, col_ptr, row_ind, val);
vcp::tcsr_to_coo(rows, cols, nnz, row_ptr, col_ind, val, row, col, value);
vcp::tcsr_to_csc(rows, cols, nnz, row_ptr, col_ind, val, col_ptr, row_ind, val2);
```

| 関数 | 動作 |
| --- | --- |
| `tcoo_sort` | (行, 列) の辞書式順に整列する(配列をその場で並べ替える)。印の配列 `is_set` を渡す版は、同じ位置の要素の順を保つ安定な整列 |
| `tcoo_sum_duplicates` | 整列済みの COO で、同じ位置の要素の値を加算して 1 つにまとめ、新しい格納数を返す |
| `tcoo_sum_duplicates_tagged` | 上と同じだが、印の付いた要素(`set`)では、それまでの和を捨てて値を置き換える |
| `tcoo_remove_zeros` | `value[k] == T(0)` の要素を除き、新しい格納数を返す |
| `tcoo_to_csr` / `tcoo_to_csc` | 行(列)ごとの個数を数えて累積和で pointer を作り、要素を振り分ける。同じ行(列)の中の順は COO の順のまま |
| `tcsr_to_coo` | CSR を COO に展開する |
| `tcsr_to_csc` | COO を経由して CSC に変換する |

`tcoo_to_csr` / `tcoo_to_csc` は重複をまとめないので、重複のない整列済みの
COO を渡します。範囲外の添字は `vcp::index_error`、負のサイズは
`vcp::invalid_argument` を投げます。

区間型では、`==` は kv の「確実に等しい」の意味なので、
`tcoo_remove_zeros` が除くのは点区間 [0, 0] だけです。0 を含む幅のある
区間は残ります。

## 疎行列とベクトルの積

```cpp
#include <vcp/tsparse/tsparse_blas.hpp>

// y := alpha*op(A)*x + beta*y、op(A) = A('N') または Aᵀ('T')
vcp::tcsrmv(trans, rows, cols, alpha, row_ptr, col_ind, val, x, beta, y);
vcp::tcscmv(trans, rows, cols, alpha, col_ptr, row_ind, val, x, beta, y);
```

`trans` は `'N'` または `'T'` です(小文字も可)。最初に `y` の全要素に
`beta` を掛けるので、`beta = 0` でも `y` は初期化しておく必要があります。

格納形式と `trans` の組み合わせによって、計算の形が 2 通りあります。

| 組み合わせ | 計算の形 |
| --- | --- |
| CSR で `'N'`、CSC で `'T'` | 出力の 1 要素ごとに、行(列)の格納要素の内積を取る(gather) |
| CSR で `'T'`、CSC で `'N'` | 入力の 1 要素ごとに、行(列)の格納要素を出力へ足し込む(scatter) |

gather の形には並列版 `tcsrmv_gather_par` / `tcscmv_gather_par` が
あります(引数に `trans` はありません)。出力の各要素を逐次版と同じ順の
演算で計算するので、結果はスレッド数によらず逐次版と一致します。並列領域の
中で投げられた例外は、1 つだけ捕まえて並列領域の外で投げ直します。

`spmats` の `mul_vec` / `trans_mul_vec` は、gather の形で格納数が
`VCP_SPMATS_OMP_THRESHOLD` 以上のときに並列版を使い、それ以外は逐次版を
使います。

## 疎行列どうしの和と積

```cpp
#include <vcp/tsparse/tsparse_spgemm.hpp>

// C := alpha*A + beta*B(A, B, C は CSR)
vcp::tsparse_spgemm::csr_csr_linear_combination(rows, cols,
    ao, ai, av, bo, bi, bv, alpha, beta, emit);
// C := A*B
vcp::tsparse_spgemm::csr_csr_multiply(rows, inner_dim, cols,
    ao, ai, av, bo, bi, bv, emit);
```

入力は CSR の 3 つの配列(`std::vector`)です。`emit(i, j, v)` は、結果の
格納要素ごとに呼ばれる関数オブジェクトです。

どちらも、結果を 1 行ずつ作ります。長さ cols の作業配列と、その行で値が
入った列の一覧を持ち、第 i 行について

- 和では、A と B の第 i 行の要素を作業配列に足し込む。
- 積では、A の第 i 行の各要素 a_ik について、B の第 k 行の要素に a_ik を
  掛けて作業配列に足し込む。

その後、列の一覧を昇順に整列し、値が厳密零でない要素を `emit` に渡して、
作業配列を零に戻します。結果は行の昇順、各行の中は列の昇順で、重複と
厳密零を含みません。

並列版 `csr_csr_linear_combination_par` / `csr_csr_multiply_par` は、
`emit` の代わりに結果の CSR 配列 `co`、`ci`、`cv` を返します。結果の
格納数は値によって変わるので(厳密零を除くため)、2 段で作ります。

1. 行を連続したブロックに分け、ブロックごとに逐次版と同じ手順で各行を
   計算し、ブロック内の作業領域に結果と各行の個数を記録する(並列)。
2. 各行の個数の累積和で、結果の pointer を決める(逐次)。
3. ブロックごとに、作業領域の結果を決まった位置に写す(並列)。

各行の演算の順は逐次版と同じなので、結果はスレッド数によらず逐次版と
一致します。

`spmats` の和・差・積は、nnz(A) + nnz(B) が `VCP_SPMATS_OMP_THRESHOLD`
以上のときに並列版を使い、それ以外は逐次版を使います。

## OpenMP

並列版は OpenMP が有効なときだけ並列に動きます。並列化は、kernel の層
(`vcp/tsparse/`)と policy の層(`spmats`・`spmats_base`・`spimats`・
`spimats_base`)で別々に止められます。

| マクロ | 効果 |
| --- | --- |
| `VCP_SPARSE_NOMP` | kernel の層(`tcsrmv_gather_par` などの並列版)の OpenMP を止める |
| `VCP_SPMATS_NOMP` | policy の層の OpenMP を止める(並列版を呼ぶかどうかの判定を含む) |
| `VCP_NOMP` | VCP 全体の OpenMP を止める(上の 2 つを含む) |

policy の層だけを止めると、kernel の並列版は呼ばれません。kernel の層だけを
止めると、policy の層が並列版を呼んでも、並列版の中は逐次に動きます。

並列版を使う閾値 `VCP_SPMATS_OMP_THRESHOLD` の既定値は 4096 です。
コンパイル時に定義して変えられます。

## 要素型の補助 `tsparse_scalar`

`tsparse` と `spmats_base` のコードは、要素型に依存する処理を
名前空間 `vcp::tsparse_scalar` の関数で書きます。

| 名前 | 内容 |
| --- | --- |
| `real_type<T>::type` | 実部の型。`std::complex<U>` では `U`、それ以外は `T` |
| `is_complex<T>` | `std::complex<U>` かどうか |
| `real_part(x)` | 実部 |
| `abs_value(x)`、`sqrt_value(x)`、`hypot_value(x, y)` | ADL で探索した `abs`・`sqrt` による絶対値・平方根、√(x² + y²) |
| `conjugate_if_needed(x)` | 複素数なら共役、それ以外はそのまま |
| `epsilon<T>()` | `std::numeric_limits<T>::epsilon()`。`numeric_limits` が特殊化されていない型(`kv::interval` など)では `T(0)` |
| `is_finite(x)` | 浮動小数点型なら `std::isfinite`、それ以外は常に `true` |
| `decimal_power_negative<T>(e)` | 10⁻ᵉ を `T` で計算する(1 を 10 で e 回割る) |
| `real_dot`、`hermitian_dot`、`real_dot_value`、`hermitian_norm_value`、`real_norm_value` | ベクトルの内積とノルム |

オプションの許容誤差などのフィールドは `real_type<T>::type` で持ちます。

## 要素型の要件

`tsparse` と `spmats_base` のコードが要素型 `T` に要求するのは、
`vcp::mats<T>` と同じ次の演算だけです。値を `double` に変換したり、`T` の
種類で処理を分けたりはしません。

- copy・代入、`T(0)`、`T(1)` などの `int` からの構築
- 四則演算と単項 `-`
- 比較演算
- ADL で探索される `abs`、`log`、`sqrt`、`exp`(`using std::abs; abs(x);` の形)

`float`、`double`、`long double`、`kv::dd`、`kv::mpfr<N>` で使えます。

`kv::interval<U>` でもコンパイルでき、動きます。ただし、判定はすべて
「確実に成り立つと言えるときだけ通す」形で書かれています。例えば
`!(x > tol)` は「x > tol と言えない」の意味です。区間型では
`epsilon<T>()` が `T(0)` になり、0 を含む区間のピボットは割り算の前に
はじかれて、分解は `numerical_singularity` などの status で止まります。
疎 LU 分解では、M 行列のように区間 Gauss 消去が成功する行列で包含が得られます。
これは近似計算の層であり、`kv::interval<U>` で計算しても精度保証付きの
計算にはなりません。精度保証付きの計算には `vcp::spimats<T>` を使って
ください([spimats_kernel.md](spimats_kernel.md))。

## 添字型

添字型 `Index` は整数型です。疎 LU・LDLᵀ・Cholesky 分解は符号付きの添字型
だけを受け付けます(`static_assert` で検査します)。`spmats` の既定は
`int` です。
