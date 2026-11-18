# Sparse Matrix / Policy Guide

この文書では、VCP Library の疎行列クラス `vcp::spmatrix` と policy の関係を
説明します。密行列 `vcp::matrix` については [matrix.md](matrix.md) を参照して
ください。

基本形は次の通りです。

```cpp
vcp::spmatrix<ElementType, Policy>
```

`ElementType` は行列要素の型です。`double`、`kv::dd`、`kv::mpfr<N>`、
`kv::interval<double>` などを指定します。`Policy` は、格納、線形方程式、
分解、固有値計算などをどの実装で行うかを指定します。`Policy` を省略すると
`vcp::spmats<ElementType>` が使われます。

## 必要なヘッダ

区間演算や `kv::dd`・`kv::mpfr<N>` を使う例では、VCP Library とは別に
取得した kv ライブラリが include path から見える必要があります。

通常の疎行列計算(`vcp::spmats<T>`。追加依存なし):

```cpp
#include <vcp/spmatrix.hpp>
```

密行列との相互変換や、密行列を右辺に取る `lss` を使う場合は、
`matrix.hpp` も include します。

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
これらを include しません。

## policy 一覧

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

各 policy の内部で使われる下請けは [kernels/sparse_layers.md](kernels/sparse_layers.md) を参照してください。

## policy の選び方

まず疎行列計算をする場合は `vcp::spmats<T>` を使います。追加依存がなく、
`double` だけでなく `kv::dd`、`kv::mpfr<N>` などでも同じインターフェースを
使えます。

```cpp
vcp::spmatrix<double> A;
```

区間演算による包含が必要な場合は `vcp::spimats<T>` を使います。

```cpp
typedef kv::interval<double> I;
vcp::spmatrix<I, vcp::spimats<double> > A;
```

`double` の大規模な問題で、fill-in の多い直接法や固有値計算が遅い場合は
`vcp::spumar` を使います。対称正定値行列の Cholesky 分解も外部に委譲したい
場合は `vcp::spcmodumar` を使います。どちらも外部ライブラリのライセンスに
注意してください(後述の「外部ライブラリのライセンスに関する注意」)。

```cpp
vcp::spmatrix<double, vcp::spumar> A;
vcp::spmatrix<double, vcp::spcmodumar> B;
```

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
- 値が厳密に零の要素は格納しません。

`operator()` による要素の読み書きもできます。

```cpp
double c = A(1, 1);
A(0, 1) = 7.0;
A(0, 1) += 1.0;
```

よく使う初期化:

```cpp
A.zeros(n, m);    // n x m, 格納要素なし
A.eye(n);         // n x n identity matrix
A.clear();
A.reserve(nnz);   // 格納数の予約
```

サイズと格納:

```cpp
auto r   = A.rowsize();
auto c   = A.columnsize();
auto nnz = A.nnz();

A.is_finalized();
A.format();                         // COO / CSR / CSC
const auto& outer = A.outer_index();
const auto& inner = A.inner_index();
const auto& val   = A.values();
```

CSR / CSC 配列から直接構成することもできます(`assign_csr` /
`assign_csc`)。この場合、各行(列)の添字は昇順・重複なしで、厳密零の値を
含まない必要があります。

添字の型は `spmatrix<T>::index_type` です(既定は `int`)。変更する場合は
`vcp::spmats<T, long long>` のように policy の第 2 テンプレート引数で
指定します。添字型は符号付き整数である必要があります。

部分行列は `initializer_list` で取り出し、ブロックに書き込めます。

```cpp
vcp::spmatrix<double> S = A({0, 1}, {1, 2});
A({0, 1}, {1, 2}) = S;
```

内部の格納形式と演算のアルゴリズムは [kernels/tsparse.md](kernels/tsparse.md) を参照してください。

## 演算

| 種類 | 関数 |
| --- | --- |
| 四則 | `A + B`、`A - B`、`A * B`、`-A`、`A * s`、`s * A`、`A / s`、`+=`、`-=`、`*=` |
| 疎行列とベクトル | `A * x`、`x * A`、`A.mul_vec(x)`、`A.trans_mul_vec(x)` |
| 疎行列と密行列 | `A * B`、`B * A`(`B` は `vcp::matrix`、結果は密行列) |
| 構造 | `transpose(A)`、`diag(A)`、`diag(v)`、`tril(A)`、`triu(A)`、`A.horzcat(B, C)`、`A.vercat(B, C)` |
| 要素ごとの関数 | `abs`、`sqrt`、`sin`、`cos`、`exp`、`log`、`pow` |
| 集約 | `max`、`min`、`sum`、`length` |
| ノルム | `normone`、`norminf`、`normtwo` |
| その他 | `ltransmul(A)`(`Aᵀ A`)、`A.is_symmetric()`、`A.is_symmetric(tol)` |

呼び出し方と結果の形は `vcp::matrix` の同名の関数に合わせています。
集約とノルムは 1×1 などの `spmatrix` を返します。`cos`、`exp`、`log` の
ように `f(0) ≠ 0` となる関数の結果は、構造的にほぼ密になります。

`std::cout << A` は `(i,j) 値` の形で格納要素を出力します。

内部の格納形式と演算のアルゴリズムは [kernels/tsparse.md](kernels/tsparse.md) を参照してください。

## 密行列・他の型との変換

疎行列から密行列へは `to_dense()` で変換します。戻り値の型は
`spmatrix<T>::dense_matrix_type`(`spmats<T>` では
`std::vector<std::vector<T> >`)です。

密行列から疎行列へは代入で変換します。厳密に零の成分は格納されません。

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

## 線形方程式

```cpp
std::vector<double> x = A.solve(b);
auto r = A.solve_with_info(b);   // r.x, r.converged, r.residual_norm など

lss(A, b, x);                    // x = A \ b
vcp::matrix<double, vcp::mats2<double> > X = lss(A, B);   // 密な右辺 B
```

`solve` は失敗時に例外を投げます。`solve_with_info` は失敗しても例外を投げず、
結果の `converged` で成否を返します。

解法は `linear_solve_options<T>` の `method` で選びます。既定の
`auto_select` は疎 LU 直接法(`sparse_lu`)に解決されます。

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

反復法には `solve_cg(b, max_iter, tol)` などの短縮形もあります。

反復法(`conjugate_gradient`・`bicgstab`・`gmres`・`jacobi`・`gauss_seidel`)と
前処理は、まだ実験的な実装です。

直接法の `converged` は、分解と求解が成功したかどうかで決まります。
`residual_norm` は参考値であり、残差が小さいことは解が正確であることを
意味しません。精度が必要な場合は `vcp::spimats<T>` を使ってください。

同じ行列で何度も解く場合は、分解を保持するハンドルを使います。

```cpp
auto h = A.lu_factorize_with_info();
if (h.valid()) {
    auto r1 = h.solve_with_info(b1);
    auto r2 = h.solve_with_info(b2);
}
```

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

`lu` は status が `success` 以外のとき例外を投げます。`lu_with_info` は
status(`vcp::sparse_lu_extract_status`)を返します。L / U / p / q が有効なのは
`success` のときだけです。

| status | 意味 |
| --- | --- |
| `success` | 分解成功 |
| `invalid_factorization` | 分解自体が失敗(特異行列など) |
| `unsupported_options` | スケーリングを伴う指定(`equilibration = true` など) |
| `unsupported_storage` | 未対応の格納形式 |
| `internal_error` | 内部エラー |

オプション `lu_extract_options<T>` は、疎 LU のオプション
`sparse_lu_options<T>` をメンバ `slu` に持ちます。

```cpp
vcp::spmatrix<double>::lu_extract_options_type opt;
opt.slu.ordering = vcp::sparse_lu_ordering::colamd;
```

順序付け(`ordering`)は `auto_select`、`natural`、`rcm`、`amd`、`colamd`、
`nested_dissection`、`nested_dissection_ml` から選べます。

`P A Q = L U` はスケーリングを含まないため、`opt.slu.equilibration = true`
は `unsupported_options` になります。

取り出した因子を使って解く関数もあります(A⁻¹ = Q U⁻¹ L⁻¹ P)。

```cpp
std::vector<double> x   = A.lu_solve(L, U, p, q, b);        // A x = b
std::vector<double> row = A.lu_inverse_row(L, U, p, q, i);  // A^{-1} の第 i 行

std::vector<double> x2;
vcp::lu_apply_result r1 = A.lu_solve_with_info(L, U, p, q, b, x2);
```

`lu_apply_result` の status は `success`、`dimension_mismatch`、
`singular_factor`(U の対角が零)、`inconclusive_division`(区間型などで U の
対角が零でないことを確定できない)、`invalid_input`、`internal_error` です。

自由関数 `lu(A, L, U, p, q)` は `A.lu(L, U, p, q)` と同じです。

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

`ldl` は status が `success` 以外のとき例外を投げます。`ldl_with_info` は
status(`vcp::sparse_ldl_status`)を返します。L / D / p が有効なのは
`success` または `zero_pivot` のときだけです。

| status | 意味 |
| --- | --- |
| `success` | 分解完了 |
| `structural_singularity` | 構造的に空の行・列がある(`structural_empty_at`) |
| `zero_pivot` | 零ピボットを飛ばして分解を最後まで続けた(`first_zero_pivot`) |
| `inconclusive_pivot_test` | ピボットの判定を確定できなかった(区間型で 0 を含む要素など。`inconclusive_at`) |
| `not_symmetric` | 対称であることを確認できない |
| `invalid_options` / `invalid_input` | 不正なオプション・入力 |
| `internal_error` | 内部エラー |

オプション `ldl_options<T>`:

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 分解法 |
| `ordering` | `auto_select` | 順序付け。選ばれた順序付けは `ordering_used` に記録されます。将来の版で選択が変わり得るため、再現性が必要な場合は明示してください |
| `check_symmetry` | `true` | 対称性の検査。下三角だけを与える場合は `false` |
| `symmetry_tol` | 1e-12 | 対称性の許容誤差 |
| `zero_pivot_tol` | 0 | 零ピボットとみなす閾値 |

慣性(正・負・零の固有値の個数 n₊, n₋, n₀)は LDLᵀ 分解から求めます。

```cpp
auto ir  = A.inertia();             // 失敗時に例外
auto ir2 = A.inertia_with_info();   // status を返す
```

分解が `zero_pivot` の場合は失敗ではなく、零固有値として n₀ に数えます。
符号を確定できない場合は `inconclusive_sign` を返します。

A − σB(B を省略すると A − σI)の慣性を、σ を変えて繰り返し求める場合は、
σ に依存しない処理を 1 回で済ませるハンドルを使います。

```cpp
auto h = A.ldl_shift_setup_with_info(B);
if (h.valid()) {
    auto i1 = h.inertia_at(sigma1);
    auto i2 = h.inertia_at(sigma2);
}
```

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

| status | 意味 |
| --- | --- |
| `success` | 分解成功 |
| `not_positive_definite` | 正定値でない(`failure_at`。MATLAB の flag は `failure_at + 1`) |
| `inconclusive_pivot_test` | ピボットの符号を確定できなかった(`inconclusive_at`) |
| `not_symmetric` | 対称であることを確認できない |
| `invalid_input` / `invalid_options` | 不正な入力・オプション |
| `internal_error` | 内部エラー |

オプション `chol_options<T>` の `ordering` の既定は `auto_select`(AMD)
です。選ばれた順序付けは `ordering_used` に記録されます。`natural` を
指定すると置換なしで分解します。

Cholesky 分解の内部のアルゴリズムは [kernels/tsparse_ldl.md](kernels/tsparse_ldl.md) を参照してください。

## 近似逆行列(AINV・FSAI)

反復法の前処理などに使う近似逆行列を構成します。近似逆行列 R そのものは
構成せず、ベクトルへの作用 `apply` と、‖I − R A‖∞ の推定値(保証なし)を
返す関数を提供します。

| 手法 | 構成 | 因子 | 作用 |
| --- | --- | --- | --- |
| AINV | `A.ainv_with_info(Z, W, D, opt)` | Z, W は単位上三角、D は対角 | `A.ainv_apply(Z, W, D, r, z)`: z = Z D⁻¹ Wᵀ r |
| FSAI(対称正定値) | `A.fsai_with_info(U, D, perm, opt)` | U は単位上三角、D は対角 | `A.fsai_apply(U, D, perm, r, z)`: z = P U D⁻¹ Uᵀ Pᵀ r |
| 適応 FSAI | `A.fsai_adaptive_with_info(U, D, perm, opt)` | FSAI と同じ | `fsai_apply` を使う |

推定値は `ainv_residual_norm_estimate` / `fsai_residual_norm_estimate` で
求めます。いずれも status を返し、例外は投げません。

近似逆行列の内部のアルゴリズムは [kernels/sparse_approx_inverse.md](kernels/sparse_approx_inverse.md) を参照してください。

## 固有値問題

疎行列の固有値問題では、求める個数 k を指定します。

```cpp
std::vector<double> lam = A.eigs(k);
auto r = A.eigs_with_info(k, opt);   // r.eigenvalues, r.eigenvectors, r.converged

auto rg = A.eigs_with_info(B, k, opt);   // 一般化固有値問題 A x = λ B x

vcp::matrix<double, vcp::mats2<double> > E, V;
eigs(A, E, V, k, opt);                   // E: k x 1, V: n x k
```

オプション `eig_options<T>` の主なフィールド:

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 解法 |
| `target` | `smallest_algebraic` | どの固有値を求めるか(`largest_magnitude`、`smallest_magnitude`、`largest_algebraic`、`smallest_algebraic`、`target_magnitude`、`target_real`) |
| `structure` | `auto_detect` | 構造の指定(`general`、`symmetric`) |
| `use_shift` / `shift` | `false` / 0 | shift-invert を使う場合のシフト |
| `tol` | 1e-12 | 収束判定 |
| `max_iter` | 1000 | 反復回数の上限 |
| `allow_complex_pairs` | `false` | 実行列の複素共役固有値の対を返すか |

`auto_select` は次の順に解法を選びます。

| 条件 | 解法 |
| --- | --- |
| `use_shift == true` | shift-invert |
| n ≤ 300 で密行列への変換が許可されている | 密行列の解法(対称: Jacobi 法、非対称: 実 Schur 分解) |
| 確実に対称と判定できる | thick-restart Lanczos 法 |
| それ以外 | Krylov–Schur 法 |

`method` に `lanczos`、`arnoldi`、`thick_restart_lanczos`、`krylov_schur`、
`shift_invert_lanczos`、`shift_invert_arnoldi`、`dense_fallback_explicit` を
明示することもできます。

`allow_complex_pairs = true` にすると、複素共役の対を
`eigenvalues`(実部)と `eigenvalues_imag`(虚部)の隣り合う 2 つとして
返します。固有ベクトル x = u + i v の実部 u と虚部 v は
`eigenvectors` の隣り合う 2 本に入ります。

一般化固有値問題では、`structure = symmetric` を指定するか、shift-invert を
使ってください。既定のオプションでは、B が対角でない場合に空の結果を返します。

`converged == true` は、返した各固有対が残差の判定を満たすことを意味します。
Krylov 部分空間法は開始ベクトルに直交する固有空間を見られないため、
「求める側の k 個をすべて見つけた」ことは保証しません。

固有値問題の内部のアルゴリズムは [kernels/tsparse_eigs.md](kernels/tsparse_eigs.md) を参照してください。

## `spimats` の詳細

`vcp::spimats<T>` は区間疎行列の精度保証付き計算の policy です。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/spimats.hpp>
#include <vcp/spmatrix.hpp>

typedef kv::interval<double> I;

vcp::spmatrix<I, vcp::spimats<double> > A;
A.resize(n, n);
A.add(i, j, I(lo, hi));
A.finalize();

std::vector<I> b(n, I(1.0));
std::vector<I> x = A.solve(b);    // 解の包含
```

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
| `solve` / `solve_with_info` / `lss` | 解の包含 |
| `eigs` / `eigs_with_info`(一般化を含む) | 対称行列の固有値の包含(固有ベクトルは返しません) |
| `is_spd()` | 区間内のすべての対称行列が正定値であることを証明できたら `true` |

`is_spd()` の `false` は「正定値でない」という主張ではなく、十分条件が
成り立たなかったことを意味します。

`AP` と `DP` の選択は、包含の成否・包含の幅・計算速度にだけ影響し、
精度保証の正しさには影響しません。

例外の扱いは他の policy と異なります。保証に失敗したときは、
`solve_with_info` などの `_with_info` 系の関数でも例外を投げます。

| 状況 | 例外 |
| --- | --- |
| 次元の不一致などの誤用 | `vcp::dimension_error` など |
| 近似計算の失敗 | `vcp::state_error` |
| 精度保証の十分条件が成り立たない | `vcp::verification_error` |

成功したときは `converged` は常に `true` です。結果のうち保証が付いているのは、
`solve` 系の `x` と固有値問題の `eigenvalues` だけです。

`is_spd()`、保証層のオプション `spi_opt`、計算の記録 `last_verify_info` は
`spimats` 固有のメンバです。`spmatrix` は policy を protected 継承している
ため、これらは policy のオブジェクトを直接使って呼び出します。

```cpp
vcp::spimats<double> A;
A.resize(n, n);
A.add(i, j, I(lo, hi));
A.finalize();
bool ok = A.is_spd();
```

`last_verify_info` には近似計算の結果などが点型のまま記録されます。これは
計算の記録であり、対象の行列についての主張ではありません。

ノルムの上界には `norminf` / `normone` を使います。区間型では
`B.get(0, 0).upper()` が保証付きの上界です。`normtwo` は近似計算であり、
保証には使えません。

精度保証の手順と混合演算は [kernels/spimats_kernel.md](kernels/spimats_kernel.md) を参照してください。

## `spumar` の詳細

`vcp::spumar` は `spmats<double>` を継承し、疎 LU を UMFPACK に、疎固有値を
ARPACK に委譲する policy です。

```cpp
#include <vcp/spumar.hpp>
#include <vcp/spmatrix.hpp>

vcp::spmatrix<double, vcp::spumar> A;
A.resize(n, n);
/* add ... */
A.finalize();
auto r  = A.solve_with_info(b);         // UMFPACK
auto er = A.eigs_with_info(4, eopt);    // ARPACK
```

インストールとリンクの例(Ubuntu):

```bash
sudo apt install libsuitesparse-dev libarpack2-dev
g++ -std=c++11 -I/path/to/libs example.cpp -lumfpack -larpack
```

| 操作 | 委譲先 |
| --- | --- |
| `solve` / `solve_with_info` | UMFPACK(`method` の指定によらない) |
| `lu` / `lu_with_info` | UMFPACK |
| `eigs_with_info`(標準) | ARPACK。条件を満たさない指定(シフトなしの `smallest_magnitude` など)は `spmats` の実装で計算 |
| `eigs_with_info(B, ...)`(一般化) | ARPACK(shift-invert)。`use_shift` または `smallest_magnitude` のとき |
| LDLᵀ 分解・慣性・LU 因子の消費(`lu_solve` など)・行列演算 | 委譲しない(`spmats` の実装) |

UMFPACK には対応する設定がないため、`lu_extract_options::slu` の `method` /
`ordering` は適用されません。

オプションは `spumar_options` で、policy の公開メンバ `muar_opt` に
あります。

| フィールド | 既定 | 意味 |
| --- | --- | --- |
| `ncv` | 0(自動) | ARPACK の部分空間の次元 |
| `umf_irstep` | −1(UMFPACK の既定) | UMFPACK の反復改良回数 |
| `umf_pivot_tol` | −1.0(UMFPACK の既定) | UMFPACK のピボット閾値 |

`spmatrix` 経由で使う場合は、`vcp::spumar::set_default_options(o)` で
プロセス全体の既定値を設定します。この関数はスレッド安全ではありません。

UMFPACK・ARPACK への委譲の仕組みは [kernels/sparse_external.md](kernels/sparse_external.md) を参照してください。

## `spcmodumar` の詳細

`vcp::spcmodumar` は `vcp::spumar` を継承し、さらに疎 Cholesky 分解
(`chol` / `chol_with_info`)を CHOLMOD に委譲する policy です。`spumar` の
委譲(UMFPACK・ARPACK)はそのまま引き継ぎます。

```cpp
#include <vcp/spcmodumar.hpp>
#include <vcp/spmatrix.hpp>

vcp::spmatrix<double, vcp::spcmodumar> A, L;
std::vector<int> p;
A.chol(L, p);                       // CHOLMOD
auto s = A.solve_with_info(b);      // UMFPACK
```

リンクには `-lumfpack -larpack -lcholmod` が必要です(`spumar` を継承する
ため、`-lcholmod` だけでは足りません)。

分解の規約(`Pᵀ A P = L Lᵀ`、p は new→old)は `spmats` と同じです。
ただし、CHOLMOD は順序付けに後処理の並べ替えを合成するため、返る p は
`spmats` と一致するとは限りません。返った p で恒等式は常に成り立ちます。

| `ordering` | CHOLMOD の設定 |
| --- | --- |
| `natural` | 自然順序 |
| `auto_select` / `amd` | AMD |
| `rcm` | VCP の RCM を与える |
| `nested_dissection` | NESDIS |

- `pd_tol` に 0 以外を指定すると `invalid_options` になります。
- `method_used` に記録される `supernodal` は `spcmodumar` だけの値です。
- supernodal 経路の `nnz_L` は明示的に格納した零を含むため、simplicial
  経路の値と直接比べられません。

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
