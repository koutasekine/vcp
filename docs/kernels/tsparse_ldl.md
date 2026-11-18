# tsparse — 疎 LDLᵀ 分解・Cholesky 分解・慣性

この文書では、対称な疎行列の分解を説明します。

| 内容 | 場所 |
| --- | --- |
| LDLᵀ 分解 | `vcp/tsparse/tsparse_sparse_ldl.hpp`、`detail/tsparse_sparse_ldl_*_impl.hpp` |
| Cholesky 分解 | `vcp/tsparse/tsparse_sparse_chol.hpp` |
| 慣性の計算 | `vcp/spmats_base/spmats_ldl_impl.hpp` |
| shift 付き LDLᵀ の handle | `vcp/spmats_base/spmats_ldl_shift.hpp`、`detail/tsparse_sparse_ldl_shift_plan_impl.hpp` |

`spmats` の `ldl`・`inertia`・`chol`・`is_spd`・`ldl_shift_setup` から呼ばれ
ます。利用者が直接呼ぶことは想定していません。使い方は
[spmatrix.md](../spmatrix.md) の「LDLᵀ 分解と慣性」「Cholesky 分解」を参照して
ください。

## LDLᵀ 分解

### 呼び出し方

```cpp
#include <vcp/tsparse/tsparse_sparse_ldl.hpp>

vcp::sparse_ldl_options<double> opt;
vcp::sparse_ldl_result<double, int> r =
    vcp::sparse_ldl_factorize_with_info(n, col_ptr, row_ind, val, opt);
// r.L_col_ptr, r.L_row_ind, r.L_val(単位下三角、対角の 1 も格納、CSC)
// r.D_diag, r.D_sub, r.D_block2(1×1・2×2 のブロック対角)
// r.perm(new→old)
```

入力は n × n の CSC 配列です。対称と仮定し、数値分解では下三角(対角を
含む)だけを読みます。上三角は対称性の検査だけに使います。`spmats` からは
`sparse_lu_make_csc_storage` で CSC に変換して渡します。

結果は Pᵀ A P = L D Lᵀ(ベクトル形で A(p, p) = L D Lᵀ)の形です。L は単位
下三角、D は 1×1 と 2×2 のブロックからなるブロック対角です。D の表現は
次の通りです。

| 配列 | 内容 |
| --- | --- |
| `D_diag[k]` | d_kk |
| `D_sub[k]` | 2×2 ブロックの (k+1, k) 要素 |
| `D_block2[k]` | k が 2×2 ブロックの先頭なら 1 |

`sparse_ldl_factorize_with_info` は例外を投げず、結果を `status` で返します。
L・D・perm が有効なのは、`status` が `success` か `zero_pivot` のときだけです。

### 計算の流れ

1. **入力の検査**:CSC の形が正しくなければ `invalid_input`。
2. **オプションの検査**:範囲外の値は `invalid_options`。
3. **対称性の検査**(`check_symmetry = true` のとき):すべての (i, j) で
   |a_ij − a_ji| ≤ `symmetry_tol` と言えなければ `not_symmetric`。転置の
   パターンを作って列ごとに突き合わせるので、計算量は O(n + nnz)。
4. **構造的な特異性の検査**:下三角に要素が 1 つもない行・列があれば
   `structural_singularity`(分解はしない)。
5. **記号分解**:並べ替え P0 を決める。`supernodal` では、さらに消去木などを
   作る(下記)。
6. **数値分解**:方法(`method`)に応じた経路で分解する。

### 並べ替え(`ordering`)

`natural`、`rcm`、`amd`、`nested_dissection`、`nested_dissection_ml` が
あり、`auto_select` は `nested_dissection_ml` になります。中身は疎 LU と同じ
関数です([tsparse_lu.md](tsparse_lu.md) の「並べ替え」)。`colamd` は
非対称な LU のための方法なので、ありません。実際に使った並べ替えは
`ordering_used` に記録します。

### 記号分解

`supernodal` の方法では、下三角のパターンから次を作ります。要素の値は
使いません。

- 対称な下三角のパターン(対角は必ず含める)
- 消去木(Liu の方法、経路圧縮付き)と、その後順
- ピボットの交換をしない場合の L の各列の非零の個数
- fundamental supernode の分割。supernode s = [c0, c1) の各列 j の非零の
  位置は {j, j+1, …, c1−1} ∪ R_s(R_s は supernode の下の行の一覧)なので、
  R_s だけを持てば L のパターンが分かる
- 数値分解に必要な作業領域の大きさ

`ldl_amalgamation = relaxed` のときは、隣り合う親子の supernode を、次の
いずれかが成り立てばまとめます(既定は `off`)。w はまとめた後の幅、z は
まとめたことで増える零の割合(%)です。

- w ≤ n0、または w ≤ n1 かつ z < z0、または w ≤ n2 かつ z < z1、または z < z2

既定値は (n0, n1, n2) = (4, 16, 48)、(z0, z1, z2) = (80, 10, 5) です。
まとめた部分の零は数値分解の間も厳密に 0 のままで、出力からは除きます。
L・D のパターンは変わらず、和の順だけが変わります。

`baseline_dynamic` では、ピボットの交換で構造が変わるので、並べ替えだけを
決めます。

### 数値分解の方法(`method`)

| 値 | 方法 |
| --- | --- |
| `baseline_dynamic` | 列ごとの left-looking 法 |
| `supernodal` | supernode のパネルによる left-looking 法 |
| `auto_select` | `baseline_dynamic` |

#### `baseline_dynamic`

第 k 列の作業ベクトルを、A の第 k 列から、確定したピボットの更新をピボットの
順に引いて作り(left-looking)、ピボットを選んで L と D を書き出します。
ピボットの交換は、行と列の番号の対応表を書き換えて行います(確定した L の
行のデータは動かしません)。

n ≤ 64 で `pivoting = bk` のときは、密行列の分解に任せます(`dense_delegated`
が `true` になります)。疎な経路は、密行列の分解と同じピボットの列を、
要素ごとに同じ浮動小数点演算を同じ順で行って作ります。

#### `supernodal`

supernode s の列 [pa, pb)(幅 w)と下の行 R_s(h 行)を、(w + h) × w の密な
パネルとして持ちます。s を消去する前に、パネルの行に関係する完了済みの
supernode d ごとに

- W = L_d(I, :) D_d(I は R_d のうち [pa, pb) にある行)
- C(I, I) −= W L_d(I, :)ᵀ(対称な部分。`tgemm` または `tgemmtr`)
- C(J, I) −= W L_d(J, :)ᵀ(J は R_d のうち pb 以上の行。`tgemm`)

で更新し、パネルの中は 1 列ずつ right-looking で消去します。

密なブロックの演算(`vcp::tgemm` など)は、要素型が浮動小数点型
(`std::is_floating_point`)で、更新する supernode の幅が `ldl_min_block_size`
(既定 32)以上のときだけ使います。それ以外は、`baseline_dynamic` ×
`none` と同じ演算を同じ順で行うスカラーの経路を使います。対称な部分の更新に
使う演算は `diag_kernel` で選びます(`auto_select` は `gemm`)。

Bunch–Kaufman のピボット選択は、パネルの中でだけ交換します。パネルの外の
列との交換が必要になると、この経路は `pivot_out_of_panel` で止まります。
`sparse_ldl_factorize_with_info` は、そのとき同じ記号分解で
`baseline_dynamic` を使って分解し直し、`supernodal_fell_back = true` と
`fallback_from` に記録します。

### ピボット選択(`pivoting`)

| 値 | 動作 |
| --- | --- |
| `bk`(既定) | Bunch–Kaufman のピボット選択(1×1・2×2、対称な交換) |
| `none` | 交換しない。対角の 1×1 ピボットを順に使う |

`bk` では、α = (1 + √17) / 8 として、第 k 列について次の規則で選びます。

1. λ = max_{i>k} |a_ik|(最大の行を r とする)。λ = 0 なら a_kk を 1×1 の
   ピボットにする。
2. |a_kk| ≥ αλ なら、a_kk を 1×1 のピボットにする。
3. σ = max_{i≥k, i≠r} |a_ir| として、|a_kk| σ ≥ αλ² なら、a_kk を 1×1 の
   ピボットにする。
4. |a_rr| ≥ ασ なら、k と r を交換して a_rr を 1×1 のピボットにする。
5. それ以外は、(k, r) を 2×2 のピボットにする。

比較はすべて「成り立つと言えるときだけ通す」形で書かれています。区間型
などで、どの規則が成り立つかを決められないときは、`inconclusive_pivot_test`
で止まり、止まった列を `inconclusive_at` に記録します。

1×1 のピボットの絶対値、または 2×2 のピボットの行列式の絶対値が
`zero_pivot_tol`(既定 0)以下と言えるときは、その列を零のピボットとして
飛ばして数え、分解を最後まで続けます。このとき `status` は `zero_pivot`、
最初の列は `first_zero_pivot` です。`none` で、零とも零でないとも言えない
ピボットがあると `inconclusive_pivot_test` で止まります。

`bk` で `zero_pivot_tol = 0` のとき、零のピボットになるのは、残りの第 k 列が
対角も含めてすべて零のときだけです(規則 2〜5 で選ぶピボットは零になり
ません)。その列は消去しなくても分解は正しく、零の固有値 1 つに対応します。
`none` では、対角が零ならその列に非零があっても飛ばすので、そのときの
L D Lᵀ は Pᵀ A P と一致しません。

### オプション

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 数値分解の方法 |
| `pivoting` | `bk` | ピボット選択 |
| `ordering` | `auto_select` | 並べ替え |
| `check_symmetry` | `true` | 対称性の検査をする |
| `symmetry_tol` | 10⁻¹² | 対称とみなす \|a_ij − a_ji\| の上限 |
| `zero_pivot_tol` | 0 | 零のピボットとみなす上限(0 は厳密に零と言えるものだけ) |
| `pivot_threshold` | 0 | 予約。0 と言えない値は `invalid_options` |
| `ldl_min_block_size` | 32 | `supernodal` で密なブロックの演算を使う幅の下限 |
| `diag_kernel` | `auto_select` | `supernodal` の対称な部分の更新(`gemm`・`gemmtr`) |
| `ldl_amalgamation` | `off` | supernode をまとめる |
| `ldl_amalg_n0` 〜 `ldl_amalg_z2` | 上記 | まとめる規則の定数 |
| `ndml_*` | −1 | `nested_dissection_ml` のパラメータ |

### 結果

主なフィールドは次の通りです。これ以外は `supernodal` の経路の診断用です。

| フィールド | 内容 |
| --- | --- |
| `status` | 状態 |
| `n_pivots_1x1`、`n_pivots_2x2` | 使ったピボットの数(飛ばした零のピボットを除く) |
| `first_zero_pivot` | 最初の零のピボットの列(なければ −1) |
| `inconclusive_at` | 判定できずに止まった列(なければ −1) |
| `structural_empty_at` | 最初の空の行・列(なければ −1) |
| `nnz_L` | L の格納数 |
| `ordering_used`、`method_used` | 実際に使った並べ替えと方法 |
| `dense_delegated` | 密行列の分解に任せたか |
| `supernodal_fell_back`、`fallback_from` | `baseline_dynamic` で分解し直したか、その理由 |
| `growth_log2`、`growth_valid` | ⌊log₂(max\|d\| / max\|a_ii\|)⌋。`growth_valid` が `true` のときだけ有効 |

| status | 意味 |
| --- | --- |
| `success` | 成功 |
| `zero_pivot` | 零のピボットを飛ばして、分解を最後まで行った |
| `structural_singularity` | 空の行・列がある |
| `inconclusive_pivot_test` | ピボットの判定ができずに止まった |
| `pivot_out_of_panel` | `supernodal` でパネルの外との交換が必要になった |
| `not_symmetric` | 対称と言えない |
| `invalid_options`、`invalid_input` | オプション・入力の誤り |
| `internal_error` | 想定外の例外 |

`sparse_ldl_status_to_string` で文字列にできます。

## 慣性

対称行列 A の慣性(正・負・零の固有値の数)は、Sylvester の慣性法則に
より、Pᵀ A P = L D Lᵀ の D の慣性と同じです。`spmats` の `inertia` は、
LDLᵀ 分解をしてから `policy_inertia_from_block_diagonal` で D を数えます。
数え方は次の通りです。

- 1×1 ブロック d_k:|d_k| ≤ tol と言えれば零、そうでなければ符号で正・負。
- 2×2 ブロック:det = d_k d_{k+1} − e_k²、tr = d_k + d_{k+1} として、
  - |det| ≤ tol と言えれば、零が 1 つと、tr の符号の固有値が 1 つ
    (|tr| ≤ tol なら零が 2 つ)。
  - det < 0 なら、正と負が 1 つずつ。
  - det > 0 なら、tr の符号の固有値が 2 つ。
- 要素のない列は零が 1 つ。
- |i − j| > 1 の要素や、2×2 ブロックの直後に続く副対角の要素があれば
  `not_block_diagonal`。

符号を決められないときは、数えずに `inconclusive_sign` とその位置を返します。
最後に、正・負・零の数の和が n であることを確かめます。

`inertia` の結果の `ldl_status` には、内部の LDLᵀ 分解の status が入ります。
`pivoting = none` で `ldl_status` が `zero_pivot` のときは、上で述べたとおり
分解が成り立たないので、`status` が `success` でも数えた値は A の慣性と
一致するとは限りません。

## shift 付き LDLᵀ の handle

A − σB(B を省けば A − σI)の慣性を、多くの σ について調べるための仕組み
です。`spimats` の固有値の包含の検証などで使います。

`policy_ldl_shift_setup_with_info` で handle を作るとき、σ によらない部分を
1 回だけ計算します。

- `spmats` から CSC への変換と、対称性の検査
- A と B の下三角のパターンを合わせたパターンと、A・B の各要素の行き先の表
  (`sparse_ldl_shift_plan`)。B を省いたときは対角の位置も用意する。
- 並べ替えと記号分解

その後、`inertia_at(σ)` は σ に依存する部分だけを計算します。

1. 合わせたパターンの上で、A − σB の値を 1 回の走査で作る(O(nnz))。
2. 記号分解を使い回して数値分解をする。
3. D を直接走査して慣性を数える。

`factor_at(σ, L, D, p)` は、L・D・p を作って返します。結果は、合わせた
行列を `ldl_with_info` で分解したものと同じです。

handle の既定のオプション(`vcp::ldl_shift_default_options<T>()`)は
`supernodal` × `none` です。`none` ではピボットの交換がないので、記号分解を
すべての σ で使い回せます。`inertia_at` の結果(`ldl_shift_inertia_result`)
の `inertia.ldl_status` が
`zero_pivot` のときは、慣性の節で述べたとおり、数えた値が A − σB の慣性と
一致するとは限りません。handle は作った後は読み取り専用で、作業領域を
呼び出しごとに渡せば、異なる σ を並列に計算できます。ある σ での失敗は
handle を無効にしません。

## Cholesky 分解

### 呼び出し方

```cpp
#include <vcp/tsparse/tsparse_sparse_chol.hpp>

vcp::sparse_chol_options<double> opt;
vcp::sparse_chol_result<double, int> r =
    vcp::sparse_chol_factorize_with_info(n, col_ptr, row_ind, val, opt);
// r.L_col_ptr, r.L_row_ind, r.L_val(下三角、正の対角を格納、CSC)、r.perm
```

結果は Pᵀ A P = L Lᵀ(ベクトル形で A(p, p) = L Lᵀ)の形です。L の対角は
正です。ピボットの交換はしないので、perm は並べ替えそのものです。入力の
扱い(下三角だけを読む、対称性の検査)は LDLᵀ 分解と同じです。L と perm が
有効なのは `status` が `success` のときだけです(途中までの因子は返しません)。

### 計算の方法

並べ替え P0 で並べ替えた下三角の CSC C = P0ᵀ A P0 を作り、2 段で分解します
(simplicial up-looking 法)。

1. **記号分解**:C のパターンから消去木を作り、各行 k について、消去木の上で
   到達できる列(ereach)を求めて、L の構造を決める。
2. **数値分解**:決めた構造に値を入れる。第 k 行は、それまでの L による疎な
   三角行列の解法で求め、対角 d について「d > `pd_tol` と言えれば √d を
   L の対角にする」。

数値分解では ereach を計算し直しません。密行列の分解への切り替えはありません。

ピボットの判定は 3 通りです。

- d > `pd_tol` と言える:続ける。
- d ≤ `pd_tol` と言える:`not_positive_definite` で止まり、その列を
  `failure_at` に記録する(MATLAB の `chol` の flag は `failure_at + 1`)。
- どちらとも言えない:`inconclusive_pivot_test` で止まる。

`pd_tol` の既定は 0 なので、厳密に零のピボットも正定値でないとして止まり
ます。空の行・列も `not_positive_definite` とし、`structural_empty_at` に
記録します。

### オプション・結果

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | `simplicial_uplooking` だけ |
| `ordering` | `auto_select` | LDLᵀ 分解と同じ。`auto_select` は `nested_dissection_ml` |
| `check_symmetry`、`symmetry_tol` | `true`、10⁻¹² | LDLᵀ 分解と同じ |
| `pd_tol` | 0 | ピボットを正とみなす下限 |
| `ndml_*` | −1 | `nested_dissection_ml` のパラメータ |

| status | 意味 |
| --- | --- |
| `success` | 成功 |
| `not_positive_definite` | 正定値でない(零のピボット、空の行・列を含む) |
| `inconclusive_pivot_test` | ピボットの判定ができずに止まった |
| `not_symmetric` | 対称と言えない |
| `invalid_options`、`invalid_input` | オプション・入力の誤り |
| `internal_error` | 想定外の例外 |

`method_used` が `supernodal` になるのは、`vcp::spcmodumar` で CHOLMOD に
任せて、CHOLMOD が supernodal 法を使ったときだけです。

`spmats` の `is_spd` は、Cholesky 分解が成功するかどうかで判定します。
