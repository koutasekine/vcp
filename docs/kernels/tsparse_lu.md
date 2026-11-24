# tsparse — 疎 LU 分解

この文書では、`tsparse` の疎 LU 分解(`vcp/tsparse/tsparse_sparse_lu.hpp` と
`detail/tsparse_sparse_lu_*_impl.hpp`)を説明します。`spmats` の線形方程式
(直接法)、LU 分解、LU 分解の handle、shift-invert 型の固有値解法から
呼ばれます。利用者が直接呼ぶことは想定していません。使い方は
[spmatrix.md](../spmatrix.md) の「線形方程式」「LU 分解」を参照してください。

## 呼び出し方

```cpp
#include <vcp/tsparse/tsparse_sparse_lu.hpp>

vcp::sparse_lu_options<double> opt;
vcp::sparse_lu_factorization<double, int> fac =
    vcp::sparse_lu_factorize_with_info(A, opt);   // A は vcp::spmats<double>
if (fac.valid()) {
    std::vector<double> x = fac.solve(b);
    // 反復改良付きの解
    vcp::sparse_lu_refinement_info<double> ir;
    std::vector<double> y = vcp::sparse_lu_solve_refined(A, fac, b, opt, &ir);
}
```

`A` には、`rowsize()`・`columnsize()` などを持つ疎行列(`vcp::spmats<T, Index>`)
を渡します。`sparse_lu_factorize_with_info` は例外を投げず、結果を
`fac.info()` の `success` と `status` で返します(メモリ不足も
`memory_allocation_failed` として返します)。分解に失敗した `fac` で
`solve` を呼ぶと `vcp::state_error` を投げます。

## 計算の流れ

分解は次の順に進みます。

1. **入力の検査**:正方でない行列とオプションの範囲外の値は `invalid_input`。
2. **記号分解**(`sparse_lu_symbolic`):非零の位置だけを見て、列の並べ替え
   Q を決める。AᵀA の列消去木と、supernode の情報も作る。
3. **数値分解**:方法(`method`)に応じた経路で、Q で並べ替えた行列を分解する。
   行の並べ替え P は、数値分解の中のピボット選択で決まる。
4. **解法**(`solve`):前進・後退代入で解く。必要なら反復改良をする。

結果は P A Q = L U の形です。L は単位下三角、U は上三角です。

## 並べ替え(`ordering`)

列の並べ替え Q は、非零の位置だけから決めます(値は読みません)。
乱数は使わず、同順位は添字の小さい方を選ぶので、結果は決定的です。

| 値 | 方法 |
| --- | --- |
| `natural` | 並べ替えない |
| `rcm` | A + Aᵀ のグラフでの逆 Cuthill–McKee 法。開始点は George–Liu の擬似周辺点 |
| `amd` | A + Aᵀ のグラフでの近似最小次数法(Amestoy–Davis–Duff) |
| `colamd` | AᵀA の構造での列近似最小次数法(Davis–Gilbert–Larimore–Ng)。AᵀA は作らない |
| `nested_dissection` | A + Aᵀ のグラフの再帰的な 2 分割。幅優先探索の 1 つの段を分離集合にする |
| `nested_dissection_ml` | 多重レベルの nested dissection(下記) |
| `auto_select` | `nested_dissection_ml` |

`nested_dissection_ml` は、METIS などと同じ形の構成です。

1. 重い辺の matching でグラフを粗くしていく。
2. いちばん粗いグラフで、貪欲なグラフ成長法で 2 分割する。
3. 細かいグラフへ戻しながら、各段で Fiduccia–Mattheyses 法で境界を改良する。
4. 辺の切断から、貪欲な最小頂点被覆で頂点の分離集合を作る。

2 つの部分を先に(再帰的に)番号付けし、分離集合を最後にします。決めた大きさ
以下になった部分は、AMD で番号付けします。パラメータは `ndml_coarsen_stop`、
`ndml_fm_passes`、`ndml_balance_pct`、`ndml_leaf_size` で変えられます。
負の値はライブラリの既定値を表します。

## 数値分解の方法(`method`)

| 値 | 方法 |
| --- | --- |
| `baseline_gp` | 列ごとの left-looking 法(Gilbert–Peierls の方法) |
| `supernode_panel` | supernode とパネルによる left-looking 法(SuperLU の構成) |
| `supernodal` | multifrontal 法(オプトイン) |
| `auto_select` | `supernode_panel`。ただし `equilibration = true` のときは `baseline_gp` |

`baseline_gp` と `supernode_panel` の結果は、同じ形の CSC の L・U と置換で
持ち、解法・反復改良・因子の抽出は共通です。`supernodal` の結果は
supernode ごとの格納で持ち、解法と因子の抽出はその格納にも対応しています。

### `baseline_gp`

第 j 列について、次を行います。

1. A の第 j 列を、行の並べ替えに合わせて疎な作業ベクトルに散らす。
2. それまでに作った L の構造で深さ優先探索をし、第 j 列の計算に必要な
   列 k < j(到達集合)を求める。到達集合の列を位相順に使って、作業ベクトルを
   更新する(x ← x − L(:, k) x_k)。
3. ピボットを選ぶ。行を交換し、U の第 j 列と L の第 j 列を書き出す。

計算量は、実際に行う浮動小数点演算の数に比例します。

### `supernode_panel`

列を幅 `panel_size`(既定 8)のパネルに分けて処理します。L は supernode
(下三角部分の構造が同じ連続した列の集まり)ごとに密な台形で持ちます。
supernode は、各列のピボットを決めた後で動的に見つけます。各パネルについて、

1. A のパネルの列を、密な n × w の作業領域に散らす。
2. supernode のグラフで 1 回の深さ優先探索をし、パネルを更新する supernode を
   求める。グラフは Eisenstat–Liu の対称な刈り込みで小さくする。
3. 更新する supernode ごとに、U の断片を連続した領域に集め(使わない列は
   零で埋める)、supernode の対角ブロックで三角行列の解法をしてから、
   行列積 1 回でパネルを更新する(Li 2005 の方法)。更新が細いときは
   スカラーの演算で更新する。
4. パネルの中の各列を、`baseline_gp` と同じ順で 1 列ずつ分解する。

行列積には `vcp::tblas_blocked::gemm`(`panel_gemm_kernel = blocked`、既定)
または `vcp::tgemm`(`reference`)を使います。`kv::dd` で尾崎スキームの
行列積を使うには、分解を実体化する前に `<vcp/tblas/tblas_dd.hpp>` を include
します。include しない場合は汎用の行列積を使います(結果は正しく、遅く
なります)。supernode の幅の上限は `supernode_panel_maxsup`(既定 64)です。

演算の順が `baseline_gp` と違うので、結果はビット単位では一致しません。

### `supernodal`

multifrontal 法です。supernode ごとに密な frontal 行列を組み立て(A の要素と、
子の frontal 行列から送られた寄与の和)、ピボットのブロックの LU 分解、
2 回の三角行列の解法、1 回の行列積で Schur 補行列を計算し、それを親の
frontal 行列へ送ります。行の並べ替えは `baseline_gp` の分解で決めたものを
使い、frontal 行列の中では対角をピボットにします。零に近いピボットがあると、この経路の結果は使わず、
`baseline_gp` の結果を使います。

関係するオプション(`supernodal_self_symbolic`、`supernodal_inplace_frontal`、
`supernodal_native_check_residual`、`supernodal_numeric_diagnostic_leftlooking`)
は、この経路だけが読みます。

## ピボット選択(`pivoting`)

| 値 | 動作 |
| --- | --- |
| `threshold_partial`(既定) | 各列で、作業ベクトルの未確定の行のうち絶対値が最大の要素をピボットにする(部分ピボット選択) |
| `none`、`diagonal` | 対角の位置をそのままピボットにする |

ピボットの絶対値が `absolute_pivot_tolerance`(既定 0)より大きいと言えない
とき、分解は `numerical_singularity` で止まります。この判定は「大きいと
言えるときだけ通す」形なので、区間型では 0 を含むピボットは割り算の前に
はじかれます。

`threshold_partial` は、`baseline_gp` と `supernode_panel` では絶対値最大の
要素を選びます。`supernodal` の一部の経路は、対角の要素の絶対値が列の最大の
0.1 倍以上のときに対角をそのまま使います(0.1 は内部の定数で、オプションでは
変えられません)。

## 平衡化(`equilibration`)

`equilibration = true` のとき、行と列の最大値によるスケールをかけます。

- Dr[i] = 1 / max_j |a_ij|(零の行では 1)
- Dc[j] = 1 / max_i |Dr[i] a_ij|(零の列では 1)

`supernode_panel` を明示した場合は `not_implemented` です。

## 解法と反復改良

`fac.solve(b)` は、次の順に計算します。

1. Dr を掛ける。
2. 行の並べ替え P を適用する。
3. L で前進代入する(単位対角)。
4. U で後退代入する。
5. 列の並べ替え Q を戻す。
6. Dc を掛ける。

`sparse_lu_solve_refined(A, fac, b, opt, &info)` は、
`opt.iterative_refinement = true`(既定)のとき、同じ型 `T` で反復改良を
します。

```
x = solve(b)
k = 1, 2, ..., iterative_refinement_max_iterations(既定 20)について:
    r = b − A x(元の A で計算)
    x = x + solve(r)
    ‖r‖₂ ≤ tol · max(‖b‖₂, 1) なら終了(tol の既定は epsilon<T>())
    ‖r‖₂ が前回から減らなければ終了
```

残差が最小だった x を返します。結果の `info` には、実行したか(`performed`)、
収束したか(`converged`)、反復回数、最初と最後の残差が入ります。収束
しなかったときの `status` は `numerical_instability_suspected` です。
高い精度の型で残差を計算することはしません。
`iterative_refinement = false` のときは `fac.solve(b)` と同じです。

`spmats` の線形方程式(直接法)は、`iterative_refinement` が `true` のとき
この関数を使います。

## 因子の抽出

```cpp
#include <vcp/tsparse/tsparse_sparse_lu_extract.hpp>

vcp::sparse_lu_extracted<double, int> ex = vcp::sparse_lu_extract_factors(fac);
if (ex.status == vcp::sparse_lu_extract_status::success) {
    // ex.L_col_ptr, ex.L_row_ind, ex.L_val(単位下三角、対角の 1 も格納)
    // ex.U_col_ptr, ex.U_row_ind, ex.U_val
    // ex.p, ex.q(new→old、A(p, q) = L U)
}
```

`spmats` の `lu` / `lu_with_info` はこの関数を使います。規約は
P A Q = L U、置換ベクトルは new→old で、MATLAB の `[L,U,P,Q] = lu(A)` と
同じ向きです。各列の行添字は昇順で、厳密零は格納しません。
平衡化で Dr・Dc を使った分解は、`unsupported_options` を返します。

## 逆行列の行・列

```cpp
#include <vcp/tsparse/tsparse_sparse_lu_inverse_row.hpp>

vcp::sparse_lu_inverse_row_extractor<double, int> R(fac);
vcp::sparse_lu_inverse_row_result<double, int> r = R.extract(i);   // A⁻¹ の第 i 行
vcp::sparse_lu_inverse_col_extractor<double, int> C(fac);
vcp::sparse_lu_inverse_col_result<double, int> c = C.extract(j);   // A⁻¹ の第 j 列
```

A⁻¹ ≈ Dc Q U⁻¹ L⁻¹ P Dr から、1 行(1 列)だけを求めます。

- 第 i 行:Uᵀ y = e_jc(jc は第 i 列の並べ替え後の位置)、Lᵀ z = y の 2 回の
  三角行列の解法。Uᵀ と Lᵀ は extractor を作るときに 1 回だけ作る。
- 第 j 列:L y = Dr[j] e_ir、U z = y。

どちらも、右辺の非零から三角行列のグラフで到達できる成分だけを計算します
(Gilbert–Peierls の疎な三角行列の解法)。計算量は実際の演算数に比例し、
n や nnz(L) + nnz(U) にはよりません。

入出力の添字は元の(並べ替える前の)座標です。点の型(`float`、`double`、
`long double`、`kv::dd`、`kv::mpfr<N>`)だけが対象です。`supernodal` の分解は、
CSC の L・U も残っている場合だけ使えます。

## オプション

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 数値分解の方法 |
| `ordering` | `auto_select` | 列の並べ替え |
| `pivoting` | `threshold_partial` | ピボット選択 |
| `absolute_pivot_tolerance` | 0 | ピボットの絶対値の下限 |
| `zero_tolerance`、`near_zero_tolerance` | ε、√ε | `supernodal` の経路で、零・零に近いピボットを判定する閾値 |
| `equilibration` | `false` | 平衡化 |
| `iterative_refinement` | `true` | 反復改良 |
| `iterative_refinement_max_iterations` | 20 | 反復改良の最大回数 |
| `iterative_refinement_tolerance` | ε | 反復改良の目標の相対残差 |
| `panel_size` | 8 | `supernode_panel` のパネルの幅 |
| `supernode_panel_maxsup` | 64 | `supernode_panel` の supernode の幅の上限 |
| `panel_gemm_kernel` | `blocked` | `supernode_panel` の行列積 |
| `blas_min_block_size` | 16 | `supernode_panel` で行列積を使う大きさの目安 |
| `supernode_relaxation` | 0 | supernode をまとめる緩和の量(0 はまとめない) |
| `ndml_*` | −1 | `nested_dissection_ml` のパラメータ |

ε は `vcp::tsparse_scalar::epsilon<real_type>()` です。

## 結果の情報

`fac.info()`(`sparse_lu_info<T, Index>`)の主なフィールドは次の通りです。
これ以外のフィールドは、`supernodal` の経路の診断用です。

| フィールド | 内容 |
| --- | --- |
| `success`、`status` | 成否と状態 |
| `method_used` | 実際に使った方法 |
| `n`、`nnz_L`、`nnz_U` | 次数と、L・U の格納数 |
| `number_of_supernodes` | supernode の数 |
| `growth_factor` | max\|U\| / max\|A_eff\|(A_eff は分解した行列)。`supernodal` の一部の経路では計算しない |

主な `status` は次の通りです。`sparse_lu_status_to_string` で文字列に
できます。

| status | 意味 |
| --- | --- |
| `success` | 成功 |
| `invalid_input` | 正方でない、オプションの範囲外 |
| `not_implemented` | 実装していないオプションの組み合わせ |
| `numerical_singularity` | ピボットが `absolute_pivot_tolerance` より大きいと言えない |
| `structural_singularity` | 現在はどの経路も返さない |
| `memory_allocation_failed` | メモリ不足 |
| `internal_error` | 想定外の例外 |
