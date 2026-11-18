# 近似逆行列(AINV・FSAI)

この文書では、`spmats` の近似逆行列 R ≈ A⁻¹ の構成法を説明します。実装は
`tsparse` ではなく、policy の層にあります。

| 方法 | 場所 |
| --- | --- |
| AINV | `vcp/spmats_base/spmats_ainv.hpp`、`spmats_ainv_impl.hpp` |
| FSAI(静的) | `vcp/spmats_base/spmats_fsai.hpp`、`spmats_fsai_impl.hpp` |
| 適応 FSAI | `vcp/spmats_base/spmats_fsai_adaptive.hpp`、`spmats_fsai_adaptive_impl.hpp` |

利用者が直接呼ぶことは想定していません。使い方は [spmatrix.md](../spmatrix.md)
の「近似逆行列(AINV・FSAI)」を参照してください。

## 共通の性質

- R そのものは作らず、因子だけを返します。R をベクトルに掛ける `apply` と、
  ‖I − R A‖∞ を計算する `residual_norm_estimate` があります。
- 構成・`apply`・推定値の計算は例外を投げず、status(`success`・
  `invalid_input`・`internal_error`)で返します。正方でない入力や、因子の
  次元・構造の誤りは `invalid_input` です。
- 小さいピボットは、失敗にせずに置き換えます(下記の「ピボットの置き換え」)。
  置き換えた回数と最初の位置は、結果の診断に記録します。特異な行列でも
  構成は `success` で終わります。
- 厳密零の値は因子に格納しません。

### ピボットの置き換え

ピボット p の絶対値が `pivot_small_threshold`(既定 10⁻¹⁴)以下のとき、p を
符号を保って ±`pivot_lift_value`(既定 10⁻³)に置き換えます。p が厳密に零の
ときは +`pivot_lift_value` にします。

### ‖I − R A‖∞ の推定値

`residual_norm_estimate` は、R を作らずに 1 行ずつ計算します。第 i 行に
ついて

1. u = Zᵀ e_i(Z の第 i 行)、v = D⁻¹ u
2. r_i = W v
3. s_iᵀ = r_iᵀ A
4. est_i = Σ_j |s_ij − δ_ij|

として、est = max_i est_i を返します。すべての行を計算します。点の演算で
計算した値なので、真のノルムについての保証はありません。FSAI では
Z = W = U として同じ関数を使います。

## AINV

非対称な A に対して、R = Z D⁻¹ Wᵀ ≈ A⁻¹ を作ります(Benzi–Tůma, 1998 の
不完全双共役化)。Z と W は単位上三角、D は対角です。

### 構成

Z 側は A、W 側は Aᵀ について、同じ手順を 1 回ずつ行います(SDS 型の
left-looking 法)。M を A または Aᵀ として、第 i 列 z_i = e_i から始め、
確定した列 j < i について

1. p = m_jᵀ z_i(M の第 j 行と z_i の内積)
2. z_i ← z_i − (p / p_j) f_j

と更新します。更新で新しく生じた要素のうち、絶対値が `drop_tolerance`
(既定 10⁻²)未満のものは捨てます。すでにある要素(単位対角を含む)は
捨てません。列を確定するときは、

1. `max_nnz_per_column`(既定 0 = 制限なし)が正なら、対角以外で絶対値の
   大きい順にその個数だけ残す(同じ大きさなら先に入った方を残す)。
2. ピボット p_i = m_iᵀ z_i を計算し、必要なら置き換える。
3. 行の昇順に格納する。

D の対角には Z 側のピボット p_i を入れます(W 側のピボットは W の計算に
だけ使います)。

### apply

z = Z (D⁻¹ (Wᵀ r)) を計算します。Wᵀ は作らず、W の行を走査して Wᵀ r を
足し込みます。

### オプション・結果

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `drop_tolerance` | 10⁻² | 新しく生じた要素を捨てる閾値 |
| `max_nnz_per_column` | 0 | 1 列に残す要素の数の上限(単位対角を除く。0 は制限なし) |
| `pivot_small_threshold` | 10⁻¹⁴ | 置き換えるピボットの閾値 |
| `pivot_lift_value` | 10⁻³ | 置き換えたピボットの大きさ |

結果 `ainv_result` には、`status`、置き換えた回数 `n_pivot_modifications`
(Z 側と W 側の合計)、最初の位置 `first_modified_pivot`、`nnz_Z`、`nnz_W`
が入ります。

## FSAI(静的)

対称正定値の A を想定して、R = P U D̂⁻¹ Uᵀ Pᵀ ≈ A⁻¹ を作ります
(Janna ほか, 2015 の FSAIPACK の静的 FSAI)。U は単位上三角、D̂ は対角、
P は並べ替え(perm は new→old)です。平方根を使わない形にしています。

### 構成

1. **並べ替え**:A_p = Pᵀ A P を作る(`ordering`。既定は `natural`)。
2. **前処理の間引き**:|a_ij|² < τ_pre² |a_ii| |a_jj| の要素を除いた A′ の
   パターンを作る(`prefilter_tolerance` = τ_pre。0 なら除かない)。対角は
   除きません。
3. **パターン**:B₀ = I、B_i = Low(B_{i−1} A′)(Low は下三角部分)を
   `pattern_power` 回くり返し、S = B_k のパターン ∪ 対角とする。途中で
   nnz(B_i) ≥ `max_density` · nnz(A) になったら止める。
4. **各行の計算**:第 i 行のパターン P_i について、密な小行列
   A_p[P_i, P_i] を集め、A_p[P_i, P_i] g = e_m を解く。解法は、ピボットの
   交換をしない消去法(対称なら平方根を使わない LDLᵀ)で、ピボットは必要なら
   置き換える。
5. **平方根を使わない形への変換**:d = g_m として、D̂_ii = 1/d、U の第 i 列を
   g/d とする。
6. **後処理の間引き**(`postfilter_tolerance` = τ_post、
   `postfilter_max_nnz_per_row`):|u_ij| < τ_post ‖u_i‖₂ の要素を除き、さらに
   上限の個数まで絶対値の大きい順に残す。除いた要素の集合 E について
   D̂_ii ← D̂_ii (1 + εᵀ A ε) と補正する。補正の係数が正と言えないときは
   補正せず、その行数を記録する。

各行の計算は独立なので、OpenMP で並列に計算します。各行は自分の領域だけに
書くので、結果はスレッド数によらず逐次の計算と一致します。この層の並列化を
止めるマクロは、policy の層に共通の `VCP_SPMATS_NOMP`(または `VCP_NOMP`)です
([tsparse.md](tsparse.md) の「OpenMP」)。

非対称な A に対しても、ブロックを上下両方の三角部分から集めてそのまま
分解するので、計算はできます。

### オプション・結果

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `pattern_power` | 1 | パターンを作るべき乗の回数 k |
| `prefilter_tolerance` | 0 | 前処理の間引きの閾値 |
| `max_density` | 5 | パターンの密度の上限(nnz(A) に対する比。0 は制限なし) |
| `ordering` | `natural` | 並べ替え(`sparse_chol_ordering`。`auto_select` も `natural`) |
| `postfilter_tolerance` | 0 | 後処理の間引きの閾値 |
| `postfilter_max_nnz_per_row` | 0 | 1 行に残す要素の数の上限(単位対角を除く。0 は制限なし) |
| `pivot_small_threshold`、`pivot_lift_value` | 10⁻¹⁴、10⁻³ | ピボットの置き換え |

結果 `fsai_result` には、`status`、`nnz_U`、パターンの大きさ `nnz_pattern`、
置き換えた回数と最初の行、後処理で除いた要素の数、補正しなかった行数が
入ります。

## 適応 FSAI

静的 FSAI と同じ形の因子 (U, D̂, perm) を、行ごとにパターンを広げながら
作ります(Janna ほか, 2015 の適応的なパターン生成)。

第 i 行について、パターン P̄(初めは空、または与えた初期因子 U0 の行)から
始め、次を最大 `max_iterations` 回くり返します。

1. A_p[P̄, P̄] g = −A_p[P̄, i] を解く(静的 FSAI と同じ消去法)。
2. 行の質 ψ_k = G_k[i, :] A G_k[i, :]ᵀ を計算する。
3. ψ_k ≤ `exit_tolerance` · ψ_0 なら終える。
4. ψ の勾配の成分の絶対値が大きい順に、パターンにない位置を
   `entries_per_step` 個加える。
5. `inner_drop_tolerance` = τ が正なら、|g_j|² ≤ τ² ‖g‖₂² の成分を除く。

最後に静的 FSAI と同じ後処理の間引きをし、D̂_ii = ψ_k とします。ψ が小さい
ときもピボットと同じく置き換えます。`apply` と推定値は、静的 FSAI の関数を
そのまま使います。

初期因子 U0 と perm0 を与える版では、`ordering` を使わず perm0 を使います。
U0 は各列の対角に 1 を明示的に持つ単位上三角でなければなりません。

各行は、前もって決めた大きさ(|P0_i| + `max_iterations` · `entries_per_step`
+ 1)の領域にだけ書くので、並列に計算しても結果は逐次の計算と一致します。

### オプション・結果

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `max_iterations` | 10 | 1 行あたりのパターンを広げる回数の上限 |
| `entries_per_step` | 5 | 1 回に加える位置の数 |
| `exit_tolerance` | 10⁻³ | 打ち切りの閾値 ε(ψ_k ≤ ε ψ_0) |
| `inner_drop_tolerance` | 0 | 反復中の間引きの閾値 |
| `ordering`、`postfilter_*`、`pivot_*` | 静的 FSAI と同じ | 静的 FSAI と同じ意味 |

結果 `fsai_adaptive_result` には、静的 FSAI と同じ項目に加えて、打ち切りの
判定で終えた行数 `n_rows_converged` と、上限まで回した行数 `n_rows_capped`
が入ります。

## 既定値について

オプションの既定値は、いずれも暫定の値です。
