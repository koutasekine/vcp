# tsparse — 疎行列の固有値解法

この文書では、`spmats` の `eigs` / `eigs_with_info`(標準固有値問題
A x = λ x と一般化固有値問題 A x = λ B x)の内部の解法を説明します。
振り分けは `vcp/spmats_base/spmats_eigs.hpp`、解法の本体は `vcp/tsparse/` の
次の header にあります。利用者が直接呼ぶことは想定していません。使い方は
[spmatrix.md](../spmatrix.md) の「固有値問題」を参照してください。

| header | 内容 |
| --- | --- |
| `tsparse_eigs.hpp` | 列挙型(`eig_solver_method`、`eig_target`、`matrix_structure_hint` など)と定数 |
| `tsparse_eigs_types.hpp` | オプションと結果の型(`eig_options`、`eig_result`、`linear_solve_options`、`linear_solve_result` など) |
| `tsparse_thick_restart_lanczos.hpp` | thick-restart Lanczos 法(名前空間 `tsparse_experimental`) |
| `tsparse_krylov_schur.hpp` | Krylov–Schur 法(名前空間 `tsparse_experimental`) |
| `tsparse_real_schur.hpp`、`tsparse_dense_schur_driver.hpp` | 密行列の実 Schur 分解と、非対称な密行列の固有値計算 |
| `tsparse_dense_linalg.hpp`、`tsparse_dense_fallback.hpp` | 密行列の Jacobi 法など、密行列への切り替え |
| `tsparse_honest_termination.hpp`、`tsparse_restart.hpp`、`tsparse_eigen_selection.hpp` | 終了の判定、再出発の状態、固有値の選択 |
| `tsparse_generalized_shift_invert.hpp` | 一般化 shift-invert の作用素 |

`tsparse_lanczos.hpp`、`tsparse_arnoldi.hpp`、`tsparse_eigensolvers.hpp`、
`tsparse_b_inner_lanczos.hpp` は、古い解法の名残で、新しい経路からは
補助関数の一部だけが使われています。

`vcp::spumar` の固有値計算(ARPACK)は [sparse_external.md](sparse_external.md)
を参照してください。

## 振り分け(標準固有値問題)

`eig_options<T>::method` の既定は `auto_select` で、次の順に振り分けます。

1. `use_shift = true` なら shift-invert 型。対称と言えれば shift-invert Lanczos、
   そうでなければ shift-invert Arnoldi(Lanczos 型が使えない要素型でも Arnoldi)。
2. n ≤ 300(`eig_auto_dense_threshold`)で、`allow_dense_conversion = true`、
   n² ≤ `max_dense_size` なら密行列の解法。対称と言えれば Jacobi 法、
   そうでなければ実 Schur 分解。
3. 対称と言えれば thick-restart Lanczos 法、そうでなければ Krylov–Schur 法。

「対称と言える」は、`structure` が `symmetric` なら真、`general` なら偽、
`auto_detect`(既定)なら行列を見て確実に対称と言えるときに真です。
thick-restart Lanczos 法が使えない要素型(`kv::interval` など、
`numeric_limits` が特殊化されていない型)では、対称でも Krylov–Schur 法を
使います。`lanczos` を明示した場合、これらの型では解法を実行せず、対応して
いないことを結果で返します。

`method` を明示した場合は次の通りです。

| 値 | 解法 |
| --- | --- |
| `thick_restart_lanczos`、`lanczos` | thick-restart Lanczos 法 |
| `krylov_schur`、`arnoldi` | Krylov–Schur 法(`arnoldi` はその別名) |
| `shift_invert_lanczos`、`shift_invert_arnoldi` | shift-invert 型 |
| `dense_fallback_explicit` | 密行列の解法 |

前処理 M を渡す版では、shift-invert 型の内側の線形方程式の解法に M を
使います。shift を使わない `auto_select` では、前処理に対応していないことを
結果で返します。

### 尺度の正規化

`numeric_limits` が特殊化された実数型では、s = 2^round(log₂‖A‖∞) として
A/s の固有値を求め、結果を s 倍して返します(shift も σ/s にします)。2 の
べきで割るので、丸めの誤差は生じません。‖A‖∞ が 0 か有限でないときは
正規化しません。

## thick-restart Lanczos 法

実対称(または複素 Hermite)な作用素のための解法です。Wu–Simon の
thick-restart 法で、射影した行列を矢じり形(arrowhead)で持ち、収束した
Ritz 対を順に固定(locking)して、重複した固有値も見つけます。

- 開始ベクトルは決定的に作ります(`random_start = false` が既定)。
- 基底は完全に再直交化します。
- `max_iter` は行列とベクトルの積の総数の上限です。

収束(`converged = true`)とするのは、次をすべて満たすときです。

1. k 個の Ritz 対が、推定残差 ≤ tol · max(ε^(2/3) ‖A‖_est, |θ|) で固定された。
2. 終了時に、返す集合より目標側に、まだ収束していない Ritz 値が確実には
   ない(下記の「終了の判定」)。
3. 重複固有値についての条件を満たす。

結果の残差は診断用で、収束の判定には使いません。作用素が対称でないと、
結果は意味を持ちません。`auto_select` は、行列が対称と言えるときだけこの
解法を選びます。

## Krylov–Schur 法

一般の実行列のための解法です。射影した行列を実 Schur 形で持ち、隣り合う
ブロックの交換(LAPACK の dtrexc/dlaexc の形)で Schur 形を並べ替えて縮約
します。2×2 の複素ブロックは縮約の境界で分けません。交換を拒否したときは、
目標の候補を捨てずに、保持する範囲を広げます。

- `max_iter` は行列とベクトルの積の総数の厳密な上限です(最後の残差の
  評価も含みます)。上限に達したら `max_iter_exhausted`。
- 乱数は使わず、同じ入力には同じ出力を返します。

収束とするのは、次をすべて満たすときです。

1. 返す k 対すべての残差を作用素で計算し直して、残差 ≤ tol(絶対または
   相対)と言える。
2. まだ収束していない候補が、返す集合より目標側に確実にはない。
3. 返す範囲に複素共役対がない(`allow_complex_pairs = true` のときは、
   複素共役対も返せます)。
4. 判定に使った部分空間が、前回の判定の後に決定的な方向で拡げたもので
   ある(重複固有値を見落とさないため)。または基底が全空間を張る。

値の近さで重複を取り除くことはしません。

## shift-invert 型

σ に近い固有値を求めるために、(A − σI)⁻¹(一般化では (A − σB)⁻¹ B)に
Lanczos 型または Krylov–Schur 型の解法を使い、その固有値 μ から
λ = σ + 1/μ を求めます。

内側の線形方程式の解法は `shift_invert_solver` で選びます。

| 値 | 方法 |
| --- | --- |
| `sparse_lu`(既定) | A − σB を疎 LU 分解で 1 回分解し、毎回それで解く([tsparse_lu.md](tsparse_lu.md))。オプションは `shift_invert_lu` |
| `ilu0_gmres` | ILU(0) を前処理にした GMRES([tsparse_iterative.md](tsparse_iterative.md)) |

`sparse_lu` で分解に失敗したときは、`factorization_failed` を返します
(`ilu0_gmres` に自動では切り替えません)。

## 一般化固有値問題

`auto_select` は、入口で `lanczos` として扱います。

| 条件 | 解法 |
| --- | --- |
| A と B がともに対角 | 対角要素の比 a_ii / b_ii から直接求める |
| shift-invert 型を指定、または `use_shift = true` | 一般化 shift-invert。A・B が対称で Lanczos 型を指定したときは、B 内積の thick-restart Lanczos 法(ARPACK の mode 3 と同じ形。B は分解せず、半正定値でもよい)。それ以外は Krylov–Schur 法 |
| `lanczos` で `structure = symmetric` | B の疎 Cholesky 分解 B = Pᵀ L Lᵀ P を使い、L⁻¹ P A Pᵀ L⁻ᵀ に thick-restart Lanczos 法を使う |
| それ以外 | `unsupported_generalized_non_diagonal` |

尺度は、A を 2 のべき、B を 4 のべきで割って正規化します。

## 密行列の解法

- **対称**:Jacobi 法。`auto_select` から来たときは、回転の回数の上限を
  max(`max_iter`, 5n²) にします。
- **非対称**:2 のべきによる対角の平衡化 → Householder 変換による
  Hessenberg 化 → 実 Schur 分解 → 準上三角行列の後退代入による固有ベクトル。
  実 Schur 分解は LAPACK の dhseqr/dlahqr の形で、すべての副対角について
  |h_k,k−1| ≤ ε(|h_k−1,k−1| + |h_k,k|) と言えるときに分離します。反復の
  上限は 30m で、達したら失敗とします。

複素共役対は、`allow_complex_pairs = false`(既定)のときは収束した結果と
しては返さず、理由を付けて `converged = false` にします。

## 終了の判定

thick-restart Lanczos 法・shift-invert 型・Krylov–Schur 法は、収束を宣言する
直前に共通の判定(`tsparse_honest_termination.hpp`)を使います。返す予定の
集合の中で最も目標から遠い値よりも目標側に、まだ収束していない候補が
確実にあるときは、収束を宣言せずに反復を続けます。比較は「確実に成り立つ」
ときだけ真とするので、区間型では続ける側に倒れません。

`converged = true` は、返した固有対の残差が条件を満たし、目標側に見えて
いる未収束の候補がないことを意味します。求めた k 個が目標側の k 個すべてで
あることは保証しません。Krylov 部分空間法では、開始ベクトルと直交する
固有空間は見えないからです。

## オプション

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 解法 |
| `structure` | `auto_detect` | 行列の構造(`general`・`symmetric`・`hermitian`) |
| `target` | `smallest_algebraic` | 求める固有値(下記) |
| `tol` | 10⁻¹² | 収束の許容誤差 |
| `use_shift`、`shift` | `false`、0 | shift-invert 型を使う、その σ |
| `max_iter` | 1000 | 行列とベクトルの積の総数の上限 |
| `subspace_dim` | 0 | 部分空間の次元(0 は解法が決める) |
| `allow_dense_conversion`、`max_dense_size` | `true`、10⁶ | 密行列の解法を使ってよいか、その大きさの上限(n²) |
| `orthogonalization` | `modified_gram_schmidt` | 直交化の方法 |
| `random_seed`、`random_start` | 0、`false` | 開始ベクトル |
| `compute_residual_history` | `false` | 残差の履歴を記録する |
| `shift_invert_solver`、`shift_invert_lu` | `sparse_lu`、既定値 | shift-invert 型の内側の解法と、疎 LU のオプション |
| `allow_complex_pairs` | `false` | 複素共役対を返してよい |

| `target` | 求める固有値 |
| --- | --- |
| `largest_magnitude`、`smallest_magnitude` | 絶対値が大きい・小さい |
| `largest_algebraic`、`smallest_algebraic` | 値が大きい・小さい |
| `target_real` | `shift` に近い |
| `target_magnitude` | 絶対値が \|`shift`\| に近い |

## 結果

`eig_result<T>` の主なフィールドは次の通りです。

| フィールド | 内容 |
| --- | --- |
| `eigenvalues`、`eigenvectors` | 固有値と固有ベクトル |
| `eigenvalues_imag`、`complex_pair_count` | 複素共役対(`allow_complex_pairs = true` のときだけ) |
| `converged` | 収束したか(上記の意味) |
| `requested_count`、`returned_count`、`converged_count` | 求めた個数、返した個数、収束した個数 |
| `iterations`、`matrix_vector_products`、`linear_solves` | 反復回数、行列とベクトルの積の回数、線形方程式を解いた回数 |
| `residuals_absolute`、`residuals_relative` | 各固有対の残差 |
| `status`、`message`、`failure_reason` | 状態と理由(文字列) |
| `used_method` | 実際に使った解法(文字列) |
| `used_dense_fallback`、`used_shift_invert`、`used_generalized_operator` | 使った経路 |
