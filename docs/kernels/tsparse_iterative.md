# tsparse — 作用素・前処理・反復法

この文書では、疎行列の反復法と、その部品を説明します。

| 内容 | 場所 |
| --- | --- |
| `spmats` の線形方程式の反復法(Jacobi・Gauss–Seidel・CG・BiCGSTAB・GMRES) | `vcp/spmats_base/spmats_lss.hpp` |
| 残差の管理 | `vcp/tsparse/tsparse_solvers.hpp` |
| 作用素の抽象型 | `vcp/tsparse/tsparse_operator.hpp` |
| 前処理 | `vcp/tsparse/tsparse_preconditioner.hpp` |
| ILU(0) 分解と GMRES(shift-invert の内側の解法) | `vcp/tsparse/tsparse_factorization.hpp` |

利用者が直接呼ぶことは想定していません。使い方は [spmatrix.md](../spmatrix.md)
の「線形方程式」を参照してください。直接法(疎 LU 分解)は
[tsparse_lu.md](tsparse_lu.md) で説明します。

反復法(Jacobi・Gauss–Seidel・CG・BiCGSTAB・GMRES)と前処理は、まだ実験的な
実装です。

## `spmats` の線形方程式の振り分け

`linear_solve_options<T>::method` で方法を選びます。

| 値 | 方法 |
| --- | --- |
| `auto_select`(既定) | `sparse_lu` |
| `sparse_lu` | 疎 LU 分解([tsparse_lu.md](tsparse_lu.md)) |
| `jacobi`、`gauss_seidel` | Jacobi 法、Gauss–Seidel 法 |
| `conjugate_gradient` | 共役勾配法 |
| `bicgstab` | BiCGSTAB 法 |
| `gmres` | リスタート付き GMRES 法 |

反復法の本体は `spmats_lss.hpp` にあり、どれも初期値 x₀ = 0 から始めます。
`vcp::spumar` / `vcp::spcmodumar` では、指定した方法によらず UMFPACK を
使います([sparse_external.md](sparse_external.md))。

### オプション

| フィールド | 既定値 | 意味 |
| --- | --- | --- |
| `method` | `auto_select` | 方法 |
| `max_iter` | 1000 | 反復回数の上限。0 は `invalid_argument` |
| `tol` | 10⁻¹² | 許容誤差。正でなければ `invalid_argument` |
| `use_relative_residual` | `true` | 相対残差で判定する |
| `check_symmetric` | `true` | CG で対称性を確かめる |
| `restart` | 30 | GMRES のリスタートの周期。0 は `invalid_argument` |
| `preconditioner` | `none` | 前処理(`none`・`jacobi`) |
| `sparse_lu` | 既定値 | 疎 LU のオプション |

### 収束の判定

収束の閾値は `tsparse_solvers::make_residual_control` で決めます。

- `use_relative_residual = true`:閾値 = tol · max(‖b‖₂, 1)
- `use_relative_residual = false`:閾値 = tol

結果の残差のフィールド(`absolute_residual_norm`、`relative_residual_norm`、
`initial_residual_norm`)は、最後に真の残差 ‖b − A x‖₂ を計算し直して
入れます。`relative_residual_norm` は ‖b − A x‖₂ / max(‖b‖₂, 1) です。
これらの値が有効なのは `converged = true` のときだけです。

### 各方法

| 方法 | 内容 | 収束の判定に使う残差 |
| --- | --- | --- |
| Jacobi | x_i ← (b_i − Σ_{j≠i} a_ij x_j) / a_ii を全成分同時に更新 | 各反復の後の真の残差 |
| Gauss–Seidel | 同じ式を、更新した成分をすぐ使って i の順に更新 | 各反復の後の真の残差 |
| CG | 前処理付き共役勾配法 | 漸化式で更新した残差 r |
| BiCGSTAB | 左から前処理した BiCGSTAB 法 | 真の残差 |
| GMRES | 左から前処理した GMRES(m) 法。Givens 回転で最小二乗問題を解く | 射影した残差が閾値を下回ったら、解を作って真の残差で判定する |

- Jacobi・Gauss–Seidel は、対角に厳密な零があると `vcp::numerical_error` を
  投げます。前処理は使いません。
- CG は、`check_symmetric = true` のとき、`is_symmetric(10⁻¹⁰)` でなければ
  `vcp::domain_error` を投げます。pᵀ A p または rᵀ z が正で有限と言えなく
  なったら、収束しないまま止めます。
- BiCGSTAB は、次のいずれかが成り立つと言えなくなったら(breakdown)、
  収束しないまま止めます。どれも尺度によらない形の判定です。
  - |ρ| > ε ‖r̂‖ ‖r‖(ρ = r̂ᵀ r)
  - |r̂ᵀ v| > ε ‖r̂‖ ‖v‖
  - tᵀ t > 0
  - |ω| ‖t‖ > ε ‖s‖

  また、値が有限でなくなったときも止めます。
- GMRES の `max_iter` は、内側の反復の総数の上限です。射影した残差は前処理を
  かけた残差のノルムで、途中の判定の閾値は前処理をかけた後の右辺から決めます。
  途中の判定が通ったら解を作り、`converged` は真の残差と、前処理をかける前の
  b から決めた閾値で決めます。真の残差が閾値を超えていたら、その解から
  リスタートします。
- BiCGSTAB の半歩の判定も同じです。前処理をかけた残差が、前処理をかけた後の
  右辺から決めた閾値を下回ったら、半歩の解を作って真の残差を調べます。
  閾値以下ならその解で終わり、超えていたら反復を続けます。

前処理 `jacobi` は、対角要素の逆数を掛けます。対角に厳密な零があると
`vcp::numerical_error` を投げます。

## 作用素

`tsparse_operator.hpp` には、y = A x を計算する作用素の型があります。

| 型 | 内容 |
| --- | --- |
| `linear_operator<T, Index>` | `rowsize()`・`columnsize()`・`apply(x, y)` を持つ抽象型 |
| `sparse_matrix_operator<SparseMatrix>` | 疎行列を参照して、`apply` で `mul_vec` を呼ぶ |
| `function_linear_operator<T, Index>` | 関数オブジェクト `f(x, y)` を作用素として包む |

## 前処理

`tsparse_preconditioner.hpp` には、z = M⁻¹ r を計算する前処理の型があります。
どれも `apply(r, z)` を持ちます。固有値問題の前処理付きの版に渡すことが
できます([tsparse_eigs.md](tsparse_eigs.md))。

| 型 | 内容 |
| --- | --- |
| `identity_preconditioner<T>` | z = r |
| `function_preconditioner<T, Functor>` | 関数オブジェクト `f(r, z)` を包む |
| `jacobi_preconditioner<SparseMatrix>` | z_i = r_i / a_ii |
| `ilu0_preconditioner<SparseMatrix>` | ILU(0) 分解による前処理 |

`jacobi_preconditioner` と `ilu0_preconditioner` は、作るときに例外を投げず、
`valid()`・`zero_pivots()`・`diagnostics()` で状態を返します。対角(ピボット)の
絶対値が 10⁻¹⁴ 以下のときは無効になり、無効なまま `apply` を呼ぶと
`vcp::numerical_error` を投げます。

## ILU(0) 分解

`tsparse_factorization::ilu0_factorize` は、A の非零の位置だけを使う
不完全 LU 分解です(fill-in を作りません)。CSR の各行 i について、
k < i の要素 a_ik ごとに

- l_ik = a_ik / u_kk
- k より右の a_kj について、(i, j) が A の非零の位置にあれば a_ij ← a_ij − l_ik u_kj

とします。L(単位対角、対角は持たない)と U を、A と同じ位置の 1 つの配列に
持ちます。ピボット u_kk の絶対値が 10⁻¹⁴ 以下のとき、または対角の要素が
ないときは、その位置を記録して `singular_or_unstable` にします。

`ilu0_solve` は、L z = b の前進代入と U x = z の後退代入で解きます。

## GMRES(shift-invert の内側の解法)

`tsparse_factorization::gmres_solve(apply_A, apply_prec, b, max_iter, tol,
restart, relative)` は、`spmatrix` に依存しない GMRES(m) 法です。固有値問題の
shift-invert 型で `shift_invert_solver = ilu0_gmres` としたとき、内側の線形
方程式を ILU(0) を前処理にして解くのに使います。

- 初期値は 0、m = min(`restart`, n) です。`max_iter` は Arnoldi の 1 段の
  総数の上限です(リスタートをまたいで数えます)。
- 左から前処理し、前処理した残差 M⁻¹(b − A x) のノルム β で判定します。
  `relative = true` のときは、閾値を tol · β₀(β₀ は最初の β)にします。
- β が有限でないとき、または真の残差が零でないのに前処理した残差が零に
  なったときは、収束しないまま止めます。
