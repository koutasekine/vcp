# 疎行列の内部構成

この文書では、疎行列クラス `vcp::spmatrix` の内部で、各 policy がどの
下請けに計算を任せているかを説明します。利用者が直接呼ぶことは想定して
いません。使い方は [spmatrix.md](../spmatrix.md) を参照してください。

## 3 つの層

疎行列の計算は、次の 3 つの層に分かれています。

| 層 | 場所 | 役割 |
| --- | --- | --- |
| 行列クラス | `vcp/spmatrix.hpp` | 利用者向けの API。計算は policy に任せる |
| policy | `vcp/spmats.hpp`、`vcp/spumar.hpp`、`vcp/spcmodumar.hpp`、`vcp/spimats.hpp` と、それぞれの `*_base/` | 格納と、演算・分解・固有値計算の振り分け |
| 下請け | `vcp/tsparse/`、外部ライブラリ(UMFPACK・ARPACK・CHOLMOD) | 疎行列のカーネル |

policy の下請けのファイルは、次のディレクトリにあります。

| ディレクトリ | 使う policy | 内容 |
| --- | --- | --- |
| `vcp/spmats_base/` | `spmats` | 演算・線形方程式・分解・固有値計算の policy 側の実装、近似逆行列(AINV・FSAI) |
| `vcp/spumar_base/` | `spumar` | UMFPACK・ARPACK の呼び出しと、格納形式の変換 |
| `vcp/spcmodumar_base/` | `spcmodumar` | CHOLMOD の呼び出し |
| `vcp/spimats_base/` | `spimats` | 区間行列と点行列の変換、点と区間の混合演算 |
| `vcp/tsparse/` | `spmats`(と、それを継承する policy) | 疎行列のカーネル。[tsparse.md](tsparse.md) |

## policy の継承関係

```
vcp::spmats<T, Index>
 ├─ vcp::spumar                      (spmats<double> を継承)
 │   └─ vcp::spcmodumar             (spumar を継承)
 └─ vcp::spimats<T, AP, DP>          (= spmats<kv::interval<T>, AP::index_type> を継承)
```

`vcp::spmatrix<T, P>` は policy `P` を protected 継承しています。
`spmatrix` の関数は policy の非 virtual な関数 `policy_X` を呼び、
`policy_X` は入口の検査をしてから virtual な関数 `policy_X_impl` を
呼びます。派生した policy は、`policy_X_impl` を上書きすることで計算の
下請けを差し替えます。上書きしない関数は、基底の `spmats` の実装を
そのまま使います。

各 policy が上書きしている `policy_X_impl` は次の通りです。

| policy | 上書きしている関数 |
| --- | --- |
| `spumar` | `policy_lss_with_info_impl`、`policy_lss_multi_with_info_impl`、`policy_lu_with_info_impl`、`policy_eigs_with_info_impl`、`policy_generalized_eigs_with_info_impl` |
| `spcmodumar` | `policy_chol_with_info_impl` |
| `spimats` | `policy_lss_with_info_impl`、`policy_lss_multi_with_info_impl`、`policy_eigs_with_info_impl`、`policy_generalized_eigs_with_info_impl`、`policy_is_spd_impl` |

## 操作ごとの下請け

`spmatrix` の操作ごとに、各 policy で使われる下請けを示します。
「継承」は、基底の policy と同じ下請けを使うことを表します。

| 操作 | `spmats` | `spumar` | `spcmodumar` | `spimats` |
| --- | --- | --- | --- | --- |
| 和・差・積 | `spmats_product.hpp` → `tsparse_spgemm` | 継承 | 継承 | 継承(区間型の要素で同じ経路)。点と区間の混合演算は `spimats_base/spimats_mixed_ops.hpp` |
| 行列とベクトルの積 | `spmats.hpp` → `tcsrmv` / `tcscmv` | 継承 | 継承 | 継承 |
| 線形方程式(直接法) | `spmats_lss.hpp` → `sparse_lu_factorize_with_info`。反復改良を指定したときは `sparse_lu_solve_refined` | UMFPACK(`spumar_base/spumar_umfpack.hpp`) | 継承(UMFPACK) | 精度保証付き。近似解と近似逆行列の行を AP の LU 分解で求め、区間演算で包含を検証する(`spimats.hpp`) |
| 線形方程式(反復法) | `spmats_lss.hpp` の Jacobi・Gauss–Seidel・CG・BiCGSTAB・GMRES | 指定した方法によらず UMFPACK | 継承(UMFPACK) | 上と同じ |
| LU 分解(因子を返す) | `spmats_lu_extract_impl.hpp` → `sparse_lu_factorize_with_info` と `sparse_lu_extract_factors` | UMFPACK | 継承(UMFPACK) | 継承(区間型の要素で計算。精度保証なし) |
| LU 因子による解・逆行列の行 | `spmats_lu_extract_impl.hpp`(三角行列の解法) | 継承 | 継承 | 継承 |
| LU 分解の handle | `spmats_lss.hpp` → `sparse_lu_factorize_with_info` | 継承 | 継承 | 継承 |
| LDLᵀ 分解・慣性 | `spmats_ldl_impl.hpp` → `sparse_ldl_factorize_with_info` | 継承 | 継承 | 継承 |
| shift 付き LDLᵀ の handle | `spmats_ldl_shift.hpp` → `sparse_ldl_symbolic_analyze` と `sparse_ldl_factorize_numeric_with_info` | 継承 | 継承 | 継承 |
| Cholesky 分解 | `spmats_chol_impl.hpp` → `sparse_chol_factorize_with_info` | 継承 | CHOLMOD(`spcmodumar_base/spcmodumar_cholmod.hpp`) | 継承 |
| 正定値の判定 `is_spd` | Cholesky 分解が成功するか | 継承 | 継承(CHOLMOD の Cholesky 分解) | 精度保証付き。点の `spmats<T>` で近似 Cholesky 分解をし、区間演算で十分条件を検証する(`spimats.hpp`) |
| 近似逆行列(AINV・FSAI) | `spmats_base/spmats_ainv_impl.hpp`、`spmats_fsai_impl.hpp`、`spmats_fsai_adaptive_impl.hpp` | 継承 | 継承 | 継承 |
| 固有値(標準・一般化) | `spmats_eigs.hpp` → tsparse の固有値解法(Lanczos 系、Krylov–Schur、shift-invert、小さい行列は密行列の解法) | ARPACK(shift-invert の分解は UMFPACK)。ARPACK で扱わない指定は `spmats` の実装に戻す | 継承 | 精度保証付き。AP で近似固有値を求め、AP の shift 付き LDLᵀ の慣性と区間演算で包含を検証する(`spimats.hpp`) |

## `spimats` の AP と DP

`vcp::spimats<T, AP, DP>` は、精度保証の計算の途中で次の 2 つの policy を
使います。

- `AP`(既定は `vcp::spmats<T>`):点の疎行列 policy。区間行列の中点を
  `spimats_base/spimats_convert.hpp` の `mid_for_approximation` で取り出し、
  近似解・近似固有値・LU 因子・shift 付き LDLᵀ をこの policy で計算します。
  `AP = vcp::spumar` とすると、これらの計算に UMFPACK・ARPACK が使われます。
- `DP`(既定は `vcp::imats<T>`):区間の密行列 policy。固有値の包含の
  検証で、小さな区間密行列の一般化固有値問題(`eigsymge`)に使います。

近似の計算の結果は包含の成否と幅にだけ影響し、検証は区間演算で行います。
`spimats` の内部では、丸めモードを自分で切り替えません。丸めは kv の区間
演算と `vcp::convert` に任せます。

## 外部ライブラリとの境界

UMFPACK・ARPACK・CHOLMOD のヘッダを include するのは、`spumar_base/` と
`spcmodumar_base/` のファイルだけです。`spmats.hpp` などの本体のヘッダは
これらを include しないので、`spumar` / `spcmodumar` を使わない限り外部
ライブラリは不要です。

## 下請けの文書

| 文書 | 内容 |
| --- | --- |
| [tsparse.md](tsparse.md) | tsparse の基本(格納形式、変換、積、要素型の要件) |
| [tsparse_lu.md](tsparse_lu.md) | 疎 LU 分解 |
| [tsparse_ldl.md](tsparse_ldl.md) | 疎 LDLᵀ 分解・Cholesky 分解・慣性 |
| [tsparse_eigs.md](tsparse_eigs.md) | 疎行列の固有値解法 |
| [tsparse_iterative.md](tsparse_iterative.md) | 作用素・前処理・反復法 |
| [sparse_approx_inverse.md](sparse_approx_inverse.md) | 近似逆行列(AINV・FSAI) |
| [spimats_kernel.md](spimats_kernel.md) | `spimats` の精度保証の手順と混合演算 |
| [sparse_external.md](sparse_external.md) | UMFPACK・ARPACK・CHOLMOD への委譲 |
