# 疎 LU 因子抽出(lu / lu_with_info)

`vcp::spmatrix<T>` に対する疎 LU 分解の**因子抽出** API の利用ガイドです。
自前の疎 LU(SLU)で分解し、L・U・置換(P/Q または p/q)を共通出力規約で
取り出します。本 API は**近似層**に属します: 因子は浮動小数点演算による
近似であり、検証(verified)を主張しません。

## 分解の規約(共通出力規約 §C)

正方行列 A に対し、行置換 P、列置換 Q、単位下三角行列 L(対角の 1 を
明示格納)、上三角行列 U(対角明示)を

```
P A Q = L U
```

を満たすように構成します。MATLAB の 4 出力形 `[L,U,P,Q] = lu(A)` と
同じ向きです。

- 置換ベクトル p, q は **new→old** 規約: `(PAQ)(i,j) = A(p[i], q[j])`、
  すなわち MATLAB 表記で `A(p,q) = L*U`。
- 置換行列は **P(k, p[k]) = 1**、**Q(q[k], k) = 1** です。
- **LDL との向きの違いに注意**: LDL(`docs/ldl.md`)は左因子が Pᵀ の分解
  `Pᵀ A P = L D Lᵀ` のため行列要素は `P(p[k], k) = 1` でした。LU は
  左因子が P そのもの(転置なし)なので、**行側の向きが LDL と転置の
  関係**になります(列側 Q は同じ向き)。
- 値が厳密零の要素は格納しません(spmats の不変条件)。
- L の各列は行番号昇順・先頭が対角 1、U の各列は行番号昇順です。

### MATLAB との対応表

| MATLAB | 本 API |
|---|---|
| `[L,U,P,Q] = lu(A)`(4 出力・行列形) | `A.lu(L, U, P, Q)` |
| `A(p,q) = L*U`(`'vector'` 形) | `A.lu(L, U, p, q)` |
| `P*A*Q - L*U` が微小 | `P*A*Q` と `L*U` の残差(受け入れテストで機械検証) |

## API(spmatrix)

```cpp
#include <vcp/spmatrix.hpp>

vcp::spmatrix<double> A, L, U, P, Q;
std::vector<vcp::spmatrix<double>::index_type> p, q;

// strict 版: status != success で vcp::throw_error(numerical_error)
A.lu(L, U, P, Q);        // 置換行列版
A.lu(L, U, p, q);        // 置換ベクトル版(MATLAB 'vector' 形式)

// 診断版(非 throwing)
auto r  = A.lu_with_info(L, U, P, Q /*, opt*/);
auto r2 = A.lu_with_info(L, U, p, q /*, opt*/);
```

- 結果型 `lu_extract_result`(= `spmatrix<T>::lu_extract_result_type`)は
  `status`(`vcp::sparse_lu_extract_status`)、`nnz_L` / `nnz_U`、
  `method_used`(実際に使われた SLU method)、`ordering_used`
  (**指定した** ordering の記録。auto_select の解決結果は SLU が
  再輸出しないため要求値のままです)を持ちます。
- L / U / p / q(P / Q)が有効な出力であるのは status が `success` の
  ときだけです。それ以外では out 引数は空のまま返ります。
- strict 版 `lu` は **status ≠ success のすべて**で `vcp::throw_error`
  します。特異行列等を扱う用途は `lu_with_info` を使ってください。

### オプション

`lu_extract_options<T>`(= `spmatrix<T>::lu_extract_options_type`)は
SLU のオプション `sparse_lu_options<T>` を 1 メンバ `slu` で包む薄い型です:

```cpp
vcp::spmatrix<double>::lu_extract_options_type opt;
opt.slu.ordering = vcp::sparse_lu_ordering::colamd;   // 例: 列順序付け
auto r = A.lu_with_info(L, U, p, q, opt);
```

## status の意味

| status | 意味 |
|---|---|
| `success` | 抽出成功(L/U/p/q 有効) |
| `invalid_factorization` | 分解自体が失敗(特異行列等) |
| `unsupported_options` | equilibration 済み因子(下記) |
| `unsupported_storage` | 未知の格納種別(supernodal は LUX-3 で対応済み) |
| `internal_error` | 内部エラー(発生したら報告してください) |

## equilibration 非対応(重要)

共通出力規約 `P A Q = L U` は**スケーリングを含みません**。
`opt.slu.equilibration = true` の指定は入口で `unsupported_options` として
honest に拒否されます(SLU の既定は false なので既定利用では発生しません)。
また `pivoting = static_mc64` も **GP 経路では**因子に Dr/Dc スケーリングを
持つため同様に `unsupported_options` になります。スケーリング込みの拡張形
`P·Dr·A·Dc·Q = L·U` の返却は将来の純増課題です。

## supernodal 対応(LUX-3)

`method = supernodal` の因子(`sparse_lu_storage_kind::supernodal`)にも
対応しています。追加 API はなく、同じ `A.lu(...)` / `A.lu_with_info(...)` /
tsparse 自由関数がそのまま通ります(policy/spmatrix 層は無変更)。

- 抽出のディスパッチは `fac.solve()` と同じ判定を鏡映します:
  - **受理済み A_eff 起源格納**(`true_numeric_source == true`): supernodal
    パネル / U_segments を平坦化して §C 規約の L/U/p/q を返します。
  - **遷移的格納**(`true_numeric_source == false`、A_eff 数値化が
    pivot_failure 等で不成立のとき): solve と同様に **baseline CSC 因子へ
    fallback** して抽出します(パネル値は数値源ではないため)。
- どちらの経路で抽出したかは診断フィールド
  `sparse_lu_extracted::storage_kind_extracted`(success 時のみ有意)で
  分かります。
- relaxed amalgamation / dense-front model がパネルに持つ**明示零は
  コピーされません**(baseline 抽出と同じ厳密零 drop 規則。nnz は
  baseline 抽出と同水準になります)。
- **native MC64 経路**(`static_mc64` + `supernodal_self_symbolic` +
  `supernodal_inplace_frontal`)は matching-only(Dr/Dc なし)なので
  §C 抽出**可能**です(GP 経路の static_mc64 と挙動が異なる点に注意)。
- LDL(`docs/ldl.md`)には supernodal 格納は存在しないため、
  LDL 側に対応する変更はありません。

## 因子消費層(lu_solve / lu_inverse_row)

抽出済み(または §C 規約で外部から用意した)因子 (L, U, p, q) を
**引数に取る**消費 API です。§C の共通規約により、因子の供給元が自前
(`A.lu`)でも将来の外部バックエンド(spumar / UMFPACK)でも同一の
消費関数がそのまま動きます。数理は A⁻¹ = Q·U⁻¹·L⁻¹·P:

```cpp
// strict 版(status != success で vcp::throw_error)
std::vector<double> x   = A.lu_solve(L, U, p, q, b);        // A x = b の近似解
std::vector<double> row = A.lu_inverse_row(L, U, p, q, i);  // A^{-1} の第 i 行

// 診断版(非 throwing)
std::vector<double> x2, row2;
vcp::lu_apply_result r1 = A.lu_solve_with_info(L, U, p, q, b, x2);
vcp::lu_apply_result r2 = A.lu_inverse_row_with_info(L, U, p, q, i, row2);
```

- 結果型 `lu_apply_result` は status のみ: `success` /
  `dimension_mismatch`(寸法・添字範囲)/ `singular_factor`
  (U 対角が certified に零または構造的欠落)/ `inconclusive_division`
  (U 対角の非零性を保証できない — 区間型で対角が 0 を跨ぐ場合等)/
  `invalid_input`(L が単位下三角でない・U に下三角要素・p/q が置換で
  ない)/ `internal_error`。
- U 対角による除算は certified 三分岐ゲート(GT1 P1)を通ります:
  **非零を保証できたときのみ除算**します。全順序型(double 等)では
  第三分岐(`inconclusive_division`)は到達不能です。
- L は単位下三角前提(入口で対角 = 1 の構造検査)なので除算はありません。
- `lu_inverse_row` の転置三角 solve は CSC 配列をそのまま行方向に読んで
  実装されており、転置行列は構築しません。
- 置換はベクトル形式 (p, q) のみ対応です(行列 P/Q 版は将来純増)。
- 実体はポリシー層の NVI ペア(`policy_lu_solve_with_info` /
  `policy_lu_inverse_row_with_info` + virtual `_impl`)で完結し、既定
  `_impl` は**引数の因子のみ**を消費します(ポリシー内部状態に非依存。
  任意の派生ポリシーがそのまま継承可能、`_impl` の override で差し替え可)。

## ポリシー層(拡張点)

**設計原則(WFIX-2)**: ポリシーは **\*this を対象とする**(密行列側
mats/pdblas と同一の関係)。かつての「対象行列 A を第 1 引数に渡す」
スタイルは廃止され、引数は他オペランド・出力・options のみです。

実体は `vcp::spmats<T,Index>` の NVI ペアです(対象行列 = *this):

```cpp
// ベクトル形(A.policy_lu_with_info(L, U, p, q, opt) の形で呼ぶ)
lu_extract_result<T,Index> policy_lu_with_info(
    spmats& L, spmats& U, std::vector<Index>& p, std::vector<Index>& q,
    const lu_extract_options<T>& opt) const;                 // 外殻(finalize + 正方性検査)
virtual lu_extract_result<T,Index> policy_lu_with_info_impl(
    spmats& L, spmats& U, std::vector<Index>& p, std::vector<Index>& q,
    const lu_extract_options<T>& opt) const;                 // 差し替え点

// 行列形(WFIX: 行列化そのものが差し替え点)
lu_extract_result<T,Index> policy_lu_with_info(
    spmats& L, spmats& U, spmats& P, spmats& Q,
    const lu_extract_options<T>& opt) const;                 // 外殻
virtual lu_extract_result<T,Index> policy_lu_matrices_with_info_impl(
    spmats& L, spmats& U, spmats& P, spmats& Q,
    const lu_extract_options<T>& opt) const;                 // 差し替え点
```

`policy_lu_with_info_impl` を override すると分解バックエンドを
差し替えられます(外部ライブラリ委譲: spumar / UMFPACK の予定差し替え点)。
外殻は override しないでください。

## T 汎用性

抽出層は T の算術・比較を持たない(値コピーと添字操作のみ。単位対角の
T(1) 構築と厳密零の非格納判定 `v == T(0)` を除く)ため、SLU が対応する
スカラー(float / double / long double / kv::dd / kv::mpfr&lt;N&gt; /
kv::interval&lt;…&gt; / std::complex&lt;double&gt;)をそのまま継承します。
Index は符号付き整数が必要です(unsigned Index は `vcp::state_error`)。
