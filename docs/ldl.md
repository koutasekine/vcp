# 疎 LDLᵀ 分解と慣性(inertia)API

`vcp::spmatrix<T>` に対する対称 LDLᵀ 分解(MATLAB 疎 `ldl` 互換の向き)と、
その利用者としての慣性(inertia)計算の利用ガイドです。
本 API は**近似層**に属します: 結果は浮動小数点演算による近似であり、
検証(verified)を主張しません。判定はすべて certified 比較
(「性質を保証できたときのみ成功宣言」)で行われ、保証できない場合は
制御されたステータス(`inconclusive_pivot_test` / `inconclusive_sign`)で
返ります。

## 分解の規約

対称行列 A に対し、置換行列 P、単位下三角行列 L(対角の 1 を明示格納)、
1×1/2×2 ブロック対角行列 D を

```
Pᵀ A P = L D Lᵀ
```

を満たすように構成します(数値法: Bunch–Kaufman 型 1×1/2×2 対称ピボット)。

- 入力 A は対称を前提とし、数値分解が読むのは**下三角部(対角含む)のみ**です
  (MATLAB `ldl(A)` の `'lower'` と同じ)。上三角は対称性検査でのみ参照します。
- 置換ベクトル p は **new→old** 規約: `(PᵀAP)(i,j) = A(p[i], p[j])`
  (MATLAB の `A(p,p)` と一致)。置換行列は `P(p[k], k) = 1` です。
- D の 2×2 ブロックは副対角を両側((k+1,k) と (k,k+1))に格納し、
  D 単体で対称行列として完結します。**値が厳密零の要素は格納しません**
  (spmats の不変条件)。特に零ピボットの対角は非格納(構造的零)です。

## API(spmatrix)

```cpp
#include <vcp/spmatrix.hpp>

vcp::spmatrix<double> A, L, D, P;
std::vector<vcp::spmatrix<double>::index_type> p;

// strict 版: status != success で vcp::throw_error(numerical_error)
A.ldl(L, D, P);          // 置換行列版
A.ldl(L, D, p);          // 置換ベクトル版(MATLAB 'vector' 形式)

// 診断版(非 throwing)
auto r  = A.ldl_with_info(L, D, P /*, opt*/);
auto r2 = A.ldl_with_info(L, D, p /*, opt*/);

// 慣性(n₊, n₋, n₀)
auto ir  = A.inertia();            // strict: status != success で throw
auto ir2 = A.inertia_with_info();  // 非 throwing
```

- strict 版 `ldl` は **status ≠ success のすべて**(`zero_pivot` /
  `not_symmetric` / `inconclusive_pivot_test` / `structural_singularity` を含む)
  で `vcp::throw_error` します。特異行列等を扱う用途は `ldl_with_info` を
  使ってください。
- L / D / P(p) が有効な出力であるのは status が `success` または
  `zero_pivot` のときだけです。それ以外では out 引数は空のまま返ります。
- `zero_pivot` は「certified に零のピボットをスキップして分解を最後まで
  続行した」ことを意味します(零固有値計数の用途)。残差の意味では
  LDLᵀ は A の該当行/列成分を再現しません(構造的零)。

## options 契約(`ldl_options<T>` = `sparse_ldl_options<T>`)

| フィールド | 既定値 | 意味 |
|---|---|---|
| `method` | `auto_select` | 現版は `baseline_dynamic`(動的 left-looking 疎核。小規模 n ≤ 64 は dense 核へ委譲し、診断 `dense_delegated` で報告)に解決 |
| `ordering` | `auto_select` | **契約**: auto_select はライブラリが選択する。現版は **amd** に解決。選択結果は診断 `ordering_used` で必ず報告される。**将来の版で選択は変わりうる**ため、再現性が必要な場合は明示指定すること。`natural` / `rcm` / `amd` / `nested_dissection` を明示指定可能(colamd は型レベルで不在) |
| `check_symmetry` | `true` | 入口の対称性検査。合格宣言側が certified(対称と**保証できない**入力は `not_symmetric`)。下三角のみを与える場合は `false` にする |
| `symmetry_tol` | 1e-12(実型) | certified に \|a_ij − a_ji\| ≤ tol で合格 |
| `zero_pivot_tol` | 0(実型) | 0 = certified 厳密零のみ零ピボット扱い |
| `pivot_threshold` | 0(実型) | **予約フィールド**。非 0(certified に 0 と確定できない値)は `invalid_options` |

tol 系フィールドは大きさ(`real_type<T>`、ADL `abs` の戻り型)で持ちます。
実数スカラーでは `real_type<T> = T` です。

## status 一覧(`sparse_ldl_status`)

| status | 意味 |
|---|---|
| `success` | 分解完了 |
| `structural_singularity` | 構造的に空の行/列(組み立てミスの示唆)。分解に入らず即返る。診断 `structural_empty_at` に最初の該当列 |
| `zero_pivot` | certified 零ピボット(分解は最後まで続行済み)。診断 `first_zero_pivot` |
| `inconclusive_pivot_test` | ピボット判定を certified に確定できず中断(第三分岐)。診断 `inconclusive_at`。区間スカラーで 0 を跨ぐ要素がある場合の正しい挙動 |
| `not_symmetric` | 対称性を保証できない入力 |
| `invalid_options` | 予約フィールド非 0、未実装 method 等 |
| `invalid_input` | n<0、CSC 不整合等 |
| `internal_error` | 実装バグ・想定外状態(P3 例外網の発火を含む。発火は不具合として扱う) |

## スカラー型 T への要求

mats<T> と同一の演算子集合(四則・比較・ADL の abs/sqrt 等)のみです。
型による分岐・拒否は行いません。`kv::interval<TT>` を含む全対象型で
コンパイル・実行可能ですが、**数値的成功は保証されません**: 比較が
certified に確定できないスカラーでは `inconclusive_pivot_test` /
`inconclusive_sign` の制御されたステータスで返ります(「正しく死ぬ」)。
全順序型(double / kv::dd / kv::mpfr<N>)では第三分岐は到達不能であり、
従来の BK と同一の挙動です。complex Hermitian(LDLᴴ)はスコープ外です
(コンパイルは可能)。

## 慣性 API

`inertia_with_info` の既定実装は「`ldl` を呼び、得た D をブロック対角として
certified 走査する」ものです(`inertia_options<T>` = `{ real_type zero_tol;
ldl_options<T> ldl; }`。`zero_tol` は慣性の零判定用で、LDL の
`zero_pivot_tol` と独立)。

- 内部 LDL が `zero_pivot` の場合は失敗ではありません(零固有値は D の
  構造的零として n₀ に計上され、status = `success`)。
- それ以外の非 success は `factorization_failed`(内部 status を診断
  `ldl_status` に保持)。
- 計上(成功宣言)は certified 比較が成立したときのみ行われ、確定できない
  ケースは `inconclusive_sign`(位置は診断 `inconclusive_at`)で返ります。
- 検算 n₊ + n₋ + n₀ = n を内部検査します(不一致は `internal_error`)。
- D 消費層 `policy_inertia_from_block_diagonal`(spmats の非 virtual
  ポリシーメンバ)は BK 由来を仮定せず、一般の 1×1/2×2 ブロック対角行列を
  受け付けます(帯域 1 超の非零は `not_block_diagonal`)。走査対象 D が
  subject であり、`D.policy_inertia_from_block_diagonal(tol)` の形で
  呼びます(WFIX-2)。

## 設計原則: ポリシーは *this を対象とする(WFIX-2)

ポリシー(spmats とその派生)は**データを保持する主体**であり、密行列側
(mats/pdblas/matrix)と同一の関係を取ります: `inv()` / `Cholesky()` が
*this を対象とするのと同じく、疎行列側のポリシーメソッドも対象行列は
常に *this です。かつての「対象行列 A を第 1 引数に渡す」スタイル
(`policy_ldl_with_info(A, L, D, p, opt)` 等)は廃止されました。
現行署名(対象 = *this):

```cpp
ldl_result<T,Index> policy_ldl_with_info(
    spmats& L, spmats& D, std::vector<Index>& perm,
    const ldl_options<T>& opt) const;            // 外殻(A.policy_ldl_with_info(L,D,p,opt) の形)
virtual ldl_result<T,Index> policy_ldl_with_info_impl(...) const;   // 差し替え点
inertia_result<Index> policy_inertia_with_info(
    const inertia_options<T>& opt) const;        // 外殻 + virtual _impl
```
- シフト付き `inertia(A − σI)` により σ を跨ぐ固有値計数ができます
  (Sylvester の慣性則)。

## MATLAB 疎 ldl との既知の差異

| 項目 | MATLAB(MA57 系) | 本実装 | 理由 |
|---|---|---|---|
| D の零ピボット対角 | 明示零を保持 | 非格納(構造的零) | spmats の不変条件(明示零禁止)。inertia は構造走査で対応 |
| ピボット | 閾値付き(THRESH、既定 0.01) | 正統 BK(α = (1+√17)/8 固定)+ certified 三分岐 | 数値的信頼性優先。閾値は予約フィールドのみ |
| 数値核 | multifrontal | dynamic left-looking(Phase A) | 段階戦略。multifrontal は Phase B(別設計)で追随余地 |

## 実装ノート

- tsparse 層の自由関数 `vcp::sparse_ldl_factorize_with_info`(生 CSC
  インターフェース)も利用できます(`<vcp/tsparse/tsparse_sparse_ldl.hpp>`)。
- ordering は既存 SLU の pattern-only 基盤(整数のみ、A+Aᵀ グラフ)を
  include 再利用しています。SLU 側の意味論は不変です。
- 疎核は dense 核(恒久回帰基準)とピボット系列が完全一致するよう
  構成されています(差分検証は `sandbox/tests/ldl2_test_1.cpp`)。
