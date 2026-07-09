# bfem 概要

`vcp/bfem` は Bernstein 多項式を基底にした有限要素法用のヘッダ群です。
三角形・四面体メッシュ上で、連続 Lagrange 型の `P^k`、Raviart-Thomas 型の
`RT^k`、broken `P_l`、Scott-Vogelius 用のベクトル値空間、2 次元 C^1
要素を扱います。

このライブラリの特徴は、要素積分に数値積分を使わないことです。Bernstein
基底の積・次数上げ・質量行列などを有理数テーブルとして構成し、スカラー型 `T`
へ一度だけ変換して使います。`T` に区間型を使えば、行列・ベクトル・スカラー量は
厳密値の包含として計算できます。

## できること

| 分野 | 主なヘッダ | 主な機能 |
|---|---|---|
| 共通メッシュ・幾何 | `<vcp/bfem/mesh.hpp>`, `<vcp/bfem/geometry.hpp>` | 2D/3D 単体メッシュ、要素測度、重心座標勾配、退化要素検出 |
| Bernstein 多項式 | `<vcp/bfem/bpoly.hpp>`, `<vcp/bfem/poly1.hpp>`, `<vcp/bfem/refine.hpp>` | 多項式の和・積・次数上げ・微分・合成・評価・範囲包含 |
| 連続 `P^k` | `<vcp/bfem/fe_space.hpp>` | 剛性、質量、荷重、重み付き質量、対流、内積、dual、次数上げ |
| 斉次 Dirichlet | `<vcp/bfem/dirichlet.hpp>` | 境界自由度の削除、縮小系から全自由度への復元 |
| `RT^k` と broken `P_l` | `<vcp/bfem/rt/rt_space.hpp>`, `<vcp/bfem/rt/broken_space.hpp>`, `<vcp/bfem/rt/rt_assemble.hpp>` | RT 質量、div 質量、cross grad、broken load、投影誤差、flux 誤差 |
| Scott-Vogelius 部品 | `<vcp/bfem/sv/*.hpp>` | ベクトル値 `P^k`、div 行列、移流項、圧力拘束、Alfeld 分割、線形縮約 |
| 2D C^1 要素 | `<vcp/bfem/c1/c1_space.hpp>`, `<vcp/bfem/c1/c1_reduce.hpp>` | C^1 行列、補間、勾配・Hessian 評価、Laplacian 残差、境界拘束 |

## 計算の考え方

`bfem` は「要素ごとの厳密な多項式演算」と「大域自由度への決定的な散布」を行います。

bfem の説明では、話題のスコープを次の 4 種類に分けると読みやすくなります。

| スコープ | 意味 | 代表例 |
|---|---|---|
| 領域非依存 | メッシュや物理要素をまだ見ない、純粋な添字・有理数・多項式テーブルの話 | `rational`, `index_map`, `coeff_registry`, `typed_registry`, `poly1` |
| 参照単体上 | 標準三角形・標準四面体上の Bernstein 多項式の話。物理的な面積・体積はまだ掛からない | `bpoly`, `inner(bpoly,bpoly)`, `eval`, `restrict_to`, `range` |
| 物理要素上 | メッシュ中の 1 要素へ写した後の局所計算。測度、Jacobian、重心座標勾配が入る | `element_geometry`, `element_op`, `rt_element_op`, `c1_element_op` |
| 領域全体 | 全要素を回って大域自由度へ散布した後の行列・ベクトル・場の話 | `fe_space`, `rt_space`, `broken_space`, `vfe_space`, `c1_space`, `dirichlet_reduction` |

1. `mesh<D,T>` が頂点座標と要素頂点番号を保持します。
2. `element_geometry<D,T>` が各要素の測度、向き、重心座標勾配を計算します。
3. `bpoly<D,T>` が参照単体上の Bernstein 係数列を保持します。
4. `element_op<D,T,P>` や `rt_element_op`、`c1_element_op` が局所行列を作ります。
5. `fe_space`、`rt_space`、`broken_space`、`c1_space` などが大域行列・ベクトルへ組み立てます。

要素順と重複成分の結合順は固定されています。同じ入力なら同じ順序で加算されます。

## OpenMP 並列化

bfem の大域行列組立の一部は OpenMP に対応しています。利用者側では通常、
プログラムを `-fopenmp` 付きでコンパイルすれば有効になります。スレッド数を
プログラム内で制御したい場合だけ、利用者コードで `<omp.h>` を include して
`omp_set_num_threads` などを呼びます。

bfem 側の OpenMP だけを止める場合は `-DVCP_BFEM_NOMP` を付けます。
プロジェクト全体で `-DVCP_NOMP` を使う場合も、bfem では OpenMP 無効として扱われます。

並列化される主な対象は、`fe_space` の `stiffness` / `mixed_mass` /
`weighted_mass`、RT や broken 空間の大域行列、Scott-Vogelius の vector stiffness
や div/advection 系、2D C^1 の stiffness / mass / Laplacian / Hessian 系です。
詳細な使い方は `docs/bfem_user_guide.md` の OpenMP 節を参照してください。

## スカラー型 `T`

通常は `double` で近似計算できます。検証付き計算では `kv::interval<double>` などの
区間型を使います。区間型で使うときは、用途に応じて以下も include します。

```cpp
#include <vcp/bfem/bound_traits_kv.hpp>
#include <vcp/bfem/geometry_traits_kv.hpp>
#include <vcp/bfem/convert_traits_kv.hpp>
```

`T` には基本的に `T(0)`, `T(1)`, 加減乗除、コピー、比較または専用 traits が必要です。
幾何構築時には要素ごとに 1 回だけ除算が行われます。2D C^1 要素では、要素ごとの除算に加えて
辺長に関する除算も使います。

## 疎行列ポリシー

大域疎行列は `vcp::spmatrix<T, SP>` で返ります。テンプレート引数 `SP` の既定値は
`vcp::spmats<T>` です。

区間型 `T` で疎行列まで作る場合、既定の `spmats` がそのまま使えないことがあります。
その場合は、区間型に対応した sparse policy を `fe_space<D,T,P,SP>` などへ渡します。
疎行列を作らない API だけを呼ぶ場合、`SP` は実体化されないため既定のまま使える場面があります。

## 責務の境界

`bfem` は有限要素空間と行列・ベクトル・多項式量を作る部品です。線形方程式の求解、
Newton 反復、固有値計算、検証定理の適用、メッシュ生成、非多項式関数との数値求積は
利用者側の責務です。

特に Scott-Vogelius 関連は「安定化済みの完全な Navier-Stokes ソルバ」ではありません。
速度空間、圧力空間、div 行列、移流項、特異拘束、縮約などを提供する部品群です。

## ファイル構成と役割

### 共通・基礎層

| ファイル | スコープ | 役割 |
|---|---|---|
| `mesh.hpp` | 領域全体の入力 | 頂点座標と要素接続を保持する薄い値型 |
| `multi_index.hpp` | 領域非依存 | Bernstein 係数の多重添字、rank/unrank |
| `rational.hpp` | 領域非依存 | bigint と rational、有理数テーブル構築用 |
| `convert_traits.hpp` | 領域非依存 | 有理数から `T` への変換点 |
| `convert_traits_kv.hpp` | 領域非依存 | kv 型向け変換補助 |
| `bound_traits_kv.hpp` | 領域非依存 | kv interval の範囲包含用 traits |
| `geometry_traits_kv.hpp` | 物理要素上 | kv interval の幾何符号判定用 traits |
| `detail/scalar_traits.hpp` | 領域非依存 | `T` が必要な演算を満たすかの静的検査 |
| `detail/table_cache.hpp` | 領域非依存 | テーブルキャッシュ共通部品 |
| `detail/rational_la.hpp` | 領域非依存 | rational 行列の小規模厳密線形代数 |

### Bernstein 多項式・要素計算

| ファイル | スコープ | 役割 |
|---|---|---|
| `coeff_tables.hpp` | 領域非依存 | 有理数の質量・次数上げ・積・微分テーブル |
| `typed_tables.hpp` | 領域非依存 | 有理数テーブルを `T` へ変換したキャッシュ |
| `bpoly.hpp` | 参照単体上 | 参照単体上の Bernstein 多項式 |
| `poly1.hpp` | 領域非依存/参照単体上 | 1 変数多項式と `f(u)` の合成 |
| `refine.hpp` | 参照単体上 | 部分単体への制限、範囲包含、2D red refinement |
| `d3/bey_table.hpp` | 参照単体上 | 3D Bey refinement と 3D 範囲精密化 |
| `geometry.hpp` | 物理要素上 | 要素幾何、測度、重心座標勾配、退化検出 |
| `element_op.hpp` | 物理要素上 | `P^k` の局所剛性・質量・荷重・対流・内積 |
| `ref_stiffness.hpp` | 領域非依存/参照単体上 | 参照剛性テンソルのキャッシュ |

### 連続 `P^k` 空間

| ファイル | スコープ | 役割 |
|---|---|---|
| `dofmap.hpp` | 領域全体 | 2D `P^k` の大域自由度番号と境界自由度 |
| `d3/topology3.hpp` | 領域全体 | 3D メッシュの辺・面トポロジ |
| `d3/dofmap3.hpp` | 領域全体 | 3D `P^k` の大域自由度番号と境界面自由度 |
| `d3/trace3.hpp` | 参照面/参照辺 | 3D 面・辺トレースの添字補助 |
| `d3/s3_perm.hpp` | 領域非依存 | 3D 面・辺向きの置換補助 |
| `fe_function.hpp` | 領域全体 | `fe_function` と `dual_vector` |
| `fe_space.hpp` | 領域全体 | 2D/3D `P^k` 空間と大域アセンブル |
| `fe_space.hpp` 内 `detail::coo_buffer` | 領域全体 | COO triplet の決定的結合 |
| `dirichlet.hpp` | 領域全体 | 斉次 Dirichlet 縮約 |

### RT・broken 空間

| ファイル | スコープ | 役割 |
|---|---|---|
| `rt/rt_tables.hpp` | 領域非依存/参照単体上 | RT 参照基底・div・flux・質量系の有理数テーブル |
| `rt/rt_typed_tables.hpp` | 領域非依存 | RT テーブルの `T` 版キャッシュ |
| `rt/rt_backend3.hpp` | 領域全体/参照面 | 3D RT の面自由度、向き、補間補助 |
| `rt/rt_element_op.hpp` | 物理要素上 | RT 局所演算 |
| `rt/rt_space.hpp` | 領域全体 | 2D/3D `RT^k` 空間、RT 質量、補間 |
| `rt/broken_space.hpp` | 領域全体 | 要素ごとに不連続な broken `P_l` 空間 |
| `rt/rt_assemble.hpp` | 領域全体 | RT、broken、`P^k` 間の大域行列・スカラー |

### Scott-Vogelius

| ファイル | スコープ | 役割 |
|---|---|---|
| `sv/alfeld.hpp` | 領域全体 | Alfeld barycentric refinement |
| `sv/vfe_space.hpp` | 領域全体 | ベクトル値 `P^k` 空間ラッパ |
| `sv/sv_assemble.hpp` | 領域全体 | ベクトル剛性、div、移流、div ノルム |
| `sv/sv_rows.hpp` | 領域全体 | 特異頂点・特異辺に由来する圧力拘束行 |
| `sv/sv_singular.hpp` | 領域全体/幾何判定 | 特異性検出の実装 |
| `sv/sv_constraint.hpp` | 領域全体 | 整数係数の拘束行型 |
| `sv/linear_reduction.hpp` | 領域全体 | 一般線形拘束による除算フリー縮約 |

### 2D C^1 要素

| ファイル | スコープ | 役割 |
|---|---|---|
| `c1/c1_tables.hpp` | 領域非依存/参照単体上 | C^1 基底・Hermite・Laplacian 系テーブル |
| `c1/c1_geometry.hpp` | 物理要素上 | C^1 の pullback と辺長関連量 |
| `c1/c1_dofmap.hpp` | 領域全体 | C^1 自由度番号 |
| `c1/c1_element_op.hpp` | 物理要素上 | C^1 局所行列・微分演算 |
| `c1/c1_space.hpp` | 領域全体 | C^1 空間、大域行列、補間、Laplacian 残差 |
| `c1/c1_reduce.hpp` | 領域全体 | C^1 境界条件の拘束行生成 |
