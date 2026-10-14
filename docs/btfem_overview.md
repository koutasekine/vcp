# btfem 概要

`vcp/btfem` は、放物型方程式の全離散 Petrov–Galerkin 法
(Hashimoto–Kinoshita–Nakao, arXiv:1811.04261v2)のためのヘッダ群です。
空間方向は `vcp/bfem` の連続 Bernstein `P^n` 空間(`fe_space`)をそのまま使い、
時間方向を 1 次元 Bernstein 多項式で張った時間空間を接続します。

bfem と同じく、要素積分に数値求積は使いません。時間方向の積分テーブルは bfem L0 の
1 次元(`D = 1`)有理数テーブルから有理数のまま合成し、スカラー型 `T` へ一度だけ変換します。
`T` に区間型を使えば、行列・ベクトル・スカラー量は厳密値の包含として計算できます。
有理数型なら厳密等式が成り立ちます。

## できること

| 分野 | 主なヘッダ | 主な機能 |
|---|---|---|
| 時間メッシュ | `<vcp/btfem/time_mesh.hpp>` | `[t0, t1]` の `m` 等分、小区間幅 `k` と `1/k`、退化検出 |
| 時間空間 | `<vcp/btfem/time_space.hpp>` | 連続区分 Bernstein 空間 `Vbar_k^q`(左端零条件なし)と `V_k^q`(左端係数除去)、大域時間行列、局所テーブル |
| 空間時間関数 | `<vcp/btfem/st_poly.hpp>` | 係数行列(時間ブロック主)の値型 `st_poly<D,T,P>` |
| Kronecker 因子 | `<vcp/btfem/kron.hpp>` | `kron_term`(時間 `matrix` × 空間 `spmatrix`)、展開 `kron_expand`、作用 `kron_apply`、ブロック Dirichlet 縮約、時間初期条件の縮約 |
| 全離散空間 | `<vcp/btfem/pg_space.hpp>` | 行列群の因子 `A M B U W Y`、時間依存係数項、非線形スキーム素材、断面・点評価・値域、初期誤差 |

## 計算の考え方

全離散試行空間は `V_k^q ⊗ S_h^n`、基底は `psi_a(t) phi_i(x)`、試験基底は主論文どおり
`∂_t psi_a ⊗ phi_i` だけです。自由度順序は時間ブロック主 `(a, i) ↦ a·n_x + i` に固定されており、
各時間ブロックはそのまま bfem `fe_function` の係数ベクトルになります。

大域行列は「時間行列と空間行列の対」(`kron_term`)を正本とし、`mn × mn` の `spmatrix` は
`kron_expand` が作ります。時間依存係数項 `(c phi_j, ∂_t phi_i)` も対の列で返り、
その和を `kron_expand` に渡します。方程式の係数(拡散係数 ν など)は API に埋め込まず、
利用者が `kron_term::scale` に入れます。

btfem の説明では、話題のスコープを次の 5 種類に分けると読みやすくなります。

| スコープ | 意味 | 代表例 |
|---|---|---|
| 領域非依存 | メッシュを見ない、純粋な添字・有理数・テーブルの話 | `detail::time_tables`、`time_bc` |
| 参照区間上 | 参照区間 `[0,1]` 上の時間 Bernstein 多項式。小区間幅 `k` はまだ掛からない | 参照テーブル `Mref` / `Sref` / `Bref` |
| 参照単体上 | bfem と同じ、標準単体上の空間 Bernstein 多項式 | `bpoly<D,T>`(bfem) |
| 物理要素上 | 1 小区間 × 1 空間要素に写した後の局所計算。`k`、`1/k`、要素測度が入る | `time_space` の `local_*`、`detail::st_local` |
| 空間時間全体 | 全小区間・全要素を回って大域自由度へ散布した後の行列・ベクトル・場 | `time_space` の `mass` 等、`pg_space`、`st_poly`、`kron_expand` |

1. `time_mesh<T>` が `[t0, t1]` を `m` 等分し、`k` と `1/k` を保持します(除算はここで 1 回だけ)。
2. `time_space<T>` が時間自由度番号 `g(l, r)` と大域時間行列、物理化済み局所テーブルを提供します。
3. `pg_space<D,T,P,SP>` が bfem の `fe_space` を参照し、小区間 × 要素のループで局所テンソル多項式
   (`detail::st_local`)を組んで行列群の因子・ベクトル・スカラー量を作ります。
4. `kron_expand` / `kron_apply` が因子から `spmatrix` または行列ベクトル積を作ります。

加算順序は「小区間 → 要素 → 時間 rank → 局所添字」で固定され、COO の重複結合は bfem の決定的結合を
使います。同一入力・単スレッドならビット単位で同じ結果になります。

## 責務の境界

btfem は行列・ベクトル・内積値・ノルム二乗の生成までを担います。線形方程式の求解、
Newton 反復(`spnewton`)、固有値・作用素ノルム、流束 σ の決定、`u_0` から左端係数を決める
補間・射影、精度保証定数の計算は利用者側または別トラックの責務です。
時刻断面は `fe_function` で返るので、可視化は bfem `graphics.hpp` へ渡します。

RT ⊗ 時間、broken ⊗ 時間、流束系は BT-2 の範囲です。

## OpenMP 並列化

`pg_space` の時間依存係数項(`weighted_dt_mass_terms`、空間時間 `b` の
`convection_dt_terms`)の要素ループは OpenMP に対応しています。bfem と同じ
`-DVCP_BFEM_NOMP` / `-DVCP_NOMP` で無効化され、スイッチは bfem と共有です。
スレッドごとの作業領域は `pg_space` が保持し、同じ形の 2 回目以降の呼び出しでは
要素ループ内でヒープ確保が起きません。スレッド数差のビット一致は契約ではありませんが、
1 スレッドの OpenMP ビルドと NOMP ビルドは同じ順序で加算するので一致します。
小区間ループ、ベクトル・スカラー経路は逐次です。

## スカラー型 `T`

bfem と同じ契約です。`bfem::rational`(テスト用 `ratx`)、`double`、`kv::dd`、`kv::mpfr<N>`、
`kv::interval<double>`、`kv::interval<kv::dd>`、`kv::interval<kv::mpfr<N>>` で
全ヘッダがコンパイルされます。区間型・kv 点型では bfem の opt-in ヘッダも include します。

```cpp
#include <vcp/bfem/bound_traits_kv.hpp>
#include <vcp/bfem/geometry_traits_kv.hpp>
#include <vcp/bfem/convert_traits_kv.hpp>
```

`range` / `range_composed` を kv の点型(`kv::dd`、`kv::mpfr<N>`)や有理数型で使う場合は、
bfem `refine.hpp` の契約どおり `bound_traits<T>` を `point_bound_traits<T>` から派生させる
特殊化を利用者側で与えます。

除算の回数はプログラム全体で「空間要素数(bfem 幾何)+ 時間メッシュ数(`time_mesh::uniform`)」
回です。有理数から `T` への enclose-once 変換(`1/m`、時間テーブル)は bfem L0 と同じ
変換予算に属し、プロセス共有キャッシュで 1 回限りです。

## 疎行列ポリシー

空間因子と展開結果は `vcp::spmatrix<T, SP>` で返ります。`SP` の既定は `vcp::spmats<T>` で、
区間型では bfem と同様に区間対応の sparse policy を `pg_space<D,T,P,SP>` へ渡します。
btfem は `spmatrix` を戻り値とデータ保持にだけ使い、`mul_vec` 等の演算 API には依存しません
(`kron_apply` は成分列挙で計算します)。返る `spmatrix` は triplet 挿入済み・`finalize` 未実施です。

## ファイル構成と役割

| ファイル | スコープ | 役割 |
|---|---|---|
| `time_mesh.hpp` | 空間時間全体の入力 | 等間隔時間メッシュ、`degenerate_interval`、唯一の除算 |
| `detail/time_tables.hpp` | 領域非依存 → 参照区間上 | L0 の `D = 1` テーブルからの有理数合成と `T` 段(プロセス共有) |
| `time_space.hpp` | 参照区間上 → 空間時間全体 | 時間自由度番号、大域時間行列、物理化済み局所テーブル |
| `st_poly.hpp` | 物理要素上・空間時間全体 | `st_poly` 値型、`detail::st_local` の演算、gather |
| `kron.hpp` | 空間時間全体 | `kron_term`、展開・作用、ブロック縮約、`initial_reduction` |
| `pg_space.hpp` | 空間時間全体 | 行列群の因子、時間依存係数項、非線形素材、断面、初期誤差、値域 |

依存方向は `time_mesh ← time_space ← st_poly ← pg_space`、`kron.hpp` は `pg_space.hpp` と
独立に利用できます。`vcp/bfem/` は無改造で、`detail::coo_buffer` / `detail::spm_adapter` /
`detail::scatter_matrix` / `detail::bpoly_access` を呼びます。

テストは `sandbox/tests/btfem_*.cpp`、Makefile ターゲットは `sandbox/Makefile` の
`btfem_*` / `run_btfem_*`、全スイート実行は `sandbox/scripts/run_btfem_all_suites.sh` です。
API の詳細は `docs/btfem_api_reference.md` を参照してください。

## BT-2: RT ⊗ 時間、broken ⊗ 時間、流束系

BT-2 は、放物型の残差を `L²(J; H⁻¹)` で扱うための**流束系の素材**を追加します。
bfem の `rt_space` / `broken_space` / `rt_assemble.hpp` を無改造で使い、
時間方向を BT-1 と同じ Bernstein 空間 `Vbar_k^q` / `V_k^q` で張ります。
どの不等式で組み合わせるかは理論側(オーナー)の責務で、本層は量だけを出します。

| 分野 | 主なヘッダ | 主な機能 |
|---|---|---|
| RT ⊗ 時間 | `<vcp/btfem/st_rt_space.hpp>` | `st_rt_field`(係数行列、時間ブロック主、列は bfem `rt_dofmap` の符号付き大域順)、`st_rt_space`(ファクトリ、断面、質量因子 `U_t ⊗ M_RT`) |
| broken ⊗ 時間 | `<vcp/btfem/st_broken_space.hpp>` | `st_broken_field`、`st_broken_space`(ファクトリ、断面、質量因子 `U_t ⊗ M_broken`、空間時間内積 `inner`) |
| 流束系 | `<vcp/btfem/st_rt_assemble.hpp>` | bfem `rt_assemble.hpp` 7 関数の時間版: `st_div_field`、因子 `st_assemble_div_mass` / `st_assemble_cross_grad` / `st_assemble_mixed_mass`、ノルム二乗 `st_flux_error_sq` / `st_div_residual_sq` / `st_projection_error_sq`、ベクトル `st_broken_load` |

### 構成規則(追加の余地)

- 時間に依存しない双一次形式 `b_x(·,·)` の時間版は `kron_term{ T(1), Vt_row.mixed_mass<P>(Vt_col), B_x }` の 1 行で追加できます(`B_x` は bfem の既存アセンブラ)。S4〜S6 がこの形です。
- 時間に依存する係数つきの形式は BT-1 `pg_space::weighted_dt_mass_terms` と同じ「時間係数分解 → `(l, σ)` ごとの `kron_term` 列」で追加します(初版未実装)。
- ノルム二乗は「小区間 × 要素の局所テンソル多項式を作る → 時間局所質量 × 空間局所内積の二重和」で追加します。共通部品は `detail::st_quadratic_sq` で、S7〜S10 がこの形です。
- 空間側に新しい局所演算が要る場合は bfem 側のトラックになります。

### 残差多項式 `w = g(ū) − ∂_t ū − b·∇ū`

S8〜S10 の `w` は小区間 × 要素ごとに局所テンソル(`detail::st_local`)として組み、大域化しません。
`g(ū)` は BT-1 の合成、`∂_t ū` は時間微分、`b·∇ū` は定数 `b`(`std::array<T, D>`)または
空間時間 `b`(`std::array<st_poly, D>`、全成分が同じ `(q_b, bc)`)で、`b` なしの overload は
bfem の `f(u)` 形に対応します。3 項は時間次数を `max` に揃えてから(L0 の `D = 1` 次数上げ
テーブル)、空間次数は bfem `sub_into` が揃えて和を取ります。

### 符号と加算順

RT 自由度は辺・面の向きで符号を持ちます(bfem Y2)。時間方向はテンソル積なので、
符号は空間自由度ごとに決まり全時間 rank に同じ符号が掛かります(`detail::st_gather_rt`)。
`to_vector` / `field_from_vector` は符号を変えません。加算順は二次形式が
「小区間 → 要素 → `(r, r')` 行優先」、ベクトルが「`l → e → s → a → i`」で固定され、
同一入力・単スレッドでビット再現します。BT-2 の自由関数には並列区間がなく、
OpenMP は因子 S3〜S6 の内側にある bfem アセンブラのものだけです。

### 係数 ν の扱い

S7 は bfem `flux_error_sq` と同じ `‖∇ū − σ‖²` を返し、ν は API に持ちません。
`σ ≈ ν∇ū` の流束を使う場合は `σ/ν` を渡すか結果に `ν²` を掛けます。因子の係数は
`kron_term::scale` に入れます。

### ファイル構成(BT-2 追加分)

| ファイル | スコープ | 役割 |
|---|---|---|
| `st_rt_space.hpp` | 空間時間全体 | `st_rt_field` / `st_rt_space`、`detail::st_rt_local` / `st_gather_rt` |
| `st_broken_space.hpp` | 空間時間全体 | `st_broken_field` / `st_broken_space`、`detail::st_broken_dofmap_view` |
| `st_rt_assemble.hpp` | 物理要素上・空間時間全体 | `st_div_field`、因子 S4〜S6、残差 `detail::st_residual_into`、二次形式 `detail::st_quadratic_sq`、S7〜S10、`detail::st_sub_into` / `st_elevate_time_into` |

依存方向は `st_rt_space` / `st_broken_space` ← `st_rt_assemble`(`pg_space.hpp` にも依存)です。
BT-1 の 6 ヘッダと `vcp/bfem/` は無改造です。テストは `sandbox/tests/btfem_rt_*.cpp`、
Makefile ターゲットは `btfem_rt_*` / `run_btfem_rt_*`(`run_btfem_rt_all`)で、
`run_btfem_all_suites.sh` にも含まれます。

## BT-3: 定数層(全離散 Petrov–Galerkin 法の誤差定数)

`vcp/btfem/constants/` は、線形熱方程式 ∂_t u − νΔu = f(u(a) = 0)の全離散 Petrov–Galerkin
近似(主論文 (2.3) / 定数論文 (10)。試行 V_h ⊗ {p ∈ V_τ^k : p(a) = 0}、試験 V_h ⊗ DV_τ^{k−1})の
3 つの誤差定数を、2 つの射影誤差定数と 1 つの離散量だけから組み立てます。
名前空間は `vcp::btfem::constants`、スカラー `T` は kv の区間型に限定されます
(`vcp::bfem::constants::detail::interval_scalar_contract`。`double` や有理数型はコンパイルエラー)。

| 分野 | 主なヘッダ | 主な機能 |
|---|---|---|
| 時間射影誤差定数 | `<vcp/btfem/constants/time_projection_constants.hpp>` | `time_projection_error_constant(d, tau)` = C_τ(d) = c_d·τ。閉じた式のみ(整数・π・平方根)、d に上限なし |
| 逆不等式 λ_h | `<vcp/btfem/constants/inverse_inequality.hpp>` | `inverse_inequality_h01_bound<D>(Th, n)`: 縮約後の K_x v = λ U_x v の最大固有値の上界(bfem κ_h² と同じ密 `eigsymge` 経路) |
| 3 定数 | `<vcp/btfem/constants/heat_pg_constants.hpp>` | 材料 (C_Ω², C_τ, λ_h, ν) から K = C̃₁²、L = c̃₀²、C̃₁、c̃₀、C̃₀ を組む純関数、空間側材料の一括取得と再利用、`time_mesh` 版・一括版 |

計算の考え方は bfem 定数層と同じです。数値は閉じた式(整数、`kv::constants<T>::pi()`、`sqrt`)か
bfem の辞書からのみ入り、小数リテラルは書きません。空間方向の C_Ω(h)² は bfem の
`ritz_projection_error_constant_h01_sq(Th, n)` をそのまま使い(2D / 3D、辞書供給。3D は
`registry_3d` 未登録クラスで bfem の例外が伝播)、新しい計算経路は作っていません。
λ_h と 3 定数は斉次 Dirichlet 空間を仮定するので `_h01` を付け、C_τ は境界条件を持たないので
印なしです(bfem VER-1 R23)。二乗量 K、L が正本で、非二乗版は平方根ラッパーです(bfem R17)。

### ファイル構成(BT-3 追加分)

| ファイル | スコープ | 役割 |
|---|---|---|
| `constants/time_projection_constants.hpp` | 参照区間上(τ でスケール) | C_τ(d) の 4 分岐(d = 0, 1: NYK 1998 最良値、d = 2: NYK 1998 定理 2、d ≥ 3: KY 1999 定理 2.3)、`detail::time_projection_coefficient_{nyk2,ky}` |
| `constants/inverse_inequality.hpp` | メッシュ全体・斉次 Dirichlet | `fe_space` → `dirichlet_reduction` → 密化 → `compsym` → `eigsymge` → 上端の最大値(点区間) |
| `constants/heat_pg_constants.hpp` | メッシュ全体・斉次 Dirichlet | `heat_pg_space_material_set_h01` / `heat_pg_space_materials_h01`、`heat_pg_error_constant_set_h01` / `heat_pg_error_constants_h01`(式のみ・材料版 2 形・一括版) |

依存方向は `time_projection_constants` / `inverse_inequality` ← `heat_pg_constants` で、
bfem(`poisson_constants.hpp`、`fe_space.hpp`、`dirichlet.hpp`)と BT-1 の `time_mesh.hpp` を
無改造で使います。`vcp/ldbase.hpp` はテストからのみ include します(G3 の突合先)。
テストは `sandbox/tests/btfem_const_*.cpp`、Makefile ターゲットは `btfem_const_*` /
`run_btfem_const_*`(`run_btfem_const_all`)で、`run_btfem_all_suites.sh` にも含まれます。
`kv::interval<kv::dd>` での出力は、kv の dd 除算の有向丸めに関する観測
(`sandbox/docs/issues/kv_dd_interval_division_rounding_issue.md`)が解消するまで
真の包含とは主張しません。
