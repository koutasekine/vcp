# bfem API リファレンス

この文書は `vcp/bfem` の利用者が直接触ることの多いクラス・関数をまとめたものです。
名前空間はすべて `vcp::bfem` です。内部実装用の `detail::` は通常使いません。

## この文書でのスコープ

API は、どこの対象を扱うかで次のように読むと整理しやすくなります。

| 表記 | どこの話か | 典型的な戻り値・状態 |
|---|---|---|
| 領域非依存 | メッシュや物理要素から独立した、有理数・添字・テーブル・1 変数多項式 | `rational`, `index_map`, `poly1`, 各種 registry |
| 参照単体上 | 標準三角形・標準四面体上の Bernstein 多項式 | `bpoly<D,T>`、参照内積係数、de Casteljau 評価 |
| 物理要素上 | メッシュ中の 1 要素に写した後の局所計算 | `element_geometry`, 局所行列 `vcp::matrix` |
| 領域全体 | 全要素と大域自由度を持つ空間・場・大域行列 | `fe_space`, `rt_space`, `spmatrix`, 大域係数ベクトル |

迷った場合は、空間クラスの `stiffness` や `mass` は「領域全体の大域行列」です。
一方、`bpoly` の `inner` は「参照単体上の内積係数」であり、物理要素上の積分値にするには
要素測度が必要です。

## 共通型

### `mesh<D,T>`

スコープ: 領域全体の入力データです。ただし、この型自身は有限要素空間や自由度番号をまだ持ちません。

include:

```cpp
#include <vcp/bfem/mesh.hpp>
```

| API | 役割 |
|---|---|
| `static mesh from_lists(vertices, elements)` | 頂点座標列と要素接続から作成。要素頂点番号の範囲を検査 |
| `int num_vertices() const` | 頂点数 |
| `int num_elements() const` | 要素数 |
| `const std::array<T,D>& vertex(int v) const` | 頂点座標 |
| `const std::array<int,D+1>& element(int e) const` | 要素の頂点番号 |

要素の向きは正規化されません。面積・体積が 0、または区間型で符号が確定しない要素は、
空間や幾何を構築した時点で `degenerate_element` になります。

### `element_geometry<D,T>`

スコープ: 物理要素上です。1 つの三角形または四面体について、参照単体から物理要素への写像を保持します。

include:

```cpp
#include <vcp/bfem/geometry.hpp>
```

| API | 役割 |
|---|---|
| `static element_geometry from_vertices(verts)` | 単体の幾何情報を作成 |
| `measure() const` | 面積または体積 |
| `inv_absdet() const` | `1 / abs(det B)` |
| `orientation() const` | det の符号 |
| `grad_lambda(i,d) const` | 重心座標 `lambda_i` の物理座標方向微分 |
| `edge_matrix(r,c) const` | アフィン写像の辺行列 `B` |
| `vertices() const` | 要素頂点座標 |

例外:

| 例外 | 条件 |
|---|---|
| `degenerate_element` | det が 0、または区間型で符号が確定しない |

## Bernstein 多項式

### `bpoly<D,T>`

スコープ: 参照単体上です。`bpoly` だけではメッシュ要素の面積・体積、物理座標、境界条件は関係しません。

include:

```cpp
#include <vcp/bfem/bpoly.hpp>
```

参照単体上の Bernstein 多項式です。係数順序は `index_map<D>` の正準順序です。

| API | 役割 |
|---|---|
| `bpoly()` | 0 次の零多項式 |
| `static constant(s)` | 定数多項式 |
| `static zero(n)` | `P^n` の零多項式 |
| `static from_coeffs(n, c)` | 係数列から作成 |
| `degree() const` | 次数 |
| `size() const` | 係数数 |
| `coeff(r)` | 係数参照 |
| `coeffs() const` | 係数列 |
| `operator+=`, `operator-=` | 同次数の in-place 加減算 |
| `operator*=`, `add_scalar(s)` | スカラー倍、定数加算 |

自由関数:

| API | 役割 |
|---|---|
| `add`, `sub`, `operator+`, `operator-` | 異次数なら低次側を自動で次数上げ |
| `mul`, `operator*` | Bernstein 多項式の積 |
| `scale` | スカラー倍 |
| `elevate(u,m)` | 次数上げ。`m >= u.degree()` |
| `dlambda(u,i)` | 重心座標方向の微分 |
| `inner(u,v)` | 参照単体上の L2 内積の係数。物理積分は要素測度を掛ける |
| `eval(u, lam)` | de Casteljau による点評価 |
| `add_into`, `sub_into`, `mul_into`, `scale_into`, `elevate_into`, `dlambda_into` | 出力先を再利用する版 |

### `poly1<T>`

スコープ: 領域非依存です。ただし `compose(f,u)` は、参照単体上の `bpoly` に対する演算です。

include:

```cpp
#include <vcp/bfem/poly1.hpp>
```

1 変数多項式 `f(x) = sum a_k x^k` です。

| API | 役割 |
|---|---|
| `static from_coeffs(a)` | `T` 係数から作成 |
| `static from_rational(num_den)` | 有理数係数から `T` へ一度だけ変換して作成 |
| `degree() const` | 次数 |
| `coeff(k) const` | 係数 |
| `derivative() const` | 形式微分 |
| `compose(f,u)` | `bpoly` へ `f(u)` を合成 |
| `compose_into(dst,f,u,workspace)` | 作業領域を再利用する合成 |

### `range`, `range_refined`, `restrict_to`

include:

```cpp
#include <vcp/bfem/refine.hpp>
#include <vcp/bfem/d3/bey_table.hpp>   // 3D range_refined を使う場合
```

| API | スコープ | 役割 |
|---|---|---|
| `range(u)` | 参照単体上 | Bernstein 係数の凸包性から値域包含を返す |
| `restrict_to(u,V)` | 参照単体上 | 親単体の重心座標で与えた部分単体へ厳密制限 |
| `range_refined(u, depth)` | 参照単体上 | 2D red refinement、3D Bey refinement で範囲を精密化 |

非算術型で `range` を使う場合は `bound_traits<T>` の特殊化が必要です。

## 連続 `P^k` 空間

### `fe_space<D,T,P,SP>`

スコープ: 領域全体です。全要素を回り、大域自由度番号に散布した行列・ベクトルを返します。
局所計算は内部で `element_op` が行います。

include:

```cpp
#include <vcp/bfem/fe_space.hpp>
```

`D == 2` または `D == 3` の連続 Bernstein `P^k` 空間です。

```cpp
fe_space<2, double> fs(msh, 2);
```

| API | 返り値 | スコープ | 役割 |
|---|---|---|---|
| `fe_space(msh,n)` |  | 領域全体 | 基準次数 `n >= 1` で構築。全要素幾何も作る |
| `base_degree() const` | `int` | 領域全体 | 構築時次数 |
| `num_elements() const` | `int` | 領域全体 | 要素数 |
| `dofs(m)` | `const dofmap<D>&` | 領域全体 | 次数 `m` の自由度写像 |
| `ndof(m)` | `int` | 領域全体 | 次数 `m` の大域自由度数 |
| `zero_function(m)` | `fe_function<D,T,P>` | 領域全体 | 零関数 |
| `function_from_coeffs(m,c)` | `fe_function<D,T,P>` | 領域全体 | 係数ベクトルから関数を作成 |
| `stiffness(m)` | `spmatrix_t` | 領域全体 | 全要素を組み立てた `(grad psi_j, grad psi_i)` |
| `mixed_mass(a,b)` | `spmatrix_t` | 領域全体 | 全要素を組み立てた `(phi_j^b, phi_i^a)`。長方形可 |
| `load(f,uh,m)` | `matrix<T,P>` | 領域全体 | 全要素を組み立てた `(f(uh), psi_i)` |
| `weighted_mass(fprime,uh,m)` | `spmatrix_t` | 領域全体 | 全要素を組み立てた `(f'(uh) psi_j, psi_i)` |
| `convection(b,m)` | `spmatrix_t` | 領域全体 | 全要素を組み立てた `(b . grad psi_j, psi_i)`。`b` は関数配列または定数配列 |
| `inner(u,v)` | `T` | 領域全体 | 全領域の `(u,v)` |
| `scalar_ff(f,uh)` | `T` | 領域全体 | 全領域の `(f(uh), f(uh))` |
| `elevate(u,m)` | `fe_function` | 領域全体 | 大域関数を変えずに次数 `m` へ上げる |
| `dual(u,m)` | `dual_vector` | 領域全体 | 大域 moment `((u,psi_i))` を直接作る |
| `eval(u,e,lam)` | `T` | 領域全体から要素上評価 | 大域関数を要素 `e` に制限し、重心座標 `lam` で評価 |

### `fe_function<D,T,P>` と `dual_vector<D,T,P>`

スコープ: 領域全体です。係数列は大域自由度順で並びます。

include:

```cpp
#include <vcp/bfem/fe_function.hpp>
```

| 型 | API | 役割 |
|---|---|---|
| `fe_function` | `degree()` | 次数 |
| `fe_function` | `coeffs()` | 大域 Bernstein 係数。直接編集可 |
| `dual_vector` | `test_degree()` | test 側次数 |
| `dual_vector` | `values()` | moment ベクトル |

`fe_function` は通常 `fe_space` の factory から作ります。`dual_vector` から
`fe_function` へ戻す API はありません。戻すには質量行列の逆が必要で、それは暗黙には行いません。

### `dofmap<D>`

スコープ: 領域全体です。局所要素番号と大域自由度番号を結びます。

`fe_space::dofs(m)` から取得します。

| API | 役割 |
|---|---|
| `ndof()` | 大域自由度数 |
| `degree()` | 次数 |
| `local_size()` | 1 要素あたり自由度数 |
| `num_elements()` | 要素数 |
| `global_dof(e,r)` | 要素 `e` の局所番号 `r` から大域自由度へ |
| `dof_sign(e,r)` | `P^k` では常に `+1` |
| `boundary_dofs()` | 全境界自由度 |
| `boundary_dofs(edges_or_faces)` | 2D では辺、3D では面を指定した境界自由度 |
| `boundary_edge_ids()` | 2D の境界辺 ID |
| `boundary_face_ids()` | 3D の境界面 ID |

## Dirichlet 縮約

### `dirichlet_reduction<T,P,SP>`

スコープ: 領域全体です。大域行列・大域ベクトルに対して自由度を削除します。

include:

```cpp
#include <vcp/bfem/dirichlet.hpp>
```

斉次 Dirichlet 条件 `u = 0` のため、指定した自由度を削除します。

| API | 役割 |
|---|---|
| `dirichlet_reduction(full_size, bdofs)` | 制約自由度リストから構築。重複・範囲外は例外 |
| `full_size()` | 元の自由度数 |
| `reduced_size()` | 縮約後自由度数 |
| `to_reduced(full)` | 元番号から縮約番号。制約自由度は `-1` |
| `to_full(reduced)` | 縮約番号から元番号 |
| `reduce(A)` | 正方行列の制約行・列を削除 |
| `reduce(v)` | ベクトルの制約成分を削除 |
| `expand(vr)` | 縮約ベクトルを全自由度へ戻す。制約成分は 0 |
| `reduce_rows(A)` | 行だけ縮約 |
| `reduce_cols(A)` | 列だけ縮約 |

## `RT^k` と broken `P_l`

### `rt_space<D,T,P,SP>`

スコープ: 領域全体です。`mass()` は大域 RT 質量行列を返します。
一方、`interpolate(provider)` の provider が返す `bpoly` は各要素上で使う物理成分です。

include:

```cpp
#include <vcp/bfem/rt/rt_space.hpp>
```

`D == 2` または `D == 3` の Raviart-Thomas 空間です。

| API | 返り値 | 役割 |
|---|---|---|
| `rt_space(msh,k)` |  | `k >= 0` の RT 空間 |
| `order() const` | `int` | RT 次数 |
| `ndof() const` | `int` | 大域自由度数 |
| `num_elements() const` | `int` | 要素数 |
| `num_vertices() const` | `int` | 頂点数 |
| `dofs() const` | `dofmap_type` | RT 自由度写像 |
| `geometry(e) const` | `element_geometry` | 要素幾何 |
| `zero_field()` | `rt_field` | 零 RT 場 |
| `field_from_coeffs(c)` | `rt_field` | 係数から RT 場 |
| `mass()` | `spmatrix_t` | `(sigma_j, tau_i)` |
| `interpolate(provider)` | `rt_field` | 物理成分で与えた場を RT DOF へ補間 |

`interpolate` の provider は、2D では `std::pair<bpoly, bpoly>`、3D では
`std::array<bpoly,3>` を要素ごとに返します。入力場は法線連続であることが前提です。

### `broken_space<D,T,P,SP>`

スコープ: 領域全体です。ただし自由度は要素ごとに独立で、係数ベクトルは要素ブロック順です。

include:

```cpp
#include <vcp/bfem/rt/broken_space.hpp>
```

要素ごとに不連続な `P_l` 空間です。

| API | 返り値 | 役割 |
|---|---|---|
| `broken_space(msh,l)` |  | `l >= 0` の broken 空間 |
| `order() const` | `int` | 次数 |
| `ndof() const` | `int` | 大域自由度数 |
| `local_size() const` | `int` | 要素ごとの係数数 |
| `dofs() const` | `broken_dofmap` | ブロック型自由度写像 |
| `geometry(e) const` | `element_geometry` | 要素幾何 |
| `zero_field()` | `broken_field` | 零場 |
| `field_from_coeffs(c)` | `broken_field` | 係数から場 |
| `mass()` | `spmatrix_t` | block diagonal mass |
| `inner(u,v) const` | `T` | broken 場の L2 内積 |

### `rt_assemble.hpp` の自由関数

スコープ: 領域全体です。複数の空間をまたいで大域行列・大域ベクトル・全領域スカラーを作ります。

include:

```cpp
#include <vcp/bfem/rt/rt_assemble.hpp>
```

| API | 役割 |
|---|---|
| `assemble_div_mass(bs, rs)` | 行 = broken、列 = RT、`(div sigma_j, q_i)` |
| `assemble_cross_grad(rs, fs, m)` | 行 = RT、列 = `P^m`、`(sigma_i, grad psi_j)` |
| `flux_error_sq(rs, sig, fs, u)` | `||grad u - sig||^2` |
| `div_residual_sq(rs, sig, f, fs, u)` | `||div sig + f(u)||^2` |
| `broken_load(bs, f, fs, uh)` | broken test での `(f(uh),q_i)` |
| `projection_error_sq(bs, f, fs, uh)` | `||f(uh) - Pi_l f(uh)||^2` |
| `assemble_mixed_mass(fs, m, bs)` | 行 = `P^m`、列 = broken、`(phi_i,q_j)` |

## Scott-Vogelius 部品

### `vfe_space<D,T,P,SP>`

スコープ: 領域全体です。`fe_space` の大域自由度を成分ごとに並べたベクトル値空間です。

include:

```cpp
#include <vcp/bfem/sv/vfe_space.hpp>
```

`fe_space` を成分ごとに並べたベクトル値空間です。`vfe_space` は `fe_space` を所有しないため、
参照先の `fe_space` は `vfe_space` より長生きさせます。

| API | 役割 |
|---|---|
| `vfe_space(msh, scalar_fe_space)` | ベクトル値空間を構築 |
| `scalar() const` | 参照先の `fe_space` |
| `ndof(m) const` | `D * scalar.ndof(m)` |
| `global_dof(d,i,m) const` | 成分 `d`、スカラー自由度 `i` の大域番号 |
| `zero_function(m) const` | 零ベクトル関数 |
| `function_from_coeffs(m,c) const` | 係数からベクトル関数 |
| `boundary_dofs(m) const` | 全成分の境界自由度 |
| `component(u,d) const` | 成分を `fe_function` としてコピー抽出 |
| `set_component(u,d,f) const` | 成分を書き戻す |

### `sv_assemble.hpp` の自由関数

スコープ: 領域全体です。局所的には要素ごとに多項式演算しますが、返り値は大域行列・大域ベクトル・全領域スカラーです。

include:

```cpp
#include <vcp/bfem/sv/sv_assemble.hpp>
```

| API | 役割 |
|---|---|
| `assemble_vector_stiffness(vs,m)` | ベクトル値剛性 `(grad u : grad v)` |
| `assemble_div_velocity(bs,vs,m)` | 行 = broken 圧力、列 = 速度、`(div u,q)` |
| `assemble_advection(vs,w,m)` | 固定速度 `w` に対する移流行列 `(w.grad u,v)` |
| `assemble_advection_derivative(vs,w,m)` | Newton 微分側 `((u.grad)w,v)` |
| `advection_vector(vs,a,b,m)` | `((a.grad)b, phi_i)` |
| `advection_scalar(vs,a,b,g)` | `((a.grad)b,g)`。`g` は `vfe_function` または `grad fe_function` |
| `div_field(bs,vs,u)` | `div u` を broken 場として返す |
| `div_norm_sq(vs,u)` | `||div u||^2` |

### 特異拘束と線形縮約

スコープ: 領域全体です。特異性検出はメッシュ幾何を見ますが、生成される拘束は大域自由度に対する行です。

include:

```cpp
#include <vcp/bfem/sv/sv_rows.hpp>
#include <vcp/bfem/sv/linear_reduction.hpp>
#include <vcp/bfem/sv/alfeld.hpp>
```

| API | 役割 |
|---|---|
| `alfeld_refine(msh)` | 各単体を重心分割 |
| `sv_pressure_constraints<D,T>(msh,n)` | Scott-Vogelius 圧力拘束行を生成 |
| `rows()` | 独立化済み整数拘束行 |
| `num_singular_vertices()` | 2D の特異頂点数 |
| `num_singular_edges()` | 3D の特異辺数 |
| `hint(provenance)` | 既知の十分条件に基づく観測情報 |
| `linear_reduction<T,P,SP>(full_size, rows)` | 一般線形拘束を消去する縮約 |
| `linear_reduction::reduce`, `reduce_rows`, `reduce_cols`, `expand` | 行列・ベクトルの縮約と復元 |

## 2D C^1 要素

### `c1_space<2,T,P,SP>`

スコープ: 領域全体です。内部では C^1 用の物理要素上局所行列を作り、大域 C^1 自由度へ散布します。
`eval_grad` や `eval_hess` は、大域関数を指定要素へ制限して評価します。

include:

```cpp
#include <vcp/bfem/c1/c1_space.hpp>
```

2D の C^1 Argyris 型空間です。次数は `k >= 5` です。

| API | 返り値 | 役割 |
|---|---|---|
| `c1_space(msh,k)` |  | 基準次数 `k >= 5` で構築 |
| `base_degree() const` | `int` | 基準次数 |
| `num_elements()`, `num_vertices()`, `num_edges()` | `int` | メッシュ情報 |
| `mesh_ref() const` | `mesh<2,T>&` | 元メッシュ |
| `geometry(e) const` | `element_geometry<2,T>` | 要素幾何 |
| `edge_inv_tsq(ed) const` | `T` | 辺長二乗の逆数 |
| `dofs(m)` | `c1_dofmap&` | C^1 自由度写像。`m >= base_degree()` |
| `ndof(m)` | `int` | 大域自由度数 |
| `zero_function(m)` | `c1_function` | 零関数 |
| `function_from_coeffs(m,c)` | `c1_function` | 係数から関数 |
| `stiffness(m)` | `spmatrix_t` | `(grad u, grad v)` |
| `mixed_mass(a,b)` | `spmatrix_t` | `(u,v)` |
| `weighted_mass(fprime,uh,m)` | `spmatrix_t` | 重み付き質量 |
| `load(f,uh,m)` | `matrix<T,P>` | 荷重 |
| `inner(u,v)` | `T` | L2 内積 |
| `scalar_ff(f,uh)` | `T` | `(f(uh),f(uh))` |
| `eval(u,e,lam)` | `T` | 値評価 |
| `eval_grad(u,e,lam)` | `std::array<T,2>` | 物理勾配 |
| `eval_hess(u,e,lam)` | `std::array<T,3>` | `xx, xy, yy` の Hessian |
| `laplacian_matrix(m)` | `spmatrix_t` | `(Delta u, Delta v)` |
| `hessian_matrix(m)` | `spmatrix_t` | `(D2 u : D2 v)` |
| `laplacian_mixed(m,l)` | `spmatrix_t` | 行 = broken `P_l`、列 = C^1 |
| `laplacian_load(f,uh,m)` | `matrix<T,P>` | `(f(uh), Delta phi_i)` |
| `laplacian_residual_sq(f,uh)` | `T` | `||Delta uh + f(uh)||^2` |
| `laplacian_field(u,bs)` | `broken_field` | `Delta u` を broken `P_{m-2}` として返す |
| `interpolate(provider,m)` | `c1_function` | 要素ごとの pullback 多項式から補間 |
| `elevate(u,m)` | `c1_function` | C^1 関数を次数上げ |

### `c1_dofmap`

スコープ: 領域全体です。頂点・辺・要素内部の C^1 自由度を大域番号で扱います。

`c1_space::dofs(m)` から取得します。

| API | 役割 |
|---|---|
| `degree()`, `ndof()`, `local_size()` | 次数・自由度数 |
| `global_dof(e,r)`, `dof_sign(e,r)` | 局所から大域への写像 |
| `vertex_dof(v,c)` | 頂点自由度。`c = 値, dx, dy, dxx, dxy, dyy` |
| `edge_trace_dof(ed,i)` | 辺上 trace 点値自由度 |
| `edge_nd_dof(ed,j)` | 辺上法線微分自由度 |
| `interior_dof(e,r)` | 要素内部自由度 |
| `boundary_edge_ids()` | 境界辺 ID |
| `edge_verts(ed)` | 辺の両端頂点 |

### C^1 境界拘束

include:

```cpp
#include <vcp/bfem/c1/c1_reduce.hpp>
```

| API | 役割 |
|---|---|
| `c1_boundary(sp,m,c1_bc_dirichlet)` | `u = 0` の本質境界拘束 |
| `c1_boundary(sp,m,c1_bc_simply_supported)` | simply-supported の本質拘束。現実装では Dirichlet と同じ |
| `c1_boundary(sp,m,c1_bc_clamped)` | `u = 0` と `du/dnu = 0` |
| `c1_boundary(sp,m,kind,edges)` | 指定境界辺だけに拘束を作る |

戻り値は rational 係数の `sv_constraint_q` 列です。`linear_reduction` に渡して縮約できます。
`c1_boundary` は `c1_coord_traits<T>` を使って、境界方向の符号判定と拘束係数用の有理数化を行います。
`double` や独自スカラー型で使う場合は、この traits を利用者側で特殊化してください。

## 代表的な例外

| 例外 | 主な条件 |
|---|---|
| `std::invalid_argument` | 次数が範囲外、係数サイズ不一致、自由度範囲外、別メッシュ混用 |
| `degenerate_element` | 要素の測度が 0、または区間型で符号不定 |
| `std::logic_error` | `element_op` で `set_geometry` 前に局所演算を呼んだ場合 |
| `sv_indeterminate_singularity` | 区間座標で SV 特異性判定が確定しない |
| `c1_indeterminate_corner` | 区間座標で C^1 境界角分類が確定しない |
