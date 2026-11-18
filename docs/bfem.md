# bfem — Bernstein 多項式による有限要素法

`vcp/bfem` は、Bernstein 多項式を基底にした有限要素法のヘッダ群です。
三角形・四面体メッシュ上で、次の空間の行列・ベクトル・スカラー量を作ります。

- 連続 Lagrange 型の `P^k`(2D・3D)
- Raviart–Thomas 型の `RT^k` と、要素ごとに不連続な broken `P_l`(2D・3D)
- Scott–Vogelius 要素のためのベクトル値空間と圧力拘束(2D・3D)
- Crouzeix–Raviart 型の CR1(2D・3D)
- 2D の C¹ 要素(Argyris 型)

要素積分に数値積分を使わないことが特徴です。Bernstein 基底の積・次数上げ・
質量行列などを有理数のテーブルとして作り、要素型 `T` へ一度だけ変換して
使います。`T` に区間型を使えば、行列・ベクトル・スカラー量は厳密値の包含として
計算されます。

この他に、多角形領域のメッシュ生成と、有限要素法の誤差定数の精度保証付き計算
(`vcp/bfem/constants/`)を含みます。

名前空間は `vcp::bfem`(誤差定数は `vcp::bfem::constants`)です。`detail::` は
内部実装です。

bfem は有限要素空間と、行列・ベクトル・多項式量を作る部品です。線形方程式の
求解、Newton 反復、固有値計算は、[spmatrix.md](spmatrix.md)、
[newton.md](newton.md) などの機能を使って利用者側で行います。Scott–Vogelius の
部品も、Navier–Stokes 方程式のソルバではなく、速度・圧力の空間と行列を作る
部品です。

## 必要なヘッダ

| 用途 | ヘッダ |
| --- | --- |
| メッシュ | `<vcp/bfem/mesh.hpp>` |
| メッシュ生成(2D 多角形) | `<vcp/bfem/meshgen.hpp>` |
| メッシュ生成(Delaunay・局所細分・3D 押し出し) | `<vcp/bfem/meshgen_delaunay.hpp>`、`<vcp/bfem/meshgen_adapt.hpp>`、`<vcp/bfem/meshgen3.hpp>` |
| メッシュの保存・読み込み・型変換 | `<vcp/bfem/meshgen_io.hpp>`、`<vcp/bfem/meshgen_convert.hpp>` |
| 連続 `P^k` | `<vcp/bfem/fe_space.hpp>` |
| 斉次 Dirichlet 縮約 | `<vcp/bfem/dirichlet.hpp>` |
| 非線形項 `f(u)` | `<vcp/bfem/poly1.hpp>` |
| `RT^k` と broken `P_l` | `<vcp/bfem/rt/rt_space.hpp>`、`<vcp/bfem/rt/broken_space.hpp>`、`<vcp/bfem/rt/rt_assemble.hpp>` |
| Scott–Vogelius | `<vcp/bfem/sv/vfe_space.hpp>`、`<vcp/bfem/sv/sv_assemble.hpp>`、`<vcp/bfem/sv/sv_rows.hpp>`、`<vcp/bfem/sv/linear_reduction.hpp>`、`<vcp/bfem/sv/alfeld.hpp>` |
| CR1 | `<vcp/bfem/cr1/cr1_space.hpp>`、`<vcp/bfem/cr1/vcr1_space.hpp>`、`<vcp/bfem/cr1/vcr1_divfree.hpp>` |
| 2D C¹ | `<vcp/bfem/c1/c1_space.hpp>`、`<vcp/bfem/c1/c1_reduce.hpp>` |
| 誤差定数 | `<vcp/bfem/constants/poisson_constants.hpp>`、`<vcp/bfem/constants/sobolev_constants.hpp>`、`<vcp/bfem/constants/stokes_constants.hpp>` |
| 図示用の出力 | `<vcp/bfem/graphics.hpp>` |

区間型で使う場合は、kv のヘッダの後に次も include します。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/bfem/bound_traits_kv.hpp>
#include <vcp/bfem/geometry_traits_kv.hpp>
#include <vcp/bfem/convert_traits_kv.hpp>
```

## 計算の考え方

bfem の API は、どの範囲の量を扱うかで次の 4 つに分けると整理しやすくなります。

| 範囲 | 意味 | 代表例 |
| --- | --- | --- |
| 領域によらない | メッシュと無関係な有理数・添字・テーブル・1 変数多項式 | `rational`、`index_map`、`poly1` |
| 参照単体上 | 標準三角形・標準四面体上の Bernstein 多項式。面積・体積はまだ掛からない | `bpoly`、`inner(bpoly, bpoly)`、`eval`、`range` |
| 物理要素上 | メッシュ中の 1 要素での局所計算。測度、Jacobian、重心座標の勾配が入る | `element_geometry`、`element_op` |
| 領域全体 | 全要素を回って大域自由度へ組み立てた行列・ベクトル・場 | `fe_space`、`rt_space`、`broken_space`、`vfe_space`、`c1_space`、`dirichlet_reduction` |

通常の利用者は「領域全体」の `fe_space` などから始めます。空間クラスの
`stiffness` や `mass` は大域行列を返します。一方、`bpoly` の `inner` は参照単体上の
内積の係数で、物理要素上の積分値にするには要素の測度を掛ける必要があります。

要素の順序と、重複した成分を足す順序は固定されています。同じ入力なら同じ順序で
加算されます。

## スカラー型と行列 policy

要素型 `T` には、`double`、`kv::interval<double>`、厳密な有理数
`vcp::bfem::rational`(`<vcp/bfem/meshgen_rational.hpp>`)などを使えます。
`T` には `T(0)`、`T(1)`、四則演算、コピー、比較(または専用の traits)が必要です。
幾何量を作るときに、要素ごとに 1 回だけ除算を行います。2D C¹ 要素では、辺の長さに
関する除算も使います。

空間クラスは、密行列の policy `P` と疎行列の policy `SP` をテンプレート引数に
取ります。

```cpp
vcp::bfem::fe_space<D, T, P = vcp::mats2<T>, SP = vcp::spmats<T> >
```

荷重ベクトルなどの密な量は `vcp::matrix<T, P>`、大域行列は `vcp::spmatrix<T, SP>`
で返ります。区間型 `T` で疎行列を作る場合は、`SP` に区間疎行列の policy
`vcp::spimats` を渡します([spmatrix.md](spmatrix.md))。

```cpp
typedef kv::interval<double> IT;
vcp::bfem::fe_space<2, IT, vcp::imats<double>, vcp::spimats<double> > fs(Th, 1);
```

疎行列を作らない API だけを呼ぶ場合は、`SP` は実体化されないので既定のままで
使えます。

多項式の係数は `poly1<T>::from_rational` を使うと、厳密な有理数の係数を区間へ
一度だけ変換して作れます。

```cpp
std::vector<std::pair<long long, long long> > a;
a.push_back(std::make_pair(1, 1));   // 1
a.push_back(std::make_pair(0, 1));   // 0 x
a.push_back(std::make_pair(1, 3));   // x^2 / 3
vcp::bfem::poly1<IT> f = vcp::bfem::poly1<IT>::from_rational(a);
```

## メッシュ

### メッシュを直接作る

`mesh<D, T>` は頂点座標と要素の頂点番号を持つ値型です。2D の単位正方形を
2 つの三角形に分ける例です。

```cpp
typedef double T;

std::vector<std::array<T, 2> > vertices;
vertices.push_back(std::array<T, 2>{{T(0), T(0)}});
vertices.push_back(std::array<T, 2>{{T(1), T(0)}});
vertices.push_back(std::array<T, 2>{{T(1), T(1)}});
vertices.push_back(std::array<T, 2>{{T(0), T(1)}});

std::vector<std::array<int, 3> > elements;
elements.push_back(std::array<int, 3>{{0, 1, 2}});
elements.push_back(std::array<int, 3>{{0, 2, 3}});

vcp::bfem::mesh<2, T> msh =
    vcp::bfem::mesh<2, T>::from_lists(vertices, elements);
```

3D の場合は `std::array<T, 3>` と `std::array<int, 4>` を使います。要素の向きは
正規化されません。面積・体積が 0 の要素や、区間型で符号を確定できない要素は、
空間や幾何を作る時点で `degenerate_element` になります。

### メッシュ生成

`generate_mesh` は、多角形領域(非凸・穴あり可)を三角形分割し、すべての辺の
長さが `h` 以下(`len^2 <= h^2`)になるまで一様に細分します。長さは常に 2 乗で
比べるので、平方根を使いません。区間型や有理数型の座標でも同じように使えます。

```cpp
#include <vcp/bfem/meshgen.hpp>

vcp::bfem::polygon_domain<T> dom;   // L 字領域
dom.outer = {{{T(0), T(0)}}, {{T(2), T(0)}}, {{T(2), T(1)}},
             {{T(1), T(1)}}, {{T(1), T(2)}}, {{T(0), T(2)}}};
// dom.holes に穴の多角形を追加できる

vcp::bfem::meshgen_status<2> st;
vcp::bfem::mesh<2, T> Th = vcp::bfem::generate_mesh(dom, T(1) / T(2), st);
```

`meshgen_status` には、境界の頂点・境界の辺・各境界辺が入力のどの線分から
来たか(`facet_source`)などが入ります。オプション `meshgen_options` で、入力の
自己交差の検査(`validate_input`、既定 `true`)と細分の上限回数(`max_refine`、
既定 30)を指定できます。

| 関数 | 内容 |
| --- | --- |
| `generate_mesh(dom, h, status[, opt])` | 2D 多角形領域の一様細分メッシュ |
| `generate_mesh_delaunay(dom, h, ...)` | 細分の前に Lawson flip を行う版(`meshgen_delaunay.hpp`) |
| `generate_mesh_adaptive(dom, h, spec, status, ...)` | newest vertex bisection による局所細分。指定した要素の細分や、角に向けた幾何級数的な細分を `meshgen_adapt_spec` で指定(`meshgen_adapt.hpp`) |
| `generate_mesh(extruded_domain, h, ...)` | 2D の断面を z 方向に押し出した 3D 領域の四面体メッシュ(`meshgen3.hpp`) |

メッシュは、機種によらない 2 進形式で保存・読み込みできます(`meshgen_io.hpp`)。
値は厳密に復元されます。要素型の変換には `convert_mesh` を使います
(`meshgen_convert.hpp`)。

```cpp
vcp::bfem::save_mesh("mesh.bin", Th);
vcp::bfem::mesh<2, T> Th2 = vcp::bfem::load_mesh_file<2, T>("mesh.bin");

vcp::bfem::mesh<2, kv::interval<double> > ITh =
    vcp::bfem::convert_mesh<kv::interval<double> >(Th);
```

## 連続 `P^k`

`fe_space<D, T>` は、2D・3D の連続 Bernstein `P^k` 空間です。

```cpp
const int k = 2;
vcp::bfem::fe_space<2, T> fs(msh, k);

int N  = fs.ndof(k);
auto A = fs.stiffness(k);       // (grad psi_j, grad psi_i)
auto M = fs.mixed_mass(k, k);   // (psi_j, psi_i)
```

非線形項 `f(u)` は `poly1<T>`(1 変数多項式)で与えます。荷重ベクトルと、
Newton 法のヤコビ行列に使う重み付き質量行列を作れます。

```cpp
auto uh = fs.zero_function(k);

// f(x) = 1 + x^2
std::vector<T> coeff;
coeff.push_back(T(1));
coeff.push_back(T(0));
coeff.push_back(T(1));
vcp::bfem::poly1<T> f = vcp::bfem::poly1<T>::from_coeffs(coeff);

auto F  = fs.load(f, uh, k);                        // (f(uh), psi_i)
auto J  = fs.weighted_mass(f.derivative(), uh, k);  // (f'(uh) psi_j, psi_i)
T   nrm = fs.scalar_ff(f, uh);                      // (f(uh), f(uh))
```

`uh.coeffs()` は `N x 1` の `vcp::matrix<T, P>` で、係数を直接代入できます。

```cpp
for (int i = 0; i < uh.coeffs().rowsize(); ++i) {
    uh.coeffs()(i, 0) = T(0);
}
```

点評価は、要素番号と重心座標で行います。同じ関数を高次の空間で扱う場合は
`elevate` を使います。

```cpp
vcp::bfem::bary_point<2, T> lam = {{T(1) / T(3), T(1) / T(3), T(1) / T(3)}};
T val = fs.eval(uh, 0, lam);

auto uh3 = fs.elevate(uh, 3);
```

## 斉次 Dirichlet 縮約

斉次 Dirichlet 条件 `u = 0` は、境界の自由度を消した縮小系を作って扱います。

```cpp
std::vector<int> bdofs = fs.dofs(k).boundary_dofs();
vcp::bfem::dirichlet_reduction<T> red(fs.ndof(k), bdofs);

auto Ar = red.reduce(A);
auto Fr = red.reduce(F);

// Ar xr = Fr を解いた後:
// auto x_full = red.expand(xr);
```

境界の一部だけを指定する場合は、2D では境界辺の ID、3D では境界面の ID を渡します。
返るリストには、その辺の内部の自由度だけでなく端点の頂点の自由度も含まれます。

```cpp
std::vector<int> edges;
edges.push_back(0);
std::vector<int> partial = fs.dofs(k).boundary_dofs(edges);
```

## `RT^k` と broken `P_l`

`rt_space` は Raviart–Thomas 空間、`broken_space` は要素ごとに不連続な `P_l`
空間です。broken 空間の自由度は要素間で共有されず、係数は要素ごとのブロック順に
並びます。フラックスや hypercircle 型の量を作るのに使います。

```cpp
vcp::bfem::rt_space<2, T> rs(msh, 1);
vcp::bfem::broken_space<2, T> bs(msh, 1);

auto P    = rs.mass();                                // (sigma_j, tau_i)
auto Nmat = vcp::bfem::assemble_div_mass(bs, rs);     // (div sigma_j, q_i)
auto rhs  = vcp::bfem::broken_load(bs, f, fs, uh);    // (f(uh), q_i)

auto sig = rs.zero_field();
T e1 = vcp::bfem::flux_error_sq(rs, sig, fs, uh);        // ||grad uh - sig||^2
T e2 = vcp::bfem::div_residual_sq(rs, sig, f, fs, uh);   // ||div sig + f(uh)||^2
T e3 = vcp::bfem::projection_error_sq(bs, f, fs, uh);    // ||f(uh) - Pi_l f(uh)||^2
```

`rt_space::interpolate(provider)` は、入力の場が法線方向に連続であることを仮定
します。連続 `P^k` の関数の勾配は一般に法線方向に連続ではないので、そのまま RT
補間する用途には注意してください。

## Scott–Vogelius

速度をベクトル値の `P^n`、圧力を broken `P_{n-1}` として組みます。

```cpp
const int n = 3;
vcp::bfem::fe_space<2, T> fs_sv(msh, n);
vcp::bfem::vfe_space<2, T> vs(msh, fs_sv);
vcp::bfem::broken_space<2, T> ps(msh, n - 1);

auto Avel = vcp::bfem::assemble_vector_stiffness(vs, n);   // (grad u : grad v)
auto B    = vcp::bfem::assemble_div_velocity(ps, vs, n);   // (div u, q)

auto u    = vs.zero_function(n);
auto divu = vcp::bfem::div_field(ps, vs, u);
T    div2 = vcp::bfem::div_norm_sq(vs, u);
```

`vfe_space` は `fe_space` を所有しないので、参照先の `fe_space` は `vfe_space` より
長く生存させてください。速度の斉次 Dirichlet 条件は、`vfe_space::boundary_dofs` と
`dirichlet_reduction` で扱います。

```cpp
vcp::bfem::dirichlet_reduction<T> vred(vs.ndof(n), vs.boundary_dofs(n));
auto A0 = vred.reduce(Avel);
auto B0 = vred.reduce_cols(B);
```

特異な頂点(2D)・辺(3D)に由来する圧力の拘束は、`sv_pressure_constraints` で
作り、`linear_reduction` で消去します。Alfeld 分割(各単体の重心分割)は
`alfeld_refine(msh)` で作れます。

```cpp
vcp::bfem::sv_pressure_constraints<2, T> cons(msh, n);
vcp::bfem::linear_reduction<T> pred(ps.ndof(), cons.rows());

auto B1 = pred.reduce_rows(B0);
```

鞍点系のブロックの組み立て、定数圧力モードの扱い、線形方程式の解法は
利用者側で行います。

## CR1

`cr1_space<D, T>` は、2D・3D の Crouzeix–Raviart 空間(1 次)です。次数は 1 に
固定なので、コンストラクタに次数の引数はありません。自由度は面(2D では辺)の
平均値です。

```cpp
vcp::bfem::cr1_space<2, T> cr(msh);

auto K  = cr.stiffness();
auto Mc = cr.mass();
std::vector<int> b = cr.boundary_dofs();
T h2 = cr.max_edge_length_sq();   // 最大の辺の長さの 2 乗
```

`max_edge_length_sq()` はどの `T` でも厳密に計算されます(区間型では真の値を
包含します)。`max_edge_length()` は平方根を取る版です。

ベクトル値の `vcr1_space<D, T>`(成分ごとの直和)と、離散的に発散が零となる部分
空間の基底を作る `vcr1_divfree.hpp` もあります。これらは Stokes 方程式の誤差定数
(`stokes_constants.hpp`)で使われます。

CR1 は `H^1` に適合しない空間です。`space_traits<S>::h1_conforming`
(`<vcp/bfem/space_traits.hpp>`)で、空間が `H^1` に適合するかを調べられます。

## 2D C¹

`c1_space<2, T>` は 2D の C¹ Argyris 型空間で、次数は `k >= 5` です。

```cpp
const int c1k = 5;
vcp::bfem::c1_space<2, T> cs(msh, c1k);

auto A = cs.stiffness(c1k);
auto M = cs.mixed_mass(c1k, c1k);

auto u  = cs.zero_function(c1k);
auto L  = cs.laplacian_matrix(c1k);   // (Delta u, Delta v)
auto H  = cs.hessian_matrix(c1k);     // (D^2 u : D^2 v)
T    r2 = cs.laplacian_residual_sq(f, u);

vcp::bfem::bary_point<2, T> lam = {{T(1) / T(3), T(1) / T(3), T(1) / T(3)}};
std::array<T, 2> g = cs.eval_grad(u, 0, lam);
std::array<T, 3> h = cs.eval_hess(u, 0, lam);   // xx, xy, yy
```

境界条件は一般の線形拘束として作り、`linear_reduction` で消去します。

```cpp
std::vector<vcp::bfem::sv_constraint_q> rows =
    vcp::bfem::c1_boundary(cs, c1k, vcp::bfem::c1_bc_dirichlet);

vcp::bfem::linear_reduction<T> cred(cs.ndof(c1k), rows);
auto Ar = cred.reduce(A);
```

`c1_boundary` は、境界の方向の分類と拘束の係数を作るために `c1_coord_traits<T>` を
使います。`vcp::bfem::rational` には既定の特殊化がありますが、`double` などでは
利用者が「符号の判定」と「座標差の有理数化」を定義してください。次の例は、座標差が
整数の場合だけに使えるものです。一般の `double` の座標では、丸めた実数をどの
有理数として扱うかを明示的に決める必要があります。

```cpp
namespace vcp {
namespace bfem {

template <>
struct c1_coord_traits<double> {
    static int sign(const double& x) { return detail::c1_sign_certified(x); }
    static detail::rational to_rational(const double& x) {
        return detail::rational(static_cast<long long>(x));
    }
};

} // namespace bfem
} // namespace vcp
```

## 誤差定数

`vcp/bfem/constants/` は、有限要素法の誤差定数・固有値の下界・Sobolev 埋め込み
定数などを精度保証付きで計算します。名前空間は `vcp::bfem::constants` です。
記号と式番号は『精度保証付き数値計算の基礎』の第 8 章に合わせています。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/bfem/bound_traits_kv.hpp>
#include <vcp/bfem/geometry_traits_kv.hpp>
#include <vcp/bfem/convert_traits_kv.hpp>
#include <vcp/bfem/meshgen.hpp>
#include <vcp/bfem/constants/poisson_constants.hpp>
#include <vcp/bfem/constants/sobolev_constants.hpp>

typedef kv::interval<double> IT;

vcp::bfem::polygon_domain<IT> dom;
dom.outer = {{{IT(0), IT(0)}}, {{IT(1), IT(0)}}, {{IT(1), IT(1)}}, {{IT(0), IT(1)}}};
vcp::bfem::meshgen_status<2> st;
vcp::bfem::mesh<2, IT> Th = vcp::bfem::generate_mesh(dom, IT(1) / IT(2), st);

IT C0h = vcp::bfem::constants::l2_projection_error_constant(Th);          // C_0 h
IT Ch  = vcp::bfem::constants::ritz_projection_error_constant_h01(Th, 1); // C_h(P^1)
IT mea = vcp::bfem::constants::domain_measure(Th);                        // |Omega|
```

要素型 `T` は kv の区間型です。関数名の約束は次の通りです。

- `_h01` が付く関数は、斉次 Dirichlet 条件の空間 `H^1_0` を仮定します。付かない関数
  (`l2_projection_*` など)は境界条件によりません。
- `_sq`・`_squared` が付く関数は 2 乗の値を返します。
- 精度保証付きの上界として使うのは、戻り値の `.upper()` です。

主な関数は次の通りです。

| 関数 | 内容 | 次元 |
| --- | --- | --- |
| `l2_projection_error_constant(Th)` | 区分定数への射影の誤差定数 C_0 h(式 (8.19)(8.20)。次数 k によらない) | 2D |
| `l2_projection_element_constant(Th, e)` | 要素 `e` ごとの同じ定数 | 2D |
| `ritz_projection_error_constant_h01(Th, k)` | `P^k` の Ritz 射影の誤差定数 C_h(定理 8.3、式 (8.25)) | 2D |
| `hypercircle_kappa_h01_squared(Th, k)` | κ_h の 2 乗(式 (8.21)) | 2D・3D |
| `ritz_projection_constants_h01(Th, k)` | 定理 8.3 の 3 つの定数を 1 回の計算でまとめて返す | 2D |
| `cr_interpolation_constant<D, T>()` | CR 補間の誤差定数 c_D | 2D・3D |
| `cr_projection_error_constant_sq(V)` | CR1 空間 `V` の射影誤差定数の 2 乗 | 2D・3D |
| `eigenvalue_lower_bounds_h01<D, T>(Th, nev, ...)` | Dirichlet Laplace 作用素の固有値の下界(CR1 空間による) | 2D・3D |
| `poincare_constant_h01_bound(lambda1_lower)` | 第 1 固有値の下界から Poincaré 定数の上界 | — |
| `aubin_talenti_constant<D, T>(pnum, pden)` | 全空間の Sobolev 不等式の最良定数(Aubin–Talenti) | — |
| `sobolev_embedding_constant_h01<D, T>(pnum, pden, mea, ...)` | `H^1_0` から `L^p` への埋め込み定数 | — |
| `domain_measure(Th)` | 領域の測度 \|Ω\| | 2D・3D |
| `stokes_*`(`stokes_constants.hpp`) | Stokes 方程式(ベクトル値・発散零の空間)の同種の定数と固有値の下界 | 2D・3D |

要素ごとの `L^2` 射影の誤差定数については、計算済みの値の辞書
(`constants/dict/`)と、その場で計算する経路があります。

内部の区間行列計算に使う密行列の policy `DP` は、既定では `vcp::imats<double>` で、
外部 BLAS を必要としません。BLAS がある場合は `DP = vcp::pidblas` を渡すと速くなります。
`vcp::imats<double, vcp::pdblas>` は点の計算だけを速くするもので、区間どうしの積は
速くならないため、同じ効果はありません。

```cpp
IT Ch = vcp::bfem::constants::ritz_projection_error_constant_h01<IT, vcp::pidblas>(Th, 1);
```

## 図示用の出力

`output_uh_for_graphics`(`<vcp/bfem/graphics.hpp>`)は、有限要素関数を要素ごとの
重心座標の格子で評価し、点の行列とセルの接続行列の組を返します。matplotlib・
MATLAB・VTK 系のツールでそのまま使える形です。

```cpp
#include <vcp/bfem/graphics.hpp>

vcp::bfem::graphics_output<2, T, vcp::mats2<T> > g =
    vcp::bfem::output_uh_for_graphics(fs, uh, 2);
std::cout << g.points;   // 行 = 点。列 = 座標 D 列 + 値 ncomp 列
std::cout << g.cells;    // 行 = 小単体。成分は points の行番号(0 始まり)
```

| メンバ | 内容 |
| --- | --- |
| `points` | `npoint x (D + ncomp)` の行列。各行がサンプル点 1 つで、物理座標 D 列と値 `ncomp` 列 |
| `cells` | `ncell x (D + 1)` の `vcp::matrix<int>`。各行が小単体 1 つで、成分は `points` の行番号(0 始まり) |
| `ncomp` | 値の成分数(スカラーは 1、Scott–Vogelius の速度は D) |
| `div` | 生成に使った分割数 |

- `fe_space`、`c1_space`、`broken_space`、`vfe_space` に対応しています。RT 空間は
  対象外です。
- 要素の番号のリスト `elems` を渡す版 `output_uh_for_graphics(space, u, elems, div)`
  もあり、行はリストの順に並びます。
- 要素あたり `C(div + D, D)` 点、`div^D` 個のセルを出力します。隣り合う要素が
  共有する辺・面の上の点は、要素ごとに重複して出力されます(不連続な場の跳びを
  保つため)。後で節点を併合する場合は、許容誤差付きで座標を照合してください。
- `div` の既定は 1(頂点のみ)です。高次要素の要素内の形(特に `c1_space`)や、
  3D の粗いメッシュを見たいときに `div` を上げます。
- MATLAB では添字に 1 を足します(`trisurf(cells + 1, x, y, u)`)。

区間型では、座標と値の両方が区間で返り、各成分は対応する厳密値を包含します。
これは「その格子点での座標と値の包含」であり、座標の区間内の任意の点で値が値の
区間に入るという主張ではありません。有理数型では座標と値は厳密です。

出力を ParaView で表示するための Python のツールが `test_PDE/bfem_viewer/` に
あります(使い方は同じディレクトリの `README_visualization.md`)。

## OpenMP

大域行列の組み立ての一部は OpenMP で並列化されています。プログラムを `-fopenmp`
付きでコンパイルすれば有効になります。スレッド数は環境変数 `OMP_NUM_THREADS`、
または利用者のコードで `<omp.h>` を include して `omp_set_num_threads` で指定します。

bfem の OpenMP だけを無効にする場合は `-DVCP_BFEM_NOMP` を付けます。`-DVCP_NOMP` を
付けた場合も、bfem の OpenMP は無効になります。

| 空間 | 並列化される主な関数 |
| --- | --- |
| 2D・3D `P^k` | `stiffness`、`mixed_mass`、`weighted_mass` |
| 2D・3D `RT^k` | `rt_space::mass`、`assemble_div_mass`、`assemble_cross_grad` |
| broken `P_l` | `broken_space::mass` |
| Scott–Vogelius | `assemble_vector_stiffness`、`assemble_div_velocity`、`assemble_advection`、`assemble_advection_derivative` |
| 2D C¹ | `stiffness`、`mixed_mass`、`weighted_mass`、`laplacian_matrix`、`hessian_matrix`、`laplacian_mixed` |

並列化は 1 回の組み立ての呼び出しの中で要素のループを分けるものです。1 つの
空間のオブジェクトを、利用者の複数のスレッドから同時に呼ぶ使い方は想定して
いません。スレッドごとの COO バッファを使い、成分を足す順序を固定しています。

## 例外

bfem が投げる例外は、すべて `vcp::error` から派生します
([error_handling.md](error_handling.md))。bfem に固有のクラスは次の通りです。

| 例外 | 基底 | 主な条件 |
| --- | --- | --- |
| `vcp::bfem::degenerate_element` | `vcp::verification_error` | 要素の測度が 0、または区間型で向きの符号を確定できない |
| `vcp::bfem::sv_indeterminate_singularity` | `vcp::verification_error` | 区間座標で Scott–Vogelius の特異性を判定できない |
| `vcp::bfem::c1_indeterminate_corner` | `vcp::verification_error` | 区間座標で C¹ の境界の角を分類できない |
| `vcp::bfem::meshgen_degenerate` | `vcp::verification_error` | メッシュ生成で、零または符号を確定できない値が出た |
| `vcp::bfem::meshgen_error` | `vcp::invalid_argument` | メッシュ生成の入力の誤り(自己交差、頂点が 3 未満の多角形、外側の境界に囲まれない穴など) |
| `vcp::bfem::meshgen_limit` | `vcp::numerical_error` | 細分の回数が `max_refine` を超えた |
| `vcp::bfem::meshgen_convert_error` | `vcp::domain_error` | メッシュの型変換で NaN・非有限値など変換できない値があった |
| `vcp::bfem::meshgen_io_error` | `vcp::io_error` | メッシュの保存・読み込みの失敗、ファイルの破損 |

この他に、次数の範囲外、係数のサイズの不一致、自由度の範囲外、異なるメッシュの
空間の混用などは `vcp::invalid_argument` になります。

## コンパイル例

bfem はヘッダだけで使えます。`vcp/` を含むディレクトリを include path に指定します
([install.md](install.md))。

```bash
g++ -std=c++11 -O3 -DNDEBUG -I/path/to/libs example.cpp -fopenmp
```

区間型・誤差定数を使う場合は、kv の include path と MPFR も必要です。

```bash
g++ -std=c++11 -O3 -DNDEBUG -DKV_FASTROUND -I/path/to/libs example.cpp -lmpfr -fopenmp
```

## よくある失敗

| 症状 | 原因と対処 |
| --- | --- |
| `degenerate_element` が投げられる | 要素がつぶれている、または区間座標で向きが確定しない。メッシュや座標の区間を確認する |
| 係数のサイズの不一致 | `function_from_coeffs` に `ndof(m) x 1` でない行列を渡している |
| 区間型で疎行列がコンパイルできない | 疎行列の policy `SP` が区間型に対応していない。`vcp::spimats` を渡す |
| OpenMP が効かない | `-fopenmp` が付いていない、または `VCP_BFEM_NOMP` / `VCP_NOMP` を定義している |
| `linear_reduction` が失敗する | 拘束の行が従属している、ピボットがない、または自由度の番号が範囲外 |
| RT 補間の後に期待した誤差が出ない | 入力の場が法線方向に連続でない可能性がある |
| C¹ で `m < k` が拒否される | `c1_space` は基準の次数以上の次数だけを扱う |

## 関数リファレンス

### `mesh<D, T>`(`<vcp/bfem/mesh.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `static mesh from_lists(vertices, elements)` | 頂点座標の列と要素の頂点番号から作る。頂点番号の範囲を検査する |
| `int num_vertices() const` | 頂点数 |
| `int num_elements() const` | 要素数 |
| `const std::array<T, D>& vertex(int v) const` | 頂点座標 |
| `const std::array<int, D + 1>& element(int e) const` | 要素の頂点番号 |

### `element_geometry<D, T>`(`<vcp/bfem/geometry.hpp>`)

1 つの三角形・四面体について、参照単体から物理要素への写像を持ちます。

| 関数 | 内容 |
| --- | --- |
| `static element_geometry from_vertices(verts)` | 単体の幾何量を作る |
| `measure() const` | 面積・体積 |
| `inv_absdet() const` | `1 / abs(det B)` |
| `orientation() const` | det の符号 |
| `grad_lambda(i, d) const` | 重心座標 `lambda_i` の物理座標方向の微分 |
| `edge_matrix(r, c) const` | アフィン写像の辺行列 `B` |
| `vertices() const` | 要素の頂点座標 |

### `bpoly<D, T>`(`<vcp/bfem/bpoly.hpp>`)

参照単体上の Bernstein 多項式です。係数の順序は `index_map<D>` の順序です。

| 関数 | 内容 |
| --- | --- |
| `bpoly()` | 0 次の零多項式 |
| `static constant(s)` | 定数多項式 |
| `static zero(n)` | `P^n` の零多項式 |
| `static from_coeffs(n, c)` | 係数の列から作る |
| `degree() const`、`size() const` | 次数、係数の数 |
| `coeff(r)`、`coeffs() const` | 係数 |
| `operator+=`、`operator-=` | 同じ次数の加減算 |
| `operator*=`、`add_scalar(s)` | スカラー倍、定数の加算 |

| 自由関数 | 内容 |
| --- | --- |
| `add`、`sub`、`operator+`、`operator-` | 次数が異なれば低い方を次数上げしてから計算 |
| `mul`、`operator*` | 積 |
| `scale` | スカラー倍 |
| `elevate(u, m)` | 次数上げ(`m >= u.degree()`) |
| `dlambda(u, i)` | 重心座標方向の微分 |
| `inner(u, v)` | 参照単体上の `L^2` 内積の係数。物理要素上の積分値は要素の測度を掛ける |
| `eval(u, lam)` | de Casteljau 法による点評価 |
| `add_into`、`sub_into`、`mul_into`、`scale_into`、`elevate_into`、`dlambda_into` | 出力先を再利用する版 |

### `poly1<T>`(`<vcp/bfem/poly1.hpp>`)

1 変数多項式 `f(x) = sum a_k x^k` です。

| 関数 | 内容 |
| --- | --- |
| `static from_coeffs(a)` | `T` の係数から作る |
| `static from_rational(num_den)` | 有理数の係数から、`T` へ一度だけ変換して作る |
| `degree() const`、`coeff(k) const` | 次数、係数 |
| `derivative() const` | 微分 |
| `compose(f, u)` | `bpoly` に `f(u)` を合成する |
| `compose_into(dst, f, u, workspace)` | 作業領域を再利用する版 |

### 値域の包含(`<vcp/bfem/refine.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `range(u)` | Bernstein 係数の凸包性による値域の包含 |
| `restrict_to(u, V)` | 親の単体の重心座標で与えた部分単体への厳密な制限 |
| `range_refined(u, depth)` | 2D は red 細分、3D は Bey 細分(`<vcp/bfem/d3/bey_table.hpp>`)で値域を精密化 |

算術型でない `T` で `range` を使うには `bound_traits<T>` の特殊化が必要です
(kv の区間型は `bound_traits_kv.hpp`)。

### `fe_space<D, T, P, SP>`(`<vcp/bfem/fe_space.hpp>`)

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `fe_space(msh, n)` | | 基準の次数 `n >= 1` で作る。全要素の幾何量も作る |
| `base_degree() const` | `int` | 基準の次数 |
| `num_elements() const` | `int` | 要素数 |
| `dofs(m)` | `const dofmap<D>&` | 次数 `m` の自由度の写像 |
| `ndof(m)` | `int` | 次数 `m` の大域自由度数 |
| `zero_function(m)` | `fe_function<D, T, P>` | 零関数 |
| `function_from_coeffs(m, c)` | `fe_function<D, T, P>` | 係数ベクトルから関数を作る |
| `stiffness(m)` | `spmatrix` | `(grad psi_j, grad psi_i)` |
| `mixed_mass(a, b)` | `spmatrix` | `(phi_j^b, phi_i^a)`。長方形も可 |
| `load(f, uh, m)` | `matrix<T, P>` | `(f(uh), psi_i)` |
| `weighted_mass(fprime, uh, m)` | `spmatrix` | `(f'(uh) psi_j, psi_i)` |
| `convection(b, m)` | `spmatrix` | `(b . grad psi_j, psi_i)`。`b` は関数の配列または定数の配列 |
| `inner(u, v)` | `T` | 全領域の `(u, v)` |
| `scalar_ff(f, uh)` | `T` | 全領域の `(f(uh), f(uh))` |
| `elevate(u, m)` | `fe_function` | 関数を変えずに次数 `m` へ上げる |
| `dual(u, m)` | `dual_vector` | 大域モーメント `(u, psi_i)` を作る |
| `eval(u, e, lam)` | `T` | 要素 `e` 上で、重心座標 `lam` の点の値 |

`fe_function`(`<vcp/bfem/fe_function.hpp>`)は `degree()` と、大域の Bernstein 係数
`coeffs()`(直接編集できる)を持ちます。`dual_vector` は `test_degree()` と
モーメントのベクトル `values()` を持ちます。`dual_vector` から `fe_function` へ
戻す関数はありません(質量行列の逆が必要になるため)。

### `dofmap<D>`

`fe_space::dofs(m)` で取得します。

| 関数 | 内容 |
| --- | --- |
| `ndof()`、`degree()`、`local_size()`、`num_elements()` | 自由度数、次数、要素あたりの自由度数、要素数 |
| `global_dof(e, r)` | 要素 `e` の局所番号 `r` の大域番号 |
| `dof_sign(e, r)` | `P^k` では常に `+1` |
| `boundary_dofs()` | 全境界の自由度 |
| `boundary_dofs(edges_or_faces)` | 2D では辺、3D では面を指定した境界の自由度 |
| `boundary_edge_ids()` / `boundary_face_ids()` | 2D の境界辺 / 3D の境界面の ID |

### `dirichlet_reduction<T, P, SP>`(`<vcp/bfem/dirichlet.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `dirichlet_reduction(full_size, bdofs)` | 拘束する自由度のリストから作る。重複・範囲外は例外 |
| `full_size()`、`reduced_size()` | 元の自由度数、縮約後の自由度数 |
| `to_reduced(full)` | 元の番号から縮約後の番号(拘束した自由度は `-1`) |
| `to_full(reduced)` | 縮約後の番号から元の番号 |
| `reduce(A)`、`reduce(v)` | 正方行列の拘束した行・列、ベクトルの拘束した成分を消す |
| `reduce_rows(A)`、`reduce_cols(A)` | 行だけ、列だけ縮約する |
| `expand(vr)` | 縮約したベクトルを全自由度に戻す(拘束した成分は 0) |

### `rt_space<D, T, P, SP>`(`<vcp/bfem/rt/rt_space.hpp>`)

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `rt_space(msh, k)` | | `k >= 0` の RT 空間 |
| `order()`、`ndof()`、`num_elements()`、`num_vertices()` | `int` | 次数、自由度数、要素数、頂点数 |
| `dofs()` | | RT の自由度の写像 |
| `geometry(e)` | `element_geometry` | 要素の幾何量 |
| `zero_field()`、`field_from_coeffs(c)` | `rt_field` | 零の場、係数からの場 |
| `mass()` | `spmatrix` | `(sigma_j, tau_i)` |
| `interpolate(provider)` | `rt_field` | 物理成分で与えた場を RT の自由度へ補間。provider は要素ごとに、2D では `std::pair<bpoly, bpoly>`、3D では `std::array<bpoly, 3>` を返す |

### `broken_space<D, T, P, SP>`(`<vcp/bfem/rt/broken_space.hpp>`)

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `broken_space(msh, l)` | | `l >= 0` の broken 空間 |
| `order()`、`ndof()`、`local_size()` | `int` | 次数、自由度数、要素あたりの係数の数 |
| `dofs()`、`geometry(e)` | | ブロック型の自由度の写像、要素の幾何量 |
| `zero_field()`、`field_from_coeffs(c)` | `broken_field` | 零の場、係数からの場 |
| `mass()` | `spmatrix` | ブロック対角の質量行列 |
| `inner(u, v)` | `T` | `L^2` 内積 |

### `rt_assemble.hpp` の自由関数(`<vcp/bfem/rt/rt_assemble.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `assemble_div_mass(bs, rs)` | 行 = broken、列 = RT。`(div sigma_j, q_i)` |
| `assemble_cross_grad(rs, fs, m)` | 行 = RT、列 = `P^m`。`(sigma_i, grad psi_j)` |
| `assemble_mixed_mass(fs, m, bs)` | 行 = `P^m`、列 = broken。`(phi_i, q_j)` |
| `broken_load(bs, f, fs, uh)` | broken の試験関数での `(f(uh), q_i)` |
| `flux_error_sq(rs, sig, fs, u)` | `‖grad u − sig‖^2` |
| `div_residual_sq(rs, sig, f, fs, u)` | `‖div sig + f(u)‖^2` |
| `projection_error_sq(bs, f, fs, uh)` | `‖f(uh) − Pi_l f(uh)‖^2` |

### `vfe_space<D, T, P, SP>`(`<vcp/bfem/sv/vfe_space.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `vfe_space(msh, scalar_fe_space)` | ベクトル値空間を作る(`fe_space` は所有しない) |
| `scalar()` | 参照先の `fe_space` |
| `ndof(m)` | `D * scalar.ndof(m)` |
| `global_dof(d, i, m)` | 成分 `d`、スカラーの自由度 `i` の大域番号 |
| `zero_function(m)`、`function_from_coeffs(m, c)` | 零関数、係数からの関数 |
| `boundary_dofs(m)` | 全成分の境界の自由度 |
| `component(u, d)`、`set_component(u, d, f)` | 成分の取り出し(コピー)と書き戻し |

### `sv_assemble.hpp` の自由関数(`<vcp/bfem/sv/sv_assemble.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `assemble_vector_stiffness(vs, m)` | `(grad u : grad v)` |
| `assemble_div_velocity(bs, vs, m)` | 行 = broken の圧力、列 = 速度。`(div u, q)` |
| `assemble_advection(vs, w, m)` | 固定した速度 `w` の移流行列 `(w . grad u, v)` |
| `assemble_advection_derivative(vs, w, m)` | Newton 法の微分の側 `((u . grad) w, v)` |
| `advection_vector(vs, a, b, m)` | `((a . grad) b, phi_i)` |
| `advection_scalar(vs, a, b, g)` | `((a . grad) b, g)`。`g` は `vfe_function` または `fe_function` の勾配 |
| `div_field(bs, vs, u)` | `div u` を broken の場として返す |
| `div_norm_sq(vs, u)` | `‖div u‖^2` |

### 圧力の拘束と線形縮約(`sv_rows.hpp`、`linear_reduction.hpp`、`alfeld.hpp`)

| 関数 | 内容 |
| --- | --- |
| `alfeld_refine(msh)` | 各単体の重心分割 |
| `sv_pressure_constraints<D, T>(msh, n)` | Scott–Vogelius の圧力の拘束の行を作る |
| `rows()` | 独立にした整数係数の拘束の行 |
| `num_singular_vertices()` / `num_singular_edges()` | 2D の特異頂点 / 3D の特異辺の数 |
| `hint(provenance)` | 既知の十分条件に基づく観測情報 |
| `linear_reduction<T, P, SP>(full_size, rows)` | 一般の線形拘束を消去する縮約 |
| `linear_reduction::reduce`、`reduce_rows`、`reduce_cols`、`expand` | 行列・ベクトルの縮約と復元 |

### `cr1_space<D, T, P, SP>`(`<vcp/bfem/cr1/cr1_space.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `cr1_space(msh)` | 1 次の CR 空間(次数の引数はない) |
| `ndof()`、`num_elements()`、`num_vertices()` | 自由度数、要素数、頂点数 |
| `dofs()`、`boundary_dofs()` | 自由度の写像、境界の自由度 |
| `stiffness()`、`mass()` | 剛性行列、質量行列 |
| `max_edge_length_sq()` | 最大の辺の長さの 2 乗(どの `T` でも厳密・包含) |
| `max_edge_length()` | その平方根 |

`vcr1_space<D, T>`(`<vcp/bfem/cr1/vcr1_space.hpp>`)はベクトル値版で、大域番号は
`d * n + i`(成分 `d`、スカラーの自由度 `i`、`n` はスカラーの自由度数)です。
`stiffness()`・`mass()` はブロック対角です。

### `c1_space<2, T, P, SP>`(`<vcp/bfem/c1/c1_space.hpp>`)

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `c1_space(msh, k)` | | 基準の次数 `k >= 5` で作る |
| `base_degree()` | `int` | 基準の次数 |
| `num_elements()`、`num_vertices()`、`num_edges()` | `int` | メッシュの情報 |
| `mesh_ref()`、`geometry(e)` | | 元のメッシュ、要素の幾何量 |
| `edge_inv_tsq(ed)` | `T` | 辺の長さの 2 乗の逆数 |
| `dofs(m)`、`ndof(m)` | | 自由度の写像(`m >= base_degree()`)、自由度数 |
| `zero_function(m)`、`function_from_coeffs(m, c)` | `c1_function` | 零関数、係数からの関数 |
| `stiffness(m)`、`mixed_mass(a, b)`、`weighted_mass(fprime, uh, m)` | `spmatrix` | 剛性、質量、重み付き質量 |
| `load(f, uh, m)` | `matrix<T, P>` | 荷重 |
| `inner(u, v)`、`scalar_ff(f, uh)` | `T` | `L^2` 内積、`(f(uh), f(uh))` |
| `eval(u, e, lam)`、`eval_grad(u, e, lam)`、`eval_hess(u, e, lam)` | | 値、物理勾配、Hessian(`xx, xy, yy`) |
| `laplacian_matrix(m)`、`hessian_matrix(m)` | `spmatrix` | `(Delta u, Delta v)`、`(D^2 u : D^2 v)` |
| `laplacian_mixed(m, l)` | `spmatrix` | 行 = broken `P_l`、列 = C¹ |
| `laplacian_load(f, uh, m)` | `matrix<T, P>` | `(f(uh), Delta phi_i)` |
| `laplacian_residual_sq(f, uh)` | `T` | `‖Delta uh + f(uh)‖^2` |
| `laplacian_field(u, bs)` | `broken_field` | `Delta u` を broken `P_{m−2}` で返す |
| `interpolate(provider, m)` | `c1_function` | 要素ごとの多項式から補間 |
| `elevate(u, m)` | `c1_function` | 次数上げ |

`c1_space::dofs(m)` で得る `c1_dofmap` は、`vertex_dof(v, c)`(`c` = 値、dx、dy、
dxx、dxy、dyy)、`edge_trace_dof(ed, i)`、`edge_nd_dof(ed, j)`、`interior_dof(e, r)`、
`boundary_edge_ids()`、`edge_verts(ed)` を持ちます。

### C¹ の境界拘束(`<vcp/bfem/c1/c1_reduce.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `c1_boundary(sp, m, c1_bc_dirichlet)` | `u = 0` |
| `c1_boundary(sp, m, c1_bc_simply_supported)` | 単純支持の本質境界条件(現在の実装では Dirichlet と同じ) |
| `c1_boundary(sp, m, c1_bc_clamped)` | `u = 0` かつ `du/dnu = 0` |
| `c1_boundary(sp, m, kind, edges)` | 指定した境界辺だけに拘束を作る |

戻り値は有理数係数の `sv_constraint_q` の列で、`linear_reduction` に渡して縮約します。
