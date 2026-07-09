# bfem ユーザガイド

この文書は、`vcp/bfem` を初めて使う人向けの最短手順です。API の詳細は
`docs/bfem_api_reference.md`、全体像は `docs/bfem_overview.md` を参照してください。

## 0. どこの話をしているか

bfem を読むときは、次の区別を常に意識してください。

| 層 | 何を扱うか | このガイドで出る例 |
|---|---|---|
| 領域非依存 | メッシュと無関係な多項式・有理数・テーブル | `poly1<T> f`, `from_rational` |
| 参照単体上 | 標準三角形・標準四面体上の Bernstein 多項式 | `bpoly`, `bary_point`, `eval` の重心座標 |
| 物理要素上 | メッシュ中の 1 つの要素での局所計算 | `geometry(e)`, `eval(u,e,lam)` の指定要素 |
| 領域全体 | 全要素を組み立てた空間・行列・ベクトル | `fe_space`, `stiffness`, `load`, `dirichlet_reduction` |

通常の利用者は「領域全体」の `fe_space` などから始めます。`bpoly` は、要素ごとの補間 provider を書く、
局所多項式を直接調べる、区間評価を細かく制御する、といった場面で使います。

## 1. 基本の include

2D/3D の連続 `P^k` 要素だけなら、まず以下で足ります。

```cpp
#include <array>
#include <vector>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/poly1.hpp>
```

RT、Scott-Vogelius、C^1 を使う場合は、それぞれ必要なヘッダを追加します。

```cpp
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_assemble.hpp>

#include <vcp/bfem/sv/vfe_space.hpp>
#include <vcp/bfem/sv/sv_assemble.hpp>
#include <vcp/bfem/sv/sv_rows.hpp>
#include <vcp/bfem/sv/linear_reduction.hpp>

#include <vcp/bfem/c1/c1_space.hpp>
#include <vcp/bfem/c1/c1_reduce.hpp>
```

## 2. メッシュを作る

この節は「領域全体の入力データ」の話です。まだ自由度番号や行列はありません。

2D の単位正方形を 2 三角形に分ける例です。

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

3D の場合は `std::array<T,3>` と `std::array<int,4>` を使います。

## 3. 連続 `P^k` で行列を作る

この節は「領域全体」の話です。`fs.stiffness(k)` や `fs.load(...)` は、すべての要素を回って
大域自由度番号へ散布した結果を返します。

```cpp
const int k = 2;
vcp::bfem::fe_space<2, T> fs(msh, k);

int N = fs.ndof(k);
auto A = fs.stiffness(k);       // 剛性行列
auto M = fs.mixed_mass(k, k);   // 質量行列
```

荷重ベクトルや Newton 用の重み付き質量には、`poly1<T>` で非線形項を渡します。

```cpp
auto uh = fs.zero_function(k);

// f(x) = 1 + x^2
std::vector<T> coeff;
coeff.push_back(T(1));
coeff.push_back(T(0));
coeff.push_back(T(1));
vcp::bfem::poly1<T> f = vcp::bfem::poly1<T>::from_coeffs(coeff);

auto F = fs.load(f, uh, k);
auto J = fs.weighted_mass(f.derivative(), uh, k);
T nrm = fs.scalar_ff(f, uh);
```

`uh.coeffs()` は `N x 1` の `vcp::matrix<T,P>` です。係数を直接代入できます。

```cpp
for (int i = 0; i < uh.coeffs().rowsize(); ++i) {
    uh.coeffs()(i, 0) = T(0);
}
```

## 4. 境界条件を入れる

この節も「領域全体」の話です。境界条件は参照単体上の多項式ではなく、大域自由度番号の集合として扱います。

斉次 Dirichlet 条件は、境界自由度を消して縮小系を作ります。

```cpp
std::vector<int> bdofs = fs.dofs(k).boundary_dofs();
vcp::bfem::dirichlet_reduction<T> red(fs.ndof(k), bdofs);

auto Ar = red.reduce(A);
auto Fr = red.reduce(F);

// 線形ソルバで Ar xr = Fr を解いた後:
// auto x_full = red.expand(xr);
```

部分境界だけを指定したいときは、2D では境界辺 ID、3D では境界面 ID を指定します。

```cpp
std::vector<int> edges;
edges.push_back(0);
std::vector<int> partial = fs.dofs(k).boundary_dofs(edges);
```

返るリストには、その辺の内点自由度だけでなく端点頂点も含まれます。不要な自由度は
利用者側でリストから取り除きます。

## 5. 点評価と次数上げ

点評価は「領域全体の関数を、指定した物理要素へ制限して評価する」操作です。
重心座標 `lam` 自体は参照単体上の座標ですが、`e` によりどの物理要素上の値かが決まります。

点評価は要素番号と重心座標で行います。

```cpp
vcp::bfem::bary_point<2, T> lam = {{T(1) / T(3), T(1) / T(3), T(1) / T(3)}};
T val = fs.eval(uh, 0, lam);
```

同じ関数を高次空間で扱いたいときは `elevate` を使います。

```cpp
auto uh3 = fs.elevate(uh, 3);
```

## 6. RT と broken 空間を使う

この節は「領域全体」の話です。RT 場も broken 場も大域係数ベクトルを持ちます。
ただし broken 空間の自由度は要素間で共有されず、係数は要素ブロック順です。

RT 空間ではフラックスや hypercircle 型の量を組みます。

```cpp
vcp::bfem::rt_space<2, T> rs(msh, 1);
vcp::bfem::broken_space<2, T> bs(msh, 1);

auto P = rs.mass();                              // (sigma, tau)
auto Nmat = vcp::bfem::assemble_div_mass(bs, rs); // (div sigma, q)
auto rhs = vcp::bfem::broken_load(bs, f, fs, uh);
```

RT 場 `sig` がある場合、残差型スカラーを直接計算できます。

```cpp
auto sig = rs.zero_field();
T e1 = vcp::bfem::flux_error_sq(rs, sig, fs, uh);
T e2 = vcp::bfem::div_residual_sq(rs, sig, f, fs, uh);
T e3 = vcp::bfem::projection_error_sq(bs, f, fs, uh);
```

`rt_space::interpolate(provider)` は、入力場が法線連続であることを仮定します。
連続 `P^k` 関数の勾配は一般に法線連続ではないため、そのまま RT 補間する用途には注意が必要です。

## 7. Scott-Vogelius 部品を使う

この節は「領域全体」の話です。特異性検出はメッシュ幾何を見ますが、拘束行は圧力の大域自由度に対して作られます。

Scott-Vogelius では、速度をベクトル値 `P^n`、圧力を broken `P_{n-1}` として組みます。

```cpp
const int n = 3;
vcp::bfem::fe_space<2, T> fs_sv(msh, n);
vcp::bfem::vfe_space<2, T> vs(msh, fs_sv);
vcp::bfem::broken_space<2, T> ps(msh, n - 1);

auto Avel = vcp::bfem::assemble_vector_stiffness(vs, n);
auto B = vcp::bfem::assemble_div_velocity(ps, vs, n);

auto u = vs.zero_function(n);
auto divu = vcp::bfem::div_field(ps, vs, u);
T div2 = vcp::bfem::div_norm_sq(vs, u);
```

速度の斉次 Dirichlet は `vfe_space::boundary_dofs` と `dirichlet_reduction` を使います。

```cpp
vcp::bfem::dirichlet_reduction<T> vred(vs.ndof(n), vs.boundary_dofs(n));
auto A0 = vred.reduce(Avel);
auto B0 = vred.reduce_cols(B);
```

圧力の特異拘束が必要な場合は `sv_pressure_constraints` と `linear_reduction` を使います。

```cpp
vcp::bfem::sv_pressure_constraints<2, T> cons(msh, n);
vcp::bfem::linear_reduction<T> pred(ps.ndof(), cons.rows());

auto B1 = pred.reduce_rows(B0);
```

この部品群は鞍点系の solve までは行いません。行列を作った後のブロック結合、定数圧力モードの処理、
線形ソルバ選択は利用者側で行います。

## 8. 2D C^1 要素を使う

この節は「領域全体」の話です。`c1_space` は C^1 自由度を大域番号で管理し、内部で物理要素上の局所行列を作って散布します。
`eval_grad` と `eval_hess` は大域関数を指定要素へ制限して評価します。

C^1 要素は 2D 専用で、次数は `k >= 5` です。

```cpp
const int c1k = 5;
vcp::bfem::c1_space<2, T> cs(msh, c1k);

auto A = cs.stiffness(c1k);
auto M = cs.mixed_mass(c1k, c1k);

auto u = cs.zero_function(c1k);
auto L = cs.laplacian_matrix(c1k);
auto H = cs.hessian_matrix(c1k);
T r2 = cs.laplacian_residual_sq(f, u);
```

勾配や Hessian の点評価もできます。

```cpp
vcp::bfem::bary_point<2, T> lam = {{T(1) / T(3), T(1) / T(3), T(1) / T(3)}};
std::array<T, 2> g = cs.eval_grad(u, 0, lam);
std::array<T, 3> h = cs.eval_hess(u, 0, lam); // xx, xy, yy
```

C^1 の境界条件は一般線形拘束として生成し、`linear_reduction` で消去します。
この API は境界の方向分類と拘束係数生成のために `c1_coord_traits<T>` を要求します。
`detail::rational` には既定特殊化がありますが、`double` や独自スカラー型では利用者側で
「符号判定」と「座標差の有理数化」を定義してください。例えば整数座標だけを使う
簡単な double スモークなら、次のように書けます。

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

この例は座標差が整数である場合の説明用です。一般の double 座標では、丸め済みの実数を
どの有理数として扱うかを明示的に決める必要があります。

```cpp
std::vector<vcp::bfem::sv_constraint_q> rows =
    vcp::bfem::c1_boundary(cs, c1k, vcp::bfem::c1_bc_dirichlet);

vcp::bfem::linear_reduction<T> red(cs.ndof(c1k), rows);
auto Ar = red.reduce(A);
```

## 9. 区間型で使うとき

この節はスカラー型と変換の話なので、主に「領域非依存」です。ただし幾何符号判定や疎行列 policy は、
物理要素上・領域全体の API を呼ぶ時点で効いてきます。

区間型では traits を明示 include します。

```cpp
#include <kv/interval.hpp>
#include <vcp/bfem/bound_traits_kv.hpp>
#include <vcp/bfem/geometry_traits_kv.hpp>
#include <vcp/bfem/convert_traits_kv.hpp>

typedef kv::interval<double> IT;
```

多項式係数は `from_rational` を使うと、厳密な有理係数を区間へ一度だけ変換できます。

```cpp
std::vector<std::pair<long long, long long> > a;
a.push_back(std::make_pair(1, 1));
a.push_back(std::make_pair(0, 1));
a.push_back(std::make_pair(1, 3));
vcp::bfem::poly1<IT> f = vcp::bfem::poly1<IT>::from_rational(a);
```

区間型で疎行列を作る場合は、`vcp::spmatrix<T,SP>` の `SP` が区間型に対応している必要があります。
対応していない場合は、ベクトル・スカラー・局所多項式 API だけを使うか、区間対応 sparse policy を渡します。

## 10. OpenMP による並列化

この節は「領域全体」の行列組立の話です。bfem では、要素ごとの局所行列を作って
大域行列へ散布する部分の一部が OpenMP で並列化されています。

OpenMP を有効にするには、通常は利用者プログラムを `-fopenmp` 付きでコンパイルします。
`_OPENMP` が定義され、かつ `VCP_BFEM_NOMP` が定義されていなければ、bfem 側の
OpenMP 経路が有効になります。

```bash
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 -fopenmp \
sandbox/tests/example.cpp \
-L${MKLROOT}/lib/intel64 \
-Wl,--no-as-needed \
-lmkl_intel_lp64 \
-lmkl_intel_thread \
-lmkl_core \
-liomp5 \
-lpthread \
-lm \
-ldl \
-lmpfr \
-o sandbox/bin/example
```

スレッド数を利用者プログラムから指定したい場合だけ、利用者側で `<omp.h>` を include します。
単に bfem 内部の OpenMP を有効にするだけなら、利用者コードに `<omp.h>` は不要です。

```cpp
#include <omp.h>

int main() {
    omp_set_num_threads(8);
    // bfem の行列組立を呼ぶ
}
```

環境変数で指定することもできます。

```bash
OMP_NUM_THREADS=8 ./sandbox/bin/example
```

bfem の OpenMP だけを無効にしたい場合は、コンパイル時に `-DVCP_BFEM_NOMP` を付けます。
プロジェクト全体の方針として `-DVCP_NOMP` を付けた場合も、bfem では
`VCP_BFEM_NOMP` として扱われます。

```bash
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 -fopenmp \
-DVCP_BFEM_NOMP \
sandbox/tests/example.cpp \
-L${MKLROOT}/lib/intel64 \
-Wl,--no-as-needed \
-lmkl_intel_lp64 \
-lmkl_intel_thread \
-lmkl_core \
-liomp5 \
-lpthread \
-lm \
-ldl \
-lmpfr \
-o sandbox/bin/example_nomp
```

主に並列化されるのは、次のような大域行列の組立です。

| 空間 | 並列化される代表 API |
|---|---|
| 2D/3D `P^k` | `stiffness`, `mixed_mass`, `weighted_mass` |
| 2D/3D `RT^k` | `rt_space::mass`, `assemble_div_mass`, `assemble_cross_grad` |
| broken `P_l` | `broken_space::mass` |
| Scott-Vogelius | `assemble_vector_stiffness`, `assemble_div_velocity`, `assemble_advection`, `assemble_advection_derivative` |
| 2D C^1 | `stiffness`, `mixed_mass`, `weighted_mass`, `laplacian_matrix`, `hessian_matrix`, `laplacian_mixed` |

`load` ベクトルや一部のスカラー評価は、現時点では逐次の API もあります。
また、1 つの `fe_space` や `c1_space` インスタンスを複数の利用者 thread から同時に呼ぶ使い方は想定していません。
bfem 内部では、1 回の組立呼び出しの中で要素ループを並列化します。

`double` では浮動小数点加算順序の問題がありますが、bfem の OpenMP 組立では
thread-local COO buffer を使い、結合順序を固定する方針です。検証付き計算で
`kv::interval<double>` などを使う場合も、区間対応 sparse policy を使う点は逐次実行時と同じです。

## 11. コンパイル例

プロジェクトの通常方針に従い、Ubuntu/WSL では g++、OpenMP、MKL、MPFR を使います。
`sandbox/tests/example.cpp` を作った場合の例です。

```bash
mkdir -p sandbox/bin

g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
sandbox/tests/example.cpp \
-L${MKLROOT}/lib/intel64 \
-Wl,--no-as-needed \
-lmkl_intel_lp64 \
-lmkl_intel_thread \
-lmkl_core \
-liomp5 \
-lpthread \
-lm \
-ldl \
-lmpfr \
-fopenmp \
-o sandbox/bin/example

./sandbox/bin/example
```

macOS ではプロジェクト方針どおり `clang++ -std=gnu++14`、OpenBLAS、MPFR、GMP、libomp を使います。

## 12. よくある失敗

| 症状 | 原因と対処 |
|---|---|
| `degenerate_element` | 要素が潰れている、または区間座標で向きが確定しない。メッシュや座標区間を確認する |
| 係数サイズ不一致 | `function_from_coeffs` へ `ndof(m) x 1` でない行列を渡している |
| 区間型で疎行列がコンパイルできない | sparse policy `SP` が区間型に対応していない |
| OpenMP が効いていない | `-fopenmp` が付いていない、または `VCP_BFEM_NOMP` / `VCP_NOMP` を定義している |
| スレッド数を変えられない | `OMP_NUM_THREADS` を設定するか、利用者コードで `<omp.h>` を include して `omp_set_num_threads` を呼ぶ |
| `linear_reduction` が失敗する | 拘束行が従属、ピボットがない、または自由度番号が範囲外 |
| RT 補間後に期待した誤差が出ない | provider が法線連続でない可能性がある |
| C^1 で `m < k` が拒否される | `c1_space` は基準次数以上の family degree だけを扱う |
