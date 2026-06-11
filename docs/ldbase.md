# Legendre 基底（ldbase）

`vcp/ldbase.hpp` は、Legendre 多項式を変形して構成した H^1_0 の基底
（Legendre 基底）を用いて、単位立方体領域 (0,1)^d 上の連立半線形楕円型
偏微分方程式

```text
-Delta u = f(u)   in (0,1)^d
       u = 0      on boundary
```

の解の精度保証付き数値計算を行うためのヘッダです。近似解の計算
（Newton 法用の行列生成）から、誤差評価に必要な各種行列・積分値・
定数の精度保証付き計算までを 1 つのクラス
`vcp::Legendre_Bases_Generator` が担当します。

基底の構成は

> M. T. Nakao and T. Kinoshita, On very accurate verification of solutions
> for boundary value problems by using spectral methods,
> JSIAM Letters, vol. 1 (2009), pp. 21--24.
> https://doi.org/10.14495/jsiaml.1.21

に基づきます。また、本ヘッダの使い方の概要は

> NOLTA, IEICE, vol. 12, no. 1 (2021), pp. 41--74, Section 7.
> https://doi.org/10.1587/nolta.12.41

に記載されています。実際の利用例は `test_PDE/` 以下
（`test_Emden.cpp`、`test_Gray_Scott.cpp` など）にあります。

## Legendre 基底の定義

区間 (0,1) 上の（シフトされた）Legendre 多項式 `P_n` に対して、
1 次元の基底関数を

```text
phi_n(x) = x (1 - x) / (n (n - 1)) * dP_{n-1}(x)/dx ,   n >= 2
```

と定義します。`phi_n(0) = phi_n(1) = 0` を満たすので H^1_0(0,1) の
基底になります。d 次元（d 変数）の基底はテンソル積

```text
Phi_{i1,...,id}(x1,...,xd) = phi_{i1}(x1) * ... * phi_{id}(xd)
```

で構成します。クラス内部では多重指数 `(i1,...,id)` を 1 列に並べた
「リスト」（`output_list()` で取得可能）で基底を管理し、近似解は

```text
uh(x) = sum_i  uh_i  Phi_i(x)
```

の係数ベクトル（`vcp::matrix`、変数が複数なら列が変数）で表します。

内部の数値積分には Gauss-Legendre 積分を使い、積分点・重みや基底の
点値は高精度型（後述の `_T`）で生成してから計算用の型に変換します。
積分次数はモードと多項式次数から自動で決まり、被積分関数が多項式の
範囲では打ち切り誤差なしで（区間演算なら包含として）計算できます。

## Include

```cpp
#include <kv/interval.hpp>   // 事前に interval.hpp が必要
#include <kv/rdouble.hpp>
#include <kv/mpfr.hpp>       // kv の MPFR 型
#include <kv/rmpfr.hpp>      // MPFR の丸め方向制御（rop 版）
#include <vcp/matrix.hpp>
#include <vcp/ldbase.hpp>
```

`ldbase.hpp` は `interval.hpp` を include 済みであることを要求します
（未 include の場合はコンパイルエラーになります）。

後述の `_T` には実質的に `kv::interval< kv::mpfr<N> >` などの MPFR
ベースの高精度区間型が必須です。基底（高次の Legendre 多項式）や
Gauss-Legendre 積分点の生成では係数が極端に大きく・小さくなり、
`double` 程度の精度・指数範囲ではオーバーフローや精度劣化が起こる
ためです。そのため `kv/mpfr.hpp` と、区間演算で丸め方向を制御する
`kv/rmpfr.hpp` の両方を include してください。リンク時には MPFR
ライブラリ（`-lmpfr`）が必要です。

## クラステンプレート

```cpp
vcp::Legendre_Bases_Generator< _T, _TM, _PM >
```

| 引数 | 意味 |
| --- | --- |
| `_T` | 基底・積分点の生成に使う高精度型。オーバーフロー防止のため MPFR ベースの区間型がほぼ必須（例: `kv::interval< kv::mpfr<500> >`） |
| `_TM` | 行列・積分値の計算に使う型（例: `double`, `kv::dd`, `kv::interval<double>`） |
| `_PM` | `vcp::matrix` の計算 policy（省略時 `vcp::mats<_TM>`） |

典型的な typedef は次の通りです。

```cpp
typedef kv::interval< kv::mpfr<500> > DataType; // 基底生成用（高精度）

// 近似計算用
typedef kv::dd AppData;
typedef vcp::mats< AppData > POLICY;

// 精度保証用（区間演算）
typedef kv::interval< double > VData;
typedef vcp::pidblas VPOLICY;

vcp::Legendre_Bases_Generator< DataType, AppData, POLICY > Approximate_Generator;
vcp::Legendre_Bases_Generator< DataType, VData, VPOLICY > Verification_Generator;
```

精度保証で使う場合は `_TM` を区間型にしてください。丸め誤差まで含めた
包含が得られます。

## 計算の流れ

1. `setting()` でモード・次数などを設定する
2. `setting_list()` または `setting_evenlist()` で基底リストを生成する
3. `setting_uh()` で近似解の係数を設定する（必要な場合）
4. 行列生成関数・積分関数・定数計算関数を呼ぶ
5. 使い終わったら `clear()`

## setting

```cpp
void setting(int m, int p, int dim, int vs = 1, int k = 2, int uh_o = -1);
```

| 引数 | 意味 |
| --- | --- |
| `m` | Legendre 基底の次数（1 次元あたりの基底数、`m >= 2`）。多項式次数は `m + 2` |
| `p` | 非線形項の最大べき（`-Delta u = u^p` の `p` のイメージ、`p >= 1`） |
| `dim` | 領域の次元 `d`（変数の個数、`dim >= 1`） |
| `vs` | 連立する未知関数の個数（例: `u, v` の 2 成分系なら `2`） |
| `k` | モード（下表） |
| `uh_o` | `uh`（または `cx`）の基底次数。基底次数 `m` と異なる次数の `uh` を与えるときに指定（`k == 1` または `k == 0` のときのみ有効） |

モード `k` は生成する内部データと使える関数を決めます。

| `k` | モード | 主な用途 |
| --- | --- | --- |
| `1` | Inverse mode | 逆作用素ノルム評価用。`(uh^p phi, phi)` 型の行列生成 |
| `2` | Residual mode | 残差ノルム評価用。`(Laplace uh, Laplace uh)` などの積分値 |
| `0` | Goerisch mode | Goerisch の定理による固有値評価用 |
| `k > 2` | Approximation mode | 近似解の計算。Gauss-Legendre 積分次数は `10k` |

近似モードでは `k` を大きくするほど積分点が増えます（例: `k = 5` で
積分次数 50）。検証モード（`k = 0, 1, 2`）では、必要な積分次数が
`m`, `p`, `uh_o` から自動で決まります。

## 基底リストの生成

```cpp
void setting_list();      // すべての基底を使う
void setting_evenlist();  // 偶数番目（x = 1/2 で対称）の基底のみを使う
```

`setting()` の直後に必ずどちらかを呼びます。解が領域の中心に関して
対称な場合は `setting_evenlist()` を使うと、1 次元あたりの基底数が
約半分になり次元数の累乗で効きます。リストは

```cpp
vcp::matrix< int > list_uh = Generator.output_list();
```

で取得でき、行数 `list_uh.rowsize()` が全基底数（= 係数ベクトル `uh`
の行数）になります。

## 近似解 uh の設定

```cpp
void setting_uh();                              // uh をすべて 1 で初期化
void setting_uh(const vcp::matrix<_TM,_PM>& uh);  // 係数を直接与える
void setting_uh(const vcp::matrix<_TM,_PM>& uh,
                const vcp::matrix<int>& list_uh, int ind);
```

- 1, 2 番目は `setting()` で `uh_o` を指定していない場合（`uh` の次数 =
  基底の次数）に使います。`uh` のサイズは（基底数, 変数の個数）です。
- 3 番目は `uh` の次数が基底の次数と異なる場合（`setting()` の `uh_o`
  指定時）に使います。`list_uh` は `uh` を計算したときの Generator の
  `output_list()` の値、`ind` はそのリストの種類で
  full list なら `1`、even list なら `2` を渡します。

典型的には、近似モードの Generator で求めた `uh` を、検証モードの
Generator に渡します。

```cpp
// 近似（基底次数 uh_m）で uh を計算した後、
// 検証用 Generator（基底次数 m, モード 1）に uh を設定する
Verification_Generator.setting(m, p, dim, vs, 1, uh_m);
Verification_Generator.setting_evenlist();
vcp::matrix< VData, VPOLICY > uhi;
vcp::convert(uh, uhi);
Verification_Generator.setting_uh(uhi, list_uh, 2); // even list なので 2
```

Goerisch mode（`k = 0`）では、同様の役割の関数

```cpp
void setting_cx(const vcp::matrix<_TM,_PM>& cx,
                const vcp::matrix<int>& list_cx, int ind);
```

で係数関数 `cx` を設定します（`setting()` の `uh_o` に `cx` の次数を
指定しておきます）。

## 行列生成関数（モード 1 と近似モード）

Newton 法や固有値問題で使う行列・ベクトルを生成します。
Residual mode（`k = 2`）と Goerisch mode（`k = 0`）では使えません。

| 関数 | 返り値 |
| --- | --- |
| `dphidphi()` | 行列 `( (nabla Phi_i, nabla Phi_j)_{L^2} )_{i,j}` |
| `phiphi()` | 行列 `( (Phi_i, Phi_j)_{L^2} )_{i,j}` |
| `uhphiphi(p1, ..., pvs)` | 行列 `( (u1^p1 ... uvs^pvs Phi_i, Phi_j)_{L^2} )_{i,j}` |
| `uhphi(p1, ..., pvs)` | ベクトル `( (u1^p1 ... uvs^pvs, Phi_j)_{L^2} )_j` |
| `fphi()` | ベクトル `( (1, Phi_j)_{L^2} )_j` |

`uhphiphi` と `uhphi` の引数は変数ごとのべき指数で、引数の個数は
`setting()` の `vs` と一致させます。べきの合計は `p` 以下が必要です。
`dphidphi` と `phiphi` は Legendre 基底の直交性に基づく厳密な公式で
計算されるため、数値積分誤差はありません。

例（2 変数系、Gray-Scott の場合）:

```cpp
DL = Generator.dphidphi();          // (nabla phi, nabla phi)
L  = Generator.phiphi();            // (phi, phi)
uh2vhphi  = Generator.uhphi(2, 1);    // ベクトル (u^2 v, phi)
uh2phiphi = Generator.uhphiphi(2, 0); // 行列 (u^2 phi, phi)
uhvhphiphi = Generator.uhphiphi(1, 1);// 行列 (u v phi, phi)
```

例（1 変数 Emden 方程式 `-Delta u = u^2` の Newton 法）:

```cpp
while (true) {
    Generator.setting_uh(uh);
    uh2phi   = Generator.uhphi(2);     // (uh^2, phi)
    uhphiphi = Generator.uhphiphi(1);  // (uh phi, phi)
    DF = DL - 2 * uhphiphi;
    F  = DL * uh - uh2phi;
    uh = uh - lss(DF, F);              // 連立一次方程式を解いて更新
    // 収束判定 ...
}
```

## 積分値の計算（モード 2 と 0）

残差ノルム `|| Laplace uh + f(uh) ||_{L^2}` の評価などに使う積分値を
返します。`L = Laplace` と書きます。

| 関数 | 返り値 | 使用可能モード |
| --- | --- | --- |
| `integral_uh(p1, ..., pvs)` | `int u1^p1 ... uvs^pvs dx` | 2 |
| `integral_LuhLuh(v)` | `int (L u_v)^2 dx` | 2 |
| `integral_LuhLuh(v1, v2)` | `int (L u_{v1}) (L u_{v2}) dx` | 0, 2 |
| `integral_Luhuh(v, p1, ..., pvs)` | `int u1^p1 ... uvs^pvs (L u_v) dx` | 2 |
| `integral_Luhcxuh(vL, vu, p1, ..., pcx)` | `int cx1^p1 ... (L u_{vL}) u_{vu} dx` | 0 |
| `integral_cxuhuh(v1, v2, p1, ..., pcx)` | `int cx1^p1 ... u_{v1} u_{v2} dx` | 0 |

変数番号 `v` は 0 始まりです。例（Emden 方程式の残差）:

```cpp
Generator.setting(uh_m, p, dim, vs, 2);  // Residual mode
Generator.setting_evenlist();
Generator.setting_uh(uhi);               // uhi は区間行列

VResData uh4    = Generator.integral_uh(4);       // int uh^4 dx
VResData LuhLuh = Generator.integral_LuhLuh(0);   // int (L uh)^2 dx
VResData Luhuh2 = Generator.integral_Luhuh(0, 2); // int uh^2 (L uh) dx

// || L uh + uh^2 ||_{L^2}^2 = LuhLuh + 2 Luhuh2 + uh4
Res = sqrt(abs(LuhLuh + 2 * Luhuh2 + uh4));
```

## 定数の精度保証付き計算

誤差評価に現れる各種定数を返します。`_TC` には区間型を指定します
（領域は (0,1)^d 固定です）。

| 関数 | 評価式 |
| --- | --- |
| `Ritz_projection_error<_TC>(N = -1)` | `\|\| nabla(u - Rh u) \|\|_{L^2} <= CN \|\| L u \|\|_{L^2}` の `CN` |
| `weighted_Ritz_projection_error<_TC>(w, N = -1)` | `sqrt(\|\| nabla(u - Rh u) \|\|^2 + \|\| u \|\|^2) <= CNw \|\| L u + w u \|\|_{L^2}` の `CNw` |
| `Poincare_constant<_TC>()` | `\|\| u \|\|_{L^2} <= Cp \|\| nabla u \|\|_{L^2}` の `Cp` |
| `weighted_Poincare_constant<_TC>(w)` | `\|\| u \|\|_{L^2} <= Cpw sqrt(\|\| nabla u \|\|^2 + w \|\| u \|\|^2)` の `Cpw` |
| `Sobolev_constant<_TC>(q)` | `\|\| u \|\|_{L^q} <= Cs \|\| nabla u \|\|_{L^2}` の `Cs` |
| `weighted_Sobolev_constant<_TC>(q, w)` | `\|\| u \|\|_{L^q} <= Csw sqrt(\|\| nabla u \|\|^2 + w \|\| u \|\|^2)` の `Csw` |

`Ritz_projection_error` の `N` を省略すると `setting()` の基底次数 `m`
が使われます。`Sobolev_constant` の `q` は埋め込みが成り立つ範囲
（`d = 1` なら `q > 0`、`d = 2` なら `q > 2`、`d >= 3` なら
`q <= 2d/(d-2)` など）で指定してください。

```cpp
Cp  = Generator.Poincare_constant< VData >();
Cs3 = Generator.Sobolev_constant< VData >(3);
CN  = Generator.Ritz_projection_error< VData >();
```

## 近似解の評価・出力

| 関数 | 意味 |
| --- | --- |
| `func<_TC>(x)` | 点 `x`（`std::vector`、サイズ `dim`）での `(x, u1(x), ..., uvs(x))` を行列で返す |
| `func<_TC>(x, v)` | 変数 `v` の点 `x` での値 |
| `global_min(x, mesh)` | 区間ベクトル `x` 上での各変数の下限値（精度保証付き、モード 2 以外） |
| `global_max(x, mesh)` | 区間ベクトル `x` 上での各変数の上限値（精度保証付き、モード 2 以外） |
| `output_uh_for_graphics(Div)` | 各軸 `Div` 分割の格子点での `(x, uh(x))` を行列で返す（近似モードのみ） |
| `output_list()` | 基底の多重指数リスト |
| `clear()` | 内部データの消去 |

`global_min` / `global_max` は区間の分割（branch and bound）により
`uh` の最小値・最大値を包含します。`mesh` は分割を打ち切る最小幅
（省略時 `2^{-11}`）です。`uh` が対称な場合は対称性を使って探索領域を
狭められます。

```cpp
std::vector< kv::interval<double> > x(Dimension);
for (int d = 0; d < Dimension; d++) x[d] = kv::interval<double>(0, 0.5);

std::vector<double> uh_min = Generator.global_min(x, std::pow(2.0, -9));
std::vector<double> uh_max = Generator.global_max(x, std::pow(2.0, -9));
```

`output_uh_for_graphics` の返り値は各行が
`(x1, ..., xd, u1, ..., uvs)` で、gnuplot 等での可視化に使えます。

## マクロ

| マクロ | 効果 |
| --- | --- |
| `VCP_LEGENDRE_DEBUG`（または `VCP_DEBUG`） | 生成過程の詳細ログを標準出力に表示 |
| `VCP_LEGENDRE_NOMP`（または `VCP_NOMP`） | OpenMP 並列化を無効化 |

行列生成・積分点評価・`global_min` などは OpenMP で並列化されて
います。コンパイル時に `-fopenmp` を付けると有効になります。

## コンパイル例

Ubuntu (WSL) + MKL の場合:

```bash
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 test_PDE/test_Emden.cpp \
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
-fopenmp
```

## 関連ファイル

| ファイル | 役割 |
| --- | --- |
| `vcp/ldbase.hpp` | 本体（`Legendre_Bases_Generator`） |
| `vcp/ldbase_assist.hpp` | 補助ヘッダ（`ldbase.hpp` の後に include） |
| `vcp/legendre_integral.hpp` | Gauss-Legendre 積分点・重みの精度保証付き生成 |
| `vcp/matrix.hpp` | 行列クラス `vcp::matrix`（[matrix.md](matrix.md)） |
| `test_PDE/test_Emden.cpp` | 1 変数（Emden 方程式）の検証例 |
| `test_PDE/test_Gray_Scott.cpp` | 2 変数連立系（Gray-Scott）の検証例 |
| `test_PDE/test_Laplace_eigenvalue.cpp` | 固有値評価の例 |
