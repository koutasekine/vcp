# btfem — 放物型方程式の全離散 Petrov–Galerkin 法

`vcp/btfem` は、放物型方程式の全離散 Petrov–Galerkin 法の行列・ベクトル・
スカラー量を作るヘッダ群です。方法は次の論文に従います。

- K. Hashimoto, T. Kinoshita, M. T. Nakao, "Numerical verification of solutions for
  nonlinear parabolic problems", arXiv:1811.04261v2

空間方向は bfem の連続 Bernstein `P^n` 空間(`fe_space`)をそのまま使い
([bfem.md](bfem.md))、時間方向は 1 次元の Bernstein 多項式で張った時間空間を
組み合わせます。bfem と同じく、要素積分に数値積分は使いません。時間方向の積分の
テーブルは有理数のまま作り、要素型 `T` へ一度だけ変換します。`T` に区間型を使えば
行列・ベクトル・スカラー量は厳密値の包含として、有理数型を使えば厳密な値として
計算されます。

名前空間は `vcp::btfem`(誤差定数は `vcp::btfem::constants`)です。空間側の型
(`mesh`、`fe_space`、`fe_function`、`poly_field`、`poly1`、`bpoly`、`bary_point`、
`range_pair`、`dirichlet_reduction`、`rt_space`、`broken_space` など)は `vcp::bfem` の
ものを使います。

btfem は行列・ベクトル・内積・ノルムの 2 乗を作るところまでを担います。線形方程式の
求解、Newton 反復([newton.md](newton.md) の `SpNewton`)、固有値・作用素ノルムの計算、
流束 σ の決定、初期値 `u_0` から左端の係数を決める補間・射影は、利用者側で行います。
時刻の断面は `fe_function` で返るので、図示には bfem の `output_uh_for_graphics` を
使えます。

## 必要なヘッダ

| 用途 | ヘッダ |
| --- | --- |
| 時間メッシュ | `<vcp/btfem/time_mesh.hpp>` |
| 時間空間 | `<vcp/btfem/time_space.hpp>` |
| 空間時間関数 | `<vcp/btfem/st_poly.hpp>` |
| Kronecker 因子と縮約 | `<vcp/btfem/kron.hpp>` |
| 全離散 Petrov–Galerkin 空間 | `<vcp/btfem/pg_space.hpp>` |
| RT ⊗ 時間、broken ⊗ 時間 | `<vcp/btfem/st_rt_space.hpp>`、`<vcp/btfem/st_broken_space.hpp>` |
| 流束系の行列・ノルム | `<vcp/btfem/st_rt_assemble.hpp>` |
| 誤差定数 | `<vcp/btfem/constants/time_projection_constants.hpp>`、`<vcp/btfem/constants/inverse_inequality.hpp>`、`<vcp/btfem/constants/heat_pg_constants.hpp>` |

`kron.hpp` は `pg_space.hpp` と独立に使えます。

## 計算の考え方

全離散の試行空間は `V_k^q ⊗ S_h^n`、基底は `psi_a(t) phi_i(x)`、試験関数の基底は
`∂_t psi_a ⊗ phi_i` です。自由度の順序は時間ブロックを主とする
`(a, i) ↦ a · n_x + i` に固定されていて、各時間ブロックはそのまま bfem の
`fe_function` の係数ベクトルになります。

大域行列は、時間の行列と空間の行列の組 `kron_term` を基本として扱います。
`mn x mn` の `spmatrix` は `kron_expand` で作り、行列を作らずに積だけを計算する場合は
`kron_apply` を使います。時間に依存する係数の項 `(c phi_j, ∂_t phi_i)` も `kron_term` の
列で返るので、その和を `kron_expand` に渡します。方程式の係数(拡散係数 ν など)は
API に埋め込まず、利用者が `kron_term::scale` に入れます。

計算は次の順に進みます。

1. `time_mesh<T>` が `[t0, t1]` を `m` 等分し、小区間の幅 `k` と `1/k` を持ちます
   (除算はここで 1 回だけ)。
2. `time_space<T>` が時間の自由度番号、大域の時間行列、小区間ごとの局所テーブルを
   作ります。
3. `pg_space<D, T, P, SP>` が bfem の `fe_space` を参照し、小区間 × 要素のループで
   行列の因子・ベクトル・スカラー量を作ります。
4. `kron_expand` / `kron_apply` が因子から `spmatrix` または行列ベクトル積を作ります。

加算の順序は「小区間 → 要素 → 時間の局所番号 → 空間の局所番号」で固定されています。
同じ入力・1 スレッドなら、ビット単位で同じ結果になります。

## スカラー型と行列 policy

要素型 `T` の条件は bfem と同じです。`vcp::bfem::rational`、`double`、`kv::dd`、
`kv::mpfr<N>`、`kv::interval<double>`、`kv::interval<kv::dd>`、
`kv::interval<kv::mpfr<N>>` で使えます。区間型・kv の点の型では、bfem と同じく
次も include します。

```cpp
#include <vcp/bfem/bound_traits_kv.hpp>
#include <vcp/bfem/geometry_traits_kv.hpp>
#include <vcp/bfem/convert_traits_kv.hpp>
```

`range` / `range_composed` を kv の点の型(`kv::dd`、`kv::mpfr<N>`)や有理数型で
使う場合は、`bound_traits<T>` を `point_bound_traits<T>` から派生させる特殊化を
利用者側で与えます(bfem の `refine.hpp` と同じ)。

除算はプログラム全体で「空間の要素数(bfem の幾何量)+ 時間メッシュの数
(`time_mesh::uniform`)」回だけです。

空間の因子と展開した行列は `vcp::spmatrix<T, SP>` で返ります。`SP` の既定は
`vcp::spmats<T>` で、区間型では bfem と同じく区間に対応した疎行列の policy
(`vcp::spimats`)を `pg_space<D, T, P, SP>` に渡します。btfem は `spmatrix` を
戻り値とデータの保持にだけ使い、`mul_vec` などの演算は使いません(`kron_apply` は
成分を列挙して計算します)。返る `spmatrix` は `finalize()` していません。

## 時間メッシュと時間空間

```cpp
#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/time_space.hpp>

vcp::btfem::time_mesh<double> tm = vcp::btfem::time_mesh<double>::uniform(0.0, 1.0, 4);
vcp::btfem::time_space<double> Vt(tm, 2, vcp::btfem::time_bc_zero_left);

vcp::matrix<double, vcp::mats2<double> > Ut = Vt.mass<vcp::mats2<double> >();   // (psi_b, psi_a)
```

`time_mesh::uniform(t0, t1, m)` は `[t0, t1]` を `m` 等分します。次の時間区間は、
利用者が `uniform(t1, t2, m')` で別に作ります。等間隔でない分割はありません。

時間空間の境界条件 `time_bc` は次の 2 つです。

| `time_bc` | 空間 | 自由度数 |
| --- | --- | --- |
| `time_bc_none` | `Vbar_k^q`(左端に初期値を持つ) | `m q + 1` |
| `time_bc_zero_left` | `V_k^q`(左端の係数を除く) | `m q` |

`V_k^q` では `dt_mass + dt_mass^T == y y^T`、`Vbar_k^q` では
`dt_mass + dt_mass^T == y y^T - e_0 e_0^T` が成り立ちます(`y` は `end_values`)。

## 全離散 Petrov–Galerkin 空間

```cpp
#include <vcp/btfem/pg_space.hpp>

vcp::bfem::fe_space<2, double> Vh(msh, 1);
vcp::btfem::pg_space<2, double> S(Vh, msh, tm);

auto A = S.A(1, 2);                                     // A_t(q = 1) ⊗ U_x(n = 2)
vcp::spmatrix<double> K = vcp::btfem::kron_expand(A);
```

`pg_space` は bfem の `fe_space` を所有しないので、`fe_space` は `pg_space` より長く
生存させてください。時間メッシュはコピーして持ちます。試行・試験の次数 `(q, n)` は
毎回明示します(既定値はありません)。空間の因子は境界条件を入れていない全自由度の
系で、Dirichlet 条件の縮約は `kron.hpp` の `block_reduce`、または因子の段階で
`kron_term::space` に `R.reduce(...)` を適用して行います。両者は厳密に一致します。

`Vbar_k^q` と `V_k^q` の添字の相互変換(左端のブロックの除去・挿入)には
`initial_reduction` を使います。

### `SpNewton` で使う場合の手順

1. `fe_space`、`mesh`、`time_mesh` から `pg_space` を作ります。
2. 未知数は `V_k^{qbar} ⊗ S_h^{nbar}` の係数(左端は固定)です。
   `initial_reduction::expand` で `Vbar` の `ubar` を復元し、`load_dt` /
   `convection_dt` と `kron_apply(A_mixed / B_mixed)` で残差を作ります。
3. ヤコビ行列は、`A`、`B`(`scale = ν`)、`weighted_dt_mass_terms`、
   `convection_dt_terms` を `kron_expand` した正方行列です。空間の Dirichlet 条件は
   `block_reduce` または因子の段階の `R.reduce` で入れます。
4. 時間発展は、`section_at_node(ubar, m)` を次の区間の初期値にして、
   `uniform(t1, t2, m')` で続けます。

物理時刻 `t` から小区間の番号と局所座標 `(l, tau)` を求める処理はありません
(bfem の `eval(u, e, lam)` と同じ方針)。

## RT ⊗ 時間と broken ⊗ 時間

放物型の残差を `L^2(J; H^{-1})` で扱うための流束系の量を作ります。bfem の
`rt_space` / `broken_space` / `rt_assemble.hpp` を空間側にそのまま使い、時間方向を
`Vbar_k^q` / `V_k^q` で張ります。どの不等式で組み合わせるかは利用者が決めます。

```cpp
#include <vcp/btfem/st_rt_space.hpp>
#include <vcp/btfem/st_broken_space.hpp>
#include <vcp/btfem/st_rt_assemble.hpp>

vcp::bfem::rt_space<2, double> rs(msh, 1);
vcp::btfem::st_rt_space<2, double> srs(rs, tm);
auto sig = srs.zero_field(2, vcp::btfem::time_bc_none);   // Vbar_k^2 ⊗ RT^1
auto M   = srs.mass(2, vcp::btfem::time_bc_none);         // U_t ⊗ M_RT

auto ub = S.zero_function(1, vcp::btfem::time_bc_none, 2);
double e = vcp::btfem::st_flux_error_sq(srs, sig, S, ub);  // ||grad ub - sig||^2
```

時間空間 `(q, bc)` は、空間を作るときではなく、場や行列を作るたびに指定します。
`rt_space` / `broken_space` は所有しないので、`st_rt_space` / `st_broken_space` より
長く生存させてください。

残差 `w = g(ū) − ∂_t ū − b·∇ū` は、小区間 × 要素ごとに局所的に作り、大域には
作りません。`g` は `poly1<T>`、`b` は定数(`std::array<T, D>`)または空間時間の場
(`std::array<st_poly, D>`、全成分が同じ `(q_b, bc)`)です。`b` を渡さない版は、
bfem の `f(u)` の形に対応します。

RT の自由度は、辺・面の向きによる符号を持ちます(bfem と同じ)。時間方向は
テンソル積なので、符号は空間の自由度ごとに決まり、すべての時間の自由度に同じ符号が
掛かります。`to_vector` / `field_from_vector` は符号を変えません。

`st_flux_error_sq` は bfem の `flux_error_sq` と同じく `‖∇ū − σ‖^2` を返し、
係数 ν を持ちません。流束 `σ ≈ ν∇ū` を使う場合は、`σ/ν` を渡すか結果に `ν^2` を
掛けてください。

時間に依存しない双一次形式 `b_x(·, ·)` の時間版は、
`kron_term{ T(1), Vt_row.mixed_mass<P>(Vt_col), B_x }` の 1 行で作れます
(`B_x` は bfem の行列)。

### bfem との対応

| bfem(空間) | btfem(空間時間) | 戻り値 |
| --- | --- | --- |
| `rt_space`、`rt_field` | `st_rt_space`、`st_rt_field` | — |
| `broken_space`、`broken_field` | `st_broken_space`、`st_broken_field` | — |
| `rt_space::mass()` | `st_rt_space::mass(q, bc)` | `kron_term` |
| `broken_space::mass()` / `inner(u, v)` | `st_broken_space::mass(q, bc)` / `inner(u, v)` | `kron_term` / `T` |
| (bfem にはない) | `st_div_field(srs, σ, sbs)` | `st_broken_field` |
| `assemble_div_mass(bs, rs)` | `st_assemble_div_mass(sbs, qb, bcb, srs, qs, bcs)` | `kron_term` |
| `assemble_cross_grad(rs, fs, m)` | `st_assemble_cross_grad(srs, qs, bcs, pg, q, n)` | `kron_term` |
| `assemble_mixed_mass(fs, m, bs)` | `st_assemble_mixed_mass(pg, q, n, sbs, qb, bcb)` | `kron_term` |
| `flux_error_sq(rs, σ, fs, u)` | `st_flux_error_sq(srs, σ, pg, ū)` | `T` |
| `div_residual_sq(rs, σ, f, fs, u)` | `st_div_residual_sq(srs, σ, g[, b], pg, ū)` | `T` |
| `broken_load(bs, f, fs, u)` | `st_broken_load(sbs, qb, bcb, g[, b], pg, ū)` | `matrix` |
| `projection_error_sq(bs, f, fs, u)` | `st_projection_error_sq(sbs, g[, b], pg, ū)` | `T` |

bfem の `f(u)` に対して、btfem では `w = g(ū) − ∂_t ū − b·∇ū` を使う点だけが違います。
ū と σ が時間によらなければ、各量は bfem の対応する量の `|J|` 倍に一致します。

## 誤差定数

`vcp/btfem/constants/` は、線形熱方程式 `∂_t u − νΔu = f`(`u(a) = 0`)の全離散
Petrov–Galerkin 近似の 3 つの誤差定数を、2 つの射影誤差定数 C_Ω(h)、C_τ と離散量
λ_h から組み立てます。式と記号は次の論文に従います。

- M. Mizuguchi, M. T. Nakao, K. Hashimoto, K. Sekine, S. Oishi, "Error estimates for
  the fully-discrete Petrov-Galerkin method of a linear parabolic problem",
  Japan J. Indust. Appl. Math. (2025)

名前空間は `vcp::btfem::constants` で、要素型 `T` は kv の区間型
(`kv::interval<double>`、`kv::interval<kv::dd>`、`kv::interval<kv::mpfr<N>>`)に
限ります。`double` や有理数型ではコンパイルエラーになります。

```cpp
#include <vcp/btfem/constants/heat_pg_constants.hpp>

typedef kv::interval<double> IT;
// Th は vcp::bfem::mesh<2, IT>

auto sm = vcp::btfem::constants::heat_pg_space_materials_h01<2, IT>(Th, 1);   // C_Omega^2, lambda_h
auto c  = vcp::btfem::constants::heat_pg_error_constants_h01<IT>(sm, IT(1) / IT(4), 1, IT(1));
// c.K, c.L, c.c1_tilde, c.c0_small, c.c0_tilde
```

初期誤差を 0(`u_0 = û_0 = 0`)、`k = q`、`C_τ := C_τ(q − 1)` として、論文の
定理 1〜3 から次の量を作ります。

| 量 | 式 | 不等式 |
| --- | --- | --- |
| K = C̃₁² | (C_Ω² + ν² λ_h C_τ² + 4ν C_τ) / ν² | ‖u − u_{h,τ}^k‖²_{L²H¹₀} ≤ K ‖f‖² |
| L = c̃₀² | 4 (C_Ω² / ν + C_τ) | ‖(u − u_{h,τ}^k)(b)‖²_{L²} ≤ L ‖f‖² |
| C̃₀ | 4 C_τ + 2 C_Ω C̃₁ | ‖u − u_{h,τ}^k‖_{L²L²} ≤ C̃₀ ‖f‖ |

| 関数 | ヘッダ | 戻り値 |
| --- | --- | --- |
| `time_projection_error_constant<T>(d, tau)` | `time_projection_constants.hpp` | C_τ(d) = c_d · τ |
| `inverse_inequality_h01_bound<D, T, DP, SP>(Th, n)` | `inverse_inequality.hpp` | λ_h の上界(点区間。下端は下界ではない) |
| `heat_pg_space_materials_h01<D, T, DP, SP>(Th, n[, rep])` | `heat_pg_constants.hpp` | 空間側の材料 {`n`、`c_omega_sq`、`lambda_h`} |
| `heat_pg_error_constants_h01<T>(c_omega_sq, c_tau, lambda_h, nu)` | 同上 | 3 定数(式だけを計算する版。材料を差し替えられる) |
| `heat_pg_error_constants_h01<T>(sm, tau, q, nu)` | 同上 | 同上(`c_tau = time_projection_error_constant(q − 1, tau)`) |
| `heat_pg_error_constants_h01<T>(sm, tm, q, nu)` | 同上 | 同上(`tau = tm.width()`) |
| `heat_pg_error_constants_h01<D, T, DP, SP>(Th, n, tau, q, nu[, rep])` | 同上 | 同上(材料を内部で作る一括版。材料版とビット単位で一致) |

空間方向の C_Ω(h)² には、bfem の `ritz_projection_error_constant_h01_sq(Th, n)` を
そのまま使います。`_h01` は斉次 Dirichlet 条件の空間を仮定する印で、C_τ は境界条件に
よらないので付きません。2 乗の量 K、L が基本で、2 乗でない版は平方根を取ります。

C_τ(d) = c_d · τ の c_d の出典は次の通りです。

| d = q − 1 | c_d | 出典 | 種別 |
| --- | --- | --- | --- |
| 0 | 1/π | M. T. Nakao, N. Yamamoto, S. Kimura, J. Approx. Theory 93 (1998) の 1 章(Schultz, Spline Analysis の定理 2.5 (2.17)) | 最良値 |
| 1 | 1/(2π) | 同論文の定理 1 | 最良値 |
| 2 | (1/(3√5π³))·√(−2160 + 25π⁴ + 4√5·√(−174960 − 1080π⁴ + 243π⁶ + 5π⁸)) ≈ 1/8.92337 | 同論文の定理 2 (14) | 上界 |
| ≥ 3 | √max(L1, L2)、L1 = 1/(2(2d+1)(2d+5)) + 1/(4(2d+5)√(2d+3)√(2d+7))、L2 = 1/(4(2d+5)√(2d+3)√(2d+7)) + 1/(2(2d+5)(2d+9)) + 1/(4(2d+9)√(2d+7)√(2d+11)) | S. Kimura, N. Yamamoto, Bull. Inform. Cybernet. 31 (1999) の定理 2.3 (2.3)(N = d) | 上界 |

C_τ は閉じた式の区間評価で、`kv::interval<double>` と `kv::interval<kv::mpfr<N>>` では
真の値の包含です。λ_h と bfem の κ_h² は上端だけに意味がある点区間なので、それらを
材料に含む K、L、C̃₁、c̃₀、C̃₀ も上端(上界)として読みます。

内部の区間行列計算に使う `DP` / `SP` の既定は、bfem の誤差定数と同じ
`vcp::imats<typename T::base_type>` / `vcp::spimats<typename T::base_type>` です。

## OpenMP

`pg_space` の時間に依存する係数の項(`weighted_dt_mass_terms`、空間時間の `b` の
`convection_dt_terms`)の要素のループは OpenMP で並列化されています。無効にする方法は
bfem と同じ(`-DVCP_BFEM_NOMP` / `-DVCP_NOMP`)です。スレッドごとの作業領域は
`pg_space` が持ち、同じ形の 2 回目以降の呼び出しでは要素のループの中でメモリを
確保しません。スレッド数が違うときのビット一致は保証しませんが、1 スレッドの OpenMP
ビルドと OpenMP なしのビルドは同じ順序で加算するので一致します。小区間のループと、
ベクトル・スカラーの計算は逐次です。

## 例外

btfem が投げる例外は、すべて `vcp::error` から派生します
([error_handling.md](error_handling.md))。

| 例外 | 主な条件 |
| --- | --- |
| `vcp::btfem::degenerate_interval`(`vcp::verification_error` の派生) | `time_mesh::uniform` の区間の長さが 0 以下、または区間型で符号を確定できない |
| `vcp::invalid_argument` | `m < 1`、`q < 1`、`n < 1`、`node(l)` の範囲外、サイズ・次数の条件違反、別の `pg_space` / `st_rt_space` / `st_broken_space` の場、`A_mixed` / `B_mixed` で `nbar != n`、`kron_expand` の項の形の不一致、`initial_reduction` の `full` / `reduced` の不整合、要素数の不一致、空間・時間メッシュの不一致 |
| `vcp::invalid_argument`(誤差定数) | `d < 0`、`tau` の下端 ≤ 0、縮約後の自由度が 0、`nu` の下端 ≤ 0、`c_tau` / `lambda_h` / `c_omega_sq` の下端 < 0 |
| `vcp::verification_error`(誤差定数) | 固有値の包含の証明に失敗、上端の最大値が負など |
| bfem の例外 | 空間側の退化・不整合はそのまま伝わる |

## 関数リファレンス

### `time_mesh<T>`(`<vcp/btfem/time_mesh.hpp>`)

| 関数 | 内容 |
| --- | --- |
| `static time_mesh uniform(t0, t1, m)` | `[t0, t1]` を `m >= 1` 等分。`k = (t1 - t0) * (1/m)`、`1/k = m / (t1 - t0)` |
| `t_begin()`、`t_end()` | `t0`、`t1` |
| `num_intervals()` | `m` |
| `width()` | 小区間の幅 `k` |
| `inv_width()` | `1/k`(`k` とは独立に計算した値。区間型では `k * inv_k` は 1 を含むだけ) |
| `node(l)` | `tau_l = t0 + l k`(`l = 0..m`) |
| `same_as(other)` | 同じメッシュか(`m` と両端。区間型では重なりで判定) |

### `time_space<T>`(`<vcp/btfem/time_space.hpp>`)

局所番号 `r = 0..q` は、`r = 0` が左端、`r = q` が右端です。大域番号は
`g(l, r) = (l-1) q + r`(`time_bc_none`)、`(l-1) q + r - 1`(`time_bc_zero_left`、
`(1, 0)` は `-1`)です。

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `time_space(tm, q, bc)` | | `q >= 1` |
| `degree()`、`bc()`、`num_intervals()`、`local_size()`、`ndof()` | | `q`、`bc`、`m`、`q + 1`、自由度数 |
| `global_dof(l, r)` | `int` | 大域番号。除いた自由度は `-1` |
| `mesh_ref()` | `const time_mesh<T>&` | 時間メッシュ |
| `mass<P>()` | `matrix<T, P>` | `U_t = (psi_b, psi_a)` |
| `stiffness<P>()` | `matrix<T, P>` | `A_t = (psi_b', psi_a')` |
| `dt_mass<P>()` | `matrix<T, P>` | `B_t = (psi_b, psi_a')`。行 `a` が試験(微分)の側 |
| `end_values<P>()` | `matrix<T, P>` | `y_a = psi_a(t_end)` |
| `mixed_mass<P>(other)`、`mixed_stiffness<P>(other)`、`mixed_dt_mass<P>(other)` | `matrix<T, P>` | 行 = この空間、列 = `other`(同じ時間メッシュに限る) |
| `local_mass()`、`local_stiffness()`、`local_dt_mass()` | | 1 小区間の局所テーブル(`k · Mref`、`(1/k) · Sref`、`Bref`) |
| `local_dt_load(sdeg)`、`local_weighted_dt_mass(sdeg, sigma)`、`local_weighted_mass(sdeg, sigma)` | | 時間に依存する係数の局所テーブル |
| `same_mesh(other)` | `bool` | 同じ時間メッシュか |

### `st_poly<D, T, P>`(`<vcp/btfem/st_poly.hpp>`)

`u(t, x) = Σ_a Σ_i C[a][i] psi_a(t) phi_i(x)` の係数行列 `C`(行 = 時間の自由度、
列 = 空間の自由度)を持ちます。行優先に並べたものが時間ブロック主のベクトル
`v[a n_x + i]` です。

| 関数 | 内容 |
| --- | --- |
| `time_degree()`、`time_bc_kind()`、`space_degree()` | `q`、`bc`、`n` |
| `num_time_dofs()`、`num_space_dofs()` | 時間・空間の自由度数 |
| `coeffs()` | 係数行列(直接編集できる) |

`st_poly` は `pg_space` のファクトリでだけ作れます。作った `pg_space` を記憶し、別の
`pg_space` に渡すと例外になります。

### Kronecker 因子と縮約(`<vcp/btfem/kron.hpp>`)

```cpp
struct kron_term {
    T                      scale;   // 既定 T(1)。ν などはここに入れる
    vcp::matrix<T, P>      time;    // 時間の行列(密)
    vcp::spmatrix<T, SP>   space;   // 空間の行列
};
```

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `kron_expand(term)` | `spmatrix` | 時間ブロック主で展開。行 `a n_x + i`、列 `b n_x' + j`、値 `scale · time(a, b) · space(i, j)` |
| `kron_expand(terms)` | `spmatrix` | 項の和 |
| `kron_apply(terms, x)` | `matrix<T, P>` | 行列を作らずに `y = Σ scale (time ⊗ space) x` |
| `block_reduce(R, nblocks, v)`、`block_expand(R, nblocks, vr)` | `matrix<T, P>` | 時間ブロックごとに `dirichlet_reduction R` を適用、その逆 |
| `block_reduce(R, nblocks, A)` | `spmatrix` | 展開した行列の行・列をブロックごとに縮約 |

`initial_reduction<T, P>` は `Vbar_k^q` と `V_k^q` の添字を相互に変換します。

| 関数 | 内容 |
| --- | --- |
| `initial_reduction(full, reduced, n_x)` | `full` は `time_bc_none`、`reduced` は同じ次数・同じメッシュの `time_bc_zero_left` |
| `full_size()`、`reduced_size()`、`space_size()` | `(m q + 1) n_x`、`m q n_x`、`n_x` |
| `reduce(v_full)` | 先頭の `n_x` 成分(時間の自由度 0 のブロック)を除く |
| `expand<D>(v_reduced, u_init)` | 先頭のブロックに `u_init.coeffs()` を置いて `Vbar` の添字に戻す |

### `pg_space<D, T, P, SP>`(`<vcp/btfem/pg_space.hpp>`)

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `pg_space(Vh, msh, tm)` | | 要素数の一致を検査する |
| `space()`、`mesh_ref()`、`time_mesh_ref()` | | 参照先の bfem 空間、メッシュ、時間メッシュ |
| `time_dofs(q, bc)` | `const time_space<T>&` | `(q, bc)` ごとに保持する時間空間 |
| `zero_function(q, bc, n)`、`function_from_coeffs(q, bc, n, C)` | `st_poly` | 零関数、係数行列からの関数 |
| `function_from_vector(q, bc, n, v)`、`to_vector(u)` | | 時間ブロック主のベクトルとの相互変換 |
| `A(q, n)`、`M(q, n)`、`B(q, n)`、`U(q, n)`、`W(q, n)`、`Y(q, n)` | `kron_term` | `A_t⊗U_x`、`U_t⊗K_x`、`B_t⊗K_x`、`U_t⊗U_x`、`A_t⊗K_x`、`y yᵀ⊗U_x`(試行・試験とも `V_k^q⊗S_h^n`) |
| `A_mixed(q, n, ubar)`、`B_mixed(q, n, ubar)` | `kron_term` | 行 = `V_k^q`、列 = `ubar` の時間空間。`nbar == n` が必要 |
| `weighted_dt_mass_terms(fprime, ubar, q, n)` | `std::vector<kron_term>` | `(f'[ubar] phi_j, ∂_t phi_i)` の因子の列 |
| `convection_dt_terms(b, q, n)` | `std::vector<kron_term>` | `(b·∇phi_j, ∂_t phi_i)`。`b` は定数または `std::array<st_poly, D>` |
| `load_dt(g, ubar, q, n)` | `matrix<T, P>` | `(g(ubar), ∂_t phi_i)` |
| `convection_dt(b, ubar, q, n)` | `matrix<T, P>` | `(b·∇ubar, ∂_t phi_i)` |
| `section_at_node(u, l)` | `fe_function` | 時刻 `tau_l` の断面(`time_bc_zero_left` で `l = 0` は零関数) |
| `section_at(u, l, tau)` | `fe_function` | 小区間 `l` の局所重心座標 `tau = (1-t, t)` での断面。`tau` が幅を持てば値の包含 |
| `time_trace(u, e, lam)` | `std::vector<bpoly<1, T> >` | 空間の点 `(e, lam)` での小区間ごとの時間多項式 |
| `eval(u, l, tau, e, lam)` | `T` | 点評価 |
| `range(u)`、`range_composed(f, ubar)` | `range_pair<T>` | 値域の包含(`f(ubar)` の版は局所の係数から) |
| `elevate_space(u, n2)` | `st_poly` | 各時間ブロックを空間の次数 `n2` に上げる |
| `initial_error_sq_l2(u0, ubar)`、`initial_error_sq_h01(u0, ubar)` | `T` | `‖u_0 − ubar(t_0)‖^2`、`‖∇(u_0 − ubar(t_0))‖^2`。`u0` は `fe_function` または `poly_field` |

線形部分 `(∂_t ubar, ∂_t phi_i) + ν(∇ubar, ∇∂_t phi_i)` は、`ubar` が `V_k^q` なら
`kron_apply({A(q, n)}, to_vector(ubar))` と `kron_apply({B(q, n)}, ...)`、`ubar` が
`time_bc_none` なら `A_mixed` / `B_mixed` で得られます。

### `st_rt_space<D, T, P, SP>`(`<vcp/btfem/st_rt_space.hpp>`)

`st_rt_field` は係数行列(行 = 時間の自由度、列 = RT の大域自由度。符号は bfem の規約)を
持ちます。

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `st_rt_space(rs, tm)` | | |
| `space()`、`time_mesh_ref()`、`time_dofs(q, bc)` | | bfem の RT 空間、時間メッシュ、時間空間 |
| `zero_field(q, bc)`、`field_from_coeffs(q, bc, C)` | `st_rt_field` | 零の場、係数行列からの場 |
| `field_from_vector(q, bc, v)`、`to_vector(s)` | | 時間ブロック主のベクトル `(a, j) ↦ a · ndof_rt + j` との相互変換 |
| `section_at_node(s, l)`、`section_at(s, l, tau)` | `rt_field` | 時刻の断面 |
| `mass(q, bc)` | `kron_term` | `U_t(q) ⊗ M_RT` |
| `validate(s)` | `void` | 場の所属・次数・サイズの検査 |

`st_rt_field` は `order()`、`time_degree()`、`time_bc_kind()`、`num_time_dofs()`、
`num_space_dofs()`、`coeffs()` を持ちます。

### `st_broken_space<D, T, P, SP>`(`<vcp/btfem/st_broken_space.hpp>`)

`st_rt_space` と同じ形で、列は bfem の broken の要素ブロック順 `e · nloc + r`、
符号は `+1` です。

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `st_broken_space(bs, tm)` | | |
| `space()`、`time_mesh_ref()`、`time_dofs(q, bc)` | | `st_rt_space` と同じ |
| `zero_field`、`field_from_coeffs`、`field_from_vector`、`to_vector` | | `st_rt_space` と同じ |
| `section_at_node(u, l)`、`section_at(u, l, tau)` | `broken_field` | 時刻の断面 |
| `mass(q, bc)` | `kron_term` | `U_t(q) ⊗ M_broken`(ブロック対角) |
| `inner(u, v)` | `T` | `(u, v)_{L^2(J; L^2)}`。`u`、`v` の `(q, bc)` は異なってよい |
| `validate(u)` | `void` | 場の検査 |

### `st_rt_assemble.hpp` の自由関数(`<vcp/btfem/st_rt_assemble.hpp>`)

`pg` の試行空間は `V_k^q`(`time_bc_zero_left`)で、`(q, n)` を取ります。`ū` は
`pg_space` の `st_poly`(`time_bc_none` / `time_bc_zero_left` のどちらでもよい)です。

| 関数 | 戻り値 | 内容 |
| --- | --- | --- |
| `st_div_field(srs, σ, sbs)` | `st_broken_field` | `div σ`(σ と同じ `(q, bc)`)。broken の次数と RT の次数が等しいことが必要 |
| `st_assemble_div_mass(sbs, qb, bcb, srs, qs, bcs)` | `kron_term` | 行 = broken ⊗ `(qb, bcb)`、列 = RT ⊗ `(qs, bcs)`。`(div σ_j, q_i)` |
| `st_assemble_cross_grad(srs, qs, bcs, pg, q, n)` | `kron_term` | 行 = RT ⊗ `(qs, bcs)`、列 = `V_k^q ⊗ S_h^n`。`(σ_i, ∇φ_j)` |
| `st_assemble_mixed_mass(pg, q, n, sbs, qb, bcb)` | `kron_term` | 行 = `V_k^q ⊗ S_h^n`、列 = broken ⊗ `(qb, bcb)`。`(φ_i, q_j)` |
| `st_flux_error_sq(srs, σ, pg, ū)` | `T` | `‖∇ū − σ‖^2_{L^2(J; L^2)}`(`n̄ >= 1`) |
| `st_div_residual_sq(srs, σ, g[, b], pg, ū)` | `T` | `‖div σ + w‖^2_{L^2(J; L^2)}` |
| `st_broken_load(sbs, qb, bcb, g[, b], pg, ū)` | `matrix<T, P>` | `(w, q_i)_{L^2(J; L^2)}` |
| `st_projection_error_sq(sbs, g[, b], pg, ū)` | `T` | `‖(I − Π_l) w‖^2_{L^2(J; L^2)}`。`Π_l` は空間方向の要素ごとの `L^2` 射影(時間方向は射影しない) |
