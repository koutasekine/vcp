# Fourier Series

`vcp/fourier_series.hpp` は、1 変数の実 Fourier 級数を係数列として扱うための
ヘッダです。PDE や周期解の計算で、非線形項を Fourier 係数上で組み立てる
ときに使います。

基本形は次の通りです。

```text
u(t) = a0 + sum_{m=1}^n (am sin(mt) + bm cos(mt))
```

`vcp::fourier_series<T>` は、この `a0` と、各次数 `m` の `am`, `bm` を
保持します。内部では `vcp::sincos<T>` が `am sin(mt) + bm cos(mt)` の
係数対を表します。

`fourier_series` 自体は定数項を `a0` として持ちます。一方で
`fourier_basis.hpp` や `fourier_cos_basis.hpp` の行列ベクトル表現では、
基底の先頭を `1/2` として扱います。この変換の係数調整は
`vec_to_fourier_series` や `fphi` が行います。

## Include

```cpp
#include <vcp/fourier_series.hpp>
```

`T` には `double`、`kv::interval<double>`、`kv::dd`、`kv::mpfr<N>` など、
四則演算、`sin`、`cos`、`sqrt` などが使える型を想定します。

## 初期化と係数設定

`zeros(n)` は、0 次から `n` 次までを扱う Fourier 級数を 0 で初期化します。
`n` は最高次数です。`set_sinm(a, m)` と `set_cosm(b, m)` は、`sin(mt)` と
`cos(mt)` の係数を設定します。

```cpp
vcp::fourier_series<double> u;

u.zeros(3);
u.set_a0(1.0);
u.set_sinm(2.0, 1);  // 2 sin(t)
u.set_cosm(3.0, 2);  // 3 cos(2t)

std::cout << u << std::endl;
```

この例は、おおよそ次の級数を表します。

```text
u(t) = 1 + 2 sin(t) + 3 cos(2t)
```

`set_sinm` と `set_cosm` は、あらかじめ `zeros(n)` などで必要な次数まで
確保してから使ってください。例えば `set_cosm(a, 5)` を使う場合は、
`zeros(5)` 以上の初期化が必要です。

係数を読むときは次を使います。

| 関数 | 意味 |
| --- | --- |
| `get_a0()` | 定数項 `a0` |
| `get_sin(m)` | `sin(mt)` の係数 |
| `get_cos(m)` | `cos(mt)` の係数 |
| `get_n()` | 保持している最高次数 |

`get_sin(m)` と `get_cos(m)` は、保持次数を超える `m` に対しては 0 を返します。
次数 `m` は 1 以上として使ってください。

## 演算

`fourier_series` は、Fourier 級数同士の加減算、乗算、スカラー倍を
演算子で扱えます。

```cpp
vcp::fourier_series<double> u, v;

u.zeros(2);
u.set_a0(1.0);
u.set_cosm(1.0, 1);

v.zeros(2);
v.set_sinm(2.0, 1);

auto f = u * v + 3.0 * u - 1.0;
```

乗算では、三角関数の積和公式により必要な次数まで係数を展開します。
例えば、`n` 次と `m` 次の級数を掛けると、結果は最大で `n + m` 次まで
持ちます。非線形項を Fourier 係数上で構成するときに重要な機能です。

## よく使う関数

| 関数 | 意味 |
| --- | --- |
| `delay(tau)` | `u(t + tau)` に対応する係数へ変換 |
| `mul_sin(N)` | `u(t) sin(Nt)` を Fourier 級数として返す |
| `mul_cos(N)` | `u(t) cos(Nt)` を Fourier 級数として返す |
| `diff()` | `du/dt` を Fourier 級数として返す |
| `Pm(m)` | `m` 次までの射影 |
| `I_minus_Pm(m)` | `m` 次以下を 0 にした高周波側 |
| `value(t)` | 点 `t` での値 |
| `L2norm()` | `(1/pi int_0^{2pi} |u(t)|^2 dt)^{1/2}` |
| `Linfnorm(N)` | 格子評価と 2 階微分係数による `L^\infty` 評価 |
| `clear()` | 係数を消去 |

`Pm(m)` と `I_minus_Pm(m)` は、現在保持している次数以下の `m` で使ってください。

## Fourier 基底クラスとの関係

`fourier_series.hpp` は Fourier 級数の係数演算を担当します。一方で、
Newton 法や PDE の離散化で必要になる行列・ベクトルへの変換は、次のヘッダが
担当します。

| ヘッダ | 役割 |
| --- | --- |
| `vcp/fourier_basis.hpp` | `1/2, sin(t), cos(t), ...` を使う一般 Fourier 基底 |
| `vcp/fourier_cos_basis.hpp` | Neumann 条件などで使う cosine 基底 |
| `vcp/fourier_quadrature.hpp` | 関数を数値積分して Fourier 係数や Jacobi 行列を作る |

例えば `fourier_cos_basis` では、行列ベクトルから Fourier 級数を作り、
非線形項を `fourier_series` の演算で構成してから、再びベクトルへ戻します。

```cpp
vcp::fourier_cos_basis<double, vcp::pdblas> basis;
basis.setting_m(100);

vcp::matrix<double, vcp::pdblas> z;
z.zeros(101, 1);
z(0) = 2.0;  // 1/2 基底の係数。fourier_series では a0 = 1 になる
z(1) = 1.0;  // cos(t) の係数

vcp::fourier_series<double> u = basis.vec_to_fourier_series(z);

auto f = u * u * u - u;
vcp::matrix<double, vcp::pdblas> coeff = basis.fphi(f);
```

実際の利用例としては、`test_PDE/test_1d_Neumann_Gray_Scott.cpp` が
`fourier_series` と `fourier_cos_basis` を使っています。

## OpenMP

OpenMP が有効な場合、`fourier_series` 内部の一部ループは並列化されます。
利用者側で外側のループを並列化するなど、Fourier 関連の内部 OpenMP を止めたい
場合は、次を定義してください。

```bash
-DVCP_FOURIER_NOMP
```

全体の内部 OpenMP を止める場合は、次も使えます。

```bash
-DVCP_NOMP
```
