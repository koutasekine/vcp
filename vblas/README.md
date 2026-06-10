# vblas

`vblas/` は、丸め方向を指定して double 行列積を計算するための実験的な
BLAS 風カーネル群です。外部 BLAS がハードウェア丸めモードの変更を
精度保証付き数値計算に必要な形で反映しない場合に、代替 backend として
使えるかを検討する目的で置かれています。

通常の VCP Library の安定 API ではありません。`vcp::matrix`、
`vcp::pdblas`、`vcp::pidblas` から自動的に使われるものではなく、
利用側が `vblas/rmatmul.hpp` などを明示的に include して呼び出します。

## 主なファイル

- `rmatmul.hpp`
  - 丸めモードを指定して 1 回の行列積 `A * B` を計算する入口です。
  - AVX-512、AVX2 + FMA、AArch64 NEON、no-SIMD fallback の各実装を
    compile 時の feature macro により選択します。
- `rmatmul_avx512.hpp`
  - AVX-512 実装です。
  - `_mm512_fmadd_round_pd` などの embedded rounding 付き intrinsic を
    使います。
- `rmatmul_avx2.hpp`
  - AVX2 + FMA 実装です。
  - FMA 命令が現在の floating-point environment の丸めモードに従うため、
    `std::fesetround()` により計算中の丸め方向を切り替えます。
- `rmatmul_neon.hpp`
  - AArch64 NEON 実装です。
  - AVX2 実装と同様に `std::fesetround()` で計算中の丸め方向を
    切り替えます。
- `rmatmul_nosimd.hpp`
  - 対応する SIMD backend が有効でない場合に使う scalar 実装です。
  - AVX2 実装と同様に block 化し、OpenMP の各 thread で丸め方向を
    一時的に切り替えます。
  - AVX2 / AVX-512 の FMA backend と同じ丸め付き積和の性質を保つため、
    scalar 計算でも `std::fma()` を使います。そのため、環境によっては
    software 実装や関数呼び出しのコストが支配的になり、SIMD backend より
    大幅に遅くなります。
- `udmatmul.hpp`
  - AVX-512 専用の古い実験用実装です。
  - 上向き丸め結果 `CU` と下向き丸め結果 `CD` を 1 回の呼び出しで
    同時に計算します。
- `Check_rmatmul_rounding.cpp`、`Check_udmatmul_rounding.cpp`
  - 丸め方向の違いが結果に反映されるかを確認するための検証コードです。
- `main.cpp`
  - 大きめの行列で `rmatmul` と `vcp::pdblas` の結果を比較する
    実験用 driver です。

## `rmatmul`

`rmatmul` の呼び出し形式は次の通りです。

```cpp
#include "rmatmul.hpp"

rmatmul(m, n, k, A, B, C, rounding_mode);
```

行列は column-major です。

- `A`: `m * n` 行列
- `B`: `n * k` 行列
- `C`: `m * k` 行列
- `C` は `A * B` の結果で上書きされます。

`rounding_mode` の意味は次の通りです。

- `1`: 上向き丸め
- `-1`: 下向き丸め
- `0`: 最近点丸め
- 上記以外の値: 最近点丸めとして扱います。

上向き丸めと下向き丸めの両方が必要な場合は、出力先を分けて
2 回呼び出します。

```cpp
rmatmul(m, n, k, A.data(), B.data(), CU.data(),  1);
rmatmul(m, n, k, A.data(), B.data(), CD.data(), -1);
```

## backend 選択

`rmatmul.hpp` は compile 時に利用可能な SIMD 実装を選びます。
優先順位は次の通りです。

1. `__AVX512F__` が定義されている場合: AVX-512
2. `__AVX2__` と `__FMA__` が定義されている場合: AVX2 + FMA
3. AArch64 NEON が有効な場合: NEON
4. 上記のいずれも有効でない場合: no-SIMD fallback

AVX-512 対応 CPU でも AVX2 実装を確認したい場合は、`-mavx512f` を付けずに
`-mavx2 -mfma` を付けてコンパイルしてください。

no-SIMD fallback は移植性のための実装であり、高速化を目的とした backend では
ありません。特に no-SIMD fallback では丸め付き積和を維持するために
`std::fma()` を使うため、AVX2 / AVX-512 の hardware FMA を使う実装より
遅くなることがあります。性能が必要な場合は、可能な限り対象 CPU に合わせて
AVX2 + FMA、AVX-512、または NEON を有効にしてコンパイルしてください。

## 丸めモードの取り扱い

`rmatmul` は、呼び出し元から指定された `rounding_mode` で行列積を計算します。
一方で、呼び出し元が `rmatmul` 実行前に設定していた floating-point
environment の丸めモードは、`rmatmul` の終了後も維持されるようにしています。

実装ごとの扱いは次の通りです。

- AVX-512 実装
  - embedded rounding 付き AVX-512 intrinsic により、命令ごとに丸め方向を
    指定します。
  - global な `cfenv` の丸めモードは変更しません。
- AVX2 実装
  - 計算前に `std::fegetround()` で現在の丸めモードを保存します。
  - OpenMP の各 worker thread では、計算用に `std::fesetround()` で
    `FE_UPWARD`、`FE_DOWNWARD`、または `FE_TONEAREST` を設定します。
  - 計算終了後、保存していた丸めモードへ `std::fesetround()` で戻します。
- NEON 実装
  - AVX2 実装と同じ方針で、`std::fegetround()` による保存と
    `std::fesetround()` による復元を行います。
- no-SIMD 実装
  - AVX2 実装と同じ方針で、OpenMP の各 worker thread ごとに丸めモードを
    保存、変更、復元します。

つまり、AVX2 / NEON / no-SIMD では計算中に各スレッドの丸めモードを
一時的に変更しますが、`rmatmul` から戻った後の呼び出し元スレッドの
丸めモードは、呼び出し前の状態に戻ります。OpenMP worker thread についても、
各 parallel region に入った時点の丸めモードへ戻します。

## `udmatmul`

`udmatmul` は AVX-512 の embedded rounding を使い、上向き丸め結果と
下向き丸め結果を同時に計算する実験用関数です。

```cpp
#include "udmatmul.hpp"

udmatmul(m, n, k, A, B, CU, CD);
```

- `CU`: 上向き丸めの結果
- `CD`: 下向き丸めの結果
- `CU` と `CD` は計算結果で上書きされます。
- AVX-512 非対応環境では利用できません。

現在は、より汎用的な入口として `rmatmul` を使用することを想定しています。

## コンパイル例

プロジェクトルートから実行する場合の AVX2 実装の確認例です。

```bash
mkdir -p sandbox/bin

g++ -I. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
-mavx2 -mfma \
vblas/Check_rmatmul_rounding.cpp \
-lpthread \
-lm \
-ldl \
-fopenmp \
-o sandbox/bin/check_rmatmul_rounding

./sandbox/bin/check_rmatmul_rounding
```

AVX-512 実装を使う場合は、AVX-512 用の compile option を追加します。

```bash
g++ -I. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
-mavx512f -DVCP_USE_AVX512 \
vblas/Check_rmatmul_rounding.cpp \
-lpthread \
-lm \
-ldl \
-fopenmp \
-o sandbox/bin/check_rmatmul_rounding_avx512
```

MKL や OpenBLAS、MPFR を含む既存の VCP 検証コードと一緒にリンクする場合は、
プロジェクトルートの `AGENTS.md` にある通常の compile option に従ってください。

## 注意

- 対象は点 double 行列です。一般の区間行列 backend ではありません。
- まだ実験段階のため、通常の VCP Library 利用では `vcp/` と `docs/` 以下の
  文書を参照してください。
- 将来 VCP の公開 API に昇格する場合は、CPU feature check、丸め方向付き計算の
  正しさテスト、性能評価、policy としてのインターフェース設計が必要です。
