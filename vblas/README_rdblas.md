# rdblas — 丸めモード指定付き double BLAS

`rdblas` は，LAPACK を構築する前提で必要となる double 精度 BLAS の全 routine を，
丸めモードを指定して計算できる形で提供する header-only library です．
`rmatmul` と同様に実験的な位置付けで，`vcp::matrix` からは自動では使われません．
この上に構築した丸めモード指定付き LAPACK (rdlapack) は兄弟 directory の
`vlapack/` にあります (`vlapack/README_rdlapack.md` 参照)．

```cpp
#include "vblas/rdblas.hpp"

// C := alpha*A*B + beta*C を上向き丸めで計算
vcp::rdgemm('N', 'N', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, 1);
```

あるいはファイル先頭で `using namespace vcp;` を宣言すれば，
`vcp::` 修飾なしで呼び出せます:

```cpp
#include "vblas/rdblas.hpp"
using namespace vcp;

rdgemm('N', 'N', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc, 1);
```

## 名前空間

**rdblas の全関数 (Level 1 / 2 / 3) は `vcp` 名前空間に収められています．**
`vcp::rdgemm`, `vcp::rdscal`, `vcp::dcopy` のように修飾するか，
呼び出し側ファイルで `using namespace vcp;` を宣言してください．

`dblas.hpp` の Fortran 互換 wrapper (`dgemm_` など) は `extern "C"` の
グローバルシンボルであり，名前空間に入っていません (変更なし)．

## 命名規約と引数

- FP 演算を含む routine は BLAS 名の先頭に `r` を付け，**reference BLAS と
  同順の引数の末尾に `rounding_mode` を 1 つ追加**します．
  - `rounding_mode = 1`: 上向き丸め (`FE_UPWARD`)
  - `rounding_mode = -1`: 下向き丸め (`FE_DOWNWARD`)
  - それ以外 (`0` を推奨): 最近点丸め (`FE_TONEAREST`)
- FP 演算を含まない routine は BLAS と同名のままです:
  `dcopy`, `dswap`, `idamax` (丸めモード引数なし)．いずれも `vcp` 名前空間内です．
- 行列は column-major，index は 0-based です．
  `idamax` の返り値も 0-based index です (Fortran BLAS の 1-based と異なる
  ので注意．`n <= 0` などでは -1 を返します)．
- 不正な引数 (`lda` 不足など) は `std::invalid_argument` を投げます
  (例外無効時は `abort`)．

## 提供 routine

| Level | routine |
|---|---|
| 1 | `rdscal` `rdaxpy` `rddot` `rdasum` `rdnrm2` `rdrot` `rdrotg` `rdrotm` `rdrotmg` / `dcopy` `dswap` `idamax` |
| 2 | `rdgemv` `rdgbmv` `rdsymv` `rdsbmv` `rdspmv` `rdtrmv` `rdtbmv` `rdtpmv` `rdtrsv` `rdtbsv` `rdtpsv` `rdger` `rdsyr` `rdspr` `rdsyr2` `rdspr2` |
| 3 | `rdgemm` `rdgemmtr` `rdsymm` `rdsyrk` `rdsyr2k` `rdtrmm` `rdtrsm` |

double 精度実数 LAPACK が参照する BLAS routine を全て含みます
(`rdgemmtr` は LAPACK 3.12.1 で reference BLAS に追加された GEMMTR
(C の uplo 三角部分のみ更新する GEMM) に対応するもので，
LAPACK 3.12.1 以降の dsytrf 系 (dlasyf) が必要とします)．

## dblas — 通常の BLAS interface (rdblas の wrapper)

`dblas.hpp` は rdblas (`vcp` 名前空間) を呼び出して構築した「通常の」double BLAS です．
関数名・引数は Fortran BLAS の C からの呼び出し規約に一致します
(末尾 underscore，全引数 pointer 渡し，INTEGER は 32bit `int` / LP64)．

```cpp
#include "vblas/dblas.hpp"

std::fesetround(FE_UPWARD); // 以後の dblas 呼び出しは上向き丸めで計算される
dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
std::fesetround(FE_TONEAREST);
```

- 各関数は呼び出し時点の丸めモードを取得し，
  `FE_DOWNWARD` なら `-1`，`FE_UPWARD` なら `1`，それ以外なら `0` を
  `vcp::rdblas` の末尾引数 `rounding_mode` として渡します．
  つまり**現在の丸めモードを尊重して計算する BLAS** になります．
- `KV_FASTROUND` 有効時は，`kv::hwround::roundup()` / `rounddown()` が直接変更する
  制御レジスタ (x86 では MXCSR，AArch64/ARM では FPCR/FPSCR) から丸め方向を読みます．
  `KV_FASTROUND` 無効時は通常の `std::fegetround()` を使います．
- 通常の BLAS と異なり，OpenMP の全 worker thread でも同じ丸めが保証され，
  計算後は各 thread の丸めモードが呼び出し前の状態へ復元されます．
- 丸めモードの変更を行わない関数 (`dcopy_` `dswap_` `idamax_`) は
  rdblas の同名関数をそのまま呼びます．`idamax_` の返り値だけは
  Fortran BLAS に合わせて 1-based index に直しています (`n <= 0` 等では 0)．
- LAPACK 3.12.1 で追加された `dgemmtr_` も提供します (全 35 関数)．
  reference LAPACK 3.12.1 以降と link する場合に dsytrf 系が必要とします．
- `lsame_` / `xerbla_` は提供しません．引数 error は rdblas に従い
  `std::invalid_argument` を投げます．
- Fortran object と link する library を作る場合は，1 つの翻訳単位で
  `VBLAS_DBLAS_EMIT_SYMBOLS` を定義して include すると全 symbol が発行されます
  (MKL 等と同名 symbol が衝突しないよう link 構成に注意)．

注意 (libiomp5 の FP 制御継承): Intel OpenMP runtime は既定
(`KMP_INHERIT_FP_CONTROL=true`) で，並列領域の入口に master の FP 制御状態を
各 worker へ継承させ，join 時に master の状態を入口時点へ戻します．
そのため `std::fesetround()` や `kv::hwround::roundup()` / `rounddown()` は
**並列領域の外で** (master thread に対して) 実行してください．
並列領域の内側で設定した丸めモードは join 後の master には残りません．

## ファイル構成

- `rdblas.hpp` — 入口 (これだけ include すればよい)
- `rdblas_common.hpp` — 丸めモード guard (RAII)・thread 数決定・
  alpha/beta epilogue・pack (copy / transpose copy) などの共通基盤
- `rdblas_level1.hpp` / `rdblas_level2.hpp` / `rdblas_level3.hpp` — 各 level の実装
- `dblas.hpp` — 通常の BLAS interface (Fortran 互換 wrapper，上記参照)

`rmatmul_common.hpp` は rdblas/dblas と同じ丸めモード保存・復元 helper を提供します．
`rdblas_level*.hpp` の公開関数は `namespace vcp { }` で囲まれています．

## 丸めモードの取り扱い (重要)

全 routine について次を保証します．

1. **計算前に丸めモードを保存**: 呼び出し thread と，OpenMP parallel region に
   入った各 worker thread のそれぞれが，自分の現在の丸めモードを保存します
   (`RoundingGuard` による RAII)．
   保存時の読み取りも dblas と同じく `KV_FASTROUND` 対応です．
2. **全 thread で丸めモードを変更**: 丸めモードは thread ごとの状態
   (x86 では MXCSR / x87 CW) なので，並列計算する全ての worker thread が
   parallel region の内側で `fesetround` を実行します．AVX-512 backend の
   rmatmul は embedded rounding (命令単位の丸め指定) を使うため
   global 状態を変更しませんが，防御的に同じ保存・復元を行います．
3. **計算後に全 thread の丸めモードを復元**: 各 thread が自分の保存値へ戻します．
   呼び出し側が事前に設定していた丸めモード (例えば `FE_DOWNWARD`) は
   呼び出し後もそのまま維持されます．

`sandbox/tests/test_rdblas.cpp` の `test_rounding_restore` と
「1/3 データで全要素 up != down」の test がこの 2 点・3 点を検証します
(1 thread でも切り替えに失敗すると担当要素が up == down になり検出されます)．

## 実装方針 (速度)

- **rdgemm**: op(A)/op(B) を連続 buffer に詰め (FP 演算なしの copy)，
  O(N^3) 本体は `rmatmul` (AVX-512 / AVX2+FMA / NEON / no-SIMD を compile 時に
  自動選択) で計算します．`alpha != 1 || beta != 0` の場合は temp に積を作り，
  `C := fma(beta, C, alpha*T)` の O(N^2) epilogue を指定丸めで適用します．
- **rdgemmtr**: 列 block (128 列) ごとに，三角の長方形部分は rdgemm で，
  対角 block は temp に積を作って三角部分のみ指定丸めの epilogue で更新します．
- **rdsymm**: 対称三角格納を full dense に展開 (copy のみ) して rdgemm に帰着．
- **rdsyrk / rdsyr2k**: `op(A)^T` を一度だけ連続 buffer に作り，行 band ごとに
  1 回の `rmatmul` + 三角 epilogue で更新します (band 方式)．
  band ごとの再 pack は左オペランド (band x k) のみで，余分な flops は
  band 対角 tile の半分 (全体の数 %) だけです．
- **rdtrmm / rdtrsm**: 再帰 2 分割で off-diagonal block を rdgemm に帰着します．
  trmm の leaf は三角 block を 0 詰め dense に展開して rmatmul で計算
  (0 との fma は結果を変えないため正当)，trsm の leaf は OpenMP 並列の
  scalar 代入 solve です．
- Level 1 / 2 は O(N^2) 以下で律速にならないため，`rdgemv` `rdger` `rdsymv`
  `rdsyr` `rdsyr2` のみ OpenMP 並列の単純 loop，それ以外は reference BLAS の
  loop 構造を移植した逐次実装です (auto-vectorize 可能な形)．

block size は compile 時 macro で調整できます (通常は変更不要):

```
-DVBLAS_RDBLAS_SYRK_BAND=256   (既定 0 = n に応じて自動)
-DVBLAS_RDBLAS_L3_LEAF=128     (trmm の leaf size)
-DVBLAS_RDBLAS_TRSM_LEAF=64    (trsm の leaf size)
```

### 性能の目安

AVX-512 (16 threads, n = 2000, 同一機で MKL と比較，値は環境依存):

| routine | rdblas (upward) | MKL (nearest) |
|---|---|---|
| dgemm | ~380 GFLOPS | ~254 GFLOPS |
| dsyrk | ~195 GFLOPS | ~222 GFLOPS |
| dtrsm | ~158 GFLOPS | ~304 GFLOPS |
| dtrmm | ~144 GFLOPS | ~276 GFLOPS |

rdgemm は丸めモード指定付きでも MKL dgemm より高速です．trmm/trsm 系は
温存 buffer 経由の更新が増えるため MKL の 5〜6 割程度ですが，
丸めモードの保証 (MKL にはない) を持ちます．

## 丸め方向の意味について

各 routine は「全ての浮動小数点演算を指定丸めで実行した結果」を返します．
積和だけからなる routine (dot, gemv, gemm, syrk など) では，同じ入力に対して
下向き丸めの結果 <= 上向き丸めの結果が要素ごとに成り立ちます．
一方，除算や負係数の混じる routine (trsv, trsm や alpha < 0 の更新) では
上向き/下向きの結果は真値を挟む包含にはならないので，
区間演算としての利用は呼び出し側で定式化してください．

## テストと benchmark

```bash
mkdir -p sandbox/bin

# 検証 (AVX-512)
g++ -I. -I${MKLROOT}/include -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
-mavx512f -DVCP_USE_AVX512 \
sandbox/tests/test_rdblas.cpp \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 \
-lpthread -lm -ldl -lmpfr -fopenmp \
-o sandbox/bin/test_rdblas

./sandbox/bin/test_rdblas

# 性能比較 (rdgemm/rdsyrk/rdtrsm/rdtrmm vs MKL)
g++ -I. -I${MKLROOT}/include -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
-mavx512f -DVCP_USE_AVX512 \
sandbox/tests/bench_rdblas.cpp \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 \
-lpthread -lm -ldl -lmpfr -fopenmp \
-o sandbox/bin/bench_rdblas

./sandbox/bin/bench_rdblas 2000
```

AVX2 backend は `-mavx512f -DVCP_USE_AVX512` の代わりに `-mavx2 -mfma`，
no-SIMD fallback はどちらも付けずにコンパイルします (3 backend とも
test_rdblas が PASS することを確認済み)．

dblas wrapper は `sandbox/tests/test_dblas.cpp` で検証します
(全 34 関数について，4 つの丸めモード下で wrapper の結果が rdblas の
直接呼び出しと bit 一致すること，丸めモードが保存されることを確認):

```bash
g++ -I. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
-mavx512f -DVCP_USE_AVX512 \
sandbox/tests/test_dblas.cpp \
-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed \
-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 \
-lpthread -lm -ldl -lmpfr -fopenmp \
-o sandbox/bin/test_dblas

./sandbox/bin/test_dblas
```

## License

`rdblas_level1.hpp` / `rdblas_level2.hpp` / `rdblas_level3.hpp` は
reference LAPACK 3.12.1 配布物に含まれる reference BLAS
(Copyright (c) The University of Tennessee 他, modified BSD license) の
派生物です．license 全文と改変内容は `vlapack/LICENSE_LAPACK.txt` を，
由来は各 header 冒頭の notice を参照してください．

## 注意

- 対象は点 double 行列です (区間行列 backend ではありません)．
- rdgemm (alpha != 1 / beta != 0)・rdsyrk・rdtrmm などは m x n 程度の
  temp buffer を確保します．巨大行列ではメモリ使用量に注意してください．
- 丸めモード変更がコンパイラ最適化で無効化されないことは
  `test_rdblas` の up/down 全要素差 test で backend ごとに確認しています．
