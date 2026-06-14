# rdlapack — 丸めモード指定付き double LAPACK

`rdlapack` は，reference LAPACK 3.12.1 (BSD license, <https://www.netlib.org/lapack/>)
の double 精度実数 routine を，BLAS 呼び出しを `rdblas` に置き換えながら C++ へ
移植した header-only library です．`rmatmul` / `rdblas` と同様に実験的な位置付けで，
`vcp::matrix` からは自動では使われません．

```cpp
#include <vcp/vlapack/rdlapack.hpp>

// A*X = B を全演算上向き丸めで解く
int info = vcp::rdgesv(n, nrhs, A, lda, ipiv, B, ldb, 1);
```

あるいはファイル先頭で `using namespace vcp;` を宣言すれば，
`vcp::` 修飾なしで呼び出せます:

```cpp
#include <vcp/vlapack/rdlapack.hpp>
using namespace vcp;

int info = rdgesv(n, nrhs, A, lda, ipiv, B, ldb, 1);
```

現在の丸めモードを読み込んで計算する Fortran 互換 wrapper (dlapack) は
`vcp/vlapack/dlapack.hpp` にあります (`docs/dlapack.md` 参照)．

## 名前空間

**rdlapack の全関数は `vcp` 名前空間に収められています．**
`vcp::rdgesv`, `vcp::rdsyev`, `vcp::dlamch` のように修飾するか，
呼び出し側ファイルで `using namespace vcp;` を宣言してください．

FP 演算を含まない補助関数 (`dlamch`, `dlaset`, `dlacpy`, `dlaswp`,
`dlasrt`, `iladlr`, `iladlc`) も同様に `vcp` 名前空間内です．

`dlapack.hpp` の Fortran 互換 wrapper (`dgesv_` など) は `extern "C"` の
グローバルシンボルであり，名前空間に入っていません (変更なし)．

## 重要な注意

1. **全ての浮動小数点演算が指定した同一の丸めモードで計算されます．**
   各 routine 内の浮動小数点演算 (四則演算と平方根) は，BLAS 部分 (rdblas) も
   LAPACK 固有の scalar 演算部分も含めて，すべて引数 `rounding_mode` で指定した
   丸めモード (1: 上向き `FE_UPWARD`, -1: 下向き `FE_DOWNWARD`, それ以外: 最近点
   `FE_TONEAREST`) で実行されます．計算終了後は，呼び出し thread (および rdblas が
   使用した OpenMP の全 worker thread) の丸めモードが呼び出し前の状態へ復元されます．
2. **これは精度保証付き数値計算ではありません．**
   rdlapack は「すべての演算を同一の丸めモードで計算する」だけです．
   上向き丸めと下向き丸めの結果が真値を挟むこと (包含) は，gesv の解や固有値の
   ような一般の計算では成り立ちません (途中に減算・除算・反復が含まれるため)．
   区間演算による厳密な誤差評価 (精度保証付き数値計算) が必要な場合は，
   別途その定式化のもとで利用してください．

## 命名規約と引数

- FP 演算を含む routine は LAPACK 名の先頭に `r` を付け，**reference LAPACK と
  同順の引数の末尾に `rounding_mode` を 1 つ追加**します．
- FP 演算を含まない routine は LAPACK と同名のままです:
  `dlamch`, `dlaset`, `dlacpy`, `dlaswp`, `dlasrt`, `iladlr`, `iladlc`．
- 行列は column-major，index は 0-based です．例外 (LAPACK と同じ流儀のもの):
  - `rdsytf2`/`rdsytrf`/`rdsytrs`/`rdsysv` の `ipiv` は **1-based 符号付き**
    (負値は 2x2 block)．
  - `rdgebal`/`rdgebak`/`rdgehrd`/`rdorghr`/`rdhseqr`/`rdlahqr` の `ilo`/`ihi` は
    **1-based**．
  - `rdgesv`/`rdgetrf`/`rdgbsv` 系の `ipiv` と `dlaswp` は 0-based．
- LAPACK の `INFO` は **返り値 (int)** です (0: 成功，> 0: LAPACK と同じ意味の
  1-based 位置/個数)．`INFO < 0` に相当する引数 error は `std::invalid_argument`
  を投げます (例外無効時は `abort`)．
- LAPACK の `WORK`/`LWORK` 引数は持ちません．作業領域は内部で `std::vector` により
  最適 size を確保します (workspace query は不要)．

## 提供 routine

| 分類 | driver | 計算 routine |
|---|---|---|
| 連立一次方程式 (一般) | `rdgesv` | `rdgetf2` `rdgetrf2` `rdgetrf` `rdgetrs` `rdgetri` |
| 連立一次方程式 (対称正定値) | `rdposv` | `rdpotf2` `rdpotrf2` `rdpotrf` `rdpotrs` `rdpotri` |
| 連立一次方程式 (対称不定値) | `rdsysv` | `rdsytf2` `rdlasyf` `rdsytrf` `rdsytrs` |
| 連立一次方程式 (帯) | `rdgbsv` `rdpbsv` | `rdgbtf2` `rdgbtrf` `rdgbtrs` `rdpbtf2` `rdpbtrf` `rdpbtrs` |
| 三角行列 | — | `rdtrtrs` `rdtrti2` `rdtrtri` `rdlauu2` `rdlauum` |
| 最小二乗 | `rdgels` | QR/LQ 一式 |
| QR / LQ / QL | — | `rdgeqr2` `rdgeqrf` `rdgelq2` `rdgelqf` `rdorg2r` `rdorgqr` `rdorgl2` `rdorglq` `rdorg2l` `rdorgql` `rdorm2r` `rdormqr` `rdorml2` `rdormlq` `rdorm2l` `rdormql` |
| 対称固有値 | `rdsyev` | `rdsytd2` `rdlatrd` `rdsytrd` `rdorgtr` `rdormtr` `rdsterf` `rdsteqr` |
| 一般化対称固有値 | `rdsygv` | `rdsygs2` `rdsygst` |
| 非対称固有値 | `rdgeev` | `rdgebal` `rdgebak` `rdgehd2` `rdlahr2` `rdgehrd` `rdorghr` `rdlahqr` `rdhseqr` `rdlaln2` `rdtrevc` |
| 特異値分解 | `rdgesvd` | `rdgebd2` `rdlabrd` `rdgebrd` `rdorgbr` `rdormbr` `rdbdsqr` |
| 補助 | — | `rdlassq` `rdlapy2` `rdlapy3` `rdladiv` `rdlartg` `rdlarfg` `rdlarf` `rdlarf1f` `rdlarf1l` `rdlarft` `rdlarfb` `rdlascl` `rdlasr` `rdlae2` `rdlaev2` `rdlas2` `rdlasv2` `rdlanv2` `rdlange` `rdlangb` `rdlansy` `rdlansb` `rdlanst` / `dlamch` `dlaset` `dlacpy` `dlaswp` `dlasrt` `iladlr` `iladlc` |

block size は reference の `ilaenv` の既定値 (GETRF/POTRF/SYTRF 64, QR/SYTRD/
GEBRD/GEHRD 系 32 など) を `vblas_rdlapack_detail::ilaenv` に再現しています．
O(N^3) 部分は rdblas (内部は rmatmul: AVX-512 / AVX2+FMA / NEON / no-SIMD を
compile 時に自動選択) が実行します．

## reference LAPACK との相違点

結果の数学的意味は同じですが，以下は実装が reference と異なります
(丸めの列が変わるため bit 単位では reference/MKL と一致しません):

- `rdbdsqr` は特異値のみ (`ncvt = nru = ncc = 0`) の場合も dqds (`dlasq1` 系)
  ではなく QR 反復で計算します (dlasq 系は未移植)．
- `rdgesvd` は基本経路 (A を直接二重対角化) のみで，m ≫ n / n ≫ m 用の
  QR/LQ 前処理経路 (`mnthr` 分岐) は持ちません (演算量のみ異なる)．
- `rdhseqr` は全 size で `rdlahqr` (double-shift QR) を使います
  (`dlaqr0` 系 multishift は未移植)．大規模行列では reference より遅くなります．
- `rdtrevc` は `howmny = 'A' / 'B'` のみ対応 (`SELECT` 指定は未対応)．
- `rdlasyf` の `DGEMMTR` (LAPACK 3.12 で導入) は LAPACK 3.11 の dlasyf と同じ
  blocked 更新で代替しています．
- `rdsysv` は `dsytrs2` ではなく常に `rdsytrs` を使います．
- 条件数推定 (`*con`)・反復改良 (`*rfs`)・expert driver (`*svx` 等) は未移植です．

## ファイル構成

rdlapack (LAPACK 層) は `vcp/vlapack/` に，依存する rdblas / dblas / rmatmul
(BLAS 層) は兄弟 directory の `vcp/vblas/` に置かれています:

```text
vblas/     BLAS 層 (rmatmul*.hpp, rdblas*.hpp, dblas.hpp)
vlapack/   LAPACK 層 (本 library)
```

- `rdlapack.hpp` — 入口 (これだけ include すればよい)
- `rdlapack_common.hpp` — 共通基盤 (`dlamch`, `ilaenv` 相当, RoundingGuard 再利用．
  `<vcp/vblas/rdblas.hpp>` を include する)
- `rdlapack_aux.hpp` — 補助 routine (Householder/Givens/norm/scale/2x2 固有・特異値)
- `rdlapack_lu.hpp` / `rdlapack_chol.hpp` / `rdlapack_qr.hpp` — LU / Cholesky / QR 系
- `rdlapack_eig.hpp` / `rdlapack_svd.hpp` / `rdlapack_geev.hpp` — 固有値 / SVD / 非対称固有値
- `rdlapack_sy.hpp` / `rdlapack_band.hpp` — Bunch-Kaufman・一般化固有値 / 帯行列

vblas 側の `rdblas*.hpp` の公開関数はすべて `namespace vcp { }` で囲まれています
(`rmatmul*.hpp` は内部実装であり名前空間の外に置かれています)．

## テスト

主要 driver について
(1) 最近点丸めでの残差と MKL LAPACKE との比較，
(2) 上向き/下向きで結果が要素レベルで異なること (丸め変更が全演算に効いている)，
(3) 呼び出し時の丸めモードに依らず同じ `rounding_mode` 指定で bit 同一の結果に
なること (ambient 独立性) と呼び出し後の丸めモード復元，を検証済み．

no-SIMD fallback は `-mavx512f -DVCP_USE_AVX512` を外して compile します
(AVX-512 / no-SIMD の両 backend で ALL PASS を確認済み)．

## License

本 directory の `rdlapack*.hpp` は reference LAPACK 3.12.1
(Copyright (c) 1992-2023 The University of Tennessee and The University of
Tennessee Research Foundation / 2000-2023 The University of California
Berkeley / 2006-2023 The University of Colorado Denver, modified BSD license)
の派生物です．LAPACK license の全文・条件一覧・免責条項と，改変内容の一覧は
**`vcp/vlapack/LICENSE_LAPACK.txt`** に収録しています (license の条件である
「ソース再配布時の著作権表示・条件・免責の保持」はこのファイルと各 header
冒頭の由来表示で満たします)．各 header には由来と改変者を明記した
notice を付しています．本移植自体は VCP Library の BSD 3-clause license に
従います (LAPACK license と互換)．
