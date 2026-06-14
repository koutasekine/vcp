# dlapack — 通常の LAPACK interface (rdlapack の wrapper)

`dlapack.hpp` は `vcp::rdlapack` を呼び出して構築した「通常の」double LAPACK です．
関数名・引数は Fortran LAPACK の C からの呼び出し規約に一致します
(末尾 underscore，全引数 pointer 渡し，INTEGER は 32bit `int` / LP64)．
dblas (`vcp/vblas/dblas.hpp`) の LAPACK 版にあたり，階層は完全に対称です:

```text
rmatmul ── vcp::rdblas ──┬── dblas    (BLAS  の Fortran 互換 wrapper，現在の丸めモードで計算)
                         └── vcp::rdlapack ── dlapack (LAPACK の Fortran 互換 wrapper，現在の丸めモードで計算)
```

`vcp::rdblas` / `vcp::rdlapack` は `vcp` 名前空間に収められた C++ 関数です．
`dblas` / `dlapack` は Fortran 互換の `extern "C"` グローバルシンボルであり，
名前空間には入っていません．
`vcp::pdblas`・`vcp::pidblas`・`vcp::pddblas` から使う場合は，内部の共通ヘッダー
`vcp/dblas_dlapack.hpp` が `USE_VCP_LAPACK` の有無に応じて `dlapack.hpp` と
外部 LAPACK 宣言を切り替えます．

```cpp
#include <vcp/vlapack/dlapack.hpp>

std::fesetround(FE_UPWARD);  // 以後の dlapack 呼び出しは全演算が上向き丸めで計算される
int info;
dgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
std::fesetround(FE_TONEAREST);
```

- 各関数は呼び出し時点の丸めモードを取得し，
  `FE_DOWNWARD` なら `-1`，`FE_UPWARD` なら `1`，それ以外なら `0` を
  `vcp::rdlapack` の末尾引数 `rounding_mode` として渡します (dblas と同じ規約)．
- `KV_FASTROUND` 有効時は，`kv::hwround::roundup()` / `rounddown()` が直接変更する
  制御レジスタ (x86 では MXCSR，AArch64/ARM では FPCR/FPSCR) から丸め方向を読みます．
  `KV_FASTROUND` 無効時は通常の `std::fegetround()` を使います．
- `vcp::rdlapack` の保証がそのまま引き継がれます: 全浮動小数点演算 (四則演算と平方根，
  BLAS 部分も scalar 部分も，OpenMP の全 worker thread も) が同じ丸めモードで
  実行され，終了後に各 thread の丸めモードは呼び出し前の状態へ復元されます．
- `vcp::rdlapack` と同様，これは「全演算を同一丸めモードで計算する」だけであり，
  **精度保証付き数値計算ではありません** (`docs/rdlapack.md` の注意を参照)．

## 既存の BLAS の取り扱い (重要)

### dlapack は既存の BLAS と干渉しない

dlapack の内部は `dlapack → vcp::rdlapack → vcp::rdblas` という **C++ の直接呼び出し**
(マングルされた C++ 関数名，多くは inline 展開) で完結しており，`dgemm_` などの
**BLAS シンボルを参照も提供もしません**．リンカのシンボル解決を経由する箇所が
ないため:

- MKL や libblas をリンクしても **dlapack の内部計算には一切影響しません**
  (使われない，ではなく「呼ぶ経路が存在しない」)．
- 同一 binary 内で，ユーザコードが直接呼ぶ `dgemm_` / `cblas_dgemm` は
  従来どおり既存 BLAS に解決されます (dblas を emit しない限り)．
  この共存は `test_dlapack` で検証済みです．

```text
                 ┌─ dlapack (dgesv_ …) ─→ vcp::rdlapack ─→ vcp::rdblas   ← C++ 直接呼び出し (差し替え不能)
あなたのプログラム ┤
                 └─ dgemm_ などの直接呼び出し ─→ シンボル解決で決まる:
                       ├─ dblas を emit していれば → dblas (→ vcp::rdblas)
                       └─ していなければ → リンクした既存 BLAS (MKL 等)
```

### 既存の BLAS の速度の恩恵は受けない

上記の構造の帰結として，dlapack の O(N^3) 部分は常に rdblas (rmatmul 核) で
実行され，MKL 等の最適化を取り込む手段はありません．目安 (本環境での実測):
gemm 系は rdblas の方が速い (丸め制御付きで MKL の約 1.5 倍) 一方，
trmm/trsm 系は MKL の 5〜6 割，非対称固有値 (dgeev) は algorithm 簡略化
(dlahqr のみ) のため大規模では明確に遅くなります．

### 使い分けの指針

| 目的 | 構成 |
|---|---|
| 全演算の丸めモード保証 (LAPACK 演算) | **dlapack** (または rounding_mode を明示する rdlapack) |
| 既存 LAPACK + BLAS 部分のみ丸め保証 | 既存 liblapack + **dblas** (検証済み，docs/rdblas.md 参照) |
| 性能優先 (丸め保証なし) | 既存 liblapack + MKL 等 (従来どおり) |
| 丸め保証付き LAPACK と高速 BLAS の同居 | dlapack + MKL を同一 binary にリンク (dblas は emit しない) |

### LAPACK シンボル側の衝突に注意

dlapack が export するのは LAPACK シンボル (`dgesv_` 等) です．liblapack や
MKL (LAPACK シンボルも含む) と同時リンクする場合，同名シンボルは
リンク構成で解決が決まります (既定では実行ファイル側 = dlapack の weak symbol
が優先)．また dlapack が提供しない routine (複素数/単精度，`*con`/`*rfs` 等)
は liblapack 側に解決されるため，**1 つのプログラム内で dlapack の `dgesv_` と
liblapack の `dgecon_` が混在する**構成になり得ます．意図した構成かどうか
`LD_DEBUG=bindings` 等で確認することを推奨します．

## LAPACK との互換性 (実装方法)

wrapper は計算を持たず，次の変換のみを行います:

1. **丸めモード**: 現在の丸めモード → rdlapack の `rounding_mode` (上記規約)
2. **Fortran ABI**: 全引数 pointer 渡し，`INFO` は出力引数に格納
3. **ipiv の 1-based 変換**: `dgetrf_`/`dgesv_`/`dgbtrf_`/`dgbsv_` 系は
   rdlapack の 0-based ipiv を wrapper が 1-based に変換して入出力する
   (LAPACK と完全互換)．`dsytrf_` 系は rdlapack 自体が LAPACK 互換
   (1-based 符号付き) のため pass-through．`dlaswp_` の k1/k2/ipiv も
   1-based のまま受け取り内部で変換する
4. **WORK/LWORK**: 作業領域は rdlapack が内部確保するため WORK は使用しない．
   `LWORK = -1` (workspace query) のときは計算せず推奨 size を `WORK(1)` に
   返す．それ以外の LWORK の値は**検査せず無視**する (LAPACK は不足時に
   INFO < 0 を返すが，dlapack では不足という概念がない)
5. **引数 error**: rdlapack が投げる `std::invalid_argument` を catch して
   `INFO = -1` に変換する (LAPACK と異なり，どの引数が不正かは特定しない)．
   `INFO > 0` (特異・非収束等) は LAPACK と同じ意味
6. **symbol の強制発行**: Fortran object と link する場合は 1 つの翻訳単位で
   `VLAPACK_DLAPACK_EMIT_SYMBOLS` を定義して include する
   (`__attribute__((used))` 付きの table で全 symbol を weak 発行)

## 提供 routine (78)

rdlapack が提供する routine の Fortran 互換版です:

| 分類 | routine |
|---|---|
| LU / 三角 | `dgesv_` `dgetf2_` `dgetrf_` `dgetrs_` `dgetri_` `dtrti2_` `dtrtri_` `dtrtrs_` `dlauu2_` `dlauum_` |
| Cholesky | `dposv_` `dpotf2_` `dpotrf_` `dpotrs_` `dpotri_` |
| 対称不定値 | `dsysv_` `dsytf2_` `dsytrf_` `dsytrs_` |
| 帯 | `dgbsv_` `dgbtf2_` `dgbtrf_` `dgbtrs_` `dpbsv_` `dpbtf2_` `dpbtrf_` `dpbtrs_` |
| QR/LQ/QL/最小二乗 | `dgels_` `dgeqr2_` `dgeqrf_` `dgelq2_` `dgelqf_` `dorgqr_` `dorglq_` `dorgql_` `dormqr_` `dormlq_` `dormql_` |
| 対称固有値 | `dsyev_` `dsytd2_` `dsytrd_` `dorgtr_` `dormtr_` `dsterf_` `dsteqr_` `dsygs2_` `dsygst_` `dsygv_` |
| SVD | `dgesvd_` `dgebd2_` `dgebrd_` `dorgbr_` `dormbr_` `dbdsqr_` |
| 非対称固有値 | `dgeev_` `dgebal_` `dgebak_` `dgehd2_` `dgehrd_` `dorghr_` `dhseqr_` `dtrevc_` (HOWMNY='A'/'B' のみ) |
| 補助 | `dlamch_` `dlapy2_` `dlapy3_` `dlassq_` `dlartg_` `dlarfg_` `dlascl_` `dlaset_` `dlacpy_` `dlasrt_` `dlaswp_` `dlange_` `dlangb_` `dlansy_` `dlansb_` `dlanst_` |

複素数・単精度，条件数推定 (`*con`)・反復改良 (`*rfs`)・expert driver は
rdlapack 同様に提供しません．

## テスト

4 つの丸めモード下で主要 wrapper の結果が rdlapack の直接呼び出しと bit 一致すること，
ipiv の 1-based 変換，workspace query (LWORK = -1)，丸めモードの保存，
および同一 binary 内での MKL (`cblas_dgemm`) との共存を確認済み．

## License

dlapack 自体は wrapper ですが，呼び出す rdlapack は reference LAPACK 3.12.1 の
派生物です．`vcp/vlapack/LICENSE_LAPACK.txt` を参照してください．
