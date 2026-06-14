# tblas — テンプレート版 BLAS (T 型要件と実装)

`tblas` は `vblas/rdblas` と同じ routine 一覧を template 化した header-only library
です．本文書は要素型 `T` に要求する演算 (concept 相当) と実装方針を定めます．
**全関数の引数仕様は `tblas/REFERENCE_tblas.md` を参照．**
**丸めモードは一切意識しません** (rounding_mode 引数なし，`fesetround` /
RoundingGuard を使わない)．丸めモード指定が必要な場合は `vblas/rdblas` を
使ってください．

```cpp
#include "tblas/tblas.hpp"
#include <kv/dd.hpp>

// C := alpha*A*B + beta*C を kv::dd で計算 (column-major, 0-based)
tgemm('N', 'N', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);   // T は引数から推論
tgemm<kv::dd>('N', 'N', m, n, k, alpha, A, lda, B, ldb, beta, C, ldc); // 明示も可
```

## file 構成

| file | 内容 |
|---|---|
| `tblas.hpp` | umbrella header (これだけ include すればよい) |
| `tblas_common.hpp` | 共通基盤 (`tblas_detail`: 引数 check / scale / 展開 copy / OpenMP shim) |
| `tblas_level1.hpp` | Level 1 (`tcopy` `tswap` `itamax` `tscal` `taxpy` `tdot` `tasum` `tnrm2` `trot` `trotg` `trotm` `trotmg`) |
| `tblas_level2.hpp` | Level 2 (`tgemv` `tgbmv` `tsymv` `tsbmv` `tspmv` `ttrmv` `ttbmv` `ttpmv` `ttrsv` `ttbsv` `ttpsv` `tger` `tsyr` `tspr` `tsyr2` `tspr2`) |
| `tblas_level3.hpp` | Level 3 (`tgemm` `tgemmtr` `tsymm` `tsyrk` `tsyr2k` `ttrmm` `ttrsm`) |

## 提供予定 routine (rdblas と同一集合，d → t に読み替え)

| Level | routine |
|---|---|
| 1 | `tscal` `taxpy` `tdot` `tasum` `tnrm2` `trot` `trotg` `trotm` `trotmg` `tcopy` `tswap` `itamax` |
| 2 | `tgemv` `tgbmv` `tsymv` `tsbmv` `tspmv` `ttrmv` `ttbmv` `ttpmv` `ttrsv` `ttbsv` `ttpsv` `tger` `tsyr` `tspr` `tsyr2` `tspr2` |
| 3 | `tgemm` `tgemmtr` `tsymm` `tsyrk` `tsyr2k` `ttrmm` `ttrsm` |

全て `template <typename T>` の関数 template (例: `tgemm<T>`)．
引数順は rdblas から末尾の `rounding_mode` を除いたもの (= reference BLAS と同順)．
column-major / 0-based / `int` index，不正引数は `std::invalid_argument`，
scalar 引数は `const T&` 渡し (mpfr 等の copy cost 回避) とします．

## T 型要件

要件は 4 段階に分け，routine ごとに必要な段階のみを要求します．
「下位段階を含む」累積構造です (R2 は R1 を含む，など)。

### R1: 基本要件 — 全 routine 共通

| # | 要件 | 用途 |
|---|---|---|
| 1 | copy 構築・copy 代入・破棄 (値セマンティクス) | 全 routine．`memcpy`/`memset` は使わず代入 loop で実装するため trivially copyable は不要 |
| 2 | `int` からの構築 `T(0)`, `T(1)`, `T(4096)` 等 | zero clear，単位元，trotmg の定数．**double literal からの構築は要求しない** (定数は全て `T(int)` と四則で構成する) |
| 3 | 二項 `+`, `-`, `*` (T × T → T) | 積和の本体 |
| 4 | 単項 `-` | trot / trotg / trsv 等の符号反転 |
| 5 | `==`, `!=` (T × T → bool) | `alpha == T(0)`, `beta == T(1)` の fast path 判定，trotm の flag 判定 |

- 複合代入 `+=` `*=` 等は**要求しない** (二項演算 + 代入のみで実装する)．
- `std::fma` 相当も**要求しない**．rdblas で `std::fma(a, b, c)` の箇所は
  `a * b + c` に置き換える (double で結果が rdblas と bit 一致しない場合がある
  ことを許容する)．
- T と `double`/`int` の混合演算は**要求しない**．literal は必ず `T(...)` で包む．

### R2: 除算 `/` (T × T → T)

| routine | 用途 |
|---|---|
| `ttrsv` `ttbsv` `ttpsv` `ttrsm` | 三角解 (non-unit diag のとき対角での除算) |
| `tnrm2` | overflow/underflow 回避の scaling (`v / scale`) |
| `trotg` `trotmg` | 回転生成 |

### R3: 順序比較 `<`, `<=`, `>`, `>=` (T × T → bool)

| routine | 用途 |
|---|---|
| `itamax` | 最大絶対値要素の探索 |
| `tnrm2` | `scale < v` の scaling 分岐 |
| `trotg` | `fabs(a) > fabs(b)`，`roe < 0` 等の分岐 |
| `trotmg` | 大小比較・rescaling loop |

### R4: 数学関数 — **ADL で探索** (`using std::abs; abs(x);` 形式)

| 関数 | 要件 | routine |
|---|---|---|
| `abs(T)` → T | 絶対値 | `itamax` `tasum` `tnrm2` `trotg` `trotmg` |
| `sqrt(T)` → T | 平方根 | `tnrm2` `trotg` |

呼び出しは必ず

```cpp
using std::abs;   // using std::sqrt;
abs(x);           // 修飾なし → 組み込み型は std::，user 型は ADL で解決
```

の形に統一し，`std::abs(x)` と直接書かない．これにより `kv::dd` /
`kv::mpfr<N>` / `kv::interval<U>` が `namespace kv` に持つ `abs`, `sqrt` が
そのまま使われる．

## routine × 要件 対応表

| routine | R1 | R2 (`/`) | R3 (`<`) | R4 (`abs`/`sqrt`) |
|---|---|---|---|---|
| `tcopy` `tswap` | ✓ | | | |
| `tscal` `taxpy` `tdot` `trot` `trotm` | ✓ | | | |
| `tgemv` `tgbmv` `tsymv` `tsbmv` `tspmv` `ttrmv` `ttbmv` `ttpmv` | ✓ | | | |
| `tger` `tsyr` `tspr` `tsyr2` `tspr2` | ✓ | | | |
| `tgemm` `tgemmtr` `tsymm` `tsyrk` `tsyr2k` `ttrmm` | ✓ | | | |
| `ttrsv` `ttbsv` `ttpsv` `ttrsm` | ✓ | ✓ | | |
| `tasum` | ✓ | | | ✓ (`abs`) |
| `itamax` | ✓ | | ✓ | ✓ (`abs`) |
| `tnrm2` | ✓ | ✓ | ✓ | ✓ (`abs` `sqrt`) |
| `trotg` | ✓ | ✓ | ✓ | ✓ (`abs` `sqrt`) |
| `trotmg` | ✓ | ✓ | ✓ | ✓ (`abs`) |

## 想定する具体型と適合性

| 型 | R1 | R2 | R3 | R4 | 備考 |
|---|---|---|---|---|---|
| `float` `double` `long double` | ✓ | ✓ | ✓ | ✓ | `std::abs` / `std::sqrt` |
| `kv::dd` | ✓ | ✓ | ✓ | ✓ | `kv` namespace の `abs`/`sqrt` を ADL で解決 |
| `kv::mpfr<N>` | ✓ | ✓ | ✓ | ✓ | 同上 |
| `kv::interval<U>` | ✓ | ✓ | △ | ✓ | 四則・`abs`・`sqrt` は包含保証付き．順序比較は区間が重なると分岐が数学的に不定なため，**R3/R4 を要する分岐系 routine (`itamax` `tnrm2` `trotg` `trotmg`) は区間型では非推奨** (R1/R2 のみの routine は安全に使用可) |

区間型の対応状況 (テスト済み):

- **`kv::interval<double>`** と **`kv::interval<kv::dd>`** の両方で，R1/R2 のみの
  routine (`tgemm` `tgemv` `ttrsv` `ttrsm` `taxpy` `tdot` など) と `tasum`
  (abs のみで分岐なし，1-norm の包含を返す) が動作することを
  `sandbox/tests/test_tblas_types.cpp` で検証済み
  (整数成分での幅 0 の exact 計算と，trsv/trsm の residual の 0 包含)．
- `kv::interval<kv::dd>` を使う場合は `<kv/rdd.hpp>` (dd の方向丸め) の include
  が必要 (`interval<double>` は `<kv/rdouble.hpp>`)．
- 分岐系 4 routine (`itamax` `tnrm2` `trotg` `trotmg`) は区間型では使わないこと．
  kv の区間比較は「確実に小さい」(`x.sup < y.inf`) の意味論なので重なった区間では
  分岐が選別できず，`tnrm2` は scale が 0 を含む区間になると内部の除算で
  `std::domain_error` を投げうる．
- tlapack は区間型に**対応しない** (理由は tlapack/README_tlapack.md の
  「区間型は対象外」を参照)．

## 設計上の注意 (要件から導かれる実装方針)

1. **定数の構成**: `trotmg` の `gam = 4096`, `gamsq = gam * gam`,
   `rgamsq = T(1) / gamsq` は全て `T(int)` と四則から作る
   (reference の `5.9604645e-8` ≒ `1/4096^2` を除算で代替)．
   これらの rescaling 閾値は本来 double 向け tuning なので，
   多倍長型では rescaling の発動頻度が変わるが結果の正しさには影響しない．
2. **fast path**: `alpha == T(0)` 等の等値判定は R1 の `==` のみで行い，
   区間型でも問題ない (区間の `==` は両端一致の判定)．
3. **一時変数**: `T tmp = T(0);` のように必ず明示的に初期化し，
   default 構築 `T()` の値には依存しない．
4. **swap**: `tswap` は `std::swap<T>` または 3 代入で実装 (memcpy 不可)．
5. **戻り値**: `tdot` `tasum` `tnrm2` は `T` を値返し，`itamax` は
   rdblas 同様 0-based `int` (空のとき -1)．
6. **OpenMP**: 並列化する場合，reduction は `T` の `+` のみで書ける形
   (手動 reduction) にする．`#pragma omp declare reduction` は使わない．
7. **namespace / file 構成**: rdblas に合わせて関数は global scope，
   補助は `tblas_detail` namespace．
   `tblas/tblas.hpp` を umbrella とし `tblas_common.hpp` /
   `tblas_level1.hpp` / `tblas_level2.hpp` / `tblas_level3.hpp` に分割．

## 実装の詳細 (rdblas との対応)

- **Level 1 / Level 2**: rdblas の loop の忠実な移植．`std::fma(a, b, c)` は
  `a * b + c` に展開し，literal は全て `T(int)` から構成する
  (double で rdblas と bit 一致しない routine がある)．
- **Level 3**:
  - `tgemm` は reference dgemm と同じ列方向 algorithm (rmatmul のような
    SIMD / blocking はしない)．
  - `tsymm` / `ttrmm` は対称・三角格納を full dense に展開して (copy のみ，
    unit diag は `T(1)`，三角の外は exact `T(0)`) `tgemm` に帰着する
    (rdblas と同方式．exact 0 との積和なので結果に影響せず，区間型でも
    幅は増えない)．
  - `tsyrk` / `tsyr2k` / `tgemmtr` は uplo 三角部分のみを dot 形式で直接計算．
  - `ttrsm` は rdblas の leaf solver (left: 列ごとの独立 solve，right: 列 sweep)
    の非再帰移植．
- **OpenMP**: 列方向に独立な loop のみ `if` 句付きで並列化する
  (`tgemv` `tsymv` `tger` `tsyr` `tsyr2` と Level 3 全 routine の主 loop，
  `ttrsm` right side は列間に依存があるため逐次)．閾値はスカラー演算回数で
  `TBLAS_OMP_THRESHOLD` (既定 2e4)，`-fopenmp` なしでも compile 可能．
- **引数規約**: reference BLAS と同順で rounding_mode なし，scalar は
  `const T&` 渡し，column-major / 0-based / `int` index，不正引数は
  `std::invalid_argument` (例外無効時は `abort`)．

## テスト

`sandbox/tests/test_tblas.cpp` (T=double を rdblas と比較，全 routine・全
option 組合せ・inc 正負・alpha/beta 特殊値で 282 checks) と
`sandbox/tests/test_tblas_types.cpp` (kv::dd / kv::mpfr<106>: double との一致と
residual の高精度性，kv::interval<double>: dd 参照値の包含と residual の 0 包含)
で検証済み．

```bash
make -C sandbox run_tblas   # build + 実行 (test_tblas, test_tblas_types)
```
