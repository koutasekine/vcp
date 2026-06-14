# tblas — テンプレート版 BLAS (T 型要件と実装)

`tblas` は `vcp/vblas/rdblas` と同じ routine 一覧を template 化した header-only library
です．本文書は要素型 `T` に要求する演算 (concept 相当) と実装方針を定めます．
**全関数の引数仕様は このドキュメントの「関数リファレンス」節 を参照．**
**丸めモードは一切意識しません** (rounding_mode 引数なし，`fesetround` /
RoundingGuard を使わない)．丸めモード指定が必要な場合は `vcp/vblas/rdblas` を
使ってください．

```cpp
#include <vcp/tblas/tblas.hpp>
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
   `vcp/tblas/tblas.hpp` を umbrella とし `tblas_common.hpp` /
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

以下の 2 点で検証済み:

- **T=double を rdblas と比較** (282 checks)：全 routine・全 option 組合せ・inc 正負・alpha/beta 特殊値．
- **kv::dd / kv::mpfr<106>**：double との一致と residual の高精度性．
  **kv::interval<double>**：dd 参照値の包含と residual の 0 包含．


---

# 関数リファレンス


`vcp/tblas/tblas.hpp` が提供する全 37 routine の仕様．設計方針と要素型 T の要件
(R1–R4) は このドキュメントの「設計方針」節 を参照．

## 共通仕様

- 全関数は `template <typename T>` の関数 template．`T` は引数から推論される
  (`tgemm(...)` のように呼べばよい．`tgemm<kv::dd>(...)` と明示してもよい)．
- **行列は column-major**: `A(i,j)` は `A[i + lda * j]`．leading dimension
  `lda` は列の先頭間隔で，`lda >= max(1, 行数)` が必要．
- **index は 0-based**．`itamax` の返り値も 0-based (Fortran BLAS と異なる)．
- **増分 (incx / incy)**: vector x の第 i 要素は `x[start + i*incx]`，
  `start = incx > 0 ? 0 : (n-1)*(-incx)` (BLAS の負増分規約)．
  `incx == 0` は不正引数 (`tscal` `itamax` `tasum` `tnrm2` は `incx <= 0` のとき
  何もしない / 0 / -1 を返す BLAS 規約)．
- **scalar 引数は `const T&` 渡し**，vector/行列は pointer 渡し．
- **不正引数**は `std::invalid_argument` を投げる (例外無効時は `abort`)．
- **丸めモードは扱わない**．丸めモード指定が必要なら `vcp/vblas/rdblas` を使う．
- `std::fma` は使わない (`a*b + c` は乗算 1 丸め + 加算 1 丸め)．
- BLAS 規約の quick return: `beta == T(0)` のとき出力は読まずに 0 埋め
  (NaN があっても伝播しない)，`beta == T(1)` のとき scaling を省略，
  `alpha == T(0)` のとき入力行列を読まない．
- 演算量が `TBLAS_OMP_THRESHOLD` (既定 2e4) を超える独立 loop は OpenMP 並列．

### option 文字 (大文字小文字どちらも可)

| 引数 | 値 | 意味 |
|---|---|---|
| `trans` | `'N'` / `'T'`,`'C'` | op(A) = A / op(A) = A^T (実数なので 'C' = 'T') |
| `uplo` | `'U'` / `'L'` | 上三角 (または上三角部のみ参照) / 下三角 |
| `diag` | `'N'` / `'U'` | 対角は格納値 / 単位対角 (対角は参照されず 1 とみなす) |
| `side` | `'L'` / `'R'` | 左から掛ける / 右から掛ける |

### 格納形式

- **band 格納** (tgbmv): `AB[ku + i - j + ldab*j] = A(i,j)`
  (`max(0,j-ku) <= i <= min(m-1,j+kl)`，`ldab >= kl+ku+1`)．
- **対称/三角 band 格納** (tsbmv / ttbmv / ttbsv): uplo='U' なら
  `AB[k + i - j + ldab*j] = A(i,j)` (`j-k <= i <= j`)，uplo='L' なら
  `AB[i - j + ldab*j] = A(i,j)` (`j <= i <= j+k`)，`ldab >= k+1`．
- **packed 格納** (tspmv / tspr / tspr2 / ttpmv / ttpsv): uplo='U' なら
  `AP[i + j(j+1)/2] = A(i,j)` (`i <= j`)，uplo='L' なら
  `AP[i - j + jn - j(j-1)/2] = A(i,j)` (`i >= j`)．長さ n(n+1)/2．

---

## Level 1

### tcopy — vector の copy
```cpp
template <typename T>
void tcopy(const int n, const T* x, const int incx, T* y, const int incy);
```
y := x．FP 演算なし (要件 R1 のみ)．`n <= 0` のとき何もしない．

### tswap — vector の交換
```cpp
template <typename T>
void tswap(const int n, T* x, const int incx, T* y, const int incy);
```
x <-> y．FP 演算なし (R1)．

### itamax — 最大絶対値要素の index
```cpp
template <typename T>
int itamax(const int n, const T* x, const int incx);
```
|x_i| が最大となる最小の **0-based** index を返す．
`n <= 0` または `incx <= 0` のとき **-1** を返す．要件: R1+R3+R4 (`abs`, `>`)．
区間型では非推奨 (比較分岐のため)．

### tscal — scalar 倍
```cpp
template <typename T>
void tscal(const int n, const T& alpha, T* x, const int incx);
```
x := alpha\*x．`incx <= 0` または `alpha == T(1)` のとき何もしない (R1)．

### taxpy — scalar 倍加算
```cpp
template <typename T>
void taxpy(const int n, const T& alpha, const T* x, const int incx, T* y, const int incy);
```
y := alpha\*x + y．`alpha == T(0)` のとき何もしない (R1)．

### tdot — 内積
```cpp
template <typename T>
T tdot(const int n, const T* x, const int incx, const T* y, const int incy);
```
x^T y を逐次和で返す．`n <= 0` のとき T(0) (R1)．

### tasum — 絶対値和
```cpp
template <typename T>
T tasum(const int n, const T* x, const int incx);
```
sum |x_i| を返す．`n <= 0` または `incx <= 0` のとき T(0)．要件: R1+R4 (`abs`)．

### tnrm2 — Euclid norm
```cpp
template <typename T>
T tnrm2(const int n, const T* x, const int incx);
```
||x||_2 を返す．overflow/underflow 回避の scale/ssq 方式 (reference dnrm2 と同方式)．
要件: R1–R4 全て．区間型では非推奨．

### trot — 平面回転の適用
```cpp
template <typename T>
void trot(const int n, T* x, const int incx, T* y, const int incy, const T& c, const T& s);
```
(x_i, y_i) := (c\*x_i + s\*y_i, c\*y_i - s\*x_i) (R1)．

### trotg — Givens 回転の生成
```cpp
template <typename T>
void trotg(T* a, T* b, T* c, T* s);
```
[c s; -s c] [a; b] = [r; 0] となる c, s を生成する (reference drotg と同方式)．
出力: `*a` = r，`*b` = 再構成用の z，`*c`, `*s`．要件: R1–R4．区間型では非推奨．

### trotm — modified Givens 回転の適用
```cpp
template <typename T>
void trotm(const int n, T* x, const int incx, T* y, const int incy, const T* param);
```
(x_i, y_i) := H (x_i, y_i)^T．`param[0]` = flag (−2: 恒等で何もしない，
−1: H 全要素，0: 対角 1 + `param[2]`,`param[3]`，1: 反対角 ±1 + `param[1]`,`param[4]`)，
`param[1..4]` = h11, h21, h12, h22．flag は `trotmg` が生成する −2/−1/0/1 を仮定 (R1)．

### trotmg — modified Givens 回転の生成
```cpp
template <typename T>
void trotmg(T* d1, T* d2, T* x1, const T& y1, T* param);
```
sqrt(d1)\*x1, sqrt(d2)\*y1 の回転消去を行う H と更新後の d1, d2, x1 を生成する
(reference drotmg と同方式)．`param` の意味は trotm と同じ．
rescaling 閾値は gam = 4096 から `T` の演算で構成 (reference の double 定数
`5.9604645e-8` とは末位が僅かに異なる)．要件: R1–R4．区間型では非推奨．

---

## Level 2

### tgemv — 一般行列 vector 積
```cpp
template <typename T>
void tgemv(const char trans, const int m, const int n,
	const T& alpha, const T* A, const int lda, const T* x, const int incx,
	const T& beta, T* y, const int incy);
```
y := alpha\*op(A)\*x + beta\*y．A は m×n．x の長さは op(A) の列数
(trans='N' なら n，'T' なら m)，y の長さは行数 (R1)．

### tgbmv — band 行列 vector 積
```cpp
template <typename T>
void tgbmv(const char trans, const int m, const int n, const int kl, const int ku,
	const T& alpha, const T* A, const int lda, const T* x, const int incx,
	const T& beta, T* y, const int incy);
```
y := alpha\*op(A)\*x + beta\*y．A は band 格納 (下帯域 kl，上帯域 ku，
`lda >= kl+ku+1`) (R1)．

### tsymv / tsbmv / tspmv — 対称行列 vector 積
```cpp
template <typename T>
void tsymv(const char uplo, const int n, const T& alpha, const T* A, const int lda,
	const T* x, const int incx, const T& beta, T* y, const int incy);
template <typename T>
void tsbmv(const char uplo, const int n, const int k, const T& alpha, const T* A, const int lda,
	const T* x, const int incx, const T& beta, T* y, const int incy);
template <typename T>
void tspmv(const char uplo, const int n, const T& alpha, const T* AP,
	const T* x, const int incx, const T& beta, T* y, const int incy);
```
y := alpha\*A\*x + beta\*y．A は n×n 対称で uplo の三角のみ参照する．
tsbmv は対称 band (帯域 k，`lda >= k+1`)，tspmv は packed 格納 (R1)．

### ttrmv / ttbmv / ttpmv — 三角行列 vector 積
```cpp
template <typename T>
void ttrmv(const char uplo, const char trans, const char diag, const int n,
	const T* A, const int lda, T* x, const int incx);
template <typename T>
void ttbmv(const char uplo, const char trans, const char diag, const int n, const int k,
	const T* A, const int lda, T* x, const int incx);
template <typename T>
void ttpmv(const char uplo, const char trans, const char diag, const int n,
	const T* AP, T* x, const int incx);
```
x := op(A)\*x (in-place)．A は n×n 三角 (ttbmv は三角 band，ttpmv は packed) (R1)．

### ttrsv / ttbsv / ttpsv — 三角方程式の求解
```cpp
template <typename T>
void ttrsv(const char uplo, const char trans, const char diag, const int n,
	const T* A, const int lda, T* x, const int incx);
template <typename T>
void ttbsv(const char uplo, const char trans, const char diag, const int n, const int k,
	const T* A, const int lda, T* x, const int incx);
template <typename T>
void ttpsv(const char uplo, const char trans, const char diag, const int n,
	const T* AP, T* x, const int incx);
```
op(A)\*x = b を解き x を解で上書きする (入力の x が b)．
特異性検査は行わない (`diag='N'` で対角 0 なら 0 除算)．要件: R1+R2 (`/`)．

### tger — rank-1 更新 (一般)
```cpp
template <typename T>
void tger(const int m, const int n, const T& alpha, const T* x, const int incx,
	const T* y, const int incy, T* A, const int lda);
```
A := alpha\*x\*y^T + A．A は m×n (R1)．

### tsyr / tspr — 対称 rank-1 更新
```cpp
template <typename T>
void tsyr(const char uplo, const int n, const T& alpha, const T* x, const int incx,
	T* A, const int lda);
template <typename T>
void tspr(const char uplo, const int n, const T& alpha, const T* x, const int incx, T* AP);
```
A := alpha\*x\*x^T + A．uplo の三角のみ更新する (R1)．

### tsyr2 / tspr2 — 対称 rank-2 更新
```cpp
template <typename T>
void tsyr2(const char uplo, const int n, const T& alpha, const T* x, const int incx,
	const T* y, const int incy, T* A, const int lda);
template <typename T>
void tspr2(const char uplo, const int n, const T& alpha, const T* x, const int incx,
	const T* y, const int incy, T* AP);
```
A := alpha\*x\*y^T + alpha\*y\*x^T + A．uplo の三角のみ更新する (R1)．

---

## Level 3

### tgemm — 一般行列積
```cpp
template <typename T>
void tgemm(const char transa, const char transb, const int m, const int n, const int k,
	const T& alpha, const T* A, const int lda, const T* B, const int ldb,
	const T& beta, T* C, const int ldc);
```
C := alpha\*op(A)\*op(B) + beta\*C．op(A): m×k，op(B): k×n，C: m×n．
`lda >= max(1, transa='N' ? m : k)` (B も同様)．reference dgemm と同じ列方向
algorithm (SIMD/blocking なし)．列方向に OpenMP 並列 (R1)．

### tgemmtr — 三角部分のみの一般行列積
```cpp
template <typename T>
void tgemmtr(const char uplo, const char transa, const char transb, const int n, const int k,
	const T& alpha, const T* A, const int lda, const T* B, const int ldb,
	const T& beta, T* C, const int ldc);
```
C の uplo 三角部分のみ := alpha\*op(A)\*op(B) + beta\*C．C: n×n，op(A): n×k，
op(B): k×n．reference BLAS の GEMMTR (LAPACK 3.12.1 で追加) に対応 (R1)．

### tsymm — 対称行列積
```cpp
template <typename T>
void tsymm(const char side, const char uplo, const int m, const int n,
	const T& alpha, const T* A, const int lda, const T* B, const int ldb,
	const T& beta, T* C, const int ldc);
```
C := alpha\*A\*B + beta\*C (side='L') / C := alpha\*B\*A + beta\*C (side='R')．
A は対称 (side='L' なら m×m，'R' なら n×n) で uplo の三角のみ参照．
実装は対称三角格納を full dense に展開 (copy のみ) して tgemm に帰着 (R1)．

### tsyrk — 対称 rank-k 更新
```cpp
template <typename T>
void tsyrk(const char uplo, const char trans, const int n, const int k,
	const T& alpha, const T* A, const int lda, const T& beta, T* C, const int ldc);
```
C := alpha\*op(A)\*op(A)^T + beta\*C の uplo 三角部分のみ更新．
trans='N': op(A) = A (n×k)，'T': op(A) = A^T (A は k×n) (R1)．

### tsyr2k — 対称 rank-2k 更新
```cpp
template <typename T>
void tsyr2k(const char uplo, const char trans, const int n, const int k,
	const T& alpha, const T* A, const int lda, const T* B, const int ldb,
	const T& beta, T* C, const int ldc);
```
C := alpha\*op(A)\*op(B)^T + alpha\*op(B)\*op(A)^T + beta\*C の uplo 三角のみ更新 (R1)．

### ttrmm — 三角行列積
```cpp
template <typename T>
void ttrmm(const char side, const char uplo, const char transa, const char diag,
	const int m, const int n, const T& alpha, const T* A, const int lda, T* B, const int ldb);
```
B := alpha\*op(A)\*B (side='L') / B := alpha\*B\*op(A) (side='R')．
A は三角 (side='L' なら m×m，'R' なら n×n)，B は m×n を in-place 更新．
実装は三角を 0 詰め dense に展開して tgemm に帰着 (exact 0 との積和なので
結果に影響せず，区間型でも幅は増えない)．`alpha == T(0)` のとき B := 0 (R1)．

### ttrsm — 三角行列方程式の求解
```cpp
template <typename T>
void ttrsm(const char side, const char uplo, const char transa, const char diag,
	const int m, const int n, const T& alpha, const T* A, const int lda, T* B, const int ldb);
```
op(A)\*X = alpha\*B (side='L') / X\*op(A) = alpha\*B (side='R') を解き，
B を解 X で上書きする．特異性検査は行わない (diag='N' の対角は正則を仮定)．
side='L' は列ごとに独立に solve (OpenMP 並列)，side='R' は列 sweep (逐次)．
要件: R1+R2 (`/`)．
