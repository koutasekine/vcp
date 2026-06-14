# tblas 関数リファレンス

`tblas/tblas.hpp` が提供する全 37 routine の仕様．設計方針と要素型 T の要件
(R1–R4) は `tblas/README_tblas.md` を参照．

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
- **丸めモードは扱わない**．丸めモード指定が必要なら `vblas/rdblas` を使う．
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
