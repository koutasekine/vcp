# tlapack 関数リファレンス

`tlapack/tlapack.hpp` が提供する全 routine の仕様．設計方針と要素型 T の要件
(R1–R6) は `tlapack/README_tlapack.md` を参照．アルゴリズムの数学的内容は
reference LAPACK 3.12.1 の同名 routine (`t` → `d`) に準ずる．

## 共通仕様

- 全関数は `template <typename T>` の関数 template (`tlamch<T>` 以外は引数から
  T が推論される)．
- **行列は column-major**: `A(i,j)` は `A[i + lda * j]`，`lda >= max(1, 行数)`．
- **返り値は LAPACK の INFO** (int)．`0`: 成功，`> 0`: LAPACK と同じ意味
  (位置・個数は **1-based**)．LAPACK で `INFO < 0` になる引数 error は返り値では
  なく `std::invalid_argument` を投げる (例外無効時は `abort`)．
- **WORK / LWORK 引数はない**．作業領域は内部で `std::vector<T>` により
  最適 size を確保する (workspace query 不要)．
- **index の規約** (rdlapack と同一):

  | 対象 | 規約 |
  |---|---|
  | `tgetrf` / `tgesv` / `tgbtrf` / `tgbsv` 系の `ipiv`，`tlaswp` | **0-based** |
  | `tsytf2` / `tsytrf` / `tsytrs` / `tsysv` の `ipiv` | **1-based 符号付き** (LAPACK と同じ，負値は 2×2 block) |
  | `tgebal` / `tgebak` / `tgehrd` / `torghr` / `thseqr` / `tlahqr` の `ilo` / `ihi` | **1-based** |
  | 返り値 INFO 中の位置・個数 | **1-based** (LAPACK と同じ) |
- **丸めモードは扱わない**．丸めモード指定が必要なら `vlapack/rdlapack` を使う．
- option 文字は大文字小文字どちらも可．
- scalar 引数は `const T&`，出力 scalar は `T&` 渡し．
- `namespace tlapack_detail` 内の関数 (`tabs` `tsqrt` `tisnan` `f_sign` `pow2i`
  `ilaenv` など) は内部実装用であり，この文書では公開 API のみ記載する．

### reference LAPACK との実装上の相違 (結果の意味は同じ)

- `tlassq` は Blue's algorithm でなく古典的な scale/ssq 更新．
- `tlanv2` の scaling 閾値は 2^-485 でなく `sqrt(safmin/eps)`．
- `tbdsqr` は特異値のみの場合も dqds (dlasq1) でなく QR 反復．
- `tgesvd` は基本経路 (直接二重対角化) のみ (m≫n / n≫m 用の QR/LQ 前処理なし)．
- `f_sign(a, b)` は b = -0.0 を正と扱う (Householder の符号選択にのみ影響し，
  結果は同値な分解になる)．
- fma を使わないため，double で rdlapack / reference と bit 単位では一致しない．

---

## 機械定数

### tlamch — machine parameter
```cpp
template <typename T>
T tlamch(const char cmach);   // 例: tlamch<kv::dd>('E')
```
`std::numeric_limits<T>` から machine parameter を構成して返す (要件 R5)．
T は推論できないため明示する．

| cmach | 値 | 意味 |
|---|---|---|
| `'E'` | epsilon()/2 | 相対 machine epsilon (丸め単位，double で 2^-53) |
| `'S'` | safe minimum | 1/sfmin が overflow しない最小の正数 |
| `'B'` | radix | 基数 (通常 2) |
| `'P'` | epsilon() | eps \* base (double で 2^-52) |
| `'N'` | digits | 仮数部 bit 数 |
| `'R'` | 1 | rounding ('R' = 丸めあり) |
| `'M'` | min_exponent | emin |
| `'U'` | min() | underflow threshold |
| `'L'` | max_exponent | emax |
| `'O'` | max() | overflow threshold |

---

## 連立一次方程式 (一般行列, LU)

### tgesv — driver: A\*X = B
```cpp
template <typename T>
int tgesv(const int n, const int nrhs, T* A, const int lda, int* ipiv, T* B, const int ldb);
```
部分 pivot 付き LU 分解で A\*X = B を解く．
出力: A は L, U で上書き，`ipiv[i]` (0-based) は行 i と交換した行，
B は解 X で上書き．INFO > 0: U(i,i) = 0 で特異 (i は 1-based)．

### tgetrf / tgetrf2 / tgetf2 — LU 分解
```cpp
template <typename T>
int tgetrf(const int m, const int n, T* A, const int lda, int* ipiv);   // blocked
int tgetrf2(const int m, const int n, T* A, const int lda, int* ipiv);  // 再帰版
int tgetf2(const int m, const int n, T* A, const int lda, int* ipiv);   // unblocked
```
A = P\*L\*U (部分 pivot，L は単位下三角)．A を L (対角より下) と U (上三角) で
上書きする．`ipiv` は 0-based (長さ min(m,n))．INFO > 0: U(i,i) = 0．

### tgetrs — LU 分解による求解
```cpp
template <typename T>
int tgetrs(const char trans, const int n, const int nrhs,
	const T* A, const int lda, const int* ipiv, T* B, const int ldb);
```
tgetrf の結果を使い op(A)\*X = B (trans = 'N' / 'T','C') を解く．B を X で上書き．

### tgetri — LU 分解からの逆行列
```cpp
template <typename T>
int tgetri(const int n, T* A, const int lda, const int* ipiv);
```
tgetrf の出力 A, ipiv から inv(A) を計算し A を上書きする．INFO > 0: 特異．

---

## 三角行列

### ttrtrs — 三角方程式の求解 (特異性検査付き)
```cpp
template <typename T>
int ttrtrs(const char uplo, const char trans, const char diag, const int n, const int nrhs,
	const T* A, const int lda, T* B, const int ldb);
```
op(A)\*X = B を解き B を上書き．INFO > 0: A(i,i) = 0 で特異 (解かずに返る)．
(tblas の `ttrsm` と異なり対角 0 を検査する．)

### ttrtri / ttrti2 — 三角行列の逆行列
```cpp
template <typename T>
int ttrtri(const char uplo, const char diag, const int n, T* A, const int lda);  // blocked
int ttrti2(const char uplo, const char diag, const int n, T* A, const int lda);  // unblocked
```
inv(A) を計算し A を上書き．INFO > 0: A(i,i) = 0 で特異．

### tlauum / tlauu2 — 三角因子の積
```cpp
template <typename T>
int tlauum(const char uplo, const int n, T* A, const int lda);  // blocked
int tlauu2(const char uplo, const int n, T* A, const int lda);  // unblocked
```
uplo='U': A の上三角 U に対し U\*U^T，'L': L^T\*L を計算し，uplo の三角に上書き
(tpotri の内部で使用)．

---

## 連立一次方程式 (対称正定値, Cholesky)

### tposv — driver: A\*X = B (A 対称正定値)
```cpp
template <typename T>
int tposv(const char uplo, const int n, const int nrhs, T* A, const int lda, T* B, const int ldb);
```
Cholesky 分解 A = U^T\*U ('U') / L\*L^T ('L') で解く．A は因子で，B は解で上書き．
INFO > 0: 第 i 主小行列が正定値でない (i は 1-based)．

### tpotrf / tpotrf2 / tpotf2 — Cholesky 分解
```cpp
template <typename T>
int tpotrf(const char uplo, const int n, T* A, const int lda);   // blocked
int tpotrf2(const char uplo, const int n, T* A, const int lda);  // 再帰版
int tpotf2(const char uplo, const int n, T* A, const int lda);   // unblocked
```
A = U^T\*U / L\*L^T．uplo の三角を因子で上書き (反対側は変更しない)．
INFO > 0: 正定値でない (対角が 0 以下，または NaN — 要件 R6 の `tisnan` で検出)．

### tpotrs — Cholesky 分解による求解
```cpp
template <typename T>
int tpotrs(const char uplo, const int n, const int nrhs,
	const T* A, const int lda, T* B, const int ldb);
```
tpotrf の因子で A\*X = B を解き B を上書き．

### tpotri — Cholesky 分解からの逆行列
```cpp
template <typename T>
int tpotri(const char uplo, const int n, T* A, const int lda);
```
tpotrf の因子から inv(A) の uplo 三角を計算し上書き．INFO > 0: 因子が特異．

---

## 連立一次方程式 (対称不定値, Bunch-Kaufman)

### tsysv — driver: A\*X = B (A 対称)
```cpp
template <typename T>
int tsysv(const char uplo, const int n, const int nrhs, T* A, const int lda,
	int* ipiv, T* B, const int ldb);
```
Bunch-Kaufman 分解 A = U\*D\*U^T / L\*D\*L^T で解く．INFO > 0: D(i,i) = 0．

### tsytrf / tsytf2 / tlasyf — Bunch-Kaufman 分解
```cpp
template <typename T>
int tsytrf(const char uplo, const int n, T* A, const int lda, int* ipiv);  // blocked
int tsytf2(const char uplo, const int n, T* A, const int lda, int* ipiv);  // unblocked
int tlasyf(const char uplo, const int n, const int nb, int& kb,
	T* A, const int lda, int* ipiv, T* W, const int ldw);  // 部分分解 (tsytrf の内部)
```
A = U\*D\*U^T / L\*D\*L^T (D は 1×1 / 2×2 block 対角)．
`ipiv` は **1-based 符号付き** (LAPACK と同じ): `ipiv[k] > 0` なら 1×1 block で
行 k+1 (1-based) と `ipiv[k]` を交換，`ipiv[k] = ipiv[k+1] < 0` なら 2×2 block．
`tlasyf` は先頭 (uplo='L') / 末尾 ('U') の nb 列を分解し，実際に処理した列数を
`kb` に返す (W は n×nb の作業領域)．

### tsytrs — Bunch-Kaufman 分解による求解
```cpp
template <typename T>
int tsytrs(const char uplo, const int n, const int nrhs,
	const T* A, const int lda, const int* ipiv, T* B, const int ldb);
```
tsytrf の結果で A\*X = B を解き B を上書き．

---

## 連立一次方程式 (帯行列)

band 格納: 一般 band は `AB[kl + ku + i - j + ldab*j] = A(i,j)`，
`ldab >= 2*kl + ku + 1` (LU の fill-in 用に kl 行余分に必要)．
対称 band は uplo='U' なら `AB[kd + i - j + ldab*j]` (i <= j)，
'L' なら `AB[i - j + ldab*j]` (i >= j)，`ldab >= kd + 1`．

### tgbsv — driver: 帯行列の A\*X = B
```cpp
template <typename T>
int tgbsv(const int n, const int kl, const int ku, const int nrhs,
	T* AB, const int ldab, int* ipiv, T* B, const int ldb);
```
部分 pivot 付き band LU で解く．`ipiv` は 0-based．INFO > 0: U(i,i) = 0．

### tgbtrf / tgbtf2 / tgbtrs — band LU 分解と求解
```cpp
template <typename T>
int tgbtrf(const int m, const int n, const int kl, const int ku,
	T* AB, const int ldab, int* ipiv);   // blocked
int tgbtf2(const int m, const int n, const int kl, const int ku,
	T* AB, const int ldab, int* ipiv);   // unblocked
int tgbtrs(const char trans, const int n, const int kl, const int ku, const int nrhs,
	const T* AB, const int ldab, const int* ipiv, T* B, const int ldb);
```

### tpbsv — driver: 対称正定値帯行列の A\*X = B
```cpp
template <typename T>
int tpbsv(const char uplo, const int n, const int kd, const int nrhs,
	T* AB, const int ldab, T* B, const int ldb);
```
band Cholesky で解く．INFO > 0: 正定値でない．

### tpbtrf / tpbtf2 / tpbtrs — band Cholesky 分解と求解
```cpp
template <typename T>
int tpbtrf(const char uplo, const int n, const int kd, T* AB, const int ldab);  // blocked
int tpbtf2(const char uplo, const int n, const int kd, T* AB, const int ldab);  // unblocked
int tpbtrs(const char uplo, const int n, const int kd, const int nrhs,
	const T* AB, const int ldab, T* B, const int ldb);
```

---

## QR / LQ / QL 分解と最小二乗

Householder 表現 (LAPACK と同じ): Q = H(0)\*H(1)\*...\*H(k-1)，
H(i) = I - tau[i]\*v\*v^T．QR では v は A の対角より下に格納され v の先頭成分 1 は
格納されない．R は A の上三角に残る．

### tgels — driver: 最小二乗 / 最小 norm 解
```cpp
template <typename T>
int tgels(const char trans, const int m, const int n, const int nrhs,
	T* A, const int lda, T* B, const int ldb);
```
A が full rank と仮定し，trans='N' かつ m >= n なら min ||B - A\*X|| (QR)，
m < n なら A\*X = B の最小 norm 解 (LQ)．trans='T' は A^T について同様．
B は max(m,n)×nrhs (`ldb >= max(1,m,n)`) で，解 X は先頭 n 行 (trans='N')
または m 行 ('T') に上書きされる．優決定の場合，残りの行には residual の
情報が残る (各列の要素平方和が residual の 2 乗和)．
INFO > 0: 三角因子の対角 i が 0 (full rank でない)．

### tgeqrf / tgeqr2 — QR 分解
```cpp
template <typename T>
int tgeqrf(const int m, const int n, T* A, const int lda, T* tau);  // blocked
int tgeqr2(const int m, const int n, T* A, const int lda, T* tau);  // unblocked
```
A = Q\*R．A の上三角に R，対角より下に Householder vector，`tau` (長さ min(m,n))
に係数を格納する．

### tgelqf / tgelq2 — LQ 分解
```cpp
template <typename T>
int tgelqf(const int m, const int n, T* A, const int lda, T* tau);  // blocked
int tgelq2(const int m, const int n, T* A, const int lda, T* tau);  // unblocked
```
A = L\*Q．A の下三角に L，対角より右に Householder vector を格納する．

### torgqr / torg2r，torglq / torgl2，torgql / torg2l — 直交行列の生成
```cpp
template <typename T>
int torgqr(const int m, const int n, const int k, T* A, const int lda, const T* tau);
int torglq(const int m, const int n, const int k, T* A, const int lda, const T* tau);
int torgql(const int m, const int n, const int k, T* A, const int lda, const T* tau);
// torg2r / torgl2 / torg2l は各 unblocked 版 (同じ引数)
```
tgeqrf / tgelqf / (QL 分解) の出力から，Q の先頭 n 列 (orgqr)・先頭 m 行 (orglq)・
末尾 n 列 (orgql) を陽に生成して A に上書きする．k は反射の個数
(`m >= n >= k` (orgqr) 等)．

### tormqr / torm2r，tormlq / torml2，tormql / torm2l — 直交行列の適用
```cpp
template <typename T>
int tormqr(const char side, const char trans, const int m, const int n, const int k,
	const T* A, const int lda, const T* tau, T* C, const int ldc);
// tormlq / tormql，および torm2r / torml2 / torm2l (unblocked) も同じ引数
```
C := Q\*C ('L','N')，Q^T\*C ('L','T')，C\*Q ('R','N')，C\*Q^T ('R','T')．
Q は分解の出力 (A, tau) の暗黙表現のまま適用する (生成より速い)．C は m×n．

### tlarfg — Householder 反射の生成
```cpp
template <typename T>
void tlarfg(const int n, T& alpha, T* x, const int incx, T& tau);
```
H\*(alpha; x) = (beta; 0)，H = I - tau\*v\*v^T となる v, tau, beta を生成する．
入出力: `alpha` は beta で，`x` (長さ n-1) は v(1:n-1) で上書き (v(0) = 1 は暗黙)．

### tlarf / tlarf1f / tlarf1l — 反射の適用 (unblocked)
```cpp
template <typename T>
void tlarf(const char side, const int m, const int n, const T* v, const int incv,
	const T& tau, T* C, const int ldc);    // v は全成分を陽に持つ
void tlarf1f(...同じ引数...);              // v(0) = 1 を暗黙とする (先頭)
void tlarf1l(...同じ引数...);              // v(last) = 1 を暗黙とする (末尾)
```
C := H\*C (side='L') / C\*H (side='R')，H = I - tau\*v\*v^T．

### tlarft / tlarfb — block 反射
```cpp
template <typename T>
void tlarft(const char direct, const char storev, const int n, const int k,
	const T* V, const int ldv, const T* tau, T* Tf, const int ldt);
template <typename T>
void tlarfb(const char side, const char trans, const char direct, const char storev,
	const int m, const int n, const int k, const T* V, const int ldv,
	const T* Tf, const int ldt, T* C, const int ldc);
```
k 個の反射の積 H = I - V\*Tf\*V^T の三角因子 Tf (k×k) を作る (tlarft) /
C に適用する (tlarfb)．direct: 'F' (H(0)..H(k-1)) / 'B'，storev: 'C' (列格納) / 'R'．
**LAPACK で T と呼ばれる引数は template parameter と衝突するため `Tf` に改名**
している (機能は同じ)．

---

## 対称固有値問題

### tsyev — driver: 対称行列の全固有値・固有 vector
```cpp
template <typename T>
int tsyev(const char jobz, const char uplo, const int n, T* A, const int lda, T* w);
```
jobz='N': 固有値のみ (tsterf)，'V': 固有 vector も (tsteqr)．
出力: `w` (長さ n) に固有値が**昇順**，jobz='V' なら A に正規直交固有 vector
(列が w[j] に対応)．INFO > 0: 収束しなかった非対角要素の個数．

### tsygv — driver: 一般化対称固有値問題
```cpp
template <typename T>
int tsygv(const int itype, const char jobz, const char uplo, const int n,
	T* A, const int lda, T* B, const int ldb, T* w);
```
itype=1: A\*z = λ\*B\*z，2: A\*B\*z = λ\*z，3: B\*A\*z = λ\*z (B は対称正定値)．
出力: `w` に固有値昇順，jobz='V' なら A に固有 vector
(itype=1,2: z^T\*B\*z = I，itype=3: z^T\*inv(B)\*z = I に正規化)，
B は Cholesky 因子で上書き．INFO > 0: n 以下なら tsyev の失敗，
n より大なら B が正定値でない (INFO - n が主小行列の位置)．

### tsytrd / tsytd2 / tlatrd — 三重対角化
```cpp
template <typename T>
int tsytrd(const char uplo, const int n, T* A, const int lda, T* d, T* e, T* tau);  // blocked
int tsytd2(const char uplo, const int n, T* A, const int lda, T* d, T* e, T* tau);  // unblocked
void tlatrd(const char uplo, const int n, const int nb, T* A, const int lda,
	T* e, T* tau, T* W, const int ldw);  // 部分変換 (tsytrd の内部)
```
Q^T\*A\*Q = T (三重対角)．`d` (長さ n) に対角，`e` (長さ n-1) に副対角，
A に Householder vector，`tau` に係数を格納する．

### torgtr / tormtr — 三重対角化の Q の生成・適用
```cpp
template <typename T>
int torgtr(const char uplo, const int n, T* A, const int lda, const T* tau);
int tormtr(const char side, const char uplo, const char trans, const int m, const int n,
	const T* A, const int lda, const T* tau, T* C, const int ldc);
```
tsytrd の出力から Q を生成 / C に適用する (uplo は tsytrd に渡したものと同じ)．

### tsterf — 三重対角行列の固有値 (Pal-Walker-Kahan QL/QR)
```cpp
template <typename T>
int tsterf(const int n, T* d, T* e);
```
`d` (対角)，`e` (副対角) の対称三重対角行列の固有値を計算し，`d` に昇順で返す
(`e` は破壊される)．INFO > 0: 収束しなかった要素数．

### tsteqr — 三重対角行列の固有値・固有 vector (implicit QL/QR)
```cpp
template <typename T>
int tsteqr(const char compz, const int n, T* d, T* e, T* Z, const int ldz);
```
compz='N': 固有値のみ，'V': Z に三重対角化の Q を渡して元の対称行列の
固有 vector を計算，'I': Z を単位行列に初期化して三重対角行列自身の
固有 vector を計算．`d` に固有値昇順，Z (n×n) に固有 vector．

### tsygst / tsygs2 — 一般化問題の標準形への変換
```cpp
template <typename T>
int tsygst(const int itype, const char uplo, const int n, T* A, const int lda,
	const T* B, const int ldb);  // blocked
int tsygs2(...同じ引数...);       // unblocked
```
tpotrf 済みの B の因子を使い，一般化固有値問題を標準形に変換する
(itype=1: inv(U^T)\*A\*inv(U) 等)．A を上書き．

---

## 非対称固有値問題

### tgeev — driver: 一般行列の固有値・固有 vector
```cpp
template <typename T>
int tgeev(const char jobvl, const char jobvr, const int n, T* A, const int lda,
	T* wr, T* wi, T* VL, const int ldvl, T* VR, const int ldvr);
```
固有値 `wr[j] + i*wi[j]` と，jobvl/jobvr='V' なら左/右固有 vector を計算する
(`'N'` なら VL/VR は参照されない)．複素共役対は隣接し wi > 0 が先．
実固有値の vector は VR の列 j，複素対は v = VR(:,j) ± i\*VR(:,j+1)
(VL も同様)．vector は Euclid norm 1 に正規化され，複素対は絶対値最大の成分が
実になるよう回転される (reference DGEEV と同じ)．A は破壊される．
INFO > 0: QR 反復が収束せず (wr/wi の INFO+1 .. n 番目 (1-based) は正しい)．

### tgebal / tgebak — balancing とその逆変換
```cpp
template <typename T>
int tgebal(const char job, const int n, T* A, const int lda, int& ilo, int& ihi, T* scale);
int tgebak(const char job, const char side, const int n, const int ilo, const int ihi,
	const T* scale, const int m, T* V, const int ldv);
```
job='N': 何もしない，'P': 置換のみ，'S': scaling のみ，'B': 両方．
`ilo`/`ihi` は **1-based** で，A(i,j) = 0 (i > j, j < ilo or i > ihi) の形に
整える．`scale` (長さ n) に置換と scaling の情報を格納．tgebak は固有 vector
V (n×m) に逆変換を適用する (side='L'/'R' は左/右固有 vector)．

### tgehrd / tgehd2 / tlahr2 — Hessenberg 化
```cpp
template <typename T>
int tgehrd(const int n, const int ilo, const int ihi, T* A, const int lda, T* tau);  // blocked
int tgehd2(const int n, const int ilo, const int ihi, T* A, const int lda, T* tau);  // unblocked
void tlahr2(const int n, const int k, const int nb, T* A, const int lda,
	T* tau, T* Tf, const int ldt, T* Y, const int ldy);  // 部分変換 (tgehrd の内部)
```
Q^T\*A\*Q = H (上 Hessenberg)．ilo/ihi は 1-based (tgebal の出力，
balancing なしなら ilo=1, ihi=n)．

### torghr — Hessenberg 化の Q の生成
```cpp
template <typename T>
int torghr(const int n, const int ilo, const int ihi, T* A, const int lda, const T* tau);
```

### thseqr / tlahqr — Hessenberg 行列の固有値 (実 Schur 形)
```cpp
template <typename T>
int thseqr(const char job, const char compz, const int n, const int ilo, const int ihi,
	T* H, const int ldh, T* wr, T* wi, T* Z, const int ldz);
template <typename T>
int tlahqr(const bool wantt, const bool wantz, const int n, const int ilo, const int ihi,
	T* H, const int ldh, T* wr, T* wi, const int iloz, const int ihiz, T* Z, const int ldz);
```
double shift QR 反復．job='E': 固有値のみ，'S': H を実 Schur 形 T に上書き．
compz='N': Z なし，'I': Z を単位行列から (Schur vector)，'V': 与えた Z に累積．
INFO > 0: 収束せず．`tlahqr` は thseqr の内部で使う unblocked 版．

### ttrevc — 実 Schur 形の固有 vector
```cpp
template <typename T>
int ttrevc(const char side, const char howmny, const int n, const T* Tf, const int ldt,
	T* VL, const int ldvl, T* VR, const int ldvr, const int mm, int& m);
```
実 Schur 形 Tf の右/左固有 vector (side='R'/'L'/'B') を計算する．
howmny='A': Schur 形の固有 vector，'B': VL/VR に与えた Schur vector で
back-transform して元行列の固有 vector (**SELECT 指定 ('S') は未対応**)．
mm は VL/VR の列数，m に実際に使った列数が返る．

### tlaln2 — 小さい線形系 (1×1 / 2×2，複素 shift 付き) の求解
```cpp
template <typename T>
int tlaln2(const bool ltrans, const int na, const int nw, const T& smin, const T& ca,
	const T* A, const int lda, const T& d1, const T& d2, const T* B, const int ldb,
	const T& wr, const T& wi, T* X, const int ldx, T& scale, T& xnorm);
```
(ca\*A - w\*D)\*X = s\*B (na = 1,2，nw = 1,2) を overflow を避けて解く
(ttrevc の内部用)．返り値 1: 摂動して解いた．

---

## 特異値分解

### tgesvd — driver: 特異値分解 A = U\*S\*V^T
```cpp
template <typename T>
int tgesvd(const char jobu, const char jobvt, const int m, const int n,
	T* A, const int lda, T* s, T* U, const int ldu, T* VT, const int ldvt);
```
`s` (長さ min(m,n)) に特異値を**降順**で返す．
jobu = 'A': U (m×m) 全列 / 'S': 先頭 min(m,n) 列 / 'O': U の先頭列で A を上書き /
'N': 計算しない．jobvt も同様に V^T の行について 'A' (n×n) / 'S' / 'O' / 'N'
(jobu と jobvt の両方を 'O' にはできない)．
INFO > 0: tbdsqr が収束しなかった副対角の個数．

### tgebrd / tgebd2 / tlabrd — 二重対角化
```cpp
template <typename T>
int tgebrd(const int m, const int n, T* A, const int lda,
	T* d, T* e, T* tauq, T* taup);  // blocked
int tgebd2(...同じ引数...);          // unblocked
void tlabrd(const int m, const int n, const int nb, T* A, const int lda, T* d, T* e,
	T* tauq, T* taup, T* X, const int ldx, T* Y, const int ldy);  // 部分変換
```
Q^T\*A\*P = B (二重対角)．`d` (長さ min(m,n)) に対角，`e` に副対角
(m >= n なら上二重対角)，A に Q と P の Householder vector，
`tauq`/`taup` に係数を格納する．

### torgbr / tormbr — 二重対角化の Q, P^T の生成・適用
```cpp
template <typename T>
int torgbr(const char vect, const int m, const int n, const int k,
	T* A, const int lda, const T* tau);
int tormbr(const char vect, const char side, const char trans, const int m, const int n,
	const int k, const T* A, const int lda, const T* tau, T* C, const int ldc);
```
vect='Q': Q を生成/適用，'P': P^T を生成/適用 (tau は tauq / taup)．

### tbdsqr — 二重対角行列の SVD (QR 反復)
```cpp
template <typename T>
int tbdsqr(const char uplo, const int n, const int ncvt, const int nru, const int ncc,
	T* d, T* e, T* VT, const int ldvt, T* U, const int ldu, T* C, const int ldc);
```
二重対角行列 (uplo='U' 上 / 'L' 下) の特異値を `d` に降順で返し，
ncvt > 0 なら VT (n×ncvt) に V^T 行列を，nru > 0 なら U (nru×n) を，
ncc > 0 なら C (n×ncc) を右/左から回転で更新する．
特異値のみ (ncvt = nru = ncc = 0) の場合も dqds でなく QR 反復で計算する．
INFO > 0: 収束しなかった副対角の個数．

---

## 補助 routine

### tlaset — 行列の定数初期化
```cpp
template <typename T>
void tlaset(const char uplo, const int m, const int n, const T& alpha, const T& beta,
	T* A, const int lda);
```
非対角を alpha，対角を beta に設定する (uplo='U'/'L' は三角部分のみ，
それ以外は全体)．

### tlacpy — 行列の copy
```cpp
template <typename T>
void tlacpy(const char uplo, const int m, const int n, const T* A, const int lda,
	T* B, const int ldb);
```
B := A (uplo='U'/'L' は三角部分のみ，それ以外は全体)．

### tlaswp — 行交換の適用
```cpp
template <typename T>
void tlaswp(const int n, T* A, const int lda, const int k1, const int k2,
	const int* ipiv, const int incx);
```
`i = k1, ..., k2-1` (**0-based 半開区間**) について行 i と行 `ipiv[i]` を交換する
(incx < 0 なら逆順)．tgetrf 系の 0-based ipiv に対応．

### tlasrt — sort
```cpp
template <typename T>
void tlasrt(const char id, const int n, T* d);
```
d を昇順 ('I') / 降順 ('D') に sort する (比較のみ，要件 R3)．

### ilatlr / ilatlc — 非零部分の大きさ
```cpp
template <typename T>
int ilatlr(const int m, const int n, const T* A, const int lda);  // 行
int ilatlc(const int m, const int n, const T* A, const int lda);  // 列
```
非零要素を含む最後の行/列までの個数 (= 1-based index，全零なら 0) を返す．

### tlassq — scaled sum of squares
```cpp
template <typename T>
void tlassq(const int n, const T* x, const int incx, T& scale, T& sumsq);
```
scale^2\*sumsq += sum x_i^2 となるよう scale, sumsq を更新する
(overflow/underflow 回避)．古典的 scale/ssq 更新で実装
(reference 3.12 の Blue's algorithm と結果が末位で異なりうる)．

### tlapy2 / tlapy3 — 安全な平方和の平方根
```cpp
template <typename T>
T tlapy2(const T& x, const T& y);                // sqrt(x^2 + y^2)
T tlapy3(const T& x, const T& y, const T& z);    // sqrt(x^2 + y^2 + z^2)
```
overflow を避けて計算する．

### tladiv — 複素除算
```cpp
template <typename T>
void tladiv(const T& a, const T& b, const T& c, const T& d, T& p, T& q);
```
(a + bi)/(c + di) = p + qi を overflow を避けて計算する (Baudin-Smith 方式)．

### tlartg — Givens 回転の生成 (LAPACK 版)
```cpp
template <typename T>
void tlartg(const T& f, const T& g, T& cs, T& sn, T& r);
```
[cs sn; -sn cs]\*[f; g] = [r; 0]．`tblas` の trotg と異なり z を返さない
LAPACK (3.10 以降) 規約: g = 0 なら cs = 1, r = f (f = g = 0 を含む)，
f = 0 なら cs = 0, sn = sign(g), r = |g|，それ以外は r の符号を f に合わせる．

### tlae2 / tlaev2 — 2×2 対称行列の固有値
```cpp
template <typename T>
void tlae2(const T& a, const T& b, const T& c, T& rt1, T& rt2);
void tlaev2(const T& a, const T& b, const T& c, T& rt1, T& rt2, T& cs1, T& sn1);
```
[a b; b c] の固有値 rt1 (絶対値大), rt2 を計算する (tlaev2 は固有 vector
(cs1, sn1) も)．

### tlas2 / tlasv2 — 2×2 三角行列の特異値
```cpp
template <typename T>
void tlas2(const T& f, const T& g, const T& h, T& ssmin, T& ssmax);
void tlasv2(const T& f, const T& g, const T& h, T& ssmin, T& ssmax,
	T& snr, T& csr, T& snl, T& csl);
```
[f g; 0 h] の特異値 ssmin, ssmax (tlasv2 は左右の回転も) を計算する．

### tlanv2 — 2×2 実 Schur 標準形
```cpp
template <typename T>
void tlanv2(T& a, T& b, T& c, T& d, T& rt1r, T& rt1i, T& rt2r, T& rt2i, T& cs, T& sn);
```
[a b; c d] を実 Schur 標準形に変換する回転 (cs, sn) と固有値を計算する
(in-place)．scaling 閾値は `sqrt(safmin/eps)` (reference の 2^-485 と僅かに異なる)．

### tlascl — 安全な scalar 倍
```cpp
template <typename T>
void tlascl(const char type, const int kl, const int ku, const T& cfrom, const T& cto,
	const int m, const int n, T* A, const int lda);
```
A := A\*(cto/cfrom) を overflow/underflow を起こさず段階的に計算する．
type: 'G' 全体 / 'L' 下三角 / 'U' 上三角 / 'H' 上 Hessenberg /
'B','Q','Z' band 格納 (kl, ku は band 幅)．cfrom = 0 や NaN は invalid_argument．

### tlasr — 平面回転列の適用
```cpp
template <typename T>
void tlasr(const char side, const char pivot, const char direct, const int m, const int n,
	const T* c, const T* s, T* A, const int lda);
```
回転列 P を A に適用する (side='L': A := P\*A，'R': A := A\*P^T)．
pivot: 'V' 隣接 / 'T' 先頭固定 / 'B' 末尾固定，direct: 'F' / 'B'．

### tlange / tlansy / tlanst / tlangb / tlansb — 行列 norm
```cpp
template <typename T>
T tlange(const char norm, const int m, const int n, const T* A, const int lda);
T tlansy(const char norm, const char uplo, const int n, const T* A, const int lda);
T tlanst(const char norm, const int n, const T* d, const T* e);       // 対称三重対角
T tlangb(const char norm, const int n, const int kl, const int ku,
	const T* AB, const int ldab);                                      // band
T tlansb(const char norm, const char uplo, const int n, const int k,
	const T* AB, const int ldab);                                      // 対称 band
```
norm = 'M': max |a_ij|，'1'/'O': 1-norm，'I': inf-norm，'F'/'E': Frobenius norm．
max 系は NaN を伝播する (要件 R6)．
