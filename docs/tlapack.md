# tlapack — テンプレート版 LAPACK (T 型要件と実装)

`tlapack` は `vcp/vlapack/rdlapack` (reference LAPACK 3.12.1 由来) と同じ routine
一覧を template 化した header-only library です．本文書は要素型 `T` に要求する
演算 (concept 相当) と実装・テストの要点を定めます．BLAS 部分は `vcp/tblas/` を
呼び出し，**丸めモードは一切意識しません** (rounding_mode 引数なし)．
**全関数の引数仕様は このドキュメントの「関数リファレンス」節 を参照．**

```cpp
#include <vcp/tlapack/tlapack.hpp>
#include <kv/dd.hpp>

// A*X = B を kv::dd で解く (column-major, ipiv は 0-based)
int info = vcp::tgesv(n, nrhs, A, lda, ipiv, B, ldb);   // T は引数から推論される
int info2 = vcp::tsyev('V', 'U', n, A, lda, w);          // 対称固有値 (kv::dd のまま高精度)

// using namespace vcp; を使えば vcp:: を省略できる
using namespace vcp;
int info3 = tgesv(n, nrhs, A, lda, ipiv, B, ldb);
```

## file 構成

| file | 内容 |
|---|---|
| `tlapack.hpp` | umbrella header (これだけ include すればよい) |
| `tlapack_common.hpp` | 共通基盤 (`tlamch<T>` `tisnan` `f_sign` `pow2i` `ilaenv`) |
| `tlapack_aux.hpp` | 補助 routine (`tlassq` `tlartg` `tlarfg` `tlarf*` `tlange` 系など) |
| `tlapack_lu.hpp` `_chol` `_band` `_qr` `_sy` `_eig` `_geev` `_svd` | rdlapack と同じ分割 |
| `tlapack_double.hpp` | `T=double` の明示的特殊化．`dblas_dlapack.hpp` 経由で double LAPACK/BLAS を呼ぶ |

## 提供予定 routine (rdlapack と同一集合，rd → t に読み替え)

| 分類 | driver | 計算 routine |
|---|---|---|
| 連立一次方程式 (一般) | `tgesv` | `tgetf2` `tgetrf2` `tgetrf` `tgetrs` `tgetri` |
| 連立一次方程式 (対称正定値) | `tposv` | `tpotf2` `tpotrf2` `tpotrf` `tpotrs` `tpotri` |
| 連立一次方程式 (対称不定値) | `tsysv` | `tsytf2` `tlasyf` `tsytrf` `tsytrs` |
| 連立一次方程式 (帯) | `tgbsv` `tpbsv` | `tgbtf2` `tgbtrf` `tgbtrs` `tpbtf2` `tpbtrf` `tpbtrs` |
| 三角行列 | — | `ttrtrs` `ttrti2` `ttrtri` `tlauu2` `tlauum` |
| 最小二乗 | `tgels` | QR/LQ 一式 |
| QR / LQ / QL | — | `tgeqr2` `tgeqrf` `tgelq2` `tgelqf` `torg2r` `torgqr` `torgl2` `torglq` `torg2l` `torgql` `torm2r` `tormqr` `torml2` `tormlq` `torm2l` `tormql` |
| 対称固有値 | `tsyev` | `tsytd2` `tlatrd` `tsytrd` `torgtr` `tormtr` `tsterf` `tsteqr` |
| 一般化対称固有値 | `tsygv` | `tsygs2` `tsygst` |
| 非対称固有値 | `tgeev` | `tgebal` `tgebak` `tgehd2` `tlahr2` `tgehrd` `torghr` `tlahqr` `thseqr` `tlaln2` `ttrevc` |
| 特異値分解 | `tgesvd` | `tgebd2` `tlabrd` `tgebrd` `torgbr` `tormbr` `tbdsqr` |
| 補助 | — | `tlassq` `tlapy2` `tlapy3` `tladiv` `tlartg` `tlarfg` `tlarf` `tlarf1f` `tlarf1l` `tlarft` `tlarfb` `tlascl` `tlasr` `tlae2` `tlaev2` `tlas2` `tlasv2` `tlanv2` `tlange` `tlangb` `tlansy` `tlansb` `tlanst` `tlamch` `tlaset` `tlacpy` `tlaswp` `tlasrt` `ilatlr` `ilatlc` |

引数規約は rdlapack と同じ (reference LAPACK と同順で rounding_mode なし，
INFO は返り値，引数 error は `std::invalid_argument`，WORK/LWORK なし，
ipiv / ilo / ihi の 0-based / 1-based 規約も rdlapack と同一)．
scalar は tblas と同様 `const T&` 渡しとします．

## T=double の LAPACK/BLAS 委譲特殊化

`T=double` については `vcp/tlapack/tlapack_double.hpp` を include すると，
`dblas_dlapack.hpp` に宣言済みの double LAPACK routine を使って明示的特殊化します．
`tlapack_double.hpp` は `vcp/tblas/tblas_double.hpp` も include するため，
特殊化されていない補助 routine が内部で BLAS を呼ぶ場合も `tgemm<double>` などは
double BLAS へ委譲されます．

```cpp
#include <vcp/tlapack/tlapack_double.hpp>

int info = vcp::tgesv<double>(n, nrhs, A, lda, ipiv, B, ldb);
```

特殊化の範囲は次の方針です．

- `dgesv_`，`dgetrf_`，`dgetrs_`，`dgetri_`，`dpotrf_`，`dsyev_`，
  `dgesvd_`，`dgeev_` など，`dblas_dlapack.hpp` に既にある routine は
  `USE_VCP_LAPACK` の有無に関係なく特殊化します．`USE_VCP_LAPACK` がある場合は
  `dblas_dlapack.hpp` 経由で `vcp/vlapack/dlapack.hpp` の wrapper を呼び，
  ない場合は外部 LAPACK (MKL / OpenBLAS / reference LAPACK) の Fortran symbol を
  呼びます．
- `dgeqrf_`，`dsyev_`，`dgesvd_`，`dgeev_` など LAPACK では `WORK/LWORK` を持つ
  routine は，特殊化側で workspace query (`LWORK = -1`) を行い，
  `std::vector<double>` で作業領域を確保してから本計算を呼びます．
- `dlae2_`，`dlaev2_`，`dlas2_`，`dlasv2_`，`dlanv2_`，`dlasr_`，
  `dlarf_`，`dlarft_`，`dlarfb_`，`dlasyf_`，`dlatrd_`，`dlahr2_`，
  `dlabrd_`，`dlahqr_`，`dlaln2_` など，`dblas_dlapack.hpp` にない補助 routine は
  `#ifndef USE_VCP_LAPACK` の場合だけ `extern "C"` 宣言して特殊化します．
  `USE_VCP_LAPACK` がある場合，`vcp/vlapack/dlapack.hpp` はこれらの補助 symbol を
  提供しないため，特殊化せず既存の template<double> 実装にフォールバックします．
- `tgetrf2<double>` / `tpotrf2<double>` も同様に，外部 LAPACK に symbol がある場合だけ
  (`USE_VCP_LAPACK` なし) 特殊化し，`USE_VCP_LAPACK` 時は既存 template 実装を使います．
- `tlapack_double.hpp` は `t...<double>` を最初に使う前に include してください．
  明示的特殊化は，対象 template が同じ翻訳単位で既に暗黙実体化された後には
  追加できません．

### double 特殊化と ipiv

公開 API の `ipiv` 規約は `tlapack.hpp` の template 実装と同じです．
ただし外部 LAPACK/`dlapack.hpp` の Fortran ABI は LU 系 pivot を 1-based で扱うため，
`tlapack_double.hpp` の特殊化内部で相互変換します．

| 対象 | `tlapack` 公開 API | double 特殊化内部の処理 |
|---|---|---|
| `tgetf2` / `tgetrf` / `tgesv` / `tgbtf2` / `tgbtrf` / `tgbsv` の `ipiv` | 0-based | LAPACK 呼び出し後に 1-based → 0-based へ変換して返す |
| `tgetrs` / `tgetri` / `tgbtrs` に入力する `ipiv` | 0-based | LAPACK 呼び出し前に一時配列で 0-based → 1-based へ変換する |
| `tlaswp` の `k1` / `k2` / `ipiv` | `k1`,`k2` は 0-based 半開区間，`ipiv` は 0-based | LAPACK の `dlaswp_` 用に `k1+1`,`k2` と 1-based pivot へ変換する |
| `tsytf2` / `tsytrf` / `tsytrs` / `tsysv` の `ipiv` | 1-based 符号付き (LAPACK と同じ) | 変換せず pass-through |

このため，利用側は `T=double` でも `kv::dd` 等でも同じ `ipiv` 規約で呼べます．

## T 型要件

**tblas の R1–R4 はすべてそのまま必要です** (LAPACK は BLAS の全演算に加えて
scalar 演算でも abs / sqrt / 順序比較 / 除算を多用するため)．

- R1: 値セマンティクス，`T(int)` 構築，`+` `-` `*` 単項 `-`，`==` `!=`
- R2: 除算 `/`
- R3: 順序比較 `<` `<=` `>` `>=`
- R4: ADL の `abs` / `sqrt`

これに加えて tlapack では以下の **R5・R6 が新たに増えます**．
R7 は新しい演算ではなく意味論上の注意です．

### R5: 機械定数 — `std::numeric_limits<T>` の特殊化 (新規・最重要)

LAPACK の反復 algorithm は「収束した」「scaling が必要」を機械定数との比較で
判定するため，以下が意味のある値を返すことを要求します:

| 関数/定数 | 意味 | tlapack での用途 (`tlamch<T>` 経由) |
|---|---|---|
| `epsilon()` | 1 と次の表現可能数の差 | 収束判定 (`tsterf` `tsteqr` `tbdsqr` `tlahqr`)，rank 判定 (`tgels`) |
| `min()` | 最小の正規化数 (safmin) | 除算 guard (`tgetrf` の sfmin，`tlarfg` `tlartg` `tgebal` `tlaln2` `tlascl`) |
| `max()` | 最大の有限数 | overflow guard (`tlapy2` `tladiv` `ttrevc`) |
| `radix` / `digits` | 基数と仮数桁数 | `tlamch('B')` ほか派生定数 |

`tlamch<T>(cmach)` を `std::numeric_limits<T>` から構成して提供します
(`'E'` = epsilon()/2 (丸め単位)，`'S'` = min()，`'O'` = max() など，
dlamch と同じ文字規約)．派生定数 (smlnum = safmin/eps，bignum = 1/smlnum，
rtmin = sqrt(safmin) など) はすべて `tlamch<T>` の値から T の四則と sqrt で
計算するので，追加の要件は生じません．

- `double` / `float` / `long double`: 標準で充足
- `kv::dd`: kv が特殊化済み (epsilon = 2^-105，min = DBL_MIN，max ≈ DBL_MAX)
- `kv::mpfr<N>`: kv が特殊化済み (epsilon = 2^(1-N)，min/max は mpfr の emin/emax 由来)
- それ以外の型: 利用者が `std::numeric_limits<T>` を特殊化すれば使用可能

注意 (kv::dd): dd の数値モデルでは min() = DBL_MIN 付近で下位語が underflow
するため「safmin から max まで一様に相対精度 eps」という LAPACK の前提は
厳密には成り立ちません．極端な scale の入力 (|x| < 2^-916 程度) では精度が
落ちる可能性がありますが，通常 scale の入力では問題になりません．

### R6: NaN 判定 — `x == x` が NaN で false を返すこと (新規)

`tpotrf` (正定値性破綻の検出)，`tlange` 系の max norm (NaN 伝播)，
`tlassq` `tlapy2` `tlascl` が NaN 判定を必要とします．kv の型には `isnan` が
ないため，tlapack は既定で

```cpp
template <typename T> bool tisnan(const T& x) { return !(x == x); }
```

を使います (double は `std::isnan` に委譲)．したがって T への要件は
「**NaN を保持しうる型では，NaN との `==` 比較が false を返すこと**」です．
kv::dd (成分 double の比較) と kv::mpfr (mpfr_cmp は NaN で不等) はこれを
満たします．NaN が存在しない型ではこの要件は自動的に満たされます (常に false
で問題ない)．より適切な判定を持つ型は ADL の `isnan` で override できます．

### R7: 符号転送 sign(a, b) — 新しい演算は不要 (意味論上の注意のみ)

Fortran の `SIGN(A, B)` (rdlapack では `std::signbit` で実装) は

```cpp
template <typename T> T tsign(const T& a, const T& b) {
	using std::abs;
	return (b < T(0)) ? -abs(a) : abs(a);
}
```

と R1/R3/R4 のみで構成します．`std::signbit` と異なり **b = -0.0 を正と扱う**
点だけが rdlapack と相違しますが，sign は `tlarfg` (Householder の beta の
符号選択) 等で使われるだけなので，-0.0 で逆符号になっても結果は同値な分解に
なるだけで正しさには影響しません．

### 要件に「ならない」もの (実装側で回避する)

rdlapack が double 専用で使っている以下は，T への要件にしません:

1. **`std::ldexp` (2^k 定数)**: `rdlassq` の Blue's algorithm 定数
   (2^-511 など) と `rdlartg` の safmn2 = 2^-485 は double 専用 tuning．
   tlapack では
   - `tlassq` は Blue's algorithm をやめ，`tnrm2` と同じ古典的な
     scale/ssq 更新 (abs / 比較 / 除算のみ) で実装する
   - `tlartg` の scaling 定数は `sqrt(tlamch('S'))` 等の traits 由来の値を
     T の演算で構成する
2. **`std::pow`**: `rdbdsqr` の tolmul = eps^(-1/8) のみ．
   `1 / sqrt(sqrt(sqrt(eps)))` (sqrt 3 回 + 除算) で代替する．
3. **`std::sort` の比較**: `tlasrt` は `std::sort` + T の `<` (R3) で済む．
4. **`std::fma`**: tblas と同様使わない (`a * b + c` に展開)．

### 区間型は対象外

LAPACK は pivot 選択・収束判定・shift 選択など，ほぼ全 routine が順序比較で
分岐するため，`kv::interval` 系を tlapack に入れるのは数学的に不適切です
(分岐が不定になり，包含も保証されない)．tlapack の対象は点の型
(`double` `float` `long double` `kv::dd` `kv::mpfr<N>` など) に限定します．
区間演算による精度保証は，tlapack の近似解を入力とする検証法 (vcp 本体の
定式化) 側で行ってください．

具体的には `kv::interval<double>` は **R5 を満たしません**: kv は
`std::numeric_limits<kv::interval<T>>` を特殊化していないため primary template
に fallback し，`epsilon()` 等が [0,0] を返します．その結果 `tlamch<itv>('S')`
の内部の `1/max()` が区間の 0 除算となり，`tlamch` を使う全 routine (tgesv を
含む) は実行時に `std::domain_error` を投げます (compile は通るが使えない)．
仮に numeric_limits を自前で特殊化しても，kv の区間比較は「確実に小さい」
(`x.sup < y.inf`) の意味論なので重なった区間では分岐が機能せず，反復法の
収束判定は満たされず，deflation (要素への exact 0 代入) で包含も壊れるため，
意味のある結果にはなりません．

## routine 群と必要要件の対応 (増分のみ)

| routine 群 | R5 (tlamch) | R6 (NaN) | 備考 |
|---|---|---|---|
| `tlaset` `tlacpy` `tlaswp` `tlasrt` `ilatlr` `ilatlc` | | | R1/R3 のみ |
| `tgetrs` `tpotrs` `ttrtrs` `ttrtri` `tlauum` `tgetri` `tsytrs` | | | R1–R4 のみ (tblas 相当) |
| `tgetf2/tgetrf` 系，`tgbtrf` 系 | ✓ (sfmin) | | pivot は R3/R4 |
| `tpotf2/tpotrf` 系，`tpbtrf` 系 | | ✓ | 対角 <= 0 か NaN で info |
| `tlarfg` `tlartg` と QR/LQ/QL・`tgels` | ✓ | | scaling guard |
| norm 系 (`tlange` `tlansy` ...) `tlassq` `tlapy2/3` | ✓ | ✓ | max norm の NaN 伝播 |
| `tsyev` / `tsygv` 系 (`tsterf` `tsteqr` ...) | ✓ | | 収束判定 |
| `tgeev` 系 (`tgebal` `tlahqr` `tlaln2` `ttrevc` ...) | ✓ | | 収束判定 + overflow guard |
| `tgesvd` 系 (`tbdsqr` ...) | ✓ | | 収束判定 (pow は sqrt で代替) |
| `tsytf2/tsytrf` 系 | | | Bunch-Kaufman の alpha = (1+sqrt(17))/8 は T の演算で構成 |

## 実装方針 (要件から導かれるもの)

1. **`tlamch<T>` を单一の extension point にする**: 機械定数の参照は全て
   `tlamch<T>` 経由とし，`std::numeric_limits<T>` に直接触る箇所を
   `tlamch` の実装 1 箇所に集約する．
2. **派生定数は遅延計算しない**: 各 routine の冒頭で `const T eps = ...` と
   して T の演算で構成する (rdlapack と同じ形)．
3. **BLAS は全て tblas**: rdblas 呼び出し (`rdgemm(..., rm)`) を tblas
   (`vcp::tgemm(...)`) に置き換える．`idamax` → `itamax` など命名も対応．
4. **block size は rdlapack の `ilaenv` 既定値をそのまま使う** (GETRF/POTRF/
   SYTRF 64，QR/SYTRD/GEBRD/GEHRD 系 32)．多倍長型では演算 cost が支配的で
   block 化の意義は cache でなく BLAS3 化にあるため，同じ値で問題ない．
5. **file 構成は rdlapack を踏襲**: `tlapack.hpp` (umbrella) +
   `tlapack_common.hpp` / `_aux` / `_lu` / `_chol` / `_band` / `_qr` /
   `_sy` / `_eig` / `_geev` / `_svd`．
6. **rdlapack の reference との相違点を引き継ぐ**: `tbdsqr` は dqds でなく
   QR 反復，`tgesvd` は基本経路のみ，など (docs/rdlapack.md 参照)．

## まとめ: tblas からの要件の増分

| | 内容 | 充足状況 |
|---|---|---|
| **R5 (新規)** | `std::numeric_limits<T>` の `epsilon()` `min()` `max()` `radix` `digits` | double 系は標準，kv::dd / kv::mpfr は kv が特殊化済み |
| **R6 (新規)** | NaN との `==` が false (既定の `tisnan` = `!(x == x)` が機能すること)．ADL `isnan` で override 可 | double 系・kv::dd・kv::mpfr とも充足 |
| **R7 (注意のみ)** | sign(a,b) は比較 + abs で構成 (b = -0.0 の扱いだけ rdlapack と相違，結果の正しさに影響なし) | 追加要件なし |
| ldexp / pow / fma | 実装側で回避するため **T への要件にしない** | — |

つまり，**kv::dd と kv::mpfr<N> は追加の作業なしで tlapack の対象型になれます**．
新しい型を持ち込む場合のみ `std::numeric_limits` の特殊化が追加で必要です．

## 実装の詳細 (rdlapack からの移植)

vlapack/rdlapack (約 9,600 行) を以下の規則で機械変換し，特殊箇所を手で書き直した:

- rounding_mode 引数と RoundingGuard を全て除去，BLAS 呼び出しは tblas に置換
- `double` → `T`，FP literal は `T(int)` と四則から構成 (`0.95` → `T(95)/T(100)` 等)
- `DBL_EPSILON*0.5 / DBL_EPSILON / DBL_MIN / DBL_MAX` →
  `tlamch<T>('E' / 'P' / 'S' / 'O')`
- `std::fabs` / `std::sqrt` / `std::isnan` →
  `tlapack_detail::tabs / tsqrt / tisnan` (ADL 解決)
- LAPACK で行列名に使われる変数 `T` (larft/larfb の三角因子，trevc の Schur 形)
  は template parameter と衝突するため `Tf` に改名
- 作業領域の `std::vector<T>` は全て `T(0)` で fill
  (kv::mpfr の default constructor は NaN のため)
- **tlassq**: Blue's algorithm (2^k 定数が型依存) をやめ，古典的な scale/ssq
  更新で再実装 (呼び出し規約と scale^2*sumsq の不変量は同じ)
- **tlanv2**: scaling 閾値 safmn2 = 2^-485 (ldexp) を値ベースの
  `sqrt(tlamch('S')/eps)` に変更 (極端な scale の入力でのみ丸めが 1 回増える)
- **tbdsqr**: `pow(eps, -1/8)` を `1/sqrt(sqrt(sqrt(eps)))` に変更
- rdlapack が reference と異なる点 (tbdsqr は dqds でなく QR 反復，tgesvd は
  基本経路のみ等) はそのまま引き継ぐ

## テスト

2 本のテストで検証済み:

1. **T=double を rdlapack (rounding_mode=0) と比較** (86 checks PASS)．
   gesv/getri/posv/potri/sysv/gbsv/pbsv/trtrs/trtri/geqrf/orgqr/gels/syev/sygv/geev/gesvd
   について解・固有値・特異値の一致 (tlapack は fma を使わないため bit 一致ではなく
   許容誤差比較)，LU の pivot 一致，residual・直交性・再構成誤差を確認．
2. **kv::dd / kv::mpfr<106>** (42 checks PASS)．同じ問題を double と T で解き，
   解・固有値 (sort 済)・特異値・pivot 列が double の結果と一致すること
   (型を変えても結果が変わらない)，かつ residual / V^T V - I / A - U S VT が
   T の精度水準 (< 1e-24) まで落ちること (実際に高精度で解けている) を確認．
3. **T=double 特殊化の smoke test**:
   `sandbox/tests/test_tblas_tlapack_double_specialization.cpp` で
   `vcp/tblas/tblas_double.hpp` / `vcp/tlapack/tlapack_double.hpp` の include，
   `taxpy<double>`，`tgesv<double>`，および `USE_VCP_LAPACK` なしの場合の
   不足補助 routine (`tlae2<double>` → `dlae2_`) の実体化を確認．
   さらに `-DUSE_VCP_BLAS -DUSE_VCP_LAPACK` 付きでは，不足補助 routine が
   template<double> 実装へフォールバックして compile できることを確認．

---

# 関数リファレンス


`vcp/tlapack/tlapack.hpp` が提供する全 routine の仕様．設計方針と要素型 T の要件
(R1–R6) は このドキュメントの「設計方針」節 を参照．アルゴリズムの数学的内容は
reference LAPACK 3.12.1 の同名 routine (`t` → `d`) に準ずる．

## 共通仕様

- 全関数は `namespace vcp` に属する `template <typename T>` の関数 template
  (`tlamch<T>` 以外は引数から T が推論される)．`using namespace vcp;` を
  使えば `vcp::` を省略できる．
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
- **丸めモードは扱わない**．丸めモード指定が必要なら `vcp/vlapack/rdlapack` を使う．
- option 文字は大文字小文字どちらも可．
- scalar 引数は `const T&`，出力 scalar は `T&` 渡し．
- `namespace tlapack_detail` 内の関数 (`tabs` `tsqrt` `tisnan` `f_sign` `pow2i`
  `ilaenv` など) は内部実装用であり，この文書では公開 API のみ記載する．
  (内部 detail 名前空間はグローバルスコープに置かれ，`namespace vcp` には属さない．)

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
T vcp::tlamch(const char cmach);   // 例: vcp::tlamch<kv::dd>('E')
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
`ipiv` は `tgetrf` が返した 0-based pivot を渡す．

### tgetri — LU 分解からの逆行列
```cpp
template <typename T>
int tgetri(const int n, T* A, const int lda, const int* ipiv);
```
tgetrf の出力 A, ipiv から inv(A) を計算し A を上書きする．INFO > 0: 特異．
`ipiv` は `tgetrf` が返した 0-based pivot を渡す．

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
`tgbtf2` / `tgbtrf` / `tgbtrs` の `ipiv` は 0-based．

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
