# tlapack — テンプレート版 LAPACK (T 型要件と実装)

`tlapack` は `vlapack/rdlapack` (reference LAPACK 3.12.1 由来) と同じ routine
一覧を template 化した header-only library です．本文書は要素型 `T` に要求する
演算 (concept 相当) と実装・テストの要点を定めます．BLAS 部分は `tblas/` を
呼び出し，**丸めモードは一切意識しません** (rounding_mode 引数なし)．
**全関数の引数仕様は `tlapack/REFERENCE_tlapack.md` を参照．**

```cpp
#include "tlapack/tlapack.hpp"
#include <kv/dd.hpp>

// A*X = B を kv::dd で解く (column-major, ipiv は 0-based)
int info = tgesv(n, nrhs, A, lda, ipiv, B, ldb);   // T は引数から推論される
int info2 = tsyev('V', 'U', n, A, lda, w);          // 対称固有値 (kv::dd のまま高精度)
```

## file 構成

| file | 内容 |
|---|---|
| `tlapack.hpp` | umbrella header (これだけ include すればよい) |
| `tlapack_common.hpp` | 共通基盤 (`tlamch<T>` `tisnan` `f_sign` `pow2i` `ilaenv`) |
| `tlapack_aux.hpp` | 補助 routine (`tlassq` `tlartg` `tlarfg` `tlarf*` `tlange` 系など) |
| `tlapack_lu.hpp` `_chol` `_band` `_qr` `_sy` `_eig` `_geev` `_svd` | rdlapack と同じ分割 |

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
   (`tgemm(...)`) に置き換える．`idamax` → `itamax` など命名も対応．
4. **block size は rdlapack の `ilaenv` 既定値をそのまま使う** (GETRF/POTRF/
   SYTRF 64，QR/SYTRD/GEBRD/GEHRD 系 32)．多倍長型では演算 cost が支配的で
   block 化の意義は cache でなく BLAS3 化にあるため，同じ値で問題ない．
5. **file 構成は rdlapack を踏襲**: `tlapack.hpp` (umbrella) +
   `tlapack_common.hpp` / `_aux` / `_lu` / `_chol` / `_band` / `_qr` /
   `_sy` / `_eig` / `_geev` / `_svd`．
6. **rdlapack の reference との相違点を引き継ぐ**: `tbdsqr` は dqds でなく
   QR 反復，`tgesvd` は基本経路のみ，など (vlapack/README_rdlapack.md 参照)．

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

vlapack/rdlapack (約 9,600 行) を以下の規則で機械変換し
(`sandbox/tmp/convert_tlapack.py`)，特殊箇所を手で書き直した:

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

`make -C sandbox run_tlapack` で 2 本のテストを build + 実行する:

1. `sandbox/tests/test_tlapack.cpp` — **T=double を rdlapack (rounding_mode=0)
   と比較** (86 checks PASS)．gesv/getri/posv/potri/sysv/gbsv/pbsv/trtrs/trtri/
   geqrf/orgqr/gels/syev/sygv/geev/gesvd について解・固有値・特異値の一致
   (tlapack は fma を使わないため bit 一致ではなく許容誤差比較)，LU の pivot 一致，
   residual・直交性・再構成誤差を確認．
2. `sandbox/tests/test_tlapack_types.cpp` — **kv::dd / kv::mpfr<106>**
   (42 checks PASS)．同じ問題を double と T で解き，解・固有値 (sort 済)・
   特異値・pivot 列が double の結果と一致すること (型を変えても結果が変わらない)，
   かつ residual / V^T V - I / A - U S VT が T の精度水準 (< 1e-24) まで
   落ちること (実際に高精度で解けている) を確認．
