# 外部ライブラリへの委譲(UMFPACK・ARPACK・CHOLMOD)

この文書では、`vcp::spumar` と `vcp::spcmodumar` が計算を外部ライブラリに
任せる仕組みを説明します。

| 内容 | 場所 |
| --- | --- |
| `spumar`(UMFPACK・ARPACK) | `vcp/spumar.hpp` |
| UMFPACK の呼び出し | `vcp/spumar_base/spumar_umfpack.hpp` |
| ARPACK の呼び出し | `vcp/spumar_base/spumar_arpack.hpp` |
| 格納形式の変換 | `vcp/spumar_base/spumar_convert.hpp` |
| `spcmodumar`(CHOLMOD) | `vcp/spcmodumar.hpp`、`vcp/spcmodumar_base/spcmodumar_cholmod.hpp` |

利用者が直接呼ぶことは想定していません。使い方とライセンスの注意は
[spmatrix.md](../spmatrix.md) の「`spumar` の詳細」「`spcmodumar` の詳細」
「外部ライブラリのライセンスに関する注意」を参照してください。

## 共通の性質

- **要素型**:`double` だけです。添字型は `std::int64_t` です(UMFPACK・CHOLMOD の
  64 ビットの添字の API、`umfpack_dl_*`・`cholmod_l_*` を使います)。ARPACK に
  渡す次元は `int` で、`int` に収まらない次元の固有値は委譲しません。
- **依存の切り離し**:外部ライブラリのヘッダを include するのは
  `spumar_base/spumar_umfpack.hpp`(`<suitesparse/umfpack.h>`)と
  `spcmodumar_base/spcmodumar_cholmod.hpp`(`<suitesparse/cholmod.h>`)だけ
  です。ARPACK は、`spumar_arpack.hpp` が Fortran の関数を宣言して使います。
  `spmats.hpp` などの本体のヘッダはこれらを include しないので、`spumar` /
  `spcmodumar` を使わないコードは外部ライブラリなしでビルドできます。
- **上書きする関数**:基底の `policy_X_impl` のうち、委譲するものだけを
  上書きします([sparse_layers.md](sparse_layers.md))。上書きしない計算
  (LDLᵀ 分解、慣性、LU 因子による解法、積など)は `spmats<double>` の実装を
  そのまま使います。
- **失敗の返し方**:外部ライブラリの失敗は、基底と同じく status で返します。
  外部ライブラリが返した値(UMFPACK の状態と Info、ARPACK の info・ierr・
  iparam など)は、`last_delegate_info` にそのまま記録します。委譲せずに
  基底の実装に戻したときは、その理由も記録します。
- **スレッド**:`last_delegate_info` は const の関数からも書き込むので、同じ
  オブジェクトを複数のスレッドから使うことはできません。
  `set_default_options()` はプロセス全体の既定値を変えます。

## `spumar` の線形方程式と LU 分解

### UMFPACK の包み方

`spumar_umfpack.hpp` は、UMFPACK の dl(64 ビットの添字・double)の API を包みます。
記号分解と数値分解の handle を持つ RAII のクラスで、1 回分解して何度も
解けるようにしています。`spumar_convert.hpp` は、`spmats` の CSC を
そのまま UMFPACK に渡します(どちらも 0 始まりの CSC で、行の添字は昇順・
重複なしなので、変換は CSC のコピーだけです)。

### 線形方程式

`policy_lss_with_info_impl` と `policy_lss_multi_with_info_impl` は、
`linear_solve_options` の `method` によらず、UMFPACK で分解して解きます
(`method` の値は結果にそのまま記録します)。右辺が複数のときは、1 回の分解を
共有します。

- `converged` は、分解と解法の status で決めます。
- 残差は、元の行列で計算し直して結果に入れます(参考値)。
- `muar_opt.umf_irstep` を 0 以上にすると、UMFPACK の反復改良の回数
  (`Control[UMFPACK_IRSTEP]`)を指定できます。`muar_opt.umf_pivot_tol` で
  ピボットの閾値を指定できます。負の値は UMFPACK の既定値です。

### LU 分解(因子を返す)

`policy_lu_with_info_impl` は、UMFPACK の分解から `umfpack_dl_get_numeric`
で因子を取り出し、P A Q = L U(p, q は new→old、L は対角の 1 も格納する
単位下三角)の形にして返します。

- 尺度は使いません(`Control[UMFPACK_SCALE] = UMFPACK_SCALE_NONE`)。
- `sparse_lu` のオプションの `method` と `ordering` は、UMFPACK に対応する
  ものがないので使いません。結果の `method_used`・`ordering_used` には要求した
  値が入り、使わなかったことを `last_delegate_info` に記録します。
- `equilibration = true` は、分解せずに `unsupported_options` を返します。
- UMFPACK が特異と判定したときは `invalid_factorization` を返します。
- UMFPACK は L を行ごとの形で返すので、CSC に転置します。厳密零は格納
  しません。

## `spumar` の固有値問題

`spumar_arpack.hpp` は、ARPACK-NG の逆通信の反復を 1 つの関数で回します。
作用素 OP と B の積は関数として渡します。

- 対称な問題は `dsaupd` / `dseupd`、非対称な問題は `dnaupd` / `dneupd` を
  使います。対称かどうかは `structure` で決めます(`auto_detect` のときは
  `is_symmetric()`、一般化では B も対称なとき)。
- 開始ベクトルは決まった式(resid_i = 1 + 0.01 · ((17 i) mod 101))で作るので、
  同じ入力には同じ結果を返します。
- 部分空間の次元 ncv は、既定で min(n, max(2k + 1, 20)) です。
  `muar_opt.ncv` が 0 でなければ min(n, max(ncv, k + 2)) です。
- `converged` は ARPACK の判定(info = 0 で、要求した個数が収束)です。
  固有対の残差は元の行列で計算し直して、診断として入れます。
- info = 1(反復の上限に達した)のときは、`converged = false` とし、収束した
  分だけを返します。
- 複素共役対は、基底と同じく `allow_complex_pairs = true` のときだけ返します
  (一般化問題では返しません)。

### 標準固有値問題の振り分け

| 条件 | ARPACK の使い方 |
| --- | --- |
| `use_shift = true` | mode 3。OP = (A − σI)⁻¹、which = LM。A − σI は UMFPACK で 1 回だけ分解する |
| `use_shift = false`、`target` が `smallest_magnitude` | mode 3、σ = 0 |
| `use_shift = false`、`target` が `target_real` | mode 3、σ = `shift` |
| `use_shift = false`、`target` が `largest_magnitude`・`largest_algebraic`・`smallest_algebraic` | mode 1(OP = A)。which は LM・LA・SA(非対称では LM・LR・SR) |
| `target_magnitude`、`dense_fallback_explicit`、k > n − 2、k = 0、n = 0 | `spmats<double>` の実装に戻す |

### 一般化固有値問題の振り分け

| 条件 | ARPACK の使い方 |
| --- | --- |
| `use_shift = true` | mode 3。OP = (A − σB)⁻¹ B。A − σB は UMFPACK で 1 回だけ分解する(特異なら `factorization_failed`) |
| `use_shift = false`、`target` が `smallest_magnitude` | mode 3、σ = 0 |
| `use_shift = false`、`target` が `largest_magnitude`・`largest_algebraic`・`smallest_algebraic` | mode 2。OP = B⁻¹ A。B を UMFPACK で 1 回だけ分解する |
| それ以外 | `spmats<double>` の実装に戻す |

### 前処理付きの版

前処理 M を渡す版では、shift-invert の内側の線形方程式を、UMFPACK ではなく
M を使った GMRES(相対許容誤差 tol/10)で解きます。内側の GMRES が収束
しなくても ARPACK の反復は止められないので、失敗した回数を数え、反復の後で
`converged = false`、`status = "inner_solve_failed"` として固有値を返しません。
UMFPACK には自動では切り替えません。

## `spcmodumar` の Cholesky 分解

`spcmodumar` は `spumar` を継承し、`policy_chol_with_info_impl` だけを
上書きして、疎 Cholesky 分解を CHOLMOD に任せます。`is_spd` も、この
Cholesky 分解を使います。

分解は次の決まった順に進みます。

1. **オプションの検査**:`pd_tol` が 0 でなければ `invalid_options`
   (CHOLMOD に対応する指定がないため)。`method` は `auto_select`(CHOLMOD の
   AUTO)と `simplicial_uplooking`(CHOLMOD の SIMPLICIAL)だけを受け付けます。
2. **入力の検査**:CSC の形と、非零の数が添字の型に収まるか。
3. **対称性の検査**:CHOLMOD に渡す前に、VCP 側で対称と言えるか確かめる
   (`check_symmetry`、`symmetry_tol`)。
4. **空の行・列の検査**:あれば `not_positive_definite`(`structural_empty_at`
   に記録)。
5. **並べ替えの準備**:

   | `ordering` | CHOLMOD |
   | --- | --- |
   | `natural` | NATURAL |
   | `amd`、`auto_select` | AMD |
   | `rcm` | VCP の RCM の並べ替えを与える(GIVEN) |
   | `nested_dissection` | NESDIS。CHOLMOD に Partition モジュールがなければ `invalid_options` |

6. **CHOLMOD での分解**:下三角だけを渡す。`method_used` には、実際に使った
   方法(CHOLMOD が supernodal 法を使えば `supernodal`)を記録する。
7. **status の対応**:成功なら `success`、正定値でなければ
   `not_positive_definite`(失敗した列を `failure_at` に記録)、それ以外は
   `internal_error`。`inconclusive_pivot_test` は返さない。
8. **結果の取り出し**:並べ替えを読み戻して全単射か確かめ、因子を下三角・
   行の昇順になっているか確かめて CSC に写す。厳密零は格納しない。

CHOLMOD の作業領域・因子・疎行列は RAII で管理し、どの経路でも解放します。
