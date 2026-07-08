# bfem rational 厳密値と interval 包含テストの使い方

この文書は、`vcp/bfem` の「rational で作った厳密値」と「`kv::interval<double>` による包含」を再確認するための既存テストの使い方をまとめたものです。

ここで扱うテストは、新しい数値解法を解くためのテストではなく、同じ問題を

- `bfem_test::ratx` / `vcp::bfem::detail::rational` で厳密に組み立てる経路
- `kv::interval<double>` で区間として組み立てる経路

の2本で計算し、区間結果が rational 厳密値を含むかを確認する検証です。

## 対象ファイル

| 分野 | テストファイル | 既存バイナリ | 主に確認する範囲 |
|---|---|---|---|
| 2D RT E2E | `sandbox/tests/bfem_rt_e2e_smoke.cpp` | `sandbox/bin/bfem_rt_e2e_smoke_release` | 2D RT の大域行列、大域ベクトル、KKT 厳密解、全領域スカラー |
| 3D RT E2E | `sandbox/tests/bfem_d3c_e2e_smoke.cpp` | `sandbox/bin/bfem_d3c_e2e_smoke_release` | 3D RT の大域行列、大域ベクトル、KKT 厳密解、全領域スカラー |
| Scott-Vogelius | `sandbox/tests/bfem_sv_interval_smoke.cpp` | `sandbox/bin/bfem_sv_interval_smoke_release` | SV の大域行列、ベクトル、div 場、全領域スカラー |
| C1 要素 | `sandbox/tests/bfem_c1_e2e_tests.cpp` | `sandbox/bin/bfem_c1_e2e_tests_release` | C1 の行列、ベクトル、場、点評価、次数上げ、全領域スカラー |

補助ヘッダ:

- `sandbox/tests/bfem_ratx.hpp`: テスト用の rational スカラー `ratx`。
- `sandbox/tests/bfem_sparse_stub_policy.hpp`: `kv::interval<double>` でも疎行列の組み立てだけを確認できる storage-only policy。
- `sandbox/tests/bfem_test_framework.hpp`: `BFEM_CHECK` と結果表示。

## 何を「厳密値」と呼んでいるか

このテスト群での rational 厳密値は、参照単体上だけに限定されません。各テストは、rational 座標・rational 係数から物理要素上の局所量を作り、それを大域自由度へ組み立てた値まで比較しています。

| スコープ | このテスト群での扱い |
|---|---|
| 参照単体上 | Bernstein 係数、質量・微分・RT/C1 テーブルなどの内部材料として使われる |
| 物理要素上 | 要素幾何を含む局所行列・局所ベクトル・局所スカラーの材料として使われる |
| 領域全体 | 大域行列、大域ベクトル、KKT 解、div 場、誤差ノルムなどを rational と interval で比較する |

つまり、下記の smoke / E2E テストで確認している主対象は「要素内の式だけ」ではなく、代表メッシュ上で組み立てた領域全体の量です。ただし、すべての次数・すべてのメッシュ・すべての API を網羅した数学的証明ではありません。

## 実行方法

推奨は Makefile の既存ターゲットを使う方法です。プロジェクトルートから実行します。

```bash
make -C sandbox bfem_rt_e2e_release
make -C sandbox bfem_d3c_release
make -C sandbox bfem_sv_release
make -C sandbox bfem_c1_release
```

上のターゲットは各分野の関連テストもまとめてビルドします。ビルド後、包含確認に直接関係する release バイナリだけを実行するには次を使います。

```bash
./sandbox/bin/bfem_rt_e2e_smoke_release
./sandbox/bin/bfem_d3c_e2e_smoke_release
./sandbox/bin/bfem_sv_interval_smoke_release
./sandbox/bin/bfem_c1_e2e_tests_release
```

関連テスト全体を debug / release の両方で実行したい場合は、次の Makefile ターゲットを使います。

```bash
make -C sandbox run_bfem_rt_e2e
make -C sandbox run_bfem_d3c
make -C sandbox run_bfem_sv
make -C sandbox run_bfem_c1
```

## 個別コンパイル例

Makefile を使わず、対象テストだけを直接コンパイルする場合の例です。`MKLROOT` が設定されている Ubuntu / WSL 環境を想定しています。

```bash
mkdir -p sandbox/bin

g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 \
sandbox/tests/bfem_rt_e2e_smoke.cpp \
-L${MKLROOT}/lib/intel64 \
-Wl,--no-as-needed \
-lmkl_intel_lp64 \
-lmkl_intel_thread \
-lmkl_core \
-liomp5 \
-lpthread \
-lm \
-ldl \
-lmpfr \
-fopenmp \
-o sandbox/bin/bfem_rt_e2e_smoke_release
```

他の3本も、入力ファイル名と出力ファイル名を置き換えれば同じ形式でコンパイルできます。

```text
sandbox/tests/bfem_d3c_e2e_smoke.cpp      -> sandbox/bin/bfem_d3c_e2e_smoke_release
sandbox/tests/bfem_sv_interval_smoke.cpp  -> sandbox/bin/bfem_sv_interval_smoke_release
sandbox/tests/bfem_c1_e2e_tests.cpp       -> sandbox/bin/bfem_c1_e2e_tests_release
```

## PASS の見方

各テストは最後に `PASS` を含む集計行を出します。代表的な成功時の形は次の通りです。

```text
[bfem_rt_e2e_smoke] 288 checks, 288 passed, 0 failed => PASS
[bfem_d3c_e2e_smoke] 5 checks, 5 passed, 0 failed => PASS
[bfem_sv_interval_smoke] 8 checks, 8 passed, 0 failed => PASS
[bfem_c1_e2e_tests] 14 checks, 14 passed, 0 failed => PASS
```

追加で、区間幅を確認するための表示が出ます。これは幅を評価するための情報であり、判定自体は epsilon ではなく包含の成否です。

```text
(73) 3D term 1 ||grad u - p||^2  in [...]
(73) 3D term 2 ||f - Pi f||^2    in [...]
div_norm_sq interval = [...] width ...
C-T8 zero-case interval residual: [...]
```

## 各テストの確認内容

### `bfem_rt_e2e_smoke.cpp`

2次元の RT E2E 検証です。

- rational 側で `mesh<2, ratx>`、`fe_space<2, ratx>`、`rt_space<2, ratx>`、`broken_space<2, ratx>` を作る。
- RT 質量行列 `P`、div 行列 `N`、broken load `f_v` を rational で組み立てる。
- `rational_la::solve_exact` で KKT 系を厳密に解く。
- `flux_error_sq` と `projection_error_sq` の rational 厳密値を作る。
- interval 側で同じ量を `kv::interval<double>` と `sparse_stub` で組み立てる。
- `P`、`N`、`f_v`、`flux_error_sq`、`projection_error_sq` の各 interval が rational 厳密値を含むことを確認する。

### `bfem_d3c_e2e_smoke.cpp`

3次元の RT E2E 検証です。

- 3D の2四面体メッシュで rational 側と interval 側を構成する。
- RT 質量行列 `P`、div 行列 `N`、broken load `f_v` を比較する。
- rational 側で KKT 系を厳密に解き、その解係数を interval に一度だけ包含変換する。
- `flux_error_sq` と `projection_error_sq` の interval が rational 厳密値を含むことを確認する。

### `bfem_sv_interval_smoke.cpp`

Scott-Vogelius 関連の interval 包含検証です。

- rational 側で criss-cross メッシュと `vfe_space` / `broken_space` を作る。
- interval 側で同じメッシュを dyadic 座標として作る。
- `assemble_vector_stiffness`、`assemble_div_velocity`、`assemble_advection`、`assemble_advection_derivative` の全成分包含を確認する。
- `advection_vector`、`advection_scalar`、`div_norm_sq`、`div_field` の包含を確認する。

### `bfem_c1_e2e_tests.cpp`

2次元 C1 要素の interval E2E 検証です。

- `c1_space<2, ratx>` と `c1_space<2, kv::interval<double>>` を同じ dyadic メッシュ上で作る。
- `stiffness`、`laplacian_matrix`、`hessian_matrix`、`mixed_mass`、`laplacian_mixed` の全成分包含を確認する。
- `load`、`laplacian_load`、`inner`、`scalar_ff`、`laplacian_residual_sq` の包含を確認する。
- `laplacian_field` の係数包含を確認する。
- `eval`、`eval_grad`、`eval_hess` の点評価包含を確認する。
- `elevate` の係数包含を確認する。
- harmonic 構成のゼロケースで、区間 residual が rational 厳密ゼロを含むことを確認する。

## 判定が epsilon でないことの確認ポイント

テスト内の包含判定は、おおむね次の形です。

1. `kv::interval<double>` の下端・上端を double の2進小数として rational に戻す。
2. `lower_as_rational <= exact_rational <= upper_as_rational` を比較する。
3. 失敗した場合は `BFEM_CHECK` が failed を数える。

したがって、表示される区間幅は参考情報ですが、PASS/FAIL 判定は「近いかどうか」ではなく「厳密値を含むかどうか」です。

## 注意点

- `kv::interval<double>` で疎行列を扱うため、行列組み立て側ではテスト専用の `bfem_test::sparse_stub` を使います。これは疎ソルバではなく、COO 形式の storage-only policy です。
- KKT 系の解法そのものは interval ソルバで解いているのではなく、rational 側で `rational_la::solve_exact` により厳密解を作り、その係数を interval に包含変換して後段の量を確認しています。
- `ratx` はテスト用スカラーです。bfem の厳密有理数 `vcp::bfem::detail::rational` を包み、テストコードが必要とする最小限の変換を追加しています。
- smoke と名が付くテストは、代表ケースでの E2E 包含確認です。全次数・全メッシュ・全 API の網羅を意味しません。

