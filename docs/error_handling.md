# Error Handling

この文書では、VCP Library のエラーの知らせ方を説明します。VCP Library の
エラーは例外で知らせます。疎行列の分解などの一部の関数は、例外の代わりに
status を返す版も持ちます。

## 例外の階層

VCP 由来の例外は、`<vcp/error.hpp>` で定義された `vcp::error` から派生します。
`vcp::error` は `std::runtime_error` を継承しています。

```
std::runtime_error
└─ vcp::error
   ├─ vcp::invalid_argument
   │  └─ vcp::bfem::meshgen_error
   ├─ vcp::dimension_error
   ├─ vcp::index_error
   ├─ vcp::state_error
   ├─ vcp::domain_error
   │  └─ vcp::bfem::meshgen_convert_error
   ├─ vcp::numerical_error
   │  ├─ vcp::verification_error
   │  │  ├─ vcp::bfem::degenerate_element
   │  │  ├─ vcp::bfem::meshgen_degenerate
   │  │  ├─ vcp::bfem::sv_indeterminate_singularity
   │  │  ├─ vcp::bfem::c1_indeterminate_corner
   │  │  └─ vcp::btfem::degenerate_interval
   │  ├─ vcp::backend_error
   │  │  ├─ vcp::blas_error
   │  │  └─ vcp::lapack_error
   │  └─ vcp::bfem::meshgen_limit
   └─ vcp::io_error
      └─ vcp::bfem::meshgen_io_error
```

`vcp::bfem::` と `vcp::btfem::` のクラスは、bfem・btfem に固有の例外です。

| 例外 | 意味 |
| --- | --- |
| `vcp::invalid_argument` | 不正な引数 |
| `vcp::dimension_error` | 行列サイズの不一致 |
| `vcp::index_error` | 添字や部分行列範囲の誤り |
| `vcp::state_error` | 対応していない行列状態や policy 状態。反復法などが収束しなかった場合 |
| `vcp::domain_error` | 対称性など数学的仮定の破れ |
| `vcp::numerical_error` | 数値計算の失敗(分解の失敗など) |
| `vcp::verification_error` | 精度保証の検証条件の失敗 |
| `vcp::backend_error` | 外部ライブラリの routine の失敗。`backend()`、`routine()`、`code()` で詳細を取得できます |
| `vcp::blas_error` | BLAS routine の失敗 |
| `vcp::lapack_error` | LAPACK routine が非ゼロの `info` を返した。`info()` で値を取得できます |
| `vcp::io_error` | ファイル入出力やバイナリ形式のエラー |

`vcp::verification_error` は、精度保証の十分条件が成り立たなかったことを
意味します。対象についての主張(たとえば「解が存在しない」)ではありません。

## 例外の捕まえ方

```cpp
#include <vcp/error.hpp>

try {
    x = lss(A, b);
} catch (const vcp::verification_error& e) {
    std::cerr << "verification failed: " << e.what() << std::endl;
} catch (const vcp::lapack_error& e) {
    std::cerr << e.routine() << ": info = " << e.info() << std::endl;
} catch (const vcp::error& e) {
    std::cerr << "VCP error: " << e.what() << std::endl;
}
```

kv など外部ライブラリが投げる例外も含めてすべて捕まえる場合は、
`std::exception` で捕まえます。

```cpp
try {
    /* ... */
} catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
}
```

## 例外を投げる関数と status を返す関数

疎行列 `vcp::spmatrix` の多くの操作には、例外を投げる版と、名前が
`_with_info` で終わる版があります。

| 版 | 失敗したとき |
| --- | --- |
| 例外を投げる版(`solve`、`ldl`、`lu`、`chol`、`inertia`、`eigs` など) | 例外を投げます |
| `_with_info` 版(`solve_with_info`、`ldl_with_info` など) | 例外を投げず、結果の `status` または `converged` で知らせます |

例外を投げる版が投げる例外は次の通りです。

| 操作 | 例外 |
| --- | --- |
| `solve`、`eigs` が収束しない | `vcp::state_error` |
| `ldl`、`lu`、`chol`、`inertia`、`lu_solve` などで status が `success` 以外 | `vcp::numerical_error` |
| `lu_solve`、`lu_inverse_row` で、因子・置換・右辺のサイズ(または行番号)が合わない | `vcp::dimension_error` |

`_with_info` 版でも、行列が正方でない、右辺のサイズが合わないといった
誤った使い方には例外を投げます(`vcp::dimension_error` など)。

区間疎行列の policy `vcp::spimats<T>` は例外です。精度保証に失敗したときは、
`_with_info` 版でも `vcp::verification_error` を投げます。詳しくは
[spmatrix.md](spmatrix.md) の「`spimats` の詳細」を参照してください。

status の値と意味は、それぞれの関数の説明を参照してください。

## 破壊的な操作と例外

BLAS/LAPACK を使う一部の操作は、大きな行列で余分なメモリを使わないため
破壊的に実行されます。そのため、例外が投げられた時点で入力行列の内容が
backend routine によって上書き済みの場合があります。
