# File I/O Guide

この文書では、VCP Library の行列ファイル入出力を説明します。

VCP には 2 種類の保存形式があります。

| API | 形式 | 主な用途 |
| --- | --- | --- |
| `save` / `load` | 従来形式 | 既存ファイルとの互換性を保つ |
| `save_portable` / `load_portable` | portable binary 形式 | 異なる OS、CPU、実行環境間でのデータ移行 |

どちらも text 形式ではありません。数値を 10 進文字列に変換しないため、
文字列化による丸め誤差は入りません。

## 必要なヘッダ

通常は `matrix.hpp` の後に `matrix_assist.hpp` を include します。

```cpp
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

`matrix_assist.hpp` から、従来形式の `vcp_fio.hpp` と portable 形式の
`vcp_portable_fio.hpp` が読み込まれます。

portable 形式だけを明示的に使う場合は、次の include でも使えます。

```cpp
#include <vcp/matrix.hpp>
#include <vcp/vcp_portable_fio.hpp>
```

kv の型を保存する場合は、先に必要な kv ヘッダと VCP の policy を
include してください。

```cpp
#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>

#include <vcp/imats.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
```

## 従来形式

従来形式は、これまでの VCP のファイルと互換性を持つ保存形式です。

```cpp
vcp::save(A, "A");
vcp::load(B, "A");
```

実際のファイル名は、行列の要素型によって変わります。

| 行列型 | ファイル名 |
| --- | --- |
| `vcp::matrix<int, P>` | `A.matrix_int` |
| `vcp::matrix<double, P>` | `A.matrix_d` |
| `vcp::matrix<kv::interval<double>, P>` | `A.matrix_id` |
| `vcp::matrix<kv::dd, P>` | `A.matrix_kvdd` |
| `vcp::matrix<kv::interval<kv::dd>, P>` | `A.matrix_ikvdd` |

従来形式は既存データを読むために残しています。新しく環境間で大きな
データを移動する場合は、portable 形式を推奨します。

## Portable 形式

portable 形式は、環境間の移行を想定した binary 形式です。

```cpp
vcp::save_portable(A, "A");
vcp::load_portable(B, "A");
```

実際のファイル名は `A.vcpmat` です。従来形式の `.matrix_*` ファイルとは
別ファイルになるため、既存ファイルは変更されません。

portable 形式は次の情報を持ちます。

```text
magic
version
type id
layout
rows
columns
payload byte size
payload
CRC64
```

データは列優先で保存します。`double` は IEEE754 binary64 の bit pattern を
little-endian で保存します。`kv::dd` は `a1` と `a2` をそれぞれ binary64
として保存します。

対応している型は次の 5 種類です。

| 行列型 | 保存内容 |
| --- | --- |
| `vcp::matrix<int, P>` | int32 little-endian |
| `vcp::matrix<double, P>` | binary64 |
| `vcp::matrix<kv::interval<double>, P>` | `lower`, `upper` の binary64 |
| `vcp::matrix<kv::dd, P>` | `a1`, `a2` の binary64 |
| `vcp::matrix<kv::interval<kv::dd>, P>` | `lower.a1`, `lower.a2`, `upper.a1`, `upper.a2` の binary64 |

`kv::mpfr<N>`、`kv::interval<kv::mpfr<N> >`、`kv::ddx` は portable 形式の
対象外です。

## 保存時の安全性

portable 形式では、保存時にまず一時ファイルへ書き込みます。

```text
A.vcpmat.tmp
```

書き込み、CRC64 の保存、`close()` まで成功した場合だけ、
`A.vcpmat` に置き換えます。保存途中で失敗した場合に、既存の `A.vcpmat` を
壊しにくくするためです。

読み込み時には次を確認します。

- magic と version が正しいか
- 読み込み先の型とファイル内の型が一致するか
- 行数、列数、payload byte size が妥当か
- CRC64 が一致するか
- CRC64 の後ろに余分なデータがないか

失敗した場合は `vcp::io_error` を投げます。

## Round-Trip Test

同じ環境で保存して読み戻す最小確認は次の形です。

```cpp
vcp::save_portable(A, "A");
vcp::load_portable(B, "A");

for (int i = 0; i < A.rowsize(); i++) {
    for (int j = 0; j < A.columnsize(); j++) {
        if (A(i, j) != B(i, j)) {
            std::cout << "file I/O mismatch" << std::endl;
        }
    }
}
```

`kv::dd` や `kv::interval<kv::dd>` で bit 単位の復元を確認したい場合は、
内部成分も確認できます。

```cpp
if (A(i, j).a1 != B(i, j).a1 || A(i, j).a2 != B(i, j).a2) {
    std::cout << "file I/O mismatch" << std::endl;
}
```

```cpp
if (A(i, j).lower().a1 != B(i, j).lower().a1 ||
    A(i, j).lower().a2 != B(i, j).lower().a2 ||
    A(i, j).upper().a1 != B(i, j).upper().a1 ||
    A(i, j).upper().a2 != B(i, j).upper().a2) {
    std::cout << "file I/O mismatch" << std::endl;
}
```

異なる環境へ巨大ファイルを移動する場合は、ファイル転送後にファイルサイズや
checksum も確認すると原因切り分けがしやすくなります。

```bash
ls -lh A.vcpmat
sha256sum A.vcpmat
```

portable 形式の CRC64 は読み込み時の破損検出に使われます。`sha256sum` は、
転送前後のファイルが同一かを外側から確認するために使います。

## 使い分け

既存の `.matrix_*` ファイルを読む場合は `load` を使います。

```cpp
vcp::load(A, "A");
```

新しく保存し、別環境へ移動する可能性がある場合は `save_portable` を
推奨します。

```cpp
vcp::save_portable(A, "A");
```

従来形式と portable 形式はファイル名が違うため、併用できます。

```cpp
vcp::save(A, "A");           // A.matrix_d など
vcp::save_portable(A, "A");  // A.vcpmat
```
