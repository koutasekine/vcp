# vblas

`vblas/` は実験用のディレクトリであり、安定版の VCP Library API では
ありません。

現在のコードは、丸め方向を命令ごとに明示した低レベル BLAS 風カーネルを
試すためのものです。外部 BLAS ライブラリが、精度保証付き数値計算に必要な形で
ハードウェア丸めモードの変更を反映しない場合に、代替 backend として使えるかを
検討する意図で作成されています。

現時点での範囲は限定的です。

- `udmatmul.hpp` は `_mm512_fmadd_round_pd` などの AVX-512 intrinsic を使います。
- 上向き丸めと下向き丸めの double 行列積を、別々の出力行列へ同時に計算します。
- 対象は点 double 行列であり、一般の区間行列 backend ではありません。
- `vcp::matrix`、`vcp::pdblas`、`vcp::pidblas` から自動的に使われるものでは
  ありません。
- AVX-512 に対応していない CPU やコンパイラでは利用できません。

このディレクトリのファイルは、実験、ベンチマーク、将来の backend 設計のための
ものです。通常の VCP Library 利用では、`vcp/` と `docs/` 以下の文書を
参照してください。`vblas/` を、現時点でサポート済みの BLAS 代替実装として
扱わないでください。

`main.cpp` は `udmatmul.hpp` の小さな実験用 driver です。丸め方向指定カーネルの
結果と、ハードウェア丸めモードを変更した状態で `vcp::pdblas` により計算した
行列積を比較します。

`vblas/` のコードを将来 VCP の公開 API に昇格する場合は、少なくとも
compile guard、CPU feature check、丸め方向付き計算の正しさテスト、
policy としてのインターフェース設計が必要です。
