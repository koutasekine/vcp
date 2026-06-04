# Legacy Installer Notes

以前の VCP Library では、`installer/` に Ubuntu や CentOS 向けの補助
インストールスクリプトを置いていました。これらのスクリプトは、開発環境、
Boost、GMP、MPFR、Intel MKL/oneAPI などの導入と、`kv/`、`vcp/`、
`test_matrix/`、`test_PDE/` の展開をまとめて行うものでした。

しかし、Linux distribution、Intel oneAPI/MKL、パッケージリポジトリ、
GPG key の扱いは時間とともに変わるため、OS 別の実行可能 installer を
現行ツリーに残すと、現在の推奨手順として誤解される可能性があります。
そのため、OS や MKL に依存する installer script は廃止しました。

現在の基本方針は次の通りです。

- VCP Library 本体はヘッダを配置するだけで利用します。
- `vcp/` と `kv/` は同じ親ディレクトリに置く配置を推奨します。
- 開発環境、Boost、GMP、MPFR は distribution の package manager で
  導入します。
- BLAS/LAPACK、OpenBLAS、Intel MKL、OpenMP の設定は
  [blas_build.md](blas_build.md) を参照します。
- GitHub archive や tar の展開に慣れていない場合は、
  [../tools/download_vcp_kv.sh](../tools/download_vcp_kv.sh) を使えます。

旧 installer script で行っていた有用な処理のうち、`kv/` と `vcp/` を
同じ親ディレクトリに展開する部分だけを `tools/download_vcp_kv.sh` として
残しています。この補助ツールは `sudo`、`apt`、`dnf`、`.bashrc` の編集、
Intel MKL/oneAPI の導入を行いません。
kv ライブラリは [mskashi/kv](https://github.com/mskashi/kv) の GitHub archive から
取得します。

過去の OS 別 installer script が必要な場合は、Git の履歴から参照してください。
