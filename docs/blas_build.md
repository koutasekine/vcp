# BLAS / LAPACK Build Notes

この文書では、VCP Library を利用する際のコンパイラ、BLAS/LAPACK、
OpenMP、Apple Silicon Mac での OpenBLAS 設定例をまとめます。

## 対応コンパイラと検証環境

VCP Library は C++11 以降のコンパイラで利用できます。現在、C++11 から
C++20 までのコンパイルで動作を確認しています。gcc と clang の両方で
検証しています。

C++23 についても Ubuntu clang version 18.1.3 と gcc version 13.3.0 で
コンパイルと実行を確認しています。ただし、これらのコンパイラでは C++23 対応に
実験的な側面が残るため、ここでは主な確認済み範囲を C++11 から C++20 として
記載しています。

検証している環境には、Linux、macOS、Apple Silicon Mac、AMD CPU 環境が
含まれます。BLAS/LAPACK を使う policy では、利用環境に応じて OpenBLAS、
Intel MKL、Apple Accelerate BLAS/LAPACK などを選択できます。ただし、
精度保証付き数値計算で BLAS の丸めモード変更を使う場合は、後述の
丸めモード確認を必ず実行してください。

VCP Library 本体は Intel 専用ではありません。Intel MKL は選択肢の一つです。
現在の `tools/download_vcp_kv.sh` は、Intel MKL/oneAPI の導入や shell 設定を
行いません。

開発環境、Boost、GMP、MPFR などの基本パッケージの導入例は
[install.md](install.md) を参照してください。

## 基本のコンパイル例

以下の例では、最小対応規格として `-std=c++11` を指定しています。
検証済みの環境では、`-std=c++14`、`-std=c++17`、`-std=c++20` も利用できます。

`vcp/` ディレクトリを含む VCP Library の展開先を include path に追加します。
`vcp/` そのものを `-I` に指定しないでください。配置と include path の考え方は
[install.md](install.md) を参照してください。

```bash
g++ -std=c++11 -I/path/to/libs example.cpp
```

区間演算を使う場合は、別途取得した kv ライブラリの include path も
必要です。推奨配置は `vcp/` と `kv/` を同じ親ディレクトリに置く形で、
この場合は include path を 1 つ指定するだけで済みます。

```bash
g++ -std=c++11 -I/path/to/libs example.cpp
```

BLAS/LAPACK を使う場合:

```bash
g++ -std=c++11 -I/path/to/libs example.cpp -llapack -lblas
```

OpenMP、MPFR、GMP も使う場合:

```bash
g++ -std=c++11 -I/path/to/libs -O3 -DNDEBUG -DKV_FASTROUND example.cpp -llapack -lblas -lmpfr -lgmp -fopenmp
```

## Intel MKL を使う場合

Intel MKL を利用する場合、MKL を含む Intel oneAPI のインストール方法は
Intel 公式の
[Intel oneAPI Toolkits Installation Guide for Linux](https://www.intel.com/content/www/us/en/docs/oneapi-toolkit/installation-guide-linux/latest/overview.html)
を参照してください。VCP Library のドキュメントでは、インストール後の
リンク方法と丸めモード確認を中心に説明します。
APT を使う場合の具体的な導入手順は Intel 公式の
[Install with APT](https://www.intel.com/content/www/us/en/docs/oneapi-toolkit/installation-guide-linux/latest/install-oneapi-toolkit-with-apt.html)
を参照してください。

Intel oneAPI を package manager でインストールしただけでは、現在開いている
shell に `MKLROOT`、include path、library path などが設定されていない場合が
あります。コンパイル前に、oneAPI の環境設定スクリプトを source してください。

Component Directory Layout の典型例:

```bash
source /opt/intel/oneapi/setvars.sh
```

同じ意味で、次のようにも書けます。

```bash
. /opt/intel/oneapi/setvars.sh
```

Unified Directory Layout の場合は、toolkit version ごとの
`oneapi-vars.sh` を使います。

```bash
source /opt/intel/oneapi/<toolkit-version>/oneapi-vars.sh
```

この設定は、基本的には現在の terminal session にだけ有効です。新しい端末を
開いた場合は再度 source してください。毎回自動で設定したい場合は、環境に
応じて `~/.bashrc`、`~/.profile`、modulefiles などを使います。

設定できているかは、例えば次で確認できます。

```bash
echo ${MKLROOT}
ls ${MKLROOT}/include
```

`MKLROOT` が空の場合は、まだ現在の shell で MKL の環境設定が有効に
なっていません。

Intel の APT インストール手順は、主に package の導入手順を説明しています。
開発用 shell の環境初期化については、Intel の
[Configure Your CPU or GPU System](https://www.intel.com/content/www/us/en/docs/oneapi-base-toolkit/get-started-guide-linux/2025-2/before-you-begin.html)
や
[Use the setvars and oneapi-vars Scripts with Linux](https://www.intel.com/content/www/us/en/docs/oneapi/programming-guide/current/use-the-setvars-and-oneapi-vars-scripts-with-linux.html)
も参照してください。

MKL を精度保証付き計算で使う場合は、丸めモードの変更が有効に働くことが
重要です。環境によっては `-lblas` や `-llapack` を使うと、意図した MKL ではなく
別の BLAS/LAPACK に解決される場合があります。そのため、MKL を使う場合は
`-lmkl_intel_lp64`、`-lmkl_intel_thread`、`-lmkl_core` などを明示して
リンクする構成を推奨します。

推奨する g++ の例:

```bash
source /opt/intel/oneapi/setvars.sh
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 <filename.cpp> -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp
```

この例では `-lblas` と `-llapack` を指定していません。MKL を使う意図を
リンクオプション上でも明確にし、丸めモード変更を利用する計算で
想定外の backend が混ざることを避けるためです。

## インストール後の丸めモード確認

精度保証付き数値計算で BLAS を使う場合、BLAS の演算が丸めモード変更を
反映することが重要です。インストール後は、次の確認プログラムを実行することを
推奨します。

| ファイル | 確認内容 |
| --- | --- |
| `test_matrix/Check_pdblas_rounding.cpp` | `vcp::pdblas` で BLAS の丸めモード変更が効くか |
| `test_matrix/Check_pidblas_rounding.cpp` | `vcp::pidblas` で利用する BLAS 演算の丸めモード変更が効くか |

OpenBLAS や Apple LAPACK を使う例:

```bash
cd test_matrix
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 Check_pdblas_rounding.cpp -llapack -lopenblas -lmpfr -lgmp -fopenmp
./a.out

g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 Check_pidblas_rounding.cpp -llapack -lopenblas -lmpfr -lgmp -fopenmp
./a.out
```

MKL を使う場合は、上の MKL 用リンクオプションに
`Check_pdblas_rounding.cpp` または `Check_pidblas_rounding.cpp` を指定して
実行してください。`-lblas` や `-llapack` に置き換えず、MKL のライブラリを
明示的にリンクすることを推奨します。

Intel MKL を使う場合の確認例:

```bash
echo "Check for BLAS rounding mode changes: Please wait..."
source /opt/intel/oneapi/setvars.sh
cd test_matrix
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 Check_pdblas_rounding.cpp -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp && ./a.out
g++ -I.. -std=c++11 -DNDEBUG -DKV_FASTROUND -O3 -m64 Check_pidblas_rounding.cpp -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lmpfr -fopenmp && ./a.out
```

成功すると、丸めモード変更可能であることを示すメッセージが表示されます。
失敗した場合は、リンクされている BLAS/LAPACK の実体や、OpenBLAS の
`CONSISTENT_FPCSR` 設定、MKL のリンクオプションを確認してください。

## Linux x86_64 で OpenBLAS を使う場合

AMD や Intel の Linux 環境で OpenBLAS を利用し、BLAS の丸めモード変更を
反映させたい場合は、OpenBLAS を source build して
`CONSISTENT_FPCSR = 1` を有効にします。

Ubuntu 系環境での設定例:

```bash
sudo apt update
sudo apt upgrade -y
sudo apt install -y build-essential clang libboost-all-dev libgmp-dev libmpfrc++-dev git checkinstall

git clone https://github.com/xianyi/OpenBLAS.git
cd OpenBLAS/
vim Makefile.rule
make
sudo checkinstall --install=no
sudo dpkg -i openblas_xxxxxxxx_amd64.deb
```

`vim Makefile.rule` で編集する箇所は、次のコメントアウトを外す作業です。

```make
# If you need to synchronize FP CSR between threads (for x86/x86_64 and aarch64 only).
# CONSISTENT_FPCSR = 1
```

次のようにします。

```make
CONSISTENT_FPCSR = 1
```

インストール後、OpenBLAS の library path を通します。

```bash
sudo vim /etc/bash.bashrc
```

追記例:

```bash
export LIBRARY_PATH=/opt/OpenBLAS/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/OpenBLAS/lib:$LD_LIBRARY_PATH
```

dynamic linker にも OpenBLAS の場所を登録します。

```bash
sudo vim /etc/ld.so.conf.d/openblas.conf
```

内容:

```text
/opt/OpenBLAS/lib
```

設定後、必要に応じて `ldconfig` を更新し、OpenBLAS が見えているか確認します。

```bash
sudo ldconfig
ldconfig -p | grep libopenblas
env | grep PATH
```

最後に、`test_matrix/Check_pdblas_rounding.cpp` と
`test_matrix/Check_pidblas_rounding.cpp` を実行し、BLAS の丸めモード変更が
反映されることを確認してください。

## Apple Accelerate BLAS の注意

macOS では、Xcode Command Line Tools などの Developer Tools を導入すると、
Accelerate framework の BLAS/LAPACK を利用できる構成になります。
Apple 公式ドキュメントでも、Accelerate framework が BLAS/LAPACK を提供する
こと、および BLAS/LAPACK の threading model を指定できることは説明されています。

一方で、VCP Library の精度保証付き数値計算で重要になる
「BLAS 演算がハードウェア丸めモードの変更を反映すること」については、
Apple Accelerate BLAS で利用できる設定方法をまだ確認できていません。
少なくともこの文書の作成時点では、Accelerate BLAS に対して OpenBLAS の
`CONSISTENT_FPCSR` のような丸めモード同期・反映の設定方法は把握していません。

そのため、Apple Accelerate BLAS を使う場合は、精度保証付き数値計算に使う前に
必ず `test_matrix/Check_pdblas_rounding.cpp` と
`test_matrix/Check_pidblas_rounding.cpp` を実行してください。これらの確認に
失敗する場合、Accelerate BLAS を丸めモード変更に依存する精度保証付き計算の
BLAS backend として使う方法は、現時点では未確認です。

Apple Silicon Mac で BLAS の丸めモード変更が必要な場合は、次節のように
`CONSISTENT_FPCSR` を有効にした OpenBLAS を明示的に `-lopenblas` でリンクする
構成を推奨します。

参考:

- [Apple Accelerate](https://developer.apple.com/accelerate/)
- [Apple Accelerate BLAS](https://developer.apple.com/documentation/accelerate/blas-library)
- [Apple BLAS threading model](https://developer.apple.com/documentation/accelerate/blas_threading)

## Apple Silicon Mac での OpenBLAS

以下は、Apple Silicon Mac で丸めモード変更可能な OpenBLAS を使うための
設定例です。
`vcp::pidblas` など、丸めモード変更を利用する検証計算で OpenBLAS を使う場合に
参考になります。

設定例として記録した環境:

```text
MacBook Pro
チップ: Apple M3
OS: macOS Sonoma v14.1
Homebrew 4.2.21-84-g1c28136
OpenBLAS 0.3.27
```

同様の構成は M1、M3、M4 などの Apple Silicon Mac でも検証しています。

Homebrew の挙動を固定するため、`~/.zshrc` に次を追加します。

```zsh
export HOMEBREW_NO_AUTO_UPDATE=1
export HOMEBREW_NO_INSTALL_FROM_API=1
export EDITOR=vim
```

`HOMEBREW_NO_AUTO_UPDATE` は Homebrew の自動更新を止める設定です。
Homebrew を日常的に使っている環境では影響に注意してください。
`HOMEBREW_NO_INSTALL_FROM_API` は、formula の `.rb` ファイルを編集して
再インストールする際に必要です。`EDITOR` は `brew edit` で使うエディタです。

OpenBLAS をインストールし、formula を編集します。

```bash
brew install openblas
brew edit openblas
```

`brew edit` で開かれる formula は、典型的には次の場所にあります。

```text
/opt/homebrew/Library/Taps/homebrew/homebrew-core/Formula/o/openblas.rb
```

`def install` の冒頭付近で、次の 3 行を追加します。

```ruby
def install
  ENV.runtime_cpu_detection
  ENV.deparallelize

  ENV["CONSISTENT_FPCSR"] = "1"
  ENV["NO_LAPACK"] = "1"
  ENV["NO_LAPACKE"] = "1"
```

`CONSISTENT_FPCSR` は、丸めモード変更に関係する OpenBLAS 側のフラグです。
`NO_LAPACK` と `NO_LAPACKE` は、macOS に標準で存在する Apple LAPACK との
競合を避け、OpenBLAS は BLAS として使う意図の設定です。

さらに、`libblas` と `liblapack` の symlink が作られないように、次の行を
コメントアウトします。

```ruby
# lib.install_symlink shared_library("libopenblas") => shared_library("libblas")
# lib.install_symlink shared_library("libopenblas") => shared_library("liblapack")
```

これにより、OpenBLAS は `-lopenblas` で明示的に使い、`-lblas` や
`-llapack` では Apple BLAS/LAPACK を使う構成にしやすくなります。

編集後、OpenBLAS を source build で再インストールします。

```zsh
brew uninstall --ignore-dependencies openblas
brew reinstall --build-from-source openblas
```

OpenBLAS を Make するため、再インストールには時間がかかります。

最後に、`~/.zshrc` に OpenBLAS の include path と library path を追加します。

```zsh
export CPATH=$CPATH:/opt/homebrew/opt/openblas/include/
export LIBRARY_PATH=$LIBRARY_PATH:/opt/homebrew/opt/openblas/lib/
```

## Apple Silicon Mac でのコンパイル例

OpenBLAS、Apple LAPACK、MPFR、GMP、OpenMP を使う場合の例です。

Apple Clang:

```bash
clang++ -I.. -std=gnu++14 -O3 -Xpreprocessor -fopenmp Check_pdblas_rounding.cpp -llapack -lopenblas -lmpfr -lgmp -lomp
```

Homebrew GCC:

```bash
g++-14 -I.. -O3 -fopenmp Check_pdblas_rounding.cpp -llapack -lopenblas -lmpfr -lgmp
```

この構成では、`-lopenblas` で丸めモード変更可能な OpenBLAS を明示的に
使います。`-lopenblas` を `-lblas` に変えると Apple BLAS が使われ、
高速になる場合がありますが、丸めモード変更が効かないことがあります。

## OpenMP の注意

Apple Clang で OpenMP を使う場合は、Homebrew の `libomp` を使い、
コンパイル時に `-Xpreprocessor -fopenmp`、リンク時に `-lomp` を指定する
構成が典型的です。

Homebrew GCC では、通常の `-fopenmp` で OpenMP を有効にできます。

利用者側で外側のループを OpenMP 並列化する場合、VCP 内部の OpenMP を
止めるために次を定義できます。

```bash
-DVCP_NOMP
```

必要に応じて、より細かい単位で次も利用できます。

```bash
-DVCP_MATS_NOMP
-DVCP_FOURIER_NOMP
-DVCP_LEGENDRE_NOMP
```
