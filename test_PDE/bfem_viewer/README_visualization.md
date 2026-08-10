# bfem 近似解の可視化 (GRF-3)

`test_PDE/test_spnewton_*.cpp` が書き出す図示用データ (`.dat`) を ParaView で
表示するための道具です。C++ 側 (VCP 本体) には何も追加しません。

```
テスト実行          変換                      表示
*.cpp  ──▶  cells.dat / points.dat  ──▶  *.vtu  ──▶  ParaView
```

---

## 1. 必要なもの

| 用途 | 必要なもの | 確認方法 |
|---|---|---|
| 変換 (`.dat` → `.vtu`) | Python 3 と **numpy** | `python3 -c "import numpy; print(numpy.__version__)"` |
| 表示 (`.vtu` → 画面) | **ParaView(`pvpython` 込み)** | `pvpython --version` |

### numpy

多くの環境に既に入っています。無い場合:

| 環境 | コマンド |
|---|---|
| Debian / Ubuntu | `sudo apt install python3-numpy` |
| pip (どの OS でも) | `pip install numpy` |

### ParaView

**`pvpython` が必要です。** Debian / Ubuntu では `paraview` パッケージに
`pvpython` は含まれておらず、`python3-paraview` の側に入っています
(`python3-paraview` は `paraview` を依存に持つので、こちらを指定すれば GUI ごと
入ります)。

```sh
sudo apt update
sudo apt install python3-paraview
```

Ubuntu 24.04 では ParaView 5.11.2 が入ります。依存が 155 パッケージ、展開
約 500 MB と大きいので、それを避けたい場合は Kitware 公式のバイナリ
(`ParaView-*-MPI-Linux-x86_64.tar.gz`) を展開する方法もあります。こちらは
`pvpython` が同梱で、apt を汚しません。ただし版が新しいと `paraview.simple` の
API 名が変わっている可能性があり、動作確認済みなのは 5.11.2 です。

Windows / macOS でも同じスクリプトが動く作りですが、**確認していません**。
`pvpython` が見つからない場合は、変換だけ済ませて「`.vtu` は生成済み」と案内して
終了します。その `.vtu` を ParaView の GUI で開けば表示できます。

### 検証用 (通常は不要)

`sandbox/tests/grf3_vtu_check.py` を走らせる場合のみ `meshio` が要ります
(`pip install meshio`)。書き出した `.vtu` を独立した実装で読み戻して照合する
ために使うもので、道具そのものには不要です。

---

## 2. 置き場所

道具は 3 本です。**この 3 本が同じフォルダにあること**だけが条件で、フォルダの
場所は任意です。`test_PDE/` はリポジトリ内の正本の置き場所にすぎず、実行時には
関係しません。

```
<任意のフォルダ>/vcp_bfem_dat_to_vtu.py     cells.dat + points.dat -> .vtu
<任意のフォルダ>/vcp_bfem_view_vtu.py       .vtu -> 画面
<任意のフォルダ>/vcp_bfem_view_dat.py       上の 2 本をまとめて呼ぶ
```

`.dat` も同じフォルダに置くと、引数以外に何も打たずに済みます。

```
<任意のフォルダ>/vcp_bfem_dat_to_vtu.py
<任意のフォルダ>/vcp_bfem_view_vtu.py
<任意のフォルダ>/vcp_bfem_view_dat.py
<任意のフォルダ>/ns_3dsv_lsc_cells.dat
<任意のフォルダ>/ns_3dsv_lsc_velocity_points.dat
```


---

## 3. 使い方

`.dat` があるフォルダで実行します。**引数は `cells` が先、`points` が後**です。

```sh
python3 vcp_bfem_view_dat.py ns_3dsv_lsc_cells.dat ns_3dsv_lsc_velocity_points.dat
```

段階的に実行することもできます。

```sh
python3 vcp_bfem_dat_to_vtu.py ns_3dsv_lsc_cells.dat ns_3dsv_lsc_velocity_points.dat
pvpython  vcp_bfem_view_vtu.py  ns_3dsv_lsc_velocity.vtu
```

ウィンドウは閉じるまで残ります (`q` でも終了)。画像として残す場合:

```sh
python3 vcp_bfem_view_dat.py cells.dat points.dat --screenshot out.png
```

### テストごとの例

```sh
# 2 次元 Emden
python3 vcp_bfem_view_dat.py emden_2dfem_cells.dat emden_2dfem_points.dat

# 3 次元 Emden
python3 vcp_bfem_view_dat.py emden_3dfem_lsc_cells.dat emden_3dfem_lsc_points.dat

# Scott-Vogelius Navier-Stokes: 速度と圧力で 2 回。cells は共通
python3 vcp_bfem_view_dat.py ns_3dsv_lsc_cells.dat ns_3dsv_lsc_velocity_points.dat
python3 vcp_bfem_view_dat.py ns_3dsv_lsc_cells.dat ns_3dsv_lsc_pressure_points.dat
```

2 次元・3 次元・ベクトル場の区別は列数から自動判別されるので、コマンドの形は
どれも同じです。

### 主な引数

| 引数 | 意味 |
|---|---|
| `-n <名前>` | データ配列名と `.vtu` の名前 (既定は points のファイル名から) |
| `--glyph magnitude` | ベクトルの矢印の長さを大きさに比例させる (既定は一定長) |
| `--value lower` / `upper` / `mid` | 区間・有理数入力でどの値を表示するか (既定 `mid`) |
| `--merge always` / `never` | 節点併合の強制 (既定は自動判定) |
| `--screenshot <path>.png` | 画面表示の代わりに画像を保存 |

---

## 4. 既定の表示

| 対象 | 表示 |
|---|---|
| スカラー (Emden の解、SV の圧力) | 色付き曲面 |
| ベクトル (SV の速度) | **長さ一定の矢印**、色 = 大きさ、領域の稜線を重ねる |

矢印を長さ一定にしているのは、大きさに比例させると大半が見えなくなるためです
(SV 速度で実測: 領域幅の 2 % より短い矢印が 91 %)。大きさは色で読みます。

節点の併合は測って決めます。同じ座標の点で値が一致すれば連続な場と判断して
併合し、食い違えば不連続な場と判断して併合しません。SV の圧力は跳びが場の大きさの
77 % に達するため併合されません (跳びを平均で消さないため)。

---

## 5. うまくいかないとき

| 症状 | 原因と対処 |
|---|---|
| `paraview.simple is not importable` | `vcp_bfem_view_vtu.py` を `python3` で実行している。`pvpython` で実行するか、`vcp_bfem_view_dat.py` を使う |
| `pvpython was not found on PATH` | `python3-paraview` が入っていない。`.vtu` は生成済みなので GUI で開いてもよい |
| `cells must be non-negative integers` | 引数が逆。`cells` が先 |
| `point rows ... are not consistent` | `cells` と `points` が別のメッシュのもの。同じ実行で出た対を渡す |
| `... is missing: all three scripts must sit in the same folder` | 3 本が揃っていない |
| 絵が一瞬出て消える | 古い生成物 (`*_view.py`) を実行している。それは廃止された。`vcp_bfem_view_vtu.py` を使う |

---

## 6. 精度保証について

**表示される絵は精度保証ではありません。**

保証されるのは 2 段目までです。

| 段 | 保証 |
|---|---|
| `.dat` の値 | **あり**。区間型なら厳密値を包含する |
| `.dat` → `.vtu` | **あり**。区間の境界は double へ厳密に往復し、有理数は外向きに丸めて厳密値を挟む |
| ParaView での描画 | **なし** |

3 段目が保証できない理由は次の 3 つで、いずれも回避できません。

1. 中心や片方の境界を選んだ時点で包含という主張が消える
2. ParaView がセルの内部を線形補間する。補間値は元の有限要素関数の包含ではない
3. カラーマップの量子化とラスタライズで区間の幅は絵に残らない

区間型の入力では、値に加えて**半径**の配列 (`<名前>_rad`) も書き出します。
「どこで包含が広いか」は 2 段目までの保証で閉じているので、こちらは判断材料に
なります。
