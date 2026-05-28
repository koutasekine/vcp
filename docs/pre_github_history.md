# Pre-GitHub History

この文書は、GitHub で履歴管理を始める前に使っていた旧バージョン管理メモを
まとめたものです。現在の変更履歴は GitHub の commit history や release を
参照してください。

## vera0.1.1

fourier_basis.hpp: 1次元フーリエ基底を張るためのクラス
fourier_series.hpp: フーリエ級数の演算クラス
sincos_class.hpp: fourier_series.hppの下請け
ldbase.hpp: 拡張: 区間係数を持つ解に対しても最大値・最小値を求められるように，funcメソッドをイコール代入からconvert代入に変更，Mode=1でも利用可能にした.
minimatrix.hpp: psa型などの符号判定を持たない型用の行列ポリシーminimats< _T >.matrixクラスの部分特殊化.

bugfix: vera0.1.0においてmatrix_assist.hppにてkvのintervalを要求していたバグを解消
その他: test_PDEのファイルを論文の結果に変更

test_Emden_Numerische_Mathematik2020.cpp:
Kouta Sekine, Mitsuhiro T. Nakao, and Shin'ichi Oishi: “A new formulation using the Schur complement for the numerical existence proof of solutions to elliptic problems: without direct estimation for an inverse of the linearized operator”, Numer. Math., 146, 907-926, Oct., 2020. doi.org/10.1007/s00211-020-01155-7.

test_Lotka_Volterra2_NOLTA2021.cpp:
Kouta Sekine, Mitsuhiro T. Nakao, and Shin'ichi Oishi: “Numerical verification methods for a system of elliptic PDEs, and their software library”, to appear in NOLTA, IEICE, Vol.12, No.1, Jan., 2021.

## vera0.1.0

mats.hpp: normtwoメソッドを追加.2ノルムの計算をする.全固有値求めている.multmmとeigsymを使って実装しているため，BLASを使うポリシでは自動で高速に.(normone,norminfと違う実装法)
matrix.hpp: normtwo関数を追加.
imats_assist.hpp: 区間行列に対するmag関数を追加.mag(A, B), Bはinterval型ではない.
imats_assist.hpp: 区間行列に対するintervalmag関数を追加.B = intervalmag(A), Bはinterval型.
ldbase.hpp: 拡張: 区間係数を持つ解に対しても最大値・最小値を求められるように，funcメソッドをイコール代入からconvert代入に変更，Mode=1でも利用可能にした.
minimatrix.hpp: psa型などの符号判定を持たない型用の行列ポリシーminimats< _T >.matrixクラスの部分特殊化.

irk_legendre_parameter.hpp: Implicit Runge-Kuttaのルジャンドル法のブッチャー表を作成するクラス
implicit_rk.hpp: Implicit Runge-Kuttaをブッチャー表が与えられたら実行するクラス．デフォルトはルジャンドル法．
test_IRK_Legendre_approximate.cpp: implicit_rk.hppの使用例1．単純な2変数のODEの近似解を求める．
test_parabolic_Emden_approximate.cpp: ldbase.hppとimplicit_rk.hppを使用してEmden方程式の初期値-境界値問題の近似解を求める．

## vera0.0.6

mats.hpp: addmmなどのsubsmmA, subsmmBの行列サイズの判定条件のミスを修正(Bug fix)
mats.hpp: 表示関数operator<<のconstつけ忘れを修正(Bug fix)
matrix.hpp: 表示関数operator<<のconstつけ忘れを修正(Bug fix)
imats.hpp: linearsolveの無駄な計算が行われていた箇所を修正(Bug fix)
imats.hpp: linearsolveにR(Ax - b)の計算ミスを修正(Bug fix)
imats.hpp: linearsolveのメモリの節約(Rとmbのclear)
imats.hpp: eigsymgeの精度保証アルゴリズムを変更
pidblas.hpp: 内部で使用するmul_m_imとmul_im_mにlCとuCのミスを修正(Bug fix)
legendre_integral.hpp: Legendre積分の分点と重みのクラスinterval_ld_weightpoint(コードの整理)
ldbase.hpp: Legendre積分をlegendre_integral.hppに移動(コードの整理)
ldbase.hpp: LegendrePointFuncを並列化(コード最適化)
newton.hpp: Newton法のクラス．クラスを継承してメソッドfとdfをオーバーライドで使用
test_newton1.cpp: newton.hppの使用例1．単純な1変数の非線形方程式の近似解を求める．
test_newton2_Emden.cpp: newton.hppの使用例2．定常Emden方程式の近似解をnewton.hppで求める．

## vera0.0.5

ldbase.hpp: 残差モードを追加，Ritz射影の誤差定数の追加
test_Emden.cpp: Ritz射影の誤差評価の使い方を追加，残差とLipschitz定数の計算を追加, Newton-Kantorovichの定理による精度保証の追加
test_Gray_Scott.cpp: 残差とLipschitz定数の計算を追加
pidblas.hpp: kvのhwround.hppをinclude(kvのAVX512機能利用時のため)
imats_assist.hpp: 区間行列の強制対称化関数compsymを追加
matrix_assist.hpp: vcp_metafunction.cppをinclude

## vera0.0.4

test_Emden.cpp : Newton法の終了後の近似解uhをGeneratorに差し戻し忘れ．最大値，最小値に影響(Bug fixes)
test_Gray_Scott.cpp : Newton法の終了後の近似解uhをGeneratorに差し戻し忘れ．最大値，最小値に影響(Bug fixes)
vcp_timer.hpp : timer用のクラスを作成．C++11のchronoを使用．

## vera0.0.3

mats : rand関数の自動OpenMP化を無効
vcp_metafunc.hpp : 新規作成．is_round_control，is_interval，is_round_controlなどのメタ関数を追加
License : BSD 3-clause "New" or "Revised" Licenseを付加
