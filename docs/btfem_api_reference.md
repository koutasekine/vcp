# btfem API リファレンス

この文書は `vcp/btfem` の利用者が直接触るクラス・関数をまとめたものです。
名前空間はすべて `vcp::btfem` です。内部実装用の `detail::` は通常使いません。
空間側の型(`mesh`, `fe_space`, `fe_function`, `poly_field`, `poly1`, `bpoly`,
`bary_point`, `range_pair`, `dirichlet_reduction`)は `vcp::bfem` のものをそのまま使います。

## この文書でのスコープ

| 表記 | どこの話か | 典型的な戻り値・状態 |
|---|---|---|
| 領域非依存 | メッシュから独立した添字・有理数・テーブル | `time_bc`、参照テーブル(有理数段) |
| 参照区間上 | 参照区間 `[0,1]` 上の時間 Bernstein 多項式。`k` はまだ掛からない | `Mref`, `Sref`, `Bref`(`detail`) |
| 参照単体上 | 標準単体上の空間 Bernstein 多項式(bfem) | `bpoly<D,T>` |
| 物理要素上 | 1 小区間(幅 `k`)× 1 空間要素に写した局所計算 | `time_space::local_*`、局所テンソル |
| 空間時間全体 | 全小区間・全要素と大域自由度を持つ空間・場・行列 | `pg_space`, `st_poly`, `kron_term`, `spmatrix`, 大域ベクトル |

迷った場合は、`pg_space` の `A` ... `Y` と `*_terms` は「空間時間全体の因子」、
`time_space` の `local_*` は「物理要素上(1 小区間)の局所テーブル」です。

## 時間メッシュ

### `time_mesh<T>`

スコープ: 空間時間全体の入力データです。`J_i = [t0, t1]` を `m` 等分します。

include:

```cpp
#include <vcp/btfem/time_mesh.hpp>
```

| API | 役割 |
|---|---|
| `static time_mesh uniform(t0, t1, m)` | `[t0, t1]` を `m >= 1` 等分。`k = (t1 - t0) * (1/m)`、`1/k = m / (t1 - t0)`(構築時の唯一の除算) |
| `t_begin() const`, `t_end() const` | `t0`, `t1` |
| `num_intervals() const` | `m` |
| `width() const` | 小区間幅 `k`(全小区間で共通) |
| `inv_width() const` | `1/k`。`k` とは独立に計算した値(区間型では `k * inv_k` は 1 を含むだけ) |
| `node(l) const` | `tau_l = t0 + l k`(`l = 0..m`)。乗算 1 回・加算 1 回 |
| `same_as(other) const` | 値としての同一メッシュ判定(`m` と両端。区間型は重なり判定) |

例外:

| 例外 | 条件 |
|---|---|
| `std::invalid_argument` | `m < 1`、`node(l)` の範囲外 |
| `degenerate_interval` | `t1 - t0` が 0・負、または区間型で符号が確定しない |

次の時間区間 `J_{i+1}` は利用者が `uniform(t1, t2, m')` を別に作ります。非等間隔分割はありません。

## 時間空間

### `time_bc`

```cpp
enum time_bc { time_bc_none, time_bc_zero_left };
```

`time_bc_none` は `Vbar_k^q`(左端に初期値を持つ、`ndof = m q + 1`)、
`time_bc_zero_left` は `V_k^q`(左端係数を除去、`ndof = m q`)です。

### `time_space<T>`

スコープ: 参照区間上の局所テーブルを物理化し、空間時間全体の時間自由度番号と大域時間行列を提供します。

include:

```cpp
#include <vcp/btfem/time_space.hpp>
```

```cpp
time_space<double> Vt(tm, 2, time_bc_zero_left);
```

局所 rank `r = 0..q` は bfem `index_map<1>(q)` の正準順序(`r = 0` が左端 `lambda_0^q`、
`r = q` が右端 `lambda_1^q`)です。大域自由度は `g(l, r) = (l-1) q + r`(`bc_none`)、
`(l-1) q + r - 1` で `(1, 0)` は `-1`(`bc_zero_left`)です。

| API | 返り値 | スコープ | 役割 |
|---|---|---|---|
| `time_space(tm, q, bc)` |  | 空間時間全体 | `q >= 1` |
| `degree()`, `bc()`, `num_intervals()`, `local_size()`, `ndof()` | `int` 等 | 空間時間全体 | `q`, `bc`, `m`, `q + 1`, 自由度数 |
| `global_dof(l, r)` | `int` | 空間時間全体 | `g(l, r)`。除去自由度は `-1` |
| `mesh_ref()` | `const time_mesh<T>&` | 空間時間全体 | 保持する時間メッシュ |
| `mass<P>()` | `matrix<T,P>` | 空間時間全体 | `U_t = (psi_b, psi_a)`(密、`ndof x ndof`) |
| `stiffness<P>()` | `matrix<T,P>` | 空間時間全体 | `A_t = (psi_b', psi_a')` |
| `dt_mass<P>()` | `matrix<T,P>` | 空間時間全体 | `B_t = (psi_b, psi_a')`。行 `a` が試験(微分)側 |
| `end_values<P>()` | `matrix<T,P>` | 空間時間全体 | `y_a = psi_a(t_end)`(`ndof x 1`) |
| `mixed_mass<P>(other)` | `matrix<T,P>` | 空間時間全体 | 行 = 本空間、列 = `other`。`(psi'_b, psi_a)` |
| `mixed_stiffness<P>(other)` | `matrix<T,P>` | 空間時間全体 | `(∂_t psi'_b, ∂_t psi_a)` |
| `mixed_dt_mass<P>(other)` | `matrix<T,P>` | 空間時間全体 | `(psi'_b, ∂_t psi_a)` |
| `local_mass()` | `const matrix<T,mats<T>>&` | 物理要素上 | `k · Mref(q,q)` |
| `local_stiffness()` | 同上 | 物理要素上 | `(1/k) · Sref(q)` |
| `local_dt_mass()` | 同上 | 物理要素上 | `Bref(q)`(`k` に依らない) |
| `local_dt_load(sdeg)` | 同上 | 物理要素上 | 行 `sigma`、列 `a`: `(beta^{sdeg}_sigma, psi_a')` |
| `local_weighted_dt_mass(sdeg, sigma)` | 同上 | 物理要素上 | `(beta^{sdeg}_sigma psi_b, psi_a')` |
| `local_weighted_mass(sdeg, sigma)` | 同上 | 物理要素上 | `k · (beta^{sdeg}_sigma psi_b, psi_a)` |
| `same_mesh(other)` | `bool` | 空間時間全体 | 同一 `time_mesh` 判定 |

`mixed_*` は同一 `time_mesh` 上の時間空間間でのみ定義され、異なるメッシュは `std::invalid_argument` です。
局所テーブルは参照区間上の有理数テーブルを `T` に一度だけ変換(プロセス共有キャッシュ)し、
`k` または `1/k` を乗算したものです。試験側 `psi_a'` は微分写像 `D^q` を乗じた形で保持されるので、
利用者が微分写像を意識する必要はありません。

恒等式(テストで確認): `V_k^q` で `dt_mass + dt_mass^T == y y^T`、`Vbar_k^q` で `== y y^T - e_0 e_0^T`。

## 空間時間関数

### `st_poly<D,T,P>`

スコープ: 空間時間全体です。`u(t, x) = Σ_a Σ_i C[a][i] psi_a(t) phi_i(x)` の係数行列
`C`(`ndof_t x ndof_x`、行 = 時間自由度、列 = 空間自由度)を持ちます。
行優先に並べたものが時間ブロック主ベクトル `v[a n_x + i]` です。

include:

```cpp
#include <vcp/btfem/st_poly.hpp>
```

| API | 役割 |
|---|---|
| `time_degree()`, `time_bc_kind()`, `space_degree()` | `q`, `bc`, `n` |
| `num_time_dofs()`, `num_space_dofs()` | `ndof_t`, `ndof_x` |
| `coeffs()` | 係数行列。直接編集可(`fe_function` と同じ) |

生成は `pg_space` のファクトリ経由のみです(既定構築子はなく、`std::array<st_poly, D>` は
波括弧初期化で作ります)。`st_poly` は作成元の `pg_space`(`fe_space` と時間メッシュのコピー)を記憶し、
別の `pg_space` に渡すと `std::invalid_argument` です。合成結果 `f(ubar)` の大域 `st_poly` は存在しません。

## Kronecker 因子と縮約

include:

```cpp
#include <vcp/btfem/kron.hpp>
```

### `kron_term<T,P,SP>`

スコープ: 空間時間全体の因子です。

```cpp
struct kron_term {
    T                      scale;   // 既定 T(1)。ν などはここに入れる
    vcp::matrix<T, P>      time;    // ndof_t x ndof_t'(密)
    vcp::spmatrix<T, SP>   space;   // n_x x n_x'
};
```

### 自由関数

| API | 返り値 | 役割 |
|---|---|---|
| `kron_expand(term)` | `spmatrix<T,SP>` | 時間ブロック主で展開。行 `a n_x + i`、列 `b n_x' + j`、値 `scale · time(a,b) · space(i,j)` |
| `kron_expand(terms)` | `spmatrix<T,SP>` | 項の和。全項の形が一致しないと `std::invalid_argument` |
| `kron_apply(terms, x)` | `matrix<T,P>` | 行列を作らず `y = Σ_t scale_t (time_t ⊗ space_t) x` を計算 |
| `block_reduce(R, nblocks, v)` | `matrix<T,P>` | 時間ブロック主ベクトルの各ブロックに `dirichlet_reduction R` を適用 |
| `block_expand(R, nblocks, vr)` | `matrix<T,P>` | 拘束成分に `T(0)` を戻す逆変換 |
| `block_reduce(R, nblocks, A)` | `spmatrix<T,SP>` | 展開済み行列の行・列をブロックごとに縮約 |

展開の加算順は「項 → 時間成分 `(a, b)` 行優先 → 空間成分の列挙順」で固定され、重複は bfem の
`coo_buffer::combine` で決定的に結合します。時間行列の厳密零(`== T(0)`。区間型では幅 0 の零のみ)は
読み飛ばします。返る `spmatrix` は `finalize` 未実施です。
`kron_apply` と `block_reduce` は `spmatrix` の成分列挙だけを使い、演算 API には依存しません。

因子レベルで縮約する場合は `kron_term::space` に `R.reduce(...)` を適用してから展開します。
両経路は厳密に一致します。

### `initial_reduction<T,P>`

スコープ: 空間時間全体です。`Vbar_k^q`(`bc_none`)添字と `V_k^q`(`bc_zero_left`)添字の
相互変換(左端ブロックの除去・挿入)を行います。

| API | 役割 |
|---|---|
| `initial_reduction(full, reduced, n_x)` | `full` は `bc_none`、`reduced` は同一次数・同一メッシュの `bc_zero_left` |
| `full_size()`, `reduced_size()`, `space_size()` | `(m q + 1) n_x`, `m q n_x`, `n_x` |
| `reduce(v_full)` | 先頭 `n_x` 成分(時間自由度 0 のブロック)を除去 |
| `expand<D>(v_reduced, u_init)` | 先頭ブロックに `u_init.coeffs()` を置いて `Vbar` 添字へ戻す |

## 全離散 Petrov–Galerkin 空間

### `pg_space<D,T,P,SP>`

スコープ: 空間時間全体です。bfem の `fe_space`(所有しない。`pg_space` より長生きさせる)と
そのメッシュ、時間メッシュ(コピー)から構築します。

include:

```cpp
#include <vcp/btfem/pg_space.hpp>
```

```cpp
vcp::bfem::fe_space<2, double> Vh(msh, 1);
vcp::btfem::time_mesh<double> tm = vcp::btfem::time_mesh<double>::uniform(0.0, 1.0, 4);
vcp::btfem::pg_space<2, double> S(Vh, msh, tm);
auto A = S.A(1, 2);                           // A_t(q = 1) ⊗ U_x(n = 2)
vcp::spmatrix<double> K = vcp::btfem::kron_expand(A);
```

試行・試験の次数 `(q, n)` は毎回明示します(既定値なし)。`ubar` の次数 `(qbar, nbar)` と `(q, n)` は
独立で、`A_mixed` / `B_mixed` だけが `nbar == n` を要求します。空間因子は全自由度系(境界条件未適用)で、
Dirichlet 縮約は `kron.hpp` の経路で行います。生成系メソッドは内部バッファを変更するため非 `const` です。

| API | 返り値 | スコープ | 役割 |
|---|---|---|---|
| `pg_space(Vh, msh, tm)` |  | 空間時間全体 | `msh.num_elements() == Vh.num_elements()` を検査 |
| `space()` | `fe_space&` | 空間時間全体 | 参照先 bfem 空間 |
| `mesh_ref()`, `time_mesh_ref()` | 参照 | 空間時間全体 | 構築時のメッシュ・時間メッシュ |
| `time_dofs(q, bc)` | `const time_space<T>&` | 空間時間全体 | `(q, bc)` でキャッシュされた時間空間 |
| `zero_function(q, bc, n)` | `st_poly` | 空間時間全体 | 全係数 `T(0)` |
| `function_from_coeffs(q, bc, n, C)` | `st_poly` | 空間時間全体 | 係数行列から。サイズ検査 |
| `function_from_vector(q, bc, n, v)` / `to_vector(u)` | `st_poly` / `matrix<T,P>` | 空間時間全体 | 時間ブロック主ベクトルとの相互変換 |
| `A(q,n)`, `M(q,n)`, `B(q,n)`, `U(q,n)`, `W(q,n)`, `Y(q,n)` | `kron_term` | 空間時間全体 | `A_t⊗U_x`, `U_t⊗K_x`, `B_t⊗K_x`, `U_t⊗U_x`, `A_t⊗K_x`, `y yᵀ⊗U_x`(試行・試験とも `V_k^q⊗S_h^n`、`scale = T(1)`) |
| `A_mixed(q,n,ubar)`, `B_mixed(q,n,ubar)` | `kron_term` | 空間時間全体 | 行 = `V_k^q`、列 = `ubar` の時間空間。`(psibar_b', psi_a')⊗U_x`, `(psibar_b, psi_a')⊗K_x` |
| `weighted_dt_mass_terms(fprime, ubar, q, n)` | `std::vector<kron_term>` | 空間時間全体 | `(f'[ubar] phi_j, ∂_t phi_i)` の因子列。項数 `m (deg f' · qbar + 1)`、`(l, sigma)` 順 |
| `convection_dt_terms(b_const, q, n)` | `std::vector<kron_term>` | 空間時間全体 | 単一項 `B_t ⊗ C_x(b)` |
| `convection_dt_terms(b, q, n)` | `std::vector<kron_term>` | 空間時間全体 | `b` は `std::array<st_poly, D>`(全成分が同じ `(q_b, bc)`)。項数 `m (q_b + 1)` |
| `load_dt(g, ubar, q, n)` | `matrix<T,P>` | 空間時間全体 | `(g(ubar), ∂_t phi_i)` |
| `convection_dt(b_const, ubar, q, n)` / `convection_dt(b, ubar, q, n)` | `matrix<T,P>` | 空間時間全体 | `(b·∇ubar, ∂_t phi_i)` |
| `section_at_node(u, l)` | `fe_function` | 空間時間全体 → 空間断面 | 時刻 `tau_l`(`l = 0..m`)の断面。ブロック切り出しのみ。`bc_zero_left` で `l = 0` は零関数 |
| `section_at(u, l, tau)` | `fe_function` | 空間時間全体 → 空間断面 | 小区間 `l`(`1..m`)の局所重心座標 `tau = (1-t, t)` での断面(de Casteljau)。`tau` が幅を持てば値の包含 |
| `time_trace(u, e, lam)` | `std::vector<bpoly<1,T>>` | 空間時間全体 → 時間方向 | 空間点 `(e, lam)` での小区間ごとの時間多項式 |
| `eval(u, l, tau, e, lam)` | `T` | 点評価 | 小区間 `l` の `tau` と要素 `e` の `lam` で評価 |
| `range(u)` | `range_pair<T>` | 空間時間全体 | 全係数の凸包(値域の包含) |
| `range_composed(f, ubar)` | `range_pair<T>` | 空間時間全体 | `f(ubar)` の局所テンソル係数の凸包 |
| `elevate_space(u, n2)` | `st_poly` | 空間時間全体 | 各時間ブロックを bfem `elevate` で次数 `n2` へ |
| `initial_error_sq_l2(u0, ubar)` | `T` | 空間時間全体 | `‖u_0 - ubar(t_0)‖²_{L²}`。`u0` は `fe_function` または `poly_field` |
| `initial_error_sq_h01(u0, ubar)` | `T` | 空間時間全体 | `‖∇(u_0 - ubar(t_0))‖²_{L²}` |

線形部 `(∂_t ubar, ∂_t phi_i) + ν(∇ubar, ∇∂_t phi_i)` は、`ubar` が `V_k^q` なら
`kron_apply({A(q,n)}, to_vector(ubar))` と `kron_apply({B(q,n)}, ...)`、`ubar` が `bc_none` なら
`A_mixed` / `B_mixed` で得ます。

物理時刻 `t` から `(l, tau)` を求める処理は本層にありません(bfem `eval(u, e, lam)` と同じ流儀)。

### 例外

| 例外 | 条件 |
|---|---|
| `degenerate_interval` | `time_mesh::uniform` の区間長が符号不定・零 |
| `std::invalid_argument` | `m < 1`, `q < 1`, `n < 1`、サイズ・次数の事前条件違反、別 `pg_space` の `st_poly`、`A_mixed` / `B_mixed` の `nbar != n`、`kron_expand` の項間形不一致、`initial_reduction` の `full` / `reduced` 不整合、要素数不一致 |
| bfem の例外 | 空間側の退化・不整合はそのまま伝播 |

### 区間型・有理数型での意味

有理数型では全 API が厳密値を返します(テストで 2 経路・製造解との厳密一致を確認)。
区間型では、各成分が対応する厳密値を包含し、`range` / `range_composed` の `lo`, `hi` は
それぞれ真の最小・最大の包含です。`section_at` の `tau` に幅を持たせると、その時間範囲の
値の包含が返ります。

## 使い方の骨子(`spnewton` 向け)

1. `fe_space` と `mesh`、`time_mesh` から `pg_space` を構築する。
2. 未知数は `V_k^{qbar} ⊗ S_h^{nbar}` の係数(左端固定)。`initial_reduction::expand` で
   `Vbar` の `ubar` を復元し、`load_dt` / `convection_dt` と `kron_apply(A_mixed / B_mixed)` で残差を組む。
3. Jacobian は `A`, `B`(`scale = ν`), `weighted_dt_mass_terms`, `convection_dt_terms` を
   `kron_expand` した正方行列。空間 Dirichlet は `block_reduce` または因子レベルの `R.reduce`。
4. 時間発展は `section_at_node(ubar, m)` を次区間の初期値として `uniform(t1, t2, m')` で続ける。

## BT-2: RT ⊗ 時間と broken ⊗ 時間

空間側の型(`rt_space`, `rt_field`, `broken_space`, `broken_field`, `rt_element_op`)は
`vcp::bfem` のものをそのまま使います。時間空間 `(q, bc)` は構築時ではなく、
場や行列を作るたびに明示します(既定値・自動推定なし)。各空間は `(q, bc)` ごとの
`time_space` を自前でキャッシュします(`pg_space` と同じ値、共有はしない)。

### `st_rt_field<D,T,P>`

スコープ: 空間時間全体です。`σ(t, x) = Σ_a Σ_j C[a][j] psi_a(t) tau_j(x)` の係数行列
`C`(`ndof_t × ndof_rt`、行 = 時間自由度、列 = RT 大域自由度。符号は bfem `rt_dofmap` の規約)を持ちます。

| API | 役割 |
|---|---|
| `order()`, `time_degree()`, `time_bc_kind()` | RT 次数 `l`、`q_σ`、`bc_σ` |
| `num_time_dofs()`, `num_space_dofs()` | `m_t`、`rt_space::ndof()` |
| `coeffs()` | 係数行列。直接編集可 |

生成は `st_rt_space` のファクトリ経由のみで、作成元(`rt_space` と時間メッシュのコピー)を記憶します。

### `st_rt_space<D,T,P,SP>`

スコープ: 空間時間全体です。bfem の `rt_space`(所有しない。`st_rt_space` より長生きさせる)と時間メッシュ(コピー)から構築します。

include:

```cpp
#include <vcp/btfem/st_rt_space.hpp>
```

```cpp
vcp::bfem::rt_space<2, double> rs(msh, 1);
vcp::btfem::st_rt_space<2, double> srs(rs, tm);
auto sig = srs.zero_field(2, vcp::btfem::time_bc_none);      // Vbar_k^2 ⊗ RT^1
auto M = srs.mass(2, vcp::btfem::time_bc_none);             // U_t ⊗ M_RT
```

| API | 返り値 | スコープ | 役割 |
|---|---|---|---|
| `st_rt_space(rs, tm)` |  | 空間時間全体 | |
| `space()` | `rt_space&` | 空間時間全体 | 参照先 bfem RT 空間 |
| `time_mesh_ref()` | `const time_mesh<T>&` | 空間時間全体 | 時間メッシュ(コピー) |
| `time_dofs(q, bc)` | `const time_space<T>&` | 空間時間全体 | `(q, bc)` でキャッシュされた時間空間 |
| `zero_field(q, bc)` | `st_rt_field` | 空間時間全体 | 全係数 `T(0)` |
| `field_from_coeffs(q, bc, C)` | `st_rt_field` | 空間時間全体 | 係数行列から。サイズ検査 |
| `field_from_vector(q, bc, v)` / `to_vector(s)` | `st_rt_field` / `matrix<T,P>` | 空間時間全体 | 時間ブロック主ベクトル `(a, j) ↦ a·ndof_rt + j` との相互変換(符号不変) |
| `section_at_node(s, l)` | `rt_field` | 空間断面 | 時刻 `tau_l`(`l = 0..m`)の RT 場(ブロック切り出し)。`bc_zero_left` で `l = 0` は零場 |
| `section_at(s, l, tau)` | `rt_field` | 空間断面 | 小区間 `l` の局所重心座標 `tau = (1-t, t)` での RT 場(de Casteljau) |
| `mass(q, bc)` | `kron_term` | 空間時間全体 | S3: `(σ_j, τ_i)_{L²(J;L²)} = U_t(q) ⊗ M_RT`、`scale = T(1)` |
| `validate(s)` | `void` | 空間時間全体 | 場の所属・次数・サイズ検査(`std::invalid_argument`)。自由関数も使う |

### `st_broken_field<D,T,P>` / `st_broken_space<D,T,P,SP>`

スコープ: 空間時間全体です。`st_rt_*` と同形で、列は bfem `broken_dofmap` の要素ブロック順
`e·nloc + r`、符号は `+1` です。

include:

```cpp
#include <vcp/btfem/st_broken_space.hpp>
```

| API | 返り値 | 役割 |
|---|---|---|
| `st_broken_space(bs, tm)` |  | bfem `broken_space`(所有しない)と時間メッシュ(コピー) |
| `space()`, `time_mesh_ref()`, `time_dofs(q, bc)` | 参照 | `st_rt_space` と同じ |
| `dofmap_view()` | `detail::st_broken_dofmap_view` | BT-1 `st_gather` に渡せる broken 自由度写像(`degree()` 付き) |
| `zero_field(q, bc)`, `field_from_coeffs(q, bc, C)`, `field_from_vector(q, bc, v)`, `to_vector(u)` | | `st_rt_space` と同じ |
| `section_at_node(u, l)`, `section_at(u, l, tau)` | `broken_field` | 断面(意味論は `pg_space` と同じ) |
| `mass(q, bc)` | `kron_term` | `(q_j, q_i)_{L²(J;L²)} = U_t(q) ⊗ M_broken`(ブロック対角) |
| `inner(u, v)` | `T` | `(u, v)_{L²(J;L²)}`。`u`, `v` の `(q, bc)` は異なってよい。加算順は小区間 → 要素 → `(r, r')` |
| `validate(u)` | `void` | 場の所属・次数・サイズ検査 |

### `st_rt_assemble.hpp` の自由関数

スコープ: 空間時間全体です。bfem `rt_assemble.hpp` の 7 関数に 1 対 1 で対応し、
引数の並びも bfem と同じ(行側の空間 → 列側の空間 → 補助引数)です。
`pg` 側の試行空間は `V_k^q`(`bc_zero_left`)固定で `(q, n)` を取ります(E-5)。
`ū` は `pg_space` の `st_poly`(`bc_none` / `bc_zero_left` どちらでも)、
`g` は `poly1<T>`、`b` は `std::array<T, D>`(定数)または `std::array<st_poly, D>`(空間時間、
全成分が同じ `(q_b, bc)`)で、`b` なしの overload が bfem の形に対応します。
`w = g(ū) − ∂_t ū − b·∇ū` です。

include:

```cpp
#include <vcp/btfem/st_rt_assemble.hpp>
```

| API | 返り値 | 役割 |
|---|---|---|
| `st_div_field(srs, σ, sbs)` | `st_broken_field` | `div σ` を broken ⊗ 時間(σ と同じ `(q, bc)`)の場として返す。`sbs.space().order() == srs.space().order()` を要求 |
| `st_assemble_div_mass(sbs, qb, bcb, srs, qs, bcs)` | `kron_term` | S4 (K3): 行 = broken ⊗ `(qb, bcb)`、列 = RT ⊗ `(qs, bcs)`、`(div σ_j, q_i)_{L²(J;L²)}` |
| `st_assemble_cross_grad(srs, qs, bcs, pg, q, n)` | `kron_term` | S5 (K4): 行 = RT ⊗ `(qs, bcs)`、列 = `V_k^q ⊗ S_h^n`、`(σ_i, ∇φ_j)` |
| `st_assemble_mixed_mass(pg, q, n, sbs, qb, bcb)` | `kron_term` | S6 (K9): 行 = `V_k^q ⊗ S_h^n`、列 = broken ⊗ `(qb, bcb)`、`(φ_i, q_j)` |
| `st_flux_error_sq(srs, σ, pg, ū)` | `T` | S7 (K5): `‖∇ū − σ‖²_{L²(J;L²)}`(ν は API に持たない、E-1)。`n̄ ≥ 1` |
| `st_div_residual_sq(srs, σ, g[, b], pg, ū)` | `T` | S8 (K5): `‖div σ + w‖²_{L²(J;L²)}` |
| `st_broken_load(sbs, qb, bcb, g[, b], pg, ū)` | `matrix<T,P>` | S9 (K7): `(w, q_i)_{L²(J;L²)}`、行 = broken ⊗ `(qb, bcb)`(時間ブロック主) |
| `st_projection_error_sq(sbs, g[, b], pg, ū)` | `T` | S10 (K8): `‖(I − Π_l) w‖²_{L²(J;L²)}`、`Π_l` は空間方向の要素ごと `L²` 射影(時間方向は射影しない) |

因子(S3〜S6)の `kron_term` は `scale = T(1)`、`time` は行側と列側の時間空間の
混合質量(長方形)、`space` は bfem のアセンブラの結果です。展開は `kron_expand` で、
行 `a·n_x + i`、列 `b·n_x' + j` です。

### bfem ↔ btfem 対応表

| bfem(空間) | btfem(空間時間) | 返り値 |
|---|---|---|
| `rt_space`, `rt_field` | `st_rt_space`, `st_rt_field` | — |
| `broken_space`, `broken_field` | `st_broken_space`, `st_broken_field` | — |
| `rt_space::mass()` | `st_rt_space::mass(q, bc)` | `kron_term` |
| `broken_space::mass()` / `inner(u,v)` | `st_broken_space::mass(q, bc)` / `inner(u,v)` | `kron_term` / `T` |
| (RT 版 `div_field` は bfem に無い) | `st_div_field(srs, σ, sbs)` | `st_broken_field` |
| `assemble_div_mass(bs, rs)` | `st_assemble_div_mass(sbs, qb, bcb, srs, qs, bcs)` | `kron_term` |
| `assemble_cross_grad(rs, fs, m)` | `st_assemble_cross_grad(srs, qs, bcs, pg, q, n)` | `kron_term` |
| `assemble_mixed_mass(fs, m, bs)` | `st_assemble_mixed_mass(pg, q, n, sbs, qb, bcb)` | `kron_term` |
| `flux_error_sq(rs, σ, fs, u)` | `st_flux_error_sq(srs, σ, pg, ū)` | `T` |
| `div_residual_sq(rs, σ, f, fs, u)` | `st_div_residual_sq(srs, σ, g[, b], pg, ū)` | `T` |
| `broken_load(bs, f, fs, u)` | `st_broken_load(sbs, qb, bcb, g[, b], pg, ū)` | `matrix` |
| `projection_error_sq(bs, f, fs, u)` | `st_projection_error_sq(sbs, g[, b], pg, ū)` | `T` |

bfem の `f(u)` に対し btfem は `w = g(ū) − ∂_t ū − b·∇ū` を取る点だけが対応からずれます。
ū と σ が時間定数なら、各量は `|J| ×` bfem の対応量に一致します(テスト Q-R7)。

### 例外(BT-2 追加分)

| 例外 | 条件 |
|---|---|
| `std::invalid_argument` | 空間側のメッシュ不一致(bfem `rt_require_same_mesh` 伝播)、時間メッシュ不一致(`time_mesh::same_as` 偽)、別の `st_rt_space` / `st_broken_space` の場、場の次数・サイズ不一致、`st_div_field` の broken 次数 ≠ RT 次数、`q < 1`、`ū` の形が `pg_space` と合わない、`b` の成分の時間空間不一致、S7 で `n̄ < 1` |
| bfem の例外 | 空間側の退化・不整合はそのまま伝播 |

### 区間型・有理数型での意味(BT-2)

有理数型では全 API が厳密値を返します(テスト G1〜G8 で 2 経路・独立経路との厳密一致を確認)。
区間型では各成分・各スカラーが対応する厳密値を包含します(テスト G9)。
`S10` の `M_l⁻¹` は有理数段で逆行列化した `T` 像(bfem `typed_rt_registry::inv_mass`)で、
BT-2 に新しい除算はありません。

## BT-3: 定数層 `vcp/btfem/constants/`

名前空間 `vcp::btfem::constants`。`T` は kv の区間型(`kv::interval<double>`、`kv::interval<kv::dd>`、
`kv::interval<kv::mpfr<N>>`)に限定され、`DP` / `SP` の既定は bfem 定数層と同じ
`vcp::imats<typename T::base_type>` / `vcp::spimats<typename T::base_type>` です。

### 対象の式(定数論文 定理 1〜3、補題 3、u₀ = û₀ = 0、k = q、C_τ := C_τ(q − 1))

| 量 | 式 | 不等式 |
|---|---|---|
| K = C̃₁² | (C_Ω² + ν² λ_h C_τ² + 4ν C_τ) / ν² | ‖u − u_{h,τ}^k‖²_{L²H¹₀} ≤ K ‖f‖² |
| L = c̃₀² | 4 (C_Ω² / ν + C_τ) | ‖(u − u_{h,τ}^k)(b)‖²_{L²} ≤ L ‖f‖² |
| C̃₀ | 4 C_τ + 2 C_Ω C̃₁ | ‖u − u_{h,τ}^k‖_{L²L²} ≤ C̃₀ ‖f‖ |

### 関数表

| 関数 | ヘッダ | 返り値 | スコープ |
|---|---|---|---|
| `time_projection_error_constant<T>(int d, const T& tau)` | `time_projection_constants.hpp` | C_τ(d) = c_d·τ | 参照区間上(τ でスケール) |
| `inverse_inequality_h01_bound<D, T, DP, SP>(const bfem::mesh<D,T>& Th, int n)` | `inverse_inequality.hpp` | λ_h の上界(点区間 `T(λ_up)`。下端は下界ではない) | メッシュ全体・斉次 Dirichlet |
| `heat_pg_space_materials_h01<D, T, DP, SP>(Th, n, mesh_resolution_report* rep = nullptr)` | `heat_pg_constants.hpp` | `heat_pg_space_material_set_h01<T>` {`n`, `c_omega_sq`, `lambda_h`}(C_Ω² を先に、λ_h を後に計算) | 同上 |
| `heat_pg_error_constants_h01<T>(c_omega_sq, c_tau, lambda_h, nu)` | 同上 | `heat_pg_error_constant_set_h01<T>`(材料 4 + K, L, c1_tilde, c0_small, c0_tilde)。式のみ、材料は差し替え可 | 同上 |
| `heat_pg_error_constants_h01<T>(const heat_pg_space_material_set_h01<T>& sm, const T& tau, int q, const T& nu)` | 同上 | 同上(c_tau = `time_projection_error_constant(q − 1, tau)`) | 同上 |
| `heat_pg_error_constants_h01<T>(sm, const time_mesh<T>& tm, int q, const T& nu)` | 同上 | 同上(tau = `tm.width()`) | 同上 |
| `heat_pg_error_constants_h01<D, T, DP, SP>(Th, n, tau, q, nu, rep = nullptr)` | 同上 | 同上(材料構造体を内部で作る一括版。材料版とビット一致) | 同上 |

計算順は K → L → √K, √L → C̃₀ に固定され、平方根は √K、√L、√C_Ω² の 3 箇所だけです。
薄い射影入口(1 成分だけ返す関数)はなく、取り出しは構造体メンバで行います。

### 出典表(C_τ(d) = c_d·τ)

| d = q − 1 | c_d | 出典 | 種別 |
|---|---|---|---|
| 0 | 1/π | Nakao–Yamamoto–Kimura, J. Approx. Theory 93 (1998) 1 章(Schultz, Spline Analysis 定理 2.5 (2.17)) | 最良値 |
| 1 | 1/(2π) | NYK 1998 定理 1 | 最良値 |
| 2 | (1/(3√5π³))·√(−2160 + 25π⁴ + 4√5·√(−174960 − 1080π⁴ + 243π⁶ + 5π⁸)) ≈ 1/8.92337 | NYK 1998 定理 2 (14) | 上界 |
| ≥ 3 | √max(L1, L2)、L1 = 1/(2(2d+1)(2d+5)) + 1/(4(2d+5)√(2d+3)√(2d+7))、L2 = 1/(4(2d+5)√(2d+3)√(2d+7)) + 1/(2(2d+5)(2d+9)) + 1/(4(2d+9)√(2d+7)√(2d+11)) | Kimura–Yamamoto, Bull. Inform. Cybernet. 31 (1999) 定理 2.3 (2.3)、N = d | 上界 |

添字の対応: KY の N(Q_N の次数は N + 1)は d、`vcp/ldbase.hpp` の `Ritz_projection_error(N)` の N も d
(項順も同じ)、Nakao–Kinoshita (2009) 定理 2 の N(P_N^{1,0} の次数は N)は d + 1(テストの交差検査にのみ使用)。

**帰着 [自己導出、有理数ゲート G1 で機械確認]**: 小区間 [0, τ] 上で v ∈ H¹ に対し V(t) := ∫₀ᵗ v、
R を「次数 d + 1、両端で V を補間し、(V − RV)′ ⊥ {q′ : q ∈ P_{d+1}, q(0) = q(τ) = 0}」で定まる
H¹₀ 型射影とすると (V − RV)′ = (I − P_d)v であり、u := V − (両端の線形補間) は H¹₀ ∩ H² に入って
u − π_{d+1}u = V − RV、|u|_{H²} = ‖v′‖ となる。よって次数 d + 1 の 1 次元 H¹₀ 射影の誤差定数に τ を
掛けたものが C_τ(d) の上界になる。

λ_h の経路は bfem `hypercircle_kappa_h01_squared`(`assemble_dense_materials` → `hypercircle_kappa_core`)と
同じ呼び出し列(`fe_space` → `dirichlet_reduction` → 疎 `finalize` → 縮約 → 密化 → `compsym` → `eigsymge` →
`diag` → 上端の最大値)で、新しい密化・対称化・固有値経路はありません。参照値: 参照三角形 36、参照四面体 80
(n = 1、縮約なし。有理数で特性多項式の根かつ最大固有値であることを確認済み)。

### 例外(BT-3 追加分)

| 例外 | 条件 |
|---|---|
| `vcp::invalid_argument` | `d < 0`、`tau` の下端 ≤ 0、`n < 1`、縮約後の自由度 0(λ_h が定義されない)、`q < 1`、`nu` の下端 ≤ 0、`c_tau` / `lambda_h` / `c_omega_sq` の下端 < 0 |
| `vcp::verification_error` | d = 2 の √ の引数の下端が正でない(転記誤りの検出)、`compsym` の交わりが空、`eigsymge` の証明失敗、上端の最大値が負 |
| `vcp::bfem::degenerate_element` | 退化要素(空間構築から伝播、`vcp::error` 階層外) |
| bfem の例外 | `ritz_projection_error_constant_h01_sq` の例外(3D 辞書未登録を含む)はそのまま伝播 |

### 区間型での意味(BT-3)

C_τ は閉じた式の区間評価で、`kv::interval<double>` と `kv::interval<kv::mpfr<N>>` では真値の包含です。
λ_h と bfem の κ_h² は上端だけに意味がある点区間なので、それらを材料に含む K、L、C̃₁、c̃₀、C̃₀ も
上端(上界)として読みます。同一入力・単スレッド・同一 DP でビット一致します(DP の違いによる
ビット差は契約外、bfem CM-1T と同じ)。
