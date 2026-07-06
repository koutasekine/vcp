// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// vcp/tsparse/tsparse_honest_termination.hpp
//
// EIG-1 F-3/F-4/F-5: honest termination の唯一のコア実装(EIG-0 C-1/C-2)。
// TRL・公開 lanczos・shift-invert 経路はこのコアを共用する(コピー禁止)。
// tsparse_restart.hpp には ritz_pair/locked_pair 版の thin overload があり、
// 本ヘッダのコアに委譲する。
//
// 依存は tsparse_scalar / tsparse_eigs / tsparse_eigen_selection のみ
// (tsparse_lanczos.hpp から include しても spmatrix.hpp への循環を作らない)。

#pragma once

#ifndef VCP_TSPARSE_HONEST_TERMINATION_HPP
#define VCP_TSPARSE_HONEST_TERMINATION_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>

namespace vcp {
namespace tsparse {

// ---------------------------------------------------------------------------
// honest_termination_check_   (EIG-1 F-3-1 / F-4; EIG-0 C-2 の実装コア)
//
// converged 宣言の直前に呼ぶ、唯一の honest stopping rule 実装。
//
// 返却予定集合(locked_values の target 順 prefix k)の最悪値より target 側
// (内側)に、未収束の active 候補が「certainly に」存在する場合 false を返す
// (終了保留)。それ以外は true(見えている矛盾はない — 正直に終了してよい)。
//
// 比較規約(GT1 P1 / EIG-1 B-16): T の素の比較演算は interval では certainly
// 比較である。継続(false)側に倒すのは「内側と保証できる未収束候補がある」
// 場合のみで、indeterminate は終了許可(true)側に落ちる。converged の最終宣言
// 自体は C-1 の終了時厳密残差検査(residual_acceptance_check_, certified-≤)が
// 別途守るため、これで嘘側に倒れることはない。
//
// 注意: これは EIG-0 C-4 の通り完全性の保証ではない。「現在の部分空間に見えて
// いる矛盾を無視しない」という規則である。呼び出し側は freshness(候補集合が
// 現在の locked prefix を知る部分空間で計算されたこと)を別途保証すること。
// ---------------------------------------------------------------------------
template <class T>
bool honest_termination_check_(
	const std::vector<T>& active_values,
	const std::vector<bool>& active_converged,
	const std::vector<T>& locked_values,
	const std::size_t k,
	const eig_target target,
	const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
	typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
	typedef std::complex<real_type> complex_type;

	if (k == 0) return true;
	// 返却数が k に満たないままの converged 宣言は正直でない(呼び出し側の
	// 数え間違いに対する防波堤)。
	if (locked_values.size() < k) return false;

	// 返却予定 prefix k の最悪値: locked を target 順に並べた k 番目。
	// 順序は返却時と同じ選択関数(select_eigen_indices_from_real)を使う。
	const std::vector<std::size_t> order =
		vcp::tsparse_eigen_selection::select_eigen_indices_from_real(
			locked_values, locked_values.size(), target, shift);
	if (order.size() < k) return false;
	const T worst_returned = locked_values[order[k - 1]];

	const real_type worst_key = vcp::tsparse_eigen_selection::target_distance(
		complex_type(vcp::tsparse_scalar::real_part(worst_returned), real_type(0)),
		target, shift);
	const bool prefer_large =
		(target == eig_target::largest_magnitude ||
		 target == eig_target::largest_algebraic);

	const std::size_t n_active =
		(active_values.size() < active_converged.size())
			? active_values.size() : active_converged.size();
	for (std::size_t i = 0; i < n_active; i++) {
		if (active_converged[i]) continue;   // 収束済み候補はロック側で扱う(F-1)
		const real_type key = vcp::tsparse_eigen_selection::target_distance(
			complex_type(vcp::tsparse_scalar::real_part(active_values[i]), real_type(0)),
			target, shift);
		// certainly-inner のときのみ継続側(false)に倒す
		const bool certainly_inner = prefer_large ? (key > worst_key)
		                                          : (key < worst_key);
		if (certainly_inner) return false;
	}
	return true;
}

// ---------------------------------------------------------------------------
// honest_termination_check_complex_   (EIG-3 D3-3。B-22 例外 (ii) の「追加」)
//
// 複素候補対応の C-2 検査。実装は上の real 版と同一の certainly-inner 規則で、
// 候補の target 計量のみ複素値 (re, im) で取る(algebraic 系 = 実部、
// magnitude 系 = hypot。tsparse_eigen_selection::target_distance に委譲)。
// 返却集合(locked_values)は実のみ(EIG-3 D3-2: 複素対は converged で返さない)。
//
// 挙動不変の保証: active_imag が全ゼロのとき本関数は real 版
// honest_termination_check_ と同一判定を返す(ks_units 群 [4a] が機械証明)。
// 既存 real 版は無変更であり、EIG-1 完了経路(TRL/lanczos/si_lanczos)は
// 本オーバーロードを呼ばない。
// ---------------------------------------------------------------------------
template <class R>
bool honest_termination_check_complex_(
	const std::vector<R>& active_real,
	const std::vector<R>& active_imag,
	const std::vector<bool>& active_converged,
	const std::vector<R>& locked_values,
	const std::size_t k,
	const eig_target target,
	const R& shift)
{
	typedef std::complex<R> complex_type;

	if (k == 0) return true;
	if (locked_values.size() < k) return false;

	const std::vector<std::size_t> order =
		vcp::tsparse_eigen_selection::select_eigen_indices_from_real(
			locked_values, locked_values.size(), target, shift);
	if (order.size() < k) return false;
	const R worst_returned = locked_values[order[k - 1]];

	const R worst_key = vcp::tsparse_eigen_selection::target_distance(
		complex_type(worst_returned, R(0)), target, shift);
	const bool prefer_large =
		(target == eig_target::largest_magnitude ||
		 target == eig_target::largest_algebraic);

	std::size_t n_active = active_real.size();
	if (active_imag.size() < n_active) n_active = active_imag.size();
	if (active_converged.size() < n_active) n_active = active_converged.size();
	for (std::size_t i = 0; i < n_active; i++) {
		if (active_converged[i]) continue;
		const R key = vcp::tsparse_eigen_selection::target_distance(
			complex_type(active_real[i], active_imag[i]), target, shift);
		const bool certainly_inner = prefer_large ? (key > worst_key)
		                                          : (key < worst_key);
		if (certainly_inner) return false;
	}
	return true;
}

// ---------------------------------------------------------------------------
// honest_termination_check_complex_pairs_   (EIG-4 T-3。純追加 — 既存 2 関数は不変)
//
// 返却集合(locked)に複素共役対を含む場合の C-2 検査。locked の各項目は
// 実固有値(im = 0、1 スロット)または複素対(im > 0 の代表 1 エントリ、
// 2 スロット)。worst key は「target 順に並べた locked をスロット数 k_slots
// まで埋めたときの最後の項目」の target_distance(complex(re, im))。
// active 側の certainly-inner 規則は honest_termination_check_complex_ と同一。
//
// 順序付けは半順序許容の自前選択ループで行う(P6: モジュールスカラー R を
// std::sort 等の SWO 前提アルゴリズムに渡さない)。比較は certainly であり、
// indeterminate は「より優先とは言えない」側に落ちる(選択が保守化するだけで
// 嘘側には倒れない — 最終防衛は C-1 の厳密残差)。
// ---------------------------------------------------------------------------
template <class R>
bool honest_termination_check_complex_pairs_(
	const std::vector<R>& active_real,
	const std::vector<R>& active_imag,
	const std::vector<bool>& active_converged,
	const std::vector<R>& locked_real,
	const std::vector<R>& locked_imag,
	const std::size_t k_slots,
	const eig_target target,
	const R& shift)
{
	typedef std::complex<R> complex_type;

	if (k_slots == 0) return true;
	if (locked_real.size() != locked_imag.size()) return false;

	const bool prefer_large =
		(target == eig_target::largest_magnitude ||
		 target == eig_target::largest_algebraic);

	// locked 項目の key とスロット数
	const std::size_t m = locked_real.size();
	std::vector<R> key(m, R(0));
	std::vector<std::size_t> slots(m, 1);
	std::size_t total_slots = 0;
	for (std::size_t i = 0; i < m; i++) {
		key[i] = vcp::tsparse_eigen_selection::target_distance(
			complex_type(locked_real[i], locked_imag[i]), target, shift);
		if (locked_imag[i] > R(0)) slots[i] = 2;
		total_slots += slots[i];
	}
	if (total_slots < k_slots) return false;   // 返却不足のままの宣言は不正直

	// target 順の自前選択(P6 準拠の線形走査。同格は最小インデックス)
	std::vector<bool> used(m, false);
	std::size_t filled = 0;
	R worst_key = R(0);
	bool have_worst = false;
	while (filled < k_slots) {
		std::size_t best = m;
		for (std::size_t i = 0; i < m; i++) {
			if (used[i]) continue;
			if (best == m) { best = i; continue; }
			const bool better = prefer_large ? (key[i] > key[best])
			                                 : (key[i] < key[best]);
			if (better) best = i;
		}
		if (best == m) break;   // 全消費(total_slots 検査済みなので到達しない)
		used[best] = true;
		filled += slots[best];
		worst_key = key[best];
		have_worst = true;
	}
	if (!have_worst) return false;

	std::size_t n_active = active_real.size();
	if (active_imag.size() < n_active) n_active = active_imag.size();
	if (active_converged.size() < n_active) n_active = active_converged.size();
	for (std::size_t i = 0; i < n_active; i++) {
		if (active_converged[i]) continue;
		const R akey = vcp::tsparse_eigen_selection::target_distance(
			complex_type(active_real[i], active_imag[i]), target, shift);
		const bool certainly_inner = prefer_large ? (akey > worst_key)
		                                          : (akey < worst_key);
		if (certainly_inner) return false;
	}
	return true;
}

// ---------------------------------------------------------------------------
// locked_prefix_indices_   (EIG-1 F-3-1 / F-4 の freshness ガード用ヘルパ)
//
// 返却予定の target 順 prefix k を構成する locked インデックス(昇順)。
// 「現在の部分空間が最後の prefix 変化より後に構築されたか」の判定に使う。
// prefix 外のロックはこの集合を変えないため、外側候補の逐次ロックが
// 検証を永遠に再要求する事態(暴走)を防ぐ。
// ---------------------------------------------------------------------------
template <class T>
std::vector<std::size_t> locked_prefix_indices_(
	const std::vector<T>& locked_values,
	const std::size_t k,
	const eig_target target,
	const typename vcp::tsparse_scalar::real_type<T>::type& shift)
{
	std::vector<std::size_t> order =
		vcp::tsparse_eigen_selection::select_eigen_indices_from_real(
			locked_values, locked_values.size(), target, shift);
	if (order.size() > k) order.resize(k);
	std::sort(order.begin(), order.end());
	return order;
}

// ---------------------------------------------------------------------------
// residual_acceptance_check_   (EIG-1 F-3-2 / F-4 / F-5; EIG-0 C-1 の実装)
//
// 終了時に再評価した「厳密残差」(内部推定値ではない)が全返却対で受理条件
//   res_abs <= tol  または  res_rel <= tol(scale = 各ソルバーの既存定義)
// を満たすときのみ true。converged=true はこの検査を通過した場合に限る。
//
// 成功宣言側は certified-≤(GT1 P1 / B-16): interval では certainly-≤ のみ
// 受理し、indeterminate は失敗側(第三分岐 = 正直な非収束)へ落ちる。
// ---------------------------------------------------------------------------
template <class R>
bool residual_acceptance_check_(
	const std::vector<R>& residuals_absolute,
	const std::vector<R>& residuals_relative,
	const R& tol)
{
	const std::size_t m = residuals_absolute.size();
	if (m == 0) return false;                       // 返却なしに受理なし
	if (residuals_relative.size() != m) return false;  // 診断の不整合は受理しない
	for (std::size_t i = 0; i < m; i++) {
		const bool accepted = (residuals_absolute[i] <= tol) ||
		                      (residuals_relative[i] <= tol);
		if (!accepted) return false;
	}
	return true;
}

// ---------------------------------------------------------------------------
// EIG-4 T-4 (D4-4、R-1 改訂式): 改訂 C-1 受理スケールの唯一の定義(B-29:
// scale 定義の分散・重複を禁止。ソルバー側での再定義は不可)。
//
//   scale(|θ|, ‖A‖_est) = max(1 + |θ|, ‖A‖_est)
//
// 受理は「既存受理(res_abs ≤ tol ∨ res_rel ≤ tol)∨ res_abs ≤ tol·scale」。
// 既存受理式を 1 ビットも変えない OR 追加なので厳密に広義緩和であり、
// 許される v1 遷移は「正直 NONCONV → 正解 OK」のみ(遷移許容集合 = {t1b})。
//
// ‖A‖_est の経路別定義(この 1 箇所に集約・文書化):
//   - 行列を所持する dispatch 層の最終ゲート(lanczos_package_to_result_ /
//     lambda_c1_acceptance_gate_): ‖A‖∞ の厳密値(O(nnz)、決定的)
//   - TRL(apply 抽象、行列非所持): 実行全体で観測した射影 3 重対角の
//     ∞ ノルム走行最大 max_i(|β_{i-1}| + |α_i| + |β_i|)。|α_i|,|β_i| ≤ ‖A‖₂
//     より高々 3‖A‖₂(対称では ≤ 3‖A‖∞)の決定的推定であり、‖A‖∞ の厳密値
//     とは一致しない(経路間で scale 値は同一でない — G-1.1 R-1 付記)。
//     追加 apply ゼロ・実行軌跡不変。
//   - KS ドライバの実行中受理は従来 scale(1 + |θ|)のまま(契約より厳しい
//     側は C-1 の含意を破らない。実行中受理の緩和は実行軌跡と mv を変え
//     v1 に宣言外差分を作るため行わない — 遷移許容集合文書 §3)
//   - dense 経路は従来の受理スケール(max(tol, tol·10n))を維持(既存定義の
//     踏襲。EIG-0 C-1 の「各ソルバーの既存定義を文書化」条項)
// ---------------------------------------------------------------------------
template <class R>
R c1_revised_scale_(const R& theta_abs, const R& anorm_est)
{
	R s = R(1) + theta_abs;
	if (anorm_est > s) s = anorm_est;
	return s;
}

// 最終 verdict ゲート用の改訂受理(成功宣言側 = certified-≤、GT1 P1/B-16。
// interval では indeterminate が失敗側 = 正直な非収束に落ちる)。
template <class R>
bool residual_acceptance_check_scaled_(
	const std::vector<R>& residuals_absolute,
	const std::vector<R>& residuals_relative,
	const R& tol,
	const std::vector<R>& theta_abs,
	const R& anorm_est)
{
	const std::size_t m = residuals_absolute.size();
	if (m == 0) return false;
	if (residuals_relative.size() != m) return false;
	if (theta_abs.size() != m) return false;
	for (std::size_t i = 0; i < m; i++) {
		const bool legacy = (residuals_absolute[i] <= tol) ||
		                    (residuals_relative[i] <= tol);
		const bool revised = (residuals_absolute[i]
		    <= tol * c1_revised_scale_<R>(theta_abs[i], anorm_est));
		if (!(legacy || revised)) return false;
	}
	return true;
}

} // namespace tsparse
} // namespace vcp

#endif // VCP_TSPARSE_HONEST_TERMINATION_HPP
