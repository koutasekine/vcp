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

} // namespace tsparse
} // namespace vcp

#endif // VCP_TSPARSE_HONEST_TERMINATION_HPP
