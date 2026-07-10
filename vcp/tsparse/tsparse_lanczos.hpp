// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_LANCZOS_HPP
#define VCP_TSPARSE_LANCZOS_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_eigensolvers.hpp>
#include <vcp/tsparse/tsparse_honest_termination.hpp>
#include <vcp/tsparse/tsparse_projected_eigensolver.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_lanczos {

// ---------------------------------------------------------------------------
// Generate a deterministic starting vector (based on indices)
// ---------------------------------------------------------------------------
template <typename T>
std::vector<T> deterministic_start_vector(const std::size_t n, const unsigned int seed) {
	typedef typename tsparse_scalar::real_type<T>::type R;
	std::vector<T> v(n, T(0));
	// LCG-based pseudorandom (simple, deterministic, no <random> dependency issue)
	unsigned int state = seed + 1u;
	for (std::size_t i = 0; i < n; i++) {
		state = state * 1664525u + 1013904223u;
		// Map to [-1, 1]
		const R val = R(static_cast<int>(state >> 16)) / R(32768) - R(1);
		v[i] = T(val);
	}
	const R nv = tsparse_scalar::real_norm_value(v);
	if (nv > R(0)) {
		for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	}
	return v;
}

// ---------------------------------------------------------------------------
// Orthogonalize v against a list of vectors (Gram-Schmidt, full)
// Returns false if resulting vector is (near) zero.
// ---------------------------------------------------------------------------
template <typename T>
bool orthogonalize_against(std::vector<T>& v,
                            const std::vector<std::vector<T> >& vecs,
                            const typename tsparse_scalar::real_type<T>::type& tol)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	const std::size_t n = v.size();
	// Two passes for numerical stability
	for (int pass = 0; pass < 2; pass++) {
		for (std::size_t j = 0; j < vecs.size(); j++) {
			const R c = tsparse_scalar::real_dot_value(vecs[j], v);
			for (std::size_t i = 0; i < n; i++) v[i] -= T(c) * vecs[j][i];
		}
	}
	const R nv = tsparse_scalar::real_norm_value(v);
	if (nv <= tol) return false;
	for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
	return true;
}

// ---------------------------------------------------------------------------
// Symmetric tridiagonal eigensolver (Jacobi-based, small m)
// Returns eigenvalues and eigenvectors of tridiagonal T(alpha, beta)
// Phase 8.1: now routes through solve_real_symmetric_projected for
// validated, dimension-safe projected solve.
// ---------------------------------------------------------------------------
template <typename T>
tsparse_dense_linalg::dense_eigen_result<T> solve_tridiag(
	const std::vector<T>& alpha,
	const std::vector<T>& beta,
	const std::size_t max_iter,
	const typename tsparse_scalar::real_type<T>::type& tol)
{
	const std::size_t m = alpha.size();
	std::vector<std::vector<T> > Tm(m, std::vector<T>(m, T(0)));
	for (std::size_t i = 0; i < m; i++) {
		Tm[i][i] = alpha[i];
		if (i + 1 < m) { Tm[i][i+1] = beta[i]; Tm[i+1][i] = beta[i]; }
	}
	vcp::tsparse_projected::projected_eigensolver_result<T> proj =
		vcp::tsparse_projected::solve_real_symmetric_projected(Tm, max_iter, tol);
	tsparse_dense_linalg::dense_eigen_result<T> result;
	result.converged = proj.success && proj.converged;
	result.iterations = proj.iterations;
	if (!proj.success) {
		return result;
	}
	result.eigenvalues.reserve(proj.pairs.size());
	result.eigenvectors.reserve(proj.pairs.size());
	for (std::size_t i = 0; i < proj.pairs.size(); i++) {
		result.eigenvalues.push_back(proj.pairs[i].value);
		result.eigenvectors.push_back(proj.pairs[i].vector);
	}
	return result;
}

// ---------------------------------------------------------------------------
// EIG-6 F-2' (G-2.1 承認): si_lanczos 経路の opt-in リスタート制御。
//
// 既定(si_spin_retention = false)は従来挙動と完全同一(公開 lanczos
// 直接経路は無変更 — v1 凍結スイートで機械検証)。opt-in 時のみ:
//  (b') スピン誘発 warm start — 「ロックなし(かつ locked<k)で終わった
//       リスタート」の直後のみ、次リスタートの開始ベクトルを直前リスタートの
//       target 側先頭候補の Ritz ベクトル(locked 直交化済み)にする。
//       ロック成立で解除(次は従来の新規決定的ベクトル = 1 リスタート
//       1 ロックの多重度発見機構は無変更)。決定的・乱数なし。
//  (e)  λ 形式 lock 併記 — μ 空間 lock 基準(現行・無変更)に加え、
//       呼び出し側が渡す λ 空間ゲート(契約 C-1 の改訂 scale 受理式そのもの。
//       共有ヘルパ再利用は呼び出し側 = spmats_eigs.hpp 側の責務)での
//       lock を許す。適用は走査先頭 min(k_remaining+2, sdim/2) 候補のみ。
//       判定用 A·x はゲート内で mv に 1:1 計上し、lambda_gate_products にも
//       別建て計上する(B-38 / (e-2))。
//  (d') 走査キャップ — lock 適格性走査を min(k_remaining+2, sdim/2) で
//       打ち切る。C-2/honest termination の候補全景(cand_vals)は sel 全体
//       から構築され、キャップの影響を受けない((d') 承認条件)。
//       リスタートあたり最悪 mv = sdim + 2·min(k+2, sdim/2) ≤ 2·sdim で
//       B-1 の機械照合上限 mv ≤ max_iter + 2·sdim は不変。
//
// 背景(D-17b 機構、EIG-6_G2.1_submission.md):
//  b-1: plain restart の部分空間破棄 → 単発 sdim step の残差床が tol に
//       届かずロック飢餓(grid_large: 47/47 リスタート locked=0)。
//  b-2: 高条件数では再適用 μ 残差に LU forward-error 床(~eps·κ(A)·μ)が乗り、
//       μ 空間絶対 lock 基準が構造的に到達不能(bcsstk18 模擬: 床 4e-6)。
//       λ 空間残差の床は backward stability により ~eps·‖A‖ で契約 C-1 は
//       満たせる — 契約より厳しすぎる内部基準のみを λ 形式で併記救済する。
// ---------------------------------------------------------------------------
struct si_lambda_lock_gate_none {
	template <typename VEC, typename MU>
	bool operator()(const VEC&, const MU&, std::size_t&) const { return false; }
};

// ---------------------------------------------------------------------------
// Full-reorthogonalized standard Lanczos eigensolver
//
// ApplyA: void(const vector<T>& x, vector<T>& y)  →  y = A x
// Returns eig_result<T> with all required fields populated.
// ---------------------------------------------------------------------------
template <typename T, class ApplyA>
struct lanczos_result_package {
	std::vector<T> eigenvalues;
	std::vector<std::vector<T> > eigenvectors;
		std::vector<typename tsparse_scalar::real_type<T>::type> residuals_abs;
		std::vector<typename tsparse_scalar::real_type<T>::type> residuals_rel;
		std::vector<typename tsparse_scalar::real_type<T>::type> history_abs;
		std::vector<typename tsparse_scalar::real_type<T>::type> history_rel;
	bool converged;
	std::size_t converged_count;
	std::size_t returned_count;
	std::size_t iterations;
	std::size_t restarts;
	std::size_t mv_count;
	std::string breakdown_reason;
	std::string failure_reason;
	std::string used_method;
};

template <typename T, class ApplyA, class ApplyNorm,
          class LambdaLockGate = si_lambda_lock_gate_none>
lanczos_result_package<T, ApplyA> lanczos_eigs(
	const std::size_t n,           // matrix size
	const std::size_t k,           // number requested
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const unsigned int random_seed,
	const bool random_start,
	const eig_target target,
	const typename tsparse_scalar::real_type<T>::type& shift_val,
	const bool compute_residual_history,
	ApplyA apply,                  // y = A x
	ApplyNorm apply_norm,          // ||v||_B (for standard: standard norm)
	// EIG-6 F-2'(opt-in、既定 = 従来挙動。ヘッダ先頭の設計コメント参照)
	const bool si_spin_retention = false,
	LambdaLockGate lambda_lock_gate = LambdaLockGate(),
	const bool use_lambda_lock_gate = false,
	std::size_t* lambda_gate_products = 0)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	lanczos_result_package<T, ApplyA> res;
	res.converged = false;
	res.converged_count = 0;
	res.returned_count = 0;
	res.iterations = 0;
	res.restarts = 0;
	res.mv_count = 0;
	res.used_method = "lanczos";

	if (n == 0 || k == 0) return res;

	const std::size_t m_limit = std::min(n, std::max(subspace_dim, k + std::size_t(3)));
	const R small_tol = tol / R(10);

	// locked eigenvalues and their eigenvectors
	std::vector<T> locked_vals;
	std::vector<std::vector<T> > locked_vecs;
	std::vector<R> locked_res;
	std::vector<T> best_vals;
	std::vector<std::vector<T> > best_vecs;
	std::vector<R> best_res;

	unsigned int seed_counter = random_start ? random_seed : 0u;

	// EIG-6 F-2' (b'): warm start 状態(si_spin_retention 時のみ使用)
	std::vector<T> warm_vec;
	bool warm_valid = false;

	for (std::size_t restart = 0; restart <= max_restarts; restart++) {
		// EIG-1 F-4: 終了は本ループ末尾の honest termination ゲート
		// (C-2 + freshness)だけが宣言する。locked >= k でも検証未了なら
		// 継続する(max_restarts 到達で正直な非収束)。
		res.restarts = restart;
		const bool warm_this_restart = si_spin_retention && warm_valid;

		// freshness ガードの基準: このリスタートの部分空間は「今の locked
		// 集合」に直交して構築される。その時点の返却予定 prefix を記録し、
		// ゲートでは「このリスタート中に prefix が変化していない」ことを
		// 要求する(prefix 外のロックは対象外 — 外側候補の逐次ロックが
		// 検証を永遠に再要求する暴走を防ぐ)。
		const std::vector<std::size_t> prefix_at_restart_start =
			vcp::tsparse::locked_prefix_indices_(locked_vals, k, target, shift_val);

		// Build starting vector orthogonal to locked vectors.
		// EIG-6 F-2' (b'): warm リスタートでは直前リスタートの target 側先頭
		// 候補の Ritz ベクトルから開始(直交化退化時は従来の決定的系列へ
		// フォールバック)。それ以外は従来どおり新規決定的ベクトル。
		std::vector<T> v0;
		bool v0_ready = false;
		if (warm_this_restart) {
			v0 = warm_vec;
			v0_ready = orthogonalize_against(v0, locked_vecs, small_tol);
		}
		if (!v0_ready) {
			v0 = deterministic_start_vector<T>(n, seed_counter);
			seed_counter++;
			if (!orthogonalize_against(v0, locked_vecs, small_tol)) {
				// Try a few more seeds
				bool ok = false;
				for (int attempt = 0; attempt < 20 && !ok; attempt++) {
					v0 = deterministic_start_vector<T>(n, seed_counter++);
					ok = orthogonalize_against(v0, locked_vecs, small_tol);
				}
				if (!ok) { res.failure_reason = "exhausted orthogonal starting vectors"; break; }
			}
		}

		// Run Lanczos from v0
		const std::size_t m_active = std::min(
			(n > locked_vals.size()) ? (n - locked_vals.size()) : std::size_t(0),
			m_limit);
		if (m_active == 0) {
			// locked が全空間を張った: 部分空間に候補は存在し得ないので
			// C-2 は自明に成立する(空の active 集合)。k 個そろっていれば
			// 正直に converged としてよい(C-1 はロック時の厳密残差 +
			// 呼び出し側の λ 空間受理検査が守る)。
			if (locked_vals.size() >= k) res.converged = true;
			break;
		}

		std::vector<std::vector<T> > basis;
		basis.reserve(m_active + 1);
		basis.push_back(v0);

		std::vector<T> alpha_vec, beta_vec;
		std::vector<T> v_prev(n, T(0));
		R beta_prev = R(0);
		std::string bd_reason;

		for (std::size_t j = 0; j < m_active; j++) {
			// z = A * basis[j]
			std::vector<T> z;
			apply(basis[j], z);
			res.mv_count++;

			if (j > 0) {
				for (std::size_t i = 0; i < n; i++) z[i] -= T(beta_prev) * v_prev[i];
			}

			// Orthogonalize z against locked vectors
			for (std::size_t li = 0; li < locked_vecs.size(); li++) {
				const R c = tsparse_scalar::real_dot_value(locked_vecs[li], z);
				for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * locked_vecs[li][i];
			}

			const R a = tsparse_scalar::real_dot_value(basis[j], z);
			alpha_vec.push_back(T(a));
			for (std::size_t i = 0; i < n; i++) z[i] -= T(a) * basis[j][i];

			for (int pass = 0; pass < 2; pass++) {
				for (std::size_t bi = 0; bi < basis.size(); bi++) {
					const R c = tsparse_scalar::real_dot_value(basis[bi], z);
					for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * basis[bi][i];
				}
				for (std::size_t li = 0; li < locked_vecs.size(); li++) {
					const R c = tsparse_scalar::real_dot_value(locked_vecs[li], z);
					for (std::size_t i = 0; i < n; i++) z[i] -= T(c) * locked_vecs[li][i];
				}
			}

			const R bn = apply_norm(z);
			if (compute_residual_history) {
				res.history_abs.push_back(bn);
				res.history_rel.push_back(bn / (R(1) + tsparse_scalar::abs_value(a)));
			}

			if (bn <= small_tol) {
				bd_reason = "happy breakdown";
				break;
			}

			beta_vec.push_back(T(bn));
			v_prev = basis[j];
			beta_prev = bn;
			std::vector<T> vnew(n);
			for (std::size_t i = 0; i < n; i++) vnew[i] = z[i] / T(bn);
			basis.push_back(vnew);
		}

		res.iterations += alpha_vec.size();

		if (alpha_vec.empty()) { res.breakdown_reason = "empty Lanczos basis"; break; }

		// Solve small tridiagonal eigenproblem
		const std::size_t m = alpha_vec.size();
		const std::size_t proj_iter = std::max(m * m * 100, std::size_t(1000));
		tsparse_dense_linalg::dense_eigen_result<T> small =
			solve_tridiag(alpha_vec, beta_vec, proj_iter, small_tol);

		if (!small.converged) {
			res.breakdown_reason = "projected_eigensolver_failed";
			res.failure_reason = "projected eigensolver failed or did not converge";
			break;
		}
		if (small.eigenvalues.empty()) {
			res.breakdown_reason = "projected_eigensolver_empty_result";
			res.failure_reason = "projected eigensolver returned empty result";
			break;
		}

		// Select target eigenvalues
		std::vector<std::complex<R> > ceigs;
		ceigs.reserve(small.eigenvalues.size());
		for (std::size_t i = 0; i < small.eigenvalues.size(); i++) {
			ceigs.push_back(std::complex<R>(tsparse_scalar::real_part(small.eigenvalues[i]), R(0)));
		}
		// EIG-1 F-4: locked >= k の検証リスタートでも 1 ロック可能にする
		// (size_t underflow 対策込み)。1 リスタート 1 ロックは多重度発見
		// 機構として維持する。
		const std::size_t k_remaining = (locked_vals.size() < k)
			? (k - locked_vals.size()) : std::size_t(1);
		std::vector<std::size_t> sel = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs, std::min(m, k_remaining + m), target, shift_val);

		// Check convergence for each selected Ritz pair.
		// Lock at most ONE new pair per restart so that repeated eigenvalues
		// are found on subsequent restarts (each with a new deflated start vector).
		// EIG-6 F-2' (d'): opt-in 時は lock 適格性走査を
		// min(k_remaining+2, sdim/2) で打ち切る(非ロック時の全数残差評価の
		// 浪費を抑制。C-2 の候補全景 cand_vals は sel 全体から構築するため
		// 影響しない)。λ 形式 lock (e) の適用窓も同じキャップ。
		const std::size_t si_scan_cap = si_spin_retention
			? std::min(k_remaining + std::size_t(2),
			           std::max(std::size_t(1), m_limit / 2))
			: sel.size();
		std::size_t si_scanned = 0;
		bool locked_this_restart = false;
		std::vector<T> si_first_cand_vec;
		bool si_first_cand_set = false;
		for (std::size_t si = 0; si < sel.size(); si++) {
			if (si_spin_retention && si_scanned >= si_scan_cap) break;
			const std::size_t idx = sel[si];
			if (idx >= small.eigenvectors.size()) continue;

			// Lift Ritz vector
			const std::vector<T>& y = small.eigenvectors[idx];
			std::vector<T> ritz_vec(n, T(0));
			for (std::size_t j = 0; j < y.size() && j < basis.size(); j++) {
				for (std::size_t i = 0; i < n; i++) ritz_vec[i] += basis[j][i] * y[j];
			}
			// Normalize and orthogonalize against already-locked vectors
			for (std::size_t li = 0; li < locked_vecs.size(); li++) {
				const R c = tsparse_scalar::real_dot_value(locked_vecs[li], ritz_vec);
				for (std::size_t i = 0; i < n; i++) ritz_vec[i] -= T(c) * locked_vecs[li][i];
			}
			const R ritz_nrm = apply_norm(ritz_vec);
			if (ritz_nrm <= small_tol) continue;
			for (std::size_t i = 0; i < n; i++) ritz_vec[i] /= T(ritz_nrm);

			// Compute residual ||A v - mu v||
			std::vector<T> Av;
			apply(ritz_vec, Av);
			res.mv_count++;
			const T mu = small.eigenvalues[idx];
			std::vector<T> r(n);
			for (std::size_t i = 0; i < n; i++) r[i] = Av[i] - mu * ritz_vec[i];
			const R res_abs = tsparse_scalar::real_norm_value(r);
			const R res_rel = res_abs / (R(1) + tsparse_scalar::abs_value(mu));

				if (compute_residual_history) {
					res.history_abs.push_back(res_abs);
					res.history_rel.push_back(res_rel);
				}
			si_scanned++;
			if (si_spin_retention && !si_first_cand_set) {
				si_first_cand_vec = ritz_vec;   // target 側先頭候補(次回 warm 用)
				si_first_cand_set = true;
			}
			if (best_vals.size() < k) {
				best_vals.push_back(mu);
				best_vecs.push_back(ritz_vec);
				best_res.push_back(res_abs);
			}

			// EIG-6 F-2' (e): μ 空間 lock 基準(現行・無変更)∨ λ 形式ゲート
			// (μ 不合格の窓内候補のみ消費。契約 C-1 の改訂 scale 受理式は
			// 呼び出し側 = spmats_eigs.hpp が共有ヘルパで実装)。
			// λ ゲートの A·x はゲート内で mv に 1:1 計上 + 別建てカウント。
			// 注記: warm ロックへの λ 認証必須化は検討・実測の結果、単発 warm
			// 反復の決定的固定点(r_μ 床)により不能と確定 — 品質の最終確定は
			// back-half の θ シフト磨き(spmats_eigs.hpp)が担う。
			bool lock_ok = (res_abs <= tol || res_rel <= tol);
			if (!lock_ok && use_lambda_lock_gate && si_scanned <= si_scan_cap) {
				const std::size_t mv_before = res.mv_count;
				lock_ok = lambda_lock_gate(ritz_vec, mu, res.mv_count);
				if (lambda_gate_products != 0)
					*lambda_gate_products += (res.mv_count - mv_before);
			}
			if (lock_ok) {
				locked_vals.push_back(mu);
				locked_vecs.push_back(ritz_vec);
				locked_res.push_back(res_abs);
				res.converged_count++;
				locked_this_restart = true;
				break;  // one lock per restart → next restart deflates this vector out
			}
		}

		// EIG-6 F-2' (b'): スピン誘発 warm start の発火/解除。
		// 発火 = ロックなし(かつ locked<k)で終わった直後のみ。
		// 解除 = ロック成立(次リスタートは従来の新規決定的ベクトル —
		// 多重度発見の fresh 復帰)。
		if (si_spin_retention) {
			if (!locked_this_restart && locked_vals.size() < k && si_first_cand_set) {
				warm_vec = si_first_cand_vec;
				warm_valid = true;
			} else {
				warm_valid = false;
			}
		}

		if (!bd_reason.empty() && locked_vals.size() < k) {
			// Happy breakdown with insufficient locked pairs → need restart
			// Keep going to next restart
		}

		if (locked_vals.size() >= k) {
			// ---------------------------------------------------------------
			// EIG-1 F-4 (EIG-0 C-2 + freshness): converged は次の両方が成立
			// する場合にのみ宣言する。
			//  (i)  freshness — このリスタート中に返却予定 prefix k が変化して
			//       いない。公開 lanczos は毎リスタートが locked 直交の新規
			//       ランダム開始なので、prefix 不変のリスタートの候補 =
			//       現在の返却集合を知る新鮮な部分空間からの証拠である。
			//       prefix を変えるロック(新規参入)が起きた直後の候補は
			//       多重度コピーの盲点を持ち得る(TRL の t1 機構と同じ)ため
			//       次の検証リスタートを待つ。prefix 外のロックは検証を
			//       再要求しない(外側候補の逐次ロックによる暴走防止)。
			//  (ii) 現部分空間の候補に、返却予定 prefix k の最悪値より内側の
			//       未収束候補が certainly 存在しない(共通コア検査)。
			// 残差未評価の候補(ロック break 後の未走査分)も未収束扱い
			// (保守側)。不成立なら検証リスタートを継続し、max_restarts
			// 到達で正直な非収束になる。
			// ---------------------------------------------------------------
			const std::vector<std::size_t> prefix_now =
				vcp::tsparse::locked_prefix_indices_(locked_vals, k, target, shift_val);
			const bool prefix_fresh = (prefix_now == prefix_at_restart_start);
			if (prefix_fresh) {
				std::vector<T> cand_vals;
				std::vector<bool> cand_conv;
				cand_vals.reserve(sel.size());
				cand_conv.reserve(sel.size());
				for (std::size_t si = 0; si < sel.size(); si++) {
					const std::size_t idx = sel[si];
					if (idx >= small.eigenvalues.size()) continue;
					cand_vals.push_back(small.eigenvalues[idx]);
					cand_conv.push_back(false);
				}
				if (vcp::tsparse::honest_termination_check_(
						cand_vals, cand_conv, locked_vals, k, target, shift_val)) {
					res.converged = true;
					res.breakdown_reason = bd_reason;
					break;
				}
			}
		}
	}

	// Sort locked eigenvalues by target
	if (!locked_vals.empty()) {
		std::vector<std::complex<R> > ceigs_locked;
		ceigs_locked.reserve(locked_vals.size());
		for (std::size_t i = 0; i < locked_vals.size(); i++) {
			ceigs_locked.push_back(std::complex<R>(tsparse_scalar::real_part(locked_vals[i]), R(0)));
		}
		const std::size_t ksorted = std::min(k, locked_vals.size());
		std::vector<std::size_t> order = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs_locked, ksorted, target, shift_val);

		std::vector<T> sorted_vals;
		std::vector<std::vector<T> > sorted_vecs;
		std::vector<R> sorted_res;
		for (std::size_t i = 0; i < order.size(); i++) {
			sorted_vals.push_back(locked_vals[order[i]]);
			sorted_vecs.push_back(locked_vecs[order[i]]);
			sorted_res.push_back(locked_res[order[i]]);
		}
		res.eigenvalues = sorted_vals;
		res.eigenvectors = sorted_vecs;
		res.residuals_abs = sorted_res;
		res.residuals_rel.resize(sorted_res.size());
		for (std::size_t i = 0; i < sorted_vals.size(); i++) {
			res.residuals_rel[i] = sorted_res[i] / (R(1) + tsparse_scalar::abs_value(sorted_vals[i]));
		}
		res.returned_count = sorted_vals.size();
	}

	if (!res.converged && res.failure_reason.empty()) {
		res.failure_reason = "maximum restarts reached without full convergence";
	}
	if (res.eigenvalues.empty() && !best_vals.empty()) {
		std::vector<std::complex<R> > ceigs_best;
		ceigs_best.reserve(best_vals.size());
		for (std::size_t i = 0; i < best_vals.size(); i++) {
			ceigs_best.push_back(std::complex<R>(tsparse_scalar::real_part(best_vals[i]), R(0)));
		}
		const std::vector<std::size_t> order = tsparse_eigensolvers::select_ritz_indices<T>(
			ceigs_best, std::min(k, best_vals.size()), target, shift_val);
		for (std::size_t i = 0; i < order.size(); i++) {
			res.eigenvalues.push_back(best_vals[order[i]]);
			res.eigenvectors.push_back(best_vecs[order[i]]);
			res.residuals_abs.push_back(best_res[order[i]]);
		}
		res.residuals_rel.resize(res.residuals_abs.size());
		for (std::size_t i = 0; i < res.eigenvalues.size(); i++) {
			res.residuals_rel[i] = res.residuals_abs[i] / (R(1) + tsparse_scalar::abs_value(res.eigenvalues[i]));
		}
		res.returned_count = res.eigenvalues.size();
	}

	return res;
}

// Convenience overload with standard L2 norm
template <typename T, class ApplyA>
lanczos_result_package<T, ApplyA> lanczos_eigs_standard(
	const std::size_t n,
	const std::size_t k,
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const unsigned int random_seed,
	const bool random_start,
	const eig_target target,
	const typename tsparse_scalar::real_type<T>::type& shift_val,
	const bool compute_residual_history,
	ApplyA apply)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	struct norm_fn {
		R operator()(const std::vector<T>& v) const {
			return tsparse_scalar::real_norm_value(v);
		}
	} norm_func;
	return lanczos_eigs<T, ApplyA, norm_fn>(
		n, k, subspace_dim, max_restarts, tol,
		random_seed, random_start, target, shift_val,
		compute_residual_history, apply, norm_func);
}

// EIG-6 F-2': si_lanczos back-half 専用の opt-in 版(標準 L2 ノルム)。
// 公開 lanczos の既存 wrapper(上)は無変更。
template <typename T, class ApplyA, class LambdaLockGate>
lanczos_result_package<T, ApplyA> lanczos_eigs_standard_si_(
	const std::size_t n,
	const std::size_t k,
	const std::size_t subspace_dim,
	const std::size_t max_restarts,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const unsigned int random_seed,
	const bool random_start,
	const eig_target target,
	const typename tsparse_scalar::real_type<T>::type& shift_val,
	const bool compute_residual_history,
	ApplyA apply,
	LambdaLockGate lambda_lock_gate,
	std::size_t* lambda_gate_products)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	struct norm_fn {
		R operator()(const std::vector<T>& v) const {
			return tsparse_scalar::real_norm_value(v);
		}
	} norm_func;
	return lanczos_eigs<T, ApplyA, norm_fn, LambdaLockGate>(
		n, k, subspace_dim, max_restarts, tol,
		random_seed, random_start, target, shift_val,
		compute_residual_history, apply, norm_func,
		/*si_spin_retention=*/true, lambda_lock_gate,
		/*use_lambda_lock_gate=*/true, lambda_gate_products);
}

} // namespace tsparse_lanczos
} // namespace vcp

#endif
