// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spimats.hpp
// SPI-1: 疎行列の区間型(精度保証)ポリシー spimats の骨組み。
// 設計書: sandbox/docs/design/spimats_design_v0.md(v0、オーナー承認済み)
// 実装指示書: sandbox/docs/plans/SPI-1_implementation_directive.md
//
// 位置づけ(設計 G1): spimats<T, AP, DP> : public spmats<kv::interval<T>>。
//   vcp::spmatrix<kv::interval<T>, vcp::spimats<T, AP, DP>> として既存 API
//   (solve / eig / eigs / generalized_eig(s))経由で保証付き結果を得る。
//   override は NVI virtual _impl 3 本のみ(§3)。
//
// 意味論の芯(設計 G3/G4、二段 throw = §4):
//   - 成功時のみ result が返る。保証の失敗は「十分条件が満たせなかった」
//     だけで対象について何の主張も生まない(無主張)。部分結果・NaN 埋め・
//     近似値の区間化返却は行わない。
//   - 入口検査の misuse は従来どおり vcp::error 系 throw。
//     AP(近似層)の失敗は vcp::state_error。保証の失敗は
//     vcp::verification_error。
//   - 【既存規約からの意図的逸脱】spmats / spumar / spcmodumar の
//     「runtime failure は status、throw しない」規約に対し、本ポリシーは
//     _with_info 系でも throw する。既存規約の対象は近似計算の runtime
//     failure(失敗も正常な結果の一つ)だが、保証の失敗は無主張であり
//     返す値が存在しないため(設計 §4.1 で確定)。
//   - 点の近似診断(converged / residuals / iterations)は区間の保証の器に
//     流用しない(G4)。近似層のスナップショットは last_verify_info(点型)
//     にのみ置く。result のフィールド有効性は設計 §5.3 の表に従う。
//
// 丸め規律(G5): 本ヘッダおよび spimats_base/ は FP 環境の丸めモード切替
//   (C の丸め切替関数・VCP の丸めガード)や自前の方向丸め算術を一切使用
//   しない。丸めは kv(mid / rad / midrad / 区間演算)と vcp::convert に
//   委譲する。
//
// 区間比較規約(設計 §8): 区間値の同一性判定に kv の == / != を使わない
//   (== は certainly-equal [確認済み kv/interval.hpp L805]、!= は
//   no-overlap [L828] であり端点一致・その否定ではない)。同一性は
//   .lower()/.upper() の端点スカラー比較、包含は subset、重なりは overlap。
//   例外:「厳密に点 0 か」の判定に限り x == _T(0) は正しい。
//
// 区画体系(設計 §11):
//   USER-IMPLEMENTS [SPI-Uk]: 精度保証の実装区画。正しさはオーナーが実装・
//     保証。Claude の想定コードコメントは補助であり拘束しない。
//   REVIEW-REQUIRED [SPI-Rk]: Claude 実装だが保証の連鎖に乗るため
//     オーナーレビュー承認(STATUS 書換はオーナーのみ)が編入条件。
//   リリース条件: 全 R 区画の STATUS がオーナー承認で更新されている
//   (未レビュー状態の区画ヘッダが grep で 1 件も出ない)こと。
//
// アクセス経路(設計 §2.3): spmatrix<_T,_P> は protected _P 継承のため、
//   spimats 固有の公開メンバ(is_spd / spi_opt / last_verify_info / 混合
//   演算)は spmatrix ラッパー越しには呼べない。solve / eig 系は spmatrix の
//   既存 API から override に到達する。固有メンバはポリシーオブジェクトの
//   直接利用(vcp::spimats<double> A; A.add(...); A.finalize(); A.is_spd();)
//   で使う(spumar の muar_opt と同方式)。詳細は
//   sandbox/docs/usage/spimats_usage.md。

#ifndef VCP_SPIMATS_HPP
#define VCP_SPIMATS_HPP

#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <kv/interval.hpp>

#include <vcp/error.hpp>
#include <vcp/imats.hpp>       // 区間密ポリシー(_DP 既定 = imats<_T>。mats.hpp 同梱)
#include <vcp/spmats.hpp>
#include <vcp/spimats_base/spimats_convert.hpp>
#include <vcp/spimats_base/spimats_mixed_ops.hpp>

namespace vcp {

namespace spimats_detail {

	// policy_is_symmetric(tol 付き)の static_assert 封鎖用(依存 false)。
	template <typename _Any>
	struct always_false_ : std::false_type {};

} // namespace spimats_detail

// ---------------------------------------------------------------------------
// spimats_options: 保証層のノブ(spumar の muar_opt 方式)。近似層のノブ
// (linear_solve_options / eig_options)とは混在させない(設計 §5.1)。
// ---------------------------------------------------------------------------
struct spimats_options {
	// 近似固有ベクトルのスナップショット opt-in(設計 §5.2: 既定では保持
	// しない --- 全固有値時に密 n×n が乗るため。throw 後の検死用)。
	bool keep_ap_eigenvectors;
	// [SPI-F4 F4-1] veigs: smallest_magnitude / largest_magnitude / target_* の
	// 要求順位(最も近い / 最も大きい)を慣性で検査するか。
	//   off     = 検査しない(last_verify_info.veigs_rank_certified = 0 = 未検査)
	//   record  = 検査して記録するだけ(既定。返却・無主張は SPI-F3 と同一)
	//   require = 証明できなければ verification_error(無主張)
	// sa / la と Rayleigh 退路は F3-4 で順位証明済みのためノブに関わらず 1 を記録。
	enum class target_rank_check_mode { off, record, require };
	target_rank_check_mode target_rank_check;

	spimats_options()
		: keep_ap_eigenvectors(false),
		  target_rank_check(target_rank_check_mode::record) {}
};

// ---------------------------------------------------------------------------
// spd_options<_T>: is_spd のノブ。INTLAB isspd.m の mmd / nosymcheck 相当
// (設計 §7.2。accurate は SPI-G1 D-4 で削除 --- 証明経路は 1 本)。
// ---------------------------------------------------------------------------
template <typename _T>
struct spd_options {
	// 区間版 is_spd の厳密対称入口検査(SPI-R4)を省略する(isspd.m の
	// nosymcheck 相当)。省略時の対称性の責任は呼び出し側。
	bool nosymcheck;
	// isspd.m の mmd 相当: 既存 chol ordering 系へ接続(sparse_chol_ordering
	// [tsparse_sparse_chol.hpp])。順列は is_spd 内で行列に先に適用し(値移動
	// のみ)、Cholesky は常に基底 spmats の自前 chol(演算順が管理下)を
	// natural で実行する。外部委譲(CHOLMOD 等)は主張の正当化を文書化して
	// から有効化するオプションの口のみ。
	vcp::sparse_chol_ordering ordering;
	// 'accurate' removed (SPI-G1 D-4): there is exactly one certification route
	// (Rump 1993 Thm 2.1 / Alg. 2.2 type a-posteriori evaluation, T-generic).
	// Keeping the field and ignoring it would be a silent fallback of an
	// explicitly requested mode. The double-only fast route (Rump 2006) returns
	// in the derived policy spidmats (separate track).

	spd_options()
		: nosymcheck(false),
		  ordering(vcp::sparse_chol_ordering::auto_select) {}
};

// ---------------------------------------------------------------------------
// spimats_verify_info<_T>: フライトレコーダ(設計 §5.4)。
// これは【計算の記録】であって対象行列についての主張ではない(G3 と矛盾
// しない)。mutable メンバとして検証の進行に沿って逐次書き込まれ、throw 後に
// ポリシーオブジェクトから検死できる。
// 近似値は点型 std::vector<_T> のまま置く: 区間型 std::vector<interval<_T>>
// を期待する保証側コードへ型レベルで流入できない(設計 §5.2)。
// ---------------------------------------------------------------------------
template <typename _T>
struct spimats_verify_info {
	// --- 近似層スナップショット(点型のまま) ---
	linear_solve_result<_T> ap_lss_result;          // lss: AP 結果を丸ごと保持
	std::vector<_T>         ap_eigenvalues;         // eig: 近似固有値のみ(既定)
	std::vector<std::vector<_T> > ap_eigenvectors;  // opt-in(spimats_options::keep_ap_eigenvectors)
	// --- 保証層の計算記録 ---
	//   stage(どこまで進んだか)、検証量(‖RA−I‖ 等)、破綻理由、
	//   is_spd の試行記録 …
	// >>> USER-IMPLEMENTS [SPI-U0: options/info 拡張(任意)] <<<
	//     正しさはオーナーが実装・保証。Claude の想定コードは補助であり
	//     拘束しない。spimats_options / spd_options 側のノブ拡張も本区画の
	//     扱いに準ずる(設計 §11: U0 = options / verify_info の拡張)。
	//
	// [SPI-D1 v1.1 実装] 設計 SPI-D1_design.md §3 の全数フィールド。
	// これは【計算の記録】であり対象行列への主張ではない(G3。判定の
	// 成否は is_spd の返値 / lss・veigs の result だけが担う)。
	// 型の規約: 記録は設計 §3 の指定どおり double 固定(_T ≠ double でも
	// double へ変換して記録する --- フライトレコーダの精度は保証に関与
	// しない)。書き込み側での変換は kv::dd / kv::mpfr の operator double
	// [確認済み kv/dd.hpp L728 / kv/mpfr.hpp L415]。

	// 進行段階: 0=未着手 / 1=前処理 / 2=ノルム評価 / 3=判定 / 4=成功
	int stage;

	// --- is_spd の試行記録(設計 §3)---
	double spd_shift;               // 採用したシフト量(fast: c / accurate: 0.8(s−M))
	double spd_rho_rad;             // 区間版: Collatz による ρ(rad) 上界 M_ρ
	std::vector<int> spd_scale_exp; // van der Sluis 2 冪スケーリング指数(適用時のみ非空)
	int spd_collatz_iters;          // 区間版: Collatz 反復回数
	bool spd_accurate_path;         // accurate 経路を通ったか
	// [SPI-F3 F3-9] 2 冪スケーリングの往復検査が不成立で scaling 全体を断念したか
	//(断念時は spd_scale_exp は空のまま)
	bool spd_scale_abandoned;

	// --- lss の記録(設計 §3)---
	double lss_g_norm_sup;          // ‖G‖∞ = Gnorm の上端
	double lss_t_norm_sup;          // 残差ノルム = Tnorm の上端
	std::size_t lss_rows_done;      // 行ストリーミングの完了行数

	// --- veigs の記録(設計 §3)---
	std::size_t veigs_bisection_steps;
	std::size_t veigs_clusters;
	// [SPI-F3 F3-6] 入口 2 冪正規化の記録: {a, b}(A ← 2^-a·A、B ← 2^-b·B。
	// 断念時は {0, 0} と abandoned = true)。ap_eigenvalues は正規化後の値のまま。
	std::vector<int> veigs_scale_exp;
	bool veigs_scale_abandoned;
	// [SPI-F3 F3-4 / F3-5] 返した各クラスタの確定番号(1 始まり両端含む、出力順。
	// veigs.m の ind_range に対応)と Rayleigh 退路の記録。
	std::vector<std::size_t> veigs_index_lo, veigs_index_hi;
	bool veigs_rayleigh_fallback;
	// [SPI-F3 F3-7 梯子] AP が収束した要求個数(k_req → min(n,k+2) → k の梯子の
	// どの段か。G-3 裁定で記録を要求。全段不収束なら 0)
	std::size_t veigs_ap_request_count;
	// [SPI-F4 F4-2] 要求順位の検査結果(0 = 未検査 / 1 = 証明済み / −1 = 証明できず)と、
	// 検査で追加した numeric 分解(sep_eval 呼び出し)の回数
	int veigs_rank_certified;
	std::size_t veigs_rank_check_factorizations;
	// >>> END [SPI-U0] <<<

	void clear() {
		ap_lss_result = linear_solve_result<_T>();
		ap_eigenvalues.clear();
		ap_eigenvectors.clear();
		// [SPI-D1 v1.1] SPI-U0 拡張フィールドの初期化。clear() 本体は U0
		// 区画外だが、初期化行の追加は「U0 拡張フィールドの扱いに準ずる」
		// --- 骨組み U0 注記の範囲内であり、これが唯一の区画外接触
		// (設計 SPI-D1_design.md §3 の許容範囲注記の転記)。
		stage = 0;
		spd_shift = 0.0;
		spd_rho_rad = 0.0;
		spd_scale_exp.clear();
		spd_collatz_iters = 0;
		spd_accurate_path = false;
		spd_scale_abandoned = false;
		lss_g_norm_sup = 0.0;
		lss_t_norm_sup = 0.0;
		lss_rows_done = 0;
		veigs_bisection_steps = 0;
		veigs_clusters = 0;
		veigs_scale_exp.clear();
		veigs_scale_abandoned = false;
		veigs_index_lo.clear();
		veigs_index_hi.clear();
		veigs_rayleigh_fallback = false;
		veigs_ap_request_count = 0;
		veigs_rank_certified = 0;
		veigs_rank_check_factorizations = 0;
	}
};

// ---------------------------------------------------------------------------
// spimats<_T, _AP, _DP>: 疎行列の区間型(精度保証)ポリシー(設計 §2.1)。
//   _T  : 点スカラー型(double / kv::dd / kv::mpfr<N> 等。double 限定でない)
//   _AP : 近似(点)疎ポリシー(spmats<_T> 派生。spumar は T = double 限定 ---
//         spumar 側の static_assert が要素挿入時に弾く [確認済み spumar.hpp
//         L253 相当])
//   _DP : 区間密ポリシー(要素 kv::interval<_T>、eigsymge を提供。既定
//         vcp::imats<_T> = imats<_T, mats<_T>>。vcp::pidblas で BLAS 高速版に
//         切り替え --- 派生側オーバーライド込み)。区間密行列は
//         matrix<kv::interval<_T>, _DP> を内部で構成(DD-11、v0.2)
// index_type は _AP から導出する(_Index を変更したい利用者は
// AP = spmats<T, long long> を渡す)。
// デフォルト引数は左から埋まるため(C++ の規則)、DP のみ変えたい場合も
// AP の明示が必要: spimats<double, vcp::spmats<double>, vcp::pdblas>。
// 利便のため alias template spimats_dp を併設(ファイル末尾)。
// ---------------------------------------------------------------------------
template <typename _T,
          class _AP = vcp::spmats<_T>,             // 近似(点)疎ポリシー
          class _DP = vcp::imats<_T> >             // 区間密ポリシー(DD-11)
class spimats : public vcp::spmats<kv::interval<_T>,
                                   typename _AP::index_type> {
public:
	typedef typename _AP::index_type index_type;   // [確認済み spmats.hpp L32]
	typedef vcp::spmats<kv::interval<_T>, index_type> base_type;
	typedef kv::interval<_T> interval_type;

	static_assert(std::is_base_of<vcp::spmats<_T, typename _AP::index_type>, _AP>::value,
		"vcp::spimats: _AP must derive from vcp::spmats<_T, _AP::index_type> (point sparse policy over the same scalar)");

	// 保証層のノブ(spumar の muar_opt 方式。spmatrix ラッパー越しには
	// 見えない --- ポリシーオブジェクト直接利用で設定する)
	spimats_options spi_opt;

	// フライトレコーダ(設計 §5.4。const _impl override から書くため
	// mutable。thread-safe ではない)
	mutable spimats_verify_info<_T> last_verify_info;

	spimats() : base_type(), spi_opt(), last_verify_info() {}

	// spmatrix の算術演算子(add / sub / mul / スカラー演算)は、結果の
	// ポリシーへ `static_cast<_P&>(C) = policy_add(...)` の形で基底
	// spmats<interval<_T>> 値を代入する [確認済み spmatrix.hpp L446 相当]。
	// 派生ポリシーはこの基底→派生代入を受ける口が必要なため、基底部分
	// オブジェクトへの代入として定義する(spi_opt / last_verify_info は
	// 変更しない --- spmats::normtwo 判断C「基底部分オブジェクト代入で
	// 動的型を保存」と同じ規律。通常のコピー代入は暗黙定義のまま)。
	spimats& operator=(const base_type& A) {
		base_type::operator=(A);
		return *this;
	}

	// -------------------------------------------------------------------
	// [SPI-F3 F3-1 / F3-2] 前処理付き固有値入口(設計 SPI-F3_design_v1.0.md §1-1)
	//
	// 基底 spmats の M 付き 4 本は非 virtual テンプレート [確認済み
	// spmats_eigs.hpp L4517-4526 / L4551-4562 / L4634-4646 / L4698-] で、
	// 派生の virtual _impl を通らず区間 Lanczos の値をそのまま成功で返す
	// (novel レビュー N1: 保証の迂回)。ここで同名テンプレートを定義して
	// 基底版を隠し、M は【近似段(AP)にだけ】渡して保証段は M なし版と
	// 同一のコアへ合流させる。M なし 4 本は using で再公開する。
	// (1) using 4 行
	// -------------------------------------------------------------------
	using base_type::policy_eigs_with_info;
	using base_type::policy_eigs;
	using base_type::policy_generalized_eigs_with_info;
	using base_type::policy_generalized_eigs;

	// (2) 前処理の型判定と adapter(F3-2)。
	// spi_prec_accepts_point_<Prec>: M.apply(const vector<_T>&, vector<_T>&) が
	//   式として成立すれば true(C++11 式 SFINAE)。点型 M(jacobi 等)は true、
	//   区間ベクトル用 apply しか持たない M は false。
	// spi_prec_adapter_<Prec>: AP(点疎)に渡す薄い包み。点型 M はそのまま、
	//   区間型 M は r を幅 0 の区間ベクトルに昇格して M.apply を呼び kv::mid で
	//   点に落とす。いずれも最後に z を scale(= 2^a、F3-6 の正規化に対する
	//   M の尺度合わせ。正規化前は 1)倍する。M は近似段の前処理でしかなく
	//   保証には一切関与しない(adapter の算術は点演算でよい)。
	template <class Prec>
	struct spi_prec_accepts_point_ {
		template <class P>
		static auto test_(int) -> decltype(
			std::declval<const P&>().apply(std::declval<const std::vector<_T>&>(),
			                               std::declval<std::vector<_T>&>()),
			std::true_type());
		template <class P>
		static std::false_type test_(...);
		static const bool value = decltype(test_<Prec>(0))::value;
	};
	template <class Prec>
	struct spi_prec_adapter_ {
		const Prec* M;
		_T scale;
		void apply_(const std::vector<_T>& r, std::vector<_T>& z, std::true_type) const {
			M->apply(r, z);
		}
		void apply_(const std::vector<_T>& r, std::vector<_T>& z, std::false_type) const {
			std::vector<interval_type> ri(r.size()), zi;
			for (std::size_t i = 0; i < r.size(); i++) ri[i] = interval_type(r[i]);
			M->apply(ri, zi);
			z.resize(zi.size());
			for (std::size_t i = 0; i < zi.size(); i++) z[i] = mid(zi[i]);
		}
		void apply(const std::vector<_T>& r, std::vector<_T>& z) const {
			apply_(r, z, std::integral_constant<bool, spi_prec_accepts_point_<Prec>::value>());
			for (std::size_t i = 0; i < z.size(); i++) z[i] = z[i] * scale;
		}
	};

	// (3) M 付き 4 本。with_info 2 本は入口検査の後に spi_veigs_dispatch_ で
	//   保証コアへ合流する。strict 2 本は基底 spmats_eigs.hpp L4635-4646 /
	//   L4698- と同じ converged・個数の検査(文言のみ spimats)。
	// spi_veigs_dispatch_: 本体 vcp/spimats.hpp との差はこの関数本体の 1 行
	//   (本体は spi_veigs_core_ を持たないので (void)M; で受けて
	//   policy_generalized_eigs_with_info_impl(B, k, opt) を呼ぶ。U5 昇格時に
	//   debug と同文へ置き換える)。
	template <class Prec>
	eig_result<interval_type> spi_veigs_dispatch_(const base_type& B, std::size_t k,
	                                              const eig_options<interval_type>& opt,
	                                              const Prec* M) const {
		spi_veigs_entry_checks_(B, k, opt); return spi_veigs_core_(B, k, opt, M);
	}
	template <class Prec>
	eig_result<interval_type> policy_eigs_with_info(std::size_t k,
	                                                const eig_options<interval_type>& opt,
	                                                const Prec& M) const {
		if (this->rowsize() != this->columnsize()) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_eigs_with_info: matrix must be square");
		}
		if (!this->is_finalized()) this->finalize();
		// B = I を厳密構成(SPI-R3)して一般化保証コアへ(override #2 と同文)
		base_type IB;
		vcp::identity_interval(this->rowsize(), IB);
		return spi_veigs_dispatch_(IB, k, opt, &M);
	}
	template <class Prec>
	eig_result<interval_type> policy_generalized_eigs_with_info(const base_type& B, std::size_t k,
	                                                            const eig_options<interval_type>& opt,
	                                                            const Prec& M) const {
		if (!this->is_finalized()) this->finalize();
		if (!B.is_finalized()) B.finalize();   // 防御(入口検査でも finalize する)
		return spi_veigs_dispatch_(B, k, opt, &M);
	}
	template <class Prec>
	std::vector<interval_type> policy_eigs(std::size_t k, const eig_options<interval_type>& opt,
	                                       const Prec& M) const {
		const eig_result<interval_type> result = policy_eigs_with_info(k, opt, M);
		if (!result.converged)
			vcp::throw_error<vcp::state_error>("spimats::policy_eigs(with preconditioner): eigensolver did not converge");
		if (result.returned_real_count < k)
			vcp::throw_error<vcp::state_error>("spimats::policy_eigs(with preconditioner): insufficient eigenvalues returned");
		return result.eigenvalues;
	}
	template <class Prec>
	std::vector<interval_type> policy_generalized_eigs(const base_type& B, std::size_t k,
	                                                   const eig_options<interval_type>& opt,
	                                                   const Prec& M) const {
		const eig_result<interval_type> result = policy_generalized_eigs_with_info(B, k, opt, M);
		if (!result.converged)
			vcp::throw_error<vcp::state_error>("spimats::policy_generalized_eigs(with preconditioner): eigensolver did not converge");
		const std::size_t n = static_cast<std::size_t>(this->rowsize());
		const std::size_t k_eff = (k < n) ? k : n;
		if (result.returned_real_count < k_eff)
			vcp::throw_error<vcp::state_error>("spimats::policy_generalized_eigs(with preconditioner): insufficient eigenvalues returned");
		return result.eigenvalues;
	}
	// >>> END [SPI-F3 §1-1] <<<

	// -------------------------------------------------------------------
	// is_spd(SPD 証明述語、設計 §7)
	//
	// bool のみを返す(Rump isspd.m 完全パリティ)。最小固有値の下界は
	// 返さない: 下界が欲しい利用者は、保証付き区間演算で A − cI / A − cB を
	// 構成して is_spd に問い合わせる(二分法の原始述語)。veigs(SPI-U5)の
	// シフト選択もこの二分法を使う想定。
	//
	// 非対称な意味論(ヘッダ契約):
	//   true  = 区間内のすべての対称行列が正定値であることが証明された
	//   false = 十分条件不成立。「SPD でない」という主張ではない(無主張)
	// 述語であり、失敗(false)はシフト探索など正常フローの一部なので
	// throw しない。misuse(非正方・区間版での厳密非対称 [nosymcheck で
	// 回避可])は従来どおり throw。
	// -------------------------------------------------------------------

	// 区間版: *this が対象(NVI 外殻 --- SPI-G1 D-5。判定対象を引数に取る
	// オーバーロードへ転送する)。
	bool is_spd(const spd_options<_T>& opt = spd_options<_T>()) const {
		return is_spd(static_cast<const base_type&>(*this), opt);
	}
	// 区間版: 判定対象 IA を引数に取る(新設 --- veigs の check_spd が使う)。
	// 外殻(非仮想)の責務: 正方検査(throw dimension_error)・finalize・
	// is_spd 所有フィールドのクリア・厳密対称入口検査(SPI-R4、nosymcheck で
	// 回避、throw domain_error)。判定本体は protected virtual の _impl
	// (派生ポリシー spidmats が override する口 --- veigs からも *this 経由で
	// 仮想呼び出しされる)。
	bool is_spd(const base_type& IA, const spd_options<_T>& opt = spd_options<_T>()) const {
		if (IA.rowsize() != IA.columnsize())
			vcp::throw_error<vcp::dimension_error>("spimats::is_spd: matrix must be square");
		if (!IA.is_finalized()) IA.finalize();
		spi_clear_spd_fields_();
		if (IA.rowsize() == 0) {                 // [SPI-F1 F1-1] n = 0: 空行列は空虚に正定値
			last_verify_info.stage = 4;
			last_verify_info.spd_accurate_path = true;
			return true;
		}
		if (!opt.nosymcheck && !policy_is_symmetric(IA))
			vcp::throw_error<vcp::domain_error>("spimats::is_spd: matrix is not exactly symmetric as an interval matrix (endpoint comparison; set spd_options::nosymcheck to skip)");
		return is_spd_interval_impl_(IA, opt);
	}
	// 点版: 点疎行列が対象(_AP は spmats<_T> 派生なのでそのまま渡せる)。
	bool is_spd(const spmats<_T, index_type>& A, const spd_options<_T>& opt = spd_options<_T>()) const {
		if (A.rowsize() != A.columnsize())
			vcp::throw_error<vcp::dimension_error>("spimats::is_spd: matrix must be square");
		if (!A.is_finalized()) A.finalize();
		spi_clear_spd_fields_();
		if (A.rowsize() == 0) {                  // [SPI-F1 F1-1] n = 0: 空行列は空虚に正定値
			last_verify_info.stage = 4;
			last_verify_info.spd_accurate_path = true;
			return true;
		}
		return is_spd_point_impl_(A, opt);
	}

protected:
	// -------------------------------------------------------------------
	// is_spd_interval_impl_: 区間版の判定本体(protected virtual、SPI-G1 D-5)。
	// 前処理(midrad 分解・対角検査・対角置換)を行って中核 spi_is_spd_core_
	// を呼ぶ(設計 SPI-G1 §2.2〜2.3)。派生 spidmats は本関数を override して
	// double 専用の fast 経路を差し込む。
	// -------------------------------------------------------------------
	virtual bool is_spd_interval_impl_(const base_type& IA, const spd_options<_T>& opt) const {
		// >>> USER-IMPLEMENTS [SPI-U1: is_spd 区間版コア(スケーリング・ρ(rad) 上界・シフト・判定を含む)] <<<
		//     正しさはオーナーが実装・保証。
		//
		// ===============================================================
		// [SPI-G1 実装] U1: is_spd 区間版 --- T 汎用の事後評価(accurate)に一本化
		//   (D-1 / D-2)。前処理(本区画)→ 中核 spi_is_spd_core_(点版と共通、
		//   D-6)。旧 fast 経路(Rump 2006 の事前誤差解析 ε = 2⁻⁵³ / η = 2⁻¹⁰⁷⁴、
		//   double 限定)と R5 double 帰着(D-7)は削除。fast は派生 spidmats
		//   (別トラック)で復活させる。
		// 規範出典(設計 SPI-G1_design_v1.1.md §1.6 の表を転記):
		//   - 近似 Cholesky G の対角が非負なら GGᵀ ⪰ 0、E := LLᵀ − λ̃²I − GGᵀ の
		//     評価で λ(LLᵀ) ≥ λ̃² − ρ(E)、ρ(E) ≤ ‖E‖∞ --- Rump 1993, Computing
		//     Suppl. 9, §2 Theorem 2.1 とその直前の議論(Golub–Van Loan Cor. 8.1.3
		//     経由。Cor. 8.1.3 本文は [UNVERIFIED: 逐語未確認])。
		//   - SPD 行列への直接適用 A ≈ G̃G̃ᵀ、σ̃ ≈ σn(A)、λ̃ = 0.9σ̃ で
		//     Algorithm 2.1 → σn(A) ≥ (λ̃² − emax)^{1/2} --- Rump 1993 §2 II)。
		//   - 区間系: ΔA := |[A] − Ã|、Algorithm 2.2 手順 6 の µ := σ₁σ₂ −
		//     半径ノルム > 0 --- Rump 1993 Theorem 2.4 / Algorithm 2.2。
		//   - 同定理の 1999 年版(s := 0.9t、s − δ > 0 ⟹ PD、σn(A) ≥ s − ‖Δ‖)
		//     --- Rump 1999, Bulgak–Zenger 編, §4 Algorithm 4.1。
		//   - 対称 A′ ∈ ⟨A,R⟩ で λν(A′) ≥ λν(A) − ‖R‖、‖R‖ = ρ(R)、Collatz
		//     ρ(R) ≤ max_i (Rx)_i/x_i --- Rump 2006, BIT 46, p.440 (2.12) と Lemma 2.6。
		//   - 対角を下端 inf に固定した A_ が SPD なら [A] の全対称行列が SPD
		//     (A′ = A″ + D₊、D₊ ⪰ 0 対角)--- isspd.m L81-82 コメント(自明)。
		//   - 2 冪対角合同 DAD・対称置換 PᵀAP での正定値性不変 --- 標準結果(自明)。
		//   - 係数 0.8(isspd.m L211)/ 逆反復停止 1e-4(L208)/ 再試行 3 回
		//     (Rump 1993 Alg. 2.2)--- 自由パラメータ(成功率のみ、主張に無関係)。
		//   - isspd.m の accurate + 区間分岐(L210-217、判定に M を含めない)は
		//     移植しない(D-3。sandbox/docs/issues/SPI-G1_isspd_interval_accurate_issue.md)。
		// 丸め規律: 丸めモード操作は書かない。isspd.m の setround(±1) は kv 区間
		//   演算 + 端点取り出しで表現(現行規律のまま)。
		// ===============================================================
		// ---- accurate 経路(SPI-G1 D-2 / D-3 / D-7。前処理は設計 §2.2)----
		// [A] = <M, R>(R1 midrad_split、kv::midrad のみ)→ 対角検査 → M の
		// 対角を diag([A]).inf に置換・R の対角を 0 に【対応: isspd.m L110-111。
		// 根拠(L81-82 コメント): 任意の対称 A′ ∈ [A] は A′ = A″ + D₊(A″ は
		// 対角を下端 inf に固定した [A_] の元、D₊ ⪰ 0 対角)なので
		// isspd(A_) ⟹ isspd(A)】→ n == 1 → 中核 spi_is_spd_core_(M, &R, dlow)。
		// _T のまま判定する(R5 double 帰着なし)。
		{
			typedef spmats<_T, index_type> PT_;
			const index_type n = IA.rowsize();
			last_verify_info.stage = 1;                      // 前処理段
			last_verify_info.spd_accurate_path = true;
			PT_ Mm, Rr;
			vcp::midrad_split(IA, Mm, Rr);
			// 対角検査【対応: isspd.m L76-85】: mid(d_i) ≤ 0 → false【L77-80】/
			// mig(d_i) ≤ 0 → false【L81-85】。いずれも失敗側ゲート !(x > 0) で
			// NaN は失敗側へ。非格納対角 = [0,0] は mid = 0 で false。midrad /
			// mig は kv::interval の friend(ADL で解決 --- 修飾名は不可
			// [確認済み kv/interval.hpp L622/L644])。
			std::vector<_T> dlow(static_cast<std::size_t>(n), _T(0));
			{
				std::vector<interval_type> dv(static_cast<std::size_t>(n), interval_type(_T(0)));
				const base_type Ac = IA.as_csr();
				const std::vector<index_type>& outerv = Ac.outer_index();
				const std::vector<index_type>& innerv = Ac.inner_index();
				const std::vector<interval_type>& valv = Ac.values();
				for (index_type i = 0; i < n; i++) {
					for (index_type p = outerv[static_cast<std::size_t>(i)];
					     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
						if (innerv[static_cast<std::size_t>(p)] == i) {
							dv[static_cast<std::size_t>(i)] = valv[static_cast<std::size_t>(p)];
						}
					}
				}
				for (index_type i = 0; i < n; i++) {
					_T dm, drr;
					midrad(dv[static_cast<std::size_t>(i)], dm, drr);
					if (!(dm > _T(0))) return false;         // 【対応: isspd.m L77-80】
					if (!(mig(dv[static_cast<std::size_t>(i)]) > _T(0))) return false;   // 【L81-85】
					dlow[static_cast<std::size_t>(i)] = dv[static_cast<std::size_t>(i)].lower();
				}
			}
			// M の対角を dlow に置換(dlow > 0 なので必ず格納)、R の対角を落とす
			// (radius 0 = 非格納)。値移動のみ・算術なし。
			{
				PT_ M2, R2;
				M2.resize(n, n);
				R2.resize(n, n);
				{
					const PT_ Mc = Mm.as_csr();
					const std::vector<index_type>& outerv = Mc.outer_index();
					const std::vector<index_type>& innerv = Mc.inner_index();
					const std::vector<_T>& valv = Mc.values();
					for (index_type i = 0; i < n; i++) {
						for (index_type p = outerv[static_cast<std::size_t>(i)];
						     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
							const index_type j = innerv[static_cast<std::size_t>(p)];
							if (i != j) M2.add(i, j, valv[static_cast<std::size_t>(p)]);
						}
					}
				}
				for (index_type i = 0; i < n; i++) {
					M2.add(i, i, dlow[static_cast<std::size_t>(i)]);
				}
				{
					const PT_ Rc = Rr.as_csr();
					const std::vector<index_type>& outerv = Rc.outer_index();
					const std::vector<index_type>& innerv = Rc.inner_index();
					const std::vector<_T>& valv = Rc.values();
					for (index_type i = 0; i < n; i++) {
						for (index_type p = outerv[static_cast<std::size_t>(i)];
						     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
							const index_type j = innerv[static_cast<std::size_t>(p)];
							if (i != j) R2.add(i, j, valv[static_cast<std::size_t>(p)]);
						}
					}
				}
				M2.finalize();
				R2.finalize();
				Mm = M2;
				Rr = R2;
			}
			// [SPI-F1 F1-2] 非有限端点ガード: M・R の非対角格納値に NaN / ±∞ が
			// 1 つでもあれば前処理段のまま false(無主張)。対角は上の対角検査
			// (mid・mig > 0 の失敗側ゲート)に任せ、diag [v, +∞] (v > 0) の true を
			// 保つ(G-0 裁定: 走査位置は midrad_split 直後ではなく対角置換の直後)。
			// 判定は比較演算のみ(!(v − v == 0) で NaN・±∞ を失敗側へ。T 汎用、
			// double cast / isfinite は書かない)。
			if (!spi_offdiag_all_finite_(Mm) || !spi_offdiag_all_finite_(Rr)) return false;
			// n == 1【対応: isspd.m L88-92】: 対角検査通過 = 下端 inf > 0 なので
			// [A] 内の全(1×1 対称)行列が正定値。
			if (n == 1) {
				last_verify_info.stage = 4;
				return true;
			}
			return spi_is_spd_core_(Mm, &Rr, dlow, opt);
		}
		// >>> END [SPI-U1] <<<
		return false;   // not reached
	}

	// -------------------------------------------------------------------
	// is_spd_point_impl_: 点版の判定本体(protected virtual、SPI-G1 D-5)。
	// 対角抽出・対角検査を行って中核 spi_is_spd_core_(R = nullptr)を呼ぶ。
	// 型ガードなし: _T = double / kv::dd / kv::mpfr<N> で同じ経路。
	// -------------------------------------------------------------------
	virtual bool is_spd_point_impl_(const spmats<_T, index_type>& A, const spd_options<_T>& opt) const {
		// >>> USER-IMPLEMENTS [SPI-U2: is_spd 点版コア(fast / accurate)] <<<
		//     正しさはオーナーが実装・保証。
		//
		// ===============================================================
		// [SPI-G1 実装] U2: is_spd 点版 --- T 汎用の事後評価(accurate)に一本化
		//   (D-1 / D-2)。前処理(本区画: 対角検査・n == 1)→ 中核
		//   spi_is_spd_core_(区間版と共通、R = nullptr で M_ρ = 0、D-6)。
		//   旧 fast 経路(Rump 2006 の事前誤差解析: t プロファイル計数・
		//   α = (t+3)ε・η 項・次元閾値の定数一式、double 限定
		//   の型ガード)は削除。fast は派生 spidmats(別トラック)で復活させる。
		// 規範出典: U1 区画冒頭の表と同一(Rump 1993 Thm 2.1 / Alg. 2.2、Rump 1999
		//   Alg. 4.1。点版は Δ = 0、M_ρ = 0 の特殊化)。isspd.m の accurate 分岐
		//   L192-226 に対応(判定 L217 の s > ‖CᵀC − A‖ と同値)。
		// 丸め規律・比較規約: 成功宣言側は certified 比較(x > c の成立)、失敗側
		//   ゲートは !(x > c) イディオム(NaN・∞ は失敗側に落ちる)。丸めモード
		//   操作は書かない。
		// ===============================================================
		{
			// ---- accurate 経路(SPI-G1 D-2。前処理は設計 §2.2 点版)----
			// 対角 d_i = A_ii(非格納は厳密 0)。!(d_i > 0) → false(無主張。NaN も
			// 失敗側)【対応: isspd.m L76-85(点版: mid = mig = 値自身)】→ n == 1
			// → 中核 spi_is_spd_core_(A, nullptr, dg)(R なし = M_ρ = 0)。
			const index_type n = A.rowsize();
			last_verify_info.stage = 1;                    // 前処理段
			last_verify_info.spd_accurate_path = true;
			std::vector<_T> dg(static_cast<std::size_t>(n), _T(0));
			{
				const spmats<_T, index_type> Ac = A.as_csr();
				const std::vector<index_type>& outerv = Ac.outer_index();
				const std::vector<index_type>& innerv = Ac.inner_index();
				const std::vector<_T>& valv = Ac.values();
				for (index_type i = 0; i < n; i++) {
					for (index_type p = outerv[static_cast<std::size_t>(i)];
					     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
						if (innerv[static_cast<std::size_t>(p)] == i) {
							dg[static_cast<std::size_t>(i)] = valv[static_cast<std::size_t>(p)];
						}
					}
				}
			}
			for (index_type i = 0; i < n; i++) {
				if (!(dg[static_cast<std::size_t>(i)] > _T(0))) return false;
			}
			// [SPI-F1 F1-2] 非有限ガード(点版): 非対角格納値に NaN / ±∞ があれば
			// 前処理段のまま false(無主張)。対角は上の !(d_i > 0) に任せる。
			if (!spi_offdiag_all_finite_(A)) return false;
			// n == 1【対応: isspd.m L88-92】: A_11 > 0(厳密比較)なので SPD。
			if (n == 1) {
				last_verify_info.stage = 4;                // 成功
				return true;
			}
			return spi_is_spd_core_(A, nullptr, dg, opt);
		}
		// >>> END [SPI-U2] <<<
		return false;   // not reached
	}

	// >>> USER-IMPLEMENTS [SPI-U2c: is_spd 共通中核(点・区間共通、spi_is_spd_core_)] <<<
	// -------------------------------------------------------------------
	// [SPI-F1] 静的補助 2 本(判定式・シフト・残差には関与しない防御ガード)。
	// spi_offdiag_all_finite_: 点疎行列 X の非対角格納値がすべて有限か
	//   (F1-2。!(v − v == 0) が成立する NaN / ±∞ を 1 つでも含めば false。
	//   比較演算のみ、T 汎用)。
	// spi_perm_is_bijection_: 順序器の置換が {0..n−1} の全単射か(F1-4。長さ n・
	//   値域 [0, n)・重複なし。chol カーネルの sparse_chol_verify_permutation_ と
	//   同基準)。
	// -------------------------------------------------------------------
	static bool spi_offdiag_all_finite_(const spmats<_T, index_type>& X) {
		const spmats<_T, index_type> Xc = X.as_csr();
		const std::vector<index_type>& outerv = Xc.outer_index();
		const std::vector<index_type>& innerv = Xc.inner_index();
		const std::vector<_T>& valv = Xc.values();
		const index_type n = Xc.rowsize();
		for (index_type i = 0; i < n; i++) {
			for (index_type p = outerv[static_cast<std::size_t>(i)];
			     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
				if (innerv[static_cast<std::size_t>(p)] == i) continue;
				const _T& v = valv[static_cast<std::size_t>(p)];
				if (!(v - v == _T(0))) return false;
			}
		}
		return true;
	}
	static bool spi_perm_is_bijection_(const std::vector<index_type>& perm, index_type n) {
		if (n < index_type(0)) return false;
		const std::size_t un = static_cast<std::size_t>(n);
		if (perm.size() != un) return false;
		std::vector<char> seen(un, char(0));
		for (std::size_t i = 0; i < un; i++) {
			const index_type q = perm[i];
			if (q < index_type(0) || q >= n) return false;                  // 値域外
			if (seen[static_cast<std::size_t>(q)] != char(0)) return false;   // 重複
			seen[static_cast<std::size_t>(q)] = char(1);
		}
		return true;
	}

	// -------------------------------------------------------------------
	// spi_is_spd_core_: 点・区間共通の証明中核(非仮想、SPI-G1 D-2 / D-6)。
	// R == nullptr が点版(M_ρ = 0 の特殊化)。
	// 規範出典: Rump 1993 (Computing Suppl. 9, §2) Theorem 2.1 / Algorithm 2.2、
	//   Rump 1999 (Bulgak–Zenger 編, §4) Algorithm 4.1、Rump 2006 (BIT 46,
	//   p.440) (2.12) と Lemma 2.6(Collatz)。INTLAB isspd.m
	//   (sandbox/reference/isspd.m)の区間前処理 L107-190 と accurate 分岐
	//   L192-226 に対応(区間分岐の判定式は移植しない --- D-3)。
	// 主張の連鎖(設計 SPI-G1 §2.3): 任意の対称 A′ ∈ [A] について
	//   D P A′ Pᵀ D = M + Δ、|Δ| ≤ R(スケール・置換後)。‖Δ‖₂ ≤ ρ(|Δ|) ≤
	//   ρ(R) ≤ M_ρ(Perron–Frobenius の単調性 + Collatz)。M = A_s + s·I + N
	//   (N ⪰ 0 対角)、A_s = L₂L₂ᵀ − E、L₂L₂ᵀ ⪰ 0 ⟹ λmin(M) ≥ s − ‖E‖₂ ≥
	//   s − ‖E‖∞(E は端点まで対称)。判定 s − ‖E‖∞ > M_ρ ⟹ λmin(A′) ≥
	//   λmin(M) − M_ρ > 0。2 冪対角合同と対称置換は正定値性を保つ。点版は
	//   Δ = 0、M_ρ = 0 で同じ連鎖(Rump 1993 Thm 2.1 / 1999 Alg. 4.1)。
	// T 汎用性: 数値計算に依存するのは「E ⊇ L₂L₂ᵀ − A_s」「M_ρ ≥ ρ(R)」
	//   「A_s の対角 ≤ M_ii − s」の 3 点のみで、いずれも interval<_T> の演算の
	//   包含性しか使わない(誤差モデル不使用)。L, L₂, t, s は近似で、精度に
	//   関する仮定を置かない。
	// 丸め規律: 丸めモード操作は書かない。isspd.m の setround(±1) は kv 区間
	//   演算 + 端点取り出し(.upper() / .lower())で表現する。保証に関与する
	//   式に double キャスト・isfinite・fabs(double) を書かない。近似段
	//   (スケーリング指数選択)のみ double 経由(行末 approximation only)。
	// 比較規約: 成功宣言側は certified 比較(x > c の成立)、失敗側ゲートは
	//   !(x > c) イディオム(NaN・∞ は失敗側に落ちる)。
	// -------------------------------------------------------------------
	bool spi_is_spd_core_(const spmats<_T, index_type>& M,
	                      const spmats<_T, index_type>* R,
	                      const std::vector<_T>& diag_low,
	                      const spd_options<_T>& opt) const {
		typedef kv::interval<_T> IT_;
		typedef spmats<_T, index_type> PT_;
		using std::abs;                                    // kv::dd / kv::mpfr は ADL の abs
		const index_type n = M.rowsize();
		last_verify_info.stage = 1;                        // 前処理段
		last_verify_info.spd_accurate_path = true;         // D-4: 証明経路は常に accurate
		std::vector<_T> dg(diag_low);                      // 対角(スケール・置換に追従)
		PT_ Awork = M.as_csr();                            // 以降の判定対象(変形はコピーに)
		const bool hasR = (R != nullptr);
		PT_ Rwork;
		if (hasR) Rwork = R->as_csr();

		// --- 1. van der Sluis スケーリング【対応: isspd.m L116-136】(D-9)---
		// 条件 maxdiag/mindiag > sqrt(n)、ガード maxdiag > 1e100 / mindiag <
		// 1e-100 では適用しない【L119】。適用判定と指数選択は double 近似で
		// よい(成功率のみに影響し、主張は後段が担う)。ただし D_ii = 2^e の
		// _T 表現と積 a_ij·D_ii·D_jj の厳密性は主張の前提(2 冪倍はガード
		// 範囲で丸めなし: dd は上位語 2^e・下位語 0、mpfr は厳密。受け入れ
		// テスト T-6 で dd の成分ビット検証)。[SPI-F3 F3-9] 非対角成分が
		// subnormal に落ちる・0 に落ちる・inf になる場合は積が厳密でなく、
		// R が真の半径の上界でなくなり得る(audit 2026-10-15 §5: r = 2^-1074、
		// 対角 [16, 1/16] で R の成分が 0 に落ちる)。そこで各積 sv = x·D_ii·D_jj
		// に往復検査 (sv·D_jj⁻¹)·D_ii⁻¹ == x(2 冪の除算は逆元 2^-e との点乗算)
		// を課し、1 成分でも不成立なら scaling 全体を断念して Awork / Rwork /
		// dg を元のまま使う(spd_scale_exp は空、spd_scale_abandoned = true)。
		// 成立すれば全成分が厳密な 2 冪倍であり、主張の連鎖は無変更。R があれば
		// R ← D R D(非負 × 正の 2 冪で非負性保存【L133-135】)。dlow も追従(厳密)。
		{
			double maxdiag = static_cast<double>(dg[0]);   // approximation only
			double mindiag = maxdiag;                      // approximation only
			for (index_type i = 1; i < n; i++) {
				const double di = static_cast<double>(dg[static_cast<std::size_t>(i)]);   // approximation only
				if (di > maxdiag) maxdiag = di;
				if (di < mindiag) mindiag = di;
			}
			if ((maxdiag / mindiag > std::sqrt(static_cast<double>(n)))
			    && !((maxdiag > 1e100) || (mindiag < 1e-100))) {   // approximation only
				// D_ii = 2^(-ceil(0.5*log2 d_i))【L127, L130】。指数を int に確定
				// してから ldexp で 2 冪を厳密構成する(log2/ceil の近似誤差は
				// 指数選択を ±1 ずらすだけで「D が 2 冪 = 乗算厳密」は崩れない)。
				std::vector<int> sexp(static_cast<std::size_t>(n));
				std::vector<_T> dsc(static_cast<std::size_t>(n)), dinv(static_cast<std::size_t>(n));
				for (index_type i = 0; i < n; i++) {
					const double di = static_cast<double>(dg[static_cast<std::size_t>(i)]);   // approximation only
					const int e = -static_cast<int>(std::ceil(0.5 * std::log2(di)));         // approximation only
					sexp[static_cast<std::size_t>(i)] = e;
					dsc[static_cast<std::size_t>(i)] = _T(std::ldexp(1.0, e));    // 2 冪の厳密表現
					dinv[static_cast<std::size_t>(i)] = _T(std::ldexp(1.0, -e));  // 逆元 2^-e(厳密)
				}
				// X ← D X D【L131 / L134】: 格納値 x_ij·D_ii·D_jj。[SPI-F3 F3-9] 各積に
				// 往復検査を課し、不成立(subnormal 落ち・0 落ち・inf)なら false
				// (X は更新しない --- Y を作ってから呼び出し側で代入する)。
				auto scale = [&](const PT_& X, PT_& Y) -> bool {
					const PT_ Xc = X.as_csr();
					const std::vector<index_type>& outerv = Xc.outer_index();
					const std::vector<index_type>& innerv = Xc.inner_index();
					const std::vector<_T>& valv = Xc.values();
					Y.resize(n, n);
					for (index_type i = 0; i < n; i++) {
						for (index_type p = outerv[static_cast<std::size_t>(i)];
						     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
							const index_type j = innerv[static_cast<std::size_t>(p)];
							const _T& x = valv[static_cast<std::size_t>(p)];
							const _T sv = x * dsc[static_cast<std::size_t>(i)]
								* dsc[static_cast<std::size_t>(j)];
							if (!((sv * dinv[static_cast<std::size_t>(j)]) * dinv[static_cast<std::size_t>(i)] == x)) return false;
							if (!(sv == _T(0))) Y.add(i, j, sv);
						}
					}
					Y.finalize();
					return true;
				};
				PT_ Ya, Yr;
				std::vector<_T> dg2(dg);
				bool sc_ok = scale(Awork, Ya);
				if (sc_ok && hasR) sc_ok = scale(Rwork, Yr);
				for (index_type i = 0; sc_ok && i < n; i++) {   // 対角も追従(往復検査つき)
					const _T x = dg[static_cast<std::size_t>(i)];
					const _T sv = x * dsc[static_cast<std::size_t>(i)] * dsc[static_cast<std::size_t>(i)];
					if (!((sv * dinv[static_cast<std::size_t>(i)]) * dinv[static_cast<std::size_t>(i)] == x)) sc_ok = false;
					dg2[static_cast<std::size_t>(i)] = sv;
				}
				if (sc_ok) {
					last_verify_info.spd_scale_exp = sexp;     // 試行記録(設計 §3)
					Awork = Ya;
					if (hasR) Rwork = Yr;
					dg = dg2;
				}
				else {
					last_verify_info.spd_scale_abandoned = true;   // 断念: 元の Awork / Rwork / dg のまま
				}
			}
		}

		// --- 2. ordering 前処理【対応: isspd.m L138-153(A と Arad の両方を
		//  置換)】---
		// 順列は先に行列へ適用し(値の移動のみ・算術なし = 厳密)、chol は常に
		// natural で実行する。M と R に同一の対称置換(SPD は合同変換 PᵀAP で
		// 不変、ρ(R) も対称置換で不変)。auto_select は基底 chol の規約
		// (auto → amd [確認済み tsparse_sparse_chol.hpp L81-82])に合わせる。
		{
			vcp::sparse_chol_ordering ord = opt.ordering;
			if (ord == vcp::sparse_chol_ordering::auto_select) {
				ord = vcp::sparse_chol_ordering::amd;
			}
			if (ord != vcp::sparse_chol_ordering::natural) {
				// 下三角(r >= c)CSC パターンを ordering 関数へ渡す [確認済み
				// tsparse_sparse_chol.hpp L363-371: 辺 {r,c} は内部で対称化]
				const PT_ Acc = Awork.as_csc();
				const std::vector<index_type>& outc = Acc.outer_index();
				const std::vector<index_type>& innc = Acc.inner_index();
				std::vector<index_type> lptr(static_cast<std::size_t>(n) + 1, index_type(0));
				std::vector<index_type> lind;
				lind.reserve(innc.size());
				for (index_type j = 0; j < n; j++) {
					lptr[static_cast<std::size_t>(j)] = static_cast<index_type>(lind.size());
					for (index_type p = outc[static_cast<std::size_t>(j)];
					     p < outc[static_cast<std::size_t>(j) + 1]; p++) {
						if (innc[static_cast<std::size_t>(p)] >= j) {
							lind.push_back(innc[static_cast<std::size_t>(p)]);
						}
					}
				}
				lptr[static_cast<std::size_t>(n)] = static_cast<index_type>(lind.size());
				std::vector<index_type> perm0;
				if (!vcp::sparse_chol_detail::sparse_chol_compute_ordering_(
					n, lptr, lind, ord, perm0)) {
					return false;   // 未解決 enum(基底と同じ防御)= 無主張
				}
				// [SPI-F1 F1-4] 置換の全単射検査(順序器の契約違反 = 内部状態の
				// 異常なので無主張ではなく state_error)。
				if (!spi_perm_is_bijection_(perm0, n)) {
					vcp::throw_error<vcp::state_error>("spimats::is_spd: ordering returned a non-bijective permutation (internal state error)");
				}
				std::vector<index_type> pinv(static_cast<std::size_t>(n));
				for (index_type k = 0; k < n; k++) {       // perm0 は new->old
					pinv[static_cast<std::size_t>(perm0[static_cast<std::size_t>(k)])] = k;
				}
				// Y(pinv[r], pinv[c]) = X(r, c)(値移動のみ)
				auto permute = [&](PT_& X) {
					const PT_ Xc = X.as_csc();
					const std::vector<index_type>& outv = Xc.outer_index();
					const std::vector<index_type>& innv = Xc.inner_index();
					const std::vector<_T>& valv = Xc.values();
					PT_ Y;
					Y.resize(n, n);
					for (index_type j = 0; j < n; j++) {
						for (index_type p = outv[static_cast<std::size_t>(j)];
						     p < outv[static_cast<std::size_t>(j) + 1]; p++) {
							const index_type r = innv[static_cast<std::size_t>(p)];
							Y.add(pinv[static_cast<std::size_t>(r)],
							      pinv[static_cast<std::size_t>(j)],
							      valv[static_cast<std::size_t>(p)]);
						}
					}
					Y.finalize();
					X = Y;
				};
				permute(Awork);
				if (hasR) permute(Rwork);
				std::vector<_T> dgp(static_cast<std::size_t>(n));
				for (index_type i = 0; i < n; i++) {       // 対角も追従
					dgp[static_cast<std::size_t>(pinv[static_cast<std::size_t>(i)])] =
						dg[static_cast<std::size_t>(i)];
				}
				dg.swap(dgp);
			}
		}

		// --- 3. M_ρ: ρ(R) の上界(Collatz 反復)【対応: isspd.m L155-190】---
		// 非負行列 R と正ベクトル x に対し ρ(R) ≤ max_i (Rx)_i / x_i(Collatz /
		// Perron–Frobenius --- Rump 2006 Lemma 2.6)。上界性は毎反復成立し、
		// 反復は上界を締めるだけなので停止条件(|M−m|/(m+M) < 0.1・上限 10
		// 【L184-186】)は品質のみ。y = R x を interval<_T> で評価し .upper() を
		// 採る(L161 setround(1) の区間置換)。比 y_i/x_i も区間除算の .upper()。
		// ゼロ行の除去【L165-174】: y_i = 0(上端厳密 0)はアクティブな列上で
		// 行 i が空 ⟹ 行/列 i を除去しても置換でブロック上三角になり、
		// スペクトルは {R′ の固有値} ∪ {0} --- ρ は保存される。R なし → M_ρ = 0。
		_T Mrho = _T(0);
		int citer = 0;
		if (hasR) {
			Mrho = _T(2);                                  // 【L159: M = 2】
			const PT_ Rcsr = Rwork.as_csr();
			const std::vector<index_type>& outerv = Rcsr.outer_index();
			const std::vector<index_type>& innerv = Rcsr.inner_index();
			const std::vector<_T>& valv = Rcsr.values();
			std::vector<char> act(static_cast<std::size_t>(n), 1);
			std::vector<_T> x(static_cast<std::size_t>(n), _T(1));   // 【L157: x = ones】
			std::vector<_T> y(static_cast<std::size_t>(n), _T(0));
			while (true) {                                 // 【L162: while 1】
				citer++;
				// y = R x の上界(アクティブ行×アクティブ列のみ)【L164】
				for (index_type i = 0; i < n; i++) {
					if (!act[static_cast<std::size_t>(i)]) continue;
					IT_ acc(_T(0));
					for (index_type p = outerv[static_cast<std::size_t>(i)];
					     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
						const index_type j = innerv[static_cast<std::size_t>(p)];
						if (!act[static_cast<std::size_t>(j)]) continue;
						acc += IT_(valv[static_cast<std::size_t>(p)])
							* IT_(x[static_cast<std::size_t>(j)]);
					}
					y[static_cast<std::size_t>(i)] = acc.upper();
				}
				// ゼロ成分の行/列除去【L165-174】
				bool removed = false;
				for (index_type i = 0; i < n; i++) {
					if (act[static_cast<std::size_t>(i)]
					    && y[static_cast<std::size_t>(i)] == _T(0)) {   // 厳密 0 判定(許容イディオム)
						act[static_cast<std::size_t>(i)] = 0;
						removed = true;
					}
				}
				bool any_active = false;
				for (index_type i = 0; i < n; i++) {
					if (act[static_cast<std::size_t>(i)]) { any_active = true; break; }
				}
				if (!any_active) {                         // 【L180-183: 空なら M = 0】
					Mrho = _T(0);
					break;
				}
				// m(下側の観測値)は停止判定専用の統計であり上界性に関与
				// しない(点計算でよい)【L171/L173/L177】。除去があれば m = 0。
				_T mlo = _T(0);
				bool mlo_set = removed;
				_T Mup = _T(0);
				_T ymax = _T(0);
				for (index_type i = 0; i < n; i++) {
					if (!act[static_cast<std::size_t>(i)]) continue;
					const std::size_t ui = static_cast<std::size_t>(i);
					// 比の上界(区間除算 + .upper())【L176/L178: x = y./x】
					const _T rup = (IT_(y[ui]) / IT_(x[ui])).upper();
					// 失敗側ゲート: 非負行列の比は ≥ 0。NaN は上界として採れない
					// ので無主張(∞ は M_ρ = ∞ として後段の !(t > M_ρ) で落ちる)
					if (!(rup >= _T(0))) return false;
					if (rup > Mup) Mup = rup;
					const _T rlo = y[ui] / x[ui];          // 統計(停止判定用)
					if (!mlo_set || rlo < mlo) { mlo = rlo; mlo_set = true; }
					if (y[ui] > ymax) ymax = y[ui];
				}
				Mrho = Mup;                                // 上界は毎反復有効
				// x = y / max(y)(正規化は近似でよい: x > 0 が保たれれば次反復も
				// Collatz 上界として成立)【L179】。ymax は除去後のアクティブ
				// 成分が非零・非負なので正(防御は失敗側)。
				if (!(ymax > _T(0))) return false;
				for (index_type i = 0; i < n; i++) {
					if (act[static_cast<std::size_t>(i)]) {
						x[static_cast<std::size_t>(i)] =
							y[static_cast<std::size_t>(i)] / ymax;
					}
				}
				// [SPI-F1 F1-3] 正規化の失敗側ゲート: アクティブ成分で !(x_i > 0)
				// (アンダーフローで 0、または NaN)なら次反復の区間除算が 0 除算に
				// なる経路なので無主張 false(Collatz の前提 x > 0 が崩れた)。
				for (index_type i = 0; i < n; i++) {
					if (act[static_cast<std::size_t>(i)]
					    && !(x[static_cast<std::size_t>(i)] > _T(0))) return false;
				}
				// 停止条件【L184-186: |M−m|/(m+M) < 0.1 または iter > 10】
				if ((abs(Mup - mlo) / (mlo + Mup) < _T(0.1)) || (citer > 10)) {
					break;
				}
			}
		}
		last_verify_info.spd_rho_rad = static_cast<double>(Mrho);   // 記録専用の double 変換(設計 §3)
		last_verify_info.spd_collatz_iters = citer;

		// --- chol オプション(近似 Cholesky・シフト後 Cholesky 共用)---
		// 実行主体は基底 spmats の自前 policy_chol_with_info(演算順管理下
		// [確認済み spmats.hpp L1955 相当])、ordering = natural(前段で置換済み)。
		// check_symmetry は厳密(tol = 0): chol カーネルは下三角のみを読む
		// [確認済み tsparse_sparse_chol.hpp 冒頭契約] ため、非対称入力を黙って
		// 通すと「下三角の対称化 = 別の行列が SPD」への主張のすり替えが起きる。
		// pd_tol は既定 0(certified d > 0 のみ受理)。
		chol_options<_T> copt;
		copt.method = vcp::sparse_chol_method::auto_select;
		copt.ordering = vcp::sparse_chol_ordering::natural;
		copt.check_symmetry = true;
		copt.symmetry_tol = _T(0);

		// --- A_s = Awork − shift·I の外側(下端)構成【対応: isspd.m L212-214 の
		// setround(-1)】。対角 A_jj − shift を kv 区間で計算し .lower() を採る =
		// A_s ⪯ Awork − shift·I を成分ごとに保証する向き。差が厳密 0 の対角は
		// 格納しない(spmats invariant。対角欠落は chol の certified 0 ピボット
		// = not_positive_definite で正しく false 側に倒れる)。
		auto build_shifted = [&](const _T& shift, PT_& As) {
			const PT_ Ac = Awork.as_csr();
			const std::vector<index_type>& outerv = Ac.outer_index();
			const std::vector<index_type>& innerv = Ac.inner_index();
			const std::vector<_T>& valv = Ac.values();
			const IT_ ish(shift);
			As.clear();
			As.resize(n, n);
			for (index_type i = 0; i < n; i++) {
				for (index_type p = outerv[static_cast<std::size_t>(i)];
				     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
					const index_type j = innerv[static_cast<std::size_t>(p)];
					if (i == j) {
						const _T v = (IT_(valv[static_cast<std::size_t>(p)]) - ish).lower();
						if (!(v == _T(0))) As.add(i, j, v);
					}
					else {
						As.add(i, j, valv[static_cast<std::size_t>(p)]);
					}
				}
			}
			As.finalize();
		};

		// --- 4. 近似 Cholesky chol(Awork) = L【対応: isspd.m L196-197】。
		// 失敗なら false【L224-226】(近似段だが、失敗 = 以降の近似が作れない)。
		PT_ Lf;
		std::vector<index_type> cperm;
		{
			const chol_result<_T, index_type> cr0 =
				Awork.policy_chol_with_info(Lf, cperm, copt);
			if (cr0.status != vcp::sparse_chol_status::success) return false;
		}

		// --- 5. 逆反復 t ≈ λmin(Awork)【対応: isspd.m L199-209】---
		// この段は全て点計算の近似で保証に使わない(意図コメント必須): シフト
		// の値の良し悪しは成功率のみに影響し、主張は 7.〜9. の「A_s の chol
		// 完走 + certified 判定」だけが担う。MATLAB の C は上三角(A = CᵀC)、
		// VCP は下三角 L(A = L Lᵀ)なので C = Lᵀ: C\v = Lᵀ\v(後退代入)、
		// C'\v = L\v(前進代入)。既存の chol 因子 solve は無いため区画内の
		// 代入ループで実装。
		const PT_ Lc = Lf.as_csc();
		const std::vector<index_type>& lptr = Lc.outer_index();
		const std::vector<index_type>& lind = Lc.inner_index();
		const std::vector<_T>& lval = Lc.values();
		// L の対角が各列先頭に格納されていることの検証(chol 契約: 正対角の
		// 非単位下三角。CSC 列内は行昇順)。欠落は防御的 false。
		for (index_type j = 0; j < n; j++) {
			if (lptr[static_cast<std::size_t>(j)] >= lptr[static_cast<std::size_t>(j) + 1]) return false;
			if (lind[static_cast<std::size_t>(lptr[static_cast<std::size_t>(j)])] != j) return false;
		}
		// 前進代入 b ← L⁻¹ b(列指向)【対応: isspd.m の C'\ 】
		auto fwd = [&](std::vector<_T>& b) {
			for (index_type j = 0; j < n; j++) {
				const index_type p0 = lptr[static_cast<std::size_t>(j)];
				b[static_cast<std::size_t>(j)] /= lval[static_cast<std::size_t>(p0)];
				for (index_type p = p0 + 1; p < lptr[static_cast<std::size_t>(j) + 1]; p++) {
					b[static_cast<std::size_t>(lind[static_cast<std::size_t>(p)])] -=
						lval[static_cast<std::size_t>(p)] * b[static_cast<std::size_t>(j)];
				}
			}
		};
		// 後退代入 b ← L⁻ᵀ b(列 j と b の内積)【対応: isspd.m の C\ 】
		auto bwd = [&](std::vector<_T>& b) {
			for (index_type j = n - 1; j >= 0; j--) {
				const index_type p0 = lptr[static_cast<std::size_t>(j)];
				_T acc = b[static_cast<std::size_t>(j)];
				for (index_type p = p0 + 1; p < lptr[static_cast<std::size_t>(j) + 1]; p++) {
					acc -= lval[static_cast<std::size_t>(p)]
						* b[static_cast<std::size_t>(lind[static_cast<std::size_t>(p)])];
				}
				b[static_cast<std::size_t>(j)] = acc / lval[static_cast<std::size_t>(p0)];
				if (j == 0) break;   // index_type が unsigned でも安全な脱出
			}
		};
		// 開始ベクトル【対応: isspd.m L199: y = C\(sum(abs(C),1))'】
		// sum(abs(C),1)' の第 i 成分 = |C| の列 i 和 = |L| の行 i 和(C = Lᵀ)。
		std::vector<_T> y(static_cast<std::size_t>(n), _T(0));
		for (index_type j = 0; j < n; j++) {
			for (index_type p = lptr[static_cast<std::size_t>(j)];
			     p < lptr[static_cast<std::size_t>(j) + 1]; p++) {
				y[static_cast<std::size_t>(lind[static_cast<std::size_t>(p)])] +=
					abs(lval[static_cast<std::size_t>(p)]);
			}
		}
		bwd(y);
		// t = n/(yᵀy)【対応: isspd.m L200】。退化(yᵀy が正でない = 0/NaN)は
		// 無主張 false(失敗側ゲート)。n の double 経由は近似段(厳密性不要)。
		_T yy = _T(0);
		for (index_type i = 0; i < n; i++) {
			yy += y[static_cast<std::size_t>(i)] * y[static_cast<std::size_t>(i)];
		}
		if (!(yy > _T(0))) return false;
		_T t = _T(static_cast<double>(n)) / yy;            // approximation only
		// 逆反復【対応: isspd.m L201-209: 3 回以上 + 相対変化 1e-4】。反復上限
		// 50 は区画内定数(isspd.m に無い発散・NaN 停留防御。近似段なので
		// 打ち切っても主張には影響しない)。停止判定は _T の abs で評価。
		{
			int it = 0;
			bool notready = true;
			std::vector<_T> z(static_cast<std::size_t>(n));
			while ((notready || it < 3) && it < 50) {
				it++;
				const _T told = t;
				z = y;
				fwd(z);                                    // z = (C'\y)【対応: isspd.m L204】
				for (index_type i = 0; i < n; i++) {
					z[static_cast<std::size_t>(i)] *= t;   // ·s【同 L204】
				}
				y = z;
				bwd(y);                                    // y = C\z【対応: isspd.m L205】
				_T zz = _T(0), yy2 = _T(0);
				for (index_type i = 0; i < n; i++) {
					zz  += z[static_cast<std::size_t>(i)] * z[static_cast<std::size_t>(i)];
					yy2 += y[static_cast<std::size_t>(i)] * y[static_cast<std::size_t>(i)];
				}
				if (!(yy2 > _T(0))) return false;          // 退化は無主張
				t = zz / yy2;                              // 【対応: isspd.m L206】
				notready = (abs(t - told) > _T(1e-4) * abs(t + told));   // 【同 L207-208】
			}
		}
		// Rump 1993 Alg. 2.2 の「σ > 半径」【対応: isspd.m L210 / L221-222:
		// s−M not positive】: 成立しなければ無主張 false(NaN も失敗側)。
		if (!(t > Mrho)) return false;

		// --- 6. シフト s = M_ρ + 0.8·(t − M_ρ)(D-3。点版は M_ρ = 0 で 0.8 t =
		// isspd.m L211 と同値)。係数 0.8 は自由パラメータ(成功率のみ)。
		_T s = Mrho + _T(0.8) * (t - Mrho);
		last_verify_info.stage = 2;                        // ノルム/シフト評価段
		last_verify_info.spd_shift = static_cast<double>(s);   // 記録専用の double 変換(設計 §3)

		// --- 7.〜10. A_s 構成 → chol(A_s) = L₂ → E ⊇ L₂L₂ᵀ − A_s → 判定。
		// 不成立なら s ← M_ρ + (s − M_ρ)/2 で再試行(D-8)。
		static const int kSpiAccurateMaxRetry = 3;   // success-rate knob only; no claim depends on it
		//   (Rump 1993, Alg. 2.2: "If algorithm 2.1 does not finish successfully try smaller values")
		for (int attempt = 0; attempt <= kSpiAccurateMaxRetry; attempt++) {
			if (attempt > 0) s = Mrho + (s - Mrho) / _T(2);
			last_verify_info.spd_shift = static_cast<double>(s);   // 記録専用(設計 §3)
			// 7. A_s = Awork − s·I(対角下端)【対応: isspd.m L212-214】。
			//    よって Awork = A_s + s·I + N、N ⪰ 0 対角。
			PT_ As;
			build_shifted(s, As);
			// 8. chol(A_s) = L₂【対応: isspd.m L215-216】。失敗なら再試行【L218-220】
			PT_ L2;
			std::vector<index_type> perm2;
			const chol_result<_T, index_type> cr2 =
				As.policy_chol_with_info(L2, perm2, copt);
			if (cr2.status != vcp::sparse_chol_status::success) continue;
			// 9. 残差 E ⊇ L₂·L₂ᵀ − A_s【対応: isspd.m L217: norm(C'*intval(C)−A)】
			//    を片側昇格の区間累積ループで評価(R-11)。点×点の積は最近接
			//    丸めで方向保証が残らないため、積の【片側】を区間へ昇格した
			//    acc += interval(l_ik)·l_jk(kv acceptable_n の 2 乗算経路)が
			//    誤差包含の理論最小。減算は A_s を点のまま(interval − 点)。
			//    対称性: (i,j) と (j,i) は同一の k 集合を同じ昇順で畳むため E は
			//    端点まで厳密対称 --- ‖E‖₂ ≤ ‖E‖∞ の前提が保たれる。
			base_type E;
			{
				const PT_ L2r = L2.as_csr();
				const PT_ L2cc = L2.as_csc();
				const PT_ Asr = As.as_csr();
				const std::vector<index_type>& ro = L2r.outer_index();
				const std::vector<index_type>& ri = L2r.inner_index();
				const std::vector<_T>& rv = L2r.values();
				const std::vector<index_type>& co = L2cc.outer_index();
				const std::vector<index_type>& ci = L2cc.inner_index();
				const std::vector<_T>& cv = L2cc.values();
				const std::vector<index_type>& ao = Asr.outer_index();
				const std::vector<index_type>& ai = Asr.inner_index();
				const std::vector<_T>& av = Asr.values();
				E.resize(n, n);
				// 行ごとの疎アキュムレータ(dense 作業配列 + touch リスト。
				// メモリ O(n)、行間で再利用)
				std::vector<interval_type> acc(static_cast<std::size_t>(n), interval_type(_T(0)));
				std::vector<char> touched(static_cast<std::size_t>(n), 0);
				std::vector<index_type> tlist;
				for (index_type i = 0; i < n; i++) {
					tlist.clear();
					// (L₂·L₂ᵀ)_ij = Σ_k l_ik·l_jk: 行 i の各 k について列 k
					// (= L₂ᵀ の行 k)を走査して行 i に散布
					for (index_type p = ro[static_cast<std::size_t>(i)];
					     p < ro[static_cast<std::size_t>(i) + 1]; p++) {
						const index_type kk = ri[static_cast<std::size_t>(p)];
						const _T vik = rv[static_cast<std::size_t>(p)];
						for (index_type q = co[static_cast<std::size_t>(kk)];
						     q < co[static_cast<std::size_t>(kk) + 1]; q++) {
							const index_type j = ci[static_cast<std::size_t>(q)];
							const std::size_t uj = static_cast<std::size_t>(j);
							if (!touched[uj]) {
								touched[uj] = 1;
								acc[uj] = interval_type(_T(0));
								tlist.push_back(j);
							}
							// 【必要な昇格(SR-1 対象外 --- 区別コメント)】片側のみ区間化
							acc[uj] += interval_type(vik) * cv[static_cast<std::size_t>(q)];
						}
					}
					// − A_s 行 i(点のままの混合減算 --- interval − 点)
					for (index_type p = ao[static_cast<std::size_t>(i)];
					     p < ao[static_cast<std::size_t>(i) + 1]; p++) {
						const index_type j = ai[static_cast<std::size_t>(p)];
						const std::size_t uj = static_cast<std::size_t>(j);
						if (!touched[uj]) {
							touched[uj] = 1;
							acc[uj] = interval_type(_T(0));
							tlist.push_back(j);
						}
						acc[uj] -= av[static_cast<std::size_t>(p)];
					}
					// 発行(厳密 [0,0] は格納しない --- spmats invariant)
					for (std::size_t tt = 0; tt < tlist.size(); tt++) {
						const std::size_t uj = static_cast<std::size_t>(tlist[tt]);
						if (!(acc[uj].lower() == _T(0) && acc[uj].upper() == _T(0))) {
							E.add(i, tlist[tt], acc[uj]);
						}
						touched[uj] = 0;
					}
				}
				E.finalize();
			}
			// ‖E‖∞ の区間評価(SPFN norminf。normtwo は保証に使用禁止 ---
			// 対称なら ‖E‖₂ ≤ √(‖E‖₁‖E‖∞) = ‖E‖∞)
			base_type NB;
			E.norminf(NB);
			const interval_type nrm = NB.get(0, 0);
			last_verify_info.stage = 3;                    // 判定段
			// 判定(D-3): 証明された λmin(M) の下界 s − ‖E‖∞ が半径上界 M_ρ を
			// 超えること。certified 比較: 右辺の上端が NaN / +∞ なら不成立
			// (無主張側)。点版は M_ρ = 0 で isspd.m L217 の s > ‖E‖ と同値。
			if (s > (nrm + IT_(Mrho)).upper()) {
				last_verify_info.stage = 4;                // 成功
				return true;
			}
		}
		return false;                                      // 再試行を尽くした = 無主張
	}


	// >>> END [SPI-U2c] <<<
private:
	// is_spd 外殻の器初期化: is_spd 所有フィールドのみ(設計 SPI-G1 §2.5)。
	void spi_clear_spd_fields_() const {
		// is_spd-owned fields only. veigs counters (veigs_bisection_steps / veigs_clusters)
		// and ap_* snapshots are NOT touched: is_spd is called from inside veigs (design §2.5).
		last_verify_info.spd_shift = 0.0;
		last_verify_info.spd_rho_rad = 0.0;
		last_verify_info.spd_scale_exp.clear();
		last_verify_info.spd_collatz_iters = 0;
		last_verify_info.spd_accurate_path = false;
		last_verify_info.spd_scale_abandoned = false;
		last_verify_info.stage = 0;
	}

public:

	// -------------------------------------------------------------------
	// policy_is_symmetric の名前隠蔽(SPI-R4、設計 §8)
	//
	// 基底はテンプレートメンバ [確認済み spmats.hpp L1216-1254(スナップ
	// ショット時 L865- から移動)] のため virtual 不可。名前隠蔽で差し替える
	// (imats が mats の非 virtual linearsolve に用いるのと同一機構)。
	// 名前隠蔽の常として、基底参照経由の呼び出しでは基底版(tol 既定
	// 1e-12 の近似判定)が動く。通常の spmatrix 利用では起きないが、
	// docs/usage に明記する。
	// -------------------------------------------------------------------

	// >>> REVIEW-REQUIRED [SPI-R4: policy_is_symmetric 厳密版(名前隠蔽)] <<<
	// STATUS: REVIEWED-OK (Kouta Sekine, 2026-07-26)
	// CLAIM: tol なし版は「全 (i,j) で A_ij と A_ji が区間として一致
	//   (端点スカラー一致)」のとき、かつそのときに限り true を返す。
	//   非格納要素は厳密 [0,0] として比較され、格納値は spmats invariant
	//   (finalize 済み格納列に厳密 [0,0] なし)により、片側のみ格納の要素は
	//   必ず false 側に倒れる。NaN 端点は端点 == が false になり(double /
	//   kv::dd で実測 [Phase 0 M-4])false 側 = 安全側。kv::mpfr<N> の
	//   スカラー == も現行 kv は mpfr_equal_p 系 [kv/mpfr.hpp L306 / L350] で
	//   NaN==NaN・NaN==0 とも偽(SPI-V3 kv_recheck.cpp で再確認 2026-10-14。
	//   旧記述「mpfr は NaN==NaN が真」は陳腐化 --- [SPI-F3] 訂正)。
	//   ±0 は == で同値(true で正しい: 値として同一)。
	// REDUCES-TO: 端点スカラーの == 比較のみ(区間演算・丸めなし)。走査は
	//   基底 policy_is_symmetric の CSR + 鏡像二分探索構造を踏襲し、比較行の
	//   みを置換 [確認済み spmats.hpp L1225-1253]。
	// SELF-ARITHMETIC: なし
	// TESTS: sandbox/tests/spimats_policy_test.cpp
	// -------------------------------------------------------------------

	// tol なし版 = 厳密対称判定(基底の「tol 既定 1e-12 へ委譲」を置換)。
	// Matrix の要件は基底と同じ(as_csr / outer_index / inner_index / values
	// を持ち、値型が kv::interval<_T>、index が index_type であること)。
	template <class Matrix>
	bool policy_is_symmetric(const Matrix& A) const {
		if (A.rowsize() != A.columnsize()) return false;
		const base_type Acsr = A.as_csr();   // [SPI-F3 F3-11] as_csr() は基底型を返す(Matrix = spimats でもコンパイル可)
		const std::vector<index_type>& outerv = Acsr.outer_index();
		const std::vector<index_type>& innerv = Acsr.inner_index();
		const std::vector<interval_type>& val = Acsr.values();
		for (index_type i = 0; i < Acsr.rowsize(); i++) {
			for (index_type p = outerv[static_cast<std::size_t>(i)];
			     p < outerv[static_cast<std::size_t>(i + 1)]; p++) {
				const index_type j = innerv[static_cast<std::size_t>(p)];
				if (i == j) continue;
				const index_type first = outerv[static_cast<std::size_t>(j)];
				const index_type last  = outerv[static_cast<std::size_t>(j + 1)];
				const typename std::vector<index_type>::const_iterator begin = innerv.begin() + first;
				const typename std::vector<index_type>::const_iterator end   = innerv.begin() + last;
				typename std::vector<index_type>::const_iterator it = std::lower_bound(begin, end, i);
				interval_type mirrored = interval_type(_T(0));   // 非格納 = 厳密 [0,0](基底と同じ値意味論)
				if (it != end && *it == i)
					mirrored = val[static_cast<std::size_t>(it - innerv.begin())];
				const interval_type& a = val[static_cast<std::size_t>(p)];
				// 端点スカラー一致(kv の区間 == は certainly-equal であり
				// 端点一致ではない --- 設計 §8 の比較規約)
				if (!(a.lower() == mirrored.lower() && a.upper() == mirrored.upper()))
					return false;
			}
		}
		return true;
	}

	// tol 付き版 = 封鎖(呼び出しをコンパイル時に拒否)。
	template <class Matrix, typename _Tol>
	bool policy_is_symmetric(const Matrix& A, const _Tol& tol) const {
		static_assert(spimats_detail::always_false_<Matrix>::value,
			"vcp::spimats::policy_is_symmetric(tol): a tolerance-based symmetry check has no verified meaning on an interval matrix; extract the midpoint into the point policy (mid_for_approximation) and check there, or use the exact (no-tol) overload");
		(void)A; (void)tol;
		return false;   // not reached
	}

	// >>> END [SPI-R4] <<<

protected:
	// -------------------------------------------------------------------
	// 混合演算メンバ(protected virtual、imats 同形。設計 §6)
	// カーネル実体は vcp/spimats_base/spimats_mixed_ops.hpp(SPI-R2)。
	// 本メンバは転送のみ。_AP はクラステンプレート引数で固定されるため
	// 通常メンバとして virtual にでき、将来の派生ポリシーが BLAS / 並列版に
	// 差し替える余地を残す。ユーザーコード側の演算子は作らない(確定 ---
	// 基底に異値型 policy_mul がなく、演算子化は既存ファイル変更を要する)。
	// 将来 mul_iv_m 等を追加する場合は、spimats_mixed_ops.hpp に新しい
	// REVIEW-REQUIRED 区画(SPI-R6 以降)としてカーネルを置き、ここに転送
	// メンバを追記する。
	// -------------------------------------------------------------------

	// IC = (*this) * B(区間疎 × 点疎)
	virtual void mul_im_m(const _AP& B, base_type& IC) const {
		vcp::spimats_kernel::mul_im_m(static_cast<const base_type&>(*this), B, IC);
	}
	// IC = B * (*this)(点疎 × 区間疎)
	virtual void mul_m_im(const _AP& B, base_type& IC) const {
		vcp::spimats_kernel::mul_m_im(B, static_cast<const base_type&>(*this), IC);
	}
	// IC = (*this) + B(区間疎 + 点疎。差は点側を −1 倍して加算)
	virtual void add_im_m(const _AP& B, base_type& IC) const {
		vcp::spimats_kernel::add_im_m(static_cast<const base_type&>(*this), B, IC);
	}
	// IC = B + (*this)(点疎 + 区間疎)
	virtual void add_m_im(const _AP& B, base_type& IC) const {
		vcp::spimats_kernel::add_m_im(B, static_cast<const base_type&>(*this), IC);
	}
	// iy = (*this) * x(区間疎行列 × 点ベクトル: 残差 A x~)
	virtual void mul_im_v(const std::vector<_T>& x, std::vector<interval_type>& iy) const {
		vcp::spimats_kernel::mul_im_v(static_cast<const base_type&>(*this), x, iy);
	}
	// iy = B * ix(点疎行列 × 区間ベクトル: R の適用。*this は関与しない)
	virtual void mul_m_iv(const _AP& B, const std::vector<interval_type>& ix,
	                      std::vector<interval_type>& iy) const {
		vcp::spimats_kernel::mul_m_iv(B, ix, iy);
	}
	// iy^T = x^T * (*this)(点横ベクトル × 区間疎行列: RA − I の行ごと構成)
	virtual void mul_v_im(const std::vector<_T>& x, std::vector<interval_type>& iy) const {
		vcp::spimats_kernel::mul_v_im(x, static_cast<const base_type&>(*this), iy);
	}

	// -------------------------------------------------------------------
	// options の点型抽出(設計 §5.1)
	// 利用規約: スカラーノブ(tol 等)は退化区間(点区間)で指定する。
	// 点型 _T への抽出は一律 upper()(退化区間なら mid / lower と同値。
	// 幅がある区間を渡されても tol として緩い側 = 安全側に倒れ、misuse
	// 検査は不要)。整数・enum・bool フィールドはそのまま写す。
	// 注意: フィールド一覧は HEAD(2026-07-18)時点の全数写像。基底
	// options にフィールドが追加された場合はここへの追記が必要(追記漏れは
	// 既定値で動く = ノブが黙って無視される)。
	// -------------------------------------------------------------------

	sparse_lu_options<_T> spi_extract_sparse_lu_options_(
		const sparse_lu_options<kv::interval<_T> >& o) const {
		sparse_lu_options<_T> p;
		p.method = o.method;
		p.ordering = o.ordering;
		p.pivoting = o.pivoting;
		p.pivot_threshold = o.pivot_threshold.upper();
		p.absolute_pivot_tolerance = o.absolute_pivot_tolerance.upper();
		p.zero_tolerance = o.zero_tolerance.upper();
		p.near_zero_tolerance = o.near_zero_tolerance.upper();
		p.equilibration = o.equilibration;
		p.iterative_refinement = o.iterative_refinement;
		p.blas_min_block_size = o.blas_min_block_size;
		p.supernode_relaxation = o.supernode_relaxation;
		p.panel_size = o.panel_size;
		p.growth_threshold = o.growth_threshold.upper();
		p.static_pivot_perturbation = o.static_pivot_perturbation.upper();
		p.compute_condition_estimate = o.compute_condition_estimate;
		p.iterative_refinement_max_iterations = o.iterative_refinement_max_iterations;
		p.iterative_refinement_tolerance = o.iterative_refinement_tolerance.upper();
		p.supernodal_numeric_diagnostic_leftlooking = o.supernodal_numeric_diagnostic_leftlooking;
		p.supernodal_self_symbolic = o.supernodal_self_symbolic;
		p.supernodal_inplace_frontal = o.supernodal_inplace_frontal;
		p.supernodal_native_check_residual = o.supernodal_native_check_residual;
		p.supernode_panel_maxsup = o.supernode_panel_maxsup;
		return p;
	}

	linear_solve_options<_T> spi_extract_lss_options_(
		const linear_solve_options<kv::interval<_T> >& o) const {
		linear_solve_options<_T> p;
		p.method = o.method;
		p.max_iter = o.max_iter;
		p.tol = o.tol.upper();
		p.check_symmetric = o.check_symmetric;
		p.use_relative_residual = o.use_relative_residual;
		p.restart = o.restart;
		p.preconditioner = o.preconditioner;
		p.sparse_lu = spi_extract_sparse_lu_options_(o.sparse_lu);
		return p;
	}

	// -------------------------------------------------------------------
	// >>> REVIEW-REQUIRED [SPI-R6: B=I 厳密検出] <<<
	// STATUS: REVIEWED-OK (Kouta Sekine, 2026-07-26)
	// CLAIM: true ⟺ B は厳密単位区間行列(正方、格納 nnz = n、各行
	//   ちょうど 1 要素、位置 (i,i)、値の端点が厳密に (1,1))。
	//   false は無主張(一般化経路に流れるだけで、保証コア SPI-U5 が
	//   見る問題は不変 = 検出漏れは無害)。true は近似層の経路切替
	//   (標準 eigs へ)の根拠となる厳密性の主張であり保証の連鎖に乗る。
	//   true 誤判定時も U5 の検証が誤った近似解を弾くため誤った包含は
	//   出ない(健全性は U5 が守る)が、主張である以上レビュー対象。
	// REDUCES-TO: 構造走査(CSR)+ 端点スカラー == 比較のみ
	//   (区間演算・丸めなし)
	// SELF-ARITHMETIC: なし
	// TESTS: sandbox/tests/spimats_policy_test.cpp((g) SPI-R6 ケース)
	// -------------------------------------------------------------------
	// B が厳密単位行列(正方・パターン = 対角のみ・全対角 [1,1])かの厳密
	// 判定(端点スカラー比較のみ。設計 §8 の比較規約に適合)。
	// 用途: override #3 の【AP 呼び出しの経路選択のみ】。基底 AP の
	// 一般化固有値ソルバーは【A が非対角のとき、既定 structure=auto では】
	// B が単位行列でも unsupported_generalized_non_diagonal を返す
	// [実測 2026-07-18: 対角 A なら既定でも収束、symmetric ヒントを渡せば
	// 非対角 A + B=I の一般化も収束する。失敗の駆動因は
	// 「非対角 A + 既定 auto」であり B=I そのものではない]。
	// 標準 policy_eigs_with_info は同条件(非対角 A・既定オプション)で
	// 収束する [実測同日] ため、#2 の B=I 転送で来た標準固有値問題 ---
	// 利用者は structure ヒントを渡さないのが普通 --- を近似層では標準
	// policy_eigs_with_info に戻す。保証コア(SPI-U5)が見る形は一般化
	// (B 引数付き)のままで、意味論の変更はない。
	bool spi_is_exact_identity_(const spmats<kv::interval<_T>, index_type>& B) const {
		if (B.rowsize() != B.columnsize()) return false;
		if (!B.is_finalized()) B.finalize();
		const spmats<kv::interval<_T>, index_type> C = B.as_csr();
		if (C.stored_nnz() != C.rowsize()) return false;
		const std::vector<index_type>& outer = C.outer_index();
		const std::vector<index_type>& inner = C.inner_index();
		const std::vector<kv::interval<_T> >& val = C.values();
		for (index_type i = 0; i < C.rowsize(); i++) {
			const index_type p = outer[static_cast<std::size_t>(i)];
			if (outer[static_cast<std::size_t>(i) + 1] - p != 1) return false;
			if (inner[static_cast<std::size_t>(p)] != i) return false;
			const kv::interval<_T>& x = val[static_cast<std::size_t>(p)];
			if (!(x.lower() == _T(1) && x.upper() == _T(1))) return false;
		}
		return true;
	}
	// >>> END [SPI-R6] <<<

	eig_options<_T> spi_extract_eig_options_(
		const eig_options<kv::interval<_T> >& o) const {
		eig_options<_T> p;
		p.method = o.method;
		p.structure = o.structure;
		p.target = o.target;
		p.tol = o.tol.upper();
		p.shift = o.shift.upper();
		p.use_shift = o.use_shift;
		p.max_iter = o.max_iter;
		p.subspace_dim = o.subspace_dim;
		p.allow_dense_conversion = o.allow_dense_conversion;
		p.max_dense_size = o.max_dense_size;
		p.orthogonalization = o.orthogonalization;
		p.random_seed = o.random_seed;
		p.random_start = o.random_start;
		p.compute_residual_history = o.compute_residual_history;
		p.shift_invert_solver = o.shift_invert_solver;
		p.shift_invert_lu = spi_extract_sparse_lu_options_(o.shift_invert_lu);
		p.allow_complex_pairs = o.allow_complex_pairs;
		return p;
	}

	// -------------------------------------------------------------------
	// override #1(設計 §3): 連立一次 mid → AP 近似解 → 保証(SPI-U3)。
	// 非 virtual 外殻 policy_lss_with_info(finalize 保証)から到達する。
	// strict 系(policy_lss)は override しない: 返る result は常に
	// converged == true なので基底の !converged → throw チェック
	// [確認済み spmats_lss.hpp L698-706]は発火しない(設計 §3)。
	// -------------------------------------------------------------------
	virtual linear_solve_result<kv::interval<_T> > policy_lss_with_info_impl(
		const std::vector<kv::interval<_T> >& b,
		const linear_solve_options<kv::interval<_T> >& opt) const override
	{
		// 入口検査(misuse -> vcp::error 系。spumar override の鏡映)
		if (this->rowsize() != this->columnsize()) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_lss_with_info: matrix must be square");
		}
		if (b.size() != static_cast<std::size_t>(this->rowsize())) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_lss_with_info: rhs dimension mismatch");
		}
		if (opt.max_iter == 0) {
			vcp::throw_error<vcp::invalid_argument>("spimats::policy_lss_with_info: max_iter must be positive");
		}
		if (!(opt.tol.upper() > _T(0))) {   // 失敗側ゲート(certified-only の書き方)
			vcp::throw_error<vcp::invalid_argument>("spimats::policy_lss_with_info: tol must be positive");
		}
		last_verify_info.clear();

		// mid 抽出(SPI-R1 mid_for_approximation: 近似入力専用・保証への
		// 使用禁止)→ LU 直解で近似解(DD-10)
		_AP Amid;
		vcp::mid_for_approximation(static_cast<const base_type&>(*this), Amid);
		std::vector<_T> bmid;
		vcp::mid_for_approximation(b, bmid);

		// >>> DD-10 [SPI-D1-R1 R-1: lss 前段 LU 1 回化差し替え] <<<
		// 旧経路(v1.1): spi_extract_lss_options_ + Amid.policy_lss_with_info
		// で近似解を得たのち、U3 内で R 行用に別途 LU 分解していた
		// (= 分解 2 回)。本差し替え(オーナー裁定 2026-07-24、DD-10)で
		// LUX 抽出 policy_lu_with_info を【全体で 1 回】だけ実行し、同じ
		// 因子 (Lf, Uf, pperm, qperm) で
		//   - xh: policy_lu_solve_with_info(順方向解)
		//   - U3 の Ri: policy_lu_inverse_row_with_info(転置解)
		// の両方を賄う。R 分解のノブは既定値(equilibration=false --- SSC 形
		// はスケーリング非対応 [確認済み spmats.hpp L1507-1509]。利用者の
		// opt.sparse_lu / method ノブは本経路では近似解の作り方に影響しない)。
		// xh 質注記: 旧経路の IR(反復改良)を経由しないため xh は粗いことが
		// ある --- xh の質は包含幅にのみ影響し、保証(後段の全量区間評価)には
		// 影響しない(mid・xh は任意の点で可 --- 設計 §6 前段)。
		spmats<_T, index_type> Lf, Uf;
		std::vector<index_type> pperm, qperm;
		{
			const lu_extract_options<_T> lopt;   // 既定(equilibration=false)
			const lu_extract_result<_T, index_type> lres =
				Amid.policy_lu_with_info(Lf, Uf, pperm, qperm, lopt);
			if (!(lres.status == vcp::sparse_lu_extract_status::success)) {
				// 分解が立たない = 近似層の死(xh も R も作れない)。
				// 設計 §4.1 の state_error(mats が Lapack 失敗で死ぬのと同型)
				vcp::throw_error<vcp::state_error>(
					"spimats::policy_lss_with_info: LU factorization of mid(A) failed (see last_verify_info)");
			}
		}
		std::vector<_T> xh0;
		{
			const lu_apply_result sres = Amid.policy_lu_solve_with_info(
				Lf, Uf, pperm, qperm, bmid, xh0);
			if (!(sres.status == vcp::lu_apply_status::success)) {
				vcp::throw_error<vcp::state_error>(
					"spimats::policy_lss_with_info: LU forward solve for xh failed (see last_verify_info)");
			}
		}
		// スナップショット縮退(R-1): 旧経路の linear_solve_result 丸ごと
		// 保持は LU 直解には存在しないため、ap_lss_result.x = xh のみを
		// 記録する(converged / 残差診断フィールドは既定値のまま ---
		// 「x が LU 直解として得られた」ことは本コメントと U3 到達が証す。
		// フライトレコーダは計算の記録であり主張ではない --- G3)。
		last_verify_info.ap_lss_result.x = xh0;
		// >>> END DD-10 [SPI-D1-R1 R-1] <<<

		// 成功時 result の器(設計 §5.3 有効性表):
		//   x            = 保証付き区間(唯一の保証値。SPI-U3 が埋める)
		//   converged    = true(AP 失敗・保証失敗はここに到達する前に throw)
		//   iterations / residual_norm_* = 既定値のまま(近似層の値は写さない、G4)
		//   method       = 要求値を保持(spumar S-3 と同旨。linear_solve_result に
		//                  文字列識別子フィールドはないため enum のまま)
		linear_solve_result<kv::interval<_T> > result;
		result.converged = true;
		result.method = opt.method;

		// >>> USER-IMPLEMENTS [SPI-U3: verified lss コア] <<<
	// [SPI-P1] 実装は SPI-D1 オラクル(spimats_debug.hpp、SHA f54f8894…)
	// 由来。オーナー逐行検証済み(2026-07-24)。USER-IMPLEMENTS の正しさの
	// 担保者はオーナー(SPI-1 二重マーカー規律のとおり)。
		//     正しさはオーナーが実装・保証。以下の想定コードコメントは
		//     補助であり拘束しない。定理由来箇所は UNVERIFIED マーク付き。
		//
		// 【想定フロー】imats::linearsolve [確認済み vcp/imats.hpp L293-364]
		//   の疎版。x~ = ap.x(点近似解)として:
		//   (1) R ~ inv(mid(A)) を点密で構成(imats は _P::inv() を使用
		//       [確認済み L302]。疎の逆は一般に密なので、R は _DP /
		//       matrix<kv::interval<_T>, vcp::imats<_T, _DP>> の密の枠で持つ。
		//       n² メモリが乗る --- 大規模疎向けの R の疎近似(AINV/FSAI 等)
		//       への差し替えはオーナー判断)。
		//       [v1.5 注記: _DP は区間密ポリシーに改定。本想定は DD-7 で失効済み]
		//       【SPI-D1 注記】上の「R を点密で構成」はオーナー裁定
		//       2026-07-23(SPI-D1 DD-7)で【行ストリーミング】に上書き:
		//       R は実体化せず、LU 分解 1 回から 1 行ずつ取り出す(設計
		//       SPI-D1_design.md §6)。想定コメント自体は履歴として残す。
		//   (2) 残差の区間評価: mul_im_v(x~)で [A]x~ を包含し b を引く
		//       → R を掛けて |R([A]x~ − b)| を得る(imats L307-316 と同型)。
		//   (3) [RA − I] の構成: R の第 i 行(点横ベクトル)ごとに
		//       mul_v_im で R_i・[A] を包含し、対角要素から 1 を引く
		//       (I は厳密。行ごと構成で密 n×n の区間行列を確保)。
		//       【SPI-D1 注記】上の「密 n×n の区間行列を確保」はオーナー裁定
		//       2026-07-23(SPI-D1 DD-7)で上書き: G の行は保持せず、行ごとに
		//       abs 行和 RiAe[i] へ畳み込む(メモリ O(n)/行。設計 §6.2)。
		//   (4) G = ‖RA − I‖_inf を基底 norminf(SPFN)で評価し
		//       G.get(0,0).upper() < 1 を判定 [確認済み: 区間 norminf の連鎖は
		//       kv abs → 区間 += → ADL kv max(端点ごと max、L658)で正しい
		//       包含。normtwo は policy_eigs 近似 + sqrt であり保証に使用禁止]。
		//   (5) NaN 防御(設計 §10): 判定の前に norminf 結果の端点有限性
		//       チェックを置く。Phase 0 M-3 実測では kv 区間演算は有限入力から
		//       NaN 端点を生成しなかった(±∞ 端点は生成・伝播し、判定は安全側に
		//       失敗する)が、網羅実験ではないため本防御を維持する。kv の
		//       端点 max は NaN を順序依存で黙殺しうる [SPFN 実測]。
		//   (6) 誤差評価: ‖R([A]x~ − b)‖ / (1 − ‖RA − I‖) と |RA − I| による
		//       成分別評価の合成(imats L343-361 と同型の構成であることは
		//       確認済み。評価式自体の正当性は Rump 流の定理由来
		//       [UNVERIFIED: 一次資料未照合。オーナーが保証する])。
		//   (7) 十分条件不成立(‖RA − I‖ ≥ 1 等)は
		//       vcp::throw_error<vcp::verification_error>(...)(無主張)。
		//       成立時のみ result.x に x~ + 誤差区間を格納して return result。
		//   リカバリ(精度引き上げ・R の再計算等)を行うなら本区画の内側で。
		//   尽くしてなお失敗なら throw(ラッパー層は続行判断をしない、§4.1)。
		//
		// ===============================================================
		// [SPI-D1 v1.1 実装] U3: verified lss コア(系 3.2.1・行ストリーミング)
		// 規範出典: 『系 3.2.1(オーナー提供 2026-07-23)』+『オーナー提供
		//   U3 フロー(2026-07-23、訂正合意版)』+ 設計 SPI-D1_design.md §6。
		// 規範式(系 3.2.1): G := RA − I。‖G‖∞ < 1 ならば A は正則で、
		//   |x* − x̂| ≤ |R(b − Ax̂)| + ‖R(Ax̂ − b)‖∞ / (1 − ‖G‖∞) · (|G|e)
		//   (成分別)。区間 [A], [b] へは全量を区間演算で包含評価して
		//   持ち上げる [UNVERIFIED: 点定理の区間持ち上げの定式化は
		//   オーナー照合(設計 §6/§9)]。
		// 変数名・処理順序はオーナー提供フローに一致させる(DD-9)。
		// 【訂正合意の記録(2026-07-23)--- よくある取り違え 2 点】
		//   ①「abs 済み RiA を n 次元 tmp に累積」は列和(‖G‖₁)になるため
		//     【置かない】。行和は RiAe が二役(‖G‖∞ 材料 + 第 2 項ベクトル)
		//     で担う。
		//   ② 乗数は【Tnorm / (1 − Gnorm)】(分子 = 残差ノルム、分母 =
		//     1 − ‖G‖∞)。Gnorm / (1 − Tnorm) 型は誤り。
		// 丸め規律: 丸めモード操作なし。全合成は kv 区間演算。
		// ===============================================================
		{
			const std::size_t n = static_cast<std::size_t>(this->rowsize());
			last_verify_info.stage = 1;                    // 前処理段(LU)

			// --- 6.1 変数(名前はフローに一致 --- DD-9)---
			// Am / bm: [A], [b] の近似中心。骨組み前段が Amid / bmid として
			// 作成済み(設計 6.1「Amid(以下 Am)」)のため別名で束ねる。
			// mid は【近似で良い】: R・x̂ は任意の点で可であり、保証は後段の
			// 区間評価が全て担う(mid_for_approximation は近似入力専用 ---
			// R1 契約どおり保証には使わない)。
			const _AP& Am = Amid;
			const std::vector<_T>& bm = bmid;
			(void)bm;   // bm は前段 DD-10 の xh 生成(policy_lu_solve_with_info)
			            // で消費済み --- フローの「bm を L,U で解く」の実現形
			// xh: 近似解。前段 DD-10(R-1)が【この L,U から】
			// policy_lu_solve_with_info で解いた xh0 を用いる --- オーナー
			// フローの「L,U から xh を解く」がそのまま実現された形
			// (設計 6.1 xh 行)。
			const std::vector<_T>& xh = xh0;
			if (xh.size() != n) {
				// LU solve が寸法不整合の解を返した(近似層の死の一形態)
				vcp::throw_error<vcp::state_error>(
					"spimats::policy_lss_with_info: approximate solution has wrong dimension");
			}

			// L, U(+順列 p, q): 前段 DD-10(R-1)で抽出済みの LUX 因子
			// (P Am Q = L U、SSC 規約 [確認済み spmats.hpp L1499-1517]、
			// ピボット込み --- オーナー確認済み)を【共有】する。分解は
			// 全体で 1 回(旧 v1.1 の U3 内再分解は R-1 で削除)。以後は
			// 行ごとの転置三角解のみ(DD-7)。
			// (設計 6.1 の「ハンドル」でなく明示因子なのは、Ri の指定 API
			// policy_lu_inverse_row_with_info が明示因子 (L,U,p,q) を取る
			// [確認済み spmats.hpp L1546-1560] ため --- 型欄「_AP / ハンドル」
			// の前者)

			// r = [A]·xh − [b](フロー: r。R2 mul_im_v で [A]·xh を包含し
			// [b] を区間減算。xh は点・[A] [b] は区間)
			std::vector<interval_type> r(n);
			this->mul_im_v(xh, r);                         // [確認済み R2: 残差 A x~ 用]
			for (std::size_t i = 0; i < n; i++) {
				r[i] -= b[i];
			}

			// 器(フロー: Rir / RiAe / err。RiA は行ごとの作業横ベクトル)
			std::vector<interval_type> Rir(n);             // Rir[i] = Ri·r
			std::vector<interval_type> RiAe(n);            // abs 済み行和 = (|G|e)_i
			std::vector<_T> Ri;                            // R = Am⁻¹ の第 i 行
			std::vector<interval_type> RiA;                // Ri·[A](横ベクトル、行ごとに再利用)
			// 訂正合意 ①: ここに「列和累積の n 次元 tmp」は【置かない】。
			// abs 済み RiA を成分別に別ベクトルへ累積すると Σ_i |G_ij|
			// (列和 = ‖G‖₁ の材料)になり、系 3.2.1 の ‖G‖∞(行和の最大)
			// と取り違える。行和は次の行ループ内で RiAe[i] に畳む。

			// --- 6.2 行ループ(i = 0 .. n−1。R は行ストリーミング --- DD-7)---
			last_verify_info.stage = 2;                    // ノルム評価段
#if VCP_SPARSE_USE_OPENMP
			// SPOMP-2(設計 §3.2/§4.2): 行並列版。行 i の処理は行間で
			// データを共有せず(共有は const 読みのみ)、書き込み先
			// Rir[i]/RiAe[i] は行ごとに素。成功時の演算列は逐次側と行ごとに
			// 同一 = result.x は全ビット一致(D-8)。ループ後の縮約(6.3)は
			// 無変更。
			{
				// 準備済み因子の構築(1 回・逐次・並列領域外 --- D-9 の
				// 構築相。検証・CSC コピー・qinv 構築をここで完了し、行ループ
				// 内は const 読みの行解コアのみ --- 設計 §4.3)。構築失敗
				// (構造契約違反)は逐次側の「最初の行での inverse-row 失敗」と
				// 同じ verification_error 文言に合流。
				vcp::spmats_lu_extract_detail::prepared_lu_factors_<_T, index_type> Fprep;
				{
					const vcp::lu_apply_status ps =
						vcp::spmats_lu_extract_detail::prepare_lu_factors_<_T, index_type>(
							Lf, Uf, pperm, qperm, Fprep);
					if (!(ps == vcp::lu_apply_status::success)) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_lss_with_info: inverse-row extraction failed (no claim)");
					}
				}
				// mul_v_im 生配列オーバーロード用の finalize 済み CSR
				// (コピーは全体で 1 回 --- D-7。旧経路はループ内の公開
				// mul_v_im が毎回 as_csr() コピーを払っていた)。
				const base_type Acsr = static_cast<const base_type&>(*this).as_csr();
				const std::vector<index_type>& u3_ao = Acsr.outer_index();
				const std::vector<index_type>& u3_ai = Acsr.inner_index();
				const std::vector<interval_type>& u3_av = Acsr.values();
				const index_type u3_rows = Acsr.rowsize();
				const index_type u3_cols = Acsr.columnsize();
				// 完了行数カウンタ(D-5): 並列側は atomic fetch_add(relaxed)。
				// 意味は【完了行数】(先頭からの連続到達数ではない)。
				std::atomic<std::size_t> rows_done(0);
				std::atomic<bool> caught(false);
				std::exception_ptr eptr;
				const std::ptrdiff_t pn = static_cast<std::ptrdiff_t>(n);
				// work = n×n(暫定 --- 設計 §2)、入れ子並列は抑止(D-10)
				const bool par_on =
					(n * n >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD))
					&& !omp_in_parallel();
#pragma omp parallel if(par_on)
				{
					// Ri/RiA のスレッド私有化(設計 §3.2-1: 共有作業ベクトル
					// だった Ri/RiA は並列領域内で各スレッドが実体を確保)
					std::vector<_T> Ri_t;
					std::vector<interval_type> RiA_t;
#pragma omp for schedule(static)
					for (std::ptrdiff_t ip = 0; ip < pn; ip++) {
						try {
							const std::size_t i = static_cast<std::size_t>(ip);
							// 6.2-1: Ri = (Am⁻¹)_i(準備済み因子の行解コア ---
							// 共有 const 読みのみ。転置三角解の演算列は旧
							// impl の字句的移動 --- P1 で単体一致確認済み)。
							// 失敗は verification_error 側(設計 6.2-1)。
							const vcp::lu_apply_status rst =
								vcp::spmats_lu_extract_detail::lu_inverse_row_core_<_T, index_type>(
									Fprep, static_cast<index_type>(i), Ri_t);
							if (!(rst == vcp::lu_apply_status::success)) {
								vcp::throw_error<vcp::verification_error>(
									"spimats::policy_lss_with_info: inverse-row extraction failed (no claim)");
							}
							// 6.2-2: Rir[i] = Ri·r(逐次側と同一の演算列。
							// SR-1: 点因子の明示昇格は書かない)
							{
								interval_type acc(_T(0));
								for (std::size_t j = 0; j < n; j++) {
									acc += Ri_t[j] * r[j];
								}
								Rir[i] = acc;
							}
							// 6.2-3: RiA = Ri·[A](生配列オーバーロード ---
							// 逐次実装。並列カーネルへは到達しない = D-10 の
							// 入れ子抑止の実現形)。
							vcp::spimats_kernel::mul_v_im(
								Ri_t, u3_rows, u3_cols, u3_ao, u3_ai, u3_av, RiA_t);
							// RiA[i] から【厳密な 1】を減算(SR-1、R-3)。
							RiA_t[i] -= 1;
							// 全成分に abs(kv abs は |·| の区間包含)。
							for (std::size_t j = 0; j < n; j++) {
								RiA_t[j] = abs(RiA_t[j]);
							}
							// 6.2-4: RiAe[i] = Σ_j RiA[j](abs 済み行和)。
							{
								interval_type rs(_T(0));
								for (std::size_t j = 0; j < n; j++) {
									rs += RiA_t[j];
								}
								RiAe[i] = rs;
							}
							// 6.2-5: 進行記録(検死用 --- 主張ではない)
							rows_done.fetch_add(1, std::memory_order_relaxed);
						}
						catch (...) {
							// 例外は領域外 rethrow(採用は CAS で 1 回 ---
							// どの行の例外を採るかは不定 = 診断のみの差、D-8)
							bool expected = false;
							if (caught.compare_exchange_strong(expected, true)) {
								eptr = std::current_exception();
							}
						}
					}
				}
				// ループ後にメンバへ最終値を書く(意味 =「完了行数」--- D-5)
				last_verify_info.lss_rows_done = rows_done.load();
				if (caught.load()) std::rethrow_exception(eptr);
			}
#else
			for (std::size_t i = 0; i < n; i++) {
				// 6.2-1: Ri = (Am⁻¹)_i を取り出す(転置三角解・P/Q 込み
				// [確認済み spmats.hpp L1533-1536: row_i(A^{-1})^T =
				// P^T L^{-T} U^{-T} Q^T e_i])。失敗は verification_error 側
				// (設計 6.2-1)。
				const lu_apply_result rres = Am.policy_lu_inverse_row_with_info(
					Lf, Uf, pperm, qperm, static_cast<index_type>(i), Ri);
				if (!(rres.status == vcp::lu_apply_status::success)) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_lss_with_info: inverse-row extraction failed (no claim)");
				}
				// 6.2-2: Rir[i] = Ri·r(点×区間の内積。区画内直接ループ ---
				// 設計 6.1 Rir 行)。点因子の明示昇格は書かない(SR-1:
				// kv の混合演算子 [数値スカラー×区間 = acceptable_n 専用
				// 経路・2 乗算] が方向保証を担う --- R-2)。
				{
					interval_type acc(_T(0));
					for (std::size_t j = 0; j < n; j++) {
						acc += Ri[j] * r[j];
					}
					Rir[i] = acc;
				}
				// 6.2-3: RiA = Ri·[A](R2 mul_v_im --- ヘッダコメント自体が
				// 「RA − I の行ごと構成」用と明記 [確認済み
				// spimats_mixed_ops.hpp L15])。
				this->mul_v_im(Ri, RiA);
				// RiA[i] から【厳密な 1】を減算 → G = RA − I の第 i 行
				// (I は厳密。整数 1 の混合減算は kv acceptable_n 経路で
				// 厳密 --- SR-1、R-3)。
				RiA[i] -= 1;
				// 全成分に abs(区間の絶対値包含: kv abs は |·| の区間包含)。
				for (std::size_t j = 0; j < n; j++) {
					RiA[j] = abs(RiA[j]);
				}
				// 6.2-4: RiAe[i] = Σ_j RiA[j](abs 済み行和 = (|G|e)_i)。
				// RiAe は【二役】: max_i で ‖G‖∞ の材料 + 第 2 項の係数
				// ベクトル |G|e(訂正合意 ① --- 行和は RiAe だけが担う)。
				{
					interval_type rs(_T(0));
					for (std::size_t j = 0; j < n; j++) {
						rs += RiA[j];
					}
					RiAe[i] = rs;
				}
				// 6.2-5: 進行記録(検死用 --- 主張ではない)
				last_verify_info.lss_rows_done = i + 1;
			}
#endif

			// --- 6.3 判定と組み立て(全合成は区間演算 --- 丸めモード操作なし)---
			// 事前の端点有限性チェックは置かない(R-4、SR-2: T 汎用区画に
			// double cast / isfinite を書かない)。NaN・∞ の防御は
			//   - Gnorm 側: 下の certified 判定 !(sup < 1) が NaN/∞ を
			//     そのまま失敗側に落とす
			//   - Rir/err 側: err[i] 計算直後の案 A 出口ゲート(6.3-5)
			// の比較演算のみで行う。
			last_verify_info.stage = 3;                    // 判定段
			// 6.3-2: Gnorm = max_i RiAe[i](最大は .upper() 比較で選び、
			// 【区間のまま保持】--- オーナー方針: 以降の合成を区間演算で行う。
			// 採った行の区間 [lo*, up*] は max_i lo_i ≥ lo* により真の ‖G‖∞
			// ∈ [max_i lo_i, up*] ⊆ [lo*, up*] を包含する = 広い側で健全)。
			interval_type Gnorm = RiAe[0];
			for (std::size_t i = 1; i < n; i++) {
				if (RiAe[i].upper() > Gnorm.upper()) Gnorm = RiAe[i];
			}
			// 記録は double フィールドへの暗黙変換(設計 §3 の型規約。
			// フライトレコーダの縮小であり判定には使わない --- SR-2 の
			// grep 対象 static_cast は書かない)
			last_verify_info.lss_g_norm_sup = Gnorm.upper();
			// certified 判定(失敗側ゲート !(sup < 1) --- SLU-GT1 P1。NaN も
			// 失敗側)。不成立は無主張 throw。この判定は分母 1 − Gnorm の
			// 正値保証を兼ねる(sup < 1 なら radix-2 の Sterbenz 型減算で
			// 1 − sup の下端 > 0)。
			if (!(Gnorm.upper() < _T(1))) {
				vcp::throw_error<vcp::verification_error>(
					"spimats::policy_lss_with_info: ||RA - I||_inf >= 1, verification failed (no claim; see last_verify_info.lss_g_norm_sup)");
			}
			// 6.3-3: Rir 全成分に abs → Tnorm = max_i abs(Rir[i])(選び方は
			// Gnorm と同じ .upper() 比較・区間のまま保持)。abs と max は
			// 1 ループに融合(R-5 --- 走査 1 回で二役)。
			Rir[0] = abs(Rir[0]);
			interval_type Tnorm = Rir[0];
			for (std::size_t i = 1; i < n; i++) {
				Rir[i] = abs(Rir[i]);
				if (Rir[i].upper() > Tnorm.upper()) Tnorm = Rir[i];
			}
			last_verify_info.lss_t_norm_sup = Tnorm.upper();   // 同上(暗黙変換の記録)
			// 6.3-4: 乗数 β = Tnorm / (1 − Gnorm)(【訂正合意 ②】: 分子 =
			// 残差ノルム Tnorm、分母 = 1 − ‖G‖∞。区間除算 --- 分母は 6.3-2 の
			// 判定により certified 正)。Gnorm / (1 − Tnorm) 型の式は誤り
			// (書きそうになったら STOP --- 指示書 §3-2)。整数 1 の混合
			// 減算は kv acceptable_n 経路(SR-1、R-6)。
			const interval_type beta = Tnorm / (1 - Gnorm);
			// 6.3-5: err[i] = abs(Rir[i]) + β · RiAe[i](第 1 項 = 成分別
			// 残差項 |R(Ax̂−b)|_i、第 2 項 = β·(|G|e)_i。Rir は 6.3-3 で
			// abs 済み)。対称化 err[i] = [−upper, upper](フローの
			// 「lower = −upper」の実装形: 系 3.2.1 は |x*−x̂| の成分別上界
			// なので誤差区間は上端の対称区間)。
			std::vector<interval_type> err(n);
			for (std::size_t i = 0; i < n; i++) {
				const interval_type e = Rir[i] + beta * RiAe[i];
				const _T eu = e.upper();
				// 案 A 出口ゲート(R-4、設計 §6.3-4 根拠): err の上端は
				// abs 済み量(Rir・RiAe ⪰ 0)と正の乗数 β の積和なので、
				// 真正な値なら必ず ≥ 0。NaN は !(NaN ≥ 0) で失敗側に落ち、
				// これが T 汎用(SR-2: double cast / isfinite なし)の
				// 唯一の NaN 防御(比較演算 1 本のみ)。+∞ は ≥ 0 を
				// 通過するが err = [−∞, +∞] は包含として健全(無限幅 =
				// 情報ゼロの正直な包含)。
				if (!(eu >= _T(0))) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_lss_with_info: error bound is not a nonnegative number (no claim)");
				}
				err[i] = interval_type(-eu, eu);
			}
			// 6.3-6: result.x[i] = xh[i] + err[i](x̂ は点 → 点+区間の混合
			// 加算 --- kv acceptable_n 経路が方向保証を担う。SR-1、R-7)。
			// solution も基底規約に合わせて複製
			// [確認済み spmats_lu_factor.hpp L112: solution = x]。
			result.x.resize(n);
			for (std::size_t i = 0; i < n; i++) {
				result.x[i] = xh[i] + err[i];
			}
			result.solution = result.x;
			last_verify_info.stage = 4;                    // 成功
			return result;
		}
		// >>> END [SPI-U3] <<<

		// 【将来追記の方法(オーナー向け)】SPD 経路(is_spd + エネルギー
		// ノルム評価)や H 行列経路(対角優位性による可逆性証明)を追加する
		// 場合は、SPI-U3 区画の内側の分岐として実装する(保証コアの選択は
		// U3 の内側で行い、前段・器の初期化は共有する。新しい R 区画や
		// override は増やさない)。
		return result;   // not reached
	}

	// -------------------------------------------------------------------
	// override #2(設計 §3): 標準固有値は B = I の転送層。
	// U4 相当の記入はない(設計 §11: U4 = B=I 転送、記入なし)。
	// 基底 policy_eig_impl(全固有値)は policy_eigs_with_info(rowsize) 経由で
	// ここに自動到達する [確認済み spmats_eigs.hpp L4683]。
	// -------------------------------------------------------------------
	virtual eig_result<kv::interval<_T> > policy_eigs_with_info_impl(
		std::size_t k,
		const eig_options<kv::interval<_T> >& opt) const override
	{
		if (this->rowsize() != this->columnsize()) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_eigs_with_info: matrix must be square");
		}
		// B = I を厳密構成(SPI-R3)して一般化固有値 override へ転送。
		// 近似層(AP)の呼び出しは #3 側で B=I を厳密検出して標準
		// policy_eigs_with_info に戻る(spi_is_exact_identity_ 参照)。
		base_type IB;
		vcp::identity_interval(this->rowsize(), IB);
		return policy_generalized_eigs_with_info_impl(IB, k, opt);
	}

	// -------------------------------------------------------------------
	// [SPI-F2 F2-3] spi_inertia_interval_: 1×1 / 2×2 ブロック対角 D(点)の
	// 慣性を kv::interval<_T> の演算で再計算する(非仮想・static。U5 の
	// sep_eval が共有カーネル policy_inertia_from_block_diagonal の success
	// 判定の直後に呼び、以後の計数にはこちらの (neg, zero) を使う)。
	// 動機: 共有カーネル(spmats_base/spmats_ldl_impl.hpp L268-310)は 2×2
	// ブロックの det = dk·dk1 − e² を点演算で評価し、その符号で計数する。
	// double / mpfr は乗算・減算が正しく丸められ単調なので符号は逆転しない
	// (0 に潰れるだけ)が、kv::dd は非単調で近特異ブロックの符号が逆転する
	// (SPI-F2 P0 の dd 乱数探索: 10,800,000 組中 199,147 組で反転)。
	// ブロック構造の読み方は共有カーネルと同一(対角 (k,k)、副対角 |i−j|=1 が
	// 2×2 ブロック (k,k+1) の印、非格納対角は 0、|i−j|>1 の格納値は非ブロック
	// 対角 → false)。判定(設計 F2-3): 2×2 は det=[dk][dk1]−[e]²、tr=[dk]+[dk1]
	// で det.lower()>0 ∧ tr.lower()>0 → pos 2 / det.lower()>0 ∧ tr.upper()<0 →
	// neg 2 / det.upper()<0 → pos 1 neg 1 / それ以外 → zero 2。1×1 は [d] の符号
	// が確定しなければ zero。構造的空列は zero 1。戻り false = 非ブロック対角
	// または計数不整合(呼び出し側は失敗側 return 1 に落とす)。
	// 区間演算は kv::interval<_T> の有向丸め(rdouble / rdd / rmpfr)に依る。
	// 共有カーネル・is_spd には触れない(SPI-F2 設計 非対象)。
	// -------------------------------------------------------------------
	static bool spi_inertia_interval_(const spmats<_T, index_type>& D,
	                                  index_type& n_pos, index_type& n_neg, index_type& n_zero) {
		n_pos = index_type(0); n_neg = index_type(0); n_zero = index_type(0);
		if (D.rowsize() != D.columnsize() || D.rowsize() < index_type(0)) return false;
		if (!D.is_finalized()) D.finalize();
		const std::size_t un = static_cast<std::size_t>(D.rowsize());
		std::vector<_T> dv(un, _T(0)), ev(un > 0 ? un - 1u : 0u, _T(0));
		std::vector<char> dhas(un, char(0)), ehas(un > 0 ? un - 1u : 0u, char(0));
		{
			const std::vector<index_type>& outer = D.outer_index();
			const std::vector<index_type>& inner = D.inner_index();
			const std::vector<_T>& value = D.values();
			for (std::size_t o = 0; o < un; ++o) {
				for (index_type q = outer[o]; q < outer[o + 1u]; ++q) {
					const std::size_t i = static_cast<std::size_t>(inner[static_cast<std::size_t>(q)]);
					const _T& v = value[static_cast<std::size_t>(q)];
					if (i == o) { dv[o] = v; dhas[o] = char(1); }
					// [SPI-F3 F3-10 (i)] 副対角 (i, i+1) と (i+1, i) が両方格納されて
					// いれば値の厳密一致を要求(不一致 = 非対称 D → false)
					else if (i + 1u == o) { if (ehas[i] == char(0)) { ev[i] = v; ehas[i] = char(1); } else if (!(ev[i] == v)) return false; }
					else if (o + 1u == i) { if (ehas[o] == char(0)) { ev[o] = v; ehas[o] = char(1); } else if (!(ev[o] == v)) return false; }
					else return false;                                  // 非ブロック対角
				}
			}
		}
		std::size_t k = 0;
		while (k < un) {
			if (k + 1u < un && ehas[k] != char(0)) {
				// 2×2 ブロック (k, k+1)(非格納対角は 0)。[SPI-F3 F3-10 (ii)] ブロックを
				// 使ったら e_{k+1} が格納されていてはならない(三重対角の鎖 = 非ブロック
				// 対角 → false)
				if (k + 2u < un && ehas[k + 1u] != char(0)) return false;
				const kv::interval<_T> idk(dv[k]), idk1(dv[k + 1u]), ie(ev[k]);
				const kv::interval<_T> det = idk * idk1 - ie * ie;
				const kv::interval<_T> tr = idk + idk1;
				if (det.lower() > _T(0) && tr.lower() > _T(0)) { n_pos = n_pos + index_type(2); }
				else if (det.lower() > _T(0) && tr.upper() < _T(0)) { n_neg = n_neg + index_type(2); }
				else if (det.upper() < _T(0)) { n_pos = n_pos + index_type(1); n_neg = n_neg + index_type(1); }
				else { n_zero = n_zero + index_type(2); }                // 符号未確定 = 無主張側
				k += 2;
			}
			else if (dhas[k] != char(0)) {
				const kv::interval<_T> idk(dv[k]);
				if (idk.lower() > _T(0)) { n_pos = n_pos + index_type(1); }
				else if (idk.upper() < _T(0)) { n_neg = n_neg + index_type(1); }
				else { n_zero = n_zero + index_type(1); }
				k += 1;
			}
			else { n_zero = n_zero + index_type(1); k += 1; }            // 構造的空列
		}
		return static_cast<std::size_t>(n_pos + n_neg + n_zero) == un;
	}

	// -------------------------------------------------------------------
	// override #3(設計 §3): 一般化固有値 mid → AP 近似固有対 → 保証(SPI-U5)。
	// [SPI-F3 F3-1] 入口検査を spi_veigs_entry_checks_ に、本体を
	// spi_veigs_core_<Prec>(B, k, opt, M) に切り出した。_impl(M なし)は
	// M = null で、M 付き入口(public 区画の spi_veigs_dispatch_)は M 非 null
	// で同じコアを通る(M は近似段にだけ渡る)。
	// -------------------------------------------------------------------
	void spi_veigs_entry_checks_(const base_type& B, std::size_t k,
	                             const eig_options<interval_type>& opt) const
	{
		// 入口検査(misuse -> vcp::error 系)
		if (this->rowsize() != this->columnsize()) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_generalized_eigs_with_info: matrix must be square");
		}
		if (B.rowsize() != B.columnsize()) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_generalized_eigs_with_info: B must be square");
		}
		if (B.rowsize() != this->rowsize()) {
			vcp::throw_error<vcp::dimension_error>("spimats::policy_generalized_eigs_with_info: matrix sizes must match");
		}
		if (k == 0) {
			vcp::throw_error<vcp::invalid_argument>("spimats::policy_generalized_eigs_with_info: requested count must be positive");
		}
		if (opt.max_iter == 0) {
			vcp::throw_error<vcp::invalid_argument>("spimats::policy_generalized_eigs_with_info: max_iter must be positive");
		}
		if (!(opt.tol.upper() > _T(0))) {   // 失敗側ゲート(certified-only の書き方)
			vcp::throw_error<vcp::invalid_argument>("spimats::policy_generalized_eigs_with_info: tol must be positive");
		}
		if (!B.is_finalized()) B.finalize();   // 基底外殻経由なら finalize 済み(防御)
		// [SPI-F3 F3-11] 入口の対称性検査(D-2): [A]・[B] が区間行列として厳密
		// 対称(端点一致)でなければ domain_error。veigs.m は対称を前提とし、
		// 経路によって内部 is_spd の domain_error / 黙った返却に分かれていた
		// のを入口で統一する。nosymcheck 相当のノブは設けない。
		if (!policy_is_symmetric(static_cast<const base_type&>(*this)) || !policy_is_symmetric(B)) {
			vcp::throw_error<vcp::domain_error>("spimats::policy_generalized_eigs_with_info: A (or B) is not exactly symmetric as an interval matrix");
		}
	}

	virtual eig_result<kv::interval<_T> > policy_generalized_eigs_with_info_impl(
		const spmats<kv::interval<_T>, index_type>& B,
		std::size_t k,
		const eig_options<kv::interval<_T> >& opt) const override
	{
		spi_veigs_entry_checks_(B, k, opt);
		return spi_veigs_core_(B, k, opt, (const identity_preconditioner<_T>*)0);
	}

	// -------------------------------------------------------------------
	// [SPI-F3 F3-6] 入口 2 冪正規化の補助(非仮想・static)。
	// spi_max_abs_exp_: X の全格納成分の端点の max|x| = m·2^e(m ∈ [1/2, 1))の
	//   e − 1 を返す(正規化後の max|x| ∈ [1, 2))。零行列・非有限・double の
	//   表現範囲外は 0(正規化しない)。指数の選択は近似でよい(approximation
	//   only --- 主張は往復検査と区間演算が担う)。
	// spi_scale_pow2_: Y = 2^e·X を端点ごとの点乗算で作り、各端点で往復検査
	//   (x·2^e)·2^-e == x(両方向とも点乗算)を課す。1 成分でも不成立
	//   (subnormal 落ち・0 落ち・inf)なら false(呼び出し側が正規化を断念)。
	//   成立すれば Y の各成分は X の成分の厳密な 2^e 倍であり、区間として
	//   等価な問題(固有値が厳密に 2^e 倍)である。
	// -------------------------------------------------------------------
	static int spi_max_abs_exp_(const base_type& X) {
		const base_type Xc = X.as_csr();
		const std::vector<interval_type>& valv = Xc.values();
		double mx = 0.0;                                   // approximation only
		for (std::size_t p = 0; p < valv.size(); p++) {
			const double lo = static_cast<double>(valv[p].lower());   // approximation only
			const double up = static_cast<double>(valv[p].upper());   // approximation only
			if (!(lo == lo) || !(up == up)) return 0;      // NaN → 正規化しない
			const double alo = (lo < 0.0) ? -lo : lo;
			const double aup = (up < 0.0) ? -up : up;
			if (alo > mx) mx = alo;
			if (aup > mx) mx = aup;
		}
		if (!(mx > 0.0) || !(mx - mx == 0.0)) return 0;   // 零行列 / ∞ → 正規化しない
		int e = 0;
		(void)std::frexp(mx, &e);                          // mx = m·2^e, m ∈ [1/2, 1)
		e = e - 1;
		if (e > 1000 || e < -1000) return 0;               // 2^±e の双方が正規数である範囲のみ
		return e;
	}
	static bool spi_scale_pow2_(const base_type& X, int e, base_type& Y) {
		const _T s = _T(std::ldexp(1.0, e));              // 2^e(厳密)
		const _T si = _T(std::ldexp(1.0, -e));            // 2^-e(厳密)
		const base_type Xc = X.as_csr();
		const std::vector<index_type>& outerv = Xc.outer_index();
		const std::vector<index_type>& innerv = Xc.inner_index();
		const std::vector<interval_type>& valv = Xc.values();
		const index_type n = Xc.rowsize();
		Y.resize(n, Xc.columnsize());
		for (index_type i = 0; i < n; i++) {
			for (index_type p = outerv[static_cast<std::size_t>(i)];
			     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
				const interval_type& x = valv[static_cast<std::size_t>(p)];
				const _T lo = x.lower() * s;
				const _T up = x.upper() * s;
				if (!(lo * si == x.lower()) || !(up * si == x.upper())) return false;   // 往復検査
				Y.add(i, innerv[static_cast<std::size_t>(p)], interval_type(lo, up));
			}
		}
		Y.finalize();
		return true;
	}
	static bool spi_scale_pow2_scalar_(const _T& x, int e, _T& y) {
		const _T s = _T(std::ldexp(1.0, e));
		const _T si = _T(std::ldexp(1.0, -e));
		y = x * s;
		return (y * si == x);
	}

	// [SPI-F3 F3-1] 保証コア(旧 _impl 本体)。M は AP 呼び出しにだけ使う。
	template <class Prec>
	eig_result<kv::interval<_T> > spi_veigs_core_(
		const base_type& B_in,
		std::size_t k,
		const eig_options<kv::interval<_T> >& opt_in,
		const Prec* M) const
	{
		last_verify_info.clear();

		// [SPI-F3 F3-6] 入口 2 冪正規化: A' = 2^-a·A、B' = 2^-b·B(a, b は
		// max|端点| の二進指数、零行列なら 0)。以降の全処理(mid 抽出・AP・
		// shift LDL・分離・Rough・Lehmann)は正規化後の (A', B') で行い、
		// 固有値は μ = λ·2^(b−a)(厳密)なので出力で λ = μ·2^(a−b) に戻す
		// (区間乗算 --- 常に包含)。シフト σ も μ 尺度 σ' = σ·2^(b−a) に
		// 変換して渡す(点乗算 + 往復検査)。往復検査が 1 つでも不成立なら
		// 正規化全体を断念して a = b = 0(記録 veigs_scale_abandoned)。
		// 2 冪倍は丸めなし(往復検査で確認済み)なので区間問題として等価
		// であり、主張は変わらない。
		int sc_a = spi_max_abs_exp_(static_cast<const base_type&>(*this));
		int sc_b = spi_max_abs_exp_(B_in);
		base_type IAn, IBn;
		eig_options<kv::interval<_T> > opt_n = opt_in;
		bool sc_ok = (sc_a != 0 || sc_b != 0);
		if (sc_ok) sc_ok = spi_scale_pow2_(static_cast<const base_type&>(*this), -sc_a, IAn);
		if (sc_ok) sc_ok = spi_scale_pow2_(B_in, -sc_b, IBn);
		if (sc_ok) {                                   // σ は use_shift に関わらず target_* の中心として AP が参照する
			_T lo_s, up_s;
			sc_ok = spi_scale_pow2_scalar_(opt_in.shift.lower(), sc_b - sc_a, lo_s)
			     && spi_scale_pow2_scalar_(opt_in.shift.upper(), sc_b - sc_a, up_s);
			if (sc_ok) opt_n.shift = interval_type(lo_s, up_s);
		}
		if (!sc_ok) {
			last_verify_info.veigs_scale_abandoned = (sc_a != 0 || sc_b != 0);
			sc_a = 0; sc_b = 0;
			IAn = static_cast<const base_type&>(*this).as_csr();
			IBn = B_in.as_csr();
			opt_n = opt_in;
		}
		last_verify_info.veigs_scale_exp.assign(2, 0);
		last_verify_info.veigs_scale_exp[0] = sc_a;
		last_verify_info.veigs_scale_exp[1] = sc_b;
		const base_type& B = IBn;                          // 以降、B / IA は正規化後
		const eig_options<kv::interval<_T> >& opt = opt_n;
		const interval_type restore_scale(_T(std::ldexp(1.0, sc_a - sc_b)));   // 2^(a−b)(厳密)
		const _T prec_scale = _T(std::ldexp(1.0, sc_a));                      // adapter 用 2^a(F3-2)

		// mid 抽出(SPI-R1)→ AP で近似固有対。
		// B が厳密単位行列(#2 の B=I 転送を含む)のときは近似層のみ標準
		// policy_eigs_with_info に戻す(基底 AP の一般化ソルバーは B=I でも
		// 不収束を返すため。spi_is_exact_identity_ のコメント参照。保証コア
		// SPI-U5 が見る形は一般化のまま)。
		_AP Amid;
		vcp::mid_for_approximation(IAn, Amid);
		const eig_options<_T> popt = spi_extract_eig_options_(opt);
		// [SPI-SH1] B2 シフト LDL handle(SLDL-SH)と workspace を Amid/Bmid
		// 抽出直後に関数ローカルで構築(設計 SPI-SH1 §1(a)。クラスメンバ・
		// static は追加しない)。options は B2 既定 = supernodal×none
		// 【裁定 V-3】。σ 非依存工程(変換・certified 対称性検査・パターン
		// マージ・ordering・symbolic)はここで 1 回だけ実行され、後段
		// sep_eval の λ ごとには値組み立て+numeric のみが残る。
		// B が厳密単位行列(#2 の B=I 転送を含む)のときは B 省略形
		// (A − σI)を用いる(SLDL-SH 受け入れ 2: B 省略 ≡ B=I 明示の
		// 一致は B2 で実証済み)。
		vcp::ldl_shift_handle<_T, index_type> sh;
		typename vcp::ldl_shift_handle<_T, index_type>::workspace_type sh_ws;
		eig_result<_T> ap;
		// [SPI-F3 F3-1 / F3-2] M が非 null なら adapter で包んで AP の M 付き
		// 入口へ(基底 spmats_eigs.hpp の前処理付き shift-invert 経路)。null
		// なら従来どおり。AP が M 付きで converged=false(非対応 method 等)を
		// 返した場合は下の従来の state_error に落ちる。
		spi_prec_adapter_<Prec> prec_ad;
		prec_ad.M = M;
		prec_ad.scale = prec_scale;                        // (A' − σ'B')⁻¹ = 2^a·(A − σB)⁻¹
		// [SPI-F3 F3-7] AP への要求個数【対応: veigs.m L109-113】: k = 1 なら
		// min(10, n)、それ以外は min(n, k + 2)。余分な近似は選択集合 S の外側の
		// 境界で実アンカーとして使う(人工アンカーはリストの真の端のみ)。
		// result.requested_count は k のまま。
		const std::size_t n_sz = static_cast<std::size_t>(this->rowsize());
		const std::size_t k_req = (k == 1) ? ((n_sz < 10) ? n_sz : std::size_t(10))
		                                   : ((n_sz < k + 2) ? n_sz : k + 2);
		// 要求個数の梯子 k_req → min(n, k + 2) → k(重複は飛ばす): 基底 AP の
		// Lanczos 系は要求個数が n に近い / 大きいと「exhausted orthogonal
		// starting vectors」「mv_budget exhausted」で不収束を返すことがあり
		// (veigs.m が頼る MATLAB eigs との差)、その場合は少ない個数で
		// 再試行して k で通っていた入力を失わない。converged な最初の結果を
		// 採用する(近似の個数が減るだけで主張の形は不変)。全段不収束なら
		// 従来の state_error。
		const std::size_t k_ladder[3] = { k_req, (n_sz < k + 2) ? n_sz : k + 2, k };
		const bool b_is_identity = spi_is_exact_identity_(B);
		_AP Bmid;
		if (!b_is_identity) vcp::mid_for_approximation(B, Bmid);
		for (int li = 0; li < 3; li++) {
			const std::size_t kr = k_ladder[li];
			if (li > 0 && kr == k_ladder[li - 1]) continue;
			if (li == 2 && kr == k_ladder[0]) continue;
			if (b_is_identity) {
				ap = (M != 0) ? Amid.policy_eigs_with_info(kr, popt, prec_ad)
				              : Amid.policy_eigs_with_info(kr, popt);
			}
			else {
				ap = (M != 0) ? Amid.policy_generalized_eigs_with_info(Bmid, kr, popt, prec_ad)
				              : Amid.policy_generalized_eigs_with_info(Bmid, kr, popt);
			}
			if (ap.converged) { last_verify_info.veigs_ap_request_count = kr; break; }
		}
		if (b_is_identity) {
			sh = Amid.policy_ldl_shift_setup_with_info(
				vcp::ldl_shift_default_options<_T>());
		}
		else {
			sh = Amid.policy_ldl_shift_setup_with_info(
				Bmid, vcp::ldl_shift_default_options<_T>());
		}

		// スナップショット(設計 §5.2: 既定は近似固有値のみ。近似固有
		// ベクトルは opt-in --- 全固有値時に密 n×n が乗るため)
		last_verify_info.ap_eigenvalues = ap.eigenvalues;
		if (spi_opt.keep_ap_eigenvectors) {
			last_verify_info.ap_eigenvectors = ap.eigenvectors;
		}

		if (!ap.converged) {
			vcp::throw_error<vcp::state_error>(
				"spimats::policy_generalized_eigs_with_info: approximation eigensolver did not converge (see last_verify_info.ap_eigenvalues)");
		}

		// [SPI-SH1] setup が invalid のときは本関数の既存の失敗規約(no
		// claim 系 verification_error)で終える(新しい報告形式は作らない
		// --- 指示書 §1-1)。判定を AP 収束検査の後に置き、従来の失敗
		// 優先順位(近似層の失敗 = state_error が先)を保存する。
		if (!sh.valid()) {
			vcp::throw_error<vcp::verification_error>(
				"spimats::policy_generalized_eigs_with_info: shift LDL setup failed (no claim)");
		}

		// 成功時 result の器(設計 §5.3 有効性表):
		//   eigenvalues     = 保証付き固有値区間(SPI-U5 が埋める)
		//   eigenvectors    = 空のまま(veigs は固有値のみ対象)
		//   converged       = true(失敗はここに到達する前後で throw)
		//   requested_count = k / returned_count 等の個数系は SPI-U5 が
		//                     「返した包含の個数」として正直に埋める
		//   residuals_* / iterations / 内部診断 = 既定値のまま(G4)
		//   used_method     = "spimats" 系の識別子
		eig_result<kv::interval<_T> > result;
		result.converged = true;
		result.requested_count = k;
		result.used_method = "spimats_veigs";

		// >>> USER-IMPLEMENTS [SPI-U5: veigs 一般化固有値コア(U4 = B=I 転送、記入なし)] <<<
		//     正しさはオーナーが実装・保証。以下の想定コードコメントは
		//     補助であり拘束しない。定理由来箇所は UNVERIFIED マーク付き。
		//
		// 【想定フロー】近似固有値は last_verify_info.ap_eigenvalues(点型)。
		//   (1) 個数勘定の道具: 点側 _AP の policy_ldl_with_info /
		//       policy_inertia_with_info [確認済み spmats.hpp L1428-1496 相当:
		//       certified Bunch-Kaufman LDL + 慣性] を mid(A) − s·mid(B) に
		//       適用した慣性で、シフト s 近傍の固有値個数の近似勘定 → 包含の
		//       精密化に用いる(Sylvester 慣性則の適用条件・区間化の正当化は
		//       [UNVERIFIED: 一次資料照合が必要])。
		//   (2) シフト s の保証側判定は is_spd の二分法(設計 §7.1):
		//       保証付き区間演算で [A] − s[B] を構成(add_im_m / 点側 −s 倍)し、
		//       is_spd([A] − s[B]) で固有値の存在範囲を締める。
		//   (3) 幾何平均型の残差評価: imats::eigsymge [確認済み vcp/imats.hpp
		//       L419-522。冒頭コメントに Shinya Miyajima, "Numerical enclosure
		//       for each eigenvalue in generalized eigenvalue problem", JCAM 236,
		//       pp.2545-2552 (2012) Theorem 3 の引用があることを確認。定理の
		//       内容自体は UNVERIFIED(一次資料未照合)] の疎版が参考:
		//       Y(AX − BX·Λ) と YBX − I の ∞ ノルム評価から各固有値の包含を
		//       構成する。疎では X / Y が密になる点の扱い(列ごと処理等)は
		//       オーナー設計。
		//   (4) ノルム評価は基底 norminf / normone(SPFN)を使用(normtwo は
		//       保証に使用禁止)。判定前の端点有限性チェック(NaN 防御)は
		//       SPI-U3 の (5) と同文。
		//   (5) クラスタ・重複固有値のリカバリ(クラスタまとめ・シフト再選択)
		//       を行うなら本区画の内側で。尽くしてなお失敗なら
		//       vcp::throw_error<vcp::verification_error>(...)(無主張、§4.1)。
		//       部分成功(クラスタ単位)の返却はしない(設計 §13 で棄却済み)。
		//   (6) 成立時のみ result.eigenvalues に包含区間を格納し、
		//       returned_count / converged_count 等を「返した包含の個数」で
		//       正直に埋めて return result。
		//
		// ===============================================================
		// [SPI-D1 v1.1 実装] U5-a: veigs 分離段(Yamamoto 2001 LDLᵀ 慣性法)
		// 規範出典: 【対応: veigs.m(sandbox/reference/veigs-main/、
		//   Liu–Yanagisawa、MIT)。行番号は同ファイル】+ 設計
		//   SPI-D1_design.md §7.1。U5-b(Lehmann–Behnke 精密化 = veigs.m
		//   L576-635 と veig.m)は STOP-R 承認後に本区画へ追記する。
		// 実装範囲(v1 U5-a、[SPI-F3 F3-7] 改訂): veigs.m L109-113 の余分な
		//   近似(k = 1 なら 10 個、それ以外は k + 2 個)を取り、target 順の
		//   先頭 k 個を含むクラスタの和集合 = 選択集合 S だけを評価する
		//   (veigs.m Step 3 の逐次クラスタ拡張ループ L196-254 の「目標
		//   クラスタから外へ」に対応。S の外側のクラスタは実アンカーとして
		//   使い、人工アンカー【L312/L392】はリストの真の端にだけ置く)。
		//   ind(SIGMA 起点選択 L137-155)は S の target 順で再現する。
		// 主張の形(各クラスタ c ∈ S): 「一般化固有値問題 (A,B) の固有値のうち
		//   グローバル番号 [N_c+1, N_{c+1}] のものはすべて区間
		//   [dn_c, up_{c+1}] に含まれる」。番号勘定は Yamamoto 2001 の
		//   慣性移送 [UNVERIFIED: 定理内容の一次資料照合はオーナー
		//   (設計 §9)]、区間端は is_spd(U1)の certified 判定に基づく。
		//   [SPI-F3 F3-4] smallest_algebraic / largest_algebraic は、返す
		//   クラスタが端であること(N_{c_lo} = 0 / N_{c_hi+1} = n)を慣性で
		//   証明する(満たせなければ F3-5 の Rayleigh 退路)。
		//   smallest_magnitude / largest_magnitude / target_* は番号 N_c を
		//   last_verify_info.veigs_index_lo / hi に露出し、それが要求順位
		//   (最も近い / 最も大きい)であることは [SPI-F4 F4-1〜F4-4]
		//   spi_opt.target_rank_check に従って慣性で検査し
		//   last_verify_info.veigs_rank_certified(0 未検査 / 1 証明済み /
		//   −1 証明できず)に記録する。既定 record では返却は変えない
		//   (require でのみ −1 を無主張にする)。
		//   証明される順位は非厳密(返却固有値より近い未返却固有値がないことを証明する。同距離は許す)。
		// 丸め規律: 丸めモード操作なし(上界・下界は kv 区間 + 端点)。
		// ===============================================================
		{
			typedef spmats<_T, index_type> PM_;            // 点疎(基底)
			const index_type n = this->rowsize();
			const base_type& IA = IAn;                     // [SPI-F3 F3-6] 正規化後の [A']
			// lambda_ratio【対応: veigs.m L97/L271: 1E-5】。min_dist(L98)は
			// veigs.m 本体でも未使用(L273-275 のコメjoint参照)のため移植しない。
			const _T ratio = _T(1e-5);
			const _T one_m_ratio = _T(1) - ratio;          // (1−ratio)(閾値用の点計算)
			// [MPFR-TODO] meps = 2⁻⁵²(MATLAB eps)は double 前提の探索刻み
			// 下限。MPFR 解禁時は対象型のマシンイプシロンへ差し替え
			// (探索制御専用 --- 主張には関与しない)。
			// [SPI-F2 F2-2] λ_B_min の探索(下の RoughLower(B, I) 展開)では
			// 刻みの下限と停止判定を lb0 相対(meps·lb0)に変更済み(絶対値
			// 2⁻⁵² のままだと λmin(B) < 2⁻⁵² で負のシフトを試し、is_spd が
			// 自明に通って λ_B_min < 0 が確定する --- 監査 2026-10-14)。A 側
			// RoughLower / RoughUpper(L2208-2242 相当)の meps は符号に関与
			// しないので未変更。
			const _T meps = _T(std::ldexp(1.0, -52));      // MATLAB eps 相当【L522/L526】
			last_verify_info.stage = 1;                    // 前処理段

			// --- 区画内ローカルラムダ群(設計 §2-3)---
			// is_spd 呼び出し(既定 spd_options<_T>() --- veigs.m の isspd 既定と
			// 同じ)。判定対象 C を引数に取る区間版オーバーロードを *this 経由で
			// 呼ぶ(SPI-G1 D-5: 派生ポリシーの _impl override が veigs から効く。
			// 外殻は is_spd 所有フィールドのみクリアするので veigs の計数は
			// 消えない --- 設計 SPI-G1 §2.5)。呼び出し回数は
			// veigs_bisection_steps に計上(設計 §7.1 の診断)。
			auto check_spd = [&](const base_type& C) -> bool {
				last_verify_info.veigs_bisection_steps++;
				return this->is_spd(C, spd_options<_T>());   // virtual dispatch on *this (D-5)
			};
			// [A] − λ[B](要素ごとの区間演算。A・B が厳密対称なら結果も
			// 厳密対称 --- veigs.m の hull(M,M') 対称化【L528/L557】は
			// MATLAB の intval 積の丸め非対称への防御であり、要素演算のみの
			// ここでは不要)
			// [SPI-K2] K1-a リバートに伴う復元。基底 policy_scalar_mul は
			// const _T&(区間)のみ受けるため昇格が必要。非昇格の混合スカラーが将来
			// 必要なら、点型側を無制約でテンプレート化した spimats_base 自由関数
			// (R2 カーネル mul_v_im の _TP 方式と同型・R8 候補)として spimats 系内で
			// 実装する --- 門番・kv 特性参照は置かない(オーナー方針 2026-07-25)
			auto pencil_AB = [&](const _T& lam) -> base_type {
				return this->policy_sub(IA, this->policy_scalar_mul(B, interval_type(lam)));
			};
			// λ[B] − [A](RoughUpper 用【L556】)
			auto pencil_BA = [&](const _T& lam) -> base_type {
				return this->policy_sub(this->policy_scalar_mul(B, interval_type(lam)), IA);
			};
			// hull 対称化【対応: veigs.m L346/L420/L607/L622 の hull(M,M')】:
			// S_ij = S_ji = [min(lo_ij,lo_ji), max(up_ij,up_ji)](端点の
			// min/max のみ --- 算術なしで厳密。対称メンバーの包含を保存)。
			// 区間積由来の行列(diff_m 系)は演算順の差で端点対称が崩れ、
			// U1 の R4 厳密対称検査に落ちるため、ここで対称化してから渡す。
			auto sym_hull = [&](const base_type& M) -> base_type {
				const base_type Mc = M.as_csr();
				const std::vector<index_type>& outerv = Mc.outer_index();
				const std::vector<index_type>& innerv = Mc.inner_index();
				const std::vector<interval_type>& valv = Mc.values();
				// 鏡像値の取得(非格納 = [0,0])。R4 と同じ二分探索構造。
				auto mirror_at = [&](index_type i, index_type j, interval_type& out) -> void {
					out = interval_type(_T(0));
					const index_type first = outerv[static_cast<std::size_t>(i)];
					const index_type last = outerv[static_cast<std::size_t>(i) + 1];
					typename std::vector<index_type>::const_iterator b0 = innerv.begin() + first;
					typename std::vector<index_type>::const_iterator e0 = innerv.begin() + last;
					typename std::vector<index_type>::const_iterator it = std::lower_bound(b0, e0, j);
					if (it != e0 && *it == j) {
						out = valv[static_cast<std::size_t>(it - innerv.begin())];
					}
				};
				base_type S;
				S.resize(Mc.rowsize(), Mc.columnsize());
				for (index_type i = 0; i < Mc.rowsize(); i++) {
					for (index_type p = outerv[static_cast<std::size_t>(i)];
					     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
						const index_type j = innerv[static_cast<std::size_t>(p)];
						const interval_type& a = valv[static_cast<std::size_t>(p)];
						if (i == j) {
							S.add(i, i, a);
							continue;
						}
						interval_type bm2;
						mirror_at(j, i, bm2);
						// (i,j) と (j,i) の組は「両方格納」なら i<j の走査で
						// 1 回だけ、「片側のみ格納」なら格納側の走査で 1 回
						// だけ発行する(重複 add 防止)。
						const bool mirror_stored =
							!(bm2.lower() == _T(0) && bm2.upper() == _T(0));
						if (mirror_stored && i > j) continue;
						const _T lo = (a.lower() < bm2.lower()) ? a.lower() : bm2.lower();
						const _T up = (a.upper() > bm2.upper()) ? a.upper() : bm2.upper();
						const interval_type h(lo, up);
						S.add(i, j, h);
						S.add(j, i, h);
					}
				}
				S.finalize();
				return S;
			};

			// --- 入口: B の正定値性前提の検証(veigs.m ヘッダ契約 L7-8、
			// 設計 §7.1)。false = 前提を証明できない → 無主張 throw。---
			if (!check_spd(B)) {
				vcp::throw_error<vcp::verification_error>(
					"spimats::policy_generalized_eigs_with_info: B could not be verified positive definite (no claim)");
			}

			// --- λ_B_min: λmin(B) の certified 下界【対応: veigs.m L167-168】---
			// 近似値: B = I なら厳密に 1(spi_is_exact_identity_ --- #2 転送の
			// 標準固有値問題)。それ以外は AP eigs(smallestabs)【L167】。
			// certified 化: is_spd([B] − λI) の探索【L168 → L518-541
			// veigs_RoughLower(B, speye, λ0, λ0/2) の展開移植】。
			_T lambda_B_min;
			{
				_T lb0;
				if (spi_is_exact_identity_(B)) {
					lb0 = _T(1);
				}
				else {
					_AP Bmid2;
					vcp::mid_for_approximation(B, Bmid2);
					eig_options<_T> bo;
					bo.target = vcp::eig_target::smallest_magnitude;   // 'smallestabs'【L167】
					const eig_result<_T> br = Bmid2.policy_eigs_with_info(1, bo);
					if (!br.converged || br.eigenvalues.empty()) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: lambda_min(B) approximation failed (no claim)");
					}
					lb0 = br.eigenvalues[0];
				}
				if (!(lb0 > _T(0))) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: lambda_min(B) approximation is not positive (no claim)");
				}
				// RoughLower(B, I) 展開【対応: veigs.m L518-541。探索制御は
				// 点計算で可(主張は is_spd の成功のみが担う)】
				base_type II;
				vcp::identity_interval(n, II);
				const _T tmin = lb0 / _T(2);               // eig_test_min = λ0/2【L168】
				// [SPI-F2 F2-2 (i)] 刻みの下限は lb0 相対(meps·lb0)。lb0 > 0 は
				// 検査済みなので、試す lam はすべて tmin − meps·lb0 > 0 以上に
				// なる(lb0 > 2⁻¹⁹ では旧 max(eps, ·) と同一の刻み --- ビット一致)。
				const _T meps_rel = meps * lb0;
				_T delta = (lb0 - tmin) / _T(4294967296.0);   // (test−min)/2^32【L522】
				if (!(delta > meps_rel)) delta = meps_rel;   // max(eps·lb0, ·)【L522 の相対化】
				_T step = delta;
				_T lam = lb0 - step;                       // λ = test − δ·2^{k−1}【L523】
				bool okB = false;
				int k = 1;
				while (true) {                             // 【L526-534】
					// [SPI-K2] K1-a リバートに伴う復元。基底 policy_scalar_mul は
					// const _T&(区間)のみ受けるため昇格が必要。非昇格の混合スカラーが将来
					// 必要なら、点型側を無制約でテンプレート化した spimats_base 自由関数
					// (R2 カーネル mul_v_im の _TP 方式と同型・R8 候補)として spimats 系内で
					// 実装する --- 門番・kv 特性参照は置かない(オーナー方針 2026-07-25)
					const base_type Pm = this->policy_sub(B,
						this->policy_scalar_mul(II, interval_type(lam)));
					if (check_spd(Pm)) { okB = true; break; }   // certified λmin(B) > lam
					k++;
					step = step * _T(2);
					lam = lb0 - step;
					if (!(lam >= tmin - meps_rel) || k > 64) break;   // 停止【L526 の相対化(G-0 裁定 4-2)。k 上限は区画内防御定数】
				}
				if (!okB) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: failed to certify a lower bound of lambda_min(B) (no claim)");   // 【対応: L539】
				}
				lambda_B_min = lam;
				// [SPI-F2 F2-2 (ii)] 正値性の防御ゲート。(i) の後では到達しない
				// はず(到達 = 設計の誤り → 報告)。err = ‖diff‖/λ_B_min と
				// Lehmann の Err = ‖e‖²/λ_B_min の符号はこの正値性に依存する。
				if (!(lambda_B_min > _T(0))) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: certified lower bound of lambda_min(B) is not positive (no claim)");
				}
			}

			// --- 粗い探索範囲【対応: veigs.m L169-172】(探索制御専用の
			// 点計算。主張には関与しない)---
			_T lambda_rough_min, lambda_rough_max;
			{
				base_type NB;
				IA.norminf(NB);                            // ‖sup|A|‖∞ 相当【L169】
				const _T An = NB.get(0, 0).upper();
				lambda_rough_min = -(An / lambda_B_min);   // 【L170】
				lambda_rough_max = (An / lambda_B_min);    // 【L171】
			}

			// --- 近似固有値のソート【対応: veigs.m L132: sort(diag(D))】---
			std::vector<_T> ev = ap.eigenvalues;
			if (ev.empty()) {
				vcp::throw_error<vcp::verification_error>(
					"spimats::policy_generalized_eigs_with_info: AP returned no eigenvalues to verify (no claim)");
			}
			std::sort(ev.begin(), ev.end());
			const std::size_t m = ev.size();

			// --- クラスタ分割【対応: veigs.m L306(右規則)/L387(左規則):
			// 相対ギャップ < lambda_ratio·|隣接値| なら同一クラスタ】---
			std::vector<std::size_t> clo, chi;             // 各クラスタの [lo, hi](ev 添字)
			{
				using std::abs;
				std::size_t s0 = 0;
				for (std::size_t j = 0; j + 1 <= m - 1; j++) {
					// [SPI-F3 F3-8] 同一クラスタ ⇔ ev[j] >= ev[j+1] − ratio·|ev[j+1]|
					// (等しい値 --- 厳密 0 の重複を含む --- は常に同一クラスタ)
					if (!(ev[j] >= ev[j + 1] - ratio * abs(ev[j + 1]))) {   // 【L306 の否定 = クラスタ境界】
						clo.push_back(s0);
						chi.push_back(j);
						s0 = j + 1;
					}
				}
				clo.push_back(s0);
				chi.push_back(m - 1);
			}
			const std::size_t C = clo.size();
			last_verify_info.veigs_clusters = C;

			// --- ソート順 → ap 原順の添字置換(U5-b の V 列対応と Rayleigh 退路の
			// v の取得に使う)---
			std::vector<std::size_t> ord(m);
			for (std::size_t t = 0; t < m; t++) ord[t] = t;
			std::sort(ord.begin(), ord.end(),
				[&](std::size_t x, std::size_t y) {
					return ap.eigenvalues[x] < ap.eigenvalues[y];
				});

			// --- [SPI-F3 F3-7] 選択集合 S【対応: veigs.m L137-155 の ind と
			// L196-254 の拡張】: target 順(sa: 昇順 / la: 降順 / sm: |ev| 昇順 /
			// lm: |ev| 降順 / target_*: |ev − σ| 昇順。σ は opt.shift の mid ---
			// AP の select_real_eigenpairs と同じく use_shift に関わらず参照)で
			// 先頭 k 個(m < k なら全部)を取り、それらを含むクラスタの最小番号
			// c_lo から最大番号 c_hi までを S とする(S は連続。間のクラスタも
			// 検証・返却する)。
			const bool tg_sa = (opt.target == vcp::eig_target::smallest_algebraic);
			const bool tg_la = (opt.target == vcp::eig_target::largest_algebraic);
			std::size_t c_lo = 0, c_hi = C - 1;
			{
				using std::abs;
				std::vector<std::size_t> cl_of(m, 0);          // ev 添字 → クラスタ番号
				for (std::size_t c = 0; c < C; c++)
					for (std::size_t t = clo[c]; t <= chi[c]; t++) cl_of[t] = c;
				std::vector<std::size_t> rank(m);              // ソート済み ev の添字を target 順に
				for (std::size_t t = 0; t < m; t++) rank[t] = t;
				const _T sigma_t = mid(opt.shift);
				switch (opt.target) {
				case vcp::eig_target::smallest_algebraic:
					break;                                     // 昇順のまま
				case vcp::eig_target::largest_algebraic:
					std::reverse(rank.begin(), rank.end());
					break;
				case vcp::eig_target::smallest_magnitude:
					std::stable_sort(rank.begin(), rank.end(),
						[&](std::size_t x, std::size_t y) { return abs(ev[x]) < abs(ev[y]); });
					break;
				case vcp::eig_target::largest_magnitude:
					std::stable_sort(rank.begin(), rank.end(),
						[&](std::size_t x, std::size_t y) { return abs(ev[x]) > abs(ev[y]); });
					break;
				default:                                       // target_real / target_magnitude
					std::stable_sort(rank.begin(), rank.end(),
						[&](std::size_t x, std::size_t y) { return abs(ev[x] - sigma_t) < abs(ev[y] - sigma_t); });
					break;
				}
				const std::size_t ksel = (k < m) ? k : m;
				c_lo = C; c_hi = 0;
				for (std::size_t t = 0; t < ksel; t++) {
					const std::size_t c = cl_of[rank[t]];
					if (c < c_lo) c_lo = c;
					if (c > c_hi) c_hi = c;
				}
			}
			// [SPI-F3 F3-8] 人工アンカーの幅 h: 0.01·|b| が 0 なら
			// 0.01·max(|ev[0]|, |ev[m−1]|)、それも 0 なら max(meps, (λ_rough_max −
			// λ_rough_min)·2^-32)(厳密 0 の固有値で窓が潰れないため)
			auto anchor_width = [&](const _T& b) -> _T {
				using std::abs;
				_T h = abs(b) * _T(0.01);
				if (!(h > _T(0))) {
					const _T e0 = abs(ev[0]), e1 = abs(ev[m - 1]);
					h = ((e0 > e1) ? e0 : e1) * _T(0.01);
				}
				if (!(h > _T(0))) {
					h = (lambda_rough_max - lambda_rough_min) / _T(4294967296.0);
					if (!(h > meps)) h = meps;
				}
				return h;
			};

			// --- Yamamoto 分離子評価(区画内ローカルラムダ)---
			// 【対応: veigs.m L318-375(右側)/L397-448(左側)は同一の
			//  評価器の二用途 --- 本実装は 1 本に統合し、境界ごとに
			//  「下クラスタ寄り(左式 L397)」「上クラスタ寄り(右式 L318)」の
			//  2 シフトを評価する(veigs.m の歩行が行うのと同じ 2 回)】
			// 手順: tmpM = [A] − λ[B]【L319/L398】→ mid の certified BK LDL
			// 【L321/L400 の ldl。policy_ldl_with_info = LDL-4】→ D ブロックの
			// certified 慣性【L323/L401 の GetInertia ⇔
			// policy_inertia_from_block_diagonal --- 1×1/2×2 計数の対応】→
			// diff = P L D Lᵀ Pᵀ − tmpM の区間評価【右側 L334 と同じ向き。左側
			// L411 は diff_m = tmpM − PLDLᵀPᵀ と符号が反転するが、本実装は diff を
			// 1 本に保ち、半分割の試験行列側で符号を切り替える(SPI-F2 F2-1)】→
			// err_est = ‖diff‖∞ 上端 / λ_B_min【L335/L412。初期値の無検査
			// 受理は対称 diff で ‖diff‖₂ ≤ ‖diff‖∞ と λmin(B) ≥ λ_B_min に
			// よる err·B 支配 --- [UNVERIFIED: Yamamoto 2001 の慣性移送
			// 定理の内容・適用条件はオーナー照合]】→ 半分割精密化
			// 【L340-353/L414-427: is_spd(err/2·B ∓ diff) が通る間 err を
			// 半減 --- 設計 §7.1 の「σ の二分」の一形態。符号は side 依存】。
			// [SPI-F2 F2-1] upper_side: false = 右側(λR、dn = λR − err が
			// λ_{N+1} の下界)/ true = 左側(λL、up = λL + err が λ_N の上界)。
			// 数学: K_pt = PLDLᵀPᵀ の負固有値数を N とすると、右側の主張は
			// F + err·B ⪰ K_pt ⇔ err·B − (K_pt − F) ⪰ 0、左側の主張は
			// F − err·B ⪯ K_pt ⇔ err·B + (K_pt − F) ⪰ 0(対称行列の固有値の
			// 単調性)。diff = K_pt − F なので試験行列は右側 err·B − diff、
			// 左側 err·B + diff。半分割前の初期値 err = ‖diff‖∞/λ_B_min は
			// err·B ⪰ ‖diff‖∞·I ⪰ ±diff で両向きに成立し、変更は半分割が
			// 走った左側にだけ効く。参照元 veigs.m は L334(右)/ L411(左)で
			// diff_m の符号を反転させ両側とも err·B − diff_m を検査しており、
			// C++ 移植で左側の反転が落ちていた(監査 2026-10-14 CONFIRMED)。
			// 戻り: 0 = 計数成立 / 1 = 分解・慣性の失敗 / 2 = ギャップ不足。
			auto sep_eval = [&](const _T& lam, const _T& gap,
			                    index_type& neg_out, index_type& zero_out,
			                    _T& err_out, bool upper_side) -> int {
				const base_type tm = pencil_AB(lam);   // (v) の区間残差が参照するため λ ごとの構築は残す
				PM_ Lf, Df;
				std::vector<index_type> lperm;
				// [SPI-SH1] (ii) per-λ mid 抽出+(iii) policy_ldl_with_info
				// (既定 auto→baseline_dynamic×bk)を shift handle の
				// factor_at に置換【対応: veigs.m L321/L400 の ldl --- 旧経路
				// の LDL-4 と同じ役割】。mid の線形性 mid([A]−λ[B]) =
				// Amid − λ·Bmid(設計 SPI-SH1 §0)により、λ ごとに実行される
				// のは値組み立て(O(nnz) axpy)+ supernodal numeric のみ。
				// 後段 (v) の区間残差は「得られた分解そのもの」を tm に対して
				// 評価する自己完結構造のため、分解経路の変更は主張の正当性を
				// 変えない(設計 SPI-SH1 §2)。inertia_at ではなく factor_at
				// を使う(保証の連鎖 (v) が L/D/perm の実体を要求するため)。
				ldl_result<_T, index_type> lr =
					sh.factor_at(lam, Lf, Df, lperm, sh_ws);
				if (!(lr.status == vcp::sparse_ldl_status::success
				      || lr.status == vcp::sparse_ldl_status::zero_pivot)) {
					// [SPI-SH1] 明示再試行【裁定 V-1】: D-6 の文書化済み
					// 脱出路(silent fallback ではない --- handle は自動切替
					// せず、再試行はこの呼び出し箇所に明示的に書く)。既定
					// supernodal×none が success/zero_pivot 以外を返したとき、
					// 同じ λ を baseline_dynamic×bk(従来経路と同じ組)で
					// 1 回だけ再分解し、その status で改めて判定する。
					// これで従来経路の成功包絡を完全に保存する。それも失敗
					// なら従来どおり return 1。診断カウンタは追加しない
					// 【裁定 V-2】。
					lr = sh.factor_at_with(lam,
						vcp::sparse_ldl_method::baseline_dynamic,
						vcp::sparse_ldl_pivoting::bk, Lf, Df, lperm, sh_ws);
					if (!(lr.status == vcp::sparse_ldl_status::success
					      || lr.status == vcp::sparse_ldl_status::zero_pivot)) {
						return 1;                          // 【対応: L369/L442 の error 側】
					}
				}
				const inertia_result<index_type> ir =
					Df.policy_inertia_from_block_diagonal(_T(0));   // GetInertia 対応【L323/L401】
				if (!(ir.status == vcp::inertia_status::success)) return 1;
				// [SPI-F2 F2-3] 慣性は区間で再計算した値を使う(共有カーネルの
				// (neg, zero) は success 判定にのみ用い、計数には参照しない)。
				// 非ブロック対角 / 計数不整合は従来どおり失敗側。
				{
					index_type ipos = index_type(0), ineg = index_type(0), izero = index_type(0);
					if (!spi_inertia_interval_(Df, ipos, ineg, izero)) return 1;
					neg_out = ineg;
					zero_out = izero;
				}
				// diff_m = P·(L·D·Lᵀ)·Pᵀ − tmpM【対応: L334/L411。LDL 規約は
				// PᵀAP = LDLᵀ・perm new->old [確認済み spmats.hpp L1431-1437
				// 相当] なので積の (i,j) 成分は (perm[i], perm[j]) へ写す】
				base_type LI, DI;
				vcp::point_to_interval(Lf, LI);
				vcp::point_to_interval(Df, DI);
				const base_type T1 = this->policy_mul(LI, DI);
				const base_type T2 = this->policy_mul(T1, LI.transpose());
				base_type Pm;
				Pm.resize(n, n);
				{
					const base_type T2c = T2.as_csr();
					const std::vector<index_type>& outerv = T2c.outer_index();
					const std::vector<index_type>& innerv = T2c.inner_index();
					const std::vector<interval_type>& valv = T2c.values();
					for (index_type i = 0; i < n; i++) {
						for (index_type p = outerv[static_cast<std::size_t>(i)];
						     p < outerv[static_cast<std::size_t>(i) + 1]; p++) {
							Pm.add(lperm[static_cast<std::size_t>(i)],
							       lperm[static_cast<std::size_t>(innerv[static_cast<std::size_t>(p)])],
							       valv[static_cast<std::size_t>(p)]);
						}
					}
				}
				Pm.finalize();
				const base_type diff = this->policy_sub(Pm, tm);
				// err_est = ‖diff‖∞ 上端 / λ_B_min【L335/L412】(SPFN norminf、
				// 除算は区間で上端採用)
				base_type NB;
				diff.norminf(NB);
				const _T dnrm = NB.get(0, 0).upper();
				// 案 A 型ゲート(R-9、SR-2): dnrm はノルムの上端なので真正な
				// 値なら非負必然。NaN のみ !(NaN ≥ 0) で失敗側(return 1)へ
				// 落ち、+∞ はそのまま通して後段(err > thr → return 2 の
				// 失敗側)の保守性に委ねる --- 比較演算 1 本のみで double
				// cast / isfinite を書かない(T 汎用区画)。
				if (!(dnrm >= _T(0))) return 1;
				// 【必要な昇格(SR-1 対象外 --- 区別コメント)】点÷点の除算は
				// 最近接丸めで上界の方向保証が残らないため、両辺を区間化して
				// 区間除算の .upper() を採る。これは冗長昇格ではなく
				// err_est ≥ ‖diff‖∞/λ_B_min を保証する唯一の書き方。
				_T err = (interval_type(dnrm) / interval_type(lambda_B_min)).upper();
				// 半分割精密化【L340-353/L414-427】: ギャップ不足の間、
				// is_spd(err/2·[B] − diff)(右側)/ is_spd(err/2·[B] + diff)
				// (左側)が certified で通るなら err を半減(SPI-F2 F2-1)。
				// (試験行列は区間積由来で端点対称が崩れるため hull 対称化
				// 【L346/L420】してから U1 に渡す)
				const _T thr = one_m_ratio * gap;
				if (err > thr) {
					int hk = 0;
					while (hk < 64) {                      // 上限は区画内防御定数
						hk++;
						const _T err_t = err / _T(2);
						// [SPI-K2] K1-a リバートに伴う復元。基底 policy_scalar_mul は
						// const _T&(区間)のみ受けるため昇格が必要。非昇格の混合スカラーが将来
						// 必要なら、点型側を無制約でテンプレート化した spimats_base 自由関数
						// (R2 カーネル mul_v_im の _TP 方式と同型・R8 候補)として spimats 系内で
						// 実装する --- 門番・kv 特性参照は置かない(オーナー方針 2026-07-25)
						const base_type errB = this->policy_scalar_mul(B, interval_type(err_t));
						const base_type tmpMx = upper_side
							? sym_hull(this->policy_add(errB, diff))      // 左側: err·B + diff【L411 の反転 + L420】
							: sym_hull(this->policy_sub(errB, diff));     // 右側: err·B − diff【L334 + L346】
						if (!check_spd(tmpMx)) break;      // 【L348-350 / L421-423】
						err = err_t;                       // 【L351】
						if (!(err > thr)) break;           // 十分縮んだら終了
					}
				}
				err_out = err;
				if (err > thr) return 2;                   // 実現不能【L357-359/L430-432】
				return 0;
			};

			// --- 各境界の評価 ---
			// 境界 j = 0..C: N[j] = 「境界より下の固有値の個数」の certified
			// 勘定、dn[j] = λ_{N[j]+1} の下界、up[j] = λ_{N[j]} の上界。
			// クラスタ c の包含 = [dn[c], up[c+1]]、番号 = [N[c]+1, N[c+1]]。
			std::vector<index_type> Nj(C + 1, index_type(0));
			std::vector<_T> dnj(C + 1, _T(0)), upj(C + 1, _T(0));
			std::vector<bool> have_dn(C + 1, false), have_up(C + 1, false);
			last_verify_info.stage = 2;                    // ノルム評価(分離)段
			using std::abs;
			// --- Rough 経路(certified is_spd 直接判定)を区画内ラムダに
			// (分離子ループの端境界と F3-5 の Rayleigh 退路が共用)---
			// 全固有値 > λ_RL【対応: L518-541 veigs_RoughLower(A,B)】
			auto rough_lower = [&](_T& lam_out) -> bool {
				_T delta = (ev[clo[0]] - lambda_rough_min) / _T(4294967296.0);
				if (!(delta > meps)) delta = meps;
				_T step = delta;
				_T lam = ev[clo[0]] - step;
				int kk = 1;
				while (true) {
					if (check_spd(pencil_AB(lam))) { lam_out = lam; return true; }
					kk++;
					step = step * _T(2);
					lam = ev[clo[0]] - step;
					if (!(lam >= lambda_rough_min - meps) || kk > 64) return false;
				}
			};
			// 全固有値 < λ_RU【対応: L547-570 veigs_RoughUpper(A,B)】
			auto rough_upper = [&](_T& lam_out) -> bool {
				_T delta = (lambda_rough_max - ev[chi[C - 1]]) / _T(4294967296.0);
				if (!(delta > meps)) delta = meps;
				_T step = delta;
				_T lam = ev[chi[C - 1]] + step;
				int kk = 1;
				while (true) {
					if (check_spd(pencil_BA(lam))) { lam_out = lam; return true; }
					kk++;
					step = step * _T(2);
					lam = ev[chi[C - 1]] + step;
					if (!(lam <= lambda_rough_max + meps) || kk > 64) return false;
				}
			};
			// [SPI-F3 F3-5] sa / la では分離段(境界評価・順位検査・整合検査)の
			// 無主張を捕捉して Rayleigh 退路を試みる。退路の結果はこの 2 変数。
			bool ray_fb = false;
			interval_type ray_encl(_T(0));
			try {
			// [SPI-F3 F3-7] 評価する境界は S に隣接する j ∈ [c_lo, c_hi + 1] のみ。
			// bottom / top(人工アンカー・Rough 経路)は S が端のクラスタを含む
			// ときだけ到達し、S の外側にクラスタがある境界は隣接近似値を
			// 実アンカーとして従来式で評価する。
			for (std::size_t j = c_lo; j <= c_hi + 1; j++) {
				const bool bottom = (j == 0);
				const bool top = (j == C);
				// 境界の下側アンカー a(下クラスタの上端)と上側アンカー b
				// (上クラスタの下端)。端は人工アンカー【L312/L392: ±1%、幅 h は F3-8】。
				_T a, b2;
				if (bottom) {
					b2 = ev[clo[0]];
					a = b2 - anchor_width(b2);             // eig_pre【L392】
				}
				else if (top) {
					a = ev[chi[C - 1]];
					b2 = a + anchor_width(a);              // eig_next【L312】
				}
				else {
					a = ev[chi[j - 1]];
					b2 = ev[clo[j]];
				}
				const _T gap = b2 - a;
				if (!(gap > _T(0))) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: degenerate separator gap (no claim)");
				}
				// --- 上側寄りシフト(右式【L318】: λR = ratio·a + (1−ratio)·b)
				//     → dn[j] = inf(λR − err) = λ_{N+1} の下界【L373-374】。
				//     最下端境界では省略せず評価し、neg==0 なら RoughLower
				//     経路【L403-408】に切替える。
				// --- 下側寄りシフト(左式【L397】: λL = ratio·b + (1−ratio)·a)
				//     → up[j] = sup(λL + err) = λ_N の上界【L446-447】。
				//     最上端境界では neg+zero==n のとき RoughUpper 経路
				//     【L325-331】に切替える。
				// 2 シフトの計数一致は境界の整合条件(下の consistency 検査)。
				index_type negR = 0, zeroR = 0, negL = 0, zeroL = 0;
				_T errR = _T(0), errL = _T(0);
				bool haveR = false, haveL = false;
				bool roughL = false, roughU = false;
				{
					const _T lamR = ratio * a + one_m_ratio * b2;   // 【L318】
					const int rc = sep_eval(lamR, gap, negR, zeroR, errR, false);   // 右側【L318-375】
					if (rc == 0) {
						if (bottom && negR == 0) {
							roughL = true;                 // 最小側【L403-405 相当】
						}
						else if (top && negR + zeroR == n) {
							roughU = true;                 // 最大側【L325-329】
						}
						else if (zeroR != 0) {
							vcp::throw_error<vcp::verification_error>(
								"spimats::policy_generalized_eigs_with_info: zero inertia block at separator (no claim)");
						}
						else {
							haveR = true;
							Nj[j] = negR;                  // global_s = neg_num【L373】
							// 【必要な昇格(片側のみ --- SR-1)】点−点に方向保証を
							// 持たせる区間化。もう片方は kv 混合減算に任せる(R1 手直し)
							// [SPI-C1] これは SR-1 の【必要な昇格】: lamR/errR は
							// いずれも点型であり、点−点を最近接丸めで行うと方向保証が
							// 失われる。区間化して .lower() を取ることで検証済み下界に
							// なる(冗長昇格ではない --- 旧記録の「冗長」分類は誤り)。
							dnj[j] = (interval_type(lamR) - errR).lower();   // 【L374】
							have_dn[j] = true;
						}
					}
					else if (rc == 2 && bottom) {
						roughL = true;                     // ギャップ不足 → 粗い下界【L361-367 の縮退】
					}
					else if (rc == 2 && top) {
						roughU = true;                     // 【L434-440 の縮退】
					}
					else {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: separator certification failed (no claim)");   // 【L369/L442】
					}
				}
				if (!roughL && !roughU) {
					const _T lamL = ratio * b2 + one_m_ratio * a;   // 【L397】
					const int rc = sep_eval(lamL, gap, negL, zeroL, errL, true);    // 左側【L397-448】
					if (rc == 0 && zeroL == 0) {
						haveL = true;
						// 【必要な昇格(片側のみ --- SR-1)】同上: 点+点の上端保証
						// [SPI-C1] これは SR-1 の【必要な昇格】: lamL/errL は
						// いずれも点型であり、点+点を最近接丸めで行うと方向保証が
						// 失われる。区間化して .upper() を取ることで検証済み上界に
						// なる(冗長昇格ではない --- 旧記録の「冗長」分類は誤り)。
						upj[j] = (interval_type(lamL) + errL).upper();   // global_rho = sup(λ+err)【L447】
						have_up[j] = true;
						// 2 シフトの計数一致(境界の両側評価は同じ個数を
						// 見ているはず --- 崩れたら近似が粗い = 無主張)
						if (haveR && !(negL == negR)) {
							vcp::throw_error<vcp::verification_error>(
								"spimats::policy_generalized_eigs_with_info: inconsistent inertia counts at separator (no claim)");
						}
						if (!haveR) { Nj[j] = negL; }
					}
					else {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: separator certification failed (lower-side shift, no claim)");
					}
				}
				// --- Rough 経路(certified is_spd 直接判定)---
				if (roughL) {
					_T lam;
					if (!rough_lower(lam)) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: rough lower bound search failed (no claim)");   // 【L539】
					}
					Nj[j] = 0;                             // certified: 下に固有値なし
					dnj[j] = lam;
					have_dn[j] = true;
					upj[j] = lam;                          // 下端境界の up は未使用
				}
				if (roughU) {
					_T lam;
					if (!rough_upper(lam)) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: rough upper bound search failed (no claim)");   // 【L568】
					}
					Nj[j] = n;                             // certified: 上に固有値なし
					upj[j] = lam;
					have_up[j] = true;
					dnj[j] = lam;                          // 上端境界の dn は未使用
				}
			}

			// --- [SPI-F3 F3-4] 要求順位の検査【対応: veigs.m L182-194 の端固定と
			// L459 の global_r / global_s 検査】: smallest_algebraic は最下端境界の
			// 慣性 N_{c_lo} = 0(RoughLower 経路または negR = 0 の certified)、
			// largest_algebraic は N_{c_hi+1} = n を要求する。満たさなければ
			// 「返すクラスタが端である」ことが証明されていない → 退路へ。
			last_verify_info.stage = 3;                    // 判定段
			if (tg_sa && !(Nj[c_lo] == index_type(0))) {
				vcp::throw_error<vcp::verification_error>(
					"spimats::policy_generalized_eigs_with_info: requested smallest eigenvalue but the lowest boundary is not certified as the bottom of the spectrum (no claim)");
			}
			if (tg_la && !(Nj[c_hi + 1] == n)) {
				vcp::throw_error<vcp::verification_error>(
					"spimats::policy_generalized_eigs_with_info: requested largest eigenvalue but the highest boundary is not certified as the top of the spectrum (no claim)");
			}

			// --- 整合検査【対応: veigs.m L458-491 の LDL Result check】---
			// S の各クラスタの certified 番号幅がクラスタの近似固有値個数と一致
			// することを要求(不一致 = 近似が粗い / 見落とし → 無主張)。
			for (std::size_t c = c_lo; c <= c_hi; c++) {
				const std::size_t szc = chi[c] - clo[c] + 1;
				if (!(static_cast<std::size_t>(Nj[c + 1] - Nj[c]) == szc)) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: certified index range does not match cluster size (no claim)");   // 【L485-491】
				}
				if (!have_dn[c] || !have_up[c + 1]) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: missing separator bound (no claim)");
				}
				if (!(dnj[c] <= upj[c + 1])) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: empty separation interval (no claim)");
				}
			}
			}   // try(分離段)
			catch (const vcp::verification_error&) {
				// [SPI-F3 F3-5] Rayleigh 退路【対応: veigs.m L459-485】: sa / la
				// に限り、端の固有値 1 個を次で返す。
				//   sa: [λ_RL, sup(vᵀ[A']v / vᵀ[B']v)]、v = 最小近似固有値の AP
				//       固有ベクトル(点)、λ_RL = RoughLower の certified 下界
				//   la: [inf(vᵀ[A']v / vᵀ[B']v), λ_RU](鏡像)
				// 根拠: B' ≻ 0(入口で certified)のもとで任意の v ≠ 0 の Rayleigh 商
				//   は λ_1 以上 λ_n 以下 [Rayleigh 原理、出典未逐語確認]。商は
				//   vᵀ[A']v・vᵀ[B']v を区間で評価し、分母の下端 > 0 でなければ無主張。
				//   λ_RL / λ_RU は is_spd による certified 境界。
				// 退路も失敗なら元の例外文言で再 throw。
				if (!(tg_sa || tg_la)) throw;
				if (ap.eigenvectors.empty()) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: AP eigenvectors unavailable for Rayleigh fallback (no claim)");
				}
				const std::size_t vi = tg_sa ? ord[0] : ord[m - 1];
				if (vi >= ap.eigenvectors.size()
				    || ap.eigenvectors[vi].size() != static_cast<std::size_t>(n)) throw;
				const std::vector<_T>& v = ap.eigenvectors[vi];
				std::vector<interval_type> Av, Bv;
				vcp::spimats_kernel::mul_im_v(IA, v, Av);    // [A']v
				vcp::spimats_kernel::mul_im_v(B, v, Bv);     // [B']v
				interval_type num(_T(0)), den(_T(0));
				for (std::size_t i2 = 0; i2 < static_cast<std::size_t>(n); i2++) {
					num += v[i2] * Av[i2];
					den += v[i2] * Bv[i2];
				}
				if (!(den.lower() > _T(0))) throw;         // 分母の正値性が確定しない → 元の無主張
				const interval_type rq = num / den;
				_T lam_r;
				if (tg_sa) {
					if (!rough_lower(lam_r)) throw;
					if (!(lam_r <= rq.upper())) throw;
					ray_encl = interval_type(lam_r, rq.upper());
				}
				else {
					if (!rough_upper(lam_r)) throw;
					if (!(rq.lower() <= lam_r)) throw;
					ray_encl = interval_type(rq.lower(), lam_r);
				}
				ray_fb = true;
			}
			if (ray_fb) {
				// 返却は 1 個(k > 1 でも 1 個 --- strict 版は基底の個数検査で
				// state_error になる: 契約どおり)
				result.eigenvalues.clear();
				result.eigenvalues.push_back(ray_encl * restore_scale);   // 元尺度へ(区間乗算)
				result.returned_count = 1;
				result.returned_real_count = 1;
				result.converged_count = 1;
				last_verify_info.veigs_index_lo.assign(1, tg_sa ? std::size_t(1) : static_cast<std::size_t>(n));
				last_verify_info.veigs_index_hi = last_verify_info.veigs_index_lo;
				result.index_lo = last_verify_info.veigs_index_lo;   // [VCP-C1 C1-10] 確定番号(1 個)
				result.index_hi = last_verify_info.veigs_index_hi;
				last_verify_info.veigs_rayleigh_fallback = true;
				last_verify_info.veigs_rank_certified = 1;   // [SPI-F4 F4-1] 端の固有値 1 個 = 順位証明済み
				result.rank_certified = last_verify_info.veigs_rank_certified;   // [VCP-C1 C1-10]
				last_verify_info.stage = 4;                // 成功(退路)
				return result;
			}

			// ===============================================================
			// [SPI-D1 v1.5 実装] U5-b: Lehmann–Behnke 精密化 = eigsymge 委譲
			// 規範出典: 『定理 4.1.3(オーナー提供 2026-07-24)』: X̃, Y, Λ̃ に
			//   対し G := I − YBX̃。‖G‖∞ < 1 ならば YBX̃ は正則で、
			//   r := |Y(AX̃−BX̃Λ̃)|e + ‖Y(AX̃−BX̃Λ̃)‖∞/(1−‖G‖∞)·|G|e、
			//   U_i := {z : |z−λ̃_i| ≤ r_i} の合併が全固有値を含み、単連結
			//   領域が m 個の U_i からなるとき重複度込みで m 個の固有値が
			//   その中に存在(連結成分の重複度計数 --- クラスタ処理の根拠)。
			//   U3 の系 3.2.1(G = RA−I・行和ノルム・乗数合成)と同型構造で
			//   あることに注意。実装は既存 imats::eigsymge [確認済み
			//   vcp/imats.hpp L419-522、Miyajima 2012 Th.3 引用] への委譲
			//   (DD-12。定理 4.1.3 と実装の同一性照合はオーナー --- 設計 §9)。
			// 行マップ: veigs.m L588-634(Lehmann_Behnke__)。veig.m の
			//   はしご探索・err 倍加・index ループは【移植しない】(v1.5 で
			//   破棄 --- eigsymge 委譲では定理 4.1.3 の一発評価に置換。
			//   veig.m は対照参照のみ)。
			// _DP(区間密ポリシー --- DD-11)の用い方: matrix<> ラッパーは
			//   本ヘッダの include 連鎖に無く、ラッパーの friend eigsymge は
			//   ポリシー複製+メンバ委譲の薄い転送 [確認済み vcp/matrix.hpp
			//   L872-880] のため、同一計算になるポリシー直接呼び出しで
			//   委譲する(_DP=pidblas の派生オーバーライドは virtual
			//   mul_im_m/mul_m_im 経由で有効 [確認済み pidblas.hpp L90-142])。
			// 【バックエンド指針(設計 §7.2-9、オーナー裁定 2026-07-25)】
			//   固有値用途の _DP は次を用いること:
			//     double: vcp::pidblas(点側 pdblas = LAPACK dsyev/dsygv)
			//     kv::dd: vcp::imats<kv::dd, vcp::pddblas>
			//     汎用 T(mpfr 含む): vcp::imats<_T, vcp::mats2<_T>>
			//       (mats2 の eigsym/eigsymge は tblas/tlapack の tsyev/tsygv
			//        + LAPACK info 検査つき --- 失敗は正直 throw)
			//   既定 _DP = vcp::imats<_T> の点側(作業ツリーでは mats2、HEAD では
			//   mats)について、旧記述「mats 自前の eigsym は(近)重複固有値
			//   入力で未終了」は現行では再現しない(SPI-V3 dp_hang_repro.cpp:
			//   最小再現 3 件とも正常終了、veigs の既定 _DP でもクラスタ精密化が
			//   停止 --- 2026-10-14。[SPI-F3] 訂正)。上記 3 系の指針は性能・
			//   LAPACK info 検査の観点で引き続き推奨。
			// ===============================================================
			std::vector<_T> flo(C), fhi(C);                // 精密化+交差後の端(S の範囲のみ有効)
			{
				using std::abs;
				// --- V 列の対応(ソート順 → ap 原順の添字置換 ord)は U5-a 側で
				// 構成済み([SPI-F3] 選択集合 S と Rayleigh 退路が共用)。
				// --- V の取得(設計 §7.2-1)【対応: veigs.m L157/L504 の V】---
				// 追加の固有値解なし: スコープ内の ap.eigenvectors を直接使う
				// (keep_ap_eigenvectors はスナップショット保持の別ノブで無関係)。
				// 空・寸法不整合は正直 throw(無主張)。
				if (ap.eigenvectors.size() < m) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: AP eigenvectors unavailable for refinement (no claim)");
				}
				for (std::size_t t = 0; t < m; t++) {
					if (ap.eigenvectors[t].size() != static_cast<std::size_t>(n)) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: AP eigenvector dimension mismatch (no claim)");
					}
				}
				// --- Bmid の LU(全体で 1 回、以後は列ごとの solve ---
				// 設計 §7.2-2)【対応: veigs.m L595 の mid(B)\ ...】---
				_AP Bm3;
				vcp::mid_for_approximation(B, Bm3);
				const vcp::lu_factor_handle<_T, index_type> bh =
					Bm3.policy_lu_factorize_with_info();
				if (!bh.valid()) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: LU factorization of mid(B) failed for refinement (no claim)");
				}
				// --- クラスタ 1 回分の精密化(区画内ローカルラムダ)---
				// 戻り false = S_B 非 PD / eigsymge 不成立 / 近似の死
				// (リカバリ候補があれば次候補で再試行)。
				auto refine_once = [&](std::size_t c, std::size_t kc,
				                       const _T& rho_c, const _T& sigma_c,
				                       interval_type& lo_out, interval_type& hi_out) -> bool {
					// シフト(do_shift = 1 相当)【対応: veigs.m L580-588】
					const _T shift = ev[clo[c]];               // λ_shift = eig_list(r)【L584】
					// 【必要な昇格(片側のみ --- SR-1)】点−点の端点方向保証
					const _T rho_s = (interval_type(rho_c) - shift).upper();     // 【L586】
					const _T sigma_s = (interval_type(sigma_c) - shift).lower(); // 【L587】
					// [A'] = [A] − shift·[B](区間)【L585。pencil_AB 再利用】
					const base_type IAsh = pencil_AB(shift);
					_AP Ash;
					vcp::mid_for_approximation(IAsh, Ash);     // 近似中心(nV 用)
					// --- 列データの構成【対応: L591-600】---
					std::vector<std::vector<_T> > Vc(kc), nV(kc);
					std::vector<std::vector<interval_type> > bv(kc), av(kc), anv(kc), errv(kc);
					for (std::size_t t = 0; t < kc; t++) {
						Vc[t] = ap.eigenvectors[ord[clo[c] + t]];
						// nV = Bmid⁻¹(Amid'·Ṽ)【L595】--- 点近似(保証に
						// 使わない: 残差 err_vec の区間評価が全てを担う)
						const linear_solve_result<_T> sres = bh.solve_with_info(
							Ash.policy_mul_vec(Ash, Vc[t]));
						if (!sres.converged) return false;
						nV[t] = sres.x;
						// 区間評価(R2 混合カーネル --- 疎×点列)
						vcp::spimats_kernel::mul_im_v(B, Vc[t], bv[t]);      // [B]Ṽ
						vcp::spimats_kernel::mul_im_v(IAsh, Vc[t], av[t]);   // [A']Ṽ
						vcp::spimats_kernel::mul_im_v(IAsh, nV[t], anv[t]);  // [A']nV
						vcp::spimats_kernel::mul_im_v(B, nV[t], errv[t]);    // [B]nV
						for (std::size_t i2 = 0; i2 < static_cast<std::size_t>(n); i2++) {
							errv[t][i2] -= av[t][i2];          // err_vec = [B]nV − [A']Ṽ【L596】
						}
					}
					// --- k×k 係数行列【対応: L591-600】---
					// Err = (err_vecᵀ err_vec)/c(c = λ_B_min --- 保守的
					// 残差エネルギー【L597。意図: B ⪰ c·I による粗い支配】)
					// A0 = Ṽᵀ[B]Ṽ【L591】/ A1 = Ṽᵀ[A']Ṽ【L592】/
					// A2 = Ṽᵀ[A']nV − nVᵀ([B]nV − [A']Ṽ) + Err【L600】
					// 点列×区間列の内積は区画内直接ループ(SR-1: kv 混合演算)
					std::vector<std::vector<interval_type> > A0(kc), A1(kc), A2(kc);
					for (std::size_t p = 0; p < kc; p++) {
						A0[p].assign(kc, interval_type(_T(0)));
						A1[p].assign(kc, interval_type(_T(0)));
						A2[p].assign(kc, interval_type(_T(0)));
					}
					for (std::size_t p = 0; p < kc; p++) {
						for (std::size_t q = 0; q < kc; q++) {
							interval_type s0(_T(0)), s1(_T(0)), s2a(_T(0)), s2b(_T(0)), se(_T(0));
							for (std::size_t i2 = 0; i2 < static_cast<std::size_t>(n); i2++) {
								s0 += Vc[p][i2] * bv[q][i2];
								s1 += Vc[p][i2] * av[q][i2];
								s2a += Vc[p][i2] * anv[q][i2];
								s2b += nV[p][i2] * errv[q][i2];
								se += errv[p][i2] * errv[q][i2];
							}
							A0[p][q] = s0;
							A1[p][q] = s1;
							// 【必要な昇格(SR-1 対象外)】点 λ_B_min での除算に
							// 方向保証を持たせる区間化
							A2[p][q] = s2a - s2b + se / interval_type(lambda_B_min);
						}
					}
					// --- 下界側(σ)と上界側(ρ)の 2 面【対応: L602-629】---
					interval_type low_end(_T(0)), up_end(_T(0));
					for (int side = 0; side < 2; side++) {
						// side 0【L603-605】: S_A = −A2 + σ'A1, S_B = −A1 + σ'A0
						// side 1【L617-620】: S_A = A2 − ρ'A1, S_B = A1 − ρ'A0
						std::vector<std::vector<interval_type> > SAe(kc), SBe(kc);
						for (std::size_t p = 0; p < kc; p++) {
							SAe[p].assign(kc, interval_type(_T(0)));
							SBe[p].assign(kc, interval_type(_T(0)));
							for (std::size_t q = 0; q < kc; q++) {
								if (side == 0) {
									SAe[p][q] = -A2[p][q] + sigma_s * A1[p][q];
									SBe[p][q] = -A1[p][q] + sigma_s * A0[p][q];
								}
								else {
									SAe[p][q] = A2[p][q] - rho_s * A1[p][q];
									SBe[p][q] = A1[p][q] - rho_s * A0[p][q];
								}
							}
						}
						// hull 対称化(端点 min/max --- U5-a の sym_hull と同法。
						// 区間内積の演算順差による端点非対称の防御)
						// 【対応: L607/L614/L621-622/L629 の hull(·,·')】
						for (std::size_t p = 0; p < kc; p++) {
							for (std::size_t q = p + 1; q < kc; q++) {
								for (int w = 0; w < 2; w++) {
									std::vector<std::vector<interval_type> >& M2 =
										(w == 0) ? SAe : SBe;
									const _T lo2 = (M2[p][q].lower() < M2[q][p].lower())
										? M2[p][q].lower() : M2[q][p].lower();
									const _T up2 = (M2[p][q].upper() > M2[q][p].upper())
										? M2[p][q].upper() : M2[q][p].upper();
									M2[p][q] = interval_type(lo2, up2);
									M2[q][p] = M2[p][q];
								}
							}
						}
						// [SPI-F3 F3-3 (i)] eigsymge に渡す小行列 S_A / S_B の全成分の両端点が
						// 有限でなければ無主張(return false → 呼び出し側で「Lehmann refinement
						// failed after recovery (no claim)」)。非有限量(大スケール入力の中間
						// オーバーフロー --- novel N2)を保証の根拠にしない。判定は比較演算のみ
						// (!(x − x == 0) が NaN / ±∞ を失敗側へ。T 汎用、double cast なし)。
						for (std::size_t p = 0; p < kc; p++) {
							for (std::size_t q = 0; q < kc; q++) {
								const _T& al = SAe[p][q].lower(); const _T& au = SAe[p][q].upper();
								const _T& bl = SBe[p][q].lower(); const _T& bu = SBe[p][q].upper();
								if (!(al - al == _T(0)) || !(au - au == _T(0))) return false;
								if (!(bl - bl == _T(0)) || !(bu - bu == _T(0))) return false;
							}
						}
						// S_B の正定値性前検査(設計 §7.2-4)【対応: L607/L622 の
						// isspd(hull(SB,SB'))】: k=1 は certified 正の下端、
						// k≥2 は eigsymge を (S_B, I) に適用し全包含の下端 > 0
						// (道具の再利用)。[UNVERIFIED: Lehmann 前提としての
						// S_B ≻ 0 の必要性は定理照合と併せてオーナー]
						if (kc == 1) {
							if (!(SBe[0][0].lower() > _T(0))) return false;
						}
						else {
							_DP SBc, Ik;
							SBc.zeros(static_cast<int>(kc), static_cast<int>(kc));
							Ik.zeros(static_cast<int>(kc), static_cast<int>(kc));
							for (std::size_t p = 0; p < kc; p++) {
								for (std::size_t q = 0; q < kc; q++) {
									SBc.v[p + kc * q] = SBe[p][q];
								}
								Ik.v[p + kc * p] = interval_type(_T(1));
							}
							bool pd = true;
							try {
								SBc.eigsymge(Ik);          // (S_B, I) の全固有値包含
								for (std::size_t p = 0; p < kc; p++) {
									if (!(SBc.v[p + kc * p].lower() > _T(0))) pd = false;
								}
							}
							catch (const vcp::verification_error&) { pd = false; }
							catch (const std::exception&) { pd = false; }   // [SPI-F3 F3-3 (ii)] 非 vcp 例外(例: 区間除算の 0 割)も無主張へ
							if (!pd) return false;
						}
						// 端値の取得
						interval_type endv(_T(0));
						if (kc == 1) {
							// k = 1 の縮退(区間スカラー比 --- オーナー確認)
							// 【対応: L611-612/L625-626】
							endv = SAe[0][0] / SBe[0][0];
						}
						else {
							// eigsymge((S_A, S_B))へ委譲【対応: L614/L629 の
							// veig の置き換え】。返る包含は overlap 併合済み
							// [確認済み imats.hpp L505-527] --- 連結成分の
							// 重複度計数(定理 4.1.3 後段)がクラスタの個数
							// 保証を与える。下界側は最小の下端、上界側は
							// 最大の上端を採る(設計 §7.2-6)。
							_DP SAc, SBc2;
							SAc.zeros(static_cast<int>(kc), static_cast<int>(kc));
							SBc2.zeros(static_cast<int>(kc), static_cast<int>(kc));
							for (std::size_t p = 0; p < kc; p++) {
								for (std::size_t q = 0; q < kc; q++) {
									SAc.v[p + kc * q] = SAe[p][q];
									SBc2.v[p + kc * q] = SBe[p][q];
								}
							}
							try {
								SAc.eigsymge(SBc2);
							}
							catch (const vcp::verification_error&) { return false; }
							catch (const std::exception&) { return false; }   // [SPI-F3 F3-3 (ii)]
							endv = SAc.v[0];
							for (std::size_t p = 1; p < kc; p++) {
								const interval_type& dv2 = SAc.v[p + kc * p];
								if (side == 0) {
									if (dv2.lower() < endv.lower()) endv = dv2;
								}
								else {
									if (dv2.upper() > endv.upper()) endv = dv2;
								}
							}
						}
						// [SPI-F3 F3-3 (iii)] 採用した端値の端点が有限でなければ無主張
						if (!(endv.lower() - endv.lower() == _T(0)) || !(endv.upper() - endv.upper() == _T(0))) return false;
						if (side == 0) low_end = endv; else up_end = endv;
					}
					// 合成【対応: L631-634】: shift を戻す(点+区間の混合加算
					// --- SR-1。端点は呼び出し側で採る)
					lo_out = low_end + shift;
					hi_out = up_end + shift;
					return true;
				};
				// --- クラスタごとの精密化+リカバリ(設計 §7.2-8)。[SPI-F3 F3-7]
				// S のクラスタ c ∈ [c_lo, c_hi] にだけ適用 ---
				for (std::size_t c = c_lo; c <= c_hi; c++) {
					const std::size_t kc = chi[c] - clo[c] + 1;
					// 候補 1(veigs.m 忠実): U5-a の分離値に、端クラスタでは
					// L494-499 の押し出しを適用(端の外側に固有値が無いことは
					// U5-a の Rough 判定で certified --- 押し出しは健全)。
					// 候補 2(リカバリ): 押し出し前の生の分離値(端クラスタ
					// のみ異なる --- 内側境界は押し出せず候補は 1 つ)。
					const _T rho_raw = upj[c];
					const _T sigma_raw = dnj[c + 1];
					_T rho_push = rho_raw;
					_T sigma_push = sigma_raw;
					if (c == 0 && Nj[0] == index_type(0)) {        // 【対応: L497-498】
						const _T cand = ev[clo[0]] - abs(ev[clo[0]]);
						if (cand < rho_push) rho_push = cand;
					}
					if (c == C - 1 && Nj[C] == n) {                // 【対応: L494-495】
						const _T cand = ev[chi[C - 1]] + abs(ev[chi[C - 1]]);
						if (cand > sigma_push) sigma_push = cand;
					}
					interval_type lo_i(_T(0)), hi_i(_T(0));
					bool done = refine_once(c, kc, rho_push, sigma_push, lo_i, hi_i);
					if (!done && !(rho_push == rho_raw && sigma_push == sigma_raw)) {
						// リカバリ: σ/ρ の取り直し(第 2 候補)
						done = refine_once(c, kc, rho_raw, sigma_raw, lo_i, hi_i);
					}
					if (!done) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: Lehmann refinement failed after recovery (no claim)");
					}
					// U5-a の分離区間との交差(設計 §7.2-7): 分離の個数保証を
					// 保持したまま幅だけ締める。交差が空 = 実装矛盾 → 無主張。
					const _T lo_r = lo_i.lower();
					const _T hi_r = hi_i.upper();
					flo[c] = (dnj[c] > lo_r) ? dnj[c] : lo_r;
					fhi[c] = (upj[c + 1] < hi_r) ? upj[c + 1] : hi_r;
					if (!(flo[c] <= fhi[c])) {
						vcp::throw_error<vcp::verification_error>(
							"spimats::policy_generalized_eigs_with_info: refinement/separation intersection is empty (no claim)");
					}
				}
			}

			// ===============================================================
			// [SPI-F4 F4-3 / F4-4 / F4-5] 要求順位の検査(sm / lm / target_*)。
			// 整合検査と Lehmann 精密化の後、出力復元(restore_scale)の前。S の
			// 確定番号を [N_lo+1, N_hi](N_lo = Nj[c_lo]、N_hi = Nj[c_hi+1])、返却
			// 区間の和を [L, U](正規化後の尺度 = flo / fhi)、σ' を正規化後のシフト
			// (sm は 0)とする。未返却の固有値は番号 ≤ N_lo と > N_hi の 2 群。
			// 一方向述語(F4-3。いずれも既存 sep_eval をそのまま呼ぶ --- 新しい
			// 証明器は作らない):
			//   U(x): 右側 sep_eval を λ = x + δ で評価。rc = 0・zero = 0・
			//         err ≤ inf(λ − x) なら「λ_{N+1} > λ − err ≥ x」が certified
			//         (SPI-F2 F2-1 の向き。N = neg)。
			//   L(x): 左側 sep_eval を λ = x − δ で評価。rc = 0・zero = 0・
			//         err ≤ inf(x − λ) なら「λ_N < λ + err ≤ x」が certified。
			//   δ = max(ratio·|x|, ratio·w, meps)(w = 窓幅 d / m)、gap 引数は
			//   δ/(1−ratio)(thr = δ。err > thr は return 2 = 失敗)。zero = 0 は
			//   分離子ループと同じ要求(零慣性ブロックでは N の帰属が決まらない)。
			//   err ≤ δ の代わりに inf(λ − x) と比べるのは、λ = x ± δ の点計算の
			//   丸めを吸収するため【必要な昇格(SR-1 型)】。
			// 条件(F4-4):
			//   sm / target_*: d = sup max(|L − σ'|, |U − σ'|)(返却固有値の距離の上界)。
			//     (N_lo == 0 または L(σ' − d) ≥ N_lo) かつ (N_hi == n または U(σ' + d) ≤ N_hi)。
			//     σ' ∓ d は方向丸め(inf / sup)で取る。
			//   lm: m = min |返却区間の端点|(返却固有値の |λ| の下界。どれかの区間が
			//     0 を含めば証明できず −1)。
			//     (N_lo == 0 または (U(−m) ≤ 0 かつ L(m) ≥ N_lo)) かつ
			//     (N_hi == n または (U(−m) ≤ N_hi かつ L(m) ≥ n))
			//     (U(−m) / L(m) は 1 回だけ評価して両節で共用 --- sep_eval 最大 2 回)。
			// 数学: Sylvester の慣性則で N は全メンバー共通の番号勘定、sep_eval の
			//   err は F2-1 の向きで certified。未返却の固有値が窓の外(sm / target)
			//   または内(lm)に厳密に入れば、返却固有値の距離の上界 d(または下界 m)
			//   と比較して順位が従う [独立導出、出典未逐語確認]。
			// 失敗(sep_eval の return 1 / 2、条件不成立、検査中の例外)は record では
			//   −1 を書いて続行、require では無主張。sa / la は F3-4 で証明済み(1)。
			// 返却区間・個数は off / record で変えない(F4-6)。追加した sep_eval の
			//   回数を veigs_rank_check_factorizations に記録(F4-5、T-6)。
			// ===============================================================
			last_verify_info.veigs_rank_certified = 0;
			last_verify_info.veigs_rank_check_factorizations = 0;
			if (tg_sa || tg_la) {
				last_verify_info.veigs_rank_certified = 1;   // F3-4: 端であることを慣性で証明済み
			}
			else if (spi_opt.target_rank_check != spimats_options::target_rank_check_mode::off) {
				auto rank_check_ = [&]() -> int {
					const index_type N_lo = Nj[c_lo];
					const index_type N_hi = Nj[c_hi + 1];
					const bool tg_lm = (opt.target == vcp::eig_target::largest_magnitude);
					const _T sigma_p = (opt.target == vcp::eig_target::smallest_magnitude) ? _T(0) : mid(opt.shift);
					// 返却区間の和 [L, U] と窓幅 w(sm / target: d、lm: m)
					_T Lb = flo[c_lo], Ub = fhi[c_lo];
					for (std::size_t c = c_lo + 1; c <= c_hi; c++) {
						if (flo[c] < Lb) Lb = flo[c];
						if (fhi[c] > Ub) Ub = fhi[c];
					}
					_T w = _T(0);
					if (tg_lm) {
						for (std::size_t c = c_lo; c <= c_hi; c++) {
							if (!(flo[c] > _T(0)) && !(fhi[c] < _T(0))) return -1;   // 0 ∈ [flo, fhi]
							const _T e0 = abs(flo[c]), e1 = abs(fhi[c]);
							const _T mc = (e0 < e1) ? e0 : e1;
							if (c == c_lo || mc < w) w = mc;
						}
					}
					else {
						const _T d0 = abs(interval_type(Lb) - interval_type(sigma_p)).upper();
						const _T d1 = abs(interval_type(Ub) - interval_type(sigma_p)).upper();
						w = (d0 > d1) ? d0 : d1;
					}
					auto delta_of = [&](const _T& x) -> _T {
						_T dl = ratio * abs(x);
						const _T dw = ratio * w;
						if (dw > dl) dl = dw;
						if (!(dl > meps)) dl = meps;
						return dl;
					};
					// U(x): λ_{N+1} > x を certify する N(成功で true)
					auto pred_U = [&](const _T& x, index_type& N_out) -> bool {
						const _T dl = delta_of(x);
						const _T lam = x + dl;
						const _T slack = (interval_type(lam) - interval_type(x)).lower();
						index_type neg = 0, zero = 0;
						_T err = _T(0);
						last_verify_info.veigs_rank_check_factorizations++;
						const int rc = sep_eval(lam, dl / one_m_ratio, neg, zero, err, false);
						if (rc != 0 || zero != index_type(0) || !(err <= slack)) return false;
						N_out = neg;
						return true;
					};
					// L(x): λ_N < x を certify する N
					auto pred_L = [&](const _T& x, index_type& N_out) -> bool {
						const _T dl = delta_of(x);
						const _T lam = x - dl;
						const _T slack = (interval_type(x) - interval_type(lam)).lower();
						index_type neg = 0, zero = 0;
						_T err = _T(0);
						last_verify_info.veigs_rank_check_factorizations++;
						const int rc = sep_eval(lam, dl / one_m_ratio, neg, zero, err, true);
						if (rc != 0 || zero != index_type(0) || !(err <= slack)) return false;
						N_out = neg;
						return true;
					};
					if (tg_lm) {
						bool haveU = false, haveL = false;
						index_type nU = 0, nL = 0;
						if (N_lo != index_type(0)) {
							if (!pred_U(-w, nU)) return -1;
							haveU = true;
							if (!(nU == index_type(0))) return -1;
							if (!pred_L(w, nL)) return -1;
							haveL = true;
							if (!(nL >= N_lo)) return -1;
						}
						if (N_hi != n) {
							if (!haveU && !pred_U(-w, nU)) return -1;
							if (!(nU <= N_hi)) return -1;
							if (!haveL && !pred_L(w, nL)) return -1;
							if (!(nL >= n)) return -1;
						}
					}
					else {
						if (N_lo != index_type(0)) {
							const _T x = (interval_type(sigma_p) - interval_type(w)).lower();
							index_type N = 0;
							if (!pred_L(x, N)) return -1;
							if (!(N >= N_lo)) return -1;
						}
						if (N_hi != n) {
							const _T x = (interval_type(sigma_p) + interval_type(w)).upper();
							index_type N = 0;
							if (!pred_U(x, N)) return -1;
							if (!(N <= N_hi)) return -1;
						}
					}
					return 1;
				};
				int certified = -1;
				try { certified = rank_check_(); }
				catch (const vcp::verification_error&) { certified = -1; }
				catch (const std::exception&) { certified = -1; }   // 検査中の例外も「証明できず」(返却は変えない)
				last_verify_info.veigs_rank_certified = certified;
				if (certified != 1
				    && spi_opt.target_rank_check == spimats_options::target_rank_check_mode::require) {
					vcp::throw_error<vcp::verification_error>(
						"spimats::policy_generalized_eigs_with_info: requested rank could not be certified (no claim)");
				}
			}

			// --- 出力(U5-a 分離区間 ∩ U5-b 精密化区間 --- 設計 §7.2-7)---
			// クラスタ c の全メンバー(番号 N_c+1..N_{c+1})に交差後区間
			// [flo_c, fhi_c] を返す。並びは昇順(ソート済み ev と同順)。
			// [SPI-F3 F3-7] 返すのは S のクラスタ(昇順)。各クラスタの確定番号
			// [N_c+1, N_{c+1}](1 始まり)を veigs_index_lo / hi に露出する(F3-4)。
			result.eigenvalues.clear();
			result.eigenvalues.reserve(m);
			last_verify_info.veigs_index_lo.clear();
			last_verify_info.veigs_index_hi.clear();
			std::size_t n_out = 0;
			for (std::size_t c = c_lo; c <= c_hi; c++) {
				// [SPI-F3 F3-6] 元尺度へ: λ = μ·2^(a−b)(区間乗算 --- 包含)
				const interval_type encl = interval_type(flo[c], fhi[c]) * restore_scale;
				for (std::size_t t = clo[c]; t <= chi[c]; t++) {
					result.eigenvalues.push_back(encl);
					result.index_lo.push_back(static_cast<std::size_t>(Nj[c]) + 1);   // [VCP-C1 C1-10] 各メンバーに確定番号
					result.index_hi.push_back(static_cast<std::size_t>(Nj[c + 1]));
					n_out++;
				}
				last_verify_info.veigs_index_lo.push_back(static_cast<std::size_t>(Nj[c]) + 1);
				last_verify_info.veigs_index_hi.push_back(static_cast<std::size_t>(Nj[c + 1]));
			}
			result.returned_count = n_out;                 // 「返した包含の個数」= |S| で正直に
			result.returned_real_count = n_out;
			result.converged_count = n_out;
			result.rank_certified = last_verify_info.veigs_rank_certified;   // [VCP-C1 C1-10] F4 の判定をそのまま
			last_verify_info.stage = 4;                    // 成功
			return result;
		}
		// >>> END [SPI-U5] <<<
		return result;   // not reached
	}
};

// alias template(設計 §2.1): DP のみ差し替えたい利用者向け。
//   vcp::spimats_dp<double, vcp::pdblas> == spimats<double, spmats<double>, pdblas>
template <typename _T, class _DP>
using spimats_dp = spimats<_T, vcp::spmats<_T>, _DP>;

} // namespace vcp

#endif // VCP_SPIMATS_HPP
