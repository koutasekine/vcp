// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spimats_convert.hpp
// SPI-1: 疎行列の区間型ポリシー(spimats)の型変換ヘルパ群。
// namespace vcp 直下の自由関数テンプレート:
//   SPI-R1: mid_for_approximation(行列・ベクトル)/ midrad_split(行列・ベクトル)
//   SPI-R3: point_to_interval(行列・ベクトル)/ identity_interval(厳密 B=I)
//   SPI-R5: convert_interval_matrix(区間型間の疎外側丸め変換)
//
// 丸め規律(設計書 spimats_design_v0.md G5):
//   本ヘッダは FP 環境の丸めモード切替(C の丸め切替関数・VCP の丸め
//   ガード)や自前の方向丸め算術を一切使用しない。丸めに関わる操作は
//   すべて kv(mid / midrad / 区間演算)と vcp::convert に委譲する。
//
// 疎パターン規則(実装指示書 §3 Phase 1):
//   - 入力は finalize 済み CSR 前提。未 finalize の入力は入口で finalize
//     する(finalize() は論理 const)。
//   - 非格納要素 = 厳密 0(区間では [0,0])。finalize 済み格納列に厳密 0 は
//     存在しない(spmats invariant: to_csr が落とし assign_* が拒否する)。
//     重複も finalize 済みなので存在しない。
//   - 区間の zero は [0,0] のみ(端点スカラー比較で判定。kv の区間 == は
//     certainly-equal であり同一性判定に使わない ---
//     設計書 §8 の比較規約。「厳密に点 0 か」の x == _T(0) のみ例外)。
//
// 依存: kv/interval.hpp は利用者側で include 済みであること(丸め backend
//   rdouble/rdd/rmpfr の include も kv の慣例どおり利用者責務)。
//   SPI-R5 の型ペアは vcp_converter.hpp の該当オーバーロードが可視である
//   必要がある(kv/rdouble.hpp + kv/rdd.hpp 等を include してから本ヘッダを
//   include する)。

#ifndef VCP_SPIMATS_CONVERT_HPP
#define VCP_SPIMATS_CONVERT_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

#include <kv/interval.hpp>

#include <vcp/error.hpp>
#include <vcp/spmats.hpp>
#include <vcp/vcp_converter.hpp>

namespace vcp {

namespace spimats_convert_detail {

	// 区間の厳密 zero 判定([0,0] のみ true)。spmatrix::is_strict_zero_
	// (spmatrix.hpp SPC-1)の区間オーバーロードと同一の端点スカラー比較。
	// kv の区間 == (certainly-equal) を同一性判定に使わない規約(設計書 §8)
	// に従い、端点で書く。NaN 端点は false 側に倒れる(格納維持 = 安全側)。
	template <typename _Tv>
	inline bool is_strict_zero_interval(const kv::interval<_Tv>& x) {
		return x.lower() == _Tv(0) && x.upper() == _Tv(0);
	}

} // namespace spimats_convert_detail

// >>> REVIEW-REQUIRED [SPI-R1: midrad_split / mid_for_approximation(行列・ベクトル)] <<<
// STATUS: UNREVIEWED
// CLAIM: midrad_split(IA, M, R) は全格納要素で kv::midrad(x, m, r) を呼び、
//   要素ごとに ∀a ∈ [IA]_ij : |a − M_ij| ≤ R_ij を満たす点行列対 (M, R) を
//   構成する(非格納要素は IA=[0,0]・M=0・R=0 で自明に成立。m または r が
//   厳密 0 の要素は該当行列へ格納しない --- 数学的値は 0 で同じ)。
//   mid_for_approximation は近似入力専用(非厳密)で、いかなる保証の主張の
//   一部にも使ってはならない。
// REDUCES-TO: kv::midrad(kv/interval.hpp L622: r は方向丸め sub_up の max で
//   真の中点まわり半径の上界)/ kv::mid(近似のみ)。本区画自身は丸めを行わない。
// SELF-ARITHMETIC: なし
// TESTS: sandbox/tests/spimats_convert_test.cpp
// ---------------------------------------------------------------------------

// mid_for_approximation(行列版): spmats<kv::interval<_T>> -> 任意の
// spmats<_T,_Index> 派生 _APM(AP = spumar 等の派生 add 意味論を保つため
// テンプレート)。要素 = kv::mid()。
//
// 契約: 近似入力専用(非厳密)。mid() は最近接丸めの近似中点であり、
//   精度保証の主張のいかなる部分にも使ってはならない。保証用の (M, R)
//   分解は midrad_split を使うこと(mid()+rad() の別呼び構成は禁止 ---
//   設計書 §9: rad は真の中点まわりの半径であり、計算された m の丸めずれ
//   分だけ包含が壊れる)。
// 格納規則: mid が厳密 0 になる格納要素(例: [-1,1])は M に格納しない
//   (spmats invariant「explicit zero は格納しない」と整合。数学的には
//   M_ij = 0 を意味し、近似入力としてはそのままでよい)。
template <typename _T, typename _Index, class _APM>
void mid_for_approximation(const spmats<kv::interval<_T>, _Index>& IA, _APM& M)
{
	static_assert(std::is_base_of<spmats<_T, _Index>, _APM>::value,
		"vcp::mid_for_approximation: M must derive from spmats<_T, _Index>");
	if (!IA.is_finalized()) IA.finalize();
	const spmats<kv::interval<_T>, _Index> C = IA.as_csr();
	const std::vector<_Index>& outer = C.outer_index();
	const std::vector<_Index>& inner = C.inner_index();
	const std::vector<kv::interval<_T> >& val = C.values();
	M.clear();
	M.resize(C.rowsize(), C.columnsize());
	for (_Index i = 0; i < C.rowsize(); i++) {
		for (_Index p = outer[static_cast<std::size_t>(i)];
		     p < outer[static_cast<std::size_t>(i) + 1]; p++) {
			const _T m = mid(val[static_cast<std::size_t>(p)]);
			if (!(m == _T(0))) {
				M.add(i, inner[static_cast<std::size_t>(p)], m);
			}
		}
	}
	M.finalize();
}

// mid_for_approximation(ベクトル版): 右辺 b 用。契約は行列版と同じ
// (近似入力専用・保証への使用禁止)。
template <typename _T>
void mid_for_approximation(const std::vector<kv::interval<_T> >& iv,
                           std::vector<_T>& v)
{
	v.resize(iv.size());
	for (std::size_t k = 0; k < iv.size(); k++) {
		v[k] = mid(iv[k]);
	}
}

// midrad_split(行列版): 保証用 mid/rad 分解(厳密)。全格納要素で
// kv::midrad() を使用する。
//   CLAIM: ∀A′∈[IA]: |A′ − M| ≤ R(要素ごと)。
// M / R の格納パターンは IA のパターンの部分集合で、互いに異なりうる
// (m == 0 の要素は M に、r == 0 の要素は R に格納されない。いずれも
// 数学的値 0 と等価で CLAIM は保たれる: r == 0 なら x は退化区間で
// A′ = m 厳密、m == 0 なら |A′ − 0| ≤ r)。
template <typename _T, typename _Index>
void midrad_split(const spmats<kv::interval<_T>, _Index>& IA,
                  spmats<_T, _Index>& M, spmats<_T, _Index>& R)
{
	if (!IA.is_finalized()) IA.finalize();
	const spmats<kv::interval<_T>, _Index> C = IA.as_csr();
	const std::vector<_Index>& outer = C.outer_index();
	const std::vector<_Index>& inner = C.inner_index();
	const std::vector<kv::interval<_T> >& val = C.values();
	M.clear();
	M.resize(C.rowsize(), C.columnsize());
	R.clear();
	R.resize(C.rowsize(), C.columnsize());
	for (_Index i = 0; i < C.rowsize(); i++) {
		for (_Index p = outer[static_cast<std::size_t>(i)];
		     p < outer[static_cast<std::size_t>(i) + 1]; p++) {
			_T m, r;
			midrad(val[static_cast<std::size_t>(p)], m, r);
			const _Index j = inner[static_cast<std::size_t>(p)];
			if (!(m == _T(0))) M.add(i, j, m);
			if (!(r == _T(0))) R.add(i, j, r);
		}
	}
	M.finalize();
	R.finalize();
}

// midrad_split(ベクトル版): 右辺 ib 用。密ベクトルなので零落としは行わず
// 全要素を書き込む。CLAIM は行列版と同じ(∀b′∈[ib]: |b′ − m| ≤ r)。
template <typename _T>
void midrad_split(const std::vector<kv::interval<_T> >& iv,
                  std::vector<_T>& m, std::vector<_T>& r)
{
	m.resize(iv.size());
	r.resize(iv.size());
	for (std::size_t k = 0; k < iv.size(); k++) {
		midrad(iv[k], m[k], r[k]);
	}
}

// >>> END [SPI-R1] <<<

// >>> REVIEW-REQUIRED [SPI-R3: 点→区間化(行列・ベクトル)・厳密 B=I 構成] <<<
// STATUS: UNREVIEWED
// CLAIM: point_to_interval は各格納要素 a を退化区間 [a, a] に厳密に写す
//   (kv::interval<_T>(a) は同型 _T からの端点直接構成で丸めを含まない)。
//   したがって [IA] = {A} であり A ∈ [IA]。identity_interval は対角 [1,1]
//   (厳密表現)の n×n 単位行列を構成する([IB] = {I})。
// REDUCES-TO: kv::interval<_T> の同型スカラーコンストラクタ(端点代入のみ)。
//   本区画自身は丸めを行わない。
// SELF-ARITHMETIC: なし
// TESTS: sandbox/tests/spimats_convert_test.cpp
// ---------------------------------------------------------------------------

// point_to_interval(行列版): 点疎行列(AP 派生は基底参照で受ける)->
// spmats<kv::interval<_T>> の退化区間化。格納パターンは厳密に保存される
// (格納値は spmats invariant により非零、退化区間 [a,a] (a≠0) は [0,0] に
// なり得ない)。
template <typename _T, typename _Index>
void point_to_interval(const spmats<_T, _Index>& A,
                       spmats<kv::interval<_T>, _Index>& IA)
{
	if (!A.is_finalized()) A.finalize();
	const spmats<_T, _Index> C = A.as_csr();
	const std::vector<_Index>& outer = C.outer_index();
	const std::vector<_Index>& inner = C.inner_index();
	const std::vector<_T>& val = C.values();
	IA.clear();
	IA.resize(C.rowsize(), C.columnsize());
	for (_Index i = 0; i < C.rowsize(); i++) {
		for (_Index p = outer[static_cast<std::size_t>(i)];
		     p < outer[static_cast<std::size_t>(i) + 1]; p++) {
			IA.add(i, inner[static_cast<std::size_t>(p)],
			       kv::interval<_T>(val[static_cast<std::size_t>(p)]));
		}
	}
	IA.finalize();
}

// point_to_interval(ベクトル版): 退化区間化。
template <typename _T>
void point_to_interval(const std::vector<_T>& v,
                       std::vector<kv::interval<_T> >& iv)
{
	iv.resize(v.size());
	for (std::size_t k = 0; k < v.size(); k++) {
		iv[k] = kv::interval<_T>(v[k]);
	}
}

// identity_interval: policy_eigs_with_info_impl(#2、B=I 転送層)用の厳密
// 単位行列構成。対角 [1,1] は全ての対象スカラー型で厳密表現可能。
// assign_csr により born-finalized CSR で構成する。
template <typename _T, typename _Index>
void identity_interval(const _Index n, spmats<kv::interval<_T>, _Index>& IB)
{
	if (n < 0) {
		vcp::throw_error<vcp::invalid_argument>(
			"vcp::identity_interval: size must be nonnegative");
	}
	const std::size_t un = static_cast<std::size_t>(n);
	std::vector<_Index> outer(un + 1);
	std::vector<_Index> inner(un);
	std::vector<kv::interval<_T> > val(un, kv::interval<_T>(_T(1)));
	for (std::size_t k = 0; k < un; k++) {
		outer[k] = static_cast<_Index>(k);
		inner[k] = static_cast<_Index>(k);
	}
	outer[un] = n;
	IB.clear();
	IB.assign_csr(n, n, outer, inner, val);
}

// >>> END [SPI-R3] <<<

// >>> REVIEW-REQUIRED [SPI-R5: 区間型間の疎外側丸め変換] <<<
// STATUS: UNREVIEWED
// CLAIM: convert_interval_matrix(IA, IB) は要素ごとに vcp::convert の
//   区間型間オーバーロード(vcp_converter.hpp L181-225: 下端 rnd=-1・
//   上端 rnd=+1 の外側丸め)を適用し、[IA]_src ⊆ [IB]_dst を満たす
//   (縮小変換でも包含が成立)。走査は spmatrix::convert_from_
//   (spmatrix.hpp SPC-1)と同じ単一走査・フォーマット保存(finalize 済み
//   CSC は CSC、それ以外は CSR)・strict-zero 規則(変換後 [0,0] のみ
//   落とす)。外側丸めで非零区間が [0,0] に潰れることはない(両端厳密 0 は
//   元が [0,0] のときのみ --- finalize 済み格納列には存在しない)ため
//   零落としは防御であり、パターンは実質保存される。
// REDUCES-TO: vcp::convert の区間型間オーバーロード(外側丸め)。
//   本区画自身は丸めを行わない。
// SELF-ARITHMETIC: なし
// TESTS: sandbox/tests/spimats_convert_test.cpp
// ---------------------------------------------------------------------------

// convert_interval_matrix: spmats<kv::interval<_TS>> ->
// spmats<kv::interval<_TD>>。同型ペア(_TS == _TD)は vcp::convert の
// identity オーバーロード(vcp_converter.hpp L60)で単純コピーになる。
// 型ペアに対応する vcp::convert オーバーロードの可視性は利用者の include
// 構成に依存する(ヘッダ冒頭の依存の項を参照)。
template <typename _TS, typename _TD, typename _Index>
void convert_interval_matrix(const spmats<kv::interval<_TS>, _Index>& IA,
                             spmats<kv::interval<_TD>, _Index>& IB)
{
	if (!IA.is_finalized()) IA.finalize();
	// フォーマット保存(spmatrix::convert_from_ と同一): finalize 済み CSC は
	// CSC のまま、それ以外(CSR / COO)は CSR。
	const bool use_csc = IA.is_finalized() && IA.format() == vcp::sparse_csc;
	const spmats<kv::interval<_TS>, _Index> src = use_csc ? IA.as_csc() : IA.as_csr();
	const _Index rows = src.rowsize();
	const _Index cols = src.columnsize();
	const std::size_t nouter = static_cast<std::size_t>(use_csc ? cols : rows);
	const std::vector<_Index>& src_outer = src.outer_index();
	const std::vector<_Index>& src_inner = src.inner_index();
	const std::vector<kv::interval<_TS> >& src_value = src.values();
	std::vector<_Index> outer(nouter + 1);
	std::vector<_Index> inner(src_value.size());
	std::vector<kv::interval<_TD> > val(src_value.size());
	// 単一走査: 各格納値を vcp::convert(外側丸め)し、変換後に厳密 [0,0] に
	// なった値のみ落とす(spmats invariant の防御。上記 CLAIM のとおり実際には
	// 発生しない)。outer は書き込みカーソル out で再構築する。
	std::size_t out = 0;
	outer[0] = 0;
	for (std::size_t i = 0; i < nouter; i++) {
		for (_Index k = src_outer[i]; k < src_outer[i + 1]; k++) {
			kv::interval<_TD> y;
			vcp::convert(src_value[static_cast<std::size_t>(k)], y);
			if (!spimats_convert_detail::is_strict_zero_interval(y)) {
				inner[out] = src_inner[static_cast<std::size_t>(k)];
				val[out] = y;
				out++;
			}
		}
		outer[i + 1] = static_cast<_Index>(out);
	}
	inner.resize(out);
	val.resize(out);
	IB.clear();
	if (use_csc) {
		IB.assign_csc(rows, cols, outer, inner, val);
	}
	else {
		IB.assign_csr(rows, cols, outer, inner, val);
	}
}

// >>> END [SPI-R5] <<<

// 将来の追記方法: 新たな変換(例: 区間ベクトルの型間外側丸め変換)は、
// 本ファイルに REVIEW-REQUIRED 区画(SPI-R6 以降の新番号・STATUS は未レビュー
// 状態の定型ヘッダ)を新設して追加する。既存区画への追記はレビュー済み
// STATUS を無効化するため行わない(オーナーレビュー後に区画単位で STATUS が
// 更新される運用)。

} // namespace vcp

#endif // VCP_SPIMATS_CONVERT_HPP
