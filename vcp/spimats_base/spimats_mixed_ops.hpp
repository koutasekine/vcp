// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spimats_mixed_ops.hpp
// SPI-1: 疎行列の区間型ポリシー(spimats)の混合演算カーネル 7 種(SPI-R2)。
// namespace vcp::spimats_kernel の自由関数テンプレート:
//   mul_im_m : 区間疎 × 点疎        -> 区間疎(シフト構成等)
//   mul_m_im : 点疎 × 区間疎        -> 区間疎(R·(A x~ − b) 等)
//   add_im_m : 区間疎 + 点疎        -> 区間疎(A − cB 等は点側を −c 倍して加算)
//   add_m_im : 点疎 + 区間疎        -> 区間疎(可換性により add_im_m へ委譲)
//   sub_im_m : 区間疎 − 点疎        -> 区間疎(点側の厳密符号反転 + add_im_m 委譲。SPI-K1)
//   mul_im_v : 区間疎行列 × 点ベクトル   -> 区間ベクトル(残差 A x~)
//   mul_m_iv : 点疎行列 × 区間ベクトル   -> 区間ベクトル(R の適用)
//   mul_v_im : 点横ベクトル × 区間疎行列 -> 区間横ベクトル(RA − I の行ごと構成)
//
// 型規約(設計書 §6): 値型 2 つ(区間側スカラー _T・点側スカラー _TP)で
//   テンプレート化。累積はすべて kv::interval<_T>。点側の行列・ベクトルを
//   区間型へ変換したコピーは生成しない(スカラーレベルは kv の
//   interval<_T> op _TP 自動昇格に委ねる。パターン合流で点値を区間の器に
//   格納する箇所のみ kv::interval<_T>(b) の退化区間構成 --- 端点直接代入で
//   丸めを含まない)。
// 丸め規律(設計書 G5): FP 環境の丸めモード切替(C の丸め切替関数・VCP の
//   丸めガード)や自前の方向丸め算術は一切書かない。丸めは kv の区間演算に
//   委譲する。
// 疎パターン規則: 入力は finalize 済み前提(未 finalize は入口で finalize)。
//   非格納要素 = 厳密 0。出力で厳密 [0,0] になった要素(区間演算の厳密相殺)
//   は格納しない(spmats invariant)。
// 初版は逐次実装(OpenMP 並列化は測ってから --- 設計書 §9 と同方針)。

#ifndef VCP_SPIMATS_MIXED_OPS_HPP
#define VCP_SPIMATS_MIXED_OPS_HPP

#include <algorithm>
#include <cstddef>
#include <vector>

#include <kv/interval.hpp>

#include <vcp/error.hpp>
#include <vcp/spmats.hpp>
#include <vcp/spimats_base/spimats_convert.hpp>   // spimats_convert_detail::is_strict_zero_interval

namespace vcp {
namespace spimats_kernel {

// >>> REVIEW-REQUIRED [SPI-R2: 混合演算 7 種のカーネル] <<<
// STATUS: UNREVIEWED
// CLAIM: 各カーネルの出力は、入力区間の任意の点実現(点側は退化区間)に
//   対する真の積・和を要素ごとに包含する。根拠: 出力の各要素は kv の区間
//   演算(interval<_T> op _TP の自動昇格・外側丸め)のみの合成で構成され、
//   本区画自身は端点算術・丸め操作を一切行わない。厳密 [0,0] となった
//   出力要素の非格納は値意味論(非格納 = 厳密 0)と等価。
// REDUCES-TO: kv::interval の四則(外側丸め)と自動昇格
//   interval<_T> op _TP / 退化区間コンストラクタ interval<_T>(b)(厳密)。
// SELF-ARITHMETIC: なし
// TESTS: sandbox/tests/spimats_mixed_ops_test.cpp
// ---------------------------------------------------------------------------

namespace detail {

	// 行フラッシュ(mul_im_m / mul_m_im の共有部・50 行制約による分割):
	// 触れた列をソートして厳密 [0,0] 以外を IC へ格納し、acc / touched を
	// 初期状態([0,0] / 0)へ戻す。事後条件: acc 全要素 [0,0]、touched 全 0。
	template <typename _T, typename _Index>
	inline void flush_accumulated_row_(spmats<kv::interval<_T>, _Index>& IC,
	                                   const _Index i, std::vector<_Index>& cols,
	                                   std::vector<kv::interval<_T> >& acc,
	                                   std::vector<char>& touched)
	{
		std::sort(cols.begin(), cols.end());
		for (std::size_t c = 0; c < cols.size(); c++) {
			const std::size_t j = static_cast<std::size_t>(cols[c]);
			if (!spimats_convert_detail::is_strict_zero_interval(acc[j])) {
				IC.add(i, cols[c], acc[j]);
			}
			acc[j] = kv::interval<_T>(_T(0));
			touched[j] = 0;
		}
		cols.clear();
	}

} // namespace detail

// mul_im_m: IC = IA * B(区間疎 (m×k) × 点疎 (k×n) -> 区間疎 (m×n))
template <typename _T, typename _TP, typename _Index>
void mul_im_m(const spmats<kv::interval<_T>, _Index>& IA,
              const spmats<_TP, _Index>& B,
              spmats<kv::interval<_T>, _Index>& IC)
{
	if (IA.columnsize() != B.rowsize()) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::mul_im_m: dimension mismatch");
	}
	if (!IA.is_finalized()) IA.finalize();
	if (!B.is_finalized()) B.finalize();
	const spmats<kv::interval<_T>, _Index> Ac = IA.as_csr();
	const spmats<_TP, _Index> Bc = B.as_csr();
	const _Index m = Ac.rowsize();
	const _Index n = Bc.columnsize();
	const std::vector<_Index>& ao = Ac.outer_index();
	const std::vector<_Index>& ai = Ac.inner_index();
	const std::vector<kv::interval<_T> >& av = Ac.values();
	const std::vector<_Index>& bo = Bc.outer_index();
	const std::vector<_Index>& bi = Bc.inner_index();
	const std::vector<_TP>& bv = Bc.values();
	IC.clear();
	IC.resize(m, n);
	std::vector<kv::interval<_T> > acc(static_cast<std::size_t>(n), kv::interval<_T>(_T(0)));
	std::vector<char> touched(static_cast<std::size_t>(n), 0);
	std::vector<_Index> cols;
	for (_Index i = 0; i < m; i++) {
		// LOOP INVARIANT(p ループ各周回開始時): acc[j] は行 i の部分積
		//   Σ_{既処理の p} IA(i, ai[p]) * B(ai[p], j) を保持し、touched[j] = 1
		//   ⟺ j ∈ cols ⟺ acc[j] に加算が行われた。cols 外の acc は [0,0]。
		for (_Index p = ao[static_cast<std::size_t>(i)]; p < ao[static_cast<std::size_t>(i) + 1]; p++) {
			const kv::interval<_T>& a = av[static_cast<std::size_t>(p)];
			const _Index k = ai[static_cast<std::size_t>(p)];
			for (_Index q = bo[static_cast<std::size_t>(k)]; q < bo[static_cast<std::size_t>(k) + 1]; q++) {
				const _Index j = bi[static_cast<std::size_t>(q)];
				acc[static_cast<std::size_t>(j)] += a * bv[static_cast<std::size_t>(q)];
				if (!touched[static_cast<std::size_t>(j)]) {
					touched[static_cast<std::size_t>(j)] = 1;
					cols.push_back(j);
				}
			}
		}
		detail::flush_accumulated_row_(IC, i, cols, acc, touched);
	}
	IC.finalize();
}

// mul_m_im: IC = B * IA(点疎 (m×k) × 区間疎 (k×n) -> 区間疎 (m×n))
template <typename _T, typename _TP, typename _Index>
void mul_m_im(const spmats<_TP, _Index>& B,
              const spmats<kv::interval<_T>, _Index>& IA,
              spmats<kv::interval<_T>, _Index>& IC)
{
	if (B.columnsize() != IA.rowsize()) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::mul_m_im: dimension mismatch");
	}
	if (!B.is_finalized()) B.finalize();
	if (!IA.is_finalized()) IA.finalize();
	const spmats<_TP, _Index> Bc = B.as_csr();
	const spmats<kv::interval<_T>, _Index> Ac = IA.as_csr();
	const _Index m = Bc.rowsize();
	const _Index n = Ac.columnsize();
	const std::vector<_Index>& bo = Bc.outer_index();
	const std::vector<_Index>& bi = Bc.inner_index();
	const std::vector<_TP>& bv = Bc.values();
	const std::vector<_Index>& ao = Ac.outer_index();
	const std::vector<_Index>& ai = Ac.inner_index();
	const std::vector<kv::interval<_T> >& av = Ac.values();
	IC.clear();
	IC.resize(m, n);
	std::vector<kv::interval<_T> > acc(static_cast<std::size_t>(n), kv::interval<_T>(_T(0)));
	std::vector<char> touched(static_cast<std::size_t>(n), 0);
	std::vector<_Index> cols;
	for (_Index i = 0; i < m; i++) {
		// LOOP INVARIANT(p ループ各周回開始時): acc[j] は行 i の部分積
		//   Σ_{既処理の p} B(i, bi[p]) * IA(bi[p], j) を保持し、touched[j] = 1
		//   ⟺ j ∈ cols。点値 bv は昇格のみで型変換コピーを作らない。
		for (_Index p = bo[static_cast<std::size_t>(i)]; p < bo[static_cast<std::size_t>(i) + 1]; p++) {
			const _TP& b = bv[static_cast<std::size_t>(p)];
			const _Index k = bi[static_cast<std::size_t>(p)];
			for (_Index q = ao[static_cast<std::size_t>(k)]; q < ao[static_cast<std::size_t>(k) + 1]; q++) {
				const _Index j = ai[static_cast<std::size_t>(q)];
				acc[static_cast<std::size_t>(j)] += av[static_cast<std::size_t>(q)] * b;
				if (!touched[static_cast<std::size_t>(j)]) {
					touched[static_cast<std::size_t>(j)] = 1;
					cols.push_back(j);
				}
			}
		}
		detail::flush_accumulated_row_(IC, i, cols, acc, touched);
	}
	IC.finalize();
}

// add_im_m: IC = IA + B(区間疎 + 点疎、同寸法。パターンは合併)
// A − cB 等の差は、呼び出し側で点側を −c 倍してから加算する。
template <typename _T, typename _TP, typename _Index>
void add_im_m(const spmats<kv::interval<_T>, _Index>& IA,
              const spmats<_TP, _Index>& B,
              spmats<kv::interval<_T>, _Index>& IC)
{
	if (IA.rowsize() != B.rowsize() || IA.columnsize() != B.columnsize()) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::add_im_m: dimension mismatch");
	}
	if (!IA.is_finalized()) IA.finalize();
	if (!B.is_finalized()) B.finalize();
	const spmats<kv::interval<_T>, _Index> Ac = IA.as_csr();
	const spmats<_TP, _Index> Bc = B.as_csr();
	const std::vector<_Index>& ao = Ac.outer_index();
	const std::vector<_Index>& ai = Ac.inner_index();
	const std::vector<kv::interval<_T> >& av = Ac.values();
	const std::vector<_Index>& bo = Bc.outer_index();
	const std::vector<_Index>& bi = Bc.inner_index();
	const std::vector<_TP>& bv = Bc.values();
	IC.clear();
	IC.resize(Ac.rowsize(), Ac.columnsize());
	for (_Index i = 0; i < Ac.rowsize(); i++) {
		_Index p = ao[static_cast<std::size_t>(i)], pe = ao[static_cast<std::size_t>(i) + 1];
		_Index q = bo[static_cast<std::size_t>(i)], qe = bo[static_cast<std::size_t>(i) + 1];
		// LOOP INVARIANT: 各周回開始時、行 i の出力済み列はすべて
		//   min(ai[p], bi[q])(未処理の最小列)より小さい。両 CSR は列ソート
		//   済み・重複なしなので 2 ポインタ合流で各列は高々 1 回出力される。
		while (p < pe || q < qe) {
			kv::interval<_T> y;
			_Index j;
			if (q >= qe || (p < pe && ai[static_cast<std::size_t>(p)] < bi[static_cast<std::size_t>(q)])) {
				j = ai[static_cast<std::size_t>(p)];
				y = av[static_cast<std::size_t>(p)]; p++;              // IA のみ: 値コピー(演算なし)
			}
			else if (p >= pe || bi[static_cast<std::size_t>(q)] < ai[static_cast<std::size_t>(p)]) {
				j = bi[static_cast<std::size_t>(q)];
				y = kv::interval<_T>(bv[static_cast<std::size_t>(q)]); q++;   // B のみ: 退化区間(厳密)
			}
			else {
				j = ai[static_cast<std::size_t>(p)];
				y = av[static_cast<std::size_t>(p)] + bv[static_cast<std::size_t>(q)]; p++; q++;   // 両方: 区間+点(昇格)
			}
			if (!spimats_convert_detail::is_strict_zero_interval(y)) {
				IC.add(i, j, y);   // 厳密相殺 [0,0] のみ非格納
			}
		}
	}
	IC.finalize();
}

// add_m_im: IC = B + IA。点と区間の要素ごとの和は可換(kv の interval+点 /
// 点+interval は同一の外側丸め結果)なので add_im_m へ委譲する(レビュー面積の
// 削減。ループ不変条件は add_im_m 側に記載)。
template <typename _T, typename _TP, typename _Index>
void add_m_im(const spmats<_TP, _Index>& B,
              const spmats<kv::interval<_T>, _Index>& IA,
              spmats<kv::interval<_T>, _Index>& IC)
{
	add_im_m(IA, B, IC);
}

// mul_im_v: iy = IA * x(区間疎 (m×n) × 点ベクトル (n) -> 区間ベクトル (m)、密出力)
template <typename _T, typename _TP, typename _Index>
void mul_im_v(const spmats<kv::interval<_T>, _Index>& IA,
              const std::vector<_TP>& x,
              std::vector<kv::interval<_T> >& iy)
{
	if (x.size() != static_cast<std::size_t>(IA.columnsize())) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::mul_im_v: dimension mismatch");
	}
	if (!IA.is_finalized()) IA.finalize();
	const spmats<kv::interval<_T>, _Index> Ac = IA.as_csr();
	const std::vector<_Index>& ao = Ac.outer_index();
	const std::vector<_Index>& ai = Ac.inner_index();
	const std::vector<kv::interval<_T> >& av = Ac.values();
	iy.assign(static_cast<std::size_t>(Ac.rowsize()), kv::interval<_T>(_T(0)));
	for (_Index i = 0; i < Ac.rowsize(); i++) {
		// LOOP INVARIANT(p ループ各周回開始時): iy[i] は行 i の部分和
		//   Σ_{既処理の p} IA(i, ai[p]) * x[ai[p]] を保持する(非格納要素の
		//   寄与は厳密 0)。
		for (_Index p = ao[static_cast<std::size_t>(i)]; p < ao[static_cast<std::size_t>(i) + 1]; p++) {
			iy[static_cast<std::size_t>(i)] +=
				av[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(ai[static_cast<std::size_t>(p)])];
		}
	}
}

// mul_m_iv: iy = B * ix(点疎 (m×n) × 区間ベクトル (n) -> 区間ベクトル (m)、密出力)
template <typename _T, typename _TP, typename _Index>
void mul_m_iv(const spmats<_TP, _Index>& B,
              const std::vector<kv::interval<_T> >& ix,
              std::vector<kv::interval<_T> >& iy)
{
	if (ix.size() != static_cast<std::size_t>(B.columnsize())) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::mul_m_iv: dimension mismatch");
	}
	if (!B.is_finalized()) B.finalize();
	const spmats<_TP, _Index> Bc = B.as_csr();
	const std::vector<_Index>& bo = Bc.outer_index();
	const std::vector<_Index>& bi = Bc.inner_index();
	const std::vector<_TP>& bv = Bc.values();
	iy.assign(static_cast<std::size_t>(Bc.rowsize()), kv::interval<_T>(_T(0)));
	for (_Index i = 0; i < Bc.rowsize(); i++) {
		// LOOP INVARIANT(p ループ各周回開始時): iy[i] は行 i の部分和
		//   Σ_{既処理の p} B(i, bi[p]) * ix[bi[p]] を保持する。点値は昇格のみ。
		for (_Index p = bo[static_cast<std::size_t>(i)]; p < bo[static_cast<std::size_t>(i) + 1]; p++) {
			iy[static_cast<std::size_t>(i)] +=
				ix[static_cast<std::size_t>(bi[static_cast<std::size_t>(p)])] * bv[static_cast<std::size_t>(p)];
		}
	}
}

// mul_v_im: iy^T = x^T * IA(点横ベクトル (m) × 区間疎 (m×n) -> 区間横ベクトル (n)、
// 密出力)。RA − I の行ごと構成(x = R の行)用。
template <typename _T, typename _TP, typename _Index>
void mul_v_im(const std::vector<_TP>& x,
              const spmats<kv::interval<_T>, _Index>& IA,
              std::vector<kv::interval<_T> >& iy)
{
	if (x.size() != static_cast<std::size_t>(IA.rowsize())) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::mul_v_im: dimension mismatch");
	}
	if (!IA.is_finalized()) IA.finalize();
	const spmats<kv::interval<_T>, _Index> Ac = IA.as_csr();
	const std::vector<_Index>& ao = Ac.outer_index();
	const std::vector<_Index>& ai = Ac.inner_index();
	const std::vector<kv::interval<_T> >& av = Ac.values();
	iy.assign(static_cast<std::size_t>(Ac.columnsize()), kv::interval<_T>(_T(0)));
	for (_Index i = 0; i < Ac.rowsize(); i++) {
		// LOOP INVARIANT(i ループ各周回開始時): iy[j] は部分和
		//   Σ_{i' < i} x[i'] * IA(i', j) を保持する(CSR 行走査による
		//   散布加算。各 (i', j) 格納要素はちょうど 1 回加算される)。
		for (_Index p = ao[static_cast<std::size_t>(i)]; p < ao[static_cast<std::size_t>(i) + 1]; p++) {
			iy[static_cast<std::size_t>(ai[static_cast<std::size_t>(p)])] +=
				x[static_cast<std::size_t>(i)] * av[static_cast<std::size_t>(p)];
		}
	}
}

// >>> END [SPI-R2] <<<

// >>> REVIEW-REQUIRED [SPI-R7: sub_im_m(区間疎 − 点疎)の薄い委譲(SPI-K1 K1-b)] <<<
// STATUS: REVIEWED-OK (Kouta Sekine, 2026-07-26)
// CLAIM: sub_im_m(IA, B, IC) の出力は IA − B の要素ごとの真の差を包含する。
//   根拠: −B の構成は値の単項符号反転のみ(IEEE 浮動小数点・kv::dd・
//   kv::mpfr のいずれも符号反転は厳密演算・パターン不変)で丸めを含まず、
//   包含は委譲先 add_im_m(SPI-R2)の CLAIM に帰着する(add_m_im の先例と
//   同じレビュー面積削減方式)。
// REDUCES-TO: add_im_m(SPI-R2)+ 単項符号反転の厳密性
// SELF-ARITHMETIC: なし(符号反転のみ・丸めなし)
// TESTS: sandbox/tests/spimats_k1_test.cpp
// ---------------------------------------------------------------------------

// sub_im_m: IC = IA − B(区間疎 − 点疎、同寸法。パターンは合併)
// 点側を厳密に符号反転した −B を構成して add_im_m へ委譲する。
// ループ不変条件・2 ポインタ合流規則は add_im_m 側に記載。
template <typename _T, typename _TP, typename _Index>
void sub_im_m(const spmats<kv::interval<_T>, _Index>& IA,
              const spmats<_TP, _Index>& B,
              spmats<kv::interval<_T>, _Index>& IC)
{
	if (IA.rowsize() != B.rowsize() || IA.columnsize() != B.columnsize()) {
		vcp::throw_error<vcp::dimension_error>("spimats_kernel::sub_im_m: dimension mismatch");
	}
	if (!B.is_finalized()) B.finalize();
	// Bn = −B: as_csr() のコピー(非破壊性のため必須の 1 本)に対し、
	// 非 const values()(SPI-R9)で値のみを 1 パス符号反転する。
	// 符号反転は IEEE / kv::dd / kv::mpfr のいずれも厳密演算・パターン不変
	// なので、outer/inner・finalized 状態・CSR 形式はそのまま保たれる
	// (COO 再構築・finalize 不要 — 旧実装の廃止理由。SPI-R8-b)。
	spmats<_TP, _Index> Bn = B.as_csr();
	std::vector<_TP>& v = Bn.values();
	for (std::size_t k = 0; k < v.size(); k++) {
		v[k] = -v[k];   // 単項符号反転(厳密・丸めなし)
	}
	add_im_m(IA, Bn, IC);
}

// >>> END [SPI-R7] <<<

// 将来の追記方法(設計書 §6): mul_iv_m(区間ベクトル × 点疎)等の追加カーネルは
// 本ファイルに REVIEW-REQUIRED 区画(SPI-R6 以降の新番号・STATUS は未レビュー
// 状態の定型ヘッダ)を新設して追加する。既存区画への追記は行わない
// (spimats_convert.hpp 末尾と同運用)。

} // namespace spimats_kernel
} // namespace vcp

#endif // VCP_SPIMATS_MIXED_OPS_HPP
