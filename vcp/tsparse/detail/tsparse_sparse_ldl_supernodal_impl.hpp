// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLDL-SP SP-1 / SP-2 -- left-looking supernodal LDL^T numeric kernel.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER the
// sparse_ldl_* type skeleton, the symbolic phase and the numeric workspace are
// in scope.  It has no "namespace vcp { }" wrapper; it is injected by
// tsparse_sparse_ldl.hpp.  It includes NOTHING itself -- <algorithm>,
// <chrono>, <type_traits>, <vector> and tblas.hpp are included at file scope
// by tsparse_sparse_ldl.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_ldl.hpp>
//
// ---------------------------------------------------------------------------
// Panel structure (design v1 SS3.2)
//
// A panel owns the columns [pa, pb) (width w) and a below-panel row list
// Rows[s] (height h).  It is a dense (w + h) x w column-major block:
//
//        panel row r < w        <->  matrix row/column pa + r  (diagonal block)
//        panel row w + t        <->  matrix row Rows[s][t]     (off-diagonal)
//
// Left-looking over supernodes: before eliminating s, every finished
// supernode d whose row cursor points into [pa, pb) contributes
//
//        W       = L_d(I, :) * D_d                 (I = rows of R_d in [pa,pb))
//        C(I, I) -= W * L_d(I, :)^T                (symmetric: tgemmtr or tgemm)
//        C(J, I) -= W * L_d(J, :)^T                (J = rows of R_d >= pb, tgemm)
//
// and the panel is then eliminated column by column, RIGHT-looking inside the
// panel, so every panel column is fully updated at the moment it becomes the
// pivot column -- which is exactly what the BK tests of SP-2 need.
//
// ---------------------------------------------------------------------------
// Byte-identity with baseline x none (design v1 SS6-2)
//
// The scalar mirror evaluates, for every element and in the same order as the
// baseline kernel, the single expression
//
//        C(i,c) -= L(i,j) * ( L(c,j) * d_j )
//
// with j ascending over all descendants and then over the panel's own earlier
// columns.  The baseline static path stores exactly W = L*d as its update
// coefficient, so the two systems perform the same operations in the same
// order and agree bit for bit.  The BLOCKED path deliberately does not: it
// accumulates one descendant's contribution in a temporary and subtracts it in
// one go, a different (equally valid) summation order.  That is why the
// blocked path is gated to std::is_floating_point scalars, for which the
// contract asks for pivot-sequence identity and a residual band, not bytes.
//
// ---------------------------------------------------------------------------
// SP-2: dynamic pivoting inside the panel (design v1 SS3.3, D-1, D-4)
//
// The BK test is NOT weakened: lambda is the maximum over the WHOLE column
// (the panel holds every structurally nonzero row of it), and sigma is the
// full column of the candidate row r.  Three outcomes:
//
//   * the pivot needs no exchange, or the exchange partner is a panel COLUMN
//     -> executed.  Rows of already-finished supernodes are RELABELLED (the
//     baseline label-map method, localized to the panel): no data moves and,
//     because relabelling only ever renames labels INSIDE [pa,pb) -- which are
//     columns, never members of any R_d of a not-yet-applied update -- the
//     static structure of every later supernode stays valid.
//   * a 2x2 pivot whose partner column falls one past the panel -> the
//     supernode is SPLIT at that column and the column is merged into the next
//     supernode's panel (D-4).  The split is metadata-only for the finished
//     part (the storage layout is unchanged by construction).
//   * anything else (the partner is not a panel column, or the merged panel
//     could not be formed structurally) -> status pivot_out_of_panel, no
//     partial factor is returned.  Never a silent fallback (D-1).
//
// No OpenMP pragma appears in this file: the sparse layer stays serial and all
// parallelism is left to the existing tblas guards (design v1 SS1, non-goal 2).

#ifndef VCP_TSPARSE_SPARSE_LDL_SUPERNODAL_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LDL_SUPERNODAL_IMPL_HPP

namespace sparse_ldl_detail {

// Compile-time type gate (design v1 SS3.2).  Non-floating scalars (kv::dd,
// kv::mpfr<N>, kv::interval<TT>, std::complex<TT>) always take the scalar
// mirror.  The gate is a COMPILE-TIME tag, not a runtime flag, so the dense
// block kernels are not even instantiated for those types.
template <class T>
struct sparse_ldl_blockable {
    static const bool value = std::is_floating_point<T>::value;
};

// tblas takes int extents; the gate refuses anything that would not fit.
inline bool sparse_ldl_fits_int_(long long v) {
    return v >= 0 && v <= 2147483647LL;
}

// ---------------------------------------------------------------------------
// One descendant update, blocked:  U := L_d(I u J, :) * W^T  (design v1 SS3.2)
//   gemmtr variant: the symmetric (I x I) part is computed as a triangle only
//                   (half the flops, D-3), the (J x I) part by tgemm;
//   gemm variant   : one tgemm over the whole (|I|+|J|) x |I| block.
// Selected by tag: the false_type overload exists only so that non-floating
// scalars never instantiate tblas at all.
// ---------------------------------------------------------------------------
template <class T>
inline int sparse_ldl_block_update_(
    std::true_type, const bool use_gemmtr,
    const int nall, const int ni, const int nj, const int wd,
    const T* Ld, const int md, const T* W, T* U)
{
    if (use_gemmtr) {
        vcp::tgemmtr<T>('L', 'N', 'T', ni, wd,
                        T(1), Ld, md, W, ni, T(0), U, nall);
        if (nj > 0) {
            vcp::tgemm<T>('N', 'T', nj, ni, wd,
                          T(1), Ld + ni, md, W, ni, T(0), U + ni, nall);
            return 2;
        }
        return 1;
    }
    vcp::tgemm<T>('N', 'T', nall, ni, wd,
                  T(1), Ld, md, W, ni, T(0), U, nall);
    return 1;
}

template <class T>
inline int sparse_ldl_block_update_(
    std::false_type, const bool, const int, const int, const int, const int,
    const T*, const int, const T*, T*)
{
    return 0;   // unreachable: the runtime gate below is false for these T
}

// ---------------------------------------------------------------------------
// Symmetric exchange of two PANEL COLUMNS gp < gq (global labels, both inside
// the current panel), design v1 SS3.3.
//
// The panel itself is swapped physically (positional row identity, exactly
// like the frozen dense kernel); the rows of ALREADY FINISHED supernodes are
// RELABELLED through rowloc -- the baseline label-map method, localized to the
// panel -- so no committed data moves.  Relabelling only ever renames labels
// inside the panel's column range, and those labels are never members of a
// row list that a pending update still has to consume, which is why the static
// structure of every later supernode survives an exchange unchanged.
// ---------------------------------------------------------------------------
template <class T, class Index>
inline void sparse_ldl_panel_exchange_(
    const Index gp, const Index gq, const Index pa, const std::size_t m,
    T* P, std::vector<Index>& perm,
    std::vector<Index>& cur_of_asm, std::vector<Index>& asm_of_cur,
    std::vector<std::vector<std::pair<Index, Index> > >& rowloc,
    std::vector<std::vector<Index> >& Rows)
{
    if (gp == gq) return;
    using std::swap;
    const std::size_t lp = static_cast<std::size_t>(gp - pa);
    const std::size_t lq = static_cast<std::size_t>(gq - pa);
    // already-eliminated panel columns: swap the two rows
    for (std::size_t c = 0; c < lp; ++c) swap(P[lp + c * m], P[lq + c * m]);
    // uneliminated symmetric part (lower storage)
    for (std::size_t i = lp + 1u; i < lq; ++i) swap(P[i + lp * m], P[lq + i * m]);
    for (std::size_t i = lq + 1u; i < m; ++i)  swap(P[i + lp * m], P[i + lq * m]);
    swap(P[lp + lp * m], P[lq + lq * m]);
    // permutation record and label maps
    std::swap(perm[static_cast<std::size_t>(gp)], perm[static_cast<std::size_t>(gq)]);
    const std::size_t ap = static_cast<std::size_t>(asm_of_cur[static_cast<std::size_t>(gp)]);
    const std::size_t aq = static_cast<std::size_t>(asm_of_cur[static_cast<std::size_t>(gq)]);
    std::swap(cur_of_asm[ap], cur_of_asm[aq]);
    std::swap(asm_of_cur[static_cast<std::size_t>(gp)], asm_of_cur[static_cast<std::size_t>(gq)]);
    // relabel the rows of finished supernodes
    std::vector<std::pair<Index, Index> >& pl = rowloc[static_cast<std::size_t>(gp)];
    std::vector<std::pair<Index, Index> >& ql = rowloc[static_cast<std::size_t>(gq)];
    for (std::size_t t = 0; t < pl.size(); ++t) {
        Rows[static_cast<std::size_t>(pl[t].first)][static_cast<std::size_t>(pl[t].second)] = gq;
    }
    for (std::size_t t = 0; t < ql.size(); ++t) {
        Rows[static_cast<std::size_t>(ql[t].first)][static_cast<std::size_t>(ql[t].second)] = gp;
    }
    pl.swap(ql);
}

} // namespace sparse_ldl_detail

// ---------------------------------------------------------------------------
// sparse_ldl_supernodal_factorize
//
//   sym  : symbolic result at level full.  perm0 seeds res.perm; the BK
//          exchanges of SP-2 are composed into it in place, so the result is
//          perm = P0 o P_BK (new->old, design v2 SS5.4), exactly as in the
//          baseline kernels.
//   res  : same output contract as the baseline kernels.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_ldl_supernodal_factorize(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const sparse_ldl_symbolic_result<Index>& sym,
    const sparse_ldl_options<T>& opt,
    sparse_ldl_numeric_workspace<T, Index>& ws,
    sparse_ldl_result<T, Index>& res)
{
    static_assert(std::is_signed<Index>::value, "sparse LDL Index must be signed");
    using std::abs;
    using std::sqrt;
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    const std::size_t un = static_cast<std::size_t>(n);
    const Index nsup = sym.n_supernodes;
    const R alpha = (R(1) + sqrt(R(17))) / R(8);   // BK constant, built in R (P2)

    res.D_diag.assign(un, T(0));
    res.D_sub.assign(un, T(0));
    res.D_block2.assign(un, char(0));
    res.n_pivots_1x1 = Index(0);
    res.n_pivots_2x2 = Index(0);
    res.first_zero_pivot = Index(-1);
    res.inconclusive_at = Index(-1);
    res.out_of_panel_at = Index(-1);
    res.n_zero_skips = Index(0);
    res.n_boundary_splits = Index(0);
    res.nnz_L = Index(0);
    res.L_col_ptr.clear();
    res.L_row_ind.clear();
    res.L_val.clear();
    res.perm = sym.perm0;
    res.n_supernodes = nsup;
    res.max_supernode_width = sym.max_supernode_width;
    res.nnz_L_static = sym.nnz_L_static;
    res.dense_delegated = false;

    // Resolved kernel choice; diag_kernel_used starts at "scalar" and is
    // upgraded by the first blocked call that actually happens, so it always
    // reports what RAN, not what was requested.
    const sparse_ldl_diag_kernel dk =
        sparse_ldl_detail::resolve_auto_diag_kernel(opt.diag_kernel);
    const bool use_gemmtr = (dk == sparse_ldl_diag_kernel::gemmtr);
    res.diag_kernel_used = sparse_ldl_kernel_used::scalar;
    res.gemm_call_count = 0;
    res.gemm_time_ns = 0;

    const bool static_mode = (opt.pivoting == sparse_ldl_pivoting::none);

    if (n == Index(0)) {
        res.L_col_ptr.assign(1u, Index(0));
        res.status = sparse_ldl_status::success;
        return;
    }

    // ---- full symmetric adjacency in ASSEMBLY labels (= perm0-permuted
    // labels), exactly like the baseline sparse kernel: after a BK exchange a
    // stored lower entry can become an upper one, so both directions are kept.
    std::vector<Index> adj_ptr(un + 1u, Index(0));
    {
        std::vector<Index> deg(un, Index(0));
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r >= static_cast<Index>(c)) {
                    const std::size_t a = static_cast<std::size_t>(sym.pinv0[static_cast<std::size_t>(r)]);
                    const std::size_t b = static_cast<std::size_t>(sym.pinv0[c]);
                    ++deg[b];
                    if (a != b) ++deg[a];
                }
            }
        }
        for (std::size_t i = 0; i < un; ++i) adj_ptr[i + 1u] = adj_ptr[i] + deg[i];
    }
    std::vector<Index> adj_ind(static_cast<std::size_t>(adj_ptr[un]));
    std::vector<T>     adj_val(static_cast<std::size_t>(adj_ptr[un]), T(0));
    {
        std::vector<Index> head_(adj_ptr.begin(), adj_ptr.end() - 1);
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r < static_cast<Index>(c)) continue;
                const std::size_t a = static_cast<std::size_t>(sym.pinv0[static_cast<std::size_t>(r)]);
                const std::size_t b = static_cast<std::size_t>(sym.pinv0[c]);
                const T& v = val[static_cast<std::size_t>(k)];
                std::size_t pos = static_cast<std::size_t>(head_[b]++);
                adj_ind[pos] = static_cast<Index>(a);
                adj_val[pos] = v;
                if (a != b) {
                    pos = static_cast<std::size_t>(head_[a]++);
                    adj_ind[pos] = static_cast<Index>(b);
                    adj_val[pos] = v;
                }
            }
        }
    }

    // ---- growth reference (D-5): max |a_ii| of the input diagonal, taken
    // from values that are scanned anyway (no extra pass, no extra memory).
    R growth_max_a = R(0), growth_max_d = R(0);
    bool growth_max_a_valid = false;
    for (std::size_t i = 0; i < un; ++i) {
        for (Index k = adj_ptr[i]; k < adj_ptr[i + 1u]; ++k) {
            if (adj_ind[static_cast<std::size_t>(k)] != static_cast<Index>(i)) continue;
            const R a = abs(adj_val[static_cast<std::size_t>(k)]);
            if (a > growth_max_a) { growth_max_a = a; growth_max_a_valid = true; }
            else if (a <= growth_max_a) { growth_max_a_valid = true; }
            else { /* undecidable magnitude: the diagnostic stays invalid */ }
        }
    }

    // ---- per-supernode state.  The column ranges and row lists are DYNAMIC:
    // a boundary split (D-4) moves one column from a supernode to its
    // successor, and BK exchanges relabel rows.
    const std::size_t nsu = static_cast<std::size_t>(nsup);
    std::vector<Index> sn_c0(nsu), sn_c1(nsu);
    std::vector<std::vector<Index> > Rows(nsu);
    std::vector<std::vector<T> >     Lblk(nsu);
    std::vector<Index> col_sn(un, Index(0));
    for (std::size_t s = 0; s < nsu; ++s) {
        sn_c0[s] = sym.supernode_ptr[s];
        sn_c1[s] = sym.supernode_ptr[s + 1u];
        Rows[s].assign(sym.sn_row_ind.begin() + static_cast<std::ptrdiff_t>(sym.sn_row_ptr[s]),
                       sym.sn_row_ind.begin() + static_cast<std::ptrdiff_t>(sym.sn_row_ptr[s + 1u]));
        for (Index j = sn_c0[s]; j < sn_c1[s]; ++j) {
            col_sn[static_cast<std::size_t>(j)] = static_cast<Index>(s);
        }
    }

    // ---- descendant lists: link[d] chains the supernodes whose next
    // unconsumed row falls inside the columns of head[]'s supernode.
    std::vector<Index> head(nsu, Index(-1));
    std::vector<Index> link(nsu, Index(-1));
    std::vector<Index> rptr(nsu, Index(0));      // cursor into Rows[d]

    // ---- label maps (current label <-> assembly label).  Identity in static
    // mode; permuted by BK exchanges.
    std::vector<Index> cur_of_asm(un), asm_of_cur(un);
    for (std::size_t i = 0; i < un; ++i) {
        cur_of_asm[i] = static_cast<Index>(i);
        asm_of_cur[i] = static_cast<Index>(i);
    }

    // ---- row -> (supernode, position in Rows) links, used to relabel the
    // rows of FINISHED supernodes on a BK exchange.  Built only for the
    // dynamic mode (the static mode performs no exchange at all).
    std::vector<std::vector<std::pair<Index, Index> > > rowloc;
    if (!static_mode) rowloc.resize(un);

    // ---- scratch
    ws.relind.assign(un, Index(-1));             // matrix row -> panel row
    std::vector<Index>& relind = ws.relind;
    std::vector<T>& U = ws.update;               // dense update block
    std::vector<T>& W = ws.panel_a;              // W = L_d(I,:) * D_d
    std::vector<Index> desc;                     // descendants of the panel
    std::vector<Index> carry_rows;               // carried column of a split
    std::vector<T>     carry_vals;
    bool carried = false;

    bool any_zero_pivot = false;
    bool stopped_inconclusive = false;
    bool stopped_out_of_panel = false;

    for (Index s = Index(0); s < nsup; ++s) {
        if (stopped_inconclusive || stopped_out_of_panel) break;
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index pa = sn_c0[ss];
        const Index pb = sn_c1[ss];
        if (pb <= pa) { continue; }              // fully absorbed by a split
        const std::size_t w = static_cast<std::size_t>(pb - pa);
        // COPY, not a reference: a boundary split replaces Rows[ss] while the
        // panel still needs its original row list (carry-out, relind release).
        const std::vector<Index> Rs = Rows[ss];
        const std::size_t h = Rs.size();
        const std::size_t m = w + h;

        Lblk[ss].assign(m * w, T(0));
        T* const P = &Lblk[ss][0];

        // ---- relative index map
        for (std::size_t r = 0; r < w; ++r) {
            relind[static_cast<std::size_t>(pa) + r] = static_cast<Index>(r);
        }
        for (std::size_t t = 0; t < h; ++t) {
            relind[static_cast<std::size_t>(Rs[t])] = static_cast<Index>(w + t);
        }

        // ---- gather A into the panel through the label maps.  Column pa is
        // skipped when it was carried over by a boundary split: its values are
        // already reduced and are scattered below instead.
        const std::size_t gather_from = carried ? 1u : 0u;
        for (std::size_t q = gather_from; q < w; ++q) {
            const Index j = pa + static_cast<Index>(q);
            const std::size_t alab = static_cast<std::size_t>(asm_of_cur[static_cast<std::size_t>(j)]);
            for (Index k = adj_ptr[alab]; k < adj_ptr[alab + 1u]; ++k) {
                const Index i = cur_of_asm[static_cast<std::size_t>(adj_ind[static_cast<std::size_t>(k)])];
                if (i < j) continue;             // upper triangle in current labels
                const Index pr = relind[static_cast<std::size_t>(i)];
                if (pr < Index(0)) {             // outside the static pattern
                    res.status = sparse_ldl_status::internal_error;
                    return;
                }
                P[static_cast<std::size_t>(pr) + q * m] = adj_val[static_cast<std::size_t>(k)];
            }
        }
        if (carried) {
            for (std::size_t t = 0; t < carry_rows.size(); ++t) {
                const Index pr = relind[static_cast<std::size_t>(carry_rows[t])];
                if (pr < Index(0)) { res.status = sparse_ldl_status::internal_error; return; }
                P[static_cast<std::size_t>(pr) + 0u * m] = carry_vals[t];
            }
            carried = false;
        }

        // ---- descendants, in ASCENDING order (the order the baseline applies
        // their updates in -- required for the byte-identity contract).
        desc.clear();
        for (Index d = head[ss]; d != Index(-1); ) {
            const Index nx = link[static_cast<std::size_t>(d)];
            desc.push_back(d);
            d = nx;
        }
        head[ss] = Index(-1);
        std::sort(desc.begin(), desc.end());     // integer keys (P6-safe)

        for (std::size_t di = 0; di < desc.size(); ++di) {
            const std::size_t dd = static_cast<std::size_t>(desc[di]);
            const Index d_c0 = sn_c0[dd];
            const std::size_t wd = static_cast<std::size_t>(sn_c1[dd] - d_c0);
            const std::vector<Index>& Rd = Rows[dd];
            const std::size_t d_h = Rd.size();
            const std::size_t md = wd + d_h;

            // I = rows of R_d inside [pa,pb) -- a contiguous run at the cursor,
            // J = the rest (all >= pb because R_d is ascending).
            const std::size_t p0 = static_cast<std::size_t>(rptr[dd]);
            std::size_t p1 = p0;
            while (p1 < d_h && Rd[p1] < pb) ++p1;
            const std::size_t ni = p1 - p0;
            const std::size_t nj = d_h - p1;
            const std::size_t nall = ni + nj;
            if (ni == 0u) { res.status = sparse_ldl_status::internal_error; return; }

            const T* const Ld = &Lblk[dd][0] + wd + p0;      // ld = md

            // Every updated row must be a row of THIS panel.  That is a
            // theorem about the static structure (L(i,j) != 0 and L(c,j) != 0
            // with j < c <= i imply L(i,c) != 0), verified once per descendant
            // in O(|I|+|J|) -- negligible against the O(nall*ni*wd) update it
            // guards, and it turns a structural defect into internal_error
            // instead of a wild write.
            for (std::size_t p = 0; p < nall; ++p) {
                if (relind[static_cast<std::size_t>(Rd[p0 + p])] < Index(0)) {
                    res.status = sparse_ldl_status::internal_error;
                    return;
                }
            }

            // W = L_d(I,:) * D_d   (1x1 and 2x2 blocks of D_d)
            W.assign(ni * wd, T(0));
            for (std::size_t jj = 0; jj < wd; ) {
                const std::size_t gj = static_cast<std::size_t>(d_c0) + jj;
                if (res.D_block2[gj] != char(0) && jj + 1u < wd) {
                    const T d1 = res.D_diag[gj];
                    const T e  = res.D_sub[gj];
                    const T d2 = res.D_diag[gj + 1u];
                    for (std::size_t q = 0; q < ni; ++q) {
                        const T l1 = Ld[q + jj * md];
                        const T l2 = Ld[q + (jj + 1u) * md];
                        W[q + jj * ni]        = l1 * d1 + l2 * e;
                        W[q + (jj + 1u) * ni] = l1 * e + l2 * d2;
                    }
                    jj += 2u;
                } else {
                    const T dj = res.D_diag[gj];
                    for (std::size_t q = 0; q < ni; ++q) {
                        W[q + jj * ni] = Ld[q + jj * md] * dj;
                    }
                    jj += 1u;
                }
            }

            const bool blocked =
                sparse_ldl_detail::sparse_ldl_blockable<T>::value &&
                (static_cast<long long>(wd) >= static_cast<long long>(opt.ldl_min_block_size)) &&
                sparse_ldl_detail::sparse_ldl_fits_int_(static_cast<long long>(md)) &&
                sparse_ldl_detail::sparse_ldl_fits_int_(static_cast<long long>(nall)) &&
                sparse_ldl_detail::sparse_ldl_fits_int_(static_cast<long long>(wd));

            if (!blocked) {
                // ---- scalar mirror: the baseline expression in the baseline
                // order (j ascending outermost).
                for (std::size_t jj = 0; jj < wd; ++jj) {
                    for (std::size_t q = 0; q < ni; ++q) {
                        const T wq = W[q + jj * ni];
                        const std::size_t cc = static_cast<std::size_t>(
                            relind[static_cast<std::size_t>(Rd[p0 + q])]);
                        for (std::size_t p = q; p < nall; ++p) {
                            const std::size_t pr = static_cast<std::size_t>(
                                relind[static_cast<std::size_t>(Rd[p0 + p])]);
                            P[pr + cc * m] -= Ld[p + jj * md] * wq;
                        }
                    }
                }
            } else {
                // ---- blocked: U = L_d(I u J, :) * W^T, then scatter-subtract.
                const std::chrono::steady_clock::time_point t_begin =
                    std::chrono::steady_clock::now();
                U.assign(nall * ni, T(0));
                const int calls = sparse_ldl_detail::sparse_ldl_block_update_<T>(
                    std::integral_constant<bool,
                        sparse_ldl_detail::sparse_ldl_blockable<T>::value>(),
                    use_gemmtr,
                    static_cast<int>(nall), static_cast<int>(ni),
                    static_cast<int>(nj), static_cast<int>(wd),
                    Ld, static_cast<int>(md), &W[0], &U[0]);
                res.diag_kernel_used = use_gemmtr ? sparse_ldl_kernel_used::gemmtr
                                                  : sparse_ldl_kernel_used::gemm;
                res.gemm_call_count += calls;
                for (std::size_t q = 0; q < ni; ++q) {
                    const std::size_t cc = static_cast<std::size_t>(
                        relind[static_cast<std::size_t>(Rd[p0 + q])]);
                    for (std::size_t p = q; p < nall; ++p) {
                        const std::size_t pr = static_cast<std::size_t>(
                            relind[static_cast<std::size_t>(Rd[p0 + p])]);
                        P[pr + cc * m] -= U[p + q * nall];
                    }
                }
                res.gemm_time_ns += static_cast<long long>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - t_begin).count());
            }

            // advance the cursor and relink d to its next target supernode
            rptr[dd] = static_cast<Index>(p1);
            if (p1 < d_h) {
                const Index next_sn = col_sn[static_cast<std::size_t>(Rd[p1])];
                link[dd] = head[static_cast<std::size_t>(next_sn)];
                head[static_cast<std::size_t>(next_sn)] = static_cast<Index>(dd);
            } else {
                link[dd] = Index(-1);
            }
        }

        // =================== panel elimination =========================
        // global row of panel position p
#define VCP_SLDL_GROW(p) (((p) < w) ? (pa + static_cast<Index>(p)) : Rs[(p) - w])

        std::size_t q = 0;
        bool split_here = false;
        bool want_split = false;          // BK needs the next supernode merged in
        Index split_partner = Index(-1);  // the BK candidate row that asked for it
        while (q < w) {
            const Index j = pa + static_cast<Index>(q);

            // ------------------------------------------------ static mode
            if (static_mode) {
                const T d = P[q + q * m];
                const R absd = abs(d);
                if (absd > opt.zero_pivot_tol) {
                    res.D_diag[static_cast<std::size_t>(j)] = d;
                    if (absd > growth_max_d) growth_max_d = absd;
                    P[q + q * m] = T(1);
                    for (std::size_t p = q + 1u; p < m; ++p) P[p + q * m] = P[p + q * m] / d;
                    for (std::size_t cc = q + 1u; cc < w; ++cc) {
                        const T wq = P[cc + q * m] * d;
                        for (std::size_t p = cc; p < m; ++p) {
                            P[p + cc * m] -= P[p + q * m] * wq;
                        }
                    }
                    ++res.n_pivots_1x1;
                } else if (absd <= opt.zero_pivot_tol) {
                    P[q + q * m] = T(1);
                    for (std::size_t p = q + 1u; p < m; ++p) P[p + q * m] = T(0);
                    if (!any_zero_pivot) {
                        res.first_zero_pivot = j;
                        any_zero_pivot = true;
                    }
                    ++res.n_zero_skips;
                } else {
                    res.inconclusive_at = j;
                    stopped_inconclusive = true;
                    break;
                }
                ++q;
                continue;
            }

            // ------------------------------------------------ dynamic (BK)
            // lambda = max_{i > j} |A(i,j)| over the WHOLE column: every
            // structurally nonzero row of column j is a panel row, so the
            // test is not weakened by the panel restriction.
            R lam = R(0);
            std::size_t rpos = 0;
            bool r_valid = false;
            bool search_inconclusive = false;
            for (std::size_t p = q + 1u; p < m; ++p) {
                const R a = abs(P[p + q * m]);
                if (a > lam) { lam = a; rpos = p; r_valid = true; }
                else if (a <= lam) { /* certified: no update */ }
                else { search_inconclusive = true; break; }
            }
            if (search_inconclusive) {
                res.inconclusive_at = j;
                stopped_inconclusive = true;
                break;
            }

            pivot_decision dec = pivot_decision::inconclusive;
            bool pivot_1x1 = false;
            bool swap_j_r = false;
            Index rrow = Index(-1);

            const R absakk = abs(P[q + q * m]);
            if (!r_valid) {
                dec = pivot_decision::acceptable;
                pivot_1x1 = true;
            } else {
                rrow = VCP_SLDL_GROW(rpos);
                const R alam = alpha * lam;
                if (absakk >= alam) {                        // rule 1
                    dec = pivot_decision::acceptable;
                    pivot_1x1 = true;
                } else if (absakk < alam) {
                    // Rules 2-4 all need sigma, i.e. the fully updated column
                    // r.  A left-looking panel holds it only when r is a panel
                    // COLUMN.  When r lies in the IMMEDIATELY FOLLOWING
                    // supernode the panel is extended by a boundary split
                    // (D-4) and the test is retried with r inside -- this is
                    // the same mechanism the crossing 2x2 needs, and it is
                    // what makes narrow (width-1) supernodes workable at all.
                    // Anything further away cannot be rescued: stop honestly
                    // (D-1), never weaken the test, never fall back silently.
                    if (rrow >= pb) {
                        want_split = true;
                        split_partner = rrow;
                        break;
                    }
                    const std::size_t lr = static_cast<std::size_t>(rrow - pa);
                    R sigma = R(0);
                    bool sigma_inconclusive = false;
                    for (std::size_t p = q; p < m; ++p) {
                        if (p == lr) continue;
                        const Index i = VCP_SLDL_GROW(p);
                        const R a = (i > rrow) ? abs(P[p + lr * m])
                                               : abs(P[lr + static_cast<std::size_t>(i - pa) * m]);
                        if (a > sigma) { sigma = a; }
                        else if (a <= sigma) { /* certified: no update */ }
                        else { sigma_inconclusive = true; break; }
                    }
                    if (!sigma_inconclusive) {
                        const R lhs2 = absakk * sigma;
                        const R rhs2 = alam * lam;
                        if (lhs2 >= rhs2) {                  // rule 2
                            dec = pivot_decision::acceptable;
                            pivot_1x1 = true;
                        } else if (lhs2 < rhs2) {
                            const R absarr = abs(P[lr + lr * m]);
                            const R asig = alpha * sigma;
                            if (absarr >= asig) {            // rule 3
                                dec = pivot_decision::acceptable;
                                pivot_1x1 = true;
                                swap_j_r = true;
                            } else if (absarr < asig) {      // rule 4
                                dec = pivot_decision::acceptable;
                                pivot_1x1 = false;
                            }
                        }
                    }
                }
            }

            if (dec == pivot_decision::inconclusive) {
                res.inconclusive_at = j;
                stopped_inconclusive = true;
                break;
            }

            // ---- symmetric exchange of two PANEL COLUMNS p < q (global
            // labels), design v1 SS3.3.  The panel is swapped physically; the
            // rows of finished supernodes are relabelled through rowloc (the
            // baseline label-map method), so no committed data moves and no
            // structure grows.
            if (pivot_1x1) {
                if (swap_j_r && rrow != j) {
                    sparse_ldl_detail::sparse_ldl_panel_exchange_(
                        j, rrow, pa, m, P, res.perm,
                        cur_of_asm, asm_of_cur, rowloc, Rows);
                }
                const T d = P[q + q * m];
                const R absd = abs(d);
                if (absd > opt.zero_pivot_tol) {
                    res.D_diag[static_cast<std::size_t>(j)] = d;
                    if (absd > growth_max_d) growth_max_d = absd;
                    P[q + q * m] = T(1);
                    for (std::size_t p = q + 1u; p < m; ++p) P[p + q * m] = P[p + q * m] / d;
                    for (std::size_t cc = q + 1u; cc < w; ++cc) {
                        const T wq = P[cc + q * m] * d;
                        for (std::size_t p = cc; p < m; ++p) {
                            P[p + cc * m] -= P[p + q * m] * wq;
                        }
                    }
                    ++res.n_pivots_1x1;
                } else if (absd <= opt.zero_pivot_tol) {
                    P[q + q * m] = T(1);
                    for (std::size_t p = q + 1u; p < m; ++p) P[p + q * m] = T(0);
                    if (!any_zero_pivot) { res.first_zero_pivot = j; any_zero_pivot = true; }
                    ++res.n_zero_skips;
                } else {
                    res.inconclusive_at = j;
                    stopped_inconclusive = true;
                    break;
                }
                ++q;
                continue;
            }

            // ---- rule 4: 2x2 pivot at (j, j+1) after exchanging (j+1 <-> r).
            if (q + 1u >= w) {
                // Unreachable by construction: j is the last panel column, so
                // any r > j satisfies r >= pb and the boundary split above has
                // already been taken.  Kept as an honest stop rather than an
                // assumption (P3).
                want_split = true;
                split_partner = rrow;
                break;
            }

            if (rrow != j + Index(1)) {
                sparse_ldl_detail::sparse_ldl_panel_exchange_(
                    j + Index(1), rrow, pa, m, P, res.perm,
                    cur_of_asm, asm_of_cur, rowloc, Rows);
            }
            {
                const T d1 = P[q + q * m];
                const T e  = P[q + 1u + q * m];
                const T d2 = P[q + 1u + (q + 1u) * m];
                const T det = d1 * d2 - e * e;
                const R absdet = abs(det);
                if (absdet > opt.zero_pivot_tol) {
                    res.D_diag[static_cast<std::size_t>(j)] = d1;
                    res.D_diag[static_cast<std::size_t>(j) + 1u] = d2;
                    res.D_sub[static_cast<std::size_t>(j)] = e;
                    res.D_block2[static_cast<std::size_t>(j)] = char(1);
                    { const R a1 = abs(d1), a2 = abs(d2), ae = abs(e);
                      if (a1 > growth_max_d) growth_max_d = a1;
                      if (a2 > growth_max_d) growth_max_d = a2;
                      if (ae > growth_max_d) growth_max_d = ae; }
                    for (std::size_t p = q + 2u; p < m; ++p) {
                        const T a1 = P[p + q * m];
                        const T a2 = P[p + (q + 1u) * m];
                        P[p + q * m]        = (d2 * a1 - e * a2) / det;
                        P[p + (q + 1u) * m] = (d1 * a2 - e * a1) / det;
                    }
                    P[q + q * m] = T(1);
                    P[q + 1u + q * m] = T(0);          // belongs to D, not to L
                    P[q + 1u + (q + 1u) * m] = T(1);
                    for (std::size_t cc = q + 2u; cc < w; ++cc) {
                        const T l1c = P[cc + q * m];
                        const T l2c = P[cc + (q + 1u) * m];
                        const T w1 = l1c * d1 + l2c * e;
                        const T w2 = l1c * e + l2c * d2;
                        for (std::size_t p = cc; p < m; ++p) {
                            P[p + cc * m] -= P[p + q * m] * w1 + P[p + (q + 1u) * m] * w2;
                        }
                    }
                    ++res.n_pivots_2x2;
                } else if (absdet <= opt.zero_pivot_tol) {
                    res.D_block2[static_cast<std::size_t>(j)] = char(1);
                    P[q + q * m] = T(1);
                    P[q + 1u + (q + 1u) * m] = T(1);
                    for (std::size_t p = q + 1u; p < m; ++p) P[p + q * m] = T(0);
                    for (std::size_t p = q + 2u; p < m; ++p) P[p + (q + 1u) * m] = T(0);
                    if (!any_zero_pivot) { res.first_zero_pivot = j; any_zero_pivot = true; }
                    res.n_zero_skips += Index(2);
                } else {
                    res.inconclusive_at = j;
                    stopped_inconclusive = true;
                    break;
                }
            }
            q += 2u;
        }

        // ---- boundary split (D-4): the BK test at column j needs a column of
        // the immediately following supernode (as the 2x2 partner, or simply
        // to evaluate sigma).  The supernode is cut at j and j is merged into
        // the successor's panel, where the test is retried with the candidate
        // inside.  For the finished part this is METADATA ONLY: with width
        // w-1 and row list {j} u R_s the leading dimension (w-1)+(h+1) = w+h
        // is unchanged and global row j already sits at panel position w-1.
        if (want_split) {
            const Index j = pa + static_cast<Index>(q);
            const std::size_t sn2 = ss + 1u;
            bool can_split = (sn2 < nsu) && (sn_c0[sn2] == pb) && (sn_c1[sn2] > pb);
            // the BK candidate must become a column of the merged panel
            if (can_split) {
                can_split = (split_partner >= pb) && (split_partner < sn_c1[sn2]);
            }
            // rows of this panel beyond the merged panel's columns must already
            // be rows of the successor, otherwise the merged panel would need a
            // structure the symbolic phase never allocated
            if (can_split) {
                const Index nb = sn_c1[sn2];
                const std::vector<Index>& Rn = Rows[sn2];
                std::size_t t2 = 0;
                for (std::size_t t = 0; t < h && can_split; ++t) {
                    if (Rs[t] < nb) continue;
                    while (t2 < Rn.size() && Rn[t2] < Rs[t]) ++t2;
                    if (t2 >= Rn.size() || Rn[t2] != Rs[t]) can_split = false;
                }
            }
            if (!can_split) {
                res.out_of_panel_at = j;
                stopped_out_of_panel = true;
            } else {
                carry_rows.clear();
                carry_vals.clear();
                carry_rows.push_back(j);
                carry_vals.push_back(P[q + q * m]);
                for (std::size_t t = 0; t < h; ++t) {
                    carry_rows.push_back(Rs[t]);
                    carry_vals.push_back(P[w + t + q * m]);
                }
                carried = true;

                std::vector<Index> merged;
                merged.reserve(h + Rows[sn2].size());
                const Index nb = sn_c1[sn2];
                for (std::size_t t = 0; t < h; ++t) if (Rs[t] >= nb) merged.push_back(Rs[t]);
                merged.insert(merged.end(), Rows[sn2].begin(), Rows[sn2].end());
                std::sort(merged.begin(), merged.end());
                merged.erase(std::unique(merged.begin(), merged.end()), merged.end());

                std::vector<Index> newrows;
                newrows.reserve(h + 1u);
                newrows.push_back(j);
                newrows.insert(newrows.end(), Rs.begin(), Rs.end());

                sn_c1[ss] = j;                    // finished part shrinks
                Rows[ss].swap(newrows);
                sn_c0[sn2] = j;                   // successor grows
                Rows[sn2].swap(merged);
                col_sn[static_cast<std::size_t>(j)] = static_cast<Index>(sn2);
                ++res.n_boundary_splits;
                split_here = true;
            }
        }

        if (stopped_inconclusive || stopped_out_of_panel) break;

        // ---- register this supernode as a future descendant.  After a split
        // its first row is the carried column j, whose update was already
        // applied inside this panel, so the cursor starts past it.
        {
            const std::vector<Index>& Rf = Rows[ss];
            const Index start = split_here ? Index(1) : Index(0);
            if (sn_c1[ss] > sn_c0[ss] && static_cast<std::size_t>(start) < Rf.size()) {
                const Index next_sn = col_sn[static_cast<std::size_t>(Rf[static_cast<std::size_t>(start)])];
                rptr[ss] = start;
                link[ss] = head[static_cast<std::size_t>(next_sn)];
                head[static_cast<std::size_t>(next_sn)] = static_cast<Index>(s);
            } else {
                rptr[ss] = static_cast<Index>(Rf.size());
            }
            // row -> (supernode, position) links for later relabelling
            if (!static_mode && sn_c1[ss] > sn_c0[ss]) {
                for (std::size_t t = 0; t < Rf.size(); ++t) {
                    rowloc[static_cast<std::size_t>(Rf[t])].push_back(
                        std::make_pair(static_cast<Index>(s), static_cast<Index>(t)));
                }
            }
        }

        // ---- release the relative index map (touched entries only)
        for (std::size_t r = 0; r < w; ++r) relind[static_cast<std::size_t>(pa) + r] = Index(-1);
        for (std::size_t t = 0; t < h; ++t) relind[static_cast<std::size_t>(Rs[t])] = Index(-1);
#undef VCP_SLDL_GROW
    }

    if (stopped_out_of_panel) {
        // D-1: no partial result is returned, only the integer diagnostic.
        res.L_col_ptr.clear(); res.L_row_ind.clear(); res.L_val.clear();
        res.D_diag.clear(); res.D_sub.clear(); res.D_block2.clear();
        res.perm.clear();
        res.status = sparse_ldl_status::pivot_out_of_panel;
        return;
    }
    if (stopped_inconclusive) {
        res.status = sparse_ldl_status::inconclusive_pivot_test;
        return;
    }

    // ---- L extraction (CSC, explicit unit diagonal, certified zeros dropped).
    // Row labels of finished supernodes may have been permuted by exchanges,
    // so each column is sorted by row (integer keys, P6-safe) exactly as the
    // baseline kernel does.
    res.L_col_ptr.assign(un + 1u, Index(0));
    std::vector<std::pair<Index, std::size_t> > order;
    for (Index s = Index(0); s < nsup; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        const Index pa = sn_c0[ss];
        const Index pb = sn_c1[ss];
        if (pb <= pa) continue;
        const std::size_t w = static_cast<std::size_t>(pb - pa);
        const std::vector<Index>& Rs = Rows[ss];
        const std::size_t h = Rs.size();
        const std::size_t m = w + h;
        const T* const P = &Lblk[ss][0];
        for (std::size_t q = 0; q < w; ++q) {
            const Index j = pa + static_cast<Index>(q);
            res.L_row_ind.push_back(j);
            res.L_val.push_back(T(1));
            order.clear();
            for (std::size_t p = q + 1u; p < m; ++p) {
                order.push_back(std::make_pair(
                    (p < w) ? (pa + static_cast<Index>(p)) : Rs[p - w], p));
            }
            std::sort(order.begin(), order.end());
            for (std::size_t t = 0; t < order.size(); ++t) {
                const T& x = P[order[t].second + q * m];
                if (x == T(0)) { /* certified zero: not stored */ }
                else {
                    res.L_row_ind.push_back(order[t].first);
                    res.L_val.push_back(x);
                }
            }
            res.L_col_ptr[static_cast<std::size_t>(j) + 1u] =
                static_cast<Index>(res.L_row_ind.size());
        }
    }
    res.nnz_L = static_cast<Index>(res.L_row_ind.size());

    sparse_ldl_detail::sparse_ldl_set_growth_(
        growth_max_d, growth_max_a, growth_max_a_valid, res);

    res.status = any_zero_pivot ? sparse_ldl_status::zero_pivot
                                : sparse_ldl_status::success;
}

#endif // VCP_TSPARSE_SPARSE_LDL_SUPERNODAL_IMPL_HPP
