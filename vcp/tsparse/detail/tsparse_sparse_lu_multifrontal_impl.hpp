// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-MF: Multifrontal numeric source for method=supernodal.
//
// This file owns the *production* numeric source switch for method=supernodal.
// It replaces the left-looking A_eff-origin driver
// (sparse_lu_factorize_supernodal_from_a_eff, SLU-8R.5.5) as the numeric source.
// The left-looking driver is retained as a diagnostic (opt-in via
// sparse_lu_options::supernodal_numeric_diagnostic_leftlooking).
//
// Design: sandbox/docs/design/SLU_multifrontal_design.md.
//
// Algorithm (multifrontal LU with static pivoting):
//   - Assembly tree: supernodes are processed in their natural (postorder)
//     elimination order (0..nsup-1, increasing first_col). Each supernode is a
//     "front".
//   - Frontal matrix F (dense, column-major, R x ncols):
//       rows    = supernode row_indices (diagonal block rows + off-diagonal L rows)
//       cols    = pivot columns C ++ "U-right" columns (columns to the right that
//                 carry an off-diagonal U entry in a pivot row of this front).
//     F is assembled from A_eff (sparse scatter of A_csc columns) plus the
//     blocked extend-add of children contribution blocks.
//   - Dense partial factorization (BLAS-3): block LU of the pivot block, two
//     trsm (L\A12, A21/U), one gemm (Schur complement A22 -= A21*A12). The gemm
//     is the dominant flop, which is the design-§17.2-line-894 goal.
//   - Blocked extend-add: each front emits one dense contribution block (the
//     Schur complement) routed to the ancestor front that eliminates the
//     smallest index it touches. Absorption is a block-structured add (NOT a
//     per-(k,j) scatter), which is the mechanism that removes the left-looking
//     scatter-bound bottleneck.
//   - Output: panel_values / U_segments / row_perm / Dr / Dc in the SAME
//     supernodal_lu_storage contract consumed by §18.2 storage-native solve
//     (tsparse_sparse_lu_supernodal_solve_impl.hpp). The solve is NOT modified.
//
// Static pivoting:
//   row_perm is inherited from the CSC bootstrap (baseline GP pivot order, plus
//   any MC64 static pivoting / equilibration already applied). The front pivot
//   block uses its diagonal entries; no cross-supernode row swaps are performed,
//   which keeps the prebuilt symbolic structure (row_indices / U_segments) valid.
//   A zero / near-zero diagonal pivot is a safe failure (true_numeric_source
//   stays false; CSC fallback remains in effect; valid()==false at the factor
//   level if the whole factorization cannot proceed). The post-factorization
//   sparse residual check (reused from the left-looking driver) is the final
//   acceptance gate, so any structural inconsistency degrades to CSC fallback
//   rather than producing a silent wrong answer.
//
// PROHIBITIONS (§25): dense_kernel::getrf is NOT used for pivot search.
//
// Invariants (SLU_multifrontal_design.md §4):
//   MF-1 non-destructive (baseline_gp / auto_select byte-identical; this driver
//        runs only inside method=supernodal).
//   MF-2 baseline_gp backward equivalence + original-system residual.
//   MF-3 §18.2 solve reuse (output matches the existing storage contract).
//   MF-4 safe failure preserved.
//   MF-5 T generic / point arithmetic.
//   MF-7 tblas/tlapack untouched (adapter only).
//
// This file MUST be #included from WITHIN namespace vcp, AFTER
// tsparse_sparse_lu_true_numeric_impl.hpp (it reuses a_eff_lookup and
// compute_supernodal_lu_residual_sparse).
//
// Do NOT include this file directly. Include <vcp/tsparse/tsparse_sparse_lu.hpp>.

#ifndef VCP_TSPARSE_SPARSE_LU_MULTIFRONTAL_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_MULTIFRONTAL_IMPL_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// mf_urow_index<T, Index>
//
// Transpose-of-U_segments index: for each global (permuted) row p, the list of
// off-diagonal U entries U[p, c'] stored in U_segments. Each entry records the
// absolute U_segments value index (scatter target) and the global column c'.
//
// These are precisely the "U-to-the-right" values produced by the front that
// owns p as a pivot row. Built once (O(nnz of U_segments)).
// ---------------------------------------------------------------------------
template <class Index>
struct mf_urow_index {
    std::vector<std::size_t> ptr;      // size n + 1 (CSR over rows)
    std::vector<std::size_t> abs_idx;  // U_segments value index
    std::vector<Index>       col;      // global new_col c'
};

template <class T, class Index>
void mf_build_urow_index(
    const supernodal_lu_storage<T, Index>& storage,
    Index                                  n,
    mf_urow_index<Index>&                  out)
{
    const std::size_t sn   = (n > Index(0)) ? static_cast<std::size_t>(n) : 0u;
    const std::size_t nsup = storage.supernodes.size();

    out.ptr.assign(sn + 1u, 0u);
    out.abs_idx.clear();
    out.col.clear();

    // Pass 1: count entries per row.
    for (std::size_t j = 0u; j < nsup; ++j) {
        const supernode_desc<Index>& desc = storage.supernodes[j];
        if (j >= storage.U_segments.seg_ptr.size()) continue;
        const Index seg_start = storage.U_segments.seg_ptr[j];
        if (desc.u_seg_col_ptr.size() <
                static_cast<std::size_t>(desc.num_cols + 1)) continue;

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const std::size_t sc = static_cast<std::size_t>(c);
            const Index rel_s = desc.u_seg_col_ptr[sc];
            const Index rel_e = desc.u_seg_col_ptr[sc + 1u];
            for (Index ki = rel_s; ki < rel_e; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(seg_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                if (u_row < Index(0) ||
                    static_cast<std::size_t>(u_row) >= sn) continue;
                ++out.ptr[static_cast<std::size_t>(u_row) + 1u];
            }
        }
    }

    for (std::size_t i = 0u; i < sn; ++i) out.ptr[i + 1u] += out.ptr[i];

    const std::size_t total = out.ptr[sn];
    out.abs_idx.resize(total);
    out.col.resize(total);
    std::vector<std::size_t> cursor(out.ptr.begin(), out.ptr.end() - 1);

    // Pass 2: fill.
    for (std::size_t j = 0u; j < nsup; ++j) {
        const supernode_desc<Index>& desc = storage.supernodes[j];
        const Index col_begin = desc.first_col;
        if (j >= storage.U_segments.seg_ptr.size()) continue;
        const Index seg_start = storage.U_segments.seg_ptr[j];
        if (desc.u_seg_col_ptr.size() <
                static_cast<std::size_t>(desc.num_cols + 1)) continue;

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            const std::size_t sc = static_cast<std::size_t>(c);
            const Index rel_s = desc.u_seg_col_ptr[sc];
            const Index rel_e = desc.u_seg_col_ptr[sc + 1u];
            for (Index ki = rel_s; ki < rel_e; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(seg_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                if (u_row < Index(0) ||
                    static_cast<std::size_t>(u_row) >= sn) continue;
                const std::size_t slot = cursor[static_cast<std::size_t>(u_row)]++;
                out.abs_idx[slot] = abs_idx;
                out.col[slot]     = new_col;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// mf_contribution_block<T, Index>
//
// A dense Schur-complement contribution block emitted by a factored front.
// rows/cols are global (permuted) indices (sorted ascending). vals is
// column-major (leading dimension == rows.size()).
// ---------------------------------------------------------------------------
template <class T, class Index>
struct mf_contribution_block {
    std::vector<Index> rows;
    std::vector<Index> cols;
    std::vector<T>     vals;
};


// ---------------------------------------------------------------------------
// mf_route_target
//
// Returns the supernode that eliminates the smallest index touched by a
// contribution block (rows ++ cols, both sorted ascending). That supernode is
// where the block must be extend-added next. Returns -1 if empty / out of range.
// ---------------------------------------------------------------------------
template <class Index>
Index mf_route_target(
    const std::vector<Index>&  rows,
    const std::vector<Index>&  cols,
    const std::vector<Index>&  col_to_supernode)
{
    Index mn = Index(-1);
    if (!rows.empty()) mn = rows.front();
    if (!cols.empty() && (mn < Index(0) || cols.front() < mn)) mn = cols.front();
    if (mn < Index(0) ||
        static_cast<std::size_t>(mn) >= col_to_supernode.size())
        return Index(-1);
    return col_to_supernode[static_cast<std::size_t>(mn)];
}

// ---------------------------------------------------------------------------
// mf_struct_block<Index>
//
// SLU-MF2: a SYMBOLIC contribution block (structure only, no values). rows and
// cols are global (permuted) indices, sorted ascending. This mirrors the numeric
// mf_contribution_block but carries the dense Schur-complement *pattern* O_s x Ur_s.
// ---------------------------------------------------------------------------
template <class Index>
struct mf_struct_block {
    std::vector<Index> rows;
    std::vector<Index> cols;
};

// ---------------------------------------------------------------------------
// build_supernodal_self_symbolic_storage  (SLU-MF2 core)
//
// Builds a supernodal_lu_storage whose row_indices / U_segments / panel
// allocation contain the *dense multifrontal (relaxed-front) fill*, instead of
// the GP-exact sparse fill the CSC bootstrap produces. The structure is computed
// by propagating dense contribution-block patterns up the assembly tree,
// replicating the EXACT routing used by the numeric multifrontal driver
// (sparse_lu_factorize_supernodal_multifrontal). Consequence: whenever the
// numeric driver would succeed, this structure already contains every entry the
// driver touches, so extend-add never forwards and never falls back -- wide
// (relaxed-amalgamation / AMD) fronts run native instead of degrading to CSC.
//
// The supernode PARTITION, permutations (row/col_perm), and equilibration
// (Dr/Dc) are taken verbatim from `templ` (the CSC bootstrap storage), so the
// ordering / equilibration / pivot framework is unchanged. Only the symbolic
// fill structure (row_indices / U_segments) and the panel allocation differ.
//
// The dense-front model produces a structural SUPERSET of GP-exact fill (extra
// positions are explicit zeros within fronts), so L*U still reconstructs A_eff
// and the numeric residual check is unaffected.
//
// On any structural inconsistency the result is returned with valid == false,
// which makes the caller keep the CSC bootstrap structure (safe, non-breaking).
// ---------------------------------------------------------------------------
template <class T, class Index>
supernodal_lu_storage<T, Index>
build_supernodal_self_symbolic_storage(
    const csc_storage<T, Index>&            A_csc,
    const supernodal_lu_storage<T, Index>&  templ,
    Index                                   n)
{
    static_assert(std::is_signed<Index>::value,
                  "build_supernodal_self_symbolic_storage: Index must be signed");

    supernodal_lu_storage<T, Index> R;
    R.valid                   = false;
    R.bootstrapped_from_csc   = false;
    R.source_of_truth_storage = false;
    R.true_numeric_source     = false;

    if (!templ.valid || templ.supernodes.empty() || n <= Index(0)) return R;
    const std::size_t sn   = static_cast<std::size_t>(n);
    const std::size_t nsup = templ.supernodes.size();
    if (A_csc.col_ptr.size() < sn + 1u) return R;

    // Carry over partition-independent data verbatim.
    R.row_perm     = templ.row_perm;
    R.inv_row_perm = templ.inv_row_perm;
    R.col_perm     = templ.col_perm;
    R.inv_col_perm = templ.inv_col_perm;
    R.Dr           = templ.Dr;
    R.Dc           = templ.Dc;
    R.supernodes.resize(nsup);
    // SLU-MF6: export per-front U-right column sets so the in-place driver does
    // not rebuild them from the U_segments transpose (the driver's single most
    // expensive symbolic segment). Filled per front below from ucolset (Ur_s).
    R.mf_uright_of.resize(nsup);

    // col -> supernode, and inv_row_perm[old]=new (rebuild locally; do not trust
    // a possibly-empty templ.inv_row_perm).
    std::vector<Index> col_to_supernode(sn, Index(-1));
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& d = templ.supernodes[s];
        const Index cb = d.first_col;
        const Index ce = cb + d.num_cols;
        if (cb < Index(0) || d.num_cols <= Index(0) ||
            static_cast<std::size_t>(ce) > sn) return R;
        for (Index c = cb; c < ce; ++c)
            col_to_supernode[static_cast<std::size_t>(c)] = static_cast<Index>(s);
    }
    std::vector<Index> inv_row_perm(sn, Index(-1));
    for (std::size_t nr = 0u; nr < templ.row_perm.size() && nr < sn; ++nr) {
        const Index old_row = templ.row_perm[nr];
        if (old_row >= Index(0) && static_cast<std::size_t>(old_row) < sn)
            inv_row_perm[static_cast<std::size_t>(old_row)] = static_cast<Index>(nr);
    }

    // A_eff row access (CSR over new rows): rowcols[p] = sorted new cols c with
    // A_eff[p,c] structurally present. Built once via counting sort, O(nnz).
    std::vector<std::size_t> rc_ptr(sn + 1u, 0u);
    for (std::size_t c = 0u; c < sn; ++c) {
        const Index cs = A_csc.col_ptr[c];
        const Index ce = A_csc.col_ptr[c + 1u];
        for (Index idx = cs; idx < ce; ++idx) {
            const std::size_t sidx = static_cast<std::size_t>(idx);
            if (sidx >= A_csc.row_ind.size()) continue;
            const Index old_row = A_csc.row_ind[sidx];
            if (old_row < Index(0) || static_cast<std::size_t>(old_row) >= sn) continue;
            const Index p = inv_row_perm[static_cast<std::size_t>(old_row)];
            if (p < Index(0)) continue;
            ++rc_ptr[static_cast<std::size_t>(p) + 1u];
        }
    }
    for (std::size_t i = 0u; i < sn; ++i) rc_ptr[i + 1u] += rc_ptr[i];
    std::vector<Index> rc_col(rc_ptr[sn], Index(0));
    {
        std::vector<std::size_t> cur(rc_ptr.begin(), rc_ptr.end() - 1);
        for (std::size_t c = 0u; c < sn; ++c) {
            const Index cs = A_csc.col_ptr[c];
            const Index ce = A_csc.col_ptr[c + 1u];
            for (Index idx = cs; idx < ce; ++idx) {
                const std::size_t sidx = static_cast<std::size_t>(idx);
                if (sidx >= A_csc.row_ind.size()) continue;
                const Index old_row = A_csc.row_ind[sidx];
                if (old_row < Index(0) || static_cast<std::size_t>(old_row) >= sn) continue;
                const Index p = inv_row_perm[static_cast<std::size_t>(old_row)];
                if (p < Index(0)) continue;
                rc_col[cur[static_cast<std::size_t>(p)]++] = static_cast<Index>(c);
            }
        }
    }

    // Pending symbolic contribution blocks bucketed by absorbing supernode.
    std::vector<std::vector<mf_struct_block<Index> > > pending(nsup);

    // Per-column U-above contributors: ucol_rows[c] collects the pivot rows
    // r < first_col(supernode(c)) that produce a (dense-model) U[r,c] entry.
    std::vector<std::vector<Index> > ucol_rows(sn);

    // Reusable dedup markers (touched-list cleared).
    std::vector<char> rmark(sn, 0), cmark(sn, 0);
    std::vector<Index> rowset, ucolset, touched_r, touched_c;

    // ---- Pass 1: per-front row_indices, Ur, panel allocation, propagation ----
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& td = templ.supernodes[s];
        const Index fc = td.first_col;
        const Index w  = td.num_cols;
        const Index ce = fc + w;
        const std::size_t sw = static_cast<std::size_t>(w);

        rowset.clear(); touched_r.clear();
        ucolset.clear(); touched_c.clear();
        const auto add_row = [&](Index r) {
            if (r < fc || static_cast<std::size_t>(r) >= sn) return;
            if (!rmark[static_cast<std::size_t>(r)]) {
                rmark[static_cast<std::size_t>(r)] = 1;
                rowset.push_back(r); touched_r.push_back(r);
            }
        };
        const auto add_ucol = [&](Index c) {
            if (c < ce || static_cast<std::size_t>(c) >= sn) return; // U-right only
            if (!cmark[static_cast<std::size_t>(c)]) {
                cmark[static_cast<std::size_t>(c)] = 1;
                ucolset.push_back(c); touched_c.push_back(c);
            }
        };

        // pivot rows.
        for (Index r = fc; r < ce; ++r) add_row(r);
        // A_eff L-rows (pivot columns) and U-right (pivot rows).
        for (Index j = fc; j < ce; ++j) {
            const std::size_t scol = static_cast<std::size_t>(j);
            const Index cs = A_csc.col_ptr[scol];
            const Index cee = A_csc.col_ptr[scol + 1u];
            for (Index idx = cs; idx < cee; ++idx) {
                const std::size_t sidx = static_cast<std::size_t>(idx);
                if (sidx >= A_csc.row_ind.size()) continue;
                const Index old_row = A_csc.row_ind[sidx];
                if (old_row < Index(0) || static_cast<std::size_t>(old_row) >= sn) continue;
                const Index p = inv_row_perm[static_cast<std::size_t>(old_row)];
                add_row(p); // r >= fc accepted (off-diag L + pivot block)
            }
        }
        for (Index p = fc; p < ce; ++p) {
            const std::size_t sp = static_cast<std::size_t>(p);
            for (std::size_t e = rc_ptr[sp]; e < rc_ptr[sp + 1u]; ++e)
                add_ucol(rc_col[e]); // c >= ce accepted (U-right)
        }
        // children dense contribution blocks (extend-add of structure).
        for (std::size_t ci = 0u; ci < pending[s].size(); ++ci) {
            const mf_struct_block<Index>& cb = pending[s][ci];
            for (std::size_t i = 0u; i < cb.rows.size(); ++i) add_row(cb.rows[i]);
            for (std::size_t j = 0u; j < cb.cols.size(); ++j) add_ucol(cb.cols[j]);
        }
        std::vector<mf_struct_block<Index> >().swap(pending[s]);

        std::sort(rowset.begin(),  rowset.end());
        std::sort(ucolset.begin(), ucolset.end());

        // SLU-MF6: export the sorted-unique U-right column set for this front so
        // the in-place driver reuses it instead of rebuilding from U_segments.
        R.mf_uright_of[s] = ucolset;

        // record row_indices + panel allocation.
        supernode_desc<Index>& desc = R.supernodes[s];
        desc.first_col = fc;
        desc.num_cols  = w;
        desc.row_indices = rowset; // sorted unique, all >= fc, contains [fc,ce)
        const Index row_count = static_cast<Index>(desc.row_indices.size());
        desc.values_offset     = static_cast<Index>(R.panel_values.size());
        desc.leading_dimension = row_count;
        R.panel_values.resize(R.panel_values.size() +
            static_cast<std::size_t>(row_count) * sw, T(0));

        // record U-above contributors for each U-right col (dense model: every
        // pivot row p in [fc,ce) pairs with every U-right col c).
        for (std::size_t e = 0u; e < ucolset.size(); ++e) {
            std::vector<Index>& dst = ucol_rows[static_cast<std::size_t>(ucolset[e])];
            for (Index p = fc; p < ce; ++p) dst.push_back(p);
        }

        // emit dense contribution block O_s x Ur_s; route exactly like numeric.
        std::vector<Index> O_s, Ur_s;
        O_s.reserve(rowset.size());
        for (std::size_t e = 0u; e < rowset.size(); ++e)
            if (rowset[e] >= ce) O_s.push_back(rowset[e]);
        Ur_s = ucolset; // already only >= ce, sorted
        if (!O_s.empty() && !Ur_s.empty()) {
            const Index r0 = O_s.front();
            const Index c0 = Ur_s.front();
            const Index key = (r0 < c0) ? r0 : c0;
            const Index tgt = (static_cast<std::size_t>(key) < sn)
                ? col_to_supernode[static_cast<std::size_t>(key)] : Index(-1);
            if (tgt > static_cast<Index>(s) && static_cast<std::size_t>(tgt) < nsup) {
                mf_struct_block<Index> blk;
                blk.rows.swap(O_s);
                blk.cols = Ur_s;
                pending[static_cast<std::size_t>(tgt)].push_back(std::move(blk));
            }
            // tgt <= s / out of range: do NOT propagate. The numeric driver hits
            // the same routing and safely falls back for that case; our structure
            // for front s itself remains self-consistent.
        }

        // clear markers.
        for (std::size_t e = 0u; e < touched_r.size(); ++e)
            rmark[static_cast<std::size_t>(touched_r[e])] = 0;
        for (std::size_t e = 0u; e < touched_c.size(); ++e)
            cmark[static_cast<std::size_t>(touched_c[e])] = 0;
    }

    // ---- Pass 2: assemble U_segments column-by-column per supernode ----
    R.U_segments.seg_ptr.assign(nsup + 1u, Index(0));
    for (std::size_t s = 0u; s < nsup; ++s) {
        supernode_desc<Index>& desc = R.supernodes[s];
        const Index fc = desc.first_col;
        const Index w  = desc.num_cols;
        const std::size_t sw = static_cast<std::size_t>(w);

        R.U_segments.seg_ptr[s] =
            static_cast<Index>(R.U_segments.row_ind.size());
        desc.u_segment_start = R.U_segments.seg_ptr[s];
        desc.u_seg_col_ptr.assign(sw + 1u, Index(0));

        Index col_off = Index(0);
        for (Index c = fc; c < fc + w; ++c) {
            desc.u_seg_col_ptr[static_cast<std::size_t>(c - fc)] = col_off;
            std::vector<Index>& rows = ucol_rows[static_cast<std::size_t>(c)];
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
            for (std::size_t e = 0u; e < rows.size(); ++e) {
                const Index r = rows[e];
                if (r >= fc) continue; // U-above only (defensive; should hold)
                R.U_segments.row_ind.push_back(r);
                R.U_segments.values.push_back(T(0));
                ++col_off;
            }
        }
        desc.u_seg_col_ptr[sw] = col_off;
        desc.u_segment_count = col_off;
    }
    R.U_segments.seg_ptr[nsup] =
        static_cast<Index>(R.U_segments.row_ind.size());

    // ---- structural invariant self-check (SLU-MF2 §12.5) ----
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& d = R.supernodes[s];
        const std::size_t rcnt = d.row_indices.size();
        if (static_cast<std::size_t>(d.leading_dimension) < rcnt) return R;
        if (static_cast<std::size_t>(d.values_offset) +
            static_cast<std::size_t>(d.leading_dimension) *
            static_cast<std::size_t>(d.num_cols) > R.panel_values.size()) return R;
        // row_indices sorted strictly increasing and must contain the pivot block.
        for (std::size_t i = 1u; i < rcnt; ++i)
            if (!(d.row_indices[i - 1u] < d.row_indices[i])) return R;
        if (rcnt < static_cast<std::size_t>(d.num_cols)) return R;
        for (Index k = Index(0); k < d.num_cols; ++k)
            if (d.row_indices[static_cast<std::size_t>(k)] != d.first_col + k) return R;
        // u_seg_col_ptr monotone, terminal == u_segment_count.
        if (d.u_seg_col_ptr.size() !=
            static_cast<std::size_t>(d.num_cols) + 1u) return R;
        for (std::size_t k = 1u; k < d.u_seg_col_ptr.size(); ++k)
            if (d.u_seg_col_ptr[k] < d.u_seg_col_ptr[k - 1u]) return R;
        if (d.u_seg_col_ptr[static_cast<std::size_t>(d.num_cols)] != d.u_segment_count)
            return R;
    }

    R.source_of_truth_storage = true;
    R.valid                   = true;
    return R;
}

// ---------------------------------------------------------------------------
// mf_relative_map<Index>  (SLU-MF5 core)
//
// Precomputed (symbolic) child->parent relative index map for the in-place
// multifrontal extend-add (push). For a source front s that emits a Schur
// contribution block O_s x Ur_s routed to parent front t, this stores -- one
// time, in the symbolic phase -- exactly the rmap/cmap that the MF3 numeric
// push used to rebuild with std::lower_bound on every front:
//
//   crel[j]   = parent-local column index of the j-th child U-right column
//               (== gc - parent.first_col for a parent pivot column, or
//                tsw + (index in parent U-right) for a parent U-right column).
//   rrel[i]   = parent-local row position (into parent.row_indices) of the
//               i-th child off-diagonal row.
//
// rrel is stored as CONTIGUOUS RUNS instead of a per-row array (§3.3): the
// child off-diagonal rows are a sorted subset of the parent rows, so rrel is
// monotone increasing and typically forms a few contiguous segments. Each run
// (run_src[k], run_dst[k], run_len[k]) means
//   parent positions [run_dst[k], run_dst[k]+run_len[k]) receive
//   child offsets   [run_src[k], run_src[k]+run_len[k]).
// The numeric push then adds each run as a CONTIGUOUS (vectorizable) block:
//   Ftp[run_dst+pc + l] += Fp[sc + run_src + l],  l in [0,run_len).
// No std::lower_bound runs in the numeric loop.
//
// coverage_ok == false means a child row/col was NOT found in the parent
// structure (a self-symbolic COVERAGE miss). The numeric driver treats this
// exactly like the old lower_bound miss: safe failure -> CSC fallback.
// ---------------------------------------------------------------------------
template <class Index>
struct mf_relative_map {
    Index              parent;       // routing target front t (-1 if none)
    bool               coverage_ok;  // false on any child row/col not in parent
    std::vector<Index> crel;         // size nu: parent-local column index
    std::vector<Index> run_src;      // child off-diag offset where a run starts
    std::vector<Index> run_dst;      // parent-local row position of run start
    std::vector<Index> run_len;      // contiguous run length
    mf_relative_map() : parent(Index(-1)), coverage_ok(true) {}
};

// ---------------------------------------------------------------------------
// mf_build_relative_maps  (SLU-MF5 symbolic precompute)
//
// Builds relmap[s] for every front s of a self-symbolic storage, mirroring the
// EXACT routing the MF3/MF5 numeric push performs (tgt = supernode(min(first
// off-diagonal row, first U-right column))). Replaces the per-front numeric
// std::lower_bound with sorted-merge passes done once, here, at symbolic time.
//
// uright_of[s] / ld_of[s] are the same precomputed per-front U-right column set
// and leading dimension the numeric driver already builds. col_to_supernode is
// the global column->supernode map. Fronts with no contribution block
// (m_off == 0 or nu == 0) get parent == -1 (no push).
//
// all_coverage_ok is set false if any front reports a coverage miss; the caller
// keeps INVARIANT/COVERAGE green only when it stays true for self-symbolic input.
// ---------------------------------------------------------------------------
template <class T, class Index>
void mf_build_relative_maps(
    const supernodal_lu_storage<T, Index>&  storage,
    const std::vector<std::vector<Index> >& uright_of,
    const std::vector<std::size_t>&         ld_of,
    const std::vector<Index>&               col_to_supernode,
    Index                                   n,
    std::vector<mf_relative_map<Index> >&   relmap,
    bool&                                   all_coverage_ok)
{
    const std::size_t nsup = storage.supernodes.size();
    const std::size_t sn   = (n > Index(0)) ? static_cast<std::size_t>(n) : 0u;
    relmap.assign(nsup, mf_relative_map<Index>());
    all_coverage_ok = true;
    (void)ld_of;

    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& desc = storage.supernodes[s];
        const Index w  = desc.num_cols;
        const std::size_t R  = desc.row_indices.size();
        const std::size_t sw = (w > Index(0)) ? static_cast<std::size_t>(w) : 0u;
        if (w <= Index(0) || R < sw) continue;
        const std::size_t m_off = R - sw;
        const std::vector<Index>& uright = uright_of[s];
        const std::size_t nu = uright.size();
        if (m_off == 0u || nu == 0u) continue; // no contribution block emitted

        mf_relative_map<Index>& rm = relmap[s];

        // routing target (identical to the numeric push).
        const Index r0  = desc.row_indices[sw];
        const Index c0  = uright[0];
        const Index key = (r0 < c0) ? r0 : c0;
        const Index tgt = (static_cast<std::size_t>(key) < sn)
            ? col_to_supernode[static_cast<std::size_t>(key)] : Index(-1);
        rm.parent = tgt;
        if (!(tgt > static_cast<Index>(s) &&
              static_cast<std::size_t>(tgt) < nsup))
            continue; // backward / out-of-range -> numeric path handles as fail

        const std::size_t t = static_cast<std::size_t>(tgt);
        const supernode_desc<Index>& td = storage.supernodes[t];
        const Index tfc = td.first_col;
        const Index tce = tfc + td.num_cols;
        const std::vector<Index>& trows   = td.row_indices;
        const std::vector<Index>& turight = uright_of[t];
        const std::size_t tsw = static_cast<std::size_t>(td.num_cols);

        // ---- rrel via sorted-merge (child off-diag rows subset of parent rows),
        //      emitted as contiguous runs. ----
        bool cov = true;
        std::size_t pr = 0u;
        const std::size_t tnr = trows.size();
        Index prev_dst = Index(-2);
        for (std::size_t i = 0u; i < m_off; ++i) {
            const Index gr = desc.row_indices[sw + i];
            while (pr < tnr && trows[pr] < gr) ++pr;
            if (pr >= tnr || trows[pr] != gr) { cov = false; break; }
            const Index dst = static_cast<Index>(pr);
            if (!rm.run_len.empty() && dst == prev_dst + Index(1)) {
                ++rm.run_len.back();
            } else {
                rm.run_src.push_back(static_cast<Index>(i));
                rm.run_dst.push_back(dst);
                rm.run_len.push_back(Index(1));
            }
            prev_dst = dst;
        }

        // ---- crel via sorted-merge: child U-right cols. Parent pivot columns
        //      [tfc,tce) come first (value-sorted), then parent U-right (>= tce). ----
        if (cov) {
            rm.crel.assign(nu, Index(-1));
            std::size_t ptu = 0u;
            const std::size_t tnu = turight.size();
            for (std::size_t j = 0u; j < nu; ++j) {
                const Index gc = uright[j];
                if (gc >= tfc && gc < tce) {
                    rm.crel[j] = gc - tfc; // parent pivot column
                } else {
                    while (ptu < tnu && turight[ptu] < gc) ++ptu;
                    if (ptu >= tnu || turight[ptu] != gc) { cov = false; break; }
                    rm.crel[j] = static_cast<Index>(tsw + ptu);
                }
            }
        }

        rm.coverage_ok = cov;
        if (!cov) {
            all_coverage_ok = false;
            // drop partial maps; numeric path safe-fails on coverage_ok == false.
            std::vector<Index>().swap(rm.run_src);
            std::vector<Index>().swap(rm.run_dst);
            std::vector<Index>().swap(rm.run_len);
            std::vector<Index>().swap(rm.crel);
        }
    }
}

} // namespace sparse_lu_detail

// ===========================================================================
// sparse_lu_factorize_supernodal_multifrontal  (public, in namespace vcp)
//
// SLU-MF primary driver. Same signature / return type as
// sparse_lu_factorize_supernodal_from_a_eff so it is a drop-in numeric source.
// On success: storage.true_numeric_source = true,
//             storage.numeric_source_kind = a_eff_true_numeric.
// ===========================================================================
template <class T, class Index>
sparse_lu_detail::supernodal_true_numeric_stats<
    typename vcp::tsparse_scalar::real_type<T>::type>
sparse_lu_factorize_supernodal_multifrontal(
    const csc_storage<T, Index>&           A_csc,
    const baseline_lu_storage<T, Index>&   csc_lu,
    supernodal_lu_storage<T, Index>&       storage,
    const sparse_lu_options<T>&            opt)
{
    (void)csc_lu; // row_perm / Dr / Dc already in storage (bootstrap copy)

    typedef sparse_lu_scalar_policy<T>            scalar_pol;
    typedef typename scalar_pol::real_type        real_type;
    typedef sparse_lu_detail::mf_contribution_block<T, Index> CB;

    const auto t_total_start = std::chrono::steady_clock::now();

    sparse_lu_detail::supernodal_true_numeric_stats<real_type> stats;
    stats.attempted = true;

    const auto finish_total_ticks = [&]() {
        stats.total_ticks = static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_total_start).count());
        stats.factorization_ticks =
            (stats.total_ticks > stats.residual_ticks)
                ? (stats.total_ticks - stats.residual_ticks) : 0u;
    };

    if (!storage.valid || storage.supernodes.empty()) {
        stats.status = supernodal_true_numeric_status::unsupported_structure;
        finish_total_ticks();
        return stats;
    }

    const std::size_t nsup = storage.supernodes.size();
    const Index n = static_cast<Index>(storage.row_perm.size());
    const std::size_t sn = (n > Index(0)) ? static_cast<std::size_t>(n) : 0u;

    if (A_csc.col_ptr.size() < sn + 1u) {
        stats.status = supernodal_true_numeric_status::unsupported_structure;
        finish_total_ticks();
        return stats;
    }

    // ------------------------------------------------------------------
    // Step 0: zero the output, build symbolic helpers.
    // ------------------------------------------------------------------
    std::fill(storage.panel_values.begin(),      storage.panel_values.end(),      T(0));
    std::fill(storage.U_segments.values.begin(), storage.U_segments.values.end(), T(0));
    stats.values_initialized_from_A           = true;
    stats.values_initialized_from_csc_numeric = false;

    const auto t_sym0 = std::chrono::steady_clock::now();

    // col -> supernode map.
    std::vector<Index> col_to_supernode(sn, Index(-1));
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& d = storage.supernodes[s];
        const Index cb = d.first_col;
        const Index ce = cb + d.num_cols;
        for (Index c = cb; c < ce; ++c)
            if (c >= Index(0) && static_cast<std::size_t>(c) < sn)
                col_to_supernode[static_cast<std::size_t>(c)] = static_cast<Index>(s);
    }

    // inverse row permutation: inv_row_perm[old_row] = new_row.
    std::vector<Index> inv_row_perm(sn, Index(-1));
    for (std::size_t nr = 0u; nr < storage.row_perm.size() && nr < sn; ++nr) {
        const Index old_row = storage.row_perm[nr];
        if (old_row >= Index(0) && static_cast<std::size_t>(old_row) < sn)
            inv_row_perm[static_cast<std::size_t>(old_row)] = static_cast<Index>(nr);
    }

    // transpose-of-U index (scatter targets for U-right values).
    sparse_lu_detail::mf_urow_index<Index> urow;
    sparse_lu_detail::mf_build_urow_index(storage, n, urow);

    stats.symbolic_ticks += static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t_sym0).count());

    // Pending contribution blocks bucketed by the supernode that will absorb them.
    std::vector<std::vector<CB> > pending(nsup);

    // Reusable per-front index maps (O(n)), cleared via touched lists.
    std::vector<Index> row_local(sn, Index(-1));
    std::vector<Index> col_local(sn, Index(-1));
    std::vector<Index> touched_rows;
    std::vector<Index> touched_cols;
    std::vector<Index> uright;     // U-right columns of the current front
    std::vector<T>     row_tmp;    // contiguous pivot-row scratch for ger

    // Reusable scratch for forwarded (partially-absorbed) contribution blocks.
    std::vector<Index> fwd_r, fwd_c;     // leftover entry rows / cols
    std::vector<T>     fwd_v;            // leftover entry values
    std::vector<Index> gb_rows, gb_cols; // compact forwarded block row/col sets

    bool ok = true;

    // ------------------------------------------------------------------
    // Step 1: factor each front in elimination (postorder) order.
    // ------------------------------------------------------------------
    for (std::size_t s = 0u; s < nsup && ok; ++s) {
        const supernode_desc<Index>& desc = storage.supernodes[s];
        const Index fc = desc.first_col;
        const Index w  = desc.num_cols;
        const std::size_t R  = desc.row_indices.size();
        const std::size_t sw = static_cast<std::size_t>(w);

        if (w <= Index(0) || R == 0u) {
            // Nothing to eliminate (degenerate); drop any pending blocks routed
            // here -- they cannot be absorbed, so the structure is malformed.
            if (!pending[s].empty()) {
                ok = false; break; }
            continue;
        }

        // ---- front column set: pivot columns ++ U-right columns ----
        uright.clear();
        for (std::size_t k = 0u; k < sw; ++k) {
            const Index p = fc + static_cast<Index>(k);
            if (p < Index(0) || static_cast<std::size_t>(p) >= sn) continue;
            const std::size_t a = urow.ptr[static_cast<std::size_t>(p)];
            const std::size_t b = urow.ptr[static_cast<std::size_t>(p) + 1u];
            for (std::size_t e = a; e < b; ++e) uright.push_back(urow.col[e]);
        }
        std::sort(uright.begin(), uright.end());
        uright.erase(std::unique(uright.begin(), uright.end()), uright.end());

        const std::size_t nu = uright.size();
        const std::size_t ncols = sw + nu;
        const std::size_t ldF = R;

        // ---- front local index maps ----
        touched_rows.clear();
        for (std::size_t r = 0u; r < R; ++r) {
            const Index gr = desc.row_indices[r];
            if (gr < Index(0) || static_cast<std::size_t>(gr) >= sn) {
                ok = false; break; }
            row_local[static_cast<std::size_t>(gr)] = static_cast<Index>(r);
            touched_rows.push_back(gr);
        }
        if (!ok) break;

        touched_cols.clear();
        for (std::size_t k = 0u; k < sw; ++k) {
            const Index gc = fc + static_cast<Index>(k);
            col_local[static_cast<std::size_t>(gc)] = static_cast<Index>(k);
            touched_cols.push_back(gc);
        }
        for (std::size_t e = 0u; e < nu; ++e) {
            const Index gc = uright[e];
            col_local[static_cast<std::size_t>(gc)] = static_cast<Index>(sw + e);
            touched_cols.push_back(gc);
        }

        // ---- assemble frontal matrix F (R x ncols, column-major) ----
        std::vector<T> F(R * ncols, T(0));

        // (a) scatter A_eff. Multifrontal assembly convention: only entries in a
        //     PIVOT ROW or PIVOT COLUMN of this front are loaded from the original
        //     matrix.  The contribution-block region (off-diagonal rows x U-right
        //     columns) is filled ONLY by children (extend-add) -- loading the
        //     original A there would double-count it against the owning front of
        //     the U-right column.  Pivot columns (cl < w) -> all front rows
        //     (pivot block + L). U-right columns (cl >= w) -> pivot rows only (U).
        const auto t_asm0 = std::chrono::steady_clock::now();
        for (std::size_t cl = 0u; cl < ncols; ++cl) {
            const Index gc = touched_cols[cl];
            const bool pivot_col = (cl < sw);
            const std::size_t scgc = static_cast<std::size_t>(gc);
            const Index cs = A_csc.col_ptr[scgc];
            const Index ce = A_csc.col_ptr[scgc + 1u];
            const T dc = (!storage.Dc.empty() && scgc < storage.Dc.size())
                ? storage.Dc[scgc] : T(1);
            for (Index idx = cs; idx < ce; ++idx) {
                const std::size_t sidx = static_cast<std::size_t>(idx);
                if (sidx >= A_csc.row_ind.size()) continue;
                const Index old_row = A_csc.row_ind[sidx];
                if (old_row < Index(0) ||
                    static_cast<std::size_t>(old_row) >= sn) continue;
                const Index new_row = inv_row_perm[static_cast<std::size_t>(old_row)];
                if (new_row < Index(0)) continue;
                const Index rl = row_local[static_cast<std::size_t>(new_row)];
                if (rl < Index(0)) continue;
                // U-right column: keep only pivot-row entries (U part).
                if (!pivot_col && static_cast<std::size_t>(rl) >= sw) continue;
                const T dr = (!storage.Dr.empty() &&
                              static_cast<std::size_t>(new_row) < storage.Dr.size())
                    ? storage.Dr[static_cast<std::size_t>(new_row)] : T(1);
                F[static_cast<std::size_t>(rl) + cl * ldF] +=
                    dr * A_csc.values[sidx] * dc;
            }
        }
        const auto t_ea0 = std::chrono::steady_clock::now();
        stats.mf_aeff_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                t_ea0 - t_asm0).count());

        // (b) blocked extend-add of children contribution blocks.
        // Each pending block is added into this front's dense matrix in one
        // block sweep (no per-(k,j) scatter). Entries whose (row,col) is not in
        // this front's structure are forwarded as a compact block to the
        // ancestor that first eliminates one of their indices, routed by the
        // minimum per-entry min(row,col) -- which is strictly > this front's
        // last column, so forwarding always advances (never backward).
        {
            std::vector<CB>& bucket = pending[s];
            for (std::size_t ci = 0u; ci < bucket.size() && ok; ++ci) {
                CB& cb = bucket[ci];
                const std::size_t cnr = cb.rows.size();
                const std::size_t cnc = cb.cols.size();
                fwd_r.clear(); fwd_c.clear(); fwd_v.clear();
                Index fwd_min = Index(-1);
                for (std::size_t jj = 0u; jj < cnc; ++jj) {
                    const Index gc = cb.cols[jj];
                    const Index cl = (static_cast<std::size_t>(gc) < sn)
                        ? col_local[static_cast<std::size_t>(gc)] : Index(-1);
                    for (std::size_t ii = 0u; ii < cnr; ++ii) {
                        const T v = cb.vals[ii + jj * cnr];
                        const Index gr = cb.rows[ii];
                        const Index rl = (static_cast<std::size_t>(gr) < sn)
                            ? row_local[static_cast<std::size_t>(gr)] : Index(-1);
                        if (rl >= Index(0) && cl >= Index(0)) {
                            F[static_cast<std::size_t>(rl) +
                              static_cast<std::size_t>(cl) * ldF] += v;
                        } else {
                            fwd_r.push_back(gr); fwd_c.push_back(gc); fwd_v.push_back(v);
                            const Index mn = (gr < gc) ? gr : gc;
                            if (fwd_min < Index(0) || mn < fwd_min) fwd_min = mn;
                        }
                    }
                }
                if (!fwd_v.empty()) {
                    const Index tgt = (fwd_min >= Index(0) &&
                                       static_cast<std::size_t>(fwd_min) < sn)
                        ? col_to_supernode[static_cast<std::size_t>(fwd_min)] : Index(-1);
                    if (tgt > static_cast<Index>(s) &&
                        static_cast<std::size_t>(tgt) < nsup) {
                        gb_rows.assign(fwd_r.begin(), fwd_r.end());
                        std::sort(gb_rows.begin(), gb_rows.end());
                        gb_rows.erase(std::unique(gb_rows.begin(), gb_rows.end()), gb_rows.end());
                        gb_cols.assign(fwd_c.begin(), fwd_c.end());
                        std::sort(gb_cols.begin(), gb_cols.end());
                        gb_cols.erase(std::unique(gb_cols.begin(), gb_cols.end()), gb_cols.end());
                        CB blk;
                        blk.rows = gb_rows; blk.cols = gb_cols;
                        const std::size_t bnr = blk.rows.size();
                        blk.vals.assign(bnr * blk.cols.size(), T(0));
                        for (std::size_t e = 0u; e < fwd_v.size(); ++e) {
                            const std::size_t ri = static_cast<std::size_t>(
                                std::lower_bound(blk.rows.begin(), blk.rows.end(), fwd_r[e]) - blk.rows.begin());
                            const std::size_t cj = static_cast<std::size_t>(
                                std::lower_bound(blk.cols.begin(), blk.cols.end(), fwd_c[e]) - blk.cols.begin());
                            blk.vals[ri + cj * bnr] += fwd_v[e];
                        }
                        pending[static_cast<std::size_t>(tgt)].push_back(std::move(blk));
                    } else {
                        ok = false; // unexpected backward routing -> safe fallback
                    }
                }
            }
            std::vector<CB>().swap(pending[s]); // free absorbed blocks
        }
        stats.mf_extend_add_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_ea0).count());
        stats.panel_nonkernel_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_asm0).count());
        if (!ok) break;

        // ---- dense partial factorization of the front ----
        T* Fp = F.data();

        // (1) factor pivot block A11 = F[0:w,0:w] (static diagonal pivoting).
        //     Trailing block update via ger adapter (contiguous pivot-row scratch).
        const auto t_blk0 = std::chrono::steady_clock::now();
        std::size_t blk_kernel_ticks = 0u;
        for (std::size_t k = 0u; k < sw && ok; ++k) {
            const T piv = Fp[k + k * ldF];
            const real_type apv = scalar_pol::abs_value(piv);
            // SLU-GT1 D3 (certified-only): reject the front unless the pivot
            // magnitude is certifiably positive (division by piv follows).
            if (!(apv > opt.zero_tolerance)) {
                ok = false; break;
            }
            // scale L column below the pivot (within block).
            const std::size_t below = sw - k - 1u;
            for (std::size_t i = k + 1u; i < sw; ++i)
                Fp[i + k * ldF] = Fp[i + k * ldF] / piv;
            if (below > 0u) {
                row_tmp.assign(below, T(0));
                for (std::size_t c = 0u; c < below; ++c)
                    row_tmp[c] = Fp[k + (k + 1u + c) * ldF];
                const auto tg0 = std::chrono::steady_clock::now();
                sparse_lu_dense_kernel<T>::ger(
                    below, below,
                    T(-1),
                    &Fp[(k + 1u) + k * ldF], row_tmp.data(),
                    &Fp[(k + 1u) + (k + 1u) * ldF], ldF);
                blk_kernel_ticks += static_cast<std::size_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - tg0).count());
                ++stats.ger_count;
                stats.flop_ger += 2.0 * static_cast<double>(below) *
                                  static_cast<double>(below);
            }
        }
        if (!ok) break;
        stats.dense_kernel_ticks += blk_kernel_ticks;
        stats.within_panel_dense_kernel_ticks += blk_kernel_ticks;
        {
            const std::size_t blk_total = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t_blk0).count());
            stats.within_nonkernel_ticks +=
                (blk_total > blk_kernel_ticks) ? (blk_total - blk_kernel_ticks) : 0u;
        }

        const std::size_t m_off = R - sw; // off-diagonal L rows

        // (2) A12 := L11^{-1} A12  (trsm, unit-lower).
        if (nu > 0u) {
            const auto t0 = std::chrono::steady_clock::now();
            sparse_lu_dense_kernel<T>::trsm(
                'L', 'L', 'N', 'U',
                sw, nu, T(1), Fp, ldF, &Fp[0 + sw * ldF], ldF);
            stats.dense_kernel_ticks += static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t0).count());
            ++stats.trsm_count;
            stats.flop_trsm += static_cast<double>(sw) * static_cast<double>(sw) *
                               static_cast<double>(nu);
        }

        // (3) A21 := A21 U11^{-1}  (trsm, upper, side R).
        if (m_off > 0u) {
            const auto t0 = std::chrono::steady_clock::now();
            sparse_lu_dense_kernel<T>::trsm(
                'R', 'U', 'N', 'N',
                m_off, sw, T(1), Fp, ldF, &Fp[sw + 0 * ldF], ldF);
            stats.dense_kernel_ticks += static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t0).count());
            ++stats.trsm_count;
            stats.flop_trsm += static_cast<double>(m_off) * static_cast<double>(sw) *
                               static_cast<double>(sw);
        }

        // (4) A22 := A22 - A21 * A12  (gemm -- the dominant BLAS-3 flop).
        if (m_off > 0u && nu > 0u) {
            const auto t0 = std::chrono::steady_clock::now();
            sparse_lu_dense_kernel<T>::gemm(
                m_off, nu, sw,
                T(-1), &Fp[sw + 0 * ldF], ldF,
                &Fp[0 + sw * ldF], ldF,
                T(1), &Fp[sw + sw * ldF], ldF);
            const std::size_t gtk = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t0).count());
            stats.dense_kernel_ticks          += gtk;
            stats.panel_update_dense_kernel_ticks += gtk; // extend-add / Schur portion
            ++stats.gemm_count;
            stats.flop_gemm += 2.0 * static_cast<double>(m_off) *
                               static_cast<double>(nu) * static_cast<double>(sw);
            stats.gemm_m_sum += m_off;
            stats.gemm_n_sum += nu;
            stats.gemm_k_sum += sw;
            if (m_off > stats.gemm_max_m) stats.gemm_max_m = m_off;
            if (nu    > stats.gemm_max_n) stats.gemm_max_n = nu;
            if (sw    > stats.gemm_max_k) stats.gemm_max_k = sw;
            std::size_t mind = m_off;
            if (nu < mind) mind = nu;
            if (sw < mind) mind = sw;
            ++stats.gemm_dim_hist[sparse_lu_detail::gemm_dim_bucket(mind)];
        }

        // ---- scatter results into supernodal storage ----
        const auto t_sc0 = std::chrono::steady_clock::now();
        const std::size_t off = static_cast<std::size_t>(desc.values_offset);
        const std::size_t ld  = static_cast<std::size_t>(desc.leading_dimension);

        // panel block (pivot block L\U + off-diagonal L): rows 0..R-1, cols 0..w-1.
        if (off + ld * sw <= storage.panel_values.size() && ld >= R) {
            T* panel = &storage.panel_values[off];
            for (std::size_t c = 0u; c < sw; ++c)
                for (std::size_t r = 0u; r < R; ++r)
                    panel[c * ld + r] = Fp[r + c * ldF];
        } else {
            ok = false; }

        // U-right values -> U_segments (one slot per (pivot row p, col c')).
        if (ok) {
            for (std::size_t k = 0u; k < sw; ++k) {
                const Index p = fc + static_cast<Index>(k);
                if (p < Index(0) || static_cast<std::size_t>(p) >= sn) continue;
                const std::size_t a = urow.ptr[static_cast<std::size_t>(p)];
                const std::size_t b = urow.ptr[static_cast<std::size_t>(p) + 1u];
                for (std::size_t e = a; e < b; ++e) {
                    const Index gc = urow.col[e];
                    const Index cl = col_local[static_cast<std::size_t>(gc)];
                    if (cl < Index(0)) {
                        ok = false; break; }
                    const std::size_t aidx = urow.abs_idx[e];
                    if (aidx >= storage.U_segments.values.size()) {
                        ok = false; break; }
                    storage.U_segments.values[aidx] =
                        Fp[k + static_cast<std::size_t>(cl) * ldF];
                }
                if (!ok) break;
            }
        }
        stats.init_scatter_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_sc0).count());

        // ---- emit contribution block (Schur complement A22) ----
        // One dense block (rows = off-diagonal rows O_s, cols = U-right) routed
        // to the ancestor front that first eliminates the smallest index it
        // touches: supernode(min over entries of min(row,col)) =
        // supernode(min(O_s[0], uright[0])). Children of that front add into one
        // dense matrix; entries that do not belong there are forwarded onward in
        // the absorption step above. This block extend-add (O(CB)) replaces the
        // left-looking per-(k,j) scatter.
        if (ok && m_off > 0u && nu > 0u) {
            const auto t_cb0 = std::chrono::steady_clock::now();
            CB cb;
            cb.rows.assign(desc.row_indices.begin() + static_cast<std::ptrdiff_t>(sw),
                           desc.row_indices.end());
            cb.cols = uright;
            cb.vals.assign(m_off * nu, T(0));
            for (std::size_t jj = 0u; jj < nu; ++jj)
                for (std::size_t ii = 0u; ii < m_off; ++ii)
                    cb.vals[ii + jj * m_off] = Fp[(sw + ii) + (sw + jj) * ldF];
            // smallest per-entry min(row,col) == min(first off-diag row, first
            // U-right col), since both index sets are sorted ascending.
            const Index r0 = cb.rows.front();
            const Index c0 = cb.cols.front();
            const Index key = (r0 < c0) ? r0 : c0;
            const Index tgt = (static_cast<std::size_t>(key) < sn)
                ? col_to_supernode[static_cast<std::size_t>(key)] : Index(-1);
            if (tgt > static_cast<Index>(s) && static_cast<std::size_t>(tgt) < nsup) {
                pending[static_cast<std::size_t>(tgt)].push_back(std::move(cb));
            } else {
                ok = false;
            }
            const std::size_t cb_dt = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t_cb0).count());
            stats.mf_cb_emit_ticks      += cb_dt;
            stats.panel_nonkernel_ticks += cb_dt;
        }

        // ---- clear front local maps ----
        for (std::size_t e = 0u; e < touched_rows.size(); ++e)
            row_local[static_cast<std::size_t>(touched_rows[e])] = Index(-1);
        for (std::size_t e = 0u; e < touched_cols.size(); ++e)
            col_local[static_cast<std::size_t>(touched_cols[e])] = Index(-1);

        stats.within_panel_count++;
    }

    if (!ok) {
        stats.supernodes_processed = stats.within_panel_count;
        stats.status = supernodal_true_numeric_status::pivot_failure;
        finish_total_ticks();
        return stats;
    }

    stats.supernodes_processed = nsup;

    // ------------------------------------------------------------------
    // Step 2: factorization residual check (reuse the scalable sparse path).
    //
    // SLU-MF7: opt-in residual gate (opt.supernodal_native_check_residual).
    // DEFAULT false = skip (accepted on structural success; avoids the L*U
    // reconstruction overhead on narrow fronts).  Set true to run the gate
    // byte-identical to the original: accept only when rel_res or abs_res <= 1e-6,
    // else fall back to GP.
    // ------------------------------------------------------------------
    if (opt.supernodal_native_check_residual) {
        real_type abs_res(0), rel_res(0);
        const auto t_res0 = std::chrono::steady_clock::now();
        bool checked = sparse_lu_detail::compute_supernodal_lu_residual_sparse(
            A_csc, storage, n, abs_res, rel_res);
        stats.residual_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_res0).count());

        stats.factorization_residual_abs = abs_res;
        stats.factorization_residual_rel = rel_res;

        if (!checked) {
            stats.factorization_residual_checked = false;
            stats.factorization_residual_passed  = false;
            stats.status = supernodal_true_numeric_status::residual_not_checked;
            finish_total_ticks();
            return stats;
        }
        stats.factorization_residual_checked = true;
        // Requirement-set arithmetic sanity threshold (certainly-<=; not a
        // rigorous error bound).  SLU-GT1 D2.
        const real_type tol = real_type(1e-6);
        stats.factorization_residual_passed = (rel_res <= tol) || (abs_res <= tol);
        if (!stats.factorization_residual_passed) {
            stats.status = supernodal_true_numeric_status::residual_failed;
            finish_total_ticks();
            return stats;
        }
    } else {
        // Residual gate skipped by opt-in: accept on structural success.
        stats.factorization_residual_checked = false;
        stats.factorization_residual_passed  = true;
    }

    // ------------------------------------------------------------------
    // Step 3: accept the multifrontal A_eff-origin numeric source.
    // ------------------------------------------------------------------
    storage.true_numeric_source                 = true;
    storage.values_initialized_from_A           = true;
    storage.values_initialized_from_csc_numeric = false;
    storage.numeric_source_kind =
        supernodal_numeric_source_kind::a_eff_true_numeric;

    stats.success = true;
    stats.status  = supernodal_true_numeric_status::success;
    finish_total_ticks();
    return stats;
}

// ===========================================================================
// sparse_lu_factorize_supernodal_multifrontal_inplace  (SLU-MF3)
//
// In-place frontal-update variant of the multifrontal driver. Instead of the
// MF2 "pull" model (each front emits a COMPACT contribution block, buckets it,
// and the parent later replays it with a scatter extend-add), this driver uses
// a "push" model: each front owns a persistent dense buffer Fbuf[s], and the
// moment a child finishes its dense factorization it scatters its Schur
// complement A22 DIRECTLY into its parent's buffer Fbuf[t] (in place), then
// frees its own buffer.
//
// Effect on the MF2 numeric-time breakdown (measured dominant costs):
//   - cb_emit (Schur copy-out into a fresh compact block): ELIMINATED. The
//     child reads A22 straight from its own front buffer (strided) and adds it
//     into the parent -- no intermediate compact buffer, no re-expansion.
//   - extend_add: now a single direct push-scatter (read child A22, add into
//     parent), instead of compact-write + compact-read + scatter (~2/3 the
//     contribution-path memory traffic of MF2).
//   - forwarding-replay machinery: removed. The self-symbolic structure
//     guarantees the routed target front contains every entry of O_s x Ur_s
//     (build_supernodal_self_symbolic_storage adds them), so a self-symbolic
//     structure never forwards. If any entry fails to map (non-self-symbolic
//     input, or structural inconsistency), the driver safely fails (ok=false)
//     and the caller keeps the CSC fallback.
//
// Numerically equivalent to the MF2 driver up to addition ORDER (children
// accumulate into the parent as they finish rather than being replayed at the
// parent), which is an O(eps) reordering; the post-factorization residual check
// remains the final acceptance gate. Dense kernels, FLOP, and front shapes are
// identical to MF2. Output satisfies the same supernodal_lu_storage contract
// (§18.2 solve unchanged).
//
// Same signature / return contract as sparse_lu_factorize_supernodal_multifrontal.
// ===========================================================================
template <class T, class Index>
sparse_lu_detail::supernodal_true_numeric_stats<
    typename vcp::tsparse_scalar::real_type<T>::type>
sparse_lu_factorize_supernodal_multifrontal_inplace(
    const csc_storage<T, Index>&           A_csc,
    const baseline_lu_storage<T, Index>&   csc_lu,
    supernodal_lu_storage<T, Index>&       storage,
    const sparse_lu_options<T>&            opt)
{
    (void)csc_lu; // row_perm / Dr / Dc already in storage (bootstrap copy)

    typedef sparse_lu_scalar_policy<T>      scalar_pol;
    typedef typename scalar_pol::real_type  real_type;

    const auto t_total_start = std::chrono::steady_clock::now();

    sparse_lu_detail::supernodal_true_numeric_stats<real_type> stats;
    stats.attempted = true;

    const auto finish_total_ticks = [&]() {
        stats.total_ticks = static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_total_start).count());
        stats.factorization_ticks =
            (stats.total_ticks > stats.residual_ticks)
                ? (stats.total_ticks - stats.residual_ticks) : 0u;
    };

    if (!storage.valid || storage.supernodes.empty()) {
        stats.status = supernodal_true_numeric_status::unsupported_structure;
        finish_total_ticks();
        return stats;
    }

    const std::size_t nsup = storage.supernodes.size();
    const Index n = static_cast<Index>(storage.row_perm.size());
    const std::size_t sn = (n > Index(0)) ? static_cast<std::size_t>(n) : 0u;

    if (A_csc.col_ptr.size() < sn + 1u) {
        stats.status = supernodal_true_numeric_status::unsupported_structure;
        finish_total_ticks();
        return stats;
    }

    std::fill(storage.panel_values.begin(),      storage.panel_values.end(),      T(0));
    std::fill(storage.U_segments.values.begin(), storage.U_segments.values.end(), T(0));
    stats.values_initialized_from_A           = true;
    stats.values_initialized_from_csc_numeric = false;

    const auto t_sym0 = std::chrono::steady_clock::now();

    // col -> supernode map.
    std::vector<Index> col_to_supernode(sn, Index(-1));
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& d = storage.supernodes[s];
        const Index cb = d.first_col;
        const Index ce = cb + d.num_cols;
        for (Index c = cb; c < ce; ++c)
            if (c >= Index(0) && static_cast<std::size_t>(c) < sn)
                col_to_supernode[static_cast<std::size_t>(c)] = static_cast<Index>(s);
    }

    // inverse row permutation.
    std::vector<Index> inv_row_perm(sn, Index(-1));
    for (std::size_t nr = 0u; nr < storage.row_perm.size() && nr < sn; ++nr) {
        const Index old_row = storage.row_perm[nr];
        if (old_row >= Index(0) && static_cast<std::size_t>(old_row) < sn)
            inv_row_perm[static_cast<std::size_t>(old_row)] = static_cast<Index>(nr);
    }

    // transpose-of-U index (scatter targets for U-right values).
    sparse_lu_detail::mf_urow_index<Index> urow;
    sparse_lu_detail::mf_build_urow_index(storage, n, urow);

    // Per-front U-right column set, ncols, and leading dimension. Precomputed so
    // a child can size and address its parent's buffer without re-deriving it.
    //
    // SLU-MF6: the self-symbolic builder already computes each front's U-right
    // column set (Ur_s) and now exports it as storage.mf_uright_of. Reuse it by
    // reference -- this removes the driver's single most expensive symbolic
    // segment (the per-front gather + std::sort + std::unique over the U_segments
    // transpose, ~70% of symbolic_ticks). The exported sets are byte-identical to
    // the rebuild (both sorted ascending, unique). A fallback rebuild from urow is
    // retained for any source_of_truth storage that did not export them.
    std::vector<std::vector<Index> > uright_fallback;
    const bool have_uright_export = (storage.mf_uright_of.size() == nsup);
    if (!have_uright_export) {
        uright_fallback.resize(nsup);
        for (std::size_t s = 0u; s < nsup; ++s) {
            const supernode_desc<Index>& desc = storage.supernodes[s];
            const Index fc = desc.first_col;
            const Index w  = desc.num_cols;
            const std::size_t sw = (w > Index(0)) ? static_cast<std::size_t>(w) : 0u;
            std::vector<Index>& ur = uright_fallback[s];
            for (std::size_t k = 0u; k < sw; ++k) {
                const Index p = fc + static_cast<Index>(k);
                if (p < Index(0) || static_cast<std::size_t>(p) >= sn) continue;
                const std::size_t a = urow.ptr[static_cast<std::size_t>(p)];
                const std::size_t b = urow.ptr[static_cast<std::size_t>(p) + 1u];
                for (std::size_t e = a; e < b; ++e) ur.push_back(urow.col[e]);
            }
            std::sort(ur.begin(), ur.end());
            ur.erase(std::unique(ur.begin(), ur.end()), ur.end());
        }
    }
    const std::vector<std::vector<Index> >& uright_of =
        have_uright_export ? storage.mf_uright_of : uright_fallback;

    std::vector<std::size_t> ncols_of(nsup, 0u);
    std::vector<std::size_t> ld_of(nsup, 0u);
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& desc = storage.supernodes[s];
        const Index w  = desc.num_cols;
        const std::size_t sw = (w > Index(0)) ? static_cast<std::size_t>(w) : 0u;
        ld_of[s]    = desc.row_indices.size();
        ncols_of[s] = sw + uright_of[s].size();
    }

    // SLU-MF5: precompute child->parent relative-index maps (rrel runs + crel)
    // once, here, so the numeric push is a pure direct add (no std::lower_bound).
    std::vector<sparse_lu_detail::mf_relative_map<Index> > relmap;
    bool relmap_coverage_ok = true;
    sparse_lu_detail::mf_build_relative_maps(
        storage, uright_of, ld_of, col_to_supernode, n,
        relmap, relmap_coverage_ok);

    stats.symbolic_ticks += static_cast<std::size_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t_sym0).count());

    // Persistent per-front dense buffers (column-major, ld_of[s] x ncols_of[s]).
    // Allocated lazily (the first time a child pushes into it, or when the front
    // itself is reached), and freed immediately after the front is factored.
    std::vector<std::vector<T> > Fbuf(nsup);
    std::vector<char>            Falloc(nsup, 0);
    const auto ensure_front = [&](std::size_t t) -> bool {
        if (Falloc[t]) return true;
        const std::size_t need = ld_of[t] * ncols_of[t];
        if (need == 0u) return false;
        Fbuf[t].assign(need, T(0));
        Falloc[t] = 1;
        return true;
    };

    // Reusable global local-index maps (cleared via touched lists).
    std::vector<Index> row_local(sn, Index(-1));
    std::vector<Index> col_local(sn, Index(-1));
    std::vector<Index> touched_rows, touched_cols;
    std::vector<T>     row_tmp;       // contiguous pivot-row scratch for ger

    bool ok = true;

    for (std::size_t s = 0u; s < nsup && ok; ++s) {
        const supernode_desc<Index>& desc = storage.supernodes[s];
        const Index fc = desc.first_col;
        const Index w  = desc.num_cols;
        const std::size_t R  = desc.row_indices.size();
        const std::size_t sw = static_cast<std::size_t>(w);

        if (w <= Index(0) || R == 0u) {
            if (Falloc[s]) { ok = false; break; } // pending push but nothing to absorb
            continue;
        }

        const std::vector<Index>& uright = uright_of[s];
        const std::size_t nu    = uright.size();
        const std::size_t ncols = ncols_of[s];
        const std::size_t ldF   = ld_of[s];

        if (!ensure_front(s)) { ok = false; break; }
        T* Fp = Fbuf[s].data();

        // ---- front local index maps ----
        touched_rows.clear();
        for (std::size_t r = 0u; r < R; ++r) {
            const Index gr = desc.row_indices[r];
            if (gr < Index(0) || static_cast<std::size_t>(gr) >= sn) { ok = false; break; }
            row_local[static_cast<std::size_t>(gr)] = static_cast<Index>(r);
            touched_rows.push_back(gr);
        }
        if (!ok) break;
        touched_cols.clear();
        for (std::size_t k = 0u; k < sw; ++k) {
            const Index gc = fc + static_cast<Index>(k);
            col_local[static_cast<std::size_t>(gc)] = static_cast<Index>(k);
            touched_cols.push_back(gc);
        }
        for (std::size_t e = 0u; e < nu; ++e) {
            const Index gc = uright[e];
            col_local[static_cast<std::size_t>(gc)] = static_cast<Index>(sw + e);
            touched_cols.push_back(gc);
        }

        // ---- (a) scatter A_eff into the front buffer (already holds children's
        //          pushed Schur contributions). Pivot columns -> all front rows;
        //          U-right columns -> pivot rows only. ----
        const auto t_asm0 = std::chrono::steady_clock::now();
        for (std::size_t cl = 0u; cl < ncols; ++cl) {
            const Index gc = touched_cols[cl];
            const bool pivot_col = (cl < sw);
            const std::size_t scgc = static_cast<std::size_t>(gc);
            const Index cs = A_csc.col_ptr[scgc];
            const Index ce = A_csc.col_ptr[scgc + 1u];
            const T dc = (!storage.Dc.empty() && scgc < storage.Dc.size())
                ? storage.Dc[scgc] : T(1);
            for (Index idx = cs; idx < ce; ++idx) {
                const std::size_t sidx = static_cast<std::size_t>(idx);
                if (sidx >= A_csc.row_ind.size()) continue;
                const Index old_row = A_csc.row_ind[sidx];
                if (old_row < Index(0) || static_cast<std::size_t>(old_row) >= sn) continue;
                const Index new_row = inv_row_perm[static_cast<std::size_t>(old_row)];
                if (new_row < Index(0)) continue;
                const Index rl = row_local[static_cast<std::size_t>(new_row)];
                if (rl < Index(0)) continue;
                if (!pivot_col && static_cast<std::size_t>(rl) >= sw) continue;
                const T dr = (!storage.Dr.empty() &&
                              static_cast<std::size_t>(new_row) < storage.Dr.size())
                    ? storage.Dr[static_cast<std::size_t>(new_row)] : T(1);
                Fp[static_cast<std::size_t>(rl) + cl * ldF] +=
                    dr * A_csc.values[sidx] * dc;
            }
        }
        stats.mf_aeff_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_asm0).count());
        stats.panel_nonkernel_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_asm0).count());

        // ---- dense partial factorization of the front (identical to MF2) ----
        const auto t_blk0 = std::chrono::steady_clock::now();
        std::size_t blk_kernel_ticks = 0u;
        for (std::size_t k = 0u; k < sw && ok; ++k) {
            const T piv = Fp[k + k * ldF];
            const real_type apv = scalar_pol::abs_value(piv);
            // SLU-GT1 D3 (certified-only): see the MF2 block above.
            if (!(apv > opt.zero_tolerance)) { ok = false; break; }
            for (std::size_t i = k + 1u; i < sw; ++i)
                Fp[i + k * ldF] = Fp[i + k * ldF] / piv;
            const std::size_t below = sw - k - 1u;
            if (below > 0u) {
                row_tmp.assign(below, T(0));
                for (std::size_t c = 0u; c < below; ++c)
                    row_tmp[c] = Fp[k + (k + 1u + c) * ldF];
                const auto tg0 = std::chrono::steady_clock::now();
                sparse_lu_dense_kernel<T>::ger(
                    below, below, T(-1),
                    &Fp[(k + 1u) + k * ldF], row_tmp.data(),
                    &Fp[(k + 1u) + (k + 1u) * ldF], ldF);
                blk_kernel_ticks += static_cast<std::size_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - tg0).count());
                ++stats.ger_count;
                stats.flop_ger += 2.0 * static_cast<double>(below) *
                                  static_cast<double>(below);
            }
        }
        if (!ok) break;
        stats.dense_kernel_ticks += blk_kernel_ticks;
        stats.within_panel_dense_kernel_ticks += blk_kernel_ticks;
        {
            const std::size_t blk_total = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t_blk0).count());
            stats.within_nonkernel_ticks +=
                (blk_total > blk_kernel_ticks) ? (blk_total - blk_kernel_ticks) : 0u;
        }

        const std::size_t m_off = R - sw;

        if (nu > 0u) {
            const auto t0 = std::chrono::steady_clock::now();
            sparse_lu_dense_kernel<T>::trsm(
                'L', 'L', 'N', 'U', sw, nu, T(1), Fp, ldF, &Fp[0 + sw * ldF], ldF);
            stats.dense_kernel_ticks += static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t0).count());
            ++stats.trsm_count;
            stats.flop_trsm += static_cast<double>(sw) * static_cast<double>(sw) *
                               static_cast<double>(nu);
        }
        if (m_off > 0u) {
            const auto t0 = std::chrono::steady_clock::now();
            sparse_lu_dense_kernel<T>::trsm(
                'R', 'U', 'N', 'N', m_off, sw, T(1), Fp, ldF, &Fp[sw + 0 * ldF], ldF);
            stats.dense_kernel_ticks += static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t0).count());
            ++stats.trsm_count;
            stats.flop_trsm += static_cast<double>(m_off) * static_cast<double>(sw) *
                               static_cast<double>(sw);
        }
        if (m_off > 0u && nu > 0u) {
            const auto t0 = std::chrono::steady_clock::now();
            sparse_lu_dense_kernel<T>::gemm(
                m_off, nu, sw,
                T(-1), &Fp[sw + 0 * ldF], ldF, &Fp[0 + sw * ldF], ldF,
                T(1), &Fp[sw + sw * ldF], ldF);
            const std::size_t gtk = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t0).count());
            stats.dense_kernel_ticks              += gtk;
            stats.panel_update_dense_kernel_ticks += gtk;
            ++stats.gemm_count;
            stats.flop_gemm += 2.0 * static_cast<double>(m_off) *
                               static_cast<double>(nu) * static_cast<double>(sw);
            stats.gemm_m_sum += m_off;
            stats.gemm_n_sum += nu;
            stats.gemm_k_sum += sw;
            if (m_off > stats.gemm_max_m) stats.gemm_max_m = m_off;
            if (nu    > stats.gemm_max_n) stats.gemm_max_n = nu;
            if (sw    > stats.gemm_max_k) stats.gemm_max_k = sw;
            std::size_t mind = m_off;
            if (nu < mind) mind = nu;
            if (sw < mind) mind = sw;
            ++stats.gemm_dim_hist[sparse_lu_detail::gemm_dim_bucket(mind)];
        }

        // ---- scatter results into supernodal storage ----
        const auto t_sc0 = std::chrono::steady_clock::now();
        const std::size_t off = static_cast<std::size_t>(desc.values_offset);
        const std::size_t ld  = static_cast<std::size_t>(desc.leading_dimension);
        if (off + ld * sw <= storage.panel_values.size() && ld >= R) {
            T* panel = &storage.panel_values[off];
            for (std::size_t c = 0u; c < sw; ++c)
                for (std::size_t r = 0u; r < R; ++r)
                    panel[c * ld + r] = Fp[r + c * ldF];
        } else { ok = false; }
        if (ok) {
            for (std::size_t k = 0u; k < sw; ++k) {
                const Index p = fc + static_cast<Index>(k);
                if (p < Index(0) || static_cast<std::size_t>(p) >= sn) continue;
                const std::size_t a = urow.ptr[static_cast<std::size_t>(p)];
                const std::size_t b = urow.ptr[static_cast<std::size_t>(p) + 1u];
                for (std::size_t e = a; e < b; ++e) {
                    const Index gc = urow.col[e];
                    const Index cl = col_local[static_cast<std::size_t>(gc)];
                    if (cl < Index(0)) { ok = false; break; }
                    const std::size_t aidx = urow.abs_idx[e];
                    if (aidx >= storage.U_segments.values.size()) { ok = false; break; }
                    storage.U_segments.values[aidx] =
                        Fp[k + static_cast<std::size_t>(cl) * ldF];
                }
                if (!ok) break;
            }
        }
        stats.init_scatter_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_sc0).count());
        if (!ok) break;

        // ---- in-place push of the Schur complement A22 into the parent front ----
        // A22 lives in Fbuf[s] at rows [sw,R) x cols [sw,ncols) (strided by ldF).
        // Route to t = supernode(min(first off-diag row, first U-right col)); the
        // self-symbolic structure guarantees t's row_indices/cols contain every
        // (O_s row, U-right col) pair, so the push always lands -- no forwarding.
        if (m_off > 0u && nu > 0u) {
            const auto t_push0 = std::chrono::steady_clock::now();
            // SLU-MF5: use the symbolic-precomputed relative-index map. The
            // numeric push is now a pure direct add over contiguous rrel runs --
            // NO std::lower_bound in this loop. Routing target / coverage match
            // the MF3 push exactly, so the result is byte-identical.
            const sparse_lu_detail::mf_relative_map<Index>& rm = relmap[s];
            const Index tgt = rm.parent;
            if (tgt > static_cast<Index>(s) && static_cast<std::size_t>(tgt) < nsup &&
                rm.coverage_ok && ensure_front(static_cast<std::size_t>(tgt))) {
                const std::size_t t   = static_cast<std::size_t>(tgt);
                const std::size_t tld = ld_of[t];
                T* Ftp = Fbuf[t].data();

                const Index* crel    = rm.crel.data();
                const Index* run_src = rm.run_src.data();
                const Index* run_dst = rm.run_dst.data();
                const Index* run_len = rm.run_len.data();
                const std::size_t nruns = rm.run_len.size();

                for (std::size_t j = 0u; j < nu; ++j) {
                    const std::size_t pc = static_cast<std::size_t>(crel[j]) * tld;
                    const std::size_t sc = (sw + j) * ldF + sw;
                    T* col = Ftp + pc;
                    const T* src = Fp + sc;
                    for (std::size_t rr = 0u; rr < nruns; ++rr) {
                        T*       dst = col + static_cast<std::size_t>(run_dst[rr]);
                        const T* s2  = src + static_cast<std::size_t>(run_src[rr]);
                        const std::size_t len = static_cast<std::size_t>(run_len[rr]);
                        for (std::size_t l = 0u; l < len; ++l) dst[l] += s2[l];
                    }
                }
            } else {
                ok = false; // backward / out-of-range routing / coverage miss
            }
            const std::size_t push_dt = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t_push0).count());
            stats.mf_extend_add_ticks   += push_dt;
            stats.panel_nonkernel_ticks += push_dt;
        }

        // ---- release this front's buffer; clear local maps ----
        std::vector<T>().swap(Fbuf[s]);
        Falloc[s] = 0;
        for (std::size_t e = 0u; e < touched_rows.size(); ++e)
            row_local[static_cast<std::size_t>(touched_rows[e])] = Index(-1);
        for (std::size_t e = 0u; e < touched_cols.size(); ++e)
            col_local[static_cast<std::size_t>(touched_cols[e])] = Index(-1);

        stats.within_panel_count++;
    }

    if (!ok) {
        stats.supernodes_processed = stats.within_panel_count;
        stats.status = supernodal_true_numeric_status::pivot_failure;
        finish_total_ticks();
        return stats;
    }
    stats.supernodes_processed = nsup;

    // ---- residual acceptance gate (reuse the scalable sparse path) ----
    // SLU-MF7: opt-in (opt.supernodal_native_check_residual).  DEFAULT false =
    // skip (structural acceptance).  Set true for GP-fallback safety net.
    if (opt.supernodal_native_check_residual) {
        real_type abs_res(0), rel_res(0);
        const auto t_res0 = std::chrono::steady_clock::now();
        bool checked = sparse_lu_detail::compute_supernodal_lu_residual_sparse(
            A_csc, storage, n, abs_res, rel_res);
        stats.residual_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_res0).count());
        stats.factorization_residual_abs = abs_res;
        stats.factorization_residual_rel = rel_res;
        if (!checked) {
            stats.factorization_residual_checked = false;
            stats.factorization_residual_passed  = false;
            stats.status = supernodal_true_numeric_status::residual_not_checked;
            finish_total_ticks();
            return stats;
        }
        stats.factorization_residual_checked = true;
        // Requirement-set arithmetic sanity threshold (certainly-<=; not a
        // rigorous error bound).  SLU-GT1 D2.
        const real_type tol = real_type(1e-6);
        stats.factorization_residual_passed = (rel_res <= tol) || (abs_res <= tol);
        if (!stats.factorization_residual_passed) {
            stats.status = supernodal_true_numeric_status::residual_failed;
            finish_total_ticks();
            return stats;
        }
    } else {
        stats.factorization_residual_checked = false;
        stats.factorization_residual_passed  = true;
    }

    storage.true_numeric_source                 = true;
    storage.values_initialized_from_A           = true;
    storage.values_initialized_from_csc_numeric = false;
    storage.numeric_source_kind =
        supernodal_numeric_source_kind::a_eff_true_numeric;

    stats.success = true;
    stats.status  = supernodal_true_numeric_status::success;
    finish_total_ticks();
    return stats;
}

// ===========================================================================
// sparse_lu_factorize_supernodal_numeric_source  (public dispatcher)
//
// Production numeric source for method=supernodal. Selects the multifrontal
// driver (default) or the retained left-looking diagnostic driver
// (opt.supernodal_numeric_diagnostic_leftlooking == true).
// ===========================================================================
template <class T, class Index>
sparse_lu_detail::supernodal_true_numeric_stats<
    typename vcp::tsparse_scalar::real_type<T>::type>
sparse_lu_factorize_supernodal_numeric_source(
    const csc_storage<T, Index>&           A_csc,
    const baseline_lu_storage<T, Index>&   csc_lu,
    supernodal_lu_storage<T, Index>&       storage,
    const sparse_lu_options<T>&            opt)
{
    if (opt.supernodal_numeric_diagnostic_leftlooking) {
        // Diagnostic: retained left-looking A_eff-origin driver (SLU-8R.5.5).
        return sparse_lu_factorize_supernodal_from_a_eff(
            A_csc, csc_lu, storage, opt);
    }
    if (opt.supernodal_inplace_frontal && storage.source_of_truth_storage) {
        // SLU-MF3: in-place frontal-update driver. Requires a self-symbolic
        // (source_of_truth) structure so routing is fully contained.
        return sparse_lu_factorize_supernodal_multifrontal_inplace(
            A_csc, csc_lu, storage, opt);
    }
    return sparse_lu_factorize_supernodal_multifrontal(
        A_csc, csc_lu, storage, opt);
}

#endif // VCP_TSPARSE_SPARSE_LU_MULTIFRONTAL_IMPL_HPP
