// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// SLU-SP1: opt-in left-looking supernode-panel LU factorization
// (design sandbox/docs/design/SLU-SP1_design.md §2; sequential SuperLU
// construction of Demmel et al. 1999 with the U-block treatment of
// Li 2005 §2.3 -- zero-padding + contiguous buffer + one GEMM per
// (supernode, panel) -- as the Phase 3 numeric).
//
// This file MUST be #included from WITHIN namespace vcp, AFTER the SLU-1/2/3
// helpers (csc_storage, baseline_lu_storage, sparse_lu_scalar_policy,
// sparse_lu_is_acceptable_pivot, sparse_lu_identity_permutation,
// sparse_lu_inverse_permutation) are in scope -- same injection contract as
// tsparse_sparse_lu_numeric_impl.hpp.  Do NOT include this file directly.
//
// kv::dd note (design D-4): the panel update below issues its Level-3 work
// as ONE tgemm<T> call per (supernode, panel).  The Ozaki-scheme dd
// specialization lives in <vcp/tblas/tblas_dd.hpp>, which -- following the
// existing opt-in discipline of tlapack_dd.hpp, and because tsparse headers
// are kv-agnostic -- is NOT included here.  A TU that wants the fast dd
// path must include <vcp/tblas/tblas_dd.hpp> before instantiating the
// factorization; otherwise the generic reference tgemm is used (correct,
// slower).  gemv/ger/trsm are NEVER called on this path; the supernode
// diagonal-block solve is scalar by design (width <= maxsup; measured by
// trsm_time_ns for the T-7 / STOP-6 gate).
//
// ---------------------------------------------------------------------------
// Algorithm summary (Phase 3)
//
//   Columns are processed in panels of width w = opt.panel_size over a dense
//   n x w panel workspace W (one SPA column per panel column).  L is held in
//   per-supernode dense trapezoids (supernode = maximal run of consecutive
//   columns with equal -- or, under weak relaxation, nested -- STORED
//   sub-diagonal structure, detected dynamically after each column's pivot).
//   U is emitted per column into CSC buffers with the exact GP store rules.
//
//   Per panel [jcol, jend), w_p = jend - jcol:
//     0. Scatter A[:, jcol..jend) into W via inv_row_perm (all panel
//        columns; later in-panel row swaps are applied to W wholesale,
//        which reproduces GP's scatter-after-swap values exactly).
//     1. PANEL DFS: one traversal of the supernode graph (adjacency =
//        per-supernode row lists with Eisenstat-Liu symmetric pruning) from
//        the union of the panel columns' A-pattern rows < jcol.  Yields
//        R_out, the ascending list of updating supernodes.
//     2. PANEL UPDATE (the §2.2 numeric, per supernode K in R_out):
//        a. per panel column: gather the U fragment W(K_cols, j) -- rows of
//           K's column range are pivotal, hence stable -- into a contiguous
//           nc x w_p buffer, ZERO-PADDED for inactive columns (Li 2005
//           §2.3);
//        b. scalar unit-lower triangular solve with K's diagonal block per
//           active column (D-4: allowed, measured);
//        c. ONE tgemm<T> call: C = L(ext rows, K) * U_frag  (m = |ext|,
//           n = w_p, k = |K|); a scalar fallback handles thin updates
//           (k == 1 or m*k below opt.blas_min_block_size);
//        d. scatter: write the solved fragments back into W(K_cols, j) and
//           subtract C from W(ext rows, j) for the active columns,
//           maintaining the per-column patterns.
//     3. WITHIN-PANEL FACTORIZATION, column by column (scalar; GP order):
//        in-panel updates by ascending scan with the active-column guard,
//        pivot selection / acceptability / row swap (renames applied to the
//        supernode trapezoids through a per-row reverse index, and to the
//        panel workspace rows wholesale), U store, dynamic supernode
//        decision, L store, Eisenstat-Liu pruning with replacement-path
//        verification.
//
//   Final: trapezoids are emitted into CSC L (exact zeros skipped -- the GP
//   store rule), U is already CSC; row/col permutations as in GP.  The
//   result is baseline_lu_storage, so the existing CSC solve / IR / LUX
//   consumers work unchanged.
//
// Correctness notes:
//   * The internal slot graph equals the GP STORED graph: creation/join
//     drop exact zeros exactly like the GP store (keeping active-but-zero
//     rows was measured to inflate the patterns 14x through pivot-swap
//     zeros).  Update application order per column is globally ascending
//     (out-panel supernodes ascending, then in-panel ascending), the GP
//     reach order; only the GEMM's inner accumulation order differs, which
//     is why byte identity with GP is NOT claimed (design §2.3 -- the T-4
//     gates are a 1e-12 solution band and 10x residual parity).
//   * Symmetric pruning keeps rows <= j (pivotal, stable) and additionally
//     KEEPS any row > j whose fill in column j was not actually stored
//     (exact cancellation): the Eisenstat-Liu replacement path runs through
//     column j's supernode and exists only for stored rows.  Intra-supernode
//     U triggers never prune (the replacement path would degenerate into a
//     self-loop at supernode granularity).
//   * Determinism: all data structures and loops are input-determined; the
//     GEMM/scalar dispatch depends only on structure sizes.
// ---------------------------------------------------------------------------

#ifndef VCP_TSPARSE_SPARSE_LU_SUPERNODE_PANEL_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_SUPERNODE_PANEL_IMPL_HPP

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// supernode_panel_factorize_result: output of the supernode-panel numeric.
// storage is the SAME baseline CSC form the GP numeric produces, so every
// existing consumer (solve, IR, LUX extraction) works unchanged.
// ---------------------------------------------------------------------------
template <class T, class Index>
struct supernode_panel_factorize_result {
    bool                          success;
    sparse_lu_status              status;
    baseline_lu_storage<T, Index> storage;

    Index       number_of_supernodes;
    std::size_t panel_count;
    std::size_t prune_event_count;      // diagnostics: pruning fired
    std::size_t panel_dfs_visit_count;  // diagnostics: supernodes visited by panel DFS
    std::size_t axpy_element_count;     // diagnostics: scalar update elements
                                        //   (scalar fallback + in-panel)
    std::size_t slot_element_count;     // diagnostics: total structural column slots

    // Phase 3 (T-7 / STOP-6 evidence):
    std::size_t trsm_time_ns;           // scalar diagonal-block solve time
    std::size_t gemm_time_ns;           // tgemm call time
    std::size_t gemm_call_count;        // number of tgemm calls
    std::size_t gemm_flop_count;        // 2*m*n*k summed over tgemm calls
    std::size_t scalar_update_count;    // supernode updates on the scalar fallback

    supernode_panel_factorize_result()
        : success(false),
          status(sparse_lu_status::numerical_singularity),
          number_of_supernodes(Index(0)),
          panel_count(0u),
          prune_event_count(0u),
          panel_dfs_visit_count(0u),
          axpy_element_count(0u),
          slot_element_count(0u),
          trsm_time_ns(0u),
          gemm_time_ns(0u),
          gemm_call_count(0u),
          gemm_flop_count(0u),
          scalar_update_count(0u) {}
};

// ---------------------------------------------------------------------------
// supernode_panel_lu_factorize (SLU-SP1 Phase 3: zero-padding GEMM numeric)
//
// Preconditions (caller ensures -- the dispatch in tsparse_sparse_lu.hpp):
//   - A_csc is the column-permuted effective matrix (default pivoting path);
//   - opt.pivoting is none / diagonal / threshold_partial (static_mc64 and
//     equilibration are answered not_implemented at the dispatch, D-6);
//   - col_perm / inv_col_perm are valid permutations of size n.
// ---------------------------------------------------------------------------
template <class T, class Index>
supernode_panel_factorize_result<T, Index>
supernode_panel_lu_factorize(
    const csc_storage<T, Index>& A_csc,
    Index n,
    const std::vector<Index>& col_perm,
    const std::vector<Index>& inv_col_perm,
    const sparse_lu_options<T>& opt)
{
    typedef typename sparse_lu_scalar_policy<T>::real_type real_type;
    typedef std::chrono::steady_clock clock_type;

    supernode_panel_factorize_result<T, Index> result;
    const std::size_t un = static_cast<std::size_t>(n);

    // Trivial 0x0 (same contract as the GP numeric)
    if (n == Index(0)) {
        baseline_lu_storage<T, Index>& lu = result.storage;
        lu.L.col_ptr.resize(1u, Index(0));
        lu.U.col_ptr.resize(1u, Index(0));
        lu.row_perm     = col_perm;
        lu.inv_row_perm = inv_col_perm;
        lu.col_perm     = col_perm;
        lu.inv_col_perm = inv_col_perm;
        result.success = true;
        result.status  = sparse_lu_status::success;
        return result;
    }

    const std::size_t panel_w =
        (opt.panel_size > 0u) ? opt.panel_size : 1u;
    const std::size_t maxsup =
        (opt.supernode_panel_maxsup > 0u) ? opt.supernode_panel_maxsup : 1u;
    const std::size_t relax_allow = opt.supernode_relaxation;
    // GEMM dispatch gate (the design's "2D block threshold").  Measured on
    // the reference kernels (2026-07-18, i7-11700):
    //   - double: the GEMM path wins from m*nc >= blas_min_block_size (16)
    //     upward (lap3d16 1.56x vs 1.03x all-scalar);
    //   - kv::dd: every small-block Ozaki tgemm call LOSES to the scalar
    //     update (0.46x at gate 16 -- per-call split/pack overhead), and
    //     tgemm<dd> only matches scalar dd throughput even at large sizes
    //     while tgemm<double> stays at the ~1.5 GFLOP/s reference speed.
    // Until SLU-K1 delivers a fast blocked tgemm<double>, non-fundamental
    // scalar types use a 65536x higher gate (default option value 16 ->
    // m*nc >= 1,048,576): at that setting the Ozaki dd GEMM effectively
    // never fires on present problem sizes -- even the LARGE blocks of
    // lap3d were measured 8% net-slower through tgemm<dd> at the current
    // reference tgemm<double> speed (lap3d16 dd 0.92x with a 16384 gate).
    // The §2.2 GEMM pipeline stays reachable for every T by lowering
    // blas_min_block_size (exercised by slusp1_03's forced-GEMM dd case),
    // and SLU-K1 re-tunes this constant when the kernel speed changes.
    const std::size_t gemm_min_elems =
        std::is_floating_point<T>::value
            ? opt.blas_min_block_size
            : opt.blas_min_block_size * 65536u;

    const real_type threshold        = opt.pivot_threshold;
    const real_type abs_tol          = opt.absolute_pivot_tolerance;
    const bool      do_partial_pivot =
        (opt.pivoting == sparse_lu_pivoting::threshold_partial);

    // ------------------------------------------------------------------
    // Panel workspace: W = n x w dense SPA block (column-major), one mark
    // lane per panel column.  Mark stamps are the GLOBAL column index, so
    // no per-panel reset is needed (columns strictly increase).
    // ------------------------------------------------------------------
    std::vector<T>     W(un * panel_w, T(0));
    std::vector<Index> mark_p(un * panel_w, Index(-1));
    std::vector<std::vector<Index> > patterns(panel_w);
    for (std::size_t jj = 0u; jj < panel_w; ++jj) patterns[jj].reserve(64u);

    std::vector<Index> row_perm     = sparse_lu_identity_permutation(n);
    std::vector<Index> inv_row_perm = sparse_lu_identity_permutation(n);

    // ------------------------------------------------------------------
    // Supernodal working storage for L (see the Phase 2 commit message):
    //   sup_rows[s] : STORED sub-diagonal pattern at creation (pivot-space
    //                 row ids), slot-indexed; renamed in place on row swaps
    //   sup_vals[s] : trapezoidal column-major values; local column c
    //                 occupies [c*ld + c, c*ld + ld), slot t <-> row
    //                 sup_rows[s][t]; NORMALIZED L entries (x / pivot)
    //   adj_slots[s]: DFS adjacency as a permutation of slot ids; only
    //                 [0, prune_len[s]) is traversed
    // ------------------------------------------------------------------
    std::vector<Index>                sup_first;
    std::vector<Index>                sup_ncols;
    std::vector<std::vector<Index> >  sup_rows;
    std::vector<std::vector<T> >      sup_vals;
    std::vector<Index>                col2sup(un, Index(-1));
    std::vector<std::vector<Index> >  adj_slots;
    std::vector<Index>                prune_len;

    struct row_loc { Index sup; Index slot; };
    std::vector<std::vector<row_loc> > row_locs;
    if (do_partial_pivot) {
        row_locs.assign(un, std::vector<row_loc>());
    }

    std::vector<Index> dfs_visited;   // panel DFS generation stamps
    std::vector<Index> u_sup_stamp;   // pruning trigger stamps

    // U CSC buffers (GP layout: explicit diagonal, ascending rows per column)
    std::vector<Index> U_col_ptr(un + 1u, Index(0));
    std::vector<Index> U_row_ind_buf;
    std::vector<T>     U_val_buf;

    // scratch
    std::vector<Index> dfs_stack;
    std::vector<Index> reach_sups;
    std::vector<T>     frag_B;        // nc x w_p zero-padded U fragments
    std::vector<T>     gemm_C;        // m x w_p GEMM result
    std::vector<char>  frag_active;   // per panel column activity

    const Index panel_w_i = static_cast<Index>(panel_w);

    for (Index jcol = Index(0); jcol < n; jcol += panel_w_i) {
        const Index jend = std::min<Index>(n, jcol + panel_w_i);
        const std::size_t w_p = static_cast<std::size_t>(jend - jcol);
        ++result.panel_count;

        // ==============================================================
        // 0. Scatter A[:, jcol..jend) into the panel workspace.
        //    inv_row_perm is the panel-start permutation; later in-panel
        //    swaps move W rows wholesale, reproducing GP's late-scatter
        //    values exactly.
        // ==============================================================
        for (std::size_t jj = 0u; jj < w_p; ++jj) {
            const Index j2 = jcol + static_cast<Index>(jj);
            const std::size_t sj2 = static_cast<std::size_t>(j2);
            T*     x  = &W[jj * un];
            Index* mk = &mark_p[jj * un];
            std::vector<Index>& pat = patterns[jj];
            for (Index p = A_csc.col_ptr[sj2]; p < A_csc.col_ptr[sj2 + 1u]; ++p) {
                const std::size_t sp       = static_cast<std::size_t>(p);
                const Index       orig_row = A_csc.row_ind[sp];
                const Index       piv_row  =
                    inv_row_perm[static_cast<std::size_t>(orig_row)];
                const std::size_t spiv     = static_cast<std::size_t>(piv_row);
                x[spiv] += A_csc.values[sp];
                if (mk[spiv] != j2) {
                    mk[spiv] = j2;
                    pat.push_back(piv_row);
                }
            }
        }

        // ==============================================================
        // 1. Panel DFS over the pruned supernode graph (union of the panel
        //    columns' A-pattern rows < jcol; traversal follows rows < jcol
        //    only -- those are pivotal and panel-stable).
        // ==============================================================
        reach_sups.clear();
        if (jcol > Index(0)) {
            const std::size_t nsup_now = sup_first.size();
            dfs_visited.resize(nsup_now, Index(-1));
            for (Index j = jcol; j < jend; ++j) {
                const std::size_t sj = static_cast<std::size_t>(j);
                for (Index p = A_csc.col_ptr[sj]; p < A_csc.col_ptr[sj + 1u]; ++p) {
                    const Index orig_row = A_csc.row_ind[static_cast<std::size_t>(p)];
                    const Index piv_row  =
                        inv_row_perm[static_cast<std::size_t>(orig_row)];
                    if (piv_row >= jcol) continue;
                    const Index s0 = col2sup[static_cast<std::size_t>(piv_row)];
                    if (dfs_visited[static_cast<std::size_t>(s0)] == jcol) continue;
                    dfs_stack.push_back(s0);
                    while (!dfs_stack.empty()) {
                        const Index s = dfs_stack.back();
                        dfs_stack.pop_back();
                        const std::size_t ss = static_cast<std::size_t>(s);
                        if (dfs_visited[ss] == jcol) continue;
                        dfs_visited[ss] = jcol;
                        reach_sups.push_back(s);
                        const std::vector<Index>& adj  = adj_slots[ss];
                        const std::vector<Index>& rows = sup_rows[ss];
                        const Index plen = prune_len[ss];
                        for (Index a = Index(0); a < plen; ++a) {
                            const Index r2 =
                                rows[static_cast<std::size_t>(
                                    adj[static_cast<std::size_t>(a)])];
                            if (r2 >= jcol) continue;
                            const Index s2 = col2sup[static_cast<std::size_t>(r2)];
                            if (dfs_visited[static_cast<std::size_t>(s2)] != jcol) {
                                dfs_stack.push_back(s2);
                            }
                        }
                    }
                }
            }
            std::sort(reach_sups.begin(), reach_sups.end());
            result.panel_dfs_visit_count += reach_sups.size();
        }

        // ==============================================================
        // 2. Panel update (design §2.2): per supernode K in ascending
        //    order -- gather zero-padded U fragments, scalar diagonal-block
        //    solve, ONE tgemm, scatter into the panel workspace.
        // ==============================================================
        for (std::size_t rs = 0u; rs < reach_sups.size(); ++rs) {
            const Index       s  = reach_sups[rs];
            const std::size_t ss = static_cast<std::size_t>(s);
            const Index fc = sup_first[ss];
            const Index nc_eff = std::min<Index>(sup_ncols[ss], jcol - fc);
            if (nc_eff <= Index(0)) continue;
            const std::vector<Index>& rows = sup_rows[ss];
            const std::vector<T>&     vals = sup_vals[ss];
            const Index ld = static_cast<Index>(rows.size());
            const Index m  = ld - (nc_eff - Index(1));   // ext row count
            const std::size_t snc = static_cast<std::size_t>(nc_eff);

            // Width-1 fast path: a single-column supernode needs no
            // fragment buffer, no triangular solve and no GEMM -- it is
            // exactly the GP column axpy.  This is the common case under
            // the fundamental partition (measured avg width ~1.3) and
            // avoids the fragment machinery overhead per (K, column).
            if (nc_eff == Index(1)) {
                const std::size_t sfc = static_cast<std::size_t>(fc);
                for (std::size_t jj = 0u; jj < w_p; ++jj) {
                    const Index j2 = jcol + static_cast<Index>(jj);
                    T*     x  = &W[jj * un];
                    Index* mk = &mark_p[jj * un];
                    if (mk[sfc] != j2) continue;
                    const T u_kj = x[sfc];
                    if (sparse_lu_scalar_policy<T>::is_exact_zero(u_kj))
                        continue;
                    std::vector<Index>& pat = patterns[jj];
                    result.axpy_element_count += static_cast<std::size_t>(ld);
                    for (Index t = Index(0); t < ld; ++t) {
                        const std::size_t st = static_cast<std::size_t>(t);
                        const Index       i  = rows[st];
                        const std::size_t si = static_cast<std::size_t>(i);
                        x[si] -= vals[st] * u_kj;
                        if (mk[si] != j2) {
                            mk[si] = j2;
                            pat.push_back(i);
                        }
                    }
                }
                continue;
            }

            // GEMM dispatch decided UP FRONT: below the gate the update is
            // applied as the direct member-column axpy (the exact GP order,
            // no fragment traffic -- measured to matter for kv::dd), above
            // it as the §2.2 gather / trisolve / GEMM / scatter pipeline.
            const std::size_t sm_pre = static_cast<std::size_t>(
                (m > Index(0)) ? m : Index(0));
            const bool use_gemm =
                (m > Index(0)) && (sm_pre * snc >= gemm_min_elems);
            if (!use_gemm) {
                ++result.scalar_update_count;
                for (std::size_t jj = 0u; jj < w_p; ++jj) {
                    const Index j2 = jcol + static_cast<Index>(jj);
                    T*     x  = &W[jj * un];
                    Index* mk = &mark_p[jj * un];
                    std::vector<Index>& pat = patterns[jj];
                    for (Index c = Index(0); c < nc_eff; ++c) {
                        const Index k = fc + c;
                        const std::size_t sk = static_cast<std::size_t>(k);
                        if (mk[sk] != j2) continue;
                        const T u_kj = x[sk];
                        const std::size_t base =
                            static_cast<std::size_t>(c) *
                            static_cast<std::size_t>(ld);
                        result.axpy_element_count +=
                            static_cast<std::size_t>(ld - c);
                        for (Index t = c; t < ld; ++t) {
                            const std::size_t st = static_cast<std::size_t>(t);
                            const Index       i  = rows[st];
                            const std::size_t si = static_cast<std::size_t>(i);
                            x[si] -= vals[base + st] * u_kj;
                            if (mk[si] != j2) {
                                mk[si] = j2;
                                pat.push_back(i);
                            }
                        }
                    }
                }
                continue;
            }

            // ----- 2a/2b. gather + scalar trisolve per active column -----
            frag_B.assign(snc * w_p, T(0));
            frag_active.assign(w_p, 0);
            bool any_active = false;
            const clock_type::time_point t_trsm0 = clock_type::now();
            for (std::size_t jj = 0u; jj < w_p; ++jj) {
                const Index j2 = jcol + static_cast<Index>(jj);
                T*     x  = &W[jj * un];
                Index* mk = &mark_p[jj * un];
                // activity: any MARKED row in K's column range (GP applies
                // marked columns regardless of value; unmarked rows are
                // exact zeros in W, so an all-unmarked fragment contributes
                // nothing).
                bool act = false;
                for (Index c = Index(0); c < nc_eff; ++c) {
                    if (mk[static_cast<std::size_t>(fc + c)] == j2) {
                        act = true;
                        break;
                    }
                }
                if (!act) continue;
                frag_active[jj] = 1;
                any_active = true;
                T* u = &frag_B[jj * snc];
                for (Index c = Index(0); c < nc_eff; ++c)
                    u[static_cast<std::size_t>(c)] =
                        x[static_cast<std::size_t>(fc + c)];
                // unit-lower solve with K's diagonal block:
                //   L(fc+c, fc+c') = vals[c'*ld + (c-1)]  (c > c')
                for (Index cp = Index(0); cp < nc_eff; ++cp) {
                    const T uc = u[static_cast<std::size_t>(cp)];
                    if (sparse_lu_scalar_policy<T>::is_exact_zero(uc)) continue;
                    const std::size_t base =
                        static_cast<std::size_t>(cp) * static_cast<std::size_t>(ld);
                    for (Index c = cp + Index(1); c < nc_eff; ++c) {
                        u[static_cast<std::size_t>(c)] -=
                            vals[base + static_cast<std::size_t>(c - Index(1))] * uc;
                    }
                }
                // write back the solved U segment + mark (rows are < jcol,
                // stable; stored later as U entries of column j2)
                for (Index c = Index(0); c < nc_eff; ++c) {
                    const Index r = fc + c;
                    const std::size_t sr = static_cast<std::size_t>(r);
                    x[sr] = u[static_cast<std::size_t>(c)];
                    if (mk[sr] != j2) {
                        mk[sr] = j2;
                        patterns[jj].push_back(r);
                    }
                }
            }
            result.trsm_time_ns += static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    clock_type::now() - t_trsm0).count());
            if (!any_active) continue;

            // ----- 2c. one GEMM per (supernode, panel) -----
            const std::size_t sm = static_cast<std::size_t>(m);
            {
                gemm_C.resize(sm * w_p);
                const clock_type::time_point t_g0 = clock_type::now();
                tgemm<T>('N', 'N',
                         static_cast<int>(m), static_cast<int>(w_p),
                         static_cast<int>(nc_eff),
                         T(1),
                         vals.data() + static_cast<std::size_t>(nc_eff - Index(1)),
                         static_cast<int>(ld),
                         frag_B.data(), static_cast<int>(nc_eff),
                         T(0), gemm_C.data(), static_cast<int>(m));
                result.gemm_time_ns += static_cast<std::size_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        clock_type::now() - t_g0).count());
                ++result.gemm_call_count;
                result.gemm_flop_count +=
                    2u * sm * w_p * snc;
                // ----- 2d. scatter into the panel workspace -----
                for (std::size_t jj = 0u; jj < w_p; ++jj) {
                    if (!frag_active[jj]) continue;
                    const Index j2 = jcol + static_cast<Index>(jj);
                    T*     x  = &W[jj * un];
                    Index* mk = &mark_p[jj * un];
                    std::vector<Index>& pat = patterns[jj];
                    const T* Cj = &gemm_C[jj * sm];
                    for (Index t = Index(0); t < m; ++t) {
                        const std::size_t st = static_cast<std::size_t>(t);
                        const Index r =
                            rows[static_cast<std::size_t>(nc_eff - Index(1) + t)];
                        const std::size_t sr = static_cast<std::size_t>(r);
                        x[sr] -= Cj[st];
                        if (mk[sr] != j2) {
                            mk[sr] = j2;
                            pat.push_back(r);
                        }
                    }
                }
            }
        }

        // ==============================================================
        // 3. Within-panel factorization (scalar, GP order)
        // ==============================================================
        for (Index j = jcol; j < jend; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const std::size_t jl = static_cast<std::size_t>(j - jcol);
            T*     x  = &W[jl * un];
            Index* mk = &mark_p[jl * un];
            std::vector<Index>& pattern = patterns[jl];

            // ----- 3a. in-panel updates (ascending scan, GP order) -----
            for (Index k = jcol; k < j; ++k) {
                const std::size_t sk = static_cast<std::size_t>(k);
                if (mk[sk] != j) continue;
                const T u_kj = x[sk];
                const Index       s  = col2sup[sk];
                const std::size_t ss = static_cast<std::size_t>(s);
                const Index       c  = k - sup_first[ss];
                const std::vector<Index>& rows = sup_rows[ss];
                const std::vector<T>&     vals = sup_vals[ss];
                const Index ld = static_cast<Index>(rows.size());
                const std::size_t base =
                    static_cast<std::size_t>(c) * static_cast<std::size_t>(ld);
                result.axpy_element_count += static_cast<std::size_t>(ld - c);
                for (Index t = c; t < ld; ++t) {
                    const std::size_t st = static_cast<std::size_t>(t);
                    const Index       i  = rows[st];
                    const std::size_t si = static_cast<std::size_t>(i);
                    x[si] -= vals[base + st] * u_kj;
                    if (mk[si] != j) {
                        mk[si] = j;
                        pattern.push_back(i);
                    }
                }
            }

            // ----- 3b. Pivot selection (GP step 3; stale-entry guard) -----
            Index     pivot_pos = j;
            real_type col_max   = real_type(0);
            T         pivot_val = T(0);

            if (do_partial_pivot) {
                for (std::size_t pi = 0u; pi < pattern.size(); ++pi) {
                    const Index i = pattern[pi];
                    if (i < j) continue;
                    const std::size_t si = static_cast<std::size_t>(i);
                    if (mk[si] != j) continue;   // stale after a row swap
                    const real_type ai =
                        sparse_lu_scalar_policy<T>::abs_value(x[si]);
                    if (ai > col_max) {
                        col_max   = ai;
                        pivot_pos = i;
                        pivot_val = x[si];
                    }
                }
            } else {
                pivot_val = x[sj];
                col_max   = sparse_lu_scalar_policy<T>::abs_value(x[sj]);
            }

            if (!sparse_lu_is_acceptable_pivot(pivot_val, pivot_val,
                                               threshold, abs_tol)) {
                result.success = false;
                result.status  = sparse_lu_status::numerical_singularity;
                return result;
            }

            // ----- 3c. Row swap (GP step 5) -----
            if (pivot_pos != j) {
                const std::size_t sp = static_cast<std::size_t>(pivot_pos);
                // current column: swap values only (marks handled below,
                // mirroring the GP semantics)
                std::swap(x[sj], x[sp]);
                // later panel columns: swap workspace values AND marks;
                // push newly-marked rows (duplicates are deduped at store)
                for (std::size_t jj = jl + 1u; jj < w_p; ++jj) {
                    const Index j2 = jcol + static_cast<Index>(jj);
                    T*     x2  = &W[jj * un];
                    Index* mk2 = &mark_p[jj * un];
                    std::swap(x2[sj], x2[sp]);
                    std::swap(mk2[sj], mk2[sp]);
                    if (mk2[sj] == j2) patterns[jj].push_back(j);
                    if (mk2[sp] == j2) patterns[jj].push_back(pivot_pos);
                }
                if (mk[sj] != j) {
                    mk[sj] = j;
                    pattern.push_back(j);
                }
                const Index r_j = row_perm[sj];
                const Index r_p = row_perm[sp];
                std::swap(row_perm[sj], row_perm[sp]);
                inv_row_perm[static_cast<std::size_t>(r_j)] = pivot_pos;
                inv_row_perm[static_cast<std::size_t>(r_p)] = j;

                // rename/swap rows j <-> pivot_pos in the stored supernode
                // trapezoids via the reverse index (three-case merge;
                // mirrors SLU-NQ2 -- see the Phase 2 commit message)
                std::vector<row_loc>& lj = row_locs[sj];
                std::vector<row_loc>& lp = row_locs[sp];
                std::vector<row_loc> nj, np;
                nj.reserve(lj.size() + lp.size());
                np.reserve(lj.size() + lp.size());
                std::size_t a = 0u, b = 0u;
                while (a < lj.size() || b < lp.size()) {
                    if (b >= lp.size() ||
                        (a < lj.size() && lj[a].sup < lp[b].sup)) {
                        const std::size_t s2 =
                            static_cast<std::size_t>(lj[a].sup);
                        sup_rows[s2][static_cast<std::size_t>(lj[a].slot)] =
                            pivot_pos;
                        np.push_back(lj[a]); ++a;
                    } else if (a >= lj.size() || lp[b].sup < lj[a].sup) {
                        const std::size_t s2 =
                            static_cast<std::size_t>(lp[b].sup);
                        sup_rows[s2][static_cast<std::size_t>(lp[b].slot)] = j;
                        nj.push_back(lp[b]); ++b;
                    } else {
                        const std::size_t s2 =
                            static_cast<std::size_t>(lj[a].sup);
                        std::vector<T>& vals2 = sup_vals[s2];
                        const std::size_t ld2 = sup_rows[s2].size();
                        const std::size_t ncst =
                            static_cast<std::size_t>(sup_ncols[s2]);
                        const std::size_t t1 =
                            static_cast<std::size_t>(lj[a].slot);
                        const std::size_t t2 =
                            static_cast<std::size_t>(lp[b].slot);
                        for (std::size_t c = 0u; c < ncst; ++c) {
                            std::swap(vals2[c * ld2 + t1], vals2[c * ld2 + t2]);
                        }
                        nj.push_back(lj[a]); ++a;
                        np.push_back(lp[b]); ++b;
                    }
                }
                lj.swap(nj);
                lp.swap(np);
            } else {
                if (mk[sj] != j) {
                    mk[sj] = j;
                    pattern.push_back(j);
                }
            }

            const T pivot_diag = x[sj];

            // ----- 3d. Store U (CSC, GP rules) + supernode decision + L -----
            std::sort(pattern.begin(), pattern.end());

            // ext_count counts the rows > j that will be STORED (exact
            // zeros dropped -- the GP store rule; see the Phase 2 commit
            // message for the 14x inflation this prevents).  Stale and
            // duplicate pattern entries are skipped here.
            Index u_cnt = Index(0);
            std::size_t ext_count = 0u;
            for (std::size_t pi = 0u; pi < pattern.size(); ++pi) {
                const Index i = pattern[pi];
                if (pi > 0u && i == pattern[pi - 1u]) continue;   // duplicate
                const std::size_t si = static_cast<std::size_t>(i);
                if (i < j) {
                    if (mk[si] == j &&
                        !sparse_lu_scalar_policy<T>::is_exact_zero(x[si])) {
                        U_row_ind_buf.push_back(i);
                        U_val_buf.push_back(x[si]);
                        ++u_cnt;
                    }
                } else if (i == j) {
                    U_row_ind_buf.push_back(i);
                    U_val_buf.push_back(x[si]);
                    ++u_cnt;
                } else if (mk[si] == j &&
                           !sparse_lu_scalar_policy<T>::is_exact_zero(x[si])) {
                    ++ext_count;
                }
            }
            U_col_ptr[sj + 1u] = U_col_ptr[sj] + u_cnt;

            // Dynamic supernode decision (see the file header): column j
            // joins the supernode of column j-1 when the width cap is not
            // reached, slot c-1 holds row j (diagonal-run condition), and
            // the STORED sub-diagonal pattern of j is contained in the
            // remaining slots with at most relax_allow inactive slots
            // (weak relaxation; relax_allow == 0 = fundamental partition).
            bool joined = false;
            if (j > Index(0) && col2sup[sj - 1u] != Index(-1)) {
                const Index       s  = col2sup[sj - 1u];
                const std::size_t ss = static_cast<std::size_t>(s);
                const Index c = j - sup_first[ss];
                const std::vector<Index>& rows = sup_rows[ss];
                const Index ld = static_cast<Index>(rows.size());
                if (static_cast<std::size_t>(c) < maxsup &&
                    c >= Index(1) && c - Index(1) < ld &&
                    rows[static_cast<std::size_t>(c - 1)] == j &&
                    c <= ld) {
                    std::size_t active_slots = 0u;
                    for (Index t = c; t < ld; ++t) {
                        const std::size_t si = static_cast<std::size_t>(
                            rows[static_cast<std::size_t>(t)]);
                        if (mk[si] == j &&
                            !sparse_lu_scalar_policy<T>::is_exact_zero(x[si]))
                            ++active_slots;
                    }
                    const std::size_t slot_count =
                        static_cast<std::size_t>(ld - c);
                    if (active_slots == ext_count &&
                        slot_count - active_slots <= relax_allow) {
                        std::vector<T>& vals = sup_vals[ss];
                        vals.resize(static_cast<std::size_t>(c + 1) *
                                        static_cast<std::size_t>(ld),
                                    T(0));
                        const std::size_t base =
                            static_cast<std::size_t>(c) *
                            static_cast<std::size_t>(ld);
                        for (Index t = c; t < ld; ++t) {
                            const std::size_t st = static_cast<std::size_t>(t);
                            const std::size_t si = static_cast<std::size_t>(
                                rows[st]);
                            vals[base + st] =
                                (mk[si] == j &&
                                 !sparse_lu_scalar_policy<T>::is_exact_zero(x[si]))
                                    ? x[si] / pivot_diag
                                    : T(0);
                        }
                        sup_ncols[ss] = c + Index(1);
                        col2sup[sj]   = s;
                        result.slot_element_count +=
                            static_cast<std::size_t>(ld - c);
                        joined = true;
                    }
                }
            }
            if (!joined) {
                // create a new supernode from the STORED sub-diagonal pattern
                const Index s = static_cast<Index>(sup_first.size());
                sup_first.push_back(j);
                sup_ncols.push_back(Index(1));
                sup_rows.push_back(std::vector<Index>());
                sup_vals.push_back(std::vector<T>());
                std::vector<Index>& rows = sup_rows.back();
                std::vector<T>&     vals = sup_vals.back();
                rows.reserve(ext_count);
                vals.reserve(ext_count);
                for (std::size_t pi = 0u; pi < pattern.size(); ++pi) {
                    const Index i = pattern[pi];
                    if (i <= j) continue;
                    if (pi > 0u && i == pattern[pi - 1u]) continue;
                    const std::size_t si = static_cast<std::size_t>(i);
                    if (mk[si] != j) continue;   // stale after a row swap
                    if (sparse_lu_scalar_policy<T>::is_exact_zero(x[si]))
                        continue;                // GP store rule
                    const std::size_t slot = rows.size();
                    rows.push_back(i);
                    vals.push_back(x[si] / pivot_diag);
                    if (do_partial_pivot) {
                        row_locs[si].push_back(
                            row_loc{ s, static_cast<Index>(slot) });
                    }
                }
                adj_slots.push_back(std::vector<Index>());
                std::vector<Index>& adj = adj_slots.back();
                adj.resize(rows.size());
                for (std::size_t t = 0u; t < rows.size(); ++t)
                    adj[t] = static_cast<Index>(t);
                prune_len.push_back(static_cast<Index>(rows.size()));
                col2sup[sj] = s;
                result.slot_element_count += rows.size();
            }

            // ----- 3e. Eisenstat-Liu symmetric pruning -----
            // Trigger: supernode s has a stored U entry in column j AND row
            // j appears in s's row list.  Intra-supernode triggers are
            // skipped (self-loop degeneracy).  The partition KEEPS rows
            // <= j (pivotal, stable) and any row > j whose fill in column
            // j was not stored (replacement-path verification; exact
            // cancellation keeps the edge -- the safe side).
            {
                const std::size_t nsup_now = sup_first.size();
                u_sup_stamp.resize(nsup_now, Index(-1));
                for (Index up = U_col_ptr[sj]; up < U_col_ptr[sj + 1u]; ++up) {
                    const Index k = U_row_ind_buf[static_cast<std::size_t>(up)];
                    if (k >= j) continue;
                    const Index s = col2sup[static_cast<std::size_t>(k)];
                    u_sup_stamp[static_cast<std::size_t>(s)] = j;
                }
                if (do_partial_pivot) {
                    const std::vector<row_loc>& locs = row_locs[sj];
                    for (std::size_t li = 0u; li < locs.size(); ++li) {
                        const Index       s  = locs[li].sup;
                        const std::size_t ss = static_cast<std::size_t>(s);
                        if (u_sup_stamp[ss] != j) continue;
                        if (col2sup[sj] == s) continue;   // self-supernode guard
                        std::vector<Index>& adj = adj_slots[ss];
                        const std::vector<Index>& rows = sup_rows[ss];
                        Index lo = Index(0);
                        Index hi = prune_len[ss];
                        while (lo < hi) {
                            const Index slot = adj[static_cast<std::size_t>(lo)];
                            const Index r = rows[static_cast<std::size_t>(slot)];
                            const std::size_t sr = static_cast<std::size_t>(r);
                            const bool keep =
                                (r <= j) || (mk[sr] != j) ||
                                sparse_lu_scalar_policy<T>::is_exact_zero(x[sr]);
                            if (keep) {
                                ++lo;
                            } else {
                                --hi;
                                std::swap(adj[static_cast<std::size_t>(lo)],
                                          adj[static_cast<std::size_t>(hi)]);
                            }
                        }
                        if (lo < prune_len[ss]) {
                            prune_len[ss] = lo;
                            ++result.prune_event_count;
                        }
                    }
                } else {
                    for (Index up = U_col_ptr[sj]; up < U_col_ptr[sj + 1u]; ++up) {
                        const Index k =
                            U_row_ind_buf[static_cast<std::size_t>(up)];
                        if (k >= j) continue;
                        const Index       s  = col2sup[static_cast<std::size_t>(k)];
                        const std::size_t ss = static_cast<std::size_t>(s);
                        if (u_sup_stamp[ss] != j) continue;
                        u_sup_stamp[ss] = Index(-1);   // process once
                        if (col2sup[sj] == s) continue;   // self-supernode guard
                        std::vector<Index>& adj = adj_slots[ss];
                        const std::vector<Index>& rows = sup_rows[ss];
                        bool has_j = false;
                        for (Index a = Index(0); a < prune_len[ss]; ++a) {
                            if (rows[static_cast<std::size_t>(
                                    adj[static_cast<std::size_t>(a)])] == j) {
                                has_j = true;
                                break;
                            }
                        }
                        if (!has_j) continue;
                        Index lo = Index(0);
                        Index hi = prune_len[ss];
                        while (lo < hi) {
                            const Index slot = adj[static_cast<std::size_t>(lo)];
                            const Index r = rows[static_cast<std::size_t>(slot)];
                            const std::size_t sr = static_cast<std::size_t>(r);
                            const bool keep =
                                (r <= j) || (mk[sr] != j) ||
                                sparse_lu_scalar_policy<T>::is_exact_zero(x[sr]);
                            if (keep) {
                                ++lo;
                            } else {
                                --hi;
                                std::swap(adj[static_cast<std::size_t>(lo)],
                                          adj[static_cast<std::size_t>(hi)]);
                            }
                        }
                        if (lo < prune_len[ss]) {
                            prune_len[ss] = lo;
                            ++result.prune_event_count;
                        }
                    }
                }
            }

            // ----- 3f. Clear this column of the panel workspace -----
            for (std::size_t pi = 0u; pi < pattern.size(); ++pi)
                x[static_cast<std::size_t>(pattern[pi])] = T(0);
            x[sj] = T(0);
            pattern.clear();
        }
    }

    // ------------------------------------------------------------------
    // Emit CSC L from the supernode trapezoids (exact zeros skipped -- the
    // GP store rule).  Entry order within a column is the slot order
    // (ascending at supernode creation; later renames may perturb it,
    // which the CSC solve tolerates).
    // ------------------------------------------------------------------
    baseline_lu_storage<T, Index>& lu = result.storage;
    lu.L.col_ptr.resize(un + 1u, Index(0));

    const std::size_t nsup = sup_first.size();
    {
        std::size_t nnz_L = 0u;
        for (std::size_t ss = 0u; ss < nsup; ++ss) {
            const std::vector<T>& vals = sup_vals[ss];
            const std::size_t ld = sup_rows[ss].size();
            const std::size_t nc = static_cast<std::size_t>(sup_ncols[ss]);
            for (std::size_t c = 0u; c < nc; ++c) {
                const Index j = sup_first[ss] + static_cast<Index>(c);
                std::size_t cnt = 0u;
                for (std::size_t t = c; t < ld; ++t) {
                    if (!sparse_lu_scalar_policy<T>::is_exact_zero(
                            vals[c * ld + t])) ++cnt;
                }
                lu.L.col_ptr[static_cast<std::size_t>(j) + 1u] =
                    static_cast<Index>(cnt);
                nnz_L += cnt;
            }
        }
        for (std::size_t jj = 0u; jj < un; ++jj)
            lu.L.col_ptr[jj + 1u] =
                lu.L.col_ptr[jj + 1u] + lu.L.col_ptr[jj];
        lu.L.row_ind.resize(nnz_L);
        lu.L.values.resize(nnz_L);
        for (std::size_t ss = 0u; ss < nsup; ++ss) {
            const std::vector<Index>& rows = sup_rows[ss];
            const std::vector<T>&     vals = sup_vals[ss];
            const std::size_t ld = rows.size();
            const std::size_t nc = static_cast<std::size_t>(sup_ncols[ss]);
            for (std::size_t c = 0u; c < nc; ++c) {
                const Index j = sup_first[ss] + static_cast<Index>(c);
                std::size_t out = static_cast<std::size_t>(
                    lu.L.col_ptr[static_cast<std::size_t>(j)]);
                for (std::size_t t = c; t < ld; ++t) {
                    const T& v = vals[c * ld + t];
                    if (sparse_lu_scalar_policy<T>::is_exact_zero(v)) continue;
                    lu.L.row_ind[out] = rows[t];
                    lu.L.values[out]  = v;
                    ++out;
                }
            }
        }
    }

    lu.U.col_ptr = U_col_ptr;
    lu.U.row_ind = U_row_ind_buf;
    lu.U.values  = U_val_buf;

    lu.row_perm     = row_perm;
    lu.inv_row_perm = sparse_lu_inverse_permutation(row_perm);
    lu.col_perm     = col_perm;
    lu.inv_col_perm = inv_col_perm;

    result.number_of_supernodes = static_cast<Index>(nsup);
    result.success = true;
    result.status  = sparse_lu_status::success;
    return result;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_SUPERNODE_PANEL_IMPL_HPP
