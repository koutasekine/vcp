// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// ORD-1 -- multilevel nested dissection ordering (dependency-free, METIS-class
// target; design ORD-1_design_v1 SS2, ruling ORD-D5 recorded in
// sandbox/docs/issues/ORD-1_stop1_issue.md SS0).
//
// This file MUST be #included from WITHIN namespace vcp, AFTER the ordering
// implementation header (tsparse_sparse_lu_ordering_impl.hpp) is in scope: it
// reuses sparse_lu_build_symmetric_pattern_graph (A + A^T pattern graph) and
// sparse_lu_amd_ordering (leaf ordering).  It has no "namespace vcp { }"
// wrapper; it is injected by tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse.hpp>                    (umbrella)
//
// Algorithm (design SS2, the CHOLMOD NESDIS / METIS standard construction):
// recursive vertex dissection.  Each piece larger than the leaf size is
// bisected by the multilevel scheme
//     coarsen (heavy-edge matching) -> initial bisection (deterministic
//     greedy graph growing on the coarsest graph) -> uncoarsen with FM
//     (Fiduccia-Mattheyses) boundary refinement per level
// and the edge cut is turned into a vertex separator by a greedy
// minimum-vertex-cover approximation.  The two halves are ordered first
// (recursively), the separator LAST (near the elimination-tree root, where
// it forms the wide dense fronts).  Pieces at or below the leaf size are
// ordered by the existing AMD on the induced subgraph (constrained-AMD
// simple form, design SS2-5).
//
// Safety invariants (same contract as the sibling orderings):
//   S-1 Pattern-only / integer-only.  No numeric values are read.
//   S-5 Deterministic (ruling D-2).  No randomness anywhere; every choice
//       has an explicit total order:
//         - matching: vertices visited ascending; the mate is the unmatched
//           neighbour of maximum edge weight, tie-broken by smallest index;
//         - initial bisection: start vertices enumerated ascending (strided
//           when the coarsest graph is larger than 128 vertices); the grown
//           vertex is the maximum-gain frontier vertex, tie smallest index;
//           the winning start is the lexicographically smallest
//           (cut, imbalance, start-order) triple, strict improvement only;
//         - FM: maximum-gain feasible move, tie smallest index; a pass is
//           rolled back to its lexicographically best
//           (balance-violation, cut) prefix, strict improvement only;
//         - vertex cover: maximum cut-degree vertex, tie preferring the
//           side from which fewer separator vertices have been taken so
//           far, then smallest index;
//         - every emitted leaf / separator block is ascending or AMD-ordered
//           (AMD itself is deterministic, S-5 of the ordering track).
//       Same input (and same parameters) therefore yields a byte-identical
//       permutation.  Recursion uses an explicit work stack (no call-stack
//       recursion; deep dissection trees cannot overflow the stack).
//
// Convention: returns a permutation with perm[new] = old, exactly like
// sparse_lu_rcm_ordering / sparse_lu_amd_ordering /
// sparse_lu_nested_dissection_ordering (verified against the AMD
// implementation, ordering impl line "perm[new] = old").
//
// Integer-only parameters (P4 discipline of the surrounding module): the
// balance bound beta is carried as an integer percentage (balance_pct = 40
// means beta = 0.40), so the header stays free of floating arithmetic.

#ifndef VCP_TSPARSE_ORDER_NDML_IMPL_HPP
#define VCP_TSPARSE_ORDER_NDML_IMPL_HPP

// NOTE: this file is textually injected INSIDE namespace vcp; every standard
// header it needs must already be included by the injecting header BEFORE the
// namespace opens.  The includes below are no-op guards when that discipline
// is followed (tsparse_sparse_lu.hpp includes <set> for SLU-OQ1 and the rest
// through the ordering implementation header).
#include <algorithm>
#include <cstddef>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/error.hpp>

// ---------------------------------------------------------------------------
// Parameters (design SS2-6).  Defaults are the design initial values; the
// final defaults are fixed by the OR-2 calibration (one line, four constants).
// The struct is Index-independent so the three factorizations can embed the
// same fields in their options without templating on the orderer.
// ---------------------------------------------------------------------------
struct sparse_order_ndml_params {
    long long coarsen_stop;   // stop coarsening at this many vertices
    int       fm_passes;      // FM refinement passes per uncoarsening level
    int       balance_pct;    // min side weight, percent of total (beta*100)
    long long leaf_size;      // pieces at or below this size are AMD-ordered
    // Separator acceptance bound (the NESDIS nd_oksep guard, internal --
    // deliberately NOT exposed in the factorization options, whose surface
    // is the four design SS2-6 parameters above): a piece whose best
    // refined separator exceeds max_sep_pct percent of the piece weight has
    // no useful separator (expander-shaped graphs); it is not split and is
    // ordered by the halo-augmented AMD instead.  This is what keeps the
    // no-degradation acceptance on separator-free graphs (G4).
    int       max_sep_pct;
    // Best-of-two selection (internal, like max_sep_pct): after building the
    // ND order, also build the plain AMD order and return whichever gives the
    // smaller exact symbolic fill (ties keep ND).  This is the standard
    // construction of production analyzers (CHOLMOD analyze tries its
    // ordering candidates and keeps the best): on graphs whose separators
    // are structurally unprofitable -- e.g. banded matrices, where any leaf
    // column drags its flanking separator and pays ~2x the band -- pure
    // vertex-separator ND cannot beat AMD, and the deterministic comparison
    // is what keeps the no-degradation acceptance (G5).  0 disables the
    // comparison (used by structural unit tests that pin the raw ND shape).
    int       best_of_two;

    // OR-2 calibration (3D shifted-FEM m=20, nnz_L only, 12-point grid
    // {leaf 100/200/400} x {fm 1/2} x {beta 40/45}, coarsen_stop 128):
    // winner leaf=100 / fm=2 / beta=40 (nnz_L 1,369,133 = 1.0250 x golden;
    // confirmed at m=30: 1.0438 vs 1.0592 for leaf=200).
    sparse_order_ndml_params()
        : coarsen_stop(128), fm_passes(2), balance_pct(40), leaf_size(100),
          max_sep_pct(30), best_of_two(1) {}
};

namespace sparse_ndml_detail {

// ---------------------------------------------------------------------------
// Weighted level graph.  xadj/adjncy is the CSR adjacency (neighbour lists
// ascending, no self-loops, no duplicates), ewt the matching edge weights
// (collapsed multi-edge mass), vwt the vertex weights (collapsed vertex
// mass).  The finest level of a piece has all weights 1.
// ---------------------------------------------------------------------------
template <class Index>
struct ndml_graph {
    Index nv;
    std::vector<Index> xadj, adjncy, ewt, vwt;
    std::vector<Index> cmap;   // this level's vertex -> next-coarser vertex
                               // (filled by the coarsening step)
};

// ---------------------------------------------------------------------------
// Induced subgraph of the global pattern graph on `verts` (ascending).
// Local index = rank in `verts`; since `verts` ascends and the global
// neighbour lists ascend, the local lists come out ascending without a sort.
// `loc` is a size-n scratch (-1 outside; reset before returning).
// ---------------------------------------------------------------------------
template <class Index>
void ndml_extract_subgraph(
    const std::vector<Index>& verts,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    std::vector<Index>& loc,
    ndml_graph<Index>& g)
{
    const std::size_t ns = verts.size();
    g.nv = static_cast<Index>(ns);
    for (std::size_t i = 0; i < ns; ++i) {
        loc[static_cast<std::size_t>(verts[i])] = static_cast<Index>(i);
    }
    g.xadj.assign(ns + 1u, Index(0));
    g.adjncy.clear();
    for (std::size_t i = 0; i < ns; ++i) {
        const std::size_t b = static_cast<std::size_t>(
            adj_ptr[static_cast<std::size_t>(verts[i])]);
        const std::size_t e = static_cast<std::size_t>(
            adj_ptr[static_cast<std::size_t>(verts[i]) + 1u]);
        for (std::size_t p = b; p < e; ++p) {
            const Index w = loc[static_cast<std::size_t>(adj_ind[p])];
            if (w >= Index(0)) g.adjncy.push_back(w);
        }
        g.xadj[i + 1u] = static_cast<Index>(g.adjncy.size());
    }
    g.ewt.assign(g.adjncy.size(), Index(1));
    g.vwt.assign(ns, Index(1));
    g.cmap.clear();
    for (std::size_t i = 0; i < ns; ++i) {
        loc[static_cast<std::size_t>(verts[i])] = Index(-1);
    }
}

// ---------------------------------------------------------------------------
// One coarsening step: deterministic heavy-edge matching + contraction.
// Returns false (and leaves `coarse` untouched) when the matching stalls
// (fewer than 5% of the vertices got a mate), which ends the coarsening
// loop -- the initial bisection then runs on the current graph.
// ---------------------------------------------------------------------------
template <class Index>
bool ndml_coarsen_step(
    ndml_graph<Index>& fine,
    ndml_graph<Index>& coarse)
{
    const std::size_t nv = static_cast<std::size_t>(fine.nv);
    std::vector<Index> match(nv, Index(-1));
    std::size_t matched2 = 0;   // vertices matched to a DIFFERENT vertex
    for (std::size_t u = 0; u < nv; ++u) {
        if (match[u] != Index(-1)) continue;
        Index best = Index(-1);
        Index best_w = Index(0);
        const std::size_t b = static_cast<std::size_t>(fine.xadj[u]);
        const std::size_t e = static_cast<std::size_t>(fine.xadj[u + 1u]);
        for (std::size_t p = b; p < e; ++p) {
            const Index v = fine.adjncy[p];
            if (match[static_cast<std::size_t>(v)] != Index(-1)) continue;
            const Index w = fine.ewt[p];
            // heavier edge wins; ties by smallest neighbour index -- the
            // lists ascend, so a strict > keeps the first (smallest) index.
            if (best < Index(0) || w > best_w) { best = v; best_w = w; }
        }
        if (best >= Index(0)) {
            match[u] = best;
            match[static_cast<std::size_t>(best)] = static_cast<Index>(u);
            matched2 += 2u;
        } else {
            match[u] = static_cast<Index>(u);   // stays single
        }
    }
    if (matched2 < nv / 20u + 2u) return false;   // stalled (< ~5% matched)

    // Coarse ids in ascending order of the smaller endpoint (deterministic).
    fine.cmap.assign(nv, Index(-1));
    Index nc = Index(0);
    for (std::size_t u = 0; u < nv; ++u) {
        if (fine.cmap[u] != Index(-1)) continue;
        const std::size_t v = static_cast<std::size_t>(match[u]);
        fine.cmap[u] = nc;
        fine.cmap[v] = nc;   // v == u when single
        ++nc;
    }

    // Contract.  Neighbour weights are accumulated with a stamp buffer; the
    // collected coarse neighbours are sorted ascending per coarse vertex.
    const std::size_t unc = static_cast<std::size_t>(nc);
    coarse.nv = nc;
    coarse.xadj.assign(unc + 1u, Index(0));
    coarse.adjncy.clear();
    coarse.ewt.clear();
    coarse.vwt.assign(unc, Index(0));
    coarse.cmap.clear();
    std::vector<Index> wbuf(unc, Index(0));
    std::vector<Index> stamp(unc, Index(-1));
    std::vector<Index> nbrs;
    // members of coarse vertex c, in ascending fine order: u with cmap==c.
    // Iterating fine vertices ascending and bucketing gives them directly.
    std::vector<Index> head(unc, Index(-1)), nxt(nv, Index(-1)), tail(unc, Index(-1));
    for (std::size_t u = 0; u < nv; ++u) {
        const std::size_t c = static_cast<std::size_t>(fine.cmap[u]);
        if (head[c] < Index(0)) head[c] = static_cast<Index>(u);
        else nxt[static_cast<std::size_t>(tail[c])] = static_cast<Index>(u);
        tail[c] = static_cast<Index>(u);
        coarse.vwt[c] += fine.vwt[u];
    }
    for (std::size_t c = 0; c < unc; ++c) {
        nbrs.clear();
        for (Index iu = head[c]; iu >= Index(0); iu = nxt[static_cast<std::size_t>(iu)]) {
            const std::size_t u = static_cast<std::size_t>(iu);
            const std::size_t b = static_cast<std::size_t>(fine.xadj[u]);
            const std::size_t e = static_cast<std::size_t>(fine.xadj[u + 1u]);
            for (std::size_t p = b; p < e; ++p) {
                const std::size_t d = static_cast<std::size_t>(
                    fine.cmap[static_cast<std::size_t>(fine.adjncy[p])]);
                if (d == c) continue;                   // collapsed self-loop
                if (stamp[d] != static_cast<Index>(c)) {
                    stamp[d] = static_cast<Index>(c);
                    wbuf[d] = Index(0);
                    nbrs.push_back(static_cast<Index>(d));
                }
                wbuf[d] += fine.ewt[p];
            }
        }
        std::sort(nbrs.begin(), nbrs.end());
        for (std::size_t t = 0; t < nbrs.size(); ++t) {
            coarse.adjncy.push_back(nbrs[t]);
            coarse.ewt.push_back(wbuf[static_cast<std::size_t>(nbrs[t])]);
        }
        coarse.xadj[c + 1u] = static_cast<Index>(coarse.adjncy.size());
    }
    return true;
}

// ---------------------------------------------------------------------------
// Edge cut of a partition (sum of ewt over edges with differing sides; each
// undirected edge appears twice in CSR, so the sum is halved).
// ---------------------------------------------------------------------------
template <class Index>
Index ndml_cut_of(
    const ndml_graph<Index>& g,
    const std::vector<char>& side)
{
    Index cut2 = Index(0);
    const std::size_t nv = static_cast<std::size_t>(g.nv);
    for (std::size_t u = 0; u < nv; ++u) {
        const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
        const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
        for (std::size_t p = b; p < e; ++p) {
            if (side[u] != side[static_cast<std::size_t>(g.adjncy[p])]) {
                cut2 += g.ewt[p];
            }
        }
    }
    return cut2 / Index(2);
}

// ---------------------------------------------------------------------------
// Initial bisection on the coarsest graph: deterministic greedy graph
// growing (GGGP).  For every start vertex of a deterministic start set the
// region A is grown -- always absorbing the maximum-gain frontier vertex
// (gain = 2*w(v->A) - wdeg(v); tie smallest index) and jumping to the
// smallest-index unassigned vertex when the frontier empties (disconnected
// graphs) -- until weight(A)*2 >= total weight.  The winner is the strictly
// lexicographically smallest (cut, |W - 2*wA|, start-order) triple.
// ---------------------------------------------------------------------------
template <class Index>
void ndml_initial_bisection(
    const ndml_graph<Index>& g,
    std::vector<char>& side)
{
    const std::size_t nv = static_cast<std::size_t>(g.nv);
    Index W = Index(0);
    for (std::size_t i = 0; i < nv; ++i) W += g.vwt[i];
    std::vector<Index> wdeg(nv, Index(0));
    for (std::size_t u = 0; u < nv; ++u) {
        const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
        const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
        for (std::size_t p = b; p < e; ++p) wdeg[u] += g.ewt[p];
    }

    // Deterministic start set: every vertex when nv <= 128, else an
    // ascending stride keeping at most 128 starts.
    const std::size_t stride = (nv <= 128u) ? 1u : (nv + 127u) / 128u;

    std::vector<char> cur(nv);
    std::vector<Index> conn(nv);
    bool have_best = false;
    Index best_cut = Index(0), best_imb = Index(0);
    for (std::size_t s0 = 0; s0 < nv; s0 += stride) {
        std::fill(cur.begin(), cur.end(), char(1));
        std::fill(conn.begin(), conn.end(), Index(0));
        // cand keyed (-gain, v): begin() = max gain, tie smallest index.
        std::set<std::pair<Index, Index> > cand;
        std::size_t cursor = 0;      // min-index scan position for reseeding
        Index wA = Index(0);
        Index grown = static_cast<Index>(s0);
        cur[s0] = 0;
        wA += g.vwt[s0];
        while (wA * Index(2) < W) {
            // relax the freshly grown vertex' neighbours
            {
                const std::size_t u = static_cast<std::size_t>(grown);
                const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
                const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
                for (std::size_t p = b; p < e; ++p) {
                    const Index v = g.adjncy[p];
                    const std::size_t uv = static_cast<std::size_t>(v);
                    if (cur[uv] == 0) continue;
                    if (conn[uv] > Index(0)) {
                        cand.erase(std::make_pair(
                            wdeg[uv] - Index(2) * conn[uv], v));
                    }
                    conn[uv] += g.ewt[p];
                    cand.insert(std::make_pair(
                        wdeg[uv] - Index(2) * conn[uv], v));
                }
            }
            Index nextv;
            if (!cand.empty()) {
                nextv = cand.begin()->second;
                cand.erase(cand.begin());
            } else {
                while (cursor < nv && cur[cursor] == 0) ++cursor;
                if (cursor >= nv) break;   // everything grown (defensive)
                nextv = static_cast<Index>(cursor);
            }
            const std::size_t un = static_cast<std::size_t>(nextv);
            cur[un] = 0;
            wA += g.vwt[un];
            grown = nextv;
        }
        const Index cut = ndml_cut_of(g, cur);
        const Index imb = (W >= Index(2) * wA) ? (W - Index(2) * wA)
                                               : (Index(2) * wA - W);
        if (!have_best || cut < best_cut ||
            (cut == best_cut && imb < best_imb)) {
            have_best = true;
            best_cut = cut;
            best_imb = imb;
            side = cur;
        }
    }
    if (!have_best) side.assign(nv, char(0));   // nv == 0 (defensive)
}

// ---------------------------------------------------------------------------
// FM boundary refinement (design SS2-3).  `floor_w` is the weighted balance
// bound beta*W: a move is feasible only if it leaves its source side at or
// above floor_w -- except when a side is already BELOW floor_w, in which
// case only moves INTO the light side are considered (deterministic
// rebalance).  Each pass moves every vertex at most once and is rolled back
// to its lexicographically best (violation, cut) prefix; strict improvement
// only, so a no-gain graph is left untouched.
// ---------------------------------------------------------------------------
template <class Index>
void ndml_fm_refine(
    const ndml_graph<Index>& g,
    std::vector<char>& side,
    const Index floor_w,
    const int passes)
{
    const std::size_t nv = static_cast<std::size_t>(g.nv);
    if (nv == 0u) return;

    std::vector<Index> wdeg(nv, Index(0));
    for (std::size_t u = 0; u < nv; ++u) {
        const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
        const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
        for (std::size_t p = b; p < e; ++p) wdeg[u] += g.ewt[p];
    }

    std::vector<Index> ext(nv);
    std::vector<char> locked(nv);
    std::vector<Index> moves;
    moves.reserve(nv);

    for (int pass = 0; pass < passes; ++pass) {
        Index sw[2] = { Index(0), Index(0) };
        for (std::size_t u = 0; u < nv; ++u) {
            sw[static_cast<std::size_t>(side[u])] += g.vwt[u];
        }
        std::fill(ext.begin(), ext.end(), Index(0));
        for (std::size_t u = 0; u < nv; ++u) {
            const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
            const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
            for (std::size_t p = b; p < e; ++p) {
                if (side[u] != side[static_cast<std::size_t>(g.adjncy[p])]) {
                    ext[u] += g.ewt[p];
                }
            }
        }
        std::fill(locked.begin(), locked.end(), char(0));
        // cand keyed (-gain, v) over unlocked BOUNDARY vertices (ext > 0).
        std::set<std::pair<Index, Index> > cand;
        for (std::size_t u = 0; u < nv; ++u) {
            if (ext[u] > Index(0)) {
                cand.insert(std::make_pair(
                    wdeg[u] - Index(2) * ext[u], static_cast<Index>(u)));
            }
        }
        Index cur_cut = ndml_cut_of(g, side);
        const Index viol0 =
            (floor_w > sw[0] ? floor_w - sw[0] : Index(0)) +
            (floor_w > sw[1] ? floor_w - sw[1] : Index(0));
        Index best_viol = viol0, best_cut = cur_cut;
        std::size_t best_prefix = 0u;
        moves.clear();

        while (moves.size() < nv) {
            const Index light =
                (sw[0] < floor_w) ? Index(0)
                                  : ((sw[1] < floor_w) ? Index(1) : Index(-1));
            Index pick = Index(-1);
            if (light < Index(0)) {
                // balanced: best feasible boundary move (max gain, tie min
                // index); feasible = source side stays at or above floor_w.
                // The scan over infeasible candidates is capped (256) so a
                // pathological all-infeasible prefix cannot go quadratic;
                // the cap is deterministic (it only ends the pass earlier).
                std::size_t scanned = 0u;
                for (typename std::set<std::pair<Index, Index> >::const_iterator
                         it = cand.begin();
                     it != cand.end() && scanned < 256u; ++it, ++scanned) {
                    const std::size_t v = static_cast<std::size_t>(it->second);
                    if (sw[static_cast<std::size_t>(side[v])] - g.vwt[v]
                            >= floor_w) {
                        pick = it->second;
                        break;
                    }
                }
                if (pick < Index(0)) break;   // no feasible move: pass ends
            } else {
                // rebalance: best unlocked move INTO the light side (gain
                // order among boundary candidates first, then a
                // deterministic linear scan for non-boundary vertices).
                std::size_t scanned = 0u;
                for (typename std::set<std::pair<Index, Index> >::const_iterator
                         it = cand.begin();
                     it != cand.end() && scanned < 256u; ++it, ++scanned) {
                    const std::size_t v = static_cast<std::size_t>(it->second);
                    if (static_cast<Index>(side[v]) != light) {
                        pick = it->second;
                        break;
                    }
                }
                if (pick < Index(0)) {
                    Index bg = Index(0);
                    for (std::size_t v = 0; v < nv; ++v) {
                        if (locked[v] || static_cast<Index>(side[v]) == light) {
                            continue;
                        }
                        const Index gain = Index(2) * ext[v] - wdeg[v];
                        if (pick < Index(0) || gain > bg) {
                            pick = static_cast<Index>(v);
                            bg = gain;
                        }
                    }
                    if (pick < Index(0)) break;   // heavy side all locked
                }
            }

            const std::size_t v = static_cast<std::size_t>(pick);
            if (ext[v] > Index(0)) {
                cand.erase(std::make_pair(
                    wdeg[v] - Index(2) * ext[v], pick));
            }
            const Index gain = Index(2) * ext[v] - wdeg[v];
            sw[static_cast<std::size_t>(side[v])] -= g.vwt[v];
            side[v] = static_cast<char>(1 - side[v]);
            sw[static_cast<std::size_t>(side[v])] += g.vwt[v];
            locked[v] = 1;
            cur_cut -= gain;
            ext[v] = wdeg[v] - ext[v];
            moves.push_back(pick);

            const std::size_t b = static_cast<std::size_t>(g.xadj[v]);
            const std::size_t e = static_cast<std::size_t>(g.xadj[v + 1u]);
            for (std::size_t p = b; p < e; ++p) {
                const Index u = g.adjncy[p];
                const std::size_t uu = static_cast<std::size_t>(u);
                if (locked[uu]) continue;
                if (ext[uu] > Index(0)) {
                    cand.erase(std::make_pair(
                        wdeg[uu] - Index(2) * ext[uu], u));
                }
                if (side[uu] == side[v]) ext[uu] -= g.ewt[p];
                else                     ext[uu] += g.ewt[p];
                if (ext[uu] > Index(0)) {
                    cand.insert(std::make_pair(
                        wdeg[uu] - Index(2) * ext[uu], u));
                }
            }

            const Index viol =
                (floor_w > sw[0] ? floor_w - sw[0] : Index(0)) +
                (floor_w > sw[1] ? floor_w - sw[1] : Index(0));
            if (viol < best_viol ||
                (viol == best_viol && cur_cut < best_cut)) {
                best_viol = viol;
                best_cut = cur_cut;
                best_prefix = moves.size();
            }
        }

        // roll back to the best prefix
        for (std::size_t t = moves.size(); t > best_prefix; --t) {
            const std::size_t v = static_cast<std::size_t>(moves[t - 1u]);
            side[v] = static_cast<char>(1 - side[v]);
        }
        if (best_prefix == 0u && best_viol == viol0) break;  // pass changed nothing
    }
}

// ---------------------------------------------------------------------------
// Vertex separator from the edge cut (design SS2-4): greedy minimum-vertex-
// cover approximation on the cut edges.  Repeatedly takes the vertex covering
// the most uncovered cut edges; ties prefer the side from which FEWER
// separator vertices have been taken so far (the design's "alternate the
// smaller side" balancing), then the smallest index.  Runs on the FINEST
// piece graph (unit edge weights).
// ---------------------------------------------------------------------------
template <class Index>
void ndml_vertex_separator(
    const ndml_graph<Index>& g,
    const std::vector<char>& side,
    std::vector<char>& in_sep)
{
    const std::size_t nv = static_cast<std::size_t>(g.nv);
    in_sep.assign(nv, char(0));
    std::vector<Index> cutdeg(nv, Index(0));
    for (std::size_t u = 0; u < nv; ++u) {
        const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
        const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
        for (std::size_t p = b; p < e; ++p) {
            if (side[u] != side[static_cast<std::size_t>(g.adjncy[p])]) {
                ++cutdeg[u];
            }
        }
    }
    std::set<std::pair<Index, Index> > cand;   // (-cutdeg, v)
    for (std::size_t u = 0; u < nv; ++u) {
        if (cutdeg[u] > Index(0)) {
            cand.insert(std::make_pair(-cutdeg[u], static_cast<Index>(u)));
        }
    }
    Index taken[2] = { Index(0), Index(0) };
    while (!cand.empty()) {
        // walk the run of maximum cut-degree entries; prefer the side with
        // fewer separator vertices so far, then the smallest index (the set
        // ascends in index inside a run, so the first hit per side wins).
        typename std::set<std::pair<Index, Index> >::const_iterator it =
            cand.begin();
        const Index topkey = it->first;
        Index pick = it->second;
        if (taken[0] != taken[1]) {
            const char want = (taken[0] < taken[1]) ? char(0) : char(1);
            for (; it != cand.end() && it->first == topkey; ++it) {
                if (side[static_cast<std::size_t>(it->second)] == want) {
                    pick = it->second;
                    break;
                }
            }
        }
        const std::size_t v = static_cast<std::size_t>(pick);
        cand.erase(std::make_pair(-cutdeg[v], pick));
        in_sep[v] = 1;
        ++taken[static_cast<std::size_t>(side[v])];
        cutdeg[v] = Index(0);
        const std::size_t b = static_cast<std::size_t>(g.xadj[v]);
        const std::size_t e = static_cast<std::size_t>(g.xadj[v + 1u]);
        for (std::size_t p = b; p < e; ++p) {
            const std::size_t u = static_cast<std::size_t>(g.adjncy[p]);
            if (in_sep[u] || side[u] == side[v]) continue;
            if (cutdeg[u] <= Index(0)) continue;
            cand.erase(std::make_pair(-cutdeg[u], static_cast<Index>(u)));
            --cutdeg[u];
            if (cutdeg[u] > Index(0)) {
                cand.insert(std::make_pair(-cutdeg[u], static_cast<Index>(u)));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Direct vertex-separator refinement (node-FM, Ashcraft-Liu style -- the
// separator-side counterpart of the METIS-family standard construction).
// part[v] in {0 = A, 1 = B, 2 = S}.  Moving a separator vertex v into side
// d expels it from S but pulls every neighbour of v on the OTHER side into
// S, so the separator-size gain is w(v) - w(N(v) & other side).  Passes are
// FM-shaped: max-gain feasible move (tie smallest index, then side A), each
// vertex moved at most once per pass, best-(violation, |S|)-prefix kept by
// snapshot, strict improvement only.  Runs on the FINEST piece graph (unit
// weights).  Deterministic (D-2).
// ---------------------------------------------------------------------------
template <class Index>
void ndml_separator_refine(
    const ndml_graph<Index>& g,
    std::vector<char>& part,
    const Index floor_w,
    const int passes)
{
    const std::size_t nv = static_cast<std::size_t>(g.nv);
    if (nv == 0u) return;

    std::vector<Index> cnt[2];        // neighbours of v in side 0 / side 1
    cnt[0].assign(nv, Index(0));
    cnt[1].assign(nv, Index(0));
    std::vector<char> locked(nv);
    std::vector<char> best_part;
    // key: ((-gain, v), dir) -- max gain first, tie smallest index, then A.
    typedef std::pair<std::pair<Index, Index>, Index> nkey;

    for (int pass = 0; pass < passes; ++pass) {
        Index sz[3] = { Index(0), Index(0), Index(0) };
        for (std::size_t u = 0; u < nv; ++u) {
            sz[static_cast<std::size_t>(part[u])] += g.vwt[u];
        }
        for (std::size_t u = 0; u < nv; ++u) {
            cnt[0][u] = Index(0);
            cnt[1][u] = Index(0);
            const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
            const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
            for (std::size_t p = b; p < e; ++p) {
                const std::size_t w = static_cast<std::size_t>(g.adjncy[p]);
                if (part[w] < char(2)) {
                    cnt[static_cast<std::size_t>(part[w])][u] += g.vwt[w];
                }
            }
        }
        std::fill(locked.begin(), locked.end(), char(0));
        std::set<nkey> cand;
        for (std::size_t u = 0; u < nv; ++u) {
            if (part[u] != char(2)) continue;
            for (Index d = Index(0); d < Index(2); ++d) {
                cand.insert(nkey(std::make_pair(
                    cnt[static_cast<std::size_t>(Index(1) - d)][u] - g.vwt[u],
                    static_cast<Index>(u)), d));
            }
        }
        const Index viol0 =
            (floor_w > sz[0] ? floor_w - sz[0] : Index(0)) +
            (floor_w > sz[1] ? floor_w - sz[1] : Index(0));
        Index best_viol = viol0, best_ssz = sz[2];
        best_part = part;
        bool improved = false;
        std::size_t moved = 0u;

        while (moved < nv && !cand.empty()) {
            // best feasible move: expelling N(v) on the other side from that
            // side must keep it at or above floor_w (scan capped like the
            // edge-FM: deterministic, only ends the pass earlier).
            Index pick = Index(-1), dir = Index(0);
            std::size_t scanned = 0u;
            for (typename std::set<nkey>::const_iterator it = cand.begin();
                 it != cand.end() && scanned < 256u; ++it, ++scanned) {
                const std::size_t v = static_cast<std::size_t>(it->first.second);
                const std::size_t od = static_cast<std::size_t>(Index(1) - it->second);
                if (sz[od] - cnt[od][v] >= floor_w) {
                    pick = it->first.second;
                    dir = it->second;
                    break;
                }
            }
            if (pick < Index(0)) break;
            const std::size_t v = static_cast<std::size_t>(pick);
            const std::size_t d = static_cast<std::size_t>(dir);
            const std::size_t od = 1u - d;

            // remove both directional entries of v
            cand.erase(nkey(std::make_pair(cnt[1u - 0u][v] - g.vwt[v], pick), Index(0)));
            cand.erase(nkey(std::make_pair(cnt[0u][v] - g.vwt[v], pick), Index(1)));

            // move v into side d; pull its other-side neighbours into S
            part[v] = static_cast<char>(d);
            sz[2] -= g.vwt[v];
            sz[d] += g.vwt[v];
            locked[v] = 1;
            ++moved;
            std::vector<Index> changed;   // vertices whose part changed
            changed.push_back(pick);
            {
                const std::size_t b = static_cast<std::size_t>(g.xadj[v]);
                const std::size_t e = static_cast<std::size_t>(g.xadj[v + 1u]);
                for (std::size_t p = b; p < e; ++p) {
                    const Index u = g.adjncy[p];
                    const std::size_t uu = static_cast<std::size_t>(u);
                    if (part[uu] == static_cast<char>(od)) {
                        if (part[uu] < char(2)) {
                            sz[static_cast<std::size_t>(part[uu])] -= g.vwt[uu];
                        }
                        part[uu] = 2;
                        sz[2] += g.vwt[uu];
                        changed.push_back(u);
                    }
                }
            }
            // recompute the neighbour counts and candidate entries of every
            // separator vertex adjacent to a changed vertex (and of the
            // changed vertices themselves when they are now in S).
            std::vector<Index> touch = changed;
            for (std::size_t t = 0; t < changed.size(); ++t) {
                const std::size_t u = static_cast<std::size_t>(changed[t]);
                const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
                const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
                for (std::size_t p = b; p < e; ++p) touch.push_back(g.adjncy[p]);
            }
            std::sort(touch.begin(), touch.end());
            touch.erase(std::unique(touch.begin(), touch.end()), touch.end());
            for (std::size_t t = 0; t < touch.size(); ++t) {
                const std::size_t u = static_cast<std::size_t>(touch[t]);
                // drop stale candidate entries keyed with the OLD counts
                cand.erase(nkey(std::make_pair(cnt[1][u] - g.vwt[u],
                                               touch[t]), Index(0)));
                cand.erase(nkey(std::make_pair(cnt[0][u] - g.vwt[u],
                                               touch[t]), Index(1)));
                cnt[0][u] = Index(0);
                cnt[1][u] = Index(0);
                const std::size_t b = static_cast<std::size_t>(g.xadj[u]);
                const std::size_t e = static_cast<std::size_t>(g.xadj[u + 1u]);
                for (std::size_t p = b; p < e; ++p) {
                    const std::size_t w = static_cast<std::size_t>(g.adjncy[p]);
                    if (part[w] < char(2)) {
                        cnt[static_cast<std::size_t>(part[w])][u] += g.vwt[w];
                    }
                }
                if (part[u] == char(2) && !locked[u]) {
                    cand.insert(nkey(std::make_pair(cnt[1][u] - g.vwt[u],
                                                    touch[t]), Index(0)));
                    cand.insert(nkey(std::make_pair(cnt[0][u] - g.vwt[u],
                                                    touch[t]), Index(1)));
                }
            }

            const Index viol =
                (floor_w > sz[0] ? floor_w - sz[0] : Index(0)) +
                (floor_w > sz[1] ? floor_w - sz[1] : Index(0));
            if (viol < best_viol ||
                (viol == best_viol && sz[2] < best_ssz)) {
                best_viol = viol;
                best_ssz = sz[2];
                best_part = part;
                improved = true;
            }
        }
        part = best_part;
        if (!improved) break;
    }
}

// ---------------------------------------------------------------------------
// Multilevel bisection of one piece: coarsen to at most coarsen_stop
// vertices, bisect the coarsest graph, uncoarsen with FM per level, lift the
// edge cut to a vertex separator on the finest piece graph.  Splits `verts`
// (ascending) into left / right / sep (each ascending).  Returns false when
// the piece could not be split (one half empty AND no separator) -- the
// caller then falls back to AMD-ordering the whole piece.
// ---------------------------------------------------------------------------
template <class Index>
bool ndml_bisect_piece(
    const std::vector<Index>& verts,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    std::vector<Index>& loc,
    const sparse_order_ndml_params& prm,
    std::vector<Index>& left,
    std::vector<Index>& right,
    std::vector<Index>& sep)
{
    left.clear(); right.clear(); sep.clear();

    std::vector<ndml_graph<Index> > levels(1u);
    ndml_extract_subgraph(verts, adj_ptr, adj_ind, loc, levels[0]);

    long long stop = prm.coarsen_stop;
    if (stop < 8) stop = 8;
    while (static_cast<long long>(levels.back().nv) > stop) {
        levels.push_back(ndml_graph<Index>());
        if (!ndml_coarsen_step(levels[levels.size() - 2u], levels.back())) {
            levels.pop_back();
            break;
        }
    }

    // total weight and the weighted balance floor beta*W (integer percent)
    Index W = Index(0);
    for (std::size_t i = 0; i < levels[0].vwt.size(); ++i) W += levels[0].vwt[i];
    int pct = prm.balance_pct;
    if (pct < 5) pct = 5;
    if (pct > 49) pct = 49;
    const Index floor_w = (W * static_cast<Index>(pct)) / Index(100);
    int passes = prm.fm_passes;
    if (passes < 0) passes = 0;
    if (passes > 8) passes = 8;

    std::vector<char> side;
    ndml_initial_bisection(levels.back(), side);
    ndml_fm_refine(levels.back(), side, floor_w, passes);
    for (std::size_t l = levels.size() - 1u; l > 0u; --l) {
        const ndml_graph<Index>& fine = levels[l - 1u];
        std::vector<char> fside(static_cast<std::size_t>(fine.nv));
        for (std::size_t u = 0; u < fside.size(); ++u) {
            fside[u] = side[static_cast<std::size_t>(fine.cmap[u])];
        }
        side.swap(fside);
        ndml_fm_refine(fine, side, floor_w, passes);
    }

    std::vector<char> in_sep;
    ndml_vertex_separator(levels[0], side, in_sep);

    // node-FM refinement of the vertex separator itself (S in/out moves)
    std::vector<char> part(side.size());
    for (std::size_t u = 0; u < side.size(); ++u) {
        part[u] = in_sep[u] ? char(2) : side[u];
    }
    ndml_separator_refine(levels[0], part, floor_w, passes);

    // NESDIS nd_oksep guard: a separator above max_sep_pct percent of the
    // piece has no dissection value (expander-shaped piece) -- reject the
    // split; the caller orders the whole piece with the halo-augmented AMD.
    {
        Index sep_w = Index(0);
        for (std::size_t u = 0; u < part.size(); ++u) {
            if (part[u] == char(2)) sep_w += levels[0].vwt[u];
        }
        int sp = prm.max_sep_pct;
        if (sp < 1) sp = 1;
        if (sp > 100) sp = 100;
        if (sep_w * Index(100) > W * static_cast<Index>(sp)) return false;
    }

    for (std::size_t u = 0; u < part.size(); ++u) {
        if (part[u] == char(2)) sep.push_back(verts[u]);
        else if (part[u] == char(0)) left.push_back(verts[u]);
        else right.push_back(verts[u]);
    }
    // no progress = would recurse on the identical piece: report failure
    if (sep.empty() && (left.empty() || right.empty())) return false;
    return true;
}

// ---------------------------------------------------------------------------
// AMD-order one piece and append to `order` -- the constrained-AMD simple
// form (design SS2-5), HALO-AUGMENTED: the subgraph handed to AMD is the
// piece PLUS its one-ring of outside neighbours (separator vertices of the
// enclosing dissection levels, all numbered after the piece).  AMD sees the
// halo's coupling, so piece vertices adjacent to separators keep their true
// (high) degree and are eliminated late; the emitted order is the AMD order
// RESTRICTED to the piece.  Without the halo, a banded graph degenerates:
// the bare subgraph's boundary vertices look low-degree, get eliminated
// first, and drag the separator coupling through the whole piece (measured
// 1.89x on G5 -- the OR-1 no-degradation failure this fixes).
//
// The induced adjacency is a valid CSC pattern input for the AMD entry (it
// symmetrizes and de-duplicates internally; here it is already symmetric,
// self-loop-free and ascending).
// ---------------------------------------------------------------------------
template <class Index>
void ndml_order_leaf_amd(
    const std::vector<Index>& verts,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    std::vector<Index>& loc,
    std::vector<Index>& order)
{
    if (verts.empty()) return;
    if (verts.size() == 1u) { order.push_back(verts[0]); return; }

    // halo = outside one-ring (deterministic: sorted, unique)
    for (std::size_t i = 0; i < verts.size(); ++i) {
        loc[static_cast<std::size_t>(verts[i])] = Index(1);
    }
    std::vector<Index> halo;
    for (std::size_t i = 0; i < verts.size(); ++i) {
        const std::size_t b = static_cast<std::size_t>(
            adj_ptr[static_cast<std::size_t>(verts[i])]);
        const std::size_t e = static_cast<std::size_t>(
            adj_ptr[static_cast<std::size_t>(verts[i]) + 1u]);
        for (std::size_t p = b; p < e; ++p) {
            const Index w = adj_ind[p];
            if (loc[static_cast<std::size_t>(w)] == Index(-1)) {
                loc[static_cast<std::size_t>(w)] = Index(-2);   // halo mark
                halo.push_back(w);
            }
        }
    }
    std::sort(halo.begin(), halo.end());
    for (std::size_t i = 0; i < verts.size(); ++i) {
        loc[static_cast<std::size_t>(verts[i])] = Index(-1);
    }
    for (std::size_t i = 0; i < halo.size(); ++i) {
        loc[static_cast<std::size_t>(halo[i])] = Index(-1);
    }

    std::vector<Index> aug(verts.size() + halo.size());
    std::merge(verts.begin(), verts.end(), halo.begin(), halo.end(),
               aug.begin());

    ndml_graph<Index> g;
    ndml_extract_subgraph(aug, adj_ptr, adj_ind, loc, g);
    const std::vector<Index> lperm =
        sparse_lu_amd_ordering(g.nv, g.xadj, g.adjncy);
    for (std::size_t k = 0; k < lperm.size(); ++k) {
        const Index v = aug[static_cast<std::size_t>(lperm[k])];
        // emit only piece members (binary search: verts is ascending)
        if (std::binary_search(verts.begin(), verts.end(), v)) {
            order.push_back(v);
        }
    }
}

// ---------------------------------------------------------------------------
// Exact symbolic off-diagonal fill count of the Cholesky/LDL factor for the
// symmetric pattern graph under permutation `perm` (perm[new] = old), used by
// the best-of-two selection.  Standard elimination-tree row-subtree
// traversal: for each row i (new order), walking each below-diagonal
// neighbour up the partially built etree until an already-marked node visits
// exactly the nonzero columns of row i of L, so the total work -- and the
// returned count -- is exactly nnz(L) minus the diagonal.  Deterministic.
// ---------------------------------------------------------------------------
template <class Index>
long long ndml_fill_count(Index n,
                          const std::vector<Index>& adj_ptr,
                          const std::vector<Index>& adj_ind,
                          const std::vector<Index>& perm)
{
    const std::size_t un = static_cast<std::size_t>(n);
    std::vector<Index> pinv(un), parent(un, Index(-1)), mark(un, Index(-1));
    for (std::size_t k = 0; k < un; ++k) {
        pinv[static_cast<std::size_t>(perm[k])] = static_cast<Index>(k);
    }
    long long count = 0;
    for (Index i = Index(0); i < n; ++i) {
        const Index old = perm[static_cast<std::size_t>(i)];
        for (Index p = adj_ptr[static_cast<std::size_t>(old)];
             p < adj_ptr[static_cast<std::size_t>(old) + 1]; ++p) {
            Index j = pinv[static_cast<std::size_t>(adj_ind[static_cast<std::size_t>(p)])];
            while (j < i && mark[static_cast<std::size_t>(j)] != i) {
                mark[static_cast<std::size_t>(j)] = i;
                ++count;                       // L(i, j) is a nonzero
                if (parent[static_cast<std::size_t>(j)] == Index(-1)) {
                    parent[static_cast<std::size_t>(j)] = i;
                }
                j = parent[static_cast<std::size_t>(j)];
            }
        }
    }
    return count;
}

} // namespace sparse_ndml_detail

// ---------------------------------------------------------------------------
// sparse_lu_nested_dissection_ml_ordering -- the ORD-1 public entry.
//
// Returns a permutation perm with perm[new] = old (same convention as
// rcm / amd / colamd / nested_dissection).  Pattern-only (S-1),
// deterministic (S-5 / ruling D-2), explicit work stack (no recursion).
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index> sparse_lu_nested_dissection_ml_ordering(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const sparse_order_ndml_params& prm = sparse_order_ndml_params())
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_nested_dissection_ml_ordering: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> perm(un);
    if (n <= Index(1)) {
        for (Index i = Index(0); i < n; ++i) perm[static_cast<std::size_t>(i)] = i;
        return perm;
    }

    std::vector<Index> adj_ptr, adj_ind;
    sparse_lu_build_symmetric_pattern_graph(n, col_ptr, row_ind, adj_ptr, adj_ind);

    long long leaf = prm.leaf_size;
    if (leaf < 1) leaf = 1;

    std::vector<Index> loc(un, Index(-1));   // shared global->local scratch

    // Explicit work stack.  An ORDER task dissects its piece; an EMIT task
    // appends its (pre-sorted) separator block.  Children are pushed right-
    // to-left so the pop order is left, right, separator -- separators are
    // numbered last (design SS2).
    struct task {
        std::vector<Index> verts;
        bool emit_only;
    };
    std::vector<task> stack;
    std::vector<Index> order;
    order.reserve(un);

    stack.push_back(task());
    stack.back().emit_only = false;
    stack.back().verts.resize(un);
    for (std::size_t i = 0; i < un; ++i) {
        stack.back().verts[i] = static_cast<Index>(i);
    }

    std::vector<Index> left, right, sep;
    while (!stack.empty()) {
        task t;
        t.verts.swap(stack.back().verts);
        t.emit_only = stack.back().emit_only;
        stack.pop_back();

        if (t.emit_only) {
            for (std::size_t i = 0; i < t.verts.size(); ++i) {
                order.push_back(t.verts[i]);
            }
            continue;
        }
        if (static_cast<long long>(t.verts.size()) <= leaf) {
            sparse_ndml_detail::ndml_order_leaf_amd(
                t.verts, adj_ptr, adj_ind, loc, order);
            continue;
        }
        if (!sparse_ndml_detail::ndml_bisect_piece(
                t.verts, adj_ptr, adj_ind, loc, prm, left, right, sep)) {
            sparse_ndml_detail::ndml_order_leaf_amd(
                t.verts, adj_ptr, adj_ind, loc, order);
            continue;
        }
        if (!sep.empty()) {
            stack.push_back(task());
            stack.back().emit_only = true;
            stack.back().verts.swap(sep);       // ascending by construction
        }
        if (!right.empty()) {
            stack.push_back(task());
            stack.back().emit_only = false;
            stack.back().verts.swap(right);
        }
        if (!left.empty()) {
            stack.push_back(task());
            stack.back().emit_only = false;
            stack.back().verts.swap(left);
        }
    }

    if (order.size() != un) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_nested_dissection_ml_ordering: order did not cover all vertices");
    }
    for (std::size_t i = 0; i < un; ++i) perm[i] = order[i];

    // Best-of-two selection (see sparse_order_ndml_params::best_of_two):
    // compare the exact symbolic fill of the ND order against plain AMD on
    // the same pattern graph and keep the better order; ties keep ND.
    if (prm.best_of_two != 0) {
        const std::vector<Index> amd_perm =
            sparse_lu_amd_ordering(n, adj_ptr, adj_ind);
        const long long fill_nd =
            sparse_ndml_detail::ndml_fill_count(n, adj_ptr, adj_ind, perm);
        const long long fill_amd =
            sparse_ndml_detail::ndml_fill_count(n, adj_ptr, adj_ind, amd_perm);
        if (fill_amd < fill_nd) {
            perm = amd_perm;
        }
    }
    return perm;
}

#endif // VCP_TSPARSE_ORDER_NDML_IMPL_HPP
