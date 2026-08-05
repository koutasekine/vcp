// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// Ordering Track O1 — RCM fill-reducing ordering — internal implementation.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER csc_storage<>
// and the SLU-1 conversion helpers (sparse_lu_inverse_permutation, ...) are in
// scope.  It has no "namespace vcp { }" wrapper; it is injected by
// tsparse_sparse_lu.hpp.
//
// Do NOT include this file directly.  Include one of:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>
//   <vcp/tsparse/tsparse.hpp>                    (umbrella)
//
// Safety invariants (roadmap rev2 §3.5):
//   S-1 Integer-only / pattern-only.  No numeric values are read.  The graph is
//       built from the sparsity pattern of A + A^T only, so this cannot affect
//       the point-arithmetic residual contract.
//   S-5 Deterministic.  Fixed tie-break (ascending degree, then ascending
//       index).  No randomness.
//
// The result is a COLUMN permutation Q with convention col_perm[new] = old,
// installed into the symbolic col_perm slot (roadmap rev2 §2A).  Numeric
// threshold partial pivoting and row_perm are untouched.

#ifndef VCP_TSPARSE_SPARSE_LU_ORDERING_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_ORDERING_IMPL_HPP

// NOTE: this file is textually injected INSIDE namespace vcp; every standard
// header it needs must already be included by the injecting header BEFORE the
// namespace opens (tsparse_sparse_lu.hpp includes <set> for SLU-OQ1).  The
// includes below are no-op guards when that discipline is followed.
#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/error.hpp>

// ---------------------------------------------------------------------------
// O1.1  Undirected structural graph from the pattern of A + A^T.
//
// Pattern-only: reads col_ptr / row_ind, never values (S-1).  Self-loops
// (diagonal entries) are dropped.  Parallel edges are de-duplicated.  Isolated
// vertices (empty row and column) appear as zero-degree vertices.  The result
// is a symmetric adjacency stored in CSR-like form (adj_ptr / adj_ind), with
// neighbour lists sorted ascending by index for deterministic traversal.
//
// col_ptr / row_ind describe an n x n matrix in CSC: for column c, the rows
// row_ind[col_ptr[c] .. col_ptr[c+1]) are the nonzero rows.  Edge {r, c} is
// added for every off-diagonal nonzero (this symmetrises A and A^T at once,
// because adding both endpoints of {r, c} covers the transpose pattern).
// ---------------------------------------------------------------------------
template <class Index>
void sparse_lu_build_symmetric_pattern_graph(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    std::vector<Index>& adj_ptr,
    std::vector<Index>& adj_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_build_symmetric_pattern_graph: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    // Pass 1: count upper bound on degree (each off-diagonal nonzero contributes
    // one edge to each endpoint; duplicates removed later).
    std::vector<Index> count(un, Index(0));
    for (Index c = Index(0); c < n; ++c) {
        const Index kb = col_ptr[static_cast<std::size_t>(c)];
        const Index ke = col_ptr[static_cast<std::size_t>(c) + 1u];
        for (Index k = kb; k < ke; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r == c) continue;                 // drop self-loop
            if (r < Index(0) || r >= n) continue; // defensive; convert validated
            ++count[static_cast<std::size_t>(c)];
            ++count[static_cast<std::size_t>(r)];
        }
    }

    adj_ptr.assign(un + 1u, Index(0));
    for (std::size_t i = 0; i < un; ++i) {
        adj_ptr[i + 1u] = adj_ptr[i] + count[i];
    }
    const std::size_t total = static_cast<std::size_t>(adj_ptr[un]);
    std::vector<Index> tmp_ind(total);
    std::vector<Index> head(un);
    for (std::size_t i = 0; i < un; ++i) {
        head[i] = adj_ptr[i];
    }

    // Pass 2: scatter both endpoints of every off-diagonal edge.
    for (Index c = Index(0); c < n; ++c) {
        const Index kb = col_ptr[static_cast<std::size_t>(c)];
        const Index ke = col_ptr[static_cast<std::size_t>(c) + 1u];
        for (Index k = kb; k < ke; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r == c) continue;
            if (r < Index(0) || r >= n) continue;
            tmp_ind[static_cast<std::size_t>(head[static_cast<std::size_t>(c)]++)] = r;
            tmp_ind[static_cast<std::size_t>(head[static_cast<std::size_t>(r)]++)] = c;
        }
    }

    // Pass 3: sort + de-duplicate each neighbour list (deterministic ascending).
    adj_ind.clear();
    adj_ind.reserve(total);
    std::vector<Index> new_ptr(un + 1u, Index(0));
    for (std::size_t i = 0; i < un; ++i) {
        const std::size_t b = static_cast<std::size_t>(adj_ptr[i]);
        const std::size_t e = static_cast<std::size_t>(adj_ptr[i + 1u]);
        std::sort(tmp_ind.begin() + b, tmp_ind.begin() + e);
        Index prev = Index(-1);
        for (std::size_t p = b; p < e; ++p) {
            const Index w = tmp_ind[p];
            if (w == prev) continue;  // de-duplicate parallel edges
            adj_ind.push_back(w);
            prev = w;
        }
        new_ptr[i + 1u] = static_cast<Index>(adj_ind.size());
    }
    adj_ptr.swap(new_ptr);
}

// ---------------------------------------------------------------------------
// O1.2  Rooted level structure (BFS) within the component containing root.
//
// Fills `levels` so that levels[d] lists, in ascending index order, the
// vertices at BFS distance d from root.  Returns the eccentricity (number of
// levels - 1).  `seen` is a scratch buffer of size n that is reset for every
// vertex this call touches (so it stays cheap across repeated calls).
// ---------------------------------------------------------------------------
template <class Index>
Index sparse_lu_rooted_level_structure(
    Index root,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    std::vector<char>& seen,
    std::vector<std::vector<Index> >& levels)
{
    levels.clear();
    std::vector<Index> touched;

    std::vector<Index> frontier;
    frontier.push_back(root);
    seen[static_cast<std::size_t>(root)] = 1;
    touched.push_back(root);

    while (!frontier.empty()) {
        levels.push_back(frontier);  // already ascending (root, or sorted below)
        std::vector<Index> next;
        for (std::size_t i = 0; i < frontier.size(); ++i) {
            const Index u = frontier[i];
            const std::size_t b = static_cast<std::size_t>(adj_ptr[static_cast<std::size_t>(u)]);
            const std::size_t e = static_cast<std::size_t>(adj_ptr[static_cast<std::size_t>(u) + 1u]);
            for (std::size_t p = b; p < e; ++p) {
                const Index w = adj_ind[p];
                if (!seen[static_cast<std::size_t>(w)]) {
                    seen[static_cast<std::size_t>(w)] = 1;
                    touched.push_back(w);
                    next.push_back(w);
                }
            }
        }
        std::sort(next.begin(), next.end());  // deterministic level membership
        frontier.swap(next);
    }

    // Reset only the entries this call set, keeping the scratch buffer reusable.
    for (std::size_t i = 0; i < touched.size(); ++i) {
        seen[static_cast<std::size_t>(touched[i])] = 0;
    }
    return static_cast<Index>(levels.size()) - Index(1);
}

// ---------------------------------------------------------------------------
// O1.3  Pseudo-peripheral start vertex (George-Liu), restricted to the
// component reachable from `seed`.  Deterministic: among the deepest level it
// picks the minimum-degree vertex, tie-broken by smallest index.  Falls back to
// `seed` if no improvement is found.  Bounded iteration count for safety.
// ---------------------------------------------------------------------------
template <class Index>
Index sparse_lu_pseudo_peripheral_start(
    Index seed,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    const std::vector<Index>& degree,
    std::vector<char>& seen)
{
    Index v = seed;
    std::vector<std::vector<Index> > levels;
    Index ecc = sparse_lu_rooted_level_structure(v, adj_ptr, adj_ind, seen, levels);

    const int max_iter = 20;  // George-Liu converges in very few iterations
    for (int it = 0; it < max_iter; ++it) {
        if (levels.empty()) break;
        const std::vector<Index>& last = levels.back();
        // Pick the minimum-degree vertex in the last level (tie: smallest index;
        // `last` is ascending, so first-min scan yields the smallest index).
        Index cand = last[0];
        Index cand_deg = degree[static_cast<std::size_t>(cand)];
        for (std::size_t i = 1; i < last.size(); ++i) {
            const Index w = last[i];
            const Index d = degree[static_cast<std::size_t>(w)];
            if (d < cand_deg) { cand = w; cand_deg = d; }
        }
        std::vector<std::vector<Index> > cand_levels;
        const Index cand_ecc = sparse_lu_rooted_level_structure(
            cand, adj_ptr, adj_ind, seen, cand_levels);
        if (cand_ecc > ecc) {
            v = cand;
            ecc = cand_ecc;
            levels.swap(cand_levels);
        } else {
            break;  // no eccentricity improvement: v is pseudo-peripheral
        }
    }
    return v;
}

// ---------------------------------------------------------------------------
// O1.4  Deterministic Reverse Cuthill-McKee ordering.
//
// Returns a column permutation Q with convention col_perm[new] = old (roadmap
// rev2 §2A).  Builds the A + A^T pattern graph, runs Cuthill-McKee over each
// connected component (component starts iterated by ascending index so every
// vertex is covered exactly once), then reverses the full order (the "R" of
// RCM).  For each visited vertex, unvisited neighbours are appended in
// ascending (degree, index) order (S-5).  n <= 1 returns the identity.
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index> sparse_lu_rcm_ordering(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_rcm_ordering: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> perm(un);
    if (n <= Index(1)) {
        for (Index i = Index(0); i < n; ++i) {
            perm[static_cast<std::size_t>(i)] = i;
        }
        return perm;
    }

    std::vector<Index> adj_ptr, adj_ind;
    sparse_lu_build_symmetric_pattern_graph(n, col_ptr, row_ind, adj_ptr, adj_ind);

    std::vector<Index> degree(un);
    for (std::size_t i = 0; i < un; ++i) {
        degree[i] = adj_ptr[i + 1u] - adj_ptr[i];
    }

    std::vector<char> visited(un, 0);
    std::vector<char> seen(un, 0);          // scratch for level-structure BFS
    std::vector<Index> order;               // Cuthill-McKee order (new index -> old)
    order.reserve(un);

    // Iterate component seeds by ascending index so every vertex is covered.
    for (Index s = Index(0); s < n; ++s) {
        if (visited[static_cast<std::size_t>(s)]) continue;

        const Index start =
            sparse_lu_pseudo_peripheral_start(s, adj_ptr, adj_ind, degree, seen);

        // BFS (Cuthill-McKee) within this component.  `qh` runs from the
        // component start to the current end of `order`; when the component is
        // drained qh == order.size() and we fall back to the outer loop.
        visited[static_cast<std::size_t>(start)] = 1;
        order.push_back(start);
        std::size_t qh = order.size() - 1u;
        while (qh < order.size()) {
            const Index u = order[qh];
            ++qh;
            const std::size_t b =
                static_cast<std::size_t>(adj_ptr[static_cast<std::size_t>(u)]);
            const std::size_t e =
                static_cast<std::size_t>(adj_ptr[static_cast<std::size_t>(u) + 1u]);
            // Gather unvisited neighbours of u.
            std::vector<Index> nbrs;
            nbrs.reserve(e - b);
            for (std::size_t p = b; p < e; ++p) {
                const Index w = adj_ind[p];
                if (!visited[static_cast<std::size_t>(w)]) {
                    nbrs.push_back(w);
                }
            }
            // Sort by ascending degree, tie-break ascending index (S-5).
            std::sort(nbrs.begin(), nbrs.end(),
                      [&degree](Index a, Index b2) {
                          const Index da = degree[static_cast<std::size_t>(a)];
                          const Index db = degree[static_cast<std::size_t>(b2)];
                          if (da != db) return da < db;
                          return a < b2;
                      });
            for (std::size_t i = 0; i < nbrs.size(); ++i) {
                const Index w = nbrs[i];
                if (!visited[static_cast<std::size_t>(w)]) {
                    visited[static_cast<std::size_t>(w)] = 1;
                    order.push_back(w);
                }
            }
        }
    }

    // Coverage invariant: every vertex visited exactly once.
    if (order.size() != un) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_rcm_ordering: component partition did not cover all vertices");
    }

    // Reverse the Cuthill-McKee order -> RCM.  perm[new] = old.
    for (std::size_t i = 0; i < un; ++i) {
        perm[i] = order[un - 1u - i];
    }
    return perm;
}

// ===========================================================================
// Ordering Track O2 — AMD (Approximate Minimum Degree) fill-reducing ordering.
//
// Reuses the Stage-1 A + A^T pattern graph (sparse_lu_build_symmetric_pattern_
// graph) and the same col_perm[new] = old convention; only the ordering routine
// itself is new.  AMD is a quotient-graph minimum-degree elimination on the
// symmetric pattern (Amestoy / Davis / Duff 1996) with the standard ingredients:
//   - approximate external degree (three upper bounds, take the minimum),
//   - mass elimination via supervariables,
//   - supervariable (indistinguishable-variable) detection,
//   - element absorption (both exact and aggressive |Le \ Lp| == 0).
//
// Safety invariants (roadmap rev2 §4.5 = §3.5 with opt.ordering = amd):
//   S-1 Integer-only / pattern-only.  Graph built from the sparsity pattern of
//       A + A^T only; no numeric values are read.
//   S-5 Deterministic.  Pivot = minimum approximate degree, tie-break smallest
//       index.  Supervariable principal = smallest index.  No randomness.
//   S-6 The result is a COLUMN permutation Q installed into the symbolic
//       col_perm slot; the answer is unchanged (cost/structure only).
//
// n <= 1 returns identity.  Disconnected components, empty rows/columns,
// zero-degree (isolated) vertices are handled.  A coverage check throws if the
// elimination did not cover every vertex exactly once.
// ===========================================================================

template <class Index>
inline void sparse_lu_amd_sorted_insert(std::vector<Index>& a, Index x)
{
    typename std::vector<Index>::iterator it =
        std::lower_bound(a.begin(), a.end(), x);
    if (it == a.end() || *it != x) a.insert(it, x);
}

template <class Index>
inline void sparse_lu_amd_sorted_remove_one(std::vector<Index>& a, Index x)
{
    typename std::vector<Index>::iterator it =
        std::lower_bound(a.begin(), a.end(), x);
    if (it != a.end() && *it == x) a.erase(it);
}

// a := a \ rm  (both sorted ascending, unique).
template <class Index>
inline void sparse_lu_amd_sorted_remove_set(std::vector<Index>& a,
                                            const std::vector<Index>& rm)
{
    if (rm.empty() || a.empty()) return;
    std::vector<Index> out;
    out.reserve(a.size());
    std::size_t i = 0, j = 0;
    while (i < a.size()) {
        while (j < rm.size() && rm[j] < a[i]) ++j;
        if (j < rm.size() && rm[j] == a[i]) { ++i; continue; }
        out.push_back(a[i]);
        ++i;
    }
    a.swap(out);
}

// ORD-F2 Phase A2: in-place variant of the set difference above, used by the
// amd section only (the colamd section keeps the allocating form).  Same
// two-pointer scan compacting into `a` itself (write cursor k <= read cursor
// i at all times), so the resulting content is identical while the per-call
// output allocation disappears.
template <class Index>
inline void sparse_lu_amd_sorted_remove_set_inplace_(
    std::vector<Index>& a,
    const std::vector<Index>& rm)
{
    if (rm.empty() || a.empty()) return;
    const std::size_t an = a.size(), rn = rm.size();
    std::size_t i = 0, j = 0, k = 0;
    while (i < an) {
        while (j < rn && rm[j] < a[i]) ++j;
        if (j < rn && rm[j] == a[i]) { ++i; continue; }
        a[k++] = a[i++];
    }
    a.resize(k);
}

template <class Index>
inline bool sparse_lu_amd_lists_equal(const std::vector<Index>& a,
                                      const std::vector<Index>& b)
{
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) if (a[i] != b[i]) return false;
    return true;
}

// ---------------------------------------------------------------------------
// ORD-F2 Phase A1: lazy-deletion binary-heap candidate PQ for O(log n) pivot
// selection, used by the colamd section (candq, ORD-F3 P2).  The SLU-OQ1
// std::set form it replaced was removed in F2-b Phase 5 together with the
// pre-F2-b amd body (a private verbatim copy lives in
// sandbox/probes/colprof.hpp).  Same total order as that std::set --
// minimum (degree, index) entry on top
// via std::push_heap / std::pop_heap with the "greater" comparator below
// (the ORD-F1 Phase 2 mechanism; a local functor instead of std::greater
// because this header is injected inside namespace vcp and must not pull in
// <functional> there).  Entries are pushed on every degree rewrite and never
// erased; a popped entry is valid iff
//     alive[v] && !is_element[v] && key == degree[v]
// (checked at the selection site).  A retire helper is unnecessary for the
// heap (ORD-F3 P2): eligibility-loss sites just flip alive[v] to 0, which
// stales every stored entry of v -- lazy deletion needs nothing else.
// Duplicate live (key, v) entries can exist when a degree returns to a
// former value; selection adopts exactly one entry per pivot and the adopted
// vertex immediately becomes an element (!alive && is_element), so leftover
// duplicates are invalidated -- no ORD-F1-style duplicate-scan stamp needed.
// ---------------------------------------------------------------------------
struct sparse_lu_ordering_pq_greater_ {
    template <class E>
    bool operator()(const E& a, const E& b) const { return b < a; }
};

// ORD-F3 P1: compaction guard for the lazy-deletion heap.  Rebuilds `cand`
// in place keeping ONLY the currently-valid entries
//     key == degree[v] && alive[v] && !is_element[v]
// and re-heapifies with the same comparator.  Byte-identity argument: the
// pop-time validity check already discards every stale entry, so dropping
// them here preserves the multiset of VALID entries exactly (equal-key
// duplicates of a still-valid vertex are kept) -- the sequence of adopted
// (degree, index) minima, hence the emitted permutation, is unchanged.
template <class Index>
inline void sparse_lu_ordering_heap_compact_(
    std::vector<std::pair<Index, Index> >& cand,
    const std::vector<Index>& degree,
    const std::vector<char>& alive,
    const std::vector<char>& is_element)
{
    std::size_t k = 0;
    for (std::size_t s = 0; s < cand.size(); ++s) {
        const std::size_t uv = static_cast<std::size_t>(cand[s].second);
        if (alive[uv] && !is_element[uv] && cand[s].first == degree[uv]) {
            cand[k++] = cand[s];
        }
    }
    cand.resize(k);
    std::make_heap(cand.begin(), cand.end(), sparse_lu_ordering_pq_greater_());
}

template <class Index>
inline void sparse_lu_ordering_cand_update_degree_(
    std::vector<std::pair<Index, Index> >& cand,
    std::vector<Index>& degree,
    const Index i,
    const Index new_degree,
    const std::vector<char>& alive,
    const std::vector<char>& is_element,
    const Index n_hint)
{
    degree[static_cast<std::size_t>(i)] = new_degree;
    cand.push_back(std::make_pair(new_degree, i));
    std::push_heap(cand.begin(), cand.end(), sparse_lu_ordering_pq_greater_());
    // ORD-F3 P1 compaction guard.  Internal constant 8n: 2x head-room over
    // the observed ~3.9n series peak (design §1 P1), so it never fires on
    // the normal series and exists purely as a hard O(n)-rebuild memory
    // bound.  Amortized O(1): a rebuild of size > 8n is preceded by > 4n
    // pushes since the previous one (a rebuild leaves <= one valid entry
    // per selectable vertex plus equal-key duplicates <= heap growth since
    // then).
    if (cand.size() > 8u * static_cast<std::size_t>(n_hint)) {
        sparse_lu_ordering_heap_compact_(cand, degree, alive, is_element);
    }
}

// ===========================================================================
// F2-b Phase 1 — sparse_lu_amd_ordering_v2_: ADD-96-faithful re-implementation
// of the AMD orderer (Amestoy/Davis/Duff, SIAM J. Matrix Anal. Appl. 17(4),
// 1996).  Detail name; NOT wired to any public entry until the Phase 3
// switch-over ruling.  Algorithmic content (paper sections in brackets):
//
//   - quotient graph on ONE flat Index work array `iw` [SS5]: variable lists
//     are [elements..., variables...] slices (pe/len/elen bookkeeping),
//     element lists are member-variable slices carved from the pivot's Lp;
//     new Lp goes in place of Ap when |Ep| == 0, otherwise into elbow room
//     at the tail, with DETERMINISTIC garbage collection (ascending-index
//     compaction) when the elbow is exhausted [SS5, MA27 storage].
//   - degree bucket lists dhead/dnext/dprev with a monotone-forward minimum-
//     degree pointer that retreats on insertion.  TIE-BREAK CONTRACT
//     (Phase 1 ruling 2026-08-04, superseding the design v1.1 SS4 min-index
//     pivot tie): pivot = minimum approximate degree; among equal degrees
//     the BUCKET HEAD is taken (head-insert / head-extract discipline, the
//     plain reading of ADD-96 SS5).  Supervariable principal selection,
//     merge order and every OTHER tie in this function remain smallest-index.
//     Determinism (D-2) is carried entirely by the following total order of
//     bucket operations -- no address values, no randomness:
//       (i)   initial fill scans vertices by DESCENDING index, inserting at
//             the head, so every initial bucket lists its members in
//             ascending index order;
//       (ii)  removal (unlink) never reorders the survivors;
//       (iii) every later insertion is a head insertion, and insertions
//             happen in a deterministic sequence: the degree-update pass
//             re-inserts the surviving members of Lp in Lp storage order
//             (which is itself deterministic, [V2-N2]);
//       (iv)  garbage collection moves list storage only and never touches
//             the bucket structure.
//     Two consecutive runs on the same input are verified byte-identical by
//     the unit test.
//   - approximate external degree = min of the ADD-96 THREE upper bounds,
//     eq (4), with |Le \ Lp| computed by the w(e) scan of Algorithm 2 and
//     eq (5).  All degrees/masses count ORIGINAL variables (supervariable
//     masses nv), as in the paper.
//   - mass elimination [SS3.2] and supervariable detection restricted to
//     i in Lp, via the paper's hash  Hash(i) = (sum(A_i) + sum(E_i))
//     mod (n-1) + 1  [SS5], exact list comparison inside a hash bucket,
//     pairs compared and merged in ASCENDING INDEX order (principal =
//     smallest index -- S-5; NOT affected by the pivot tie-break ruling).
//   - element absorption: natural (e in Ep) and aggressive (|Le \ Lp| == 0,
//     i.e. w(e) == 0) [SS5] -- the same rule as the F2-a orderer.
//   - NO dense-row special-casing (ruling H-4).  No randomness (D-2).
//
// Deviation ledger (documented, not silent):
//   [V2-N1] bound 1 of eq (4) is taken literally as n - k with k = the pivot
//           mass eliminated THROUGH p (the paper's "size of the active
//           submatrix" after step k).
//   [V2-N2] the rewritten element list of a variable i in Lp is ordered
//           [p, surviving old elements in stored order]; Lp itself collects
//           Ap first (stored order), then the members of each e in Ep
//           (stored order), first occurrence only.  The paper fixes no
//           order; this one is deterministic and documented.
//   [V2-N3] when the elbow room is exhausted even after a garbage
//           collection, iw grows by 3/2 (deterministic amount).  The paper
//           assumes "elbow room of size n is sufficient in practice".
// ===========================================================================
template <class Index>
std::vector<Index> sparse_lu_amd_ordering_v2_(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_amd_ordering_v2_: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> perm(un);
    if (n <= Index(1)) {
        for (Index i = Index(0); i < n; ++i) {
            perm[static_cast<std::size_t>(i)] = i;
        }
        return perm;
    }

    // ---- A + A^T pattern graph (S-1), sorted unique neighbour lists.
    std::vector<Index> adj_ptr, adj_ind;
    sparse_lu_build_symmetric_pattern_graph(n, col_ptr, row_ind, adj_ptr, adj_ind);
    const std::size_t adj_total = adj_ind.size();

    // ---- flat storage: iw holds every list; pe[i] start, len[i] entries.
    // For a live VARIABLE i the slice is [elements (elen[i])..., variables].
    // For an ELEMENT e the slice is its member variables (elen unused).
    std::vector<Index> iw(adj_total + un + 1u);
    std::vector<Index> pe(un), len(un), elen(un, Index(0));
    for (std::size_t i = 0; i < un; ++i) {
        pe[i]  = adj_ptr[i];
        len[i] = adj_ptr[i + 1u] - adj_ptr[i];
    }
    for (std::size_t k = 0; k < adj_total; ++k) iw[k] = adj_ind[k];
    std::size_t pfree = adj_total;             // first free slot (elbow room)

    std::vector<Index> nv(un, Index(1));       // supervariable mass (0 = non-principal)
    std::vector<Index> degree(un);             // approximate external degree (mass)
    std::vector<Index> esize(un, Index(0));    // element mass |Le| at formation
    std::vector<char>  is_elem(un, 0);
    std::vector<char>  dead(un, 0);            // absorbed element / non-principal var
    // Emission chains (same intrusive idiom as the F2-a orderer): the chain
    // of principal i enumerates the original variables it represents, in
    // merge-concatenation order.
    std::vector<Index> mem_head(un), mem_next(un, Index(-1)), mem_tail(un);
    for (std::size_t i = 0; i < un; ++i) {
        mem_head[i] = static_cast<Index>(i);
        mem_tail[i] = static_cast<Index>(i);
    }

    // ---- degree buckets (head-insert / head-extract discipline, see the
    // tie-break contract in the function header).
    std::vector<Index> dhead(un, Index(-1)), dnext(un, Index(-1)), dprev(un, Index(-1));
    std::vector<char>  in_bucket(un, 0);
    for (std::size_t s = un; s-- > 0u; ) {     // (i) descending index -> ascending initial lists
        const Index i = static_cast<Index>(s);
        degree[s] = Index(0);
        for (Index k = pe[s]; k < pe[s] + len[s]; ++k) {
            degree[s] += nv[static_cast<std::size_t>(iw[static_cast<std::size_t>(k)])];
        }
        const std::size_t d = static_cast<std::size_t>(degree[s]);
        dnext[s] = dhead[d];
        dprev[s] = Index(-1);
        if (dhead[d] >= Index(0)) dprev[static_cast<std::size_t>(dhead[d])] = i;
        dhead[d] = i;
        in_bucket[s] = 1;
    }

    // ---- per-round scratch
    std::vector<Index> w(un, Index(-1));       // Algorithm 2: |Le \ Lp| masses
    std::vector<char>  in_lp(un, 0);           // marks Lp members and p
    std::vector<Index> touched_w;              // elements with w set this round
    std::vector<Index> hhead(un, Index(-1)), hnext(un, Index(-1));
    std::vector<Index> touched_h;              // hash buckets used this round
    std::vector<char>  cmp_mark(un, 0);
    std::vector<Index> bucket_members;         // scratch for one hash bucket
    std::vector<Index> rebuf;                  // list-rewrite scratch (SPC-P1
    rebuf.reserve(un);                         // lesson: never let an in-place
                                               // write cursor pass unread slots)

    std::vector<Index> order;
    order.reserve(un);
    Index eliminated = Index(0);
    std::size_t mindeg = 0;

    // O(1) unlink from the degree bucket structure.
    // (local lambdas: C++11, no captures beyond references)
    struct bucket_ops {
        std::vector<Index>& dhead; std::vector<Index>& dnext;
        std::vector<Index>& dprev; std::vector<char>& in_bucket;
        const std::vector<Index>& degree;
        void remove(Index i) {
            const std::size_t ui = static_cast<std::size_t>(i);
            if (!in_bucket[ui]) return;
            const Index pv = dprev[ui], nx = dnext[ui];
            if (pv >= Index(0)) dnext[static_cast<std::size_t>(pv)] = nx;
            else dhead[static_cast<std::size_t>(degree[ui])] = nx;
            if (nx >= Index(0)) dprev[static_cast<std::size_t>(nx)] = pv;
            dprev[ui] = dnext[ui] = Index(-1);
            in_bucket[ui] = 0;
        }
        // (iii) O(1) HEAD insertion into bucket degree[i] (ruling: bucket-head
        // pivot tie-break; the caller sequence is the deterministic source of
        // the within-bucket order)
        void insert_head(Index i) {
            const std::size_t ui = static_cast<std::size_t>(i);
            const std::size_t d = static_cast<std::size_t>(degree[ui]);
            dprev[ui] = Index(-1);
            dnext[ui] = dhead[d];
            if (dhead[d] >= Index(0)) dprev[static_cast<std::size_t>(dhead[d])] = i;
            dhead[d] = i;
            in_bucket[ui] = 1;
        }
    } buckets = { dhead, dnext, dprev, in_bucket, degree };

    while (eliminated < n) {
        // ---- pivot: head of the first non-empty bucket = minimum degree,
        // smallest index (bucket lists are ascending).  The pointer only
        // moves forward here; insert sites retreat it.
        while (mindeg < un && dhead[mindeg] < Index(0)) ++mindeg;
        if (mindeg >= un) {
            vcp::throw_error<vcp::state_error>(
                "sparse_lu_amd_ordering_v2_: degree lists exhausted early");
        }
        const Index p = dhead[mindeg];
        const std::size_t up = static_cast<std::size_t>(p);
        buckets.remove(p);

        // ---- build Lp = (Ap  U  union of Le over e in Ep) minus p's own
        // supervariable, deduplicated via in_lp marks.  Storage: in place of
        // Ap when |Ep| == 0 (Lp is then a subset of Ap), else elbow room.
        in_lp[up] = 1;
        std::size_t lp_start;
        if (elen[up] == Index(0)) {
            lp_start = static_cast<std::size_t>(pe[up]);
        } else {
            // worst-case need for the deduplicated union
            std::size_t need = static_cast<std::size_t>(len[up] - elen[up]);
            for (Index k = pe[up]; k < pe[up] + elen[up]; ++k) {
                const std::size_t ue = static_cast<std::size_t>(iw[static_cast<std::size_t>(k)]);
                if (!dead[ue]) need += static_cast<std::size_t>(len[ue]);
            }
            if (pfree + need > iw.size()) {
                // deterministic garbage collection.  The slices MUST be
                // compacted in ascending STORAGE (pe) order -- not node-id
                // order: after a relocation the two orders diverge, and an
                // id-ordered forward copy lets the write cursor overrun
                // still-unread slices (in-bounds corruption).  Sorting by pe
                // is deterministic (live slices are disjoint, so keys are
                // unique).  Positions change, contents do not.
                std::vector<std::pair<Index, Index> > live;   // (pe, node)
                live.reserve(un);
                for (std::size_t v = 0; v < un; ++v) {
                    const bool live_var  = !is_elem[v] && !dead[v] && nv[v] > Index(0)
                                           && v != up;
                    const bool live_self = (v == up);
                    const bool live_elem = is_elem[v] && !dead[v];
                    if (!(live_var || live_elem || live_self)) { len[v] = Index(0); continue; }
                    live.push_back(std::make_pair(pe[v], static_cast<Index>(v)));
                }
                std::sort(live.begin(), live.end());
                std::size_t cursor = 0;
                for (std::size_t s = 0; s < live.size(); ++s) {
                    const std::size_t v = static_cast<std::size_t>(live[s].second);
                    const std::size_t b = static_cast<std::size_t>(pe[v]);
                    const std::size_t l = static_cast<std::size_t>(len[v]);
                    for (std::size_t t = 0; t < l; ++t) iw[cursor + t] = iw[b + t];
                    pe[v] = static_cast<Index>(cursor);
                    cursor += l;
                }
                pfree = cursor;
                if (pfree + need > iw.size()) {
                    iw.resize(pfree + need + (iw.size() >> 1));   // [V2-N3]
                }
            }
            lp_start = pfree;
        }

        std::size_t lp_len = 0;
        Index degp = Index(0);                 // |Lp| in original-variable mass
        {
            // Ap part (stored order) -- read BEFORE writing when in place:
            // the write cursor lp_start + lp_len never passes the read
            // cursor (elen == 0 case writes over its own prefix).
            const std::size_t ab = static_cast<std::size_t>(pe[up]) +
                                   static_cast<std::size_t>(elen[up]);
            const std::size_t al = static_cast<std::size_t>(len[up] - elen[up]);
            for (std::size_t t = 0; t < al; ++t) {
                const Index v = iw[ab + t];
                const std::size_t uv = static_cast<std::size_t>(v);
                // is_elem: a stale variable entry whose vertex has since been
                // eliminated must NOT enter Lp -- treating an element id as a
                // variable would corrupt its member list in the update pass.
                if (is_elem[uv] || dead[uv] || nv[uv] <= Index(0) || in_lp[uv]) continue;
                in_lp[uv] = 1;
                iw[lp_start + lp_len++] = v;
                degp += nv[uv];
            }
            // members of each element of Ep (stored order), first occurrence
            // only; the element itself is absorbed naturally.
            for (Index k = pe[up]; k < pe[up] + elen[up]; ++k) {
                const Index e = iw[static_cast<std::size_t>(k)];
                const std::size_t ue = static_cast<std::size_t>(e);
                if (dead[ue]) continue;
                for (Index t = pe[ue]; t < pe[ue] + len[ue]; ++t) {
                    const Index v = iw[static_cast<std::size_t>(t)];
                    const std::size_t uv = static_cast<std::size_t>(v);
                    if (is_elem[uv] || dead[uv] || nv[uv] <= Index(0) || in_lp[uv]) continue;
                    in_lp[uv] = 1;
                    iw[lp_start + lp_len++] = v;
                    degp += nv[uv];
                }
                dead[ue] = 1;                  // natural absorption
                len[ue] = Index(0);
            }
        }
        if (elen[up] != Index(0)) pfree = lp_start + lp_len;

        // every Lp member leaves the degree lists for this round
        for (std::size_t t = 0; t < lp_len; ++t) buckets.remove(iw[lp_start + t]);

        // ---- Algorithm 2: w(e) = |Le \ Lp| (mass) for every element on the
        // element list of some i in Lp.  w starts < 0 (eq 5 sentinel).
        for (std::size_t t = 0; t < lp_len; ++t) {
            const std::size_t ui = static_cast<std::size_t>(iw[lp_start + t]);
            for (Index k = pe[ui]; k < pe[ui] + elen[ui]; ++k) {
                const Index e = iw[static_cast<std::size_t>(k)];
                const std::size_t ue = static_cast<std::size_t>(e);
                if (dead[ue]) continue;
                if (w[ue] < Index(0)) {
                    w[ue] = esize[ue];
                    touched_w.push_back(e);
                }
                w[ue] -= nv[ui];
            }
        }

        // ---- degree update pass over i in Lp (stored order).  Compresses
        // both list halves in place, performs aggressive absorption, applies
        // the eq (4) three-bound minimum, and hashes i for supervariable
        // detection.
        const Index k_after = eliminated + nv[up];       // [V2-N1]
        touched_h.clear();
        for (std::size_t t = 0; t < lp_len; ++t) {
            const Index i = iw[lp_start + t];
            const std::size_t ui = static_cast<std::size_t>(i);

            // The rewritten list is assembled in `rebuf` and placed back:
            // in the same slice when it fits (the usual case -- either p was
            // in A_i or an element of E_i died, so the list shrank), else
            // RELOCATED to the elbow room (stale entries can break the
            // usual-case size argument, so fitting is checked, not assumed).
            rebuf.clear();
            // element half: p first, then surviving old elements [V2-N2]
            rebuf.push_back(p);
            Index esum = Index(0);
            unsigned long long hsum = static_cast<unsigned long long>(p);
            for (Index k = pe[ui]; k < pe[ui] + elen[ui]; ++k) {
                const Index e = iw[static_cast<std::size_t>(k)];
                const std::size_t ue = static_cast<std::size_t>(e);
                if (dead[ue]) continue;
                if (w[ue] == Index(0)) {       // aggressive absorption
                    dead[ue] = 1;
                    len[ue] = Index(0);
                    continue;
                }
                rebuf.push_back(e);
                esum += (w[ue] >= Index(0)) ? w[ue] : esize[ue];   // eq (5)
                hsum += static_cast<unsigned long long>(e);
            }
            const Index new_elen = static_cast<Index>(rebuf.size());
            // variable half: drop Lp members (now represented by p), p
            // itself, dead and non-principal entries
            Index asum = Index(0);
            for (Index k = pe[ui] + elen[ui]; k < pe[ui] + len[ui]; ++k) {
                const Index v = iw[static_cast<std::size_t>(k)];
                const std::size_t uv = static_cast<std::size_t>(v);
                if (is_elem[uv] || in_lp[uv] || dead[uv] || nv[uv] <= Index(0)) continue;
                rebuf.push_back(v);
                asum += nv[uv];
                hsum += static_cast<unsigned long long>(v);
            }
            std::size_t base = static_cast<std::size_t>(pe[ui]);
            if (rebuf.size() > static_cast<std::size_t>(len[ui])) {
                // does not fit in place: relocate to the elbow room.  No
                // garbage collection here (the Lp slice is not yet owned by
                // the not-yet-finalized element p and must not move);
                // deterministic resize instead when the elbow is short.
                if (pfree + rebuf.size() > iw.size()) {
                    iw.resize(pfree + rebuf.size() + (iw.size() >> 1));   // [V2-N3]
                }
                base = pfree;
                pe[ui] = static_cast<Index>(base);
                pfree += rebuf.size();
            }
            for (std::size_t t2 = 0; t2 < rebuf.size(); ++t2) {
                iw[base + t2] = rebuf[t2];
            }
            elen[ui] = new_elen;
            len[ui]  = static_cast<Index>(rebuf.size());

            // eq (4): min of the three upper bounds (all masses)
            const Index lp_minus_i = degp - nv[ui];
            Index d_new = n - k_after;                        // bound 1
            const Index b2 = degree[ui] + lp_minus_i;         // bound 2
            if (b2 < d_new) d_new = b2;
            const Index b3 = asum + lp_minus_i + esum;        // bound 3
            if (b3 < d_new) d_new = b3;
            degree[ui] = d_new;

            // paper hash, buckets chained in first-touch order
            const std::size_t h =
                static_cast<std::size_t>(hsum % static_cast<unsigned long long>(un - 1u)) + 1u;
            const std::size_t hslot = h - 1u;                 // store in [0, n-1)
            if (hhead[hslot] < Index(0)) touched_h.push_back(static_cast<Index>(hslot));
            hnext[ui] = hhead[hslot];
            hhead[hslot] = i;
        }

        // ---- supervariable detection inside each used hash bucket: compare
        // pairs in ascending index order; the smaller index is principal.
        for (std::size_t hb = 0; hb < touched_h.size(); ++hb) {
            const std::size_t hslot = static_cast<std::size_t>(touched_h[hb]);
            bucket_members.clear();
            for (Index c = hhead[hslot]; c >= Index(0);
                 c = hnext[static_cast<std::size_t>(c)]) {
                bucket_members.push_back(c);
            }
            hhead[hslot] = Index(-1);
            if (bucket_members.size() < 2u) continue;
            std::sort(bucket_members.begin(), bucket_members.end());
            for (std::size_t ia = 0; ia < bucket_members.size(); ++ia) {
                const Index a = bucket_members[ia];
                const std::size_t ua = static_cast<std::size_t>(a);
                if (nv[ua] <= Index(0)) continue;             // already merged away
                // mark a's list
                for (Index k = pe[ua]; k < pe[ua] + len[ua]; ++k) {
                    cmp_mark[static_cast<std::size_t>(iw[static_cast<std::size_t>(k)])] = 1;
                }
                for (std::size_t ib = ia + 1u; ib < bucket_members.size(); ++ib) {
                    const Index b = bucket_members[ib];
                    const std::size_t ub = static_cast<std::size_t>(b);
                    if (nv[ub] <= Index(0)) continue;
                    if (len[ub] != len[ua] || elen[ub] != elen[ua]) continue;
                    bool same = true;
                    for (Index k = pe[ub]; k < pe[ub] + len[ub]; ++k) {
                        if (!cmp_mark[static_cast<std::size_t>(iw[static_cast<std::size_t>(k)])]) {
                            same = false;
                            break;
                        }
                    }
                    if (!same) continue;
                    // merge b into a (principal = smaller index, S-5)
                    nv[ua] += nv[ub];
                    degree[ua] -= nv[ub];                     // Algorithm 1: d_i -= |j|
                    nv[ub] = Index(0);
                    dead[ub] = 1;
                    len[ub] = Index(0);
                    elen[ub] = Index(0);
                    mem_next[static_cast<std::size_t>(mem_tail[ua])] = mem_head[ub];
                    mem_tail[ua] = mem_tail[ub];
                }
                for (Index k = pe[ua]; k < pe[ua] + len[ua]; ++k) {
                    cmp_mark[static_cast<std::size_t>(iw[static_cast<std::size_t>(k)])] = 0;
                }
            }
        }

        // ---- surviving Lp members re-enter the degree lists at the bucket
        // head, in Lp storage order (determinism source (iii)); mindeg
        // retreats.
        for (std::size_t t = 0; t < lp_len; ++t) {
            const Index i = iw[lp_start + t];
            const std::size_t ui = static_cast<std::size_t>(i);
            if (nv[ui] <= Index(0) || dead[ui]) continue;
            if (degree[ui] < Index(0)) degree[ui] = Index(0);  // defensive clamp
            buckets.insert_head(i);
            if (static_cast<std::size_t>(degree[ui]) < mindeg) {
                mindeg = static_cast<std::size_t>(degree[ui]);
            }
        }

        // ---- finalize element p; emit its represented original variables.
        is_elem[up] = 1;
        pe[up] = static_cast<Index>(lp_start);
        len[up] = static_cast<Index>(lp_len);
        elen[up] = Index(0);
        esize[up] = degp;
        for (Index s = mem_head[up]; s >= Index(0);
             s = mem_next[static_cast<std::size_t>(s)]) {
            order.push_back(s);
        }
        eliminated += nv[up];

        // ---- reset the per-round marks (touch lists only)
        in_lp[up] = 0;
        for (std::size_t t = 0; t < lp_len; ++t) {
            in_lp[static_cast<std::size_t>(iw[lp_start + t])] = 0;
        }
        for (std::size_t t = 0; t < touched_w.size(); ++t) {
            w[static_cast<std::size_t>(touched_w[t])] = Index(-1);
        }
        touched_w.clear();
    }

    if (order.size() != un) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_amd_ordering_v2_: elimination did not cover all vertices");
    }
    for (std::size_t i = 0; i < un; ++i) {
        perm[i] = order[i];
    }
    return perm;
}

// ---------------------------------------------------------------------------
// F2-b Phase 3 (gate-2 continuation ruling 2026-08-04): the PUBLIC AMD entry.
// Forwards to the ADD-96 v2 implementation above; every amd consumer (lu /
// ldl / chol / ndml leaves / fsai) resolves through this name.  The pre-F2-b
// body (sparse_lu_amd_ordering_legacy_, zero callers) was physically removed
// in F2-b Phase 5 per the close-out ruling (H-3-iv); roll-back = git history
// (last present at commit aa5a108).
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index> sparse_lu_amd_ordering(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    return sparse_lu_amd_ordering_v2_(n, col_ptr, row_ind);
}

// ===========================================================================
// Ordering Track O3 — COLAMD (Column Approximate Minimum Degree) ordering.
//
// COLAMD (Davis / Gilbert / Larimore / Ng) computes a fill-reducing COLUMN
// ordering for UNSYMMETRIC patterns: a minimum-degree elimination on the
// column-intersection graph, i.e. the sparsity pattern of A^T A.  The key fact
// that makes this safe and cheap is that A^T A is NOT formed explicitly (it can
// be dense); instead the rows of A are treated as the INITIAL elements of a
// quotient graph (each row is a hyperedge joining all columns nonzero in it),
// exactly the structure the design doc §634 prescribes ("A^T A の tree を陽に
// 形成せず求める").  Two columns are adjacent in A^T A iff some row contains
// both — which is precisely "they share an element".
//
// Contrast with O1/O2: RCM/AMD work on the A + A^T pattern (symmetrised); for a
// strongly unsymmetric A that symmetrisation over-densifies and its order is not
// optimal for the original unsymmetric A.  COLAMD orders A^T A directly, so it
// can beat AMD on advection-dominated / circuit-like unsymmetric patterns.
//
// Quotient-graph state.  Columns are VARIABLES (ids 0..n-1); rows are the
// initial ELEMENTS (ids 0..n-1).  Eliminating a pivot column creates a NEW
// element (id >= n, appended).  For a (principal) column j: E_col[j] is its
// sorted set of adjacent element ids.  For an element e: V_elem[e] is its sorted
// set of member (live principal) columns; j in V_elem[e] iff e in E_col[j].
// The external column degree is the mass of the union of V_elem over E_col[j],
// minus j itself — recomputed exactly for the pivot's neighbours (local-exact /
// global-approximate, the standard min-degree discipline).
//
// Safety invariants (roadmap O3 §2.3 = S-1..S-6):
//   S-1 Integer-only / pattern-only.  Reads col_ptr / row_ind only; never values.
//   S-5 Deterministic.  Pivot = minimum degree, tie-break smallest index;
//       supervariable principal = smallest index; sorted lists; no randomness.
//   S-6 The result is a COLUMN permutation Q installed into the col_perm slot;
//       the answer is unchanged (cost/structure only).
//
// n <= 1 returns identity.  Disconnected patterns, empty rows/columns, and
// zero-degree (fully zero) columns are handled.  A coverage check throws if the
// elimination did not cover every column exactly once.
// ===========================================================================
template <class Index>
std::vector<Index> sparse_lu_colamd_ordering(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_colamd_ordering: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> perm(un);
    if (n <= Index(1)) {
        for (Index i = Index(0); i < n; ++i) {
            perm[static_cast<std::size_t>(i)] = i;
        }
        return perm;
    }

    // ---- build the bipartite column<->element (row) incidence WITHOUT forming
    // A^T A (S-1: pattern-only).  E_col[j] = rows containing column j (the CSC
    // column pattern, de-duplicated).  V_elem[i] = columns nonzero in row i (the
    // transpose pattern).  Both are sorted ascending and unique; they are exact
    // transposes of each other.
    std::vector<std::vector<Index> > E_col(un);  // column -> adjacent elements
    std::vector<std::vector<Index> > V_elem(un); // element(row) -> member columns
    std::vector<char> elem_alive(un, 0);

    for (Index c = Index(0); c < n; ++c) {
        const Index kb = col_ptr[static_cast<std::size_t>(c)];
        const Index ke = col_ptr[static_cast<std::size_t>(c) + 1u];
        std::vector<Index>& Ec = E_col[static_cast<std::size_t>(c)];
        for (Index k = kb; k < ke; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r < Index(0) || r >= n) continue;  // defensive; convert validated
            Ec.push_back(r);
            V_elem[static_cast<std::size_t>(r)].push_back(c);
            elem_alive[static_cast<std::size_t>(r)] = 1;
        }
        std::sort(Ec.begin(), Ec.end());
        Ec.erase(std::unique(Ec.begin(), Ec.end()), Ec.end());
    }
    for (std::size_t i = 0; i < un; ++i) {
        std::sort(V_elem[i].begin(), V_elem[i].end());
        V_elem[i].erase(std::unique(V_elem[i].begin(), V_elem[i].end()),
                        V_elem[i].end());
    }

    std::vector<Index> nv(un, Index(1));      // supervariable mass
    std::vector<std::vector<Index> > members(un);
    std::vector<Index> col_deg(un, Index(0)); // external column degree (mass)
    std::vector<char>  col_alive(un, 1);
    for (std::size_t i = 0; i < un; ++i) members[i].push_back(static_cast<Index>(i));

    // Stamp scratch for distinct-union counting over columns.
    std::vector<Index> mark(un, Index(-1));
    Index stamp = Index(0);

    // Exact external degree of a live principal column j: mass of distinct live
    // columns (other than j) reachable through j's elements.
    // (Defined as a lambda-free helper via a small loop for C++11 clarity.)
    // Initial degrees.
    for (Index j = Index(0); j < n; ++j) {
        ++stamp;
        const std::size_t uj = static_cast<std::size_t>(j);
        Index cnt = Index(0);
        const std::vector<Index>& Ej = E_col[uj];
        for (std::size_t s = 0; s < Ej.size(); ++s) {
            const std::vector<Index>& Ve = V_elem[static_cast<std::size_t>(Ej[s])];
            for (std::size_t t = 0; t < Ve.size(); ++t) {
                const Index c = Ve[t];
                const std::size_t uc = static_cast<std::size_t>(c);
                if (c == j) continue;
                if (mark[uc] == stamp) continue;
                mark[uc] = stamp;
                cnt += nv[uc];
            }
        }
        col_deg[uj] = cnt;
    }

    // SLU-OQ1: ordered candidate set (design D-1) -- since ORD-F3 P4 a
    // lazy-deletion binary heap (the ORD-F2 A1 mechanism rolled out
    // horizontally): at least one (col_deg[j], j) entry per selectable
    // column, minimum on top via pq_greater_, stale entries discarded at
    // pop time by the validity check
    //     col_alive[v] && key == col_deg[v].
    // No duplicate-scan stamp is needed, same argument as amd: exactly one
    // entry is adopted per pivot and the adopted column immediately goes
    // col_alive = 0, which stales every leftover entry of it.  Named candq
    // (not cand) because the aggressive-absorption block below has a local
    // scratch vector `cand` of its own (candidate ELEMENTS -- untouched).
    // The pop-time validity check reuses the shared update helper, whose
    // is_element argument is served by an all-zero array (columns never
    // become elements in colamd's column view; new elements live in V_elem
    // ids >= n and never enter candq).
    std::vector<std::pair<Index, Index> > candq;
    candq.reserve(un);
    for (Index j = Index(0); j < n; ++j) {
        candq.push_back(std::make_pair(col_deg[static_cast<std::size_t>(j)], j));
    }
    std::make_heap(candq.begin(), candq.end(), sparse_lu_ordering_pq_greater_());
    const std::vector<char> col_is_element(un, 0);  // always 0, see above

    std::vector<Index> order;
    order.reserve(un);
    Index next_elem = n;        // id of the next NEW element to create
    Index eliminated = Index(0);

    while (eliminated < n) {
        // ---- pivot: minimum external degree, tie-break smallest index (S-5).
        // SLU-OQ1 / ORD-F3 P4: the first VALID popped entry
        // (col_alive && key == col_deg[v]) is the lexicographic min of
        // (degree, index) among selectable columns -- identical to
        // *cand.begin() of the former std::set (design F1), so the emitted
        // permutation is bit-identical.  Adopting the popped entry replaces
        // the former cand.erase(cand.begin()) (T-C1).
        Index p = Index(-1);
        while (!candq.empty()) {
            const std::pair<Index, Index> top = candq.front();
            std::pop_heap(candq.begin(), candq.end(),
                          sparse_lu_ordering_pq_greater_());
            candq.pop_back();
            const std::size_t uv = static_cast<std::size_t>(top.second);
            if (col_alive[uv] && top.first == col_deg[uv]) {
                p = top.second;  // p goes col_alive = 0 below (T-C1)
                break;
            }
        }
        if (p < Index(0)) break;  // defensive (was: cand.empty() fallthrough)
        const std::size_t up = static_cast<std::size_t>(p);

        // ---- form Lp = union of V_elem[e] over e in E_col[p], live columns,
        // excluding p.  These columns become members of the new pivot element.
        ++stamp;
        std::vector<Index> lp;
        const std::vector<Index>& Ep = E_col[up];
        for (std::size_t s = 0; s < Ep.size(); ++s) {
            const std::vector<Index>& Ve = V_elem[static_cast<std::size_t>(Ep[s])];
            for (std::size_t t = 0; t < Ve.size(); ++t) {
                const Index c = Ve[t];
                const std::size_t uc = static_cast<std::size_t>(c);
                if (c == p) continue;
                if (!col_alive[uc]) continue;
                if (mark[uc] == stamp) continue;
                mark[uc] = stamp;
                lp.push_back(c);
            }
        }
        std::sort(lp.begin(), lp.end());

        // ---- create the new element ep absorbing every element in E_col[p].
        const Index ep = next_elem++;
        const std::vector<Index> absorbed = Ep;  // sorted; absorbed into ep
        V_elem.push_back(lp);                     // V_elem[ep] = Lp
        elem_alive.push_back(1);

        // Eliminate p: emit its represented original columns, retire it.
        col_alive[up] = 0;
        for (std::size_t s = 0; s < members[up].size(); ++s) {
            order.push_back(members[up][s]);
        }
        eliminated += nv[up];

        // ---- rewrite each column j in Lp: drop absorbed elements, attach ep.
        // Every column referencing an absorbed element is in Lp (its members are
        // a subset of Lp), so no dangling references remain.
        for (std::size_t t = 0; t < lp.size(); ++t) {
            const std::size_t uj = static_cast<std::size_t>(lp[t]);
            sparse_lu_amd_sorted_remove_set(E_col[uj], absorbed);
            sparse_lu_amd_sorted_insert(E_col[uj], ep);
        }
        // Retire absorbed elements.
        for (std::size_t s = 0; s < absorbed.size(); ++s) {
            const std::size_t ue = static_cast<std::size_t>(absorbed[s]);
            elem_alive[ue] = 0;
            V_elem[ue].clear();
        }

        // ---- aggressive element absorption: any element e' (e' != ep) all of
        // whose live members lie in Lp is now redundant (its clique is a subset
        // of ep's), so absorb it.  Mark Lp membership with a stamp, then for each
        // distinct candidate element referenced by an Lp column test set-subset.
        {
            const Index lp_stamp = ++stamp;
            for (std::size_t t = 0; t < lp.size(); ++t) {
                mark[static_cast<std::size_t>(lp[t])] = lp_stamp;  // Lp membership
            }
            std::vector<Index> cand;  // distinct candidate elements (e != ep)
            for (std::size_t t = 0; t < lp.size(); ++t) {
                const std::vector<Index>& Ej = E_col[static_cast<std::size_t>(lp[t])];
                for (std::size_t s = 0; s < Ej.size(); ++s) {
                    if (Ej[s] != ep) cand.push_back(Ej[s]);
                }
            }
            std::sort(cand.begin(), cand.end());
            cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
            std::vector<Index> absorb2;
            for (std::size_t s = 0; s < cand.size(); ++s) {
                const std::size_t ue = static_cast<std::size_t>(cand[s]);
                if (!elem_alive[ue] || V_elem[ue].empty()) continue;
                bool subset = true;
                for (std::size_t q = 0; q < V_elem[ue].size(); ++q) {
                    const std::size_t uc = static_cast<std::size_t>(V_elem[ue][q]);
                    if (!col_alive[uc]) continue;          // dead member: ignore
                    if (mark[uc] != lp_stamp) { subset = false; break; }
                }
                if (subset) absorb2.push_back(cand[s]);
            }
            if (!absorb2.empty()) {
                std::sort(absorb2.begin(), absorb2.end());
                absorb2.erase(std::unique(absorb2.begin(), absorb2.end()),
                              absorb2.end());
                for (std::size_t t = 0; t < lp.size(); ++t) {
                    sparse_lu_amd_sorted_remove_set(
                        E_col[static_cast<std::size_t>(lp[t])], absorb2);
                }
                for (std::size_t s = 0; s < absorb2.size(); ++s) {
                    const std::size_t ue = static_cast<std::size_t>(absorb2[s]);
                    elem_alive[ue] = 0;
                    V_elem[ue].clear();
                }
            }
        }

        // ---- supervariable detection within Lp: columns with identical element
        // lists are indistinguishable in A^T A; merge larger index into smaller
        // (S-5).  Bucket by a structural hash, then compare exactly within a run.
        if (lp.size() > 1u) {
            std::vector<std::pair<unsigned long long, Index> > hb;
            hb.reserve(lp.size());
            for (std::size_t t = 0; t < lp.size(); ++t) {
                const Index i = lp[t];
                const std::size_t ui = static_cast<std::size_t>(i);
                if (!col_alive[ui]) continue;
                unsigned long long h = 1469598103934665603ull;
                for (std::size_t s = 0; s < E_col[ui].size(); ++s) {
                    h = (h ^ static_cast<unsigned long long>(E_col[ui][s] + 1))
                        * 1099511628211ull;
                }
                hb.push_back(std::make_pair(h, i));
            }
            std::sort(hb.begin(), hb.end());
            for (std::size_t a = 0; a < hb.size(); ++a) {
                const std::size_t ui = static_cast<std::size_t>(hb[a].second);
                if (!col_alive[ui]) continue;
                for (std::size_t b = a + 1;
                     b < hb.size() && hb[b].first == hb[a].first; ++b) {
                    const Index j = hb[b].second;
                    const std::size_t uj = static_cast<std::size_t>(j);
                    if (!col_alive[uj]) continue;
                    if (!sparse_lu_amd_lists_equal(E_col[ui], E_col[uj])) continue;
                    // Merge j (larger index) into i (smaller index).
                    nv[ui] += nv[uj];
                    members[ui].insert(members[ui].end(),
                                       members[uj].begin(), members[uj].end());
                    members[uj].clear();
                    for (std::size_t s = 0; s < E_col[uj].size(); ++s) {
                        sparse_lu_amd_sorted_remove_one(
                            V_elem[static_cast<std::size_t>(E_col[uj][s])], j);
                    }
                    // T-C2: no candq retire needed (lazy heap, ORD-F3 P4)
                    // -- flipping col_alive[uj] below stales every stored
                    // entry of j.
                    col_alive[uj] = 0;
                    nv[uj] = Index(0);
                    E_col[uj].clear();
                }
            }
        }

        // ---- recompute exact external degree for surviving principals in Lp.
        for (std::size_t t = 0; t < lp.size(); ++t) {
            const std::size_t uj = static_cast<std::size_t>(lp[t]);
            if (!col_alive[uj]) continue;
            ++stamp;
            Index cnt = Index(0);
            const std::vector<Index>& Ej = E_col[uj];
            for (std::size_t s = 0; s < Ej.size(); ++s) {
                const std::vector<Index>& Ve =
                    V_elem[static_cast<std::size_t>(Ej[s])];
                for (std::size_t q = 0; q < Ve.size(); ++q) {
                    const Index c = Ve[q];
                    const std::size_t uc = static_cast<std::size_t>(c);
                    if (c == lp[t]) continue;
                    if (!col_alive[uc]) continue;
                    if (mark[uc] == stamp) continue;
                    mark[uc] = stamp;
                    cnt += nv[uc];
                }
            }
            sparse_lu_ordering_cand_update_degree_(candq, col_deg, lp[t], cnt,
                                                   col_alive, col_is_element,
                                                   n);  // W-C2
        }
    }

    // Coverage invariant: every column eliminated exactly once.
    if (order.size() != un) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_colamd_ordering: elimination did not cover all columns");
    }

    // Elimination order is the column permutation directly.  perm[new] = old.
    for (std::size_t i = 0; i < un; ++i) {
        perm[i] = order[i];
    }
    return perm;
}

// ---------------------------------------------------------------------------
// SLU-MF4  Nested dissection ordering (self-contained, no external library).
//
// Purpose (894 axis-1): produce large separators so the assembly tree carries
// wide dense fronts; this raises the average supernode/front width k, which is
// what makes the dense BLAS-3 work (gemm ~ 2 m n k) dominate the memory-bound
// contribution scatter (extend-add ~ m n) in the multifrontal numeric source.
//
// Method: recursive graph bisection on the pattern of A + A^T.  For each
// connected piece we build a rooted level structure from a (restricted)
// pseudo-peripheral start and take one BFS level as the vertex separator,
// recursing on the two halves and numbering the separator LAST.  Numbering
// separators last places them near the root of the elimination tree, where they
// become the wide dense fronts.
//
// Convention: returns a COLUMN permutation Q with col_perm[new] = old, exactly
// like the rcm / amd / colamd orderings.  Pattern-only (S-1), deterministic
// (S-5): fixed ascending tie-breaks, no randomness.  Leaf pieces below a small
// threshold are emitted in ascending index order.
// ---------------------------------------------------------------------------

// Restricted rooted level structure: BFS that only follows neighbours whose
// membership marker mark[w] == curtag.  Fills `levels` (each level ascending);
// resets `seen` for exactly the vertices it touched.  Returns the eccentricity
// (number of levels - 1) reached from `root` inside the marked subset.
template <class Index>
Index sparse_lu_nd_restricted_levels(
    Index root,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    const std::vector<Index>& mark,
    Index curtag,
    std::vector<char>& seen,
    std::vector<std::vector<Index> >& levels)
{
    levels.clear();
    std::vector<Index> touched;

    std::vector<Index> frontier;
    frontier.push_back(root);
    seen[static_cast<std::size_t>(root)] = 1;
    touched.push_back(root);

    while (!frontier.empty()) {
        levels.push_back(frontier);
        std::vector<Index> next;
        for (std::size_t i = 0; i < frontier.size(); ++i) {
            const Index u = frontier[i];
            const std::size_t b = static_cast<std::size_t>(adj_ptr[static_cast<std::size_t>(u)]);
            const std::size_t e = static_cast<std::size_t>(adj_ptr[static_cast<std::size_t>(u) + 1u]);
            for (std::size_t p = b; p < e; ++p) {
                const Index w = adj_ind[p];
                if (mark[static_cast<std::size_t>(w)] != curtag) continue;
                if (seen[static_cast<std::size_t>(w)]) continue;
                seen[static_cast<std::size_t>(w)] = 1;
                touched.push_back(w);
                next.push_back(w);
            }
        }
        std::sort(next.begin(), next.end());
        frontier.swap(next);
    }
    for (std::size_t i = 0; i < touched.size(); ++i) {
        seen[static_cast<std::size_t>(touched[i])] = 0;
    }
    return static_cast<Index>(levels.size()) - Index(1);
}

// Restricted pseudo-peripheral start (George-Liu), confined to the marked
// subset.  Bounded iterations; deterministic (min-degree in deepest level,
// tie-broken by smallest index).
template <class Index>
Index sparse_lu_nd_restricted_peripheral(
    Index seed,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    const std::vector<Index>& degree,
    const std::vector<Index>& mark,
    Index curtag,
    std::vector<char>& seen)
{
    Index v = seed;
    std::vector<std::vector<Index> > levels;
    Index ecc = sparse_lu_nd_restricted_levels(
        v, adj_ptr, adj_ind, mark, curtag, seen, levels);

    const int max_iter = 16;
    for (int it = 0; it < max_iter; ++it) {
        if (levels.empty()) break;
        const std::vector<Index>& last = levels.back();
        Index cand = last[0];
        Index cand_deg = degree[static_cast<std::size_t>(cand)];
        for (std::size_t i = 1; i < last.size(); ++i) {
            const Index w = last[i];
            const Index d = degree[static_cast<std::size_t>(w)];
            if (d < cand_deg) { cand = w; cand_deg = d; }
        }
        std::vector<std::vector<Index> > cand_levels;
        const Index cand_ecc = sparse_lu_nd_restricted_levels(
            cand, adj_ptr, adj_ind, mark, curtag, seen, cand_levels);
        if (cand_ecc > ecc) {
            v = cand;
            ecc = cand_ecc;
            levels.swap(cand_levels);
        } else {
            break;
        }
    }
    return v;
}

// Recursive bisection.  `sub` is the current vertex piece (any order).  Appends
// the elimination order of this piece to `order` (children first, separator
// last).  `mark`/`seen` are size-n scratch; `tag` is a strictly increasing tag
// source so each call owns a unique membership stamp.
template <class Index>
void sparse_lu_nd_recurse(
    const std::vector<Index>& sub,
    const std::vector<Index>& adj_ptr,
    const std::vector<Index>& adj_ind,
    const std::vector<Index>& degree,
    std::vector<Index>& mark,
    std::vector<char>& seen,
    Index& tag,
    std::size_t leaf,
    std::vector<Index>& order)
{
    const std::size_t ns = sub.size();
    if (ns == 0) return;
    if (ns <= leaf) {
        std::vector<Index> s = sub;
        std::sort(s.begin(), s.end());
        for (std::size_t i = 0; i < s.size(); ++i) order.push_back(s[i]);
        return;
    }

    // Stamp the current subset so the restricted BFS stays inside it.
    const Index curtag = ++tag;
    for (std::size_t i = 0; i < ns; ++i) {
        mark[static_cast<std::size_t>(sub[i])] = curtag;
    }

    // Start from the smallest-index vertex, then refine to pseudo-peripheral.
    Index start = sub[0];
    for (std::size_t i = 1; i < ns; ++i) if (sub[i] < start) start = sub[i];
    start = sparse_lu_nd_restricted_peripheral(
        start, adj_ptr, adj_ind, degree, mark, curtag, seen);

    std::vector<std::vector<Index> > levels;
    sparse_lu_nd_restricted_levels(
        start, adj_ptr, adj_ind, mark, curtag, seen, levels);

    // Flatten the component reachable from `start`.
    std::vector<Index> comp;
    comp.reserve(ns);
    for (std::size_t d = 0; d < levels.size(); ++d)
        for (std::size_t i = 0; i < levels[d].size(); ++i)
            comp.push_back(levels[d][i]);

    if (comp.size() < ns) {
        // Disconnected piece: split off this component, recurse on each part
        // independently (they share no edges, so ordering between them is free).
        for (std::size_t i = 0; i < comp.size(); ++i)
            seen[static_cast<std::size_t>(comp[i])] = 1;
        std::vector<Index> rest;
        rest.reserve(ns - comp.size());
        for (std::size_t i = 0; i < ns; ++i)
            if (!seen[static_cast<std::size_t>(sub[i])]) rest.push_back(sub[i]);
        for (std::size_t i = 0; i < comp.size(); ++i)
            seen[static_cast<std::size_t>(comp[i])] = 0;
        sparse_lu_nd_recurse(comp, adj_ptr, adj_ind, degree, mark, seen, tag, leaf, order);
        sparse_lu_nd_recurse(rest, adj_ptr, adj_ind, degree, mark, seen, tag, leaf, order);
        return;
    }

    // Need at least 3 levels for both halves to be non-empty around a separator.
    if (levels.size() < 3) {
        std::vector<Index> s = sub;
        std::sort(s.begin(), s.end());
        for (std::size_t i = 0; i < s.size(); ++i) order.push_back(s[i]);
        return;
    }

    // Choose the separator level d* (1 <= d* <= L-2) minimizing imbalance
    // |left - right|, tie-broken by smaller separator then smaller d*.
    const std::size_t L = levels.size();
    const Index total = static_cast<Index>(ns);
    Index best_d = Index(1);
    Index best_score = total + Index(1);
    Index best_sep = total + Index(1);
    Index prefix = Index(0);
    for (std::size_t d = 0; d + 1 < L; ++d) {
        if (d >= 1) {
            const Index sep = static_cast<Index>(levels[d].size());
            const Index left = prefix;                 // levels[0..d-1]
            const Index right = total - left - sep;     // levels[d+1..]
            Index score = (left > right) ? (left - right) : (right - left);
            if (score < best_score ||
                (score == best_score && sep < best_sep)) {
                best_score = score;
                best_sep = sep;
                best_d = static_cast<Index>(d);
            }
        }
        prefix += static_cast<Index>(levels[d].size());
    }

    const std::size_t ds = static_cast<std::size_t>(best_d);
    std::vector<Index> left_part, right_part, sep_part = levels[ds];
    for (std::size_t d = 0; d < ds; ++d)
        for (std::size_t i = 0; i < levels[d].size(); ++i)
            left_part.push_back(levels[d][i]);
    for (std::size_t d = ds + 1; d < L; ++d)
        for (std::size_t i = 0; i < levels[d].size(); ++i)
            right_part.push_back(levels[d][i]);

    // Number both halves first, then the separator (placed near the root).
    sparse_lu_nd_recurse(left_part, adj_ptr, adj_ind, degree, mark, seen, tag, leaf, order);
    sparse_lu_nd_recurse(right_part, adj_ptr, adj_ind, degree, mark, seen, tag, leaf, order);
    std::sort(sep_part.begin(), sep_part.end());
    for (std::size_t i = 0; i < sep_part.size(); ++i) order.push_back(sep_part[i]);
}

template <class Index>
std::vector<Index> sparse_lu_nested_dissection_ordering(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_nested_dissection_ordering: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> perm(un);
    if (n <= Index(1)) {
        for (Index i = Index(0); i < n; ++i) perm[static_cast<std::size_t>(i)] = i;
        return perm;
    }

    std::vector<Index> adj_ptr, adj_ind;
    sparse_lu_build_symmetric_pattern_graph(n, col_ptr, row_ind, adj_ptr, adj_ind);

    std::vector<Index> degree(un);
    for (std::size_t i = 0; i < un; ++i)
        degree[i] = adj_ptr[i + 1u] - adj_ptr[i];

    std::vector<Index> mark(un, Index(-1));
    std::vector<char> seen(un, 0);
    Index tag = Index(0);

    // Leaf threshold: below this, a piece is numbered in ascending index order
    // rather than bisected further.  Small enough that wide separators dominate,
    // large enough to avoid a long tail of tiny one-vertex separators.
    const std::size_t leaf = 8u;

    std::vector<Index> all(un);
    for (std::size_t i = 0; i < un; ++i) all[i] = static_cast<Index>(i);

    std::vector<Index> order;
    order.reserve(un);
    sparse_lu_nd_recurse(all, adj_ptr, adj_ind, degree, mark, seen, tag, leaf, order);

    if (order.size() != un) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_nested_dissection_ordering: order did not cover all vertices");
    }
    for (std::size_t i = 0; i < un; ++i) perm[i] = order[i];
    return perm;
}

#endif // VCP_TSPARSE_SPARSE_LU_ORDERING_IMPL_HPP
