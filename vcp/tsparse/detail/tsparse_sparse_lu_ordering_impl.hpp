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
// SLU-OQ1: ordered candidate set for O(log n) pivot selection.
//
// `cand` holds exactly one (degree[i], i) pair for every selectable vertex i
// (amd: alive && !is_element; colamd: col_alive).  Selection is *cand.begin()
// = the lexicographic minimum of (degree, index), which is IDENTICAL to the
// former ascending linear scan with strict-< update (first index attaining
// the running minimum) -- design F1, so the emitted permutation is
// bit-identical.  Every degree[] rewrite and every eligibility loss MUST go
// through these two helpers so the stored degree and the set key never
// diverge (Phase 0 inventory §2/§3 enumerates all call sites).
// ---------------------------------------------------------------------------
template <class Index>
inline void sparse_lu_ordering_cand_update_degree_(
    std::set<std::pair<Index, Index> >& cand,
    std::vector<Index>& degree,
    const Index i,
    const Index new_degree)
{
    cand.erase(std::make_pair(degree[static_cast<std::size_t>(i)], i));
    degree[static_cast<std::size_t>(i)] = new_degree;
    cand.insert(std::make_pair(new_degree, i));
}

template <class Index>
inline void sparse_lu_ordering_cand_retire_(
    std::set<std::pair<Index, Index> >& cand,
    const std::vector<Index>& degree,
    const Index i)
{
    cand.erase(std::make_pair(degree[static_cast<std::size_t>(i)], i));
}

// ---------------------------------------------------------------------------
// ORD-F2 Phase A1: lazy-deletion binary-heap variant of the update_degree
// helper above, used by the amd section only (the colamd section keeps the
// std::set form; the helper overloads on the container type).  Same total
// order as the std::set it replaces -- minimum (degree, index) entry on top
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

template <class Index>
std::vector<Index> sparse_lu_amd_ordering(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse LU Index must be signed");
    if (n < Index(0)) {
        vcp::throw_error<vcp::invalid_argument>(
            "sparse_lu_amd_ordering: negative n");
    }
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> perm(un);
    if (n <= Index(1)) {
        for (Index i = Index(0); i < n; ++i) {
            perm[static_cast<std::size_t>(i)] = i;
        }
        return perm;
    }

    // A + A^T pattern graph (S-1: pattern-only).  Neighbour lists are sorted
    // ascending and unique, which is exactly the invariant the quotient-graph
    // set operations below rely on.
    std::vector<Index> adj_ptr, adj_ind;
    sparse_lu_build_symmetric_pattern_graph(n, col_ptr, row_ind, adj_ptr, adj_ind);

    // Quotient-graph state.  For a (principal) variable i: A_var[i] is the set
    // of adjacent principal variables and E_elem[i] the set of adjacent
    // elements.  For an element e: V_elem[e] is the set of its member principal
    // variables and elem_size[e] its mass |Le| = sum of nv over its members.
    std::vector<std::vector<Index> > A_var(un);
    std::vector<std::vector<Index> > E_elem(un);
    std::vector<std::vector<Index> > V_elem(un);
    // ORD-F3 P3: the represented-original-variables list `members` (formerly
    // vector<vector<Index>>, one allocation per vertex) flattened to an
    // intrusive linked list over the vertex ids: the chain of principal i,
    // followed head -> next until -1, enumerates EXACTLY the sequence the
    // former members[i] vector held (init = singleton {i}; merge = append
    // j's whole chain at i's tail, the same concatenation order as the
    // former insert(end, begin, end)).
    std::vector<Index> mem_head(un);               // chain start (= i itself)
    std::vector<Index> mem_next(un, Index(-1));    // successor, -1 terminates
    std::vector<Index> mem_tail(un);               // last node of i's chain
    std::vector<Index> nv(un, Index(1));           // supervariable mass
    std::vector<Index> elem_size(un, Index(0));    // |Le| mass (elements only)
    std::vector<Index> degree(un, Index(0));       // approximate external degree
    std::vector<char>  is_element(un, 0);
    std::vector<char>  alive(un, 1);               // principal variable, not yet eliminated

    for (std::size_t i = 0; i < un; ++i) {
        const std::size_t b = static_cast<std::size_t>(adj_ptr[i]);
        const std::size_t e = static_cast<std::size_t>(adj_ptr[i + 1u]);
        A_var[i].assign(adj_ind.begin() + b, adj_ind.begin() + e);
        degree[i] = static_cast<Index>(A_var[i].size());  // nv == 1 initially
        mem_head[i] = static_cast<Index>(i);   // singleton chain {i}
        mem_tail[i] = static_cast<Index>(i);
    }

    // SLU-OQ1: ordered candidate set (design D-1) -- since ORD-F2 Phase A1 a
    // lazy-deletion binary heap with the same total order (see the helper
    // block above): at least one (degree[i], i) entry per selectable vertex,
    // stale entries discarded at pop time.
    std::vector<std::pair<Index, Index> > cand;
    cand.reserve(un);
    for (Index i = Index(0); i < n; ++i) {
        cand.push_back(std::make_pair(degree[static_cast<std::size_t>(i)], i));
    }
    std::make_heap(cand.begin(), cand.end(), sparse_lu_ordering_pq_greater_());

    std::vector<Index> order;
    order.reserve(un);

    // Per-pivot scratch for the |Le \ Lp| set-difference computation.
    std::vector<Index> elem_w(un, Index(0));
    std::vector<Index> elem_stamp(un, Index(-1));
    Index stamp = Index(0);

    // Per-pivot Lp membership marker (reset via lp_touched after each pivot).
    std::vector<char> in_lp(un, 0);

    // ORD-F2 Phase A2: per-pivot scratch hoisted out of the pivot loop so the
    // capacities are reused across pivots (cleared / re-assigned at the top of
    // each use; contents per pivot are identical to the former loop-local
    // vectors).
    std::vector<Index> lp;
    std::vector<Index> lp_touched;
    std::vector<Index> absorbed;
    std::vector<Index> rmv;
    std::vector<Index> touched_elems;
    std::vector<Index> agg;
    std::vector<Index> merged;
    std::vector<std::pair<unsigned long long, Index> > hb;

    Index eliminated = Index(0);

    while (eliminated < n) {
        // ---- pivot selection: min approximate degree, tie-break smallest index.
        // SLU-OQ1 / ORD-F2 Phase A1: the first VALID popped entry
        // (alive && !is_element && key == degree[v]) is the lexicographic min
        // of (degree, index) among selectable vertices -- identical to
        // *cand.begin() of the former std::set (design F1), so the emitted
        // permutation is bit-identical.  Adopting the popped entry replaces
        // the former cand.erase(cand.begin()) (T-A1).
        Index p = Index(-1);
        while (!cand.empty()) {
            const std::pair<Index, Index> top = cand.front();
            std::pop_heap(cand.begin(), cand.end(),
                          sparse_lu_ordering_pq_greater_());
            cand.pop_back();
            const std::size_t uv = static_cast<std::size_t>(top.second);
            if (alive[uv] && !is_element[uv] && top.first == degree[uv]) {
                p = top.second;  // p becomes an element below (T-A1)
                break;
            }
        }
        if (p < Index(0)) break;  // defensive (was: cand.empty() fallthrough)
        const std::size_t up = static_cast<std::size_t>(p);

        // ---- form Lp = (A_var[p]) U (U_{e in E_elem[p]} V_elem[e]) \ {p}.
        lp.clear();
        lp_touched.clear();
        for (std::size_t t = 0; t < A_var[up].size(); ++t) {
            const Index j = A_var[up][t];
            const std::size_t uj = static_cast<std::size_t>(j);
            if (alive[uj] && !is_element[uj] && j != p && !in_lp[uj]) {
                in_lp[uj] = 1; lp.push_back(j); lp_touched.push_back(j);
            }
        }
        for (std::size_t te = 0; te < E_elem[up].size(); ++te) {
            const Index e = E_elem[up][te];
            const std::vector<Index>& Ve = V_elem[static_cast<std::size_t>(e)];
            for (std::size_t t = 0; t < Ve.size(); ++t) {
                const Index j = Ve[t];
                const std::size_t uj = static_cast<std::size_t>(j);
                if (alive[uj] && !is_element[uj] && j != p && !in_lp[uj]) {
                    in_lp[uj] = 1; lp.push_back(j); lp_touched.push_back(j);
                }
            }
        }
        std::sort(lp.begin(), lp.end());

        // Elements absorbed exactly into the new element ep := p.
        absorbed.assign(E_elem[up].begin(), E_elem[up].end());  // sorted

        // ---- turn p into element ep.  Its members are already eliminated below.
        is_element[up] = 1;
        alive[up] = 0;
        V_elem[up] = lp;
        Index lp_mass = Index(0);
        for (std::size_t t = 0; t < lp.size(); ++t) {
            lp_mass += nv[static_cast<std::size_t>(lp[t])];
        }
        elem_size[up] = lp_mass;

        // Removal set for variable lists: Lp U {p}.
        rmv.assign(lp.begin(), lp.end());
        sparse_lu_amd_sorted_insert(rmv, p);

        // ---- rewrite each i in Lp: drop redundant variable edges (Lp, p), drop
        // absorbed elements, attach ep.
        for (std::size_t t = 0; t < lp.size(); ++t) {
            const std::size_t ui = static_cast<std::size_t>(lp[t]);
            sparse_lu_amd_sorted_remove_set_inplace_(A_var[ui], rmv);
            sparse_lu_amd_sorted_remove_set_inplace_(E_elem[ui], absorbed);
            sparse_lu_amd_sorted_insert(E_elem[ui], p);
        }
        // Retire absorbed elements (only Lp referenced them).
        for (std::size_t t = 0; t < absorbed.size(); ++t) {
            const std::size_t ue = static_cast<std::size_t>(absorbed[t]);
            V_elem[ue].clear();
            elem_size[ue] = Index(0);
        }

        // ---- set differences elem_w[e] = |Le \ Lp| for elements e adjacent to
        // Lp (e != ep), via the standard stamp trick.
        ++stamp;
        touched_elems.clear();
        for (std::size_t t = 0; t < lp.size(); ++t) {
            const std::size_t ui = static_cast<std::size_t>(lp[t]);
            const std::vector<Index>& Ei = E_elem[ui];
            for (std::size_t s = 0; s < Ei.size(); ++s) {
                const Index e = Ei[s];
                if (e == p) continue;
                const std::size_t ue = static_cast<std::size_t>(e);
                if (elem_stamp[ue] != stamp) {
                    elem_stamp[ue] = stamp;
                    elem_w[ue] = elem_size[ue];
                    touched_elems.push_back(e);
                }
                elem_w[ue] -= nv[ui];
            }
        }

        // ---- aggressive element absorption: |Le \ Lp| == 0  =>  Le subset of Lp
        // (all of e's variables are in ep), so e is redundant; absorb into ep.
        agg.clear();
        for (std::size_t s = 0; s < touched_elems.size(); ++s) {
            const Index e = touched_elems[s];
            if (elem_w[static_cast<std::size_t>(e)] == Index(0)) agg.push_back(e);
        }
        if (!agg.empty()) {
            std::sort(agg.begin(), agg.end());
            for (std::size_t t = 0; t < lp.size(); ++t) {
                sparse_lu_amd_sorted_remove_set_inplace_(
                    E_elem[static_cast<std::size_t>(lp[t])], agg);
            }
            for (std::size_t s = 0; s < agg.size(); ++s) {
                const std::size_t ue = static_cast<std::size_t>(agg[s]);
                V_elem[ue].clear();
                elem_size[ue] = Index(0);
            }
        }

        // ---- supervariable (indistinguishable-variable) detection within Lp.
        // Two variables with identical variable-lists AND element-lists are
        // indistinguishable; merge the larger index into the smaller (S-5).
        // Bucket by a structural hash via a sorted (hash, index) vector, then
        // compare exactly within a run of equal hashes.
        if (lp.size() > 1u) {
            hb.clear();
            hb.reserve(lp.size());
            for (std::size_t t = 0; t < lp.size(); ++t) {
                const Index i = lp[t];
                const std::size_t ui = static_cast<std::size_t>(i);
                if (!alive[ui]) continue;
                unsigned long long h = 1469598103934665603ull;
                for (std::size_t s = 0; s < A_var[ui].size(); ++s) {
                    h = (h ^ static_cast<unsigned long long>(A_var[ui][s] + 1))
                        * 1099511628211ull;
                }
                h = h * 31ull + 0x9e3779b97f4a7c15ull;
                for (std::size_t s = 0; s < E_elem[ui].size(); ++s) {
                    h = (h ^ static_cast<unsigned long long>(E_elem[ui][s] + 2))
                        * 1099511628211ull;
                }
                hb.push_back(std::make_pair(h, i));
            }
            // Sort by (hash, index): equal-hash candidates become contiguous and
            // ascending in index, so the principal (smallest index) comes first.
            std::sort(hb.begin(), hb.end());
            for (std::size_t a = 0; a < hb.size(); ++a) {
                const std::size_t ui = static_cast<std::size_t>(hb[a].second);
                if (!alive[ui]) continue;
                merged.clear();
                for (std::size_t b = a + 1;
                     b < hb.size() && hb[b].first == hb[a].first; ++b) {
                    const Index j = hb[b].second;
                    const std::size_t uj = static_cast<std::size_t>(j);
                    if (!alive[uj]) continue;
                    if (!sparse_lu_amd_lists_equal(A_var[ui], A_var[uj])) continue;
                    if (!sparse_lu_amd_lists_equal(E_elem[ui], E_elem[uj])) continue;
                    // Merge j (larger index) into i (smaller index).  The
                    // V_elem / A_var removals of j are deferred to the batch
                    // below (ORD-F2 Phase A2).
                    nv[ui] += nv[uj];
                    // ORD-F3 P3: append j's whole chain at i's tail -- the
                    // same concatenation order as the former
                    // members[ui].insert(end, begin, end).  No clear of j's
                    // chain is needed: j goes !alive here and is never
                    // selected as pivot nor as principal again, so its chain
                    // is only ever reached through i from now on.
                    mem_next[static_cast<std::size_t>(mem_tail[ui])] =
                        mem_head[uj];
                    mem_tail[ui] = mem_tail[uj];
                    merged.push_back(j);
                    // T-A2: no cand retire needed (lazy heap) -- flipping
                    // alive[uj] below stales every stored entry of j.
                    alive[uj] = 0;
                    nv[uj] = Index(0);
                    A_var[uj].clear();
                    E_elem[uj].clear();
                }
                // ORD-F2 Phase A2: batched removal of the merged j's.  Every
                // merged j had A_var[uj] == A_var[ui] and E_elem[uj] ==
                // E_elem[ui] (the merge condition), and the principal's two
                // lists never change during this run (i is not a member of
                // its own lists, so no removal above targets them).  Hence
                // one remove-set pass per target list over the ascending
                // `merged` values performs exactly the removals the former
                // per-j sorted_remove_one calls did.  Deferral cannot flip a
                // later verdict in the run: a candidate whose lists a
                // deferred removal would touch is adjacent to a merged j,
                // i.e. it appears in A_var[ui] while never containing itself
                // -- such a candidate fails lists_equal against the
                // principal in both orderings.
                if (!merged.empty()) {
                    const std::vector<Index>& Ei = E_elem[ui];
                    for (std::size_t s = 0; s < Ei.size(); ++s) {
                        sparse_lu_amd_sorted_remove_set_inplace_(
                            V_elem[static_cast<std::size_t>(Ei[s])], merged);
                    }
                    const std::vector<Index>& Ai = A_var[ui];
                    for (std::size_t s = 0; s < Ai.size(); ++s) {
                        sparse_lu_amd_sorted_remove_set_inplace_(
                            A_var[static_cast<std::size_t>(Ai[s])], merged);
                    }
                }
            }
        }

        // ---- approximate external degree for surviving principals in Lp:
        // d_i = min( n - elim_after - nv[i],          (A) remaining-variable bound
        //            deg_old[i] + (|Lp| - nv[i]),     (B) previous-degree bound
        //            |Ai| + (|Lp| - nv[i]) + sum_e |Le \ Lp| )  (C) approx external
        const Index elim_after = eliminated + nv[up];
        for (std::size_t t = 0; t < lp.size(); ++t) {
            const std::size_t ui = static_cast<std::size_t>(lp[t]);
            if (!alive[ui]) continue;
            const Index ln = elem_size[up] - nv[ui];  // |Lp \ {i}| in mass
            Index a_mass = Index(0);
            for (std::size_t s = 0; s < A_var[ui].size(); ++s) {
                a_mass += nv[static_cast<std::size_t>(A_var[ui][s])];
            }
            Index e_diff = Index(0);
            for (std::size_t s = 0; s < E_elem[ui].size(); ++s) {
                const Index e = E_elem[ui][s];
                if (e == p) continue;
                e_diff += elem_w[static_cast<std::size_t>(e)];
            }
            Index dC = a_mass + ln + e_diff;
            Index dB = degree[ui] + ln;
            Index dA = n - elim_after - nv[ui];
            if (dA < Index(0)) dA = Index(0);
            Index d = dC;
            if (dB < d) d = dB;
            if (dA < d) d = dA;
            if (d < Index(0)) d = Index(0);
            sparse_lu_ordering_cand_update_degree_(cand, degree, lp[t], d,
                                                   alive, is_element, n);  // W-A2
        }

        // ---- emit the original variables represented by p, in member order
        // (ORD-F3 P3: walk the chain head -> next; same sequence as the
        // former members[up] vector iteration).
        for (Index s = mem_head[up]; s >= Index(0);
             s = mem_next[static_cast<std::size_t>(s)]) {
            order.push_back(s);
        }
        eliminated += nv[up];

        // reset Lp markers for the next pivot.
        for (std::size_t s = 0; s < lp_touched.size(); ++s) {
            in_lp[static_cast<std::size_t>(lp_touched[s])] = 0;
        }
    }

    // Coverage invariant: every vertex eliminated exactly once.
    if (order.size() != un) {
        vcp::throw_error<vcp::state_error>(
            "sparse_lu_amd_ordering: elimination did not cover all vertices");
    }

    // Elimination order is the fill-reducing permutation directly.  perm[new] = old.
    for (std::size_t i = 0; i < un; ++i) {
        perm[i] = order[i];
    }
    return perm;
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
