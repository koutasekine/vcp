// vcp/bfem/meshgen_adapt.hpp
// MG-2 (A): conforming local refinement by newest vertex bisection (NVB)
// on the coarse (pre-uniform-refinement) mesh, design MG-2 v1.1 section 4.
//
// One-way dependency: this header includes meshgen.hpp; meshgen.hpp does
// not know about this header (the meshgen_delaunay.hpp precedent).
//
// Core properties (design 4.1):
//  - the bisection / conformity-propagation decisions are PURELY
//    combinatorial (integer ids only) -- no geometric predicate is
//    evaluated after the initial labeling, so every scalar type T behaves
//    identically once the labels are fixed;
//  - geometry enters ONLY through the initial labeling (longest edge by
//    squared-length comparison; adopted conventions below);
//  - vertices are identified by combinatorial keys only; an edge midpoint
//    is created once per undirected edge (smaller-id endpoint first in the
//    midpoint formula) and shared by both incident triangles;
//  - boundary lineage goes through the shared primitive
//    meshgen_detail::split_boundary_lineage of meshgen.hpp (MG-2a), the
//    same implementation the uniform path derives its sweeps from.
//
// Initial labeling conventions (design 4.2; recorded in the report):
//  - refinement edge = the longest edge by certain dist2 comparisons;
//  - a certain TIE is broken toward the lexicographically SMALLEST
//    unordered edge key (min id, max id);
//  - if any needed comparison is not certain (interval types), the WHOLE
//    triangle falls back to the id rule: refinement edge = the edge with
//    the lexicographically smallest unordered edge key. The fallback count
//    is recorded in nvb_state::fallback_labels (T-A7).
//
// Termination safety net: the conformity chain marks triangles on the
// current chain; revisiting one would mean a refinement-edge cycle, a
// state the design argues unreachable ([design claim, T-A7]); it raises
// meshgen_degenerate like the other internal guards of the layer.

#ifndef VCP_BFEM_MESHGEN_ADAPT_HPP
#define VCP_BFEM_MESHGEN_ADAPT_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>

#include <vcp/bfem/meshgen.hpp>
#include <vcp/bfem/meshgen_delaunay.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// NVB bookkeeping. peak[t] is the position (0..2) of the newest vertex of
// triangle t inside cm.triangles[t]; the refinement edge is the opposite
// edge (peak+1, peak+2). edge_tris is the active edge -> incident-triangle
// adjacency, maintained incrementally across bisections. The state stays
// valid across refine_* calls as long as the caller does not modify the
// coarse mesh by other means.
// ---------------------------------------------------------------------------
template <typename T>
struct nvb_state {
    std::vector<int> peak;
    std::vector<int> gen;     // bisection generation per triangle (root = 0)
    std::map<std::pair<int, int>, std::vector<int> > edge_tris;
    bool labeled;
    int fallback_labels;      // triangles labeled by the id fallback (T-A7)
    long long bisections;     // total single-triangle bisections performed
    nvb_state() : peak(), gen(), edge_tris(), labeled(false),
                  fallback_labels(0), bisections(0) {}
};

namespace meshgen_adapt_detail {

using meshgen_detail::coarse_mesh;
using meshgen_detail::edge_key;
using meshgen_detail::dist2;
using meshgen_detail::certainly_less;
using meshgen_detail::certainly_equal;

// lexicographic order on unordered edge keys (int comparison only)
inline bool edge_key_less(const std::pair<int, int>& a,
                          const std::pair<int, int>& b) {
    if (a.first != b.first) return a.first < b.first;
    return a.second < b.second;
}

// initial labeling of one triangle (conventions in the header comment).
// Returns the refinement-edge position k (edge k = (v_k, v_{k+1})) and
// sets fallback = true when the id rule was used.
template <typename T>
int initial_refedge(const coarse_mesh<T>& cm, const std::array<int, 3>& e,
                    bool& fallback) {
    fallback = false;
    T d[3] = {T(0), T(0), T(0)};
    std::pair<int, int> key[3];
    for (int k = 0; k < 3; ++k) {
        const int u = e[static_cast<std::size_t>(k)];
        const int v = e[static_cast<std::size_t>((k + 1) % 3)];
        d[k] = dist2(cm.vertices[static_cast<std::size_t>(u)],
                     cm.vertices[static_cast<std::size_t>(v)]);
        key[k] = edge_key(u, v);
    }
    int best = 0;
    for (int k = 1; k < 3; ++k) {
        if (certainly_less(d[best], d[k])) { best = k; continue; }
        if (certainly_less(d[k], d[best])) continue;
        if (certainly_equal(d[k], d[best])) {
            if (edge_key_less(key[k], key[best])) best = k;
            continue;
        }
        // not certain: id fallback for the whole triangle
        fallback = true;
        break;
    }
    if (fallback) {
        best = 0;
        for (int k = 1; k < 3; ++k)
            if (edge_key_less(key[k], key[best])) best = k;
    }
    return best;
}

// the incremental NVB engine: wraps one refine_* call. target flags follow
// the marked lineage (child1 keeps the parent's id and thus its flag;
// child2 copies it); split_round records the sweep in which an id was last
// replaced by its child1, so a round never splits the same lineage twice.
template <typename T>
struct nvb_engine {
    coarse_mesh<T>& cm;
    nvb_state<T>& nvb;
    std::vector<char> is_target;
    std::vector<int> split_round;
    std::vector<char> in_chain;
    int round;

    nvb_engine(coarse_mesh<T>& c, nvb_state<T>& s)
        : cm(c), nvb(s), is_target(), split_round(), in_chain(), round(0) {
        const std::size_t n = cm.triangles.size();
        is_target.assign(n, 0);
        split_round.assign(n, -1);
        in_chain.assign(n, 0);
    }

    void adjacency_remove(const std::pair<int, int>& key, int t) {
        std::map<std::pair<int, int>, std::vector<int> >::iterator it =
            nvb.edge_tris.find(key);
        if (it == nvb.edge_tris.end())
            throw meshgen_degenerate(
                "vcp::bfem::meshgen_adapt: internal: missing adjacency edge");
        std::vector<int>& lst = it->second;
        bool removed = false;
        for (std::size_t i = 0; i < lst.size(); ++i)
            if (lst[i] == t) {
                lst.erase(lst.begin() + static_cast<std::ptrdiff_t>(i));
                removed = true;
                break;
            }
        if (!removed)
            throw meshgen_degenerate(
                "vcp::bfem::meshgen_adapt: internal: missing adjacency entry");
        if (lst.empty()) nvb.edge_tris.erase(it);
    }

    void adjacency_add(const std::pair<int, int>& key, int t) {
        std::vector<int>& lst = nvb.edge_tris[key];
        lst.push_back(t);
        if (lst.size() > 2)
            throw meshgen_degenerate(
                "vcp::bfem::meshgen_adapt: internal: edge with more than two "
                "incident triangles");
    }

    // active triangle sharing edge `key` other than t, or -1
    int neighbor_across(const std::pair<int, int>& key, int t) const {
        std::map<std::pair<int, int>, std::vector<int> >::const_iterator it =
            nvb.edge_tris.find(key);
        if (it == nvb.edge_tris.end())
            throw meshgen_degenerate(
                "vcp::bfem::meshgen_adapt: internal: missing adjacency edge");
        const std::vector<int>& lst = it->second;
        for (std::size_t i = 0; i < lst.size(); ++i)
            if (lst[i] != t) return lst[i];
        return -1;
    }

    std::pair<int, int> refedge_key(int t) const {
        const std::array<int, 3>& e =
            cm.triangles[static_cast<std::size_t>(t)];
        const int p = nvb.peak[static_cast<std::size_t>(t)];
        return edge_key(e[static_cast<std::size_t>((p + 1) % 3)],
                        e[static_cast<std::size_t>((p + 2) % 3)]);
    }

    // split triangle x against the already-created midpoint m of its
    // refinement edge; child1 replaces x (same id), child2 is appended
    void bisect_one(int x, int m) {
        const std::array<int, 3> e = cm.triangles[static_cast<std::size_t>(x)];
        const int p = nvb.peak[static_cast<std::size_t>(x)];
        const int A = e[static_cast<std::size_t>(p)];
        const int B = e[static_cast<std::size_t>((p + 1) % 3)];
        const int C = e[static_cast<std::size_t>((p + 2) % 3)];

        adjacency_remove(edge_key(A, B), x);
        adjacency_remove(edge_key(B, C), x);
        adjacency_remove(edge_key(C, A), x);

        // child1 = (A, B, m), newest vertex m at position 2 -> refinement
        // edge (A, B); child2 = (A, m, C), newest vertex at position 1 ->
        // refinement edge (C, A): the parent's two non-refinement edges
        const int z = static_cast<int>(cm.triangles.size());
        std::array<int, 3> c1 = {{A, B, m}};
        std::array<int, 3> c2 = {{A, m, C}};
        cm.triangles[static_cast<std::size_t>(x)] = c1;
        nvb.peak[static_cast<std::size_t>(x)] = 2;
        cm.triangles.push_back(c2);
        nvb.peak.push_back(1);
        nvb.gen[static_cast<std::size_t>(x)] =
            nvb.gen[static_cast<std::size_t>(x)] + 1;
        nvb.gen.push_back(nvb.gen[static_cast<std::size_t>(x)]);
        is_target.push_back(is_target[static_cast<std::size_t>(x)]);
        split_round[static_cast<std::size_t>(x)] = round;
        split_round.push_back(round);
        in_chain.push_back(0);

        adjacency_add(edge_key(A, B), x);
        adjacency_add(edge_key(B, m), x);
        adjacency_add(edge_key(m, A), x);
        adjacency_add(edge_key(A, m), z);
        adjacency_add(edge_key(m, C), z);
        adjacency_add(edge_key(C, A), z);

        nvb.bisections = nvb.bisections + 1;
    }

    // split the refinement edge of x (and of its compatible partner, if
    // any): midpoint creation (smaller-id endpoint first), boundary
    // lineage through the shared MG-2a primitive, then both bisections
    void split_pair(int x) {
        const std::pair<int, int> ek = refedge_key(x);
        const int n = neighbor_across(ek, x);
        if (n >= 0 && refedge_key(n) != ek)
            throw meshgen_degenerate(
                "vcp::bfem::meshgen_adapt: internal: partner not compatible");
        const int u = ek.first;    // u < v by edge_key
        const int v = ek.second;
        std::array<T, 2> c;
        c[0] = (cm.vertices[static_cast<std::size_t>(u)][0]
              + cm.vertices[static_cast<std::size_t>(v)][0]) / T(2);
        c[1] = (cm.vertices[static_cast<std::size_t>(u)][1]
              + cm.vertices[static_cast<std::size_t>(v)][1]) / T(2);
        const int m = static_cast<int>(cm.vertices.size());
        cm.vertices.push_back(c);
        meshgen_detail::split_boundary_lineage(cm.segment_of, u, v, m);
        bisect_one(x, m);
        if (n >= 0) bisect_one(n, m);
    }

    // conforming bisection of element t: follow the refinement-edge chain
    // until a compatible pair (or a boundary refinement edge), split it,
    // then unwind. Purely combinatorial.
    void bisect_element(int t) {
        std::vector<int> chain;
        int cur = t;
        while (true) {
            const std::pair<int, int> ek = refedge_key(cur);
            const int n = neighbor_across(ek, cur);
            if (n < 0 || refedge_key(n) == ek) break;
            if (in_chain[static_cast<std::size_t>(n)] != 0)
                throw meshgen_degenerate(
                    "vcp::bfem::meshgen_adapt: refinement-edge cycle "
                    "(unreachable by the labeling design; report if seen)");
            in_chain[static_cast<std::size_t>(cur)] = 1;
            chain.push_back(cur);
            cur = n;
        }
        split_pair(cur);
        while (!chain.empty()) {
            const int prev = chain.back();
            chain.pop_back();
            in_chain[static_cast<std::size_t>(prev)] = 0;
            const std::pair<int, int> ek = refedge_key(prev);
            const int n = neighbor_across(ek, prev);
            if (n >= 0 && refedge_key(n) != ek)
                throw meshgen_degenerate(
                    "vcp::bfem::meshgen_adapt: internal: chain unwind not "
                    "compatible");
            split_pair(prev);
        }
    }
};

// rebuild-from-scratch adjacency + initial labeling (first refine_* call)
template <typename T>
void nvb_initialize(const coarse_mesh<T>& cm, nvb_state<T>& nvb) {
    nvb.peak.clear();
    nvb.gen.assign(cm.triangles.size(), 0);
    nvb.edge_tris.clear();
    nvb.fallback_labels = 0;
    for (std::size_t t = 0; t < cm.triangles.size(); ++t) {
        const std::array<int, 3>& e = cm.triangles[t];
        bool fb = false;
        const int k = initial_refedge<T>(cm, e, fb);
        if (fb) nvb.fallback_labels = nvb.fallback_labels + 1;
        nvb.peak.push_back((k + 2) % 3);
        for (int j = 0; j < 3; ++j) {
            std::vector<int>& lst = nvb.edge_tris[
                edge_key(e[static_cast<std::size_t>(j)],
                         e[static_cast<std::size_t>((j + 1) % 3)])];
            lst.push_back(static_cast<int>(t));
            if (lst.size() > 2)
                throw meshgen_error(
                    "vcp::bfem::meshgen_adapt: non-manifold input mesh");
        }
    }
    nvb.labeled = true;
}

} // namespace meshgen_adapt_detail

// ---------------------------------------------------------------------------
// primary marking primitive (design 4.3): bisect every element of the
// marked lineage once per round, `rounds` times, keeping the mesh
// conforming through NVB propagation. Element ids are positions in
// cm.triangles at call time; a bisected element keeps its id for one child
// and appends the other, and BOTH children stay in the marked lineage for
// the following rounds. The state is labeled lazily on the first call.
// ---------------------------------------------------------------------------
template <typename T>
void refine_marked(meshgen_detail::coarse_mesh<T>& cm, nvb_state<T>& nvb,
                   const std::vector<int>& marked_elements, int rounds) {
    if (rounds < 0)
        throw meshgen_error(
            "vcp::bfem::meshgen_adapt: rounds must be non-negative");
    if (!nvb.labeled)
        meshgen_adapt_detail::nvb_initialize(cm, nvb);
    if (nvb.peak.size() != cm.triangles.size() ||
        nvb.gen.size() != cm.triangles.size())
        throw meshgen_error(
            "vcp::bfem::meshgen_adapt: nvb_state out of sync with the mesh");

    meshgen_adapt_detail::nvb_engine<T> eng(cm, nvb);
    const int n0 = static_cast<int>(cm.triangles.size());
    for (std::size_t i = 0; i < marked_elements.size(); ++i) {
        const int t = marked_elements[i];
        if (t < 0 || t >= n0)
            throw meshgen_error(
                "vcp::bfem::meshgen_adapt: marked element id out of range");
        eng.is_target[static_cast<std::size_t>(t)] = 1;
    }

    for (int r = 0; r < rounds; ++r) {
        eng.round = r;
        std::vector<int> work;
        for (std::size_t t = 0; t < cm.triangles.size(); ++t)
            if (eng.is_target[t] != 0) work.push_back(static_cast<int>(t));
        for (std::size_t i = 0; i < work.size(); ++i) {
            const int t = work[i];
            // already replaced by its child1 in this round (either as an
            // earlier work item or by chain propagation): goal met
            if (eng.split_round[static_cast<std::size_t>(t)] == r) continue;
            eng.bisect_element(t);
        }
    }
}

namespace meshgen_adapt_detail {

// element centroid (the layer's only new geometric quantity; design C5)
template <typename T>
std::array<T, 2> centroid(const coarse_mesh<T>& cm,
                          const std::array<int, 3>& e) {
    std::array<T, 2> c;
    for (int j = 0; j < 2; ++j)
        c[static_cast<std::size_t>(j)] =
            (cm.vertices[static_cast<std::size_t>(e[0])][static_cast<std::size_t>(j)]
           + cm.vertices[static_cast<std::size_t>(e[1])][static_cast<std::size_t>(j)]
           + cm.vertices[static_cast<std::size_t>(e[2])][static_cast<std::size_t>(j)])
            / T(3);
    return c;
}

} // namespace meshgen_adapt_detail

// ---------------------------------------------------------------------------
// convenience 1 (design 4.3): size-field driven marking. h_of_x is called
// with the element centroid and returns the local target size as a T; an
// element is marked when some squared edge length CERTAINLY exceeds
// h(c)^2 (a not-certain comparison never marks -- conservative, exception
// free). At most max_rounds sweeps; stops early when nothing is marked.
// ---------------------------------------------------------------------------
template <typename T, typename SizeFn>
void refine_size_field(meshgen_detail::coarse_mesh<T>& cm, nvb_state<T>& nvb,
                       SizeFn h_of_x, int max_rounds) {
    if (max_rounds < 0)
        throw meshgen_error(
            "vcp::bfem::meshgen_adapt: max_rounds must be non-negative");
    for (int r = 0; r < max_rounds; ++r) {
        std::vector<int> marked;
        for (std::size_t t = 0; t < cm.triangles.size(); ++t) {
            const std::array<int, 3>& e = cm.triangles[t];
            const std::array<T, 2> c =
                meshgen_adapt_detail::centroid(cm, e);
            const T hx = h_of_x(c);
            const T hh = hx * hx;
            for (int k = 0; k < 3; ++k) {
                const T d2 = meshgen_adapt_detail::dist2(
                    cm.vertices[static_cast<std::size_t>(e[static_cast<std::size_t>(k)])],
                    cm.vertices[static_cast<std::size_t>(e[static_cast<std::size_t>((k + 1) % 3)])]);
                if (meshgen_adapt_detail::certainly_less(hh, d2)) {
                    marked.push_back(static_cast<int>(t));
                    break;
                }
            }
        }
        if (marked.empty()) break;
        refine_marked(cm, nvb, marked, 1);
    }
}

// ---------------------------------------------------------------------------
// convenience 2, the main entrance (design 4.3, rulings Q5/C5): geometric
// grading toward re-entrant corners. Layer k (k = 1..layers) is the ball
// dist2(centroid, corner) < sigma^(2k) around any listed corner; every
// pass bisects the current elements of that ball once, so an element in
// the k-th ball ends up bisected k times: the refinement depth grows by
// one per layer inward. Comparisons are certain-only (a not-certain
// membership never marks). sigma is the user's geometric ratio (design
// 4.3: not built in, certainly inside (0,1)).
// ---------------------------------------------------------------------------
template <typename T>
void refine_geometric_corner(meshgen_detail::coarse_mesh<T>& cm,
                             nvb_state<T>& nvb,
                             const std::vector<int>& corner_vertex_ids,
                             const T& sigma, int layers) {
    if (layers < 0)
        throw meshgen_error(
            "vcp::bfem::meshgen_adapt: layers must be non-negative");
    if (!meshgen_detail::certainly_pos(sigma) ||
        !meshgen_adapt_detail::certainly_less(sigma, T(1)))
        throw meshgen_error(
            "vcp::bfem::meshgen_adapt: sigma must be certainly inside (0,1)");
    // corner coordinates are frozen up front (ids denote input vertices)
    std::vector<std::array<T, 2> > corners;
    for (std::size_t i = 0; i < corner_vertex_ids.size(); ++i) {
        const int v = corner_vertex_ids[i];
        if (v < 0 || v >= static_cast<int>(cm.vertices.size()))
            throw meshgen_error(
                "vcp::bfem::meshgen_adapt: corner vertex id out of range");
        corners.push_back(cm.vertices[static_cast<std::size_t>(v)]);
    }
    const T s2 = sigma * sigma;
    T ball = T(1);
    for (int k = 1; k <= layers; ++k) {
        ball = ball * s2;      // sigma^(2k)
        std::vector<int> marked;
        for (std::size_t t = 0; t < cm.triangles.size(); ++t) {
            const std::array<T, 2> c =
                meshgen_adapt_detail::centroid(cm, cm.triangles[t]);
            for (std::size_t i = 0; i < corners.size(); ++i) {
                if (meshgen_adapt_detail::certainly_less(
                        meshgen_adapt_detail::dist2(c, corners[i]), ball)) {
                    marked.push_back(static_cast<int>(t));
                    break;
                }
            }
        }
        if (!marked.empty()) refine_marked(cm, nvb, marked, 1);
    }
}

// ---------------------------------------------------------------------------
// unified adaptive entrance (design 4.5): local refinement inserted BEFORE
// the uniform sweeps. An empty spec reproduces the default (or Delaunay)
// path EXACTLY -- it delegates to the very same entry points (T-A6 /
// G-MG2-0). The explicit marks are applied first, then the corner grading.
// ---------------------------------------------------------------------------
template <typename T>
struct meshgen_adapt_spec {
    bool use_delaunay;                   // Lawson flips before the NVB stage
    std::vector<int> marked_elements;    // explicit coarse-element marks
    int marked_rounds;                   // rounds for the explicit marks
    std::vector<int> corner_vertex_ids;  // grading corners (input vertex ids)
    T sigma;                             // geometric ratio, certainly in (0,1)
    int layers;                          // grading layers
    meshgen_adapt_spec()
        : use_delaunay(false), marked_elements(), marked_rounds(0),
          corner_vertex_ids(), sigma(T(0)), layers(0) {}
    bool empty() const {
        return marked_elements.empty() && corner_vertex_ids.empty();
    }
};

template <typename T>
mesh<2, T> generate_mesh_adaptive(const polygon_domain<T>& dom, const T& h,
                                  const meshgen_adapt_spec<T>& spec,
                                  meshgen_status<2>& status,
                                  const meshgen_options& opt = meshgen_options()) {
    if (spec.empty()) {
        // exactly the existing entry points (T-A6): no code of this header
        // touches the mesh
        if (spec.use_delaunay)
            return generate_mesh_delaunay(dom, h, status, opt);
        return generate_mesh(dom, h, status, opt);
    }
    if (!meshgen_detail::certainly_pos(h))
        throw meshgen_error("vcp::bfem::meshgen: h must be certainly positive");
    status = meshgen_status<2>();
    meshgen_detail::coarse_mesh<T> cm;
    bool complete = true;
    if (spec.use_delaunay) {
        meshgen_detail::delaunay_coarse<T> dc =
            meshgen_detail::delaunay_flip_coarse(dom, opt);
        cm = dc.cm;
        complete = dc.complete;
    } else {
        cm = meshgen_detail::build_coarse_mesh(dom, opt);
    }
    nvb_state<T> nvb;
    if (!spec.marked_elements.empty())
        refine_marked(cm, nvb, spec.marked_elements, spec.marked_rounds);
    if (!spec.corner_vertex_ids.empty())
        refine_geometric_corner(cm, nvb, spec.corner_vertex_ids,
                                spec.sigma, spec.layers);
    std::vector<std::array<T, 2> > vertices;
    std::vector<std::array<int, 3> > elements;
    meshgen_detail::refine_and_finalize(cm, h * h, opt, status,
                                        vertices, elements);
    status.delaunay_complete = complete;
    return mesh<2, T>::from_lists(vertices, elements);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_ADAPT_HPP
