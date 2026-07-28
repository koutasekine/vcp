// vcp/bfem/meshgen_delaunay.hpp
// MG-1 Phase B: Lawson flips on the coarse (pre-refinement) mesh, provided
// through generate_mesh_delaunay (design section 5.3, rulings Q6/Q7).
//
// One-way dependency: this header includes meshgen.hpp; meshgen.hpp does not
// know about this header.
//
// Guarantees (design 5.3): only edges with a *certain* incircle violation
// are flipped; boundary edges and bridge-origin edges are never flipped;
// termination is enforced by a flip-count cap of 8 x (internal edge count),
// reported through meshgen_status::delaunay_complete (never an exception).
// The incircle sign uses meshgen_traits<T>::sign; an indefinite sign defers
// the flip -- the catch below is the single catch permitted in the layer.

#ifndef VCP_BFEM_MESHGEN_DELAUNAY_HPP
#define VCP_BFEM_MESHGEN_DELAUNAY_HPP

#include <vcp/bfem/meshgen.hpp>

namespace vcp {
namespace bfem {

namespace meshgen_detail {

// incircle determinant, d-based translation (additions/multiplications only):
//   | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2 |
//   | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2 |
//   | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2 |
// positive iff d lies strictly inside the circumcircle of the CCW (a,b,c)
template <typename T>
T incircle_det(const std::array<T, 2>& a, const std::array<T, 2>& b,
               const std::array<T, 2>& c, const std::array<T, 2>& d) {
    const T ax = a[0] - d[0];
    const T ay = a[1] - d[1];
    const T bx = b[0] - d[0];
    const T by = b[1] - d[1];
    const T cx = c[0] - d[0];
    const T cy = c[1] - d[1];
    const T aq = ax * ax + ay * ay;
    const T bq = bx * bx + by * by;
    const T cq = cx * cx + cy * cy;
    return ax * (by * cq - bq * cy)
         - ay * (bx * cq - bq * cx)
         + aq * (bx * cy - by * cx);
}

// coarse mesh after Lawson flips, plus the constraint-edge set used
template <typename T>
struct delaunay_coarse {
    coarse_mesh<T> cm;                                // triangles after flips
    std::map<std::pair<int, int>, bool> constrained;  // input segs + bridges
    int flips;
    bool complete;
    delaunay_coarse() : cm(), constrained(), flips(0), complete(true) {}
};

const int meshgen_flip_cap_factor = 8;   // cap = 8 x (internal edge count)

template <typename T>
delaunay_coarse<T> delaunay_flip_coarse(const polygon_domain<T>& dom,
                                        const meshgen_options& opt) {
    delaunay_coarse<T> r;
    std::vector<int> walk;
    r.cm = build_coarse_mesh(dom, opt, walk);
    std::vector<std::array<int, 3> >& tris = r.cm.triangles;
    const std::vector<std::array<T, 2> >& V = r.cm.vertices;

    // constraint edges: input boundary segments ...
    for (std::map<std::pair<int, int>, meshgen_facet_source>::const_iterator
             it = r.cm.segment_of.begin(); it != r.cm.segment_of.end(); ++it)
        r.constrained[it->first] = true;
    // ... and bridge edges (vertex pairs traversed twice by the bridged walk)
    {
        std::map<std::pair<int, int>, int> wc;
        const int nw = static_cast<int>(walk.size());
        for (int i = 0; i < nw; ++i)
            wc[edge_key(walk[static_cast<std::size_t>(i)],
                        walk[static_cast<std::size_t>((i + 1) % nw)])] += 1;
        for (std::map<std::pair<int, int>, int>::const_iterator it = wc.begin();
             it != wc.end(); ++it)
            if (it->second > 1) r.constrained[it->first] = true;
    }

    // internal edge count (adjacency count == 2); constant under flips
    int ninternal = 0;
    {
        std::map<std::pair<int, int>, int> cnt;
        for (std::size_t t = 0; t < tris.size(); ++t)
            for (int k = 0; k < 3; ++k)
                cnt[edge_key(tris[t][static_cast<std::size_t>(k)],
                             tris[t][static_cast<std::size_t>((k + 1) % 3)])] += 1;
        for (std::map<std::pair<int, int>, int>::const_iterator it = cnt.begin();
             it != cnt.end(); ++it)
            if (it->second == 2) ++ninternal;
    }
    const int cap = meshgen_flip_cap_factor * ninternal;

    while (true) {
        // edge -> list of (triangle index, position of the directed edge)
        std::map<std::pair<int, int>, std::vector<std::pair<int, int> > > adj;
        for (std::size_t t = 0; t < tris.size(); ++t)
            for (int k = 0; k < 3; ++k)
                adj[edge_key(tris[t][static_cast<std::size_t>(k)],
                             tris[t][static_cast<std::size_t>((k + 1) % 3)])]
                    .push_back(std::pair<int, int>(static_cast<int>(t), k));

        bool flipped = false;
        for (std::map<std::pair<int, int>,
                      std::vector<std::pair<int, int> > >::const_iterator
                 it = adj.begin(); it != adj.end() && !flipped; ++it) {
            if (it->second.size() != 2) continue;
            if (r.constrained.find(it->first) != r.constrained.end()) continue;
            const int t1 = it->second[0].first;
            const int k1 = it->second[0].second;
            const int t2 = it->second[1].first;
            const int k2 = it->second[1].second;
            // t1 = (u,v,p) CCW with directed edge (u,v); t2 = (v,u,q) CCW
            const int u = tris[static_cast<std::size_t>(t1)]
                              [static_cast<std::size_t>(k1)];
            const int v = tris[static_cast<std::size_t>(t1)]
                              [static_cast<std::size_t>((k1 + 1) % 3)];
            const int p = tris[static_cast<std::size_t>(t1)]
                              [static_cast<std::size_t>((k1 + 2) % 3)];
            const int q = tris[static_cast<std::size_t>(t2)]
                              [static_cast<std::size_t>((k2 + 2) % 3)];
            // certain incircle violation? indefinite sign defers the flip
            // (the single permitted catch of the meshgen layer)
            int s = 0;
            bool decided = true;
            try {
                s = meshgen_traits<T>::sign(incircle_det(
                    V[static_cast<std::size_t>(u)], V[static_cast<std::size_t>(v)],
                    V[static_cast<std::size_t>(p)], V[static_cast<std::size_t>(q)]));
            } catch (const meshgen_degenerate&) {
                decided = false;
            }
            if (!decided || s != +1) continue;
            // validity guard: both replacement triangles must be certainly CCW
            const std::array<T, 2>& U = V[static_cast<std::size_t>(u)];
            const std::array<T, 2>& W = V[static_cast<std::size_t>(v)];
            const std::array<T, 2>& P = V[static_cast<std::size_t>(p)];
            const std::array<T, 2>& Q = V[static_cast<std::size_t>(q)];
            if (!certainly_pos(orient2d(U, Q, P))
             || !certainly_pos(orient2d(Q, W, P)))
                continue;
            if (r.flips == cap) {
                r.complete = false;
                return r;
            }
            std::array<int, 3> n1 = {{u, q, p}};
            std::array<int, 3> n2 = {{q, v, p}};
            tris[static_cast<std::size_t>(t1)] = n1;
            tris[static_cast<std::size_t>(t2)] = n2;
            ++r.flips;
            flipped = true;
        }
        if (!flipped) break;   // no certainly-violating flippable edge remains
    }
    return r;
}

} // namespace meshgen_detail

// ---------------------------------------------------------------------------
// public entry point, Phase B (design section 6):
// ear clipping -> Lawson flips on the coarse mesh -> the Phase A refinement
// and boundary-output path. The final (refined) mesh's CDT property is not
// claimed (Q6); status.delaunay_complete reports flip termination.
// ---------------------------------------------------------------------------
template <typename T>
mesh<2, T> generate_mesh_delaunay(const polygon_domain<T>& dom, const T& h,
                                  meshgen_status<2>& status,
                                  const meshgen_options& opt = meshgen_options()) {
    if (!meshgen_detail::certainly_pos(h))
        throw meshgen_error("vcp::bfem::meshgen: h must be certainly positive");
    status = meshgen_status<2>();
    meshgen_detail::delaunay_coarse<T> dc =
        meshgen_detail::delaunay_flip_coarse(dom, opt);
    std::vector<std::array<T, 2> > vertices;
    std::vector<std::array<int, 3> > elements;
    meshgen_detail::refine_and_finalize(dc.cm, h * h, opt, status,
                                        vertices, elements);
    status.delaunay_complete = dc.complete;
    return mesh<2, T>::from_lists(vertices, elements);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_DELAUNAY_HPP
