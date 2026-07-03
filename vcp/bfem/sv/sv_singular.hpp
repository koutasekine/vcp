// vcp/bfem/sv/sv_singular.hpp
// Phase 5d (Scott-Vogelius parts): exact detection of singular vertices (2D)
// and singular edges (3D) with deterministic fan-component decomposition
// (V2 detection half; external design section 4.1, internal design section 2).
//
// Normative rules:
//  - 2D vertex v is SINGULAR (within one fan component) iff the direction
//    lines of the edges meeting at v fall into <= 2 +-collinear classes
//    (exact 2x2 determinant signs);
//  - 3D edge e is SINGULAR (within one fan component) iff the planes of the
//    faces meeting at e fall into <= 2 +-collinear normal classes (exact
//    cross products);
//  - fan components (v0.2, B-1): the incident elements are decomposed into
//    connected components by edge (2D) / face (3D) adjacency THROUGH the
//    entity before anything else; pinch configurations are handled one
//    component at a time. The walk is deterministic: a chain starts at its
//    smallest-id endpoint element; a cycle starts at its smallest-id element
//    and proceeds toward the smaller-id neighbor.
//  - SV-4: when a sign needed by the collinearity decisions cannot be
//    certified (interval coordinates straddling zero without being the exact
//    point zero), sv_indeterminate_singularity is thrown (fail-fast; the
//    space itself cannot be defined). Exact and point types always decide.
//  - detection depends on the geometry ONLY (never on the degree n); the
//    number of constraint rows is the only n-dependent quantity (sv_rows).
//
// No division operator appears in this file (G-SV-1).

#ifndef VCP_BFEM_SV_SV_SINGULAR_HPP
#define VCP_BFEM_SV_SV_SINGULAR_HPP

#include <vector>
#include <array>
#include <map>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dofmap.hpp>            // detail::mesh_topology2
#include <vcp/bfem/d3/topology3.hpp>      // detail::mesh_topology3

namespace vcp {
namespace bfem {

// SV-4 fail-fast: the collinearity of interval coordinates is undecidable
class sv_indeterminate_singularity : public std::runtime_error {
public:
    explicit sv_indeterminate_singularity(const std::string& msg)
        : std::runtime_error(msg) {}
};

namespace detail {

// exact sign with certain comparisons: -1 / 0 / +1, or 2 when the sign is
// not certified (for intervals, x <= 0 <= x certifies the exact point zero)
template <typename T>
int sv_sign_tri(const T& x) {
    if (x < T(0)) return -1;
    if (T(0) < x) return 1;
    if (x <= T(0) && T(0) <= x) return 0;
    return 2;
}

template <typename T>
bool sv_collinear2(const T& ux, const T& uy, const T& vx, const T& vy) {
    int s = sv_sign_tri(ux * vy - uy * vx);
    if (s == 2)
        throw sv_indeterminate_singularity(
            "bfem::sv: 2D collinearity undecidable on interval input");
    return s == 0;
}

template <typename T>
bool sv_collinear3(const std::array<T, 3>& u, const std::array<T, 3>& v) {
    int sx = sv_sign_tri(u[1] * v[2] - u[2] * v[1]);
    if (sx == -1 || sx == 1) return false;
    int sy = sv_sign_tri(u[2] * v[0] - u[0] * v[2]);
    if (sy == -1 || sy == 1) return false;
    int sz = sv_sign_tri(u[0] * v[1] - u[1] * v[0]);
    if (sz == -1 || sz == 1) return false;
    if (sx == 0 && sy == 0 && sz == 0) return true;
    throw sv_indeterminate_singularity(
        "bfem::sv: 3D collinearity undecidable on interval input");
}

// ---------------------------------------------------------------------------
// one fan component around a vertex (2D) / an edge (3D)
// ---------------------------------------------------------------------------
struct sv_fan {
    std::vector<int> elems;    // deterministic fan order (chain or cycle)
    bool closed;               // true: cycle (interior fan)
    int classes;               // +-collinear direction / plane classes
};

// deterministic decomposition + walk over an incidence structure:
//   star_elems : incident element ids (ascending)
//   ent_of     : for each star element, the ids of the (exactly two)
//                through-entities (edges at v / faces at e) it contributes
inline std::vector<std::vector<int> > sv_fan_walk(
    const std::vector<int>& star_elems,
    const std::vector<std::array<int, 2> >& ent_of) {
    const std::size_t ns = star_elems.size();
    // through-entity -> star positions
    std::map<int, std::vector<std::size_t> > ent2pos;
    for (std::size_t s = 0; s < ns; ++s) {
        ent2pos[ent_of[s][0]].push_back(s);
        ent2pos[ent_of[s][1]].push_back(s);
    }
    // neighbor lists (<= 2 by manifoldness of the through-entities)
    std::vector<std::vector<std::size_t> > nbr(ns);
    for (std::map<int, std::vector<std::size_t> >::const_iterator it =
             ent2pos.begin(); it != ent2pos.end(); ++it) {
        const std::vector<std::size_t>& ps = it->second;
        if (ps.size() == 2) {
            nbr[ps[0]].push_back(ps[1]);
            nbr[ps[1]].push_back(ps[0]);
        }
    }
    // connected components (ascending seed order)
    std::vector<char> seen(ns, 0);
    std::vector<std::vector<int> > fans;
    for (std::size_t seed = 0; seed < ns; ++seed) {
        if (seen[seed]) continue;
        // collect the component
        std::vector<std::size_t> comp;
        std::vector<std::size_t> stack;
        stack.push_back(seed);
        seen[seed] = 1;
        while (!stack.empty()) {
            std::size_t cur = stack.back();
            stack.pop_back();
            comp.push_back(cur);
            for (std::size_t k = 0; k < nbr[cur].size(); ++k) {
                std::size_t nx = nbr[cur][k];
                if (!seen[nx]) { seen[nx] = 1; stack.push_back(nx); }
            }
        }
        std::sort(comp.begin(), comp.end());
        // endpoints: component members with fewer than 2 in-component
        // neighbors (through-entity unshared inside the star)
        std::vector<std::size_t> ends;
        for (std::size_t k = 0; k < comp.size(); ++k)
            if (nbr[comp[k]].size() < 2) ends.push_back(comp[k]);
        std::vector<int> order;
        if (comp.size() == 1) {
            order.push_back(star_elems[comp[0]]);
        } else if (!ends.empty()) {
            // chain: start at the smallest-id endpoint (comp is ascending,
            // star_elems ascending, so the first endpoint is the smallest)
            std::size_t cur = ends[0];
            std::size_t prev = ns;                     // sentinel
            for (;;) {
                order.push_back(star_elems[cur]);
                std::size_t nxt = ns;
                for (std::size_t k = 0; k < nbr[cur].size(); ++k)
                    if (nbr[cur][k] != prev) nxt = nbr[cur][k];
                if (nxt == ns) break;
                prev = cur;
                cur = nxt;
            }
        } else {
            // cycle: start at the smallest element, toward the smaller
            // neighbor
            std::size_t start = comp[0];
            std::size_t nxt = nbr[start][0] < nbr[start][1] ? nbr[start][0]
                                                            : nbr[start][1];
            std::size_t prev = start;
            std::size_t cur = nxt;
            order.push_back(star_elems[start]);
            while (cur != start) {
                order.push_back(star_elems[cur]);
                std::size_t n2 = (nbr[cur][0] != prev) ? nbr[cur][0] : nbr[cur][1];
                prev = cur;
                cur = n2;
            }
        }
        fans.push_back(order);
    }
    return fans;
}

// local position of global vertex v inside element e (2D)
inline int sv_local_vertex_pos2(const mesh_topology2& tp, int e, int v) {
    const std::array<int, 3>& t = tp.tri[static_cast<std::size_t>(e)];
    for (int p = 0; p < 3; ++p)
        if (t[static_cast<std::size_t>(p)] == v) return p;
    assert(false);
    return -1;
}

// ---------------------------------------------------------------------------
// 2D: fans (with direction class counts) around vertex v
// ---------------------------------------------------------------------------
template <typename T>
std::vector<sv_fan> sv_vertex_fans(const mesh<2, T>& msh,
                                   const mesh_topology2& tp,
                                   const std::vector<int>& star_elems, int v) {
    const std::size_t ns = star_elems.size();
    std::vector<std::array<int, 2> > ent_of(ns);
    for (std::size_t s = 0; s < ns; ++s) {
        int e = star_elems[s];
        int p = sv_local_vertex_pos2(tp, e, v);
        // the two local edges THROUGH v are the edges opposite the other two
        // local vertices
        int k1 = (p + 1) % 3;
        int k2 = (p + 2) % 3;
        ent_of[s][0] = tp.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(k1)];
        ent_of[s][1] = tp.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(k2)];
    }
    std::vector<std::vector<int> > orders = sv_fan_walk(star_elems, ent_of);

    std::vector<sv_fan> fans;
    for (std::size_t f = 0; f < orders.size(); ++f) {
        sv_fan fan;
        fan.elems = orders[f];
        // closed <=> every through-edge of the component is shared by two
        // component elements <=> #edges == #elems (cycle); chain has one more
        std::vector<int> eids;
        for (std::size_t s = 0; s < fan.elems.size(); ++s) {
            int e = fan.elems[s];
            int p = sv_local_vertex_pos2(tp, e, v);
            for (int dk = 1; dk <= 2; ++dk) {
                int eid = tp.tri_edge[static_cast<std::size_t>(e)]
                                     [static_cast<std::size_t>((p + dk) % 3)];
                eids.push_back(eid);
            }
        }
        std::sort(eids.begin(), eids.end());
        eids.erase(std::unique(eids.begin(), eids.end()), eids.end());
        fan.closed = (eids.size() == fan.elems.size());
        // direction classes of the through-edges
        std::vector<std::array<T, 2> > reps;
        for (std::size_t s = 0; s < eids.size(); ++s) {
            const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(eids[s])];
            int other = (ev[0] == v) ? ev[1] : ev[0];
            std::array<T, 2> dir;
            dir[0] = msh.vertex(other)[0] - msh.vertex(v)[0];
            dir[1] = msh.vertex(other)[1] - msh.vertex(v)[1];
            bool matched = false;
            for (std::size_t r = 0; r < reps.size() && !matched; ++r)
                if (sv_collinear2(reps[r][0], reps[r][1], dir[0], dir[1]))
                    matched = true;
            if (!matched) reps.push_back(dir);
        }
        fan.classes = static_cast<int>(reps.size());
        fans.push_back(fan);
    }
    return fans;
}

// local positions (p, q) of the global edge (a < b) inside tet e:
// p holds a (the smaller global id), q holds b
inline void sv_local_edge_pos3(const std::array<int, 4>& tet, int a, int b,
                               int& p, int& q) {
    p = -1;
    q = -1;
    for (int i = 0; i < 4; ++i) {
        if (tet[static_cast<std::size_t>(i)] == a) p = i;
        if (tet[static_cast<std::size_t>(i)] == b) q = i;
    }
    assert(p >= 0 && q >= 0);
}

// ---------------------------------------------------------------------------
// 3D: fans (with plane class counts) around edge ed
// ---------------------------------------------------------------------------
template <typename T>
std::vector<sv_fan> sv_edge_fans(const mesh<3, T>& msh,
                                 const mesh_topology3& tp,
                                 const std::vector<int>& star_elems, int ed) {
    const int a = tp.edges[static_cast<std::size_t>(ed)][0];
    const int b = tp.edges[static_cast<std::size_t>(ed)][1];
    const std::size_t ns = star_elems.size();
    std::vector<std::array<int, 2> > ent_of(ns);
    for (std::size_t s = 0; s < ns; ++s) {
        int e = star_elems[s];
        int p, q;
        sv_local_edge_pos3(tp.tet[static_cast<std::size_t>(e)], a, b, p, q);
        // the two local faces THROUGH the edge are the faces opposite the
        // other two local vertices
        int cnt = 0;
        for (int k = 0; k < 4; ++k) {
            if (k == p || k == q) continue;
            ent_of[s][static_cast<std::size_t>(cnt)] =
                tp.tet_face[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)];
            ++cnt;
        }
    }
    std::vector<std::vector<int> > orders = sv_fan_walk(star_elems, ent_of);

    std::vector<sv_fan> fans;
    for (std::size_t f = 0; f < orders.size(); ++f) {
        sv_fan fan;
        fan.elems = orders[f];
        std::vector<int> fids;
        for (std::size_t s = 0; s < fan.elems.size(); ++s) {
            int e = fan.elems[s];
            int p, q;
            sv_local_edge_pos3(tp.tet[static_cast<std::size_t>(e)], a, b, p, q);
            for (int k = 0; k < 4; ++k) {
                if (k == p || k == q) continue;
                fids.push_back(tp.tet_face[static_cast<std::size_t>(e)]
                                          [static_cast<std::size_t>(k)]);
            }
        }
        std::sort(fids.begin(), fids.end());
        fids.erase(std::unique(fids.begin(), fids.end()), fids.end());
        fan.closed = (fids.size() == fan.elems.size());
        // plane classes via face normals (exact cross products)
        std::vector<std::array<T, 3> > reps;
        for (std::size_t s = 0; s < fids.size(); ++s) {
            const std::array<int, 3>& fv =
                tp.faces[static_cast<std::size_t>(fids[s])];
            std::array<T, 3> u, w, n;
            for (int d = 0; d < 3; ++d) {
                u[static_cast<std::size_t>(d)] =
                    msh.vertex(fv[1])[static_cast<std::size_t>(d)]
                    - msh.vertex(fv[0])[static_cast<std::size_t>(d)];
                w[static_cast<std::size_t>(d)] =
                    msh.vertex(fv[2])[static_cast<std::size_t>(d)]
                    - msh.vertex(fv[0])[static_cast<std::size_t>(d)];
            }
            n[0] = u[1] * w[2] - u[2] * w[1];
            n[1] = u[2] * w[0] - u[0] * w[2];
            n[2] = u[0] * w[1] - u[1] * w[0];
            bool matched = false;
            for (std::size_t r = 0; r < reps.size() && !matched; ++r)
                if (sv_collinear3(reps[r], n)) matched = true;
            if (!matched) reps.push_back(n);
        }
        fan.classes = static_cast<int>(reps.size());
        fans.push_back(fan);
    }
    return fans;
}

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_SV_SINGULAR_HPP
