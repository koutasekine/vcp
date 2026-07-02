// vcp/bfem/d3/topology3.hpp
// Phase 5a (3D common infrastructure): tetrahedral mesh topology (detail) --
// edge/face tables, face orientation permutation sigma + parity, boundary
// face detection, and THE single implementation point of the face convention.
//
// Conforms to: 3D common external design v0.1 (sections 2, 3) and
//              3D common internal design v0.1 (section 2).
//
// Face convention (external design section 2, normative):
//  - canonical face  = ascending global triple (g0 < g1 < g2), faces ordered
//    lexicographically by that triple;
//  - local face k    = the face opposite local vertex k; its local vertex
//    order is the ascending local numbers excluding k (no outward-normal
//    normalization -- sigma and its parity absorb the orientation);
//  - sigma           = the S_3 permutation sending local position i to the
//    canonical position of the global vertex seen at local position i;
//    parity is published alongside (unused by P^n, consumed by 3D RT).
//
// Non-manifold check: a FACE shared by 3+ elements is rejected
// (std::invalid_argument). Edges are intentionally NOT multiplicity-checked:
// in 3D an edge may be shared by arbitrarily many elements (AX3D-3).
//
// The 2D topology (vcp/bfem/dofmap.hpp) is frozen; this is a parallel
// D = 3 implementation, no 2D code is touched.

#ifndef VCP_BFEM_D3_TOPOLOGY3_HPP
#define VCP_BFEM_D3_TOPOLOGY3_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/d3/s3_perm.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// static local tables (normative):
//  - local edges in the lexicographic order of the local vertex pairs
//    01, 02, 03, 12, 13, 23 (internal design 2.1)
//  - local face k = ascending local numbers excluding k (external design 2-2)
// ---------------------------------------------------------------------------
struct tet_local {
    static int edge_vertex(int le, int side) {           // le in [0,6), side in {0,1}
        assert(le >= 0 && le < 6 && side >= 0 && side < 2);
        static const int ev[6][2] = {
            {0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
        };
        return ev[le][side];
    }
    static int edge_of_pair(int p, int q) {              // p < q local
        assert(p >= 0 && q > p && q < 4);
        static const int pe[4][4] = {                    // pe[p][q], p < q
            {-1, 0, 1, 2}, {-1, -1, 3, 4}, {-1, -1, -1, 5}, {-1, -1, -1, -1}
        };
        return pe[p][q];
    }
    static int face_vertex(int k, int i) {               // local face k, position i
        assert(k >= 0 && k < 4 && i >= 0 && i < 3);
        static const int fv[4][3] = {
            {1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}
        };
        return fv[k][i];
    }
};

// ---------------------------------------------------------------------------
// THE single implementation point of the face convention (internal design
// section 5, S-5A-2): a degree-m 3D multi-index alpha with alpha_k == 0 is
// mapped to the face multi-index (b0, b1, b2) in the CANONICAL (ascending
// global) vertex order of the face:
//   beta_local[i] = alpha[ face_vertex(k, i) ]   (local ascending order)
//   beta_canon[sigma(i)] = beta_local[i]         (s3_perm::apply)
// Both dof_build3 (face-interior block rank) and trace_index3 (face trace)
// call this function; the convention has no second implementation.
// ---------------------------------------------------------------------------
inline std::array<int, 3> face_canonical_beta(const multi_index<3>& alpha,
                                              int k, int perm_code) {
    assert(alpha.a[static_cast<std::size_t>(k)] == 0);
    std::array<int, 3> beta_local;
    for (int i = 0; i < 3; ++i)
        beta_local[static_cast<std::size_t>(i)] =
            alpha.a[static_cast<std::size_t>(tet_local::face_vertex(k, i))];
    std::array<int, 3> beta_canon;
    s3_perm::apply(perm_code, beta_local, beta_canon);
    return beta_canon;
}

// ---------------------------------------------------------------------------
// topology (D = 3): edge and face tables and boundary detection;
// coordinates unused.
// ---------------------------------------------------------------------------
struct mesh_topology3 {
    int nv, nt;
    std::vector<std::array<int, 4> > tet;
    std::vector<std::array<int, 2> > edges;            // edge id -> (vmin, vmax), lex sorted
    std::vector<std::array<int, 6> > tet_edge;         // element -> edge ids (local edge order 01,02,03,12,13,23)
    std::vector<std::array<bool, 6> > tet_edge_fwd;    // local direction (p -> q, p < q local) == canonical?
    std::vector<std::array<int, 3> > faces;            // face id -> (g0, g1, g2) ascending, lex sorted
    std::vector<std::array<int, 4> > tet_face;         // element -> face ids (local face k)
    std::vector<std::array<int, 4> > tet_face_perm;    // sigma as s3_perm code 0..5
    std::vector<std::array<int, 4> > tet_face_parity;  // +1 / -1
    std::vector<std::array<int, 3> > face_edge;        // face id -> edge ids of (g0,g1), (g0,g2), (g1,g2)
    std::vector<int> boundary_faces;                   // faces with exactly 1 element (ascending)

    template <typename T>
    static mesh_topology3 build(const mesh<3, T>& m) {
        mesh_topology3 tp;
        tp.nv = m.num_vertices();
        tp.nt = m.num_elements();
        tp.tet.resize(static_cast<std::size_t>(tp.nt));
        for (int e = 0; e < tp.nt; ++e) {
            tp.tet[static_cast<std::size_t>(e)] = m.element(e);
            const std::array<int, 4>& t = tp.tet[static_cast<std::size_t>(e)];
            for (int a = 0; a < 4; ++a)
                for (int b = a + 1; b < 4; ++b)
                    if (t[static_cast<std::size_t>(a)] == t[static_cast<std::size_t>(b)])
                        throw std::invalid_argument(
                            "bfem::topology3: element with repeated vertex");
        }
        build_edges(tp);
        build_faces(tp);
        return tp;
    }

    int num_edges() const { return static_cast<int>(edges.size()); }
    int num_faces() const { return static_cast<int>(faces.size()); }

private:
    struct erec { int vmin, vmax, e, le; };
    struct erec_less_t {
        bool operator()(const erec& a, const erec& b) const {
            if (a.vmin != b.vmin) return a.vmin < b.vmin;
            if (a.vmax != b.vmax) return a.vmax < b.vmax;
            if (a.e != b.e) return a.e < b.e;          // determinism inside a key
            return a.le < b.le;
        }
    };
    struct frec { int g0, g1, g2, e, k; };
    struct frec_less_t {
        bool operator()(const frec& a, const frec& b) const {
            if (a.g0 != b.g0) return a.g0 < b.g0;
            if (a.g1 != b.g1) return a.g1 < b.g1;
            if (a.g2 != b.g2) return a.g2 < b.g2;
            if (a.e != b.e) return a.e < b.e;          // determinism inside a key
            return a.k < b.k;
        }
    };

    static void build_edges(mesh_topology3& tp) {
        std::vector<erec> rs;
        rs.reserve(static_cast<std::size_t>(6 * tp.nt));
        for (int e = 0; e < tp.nt; ++e) {
            const std::array<int, 4>& t = tp.tet[static_cast<std::size_t>(e)];
            for (int le = 0; le < 6; ++le) {
                int a = t[static_cast<std::size_t>(tet_local::edge_vertex(le, 0))];
                int b = t[static_cast<std::size_t>(tet_local::edge_vertex(le, 1))];
                erec r;
                r.vmin = a < b ? a : b;
                r.vmax = a < b ? b : a;
                r.e = e;
                r.le = le;
                rs.push_back(r);
            }
        }
        std::sort(rs.begin(), rs.end(), erec_less_t());
        tp.tet_edge.resize(static_cast<std::size_t>(tp.nt));
        tp.tet_edge_fwd.resize(static_cast<std::size_t>(tp.nt));
        int ne = 0;
        std::size_t i = 0;
        while (i < rs.size()) {
            std::size_t j = i;
            while (j < rs.size() && rs[j].vmin == rs[i].vmin && rs[j].vmax == rs[i].vmax) ++j;
            // NOTE: no multiplicity limit here -- in 3D an edge may be shared
            // by any number of elements (AX3D-3); only FACES carry the
            // non-manifold check (build_faces).
            std::array<int, 2> ev = { { rs[i].vmin, rs[i].vmax } };
            tp.edges.push_back(ev);
            for (std::size_t s = i; s < j; ++s) {
                int e = rs[s].e, le = rs[s].le;
                tp.tet_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(le)] = ne;
                const std::array<int, 4>& t = tp.tet[static_cast<std::size_t>(e)];
                // local positive direction: p -> q with (p, q) = local pair, p < q
                tp.tet_edge_fwd[static_cast<std::size_t>(e)][static_cast<std::size_t>(le)] =
                    (t[static_cast<std::size_t>(tet_local::edge_vertex(le, 0))]
                     < t[static_cast<std::size_t>(tet_local::edge_vertex(le, 1))]);
            }
            ++ne;
            i = j;
        }
    }

    // edge id of the global pair (a < b) by binary search over the lex-sorted
    // edge table (used to fill face_edge)
    static int edge_id_of(const mesh_topology3& tp, int a, int b) {
        assert(a < b);
        std::array<int, 2> key = { { a, b } };
        std::vector<std::array<int, 2> >::const_iterator it =
            std::lower_bound(tp.edges.begin(), tp.edges.end(), key);
        assert(it != tp.edges.end() && (*it)[0] == a && (*it)[1] == b);
        return static_cast<int>(it - tp.edges.begin());
    }

    static void build_faces(mesh_topology3& tp) {
        std::vector<frec> rs;
        rs.reserve(static_cast<std::size_t>(4 * tp.nt));
        for (int e = 0; e < tp.nt; ++e) {
            const std::array<int, 4>& t = tp.tet[static_cast<std::size_t>(e)];
            for (int k = 0; k < 4; ++k) {
                int g[3];
                for (int i = 0; i < 3; ++i)
                    g[i] = t[static_cast<std::size_t>(tet_local::face_vertex(k, i))];
                frec r;
                r.g0 = g[0]; r.g1 = g[1]; r.g2 = g[2];
                if (r.g0 > r.g1) std::swap(r.g0, r.g1);
                if (r.g1 > r.g2) std::swap(r.g1, r.g2);
                if (r.g0 > r.g1) std::swap(r.g0, r.g1);
                r.e = e;
                r.k = k;
                rs.push_back(r);
            }
        }
        std::sort(rs.begin(), rs.end(), frec_less_t());
        tp.tet_face.resize(static_cast<std::size_t>(tp.nt));
        tp.tet_face_perm.resize(static_cast<std::size_t>(tp.nt));
        tp.tet_face_parity.resize(static_cast<std::size_t>(tp.nt));
        int nf = 0;
        std::size_t i = 0;
        while (i < rs.size()) {
            std::size_t j = i;
            while (j < rs.size() && rs[j].g0 == rs[i].g0 && rs[j].g1 == rs[i].g1
                   && rs[j].g2 == rs[i].g2) ++j;
            std::size_t count = j - i;
            if (count > 2)
                throw std::invalid_argument(
                    "bfem::topology3: non-manifold input (a face shared by 3+ elements)");
            std::array<int, 3> fv = { { rs[i].g0, rs[i].g1, rs[i].g2 } };
            tp.faces.push_back(fv);
            std::array<int, 3> fe = { { edge_id_of(tp, fv[0], fv[1]),
                                        edge_id_of(tp, fv[0], fv[2]),
                                        edge_id_of(tp, fv[1], fv[2]) } };
            tp.face_edge.push_back(fe);
            if (count == 1) tp.boundary_faces.push_back(nf);
            for (std::size_t s = i; s < j; ++s) {
                int e = rs[s].e, k = rs[s].k;
                tp.tet_face[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)] = nf;
                const std::array<int, 4>& t = tp.tet[static_cast<std::size_t>(e)];
                // sigma: local position i (ascending local order) -> canonical
                // position of its global vertex = #{ j : G_j < G_i }
                int G[3];
                for (int p = 0; p < 3; ++p)
                    G[p] = t[static_cast<std::size_t>(tet_local::face_vertex(k, p))];
                int s0 = (G[1] < G[0]) + (G[2] < G[0]);
                int s1 = (G[0] < G[1]) + (G[2] < G[1]);
                int s2 = (G[0] < G[2]) + (G[1] < G[2]);
                int code = s3_perm::from_images(s0, s1, s2);
                tp.tet_face_perm[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)] = code;
                tp.tet_face_parity[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)] =
                    s3_perm::parity(code);
            }
            ++nf;
            i = j;
        }
        // boundary_faces built in ascending face id order (already sorted)
    }
};

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_D3_TOPOLOGY3_HPP
