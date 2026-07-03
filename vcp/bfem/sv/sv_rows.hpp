// vcp/bfem/sv/sv_rows.hpp
// Phase 5d (Scott-Vogelius parts): constraint row candidates, exact
// independization, and the public sv_pressure_constraints class with the
// sv_stability_hint observation API (V2 generation half; external design
// sections 2, 4.2-4.4, internal design section 3).
//
// Row generation (normative):
//  - 2D singular vertex v, fan T_0..T_{m-1}: one row
//        sum_i (-1)^i [coefficient of q at v on T_i] = 0
//    where the vertex coefficient of broken P_l is the element-block dof
//    e_i * N_l + rank(l * e_{p_i}) (p_i = local position of v in T_i);
//  - 3D singular edge (a < b), fan T_0..T_{m-1}: one row per trace position
//    t = 0..l (that is n = l+1 rows), with the edge-trace coefficient
//    alpha = (l - t) e_{p_i} + t e_{q_i} (p_i holds a, q_i holds b -- the
//    canonical direction rule: t counts toward the LARGER global vertex,
//    mirroring the dofmap convention);
//  - every candidate row is normalized deterministically (dofs ascending,
//    the coefficient of the smallest dof positive), which removes the
//    residual freedom of the fan walk (start / direction);
//  - independization: candidates are sorted by the canonical key
//    (kind, entity id, fan component, trace position) and selected GREEDILY;
//    the independence test is the exact rational rank of rational_la over
//    the compressed support columns (internal design section 3-2). The
//    output is always an independent set of the ORIGINAL integer rows, and
//    is invariant under any permutation of the candidate input (S-SV-5).
//
// SV-1 (non-rejection): sv_pressure_constraints never refuses an input mesh;
// sv_stability_hint reports the known sufficient stability conditions as an
// OBSERVATION (enum + reason string), with the mesh provenance supplied by
// the caller (v0.2, B-2 -- deciding "Alfeld or not" from a raw mesh is
// unsound, so the knowledge of having called alfeld_refine is passed in).
//
// No division operator appears in this file (G-SV-1); the exact rank runs
// inside the frozen rational_la.

#ifndef VCP_BFEM_SV_SV_ROWS_HPP
#define VCP_BFEM_SV_SV_ROWS_HPP

#include <vector>
#include <array>
#include <map>
#include <string>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/d3/topology3.hpp>
#include <vcp/bfem/rt/rational_la.hpp>
#include <vcp/bfem/sv/sv_constraint.hpp>
#include <vcp/bfem/sv/sv_singular.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// sv_stability_hint (SV-1, B-2): observation only, no claim of stability
// ---------------------------------------------------------------------------
enum sv_mesh_provenance {
    sv_provenance_general,       // nothing is known about the mesh origin
    sv_provenance_alfeld         // the caller built the mesh by alfeld_refine
};

enum sv_hint_status {
    sv_hint_stable_general,      // known sufficient condition on general meshes
    sv_hint_stable_alfeld,       // known sufficient condition on Alfeld meshes
    sv_hint_unknown              // outside the recorded sufficient conditions
};

struct sv_stability_hint {
    sv_hint_status status;
    std::string reason;
};

namespace detail {

// candidate row with its canonical sort key
struct sv_candidate {
    int kind;                    // 0 = vertex row, 1 = edge row
    int entity;                  // vertex id / edge id
    int comp;                    // fan component index (pinch decomposition)
    int t;                       // trace position (edge rows; 0 for vertices)
    sv_constraint row;
};

inline bool sv_candidate_less(const sv_candidate& a, const sv_candidate& b) {
    if (a.kind != b.kind) return a.kind < b.kind;
    if (a.entity != b.entity) return a.entity < b.entity;
    if (a.comp != b.comp) return a.comp < b.comp;
    return a.t < b.t;
}

// deterministic normalization: dofs ascending, first coefficient positive
inline void sv_normalize_row(sv_constraint& r) {
    std::vector<std::pair<int, int> > p;
    p.reserve(r.dof.size());
    for (std::size_t k = 0; k < r.dof.size(); ++k)
        p.push_back(std::make_pair(r.dof[k], r.coef[k]));
    std::sort(p.begin(), p.end());
    int flip = (!p.empty() && p[0].second < 0) ? -1 : 1;
    for (std::size_t k = 0; k < p.size(); ++k) {
        r.dof[k] = p[k].first;
        r.coef[k] = flip * p[k].second;
    }
}

// exact independence of integer rows through the frozen rational_la rank
// over the compressed union of supports
inline bool sv_extends_rank(const std::vector<const sv_constraint*>& acc,
                            const sv_constraint& cand) {
    std::map<int, int> col;
    for (std::size_t r = 0; r < acc.size(); ++r)
        for (std::size_t k = 0; k < acc[r]->dof.size(); ++k)
            col.insert(std::make_pair(acc[r]->dof[k], 0));
    for (std::size_t k = 0; k < cand.dof.size(); ++k)
        col.insert(std::make_pair(cand.dof[k], 0));
    int c = 0;
    for (std::map<int, int>::iterator it = col.begin(); it != col.end(); ++it) {
        it->second = c;
        ++c;
    }
    rmat A(static_cast<int>(acc.size()) + 1, c);
    for (std::size_t r = 0; r < acc.size(); ++r)
        for (std::size_t k = 0; k < acc[r]->dof.size(); ++k)
            A.at(static_cast<int>(r), col[acc[r]->dof[k]]) =
                rational(acc[r]->coef[k]);
    for (std::size_t k = 0; k < cand.dof.size(); ++k)
        A.at(static_cast<int>(acc.size()), col[cand.dof[k]]) =
            rational(cand.coef[k]);
    return rank_exact(A) == static_cast<int>(acc.size()) + 1;
}

// canonical-order greedy selection of an independent subset; the result does
// not depend on the input permutation (S-SV-5)
inline std::vector<sv_candidate> sv_independize(std::vector<sv_candidate> cand) {
    std::sort(cand.begin(), cand.end(), sv_candidate_less);
    std::vector<sv_candidate> out;
    std::vector<const sv_constraint*> acc;
    for (std::size_t k = 0; k < cand.size(); ++k) {
        if (sv_extends_rank(acc, cand[k].row)) {
            out.push_back(cand[k]);
            acc.push_back(&out.back().row);
            // out reallocation would invalidate acc: rebuild (small sizes)
            acc.clear();
            for (std::size_t r = 0; r < out.size(); ++r)
                acc.push_back(&out[r].row);
        }
    }
    return out;
}

// ---------------------------------------------------------------------------
// candidate generation
// ---------------------------------------------------------------------------

// 2D: all vertex-row candidates for broken P_l (l = n - 1)
template <typename T>
std::vector<sv_candidate> sv_vertex_candidates(const mesh<2, T>& msh,
                                               const mesh_topology2& tp,
                                               int l,
                                               std::vector<int>* singular_out) {
    const index_map<2>& im = coeff_registry<2>::indices(l);
    const int Nl = im.size();
    // vertex stars
    std::vector<std::vector<int> > star(static_cast<std::size_t>(tp.nv));
    for (int e = 0; e < tp.nt; ++e)
        for (int p = 0; p < 3; ++p)
            star[static_cast<std::size_t>(
                tp.tri[static_cast<std::size_t>(e)][static_cast<std::size_t>(p)])]
                .push_back(e);
    std::vector<sv_candidate> cand;
    for (int v = 0; v < tp.nv; ++v) {
        std::vector<sv_fan> fans =
            sv_vertex_fans(msh, tp, star[static_cast<std::size_t>(v)], v);
        bool any = false;
        for (std::size_t f = 0; f < fans.size(); ++f) {
            if (fans[f].classes > 2) continue;         // not singular
            any = true;
            sv_candidate c;
            c.kind = 0;
            c.entity = v;
            c.comp = static_cast<int>(f);
            c.t = 0;
            int sign = 1;
            for (std::size_t i = 0; i < fans[f].elems.size(); ++i) {
                int e = fans[f].elems[i];
                int p = sv_local_vertex_pos2(tp, e, v);
                multi_index<2> al;
                al.a[0] = 0; al.a[1] = 0; al.a[2] = 0;
                al.a[static_cast<std::size_t>(p)] = l;
                c.row.dof.push_back(e * Nl + im.rank(al));
                c.row.coef.push_back(sign);
                sign = -sign;
            }
            sv_normalize_row(c.row);
            cand.push_back(c);
        }
        if (any && singular_out) singular_out->push_back(v);
    }
    return cand;
}

// 3D: all edge-row candidates for broken P_l (l = n - 1)
template <typename T>
std::vector<sv_candidate> sv_edge_candidates(const mesh<3, T>& msh,
                                             const mesh_topology3& tp,
                                             int l,
                                             std::vector<int>* singular_out) {
    const index_map<3>& im = coeff_registry<3>::indices(l);
    const int Nl = im.size();
    // edge stars (through tet_edge: local edges 01,02,03,12,13,23)
    std::vector<std::vector<int> > star(
        static_cast<std::size_t>(tp.num_edges()));
    for (int e = 0; e < tp.nt; ++e)
        for (int le = 0; le < 6; ++le)
            star[static_cast<std::size_t>(
                tp.tet_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(le)])]
                .push_back(e);
    std::vector<sv_candidate> cand;
    for (int ed = 0; ed < tp.num_edges(); ++ed) {
        const int a = tp.edges[static_cast<std::size_t>(ed)][0];
        const int b = tp.edges[static_cast<std::size_t>(ed)][1];
        std::vector<sv_fan> fans =
            sv_edge_fans(msh, tp, star[static_cast<std::size_t>(ed)], ed);
        bool any = false;
        for (std::size_t f = 0; f < fans.size(); ++f) {
            if (fans[f].classes > 2) continue;         // not singular
            any = true;
            for (int t = 0; t <= l; ++t) {
                sv_candidate c;
                c.kind = 1;
                c.entity = ed;
                c.comp = static_cast<int>(f);
                c.t = t;
                int sign = 1;
                for (std::size_t i = 0; i < fans[f].elems.size(); ++i) {
                    int e = fans[f].elems[i];
                    int p, q;
                    sv_local_edge_pos3(tp.tet[static_cast<std::size_t>(e)],
                                       a, b, p, q);
                    multi_index<3> al;
                    al.a[0] = 0; al.a[1] = 0; al.a[2] = 0; al.a[3] = 0;
                    al.a[static_cast<std::size_t>(p)] = l - t;
                    al.a[static_cast<std::size_t>(q)] = t;
                    c.row.dof.push_back(e * Nl + im.rank(al));
                    c.row.coef.push_back(sign);
                    sign = -sign;
                }
                sv_normalize_row(c.row);
                cand.push_back(c);
            }
        }
        if (any && singular_out) singular_out->push_back(ed);
    }
    return cand;
}

// dimension dispatch for the public class
template <int D>
struct sv_backend;

template <>
struct sv_backend<2> {
    typedef mesh_topology2 topology_type;
    template <typename T>
    static std::vector<sv_candidate> candidates(const mesh<2, T>& msh, int l,
                                                std::vector<int>* sing_v,
                                                std::vector<int>*) {
        topology_type tp = topology_type::build(msh);
        return sv_vertex_candidates(msh, tp, l, sing_v);
    }
};

template <>
struct sv_backend<3> {
    typedef mesh_topology3 topology_type;
    template <typename T>
    static std::vector<sv_candidate> candidates(const mesh<3, T>& msh, int l,
                                                std::vector<int>*,
                                                std::vector<int>* sing_e) {
        topology_type tp = topology_type::build(msh);
        return sv_edge_candidates(msh, tp, l, sing_e);
    }
};

} // namespace detail

// ---------------------------------------------------------------------------
// sv_pressure_constraints<D, T> (V2)
// ---------------------------------------------------------------------------
template <int D, typename T>
class sv_pressure_constraints {
    static_assert(D == 2 || D == 3,
                  "bfem::sv_pressure_constraints: only D == 2 or D == 3");
public:
    // detection (geometry only) + generation (l = n - 1) + independization
    sv_pressure_constraints(const mesh<D, T>& msh, int n)
        : n_(n), rows_(), singular_vertices_(), singular_edges_() {
        if (n < 1)
            throw std::invalid_argument(
                "bfem::sv_pressure_constraints: n must be >= 1");
        std::vector<detail::sv_candidate> cand =
            detail::sv_backend<D>::candidates(msh, n - 1,
                                              &singular_vertices_,
                                              &singular_edges_);
        std::vector<detail::sv_candidate> sel = detail::sv_independize(cand);
        rows_.reserve(sel.size());
        for (std::size_t k = 0; k < sel.size(); ++k)
            rows_.push_back(sel[k].row);
    }

    int degree() const { return n_; }
    const std::vector<sv_constraint>& rows() const { return rows_; }

    // observation API
    int num_singular_vertices() const {
        return static_cast<int>(singular_vertices_.size());
    }
    int num_singular_edges() const {
        return static_cast<int>(singular_edges_.size());
    }
    const std::vector<int>& singular_vertices() const { return singular_vertices_; }
    const std::vector<int>& singular_edges() const { return singular_edges_; }

    // SV-1 / B-2: known sufficient stability conditions as an observation
    sv_stability_hint hint(sv_mesh_provenance p) const {
        sv_stability_hint h;
        if (p == sv_provenance_alfeld) {
            const int lo = (D == 2) ? 2 : 3;
            if (n_ >= lo) {
                h.status = sv_hint_stable_alfeld;
                h.reason = (D == 2)
                    ? "Alfeld (barycentric) split, n >= 2 (2D known sufficient condition)"
                    : "Alfeld (barycentric) split, n >= 3 (3D known sufficient condition)";
                return h;
            }
        }
        const int hi = (D == 2) ? 4 : 6;
        if (n_ >= hi) {
            h.status = sv_hint_stable_general;
            h.reason = (D == 2)
                ? "general mesh, n >= 4 (2D known sufficient condition)"
                : "general mesh, n >= 6 (3D known sufficient condition)";
            return h;
        }
        h.status = sv_hint_unknown;
        h.reason = "outside the recorded sufficient conditions (observation only)";
        return h;
    }

private:
    int n_;
    std::vector<sv_constraint> rows_;
    std::vector<int> singular_vertices_;   // 2D (empty in 3D)
    std::vector<int> singular_edges_;      // 3D (empty in 2D)
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_SV_ROWS_HPP
