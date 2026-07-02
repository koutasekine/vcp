// vcp/bfem/dofmap.hpp
// Layer 3: mesh topology (detail) and dofmap<D> (H2) -- global numbering,
// local-to-global map with edge identification, boundary dof enumeration.
//
// Conforms to: L3 external design v0.2 (section 3.2) and
//              L3 internal design v0.2 (sections 2, 3).
//
// Numbering (D = 2, normative):
//   vertex block [0, nv) -> edge block (edge id order, canonical direction
//   inside each edge) -> interior block (element order, canonical L0 order of
//   the interior indices).
// Canonical edge: global vertex (small -> large); the local-to-global map
// absorbs per-element orientation flips (the continuity implementation point,
// verified exactly by E1/K1).
//
// Local edge direction convention (B-2, normative): local edge k (opposite
// local vertex k) has local positive direction "local vertex (k+1)%3 ->
// (k+2)%3"; tri_edge_fwd[e][k] = ( global((k+1)%3) < global((k+2)%3) );
// the in-edge position t is the exponent on the terminal side ((k+2)%3), and
// the flip rule is t' = m - t.
//
// Weighted scatter/gather (Y2): the general kernels take signs
// ("global += sign_i * sign_j * local"); for the P^n family the signs are
// identically +1 and assembly dispatches at compile time to the identity
// kernels that contain no sign multiplication (structural check S-L3-1).

#ifndef VCP_BFEM_DOFMAP_HPP
#define VCP_BFEM_DOFMAP_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class fe_space;
template <int D> class dofmap;

namespace detail {

// ---------------------------------------------------------------------------
// topology (D = 2): edge table and boundary detection; coordinates unused.
// ---------------------------------------------------------------------------
struct mesh_topology2 {
    int nv, nt;
    std::vector<std::array<int, 3> > tri;
    std::vector<std::array<int, 2> > edges;          // edge id -> (v_min, v_max), lex sorted
    std::vector<std::array<int, 3> > tri_edge;       // element -> edge ids (local edge k)
    std::vector<std::array<bool, 3> > tri_edge_fwd;  // local direction == canonical?
    std::vector<int> boundary_edges;                 // edges with exactly 1 element (sorted)

    template <typename T>
    static mesh_topology2 build(const mesh<2, T>& m) {
        mesh_topology2 tp;
        tp.nv = m.num_vertices();
        tp.nt = m.num_elements();
        tp.tri.resize(static_cast<std::size_t>(tp.nt));
        for (int e = 0; e < tp.nt; ++e) tp.tri[static_cast<std::size_t>(e)] = m.element(e);
        // enumerate all local edges as (vmin, vmax, e, k), sort, unify
        struct rec { int vmin, vmax, e, k; };
        std::vector<rec> rs;
        rs.reserve(static_cast<std::size_t>(3 * tp.nt));
        for (int e = 0; e < tp.nt; ++e) {
            const std::array<int, 3>& t = tp.tri[static_cast<std::size_t>(e)];
            for (int k = 0; k < 3; ++k) {
                int a = t[static_cast<std::size_t>((k + 1) % 3)];
                int b = t[static_cast<std::size_t>((k + 2) % 3)];
                if (a == b)
                    throw std::invalid_argument("bfem::topology: element with repeated vertex");
                rec r;
                r.vmin = a < b ? a : b;
                r.vmax = a < b ? b : a;
                r.e = e;
                r.k = k;
                rs.push_back(r);
            }
        }
        std::sort(rs.begin(), rs.end(), rec_less_t());
        tp.tri_edge.resize(static_cast<std::size_t>(tp.nt));
        tp.tri_edge_fwd.resize(static_cast<std::size_t>(tp.nt));
        int ne = 0;
        std::size_t i = 0;
        while (i < rs.size()) {
            std::size_t j = i;
            while (j < rs.size() && rs[j].vmin == rs[i].vmin && rs[j].vmax == rs[i].vmax) ++j;
            std::size_t count = j - i;
            if (count > 2)
                throw std::invalid_argument(
                    "bfem::topology: non-manifold input (an edge shared by 3+ elements)");
            std::array<int, 2> ev = { { rs[i].vmin, rs[i].vmax } };
            tp.edges.push_back(ev);
            if (count == 1) tp.boundary_edges.push_back(ne);
            for (std::size_t s = i; s < j; ++s) {
                int e = rs[s].e, k = rs[s].k;
                tp.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)] = ne;
                const std::array<int, 3>& t = tp.tri[static_cast<std::size_t>(e)];
                // local positive direction: (k+1)%3 -> (k+2)%3 (B-2)
                tp.tri_edge_fwd[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)] =
                    (t[static_cast<std::size_t>((k + 1) % 3)]
                     < t[static_cast<std::size_t>((k + 2) % 3)]);
            }
            ++ne;
            i = j;
        }
        // boundary_edges built in ascending edge id order (already sorted)
        return tp;
    }

    int num_edges() const { return static_cast<int>(edges.size()); }

private:
    struct rec_less_t {
        template <typename R>
        bool operator()(const R& a, const R& b) const {
            if (a.vmin != b.vmin) return a.vmin < b.vmin;
            if (a.vmax != b.vmax) return a.vmax < b.vmax;
            if (a.e != b.e) return a.e < b.e;       // determinism inside a key
            return a.k < b.k;
        }
    };
};

struct dofmap_builder;   // fills dofmap<2> from a topology (A-2: no public ctor)

// ---------------------------------------------------------------------------
// Y2 scatter/gather kernels. General (signed) kernels accept per-dof signs;
// the P^n identity kernels contain no sign multiplication and are selected at
// compile time through the family tag of dofmap (pn_family_tag).
// ---------------------------------------------------------------------------
struct pn_family_tag {};        // signs identically +1, identity kernels
struct general_family_tag {};   // signed kernels (RT/Nedelec families later)

} // namespace detail

// ---------------------------------------------------------------------------
// dofmap<D> (public API; construction only via fe_space::dofs(m), A-2)
// ---------------------------------------------------------------------------
template <int D>
class dofmap {
    static_assert(D == 2, "bfem::dofmap: initial version supports D == 2 only");
public:
    typedef detail::pn_family_tag family_tag;   // Y2 compile-time family

    int ndof() const { return ndof_; }
    int degree() const { return m_; }

    // element e, local rank (L0 canonical order) -> global dof, O(1)
    int global_dof(int e, int local_rank) const {
        assert(e >= 0 && e < nt_);
        assert(local_rank >= 0 && local_rank < nloc_);
        return l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(local_rank)];
    }
    // local-to-global weight (sign). Identically +1 for P^n (Y2); the P^n
    // dofmap does NOT store a sign array.
    int dof_sign(int e, int local_rank) const {
        assert(e >= 0 && e < nt_);
        assert(local_rank >= 0 && local_rank < nloc_);
        return 1;
    }
    int local_size() const { return nloc_; }
    int num_elements() const { return nt_; }

    // all boundary dofs: boundary-edge interior dofs + their endpoint vertices
    std::vector<int> boundary_dofs() const {
        return boundary_dofs(boundary_edges_);
    }
    // partial boundary (B-3): the edge dofs on the given edges AND the
    // endpoint vertex dofs of those edges (edit the returned list to drop
    // endpoints if needed)
    std::vector<int> boundary_dofs(const std::vector<int>& edges) const {
        std::vector<int> out;
        std::vector<char> seen(static_cast<std::size_t>(ndof_), 0);
        for (std::size_t s = 0; s < edges.size(); ++s) {
            int ed = edges[s];
            if (ed < 0 || ed >= ne_)
                throw std::invalid_argument("bfem::dofmap::boundary_dofs: edge id out of range");
            for (int side = 0; side < 2; ++side) {
                int v = edge_verts_[static_cast<std::size_t>(ed)][static_cast<std::size_t>(side)];
                if (!seen[static_cast<std::size_t>(v)]) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    out.push_back(v);
                }
            }
            for (int t = 1; t <= m_ - 1; ++t) {
                int g = nv_ + ed * (m_ - 1) + (t - 1);
                if (!seen[static_cast<std::size_t>(g)]) {
                    seen[static_cast<std::size_t>(g)] = 1;
                    out.push_back(g);
                }
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }
    // boundary edge ids (ascending)
    const std::vector<int>& boundary_edge_ids() const { return boundary_edges_; }

private:
    dofmap() : m_(0), nv_(0), ne_(0), nt_(0), nloc_(0), ndof_(0) {}
    int m_, nv_, ne_, nt_, nloc_, ndof_;
    std::vector<int> l2g_;                            // nt x nloc
    std::vector<std::array<int, 2> > edge_verts_;     // edge id -> (vmin, vmax)
    std::vector<int> boundary_edges_;

    friend struct detail::dofmap_builder;
};

namespace detail {

struct dofmap_builder {
    // numbering per external design 3.2 / internal design 3
    static dofmap<2> build(const mesh_topology2& tp, int m) {
        if (m < 1)
            throw std::invalid_argument("bfem::dofmap: degree must be >= 1");
        dofmap<2> dm;
        dm.m_ = m;
        dm.nv_ = tp.nv;
        dm.ne_ = tp.num_edges();
        dm.nt_ = tp.nt;
        dm.edge_verts_ = tp.edges;
        dm.boundary_edges_ = tp.boundary_edges;
        const index_map<2>& im = coeff_registry<2>::indices(m);
        dm.nloc_ = im.size();
        int n_int = (m - 1) * (m - 2) / 2;
        dm.ndof_ = tp.nv + tp.num_edges() * (m - 1) + tp.nt * n_int;
        dm.l2g_.assign(static_cast<std::size_t>(tp.nt)
                       * static_cast<std::size_t>(dm.nloc_), -1);
        int int_base = tp.nv + tp.num_edges() * (m - 1);
        for (int e = 0; e < tp.nt; ++e) {
            int int_count = 0;   // interior indices in canonical rank order
            for (int r = 0; r < dm.nloc_; ++r) {
                multi_index<2> al = im.unrank(r);
                int zero_count = 0, zpos = -1, ppos = -1;
                for (int i = 0; i < 3; ++i) {
                    if (al.a[static_cast<std::size_t>(i)] == 0) { ++zero_count; zpos = i; }
                    else ppos = i;
                }
                int g;
                if (zero_count == 2) {
                    // vertex dof: alpha = m * e_p
                    g = tp.tri[static_cast<std::size_t>(e)][static_cast<std::size_t>(ppos)];
                } else if (zero_count == 1) {
                    // edge dof on local edge k = zpos; in-edge position t is
                    // the exponent on the terminal side (k+2)%3 (B-2)
                    int k = zpos;
                    int t = al.a[static_cast<std::size_t>((k + 2) % 3)];
                    int ed = tp.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)];
                    int tprime = tp.tri_edge_fwd[static_cast<std::size_t>(e)]
                                                [static_cast<std::size_t>(k)]
                                 ? t : m - t;
                    g = tp.nv + ed * (m - 1) + (tprime - 1);
                } else {
                    // interior dof: element block, canonical order
                    g = int_base + e * n_int + int_count;
                    ++int_count;
                }
                dm.l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(dm.nloc_)
                        + static_cast<std::size_t>(r)] = g;
            }
        }
        return dm;
    }
};

// ---------------------------------------------------------------------------
// Y2 kernels (matrix scatter, vector scatter, gather).
// General kernels: signed. Identity kernels: no sign multiplication at all.
// Assembly code calls the dispatch() overloads on dofmap<D>::family_tag.
// ---------------------------------------------------------------------------

// --- general (signed) matrix scatter: buf.push(gi, gj, si * sj * v) ---
template <typename T, typename Buf, typename Loc>
void scatter_matrix_general(const dofmap<2>& dm, int e, const Loc& loc, int n,
                            Buf& buf) {
    for (int a = 0; a < n; ++a) {
        int ga = dm.global_dof(e, a);
        T sa = T(dm.dof_sign(e, a));
        for (int b = 0; b < n; ++b) {
            int gb = dm.global_dof(e, b);
            T sb = T(dm.dof_sign(e, b));
            buf.push(ga, gb, sa * sb * loc(a, b));
        }
    }
}

// --- P^n identity matrix scatter: no sign multiplication (S-L3-1) ---
template <typename T, typename Buf, typename Loc>
void scatter_matrix_identity(const dofmap<2>& dm, int e, const Loc& loc, int n,
                             Buf& buf) {
    for (int a = 0; a < n; ++a) {
        int ga = dm.global_dof(e, a);
        for (int b = 0; b < n; ++b)
            buf.push(ga, dm.global_dof(e, b), loc(a, b));
    }
}

template <typename T, typename Buf, typename Loc>
void scatter_matrix(const dofmap<2>& dm, int e, const Loc& loc, int n, Buf& buf,
                    pn_family_tag) {
    scatter_matrix_identity<T>(dm, e, loc, n, buf);
}
template <typename T, typename Buf, typename Loc>
void scatter_matrix(const dofmap<2>& dm, int e, const Loc& loc, int n, Buf& buf,
                    general_family_tag) {
    scatter_matrix_general<T>(dm, e, loc, n, buf);
}

// --- gather: local = sign * global (identity: plain copy) ---
template <typename T, typename Vec>
void gather_general(const dofmap<2>& dm, int e, const Vec& g, T* loc, int n) {
    for (int r = 0; r < n; ++r)
        loc[r] = T(dm.dof_sign(e, r)) * g(dm.global_dof(e, r), 0);
}
template <typename T, typename Vec>
void gather_identity(const dofmap<2>& dm, int e, const Vec& g, T* loc, int n) {
    for (int r = 0; r < n; ++r)
        loc[r] = g(dm.global_dof(e, r), 0);
}
template <typename T, typename Vec>
void gather(const dofmap<2>& dm, int e, const Vec& g, T* loc, int n, pn_family_tag) {
    gather_identity(dm, e, g, loc, n);
}
template <typename T, typename Vec>
void gather(const dofmap<2>& dm, int e, const Vec& g, T* loc, int n, general_family_tag) {
    gather_general(dm, e, g, loc, n);
}

// --- vector scatter: global += sign * local ---
template <typename T, typename Vec>
void scatter_vector_general(const dofmap<2>& dm, int e, const T* loc, int n, Vec& g) {
    for (int r = 0; r < n; ++r)
        g(dm.global_dof(e, r), 0) += T(dm.dof_sign(e, r)) * loc[r];
}
template <typename T, typename Vec>
void scatter_vector_identity(const dofmap<2>& dm, int e, const T* loc, int n, Vec& g) {
    for (int r = 0; r < n; ++r)
        g(dm.global_dof(e, r), 0) += loc[r];
}
template <typename T, typename Vec>
void scatter_vector(const dofmap<2>& dm, int e, const T* loc, int n, Vec& g, pn_family_tag) {
    scatter_vector_identity(dm, e, loc, n, g);
}
template <typename T, typename Vec>
void scatter_vector(const dofmap<2>& dm, int e, const T* loc, int n, Vec& g, general_family_tag) {
    scatter_vector_general(dm, e, loc, n, g);
}

} // namespace detail

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_DOFMAP_HPP
