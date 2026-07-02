// vcp/bfem/d3/dofmap3.hpp
// Phase 5a (3D common infrastructure): dofmap<3> (F3D-4) -- global numbering
// for the tetrahedral P^m family, local-to-global with edge and FACE
// identification, boundary dof enumeration.
//
// Conforms to: 3D common external design v0.1 (sections 2, 4) and
//              3D common internal design v0.1 (section 4).
//
// Numbering (normative):
//   vertex block [0, nv)
//   -> edge block (edge id order, canonical direction inside each edge:
//      t' counts toward the LARGER global vertex, flip rule t' = m - t)
//   -> face block (face id order; inside a face the canonical numbering of
//      external design 2-3: face multi-index (a,b,c), a,b,c >= 1, in the
//      canonical ascending-global vertex order, ranked by
//      index_map<2>(m-3).rank((a,b,c) - (1,1,1)))
//   -> interior block (element order, canonical L0 order of the interior
//      indices).
// The orientation permutation sigma is folded into l2g at build time
// (internal design 10.1); dof_sign is identically +1 for P^n (Y2).
//
// dofmap<3> is an explicit specialization placed beside the frozen 2D
// primary template (vcp/bfem/dofmap.hpp declares dofmap<D> with a D == 2
// static_assert in the primary; the primary is never instantiated for D == 3,
// so no 2D file changes). Construction only via detail::dofmap_builder3
// (the A-2 discipline: fe_space<3> will be the public entry in phase 5b;
// until then the builder is detail + test support, like RT-L0 trace_index).

#ifndef VCP_BFEM_D3_DOFMAP3_HPP
#define VCP_BFEM_D3_DOFMAP3_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/d3/topology3.hpp>

namespace vcp {
namespace bfem {

namespace detail {
struct dofmap_builder3;   // fills dofmap<3> from a topology (A-2: no public ctor)
} // namespace detail

// ---------------------------------------------------------------------------
// dofmap<3> (public API identical to dofmap<2>; boundary entities are FACES)
// ---------------------------------------------------------------------------
template <>
class dofmap<3> {
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

    // all boundary dofs: for every boundary face, the face-interior dofs plus
    // the dofs of its 3 edges and 3 vertices
    std::vector<int> boundary_dofs() const {
        return boundary_dofs(boundary_faces_);
    }
    // partial boundary (external design section 4, the 2D B-3 rule lifted to
    // faces): for each given FACE, include its face-interior dofs, the edge
    // dofs of its 3 edges AND the vertex dofs of its 3 vertices (edit the
    // returned list to drop entities if needed)
    std::vector<int> boundary_dofs(const std::vector<int>& faces) const {
        std::vector<int> out;
        std::vector<char> seen(static_cast<std::size_t>(ndof_), 0);
        int nfl = (m_ - 1) * (m_ - 2) / 2;                // face-interior per face
        int edge_base = nv_;
        int face_base = nv_ + ne_ * (m_ - 1);
        for (std::size_t s = 0; s < faces.size(); ++s) {
            int f = faces[s];
            if (f < 0 || f >= nf_)
                throw std::invalid_argument(
                    "bfem::dofmap<3>::boundary_dofs: face id out of range");
            for (int i = 0; i < 3; ++i) {
                int v = face_verts_[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)];
                if (!seen[static_cast<std::size_t>(v)]) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    out.push_back(v);
                }
            }
            for (int i = 0; i < 3; ++i) {
                int ed = face_edges_[static_cast<std::size_t>(f)][static_cast<std::size_t>(i)];
                for (int t = 1; t <= m_ - 1; ++t) {
                    int g = edge_base + ed * (m_ - 1) + (t - 1);
                    if (!seen[static_cast<std::size_t>(g)]) {
                        seen[static_cast<std::size_t>(g)] = 1;
                        out.push_back(g);
                    }
                }
            }
            for (int r = 0; r < nfl; ++r) {
                int g = face_base + f * nfl + r;
                if (!seen[static_cast<std::size_t>(g)]) {
                    seen[static_cast<std::size_t>(g)] = 1;
                    out.push_back(g);
                }
            }
        }
        std::sort(out.begin(), out.end());
        return out;
    }
    // boundary face ids (ascending)
    const std::vector<int>& boundary_face_ids() const { return boundary_faces_; }

private:
    dofmap() : m_(0), nv_(0), ne_(0), nf_(0), nt_(0), nloc_(0), ndof_(0) {}
    int m_, nv_, ne_, nf_, nt_, nloc_, ndof_;
    std::vector<int> l2g_;                            // nt x nloc
    std::vector<std::array<int, 2> > edge_verts_;     // edge id -> (vmin, vmax)
    std::vector<std::array<int, 3> > face_verts_;     // face id -> ascending triple
    std::vector<std::array<int, 3> > face_edges_;     // face id -> 3 edge ids
    std::vector<int> boundary_faces_;

    friend struct detail::dofmap_builder3;
};

namespace detail {

struct dofmap_builder3 {
    // numbering per 3D external design section 4 / internal design section 4
    static dofmap<3> build(const mesh_topology3& tp, int m) {
        if (m < 1)
            throw std::invalid_argument("bfem::dofmap<3>: degree must be >= 1");
        dofmap<3> dm;
        dm.m_ = m;
        dm.nv_ = tp.nv;
        dm.ne_ = tp.num_edges();
        dm.nf_ = tp.num_faces();
        dm.nt_ = tp.nt;
        dm.edge_verts_ = tp.edges;
        dm.face_verts_ = tp.faces;
        dm.face_edges_ = tp.face_edge;
        dm.boundary_faces_ = tp.boundary_faces;
        const index_map<3>& im = coeff_registry<3>::indices(m);
        dm.nloc_ = im.size();
        int nfl = (m - 1) * (m - 2) / 2;                       // per face
        int nil = (m - 1) * (m - 2) * (m - 3) / 6;             // per element
        dm.ndof_ = tp.nv + tp.num_edges() * (m - 1)
                 + tp.num_faces() * nfl + tp.nt * nil;
        dm.l2g_.assign(static_cast<std::size_t>(tp.nt)
                       * static_cast<std::size_t>(dm.nloc_), -1);
        int edge_base = tp.nv;
        int face_base = tp.nv + tp.num_edges() * (m - 1);
        int int_base = face_base + tp.num_faces() * nfl;
        // rank map of the face-interior indices (only needed when m >= 3)
        const index_map<2>* im_face =
            (m >= 3) ? &coeff_registry<2>::indices(m - 3) : 0;
        for (int e = 0; e < tp.nt; ++e) {
            int int_count = 0;   // interior indices in canonical rank order
            for (int r = 0; r < dm.nloc_; ++r) {
                multi_index<3> al = im.unrank(r);
                int zero_count = 0, zpos = -1, ppos = -1;
                int nz[2]; int nnz = 0;
                for (int i = 0; i < 4; ++i) {
                    if (al.a[static_cast<std::size_t>(i)] == 0) { ++zero_count; zpos = i; }
                    else { ppos = i; if (nnz < 2) nz[nnz] = i; ++nnz; }
                }
                int g;
                if (zero_count == 3) {
                    // vertex dof: alpha = m * e_p
                    g = tp.tet[static_cast<std::size_t>(e)][static_cast<std::size_t>(ppos)];
                } else if (zero_count == 2) {
                    // edge dof on local edge (p, q), p < q local; the in-edge
                    // position t is the exponent on the larger-local side q,
                    // flip rule t' = m - t (internal design section 4)
                    int p = nz[0], q = nz[1];
                    int le = tet_local::edge_of_pair(p, q);
                    int t = al.a[static_cast<std::size_t>(q)];
                    int ed = tp.tet_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(le)];
                    int tprime = tp.tet_edge_fwd[static_cast<std::size_t>(e)]
                                                [static_cast<std::size_t>(le)]
                                 ? t : m - t;
                    g = edge_base + ed * (m - 1) + (tprime - 1);
                } else if (zero_count == 1) {
                    // face dof on local face k = zpos: canonical face
                    // multi-index via THE shared face convention
                    // (face_canonical_beta, S-5A-2), then the L0 canonical
                    // rank of (a,b,c) - (1,1,1) in degree m - 3
                    int k = zpos;
                    std::array<int, 3> bc = face_canonical_beta(
                        al, k,
                        tp.tet_face_perm[static_cast<std::size_t>(e)]
                                        [static_cast<std::size_t>(k)]);
                    multi_index<2> b2;
                    b2.a[0] = bc[0] - 1;
                    b2.a[1] = bc[1] - 1;
                    b2.a[2] = bc[2] - 1;
                    int rank2 = im_face->rank(b2);
                    int f = tp.tet_face[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)];
                    g = face_base + f * nfl + rank2;
                } else {
                    // interior dof: element block, canonical order
                    g = int_base + e * nil + int_count;
                    ++int_count;
                }
                dm.l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(dm.nloc_)
                        + static_cast<std::size_t>(r)] = g;
            }
        }
        return dm;
    }
};

} // namespace detail

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_D3_DOFMAP3_HPP
