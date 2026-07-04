// vcp/bfem/c1/c1_dofmap.hpp
// Phase 6 (2D C1 Argyris family): global numbering and the signed
// local-to-global map (identity + sign; C1-4 -- no matrix-weighted scatter).
//
// Conforms to: C1 external design v0.2 (section 2, "global identification
// sign rule") and C1 internal design v0.2 (section 5).
//
// Global numbering (normative):
//   vertex block: vertex v owns [6v, 6v+6) in the component order
//     (value, du/dx, du/dy, dxx, dxy, dyy) -- physical Cartesian components,
//     frame independent (sign +1 always);
//   edge block: edge ed (canonical direction: global vertex small -> large)
//     owns (2m-9) dofs: trace point values (parameter (i+1)/(m-4), i =
//     0..m-6) then normal-derivative point values against nu_global =
//     R_{-90}(canonical edge vector) (parameter (j+1)/(m-3), j = 0..m-5);
//   interior block: element order, canonical interior lattice order.
//
// local -> global (C1-4): vertex and interior identity (+1); edge trace:
// index reversal only when the local direction is not canonical (+1); edge
// normal derivative: index reversal AND sign = (fwd ? +1 : -1) (nu_local =
// R_{-90}(local edge vector) = +-nu_global, the RT Q1 rule transplanted).
// Index reversal is folded into l2g; signs live in a separate +-1 array;
// the map drives the FROZEN Y2 general (signed) kernels via
// general_family_tag (S-C1-5: no new scatter skeleton exists).

#ifndef VCP_BFEM_C1_C1_DOFMAP_HPP
#define VCP_BFEM_C1_C1_DOFMAP_HPP

#include <vector>
#include <array>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/dofmap.hpp>           // mesh_topology2 + Y2 kernels
#include <vcp/bfem/c1/c1_tables.hpp>     // local dof index helpers

namespace vcp {
namespace bfem {

class c1_dofmap {
public:
    typedef detail::general_family_tag family_tag;

    c1_dofmap() : m_(0), nv_(0), ne_(0), nt_(0), nloc_(0), nint_(0), ndof_(0) {}

    static c1_dofmap build(const detail::mesh_topology2& tp, int m) {
        detail::c1_check_k(m, "c1_dofmap");
        c1_dofmap dm;
        dm.m_ = m;
        dm.nv_ = tp.nv;
        dm.ne_ = tp.num_edges();
        dm.nt_ = tp.nt;
        dm.nloc_ = detail::c1_dim(m);
        dm.nint_ = detail::c1_nint(m);
        dm.ndof_ = 6 * dm.nv_ + detail::c1_nedge(m) * dm.ne_ + dm.nint_ * dm.nt_;
        dm.edge_verts_ = tp.edges;
        dm.boundary_edges_ = tp.boundary_edges;
        dm.l2g_.assign(static_cast<std::size_t>(dm.nt_)
                       * static_cast<std::size_t>(dm.nloc_), -1);
        dm.sgn_.assign(dm.l2g_.size(), static_cast<signed char>(1));
        const int ntr = detail::c1_ntrace(m);
        const int nnd = detail::c1_nnd(m);
        for (int e = 0; e < dm.nt_; ++e) {
            const std::size_t base = static_cast<std::size_t>(e)
                                     * static_cast<std::size_t>(dm.nloc_);
            // vertex blocks: identity on the 6 Cartesian components
            for (int p = 0; p < 3; ++p) {
                int v = tp.tri[static_cast<std::size_t>(e)][static_cast<std::size_t>(p)];
                for (int c = 0; c < 6; ++c)
                    dm.l2g_[base + static_cast<std::size_t>(
                        detail::c1_vertex_dof(p, c))] = dm.vertex_dof(v, c);
            }
            // edge blocks
            for (int s = 0; s < 3; ++s) {
                int ed = tp.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(s)];
                bool fwd = tp.tri_edge_fwd[static_cast<std::size_t>(e)]
                                          [static_cast<std::size_t>(s)];
                for (int i = 0; i < ntr; ++i) {
                    std::size_t idx = base + static_cast<std::size_t>(
                        detail::c1_edge_trace_dof(m, s, i));
                    dm.l2g_[idx] = dm.edge_trace_dof(ed, fwd ? i : ntr - 1 - i);
                }
                for (int j = 0; j < nnd; ++j) {
                    std::size_t idx = base + static_cast<std::size_t>(
                        detail::c1_edge_nd_dof(m, s, j));
                    dm.l2g_[idx] = dm.edge_nd_dof(ed, fwd ? j : nnd - 1 - j);
                    dm.sgn_[idx] = static_cast<signed char>(fwd ? 1 : -1);
                }
            }
            // interior block: element order, canonical lattice order
            for (int r = 0; r < dm.nint_; ++r)
                dm.l2g_[base + static_cast<std::size_t>(
                    detail::c1_interior_dof(m, r))] = dm.interior_dof(e, r);
        }
        return dm;
    }

    int degree() const { return m_; }
    int ndof() const { return ndof_; }
    int local_size() const { return nloc_; }
    int num_elements() const { return nt_; }
    int num_vertices() const { return nv_; }
    int num_edges() const { return ne_; }

    int global_dof(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(r)];
    }
    int dof_sign(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return sgn_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(r)];
    }

    // ---- global numbering accessors (boundary generation and tests) ----
    int vertex_dof(int v, int c) const {
        assert(v >= 0 && v < nv_ && c >= 0 && c < 6);
        return 6 * v + c;
    }
    int edge_trace_dof(int ed, int i) const {
        assert(ed >= 0 && ed < ne_ && i >= 0 && i < detail::c1_ntrace(m_));
        return 6 * nv_ + ed * detail::c1_nedge(m_) + i;
    }
    int edge_nd_dof(int ed, int j) const {
        assert(ed >= 0 && ed < ne_ && j >= 0 && j < detail::c1_nnd(m_));
        return 6 * nv_ + ed * detail::c1_nedge(m_) + detail::c1_ntrace(m_) + j;
    }
    int interior_dof(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nint_);
        return 6 * nv_ + detail::c1_nedge(m_) * ne_ + e * nint_ + r;
    }

    const std::vector<int>& boundary_edge_ids() const { return boundary_edges_; }
    const std::array<int, 2>& edge_verts(int ed) const {
        assert(ed >= 0 && ed < ne_);
        return edge_verts_[static_cast<std::size_t>(ed)];
    }

private:
    int m_, nv_, ne_, nt_, nloc_, nint_, ndof_;
    std::vector<int> l2g_;
    std::vector<signed char> sgn_;
    std::vector<std::array<int, 2> > edge_verts_;
    std::vector<int> boundary_edges_;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_C1_C1_DOFMAP_HPP
