// vcp/bfem/cr1/cr1_dofmap.hpp
// CR1-1 (P1 Crouzeix-Raviart nonconforming element, D = 2, 3): the global
// numbering layer -- one degree of freedom per facet (D = 2: edge,
// D = 3: face).
//
// Conforms to: CR1-1 design v1.0 (sections 2, 3.1) and
//              CR1-1 implementation directive v1.0 (phase P1).
//
// Numbering (normative): the facet id IS the id the frozen topology layer
// assigns -- detail::mesh_topology2::edges (vcp/bfem/dofmap.hpp) for D = 2 and
// detail::mesh_topology3::faces (vcp/bfem/d3/topology3.hpp) for D = 3.  Both
// tables store the facet as its ASCENDING global vertex tuple and order the
// facets lexicographically by that tuple, which is exactly the convention the
// CR1 design prescribes (design 2, "vertex pair small -> large, lexicographic
// order").  No second enumeration exists here; the rt/ headers are NOT
// included (the convention is shared, the code is not).
//
// Local numbering (normative, identical for D = 2 and D = 3): local index i
// denotes the facet OPPOSITE local vertex i.  This is precisely the local
// facet index of both topology tables (mesh_topology2::tri_edge[e][k] is the
// edge opposite local vertex k; mesh_topology3::tet_face[e][k] is the face
// opposite local vertex k), so l2g is a table lookup with no re-derivation.
//
// Family tag (design 2, "orientation"): the CR1 degree of freedom is the
// AVERAGE of a scalar over a facet.  An average is invariant under any
// re-parameterization of the facet, so the two elements sharing a facet see
// one and the same functional and the local-to-global weight is identically
// +1.  The dofmap therefore declares detail::pn_family_tag and the frozen Y2
// identity scatter/gather kernels of vcp/bfem/dofmap.hpp apply verbatim
// (they multiply by no sign at all).
//
// Boundary degrees of freedom (design 2, "boundary condition"): the facets
// carried by exactly one element, i.e. the boundary_edges / boundary_faces
// list of the topology, already ascending and duplicate free.

#ifndef VCP_BFEM_CR1_CR1_DOFMAP_HPP
#define VCP_BFEM_CR1_CR1_DOFMAP_HPP

#include <vector>
#include <array>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dofmap.hpp>         // detail::mesh_topology2, detail::pn_family_tag
#include <vcp/bfem/d3/topology3.hpp>   // detail::mesh_topology3

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// cr1_topology<D>: the ONLY dimension dispatch of the CR1 numbering layer.
// It renames the facet surface of the two frozen topology structs; it holds
// no state and derives nothing.
// ---------------------------------------------------------------------------
template <int D>
struct cr1_topology;

template <>
struct cr1_topology<2> {
    typedef mesh_topology2 topology_type;
    static int num_facets(const topology_type& tp) { return tp.num_edges(); }
    // facet opposite local vertex i of element e
    static int facet_of(const topology_type& tp, int e, int i) {
        return tp.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(i)];
    }
    static const std::vector<std::array<int, 2> >& facet_vertices(const topology_type& tp) {
        return tp.edges;
    }
    static const std::vector<int>& boundary_facets(const topology_type& tp) {
        return tp.boundary_edges;
    }
};

template <>
struct cr1_topology<3> {
    typedef mesh_topology3 topology_type;
    static int num_facets(const topology_type& tp) { return tp.num_faces(); }
    static int facet_of(const topology_type& tp, int e, int i) {
        return tp.tet_face[static_cast<std::size_t>(e)][static_cast<std::size_t>(i)];
    }
    static const std::vector<std::array<int, 3> >& facet_vertices(const topology_type& tp) {
        return tp.faces;
    }
    static const std::vector<int>& boundary_facets(const topology_type& tp) {
        return tp.boundary_faces;
    }
};

// ---------------------------------------------------------------------------
// cr1_dofmap<D>: ndof / l2g / boundary_dofs.
//
// The Y2 kernel contract of vcp/bfem/dofmap.hpp asks a dofmap for
// global_dof(e, r), dof_sign(e, r) and local_size(); the CR1 design names the
// map l2g(e, i).  Both names are published and global_dof forwards to l2g, so
// there is one implementation and the frozen kernels bind without an adapter.
// ---------------------------------------------------------------------------
template <int D>
class cr1_dofmap {
    static_assert(D == 2 || D == 3,
                  "bfem::detail::cr1_dofmap: only D == 2 or D == 3");
public:
    typedef pn_family_tag family_tag;
    typedef typename cr1_topology<D>::topology_type topology_type;

    cr1_dofmap() : ndof_(0), nt_(0), l2g_(), facets_(), bdofs_() {}

    static cr1_dofmap from_topology(const topology_type& tp) {
        cr1_dofmap dm;
        dm.nt_ = tp.nt;
        dm.ndof_ = cr1_topology<D>::num_facets(tp);
        dm.facets_ = cr1_topology<D>::facet_vertices(tp);
        dm.bdofs_ = cr1_topology<D>::boundary_facets(tp);
        dm.l2g_.assign(static_cast<std::size_t>(dm.nt_)
                       * static_cast<std::size_t>(D + 1), -1);
        for (int e = 0; e < dm.nt_; ++e)
            for (int i = 0; i <= D; ++i)
                dm.l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(D + 1)
                        + static_cast<std::size_t>(i)] = cr1_topology<D>::facet_of(tp, e, i);
        return dm;
    }

    template <typename T>
    static cr1_dofmap from_mesh(const mesh<D, T>& m) {
        return from_topology(topology_type::build(m));
    }

    int ndof() const { return ndof_; }
    int local_size() const { return D + 1; }         // one dof per facet
    int num_elements() const { return nt_; }

    // element e, local index i (the facet opposite local vertex i) -> global
    int l2g(int e, int i) const {
        assert(e >= 0 && e < nt_);
        assert(i >= 0 && i <= D);
        return l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(D + 1)
                    + static_cast<std::size_t>(i)];
    }
    // Y2 kernel contract (identical to l2g; see the note above)
    int global_dof(int e, int i) const { return l2g(e, i); }
    int dof_sign(int e, int i) const {               // identically +1
        assert(e >= 0 && e < nt_);
        assert(i >= 0 && i <= D);
        return 1;
    }

    // ascending, duplicate free (inherited from the topology tables)
    std::vector<int> boundary_dofs() const { return bdofs_; }

    // the facet as its ASCENDING global vertex tuple: THE data of the facet
    // parameterization convention (see vcp/bfem/cr1/cr1_element_op.hpp)
    const std::array<int, D>& facet_vertices(int f) const {
        if (f < 0 || f >= ndof_)
            throw std::invalid_argument(
                "bfem::detail::cr1_dofmap::facet_vertices: facet id out of range");
        return facets_[static_cast<std::size_t>(f)];
    }

private:
    int ndof_, nt_;
    std::vector<int> l2g_;                     // nt x (D + 1)
    std::vector<std::array<int, D> > facets_;  // facet id -> ascending global tuple
    std::vector<int> bdofs_;
};

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CR1_CR1_DOFMAP_HPP
