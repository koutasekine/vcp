// vcp/bfem/sv/alfeld.hpp
// Phase 5d (Scott-Vogelius parts): alfeld_refine (V4, SV-ii / SV-5) --
// barycentric (Alfeld) split of every element into D + 1 children.
//
// Conforms to: SV external design v0.2 (section 6) and
//              SV internal design v0.2 (section 7).
//
// Normative conventions (SV-5):
//  - barycenter b_e = sum_i v_i * rational_to<T>(1, D+1): multiplications and
//    additions only, so interval T yields a rigorous enclosure of the exact
//    barycenter (the enclose-once constant carries the single division of the
//    conversion; no division appears in this file);
//  - child order is canonical: parent-element major, child k = 0..D in the
//    order of the parent's faces (2D: edges), where child k keeps the parent
//    vertex tuple with position k replaced by the barycenter ("face k as
//    base, barycenter as apex"). The replacement expresses the barycenter as
//    an affine combination whose coefficient on v_k is 1/(D+1) > 0, so every
//    child inherits the parent orientation and |child| = |parent| * 1/(D+1);
//  - vertex numbering is preserved: existing vertices keep their ids, the new
//    barycenters are appended as one trailing block in element order
//    (barycenter of parent e = old_num_vertices + e).
//
// Normative note (v0.2, A-3): the interior entities created by an Alfeld
// split are NOT singular in the SV sense -- a 2D barycenter has 3 edge
// direction lines, and a 3D interior edge (barycenter, vertex) has 3 face
// planes. The stabilization mechanism of the Alfeld split is not the creation
// of singular entities. (Checked by RSV-1 on the detection side.)
//
// This header is new phase-5d code; no frozen file is touched.

#ifndef VCP_BFEM_SV_ALFELD_HPP
#define VCP_BFEM_SV_ALFELD_HPP

#include <vector>
#include <array>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/convert_traits.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T>
mesh<D, T> alfeld_refine(const mesh<D, T>& msh) {
    static_assert(D == 2 || D == 3, "bfem::alfeld_refine: only D == 2 or D == 3");
    const int nv = msh.num_vertices();
    const int nt = msh.num_elements();
    const T w = rational_to<T>(1, D + 1);          // enclose-once constant

    std::vector<std::array<T, D> > verts;
    verts.reserve(static_cast<std::size_t>(nv + nt));
    for (int v = 0; v < nv; ++v) verts.push_back(msh.vertex(v));

    // barycenters, appended in element order (id = nv + e)
    for (int e = 0; e < nt; ++e) {
        const std::array<int, D + 1>& el = msh.element(e);
        std::array<T, D> b;
        for (int d = 0; d < D; ++d) {
            T acc = msh.vertex(el[0])[static_cast<std::size_t>(d)] * w;
            for (int i = 1; i <= D; ++i)
                acc += msh.vertex(el[static_cast<std::size_t>(i)])
                           [static_cast<std::size_t>(d)] * w;
            b[static_cast<std::size_t>(d)] = acc;
        }
        verts.push_back(b);
    }

    // children: parent major, child k = parent with vertex k -> barycenter
    std::vector<std::array<int, D + 1> > elems;
    elems.reserve(static_cast<std::size_t>(nt) * static_cast<std::size_t>(D + 1));
    for (int e = 0; e < nt; ++e) {
        const std::array<int, D + 1>& el = msh.element(e);
        for (int k = 0; k <= D; ++k) {
            std::array<int, D + 1> child = el;
            child[static_cast<std::size_t>(k)] = nv + e;
            elems.push_back(child);
        }
    }

    return mesh<D, T>::from_lists(verts, elems);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_ALFELD_HPP
