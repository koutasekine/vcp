// vcp/bfem/d3/bey_table.hpp
// Phase 5a (3D common infrastructure): Bey red refinement of the reference
// tetrahedron (F3D-7) -- the 8-child vertex table with the FIXED diagonal,
// and the D = 3 overload of range_refined.
//
// Conforms to: 3D common external design v0.1 (section 7),
//              3D common internal design v0.1 (section 6) and
//              L1 internal design v0.4 (section 11.4, W3D-5 correction).
//
// FIXED diagonal (normative, W3D-5): the interior octahedron (vertices
// m01, m02, m03, m12, m13, m23) is cut along the CONVENTIONAL diagonal
// m01 - m23, statically. There is deliberately NO "shortest diagonal"
// selection logic anywhere in this header: range_refined subdivides the
// REFERENCE simplex, where a physical-shape-based choice is meaningless
// (the v0.3 "shortest diagonal" rule was retracted; see L1 internal v0.4).
//
// Children (barycentric coordinates of the parent, entries in {0, 1/2, 1};
// h = 1/2 enters once through rational_to<T> -- enclose-once):
//   corner children  C0 = (v0,  m01, m02, m03)
//                    C1 = (m01, v1,  m12, m13)
//                    C2 = (m02, m12, v2,  m23)
//                    C3 = (m03, m13, m23, v3 )
//   octahedron split C4 = (m01, m02, m03, m23)   -- all four share the
//                    C5 = (m01, m03, m13, m23)      diagonal m01 - m23;
//                    C6 = (m01, m13, m12, m23)      (m02, m03, m13, m12) is
//                    C7 = (m01, m12, m02, m23)      the equatorial cycle
//
// refine.hpp (2D, frozen) is NOT modified: the refinement algorithm is
// restrict_to-based and D independent; only this vertex table is new, and
// the D = 3 range_refined below is a parallel overload of the D = 2 one.

#ifndef VCP_BFEM_D3_BEY_TABLE_HPP
#define VCP_BFEM_D3_BEY_TABLE_HPP

#include <array>
#include <stdexcept>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/refine.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// the 8 children x 4 vertices in barycentric coordinates of the parent
template <typename T>
std::array<std::array<bary_point<3, T>, 4>, 8> bey_children_3d() {
    const T o(1);
    const T z(0);
    const T h = ::vcp::bfem::rational_to<T>(1, 2);
    // barycentric rows for the parent vertices and the edge midpoints
    // (index 0..3 = v0..v3, then m01, m02, m03, m12, m13, m23)
    //                         v0 v1 v2 v3 m01 m02 m03 m12 m13 m23
    // encoded per child below via small index tables
    static const int child_vertex[8][4] = {
        { 0, 4, 5, 6 },   // C0 = (v0,  m01, m02, m03)
        { 4, 1, 7, 8 },   // C1 = (m01, v1,  m12, m13)
        { 5, 7, 2, 9 },   // C2 = (m02, m12, v2,  m23)
        { 6, 8, 9, 3 },   // C3 = (m03, m13, m23, v3 )
        { 4, 5, 6, 9 },   // C4 = (m01, m02, m03, m23)
        { 4, 6, 8, 9 },   // C5 = (m01, m03, m13, m23)
        { 4, 8, 7, 9 },   // C6 = (m01, m13, m12, m23)
        { 4, 7, 5, 9 }    // C7 = (m01, m12, m02, m23)
    };
    static const int point_bary[10][4] = {   // 2 means h, 1 means 1, 0 means 0
        { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 },
        { 2, 2, 0, 0 }, { 2, 0, 2, 0 }, { 2, 0, 0, 2 },
        { 0, 2, 2, 0 }, { 0, 2, 0, 2 }, { 0, 0, 2, 2 }
    };
    std::array<std::array<bary_point<3, T>, 4>, 8> kids;
    for (int c = 0; c < 8; ++c) {
        for (int i = 0; i < 4; ++i) {
            const int* row = point_bary[child_vertex[c][i]];
            for (int d = 0; d < 4; ++d) {
                const T& val = (row[d] == 0) ? z : ((row[d] == 1) ? o : h);
                kids[static_cast<std::size_t>(c)][static_cast<std::size_t>(i)]
                    [static_cast<std::size_t>(d)] = val;
            }
        }
    }
    return kids;
}

template <typename T>
void range_refined_rec_3d(
        const bpoly<3, T>& u, int d,
        const std::array<std::array<bary_point<3, T>, 4>, 8>& kids,
        range_pair<T>& acc, bool& started) {
    if (d == 0) {
        range_merge(acc, started, u);
        return;
    }
    for (int c = 0; c < 8; ++c)
        range_refined_rec_3d(restrict_to(u, kids[static_cast<std::size_t>(c)]),
                             d - 1, kids, acc, started);
}

} // namespace detail

// D = 3 overload of range_refined (F10 lifted by the Bey table; the D = 2
// overload in refine.hpp is untouched). depth == 0 equals range(u).
template <typename T>
range_pair<T> range_refined(const bpoly<3, T>& u, int depth) {
    if (depth < 0)
        throw std::invalid_argument("bfem::range_refined: negative depth");
    if (depth == 0) return range(u);
    std::array<std::array<bary_point<3, T>, 4>, 8> kids =
        detail::bey_children_3d<T>();
    range_pair<T> acc;
    bool started = false;
    detail::range_refined_rec_3d(u, depth, kids, acc, started);
    return acc;
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_D3_BEY_TABLE_HPP
