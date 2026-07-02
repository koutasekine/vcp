// vcp/bfem/d3/trace3.hpp
// Phase 5a (3D common infrastructure): trace index maps (F3D-5) --
// face trace (3D coefficients -> canonically oriented 2D coefficients) and
// edge trace (the 2D trace_index mechanism lifted to the tetrahedron).
//
// Conforms to: 3D common external design v0.1 (section 5) and
//              3D common internal design v0.1 (section 5).
//
// Single implementation point (S-5A-2, normative): the face convention
// (beta construction + sigma application) is face_canonical_beta() in
// topology3.hpp. face_trace_index below iterates the 3D ranks with
// alpha_k == 0 and calls THAT function to place each coefficient at its 2D
// canonical rank; it does not re-implement the convention. dof_build3
// (dofmap3.hpp) calls the same function for the face-interior block.
//
// detail + test support (the RT-L0 trace_index treatment): consumed by the
// continuity gate ET-1 here, by the 3D RT face-moment dofs in phase 5c and
// by the boundary handling in 5d.

#ifndef VCP_BFEM_D3_TRACE3_HPP
#define VCP_BFEM_D3_TRACE3_HPP

#include <vector>
#include <array>
#include <cassert>

#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/d3/topology3.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// face trace: element e, local face k, degree m -> index list `out` of
// length N(2, m) with out[rank2] = rank3, where rank2 runs over the L0
// canonical 2D order of the face polynomial written in the CANONICAL
// (ascending global) vertex order of the face, and rank3 is the 3D rank of
// the coefficient that restricts to it (alpha_k == 0 slice, sigma applied).
// ---------------------------------------------------------------------------
inline std::vector<int> face_trace_index(const mesh_topology3& tp,
                                         int e, int k, int m) {
    assert(e >= 0 && e < tp.nt && k >= 0 && k < 4 && m >= 0);
    const index_map<3>& im3 = coeff_registry<3>::indices(m);
    const index_map<2>& im2 = coeff_registry<2>::indices(m);
    std::vector<int> out(static_cast<std::size_t>(im2.size()), -1);
    const int code = tp.tet_face_perm[static_cast<std::size_t>(e)]
                                     [static_cast<std::size_t>(k)];
    for (int r = 0; r < im3.size(); ++r) {
        multi_index<3> al = im3.unrank(r);
        if (al.a[static_cast<std::size_t>(k)] != 0) continue;   // off the face
        // THE shared face convention (S-5A-2): no second implementation here
        std::array<int, 3> bc = face_canonical_beta(al, k, code);
        multi_index<2> b2;
        b2.a[0] = bc[0];
        b2.a[1] = bc[1];
        b2.a[2] = bc[2];
        out[static_cast<std::size_t>(im2.rank(b2))] = r;
    }
    return out;
}

// ---------------------------------------------------------------------------
// edge trace: local edge le (local pair (p, q), p < q, order 01,02,03,12,
// 13,23), degree m -> index list of length m + 1 with out[t] = the 3D rank
// of alpha = (m - t) e_p + t e_q (t counts toward the larger LOCAL vertex q;
// the global flip t' = m - t is the dofmap's business, not the trace's --
// same division of labor as the 2D trace_index).
// ---------------------------------------------------------------------------
inline std::vector<int> edge_trace_index3(int le, int m) {
    assert(le >= 0 && le < 6 && m >= 0);
    const index_map<3>& im3 = coeff_registry<3>::indices(m);
    const int p = tet_local::edge_vertex(le, 0);
    const int q = tet_local::edge_vertex(le, 1);
    std::vector<int> out(static_cast<std::size_t>(m + 1), -1);
    for (int t = 0; t <= m; ++t) {
        multi_index<3> al;
        al.a[0] = al.a[1] = al.a[2] = al.a[3] = 0;
        al.a[static_cast<std::size_t>(p)] = m - t;
        al.a[static_cast<std::size_t>(q)] = t;
        out[static_cast<std::size_t>(t)] = im3.rank(al);
    }
    return out;
}

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_D3_TRACE3_HPP
