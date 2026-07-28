// vcp/bfem/meshgen3.hpp
// MG-1 Phase C: tetrahedral meshing of an extruded domain (2D cross-section
// x z-interval, common h -- ruling Q4), design section 5.4.
//
// Dependencies: meshgen.hpp (the base layer never includes this header).
// meshgen_delaunay.hpp is additionally included ONLY to provide the
// Delaunay-cross-section entry required by design section 6 / directive
// Phase 5; the base remains ignorant of both (recorded in the report).
//
// Size contract: every edge of the produced tetrahedra satisfies
// len^2 <= h^2. Since prism-splitting diagonals combine a cross-section
// edge with a layer rise, BOTH the cross-section refinement and the layer
// thickness use the squared threshold (h*h)/2, which keeps the stated
// section properties (cross-section edges^2 <= h^2, thickness^2 <= h^2)
// AND the all-edge contract of design section 3.3 (recorded in the report).
//
// facet_source encoding for D = 3 (recorded in the report):
//   side faces:   loop = 2D loop (0 = outer, 1.. = holes), segment = the 2D
//                 segment index (inherited from the cross-section boundary)
//   bottom z=z0:  loop = -1, segment = cross-section triangle index
//   top    z=z1:  loop = -2, segment = cross-section triangle index

#ifndef VCP_BFEM_MESHGEN3_HPP
#define VCP_BFEM_MESHGEN3_HPP

#include <vcp/bfem/meshgen.hpp>
#include <vcp/bfem/meshgen_delaunay.hpp>

namespace vcp {
namespace bfem {

template <typename T>
struct extruded_domain {                                   // D = 3 (Phase C)
    polygon_domain<T> base;   // cross-section (non-convex, holes allowed)
    T z0;                     // extrusion interval [z0, z1]
    T z1;
};

namespace meshgen_detail {

// sorted vertex triple (int comparisons only)
inline std::array<int, 3> face_key(int a, int b, int c) {
    int t;
    if (b < a) { t = a; a = b; b = t; }
    if (c < b) { t = b; b = c; c = t; }
    if (b < a) { t = a; a = b; b = t; }
    std::array<int, 3> k = {{a, b, c}};
    return k;
}

// extrude an already-built (and possibly Delaunay-flipped) 2D coarse mesh.
// Layer levels are dyadic bisections of [z0, z1] (midpoint form (a+b)/2,
// mirroring the 2D refinement); prisms split into 3 tetrahedra by the
// min-index rule, which is quad-consistent across neighbouring prisms.
template <typename T>
void extrude_and_finalize(coarse_mesh<T>& cm2, const T& z0, const T& z1,
                          const T& h, const meshgen_options& opt,
                          meshgen_status<3>& status,
                          std::vector<std::array<T, 3> >& vertices,
                          std::vector<std::array<int, 4> >& elements) {
    if (!certainly_pos(z1 - z0))
        throw meshgen_error(
            "vcp::bfem::meshgen: z1 must certainly exceed z0");

    const T hh2 = (h * h) / T(2);   // shared squared threshold (see header)

    // 1) refine the cross-section
    meshgen_status<2> st2;
    std::vector<std::array<T, 2> > v2;
    std::vector<std::array<int, 3> > e2;
    refine_and_finalize(cm2, hh2, opt, st2, v2, e2);
    const int nv2 = static_cast<int>(v2.size());

    // 2) layer thickness: halve until dz*dz < hh2 is certainly true
    T dz = z1 - z0;
    int zsteps = 0;
    while (!certainly_less(dz * dz, hh2)) {
        if (zsteps == opt.max_refine)
            throw meshgen_limit(
                "vcp::bfem::meshgen: max_refine exceeded (layers)");
        dz = dz / T(2);
        ++zsteps;
    }
    // z levels by repeated midpoint bisection (2^zsteps + 1 values)
    std::vector<T> zs;
    zs.push_back(z0);
    zs.push_back(z1);
    for (int s = 0; s < zsteps; ++s) {
        std::vector<T> nz;
        for (std::size_t i = 0; i + 1 < zs.size(); ++i) {
            nz.push_back(zs[i]);
            nz.push_back((zs[i] + zs[i + 1]) / T(2));
        }
        nz.push_back(zs[zs.size() - 1]);
        zs.swap(nz);
    }
    const int nlay = static_cast<int>(zs.size()) - 1;

    // 3) vertices: layer L holds ids L*nv2 .. L*nv2 + nv2-1
    vertices.clear();
    for (int L = 0; L <= nlay; ++L)
        for (int i = 0; i < nv2; ++i) {
            std::array<T, 3> p =
                {{v2[static_cast<std::size_t>(i)][0],
                  v2[static_cast<std::size_t>(i)][1],
                  zs[static_cast<std::size_t>(L)]}};
            vertices.push_back(p);
        }

    // 4) prisms -> 3 tetrahedra (min-index rule)
    elements.clear();
    for (int L = 0; L < nlay; ++L) {
        const int off = L * nv2;
        for (std::size_t t = 0; t < e2.size(); ++t) {
            const std::array<int, 3>& tri = e2[t];
            // rotate so the smallest 2D id leads (CCW order preserved)
            int r = 0;
            if (tri[1] < tri[static_cast<std::size_t>(r)]) r = 1;
            if (tri[2] < tri[static_cast<std::size_t>(r)]) r = 2;
            const int b0 = tri[static_cast<std::size_t>(r)] + off;
            const int b1 = tri[static_cast<std::size_t>((r + 1) % 3)] + off;
            const int b2 = tri[static_cast<std::size_t>((r + 2) % 3)] + off;
            const int t0 = b0 + nv2;
            const int t1 = b1 + nv2;
            const int t2 = b2 + nv2;
            if (b1 < b2) {
                std::array<int, 4> e0 = {{b0, b1, b2, t2}};
                std::array<int, 4> e1 = {{b0, b1, t2, t1}};
                std::array<int, 4> e2t = {{b0, t1, t2, t0}};
                elements.push_back(e0);
                elements.push_back(e1);
                elements.push_back(e2t);
            } else {
                std::array<int, 4> e0 = {{b0, b1, b2, t1}};
                std::array<int, 4> e1 = {{b0, b2, t2, t1}};
                std::array<int, 4> e2t = {{b0, t0, t1, t2}};
                elements.push_back(e0);
                elements.push_back(e1);
                elements.push_back(e2t);
            }
        }
    }

    // 5) boundary face origins (lineage only, no geometric decision)
    std::map<std::array<int, 3>, meshgen_facet_source> fsrc;
    for (std::size_t i = 0; i < st2.boundary_facets.size(); ++i) {
        const int u = st2.boundary_facets[i][0];
        const int v = st2.boundary_facets[i][1];
        const int m = (u < v) ? u : v;
        const int M = (u < v) ? v : u;
        for (int L = 0; L < nlay; ++L) {
            const int off = L * nv2;
            // the quad splits along the diagonal from the smaller 2D id,
            // matching the prism rule above
            fsrc[face_key(m + off, M + off, M + off + nv2)] =
                st2.facet_source[i];
            fsrc[face_key(m + off, M + off + nv2, m + off + nv2)] =
                st2.facet_source[i];
        }
    }
    for (std::size_t t = 0; t < e2.size(); ++t) {
        const std::array<int, 3>& tri = e2[t];
        fsrc[face_key(tri[0], tri[1], tri[2])] =
            meshgen_facet_source(-1, static_cast<int>(t));       // bottom
        const int off = nlay * nv2;
        fsrc[face_key(tri[0] + off, tri[1] + off, tri[2] + off)] =
            meshgen_facet_source(-2, static_cast<int>(t));       // top
    }

    // 6) boundary facets: faces with adjacency count == 1 (F2 counting)
    std::map<std::array<int, 3>, int> cnt;
    for (std::size_t t = 0; t < elements.size(); ++t) {
        const std::array<int, 4>& e = elements[t];
        cnt[face_key(e[0], e[1], e[2])] += 1;
        cnt[face_key(e[0], e[1], e[3])] += 1;
        cnt[face_key(e[0], e[2], e[3])] += 1;
        cnt[face_key(e[1], e[2], e[3])] += 1;
    }
    status.boundary_facets.clear();
    status.facet_source.clear();
    status.boundary_vertices.clear();
    std::map<int, bool> bverts;
    for (std::map<std::array<int, 3>, int>::const_iterator it = cnt.begin();
         it != cnt.end(); ++it) {
        if (it->second != 1) continue;
        std::map<std::array<int, 3>, meshgen_facet_source>::const_iterator
            src = fsrc.find(it->first);
        if (src == fsrc.end())
            throw meshgen_degenerate(
                "vcp::bfem::meshgen: internal: boundary face without origin");
        status.boundary_facets.push_back(it->first);
        status.facet_source.push_back(src->second);
        bverts[it->first[0]] = true;
        bverts[it->first[1]] = true;
        bverts[it->first[2]] = true;
    }
    for (std::map<int, bool>::const_iterator it = bverts.begin();
         it != bverts.end(); ++it)
        status.boundary_vertices.push_back(it->first);
    status.refine_steps = st2.refine_steps;   // cross-section sweeps
}

} // namespace meshgen_detail

// ---------------------------------------------------------------------------
// public entry points, Phase C (design section 6; same 3 forms as 2D plus
// the Delaunay-cross-section form)
// ---------------------------------------------------------------------------
template <typename T>
void generate_mesh_lists(const extruded_domain<T>& dom, const T& h,
                         std::vector<std::array<T, 3> >& vertices,
                         std::vector<std::array<int, 4> >& elements,
                         meshgen_status<3>& status,
                         const meshgen_options& opt = meshgen_options()) {
    if (!meshgen_detail::certainly_pos(h))
        throw meshgen_error("vcp::bfem::meshgen: h must be certainly positive");
    status = meshgen_status<3>();
    meshgen_detail::coarse_mesh<T> cm =
        meshgen_detail::build_coarse_mesh(dom.base, opt);
    meshgen_detail::extrude_and_finalize(cm, dom.z0, dom.z1, h, opt, status,
                                         vertices, elements);
}

template <typename T>
mesh<3, T> generate_mesh(const extruded_domain<T>& dom, const T& h,
                         meshgen_status<3>& status,
                         const meshgen_options& opt = meshgen_options()) {
    std::vector<std::array<T, 3> > vertices;
    std::vector<std::array<int, 4> > elements;
    generate_mesh_lists(dom, h, vertices, elements, status, opt);
    return mesh<3, T>::from_lists(vertices, elements);
}

template <typename T>
mesh<3, T> generate_mesh(const extruded_domain<T>& dom, const T& h,
                         const meshgen_options& opt = meshgen_options()) {
    meshgen_status<3> status;
    return generate_mesh(dom, h, status, opt);
}

template <typename T>
mesh<3, T> generate_mesh_delaunay(const extruded_domain<T>& dom, const T& h,
                                  meshgen_status<3>& status,
                                  const meshgen_options& opt = meshgen_options()) {
    if (!meshgen_detail::certainly_pos(h))
        throw meshgen_error("vcp::bfem::meshgen: h must be certainly positive");
    status = meshgen_status<3>();
    meshgen_detail::delaunay_coarse<T> dc =
        meshgen_detail::delaunay_flip_coarse(dom.base, opt);
    std::vector<std::array<T, 3> > vertices;
    std::vector<std::array<int, 4> > elements;
    meshgen_detail::extrude_and_finalize(dc.cm, dom.z0, dom.z1, h, opt,
                                         status, vertices, elements);
    status.delaunay_complete = dc.complete;
    return mesh<3, T>::from_lists(vertices, elements);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN3_HPP
