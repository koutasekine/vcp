// vcp/bfem/meshgen_collapse.hpp
// MG-2 (B): collapse the INTERIOR vertices of an interval-typed mesh to
// representative points (width 0), design MG-2 v1.1 section 6.
//
// Foundation (design 1.3, owner's principle): with the combinatorial
// structure fixed, every triangle certainly CCW, and the BOUNDARY vertices
// untouched, the union of the elements is the boundary polygon itself, so
// the interior vertices are free parameters and may move to any
// representative point inside their coordinate boxes.
//
// Why there is NO "restore if broken" branch (design 6.2, ruling Q3/C3):
// by inclusion monotonicity of interval arithmetic, re-evaluating orient2d
// over subsets of the input boxes yields a subset of the original interval
// determinant, so a certainly positive orientation cannot become
// indefinite or reversed. The claim is machine-verified by T-C1 instead of
// being trusted (design section 10).
//
// Representative point (design 6.3): the midpoint candidate
// (lower + upper) / 2 is checked BY COMPARISON to lie inside the box; if
// the check fails (rounding-mode effects), the lower endpoint -- trivially
// inside -- is used and the fallback is counted in the report.
//
// Point types are a no-op (the report stays zeroed except the interior
// count contract below). Boundary vertices are never written to at all,
// so they stay bit-identical.
//
// Lexical policy (MG-1 section 4.1): no decimal literals, no fp type
// tokens, no sqrt/abs/min/max, integer literals only; the interval
// endpoints are handled in the base scalar type of the interval itself.

#ifndef VCP_BFEM_MESHGEN_COLLAPSE_HPP
#define VCP_BFEM_MESHGEN_COLLAPSE_HPP

#include <array>
#include <vector>

#include <kv/interval.hpp>

#include <vcp/bfem/meshgen.hpp>

namespace vcp {
namespace bfem {

// widths are only exactly representable in the scalar family of T itself,
// and the lexical policy bans fixed fp types here, so the report carries
// the maximum boundary width as a T (a width-0 value of T for point types)
template <typename T>
struct meshgen_collapse_report {
    long long collapsed_interior;      // interior vertices given a rep point
    long long fallback_count;          // per-coordinate lower() fallbacks
    long long wide_boundary_vertices;  // boundary vertices with width > 0
    T max_boundary_width;              // largest boundary coordinate width
    meshgen_collapse_report()
        : collapsed_interior(0), fallback_count(0),
          wide_boundary_vertices(0), max_boundary_width(T(0)) {}
};

namespace meshgen_collapse_detail {

// primary template: point scalars -- nothing to collapse
template <typename T>
struct collapse_traits {
    static bool is_interval() { return false; }
    static void collapse(T&, long long&) {}
    static bool positive_width(const T&) { return false; }
    // width as a T; never called with a positive result for point types
    static T width(const T&) { return T(0); }
    static bool width_less(const T&, const T&) { return false; }
};

// kv::interval<X>: the only interval family of the workspace (the same
// scope as the meshgen_io_tags / meshgen_convert specializations)
template <typename X>
struct collapse_traits<kv::interval<X> > {
    typedef kv::interval<X> T;
    static bool is_interval() { return true; }

    // design 6.3: midpoint candidate, containment BY COMPARISON, lower()
    // fallback (counted); the assignment builds a width-0 point interval
    static void collapse(T& x, long long& fallback_count) {
        X m = (x.lower() + x.upper()) / X(2);
        if (!(x.lower() <= m && m <= x.upper())) {
            m = x.lower();
            fallback_count = fallback_count + 1;
        }
        x = T(m);
    }

    static bool positive_width(const T& x) { return x.lower() < x.upper(); }

    static T width(const T& x) { return T(x.upper() - x.lower()); }

    static bool width_less(const T& a, const T& b) {
        return a.upper() < b.upper();
    }
};

} // namespace meshgen_collapse_detail

// ---------------------------------------------------------------------------
// public API (design 6.1). Interior = vertices NOT listed in
// status.boundary_vertices. Boundary vertices are never written; interval
// interiors become width-0 representative points; point types no-op.
// ---------------------------------------------------------------------------
template <int D, typename T>
void collapse_interior_to_points(std::vector<std::array<T, D> >& vertices,
                                 const meshgen_status<D>& status,
                                 meshgen_collapse_report<T>& rep) {
    typedef meshgen_collapse_detail::collapse_traits<T> tr;
    rep = meshgen_collapse_report<T>();

    const int nv = static_cast<int>(vertices.size());
    std::vector<char> on_boundary(static_cast<std::size_t>(nv), 0);
    for (std::size_t i = 0; i < status.boundary_vertices.size(); ++i) {
        const int v = status.boundary_vertices[i];
        if (v < 0 || v >= nv)
            throw meshgen_error(
                "vcp::bfem::meshgen_collapse: boundary vertex id out of "
                "range for the vertex list");
        on_boundary[static_cast<std::size_t>(v)] = 1;
    }

    if (!tr::is_interval()) return;   // point types: no-op, zeroed report

    for (int v = 0; v < nv; ++v) {
        if (on_boundary[static_cast<std::size_t>(v)] != 0) {
            // never written; only measured for the report
            bool any_wide = false;
            for (int d = 0; d < D; ++d) {
                const T& x = vertices[static_cast<std::size_t>(v)]
                                     [static_cast<std::size_t>(d)];
                if (tr::positive_width(x)) {
                    any_wide = true;
                    const T w = tr::width(x);
                    if (tr::width_less(rep.max_boundary_width, w))
                        rep.max_boundary_width = w;
                }
            }
            if (any_wide)
                rep.wide_boundary_vertices = rep.wide_boundary_vertices + 1;
            continue;
        }
        for (int d = 0; d < D; ++d)
            tr::collapse(vertices[static_cast<std::size_t>(v)]
                                 [static_cast<std::size_t>(d)],
                         rep.fallback_count);
        rep.collapsed_interior = rep.collapsed_interior + 1;
    }
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_COLLAPSE_HPP
