// vcp/bfem/constants/detail/poisson_dict_impl.hpp
//
// CONST-C1: mesh-level consumers of the dictionary resolve layer -- the
// k dependent squared L^2 projection error constant C_{0,h}(d)^2 of a
// whole mesh, its per-element bulk form, and the squared Ritz projection
// error constant C_h(k)^2 composed from the UNCHANGED hypercircle kappa^2
// of poisson_constants.hpp.
//
// This is a DETAIL header (owner ruling on CONST-C1 STOP-1, point two):
// the one and only user entry point is
// vcp/bfem/constants/poisson_constants.hpp, whose last line includes this
// file.  Including this file directly is refused below, so the public
// include surface of the constants layer stays the single header it always
// was, and every translation unit that includes poisson_constants.hpp sees
// the whole API regardless of what else it includes and in which order.
//
// DEGREE ARGUMENT CONVENTION (CONST-C1 design section one, restated on
// every function below): d is the degree of the PROJECTION TARGET space
// P^d, i.e. || v - pi_h v || <= C_{0,h}(d) |v|_1 projects onto P^d.  In
// the conforming P^k hypercircle context the degree served is d = k - 1,
// and that k -> d conversion happens INSIDE
// ritz_projection_error_constant_h01_sq only; no k-taking alias of the l2
// functions exists (two int overloads of the same signature would be a
// misuse source).
//
// WHY THE CONDITIONAL RESOLVE INCLUDE, THE FORWARD DECLARATIONS AND THE
// Src TEMPLATE PARAMETER EXIST (include-cycle robustness).
// dict/resolve.hpp reaches back to poisson_constants.hpp through
// dict/dict_entry.hpp -> element_projection.hpp, so when a translation
// unit enters the constants cluster through one of THOSE headers instead
// of poisson_constants.hpp (the constb1 / constb2a / constb2b suites do),
// this file is parsed at a moment where the body of the entry header is
// still incomplete.  Three measures make this file well formed in EVERY
// include order:
//   (a) the resolve include below is guarded on element_projection.hpp's
//       header guard.  That guard being set here means this file is being
//       parsed from INSIDE the cluster (every cluster path into
//       poisson_constants.hpp passes through element_projection.hpp, and
//       no order completes element_projection.hpp before the tail of
//       poisson_constants.hpp runs); starting the resolve chain at such a
//       moment would re-enter the dict headers on top of an incomplete
//       dependency.  When the guard is NOT set -- the ordinary case: the
//       translation unit entered through poisson_constants.hpp itself --
//       resolve.hpp is included and the whole cluster completes here;
//   (b) the resolve layer is named only through the forward declarations
//       below, which the real definitions in resolve.hpp legally
//       complete, in whichever order the cluster was entered;
//   (c) enumerator accesses are template-dependent expressions (the Src
//       parameter defaulting to l2_source), so their lookup is deferred
//       to the instantiation point.
// The API is thus parsed and visible in EVERY translation unit that
// includes poisson_constants.hpp; no entity of this file appears or
// disappears with the include order.  One residual constraint, stated
// honestly: a translation unit that (unusually) enters the cluster
// through element_projection.hpp or dict/dict_entry.hpp AND then calls
// the mesh-level API must also include a header that defines the resolve
// layer (poisson first, or dict/resolve.hpp) -- otherwise the
// instantiation fails on the still-incomplete forward declarations.  A
// unit following the documented single entry point never sees this.
//
// Authority: sandbox/docs/design/CONST-C1_design_v1.0.md (with the STOP-1
// owner ruling: detail-header placement, this file), and
// sandbox/docs/plans/CONST-C1_implementation_directive_v1.0.md.
//
// Lexical regime: traditional (no decimal literal, bare or in string; the
// self-check lives in sandbox/tests/constc1_tests.cpp, consta-style).

#ifndef VCP_BFEM_CONSTANTS_DETAIL_POISSON_DICT_IMPL_HPP
#define VCP_BFEM_CONSTANTS_DETAIL_POISSON_DICT_IMPL_HPP

#ifndef VCP_BFEM_CONSTANTS_POISSON_CONSTANTS_HPP
#error "poisson_dict_impl.hpp is a detail header: include <vcp/bfem/constants/poisson_constants.hpp> instead"
#endif

#include <vector>
#include <array>
#include <cstddef>

// guarded on the element_projection.hpp header guard: see measure (a) in
// the header comment (set here == parsed from inside the cluster, where
// starting the resolve chain would re-enter dict headers on top of an
// incomplete dependency; unset == poisson-first, the ordinary case)
#ifndef VCP_BFEM_CONSTANTS_ELEMENT_PROJECTION_HPP
#include <vcp/bfem/constants/dict/resolve.hpp>
#endif

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// forward declarations of the resolve layer (see the header comment).  The
// real definitions live in dict/resolve.hpp; these declarations either
// precede them (cluster entered through a dict/ or engine header) or
// restate them (cluster entered here), and both orders are legal.
// ---------------------------------------------------------------------------
enum class l2_source;
template <int D, typename T> struct l2_projection_resolution;
template <int D, typename T>
l2_projection_resolution<D, T>
resolve_l2_projection_constant_sq(const T vertices[D + 1][D], int d);

// ---------------------------------------------------------------------------
// mesh_resolution_report (design section one): where the served constants
// came from, reported honestly.  Filled by the _sq entry points below.
//
//   n_elements       elements visited (== the mesh's element count)
//   count_*          elements per resolution source; their sum equals
//                    n_elements on every successful return (count_none is
//                    incremented just before the none-throw, so a caller
//                    only ever reads it nonzero from a report object it
//                    passed to a call that then threw -- rep is written on
//                    success only)
//   min_served_d     the smallest served_d over ALL elements, i.e. the
//                    coarsest degree actually served anywhere (== d when
//                    every element is served at the requested degree; a
//                    p0_closed_form element pins it to zero)
//   worst_element    the element attaining the mesh maximum
//   worst_source     that element's resolution source
//
// The default constructor VALUE-initializes worst_source (a forward
// declared scoped enum is a complete type, but its enumerators cannot be
// named here -- header comment); worst_element carrying the negative
// sentinel is the "not filled yet" marker.  Every successful _sq call
// overwrites the whole struct and resets worst_source to the none
// enumerator before its element scan.
// ---------------------------------------------------------------------------
struct mesh_resolution_report {
    long long n_elements;
    long long count_registry;
    long long count_coverage;
    long long count_envelope;
    long long count_p0;
    long long count_none;
    int min_served_d;
    long long worst_element;
    l2_source worst_source;

    mesh_resolution_report()
        : n_elements(0), count_registry(0), count_coverage(0),
          count_envelope(0), count_p0(0), count_none(0),
          min_served_d(-1), worst_element(-1), worst_source() {}
};

// ---------------------------------------------------------------------------
// l2_projection_element_constants_sq (design section one): the per-element
// vector of squared constants, entry e an upper bound of C_d(K_e)^2
// obtained as resolve(K_e, d).cd_sq_over_h_sq * h^2(K_e), in element
// order (the element-wise residual indicator form).  Entry e equals that
// direct per-element resolve composition bit for bit (gate G-E8).
//
// d is the PROJECTION TARGET degree (P^d; in the conforming P^k
// hypercircle context d = k - 1 -- that conversion lives inside
// ritz_projection_error_constant_h01_sq below and nowhere else).
//
// An element the dictionary cannot serve at all (resolution source none --
// impossible for D = 2, where the P^0 closed form always answers) throws
// vcp::verification_error naming the element.  Exceptions of resolve
// itself (d negative, genuine-interval or non-dyadic vertex, degenerate
// element on the closed-form path) propagate unchanged.
//
// The report, when asked for, carries the full census; the maximum rule
// is the file-wide one of poisson_constants.hpp (design 5.2 there): the
// LARGEST upper end wins, the comparison is made on the upper ends, ties
// keep the earlier element.
//
// Src is an implementation detail (deferred lookup, header comment); leave
// it defaulted.
// ---------------------------------------------------------------------------
template <typename T, typename Src = l2_source>
std::vector<T> l2_projection_element_constants_sq(
        const vcp::bfem::mesh<2, T>& Th, int d,
        mesh_resolution_report* rep = nullptr) {
    if (Th.num_elements() <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_constants_sq: "
            "mesh has no element");

    mesh_resolution_report r;
    r.n_elements = Th.num_elements();
    r.worst_source = Src::none;

    std::vector<T> out;
    out.reserve(static_cast<std::size_t>(Th.num_elements()));
    for (int e = 0; e < Th.num_elements(); ++e) {
        const std::array<int, 3>& el = Th.element(e);
        T verts[3][2];
        for (int i = 0; i < 3; ++i) {
            const std::array<T, 2>& p = Th.vertex(el[i]);
            verts[i][0] = p[0];
            verts[i][1] = p[1];
        }
        const l2_projection_resolution<2, T> res =
            resolve_l2_projection_constant_sq<2, T>(verts, d);
        if (!res.ok) {
            ++r.count_none;
            vcp::throw_error<vcp::verification_error>(
                "vcp::bfem::constants::l2_projection_element_constants_sq: "
                "the dictionary cannot serve element e = ", e,
                " at degree d = ", d, " (resolution source none)");
        }
        switch (res.source) {
        case Src::registry_exact:  ++r.count_registry; break;
        case Src::coverage_cell:   ++r.count_coverage; break;
        case Src::degree_envelope: ++r.count_envelope; break;
        case Src::p0_closed_form:  ++r.count_p0;       break;
        default:                   break;   // none is unreachable: thrown above
        }
        if (r.min_served_d < 0 || res.served_d < r.min_served_d)
            r.min_served_d = res.served_d;

        const T val = res.cd_sq_over_h_sq * res.h_sq;
        if (out.empty() ||
            out[static_cast<std::size_t>(r.worst_element)].upper() <
                val.upper()) {
            r.worst_element = e;
            r.worst_source = res.source;
        }
        out.push_back(val);
    }
    if (rep != nullptr) *rep = r;
    return out;
}

// ---------------------------------------------------------------------------
// l2_projection_error_constant_sq (design section one):
//
//     C_{0,h}(d)^2 = max_K [ resolve(K, d).cd_sq_over_h_sq * h^2(K) ]
//
// the mesh-level squared L^2 projection error constant onto P^d, served by
// the dictionary (registry -> coverage cell -> degree envelope -> P^0
// closed form).  SQUARED deliberately (ruling R17): no square root on this
// path, so an exact fraction scalar T closes end to end (the resolve
// strings parse exactly, h^2 is exact, the maximum is a total-order
// comparison on such a T).
//
// d is the PROJECTION TARGET degree (P^d; conforming P^k hypercircle
// context: d = k - 1, converted only inside
// ritz_projection_error_constant_h01_sq below).
//
// The k independent P^0 predecessor l2_projection_error_constant(Th) of
// poisson_constants.hpp is untouched (chapter eight reproduction path);
// this is its dictionary served, degree dependent, squared-domain
// successor.  The maximum keeps the file-wide rule: the largest upper end
// wins and the selected interval is returned unchanged.
// ---------------------------------------------------------------------------
template <typename T>
T l2_projection_error_constant_sq(const vcp::bfem::mesh<2, T>& Th, int d,
                                  mesh_resolution_report* rep = nullptr) {
    mesh_resolution_report r;
    const std::vector<T> v = l2_projection_element_constants_sq(Th, d, &r);
    if (rep != nullptr) *rep = r;
    return v[static_cast<std::size_t>(r.worst_element)];
}

// ---------------------------------------------------------------------------
// l2_projection_element_constants_sq, mesh<3, T> (CONST-C2 design section
// two): the three dimensional twin of the pair above.  Same semantics, same
// report, same maximum rule; the ONLY structural difference is the
// resolution order the dictionary can offer in three dimensions.
//
// RESOLUTION ORDER IN 3D.  resolve_l2_projection_constant_sq's layers (2)
// coverage_cell and (4) p0_closed_form are BOTH guarded on D == 2 -- the
// coverage grid is a two dimensional (a, b^2) normal form, and the P^0
// closed form l2_projection_element_bound of poisson_constants.hpp takes
// three PLANAR vertices and has no three dimensional counterpart
// (CONST-C2 P0 reconnaissance; design section two's read-through, owner
// approved 2026-08-20).  So the 3D order is
//
//     registry_exact -> degree_envelope -> none
//
// and count_coverage / count_p0 stay zero on every 3D call.  In 2D the P^0
// closed form makes "none" unreachable; in 3D it IS reachable, because the
// generated registry holds a FINITE list of similarity classes.  A mesh
// whose element is outside that list therefore throws, and the message
// names the explicit escape hatch: ondemand_l2_projection_constant_sq of
// dict/resolve.hpp computes the class on demand and returns a ledger entry
// for dict/ondemand_3d.hpp.
//
// GENERATION TIME ENVELOPE against RUNTIME ENVELOPE -- do not confuse the
// two (owner condition on the CONST-C2 census ruling, 2026-08-20).  The
// generated registry_3d.hpp ships every one of its 51 classes at
// d = 0 .. 8; the entries for d = 3 .. 8 carry a value that was INHERITED
// from d' = 2 when the table was generated (their provenance comment
// records it), but they are real, exactly keyed entries, so resolve fires
// them as registry_exact with served_d = d.  The runtime degree_envelope
// layer therefore only starts at the registry's permanently missing degree
// d = 9 (engine factorial cap, see registry_3d.hpp's missing list), where
// it inherits the d = 8 string and reports served_d = 8.  In particular a
// min_served_d of 2 does NOT show up in a 3D census, and min_served_d is
// not a measure of the information actually behind the value.
//
// Src is an implementation detail (deferred enumerator lookup, header
// comment); leave it defaulted.
// ---------------------------------------------------------------------------
template <typename T, typename Src = l2_source>
std::vector<T> l2_projection_element_constants_sq(
        const vcp::bfem::mesh<3, T>& Th, int d,
        mesh_resolution_report* rep = nullptr) {
    if (Th.num_elements() <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_constants_sq: "
            "mesh has no element");

    mesh_resolution_report r;
    r.n_elements = Th.num_elements();
    r.worst_source = Src::none;

    std::vector<T> out;
    out.reserve(static_cast<std::size_t>(Th.num_elements()));
    for (int e = 0; e < Th.num_elements(); ++e) {
        const std::array<int, 4>& el = Th.element(e);
        T verts[4][3];
        for (int i = 0; i < 4; ++i) {
            const std::array<T, 3>& p = Th.vertex(el[i]);
            verts[i][0] = p[0];
            verts[i][1] = p[1];
            verts[i][2] = p[2];
        }
        const l2_projection_resolution<3, T> res =
            resolve_l2_projection_constant_sq<3, T>(verts, d);
        if (!res.ok) {
            ++r.count_none;
            vcp::throw_error<vcp::verification_error>(
                "vcp::bfem::constants::l2_projection_element_constants_sq: "
                "the dictionary cannot serve element e = ", e,
                " at degree d = ", d, " (resolution source none; in three "
                "dimensions there is no P^0 closed-form fallback, so the "
                "element's similarity class is simply not in the registry "
                "-- compute it explicitly with "
                "ondemand_l2_projection_constant_sq<3, T> and append the "
                "returned entry to dict/ondemand_3d.hpp)");
        }
        switch (res.source) {
        case Src::registry_exact:  ++r.count_registry; break;
        case Src::coverage_cell:   ++r.count_coverage; break;
        case Src::degree_envelope: ++r.count_envelope; break;
        case Src::p0_closed_form:  ++r.count_p0;       break;
        default:                   break;   // none is unreachable: thrown above
        }
        if (r.min_served_d < 0 || res.served_d < r.min_served_d)
            r.min_served_d = res.served_d;

        const T val = res.cd_sq_over_h_sq * res.h_sq;
        if (out.empty() ||
            out[static_cast<std::size_t>(r.worst_element)].upper() <
                val.upper()) {
            r.worst_element = e;
            r.worst_source = res.source;
        }
        out.push_back(val);
    }
    if (rep != nullptr) *rep = r;
    return out;
}

// ---------------------------------------------------------------------------
// l2_projection_error_constant_sq, mesh<3, T> (CONST-C2 design section two):
//
//     C_{0,h}(d)^2 = max_K [ resolve(K, d).cd_sq_over_h_sq * h^2(K) ]
//
// over a tetrahedral mesh, served by the dictionary's 3D order
// (registry -> degree envelope; see the comment above for why the coverage
// and P^0 layers do not participate, and for the generation time against
// runtime envelope distinction).  Squared for the same reason as the 2D
// twin (ruling R17): no square root on this path, so an exact fraction
// scalar T closes end to end.  Maximum rule unchanged: the largest upper
// end wins, the comparison is on the upper ends, ties keep the earlier
// element, and the selected interval is returned unchanged.
//
// There is deliberately NO 3D ritz_projection_error_constant_h01_sq: the
// hypercircle kappa_h it composes with is mesh<2, T> only
// (hypercircle_kappa_h01_squared of poisson_constants.hpp), so the
// composition has no three dimensional meaning yet.  CONST-C2 design
// section seven keeps that as an independent track.
// ---------------------------------------------------------------------------
template <typename T>
T l2_projection_error_constant_sq(const vcp::bfem::mesh<3, T>& Th, int d,
                                  mesh_resolution_report* rep = nullptr) {
    mesh_resolution_report r;
    const std::vector<T> v = l2_projection_element_constants_sq(Th, d, &r);
    if (rep != nullptr) *rep = r;
    return v[static_cast<std::size_t>(r.worst_element)];
}

// ---------------------------------------------------------------------------
// ritz_projection_error_constant_h01_sq (design section one):
//
//     C_h(k)^2 = hypercircle_kappa_h01_squared(Th, k) + C_{0,h}(k - 1)^2
//
// the squared composition of Theorem 8.3 / (8.25) with the dictionary
// served C_{0,h} in place of the fixed P^0 constant: X_h = P^{k-1}, so the
// projection the hypercircle argument applies is served at degree
// d = k - 1 instead of degree zero, closing the slack (c) of the k >= 2
// discussion in poisson_constants.hpp.  THIS ADDITION IS THE ONLY k -> d
// CONVERSION POINT of the layer (gate G-E4 pins it).
//
// Same argument list as ritz_projection_error_constant_h01 of
// poisson_constants.hpp plus the optional report (filled from the C_{0,h}
// pass); kappa^2 is consumed UNCHANGED (its own dense pass, certification
// and exceptions as documented there), and the returned value is the
// SQUARE of C_h -- no square root is taken, matching the squared-domain
// dictionary exit.
// ---------------------------------------------------------------------------
template <typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T ritz_projection_error_constant_h01_sq(const vcp::bfem::mesh<2, T>& Th, int k,
                                    mesh_resolution_report* rep = nullptr) {
    // kappa first: it owns the k >= 1 validation of the classic dense path
    const T kappa2 = hypercircle_kappa_h01_squared<T, DP, SP>(Th, k);
    const T c0h_sq = l2_projection_error_constant_sq(Th, k - 1, rep);
    return kappa2 + c0h_sq;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DETAIL_POISSON_DICT_IMPL_HPP
