// vcp/bfem/constants/dict/dict_entry.hpp
//
// CONST-B2a: dictionary first layer -- the exact similarity-class key and the
// registry lookup for the squared polynomial L^2 projection error constant
// of one simplex element (D = 2 triangle, D = 3 tetrahedron).
//
// Key (design 1.1).  The similarity class of a rational-vertex simplex is
// decided EXACTLY through the squared-edge-length matrix:
//
//   - the squared edge lengths l_ij^2 = |v_i - v_j|^2 are exact rationals;
//   - the list is normalized by its MINIMUM element (scale elimination: the
//     similarity RATIO may be irrational, the ratio of squared lengths is
//     rational whenever both simplices have rational vertices);
//   - over all (D+1)! vertex relabelings the upper-triangular row-major list
//     (ordered pairs (i,j), i < j) is taken lexicographically MINIMAL.
//
//   The squared-edge-length matrix determines the simplex up to isometry and
//   reflection (Cayley-Menger; machine-checked by gate G-C1), and C_d is
//   invariant under isometry, reflection and scaling (sigma note corollary,
//   CONST-B1.2 G-B11), so key equality implies that the registry value is
//   shared soundly.  Fractions are compared AS FRACTIONS -- nothing is
//   decimalized on this path.
//
// Stored value (design 1.2).  A registry entry carries the dimensionless
//   Chat^2 := C_d(K_rep)^2 / h^2(K_rep)
// as a 12-significant-digit UPWARD decimal string; the consumption inequality
// is C_d(K)^2 <= Chat^2 * h^2(K), square-root free end to end (R17).  The
// string enters the requested scalar type exclusively through the CONST-A
// conversion layer (constant_from_string, literal.hpp); this header carries
// no numeric constant of its own.
//
// Registry access.  The generated tables live in dict/registry_2d.hpp and
// dict/registry_3d.hpp (regime B: decimal strings inside their own
// VCP_CONSTANTS_TABLE marker block).  The two headers are included from here
// AFTER the entry structs and the access-function declarations, so both
// inclusion orders (dict_entry first or a registry first) compile.
//
// Lookup contract.  lookup_l2_projection_class<D, T> mirrors the engine
// entry point of element_projection.hpp: vertices are POINT intervals with
// finite binary-fraction (dyadic) coordinates -- a genuine interval or a
// non-dyadic coordinate throws vcp::invalid_argument through the shared
// rationalization helper.  A not-found key returns found == false and
// nothing else; no fallback of any kind lives in this layer (B-2b).
//
// Lexical regime: traditional (no decimal literal, bare or in string).
//
// Authority: sandbox/docs/design/CONST-B2a_design_v1.0.md and
// sandbox/docs/plans/CONST-B2a_implementation_directive_v1.0.md.

#ifndef VCP_BFEM_CONSTANTS_DICT_DICT_ENTRY_HPP
#define VCP_BFEM_CONSTANTS_DICT_DICT_ENTRY_HPP

#include <array>
#include <algorithm>

#include <vcp/bfem/constants/element_projection.hpp>
#include <vcp/bfem/constants/literal.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// registry entry layout (design 1.3, R21).  The canonical key is stored as
// integer fractions (numerator / denominator per component, in the canonical
// component order); the value is the single decimal string.  Everything else
// about an entry (representative vertices, source, level, cap or envelope
// notes) lives in the provenance comment above it inside the registry file.
// ---------------------------------------------------------------------------
struct l2_projection_registry_entry_2d {
    long long key_num[3];
    long long key_den[3];
    int d;
    const char* cd_sq_over_h_sq;
};

struct l2_projection_registry_entry_3d {
    long long key_num[6];
    long long key_den[6];
    int d;
    const char* cd_sq_over_h_sq;
};

// generated tables (defined inline in the registry headers included below)
const l2_projection_registry_entry_2d*
l2_projection_registry_2d_table(int& count);
const l2_projection_registry_entry_3d*
l2_projection_registry_3d_table(int& count);

} // namespace constants
} // namespace bfem
} // namespace vcp

#include <vcp/bfem/constants/dict/registry_2d.hpp>
#include <vcp/bfem/constants/dict/registry_3d.hpp>

namespace vcp {
namespace bfem {
namespace constants {
namespace detail {

// ---------------------------------------------------------------------------
// exact similarity-class key of a rational-vertex simplex.
// dict_edge_count<D>: 3 (D = 2) / 6 (D = 3).
// ---------------------------------------------------------------------------
template <int D>
struct dict_key {
    enum { edge_count = (D * (D + 1)) / 2 };
    std::array<ep_rational, static_cast<std::size_t>(edge_count)> e;
};

template <int D>
bool dict_edge_list_less(
        const std::array<ep_rational,
                         static_cast<std::size_t>(dict_key<D>::edge_count)>& a,
        const std::array<ep_rational,
                         static_cast<std::size_t>(dict_key<D>::edge_count)>& b) {
    for (std::size_t k = 0;
         k < static_cast<std::size_t>(dict_key<D>::edge_count); ++k) {
        if (a[k] < b[k]) return true;
        if (b[k] < a[k]) return false;
    }
    return false;
}

// squared distances under a relabeling perm, upper-triangular row-major
template <int D>
void dict_edge_list_of(
        const std::array<std::array<ep_rational, D>, D + 1>& v,
        const int* perm,
        std::array<ep_rational,
                   static_cast<std::size_t>(dict_key<D>::edge_count)>& out) {
    std::size_t k = 0;
    for (int i = 0; i <= D; ++i)
        for (int j = i + 1; j <= D; ++j) {
            ep_rational s(0);
            for (int c = 0; c < D; ++c) {
                const ep_rational t =
                    v[static_cast<std::size_t>(perm[i])]
                     [static_cast<std::size_t>(c)]
                  - v[static_cast<std::size_t>(perm[j])]
                     [static_cast<std::size_t>(c)];
                s += t * t;
            }
            out[k++] = s;
        }
}

// exact h^2 = max squared edge length (identity permutation suffices)
template <int D>
ep_rational dict_max_edge_sq(
        const std::array<std::array<ep_rational, D>, D + 1>& v) {
    int perm[D + 1];
    for (int i = 0; i <= D; ++i) perm[i] = i;
    std::array<ep_rational,
               static_cast<std::size_t>(dict_key<D>::edge_count)> e;
    dict_edge_list_of<D>(v, perm, e);
    ep_rational mx = e[0];
    for (std::size_t k = 1;
         k < static_cast<std::size_t>(dict_key<D>::edge_count); ++k)
        if (mx < e[k]) mx = e[k];
    return mx;
}

// the canonical key: normalize by the minimum, then lexicographic minimum
// over all (D+1)! relabelings.  A degenerate simplex (some edge length zero,
// hence a zero minimum) is rejected before the division.
template <int D>
dict_key<D> dict_key_of(
        const std::array<std::array<ep_rational, D>, D + 1>& v) {
    int perm[D + 1];
    for (int i = 0; i <= D; ++i) perm[i] = i;
    std::array<ep_rational,
               static_cast<std::size_t>(dict_key<D>::edge_count)> raw;
    dict_edge_list_of<D>(v, perm, raw);
    ep_rational mn = raw[0];
    for (std::size_t k = 1;
         k < static_cast<std::size_t>(dict_key<D>::edge_count); ++k)
        if (raw[k] < mn) mn = raw[k];
    if (mn.is_zero())
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::lookup_l2_projection_class: two vertices "
            "coincide (a squared edge length is zero); the simplex is "
            "degenerate");

    dict_key<D> best;
    bool first = true;
    do {
        std::array<ep_rational,
                   static_cast<std::size_t>(dict_key<D>::edge_count)> cur;
        dict_edge_list_of<D>(v, perm, cur);
        for (std::size_t k = 0;
             k < static_cast<std::size_t>(dict_key<D>::edge_count); ++k)
            cur[k] = cur[k] / mn;
        if (first || dict_edge_list_less<D>(cur, best.e)) {
            best.e = cur;
            first = false;
        }
    } while (std::next_permutation(perm, perm + D + 1));
    return best;
}

// exact key == stored integer-fraction key?
template <int D, typename E>
bool dict_key_matches(const dict_key<D>& k, const E& entry) {
    for (std::size_t c = 0;
         c < static_cast<std::size_t>(dict_key<D>::edge_count); ++c) {
        const ep_rational stored(
            static_cast<long long>(entry.key_num[c]),
            static_cast<long long>(entry.key_den[c]));
        if (!(stored == k.e[c])) return false;
    }
    return true;
}

// registry selection by dimension
template <int D>
struct dict_registry_of;

template <>
struct dict_registry_of<2> {
    typedef l2_projection_registry_entry_2d entry_type;
    static const entry_type* table(int& count) {
        return l2_projection_registry_2d_table(count);
    }
};

template <>
struct dict_registry_of<3> {
    typedef l2_projection_registry_entry_3d entry_type;
    static const entry_type* table(int& count) {
        return l2_projection_registry_3d_table(count);
    }
};

} // namespace detail

// ---------------------------------------------------------------------------
// lookup result.  found == true: cd_sq_over_h_sq is the registry Chat^2 in
// the requested type (through constant_from_string -- outward enclosure for
// kv intervals, exact for rational, upward point value otherwise) and h_sq
// is an enclosure of the exact h^2(K) of the QUERIED simplex, so that
// C_d(K)^2 <= Chat^2 * h^2(K) is ready to consume.  found == false: only
// h_sq is meaningful; there is deliberately no fallback here (B-2b).
// ---------------------------------------------------------------------------
template <typename T>
struct l2_projection_class_lookup_result {
    bool found;
    T cd_sq_over_h_sq;
    T h_sq;

    l2_projection_class_lookup_result()
        : found(false), cd_sq_over_h_sq(), h_sq() {}
};

// ---------------------------------------------------------------------------
// lookup_l2_projection_class (design 4, directive P1): exact key of the
// vertices, then a linear scan of the generated registry for (key, d).
// ---------------------------------------------------------------------------
template <int D, typename T>
l2_projection_class_lookup_result<T>
lookup_l2_projection_class(const T vertices[D + 1][D], int d) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::lookup_l2_projection_class: "
                  "D must be 2 or 3");
    if (d < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::lookup_l2_projection_class: projection "
            "degree d must be >= 0 (got ", d, ")");

    std::array<std::array<detail::ep_rational, D>, D + 1> rv;
    for (int i = 0; i <= D; ++i)
        for (int c = 0; c < D; ++c)
            rv[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] =
                detail::rational_of_point_interval<T>(vertices[i][c]);

    const detail::dict_key<D> key = detail::dict_key_of<D>(rv);

    l2_projection_class_lookup_result<T> r;
    r.h_sq = detail::interval_of_exact_rational<T>(
        detail::dict_max_edge_sq<D>(rv));

    int count = 0;
    typedef typename detail::dict_registry_of<D>::entry_type entry_type;
    const entry_type* table = detail::dict_registry_of<D>::table(count);
    for (int i = 0; i < count; ++i) {
        if (table[i].d != d) continue;
        if (!detail::dict_key_matches<D>(key, table[i])) continue;
        r.found = true;
        r.cd_sq_over_h_sq =
            constant_from_string<T>::get(table[i].cd_sq_over_h_sq);
        return r;
    }
    return r;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DICT_DICT_ENTRY_HPP
