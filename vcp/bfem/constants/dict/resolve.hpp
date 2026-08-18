// vcp/bfem/constants/dict/resolve.hpp
//
// CONST-B2b: integrated lookup of the squared polynomial L^2 projection
// error constant over the whole dictionary (design section 2), plus the
// EXPLICIT ondemand computation API (design section 3).
//
// Resolution order (the fired source is always reported):
//   (1) registry_exact   exact canonical-key match at degree d, in the
//                        generated registry AND the ondemand ledgers
//                        (dict/ondemand_2d.hpp, ondemand_3d.hpp -- the
//                        ledgers are read on an equal footing, R21);
//   (2) coverage_cell    [D = 2 only] the (a, b^2) normal form lies in the
//                        covered domain (theta_min cut below) and the cell
//                        carries an entry for d (dict/coverage_2d.hpp);
//   (3) degree_envelope  the LARGEST d' < d served by (1) or (2)
//                        (monotonicity P^{d'} subset P^d => C_d <= C_{d'},
//                        the same lemma the registry inheritance uses);
//   (4) p0_closed_form   [D = 2 only] l2_projection_element_bound of
//                        poisson_constants.hpp, squared and divided by
//                        h^2(K) (valid for every d >= 0 since C_d <= C_0);
//   (5) none.
//
// Ondemand is NOT folded into resolve: a resolve call never starts an
// engine run (no implicit heavy computation).  The explicit API
// ondemand_l2_projection_constant_sq returns the value AND the ledger
// entry text; IT NEVER WRITES A FILE -- appending is the user's explicit
// operation (probe sandbox/probes/constb2b_ondemand_append.cpp).
//
// theta_min cut.  tau = 76542 / 10^7 is a rational LOWER bound of
// tan^2(5 degrees) (owner-ruled 2026-08-18; kv interval evaluation of
// tan(pi/36)^2 gives the certified lower end 0.0076542662455523309, and
// tau <= that end by exact fraction comparison -- machine-checked by gate
// G-D1t of constb2b_tests).  A class is INSIDE the covered angular domain
// iff b^2 >= tau * a^2 exactly; the smallest angle of the normal-form
// triangle (0,0), (1,0), (a, b) sits at the origin vertex with
// tan(theta) = b / a, so the cut serves a SUPERSET of { theta >= 5 deg }
// (tau is a lower bound), which only widens the coverage soundly.
//
// Normal form (design 1.1).  From the canonical similarity key of
// dict_entry.hpp (min-normalized squared edge list, permutation-minimal):
// let m be the largest component; remove ONE occurrence of m and let
// r_big >= r_small be the remaining two components divided by m.  Then
//   a = (1 + r_big - r_small) / 2,   b^2 = r_big - a^2
// (exact rationals; a in [1/2, 1], b^2 <= 1 - a^2 automatically).  Ties
// among longest edges collapse in the multiset, so the normal form is
// unique -- the same uniqueness the canonical key provides (B2a 1.1).
//
// Cell membership is decided by exact rational floor against the grid
// constants of coverage_2d.hpp; a query landing exactly on the top edge
// b^2 = 3/4 (the equilateral corner (1/2, 3/4)) or on an interior grid
// line belongs to the CLOSED cell chosen by the floor-plus-clamp rule
// below, and every such cell certifies its closure, so boundary points
// are always soundly served (owner gate addendum 2026-08-18; gates
// G-D1 / G-D3).
//
// Scalar contract: T is a kv::interval-like scalar (same contract as
// l2_projection_element_sq_bound; path (4) takes square roots and pi).
// Vertices are point intervals with dyadic coordinates -- the shared
// rationalization helper throws vcp::invalid_argument otherwise.
//
// Lexical regime: traditional (no decimal literal, bare or in string;
// every ledger string is BUILT AT RUN TIME from exact rationals).
//
// Authority: sandbox/docs/design/CONST-B2b_design_v1.0.md and
// sandbox/docs/plans/CONST-B2b_implementation_directive_v1.0.md.

#ifndef VCP_BFEM_CONSTANTS_DICT_RESOLVE_HPP
#define VCP_BFEM_CONSTANTS_DICT_RESOLVE_HPP

#include <string>
#include <array>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <unistd.h>

#include <vcp/bfem/constants/dict/dict_entry.hpp>
#include <vcp/bfem/constants/dict/coverage_2d.hpp>
#include <vcp/bfem/constants/dict/ondemand_2d.hpp>
#include <vcp/bfem/constants/dict/ondemand_3d.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// theta_min = 5 degrees: rational lower bound of tan^2(theta_min)
// (integer fraction; the ONLY new number this header introduces, ruled
// 2026-08-18 -- see the header comment for its derivation and gate)
// ---------------------------------------------------------------------------
constexpr long long l2_theta_min_tan_sq_num = 76542LL;
constexpr long long l2_theta_min_tan_sq_den = 10000000LL;

// ---------------------------------------------------------------------------
// resolution result (design 2)
// ---------------------------------------------------------------------------
enum class l2_source {
    registry_exact, coverage_cell, degree_envelope, p0_closed_form, none
};

template <int D, typename T>
struct l2_projection_resolution {
    bool ok;
    l2_source source;
    T cd_sq_over_h_sq;     // consumption: C_d(K)^2 <= cd_sq_over_h_sq * h^2(K)
    T h_sq;                // enclosure of the exact h^2(K) of the query
    int served_d;          // the d' actually served (== d except envelope;
                           // 0 on the p0_closed_form path)
    long long cell_i;      // coverage / envelope-via-coverage: cell indices
    long long cell_j;      // (-1 otherwise)

    l2_projection_resolution()
        : ok(false), source(l2_source::none), cd_sq_over_h_sq(), h_sq(),
          served_d(-1), cell_i(-1), cell_j(-1) {}
};

namespace detail {

// ---------------------------------------------------------------------------
// exact helpers
// ---------------------------------------------------------------------------

// floor of a nonnegative exact rational into long long (guarded)
inline long long b2b_floor_ll(const ep_rational& t) {
    ep_bigint q, r;
    ep_bigint::divmod(t.num(), t.den(), q, r);
    const std::string s = q.to_string();
    if (s.size() > 18 || (!s.empty() && s[0] == '-'))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::resolve_l2_projection_constant_sq: cell "
            "index out of the long long range (malformed grid arithmetic)");
    long long v = 0;
    for (std::size_t i = 0; i < s.size(); ++i)
        v = v * 10 + static_cast<long long>(s[i] - '0');
    return v;
}

// (a, b^2) normal form of a 2D canonical key (multiset rule; header note)
inline void b2b_normal_form_2d(const dict_key<2>& key,
                               ep_rational& a, ep_rational& b_sq) {
    const ep_rational* e = key.e.data();
    int im = 0;
    for (int k = 1; k < 3; ++k) if (e[im] < e[k]) im = k;
    ep_rational r1 = e[(im + 1) % 3] / e[im];
    ep_rational r2 = e[(im + 2) % 3] / e[im];
    if (r1 < r2) { const ep_rational t = r1; r1 = r2; r2 = t; }
    // r1 = r_big >= r2 = r_small
    a = (ep_rational(1) + r1 - r2) / ep_rational(2);
    b_sq = r1 - a * a;
}

// theta_min cut: b^2 >= tau * a^2, exact
inline bool b2b_inside_theta_min(const ep_rational& a,
                                 const ep_rational& b_sq) {
    const ep_rational tau(l2_theta_min_tan_sq_num, l2_theta_min_tan_sq_den);
    return !(b_sq < tau * a * a);
}

// cell of a normal-form point (floor + closed-boundary clamp).  The clamp
// is what assigns the top edge b^2 = 3/4 (equilateral corner) to the
// closed top cell -- owner gate addendum 2026-08-18.
inline void b2b_cell_of(const ep_rational& a, const ep_rational& b_sq,
                        long long& ci, long long& cj) {
    const l2_projection_coverage_grid_2d g = l2_projection_coverage_2d_grid();
    const ep_rational x0(g.x0_num, g.x0_den);
    const ep_rational xw(g.xw_num, g.xw_den);
    const ep_rational y0(g.y0_num, g.y0_den);
    const ep_rational yw(g.yw_num, g.yw_den);
    ci = b2b_floor_ll((a - x0) / xw);
    cj = b2b_floor_ll((b_sq - y0) / yw);
    if (ci >= g.nx) ci = g.nx - 1;
    if (cj >= g.ny) cj = g.ny - 1;
}

// coverage table scan for (ci, cj, d)
inline const l2_projection_coverage_entry_2d*
b2b_coverage_find(long long ci, long long cj, int d) {
    int count = 0;
    const l2_projection_coverage_entry_2d* t =
        l2_projection_coverage_2d_table(count);
    for (int k = 0; k < count; ++k)
        if (t[k].ci == ci && t[k].cj == cj && t[k].d == d) return &t[k];
    return 0;
}

// exact-key scan over ONE registry-format table
template <int D, typename E>
inline const E* b2b_key_scan(const E* table, int count,
                             const dict_key<D>& key, int d) {
    for (int k = 0; k < count; ++k) {
        if (table[k].d != d) continue;
        if (dict_key_matches<D>(key, table[k])) return &table[k];
    }
    return 0;
}

// exact-key lookup at degree dd across registry THEN ondemand ledger
template <int D>
inline const typename dict_registry_of<D>::entry_type*
b2b_exact_find(const dict_key<D>& key, int dd);

// NOTE: the table pointer is fetched in its OWN statement before the scan
// call -- table(n) writes the count into n, and folding both into one call
// expression would read n at an unspecified time (argument evaluation
// order), silently scanning 0 entries.
template <>
inline const l2_projection_registry_entry_2d*
b2b_exact_find<2>(const dict_key<2>& key, int dd) {
    int n = 0;
    const l2_projection_registry_entry_2d* t =
        l2_projection_registry_2d_table(n);
    const l2_projection_registry_entry_2d* e = b2b_key_scan<2>(t, n, key, dd);
    if (e) return e;
    n = 0;
    t = l2_projection_ondemand_2d_table(n);
    return b2b_key_scan<2>(t, n, key, dd);
}

template <>
inline const l2_projection_registry_entry_3d*
b2b_exact_find<3>(const dict_key<3>& key, int dd) {
    int n = 0;
    const l2_projection_registry_entry_3d* t =
        l2_projection_registry_3d_table(n);
    const l2_projection_registry_entry_3d* e = b2b_key_scan<3>(t, n, key, dd);
    if (e) return e;
    n = 0;
    t = l2_projection_ondemand_3d_table(n);
    return b2b_key_scan<3>(t, n, key, dd);
}

// ---------------------------------------------------------------------------
// 12-significant-digit UPWARD decimal string of an exact rational in
// (0, 1), pure bigint arithmetic (no floating point, no decimal literal).
// The B2bp formatter produced the same strings from doubles; this one is
// scalar-generic and exact by construction (ceil in the last digit).
// ---------------------------------------------------------------------------
inline std::string b2b_upward_decimal_12(const ep_rational& q) {
    const ep_rational one(1);
    if (!(ep_rational(0) < q))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
            "nonpositive value where a positive bound was expected");
    if (!(q < one))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
            "Chat^2 >= 1 (unexpected)");
    const ep_rational ten(10);
    int k = 0;
    ep_rational t = q;
    while (t < one) {
        t = t * ten;
        ++k;
        if (k > 400)
            vcp::throw_error<vcp::invalid_argument>(
                "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
                "value too small for the decimal formatter");
    }
    // t = q * 10^k in [1, 10); M = ceil(q * 10^(k+11)) in [10^11, 10^12]
    const ep_rational u =
        t * ep_rational(ep_bigint(100000000000LL), ep_bigint(1));
    ep_bigint m, r;
    ep_bigint::divmod(u.num(), u.den(), m, r);
    if (!(r == ep_bigint(0))) m = m + ep_bigint(1);
    std::string digits = m.to_string();
    if (digits.size() == 13) {          // carry: M == 10^12
        digits = digits.substr(0, 12);  // "100000000000"
        --k;
        if (k == 0)
            vcp::throw_error<vcp::invalid_argument>(
                "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
                "Chat^2 rounds up to 1 (unexpected)");
    }
    if (digits.size() != 12)
        vcp::throw_error<vcp::verification_error>(
            "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
            "mantissa is not 12 digits (formatter defect)");
    std::string out = "0.";
    for (int z = 0; z < k - 1; ++z) out += '0';
    out += digits;
    // defense in depth: the parsed string must dominate q exactly
    const ep_rational parsed = constant_from_string<ep_rational>::get(
        out.c_str());
    if (parsed < q)
        vcp::throw_error<vcp::verification_error>(
            "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
            "upward decimal string below the exact value (formatter defect)");
    return out;
}

// long long of a bigint (guarded; ledger keys must fit the entry struct)
inline long long b2b_ll_of_bigint(const ep_bigint& b) {
    const std::string s = b.to_string();
    if (s.size() > 18)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::ondemand_l2_projection_constant_sq: key "
            "fraction does not fit the long long entry field");
    long long v = 0;
    std::size_t i = 0;
    bool neg = false;
    if (i < s.size() && s[i] == '-') { neg = true; ++i; }
    for (; i < s.size(); ++i) v = v * 10 + static_cast<long long>(s[i] - '0');
    return neg ? -v : v;
}

} // namespace detail

// ---------------------------------------------------------------------------
// resolve_l2_projection_constant_sq (design 2)
// ---------------------------------------------------------------------------
template <int D, typename T>
l2_projection_resolution<D, T>
resolve_l2_projection_constant_sq(const T vertices[D + 1][D], int d) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::resolve_l2_projection_constant_sq: "
                  "D must be 2 or 3");
    detail::interval_scalar_contract<T>::require();
    if (d < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::resolve_l2_projection_constant_sq: "
            "projection degree d must be >= 0 (got ", d, ")");

    std::array<std::array<detail::ep_rational, D>, D + 1> rv;
    for (int i = 0; i <= D; ++i)
        for (int c = 0; c < D; ++c)
            rv[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] =
                detail::rational_of_point_interval<T>(vertices[i][c]);
    const detail::dict_key<D> key = detail::dict_key_of<D>(rv);

    l2_projection_resolution<D, T> res;
    res.h_sq = detail::interval_of_exact_rational<T>(
        detail::dict_max_edge_sq<D>(rv));

    // (a, b^2) normal form and coverage-domain flags (D = 2 only)
    detail::ep_rational nf_a, nf_bsq;
    bool in_theta = false;
    long long ci = -1, cj = -1;
    if (D == 2) {
        detail::dict_key<2> k2;
        for (int k = 0; k < 3; ++k)
            k2.e[static_cast<std::size_t>(k)] =
                key.e[static_cast<std::size_t>(k)];
        detail::b2b_normal_form_2d(k2, nf_a, nf_bsq);
        in_theta = detail::b2b_inside_theta_min(nf_a, nf_bsq);
        if (in_theta) detail::b2b_cell_of(nf_a, nf_bsq, ci, cj);
    }

    // (1) exact key at d: registry, then ondemand ledger
    {
        const typename detail::dict_registry_of<D>::entry_type* e =
            detail::b2b_exact_find<D>(key, d);
        if (e) {
            res.ok = true;
            res.source = l2_source::registry_exact;
            res.served_d = d;
            res.cd_sq_over_h_sq =
                constant_from_string<T>::get(e->cd_sq_over_h_sq);
            return res;
        }
    }

    // (2) coverage cell at d (D = 2, inside the theta_min cut)
    if (D == 2 && in_theta) {
        const l2_projection_coverage_entry_2d* e =
            detail::b2b_coverage_find(ci, cj, d);
        if (e) {
            res.ok = true;
            res.source = l2_source::coverage_cell;
            res.served_d = d;
            res.cell_i = ci;
            res.cell_j = cj;
            res.cd_sq_over_h_sq =
                constant_from_string<T>::get(e->cd_sq_over_h_sq);
            return res;
        }
    }

    // (3) degree envelope: the largest d' < d held by (1) or (2)
    for (int dd = d - 1; dd >= 0; --dd) {
        const typename detail::dict_registry_of<D>::entry_type* e =
            detail::b2b_exact_find<D>(key, dd);
        if (e) {
            res.ok = true;
            res.source = l2_source::degree_envelope;
            res.served_d = dd;
            res.cd_sq_over_h_sq =
                constant_from_string<T>::get(e->cd_sq_over_h_sq);
            return res;
        }
        if (D == 2 && in_theta) {
            const l2_projection_coverage_entry_2d* c =
                detail::b2b_coverage_find(ci, cj, dd);
            if (c) {
                res.ok = true;
                res.source = l2_source::degree_envelope;
                res.served_d = dd;
                res.cell_i = ci;
                res.cell_j = cj;
                res.cd_sq_over_h_sq =
                    constant_from_string<T>::get(c->cd_sq_over_h_sq);
                return res;
            }
        }
    }

    // (4) P^0 closed form (D = 2; valid for every d since C_d <= C_0)
    if (D == 2) {
        std::array<T, 2> p0, p1, p2;
        for (int c = 0; c < 2; ++c) {
            p0[static_cast<std::size_t>(c)] = vertices[0][c];
            p1[static_cast<std::size_t>(c)] = vertices[1][c];
            p2[static_cast<std::size_t>(c)] = vertices[2][c];
        }
        T c0 = l2_projection_element_bound(p0, p1, p2);
        {
            const T c1 = l2_projection_element_bound(p1, p2, p0);
            if (c1.upper() < c0.upper()) c0 = c1;
            const T c2 = l2_projection_element_bound(p2, p0, p1);
            if (c2.upper() < c0.upper()) c0 = c2;
        }
        res.ok = true;
        res.source = l2_source::p0_closed_form;
        res.served_d = 0;
        res.cd_sq_over_h_sq = (c0 * c0) / res.h_sq;
        return res;
    }

    // (5) none
    return res;
}

// ---------------------------------------------------------------------------
// ondemand computation (design 3).  Explicit engine invocation; returns the
// value AND the ledger entry text.  NO FILE IS WRITTEN HERE -- appending to
// dict/ondemand_2d.hpp / ondemand_3d.hpp is the caller's explicit act
// (probe constb2b_ondemand_append.cpp).
//
// level == 0 selects the production schedule of the generated registry
// (B2bp): 2D -> L = 3; 3D -> L = 3 (d <= 2), L = 2 (3 <= d <= 6),
// L = 1 (d >= 7).  The achieved alpha = ch^2 / mu is recorded in the
// provenance comment (12-digit upward, information only).
// ---------------------------------------------------------------------------
template <typename T>
struct l2_ondemand_result {
    bool ok;
    T cd_sq_over_h_sq;          // the value the ledger entry will serve
    T chat_sq_raw_upper;        // point interval of the pre-format binary
                                // upper end (for INDEPENDENT re-checks:
                                // string parsed exactly must dominate it)
    int d;
    int level;
    std::string decimal_string; // 12-digit upward Chat^2 string
    std::string entry_text;     // registry-format ledger entry block

    l2_ondemand_result()
        : ok(false), cd_sq_over_h_sq(), chat_sq_raw_upper(), d(-1),
          level(-1), decimal_string(), entry_text() {}
};

template <int D, typename T>
l2_ondemand_result<T>
ondemand_l2_projection_constant_sq(const T vertices[D + 1][D], int d,
                                   int level = 0) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
                  "D must be 2 or 3");
    detail::interval_scalar_contract<T>::require();
    if (level == 0)
        level = (D == 2) ? 3 : (d <= 2 ? 3 : (d <= 6 ? 2 : 1));
    if (level < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::ondemand_l2_projection_constant_sq: "
            "refinement level must be >= 1 (or 0 for the default schedule)");

    // exact key and exact h^2 of the queried simplex
    std::array<std::array<detail::ep_rational, D>, D + 1> rv;
    for (int i = 0; i <= D; ++i)
        for (int c = 0; c < D; ++c)
            rv[static_cast<std::size_t>(i)][static_cast<std::size_t>(c)] =
                detail::rational_of_point_interval<T>(vertices[i][c]);
    const detail::dict_key<D> key = detail::dict_key_of<D>(rv);
    const detail::ep_rational h_sq_exact = detail::dict_max_edge_sq<D>(rv);

    // the engine run (the ONLY heavy step; explicit by design)
    const element_projection_result<D, T> er =
        l2_projection_element_sq_bound<D, T>(vertices, d, level);

    // Chat^2 = cd_sq_upper / h^2, upper end exactified through the shared
    // dyadic rationalizer, then formatted 12-digit upward
    const T h_sq_i = detail::interval_of_exact_rational<T>(h_sq_exact);
    const T chat_sq = er.cd_sq_upper / h_sq_i;
    const detail::ep_rational chat_up_exact =
        detail::rational_of_point_interval<T>(T(chat_sq.upper()));
    l2_ondemand_result<T> out;
    out.d = d;
    out.level = level;
    out.chat_sq_raw_upper = T(chat_sq.upper());
    out.decimal_string = detail::b2b_upward_decimal_12(chat_up_exact);
    out.cd_sq_over_h_sq =
        constant_from_string<T>::get(out.decimal_string.c_str());

    // ledger entry text (registry format; provenance = ondemand + host +
    // date + level + achieved alpha, all information-only comments)
    std::string txt;
    {
        char buf[64];
        std::string keystr = "[";
        std::string nums = "{ ", dens = "{ ";
        for (std::size_t k = 0;
             k < static_cast<std::size_t>(
                     detail::dict_key<D>::edge_count); ++k) {
            const detail::ep_rational& e = key.e[k];
            if (k) { keystr += " "; nums += ", "; dens += ", "; }
            keystr += e.to_string();
            std::snprintf(buf, sizeof buf, "%lldLL",
                          detail::b2b_ll_of_bigint(e.num()));
            nums += buf;
            std::snprintf(buf, sizeof buf, "%lldLL",
                          detail::b2b_ll_of_bigint(e.den()));
            dens += buf;
        }
        keystr += "]";
        nums += " }";
        dens += " }";

        std::string vstr;
        for (int i = 0; i <= D; ++i) {
            vstr += (i ? " (" : "(");
            for (int c = 0; c < D; ++c) {
                if (c) vstr += ",";
                vstr += rv[static_cast<std::size_t>(i)]
                          [static_cast<std::size_t>(c)].to_string();
            }
            vstr += ")";
        }

        char host[256];
        if (gethostname(host, sizeof host) != 0)
            std::strcpy(host, "unknown-host");
        char datebuf[32];
        {
            const std::time_t now = std::time(0);
            std::tm tmv;
            localtime_r(&now, &tmv);
            std::snprintf(datebuf, sizeof datebuf, "%04d-%02d-%02d",
                          tmv.tm_year + 1900, tmv.tm_mon + 1, tmv.tm_mday);
        }
        // achieved alpha = ch^2 / mu, 12-digit upward (information only;
        // alpha >= 1 is reported verbally, not formatted)
        std::string alpha_str;
        {
            const T alpha = er.ch_sq / er.mu_upper;
            const detail::ep_rational alpha_up =
                detail::rational_of_point_interval<T>(T(alpha.upper()));
            if (alpha_up < detail::ep_rational(1))
                alpha_str = detail::b2b_upward_decimal_12(alpha_up);
            else
                alpha_str = "1 or more (alpha unmet)";
        }

        char head[128];
        std::snprintf(head, sizeof head,
                      "        // ondemand %dd key=", D);
        txt += head;
        txt += keystr;
        txt += "\n        //   vertices: ";
        txt += vstr;
        txt += "\n        // source: CONST-B2b ondemand on ";
        txt += host;
        txt += " ";
        txt += datebuf;
        std::snprintf(head, sizeof head,
                      ", L=%d, d=%d, alpha-achieved<=", level, d);
        txt += head;
        txt += alpha_str;
        txt += " (info)\n        { ";
        txt += nums;
        txt += ",\n          ";
        txt += dens;
        std::snprintf(head, sizeof head, ",\n          %d, \"", d);
        txt += head;
        txt += out.decimal_string;
        txt += "\" },\n";
    }
    out.entry_text = txt;
    out.ok = true;
    return out;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DICT_RESOLVE_HPP
