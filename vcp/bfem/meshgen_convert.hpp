// vcp/bfem/meshgen_convert.hpp
// MG-3 (A): mesh type conversion (design v1.4).
//
// This layer holds ONLY coordinate loops and type dispatch; every numeric
// conversion is delegated (single-source principle):
//  (i)  TS = vcp::bfem::rational  -> convert_traits<TD>::from_rational
//       (enclose-once contract, convert_traits.hpp L7-14; never restated
//       or reimplemented here),
//  (ii) pairs covered by vcp::convert (vcp/vcp_converter.hpp: the
//       {dd, mpfr<N>, and double via its guard} point and kv::interval
//       families, C9/C10) -> vcp::convert(src, dst). Interval targets keep
//       enclosure by vcp_converter's per-endpoint directed rounding;
//       interval -> point is vcp_converter's mid semantics (C11),
//  (iii) anything else -> fallback TD(x) (faithfulness not claimed).
// Detection of (ii) is C++11 expression SFINAE on vcp::convert(const TS&,
// TD&). The kv headers below are included BEFORE vcp_converter.hpp so that
// its guard-conditional overload set is fully enabled (directive v1.4).
//
// The base-2 mantissa decomposition detail (F8) stays here from Phase 1:
// it is the internal mechanism of the exact IO (meshgen_io.hpp) and of the
// exact-rational referee used by the tests (C6: it is NOT a public
// conversion path; the v1.3 decompose+rebuild FP->FP path is NOT used).
//
// Lexical policy (design section 5.4): no decimal literals, no floating
// point type tokens, no sqrt/abs/min/max, integer literals only.

#ifndef VCP_BFEM_MESHGEN_CONVERT_HPP
#define VCP_BFEM_MESHGEN_CONVERT_HPP

#include <array>
#include <vector>
#include <stdexcept>
#include <utility>
#include <type_traits>
#include <limits>

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>
#include <kv/dd.hpp>
#include <kv/rdd.hpp>
#include <kv/mpfr.hpp>
#include <kv/rmpfr.hpp>

#include <vcp/vcp_converter.hpp>

#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/mesh.hpp>

namespace vcp {
namespace bfem {

class meshgen_convert_error : public std::runtime_error {
public:
    explicit meshgen_convert_error(const char* msg) : std::runtime_error(msg) {}
};

namespace meshgen_convert_detail {

// Result of the base-2 decomposition of a nonzero finite value x:
//   x = sign * (bits[0].bits[1]bits[2]... in base 2) * 2^exp,
// bits[0] == 1 (normalized leading bit), bits.size() == nbits, and the last
// stored bit is 1 (the expansion is emitted until the remainder is exactly
// zero). sign == 0 encodes x == 0 (exp == 0, bits empty).
struct fp_decomp {
    int sign;
    long long exp;
    std::vector<int> bits;
    fp_decomp() : sign(0), exp(0), bits() {}
};

// F8 algorithm. Allowed operations on T: *T(2), /T(2), binary -, <, ==,
// construction from small ints. Exact for base-2 floating point T of any
// mantissa length (the mantissa length is discovered, not assumed).
// NaN / non-finite inputs are rejected (x==x and x-x==T(0) probes) so the
// scaling loops below cannot run unbounded on them; the exponent guards
// bound the loops for any remaining pathological input.
template <typename T>
fp_decomp decompose_base2(const T& x, long long nbits_limit = 65536) {
    fp_decomp d;
    if (x == T(0)) return d;
    if (!(x == x))
        throw meshgen_convert_error(
            "vcp::bfem::meshgen_convert: NaN cannot be decomposed");
    if (!(x - x == T(0)))
        throw meshgen_convert_error(
            "vcp::bfem::meshgen_convert: non-finite value cannot be decomposed");

    const T zero = T(0);
    const T one = T(1);
    const T two = T(2);

    T y = x;
    if (x < zero) {
        d.sign = -1;
        y = zero - x;
    } else {
        d.sign = 1;
    }

    // normalize to 1 <= y < 2, tracking the exponent of the leading bit.
    // Both directions are exact in base-2 arithmetic (doubling below 1
    // cannot overflow; halving at or above 2 cannot enter the subnormal
    // range). The guard bounds the loop far beyond any realistic exponent.
    const long long exp_guard = 16777216;
    long long e = 0;
    long long guard = 0;
    while (y < one) {
        y = y * two;
        e = e - 1;
        guard = guard + 1;
        if (guard > exp_guard)
            throw meshgen_convert_error(
                "vcp::bfem::meshgen_convert: exponent guard exceeded (small)");
    }
    guard = 0;
    while (!(y < two)) {
        y = y / two;
        e = e + 1;
        guard = guard + 1;
        if (guard > exp_guard)
            throw meshgen_convert_error(
                "vcp::bfem::meshgen_convert: exponent guard exceeded (large)");
    }
    d.exp = e;

    // emit mantissa bits until the remainder is exactly zero. Invariant:
    // r in [0, 1), so r*2 < 2 (no overflow) and the subtraction of one from
    // a value in [1, 2) is exact (Sterbenz regime; F8 measured this for
    // long double / kv::dd / kv::mpfr as well).
    d.bits.push_back(1);
    T r = y - one;
    while (!(r == zero)) {
        if (static_cast<long long>(d.bits.size()) >= nbits_limit)
            throw meshgen_convert_error(
                "vcp::bfem::meshgen_convert: mantissa expansion did not "
                "terminate within the nbits cap (non-base-2 scalar?)");
        r = r * two;
        if (r < one) {
            d.bits.push_back(0);
        } else {
            d.bits.push_back(1);
            r = r - one;
        }
    }
    return d;
}

// Exact rational value of a decomposition (referee / IO detail, C6: not a
// public conversion). Numerator by bigint Horner over the bit list, then
// the scale 2^k applied ONE BIT AT A TIME on bigint (shl1) -- unbounded
// integers cannot overflow, and no scale power is ever formed in a
// floating point type (rule 4 / F8 trap applies to T, not to bigint).
inline rational decomp_to_rational(const fp_decomp& d) {
    if (d.sign == 0) return rational();
    detail::bigint m(0);
    const detail::bigint two(2);
    for (std::size_t i = 0; i < d.bits.size(); ++i)
        m = m * two + detail::bigint(static_cast<long long>(d.bits[i]));
    const long long k = d.exp - (static_cast<long long>(d.bits.size()) - 1);
    detail::bigint den(1);
    if (k >= 0) {
        for (long long i = 0; i < k; ++i) m.shl1();
    } else {
        for (long long i = 0; i < -k; ++i) den.shl1();
    }
    if (d.sign < 0) m = -m;
    return rational(m, den);
}

// convenience referee entry: exact rational value of a base-2 scalar
template <typename T>
rational exact_rational_of(const T& x, long long nbits_limit = 65536) {
    return decomp_to_rational(decompose_base2<T>(x, nbits_limit));
}

// C++11 expression SFINAE: is vcp::convert(const TS&, TD&) callable?
// (Overload ambiguity in that call also yields false here; the dispatch
// expectations are pinned pair-by-pair in the Phase 2 test so that a
// silently changed resolution surfaces as a test failure.)
template <typename TS, typename TD>
struct has_vcp_convert {
    template <typename S, typename D>
    static auto probe(int)
        -> decltype(vcp::convert(std::declval<const S&>(),
                                 std::declval<D&>()),
                    std::true_type());
    template <typename S, typename D>
    static std::false_type probe(...);
    static const bool value = decltype(probe<TS, TD>(0))::value;
};

// The one builtin arithmetic scalar that vcp_converter targets (binary64;
// identified through numeric_limits so no type-name token appears here).
template <typename T>
struct is_binary64_point {
    static const bool value =
        std::numeric_limits<T>::is_specialized &&
        !std::numeric_limits<T>::is_integer &&
        std::numeric_limits<T>::radix == 2 &&
        std::numeric_limits<T>::digits == 53;
};

// C12 gate (measured necessity): vcp_converter's int overloads (e.g.
// convert(const int&, D&)) are reachable from ANY arithmetic source via a
// narrowing standard conversion, so raw callability would route
// float/long double through an int truncation. C12 mandates the TD(x)
// fallback for those sources; therefore an arithmetic TS enters the
// vcp::convert branch only when it is the binary64 scalar itself.
template <typename TS, typename TD>
struct use_vcp_convert {
    // is_binary64_point must not be instantiated for non-arithmetic T:
    // kv specializes std::numeric_limits for its scalars without the full
    // member set (measured), so the gate is selected lazily via conditional.
    typedef typename std::conditional<std::is_arithmetic<TS>::value,
                                      is_binary64_point<TS>,
                                      std::true_type>::type source_gate;
    static const bool value =
        has_vcp_convert<TS, TD>::value && source_gate::value;
};

// (ii)/(iii) dispatch below the rational specialization.
template <typename TS, typename TD, bool HasVcpConvert>
struct meshgen_converter_impl;

template <typename TS, typename TD>
struct meshgen_converter_impl<TS, TD, true> {
    static TD apply(const TS& x) {
        TD y;
        vcp::convert(x, y);
        return y;
    }
};

template <typename TS, typename TD>
struct meshgen_converter_impl<TS, TD, false> {
    static TD apply(const TS& x) { return TD(x); }
};

} // namespace meshgen_convert_detail

// per-coordinate converter, priority (design section 4.5-2):
//  (i) rational source -> convert_traits (partial specialization below)
//  (ii) vcp::convert pair -> vcp::convert
//  (iii) otherwise -> TD(x)
template <typename TS, typename TD>
struct meshgen_converter {
    static TD convert(const TS& x) {
        return meshgen_convert_detail::meshgen_converter_impl<
            TS, TD,
            meshgen_convert_detail::use_vcp_convert<TS, TD>::value>::apply(x);
    }
};

template <typename TD>
struct meshgen_converter<rational, TD> {
    static TD convert(const rational& x) {
        return convert_traits<TD>::from_rational(x.num(), x.den());
    }
};

// ---- public API (design section 3) ----

// elements are shared, so only coordinates are converted
template <typename TD, int D, typename TS>
void convert_mesh_lists(const std::vector<std::array<TS, D> >& vin,
                        std::vector<std::array<TD, D> >& vout) {
    vout.clear();
    vout.reserve(vin.size());
    for (std::size_t v = 0; v < vin.size(); ++v) {
        std::array<TD, D> p;
        for (int d = 0; d < D; ++d)
            p[static_cast<std::size_t>(d)] =
                meshgen_converter<TS, TD>::convert(
                    vin[v][static_cast<std::size_t>(d)]);
        vout.push_back(p);
    }
}

template <typename TD, int D, typename TS>
mesh<D, TD> convert_mesh(const mesh<D, TS>& m) {
    std::vector<std::array<TS, D> > vin;
    vin.reserve(static_cast<std::size_t>(m.num_vertices()));
    for (int v = 0; v < m.num_vertices(); ++v)
        vin.push_back(m.vertex(v));
    std::vector<std::array<TD, D> > vout;
    convert_mesh_lists<TD, D, TS>(vin, vout);
    std::vector<std::array<int, D + 1> > elems;
    elems.reserve(static_cast<std::size_t>(m.num_elements()));
    for (int e = 0; e < m.num_elements(); ++e)
        elems.push_back(m.element(e));
    return mesh<D, TD>::from_lists(vout, elems);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_CONVERT_HPP
