// vcp/bfem/convert_traits.hpp
// Layer 0: enclose-once conversion rational -> T (U2).
//
// Conforms to: L0 external design v0.3 (sections 5.1, 5.2, 5.3) and
//              L0 internal design v0.3 (section 8).
//
// Contract (external design 5.1, mandatory):
//  - interval T: the result MUST enclose the exact rational value. The default
//    implementation below achieves this using only the inclusion property of
//    the interval arithmetic of T (Horner over limbs, then one division).
//    Any change that breaks this inductive argument (e.g. pointwise
//    intermediate values) is forbidden for any optimization purpose.
//  - point T: the result is an approximation (faithful, a few ulp as a
//    quality goal); Layer 0 claims no exactness for point types.
//
// This is the only place in Layer 0 where a division of T appears.

#ifndef VCP_BFEM_CONVERT_TRAITS_HPP
#define VCP_BFEM_CONVERT_TRAITS_HPP

#include <stdexcept>

#include <vcp/bfem/rational.hpp>

namespace vcp {
namespace bfem {

namespace detail {

// Evaluate |x| (then sign) on T by Horner over base 2^32 limbs.
// Each limb is fed as T(hi16) * T(65536) + T(lo16) so that only T(int) with
// arguments < 65536 is required (no T(uint64) constructor is assumed; the
// values are exactly representable in any T with at least 16 mantissa bits,
// and for interval T every step is inclusion preserving).
template <typename T>
T horner_limbs(const bigint& x) {
    const T sixty_five_536 = T(65536);            // 2^16, exact in T
    const T base = sixty_five_536 * sixty_five_536; // 2^32 built on T
    T r = T(0);
    for (int i = x.num_limbs() - 1; i >= 0; --i) {
        unsigned long v = static_cast<unsigned long>(x.limb(i));
        int hi = static_cast<int>(v >> 16);
        int lo = static_cast<int>(v & 0xFFFFu);
        T limb_t = T(hi) * sixty_five_536 + T(lo);
        r = r * base + limb_t;
    }
    if (x.negative()) r = -r;
    return r;
}

} // namespace detail

template <typename T>
struct convert_traits {
    // default implementation: Horner over limbs on T, one division num/den.
    // If the arithmetic of T is interval arithmetic the result is
    // automatically an enclosure of num/den (external design 5.1).
    static T from_rational(const detail::bigint& num, const detail::bigint& den) {
        if (den.is_zero())
            throw std::logic_error("bfem::convert_traits: zero denominator");
        return detail::horner_limbs<T>(num) / detail::horner_limbs<T>(den);
    }
};

// Public helper (v0.3 addition, external design 5.3): converts a small exact
// rational constant with the enclose-once contract, without exposing
// detail::bigint at the public boundary.
template <typename T>
T rational_to(long long num, long long den) {
    if (den == 0)
        throw std::invalid_argument("bfem::rational_to: zero denominator");
    return convert_traits<T>::from_rational(detail::bigint(num), detail::bigint(den));
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONVERT_TRAITS_HPP
