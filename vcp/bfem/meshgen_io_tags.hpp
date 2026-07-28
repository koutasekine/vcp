// vcp/bfem/meshgen_io_tags.hpp
// MG-3 (B): meshgen_io_traits specializations and type tags (design v1.4,
// sections 5.1/5.2).
//
//  - vcp::bfem::rational, tag "vbrat": numerator and denominator each as
//    sign i8 + nlimbs u32 + u32 limb list (bigint::limb(i) written
//    directly; reading rebuilds the bigint by base-2^32 Horner because
//    bigint has no public limb-insertion API, F4).
//  - kv::interval<F>, tag "itv(" + inner tag + ")": inf then sup as two F
//    value records; reconstruction through the rounding-free two-endpoint
//    constructor interval<F>(inf, sup) (F9).
//
// This is the single header of the track exempted from the type-name
// lexical ban (design 5.4); arithmetic decimal literals remain forbidden.

#ifndef VCP_BFEM_MESHGEN_IO_TAGS_HPP
#define VCP_BFEM_MESHGEN_IO_TAGS_HPP

#include <string>
#include <vector>

#include <kv/interval.hpp>

#include <vcp/bfem/meshgen_io.hpp>
#include <vcp/bfem/rational.hpp>

namespace vcp {
namespace bfem {

template <>
struct meshgen_io_traits<rational> {
    static std::string tag() { return std::string("vbrat"); }

    static void write_bigint(meshgen_io_detail::owriter& w,
                             const detail::bigint& b) {
        w.i8(b.is_zero() ? 0 : (b.negative() ? -1 : 1));
        const int n = b.num_limbs();
        w.u32(static_cast<std::uint32_t>(n));
        for (int i = 0; i < n; ++i) w.u32(b.limb(i));
    }

    static detail::bigint read_bigint(meshgen_io_detail::ireader& r) {
        const int s = r.i8();
        const std::uint32_t n = r.u32();
        if (!(s == 0 || s == 1 || s == -1) || ((s == 0) != (n == 0)))
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt bigint record (sign)");
        std::vector<std::uint32_t> limbs;
        limbs.reserve(n);
        for (std::uint32_t i = 0; i < n; ++i) limbs.push_back(r.u32());
        if (n != 0 && limbs[n - 1] == 0u)
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt bigint record (leading zero)");
        const detail::bigint base(4294967296ll);   // 2^32
        detail::bigint b(0);
        for (std::uint32_t i = n; i > 0; --i)
            b = b * base +
                detail::bigint(static_cast<long long>(limbs[i - 1]));
        if (s < 0) b = -b;
        return b;
    }

    static void write_value(meshgen_io_detail::owriter& w, const rational& x) {
        write_bigint(w, x.num());
        write_bigint(w, x.den());
    }

    static rational read_value(meshgen_io_detail::ireader& r,
                               bool verify_roundtrip) {
        detail::bigint num = read_bigint(r);
        detail::bigint den = read_bigint(r);
        if (den.is_zero() || den.negative())
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt rational record (denominator)");
        rational v(num, den);
        if (verify_roundtrip) {
            // records are written from normalized rationals; a record that
            // renormalizes to different observers is corrupt
            if (!(v.num() == num) || !(v.den() == den))
                throw meshgen_io_error(
                    "vcp::bfem::meshgen_io: verify_roundtrip mismatch "
                    "(non-normalized rational record)");
        }
        return v;
    }
};

template <typename F>
struct meshgen_io_traits<kv::interval<F> > {
    static std::string tag() {
        return std::string("itv(") + meshgen_io_traits<F>::tag() +
               std::string(")");
    }

    static void write_value(meshgen_io_detail::owriter& w,
                            const kv::interval<F>& x) {
        meshgen_io_traits<F>::write_value(w, x.lower());
        meshgen_io_traits<F>::write_value(w, x.upper());
    }

    static kv::interval<F> read_value(meshgen_io_detail::ireader& r,
                                      bool verify_roundtrip) {
        F lo = meshgen_io_traits<F>::read_value(r, verify_roundtrip);
        F hi = meshgen_io_traits<F>::read_value(r, verify_roundtrip);
        if (hi < lo)
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt interval record (reversed "
                "endpoints)");
        return kv::interval<F>(lo, hi);   // F9: no rounding
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_IO_TAGS_HPP
