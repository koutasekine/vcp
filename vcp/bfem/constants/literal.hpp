// vcp/bfem/constants/literal.hpp
//
// CONST-A: conversion traits from the authorized decimal-string table of the
// constants layer to the scalar type T a computation actually uses.  The
// table itself lives in poisson_constants.hpp (the VCP_CONSTANTS_TABLE block);
// this header holds NO value of its own -- it is the pure string-to-type
// conversion layer, shared by the poisson side today and by the stokes side
// once that header grows a table of its own (CONST-A design, section three).
//
// The three routes (owner ruling R15 of the CONST-A design):
//
//     kv::interval<X>   kv's own string constructor, whose two ends are built
//                       by rop<X>::fromstring_down / fromstring_up, i.e. the
//                       result is an OUTWARD rounded enclosure of the decimal
//                       (verified against the kv source in CONST-A phase P0);
//     rational          exact parse: the digit string over a power of ten,
//                       reduced by the rational constructor itself.  Accepted
//                       grammar: one or more digits, one dot, one or more
//                       digits -- nothing else.  Sign, exponent, spaces and
//                       every other character throw std::invalid_argument,
//                       because the table format is fixed and anything beyond
//                       it in a table string is a defect, not an input;
//     double            UPWARD point value: the upper end of the outward kv
//                       enclosure, so the point value keeps its upper-bound
//                       status.  (The alternative of forbidding the point
//                       instantiation was considered and not adopted -- R15.)
//
// Lexical regime of THIS header: the traditional one -- no decimal literal
// anywhere, neither as a bare token nor inside a string (the digits below are
// produced by character arithmetic against '0').  The relaxed regime of R16
// applies to the TABLE block of poisson_constants.hpp only.
//
// Authority: sandbox/docs/design/CONST-A_design_v1.0.md.

#ifndef VCP_BFEM_CONSTANTS_LITERAL_HPP
#define VCP_BFEM_CONSTANTS_LITERAL_HPP

#include <string>
#include <stdexcept>
#include <cstddef>

#include <kv/interval.hpp>
// rop<double>: the directed-rounding fromstring_down / fromstring_up pair.
// Without it the generic rop<T> falls back to istringstream (round to
// nearest), which would silently break the outward-enclosure contract of the
// interval route and the upper-bound contract of the double route.
#include <kv/rdouble.hpp>

#include <vcp/bfem/rational.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// constant_from_string<T>::get(s): convert an authorized table string to T.
// The primary template is deliberately left undefined: a scalar type without
// an explicit route below is a compile error at the point of use, never a
// silent conversion through some default.
// ---------------------------------------------------------------------------
template <typename T> struct constant_from_string;

// ---- route one: kv interval, outward rounded enclosure of the decimal ----
template <typename X> struct constant_from_string< kv::interval<X> > {
    static kv::interval<X> get(const char* s) {
        return kv::interval<X>(s);
    }
};

// ---- route two: exact rational --------------------------------------------
template <> struct constant_from_string<rational> {
    static rational get(const char* s) {
        if (s == 0)
            throw std::invalid_argument(
                "vcp::bfem::constants::constant_from_string<rational>: "
                "null string");
        const std::string str(s);
        std::size_t dot = std::string::npos;
        for (std::size_t i = 0; i < str.size(); ++i) {
            const char c = str[i];
            if (c == '.') {
                if (dot != std::string::npos)
                    throw std::invalid_argument(
                        "vcp::bfem::constants::constant_from_string<rational>: "
                        "more than one dot");
                dot = i;
            } else if (!(c >= '0' && c <= '9')) {
                throw std::invalid_argument(
                    "vcp::bfem::constants::constant_from_string<rational>: "
                    "character outside the fixed table grammar "
                    "(digits, one dot, digits)");
            }
        }
        if (dot == std::string::npos || dot == 0 || dot + 1 == str.size())
            throw std::invalid_argument(
                "vcp::bfem::constants::constant_from_string<rational>: "
                "the fixed table grammar is digits, one dot, digits");

        // numerator: the digit string with the dot removed; denominator: ten
        // to the number of fractional digits.  rational(bigint, bigint)
        // normalizes (reduces by the gcd), so the result is canonical.
        // bigint is FULLY qualified: an including header may well have a
        // detail namespace of its own inside vcp::bfem::constants (poisson
        // does), and the unqualified name would resolve to that one first.
        const std::string digits = str.substr(0, dot) + str.substr(dot + 1);
        const std::size_t frac = str.size() - dot - 1;
        vcp::bfem::detail::bigint num =
            vcp::bfem::detail::bigint::from_string(digits);
        vcp::bfem::detail::bigint den(1);
        const vcp::bfem::detail::bigint ten(10);
        for (std::size_t k = 0; k < frac; ++k) den = den * ten;
        return rational(num, den);
    }
};

// ---- route three: upward point value (double) -----------------------------
template <> struct constant_from_string<double> {
    static double get(const char* s) {
        return kv::interval<double>(s).upper();
    }
};

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_LITERAL_HPP
