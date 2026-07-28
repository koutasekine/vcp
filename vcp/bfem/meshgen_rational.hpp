// vcp/bfem/meshgen_rational.hpp
// MG-3 (A0): public alias for the exact rational scalar of the meshgen /
// conversion / exact-IO track (design MG-3 v1.1, section 3, ruling R1), plus
// a string parser convenience built on bigint::from_string.
//
// Lexical policy (design section 5.4): no decimal literals, no floating
// point type tokens, integer literals only.

#ifndef VCP_BFEM_MESHGEN_RATIONAL_HPP
#define VCP_BFEM_MESHGEN_RATIONAL_HPP

#include <string>

#include <vcp/bfem/rational.hpp>

namespace vcp {
namespace bfem {

// R1: the canonical exact rational type of the exact meshing pipeline.
// (rational.hpp already exposes the same alias; restated here so that this
// header is the documented public entry point of the MG-3 track.)
typedef detail::rational rational;

// Parses "num/den" or "num" (decimal, optional leading '+'/'-' on either
// part). Malformed input raises std::invalid_argument (from bigint parsing
// or the structural checks below). A zero denominator is not a format
// error; it surfaces as the std::logic_error thrown by rational itself.
inline rational rational_from_string(const std::string& s) {
    const std::string::size_type slash = s.find('/');
    if (slash == std::string::npos)
        return rational(detail::bigint::from_string(s), detail::bigint(1));
    if (s.find('/', slash + 1) != std::string::npos)
        throw std::invalid_argument(
            "vcp::bfem::rational_from_string: more than one '/'");
    if (slash == 0 || slash + 1 == s.size())
        throw std::invalid_argument(
            "vcp::bfem::rational_from_string: empty numerator or denominator");
    detail::bigint num = detail::bigint::from_string(s.substr(0, slash));
    detail::bigint den = detail::bigint::from_string(s.substr(slash + 1));
    return rational(num, den);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_RATIONAL_HPP
