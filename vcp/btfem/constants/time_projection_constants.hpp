// vcp/btfem/constants/time_projection_constants.hpp
//
// BT-3: the time direction L^2 projection error constant C_tau(d) of the
// fully-discrete Petrov-Galerkin method of the linear heat equation.
//
// Conforms to: BT-3 design v1.0 (sections 1, 2.1, 3.2) and
//              BT-3 requirements v1.0 (S3, S4, S6, FR-1, NF-2).
//
// WHAT IS COMPUTED.  With P_d the L^2(J_i) projection onto the polynomials of
// degree d on one sub-interval J_i of width tau (constants paper (5)), the
// constant of
//
//     || v - P_d v ||_{L^2(J_i)} <= C_tau(d) || v' ||_{L^2(J_i)}   for all v in H^1(J_i)
//
// (constants paper (6), stated per sub-interval; summing the squares over
// the sub-intervals gives the J-wide form with tau = max |J_i|) is returned as
//
//     C_tau(d) = c_d * tau
//
// where c_d is the error constant of the ONE DIMENSIONAL H^1_0 projection
// onto the polynomials of degree d + 1 on the unit interval:
//
//     || u' - (pi_{d+1} u)' ||_{L^2(0,1)} <= c_d || u'' ||_{L^2(0,1)}.
//
// REDUCTION [self derived; BT-3 requirements 2.2, checked by the rational
// gate G1 of btfem_const_time_tests.cpp]: for v in H^1(0, tau) put
// V(t) := int_0^t v, let R V be the degree d + 1 polynomial that interpolates V
// at both ends and satisfies (V - RV)' orthogonal to {q' : q in P_{d+1},
// q(0) = q(tau) = 0}.  Then (V - RV)' = (I - P_d) v, and u := V - (linear
// interpolant of V) lies in H^1_0 cap H^2 with u - pi_{d+1} u = V - RV and
// |u|_{H^2} = || v' ||.  Hence the H^1_0 projection constant of degree d + 1,
// scaled by tau, is an upper bound of C_tau(d).
//
// VALUES (BT-3 requirements 2.2, table):
//
//   d      c_d                                     source               kind
//   0      1 / pi                                  NYK 1998 section 1   best
//                                                  (Schultz, Spline Analysis,
//                                                   Theorem 2.5 (2.17))
//   1      1 / (2 pi)                              NYK 1998 Theorem 1   best
//   2      (1 / (3 sqrt(5) pi^3))                  NYK 1998 Theorem 2   upper
//          * sqrt( -2160 + 25 pi^4                 equation (14)        bound
//                  + 4 sqrt(5) sqrt( -174960 - 1080 pi^4
//                                    + 243 pi^6 + 5 pi^8 ) )
//          (approximately 1 / 8.92337)
//   >= 3   sqrt( max(L1, L2) )                     KY 1999 Theorem 2.3  upper
//          L1 = 1 / (2 (2d+1)(2d+5))               equation (2.3), N = d bound
//               + 1 / (4 (2d+5) sqrt(2d+3) sqrt(2d+7))
//          L2 = 1 / (4 (2d+5) sqrt(2d+3) sqrt(2d+7))
//               + 1 / (2 (2d+5)(2d+9))
//               + 1 / (4 (2d+9) sqrt(2d+7) sqrt(2d+11))
//
//   NYK 1998 = M. T. Nakao, N. Yamamoto, S. Kimura, J. Approx. Theory 93
//              (1998) 491-500.  Their C_N is the constant of the projection
//              onto piecewise polynomials of degree N + 1, so N = d.
//   KY 1999  = S. Kimura, N. Yamamoto, Bull. Inform. Cybernet. 31 (1999)
//              109-115.  Q_N(e) has degree N + 1, so N = d; (2.3) is the
//              same term order as vcp/ldbase.hpp Ritz_projection_error(N).
//   Index correspondence: KY N = d, ldbase N = d, Nakao-Kinoshita (2009)
//   Theorem 2 N = d + 1 (the latter is used by the tests as a cross check
//   only and is not part of this header).
//
// _h01 LEDGER (bfem VER-1 R23): C_tau carries no boundary condition on the
// time interval (the H^1_0 structure appears only inside the reduction
// above), so the entry point has NO _h01 mark -- the same reading as the
// bfem l2_projection_error_constant.
//
// SCOPE (bfem four-way notation): "reference interval, scaled by tau".
//
// LEXICAL POLICY (bfem CONST-A, BT-3 requirements NF-2): no decimal literal,
// no `double` / `float` token; the only non-exact operations are the
// arithmetic of T, sqrt and kv::constants<T>::pi().  pi is taken ONCE and
// its powers are formed by multiplication (no pow), which keeps the
// interval over-estimation of the d = 2 branch small (design 3.2).
//
// TYPE CONTRACT: T is restricted to the kv interval types by
// vcp::bfem::constants::detail::interval_scalar_contract<T> (a point type
// would evaluate pi inexactly and the upper bound property would be lost).
//
// EXCEPTIONS (bfem three-way ledger):
//   d < 0                                      vcp::invalid_argument
//   tau.lower() <= 0                           vcp::invalid_argument
//   d == 2 and a sqrt argument has lower <= 0  vcp::verification_error
//                                              (transcription error detector;
//                                              the arguments are about 898.4
//                                              and 543.3 in exact arithmetic)
//
// Authority: sandbox/docs/design/BT-3_design_v1.0.md,
//            sandbox/docs/design/BT-3_requirements_v1.0.md.

#ifndef VCP_BTFEM_CONSTANTS_TIME_PROJECTION_CONSTANTS_HPP
#define VCP_BTFEM_CONSTANTS_TIME_PROJECTION_CONSTANTS_HPP

#include <algorithm>
#include <cmath>

#include <kv/interval.hpp>
#include <kv/constants.hpp>

#include <vcp/error.hpp>
#include <vcp/bfem/constants/poisson_constants.hpp>   // detail::interval_scalar_contract

namespace vcp {
namespace btfem {
namespace constants {

namespace detail {

// c_d for d >= 3: KY 1999 (2.3) with N = d, in the term order of
// vcp/ldbase.hpp Ritz_projection_error (the three square roots first, L1 as
// two terms, L2 as three terms, the max of the upper ends returned to a
// point interval before the outer sqrt).
template <typename T>
T time_projection_coefficient_ky(int d) {
    using std::sqrt;
    using std::max;
    const T sqrt2N_3  = sqrt(T(2 * d + 3));
    const T sqrt2N_7  = sqrt(T(2 * d + 7));
    const T sqrt2N_11 = sqrt(T(2 * d + 11));

    const T L1 = T(1) / T(2 * (2 * d + 1) * (2 * d + 5))
               + T(1) / (T(4 * (2 * d + 5)) * sqrt2N_3 * sqrt2N_7);
    const T L2 = T(1) / (T(4 * (2 * d + 5)) * sqrt2N_3 * sqrt2N_7)
               + T(1) / T(2 * (2 * d + 5) * (2 * d + 9))
               + T(1) / (T(4 * (2 * d + 9)) * sqrt2N_7 * sqrt2N_11);

    return sqrt(T(max(L1.upper(), L2.upper())));
}

// c_2: NYK 1998 Theorem 2, equation (14).  pi is taken once; the powers are
// products.  Both square root arguments are large-cancellation sums, so
// their lower ends are checked to be positive (a transcription error would
// show up here first).
template <typename T>
T time_projection_coefficient_nyk2() {
    using std::sqrt;
    const T pi  = kv::constants<T>::pi();
    const T pi2 = pi * pi;
    const T pi3 = pi2 * pi;
    const T pi4 = pi2 * pi2;
    const T pi6 = pi4 * pi2;
    const T pi8 = pi4 * pi4;
    const T sqrt5 = sqrt(T(5));

    const T inner = T(-174960) - T(1080) * pi4 + T(243) * pi6 + T(5) * pi8;
    if (!(inner.lower() > 0))
        vcp::throw_error<vcp::verification_error>(
            "vcp::btfem::constants::time_projection_error_constant: d = 2 inner sqrt "
            "argument is not certainly positive (lower end ", inner.lower(), ")");
    const T outer = T(-2160) + T(25) * pi4 + T(4) * sqrt5 * sqrt(inner);
    if (!(outer.lower() > 0))
        vcp::throw_error<vcp::verification_error>(
            "vcp::btfem::constants::time_projection_error_constant: d = 2 outer sqrt "
            "argument is not certainly positive (lower end ", outer.lower(), ")");
    return sqrt(outer) / (T(3) * sqrt5 * pi3);
}

} // namespace detail

// ---------------------------------------------------------------------------
// time_projection_error_constant (design 2.1):  C_tau(d) = c_d * tau.
//   d   : degree of the L^2 projection target (d = q - 1 for the btfem time
//         degree q; the constants paper writes C_tau(k - 1))
//   tau : sub-interval width (interval with positive lower end)
// ---------------------------------------------------------------------------
template <typename T>
T time_projection_error_constant(int d, const T& tau) {
    vcp::bfem::constants::detail::interval_scalar_contract<T>::require();
    if (d < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::time_projection_error_constant: d must be >= 0 (got ",
            d, ")");
    if (!(tau.lower() > 0))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::time_projection_error_constant: tau must have a "
            "positive lower end (got ", tau.lower(), ")");

    T c;
    if (d == 0) {
        c = T(1) / kv::constants<T>::pi();
    } else if (d == 1) {
        c = T(1) / (T(2) * kv::constants<T>::pi());
    } else if (d == 2) {
        c = detail::time_projection_coefficient_nyk2<T>();
    } else {
        c = detail::time_projection_coefficient_ky<T>(d);
    }
    return c * tau;
}

} // namespace constants
} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_CONSTANTS_TIME_PROJECTION_CONSTANTS_HPP
