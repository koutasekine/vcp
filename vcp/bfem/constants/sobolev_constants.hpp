// vcp/bfem/constants/sobolev_constants.hpp
//
// Verified upper bounds of the Sobolev embedding constant C_p of
//
//     || u ||_{L^p(Omega)} <= C_p || u ||_{H^1_0(Omega)},   u in H^1_0(Omega)
//
// by two independent routes, plus their min selection.  Common ground layer
// for the Poisson (2D/3D) and Stokes (componentwise) consumers.
//
// Primary source (CONST-D design section zero; every formula below was
// checked verbatim against the paper before implementation):
//
//     K. Tanaka, K. Sekine, M. Mizuguchi, S. Oishi, "Sharp numerical
//     inclusion of the best constant for embedding H^1_0(Omega) -> L^p(Omega)
//     on bounded convex domain", Journal of Computational and Applied
//     Mathematics 311 (2017) 306--313, Appendix (Theorem A.1 = the
//     Aubin--Talenti best constant T_p, Corollary A.2 = the bounded-domain
//     measure route, Theorem A.3 = Plum's spectral route).
//
// The two routes:
//
//     (A.2)  C_p = |Omega|^{(2-q)/(2q)} T_p,  q = np/(n+p): needs the domain
//            measure only, NO spectral information.  T_p is the Gamma/pi
//            expression of Theorem A.1, evaluated through kv's certified
//            gamma and pow with all exponents built as exact fractions.
//     (A.3)  from a verified LOWER bound rho_lb of the minimal spectral
//            point rho of -Delta on H^1_0(Omega) and a weight sigma >= 0
//            of the inner product (grad., grad.) + sigma (., .): C_p is
//            monotone decreasing in rho, so evaluating at a lower bound
//            yields an upper bound.  This is the FIRST entry point of the
//            layer honouring the sigma-weighted inner product (R12).
//
// Both are valid upper bounds, so their min is one as well; sigma > 0 leaves
// (A.2) valid because the sigma-weighted norm dominates || grad u ||, but
// only (A.3) benefits from sigma (design section one).
//
// BC ledger (owner ruling R23: the BC assumption is part of the NAME, the
// tag sits right after the concept name and before the form suffixes):
// every entry point below assumes the homogeneous Dirichlet space H^1_0 and
// carries _h01, except aubin_talenti_constant (a constant of the whole-space
// inequality on W^{1,q}(R^n) -- no boundary condition to tag) and
// domain_measure (pure geometry, BC-free).  H^1(Omega) general (Neumann
// type) embeddings are OUT of scope of this header (future track).
//
// Return value convention: C_p ITSELF is returned, not its square.  R17
// (prefer squared forms) is deliberately NOT applied here: a fractional pow
// is essential in every route, so avoiding the square root buys nothing
// (design section two states this exemption).
//
// The GUARANTEED UPPER BOUND IS THE RETURN VALUE'S .upper(); the lower end
// carries no claim beyond being a valid enclosure end of the computed
// expression evaluated at the supplied arguments.
//
// Lexical policy (traditional regime): no decimal literals, bare or in
// string, no double / float tokens; values never appear in comments (R14).
// sqrt, kv::constants<T>::pi(), kv::gamma and pow(interval, interval) are
// the transcendental vocabulary.  Exponents are composed as exact long long
// fractions FIRST and only then converted through vcp::bfem::rational_to.
//
// Authority: sandbox/docs/design/CONST-D_design_v1.2.md.
#ifndef VCP_BFEM_CONSTANTS_SOBOLEV_CONSTANTS_HPP
#define VCP_BFEM_CONSTANTS_SOBOLEV_CONSTANTS_HPP

#include <array>
#include <cstddef>

#include <kv/interval.hpp>
#include <kv/constants.hpp>
#include <kv/gamma.hpp>

#include <vcp/error.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/convert_traits.hpp>

namespace vcp {
namespace bfem {
namespace constants {

namespace detail {

// scalar contract: T must be a kv::interval-like type exposing T::base_type
// (the same contract as poisson_constants' interval_scalar_contract; gamma
// and pow(interval, interval) are only meaningful on such a type).
template <typename T>
struct sobolev_scalar_contract {
    typedef typename T::base_type base_type;
    static void require() {
        static_assert(sizeof(base_type) > 0,
                      "vcp::bfem::constants (sobolev): T must be a "
                      "kv::interval-like type exposing T::base_type");
    }
};

// exact fraction -> T (enclosure); num may be negative, den must be positive
// in every call site below (the sign lives in the numerator).
template <typename T>
T frac(long long num, long long den) {
    return vcp::bfem::rational_to<T>(num, den);
}

// shared validation of the rational p = pnum / pden (positivity only; the
// route-specific domain conditions live in the entry points).
inline void require_positive_p(long long pnum, long long pden,
                               const char* who) {
    if (pnum <= 0 || pden <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::", who,
            ": p must be a positive rational (pnum = ", pnum,
            ", pden = ", pden, ")");
}

} // namespace detail

// ---------------------------------------------------------------------------
// aubin_talenti_constant<D, T>(pnum, pden): the best constant T_p of the
// whole-space Sobolev inequality on W^{1,q}(R^n) (Theorem A.1, with
// p = nq/(n-q), i.e. q = np/(n+p)):
//
//     T_p = pi^{-1/2} n^{-1/q} ((q-1)/(n-q))^{1-1/q}
//           { Gamma(1+n/2) Gamma(n) / (Gamma(n/q) Gamma(1+n-n/q)) }^{1/n}
//
// Domain: 1 < q < n.  With q = np/(n+p) the upper part q < n holds for every
// positive p, so the one condition is q > 1, i.e. p > n/(n-1); outside it
// vcp::invalid_argument is thrown.  All exponents and gamma arguments are
// exact fractions of (pnum, pden) composed in long long BEFORE conversion.
// ---------------------------------------------------------------------------
template <int D, typename T>
T aubin_talenti_constant(long long pnum, long long pden) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::aubin_talenti_constant: "
                  "D must be 2 or 3");
    detail::sobolev_scalar_contract<T>::require();
    detail::require_positive_p(pnum, pden, "aubin_talenti_constant");

    // q = np/(n+p) as the exact fraction qnum/qden
    const long long n = D;
    const long long qnum = n * pnum;
    const long long qden = n * pden + pnum;

    // q > 1  <=>  p > n/(n-1)  <=>  pnum (n-1) > n pden
    if (!(pnum * (n - 1) > n * pden))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::aubin_talenti_constant: p = ", pnum, "/",
            pden, " is out of domain for n = ", n,
            " (q = np/(n+p) must exceed one, i.e. p > n/(n-1))");

    using std::pow;
    // pi^{-1/2}
    T tp = pow(kv::constants<T>::pi(), detail::frac<T>(-1, 2));
    // n^{-1/q} = n^{-qden/qnum}
    tp = tp * pow(T(static_cast<int>(n)), detail::frac<T>(-qden, qnum));
    // ((q-1)/(n-q))^{1-1/q}: base = (qnum-qden)/(n qden - qnum) exactly
    // (both parts positive: the numerator by q > 1, the denominator equals
    // n^2 pden after cancellation, hence q < n unconditionally)
    tp = tp * pow(detail::frac<T>(qnum - qden, n * qden - qnum),
                  detail::frac<T>(qnum - qden, qnum));
    // Gamma block ^{1/n}: arguments 1 + n/2, n, n/q = n qden / qnum and
    // 1 + n - n/q = ((1+n) qnum - n qden)/qnum, all positive exact fractions
    const T gnum = kv::gamma(detail::frac<T>(n + 2, 2)) *
                   kv::gamma(T(static_cast<int>(n)));
    const T gden = kv::gamma(detail::frac<T>(n * qden, qnum)) *
                   kv::gamma(detail::frac<T>((1 + n) * qnum - n * qden, qnum));
    tp = tp * pow(gnum / gden, detail::frac<T>(1, n));
    return tp;
}

// ---------------------------------------------------------------------------
// sobolev_embedding_constant_h01_measure<D, T>(pnum, pden, mea):
// route (A.2), Corollary A.2:
//
//     C_p = |Omega|^{(2-q)/(2q)} T_p,   q = np/(n+p)
//
// mea is the measure |Omega| of the domain, USER SUPPLIED (see
// domain_measure below for the mesh helper); it must be certainly positive.
// Domain: p in (n/(n-1), 2n/(n-2)] for n >= 3, p in (n/(n-1), infinity) for
// n = 2 -- in particular p = 2 is OUT for n = 2 (the corollary needs q > 1).
// At the critical exponent p = 2n/(n-2) of n >= 3 the measure exponent
// vanishes exactly and the factor is skipped (C_p = T_p, measure free).
// NO spectral information enters this route, and sigma does not either: the
// sigma-weighted H^1_0 norm dominates || grad u ||, so the bound stays valid
// under every sigma >= 0 (design section one; gate G-F5 checks it).
// ---------------------------------------------------------------------------
template <int D, typename T>
T sobolev_embedding_constant_h01_measure(long long pnum, long long pden,
                                         const T& mea) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::sobolev_embedding_constant_h01_"
                  "measure: D must be 2 or 3");
    detail::sobolev_scalar_contract<T>::require();
    detail::require_positive_p(pnum, pden,
                               "sobolev_embedding_constant_h01_measure");

    const long long n = D;
    // upper limit for n >= 3: p <= 2n/(n-2)  <=>  pnum (n-2) <= 2 n pden
    if (n >= 3 && !(pnum * (n - 2) <= 2 * n * pden))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_measure: "
            "p = ", pnum, "/", pden, " exceeds the critical exponent "
            "2n/(n-2) for n = ", n);
    // lower limit p > n/(n-1) (q > 1) is enforced by aubin_talenti_constant;
    // it is checked here first so that the message names THIS entry point
    if (!(pnum * (n - 1) > n * pden))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_measure: "
            "p = ", pnum, "/", pden, " is out of domain for n = ", n,
            " (p > n/(n-1) is required; for n = 2 note that p = 2 is OUT)");
    // failure-side gate: a measure not certainly positive cannot be used
    if (!(mea > T(0)))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_measure: "
            "the domain measure is not certainly positive");

    using std::pow;
    const long long qnum = n * pnum;
    const long long qden = n * pden + pnum;
    // (2-q)/(2q) = (2 qden - qnum)/(2 qnum); nonnegative on the domain
    // (q <= 2), zero exactly at the critical exponent of n >= 3
    const long long enum_ = 2 * qden - qnum;
    T c = aubin_talenti_constant<D, T>(pnum, pden);
    if (enum_ != 0)
        c = pow(mea, detail::frac<T>(enum_, 2 * qnum)) * c;
    return c;
}

// ---------------------------------------------------------------------------
// sobolev_embedding_constant_h01_spectrum<D, T>(pnum, pden, rho_lb, sigma):
// route (A.3), Theorem A.3 (Plum).  rho_lb is a verified LOWER bound of the
// minimal spectral point rho of -Delta on H^1_0(Omega) endowed with the
// inner product (grad., grad.) + sigma (., .): C_p is monotone decreasing
// in rho, so the value at rho_lb is an upper bound for the value at rho.
// The supply of rho_lb is the existing Poincare path (CR1 + CONST-A lower
// bound formula) -- this header only receives the value (layer separation).
//
//     (a) n = 2, p in [2, infinity), nu = floor(p/2):
//         C_p = (1/2)^{1/2 + (2 nu - 3)/p}
//               [ (p/2)(p/2 - 1) ... (p/2 - nu + 2) ]^{2/p}
//               (rho + (p/2) sigma)^{-1/p}
//         (the bracketed product is EMPTY = one when nu = 1; at p = 2 the
//         power-of-one-half exponent vanishes and C_2 = (rho+sigma)^{-1/2},
//         the sigma-weighted Poincare constant -- gate G-F4 crosses it
//         against poincare_constant_h01_sq_bound)
//     (b) n >= 3, p in [2, 2n/(n-2)], s = n(1/p - 1/2 + 1/n) in [0, 1]:
//         C_p = ((n-1)/(sqrt(n)(n-2)))^{1-s} (s/(s rho + sigma))^{s/2}
//         (at s = 0, i.e. the critical exponent, the spectral factor
//         degenerates to one and C_p is the coefficient alone; at s = 1,
//         i.e. p = 2, the coefficient factor degenerates to one -- both
//         endpoints are taken EXACTLY, no pow with a zero exponent is
//         evaluated: gate G-F8)
//
// Throws vcp::invalid_argument on: p out of [2, .) resp. [2, 2n/(n-2)],
// sigma not certainly nonnegative, rho_lb and sigma both not certainly
// positive (rho = 0 requires sigma > 0), or a spectral denominator that is
// not certainly positive.
// ---------------------------------------------------------------------------
template <int D, typename T>
T sobolev_embedding_constant_h01_spectrum(long long pnum, long long pden,
                                          const T& rho_lb, const T& sigma) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::sobolev_embedding_constant_h01_"
                  "spectrum: D must be 2 or 3");
    detail::sobolev_scalar_contract<T>::require();
    detail::require_positive_p(pnum, pden,
                               "sobolev_embedding_constant_h01_spectrum");

    const long long n = D;
    if (!(pnum >= 2 * pden))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_spectrum: "
            "p = ", pnum, "/", pden, " is below two (p >= 2 is required "
            "for every n)");
    if (n >= 3 && !(pnum * (n - 2) <= 2 * n * pden))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_spectrum: "
            "p = ", pnum, "/", pden, " exceeds the critical exponent "
            "2n/(n-2) for n = ", n);
    // failure-side gates on the spectral inputs
    if (!(sigma >= T(0)))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_spectrum: "
            "sigma is not certainly nonnegative");
    if (!(rho_lb > T(0)) && !(sigma > T(0)))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_spectrum: "
            "rho_lb is not certainly positive and sigma is not certainly "
            "positive (rho = 0 requires sigma > 0)");

    using std::pow;
    using std::sqrt;

    if (n == 2) {
        // (a): nu = floor(p/2) >= 1 on the domain
        const long long nu = pnum / (2 * pden);
        T c(1);
        // (1/2)^{1/2 + (2 nu - 3)/p}: exponent (pnum + 2(2 nu - 3) pden)
        // over (2 pnum); it vanishes exactly at p = 2 (nu = 1)
        const long long hnum = pnum + 2 * (2 * nu - 3) * pden;
        if (hnum != 0)
            c = c * pow(detail::frac<T>(1, 2),
                        detail::frac<T>(hnum, 2 * pnum));
        // [ (p/2)(p/2-1)...(p/2-nu+2) ]^{2/p}: nu - 1 factors
        // (pnum - 2 j pden)/(2 pden), j = 0 .. nu-2; empty when nu = 1
        if (nu >= 2) {
            T prod = detail::frac<T>(pnum, 2 * pden);
            for (long long j = 1; j <= nu - 2; ++j)
                prod = prod * detail::frac<T>(pnum - 2 * j * pden, 2 * pden);
            c = c * pow(prod, detail::frac<T>(2 * pden, pnum));
        }
        // (rho + (p/2) sigma)^{-1/p}
        const T base = rho_lb + detail::frac<T>(pnum, 2 * pden) * sigma;
        if (!(base > T(0)))
            vcp::throw_error<vcp::invalid_argument>(
                "vcp::bfem::constants::sobolev_embedding_constant_h01_"
                "spectrum: rho_lb + (p/2) sigma is not certainly positive "
                "(n = 2 route)");
        return c * pow(base, detail::frac<T>(-pden, pnum));
    }

    // (b) n >= 3: s = n(1/p - 1/2 + 1/n) = (2 n pden - (n-2) pnum)/(2 pnum);
    // the domain pins snum in [0, sden]
    const long long snum = 2 * n * pden - (n - 2) * pnum;
    const long long sden = 2 * pnum;
    const T coeff = T(static_cast<int>(n - 1)) /
                    (sqrt(T(static_cast<int>(n))) * T(static_cast<int>(n - 2)));
    if (snum == 0) return coeff;             // critical exponent: s = 0
    T c = (snum == sden)
              ? T(1)                          // p = 2: s = 1, coeff factor out
              : pow(coeff, detail::frac<T>(sden - snum, sden));
    const T s = detail::frac<T>(snum, sden);
    const T den2 = s * rho_lb + sigma;
    if (!(den2 > T(0)))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::sobolev_embedding_constant_h01_spectrum: "
            "s rho_lb + sigma is not certainly positive (n >= 3 route)");
    return c * pow(s / den2, detail::frac<T>(snum, 2 * sden));
}

// ---------------------------------------------------------------------------
// sobolev_embedding_constant_h01<D, T>(pnum, pden, mea, rho_lb, sigma,
// which): the min selection of the two routes.  Both are valid upper
// bounds, so the smaller UPPER end is the sharper guaranteed constant; the
// comparison is made on the upper ends and the selected interval is
// returned UNCHANGED (the l2_projection_element_constant precedent -- a
// total, deterministic order even when the enclosures overlap).  Ties keep
// route (A.2).  which, when non-null, receives zero for route (A.2)
// (measure) and one for route (A.3) (spectrum).
//
// Domain = the INTERSECTION of the two routes' domains: for n = 2 that is
// p in (2, infinity) (p = 2 is out through the (A.2) side), for n >= 3 it
// is p in [2, 2n/(n-2)].  Out-of-domain input throws from the respective
// route with its own message.
// ---------------------------------------------------------------------------
template <int D, typename T>
T sobolev_embedding_constant_h01(long long pnum, long long pden, const T& mea,
                                 const T& rho_lb, const T& sigma,
                                 int* which = nullptr) {
    const T a2 = sobolev_embedding_constant_h01_measure<D, T>(pnum, pden, mea);
    const T a3 = sobolev_embedding_constant_h01_spectrum<D, T>(pnum, pden,
                                                               rho_lb, sigma);
    if (a3.upper() < a2.upper()) {
        if (which) *which = 1;
        return a3;
    }
    if (which) *which = 0;
    return a2;
}

// ---------------------------------------------------------------------------
// domain_measure<D, T>(Th): |Omega| as the sum of the element measures,
// through element_geometry (P0 reconnaissance: measure() is |det|/D!, the
// existing per-element API; this is the one-line composition the design
// left to implementation discretion).  BC-free, pure geometry.  Throws
// vcp::invalid_argument on an empty mesh; a degenerate element throws from
// geometry_traits<T>::sign inside element_geometry.
// ---------------------------------------------------------------------------
template <int D, typename T>
T domain_measure(const vcp::bfem::mesh<D, T>& Th) {
    if (Th.num_elements() <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::domain_measure: mesh has no element");
    T sum(0);
    for (int e = 0; e < Th.num_elements(); ++e) {
        const std::array<int, D + 1>& el = Th.element(e);
        std::array<std::array<T, D>, D + 1> v;
        for (int i = 0; i <= D; ++i)
            v[static_cast<std::size_t>(i)] =
                Th.vertex(el[static_cast<std::size_t>(i)]);
        sum = sum + vcp::bfem::element_geometry<D, T>::from_vertices(v).measure();
    }
    return sum;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_SOBOLEV_CONSTANTS_HPP
