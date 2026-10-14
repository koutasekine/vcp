// vcp/btfem/constants/heat_pg_constants.hpp
//
// BT-3: the three error constants C~_1, C~_0, c~_0 of the fully-discrete
// Petrov-Galerkin approximation of the linear heat equation
// d_t u - nu Laplace u = f, u(a) = 0 (constants paper (11) with
// u_0 = u^_0 = 0, BT-3 requirements S1), assembled from the two projection
// error constants C_Omega(h) (space, bfem) and C_tau(q - 1) (time, this
// layer) and the inverse inequality quantity lambda_h.
//
// Conforms to: BT-3 design v1.0 (sections 1, 2.3, 2.4, 3.3) and
//              BT-3 requirements v1.0 (S1, S2, S6, S7, S8, S9, FR-3, FR-4).
//
// FORMULAS (BT-3 requirements 2.1; the constants paper's Theorems 1, 2, 3
// and Lemma 3 with u_0 = u^_0 = 0; k = q, C_tau := C_tau(k - 1)):
//
//     K  := ( C_Omega^2 + nu^2 lambda_h C_tau^2 + 4 nu C_tau ) / nu^2   Theorem 1
//     L  := 4 ( C_Omega^2 / nu + C_tau )                                Theorem 2
//     C~_1 = sqrt(K)        || u - u_{h,tau}^k ||_{L^2 H^1_0}  <= C~_1 || f ||
//     c~_0 = sqrt(L)        || (u - u_{h,tau}^k)(b) ||_{L^2}   <= c~_0 || f ||
//     C~_0 = 4 C_tau + 2 C_Omega C~_1                               Theorem 3
//                           || u - u_{h,tau}^k ||_{L^2 L^2}  <= C~_0 || f ||
//            (Theorem 3 with Lemma 3, || d_t (u - u_{h,tau}^k) || <= 2 || f ||,
//             and the vanishing initial term)
//
// SQUARED CANON (bfem R17, requirements S7): K = C~_1^2 and L = c~_0^2 are
// the squared quantities and are formed from C_Omega^2 WITHOUT any square
// root; the only square roots are sqrt(K), sqrt(L) and sqrt(C_Omega^2)
// inside C~_0.  The evaluation order is FIXED (design 2.3):
//
//     1. K = (c_omega_sq + nu*nu*lambda_h*c_tau*c_tau + 4*nu*c_tau) / (nu*nu)
//     2. L = 4 * (c_omega_sq / nu + c_tau)
//     3. c1_tilde = sqrt(K),  c0_small = sqrt(L)
//     4. c0_tilde = 4*c_tau + 2*sqrt(c_omega_sq)*c1_tilde
//
// MATERIALS.  The space side (C_Omega(h)^2 and lambda_h) depends on the mesh
// and the degree only, so heat_pg_space_materials_h01 computes it ONCE into
// heat_pg_space_material_set_h01 and the per time interval quantities
// (tau, q, nu) are combined afterwards without recomputation.
//   C_Omega(h)^2 is bfem ritz_projection_error_constant_h01_sq(Th, n)
//   (2D / 3D, dictionary served, requirements S2; its exceptions propagate,
//   in particular the 3D registry miss, S10).
//   lambda_h is inverse_inequality_h01_bound(Th, n) of this layer.
//   The order is c_omega_sq FIRST, then lambda_h (design 2.3): a dictionary
//   miss is reported before any space of this layer is built.
//
// NO THIN PROJECTION ENTRY (design 2.3 / 8): the constants are read from
// the members of heat_pg_error_constant_set_h01.
//
// _h01 LEDGER (bfem VER-1 R23): every entry point below assumes the
// homogeneous Dirichlet space V_h = fe_space(Th, n) cap H^1_0, hence _h01.
// SCOPE (bfem four-way notation): "whole mesh, homogeneous Dirichlet".
//
// LEXICAL POLICY (bfem CONST-A, BT-3 NF-2): no decimal literal, no
// `double` / `float` token; sqrt and the arithmetic of T only.
//
// EXCEPTIONS (bfem three-way ledger, design 2.3):
//   nu.lower() <= 0, c_tau / lambda_h / c_omega_sq lower < 0, q < 1
//                                         vcp::invalid_argument
//   from the materials                    bfem / inverse_inequality.hpp
//                                         exceptions propagate unchanged
//                                         (vcp::invalid_argument,
//                                          vcp::verification_error,
//                                          vcp::bfem::degenerate_element)
//
// Authority: sandbox/docs/design/BT-3_design_v1.0.md,
//            sandbox/docs/design/BT-3_requirements_v1.0.md.

#ifndef VCP_BTFEM_CONSTANTS_HEAT_PG_CONSTANTS_HPP
#define VCP_BTFEM_CONSTANTS_HEAT_PG_CONSTANTS_HPP

#include <cmath>

#include <kv/interval.hpp>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
#include <vcp/spmatrix.hpp>
#include <vcp/imats.hpp>
#include <vcp/spimats.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/constants/poisson_constants.hpp>   // ritz_projection_error_constant_h01_sq, mesh_resolution_report
#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/constants/time_projection_constants.hpp>
#include <vcp/btfem/constants/inverse_inequality.hpp>

namespace vcp {
namespace btfem {
namespace constants {

// ---------------------------------------------------------------------------
// heat_pg_space_material_set_h01 (design 2.3): the space side materials,
// computed once per (mesh, degree).
// ---------------------------------------------------------------------------
template <typename T>
struct heat_pg_space_material_set_h01 {
    int n;              // space degree
    T   c_omega_sq;     // C_Omega(h)^2 (bfem ritz_projection_error_constant_h01_sq)
    T   lambda_h;       // inverse_inequality_h01_bound (point interval, upper end meaningful)

    heat_pg_space_material_set_h01() : n(0), c_omega_sq(), lambda_h() {}
};

// ---------------------------------------------------------------------------
// heat_pg_error_constant_set_h01 (design 2.3): the four materials and the
// five derived quantities.
// ---------------------------------------------------------------------------
template <typename T>
struct heat_pg_error_constant_set_h01 {
    T c_omega_sq;   // material: C_Omega(h)^2
    T c_tau;        // material: C_tau(q - 1)
    T lambda_h;     // material: lambda_h upper bound
    T nu;           // material: nu
    T K;            // = c1_tilde^2   (Theorem 1)   squared canon
    T L;            // = c0_small^2   (Theorem 2)   squared canon
    T c1_tilde;     // sqrt(K)
    T c0_small;     // sqrt(L)
    T c0_tilde;     // 4 c_tau + 2 sqrt(c_omega_sq) c1_tilde   (Theorem 3 + Lemma 3)

    heat_pg_error_constant_set_h01()
        : c_omega_sq(), c_tau(), lambda_h(), nu(), K(), L(), c1_tilde(), c0_small(), c0_tilde() {}
};

// ---------------------------------------------------------------------------
// heat_pg_space_materials_h01 (design 2.3): C_Omega(h)^2 first, lambda_h second.
// ---------------------------------------------------------------------------
template <int D, typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
heat_pg_space_material_set_h01<T>
heat_pg_space_materials_h01(const vcp::bfem::mesh<D, T>& Th, int n,
                            vcp::bfem::constants::mesh_resolution_report* rep = nullptr) {
    static_assert(D == 2 || D == 3,
                  "vcp::btfem::constants::heat_pg_space_materials_h01: only D == 2 or D == 3");
    vcp::bfem::constants::detail::interval_scalar_contract<T>::require();
    if (n < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::heat_pg_space_materials_h01: n must be >= 1 (got ", n, ")");
    heat_pg_space_material_set_h01<T> sm;
    sm.n = n;
    sm.c_omega_sq = vcp::bfem::constants::ritz_projection_error_constant_h01_sq<T, DP, SP>(Th, n, rep);
    sm.lambda_h   = inverse_inequality_h01_bound<D, T, DP, SP>(Th, n);
    return sm;
}

// ---------------------------------------------------------------------------
// heat_pg_error_constants_h01, formula only (design 2.3): the materials are
// whatever the caller supplies (so they can be replaced), the evaluation
// order is the fixed 1-4 above.
// ---------------------------------------------------------------------------
template <typename T>
heat_pg_error_constant_set_h01<T>
heat_pg_error_constants_h01(const T& c_omega_sq, const T& c_tau,
                            const T& lambda_h, const T& nu) {
    using std::sqrt;
    vcp::bfem::constants::detail::interval_scalar_contract<T>::require();
    if (!(nu.lower() > 0))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::heat_pg_error_constants_h01: nu must have a positive lower end (got ",
            nu.lower(), ")");
    if (c_tau.lower() < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::heat_pg_error_constants_h01: c_tau must be nonnegative (lower end ",
            c_tau.lower(), ")");
    if (lambda_h.lower() < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::heat_pg_error_constants_h01: lambda_h must be nonnegative (lower end ",
            lambda_h.lower(), ")");
    if (c_omega_sq.lower() < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::heat_pg_error_constants_h01: c_omega_sq must be nonnegative (lower end ",
            c_omega_sq.lower(), ")");

    heat_pg_error_constant_set_h01<T> r;
    r.c_omega_sq = c_omega_sq;
    r.c_tau      = c_tau;
    r.lambda_h   = lambda_h;
    r.nu         = nu;
    // 1. K
    r.K = (c_omega_sq + nu * nu * lambda_h * c_tau * c_tau + T(4) * nu * c_tau) / (nu * nu);
    // 2. L
    r.L = T(4) * (c_omega_sq / nu + c_tau);
    // 3. square roots
    r.c1_tilde = sqrt(r.K);
    r.c0_small = sqrt(r.L);
    // 4. C~_0
    r.c0_tilde = T(4) * c_tau + T(2) * sqrt(c_omega_sq) * r.c1_tilde;
    return r;
}

// ---------------------------------------------------------------------------
// materials + time side (design 2.3): c_tau = time_projection_error_constant(q - 1, tau)
// ---------------------------------------------------------------------------
template <typename T>
heat_pg_error_constant_set_h01<T>
heat_pg_error_constants_h01(const heat_pg_space_material_set_h01<T>& sm,
                            const T& tau, int q, const T& nu) {
    if (q < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::heat_pg_error_constants_h01: q must be >= 1 (got ", q, ")");
    const T c_tau = time_projection_error_constant<T>(q - 1, tau);
    return heat_pg_error_constants_h01<T>(sm.c_omega_sq, c_tau, sm.lambda_h, nu);
}

template <typename T>
heat_pg_error_constant_set_h01<T>
heat_pg_error_constants_h01(const heat_pg_space_material_set_h01<T>& sm,
                            const vcp::btfem::time_mesh<T>& tm, int q, const T& nu) {
    return heat_pg_error_constants_h01<T>(sm, tm.width(), q, nu);
}

// ---------------------------------------------------------------------------
// one-shot convenience (design 2.3): builds the material set and calls the
// material form; same path, hence the same bits as the two-step use.
// ---------------------------------------------------------------------------
template <int D, typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
heat_pg_error_constant_set_h01<T>
heat_pg_error_constants_h01(const vcp::bfem::mesh<D, T>& Th, int n,
                            const T& tau, int q, const T& nu,
                            vcp::bfem::constants::mesh_resolution_report* rep = nullptr) {
    const heat_pg_space_material_set_h01<T> sm = heat_pg_space_materials_h01<D, T, DP, SP>(Th, n, rep);
    return heat_pg_error_constants_h01<T>(sm, tau, q, nu);
}

} // namespace constants
} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_CONSTANTS_HEAT_PG_CONSTANTS_HPP
