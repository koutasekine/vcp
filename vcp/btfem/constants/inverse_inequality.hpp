// vcp/btfem/constants/inverse_inequality.hpp
//
// BT-3: the inverse inequality quantity lambda_h of the fully-discrete
// Petrov-Galerkin method of the linear heat equation (constants paper (18)):
//
//     lambda_h := max_{v_h in V_h \ {0}}  || v_h ||^2_{H^1_0} / || v_h ||^2_{L^2}
//
// Conforms to: BT-3 design v1.0 (sections 1, 2.2, 3.1) and
//              BT-3 requirements v1.0 (S5, S6, S8, S9, FR-2, FR-6).
//
// WHAT IS COMPUTED.  V_h is the homogeneous Dirichlet reduction of the
// bfem conforming P^n space fe_space<D, T>(Th, n) (hence the _h01 mark, bfem
// VER-1 R23).  lambda_h is the largest eigenvalue of the reduced
// generalized problem
//
//     K_x v = lambda U_x v        K_x: stiffness, U_x: mass (both reduced)
//
// and is NOT a constant: it is a computed spectral quantity of the mesh and
// degree at hand.  The path is the one of the bfem kappa_h^2
// (poisson_constants.hpp assemble_dense_materials / hypercircle_kappa_core,
// bfem CM-1R / CONST-F), call for call:
//
//   fe_space -> dirichlet_reduction -> sparse assembly + finalize ->
//   sparse reduce + finalize -> densify (matrix<T, DP> = spmatrix) ->
//   vcp::compsym on both -> eigsymge(K, U, E) -> diag(E) ->
//   maximum of the upper ends -> point interval.
//
// No new densification, symmetrisation or eigenvalue route is introduced
// here; the accumulation order is the bfem assembly order and the
// computation is sequential (same input, single thread: bit reproducible;
// DP differences are not contracted, bfem CM-1T).
//
// RETURN VALUE.  As for hypercircle_kappa_h01_squared, the returned T is the
// POINT interval T(lambda_up) with lambda_up the largest upper end of the
// eigsymge enclosures.  Its lower end is therefore NOT a lower bound of
// lambda_h (bfem CM-1S U2 note); only the upper end carries meaning.
//
// SCOPE (bfem four-way notation): "whole mesh, homogeneous Dirichlet".
//
// REFERENCE VALUES (tests): on the reference triangle / tetrahedron with
// n = 1 and no reduction the largest generalized eigenvalue is 36 / 80
// (BT-3 requirements section 1; proved in exact rational arithmetic by gate
// G5 of btfem_const_lambda_tests.cpp).
//
// LEXICAL POLICY (bfem CONST-A, BT-3 NF-2): no decimal literal, no
// `double` / `float` token; the base type is reached through T::base_type.
//
// EXCEPTIONS (bfem three-way ledger, design 2.2):
//   n < 1                                   vcp::invalid_argument
//   reduced dimension 0 (max over empty)    vcp::invalid_argument
//   compsym empty intersection,             vcp::verification_error
//   eigsymge certification failure,
//   largest upper end negative
//   degenerate element                      vcp::bfem::degenerate_element
//                                           (propagated from the space
//                                           construction; outside the
//                                           vcp::error hierarchy)
//
// Authority: sandbox/docs/design/BT-3_design_v1.0.md,
//            sandbox/docs/design/BT-3_requirements_v1.0.md.

#ifndef VCP_BTFEM_CONSTANTS_INVERSE_INEQUALITY_HPP
#define VCP_BTFEM_CONSTANTS_INVERSE_INEQUALITY_HPP

#include <vector>
#include <cstddef>

#include <kv/interval.hpp>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>     // vcp::compsym (must follow matrix.hpp; bfem order)
#include <vcp/spmatrix.hpp>
#include <vcp/imats.hpp>
#include <vcp/spimats.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/constants/poisson_constants.hpp>   // detail::interval_scalar_contract

namespace vcp {
namespace btfem {
namespace constants {

// ---------------------------------------------------------------------------
// inverse_inequality_h01_bound (design 2.2 / 3.1): upper bound of lambda_h.
//   Th : bfem mesh, D in {2, 3}
//   n  : degree of the conforming P^n space (n >= 1)
// ---------------------------------------------------------------------------
template <int D, typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
T inverse_inequality_h01_bound(const vcp::bfem::mesh<D, T>& Th, int n) {
    static_assert(D == 2 || D == 3,
                  "vcp::btfem::constants::inverse_inequality_h01_bound: only D == 2 or D == 3");
    typedef typename T::base_type B;
    vcp::bfem::constants::detail::interval_scalar_contract<T>::require();
    if (n < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::inverse_inequality_h01_bound: n must be >= 1 (got ",
            n, ")");

    // 2. space and homogeneous Dirichlet reduction (bfem assemble_dense_materials)
    vcp::bfem::fe_space<D, T, DP, SP> fs(Th, n);
    const int nd = fs.ndof(n);
    vcp::bfem::dirichlet_reduction<T, DP, SP> dr(nd, fs.dofs(n).boundary_dofs());
    const int ni = dr.reduced_size();
    if (ni <= 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::btfem::constants::inverse_inequality_h01_bound: the reduced space is {0} "
            "(nd = ", nd, "); lambda_h is not defined");

    // 3. sparse assembly, then finalize
    vcp::spmatrix<T, SP> S = fs.stiffness(n);
    vcp::spmatrix<T, SP> M = fs.mixed_mass(n, n);
    S.finalize();
    M.finalize();

    // 4. sparse reduction, then finalize
    vcp::spmatrix<T, SP> Sr = dr.reduce(S);
    vcp::spmatrix<T, SP> Mr = dr.reduce(M);
    Sr.finalize();
    Mr.finalize();

    // 5. densify (matrix<T, DP>::operator=(const spmatrix&), bfem CM-1R R2)
    vcp::matrix<T, DP> K;
    vcp::matrix<T, DP> U;
    K = Sr;
    U = Mr;

    // 6. exact symmetrisation by intersection (bfem CM-1R R7); an empty
    //    intersection raises vcp::verification_error inside compsym
    vcp::compsym(K);
    vcp::compsym(U);

    // 7. enclosure of every generalized eigenvalue; U is SPD and eigsymge's
    //    own certification throws vcp::verification_error on failure
    vcp::matrix<T, DP> E;
    eigsymge(K, U, E);
    vcp::matrix<T, DP> ev = diag(E);

    B lambda_up = ev(0, 0).upper();
    for (int i = 1; i < ni; ++i) {
        if (lambda_up < ev(i, 0).upper()) lambda_up = ev(i, 0).upper();
    }
    if (lambda_up < B(0))
        vcp::throw_error<vcp::verification_error>(
            "vcp::btfem::constants::inverse_inequality_h01_bound: the largest upper end is "
            "negative (", lambda_up, "); the enclosure is not usable");

    return T(lambda_up);
}

} // namespace constants
} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_CONSTANTS_INVERSE_INEQUALITY_HPP
