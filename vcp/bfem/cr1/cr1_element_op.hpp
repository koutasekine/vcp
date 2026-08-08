// vcp/bfem/cr1/cr1_element_op.hpp
// CR1-1 (P1 Crouzeix-Raviart nonconforming element, D = 2, 3): the element
// layer -- local basis, local stiffness, local mass.  Closed forms only; this
// header contains NO numerical quadrature (design 2, mandatory).
//
// Conforms to: CR1-1 design v1.0 (sections 2, 3.2) and
//              CR1-1 implementation directive v1.0 (phase P2).
//
// ---------------------------------------------------------------------------
// Facet parameterization convention (normative, single definition point)
// ---------------------------------------------------------------------------
// The CR1 degree of freedom on a facet e is the PARAMETER-REPRESENTATION mean
// (design 0, ruling R1 -- this is the master definition, not the midpoint /
// centroid value):
//
//   D = 2 (edge, vertices p_a, p_b):
//       dof_e(v) = int_0^1 v(x(t)) dt,        x(t) = (1 - t) p_a + t p_b
//   D = 3 (face, vertices p_a, p_b, p_c):
//       dof_e(v) = 2 int_{t1 + t2 <= 1, t >= 0} v(x(t1, t2)) dt1 dt2,
//                                             x = (1 - t1 - t2) p_a + t1 p_b + t2 p_c
//
// The parameterization vertices p_a, p_b (, p_c) are ALWAYS taken in ASCENDING
// GLOBAL VERTEX ORDER.  That ordered tuple has exactly one source in this
// library -- cr1_dofmap<D>::facet_vertices(f) -- and cr1_space<D,...>
// republishes it, so the assembly, the interpolation entry points and the
// gate tests all read the same data.  A facet mean is invariant under the
// re-parameterization anyway (that is why the family tag is pn_family_tag);
// the order is fixed only so that the parameterization is REPRODUCIBLE.
//
// The parameter mean equals the arc-length / area mean int_e v ds / |e|
// because the parameterization is affine (constant Jacobian), but the
// parameter form carries no square root and therefore stays inside an exact
// rational scalar T.  That is why it is the master form.
//
// ---------------------------------------------------------------------------
// Local basis  [self-derived; verified exactly by gate R-G2]
// ---------------------------------------------------------------------------
//   phi_i = 1 - D lambda_i,  lambda_i the barycentric coordinate of the vertex
//   OPPOSITE the facet that carries local index i.
// On the facet e_j one has lambda_j = 0, and for i != j the restriction of
// lambda_i to e_j is a barycentric coordinate of that (D-1)-simplex, whose
// facet mean is 1/D.  Hence dof_{e_j}(phi_i) = 1 - D * (1/D) = 0 for i != j
// and 1 - 0 = 1 for i == j, i.e. dof_{e_j}(phi_i) = delta_ij.  Also
// sum_i phi_i = (D + 1) - D sum_i lambda_i = 1 (partition of unity).
//
// ---------------------------------------------------------------------------
// Local stiffness (design 2)
// ---------------------------------------------------------------------------
//   grad phi_i = -D grad lambda_i is constant on the element, so
//   A^K_ij = int_K grad phi_i . grad phi_j = |K| D^2 (grad lambda_i . grad lambda_j).
// grad lambda and |K| come from the frozen element_geometry<D,T> layer
// (vcp/bfem/geometry.hpp): |K| = |det B| / D! is formed by a SIGN decision
// followed by a negation, never by a square root, and the single division of
// the element is 1/det -- exactly the contract this design asks for.
//
// ---------------------------------------------------------------------------
// Local mass (closed form, no quadrature)  [self-derived; gates R-G4/R-G6]
// ---------------------------------------------------------------------------
// Barycentric monomial formula (design 2; the 3D instance is equation (32) of
// Liu-Nakao-Oishi 2022, the general D form is [source not verbatim-verified]
// and is machine-checked here by R-G4/R-G5/R-G6):
//
//   int_K lambda_0^{a_0} ... lambda_D^{a_D} dx = D! |K| (prod_i a_i!) / (D + sum_i a_i)!
//
// hence  int_K lambda_i = |K| / (D + 1)  and
//        int_K lambda_i lambda_j = |K| D! (1 + delta_ij) / (D + 2)!
//                                = |K| (1 + delta_ij) / F,   F := (D + 1)(D + 2).
// Expanding phi_i phi_j = 1 - D lambda_i - D lambda_j + D^2 lambda_i lambda_j:
//
//   M^K_ij = |K| * [ F - 2 D (D + 2) + D^2 (1 + delta_ij) ] / F
//
// which is an EXACT rational multiple of |K| depending only on (D, i == j).
// Instances: D = 2 gives |K|/3 on the diagonal and 0 off it (the 2D CR mass
// matrix is diagonal); D = 3 gives 2|K|/5 and -|K|/20.  Row sum:
// sum_j M^K_ij = |K| / (D + 1) = int_K phi_i, the identity gate R-G6 checks.

#ifndef VCP_BFEM_CR1_CR1_ELEMENT_OP_HPP
#define VCP_BFEM_CR1_CR1_ELEMENT_OP_HPP

#include <array>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/cr1/cr1_dofmap.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// cr1_local_matrix<T, N>: the (D + 1) x (D + 1) element block.  It exposes the
// operator()(i, j) the frozen Y2 scatter kernels of vcp/bfem/dofmap.hpp call,
// so no adapter sits between the element layer and the assembly.
// ---------------------------------------------------------------------------
template <typename T, int N>
struct cr1_local_matrix {
    std::array<T, static_cast<std::size_t>(N) * static_cast<std::size_t>(N)> a;

    const T& operator()(int i, int j) const {
        assert(i >= 0 && i < N && j >= 0 && j < N);
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(N)
                 + static_cast<std::size_t>(j)];
    }
    T& operator()(int i, int j) {
        assert(i >= 0 && i < N && j >= 0 && j < N);
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(N)
                 + static_cast<std::size_t>(j)];
    }
};

// ---------------------------------------------------------------------------
// cr1_element_op<D, T>: stateless -- the element geometry is an argument, not
// owned state (the CR1 element block needs nothing beyond |K| and grad lambda,
// so there is no scratch worth caching between elements).
// Both blocks are built on the upper triangle and MIRRORED, so A(i, j) and
// A(j, i) are bit identical for every scalar T (the assembled matrix then
// satisfies the strict is_symmetric() test of the dense policies).
// ---------------------------------------------------------------------------
template <int D, typename T>
struct cr1_element_op {
    static_assert(D == 2 || D == 3,
                  "bfem::detail::cr1_element_op: only D == 2 or D == 3");

    typedef cr1_local_matrix<T, D + 1> local_matrix_type;

    // A^K_ij = |K| D^2 (grad lambda_i . grad lambda_j)
    static void local_stiffness(const element_geometry<D, T>& g,
                                local_matrix_type& loc) {
        const T w = g.measure() * T(D * D);
        for (int i = 0; i <= D; ++i) {
            for (int j = i; j <= D; ++j) {
                T dot(0);
                for (int d = 0; d < D; ++d)
                    dot += g.grad_lambda(i, d) * g.grad_lambda(j, d);
                loc(i, j) = w * dot;
                if (j != i) loc(j, i) = loc(i, j);      // mirrored, not recomputed
            }
        }
    }

    // M^K_ij = |K| * mass_coefficient(i == j)
    static void local_mass(const element_geometry<D, T>& g,
                           local_matrix_type& loc) {
        const T diag = g.measure() * mass_coefficient(true);
        const T off = g.measure() * mass_coefficient(false);
        for (int i = 0; i <= D; ++i)
            for (int j = i; j <= D; ++j) {
                loc(i, j) = (i == j) ? diag : off;
                if (j != i) loc(j, i) = loc(i, j);
            }
    }

    // value of the local basis at the barycentric point lam: phi_i = 1 - D lam_i
    static T basis_value(int i, const T* lam) {
        assert(i >= 0 && i <= D);
        return T(1) - T(D) * lam[i];
    }

    // F = (D + 1)(D + 2)
    static long long mass_denominator() {
        return static_cast<long long>(D + 1) * static_cast<long long>(D + 2);
    }
    // F - 2 D (D + 2) + D^2 (1 + delta_ij)
    static long long mass_numerator(bool diagonal) {
        const long long m = diagonal ? 2 : 1;
        return mass_denominator()
             - 2 * static_cast<long long>(D) * static_cast<long long>(D + 2)
             + static_cast<long long>(D) * static_cast<long long>(D) * m;
    }
    // enclose-once conversion of the exact rational coefficient (one magic
    // static per (D, T); no per-element division is added)
    static const T& mass_coefficient(bool diagonal) {
        static const T off = rational_to<T>(mass_numerator(false), mass_denominator());
        static const T dia = rational_to<T>(mass_numerator(true), mass_denominator());
        return diagonal ? dia : off;
    }
};

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CR1_CR1_ELEMENT_OP_HPP
