// vcp/bfem/fe_function.hpp
// Layer 3: primal/dual vector types (H9, X4).
//
// Conforms to: L3 external design v0.2 (section 4).
//
// fe_function (primal, canonical state: Bernstein coefficients) is created
// only through the fe_space factories (B-1); dual_vector (moments, output
// only) is created only by fe_space::dual. There is NO dual -> primal
// conversion anywhere in this layer (the mass inverse is banned by type).

#ifndef VCP_BFEM_FE_FUNCTION_HPP
#define VCP_BFEM_FE_FUNCTION_HPP

#include <utility>

#include <vcp/matrix.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class fe_space;

template <int D, typename T, typename P = vcp::mats<T> >
class fe_function {
public:
    int degree() const { return deg_; }
    const vcp::matrix<T, P>& coeffs() const { return c_; }
    vcp::matrix<T, P>&       coeffs() { return c_; }        // direct edit allowed

private:
    explicit fe_function(int m) : deg_(m), c_() {}
    int deg_;
    vcp::matrix<T, P> c_;                                   // ndof(m) x 1
    // fe_space with any sparse policy SP creates/uses the same fe_function
    template <int DD, typename TT, typename PP, class SS> friend class fe_space;
};

template <int D, typename T, typename P = vcp::mats<T> >
class dual_vector {
public:
    int test_degree() const { return deg_; }                // degree of psi
    const vcp::matrix<T, P>& values() const { return v_; }
    // NOTE: no conversion to fe_function exists (X4, type-enforced)

private:
    explicit dual_vector(int m) : deg_(m), v_() {}
    int deg_;
    vcp::matrix<T, P> v_;
    template <int DD, typename TT, typename PP, class SS> friend class fe_space;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_FE_FUNCTION_HPP
