// vcp/bfem/sv/sv_constraint.hpp
// Phase 5d (Scott-Vogelius parts): the integer sparse constraint row
// (external design section 4.4) shared by sv_rows (generation) and
// linear_reduction (consumption).
//
// Recorded supplement to the v0.2 sketch: the design sketch writes
// "template <int D, typename T> struct sv_constraint", but neither parameter
// is used by the struct (dofs and coefficients are plain integers) and the
// consuming linear_reduction<T,P,SP> signature elides them; the struct is
// therefore defined parameter-free (noted in the deviation list of the gate
// report).

#ifndef VCP_BFEM_SV_SV_CONSTRAINT_HPP
#define VCP_BFEM_SV_SV_CONSTRAINT_HPP

#include <vector>

namespace vcp {
namespace bfem {

// one integer constraint row: sum_k coef[k] * x_{dof[k]} = 0
struct sv_constraint {
    std::vector<int> dof;
    std::vector<int> coef;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_SV_CONSTRAINT_HPP
