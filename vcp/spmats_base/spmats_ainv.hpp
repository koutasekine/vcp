// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_ainv.hpp
// Policy-layer types for the AINV approximate inverse (AINV-1): factored
// biconjugation R ~ A^{-1} = Z D^{-1} W^T after [BT98] (M. Benzi, M. Tuma,
// "A Sparse Approximate Inverse Preconditioner for Nonsymmetric Linear
// Systems", SIAM J. Sci. Comput. 19(3), 1998).  Types: ainv_status
// (+to_string) / ainv_options<T> / ainv_result<T,Index> (design
// ainv_design_v1.2.md §2.1).
// This file is included by spmats.hpp BEFORE the spmats class body (the
// types appear in policy method signatures).  The out-of-line policy method
// definitions live in spmats_base/spmats_ainv_impl.hpp (included after the
// class body, like spmats_ldl_impl.hpp).

#ifndef VCP_SPMATS_AINV_HPP
#define VCP_SPMATS_AINV_HPP

#include <cstddef>

#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {

// ainv_status: non-throwing status of the AINV construction / apply /
// residual-norm-estimate entry points (design §1.2: numerical events --
// tiny pivots, fill growth -- are NEVER a failure status; pivot lifting is
// applied silently and reported through the result diagnostics, D-4).
enum class ainv_status {
    success,        // outputs are valid (pivot modifications included)
    invalid_input,  // non-square input / dimension mismatch / non-diagonal D
    internal_error
};

inline const char* ainv_status_to_string(ainv_status s) {
    switch (s) {
    case ainv_status::success:        return "success";
    case ainv_status::invalid_input:  return "invalid_input";
    case ainv_status::internal_error: return "internal_error";
    }
    return "unknown";
}

// ainv_options<T>: dual-threshold dropping ([BT98] §8, ILUT-style) +
// pivot-lift control (design §3.3 / §4.3).  All tolerances are REAL-typed
// (they gate magnitudes) and T-dependent constants come from
// decimal_power_negative (no double literals).
template <class T>
struct ainv_options {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    real_type   drop_tolerance;          // tau: new fill-in with |value| < tau is dropped
    std::size_t max_nnz_per_column;      // 0 = unlimited; unit diagonal not counted
    real_type   pivot_small_threshold;   // |p_i| <= this -> lift
    real_type   pivot_lift_value;        // lifted pivot magnitude

    ainv_options()
        : drop_tolerance(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(2)),   // provisional (design v1.2 OPEN-2)
          max_nnz_per_column(0),                                            // provisional (design v1.2 OPEN-2)
          pivot_small_threshold(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(14)),  // provisional (design v1.2 OPEN-3)
          pivot_lift_value(
              vcp::tsparse_scalar::decimal_power_negative<real_type>(3)) {} // provisional (design v1.2 OPEN-3)
};

// ainv_result<T,Index>: status + integer diagnostics only (Z / W / D are
// out parameters of policy_ainv_with_info).  All diagnostics are advisory
// (D-4); fill is counted on the FACTORS, the same accounting as [BT98] §8
// (R itself is never materialized, D-13).
template <class T, class Index>
struct ainv_result {
    ainv_status status;

    Index n_pivot_modifications;   // lift count, Z-side p + W-side q combined
    Index first_modified_pivot;    // column index of the first lift; -1 = none
    Index nnz_Z, nnz_W;

    ainv_result()
        : status(ainv_status::internal_error),
          n_pivot_modifications(Index(0)),
          first_modified_pivot(Index(-1)),
          nnz_Z(Index(0)), nnz_W(Index(0)) {}
};

} // namespace vcp

#endif // VCP_SPMATS_AINV_HPP
