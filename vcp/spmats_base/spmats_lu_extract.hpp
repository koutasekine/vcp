// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_lu_extract.hpp
// Policy-layer types for the LU factor extraction (LUX-1) and the factor
// consuming layer (LUX-2): lu_extract_options / lu_extract_result /
// lu_apply_result (design lux_design_v0 SS2/SS2a; frozen decisions P-5,
// P-6, P-10).  This file is included by spmats.hpp BEFORE the spmats class
// body (the types appear in policy method signatures).  The out-of-line
// policy method definitions live in
// spmats_base/spmats_lu_extract_impl.hpp (included after the class body,
// like spmats_ldl_impl.hpp).

#ifndef VCP_SPMATS_LU_EXTRACT_HPP
#define VCP_SPMATS_LU_EXTRACT_HPP

#include <vcp/tsparse/tsparse_sparse_lu.hpp>
#include <vcp/tsparse/tsparse_sparse_lu_extract.hpp>

namespace vcp {

// P-5: lu_extract_options<T> wraps the SLU options in a single member.
// equilibration == true is rejected at the policy entry with
// unsupported_options (the SSC convention P.A.Q = L.U has no scaling; the
// SLU default is equilibration == false, so default usage never trips this).
template <class T>
struct lu_extract_options {
    sparse_lu_options<T> slu;

    lu_extract_options() : slu() {}
};

// P-6: lu_extract_result<T,Index>: status + integer diagnostics and the
// method/ordering usage record only (L / U / p / q are out parameters of
// policy_lu_with_info).  Field validity follows status; nnz_L / nnz_U are
// meaningful only on success.
//   method_used:   the method actually used, from sparse_lu_info
//                  (auto_select resolution recorded); when the factorization
//                  was not run (entry rejection) it echoes the request.
//   ordering_used: the ordering mode passed to SLU.  SLU does not re-export
//                  its auto_select resolution, so this records the REQUEST,
//                  not the resolved ordering.
template <class T, class Index>
struct lu_extract_result {
    sparse_lu_extract_status status;

    Index nnz_L, nnz_U;
    sparse_lu_method   method_used;
    sparse_lu_ordering ordering_used;

    lu_extract_result()
        : status(sparse_lu_extract_status::internal_error),
          nnz_L(Index(0)), nnz_U(Index(0)),
          method_used(sparse_lu_method::auto_select),
          ordering_used(sparse_lu_ordering::auto_select) {}
};

// ---------------------------------------------------------------------------
// factor consuming layer (LUX-2, design lux_design_v0 SS2a)
// ---------------------------------------------------------------------------

// P-10: lu_apply_result -- status-only result type of the factor consumers
// (policy_lu_solve_with_info / policy_lu_inverse_row_with_info).  Division
// by the U diagonal follows the certified three-branch gate (GT1 P1, same
// standard form as LDL SS3.1 with tol = 0): certified nonzero -> divide /
// certified zero (or structurally missing diagonal) -> singular_factor /
// undecidable -> inconclusive_division.  For totally ordered scalars the
// inconclusive branch is unreachable.
enum class lu_apply_status {
    success,
    dimension_mismatch,      // factor / permutation / rhs / index sizes inconsistent
    singular_factor,         // certified zero (or missing) U diagonal
    inconclusive_division,   // U diagonal nonzero-ness could not be certified (P1)
    invalid_input,           // structural contract violated (e.g. L not unit lower,
                             // p/q not permutations, U not upper triangular)
    internal_error
};

inline const char* lu_apply_status_to_string(lu_apply_status s) {
    switch (s) {
    case lu_apply_status::success:                return "success";
    case lu_apply_status::dimension_mismatch:     return "dimension_mismatch";
    case lu_apply_status::singular_factor:        return "singular_factor";
    case lu_apply_status::inconclusive_division:  return "inconclusive_division";
    case lu_apply_status::invalid_input:          return "invalid_input";
    case lu_apply_status::internal_error:         return "internal_error";
    }
    return "unknown";
}

struct lu_apply_result {
    lu_apply_status status;

    lu_apply_result() : status(lu_apply_status::internal_error) {}
};

} // namespace vcp

#endif // VCP_SPMATS_LU_EXTRACT_HPP
