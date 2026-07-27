// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_ldl.hpp
// Policy-layer types for the sparse LDL^T factorization (LDL-3) and the
// inertia API (LDL-4): ldl_options / ldl_result / inertia_options /
// inertia_result / inertia_status (design v2 SS5.1, SS7.1, D-4).
// This file is included by spmats.hpp BEFORE the spmats class body (the
// types appear in policy method signatures).  The out-of-line policy method
// definitions live in spmats_base/spmats_ldl_impl.hpp (included after the
// class body, like spmats_lss.hpp).

#ifndef VCP_SPMATS_LDL_HPP
#define VCP_SPMATS_LDL_HPP

#include <vector>

#include <vcp/tsparse/tsparse_sparse_ldl.hpp>

namespace vcp {

// D-4: the policy-layer options type is the tsparse options type itself.
template <class T>
using ldl_options = sparse_ldl_options<T>;

// ldl_result<T,Index>: status + integer diagnostics only (L / D / perm are
// out parameters of policy_ldl_with_info, design v2 SS5.1).  Field validity
// follows status; numeric sentinels are never used (P4).
template <class T, class Index>
struct ldl_result {
    sparse_ldl_status status;

    Index n_pivots_1x1, n_pivots_2x2;   // executed pivots (zero-skips excluded)
    Index first_zero_pivot;             // -1 = none
    Index inconclusive_at;              // interrupted column; -1 = none
    Index structural_empty_at;          // first structurally empty row/col; -1 = none
    Index nnz_L;
    sparse_ldl_ordering ordering_used;  // auto_select resolution, always recorded
    sparse_ldl_method   method_used;
    bool dense_delegated;               // baseline_dynamic delegated to the dense kernel

    // ---- SLDL-SP diagnostics (design v1 SS4), conducted verbatim from the
    // tsparse result.  Integers and enums only (P4); growth_log2 is valid iff
    // growth_valid, and gemm_time_ns is meaningful iff gemm_call_count > 0
    // (it is never measured for non-floating scalars -- R-6).
    sparse_ldl_pivoting    pivot_mode_used;
    sparse_ldl_kernel_used diag_kernel_used;
    Index n_supernodes;
    Index max_supernode_width;
    Index n_boundary_splits;
    Index nnz_L_static;
    Index n_zero_skips;
    Index out_of_panel_at;              // -1 = none
    long long gemm_call_count;
    long long gemm_time_ns;
    int  growth_log2;
    bool growth_valid;

    ldl_result()
        : status(sparse_ldl_status::internal_error),
          n_pivots_1x1(Index(0)), n_pivots_2x2(Index(0)),
          first_zero_pivot(Index(-1)), inconclusive_at(Index(-1)),
          structural_empty_at(Index(-1)), nnz_L(Index(0)),
          ordering_used(sparse_ldl_ordering::auto_select),
          method_used(sparse_ldl_method::auto_select),
          dense_delegated(false),
          pivot_mode_used(sparse_ldl_pivoting::bk),
          diag_kernel_used(sparse_ldl_kernel_used::not_applicable),
          n_supernodes(Index(0)), max_supernode_width(Index(0)),
          n_boundary_splits(Index(0)), nnz_L_static(Index(0)),
          n_zero_skips(Index(0)), out_of_panel_at(Index(-1)),
          gemm_call_count(0), gemm_time_ns(0),
          growth_log2(0), growth_valid(false) {}
};

// ---------------------------------------------------------------------------
// inertia (LDL-4, design v2 SS7)
// ---------------------------------------------------------------------------

enum class inertia_status {
    success,
    not_block_diagonal,        // D-consumer: nonzero beyond bandwidth 1
    inconclusive_sign,         // sign/zero decision could not be certified (P1)
    factorization_failed,      // default _impl: internal ldl not success/zero_pivot
    invalid_input,
    internal_error
};

inline const char* inertia_status_to_string(inertia_status s) {
    switch (s) {
    case inertia_status::success:              return "success";
    case inertia_status::not_block_diagonal:   return "not_block_diagonal";
    case inertia_status::inconclusive_sign:    return "inconclusive_sign";
    case inertia_status::factorization_failed: return "factorization_failed";
    case inertia_status::invalid_input:        return "invalid_input";
    case inertia_status::internal_error:       return "internal_error";
    }
    return "unknown";
}

template <class Index>
struct inertia_result {
    Index n_pos, n_neg, n_zero;
    inertia_status status;
    Index inconclusive_at;              // position of an uncertifiable sign; -1 = none
    // diagnostics of the default _impl (design v2 SS7.1): the internal LDL
    // status.  Meaningful when the default implementation was used; kept at
    // success by the D-consumer entry points.
    sparse_ldl_status ldl_status;

    inertia_result()
        : n_pos(Index(0)), n_neg(Index(0)), n_zero(Index(0)),
          status(inertia_status::internal_error),
          inconclusive_at(Index(-1)),
          ldl_status(sparse_ldl_status::success) {}
};

// D-4: inertia_options = { zero_tol, ldl } (zero_tol default 0,
// certified-only; independent of ldl.zero_pivot_tol).  zero_tol is
// REAL-typed like every other tolerance in this module (it gates |d| --
// a real_type<T> magnitude; identical to T for real scalars, and required
// for the complex instantiation of the virtual policy _impl).
template <class T>
struct inertia_options {
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
    real_type zero_tol;
    ldl_options<T> ldl;

    inertia_options() : zero_tol(real_type(0)), ldl() {}
};

} // namespace vcp

#endif // VCP_SPMATS_LDL_HPP
