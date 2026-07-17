// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License
//
// spmats_chol.hpp
// Policy-layer types for the sparse LL^T (Cholesky) factorization (CHOL-2):
// chol_options / chol_result (chol design v1 SS4.3).  This file is included
// by spmats.hpp BEFORE the spmats class body (the types appear in policy
// method signatures).  The out-of-line policy method definitions live in
// spmats_base/spmats_chol_impl.hpp (included after the class body, like
// spmats_ldl_impl.hpp).

#ifndef VCP_SPMATS_CHOL_HPP
#define VCP_SPMATS_CHOL_HPP

#include <vector>

#include <vcp/tsparse/tsparse_sparse_chol.hpp>

namespace vcp {

// the policy-layer options type is the tsparse options type itself (same
// convention as ldl_options).
template <class T>
using chol_options = sparse_chol_options<T>;

// chol_result<T,Index>: status + integer diagnostics only (L / perm are out
// parameters of policy_chol_with_info, design v1 SS4.3).  Field validity
// follows status; numeric sentinels are never used (P4).
template <class T, class Index>
struct chol_result {
    sparse_chol_status status;

    Index failure_at;            // -1 = none.  MATLAB flag = failure_at + 1
    Index inconclusive_at;       // -1 = none
    Index structural_empty_at;   // -1 = none (D-5)
    Index nnz_L;                 // symbolic count (kernel stored count; the
                                 // materialized L stores <= nnz_L after the
                                 // certified-zero drop, design SS4.3)
    sparse_chol_ordering ordering_used;  // auto_select resolution, always recorded
    sparse_chol_method   method_used;

    chol_result()
        : status(sparse_chol_status::internal_error),
          failure_at(Index(-1)), inconclusive_at(Index(-1)),
          structural_empty_at(Index(-1)), nnz_L(Index(0)),
          ordering_used(sparse_chol_ordering::auto_select),
          method_used(sparse_chol_method::auto_select) {}
};

} // namespace vcp

#endif // VCP_SPMATS_CHOL_HPP
