// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// LDL -- sparse symmetric LDL^T factorization (approximate layer).
//
// Convention (design v2 SS1.1, SS5.4, MATLAB sparse ldl orientation):
//     P^T * A * P = L * D * L^T,      A(p,p) = L * D * L^T (vector form),
// with perm p new->old ( (P^T A P)(i,j) = A(p[i], p[j]) ) and P(p[k],k) = 1.
// L is unit lower triangular (unit diagonal stored explicitly); D is a
// 1x1 / 2x2 block diagonal produced by certified Bunch-Kaufman symmetric
// pivoting (design v2 SS1.3).
//
// Input is assumed symmetric and only the LOWER triangle (diagonal included)
// is read by the numeric factorization (design v2 SS1.2); the upper triangle
// is referenced only by the optional symmetry check.
//
// Scalar contract: identical to the module scalar contract of
// tsparse_scalar.hpp (SLU-GT1 D8) -- arithmetic, certainly comparisons and
// ADL-resolved abs/sqrt only; no type-dependent branching.  For scalars
// whose comparisons cannot certify a pivot decision (e.g. kv::interval<TT>
// with zero-straddling entries) the factorization stops with the controlled
// status inconclusive_pivot_test (P1 third branch); numerical success is not
// guaranteed and not claimed.

#pragma once

#ifndef VCP_TSPARSE_SPARSE_LDL_HPP
#define VCP_TSPARSE_SPARSE_LDL_HPP

#include <cstddef>
#include <exception>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

// pivot_decision {acceptable, reject, inconclusive} is reused from the SLU
// header (design v2 SS1.3: include-only reuse; SLU itself is not modified).
// The pattern-only ordering functions (sparse_lu_rcm_ordering /
// sparse_lu_amd_ordering / sparse_lu_nested_dissection_ordering, LDL-1) also
// live in the impl injected by this header.
#include <vcp/tsparse/tsparse_sparse_lu.hpp>

namespace vcp {

// ===========================================================================
// enums / options / status / result (design v2 SS4.1-4.2)
// ===========================================================================

enum class sparse_ldl_method {
    auto_select,
    baseline_dynamic
    // Phase B (future, separate design doc): multifrontal
};

enum class sparse_ldl_ordering {
    auto_select,
    natural,
    rcm,
    amd,
    nested_dissection
    // colamd is intentionally absent at the type level (A^T A graph is for
    // nonsymmetric LU; design v2 SS3).
};

enum class sparse_ldl_status {
    success,
    structural_singularity,   // structurally empty row/column (assembly-bug hint)
    zero_pivot,               // certified zero pivot(s); factorization completed
    inconclusive_pivot_test,  // pivot decision could not be certified; stopped (P1)
    not_symmetric,            // symmetry of the input could not be certified
    invalid_options,          // colamd/pivot_threshold!=0/unimplemented method etc.
    invalid_input,            // n<0, malformed CSC, ...
    internal_error            // implementation bug / unexpected state (P3 net)
};

inline const char* sparse_ldl_status_to_string(sparse_ldl_status s) {
    switch (s) {
    case sparse_ldl_status::success:                 return "success";
    case sparse_ldl_status::structural_singularity: return "structural_singularity";
    case sparse_ldl_status::zero_pivot:              return "zero_pivot";
    case sparse_ldl_status::inconclusive_pivot_test: return "inconclusive_pivot_test";
    case sparse_ldl_status::not_symmetric:           return "not_symmetric";
    case sparse_ldl_status::invalid_options:         return "invalid_options";
    case sparse_ldl_status::invalid_input:           return "invalid_input";
    case sparse_ldl_status::internal_error:          return "internal_error";
    }
    return "unknown";
}

template <class T>
struct sparse_ldl_options {
    // Tolerances are REAL-typed (same convention as sparse_lu_options):
    // every certified gate compares magnitudes, i.e. values of
    // real_type<T>::type = the ADL abs return type.  For real scalars
    // real_type == T.
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    sparse_ldl_method   method;    // auto_select -> baseline_dynamic
    // ordering contract (design v2 SS0.1 decision 1): auto_select is resolved
    // by the library; the resolution is always reported in the diagnostic
    // ordering_used.  The resolution may change in future versions -- specify
    // an explicit ordering when reproducibility is required.
    sparse_ldl_ordering ordering;
    bool check_symmetry;
    real_type symmetry_tol;        // certified |a_ij - a_ji| <= tol accepts
    real_type zero_pivot_tol;      // 0 = only a certified exact zero
    real_type pivot_threshold;     // RESERVED (design v2 SS1.3): v0 rejects
                                   // any value not certified equal to 0 with
                                   // invalid_options.

    sparse_ldl_options()
        : method(sparse_ldl_method::auto_select),
          ordering(sparse_ldl_ordering::auto_select),
          check_symmetry(true),
          // same default policy as policy_is_symmetric (B-4: via the D4
          // customization point, not numeric_limits).
          symmetry_tol(vcp::tsparse_scalar::decimal_power_negative<real_type>(12u)),
          zero_pivot_tol(real_type(0)),
          pivot_threshold(real_type(0)) {}
};

template <class T, class Index>
struct sparse_ldl_result {
    sparse_ldl_status status;

    // L: CSC (n x n, unit diagonal stored explicitly), column rows ascending.
    // Valid only when status is success or zero_pivot.
    std::vector<Index> L_col_ptr, L_row_ind;
    std::vector<T>     L_val;

    // D: internal representation (the boundary layer converts to spmatrix;
    // value 0 entries are kept here, design v2 SS4.2).
    std::vector<T>     D_diag;     // d_kk
    std::vector<T>     D_sub;      // (k+1,k) element of a 2x2 block
    std::vector<char>  D_block2;   // D_block2[k]=1 <=> k leads a 2x2 block

    std::vector<Index> perm;       // new->old (design v2 SS5.4)

    // Diagnostics: integers only; validity is judged from status and the
    // Index sentinels below, never from numeric sentinels (P4).
    Index n_pivots_1x1, n_pivots_2x2;   // executed pivots (zero-skips excluded)
    Index first_zero_pivot;             // -1 = none
    Index inconclusive_at;              // interrupted column; -1 = none
    Index structural_empty_at;          // first structurally empty row/col; -1 = none
    Index nnz_L;
    sparse_ldl_ordering ordering_used;  // auto_select resolution, always recorded
    sparse_ldl_method   method_used;
    // LDL-2: true when baseline_dynamic delegated the numeric work to the
    // dense fallback kernel (small n).  method_used stays baseline_dynamic.
    bool dense_delegated;

    sparse_ldl_result()
        : status(sparse_ldl_status::internal_error),
          n_pivots_1x1(Index(0)), n_pivots_2x2(Index(0)),
          first_zero_pivot(Index(-1)), inconclusive_at(Index(-1)),
          structural_empty_at(Index(-1)), nnz_L(Index(0)),
          ordering_used(sparse_ldl_ordering::auto_select),
          method_used(sparse_ldl_method::auto_select),
          dense_delegated(false) {}
};

// ---------------------------------------------------------------------------
// dense fallback certified BK kernel (LDL-0; frozen regression reference)
// ---------------------------------------------------------------------------
#include <vcp/tsparse/detail/tsparse_sparse_ldl_dense_impl.hpp>

// ---------------------------------------------------------------------------
// sparse dynamic left-looking certified BK kernel (LDL-2)
// ---------------------------------------------------------------------------
#include <vcp/tsparse/detail/tsparse_sparse_ldl_sparse_impl.hpp>

// ===========================================================================
// entry: sparse_ldl_factorize_with_info (design v2 SS4.3)
// ===========================================================================

namespace sparse_ldl_detail {

// Input validation (LDL-0 Phase 2 step 1).  Returns true when the CSC triple
// is well-formed; malformed input maps to invalid_input at the entry.
template <class T, class Index>
bool sparse_ldl_validate_csc_(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val)
{
    if (n < Index(0)) return false;
    const std::size_t un = static_cast<std::size_t>(n);
    if (col_ptr.size() != un + 1u) return false;
    if (col_ptr[0] != Index(0)) return false;
    for (std::size_t c = 0; c < un; ++c) {
        if (col_ptr[c + 1u] < col_ptr[c]) return false;   // non-monotone
    }
    const std::size_t nnz = static_cast<std::size_t>(col_ptr[un]);
    if (row_ind.size() != nnz || val.size() != nnz) return false;
    for (std::size_t c = 0; c < un; ++c) {
        Index prev = Index(-1);
        for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r < Index(0) || r >= n) return false;     // out of range
            if (r <= prev) return false;                  // non-ascending / dup
            prev = r;
        }
    }
    return true;
}

// Symmetry check (LDL-0 Phase 2 step 3; D-2: transpose-pattern construction
// + per-column two-pointer merge, O(n + nnz)).  The whole input (upper
// triangle included) is compared against its transpose.  Acceptance is the
// certified side (P1): the input is declared symmetric only when
// |a_ij - a_ji| <= symmetry_tol is certified for every pair (one-sided
// entries compare against 0); both a certified mismatch and an uncertifiable
// comparison fall to not_symmetric.
template <class T, class Index>
bool sparse_ldl_symmetry_certified_(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const typename vcp::tsparse_scalar::real_type<T>::type& tol)
{
    using std::abs;
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;
    const std::size_t un = static_cast<std::size_t>(n);
    const std::size_t nnz = row_ind.size();

    // Transpose (CSC of A^T): counting sort, O(n + nnz); row lists come out
    // ascending because columns are scanned in ascending order.
    std::vector<Index> t_ptr(un + 1u, Index(0));
    for (std::size_t k = 0; k < nnz; ++k) {
        ++t_ptr[static_cast<std::size_t>(row_ind[k]) + 1u];
    }
    for (std::size_t i = 0; i < un; ++i) t_ptr[i + 1u] = t_ptr[i + 1u] + t_ptr[i];
    std::vector<Index> t_ind(nnz);
    std::vector<T>     t_val(nnz, T(0));
    {
        std::vector<Index> head(t_ptr.begin(), t_ptr.end() - 1);
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const std::size_t r = static_cast<std::size_t>(row_ind[static_cast<std::size_t>(k)]);
                const std::size_t pos = static_cast<std::size_t>(head[r]++);
                t_ind[pos] = static_cast<Index>(c);
                t_val[pos] = val[static_cast<std::size_t>(k)];
            }
        }
    }

    // Per-column two-pointer merge of A(:,j) against A^T(:,j).
    for (std::size_t j = 0; j < un; ++j) {
        Index pa = col_ptr[j];
        Index pt = t_ptr[j];
        const Index ea = col_ptr[j + 1u];
        const Index et = t_ptr[j + 1u];
        while (pa < ea || pt < et) {
            const Index ra = (pa < ea) ? row_ind[static_cast<std::size_t>(pa)] : n;
            const Index rt = (pt < et) ? t_ind[static_cast<std::size_t>(pt)] : n;
            T diff;
            if (ra == rt)      { diff = val[static_cast<std::size_t>(pa)] - t_val[static_cast<std::size_t>(pt)]; ++pa; ++pt; }
            else if (ra < rt)  { diff = val[static_cast<std::size_t>(pa)]; ++pa; }
            else               { diff = t_val[static_cast<std::size_t>(pt)]; ++pt; }
            const R adiff = abs(diff);
            if (adiff <= tol) { /* certified match: accept this pair */ }
            else { return false; }   // certified mismatch OR uncertifiable -> not symmetric
        }
    }
    return true;
}

} // namespace sparse_ldl_detail

// ---------------------------------------------------------------------------
// sparse_ldl_factorize_with_info -- non-throwing entry (P3: runtime failure
// is returned as a status; the final catch net below maps an escaped
// exception to internal_error, and its firing is treated as a defect).
//
// Processing order (fixed, LDL-0 Phase 2): input check -> options check ->
// symmetry check -> structural-singularity check -> numeric factorization.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_ldl_result<T, Index>
sparse_ldl_factorize_with_info(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const sparse_ldl_options<T>& opt)
{
    static_assert(std::is_signed<Index>::value, "sparse LDL Index must be signed");

    sparse_ldl_result<T, Index> res;
    res.ordering_used = opt.ordering;
    res.method_used   = opt.method;

    try {
        // ---- 1. input check -> invalid_input
        if (!sparse_ldl_detail::sparse_ldl_validate_csc_(n, col_ptr, row_ind, val)) {
            res.status = sparse_ldl_status::invalid_input;
            return res;
        }
        const std::size_t un = static_cast<std::size_t>(n);

        // ---- 2. options check -> invalid_options
        // pivot_threshold is reserved: accepted only when certified equal to 0.
        typedef typename vcp::tsparse_scalar::real_type<T>::type entry_real_type;
        if (opt.pivot_threshold == entry_real_type(0)) { /* certified zero: accepted */ }
        else { res.status = sparse_ldl_status::invalid_options; return res; }

        switch (opt.method) {
        case sparse_ldl_method::auto_select:
            res.method_used = sparse_ldl_method::baseline_dynamic;
            break;
        case sparse_ldl_method::baseline_dynamic:
            res.method_used = sparse_ldl_method::baseline_dynamic;
            break;
        default:
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }

        switch (opt.ordering) {
        case sparse_ldl_ordering::auto_select:
            // Contract (design v2 SS0.1 decision 1 / SS3): auto_select is
            // resolved BY THE LIBRARY; this version resolves it to amd.  The
            // resolution is always reported in ordering_used and may change
            // in future versions -- specify an explicit ordering when
            // reproducibility is required.  (The SLU auto_select semantics
            // are frozen and untouched by this.)
            res.ordering_used = sparse_ldl_ordering::amd;
            break;
        case sparse_ldl_ordering::natural:
        case sparse_ldl_ordering::rcm:
        case sparse_ldl_ordering::amd:
        case sparse_ldl_ordering::nested_dissection:
            res.ordering_used = opt.ordering;
            break;
        default:
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }

        // ---- 3. symmetry check -> not_symmetric (opt-out via check_symmetry)
        if (opt.check_symmetry) {
            if (!sparse_ldl_detail::sparse_ldl_symmetry_certified_(
                    n, col_ptr, row_ind, val, opt.symmetry_tol)) {
                res.status = sparse_ldl_status::not_symmetric;
                return res;
            }
        }

        // ---- 4. structural-singularity check -> structural_singularity
        // Effective symmetric pattern from the LOWER triangle only (i >= j):
        // entry (i,j) occupies row/col i and row/col j.  Strictly upper
        // stored entries are ignored by the factorization and therefore do
        // not count here either.  No factorization is attempted (priority
        // over zero_pivot, design v2 SS4.1).
        {
            std::vector<char> occupied(un, char(0));
            for (std::size_t c = 0; c < un; ++c) {
                for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                    const Index r = row_ind[static_cast<std::size_t>(k)];
                    if (r >= static_cast<Index>(c)) {
                        occupied[static_cast<std::size_t>(r)] = char(1);
                        occupied[c] = char(1);
                    }
                }
            }
            for (std::size_t i = 0; i < un; ++i) {
                if (occupied[i] == char(0)) {
                    // "最初の該当列を診断に" (LDL-0 Phase 2 step 4): reported
                    // in the dedicated integer diagnostic structural_empty_at
                    // (overloading first_zero_pivot / inconclusive_at would
                    // misreport their meaning).
                    res.status = sparse_ldl_status::structural_singularity;
                    res.structural_empty_at = static_cast<Index>(i);
                    return res;
                }
            }
        }

        // ---- 5. ordering (LDL-1): pattern-only pre-permutation P0
        // (perm0[new] = old) computed by the existing SLU ordering functions
        // (integer-only, A+A^T graph; the graph builder symmetrizes every
        // off-diagonal edge, so the lower-triangle CSC pattern is a valid
        // direct input).  SLU files are reused by include only, unmodified.
        std::vector<Index> perm0(un);
        for (std::size_t i = 0; i < un; ++i) perm0[i] = static_cast<Index>(i);
        switch (res.ordering_used) {
        case sparse_ldl_ordering::natural:
            break;   // identity
        case sparse_ldl_ordering::rcm:
            perm0 = sparse_lu_rcm_ordering(n, col_ptr, row_ind);
            break;
        case sparse_ldl_ordering::amd:
            perm0 = sparse_lu_amd_ordering(n, col_ptr, row_ind);
            break;
        case sparse_ldl_ordering::nested_dissection:
            perm0 = sparse_lu_nested_dissection_ordering(n, col_ptr, row_ind);
            break;
        default:
            // auto_select was already resolved above; reaching here is a bug.
            res.status = sparse_ldl_status::internal_error;
            return res;
        }

        // ---- 6. numeric factorization (LDL-2).  baseline_dynamic uses the
        // sparse dynamic left-looking kernel; small problems (n <= 64) are
        // delegated to the dense fallback kernel (reported via the
        // dense_delegated diagnostic; method_used stays baseline_dynamic).
        if (un <= 64u) {
            // CSC -> dense scatter of the lower triangle (strictly upper
            // entries are ignored, design v2 SS1.2) under P0: element (i,j)
            // lands at (pinv[i], pinv[j]) of P0^T A P0, on the lower side.
            std::vector<Index> pinv(un);   // old -> new
            for (std::size_t i = 0; i < un; ++i) pinv[static_cast<std::size_t>(perm0[i])] = static_cast<Index>(i);
            std::vector<T> work(un * un, T(0));
            for (std::size_t c = 0; c < un; ++c) {
                for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                    const Index r = row_ind[static_cast<std::size_t>(k)];
                    if (r >= static_cast<Index>(c)) {
                        const std::size_t a = static_cast<std::size_t>(pinv[static_cast<std::size_t>(r)]);
                        const std::size_t b = static_cast<std::size_t>(pinv[c]);
                        if (a >= b) work[a + b * un] = val[static_cast<std::size_t>(k)];
                        else        work[b + a * un] = val[static_cast<std::size_t>(k)];
                    }
                }
            }
            // BK exchanges are composed into perm in place, so seeding with
            // P0 yields the final perm = P0 o P_BK (new->old, design v2 SS5.4).
            res.perm = perm0;
            res.dense_delegated = true;
            sparse_ldl_dense_bk_factorize(n, work, res.perm, opt, res);
        } else {
            res.dense_delegated = false;
            sparse_ldl_sparse_bk_factorize(n, col_ptr, row_ind, val, perm0, opt, res);
        }
        return res;

    } catch (const std::exception&) {
        // P3 final protection net: reaching here means a certified gate was
        // missed somewhere upstream; treated as a defect, not a normal path.
        res.status = sparse_ldl_status::internal_error;
        return res;
    }
}

} // namespace vcp

#endif // VCP_TSPARSE_SPARSE_LDL_HPP
