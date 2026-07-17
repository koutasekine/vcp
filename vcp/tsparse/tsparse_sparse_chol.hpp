// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// CHOL -- sparse symmetric positive definite LL^T (Cholesky) factorization
// (approximate layer).
//
// Convention (chol design v1 SS0/SS1.1, MATLAB chol(A,'lower','vector')):
//     P^T * A * P = L * L^T,      A(p,p) = L * L^T (vector form),
// with perm p new->old ( (P^T A P)(i,j) = A(p[i], p[j]) ) and P(p[k],k) = 1.
// L is NON-UNIT lower triangular with positive diagonal (l_kk = sqrt of the
// certified-positive pivot); there is no D factor.  There is no pivoting, so
// perm is exactly the ordering output (static; SS5.6 difference table vs LDL).
//
// Numeric method (D-7/G2): simplicial up-looking in two separated stages --
// symbolic analysis (etree + ereach on the pattern of the EXPLICITLY
// CONSTRUCTED permuted lower CSC C = P0^T A P0) determines the structure of
// L first, then the numeric stage fills the saved pattern (ereach is NOT
// re-run on the numeric path).  There is no dense fallback (D-11).
//
// Input is assumed symmetric and only the LOWER triangle (diagonal included)
// is read by the factorization (design SS1.2); the upper triangle is
// referenced only by the optional symmetry check.
//
// Scalar contract: identical to the module scalar contract of
// tsparse_scalar.hpp (SLU-GT1 D8) -- arithmetic, certainly comparisons and
// ADL-resolved abs/sqrt only; no type-dependent branching.  The pivot gate is
// the certified three-branch (P1, design SS1.3): certified d > pd_tol ->
// sqrt and continue / certified d <= pd_tol -> not_positive_definite /
// neither certifiable -> inconclusive_pivot_test.  For scalars whose
// comparisons cannot certify (e.g. kv::interval<TT> with zero-straddling
// pivots) the controlled stop is the correct behaviour; numerical success is
// not guaranteed and not claimed.  NaN input falls to the third branch
// naturally (all comparisons false; no isnan is used, B-3 / design SS1.4).

#pragma once

#ifndef VCP_TSPARSE_SPARSE_CHOL_HPP
#define VCP_TSPARSE_SPARSE_CHOL_HPP

#include <cstddef>
#include <exception>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

// The pattern-only ordering functions (sparse_lu_rcm_ordering /
// sparse_lu_amd_ordering / sparse_lu_nested_dissection_ordering) are reused
// by include only (design SS5.5/SS7.1); SLU itself is not modified.  They are
// the ONLY SLU dependency of this kernel -- the SLU numeric / symbolic-reach
// / supernodal code paths are never entered.  ordering = natural is the
// SLU-independent escape route (design SS7.5-2).
#include <vcp/tsparse/tsparse_sparse_lu.hpp>

namespace vcp {

// ===========================================================================
// enums / options / status / result / symbolic (design SS4, SS8)
// ===========================================================================

enum class sparse_chol_method {
    auto_select,          // v1 resolves to simplicial_uplooking
    simplicial_uplooking
    // future (separate design doc): supernodal  (reserved by comment only,
    // design SS8-2; no enum value is added in v1)
    // [SPCM D-9 pure addition, 2026-07-17] the reservation above is realized
    // as a VALUE by the SPCM campaign (spcmodumar design v1.0 SS0.1 D-9):
    // supernodal is recorded in method_used by the external CHOLMOD
    // delegation policy (vcp::spcmodumar) ONLY, based on the backend's
    // L->is_super.  The native kernel NEVER sets it and NEVER accepts it as
    // an input method (an unmapped input method stays invalid_options).
    // chol design v1 SS4.1 / SS8-2 carry matching correction notes.
    , supernodal
};

enum class sparse_chol_ordering {
    auto_select,          // resolved to amd in this version (always reported
                          // in ordering_used; LDL decision-1 contract, G3)
    natural,              // identity permutation; calls no ordering function
    rcm,
    amd,
    nested_dissection
    // colamd is intentionally absent at the type level (A^T A graph is for
    // nonsymmetric LU; full parity with sparse_ldl_ordering, G3).
};

enum class sparse_chol_status {
    success,
    not_positive_definite,    // certified diagonal pivot <= pd_tol; also a
                              // structurally empty row/column (D-5, folded in
                              // and located by structural_empty_at) and a
                              // certified zero pivot (D-9, strict positive
                              // definiteness as in MATLAB chol)
    inconclusive_pivot_test,  // pivot decision could not be certified (P1)
    not_symmetric,            // symmetry of the input could not be certified
    invalid_options,          // unknown method / ordering enum values
    invalid_input,            // n<0, malformed CSC, ...
    internal_error            // implementation bug / unexpected state (P3
                              // net), and the ordering bijection-verification
                              // firewall (design SS7.5-1)
};

inline const char* sparse_chol_status_to_string(sparse_chol_status s) {
    switch (s) {
    case sparse_chol_status::success:                 return "success";
    case sparse_chol_status::not_positive_definite:   return "not_positive_definite";
    case sparse_chol_status::inconclusive_pivot_test: return "inconclusive_pivot_test";
    case sparse_chol_status::not_symmetric:           return "not_symmetric";
    case sparse_chol_status::invalid_options:         return "invalid_options";
    case sparse_chol_status::invalid_input:           return "invalid_input";
    case sparse_chol_status::internal_error:          return "internal_error";
    }
    return "unknown";
}

template <class T>
struct sparse_chol_options {
    // Tolerances are REAL-typed (same convention as sparse_ldl_options):
    // every certified gate compares values of real_type<T>::type = the ADL
    // abs return type.  For real scalars real_type == T.
    typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;

    sparse_chol_method   method;   // auto_select -> simplicial_uplooking
    // ordering contract (G3, LDL decision-1 parity): auto_select is resolved
    // by the library; the resolution is always reported in the diagnostic
    // ordering_used.  The resolution may change in future versions -- specify
    // an explicit ordering when reproducibility is required.
    sparse_chol_ordering ordering;
    bool check_symmetry;
    real_type symmetry_tol;        // certified |a_ij - a_ji| <= tol accepts
    real_type pd_tol;              // pivot acceptance: certified d > pd_tol.
                                   // default 0 = only certified d > 0 (D-9:
                                   // a certified zero pivot is rejected)

    sparse_chol_options()
        : method(sparse_chol_method::auto_select),
          ordering(sparse_chol_ordering::auto_select),
          check_symmetry(true),
          // same default policy as policy_is_symmetric (B-4: via the D4
          // customization point, not numeric_limits).
          symmetry_tol(vcp::tsparse_scalar::decimal_power_negative<real_type>(12u)),
          pd_tol(real_type(0)) {}
};

template <class T, class Index>
struct sparse_chol_result {
    sparse_chol_status status;

    // L: CSC (n x n, lower triangular, diagonal stored explicitly, rows
    // ascending within each column).  Valid ONLY when status is success
    // (D-3: no partial factor is returned).  The kernel keeps stored entries
    // whose value cancelled to exact zero (stored count == nnz_L, SS4.3);
    // the certified-zero drop happens at the spmats boundary layer only.
    std::vector<Index> L_col_ptr, L_row_ind;
    std::vector<T>     L_val;

    std::vector<Index> perm;   // new->old (= the ordering output; static,
                               // SS5.6).  Valid only when status is success.

    // Diagnostics: integers only (P4; no numeric sentinels).
    Index failure_at;            // -1 = none.  MATLAB flag = failure_at + 1
    Index inconclusive_at;       // -1 = none
    Index structural_empty_at;   // -1 = none (D-5)
    Index nnz_L;                 // symbolic count; valid as a symbolic value
                                 // even when the numeric stage stopped
    sparse_chol_ordering ordering_used;  // auto_select resolution, always recorded
    sparse_chol_method   method_used;

    sparse_chol_result()
        : status(sparse_chol_status::internal_error),
          failure_at(Index(-1)), inconclusive_at(Index(-1)),
          structural_empty_at(Index(-1)), nnz_L(Index(0)),
          ordering_used(sparse_chol_ordering::auto_select),
          method_used(sparse_chol_method::auto_select) {}
};

// ---------------------------------------------------------------------------
// sparse_chol_symbolic (D-7 / SS8): the symbolic analysis output, an
// independent struct so the numeric kernel can consume a symbolic object
// regardless of how it was produced (supernodal seed).  Required contract
// fields (design SS2.1): parent, L_col_ptr, L_row_ind, nnz_L.
// ---------------------------------------------------------------------------
template <class Index>
struct sparse_chol_symbolic {
    Index n;
    std::vector<Index> parent;      // elimination tree; -1 = root
    std::vector<Index> L_col_ptr;   // size n+1
    std::vector<Index> L_row_ind;   // diagonal first, rows ascending per column
    Index nnz_L;                    // == L_col_ptr[n]
    bool valid;                     // false = the input pattern violated the
                                    // lower-CSC contract (validity flag, not
                                    // a numeric sentinel; P4)

    sparse_chol_symbolic()
        : n(Index(0)), nnz_L(Index(0)), valid(false) {}
};

namespace sparse_chol_detail {

// ---------------------------------------------------------------------------
// Input validation (CHOL-0 step 2).  Same specification as
// sparse_ldl_detail::sparse_ldl_validate_csc_ (copied and renamed; the LDL
// detail is deliberately NOT called -- LDL is frozen, B-6).  Returns true
// when the CSC triple is well-formed; malformed input maps to invalid_input
// at the entry.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool sparse_chol_validate_csc_(
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

// ---------------------------------------------------------------------------
// Symmetry check (CHOL-0 step 3; same method as
// sparse_ldl_detail::sparse_ldl_symmetry_certified_, copied and renamed:
// transpose-pattern construction + per-column two-pointer merge, O(n + nnz)).
// The whole input (upper triangle included) is compared against its
// transpose.  Acceptance is the certified side (P1): the input is declared
// symmetric only when |a_ij - a_ji| <= symmetry_tol is certified for every
// pair (one-sided entries compare against 0); both a certified mismatch and
// an uncertifiable comparison fall to not_symmetric (design SS1.2).
// ---------------------------------------------------------------------------
template <class T, class Index>
bool sparse_chol_symmetry_certified_(
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

// ---------------------------------------------------------------------------
// Structural-empty scan (CHOL-0 step 4; same shape as the LDL entry scan,
// D-5: only the reporting differs -- the CALLER maps a hit to
// not_positive_definite with structural_empty_at, instead of LDL's
// structural_singularity).  Effective symmetric pattern from the LOWER
// triangle only (i >= j): entry (i,j) occupies row/col i and row/col j.
// Strictly upper stored entries are ignored by the factorization and
// therefore do not count here either.  Returns the first structurally empty
// row/col index, or -1 when none.
// ---------------------------------------------------------------------------
template <class Index>
Index sparse_chol_scan_structural_empty_(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    const std::size_t un = static_cast<std::size_t>(n);
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
        if (occupied[i] == char(0)) return static_cast<Index>(i);
    }
    return Index(-1);
}

// ---------------------------------------------------------------------------
// Ordering bijection-verification firewall (design SS7.5-1, D-13).
// O(n): size, range and a visited bitmap.  perm passes only when it is a
// bijection of {0..n-1}.  A failure at the entry maps to internal_error --
// the factorization's correctness is independent of the ordering QUALITY,
// so this check cuts the only path by which an upstream ordering bug could
// turn into silent corruption.
// ---------------------------------------------------------------------------
template <class Index>
bool sparse_chol_verify_permutation_(
    const Index n,
    const std::vector<Index>& perm)
{
    if (n < Index(0)) return false;
    const std::size_t un = static_cast<std::size_t>(n);
    if (perm.size() != un) return false;
    std::vector<char> seen(un, char(0));
    for (std::size_t i = 0; i < un; ++i) {
        const Index p = perm[i];
        if (p < Index(0) || p >= n) return false;      // out of range
        if (seen[static_cast<std::size_t>(p)] != char(0)) return false;  // repeat
        seen[static_cast<std::size_t>(p)] = char(1);
    }
    return true;
}

// ---------------------------------------------------------------------------
// Ordering dispatch (CHOL-0 step 5; same shape as the LDL entry switch).
// `ordering` must already be RESOLVED (auto_select decided by the caller and
// recorded in ordering_used).  natural is the identity permutation and calls
// no ordering function (SLU-independent escape route, SS7.5-2).  The SLU
// ordering functions are pattern-only / integer-only and symmetrize every
// off-diagonal edge {r,c} internally, so the lower-triangle CSC pattern is a
// valid direct input (design SS5.5).  Returns false for an unresolved /
// unknown enum (caller maps to internal_error -- auto_select must have been
// resolved before this point).
// ---------------------------------------------------------------------------
template <class Index>
bool sparse_chol_compute_ordering_(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const sparse_chol_ordering ordering,
    std::vector<Index>& perm0)
{
    const std::size_t un = static_cast<std::size_t>(n);
    perm0.resize(un);
    for (std::size_t i = 0; i < un; ++i) perm0[i] = static_cast<Index>(i);
    switch (ordering) {
    case sparse_chol_ordering::natural:
        return true;   // identity
    case sparse_chol_ordering::rcm:
        perm0 = sparse_lu_rcm_ordering(n, col_ptr, row_ind);
        return true;
    case sparse_chol_ordering::amd:
        perm0 = sparse_lu_amd_ordering(n, col_ptr, row_ind);
        return true;
    case sparse_chol_ordering::nested_dissection:
        perm0 = sparse_lu_nested_dissection_ordering(n, col_ptr, row_ind);
        return true;
    default:
        return false;  // auto_select not resolved / unknown: caller bug
    }
}

// ---------------------------------------------------------------------------
// Explicit construction of the permuted lower CSC C = P0^T A P0 (CHOL-0
// step 6; design SS2.0-4 -- an INTENTIONAL difference from LDL, which never
// materializes the permuted matrix: the chol symbolic analysis needs the
// permuted pattern, SS5.6).  Only stored LOWER entries (r >= c) of A are
// read (SS1.2); entry (r,c) lands at (max(a,b), min(a,b)) with a = pinv[r],
// b = pinv[c].  Distinct read entries map to distinct unordered pairs, so no
// duplicates can arise.  Rows come out ascending within each column by a
// two-pass stable counting sort (by row, then by column); the ascending
// invariant is re-checked defensively (false -> internal_error at the
// caller).
// ---------------------------------------------------------------------------
template <class T, class Index>
bool sparse_chol_build_permuted_lower_csc_(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const std::vector<Index>& perm0,
    std::vector<Index>& C_col_ptr,
    std::vector<Index>& C_row_ind,
    std::vector<T>&     C_val)
{
    const std::size_t un = static_cast<std::size_t>(n);

    std::vector<Index> pinv(un);   // old -> new
    for (std::size_t i = 0; i < un; ++i) pinv[static_cast<std::size_t>(perm0[i])] = static_cast<Index>(i);

    // gather the permuted lower entries (source order)
    std::vector<Index> er, ec;
    std::vector<T>     ev;
    for (std::size_t c = 0; c < un; ++c) {
        for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r >= static_cast<Index>(c)) {
                const Index a = pinv[static_cast<std::size_t>(r)];
                const Index b = pinv[c];
                er.push_back(a >= b ? a : b);
                ec.push_back(a >= b ? b : a);
                ev.push_back(val[static_cast<std::size_t>(k)]);
            }
        }
    }
    const std::size_t m = er.size();

    // stable counting sort by ROW, then stable counting sort by COLUMN:
    // the final order is by column with rows ascending inside each column.
    std::vector<std::size_t> ord(m), ord2(m);
    {
        std::vector<Index> cnt(un + 1u, Index(0));
        for (std::size_t t = 0; t < m; ++t) ++cnt[static_cast<std::size_t>(er[t]) + 1u];
        for (std::size_t i = 0; i < un; ++i) cnt[i + 1u] = cnt[i + 1u] + cnt[i];
        std::vector<Index> head(cnt.begin(), cnt.end() - 1);
        for (std::size_t t = 0; t < m; ++t) {
            ord[static_cast<std::size_t>(head[static_cast<std::size_t>(er[t])]++)] = t;
        }
    }
    C_col_ptr.assign(un + 1u, Index(0));
    {
        std::vector<Index> cnt(un + 1u, Index(0));
        for (std::size_t t = 0; t < m; ++t) ++cnt[static_cast<std::size_t>(ec[t]) + 1u];
        for (std::size_t i = 0; i < un; ++i) cnt[i + 1u] = cnt[i + 1u] + cnt[i];
        C_col_ptr = cnt;
        std::vector<Index> head(cnt.begin(), cnt.end() - 1);
        for (std::size_t t = 0; t < m; ++t) {
            const std::size_t e = ord[t];
            ord2[static_cast<std::size_t>(head[static_cast<std::size_t>(ec[e])]++)] = e;
        }
    }
    C_row_ind.resize(m);
    C_val.assign(m, T(0));
    for (std::size_t t = 0; t < m; ++t) {
        C_row_ind[t] = er[ord2[t]];
        C_val[t]     = ev[ord2[t]];
    }

    // defensive re-check of the construction invariant (strictly ascending
    // rows >= the column index inside every column; violation would mean a
    // duplicate or an implementation bug -> internal_error at the caller).
    for (std::size_t c = 0; c < un; ++c) {
        Index prev = static_cast<Index>(c) - Index(1);
        for (Index k = C_col_ptr[c]; k < C_col_ptr[c + 1u]; ++k) {
            const Index r = C_row_ind[static_cast<std::size_t>(k)];
            if (r <= prev || r >= n) return false;
            prev = r;
        }
    }
    return true;
}

} // namespace sparse_chol_detail

// ---------------------------------------------------------------------------
// sparse_chol_symbolic_analyze (CHOL-0 step 7; design SS2.1, D-7).
//
// PATTERN-ONLY (G-0.3: no value array appears in the signature).  Input:
// the pattern of a lower-triangular CSC with rows ascending inside every
// column (as produced by sparse_chol_build_permuted_lower_csc_); the
// diagonal entry may be structurally absent (a missing diagonal is a
// NUMERIC question answered by the pivot gate, not a symbolic one).
//
// Output: sparse_chol_symbolic with
//   parent    : elimination tree (parent[j] = -1 for a root), built with the
//               standard ancestor-compression traversal (each traversed node
//               has its ancestor pointer compressed to the current row k)
//               [textbook algorithm; correctness gated by G-0.2],
//   L_col_ptr / L_row_ind : the exact pattern of L (diagonal stored
//               explicitly and first in each column, rows ascending), from
//               one ereach (etree row subtree walk) per row,
//   nnz_L     : L_col_ptr[n].
// The saved pattern is REUSED by the numeric kernel (ereach is not re-run on
// the numeric path; frozen decision in the implementation plan).
//
// valid == false reports a pattern that violates the lower-CSC contract (or
// an etree walk that escaped its row -- impossible for a well-formed lower
// pattern, kept as a defensive net); the one-shot entry maps it to
// internal_error because it only passes self-constructed patterns here.
// ---------------------------------------------------------------------------
template <class Index>
sparse_chol_symbolic<Index>
sparse_chol_symbolic_analyze(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    static_assert(std::is_signed<Index>::value, "sparse CHOL Index must be signed");

    sparse_chol_symbolic<Index> sym;
    sym.valid = false;
    if (n < Index(0)) return sym;
    const std::size_t un = static_cast<std::size_t>(n);
    if (col_ptr.size() != un + 1u) return sym;
    if (un > 0 && col_ptr[0] != Index(0)) return sym;

    // light pattern validation (lower CSC, ascending rows, in range)
    for (std::size_t c = 0; c < un; ++c) {
        if (col_ptr[c + 1u] < col_ptr[c]) return sym;
        Index prev = static_cast<Index>(c) - Index(1);
        for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r <= prev || r >= n) return sym;
            prev = r;
        }
    }
    if (row_ind.size() != static_cast<std::size_t>(col_ptr[un])) return sym;

    sym.n = n;
    sym.parent.assign(un, Index(-1));
    sym.L_col_ptr.assign(un + 1u, Index(0));
    sym.L_row_ind.clear();
    sym.nnz_L = Index(0);

    // ---- row-wise view of the pattern (CSR of C): row k lists the columns
    // j <= k with C(k,j) stored, ascending (counting sort over ascending
    // column scan).  This is what both the etree and the ereach walks
    // consume (they need "the entries of row k left of the diagonal").
    std::vector<Index> r_ptr(un + 1u, Index(0));
    const std::size_t nnzC = row_ind.size();
    std::vector<Index> r_ind(nnzC);
    {
        std::vector<Index> cnt(un + 1u, Index(0));
        for (std::size_t t = 0; t < nnzC; ++t) ++cnt[static_cast<std::size_t>(row_ind[t]) + 1u];
        for (std::size_t i = 0; i < un; ++i) cnt[i + 1u] = cnt[i + 1u] + cnt[i];
        r_ptr = cnt;
        std::vector<Index> head(cnt.begin(), cnt.end() - 1);
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const std::size_t r = static_cast<std::size_t>(row_ind[static_cast<std::size_t>(k)]);
                r_ind[static_cast<std::size_t>(head[r]++)] = static_cast<Index>(c);
            }
        }
    }

    // ---- elimination tree with ancestor compression (values never read)
    {
        std::vector<Index> ancestor(un, Index(-1));
        for (std::size_t k = 0; k < un; ++k) {
            for (Index p = r_ptr[k]; p < r_ptr[k + 1u]; ++p) {
                Index i = r_ind[static_cast<std::size_t>(p)];
                while (i != Index(-1) && i < static_cast<Index>(k)) {
                    const Index inext = ancestor[static_cast<std::size_t>(i)];
                    ancestor[static_cast<std::size_t>(i)] = static_cast<Index>(k);
                    if (inext == Index(-1)) sym.parent[static_cast<std::size_t>(i)] = static_cast<Index>(k);
                    i = inext;
                }
            }
        }
    }

    // ---- ereach per row k: the strictly-lower pattern of L's row k is the
    // union of the etree paths from every entry j < k of row k up to (but
    // excluding) k.  Membership only is recorded (the numeric kernel
    // processes row patterns in ascending order, which is a valid
    // topological order because parent > child in the etree).
    std::vector<Index> rowpat_ptr(un + 1u, Index(0));
    std::vector<Index> rowpat;             // strict-lower row patterns, flat
    {
        std::vector<Index> wmark(un, Index(-1));   // stamp = k
        std::vector<Index> path;
        path.reserve(un);
        for (std::size_t k = 0; k < un; ++k) {
            wmark[k] = static_cast<Index>(k);      // mark k itself
            for (Index p = r_ptr[k]; p < r_ptr[k + 1u]; ++p) {
                Index j = r_ind[static_cast<std::size_t>(p)];
                if (j >= static_cast<Index>(k)) continue;   // diagonal
                while (j != Index(-1) && wmark[static_cast<std::size_t>(j)] != static_cast<Index>(k)) {
                    path.push_back(j);
                    wmark[static_cast<std::size_t>(j)] = static_cast<Index>(k);
                    j = sym.parent[static_cast<std::size_t>(j)];
                }
                if (j == Index(-1)) {
                    // the walk escaped the row-k subtree: impossible for a
                    // well-formed lower pattern (etree theorem); defensive.
                    sym.valid = false;
                    return sym;
                }
            }
            for (std::size_t t = 0; t < path.size(); ++t) rowpat.push_back(path[t]);
            path.clear();
            rowpat_ptr[k + 1u] = static_cast<Index>(rowpat.size());
        }
    }

    // ---- CSC structure of L: column j holds its diagonal first, then the
    // rows k > j that reference j in their row pattern, ascending (k is
    // appended in ascending scan order).
    {
        std::vector<Index> cnt(un, Index(1));      // 1 = explicit diagonal
        for (std::size_t t = 0; t < rowpat.size(); ++t) ++cnt[static_cast<std::size_t>(rowpat[t])];
        for (std::size_t j = 0; j < un; ++j) sym.L_col_ptr[j + 1u] = sym.L_col_ptr[j] + cnt[j];
        sym.nnz_L = sym.L_col_ptr[un];
        sym.L_row_ind.assign(static_cast<std::size_t>(sym.nnz_L), Index(0));
        std::vector<Index> head(un);
        for (std::size_t j = 0; j < un; ++j) {
            head[j] = sym.L_col_ptr[j];
            sym.L_row_ind[static_cast<std::size_t>(head[j]++)] = static_cast<Index>(j);   // diagonal first
        }
        for (std::size_t k = 0; k < un; ++k) {
            for (Index t = rowpat_ptr[k]; t < rowpat_ptr[k + 1u]; ++t) {
                const std::size_t j = static_cast<std::size_t>(rowpat[static_cast<std::size_t>(t)]);
                sym.L_row_ind[static_cast<std::size_t>(head[j]++)] = static_cast<Index>(k);
            }
        }
    }

    sym.valid = true;
    return sym;
}

// ---------------------------------------------------------------------------
// sparse_chol_numeric_factorize (CHOL-1 step 1; design SS2.2, D-7).
//
// Up-looking numeric stage over the SAVED symbolic pattern (ereach is NOT
// re-run here; the row-wise views below are pattern transposes of the saved
// structure, not re-analyses).  Consumes any sparse_chol_symbolic regardless
// of how it was produced (SS8-3).
//
//   sym         : valid symbolic analysis of the pattern of C
//   C_*         : permuted lower CSC C = P0^T A P0 (rows ascending; only
//                 producer-validated input is expected -- the one-shot entry
//                 passes its own construction)
//   opt         : pd_tol / (method already resolved by the caller)
//   res         : L_col_ptr / L_row_ind / L_val / status / failure_at /
//                 inconclusive_at / nnz_L are written here.  perm /
//                 ordering_used / method_used belong to the caller.
//
// Row k processing (textbook up-looking; correctness gated by G-1.1/G-1.2
// dense-reference comparison): the sparse triangular solve accumulates, for
// every j in L's row-k pattern in ascending order (a valid topological order
// of the etree paths since parent > child),
//     x_j   = C(k,j) - sum_{t<j, L(k,t)!=0} L(j,t) * L(k,t),
//     L(k,j) = x_j / l_jj,
//     d      = C(k,k) - sum_j L(k,j)^2,
// and the diagonal pivot d passes the certified three-branch (P1, SS1.3):
//   certified d >  pd_tol -> l_kk = sqrt(d) (ADL, unqualified; P2/B-2)
//   certified d <= pd_tol -> not_positive_definite, failure_at = k (D-9)
//   neither certified     -> inconclusive_pivot_test, inconclusive_at = k
// The sign gate goes through real_part (identity for real and interval
// scalars; keeps the complex instantiation of the virtual policy _impl
// compilable, same convention as the LDL/inertia layer -- LL^H itself is out
// of scope).  NaN falls to the third branch naturally (SS1.4, B-3).
//
// The kernel keeps every stored position of the symbolic pattern in L_val
// (numerically cancelled zeros INCLUDED: stored count == nnz_L, SS4.3); the
// certified-zero drop is the spmats boundary layer's job.  On a stop, L is
// cleared (D-3: no partial factor) and only the diagnostics remain.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_chol_numeric_factorize(
    const sparse_chol_symbolic<Index>& sym,
    const std::vector<Index>& C_col_ptr,
    const std::vector<Index>& C_row_ind,
    const std::vector<T>&     C_val,
    const sparse_chol_options<T>& opt,
    sparse_chol_result<T, Index>& res)
{
    static_assert(std::is_signed<Index>::value, "sparse CHOL Index must be signed");
    using std::sqrt;
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    res.failure_at      = Index(-1);
    res.inconclusive_at = Index(-1);
    res.L_col_ptr.clear();
    res.L_row_ind.clear();
    res.L_val.clear();

    // defensive: a producer-validated symbolic object and a size-consistent C
    if (!sym.valid || sym.n < Index(0)) {
        res.status = sparse_chol_status::internal_error;
        return;
    }
    const Index n = sym.n;
    const std::size_t un = static_cast<std::size_t>(n);
    if (C_col_ptr.size() != un + 1u ||
        C_row_ind.size() != static_cast<std::size_t>(C_col_ptr[un]) ||
        C_val.size()     != C_row_ind.size() ||
        sym.L_col_ptr.size() != un + 1u ||
        sym.L_row_ind.size() != static_cast<std::size_t>(sym.nnz_L)) {
        res.status = sparse_chol_status::internal_error;
        return;
    }

    res.nnz_L     = sym.nnz_L;
    res.L_col_ptr = sym.L_col_ptr;
    res.L_row_ind = sym.L_row_ind;
    res.L_val.assign(static_cast<std::size_t>(sym.nnz_L), T(0));

    // ---- row-wise view of C (values): counting sort; columns scanned
    // ascending, so each row list carries ascending column indices.
    const std::size_t nnzC = C_row_ind.size();
    std::vector<Index> cr_ptr(un + 1u, Index(0));
    std::vector<Index> cr_ind(nnzC);
    std::vector<T>     cr_val(nnzC, T(0));
    {
        std::vector<Index> cnt(un + 1u, Index(0));
        for (std::size_t t = 0; t < nnzC; ++t) ++cnt[static_cast<std::size_t>(C_row_ind[t]) + 1u];
        for (std::size_t i = 0; i < un; ++i) cnt[i + 1u] = cnt[i + 1u] + cnt[i];
        cr_ptr = cnt;
        std::vector<Index> head(cnt.begin(), cnt.end() - 1);
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = C_col_ptr[c]; k < C_col_ptr[c + 1u]; ++k) {
                const std::size_t r = static_cast<std::size_t>(C_row_ind[static_cast<std::size_t>(k)]);
                const std::size_t pos = static_cast<std::size_t>(head[r]++);
                cr_ind[pos] = static_cast<Index>(c);
                cr_val[pos] = C_val[static_cast<std::size_t>(k)];
            }
        }
    }

    // ---- row-wise pattern of the strictly-lower part of L (pattern
    // transpose of the saved structure; ascending column indices per row).
    std::vector<Index> lr_ptr(un + 1u, Index(0));
    std::vector<Index> lr_ind(static_cast<std::size_t>(sym.nnz_L) >= un
                              ? static_cast<std::size_t>(sym.nnz_L) - un : 0u);
    {
        std::vector<Index> cnt(un + 1u, Index(0));
        for (std::size_t j = 0; j < un; ++j) {
            for (Index p = sym.L_col_ptr[j] + Index(1); p < sym.L_col_ptr[j + 1u]; ++p) {
                ++cnt[static_cast<std::size_t>(sym.L_row_ind[static_cast<std::size_t>(p)]) + 1u];
            }
        }
        for (std::size_t i = 0; i < un; ++i) cnt[i + 1u] = cnt[i + 1u] + cnt[i];
        lr_ptr = cnt;
        std::vector<Index> head(cnt.begin(), cnt.end() - 1);
        for (std::size_t j = 0; j < un; ++j) {
            for (Index p = sym.L_col_ptr[j] + Index(1); p < sym.L_col_ptr[j + 1u]; ++p) {
                const std::size_t r = static_cast<std::size_t>(sym.L_row_ind[static_cast<std::size_t>(p)]);
                lr_ind[static_cast<std::size_t>(head[r]++)] = static_cast<Index>(j);
            }
        }
    }

    // ---- up-looking main loop
    std::vector<T> x(un, T(0));            // dense workspace; every touched
                                           // position is inside row k's saved
                                           // pattern and is zeroed on read
    std::vector<Index> fill(un, Index(0)); // next append position per column
    for (std::size_t j = 0; j < un; ++j) fill[j] = sym.L_col_ptr[j] + Index(1);

    bool corrupted = false;                // defensive symbolic/numeric mismatch
    for (std::size_t k = 0; k < un && !corrupted; ++k) {
        // scatter row k of C (j < k into x; the diagonal seeds d)
        T d = T(0);
        for (Index p = cr_ptr[k]; p < cr_ptr[k + 1u]; ++p) {
            const Index j = cr_ind[static_cast<std::size_t>(p)];
            if (j == static_cast<Index>(k)) d = cr_val[static_cast<std::size_t>(p)];
            else if (j < static_cast<Index>(k)) x[static_cast<std::size_t>(j)] = cr_val[static_cast<std::size_t>(p)];
            else { corrupted = true; break; }   // C not lower triangular
        }
        if (corrupted) break;

        // sparse triangular solve along the saved row pattern (ascending)
        for (Index t = lr_ptr[k]; t < lr_ptr[k + 1u] && !corrupted; ++t) {
            const std::size_t j = static_cast<std::size_t>(lr_ind[static_cast<std::size_t>(t)]);
            const T xj = x[j];
            x[j] = T(0);
            const T ljj = res.L_val[static_cast<std::size_t>(sym.L_col_ptr[j])];
            const T lkj = xj / ljj;   // l_jj passed the certified d > pd_tol
                                      // gate in row j: no zero division on
                                      // this path (SS1.4)
            for (Index p = sym.L_col_ptr[j] + Index(1); p < fill[j]; ++p) {
                x[static_cast<std::size_t>(sym.L_row_ind[static_cast<std::size_t>(p)])] -=
                    res.L_val[static_cast<std::size_t>(p)] * lkj;
            }
            if (fill[j] >= sym.L_col_ptr[j + 1u] ||
                sym.L_row_ind[static_cast<std::size_t>(fill[j])] != static_cast<Index>(k)) {
                corrupted = true;     // append slot missing / not row k's slot
                break;
            }
            res.L_val[static_cast<std::size_t>(fill[j])] = lkj;
            fill[j] = fill[j] + Index(1);
            d -= lkj * lkj;
        }
        if (corrupted) break;

        // certified three-branch on the diagonal pivot (P1, SS1.3; B-5: the
        // success side is a certified >, never a !(x > tol) idiom)
        const R rd = vcp::tsparse_scalar::real_part(d);
        if (rd > opt.pd_tol) {
            if (sym.L_row_ind[static_cast<std::size_t>(sym.L_col_ptr[k])] != static_cast<Index>(k)) {
                corrupted = true;
                break;
            }
            res.L_val[static_cast<std::size_t>(sym.L_col_ptr[k])] = sqrt(d);
        } else if (rd <= opt.pd_tol) {
            res.status = sparse_chol_status::not_positive_definite;
            res.failure_at = static_cast<Index>(k);
            res.L_col_ptr.clear(); res.L_row_ind.clear(); res.L_val.clear();
            return;
        } else {
            res.status = sparse_chol_status::inconclusive_pivot_test;
            res.inconclusive_at = static_cast<Index>(k);
            res.L_col_ptr.clear(); res.L_row_ind.clear(); res.L_val.clear();
            return;
        }
    }

    if (corrupted) {
        res.status = sparse_chol_status::internal_error;
        res.L_col_ptr.clear(); res.L_row_ind.clear(); res.L_val.clear();
        return;
    }
    res.status = sparse_chol_status::success;
}

// ---------------------------------------------------------------------------
// sparse_chol_factorize_with_info -- non-throwing one-shot entry (CHOL-1
// step 2; P3: runtime failure is returned as a status; the final catch net
// maps an escaped exception to internal_error, and its firing is treated as
// a defect).  This is the only pipeline the policy layer publishes (D-7).
//
// Processing order (fixed): input check -> options check -> symmetry check
// -> structural-empty check (D-5) -> ordering (+ bijection firewall,
// SS7.5-1) -> permuted lower CSC -> symbolic analysis -> numeric stage.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_chol_result<T, Index>
sparse_chol_factorize_with_info(
    Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const sparse_chol_options<T>& opt)
{
    static_assert(std::is_signed<Index>::value, "sparse CHOL Index must be signed");

    sparse_chol_result<T, Index> res;
    res.ordering_used = opt.ordering;
    res.method_used   = opt.method;

    try {
        // ---- 1. input check -> invalid_input
        if (!sparse_chol_detail::sparse_chol_validate_csc_(n, col_ptr, row_ind, val)) {
            res.status = sparse_chol_status::invalid_input;
            return res;
        }

        // ---- 2. options check -> invalid_options
        switch (opt.method) {
        case sparse_chol_method::auto_select:
            res.method_used = sparse_chol_method::simplicial_uplooking;
            break;
        case sparse_chol_method::simplicial_uplooking:
            res.method_used = sparse_chol_method::simplicial_uplooking;
            break;
        default:
            res.status = sparse_chol_status::invalid_options;
            return res;
        }

        switch (opt.ordering) {
        case sparse_chol_ordering::auto_select:
            // Contract (G3, LDL decision-1 parity): auto_select is resolved
            // BY THE LIBRARY; this version resolves it to amd.  The
            // resolution is always reported in ordering_used and may change
            // in future versions -- specify an explicit ordering when
            // reproducibility is required.  (The SLU auto_select semantics
            // are frozen and untouched by this.)
            res.ordering_used = sparse_chol_ordering::amd;
            break;
        case sparse_chol_ordering::natural:
        case sparse_chol_ordering::rcm:
        case sparse_chol_ordering::amd:
        case sparse_chol_ordering::nested_dissection:
            res.ordering_used = opt.ordering;
            break;
        default:
            res.status = sparse_chol_status::invalid_options;
            return res;
        }

        // ---- 3. symmetry check -> not_symmetric (opt-out via check_symmetry)
        if (opt.check_symmetry) {
            if (!sparse_chol_detail::sparse_chol_symmetry_certified_(
                    n, col_ptr, row_ind, val, opt.symmetry_tol)) {
                res.status = sparse_chol_status::not_symmetric;
                return res;
            }
        }

        // ---- 4. structural-empty check -> not_positive_definite (D-5: the
        // LDL structural_singularity is folded into not_positive_definite;
        // the location is reported in the dedicated integer diagnostic
        // structural_empty_at; failure_at stays -1 -- no pivot was tested).
        {
            const Index empty_at = sparse_chol_detail::sparse_chol_scan_structural_empty_(
                n, col_ptr, row_ind);
            if (empty_at >= Index(0)) {
                res.status = sparse_chol_status::not_positive_definite;
                res.structural_empty_at = empty_at;
                return res;
            }
        }

        // ---- 5. ordering (pattern-only pre-permutation P0, perm0[new] =
        // old) + bijection-verification firewall (D-13/SS7.5-1): a non-
        // bijective return from the ordering layer is stopped here as
        // internal_error before it can touch the factorization.
        std::vector<Index> perm0;
        if (!sparse_chol_detail::sparse_chol_compute_ordering_(
                n, col_ptr, row_ind, res.ordering_used, perm0)) {
            res.status = sparse_chol_status::internal_error;
            return res;
        }
        if (!sparse_chol_detail::sparse_chol_verify_permutation_(n, perm0)) {
            res.status = sparse_chol_status::internal_error;
            return res;
        }

        // ---- 6. explicit permuted lower CSC C = P0^T A P0 (SS2.0-4)
        std::vector<Index> C_col_ptr, C_row_ind;
        std::vector<T>     C_val;
        if (!sparse_chol_detail::sparse_chol_build_permuted_lower_csc_(
                n, col_ptr, row_ind, val, perm0, C_col_ptr, C_row_ind, C_val)) {
            res.status = sparse_chol_status::internal_error;
            return res;
        }

        // ---- 7. symbolic analysis (pattern-only)
        const sparse_chol_symbolic<Index> sym =
            sparse_chol_symbolic_analyze<Index>(n, C_col_ptr, C_row_ind);
        if (!sym.valid) {
            res.status = sparse_chol_status::internal_error;
            return res;
        }
        res.nnz_L = sym.nnz_L;   // symbolic value; stays valid on a stop

        // ---- 8. numeric stage (saved pattern reused; no re-ereach)
        sparse_chol_numeric_factorize<T, Index>(sym, C_col_ptr, C_row_ind, C_val, opt, res);
        if (res.status == sparse_chol_status::success) {
            res.perm = perm0;   // static: the ordering output itself (SS5.6)
        } else {
            res.perm.clear();   // D-3: no valid perm without a valid L
        }
        return res;

    } catch (const std::exception&) {
        // P3 final protection net: reaching here means a certified gate was
        // missed somewhere upstream; treated as a defect, not a normal path.
        res.status = sparse_chol_status::internal_error;
        return res;
    }
}

} // namespace vcp

#endif // VCP_TSPARSE_SPARSE_CHOL_HPP
