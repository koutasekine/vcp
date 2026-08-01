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

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <type_traits>
#include <utility>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>
// Dense block kernels of the supernodal numeric phase (SLDL-SP SP-1).  The
// include MUST stay here, at file scope: the detail impls below are injected
// from INSIDE namespace vcp and tblas.hpp opens namespace vcp itself.
#include <vcp/tblas/tblas.hpp>

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
    baseline_dynamic,
    supernodal      // SLDL-SP (Phase B / B1): left-looking supernodal panels
    // still open: multifrontal
};

// SLDL-SP design v1 SS3.3: the pivot strategy is ORTHOGONAL to the method.
//   bk   : dynamic Bunch-Kaufman (1x1 / 2x2, symmetric exchanges) -- default,
//          i.e. the pre-SLDL-SP behaviour of both baseline and supernodal.
//   none : static, exchange-free.  Diagonal 1x1 pivots in their natural
//          order; a certified zero pivot is skipped and counted (the
//          zero-eigenvalue counting use case); a pivot that can be certified
//          neither zero nor nonzero stops with inconclusive_pivot_test.
// none exists for both methods: baseline x none is the scalar reference
// system that supernodal x none is verified against (design v1 SS6-2).
enum class sparse_ldl_pivoting {
    bk,
    none
};

// SLDL-SP design v1 D-3: which dense kernel computes the symmetric diagonal
// block update C -= W * L1^T of a supernode panel.  auto_select is resolved
// by the library and always reported in diag_kernel_used.
enum class sparse_ldl_diag_kernel {
    auto_select,   // resolved to gemmtr in v1 (SP-3 calibrates the default)
    gemmtr,        // tgemmtr: triangular part only, half the flops
    gemm           // tgemm: full block, the faster double kernel
};

// SLDL-AM design v1.1 SS2: relaxed amalgamation of the supernodal SYMBOLIC
// partition.  off (the default until ruling AM-D2) keeps the fundamental
// partition and is byte-identical to the pre-SLDL-AM behaviour: the merge
// pass is not even entered.  relaxed merges adjacent parent-child supernode
// pairs under the staged rule of SS2 (constants below); the merged panels
// carry EXPLICIT ZEROS that are exact 0.0 through the numeric phase and are
// dropped by the output emission, so L/D patterns and nnz_L are unchanged --
// only the summation order (and hence ulps) of true-pattern values moves.
enum class sparse_ldl_amalgamation {
    off,
    relaxed
};

// What the numeric phase actually used (result diagnostic).
enum class sparse_ldl_kernel_used {
    not_applicable,   // no supernodal numeric work was performed
    scalar,           // scalar mirror (type gate not satisfied)
    gemmtr,
    gemm
};

enum class sparse_ldl_ordering {
    auto_select,
    natural,
    rcm,
    amd,
    nested_dissection,
    // ORD-1 (pure addition, ruling D-1): multilevel nested dissection
    // (dependency-free, METIS-class target; tsparse_order_ndml_impl.hpp).
    // The existing nested_dissection is kept unchanged for reproducibility.
    nested_dissection_ml
    // colamd is intentionally absent at the type level (A^T A graph is for
    // nonsymmetric LU; design v2 SS3).
};

enum class sparse_ldl_status {
    success,
    structural_singularity,   // structurally empty row/column (assembly-bug hint)
    zero_pivot,               // certified zero pivot(s); factorization completed
    inconclusive_pivot_test,  // pivot decision could not be certified; stopped (P1)
    pivot_out_of_panel,       // SLDL-SP D-1: the supernodal panel needed a
                              // symmetric exchange with a column OUTSIDE the
                              // panel.  The BK test is never weakened and no
                              // silent fallback is taken: the factorization
                              // stops and returns diagnostics only.  Callers
                              // may re-run with method = baseline_dynamic,
                              // which handles exchanges at any distance.
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
    case sparse_ldl_status::pivot_out_of_panel:      return "pivot_out_of_panel";
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
    sparse_ldl_pivoting pivoting;  // default bk (the pre-SLDL-SP behaviour)
    // Type gate of the supernodal blocked kernels (design v1 SS3.2): dense
    // block operations are used only for std::is_floating_point scalars AND
    // only when the descendant supernode is at least this wide; every other
    // case runs the scalar mirror.  The default is provisional and is
    // calibrated in SP-3.
    int ldl_min_block_size;
    sparse_ldl_diag_kernel diag_kernel;
    // SLDL-AM: relaxed amalgamation switch + staged-rule constants (design
    // v1.1 SS2, ruling AM-D1).  A merged candidate of width w, panel rows m
    // and zero PERCENTAGE z (= 100 * padded zeros / (m*w), integers only --
    // P4) is merged iff
    //     w <= n0                        (unconditional small-supernode rescue)
    //  or (w <= n1 and z < z0)
    //  or (w <= n2 and z < z1)
    //  or (z < z2).
    // The default constants are the CHOLMOD-equivalent (4,16,48 / 80,10,5).
    // ldl_amalgamation stays off until ruling AM-D2 (default-ON is a separate
    // one-line owner-decided commit).
    sparse_ldl_amalgamation ldl_amalgamation;
    int ldl_amalg_n0, ldl_amalg_n1, ldl_amalg_n2;   // stage widths
    int ldl_amalg_z0, ldl_amalg_z1, ldl_amalg_z2;   // integer percentages
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
    // ORD-1 (appended; additive): parameters of ordering =
    // nested_dissection_ml, read ONLY on that ordering.  A negative value
    // means "use the library default" (sparse_order_ndml_params, the single
    // OR-2 calibration source); every other ordering ignores these fields.
    long long ndml_coarsen_stop;
    int       ndml_fm_passes;
    int       ndml_balance_pct;
    long long ndml_leaf_size;

    sparse_ldl_options()
        : method(sparse_ldl_method::auto_select),
          pivoting(sparse_ldl_pivoting::bk),
          // provisional default (implementation guide SS2-2): the LU value 16
          // is a syrk-shaped threshold and is deliberately NOT inherited;
          // SP-3 calibrates this on four machines.
          ldl_min_block_size(32),
          diag_kernel(sparse_ldl_diag_kernel::auto_select),
          ldl_amalgamation(sparse_ldl_amalgamation::off),
          ldl_amalg_n0(4), ldl_amalg_n1(16), ldl_amalg_n2(48),
          ldl_amalg_z0(80), ldl_amalg_z1(10), ldl_amalg_z2(5),
          ordering(sparse_ldl_ordering::auto_select),
          check_symmetry(true),
          // same default policy as policy_is_symmetric (B-4: via the D4
          // customization point, not numeric_limits).
          symmetry_tol(vcp::tsparse_scalar::decimal_power_negative<real_type>(12u)),
          zero_pivot_tol(real_type(0)),
          pivot_threshold(real_type(0)),
          ndml_coarsen_stop(-1), ndml_fm_passes(-1),
          ndml_balance_pct(-1), ndml_leaf_size(-1) {}
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

    // ---- SLDL-SP diagnostics (design v1 SS4; integers and enums only, P4).
    sparse_ldl_pivoting    pivot_mode_used;
    sparse_ldl_kernel_used diag_kernel_used;   // not_applicable outside supernodal
    Index n_supernodes;          // supernodes of the symbolic partition
    Index max_supernode_width;
    Index n_boundary_splits;     // 2x2 pivots that split a supernode (D-4)
    Index nnz_L_static;          // entries of the static (exchange-free) factor
    // Certified zero pivots that were skipped and counted.  Maintained by the
    // static mode and by the supernodal kernel in both pivot modes.  The
    // FROZEN baseline x bk path does not maintain it -- its code path is
    // unchanged by contract (implementation guide SS0) -- so read
    // first_zero_pivot / status there instead of this counter.
    Index n_zero_skips;
    Index out_of_panel_at;       // column that needed an out-of-panel exchange; -1 = none
    long long gemm_call_count;   // dense block calls issued by the numeric phase
    long long gemm_time_ns;      // time spent in them; valid only when
                                 // gemm_call_count > 0 (never measured for
                                 // non-floating scalars -- R-6)
    // Growth diagnostic (D-5): floor(log2(max|d| / max|a_ii|)), clipped to
    // +-1024.  growth_valid is the ONLY validity criterion (P4: no numeric
    // sentinel); it stays false when the input diagonal is not certified
    // nonzero or when the mode does not track growth.
    int  growth_log2;
    bool growth_valid;

    // ---- SLDL-AM diagnostics (design v1.1 SS2.1 / SS5-7; integers only, P4).
    // Filled from the symbolic analysis by the numeric driver for the
    // supernodal method; they keep their defaults on the baseline paths.
    sparse_ldl_amalgamation amalgamation_used;
    Index     n_amalgamations;            // merges performed by the relax pass
    long long n_padded_zeros;             // explicit zeros of the padded pattern
    long long mean_supernode_width_x100;  // 100*n / n_supernodes (integer)
    // Exact flop estimates Sum_snode h*w*(w+h) (h = below-diagonal rows,
    // w = width), as 128-bit values split into (hi,lo) 64-bit halves
    // (SS2.1: the >64-bit pair).  flops_true is the sum over the FUNDAMENTAL
    // partition (= the relax=off work); flops_padded is the sum over the
    // partition the numeric phase actually runs (merged when relax=on, equal
    // to flops_true when off).  flop_pad_ratio = flops_padded / flops_true.
    unsigned long long flops_padded_lo, flops_padded_hi;
    unsigned long long flops_true_lo,   flops_true_hi;

    sparse_ldl_result()
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
          growth_log2(0), growth_valid(false),
          amalgamation_used(sparse_ldl_amalgamation::off),
          n_amalgamations(Index(0)), n_padded_zeros(0),
          mean_supernode_width_x100(0),
          flops_padded_lo(0u), flops_padded_hi(0u),
          flops_true_lo(0u), flops_true_hi(0u) {}
};

// ---------------------------------------------------------------------------
// symbolic phase (SLDL-SP SP-0): ordering, pattern, etree, postorder, column
// counts, fundamental supernodes, workspace sizes.  T-independent.
// ---------------------------------------------------------------------------
#include <vcp/tsparse/detail/tsparse_sparse_ldl_symbolic_impl.hpp>

// ===========================================================================
// numeric phase (SLDL-SP SP-0): workspace, kernels, and the symbolic-driven
// numeric entry
// ===========================================================================

// Reusable numeric scratch (design v1 SS4 / B2 boundary requirement 3).  The
// buffers are owned by the caller so that a B2-style shift loop can factor
// A - sigma*B repeatedly without reallocating.  The baseline kernels keep
// their own internal state and ignore this object; it is used and sized by
// the supernodal kernel.
template <class T, class Index>
struct sparse_ldl_numeric_workspace {
    static_assert(std::is_signed<Index>::value,
                  "sparse LDL Index must be signed");

    std::vector<T>     panel;      // dense supernode panel (column-major)
    std::vector<T>     panel_a;    // W = L*D of the current descendant update
    std::vector<T>     update;     // dense update block from one descendant
    std::vector<Index> relind;     // matrix row -> position inside the panel
    std::vector<Index> iwork;      // general integer scratch

    void clear_buffers() {
        panel.clear(); panel_a.clear(); update.clear();
        relind.clear(); iwork.clear();
    }
};

// ---------------------------------------------------------------------------
// dense fallback certified BK kernel (LDL-0; frozen regression reference)
// ---------------------------------------------------------------------------
#include <vcp/tsparse/detail/tsparse_sparse_ldl_dense_impl.hpp>

// ---------------------------------------------------------------------------
// sparse dynamic left-looking certified BK kernel (LDL-2) + the static
// (pivoting = none) entry branch added by SLDL-SP SP-1
// ---------------------------------------------------------------------------
#include <vcp/tsparse/detail/tsparse_sparse_ldl_sparse_impl.hpp>

// ---------------------------------------------------------------------------
// left-looking supernodal numeric kernel (SLDL-SP SP-1)
// ---------------------------------------------------------------------------
#include <vcp/tsparse/detail/tsparse_sparse_ldl_supernodal_impl.hpp>

// ---------------------------------------------------------------------------
// sparse_ldl_factorize_numeric_with_info -- numeric phase driven by a symbolic
// result (design v1 SS4: the one-shot API is a DELEGATION to this pair, never
// an independent implementation).
//
// Preconditions (checked, reported through res.status):
//   - sym.status == success and sym.n == n,
//   - the CSC triple is the same one the symbolic phase analysed.
// res.method_used / res.ordering_used must already be resolved by the caller;
// this function does not re-resolve auto_select.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_ldl_factorize_numeric_with_info(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const sparse_ldl_symbolic_result<Index>& sym,
    const sparse_ldl_options<T>& opt,
    sparse_ldl_numeric_workspace<T, Index>& ws,
    sparse_ldl_result<T, Index>& res)
{
    static_assert(std::is_signed<Index>::value, "sparse LDL Index must be signed");

    if (sym.status != sparse_ldl_symbolic_status::success || sym.n != n ||
        sym.perm0.size() != static_cast<std::size_t>(n)) {
        res.status = sparse_ldl_status::internal_error;
        return;
    }
    const std::size_t un = static_cast<std::size_t>(n);
    res.pivot_mode_used = opt.pivoting;

    // ---- supernodal (SLDL-SP SP-1): needs the full symbolic analysis.
    if (res.method_used == sparse_ldl_method::supernodal) {
        if (sym.level != sparse_ldl_symbolic_level::full) {
            res.status = sparse_ldl_status::internal_error;
            return;
        }
        // ---- SLDL-AM bridge (design v2 SS3-4 / SS2.0).  The B2 shift
        // handle builds its symbolic analysis with default (off)
        // amalgamation options -- its setup code is frozen and copies only
        // ordering and level -- so an explicit relax=on request must be
        // honoured HERE: the relaxed analysis (postorder relabelling BEFORE
        // fundamental detection, then the merge pass) is recomputed locally
        // by the same deterministic sparse_ldl_symbolic_analyze the one-shot
        // entry runs, on the same CSC pattern and resolved ordering.  Both
        // sides therefore hand the kernel identical arrays, which is what
        // makes acceptance 4 (factor_at == one-shot, byte identity) hold.
        // Cost note: this path pays one full symbolic analysis (ordering
        // included) per numeric call; it exists only for the off-built-
        // analysis + relax=on combination.  The off path below remains a
        // plain pass-through that enters neither relabelling nor merging.
        if (opt.ldl_amalgamation == sparse_ldl_amalgamation::relaxed &&
            sym.amalgamation_used == sparse_ldl_amalgamation::off) {
            sparse_ldl_symbolic_options sopt2;
            sopt2.ordering = sym.ordering_used;   // already resolved
            sopt2.level = sparse_ldl_symbolic_level::full;
            sopt2.amalgamation = sparse_ldl_amalgamation::relaxed;
            sopt2.amalg_n0 = opt.ldl_amalg_n0;
            sopt2.amalg_n1 = opt.ldl_amalg_n1;
            sopt2.amalg_n2 = opt.ldl_amalg_n2;
            sopt2.amalg_z0 = opt.ldl_amalg_z0;
            sopt2.amalg_z1 = opt.ldl_amalg_z1;
            sopt2.amalg_z2 = opt.ldl_amalg_z2;
            // ORD-1: keep the ndml parameters aligned with the one-shot
            // analysis so the recomputed ordering is identical.
            sopt2.ndml_coarsen_stop = opt.ndml_coarsen_stop;
            sopt2.ndml_fm_passes    = opt.ndml_fm_passes;
            sopt2.ndml_balance_pct  = opt.ndml_balance_pct;
            sopt2.ndml_leaf_size    = opt.ndml_leaf_size;
            const sparse_ldl_symbolic_result<Index> relaxed_sym =
                sparse_ldl_symbolic_analyze(n, col_ptr, row_ind, sopt2);
            if (relaxed_sym.status != sparse_ldl_symbolic_status::success ||
                relaxed_sym.ordering_used != sym.ordering_used) {
                res.status = sparse_ldl_status::internal_error;
                return;
            }
            sparse_ldl_supernodal_factorize(
                n, col_ptr, row_ind, val, relaxed_sym, opt, ws, res);
            res.amalgamation_used         = relaxed_sym.amalgamation_used;
            res.n_amalgamations           = relaxed_sym.n_amalgamations;
            res.n_padded_zeros            = relaxed_sym.n_padded_zeros;
            res.mean_supernode_width_x100 = relaxed_sym.mean_supernode_width_x100;
            res.flops_padded_lo = relaxed_sym.flops_padded_lo;
            res.flops_padded_hi = relaxed_sym.flops_padded_hi;
            res.flops_true_lo   = relaxed_sym.flops_true_lo;
            res.flops_true_hi   = relaxed_sym.flops_true_hi;
            return;
        }
        if (opt.ldl_amalgamation == sparse_ldl_amalgamation::off &&
            sym.amalgamation_used == sparse_ldl_amalgamation::relaxed) {
            // A merged analysis cannot be un-merged (the fundamental
            // partition is gone); honest death instead of a silently
            // different structure.
            res.status = sparse_ldl_status::internal_error;
            return;
        }
        sparse_ldl_supernodal_factorize(n, col_ptr, row_ind, val, sym, opt, ws, res);
        res.amalgamation_used         = sym.amalgamation_used;
        res.n_amalgamations           = sym.n_amalgamations;
        res.n_padded_zeros            = sym.n_padded_zeros;
        res.mean_supernode_width_x100 = sym.mean_supernode_width_x100;
        res.flops_padded_lo = sym.flops_padded_lo;
        res.flops_padded_hi = sym.flops_padded_hi;
        res.flops_true_lo   = sym.flops_true_lo;
        res.flops_true_hi   = sym.flops_true_hi;
        return;
    }

    // ---- baseline_dynamic: the frozen path.
    // The dense fallback kernel is FROZEN (LDL-0 design decision D-3) and
    // implements the dynamic BK rules only, so the small-n delegation applies
    // to pivoting = bk alone; the static mode always runs the sparse kernel,
    // which reports dense_delegated = false honestly.  bk therefore keeps its
    // byte-identical behaviour at every n.
    if (un <= 64u && opt.pivoting == sparse_ldl_pivoting::bk) {
        // CSC -> dense scatter of the lower triangle (strictly upper entries
        // are ignored, design v2 SS1.2) under P0: element (i,j) lands at
        // (pinv[i], pinv[j]) of P0^T A P0, on the lower side.
        std::vector<T> work(un * un, T(0));
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r >= static_cast<Index>(c)) {
                    const std::size_t a = static_cast<std::size_t>(sym.pinv0[static_cast<std::size_t>(r)]);
                    const std::size_t b = static_cast<std::size_t>(sym.pinv0[c]);
                    if (a >= b) work[a + b * un] = val[static_cast<std::size_t>(k)];
                    else        work[b + a * un] = val[static_cast<std::size_t>(k)];
                }
            }
        }
        // BK exchanges are composed into perm in place, so seeding with P0
        // yields the final perm = P0 o P_BK (new->old, design v2 SS5.4).
        res.perm = sym.perm0;
        res.dense_delegated = true;
        sparse_ldl_dense_bk_factorize(n, work, res.perm, opt, res);
    } else {
        res.dense_delegated = false;
        sparse_ldl_sparse_bk_factorize(n, col_ptr, row_ind, val, sym.perm0, opt, res);
    }
}

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
        case sparse_ldl_method::baseline_dynamic:
        case sparse_ldl_method::supernodal:
            // auto_select resolution is centralised in resolve_auto_method
            // (design v1 SS4); v1 keeps auto -> baseline_dynamic (switching
            // it is the separate SP-D1 track).
            res.method_used = sparse_ldl_detail::resolve_auto_method(opt.method);
            break;
        default:
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }
        switch (opt.pivoting) {
        case sparse_ldl_pivoting::bk:
        case sparse_ldl_pivoting::none:
            res.pivot_mode_used = opt.pivoting;
            break;
        default:
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }
        switch (opt.diag_kernel) {
        case sparse_ldl_diag_kernel::auto_select:
        case sparse_ldl_diag_kernel::gemmtr:
        case sparse_ldl_diag_kernel::gemm:
            break;
        default:
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }
        if (opt.ldl_min_block_size < 1) {
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }
        switch (opt.ldl_amalgamation) {
        case sparse_ldl_amalgamation::off:
        case sparse_ldl_amalgamation::relaxed:
            break;
        default:
            res.status = sparse_ldl_status::invalid_options;
            return res;
        }
        if (opt.ldl_amalg_n0 < 0 || opt.ldl_amalg_n1 < 0 || opt.ldl_amalg_n2 < 0 ||
            opt.ldl_amalg_z0 < 0 || opt.ldl_amalg_z0 > 100 ||
            opt.ldl_amalg_z1 < 0 || opt.ldl_amalg_z1 > 100 ||
            opt.ldl_amalg_z2 < 0 || opt.ldl_amalg_z2 > 100) {
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
        case sparse_ldl_ordering::nested_dissection_ml:
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

        // ---- 5. symbolic phase (SLDL-SP SP-0): the pattern-only
        // pre-permutation P0 (perm0[new] = old) and -- at analysis level
        // full -- the etree / column counts / supernode structure.  The
        // ordering functions themselves are the existing SLU ones, reused by
        // include only and unmodified.
        sparse_ldl_symbolic_options sopt;
        sopt.ordering = opt.ordering;
        // baseline_dynamic exchanges rows and columns dynamically, so the
        // static structure has no meaning for it: analysing only the ordering
        // keeps the frozen path free of the O(nnz_L) symbolic cost.  The
        // supernodal method needs the full analysis.
        sopt.level = (res.method_used == sparse_ldl_method::supernodal)
                   ? sparse_ldl_symbolic_level::full
                   : sparse_ldl_symbolic_level::ordering_only;
        // SLDL-AM: the one-shot entry merges at the SYMBOLIC side; the
        // numeric driver's bridge then sees a matching analysis and passes
        // it through unchanged (effective only at level full).
        sopt.amalgamation = opt.ldl_amalgamation;
        sopt.amalg_n0 = opt.ldl_amalg_n0;
        sopt.amalg_n1 = opt.ldl_amalg_n1;
        sopt.amalg_n2 = opt.ldl_amalg_n2;
        sopt.amalg_z0 = opt.ldl_amalg_z0;
        sopt.amalg_z1 = opt.ldl_amalg_z1;
        sopt.amalg_z2 = opt.ldl_amalg_z2;
        // ORD-1: nested_dissection_ml parameters (ignored by every other
        // ordering; negative = library default).
        sopt.ndml_coarsen_stop = opt.ndml_coarsen_stop;
        sopt.ndml_fm_passes    = opt.ndml_fm_passes;
        sopt.ndml_balance_pct  = opt.ndml_balance_pct;
        sopt.ndml_leaf_size    = opt.ndml_leaf_size;
        const sparse_ldl_symbolic_result<Index> sym =
            sparse_ldl_symbolic_analyze(n, col_ptr, row_ind, sopt);
        if (sym.status != sparse_ldl_symbolic_status::success) {
            res.status = (sym.status == sparse_ldl_symbolic_status::invalid_options)
                       ? sparse_ldl_status::invalid_options
                       : ((sym.status == sparse_ldl_symbolic_status::invalid_input)
                          ? sparse_ldl_status::invalid_input
                          : sparse_ldl_status::internal_error);
            return res;
        }
        // The ordering resolution lives in the symbolic phase now; the entry
        // check above rejected out-of-enum values, so the two must agree.
        if (sym.ordering_used != res.ordering_used) {
            res.status = sparse_ldl_status::internal_error;
            return res;
        }

        // ---- 6. numeric phase (LDL-2 kernels), by delegation.
        sparse_ldl_numeric_workspace<T, Index> ws;
        sparse_ldl_factorize_numeric_with_info(
            n, col_ptr, row_ind, val, sym, opt, ws, res);
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
