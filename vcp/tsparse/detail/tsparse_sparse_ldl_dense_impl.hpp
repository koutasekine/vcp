// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// LDL-0 -- dense fallback certified Bunch-Kaufman LDL^T kernel.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER the
// sparse_ldl_* type skeleton (options / result / status) and the SLU
// pivot_decision enum are in scope.  It has no "namespace vcp { }" wrapper;
// it is injected by tsparse_sparse_ldl.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_ldl.hpp>
//
// Role (design v2 SS8, LDL-0): permanent regression-reference kernel.
// LDL-2 verifies the sparse dynamic left-looking kernel against this one
// (identical pivot sequence, D-6), so this kernel and its signature are
// FROZEN after LDL-0 (design decision D-3).
//
// Layout (D-3): full dense column-major n x n work buffer; only the lower
// triangle (i >= j) is used; the strict upper triangle is dead storage.
//
// All gates are written in the certified three-way form (LDL-0 SS3.1,
// SLU-GT1 P1): the success-declaring branch keeps a certified comparison,
// the opposite branch requires the certified negation, and when neither can
// be certified the factorization stops with inconclusive_pivot_test
// (pivot_decision::inconclusive; division by an uncertified pivot is
// forbidden).  For totally ordered scalars (double / kv::dd / kv::mpfr<N>)
// the third branch is unreachable and the kernel behaves as classic BK.

#ifndef VCP_TSPARSE_SPARSE_LDL_DENSE_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LDL_DENSE_IMPL_HPP

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Dense certified Bunch-Kaufman factorization (lower, partial pivoting;
// LAPACK dsytrf family; design v2 SS1.3).
//
//   work : n*n column-major; on entry the lower triangle holds the values of
//          the (pre-permuted) symmetric input; on exit, columns of finalized
//          1x1 / 2x2 pivots hold the L factor below the pivot block.
//   perm : new->old permutation (design v2 SS5.4); BK symmetric exchanges
//          are composed into it in place, so passing the ordering
//          pre-permutation P0 yields the final perm = P0 o P_BK.
//   res  : receives L (CSC, explicit unit diagonal), the D triple
//          (D_diag / D_sub / D_block2), pivot counters, first_zero_pivot,
//          inconclusive_at, nnz_L and status (success / zero_pivot /
//          inconclusive_pivot_test only; the caller owns every other status).
//
// L / D / perm are valid only when the returned status is success or
// zero_pivot (zero pivots skip their column but the factorization runs to
// completion, design v2 SS1.3).  On inconclusive_pivot_test the
// factorization stops at column inconclusive_at and L is left empty.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_ldl_dense_bk_factorize(
    const Index n,
    std::vector<T>& work,
    std::vector<Index>& perm,
    const sparse_ldl_options<T>& opt,
    sparse_ldl_result<T, Index>& res)
{
    static_assert(std::is_signed<Index>::value, "sparse LDL Index must be signed");
    using std::abs;
    using std::sqrt;
    using std::swap;

    // All certified gates compare MAGNITUDES: R = real_type<T> is the ADL
    // abs return type (R == T for real scalars; double for complex, whose
    // instantiation must compile through the virtual policy _impl even
    // though LDL^H is out of scope).
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    const std::size_t un = static_cast<std::size_t>(n);

    // alpha = (1 + sqrt(17)) / 8, built in R (P2: unqualified ADL sqrt).
    const R alpha = (R(1) + sqrt(R(17))) / R(8);

    res.D_diag.assign(un, T(0));
    res.D_sub.assign(un, T(0));
    res.D_block2.assign(un, char(0));
    res.n_pivots_1x1 = Index(0);
    res.n_pivots_2x2 = Index(0);
    res.first_zero_pivot = Index(-1);
    res.inconclusive_at = Index(-1);
    res.nnz_L = Index(0);
    res.L_col_ptr.clear();
    res.L_row_ind.clear();
    res.L_val.clear();

    bool any_zero_pivot = false;
    bool stopped_inconclusive = false;

    // Scratch buffers for the original pivot-column values (allocated once).
    std::vector<T> wbuf1(un, T(0));
    std::vector<T> wbuf2(un, T(0));

#define VCP_LDL_W(i, j) work[(i) + (j) * un]

    // Symmetric row/column exchange p <-> q (p < q) in lower-triangular
    // storage, including the row relabel of already-finalized L columns
    // (uniform: for every column j < p both (p,j) and (q,j) are strictly
    // lower entries), plus the perm record.  W(q,p) is invariant.
    // (LDL-0 SS3.2 "対称交換".)

    std::size_t k = 0;
    while (k < un) {
        // ------ lambda = max_{i>k} |W(i,k)| : certified explicit loop (P6).
        // Update only when the candidate is certainly greater; when neither
        // "greater" nor "not greater" can be certified the whole search is
        // inconclusive and the factorization stops (design v2 SS1.3).
        R lam = R(0);
        std::size_t r = 0;
        bool r_valid = false;              // r_valid == false  <=>  lambda certified zero
        bool search_inconclusive = false;
        for (std::size_t i = k + 1; i < un; ++i) {
            const R a = abs(VCP_LDL_W(i, k));
            if (a > lam) { lam = a; r = i; r_valid = true; }
            else if (a <= lam) { /* certified: no update */ }
            else { search_inconclusive = true; break; }
        }
        if (search_inconclusive) {
            res.inconclusive_at = static_cast<Index>(k);
            stopped_inconclusive = true;
            break;
        }

        // ------ Bunch-Kaufman pivot decision (design v2 SS1.3 rules 1-4),
        // every rule in the certified three-way standard form (LDL-0 SS3.1).
        pivot_decision dec = pivot_decision::inconclusive;
        bool pivot_1x1 = false;   // valid when dec == acceptable
        bool swap_k_r = false;    // rule 3: symmetric exchange (k <-> r)

        const R absakk = abs(VCP_LDL_W(k, k));
        if (!r_valid) {
            // lambda certified zero -> rule 1 degenerates to |a_kk| >= 0,
            // which is certified for any abs value; take the 1x1 pivot and
            // let the division gate decide zero vs nonzero.
            dec = pivot_decision::acceptable;
            pivot_1x1 = true;
        } else {
            const R alam = alpha * lam;
            if (absakk >= alam) {                          // rule 1
                dec = pivot_decision::acceptable;
                pivot_1x1 = true;
            } else if (absakk < alam) {
                // sigma = max_{i in [k,n), i != r} |A(i,r)| (row r part read
                // transposed from lower storage): certified explicit loop.
                R sigma = R(0);
                bool sigma_inconclusive = false;
                for (std::size_t i = k; i < un; ++i) {
                    if (i == r) continue;
                    const R a = (i > r) ? abs(VCP_LDL_W(i, r)) : abs(VCP_LDL_W(r, i));
                    if (a > sigma) { sigma = a; }
                    else if (a <= sigma) { /* certified: no update */ }
                    else { sigma_inconclusive = true; break; }
                }
                if (!sigma_inconclusive) {
                    const R lhs2 = absakk * sigma;
                    const R rhs2 = alam * lam;
                    if (lhs2 >= rhs2) {                    // rule 2
                        dec = pivot_decision::acceptable;
                        pivot_1x1 = true;
                    } else if (lhs2 < rhs2) {
                        const R absarr = abs(VCP_LDL_W(r, r));
                        const R asig = alpha * sigma;
                        if (absarr >= asig) {              // rule 3
                            dec = pivot_decision::acceptable;
                            pivot_1x1 = true;
                            swap_k_r = true;
                        } else if (absarr < asig) {        // rule 4 (negations
                            dec = pivot_decision::acceptable;  //  1-3 certified)
                            pivot_1x1 = false;
                        }
                        // else: dec stays inconclusive
                    }
                    // else: dec stays inconclusive
                }
                // else: dec stays inconclusive
            }
            // else: dec stays inconclusive
        }

        if (dec == pivot_decision::inconclusive) {
            res.inconclusive_at = static_cast<Index>(k);
            stopped_inconclusive = true;
            break;
        }

        if (pivot_1x1) {
            // ------ optional rule-3 symmetric exchange (k <-> r).
            if (swap_k_r && r != k) {
                const std::size_t p = k, q = r;
                for (std::size_t j = 0; j < p; ++j) swap(VCP_LDL_W(p, j), VCP_LDL_W(q, j));
                for (std::size_t i = p + 1; i < q; ++i) swap(VCP_LDL_W(i, p), VCP_LDL_W(q, i));
                for (std::size_t i = q + 1; i < un; ++i) swap(VCP_LDL_W(i, p), VCP_LDL_W(i, q));
                swap(VCP_LDL_W(p, p), VCP_LDL_W(q, q));
                std::swap(perm[p], perm[q]);
            }

            // ------ 1x1 division gate (certified three-way; LDL-0 SS3.2).
            const T d = VCP_LDL_W(k, k);
            const R absd = abs(d);
            if (absd > opt.zero_pivot_tol) {
                // certified nonzero -> eliminate.
                res.D_diag[k] = d;
                for (std::size_t i = k + 1; i < un; ++i) {
                    wbuf1[i] = VCP_LDL_W(i, k);            // original a_ik
                    VCP_LDL_W(i, k) = wbuf1[i] / d;        // L(i,k)
                }
                for (std::size_t j = k + 1; j < un; ++j) {
                    const T wj = wbuf1[j];
                    for (std::size_t i = j; i < un; ++i) {
                        VCP_LDL_W(i, j) -= VCP_LDL_W(i, k) * wj;   // A(i,j) -= l_ik * a_jk
                    }
                }
                ++res.n_pivots_1x1;
            } else if (absd <= opt.zero_pivot_tol) {
                // certified zero pivot: D stays a structural zero, the column
                // elimination is skipped, factorization continues
                // (zero-eigenvalue counting use case, design v2 SS1.3).
                for (std::size_t i = k + 1; i < un; ++i) VCP_LDL_W(i, k) = T(0);
                if (!any_zero_pivot) {
                    res.first_zero_pivot = static_cast<Index>(k);
                    any_zero_pivot = true;
                }
            } else {
                // neither nonzero nor zero can be certified -> stop before
                // any division (P1/D3: fall before dividing).
                res.inconclusive_at = static_cast<Index>(k);
                stopped_inconclusive = true;
                break;
            }
            ++k;
        } else {
            // ------ rule 4: 2x2 pivot after symmetric exchange (k+1 <-> r).
            if (r != k + 1) {
                const std::size_t p = k + 1, q = r;
                for (std::size_t j = 0; j < p; ++j) swap(VCP_LDL_W(p, j), VCP_LDL_W(q, j));
                for (std::size_t i = p + 1; i < q; ++i) swap(VCP_LDL_W(i, p), VCP_LDL_W(q, i));
                for (std::size_t i = q + 1; i < un; ++i) swap(VCP_LDL_W(i, p), VCP_LDL_W(i, q));
                swap(VCP_LDL_W(p, p), VCP_LDL_W(q, q));
                std::swap(perm[p], perm[q]);
            }

            const T d1 = VCP_LDL_W(k, k);
            const T e  = VCP_LDL_W(k + 1, k);
            const T d2 = VCP_LDL_W(k + 1, k + 1);
            const T det = d1 * d2 - e * e;
            const R absdet = abs(det);
            // 2x2 division gate: tol2 = opt.zero_pivot_tol (D-1; default 0 =
            // only a certified exact zero is treated as a zero block).
            if (absdet > opt.zero_pivot_tol) {
                res.D_diag[k] = d1;
                res.D_diag[k + 1] = d2;
                res.D_sub[k] = e;
                res.D_block2[k] = char(1);
                for (std::size_t i = k + 2; i < un; ++i) {
                    wbuf1[i] = VCP_LDL_W(i, k);            // original a_ik
                    wbuf2[i] = VCP_LDL_W(i, k + 1);        // original a_i,k+1
                    VCP_LDL_W(i, k)     = (d2 * wbuf1[i] - e * wbuf2[i]) / det;
                    VCP_LDL_W(i, k + 1) = (d1 * wbuf2[i] - e * wbuf1[i]) / det;
                }
                for (std::size_t j = k + 2; j < un; ++j) {
                    const T wj1 = wbuf1[j];
                    const T wj2 = wbuf2[j];
                    for (std::size_t i = j; i < un; ++i) {
                        // A(i,j) -= [l_i1 l_i2] * B * [l_j1; l_j2]
                        //         = l_i1 * a_jk + l_i2 * a_j,k+1
                        VCP_LDL_W(i, j) -= VCP_LDL_W(i, k) * wj1 + VCP_LDL_W(i, k + 1) * wj2;
                    }
                }
                ++res.n_pivots_2x2;
            } else if (absdet <= opt.zero_pivot_tol) {
                // certified zero 2x2 block: skip both columns as a structural
                // zero block and continue (status zero_pivot).  Unreachable
                // for exact arithmetic on a genuine rule-4 block
                // (|det| >= (1 - alpha^2) * lambda^2 > 0) but kept honest.
                for (std::size_t i = k + 1; i < un; ++i) VCP_LDL_W(i, k) = T(0);
                for (std::size_t i = k + 2; i < un; ++i) VCP_LDL_W(i, k + 1) = T(0);
                res.D_block2[k] = char(1);
                if (!any_zero_pivot) {
                    res.first_zero_pivot = static_cast<Index>(k);
                    any_zero_pivot = true;
                }
            } else {
                res.inconclusive_at = static_cast<Index>(k);
                stopped_inconclusive = true;
                break;
            }
            k += 2;
        }
    }

    if (stopped_inconclusive) {
        // Interrupted factorization: L is deliberately left empty; L / D /
        // perm are not valid outputs under inconclusive_pivot_test.
        res.status = sparse_ldl_status::inconclusive_pivot_test;
        return;
    }

    // ------ L extraction: CSC with explicit unit diagonal (design v2 SS5.2).
    // A strictly-lower work entry is dropped only when it is CERTIFIED equal
    // to T(0); a value whose zero-ness cannot be certified is stored (LDL-0
    // SS3.3).  For a 2x2 pivot block the (k+1,k) work slot holds the D
    // subdiagonal e_k, which belongs to D, not to L: L(k+1,k) = 0
    // structurally, so extraction of a block-leading column starts at k+2.
    res.L_col_ptr.assign(un + 1u, Index(0));
    for (std::size_t j = 0; j < un; ++j) {
        res.L_row_ind.push_back(static_cast<Index>(j));   // explicit unit diagonal
        res.L_val.push_back(T(1));
        const std::size_t first = (res.D_block2[j] != char(0)) ? (j + 2) : (j + 1);
        for (std::size_t i = first; i < un; ++i) {
            const T x = VCP_LDL_W(i, j);
            if (x == T(0)) { /* certified zero: not stored */ }
            else {
                res.L_row_ind.push_back(static_cast<Index>(i));
                res.L_val.push_back(x);
            }
        }
        res.L_col_ptr[j + 1u] = static_cast<Index>(res.L_row_ind.size());
    }
    res.nnz_L = static_cast<Index>(res.L_row_ind.size());

    res.status = any_zero_pivot ? sparse_ldl_status::zero_pivot
                                : sparse_ldl_status::success;
#undef VCP_LDL_W
}

#endif // VCP_TSPARSE_SPARSE_LDL_DENSE_IMPL_HPP
