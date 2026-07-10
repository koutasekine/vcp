// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// LDL-2 -- sparse dynamic left-looking certified Bunch-Kaufman LDL^T kernel.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER the
// sparse_ldl_* type skeleton and the dense kernel (which defines the frozen
// pivot rules this kernel replicates) are in scope.  It has no
// "namespace vcp { }" wrapper; it is injected by tsparse_sparse_ldl.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_ldl.hpp>
//
// Contract (D-6): this kernel reproduces the dense kernel's PIVOT SEQUENCE
// (1x1/2x2 kind sequence, exchange columns, perm) exactly.  The
// implementation goes further and performs, for every matrix element, the
// same floating-point operations in the same order as the dense kernel:
//  - working column c is built left-looking as A(:,c) minus the committed
//    pivot updates in ascending pivot order (the order the dense
//    right-looking kernel applied them);
//  - each 1x1 update term is L(i,j) * abar_cj with abar the PRE-division
//    reduced value (the dense kernel's wbuf coefficient), never
//    L(i,j)*d*L(c,j) (which would round differently);
//  - rows ABOVE the target column (needed for the sigma search / d2) use the
//    transposed form L(c,j) * abar_ij, matching the dense element (c,i);
//  - 2x2 updates subtract l1*w1 + l2*w2 as one expression, as dense does;
//  - lambda/sigma searches iterate the working pattern in ascending row
//    order (certified explicit loops, P6), so certified ties resolve to the
//    same row as the dense ascending scan; absent (structurally zero)
//    entries can never update a max whose candidates are certified >= 0.
//
// Occurrence-position stacks (SLU-SNQ lesson, G-L2.2): the dense working
// vectors are reset through their touched-row stacks only; there is no O(n)
// clear or O(n) scan per column anywhere in this kernel.
//
// Dynamic symmetric exchange: label maps (current <-> assembly labels) make
// the uncommitted-part exchange O(1), and the committed-L row relabel walks
// only the two row lists ((block,pos) links), never a full row copy.

#ifndef VCP_TSPARSE_SPARSE_LDL_SPARSE_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LDL_SPARSE_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace sparse_ldl_detail {

// One committed pivot block (width 1 or 2).  Column labels are pivot
// positions and never change after commit; row labels are CURRENT labels and
// are relabeled through row lists on symmetric exchange.
template <class T, class Index>
struct sparse_ldl_committed_block {
    Index jcol;               // leading pivot column
    int   width;              // 1 or 2
    std::vector<Index> rows;  // strictly-below-block rows (current labels), ascending at commit
    std::vector<T> l1, a1;    // L values and pre-division reduced values, column jcol
    std::vector<T> l2, a2;    // same for column jcol+1 (width 2 only)
};

} // namespace sparse_ldl_detail

// ---------------------------------------------------------------------------
// Sparse dynamic left-looking certified BK factorization.
//
//   col_ptr/row_ind/val : validated CSC; only the lower triangle (i >= j) is
//                         read (design v2 SS1.2).
//   perm0               : new->old pre-permutation P0 (ordering).  BK
//                         exchanges are composed on top; res.perm receives
//                         the final perm = P0 o P_BK.
//   res                 : same output contract as sparse_ldl_dense_bk_factorize.
// ---------------------------------------------------------------------------
template <class T, class Index>
void sparse_ldl_sparse_bk_factorize(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const std::vector<T>&     val,
    const std::vector<Index>& perm0,
    const sparse_ldl_options<T>& opt,
    sparse_ldl_result<T, Index>& res)
{
    static_assert(std::is_signed<Index>::value, "sparse LDL Index must be signed");
    using std::abs;
    using std::sqrt;

    typedef sparse_ldl_detail::sparse_ldl_committed_block<T, Index> block_t;
    // R = real_type<T>: magnitude type of the certified gates (see the dense
    // kernel note; identical rule so the pivot decisions stay in lockstep).
    typedef typename vcp::tsparse_scalar::real_type<T>::type R;

    const std::size_t un = static_cast<std::size_t>(n);
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
    res.perm = perm0;

    // ---- full symmetric adjacency (values) in assembly labels = P0-permuted
    // labels.  Each stored lower entry (r,c,v) appears as (a->(b,v)) and, off
    // the diagonal, (b->(a,v)) with a = pinv[r], b = pinv[c].
    std::vector<Index> pinv(un);
    for (std::size_t i = 0; i < un; ++i) pinv[static_cast<std::size_t>(perm0[i])] = static_cast<Index>(i);
    std::vector<Index> adj_ptr(un + 1u, Index(0));
    {
        std::vector<Index> deg(un, Index(0));
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r >= static_cast<Index>(c)) {
                    const std::size_t a = static_cast<std::size_t>(pinv[static_cast<std::size_t>(r)]);
                    const std::size_t b = static_cast<std::size_t>(pinv[c]);
                    ++deg[b];
                    if (a != b) ++deg[a];
                }
            }
        }
        for (std::size_t i = 0; i < un; ++i) adj_ptr[i + 1u] = adj_ptr[i] + deg[i];
    }
    const std::size_t adj_nnz = static_cast<std::size_t>(adj_ptr[un]);
    std::vector<Index> adj_ind(adj_nnz);
    std::vector<T>     adj_val(adj_nnz, T(0));
    {
        std::vector<Index> head(adj_ptr.begin(), adj_ptr.end() - 1);
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r >= static_cast<Index>(c)) {
                    const std::size_t a = static_cast<std::size_t>(pinv[static_cast<std::size_t>(r)]);
                    const std::size_t b = static_cast<std::size_t>(pinv[c]);
                    const T& v = val[static_cast<std::size_t>(k)];
                    std::size_t pos = static_cast<std::size_t>(head[b]++);
                    adj_ind[pos] = static_cast<Index>(a);
                    adj_val[pos] = v;
                    if (a != b) {
                        pos = static_cast<std::size_t>(head[a]++);
                        adj_ind[pos] = static_cast<Index>(b);
                        adj_val[pos] = v;
                    }
                }
            }
        }
    }

    // ---- label maps: assembly label <-> current label (O(1) exchange of the
    // uncommitted part).
    std::vector<Index> cur_of_asm(un), asm_of_cur(un);
    for (std::size_t i = 0; i < un; ++i) { cur_of_asm[i] = static_cast<Index>(i); asm_of_cur[i] = static_cast<Index>(i); }

    // ---- committed blocks + row lists ((block id, position) links)
    std::vector<block_t> blocks;
    blocks.reserve(un);
    std::vector<std::vector<std::pair<Index, Index> > > rowlist(un);
    // col_block_of[j]: block id covering pivot column j, -1 = zero-skip column
    std::vector<Index> col_block_of(un, Index(-1));

    // ---- two working columns with occurrence-position stacks
    std::vector<T> wa(un, T(0)), wb(un, T(0));
    std::vector<char> ma(un, char(0)), mb(un, char(0));
    std::vector<Index> sa, sb;
    sa.reserve(un); sb.reserve(un);

    bool any_zero_pivot = false;
    bool stopped_inconclusive = false;

    // ------ helpers (C++11 lambdas; integer sorts only, P6) ------
    // reset a working column through its stack (no O(n) clear)
#define VCP_LDL_RESET(w, m, s) \
    do { for (std::size_t t_ = 0; t_ < (s).size(); ++t_) { (m)[static_cast<std::size_t>((s)[t_])] = char(0); (w)[static_cast<std::size_t>((s)[t_])] = T(0); } (s).clear(); } while (0)

#define VCP_LDL_TOUCH(w, m, s, i) \
    do { if (!(m)[i]) { (m)[i] = char(1); (w)[i] = T(0); (s).push_back(static_cast<Index>(i)); } } while (0)

    // Build the working column for current column c_target, rows >= cutoff:
    //   w = A(:,c_target) (reduced) with, for rows i < c_target, the
    //   transposed element (c_target, i) -- bitwise the dense kernel values.
    // Pattern stack is left SORTED ascending.
    // (declared as a lambda so it can capture the kernel state)
    const std::size_t UN = un;
    auto build_column = [&](const std::size_t c_target, const std::size_t cutoff,
                            std::vector<T>& w, std::vector<char>& m, std::vector<Index>& s) {
        VCP_LDL_RESET(w, m, s);
        // gather A entries
        const std::size_t alab = static_cast<std::size_t>(asm_of_cur[c_target]);
        for (Index k = adj_ptr[alab]; k < adj_ptr[alab + 1u]; ++k) {
            const std::size_t i = static_cast<std::size_t>(cur_of_asm[static_cast<std::size_t>(adj_ind[static_cast<std::size_t>(k)])]);
            if (i >= cutoff) {
                VCP_LDL_TOUCH(w, m, s, i);
                w[i] = adj_val[static_cast<std::size_t>(k)];
            }
        }
        // committed updates, ascending pivot order (rowlist is append-ordered)
        const std::vector<std::pair<Index, Index> >& rl = rowlist[c_target];
        for (std::size_t q = 0; q < rl.size(); ++q) {
            const block_t& B = blocks[static_cast<std::size_t>(rl[q].first)];
            const std::size_t pos = static_cast<std::size_t>(rl[q].second);
            if (B.width == 1) {
                const T lr = B.l1[pos];
                const T ar = B.a1[pos];
                for (std::size_t t = 0; t < B.rows.size(); ++t) {
                    const std::size_t i = static_cast<std::size_t>(B.rows[t]);
                    if (i < cutoff) continue;
                    VCP_LDL_TOUCH(w, m, s, i);
                    if (i >= c_target) w[i] -= B.l1[t] * ar;   // dense: W(i,c) -= L(i,j)*abar_cj
                    else               w[i] -= lr * B.a1[t];   // dense: W(c,i) -= L(c,j)*abar_ij
                }
            } else {
                const T lr1 = B.l1[pos], lr2 = B.l2[pos];
                const T ar1 = B.a1[pos], ar2 = B.a2[pos];
                for (std::size_t t = 0; t < B.rows.size(); ++t) {
                    const std::size_t i = static_cast<std::size_t>(B.rows[t]);
                    if (i < cutoff) continue;
                    VCP_LDL_TOUCH(w, m, s, i);
                    if (i >= c_target) w[i] -= B.l1[t] * ar1 + B.l2[t] * ar2;
                    else               w[i] -= lr1 * B.a1[t] + lr2 * B.a2[t];
                }
            }
        }
        std::sort(s.begin(), s.end());   // integer keys (P6-safe)
        (void)UN;
    };

    // symmetric exchange of CURRENT labels p <-> q: perm, label maps and the
    // committed-L row relabel (row lists only; no full-row copy).
    auto exchange_labels = [&](const std::size_t p, const std::size_t q) {
        std::swap(res.perm[p], res.perm[q]);
        const std::size_t ap = static_cast<std::size_t>(asm_of_cur[p]);
        const std::size_t aq = static_cast<std::size_t>(asm_of_cur[q]);
        std::swap(cur_of_asm[ap], cur_of_asm[aq]);
        std::swap(asm_of_cur[p], asm_of_cur[q]);
        for (std::size_t t = 0; t < rowlist[p].size(); ++t) {
            blocks[static_cast<std::size_t>(rowlist[p][t].first)]
                .rows[static_cast<std::size_t>(rowlist[p][t].second)] = static_cast<Index>(q);
        }
        for (std::size_t t = 0; t < rowlist[q].size(); ++t) {
            blocks[static_cast<std::size_t>(rowlist[q][t].first)]
                .rows[static_cast<std::size_t>(rowlist[q][t].second)] = static_cast<Index>(p);
        }
        std::swap(rowlist[p], rowlist[q]);
    };

    // apply a label exchange to an already-built working column
    auto exchange_w = [&](std::vector<T>& w, std::vector<char>& m, std::vector<Index>& s,
                          const std::size_t p, const std::size_t q) {
        const T vp = m[p] ? w[p] : T(0);
        const T vq = m[q] ? w[q] : T(0);
        VCP_LDL_TOUCH(w, m, s, p);
        VCP_LDL_TOUCH(w, m, s, q);
        w[p] = vq;
        w[q] = vp;
        std::sort(s.begin(), s.end());
    };

    // ------ main pivot loop ------
    std::size_t k = 0;
    while (k < un) {
        build_column(k, k, wa, ma, sa);

        // lambda = max_{i>k} |w_k(i)|: certified explicit ascending loop.
        R lam = R(0);
        std::size_t r = 0;
        bool r_valid = false;
        bool search_inconclusive = false;
        for (std::size_t t = 0; t < sa.size(); ++t) {
            const std::size_t i = static_cast<std::size_t>(sa[t]);
            if (i <= k) continue;
            const R a = abs(wa[i]);
            if (a > lam) { lam = a; r = i; r_valid = true; }
            else if (a <= lam) { /* certified: no update */ }
            else { search_inconclusive = true; break; }
        }
        if (search_inconclusive) {
            res.inconclusive_at = static_cast<Index>(k);
            stopped_inconclusive = true;
            break;
        }

        pivot_decision dec = pivot_decision::inconclusive;
        bool pivot_1x1 = false;
        bool swap_k_r = false;

        const R absakk = abs(ma[k] ? wa[k] : T(0));
        if (!r_valid) {
            dec = pivot_decision::acceptable;
            pivot_1x1 = true;
        } else {
            const R alam = alpha * lam;
            if (absakk >= alam) {                          // rule 1
                dec = pivot_decision::acceptable;
                pivot_1x1 = true;
            } else if (absakk < alam) {
                // sigma over current column r, rows [k, n), i != r
                build_column(r, k, wb, mb, sb);
                R sigma = R(0);
                bool sigma_inconclusive = false;
                for (std::size_t t = 0; t < sb.size(); ++t) {
                    const std::size_t i = static_cast<std::size_t>(sb[t]);
                    if (i < k || i == r) continue;
                    const R a = abs(wb[i]);
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
                        const R absarr = abs(mb[r] ? wb[r] : T(0));
                        const R asig = alpha * sigma;
                        if (absarr >= asig) {              // rule 3
                            dec = pivot_decision::acceptable;
                            pivot_1x1 = true;
                            swap_k_r = true;
                        } else if (absarr < asig) {        // rule 4
                            dec = pivot_decision::acceptable;
                            pivot_1x1 = false;
                        }
                    }
                }
            }
        }

        if (dec == pivot_decision::inconclusive) {
            res.inconclusive_at = static_cast<Index>(k);
            stopped_inconclusive = true;
            break;
        }

        if (pivot_1x1) {
            // pivot column values: w_k, or (rule 3) exchanged w_r
            std::vector<T>* w = &wa; std::vector<char>* m = &ma; std::vector<Index>* s = &sa;
            if (swap_k_r && r != k) {
                exchange_labels(k, r);
                exchange_w(wb, mb, sb, k, r);
                w = &wb; m = &mb; s = &sb;
            }
            const T d = ((*m)[k] ? (*w)[k] : T(0));
            const R absd = abs(d);
            if (absd > opt.zero_pivot_tol) {
                res.D_diag[k] = d;
                blocks.push_back(block_t());
                block_t& B = blocks.back();
                const Index bid = static_cast<Index>(blocks.size() - 1u);
                B.jcol = static_cast<Index>(k);
                B.width = 1;
                for (std::size_t t = 0; t < s->size(); ++t) {
                    const std::size_t i = static_cast<std::size_t>((*s)[t]);
                    if (i <= k) continue;
                    const T av = (*w)[i];
                    B.rows.push_back(static_cast<Index>(i));
                    B.a1.push_back(av);
                    B.l1.push_back(av / d);
                    rowlist[i].push_back(std::make_pair(bid, static_cast<Index>(B.rows.size() - 1u)));
                }
                col_block_of[k] = bid;
                ++res.n_pivots_1x1;
            } else if (absd <= opt.zero_pivot_tol) {
                // certified zero pivot: structural zero, skip, continue
                if (!any_zero_pivot) {
                    res.first_zero_pivot = static_cast<Index>(k);
                    any_zero_pivot = true;
                }
            } else {
                res.inconclusive_at = static_cast<Index>(k);
                stopped_inconclusive = true;
                break;
            }
            ++k;
        } else {
            // rule 4: 2x2 pivot after exchange (k+1 <-> r)
            if (r != k + 1) {
                exchange_labels(k + 1, r);
                exchange_w(wa, ma, sa, k + 1, r);
                exchange_w(wb, mb, sb, k + 1, r);
            }
            const T d1 = (ma[k]      ? wa[k]      : T(0));
            const T e  = (ma[k + 1u] ? wa[k + 1u] : T(0));
            const T d2 = (mb[k + 1u] ? wb[k + 1u] : T(0));
            const T det = d1 * d2 - e * e;
            const R absdet = abs(det);
            if (absdet > opt.zero_pivot_tol) {             // D-1: tol2 = zero_pivot_tol
                res.D_diag[k] = d1;
                res.D_diag[k + 1u] = d2;
                res.D_sub[k] = e;
                res.D_block2[k] = char(1);
                blocks.push_back(block_t());
                block_t& B = blocks.back();
                const Index bid = static_cast<Index>(blocks.size() - 1u);
                B.jcol = static_cast<Index>(k);
                B.width = 2;
                // union of the two sorted patterns, rows >= k+2
                std::size_t ta = 0, tb = 0;
                while (ta < sa.size() || tb < sb.size()) {
                    const std::size_t ia = (ta < sa.size()) ? static_cast<std::size_t>(sa[ta]) : UN;
                    const std::size_t ib = (tb < sb.size()) ? static_cast<std::size_t>(sb[tb]) : UN;
                    std::size_t i;
                    T w1, w2;
                    if (ia == ib)      { i = ia; w1 = wa[ia]; w2 = wb[ib]; ++ta; ++tb; }
                    else if (ia < ib)  { i = ia; w1 = wa[ia]; w2 = T(0);  ++ta; }
                    else               { i = ib; w1 = T(0);  w2 = wb[ib]; ++tb; }
                    if (i < k + 2u) continue;
                    const T l1 = (d2 * w1 - e * w2) / det;
                    const T l2 = (d1 * w2 - e * w1) / det;
                    B.rows.push_back(static_cast<Index>(i));
                    B.a1.push_back(w1);
                    B.a2.push_back(w2);
                    B.l1.push_back(l1);
                    B.l2.push_back(l2);
                    rowlist[i].push_back(std::make_pair(bid, static_cast<Index>(B.rows.size() - 1u)));
                }
                col_block_of[k] = bid;
                col_block_of[k + 1u] = bid;
                ++res.n_pivots_2x2;
            } else if (absdet <= opt.zero_pivot_tol) {
                // certified zero 2x2 block: skip both columns, continue
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
        res.status = sparse_ldl_status::inconclusive_pivot_test;
        return;
    }

    // ------ L extraction (CSC, explicit unit diagonal; same drop rule as the
    // dense kernel: a value is dropped only when CERTIFIED equal to T(0)).
    // Committed block rows carry final current labels; each pattern was
    // sorted at commit but relabeling exchanges may have perturbed the
    // order, so sort by row per column (integer keys).
    res.L_col_ptr.assign(un + 1u, Index(0));
    std::vector<std::pair<Index, std::size_t> > order;   // (row, position in block)
    for (std::size_t j = 0; j < un; ++j) {
        res.L_row_ind.push_back(static_cast<Index>(j));
        res.L_val.push_back(T(1));
        const Index bid = col_block_of[j];
        if (bid >= Index(0)) {
            const block_t& B = blocks[static_cast<std::size_t>(bid)];
            const bool second = (B.width == 2 && static_cast<std::size_t>(B.jcol) + 1u == j);
            const std::vector<T>& lv = second ? B.l2 : B.l1;
            order.clear();
            for (std::size_t t = 0; t < B.rows.size(); ++t) {
                order.push_back(std::make_pair(B.rows[t], t));
            }
            std::sort(order.begin(), order.end());       // integer keys (P6-safe)
            for (std::size_t t = 0; t < order.size(); ++t) {
                const T& x = lv[order[t].second];
                if (x == T(0)) { /* certified zero: not stored */ }
                else {
                    res.L_row_ind.push_back(order[t].first);
                    res.L_val.push_back(x);
                }
            }
        }
        res.L_col_ptr[j + 1u] = static_cast<Index>(res.L_row_ind.size());
    }
    res.nnz_L = static_cast<Index>(res.L_row_ind.size());

    res.status = any_zero_pivot ? sparse_ldl_status::zero_pivot
                                : sparse_ldl_status::success;

#undef VCP_LDL_RESET
#undef VCP_LDL_TOUCH
}

#endif // VCP_TSPARSE_SPARSE_LDL_SPARSE_IMPL_HPP
