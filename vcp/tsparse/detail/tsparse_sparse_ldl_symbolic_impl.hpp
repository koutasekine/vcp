// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLDL-SP SP-0 -- symbolic phase of the sparse LDL^T factorization.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER the
// sparse_ldl_* type skeleton (enums / options / result) and the SLU header
// (pattern-only ordering functions) are in scope.  It has no
// "namespace vcp { }" wrapper; it is injected by tsparse_sparse_ldl.hpp.
//
// Do NOT include this file directly.  Include:
//   <vcp/tsparse/tsparse_sparse_ldl.hpp>
//
// Role (SLDL-SP design v1 SS3.1): the symbolic phase is T-INDEPENDENT and
// self-contained (B2 boundary requirement 1).  From the lower triangle of the
// input pattern it computes, in the ordering-permuted labels:
//   - the pre-permutation perm0 (new->old) and its inverse,
//   - the symmetric lower pattern (diagonal always present),
//   - the elimination tree (Liu, with path compression),
//   - a postorder of that tree,
//   - the column counts |L(:,j)| of the STATIC (exchange-free) factor,
//   - the fundamental supernode partition,
//   - the below-supernode row pattern of every supernode,
//   - the workspace sizes the numeric phase needs.
// The result object is READ-ONLY after construction (no mutable members, no
// in-place refinement by the numeric phase).
//
// Column patterns are NOT stored one by one: inside supernode s = [c0, c1)
// every column j has the exact static pattern
//     L(:,j) = {j, j+1, ..., c1-1} U R_s
// (R_s = the stored below-supernode row list), which is the defining property
// of a fundamental supernode.  Storing R_s per supernode instead of the full
// L pattern is therefore lossless and costs supernodal, not nnz_L, memory.
//
// Analysis level: baseline_dynamic performs dynamic symmetric exchanges, so
// the static structure is meaningless for it and computing it would cost
// O(nnz_L) time and memory for nothing.  sparse_ldl_symbolic_level::
// ordering_only therefore stops after perm0; the frozen baseline path uses
// exactly that level and is unaffected by everything below it.

#ifndef VCP_TSPARSE_SPARSE_LDL_SYMBOLIC_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LDL_SYMBOLIC_IMPL_HPP

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>

// ---------------------------------------------------------------------------
// Analysis level (see the header comment).
// ---------------------------------------------------------------------------
enum class sparse_ldl_symbolic_level {
    ordering_only,   // perm0 / pinv0 / ordering_used only
    full             // + pattern, etree, postorder, col_count, supernodes
};

enum class sparse_ldl_symbolic_status {
    success,
    invalid_input,     // n < 0 or malformed CSC
    invalid_options,   // ordering value outside the enum
    internal_error
};

struct sparse_ldl_symbolic_options {
    sparse_ldl_ordering      ordering;
    sparse_ldl_symbolic_level level;

    sparse_ldl_symbolic_options()
        : ordering(sparse_ldl_ordering::auto_select),
          level(sparse_ldl_symbolic_level::full) {}
};

// ---------------------------------------------------------------------------
// Symbolic result.  Integer-only (P4); every field is filled by
// sparse_ldl_symbolic_analyze and never modified afterwards.
// ---------------------------------------------------------------------------
template <class Index>
struct sparse_ldl_symbolic_result {
    static_assert(std::is_signed<Index>::value,
                  "sparse LDL Index must be signed");

    sparse_ldl_symbolic_status status;
    sparse_ldl_symbolic_level  level;
    sparse_ldl_ordering        ordering_used;   // auto_select resolution
    Index n;

    // ---- always valid on success
    std::vector<Index> perm0;    // new -> old (the ordering pre-permutation)
    std::vector<Index> pinv0;    // old -> new

    // ---- valid only when level == full
    // Permuted symmetric lower pattern (rows ascending per column, diagonal
    // always present even when structurally absent from the input).
    std::vector<Index> A_col_ptr, A_row_ind;
    std::vector<Index> parent;         // etree; -1 = root, else j < parent[j] < n
    std::vector<Index> postorder;      // postorder[k] = column visited k-th
    std::vector<Index> col_count;      // |L(:,j)| of the static factor, diagonal included

    // Fundamental supernode partition: supernode s owns the columns
    // [supernode_ptr[s], supernode_ptr[s+1]).
    std::vector<Index> supernode_ptr;        // size n_supernodes + 1
    std::vector<Index> column_to_supernode;  // size n

    // Below-supernode row pattern: supernode s occupies
    // sn_row_ind[sn_row_ptr[s] .. sn_row_ptr[s+1]), sorted ascending, every
    // entry >= supernode_ptr[s+1].
    std::vector<Index> sn_row_ptr, sn_row_ind;

    Index nnz_L_static;          // entries of the static L (unit diagonal included)
    Index n_supernodes;
    Index max_supernode_width;
    Index max_panel_rows;        // max over s of (width + |R_s|)
    Index max_update_rows;       // max over s of |R_s| (update source height)

    sparse_ldl_symbolic_result()
        : status(sparse_ldl_symbolic_status::internal_error),
          level(sparse_ldl_symbolic_level::ordering_only),
          ordering_used(sparse_ldl_ordering::auto_select),
          n(Index(0)),
          nnz_L_static(Index(0)), n_supernodes(Index(0)),
          max_supernode_width(Index(0)), max_panel_rows(Index(0)),
          max_update_rows(Index(0)) {}
};

namespace sparse_ldl_detail {

// Pattern-only CSC validation (values are irrelevant here).  Mirrors the rules
// of sparse_ldl_validate_csc_ minus the value-array length checks.
template <class Index>
bool sparse_ldl_validate_pattern_(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind)
{
    if (n < Index(0)) return false;
    const std::size_t un = static_cast<std::size_t>(n);
    if (col_ptr.size() != un + 1u) return false;
    if (col_ptr[0] != Index(0)) return false;
    for (std::size_t c = 0; c < un; ++c) {
        if (col_ptr[c + 1u] < col_ptr[c]) return false;
    }
    if (row_ind.size() != static_cast<std::size_t>(col_ptr[un])) return false;
    for (std::size_t c = 0; c < un; ++c) {
        Index prev = Index(-1);
        for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
            const Index r = row_ind[static_cast<std::size_t>(k)];
            if (r < Index(0) || r >= n) return false;
            if (r <= prev) return false;
            prev = r;
        }
    }
    return true;
}

// auto_select resolution for the ordering, in ONE place (design v1 SS4).
// Contract (ldl_design_v2 SS0.1 decision 1): the library chooses, and the
// choice is always reported in ordering_used.
inline sparse_ldl_ordering resolve_auto_ordering(sparse_ldl_ordering o) {
    return (o == sparse_ldl_ordering::auto_select) ? sparse_ldl_ordering::amd : o;
}

// auto_select resolution for the method, in ONE place (design v1 SS4).
// v1 keeps auto_select -> baseline_dynamic; changing that resolution is the
// separate SP-D1 track and must not happen here as a side effect.
inline sparse_ldl_method resolve_auto_method(sparse_ldl_method m) {
    return (m == sparse_ldl_method::auto_select) ? sparse_ldl_method::baseline_dynamic : m;
}

// auto_select resolution for the diagonal-block kernel (design v1 D-3).
// auto resolves to gemm: the final default fixed by SP-3
// from the four-machine calibration.
inline sparse_ldl_diag_kernel resolve_auto_diag_kernel(sparse_ldl_diag_kernel k) {
    return (k == sparse_ldl_diag_kernel::auto_select) ? sparse_ldl_diag_kernel::gemm : k;
}

// ---------------------------------------------------------------------------
// Growth diagnostic (design v1 D-5, P4-compliant integer representation).
//
// growth_log2 = floor(log2( max|d| / max|a_ii| )) computed WITHOUT a logarithm
// (the scalar contract has no log requirement for this module and a log of a
// magnitude would be one more requirement on T): the ratio is normalised into
// [1,2) by exact halving/doubling, and the exponent is counted.  Both loops
// use certified comparisons, so for a scalar whose comparisons cannot decide,
// the loop stops early and the reported exponent is a certified LOWER bound
// -- never an overstatement of the growth.
//
// Validity is the bool flag alone (P4): it stays false when the reference
// max|a_ii| or the observed max|d| cannot be certified positive.
// ---------------------------------------------------------------------------
template <class R, class T, class Index>
inline void sparse_ldl_set_growth_(
    const R& max_d, const R& max_a, const bool max_a_valid,
    sparse_ldl_result<T, Index>& res)
{
    res.growth_log2 = 0;
    res.growth_valid = false;
    if (!max_a_valid) return;
    if (!(max_a > R(0))) return;    // no certified positive reference
    if (!(max_d > R(0))) return;    // degenerate factorization (all pivots zero)

    const R two = R(2);
    R r = max_d / max_a;
    int e = 0;
    while (e < 1024) {
        if (r >= two) { r = r / two; ++e; }   // exact in any binary format
        else break;                            // certified < 2, or undecidable
    }
    while (e > -1024) {
        if (r < R(1)) { r = r * two; --e; }
        else break;                            // certified >= 1, or undecidable
    }
    res.growth_log2 = e;
    res.growth_valid = true;
}

} // namespace sparse_ldl_detail

// ---------------------------------------------------------------------------
// sparse_ldl_symbolic_analyze -- non-throwing symbolic phase.
//
//   col_ptr/row_ind : CSC pattern of the input; only the LOWER triangle
//                     (i >= j) is read, exactly like the numeric kernels.
// ---------------------------------------------------------------------------
template <class Index>
sparse_ldl_symbolic_result<Index>
sparse_ldl_symbolic_analyze(
    const Index n,
    const std::vector<Index>& col_ptr,
    const std::vector<Index>& row_ind,
    const sparse_ldl_symbolic_options& sopt)
{
    static_assert(std::is_signed<Index>::value,
                  "sparse LDL Index must be signed");

    sparse_ldl_symbolic_result<Index> sym;
    sym.level = sopt.level;
    sym.n = n;

    if (!sparse_ldl_detail::sparse_ldl_validate_pattern_(n, col_ptr, row_ind)) {
        sym.status = sparse_ldl_symbolic_status::invalid_input;
        return sym;
    }
    switch (sopt.ordering) {
    case sparse_ldl_ordering::auto_select:
    case sparse_ldl_ordering::natural:
    case sparse_ldl_ordering::rcm:
    case sparse_ldl_ordering::amd:
    case sparse_ldl_ordering::nested_dissection:
        break;
    default:
        sym.status = sparse_ldl_symbolic_status::invalid_options;
        return sym;
    }

    const std::size_t un = static_cast<std::size_t>(n);
    sym.ordering_used = sparse_ldl_detail::resolve_auto_ordering(sopt.ordering);

    // ---- 1. ordering (pattern-only, integer-only; the SLU ordering functions
    // are reused by include, unmodified: the graph builder symmetrizes every
    // off-diagonal edge, so the lower-triangle CSC pattern is a valid input).
    sym.perm0.resize(un);
    for (std::size_t i = 0; i < un; ++i) sym.perm0[i] = static_cast<Index>(i);
    switch (sym.ordering_used) {
    case sparse_ldl_ordering::natural:
        break;
    case sparse_ldl_ordering::rcm:
        sym.perm0 = sparse_lu_rcm_ordering(n, col_ptr, row_ind);
        break;
    case sparse_ldl_ordering::amd:
        sym.perm0 = sparse_lu_amd_ordering(n, col_ptr, row_ind);
        break;
    case sparse_ldl_ordering::nested_dissection:
        sym.perm0 = sparse_lu_nested_dissection_ordering(n, col_ptr, row_ind);
        break;
    default:
        sym.status = sparse_ldl_symbolic_status::internal_error;
        return sym;
    }
    if (sym.perm0.size() != un) {
        sym.status = sparse_ldl_symbolic_status::internal_error;
        return sym;
    }
    sym.pinv0.assign(un, Index(-1));
    for (std::size_t i = 0; i < un; ++i) {
        const Index o = sym.perm0[i];
        if (o < Index(0) || o >= n || sym.pinv0[static_cast<std::size_t>(o)] != Index(-1)) {
            sym.status = sparse_ldl_symbolic_status::internal_error;   // not a permutation
            return sym;
        }
        sym.pinv0[static_cast<std::size_t>(o)] = static_cast<Index>(i);
    }

    if (sopt.level == sparse_ldl_symbolic_level::ordering_only) {
        sym.status = sparse_ldl_symbolic_status::success;
        return sym;
    }

    // ---- 2. permuted symmetric lower pattern (diagonal always present).
    // A stored lower entry (r,c) becomes (max(a,b), min(a,b)) with a=pinv0[r],
    // b=pinv0[c]; strictly upper stored entries are ignored, exactly as the
    // numeric kernels ignore them (design v2 SS1.2).
    {
        std::vector<Index> deg(un, Index(1));   // the diagonal of every column
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r < static_cast<Index>(c)) continue;             // strictly upper: ignored
                const Index a = sym.pinv0[static_cast<std::size_t>(r)];
                const Index b = sym.pinv0[c];
                if (a == b) continue;                                 // diagonal already counted
                ++deg[static_cast<std::size_t>(a < b ? a : b)];
            }
        }
        sym.A_col_ptr.assign(un + 1u, Index(0));
        for (std::size_t j = 0; j < un; ++j) {
            sym.A_col_ptr[j + 1u] = sym.A_col_ptr[j] + deg[j];
        }
        sym.A_row_ind.assign(static_cast<std::size_t>(sym.A_col_ptr[un]), Index(0));
        std::vector<Index> head(sym.A_col_ptr.begin(), sym.A_col_ptr.end() - 1);
        for (std::size_t j = 0; j < un; ++j) {
            sym.A_row_ind[static_cast<std::size_t>(head[j]++)] = static_cast<Index>(j);
        }
        for (std::size_t c = 0; c < un; ++c) {
            for (Index k = col_ptr[c]; k < col_ptr[c + 1u]; ++k) {
                const Index r = row_ind[static_cast<std::size_t>(k)];
                if (r < static_cast<Index>(c)) continue;
                const Index a = sym.pinv0[static_cast<std::size_t>(r)];
                const Index b = sym.pinv0[c];
                if (a == b) continue;
                const Index lo = (a < b) ? a : b;
                const Index hi = (a < b) ? b : a;
                sym.A_row_ind[static_cast<std::size_t>(head[static_cast<std::size_t>(lo)]++)] = hi;
            }
        }
        for (std::size_t j = 0; j < un; ++j) {
            std::sort(sym.A_row_ind.begin() + static_cast<std::ptrdiff_t>(sym.A_col_ptr[j]),
                      sym.A_row_ind.begin() + static_cast<std::ptrdiff_t>(sym.A_col_ptr[j + 1u]));
        }
        // Duplicate lower entries cannot occur: the input CSC is validated to
        // have strictly ascending rows per column, and (r,c) -> (max,min) is
        // injective on the lower triangle.
    }

    // ---- 3. elimination tree (Liu's algorithm with path compression).
    // Needs, for every column k, the entries A(k, i) with i < k -- i.e. the
    // rows k of the earlier columns i.  rowlist[k] collects them in ascending
    // i because the columns are scanned in ascending order.
    std::vector<std::vector<Index> > rowlist(un);
    for (std::size_t i = 0; i < un; ++i) {
        for (Index p = sym.A_col_ptr[i]; p < sym.A_col_ptr[i + 1u]; ++p) {
            const Index r = sym.A_row_ind[static_cast<std::size_t>(p)];
            if (r > static_cast<Index>(i)) {
                rowlist[static_cast<std::size_t>(r)].push_back(static_cast<Index>(i));
            }
        }
    }
    sym.parent.assign(un, Index(-1));
    {
        std::vector<Index> ancestor(un, Index(-1));
        for (std::size_t k = 0; k < un; ++k) {
            const std::vector<Index>& rl = rowlist[k];
            for (std::size_t t = 0; t < rl.size(); ++t) {
                Index i = rl[t];
                while (i != Index(-1) && i < static_cast<Index>(k)) {
                    const Index inext = ancestor[static_cast<std::size_t>(i)];
                    ancestor[static_cast<std::size_t>(i)] = static_cast<Index>(k);
                    if (inext == Index(-1)) {
                        sym.parent[static_cast<std::size_t>(i)] = static_cast<Index>(k);
                    }
                    i = inext;
                }
            }
        }
    }

    // ---- 4. postorder of the elimination forest (explicit stack, no
    // recursion).  Children are visited in ascending column order.
    sym.postorder.clear();
    sym.postorder.reserve(un);
    {
        // child lists (counting sort by parent, ascending child order)
        std::vector<Index> child_ptr(un + 2u, Index(0));
        for (std::size_t j = 0; j < un; ++j) {
            const Index p = sym.parent[j];
            const std::size_t slot = (p == Index(-1)) ? un : static_cast<std::size_t>(p);
            ++child_ptr[slot + 1u];
        }
        for (std::size_t s = 0; s + 1u < child_ptr.size(); ++s) {
            child_ptr[s + 1u] = child_ptr[s + 1u] + child_ptr[s];
        }
        std::vector<Index> child_ind(un, Index(0));
        {
            std::vector<Index> head(child_ptr.begin(), child_ptr.end() - 1);
            for (std::size_t j = 0; j < un; ++j) {
                const Index p = sym.parent[j];
                const std::size_t slot = (p == Index(-1)) ? un : static_cast<std::size_t>(p);
                child_ind[static_cast<std::size_t>(head[slot]++)] = static_cast<Index>(j);
            }
        }
        // iterative DFS: next_child[v] = index of the next unvisited child
        std::vector<Index> next_child(un, Index(0));
        for (std::size_t j = 0; j < un; ++j) next_child[j] = child_ptr[j];
        std::vector<Index> stack;
        stack.reserve(un);
        for (Index rt = child_ptr[un]; rt < child_ptr[un + 1u]; ++rt) {
            stack.push_back(child_ind[static_cast<std::size_t>(rt)]);
            while (!stack.empty()) {
                const Index v = stack.back();
                const std::size_t sv = static_cast<std::size_t>(v);
                if (next_child[sv] < child_ptr[sv + 1u]) {
                    const Index c = child_ind[static_cast<std::size_t>(next_child[sv]++)];
                    stack.push_back(c);
                } else {
                    sym.postorder.push_back(v);
                    stack.pop_back();
                }
            }
        }
        if (sym.postorder.size() != un) {
            sym.status = sparse_ldl_symbolic_status::internal_error;
            return sym;
        }
    }

    // ---- 5. column counts of the STATIC factor.
    // Row-oriented symbolic factorization: for row k, every seed i < k with
    // A(k,i) != 0 walks up the etree until a column already marked for row k
    // is met; each newly marked column j gains the entry L(k,j).  This is the
    // exact static structure (no exchanges), computed in O(nnz_L) time and
    // O(n) memory -- the row indices themselves are not stored here (see the
    // header comment: they are recovered per supernode in step 7).
    sym.col_count.assign(un, Index(1));   // the unit diagonal
    {
        std::vector<Index> mark(un, Index(-1));
        for (std::size_t k = 0; k < un; ++k) {
            mark[k] = static_cast<Index>(k);
            const std::vector<Index>& rl = rowlist[k];
            for (std::size_t t = 0; t < rl.size(); ++t) {
                Index j = rl[t];
                while (j != Index(-1) && mark[static_cast<std::size_t>(j)] != static_cast<Index>(k)) {
                    mark[static_cast<std::size_t>(j)] = static_cast<Index>(k);
                    ++sym.col_count[static_cast<std::size_t>(j)];
                    j = sym.parent[static_cast<std::size_t>(j)];
                }
            }
        }
    }

    // ---- 6. fundamental supernode partition (implementation guide SS1-2):
    // column j joins j-1  <=>  parent(j-1) == j  AND  col_count(j) == col_count(j-1) - 1.
    sym.supernode_ptr.clear();
    sym.column_to_supernode.assign(un, Index(0));
    if (un > 0u) {
        sym.supernode_ptr.push_back(Index(0));
        for (std::size_t j = 1; j < un; ++j) {
            const bool same =
                (sym.parent[j - 1u] == static_cast<Index>(j)) &&
                (sym.col_count[j] == sym.col_count[j - 1u] - Index(1));
            if (!same) sym.supernode_ptr.push_back(static_cast<Index>(j));
        }
        sym.supernode_ptr.push_back(static_cast<Index>(n));
    } else {
        sym.supernode_ptr.push_back(Index(0));
    }
    sym.n_supernodes = static_cast<Index>(sym.supernode_ptr.size()) - Index(1);
    for (Index s = Index(0); s < sym.n_supernodes; ++s) {
        const std::size_t ss = static_cast<std::size_t>(s);
        for (Index j = sym.supernode_ptr[ss]; j < sym.supernode_ptr[ss + 1u]; ++j) {
            sym.column_to_supernode[static_cast<std::size_t>(j)] = s;
        }
        const Index w = sym.supernode_ptr[ss + 1u] - sym.supernode_ptr[ss];
        if (w > sym.max_supernode_width) sym.max_supernode_width = w;
    }

    // ---- 7. below-supernode row patterns, bottom-up over the supernode tree.
    //   R_s = { rows >= c1 of A(:, [c0,c1)) }  U  ( U_{child t of s} R_t \ [0,c1) )
    // Supernodes are processed in ascending index; every child supernode has a
    // smaller index than its parent (parent columns are larger), so its R is
    // already final when the parent is reached.
    {
        std::vector<std::vector<Index> > R(static_cast<std::size_t>(sym.n_supernodes));
        std::vector<Index> mark(un, Index(-1));
        std::vector<Index> buf;
        // supernode children: parent supernode of s = supernode of parent(last col)
        std::vector<Index> sn_parent(static_cast<std::size_t>(sym.n_supernodes), Index(-1));
        for (Index s = Index(0); s < sym.n_supernodes; ++s) {
            const std::size_t ss = static_cast<std::size_t>(s);
            const Index c1 = sym.supernode_ptr[ss + 1u];
            const Index plast = sym.parent[static_cast<std::size_t>(c1 - Index(1))];
            sn_parent[ss] = (plast == Index(-1))
                          ? Index(-1)
                          : sym.column_to_supernode[static_cast<std::size_t>(plast)];
        }
        std::vector<std::vector<Index> > sn_child(static_cast<std::size_t>(sym.n_supernodes));
        for (Index s = Index(0); s < sym.n_supernodes; ++s) {
            const Index p = sn_parent[static_cast<std::size_t>(s)];
            if (p != Index(-1)) sn_child[static_cast<std::size_t>(p)].push_back(s);
        }

        sym.sn_row_ptr.assign(static_cast<std::size_t>(sym.n_supernodes) + 1u, Index(0));
        sym.sn_row_ind.clear();
        for (Index s = Index(0); s < sym.n_supernodes; ++s) {
            const std::size_t ss = static_cast<std::size_t>(s);
            const Index c0 = sym.supernode_ptr[ss];
            const Index c1 = sym.supernode_ptr[ss + 1u];
            buf.clear();
            for (Index j = c0; j < c1; ++j) {
                for (Index p = sym.A_col_ptr[static_cast<std::size_t>(j)];
                     p < sym.A_col_ptr[static_cast<std::size_t>(j) + 1u]; ++p) {
                    const Index r = sym.A_row_ind[static_cast<std::size_t>(p)];
                    if (r >= c1 && mark[static_cast<std::size_t>(r)] != s) {
                        mark[static_cast<std::size_t>(r)] = s;
                        buf.push_back(r);
                    }
                }
            }
            const std::vector<Index>& ch = sn_child[ss];
            for (std::size_t t = 0; t < ch.size(); ++t) {
                const std::vector<Index>& Rt = R[static_cast<std::size_t>(ch[t])];
                for (std::size_t q = 0; q < Rt.size(); ++q) {
                    const Index r = Rt[q];
                    if (r >= c1 && mark[static_cast<std::size_t>(r)] != s) {
                        mark[static_cast<std::size_t>(r)] = s;
                        buf.push_back(r);
                    }
                }
                // a child's rows are consumed once; free them eagerly
                std::vector<Index>().swap(R[static_cast<std::size_t>(ch[t])]);
            }
            std::sort(buf.begin(), buf.end());   // integer keys (P6-safe)
            R[ss] = buf;
            // R[ss] stays alive until its parent supernode (always a LARGER
            // index) consumes it; the flat output is appended here, in
            // ascending supernode order, so the eager child release above
            // cannot lose anything.
            sym.sn_row_ind.insert(sym.sn_row_ind.end(), buf.begin(), buf.end());
            const Index w = c1 - c0;
            const Index h = static_cast<Index>(buf.size());
            if (w + h > sym.max_panel_rows) sym.max_panel_rows = w + h;
            if (h > sym.max_update_rows) sym.max_update_rows = h;
            sym.nnz_L_static = sym.nnz_L_static + (w * (w + Index(1))) / Index(2) + w * h;
            sym.sn_row_ptr[ss + 1u] = sym.sn_row_ptr[ss] + h;
        }

        // Cross-check against the independently computed column counts: for a
        // fundamental supernode, col_count(c0) = width + |R_s|.  A mismatch
        // means the partition and the counts disagree -> internal_error rather
        // than a silently wrong structure.
        for (Index s = Index(0); s < sym.n_supernodes; ++s) {
            const std::size_t ss = static_cast<std::size_t>(s);
            const Index c0 = sym.supernode_ptr[ss];
            const Index c1 = sym.supernode_ptr[ss + 1u];
            const Index h = sym.sn_row_ptr[ss + 1u] - sym.sn_row_ptr[ss];
            if (sym.col_count[static_cast<std::size_t>(c0)] != (c1 - c0) + h) {
                sym.status = sparse_ldl_symbolic_status::internal_error;
                return sym;
            }
        }
    }

    sym.status = sparse_ldl_symbolic_status::success;
    return sym;
}

#endif // VCP_TSPARSE_SPARSE_LDL_SYMBOLIC_IMPL_HPP
