// vcp/bfem/rt/rational_la.hpp
// RT Layer 0 detail: exact rational dense matrices (rmat) and exact Gauss
// elimination -- solve_exact / rank_exact / solve_consistent.
//
// Conforms to: RT-L0 external design v0.3 (section 2 of the internal design)
//              RT-L0 internal design v0.3 (section 2, B-3).
//
// L0 files are frozen; this is a NEW detail module under vcp/bfem/rt/ that
// only consumes the frozen detail::rational type. No floating point appears
// here. rmat is generation-time only (off the hot path; the zero-allocation
// contract does not apply, internal design section 2).
//
// Pivot rule (settled): among the nonzero pivot candidates pick the row whose
// pivot entry has the smallest total bigint limb count (numerator +
// denominator). Exactness holds for ANY nonzero pivot; the choice only
// dampens intermediate bigint growth (internal design section 2). Bareiss is
// the recorded fallback if growth is ever measured to be a problem (10.2).

#ifndef VCP_BFEM_RT_RATIONAL_LA_HPP
#define VCP_BFEM_RT_RATIONAL_LA_HPP

#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/rational.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// rmat: row-major dense matrix of detail::rational.
// ---------------------------------------------------------------------------
struct rmat {
    int rows, cols;
    std::vector<rational> a;

    rmat() : rows(0), cols(0), a() {}
    rmat(int r, int c)
        : rows(r), cols(c),
          a(static_cast<std::size_t>(r) * static_cast<std::size_t>(c)) {
        if (r < 0 || c < 0)
            throw std::invalid_argument("bfem::rmat: negative size");
    }

    rational& at(int i, int j) {
        assert(i >= 0 && i < rows && j >= 0 && j < cols);
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)
                 + static_cast<std::size_t>(j)];
    }
    const rational& at(int i, int j) const {
        assert(i >= 0 && i < rows && j >= 0 && j < cols);
        return a[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)
                 + static_cast<std::size_t>(j)];
    }

    static rmat identity(int n) {
        rmat r(n, n);
        for (int i = 0; i < n; ++i) r.at(i, i) = rational(1);
        return r;
    }
};

inline rmat mul(const rmat& A, const rmat& B) {
    if (A.cols != B.rows)
        throw std::invalid_argument("bfem::rmat::mul: size mismatch");
    rmat C(A.rows, B.cols);
    for (int i = 0; i < A.rows; ++i) {
        for (int k = 0; k < A.cols; ++k) {
            const rational& aik = A.at(i, k);
            if (aik.is_zero()) continue;
            for (int j = 0; j < B.cols; ++j) {
                const rational& bkj = B.at(k, j);
                if (bkj.is_zero()) continue;
                C.at(i, j) += aik * bkj;
            }
        }
    }
    return C;
}

inline rmat transpose(const rmat& A) {
    rmat T(A.cols, A.rows);
    for (int i = 0; i < A.rows; ++i)
        for (int j = 0; j < A.cols; ++j)
            T.at(j, i) = A.at(i, j);
    return T;
}

// ---------------------------------------------------------------------------
// elimination statistics (RTX-2 / internal design 6: bigint growth watch).
// Passed as an optional out-parameter; a null pointer disables recording.
// ---------------------------------------------------------------------------
struct la_stats {
    int max_limbs;                    // max (num + den) limbs of any entry seen
    unsigned long long ops;           // rational multiply/divide count
    la_stats() : max_limbs(0), ops(0) {}
};

namespace la_impl {

inline int limbs_of(const rational& q) {
    return q.num().num_limbs() + q.den().num_limbs();
}

inline void record(la_stats* st, const rational& q) {
    if (!st) return;
    int l = limbs_of(q);
    if (l > st->max_limbs) st->max_limbs = l;
    ++st->ops;
}

// pick the pivot row in rows [from, rows) of column col: nonzero entry with
// the smallest total limb count; -1 when the whole column tail is zero.
inline int pick_pivot(const rmat& A, int col, int from) {
    int best = -1;
    int best_limbs = 0;
    for (int r = from; r < A.rows; ++r) {
        const rational& v = A.at(r, col);
        if (v.is_zero()) continue;
        int l = limbs_of(v);
        if (best < 0 || l < best_limbs) {
            best = r;
            best_limbs = l;
        }
    }
    return best;
}

inline void swap_rows(rmat& A, int r1, int r2) {
    if (r1 == r2) return;
    for (int j = 0; j < A.cols; ++j) {
        rational t = A.at(r1, j);
        A.at(r1, j) = A.at(r2, j);
        A.at(r2, j) = t;
    }
}

// forward elimination of the augmented pair (A | B), pivoting over the first
// min(rows, cols_of_A) columns. Returns the pivot count (= exact rank of A).
// Row swaps and eliminations are applied to B in lockstep.
inline int forward_eliminate(rmat& A, rmat* B, la_stats* st) {
    int piv = 0;
    for (int col = 0; col < A.cols && piv < A.rows; ++col) {
        int p = pick_pivot(A, col, piv);
        if (p < 0) continue;                     // zero column tail: no pivot here
        swap_rows(A, piv, p);
        if (B) swap_rows(*B, piv, p);
        const rational pivval = A.at(piv, col);  // copy: row content mutates below
        for (int r = piv + 1; r < A.rows; ++r) {
            if (A.at(r, col).is_zero()) continue;
            rational f = A.at(r, col) / pivval;
            record(st, f);
            A.at(r, col) = rational(0);
            for (int j = col + 1; j < A.cols; ++j) {
                if (A.at(piv, j).is_zero()) continue;
                A.at(r, j) -= f * A.at(piv, j);
                record(st, A.at(r, j));
            }
            if (B) {
                for (int j = 0; j < B->cols; ++j) {
                    if (B->at(piv, j).is_zero()) continue;
                    B->at(r, j) -= f * B->at(piv, j);
                    record(st, B->at(r, j));
                }
            }
        }
        ++piv;
    }
    return piv;
}

} // namespace la_impl

// ---------------------------------------------------------------------------
// solve_exact: A x = B for square regular A (multiple right hand sides; the
// inverse is solve_exact(A, identity)). Singular input raises logic_error --
// the RT Vandermonde is regular by theory, so singularity is an
// implementation bug (internal design section 2).
// ---------------------------------------------------------------------------
inline rmat solve_exact(rmat A, rmat B, la_stats* st = 0) {
    if (A.rows != A.cols)
        throw std::invalid_argument("bfem::solve_exact: A not square");
    if (B.rows != A.rows)
        throw std::invalid_argument("bfem::solve_exact: B row mismatch");
    const int n = A.rows;
    int piv = la_impl::forward_eliminate(A, &B, st);
    if (piv < n)
        throw std::logic_error("bfem::solve_exact: singular matrix");
    // back substitution (A is now upper triangular with nonzero diagonal)
    rmat X(n, B.cols);
    for (int i = n - 1; i >= 0; --i) {
        for (int j = 0; j < B.cols; ++j) {
            rational acc = B.at(i, j);
            for (int k = i + 1; k < n; ++k) {
                if (A.at(i, k).is_zero() || X.at(k, j).is_zero()) continue;
                acc -= A.at(i, k) * X.at(k, j);
            }
            if (!acc.is_zero()) acc /= A.at(i, i);
            la_impl::record(st, acc);
            X.at(i, j) = acc;
        }
    }
    return X;
}

// ---------------------------------------------------------------------------
// rank_exact (v0.3, B-3): exact rank = pivot count of the elimination.
// ---------------------------------------------------------------------------
inline int rank_exact(rmat A, la_stats* st = 0) {
    return la_impl::forward_eliminate(A, 0, st);
}

// ---------------------------------------------------------------------------
// solve_consistent (v0.3, B-3): A x = B for a tall consistent system
// (A: r x c, r >= c, full column rank). After elimination every residual row
// must be EXACTLY zero on the B side; otherwise logic_error (this implements
// the generation-time check of internal design 5.2: the RT edge flux MUST
// reduce to degree k, a nonzero residual is an implementation bug).
// Rank-deficient A also raises logic_error.
// ---------------------------------------------------------------------------
inline rmat solve_consistent(rmat A, rmat B, la_stats* st = 0) {
    if (A.rows < A.cols)
        throw std::invalid_argument("bfem::solve_consistent: A has rows < cols");
    if (B.rows != A.rows)
        throw std::invalid_argument("bfem::solve_consistent: B row mismatch");
    const int c = A.cols;
    int piv = la_impl::forward_eliminate(A, &B, st);
    if (piv < c)
        throw std::logic_error("bfem::solve_consistent: rank-deficient matrix");
    // consistency: rows c..r-1 of the eliminated B must be exactly zero
    for (int r = c; r < A.rows; ++r)
        for (int j = 0; j < B.cols; ++j)
            if (!B.at(r, j).is_zero())
                throw std::logic_error(
                    "bfem::solve_consistent: inconsistent right hand side "
                    "(nonzero elimination residual)");
    // back substitution on the top c x c triangle
    rmat X(c, B.cols);
    for (int i = c - 1; i >= 0; --i) {
        for (int j = 0; j < B.cols; ++j) {
            rational acc = B.at(i, j);
            for (int k = i + 1; k < c; ++k) {
                if (A.at(i, k).is_zero() || X.at(k, j).is_zero()) continue;
                acc -= A.at(i, k) * X.at(k, j);
            }
            if (!acc.is_zero()) acc /= A.at(i, i);
            la_impl::record(st, acc);
            X.at(i, j) = acc;
        }
    }
    return X;
}

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RATIONAL_LA_HPP
