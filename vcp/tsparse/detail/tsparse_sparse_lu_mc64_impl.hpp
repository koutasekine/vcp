// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// O4 -- MC64 / static pivoting (zero-free diagonal; saddle-point solvability).
//
// Implements the value-reading maximum-weight bipartite matching (MC64-style)
// plus the accompanying Dr/Dc equilibration scaling that places large entries on
// the diagonal so that the existing threshold partial-pivoting factorization can
// factor saddle-point systems (RT mixed / incompressible Stokes-like
// [ A B^T ; B 0 ]) whose zero (2,2) block defeats diagonal-only pivoting.
//
// SCOPE (O4, strict): matching + Dr/Dc scaling + the existing pivoting only.
// No near-zero pivot perturbation and no iterative refinement (those are O4.1).
// The system is NOT perturbed: A x = b is transformed to a tracked, exactly
// invertible form A -> P_static * Dr * A * Dc * Qc, and the solve un-applies
// P_static, Dr, Dc, Qc so the residual is measured on the ORIGINAL A.
//
// Invariants (roadmap §3.3 S-A..S-E):
//   S-A  Reads values; transforms the system; solve un-applies the transform.
//   S-B  Dr/Dc/P_static are explicit, tracked transformations stored in the
//        existing storage slots (not silently folded into the factor).
//   S-C  Default pivoting (threshold_partial) is untouched; static_mc64 is opt-in.
//   S-D  No perfect matching / structurally singular -> safe failure
//        (structural_singularity); numerically singular B -> numeric stage fails
//        safely.  Never a silent wrong answer.
//   S-E  Deterministic matching (fixed index tie-break; no randomness).
//
// The matching/scaling is computed in plain double scalar arithmetic and is
// T-INDEPENDENT: magnitudes are converted to double via static_cast (every
// supported real_type, incl. kv::dd / kv::mpfr, provides operator double()), and
// the resulting scale factors are converted back with static_cast<T>(...).  No
// branch depends on what T is.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER csc_storage,
// sparse_lu_scalar_policy, sparse_lu_identity_permutation and
// sparse_lu_inverse_permutation are in scope.  It has no namespace wrapper; it is
// injected by tsparse_sparse_lu.hpp.  Do NOT include this file directly.

#ifndef VCP_TSPARSE_SPARSE_LU_MC64_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_MC64_IMPL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include <vcp/error.hpp>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// sparse_lu_mc64_transform: result of the MC64 static-pivoting transform.
//
//   success    : a perfect (zero-free-diagonal) matching was found.
//   status     : success or structural_singularity (no perfect matching).
//   B          : effective matrix = P_static * Dr * A * Dc * Qc (CSC), the matrix
//                the numeric factorization runs on.  Its diagonal is zero-free.
//   p_static   : p_static[B_row] = original row brought to B_row by the matching
//                composed with the column ordering Qc.  Used to compose the final
//                row permutation with the numeric dynamic pivoting.
//   Dr, Dc     : size-n equilibration scalings in ORIGINAL row / column
//                coordinates (T-typed; solve multiplies b by Dr and x by Dc).
// ---------------------------------------------------------------------------
template <class T, class Index>
struct sparse_lu_mc64_transform {
    bool                  success;
    sparse_lu_status      status;
    csc_storage<T, Index> B;
    std::vector<Index>    p_static;
    std::vector<T>        Dr;
    std::vector<T>        Dc;

    sparse_lu_mc64_transform()
        : success(false),
          status(sparse_lu_status::structural_singularity) {}
};

// ---------------------------------------------------------------------------
// sparse_lu_mc64_assignment: minimum-cost perfect bipartite assignment.
//
// Deterministic Jonker-Volgenant / successive-shortest-augmenting-path
// (the classic Kuhn-Munkres "Hungarian" potential form), operating on a dense
// n*n double cost matrix in row-major (cost[i*n + j]).  Rows are workers, columns
// are jobs.  Minimizes sum_i cost[i][assign[i]].
//
// On return:
//   match_row_for_col[j] = row matched to column j   (perfect assignment)
//   u[i] (rows), v[j] (cols): LP dual potentials with u[i]+v[j] <= cost[i][j],
//   equality on matched pairs.  These are the MC64 scaling potentials.
//
// The graph is complete (forbidden edges carry a large finite cost `big`), so a
// perfect assignment always exists; the caller checks matched edges against `big`
// to detect structural singularity.  Deterministic: ties resolved by lowest
// index (strict `<` comparisons), no randomness (S-E).
// ---------------------------------------------------------------------------
inline void sparse_lu_mc64_assignment(
    const std::vector<double>& cost,           // size n*n, row-major
    int                        n,
    std::vector<int>&          match_row_for_col,
    std::vector<double>&       u,               // size n (row potentials)
    std::vector<double>&       v)               // size n (col potentials)
{
    const double INF = std::numeric_limits<double>::max() / 4.0;

    // 1-indexed work arrays (e-maxx convention) to keep the algorithm transparent.
    std::vector<double> U(static_cast<std::size_t>(n) + 1u, 0.0);
    std::vector<double> V(static_cast<std::size_t>(n) + 1u, 0.0);
    std::vector<int>    p(static_cast<std::size_t>(n) + 1u, 0);   // p[j] = row of col j
    std::vector<int>    way(static_cast<std::size_t>(n) + 1u, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(static_cast<std::size_t>(n) + 1u, INF);
        std::vector<char>   used(static_cast<std::size_t>(n) + 1u, 0);
        do {
            used[static_cast<std::size_t>(j0)] = 1;
            const int i0 = p[static_cast<std::size_t>(j0)];
            double delta = INF;
            int    j1    = -1;
            for (int j = 1; j <= n; ++j) {
                if (used[static_cast<std::size_t>(j)]) continue;
                const double cur =
                    cost[static_cast<std::size_t>(i0 - 1) * static_cast<std::size_t>(n)
                         + static_cast<std::size_t>(j - 1)]
                    - U[static_cast<std::size_t>(i0)]
                    - V[static_cast<std::size_t>(j)];
                if (cur < minv[static_cast<std::size_t>(j)]) {
                    minv[static_cast<std::size_t>(j)] = cur;
                    way[static_cast<std::size_t>(j)]  = j0;
                }
                if (minv[static_cast<std::size_t>(j)] < delta) {
                    delta = minv[static_cast<std::size_t>(j)];
                    j1    = j;
                }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[static_cast<std::size_t>(j)]) {
                    U[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] += delta;
                    V[static_cast<std::size_t>(j)]                              -= delta;
                } else {
                    minv[static_cast<std::size_t>(j)] -= delta;
                }
            }
            j0 = j1;
        } while (p[static_cast<std::size_t>(j0)] != 0);
        do {
            const int j1 = way[static_cast<std::size_t>(j0)];
            p[static_cast<std::size_t>(j0)] = p[static_cast<std::size_t>(j1)];
            j0 = j1;
        } while (j0);
    }

    match_row_for_col.assign(static_cast<std::size_t>(n), 0);
    u.assign(static_cast<std::size_t>(n), 0.0);
    v.assign(static_cast<std::size_t>(n), 0.0);
    for (int j = 1; j <= n; ++j) {
        match_row_for_col[static_cast<std::size_t>(j - 1)] =
            p[static_cast<std::size_t>(j)] - 1;       // 0-indexed row matched to col j-1
        v[static_cast<std::size_t>(j - 1)] = V[static_cast<std::size_t>(j)];
    }
    for (int i = 1; i <= n; ++i) {
        u[static_cast<std::size_t>(i - 1)] = U[static_cast<std::size_t>(i)];
    }
}

// ---------------------------------------------------------------------------
// sparse_lu_mc64_make_transform
//
// Computes the MC64 maximum-weight matching + Dr/Dc scaling for A_nat (natural
// CSC, original coordinates) and assembles the effective matrix
//   B = P_static * Dr * A * Dc * Qc
// where Qc is the column ordering passed in col_perm (col_perm[new]=old).
//
// Maximum-PRODUCT matching <=> minimum-cost assignment with
//   cost[i][j] = log(maxcol[j]) - log(|a_ij|)   (>= 0, for a_ij != 0)
//   cost[i][j] = big                            (forbidden: no entry / exact 0)
// The dual potentials u_i, v_j give the standard MC64 scaling
//   Dr[i] = exp(u_i),   Dc[j] = exp(v_j) / maxcol[j]
// so matched entries scale to magnitude 1 and all others to <= 1.
//
// Robustness guard: if any scaling factor is non-finite or non-positive (e.g.
// extreme dynamic range), the scaling falls back to identity.  The matching alone
// still yields a zero-free diagonal, preserving solvability; only equilibration
// is forgone.  This guard never silently changes the answer (solve un-applies
// whatever Dr/Dc are stored).
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_mc64_transform<T, Index>
sparse_lu_mc64_make_transform(
    const csc_storage<T, Index>& A_nat,
    Index                        n,
    const std::vector<Index>&    col_perm)
{
    typedef typename sparse_lu_scalar_policy<T>::real_type real_type;

    sparse_lu_mc64_transform<T, Index> out;
    const std::size_t un = static_cast<std::size_t>(n);

    if (n == Index(0)) {
        out.success = true;
        out.status  = sparse_lu_status::success;
        out.B.col_ptr.assign(1u, Index(0));
        return out;
    }

    const int ni = static_cast<int>(n);

    // ------------------------------------------------------------------
    // 1. Dense magnitude matrix |a_ij| (double) and per-column maxima.
    //    aval[i*n + j] = |A[i,j]| as double (0 if absent / exact zero).
    // ------------------------------------------------------------------
    std::vector<double> aval(un * un, 0.0);
    std::vector<char>   exists(un * un, 0);
    std::vector<double> maxcol(un, 0.0);

    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index k = A_nat.col_ptr[sj]; k < A_nat.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index       i  = A_nat.row_ind[sk];
            if (sparse_lu_scalar_policy<T>::is_exact_zero(A_nat.values[sk])) continue;
            const real_type ar = sparse_lu_scalar_policy<T>::abs_value(A_nat.values[sk]);
            const double    av = static_cast<double>(ar);   // T-independent conversion
            const std::size_t idx =
                static_cast<std::size_t>(i) * un + sj;
            if (av > aval[idx]) {        // keep the largest magnitude if duplicated
                aval[idx] = av;
            }
            exists[idx] = 1;
            if (av > maxcol[sj]) maxcol[sj] = av;
        }
    }

    // ------------------------------------------------------------------
    // 2. Cost matrix with a large finite "forbidden" cost so the assignment
    //    always completes; matched forbidden edges flag structural singularity.
    // ------------------------------------------------------------------
    double finite_max = 0.0;
    for (std::size_t i = 0; i < un; ++i) {
        for (std::size_t j = 0; j < un; ++j) {
            const std::size_t idx = i * un + j;
            if (exists[idx]) {
                const double c = std::log(maxcol[j]) - std::log(aval[idx]);
                if (c > finite_max) finite_max = c;
            }
        }
    }
    const double big = (static_cast<double>(ni) + 1.0) * (finite_max + 1.0) + 1.0;

    std::vector<double> cost(un * un, big);
    for (std::size_t i = 0; i < un; ++i) {
        for (std::size_t j = 0; j < un; ++j) {
            const std::size_t idx = i * un + j;
            if (exists[idx]) {
                cost[idx] = std::log(maxcol[j]) - std::log(aval[idx]);
                if (cost[idx] < 0.0) cost[idx] = 0.0;   // guard tiny negative rounding
            }
        }
    }

    // ------------------------------------------------------------------
    // 3. Minimum-cost perfect assignment (deterministic).
    // ------------------------------------------------------------------
    std::vector<int>    match_row_for_col;
    std::vector<double> u, v;
    sparse_lu_mc64_assignment(cost, ni, match_row_for_col, u, v);

    // Structural singularity: any matched edge is a forbidden (no-entry) slot.
    for (std::size_t j = 0; j < un; ++j) {
        const std::size_t i = static_cast<std::size_t>(match_row_for_col[j]);
        if (!exists[i * un + j]) {
            out.success = false;
            out.status  = sparse_lu_status::structural_singularity;
            return out;
        }
    }

    // ------------------------------------------------------------------
    // 4. MC64 scaling from the dual potentials (original coordinates).
    //    Dr[i] = exp(u_i); Dc[j] = exp(v_j) / maxcol[j].
    //    Fall back to identity if any factor is degenerate.
    // ------------------------------------------------------------------
    std::vector<double> Dr_r(un, 1.0);
    std::vector<double> Dc_r(un, 1.0);
    bool scaling_ok = true;
    for (std::size_t i = 0; i < un && scaling_ok; ++i) {
        const double d = std::exp(u[i]);
        if (!(d > 0.0) || !std::isfinite(d)) scaling_ok = false;
        Dr_r[i] = d;
    }
    for (std::size_t j = 0; j < un && scaling_ok; ++j) {
        if (!(maxcol[j] > 0.0)) { scaling_ok = false; break; }
        const double d = std::exp(v[j]) / maxcol[j];
        if (!(d > 0.0) || !std::isfinite(d)) scaling_ok = false;
        Dc_r[j] = d;
    }
    if (!scaling_ok) {
        std::fill(Dr_r.begin(), Dr_r.end(), 1.0);
        std::fill(Dc_r.begin(), Dc_r.end(), 1.0);
    }

    // ------------------------------------------------------------------
    // 5. Static row permutation composed with the column ordering Qc:
    //    B-row jc holds original row match_row_for_col[col_perm[jc]].
    //    p_static[jc] = that original row;  inv_p_static is its inverse.
    // ------------------------------------------------------------------
    out.p_static.assign(un, Index(0));
    std::vector<Index> inv_p_static(un, Index(0));
    for (Index jc = Index(0); jc < n; ++jc) {
        const std::size_t sjc = static_cast<std::size_t>(jc);
        const Index       oc  = col_perm[sjc];                 // Qc: orig column
        const Index       mr  =
            static_cast<Index>(match_row_for_col[static_cast<std::size_t>(oc)]);
        out.p_static[sjc] = mr;
        inv_p_static[static_cast<std::size_t>(mr)] = jc;
    }

    // Store Dr/Dc as T (original coordinates).
    out.Dr.assign(un, T(0));
    out.Dc.assign(un, T(0));
    for (std::size_t i = 0; i < un; ++i) out.Dr[i] = static_cast<T>(Dr_r[i]);
    for (std::size_t j = 0; j < un; ++j) out.Dc[j] = static_cast<T>(Dc_r[j]);

    // ------------------------------------------------------------------
    // 6. Assemble B = P_static * Dr * A * Dc * Qc (CSC, ascending rows/col).
    //    B[:, jc] = scaled entries of A[:, col_perm[jc]], rows mapped by
    //    inv_p_static.  Diagonal is guaranteed nonzero (matched, scaled entry).
    // ------------------------------------------------------------------
    csc_storage<T, Index>& B = out.B;
    B.col_ptr.assign(un + 1u, Index(0));
    for (Index jc = Index(0); jc < n; ++jc) {
        const std::size_t sjc = static_cast<std::size_t>(jc);
        const Index       oc  = col_perm[sjc];
        const std::size_t soc = static_cast<std::size_t>(oc);
        B.col_ptr[sjc + 1u] = A_nat.col_ptr[soc + 1u] - A_nat.col_ptr[soc];
    }
    for (std::size_t j = 0; j < un; ++j) B.col_ptr[j + 1u] += B.col_ptr[j];

    const Index nnz = B.col_ptr[un];
    B.row_ind.assign(static_cast<std::size_t>(nnz), Index(0));
    B.values.assign(static_cast<std::size_t>(nnz), T(0));

    std::vector<std::pair<Index, T> > colbuf;
    for (Index jc = Index(0); jc < n; ++jc) {
        const std::size_t sjc = static_cast<std::size_t>(jc);
        const Index       oc  = col_perm[sjc];
        const std::size_t soc = static_cast<std::size_t>(oc);
        colbuf.clear();
        for (Index k = A_nat.col_ptr[soc]; k < A_nat.col_ptr[soc + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index       r  = A_nat.row_ind[sk];
            const std::size_t sr = static_cast<std::size_t>(r);
            const Index       br = inv_p_static[sr];           // B-row index
            const T scaled = out.Dr[sr] * A_nat.values[sk] * out.Dc[soc];
            colbuf.push_back(std::make_pair(br, scaled));
        }
        std::sort(colbuf.begin(), colbuf.end(),
                  [](const std::pair<Index, T>& a, const std::pair<Index, T>& b) {
                      return a.first < b.first;
                  });
        Index out_pos = B.col_ptr[sjc];
        for (std::size_t t = 0; t < colbuf.size(); ++t) {
            const std::size_t sp = static_cast<std::size_t>(out_pos);
            B.row_ind[sp] = colbuf[t].first;
            B.values[sp]  = colbuf[t].second;
            ++out_pos;
        }
    }

    out.success = true;
    out.status  = sparse_lu_status::success;
    return out;
}

// ---------------------------------------------------------------------------
// sparse_lu_mc64_matching: matching-only result for the GP-less native path.
//
//   success   : a perfect (zero-free-diagonal) matching was found.
//   status    : success or structural_singularity (no perfect matching).
//   p_static  : p_static[B_row] = original row brought to B_row by the matching
//               composed with the column ordering Qc (== make_transform's
//               out.p_static).  This is ALL the native multifrontal path consumes
//               (Dr/Dc and the assembled B are unused on the native path, MF7
//               §1.3), so this struct deliberately omits them.
// ---------------------------------------------------------------------------
template <class Index>
struct sparse_lu_mc64_matching {
    bool               success;
    sparse_lu_status   status;
    std::vector<Index> p_static;

    sparse_lu_mc64_matching()
        : success(false),
          status(sparse_lu_status::structural_singularity) {}
};

// ---------------------------------------------------------------------------
// sparse_lu_mc64_match_native (SLU-MF8)
//
// Sparse, matching-only MC64 for the native multifrontal path.  Computes EXACTLY
// the same maximum-product matching as sparse_lu_mc64_make_transform but WITHOUT
// ever materializing the dense n*n magnitude/cost matrices: it works on the
// nonzero adjacency only and finds the minimum-cost perfect assignment with a
// successive-shortest-augmenting-path search (Dijkstra + dual potentials, binary
// heap), relaxing real (nonzero) edges only.  Memory is O(n + nnz) instead of the
// dense path's O(n^2), and the per-augmentation work scales with incident
// nonzeros instead of the full n-column rescan.
//
// Matching mathematics are unchanged (maximum-product == minimum-cost assignment
// with cost[i][j] = log(maxcol[j]) - log(|a_ij|), clipped at 0).  Determinism:
// each row's incident columns are visited in ascending column index, and the heap
// breaks distance ties by lower column index, so the matching is reproducible
// (S-E).  On these matrix classes the optimum is unique, so the returned p_static
// is identical to the dense path's (verified by the E4 harness); when a true tie
// admits another equally-optimal matching the residual acceptance gate in
// try_build_native_mc64_supernodal still guarantees full accuracy.
//
// No perfect matching on the real-edge graph -> success=false /
// structural_singularity (S-D).  Magnitudes are read through the scalar policy and
// converted to double exactly as the dense path (T-independent).
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_mc64_matching<Index>
sparse_lu_mc64_match_native(
    const csc_storage<T, Index>& A_nat,
    Index                        n,
    const std::vector<Index>&    col_perm)
{
    sparse_lu_mc64_matching<Index> out;
    const std::size_t un = static_cast<std::size_t>(n);

    if (n == Index(0)) {
        out.success = true;
        out.status  = sparse_lu_status::success;
        return out;
    }

    const int ni = static_cast<int>(n);

    // ------------------------------------------------------------------
    // 1. Per-column maxima (O(nnz), col_ptr scan; same as the dense path).
    // ------------------------------------------------------------------
    std::vector<double> maxcol(un, 0.0);
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index k = A_nat.col_ptr[sj]; k < A_nat.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            if (sparse_lu_scalar_policy<T>::is_exact_zero(A_nat.values[sk])) continue;
            const double av = static_cast<double>(
                sparse_lu_scalar_policy<T>::abs_value(A_nat.values[sk]));
            if (av > maxcol[sj]) maxcol[sj] = av;
        }
    }

    // ------------------------------------------------------------------
    // 2. Row-major nonzero adjacency: row i -> (col j, cost), real edges only.
    //    cost[i][j] = log(maxcol[j]) - log(|a_ij|) >= 0 (clip tiny negatives).
    //    Sorted by column index so the SSP expansion is deterministic.
    // ------------------------------------------------------------------
    std::vector<std::vector<std::pair<int, double> > > adj(un);
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        if (!(maxcol[sj] > 0.0)) continue;          // empty column: no real edge
        const double logmax = std::log(maxcol[sj]);
        for (Index k = A_nat.col_ptr[sj]; k < A_nat.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            if (sparse_lu_scalar_policy<T>::is_exact_zero(A_nat.values[sk])) continue;
            const double av = static_cast<double>(
                sparse_lu_scalar_policy<T>::abs_value(A_nat.values[sk]));
            if (av == 0.0) continue;
            double c = logmax - std::log(av);
            if (c < 0.0) c = 0.0;
            const int i = static_cast<int>(A_nat.row_ind[sk]);
            adj[static_cast<std::size_t>(i)].push_back(
                std::make_pair(static_cast<int>(j), c));
        }
    }
    for (std::size_t i = 0; i < un; ++i) {
        std::sort(adj[i].begin(), adj[i].end(),
                  [](const std::pair<int, double>& a,
                     const std::pair<int, double>& b) { return a.first < b.first; });
    }

    // ------------------------------------------------------------------
    // 3. Minimum-cost perfect assignment via successive shortest augmenting
    //    paths (Dijkstra + column dual potentials, real edges only).
    //    col_match[j] = row matched to column j (or -1).
    // ------------------------------------------------------------------
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> pi_col(un, 0.0);            // column dual potentials
    std::vector<int>    col_match(un, -1);
    typedef std::pair<double, int> QE;              // (reduced distance, col)

    // Binary min-heap over QE: lowest distance first, ties broken by lowest column
    // index (std::pair's lexicographic operator>).  Uses push_heap/pop_heap from
    // <algorithm> (already in scope) -- this file is injected inside namespace vcp,
    // so <queue> must NOT be #included here.
    struct qe_greater {
        bool operator()(const QE& a, const QE& b) const { return a > b; }
    };
    const qe_greater heap_gt = qe_greater();
    std::vector<QE> heap;

    for (int r = 0; r < ni; ++r) {
        std::vector<double> dist(un, INF);
        std::vector<char>   done(un, 0);
        std::vector<int>    par_col(un, -1);        // predecessor column on the path
        std::vector<int>    col_via_row(un, -1);    // row used to settle this column
        heap.clear();

        const std::vector<std::pair<int, double> >& r_adj =
            adj[static_cast<std::size_t>(r)];
        for (std::size_t e = 0; e < r_adj.size(); ++e) {
            const int    j  = r_adj[e].first;
            const double rc = r_adj[e].second - pi_col[static_cast<std::size_t>(j)];
            if (rc < dist[static_cast<std::size_t>(j)]) {
                dist[static_cast<std::size_t>(j)]        = rc;
                col_via_row[static_cast<std::size_t>(j)] = r;
                par_col[static_cast<std::size_t>(j)]     = -1;
                heap.push_back(QE(rc, j));
                std::push_heap(heap.begin(), heap.end(), heap_gt);
            }
        }

        int    free_col = -1;
        double free_d   = 0.0;
        while (!heap.empty()) {
            const QE     top = heap.front();
            std::pop_heap(heap.begin(), heap.end(), heap_gt);
            heap.pop_back();
            const int    j = top.second;
            const double d = top.first;
            if (done[static_cast<std::size_t>(j)]) continue;
            if (d > dist[static_cast<std::size_t>(j)] + 1e-300) continue;
            done[static_cast<std::size_t>(j)] = 1;
            if (col_match[static_cast<std::size_t>(j)] < 0) {
                free_col = j;
                free_d   = d;
                break;
            }
            const int mi = col_match[static_cast<std::size_t>(j)];
            const std::vector<std::pair<int, double> >& m_adj =
                adj[static_cast<std::size_t>(mi)];
            for (std::size_t e = 0; e < m_adj.size(); ++e) {
                const int jj = m_adj[e].first;
                if (done[static_cast<std::size_t>(jj)]) continue;
                const double rc = m_adj[e].second - pi_col[static_cast<std::size_t>(jj)];
                const double nd = d + rc;
                if (nd < dist[static_cast<std::size_t>(jj)]) {
                    dist[static_cast<std::size_t>(jj)]        = nd;
                    col_via_row[static_cast<std::size_t>(jj)] = mi;
                    par_col[static_cast<std::size_t>(jj)]     = j;
                    heap.push_back(QE(nd, jj));
                    std::push_heap(heap.begin(), heap.end(), heap_gt);
                }
            }
        }

        if (free_col < 0) {
            // No augmenting path on the real-edge graph: structurally singular.
            out.success = false;
            out.status  = sparse_lu_status::structural_singularity;
            return out;
        }

        // Dual maintenance: shift settled columns' potentials by their slack.
        for (int j = 0; j < ni; ++j) {
            if (done[static_cast<std::size_t>(j)]) {
                pi_col[static_cast<std::size_t>(j)] +=
                    dist[static_cast<std::size_t>(j)] - free_d;
            }
        }

        // Augment along the path free_col <- ... <- r.
        int j = free_col;
        while (j != -1) {
            const int i  = col_via_row[static_cast<std::size_t>(j)];
            const int pj = par_col[static_cast<std::size_t>(j)];
            col_match[static_cast<std::size_t>(j)] = i;
            j = pj;
        }
    }

    // ------------------------------------------------------------------
    // 4. Compose the static row permutation with the column ordering Qc:
    //    B-row jc holds original row col_match[col_perm[jc]] (== dense path step 5).
    // ------------------------------------------------------------------
    out.p_static.assign(un, Index(0));
    for (Index jc = Index(0); jc < n; ++jc) {
        const std::size_t sjc = static_cast<std::size_t>(jc);
        const Index       oc  = col_perm[sjc];          // Qc: original column
        out.p_static[sjc] =
            static_cast<Index>(col_match[static_cast<std::size_t>(oc)]);
    }

    out.success = true;
    out.status  = sparse_lu_status::success;
    return out;
}

// ---------------------------------------------------------------------------
// sparse_lu_mc64_compose_row_perm
//
// Composes the static MC64 row permutation with the numeric dynamic pivoting
// permutation produced while factorizing B:
//   final_row_perm[new] = p_static[ dynamic_row_perm[new] ]   (-> original row)
// Both arrays are forward permutations (perm[new]=old in their own space).
// ---------------------------------------------------------------------------
template <class Index>
std::vector<Index>
sparse_lu_mc64_compose_row_perm(
    const std::vector<Index>& p_static,
    const std::vector<Index>& dynamic_row_perm)
{
    const std::size_t un = dynamic_row_perm.size();
    std::vector<Index> composed(un, Index(0));
    for (std::size_t i = 0; i < un; ++i) {
        const std::size_t d = static_cast<std::size_t>(dynamic_row_perm[i]);
        composed[i] = p_static[d];
    }
    return composed;
}

} // namespace sparse_lu_detail

#endif // VCP_TSPARSE_SPARSE_LU_MC64_IMPL_HPP
