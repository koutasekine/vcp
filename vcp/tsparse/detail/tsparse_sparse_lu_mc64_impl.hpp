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
//   S-E  Deterministic matching PER SCALAR TYPE T (fixed index tie-break; no
//        randomness).  Cross-type identity of the permutation is NOT guaranteed.
//
// SLU-GT1 D1: the matching/scaling is computed within the requirement set of
// the module scalar contract -- R = real_type<T> arithmetic (+,-,*,/),
// certainly comparisons, and ADL-resolved log/exp only.  No conversion of
// T/R-dependent values to double, and no branch depends on what T is.
// Components whose magnitude cannot be certified positive (!(|a_ij| > 0),
// e.g. an interval containing 0) are treated as STRUCTURALLY ABSENT for the
// matching: this keeps log() arguments certifiably positive and prevents a
// silently degenerate matching (S-D).  "Unreached / not yet computed" states
// in the shortest-path searches are represented by explicit flags, never by
// numeric infinity sentinels (SLU-GT1 P4).
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
// index (strict `<` comparisons), no randomness (S-E, per scalar type).
//
// SLU-GT1 D1/P4: R-generic (R = real_type of the module scalar).  The
// "not yet computed" state of minv[] is an explicit flag (minv_set[]), not a
// numeric infinity sentinel.  For totally ordered R this is behaviorally
// identical to the classic INF initialization (the first relaxation always
// records).  For interval R the certainly-< comparisons may keep a
// suboptimal candidate, but every do-while pass still settles exactly one
// column (the graph is complete, so after the first scan every unused column
// has minv_set), preserving termination and per-T determinism.
// ---------------------------------------------------------------------------
template <class R>
inline void sparse_lu_mc64_assignment(
    const std::vector<R>&      cost,           // size n*n, row-major
    int                        n,
    std::vector<int>&          match_row_for_col,
    std::vector<R>&            u,               // size n (row potentials)
    std::vector<R>&            v)               // size n (col potentials)
{
    // 1-indexed work arrays (e-maxx convention) to keep the algorithm transparent.
    std::vector<R>   U(static_cast<std::size_t>(n) + 1u, R(0));
    std::vector<R>   V(static_cast<std::size_t>(n) + 1u, R(0));
    std::vector<int> p(static_cast<std::size_t>(n) + 1u, 0);   // p[j] = row of col j
    std::vector<int> way(static_cast<std::size_t>(n) + 1u, 0);

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<R>    minv(static_cast<std::size_t>(n) + 1u, R(0));
        std::vector<char> minv_set(static_cast<std::size_t>(n) + 1u, 0);
        std::vector<char> used(static_cast<std::size_t>(n) + 1u, 0);
        do {
            used[static_cast<std::size_t>(j0)] = 1;
            const int i0 = p[static_cast<std::size_t>(j0)];
            R    delta     = R(0);
            bool delta_set = false;
            int  j1        = -1;
            for (int j = 1; j <= n; ++j) {
                if (used[static_cast<std::size_t>(j)]) continue;
                const R cur =
                    cost[static_cast<std::size_t>(i0 - 1) * static_cast<std::size_t>(n)
                         + static_cast<std::size_t>(j - 1)]
                    - U[static_cast<std::size_t>(i0)]
                    - V[static_cast<std::size_t>(j)];
                if (!minv_set[static_cast<std::size_t>(j)] ||
                    cur < minv[static_cast<std::size_t>(j)]) {
                    minv[static_cast<std::size_t>(j)]     = cur;
                    minv_set[static_cast<std::size_t>(j)] = 1;
                    way[static_cast<std::size_t>(j)]      = j0;
                }
                if (minv_set[static_cast<std::size_t>(j)] &&
                    (!delta_set || minv[static_cast<std::size_t>(j)] < delta)) {
                    delta     = minv[static_cast<std::size_t>(j)];
                    delta_set = true;
                    j1        = j;
                }
            }
            for (int j = 0; j <= n; ++j) {
                if (used[static_cast<std::size_t>(j)]) {
                    U[static_cast<std::size_t>(p[static_cast<std::size_t>(j)])] += delta;
                    V[static_cast<std::size_t>(j)]                              -= delta;
                } else if (minv_set[static_cast<std::size_t>(j)]) {
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
    u.assign(static_cast<std::size_t>(n), R(0));
    v.assign(static_cast<std::size_t>(n), R(0));
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
    using std::log;
    using std::exp;

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
    // 1. Dense magnitude matrix |a_ij| (R = real_type) and per-column maxima.
    //    aval[i*n + j] = |A[i,j]| (0 if absent / magnitude not certifiably
    //    positive).  SLU-GT1 D1: components with !(|a_ij| > 0) are treated as
    //    structurally absent, so every stored aval/maxcol entry is certifiably
    //    positive and every log() argument below is certified > 0.
    // ------------------------------------------------------------------
    std::vector<real_type> aval(un * un, real_type(0));
    std::vector<char>      exists(un * un, 0);
    std::vector<real_type> maxcol(un, real_type(0));

    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index k = A_nat.col_ptr[sj]; k < A_nat.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const Index       i  = A_nat.row_ind[sk];
            const real_type ar = sparse_lu_scalar_policy<T>::abs_value(A_nat.values[sk]);
            if (!(ar > real_type(0))) continue;   // not certifiably nonzero: absent
            const std::size_t idx =
                static_cast<std::size_t>(i) * un + sj;
            if (ar > aval[idx]) {        // keep the largest magnitude if duplicated
                aval[idx] = ar;
            }
            exists[idx] = 1;
            if (ar > maxcol[sj]) maxcol[sj] = ar;
        }
    }

    // ------------------------------------------------------------------
    // 2. Cost matrix with a large finite "forbidden" cost so the assignment
    //    always completes; matched forbidden edges flag structural singularity.
    // ------------------------------------------------------------------
    real_type finite_max(0);
    for (std::size_t i = 0; i < un; ++i) {
        for (std::size_t j = 0; j < un; ++j) {
            const std::size_t idx = i * un + j;
            if (exists[idx]) {
                const real_type c = log(maxcol[j]) - log(aval[idx]);
                if (c > finite_max) finite_max = c;
            }
        }
    }
    const real_type big =
        (real_type(ni) + real_type(1)) * (finite_max + real_type(1)) + real_type(1);

    std::vector<real_type> cost(un * un, big);
    for (std::size_t i = 0; i < un; ++i) {
        for (std::size_t j = 0; j < un; ++j) {
            const std::size_t idx = i * un + j;
            if (exists[idx]) {
                cost[idx] = log(maxcol[j]) - log(aval[idx]);
                if (cost[idx] < real_type(0)) cost[idx] = real_type(0);   // guard tiny negative rounding
            }
        }
    }

    // ------------------------------------------------------------------
    // 3. Minimum-cost perfect assignment (deterministic per scalar type).
    // ------------------------------------------------------------------
    std::vector<int>       match_row_for_col;
    std::vector<real_type> u, v;
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
    std::vector<real_type> Dr_r(un, real_type(1));
    std::vector<real_type> Dc_r(un, real_type(1));
    bool scaling_ok = true;
    for (std::size_t i = 0; i < un && scaling_ok; ++i) {
        const real_type d = exp(u[i]);
        if (!(d > real_type(0)) || !vcp::tsparse_scalar::is_finite(d)) scaling_ok = false;
        Dr_r[i] = d;
    }
    for (std::size_t j = 0; j < un && scaling_ok; ++j) {
        if (!(maxcol[j] > real_type(0))) { scaling_ok = false; break; }
        const real_type d = exp(v[j]) / maxcol[j];
        if (!(d > real_type(0)) || !vcp::tsparse_scalar::is_finite(d)) scaling_ok = false;
        Dc_r[j] = d;
    }
    if (!scaling_ok) {
        std::fill(Dr_r.begin(), Dr_r.end(), real_type(1));
        std::fill(Dc_r.begin(), Dc_r.end(), real_type(1));
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

    // Store Dr/Dc as T (original coordinates).  T(real_type) construction:
    // identity for real T, real-part construction for complex T (SLU-GT1 D1).
    out.Dr.assign(un, T(0));
    out.Dc.assign(un, T(0));
    for (std::size_t i = 0; i < un; ++i) out.Dr[i] = T(Dr_r[i]);
    for (std::size_t j = 0; j < un; ++j) out.Dc[j] = T(Dc_r[j]);

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
// structural_singularity (S-D).  SLU-GT1 D1: magnitudes are read through the
// scalar policy and kept in R = real_type<T> (requirement-set arithmetic,
// ADL log; no double conversion).  Components whose magnitude cannot be
// certified positive are treated as structurally absent, exactly as the dense
// path.  Unreached dist[] states are explicit flags, not infinity sentinels
// (SLU-GT1 P4).  Matching is deterministic PER SCALAR TYPE T (S-E).
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_mc64_matching<Index>
sparse_lu_mc64_match_native(
    const csc_storage<T, Index>& A_nat,
    Index                        n,
    const std::vector<Index>&    col_perm)
{
    typedef typename sparse_lu_scalar_policy<T>::real_type real_type;
    using std::log;

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
    //    SLU-GT1 D1: !(|a_ij| > 0) entries are structurally absent, so every
    //    stored maxcol entry is certifiably positive (log-safe below).
    // ------------------------------------------------------------------
    std::vector<real_type> maxcol(un, real_type(0));
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        for (Index k = A_nat.col_ptr[sj]; k < A_nat.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const real_type ar =
                sparse_lu_scalar_policy<T>::abs_value(A_nat.values[sk]);
            if (!(ar > real_type(0))) continue;   // not certifiably nonzero: absent
            if (ar > maxcol[sj]) maxcol[sj] = ar;
        }
    }

    // ------------------------------------------------------------------
    // 2. Row-major nonzero adjacency: row i -> (col j, cost), real edges only.
    //    cost[i][j] = log(maxcol[j]) - log(|a_ij|) >= 0 (clip tiny negatives).
    //    Sorted by column index so the SSP expansion is deterministic.
    // ------------------------------------------------------------------
    std::vector<std::vector<std::pair<int, real_type> > > adj(un);
    for (Index j = Index(0); j < n; ++j) {
        const std::size_t sj = static_cast<std::size_t>(j);
        if (!(maxcol[sj] > real_type(0))) continue; // empty column: no real edge
        const real_type logmax = log(maxcol[sj]);
        for (Index k = A_nat.col_ptr[sj]; k < A_nat.col_ptr[sj + 1u]; ++k) {
            const std::size_t sk = static_cast<std::size_t>(k);
            const real_type ar =
                sparse_lu_scalar_policy<T>::abs_value(A_nat.values[sk]);
            if (!(ar > real_type(0))) continue;   // certified-positive guard before log
            real_type c = logmax - log(ar);
            if (c < real_type(0)) c = real_type(0);
            const int i = static_cast<int>(A_nat.row_ind[sk]);
            adj[static_cast<std::size_t>(i)].push_back(
                std::make_pair(static_cast<int>(j), c));
        }
    }
    for (std::size_t i = 0; i < un; ++i) {
        std::sort(adj[i].begin(), adj[i].end(),
                  [](const std::pair<int, real_type>& a,
                     const std::pair<int, real_type>& b) { return a.first < b.first; });
    }

    // ------------------------------------------------------------------
    // 3. Minimum-cost perfect assignment via successive shortest augmenting
    //    paths (Dijkstra + column dual potentials, real edges only).
    //    col_match[j] = row matched to column j (or -1).
    //    SLU-GT1 P4: "column not reached yet" is the explicit flag
    //    dist_set[j] == 0, not an infinity sentinel.  For totally ordered R
    //    this is behaviorally identical (the first relaxation always records);
    //    for interval R indeterminate certainly-< keeps the incumbent, and
    //    termination is preserved (each pop settles at most one column).
    // ------------------------------------------------------------------
    std::vector<real_type> pi_col(un, real_type(0)); // column dual potentials
    std::vector<int>       col_match(un, -1);
    typedef std::pair<real_type, int> QE;            // (reduced distance, col)

    // Hand-rolled binary min-heap over QE: lowest distance first, ties broken by
    // lowest column index.  The order predicate is the same lexicographic pair
    // comparison the previous std::push_heap/pop_heap usage induced:
    //   qe_before(a, b) = a.first < b.first
    //                  || (!(b.first < a.first) && a.second < b.second)
    // It deliberately falls through to the column tie-break when the distances
    // are equal OR certainly-incomparable, so interval distances that cannot be
    // ordered are still decided deterministically by column index (S-E: per-T
    // determinism via fixed tie-break is preserved).
    //
    // SLU-GT1.1 F-1 (B-7): this hand-rolled heap does NOT require the predicate
    // to be a strict weak ordering.  Under R's certainly comparisons (a partial
    // order: incomparability is not transitive) the extraction order is
    // best-effort, and any ordering degradation affects only matching quality,
    // never correctness.  std::push_heap/pop_heap (and other std algorithms
    // with comparator requirements) are avoided because passing a non-SWO
    // predicate to them is formally UB.  This file is injected inside
    // namespace vcp, so <queue> must NOT be #included here either.
    std::vector<QE> heap;
    const auto qe_before = [](const QE& a, const QE& b) -> bool {
        return a.first < b.first ||
               (!(b.first < a.first) && a.second < b.second);
    };
    // sift-up (push): iterative, no recursion.
    const auto heap_sift_up = [&heap, &qe_before](std::size_t c) {
        while (c > 0u) {
            const std::size_t p = (c - 1u) / 2u;
            if (!qe_before(heap[c], heap[p])) break;
            const QE tmp = heap[c]; heap[c] = heap[p]; heap[p] = tmp;
            c = p;
        }
    };
    // sift-down (pop): iterative, no recursion.
    const auto heap_sift_down = [&heap, &qe_before](std::size_t c) {
        const std::size_t sz = heap.size();
        for (;;) {
            const std::size_t l = 2u * c + 1u;
            if (l >= sz) break;
            std::size_t m = l;
            const std::size_t rt = l + 1u;
            if (rt < sz && qe_before(heap[rt], heap[l])) m = rt;
            if (!qe_before(heap[m], heap[c])) break;
            const QE tmp = heap[c]; heap[c] = heap[m]; heap[m] = tmp;
            c = m;
        }
    };

    for (int r = 0; r < ni; ++r) {
        std::vector<real_type> dist(un, real_type(0));
        std::vector<char>      dist_set(un, 0);     // P4 flag: dist[j] valid iff set
        std::vector<char>      done(un, 0);
        std::vector<int>       par_col(un, -1);     // predecessor column on the path
        std::vector<int>       col_via_row(un, -1); // row used to settle this column
        heap.clear();

        const std::vector<std::pair<int, real_type> >& r_adj =
            adj[static_cast<std::size_t>(r)];
        for (std::size_t e = 0; e < r_adj.size(); ++e) {
            const int       j  = r_adj[e].first;
            const real_type rc = r_adj[e].second - pi_col[static_cast<std::size_t>(j)];
            if (!dist_set[static_cast<std::size_t>(j)] ||
                rc < dist[static_cast<std::size_t>(j)]) {
                dist[static_cast<std::size_t>(j)]        = rc;
                dist_set[static_cast<std::size_t>(j)]    = 1;
                col_via_row[static_cast<std::size_t>(j)] = r;
                par_col[static_cast<std::size_t>(j)]     = -1;
                heap.push_back(QE(rc, j));
                heap_sift_up(heap.size() - 1u);
            }
        }

        int       free_col = -1;
        real_type free_d(0);
        while (!heap.empty()) {
            const QE        top = heap.front();
            if (heap.size() > 1u) heap.front() = heap.back();
            heap.pop_back();
            if (!heap.empty()) heap_sift_down(0u);
            const int       j = top.second;
            const real_type d = top.first;
            if (done[static_cast<std::size_t>(j)]) continue;
            // Stale heap entry (superseded by a later, better relaxation).
            // Heap entries are only pushed with dist_set[j] == 1.
            if (dist_set[static_cast<std::size_t>(j)] &&
                (d > dist[static_cast<std::size_t>(j)] + real_type(1e-300))) continue;
            done[static_cast<std::size_t>(j)] = 1;
            if (col_match[static_cast<std::size_t>(j)] < 0) {
                free_col = j;
                free_d   = d;
                break;
            }
            const int mi = col_match[static_cast<std::size_t>(j)];
            const std::vector<std::pair<int, real_type> >& m_adj =
                adj[static_cast<std::size_t>(mi)];
            for (std::size_t e = 0; e < m_adj.size(); ++e) {
                const int jj = m_adj[e].first;
                if (done[static_cast<std::size_t>(jj)]) continue;
                const real_type rc = m_adj[e].second - pi_col[static_cast<std::size_t>(jj)];
                const real_type nd = d + rc;
                if (!dist_set[static_cast<std::size_t>(jj)] ||
                    nd < dist[static_cast<std::size_t>(jj)]) {
                    dist[static_cast<std::size_t>(jj)]        = nd;
                    dist_set[static_cast<std::size_t>(jj)]    = 1;
                    col_via_row[static_cast<std::size_t>(jj)] = mi;
                    par_col[static_cast<std::size_t>(jj)]     = j;
                    heap.push_back(QE(nd, jj));
                    heap_sift_up(heap.size() - 1u);
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
