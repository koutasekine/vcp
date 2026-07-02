// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLU-8R.5.5: A_eff-origin true numeric supernodal factorization -- internal implementation.
//
// This file owns the numeric source switch for method=supernodal.
//
// Contract:
//   panel_values / U_segments are initialized from A_eff (NOT from CSC L/U).
//   §17.2(A) and §17.2(B) are run in topological (per-supernode) interleaved order.
//   On success: storage.true_numeric_source = true.
//               storage.numeric_source_kind = a_eff_true_numeric.
//   On failure: panel_values / U_segments.values may be modified in-place
//               (zeroed and scatter-filled from A_eff before the failure point).
//               However, true_numeric_source remains false, and native solve
//               is disabled (CSC fallback remains in effect).
//               This is an in-place may-modify contract (SLU-8R.5.5.2 Option B).
//               See sandbox/docs/reviews/SLU8R55_numeric_source_switch_review.md.
//
// PROHIBITIONS (§25):
//   dense_kernel::getrf MUST NOT be called for pivot search in §17.2(B).
//   within_panel_factorization uses explicit select_threshold_pivot (no getrf).
//
// Gate 6 status: PASS candidate only.
//   Full Gate 6 confirmation requires SLU-8R.6 review.
//   issue_SLU8_contract_violation.md: OPEN.
//   SLU-8 full conformance: NOT claimed.
//
// This file MUST be #included from WITHIN namespace vcp, AFTER:
//   - supernodal_lu_storage<T, Index>, supernode_desc<Index>, csc_storage<T, Index>
//   - sparse_lu_detail::csc_lookup_value (from bootstrap_impl.hpp)
//   - sparse_lu_detail::compute_panel_update_set, gather_panel_workspace,
//     apply_supernode_panel_update, scatter_panel_workspace (from panel_update_impl.hpp)
//   - sparse_lu_detail::gather_within_panel_workspace, select_threshold_pivot,
//     apply_within_panel_row_swap, scale_pivot_column, update_trailing_columns,
//     scatter_within_panel_factor (from within_panel_factor_impl.hpp)
//   - sparse_lu_detail::supernodal_true_numeric_stats, within_panel_factor_stats,
//     supernode_panel_update_stats (from tsparse_sparse_lu.hpp)
//   - sparse_lu_scalar_policy<T>, sparse_lu_options<T>
//
// Do NOT include this file directly. Include:
//   <vcp/tsparse/tsparse_sparse_lu.hpp>

#ifndef VCP_TSPARSE_SPARSE_LU_TRUE_NUMERIC_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LU_TRUE_NUMERIC_IMPL_HPP

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace sparse_lu_detail {

// ---------------------------------------------------------------------------
// a_eff_lookup
//
// Returns A_eff[new_row, new_col] = Dr[new_row] * A_csc[old_row, new_col] * Dc[new_col]
//   where old_row = storage.row_perm[new_row]
//   A_csc has column permutation already applied (col_perm baked in).
//
// Returns T(0) if the entry is structurally zero in A_csc.
// ---------------------------------------------------------------------------
template <class T, class Index>
T a_eff_lookup(
    const csc_storage<T, Index>&           A_csc,
    const supernodal_lu_storage<T, Index>& storage,
    Index new_row,
    Index new_col)
{
    const std::size_t snr = static_cast<std::size_t>(new_row);
    const std::size_t snc = static_cast<std::size_t>(new_col);

    if (snr >= storage.row_perm.size()) return T(0);
    const Index old_row = storage.row_perm[snr];

    // A_csc[old_row, new_col] -- binary search on sorted row_ind
    const T a_val = csc_lookup_value(
        A_csc.col_ptr, A_csc.row_ind, A_csc.values, new_col, old_row);

    const T dr = (!storage.Dr.empty() && snr < storage.Dr.size())
        ? storage.Dr[snr] : T(1);
    const T dc = (!storage.Dc.empty() && snc < storage.Dc.size())
        ? storage.Dc[snc] : T(1);

    return dr * a_val * dc;
}

// ---------------------------------------------------------------------------
// initialize_supernodal_from_a_eff
//
// Resets all panel_values and U_segments.values to zero, then scatters
// A_eff values into the positions dictated by the supernode structure.
//
// Uses storage.row_perm (initial, from CSC bootstrap) for old_row lookup.
// Dr / Dc (if non-empty) are applied.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool initialize_supernodal_from_a_eff(
    const csc_storage<T, Index>&           A_csc,
    supernodal_lu_storage<T, Index>&       storage,
    supernodal_true_numeric_stats&         stats)
{
    std::fill(storage.panel_values.begin(),     storage.panel_values.end(),     T(0));
    std::fill(storage.U_segments.values.begin(), storage.U_segments.values.end(), T(0));

    const std::size_t nsup = storage.supernodes.size();
    std::size_t scatter_count = 0u;

    for (std::size_t j = 0u; j < nsup; ++j) {
        const supernode_desc<Index>& desc = storage.supernodes[j];
        const Index col_begin = desc.first_col;
        const Index ld        = desc.leading_dimension;
        const Index R         = static_cast<Index>(desc.row_indices.size());
        const std::size_t offset = static_cast<std::size_t>(desc.values_offset);

        if (offset + static_cast<std::size_t>(ld) *
                     static_cast<std::size_t>(desc.num_cols) >
                     storage.panel_values.size()) continue;

        T* panel = &storage.panel_values[offset];

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            const std::size_t sc = static_cast<std::size_t>(c);
            T* col_data = panel + sc * static_cast<std::size_t>(ld);

            for (Index r = Index(0); r < R; ++r) {
                const Index new_row = desc.row_indices[static_cast<std::size_t>(r)];
                col_data[static_cast<std::size_t>(r)] =
                    a_eff_lookup(A_csc, storage, new_row, new_col);
                ++scatter_count;
            }
        }

        // U_segments: off-diagonal U positions (u_row < col_begin)
        if (j >= storage.U_segments.seg_ptr.size()) continue;
        const Index j_seg_start = storage.U_segments.seg_ptr[j];
        if (desc.u_seg_col_ptr.size() <
                static_cast<std::size_t>(desc.num_cols + 1)) continue;

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            const std::size_t sc = static_cast<std::size_t>(c);
            const Index rel_s = desc.u_seg_col_ptr[sc];
            const Index rel_e = desc.u_seg_col_ptr[sc + 1u];

            for (Index ki = rel_s; ki < rel_e; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(j_seg_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                storage.U_segments.values[abs_idx] =
                    a_eff_lookup(A_csc, storage, u_row, new_col);
                ++scatter_count;
            }
        }
    }

    stats.a_entry_scatter_count = scatter_count;
    return true;
}

// ---------------------------------------------------------------------------
// factorize_within_panel_single
//
// Applies §17.2(B) to a single supernode j.
// Reuses the building-block functions from within_panel_factor_impl.hpp:
//   gather_within_panel_workspace, select_threshold_pivot,
//   apply_within_panel_row_swap, scale_pivot_column,
//   update_trailing_columns, scatter_within_panel_factor.
//
// Returns false if a zero pivot is encountered.
// stats: accumulates pivot/ger counts.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool factorize_within_panel_single(
    supernodal_lu_storage<T, Index>& storage,
    Index j,
    const sparse_lu_options<T>&      opt,
    within_panel_factor_stats&       stats)
{
    typedef sparse_lu_scalar_policy<T> scalar_pol;
    typedef typename scalar_pol::real_type real_type;

    const std::size_t sj = static_cast<std::size_t>(j);
    if (sj >= storage.supernodes.size()) return false;

    const supernode_desc<Index>& desc_j = storage.supernodes[sj];
    const Index w = desc_j.num_cols;
    const Index R = static_cast<Index>(desc_j.row_indices.size());

    if (w <= Index(0) || R <= Index(0)) return true;

    within_panel_workspace<T, Index> work;
    gather_within_panel_workspace(storage, j, work);
    if (!work.valid) return false;

    bool completed = true;
    for (Index c = Index(0); c < w; ++c) {
        stats.pivot_search_count++;

        Index pivot_row;
        T     pivot_value;
        pivot_decision decision =
            select_threshold_pivot(work, c, opt, pivot_row, pivot_value);

        if (decision == pivot_decision::acceptable) {
            stats.pivot_accept_count++;
        } else if (decision == pivot_decision::reject) {
            stats.pivot_reject_count++;
        } else {
            stats.inconclusive_pivot_count++;
        }

        if (pivot_row != c) {
            apply_within_panel_row_swap(work, c, pivot_row, storage, stats);
        }

        const T& actual_pivot =
            work.values[static_cast<std::size_t>(c * work.ld + c)];
        const real_type abs_pv = scalar_pol::abs_value(actual_pivot);

        if (abs_pv <= opt.zero_tolerance) {
            stats.zero_pivot_count++;
            completed = false;
            continue;
        }
        if (abs_pv <= opt.near_zero_tolerance) {
            stats.near_zero_pivot_count++;
        }

        scale_pivot_column(work, c, actual_pivot, stats);
        update_trailing_columns(work, c, stats);
    }

    scatter_within_panel_factor(storage, j, work, stats);
    stats.panel_count++;

    return completed;
}

// ---------------------------------------------------------------------------
// compute_supernodal_lu_residual  (DENSE REFERENCE)
//
// Reconstructs dense L (unit lower) and U (upper) from supernodal storage and
// computes the point-arithmetic Frobenius residual:
//   abs_res = ||L*U - A_eff||_F
//   rel_res = abs_res / max(1.0, ||A_eff||_F)
//
// SLU-8R.8: this dense O(n^3) path is RETAINED ONLY as a small-n correctness
// reference for the dense-vs-sparse cross-check (see
// compute_supernodal_lu_residual_sparse). It is NOT used on the production
// acceptance path. Production residual checking is scalable and never skipped
// (see the sparse variant below). Because it is reference-only, the dense
// reconstruction is capped at sn_residual_dense_reference_n_max; for larger n
// it returns false and the caller MUST use the sparse variant.
//
// NOTE: this is a point-arithmetic residual sanity check (Frobenius norm).
// ---------------------------------------------------------------------------
static const std::size_t sn_residual_dense_reference_n_max = 200u;

template <class T, class Index>
bool compute_supernodal_lu_residual(
    const csc_storage<T, Index>&           A_csc,
    const supernodal_lu_storage<T, Index>& storage,
    Index n,
    double& abs_res,
    double& rel_res)
{
    typedef sparse_lu_scalar_policy<T> scalar_pol;
    abs_res = 0.0; rel_res = 0.0;

    if (static_cast<std::size_t>(n) > sn_residual_dense_reference_n_max ||
        n <= Index(0))
        return false;

    const std::size_t sn   = static_cast<std::size_t>(n);
    const std::size_t nsup = storage.supernodes.size();

    // Dense L (unit lower) and U (upper), column-major.
    std::vector<T> L(sn * sn, T(0));
    std::vector<T> U(sn * sn, T(0));

    // L diagonal = identity (unit lower)
    for (std::size_t i = 0u; i < sn; ++i) L[i * sn + i] = T(1);

    for (std::size_t j = 0u; j < nsup; ++j) {
        const supernode_desc<Index>& desc = storage.supernodes[j];
        const Index col_begin = desc.first_col;
        const Index ld        = desc.leading_dimension;
        const std::size_t R   = desc.row_indices.size();
        const std::size_t off = static_cast<std::size_t>(desc.values_offset);

        if (off + static_cast<std::size_t>(ld) *
                  static_cast<std::size_t>(desc.num_cols) >
                  storage.panel_values.size()) continue;

        const T* panel = &storage.panel_values[off];

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            const std::size_t sc  = static_cast<std::size_t>(c);
            const std::size_t snc = static_cast<std::size_t>(new_col);
            const T* col_data = panel + sc * static_cast<std::size_t>(ld);

            for (std::size_t r = 0u; r < R; ++r) {
                const Index new_row = desc.row_indices[r];
                const std::size_t snr = static_cast<std::size_t>(new_row);
                const T val = col_data[r];

                if (new_row <= new_col) {
                    // On or above diagonal → U entry
                    U[snc * sn + snr] = val;
                } else {
                    // Below diagonal → L multiplier
                    L[snc * sn + snr] = val;
                }
            }
        }

        // U_segments: off-diagonal U (rows < col_begin)
        if (j >= storage.U_segments.seg_ptr.size()) continue;
        const Index seg_start = storage.U_segments.seg_ptr[j];
        if (desc.u_seg_col_ptr.size() <
                static_cast<std::size_t>(desc.num_cols + 1)) continue;

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            const std::size_t snc = static_cast<std::size_t>(new_col);
            const std::size_t sc  = static_cast<std::size_t>(c);
            const Index rel_s = desc.u_seg_col_ptr[sc];
            const Index rel_e = desc.u_seg_col_ptr[sc + 1u];

            for (Index ki = rel_s; ki < rel_e; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(seg_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                const std::size_t snr = static_cast<std::size_t>(u_row);
                U[snc * sn + snr] = storage.U_segments.values[abs_idx];
            }
        }
    }

    // Compute L*U (column-major: (L*U)[:,c] = sum_k L[:,k] * U[k,c])
    std::vector<T> LU(sn * sn, T(0));
    for (std::size_t c = 0u; c < sn; ++c) {
        for (std::size_t r = 0u; r < sn; ++r) {
            T sum = T(0);
            for (std::size_t k = 0u; k < sn; ++k) {
                sum = sum + L[k * sn + r] * U[c * sn + k];
            }
            LU[c * sn + r] = sum;
        }
    }

    // Compute Frobenius residual against A_eff
    double res_sq  = 0.0;
    double aeff_sq = 0.0;
    for (std::size_t new_col = 0u; new_col < sn; ++new_col) {
        for (std::size_t new_row = 0u; new_row < sn; ++new_row) {
            const T aeff = a_eff_lookup(
                A_csc, storage,
                static_cast<Index>(new_row),
                static_cast<Index>(new_col));
            const T diff = LU[new_col * sn + new_row] - aeff;

            const double da = static_cast<double>(scalar_pol::abs_value(aeff));
            const double dd = static_cast<double>(scalar_pol::abs_value(diff));
            aeff_sq += da * da;
            res_sq  += dd * dd;
        }
    }

    abs_res = std::sqrt(res_sq);
    const double aeff_norm = std::sqrt(aeff_sq);
    rel_res = abs_res / ((aeff_norm > 1.0) ? aeff_norm : 1.0);
    return true;
}

// ---------------------------------------------------------------------------
// compute_supernodal_lu_residual_sparse  (PRODUCTION, SCALABLE)
//
// SLU-8R.8 B2: scalable sparse/column-accumulator residual checking that does
// NOT skip for large n. This is the production acceptance-path residual checker.
//
// Algorithm (column accumulator):
//   For each global column j:
//     acc := A_eff(:, j)                      (sparse, only A nonzeros)
//     For each k with U(k, j) != 0:
//         acc -= U(k, j) * L(:, k)            (L(k,k) = 1, plus strictly-lower)
//     res_sq  += ||acc||_2^2
//     aeff_sq += ||A_eff(:, j)||_2^2
//
//   abs_res = sqrt(res_sq) = ||A_eff - L*U||_F
//   rel_res = abs_res / max(1.0, sqrt(aeff_sq)) = abs_res / max(1, ||A_eff||_F)
//
// A_eff semantics: identical to a_eff_lookup. A_eff(:,j) nonzeros are produced
// from A_csc column j (col_perm already baked in) mapped through the SAME
// final post-pivot row_perm / Dr / Dc carried in `storage`:
//   A_eff[new_row, j] = Dr[new_row] * A_csc[old_row, j] * Dc[j],
//   old_row = storage.row_perm[new_row]  (inverse used to scatter by column).
// L and U are read from the SAME final supernodal storage (panel_values +
// U_segments), so A_eff and L*U share one storage/permutation/scaling source.
// This is the SLU-8R.7/R.7.1 invariant that prevents B1-type mismatch.
//
// Workspace is O(n): a dense `work` accumulator cleared via a `touched` list,
// no dense n x n matrix is ever formed.
//
// NOTE: point-arithmetic Frobenius residual sanity check.
// ---------------------------------------------------------------------------
template <class T, class Index>
bool compute_supernodal_lu_residual_sparse(
    const csc_storage<T, Index>&           A_csc,
    const supernodal_lu_storage<T, Index>& storage,
    Index n,
    double& abs_res,
    double& rel_res)
{
    typedef sparse_lu_scalar_policy<T> scalar_pol;
    abs_res = 0.0; rel_res = 0.0;

    if (n <= Index(0)) return false;
    if (A_csc.col_ptr.size() < static_cast<std::size_t>(n) + 1u) return false;

    const std::size_t sn   = static_cast<std::size_t>(n);
    const std::size_t nsup = storage.supernodes.size();

    typedef std::pair<Index, T> entry_type;

    // Per-column sparse L (strictly lower, unit diagonal implicit) and per-column
    // sparse U (rows <= col, includes the pivot diagonal). Built once from the
    // final supernodal storage; total size is O(nnz of L+U).
    std::vector<std::vector<entry_type> > Lcol(sn);
    std::vector<std::vector<entry_type> > Ucol(sn);

    for (std::size_t j = 0u; j < nsup; ++j) {
        const supernode_desc<Index>& desc = storage.supernodes[j];
        const Index col_begin = desc.first_col;
        const Index ld        = desc.leading_dimension;
        const std::size_t R   = desc.row_indices.size();
        const std::size_t off = static_cast<std::size_t>(desc.values_offset);

        if (off + static_cast<std::size_t>(ld) *
                  static_cast<std::size_t>(desc.num_cols) >
                  storage.panel_values.size()) continue;

        const T* panel = &storage.panel_values[off];

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            if (new_col < Index(0) ||
                static_cast<std::size_t>(new_col) >= sn) continue;
            const std::size_t sc = static_cast<std::size_t>(c);
            const T* col_data = panel + sc * static_cast<std::size_t>(ld);

            for (std::size_t r = 0u; r < R; ++r) {
                const Index new_row = desc.row_indices[r];
                if (new_row < Index(0) ||
                    static_cast<std::size_t>(new_row) >= sn) continue;
                const T val = col_data[r];
                if (new_row <= new_col) {
                    Ucol[static_cast<std::size_t>(new_col)].push_back(
                        entry_type(new_row, val));
                } else {
                    Lcol[static_cast<std::size_t>(new_col)].push_back(
                        entry_type(new_row, val));
                }
            }
        }

        // U_segments: off-diagonal U (rows < col_begin)
        if (j >= storage.U_segments.seg_ptr.size()) continue;
        const Index seg_start = storage.U_segments.seg_ptr[j];
        if (desc.u_seg_col_ptr.size() <
                static_cast<std::size_t>(desc.num_cols + 1)) continue;

        for (Index c = Index(0); c < desc.num_cols; ++c) {
            const Index new_col = col_begin + c;
            if (new_col < Index(0) ||
                static_cast<std::size_t>(new_col) >= sn) continue;
            const std::size_t sc = static_cast<std::size_t>(c);
            const Index rel_s = desc.u_seg_col_ptr[sc];
            const Index rel_e = desc.u_seg_col_ptr[sc + 1u];

            for (Index ki = rel_s; ki < rel_e; ++ki) {
                const std::size_t abs_idx =
                    static_cast<std::size_t>(seg_start + ki);
                if (abs_idx >= storage.U_segments.row_ind.size()) continue;
                const Index u_row = storage.U_segments.row_ind[abs_idx];
                if (u_row < Index(0) ||
                    static_cast<std::size_t>(u_row) >= sn) continue;
                Ucol[static_cast<std::size_t>(new_col)].push_back(
                    entry_type(u_row, storage.U_segments.values[abs_idx]));
            }
        }
    }

    // inverse row permutation: inv_row_perm[old_row] = new_row.
    std::vector<Index> inv_row_perm(sn, Index(-1));
    for (std::size_t nr = 0u;
         nr < storage.row_perm.size() && nr < sn; ++nr) {
        const Index old_row = storage.row_perm[nr];
        if (old_row >= Index(0) && static_cast<std::size_t>(old_row) < sn)
            inv_row_perm[static_cast<std::size_t>(old_row)] =
                static_cast<Index>(nr);
    }

    // Column accumulator workspace (O(n)).
    std::vector<T>     work(sn, T(0));
    std::vector<char>  marked(sn, 0);
    std::vector<Index> touched;
    touched.reserve(64);

    double res_sq  = 0.0;
    double aeff_sq = 0.0;

    for (std::size_t jcol = 0u; jcol < sn; ++jcol) {
        touched.clear();

        // acc := A_eff(:, jcol)
        const Index cs = A_csc.col_ptr[jcol];
        const Index ce = A_csc.col_ptr[jcol + 1u];
        const T dc = (!storage.Dc.empty() && jcol < storage.Dc.size())
            ? storage.Dc[jcol] : T(1);

        for (Index idx = cs; idx < ce; ++idx) {
            const std::size_t sidx = static_cast<std::size_t>(idx);
            if (sidx >= A_csc.row_ind.size()) continue;
            const Index old_row = A_csc.row_ind[sidx];
            if (old_row < Index(0) ||
                static_cast<std::size_t>(old_row) >= sn) continue;
            const Index new_row = inv_row_perm[static_cast<std::size_t>(old_row)];
            if (new_row < Index(0)) continue;
            const std::size_t snr = static_cast<std::size_t>(new_row);
            const T dr = (!storage.Dr.empty() && snr < storage.Dr.size())
                ? storage.Dr[snr] : T(1);
            const T aval = dr * A_csc.values[sidx] * dc;

            if (!marked[snr]) { marked[snr] = 1; touched.push_back(new_row); }
            work[snr] = work[snr] + aval;

            const double da = static_cast<double>(scalar_pol::abs_value(aval));
            aeff_sq += da * da;
        }

        // acc -= sum_k U(k, jcol) * L(:, k)   with L(k,k) = 1
        const std::vector<entry_type>& ucolj = Ucol[jcol];
        for (std::size_t t = 0u; t < ucolj.size(); ++t) {
            const Index k    = ucolj[t].first;
            const T     u_kj = ucolj[t].second;
            const std::size_t sk = static_cast<std::size_t>(k);

            // diagonal contribution L(k,k) = 1
            if (!marked[sk]) { marked[sk] = 1; touched.push_back(k); }
            work[sk] = work[sk] - u_kj;

            // strictly-lower L(:, k)
            const std::vector<entry_type>& lcolk = Lcol[sk];
            for (std::size_t s = 0u; s < lcolk.size(); ++s) {
                const Index r = lcolk[s].first;
                const std::size_t sr = static_cast<std::size_t>(r);
                if (!marked[sr]) { marked[sr] = 1; touched.push_back(r); }
                work[sr] = work[sr] - u_kj * lcolk[s].second;
            }
        }

        // accumulate ||acc||^2 and clear workspace via touched list
        for (std::size_t t = 0u; t < touched.size(); ++t) {
            const std::size_t sr = static_cast<std::size_t>(touched[t]);
            const double dv =
                static_cast<double>(scalar_pol::abs_value(work[sr]));
            res_sq += dv * dv;
            work[sr]   = T(0);
            marked[sr] = 0;
        }
    }

    abs_res = std::sqrt(res_sq);
    const double aeff_norm = std::sqrt(aeff_sq);
    rel_res = abs_res / ((aeff_norm > 1.0) ? aeff_norm : 1.0);
    return true;
}

} // namespace sparse_lu_detail

// ---------------------------------------------------------------------------
// sparse_lu_factorize_supernodal_from_a_eff  (public, in namespace vcp)
//
// SLU-8R.5.5 primary driver: initializes supernodal storage from A_eff and
// runs the interleaved §17.2(A)/(B) factorization.
//
// Returns supernodal_true_numeric_stats describing outcome.
// On success: storage.true_numeric_source = true.
// On failure: panel_values / U_segments.values may have been modified in-place
//   (zeroed and scatter-filled from A_eff before the failure point).
//   In-place may-modify contract (SLU-8R.5.5.2 Option B):
//     - true_numeric_source remains false
//     - native solve remains disabled
//     - CSC fallback continues to be used
//     - caller-visible production fac.solve() correctness is maintained via fallback
//
// GATE 6: PASS candidate.
//   Full Gate 6 requires SLU-8R.6 conformance review.
//   issue_SLU8_contract_violation.md: OPEN.
// ---------------------------------------------------------------------------
template <class T, class Index>
sparse_lu_detail::supernodal_true_numeric_stats
sparse_lu_factorize_supernodal_from_a_eff(
    const csc_storage<T, Index>&           A_csc,
    const baseline_lu_storage<T, Index>&   csc_lu,
    supernodal_lu_storage<T, Index>&       storage,
    const sparse_lu_options<T>&            opt)
{
    (void)csc_lu; // row_perm/Dr/Dc already in storage (copied from csc_lu at bootstrap)

    const auto t_total_start = std::chrono::steady_clock::now();

    sparse_lu_detail::supernodal_true_numeric_stats stats;
    stats.attempted = true;

    // SLU-8R.6.2: helper to set total_ticks on every return path without repetition.
    // Zero is legitimate (clock resolution); no max(1,...) forced nonzero.
    const auto finish_total_ticks = [&]() {
        stats.total_ticks = static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_total_start).count());
        // SLU-SN-OPT: pure factorization time excludes the post-factorization
        // residual verification phase (design §17.2 line 894: "分解時間").
        stats.factorization_ticks =
            (stats.total_ticks > stats.residual_ticks)
                ? (stats.total_ticks - stats.residual_ticks) : 0u;
    };

    if (!storage.valid || storage.supernodes.empty()) {
        stats.status =
            supernodal_true_numeric_status::unsupported_structure;
        finish_total_ticks();
        return stats;
    }

    const std::size_t nsup = storage.supernodes.size();
    const Index n = static_cast<Index>(storage.row_perm.size());

    // SLU-SN-OPT: precompute col -> supernode map once (O(n)). Used to replace the
    // O(nsup^2) per-supernode scan in compute_panel_update_set with an O(#U-rows)
    // lookup. The resulting update set (and its ascending order) is identical, so
    // the numeric factorization is byte-identical.
    std::vector<Index> col_to_supernode(static_cast<std::size_t>(n < Index(0) ? Index(0) : n), Index(-1));
    for (std::size_t s = 0u; s < nsup; ++s) {
        const supernode_desc<Index>& d = storage.supernodes[s];
        const Index cb = d.first_col;
        const Index ce = cb + d.num_cols;
        for (Index c = cb; c < ce; ++c) {
            if (c >= Index(0) && static_cast<std::size_t>(c) < col_to_supernode.size())
                col_to_supernode[static_cast<std::size_t>(c)] = static_cast<Index>(s);
        }
    }

    // ------------------------------------------------------------------
    // Step 1: Initialize panel_values / U_segments from A_eff
    // ------------------------------------------------------------------
    {
        const auto t_init0 = std::chrono::steady_clock::now();
        bool init_ok = sparse_lu_detail::initialize_supernodal_from_a_eff(
            A_csc, storage, stats);
        stats.init_scatter_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_init0).count());
        if (!init_ok) {
            stats.status =
                supernodal_true_numeric_status::missing_a_entry;
            finish_total_ticks();
            return stats;
        }
    }
    stats.values_initialized_from_A           = true;
    stats.values_initialized_from_csc_numeric = false;

    // ------------------------------------------------------------------
    // Step 2: Interleaved §17.2(A) + §17.2(B) per supernode j
    // ------------------------------------------------------------------
    // SLU-SN-OPT: one reusable scratch shared across all (k,j) panel updates,
    // eliminating per-call heap allocation in apply_supernode_panel_update.
    // panel_work is also hoisted so its row_pos map (sized to n) is reused across
    // supernodes instead of reallocated per panel.
    sparse_lu_detail::panel_apply_scratch<T, Index> apply_scr;
    sparse_lu_detail::supernode_panel_workspace<T, Index> panel_work;
    for (std::size_t j = 0u; j < nsup; ++j) {
        const Index j_idx = static_cast<Index>(j);

        // §17.2(A): panel update from completed k < j
        const auto t_sym0 = std::chrono::steady_clock::now();
        sparse_lu_detail::supernode_update_set<Index> uset =
            sparse_lu_detail::compute_panel_update_set(storage, j_idx, col_to_supernode);
        stats.symbolic_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_sym0).count());

        if (!uset.updating_supernodes.empty()) {
            const auto t_pa0 = std::chrono::steady_clock::now();
            sparse_lu_detail::supernode_panel_workspace<T, Index>& work = panel_work;
            sparse_lu_detail::gather_panel_workspace(storage, j_idx, work);

            if (!work.valid) {
                stats.status = supernodal_true_numeric_status::failed;
                finish_total_ticks(); // SLU-8R.6.2: was missing — all return paths must set total_ticks
                return stats;
            }

            sparse_lu_detail::supernode_panel_update_stats pu_stats;
            for (std::size_t ki = 0u;
                 ki < uset.updating_supernodes.size(); ++ki) {
                sparse_lu_detail::apply_supernode_panel_update(
                    storage,
                    uset.updating_supernodes[ki],
                    j_idx,
                    work,
                    pu_stats,
                    apply_scr);
            }
            sparse_lu_detail::scatter_panel_workspace(
                storage, j_idx, work, pu_stats);

            // SLU-SN-OPT: §17.2(A) non-kernel time = block wall-clock minus the
            // dense-kernel adapter time accumulated inside the block.
            {
                const std::size_t pa_total = static_cast<std::size_t>(
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - t_pa0).count());
                stats.panel_nonkernel_ticks +=
                    (pa_total > pu_stats.dense_kernel_ticks)
                        ? (pa_total - pu_stats.dense_kernel_ticks) : 0u;
            }

            stats.panel_update_count              += pu_stats.update_count;
            stats.trsm_count                      += pu_stats.trsm_count;
            stats.gemm_count                      += pu_stats.gemm_count;
            stats.gemv_count                      += pu_stats.gemv_count;
            // SLU-8R.6.1: accumulate §17.2(A) dense kernel timing into true-numeric stats.
            stats.panel_update_dense_kernel_ticks += pu_stats.dense_kernel_ticks;
            stats.dense_kernel_ticks              += pu_stats.dense_kernel_ticks;
            // SLU-PERF: accumulate §17.2(A) dense-kernel FLOP + gemm shape.
            stats.flop_trsm  += pu_stats.flop_trsm;
            stats.flop_gemm  += pu_stats.flop_gemm;
            stats.flop_gemv  += pu_stats.flop_gemv;
            stats.gemm_m_sum += pu_stats.gemm_m_sum;
            stats.gemm_n_sum += pu_stats.gemm_n_sum;
            stats.gemm_k_sum += pu_stats.gemm_k_sum;
            if (pu_stats.gemm_max_m > stats.gemm_max_m) stats.gemm_max_m = pu_stats.gemm_max_m;
            if (pu_stats.gemm_max_n > stats.gemm_max_n) stats.gemm_max_n = pu_stats.gemm_max_n;
            if (pu_stats.gemm_max_k > stats.gemm_max_k) stats.gemm_max_k = pu_stats.gemm_max_k;
            for (int hi = 0; hi < 8; ++hi)
                stats.gemm_dim_hist[hi] += pu_stats.gemm_dim_hist[hi];
        }

        // §17.2(B): within-panel factorization for supernode j
        const auto t_wp0 = std::chrono::steady_clock::now();
        sparse_lu_detail::within_panel_factor_stats wpf_j;
        bool completed_j = sparse_lu_detail::factorize_within_panel_single(
            storage, j_idx, opt, wpf_j);
        {
            const std::size_t wp_total = static_cast<std::size_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - t_wp0).count());
            stats.within_nonkernel_ticks +=
                (wp_total > wpf_j.dense_kernel_ticks)
                    ? (wp_total - wpf_j.dense_kernel_ticks) : 0u;
        }

        stats.ger_count += wpf_j.ger_count;
        // SLU-8R.6.1: accumulate §17.2(B) dense kernel timing into true-numeric stats.
        stats.within_panel_dense_kernel_ticks += wpf_j.dense_kernel_ticks;
        stats.dense_kernel_ticks              += wpf_j.dense_kernel_ticks;
        // SLU-PERF: accumulate §17.2(B) ger FLOP.
        stats.flop_ger += wpf_j.flop_ger;
        stats.within_panel_count++;

        if (!completed_j) {
            stats.supernodes_processed = j + 1u;
            stats.status = supernodal_true_numeric_status::pivot_failure;
            // Leave storage without true_numeric_source on pivot failure.
            finish_total_ticks();
            return stats;
        }
    }

    stats.supernodes_processed = nsup;

    // ------------------------------------------------------------------
    // Step 3: Factorization residual check (SLU-8R.8 B2)
    //
    // Scalable sparse/column-accumulator residual checking. This path is NEVER
    // skipped (no n_max cutoff): large n is checked just like small n. An
    // accepted true-numeric source therefore ALWAYS has residual_checked == true
    // AND residual_passed == true. There is no skip-accept.
    //
    // Residual type: point-arithmetic Frobenius sanity check.
    // ------------------------------------------------------------------
    {
        double abs_res = 0.0, rel_res = 0.0;
        const auto t_res0 = std::chrono::steady_clock::now();
        bool checked = sparse_lu_detail::compute_supernodal_lu_residual_sparse(
            A_csc, storage, n, abs_res, rel_res);
        stats.residual_ticks += static_cast<std::size_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t_res0).count());

        stats.factorization_residual_abs = abs_res;
        stats.factorization_residual_rel = rel_res;

        if (!checked) {
            // Sparse residual could not run (degenerate n / malformed CSC).
            // Per SLU-8R.8: an unchecked residual MUST NOT be accepted as a
            // residual-checked true-numeric source. true_numeric_source stays
            // false; CSC fallback remains in effect.
            stats.factorization_residual_checked = false;
            stats.factorization_residual_passed  = false;
            stats.status =
                supernodal_true_numeric_status::residual_not_checked;
            finish_total_ticks();
            return stats;
        }

        stats.factorization_residual_checked = true;
        // Point-arithmetic sanity threshold for well-conditioned deterministic
        // matrices: 1e-6 relative (or absolute). This is a residual sanity gate,
        // NOT a rigorous error bound.
        const double tol = 1e-6;
        stats.factorization_residual_passed =
            (rel_res <= tol || abs_res <= tol);

        if (!stats.factorization_residual_passed) {
            stats.status =
                supernodal_true_numeric_status::residual_failed;
            // Do NOT set true_numeric_source on residual failure.
            finish_total_ticks();
            return stats;
        }
    }

    // ------------------------------------------------------------------
    // Step 4: Accept true-numeric source
    // ------------------------------------------------------------------
    storage.true_numeric_source            = true;
    storage.values_initialized_from_A      = true;
    storage.values_initialized_from_csc_numeric = false;
    storage.numeric_source_kind =
        supernodal_numeric_source_kind::a_eff_true_numeric;

    stats.success = true;
    stats.status  = supernodal_true_numeric_status::success;
    finish_total_ticks();
    return stats;
}

#endif // VCP_TSPARSE_SPARSE_LU_TRUE_NUMERIC_IMPL_HPP
