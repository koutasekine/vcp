// VCP Library
// http://verified.computation.jp
//
// VCP Library is licensed under the BSD 3-clause "New" or "Revised" License

// SLDL-SH SH-0 -- shift assembly plan for the A - sigma*B iteration (B2).
//
// sparse_ldl_shift_plan<T, Index> precomputes, ONCE per (pattern(A),
// pattern(B)) pair, everything the sigma loop of ldl_shift_handle needs to
// assemble the LOWER-triangle CSC values of A - sigma*B in a single O(nnz)
// pass (design SLDL-SH v1 SS2.1):
//
//   - the merged lower-triangle CSC pattern  pattern(A) U pattern(B)
//     (rows ascending per column, no duplicates).  When B is omitted the
//     shift target is A - sigma*I and the merged pattern is
//     pattern(A) U diagonal: a column whose diagonal is structurally absent
//     from A gets it added, so the -sigma contribution always has a slot.
//   - index maps mapA / mapB ALIGNED WITH THE INPUT val ARRAYS: entry k of
//     the input CSC lands at merged position mapA[k] (resp. mapB[k]).
//     A strictly-upper stored entry maps to -1 and is ignored, mirroring the
//     factorization's lower-triangle-only read (ldl design v2 SS1.2).
//   - diag_pos[j] (B omitted only): the merged position of (j, j).
//
// build_values(sigma, valA, valB, out) then performs the P2-compliant
// T-generic one-pass assembly
//
//     out = 0;  out[mapA[k]] += valA[k];  out[mapB[k]] -= sigma * valB[k];
//     (B omitted:                         out[diag_pos[j]] -= sigma;)
//
// out is (re)sized with assign(), so a caller that reuses the same vector
// performs NO dynamic allocation from the second call on (STOP-1 ruling A-1).
// Cancellation may leave exact zeros in out; that is intentional -- the
// numeric kernels accept explicit zeros and the certified-zero drop happens
// at their output boundary, exactly as in the one-shot path.
//
// Unlike the other tsparse_sparse_ldl_*_impl.hpp files this header is NOT
// injected by tsparse_sparse_ldl.hpp (which is frozen for this track): it is
// SELF-CONTAINED and is included at file scope by the B2 policy layer
// (vcp/spmats_base/spmats_ldl_shift.hpp).  Including it directly is fine.
//
// Scalar contract: module scalar contract of tsparse_scalar.hpp (SLU-GT1 D8)
// -- the value pass uses T arithmetic only (no comparisons, no abs); the
// pattern part is T-independent integer code (P4/P6-safe).

#ifndef VCP_TSPARSE_SPARSE_LDL_SHIFT_PLAN_IMPL_HPP
#define VCP_TSPARSE_SPARSE_LDL_SHIFT_PLAN_IMPL_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

#include <vcp/tsparse/tsparse_sparse_ldl.hpp>

namespace vcp {

enum class sparse_ldl_shift_plan_status {
    not_built,      // default-constructed; no pattern attached
    success,
    invalid_input   // malformed CSC pattern (either input)
};

inline const char* sparse_ldl_shift_plan_status_to_string(
    sparse_ldl_shift_plan_status s)
{
    switch (s) {
    case sparse_ldl_shift_plan_status::not_built:     return "not_built";
    case sparse_ldl_shift_plan_status::success:       return "success";
    case sparse_ldl_shift_plan_status::invalid_input: return "invalid_input";
    }
    return "unknown";
}

template <class T, class Index>
struct sparse_ldl_shift_plan {
    static_assert(std::is_signed<Index>::value,
                  "sparse LDL Index must be signed");

    // Every field is filled by the constructor and never modified afterwards
    // (read-only after construction; the B2 handle shares plans as const).
    sparse_ldl_shift_plan_status status;
    Index n;
    bool has_B;

    // Merged lower-triangle CSC pattern (rows ascending per column).
    std::vector<Index> col_ptr, row_ind;

    // Maps aligned with the INPUT val arrays; -1 = strictly-upper stored
    // entry, ignored by the assembly (never a valid merged position).
    std::vector<Index> mapA;        // size = nnz of the A input
    std::vector<Index> mapB;        // size = nnz of the B input; empty when !has_B
    std::vector<Index> diag_pos;    // size n when !has_B: merged position of (j,j)

    sparse_ldl_shift_plan()
        : status(sparse_ldl_shift_plan_status::not_built),
          n(Index(0)), has_B(false) {}

    // A - sigma*I (B omitted)
    sparse_ldl_shift_plan(
        const Index n_in,
        const std::vector<Index>& colA,
        const std::vector<Index>& rowA)
        : status(sparse_ldl_shift_plan_status::not_built),
          n(n_in), has_B(false)
    {
        build_(colA, rowA, static_cast<const std::vector<Index>*>(0),
               static_cast<const std::vector<Index>*>(0));
    }

    // A - sigma*B
    sparse_ldl_shift_plan(
        const Index n_in,
        const std::vector<Index>& colA,
        const std::vector<Index>& rowA,
        const std::vector<Index>& colB,
        const std::vector<Index>& rowB)
        : status(sparse_ldl_shift_plan_status::not_built),
          n(n_in), has_B(true)
    {
        build_(colA, rowA, &colB, &rowB);
    }

    Index nnz() const { return static_cast<Index>(row_ind.size()); }

    // -----------------------------------------------------------------------
    // build_values -- the sigma-dependent O(nnz) assembly pass (T-generic).
    // Returns false (out contents unspecified) when the plan is not built or
    // an input val array does not match the pattern the plan was built from;
    // never throws (P3).  valB must be empty when the plan has no B.
    // -----------------------------------------------------------------------
    bool build_values(const T& sigma,
                      const std::vector<T>& valA,
                      const std::vector<T>& valB,
                      std::vector<T>& out) const
    {
        if (status != sparse_ldl_shift_plan_status::success) return false;
        if (valA.size() != mapA.size()) return false;
        if (has_B) { if (valB.size() != mapB.size()) return false; }
        else       { if (!valB.empty()) return false; }

        out.assign(row_ind.size(), T(0));
        for (std::size_t k = 0; k < mapA.size(); ++k) {
            const Index m = mapA[k];
            if (m >= Index(0)) out[static_cast<std::size_t>(m)] += valA[k];
        }
        if (has_B) {
            for (std::size_t k = 0; k < mapB.size(); ++k) {
                const Index m = mapB[k];
                if (m >= Index(0)) out[static_cast<std::size_t>(m)] -= sigma * valB[k];
            }
        } else {
            for (std::size_t j = 0; j < diag_pos.size(); ++j) {
                out[static_cast<std::size_t>(diag_pos[j])] -= sigma;
            }
        }
        return true;
    }

    // B-omitted convenience overload (A - sigma*I).
    bool build_values(const T& sigma,
                      const std::vector<T>& valA,
                      std::vector<T>& out) const
    {
        const std::vector<T> empty_valB;
        return build_values(sigma, valA, empty_valB, out);
    }

private:
    // Pattern merge (setup-time; allocation here is fine -- the sigma loop
    // never re-enters this).  colB/rowB are null when B is omitted.
    void build_(const std::vector<Index>& colA,
                const std::vector<Index>& rowA,
                const std::vector<Index>* colB,
                const std::vector<Index>* rowB)
    {
        if (!sparse_ldl_detail::sparse_ldl_validate_pattern_(n, colA, rowA)) {
            status = sparse_ldl_shift_plan_status::invalid_input;
            return;
        }
        if (has_B &&
            !sparse_ldl_detail::sparse_ldl_validate_pattern_(n, *colB, *rowB)) {
            status = sparse_ldl_shift_plan_status::invalid_input;
            return;
        }

        const std::size_t un = static_cast<std::size_t>(n);
        col_ptr.assign(un + 1u, Index(0));
        row_ind.clear();
        row_ind.reserve(rowA.size() + (has_B ? rowB->size() : un));
        mapA.assign(rowA.size(), Index(-1));
        if (has_B) mapB.assign(rowB->size(), Index(-1));
        else       mapB.clear();
        if (!has_B) diag_pos.assign(un, Index(-1));
        else        diag_pos.clear();

        for (std::size_t j = 0; j < un; ++j) {
            const Index cj = static_cast<Index>(j);

            // Cursor over the lower entries of A(:,j); strictly-upper stored
            // entries are consumed here (map -1) so both cursors only ever
            // look at lower rows during the merge.
            Index ka = colA[j];
            const Index ea = colA[j + 1u];
            while (ka < ea && rowA[static_cast<std::size_t>(ka)] < cj) {
                mapA[static_cast<std::size_t>(ka)] = Index(-1);
                ++ka;
            }
            // Second source: lower entries of B(:,j), or the singleton {j}.
            Index kb = has_B ? (*colB)[j] : Index(0);
            const Index eb = has_B ? (*colB)[j + 1u] : Index(1);
            if (has_B) {
                while (kb < eb && (*rowB)[static_cast<std::size_t>(kb)] < cj) {
                    mapB[static_cast<std::size_t>(kb)] = Index(-1);
                    ++kb;
                }
            }

            while (ka < ea || kb < eb) {
                const Index ra = (ka < ea) ? rowA[static_cast<std::size_t>(ka)] : n;
                const Index rb = (kb < eb)
                               ? (has_B ? (*rowB)[static_cast<std::size_t>(kb)] : cj)
                               : n;
                const Index r = (ra < rb) ? ra : rb;
                const Index pos = static_cast<Index>(row_ind.size());
                row_ind.push_back(r);
                if (ra == r) { mapA[static_cast<std::size_t>(ka)] = pos; ++ka; }
                if (rb == r) {
                    if (has_B) { mapB[static_cast<std::size_t>(kb)] = pos; }
                    ++kb;
                }
                if (!has_B && r == cj) diag_pos[j] = pos;
            }
            col_ptr[j + 1u] = static_cast<Index>(row_ind.size());
        }

        status = sparse_ldl_shift_plan_status::success;
    }
};

} // namespace vcp

#endif // VCP_TSPARSE_SPARSE_LDL_SHIFT_PLAN_IMPL_HPP
