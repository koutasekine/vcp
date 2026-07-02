// vcp/bfem/dirichlet.hpp
// Layer 3: homogeneous Dirichlet reduction (H10, X3) -- fully independent of
// the assembly: matrices/vectors + boundary dof list in, reduced system out.
//
// Conforms to: L3 external design v0.2 (section 6) and
//              L3 internal design v0.2 (section 8).
//
// reduce(A) uses the MAIN plan of the reconciliation (component enumeration
// via spm_adapter::for_each_entry; see fe_space.hpp header note): rows and
// columns of constrained dofs are removed (homogeneous condition, symmetry
// preserved), no right-hand-side lifting occurs.

#ifndef VCP_BFEM_DIRICHLET_HPP
#define VCP_BFEM_DIRICHLET_HPP

#include <vector>
#include <algorithm>
#include <utility>
#include <stdexcept>

#include <vcp/matrix.hpp>

#include <vcp/bfem/fe_space.hpp>   // detail::coo_buffer / detail::spm_adapter

namespace vcp {
namespace bfem {

// SP: sparse policy, forwarded into spmatrix_t (default: the approximate
// spmats<T>; see the policy note at the top of fe_space.hpp)
template <typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class dirichlet_reduction {
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;

    // bdofs: constrained global dofs; duplicates or out-of-range entries are
    // an error (duplicates are NOT silently removed -- external design 8)
    dirichlet_reduction(int full_size, std::vector<int> bdofs)
        : full_(full_size), to_reduced_(), to_full_() {
        if (full_size < 0)
            throw std::invalid_argument("bfem::dirichlet_reduction: negative size");
        std::sort(bdofs.begin(), bdofs.end());
        for (std::size_t k = 0; k < bdofs.size(); ++k) {
            if (bdofs[k] < 0 || bdofs[k] >= full_size)
                throw std::invalid_argument("bfem::dirichlet_reduction: bdof out of range");
            if (k > 0 && bdofs[k] == bdofs[k - 1])
                throw std::invalid_argument("bfem::dirichlet_reduction: duplicate bdof");
        }
        to_reduced_.assign(static_cast<std::size_t>(full_size), 0);
        for (std::size_t k = 0; k < bdofs.size(); ++k)
            to_reduced_[static_cast<std::size_t>(bdofs[k])] = -1;
        to_full_.reserve(static_cast<std::size_t>(full_size) - bdofs.size());
        int r = 0;
        for (int i = 0; i < full_size; ++i) {
            if (to_reduced_[static_cast<std::size_t>(i)] == -1) continue;
            to_reduced_[static_cast<std::size_t>(i)] = r;
            to_full_.push_back(i);
            ++r;
        }
    }

    int full_size() const { return full_; }
    int reduced_size() const { return static_cast<int>(to_full_.size()); }
    int to_reduced(int full_dof) const {          // constrained dofs map to -1
        if (full_dof < 0 || full_dof >= full_)
            throw std::invalid_argument("bfem::dirichlet_reduction::to_reduced: out of range");
        return to_reduced_[static_cast<std::size_t>(full_dof)];
    }
    int to_full(int reduced_dof) const {
        if (reduced_dof < 0 || reduced_dof >= reduced_size())
            throw std::invalid_argument("bfem::dirichlet_reduction::to_full: out of range");
        return to_full_[static_cast<std::size_t>(reduced_dof)];
    }

    // remove constrained rows AND columns (u = 0; keeps symmetry / positive
    // definiteness when present)
    spmatrix_t reduce(const spmatrix_t& A) const {
        if (A.rowsize() != full_ || A.columnsize() != full_)
            throw std::invalid_argument("bfem::dirichlet_reduction::reduce: size mismatch");
        detail::coo_buffer<T> buf;
        collector col = { this, &buf };
        detail::spm_adapter<T, SP>::for_each_entry(A, col);
        buf.combine();
        return detail::spm_adapter<T, SP>::build(reduced_size(), reduced_size(), buf);
    }

    vcp::matrix<T, P> reduce(const vcp::matrix<T, P>& v) const {
        if (v.rowsize() != full_ || v.columnsize() != 1)
            throw std::invalid_argument("bfem::dirichlet_reduction::reduce: vector size mismatch");
        vcp::matrix<T, P> out;
        out.zeros(reduced_size(), 1);
        for (int r = 0; r < reduced_size(); ++r)
            out(r, 0) = v(to_full_[static_cast<std::size_t>(r)], 0);
        return out;
    }

    // reinsert T(0) at the constrained dofs
    vcp::matrix<T, P> expand(const vcp::matrix<T, P>& v_reduced) const {
        if (v_reduced.rowsize() != reduced_size() || v_reduced.columnsize() != 1)
            throw std::invalid_argument("bfem::dirichlet_reduction::expand: vector size mismatch");
        vcp::matrix<T, P> out;
        out.zeros(full_, 1);
        for (int r = 0; r < reduced_size(); ++r)
            out(to_full_[static_cast<std::size_t>(r)], 0) = v_reduced(r, 0);
        return out;
    }

private:
    int full_;
    std::vector<int> to_reduced_;
    std::vector<int> to_full_;

    struct collector {
        const dirichlet_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            int ri = self->to_reduced_[static_cast<std::size_t>(i)];
            int rj = self->to_reduced_[static_cast<std::size_t>(j)];
            if (ri >= 0 && rj >= 0) buf->push(ri, rj, v);
        }
    };
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_DIRICHLET_HPP
