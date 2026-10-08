// vcp/btfem/kron.hpp
// btfem T5 / T11: Kronecker factors kron_term<T,P,SP>, the expansion helpers
// kron_expand (single term / sum of terms), the matrix-free action kron_apply,
// the block-wise Dirichlet reduction block_reduce / block_expand, and the time
// initial-condition reduction initial_reduction<T,P>.
//
// Conforms to: BT-1 external design v1.1 (sections 7.1, 7.2, 7.3) and
//              BT-1 internal design v1.0 (sections 7.1-7.4).
//
// spmatrix is used for data holding and as a return type only; all component
// access goes through bfem detail::spm_adapter::for_each_entry and all
// accumulation through detail::coo_buffer::combine (D-1, D-10). Time-block
// major indexing: row a n_x + i, column b n_x' + j (R12).
// Accumulation order (external 7.1): term order -> time (a, b) row major ->
// space enumeration order; duplicates combined by the deterministic
// coo_buffer::combine. Exact zeros of the time factor (tv == T(0); for
// interval T only the degenerate zero) are skipped.

#ifndef VCP_BTFEM_KRON_HPP
#define VCP_BTFEM_KRON_HPP

#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/fe_space.hpp>      // detail::coo_buffer / detail::spm_adapter
#include <vcp/bfem/fe_function.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/btfem/time_space.hpp>

namespace vcp {
namespace btfem {

// ---------------------------------------------------------------------------
// kron_term: scale * (time (x) space)
// ---------------------------------------------------------------------------
template <typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
struct kron_term {
    T scale;                                    // scalar factor (nu etc.), default T(1)
    vcp::matrix<T, P> time;                     // ndof_t x ndof_t' (dense)
    vcp::spmatrix<T, SP> space;                 // n_x x n_x'

    kron_term() : scale(1), time(), space() {}
    kron_term(const T& s, vcp::matrix<T, P> t, vcp::spmatrix<T, SP> sp)
        : scale(s), time(std::move(t)), space(std::move(sp)) {}

    int time_rows() const { return time.rowsize(); }
    int time_cols() const { return time.columnsize(); }
    int space_rows() const { return static_cast<int>(space.rowsize()); }
    int space_cols() const { return static_cast<int>(space.columnsize()); }
};

namespace detail {

template <typename T>
struct kron_pusher {
    bfem::detail::coo_buffer<T>* buf;
    int row_off, col_off;
    const T* s;
    template <typename I>
    void operator()(I i, I j, const T& v) const {
        buf->push(row_off + static_cast<int>(i), col_off + static_cast<int>(j), (*s) * v);
    }
};

template <typename T>
struct kron_counter {
    std::size_t* n;
    template <typename I>
    void operator()(I, I, const T&) const { ++(*n); }
};

template <typename T, typename P, class SP>
void kron_check_shapes(const std::vector<kron_term<T, P, SP> >& terms, const char* where) {
    if (terms.empty()) {
        std::string msg("btfem::");
        msg += where;
        msg += ": empty term list";
        throw std::invalid_argument(msg);
    }
    for (std::size_t t = 1; t < terms.size(); ++t) {
        if (terms[t].time_rows() != terms[0].time_rows() || terms[t].time_cols() != terms[0].time_cols()
            || terms[t].space_rows() != terms[0].space_rows()
            || terms[t].space_cols() != terms[0].space_cols()) {
            std::string msg("btfem::");
            msg += where;
            msg += ": term shapes differ";
            throw std::invalid_argument(msg);
        }
    }
}

template <typename T, typename P, class SP>
void kron_expand_into(bfem::detail::coo_buffer<T>& buf, const kron_term<T, P, SP>& t) {
    const int nx = t.space_rows();
    const int nxp = t.space_cols();
    const T zero(0);
    for (int a = 0; a < t.time_rows(); ++a) {
        for (int b = 0; b < t.time_cols(); ++b) {
            const T& tv = t.time(a, b);
            if (tv == zero) continue;                         // exact zero only
            const T s = t.scale * tv;
            kron_pusher<T> push = { &buf, a * nx, b * nxp, &s };
            bfem::detail::spm_adapter<T, SP>::for_each_entry(t.space, push);
        }
    }
}

template <typename T, typename P, class SP>
std::size_t kron_nnz_estimate(const kron_term<T, P, SP>& t) {
    std::size_t ns = 0;
    kron_counter<T> cnt = { &ns };
    bfem::detail::spm_adapter<T, SP>::for_each_entry(t.space, cnt);
    std::size_t nt = 0;
    const T zero(0);
    for (int a = 0; a < t.time_rows(); ++a)
        for (int b = 0; b < t.time_cols(); ++b)
            if (!(t.time(a, b) == zero)) ++nt;
    return ns * nt;
}

} // namespace detail

// ---------------------------------------------------------------------------
// kron_expand: time-block major expansion to a (rows_t n_x) x (cols_t n_x')
// spmatrix (triplets inserted, finalize not called -- bfem X5)
// ---------------------------------------------------------------------------
template <typename T, typename P, class SP>
vcp::spmatrix<T, SP> kron_expand(const std::vector<kron_term<T, P, SP> >& terms) {
    detail::kron_check_shapes(terms, "kron_expand");
    bfem::detail::coo_buffer<T> buf;
    std::size_t est = 0;
    for (std::size_t t = 0; t < terms.size(); ++t) est += detail::kron_nnz_estimate(terms[t]);
    buf.reserve(est);
    for (std::size_t t = 0; t < terms.size(); ++t) detail::kron_expand_into(buf, terms[t]);
    buf.combine();
    const kron_term<T, P, SP>& t0 = terms[0];
    return bfem::detail::spm_adapter<T, SP>::build(t0.time_rows() * t0.space_rows(),
                                                   t0.time_cols() * t0.space_cols(), buf);
}

template <typename T, typename P, class SP>
vcp::spmatrix<T, SP> kron_expand(const kron_term<T, P, SP>& t) {
    std::vector<kron_term<T, P, SP> > one(1, t);
    return kron_expand(one);
}

// ---------------------------------------------------------------------------
// kron_apply: y = sum_t scale_t (time_t (x) space_t) x without expansion.
// Order: term -> time column b -> space enumeration -> time row a (internal 7.2).
// ---------------------------------------------------------------------------
namespace detail {
template <typename T>
struct kron_spmv {
    std::vector<T>* z;
    const T* xcol;                       // x block (pointer to x(b n_x', 0))
    template <typename I>
    void operator()(I i, I j, const T& v) const {
        (*z)[static_cast<std::size_t>(i)] += v * xcol[static_cast<std::size_t>(j)];
    }
};
} // namespace detail

template <typename T, typename P, class SP>
vcp::matrix<T, P> kron_apply(const std::vector<kron_term<T, P, SP> >& terms,
                             const vcp::matrix<T, P>& x) {
    detail::kron_check_shapes(terms, "kron_apply");
    const kron_term<T, P, SP>& t0 = terms[0];
    const int nx = t0.space_rows(), nxp = t0.space_cols();
    const int rt = t0.time_rows(), ct = t0.time_cols();
    if (x.rowsize() != ct * nxp || x.columnsize() != 1)
        throw std::invalid_argument("btfem::kron_apply: x size != cols_t * n_x' x 1");
    vcp::matrix<T, P> y;
    y.zeros(rt * nx, 1);
    std::vector<T> z(static_cast<std::size_t>(nx));
    std::vector<T> xb(static_cast<std::size_t>(nxp));
    const T zero(0);
    for (std::size_t t = 0; t < terms.size(); ++t) {
        const kron_term<T, P, SP>& k = terms[t];
        for (int b = 0; b < ct; ++b) {
            for (int j = 0; j < nxp; ++j) xb[static_cast<std::size_t>(j)] = x(b * nxp + j, 0);
            for (int i = 0; i < nx; ++i) z[static_cast<std::size_t>(i)] = zero;
            detail::kron_spmv<T> mv = { &z, xb.data() };
            bfem::detail::spm_adapter<T, SP>::for_each_entry(k.space, mv);
            for (int a = 0; a < rt; ++a) {
                const T& tv = k.time(a, b);
                if (tv == zero) continue;
                const T s = k.scale * tv;
                for (int i = 0; i < nx; ++i) y(a * nx + i, 0) += s * z[static_cast<std::size_t>(i)];
            }
        }
    }
    return y;
}

// ---------------------------------------------------------------------------
// block-wise spatial Dirichlet reduction (external 7.2)
// ---------------------------------------------------------------------------
template <typename T, typename P, class SP>
vcp::matrix<T, P> block_reduce(const bfem::dirichlet_reduction<T, P, SP>& R, int nblocks,
                               const vcp::matrix<T, P>& v) {
    const int nf = R.full_size(), nr = R.reduced_size();
    if (nblocks < 0 || v.rowsize() != nblocks * nf || v.columnsize() != 1)
        throw std::invalid_argument("btfem::block_reduce: vector size != nblocks * full_size x 1");
    vcp::matrix<T, P> out;
    out.zeros(nblocks * nr, 1);
    for (int a = 0; a < nblocks; ++a)
        for (int r = 0; r < nr; ++r)
            out(a * nr + r, 0) = v(a * nf + R.to_full(r), 0);
    return out;
}

template <typename T, typename P, class SP>
vcp::matrix<T, P> block_expand(const bfem::dirichlet_reduction<T, P, SP>& R, int nblocks,
                               const vcp::matrix<T, P>& vr) {
    const int nf = R.full_size(), nr = R.reduced_size();
    if (nblocks < 0 || vr.rowsize() != nblocks * nr || vr.columnsize() != 1)
        throw std::invalid_argument("btfem::block_expand: vector size != nblocks * reduced_size x 1");
    vcp::matrix<T, P> out;
    out.zeros(nblocks * nf, 1);
    for (int a = 0; a < nblocks; ++a)
        for (int r = 0; r < nr; ++r)
            out(a * nf + R.to_full(r), 0) = vr(a * nr + r, 0);
    return out;
}

namespace detail {
template <typename T, typename P, class SP>
struct block_collector {
    const bfem::dirichlet_reduction<T, P, SP>* R;
    bfem::detail::coo_buffer<T>* buf;
    int nf, nr;
    template <typename I>
    void operator()(I i, I j, const T& v) const {
        const int ii = static_cast<int>(i), jj = static_cast<int>(j);
        const int a = ii / nf, ia = ii % nf;
        const int b = jj / nf, jb = jj % nf;
        const int ri = R->to_reduced(ia), rj = R->to_reduced(jb);
        if (ri >= 0 && rj >= 0) buf->push(a * nr + ri, b * nr + rj, v);
    }
};
} // namespace detail

template <typename T, typename P, class SP>
vcp::spmatrix<T, SP> block_reduce(const bfem::dirichlet_reduction<T, P, SP>& R, int nblocks,
                                  const vcp::spmatrix<T, SP>& A) {
    const int nf = R.full_size(), nr = R.reduced_size();
    if (nblocks < 0 || static_cast<int>(A.rowsize()) != nblocks * nf
        || static_cast<int>(A.columnsize()) != nblocks * nf)
        throw std::invalid_argument("btfem::block_reduce: matrix size != (nblocks * full_size)^2");
    bfem::detail::coo_buffer<T> buf;
    detail::block_collector<T, P, SP> col = { &R, &buf, nf, nr };
    bfem::detail::spm_adapter<T, SP>::for_each_entry(A, col);
    buf.combine();
    return bfem::detail::spm_adapter<T, SP>::build(nblocks * nr, nblocks * nr, buf);
}

// ---------------------------------------------------------------------------
// initial_reduction<T, P> (T11): Vbar^q_k (bc_none, m q + 1 time dofs) <->
// V^q_k (bc_zero_left, m q time dofs) on the same mesh; the left-end block
// (time dof 0, n_x components) is removed / reinserted from u_init.
// ---------------------------------------------------------------------------
template <typename T, typename P = vcp::mats<T> >
class initial_reduction {
public:
    initial_reduction(const time_space<T>& full, const time_space<T>& reduced, int n_x)
        : mt_full_(full.ndof()), mt_red_(reduced.ndof()), nx_(n_x) {
        if (full.degree() != reduced.degree())
            throw std::invalid_argument("btfem::initial_reduction: degrees differ");
        if (full.bc() != time_bc_none || reduced.bc() != time_bc_zero_left)
            throw std::invalid_argument("btfem::initial_reduction: full must be bc_none, reduced bc_zero_left");
        if (!full.same_mesh(reduced))
            throw std::invalid_argument("btfem::initial_reduction: different time meshes");
        if (n_x < 1)
            throw std::invalid_argument("btfem::initial_reduction: n_x < 1");
        if (mt_full_ != mt_red_ + 1)
            throw std::invalid_argument("btfem::initial_reduction: inconsistent dof counts");
    }

    int full_size() const { return mt_full_ * nx_; }
    int reduced_size() const { return mt_red_ * nx_; }
    int space_size() const { return nx_; }

    vcp::matrix<T, P> reduce(const vcp::matrix<T, P>& v_full) const {
        if (v_full.rowsize() != full_size() || v_full.columnsize() != 1)
            throw std::invalid_argument("btfem::initial_reduction::reduce: size mismatch");
        vcp::matrix<T, P> out;
        out.zeros(reduced_size(), 1);
        for (int i = 0; i < reduced_size(); ++i) out(i, 0) = v_full(nx_ + i, 0);
        return out;
    }
    template <int D>
    vcp::matrix<T, P> expand(const vcp::matrix<T, P>& v_reduced,
                             const bfem::fe_function<D, T, P>& u_init) const {
        if (v_reduced.rowsize() != reduced_size() || v_reduced.columnsize() != 1)
            throw std::invalid_argument("btfem::initial_reduction::expand: size mismatch");
        if (u_init.coeffs().rowsize() != nx_ || u_init.coeffs().columnsize() != 1)
            throw std::invalid_argument("btfem::initial_reduction::expand: u_init size != n_x");
        vcp::matrix<T, P> out;
        out.zeros(full_size(), 1);
        for (int i = 0; i < nx_; ++i) out(i, 0) = u_init.coeffs()(i, 0);
        for (int i = 0; i < reduced_size(); ++i) out(nx_ + i, 0) = v_reduced(i, 0);
        return out;
    }

private:
    int mt_full_, mt_red_, nx_;
};

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_KRON_HPP
