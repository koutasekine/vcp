// vcp/bfem/rt/broken_space.hpp
// RT Layer 3: broken_dofmap (detail), broken_field and
// broken_space<D,T,P,SP> (K2 -- Q3).
//
// Conforms to: RT-L3 external design v0.2 (sections 3.1, 4) and
//              RT-L3 internal design v0.2 (section 3).
//
// Element-wise discontinuous P_l: the dof map is the block arithmetic
// gdof(e, r) = e N(D,l) + r ONLY (no edge identification, no signs; the
// broken_dofmap holds no array). Broken fields are barycentric-defined
// (W-RT9). The mass matrix is block diagonal (|T| M^{(l,l)} per element);
// entries outside the blocks are never emitted (RE-5 / RTX-15 checks the
// ABSENCE of triplets, not zero values).

#ifndef VCP_BFEM_RT_BROKEN_SPACE_HPP
#define VCP_BFEM_RT_BROKEN_SPACE_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/poly_field.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/fe_space.hpp>          // detail::coo_buffer / spm_adapter
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/element_op.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/detail/scalar_traits.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class broken_space;

namespace detail {

// pure arithmetic, no storage (internal design section 3); provides the
// same surface as the other dofmaps so the generalized Y2 kernels accept it
class broken_dofmap {
public:
    typedef pn_family_tag family_tag;             // identity (signs are +1)

    broken_dofmap() : nloc_(0), nt_(0) {}
    broken_dofmap(int nloc, int nt) : nloc_(nloc), nt_(nt) {}

    int local_size() const { return nloc_; }
    int num_elements() const { return nt_; }
    int ndof() const { return nloc_ * nt_; }
    int global_dof(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return e * nloc_ + r;
    }
    int dof_sign(int, int) const { return 1; }

private:
    int nloc_, nt_;
};

} // namespace detail

// ---------------------------------------------------------------------------
// broken_field<D,T,P> (B-2): degree l + global coefficients (element blocks,
// each block in the canonical L0 order, barycentric definition).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class broken_field {
public:
    int order() const { return l_; }
    const vcp::matrix<T, P>& coeffs() const { return c_; }
    vcp::matrix<T, P>&       coeffs() { return c_; }

private:
    explicit broken_field(int l) : l_(l), c_() {}
    int l_;
    vcp::matrix<T, P> c_;
    template <int DD, typename TT, typename PP, class SS> friend class broken_space;
};

// ---------------------------------------------------------------------------
// broken_space<D,T,P,SP> (K2)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class broken_space {
    // phase 5c (D5C-5): broken_space is dimension uniform; D == 3 is enabled
    // by relaxing this guard only (no other change -- see
    // sandbox/docs/issues/bfem_d3c_broken_space_allowlist_issue.md)
    static_assert(D == 2 || D == 3, "bfem::broken_space: only D == 2 or D == 3");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef broken_field<D, T, P> field_type;

    broken_space(const mesh<D, T>& msh, int l)
        : l_(l), nv_(msh.num_vertices()), nt_(msh.num_elements()),
          nloc_(0), dm_(), geom_(), op_(), buf_(), loc_(), ub_(), vb_() {
        bfem_scalar_traits<T>::require();   // C-1 contract (L4, additive)
        if (l < 0)
            throw std::invalid_argument("bfem::broken_space: l must be >= 0");
        nloc_ = coeff_registry<D>::indices(l).size();
        dm_ = detail::broken_dofmap(nloc_, nt_);
        geom_.reserve(static_cast<std::size_t>(nt_));
        for (int e = 0; e < nt_; ++e) {
            std::array<std::array<T, D>, D + 1> vv;
            for (int c = 0; c <= D; ++c)
                vv[static_cast<std::size_t>(c)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(c)]);
            geom_.push_back(element_geometry<D, T>::from_vertices(vv));
        }
    }

    int order() const { return l_; }
    int ndof() const { return nt_ * nloc_; }
    int local_size() const { return nloc_; }
    int num_elements() const { return nt_; }
    int num_vertices() const { return nv_; }

    const detail::broken_dofmap& dofs() const { return dm_; }
    const element_geometry<D, T>& geometry(int e) const {
        assert(e >= 0 && e < nt_);
        return geom_[static_cast<std::size_t>(e)];
    }

    // ---- field factories ----
    field_type zero_field() {
        field_type f(l_);
        f.c_.zeros(ndof(), 1);
        return f;
    }
    field_type field_from_coeffs(vcp::matrix<T, P> c) {
        if (c.rowsize() != ndof() || c.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::broken_space::field_from_coeffs: size != ndof x 1");
        field_type f(l_);
        f.c_ = std::move(c);
        return f;
    }

    // ---- K2: block diagonal mass (|T| M^{(l,l)} per element) ----
    spmatrix_t mass() {
        buf_.clear();
        buf_.reserve(static_cast<std::size_t>(nt_)
                     * static_cast<std::size_t>(nloc_)
                     * static_cast<std::size_t>(nloc_));
#if VCP_BFEM_USE_OPENMP
        const int nt = nt_;
        int nrun = 1;
        std::vector<detail::coo_buffer<T> > tbuf;
#pragma omp parallel
        {
#pragma omp single
            {
                nrun = omp_get_num_threads();
                tbuf.resize(static_cast<std::size_t>(nrun));
            }
#pragma omp barrier
            const int tid = omp_get_thread_num();
            const int e0 = static_cast<int>((static_cast<long long>(nt) * tid) / nrun);
            const int e1 = static_cast<int>((static_cast<long long>(nt) * (tid + 1)) / nrun);
            element_op<D, T, P> op_l;
            vcp::matrix<T, P> loc_l;
            detail::coo_buffer<T>& b = tbuf[static_cast<std::size_t>(tid)];
            b.reserve(static_cast<std::size_t>(e1 - e0)
                      * static_cast<std::size_t>(nloc_)
                      * static_cast<std::size_t>(nloc_));
            for (int e = e0; e < e1; ++e) {
                op_l.set_geometry(geom_[static_cast<std::size_t>(e)]);
                op_l.local_mass(l_, l_, loc_l);
                detail::scatter_matrix<T>(dm_, e, loc_l, nloc_, b,
                                          detail::broken_dofmap::family_tag());
            }
        }
        buf_.append_all(tbuf, nrun);
#else
        for (int e = 0; e < nt_; ++e) {              // element order (X9)
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            op_.local_mass(l_, l_, loc_);
            // block-contiguous identity scatter (no cross-block entry is
            // ever pushed -- RE-5's block diagonality by construction)
            detail::scatter_matrix<T>(dm_, e, loc_, nloc_, buf_,
                                      detail::broken_dofmap::family_tag());
        }
#endif
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(ndof(), ndof(), buf_);
    }

    // ---- PF-1 L2-2: load (f, q_i) for a coordinate polynomial field ----
    //
    // L2-2 の用途(重要)(PF-1 設計書 §4 より全文転記):
    // broken_space::load は荷重ベクトル (f, q_i) の計算を提供するだけであり、
    // NS の非圧縮制約の右辺は 0 のままなので NS 本体では呼ばれない。SV の
    // 厳密 div-free 性は制約側((div u_h, q_h) = 0 ∀q_h、div V_h ⊆ Q_h)の
    // 性質であって右辺の API とは独立に保たれる。L2-2 の実際の用途は:
    //
    //   (i)   hypercircle の制約右辺 (ϵ div p_h − g, q_h) = 0 における
    //         (g, q_h) の計算
    //   (ii)  div u = g 型(圧縮性・湧き出しあり)の問題を意図的に解く場合
    //   (iii) 双対問題・誤差評価での汎関数計算
    //
    // 誤って NS の圧力ブロック右辺に置けば div u_h = g を課したことになる。
    // その帰結は API のバグではなく問題設定の変更だが、誤用防止のため
    // ここに明記する。
    //
    // The polynomial degree of this space is fixed at construction, so no
    // degree argument is taken (design section 2, rule 4); the integrand is
    // restricted at n = max(f.total_degree(), order()). Element loop /
    // local_load kernel / identity scatter mirror
    // fe_space::load(poly_field, m). broken_space stores no mesh member, so
    // the mesh view required by restrict_to is rebuilt once per call from
    // the stored element geometries -- vertex order is exactly the element
    // vertex order used at construction, hence the restriction is identical.
    vcp::matrix<T, P> load(const poly_field<D, T>& f) {
        vcp::matrix<T, P> F;
        F.zeros(ndof(), 1);
        if (f.is_zero()) return F;
        const int n = f.total_degree() > l_ ? f.total_degree() : l_;
        std::vector<std::array<T, D> > verts;
        verts.reserve(static_cast<std::size_t>(nt_)
                      * static_cast<std::size_t>(D + 1));
        std::vector<std::array<int, D + 1> > elems(
            static_cast<std::size_t>(nt_));
        for (int e = 0; e < nt_; ++e) {
            const std::array<std::array<T, D>, D + 1>& vv =
                geom_[static_cast<std::size_t>(e)].vertices();
            for (int k = 0; k <= D; ++k) {
                elems[static_cast<std::size_t>(e)][static_cast<std::size_t>(k)] =
                    static_cast<int>(verts.size());
                verts.push_back(vv[static_cast<std::size_t>(k)]);
            }
        }
        mesh<D, T> msh =
            mesh<D, T>::from_lists(std::move(verts), std::move(elems));
        bpoly<D, T> w;
        for (int e = 0; e < nt_; ++e) {              // element order (X9)
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            w = f.restrict_to(msh, e, n);
            op_.local_load(w, l_, loc_);
            detail::scatter_vector(dm_, e, loc_.data(), nloc_, F,
                                   detail::broken_dofmap::family_tag());
        }
        return F;
    }

    // ---- L2 scalar: (u, v)_{L2(Omega)} of two broken fields ----
    // const per external design section 4; the gather buffers are mutable
    // (the instance stays externally immutable; not thread safe, as usual)
    // ---- GRF-1: point evaluation on one element ----
    // Barycentric-defined block (W-RT9): the element block IS the bpoly
    // coefficient vector in canonical L0 order, so evaluation is de
    // Casteljau on the gathered block. const with mutable gather buffer,
    // exactly the inner() pattern (externally immutable, not thread safe).
    T eval(const field_type& u, int e, const bary_point<D, T>& lam) const {
        validate(u);
        if (e < 0 || e >= nt_)
            throw std::invalid_argument(
                "bfem::broken_space::eval: element out of range");
        gather_block(u, e, ub_);
        return ::vcp::bfem::eval(ub_, lam);
    }

    T inner(const field_type& u, const field_type& v) const {
        validate(u);
        validate(v);
        T acc(0);
        for (int e = 0; e < nt_; ++e) {
            gather_block(u, e, ub_);
            gather_block(v, e, vb_);
            acc += geom_[static_cast<std::size_t>(e)].measure()
                   * ::vcp::bfem::inner(ub_, vb_);
        }
        return acc;
    }

private:
    int l_, nv_, nt_, nloc_;
    detail::broken_dofmap dm_;
    std::vector<element_geometry<D, T> > geom_;
    element_op<D, T, P> op_;
    detail::coo_buffer<T> buf_;
    vcp::matrix<T, P> loc_;
    mutable bpoly<D, T> ub_, vb_;

    void validate(const field_type& u) const {
        if (u.order() != l_ || u.coeffs().rowsize() != ndof()
            || u.coeffs().columnsize() != 1)
            throw std::invalid_argument(
                "bfem::broken_space: field does not match this space");
    }
    void gather_block(const field_type& u, int e, bpoly<D, T>& dst) const {
        detail::bpoly_access::prepare(dst, l_, false);
        std::vector<T>& d = detail::bpoly_access::vec(dst);
        for (int r = 0; r < nloc_; ++r)
            d[static_cast<std::size_t>(r)] = u.coeffs()(e * nloc_ + r, 0);
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_BROKEN_SPACE_HPP
