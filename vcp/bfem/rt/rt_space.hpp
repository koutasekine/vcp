// vcp/bfem/rt/rt_space.hpp
// RT Layer 3: rt_dofmap (detail), rt_field and rt_space<D,T,P,SP> (K1, K6).
//
// Conforms to: RT-L3 external design v0.2 (sections 2, 3, 6) and
//              RT-L3 internal design v0.2 (sections 2, 4, 5).
//
// Sign rule (Q1, normative -- external design section 2):
//   global edge DOF = moment against nu_global = R_{-90}(canonical physical
//   edge vector), canonical direction = global vertex small -> large.
//   local -> global: fwd ? (sign +1, index j) : (sign -1, index k - j).
//   ORIENT-FREE: the det sign never appears here (it is consumed by the
//   RT-L2 f_div factor -- the inter-layer sign division of labor).
//   Interior DOFs: element private, sign +1.
// Numbering: edge block (edge id order, j = 0..k in canonical direction)
// then interior block (element order, RT-L0 canonical interior order).
// Index reversal is folded into l2g; signs are a separate +-1 array
// (internal design 2.1).
//
// rt_interp (K6, Q5 test-support scope): pullback via adj(B) (sign
// permutation of B entries, multiplications only -- zero added divisions),
// reference DOF application with the typed L0 1D/2D mass tables (the SAME
// convention as the RT-L0 rational-stage DOF construction, including the
// |T_hat| = 1/D! factor on interior moments), then the sign rule above,
// write-once in element order. Precondition: the input field is normal-
// continuous (H(div)-conforming); the gradient of a conforming P^n function
// is NOT (tangential continuity only) -- the classic trap (A-1).
// Interval T: both sides of a shared edge enclose the same exact value, so
// last-write-wins stays sound (C-1).

#ifndef VCP_BFEM_RT_RT_SPACE_HPP
#define VCP_BFEM_RT_RT_SPACE_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/fe_space.hpp>          // detail::coo_buffer / spm_adapter
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/rt/rt_tables.hpp>
#include <vcp/bfem/rt/rt_typed_tables.hpp>
#include <vcp/bfem/rt/rt_element_op.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class rt_space;

namespace detail {

// ---------------------------------------------------------------------------
// rt_dofmap: RT^k local-to-global map (index reversal folded in) plus the
// +-1 sign array. Drives the GENERAL (signed) Y2 kernels of dofmap.hpp via
// general_family_tag.
// ---------------------------------------------------------------------------
class rt_dofmap {
public:
    typedef general_family_tag family_tag;

    rt_dofmap() : k_(0), ne_(0), nt_(0), nloc_(0), nint_(0), ndof_(0) {}

    static rt_dofmap build(const mesh_topology2& tp, int k) {
        if (k < 0)
            throw std::invalid_argument("bfem::rt_dofmap: k < 0");
        rt_dofmap dm;
        dm.k_ = k;
        dm.ne_ = tp.num_edges();
        dm.nt_ = tp.nt;
        dm.nint_ = k * (k + 1);                      // 2 N(2, k-1)
        dm.nloc_ = 3 * (k + 1) + dm.nint_;           // == dim(k)
        dm.ndof_ = dm.ne_ * (k + 1) + dm.nt_ * dm.nint_;
        dm.l2g_.assign(static_cast<std::size_t>(dm.nt_)
                       * static_cast<std::size_t>(dm.nloc_), -1);
        dm.sgn_.assign(dm.l2g_.size(), static_cast<signed char>(1));
        const int int_base = dm.ne_ * (k + 1);
        for (int e = 0; e < dm.nt_; ++e) {
            for (int s = 0; s < 3; ++s) {
                int ed = tp.tri_edge[static_cast<std::size_t>(e)]
                                    [static_cast<std::size_t>(s)];
                bool fwd = tp.tri_edge_fwd[static_cast<std::size_t>(e)]
                                          [static_cast<std::size_t>(s)];
                for (int j = 0; j <= k; ++j) {
                    int r = s * (k + 1) + j;
                    std::size_t idx = static_cast<std::size_t>(e)
                                      * static_cast<std::size_t>(dm.nloc_)
                                      + static_cast<std::size_t>(r);
                    dm.l2g_[idx] = ed * (k + 1) + (fwd ? j : k - j);
                    dm.sgn_[idx] = static_cast<signed char>(fwd ? 1 : -1);
                }
            }
            for (int r = 0; r < dm.nint_; ++r) {
                std::size_t idx = static_cast<std::size_t>(e)
                                  * static_cast<std::size_t>(dm.nloc_)
                                  + static_cast<std::size_t>(3 * (k + 1) + r);
                dm.l2g_[idx] = int_base + e * dm.nint_ + r;
            }
        }
        return dm;
    }

    int order() const { return k_; }
    int ndof() const { return ndof_; }
    int local_size() const { return nloc_; }
    int num_elements() const { return nt_; }
    int num_edge_dofs() const { return ne_ * (k_ + 1); }

    int global_dof(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(r)];
    }
    int dof_sign(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return sgn_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(r)];
    }

private:
    int k_, ne_, nt_, nloc_, nint_, ndof_;
    std::vector<int> l2g_;
    std::vector<signed char> sgn_;
};

} // namespace detail

// ---------------------------------------------------------------------------
// rt_field<D,T,P> (B-2): degree + global coefficients in the GLOBAL
// convention of external design 2.2 (signs/index reversal absorbed by the
// dofmap). Factory-created; coefficients are directly editable.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class rt_field {
public:
    int order() const { return k_; }
    const vcp::matrix<T, P>& coeffs() const { return c_; }
    vcp::matrix<T, P>&       coeffs() { return c_; }

private:
    explicit rt_field(int k) : k_(k), c_() {}
    int k_;
    vcp::matrix<T, P> c_;
    template <int DD, typename TT, typename PP, class SS> friend class rt_space;
};

// ---------------------------------------------------------------------------
// rt_space<D,T,P,SP> (K1). Contract inherited from fe_space section 5.1:
// all-dof output, non-const generation members, one instance NOT thread
// safe, element-order accumulation, deterministic duplicate combination,
// bit reproducibility, SP lazy realization rules. Space members reuse
// instance buffers (allocation-free from the second call on, A-2).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class rt_space {
    static_assert(D == 2, "bfem::rt_space: initial version supports D == 2 only");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef rt_field<D, T, P> field_type;

    // geometries are computed at construction (degeneracy fails fast, as in
    // fe_space)
    rt_space(const mesh<D, T>& msh, int k)
        : k_(k), topo_(), dm_(), geom_(), op_(), buf_(), sloc_(), loc_() {
        if (k < 0)
            throw std::invalid_argument("bfem::rt_space: k must be >= 0");
        topo_ = detail::mesh_topology2::build(msh);
        dm_ = detail::rt_dofmap::build(topo_, k);
        geom_.reserve(static_cast<std::size_t>(topo_.nt));
        for (int e = 0; e < topo_.nt; ++e) {
            std::array<std::array<T, D>, D + 1> vv;
            for (int c = 0; c <= D; ++c)
                vv[static_cast<std::size_t>(c)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(c)]);
            geom_.push_back(element_geometry<D, T>::from_vertices(vv));
        }
        sloc_.k = k;
        sloc_.c.assign(static_cast<std::size_t>(dm_.local_size()), T(0));
    }

    int order() const { return k_; }
    int ndof() const { return dm_.ndof(); }
    int num_elements() const { return topo_.nt; }
    int num_vertices() const { return topo_.nv; }

    // accessors for the cross-family free functions and the tests (internal
    // design 4: free functions "borrow dofmap and geometry from the spaces")
    const detail::rt_dofmap& dofs() const { return dm_; }
    const element_geometry<D, T>& geometry(int e) const {
        assert(e >= 0 && e < topo_.nt);
        return geom_[static_cast<std::size_t>(e)];
    }

    // ---- field factories (B-2) ----
    field_type zero_field() {
        field_type f(k_);
        f.c_.zeros(ndof(), 1);
        return f;
    }
    field_type field_from_coeffs(vcp::matrix<T, P> c) {
        if (c.rowsize() != ndof() || c.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::rt_space::field_from_coeffs: size != ndof x 1");
        field_type f(k_);
        f.c_ = std::move(c);
        return f;
    }

    // ---- K1: global RT mass (H1) ----
    spmatrix_t mass() {
        const int nloc = dm_.local_size();
        buf_.clear();
        buf_.reserve(static_cast<std::size_t>(topo_.nt)
                     * static_cast<std::size_t>(nloc) * static_cast<std::size_t>(nloc));
        for (int e = 0; e < topo_.nt; ++e) {             // element order (X9)
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            op_.local_rt_mass(k_, loc_);
            // Y2 GENERAL kernel: signed scatter (the promotion's production
            // entry point -- S-RT3-1(ii))
            detail::scatter_matrix<T>(dm_, e, loc_, nloc, buf_,
                                      detail::rt_dofmap::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm_.ndof(), dm_.ndof(), buf_);
    }

    // ---- K6: RT interpolation (Q5, test-support scope) ----
    // Provider: e -> std::pair<bpoly σ_x, bpoly σ_y> (PHYSICAL components,
    // barycentric bpoly). Precondition: normal-continuous input (see header
    // note); shared-edge DOFs are written last-write-wins in element order.
    template <typename Provider>
    field_type interpolate(Provider f) {
        field_type out = zero_field();
        const int k = k_;
        const int nint2 = (k >= 1) ? coeff_registry<D>::indices(k - 1).size() : 0;
        for (int e = 0; e < topo_.nt; ++e) {
            std::pair<bpoly<D, T>, bpoly<D, T> > s = f(e);
            const element_geometry<D, T>& g = geom_[static_cast<std::size_t>(e)];
            // pullback: sigma_hat = adj(B) (sigma o F); adj(B) is a signed
            // permutation of the B entries (multiplications only)
            const T a00 = g.edge_matrix(1, 1), a01 = -g.edge_matrix(0, 1);
            const T a10 = -g.edge_matrix(1, 0), a11 = g.edge_matrix(0, 0);
            scale_into(ibuf_[0], s.first, a00);
            scale_into(ibuf_[1], s.second, a01);
            add_into(hat_[0], ibuf_[0], ibuf_[1]);
            scale_into(ibuf_[0], s.first, a10);
            scale_into(ibuf_[1], s.second, a11);
            add_into(hat_[1], ibuf_[0], ibuf_[1]);
            // reference DOF application (RT-L0 section 4 in T arithmetic)
            for (int sdg = 0; sdg < 3; ++sdg) {
                bool fwd = topo_.tri_edge_fwd[static_cast<std::size_t>(e)]
                                             [static_cast<std::size_t>(sdg)];
                int ed = topo_.tri_edge[static_cast<std::size_t>(e)]
                                       [static_cast<std::size_t>(sdg)];
                for (int j = 0; j <= k; ++j) {
                    T ell(0);
                    for (int d = 0; d < 2; ++d) {
                        int nu = detail::rt_edge_normal(sdg, d);
                        if (nu == 0) continue;
                        const bpoly<D, T>& hd = hat_[static_cast<std::size_t>(d)];
                        const int deg = hd.degree();
                        const typed_mass_table<1, T>& M1 =
                            typed_registry<1, T>::mass(k, deg);
                        T acc(0);
                        for (int t = 0; t <= deg; ++t)
                            acc += M1.at(j, t)
                                   * hd.coeff(detail::rt_trace_index(sdg, deg, t));
                        ell += (nu > 0) ? acc : -acc;    // nonzero nu is +-1
                    }
                    // global convention (external 2.2 inverse): fwd keeps
                    // (index, sign); reversed flips both
                    int gidx = ed * (k + 1) + (fwd ? j : k - j);
                    out.c_(gidx, 0) = fwd ? ell : -ell;  // write-once assignment
                }
            }
            if (k >= 1) {
                const int int_base = dm_.num_edge_dofs() + e * (2 * nint2);
                for (int d = 0; d < 2; ++d) {
                    const bpoly<D, T>& hd = hat_[static_cast<std::size_t>(d)];
                    const typed_mass_table<D, T>& M2 =
                        typed_registry<D, T>::mass(k - 1, hd.degree());
                    for (int ar = 0; ar < nint2; ++ar) {
                        T acc = M2.at(ar, 0) * hd.coeff(0);
                        for (int b = 1; b < M2.cols(); ++b)
                            acc += M2.at(ar, b) * hd.coeff(b);
                        // interior moment includes |T_hat| = 1/D! (the RT-L0
                        // normative DOF convention)
                        out.c_(int_base + d * nint2 + ar, 0) = that_const() * acc;
                    }
                }
            }
        }
        return out;
    }

private:
    int k_;
    detail::mesh_topology2 topo_;
    detail::rt_dofmap dm_;
    std::vector<element_geometry<D, T> > geom_;
    rt_element_op<D, T, P> op_;
    detail::coo_buffer<T> buf_;
    rt_local_coeffs<T> sloc_;
    vcp::matrix<T, P> loc_;
    bpoly<D, T> ibuf_[2], hat_[2];

    static const T& that_const() {                       // |T_hat| = 1/D!
        static const T c = rational_to<T>(1, detail::factorial_of<D>::value);
        return c;
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_SPACE_HPP
