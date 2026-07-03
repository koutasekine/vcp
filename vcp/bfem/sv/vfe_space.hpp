// vcp/bfem/sv/vfe_space.hpp
// Phase 5d (Scott-Vogelius parts): vfe_function and vfe_space<D,T,P,SP>
// (V1, SV-2) -- the component-major vector-valued wrapper of fe_space, plus
// the detail scatter plumbing (vfe_scatter) reused by sv_assemble.
//
// Conforms to: SV external design v0.2 (section 3) and
//              SV internal design v0.2 (section 5).
//
// Component-major numbering (SV-2, normative): with N = scalar ndof(m), the
// global dof of component d, scalar dof i is d * N + i.
//
// Lifetime / thread contract (v0.2, B-4): vfe_space wraps the scalar fe_space
// BY REFERENCE and does not own it; the referee must outlive the vfe_space
// (caller's responsibility). The thread model is one with the referee: no
// concurrent calls across vfe_space instances wrapping the SAME fe_space
// (instance duplication must duplicate the fe_space as well).
//
// Geometry note (recorded supplement to the v0.2 sketch): the frozen fe_space
// exposes no per-element geometry, so vfe_space -- exactly like the frozen
// broken_space -- receives the mesh at construction and computes its own
// element geometries once (one division per element inside
// element_geometry::from_vertices; no division appears in this file). The
// K-style free assemblers of sv_assemble.hpp borrow vfe_space::geometry(e)
// the same way the RT K functions borrow broken_space::geometry(e).
//
// This header is new phase-5d code; no frozen file is touched. The scatter
// plumbing below reuses the frozen Y2 kernels and coo_buffer; it adds no new
// scatter skeleton (S-SV-3):
//   - sv_offset_dofmap<DM> is a shifted VIEW of a scalar dofmap; the actual
//     matrix/vector scatter runs through the frozen detail::scatter_matrix /
//     detail::scatter_vector identity kernels (dofmap.hpp, pn_family_tag);
//   - the rectangular two-dofmap identity scatter (broken rows x fe columns
//     and the (c,d) block form) uses plain buf.push exactly like the frozen
//     K9 assemble_mixed_mass (rt_assemble.hpp) -- no sign multiplication.

#ifndef VCP_BFEM_SV_VFE_SPACE_HPP
#define VCP_BFEM_SV_VFE_SPACE_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/d3/dofmap3.hpp>
#include <vcp/bfem/fe_function.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/geometry.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class vfe_space;

// ---------------------------------------------------------------------------
// vfe_function<D,T,P>: degree + coefficient column of length D * N(m)
// (component major). Created only through the vfe_space factories; component
// extraction / write-back go through vfe_space (which owns the scalar-space
// linkage, per the B-4 lifetime contract).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class vfe_function {
public:
    int degree() const { return deg_; }
    const vcp::matrix<T, P>& coeffs() const { return c_; }
    vcp::matrix<T, P>&       coeffs() { return c_; }        // direct edit allowed

private:
    explicit vfe_function(int m) : deg_(m), c_() {}
    int deg_;
    vcp::matrix<T, P> c_;                                   // (D * ndof(m)) x 1
    template <int DD, typename TT, typename PP, class SS> friend class vfe_space;
};

namespace detail {

// shifted view of a scalar dofmap (vfe_scatter): drives the frozen Y2
// kernels; the family tag (and therefore the identity/no-sign dispatch) is
// inherited from the wrapped dofmap.
template <typename DM>
struct sv_offset_dofmap {
    typedef typename DM::family_tag family_tag;
    const DM* dm;
    int off;
    sv_offset_dofmap(const DM& d, int o) : dm(&d), off(o) {}
    int global_dof(int e, int r) const { return off + dm->global_dof(e, r); }
    int dof_sign(int e, int r) const { return dm->dof_sign(e, r); }
};

// rectangular identity x identity scatter with row/column offsets (the K9
// plain-push form; both families are P^n / broken identity, so no sign
// multiplication appears)
template <typename T, typename Buf, typename Loc, typename DMR, typename DMC>
void sv_scatter_rect_identity(const DMR& dmr, int row_off,
                              const DMC& dmc, int col_off,
                              int e, const Loc& loc, int nr, int nc, Buf& buf) {
    for (int i = 0; i < nr; ++i) {
        int gi = row_off + dmr.global_dof(e, i);
        for (int j = 0; j < nc; ++j)
            buf.push(gi, col_off + dmc.global_dof(e, j), loc(i, j));
    }
}

} // namespace detail

// ---------------------------------------------------------------------------
// vfe_space<D, T, P, SP> (V1)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class vfe_space {
    static_assert(D == 2 || D == 3, "bfem::vfe_space: only D == 2 or D == 3");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef vfe_function<D, T, P> function_type;
    typedef fe_function<D, T, P> scalar_function_type;

    vfe_space(const mesh<D, T>& msh, fe_space<D, T, P, SP>& scalar)
        : scalar_(&scalar), nv_(msh.num_vertices()), nt_(msh.num_elements()),
          geom_() {
        if (scalar.num_elements() != msh.num_elements())
            throw std::invalid_argument(
                "bfem::vfe_space: fe_space is not built from this mesh");
        geom_.reserve(static_cast<std::size_t>(nt_));
        for (int e = 0; e < nt_; ++e) {
            std::array<std::array<T, D>, D + 1> vv;
            for (int k = 0; k <= D; ++k)
                vv[static_cast<std::size_t>(k)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(k)]);
            geom_.push_back(element_geometry<D, T>::from_vertices(vv));
        }
    }

    // the wrapped scalar space (non-owning reference, B-4)
    fe_space<D, T, P, SP>& scalar() const { return *scalar_; }
    int num_elements() const { return nt_; }
    int num_vertices() const { return nv_; }
    const element_geometry<D, T>& geometry(int e) const {
        assert(e >= 0 && e < nt_);
        return geom_[static_cast<std::size_t>(e)];
    }

    int ndof(int m) const { return D * scalar_->ndof(m); }

    // component-major global numbering (SV-2): d * N + i
    int global_dof(int d, int i, int m) const {
        const int N = scalar_->ndof(m);
        if (d < 0 || d >= D)
            throw std::invalid_argument("bfem::vfe_space::global_dof: bad component");
        if (i < 0 || i >= N)
            throw std::invalid_argument("bfem::vfe_space::global_dof: bad scalar dof");
        return d * N + i;
    }

    // ---- factories ----
    function_type zero_function(int m) const {
        function_type f(m);
        f.c_.zeros(ndof(m), 1);
        return f;
    }
    function_type function_from_coeffs(int m, vcp::matrix<T, P> c) const {
        if (c.rowsize() != ndof(m) || c.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::vfe_space::function_from_coeffs: size != D * ndof(m) x 1");
        function_type f(m);
        f.c_ = std::move(c);
        return f;
    }

    // ---- boundary dofs: replicate the scalar list into every component ----
    std::vector<int> boundary_dofs(int m) const {
        std::vector<int> s = scalar_->dofs(m).boundary_dofs();
        const int N = scalar_->ndof(m);
        std::vector<int> out;
        out.reserve(static_cast<std::size_t>(D) * s.size());
        for (int d = 0; d < D; ++d)
            for (std::size_t k = 0; k < s.size(); ++k)
                out.push_back(d * N + s[k]);
        return out;
    }

    // ---- component view (copy extraction) and write-back ----
    scalar_function_type component(const function_type& u, int d) const {
        validate(u);
        if (d < 0 || d >= D)
            throw std::invalid_argument("bfem::vfe_space::component: bad component");
        const int N = scalar_->ndof(u.degree());
        vcp::matrix<T, P> c;
        c.zeros(N, 1);
        for (int i = 0; i < N; ++i) c(i, 0) = u.coeffs()(d * N + i, 0);
        return scalar_->function_from_coeffs(u.degree(), std::move(c));
    }
    void set_component(function_type& u, int d,
                       const scalar_function_type& f) const {
        validate(u);
        if (d < 0 || d >= D)
            throw std::invalid_argument("bfem::vfe_space::set_component: bad component");
        if (f.degree() != u.degree())
            throw std::invalid_argument("bfem::vfe_space::set_component: degree mismatch");
        const int N = scalar_->ndof(u.degree());
        if (f.coeffs().rowsize() != N || f.coeffs().columnsize() != 1)
            throw std::invalid_argument("bfem::vfe_space::set_component: size mismatch");
        for (int i = 0; i < N; ++i) u.coeffs()(d * N + i, 0) = f.coeffs()(i, 0);
    }

    void validate(const function_type& u) const {
        if (u.coeffs().columnsize() != 1
            || u.coeffs().rowsize() != ndof(u.degree()))
            throw std::invalid_argument(
                "bfem::vfe_space: vfe_function does not match this space");
    }

private:
    fe_space<D, T, P, SP>* scalar_;                          // non-owning (B-4)
    int nv_, nt_;
    std::vector<element_geometry<D, T> > geom_;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_VFE_SPACE_HPP
