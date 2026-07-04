// vcp/bfem/c1/c1_space.hpp
// Phase 6 (2D C1 Argyris family): the global layer c1_space<D=2,T,P,SP> --
// H-A (full P^k vocabulary + eval_grad / eval_hess), H-B (the Delta residual
// system), c1 interpolation and c1_elevate (DOF re-interpolation).
//
// Conforms to: C1 external design v0.2 (sections 4, 5, 6) and
//              C1 internal design v0.2 (section 6).
//
// Construction order (B-4): topology -> edge cache (fail-fast BEFORE the
// per-edge division on uncertified |t_e|^2) -> element geometries (fail-fast
// on degenerate elements; one division per element). The per-program T
// division count of a c1_space therefore settles at n_t + n_e (S-C1-2).
//
// Assembly: the fe_space skeleton is reused verbatim -- element-order
// accumulation, detail::coo_buffer deterministic combination,
// detail::spm_adapter build, SP lazy realization; scatter/gather run the
// FROZEN Y2 general (signed) kernels driven by c1_dofmap
// (general_family_tag). No matrix-weighted scatter exists (C1-4 / S-C1-5);
// all geometry mixing is closed inside the element layer (c1_element_op).
//
// H-B: laplacian_field supplies Delta u_h as a broken P_{k-2} field -- the
// splice point for the user's correction-term quadrature (the bfem side
// integrates polynomials only; see the responsibility split of external
// design section 1). laplacian_residual_sq is the K8/RG6-style convenience
// scalar || Delta u_h + f(u_h) ||^2.
//
// c1_elevate (C1-10): vertex DOFs are degree-independent functionals and
// are COPIED (never re-evaluated -- verified by counting in C-T7); edge and
// interior DOFs are exactly re-interpolated per element through the
// partial closed-form inverse of M_T.

#ifndef VCP_BFEM_C1_C1_SPACE_HPP
#define VCP_BFEM_C1_C1_SPACE_HPP

#include <vector>
#include <array>
#include <map>
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
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_assemble.hpp>    // detail::scatter_matrix_general2
#include <vcp/bfem/c1/c1_tables.hpp>
#include <vcp/bfem/c1/c1_geometry.hpp>
#include <vcp/bfem/c1/c1_dofmap.hpp>
#include <vcp/bfem/c1/c1_element_op.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class c1_space;

// ---------------------------------------------------------------------------
// c1_function<D,T,P>: degree + global C1 coefficients (factory-created)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class c1_function {
public:
    int degree() const { return deg_; }
    const vcp::matrix<T, P>& coeffs() const { return c_; }
    vcp::matrix<T, P>&       coeffs() { return c_; }

private:
    explicit c1_function(int m) : deg_(m), c_() {}
    int deg_;
    vcp::matrix<T, P> c_;
    template <int DD, typename TT, typename PP, class SS> friend class c1_space;
};

// ---------------------------------------------------------------------------
// c1_space<D = 2, T, P, SP>
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class c1_space {
    static_assert(D == 2, "bfem::c1_space: D == 2 only (C1 scope)");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef c1_function<D, T, P> function_type;

    c1_space(const mesh<2, T>& msh, int k)
        : mesh_(msh), k_(k), topo_(), edge_cache_(), dmaps_(), geom_(),
          op_(), buf_(), loc_(), uloc_(), vloc_(), wloc_(), rloc_(), cws_(),
          gbuf_(), hatbuf_() {
        detail::c1_check_k(k, "c1_space");
        topo_ = detail::mesh_topology2::build(msh);
        // edge cache FIRST (B-4: |t_e|^2 certification precedes everything;
        // one division per edge)
        edge_cache_ = detail::c1_edge_inv_tsq(msh, topo_);
        geom_.reserve(static_cast<std::size_t>(topo_.nt));
        for (int e = 0; e < topo_.nt; ++e) {
            std::array<std::array<T, 2>, 3> vv;
            for (int c = 0; c < 3; ++c)
                vv[static_cast<std::size_t>(c)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(c)]);
            geom_.push_back(element_geometry<2, T>::from_vertices(vv));
        }
    }

    // ---- observers ----
    int base_degree() const { return k_; }
    int num_elements() const { return topo_.nt; }
    int num_vertices() const { return topo_.nv; }
    int num_edges() const { return topo_.num_edges(); }
    const mesh<2, T>& mesh_ref() const { return mesh_; }
    const element_geometry<2, T>& geometry(int e) const {
        assert(e >= 0 && e < topo_.nt);
        return geom_[static_cast<std::size_t>(e)];
    }
    const T& edge_inv_tsq(int ed) const {
        assert(ed >= 0 && ed < topo_.num_edges());
        return edge_cache_[static_cast<std::size_t>(ed)];
    }
    // family degrees m are >= the base degree (AX-C1-3: m < k is rejected)
    const c1_dofmap& dofs(int m) {
        if (m < k_)
            throw std::invalid_argument("bfem::c1_space::dofs: m < base degree");
        typename std::map<int, c1_dofmap>::iterator it = dmaps_.find(m);
        if (it != dmaps_.end()) return it->second;
        return dmaps_.insert(std::make_pair(m, c1_dofmap::build(topo_, m)))
            .first->second;
    }
    int ndof(int m) { return dofs(m).ndof(); }

    // ---- function factories ----
    function_type zero_function(int m) {
        function_type f(m);
        f.c_.zeros(ndof(m), 1);
        return f;
    }
    function_type function_from_coeffs(int m, vcp::matrix<T, P> c) {
        if (c.rowsize() != ndof(m) || c.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::c1_space::function_from_coeffs: size != ndof(m) x 1");
        function_type f(m);
        f.c_ = std::move(c);
        return f;
    }

    // ---- H-A: global matrices ----
    spmatrix_t stiffness(int m) {
        const c1_dofmap& dm = dofs(m);
        begin_matrix(dm, dm);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            op_.local_stiffness(m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename c1_dofmap::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }
    spmatrix_t mixed_mass(int a, int b) {
        const c1_dofmap& dma = dofs(a);
        const c1_dofmap& dmb = dofs(b);
        begin_matrix(dma, dmb);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            op_.local_mass(a, b, loc_);
            detail::scatter_matrix_general2<T>(dma, dmb, e, loc_,
                                               dma.local_size(),
                                               dmb.local_size(), buf_);
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dma.ndof(), dmb.ndof(), buf_);
    }
    spmatrix_t weighted_mass(const poly1<T>& fprime, const function_type& uh,
                             int m) {
        validate_fn(uh);
        const c1_dofmap& dmu = dofs(uh.degree());
        const c1_dofmap& dm = dofs(m);
        begin_matrix(dm, dm);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, fprime, uloc_, cws_);
            op_.local_weighted_mass(wloc_, m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename c1_dofmap::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }
    vcp::matrix<T, P> load(const poly1<T>& f, const function_type& uh, int m) {
        validate_fn(uh);
        const c1_dofmap& dmu = dofs(uh.degree());
        const c1_dofmap& dm = dofs(m);
        vcp::matrix<T, P> F;
        F.zeros(dm.ndof(), 1);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, f, uloc_, cws_);
            op_.local_load(wloc_, m, loc_);
            detail::scatter_vector(dm, e, loc_.data(), dm.local_size(), F,
                                   typename c1_dofmap::family_tag());
        }
        return F;
    }

    // ---- H-A: scalars and point evaluation ----
    T inner(const function_type& u, const function_type& v) {
        validate_fn(u);
        validate_fn(v);
        const c1_dofmap& dmu = dofs(u.degree());
        const c1_dofmap& dmv = dofs(v.degree());
        T acc(0);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, u, uloc_);
            gather_hat_fn(dmv, e, v, vloc_);
            acc += op_.local_inner(uloc_, vloc_);
        }
        return acc;
    }
    T scalar_ff(const poly1<T>& f, const function_type& uh) {
        validate_fn(uh);
        const c1_dofmap& dmu = dofs(uh.degree());
        T acc(0);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, f, uloc_, cws_);
            acc += op_.local_inner(wloc_, wloc_);
        }
        return acc;
    }
    T eval(const function_type& u, int e, const bary_point<2, T>& lam) {
        validate_fn(u);
        check_elem(e, "eval");
        set_elem(e);
        gather_hat_fn(dofs(u.degree()), e, u, uloc_);
        return ::vcp::bfem::eval(uloc_, lam);
    }
    std::array<T, 2> eval_grad(const function_type& u, int e,
                               const bary_point<2, T>& lam) {
        validate_fn(u);
        check_elem(e, "eval_grad");
        set_elem(e);
        gather_hat_fn(dofs(u.degree()), e, u, uloc_);
        std::array<T, 2> out;
        for (int d = 0; d < 2; ++d) {
            op_.grad_component(uloc_, d, rloc_);
            out[static_cast<std::size_t>(d)] = ::vcp::bfem::eval(rloc_, lam);
        }
        return out;
    }
    std::array<T, 3> eval_hess(const function_type& u, int e,
                               const bary_point<2, T>& lam) {
        validate_fn(u);
        check_elem(e, "eval_hess");
        set_elem(e);
        gather_hat_fn(dofs(u.degree()), e, u, uloc_);
        std::array<T, 3> out;
        op_.grad_component(uloc_, 0, vloc_);          // du/dx
        op_.grad_component(vloc_, 0, rloc_);
        out[0] = ::vcp::bfem::eval(rloc_, lam);       // dxx
        op_.grad_component(vloc_, 1, rloc_);
        out[1] = ::vcp::bfem::eval(rloc_, lam);       // dxy
        op_.grad_component(uloc_, 1, vloc_);          // du/dy
        op_.grad_component(vloc_, 1, rloc_);
        out[2] = ::vcp::bfem::eval(rloc_, lam);       // dyy
        return out;
    }

    // ---- H-B: the Delta residual system ----
    spmatrix_t laplacian_matrix(int m) {
        const c1_dofmap& dm = dofs(m);
        begin_matrix(dm, dm);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            op_.local_laplacian(m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename c1_dofmap::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }
    spmatrix_t hessian_matrix(int m) {                 // (D^2 u : D^2 v)
        const c1_dofmap& dm = dofs(m);
        begin_matrix(dm, dm);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            op_.local_hessian(m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename c1_dofmap::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }
    // rows: broken P_l (block numbering e N_l + r), cols: V^{C1,m}
    spmatrix_t laplacian_mixed(int m, int l) {
        if (l < 0)
            throw std::invalid_argument("bfem::c1_space::laplacian_mixed: l < 0");
        const c1_dofmap& dm = dofs(m);
        const int nl = coeff_registry<2>::indices(l).size();
        detail::broken_dofmap bdm(nl, topo_.nt);
        buf_.clear();
        buf_.reserve(static_cast<std::size_t>(topo_.nt)
                     * static_cast<std::size_t>(nl)
                     * static_cast<std::size_t>(dm.local_size()));
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            op_.local_laplacian_mixed(m, l, loc_);
            detail::scatter_matrix_general2<T>(bdm, dm, e, loc_, nl,
                                               dm.local_size(), buf_);
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(bdm.ndof(), dm.ndof(), buf_);
    }
    vcp::matrix<T, P> laplacian_load(const poly1<T>& f, const function_type& uh,
                                     int m) {
        validate_fn(uh);
        const c1_dofmap& dmu = dofs(uh.degree());
        const c1_dofmap& dm = dofs(m);
        vcp::matrix<T, P> F;
        F.zeros(dm.ndof(), 1);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, f, uloc_, cws_);
            op_.local_laplacian_load(wloc_, m, loc_);
            detail::scatter_vector(dm, e, loc_.data(), dm.local_size(), F,
                                   typename c1_dofmap::family_tag());
        }
        return F;
    }
    // || Delta u_h + f(u_h) ||^2_{L2(Omega)} (element-order accumulation)
    T laplacian_residual_sq(const poly1<T>& f, const function_type& uh) {
        validate_fn(uh);
        const c1_dofmap& dmu = dofs(uh.degree());
        T acc(0);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, f, uloc_, cws_);
            op_.laplacian_coeffs(uloc_, vloc_);
            add_into(rloc_, vloc_, wloc_);
            acc += op_.local_inner(rloc_, rloc_);
        }
        return acc;
    }
    // Delta u_h as a broken P_{k-2} field: the correction-term splice point
    template <class SPB>
    broken_field<2, T, P> laplacian_field(const function_type& u,
                                          broken_space<2, T, P, SPB>& bs) {
        validate_fn(u);
        if (bs.order() != u.degree() - 2)
            throw std::invalid_argument(
                "bfem::c1_space::laplacian_field: broken order != degree - 2");
        if (bs.num_elements() != topo_.nt || bs.num_vertices() != topo_.nv)
            throw std::invalid_argument(
                "bfem::c1_space::laplacian_field: different mesh");
        const c1_dofmap& dmu = dofs(u.degree());
        broken_field<2, T, P> out = bs.zero_field();
        const int nl = bs.local_size();
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmu, e, u, uloc_);
            op_.laplacian_coeffs(uloc_, vloc_);
            for (int r = 0; r < nl; ++r)
                out.coeffs()(e * nl + r, 0) = vloc_.coeff(r);
        }
        return out;
    }

    // ---- c1 interpolation: provider(e) -> local pullback polynomial ----
    // (degree <= m, elevated internally; input must be C1 conforming --
    // shared DOFs are written last-write-wins, exact for exact input)
    template <typename Provider>
    function_type interpolate(Provider f, int m) {
        const c1_dofmap& dm = dofs(m);
        function_type out = zero_function(m);
        const typename c1_typed_registry<T>::mat_table& L =
            c1_typed_registry<T>::dof_matrix(m);
        const int dim = dm.local_size();
        grow(gbuf_, dim);
        grow(hatbuf_, dim);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            elevate_into(uloc_, f(e), m);
            // full DOF application hat_c = L b, then c = P^{-1} hat_c
            const std::vector<T>& b = detail::bpoly_access::vec(uloc_);
            for (int r = 0; r < dim; ++r) {
                T acc = L.at(r, 0) * b[0];
                for (int j = 1; j < L.cols(); ++j)
                    acc += L.at(r, j) * b[static_cast<std::size_t>(j)];
                hatbuf_[static_cast<std::size_t>(r)] = acc;
            }
            op_.pull(m).apply_inverse(hatbuf_.data(), gbuf_.data());
            scatter_assign(dm, e, out);
        }
        return out;
    }

    // ---- c1_elevate (C1-10): vertex DOFs copied, the rest re-interpolated --
    function_type elevate(const function_type& u, int m) {
        validate_fn(u);
        const int k = u.degree();
        if (m < k)
            throw std::invalid_argument("bfem::c1_space::elevate: m < degree(u)");
        const c1_dofmap& dmk = dofs(k);
        const c1_dofmap& dmm = dofs(m);
        function_type out = zero_function(m);
        // vertex blocks are the SAME functionals at every degree: copy
        // (identical global indices 6v + c; no evaluation happens here)
        for (int i = 0; i < 6 * topo_.nv; ++i)
            out.c_(i, 0) = u.coeffs()(i, 0);
        if (m == k) {                                  // copy the rest as well
            for (int i = 6 * topo_.nv; i < dmk.ndof(); ++i)
                out.c_(i, 0) = u.coeffs()(i, 0);
            return out;
        }
        const typename c1_typed_registry<T>::mat_table& L =
            c1_typed_registry<T>::dof_matrix(m);
        const int dim = dmm.local_size();
        const int ntr = detail::c1_ntrace(m);
        const int nnd = detail::c1_nnd(m);
        const int nint = detail::c1_nint(m);
        grow(gbuf_, dim);
        grow(hatbuf_, dim);
        for (int e = 0; e < topo_.nt; ++e) {
            set_elem(e);
            gather_hat_fn(dmk, e, u, vloc_);
            elevate_into(uloc_, vloc_, m);
            const std::vector<T>& b = detail::bpoly_access::vec(uloc_);
            // local physical vertex entries: copies of the (already copied)
            // global values -- gathered, not evaluated
            for (int p = 0; p < 3; ++p) {
                int v = topo_.tri[static_cast<std::size_t>(e)][static_cast<std::size_t>(p)];
                for (int c = 0; c < 6; ++c)
                    gbuf_[static_cast<std::size_t>(detail::c1_vertex_dof(p, c))] =
                        out.c_(6 * v + c, 0);
            }
            // edge trace + interior: value rows of L (identity under P^{-1})
            for (int s = 0; s < 3; ++s)
                for (int i = 0; i < ntr; ++i)
                    row_into(L, detail::c1_edge_trace_dof(m, s, i), b);
            for (int r = 0; r < nint; ++r)
                row_into(L, detail::c1_interior_dof(m, r), b);
            // ND rows: hat values then the partial closed-form inverse
            const detail::c1_pullback<T>& Pm = op_.pull(m);
            for (int s = 0; s < 3; ++s) {
                for (int j = 0; j < nnd; ++j) {
                    int r = detail::c1_edge_nd_dof(m, s, j);
                    T acc = L.at(r, 0) * b[0];
                    for (int q = 1; q < L.cols(); ++q)
                        acc += L.at(r, q) * b[static_cast<std::size_t>(q)];
                    hatbuf_[static_cast<std::size_t>(j)] = acc;
                }
                Pm.invert_nd_edge(s, hatbuf_.data(), gbuf_.data());
            }
            // write edge + interior dofs (vertex dofs are already copied)
            scatter_assign_nonvertex(dmm, e, out);
        }
        return out;
    }

private:
    mesh<2, T> mesh_;
    int k_;
    detail::mesh_topology2 topo_;
    std::vector<T> edge_cache_;
    std::map<int, c1_dofmap> dmaps_;
    std::vector<element_geometry<2, T> > geom_;
    c1_element_op<T, P> op_;
    detail::coo_buffer<T> buf_;
    vcp::matrix<T, P> loc_;
    bpoly<2, T> uloc_, vloc_, wloc_, rloc_;
    compose_workspace<2, T> cws_;
    std::vector<T> gbuf_, hatbuf_;

    void check_elem(int e, const char* where) const {
        if (e < 0 || e >= topo_.nt) {
            std::string msg("bfem::c1_space::");
            msg += where;
            msg += ": element out of range";
            throw std::invalid_argument(msg);
        }
    }
    void validate_fn(const function_type& u) {
        if (u.coeffs().columnsize() != 1
            || u.coeffs().rowsize() != ndof(u.degree()))
            throw std::invalid_argument(
                "bfem::c1_space: c1_function does not match this space");
    }
    void set_elem(int e) {
        std::array<T, 3> lts;
        for (int s = 0; s < 3; ++s)
            lts[static_cast<std::size_t>(s)] = edge_cache_[static_cast<std::size_t>(
                topo_.tri_edge[static_cast<std::size_t>(e)][static_cast<std::size_t>(s)])];
        op_.set_geometry(geom_[static_cast<std::size_t>(e)], lts);
    }
    void gather_hat_fn(const c1_dofmap& dm, int e, const function_type& u,
                       bpoly<2, T>& dst) {
        grow(gbuf_, dm.local_size());
        detail::gather(dm, e, u.coeffs(), gbuf_.data(), dm.local_size(),
                       typename c1_dofmap::family_tag());
        op_.gather_hat(dm.degree(), gbuf_.data(), dst);
    }
    void begin_matrix(const c1_dofmap& dma, const c1_dofmap& dmb) {
        buf_.clear();
        buf_.reserve(static_cast<std::size_t>(topo_.nt)
                     * static_cast<std::size_t>(dma.local_size())
                     * static_cast<std::size_t>(dmb.local_size()));
    }
    static void grow(std::vector<T>& v, int n) {
        if (v.size() < static_cast<std::size_t>(n))
            v.resize(static_cast<std::size_t>(n));
    }
    void row_into(const typename c1_typed_registry<T>::mat_table& L, int r,
                  const std::vector<T>& b) {
        T acc = L.at(r, 0) * b[0];
        for (int q = 1; q < L.cols(); ++q)
            acc += L.at(r, q) * b[static_cast<std::size_t>(q)];
        gbuf_[static_cast<std::size_t>(r)] = acc;
    }
    // write-once assignment scatter: global = sign x local
    void scatter_assign(const c1_dofmap& dm, int e, function_type& out) {
        for (int r = 0; r < dm.local_size(); ++r) {
            const T& v = gbuf_[static_cast<std::size_t>(r)];
            out.c_(dm.global_dof(e, r), 0) = (dm.dof_sign(e, r) > 0) ? v : -v;
        }
    }
    void scatter_assign_nonvertex(const c1_dofmap& dm, int e, function_type& out) {
        for (int r = 18; r < dm.local_size(); ++r) {
            const T& v = gbuf_[static_cast<std::size_t>(r)];
            out.c_(dm.global_dof(e, r), 0) = (dm.dof_sign(e, r) > 0) ? v : -v;
        }
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_C1_C1_SPACE_HPP
