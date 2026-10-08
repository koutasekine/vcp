// vcp/btfem/st_broken_space.hpp
// btfem BT-2 S2: the space-time broken space
//   Vbar^{q_b}_k (x) broken P_l   (st_broken_space<D,T,P,SP>)
// and its field value type st_broken_field<D,T,P> (coefficient matrix,
// time-block major: rows = time dofs, columns = broken global dofs in the
// bfem broken_dofmap element-block order e nloc + r), plus the dofmap view
// detail::st_broken_dofmap_view that lets BT-1 detail::st_gather read a
// broken block.
//
// Conforms to: BT-2 external design v1.0 (sections 1.3, 3.2, 8) and
//              BT-2 internal design v1.0 (sections 2.1, 2.2, 3.2).
//
// bfem broken_space is used unchanged and is NOT owned (E-3 / D-6). The time
// mesh is copied; time spaces are created per (q, bc) on demand and cached in
// this instance (J-4). Broken dofs carry no sign (+1). Accumulation order of
// inner(): sub-interval l -> element e -> (r, r') row major. No division
// appears in this header.

#ifndef VCP_BTFEM_ST_BROKEN_SPACE_HPP
#define VCP_BTFEM_ST_BROKEN_SPACE_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/element_op.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/detail/scalar_traits.hpp>

#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/time_space.hpp>
#include <vcp/btfem/detail/time_tables.hpp>
#include <vcp/btfem/st_poly.hpp>
#include <vcp/btfem/kron.hpp>

namespace vcp {
namespace btfem {

template <int D, typename T, typename P, class SP> class st_broken_space;

// ---------------------------------------------------------------------------
// st_broken_field<D, T, P>: u(t, x) = sum_a sum_i C[a][i] psi_a(t) q_i(x),
// C is m_t x ndof_broken. Created only through the st_broken_space factories.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class st_broken_field {
public:
    int order() const { return l_; }
    int time_degree() const { return q_; }
    time_bc time_bc_kind() const { return bc_; }
    int num_time_dofs() const { return mt_; }
    int num_space_dofs() const { return nx_; }
    const vcp::matrix<T, P>& coeffs() const { return C_; }
    vcp::matrix<T, P>&       coeffs() { return C_; }

private:
    st_broken_field(int l, int q, time_bc bc, int mt, int nx, const void* sid, const void* tid)
        : l_(l), q_(q), bc_(bc), mt_(mt), nx_(nx), C_(), space_id_(sid), tmesh_id_(tid) {}
    int l_, q_;
    time_bc bc_;
    int mt_, nx_;
    vcp::matrix<T, P> C_;                       // mt_ x nx_
    const void* space_id_;                      // identity of the bfem broken_space
    const void* tmesh_id_;                      // identity of the st_broken_space time_mesh copy
    template <int DD, typename TT, typename PP, class SS> friend class st_broken_space;
};

namespace detail {

// dofmap view of bfem detail::broken_dofmap for BT-1 st_gather, which also
// requires degree() (the broken_dofmap holds no degree; internal 3.2)
struct st_broken_dofmap_view {
    const bfem::detail::broken_dofmap* dm;
    int l;
    st_broken_dofmap_view(const bfem::detail::broken_dofmap& d, int l_) : dm(&d), l(l_) {}
    int degree() const { return l; }
    int local_size() const { return dm->local_size(); }
    int num_elements() const { return dm->num_elements(); }
    int ndof() const { return dm->ndof(); }
    int global_dof(int e, int r) const { return dm->global_dof(e, r); }
    int dof_sign(int, int) const { return 1; }
};

} // namespace detail

// ---------------------------------------------------------------------------
// st_broken_space<D, T, P, SP>
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class st_broken_space {
    static_assert(D == 2 || D == 3, "btfem::st_broken_space: only D == 2 or D == 3");
public:
    typedef st_broken_field<D, T, P>            field_type;
    typedef bfem::broken_space<D, T, P, SP>     space_type;
    typedef kron_term<T, P, SP>                 term_type;
    typedef bfem::broken_field<D, T, P>         section_type;

    // bs: the bfem broken space (not owned; must outlive this object, E-3 / D-6)
    // tm: the time mesh (copied)
    st_broken_space(space_type& bs, const time_mesh<T>& tm)
        : bs_(&bs), tm_(tm), tspaces_(), op_(), u_loc_(), v_loc_() {
        bfem::bfem_scalar_traits<T>::require();
    }

    space_type& space() { return *bs_; }
    const space_type& space() const { return *bs_; }
    const time_mesh<T>& time_mesh_ref() const { return tm_; }
    const time_space<T>& time_dofs(int q, time_bc bc) {
        std::pair<int, int> key(q, static_cast<int>(bc));
        typename std::map<std::pair<int, int>, time_space<T> >::iterator it = tspaces_.find(key);
        if (it != tspaces_.end()) return it->second;
        return tspaces_.insert(std::make_pair(key, time_space<T>(tm_, q, bc))).first->second;
    }
    // dofmap view for BT-1 st_gather (broken order l)
    detail::st_broken_dofmap_view dofmap_view() const {
        return detail::st_broken_dofmap_view(bs_->dofs(), bs_->order());
    }

    // ---- factories and vector conversion (external 3.2) ----
    field_type zero_field(int q, time_bc bc) {
        field_type u = make_field(q, bc);
        u.C_.zeros(u.mt_, u.nx_);
        return u;
    }
    field_type field_from_coeffs(int q, time_bc bc, vcp::matrix<T, P> C) {
        field_type u = make_field(q, bc);
        if (C.rowsize() != u.mt_ || C.columnsize() != u.nx_)
            throw std::invalid_argument("btfem::st_broken_space::field_from_coeffs: size != ndof_t x ndof_broken");
        u.C_ = std::move(C);
        return u;
    }
    field_type field_from_vector(int q, time_bc bc, const vcp::matrix<T, P>& v) {
        field_type u = make_field(q, bc);
        if (v.rowsize() != u.mt_ * u.nx_ || v.columnsize() != 1)
            throw std::invalid_argument("btfem::st_broken_space::field_from_vector: size != ndof_t ndof_broken x 1");
        u.C_.zeros(u.mt_, u.nx_);
        for (int a = 0; a < u.mt_; ++a)
            for (int i = 0; i < u.nx_; ++i)
                u.C_(a, i) = v(a * u.nx_ + i, 0);
        return u;
    }
    vcp::matrix<T, P> to_vector(const field_type& u) const {
        validate(u);
        vcp::matrix<T, P> v;
        v.zeros(u.mt_ * u.nx_, 1);
        for (int a = 0; a < u.mt_; ++a)
            for (int i = 0; i < u.nx_; ++i)
                v(a * u.nx_ + i, 0) = u.C_(a, i);
        return v;
    }

    // ---- sections (BT-1 D-11 semantics) ----
    section_type section_at_node(const field_type& u, int l) {
        validate(u);
        const int m = tm_.num_intervals();
        if (l < 0 || l > m)
            throw std::invalid_argument("btfem::st_broken_space::section_at_node: l out of range");
        const time_space<T>& Vt = time_dofs(u.q_, u.bc_);
        const int gt = (l == 0) ? Vt.global_dof(1, 0) : Vt.global_dof(l, u.q_);
        if (gt < 0) return bs_->zero_field();
        vcp::matrix<T, P> c;
        c.zeros(u.nx_, 1);
        for (int i = 0; i < u.nx_; ++i) c(i, 0) = u.C_(gt, i);
        return bs_->field_from_coeffs(c);
    }
    section_type section_at(const field_type& u, int l, const bfem::bary_point<1, T>& tau) {
        validate(u);
        const int m = tm_.num_intervals();
        if (l < 1 || l > m)
            throw std::invalid_argument("btfem::st_broken_space::section_at: l out of range");
        const time_space<T>& Vt = time_dofs(u.q_, u.bc_);
        vcp::matrix<T, P> c;
        c.zeros(u.nx_, 1);
        std::vector<T> unit(static_cast<std::size_t>(u.q_ + 1), T(0));
        for (int r = 0; r <= u.q_; ++r) {
            const int gt = Vt.global_dof(l, r);
            if (gt < 0) continue;
            unit[static_cast<std::size_t>(r)] = T(1);
            const T w = bfem::eval(bfem::bpoly<1, T>::from_coeffs(u.q_, unit), tau);
            unit[static_cast<std::size_t>(r)] = T(0);
            for (int i = 0; i < u.nx_; ++i) c(i, 0) += w * u.C_(gt, i);
        }
        return bs_->field_from_coeffs(c);
    }

    // ---- (q_j, q_i)_{L2(J;L2)} = U_t(q_b) (x) M_broken (block diagonal) ----
    term_type mass(int q, time_bc bc) {
        return term_type(T(1), time_dofs(q, bc).template mass<P>(), bs_->mass());
    }

    // ---- (u, v)_{L2(J;L2)} of two fields of this space (time spaces may
    //      differ in (q, bc); they share the time mesh by construction) ----
    // sum_l sum_e sum_{r, r'} k Mref(q_u, q_v)(r, r') * |e| (u_r, v_r')_ref
    T inner(const field_type& u, const field_type& v) {
        validate(u);
        validate(v);
        const time_space<T>& Vtu = time_dofs(u.q_, u.bc_);
        const time_space<T>& Vtv = time_dofs(v.q_, v.bc_);
        const typename time_space<T>::table_type& M =
            detail::time_tables_typed<T>::mass_ref(u.q_, v.q_);
        const T& k = tm_.width();
        const detail::st_broken_dofmap_view dm = dofmap_view();
        const int m = tm_.num_intervals();
        const int nt = bs_->num_elements();
        T acc(0);
        for (int l = 1; l <= m; ++l) {
            for (int e = 0; e < nt; ++e) {
                op_.set_geometry(bs_->geometry(e));
                detail::st_gather(u_loc_, u.C_, Vtu, l, dm, e);
                detail::st_gather(v_loc_, v.C_, Vtv, l, dm, e);
                for (int r = 0; r <= u.q_; ++r) {
                    for (int rp = 0; rp <= v.q_; ++rp) {
                        const T tv = k * M(r, rp);
                        acc += tv * op_.local_inner(u_loc_.r(r), v_loc_.r(rp));
                    }
                }
            }
        }
        return acc;
    }

    // membership / shape check (public: used by the st_rt_assemble free functions)
    void validate(const field_type& u) const {
        if (u.space_id_ != static_cast<const void*>(bs_) || u.tmesh_id_ != static_cast<const void*>(&tm_))
            throw std::invalid_argument("btfem::st_broken_space: st_broken_field does not belong to this space");
        if (u.l_ != bs_->order())
            throw std::invalid_argument("btfem::st_broken_space: st_broken_field order != broken_space::order()");
        if (u.C_.rowsize() != u.mt_ || u.C_.columnsize() != u.nx_)
            throw std::invalid_argument("btfem::st_broken_space: st_broken_field coefficient matrix size changed");
    }

private:
    space_type* bs_;
    time_mesh<T> tm_;
    std::map<std::pair<int, int>, time_space<T> > tspaces_;
    bfem::element_op<D, T, P> op_;
    detail::st_local<D, T> u_loc_, v_loc_;

    field_type make_field(int q, time_bc bc) {
        const time_space<T>& Vt = time_dofs(q, bc);
        return field_type(bs_->order(), q, bc, Vt.ndof(), bs_->ndof(),
                          static_cast<const void*>(bs_), static_cast<const void*>(&tm_));
    }
};

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_ST_BROKEN_SPACE_HPP
