// vcp/btfem/st_rt_space.hpp
// btfem BT-2 S1: the space-time Raviart-Thomas space
//   Sigma_h^k = Vbar^{q_sigma}_k (x) RT^l   (st_rt_space<D,T,P,SP>)
// and its field value type st_rt_field<D,T,P> (coefficient matrix, time-block
// major: rows = time dofs, columns = RT global dofs in the bfem rt_dofmap
// convention), the signed local gather detail::st_gather_rt and the local
// time-rank list detail::st_rt_local.
//
// Conforms to: BT-2 external design v1.0 (sections 1.3, 3.1, 3.3 (div is in
//              st_rt_assemble.hpp), 8) and
//              BT-2 internal design v1.0 (sections 2.1, 2.2, 3.1).
//
// bfem rt_space is used unchanged and is NOT owned (non-const reference,
// E-3; it must outlive this object, the D-6 lifetime rule of pg_space). The
// time mesh is copied. The time space (q, bc) is never stored in the space:
// every field / matrix is created with an explicit (q, bc) (R10) and the
// time_space objects are cached per (q, bc) in this instance (J-4).
// Sign rule (external 1.3 principle 5): RT dofs carry the bfem rt_dofmap sign
// (edge / face orientation); the time direction is a tensor product, so the
// same sign multiplies every time rank. No division appears in this header.

#ifndef VCP_BTFEM_ST_RT_SPACE_HPP
#define VCP_BTFEM_ST_RT_SPACE_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/rt_element_op.hpp>
#include <vcp/bfem/detail/scalar_traits.hpp>

#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/time_space.hpp>
#include <vcp/btfem/kron.hpp>

namespace vcp {
namespace btfem {

template <int D, typename T, typename P, class SP> class st_rt_space;

// ---------------------------------------------------------------------------
// st_rt_field<D, T, P>: sigma(t, x) = sum_a sum_j C[a][j] psi_a(t) tau_j(x),
// C is m_t x ndof_rt (time-block major flattening v[a ndof_rt + j]).
// Created only through the st_rt_space factories (internal 2.1).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class st_rt_field {
public:
    int order() const { return k_; }
    int time_degree() const { return q_; }
    time_bc time_bc_kind() const { return bc_; }
    int num_time_dofs() const { return mt_; }
    int num_space_dofs() const { return nx_; }
    const vcp::matrix<T, P>& coeffs() const { return C_; }
    vcp::matrix<T, P>&       coeffs() { return C_; }

private:
    st_rt_field(int k, int q, time_bc bc, int mt, int nx, const void* sid, const void* tid)
        : k_(k), q_(q), bc_(bc), mt_(mt), nx_(nx), C_(), space_id_(sid), tmesh_id_(tid) {}
    int k_, q_;
    time_bc bc_;
    int mt_, nx_;
    vcp::matrix<T, P> C_;                       // mt_ x nx_
    const void* space_id_;                      // identity of the bfem rt_space
    const void* tmesh_id_;                      // identity of the st_rt_space time_mesh copy
    template <int DD, typename TT, typename PP, class SS> friend class st_rt_space;
};

namespace detail {

// ---------------------------------------------------------------------------
// st_rt_local<T>: local RT coefficient rows of one sub-interval x element,
// row[r] (r = 0..qt) in the bfem rt_local_coeffs convention (reference DOF
// basis, canonical order, signs already applied by st_gather_rt).
// ---------------------------------------------------------------------------
template <typename T>
struct st_rt_local {
    int qt;                                             // time degree
    int k;                                              // RT order
    std::vector<bfem::rt_local_coeffs<T> > row;         // row[r], r = 0..qt

    st_rt_local() : qt(0), k(0), row(1) { row[0].k = 0; }

    void prepare(int qt_, int k_, int dim) {
        qt = qt_;
        k = k_;
        if (row.size() < static_cast<std::size_t>(qt_ + 1))
            row.resize(static_cast<std::size_t>(qt_ + 1));
        for (int r = 0; r <= qt_; ++r) {
            bfem::rt_local_coeffs<T>& s = row[static_cast<std::size_t>(r)];
            s.k = k_;
            if (s.c.size() != static_cast<std::size_t>(dim))
                s.c.assign(static_cast<std::size_t>(dim), T(0));
        }
    }
    bfem::rt_local_coeffs<T>& r(int i) {
        assert(i >= 0 && i <= qt);
        return row[static_cast<std::size_t>(i)];
    }
    const bfem::rt_local_coeffs<T>& r(int i) const {
        assert(i >= 0 && i <= qt);
        return row[static_cast<std::size_t>(i)];
    }
};

// signed gather (internal 3.1): row[r].c[j] = sign(e, j) * C(g_t(l, r), g_rt(e, j));
// a removed time dof (g_t < 0) reads as T(0). The same sign multiplies every
// time rank (external 1.3 principle 5). The multiplication order
// T(sign) * value is that of bfem detail::gather_general (Q-R1).
template <typename T, typename P, typename DM>
void st_gather_rt(st_rt_local<T>& out, const vcp::matrix<T, P>& C, const time_space<T>& Vt,
                  int l, const DM& dm, int e, int k) {
    const int qt = Vt.degree();
    const int dim = dm.local_size();
    out.prepare(qt, k, dim);
    for (int r = 0; r <= qt; ++r) {
        const int gt = Vt.global_dof(l, r);
        std::vector<T>& c = out.row[static_cast<std::size_t>(r)].c;
        if (gt < 0) {
            const T z(0);
            for (int j = 0; j < dim; ++j) c[static_cast<std::size_t>(j)] = z;
        } else {
            for (int j = 0; j < dim; ++j)
                c[static_cast<std::size_t>(j)] = T(dm.dof_sign(e, j)) * C(gt, dm.global_dof(e, j));
        }
    }
}

} // namespace detail

// ---------------------------------------------------------------------------
// st_rt_space<D, T, P, SP>
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class st_rt_space {
    static_assert(D == 2 || D == 3, "btfem::st_rt_space: only D == 2 or D == 3");
public:
    typedef st_rt_field<D, T, P>            field_type;
    typedef bfem::rt_space<D, T, P, SP>     space_type;
    typedef kron_term<T, P, SP>             term_type;
    typedef bfem::rt_field<D, T, P>         section_type;

    // rs: the bfem RT space (not owned; must outlive this object, E-3 / D-6)
    // tm: the time mesh (copied)
    st_rt_space(space_type& rs, const time_mesh<T>& tm)
        : rs_(&rs), tm_(tm), tspaces_() {
        bfem::bfem_scalar_traits<T>::require();
    }

    space_type& space() { return *rs_; }
    const space_type& space() const { return *rs_; }
    const time_mesh<T>& time_mesh_ref() const { return tm_; }
    // (q, bc) -> time_space, cached in this instance (J-4; same logic as pg_space)
    const time_space<T>& time_dofs(int q, time_bc bc) {
        std::pair<int, int> key(q, static_cast<int>(bc));
        typename std::map<std::pair<int, int>, time_space<T> >::iterator it = tspaces_.find(key);
        if (it != tspaces_.end()) return it->second;
        return tspaces_.insert(std::make_pair(key, time_space<T>(tm_, q, bc))).first->second;
    }

    // ---- factories and vector conversion (external 3.1) ----
    field_type zero_field(int q, time_bc bc) {
        field_type s = make_field(q, bc);
        s.C_.zeros(s.mt_, s.nx_);
        return s;
    }
    field_type field_from_coeffs(int q, time_bc bc, vcp::matrix<T, P> C) {
        field_type s = make_field(q, bc);
        if (C.rowsize() != s.mt_ || C.columnsize() != s.nx_)
            throw std::invalid_argument("btfem::st_rt_space::field_from_coeffs: size != ndof_t x ndof_rt");
        s.C_ = std::move(C);
        return s;
    }
    field_type field_from_vector(int q, time_bc bc, const vcp::matrix<T, P>& v) {
        field_type s = make_field(q, bc);
        if (v.rowsize() != s.mt_ * s.nx_ || v.columnsize() != 1)
            throw std::invalid_argument("btfem::st_rt_space::field_from_vector: size != ndof_t ndof_rt x 1");
        s.C_.zeros(s.mt_, s.nx_);
        for (int a = 0; a < s.mt_; ++a)
            for (int i = 0; i < s.nx_; ++i)
                s.C_(a, i) = v(a * s.nx_ + i, 0);
        return s;
    }
    // time-block major (a, i) -> a ndof_rt + i (signs unchanged)
    vcp::matrix<T, P> to_vector(const field_type& s) const {
        validate(s);
        vcp::matrix<T, P> v;
        v.zeros(s.mt_ * s.nx_, 1);
        for (int a = 0; a < s.mt_; ++a)
            for (int i = 0; i < s.nx_; ++i)
                v(a * s.nx_ + i, 0) = s.C_(a, i);
        return v;
    }

    // ---- sections (BT-1 external 5.2 / D-11 semantics) ----
    // RT field at the node tau_l (l = 0..m): block row copy, no evaluation;
    // bc_zero_left and l = 0 give the zero field
    section_type section_at_node(const field_type& s, int l) {
        validate(s);
        const int m = tm_.num_intervals();
        if (l < 0 || l > m)
            throw std::invalid_argument("btfem::st_rt_space::section_at_node: l out of range");
        const time_space<T>& Vt = time_dofs(s.q_, s.bc_);
        const int gt = (l == 0) ? Vt.global_dof(1, 0) : Vt.global_dof(l, s.q_);
        if (gt < 0) return rs_->zero_field();
        vcp::matrix<T, P> c;
        c.zeros(s.nx_, 1);
        for (int i = 0; i < s.nx_; ++i) c(i, 0) = s.C_(gt, i);
        return rs_->field_from_coeffs(c);
    }
    // RT field inside sub-interval l (1..m) at local barycentric tau = (1 - t, t):
    // sum_r B^q_r(tau) C[g(l, r), .] (de Casteljau weights, r ascending)
    section_type section_at(const field_type& s, int l, const bfem::bary_point<1, T>& tau) {
        validate(s);
        const int m = tm_.num_intervals();
        if (l < 1 || l > m)
            throw std::invalid_argument("btfem::st_rt_space::section_at: l out of range");
        const time_space<T>& Vt = time_dofs(s.q_, s.bc_);
        vcp::matrix<T, P> c;
        c.zeros(s.nx_, 1);
        std::vector<T> unit(static_cast<std::size_t>(s.q_ + 1), T(0));
        for (int r = 0; r <= s.q_; ++r) {
            const int gt = Vt.global_dof(l, r);
            if (gt < 0) continue;
            unit[static_cast<std::size_t>(r)] = T(1);
            const T w = bfem::eval(bfem::bpoly<1, T>::from_coeffs(s.q_, unit), tau);
            unit[static_cast<std::size_t>(r)] = T(0);
            for (int i = 0; i < s.nx_; ++i) c(i, 0) += w * s.C_(gt, i);
        }
        return rs_->field_from_coeffs(c);
    }

    // ---- S3: (sigma_j, tau_i)_{L2(J;L2)} = U_t(q_sigma) (x) M_RT ----
    term_type mass(int q, time_bc bc) {
        return term_type(T(1), time_dofs(q, bc).template mass<P>(), rs_->mass());
    }

    // membership / shape check of a field (public so that the free functions
    // of st_rt_assemble.hpp can use it; std::invalid_argument on violation)
    void validate(const field_type& s) const {
        if (s.space_id_ != static_cast<const void*>(rs_) || s.tmesh_id_ != static_cast<const void*>(&tm_))
            throw std::invalid_argument("btfem::st_rt_space: st_rt_field does not belong to this space");
        if (s.k_ != rs_->order())
            throw std::invalid_argument("btfem::st_rt_space: st_rt_field RT order != rt_space::order()");
        if (s.C_.rowsize() != s.mt_ || s.C_.columnsize() != s.nx_)
            throw std::invalid_argument("btfem::st_rt_space: st_rt_field coefficient matrix size changed");
    }

private:
    space_type* rs_;
    time_mesh<T> tm_;
    std::map<std::pair<int, int>, time_space<T> > tspaces_;

    field_type make_field(int q, time_bc bc) {
        const time_space<T>& Vt = time_dofs(q, bc);
        return field_type(rs_->order(), q, bc, Vt.ndof(), rs_->ndof(),
                          static_cast<const void*>(rs_), static_cast<const void*>(&tm_));
    }
};

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_ST_RT_SPACE_HPP
