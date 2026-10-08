// vcp/btfem/st_poly.hpp
// btfem T4: the space-time function value type st_poly<D,T,P> (coefficient
// matrix, time-block major) and the local tensor polynomial machinery of the
// sub-interval x element loops (detail::st_local, st_workspace, st_gather).
//
// Conforms to: BT-1 external design v1.1 (section 5.1) and
//              BT-1 internal design v1.0 (sections 5.1, 5.2, 5.3).
//
// A local tensor polynomial u(tau, lambda) = sum_r sum_alpha c_{r,alpha}
// B^{qt}_r(tau) B^n_alpha(lambda) is stored as the list of its time-rank rows,
// each a bfem bpoly<D,T> of degree n. All operations are *_into forms with
// caller-owned outputs and workspace (bfem L1 contract); loops run in fixed
// order (time rank r, then r') so that the accumulation order is
// deterministic. No division appears in this header.

#ifndef VCP_BTFEM_ST_POLY_HPP
#define VCP_BTFEM_ST_POLY_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/refine.hpp>
#include <vcp/bfem/element_op.hpp>
#include <vcp/btfem/time_space.hpp>
#include <vcp/btfem/detail/time_tables.hpp>

namespace vcp {
namespace btfem {

template <int D, typename T, typename P, class SP> class pg_space;

// ---------------------------------------------------------------------------
// st_poly<D, T, P>: u(t, x) = sum_a sum_i C[a][i] psi_a(t) phi_i(x) in
// V^q_k (or Vbar^q_k) x S_h^n. Rows = time dofs a, columns = space dofs i;
// row-major flattening equals the time-block-major vector v[a n_x + i].
// Created only through the pg_space factories.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
class st_poly {
public:
    int time_degree() const { return q_; }
    time_bc time_bc_kind() const { return bc_; }
    int space_degree() const { return n_; }
    int num_time_dofs() const { return mt_; }
    int num_space_dofs() const { return nx_; }
    const vcp::matrix<T, P>& coeffs() const { return C_; }
    vcp::matrix<T, P>& coeffs() { return C_; }

private:
    st_poly(int q, time_bc bc, int n, int mt, int nx, const void* sid, const void* tid)
        : q_(q), n_(n), bc_(bc), mt_(mt), nx_(nx), C_(), space_id_(sid), tmesh_id_(tid) {}
    int q_, n_;
    time_bc bc_;
    int mt_, nx_;
    vcp::matrix<T, P> C_;                       // mt_ x nx_
    const void* space_id_;                      // identity of the owning fe_space
    const void* tmesh_id_;                      // identity of the pg_space time_mesh copy
    template <int DD, typename TT, typename PP, class SS> friend class pg_space;
};

namespace detail {

// ---------------------------------------------------------------------------
// st_local<D, T>: local tensor polynomial on one sub-interval x element.
// row is grow-only; rows [0, qt] are the live ones.
// ---------------------------------------------------------------------------
template <int D, typename T>
struct st_local {
    int qt;                                     // time degree
    int n;                                      // space degree
    std::vector<bfem::bpoly<D, T> > row;        // row[r], r = 0..qt (degree n each)

    st_local() : qt(0), n(0), row(1) {}

    // size the rows (no shrink; reuse of bpoly storage via bpoly_access)
    void prepare(int qt_, int n_, bool zero_fill) {
        qt = qt_;
        n = n_;
        if (row.size() < static_cast<std::size_t>(qt_ + 1))
            row.resize(static_cast<std::size_t>(qt_ + 1));
        for (int r = 0; r <= qt_; ++r)
            bfem::detail::bpoly_access::prepare(row[static_cast<std::size_t>(r)], n_, zero_fill);
    }
    bfem::bpoly<D, T>& r(int i) {
        assert(i >= 0 && i <= qt);
        return row[static_cast<std::size_t>(i)];
    }
    const bfem::bpoly<D, T>& r(int i) const {
        assert(i >= 0 && i <= qt);
        return row[static_cast<std::size_t>(i)];
    }
    void add_scalar(const T& s) {               // partition of unity in tau and lambda
        for (int i = 0; i <= qt; ++i) row[static_cast<std::size_t>(i)].add_scalar(s);
    }
};

template <int D, typename T>
struct st_workspace {
    st_local<D, T> acc, tmp, tmp2, g;           // Horner ping-pong / gradient scratch
    bfem::bpoly<D, T> sb0, sb1, sb2;            // bpoly scratch
    std::vector<T> tvals;                       // eval scratch (time direction)
};

// copy (same shape, assignment only)
template <int D, typename T>
void st_copy(st_local<D, T>& dst, const st_local<D, T>& src) {
    dst.prepare(src.qt, src.n, false);
    for (int r = 0; r <= src.qt; ++r) {
        std::vector<T>& d = bfem::detail::bpoly_access::vec(dst.r(r));
        const std::vector<T>& s = bfem::detail::bpoly_access::vec(src.r(r));
        for (std::size_t i = 0; i < d.size(); ++i) d[i] = s[i];
    }
}

// out = u * v: time product coefficients from typed_registry<1,T>::product
// (target rank r + r'), space product by bfem mul_into. Loop order r -> r'.
template <int D, typename T>
void st_mul_into(st_local<D, T>& out, const st_local<D, T>& u, const st_local<D, T>& v,
                 st_workspace<D, T>& ws) {
    assert(&out != &u && &out != &v);
    const bfem::typed_product_table<1, T>& Pt =
        bfem::typed_registry<1, T>::product(u.qt, v.qt);
    out.prepare(u.qt + v.qt, u.n + v.n, true);
    for (int r = 0; r <= u.qt; ++r) {
        for (int rp = 0; rp <= v.qt; ++rp) {
            bfem::mul_into(ws.sb0, u.r(r), v.r(rp));                 // degree nu + nv
            bfem::scale_into(ws.sb1, ws.sb0, Pt.coeff(r, rp));
            out.r(Pt.target_rank(r, rp)) += ws.sb1;                  // same degree
        }
    }
}

// out = u + v (equal time degree required; space degrees lifted by bfem add_into)
template <int D, typename T>
void st_add_into(st_local<D, T>& out, const st_local<D, T>& u, const st_local<D, T>& v) {
    assert(&out != &u && &out != &v);
    if (u.qt != v.qt)
        throw std::invalid_argument("btfem::st_add_into: time degree mismatch");
    const int n = u.n >= v.n ? u.n : v.n;
    out.prepare(u.qt, n, false);
    for (int r = 0; r <= u.qt; ++r) bfem::add_into(out.r(r), u.r(r), v.r(r));
}

// out = f(u) by Horner on the tensor Bernstein algebra: acc = a_K;
// acc = acc * u + a_k (k = K-1 .. 0). Degrees grow to (K qu, K nu).
template <int D, typename T>
void st_compose_into(st_local<D, T>& out, const bfem::poly1<T>& f, const st_local<D, T>& u,
                     st_workspace<D, T>& ws) {
    assert(&out != &u);
    const int K = f.degree();
    if (K == 0) {
        out.prepare(0, 0, false);
        bfem::detail::bpoly_access::vec(out.r(0))[0] = f.coeff(0);
        return;
    }
    st_local<D, T>* r = &ws.acc;
    st_local<D, T>* t = &ws.tmp;
    r->prepare(0, 0, false);
    bfem::detail::bpoly_access::vec(r->r(0))[0] = f.coeff(K);
    for (int k = K - 1; k >= 0; --k) {
        st_mul_into(*t, *r, u, ws);
        t->add_scalar(f.coeff(k));
        st_local<D, T>* sw = r;
        r = t;
        t = sw;
    }
    st_copy(out, *r);
}

// time derivative: out.row[rp] = (qt * inv_k) * (u.row[rp+1] - u.row[rp]),
// time degree qt - 1 (qt == 0 gives the zero of time degree 0)
template <int D, typename T>
void st_dt_into(st_local<D, T>& out, const st_local<D, T>& u, const T& inv_k,
                st_workspace<D, T>& ws) {
    assert(&out != &u);
    if (u.qt == 0) {
        out.prepare(0, u.n, true);
        return;
    }
    const deriv_matrix& Dq = time_tables_rational::deriv(u.qt);
    const T fac = T(Dq.q) * inv_k;
    out.prepare(u.qt - 1, u.n, false);
    for (int rp = 0; rp < u.qt; ++rp) {
        // D[rp][rp+1] = +q, D[rp][rp] = -q (verified structure of derivative_map<1>)
        assert(Dq.at(rp, rp + 1) == Dq.q && Dq.at(rp, rp) == -Dq.q);
        bfem::sub_into(ws.sb0, u.r(rp + 1), u.r(rp));
        bfem::scale_into(out.r(rp), ws.sb0, fac);
    }
}

// spatial gradient component d: each row through element_op::grad_component
template <int D, typename T, typename P>
void st_grad_into(st_local<D, T>& out, const st_local<D, T>& u, int d,
                  bfem::element_op<D, T, P>& op) {
    assert(&out != &u);
    out.prepare(u.qt, u.n >= 1 ? u.n - 1 : 0, false);
    for (int r = 0; r <= u.qt; ++r) op.grad_component(u.r(r), d, out.r(r));
}

// out = sum_d b_d * (grad u)_d  (b points to D components of the same time
// degree; a pointer avoids the int / size_t deduction clash of std::array<_, D>)
template <int D, typename T, typename P>
void st_b_dot_grad_into(st_local<D, T>& out, const st_local<D, T>* b,
                        const st_local<D, T>& u, st_workspace<D, T>& ws,
                        bfem::element_op<D, T, P>& op) {
    assert(&out != &u);
    for (int d = 0; d < D; ++d) {
        st_grad_into(ws.g, u, d, op);
        if (d == 0) {
            st_mul_into(out, b[0], ws.g, ws);
        } else {
            st_mul_into(ws.tmp2, b[d], ws.g, ws);
            st_add_into(ws.acc, out, ws.tmp2);
            st_copy(out, ws.acc);
        }
    }
}

// point evaluation: space by bfem eval per row, then 1D de Casteljau in tau
// (same operation order as bfem::eval for D = 1)
template <int D, typename T>
T st_eval(const st_local<D, T>& u, const bfem::bary_point<1, T>& tau,
          const bfem::bary_point<D, T>& lam, st_workspace<D, T>& ws) {
    if (ws.tvals.size() < static_cast<std::size_t>(u.qt + 1))
        ws.tvals.resize(static_cast<std::size_t>(u.qt + 1));
    for (int r = 0; r <= u.qt; ++r)
        ws.tvals[static_cast<std::size_t>(r)] = bfem::eval(u.r(r), lam);
    for (int k = u.qt; k >= 1; --k) {
        for (int r = 0; r < k; ++r) {
            T acc(0);
            acc += tau[0] * ws.tvals[static_cast<std::size_t>(r)];
            acc += tau[1] * ws.tvals[static_cast<std::size_t>(r + 1)];
            ws.tvals[static_cast<std::size_t>(r)] = acc;
        }
    }
    return ws.tvals[0];
}

// coefficient range merge over all rows (bound_traits)
template <int D, typename T>
void st_range_merge(bfem::range_pair<T>& acc, bool& started, const st_local<D, T>& u) {
    for (int r = 0; r <= u.qt; ++r) bfem::detail::range_merge(acc, started, u.r(r));
}

// gather: coefficient matrix C (time dofs x space dofs) -> local tensor on
// sub-interval l, element e. Removed time dof (g < 0) reads as T(0).
template <int D, typename T, typename P, typename DM>
void st_gather(st_local<D, T>& out, const vcp::matrix<T, P>& C, const time_space<T>& Vt,
               int l, const DM& dm, int e) {
    const int qt = Vt.degree();
    const int N = dm.local_size();
    out.prepare(qt, dm.degree(), false);
    for (int r = 0; r <= qt; ++r) {
        int gt = Vt.global_dof(l, r);
        std::vector<T>& c = bfem::detail::bpoly_access::vec(out.r(r));
        if (gt < 0) {
            const T z(0);
            for (int a = 0; a < N; ++a) c[static_cast<std::size_t>(a)] = z;
        } else {
            for (int a = 0; a < N; ++a) c[static_cast<std::size_t>(a)] = C(gt, dm.global_dof(e, a));
        }
    }
}

} // namespace detail
} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_ST_POLY_HPP
