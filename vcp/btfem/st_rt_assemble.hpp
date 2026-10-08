// vcp/btfem/st_rt_assemble.hpp
// btfem BT-2 S3-S11: the flux-system free functions of the space-time RT /
// broken spaces -- the time versions of the seven bfem rt_assemble.hpp
// functions (K3-K9) and the three squared norms -- plus the BT-2 detail
// parts that BT-1's st_poly.hpp does not provide (J-1):
//   detail::st_sub_into             row-wise bpoly subtraction (equal time degree)
//   detail::st_elevate_time_into    time-degree elevation of an st_local (L0 D = 1 table)
//   detail::st_elevate_time_rt_into the same elevation of an st_rt_local
//
// Conforms to: BT-2 external design v1.0 (sections 1.4, 3.3, 4, 5, 6, 8) and
//              BT-2 internal design v1.0 (sections 4, 5, 6, 7, 8).
//
// Accumulation order (normative): sub-interval l -> element e -> (r, r') row
// major for the quadratic forms; l -> e -> s -> a -> i for the vectors. Same
// input, single thread: bit reproducible. No parallel region and no
// allocation-free contract in this header (J-3; bfem K5-K8 model). No
// division appears in this header.

#ifndef VCP_BTFEM_ST_RT_ASSEMBLE_HPP
#define VCP_BTFEM_ST_RT_ASSEMBLE_HPP

#include <vector>
#include <array>
#include <utility>
#include <string>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/typed_tables.hpp>
#include <vcp/bfem/element_op.hpp>
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_element_op.hpp>
#include <vcp/bfem/rt/rt_typed_tables.hpp>
#include <vcp/bfem/rt/rt_assemble.hpp>

#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/time_space.hpp>
#include <vcp/btfem/detail/time_tables.hpp>
#include <vcp/btfem/st_poly.hpp>
#include <vcp/btfem/kron.hpp>
#include <vcp/btfem/pg_space.hpp>
#include <vcp/btfem/st_rt_space.hpp>
#include <vcp/btfem/st_broken_space.hpp>

namespace vcp {
namespace btfem {
namespace detail {

// ---------------------------------------------------------------------------
// out = u - v (equal time degree required; space degrees lifted by bfem
// sub_into) -- the signed counterpart of BT-1 st_add_into (internal 4.1, J-1)
// ---------------------------------------------------------------------------
template <int D, typename T>
void st_sub_into(st_local<D, T>& out, const st_local<D, T>& u, const st_local<D, T>& v) {
    assert(&out != &u && &out != &v);
    if (u.qt != v.qt)
        throw std::invalid_argument("btfem::st_sub_into: time degree mismatch");
    const int n = u.n >= v.n ? u.n : v.n;
    out.prepare(u.qt, n, false);
    for (int r = 0; r <= u.qt; ++r) bfem::sub_into(out.r(r), u.r(r), v.r(r));
}

// ---------------------------------------------------------------------------
// time-degree elevation (internal 4.2): out.row[r'] = sum_r E[r'][r] u.row[r]
// with E = typed_registry<1,T>::elevation(q_from, q_to) (the L0 elevation
// table for D = 1; enclose-once T image of the rational table). The sparse
// table is stored per source rank r; iterating r ascending and accumulating
// into the targets adds, for every target r', its terms in ascending r --
// the same order as the dense r'-outer / r-inner formulation. Space degree
// and the point values in tau are unchanged (Q-R3). q_to == q_from copies.
// ---------------------------------------------------------------------------
template <int D, typename T>
void st_elevate_time_into(st_local<D, T>& out, const st_local<D, T>& u, int q_to,
                          st_workspace<D, T>& ws) {
    assert(&out != &u);
    if (q_to < u.qt)
        throw std::invalid_argument("btfem::st_elevate_time_into: q_to < time degree");
    if (q_to == u.qt) {
        st_copy(out, u);
        return;
    }
    const bfem::typed_elevation_table<1, T>& E = bfem::typed_registry<1, T>::elevation(u.qt, q_to);
    out.prepare(q_to, u.n, true);
    for (int r = 0; r <= u.qt; ++r) {
        typename bfem::typed_elevation_table<1, T>::entry_range rr = E.row(r);
        for (typename bfem::typed_elevation_table<1, T>::entry_iterator p = rr.begin(); p != rr.end(); ++p) {
            typename bfem::typed_elevation_table<1, T>::entry en = *p;
            bfem::scale_into(ws.sb0, u.r(r), *en.coeff);
            out.r(en.target_rank) += ws.sb0;                     // same space degree
        }
    }
}

// the same elevation of the RT coefficient rows (coefficient-wise linear
// combination with the same E; RT order and dim unchanged)
template <typename T>
void st_elevate_time_rt_into(st_rt_local<T>& out, const st_rt_local<T>& u, int q_to) {
    assert(&out != &u);
    if (q_to < u.qt)
        throw std::invalid_argument("btfem::st_elevate_time_rt_into: q_to < time degree");
    const int dim = static_cast<int>(u.row[0].c.size());
    if (q_to == u.qt) {
        out.prepare(u.qt, u.k, dim);
        for (int r = 0; r <= u.qt; ++r)
            for (int j = 0; j < dim; ++j)
                out.row[static_cast<std::size_t>(r)].c[static_cast<std::size_t>(j)] =
                    u.row[static_cast<std::size_t>(r)].c[static_cast<std::size_t>(j)];
        return;
    }
    const bfem::typed_elevation_table<1, T>& E = bfem::typed_registry<1, T>::elevation(u.qt, q_to);
    out.prepare(q_to, u.k, dim);
    const T z(0);
    for (int rp = 0; rp <= q_to; ++rp)
        for (int j = 0; j < dim; ++j)
            out.row[static_cast<std::size_t>(rp)].c[static_cast<std::size_t>(j)] = z;
    for (int r = 0; r <= u.qt; ++r) {
        const std::vector<T>& src = u.row[static_cast<std::size_t>(r)].c;
        typename bfem::typed_elevation_table<1, T>::entry_range rr = E.row(r);
        for (typename bfem::typed_elevation_table<1, T>::entry_iterator p = rr.begin(); p != rr.end(); ++p) {
            typename bfem::typed_elevation_table<1, T>::entry en = *p;
            std::vector<T>& dst = out.row[static_cast<std::size_t>(en.target_rank)].c;
            for (int j = 0; j < dim; ++j)
                dst[static_cast<std::size_t>(j)] += (*en.coeff) * src[static_cast<std::size_t>(j)];
        }
    }
}

inline void st_require_same_time_mesh(bool same, const char* where) {
    if (!same) {
        std::string msg("btfem::");
        msg += where;
        msg += ": spaces live on different time meshes";
        throw std::invalid_argument(msg);
    }
}

} // namespace detail

// ===========================================================================
// S1 operation: div sigma -> broken (x) time (external 3.3, internal 7)
// ===========================================================================
// div sigma as a field of sbs with the time space (q_sigma, bc_sigma) of
// sigma; requires sbs.space().order() == srs.space().order() (RT^k has
// div in broken P_k). Time rank by time rank through rt_element_op::
// div_as_bpoly; broken dofs are element private (write-once, no clash); a
// time dof shared by two sub-intervals receives the same value twice.
template <int D, typename T, typename P, class SPR, class SPB>
typename st_broken_space<D, T, P, SPB>::field_type
st_div_field(st_rt_space<D, T, P, SPR>& srs, const typename st_rt_space<D, T, P, SPR>::field_type& sig,
             st_broken_space<D, T, P, SPB>& sbs) {
    srs.validate(sig);
    detail::st_require_same_time_mesh(srs.time_mesh_ref().same_as(sbs.time_mesh_ref()), "st_div_field");
    bfem::rt_space<D, T, P, SPR>& rs = srs.space();
    bfem::broken_space<D, T, P, SPB>& bs = sbs.space();
    bfem::detail::rt_require_same_mesh(bs.num_elements(), bs.num_vertices(),
                                       rs.num_elements(), rs.num_vertices(), "st_div_field");
    if (bs.order() != rs.order())
        throw std::invalid_argument("btfem::st_div_field: broken order != RT order");
    const int k = rs.order();
    const int nloc = bs.local_size();
    const time_space<T>& Vt = srs.time_dofs(sig.time_degree(), sig.time_bc_kind());
    typename st_broken_space<D, T, P, SPB>::field_type out = sbs.zero_field(sig.time_degree(), sig.time_bc_kind());
    bfem::rt_element_op<D, T, P> op;
    detail::st_rt_local<T> sloc;
    bfem::bpoly<D, T> d;
    const int m = srs.time_mesh_ref().num_intervals();
    const int nt = rs.num_elements();
    for (int l = 1; l <= m; ++l) {
        for (int e = 0; e < nt; ++e) {
            op.set_geometry(rs.geometry(e));
            detail::st_gather_rt(sloc, sig.coeffs(), Vt, l, rs.dofs(), e, k);
            for (int r = 0; r <= sig.time_degree(); ++r) {
                const int gt = Vt.global_dof(l, r);
                if (gt < 0) continue;
                op.div_as_bpoly(sloc.r(r), d);
                for (int i = 0; i < nloc; ++i) out.coeffs()(gt, e * nloc + i) = d.coeff(i);
            }
        }
    }
    return out;
}

// ===========================================================================
// S4-S6: time-independent bilinear forms (external 4, internal 5) -- one
// kron_term each: scale T(1), time = mixed time mass (rows = row space,
// cols = column space), space = the bfem assembler. The same-mesh checks are
// those of the bfem assembler (space) and time_space::mixed_mass (time).
// ===========================================================================
// S4 (K3): rows = broken (x) (qb, bcb), cols = RT (x) (qs, bcs); (div sigma_j, q_i)
template <int D, typename T, typename P, class SPB, class SPR>
kron_term<T, P, SPB> st_assemble_div_mass(st_broken_space<D, T, P, SPB>& sbs, int qb, time_bc bcb,
                                          st_rt_space<D, T, P, SPR>& srs, int qs, time_bc bcs) {
    detail::st_require_same_time_mesh(sbs.time_mesh_ref().same_as(srs.time_mesh_ref()), "st_assemble_div_mass");
    return kron_term<T, P, SPB>(T(1), sbs.time_dofs(qb, bcb).template mixed_mass<P>(srs.time_dofs(qs, bcs)),
                                bfem::assemble_div_mass(sbs.space(), srs.space()));
}

// S5 (K4): rows = RT (x) (qs, bcs), cols = V_k^q (x) S_h^n (pg_space trial side,
// bc_zero_left fixed -- E-5); (sigma_i, grad phi_j)
template <int D, typename T, typename P, class SPR, class SPF>
kron_term<T, P, SPR> st_assemble_cross_grad(st_rt_space<D, T, P, SPR>& srs, int qs, time_bc bcs,
                                            pg_space<D, T, P, SPF>& pg, int q, int n) {
    detail::st_require_same_time_mesh(srs.time_mesh_ref().same_as(pg.time_mesh_ref()), "st_assemble_cross_grad");
    return kron_term<T, P, SPR>(T(1), srs.time_dofs(qs, bcs).template mixed_mass<P>(pg.time_dofs(q, time_bc_zero_left)),
                                bfem::assemble_cross_grad(srs.space(), pg.space(), n));
}

// S6 (K9): rows = V_k^q (x) S_h^n, cols = broken (x) (qb, bcb); (phi_i, q_j)
template <int D, typename T, typename P, class SPF, class SPB>
kron_term<T, P, SPF> st_assemble_mixed_mass(pg_space<D, T, P, SPF>& pg, int q, int n,
                                            st_broken_space<D, T, P, SPB>& sbs, int qb, time_bc bcb) {
    detail::st_require_same_time_mesh(pg.time_mesh_ref().same_as(sbs.time_mesh_ref()), "st_assemble_mixed_mass");
    return kron_term<T, P, SPF>(T(1), pg.time_dofs(q, time_bc_zero_left).template mixed_mass<P>(sbs.time_dofs(qb, bcb)),
                                bfem::assemble_mixed_mass(pg.space(), n, sbs.space()));
}

namespace detail {

// ===========================================================================
// S11: the residual polynomial w = g(ubar) - dt ubar - b . grad ubar on one
// sub-interval x element (internal 4.1). Local tensors only; never global.
// ===========================================================================
template <int D, typename T>
struct st_residual_scratch {
    st_local<D, T> g_loc, dt_loc, bg_loc, e1, e2, tmp;
};

// w = g(u) - dt u - b . grad u with the three terms lifted to the common time
// degree q_w = max(deg g * q_u, q_u - 1, q_u + q_b) by st_elevate_time_into;
// b == 0 (null) omits the convection term. Space degrees are lifted by the
// row-wise bfem sub_into. The sum is taken as (g - dt) - bg (fixed order).
template <int D, typename T, typename P>
void st_residual_into(st_local<D, T>& w, const bfem::poly1<T>& g, const st_local<D, T>& u,
                      const T& inv_k, const st_local<D, T>* b, st_workspace<D, T>& ws,
                      st_residual_scratch<D, T>& sc, bfem::element_op<D, T, P>& op) {
    assert(&w != &u);
    st_compose_into(sc.g_loc, g, u, ws);
    st_dt_into(sc.dt_loc, u, inv_k, ws);
    int qw = sc.g_loc.qt >= sc.dt_loc.qt ? sc.g_loc.qt : sc.dt_loc.qt;
    if (b) {
        st_b_dot_grad_into(sc.bg_loc, b, u, ws, op);
        if (sc.bg_loc.qt > qw) qw = sc.bg_loc.qt;
    }
    st_elevate_time_into(sc.e1, sc.g_loc, qw, ws);
    st_elevate_time_into(sc.e2, sc.dt_loc, qw, ws);
    if (!b) {
        st_sub_into(w, sc.e1, sc.e2);
        return;
    }
    st_sub_into(sc.tmp, sc.e1, sc.e2);
    st_elevate_time_into(sc.e1, sc.bg_loc, qw, ws);
    st_sub_into(w, sc.tmp, sc.e1);
}

// ---------------------------------------------------------------------------
// validated access to ubar / b of a pg_space (shape checks: pg_space keeps
// its identity check private; the shapes are what the loops rely on)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPF>
void st_check_ubar(pg_space<D, T, P, SPF>& pg, const st_poly<D, T, P>& u, const char* where) {
    const time_space<T>& Vt = pg.time_dofs(u.time_degree(), u.time_bc_kind());
    if (u.num_time_dofs() != Vt.ndof() || u.num_space_dofs() != pg.space().ndof(u.space_degree())
        || u.coeffs().rowsize() != u.num_time_dofs() || u.coeffs().columnsize() != u.num_space_dofs()) {
        std::string msg("btfem::");
        msg += where;
        msg += ": st_poly does not match the pg_space (shape)";
        throw std::invalid_argument(msg);
    }
}

// the residual source: ubar (+ b) of a pg_space, g, and the per-(l, e)
// construction of w. bkind: 0 none, 1 constant b, 2 space-time b.
template <int D, typename T, typename P, class SPF>
struct st_residual_source {
    typedef st_poly<D, T, P> fn_t;

    pg_space<D, T, P, SPF>* pg;
    const fn_t* ubar;
    const bfem::poly1<T>* g;
    int bkind;
    const time_space<T>* Vtu;
    const time_space<T>* Vtb;
    const bfem::dofmap<D>* dmu;
    std::array<const bfem::dofmap<D>*, static_cast<std::size_t>(D)> dmb;
    const std::array<fn_t, static_cast<std::size_t>(D)>* bf;
    int qw;                                     // time degree of w
    int nw;                                     // space degree of w
    st_workspace<D, T> ws;
    st_local<D, T> u_loc;
    std::array<st_local<D, T>, static_cast<std::size_t>(D)> b_loc;
    st_residual_scratch<D, T> sc;

    st_residual_source(pg_space<D, T, P, SPF>& pg_, const bfem::poly1<T>& g_, const fn_t& u_, const char* where)
        : pg(&pg_), ubar(&u_), g(&g_), bkind(0), Vtu(0), Vtb(0), dmu(0), dmb(), bf(0), qw(0), nw(0),
          ws(), u_loc(), b_loc(), sc() {
        st_check_ubar(pg_, u_, where);
        Vtu = &pg_.time_dofs(u_.time_degree(), u_.time_bc_kind());
        dmu = &pg_.space().dofs(u_.space_degree());
        for (int d = 0; d < D; ++d) dmb[static_cast<std::size_t>(d)] = 0;
        finish_degrees(0, 0);
    }
    void set_b(const std::array<T, static_cast<std::size_t>(D)>& bc, const char* where) {
        require_grad(where);
        bkind = 1;
        for (int d = 0; d < D; ++d) {
            b_loc[static_cast<std::size_t>(d)].prepare(0, 0, false);
            bfem::detail::bpoly_access::vec(b_loc[static_cast<std::size_t>(d)].r(0))[0] = bc[static_cast<std::size_t>(d)];
        }
        finish_degrees(0, 0);
    }
    void set_b(const std::array<fn_t, static_cast<std::size_t>(D)>& b, const char* where) {
        require_grad(where);
        for (int d = 0; d < D; ++d) {
            st_check_ubar(*pg, b[static_cast<std::size_t>(d)], where);
            if (b[static_cast<std::size_t>(d)].time_degree() != b[0].time_degree()
                || b[static_cast<std::size_t>(d)].time_bc_kind() != b[0].time_bc_kind()) {
                std::string msg("btfem::");
                msg += where;
                msg += ": components of b must share the time space (q, bc)";
                throw std::invalid_argument(msg);
            }
        }
        bkind = 2;
        bf = &b;
        Vtb = &pg->time_dofs(b[0].time_degree(), b[0].time_bc_kind());
        int nbmax = 0;
        for (int d = 0; d < D; ++d) {
            dmb[static_cast<std::size_t>(d)] = &pg->space().dofs(b[static_cast<std::size_t>(d)].space_degree());
            if (b[static_cast<std::size_t>(d)].space_degree() > nbmax) nbmax = b[static_cast<std::size_t>(d)].space_degree();
        }
        finish_degrees(b[0].time_degree(), nbmax);
    }
    // gather ubar (and b) on (l, e) and build w (op geometry must be set)
    void build(int l, int e, bfem::element_op<D, T, P>& op, st_local<D, T>& w) {
        st_gather(u_loc, ubar->coeffs(), *Vtu, l, *dmu, e);
        if (bkind == 2)
            for (int d = 0; d < D; ++d)
                st_gather(b_loc[static_cast<std::size_t>(d)], (*bf)[static_cast<std::size_t>(d)].coeffs(),
                          *Vtb, l, *dmb[static_cast<std::size_t>(d)], e);
        st_residual_into(w, *g, u_loc, pg->time_mesh_ref().inv_width(),
                         bkind == 0 ? static_cast<const st_local<D, T>*>(0) : b_loc.data(), ws, sc, op);
        assert(w.qt == qw);
    }

private:
    void require_grad(const char* where) const {
        if (ubar->space_degree() < 1) {
            std::string msg("btfem::");
            msg += where;
            msg += ": ubar space degree < 1 (convection term needs grad)";
            throw std::invalid_argument(msg);
        }
    }
    void finish_degrees(int qb, int nb) {
        const int qu = ubar->time_degree(), nu = ubar->space_degree();
        const int K = g->degree();
        qw = K * qu;
        if (qu - 1 > qw) qw = qu - 1;
        nw = K * nu;
        if (nu > nw) nw = nu;
        if (bkind != 0) {
            if (qu + qb > qw) qw = qu + qb;
            if (nu - 1 + nb > nw) nw = nu - 1 + nb;
        }
    }
};

// ===========================================================================
// the common quadratic form (internal 6.1):
//   sum_l sum_e sum_{r, r'} Mt(r, r') * inner.value(r, r'),  Mt = k Mref(qz, qz)
// inner.prepare(l, e) readies the local quantities of (l, e); the loop order
// l -> e -> r -> r' is fixed; no symmetry shortcut.
// ===========================================================================
template <typename T, typename Inner>
T st_quadratic_sq(const time_mesh<T>& tm, int nt, int qz, Inner& inner) {
    const typename time_space<T>::table_type& M = time_tables_typed<T>::mass_ref(qz, qz);
    const T& k = tm.width();
    std::vector<T> Mt(static_cast<std::size_t>((qz + 1) * (qz + 1)));
    for (int r = 0; r <= qz; ++r)
        for (int rp = 0; rp <= qz; ++rp)
            Mt[static_cast<std::size_t>(r * (qz + 1) + rp)] = k * M(r, rp);
    const int m = tm.num_intervals();
    T acc(0);
    for (int l = 1; l <= m; ++l) {
        for (int e = 0; e < nt; ++e) {
            inner.prepare(l, e);
            for (int r = 0; r <= qz; ++r)
                for (int rp = 0; rp <= qz; ++rp)
                    acc += Mt[static_cast<std::size_t>(r * (qz + 1) + rp)] * inner.value(r, rp);
        }
    }
    return acc;
}

// ---------------------------------------------------------------------------
// S7 inner: B_e(r, r') = u_r^T S u_r' - u_r^T X^T s_r' - s_r^T X u_r' + s_r^T M s_r'
//   S = local_stiffness(nbar), X = local_cross_grad(k, nbar), M = local_rt_mass(k)
// (internal 6.2, J-2). sigma and ubar lifted to q_z = max(q_sigma, q_ubar).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPR, class SPF>
struct st_flux_inner {
    st_rt_space<D, T, P, SPR>* srs;
    const st_rt_field<D, T, P>* sig;
    pg_space<D, T, P, SPF>* pg;
    const st_poly<D, T, P>* ubar;
    const time_space<T>* Vts;
    const time_space<T>* Vtu;
    const bfem::dofmap<D>* dmu;
    int k, n, qz;
    bfem::rt_element_op<D, T, P> rop;
    bfem::element_op<D, T, P> op;
    vcp::matrix<T, P> S, X, M;
    st_rt_local<T> sraw, sel;
    st_local<D, T> uraw, uel;
    st_workspace<D, T> ws;

    st_flux_inner(st_rt_space<D, T, P, SPR>& srs_, const st_rt_field<D, T, P>& sig_,
                  pg_space<D, T, P, SPF>& pg_, const st_poly<D, T, P>& u_, const char* where)
        : srs(&srs_), sig(&sig_), pg(&pg_), ubar(&u_), Vts(0), Vtu(0), dmu(0), k(0), n(0), qz(0),
          rop(), op(), S(), X(), M(), sraw(), sel(), uraw(), uel(), ws() {
        srs_.validate(sig_);
        st_check_ubar(pg_, u_, where);
        st_require_same_time_mesh(srs_.time_mesh_ref().same_as(pg_.time_mesh_ref()), where);
        bfem::detail::rt_require_same_mesh(srs_.space().num_elements(), srs_.space().num_vertices(),
                                           pg_.space().num_elements(), pg_.space().ndof(1), where);
        if (u_.space_degree() < 1) {
            std::string msg("btfem::");
            msg += where;
            msg += ": ubar space degree < 1";
            throw std::invalid_argument(msg);
        }
        Vts = &srs_.time_dofs(sig_.time_degree(), sig_.time_bc_kind());
        Vtu = &pg_.time_dofs(u_.time_degree(), u_.time_bc_kind());
        dmu = &pg_.space().dofs(u_.space_degree());
        k = srs_.space().order();
        n = u_.space_degree();
        qz = sig_.time_degree() >= u_.time_degree() ? sig_.time_degree() : u_.time_degree();
    }
    void prepare(int l, int e) {
        rop.set_geometry(srs->space().geometry(e));
        op.set_geometry(pg->space().geometry(e));
        st_gather_rt(sraw, sig->coeffs(), *Vts, l, srs->space().dofs(), e, k);
        st_gather(uraw, ubar->coeffs(), *Vtu, l, *dmu, e);
        st_elevate_time_rt_into(sel, sraw, qz);
        st_elevate_time_into(uel, uraw, qz, ws);
        op.local_stiffness(n, S);
        rop.local_cross_grad(k, n, X);
        rop.local_rt_mass(k, M);
    }
    T value(int r, int rp) const {
        const std::vector<T>& ur = bfem::detail::bpoly_access::vec(uel.r(r));
        const std::vector<T>& urp = bfem::detail::bpoly_access::vec(uel.r(rp));
        const std::vector<T>& sr = sel.r(r).c;
        const std::vector<T>& srp = sel.r(rp).c;
        const int N = S.rowsize();
        const int dim = M.rowsize();
        T uSu(0);
        for (int a = 0; a < N; ++a) {
            T row = S(a, 0) * urp[0];
            for (int b = 1; b < N; ++b) row += S(a, b) * urp[static_cast<std::size_t>(b)];
            uSu += ur[static_cast<std::size_t>(a)] * row;
        }
        T uXs(0);
        for (int a = 0; a < N; ++a) {
            T row = X(0, a) * srp[0];
            for (int j = 1; j < dim; ++j) row += X(j, a) * srp[static_cast<std::size_t>(j)];
            uXs += ur[static_cast<std::size_t>(a)] * row;
        }
        T sXu(0);
        for (int j = 0; j < dim; ++j) {
            T row = X(j, 0) * urp[0];
            for (int a = 1; a < N; ++a) row += X(j, a) * urp[static_cast<std::size_t>(a)];
            sXu += sr[static_cast<std::size_t>(j)] * row;
        }
        T sMs(0);
        for (int i = 0; i < dim; ++i) {
            T row = M(i, 0) * srp[0];
            for (int j = 1; j < dim; ++j) row += M(i, j) * srp[static_cast<std::size_t>(j)];
            sMs += sr[static_cast<std::size_t>(i)] * row;
        }
        return uSu - uXs - sXu + sMs;
    }
};

// ---------------------------------------------------------------------------
// S8 inner: z = div sigma + w, B_e(r, r') = local_inner(z_r, z_r') (internal 6.3)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPR, class SPF>
struct st_divres_inner {
    st_rt_space<D, T, P, SPR>* srs;
    const st_rt_field<D, T, P>* sig;
    st_residual_source<D, T, P, SPF> src;
    const time_space<T>* Vts;
    int k, qz;
    bfem::rt_element_op<D, T, P> rop;
    bfem::element_op<D, T, P> op;
    st_rt_local<T> sraw;
    st_local<D, T> dv, w, dvel, wel, z;
    st_workspace<D, T> ws;

    st_divres_inner(st_rt_space<D, T, P, SPR>& srs_, const st_rt_field<D, T, P>& sig_,
                    const bfem::poly1<T>& g, pg_space<D, T, P, SPF>& pg_, const st_poly<D, T, P>& u_,
                    const char* where)
        : srs(&srs_), sig(&sig_), src(pg_, g, u_, where), Vts(0), k(0), qz(0),
          rop(), op(), sraw(), dv(), w(), dvel(), wel(), z(), ws() {
        srs_.validate(sig_);
        st_require_same_time_mesh(srs_.time_mesh_ref().same_as(pg_.time_mesh_ref()), where);
        bfem::detail::rt_require_same_mesh(srs_.space().num_elements(), srs_.space().num_vertices(),
                                           pg_.space().num_elements(), pg_.space().ndof(1), where);
        Vts = &srs_.time_dofs(sig_.time_degree(), sig_.time_bc_kind());
        k = srs_.space().order();
        update_qz();
    }
    void update_qz() { qz = sig->time_degree() >= src.qw ? sig->time_degree() : src.qw; }
    void prepare(int l, int e) {
        rop.set_geometry(srs->space().geometry(e));
        op.set_geometry(src.pg->space().geometry(e));
        st_gather_rt(sraw, sig->coeffs(), *Vts, l, srs->space().dofs(), e, k);
        dv.prepare(sraw.qt, k, false);
        for (int r = 0; r <= sraw.qt; ++r) rop.div_as_bpoly(sraw.r(r), dv.r(r));
        src.build(l, e, op, w);
        st_elevate_time_into(dvel, dv, qz, ws);
        st_elevate_time_into(wel, w, qz, ws);
        st_add_into(z, dvel, wel);
    }
    T value(int r, int rp) const { return op.local_inner(z.r(r), z.r(rp)); }
};

// ---------------------------------------------------------------------------
// S10 inner: B_e(r, r') = local_inner(w_r, w_r') - D! inv_absdet f_r^T Mi f_r'
//   f_r = local_load(w_r, l_broken), Mi = typed_rt_registry<D,T>::inv_mass(l)
// (internal 6.5; bfem K8 bilinearized, E-2)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPB, class SPF>
struct st_projres_inner {
    st_broken_space<D, T, P, SPB>* sbs;
    st_residual_source<D, T, P, SPF> src;
    int lb, nloc, qz;
    const typename bfem::typed_rt_registry<D, T>::mat_table* Mi;
    T dfact;
    T scale;                                     // D! inv_absdet of the current element
    bfem::element_op<D, T, P> op;
    st_local<D, T> w;
    std::vector<vcp::matrix<T, P> > fv;

    st_projres_inner(st_broken_space<D, T, P, SPB>& sbs_, const bfem::poly1<T>& g,
                     pg_space<D, T, P, SPF>& pg_, const st_poly<D, T, P>& u_, const char* where)
        : sbs(&sbs_), src(pg_, g, u_, where), lb(0), nloc(0), qz(0), Mi(0),
          dfact(static_cast<int>(bfem::detail::factorial_of<D>::value)), scale(), op(), w(), fv() {
        st_require_same_time_mesh(sbs_.time_mesh_ref().same_as(pg_.time_mesh_ref()), where);
        bfem::detail::rt_require_same_mesh(sbs_.space().num_elements(), sbs_.space().num_vertices(),
                                           pg_.space().num_elements(), pg_.space().ndof(1), where);
        lb = sbs_.space().order();
        nloc = sbs_.space().local_size();
        Mi = &bfem::typed_rt_registry<D, T>::inv_mass(lb);
        update_qz();
    }
    void update_qz() {
        qz = src.qw;
        if (fv.size() < static_cast<std::size_t>(qz + 1)) fv.resize(static_cast<std::size_t>(qz + 1));
    }
    void prepare(int l, int e) {
        op.set_geometry(sbs->space().geometry(e));
        src.build(l, e, op, w);
        for (int r = 0; r <= qz; ++r) op.local_load(w.r(r), lb, fv[static_cast<std::size_t>(r)]);
        scale = dfact * sbs->space().geometry(e).inv_absdet();
    }
    T value(int r, int rp) const {
        const vcp::matrix<T, P>& fr = fv[static_cast<std::size_t>(r)];
        const vcp::matrix<T, P>& frp = fv[static_cast<std::size_t>(rp)];
        T q(0);
        for (int i = 0; i < nloc; ++i) {
            T row = Mi->at(i, 0) * frp(0, 0);
            for (int j = 1; j < nloc; ++j) row += Mi->at(i, j) * frp(j, 0);
            q += fr(i, 0) * row;
        }
        return op.local_inner(w.r(r), w.r(rp)) - scale * q;
    }
};

// S9 core: F(gt ndof_b + e nloc + i) += Lt(a, s) fv_s(i), Lt = k Mref(qb, qw)
// (internal 6.4; order l -> e -> s -> a -> i)
template <int D, typename T, typename P, class SPB, class SPF>
vcp::matrix<T, P> st_broken_load_core(st_broken_space<D, T, P, SPB>& sbs, int qb, time_bc bcb,
                                      st_residual_source<D, T, P, SPF>& src, const char* where) {
    pg_space<D, T, P, SPF>& pg = *src.pg;
    st_require_same_time_mesh(sbs.time_mesh_ref().same_as(pg.time_mesh_ref()), where);
    bfem::detail::rt_require_same_mesh(sbs.space().num_elements(), sbs.space().num_vertices(),
                                       pg.space().num_elements(), pg.space().ndof(1), where);
    const time_space<T>& Vtb = sbs.time_dofs(qb, bcb);
    const int lb = sbs.space().order();
    const int nloc = sbs.space().local_size();
    const int nb = sbs.space().ndof();
    const int qw = src.qw;
    const typename time_space<T>::table_type& M = time_tables_typed<T>::mass_ref(qb, qw);
    const T& k = sbs.time_mesh_ref().width();
    vcp::matrix<T, P> F;
    F.zeros(Vtb.ndof() * nb, 1);
    bfem::element_op<D, T, P> op;
    st_local<D, T> w;
    vcp::matrix<T, P> fv;
    const int m = sbs.time_mesh_ref().num_intervals();
    const int nt = sbs.space().num_elements();
    for (int l = 1; l <= m; ++l) {
        for (int e = 0; e < nt; ++e) {
            op.set_geometry(sbs.space().geometry(e));
            src.build(l, e, op, w);
            for (int s = 0; s <= qw; ++s) {
                op.local_load(w.r(s), lb, fv);
                for (int a = 0; a <= qb; ++a) {
                    const int gt = Vtb.global_dof(l, a);
                    if (gt < 0) continue;
                    const T tc = k * M(a, s);
                    for (int i = 0; i < nloc; ++i)
                        F(gt * nb + e * nloc + i, 0) += tc * fv(i, 0);
                }
            }
        }
    }
    return F;
}

} // namespace detail

// ===========================================================================
// S7-S10 (external 6)
// ===========================================================================
// S7 (K5 flux_error_sq): || grad ubar - sigma ||^2_{L2(J;L2)} (nu is the
// user's business, E-1)
template <int D, typename T, typename P, class SPR, class SPF>
T st_flux_error_sq(st_rt_space<D, T, P, SPR>& srs, const typename st_rt_space<D, T, P, SPR>::field_type& sig,
                   pg_space<D, T, P, SPF>& pg, const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_flux_inner<D, T, P, SPR, SPF> inner(srs, sig, pg, ubar, "st_flux_error_sq");
    return detail::st_quadratic_sq(srs.time_mesh_ref(), srs.space().num_elements(), inner.qz, inner);
}

// S8 (K5 div_residual_sq): || div sigma + w ||^2_{L2(J;L2)}, w = g(ubar) - dt ubar [- b . grad ubar]
template <int D, typename T, typename P, class SPR, class SPF>
T st_div_residual_sq(st_rt_space<D, T, P, SPR>& srs, const typename st_rt_space<D, T, P, SPR>::field_type& sig,
                     const bfem::poly1<T>& g, pg_space<D, T, P, SPF>& pg,
                     const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_divres_inner<D, T, P, SPR, SPF> inner(srs, sig, g, pg, ubar, "st_div_residual_sq");
    return detail::st_quadratic_sq(srs.time_mesh_ref(), srs.space().num_elements(), inner.qz, inner);
}
template <int D, typename T, typename P, class SPR, class SPF>
T st_div_residual_sq(st_rt_space<D, T, P, SPR>& srs, const typename st_rt_space<D, T, P, SPR>::field_type& sig,
                     const bfem::poly1<T>& g, const std::array<T, static_cast<std::size_t>(D)>& b_const,
                     pg_space<D, T, P, SPF>& pg, const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_divres_inner<D, T, P, SPR, SPF> inner(srs, sig, g, pg, ubar, "st_div_residual_sq");
    inner.src.set_b(b_const, "st_div_residual_sq");
    inner.update_qz();
    return detail::st_quadratic_sq(srs.time_mesh_ref(), srs.space().num_elements(), inner.qz, inner);
}
template <int D, typename T, typename P, class SPR, class SPF>
T st_div_residual_sq(st_rt_space<D, T, P, SPR>& srs, const typename st_rt_space<D, T, P, SPR>::field_type& sig,
                     const bfem::poly1<T>& g,
                     const std::array<typename pg_space<D, T, P, SPF>::function_type, static_cast<std::size_t>(D)>& b,
                     pg_space<D, T, P, SPF>& pg, const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_divres_inner<D, T, P, SPR, SPF> inner(srs, sig, g, pg, ubar, "st_div_residual_sq");
    inner.src.set_b(b, "st_div_residual_sq");
    inner.update_qz();
    return detail::st_quadratic_sq(srs.time_mesh_ref(), srs.space().num_elements(), inner.qz, inner);
}

// S9 (K7 broken_load): (w, q_i)_{L2(J;L2)}, rows = broken (x) (qb, bcb), time-block major
template <int D, typename T, typename P, class SPB, class SPF>
vcp::matrix<T, P> st_broken_load(st_broken_space<D, T, P, SPB>& sbs, int qb, time_bc bcb,
                                 const bfem::poly1<T>& g, pg_space<D, T, P, SPF>& pg,
                                 const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_residual_source<D, T, P, SPF> src(pg, g, ubar, "st_broken_load");
    return detail::st_broken_load_core(sbs, qb, bcb, src, "st_broken_load");
}
template <int D, typename T, typename P, class SPB, class SPF>
vcp::matrix<T, P> st_broken_load(st_broken_space<D, T, P, SPB>& sbs, int qb, time_bc bcb,
                                 const bfem::poly1<T>& g, const std::array<T, static_cast<std::size_t>(D)>& b_const,
                                 pg_space<D, T, P, SPF>& pg,
                                 const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_residual_source<D, T, P, SPF> src(pg, g, ubar, "st_broken_load");
    src.set_b(b_const, "st_broken_load");
    return detail::st_broken_load_core(sbs, qb, bcb, src, "st_broken_load");
}
template <int D, typename T, typename P, class SPB, class SPF>
vcp::matrix<T, P> st_broken_load(st_broken_space<D, T, P, SPB>& sbs, int qb, time_bc bcb,
                                 const bfem::poly1<T>& g,
                                 const std::array<typename pg_space<D, T, P, SPF>::function_type, static_cast<std::size_t>(D)>& b,
                                 pg_space<D, T, P, SPF>& pg,
                                 const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_residual_source<D, T, P, SPF> src(pg, g, ubar, "st_broken_load");
    src.set_b(b, "st_broken_load");
    return detail::st_broken_load_core(sbs, qb, bcb, src, "st_broken_load");
}

// S10 (K8 projection_error_sq): || (I - Pi_l) w ||^2_{L2(J;L2)}, Pi_l the
// element-wise spatial L2 projection onto broken P_l (time untouched)
template <int D, typename T, typename P, class SPB, class SPF>
T st_projection_error_sq(st_broken_space<D, T, P, SPB>& sbs, const bfem::poly1<T>& g,
                         pg_space<D, T, P, SPF>& pg, const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_projres_inner<D, T, P, SPB, SPF> inner(sbs, g, pg, ubar, "st_projection_error_sq");
    return detail::st_quadratic_sq(sbs.time_mesh_ref(), sbs.space().num_elements(), inner.qz, inner);
}
template <int D, typename T, typename P, class SPB, class SPF>
T st_projection_error_sq(st_broken_space<D, T, P, SPB>& sbs, const bfem::poly1<T>& g,
                         const std::array<T, static_cast<std::size_t>(D)>& b_const,
                         pg_space<D, T, P, SPF>& pg, const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_projres_inner<D, T, P, SPB, SPF> inner(sbs, g, pg, ubar, "st_projection_error_sq");
    inner.src.set_b(b_const, "st_projection_error_sq");
    inner.update_qz();
    return detail::st_quadratic_sq(sbs.time_mesh_ref(), sbs.space().num_elements(), inner.qz, inner);
}
template <int D, typename T, typename P, class SPB, class SPF>
T st_projection_error_sq(st_broken_space<D, T, P, SPB>& sbs, const bfem::poly1<T>& g,
                         const std::array<typename pg_space<D, T, P, SPF>::function_type, static_cast<std::size_t>(D)>& b,
                         pg_space<D, T, P, SPF>& pg, const typename pg_space<D, T, P, SPF>::function_type& ubar) {
    detail::st_projres_inner<D, T, P, SPB, SPF> inner(sbs, g, pg, ubar, "st_projection_error_sq");
    inner.src.set_b(b, "st_projection_error_sq");
    inner.update_qz();
    return detail::st_quadratic_sq(sbs.time_mesh_ref(), sbs.space().num_elements(), inner.qz, inner);
}

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_ST_RT_ASSEMBLE_HPP
