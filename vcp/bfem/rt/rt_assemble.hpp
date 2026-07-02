// vcp/bfem/rt/rt_assemble.hpp
// RT Layer 3: cross-family global assemblers and global scalars
// (K3, K4, K5, K7, K8, K9 -- free functions, Q4/Q8).
//
// Conforms to: RT-L3 external design v0.2 (section 5) and
//              RT-L3 internal design v0.2 (section 4).
//
// Contracts:
//  - same-mesh: both spaces must be built from the same mesh; checked by
//    element/vertex counts (identity itself is a documented contract, Q2);
//    violation raises std::invalid_argument.
//  - the free functions borrow dofmaps and geometry from the spaces and run
//    their OWN element op and buffers ("do not mutate the space", external
//    design section 5). Per-call setup allocation is allowed; the element
//    loop is allocation free from the second element on (A-2 / RK-5).
//  - element-order accumulation, deterministic combination, SP rules as in
//    the single-family assemblers.
//  - scatter: the RT side runs SIGNED (Y2 general kernels / the two-dofmap
//    mixed form below); fe and broken sides are identity.
//  - the returned sparse matrix uses the ROW side space's SP policy.

#ifndef VCP_BFEM_RT_RT_ASSEMBLE_HPP
#define VCP_BFEM_RT_RT_ASSEMBLE_HPP

#include <vector>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_element_op.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// two-dofmap (row x column) general scatter -- the Y2 mixed form of the
// promotion (RT-L3 internal design 2.3). Signs come from BOTH sides'
// dof_sign; identity-family maps return the constant +1.
template <typename T, typename Buf, typename Loc, typename DMR, typename DMC>
void scatter_matrix_general2(const DMR& dmr, const DMC& dmc, int e,
                             const Loc& loc, int nr, int nc, Buf& buf) {
    for (int a = 0; a < nr; ++a) {
        int ga = dmr.global_dof(e, a);
        T sa = T(dmr.dof_sign(e, a));
        for (int b = 0; b < nc; ++b)
            buf.push(ga, dmc.global_dof(e, b),
                     sa * T(dmc.dof_sign(e, b)) * loc(a, b));
    }
}

inline void rt_require_same_mesh(int nt_a, int nv_a, int nt_b, int nv_b,
                                 const char* where) {
    if (nt_a != nt_b || nv_a != nv_b) {
        std::string msg("bfem::");
        msg += where;
        msg += ": spaces are not built from the same mesh";
        throw std::invalid_argument(msg);
    }
}

// gather an fe_function's local block as a degree-m bpoly (identity l2g via
// the public dofmap API -- fe_space privates are never touched)
template <int D, typename T, typename P>
void rt_gather_fe(const dofmap<D>& dm, int e, const vcp::matrix<T, P>& g,
                  bpoly<D, T>& dst) {
    bpoly_access::prepare(dst, dm.degree(), false);
    gather_identity(dm, e, g, bpoly_access::vec(dst).data(), dm.local_size());
}

} // namespace detail

// ---------------------------------------------------------------------------
// K3 (H2): rows = broken P_l, cols = RT^k; entries (div sigma_j, q_i)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPB, class SPR>
vcp::spmatrix<T, SPB> assemble_div_mass(broken_space<D, T, P, SPB>& bs,
                                        rt_space<D, T, P, SPR>& rs) {
    detail::rt_require_same_mesh(bs.num_elements(), bs.num_vertices(),
                                 rs.num_elements(), rs.num_vertices(),
                                 "assemble_div_mass");
    const int k = rs.order();
    const int l = bs.order();
    rt_element_op<D, T, P> op;
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    const int nr = bs.local_size();
    const int nc = rs.dofs().local_size();
    buf.reserve(static_cast<std::size_t>(bs.num_elements())
                * static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc));
    for (int e = 0; e < bs.num_elements(); ++e) {
        op.set_geometry(rs.geometry(e));
        op.local_div_mass(k, l, loc);
        detail::scatter_matrix_general2<T>(bs.dofs(), rs.dofs(), e, loc, nr, nc, buf);
    }
    buf.combine();
    return detail::spm_adapter<T, SPB>::build(bs.ndof(), rs.ndof(), buf);
}

// ---------------------------------------------------------------------------
// K4 (H3): rows = RT^k, cols = V_h^m (m >= 1); entries (sigma_j, grad psi_a)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPR, class SPF>
vcp::spmatrix<T, SPR> assemble_cross_grad(rt_space<D, T, P, SPR>& rs,
                                          fe_space<D, T, P, SPF>& fs, int m) {
    if (m < 1)
        throw std::invalid_argument("bfem::assemble_cross_grad: m < 1");
    detail::rt_require_same_mesh(rs.num_elements(), rs.num_vertices(),
                                 fs.num_elements(), fs.ndof(1),
                                 "assemble_cross_grad");
    const int k = rs.order();
    const dofmap<D>& fdm = fs.dofs(m);
    rt_element_op<D, T, P> op;
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    const int nr = rs.dofs().local_size();
    const int nc = fdm.local_size();
    buf.reserve(static_cast<std::size_t>(rs.num_elements())
                * static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc));
    for (int e = 0; e < rs.num_elements(); ++e) {
        op.set_geometry(rs.geometry(e));
        op.local_cross_grad(k, m, loc);
        detail::scatter_matrix_general2<T>(rs.dofs(), fdm, e, loc, nr, nc, buf);
    }
    buf.combine();
    return detail::spm_adapter<T, SPR>::build(rs.ndof(), fs.ndof(m), buf);
}

// ---------------------------------------------------------------------------
// K5 (H4): global scalars
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPR, class SPF>
T flux_error_sq(rt_space<D, T, P, SPR>& rs, const rt_field<D, T, P>& sig,
                fe_space<D, T, P, SPF>& fs, const fe_function<D, T, P>& u) {
    detail::rt_require_same_mesh(rs.num_elements(), rs.num_vertices(),
                                 fs.num_elements(), fs.ndof(1), "flux_error_sq");
    if (sig.order() != rs.order()
        || sig.coeffs().rowsize() != rs.ndof() || sig.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::flux_error_sq: rt_field mismatch");
    const dofmap<D>& fdm = fs.dofs(u.degree());
    if (u.coeffs().rowsize() != fdm.ndof() || u.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::flux_error_sq: fe_function mismatch");
    rt_element_op<D, T, P> op;
    rt_local_coeffs<T> sloc;
    sloc.k = rs.order();
    sloc.c.assign(static_cast<std::size_t>(rs.dofs().local_size()), T(0));
    bpoly<D, T> uloc;
    T acc(0);
    for (int e = 0; e < rs.num_elements(); ++e) {        // element order (X9)
        op.set_geometry(rs.geometry(e));
        detail::gather_general(rs.dofs(), e, sig.coeffs(), sloc.c.data(),
                               rs.dofs().local_size());
        detail::rt_gather_fe(fdm, e, u.coeffs(), uloc);
        acc += op.local_flux_error_sq(sloc, uloc);
    }
    return acc;
}

// || div sigma + f(u_h) ||^2_{L2(Omega)} (w = f(u_h) composed per element)
template <int D, typename T, typename P, class SPR, class SPF>
T div_residual_sq(rt_space<D, T, P, SPR>& rs, const rt_field<D, T, P>& sig,
                  const poly1<T>& f, fe_space<D, T, P, SPF>& fs,
                  const fe_function<D, T, P>& u) {
    detail::rt_require_same_mesh(rs.num_elements(), rs.num_vertices(),
                                 fs.num_elements(), fs.ndof(1), "div_residual_sq");
    if (sig.order() != rs.order()
        || sig.coeffs().rowsize() != rs.ndof() || sig.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::div_residual_sq: rt_field mismatch");
    const dofmap<D>& fdm = fs.dofs(u.degree());
    if (u.coeffs().rowsize() != fdm.ndof() || u.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::div_residual_sq: fe_function mismatch");
    rt_element_op<D, T, P> op;
    rt_local_coeffs<T> sloc;
    sloc.k = rs.order();
    sloc.c.assign(static_cast<std::size_t>(rs.dofs().local_size()), T(0));
    bpoly<D, T> uloc, w;
    compose_workspace<D, T> cws;
    T acc(0);
    for (int e = 0; e < rs.num_elements(); ++e) {
        op.set_geometry(rs.geometry(e));
        detail::gather_general(rs.dofs(), e, sig.coeffs(), sloc.c.data(),
                               rs.dofs().local_size());
        detail::rt_gather_fe(fdm, e, u.coeffs(), uloc);
        compose_into(w, f, uloc, cws);
        acc += op.local_div_residual_sq(sloc, w);
    }
    return acc;
}

// ---------------------------------------------------------------------------
// K7 (v0.2): KKT right hand side f_v -- g_i = (f(u_h), q_i), q in broken P_l
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPB, class SPF>
vcp::matrix<T, P> broken_load(broken_space<D, T, P, SPB>& bs, const poly1<T>& f,
                              fe_space<D, T, P, SPF>& fs,
                              const fe_function<D, T, P>& uh) {
    detail::rt_require_same_mesh(bs.num_elements(), bs.num_vertices(),
                                 fs.num_elements(), fs.ndof(1), "broken_load");
    const dofmap<D>& fdm = fs.dofs(uh.degree());
    if (uh.coeffs().rowsize() != fdm.ndof() || uh.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::broken_load: fe_function mismatch");
    const int l = bs.order();
    const int nloc = bs.local_size();
    element_op<D, T, P> op;
    bpoly<D, T> uloc, w;
    compose_workspace<D, T> cws;
    vcp::matrix<T, P> loc;
    vcp::matrix<T, P> F;
    F.zeros(bs.ndof(), 1);
    for (int e = 0; e < bs.num_elements(); ++e) {
        op.set_geometry(bs.geometry(e));
        detail::rt_gather_fe(fdm, e, uh.coeffs(), uloc);
        compose_into(w, f, uloc, cws);
        op.local_load(w, l, loc);
        for (int i = 0; i < nloc; ++i)                   // block-contiguous write
            F(e * nloc + i, 0) = loc(i, 0);
    }
    return F;
}

// ---------------------------------------------------------------------------
// K8 (v0.2): || f(u_h) - Pi_{M_l} f(u_h) ||^2_{L2(Omega)}.
// Element contribution (external design section 5):
//   (w, w)_T - D! inv_absdet (f_{v,T}^T Mhat^{-1} f_{v,T}),   w = f(u_h),
// with Mhat^{-1} the EXACT inverse reference mass (typed T-R7). ZERO added
// divisions (RK-6); interval T yields a rigorous enclosure.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPB, class SPF>
T projection_error_sq(broken_space<D, T, P, SPB>& bs, const poly1<T>& f,
                      fe_space<D, T, P, SPF>& fs, const fe_function<D, T, P>& uh) {
    detail::rt_require_same_mesh(bs.num_elements(), bs.num_vertices(),
                                 fs.num_elements(), fs.ndof(1),
                                 "projection_error_sq");
    const dofmap<D>& fdm = fs.dofs(uh.degree());
    if (uh.coeffs().rowsize() != fdm.ndof() || uh.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::projection_error_sq: fe_function mismatch");
    const int l = bs.order();
    const int nloc = bs.local_size();
    const typename typed_rt_registry<D, T>::mat_table& Mi =
        typed_rt_registry<D, T>::inv_mass(l);
    const T dfact = T(static_cast<int>(detail::factorial_of<D>::value));
    element_op<D, T, P> op;
    bpoly<D, T> uloc, w;
    compose_workspace<D, T> cws;
    vcp::matrix<T, P> fv;
    T acc(0);
    for (int e = 0; e < bs.num_elements(); ++e) {
        op.set_geometry(bs.geometry(e));
        detail::rt_gather_fe(fdm, e, uh.coeffs(), uloc);
        compose_into(w, f, uloc, cws);
        // (w, w)_T
        acc += op.local_inner(w, w);
        // f_{v,T} = |T| M^{(l, deg w)} coeffs(w)  (the local load vector)
        op.local_load(w, l, fv);
        // quadratic form with the exact inverse reference mass
        T q(0);
        for (int i = 0; i < nloc; ++i) {
            T row = Mi.at(i, 0) * fv(0, 0);
            for (int j = 1; j < nloc; ++j)
                row += Mi.at(i, j) * fv(j, 0);
            q += fv(i, 0) * row;
        }
        acc -= (dfact * bs.geometry(e).inv_absdet()) * q;
    }
    return acc;
}

// ---------------------------------------------------------------------------
// K9 (v0.2): B_{ij} = (phi_i, q_j); rows = fe (degree m), cols = broken P_l
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPF, class SPB>
vcp::spmatrix<T, SPF> assemble_mixed_mass(fe_space<D, T, P, SPF>& fs, int m,
                                          broken_space<D, T, P, SPB>& bs) {
    if (m < 1)
        throw std::invalid_argument("bfem::assemble_mixed_mass: m < 1");
    detail::rt_require_same_mesh(fs.num_elements(), fs.ndof(1),
                                 bs.num_elements(), bs.num_vertices(),
                                 "assemble_mixed_mass");
    const dofmap<D>& fdm = fs.dofs(m);
    const int l = bs.order();
    element_op<D, T, P> op;
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    const int nr = fdm.local_size();
    const int nc = bs.local_size();
    buf.reserve(static_cast<std::size_t>(bs.num_elements())
                * static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc));
    for (int e = 0; e < bs.num_elements(); ++e) {
        op.set_geometry(bs.geometry(e));
        op.local_mass(m, l, loc);
        // identity x identity: plain pushes (no sign multiplication)
        for (int i = 0; i < nr; ++i) {
            int gi = fdm.global_dof(e, i);
            for (int j = 0; j < nc; ++j)
                buf.push(gi, bs.dofs().global_dof(e, j), loc(i, j));
        }
    }
    buf.combine();
    return detail::spm_adapter<T, SPF>::build(fs.ndof(m), bs.ndof(), buf);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_ASSEMBLE_HPP
