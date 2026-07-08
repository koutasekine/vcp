// vcp/bfem/sv/sv_assemble.hpp
// Phase 5d (Scott-Vogelius parts): the NS matrix / vector / scalar family
// (V5, external design section 7) and the ONE new kernel of the phase,
// local_grad_mixed_mass (v0.2, A-1 / SV-6).
//
// Conforms to: SV external design v0.2 (section 7) and
//              SV internal design v0.2 (section 6, 6.1).
//
// K-style free functions (Q4/Q8 discipline of rt_assemble.hpp): the
// assemblers borrow dofmaps from the wrapped scalar fe_space and geometry
// from vfe_space (broken_space on the pressure side) and run their OWN
// element op and buffers. Per-call setup allocation is allowed; the element
// loop is allocation free from the second element on. Element-order
// accumulation, deterministic combination and the SP rules are inherited
// from the frozen assemblers. The returned sparse matrix uses the ROW side
// space's SP policy.
//
// Division discipline (normative): NO division operator appears in this file
// (G-SV-1). Geometry factors (measure, grad_lambda) are the precomputed
// members of element_geometry; every kernel below is "table reads +
// multiplications + additions", so interval T yields rigorous enclosures by
// the standard three-step argument (L0 table enclosure + ring ops).
//
// Wiring (SV-6: everything except the one new kernel is a re-wiring of
// frozen kernels):
//   vector_stiffness        local_stiffness(m) computed ONCE per element,
//                           placed D times (offset views of the Y2 identity
//                           scatter -- RSV-4 counts this)
//   div_velocity            local_grad_mixed_mass(d, m, l) per component,
//                           broken rows x fe columns (K9-style plain push)
//   advection C(w)          gather w -> local_convection(w, m), block diag
//   advection_derivative    grad_component(w_c, d) -> local_weighted_mass,
//   D(w)                    (c, d) dense blocks in one element loop
//   advection_vector        b_dot_grad(a, b_c) -> local_load, per component
//   advection_scalar        b_dot_grad(a, b_c) + component L1 inner
//   div_field / div_norm_sq grad_component sum (exact), broken coefficients

#ifndef VCP_BFEM_SV_SV_ASSEMBLE_HPP
#define VCP_BFEM_SV_SV_ASSEMBLE_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/element_op.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/typed_tables.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/sv/vfe_space.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// local_grad_mixed_mass (internal design 6.1 -- the only new kernel):
//   out(i, j) = (d/dx_d phi_j^{(m)}, q_i^{(l)})_T
//             = |T| sum_{i'} gl(i', d) * m * M^{(l, m-1)}[i, shift(j, i')]
// with shift = derivative_map<D>(m) (vanishing terms contribute nothing) and
// gl(i', d) = component d of grad lambda_{i'} (precomputed, multiplications
// only). Counting profile (S-SV-2, exact): T multiplications =
//   2 (D+1)                          (the weights w = measure * gl * m)
// + N_l * #{(j, i') : alpha_j has alpha_{i'} > 0}
// which is bounded by (D+1) N_l N_m + O(N_m); divisions = 0.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P>
void sv_local_grad_mixed_mass(const element_geometry<D, T>& g, int d,
                              int m, int l, vcp::matrix<T, P>& out) {
    if (m < 1)
        throw std::invalid_argument("bfem::sv_local_grad_mixed_mass: m < 1");
    if (l < 0)
        throw std::invalid_argument("bfem::sv_local_grad_mixed_mass: l < 0");
    if (d < 0 || d >= D)
        throw std::invalid_argument("bfem::sv_local_grad_mixed_mass: bad d");
    const typed_mass_table<D, T>& M = typed_registry<D, T>::mass(l, m - 1);
    const derivative_map<D>& dm = detail::deriv_cache<D>::get(m);
    const int Nl = M.rows();
    const int Nm = dm.source_size();
    out.zeros(Nl, Nm);
    const T tm = T(m);
    for (int ip = 0; ip <= D; ++ip) {
        const T w = g.measure() * g.grad_lambda(ip, d) * tm;
        for (int j = 0; j < Nm; ++j) {
            int r = dm.target(j, ip);
            if (r < 0) continue;                     // alpha_{i'} == 0: vanishes
            for (int i = 0; i < Nl; ++i)
                out(i, j) += w * M.at(i, r);
        }
    }
}

template <int D>
inline void sv_require_same_mesh(int nt_a, int nv_a, int nt_b, int nv_b,
                                 const char* where) {
    if (nt_a != nt_b || nv_a != nv_b) {
        std::string msg("bfem::");
        msg += where;
        msg += ": spaces are not built from the same mesh";
        throw std::invalid_argument(msg);
    }
}

// gather component d of a vfe_function on element e as a degree-n bpoly
// (identity Y2 gather through an offset view of the scalar dofmap)
template <int D, typename T, typename P>
void sv_gather_component(const dofmap<D>& dm, int N, int e,
                         const vfe_function<D, T, P>& u, int d,
                         bpoly<D, T>& dst) {
    bpoly_access::prepare(dst, dm.degree(), false);
    sv_offset_dofmap<dofmap<D> > odm(dm, d * N);
    gather(odm, e, u.coeffs(), bpoly_access::vec(dst).data(), dm.local_size(),
           typename dofmap<D>::family_tag());
}

} // namespace detail

// ---------------------------------------------------------------------------
// vector stiffness (grad u : grad v): blkdiag(S), computed once per element
// and placed D times (RSV-4)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
typename vfe_space<D, T, P, SP>::spmatrix_t
assemble_vector_stiffness(vfe_space<D, T, P, SP>& vs, int m) {
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dm = fs.dofs(m);
    const int N = dm.ndof();
    const int nloc = dm.local_size();
    element_op<D, T, P> op;
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    buf.reserve(static_cast<std::size_t>(vs.num_elements())
                * static_cast<std::size_t>(D)
                * static_cast<std::size_t>(nloc) * static_cast<std::size_t>(nloc));
#if VCP_BFEM_USE_OPENMP
    const int nt = vs.num_elements();
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
                  * static_cast<std::size_t>(D)
                  * static_cast<std::size_t>(nloc)
                  * static_cast<std::size_t>(nloc));
        for (int e = e0; e < e1; ++e) {
            op_l.set_geometry(vs.geometry(e));
            op_l.local_stiffness(m, loc_l);
            for (int d = 0; d < D; ++d) {
                detail::sv_offset_dofmap<dofmap<D> > odm(dm, d * N);
                detail::scatter_matrix<T>(odm, e, loc_l, nloc, b,
                                          typename dofmap<D>::family_tag());
            }
        }
    }
    buf.append_all(tbuf, nrun);
#else
    for (int e = 0; e < vs.num_elements(); ++e) {            // element order (X9)
        op.set_geometry(vs.geometry(e));
        op.local_stiffness(m, loc);                          // computed ONCE
        for (int d = 0; d < D; ++d) {                        // placed D times
            detail::sv_offset_dofmap<dofmap<D> > odm(dm, d * N);
            detail::scatter_matrix<T>(odm, e, loc, nloc, buf,
                                      typename dofmap<D>::family_tag());
        }
    }
#endif
    buf.combine();
    return detail::spm_adapter<T, SP>::build(D * N, D * N, buf);
}

// ---------------------------------------------------------------------------
// div matrix (div u, q): rows = broken P_l, cols = (V_h^m)^D; component
// block d = local_grad_mixed_mass(d, m, l)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPB, class SPV>
vcp::spmatrix<T, SPB>
assemble_div_velocity(broken_space<D, T, P, SPB>& bs,
                      vfe_space<D, T, P, SPV>& vs, int m) {
    detail::sv_require_same_mesh<D>(bs.num_elements(), bs.num_vertices(),
                                    vs.num_elements(), vs.num_vertices(),
                                    "assemble_div_velocity");
    fe_space<D, T, P, SPV>& fs = vs.scalar();
    const dofmap<D>& dm = fs.dofs(m);
    const int N = dm.ndof();
    const int nr = bs.local_size();
    const int nc = dm.local_size();
    const int l = bs.order();
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    buf.reserve(static_cast<std::size_t>(vs.num_elements())
                * static_cast<std::size_t>(D)
                * static_cast<std::size_t>(nr) * static_cast<std::size_t>(nc));
#if VCP_BFEM_USE_OPENMP
    const int nt = vs.num_elements();
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
        vcp::matrix<T, P> loc_l;
        detail::coo_buffer<T>& b = tbuf[static_cast<std::size_t>(tid)];
        b.reserve(static_cast<std::size_t>(e1 - e0)
                  * static_cast<std::size_t>(D)
                  * static_cast<std::size_t>(nr)
                  * static_cast<std::size_t>(nc));
        for (int e = e0; e < e1; ++e) {
            for (int d = 0; d < D; ++d) {
                detail::sv_local_grad_mixed_mass<D, T, P>(vs.geometry(e), d, m, l, loc_l);
                detail::sv_scatter_rect_identity<T>(bs.dofs(), 0, dm, d * N,
                                                    e, loc_l, nr, nc, b);
            }
        }
    }
    buf.append_all(tbuf, nrun);
#else
    for (int e = 0; e < vs.num_elements(); ++e) {
        for (int d = 0; d < D; ++d) {
            detail::sv_local_grad_mixed_mass<D, T, P>(vs.geometry(e), d, m, l, loc);
            detail::sv_scatter_rect_identity<T>(bs.dofs(), 0, dm, d * N,
                                                e, loc, nr, nc, buf);
        }
    }
#endif
    buf.combine();
    return detail::spm_adapter<T, SPB>::build(bs.ndof(), D * N, buf);
}

// ---------------------------------------------------------------------------
// advection C(w): (w . grad u, v), block diagonal of the frozen convection
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
typename vfe_space<D, T, P, SP>::spmatrix_t
assemble_advection(vfe_space<D, T, P, SP>& vs,
                   const vfe_function<D, T, P>& w, int m) {
    vs.validate(w);
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dm = fs.dofs(m);
    const dofmap<D>& dmw = fs.dofs(w.degree());
    const int N = dm.ndof();
    const int Nw = dmw.ndof();
    const int nloc = dm.local_size();
    element_op<D, T, P> op;
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    std::array<bpoly<D, T>, D> wloc;
    buf.reserve(static_cast<std::size_t>(vs.num_elements())
                * static_cast<std::size_t>(D)
                * static_cast<std::size_t>(nloc) * static_cast<std::size_t>(nloc));
#if VCP_BFEM_USE_OPENMP
    const int nt = vs.num_elements();
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
        std::array<bpoly<D, T>, D> wloc_l;
        detail::coo_buffer<T>& b = tbuf[static_cast<std::size_t>(tid)];
        b.reserve(static_cast<std::size_t>(e1 - e0)
                  * static_cast<std::size_t>(D)
                  * static_cast<std::size_t>(nloc)
                  * static_cast<std::size_t>(nloc));
        for (int e = e0; e < e1; ++e) {
            op_l.set_geometry(vs.geometry(e));
            for (int d = 0; d < D; ++d)
                detail::sv_gather_component(dmw, Nw, e, w, d,
                                            wloc_l[static_cast<std::size_t>(d)]);
            op_l.local_convection(wloc_l, m, loc_l);
            for (int d = 0; d < D; ++d) {
                detail::sv_offset_dofmap<dofmap<D> > odm(dm, d * N);
                detail::scatter_matrix<T>(odm, e, loc_l, nloc, b,
                                          typename dofmap<D>::family_tag());
            }
        }
    }
    buf.append_all(tbuf, nrun);
#else
    for (int e = 0; e < vs.num_elements(); ++e) {
        op.set_geometry(vs.geometry(e));
        for (int d = 0; d < D; ++d)
            detail::sv_gather_component(dmw, Nw, e, w, d,
                                        wloc[static_cast<std::size_t>(d)]);
        op.local_convection(wloc, m, loc);                   // computed ONCE
        for (int d = 0; d < D; ++d) {                        // placed D times
            detail::sv_offset_dofmap<dofmap<D> > odm(dm, d * N);
            detail::scatter_matrix<T>(odm, e, loc, nloc, buf,
                                      typename dofmap<D>::family_tag());
        }
    }
#endif
    buf.combine();
    return detail::spm_adapter<T, SP>::build(D * N, D * N, buf);
}

// ---------------------------------------------------------------------------
// advection derivative D(w): ((u . grad) w, v); block (c, d) is the weighted
// mass with weight d w_c / d x_d = grad_component(w_c, d)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
typename vfe_space<D, T, P, SP>::spmatrix_t
assemble_advection_derivative(vfe_space<D, T, P, SP>& vs,
                              const vfe_function<D, T, P>& w, int m) {
    vs.validate(w);
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dm = fs.dofs(m);
    const dofmap<D>& dmw = fs.dofs(w.degree());
    const int N = dm.ndof();
    const int Nw = dmw.ndof();
    const int nloc = dm.local_size();
    element_op<D, T, P> op;
    detail::coo_buffer<T> buf;
    vcp::matrix<T, P> loc;
    bpoly<D, T> wc, gd;
    buf.reserve(static_cast<std::size_t>(vs.num_elements())
                * static_cast<std::size_t>(D) * static_cast<std::size_t>(D)
                * static_cast<std::size_t>(nloc) * static_cast<std::size_t>(nloc));
#if VCP_BFEM_USE_OPENMP
    const int nt = vs.num_elements();
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
        bpoly<D, T> wc_l, gd_l;
        detail::coo_buffer<T>& b = tbuf[static_cast<std::size_t>(tid)];
        b.reserve(static_cast<std::size_t>(e1 - e0)
                  * static_cast<std::size_t>(D)
                  * static_cast<std::size_t>(D)
                  * static_cast<std::size_t>(nloc)
                  * static_cast<std::size_t>(nloc));
        for (int e = e0; e < e1; ++e) {
            op_l.set_geometry(vs.geometry(e));
            for (int c = 0; c < D; ++c) {
                detail::sv_gather_component(dmw, Nw, e, w, c, wc_l);
                for (int d = 0; d < D; ++d) {
                    op_l.grad_component(wc_l, d, gd_l);
                    op_l.local_weighted_mass(gd_l, m, loc_l);
                    detail::sv_scatter_rect_identity<T>(dm, c * N, dm, d * N,
                                                        e, loc_l, nloc, nloc, b);
                }
            }
        }
    }
    buf.append_all(tbuf, nrun);
#else
    for (int e = 0; e < vs.num_elements(); ++e) {
        op.set_geometry(vs.geometry(e));
        for (int c = 0; c < D; ++c) {
            detail::sv_gather_component(dmw, Nw, e, w, c, wc);
            for (int d = 0; d < D; ++d) {
                op.grad_component(wc, d, gd);
                op.local_weighted_mass(gd, m, loc);
                detail::sv_scatter_rect_identity<T>(dm, c * N, dm, d * N,
                                                    e, loc, nloc, nloc, buf);
            }
        }
    }
#endif
    buf.combine();
    return detail::spm_adapter<T, SP>::build(D * N, D * N, buf);
}

// ---------------------------------------------------------------------------
// advection vector ((a . grad) b, phi_i): the load form (no matrix is built)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
vcp::matrix<T, P> advection_vector(vfe_space<D, T, P, SP>& vs,
                                   const vfe_function<D, T, P>& a,
                                   const vfe_function<D, T, P>& b, int m) {
    vs.validate(a);
    vs.validate(b);
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dm = fs.dofs(m);
    const dofmap<D>& dma = fs.dofs(a.degree());
    const dofmap<D>& dmb = fs.dofs(b.degree());
    const int N = dm.ndof();
    const int Na = dma.ndof();
    const int Nb = dmb.ndof();
    element_op<D, T, P> op;
    vcp::matrix<T, P> loc;
    std::array<bpoly<D, T>, D> aloc;
    bpoly<D, T> bc, wb;
    vcp::matrix<T, P> F;
    F.zeros(D * N, 1);
    for (int e = 0; e < vs.num_elements(); ++e) {            // element order (X9)
        op.set_geometry(vs.geometry(e));
        for (int d = 0; d < D; ++d)
            detail::sv_gather_component(dma, Na, e, a, d,
                                        aloc[static_cast<std::size_t>(d)]);
        for (int c = 0; c < D; ++c) {
            detail::sv_gather_component(dmb, Nb, e, b, c, bc);
            op.b_dot_grad(aloc, bc, wb);
            op.local_load(wb, m, loc);
            detail::sv_offset_dofmap<dofmap<D> > odm(dm, c * N);
            detail::scatter_vector(odm, e, loc.data(), dm.local_size(), F,
                                   typename dofmap<D>::family_tag());
        }
    }
    return F;
}

// ---------------------------------------------------------------------------
// advection scalar ((a . grad) b, g): g a vfe_function, or g = grad u_h of a
// scalar fe_function (the T-SV-7 gradient-field path)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
T advection_scalar(vfe_space<D, T, P, SP>& vs,
                   const vfe_function<D, T, P>& a,
                   const vfe_function<D, T, P>& b,
                   const vfe_function<D, T, P>& g) {
    vs.validate(a);
    vs.validate(b);
    vs.validate(g);
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dma = fs.dofs(a.degree());
    const dofmap<D>& dmb = fs.dofs(b.degree());
    const dofmap<D>& dmg = fs.dofs(g.degree());
    const int Na = dma.ndof();
    const int Nb = dmb.ndof();
    const int Ng = dmg.ndof();
    element_op<D, T, P> op;
    std::array<bpoly<D, T>, D> aloc;
    bpoly<D, T> bc, wb, gc;
    T acc(0);
    for (int e = 0; e < vs.num_elements(); ++e) {            // element order (X9)
        op.set_geometry(vs.geometry(e));
        for (int d = 0; d < D; ++d)
            detail::sv_gather_component(dma, Na, e, a, d,
                                        aloc[static_cast<std::size_t>(d)]);
        for (int c = 0; c < D; ++c) {
            detail::sv_gather_component(dmb, Nb, e, b, c, bc);
            op.b_dot_grad(aloc, bc, wb);
            detail::sv_gather_component(dmg, Ng, e, g, c, gc);
            acc += op.local_inner(wb, gc);
        }
    }
    return acc;
}

template <int D, typename T, typename P, class SP>
T advection_scalar(vfe_space<D, T, P, SP>& vs,
                   const vfe_function<D, T, P>& a,
                   const vfe_function<D, T, P>& b,
                   const fe_function<D, T, P>& u) {           // g = grad u_h
    vs.validate(a);
    vs.validate(b);
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dma = fs.dofs(a.degree());
    const dofmap<D>& dmb = fs.dofs(b.degree());
    const dofmap<D>& dmu = fs.dofs(u.degree());
    if (u.coeffs().rowsize() != dmu.ndof() || u.coeffs().columnsize() != 1)
        throw std::invalid_argument("bfem::advection_scalar: fe_function mismatch");
    const int Na = dma.ndof();
    const int Nb = dmb.ndof();
    element_op<D, T, P> op;
    std::array<bpoly<D, T>, D> aloc;
    bpoly<D, T> bc, wb, ul, gc;
    T acc(0);
    for (int e = 0; e < vs.num_elements(); ++e) {
        op.set_geometry(vs.geometry(e));
        for (int d = 0; d < D; ++d)
            detail::sv_gather_component(dma, Na, e, a, d,
                                        aloc[static_cast<std::size_t>(d)]);
        detail::bpoly_access::prepare(ul, dmu.degree(), false);
        detail::gather(dmu, e, u.coeffs(),
                       detail::bpoly_access::vec(ul).data(), dmu.local_size(),
                       typename dofmap<D>::family_tag());
        for (int c = 0; c < D; ++c) {
            detail::sv_gather_component(dmb, Nb, e, b, c, bc);
            op.b_dot_grad(aloc, bc, wb);
            op.grad_component(ul, c, gc);
            acc += op.local_inner(wb, gc);
        }
    }
    return acc;
}

// ---------------------------------------------------------------------------
// div_field: div u_h in broken P_{m-1}, exact (grad_component sum); the
// coefficients are written block by block (assignment, deterministic)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SPB, class SPV>
broken_field<D, T, P> div_field(broken_space<D, T, P, SPB>& bs,
                                vfe_space<D, T, P, SPV>& vs,
                                const vfe_function<D, T, P>& u) {
    vs.validate(u);
    detail::sv_require_same_mesh<D>(bs.num_elements(), bs.num_vertices(),
                                    vs.num_elements(), vs.num_vertices(),
                                    "div_field");
    if (bs.order() != u.degree() - 1)
        throw std::invalid_argument("bfem::div_field: broken order != degree - 1");
    fe_space<D, T, P, SPV>& fs = vs.scalar();
    const dofmap<D>& dmu = fs.dofs(u.degree());
    const int Nu = dmu.ndof();
    const int nloc = bs.local_size();
    element_op<D, T, P> op;
    bpoly<D, T> uc, gd, acc;
    broken_field<D, T, P> out = bs.zero_field();
    for (int e = 0; e < vs.num_elements(); ++e) {
        op.set_geometry(vs.geometry(e));
        detail::sv_gather_component(dmu, Nu, e, u, 0, uc);
        op.grad_component(uc, 0, acc);
        for (int d = 1; d < D; ++d) {
            detail::sv_gather_component(dmu, Nu, e, u, d, uc);
            op.grad_component(uc, d, gd);
            acc += gd;                                       // same degree m-1
        }
        for (int r = 0; r < nloc; ++r)
            out.coeffs()(e * nloc + r, 0) = acc.coeff(r);    // assignment
    }
    return out;
}

// ---------------------------------------------------------------------------
// div_norm_sq: || div u_h ||^2_{L2(Omega)} (no broken_space needed; interval
// T yields a rigorous enclosure)
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
T div_norm_sq(vfe_space<D, T, P, SP>& vs, const vfe_function<D, T, P>& u) {
    vs.validate(u);
    fe_space<D, T, P, SP>& fs = vs.scalar();
    const dofmap<D>& dmu = fs.dofs(u.degree());
    const int Nu = dmu.ndof();
    element_op<D, T, P> op;
    bpoly<D, T> uc, gd, acc;
    T s(0);
    for (int e = 0; e < vs.num_elements(); ++e) {            // element order (X9)
        op.set_geometry(vs.geometry(e));
        detail::sv_gather_component(dmu, Nu, e, u, 0, uc);
        op.grad_component(uc, 0, acc);
        for (int d = 1; d < D; ++d) {
            detail::sv_gather_component(dmu, Nu, e, u, d, uc);
            op.grad_component(uc, d, gd);
            acc += gd;
        }
        s += op.local_inner(acc, acc);
    }
    return s;
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_SV_ASSEMBLE_HPP
