// vcp/btfem/pg_space.hpp
// btfem T4-T12: pg_space<D,T,P,SP> -- the fully discrete Petrov-Galerkin space
// V_k^q (x) S_h^n on a bfem fe_space: st_poly factories, matrix-group factors
// (T6), time-dependent coefficient terms (T7), nonlinear scheme material (T8),
// initial errors (T9), range enclosures (T10), sections / evaluation (T4).
//
// Conforms to: BT-1 external design v1.1 (sections 5.2, 6, 6.1-6.5) and
//              BT-1 internal design v1.0 (sections 5.3, 5.4, 6, 8).
//
// Accumulation order (normative): sub-interval l -> element e -> time rank
// sigma -> local indices; COO duplicates combined by bfem coo_buffer::combine.
// Same input, single thread: bit reproducible. The OpenMP path (bfem L3
// model: persistent per-thread element_op + workspace + coo_buffers,
// thread-order concatenation) parallelizes the element loop of the T7 matrix
// terms (weighted_dt_mass_terms, space-time convection_dt_terms); the
// sub-interval and sigma loops and the vector / scalar paths are sequential.
// No division appears in this header.

#ifndef VCP_BTFEM_PG_SPACE_HPP
#define VCP_BTFEM_PG_SPACE_HPP

#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/fe_function.hpp>
#include <vcp/bfem/poly_field.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/refine.hpp>
#include <vcp/bfem/detail/scalar_traits.hpp>

#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/time_space.hpp>
#include <vcp/btfem/st_poly.hpp>
#include <vcp/btfem/kron.hpp>

// VCP_BFEM_USE_OPENMP is defined by vcp/bfem/fe_space.hpp (shared switch:
// -DVCP_BFEM_NOMP / -DVCP_NOMP disable bfem and btfem together)
#if VCP_BFEM_USE_OPENMP
#  include <omp.h>
#endif

namespace vcp {
namespace btfem {

namespace detail {
// per-thread assembly workspace of the OpenMP element loops (bfem L3 model:
// one element_op + local scratch + COO buffers per thread; persistent across
// calls so that repeated assembly allocates nothing from the second call on)
template <int D, typename T, typename P>
struct pg_thread_ws {
    bfem::element_op<D, T, P> op;
    st_workspace<D, T> ws;
    st_local<D, T> u_loc, w_loc;
    std::array<st_local<D, T>, static_cast<std::size_t>(D)> b_loc;
    std::array<bfem::bpoly<D, T>, static_cast<std::size_t>(D)> bsig;
    vcp::matrix<T, P> loc;
    std::vector<bfem::detail::coo_buffer<T> > bufs;     // one per sigma
};
} // namespace detail

template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class pg_space {
    static_assert(D == 2 || D == 3, "btfem::pg_space: only D == 2 or D == 3");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef st_poly<D, T, P> function_type;
    typedef kron_term<T, P, SP> term_type;
    typedef bfem::fe_space<D, T, P, SP> space_type;
    typedef bfem::fe_function<D, T, P> section_type;

    // Vh: the spatial space (not owned; must outlive the pg_space, D-6)
    // msh: the mesh Vh was built from (not owned; poly_field::restrict_to, D-7)
    // tm : the time mesh (copied)
    pg_space(space_type& Vh, const bfem::mesh<D, T>& msh, const time_mesh<T>& tm)
        : Vh_(&Vh), msh_(&msh), tm_(tm), tspaces_(), op_(), ws_(),
          u_loc_(), w_loc_(), b_loc_(), bsig_(), loc_(), buf_() {
        bfem::bfem_scalar_traits<T>::require();
        if (msh.num_elements() != Vh.num_elements())
            throw std::invalid_argument("btfem::pg_space: mesh and fe_space element counts differ");
    }

    // ---- T12 accessors ----
    space_type& space() { return *Vh_; }
    const bfem::mesh<D, T>& mesh_ref() const { return *msh_; }
    const time_mesh<T>& time_mesh_ref() const { return tm_; }
    const time_space<T>& time_dofs(int q, time_bc bc) {
        std::pair<int, int> key(q, static_cast<int>(bc));
        typename std::map<std::pair<int, int>, time_space<T> >::iterator it = tspaces_.find(key);
        if (it != tspaces_.end()) return it->second;
        return tspaces_.insert(std::make_pair(key, time_space<T>(tm_, q, bc))).first->second;
    }

    // ---- factories and vector conversion (external 5.2) ----
    function_type zero_function(int q, time_bc bc, int n) {
        function_type u = make_fn(q, bc, n);
        u.C_.zeros(u.mt_, u.nx_);
        return u;
    }
    function_type function_from_coeffs(int q, time_bc bc, int n, vcp::matrix<T, P> C) {
        function_type u = make_fn(q, bc, n);
        if (C.rowsize() != u.mt_ || C.columnsize() != u.nx_)
            throw std::invalid_argument("btfem::pg_space::function_from_coeffs: size != ndof_t x ndof_x");
        u.C_ = std::move(C);
        return u;
    }
    function_type function_from_vector(int q, time_bc bc, int n, const vcp::matrix<T, P>& v) {
        function_type u = make_fn(q, bc, n);
        if (v.rowsize() != u.mt_ * u.nx_ || v.columnsize() != 1)
            throw std::invalid_argument("btfem::pg_space::function_from_vector: size != ndof_t ndof_x x 1");
        u.C_.zeros(u.mt_, u.nx_);
        for (int a = 0; a < u.mt_; ++a)
            for (int i = 0; i < u.nx_; ++i)
                u.C_(a, i) = v(a * u.nx_ + i, 0);
        return u;
    }
    vcp::matrix<T, P> to_vector(const function_type& u) const {
        validate(u);
        vcp::matrix<T, P> v;
        v.zeros(u.mt_ * u.nx_, 1);
        for (int a = 0; a < u.mt_; ++a)
            for (int i = 0; i < u.nx_; ++i)
                v(a * u.nx_ + i, 0) = u.C_(a, i);
        return v;
    }

    // ---- T6: matrix-group factors, trial = test = V_k^q (x) S_h^n ----
    term_type A(int q, int n) { return term_type(T(1), time_dofs(q, time_bc_zero_left).template stiffness<P>(), Vh_->mixed_mass(n, n)); }
    term_type M(int q, int n) { return term_type(T(1), time_dofs(q, time_bc_zero_left).template mass<P>(), Vh_->stiffness(n)); }
    term_type B(int q, int n) { return term_type(T(1), time_dofs(q, time_bc_zero_left).template dt_mass<P>(), Vh_->stiffness(n)); }
    term_type U(int q, int n) { return term_type(T(1), time_dofs(q, time_bc_zero_left).template mass<P>(), Vh_->mixed_mass(n, n)); }
    term_type W(int q, int n) { return term_type(T(1), time_dofs(q, time_bc_zero_left).template stiffness<P>(), Vh_->stiffness(n)); }
    term_type Y(int q, int n) {
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        vcp::matrix<T, P> y = Vt.template end_values<P>();
        vcp::matrix<T, P> yy;
        yy.zeros(Vt.ndof(), Vt.ndof());
        for (int a = 0; a < Vt.ndof(); ++a)
            for (int b = 0; b < Vt.ndof(); ++b)
                yy(a, b) = y(a, 0) * y(b, 0);
        return term_type(T(1), yy, Vh_->mixed_mass(n, n));
    }
    // rows = V_k^q (bc_zero_left), cols = the time space of ubar; nbar == n (D-3)
    term_type A_mixed(int q, int n, const function_type& ubar) {
        validate(ubar);
        require_same_space_degree(ubar, n, "A_mixed");
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        return term_type(T(1), Vt.template mixed_stiffness<P>(Vtb), Vh_->mixed_mass(n, n));
    }
    term_type B_mixed(int q, int n, const function_type& ubar) {
        validate(ubar);
        require_same_space_degree(ubar, n, "B_mixed");
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        return term_type(T(1), Vt.template mixed_dt_mass<P>(Vtb), Vh_->stiffness(n));
    }

    // ---- T7: (w phi_j, dt phi_i), w = fprime[ubar]; m (sdeg + 1) terms in (l, sigma) order ----
    std::vector<term_type> weighted_dt_mass_terms(const bfem::poly1<T>& fprime,
                                                  const function_type& ubar, int q, int n) {
        validate(ubar);
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        const bfem::dofmap<D>& dm = Vh_->dofs(n);
        const bfem::dofmap<D>& dmb = Vh_->dofs(ubar.n_);
        const int sdeg = fprime.degree() * ubar.q_;
        const int m = tm_.num_intervals();
        const int nt = Vh_->num_elements();
        const int N = dm.local_size();
        const int nx = dm.ndof();
        std::vector<term_type> terms;
        terms.reserve(static_cast<std::size_t>(m) * static_cast<std::size_t>(sdeg + 1));
        prepare_buffers(sdeg + 1, nt, N);
        for (int l = 1; l <= m; ++l) {
            for (int s = 0; s <= sdeg; ++s) buf_[static_cast<std::size_t>(s)].clear();
            weighted_rows_loop(fprime, ubar, Vtb, dmb, dm, l, n, nt, N, sdeg);
            for (int s = 0; s <= sdeg; ++s) {
                bfem::detail::coo_buffer<T>& b = buf_[static_cast<std::size_t>(s)];
                b.combine();
                terms.push_back(term_type(T(1), scatter_time(Vt.local_weighted_dt_mass(sdeg, s), Vt, l),
                                          bfem::detail::spm_adapter<T, SP>::build(nx, nx, b)));
            }
        }
        return terms;
    }
    // constant b: the single term B_t (x) C_x(b)
    std::vector<term_type> convection_dt_terms(const std::array<T, D>& b_const, int q, int n) {
        std::vector<term_type> terms;
        terms.push_back(term_type(T(1), time_dofs(q, time_bc_zero_left).template dt_mass<P>(),
                                  Vh_->convection(b_const, n)));
        return terms;
    }
    // space-time b: m (q_b + 1) terms, time factor local_weighted_dt_mass(q_b, sigma)
    std::vector<term_type> convection_dt_terms(const std::array<function_type, D>& b, int q, int n) {
        validate_b(b);
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(b[0].q_, b[0].bc_);
        const bfem::dofmap<D>& dm = Vh_->dofs(n);
        const int qb = b[0].q_;
        const int m = tm_.num_intervals();
        const int nt = Vh_->num_elements();
        const int N = dm.local_size();
        const int nx = dm.ndof();
        std::vector<term_type> terms;
        terms.reserve(static_cast<std::size_t>(m) * static_cast<std::size_t>(qb + 1));
        prepare_buffers(qb + 1, nt, N);
        for (int l = 1; l <= m; ++l) {
            for (int s = 0; s <= qb; ++s) buf_[static_cast<std::size_t>(s)].clear();
            convection_rows_loop(b, Vtb, dm, l, n, nt, N, qb);
            for (int s = 0; s <= qb; ++s) {
                bfem::detail::coo_buffer<T>& bb = buf_[static_cast<std::size_t>(s)];
                bb.combine();
                terms.push_back(term_type(T(1), scatter_time(Vt.local_weighted_dt_mass(qb, s), Vt, l),
                                          bfem::detail::spm_adapter<T, SP>::build(nx, nx, bb)));
            }
        }
        return terms;
    }

    // ---- T8: vectors over the test dofs of V_k^q (x) S_h^n ----
    // (g(ubar), dt phi_i)
    vcp::matrix<T, P> load_dt(const bfem::poly1<T>& g, const function_type& ubar, int q, int n) {
        validate(ubar);
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        const bfem::dofmap<D>& dm = Vh_->dofs(n);
        const bfem::dofmap<D>& dmb = Vh_->dofs(ubar.n_);
        const int sdeg = g.degree() * ubar.q_;
        vcp::matrix<T, P> F;
        F.zeros(Vt.ndof() * dm.ndof(), 1);
        const int m = tm_.num_intervals();
        const int nt = Vh_->num_elements();
        for (int l = 1; l <= m; ++l) {
            for (int e = 0; e < nt; ++e) {
                op_.set_geometry(Vh_->geometry(e));
                detail::st_gather(u_loc_, ubar.C_, Vtb, l, dmb, e);
                detail::st_compose_into(w_loc_, g, u_loc_, ws_);
                accumulate_dt_load(F, Vt, dm, l, e, n, sdeg);
            }
        }
        return F;
    }
    // (b . grad ubar, dt phi_i), constant b
    vcp::matrix<T, P> convection_dt(const std::array<T, D>& b_const, const function_type& ubar, int q, int n) {
        validate(ubar);
        if (ubar.n_ < 1)
            throw std::invalid_argument("btfem::pg_space::convection_dt: ubar space degree < 1");
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        const bfem::dofmap<D>& dm = Vh_->dofs(n);
        const bfem::dofmap<D>& dmb = Vh_->dofs(ubar.n_);
        const int sdeg = ubar.q_;
        for (int d = 0; d < D; ++d) {
            b_loc_[static_cast<std::size_t>(d)].prepare(0, 0, false);
            bfem::detail::bpoly_access::vec(b_loc_[static_cast<std::size_t>(d)].r(0))[0] =
                b_const[static_cast<std::size_t>(d)];
        }
        vcp::matrix<T, P> F;
        F.zeros(Vt.ndof() * dm.ndof(), 1);
        const int m = tm_.num_intervals();
        const int nt = Vh_->num_elements();
        for (int l = 1; l <= m; ++l) {
            for (int e = 0; e < nt; ++e) {
                op_.set_geometry(Vh_->geometry(e));
                detail::st_gather(u_loc_, ubar.C_, Vtb, l, dmb, e);
                detail::st_b_dot_grad_into(w_loc_, b_loc_.data(), u_loc_, ws_, op_);
                accumulate_dt_load(F, Vt, dm, l, e, n, sdeg);
            }
        }
        return F;
    }
    // (b . grad ubar, dt phi_i), space-time b
    vcp::matrix<T, P> convection_dt(const std::array<function_type, D>& b, const function_type& ubar, int q, int n) {
        validate(ubar);
        validate_b(b);
        if (ubar.n_ < 1)
            throw std::invalid_argument("btfem::pg_space::convection_dt: ubar space degree < 1");
        const time_space<T>& Vt = time_dofs(q, time_bc_zero_left);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        const time_space<T>& Vtbb = time_dofs(b[0].q_, b[0].bc_);
        const bfem::dofmap<D>& dm = Vh_->dofs(n);
        const bfem::dofmap<D>& dmb = Vh_->dofs(ubar.n_);
        const int sdeg = b[0].q_ + ubar.q_;
        vcp::matrix<T, P> F;
        F.zeros(Vt.ndof() * dm.ndof(), 1);
        const int m = tm_.num_intervals();
        const int nt = Vh_->num_elements();
        for (int l = 1; l <= m; ++l) {
            for (int e = 0; e < nt; ++e) {
                op_.set_geometry(Vh_->geometry(e));
                detail::st_gather(u_loc_, ubar.C_, Vtb, l, dmb, e);
                for (int d = 0; d < D; ++d)
                    detail::st_gather(b_loc_[static_cast<std::size_t>(d)], b[static_cast<std::size_t>(d)].C_,
                                      Vtbb, l, Vh_->dofs(b[static_cast<std::size_t>(d)].n_), e);
                detail::st_b_dot_grad_into(w_loc_, b_loc_.data(), u_loc_, ws_, op_);
                accumulate_dt_load(F, Vt, dm, l, e, n, sdeg);
            }
        }
        return F;
    }

    // ---- T4 operations (external 5.2) ----
    // spatial section at the node tau_l (l = 0..m): a block row copy, no evaluation
    section_type section_at_node(const function_type& u, int l) {
        validate(u);
        const int m = tm_.num_intervals();
        if (l < 0 || l > m)
            throw std::invalid_argument("btfem::pg_space::section_at_node: l out of range");
        const time_space<T>& Vt = time_dofs(u.q_, u.bc_);
        int gt = (l == 0) ? Vt.global_dof(1, 0) : Vt.global_dof(l, u.q_);
        if (gt < 0) return Vh_->zero_function(u.n_);            // bc_zero_left, l = 0
        vcp::matrix<T, P> c;
        c.zeros(u.nx_, 1);
        for (int i = 0; i < u.nx_; ++i) c(i, 0) = u.C_(gt, i);
        return Vh_->function_from_coeffs(u.n_, c);
    }
    // spatial section inside sub-interval l (1..m) at local barycentric tau =
    // (1 - t, t): sum_r B^q_r(tau) C[g(l, r), .] (de Casteljau weights, r ascending)
    section_type section_at(const function_type& u, int l, const bfem::bary_point<1, T>& tau) {
        validate(u);
        const int m = tm_.num_intervals();
        if (l < 1 || l > m)
            throw std::invalid_argument("btfem::pg_space::section_at: l out of range");
        const time_space<T>& Vt = time_dofs(u.q_, u.bc_);
        vcp::matrix<T, P> c;
        c.zeros(u.nx_, 1);
        std::vector<T> unit(static_cast<std::size_t>(u.q_ + 1), T(0));
        for (int r = 0; r <= u.q_; ++r) {
            int gt = Vt.global_dof(l, r);
            if (gt < 0) continue;
            unit[static_cast<std::size_t>(r)] = T(1);
            const T w = bfem::eval(bfem::bpoly<1, T>::from_coeffs(u.q_, unit), tau);
            unit[static_cast<std::size_t>(r)] = T(0);
            for (int i = 0; i < u.nx_; ++i) c(i, 0) += w * u.C_(gt, i);
        }
        return Vh_->function_from_coeffs(u.n_, c);
    }
    // time trace at the spatial point (e, lam): one bpoly<1,T> per sub-interval
    std::vector<bfem::bpoly<1, T> > time_trace(const function_type& u, int e, const bfem::bary_point<D, T>& lam) {
        validate(u);
        if (e < 0 || e >= Vh_->num_elements())
            throw std::invalid_argument("btfem::pg_space::time_trace: element out of range");
        const time_space<T>& Vtb = time_dofs(u.q_, u.bc_);
        const bfem::dofmap<D>& dmb = Vh_->dofs(u.n_);
        const int m = tm_.num_intervals();
        std::vector<bfem::bpoly<1, T> > out;
        out.reserve(static_cast<std::size_t>(m));
        std::vector<T> vals(static_cast<std::size_t>(u.q_ + 1));
        for (int l = 1; l <= m; ++l) {
            detail::st_gather(u_loc_, u.C_, Vtb, l, dmb, e);
            for (int r = 0; r <= u.q_; ++r) vals[static_cast<std::size_t>(r)] = bfem::eval(u_loc_.r(r), lam);
            out.push_back(bfem::bpoly<1, T>::from_coeffs(u.q_, vals));
        }
        return out;
    }
    // point evaluation (sub-interval l, tau; element e, lam)
    T eval(const function_type& u, int l, const bfem::bary_point<1, T>& tau, int e, const bfem::bary_point<D, T>& lam) {
        validate(u);
        if (l < 1 || l > tm_.num_intervals())
            throw std::invalid_argument("btfem::pg_space::eval: l out of range");
        if (e < 0 || e >= Vh_->num_elements())
            throw std::invalid_argument("btfem::pg_space::eval: element out of range");
        const time_space<T>& Vtb = time_dofs(u.q_, u.bc_);
        detail::st_gather(u_loc_, u.C_, Vtb, l, Vh_->dofs(u.n_), e);
        return detail::st_eval(u_loc_, tau, lam, ws_);
    }
    // T10: coefficient hull (enclosure of the range; bound_traits<T>)
    bfem::range_pair<T> range(const function_type& u) const {
        validate(u);
        bfem::range_pair<T> acc;
        acc.lo = u.C_(0, 0);
        acc.hi = u.C_(0, 0);
        for (int a = 0; a < u.mt_; ++a)
            for (int i = 0; i < u.nx_; ++i) {
                if (a == 0 && i == 0) continue;
                bfem::bound_traits<T>::min_update(acc.lo, u.C_(a, i));
                bfem::bound_traits<T>::max_update(acc.hi, u.C_(a, i));
            }
        return acc;
    }
    // T10: hull of the local tensor coefficients of f(ubar) over all (l, e)
    bfem::range_pair<T> range_composed(const bfem::poly1<T>& f, const function_type& ubar) {
        validate(ubar);
        const time_space<T>& Vtb = time_dofs(ubar.q_, ubar.bc_);
        const bfem::dofmap<D>& dmb = Vh_->dofs(ubar.n_);
        bfem::range_pair<T> acc;
        bool started = false;
        const int m = tm_.num_intervals();
        const int nt = Vh_->num_elements();
        for (int l = 1; l <= m; ++l)
            for (int e = 0; e < nt; ++e) {
                detail::st_gather(u_loc_, ubar.C_, Vtb, l, dmb, e);
                detail::st_compose_into(w_loc_, f, u_loc_, ws_);
                detail::st_range_merge(acc, started, w_loc_);
            }
        return acc;
    }
    // spatial degree elevation of every time block (bfem fe_space::elevate)
    function_type elevate_space(const function_type& u, int n2) {
        validate(u);
        if (n2 < u.n_)
            throw std::invalid_argument("btfem::pg_space::elevate_space: n2 < space degree");
        function_type out = zero_function(u.q_, u.bc_, n2);
        vcp::matrix<T, P> c;
        for (int a = 0; a < u.mt_; ++a) {
            c.zeros(u.nx_, 1);
            for (int i = 0; i < u.nx_; ++i) c(i, 0) = u.C_(a, i);
            section_type blk = Vh_->elevate(Vh_->function_from_coeffs(u.n_, c), n2);
            for (int i = 0; i < out.nx_; ++i) out.C_(a, i) = blk.coeffs()(i, 0);
        }
        return out;
    }

    // ---- T9: initial errors ||u0 - ubar(t_0)||^2 (L2 and H_0^1 seminorm) ----
    T initial_error_sq_l2(const bfem::fe_function<D, T, P>& u0, const function_type& ubar) {
        section_type diff = initial_diff(u0, ubar);
        return Vh_->inner(diff, diff);
    }
    T initial_error_sq_h01(const bfem::fe_function<D, T, P>& u0, const function_type& ubar) {
        section_type diff = initial_diff(u0, ubar);
        const bfem::dofmap<D>& dm = Vh_->dofs(diff.degree());
        T acc(0);
        const int nt = Vh_->num_elements();
        for (int e = 0; e < nt; ++e) {
            op_.set_geometry(Vh_->geometry(e));
            gather_fe(u_loc_.r(0), diff, dm, e);
            for (int d = 0; d < D; ++d) {
                op_.grad_component(u_loc_.r(0), d, ws_.sb2);
                acc += op_.local_inner(ws_.sb2, ws_.sb2);
            }
        }
        return acc;
    }
    T initial_error_sq_l2(const bfem::poly_field<D, T>& u0, const function_type& ubar) {
        return initial_error_pf(u0, ubar, false);
    }
    T initial_error_sq_h01(const bfem::poly_field<D, T>& u0, const function_type& ubar) {
        return initial_error_pf(u0, ubar, true);
    }

private:
    space_type* Vh_;
    const bfem::mesh<D, T>* msh_;
    time_mesh<T> tm_;
    std::map<std::pair<int, int>, time_space<T> > tspaces_;
    bfem::element_op<D, T, P> op_;
    detail::st_workspace<D, T> ws_;
    detail::st_local<D, T> u_loc_, w_loc_;
    std::array<detail::st_local<D, T>, static_cast<std::size_t>(D)> b_loc_;
    std::array<bfem::bpoly<D, T>, static_cast<std::size_t>(D)> bsig_;
    vcp::matrix<T, P> loc_;
    std::vector<bfem::detail::coo_buffer<T> > buf_;

    function_type make_fn(int q, time_bc bc, int n) {
        const time_space<T>& Vt = time_dofs(q, bc);
        const int nx = Vh_->ndof(n);
        return function_type(q, bc, n, Vt.ndof(), nx, static_cast<const void*>(Vh_),
                             static_cast<const void*>(&tm_));
    }
    void validate(const function_type& u) const {
        if (u.space_id_ != static_cast<const void*>(Vh_) || u.tmesh_id_ != static_cast<const void*>(&tm_))
            throw std::invalid_argument("btfem::pg_space: st_poly does not belong to this pg_space");
        if (u.C_.rowsize() != u.mt_ || u.C_.columnsize() != u.nx_)
            throw std::invalid_argument("btfem::pg_space: st_poly coefficient matrix size changed");
    }
    void validate_b(const std::array<function_type, D>& b) const {
        for (int d = 0; d < D; ++d) {
            validate(b[static_cast<std::size_t>(d)]);
            if (b[static_cast<std::size_t>(d)].q_ != b[0].q_ || b[static_cast<std::size_t>(d)].bc_ != b[0].bc_)
                throw std::invalid_argument("btfem::pg_space: components of b must share the time space (q, bc)");
        }
    }
    void require_same_space_degree(const function_type& ubar, int n, const char* where) const {
        if (ubar.n_ != n) {
            std::string msg("btfem::pg_space::");
            msg += where;
            msg += ": ubar space degree != n (D-3; use elevate_space)";
            throw std::invalid_argument(msg);
        }
    }
    void prepare_buffers(int count, int nt, int N) {
        if (buf_.size() < static_cast<std::size_t>(count)) buf_.resize(static_cast<std::size_t>(count));
        for (int s = 0; s < count; ++s)
            buf_[static_cast<std::size_t>(s)].reserve(static_cast<std::size_t>(nt)
                                                      * static_cast<std::size_t>(N)
                                                      * static_cast<std::size_t>(N));
    }
    // dense ndof_t x ndof_t time factor of ONE sub-interval l (assignment only)
    vcp::matrix<T, P> scatter_time(const typename time_space<T>::table_type& loc,
                                   const time_space<T>& Vt, int l) const {
        vcp::matrix<T, P> out;
        out.zeros(Vt.ndof(), Vt.ndof());
        for (int a = 0; a <= Vt.degree(); ++a) {
            int ga = Vt.global_dof(l, a);
            if (ga < 0) continue;
            for (int b = 0; b <= Vt.degree(); ++b) {
                int gb = Vt.global_dof(l, b);
                if (gb < 0) continue;
                out(ga, gb) = loc(a, b);
            }
        }
        return out;
    }
    // element loop of weighted_dt_mass_terms for one l: compose once per
    // element, scatter row sigma into buf_[sigma] (internal 6.1)
    void weighted_rows_loop(const bfem::poly1<T>& fprime, const function_type& ubar,
                            const time_space<T>& Vtb, const bfem::dofmap<D>& dmb,
                            const bfem::dofmap<D>& dm, int l, int n, int nt, int N, int sdeg) {
#if VCP_BFEM_USE_OPENMP
        int nrun = 1;
#pragma omp parallel
        {
#pragma omp single
            {
                nrun = omp_get_num_threads();
                prepare_threads(nrun, sdeg + 1, nt / nrun, N);
            }
#pragma omp barrier
            const int tid = omp_get_thread_num();
            const int e0 = static_cast<int>((static_cast<long long>(nt) * tid) / nrun);
            const int e1 = static_cast<int>((static_cast<long long>(nt) * (tid + 1)) / nrun);
            detail::pg_thread_ws<D, T, P>& w = tws_[static_cast<std::size_t>(tid)];
            for (int e = e0; e < e1; ++e) {
                w.op.set_geometry(Vh_->geometry(e));
                detail::st_gather(w.u_loc, ubar.C_, Vtb, l, dmb, e);
                detail::st_compose_into(w.w_loc, fprime, w.u_loc, w.ws);
                for (int s = 0; s <= sdeg; ++s) {
                    w.op.local_weighted_mass(w.w_loc.r(s), n, w.loc);
                    bfem::detail::scatter_matrix<T>(dm, e, w.loc, N, w.bufs[static_cast<std::size_t>(s)],
                                                    typename bfem::dofmap<D>::family_tag());
                }
            }
        }
        merge_threads(nrun, sdeg + 1);
#else
        for (int e = 0; e < nt; ++e) {
            op_.set_geometry(Vh_->geometry(e));
            detail::st_gather(u_loc_, ubar.C_, Vtb, l, dmb, e);
            detail::st_compose_into(w_loc_, fprime, u_loc_, ws_);
            assert(w_loc_.qt == sdeg);
            for (int s = 0; s <= sdeg; ++s) {
                op_.local_weighted_mass(w_loc_.r(s), n, loc_);
                bfem::detail::scatter_matrix<T>(dm, e, loc_, N, buf_[static_cast<std::size_t>(s)],
                                                typename bfem::dofmap<D>::family_tag());
            }
        }
#endif
    }
    void convection_rows_loop(const std::array<function_type, D>& b, const time_space<T>& Vtb,
                              const bfem::dofmap<D>& dm, int l, int n, int nt, int N, int qb) {
        // dofmaps of the b components are materialized before the parallel region
        std::array<const bfem::dofmap<D>*, static_cast<std::size_t>(D)> dmbs;
        for (int d = 0; d < D; ++d) dmbs[static_cast<std::size_t>(d)] = &Vh_->dofs(b[static_cast<std::size_t>(d)].n_);
#if VCP_BFEM_USE_OPENMP
        int nrun = 1;
#pragma omp parallel
        {
#pragma omp single
            {
                nrun = omp_get_num_threads();
                prepare_threads(nrun, qb + 1, nt / nrun, N);
            }
#pragma omp barrier
            const int tid = omp_get_thread_num();
            const int e0 = static_cast<int>((static_cast<long long>(nt) * tid) / nrun);
            const int e1 = static_cast<int>((static_cast<long long>(nt) * (tid + 1)) / nrun);
            detail::pg_thread_ws<D, T, P>& w = tws_[static_cast<std::size_t>(tid)];
            for (int e = e0; e < e1; ++e) {
                w.op.set_geometry(Vh_->geometry(e));
                for (int d = 0; d < D; ++d)
                    detail::st_gather(w.b_loc[static_cast<std::size_t>(d)], b[static_cast<std::size_t>(d)].C_,
                                      Vtb, l, *dmbs[static_cast<std::size_t>(d)], e);
                for (int s = 0; s <= qb; ++s) {
                    for (int d = 0; d < D; ++d)
                        w.bsig[static_cast<std::size_t>(d)] = w.b_loc[static_cast<std::size_t>(d)].r(s);
                    w.op.local_convection(w.bsig, n, w.loc);
                    bfem::detail::scatter_matrix<T>(dm, e, w.loc, N, w.bufs[static_cast<std::size_t>(s)],
                                                    typename bfem::dofmap<D>::family_tag());
                }
            }
        }
        merge_threads(nrun, qb + 1);
#else
        for (int e = 0; e < nt; ++e) {
            op_.set_geometry(Vh_->geometry(e));
            for (int d = 0; d < D; ++d)
                detail::st_gather(b_loc_[static_cast<std::size_t>(d)], b[static_cast<std::size_t>(d)].C_,
                                  Vtb, l, *dmbs[static_cast<std::size_t>(d)], e);
            for (int s = 0; s <= qb; ++s) {
                for (int d = 0; d < D; ++d)
                    bsig_[static_cast<std::size_t>(d)] = b_loc_[static_cast<std::size_t>(d)].r(s);
                op_.local_convection(bsig_, n, loc_);
                bfem::detail::scatter_matrix<T>(dm, e, loc_, N, buf_[static_cast<std::size_t>(s)],
                                                typename bfem::dofmap<D>::family_tag());
            }
        }
#endif
    }
    // F(gt n_x + gdof) += L(sigma, a) * (w_sigma, phi_alpha)_e  (order sigma -> a -> alpha)
    void accumulate_dt_load(vcp::matrix<T, P>& F, const time_space<T>& Vt, const bfem::dofmap<D>& dm,
                            int l, int e, int n, int sdeg) {
        assert(w_loc_.qt == sdeg);
        const typename time_space<T>::table_type& L = Vt.local_dt_load(sdeg);
        const int nx = dm.ndof();
        const int N = dm.local_size();
        for (int s = 0; s <= sdeg; ++s) {
            op_.local_load(w_loc_.r(s), n, loc_);
            for (int a = 0; a <= Vt.degree(); ++a) {
                int gt = Vt.global_dof(l, a);
                if (gt < 0) continue;
                const T& tcoef = L(s, a);
                for (int al = 0; al < N; ++al)
                    F(gt * nx + dm.global_dof(e, al), 0) += tcoef * loc_(al, 0);
            }
        }
    }
    // gather of a bfem fe_function into a bpoly (identity P^n map, public dofmap API)
    void gather_fe(bfem::bpoly<D, T>& dst, const bfem::fe_function<D, T, P>& u, const bfem::dofmap<D>& dm, int e) const {
        bfem::detail::bpoly_access::prepare(dst, dm.degree(), false);
        std::vector<T>& c = bfem::detail::bpoly_access::vec(dst);
        for (int a = 0; a < dm.local_size(); ++a) c[static_cast<std::size_t>(a)] = u.coeffs()(dm.global_dof(e, a), 0);
    }
    // u0 - ubar(t_0) as an fe_function of degree max(deg u0, nbar)
    section_type initial_diff(const bfem::fe_function<D, T, P>& u0, const function_type& ubar) {
        validate(ubar);
        if (u0.coeffs().columnsize() != 1 || u0.coeffs().rowsize() != Vh_->ndof(u0.degree()))
            throw std::invalid_argument("btfem::pg_space::initial_error: u0 does not match the fe_space");
        section_type e0 = section_at_node(ubar, 0);
        const int nn = u0.degree() > e0.degree() ? u0.degree() : e0.degree();
        section_type a = (u0.degree() == nn) ? u0 : Vh_->elevate(u0, nn);
        section_type b = (e0.degree() == nn) ? e0 : Vh_->elevate(e0, nn);
        vcp::matrix<T, P> c;
        c.zeros(Vh_->ndof(nn), 1);
        for (int i = 0; i < c.rowsize(); ++i) c(i, 0) = a.coeffs()(i, 0) - b.coeffs()(i, 0);
        return Vh_->function_from_coeffs(nn, c);
    }
    // poly_field route: element-wise restrict_to vs the gathered section block
    T initial_error_pf(const bfem::poly_field<D, T>& u0, const function_type& ubar, bool h01) {
        validate(ubar);
        section_type e0 = section_at_node(ubar, 0);
        const bfem::dofmap<D>& dm = Vh_->dofs(e0.degree());
        const int nn = u0.total_degree() > e0.degree() ? u0.total_degree() : e0.degree();
        T acc(0);
        const int nt = Vh_->num_elements();
        for (int e = 0; e < nt; ++e) {
            op_.set_geometry(Vh_->geometry(e));
            bfem::bpoly<D, T> u0loc = u0.restrict_to(*msh_, e, nn);
            gather_fe(ws_.sb0, e0, dm, e);
            bfem::sub_into(ws_.sb1, u0loc, ws_.sb0);            // degree nn (lifts e0)
            if (!h01) {
                acc += op_.local_inner(ws_.sb1, ws_.sb1);
            } else {
                for (int d = 0; d < D; ++d) {
                    op_.grad_component(ws_.sb1, d, ws_.sb2);
                    acc += op_.local_inner(ws_.sb2, ws_.sb2);
                }
            }
        }
        return acc;
    }
#if VCP_BFEM_USE_OPENMP
    std::vector<detail::pg_thread_ws<D, T, P> > tws_;      // persistent per-thread workspaces
    void prepare_threads(int nthreads, int count, int e_per, int N) {
        if (tws_.size() < static_cast<std::size_t>(nthreads)) tws_.resize(static_cast<std::size_t>(nthreads));
        for (int t = 0; t < nthreads; ++t) {
            detail::pg_thread_ws<D, T, P>& w = tws_[static_cast<std::size_t>(t)];
            if (w.bufs.size() < static_cast<std::size_t>(count)) w.bufs.resize(static_cast<std::size_t>(count));
            for (int s = 0; s < count; ++s) {
                w.bufs[static_cast<std::size_t>(s)].clear();
                w.bufs[static_cast<std::size_t>(s)].reserve(static_cast<std::size_t>(e_per + 1)
                                                            * static_cast<std::size_t>(N)
                                                            * static_cast<std::size_t>(N));
            }
        }
    }
    // thread-order concatenation into buf_[s] (bfem append_all semantics)
    void merge_threads(int nthreads, int count) {
        for (int s = 0; s < count; ++s) {
            bfem::detail::coo_buffer<T>& dst = buf_[static_cast<std::size_t>(s)];
            std::size_t total = dst.a.size();
            for (int t = 0; t < nthreads; ++t) total += tws_[static_cast<std::size_t>(t)].bufs[static_cast<std::size_t>(s)].a.size();
            dst.reserve(total);
            for (int t = 0; t < nthreads; ++t) dst.append(tws_[static_cast<std::size_t>(t)].bufs[static_cast<std::size_t>(s)]);
        }
    }
#endif
};

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_PG_SPACE_HPP
