// vcp/btfem/time_space.hpp
// btfem T2/T3: time_space<T> -- the continuous piecewise Bernstein spaces
// Vbar_k^q (time_bc_none) and V_k^q (time_bc_zero_left) on a time_mesh, their
// global dof numbering, the global time matrices and the physicalized local
// tables used by pg_space.
//
// Conforms to: BT-1 external design v1.1 (sections 4.1, 4.2) and
//              BT-1 internal design v1.0 (sections 3, 4.2, 4.3).
//
// Division count: zero (the reference tables are shared T-stage tables of
// detail/time_tables.hpp; physicalization multiplies by k or 1/k of the mesh).

#ifndef VCP_BTFEM_TIME_SPACE_HPP
#define VCP_BTFEM_TIME_SPACE_HPP

#include <vector>
#include <map>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>

#include <vcp/btfem/time_mesh.hpp>
#include <vcp/btfem/detail/time_tables.hpp>

namespace vcp {
namespace btfem {

enum time_bc { time_bc_none, time_bc_zero_left };

template <typename T>
class time_space {
public:
    typedef vcp::matrix<T, vcp::mats<T> > table_type;

    time_space(const time_mesh<T>& tm, int q, time_bc bc)
        : tm_(tm), q_(q), bc_(bc), m_(tm.num_intervals()), ndof_(0), l2g_(),
          local_mass_(), local_stiffness_(), wm_cache_() {
        if (q < 1)
            throw std::invalid_argument("btfem::time_space: q < 1");
        // l2g_t (internal 3.1)
        l2g_.assign(static_cast<std::size_t>(m_) * static_cast<std::size_t>(q + 1), -1);
        for (int l = 1; l <= m_; ++l) {
            for (int r = 0; r <= q; ++r) {
                int g = (l - 1) * q + r;
                if (bc == time_bc_zero_left) g -= 1;    // (1,0) -> -1
                l2g_[static_cast<std::size_t>(l - 1) * static_cast<std::size_t>(q + 1)
                     + static_cast<std::size_t>(r)] = g;
            }
        }
        ndof_ = (bc == time_bc_none) ? m_ * q + 1 : m_ * q;
        // physicalized local tables (internal 4.2): k * Mref, (1/k) * Sref
        const table_type& M = detail::time_tables_typed<T>::mass_ref(q, q);
        const table_type& S = detail::time_tables_typed<T>::stiffness_ref(q);
        local_mass_.zeros(q + 1, q + 1);
        local_stiffness_.zeros(q + 1, q + 1);
        const T& k = tm_.width();
        const T& ik = tm_.inv_width();
        for (int a = 0; a <= q; ++a) {
            for (int b = 0; b <= q; ++b) {
                local_mass_(a, b) = k * M(a, b);
                local_stiffness_(a, b) = ik * S(a, b);
            }
        }
    }

    int degree() const { return q_; }
    time_bc bc() const { return bc_; }
    int num_intervals() const { return m_; }
    int local_size() const { return q_ + 1; }
    int ndof() const { return ndof_; }
    // g(l, r), l = 1..m, r = 0..q; removed dof -> -1
    int global_dof(int l, int r) const {
        if (l < 1 || l > m_ || r < 0 || r > q_)
            throw std::invalid_argument("btfem::time_space::global_dof: (l, r) out of range");
        return l2g_[static_cast<std::size_t>(l - 1) * static_cast<std::size_t>(q_ + 1)
                    + static_cast<std::size_t>(r)];
    }
    const time_mesh<T>& mesh_ref() const { return tm_; }

    // ---- global time matrices (dense, ndof x ndof) ----
    template <typename P> vcp::matrix<T, P> mass() const { return scatter_square<P>(local_mass_); }
    template <typename P> vcp::matrix<T, P> stiffness() const { return scatter_square<P>(local_stiffness_); }
    template <typename P> vcp::matrix<T, P> dt_mass() const { return scatter_square<P>(local_dt_mass()); }
    // y_a = psi_a(t_end): T(1) at g(m, q), ndof x 1
    template <typename P> vcp::matrix<T, P> end_values() const {
        vcp::matrix<T, P> y;
        y.zeros(ndof_, 1);
        y(global_dof(m_, q_), 0) = T(1);
        return y;
    }
    // rows = this space, cols = other (same time_mesh required; D-9)
    template <typename P> vcp::matrix<T, P> mixed_mass(const time_space& other) const {
        require_same_mesh(other, "mixed_mass");
        const table_type& M = detail::time_tables_typed<T>::mass_ref(q_, other.q_);
        table_type loc;
        loc.zeros(q_ + 1, other.q_ + 1);
        const T& k = tm_.width();
        for (int a = 0; a <= q_; ++a)
            for (int b = 0; b <= other.q_; ++b)
                loc(a, b) = k * M(a, b);
        return scatter_mixed<P>(loc, other);
    }
    template <typename P> vcp::matrix<T, P> mixed_stiffness(const time_space& other) const {
        require_same_mesh(other, "mixed_stiffness");
        const table_type& S = detail::time_tables_typed<T>::mixed_stiffness_ref(q_, other.q_);
        table_type loc;
        loc.zeros(q_ + 1, other.q_ + 1);
        const T& ik = tm_.inv_width();
        for (int a = 0; a <= q_; ++a)
            for (int b = 0; b <= other.q_; ++b)
                loc(a, b) = ik * S(a, b);
        return scatter_mixed<P>(loc, other);
    }
    template <typename P> vcp::matrix<T, P> mixed_dt_mass(const time_space& other) const {
        require_same_mesh(other, "mixed_dt_mass");
        return scatter_mixed<P>(detail::time_tables_typed<T>::mixed_dt_mass_ref(q_, other.q_), other);
    }

    // ---- local tables on one sub-interval (physical factors applied) ----
    const table_type& local_mass() const { return local_mass_; }             // k Mref(q,q)
    const table_type& local_stiffness() const { return local_stiffness_; }   // (1/k) Sref(q)
    const table_type& local_dt_mass() const {                                // Bref(q)
        return detail::time_tables_typed<T>::dt_mass_ref(q_);
    }
    // rows sigma = 0..sdeg, cols a: (beta^{sdeg}_sigma, psi_a')_{J_l}
    const table_type& local_dt_load(int sdeg) const {
        return detail::time_tables_typed<T>::dt_load_ref(sdeg, q_);
    }
    // (beta^{sdeg}_sigma psi_b, psi_a')_{J_l}, (q+1) x (q+1), row a col b
    const table_type& local_weighted_dt_mass(int sdeg, int sigma) const {
        const std::vector<table_type>& W =
            detail::time_tables_typed<T>::weighted_dt_mass_ref(sdeg, q_);
        if (sigma < 0 || sigma > sdeg)
            throw std::invalid_argument("btfem::time_space::local_weighted_dt_mass: sigma out of range");
        return W[static_cast<std::size_t>(sigma)];
    }
    // k (beta^{sdeg}_sigma psi_b, psi_a)_{J_l}
    const table_type& local_weighted_mass(int sdeg, int sigma) const {
        if (sdeg < 0)
            throw std::invalid_argument("btfem::time_space::local_weighted_mass: sdeg < 0");
        if (sigma < 0 || sigma > sdeg)
            throw std::invalid_argument("btfem::time_space::local_weighted_mass: sigma out of range");
        typename std::map<int, std::vector<table_type> >::iterator it = wm_cache_.find(sdeg);
        if (it == wm_cache_.end()) {
            const std::vector<table_type>& U =
                detail::time_tables_typed<T>::weighted_mass_ref(sdeg, q_);
            std::vector<table_type> phys;
            phys.reserve(U.size());
            const T& k = tm_.width();
            for (std::size_t s = 0; s < U.size(); ++s) {
                table_type t;
                t.zeros(q_ + 1, q_ + 1);
                for (int a = 0; a <= q_; ++a)
                    for (int b = 0; b <= q_; ++b)
                        t(a, b) = k * U[s](a, b);
                phys.push_back(t);
            }
            it = wm_cache_.insert(std::make_pair(sdeg, phys)).first;
        }
        return it->second[static_cast<std::size_t>(sigma)];
    }

    bool same_mesh(const time_space& other) const { return tm_.same_as(other.tm_); }

private:
    time_mesh<T> tm_;
    int q_;
    time_bc bc_;
    int m_;
    int ndof_;
    std::vector<int> l2g_;
    table_type local_mass_;
    table_type local_stiffness_;
    mutable std::map<int, std::vector<table_type> > wm_cache_;   // sdeg -> k * Uref[sigma]

    void require_same_mesh(const time_space& other, const char* where) const {
        if (!same_mesh(other)) {
            std::string msg("btfem::time_space::");
            msg += where;
            msg += ": time spaces live on different time meshes";
            throw std::invalid_argument(msg);
        }
    }
    // scatter a (q+1)x(q+1) local table with l2g_t; order: l ascending, (a,b) row major
    template <typename P>
    vcp::matrix<T, P> scatter_square(const table_type& loc) const {
        vcp::matrix<T, P> out;
        out.zeros(ndof_, ndof_);
        for (int l = 1; l <= m_; ++l) {
            for (int a = 0; a <= q_; ++a) {
                int ga = global_dof(l, a);
                if (ga < 0) continue;
                for (int b = 0; b <= q_; ++b) {
                    int gb = global_dof(l, b);
                    if (gb < 0) continue;
                    out(ga, gb) += loc(a, b);
                }
            }
        }
        return out;
    }
    template <typename P>
    vcp::matrix<T, P> scatter_mixed(const table_type& loc, const time_space& other) const {
        vcp::matrix<T, P> out;
        out.zeros(ndof_, other.ndof_);
        for (int l = 1; l <= m_; ++l) {
            for (int a = 0; a <= q_; ++a) {
                int ga = global_dof(l, a);
                if (ga < 0) continue;
                for (int b = 0; b <= other.q_; ++b) {
                    int gb = other.global_dof(l, b);
                    if (gb < 0) continue;
                    out(ga, gb) += loc(a, b);
                }
            }
        }
        return out;
    }
};

} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_TIME_SPACE_HPP
