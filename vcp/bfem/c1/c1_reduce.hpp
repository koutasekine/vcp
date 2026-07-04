// vcp/bfem/c1/c1_reduce.hpp
// Phase 6 (2D C1 Argyris family): boundary constraint generation
// c1_boundary(kind) with the normative corner classification (v0.2, B-3),
// feeding the rational-coefficient linear_reduction overload (C1-6).
//
// Conforms to: C1 external design v0.2 (section 7) and
//              C1 internal design v0.2 (section 7).
//
// Corner classification (normative table of external design 7): a boundary
// vertex v with incident LISTED edges is classified by the exact
// collinearity of the edge tangents (certified sign of the 2 x 2 cross
// product, the sv_singular family of decisions):
//   1 direction class ("on a straight boundary"):
//       gradient: 1 row  t . grad u = 0           (rational coefficients)
//       Hessian:  1 row  t^T H t = 0
//   2 direction classes ("true corner"):
//       gradient: 2 rows == grad u = 0             (degenerates to simple
//                                                   eliminations -- B-3)
//       Hessian:  2 rows t_i^T H t_i = 0 (i = 1, 2; independent by the
//                 Veronese image of distinct directions)
//   >= 3 direction classes: gradient as the corner; Hessian: H = 0 (3
//       simple rows -- the t^T H t system reaches full rank);
//   interval-indeterminate collinearity: fail-fast with
//       c1_indeterminate_corner (SV-4 contract).
// clamped (u = du/dnu = 0) adds the edge ND rows, forces grad u = 0 at every
// listed vertex and extends the straight-boundary Hessian rows to
// { t^T H t = 0, t^T H nu = 0 } (corners: H = 0). simply-supported shares
// the ESSENTIAL row set with the homogeneous Dirichlet kind (the moment
// condition is natural, not essential -- recorded in the usage document).
//
// Coordinates: the classification decides on T through
// c1_coord_traits<T>::sign (certified; interval indeterminate throws); the
// row coefficients are EXACT rationals obtained through
// c1_coord_traits<T>::to_rational (coordinate differences -- the C1-6
// requirement basis). The trait is specialized here for detail::rational;
// other scalars opt in on the user/test side (the convert_traits pattern).
//
// No division operator appears in this file (G-C1-1); the rational-stage
// divisions of the elimination run inside sv/linear_reduction.hpp.

#ifndef VCP_BFEM_C1_C1_REDUCE_HPP
#define VCP_BFEM_C1_C1_REDUCE_HPP

#include <vector>
#include <map>
#include <algorithm>
#include <string>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/sv/linear_reduction.hpp>   // sv_constraint_q + reduction
#include <vcp/bfem/c1/c1_space.hpp>

namespace vcp {
namespace bfem {

// SV-4 contract transplanted: interval corner classification undecidable
class c1_indeterminate_corner : public std::runtime_error {
public:
    explicit c1_indeterminate_corner(const std::string& msg)
        : std::runtime_error(msg) {}
};

namespace detail {

// certified sign: -1 / 0 / +1, throw when not certified (interval
// straddling zero without being the exact point zero)
template <typename T>
int c1_sign_certified(const T& x) {
    if (x < T(0)) return -1;
    if (T(0) < x) return 1;
    if (x <= T(0) && T(0) <= x) return 0;
    throw c1_indeterminate_corner(
        "bfem::c1_boundary: collinearity undecidable on interval input");
}

} // namespace detail

// opt-in coordinate trait: sign certification on T + exact rational image
template <typename T>
struct c1_coord_traits;

template <>
struct c1_coord_traits<detail::rational> {
    static int sign(const detail::rational& x) {
        return detail::c1_sign_certified(x);
    }
    static const detail::rational& to_rational(const detail::rational& x) {
        return x;
    }
};

enum c1_boundary_kind {
    c1_bc_dirichlet,           // u = 0 (essential rows of the 2nd order case)
    c1_bc_simply_supported,    // same ESSENTIAL rows as dirichlet (recorded)
    c1_bc_clamped              // u = 0 and du/dnu = 0
};

namespace detail {

// push a rational row, pruning exact zero coefficients
inline void c1_push_row(std::vector<sv_constraint_q>& out,
                        const std::vector<int>& dof,
                        const std::vector<rational>& coef) {
    sv_constraint_q r;
    for (std::size_t i = 0; i < dof.size(); ++i) {
        if (coef[i].is_zero()) continue;
        r.dof.push_back(dof[i]);
        r.coef.push_back(coef[i]);
    }
    assert(!r.dof.empty());
    out.push_back(r);
}

inline void c1_push_single(std::vector<sv_constraint_q>& out, int dof) {
    sv_constraint_q r;
    r.dof.push_back(dof);
    r.coef.push_back(rational(1));
    out.push_back(r);
}

} // namespace detail

// ---------------------------------------------------------------------------
// c1_boundary: constraint rows of the given kind on the given edge list
// (default overload: the whole domain boundary). Deterministic emission
// order: edges ascending (trace rows, then ND rows for clamped), then
// vertices ascending (value, gradient rows, Hessian rows).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
std::vector<sv_constraint_q> c1_boundary(c1_space<D, T, P, SP>& sp, int m,
                                         c1_boundary_kind kind,
                                         const std::vector<int>& edges) {
    typedef detail::rational rat;
    typedef c1_coord_traits<T> CT;
    const c1_dofmap& dm = sp.dofs(m);
    const mesh<2, T>& msh = sp.mesh_ref();
    // deterministic, duplicate-free ascending edge list
    std::vector<int> eds = edges;
    std::sort(eds.begin(), eds.end());
    eds.erase(std::unique(eds.begin(), eds.end()), eds.end());
    for (std::size_t s = 0; s < eds.size(); ++s)
        if (eds[s] < 0 || eds[s] >= dm.num_edges())
            throw std::invalid_argument("bfem::c1_boundary: edge id out of range");
    std::vector<sv_constraint_q> out;
    // ---- edge rows ----
    for (std::size_t s = 0; s < eds.size(); ++s) {
        int ed = eds[s];
        for (int i = 0; i < detail::c1_ntrace(m); ++i)
            detail::c1_push_single(out, dm.edge_trace_dof(ed, i));
        if (kind == c1_bc_clamped)
            for (int j = 0; j < detail::c1_nnd(m); ++j)
                detail::c1_push_single(out, dm.edge_nd_dof(ed, j));
    }
    // ---- vertex incidence of the listed edges ----
    std::map<int, std::vector<int> > vinc;   // vertex -> incident listed edges
    for (std::size_t s = 0; s < eds.size(); ++s) {
        const std::array<int, 2>& vv = dm.edge_verts(eds[s]);
        vinc[vv[0]].push_back(eds[s]);
        vinc[vv[1]].push_back(eds[s]);
    }
    for (std::map<int, std::vector<int> >::const_iterator it = vinc.begin();
         it != vinc.end(); ++it) {
        const int v = it->first;
        // tangents of the incident listed edges (T for classification,
        // rational for the coefficients)
        std::vector<std::array<T, 2> > tt;
        std::vector<std::array<rat, 2> > tq;
        for (std::size_t s = 0; s < it->second.size(); ++s) {
            const std::array<int, 2>& vv = dm.edge_verts(it->second[s]);
            const std::array<T, 2>& a = msh.vertex(vv[0]);
            const std::array<T, 2>& b = msh.vertex(vv[1]);
            std::array<T, 2> t;
            t[0] = b[0] - a[0];
            t[1] = b[1] - a[1];
            tt.push_back(t);
        }
        // direction classes by certified collinearity (greedy over reps)
        std::vector<std::size_t> reps;
        for (std::size_t s = 0; s < tt.size(); ++s) {
            bool found = false;
            for (std::size_t r2 = 0; r2 < reps.size() && !found; ++r2) {
                const std::array<T, 2>& u = tt[reps[r2]];
                T cr = u[0] * tt[s][1] - u[1] * tt[s][0];
                if (CT::sign(cr) == 0) found = true;
            }
            if (!found) reps.push_back(s);
        }
        for (std::size_t r2 = 0; r2 < reps.size(); ++r2) {
            std::array<rat, 2> q;
            q[0] = CT::to_rational(tt[reps[r2]][0]);
            q[1] = CT::to_rational(tt[reps[r2]][1]);
            tq.push_back(q);
        }
        const int classes = static_cast<int>(tq.size());
        // value row
        detail::c1_push_single(out, dm.vertex_dof(v, 0));
        // gradient rows
        if (kind == c1_bc_clamped || classes >= 2) {
            // grad u = 0: the simple-elimination degeneration (B-3)
            detail::c1_push_single(out, dm.vertex_dof(v, 1));
            detail::c1_push_single(out, dm.vertex_dof(v, 2));
        } else {
            std::vector<int> dof;
            std::vector<rat> coef;
            dof.push_back(dm.vertex_dof(v, 1));
            dof.push_back(dm.vertex_dof(v, 2));
            coef.push_back(tq[0][0]);
            coef.push_back(tq[0][1]);
            detail::c1_push_row(out, dof, coef);
        }
        // Hessian rows
        const bool full_h = (kind == c1_bc_clamped) ? (classes >= 2)
                                                    : (classes >= 3);
        if (full_h) {
            detail::c1_push_single(out, dm.vertex_dof(v, 3));
            detail::c1_push_single(out, dm.vertex_dof(v, 4));
            detail::c1_push_single(out, dm.vertex_dof(v, 5));
        } else if (kind == c1_bc_clamped) {
            // straight boundary, clamped: t^T H t = 0 and t^T H nu = 0
            std::vector<int> dof;
            dof.push_back(dm.vertex_dof(v, 3));
            dof.push_back(dm.vertex_dof(v, 4));
            dof.push_back(dm.vertex_dof(v, 5));
            std::vector<rat> c1v;
            c1v.push_back(tq[0][0] * tq[0][0]);
            c1v.push_back(tq[0][0] * tq[0][1] + tq[0][0] * tq[0][1]);
            c1v.push_back(tq[0][1] * tq[0][1]);
            detail::c1_push_row(out, dof, c1v);
            // nu = R_{-90} t = (t1, -t0):
            // t^T H nu = t0 t1 Hxx + (t1^2 - t0^2) Hxy - t0 t1 Hyy
            std::vector<rat> c2v;
            c2v.push_back(tq[0][0] * tq[0][1]);
            c2v.push_back(tq[0][1] * tq[0][1] - tq[0][0] * tq[0][0]);
            c2v.push_back(-(tq[0][0] * tq[0][1]));
            detail::c1_push_row(out, dof, c2v);
        } else {
            // dirichlet family: t_i^T H t_i = 0 per direction class
            for (int r2 = 0; r2 < classes; ++r2) {
                std::vector<int> dof;
                dof.push_back(dm.vertex_dof(v, 3));
                dof.push_back(dm.vertex_dof(v, 4));
                dof.push_back(dm.vertex_dof(v, 5));
                std::vector<rat> cv;
                cv.push_back(tq[static_cast<std::size_t>(r2)][0]
                             * tq[static_cast<std::size_t>(r2)][0]);
                cv.push_back(tq[static_cast<std::size_t>(r2)][0]
                                 * tq[static_cast<std::size_t>(r2)][1]
                             + tq[static_cast<std::size_t>(r2)][0]
                                   * tq[static_cast<std::size_t>(r2)][1]);
                cv.push_back(tq[static_cast<std::size_t>(r2)][1]
                             * tq[static_cast<std::size_t>(r2)][1]);
                detail::c1_push_row(out, dof, cv);
            }
        }
    }
    return out;
}

template <int D, typename T, typename P, class SP>
std::vector<sv_constraint_q> c1_boundary(c1_space<D, T, P, SP>& sp, int m,
                                         c1_boundary_kind kind) {
    return c1_boundary(sp, m, kind, sp.dofs(m).boundary_edge_ids());
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_C1_C1_REDUCE_HPP
