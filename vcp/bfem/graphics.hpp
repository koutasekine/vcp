// vcp/bfem/graphics.hpp
// GRF-1: lattice sampling of finite element functions for graphics output.
//
// Conforms to: GRF-1 external design v1.0.
//
// One free-function family output_uh_for_graphics(space, u [, elems], div)
// covering the five plottable spaces (fe_space<2>/<3>, c1_space<2>,
// vfe_space, broken_space; rt_space is intentionally NOT covered -- GRF-1
// design section 2). Each call returns a graphics_output package:
//
//   points : (nsel * npt) x (D + ncomp)   rows = sample points, columns =
//            physical coordinates then value components
//   cells  : (nsel * ncell1) x (D + 1)    0-based row indices into points
//
// npt = C(div + D, D) lattice points per element, ncell1 = div^D sub-cells
// per element. Points shared by neighbouring elements are emitted once per
// incident element ON PURPOSE (discontinuous fields keep their jumps; the
// consumer may merge with a tolerance -- design section 6).
//
// Interval scalars: every emitted interval encloses the corresponding exact
// value (plain interval arithmetic; verified against the exact rational run
// in the GRF-1 gate). A row (x-box, u-box) encloses the sample AT ITS
// LATTICE POINT; it is NOT a claim about all points inside the x-box.

#ifndef VCP_BFEM_GRAPHICS_HPP
#define VCP_BFEM_GRAPHICS_HPP

#include <vector>
#include <array>
#include <stdexcept>

#include <vcp/matrix.hpp>

#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/c1/c1_space.hpp>
#include <vcp/bfem/sv/vfe_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// graphics_output<D,T,P>: the (points, cells) package
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T> >
struct graphics_output {
    vcp::matrix<T, P> points;   // npoint x (D + ncomp)
    vcp::matrix<int>  cells;    // ncell x (D + 1), 0-based rows of points
    int ncomp;                  // value components (scalar 1, vector D)
    int div;                    // lattice subdivision used

    graphics_output() : points(), cells(), ncomp(0), div(0) {}

    int num_points() const { return points.rowsize(); }
    int num_cells()  const { return cells.rowsize(); }
};

namespace detail {

// ---- per-element lattice cells (rank indices in index_map<D>(div)) --------
//
// D == 2: div^2 triangles: upward (alpha+e0, alpha+e1, alpha+e2) for
// |alpha| = div-1 and downward (alpha+e0+e1, alpha+e1+e2, alpha+e0+e2) for
// |alpha| = div-2.
//
// D == 3: div^3 tetrahedra by the principal-lattice decomposition:
//   corner   (|alpha| = div-1): alpha + {e0, e1, e2, e3}
//   octahedron (|alpha| = div-2): 6 vertices alpha+ei+ej split by the
//       diagonal (alpha+e0+e1, alpha+e2+e3) into 4 tetrahedra
//   inverted (|alpha| = div-3): alpha + {f0, f1, f2, f3}, fi = 1 - ei
// [出典未逐語確認] -- accepted through the machine audit in the GRF-1 gate
// (exact rational volume identity + face incidence), not through a citation.

inline void lattice_cells(const index_map<2>& im, int div,
                          std::vector<std::array<int, 3> >& out) {
    out.clear();
    out.reserve(static_cast<std::size_t>(div) * static_cast<std::size_t>(div));
    std::array<int, 3> c;
    {
        const index_map<2> base(div - 1);
        for (int r = 0; r < base.size(); ++r) {
            const multi_index<2> al = base.unrank(r);
            for (int k = 0; k < 3; ++k) {
                multi_index<2> v = al;
                v.a[static_cast<std::size_t>(k)] += 1;
                c[static_cast<std::size_t>(k)] = im.rank(v);
            }
            out.push_back(c);
        }
    }
    if (div >= 2) {
        const index_map<2> base(div - 2);
        const int pair_i[3] = { 0, 1, 0 };
        const int pair_j[3] = { 1, 2, 2 };
        for (int r = 0; r < base.size(); ++r) {
            const multi_index<2> al = base.unrank(r);
            for (int k = 0; k < 3; ++k) {
                multi_index<2> v = al;
                v.a[static_cast<std::size_t>(pair_i[k])] += 1;
                v.a[static_cast<std::size_t>(pair_j[k])] += 1;
                c[static_cast<std::size_t>(k)] = im.rank(v);
            }
            out.push_back(c);
        }
    }
}

inline void lattice_cells(const index_map<3>& im, int div,
                          std::vector<std::array<int, 4> >& out) {
    out.clear();
    out.reserve(static_cast<std::size_t>(div) * static_cast<std::size_t>(div)
                * static_cast<std::size_t>(div));
    std::array<int, 4> c;
    {   // corner tetrahedra
        const index_map<3> base(div - 1);
        for (int r = 0; r < base.size(); ++r) {
            const multi_index<3> al = base.unrank(r);
            for (int k = 0; k < 4; ++k) {
                multi_index<3> v = al;
                v.a[static_cast<std::size_t>(k)] += 1;
                c[static_cast<std::size_t>(k)] = im.rank(v);
            }
            out.push_back(c);
        }
    }
    if (div >= 2) {   // octahedra -> 4 tetrahedra each
        const index_map<3> base(div - 2);
        const int pi[6] = { 0, 2, 0, 1, 1, 0 };   // p01, p23, p02, p12, p13, p03
        const int pj[6] = { 1, 3, 2, 2, 3, 3 };
        int q[6];
        for (int r = 0; r < base.size(); ++r) {
            const multi_index<3> al = base.unrank(r);
            for (int k = 0; k < 6; ++k) {
                multi_index<3> v = al;
                v.a[static_cast<std::size_t>(pi[k])] += 1;
                v.a[static_cast<std::size_t>(pj[k])] += 1;
                q[k] = im.rank(v);
            }
            // equator cycle p02 - p12 - p13 - p03 around diagonal p01 - p23
            const int eq[4] = { q[2], q[3], q[4], q[5] };
            for (int k = 0; k < 4; ++k) {
                c[0] = q[0]; c[1] = q[1];
                c[2] = eq[k]; c[3] = eq[(k + 1) % 4];
                out.push_back(c);
            }
        }
    }
    if (div >= 3) {   // inverted tetrahedra
        const index_map<3> base(div - 3);
        for (int r = 0; r < base.size(); ++r) {
            const multi_index<3> al = base.unrank(r);
            for (int k = 0; k < 4; ++k) {
                multi_index<3> v = al;
                for (int t = 0; t < 4; ++t)
                    if (t != k) v.a[static_cast<std::size_t>(t)] += 1;
                c[static_cast<std::size_t>(k)] = im.rank(v);
            }
            out.push_back(c);
        }
    }
}

// ---- shared engine --------------------------------------------------------
// Space duck type: int num_elements() const,
//                  const element_geometry<D,T>& geometry(int) const.
// Eval duck type:  void operator()(int e, const bary_point<D,T>& lam, T* val)
//                  writing ncomp values.
template <int D, typename T, typename P, class Space, class Eval>
graphics_output<D, T, P> sample_lattice(const Space& sp,
                                        const std::vector<int>& elems,
                                        int div, int ncomp, Eval& ev) {
    if (div < 1)
        throw std::invalid_argument(
            "bfem::output_uh_for_graphics: div must be >= 1");
    const int nt = sp.num_elements();
    for (std::size_t s = 0; s < elems.size(); ++s)
        if (elems[s] < 0 || elems[s] >= nt)
            throw std::invalid_argument(
                "bfem::output_uh_for_graphics: element index out of range");

    const index_map<D> im(div);
    const int npt = im.size();
    std::vector<std::array<int, D + 1> > cell1;
    lattice_cells(im, div, cell1);
    const int nc1 = static_cast<int>(cell1.size());
    const int nsel = static_cast<int>(elems.size());

    graphics_output<D, T, P> g;
    g.ncomp = ncomp;
    g.div = div;
    g.points.zeros(nsel * npt, D + ncomp);
    g.cells.zeros(nsel * nc1, D + 1);

    std::vector<T> val(static_cast<std::size_t>(ncomp));
    for (int s = 0; s < nsel; ++s) {
        const int e = elems[static_cast<std::size_t>(s)];
        const std::array<std::array<T, D>, D + 1>& v =
            sp.geometry(e).vertices();
        const int prow = s * npt;
        for (int r = 0; r < npt; ++r) {
            const multi_index<D> al = im.unrank(r);
            bary_point<D, T> lam;
            for (int k = 0; k <= D; ++k)
                lam[static_cast<std::size_t>(k)] =
                    T(al.a[static_cast<std::size_t>(k)]) / T(div);
            for (int d = 0; d < D; ++d) {
                T x(0);
                for (int k = 0; k <= D; ++k)
                    x += lam[static_cast<std::size_t>(k)]
                         * v[static_cast<std::size_t>(k)][static_cast<std::size_t>(d)];
                g.points(prow + r, d) = x;
            }
            ev(e, lam, &val[0]);
            for (int c = 0; c < ncomp; ++c)
                g.points(prow + r, D + c) = val[static_cast<std::size_t>(c)];
        }
        const int crow = s * nc1;
        for (int c = 0; c < nc1; ++c)
            for (int k = 0; k <= D; ++k)
                g.cells(crow + c, k) =
                    prow + cell1[static_cast<std::size_t>(c)]
                                [static_cast<std::size_t>(k)];
    }
    return g;
}

inline std::vector<int> all_elements(int nt) {
    std::vector<int> e(static_cast<std::size_t>(nt));
    for (int i = 0; i < nt; ++i) e[static_cast<std::size_t>(i)] = i;
    return e;
}

// ---- eval adapters --------------------------------------------------------
template <int D, typename T, typename P, class SP>
struct fe_eval {
    fe_space<D, T, P, SP>* sp;
    const fe_function<D, T, P>* u;
    void operator()(int e, const bary_point<D, T>& lam, T* val) {
        val[0] = sp->eval(*u, e, lam);
    }
};

template <typename T, typename P, class SP>
struct c1_eval {
    c1_space<2, T, P, SP>* sp;
    const c1_function<2, T, P>* u;
    void operator()(int e, const bary_point<2, T>& lam, T* val) {
        val[0] = sp->eval(*u, e, lam);
    }
};

template <int D, typename T, typename P, class SP>
struct broken_eval {
    const broken_space<D, T, P, SP>* sp;
    const broken_field<D, T, P>* u;
    void operator()(int e, const bary_point<D, T>& lam, T* val) {
        val[0] = sp->eval(*u, e, lam);
    }
};

template <int D, typename T, typename P, class SP>
struct vfe_eval {
    fe_space<D, T, P, SP>* scalar;
    std::vector<fe_function<D, T, P> >* comp;   // extracted once per call
    void operator()(int e, const bary_point<D, T>& lam, T* val) {
        for (int d = 0; d < D; ++d)
            val[d] = scalar->eval((*comp)[static_cast<std::size_t>(d)], e, lam);
    }
};

} // namespace detail

// ---------------------------------------------------------------------------
// output_uh_for_graphics: the public overload family. The (space, u, div)
// overloads sample every element; the (space, u, elems, div) overloads
// sample exactly the listed elements in list order.
// ---------------------------------------------------------------------------

// ---- fe_space<D> (P^k, scalar, ncomp = 1) ----
template <int D, typename T, typename P, class SP>
graphics_output<D, T, P>
output_uh_for_graphics(fe_space<D, T, P, SP>& Vh,
                       const fe_function<D, T, P>& u,
                       const std::vector<int>& elems, int div = 1) {
    detail::fe_eval<D, T, P, SP> ev = { &Vh, &u };
    return detail::sample_lattice<D, T, P>(Vh, elems, div, 1, ev);
}
template <int D, typename T, typename P, class SP>
graphics_output<D, T, P>
output_uh_for_graphics(fe_space<D, T, P, SP>& Vh,
                       const fe_function<D, T, P>& u, int div = 1) {
    return output_uh_for_graphics(
        Vh, u, detail::all_elements(Vh.num_elements()), div);
}

// ---- c1_space<2> (Argyris-type, scalar, ncomp = 1) ----
template <typename T, typename P, class SP>
graphics_output<2, T, P>
output_uh_for_graphics(c1_space<2, T, P, SP>& Vh,
                       const c1_function<2, T, P>& u,
                       const std::vector<int>& elems, int div = 1) {
    detail::c1_eval<T, P, SP> ev = { &Vh, &u };
    return detail::sample_lattice<2, T, P>(Vh, elems, div, 1, ev);
}
template <typename T, typename P, class SP>
graphics_output<2, T, P>
output_uh_for_graphics(c1_space<2, T, P, SP>& Vh,
                       const c1_function<2, T, P>& u, int div = 1) {
    return output_uh_for_graphics(
        Vh, u, detail::all_elements(Vh.num_elements()), div);
}

// ---- broken_space<D> (discontinuous P_l, scalar, ncomp = 1) ----
template <int D, typename T, typename P, class SP>
graphics_output<D, T, P>
output_uh_for_graphics(const broken_space<D, T, P, SP>& Qh,
                       const broken_field<D, T, P>& p,
                       const std::vector<int>& elems, int div = 1) {
    detail::broken_eval<D, T, P, SP> ev = { &Qh, &p };
    return detail::sample_lattice<D, T, P>(Qh, elems, div, 1, ev);
}
template <int D, typename T, typename P, class SP>
graphics_output<D, T, P>
output_uh_for_graphics(const broken_space<D, T, P, SP>& Qh,
                       const broken_field<D, T, P>& p, int div = 1) {
    return output_uh_for_graphics(
        Qh, p, detail::all_elements(Qh.num_elements()), div);
}

// ---- vfe_space<D> (SV velocity, vector, ncomp = D) ----
// Components are copy-extracted once per call (vfe_space::component), then
// evaluated through the wrapped scalar space; shares the referee's thread
// contract (B-4).
template <int D, typename T, typename P, class SP>
graphics_output<D, T, P>
output_uh_for_graphics(vfe_space<D, T, P, SP>& Vh,
                       const vfe_function<D, T, P>& u,
                       const std::vector<int>& elems, int div = 1) {
    std::vector<fe_function<D, T, P> > comp;
    comp.reserve(static_cast<std::size_t>(D));
    for (int d = 0; d < D; ++d) comp.push_back(Vh.component(u, d));
    detail::vfe_eval<D, T, P, SP> ev = { &Vh.scalar(), &comp };
    return detail::sample_lattice<D, T, P>(Vh, elems, div, D, ev);
}
template <int D, typename T, typename P, class SP>
graphics_output<D, T, P>
output_uh_for_graphics(vfe_space<D, T, P, SP>& Vh,
                       const vfe_function<D, T, P>& u, int div = 1) {
    return output_uh_for_graphics(
        Vh, u, detail::all_elements(Vh.num_elements()), div);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_GRAPHICS_HPP
