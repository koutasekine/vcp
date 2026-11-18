// vcp/bfem/cr1/vcr1_divfree_local.hpp
// STK-S1 (D-1): a LOCAL basis of the discrete divergence-free subspace of the
// vector Crouzeix-Raviart space, written down element by element as rational
// expressions of the vertex coordinates (two and three dimensions).
//
// Conforms to: sandbox/handoff/stokes_sparse/stokes_sparse_method.md (the
//              method note, chapters 3 and 4: the mathematics) and
//              sandbox/docs/design/STK-S1_design_v1.0.md (decision D-1: how
//              it enters VCP).
//
// ---------------------------------------------------------------------------
// (i) What this header is and is not
// ---------------------------------------------------------------------------
// The space is the one of vcr1_divfree.hpp, section (i):
//
//     V_{0,h} = { v_h in (CR)^D, vanishing boundary facet means :
//                 int_K div v_h = 0 on every element K },
//
// i.e. the kernel of the constraint matrix C (n_elements x D n_r) on the
// Dirichlet-reduced vector CR space (reduced vector dof d * n_r + r, the
// numbering of vcr1_space.hpp, section (ii)).  vcr1_divfree.hpp builds an
// EXACT basis of that kernel by an elimination over exact rational vertices;
// stokes_constants.hpp (CONST-K) builds an interval enclosure of the basis
// [-C1^{-1} C2; I] by a verified dense solve.  Both are DENSE.  This header
// builds a SPARSE basis Z whose columns have a support fixed by the mesh
// connectivity alone, without any elimination, and whose values are interval
// enclosures that hold for every mesh inside the vertex boxes (interval
// vertices are accepted; nothing is rationalized).
//
// It is NOT a constants header: no constant, no bound, no eigenvalue here.
// The reduced pencil Z^T blockdiag(K) Z and the eigenvalue call live in the
// constants layer (stokes_constants.hpp, detail::
// stokes_divfree_eigen_enclosures_sparse).
//
// ---------------------------------------------------------------------------
// (ii) The functions (method note, chapter 4)
// ---------------------------------------------------------------------------
// A CR function is given by its facet values (the value at the facet
// barycentre = the facet mean); the divergence integral over an element K is
// sum_{f in dK} |f| v(m_f) . n_{K,f}.
//
// D = 2.
//   tangential, one per interior edge e = (a, b) (ascending vertex pair): the
//     value b - a on e, zero on every other edge.  The flux through e is zero
//     from either side.
//   rotational, one per interior vertex a: on every edge e = (a, b) at a the
//     value rot(b - a) / |b - a|^2, rot = the counterclockwise quarter turn
//     (x, y) -> (-y, x).  rot(b - a) / |b - a| is a unit normal of e, so the
//     flux through e is +-1; in every triangle at a the two edges through a
//     carry +1 and -1.
// D = 3.
//   tangential, two per interior face f = (p0, p1, p2) (ascending): the
//     values p1 - p0 and p2 - p0 on f.
//   rotational, one per interior edge e = (a, b) (ascending) that is NOT a
//     tree edge: on every face f = (a, b, c) around e the value
//     2 N_f / |N_f|^2, N_f = (b - a) x (c - a).  |f| = |N_f| / 2, so the flux
//     through f is +-1; in every tetrahedron around e the two faces through e
//     carry +1 and -1.
//   tree edges: the rotational functions of all interior edges are linearly
//     dependent, one relation per interior vertex (measured by the method
//     note, chapter 7; not proved there).  The graph whose nodes are the
//     interior vertices plus ONE node standing for the whole boundary and
//     whose edges are the interior edges gets a spanning forest by a
//     union-find sweep over the edges in ascending edge id; the functions of
//     the forest edges are left out (tree-cotree).
//
// "Interior" is read off the frozen topology tables (mesh_topology2 /
// mesh_topology3) and nothing else: a facet is interior when the Dirichlet
// reduction keeps it (checked here to be the same as "not in boundary_edges
// / boundary_faces"), a vertex is interior when no boundary facet contains
// it, and an edge of a three dimensional mesh is interior when it is not a
// face_edge of any boundary face.  Every facet at an interior vertex (2D) /
// around an interior edge (3D) is therefore an interior facet.
//
// All values are polynomials in the vertex coordinates divided by one sum of
// squares (|b - a|^2 or |N_f|^2): no square root.  They are evaluated in the
// interval scalar as they stand.
//
// ---------------------------------------------------------------------------
// (iii) Why the columns are a basis (method note, section 4.3)
// ---------------------------------------------------------------------------
//   [in the kernel]  the flux cancellation above holds for every actual mesh.
//       In addition the interval product C Z is checked here to contain zero
//       for every column and every element that meets the support of the
//       column (the other entries are structurally zero).  C is the matrix of
//       stokes_divfree_build_basis, step 1: row e carries |K| (-D) d/dx_d
//       lambda_i at the reduced dof of component d and of the facet opposite
//       local vertex i.
//   [dimension]      dim V_{0,h} = D n_r - N_T + c with c the number of
//       connected components of the dual graph (elements joined through
//       interior facets; method note, chapter 3).  c == 1 is REQUIRED and the
//       number of columns built must be D n_r - N_T + 1, both as integer
//       checks; otherwise vcp::verification_error.  The count differs from
//       the formula by chi(Omega) - 1, so a domain that is not simply
//       connected is refused here (owner ruling 3 of the STK-S1 design: no
//       global function is added).
//   [independence]   NOT checked in this header.  The consumer hands
//       Z^T blockdiag(M) Z to the verified sparse eigensolver, whose entry
//       proves that interval matrix positive definite or refuses with no
//       claim; M is positive definite, so a positive definite Z^T M Z means
//       independent columns for every mesh inside the vertex boxes.
// The correctness of what a consumer certifies rests on these three checks
// only; no theorem of the literature on local divergence-free bases is used.
//
// ---------------------------------------------------------------------------
// (iv) Column order and scaling (normative)
// ---------------------------------------------------------------------------
// Columns: the tangential functions in ascending facet id (D - 1 per facet),
// then the rotational ones (D = 2: ascending interior vertex id, D = 3:
// ascending interior edge id, tree edges skipped).
// Each column stores, facet by facet of its support (ascending reduced facet
// id is NOT promised; the order is the order of construction), the D
// component values at the reduced vector dofs d * n_r + r, d = 0 .. D - 1.
// All D components of a support facet are stored even when a value is zero:
// the pattern is the support, never the floating point values.
// Scaling (argument): none -- the expressions as they stand; pow2 -- every
// column multiplied by the power of two nearest (on the logarithmic scale) to
// the reciprocal of its largest midpoint magnitude.  A column scaling by a
// positive number changes neither the span nor any of the three checks.
// The default of the argument is pow2, the value the Stokes entry passes
// (STK-S1 G-2 ruling: the better conditioned side by default; not a change of
// what is certified).
//
// Lexical policy (the constants layer's, STK-S1 decision D-6): no decimal
// literal and no floating point type keyword in the code of this file.
// No OpenMP pragma.

#ifndef VCP_BFEM_CR1_VCR1_DIVFREE_LOCAL_HPP
#define VCP_BFEM_CR1_VCR1_DIVFREE_LOCAL_HPP

#include <vector>
#include <array>
#include <utility>
#include <algorithm>
#include <cstddef>
#include <string>
#include <stdexcept>

#include <vcp/error.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dofmap.hpp>          // detail::mesh_topology2
#include <vcp/bfem/d3/topology3.hpp>    // detail::mesh_topology3
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/cr1/cr1_dofmap.hpp>  // detail::cr1_topology
#include <vcp/bfem/cr1/vcr1_space.hpp>

namespace vcp {
namespace bfem {

// column scaling of the local basis, section (iv)
enum vcr1_local_scaling {
    vcr1_local_scaling_none = 0,
    vcr1_local_scaling_pow2 = 1
};

// ---------------------------------------------------------------------------
// vcr1_divfree_local_basis<T>: the basis Z, column by column, and the counts
// the consumer reports.
// ---------------------------------------------------------------------------
template <typename T>
struct vcr1_divfree_local_basis {
    int n_elements;                 // N_T
    int n_reduced;                  // n_r (interior facets)
    int columns;                    // D n_r, the reduced vector dofs
    int rank;                       // N_T - 1 (c == 1, section (iii))
    int dim;                        // D n_r - N_T + 1 == number of columns of Z
    int n_tangential;               // tangential functions
    int n_rotational;               // rotational functions kept
    int n_tree_removed;             // rotational functions left out (D = 3)
    int n_components;               // c of the dual graph (1 on return)
    int scaling;                    // the vcr1_local_scaling used
    long long pairs_checked;        // entries of C Z verified to contain zero
    int facets_per_element;         // D + 1
    // element e, local index i -> reduced facet id or -1 (boundary facet):
    // the connectivity a consumer needs for the pattern of Z^T K Z
    std::vector<int> element_facets;
    // column -> the reduced facets of its support, in the order of cols
    std::vector<std::vector<int> > support;
    // column -> (reduced vector dof d * n_r + r, value); for the q-th support
    // facet the entries q * D + d, d = 0 .. D - 1
    std::vector<std::vector<std::pair<int, T> > > cols;
};

namespace detail {

// one support facet of a column before the final layout
template <int D, typename T>
struct vcr1_local_entry {
    int r;                          // reduced facet id
    std::array<T, D> v;             // the vector value on the facet
};

inline int vcr1_local_find(std::vector<int>& parent, int x) {
    while (parent[static_cast<std::size_t>(x)] != x) {
        parent[static_cast<std::size_t>(x)] =
            parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(x)])];
        x = parent[static_cast<std::size_t>(x)];
    }
    return x;
}

// the power of two p with (p m)^2 in [1/2, 2), m > 0; 1 when m is not a
// positive finite number or the search leaves the finite range.  Only the
// choice of the scale depends on the rounding of these products; a column
// multiplied by any positive number stays a valid column.
template <typename B>
B vcr1_local_pow2_scale(const B& m) {
    const B one(1), two(2), zero(0);
    if (!(m > zero)) return one;
    B p(1);
    const int cap = 4096;
    int guard = 0;
    while (guard < cap) {
        const B q = p * m;
        if (!(q * q >= two)) break;
        p = p / two;
        ++guard;
    }
    while (guard < cap) {
        const B q = p * m;
        if (!(q * q * two < one)) break;
        p = p * two;
        ++guard;
    }
    if (guard >= cap || !(p > zero) || !(p - p == zero)) return one;
    return p;
}

// the squared length must be certainly positive before it divides
template <typename T>
void vcr1_local_require_positive(const T& s, const char* who, const char* what, int id) {
    typedef typename T::base_type B;
    if (!(s.lower() > B(0)))
        vcp::throw_error<vcp::verification_error>(
            who, ": ", what, " ", id,
            " is not certainly nondegenerate (its squared measure is not certainly positive)");
}

template <int D> struct vcr1_local_generator;

// ---------------------------------------------------------------------------
// D = 2: tangential per interior edge, rotational per interior vertex
// ---------------------------------------------------------------------------
template <>
struct vcr1_local_generator<2> {
    template <typename T, typename P, class SP>
    static void build(const mesh<2, T>& msh, const mesh_topology2& tp,
                      const dirichlet_reduction<T, P, SP>& red, const char* who,
                      std::vector<std::vector<vcr1_local_entry<2, T> > >& out,
                      int& n_tan, int& n_rot, int& n_tree) {
        const int ne = tp.num_edges();
        std::vector<int> bvert(static_cast<std::size_t>(tp.nv), 0);
        for (std::size_t k = 0; k < tp.boundary_edges.size(); ++k) {
            const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(tp.boundary_edges[k])];
            bvert[static_cast<std::size_t>(ev[0])] = 1;
            bvert[static_cast<std::size_t>(ev[1])] = 1;
        }
        n_tan = 0;
        n_rot = 0;
        n_tree = 0;
        // tangential: b - a on the interior edge (a, b)
        for (int f = 0; f < ne; ++f) {
            const int r = red.to_reduced(f);
            if (r < 0) continue;
            const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(f)];
            const std::array<T, 2>& pa = msh.vertex(ev[0]);
            const std::array<T, 2>& pb = msh.vertex(ev[1]);
            vcr1_local_entry<2, T> en;
            en.r = r;
            en.v[0] = pb[0] - pa[0];
            en.v[1] = pb[1] - pa[1];
            out.push_back(std::vector<vcr1_local_entry<2, T> >(1, en));
            ++n_tan;
        }
        // the edges at every vertex, ascending edge id
        std::vector<std::vector<int> > vedges(static_cast<std::size_t>(tp.nv));
        for (int f = 0; f < ne; ++f) {
            const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(f)];
            vedges[static_cast<std::size_t>(ev[0])].push_back(f);
            vedges[static_cast<std::size_t>(ev[1])].push_back(f);
        }
        // rotational: rot(b - a) / |b - a|^2 on every edge (a, b) at the
        // interior vertex a
        for (int a = 0; a < tp.nv; ++a) {
            if (bvert[static_cast<std::size_t>(a)]) continue;
            const std::vector<int>& es = vedges[static_cast<std::size_t>(a)];
            if (es.empty()) continue;            // a vertex of no element
            std::vector<vcr1_local_entry<2, T> > col;
            const std::array<T, 2>& pa = msh.vertex(a);
            for (std::size_t k = 0; k < es.size(); ++k) {
                const int f = es[k];
                const int r = red.to_reduced(f);
                if (r < 0)
                    vcp::throw_error<vcp::verification_error>(
                        who, ": a boundary edge at an interior vertex (edge ", f,
                        ", vertex ", a, "); the topology tables are inconsistent");
                const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(f)];
                const int b = (ev[0] == a) ? ev[1] : ev[0];
                const std::array<T, 2>& pb = msh.vertex(b);
                const T t0 = pb[0] - pa[0];
                const T t1 = pb[1] - pa[1];
                const T len2 = t0 * t0 + t1 * t1;
                vcr1_local_require_positive(len2, who, "edge", f);
                vcr1_local_entry<2, T> en;
                en.r = r;
                en.v[0] = -t1 / len2;
                en.v[1] = t0 / len2;
                col.push_back(en);
            }
            out.push_back(col);
            ++n_rot;
        }
    }
};

// ---------------------------------------------------------------------------
// D = 3: two tangential per interior face, rotational per interior edge off
// the spanning forest
// ---------------------------------------------------------------------------
template <>
struct vcr1_local_generator<3> {
    template <typename T, typename P, class SP>
    static void build(const mesh<3, T>& msh, const mesh_topology3& tp,
                      const dirichlet_reduction<T, P, SP>& red, const char* who,
                      std::vector<std::vector<vcr1_local_entry<3, T> > >& out,
                      int& n_tan, int& n_rot, int& n_tree) {
        const int nf = tp.num_faces();
        const int ne = tp.num_edges();
        std::vector<int> bvert(static_cast<std::size_t>(tp.nv), 0);
        std::vector<int> bedge(static_cast<std::size_t>(ne), 0);
        for (std::size_t k = 0; k < tp.boundary_faces.size(); ++k) {
            const std::size_t f = static_cast<std::size_t>(tp.boundary_faces[k]);
            for (int i = 0; i < 3; ++i) {
                bvert[static_cast<std::size_t>(tp.faces[f][static_cast<std::size_t>(i)])] = 1;
                bedge[static_cast<std::size_t>(tp.face_edge[f][static_cast<std::size_t>(i)])] = 1;
            }
        }
        n_tan = 0;
        n_rot = 0;
        n_tree = 0;
        // tangential: p1 - p0 and p2 - p0 on the interior face (p0, p1, p2)
        for (int f = 0; f < nf; ++f) {
            const int r = red.to_reduced(f);
            if (r < 0) continue;
            const std::array<int, 3>& fv = tp.faces[static_cast<std::size_t>(f)];
            const std::array<T, 3>& p0 = msh.vertex(fv[0]);
            for (int j = 1; j <= 2; ++j) {
                const std::array<T, 3>& pj = msh.vertex(fv[static_cast<std::size_t>(j)]);
                vcr1_local_entry<3, T> en;
                en.r = r;
                for (int d = 0; d < 3; ++d)
                    en.v[static_cast<std::size_t>(d)] =
                        pj[static_cast<std::size_t>(d)] - p0[static_cast<std::size_t>(d)];
                out.push_back(std::vector<vcr1_local_entry<3, T> >(1, en));
                ++n_tan;
            }
        }
        // the faces around every edge with the vertex opposite the edge,
        // ascending face id.  face_edge[f] lists the edges (g0, g1), (g0, g2),
        // (g1, g2), so the vertex opposite entry k is g_{2 - k}.
        std::vector<std::vector<std::pair<int, int> > > efaces(static_cast<std::size_t>(ne));
        for (int f = 0; f < nf; ++f)
            for (int k = 0; k < 3; ++k)
                efaces[static_cast<std::size_t>(
                           tp.face_edge[static_cast<std::size_t>(f)][static_cast<std::size_t>(k)])]
                    .push_back(std::make_pair(
                        f, tp.faces[static_cast<std::size_t>(f)][static_cast<std::size_t>(2 - k)]));
        // the spanning forest: nodes = the interior vertices and the node nv
        // for the whole boundary, edges = the interior edges, ascending id
        std::vector<int> parent(static_cast<std::size_t>(tp.nv) + 1);
        for (int v = 0; v <= tp.nv; ++v) parent[static_cast<std::size_t>(v)] = v;
        std::vector<int> tree(static_cast<std::size_t>(ne), 0);
        for (int e = 0; e < ne; ++e) {
            if (bedge[static_cast<std::size_t>(e)]) continue;
            if (efaces[static_cast<std::size_t>(e)].empty()) continue;
            const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(e)];
            const int u = bvert[static_cast<std::size_t>(ev[0])] ? tp.nv : ev[0];
            const int w = bvert[static_cast<std::size_t>(ev[1])] ? tp.nv : ev[1];
            const int ru = vcr1_local_find(parent, u);
            const int rw = vcr1_local_find(parent, w);
            if (ru != rw) {
                parent[static_cast<std::size_t>(ru)] = rw;
                tree[static_cast<std::size_t>(e)] = 1;
                ++n_tree;
            }
        }
        // rotational: 2 N_f / |N_f|^2, N_f = (b - a) x (c - a), on every face
        // (a, b, c) around the interior edge (a, b) off the forest
        for (int e = 0; e < ne; ++e) {
            if (bedge[static_cast<std::size_t>(e)]) continue;
            if (tree[static_cast<std::size_t>(e)]) continue;
            const std::vector<std::pair<int, int> >& fs = efaces[static_cast<std::size_t>(e)];
            if (fs.empty()) continue;
            const std::array<int, 2>& ev = tp.edges[static_cast<std::size_t>(e)];
            const std::array<T, 3>& pa = msh.vertex(ev[0]);
            const std::array<T, 3>& pb = msh.vertex(ev[1]);
            std::array<T, 3> u;
            for (int d = 0; d < 3; ++d)
                u[static_cast<std::size_t>(d)] =
                    pb[static_cast<std::size_t>(d)] - pa[static_cast<std::size_t>(d)];
            std::vector<vcr1_local_entry<3, T> > col;
            for (std::size_t k = 0; k < fs.size(); ++k) {
                const int f = fs[k].first;
                const int r = red.to_reduced(f);
                if (r < 0)
                    vcp::throw_error<vcp::verification_error>(
                        who, ": a boundary face around an interior edge (face ", f,
                        ", edge ", e, "); the topology tables are inconsistent");
                const std::array<T, 3>& pc = msh.vertex(fs[k].second);
                std::array<T, 3> w;
                for (int d = 0; d < 3; ++d)
                    w[static_cast<std::size_t>(d)] =
                        pc[static_cast<std::size_t>(d)] - pa[static_cast<std::size_t>(d)];
                std::array<T, 3> nrm;
                nrm[0] = u[1] * w[2] - u[2] * w[1];
                nrm[1] = u[2] * w[0] - u[0] * w[2];
                nrm[2] = u[0] * w[1] - u[1] * w[0];
                const T n2 = nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2];
                vcr1_local_require_positive(n2, who, "face", f);
                vcr1_local_entry<3, T> en;
                en.r = r;
                for (int d = 0; d < 3; ++d)
                    en.v[static_cast<std::size_t>(d)] =
                        T(2) * nrm[static_cast<std::size_t>(d)] / n2;
                col.push_back(en);
            }
            out.push_back(col);
            ++n_rot;
        }
    }
};

} // namespace detail

// ---------------------------------------------------------------------------
// build_vcr1_divfree_local_basis<D, T, P, SP>(msh, V, red, scaling, who): the
// local basis of section (ii) on the mesh msh, for the vector CR space V of
// that mesh and the SCALAR Dirichlet reduction red of its boundary facets.
//
// Steps (STK-S1 design, D-1):
//   1. c = the number of connected components of the dual graph; c != 1 is
//      vcp::verification_error.
//   2. the tangential and the rotational functions, section (ii).
//   3. D = 3: the spanning forest, its rotational functions left out.
//   4. the number of columns must be D n_r - N_T + 1 (integer check),
//      otherwise vcp::verification_error with both numbers.
//   5. the interval product C Z contains zero for every column and every
//      element meeting its support, otherwise vcp::verification_error.
// The scaling of section (iv) is applied between steps 4 and 5, so step 5
// checks the columns that are returned.
//
// T: a kv interval scalar (T::base_type, lower(), upper(), mid).
// `who` prefixes every exception text (the qualified name of the entry point
// that calls this function).
// Exceptions: vcp::invalid_argument (msh / V / red do not belong together),
// vcp::verification_error (steps 1, 4, 5; an edge or a face whose squared
// measure is not certainly positive).
// ---------------------------------------------------------------------------
template <int D, typename T, typename P, class SP>
vcr1_divfree_local_basis<T> build_vcr1_divfree_local_basis(
        const mesh<D, T>& msh,
        const vcr1_space<D, T, P, SP>& V,
        const dirichlet_reduction<T, P, SP>& red,
        vcr1_local_scaling scaling = vcr1_local_scaling_pow2,
        const char* who = "vcp::bfem::build_vcr1_divfree_local_basis") {
    static_assert(D == 2 || D == 3,
                  "bfem::build_vcr1_divfree_local_basis: only D == 2 or D == 3");
    typedef typename T::base_type B;
    typedef typename detail::cr1_topology<D>::topology_type topology_type;

    const topology_type tp = topology_type::build(msh);
    const int nt = tp.nt;
    const int nfac = detail::cr1_topology<D>::num_facets(tp);
    const int nr = red.reduced_size();
    if (nt != V.num_elements() || nfac != V.scalar_ndof() || nfac != red.full_size())
        vcp::throw_error<vcp::invalid_argument>(
            std::string(who) + ": the mesh, the vector CR space and the Dirichlet "
            "reduction do not belong together (sizes)");
    {
        // the reduction must remove exactly the boundary facets
        std::vector<int> isb(static_cast<std::size_t>(nfac), 0);
        const std::vector<int>& bf = detail::cr1_topology<D>::boundary_facets(tp);
        for (std::size_t k = 0; k < bf.size(); ++k) isb[static_cast<std::size_t>(bf[k])] = 1;
        for (int f = 0; f < nfac; ++f)
            if ((red.to_reduced(f) < 0) != (isb[static_cast<std::size_t>(f)] != 0))
                vcp::throw_error<vcp::invalid_argument>(
                    std::string(who) + ": the Dirichlet reduction is not the reduction "
                    "of the boundary facets of the mesh");
    }

    vcr1_divfree_local_basis<T> out;
    out.n_elements = nt;
    out.n_reduced = nr;
    out.columns = D * nr;
    out.rank = nt - 1;
    out.dim = 0;
    out.n_tangential = 0;
    out.n_rotational = 0;
    out.n_tree_removed = 0;
    out.n_components = 0;
    out.scaling = static_cast<int>(scaling);
    out.pairs_checked = 0;
    out.facets_per_element = D + 1;

    // the connectivity: element -> reduced facets, reduced facet -> elements
    out.element_facets.assign(static_cast<std::size_t>(nt) * static_cast<std::size_t>(D + 1), -1);
    std::vector<int> facet_elem(static_cast<std::size_t>(nr) * 2, -1);
    for (int e = 0; e < nt; ++e)
        for (int i = 0; i <= D; ++i) {
            const int r = red.to_reduced(V.scalar().dofs().l2g(e, i));
            out.element_facets[static_cast<std::size_t>(e) * static_cast<std::size_t>(D + 1)
                               + static_cast<std::size_t>(i)] = r;
            if (r < 0) continue;
            const std::size_t s = static_cast<std::size_t>(r) * 2;
            if (facet_elem[s] < 0) facet_elem[s] = e;
            else facet_elem[s + 1] = e;
        }

    // 1. the dual graph
    {
        std::vector<int> parent(static_cast<std::size_t>(nt));
        for (int e = 0; e < nt; ++e) parent[static_cast<std::size_t>(e)] = e;
        int c = nt;
        for (int r = 0; r < nr; ++r) {
            const int e0 = facet_elem[static_cast<std::size_t>(r) * 2];
            const int e1 = facet_elem[static_cast<std::size_t>(r) * 2 + 1];
            if (e0 < 0 || e1 < 0)
                vcp::throw_error<vcp::verification_error>(
                    who, ": an interior facet with fewer than two elements (reduced facet ",
                    r, "); the topology tables are inconsistent");
            const int r0 = detail::vcr1_local_find(parent, e0);
            const int r1 = detail::vcr1_local_find(parent, e1);
            if (r0 != r1) {
                parent[static_cast<std::size_t>(r0)] = r1;
                --c;
            }
        }
        out.n_components = c;
        if (c != 1)
            vcp::throw_error<vcp::verification_error>(
                who, ": the dual graph of the mesh (elements joined through interior "
                "facets) has ", c, " connected components; exactly one is required");
    }

    // 2., 3. the functions
    std::vector<std::vector<detail::vcr1_local_entry<D, T> > > raw;
    detail::vcr1_local_generator<D>::build(msh, tp, red, who, raw, out.n_tangential,
                                           out.n_rotational, out.n_tree_removed);

    // 4. the count
    {
        const long long have = static_cast<long long>(raw.size());
        const long long want = static_cast<long long>(D) * nr - nt + 1;
        if (have != want)
            vcp::throw_error<vcp::verification_error>(
                who, ": the number of local basis functions (", have, " = ",
                out.n_tangential, " tangential + ", out.n_rotational,
                " rotational) differs from D n_r - N_T + 1 = ", want,
                " (dimension check; the domain is not simply connected, or the mesh is "
                "outside the scope of the local basis)");
        out.dim = static_cast<int>(have);
    }

    // scaling and the final layout
    const std::size_t ncol = raw.size();
    out.support.resize(ncol);
    out.cols.resize(ncol);
    for (std::size_t j = 0; j < ncol; ++j) {
        std::vector<detail::vcr1_local_entry<D, T> >& col = raw[j];
        if (scaling == vcr1_local_scaling_pow2) {
            B m(0);
            for (std::size_t q = 0; q < col.size(); ++q)
                for (int d = 0; d < D; ++d) {
                    B a = mid(col[q].v[static_cast<std::size_t>(d)]);
                    if (a < B(0)) a = -a;
                    if (a > m) m = a;
                }
            const T p(detail::vcr1_local_pow2_scale(m));
            for (std::size_t q = 0; q < col.size(); ++q)
                for (int d = 0; d < D; ++d)
                    col[q].v[static_cast<std::size_t>(d)] = p * col[q].v[static_cast<std::size_t>(d)];
        }
        out.support[j].reserve(col.size());
        out.cols[j].reserve(col.size() * static_cast<std::size_t>(D));
        for (std::size_t q = 0; q < col.size(); ++q) {
            out.support[j].push_back(col[q].r);
            for (int d = 0; d < D; ++d)
                out.cols[j].push_back(std::make_pair(d * nr + col[q].r,
                                                     col[q].v[static_cast<std::size_t>(d)]));
        }
    }

    // 5. C Z contains zero: every column, every element meeting its support
    {
        std::vector<int> where(static_cast<std::size_t>(nr), -1);   // reduced facet -> position in the column
        std::vector<int> elems;
        for (std::size_t j = 0; j < ncol; ++j) {
            const std::vector<int>& sup = out.support[j];
            elems.clear();
            for (std::size_t q = 0; q < sup.size(); ++q) {
                where[static_cast<std::size_t>(sup[q])] = static_cast<int>(q);
                elems.push_back(facet_elem[static_cast<std::size_t>(sup[q]) * 2]);
                elems.push_back(facet_elem[static_cast<std::size_t>(sup[q]) * 2 + 1]);
            }
            std::sort(elems.begin(), elems.end());
            elems.erase(std::unique(elems.begin(), elems.end()), elems.end());
            for (std::size_t k = 0; k < elems.size(); ++k) {
                const int e = elems[k];
                const element_geometry<D, T>& g = V.scalar().geometry(e);
                const T w = g.measure() * T(-D);
                T sum(0);
                for (int i = 0; i <= D; ++i) {
                    const int r = out.element_facets[static_cast<std::size_t>(e)
                                                     * static_cast<std::size_t>(D + 1)
                                                     + static_cast<std::size_t>(i)];
                    if (r < 0) continue;
                    const int q = where[static_cast<std::size_t>(r)];
                    if (q < 0) continue;
                    for (int d = 0; d < D; ++d)
                        sum += (w * g.grad_lambda(i, d))
                               * out.cols[j][static_cast<std::size_t>(q) * static_cast<std::size_t>(D)
                                             + static_cast<std::size_t>(d)].second;
                }
                ++out.pairs_checked;
                if (!(sum.lower() <= B(0) && B(0) <= sum.upper()))
                    vcp::throw_error<vcp::verification_error>(
                        who, ": the interval product C Z does not contain zero (element ",
                        e, ", column ", static_cast<long long>(j), ")");
            }
            for (std::size_t q = 0; q < sup.size(); ++q)
                where[static_cast<std::size_t>(sup[q])] = -1;
        }
    }
    return out;
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CR1_VCR1_DIVFREE_LOCAL_HPP
