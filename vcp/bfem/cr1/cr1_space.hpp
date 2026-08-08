// vcp/bfem/cr1/cr1_space.hpp
// CR1-1 (P1 Crouzeix-Raviart nonconforming element, D = 2, 3): the global
// layer cr1_space<D,T,P,SP> -- broken H^1 stiffness M_h(.,.), mass N(.,.),
// the facet-mean interpolation entry points and the mesh size h^2.
//
// Conforms to: CR1-1 design v1.0 (sections 2, 3.3) and
//              CR1-1 implementation directive v1.0 (phase P3).
//
// ---------------------------------------------------------------------------
// (i) The design intent that runs through the intended uses
// ---------------------------------------------------------------------------
// One space serves three uses: the 2D eigenvalue lower bound, the 3D
// eigenvalue lower bound and (later) Stokes.  What they share is the SPACE,
// not the constants: every one of them needs the same broken-gradient bilinear
// form, the same facet-mean degrees of freedom and the same Dirichlet facet
// list, but each attaches its own constant estimate on top.  This header
// therefore stops at the space.  It is deliberately D-generic (the only
// dimension dispatch of the whole cr1/ directory is cr1_topology<D> in
// cr1_dofmap.hpp), it owns no constant, and it never assumes conformity:
// V_h is NOT a subspace of H^1, which is recorded machine-readably in
// vcp/bfem/space_traits.hpp (h1_conforming == false).
//
// ---------------------------------------------------------------------------
// (ii) Rulings this header implements (design 0)
// ---------------------------------------------------------------------------
//  R1  the master degree of freedom is the facet PARAMETER MEAN (integral
//      form); the midpoint / centroid value is only an equivalent expression
//      ON P1 and appears nowhere in the code.
//  R2  space_traits is a separate, centrally placed header; no static_assert
//      is planted in any consumer here (that belongs to the constants track).
//  R3  the name fixes the degree: cr1_space, and the constructor takes NO
//      degree argument.
//  R4  max_edge_length_sq() (h^2) is the master API and is exact for every T
//      (made sound for EVERY T in CR1-1.1: the running maximum folds with the
//      max of T, which for kv::interval is the endpoint max, so the result
//      encloses the true maximum even when candidate enclosures overlap);
//      max_edge_length() applies a square root and is a convenience for
//      floating point / interval T only.  Being an ordinary non-virtual
//      member it is instantiated only when called, so an exact rational T
//      without a usable square root still compiles.
//  R5  no constant expression, no assist layer, no vector-valued CR, no local
//      divergence-free constraint: all of that is a later track.
//  R6  bfem gains this directory and nothing else -- not one existing header
//      is modified.
//
// ---------------------------------------------------------------------------
// (iii) Integral master form and its pointwise equivalent
// ---------------------------------------------------------------------------
// dof_e(v) is defined by the parameter-representation mean written out in
// vcp/bfem/cr1/cr1_element_op.hpp.  For v in P1 that number coincides with the
// midpoint value (D = 2) resp. the centroid value (D = 3) of the facet, and it
// also coincides with the arc-length / area mean int_e v ds / |e| for any v.
// The parameter form is the master because it carries no square root and stays
// inside an exact rational scalar.  The equivalences are stated here as a
// reading aid; the code implements the integral form only.
//
// ---------------------------------------------------------------------------
// (iv) Sources (design 1; verbatim-verified against the primary sources)
// ---------------------------------------------------------------------------
//  - degree of freedom int_{S_i} (Pi_h u - u) ds = 0: Liu 2015 equation (12),
//    Liu 2020 equation (15);
//  - "projection == interpolation" (Delta v_h = 0 inside an element,
//    d v_h / d n constant on a facet, facet integral zero, hence M_h
//    orthogonality): Liu 2015 equations (12)(13) and the proof of section 3.1;
//  - boundary condition as a vanishing facet integral: Liu-Nakao-Oishi 2022
//    section 3.2;
//  - h = the longest edge of the simplicial subdivision: Liu 2015 section 5,
//    theorem 3.4.
//  - the local basis phi_i = 1 - D lambda_i is [self-derived]; gate R-G2
//    verifies dof_{e_j}(phi_i) == delta_ij exactly.
//
// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------
// Element order accumulation into detail::coo_buffer, deterministic duplicate
// combination, detail::spm_adapter build -- the same vessel fe_space uses, reached
// through the FROZEN Y2 identity scatter kernels (family tag pn_family_tag;
// the CR1 local-to-global weight is identically +1 because a facet mean does
// not depend on the facet parameterization).  No hand-written sparse
// construction and no OpenMP pragma live in this header.

#ifndef VCP_BFEM_CR1_CR1_SPACE_HPP
#define VCP_BFEM_CR1_CR1_SPACE_HPP

#include <vector>
#include <array>
#include <algorithm>            // std::max (the totally ordered fold of R4)
#include <cmath>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/fe_space.hpp>        // detail::coo_buffer / detail::spm_adapter
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/cr1/cr1_dofmap.hpp>
#include <vcp/bfem/cr1/cr1_element_op.hpp>
#include <vcp/bfem/detail/scalar_traits.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P, class SP> class cr1_space;

// ---------------------------------------------------------------------------
// cr1_space<D, T, P, SP>
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class cr1_space {
    static_assert(D == 2 || D == 3, "bfem::cr1_space: only D == 2 or D == 3");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef detail::cr1_dofmap<D> dofmap_type;
    typedef typename dofmap_type::topology_type topology_type;

    // R3: no degree argument.  All element geometries are built here, so a
    // degenerate element fails fast and the per-program division count settles
    // at exactly one per element (the element_geometry contract).
    explicit cr1_space(const mesh<D, T>& msh)
        : mesh_(msh), topo_(), dm_(), geom_() {
        bfem_scalar_traits<T>::require();       // C-1 scalar contract
        topo_ = topology_type::build(msh);
        dm_ = dofmap_type::from_topology(topo_);
        geom_.reserve(static_cast<std::size_t>(topo_.nt));
        for (int e = 0; e < topo_.nt; ++e) {
            std::array<std::array<T, D>, D + 1> vv;
            for (int k = 0; k <= D; ++k)
                vv[static_cast<std::size_t>(k)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(k)]);
            geom_.push_back(element_geometry<D, T>::from_vertices(vv));
        }
    }

    // ---- observers ----
    int ndof() const { return dm_.ndof(); }
    int num_elements() const { return topo_.nt; }
    int num_vertices() const { return topo_.nv; }
    const dofmap_type& dofs() const { return dm_; }
    std::vector<int> boundary_dofs() const { return dm_.boundary_dofs(); }
    const element_geometry<D, T>& geometry(int e) const {
        if (e < 0 || e >= topo_.nt)
            throw std::invalid_argument("bfem::cr1_space::geometry: element out of range");
        return geom_[static_cast<std::size_t>(e)];
    }
    // the facet as its ascending global vertex tuple: the data of the facet
    // parameterization convention, republished so that a user-supplied AvgFn
    // can build exactly the parameterization the degree of freedom is defined
    // by (see cr1_element_op.hpp)
    const std::array<int, D>& facet_vertices(int f) const {
        return dm_.facet_vertices(f);
    }

    // ---- bilinear forms ----
    // M_h(u, v) = sum_K int_K grad u . grad v  (broken H^1)
    spmatrix_t stiffness() const {
        detail::coo_buffer<T> buf;
        local_matrix_type loc;
        reserve_(buf);
        for (int e = 0; e < topo_.nt; ++e) {    // element order
            element_op_type::local_stiffness(geom_[static_cast<std::size_t>(e)], loc);
            detail::scatter_matrix<T>(dm_, e, loc, D + 1, buf,
                                      typename dofmap_type::family_tag());
        }
        buf.combine();
        return detail::spm_adapter<T, SP>::build(dm_.ndof(), dm_.ndof(), buf);
    }

    // N(u, v) = int_Omega u v
    spmatrix_t mass() const {
        detail::coo_buffer<T> buf;
        local_matrix_type loc;
        reserve_(buf);
        for (int e = 0; e < topo_.nt; ++e) {
            element_op_type::local_mass(geom_[static_cast<std::size_t>(e)], loc);
            detail::scatter_matrix<T>(dm_, e, loc, D + 1, buf,
                                      typename dofmap_type::family_tag());
        }
        buf.combine();
        return detail::spm_adapter<T, SP>::build(dm_.ndof(), dm_.ndof(), buf);
    }

    // ---- mesh size (R4) ----
    // h^2 = max over all elements and all vertex pairs of the squared distance
    // (for a simplex that set IS the set of element edges).  Master API: no
    // square root, hence exact for an exact rational T.
    //
    // The running maximum folds with max(), NOT with a comparison, and the
    // result therefore ENCLOSES the true maximum for EVERY T.
    //
    // For an interval T the distinction is the whole point.  The comparison
    // operators of an interval type are CERTAIN comparisons, so a running
    // maximum written as "if (best < cand) best = cand;" keeps whichever
    // candidate came first whenever two candidate enclosures overlap, and the
    // kept enclosure need not contain the exact maximum at all.  kv documents
    // exactly this trap and supplies the remedy: max([a,b],[c,d]) =
    // [max(a,c), max(b,d)] (kv manual section 5.11; kv::interval's friend max).
    // The endpoint max is a valid enclosure of the maximum of the true values
    // because max is monotone in each argument [self-derived; gate I-G6
    // verifies it against an independently folded endpoint max on a mesh whose
    // two longest candidates provably overlap].
    //
    // For a totally ordered T (an exact rational, a floating point type) the
    // unqualified call resolves to std::max and the returned VALUE is exactly
    // what the previous comparison form produced -- the exact maximum.  That
    // equivalence is a contract: gate R-G8 is unchanged and must stay green.
    T max_edge_length_sq() const {
        T best(0);
        for (int e = 0; e < topo_.nt; ++e) {
            const std::array<int, D + 1>& el = mesh_.element(e);
            for (int a = 0; a <= D; ++a)
                for (int b = a + 1; b <= D; ++b) {
                    const std::array<T, D>& pa =
                        mesh_.vertex(el[static_cast<std::size_t>(a)]);
                    const std::array<T, D>& pb =
                        mesh_.vertex(el[static_cast<std::size_t>(b)]);
                    T s(0);
                    for (int d = 0; d < D; ++d) {
                        const T t = pa[static_cast<std::size_t>(d)]
                                  - pb[static_cast<std::size_t>(d)];
                        s += t * t;
                    }
                    using std::max;         // ADL: kv endpoint max for interval T
                    best = max(best, s);
                }
        }
        return best;
    }
    // convenience only (R4): instantiated exclusively when called, so a T
    // without a usable square root is fine as long as this member is unused
    T max_edge_length() const {
        using std::sqrt;
        return sqrt(max_edge_length_sq());
    }

    // ---- interpolation (R1: the degrees of freedom ARE the facet means) ----
    // The exact-integration machinery for restricting and integrating a
    // polynomial deliberately does NOT live in bfem (layer separation), so the
    // entry points take the facet means themselves.
    std::vector<T> interpolate_from_averages(const std::vector<T>& avg) const {
        if (static_cast<int>(avg.size()) != dm_.ndof())
            throw std::invalid_argument(
                "bfem::cr1_space::interpolate_from_averages: size != ndof");
        return avg;
    }
    // AvgFn: (global facet id) -> T, the parameter mean over that facet
    template <class AvgFn>
    std::vector<T> interpolate(AvgFn f) const {
        std::vector<T> out;
        out.reserve(static_cast<std::size_t>(dm_.ndof()));
        for (int g = 0; g < dm_.ndof(); ++g) out.push_back(f(g));
        return out;
    }

    // ---- local evaluation: sum_i dof_{l2g(e,i)} (1 - D lam_i) ----
    T value_local(int e, const std::vector<T>& dofvec, const T* lam) const {
        if (e < 0 || e >= topo_.nt)
            throw std::invalid_argument("bfem::cr1_space::value_local: element out of range");
        if (static_cast<int>(dofvec.size()) != dm_.ndof())
            throw std::invalid_argument("bfem::cr1_space::value_local: size != ndof");
        T r(0);
        for (int i = 0; i <= D; ++i)
            r += dofvec[static_cast<std::size_t>(dm_.l2g(e, i))]
                 * element_op_type::basis_value(i, lam);
        return r;
    }

private:
    typedef detail::cr1_element_op<D, T> element_op_type;
    typedef typename element_op_type::local_matrix_type local_matrix_type;

    void reserve_(detail::coo_buffer<T>& buf) const {
        buf.reserve(static_cast<std::size_t>(topo_.nt)
                    * static_cast<std::size_t>(D + 1)
                    * static_cast<std::size_t>(D + 1));
    }

    mesh<D, T> mesh_;
    topology_type topo_;
    dofmap_type dm_;
    std::vector<element_geometry<D, T> > geom_;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CR1_CR1_SPACE_HPP
