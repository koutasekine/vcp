// vcp/bfem/rt/rt_backend3.hpp
// RT phase 5c (3D RT^k): the D == 3 side of the RT dimension dispatch --
// reference facet constants, the reference face trace, rt_dof_backend<3>
// (DOF-row construction backend of the RT-L0 tables), rt_dofmap3 (RT-L3
// local-to-global with parity signs) and the 3D face-moment application of
// rt_space::interpolate.
//
// Conforms to: bfem_d3c_plan.md section 1 (normative design deltas D5C-1,
// D5C-2, D5C-4 on top of RT-L0 v0.3 / RT-L3 v0.2) and the 3D common
// infrastructure (phase 5a: topology3 / s3_perm are the single authorities
// for the face convention and the permutation sigma).
//
// Normative conventions implemented here (D5C-1):
//  - local face k = ascending local vertices excluding k (5a face convention);
//    parameterization gamma(s,t) = P0 + s(P1-P0) + t(P2-P0) over the
//    canonical face vertex order (P0, P1, P2);
//  - nu_k = (P1-P0) x (P2-P0): cross product = normal x area absorption (the
//    face version of the 2D "n ds = nu dt"; no unit-normal normalization);
//  - reference nu_hat table (RD-1-3D, hardcoded below):
//      k=0: (1,1,1)  k=1: (1,0,0)  k=2: (0,-1,0)  k=3: (0,0,1)
//    OUTWARD normal x area = (-1)^k nu_hat_k -- the outwardness ALTERNATES
//    with the face number under the ascending-order convention (unlike 2D,
//    where all three edge normals are outward). Green identities must carry
//    the (-1)^k factor (checked by RA-3-3D / RC-4-3D / AXN-5C).
//  - face moments ell^F_beta(sigma) = int_{T2hat} sigma(gamma) . nu B^k_beta,
//    PLAIN parameter integrals: the L0 2D mass tables are "coefficients of
//    |T|", so the plain value carries the extra factor |T2hat| = 1/2 (D5C-6);
//  - sign/permutation composition (local -> canonical face vertex order):
//    the parameter-change Jacobian has |det| = 1, the moment basis permutes
//    beta -> sigma(beta), nu flips by parity(sigma). Hence
//      ell^local_beta = parity(sigma) * ell^canonical_{sigma(beta)}:
//    the index permutation sigma(beta) is folded into l2g and the sign array
//    is parity (D5C-4); det/orient never appear here (consumed by the RT-L2
//    f_div factor -- the same inter-layer sign division of labor as 2D);
//  - 3D RT has NO edge degrees of freedom (R3D-5): dim(k) =
//    4 (k+1)(k+2)/2 + 3 N(3,k-1) = (k+1)(k+2)(k+4)/2.
//
// No floating point appears in this header; the only scalar constants are
// exact rationals through convert_traits (enclose-once).

#ifndef VCP_BFEM_RT_RT_BACKEND3_HPP
#define VCP_BFEM_RT_RT_BACKEND3_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/typed_tables.hpp>
#include <vcp/bfem/dofmap.hpp>            // general_family_tag
#include <vcp/bfem/d3/s3_perm.hpp>
#include <vcp/bfem/d3/topology3.hpp>      // tet_local / face_canonical_beta / mesh_topology3
#include <vcp/bfem/detail/rational_la.hpp>    // rmat (generation time only)

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// rt_dof_backend<D>: the DOF-row construction backend of the RT-L0 reference
// tables (design delta D5C-3). The dimension-uniform mathematics (spanning
// set, divergence, contractions, cache) is D-generalized in rt_tables.hpp;
// ONLY the facet DOF geometry (D = 2: edges / D = 3: faces) dispatches here.
// The D == 2 specialization lives beside its primitives in rt_tables.hpp;
// this header declares the primary and defines the D == 3 side.
// ---------------------------------------------------------------------------
template <int D>
struct rt_dof_backend;

// ---------------------------------------------------------------------------
// D = 3 reference facet constants (RD-1-3D normative table; header note)
// ---------------------------------------------------------------------------
inline int rt_face_normal3(int f, int d) {
    assert(f >= 0 && f < 4 && d >= 0 && d < 3);
    static const int nu[4][3] = {
        { 1, 1, 1 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 0, 1 }
    };
    return nu[f][d];
}

// outward normal x area = (-1)^f nu_hat_f (the alternation trap; RD-1-3D)
inline int rt_face_outward3(int f) {
    assert(f >= 0 && f < 4);
    return (f % 2 == 0) ? 1 : -1;
}

// reference gradients grad lambda_i, D = 3 (closed form; the 3D row of the
// same table rt_grad_lambda hardcodes for D = 2)
inline int rt_grad_lambda3(int i, int d) {
    assert(i >= 0 && i <= 3 && d >= 0 && d < 3);
    return (i == 0) ? -1 : ((i - 1 == d) ? 1 : 0);
}

// dim RT_k (D = 3) = (k+1)(k+2)(k+4)/2 (D5C-1; no edge block, R3D-5)
inline int rt_dim3(int k) { return (k + 1) * (k + 2) * (k + 4) / 2; }

// ---------------------------------------------------------------------------
// reference face trace (identity permutation): out[rank2] = rank3 over the
// alpha_f == 0 slice of degree m, the face multi-index written in the LOCAL
// ascending vertex order of face f. Uses THE single face-convention
// implementation face_canonical_beta (S-5A-2) with the identity code 0 --
// the forward function only; no inverse re-implementation (the trace3.hpp
// discipline).
// ---------------------------------------------------------------------------
inline std::vector<int> rt_ref_face_trace(int f, int m) {
    assert(f >= 0 && f < 4 && m >= 0);
    const index_map<3>& im3 = coeff_registry<3>::indices(m);
    const index_map<2>& im2 = coeff_registry<2>::indices(m);
    std::vector<int> out(static_cast<std::size_t>(im2.size()), -1);
    for (int r = 0; r < im3.size(); ++r) {
        multi_index<3> al = im3.unrank(r);
        if (al.a[static_cast<std::size_t>(f)] != 0) continue;   // off the face
        std::array<int, 3> bc = face_canonical_beta(al, f, 0);  // identity sigma
        multi_index<2> b2;
        b2.a[0] = bc[0];
        b2.a[1] = bc[1];
        b2.a[2] = bc[2];
        out[static_cast<std::size_t>(im2.rank(b2))] = r;
    }
    return out;
}

// ---------------------------------------------------------------------------
// rt_dof_backend<3> (generation time; rmat values are exact rationals)
// ---------------------------------------------------------------------------
template <>
struct rt_dof_backend<3> {
    static const int n_facets = 4;

    static int dim(int k) { return rt_dim3(k); }

    // moments per facet: |{beta : |beta| = k}| = (k+1)(k+2)/2
    static int facet_moments(int k) { return (k + 1) * (k + 2) / 2; }

    // facet trace coefficients of a degree-m field: N(2, m)
    static int facet_coeffs(int m) { return (m + 1) * (m + 2) / 2; }

    static int facet_normal(int f, int d) { return rt_face_normal3(f, d); }

    static int grad_lambda(int i, int d) { return rt_grad_lambda3(i, d); }

    static std::vector<int> facet_trace(int f, int m) {
        return rt_ref_face_trace(f, m);
    }

    // plain parameter integrals int_{T2hat} B^k_j B^m_t ds dt =
    // |T2hat| * (2D mass "coefficient of |T|") with |T2hat| = 1/2 (D5C-6)
    static rmat facet_mass_plain(int k, int m) {
        const mass_table<2>& M = coeff_registry<2>::mass(k, m);
        const rational half(1, 2);
        rmat R(M.rows(), M.cols());
        for (int i = 0; i < M.rows(); ++i)
            for (int j = 0; j < M.cols(); ++j) {
                const rational& v = M.at(i, j);
                if (v.is_zero()) continue;
                R.at(i, j) = half * v;
            }
        return R;
    }

    // facet elevation matrix E^{k -> m} as a tall N(2,m) x N(2,k) rmat
    static rmat facet_elevation(int k, int m) {
        const elevation_table<2>& E = coeff_registry<2>::elevation(k, m);
        const int nk = coeff_registry<2>::indices(k).size();
        const int nm = coeff_registry<2>::indices(m).size();
        rmat E2(nm, nk);
        for (int src = 0; src < nk; ++src) {
            elevation_table<2>::entry_range rr = E.row(src);
            for (elevation_table<2>::entry_iterator p = rr.begin(); p != rr.end(); ++p) {
                elevation_table<2>::entry en = *p;
                E2.at(en.target_rank, src) = *en.coeff;
            }
        }
        return E2;
    }
};

// ---------------------------------------------------------------------------
// rt_dofmap3: RT^k(3D) local-to-global map (D5C-4). Block order: face block
// (face id order x in-face beta in the canonical 2D order of index_map<2>(k))
// then interior block (element order x component major x index_map<3>(k-1)).
// The index permutation sigma(beta) is folded into l2g; the sign array is
// parity(tet_face_perm). Drives the SAME general (signed) Y2 kernels as the
// 2D rt_dofmap (general_family_tag; the kernels are unchanged -- they just
// receive parity-derived +-1 instead of the 2D fwd-derived +-1).
// ---------------------------------------------------------------------------
class rt_dofmap3 {
public:
    typedef general_family_tag family_tag;

    rt_dofmap3() : k_(0), nf_(0), nt_(0), nfl_(0), nloc_(0), nint_(0), ndof_(0) {}

    static rt_dofmap3 build(const mesh_topology3& tp, int k) {
        if (k < 0)
            throw std::invalid_argument("bfem::rt_dofmap3: k < 0");
        rt_dofmap3 dm;
        dm.k_ = k;
        dm.nf_ = tp.num_faces();
        dm.nt_ = tp.nt;
        const index_map<2>& imk = coeff_registry<2>::indices(k);
        dm.nfl_ = imk.size();                                    // (k+1)(k+2)/2
        dm.nint_ = (k >= 1) ? 3 * coeff_registry<3>::indices(k - 1).size() : 0;
        dm.nloc_ = 4 * dm.nfl_ + dm.nint_;                       // == rt_dim3(k)
        dm.ndof_ = dm.nf_ * dm.nfl_ + dm.nt_ * dm.nint_;
        dm.l2g_.assign(static_cast<std::size_t>(dm.nt_)
                       * static_cast<std::size_t>(dm.nloc_), -1);
        dm.sgn_.assign(dm.l2g_.size(), static_cast<signed char>(1));
        const int int_base = dm.nf_ * dm.nfl_;
        for (int e = 0; e < dm.nt_; ++e) {
            for (int kf = 0; kf < 4; ++kf) {
                const int f = tp.tet_face[static_cast<std::size_t>(e)]
                                         [static_cast<std::size_t>(kf)];
                const int code = tp.tet_face_perm[static_cast<std::size_t>(e)]
                                                 [static_cast<std::size_t>(kf)];
                const int par = tp.tet_face_parity[static_cast<std::size_t>(e)]
                                                  [static_cast<std::size_t>(kf)];
                for (int j = 0; j < dm.nfl_; ++j) {
                    multi_index<2> b = imk.unrank(j);
                    std::array<int, 3> bl = { { b.a[0], b.a[1], b.a[2] } };
                    std::array<int, 3> bc;
                    s3_perm::apply(code, bl, bc);
                    multi_index<2> b2;
                    b2.a[0] = bc[0];
                    b2.a[1] = bc[1];
                    b2.a[2] = bc[2];
                    std::size_t idx = static_cast<std::size_t>(e)
                                      * static_cast<std::size_t>(dm.nloc_)
                                      + static_cast<std::size_t>(kf * dm.nfl_ + j);
                    dm.l2g_[idx] = f * dm.nfl_ + imk.rank(b2);
                    dm.sgn_[idx] = static_cast<signed char>(par);
                }
            }
            for (int r = 0; r < dm.nint_; ++r) {
                std::size_t idx = static_cast<std::size_t>(e)
                                  * static_cast<std::size_t>(dm.nloc_)
                                  + static_cast<std::size_t>(4 * dm.nfl_ + r);
                dm.l2g_[idx] = int_base + e * dm.nint_ + r;
            }
        }
        return dm;
    }

    int order() const { return k_; }
    int ndof() const { return ndof_; }
    int local_size() const { return nloc_; }
    int num_elements() const { return nt_; }
    int num_face_dofs() const { return nf_ * nfl_; }
    int num_faces() const { return nf_; }
    int face_moments() const { return nfl_; }

    int global_dof(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return l2g_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(r)];
    }
    int dof_sign(int e, int r) const {
        assert(e >= 0 && e < nt_ && r >= 0 && r < nloc_);
        return sgn_[static_cast<std::size_t>(e) * static_cast<std::size_t>(nloc_)
                    + static_cast<std::size_t>(r)];
    }

private:
    int k_, nf_, nt_, nfl_, nloc_, nint_, ndof_;
    std::vector<int> l2g_;
    std::vector<signed char> sgn_;
};

// ---------------------------------------------------------------------------
// 3D face-moment application of rt_space::interpolate (the D = 3 dispatch of
// the reference DOF application; the pullback and the interior moments stay
// dimension-uniform in rt_space.hpp). hat points to the 3 components of
// sigma_hat = adj(B)(sigma o F) as bpoly<3, T>; out is the ndof x 1 global
// coefficient column. Write-once assignment in element order (K6 contract).
//   local moment  = plain face parameter integral (|T2hat| = 1/2 included);
//   global column = parity * local, at index f * nfl + rank(sigma(beta)).
// ---------------------------------------------------------------------------
template <typename T, typename Poly, typename Mat>
void rt_interp_faces3(const mesh_topology3& tp, int e, int k,
                      const Poly* const* hat, Mat& out) {
    const index_map<2>& imk = coeff_registry<2>::indices(k);
    const int nfl = imk.size();
    static const T half = rational_to<T>(1, 2);          // |T2hat|, enclose-once
    for (int kf = 0; kf < 4; ++kf) {
        const int f = tp.tet_face[static_cast<std::size_t>(e)]
                                 [static_cast<std::size_t>(kf)];
        const int code = tp.tet_face_perm[static_cast<std::size_t>(e)]
                                         [static_cast<std::size_t>(kf)];
        const int par = tp.tet_face_parity[static_cast<std::size_t>(e)]
                                          [static_cast<std::size_t>(kf)];
        // per-component face trace of the hat degrees (allocation is fine:
        // interpolate is K6 test-support scope, as in 2D)
        std::vector<int> tr[3];
        for (int d = 0; d < 3; ++d) {
            if (rt_face_normal3(kf, d) == 0) continue;
            tr[d] = rt_ref_face_trace(kf, hat[d]->degree());
        }
        for (int j = 0; j < nfl; ++j) {
            T ell(0);
            for (int d = 0; d < 3; ++d) {
                const int nu = rt_face_normal3(kf, d);
                if (nu == 0) continue;
                const Poly& hd = *hat[d];
                const int deg = hd.degree();
                const typed_mass_table<2, T>& M2 = typed_registry<2, T>::mass(k, deg);
                const std::vector<int>& trd = tr[d];
                T acc(0);
                for (int t = 0; t < M2.cols(); ++t)
                    acc += M2.at(j, t)
                           * hd.coeff(trd[static_cast<std::size_t>(t)]);
                (void)deg;
                ell += (nu > 0) ? acc : -acc;            // nonzero nu is +-1
            }
            ell = half * ell;                            // plain T2hat integral
            multi_index<2> b = imk.unrank(j);
            std::array<int, 3> bl = { { b.a[0], b.a[1], b.a[2] } };
            std::array<int, 3> bc;
            s3_perm::apply(code, bl, bc);
            multi_index<2> b2;
            b2.a[0] = bc[0];
            b2.a[1] = bc[1];
            b2.a[2] = bc[2];
            const int gidx = f * nfl + imk.rank(b2);
            out(gidx, 0) = (par > 0) ? ell : -ell;       // write-once assignment
        }
    }
}

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_BACKEND3_HPP
