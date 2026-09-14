// vcp/bfem/constants/detail/stokes_hypercircle_impl.hpp
//
// CONST-H: the hypercircle quantity kappa_h^S of the Stokes projection
// (Liu, Nakao, Oishi, Commun. Nonlinear Sci. Numer. Simul. 108 (2022)
// 106223, section 3-1) as a verified generalized eigenvalue -- the Stokes
// twin of detail::hypercircle_kappa_core of poisson_constants.hpp.
//
// Primary source, verbatim (CONST-H P0 correspondence table):
//     kappa_h = max_{f_h in X_h} min_{p_h in RT_h, v_h in V_h}
//               ||p_h - grad v_h|| / ||f_h||,
//     "where the minimization respect to p_h is subject to the condition
//      div p_h + grad phi_h + f_h = 0 for certain phi_h in U_h",
//     "The computation of kappa_h requires solving a matrix eigenvalue
//      problem."
// with the space selection d = m = k - 1 of the paper: V_h the
// Scott-Vogelius velocity space (vector P^k on the barycentric refinement
// T_SV of the mesh T_h, exactly divergence-free, homogeneous Dirichlet),
// RT_h = (RT_m)^D on T_h, X_h = (broken P^d)^D on T_h, U_h = P^k on T_h
// (phi_h modulo constants).  The two meshes are INDEPENDENT (paper section 5-1,
// "independent mesh selection ... NConf ~ NMix"); T_SV = alfeld_refine(T_h).
//
// ---------------------------------------------------------------------------
// DERIVATION RECORD [self-derived, based on Lemma 6 / 自己導出・Lemma 6 に基づく]
// (owner ruling of gate 1, item (ii); the paper's [13] was not at hand).
//
// Three premises, each checked against the paper's text:
//   (P1) V_h subset V:  v_h in V_h is exactly divergence-free and vanishes
//        on the boundary ("Vh (subset V) ... strictly divergence-free over
//        the domain", section 3; the SV constraint (div v_h, q_h) = 0 for all
//        q_h in X_{h,0} together with int div v_h = 0 gives div v_h == 0 in
//        the broken space that contains it).
//   (P2) phi_h in U_h subset H^1(Omega): the multiplier of the constraint is
//        a conforming P^k function on T_h, so grad phi_h is a legitimate
//        H(div)-free term of Lemma 6's "div p + grad phi + f = 0, phi in H^1".
//   (P3) the same f_h in both minimizations: the v_h-problem and the
//        p_h-problem are coupled ONLY through f_h (paper: "for each fh in
//        Xh, the minimization of ||ph - grad vh|| is to find uh and ph that
//        minimize ||grad u - uh|| and ||grad u - ph||, respectively").
//
// The three lines.  Fix f_h.  For any admissible p_h (constraint) and any
// v_h in V_h, integrating by parts (p_h in H(div), v_h in H^1_0) and using
// the constraint, then (P1):
//   (1) (p_h, grad v_h) = -(div p_h, v_h) = (grad phi_h + f_h, v_h)
//                       = -(phi_h, div v_h) + (f_h, v_h) = (f_h, v_h),
// so ||p_h - grad v_h||^2 = ||p_h||^2 - 2 (f_h, v_h) + ||grad v_h||^2 and the
// two minimizations DECOUPLE:
//   (2) min over v_h of ||grad v_h||^2 - 2 (f_h, v_h) is attained by the
//       Galerkin solution v_h* of (grad v_h, grad w_h) = (f_h, w_h) for all
//       w_h in V_h (= P_h applied to the exact Stokes solution for f_h, the
//       P_h minimality), with value -||grad v_h*||^2;
//       min over (p_h, phi_h) of ||p_h||^2 subject to the constraint is the
//       equality constrained least squares problem (p_h*, phi_h*).
//   (3) hence min ||p_h - grad v_h||^2 = ||p_h*||^2 - ||grad v_h*||^2, and
//       kappa_h^2 = max_{f_h} (||p_h*||^2 - ||grad v_h*||^2) / ||f_h||^2 is
//       the largest generalized eigenvalue of (Q, M3) with
//           Q  = H^T P3 H - K^T S K,   K^T S K = Bmix^T K,
//       K the solution operator of the SV saddle system, H that of the
//       constrained least squares KKT system (both f_h -> coefficients).
//       Lemma 6 with v = v_h*, p = p_h* (admissible by (P2), Lemma 6 needs
//       phi in H^1 only) then gives ||grad(u - v_h*)|| <= ||p_h* - grad
//       v_h*|| <= kappa_h ||f_h||, which is (15) of the paper for the
//       projection P_h (P_h minimality: ||grad(u - P_h u)|| <= ||grad(u - v_h)||
//       for every v_h in V_h).
//
// What is verified where.  (P1) is enforced by the saddle system (the
// pressure test space is the whole X_{h,0} after the SV rows and one pinned
// dof; the CONST-H tests check ||div v_h||^2 of the Galerkin columns and the
// exact rank identity on the rational scalar).  The cross term identity (1)
// is NOT used by the tests' direct check G-T3a, which recomputes
// ||p_h - grad v_h||^2 from the fields with (p_h, grad v_h) = -(div p_h, v_h)
// only (integration by parts, RT H(div) conformity), so an error in (1)-(3)
// would surface there.
// ---------------------------------------------------------------------------
//
// MATRIX FORM (this file; component c = 0 .. D-1 is component-major as in
// vfe_space, the X_h side likewise c * nb + j):
//   K step  [ Sr  B^T ] [ v ]   [ Bmix f ]      NConf = nU + nP
//           [ B    0  ] [ p ] = [   0    ]
//   H step  [ P3   0    N3^T ] [ p   ]   [  0   ]   NMix = D nr + nphi_r + D nb
//           [ 0    0    G^T  ] [ phi ] = [  0   ]
//           [ N3   G    0    ] [ mu  ]   [ -M3  ]
//   Q = H^T P3 H - Bmix^T K, symmetric intersection (compsym), then
//   eigsymge(Q, M3); kappa_h^2 = the largest upper end.
// Sr: no-slip reduced vector stiffness on T_SV; B: (div v, q) on T_SV with
// the SV constraint rows, one pinned pressure dof and the no-slip columns
// removed; Bmix = (w_i, f_j) between the reduced velocity dofs on T_SV and
// X_h on T_h; P3 = blockdiag(RT mass); N3 = blockdiag((div p, q)); G row
// block c = (d_c phi, q) (column block c of assemble_div_velocity on T_h
// with the pinned phi column dropped); M3 = blockdiag(X_h mass).
//
// CROSS-MESH COUPLING.  The only quantity that mixes T_h and T_SV is Bmix.
// It is assembled as Bmix_SV * I with Bmix_SV the same-mesh mixed mass on
// T_SV and I the EXACT prolongation X_h(T_h) -> X_h(T_SV): on child c of the
// barycentric split the parent barycentric coordinates of a child point are
// lambda_j = mu_j + mu_c / (D+1) (j != c), lambda_c = mu_c / (D+1), so the
// child Bernstein coefficients of a parent Bernstein basis function follow
// from an exact rational Bernstein-Lagrange solve on the child domain
// points (D + 1 blocks, shared by every parent; rational_la::solve_exact),
// converted once to T (convert_traits, enclose-once).  The identity
// I^T M_SV I == M_h is checked entry by entry (interval overlap) and the
// counts go to the report; a violation is a vcp::verification_error.
//
// CERTIFICATION.  lss certifies both saddle systems (||I - RA|| < 1: this
// also discharges the full row rank of B, i.e. the non-singularity of the
// SV pair on the mesh at hand, and of the KKT matrix), H^T P3 H is crossed
// with M3 L (the Lagrange identity of the KKT system; interval overlap of
// every entry), Q is intersected with its transpose (compsym), and
// eigsymge encloses every generalized eigenvalue.  Nothing of kappa is a
// literature constant: the paper's printed value is a TEST-SIDE literature
// value only (gate G-T3b, direction kappa(ours) <= paper).
//
// Scope: the dimension generic core below is exercised and gated for D = 3
// only (owner ruling: 3D direct; the SV verified construction class used is
// "Alfeld 3D, k >= 3").  The public entry points of stokes_constants.hpp
// reject D != 3 at compile time.
//
// Lexical regime: no decimal literal, no `double` / `float` token (gate
// G-T6, the regime of poisson_constants.hpp).

#ifndef VCP_BFEM_CONSTANTS_DETAIL_STOKES_HYPERCIRCLE_IMPL_HPP
#define VCP_BFEM_CONSTANTS_DETAIL_STOKES_HYPERCIRCLE_IMPL_HPP

#include <vector>
#include <array>
#include <string>
#include <cstddef>
#include <stdexcept>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
#include <vcp/spmatrix.hpp>
#include <vcp/imats.hpp>
#include <vcp/spimats.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/detail/rational_la.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/rt/rt_space.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/rt/rt_assemble.hpp>
#include <vcp/bfem/sv/alfeld.hpp>
#include <vcp/bfem/sv/vfe_space.hpp>
#include <vcp/bfem/sv/sv_assemble.hpp>
#include <vcp/bfem/sv/sv_rows.hpp>
#include <vcp/bfem/sv/linear_reduction.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// stokes_hypercircle_report<T>: sizes, certificates and the spectrum of one
// call of the core.  kappa2 is the point interval of the largest upper end
// (the convention of ritz_projection_constant_set_h01::kappa2); lambda keeps
// every enclosure so a caller can see the width.
// ---------------------------------------------------------------------------
template <typename T>
struct stokes_hypercircle_report {
    int degree_k;
    int n_coarse_elements;        // T_h
    int n_sv_elements;            // T_SV = (D + 1) n_coarse_elements
    int n_rt_scalar;              // nr: one RT_m component on T_h
    int n_x_scalar;               // nb: one broken P^d component on T_h
    int n_phi;                    // U_h = P^k on T_h
    int n_phi_reduced;            // n_phi - 1
    int n_velocity;               // nU: no-slip reduced vector P^k on T_SV
    int n_pressure;               // nP: X_{h,0} after SV rows and one pin
    int dim_divfree;              // nU - nP = dim V_h
    int n_conf;                   // nU + nP
    int n_mix;                    // D nr + n_phi_reduced + D nb
    int sv_constraint_rows;
    int sv_singular_entities;     // singular edges (D = 3) / vertices (D = 2)
    std::string sv_hint;
    long long prolongation_pairs_checked;   // entries of I^T M_SV I vs M_h
    long long prolongation_violations;
    long long lagrange_pairs_checked;       // entries of H^T P3 H vs M3 L
    long long lagrange_violations;
    long long sym_pairs_checked;            // compsym (i < j)
    std::vector<T> lambda;                  // eigsymge enclosures, ascending
    T kappa2;                               // point interval, largest upper end
    T min_lambda_lower;                     // information: Q is PSD in exact arithmetic

    stokes_hypercircle_report()
        : degree_k(0), n_coarse_elements(0), n_sv_elements(0), n_rt_scalar(0),
          n_x_scalar(0), n_phi(0), n_phi_reduced(0), n_velocity(0),
          n_pressure(0), dim_divfree(0), n_conf(0), n_mix(0),
          sv_constraint_rows(0), sv_singular_entities(0), sv_hint(),
          prolongation_pairs_checked(0), prolongation_violations(0),
          lagrange_pairs_checked(0), lagrange_violations(0),
          sym_pairs_checked(0), lambda(), kappa2(), min_lambda_lower() {}
};

namespace detail {

// ---------------------------------------------------------------------------
// stokes_hypercircle_materials<T, DP>: the dense reduced matrices of the two
// saddle systems and of Q, kept only when a caller asks for them (the
// CONST-H direct check G-T3a rebuilds the constrained minimization for
// concrete f_h from these).  Everything is component-major.
// ---------------------------------------------------------------------------
template <typename T, class DP>
struct stokes_hypercircle_materials {
    vcp::matrix<T, DP> Sr;        // nU x nU
    vcp::matrix<T, DP> B;         // nP x nU
    vcp::matrix<T, DP> Bmix;      // nU x (D nb)
    vcp::matrix<T, DP> P3;        // (D nr) x (D nr)
    vcp::matrix<T, DP> N3;        // (D nb) x (D nr)
    vcp::matrix<T, DP> G;         // (D nb) x nphi_r
    vcp::matrix<T, DP> M3;        // (D nb) x (D nb)
    vcp::matrix<T, DP> Q;         // (D nb) x (D nb), after compsym
    std::vector<int> velocity_to_full;   // reduced velocity dof -> full dof (D * N_SV)
};

// ---------------------------------------------------------------------------
// exact prolongation blocks of the barycentric split, degree d: block c
// (child c = parent with vertex c replaced by the barycenter, the child
// order of alfeld_refine) maps parent Bernstein coefficients to child
// Bernstein coefficients, C(beta, alpha) = coefficient of child basis beta
// in the restriction of parent basis alpha.
// ---------------------------------------------------------------------------
template <int D>
std::vector< ::vcp::bfem::detail::rmat > stokes_prolongation_blocks(int d) {
    typedef ::vcp::bfem::detail::rational rat;
    typedef ::vcp::bfem::detail::rmat rmat;
    typedef bpoly<D, rat> rpoly;
    const index_map<D>& im = coeff_registry<D>::indices(d);
    const int n = im.size();
    std::vector<rmat> out;
    out.reserve(static_cast<std::size_t>(D + 1));
    for (int c = 0; c <= D; ++c) {
        rmat V(n, n), F(n, n);
        for (int r = 0; r < n; ++r) {
            const multi_index<D> beta = im.unrank(r);
            bary_point<D, rat> mu, lam;
            for (int j = 0; j <= D; ++j)
                mu[static_cast<std::size_t>(j)] = rat(beta.a[static_cast<std::size_t>(j)], d);
            const rat q = mu[static_cast<std::size_t>(c)] / rat(D + 1);
            for (int j = 0; j <= D; ++j)
                lam[static_cast<std::size_t>(j)] = (j == c) ? q : mu[static_cast<std::size_t>(j)] + q;
            for (int s = 0; s < n; ++s) {
                std::vector<rat> e(static_cast<std::size_t>(n), rat(0));
                e[static_cast<std::size_t>(s)] = rat(1);
                const rpoly bs = rpoly::from_coeffs(d, e);
                V.at(r, s) = eval(bs, mu);        // child basis s at child point r
                F.at(r, s) = eval(bs, lam);       // parent basis s at the same point
            }
        }
        out.push_back(::vcp::bfem::detail::solve_exact(V, F));   // V C = F
    }
    return out;
}

template <typename T, class DP>
vcp::matrix<T, DP> stokes_rmat_to_dense(const ::vcp::bfem::detail::rmat& R) {
    vcp::matrix<T, DP> M;
    M.zeros(R.rows, R.cols);
    for (int i = 0; i < R.rows; ++i)
        for (int j = 0; j < R.cols; ++j) {
            const ::vcp::bfem::detail::rational& q = R.at(i, j);
            if (!q.is_zero())
                M(i, j) = convert_traits<T>::from_rational(q.num(), q.den());
        }
    return M;
}

template <typename T>
bool stokes_intervals_overlap(const T& a, const T& b) {
    return a.lower() <= b.upper() && b.lower() <= a.upper();
}

// ---------------------------------------------------------------------------
// stokes_hypercircle_kappa_core<D, T, DP, SP>(Th, k, rep, mat): the pipeline
// of the file head.  Th is the COARSE mesh T_h; the SV mesh is built here.
// rep is always filled; mat only when non-null.
//
// Exceptions: vcp::invalid_argument (k < 1, empty mesh), vcp::verification_error
// (prolongation identity, Lagrange identity, compsym, lss or eigsymge could
// not certify), std::invalid_argument / vcp::bfem exceptions from the space
// layer (degenerate element, undecidable SV singularity on interval input).
// ---------------------------------------------------------------------------
template <int D, typename T, class DP, class SP>
void stokes_hypercircle_kappa_core(const vcp::bfem::mesh<D, T>& Th, int k,
                                   stokes_hypercircle_report<T>& rep,
                                   stokes_hypercircle_materials<T, DP>* mat) {
    typedef typename T::base_type B;
    typedef vcp::matrix<T, DP> dmat;
    typedef vcp::spmatrix<T, SP> smat;
    if (k < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::stokes_hypercircle_kappa_core: k must be >= 1 (got ", k, ")");
    if (Th.num_elements() < 1)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::stokes_hypercircle_kappa_core: empty mesh");
    const int d = k - 1, m = k - 1;
    rep = stokes_hypercircle_report<T>();
    rep.degree_k = k;

    // ---- T_SV ----
    const vcp::bfem::mesh<D, T> Tsv = alfeld_refine(Th);
    rep.n_coarse_elements = Th.num_elements();
    rep.n_sv_elements = Tsv.num_elements();

    // ---- coarse spaces: RT_m, X_h, U_h on T_h ----
    fe_space<D, T, DP, SP> fs_c(Th, k);
    vfe_space<D, T, DP, SP> vfe_c(Th, fs_c);
    broken_space<D, T, DP, SP> bs_c(Th, d);
    rt_space<D, T, DP, SP> rs_c(Th, m);
    const int nr = rs_c.ndof(), nb = bs_c.ndof(), nphi = fs_c.ndof(k);
    const int nloc = bs_c.local_size();
    smat Pm = rs_c.mass();                               Pm.finalize();
    smat Mm = bs_c.mass();                               Mm.finalize();
    smat Nm = assemble_div_mass(bs_c, rs_c);             Nm.finalize();
    smat Dm = assemble_div_velocity(bs_c, vfe_c, k);     Dm.finalize();
    dmat Pd, Md, Nd, Dd;
    Pd = Pm; Md = Mm; Nd = Nm; Dd = Dm;
    std::vector<int> pin1(1, 0);
    dirichlet_reduction<T, DP, SP> pin_phi(nphi, pin1);
    const int nphir = pin_phi.reduced_size();
    rep.n_rt_scalar = nr; rep.n_x_scalar = nb; rep.n_phi = nphi; rep.n_phi_reduced = nphir;

    // ---- SV spaces on T_SV ----
    fe_space<D, T, DP, SP> fs_s(Tsv, k);
    vfe_space<D, T, DP, SP> vfe_s(Tsv, fs_s);
    broken_space<D, T, DP, SP> bs_s(Tsv, d);
    const int Ns = fs_s.ndof(k), nbs = bs_s.ndof();
    dirichlet_reduction<T, DP, SP> noslip(vfe_s.ndof(k), vfe_s.boundary_dofs(k));
    const int nU = noslip.reduced_size();
    smat Sf = assemble_vector_stiffness(vfe_s, k);       Sf.finalize();
    smat Sr = noslip.reduce(Sf);                         Sr.finalize();
    smat Bf = assemble_div_velocity(bs_s, vfe_s, k);     Bf.finalize();
    sv_pressure_constraints<D, T> svc(Tsv, k);
    const sv_stability_hint hint = svc.hint(sv_provenance_alfeld);
    linear_reduction<T, DP, SP> svrows(nbs, svc.rows());
    smat B1 = svrows.reduce_rows(Bf);                    B1.finalize();
    dirichlet_reduction<T, DP, SP> pin_p(svrows.reduced_size(), pin1);
    smat B2 = pin_p.reduce_rows(B1);                     B2.finalize();
    smat Bsp = noslip.reduce_cols(B2);                   Bsp.finalize();
    const int nP = pin_p.reduced_size();
    smat Bmix_s = assemble_mixed_mass(fs_s, k, bs_s);    Bmix_s.finalize();   // Ns x nbs
    smat Mss = bs_s.mass();                              Mss.finalize();
    dmat Srd, Bd, Bmix_sd, Ms;
    Srd = Sr; Bd = Bsp; Bmix_sd = Bmix_s; Ms = Mss;
    rep.n_velocity = nU; rep.n_pressure = nP; rep.dim_divfree = nU - nP;
    rep.n_conf = nU + nP;
    rep.sv_constraint_rows = static_cast<int>(svc.rows().size());
    rep.sv_singular_entities = (D == 3) ? svc.num_singular_edges() : svc.num_singular_vertices();
    rep.sv_hint = hint.reason;

    // ---- exact prolongation I and the identity I^T M_SV I == M_h ----
    const std::vector< ::vcp::bfem::detail::rmat > Ic = stokes_prolongation_blocks<D>(d);
    std::vector<dmat> Icd;
    for (int c = 0; c <= D; ++c) Icd.push_back(stokes_rmat_to_dense<T, DP>(Ic[static_cast<std::size_t>(c)]));
    const int nt_h = Th.num_elements();
    long long pairs = 0, bad = 0;
    dmat Bmix_c;
    Bmix_c.zeros(Ns, nb);
    for (int e = 0; e < nt_h; ++e) {
        dmat Me;
        Me.zeros(nloc, nloc);
        for (int c = 0; c <= D; ++c) {
            const int ce = (D + 1) * e + c;
            const dmat& I = Icd[static_cast<std::size_t>(c)];
            dmat Msb = Ms({ce * nloc, (ce + 1) * nloc - 1}, {ce * nloc, (ce + 1) * nloc - 1});
            Me = Me + transpose(I) * (Msb * I);
            dmat Bb = Bmix_sd({0, Ns - 1}, {ce * nloc, (ce + 1) * nloc - 1});
            dmat Bc = Bb * I;
            for (int i = 0; i < Ns; ++i)
                for (int j = 0; j < nloc; ++j) Bmix_c(i, e * nloc + j) += Bc(i, j);
        }
        for (int i = 0; i < nloc; ++i)
            for (int j = 0; j < nloc; ++j) {
                ++pairs;
                if (!stokes_intervals_overlap(Me(i, j), Md(e * nloc + i, e * nloc + j))) ++bad;
            }
    }
    rep.prolongation_pairs_checked = pairs;
    rep.prolongation_violations = bad;
    if (bad != 0)
        vcp::throw_error<vcp::verification_error>(
            "vcp::bfem::constants::stokes_hypercircle_kappa_core: the prolongation identity "
            "I^T M_SV I == M_h failed on ", bad, " of ", pairs, " entries");

    // Bmix (nU x D nb): blockdiag(Bmix_c) with the no-slip rows removed
    const int nb3 = D * nb;
    dmat Bmix;
    Bmix.zeros(nU, nb3);
    std::vector<int> v2f;
    v2f.reserve(static_cast<std::size_t>(nU));
    for (int c = 0; c < D; ++c)
        for (int i = 0; i < Ns; ++i) {
            const int r = noslip.to_reduced(c * Ns + i);
            if (r < 0) continue;
            for (int j = 0; j < nb; ++j) Bmix(r, c * nb + j) = Bmix_c(i, j);
        }
    for (int r = 0; r < nU; ++r) v2f.push_back(noslip.to_full(r));

    // ---- K step: the SV saddle system, D nb right hand sides ----
    dmat QK;
    {
        dmat A, rhs;
        A.zeros(nU + nP, nU + nP);
        A({0, nU - 1}, {0, nU - 1}) = Srd;
        A({0, nU - 1}, {nU, nU + nP - 1}) = transpose(Bd);
        A({nU, nU + nP - 1}, {0, nU - 1}) = Bd;
        rhs.zeros(nU + nP, nb3);
        rhs({0, nU - 1}, {}) = Bmix;
        dmat sol = lss(A, rhs);
        dmat K = sol({0, nU - 1}, {});
        QK = transpose(Bmix) * K;
    }

    // ---- H step: the constrained least squares KKT system ----
    const int n3r = D * nr, nk = n3r + nphir + nb3;
    dmat P3, N3, G, M3;
    P3.zeros(n3r, n3r);
    N3.zeros(nb3, n3r);
    G.zeros(nb3, nphir);
    M3.zeros(nb3, nb3);
    for (int c = 0; c < D; ++c) {
        P3({c * nr, (c + 1) * nr - 1}, {c * nr, (c + 1) * nr - 1}) = Pd;
        N3({c * nb, (c + 1) * nb - 1}, {c * nr, (c + 1) * nr - 1}) = Nd;
        M3({c * nb, (c + 1) * nb - 1}, {c * nb, (c + 1) * nb - 1}) = Md;
        for (int i = 0; i < nb; ++i)
            for (int j = 0; j < nphi; ++j) {
                const int r = pin_phi.to_reduced(j);
                if (r >= 0) G(c * nb + i, r) = Dd(i, c * nphi + j);
            }
    }
    rep.n_mix = nk;
    dmat QH;
    {
        dmat A, rhs;
        A.zeros(nk, nk);
        A({0, n3r - 1}, {0, n3r - 1}) = P3;
        A({n3r + nphir, nk - 1}, {0, n3r - 1}) = N3;
        A({0, n3r - 1}, {n3r + nphir, nk - 1}) = transpose(N3);
        A({n3r + nphir, nk - 1}, {n3r, n3r + nphir - 1}) = G;
        A({n3r, n3r + nphir - 1}, {n3r + nphir, nk - 1}) = transpose(G);
        rhs.zeros(nk, nb3);
        rhs({n3r + nphir, nk - 1}, {}) = -M3;
        dmat sol = lss(A, rhs);
        dmat H = sol({0, n3r - 1}, {});
        dmat L = sol({n3r + nphir, nk - 1}, {});
        QH = transpose(H) * (P3 * H);
        dmat QH2 = M3 * L;                       // the Lagrange identity
        long long lp = 0, lb = 0;
        for (int i = 0; i < nb3; ++i)
            for (int j = 0; j < nb3; ++j) {
                ++lp;
                if (!stokes_intervals_overlap(QH(i, j), QH2(i, j))) ++lb;
            }
        rep.lagrange_pairs_checked = lp;
        rep.lagrange_violations = lb;
        if (lb != 0)
            vcp::throw_error<vcp::verification_error>(
                "vcp::bfem::constants::stokes_hypercircle_kappa_core: the Lagrange identity "
                "H^T P3 H == M3 L failed on ", lb, " of ", lp, " entries");
    }

    // ---- Q, symmetric intersection, generalized eigenvalues ----
    dmat Q = QH - QK;
    vcp::compsym(Q);
    rep.sym_pairs_checked = static_cast<long long>(nb3) * static_cast<long long>(nb3 - 1) / 2;
    dmat E;
    eigsymge(Q, M3, E);
    dmat ev = diag(E);
    rep.lambda.resize(static_cast<std::size_t>(nb3));
    B kappa2u = ev(0, 0).upper();
    T minlo = ev(0, 0);
    for (int i = 0; i < nb3; ++i) {
        rep.lambda[static_cast<std::size_t>(i)] = ev(i, 0);
        if (kappa2u < ev(i, 0).upper()) kappa2u = ev(i, 0).upper();
        if (ev(i, 0).lower() < minlo.lower()) minlo = ev(i, 0);
    }
    if (kappa2u < B(0))
        vcp::throw_error<vcp::verification_error>(
            "vcp::bfem::constants::stokes_hypercircle_kappa_core: kappa_h^2 upper bound is negative (",
            kappa2u, "); the enclosure is not usable");
    rep.kappa2 = T(kappa2u);
    rep.min_lambda_lower = T(minlo.lower());

    if (mat != nullptr) {
        mat->Sr = Srd; mat->B = Bd; mat->Bmix = Bmix;
        mat->P3 = P3; mat->N3 = N3; mat->G = G; mat->M3 = M3; mat->Q = Q;
        mat->velocity_to_full = v2f;
    }
}

} // namespace detail
} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_DETAIL_STOKES_HYPERCIRCLE_IMPL_HPP
