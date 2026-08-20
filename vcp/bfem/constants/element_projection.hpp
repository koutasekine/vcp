// vcp/bfem/constants/element_projection.hpp
//
// CONST-B1: guaranteed upper bound of the squared polynomial L^2 projection
// error constant C_d(K)^2 of one simplex element K (D = 2 triangle,
// D = 3 tetrahedron), Liu 2020 model problem 2 restricted to one element:
//
//     C_d(K) = sup_{u in H^1(K), |u|_1 > 0} ||u - P_d u||_{0,K} / |u|_{1,K}
//
// with P_d the L^2(K)-orthogonal projection onto P^d(K).  C_d(K)^2 = 1/lambda_1
// of the eigenproblem (grad u, grad v) = lambda ((I-P_d)u, (I-P_d)v).
//
// Method (design 1.2-1.4).  The constraint defining the quotient-transversal
// subspace V^h is (Q): the UNWEIGHTED sum of the boundary-facet parameter
// means vanishes -- rational coefficients, preserved by the CR interpolation
// (the quotient space argument; machine-checked by gate G-B4).  On the level-L
// uniform red refinement of K the CR1 space with (Q) eliminated gives the
// INVERTED pencil N_h x = mu A_h x (A_h positive definite after the
// reduction, N_h positive semi-definite); with an upper bound U of mu_max and
// C_h = c_D h_L from the CONST-A CR constants,
//
//     C_d(K)^2 <= U + (c_D h_L)^2                                  ... (star)
//
// No square root appears anywhere on this path (ruling R17).
//
// Exactness split (design 1.4).  Everything up to and including A_h and
// N_h = Z^T (M_mass - B^T A_p^{-1} B) Z is assembled in EXACT rational
// arithmetic: the sub-simplex integrals of monomials against the CR basis are
// evaluated by the barycentric factorial formula, and A_p^{-1} B goes through
// detail solve_exact (the ill conditioning of the monomial Gram A_p is thereby
// neutralized -- it only costs time, never accuracy).  Only after the
// assembly are the two matrices converted (enclose-once, entry by entry) to
// the interval scalar T, and the single verified computation is the
// generalized eigenvalue enclosure of the inverted pencil.
//
// Rational linear algebra dependency (owner-approved P0-3 option A): the
// shared seam vcp/bfem/detail/rational_la.hpp is consumed read-only, exactly
// as the c1/ and sv/ layers already do.  The rt/rational_la.hpp twin is a
// byte-identical leftover copy and is deliberately NOT the include path.
//
// Constants policy: this header introduces NO numeric constant of its own.
// c_D enters exclusively through cr_interpolation_constant /
// cr_projection_error_constant_sq (the CONST-A string table); every other
// number in this file is an exact small-integer combinatorial quantity.
// Lexical regime: traditional (no decimal literal, bare or in string).
//
// Authority: sandbox/docs/design/CONST-B1_design_v1.0.md and
// sandbox/docs/plans/CONST-B1_implementation_directive_v1.0.md.

#ifndef VCP_BFEM_CONSTANTS_ELEMENT_PROJECTION_HPP
#define VCP_BFEM_CONSTANTS_ELEMENT_PROJECTION_HPP

#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/constants/poisson_constants.hpp>
#include <vcp/bfem/detail/rational_la.hpp>
#include <vcp/bfem/refine.hpp>
#include <vcp/bfem/d3/bey_table.hpp>
#include <vcp/bfem/cr1/cr1_space.hpp>

namespace vcp {
namespace bfem {
namespace constants {

// ---------------------------------------------------------------------------
// element_projection_result (design 2): the (star) bound and its diagnostics.
// The GUARANTEED UPPER BOUND is cd_sq_upper.upper(); the lower ends carry no
// claim (same reading as ritz_projection_error_constant_h01, CM-1S U2).
// ---------------------------------------------------------------------------
template <int D, typename T>
struct element_projection_result {
    T cd_sq_upper;        // (star): C_d(K)^2 <= cd_sq_upper.upper()
    T mu_upper;           // U, upper bound of mu_max of the inverted pencil
    T ch_sq;              // (c_D h_L)^2 through the CONST-A entry point
    int kernel_dim;       // exact rational kernel dimension of the reduced N_h
    int level;            // refinement level L actually used
    long long n_reduced;  // degrees of freedom after the (Q) elimination

    element_projection_result()
        : cd_sq_upper(), mu_upper(), ch_sq(),
          kernel_dim(0), level(0), n_reduced(0) {}
};

namespace detail {

// The frozen exact-arithmetic layer, root-qualified once.  (The nested
// namespace constants::detail SHADOWS vcp::bfem::detail for qualified lookup,
// so every consumed name is pulled in explicitly here -- the CONST-A trap.)
typedef ::vcp::bfem::detail::bigint    ep_bigint;
typedef ::vcp::bfem::detail::rational  ep_rational;
typedef ::vcp::bfem::detail::rmat      ep_rmat;
using ::vcp::bfem::detail::solve_exact;
using ::vcp::bfem::detail::rank_exact;
using ::vcp::bfem::detail::mul;
using ::vcp::bfem::detail::transpose;

// ---------------------------------------------------------------------------
// exact rationalization of one point-interval coordinate (design 2, second
// bullet): the two endpoints must COINCIDE, and the point value must be a
// finite binary fraction of the base type (every representable floating point
// value is).  A genuine interval, or a value that fails to terminate within
// the bit cap, is rejected -- the exact pipeline has no way to carry it.
// ---------------------------------------------------------------------------
enum { ep_max_dyadic_bits = 16384 };

template <typename T>
ep_rational rational_of_point_interval(const T& x) {
    typedef typename T::base_type B;
    if (!(x.lower() == x.upper()))
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_sq_bound: vertex "
            "coordinate is a genuine interval (endpoints differ); the exact "
            "rational pipeline requires point-interval vertices");
    B v = x.lower();
    const B zero(0);
    if (v == zero) return ep_rational(0);
    const bool neg = (v < zero);
    B a = neg ? B(-v) : v;
    const B two(2);
    // highest exponent k with 2^k <= a
    B p(1);
    long long k = 0;
    while (p + p <= a) {
        p = p + p;
        ++k;
        if (k > static_cast<long long>(ep_max_dyadic_bits))
            vcp::throw_error<vcp::invalid_argument>(
                "vcp::bfem::constants::l2_projection_element_sq_bound: vertex "
                "coordinate magnitude exceeds the dyadic bit cap");
    }
    // greedy binary expansion a = sum bit_e 2^e, e = k, k-1, ...
    ep_bigint num(0);
    B r = a;
    long long e = k;
    while (true) {
        num = num + num;
        if (p <= r) {
            r = r - p;
            num = num + ep_bigint(1);
        }
        if (r == zero) break;
        p = p / two;
        --e;
        if (k - e > static_cast<long long>(ep_max_dyadic_bits))
            vcp::throw_error<vcp::invalid_argument>(
                "vcp::bfem::constants::l2_projection_element_sq_bound: vertex "
                "coordinate is not a finite binary fraction (bit cap hit); "
                "cannot rationalize exactly");
    }
    ep_bigint den(1);
    if (e < 0) {
        for (long long i = 0; i < -e; ++i) den.shl1();
    } else {
        for (long long i = 0; i < e; ++i) num.shl1();
    }
    ep_rational q(num, den);
    if (neg) q = -q;
    return q;
}

// ---------------------------------------------------------------------------
// red-refinement child tables in exact rational barycentric coordinates.
// D = 2: the F10 table of refine.hpp (4 children); D = 3: the Bey table of
// d3/bey_table.hpp with the normative FIXED diagonal m01-m23 (8 children).
// The tables are the frozen normative source; this header applies them to
// PHYSICAL vertices (child vertex = sum_k bary_k * parent_k, exact since the
// entries are in {0, 1/2, 1}) and adds nothing of its own.
// ---------------------------------------------------------------------------
template <int D>
struct red_child_table;

template <>
struct red_child_table<2> {
    enum { n_children = 4 };
    typedef std::array<std::array< ::vcp::bfem::bary_point<2, ep_rational>, 3>, 4>
        table_type;
    static table_type get() {
        return ::vcp::bfem::detail::red_children_2d<ep_rational>();
    }
};

template <>
struct red_child_table<3> {
    enum { n_children = 8 };
    typedef std::array<std::array< ::vcp::bfem::bary_point<3, ep_rational>, 4>, 8>
        table_type;
    static table_type get() {
        return ::vcp::bfem::detail::bey_children_3d<ep_rational>();
    }
};

// ---------------------------------------------------------------------------
// level-L uniform red refinement of ONE simplex, exact rational vertex lists.
// Vertex deduplication is exact: detail::rational is kept normalized
// (gcd == 1, positive denominator), so value equality is representation
// equality and the lexicographic std::map key needs nothing further.
// ---------------------------------------------------------------------------
template <int D>
struct refined_simplex_lists {
    std::vector<std::array<ep_rational, D> > vertices;
    std::vector<std::array<int, D + 1> > elements;
};

template <int D>
void refine_simplex_uniform_red(
        const std::array<std::array<ep_rational, D>, D + 1>& simplex_vertices,
        int level,
        refined_simplex_lists<D>& out) {
    if (level < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_sq_bound: negative "
            "refinement level");
    out.vertices.clear();
    out.elements.clear();
    for (int i = 0; i <= D; ++i)
        out.vertices.push_back(simplex_vertices[static_cast<std::size_t>(i)]);
    {
        std::array<int, D + 1> e0;
        for (int i = 0; i <= D; ++i) e0[static_cast<std::size_t>(i)] = i;
        out.elements.push_back(e0);
    }
    const typename red_child_table<D>::table_type kids = red_child_table<D>::get();
    for (int lv = 0; lv < level; ++lv) {
        std::vector<std::array<int, D + 1> > next_elements;
        std::map<std::array<ep_rational, D>, int> index_of;
        std::vector<std::array<ep_rational, D> > next_vertices;
        for (std::size_t e = 0; e < out.elements.size(); ++e) {
            std::array<std::array<ep_rational, D>, D + 1> pv;
            for (int i = 0; i <= D; ++i)
                pv[static_cast<std::size_t>(i)] =
                    out.vertices[static_cast<std::size_t>(
                        out.elements[e][static_cast<std::size_t>(i)])];
            for (int c = 0; c < static_cast<int>(red_child_table<D>::n_children); ++c) {
                std::array<int, D + 1> child;
                for (int i = 0; i <= D; ++i) {
                    std::array<ep_rational, D> phys;
                    for (int dcoord = 0; dcoord < D; ++dcoord) {
                        ep_rational acc(0);
                        for (int kvert = 0; kvert <= D; ++kvert)
                            acc += kids[static_cast<std::size_t>(c)]
                                       [static_cast<std::size_t>(i)]
                                       [static_cast<std::size_t>(kvert)]
                                 * pv[static_cast<std::size_t>(kvert)]
                                     [static_cast<std::size_t>(dcoord)];
                        phys[static_cast<std::size_t>(dcoord)] = acc;
                    }
                    typename std::map<std::array<ep_rational, D>, int>::iterator
                        it = index_of.find(phys);
                    int idx;
                    if (it == index_of.end()) {
                        idx = static_cast<int>(next_vertices.size());
                        next_vertices.push_back(phys);
                        index_of.insert(std::make_pair(phys, idx));
                    } else {
                        idx = it->second;
                    }
                    child[static_cast<std::size_t>(i)] = idx;
                }
                next_elements.push_back(child);
            }
        }
        out.vertices.swap(next_vertices);
        out.elements.swap(next_elements);
    }
}

// ---------------------------------------------------------------------------
// factorials for the barycentric integral formula (long long is exact far
// beyond the degrees this engine meets: the largest factorial taken is
// (2 d + D)! of the monomial Gram, and the public entry point caps the
// request so that this never overflows).
// ---------------------------------------------------------------------------
inline long long ep_factorial(int n) {
    assert(n >= 0 && n <= 20);
    long long r = 1;
    for (int i = 2; i <= n; ++i) r *= static_cast<long long>(i);
    return r;
}

// ---------------------------------------------------------------------------
// sparse polynomial in the barycentric coordinates (lambda_0 .. lambda_D) of
// ONE element, exact rational coefficients.  The exponent tuples are NOT
// required to be homogeneous; the closed integral formula below holds for
// every monomial:
//
//     int_S lambda^beta dx = |S| * D! * prod_i beta_i! / (|beta| + D)!
// ---------------------------------------------------------------------------
template <int D>
struct bary_polynomial {
    typedef std::array<int, D + 1> exponent_type;
    typedef std::map<exponent_type, ep_rational> term_map;
    term_map terms;

    void add_term(const exponent_type& e, const ep_rational& c) {
        if (c.is_zero()) return;
        typename term_map::iterator it = terms.find(e);
        if (it == terms.end()) {
            terms.insert(std::make_pair(e, c));
        } else {
            it->second += c;
            if (it->second.is_zero()) terms.erase(it);
        }
    }
};

template <int D>
bary_polynomial<D> bary_polynomial_product(const bary_polynomial<D>& a,
                                           const bary_polynomial<D>& b) {
    bary_polynomial<D> r;
    typedef typename bary_polynomial<D>::term_map::const_iterator cit;
    for (cit ia = a.terms.begin(); ia != a.terms.end(); ++ia) {
        for (cit ib = b.terms.begin(); ib != b.terms.end(); ++ib) {
            typename bary_polynomial<D>::exponent_type e;
            for (int i = 0; i <= D; ++i)
                e[static_cast<std::size_t>(i)] =
                    ia->first[static_cast<std::size_t>(i)]
                  + ib->first[static_cast<std::size_t>(i)];
            r.add_term(e, ia->second * ib->second);
        }
    }
    return r;
}

template <int D>
ep_rational integrate_bary_polynomial(const bary_polynomial<D>& p,
                                      const ep_rational& element_measure) {
    ep_rational acc(0);
    typedef typename bary_polynomial<D>::term_map::const_iterator cit;
    for (cit it = p.terms.begin(); it != p.terms.end(); ++it) {
        int total = 0;
        long long numer = ep_factorial(D);
        for (int i = 0; i <= D; ++i) {
            const int b = it->first[static_cast<std::size_t>(i)];
            total += b;
            numer *= ep_factorial(b);
        }
        acc += it->second * ep_rational(numer, ep_factorial(total + D));
    }
    return acc * element_measure;
}

// ---------------------------------------------------------------------------
// monomial exponent list of P^d in D physical variables, graded
// lexicographic, deterministic.  size == binom(d + D, D).
// ---------------------------------------------------------------------------
template <int D>
void monomial_exponents_up_to_degree(int d, std::vector<std::array<int, D> >& out);

template <>
inline void monomial_exponents_up_to_degree<2>(
        int d, std::vector<std::array<int, 2> >& out) {
    out.clear();
    for (int t = 0; t <= d; ++t)
        for (int i = 0; i <= t; ++i) {
            std::array<int, 2> a;
            a[0] = t - i;
            a[1] = i;
            out.push_back(a);
        }
}

template <>
inline void monomial_exponents_up_to_degree<3>(
        int d, std::vector<std::array<int, 3> >& out) {
    out.clear();
    for (int t = 0; t <= d; ++t)
        for (int i = 0; i <= t; ++i)
            for (int j = 0; j <= t - i; ++j) {
                std::array<int, 3> a;
                a[0] = t - i - j;
                a[1] = i;
                a[2] = j;
                out.push_back(a);
            }
}

// ---------------------------------------------------------------------------
// per-element monomial power table: the physical coordinate x_c restricted to
// the element is the barycentric linear form sum_k v_k[c] lambda_k; the table
// holds its powers 0..dmax so each monomial is a lookup-and-multiply.
// ---------------------------------------------------------------------------
template <int D>
struct element_monomial_powers {
    std::array<std::vector<bary_polynomial<D> >, D> powers;

    void build(const std::array<std::array<ep_rational, D>, D + 1>& verts,
               int dmax) {
        for (int c = 0; c < D; ++c) {
            std::vector<bary_polynomial<D> >& pw =
                powers[static_cast<std::size_t>(c)];
            pw.clear();
            bary_polynomial<D> one;
            typename bary_polynomial<D>::exponent_type z;
            for (int i = 0; i <= D; ++i) z[static_cast<std::size_t>(i)] = 0;
            one.add_term(z, ep_rational(1));
            pw.push_back(one);
            bary_polynomial<D> lin;
            for (int k = 0; k <= D; ++k) {
                typename bary_polynomial<D>::exponent_type ek = z;
                ek[static_cast<std::size_t>(k)] = 1;
                lin.add_term(ek, verts[static_cast<std::size_t>(k)]
                                      [static_cast<std::size_t>(c)]);
            }
            for (int j = 1; j <= dmax; ++j)
                pw.push_back(bary_polynomial_product<D>(
                    pw[static_cast<std::size_t>(j - 1)], lin));
        }
    }

    bary_polynomial<D> monomial(const std::array<int, D>& alpha) const {
        bary_polynomial<D> m =
            powers[0][static_cast<std::size_t>(alpha[0])];
        for (int c = 1; c < D; ++c)
            m = bary_polynomial_product<D>(
                m, powers[static_cast<std::size_t>(c)]
                         [static_cast<std::size_t>(alpha[static_cast<std::size_t>(c)])]);
        return m;
    }
};

// ---------------------------------------------------------------------------
// element_projection_core: every exact rational object of the pipeline, kept
// for the gates (G-B1 compares N_h against an independent construction; G-B4
// re-reduces the UNREDUCED pair with the measure-weighted constraint).
// ---------------------------------------------------------------------------
template <int D>
struct element_projection_core {
    refined_simplex_lists<D> mesh_lists;      // level-L rational lists
    int ndof;                                 // CR dofs of the refined mesh
    std::vector<int> boundary_dofs;           // ascending (topology order)
    int eliminated_dof;                       // the (Q)-pivot: boundary_dofs[0]
    ep_rmat stiffness_full;                   // CR broken stiffness, unreduced
    ep_rmat mass_full;                        // CR mass, unreduced
    ep_rmat monomial_gram;                    // A_p (n_p x n_p) on K itself
    ep_rmat monomial_cr_mixed;                // B (n_p x ndof)
    ep_rmat complement_full;                  // M - B^T A_p^{-1} B, unreduced
    ep_rmat stiffness_reduced;                // A_h = Z^T A Z
    ep_rmat complement_reduced;               // N_h = Z^T (M - B^T X) Z
    ep_rational ch_sq_exact;                  // (c_D h_L)^2, exact rational
    int kernel_dim;                           // n_reduced - rank_exact(N_h)
    int n_reduced;                            // ndof - 1

    element_projection_core()
        : mesh_lists(), ndof(0), boundary_dofs(), eliminated_dof(-1),
          stiffness_full(), mass_full(), monomial_gram(), monomial_cr_mixed(),
          complement_full(), stiffness_reduced(), complement_reduced(),
          ch_sq_exact(), kernel_dim(0), n_reduced(0) {}
};

inline ep_rmat rmat_difference(const ep_rmat& a, const ep_rmat& b) {
    if (a.rows != b.rows || a.cols != b.cols)
        vcp::throw_error<vcp::dimension_error>(
            "vcp::bfem::constants::element_projection: rmat size mismatch in "
            "difference");
    ep_rmat r(a.rows, a.cols);
    for (int i = 0; i < a.rows; ++i)
        for (int j = 0; j < a.cols; ++j)
            r.at(i, j) = a.at(i, j) - b.at(i, j);
    return r;
}

// Z^T M Z for the (Q) elimination matrix Z (n x (n-1)): the reduced unknowns
// are all dofs except the pivot j0 = boundary_dofs[0], in ascending order, and
// the pivot row of Z carries -1 at every OTHER boundary dof column.
inline ep_rmat constraint_elimination_matrix(int ndof,
                                             const std::vector<int>& bdofs) {
    if (ndof < 2 || bdofs.empty())
        vcp::throw_error<vcp::dimension_error>(
            "vcp::bfem::constants::element_projection: the refined mesh must "
            "carry at least two dofs and one boundary facet");
    const int j0 = bdofs[0];
    ep_rmat Z(ndof, ndof - 1);
    for (int i = 0; i < ndof; ++i) {
        if (i == j0) continue;
        Z.at(i, i < j0 ? i : i - 1) = ep_rational(1);
    }
    for (std::size_t k = 1; k < bdofs.size(); ++k) {
        const int b = bdofs[k];
        Z.at(j0, b < j0 ? b : b - 1) = ep_rational(-1);
    }
    return Z;
}

template <int D>
void build_element_projection_core(
        const std::array<std::array<ep_rational, D>, D + 1>& simplex_vertices,
        int d, int level,
        element_projection_core<D>& out) {
    if (d < 0)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_sq_bound: projection "
            "degree d must be >= 0 (got ", d, ")");
    // (2 d + D)! must stay inside long long (ep_factorial cap 20): the Gram
    // integrand degree is 2 d, the mixed integrand degree d + 1.
    if (2 * d + D > 20)
        vcp::throw_error<vcp::invalid_argument>(
            "vcp::bfem::constants::l2_projection_element_sq_bound: projection "
            "degree d too large for the exact factorial table (2 d + D must "
            "be <= 20; got d = ", d, ")");

    // ---- geometry of K itself (throws degenerate_element on a flat input,
    //      BEFORE any refinement work) ----
    const element_geometry<D, ep_rational> parent_geometry =
        element_geometry<D, ep_rational>::from_vertices(simplex_vertices);

    // ---- level-L red refinement, exact lists, CR numbering ----
    refine_simplex_uniform_red<D>(simplex_vertices, level, out.mesh_lists);
    const mesh<D, ep_rational> refined_mesh =
        mesh<D, ep_rational>::from_lists(out.mesh_lists.vertices,
                                         out.mesh_lists.elements);
    typedef ::vcp::bfem::detail::cr1_dofmap<D> dofmap_type;
    const dofmap_type dm = dofmap_type::from_mesh(refined_mesh);
    const int n = dm.ndof();
    out.ndof = n;
    out.boundary_dofs = dm.boundary_dofs();

    // ---- exact CR broken stiffness and mass (the same frozen local blocks
    //      cr1_space assembles; the sparse vessel is bypassed because the
    //      spmats solver stack does not compile on an exact rational scalar,
    //      which is precisely why the CR1 tests carry the ratx wrapper) ----
    typedef ::vcp::bfem::detail::cr1_element_op<D, ep_rational> element_op_type;
    typename element_op_type::local_matrix_type loc;
    out.stiffness_full = ep_rmat(n, n);
    out.mass_full = ep_rmat(n, n);
    std::vector<element_geometry<D, ep_rational> > geometries;
    geometries.reserve(out.mesh_lists.elements.size());
    for (std::size_t e = 0; e < out.mesh_lists.elements.size(); ++e) {
        std::array<std::array<ep_rational, D>, D + 1> vv;
        for (int i = 0; i <= D; ++i)
            vv[static_cast<std::size_t>(i)] =
                out.mesh_lists.vertices[static_cast<std::size_t>(
                    out.mesh_lists.elements[e][static_cast<std::size_t>(i)])];
        geometries.push_back(element_geometry<D, ep_rational>::from_vertices(vv));
    }
    for (std::size_t e = 0; e < out.mesh_lists.elements.size(); ++e) {
        const int ei = static_cast<int>(e);
        element_op_type::local_stiffness(geometries[e], loc);
        for (int i = 0; i <= D; ++i)
            for (int j = 0; j <= D; ++j)
                out.stiffness_full.at(dm.l2g(ei, i), dm.l2g(ei, j)) += loc(i, j);
        element_op_type::local_mass(geometries[e], loc);
        for (int i = 0; i <= D; ++i)
            for (int j = 0; j <= D; ++j)
                out.mass_full.at(dm.l2g(ei, i), dm.l2g(ei, j)) += loc(i, j);
    }

    // ---- monomial layer: A_p on K itself, B against the CR basis ----
    std::vector<std::array<int, D> > exponents;
    monomial_exponents_up_to_degree<D>(d, exponents);
    const int np = static_cast<int>(exponents.size());

    {
        element_monomial_powers<D> parent_powers;
        parent_powers.build(simplex_vertices, d);
        std::vector<bary_polynomial<D> > parent_monomials;
        parent_monomials.reserve(static_cast<std::size_t>(np));
        for (int a = 0; a < np; ++a)
            parent_monomials.push_back(
                parent_powers.monomial(exponents[static_cast<std::size_t>(a)]));
        out.monomial_gram = ep_rmat(np, np);
        for (int a = 0; a < np; ++a)
            for (int b = a; b < np; ++b) {
                const ep_rational v = integrate_bary_polynomial<D>(
                    bary_polynomial_product<D>(
                        parent_monomials[static_cast<std::size_t>(a)],
                        parent_monomials[static_cast<std::size_t>(b)]),
                    parent_geometry.measure());
                out.monomial_gram.at(a, b) = v;
                if (b != a) out.monomial_gram.at(b, a) = v;
            }
    }

    out.monomial_cr_mixed = ep_rmat(np, n);
    for (std::size_t e = 0; e < out.mesh_lists.elements.size(); ++e) {
        const int ei = static_cast<int>(e);
        std::array<std::array<ep_rational, D>, D + 1> vv;
        for (int i = 0; i <= D; ++i)
            vv[static_cast<std::size_t>(i)] =
                out.mesh_lists.vertices[static_cast<std::size_t>(
                    out.mesh_lists.elements[e][static_cast<std::size_t>(i)])];
        element_monomial_powers<D> pw;
        pw.build(vv, d);
        // phi_i = 1 - D lambda_i on this element
        for (int a = 0; a < np; ++a) {
            const bary_polynomial<D> mono =
                pw.monomial(exponents[static_cast<std::size_t>(a)]);
            for (int i = 0; i <= D; ++i) {
                bary_polynomial<D> phi;
                typename bary_polynomial<D>::exponent_type z;
                for (int t = 0; t <= D; ++t) z[static_cast<std::size_t>(t)] = 0;
                phi.add_term(z, ep_rational(1));
                typename bary_polynomial<D>::exponent_type eli = z;
                eli[static_cast<std::size_t>(i)] = 1;
                phi.add_term(eli, ep_rational(-D));
                const ep_rational v = integrate_bary_polynomial<D>(
                    bary_polynomial_product<D>(mono, phi),
                    geometries[e].measure());
                out.monomial_cr_mixed.at(a, dm.l2g(ei, i)) += v;
            }
        }
    }

    // ---- N = M - B^T (A_p^{-1} B), exact (the Gram ill conditioning is
    //      neutralized by the exact solve; design 1.4) ----
    {
        const ep_rmat X = solve_exact(out.monomial_gram, out.monomial_cr_mixed);
        out.complement_full = rmat_difference(
            out.mass_full, mul(transpose(out.monomial_cr_mixed), X));
    }

    // ---- (Q) elimination and the reduced pair ----
    const ep_rmat Z = constraint_elimination_matrix(n, out.boundary_dofs);
    out.eliminated_dof = out.boundary_dofs[0];
    out.stiffness_reduced = mul(transpose(Z), mul(out.stiffness_full, Z));
    out.complement_reduced = mul(transpose(Z), mul(out.complement_full, Z));
    out.n_reduced = n - 1;

    // ---- diagnostics: exact kernel dimension, exact (c_D h_L)^2 ----
    out.kernel_dim = out.n_reduced - rank_exact(out.complement_reduced);
    {
        // exact h_L^2 = max squared edge length over the refined elements
        // (same quantity cr1_space::max_edge_length_sq computes; rational is
        // totally ordered so the plain comparison fold is the exact maximum)
        ep_rational h2(0);
        for (std::size_t e = 0; e < out.mesh_lists.elements.size(); ++e) {
            for (int a = 0; a <= D; ++a)
                for (int b = a + 1; b <= D; ++b) {
                    ep_rational s(0);
                    for (int c = 0; c < D; ++c) {
                        const ep_rational t =
                            out.mesh_lists.vertices[static_cast<std::size_t>(
                                out.mesh_lists.elements[e][static_cast<std::size_t>(a)])]
                                [static_cast<std::size_t>(c)]
                          - out.mesh_lists.vertices[static_cast<std::size_t>(
                                out.mesh_lists.elements[e][static_cast<std::size_t>(b)])]
                                [static_cast<std::size_t>(c)];
                        s += t * t;
                    }
                    if (h2 < s) h2 = s;
                }
        }
        const ep_rational cd = constant_from_string<ep_rational>::get(
            D == 2 ? cr_interpolation_constant_2d_str
                   : cr_interpolation_constant_3d_str);
        out.ch_sq_exact = cd * cd * h2;
    }
}

// ---------------------------------------------------------------------------
// enclose-once conversion of the exact rational matrices to the interval
// scalar T (design 1.4: "assemble rational, then convert").
// ---------------------------------------------------------------------------
template <typename T>
T interval_of_exact_rational(const ep_rational& q) {
    return ::vcp::bfem::convert_traits<T>::from_rational(q.num(), q.den());
}

template <typename T, class DP>
void interval_matrix_of_rmat(const ep_rmat& r, vcp::matrix<T, DP>& out) {
    out.zeros(r.rows, r.cols);
    for (int i = 0; i < r.rows; ++i)
        for (int j = 0; j < r.cols; ++j)
            out(i, j) = interval_of_exact_rational<T>(r.at(i, j));
}

} // namespace detail

// ---------------------------------------------------------------------------
// l2_projection_element_sq_bound (design 2): vertex coordinates (D + 1
// point-interval points), projection degree d, refinement level L
// -> the (star) bound C_d(K)^2 <= U + (c_D h_L)^2 with diagnostics.
//
// T is a kv::interval-like scalar (T::base_type required; an exact rational T
// cannot carry the verified eigenvalue step and is rejected at compile time
// by the same contract the other entry points of this layer use).  DP / SP
// follow the dense/sparse policy convention of this header's older entry
// points; the design signature <D, T> stays callable through the defaults.
//
// The verified chain: mu_max enclosure of the INVERTED pencil
// N_h x = mu A_h x through eigsymge (A side N_h semi-definite is admissible;
// the certification ||Y B X - I|| < 1 sits on the A_h side, which is positive
// definite after the (Q) elimination -- the same shape as the CM-1 P7 call
// eigsymge(Q, Md, E) above).  U is the largest upper end of the returned
// enclosures; the kernel of N_h maps to mu = 0 and cannot disturb it.
//
// Exceptions: vcp::invalid_argument (d < 0 / degree cap / negative level /
// genuine interval vertex / non-dyadic coordinate), vcp::verification_error
// (eigsymge certification failed, or the mu_max enclosure is certainly
// negative although N_h is positive semi-definite -- inclusion broken),
// vcp::bfem::degenerate_element (flat input simplex; NOT in the vcp::error
// hierarchy, same caveat as ritz_projection_error_constant_h01),
// std::logic_error from solve_exact (the exact monomial Gram found singular
// -- an implementation bug by theory, never an input condition).
// ---------------------------------------------------------------------------
template <int D, typename T,
          class DP = vcp::imats<typename T::base_type>,
          class SP = vcp::spimats<typename T::base_type> >
element_projection_result<D, T>
l2_projection_element_sq_bound(const T vertices[D + 1][D], int d, int level) {
    static_assert(D == 2 || D == 3,
                  "vcp::bfem::constants::l2_projection_element_sq_bound: "
                  "D must be 2 or 3");
    detail::interval_scalar_contract<T>::require();
    typedef typename T::base_type B;

    // ---- exact rational pipeline (P1) ----
    std::array<std::array<detail::ep_rational, D>, D + 1> rational_vertices;
    for (int i = 0; i <= D; ++i)
        for (int c = 0; c < D; ++c)
            rational_vertices[static_cast<std::size_t>(i)]
                             [static_cast<std::size_t>(c)] =
                detail::rational_of_point_interval<T>(
                    vertices[i][c]);
    detail::element_projection_core<D> core;
    detail::build_element_projection_core<D>(rational_vertices, d, level, core);

    // ---- verified eigenvalue step on the inverted pencil ----
    vcp::matrix<T, DP> stiffness_i, complement_i, enclosure;
    detail::interval_matrix_of_rmat<T, DP>(core.stiffness_reduced, stiffness_i);
    detail::interval_matrix_of_rmat<T, DP>(core.complement_reduced, complement_i);
    eigsymge(complement_i, stiffness_i, enclosure);

    B mu_up;
    {
        const bool diag_column = (enclosure.columnsize() == 1);
        mu_up = enclosure(0, 0).upper();
        for (int i = 1; i < core.n_reduced; ++i) {
            const T ei = diag_column ? enclosure(i, 0) : enclosure(i, i);
            if (mu_up < ei.upper()) mu_up = ei.upper();
        }
    }
    if (mu_up < B(0))
        vcp::throw_error<vcp::verification_error>(
            "vcp::bfem::constants::l2_projection_element_sq_bound: the mu_max "
            "enclosure is certainly negative although N_h is positive "
            "semi-definite; the inclusion is not usable");

    // ---- (c_D h_L)^2 through the CONST-A entry point on the refined mesh ----
    element_projection_result<D, T> result;
    {
        std::vector<std::array<T, D> > vi;
        vi.reserve(core.mesh_lists.vertices.size());
        for (std::size_t v = 0; v < core.mesh_lists.vertices.size(); ++v) {
            std::array<T, D> p;
            for (int c = 0; c < D; ++c)
                p[static_cast<std::size_t>(c)] =
                    detail::interval_of_exact_rational<T>(
                        core.mesh_lists.vertices[v][static_cast<std::size_t>(c)]);
            vi.push_back(p);
        }
        const mesh<D, T> refined_interval_mesh =
            mesh<D, T>::from_lists(vi, core.mesh_lists.elements);
        const cr1_space<D, T, DP, SP> refined_space(refined_interval_mesh);
        result.ch_sq = cr_projection_error_constant_sq(refined_space);
    }

    // ---- (star) ----
    result.mu_upper = T(mu_up);
    result.cd_sq_upper = result.mu_upper + result.ch_sq;
    result.kernel_dim = core.kernel_dim;
    result.level = level;
    result.n_reduced = static_cast<long long>(core.n_reduced);
    return result;
}

} // namespace constants
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CONSTANTS_ELEMENT_PROJECTION_HPP
