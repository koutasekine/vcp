// vcp/bfem/fe_space.hpp
// Layer 3: fe_space<D,T,P> (H3-H8), the assembly core, plus the COO buffer
// and the spmatrix adapter.
//
// Conforms to: L3 external design v0.2 (sections 5) and
//              L3 internal design v0.2 (sections 4, 5, 6, 7).
//
// D = 3 enablement (phase 5b, D5B-1/D5B-2 of
// sandbox/docs/plans/bfem_d3b_plan.md): fe_space is generalized over D
// instead of being duplicated -- the ONLY dimension-dependent pieces are the
// topology type and the dofmap builder, both selected at compile time by
// detail::fe_space_backend<D> below. The external API, the SP lazy
// instantiation, the thread model and the determinism contract are identical
// for D = 2 and D = 3. fe_space<D>::dofs(m) is the single public
// construction entry of dofmap<D> (the 2D A-2 discipline, lifted to 3D:
// detail::dofmap_builder3 stays detail).
//
// spmatrix reconciliation (external design 11.1, resolved 2026-07-02, see
// sandbox/docs/reviews/bfem_api_reconciliation_review.md):
//  (i)  sequential (i, j, value) insertion exists: spmatrix::add(i, j, v);
//  (ii) component enumeration is public: coo_rows()/coo_columns()/coo_values()
//       in the COO stage and outer_index()/inner_index()/values() when
//       finalized -- the MAIN plan for dirichlet reduce(A) is adopted (the
//       triplet-refilter fallback is not needed).
//  spmatrix_t = vcp::spmatrix<T, SP>. finalize() is never called here (X5:
//  consumer-triggered finalize semantics).
//
// Sparse policy parameter SP (default vcp::spmats<T>): like the dense pair
// mats/imats, the default spmats<T> is the APPROXIMATE policy and by design
// does not accept interval scalars; interval (verified) use plugs a different
// policy into SP. The assembly below only touches the storage surface of the
// policy (resize/reserve/add + COO/CSR/CSC enumeration), so any policy that
// provides it works; a production interval sparse policy (the imats analog)
// does not exist in VCP yet and is a separate task.
//
// Thread model (A-3): one fe_space instance is NOT thread safe (it owns the
// element_op, the COO buffer and the L1 scratch); shared caches are safe.
// All accumulation happens in element order (X9) and duplicate combination is
// deterministic (stable sort + insertion-order addition): E6 bit
// reproducibility.

#ifndef VCP_BFEM_FE_SPACE_HPP
#define VCP_BFEM_FE_SPACE_HPP

#ifdef VCP_NOMP
#  ifndef VCP_BFEM_NOMP
#    define VCP_BFEM_NOMP
#  endif
#endif

#if defined(_OPENMP) && !defined(VCP_BFEM_NOMP)
#  define VCP_BFEM_USE_OPENMP 1
#  include <omp.h>
#else
#  define VCP_BFEM_USE_OPENMP 0
#endif

#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/poly_field.hpp>
#include <vcp/bfem/dofmap.hpp>
#include <vcp/bfem/d3/dofmap3.hpp>
#include <vcp/bfem/fe_function.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/element_op.hpp>
#include <vcp/bfem/detail/scalar_traits.hpp>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// coo_buffer: triplet accumulation and deterministic duplicate combination
// (X5, X9). stable_sort on (i, j) keeps the insertion (= element) order
// within a key; addition then runs front to back: the "element-order sum"
// survives sorting, the root of E6 (bit reproducibility, interval bounds
// included).
// ---------------------------------------------------------------------------
template <typename T>
struct coo_buffer {
    struct trip {
        int i, j;
        T v;
    };
    std::vector<trip> a;

    void clear() { a.clear(); }              // capacity kept (reuse across calls)
    void reserve(std::size_t n) { if (a.capacity() < n) a.reserve(n); }
    void append(const coo_buffer& o) {
        a.insert(a.end(), o.a.begin(), o.a.end());
    }
    void append_all(const std::vector<coo_buffer>& src, int nsrc) {
        std::size_t add = 0;
        for (int t = 0; t < nsrc; ++t)
            add += src[static_cast<std::size_t>(t)].a.size();

        const std::size_t old = a.size();
        a.resize(old + add);
        std::size_t pos = old;
        for (int t = 0; t < nsrc; ++t) {
            const std::vector<trip>& v = src[static_cast<std::size_t>(t)].a;
            std::copy(v.begin(), v.end(), a.begin() + pos);
            pos += v.size();
        }
    }
    void push(int i, int j, const T& v) {
        trip t;
        t.i = i;
        t.j = j;
        t.v = v;
        a.push_back(t);
    }
    static bool ij_less(const trip& x, const trip& y) {
        if (x.i != y.i) return x.i < y.i;
        return x.j < y.j;
    }
    void combine() {
        stable_sort_ij();
        std::size_t out = 0;
        std::size_t k = 0;
        while (k < a.size()) {
            trip acc = a[k];
            std::size_t j = k + 1;
            while (j < a.size() && a[j].i == acc.i && a[j].j == acc.j) {
                acc.v += a[j].v;             // insertion-order accumulation
                ++j;
            }
            a[out] = acc;
            ++out;
            k = j;
        }
        a.resize(out);
    }

private:
    std::vector<trip> tmp_;                  // merge buffer, reused (K5)

    // bottom-up stable merge sort with a persistent temporary buffer
    // (std::stable_sort allocates a fresh temporary on every call, which
    // would break the global zero-allocation contract in repeated assembly)
    void stable_sort_ij() {
        std::size_t n = a.size();
        if (n < 2) return;
        if (tmp_.size() < n) tmp_.resize(n);
        std::vector<trip>* src = &a;
        std::vector<trip>* dst = &tmp_;
        bool swapped = false;
        for (std::size_t width = 1; width < n; width *= 2) {
            const std::size_t span = 2 * width;
            const long long nruns = static_cast<long long>((n + span - 1) / span);
#if VCP_BFEM_USE_OPENMP
#pragma omp parallel for if(n >= static_cast<std::size_t>(4096))
#endif
            for (long long r = 0; r < nruns; ++r) {
                const std::size_t lo = static_cast<std::size_t>(r) * span;
                std::size_t mid = lo + width < n ? lo + width : n;
                std::size_t hi = lo + span < n ? lo + span : n;
                std::size_t p = lo, q = mid, o = lo;
                while (p < mid && q < hi) {
                    // stability: take from the left run on ties
                    if (ij_less((*src)[q], (*src)[p])) (*dst)[o++] = (*src)[q++];
                    else                               (*dst)[o++] = (*src)[p++];
                }
                while (p < mid) (*dst)[o++] = (*src)[p++];
                while (q < hi)  (*dst)[o++] = (*src)[q++];
            }
            std::vector<trip>* t = src;
            src = dst;
            dst = t;
            swapped = !swapped;
        }
        if (swapped) {
            a.swap(tmp_);                    // sorted data lives in tmp_
            a.resize(n);                     // tmp_ may be larger from earlier use
        }
    }
};

// ---------------------------------------------------------------------------
// spm_adapter: the single absorption point of the spmatrix API (X5;
// L3 internal design 5). build() inserts the combined triplets in (i, j)
// order and never calls finalize(). SP is the sparse policy (see the header
// note above).
// ---------------------------------------------------------------------------
template <typename T, class SP = vcp::spmats<T> >
struct spm_adapter {
    typedef vcp::spmatrix<T, SP> spmatrix_t;

    static spmatrix_t build(int rows, int cols, const coo_buffer<T>& c) {
        spmatrix_t A(rows, cols);
        A.reserve(static_cast<int>(c.a.size()));
        for (std::size_t k = 0; k < c.a.size(); ++k)
            A.add(c.a[k].i, c.a[k].j, c.a[k].v);
        return A;                            // finalize is NOT called (X5)
    }

    typedef typename spmatrix_t::index_type index_type;

    template <typename F>
    static void for_each_entry(const spmatrix_t& A, F f) {
        if (!A.is_finalized()) {
            const std::vector<index_type>& r = A.coo_rows();
            const std::vector<index_type>& c = A.coo_columns();
            const std::vector<T>& v = A.coo_values();
            for (std::size_t k = 0; k < v.size(); ++k) f(r[k], c[k], v[k]);
        } else if (A.format() == vcp::sparse_csr) {
            const std::vector<index_type>& outer = A.outer_index();
            const std::vector<index_type>& inner = A.inner_index();
            const std::vector<T>& v = A.values();
            for (index_type i = 0; i < A.rowsize(); ++i)
                for (index_type p = outer[static_cast<std::size_t>(i)];
                     p < outer[static_cast<std::size_t>(i) + 1]; ++p)
                    f(i, inner[static_cast<std::size_t>(p)], v[static_cast<std::size_t>(p)]);
        } else {                             // sparse_csc
            const std::vector<index_type>& outer = A.outer_index();
            const std::vector<index_type>& inner = A.inner_index();
            const std::vector<T>& v = A.values();
            for (index_type j = 0; j < A.columnsize(); ++j)
                for (index_type p = outer[static_cast<std::size_t>(j)];
                     p < outer[static_cast<std::size_t>(j) + 1]; ++p)
                    f(inner[static_cast<std::size_t>(p)], j, v[static_cast<std::size_t>(p)]);
        }
    }
};

// ---------------------------------------------------------------------------
// fe_space_backend<D> (D5B-1/D5B-2): the single dimension dispatch of the
// assembly layer. Everything else in fe_space is D-generic; only the
// topology type and the dofmap builder differ between D = 2 and D = 3.
// ---------------------------------------------------------------------------
template <int D>
struct fe_space_backend;

template <>
struct fe_space_backend<2> {
    typedef mesh_topology2 topology_type;
    static dofmap<2> build_dofmap(const topology_type& tp, int m) {
        return dofmap_builder::build(tp, m);
    }
};

template <>
struct fe_space_backend<3> {
    typedef mesh_topology3 topology_type;
    static dofmap<3> build_dofmap(const topology_type& tp, int m) {
        return dofmap_builder3::build(tp, m);
    }
};

} // namespace detail

// ---------------------------------------------------------------------------
// fe_space<D, T, P> (requirement 17). Generation methods are non-const (A-3):
// they mutate the owned buffers, which is exactly how the zero-allocation
// contract is realized.
// ---------------------------------------------------------------------------
template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class fe_space {
    static_assert(D == 2 || D == 3, "bfem::fe_space: only D == 2 or D == 3");
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;
    typedef fe_function<D, T, P> function_type;
    typedef dual_vector<D, T, P> dual_type;

    // computes ALL element geometries at construction (degeneracy fails fast;
    // the per-program division count settles at exactly one per element)
    fe_space(const mesh<D, T>& msh, int n)
        : mesh_(msh), n_(n), topo_(), dmaps_(), geom_(), op_(),
          buf_(), uloc_(), vloc_(), wloc_(), cws_(), loc_() {
        bfem_scalar_traits<T>::require();   // C-1 contract (L4, additive)
        if (n < 1)
            throw std::invalid_argument("bfem::fe_space: base degree must be >= 1");
        topo_ = detail::fe_space_backend<D>::topology_type::build(msh);
        geom_.reserve(static_cast<std::size_t>(topo_.nt));
        for (int e = 0; e < topo_.nt; ++e) {
            std::array<std::array<T, D>, D + 1> vv;
            for (int k = 0; k <= D; ++k)
                vv[static_cast<std::size_t>(k)] =
                    msh.vertex(msh.element(e)[static_cast<std::size_t>(k)]);
            geom_.push_back(element_geometry<D, T>::from_vertices(vv));
        }
    }

    int base_degree() const { return n_; }
    int ndof(int m) { return dofs(m).ndof(); }
    const dofmap<D>& dofs(int m) {
        if (m < 1)
            throw std::invalid_argument("bfem::fe_space::dofs: m must be >= 1");
        typename std::map<int, dofmap<D> >::iterator it = dmaps_.find(m);
        if (it != dmaps_.end()) return it->second;
        return dmaps_.insert(
                   std::make_pair(m, detail::fe_space_backend<D>::build_dofmap(topo_, m)))
            .first->second;
    }
    int num_elements() const { return topo_.nt; }
    // GRF-1: read-only element geometry, mirroring rt_space / broken_space /
    // vfe_space / c1_space (additive; fe_space previously exposed no
    // per-element geometry)
    const element_geometry<D, T>& geometry(int e) const {
        assert(e >= 0 && e < topo_.nt);
        return geom_[static_cast<std::size_t>(e)];
    }

    // ---- B-1: fe_function factories ----
    function_type zero_function(int m) {
        function_type f(m);
        f.c_.zeros(ndof(m), 1);
        return f;
    }
    function_type function_from_coeffs(int m, vcp::matrix<T, P> c) {
        if (c.rowsize() != ndof(m) || c.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::fe_space::function_from_coeffs: size != ndof(m) x 1");
        function_type f(m);
        f.c_ = std::move(c);
        return f;
    }

    // ---- H3: global stiffness ----
    spmatrix_t stiffness(int m) {
        const dofmap<D>& dm = dofs(m);
        begin_matrix(dm, dm);
#if VCP_BFEM_USE_OPENMP
        const int nt = topo_.nt;
        const int nloc = dm.local_size();
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
                      * static_cast<std::size_t>(nloc)
                      * static_cast<std::size_t>(nloc));
            for (int e = e0; e < e1; ++e) {
                op_l.set_geometry(geom_[static_cast<std::size_t>(e)]);
                op_l.local_stiffness(m, loc_l);
                detail::scatter_matrix<T>(dm, e, loc_l, nloc, b,
                                          typename dofmap<D>::family_tag());
            }
        }
        buf_.append_all(tbuf, nrun);
#else
        for (int e = 0; e < topo_.nt; ++e) {                 // element order (X9)
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            op_.local_stiffness(m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename dofmap<D>::family_tag());
        }
#endif
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }

    // ---- H7: global mixed mass (rows: degree a test side, cols: degree b) ----
    spmatrix_t mixed_mass(int a, int b) {
        const dofmap<D>& dma = dofs(a);
        const dofmap<D>& dmb = dofs(b);
        begin_matrix(dma, dmb);
#if VCP_BFEM_USE_OPENMP
        const int nt = topo_.nt;
        const int nra = dma.local_size();
        const int ncb = dmb.local_size();
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
            detail::coo_buffer<T>& bl = tbuf[static_cast<std::size_t>(tid)];
            bl.reserve(static_cast<std::size_t>(e1 - e0)
                       * static_cast<std::size_t>(nra)
                       * static_cast<std::size_t>(ncb));
            for (int e = e0; e < e1; ++e) {
                op_l.set_geometry(geom_[static_cast<std::size_t>(e)]);
                op_l.local_mass(a, b, loc_l);
                for (int i = 0; i < nra; ++i) {
                    int gi = dma.global_dof(e, i);
                    for (int j = 0; j < ncb; ++j)
                        bl.push(gi, dmb.global_dof(e, j), loc_l(i, j));
                }
            }
        }
        buf_.append_all(tbuf, nrun);
#else
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            op_.local_mass(a, b, loc_);
            for (int i = 0; i < dma.local_size(); ++i) {
                int gi = dma.global_dof(e, i);
                for (int j = 0; j < dmb.local_size(); ++j)
                    buf_.push(gi, dmb.global_dof(e, j), loc_(i, j));
            }
        }
#endif
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dma.ndof(), dmb.ndof(), buf_);
    }

    // ---- H4: load (f(u_h), psi_i), psi in V_h^m ----
    vcp::matrix<T, P> load(const poly1<T>& f, const function_type& uh, int m) {
        validate_fn(uh);
        const dofmap<D>& dmu = dofs(uh.degree());
        const dofmap<D>& dmm = dofs(m);
        vcp::matrix<T, P> F;
        F.zeros(dmm.ndof(), 1);
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            gather_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, f, uloc_, cws_);             // B-5: fe_space buffers
            op_.local_load(wloc_, m, loc_);
            detail::scatter_vector(dmm, e, loc_.data(), dmm.local_size(), F,
                                   typename dofmap<D>::family_tag());
        }
        return F;
    }

    // ---- PF-1 L2-1: load (f, psi_i) for a coordinate polynomial field ----
    // Same element loop / kernel / scatter as load(poly1, uh, m) above; only
    // the composed bpoly is replaced by f.restrict_to. The integrand degree
    // n = max(f.total_degree(), m) (>= total_degree, so restrict_to is
    // exact); local_load's degree argument is the TEST BASIS degree m
    // (element_op.hpp G4: F_alpha = (w, phi^m_alpha)_T, exact for any deg w).
    vcp::matrix<T, P> load(const poly_field<D, T>& f, int m) {
        const dofmap<D>& dmm = dofs(m);
        vcp::matrix<T, P> F;
        F.zeros(dmm.ndof(), 1);
        if (f.is_zero()) return F;
        const int n = f.total_degree() > m ? f.total_degree() : m;
        for (int e = 0; e < topo_.nt; ++e) {                 // element order (X9)
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            wloc_ = f.restrict_to(mesh_, e, n);
            op_.local_load(wloc_, m, loc_);
            detail::scatter_vector(dmm, e, loc_.data(), dmm.local_size(), F,
                                   typename dofmap<D>::family_tag());
        }
        return F;
    }

    // ---- H5: (w psi_j, psi_i), w = f'(u_h) (requirement 16: call with m) ----
    spmatrix_t weighted_mass(const poly1<T>& fprime, const function_type& uh, int m) {
        validate_fn(uh);
        const dofmap<D>& dmu = dofs(uh.degree());
        const dofmap<D>& dm = dofs(m);
        begin_matrix(dm, dm);
#if VCP_BFEM_USE_OPENMP
        const int nt = topo_.nt;
        const int nloc = dm.local_size();
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
            bpoly<D, T> uloc_l, wloc_l;
            compose_workspace<D, T> cws_l;
            detail::coo_buffer<T>& b = tbuf[static_cast<std::size_t>(tid)];
            b.reserve(static_cast<std::size_t>(e1 - e0)
                      * static_cast<std::size_t>(nloc)
                      * static_cast<std::size_t>(nloc));
            for (int e = e0; e < e1; ++e) {
                op_l.set_geometry(geom_[static_cast<std::size_t>(e)]);
                gather_fn(dmu, e, uh, uloc_l);
                compose_into(wloc_l, fprime, uloc_l, cws_l);
                op_l.local_weighted_mass(wloc_l, m, loc_l);
                detail::scatter_matrix<T>(dm, e, loc_l, nloc, b,
                                          typename dofmap<D>::family_tag());
            }
        }
        buf_.append_all(tbuf, nrun);
#else
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            gather_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, fprime, uloc_, cws_);
            op_.local_weighted_mass(wloc_, m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename dofmap<D>::family_tag());
        }
#endif
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }

    // ---- H6: convection (b . grad psi_j, psi_i) ----
    spmatrix_t convection(const std::array<function_type, D>& b, int m) {
        for (int d = 0; d < D; ++d) validate_fn(b[static_cast<std::size_t>(d)]);
        const dofmap<D>& dm = dofs(m);
        begin_matrix(dm, dm);
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            for (int d = 0; d < D; ++d)
                gather_fn(dofs(b[static_cast<std::size_t>(d)].degree()), e,
                          b[static_cast<std::size_t>(d)],
                          bloc_[static_cast<std::size_t>(d)]);
            op_.local_convection(bloc_, m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename dofmap<D>::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }
    spmatrix_t convection(const std::array<T, D>& b_const, int m) {
        const dofmap<D>& dm = dofs(m);
        begin_matrix(dm, dm);
        for (int d = 0; d < D; ++d)
            bloc_[static_cast<std::size_t>(d)] =
                bpoly<D, T>::constant(b_const[static_cast<std::size_t>(d)]);
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            op_.local_convection(bloc_, m, loc_);
            detail::scatter_matrix<T>(dm, e, loc_, dm.local_size(), buf_,
                                      typename dofmap<D>::family_tag());
        }
        buf_.combine();
        return detail::spm_adapter<T, SP>::build(dm.ndof(), dm.ndof(), buf_);
    }

    // ---- H8: global scalars (element-order accumulation, X9) ----
    T inner(const function_type& u, const function_type& v) {
        validate_fn(u);
        validate_fn(v);
        const dofmap<D>& dmu = dofs(u.degree());
        const dofmap<D>& dmv = dofs(v.degree());
        T acc(0);
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            gather_fn(dmu, e, u, uloc_);
            gather_fn(dmv, e, v, vloc_);
            acc += op_.local_inner(uloc_, vloc_);
        }
        return acc;
    }
    T scalar_ff(const poly1<T>& f, const function_type& uh) {
        validate_fn(uh);
        const dofmap<D>& dmu = dofs(uh.degree());
        T acc(0);
        for (int e = 0; e < topo_.nt; ++e) {                 // fixed element order
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            gather_fn(dmu, e, uh, uloc_);
            compose_into(wloc_, f, uloc_, cws_);
            acc += op_.local_inner(wloc_, wloc_);
        }
        return acc;
    }

    // ---- H9 operations ----
    // global elevation with the write-once rule (deterministic overwrite in
    // element order; all writes agree exactly for exact T -- K3)
    function_type elevate(const function_type& u, int m) {
        validate_fn(u);
        if (m < u.degree())
            throw std::invalid_argument("bfem::fe_space::elevate: m < degree(u)");
        const dofmap<D>& dmu = dofs(u.degree());
        const dofmap<D>& dmm = dofs(m);
        function_type out = zero_function(m);
        for (int e = 0; e < topo_.nt; ++e) {
            gather_fn(dmu, e, u, uloc_);
            elevate_into(vloc_, uloc_, m);
            for (int r = 0; r < dmm.local_size(); ++r)
                out.c_(dmm.global_dof(e, r), 0) = vloc_.coeff(r);   // assignment
        }
        return out;
    }

    // dual (requirement 15): u'_i = (u_h, psi_i), element-wise accumulation,
    // no global matrix is formed (M^{-1}-free by construction)
    dual_type dual(const function_type& u, int m) {
        validate_fn(u);
        const dofmap<D>& dmu = dofs(u.degree());
        const dofmap<D>& dmm = dofs(m);
        dual_type d(m);
        d.v_.zeros(dmm.ndof(), 1);
        for (int e = 0; e < topo_.nt; ++e) {
            op_.set_geometry(geom_[static_cast<std::size_t>(e)]);
            gather_fn(dmu, e, u, uloc_);
            op_.local_mass(m, u.degree(), loc_);
            for (int i = 0; i < dmm.local_size(); ++i) {
                T acc(0);
                for (int j = 0; j < dmu.local_size(); ++j)
                    acc += loc_(i, j) * uloc_.coeff(j);
                d.v_(dmm.global_dof(e, i), 0) += acc;
            }
        }
        return d;
    }

    T eval(const function_type& u, int e, const bary_point<D, T>& lam) {
        validate_fn(u);
        if (e < 0 || e >= topo_.nt)
            throw std::invalid_argument("bfem::fe_space::eval: element out of range");
        gather_fn(dofs(u.degree()), e, u, uloc_);
        return ::vcp::bfem::eval(uloc_, lam);
    }

private:
    mesh<D, T> mesh_;
    int n_;
    typename detail::fe_space_backend<D>::topology_type topo_;
    std::map<int, dofmap<D> > dmaps_;
    std::vector<element_geometry<D, T> > geom_;
    element_op<D, T, P> op_;
    detail::coo_buffer<T> buf_;
    bpoly<D, T> uloc_, vloc_, wloc_;
    std::array<bpoly<D, T>, D> bloc_;
    compose_workspace<D, T> cws_;                            // B-5
    vcp::matrix<T, P> loc_;
    friend class fe_function<D, T, P>;

    void validate_fn(const function_type& u) {
        if (u.coeffs().columnsize() != 1 || u.coeffs().rowsize() != ndof(u.degree()))
            throw std::invalid_argument(
                "bfem::fe_space: fe_function does not match this space");
    }
    void gather_fn(const dofmap<D>& dm, int e, const function_type& u,
                   bpoly<D, T>& dst) {
        detail::bpoly_access::prepare(dst, dm.degree(), false);
        detail::gather(dm, e, u.coeffs(),
                       detail::bpoly_access::vec(dst).data(), dm.local_size(),
                       typename dofmap<D>::family_tag());
    }
    void begin_matrix(const dofmap<D>& dma, const dofmap<D>& dmb) {
        buf_.clear();
        buf_.reserve(static_cast<std::size_t>(topo_.nt)
                     * static_cast<std::size_t>(dma.local_size())
                     * static_cast<std::size_t>(dmb.local_size()));
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_FE_SPACE_HPP
