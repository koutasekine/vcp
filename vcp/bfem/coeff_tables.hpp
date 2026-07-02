// vcp/bfem/coeff_tables.hpp
// Layer 0: exact rational coefficient tables (T1-T5) and the shared,
// type independent registry coeff_registry<D>.
//
// Conforms to: L0 external design v0.3 (sections 3.4, 4) and
//              L0 internal design v0.3 (sections 5, 7).
//
// All tables are immutable after construction and normalized as coefficients
// of |T| (reference simplex values). No floating point appears here.

#ifndef VCP_BFEM_COEFF_TABLES_HPP
#define VCP_BFEM_COEFF_TABLES_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/multi_index.hpp>

namespace vcp {
namespace bfem {

template <int D> class coeff_registry;

// ---------------------------------------------------------------------------
// T2: mixed degree mass table (dense).  (M^{a,b})_{ij} = integral of
// B^a_alpha B^b_beta over the reference simplex, divided by |T|.
// ---------------------------------------------------------------------------
template <int D>
class mass_table {
public:
    int rows() const { return rows_; }                 // N(D, a)
    int cols() const { return cols_; }                 // N(D, b)
    const rational& at(int i, int j) const {
        assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
        return v_[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)
                  + static_cast<std::size_t>(j)];
    }
    int degree_row() const { return a_; }
    int degree_col() const { return b_; }

private:
    friend class coeff_registry<D>;
    mass_table(int a, int b, int rows, int cols, std::vector<rational> v)
        : a_(a), b_(b), rows_(rows), cols_(cols), v_(std::move(v)) {}
    int a_, b_, rows_, cols_;
    std::vector<rational> v_;
};

// ---------------------------------------------------------------------------
// T3: degree elevation table (sparse, scatter rows per source alpha).
// Application semantics: c^m = E c^n, (c^m)_beta = sum_alpha E_{beta,alpha} (c^n)_alpha.
// Row i (source rank) lists the C(m-n+D, D) targets it contributes to.
// Row pointers are implicit: row length is constant (internal design 5.3).
// ---------------------------------------------------------------------------
template <int D>
class elevation_table {
public:
    struct entry {
        int target_rank;
        const rational* coeff;
    };

    class entry_iterator {
    public:
        entry_iterator(const int* t, const rational* c) : t_(t), c_(c) {}
        entry operator*() const { entry e = { *t_, c_ }; return e; }
        entry_iterator& operator++() { ++t_; ++c_; return *this; }
        bool operator!=(const entry_iterator& o) const { return t_ != o.t_; }
        bool operator==(const entry_iterator& o) const { return t_ == o.t_; }
    private:
        const int* t_;
        const rational* c_;
    };

    struct entry_range {
        entry_iterator b, e;
        entry_iterator begin() const { return b; }
        entry_iterator end() const { return e; }
    };

    int n() const { return n_; }
    int m() const { return m_; }
    int source_size() const { return source_size_; }   // N(D, n)
    int target_size() const { return target_size_; }   // N(D, m)
    int row_length() const { return row_len_; }        // C(m-n+D, D)

    // contributions of source rank i (row length entries)
    entry_range row(int i) const {
        assert(i >= 0 && i < source_size_);
        std::size_t off = static_cast<std::size_t>(i) * static_cast<std::size_t>(row_len_);
        entry_range r = {
            entry_iterator(&targets_[off], &coeffs_[off]),
            entry_iterator(&targets_[off] + row_len_, &coeffs_[off] + row_len_)
        };
        return r;
    }

private:
    friend class coeff_registry<D>;
    elevation_table(int n, int m, int source_size, int target_size, int row_len,
                    std::vector<int> targets, std::vector<rational> coeffs)
        : n_(n), m_(m), source_size_(source_size), target_size_(target_size),
          row_len_(row_len), targets_(std::move(targets)), coeffs_(std::move(coeffs)) {}
    int n_, m_, source_size_, target_size_, row_len_;
    std::vector<int> targets_;
    std::vector<rational> coeffs_;
};

// ---------------------------------------------------------------------------
// T4: product coefficient table (dense coefficients + precomputed target
// ranks).  B^a_alpha * B^b_beta = c(alpha,beta) * B^{a+b}_{alpha+beta}.
// The precomputed int targets remove all rank computation from the Layer 1
// convolution loop (internal design 5.4).
// ---------------------------------------------------------------------------
template <int D>
class product_table {
public:
    int a() const { return a_; }
    int b() const { return b_; }
    int rows() const { return rows_; }                 // N(D, a)
    int cols() const { return cols_; }                 // N(D, b)
    const rational& coeff(int i, int j) const {
        assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
        return c_[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)
                  + static_cast<std::size_t>(j)];
    }
    int target_rank(int i, int j) const {
        assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
        return t_[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)
                  + static_cast<std::size_t>(j)];
    }

private:
    friend class coeff_registry<D>;
    product_table(int a, int b, int rows, int cols,
                  std::vector<rational> c, std::vector<int> t)
        : a_(a), b_(b), rows_(rows), cols_(cols),
          c_(std::move(c)), t_(std::move(t)) {}
    int a_, b_, rows_, cols_;
    std::vector<rational> c_;
    std::vector<int> t_;
};

// ---------------------------------------------------------------------------
// T5: derivative structure map. dB^n_alpha/dlambda_i = n * B^{n-1}_{alpha-e_i};
// only the integer factor n and the index shift are needed (no table of
// rationals). target(r, i) < 0 means alpha_i == 0 (the term vanishes).
// All (r, i) pairs are precomputed to an int array at construction.
// ---------------------------------------------------------------------------
template <int D>
class derivative_map {
public:
    explicit derivative_map(int n) : n_(n), size_(0), t_() {
        if (n < 0)
            throw std::invalid_argument("bfem::derivative_map: negative degree");
        const index_map<D>& im_n = coeff_registry<D>::indices(n);
        size_ = im_n.size();
        t_.assign(static_cast<std::size_t>(size_) * static_cast<std::size_t>(D + 1), -1);
        if (n >= 1) {
            const index_map<D>& im_lo = coeff_registry<D>::indices(n - 1);
            for (int r = 0; r < size_; ++r) {
                multi_index<D> alpha = im_n.unrank(r);
                for (int i = 0; i <= D; ++i) {
                    std::pair<multi_index<D>, bool> me = alpha.minus_e(i);
                    if (me.second)
                        t_[static_cast<std::size_t>(r) * static_cast<std::size_t>(D + 1)
                           + static_cast<std::size_t>(i)] = im_lo.rank(me.first);
                }
            }
        }
    }

    int source_size() const { return size_; }          // N(D, n)
    // < 0 means the term vanishes (alpha_i == 0)
    int target(int source_rank, int i) const {
        assert(source_rank >= 0 && source_rank < size_ && i >= 0 && i <= D);
        return t_[static_cast<std::size_t>(source_rank) * static_cast<std::size_t>(D + 1)
                  + static_cast<std::size_t>(i)];
    }
    int factor() const { return n_; }                  // the integer n

private:
    int n_;
    int size_;
    std::vector<int> t_;
};

// ---------------------------------------------------------------------------
// coeff_registry<D>: type independent, process shared registry (rational
// stage of the two stage cache). All members are lazily generated and cached;
// returned references stay valid until program termination.
// Thread safety: one mutex per D; generation happens while holding the lock
// (internal design 7.1).
// ---------------------------------------------------------------------------
template <int D>
class coeff_registry {
public:
    // T1: w(D, n) = 1 / C(n + D, D)
    static const rational& basis_integral(int n) {
        check_deg(n, "basis_integral");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        typename std::map<int, rational>::iterator it = s.integrals.find(n);
        if (it != s.integrals.end()) return it->second;
        rational w(detail::bigint(1), detail::binomial_cache::binom(n + D, D));
        return s.integrals.insert(std::make_pair(n, std::move(w))).first->second;
    }

    // T2: precondition a >= 0, b >= 0
    static const mass_table<D>& mass(int a, int b) {
        check_deg(a, "mass");
        check_deg(b, "mass");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        typename std::map<std::pair<int, int>, mass_table<D> >::iterator it =
            s.mass.find(key);
        if (it != s.mass.end()) return it->second;
        mass_table<D> tbl = build_mass(s, a, b);
        return s.mass.insert(std::make_pair(key, std::move(tbl))).first->second;
    }

    // T3: precondition m >= n >= 0 (m == n returns the identity table)
    static const elevation_table<D>& elevation(int n, int m) {
        check_deg(n, "elevation");
        if (m < n)
            throw std::invalid_argument("bfem::coeff_registry::elevation: m < n");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(n, m);
        typename std::map<std::pair<int, int>, elevation_table<D> >::iterator it =
            s.elev.find(key);
        if (it != s.elev.end()) return it->second;
        elevation_table<D> tbl = build_elevation(s, n, m);
        return s.elev.insert(std::make_pair(key, std::move(tbl))).first->second;
    }

    // T4: precondition a >= 0, b >= 0
    static const product_table<D>& product(int a, int b) {
        check_deg(a, "product");
        check_deg(b, "product");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        typename std::map<std::pair<int, int>, product_table<D> >::iterator it =
            s.prod.find(key);
        if (it != s.prod.end()) return it->second;
        product_table<D> tbl = build_product(s, a, b);
        return s.prod.insert(std::make_pair(key, std::move(tbl))).first->second;
    }

    // U1: shared index_map cache (every layer references the same instance)
    static const index_map<D>& indices(int n) {
        check_deg(n, "indices");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return indices_locked(s, n);
    }

private:
    struct state {
        std::mutex mtx;
        std::map<int, index_map<D> > imaps;
        std::map<int, rational> integrals;
        std::map<std::pair<int, int>, mass_table<D> > mass;
        std::map<std::pair<int, int>, elevation_table<D> > elev;
        std::map<std::pair<int, int>, product_table<D> > prod;
    };
    static state& st() {
        static state s;   // C++11 magic static: thread safe initialization
        return s;
    }
    static void check_deg(int n, const char* where) {
        if (n < 0) {
            std::string msg("bfem::coeff_registry::");
            msg += where;
            msg += ": negative degree";
            throw std::invalid_argument(msg);
        }
    }

    // callers must hold s.mtx
    static const index_map<D>& indices_locked(state& s, int n) {
        typename std::map<int, index_map<D> >::iterator it = s.imaps.find(n);
        if (it != s.imaps.end()) return it->second;
        return s.imaps.insert(std::make_pair(n, index_map<D>(n))).first->second;
    }

    // unrank all multi indices of degree n once (internal design 5.2:
    // repeated unrank calls inside generation loops are forbidden)
    static std::vector<multi_index<D> > all_indices(const index_map<D>& im) {
        std::vector<multi_index<D> > v;
        v.reserve(static_cast<std::size_t>(im.size()));
        for (int r = 0; r < im.size(); ++r) v.push_back(im.unrank(r));
        return v;
    }

    static detail::bigint mult(int n, const multi_index<D>& alpha) {
        return detail::multinomial(n, alpha.a.data(), D + 1);
    }

    // (M^{a,b})_{alpha,beta} = mult(a,alpha) * mult(b,beta)
    //                          / ( mult(a+b, alpha+beta) * C(a+b+D, D) )
    static mass_table<D> build_mass(state& s, int a, int b) {
        const index_map<D>& im_a = indices_locked(s, a);
        const index_map<D>& im_b = indices_locked(s, b);
        std::vector<multi_index<D> > as = all_indices(im_a);
        std::vector<multi_index<D> > bs = all_indices(im_b);
        std::vector<detail::bigint> mb(bs.size());
        for (std::size_t j = 0; j < bs.size(); ++j) mb[j] = mult(b, bs[j]);
        const detail::bigint& cbd = detail::binomial_cache::binom(a + b + D, D);
        std::vector<rational> v;
        v.reserve(as.size() * bs.size());
        for (std::size_t i = 0; i < as.size(); ++i) {
            detail::bigint ma = mult(a, as[i]);
            for (std::size_t j = 0; j < bs.size(); ++j) {
                multi_index<D> ab = as[i] + bs[j];
                v.push_back(rational(ma * mb[j], mult(a + b, ab) * cbd));
            }
        }
        return mass_table<D>(a, b, im_a.size(), im_b.size(), std::move(v));
    }

    // E_{beta,alpha} = mult(n,alpha) * mult(m-n, beta-alpha) / mult(m,beta)
    static elevation_table<D> build_elevation(state& s, int n, int m) {
        const index_map<D>& im_n = indices_locked(s, n);
        const index_map<D>& im_m = indices_locked(s, m);
        int r = m - n;
        const index_map<D>& im_r = indices_locked(s, r);
        std::vector<multi_index<D> > ns = all_indices(im_n);
        std::vector<multi_index<D> > ds = all_indices(im_r);
        std::vector<detail::bigint> md(ds.size());
        for (std::size_t k = 0; k < ds.size(); ++k) md[k] = mult(r, ds[k]);
        int row_len = im_r.size();                   // C(m-n+D, D)
        std::vector<int> targets;
        std::vector<rational> coeffs;
        targets.reserve(ns.size() * static_cast<std::size_t>(row_len));
        coeffs.reserve(ns.size() * static_cast<std::size_t>(row_len));
        for (std::size_t i = 0; i < ns.size(); ++i) {
            detail::bigint mn = mult(n, ns[i]);
            for (std::size_t k = 0; k < ds.size(); ++k) {
                multi_index<D> beta = ns[i] + ds[k];
                targets.push_back(im_m.rank(beta));
                coeffs.push_back(rational(mn * md[k], mult(m, beta)));
            }
        }
        return elevation_table<D>(n, m, im_n.size(), im_m.size(), row_len,
                                  std::move(targets), std::move(coeffs));
    }

    // c(alpha,beta) = mult(a,alpha) * mult(b,beta) / mult(a+b, alpha+beta)
    static product_table<D> build_product(state& s, int a, int b) {
        const index_map<D>& im_a = indices_locked(s, a);
        const index_map<D>& im_b = indices_locked(s, b);
        const index_map<D>& im_ab = indices_locked(s, a + b);
        std::vector<multi_index<D> > as = all_indices(im_a);
        std::vector<multi_index<D> > bs = all_indices(im_b);
        std::vector<detail::bigint> mb(bs.size());
        for (std::size_t j = 0; j < bs.size(); ++j) mb[j] = mult(b, bs[j]);
        std::vector<rational> c;
        std::vector<int> t;
        c.reserve(as.size() * bs.size());
        t.reserve(as.size() * bs.size());
        for (std::size_t i = 0; i < as.size(); ++i) {
            detail::bigint ma = mult(a, as[i]);
            for (std::size_t j = 0; j < bs.size(); ++j) {
                multi_index<D> ab = as[i] + bs[j];
                c.push_back(rational(ma * mb[j], mult(a + b, ab)));
                t.push_back(im_ab.rank(ab));
            }
        }
        return product_table<D>(a, b, im_a.size(), im_b.size(),
                                std::move(c), std::move(t));
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_COEFF_TABLES_HPP
