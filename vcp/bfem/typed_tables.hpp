// vcp/bfem/typed_tables.hpp
// Layer 0: T stage of the two stage cache. Holds the image of the rational
// tables converted to T exactly once via convert_traits<T> (enclose-once).
//
// Conforms to: L0 external design v0.3 (section 5) and
//              L0 internal design v0.3 (section 7.2).
//
// Layout mirrors the rational stage with rational replaced by T; accessor
// names are identical.

#ifndef VCP_BFEM_TYPED_TABLES_HPP
#define VCP_BFEM_TYPED_TABLES_HPP

#include <vector>
#include <map>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/detail/table_cache.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T> class typed_registry;

// ---------------------------------------------------------------------------
template <int D, typename T>
class typed_mass_table {
public:
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    const T& at(int i, int j) const {
        assert(i >= 0 && i < rows_ && j >= 0 && j < cols_);
        return v_[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols_)
                  + static_cast<std::size_t>(j)];
    }
    int degree_row() const { return a_; }
    int degree_col() const { return b_; }

private:
    friend class typed_registry<D, T>;
    typed_mass_table(int a, int b, int rows, int cols, std::vector<T> v)
        : a_(a), b_(b), rows_(rows), cols_(cols), v_(std::move(v)) {}
    int a_, b_, rows_, cols_;
    std::vector<T> v_;
};

// ---------------------------------------------------------------------------
template <int D, typename T>
class typed_elevation_table {
public:
    struct entry {
        int target_rank;
        const T* coeff;
    };

    class entry_iterator {
    public:
        entry_iterator(const int* t, const T* c) : t_(t), c_(c) {}
        entry operator*() const { entry e = { *t_, c_ }; return e; }
        entry_iterator& operator++() { ++t_; ++c_; return *this; }
        bool operator!=(const entry_iterator& o) const { return t_ != o.t_; }
        bool operator==(const entry_iterator& o) const { return t_ == o.t_; }
    private:
        const int* t_;
        const T* c_;
    };

    struct entry_range {
        entry_iterator b, e;
        entry_iterator begin() const { return b; }
        entry_iterator end() const { return e; }
    };

    int n() const { return n_; }
    int m() const { return m_; }
    int source_size() const { return source_size_; }
    int target_size() const { return target_size_; }
    int row_length() const { return row_len_; }

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
    friend class typed_registry<D, T>;
    typed_elevation_table(int n, int m, int source_size, int target_size,
                          int row_len, std::vector<int> targets, std::vector<T> coeffs)
        : n_(n), m_(m), source_size_(source_size), target_size_(target_size),
          row_len_(row_len), targets_(std::move(targets)), coeffs_(std::move(coeffs)) {}
    int n_, m_, source_size_, target_size_, row_len_;
    std::vector<int> targets_;
    std::vector<T> coeffs_;
};

// ---------------------------------------------------------------------------
template <int D, typename T>
class typed_product_table {
public:
    int a() const { return a_; }
    int b() const { return b_; }
    int rows() const { return rows_; }
    int cols() const { return cols_; }
    const T& coeff(int i, int j) const {
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
    friend class typed_registry<D, T>;
    typed_product_table(int a, int b, int rows, int cols,
                        std::vector<T> c, std::vector<int> t)
        : a_(a), b_(b), rows_(rows), cols_(cols),
          c_(std::move(c)), t_(std::move(t)) {}
    int a_, b_, rows_, cols_;
    std::vector<T> c_;
    std::vector<int> t_;
};

// ---------------------------------------------------------------------------
// typed_registry<D, T>: converts rational stage tables to T exactly once per
// (D, T, table key) and caches the result. If the rational table already
// exists no rational arithmetic happens here (enclose-once implementation
// point). Each T has its own mutex; different T initialize concurrently.
// ---------------------------------------------------------------------------
template <int D, typename T>
class typed_registry {
public:
    static const T& basis_integral(int n) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.integrals.get_or_build(n, [&]() -> T {
            const rational& w = coeff_registry<D>::basis_integral(n);
            T v = conv(w);
            return v;
        });
    }

    static const typed_mass_table<D, T>& mass(int a, int b) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        return s.mass.get_or_build(key, [&]() -> typed_mass_table<D, T> {
            const mass_table<D>& src = coeff_registry<D>::mass(a, b);
            std::vector<T> v;
            v.reserve(static_cast<std::size_t>(src.rows()) * static_cast<std::size_t>(src.cols()));
            for (int i = 0; i < src.rows(); ++i)
                for (int j = 0; j < src.cols(); ++j)
                    v.push_back(conv(src.at(i, j)));
            typed_mass_table<D, T> tbl(a, b, src.rows(), src.cols(), std::move(v));
            return tbl;
        });
    }

    static const typed_elevation_table<D, T>& elevation(int n, int m) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(n, m);
        return s.elev.get_or_build(key, [&]() -> typed_elevation_table<D, T> {
            const elevation_table<D>& src = coeff_registry<D>::elevation(n, m);
            std::vector<int> targets;
            std::vector<T> coeffs;
            std::size_t total = static_cast<std::size_t>(src.source_size())
                                * static_cast<std::size_t>(src.row_length());
            targets.reserve(total);
            coeffs.reserve(total);
            for (int i = 0; i < src.source_size(); ++i) {
                typename elevation_table<D>::entry_range rr = src.row(i);
                for (typename elevation_table<D>::entry_iterator p = rr.begin();
                     p != rr.end(); ++p) {
                    typename elevation_table<D>::entry e = *p;
                    targets.push_back(e.target_rank);
                    coeffs.push_back(conv(*e.coeff));
                }
            }
            typed_elevation_table<D, T> tbl(n, m, src.source_size(), src.target_size(),
                                            src.row_length(),
                                            std::move(targets), std::move(coeffs));
            return tbl;
        });
    }

    static const typed_product_table<D, T>& product(int a, int b) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        return s.prod.get_or_build(key, [&]() -> typed_product_table<D, T> {
            const product_table<D>& src = coeff_registry<D>::product(a, b);
            std::vector<T> c;
            std::vector<int> t;
            std::size_t total = static_cast<std::size_t>(src.rows())
                                * static_cast<std::size_t>(src.cols());
            c.reserve(total);
            t.reserve(total);
            for (int i = 0; i < src.rows(); ++i) {
                for (int j = 0; j < src.cols(); ++j) {
                    c.push_back(conv(src.coeff(i, j)));
                    t.push_back(src.target_rank(i, j));
                }
            }
            typed_product_table<D, T> tbl(a, b, src.rows(), src.cols(),
                                          std::move(c), std::move(t));
            return tbl;
        });
    }

private:
    // L4: shared detail::table_cache vessel (one mutex per D x T, unchanged)
    struct maps {
        detail::table_cache<int, T> integrals;
        detail::table_cache<std::pair<int, int>, typed_mass_table<D, T> > mass;
        detail::table_cache<std::pair<int, int>, typed_elevation_table<D, T> > elev;
        detail::table_cache<std::pair<int, int>, typed_product_table<D, T> > prod;
    };
    typedef detail::table_cache_state<maps> state;
    static state& st() {
        return detail::table_cache_instance<state>();
    }
    static T conv(const rational& r) {
        return convert_traits<T>::from_rational(r.num(), r.den());
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_TYPED_TABLES_HPP
