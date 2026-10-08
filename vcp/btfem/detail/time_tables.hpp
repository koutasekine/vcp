// vcp/btfem/detail/time_tables.hpp
// btfem T3 (internal): reference-interval [0,1] time tables synthesized from
// the bfem L0 D = 1 registry (coeff_registry<1>: mass / product,
// derivative_map<1>), in exact rationals, then converted to T exactly once.
//
// Conforms to: BT-1 external design v1.1 (section 4.2) and
//              BT-1 internal design v1.0 (section 4).
//
// No new rational formula is written here: every entry is a finite sum of
// products of L0 table entries (mass, product) and the integer derivative
// matrix D^q obtained from derivative_map<1>.
//
// Two shared stages (process wide, lazily generated, bfem table_cache vessel):
//   rational stage  time_tables_rational          keyed by q / (s, q) / (q, qb)
//   T stage         time_tables_typed<T>          same keys, vcp::matrix<T, mats<T>>
// The physical factors k and 1/k (internal 4.2) are applied by time_space on
// its own copies (multiplications only); the T-stage conversion divisions
// belong to the L0 conversion budget (same convention as typed_registry).
//
// Conventions (external 4.1, verified in P0): local rank r <-> alpha = (q-r, r),
// B^q_r = C(q,r) (1-tau)^{q-r} tau^r, tau = lambda_1;
// d/dtau B^q_r = q (B^{q-1}_{r-1} - B^{q-1}_r) (vanishing terms dropped).
// All (a, b)-indexed tables below are stored with the TEST (derivative) index
// as the ROW a and the trial index as the COLUMN b.

#ifndef VCP_BTFEM_DETAIL_TIME_TABLES_HPP
#define VCP_BTFEM_DETAIL_TIME_TABLES_HPP

#include <vector>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cassert>

#include <vcp/matrix.hpp>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/coeff_tables.hpp>
#include <vcp/bfem/convert_traits.hpp>
#include <vcp/bfem/detail/table_cache.hpp>

namespace vcp {
namespace btfem {
namespace detail {

// ---------------------------------------------------------------------------
// small dense rational matrix (row major)
// ---------------------------------------------------------------------------
struct rmat {
    int rows, cols;
    std::vector<bfem::rational> v;

    rmat() : rows(0), cols(0), v() {}
    rmat(int r, int c) : rows(r), cols(c),
        v(static_cast<std::size_t>(r) * static_cast<std::size_t>(c), bfem::rational(0)) {}
    bfem::rational& at(int i, int j) {
        assert(i >= 0 && i < rows && j >= 0 && j < cols);
        return v[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)
                 + static_cast<std::size_t>(j)];
    }
    const bfem::rational& at(int i, int j) const {
        assert(i >= 0 && i < rows && j >= 0 && j < cols);
        return v[static_cast<std::size_t>(i) * static_cast<std::size_t>(cols)
                 + static_cast<std::size_t>(j)];
    }
};

// ---------------------------------------------------------------------------
// integer derivative matrix D^q: (q) x (q+1); D[rp][r] = coefficient of
// B^{q-1}_{rp} in d/dtau B^q_r. Built from derivative_map<1>(q): direction
// lambda_1 contributes +q at target(r,1), direction lambda_0 contributes -q at
// target(r,0) (d/dtau = d/dlambda_1 - d/dlambda_0).
// ---------------------------------------------------------------------------
struct deriv_matrix {
    int q;
    std::vector<int> d;      // q x (q+1), row major

    deriv_matrix() : q(0), d() {}
    explicit deriv_matrix(int qq) : q(qq), d() {
        if (qq < 1)
            throw std::invalid_argument("btfem::deriv_matrix: q < 1");
        d.assign(static_cast<std::size_t>(qq) * static_cast<std::size_t>(qq + 1), 0);
        bfem::derivative_map<1> dm(qq);
        for (int r = 0; r <= qq; ++r) {
            int t1 = dm.target(r, 1);
            int t0 = dm.target(r, 0);
            if (t1 >= 0) at(t1, r) += dm.factor();
            if (t0 >= 0) at(t0, r) -= dm.factor();
        }
    }
    int& at(int rp, int r) {
        assert(rp >= 0 && rp < q && r >= 0 && r <= q);
        return d[static_cast<std::size_t>(rp) * static_cast<std::size_t>(q + 1)
                 + static_cast<std::size_t>(r)];
    }
    int at(int rp, int r) const {
        assert(rp >= 0 && rp < q && r >= 0 && r <= q);
        return d[static_cast<std::size_t>(rp) * static_cast<std::size_t>(q + 1)
                 + static_cast<std::size_t>(r)];
    }
};

// ---------------------------------------------------------------------------
// rational stage (process shared)
// ---------------------------------------------------------------------------
class time_tables_rational {
public:
    // D^q (q >= 1)
    static const deriv_matrix& deriv(int q) {
        check_q(q, "deriv");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.deriv.get_or_build(q, [&]() -> deriv_matrix { return deriv_matrix(q); });
    }
    // Mref(a, b)[i][j] = int_0^1 B^a_i B^b_j : the L0 table itself
    static const bfem::mass_table<1>& mass_ref(int a, int b) {
        return bfem::coeff_registry<1>::mass(a, b);
    }
    // Sref(q)[a][b] = int (B^q_b)' (B^q_a)'  (q >= 1)
    static const rmat& stiffness_ref(int q) {
        check_q(q, "stiffness_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.stiff.get_or_build(q, [&]() -> rmat { return build_mixed_stiffness(q, q); });
    }
    // Bref(q)[a][b] = int B^q_b (B^q_a)'  (q >= 1)
    static const rmat& dt_mass_ref(int q) {
        check_q(q, "dt_mass_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.dtmass.get_or_build(q, [&]() -> rmat { return build_mixed_dt_mass(q, q); });
    }
    // Lref(s, q)[sigma][a] = int B^s_sigma (B^q_a)'  ((s+1) x (q+1), s >= 0, q >= 1)
    static const rmat& dt_load_ref(int sdeg, int q) {
        check_s(sdeg, "dt_load_ref");
        check_q(q, "dt_load_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(sdeg, q);
        return s.dtload.get_or_build(key, [&]() -> rmat { return build_dt_load(sdeg, q); });
    }
    // Wref(s, q)[sigma][a][b] = int B^s_sigma B^q_b (B^q_a)'
    static const std::vector<rmat>& weighted_dt_mass_ref(int sdeg, int q) {
        check_s(sdeg, "weighted_dt_mass_ref");
        check_q(q, "weighted_dt_mass_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(sdeg, q);
        return s.wdt.get_or_build(key, [&]() -> std::vector<rmat> {
            return build_weighted_dt_mass(sdeg, q);
        });
    }
    // Uref(s, q)[sigma][a][b] = int B^s_sigma B^q_b B^q_a
    static const std::vector<rmat>& weighted_mass_ref(int sdeg, int q) {
        check_s(sdeg, "weighted_mass_ref");
        check_q(q, "weighted_mass_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(sdeg, q);
        return s.wm.get_or_build(key, [&]() -> std::vector<rmat> {
            return build_weighted_mass(sdeg, q);
        });
    }
    // mixed Sref(q, qb)[a][b] = int (B^qb_b)' (B^q_a)'   ((q+1) x (qb+1))
    static const rmat& mixed_stiffness_ref(int q, int qb) {
        check_q(q, "mixed_stiffness_ref");
        check_q(qb, "mixed_stiffness_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(q, qb);
        return s.mstiff.get_or_build(key, [&]() -> rmat { return build_mixed_stiffness(q, qb); });
    }
    // mixed Bref(q, qb)[a][b] = int B^qb_b (B^q_a)'   ((q+1) x (qb+1))
    static const rmat& mixed_dt_mass_ref(int q, int qb) {
        check_q(q, "mixed_dt_mass_ref");
        check_q(qb, "mixed_dt_mass_ref");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(q, qb);
        return s.mdtmass.get_or_build(key, [&]() -> rmat { return build_mixed_dt_mass(q, qb); });
    }

private:
    struct maps {
        bfem::detail::table_cache<int, deriv_matrix> deriv;
        bfem::detail::table_cache<int, rmat> stiff;
        bfem::detail::table_cache<int, rmat> dtmass;
        bfem::detail::table_cache<std::pair<int, int>, rmat> dtload;
        bfem::detail::table_cache<std::pair<int, int>, std::vector<rmat> > wdt;
        bfem::detail::table_cache<std::pair<int, int>, std::vector<rmat> > wm;
        bfem::detail::table_cache<std::pair<int, int>, rmat> mstiff;
        bfem::detail::table_cache<std::pair<int, int>, rmat> mdtmass;
    };
    typedef bfem::detail::table_cache_state<maps> state;
    static state& st() { return bfem::detail::table_cache_instance<state>(); }

    static void check_q(int q, const char* where) {
        if (q < 1) {
            std::string msg("btfem::time_tables_rational::");
            msg += where;
            msg += ": q < 1";
            throw std::invalid_argument(msg);
        }
    }
    static void check_s(int s, const char* where) {
        if (s < 0) {
            std::string msg("btfem::time_tables_rational::");
            msg += where;
            msg += ": sdeg < 0";
            throw std::invalid_argument(msg);
        }
    }

    // builders (called with the stage lock held; they only read the L0
    // registry, whose own lock is independent)
    static rmat build_mixed_stiffness(int q, int qb) {
        deriv_matrix Dq(q), Db(qb);
        const bfem::mass_table<1>& M = mass_ref(q - 1, qb - 1);
        rmat out(q + 1, qb + 1);
        for (int a = 0; a <= q; ++a) {
            for (int b = 0; b <= qb; ++b) {
                bfem::rational acc(0);
                for (int rp = 0; rp < q; ++rp) {
                    int da = Dq.at(rp, a);
                    if (da == 0) continue;
                    for (int rpp = 0; rpp < qb; ++rpp) {
                        int db = Db.at(rpp, b);
                        if (db == 0) continue;
                        acc += bfem::rational(static_cast<long long>(da) * db) * M.at(rp, rpp);
                    }
                }
                out.at(a, b) = acc;
            }
        }
        return out;
    }
    static rmat build_mixed_dt_mass(int q, int qb) {
        deriv_matrix Dq(q);
        const bfem::mass_table<1>& M = mass_ref(q - 1, qb);
        rmat out(q + 1, qb + 1);
        for (int a = 0; a <= q; ++a) {
            for (int b = 0; b <= qb; ++b) {
                bfem::rational acc(0);
                for (int rp = 0; rp < q; ++rp) {
                    int da = Dq.at(rp, a);
                    if (da == 0) continue;
                    acc += bfem::rational(static_cast<long long>(da)) * M.at(rp, b);
                }
                out.at(a, b) = acc;
            }
        }
        return out;
    }
    static rmat build_dt_load(int sdeg, int q) {
        deriv_matrix Dq(q);
        const bfem::mass_table<1>& M = mass_ref(q - 1, sdeg);
        rmat out(sdeg + 1, q + 1);
        for (int sg = 0; sg <= sdeg; ++sg) {
            for (int a = 0; a <= q; ++a) {
                bfem::rational acc(0);
                for (int rp = 0; rp < q; ++rp) {
                    int da = Dq.at(rp, a);
                    if (da == 0) continue;
                    acc += bfem::rational(static_cast<long long>(da)) * M.at(rp, sg);
                }
                out.at(sg, a) = acc;
            }
        }
        return out;
    }
    static std::vector<rmat> build_weighted_dt_mass(int sdeg, int q) {
        deriv_matrix Dq(q);
        const bfem::product_table<1>& Pt = bfem::coeff_registry<1>::product(sdeg, q);
        const bfem::mass_table<1>& M = mass_ref(q - 1, sdeg + q);
        std::vector<rmat> out(static_cast<std::size_t>(sdeg + 1), rmat(q + 1, q + 1));
        for (int sg = 0; sg <= sdeg; ++sg) {
            rmat& W = out[static_cast<std::size_t>(sg)];
            for (int a = 0; a <= q; ++a) {
                for (int b = 0; b <= q; ++b) {
                    const bfem::rational& c = Pt.coeff(sg, b);
                    int tgt = Pt.target_rank(sg, b);
                    bfem::rational acc(0);
                    for (int rp = 0; rp < q; ++rp) {
                        int da = Dq.at(rp, a);
                        if (da == 0) continue;
                        acc += bfem::rational(static_cast<long long>(da)) * M.at(rp, tgt);
                    }
                    W.at(a, b) = c * acc;
                }
            }
        }
        return out;
    }
    static std::vector<rmat> build_weighted_mass(int sdeg, int q) {
        const bfem::product_table<1>& Pt = bfem::coeff_registry<1>::product(sdeg, q);
        const bfem::mass_table<1>& M = mass_ref(q, sdeg + q);
        std::vector<rmat> out(static_cast<std::size_t>(sdeg + 1), rmat(q + 1, q + 1));
        for (int sg = 0; sg <= sdeg; ++sg) {
            rmat& U = out[static_cast<std::size_t>(sg)];
            for (int a = 0; a <= q; ++a)
                for (int b = 0; b <= q; ++b)
                    U.at(a, b) = Pt.coeff(sg, b) * M.at(a, Pt.target_rank(sg, b));
        }
        return out;
    }
};

// ---------------------------------------------------------------------------
// T stage (process shared per T): each rational entry converted exactly once
// through convert_traits<T>::from_rational (enclose-once). Reference interval
// values -- no k factor yet.
// ---------------------------------------------------------------------------
template <typename T>
class time_tables_typed {
public:
    typedef vcp::matrix<T, vcp::mats<T> > table_type;

    static const table_type& mass_ref(int a, int b) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(a, b);
        return s.mass.get_or_build(key, [&]() -> table_type {
            return conv_mass(time_tables_rational::mass_ref(a, b));
        });
    }
    static const table_type& stiffness_ref(int q) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.stiff.get_or_build(q, [&]() -> table_type {
            return conv(time_tables_rational::stiffness_ref(q));
        });
    }
    static const table_type& dt_mass_ref(int q) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return s.dtmass.get_or_build(q, [&]() -> table_type {
            return conv(time_tables_rational::dt_mass_ref(q));
        });
    }
    static const table_type& dt_load_ref(int sdeg, int q) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(sdeg, q);
        return s.dtload.get_or_build(key, [&]() -> table_type {
            return conv(time_tables_rational::dt_load_ref(sdeg, q));
        });
    }
    static const std::vector<table_type>& weighted_dt_mass_ref(int sdeg, int q) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(sdeg, q);
        return s.wdt.get_or_build(key, [&]() -> std::vector<table_type> {
            return conv_vec(time_tables_rational::weighted_dt_mass_ref(sdeg, q));
        });
    }
    static const std::vector<table_type>& weighted_mass_ref(int sdeg, int q) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(sdeg, q);
        return s.wm.get_or_build(key, [&]() -> std::vector<table_type> {
            return conv_vec(time_tables_rational::weighted_mass_ref(sdeg, q));
        });
    }
    static const table_type& mixed_stiffness_ref(int q, int qb) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(q, qb);
        return s.mstiff.get_or_build(key, [&]() -> table_type {
            return conv(time_tables_rational::mixed_stiffness_ref(q, qb));
        });
    }
    static const table_type& mixed_dt_mass_ref(int q, int qb) {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        std::pair<int, int> key(q, qb);
        return s.mdtmass.get_or_build(key, [&]() -> table_type {
            return conv(time_tables_rational::mixed_dt_mass_ref(q, qb));
        });
    }

    // the enclose-once conversion of one rational
    static T conv1(const bfem::rational& r) {
        return bfem::convert_traits<T>::from_rational(r.num(), r.den());
    }

private:
    struct maps {
        bfem::detail::table_cache<std::pair<int, int>, table_type> mass;
        bfem::detail::table_cache<int, table_type> stiff;
        bfem::detail::table_cache<int, table_type> dtmass;
        bfem::detail::table_cache<std::pair<int, int>, table_type> dtload;
        bfem::detail::table_cache<std::pair<int, int>, std::vector<table_type> > wdt;
        bfem::detail::table_cache<std::pair<int, int>, std::vector<table_type> > wm;
        bfem::detail::table_cache<std::pair<int, int>, table_type> mstiff;
        bfem::detail::table_cache<std::pair<int, int>, table_type> mdtmass;
    };
    typedef bfem::detail::table_cache_state<maps> state;
    static state& st() { return bfem::detail::table_cache_instance<state>(); }

    static table_type conv(const rmat& R) {
        table_type out;
        out.zeros(R.rows, R.cols);
        for (int i = 0; i < R.rows; ++i)
            for (int j = 0; j < R.cols; ++j)
                out(i, j) = conv1(R.at(i, j));
        return out;
    }
    static table_type conv_mass(const bfem::mass_table<1>& M) {
        table_type out;
        out.zeros(M.rows(), M.cols());
        for (int i = 0; i < M.rows(); ++i)
            for (int j = 0; j < M.cols(); ++j)
                out(i, j) = conv1(M.at(i, j));
        return out;
    }
    static std::vector<table_type> conv_vec(const std::vector<rmat>& V) {
        std::vector<table_type> out;
        out.reserve(V.size());
        for (std::size_t k = 0; k < V.size(); ++k) out.push_back(conv(V[k]));
        return out;
    }
};

} // namespace detail
} // namespace btfem
} // namespace vcp

#endif // VCP_BTFEM_DETAIL_TIME_TABLES_HPP
