// vcp/bfem/sv/linear_reduction.hpp
// Phase 5d (Scott-Vogelius parts): linear_reduction<T,P,SP> (V3, SV-3/SV-8)
// -- division-free congruence reduction by an independent set of integer
// constraint rows (the generalization of dirichlet_reduction).
//
// Conforms to: SV external design v0.2 (section 5) and
//              SV internal design v0.2 (section 4).
//
// Construction (internal design section 4):
//  1. every row is validated (index range, nonzero integer coefficients,
//     no duplicate dof inside a row);
//  2. integer Gauss-Jordan with +-1 pivots brings the rows to the normal
//     form "each pivot appears in no other row". The pivot of a row is the
//     LARGEST dof index among its +-1 coefficients (deterministic). Row
//     updates are integer multiply-subtract only (the pivot is +-1, so no
//     division exists anywhere in this file -- G-SV-1);
//  3. a row that becomes zero is DEPENDENT, a row without any +-1 entry has
//     no admissible pivot: both raise std::invalid_argument;
//  4. the embedding T (full x reduced, int entries) is stored row-wise:
//     kept dof i        -> the single pair (reduced(i), 1)
//     pivot dof p (row r with pivot coefficient s = +-1, entries c_j on the
//     kept dofs j)      -> pairs (reduced(j), -s * c_j)
//     so x = T x_red satisfies every constraint EXACTLY (expand roundtrip).
//
// reduce(A) = T^t A T, reduce_rows(A) = T^t A, reduce_cols(A) = A T
// (v0.2, A-2), reduce(b) = T^t b, expand(x) = T x. All matrix work is COO
// enumeration through the frozen spm_adapter + coo_buffer (deterministic
// combination); the only T operations are multiplications by the integer
// weights (skipped when the weight is one) and the additions of combine.
// Acceptance identity (T-SV-4): reduce(A) == reduce_rows(reduce_cols(A))
// exactly for exact T.
//
// RATIONAL-COEFFICIENT EXTENSION (phase 6, C1-6 -- ADDITIVE ONLY, the
// integer path above is textually and behaviourally unchanged):
//  - sv_constraint_q carries exact rational coefficients (the C1 corner
//    constraints t . grad u = 0 / t^T H t = 0 mix dofs with coordinate-
//    difference coefficients);
//  - the overloaded constructor eliminates the rows EXACTLY at the rational
//    stage (Gauss-Jordan on detail::rational; pivot = the first nonzero in
//    canonical ascending-dof order -- exactness holds for any nonzero pivot,
//    determinism is the criterion; rational division stays confined to this
//    generation stage);
//  - the embedding weights are converted to T exactly ONCE through
//    convert_traits (enclose-once; verified bitwise by CI-6): afterwards
//    reduce/reduce_rows/reduce_cols/expand are multiplications only (weight
//    one is skipped as in the integer path; T division count stays ZERO);
//  - the integer +-1 path remains the special case served by the original
//    constructor (same public API; SV regressions stay bit-invariant).

#ifndef VCP_BFEM_SV_LINEAR_REDUCTION_HPP
#define VCP_BFEM_SV_LINEAR_REDUCTION_HPP

#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/fe_space.hpp>          // detail::coo_buffer / spm_adapter
#include <vcp/bfem/rational.hpp>          // detail::rational (q extension)
#include <vcp/bfem/convert_traits.hpp>    // enclose-once weights (q extension)
#include <vcp/bfem/sv/sv_constraint.hpp>

namespace vcp {
namespace bfem {

// one exact rational constraint row: sum_k coef[k] * x_{dof[k]} = 0
// (phase 6 additive extension; the integer sv_constraint stays untouched)
struct sv_constraint_q {
    std::vector<int> dof;
    std::vector<detail::rational> coef;
};

template <typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class linear_reduction {
public:
    typedef vcp::spmatrix<T, SP> spmatrix_t;

    linear_reduction(int full_size, std::vector<sv_constraint> rows)
        : full_(full_size), to_reduced_(), to_full_(), weights_(),
          qmode_(false), qweights_() {
        if (full_size < 0)
            throw std::invalid_argument("bfem::linear_reduction: negative size");

        // working rows as sorted maps dof -> integer coefficient
        typedef std::map<int, long long> irow;
        std::vector<irow> w;
        w.reserve(rows.size());
        for (std::size_t k = 0; k < rows.size(); ++k) {
            const sv_constraint& r = rows[k];
            if (r.dof.size() != r.coef.size() || r.dof.empty())
                throw std::invalid_argument(
                    "bfem::linear_reduction: malformed constraint row");
            irow m;
            for (std::size_t s = 0; s < r.dof.size(); ++s) {
                if (r.dof[s] < 0 || r.dof[s] >= full_size)
                    throw std::invalid_argument(
                        "bfem::linear_reduction: dof out of range");
                if (r.coef[s] == 0)
                    throw std::invalid_argument(
                        "bfem::linear_reduction: zero coefficient");
                if (!m.insert(std::make_pair(r.dof[s],
                        static_cast<long long>(r.coef[s]))).second)
                    throw std::invalid_argument(
                        "bfem::linear_reduction: duplicate dof in a row");
            }
            w.push_back(m);
        }

        // integer Gauss-Jordan with +-1 pivots (largest-index rule)
        std::vector<int> piv(w.size(), -1);
        std::vector<long long> pivc(w.size(), 0);
        for (std::size_t k = 0; k < w.size(); ++k) {
            // eliminate the already-chosen pivots from row k
            for (std::size_t j = 0; j < k; ++j) {
                irow::iterator it = w[k].find(piv[j]);
                if (it == w[k].end()) continue;
                long long f = it->second * pivc[j];    // c * s (s * s == 1)
                for (irow::const_iterator p = w[j].begin(); p != w[j].end(); ++p) {
                    long long& e = w[k][p->first];
                    e -= f * p->second;
                    if (e == 0) w[k].erase(p->first);
                }
            }
            if (w[k].empty())
                throw std::invalid_argument(
                    "bfem::linear_reduction: dependent constraint row");
            // pivot: largest dof index among the +-1 coefficients
            int pk = -1;
            long long sk = 0;
            for (irow::const_iterator p = w[k].begin(); p != w[k].end(); ++p) {
                if (p->second == 1 || p->second == -1) { pk = p->first; sk = p->second; }
            }
            if (pk < 0)
                throw std::invalid_argument(
                    "bfem::linear_reduction: no +-1 pivot in a constraint row");
            piv[static_cast<std::size_t>(k)] = pk;
            pivc[static_cast<std::size_t>(k)] = sk;
            // eliminate the new pivot from the earlier rows (Jordan step)
            for (std::size_t j = 0; j < k; ++j) {
                irow::iterator it = w[j].find(pk);
                if (it == w[j].end()) continue;
                long long f = it->second * sk;
                for (irow::const_iterator p = w[k].begin(); p != w[k].end(); ++p) {
                    long long& e = w[j][p->first];
                    e -= f * p->second;
                    if (e == 0) w[j].erase(p->first);
                }
            }
        }

        // reduced numbering over the kept (non-pivot) dofs
        to_reduced_.assign(static_cast<std::size_t>(full_size), 0);
        for (std::size_t k = 0; k < w.size(); ++k)
            to_reduced_[static_cast<std::size_t>(piv[k])] = -1;
        int r = 0;
        to_full_.reserve(static_cast<std::size_t>(full_size) - w.size());
        for (int i = 0; i < full_size; ++i) {
            if (to_reduced_[static_cast<std::size_t>(i)] == -1) continue;
            to_reduced_[static_cast<std::size_t>(i)] = r;
            to_full_.push_back(i);
            ++r;
        }

        // embedding rows of T
        weights_.assign(static_cast<std::size_t>(full_size),
                        std::vector<std::pair<int, long long> >());
        for (int i = 0; i < full_size; ++i) {
            int ri = to_reduced_[static_cast<std::size_t>(i)];
            if (ri >= 0)
                weights_[static_cast<std::size_t>(i)].push_back(
                    std::make_pair(ri, 1LL));
        }
        for (std::size_t k = 0; k < w.size(); ++k) {
            std::vector<std::pair<int, long long> >& row =
                weights_[static_cast<std::size_t>(piv[k])];
            for (irow::const_iterator p = w[k].begin(); p != w[k].end(); ++p) {
                if (p->first == piv[k]) continue;
                // x_p = -s * sum c_j x_j over kept dofs j (normal form)
                row.push_back(std::make_pair(
                    to_reduced_[static_cast<std::size_t>(p->first)],
                    -pivc[k] * p->second));
            }
        }
    }

    // ---- rational-coefficient overload (phase 6 additive extension) ----
    // Exact rational Gauss-Jordan normal form; pivot = the FIRST nonzero in
    // canonical ascending-dof order (deterministic; any nonzero is exact).
    // The embedding weights pass through convert_traits exactly once.
    linear_reduction(int full_size, const std::vector<sv_constraint_q>& rows)
        : full_(full_size), to_reduced_(), to_full_(), weights_(),
          qmode_(true), qweights_() {
        if (full_size < 0)
            throw std::invalid_argument("bfem::linear_reduction: negative size");
        typedef detail::rational rat;
        typedef std::map<int, rat> qrow;
        std::vector<qrow> w;
        w.reserve(rows.size());
        for (std::size_t k = 0; k < rows.size(); ++k) {
            const sv_constraint_q& r = rows[k];
            if (r.dof.size() != r.coef.size() || r.dof.empty())
                throw std::invalid_argument(
                    "bfem::linear_reduction: malformed constraint row");
            qrow m;
            for (std::size_t s = 0; s < r.dof.size(); ++s) {
                if (r.dof[s] < 0 || r.dof[s] >= full_size)
                    throw std::invalid_argument(
                        "bfem::linear_reduction: dof out of range");
                if (r.coef[s].is_zero())
                    throw std::invalid_argument(
                        "bfem::linear_reduction: zero coefficient");
                if (!m.insert(std::make_pair(r.dof[s], r.coef[s])).second)
                    throw std::invalid_argument(
                        "bfem::linear_reduction: duplicate dof in a row");
            }
            w.push_back(m);
        }

        // exact Gauss-Jordan with canonical-first pivots (rational stage;
        // the ONLY divisions of the extension happen here, on rationals)
        std::vector<int> piv(w.size(), -1);
        std::vector<rat> pivc(w.size());
        for (std::size_t k = 0; k < w.size(); ++k) {
            for (std::size_t j = 0; j < k; ++j) {
                typename qrow::iterator it = w[k].find(piv[j]);
                if (it == w[k].end()) continue;
                rat f = it->second / pivc[j];
                for (typename qrow::const_iterator p = w[j].begin();
                     p != w[j].end(); ++p) {
                    if (p->first == piv[j]) {
                        w[k].erase(p->first);
                        continue;
                    }
                    rat& e = w[k][p->first];
                    e -= f * p->second;
                    if (e.is_zero()) w[k].erase(p->first);
                }
            }
            if (w[k].empty())
                throw std::invalid_argument(
                    "bfem::linear_reduction: dependent constraint row");
            // pivot: first nonzero in canonical (ascending dof) order
            piv[static_cast<std::size_t>(k)] = w[k].begin()->first;
            pivc[static_cast<std::size_t>(k)] = w[k].begin()->second;
            for (std::size_t j = 0; j < k; ++j) {
                typename qrow::iterator it = w[j].find(piv[k]);
                if (it == w[j].end()) continue;
                rat f = it->second / pivc[k];
                for (typename qrow::const_iterator p = w[k].begin();
                     p != w[k].end(); ++p) {
                    if (p->first == piv[k]) {
                        w[j].erase(p->first);
                        continue;
                    }
                    rat& e = w[j][p->first];
                    e -= f * p->second;
                    if (e.is_zero()) w[j].erase(p->first);
                }
            }
        }

        // reduced numbering over the kept dofs (identical bookkeeping)
        to_reduced_.assign(static_cast<std::size_t>(full_size), 0);
        for (std::size_t k = 0; k < w.size(); ++k)
            to_reduced_[static_cast<std::size_t>(piv[k])] = -1;
        int r = 0;
        to_full_.reserve(static_cast<std::size_t>(full_size) - w.size());
        for (int i = 0; i < full_size; ++i) {
            if (to_reduced_[static_cast<std::size_t>(i)] == -1) continue;
            to_reduced_[static_cast<std::size_t>(i)] = r;
            to_full_.push_back(i);
            ++r;
        }

        // embedding rows: kept dof -> weight one; pivot dof p of row k ->
        // x_p = -(1/c_p) sum_{j != p} c_j x_j, each weight converted ONCE
        qweights_.assign(static_cast<std::size_t>(full_size),
                         std::vector<qentry>());
        for (int i = 0; i < full_size; ++i) {
            int ri = to_reduced_[static_cast<std::size_t>(i)];
            if (ri >= 0) {
                qentry e;
                e.col = ri;
                e.w = T(1);
                e.one = true;
                qweights_[static_cast<std::size_t>(i)].push_back(e);
            }
        }
        for (std::size_t k = 0; k < w.size(); ++k) {
            std::vector<qentry>& row = qweights_[static_cast<std::size_t>(piv[k])];
            for (typename qrow::const_iterator p = w[k].begin(); p != w[k].end(); ++p) {
                if (p->first == piv[k]) continue;
                rat wq = -(p->second / pivc[k]);
                qentry e;
                e.col = to_reduced_[static_cast<std::size_t>(p->first)];
                e.w = convert_traits<T>::from_rational(wq.num(), wq.den());
                e.one = (wq == rat(1));
                row.push_back(e);
            }
        }
    }

    int full_size() const { return full_; }
    int reduced_size() const { return static_cast<int>(to_full_.size()); }
    int num_constraints() const { return full_ - reduced_size(); }
    int to_reduced(int full_dof) const {          // pivot dofs map to -1
        if (full_dof < 0 || full_dof >= full_)
            throw std::invalid_argument(
                "bfem::linear_reduction::to_reduced: out of range");
        return to_reduced_[static_cast<std::size_t>(full_dof)];
    }
    int to_full(int reduced_dof) const {
        if (reduced_dof < 0 || reduced_dof >= reduced_size())
            throw std::invalid_argument(
                "bfem::linear_reduction::to_full: out of range");
        return to_full_[static_cast<std::size_t>(reduced_dof)];
    }

    // ---- congruence reduction T^t A T ----
    spmatrix_t reduce(const spmatrix_t& A) const {
        if (A.rowsize() != full_ || A.columnsize() != full_)
            throw std::invalid_argument("bfem::linear_reduction::reduce: size mismatch");
        detail::coo_buffer<T> buf;
        if (qmode_) {
            q_both_collector col = { this, &buf };
            detail::spm_adapter<T, SP>::for_each_entry(A, col);
        } else {
            both_collector col = { this, &buf };
            detail::spm_adapter<T, SP>::for_each_entry(A, col);
        }
        buf.combine();
        return detail::spm_adapter<T, SP>::build(reduced_size(), reduced_size(), buf);
    }

    // ---- one-sided reductions (v0.2, A-2) ----
    spmatrix_t reduce_rows(const spmatrix_t& A) const {
        if (A.rowsize() != full_)
            throw std::invalid_argument(
                "bfem::linear_reduction::reduce_rows: row size mismatch");
        detail::coo_buffer<T> buf;
        if (qmode_) {
            q_row_collector col = { this, &buf };
            detail::spm_adapter<T, SP>::for_each_entry(A, col);
        } else {
            row_collector col = { this, &buf };
            detail::spm_adapter<T, SP>::for_each_entry(A, col);
        }
        buf.combine();
        return detail::spm_adapter<T, SP>::build(reduced_size(), A.columnsize(), buf);
    }
    spmatrix_t reduce_cols(const spmatrix_t& A) const {
        if (A.columnsize() != full_)
            throw std::invalid_argument(
                "bfem::linear_reduction::reduce_cols: column size mismatch");
        detail::coo_buffer<T> buf;
        if (qmode_) {
            q_col_collector col = { this, &buf };
            detail::spm_adapter<T, SP>::for_each_entry(A, col);
        } else {
            col_collector col = { this, &buf };
            detail::spm_adapter<T, SP>::for_each_entry(A, col);
        }
        buf.combine();
        return detail::spm_adapter<T, SP>::build(A.rowsize(), reduced_size(), buf);
    }

    // ---- vectors ----
    vcp::matrix<T, P> reduce(const vcp::matrix<T, P>& v) const {   // T^t v
        if (v.rowsize() != full_ || v.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::linear_reduction::reduce: vector size mismatch");
        vcp::matrix<T, P> out;
        out.zeros(reduced_size(), 1);
        if (qmode_) {
            for (int i = 0; i < full_; ++i) {
                const std::vector<qentry>& wr =
                    qweights_[static_cast<std::size_t>(i)];
                for (std::size_t s = 0; s < wr.size(); ++s) {
                    if (wr[s].one) out(wr[s].col, 0) += v(i, 0);
                    else out(wr[s].col, 0) += v(i, 0) * wr[s].w;
                }
            }
            return out;
        }
        for (int i = 0; i < full_; ++i) {
            const std::vector<std::pair<int, long long> >& wr =
                weights_[static_cast<std::size_t>(i)];
            for (std::size_t s = 0; s < wr.size(); ++s) {
                if (wr[s].second == 1LL) out(wr[s].first, 0) += v(i, 0);
                else out(wr[s].first, 0) += v(i, 0) * wfac(wr[s].second);
            }
        }
        return out;
    }
    vcp::matrix<T, P> expand(const vcp::matrix<T, P>& v_reduced) const {  // T x
        if (v_reduced.rowsize() != reduced_size() || v_reduced.columnsize() != 1)
            throw std::invalid_argument(
                "bfem::linear_reduction::expand: vector size mismatch");
        vcp::matrix<T, P> out;
        out.zeros(full_, 1);
        if (qmode_) {
            for (int i = 0; i < full_; ++i) {
                const std::vector<qentry>& wr =
                    qweights_[static_cast<std::size_t>(i)];
                for (std::size_t s = 0; s < wr.size(); ++s) {
                    if (wr[s].one) out(i, 0) += v_reduced(wr[s].col, 0);
                    else out(i, 0) += v_reduced(wr[s].col, 0) * wr[s].w;
                }
            }
            return out;
        }
        for (int i = 0; i < full_; ++i) {
            const std::vector<std::pair<int, long long> >& wr =
                weights_[static_cast<std::size_t>(i)];
            for (std::size_t s = 0; s < wr.size(); ++s) {
                if (wr[s].second == 1LL) out(i, 0) += v_reduced(wr[s].first, 0);
                else out(i, 0) += v_reduced(wr[s].first, 0) * wfac(wr[s].second);
            }
        }
        return out;
    }

private:
    int full_;
    std::vector<int> to_reduced_;
    std::vector<int> to_full_;
    // row i of the embedding T: pairs (reduced column, integer weight)
    std::vector<std::vector<std::pair<int, long long> > > weights_;
    // rational-mode embedding (phase 6 additive extension): weights are the
    // enclose-once T images of the exact rational elimination
    struct qentry {
        int col;
        T w;
        bool one;                    // exact rational weight == 1: skip mult
    };
    bool qmode_;
    std::vector<std::vector<qentry> > qweights_;

    static T wfac(long long w) { return T(static_cast<int>(w)); }

    struct both_collector {
        const linear_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            const std::vector<std::pair<int, long long> >& wi =
                self->weights_[static_cast<std::size_t>(i)];
            const std::vector<std::pair<int, long long> >& wj =
                self->weights_[static_cast<std::size_t>(j)];
            for (std::size_t a = 0; a < wi.size(); ++a) {
                for (std::size_t b = 0; b < wj.size(); ++b) {
                    long long w = wi[a].second * wj[b].second;
                    if (w == 1LL) buf->push(wi[a].first, wj[b].first, v);
                    else buf->push(wi[a].first, wj[b].first, v * wfac(w));
                }
            }
        }
    };
    struct row_collector {
        const linear_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            const std::vector<std::pair<int, long long> >& wi =
                self->weights_[static_cast<std::size_t>(i)];
            for (std::size_t a = 0; a < wi.size(); ++a) {
                if (wi[a].second == 1LL) buf->push(wi[a].first, j, v);
                else buf->push(wi[a].first, j, v * wfac(wi[a].second));
            }
        }
    };
    struct col_collector {
        const linear_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            const std::vector<std::pair<int, long long> >& wj =
                self->weights_[static_cast<std::size_t>(j)];
            for (std::size_t b = 0; b < wj.size(); ++b) {
                if (wj[b].second == 1LL) buf->push(i, wj[b].first, v);
                else buf->push(i, wj[b].first, v * wfac(wj[b].second));
            }
        }
    };

    // ---- rational-mode collectors (phase 6 additive extension) ----
    struct q_both_collector {
        const linear_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            const std::vector<qentry>& wi =
                self->qweights_[static_cast<std::size_t>(i)];
            const std::vector<qentry>& wj =
                self->qweights_[static_cast<std::size_t>(j)];
            for (std::size_t a = 0; a < wi.size(); ++a) {
                for (std::size_t b = 0; b < wj.size(); ++b) {
                    if (wi[a].one && wj[b].one)
                        buf->push(wi[a].col, wj[b].col, v);
                    else if (wi[a].one)
                        buf->push(wi[a].col, wj[b].col, v * wj[b].w);
                    else if (wj[b].one)
                        buf->push(wi[a].col, wj[b].col, v * wi[a].w);
                    else
                        buf->push(wi[a].col, wj[b].col, v * wi[a].w * wj[b].w);
                }
            }
        }
    };
    struct q_row_collector {
        const linear_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            const std::vector<qentry>& wi =
                self->qweights_[static_cast<std::size_t>(i)];
            for (std::size_t a = 0; a < wi.size(); ++a) {
                if (wi[a].one) buf->push(wi[a].col, j, v);
                else buf->push(wi[a].col, j, v * wi[a].w);
            }
        }
    };
    struct q_col_collector {
        const linear_reduction* self;
        detail::coo_buffer<T>* buf;
        void operator()(int i, int j, const T& v) const {
            const std::vector<qentry>& wj =
                self->qweights_[static_cast<std::size_t>(j)];
            for (std::size_t b = 0; b < wj.size(); ++b) {
                if (wj[b].one) buf->push(i, wj[b].col, v);
                else buf->push(i, wj[b].col, v * wj[b].w);
            }
        }
    };
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_SV_LINEAR_REDUCTION_HPP
