// vcp/bfem/rational.hpp
// Layer 0: exact rational arithmetic (detail::bigint / detail::rational)
// and the binomial coefficient cache (detail::binomial_cache).
//
// Conforms to: L0 external design v0.3 (section 3.3) and
//              L0 internal design v0.3 (sections 2, 3, 4).
//
// No floating point type appears anywhere in this header (requirement 12).
// Dependencies: C++11 standard library only.

#ifndef VCP_BFEM_RATIONAL_HPP
#define VCP_BFEM_RATIONAL_HPP

#include <vector>
#include <deque>
#include <string>
#include <utility>
#include <stdexcept>
#include <mutex>
#include <cstdint>
#include <cassert>

namespace vcp {
namespace bfem {
namespace detail {

// ---------------------------------------------------------------------------
// bigint: arbitrary precision signed integer.
//   representation: base 2^32 limbs, little endian, no leading zero limbs,
//   sign flag (zero is normalized to neg_ == false).
//   All operations give the strong exception guarantee (internal design 2.3).
// ---------------------------------------------------------------------------
class bigint {
public:
    bigint() : limb_(), neg_(false) {}

    bigint(long long v) : limb_(), neg_(false) {
        if (v < 0) {
            neg_ = true;
            // avoid overflow on LLONG_MIN: negate via unsigned
            std::uint64_t u = ~static_cast<std::uint64_t>(v) + 1u;
            push_u64(u);
        } else {
            push_u64(static_cast<std::uint64_t>(v));
        }
    }

    bool is_zero() const { return limb_.empty(); }
    bool negative() const { return neg_; }
    int  num_limbs() const { return static_cast<int>(limb_.size()); }
    std::uint32_t limb(int i) const {
        assert(i >= 0 && i < num_limbs());
        return limb_[static_cast<std::size_t>(i)];
    }
    bool is_even() const { return limb_.empty() || (limb_[0] & 1u) == 0u; }

    // true iff |x| fits in a signed 64-bit value (used by index_map's
    // bigint -> int table baking, internal design section 6)
    bool fits_int64() const {
        if (limb_.size() > 2) return false;
        if (limb_.size() == 2 && (limb_[1] & 0x80000000u) != 0u) return false;
        return true;
    }
    long long to_int64() const {
        if (!fits_int64()) throw std::overflow_error("bfem::bigint: to_int64 overflow");
        std::uint64_t u = 0;
        if (limb_.size() >= 1) u |= static_cast<std::uint64_t>(limb_[0]);
        if (limb_.size() == 2) u |= static_cast<std::uint64_t>(limb_[1]) << 32;
        long long r = static_cast<long long>(u);
        return neg_ ? -r : r;
    }

    // ---- comparison ----
    // returns -1 / 0 / +1
    static int cmp(const bigint& a, const bigint& b) {
        if (a.neg_ != b.neg_) return a.neg_ ? -1 : 1;
        int c = cmp_abs(a.limb_, b.limb_);
        return a.neg_ ? -c : c;
    }
    bool operator==(const bigint& o) const { return cmp(*this, o) == 0; }
    bool operator!=(const bigint& o) const { return cmp(*this, o) != 0; }
    bool operator<(const bigint& o) const  { return cmp(*this, o) < 0; }
    bool operator>(const bigint& o) const  { return cmp(*this, o) > 0; }
    bool operator<=(const bigint& o) const { return cmp(*this, o) <= 0; }
    bool operator>=(const bigint& o) const { return cmp(*this, o) >= 0; }

    // ---- arithmetic ----
    friend bigint operator+(const bigint& a, const bigint& b) {
        bigint r;
        if (a.neg_ == b.neg_) {
            r.limb_ = add_abs(a.limb_, b.limb_);
            r.neg_ = a.neg_;
        } else {
            int c = cmp_abs(a.limb_, b.limb_);
            if (c == 0) return r; // zero
            if (c > 0) { r.limb_ = sub_abs(a.limb_, b.limb_); r.neg_ = a.neg_; }
            else       { r.limb_ = sub_abs(b.limb_, a.limb_); r.neg_ = b.neg_; }
        }
        r.normalize_zero();
        return r;
    }
    friend bigint operator-(const bigint& a, const bigint& b) {
        bigint nb = b;
        if (!nb.is_zero()) nb.neg_ = !nb.neg_;
        return a + nb;
    }
    bigint operator-() const {
        bigint r = *this;
        if (!r.is_zero()) r.neg_ = !r.neg_;
        return r;
    }
    friend bigint operator*(const bigint& a, const bigint& b) {
        bigint r;
        if (a.is_zero() || b.is_zero()) return r;
        r.limb_ = mul_abs(a.limb_, b.limb_);
        r.neg_ = (a.neg_ != b.neg_);
        r.normalize_zero();
        return r;
    }

    // truncated division: q = trunc(a/b), r = a - q*b (sign of r = sign of a).
    // identity q*b + r == a is the I1 test invariant.
    static void divmod(const bigint& a, const bigint& b, bigint& q, bigint& r) {
        if (b.is_zero()) throw std::logic_error("bfem::bigint: division by zero");
        divmod_abs(a.limb_, b.limb_, q.limb_, r.limb_);
        q.neg_ = !q.limb_.empty() && (a.neg_ != b.neg_);
        r.neg_ = !r.limb_.empty() && a.neg_;
        q.normalize_zero();
        r.normalize_zero();
    }

    // gcd(|a|, |b|), binary algorithm (no divmod dependency).
    static bigint gcd(bigint a, bigint b) {
        a.neg_ = false;
        b.neg_ = false;
        if (a.is_zero()) return b;
        if (b.is_zero()) return a;
        int shift = 0;
        while (a.is_even() && b.is_even()) { a.shr1(); b.shr1(); ++shift; }
        while (a.is_even()) a.shr1();
        // invariant: a is odd
        while (!b.is_zero()) {
            while (b.is_even()) b.shr1();
            if (cmp_abs(a.limb_, b.limb_) > 0) a.limb_.swap(b.limb_);
            b.limb_ = sub_abs(b.limb_, a.limb_);
            bigint::trim(b.limb_);
        }
        while (shift-- > 0) a.shl1();
        return a;
    }

    // ---- string I/O (debug / test vectors; base 10^9 chunks) ----
    std::string to_string() const {
        if (is_zero()) return "0";
        std::vector<std::uint32_t> work = limb_;
        std::string digits;
        while (!work.empty()) {
            std::uint64_t rem = 0;
            for (std::size_t i = work.size(); i-- > 0;) {
                std::uint64_t cur = (rem << 32) | work[i];
                work[i] = static_cast<std::uint32_t>(cur / 1000000000u);
                rem = cur % 1000000000u;
            }
            while (!work.empty() && work.back() == 0) work.pop_back();
            for (int d = 0; d < 9; ++d) {
                digits.push_back(static_cast<char>('0' + rem % 10));
                rem /= 10;
            }
        }
        while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
        if (neg_) digits.push_back('-');
        return std::string(digits.rbegin(), digits.rend());
    }

    static bigint from_string(const std::string& s) {
        bigint r;
        std::size_t pos = 0;
        bool neg = false;
        if (pos < s.size() && (s[pos] == '-' || s[pos] == '+')) {
            neg = (s[pos] == '-');
            ++pos;
        }
        if (pos >= s.size()) throw std::invalid_argument("bfem::bigint: empty numeral");
        const bigint chunk_base(1000000000LL);
        // process leading partial chunk, then 9-digit chunks
        std::size_t ndig = s.size() - pos;
        std::size_t first = ndig % 9;
        if (first == 0) first = 9;
        std::size_t end = pos + first;
        while (pos < s.size()) {
            long long v = 0;
            for (; pos < end; ++pos) {
                if (s[pos] < '0' || s[pos] > '9')
                    throw std::invalid_argument("bfem::bigint: bad digit");
                v = v * 10 + (s[pos] - '0');
            }
            r = r * chunk_base + bigint(v);
            end = pos + 9;
        }
        if (neg && !r.is_zero()) r.neg_ = true;
        return r;
    }

    // left shift by one bit (used by tests and gcd)
    void shl1() {
        std::uint32_t carry = 0;
        for (std::size_t i = 0; i < limb_.size(); ++i) {
            std::uint32_t nc = limb_[i] >> 31;
            limb_[i] = (limb_[i] << 1) | carry;
            carry = nc;
        }
        if (carry) limb_.push_back(carry);
    }
    void shr1() {
        std::uint32_t carry = 0;
        for (std::size_t i = limb_.size(); i-- > 0;) {
            std::uint32_t nc = limb_[i] & 1u;
            limb_[i] = (limb_[i] >> 1) | (carry << 31);
            carry = nc;
        }
        trim(limb_);
        normalize_zero();
    }

private:
    std::vector<std::uint32_t> limb_;
    bool neg_;

    void push_u64(std::uint64_t u) {
        if (u == 0) return;
        limb_.push_back(static_cast<std::uint32_t>(u & 0xFFFFFFFFu));
        std::uint32_t hi = static_cast<std::uint32_t>(u >> 32);
        if (hi) limb_.push_back(hi);
    }
    void normalize_zero() {
        if (limb_.empty()) neg_ = false;
    }
    static void trim(std::vector<std::uint32_t>& v) {
        while (!v.empty() && v.back() == 0) v.pop_back();
    }
    static int cmp_abs(const std::vector<std::uint32_t>& a,
                       const std::vector<std::uint32_t>& b) {
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (std::size_t i = a.size(); i-- > 0;) {
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        }
        return 0;
    }
    static std::vector<std::uint32_t> add_abs(const std::vector<std::uint32_t>& a,
                                              const std::vector<std::uint32_t>& b) {
        const std::vector<std::uint32_t>& x = a.size() >= b.size() ? a : b;
        const std::vector<std::uint32_t>& y = a.size() >= b.size() ? b : a;
        std::vector<std::uint32_t> r;
        r.reserve(x.size() + 1);
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            std::uint64_t s = carry + x[i] + (i < y.size() ? y[i] : 0u);
            r.push_back(static_cast<std::uint32_t>(s & 0xFFFFFFFFu));
            carry = s >> 32;
        }
        if (carry) r.push_back(static_cast<std::uint32_t>(carry));
        return r;
    }
    // precondition: |a| >= |b|
    static std::vector<std::uint32_t> sub_abs(const std::vector<std::uint32_t>& a,
                                              const std::vector<std::uint32_t>& b) {
        std::vector<std::uint32_t> r;
        r.reserve(a.size());
        std::int64_t borrow = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            std::int64_t s = static_cast<std::int64_t>(a[i]) - borrow
                             - (i < b.size() ? static_cast<std::int64_t>(b[i]) : 0);
            if (s < 0) { s += (static_cast<std::int64_t>(1) << 32); borrow = 1; }
            else borrow = 0;
            r.push_back(static_cast<std::uint32_t>(s));
        }
        assert(borrow == 0);
        trim(r);
        return r;
    }
    static std::vector<std::uint32_t> mul_abs(const std::vector<std::uint32_t>& a,
                                              const std::vector<std::uint32_t>& b) {
        std::vector<std::uint32_t> r(a.size() + b.size(), 0u);
        for (std::size_t i = 0; i < a.size(); ++i) {
            std::uint64_t carry = 0;
            std::uint64_t ai = a[i];
            for (std::size_t j = 0; j < b.size(); ++j) {
                std::uint64_t cur = static_cast<std::uint64_t>(r[i + j])
                                    + ai * b[j] + carry;
                r[i + j] = static_cast<std::uint32_t>(cur & 0xFFFFFFFFu);
                carry = cur >> 32;
            }
            std::size_t k = i + b.size();
            while (carry) {
                std::uint64_t cur = static_cast<std::uint64_t>(r[k]) + carry;
                r[k] = static_cast<std::uint32_t>(cur & 0xFFFFFFFFu);
                carry = cur >> 32;
                ++k;
            }
        }
        trim(r);
        return r;
    }
    // binary long division on magnitudes: a = q*b + r, 0 <= r < b
    static void divmod_abs(const std::vector<std::uint32_t>& a,
                           const std::vector<std::uint32_t>& b,
                           std::vector<std::uint32_t>& q,
                           std::vector<std::uint32_t>& r) {
        assert(!b.empty());
        std::vector<std::uint32_t> quot(a.size(), 0u);
        std::vector<std::uint32_t> rem;
        int total_bits = static_cast<int>(a.size()) * 32;
        for (int bit = total_bits - 1; bit >= 0; --bit) {
            // rem = rem*2 + bit_of_a
            std::uint32_t carry = 0;
            for (std::size_t i = 0; i < rem.size(); ++i) {
                std::uint32_t nc = rem[i] >> 31;
                rem[i] = (rem[i] << 1) | carry;
                carry = nc;
            }
            if (carry) rem.push_back(carry);
            std::uint32_t abit = (a[static_cast<std::size_t>(bit / 32)] >> (bit % 32)) & 1u;
            if (abit) {
                if (rem.empty()) rem.push_back(1u);
                else rem[0] |= 1u;
            }
            if (cmp_abs(rem, b) >= 0) {
                rem = sub_abs(rem, b);
                quot[static_cast<std::size_t>(bit / 32)] |=
                    (static_cast<std::uint32_t>(1) << (bit % 32));
            }
        }
        trim(quot);
        trim(rem);
        q = quot;
        r = rem;
    }
};

// ---------------------------------------------------------------------------
// rational: exact rational number.
//   invariants (internal design 3.1): den > 0, gcd(|num|, den) == 1, zero is 0/1.
//   Publicly the L0 external design exposes read-only observers; the full
//   arithmetic lives here in detail so that rational satisfies the T type
//   requirements of L1 (external design section 10 + total order, C-4) and can
//   be used as the exact test scalar (bpoly<D, rational> instantiation).
// ---------------------------------------------------------------------------
class rational {
public:
    rational() : num_(0), den_(1) {}
    rational(int v) : num_(v), den_(1) {}
    rational(long long v) : num_(v), den_(1) {}
    rational(long long num, long long den) : num_(num), den_(den) { normalize(); }
    rational(bigint num, bigint den)
        : num_(std::move(num)), den_(std::move(den)) { normalize(); }

    // ---- observers (the public L0 API surface) ----
    const bigint& num() const { return num_; }   // normalized: gcd == 1
    const bigint& den() const { return den_; }   // normalized: den > 0
    bool is_zero() const { return num_.is_zero(); }
    std::string to_string() const { return num_.to_string() + "/" + den_.to_string(); }

    // ---- arithmetic (detail; T requirements of L1 section 10 + C-4) ----
    rational& operator+=(const rational& o) {
        num_ = num_ * o.den_ + o.num_ * den_;
        den_ = den_ * o.den_;
        normalize();
        return *this;
    }
    rational& operator-=(const rational& o) {
        num_ = num_ * o.den_ - o.num_ * den_;
        den_ = den_ * o.den_;
        normalize();
        return *this;
    }
    rational& operator*=(const rational& o) {
        num_ = num_ * o.num_;
        den_ = den_ * o.den_;
        normalize();
        return *this;
    }
    rational& operator/=(const rational& o) {
        if (o.num_.is_zero())
            throw std::logic_error("bfem::rational: division by zero");
        num_ = num_ * o.den_;
        den_ = den_ * o.num_;
        normalize();
        return *this;
    }
    friend rational operator+(rational a, const rational& b) { a += b; return a; }
    friend rational operator-(rational a, const rational& b) { a -= b; return a; }
    friend rational operator*(rational a, const rational& b) { a *= b; return a; }
    friend rational operator/(rational a, const rational& b) { a /= b; return a; }
    rational operator-() const {
        rational r = *this;
        r.num_ = -r.num_;
        return r;
    }

    // total order (den > 0 always, so cross multiplication preserves order)
    friend bool operator==(const rational& a, const rational& b) {
        return a.num_ == b.num_ && a.den_ == b.den_;
    }
    friend bool operator!=(const rational& a, const rational& b) { return !(a == b); }
    friend bool operator<(const rational& a, const rational& b) {
        return bigint::cmp(a.num_ * b.den_, b.num_ * a.den_) < 0;
    }
    friend bool operator>(const rational& a, const rational& b)  { return b < a; }
    friend bool operator<=(const rational& a, const rational& b) { return !(b < a); }
    friend bool operator>=(const rational& a, const rational& b) { return !(a < b); }

private:
    bigint num_;
    bigint den_;

    void normalize() {
        if (den_.is_zero())
            throw std::logic_error("bfem::rational: zero denominator");
        if (den_.negative()) { den_ = -den_; num_ = -num_; }
        if (num_.is_zero()) { den_ = bigint(1); return; }
        bigint g = bigint::gcd(num_, den_);
        if (!(g == bigint(1))) {
            bigint q, r;
            bigint::divmod(num_, g, q, r);
            assert(r.is_zero());
            num_ = q;
            bigint::divmod(den_, g, q, r);
            assert(r.is_zero());
            den_ = q;
        }
    }
};

// ---------------------------------------------------------------------------
// binomial_cache: Pascal triangle of bigint, grown row by row on demand.
//   Storage: deque<vector<bigint>> (row append never invalidates references),
//   growth protected by a dedicated mutex (C-1).
//   Access is restricted by design to table generation (under the registry
//   lock) and index_map construction; runtime rank/unrank paths never touch it.
// ---------------------------------------------------------------------------
class binomial_cache {
public:
    // C(n, k). precondition: n >= 0, 0 <= k <= n.
    static const bigint& binom(int n, int k) {
        if (n < 0 || k < 0 || k > n)
            throw std::invalid_argument("bfem::binomial_cache: bad (n, k)");
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        while (static_cast<int>(s.rows.size()) <= n) {
            int r = static_cast<int>(s.rows.size());
            std::vector<bigint> row(static_cast<std::size_t>(r) + 1);
            row[0] = bigint(1);
            row[static_cast<std::size_t>(r)] = bigint(1);
            for (int j = 1; j < r; ++j) {
                row[static_cast<std::size_t>(j)] =
                    s.rows[static_cast<std::size_t>(r - 1)][static_cast<std::size_t>(j - 1)]
                    + s.rows[static_cast<std::size_t>(r - 1)][static_cast<std::size_t>(j)];
            }
            s.rows.push_back(std::move(row));
        }
        return s.rows[static_cast<std::size_t>(n)][static_cast<std::size_t>(k)];
    }

    // number of rows currently cached (test hook, I3)
    static int cached_rows() {
        state& s = st();
        std::lock_guard<std::mutex> lk(s.mtx);
        return static_cast<int>(s.rows.size());
    }

private:
    struct state {
        std::mutex mtx;
        std::deque<std::vector<bigint> > rows;
    };
    static state& st() {
        static state s;
        return s;
    }
};

// multinomial coefficient C(n, alpha) = n! / prod alpha_i!, constructed as a
// product of binomials (internal design section 4; factorials are never built):
//   C(n, alpha) = prod_{i=1..D} C(alpha_0 + ... + alpha_i, alpha_i)
// alpha is passed as a pointer range to stay independent of multi_index<D>.
inline bigint multinomial(int n, const int* alpha, int count) {
    bigint r(1);
    int partial = (count > 0) ? alpha[0] : 0;
    for (int i = 1; i < count; ++i) {
        partial += alpha[i];
        r = r * binomial_cache::binom(partial, alpha[i]);
    }
    assert(partial == n);
    (void)n;
    return r;
}

} // namespace detail

// Public exposure of the exact rational type (read-only usage is the
// documented contract for upper layers; see gate report for the recorded
// supplementation regarding operator visibility).
typedef detail::rational rational;

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RATIONAL_HPP
