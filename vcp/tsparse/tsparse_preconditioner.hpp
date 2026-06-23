// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_PRECONDITIONER_HPP
#define VCP_TSPARSE_PRECONDITIONER_HPP

#include <cstddef>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_factorization.hpp>

namespace vcp {

// -----------------------------------------------------------------------
// identity_preconditioner<T>
// -----------------------------------------------------------------------
template <typename T>
struct identity_preconditioner {
	typedef T value_type;

	void apply(const std::vector<T>& r, std::vector<T>& z) const {
		if (r.size() != z.size()) z.resize(r.size());
		z = r;
	}
};

// -----------------------------------------------------------------------
// preconditioner_operator<T> – abstract base class interface
// -----------------------------------------------------------------------
template <typename T>
class preconditioner_operator {
public:
	virtual ~preconditioner_operator() {}
	virtual void apply(const std::vector<T>& r, std::vector<T>& z) const = 0;
};

// -----------------------------------------------------------------------
// function_preconditioner<T, Functor>
// Wraps any callable f(r, z) as a preconditioner.
// -----------------------------------------------------------------------
template <typename T, class Functor>
class function_preconditioner {
public:
	typedef T value_type;

	explicit function_preconditioner(const Functor& f) : f_(f) {}

	void apply(const std::vector<T>& r, std::vector<T>& z) const {
		f_(r, z);
	}

private:
	Functor f_;
};

// -----------------------------------------------------------------------
// jacobi_preconditioner<SparseMatrix>
//
// Non-throwing construction.  Stores valid() / zero_pivots() / diagnostics().
// apply() throws vcp::numerical_error if !valid().
// apply() throws vcp::dimension_error on size mismatch.
//
// Near-zero diagonal detection threshold:  abs(d_i) <= pivot_tol
// where pivot_tol = decimal_power_negative<real_type>(14).
// -----------------------------------------------------------------------
template <class SparseMatrix>
class jacobi_preconditioner {
public:
	typedef typename SparseMatrix::value_type value_type;
	typedef typename SparseMatrix::index_type index_type;
	typedef typename vcp::tsparse_scalar::real_type<value_type>::type real_type;

	explicit jacobi_preconditioner(const SparseMatrix& K)
		: valid_(true), zero_pivots_(0), first_zero_pivot_(0)
	{
		const index_type n = K.rowsize();
		inv_diag_.assign(static_cast<std::size_t>(n), value_type(0));
		const real_type pivot_tol =
			vcp::tsparse_scalar::decimal_power_negative<real_type>(14);

		for (index_type i = 0; i < n; i++) {
			const value_type diag = K.get(i, i);
			const real_type abs_d = vcp::tsparse_scalar::abs_value(diag);
			if (abs_d <= pivot_tol) {
				if (zero_pivots_ == 0)
					first_zero_pivot_ = static_cast<std::size_t>(i);
				zero_pivots_++;
				valid_ = false;
				inv_diag_[static_cast<std::size_t>(i)] = value_type(0);
			} else {
				inv_diag_[static_cast<std::size_t>(i)] =
					value_type(1) / diag;
			}
		}
		build_diagnostics();
	}

	bool        valid()       const { return valid_; }
	std::size_t zero_pivots() const { return zero_pivots_; }

	const std::string& diagnostics() const { return diagnostics_; }

	void apply(const std::vector<value_type>& r, std::vector<value_type>& z) const {
		if (!valid_) {
			vcp::throw_error<vcp::numerical_error>(
				"jacobi_preconditioner::apply: invalid preconditioner; ",
				diagnostics_);
		}
		if (r.size() != inv_diag_.size()) {
			vcp::throw_error<vcp::dimension_error>(
				"jacobi_preconditioner::apply: dimension mismatch");
		}
		z.resize(r.size());
		for (std::size_t i = 0; i < r.size(); i++)
			z[i] = inv_diag_[i] * r[i];
	}

private:
	void build_diagnostics() {
		std::ostringstream os;
		os << "jacobi_preconditioner: zero_pivots=" << zero_pivots_;
		if (zero_pivots_ > 0)
			os << " first_zero_pivot=" << first_zero_pivot_;
		diagnostics_ = os.str();
	}

	bool        valid_;
	std::size_t zero_pivots_;
	std::size_t first_zero_pivot_;
	std::vector<value_type> inv_diag_;
	std::string diagnostics_;
};

// -----------------------------------------------------------------------
// ilu0_preconditioner<SparseMatrix>
//
// Builds ILU(0) factorization on construction (non-throwing).
// valid() / zero_pivots() / diagnostics() reflect factorization state.
// apply() throws vcp::numerical_error if factorization failed.
// -----------------------------------------------------------------------
template <class SparseMatrix>
class ilu0_preconditioner {
public:
	typedef typename SparseMatrix::value_type value_type;
	typedef typename SparseMatrix::index_type index_type;
	typedef typename vcp::tsparse_scalar::real_type<value_type>::type real_type;

	explicit ilu0_preconditioner(const SparseMatrix& K)
		: valid_(false), zero_pivots_(0)
	{
		SparseMatrix K_csr = K.as_csr();
		data_ = vcp::tsparse_factorization::ilu0_factorize<value_type, index_type>(
			K_csr.outer_index(),
			K_csr.inner_index(),
			K_csr.values(),
			static_cast<std::size_t>(K_csr.rowsize()),
			vcp::tsparse_scalar::decimal_power_negative<real_type>(14));
		valid_ = data_.factorized && !data_.singular_or_unstable;
		zero_pivots_ = data_.zero_pivots;
		diagnostics_ = data_.diagnostics;
	}

	bool        valid()       const { return valid_; }
	std::size_t zero_pivots() const { return zero_pivots_; }

	const std::string& diagnostics() const { return diagnostics_; }

	void apply(const std::vector<value_type>& r, std::vector<value_type>& z) const {
		if (!valid_) {
			vcp::throw_error<vcp::numerical_error>(
				"ilu0_preconditioner::apply: invalid factorization; ",
				diagnostics_);
		}
		if (r.size() != data_.n) {
			vcp::throw_error<vcp::dimension_error>(
				"ilu0_preconditioner::apply: dimension mismatch");
		}
		z = vcp::tsparse_factorization::ilu0_solve(data_, r);
	}

private:
	vcp::tsparse_factorization::ilu0_data<value_type, index_type> data_;
	bool        valid_;
	std::size_t zero_pivots_;
	std::string diagnostics_;
};

// -----------------------------------------------------------------------
// Trait helpers (SFINAE) for accessing optional preconditioner members.
// These are used by spmatrix to generically read diagnostics/validity.
// -----------------------------------------------------------------------
namespace tsparse_prec_traits {

template <class P>
std::string get_diagnostics_impl(const P& p, int,
    typename std::enable_if<!std::is_void<
        decltype(p.diagnostics())>::value>::type* = 0)
{ return p.diagnostics(); }

template <class P>
std::string get_diagnostics_impl(const P&, ...) { return "user_preconditioner"; }

template <class P>
std::string get_diagnostics(const P& p) {
	return get_diagnostics_impl(p, 0);
}

template <class P>
bool get_valid_impl(const P& p, int,
    typename std::enable_if<!std::is_void<
        decltype(p.valid())>::value>::type* = 0)
{ return p.valid(); }

template <class P>
bool get_valid_impl(const P&, ...) { return true; }

template <class P>
bool get_valid(const P& p) {
	return get_valid_impl(p, 0);
}

template <class P>
std::size_t get_zero_pivots_impl(const P& p, int,
    typename std::enable_if<!std::is_void<
        decltype(p.zero_pivots())>::value>::type* = 0)
{ return p.zero_pivots(); }

template <class P>
std::size_t get_zero_pivots_impl(const P&, ...) { return 0; }

template <class P>
std::size_t get_zero_pivots(const P& p) {
	return get_zero_pivots_impl(p, 0);
}

} // namespace tsparse_prec_traits

} // namespace vcp

#endif
