// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATRIX_HPP
#define VCP_SPMATRIX_HPP

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <complex>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <vcp/spmats.hpp>
#include <vcp/tsparse/tsparse_lanczos.hpp>
#include <vcp/tsparse/tsparse_arnoldi.hpp>
#include <vcp/tsparse/tsparse_factorization.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>

namespace vcp {

	template <typename T> struct spmatrix_is_complex : vcp::tsparse_scalar::is_complex<T> {};

	template <typename T> struct spmatrix_real_type {
		typedef typename vcp::tsparse_scalar::real_type<T>::type type;
	};
	template <typename T> struct spmatrix_real_type<std::complex<T> > {
		typedef typename vcp::tsparse_scalar::real_type<std::complex<T> >::type type;
	};

	template <typename T> struct eig_value_traits {
		typedef T real_type;
		typedef std::complex<T> complex_type;
	};

	template <typename T> struct eig_value_traits<std::complex<T> > {
		typedef T real_type;
		typedef std::complex<T> complex_type;
	};

	enum class linear_solver_method {
		jacobi,
		gauss_seidel,
		conjugate_gradient,
		bicgstab,
		gmres
	};

	enum class preconditioner_type {
		none,
		jacobi
	};

	// eig_method, eigs_target, generalized_eig_method, matrix_structure_hint,
	// orthogonalization_method, eig_solver_method, eig_target
	// are all defined in <vcp/tsparse/tsparse_eigs.hpp> (included via tsparse.hpp).

	template <typename T>
	struct linear_solve_options {
		typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
		linear_solver_method method;
		std::size_t max_iter;
		real_type tol;
		bool check_symmetric;
		bool use_relative_residual;
		std::size_t restart;
		preconditioner_type preconditioner;

		linear_solve_options()
			: method(linear_solver_method::conjugate_gradient), max_iter(1000),
			  tol(vcp::tsparse_scalar::decimal_power_negative<real_type>(12)),
			  check_symmetric(true), use_relative_residual(true), restart(30),
			  preconditioner(preconditioner_type::none) {}
	};

	template <typename T>
	struct eig_options {
		typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
		eig_solver_method method;
		matrix_structure_hint structure;
		eig_target target;
		real_type tol;
		real_type shift;
		bool use_shift;
		std::size_t max_iter;
		std::size_t subspace_dim;
		bool allow_dense_conversion;
		std::size_t max_dense_size;
		orthogonalization_method orthogonalization;
		unsigned int random_seed;
		bool random_start;
		bool compute_residual_history;

		eig_options()
			: method(eig_solver_method::lanczos),
			  structure(matrix_structure_hint::auto_detect),
			  target(eig_target::smallest_algebraic),
			  tol(vcp::tsparse_scalar::decimal_power_negative<real_type>(12)),
			  shift(real_type(0)), use_shift(false),
			  max_iter(1000), subspace_dim(0),
			  allow_dense_conversion(true), max_dense_size(1000000),
			  orthogonalization(orthogonalization_method::modified_gram_schmidt),
			  random_seed(0), random_start(false),
			  compute_residual_history(false) {}
	};

	template <class T> struct linear_solve_result {
		typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
		std::vector<T> x;
		std::vector<T> solution;
		bool converged;
		std::size_t iterations;
		real_type residual_norm;
		real_type absolute_residual_norm;
		real_type relative_residual_norm;
		real_type initial_residual_norm;
		linear_solver_method method;

		linear_solve_result()
			: converged(false), iterations(0), residual_norm((std::numeric_limits<real_type>::infinity)()),
			  absolute_residual_norm((std::numeric_limits<real_type>::infinity)()),
			  relative_residual_norm((std::numeric_limits<real_type>::infinity)()),
			  initial_residual_norm(real_type(0)), method(linear_solver_method::conjugate_gradient) {}
	};

	template <class T> struct eig_result {
		typedef typename vcp::tsparse_scalar::real_type<T>::type real_type;
		typedef typename eig_value_traits<T>::complex_type eigenvalue_type;
		std::vector<T> eigenvalues;
		std::vector<eigenvalue_type> complex_eigenvalues;
		std::vector<std::vector<T> > eigenvectors;
		bool converged;
		std::size_t requested_count;
		std::size_t returned_count;
		std::size_t returned_real_count;
		std::size_t returned_complex_count;
		std::size_t converged_count;
		std::size_t iterations;
		std::size_t matrix_vector_products;
		std::size_t linear_solves;
		std::vector<real_type> residuals_absolute;
		std::vector<real_type> residuals_relative;
		real_type residual_norm_absolute;
		real_type residual_norm_relative;
		std::vector<real_type> residual_history_absolute;
		std::vector<real_type> residual_history_relative;
		eig_solver_method method;             // eig_method (backward compat) field
		std::string status;
		std::string message;
		std::string failure_reason;
		std::string breakdown_reason;
		std::string used_method;              // string version (changed from eig_method)
		std::string used_orthogonalization;
		bool used_dense_fallback;
		bool used_shift_invert;
		bool used_generalized_operator;
		std::size_t used_subspace_dim;
		std::size_t inner_iterations;
		std::size_t inner_failure_count;
		real_type inner_residual_norm;
		std::string inner_failure_reason;
		std::string factorization_diagnostics;
		std::size_t factorization_zero_pivots;

		eig_result()
			: eigenvalues(), complex_eigenvalues(), eigenvectors(),
			  converged(false), requested_count(0), returned_count(0),
			  returned_real_count(0), returned_complex_count(0), converged_count(0),
			  iterations(0), matrix_vector_products(0), linear_solves(0),
			  residuals_absolute(), residuals_relative(),
			  residual_norm_absolute((std::numeric_limits<real_type>::infinity)()),
			  residual_norm_relative((std::numeric_limits<real_type>::infinity)()),
			  residual_history_absolute(), residual_history_relative(),
			  method(eig_solver_method::lanczos), status(), message(), failure_reason(),
			  breakdown_reason(), used_method(), used_orthogonalization(),
			  used_dense_fallback(false), used_shift_invert(false), used_generalized_operator(false),
			  used_subspace_dim(0), inner_iterations(0), inner_failure_count(0),
			  inner_residual_norm(real_type(0)), inner_failure_reason(),
			  factorization_diagnostics(), factorization_zero_pivots(0) {}
	};

	template <typename _T, class _P = spmats<_T> > class spmatrix : protected _P {
	public:
		typedef _T value_type;
		typedef _P policy_type;
		typedef typename _P::index_type index_type;
		typedef typename _P::format_type format_type;
		typedef typename spmatrix_real_type<_T>::type scalar_real_type;
		typedef vcp::linear_solve_options<_T> linear_solve_options_type;
		typedef vcp::eig_options<_T> eig_options_type;

		spmatrix() : _P() {}
		spmatrix(const index_type rows, const index_type cols) : _P() { this->resize(rows, cols); }
		~spmatrix() = default;
		spmatrix(const spmatrix&) = default;
		spmatrix(spmatrix&&) = default;
		spmatrix& operator=(const spmatrix&) = default;
		spmatrix& operator=(spmatrix&&) = default;

		index_type rowsize() const { return _P::rowsize(); }
		index_type columnsize() const { return _P::columnsize(); }
		index_type nnz() const { return _P::nnz(); }
		index_type stored_nnz() const { return _P::stored_nnz(); }
		bool is_finalized() const { return _P::is_finalized(); }
		bool is_sorted() const { return _P::is_sorted(); }
		bool is_unique() const { return _P::is_unique(); }
		format_type format() const { return _P::format(); }

		void resize(const index_type rows, const index_type cols) { _P::resize(rows, cols); }
		void clear() { _P::clear(); }
		void reserve(const index_type n) { _P::reserve(n); }

		void add(const index_type i, const index_type j, const _T& value) { _P::add(i, j, value); }
		void set(const index_type i, const index_type j, const _T& value) { _P::set(i, j, value); }
		_T get(const index_type i, const index_type j) const { return _P::get(i, j); }

		void sort_coo() { _P::sort_coo(); }
		void normalize_coo() { _P::normalize_coo(); }
		void finalize() { _P::finalize(); }
		void to_csr() { _P::to_csr(); }
		void to_csc() { _P::to_csc(); }

		spmatrix as_csr() const {
			spmatrix B;
			static_cast<_P&>(B) = _P::as_csr();
			return B;
		}

		spmatrix as_csc() const {
			spmatrix B;
			static_cast<_P&>(B) = _P::as_csc();
			return B;
		}

		std::vector<_T> mul_vec(const std::vector<_T>& x) const { return _P::mul_vec(x); }
		void mul_vec(const _T* x, _T* y) const { _P::mul_vec(x, y); }
		std::vector<_T> trans_mul_vec(const std::vector<_T>& x) const { return _P::trans_mul_vec(x); }
		void trans_mul_vec(const _T* x, _T* y) const { _P::trans_mul_vec(x, y); }

		std::vector<std::vector<_T> > to_dense() const {
			std::vector<std::vector<_T> > dense(static_cast<std::size_t>(rowsize()), std::vector<_T>(static_cast<std::size_t>(columnsize()), _T(0)));
			spmatrix A = this->as_csr();
			const std::vector<index_type>& outer = A.outer_index();
			const std::vector<index_type>& inner = A.inner_index();
			const std::vector<_T>& val = A.values();
			for (index_type i = 0; i < A.rowsize(); i++) {
				for (index_type p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
					dense[static_cast<std::size_t>(i)][static_cast<std::size_t>(inner[static_cast<std::size_t>(p)])] = val[static_cast<std::size_t>(p)];
				}
			}
			return dense;
		}

		bool is_symmetric() const {
			return is_symmetric(vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(12));
		}

		bool is_symmetric(const scalar_real_type& tol) const {
			if (tol <= scalar_real_type(0)) vcp::throw_error<vcp::invalid_argument>("spmatrix::is_symmetric: tol must be positive");
			if (rowsize() != columnsize()) return false;
			spmatrix A = this->as_csr();
			const std::vector<index_type>& outer = A.outer_index();
			const std::vector<index_type>& inner = A.inner_index();
			const std::vector<_T>& val = A.values();
			for (index_type i = 0; i < A.rowsize(); i++) {
				for (index_type p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
					const index_type j = inner[static_cast<std::size_t>(p)];
					if (i == j) continue;
					const index_type first = outer[static_cast<std::size_t>(j)];
					const index_type last = outer[static_cast<std::size_t>(j + 1)];
					const typename std::vector<index_type>::const_iterator begin = inner.begin() + first;
					const typename std::vector<index_type>::const_iterator end = inner.begin() + last;
					typename std::vector<index_type>::const_iterator it = std::lower_bound(begin, end, i);
					_T mirrored = _T(0);
					if (it != end && *it == i) {
						mirrored = val[static_cast<std::size_t>(it - inner.begin())];
					}
					if (abs_value(val[static_cast<std::size_t>(p)] - mirrored) > tol) return false;
				}
			}
			return true;
		}

		spmatrix operator+(const spmatrix& rhs) const { return add(rhs); }
		spmatrix operator-(const spmatrix& rhs) const { return sub(rhs); }
		spmatrix operator*(const spmatrix& rhs) const { return matmul(rhs); }

		spmatrix add(const spmatrix& rhs) const {
			if (rowsize() != rhs.rowsize() || columnsize() != rhs.columnsize()) {
				vcp::throw_error<vcp::dimension_error>("spmatrix::add: dimension mismatch");
			}
			spmatrix A = this->as_csr();
			spmatrix B = rhs.as_csr();
			spmatrix C(rowsize(), columnsize());
			C.reserve(A.nnz() + B.nnz());
			vcp::tsparse_spgemm::csr_csr_linear_combination(A.rowsize(), A.columnsize(),
				A.outer_index(), A.inner_index(), A.values(),
				B.outer_index(), B.inner_index(), B.values(),
				_T(1), _T(1),
				[&](const index_type i, const index_type j, const _T& value) {
					C.add(i, j, value);
				});
			C.finalize();
			return C;
		}

		spmatrix sub(const spmatrix& rhs) const {
			if (rowsize() != rhs.rowsize() || columnsize() != rhs.columnsize()) {
				vcp::throw_error<vcp::dimension_error>("spmatrix::sub: dimension mismatch");
			}
			spmatrix A = this->as_csr();
			spmatrix B = rhs.as_csr();
			spmatrix C(rowsize(), columnsize());
			C.reserve(A.nnz() + B.nnz());
			vcp::tsparse_spgemm::csr_csr_linear_combination(A.rowsize(), A.columnsize(),
				A.outer_index(), A.inner_index(), A.values(),
				B.outer_index(), B.inner_index(), B.values(),
				_T(1), _T(-1),
				[&](const index_type i, const index_type j, const _T& value) {
					C.add(i, j, value);
				});
			C.finalize();
			return C;
		}

		spmatrix matmul(const spmatrix& rhs) const {
			if (columnsize() != rhs.rowsize()) {
				vcp::throw_error<vcp::dimension_error>("spmatrix::matmul: dimension mismatch");
			}
			spmatrix A = this->as_csr();
			spmatrix B = rhs.as_csr();
			spmatrix C(rowsize(), rhs.columnsize());
			const std::vector<index_type>& ao = A.outer_index();
			const std::vector<index_type>& ai = A.inner_index();
			const std::vector<_T>& av = A.values();
			const std::vector<index_type>& bo = B.outer_index();
			const std::vector<index_type>& bi = B.inner_index();
			const std::vector<_T>& bv = B.values();
			vcp::tsparse_spgemm::csr_csr_multiply(A.rowsize(), A.columnsize(), B.columnsize(), ao, ai, av, bo, bi, bv,
				[&](const index_type i, const index_type j, const _T& value) {
					C.add(i, j, value);
				});
			C.finalize();
			return C;
		}

		std::vector<_T> solve(const std::vector<_T>& b, const linear_solve_options_type& options = linear_solve_options_type()) const {
			require_real_scalar("spmatrix::solve");
			linear_solve_result<_T> result = solve_with_info(b, options);
			if (!result.converged) {
				vcp::throw_error<vcp::state_error>("spmatrix::solve: iterative solver did not converge");
			}
			return result.x;
		}

		linear_solve_result<_T> solve_with_info(const std::vector<_T>& b, const linear_solve_options_type& options = linear_solve_options_type()) const {
			require_real_scalar("spmatrix::solve_with_info");
			switch (options.method) {
			case linear_solver_method::jacobi:
				return solve_jacobi_with_info(b, options.max_iter, options.tol, options.use_relative_residual);
			case linear_solver_method::gauss_seidel:
				return solve_gauss_seidel_with_info(b, options.max_iter, options.tol, options.use_relative_residual);
			case linear_solver_method::conjugate_gradient:
				return solve_cg_with_info(b, options.max_iter, options.tol, options.check_symmetric, options.preconditioner, options.use_relative_residual);
			case linear_solver_method::bicgstab:
				return solve_bicgstab_with_info(b, options.max_iter, options.tol, options.use_relative_residual, options.preconditioner);
			case linear_solver_method::gmres:
				return solve_gmres_with_info(b, options.max_iter, options.tol, options.restart, options.use_relative_residual, options.preconditioner);
			}
			vcp::throw_error<vcp::state_error>("spmatrix::solve_with_info: unknown method");
			return linear_solve_result<_T>();
		}

		std::vector<_T> solve_jacobi(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			require_real_scalar("spmatrix::solve_jacobi");
			return checked_solve_result(solve_jacobi_with_info(b, max_iter, tol, true), "spmatrix::solve_jacobi");
		}

		std::vector<_T> solve_gauss_seidel(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			require_real_scalar("spmatrix::solve_gauss_seidel");
			return checked_solve_result(solve_gauss_seidel_with_info(b, max_iter, tol, true), "spmatrix::solve_gauss_seidel");
		}

		std::vector<_T> solve_cg(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			require_real_scalar("spmatrix::solve_cg");
			return checked_solve_result(solve_cg_with_info(b, max_iter, tol, true, preconditioner_type::none, true), "spmatrix::solve_cg");
		}

		std::vector<_T> solve_bicgstab(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			require_real_scalar("spmatrix::solve_bicgstab");
			return checked_solve_result(solve_bicgstab_with_info(b, max_iter, tol, true), "spmatrix::solve_bicgstab");
		}

		std::vector<_T> solve_gmres(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol) const {
			require_real_scalar("spmatrix::solve_gmres");
			return checked_solve_result(solve_gmres_with_info(b, max_iter, tol, 30, true), "spmatrix::solve_gmres");
		}

		linear_solve_result<_T> solve_jacobi_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                               const bool use_relative_residual = true) const {
			require_real_scalar("spmatrix::solve_jacobi");
			validate_solve_input(b, max_iter, tol, "spmatrix::solve_jacobi");
			const index_type n = rowsize();
			std::vector<_T> diag(static_cast<std::size_t>(n), _T(0));
			for (index_type i = 0; i < n; i++) {
				diag[static_cast<std::size_t>(i)] = get(i, i);
				if (diag[static_cast<std::size_t>(i)] == _T(0)) {
					vcp::throw_error<vcp::numerical_error>("spmatrix::solve_jacobi: zero diagonal");
				}
			}
			spmatrix A = this->as_csr();
			std::vector<_T> x(static_cast<std::size_t>(n), _T(0));
			std::vector<_T> next(static_cast<std::size_t>(n), _T(0));
			linear_solve_result<_T> result;
			result.converged = false;
			result.iterations = 0;
			result.method = linear_solver_method::jacobi;
			const vcp::tsparse_solvers::residual_control<_T> control =
				vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
			scalar_real_type current_residual = residual_norm_value(A, x, b);
			result.residual_norm = current_residual;
			result.converged = current_residual <= control.threshold;
			for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
				for (index_type i = 0; i < n; i++) {
					_T sigma = _T(0);
					const std::vector<index_type>& outer = A.outer_index();
					const std::vector<index_type>& inner = A.inner_index();
					const std::vector<_T>& val = A.values();
					for (index_type p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
						const index_type j = inner[static_cast<std::size_t>(p)];
						if (j != i) sigma += val[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(j)];
					}
					next[static_cast<std::size_t>(i)] = (b[static_cast<std::size_t>(i)] - sigma) / diag[static_cast<std::size_t>(i)];
				}
				x.swap(next);
				current_residual = residual_norm_value(A, x, b);
				result.residual_norm = current_residual;
				result.iterations = iter;
				if (current_residual <= control.threshold) {
					result.converged = true;
					break;
				}
			}
			result.x = x;
			set_linear_residual_fields(result, A, b);
			return result;
		}

		linear_solve_result<_T> solve_gauss_seidel_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                                     const bool use_relative_residual = true) const {
			require_real_scalar("spmatrix::solve_gauss_seidel");
			validate_solve_input(b, max_iter, tol, "spmatrix::solve_gauss_seidel");
			const index_type n = rowsize();
			for (index_type i = 0; i < n; i++) {
				if (get(i, i) == _T(0)) {
					vcp::throw_error<vcp::numerical_error>("spmatrix::solve_gauss_seidel: zero diagonal");
				}
			}
			spmatrix A = this->as_csr();
			const std::vector<index_type>& outer = A.outer_index();
			const std::vector<index_type>& inner = A.inner_index();
			const std::vector<_T>& val = A.values();
			std::vector<_T> x(static_cast<std::size_t>(n), _T(0));
			linear_solve_result<_T> result;
			result.converged = false;
			result.iterations = 0;
			result.method = linear_solver_method::gauss_seidel;
			const vcp::tsparse_solvers::residual_control<_T> control =
				vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
			scalar_real_type current_residual = residual_norm_value(A, x, b);
			result.residual_norm = current_residual;
			result.converged = current_residual <= control.threshold;
			for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
				for (index_type i = 0; i < n; i++) {
					_T sigma = _T(0);
					_T diag = _T(0);
					for (index_type p = outer[static_cast<std::size_t>(i)]; p < outer[static_cast<std::size_t>(i + 1)]; p++) {
						const index_type j = inner[static_cast<std::size_t>(p)];
						if (j == i) diag = val[static_cast<std::size_t>(p)];
						else sigma += val[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(j)];
					}
					x[static_cast<std::size_t>(i)] = (b[static_cast<std::size_t>(i)] - sigma) / diag;
				}
				current_residual = residual_norm_value(A, x, b);
				result.residual_norm = current_residual;
				result.iterations = iter;
				if (current_residual <= control.threshold) {
					result.converged = true;
					break;
				}
			}
			result.x = x;
			set_linear_residual_fields(result, A, b);
			return result;
		}

		linear_solve_result<_T> solve_cg_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                           const bool check_symmetric = true,
		                                           const preconditioner_type preconditioner = preconditioner_type::none,
		                                           const bool use_relative_residual = true) const {
			require_real_scalar("spmatrix::solve_cg");
			validate_solve_input(b, max_iter, tol, "spmatrix::solve_cg");
			if (check_symmetric && !is_symmetric(vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(10))) {
				vcp::throw_error<vcp::domain_error>("spmatrix::solve_cg: matrix must be symmetric");
			}
			spmatrix A = this->as_csr();
			const std::size_t n = b.size();
			std::vector<_T> x(n, _T(0));
			std::vector<_T> r = b;
			const std::vector<_T> inv_diag = make_jacobi_inverse_diagonal(A, preconditioner, "spmatrix::solve_cg");
			std::vector<_T> z = apply_left_preconditioner(r, inv_diag, preconditioner);
			std::vector<_T> p = z;
			scalar_real_type rzold = dot_value(r, z);
			const vcp::tsparse_solvers::residual_control<_T> control =
				vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
			linear_solve_result<_T> result;
			scalar_real_type current_residual = norm_value(r);
			result.converged = current_residual <= control.threshold;
			result.iterations = 0;
			result.residual_norm = current_residual;
			result.method = linear_solver_method::conjugate_gradient;
			for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
				const std::vector<_T> Ap = A.mul_vec(p);
				const scalar_real_type denom = dot_value(p, Ap);
				if (!(denom > scalar_real_type(0)) || !vcp::tsparse_scalar::is_finite(denom)) {
					result.converged = false;
					break;
				}
				const _T alpha = _T(rzold / denom);
				for (std::size_t i = 0; i < n; i++) {
					x[i] += alpha * p[i];
					r[i] -= alpha * Ap[i];
				}
				current_residual = norm_value(r);
				result.residual_norm = current_residual;
				result.iterations = iter;
				if (current_residual <= control.threshold) {
					result.converged = true;
					break;
				}
				z = apply_left_preconditioner(r, inv_diag, preconditioner);
				const scalar_real_type rznew = dot_value(r, z);
				if (!(rznew > scalar_real_type(0)) || !vcp::tsparse_scalar::is_finite(rznew)) {
					result.converged = false;
					break;
				}
				const _T beta = _T(rznew / rzold);
				for (std::size_t i = 0; i < n; i++) p[i] = z[i] + beta * p[i];
				rzold = rznew;
			}
			result.x = x;
			set_linear_residual_fields(result, A, b);
			return result;
		}

		linear_solve_result<_T> solve_bicgstab_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                                 const bool use_relative_residual = true,
		                                                 const preconditioner_type preconditioner = preconditioner_type::none) const {
			require_real_scalar("spmatrix::solve_bicgstab");
			validate_solve_input(b, max_iter, tol, "spmatrix::solve_bicgstab");
			spmatrix A = this->as_csr();
			const std::size_t n = b.size();
			const std::vector<_T> inv_diag = make_jacobi_inverse_diagonal(A, preconditioner, "spmatrix::solve_bicgstab");
			const std::vector<_T> pb = apply_left_preconditioner(b, inv_diag, preconditioner);
			const vcp::tsparse_solvers::residual_control<_T> control =
				vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
			std::vector<_T> x(n, _T(0));
			std::vector<_T> r = pb;
			std::vector<_T> r_hat = r;
			std::vector<_T> p(n, _T(0));
			std::vector<_T> v(n, _T(0));
			scalar_real_type rho_old = scalar_real_type(1);
			scalar_real_type alpha = scalar_real_type(1);
			scalar_real_type omega = scalar_real_type(1);

			linear_solve_result<_T> result;
			result.x = x;
			scalar_real_type current_residual = residual_norm_value(A, x, b);
			result.residual_norm = current_residual;
			result.converged = current_residual <= control.threshold;
			result.iterations = 0;
			result.method = linear_solver_method::bicgstab;
			const scalar_real_type eps = std::numeric_limits<scalar_real_type>::epsilon();
			for (std::size_t iter = 1; iter <= max_iter && !result.converged; iter++) {
				const scalar_real_type rho = dot_value(r_hat, r);
				if (vcp::tsparse_scalar::abs_value(rho) <= eps || !vcp::tsparse_scalar::is_finite(rho)) break;
				const scalar_real_type beta = (rho / rho_old) * (alpha / omega);
				for (std::size_t i = 0; i < n; i++) p[i] = r[i] + _T(beta) * (p[i] - _T(omega) * v[i]);
				v = apply_preconditioned_operator(A, p, inv_diag, preconditioner);
				const scalar_real_type denom = dot_value(r_hat, v);
				if (vcp::tsparse_scalar::abs_value(denom) <= eps || !vcp::tsparse_scalar::is_finite(denom)) break;
				alpha = rho / denom;
				std::vector<_T> s(n);
				for (std::size_t i = 0; i < n; i++) s[i] = r[i] - _T(alpha) * v[i];
				const scalar_real_type snorm = norm_value(s);
				if (snorm <= control.threshold) {
					for (std::size_t i = 0; i < n; i++) x[i] += _T(alpha) * p[i];
					current_residual = residual_norm_value(A, x, b);
					result.residual_norm = current_residual;
					result.iterations = iter;
					result.converged = current_residual <= control.threshold;
					break;
				}
				const std::vector<_T> t = apply_preconditioned_operator(A, s, inv_diag, preconditioner);
				const scalar_real_type tt = dot_value(t, t);
				if (tt <= eps || !vcp::tsparse_scalar::is_finite(tt)) break;
				omega = dot_value(t, s) / tt;
				if (vcp::tsparse_scalar::abs_value(omega) <= eps || !vcp::tsparse_scalar::is_finite(omega)) break;
				for (std::size_t i = 0; i < n; i++) {
					x[i] += _T(alpha) * p[i] + _T(omega) * s[i];
					r[i] = s[i] - _T(omega) * t[i];
				}
				current_residual = residual_norm_value(A, x, b);
				result.residual_norm = current_residual;
				result.iterations = iter;
				result.converged = current_residual <= control.threshold;
				rho_old = rho;
			}
			result.x = x;
			set_linear_residual_fields(result, A, b);
			return result;
		}

		linear_solve_result<_T> solve_gmres_with_info(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol,
		                                             const std::size_t restart, const bool use_relative_residual = true,
		                                             const preconditioner_type preconditioner = preconditioner_type::none) const {
			require_real_scalar("spmatrix::solve_gmres");
			validate_solve_input(b, max_iter, tol, "spmatrix::solve_gmres");
			if (restart == 0) vcp::throw_error<vcp::invalid_argument>("spmatrix::solve_gmres: restart must be positive");
			spmatrix A = this->as_csr();
			const std::size_t n = b.size();
			const std::vector<_T> inv_diag = make_jacobi_inverse_diagonal(A, preconditioner, "spmatrix::solve_gmres");
			const std::vector<_T> pb = apply_left_preconditioner(b, inv_diag, preconditioner);
			const vcp::tsparse_solvers::residual_control<_T> control =
				vcp::tsparse_solvers::make_residual_control(b, tol, use_relative_residual);
			std::vector<_T> x(n, _T(0));
			linear_solve_result<_T> result;
			result.x = x;
			result.converged = false;
			result.iterations = 0;
			scalar_real_type current_residual = residual_norm_value(A, x, b);
			result.residual_norm = current_residual;
			result.method = linear_solver_method::gmres;
			if (current_residual <= control.threshold) {
				result.converged = true;
				set_linear_residual_fields(result, A, b);
				return result;
			}
			const std::size_t mmax = std::min(restart, n);
			while (result.iterations < max_iter && !result.converged) {
				std::vector<_T> Ax = A.mul_vec(x);
				const std::vector<_T> pAx = apply_left_preconditioner(Ax, inv_diag, preconditioner);
				std::vector<_T> r(n);
				for (std::size_t i = 0; i < n; i++) r[i] = pb[i] - pAx[i];
				const scalar_real_type beta = norm_value(r);
				current_residual = residual_norm_value(A, x, b);
				result.residual_norm = current_residual;
				if (current_residual <= control.threshold) {
					result.converged = true;
					break;
				}
				std::vector<std::vector<_T> > V(mmax + 1, std::vector<_T>(n, _T(0)));
				std::vector<std::vector<_T> > H(mmax + 1, std::vector<_T>(mmax, _T(0)));
				std::vector<scalar_real_type> cs(mmax, scalar_real_type(0));
				std::vector<scalar_real_type> sn(mmax, scalar_real_type(0));
				std::vector<_T> g(mmax + 1, _T(0));
				g[0] = _T(beta);
				for (std::size_t i = 0; i < n; i++) V[0][i] = r[i] / _T(beta);
				std::size_t m = 0;
				for (; m < mmax && result.iterations < max_iter; m++) {
					std::vector<_T> w = apply_preconditioned_operator(A, V[m], inv_diag, preconditioner);
					for (std::size_t j = 0; j <= m; j++) {
						H[j][m] = _T(dot_value(w, V[j]));
						for (std::size_t i = 0; i < n; i++) w[i] -= H[j][m] * V[j][i];
					}
					const scalar_real_type hnext = norm_value(w);
					H[m + 1][m] = _T(hnext);
					if (hnext > std::numeric_limits<scalar_real_type>::epsilon()) {
						for (std::size_t i = 0; i < n; i++) V[m + 1][i] = w[i] / _T(hnext);
					}
					for (std::size_t j = 0; j < m; j++) {
						const _T hij = H[j][m];
						const _T hip1j = H[j + 1][m];
						H[j][m] = _T(cs[j]) * hij + _T(sn[j]) * hip1j;
						H[j + 1][m] = _T(-sn[j]) * hij + _T(cs[j]) * hip1j;
					}
					const scalar_real_type h0 = abs_value(H[m][m]);
					const scalar_real_type h1 = abs_value(H[m + 1][m]);
					const scalar_real_type rho = vcp::tsparse_scalar::hypot_value(h0, h1);
					if (rho <= std::numeric_limits<scalar_real_type>::epsilon()) {
						cs[m] = scalar_real_type(1);
						sn[m] = scalar_real_type(0);
					}
					else {
						cs[m] = vcp::tsparse_scalar::real_part(H[m][m]) / rho;
						sn[m] = vcp::tsparse_scalar::real_part(H[m + 1][m]) / rho;
					}
					const _T hmm = H[m][m];
					const _T hp1m = H[m + 1][m];
					H[m][m] = _T(cs[m]) * hmm + _T(sn[m]) * hp1m;
					H[m + 1][m] = _T(0);
					const _T gm = g[m];
					g[m] = _T(cs[m]) * gm;
					g[m + 1] = _T(-sn[m]) * gm;
					result.iterations++;
					const scalar_real_type projected_residual = abs_value(g[m + 1]);
					result.residual_norm = projected_residual;
					if (projected_residual <= control.threshold) {
						m++;
						std::vector<_T> y = solve_upper_triangular(H, g, m);
						for (std::size_t j = 0; j < m; j++) {
							for (std::size_t i = 0; i < n; i++) x[i] += V[j][i] * y[j];
						}
						current_residual = residual_norm_value(A, x, b);
						result.residual_norm = current_residual;
						result.converged = current_residual <= control.threshold;
						break;
					}
					if (hnext <= std::numeric_limits<scalar_real_type>::epsilon()) {
						m++;
						break;
					}
				}
				if (m == 0) break;
				if (!result.converged) {
					std::vector<_T> y = solve_upper_triangular(H, g, m);
					for (std::size_t j = 0; j < m; j++) {
						for (std::size_t i = 0; i < n; i++) x[i] += V[j][i] * y[j];
					}
					current_residual = residual_norm_value(A, x, b);
					result.residual_norm = current_residual;
					result.converged = current_residual <= control.threshold;
				}
			}
			result.x = x;
			set_linear_residual_fields(result, A, b);
			return result;
		}

		eig_result<_T> eig(const eig_options_type& options) const {
			require_real_scalar("spmatrix::eig");
			eig_result<_T> result = eig_with_info(options);
			if (!result.converged) vcp::throw_error<vcp::state_error>("spmatrix::eig: eigensolver did not converge");
			return result;
		}

		eig_result<_T> eig_with_info(const eig_options_type& options) const {
			require_real_scalar("spmatrix::eig");
			validate_eig_input("spmatrix::eig");
			if (options.max_iter == 0 || options.tol <= scalar_real_type(0)) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eig: invalid iteration option");
			}
			if (options.method != eig_solver_method::dense_fallback_explicit) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eig: full dense eig requires dense_fallback_explicit");
			}
			check_dense_allowed(options, "spmatrix::eig");
			std::vector<std::vector<_T> > dense = to_dense();
			eig_result<_T> result = dense_eig(dense, options);
			result.used_dense_fallback = true;
			set_result_counts(result, static_cast<std::size_t>(rowsize()));
			return result;
		}

		eig_options_type default_eigs_options() const {
			eig_options_type options;
			options.structure = matrix_structure_hint::auto_detect;
			return options;
		}

		static eig_options_type dense_eigs_options() {
			eig_options_type options;
			options.method = eig_solver_method::dense_fallback_explicit;
			return options;
		}

		std::vector<_T> eigvals(const eig_options_type& options) const {
			return eig(options).eigenvalues;
		}

		std::vector<_T> eigs(const std::size_t k) const {
			require_real_scalar("spmatrix::eigs");
			return eigs(k, default_eigs_options());
		}

		std::vector<_T> eigs(const std::size_t k, const eig_options_type& options) const {
			require_real_scalar("spmatrix::eigs");
			eig_result<_T> result = eigs_with_info(k, options);
			if (!result.converged) vcp::throw_error<vcp::state_error>("spmatrix::eigs: eigensolver did not converge");
			if (result.returned_real_count < k) vcp::throw_error<vcp::state_error>("spmatrix::eigs: insufficient real eigenvalues returned");
			return result.eigenvalues;
		}

		eig_result<_T> eigs_with_info(const std::size_t k) const {
			require_real_scalar("spmatrix::eigs");
			return eigs_with_info(k, default_eigs_options());
		}

		eig_result<_T> eigs_with_info(const std::size_t k, const eig_options_type& options) const {
			require_real_scalar("spmatrix::eigs");
			eig_options_type active = resolve_eigs_options(options);
			validate_eigs_input(k, "spmatrix::eigs", active.method);
			if (active.max_iter == 0 || active.tol <= scalar_real_type(0)) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eigs: invalid iteration option");
			}
			if (active.method == eig_solver_method::lanczos
			 || active.method == eig_solver_method::shift_invert_lanczos)
				return lanczos_eigs_new(k, active);
			if (active.method == eig_solver_method::arnoldi
			 || active.method == eig_solver_method::shift_invert_arnoldi)
				return arnoldi_eigs_new(k, active);
			if (active.method != eig_solver_method::dense_fallback_explicit) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eigs: unknown eigensolver method");
			}
			eig_result<_T> result = eig_with_info(active);
			select_eigenpairs(result, k, active.target, active.shift);
			set_result_counts(result, k);
			return result;
		}

		std::vector<_T> eig(const spmatrix& B) const {
			require_real_scalar("spmatrix::eig(A,B)");
			return eig(B, eig_options_type()).eigenvalues;
		}

		eig_result<_T> eig(const spmatrix& B, const eig_options_type& options) const {
			require_real_scalar("spmatrix::eig(A,B)");
			eig_result<_T> result = eigs_with_info(B, static_cast<std::size_t>(rowsize()), options);
			if (!result.converged) vcp::throw_error<vcp::state_error>("spmatrix::eig(A,B): eigensolver did not converge");
			return result;
		}

		std::vector<_T> eigs(const spmatrix& B, const std::size_t k, const eig_options_type& options = eig_options_type()) const {
			require_real_scalar("spmatrix::eigs(A,B)");
			eig_result<_T> result = eigs_with_info(B, k, options);
			if (!result.converged) vcp::throw_error<vcp::state_error>("spmatrix::eigs(A,B): eigensolver did not converge");
			if (result.returned_real_count < k) vcp::throw_error<vcp::state_error>("spmatrix::eigs(A,B): insufficient real eigenvalues returned");
			return result.eigenvalues;
		}

		eig_result<_T> eigs_with_info(const spmatrix& B, const std::size_t k, const eig_options_type& options = eig_options_type()) const {
			require_real_scalar("spmatrix::eigs(A,B)");
			validate_generalized_eig_input(B, "spmatrix::eigs(A,B)");
			if (k == 0 || k > static_cast<std::size_t>(rowsize())) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eigs(A,B): invalid k");
			}
			if (options.max_iter == 0 || options.tol <= scalar_real_type(0)) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix::eigs(A,B): invalid iteration option");
			}
			if (is_diagonal_matrix(*this) && is_diagonal_matrix(B)) {
				return generalized_diagonal_eigs(B, k, options);
			}
			eig_result<_T> result;
			result.requested_count = k;
						result.converged = false;
			result.status = "unsupported_generalized_non_diagonal";
			result.used_generalized_operator = true;
			result.used_dense_fallback = false;
			result.method = options.method;
			result.used_method = eig_method_to_string(options.method);
			result.failure_reason = "non-diagonal generalized sparse eigs is outside the merge API; use diagonal generalized path or future generalized shift-invert";
			result.message = result.failure_reason;
			set_result_counts(result, k);
			return result;
		}

		spmatrix transpose() const {
			spmatrix B;
			static_cast<_P&>(B) = _P::transpose();
			return B;
		}

		void assign_csr(const index_type rows, const index_type cols,
		                const std::vector<index_type>& row_ptr,
		                const std::vector<index_type>& col_ind,
		                const std::vector<_T>& value) {
			_P::assign_csr(rows, cols, row_ptr, col_ind, value);
		}

		void assign_csc(const index_type rows, const index_type cols,
		                const std::vector<index_type>& col_ptr,
		                const std::vector<index_type>& row_ind,
		                const std::vector<_T>& value) {
			_P::assign_csc(rows, cols, col_ptr, row_ind, value);
		}

		const std::vector<index_type>& outer_index() const { return _P::outer_index(); }
		const std::vector<index_type>& inner_index() const { return _P::inner_index(); }
		const std::vector<_T>& values() const { return _P::values(); }
		const std::vector<index_type>& coo_rows() const { return _P::coo_rows(); }
		const std::vector<index_type>& coo_columns() const { return _P::coo_columns(); }
		const std::vector<_T>& coo_values() const { return _P::coo_values(); }

		friend std::vector<_T> operator*(const spmatrix& A, const std::vector<_T>& x) {
			return A.mul_vec(x);
		}

		friend spmatrix transpose(const spmatrix& A) {
			return A.transpose();
		}

	private:
		static std::vector<_T> checked_solve_result(const linear_solve_result<_T>& result, const char* routine) {
			if (!result.converged) vcp::throw_error<vcp::state_error>(routine, ": iterative solver did not converge");
			return result.x;
		}

		static void require_real_scalar(const char* routine) {
			if (spmatrix_is_complex<_T>::value) {
				vcp::throw_error<vcp::domain_error>(routine, ": std::complex scalar is supported only for basic sparse operations");
			}
		}

		static scalar_real_type abs_value(const _T& x) {
			return vcp::tsparse_scalar::abs_value(x);
		}

		void validate_solve_input(const std::vector<_T>& b, const std::size_t max_iter, const scalar_real_type& tol, const char* routine) const {
			if (rowsize() != columnsize()) vcp::throw_error<vcp::dimension_error>(routine, ": matrix must be square");
			if (b.size() != static_cast<std::size_t>(rowsize())) vcp::throw_error<vcp::dimension_error>(routine, ": right hand side dimension mismatch");
			if (max_iter == 0) vcp::throw_error<vcp::invalid_argument>(routine, ": max_iter must be positive");
			if (tol <= scalar_real_type(0)) vcp::throw_error<vcp::invalid_argument>(routine, ": tol must be positive");
		}

		void validate_eig_input(const char* routine) const {
			if (rowsize() != columnsize()) vcp::throw_error<vcp::dimension_error>(routine, ": matrix must be square");
		}

		void validate_eigs_input(const std::size_t k, const char* routine, const eig_method = eig_solver_method::lanczos) const {
			validate_eig_input(routine);
			if (k == 0) vcp::throw_error<vcp::invalid_argument>(routine, ": k must be positive");
			if (k > static_cast<std::size_t>(rowsize())) vcp::throw_error<vcp::invalid_argument>(routine, ": k is larger than matrix size");
		}

		void validate_generalized_eig_input(const spmatrix& B, const char* routine) const {
			if (rowsize() != columnsize() || B.rowsize() != B.columnsize()) {
				vcp::throw_error<vcp::dimension_error>(routine, ": matrices must be square");
			}
			if (rowsize() != B.rowsize()) {
				vcp::throw_error<vcp::dimension_error>(routine, ": matrix sizes must match");
			}
		}

		static bool is_diagonal_matrix(const spmatrix& A) {
			spmatrix C = A.as_csr();
			const std::vector<index_type>& outer = C.outer_index();
			const std::vector<index_type>& inner = C.inner_index();
			const std::vector<_T>& val = C.values();
			for (std::size_t i = 0; i < static_cast<std::size_t>(C.rowsize()); i++) {
				for (index_type p = outer[i]; p < outer[i + 1]; p++) {
					const std::size_t j = static_cast<std::size_t>(inner[static_cast<std::size_t>(p)]);
					if (i != j && tsparse_scalar::abs_value(val[static_cast<std::size_t>(p)]) > scalar_real_type(0)) return false;
				}
			}
			return true;
		}

		eig_result<_T> generalized_diagonal_eigs(const spmatrix& B, const std::size_t k,
		                                        const eig_options_type& options) const {
			const std::size_t n = static_cast<std::size_t>(rowsize());
			eig_result<_T> result;
			result.requested_count = k;
						result.method = options.method;
			result.used_method = eig_method_to_string(options.method);
			result.used_generalized_operator = true;
			result.used_dense_fallback = false;
			std::vector<_T> values;
			std::vector<std::vector<_T> > vectors;
			values.reserve(n);
			vectors.reserve(n);
			const scalar_real_type pivot_tol = std::numeric_limits<scalar_real_type>::epsilon();
			for (std::size_t i = 0; i < n; i++) {
				const _T bdiag = B.get(static_cast<index_type>(i), static_cast<index_type>(i));
				if (tsparse_scalar::abs_value(bdiag) <= pivot_tol) {
					vcp::throw_error<vcp::numerical_error>("spmatrix::eigs(A,B): zero diagonal in B");
				}
				values.push_back(get(static_cast<index_type>(i), static_cast<index_type>(i)) / bdiag);
				std::vector<_T> e(n, _T(0));
				e[i] = _T(1);
				vectors.push_back(e);
			}
			const std::vector<std::size_t> order =
				vcp::tsparse_eigen_selection::select_real_eigenpairs(values, k, options.target, options.shift);
			for (std::size_t i = 0; i < order.size(); i++) {
				result.eigenvalues.push_back(values[order[i]]);
				result.eigenvectors.push_back(vectors[order[i]]);
			}
			result.residuals_absolute = generalized_eigenpair_residuals(*this, B, result.eigenvalues, result.eigenvectors);
			result.residuals_relative = generalized_eigenpair_relative_residuals(*this, B, result.eigenvalues, result.eigenvectors);
			result.converged = true;
			result.status = "converged";
			result.message = "converged";
			populate_real_complex_eigenvalues(result);
			set_result_counts(result, k);
			return result;
		}

		void check_dense_allowed(const eig_options_type& options, const char*) const {
			vcp::tsparse_dense_fallback::check_dense_size(static_cast<std::size_t>(rowsize()),
				static_cast<std::size_t>(columnsize()), options.max_dense_size, options.allow_dense_conversion);
		}

		eig_options_type resolve_eigs_options(const eig_options_type& options) const {
			eig_options_type active = options;
			if (active.structure == matrix_structure_hint::hermitian) {
				vcp::throw_error<vcp::state_error>("spmatrix::eigs: hermitian structure hint is reserved for future complex support");
				return active;
			}
			if (active.structure == matrix_structure_hint::symmetric) {
				if (rowsize() != columnsize()) {
					vcp::throw_error<vcp::dimension_error>("spmatrix::eigs: symmetric structure hint requires a square matrix");
					return active;
				}
				if (active.method == eig_solver_method::arnoldi) active.method = eig_solver_method::lanczos;
			}
			else if (active.structure == matrix_structure_hint::general) {
				if (active.method == eig_solver_method::lanczos) active.method = eig_solver_method::arnoldi;
			}
			if (active.use_shift) {
				if (active.method == eig_solver_method::lanczos) active.method = eig_solver_method::shift_invert_lanczos;
				else if (active.method == eig_solver_method::arnoldi) active.method = eig_solver_method::shift_invert_arnoldi;
				else if (active.method == eig_solver_method::dense_fallback_explicit) {
					vcp::throw_error<vcp::invalid_argument>("spmatrix::eigs: shift is not supported by dense_fallback_explicit");
				}
			}
			if (active.structure == matrix_structure_hint::auto_detect
			 && active.method == eig_solver_method::lanczos
			 && !is_symmetric(vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(10))) {
				vcp::throw_error<vcp::domain_error>("spmatrix::eigs: Lanczos requires a symmetric matrix");
			}
			return active;
		}

		static scalar_real_type dot_value(const std::vector<_T>& a, const std::vector<_T>& b) {
			return vcp::tsparse_scalar::real_dot_value(a, b);
		}

		static scalar_real_type norm_value(const std::vector<_T>& a) {
			return vcp::tsparse_scalar::real_norm_value(a);
		}

		static scalar_real_type residual_norm_value(const spmatrix& A, const std::vector<_T>& x, const std::vector<_T>& b) {
			return vcp::tsparse_solvers::residual_norm_value(A, x, b);
		}

		static void set_linear_residual_fields(linear_solve_result<_T>& result, const spmatrix& A, const std::vector<_T>& b) {
			vcp::tsparse_solvers::set_linear_residual_fields(result, A, b);
		}

		static std::vector<_T> make_jacobi_inverse_diagonal(const spmatrix& A, const preconditioner_type preconditioner, const char* routine) {
			std::vector<_T> inv_diag;
			if (preconditioner == preconditioner_type::none) return inv_diag;
			if (preconditioner != preconditioner_type::jacobi) {
				vcp::throw_error<vcp::state_error>(routine, ": unknown preconditioner");
				return inv_diag;
			}
			const std::size_t n = static_cast<std::size_t>(A.rowsize());
			inv_diag.assign(n, _T(0));
			for (std::size_t i = 0; i < n; i++) {
				const _T diag = A.get(static_cast<index_type>(i), static_cast<index_type>(i));
				if (diag == _T(0)) {
					vcp::throw_error<vcp::numerical_error>(routine, ": Jacobi preconditioner has zero diagonal");
				}
				inv_diag[i] = _T(1) / diag;
			}
			return inv_diag;
		}

		static std::vector<_T> apply_left_preconditioner(const std::vector<_T>& r, const std::vector<_T>& inv_diag,
		                                                 const preconditioner_type preconditioner) {
			if (preconditioner == preconditioner_type::none) return r;
			std::vector<_T> z(r.size(), _T(0));
			for (std::size_t i = 0; i < r.size(); i++) z[i] = inv_diag[i] * r[i];
			return z;
		}

		static std::vector<_T> apply_preconditioned_operator(const spmatrix& A, const std::vector<_T>& x,
		                                                     const std::vector<_T>& inv_diag,
		                                                     const preconditioner_type preconditioner) {
			const std::vector<_T> ax = A.mul_vec(x);
			return apply_left_preconditioner(ax, inv_diag, preconditioner);
		}

		static std::string eig_method_to_string(const eig_solver_method m) {
			switch (m) {
			case eig_solver_method::lanczos: return "lanczos";
			case eig_solver_method::arnoldi: return "arnoldi";
			case eig_solver_method::shift_invert_lanczos: return "shift_invert_lanczos";
			case eig_solver_method::shift_invert_arnoldi: return "shift_invert_arnoldi";
			case eig_solver_method::dense_fallback_explicit: return "dense_fallback_explicit";
			}
			return "unknown";
		}

		static void set_eig_diagnostics(eig_result<_T>& result, const eig_solver_method method,
		                                const std::size_t subspace_dim, const std::size_t matvec_count,
		                                const std::string& breakdown, const std::string& failure) {
			result.method = method;
			result.used_method = eig_method_to_string(method);
			result.used_subspace_dim = subspace_dim;
			result.matrix_vector_products = matvec_count;
			result.breakdown_reason = breakdown;
			if (result.converged) {
				result.status = "converged";
				result.message = breakdown.empty() ? "converged" : breakdown;
				result.failure_reason.clear();
			}
			else {
				result.status = "not_converged";
				result.failure_reason = failure.empty() ? "residual tolerance not reached" : failure;
				result.message = result.failure_reason;
			}
			set_result_counts(result, result.requested_count);
		}

		static eig_result<_T> convert_dense_result(const vcp::tsparse_dense_linalg::dense_eigen_result<_T>& source,
		                                          const eig_solver_method method) {
			eig_result<_T> result;
			result.eigenvalues = source.eigenvalues;
			result.complex_eigenvalues = source.complex_eigenvalues;
			result.eigenvectors = source.eigenvectors;
			result.converged = source.converged;
			result.iterations = source.iterations;
			result.residuals_absolute = source.residuals;
			result.residual_norm_absolute = source.residual_norm;
			set_eig_diagnostics(result, method, source.eigenvalues.size(), 0, std::string(), std::string());
			return result;
		}

		static void set_result_counts(eig_result<_T>& result, const std::size_t requested) {
			result.requested_count = requested;
			result.returned_real_count = result.eigenvalues.size();
			result.returned_complex_count = 0;
			for (std::size_t i = 0; i < result.complex_eigenvalues.size(); i++) {
				if (tsparse_scalar::abs_value(result.complex_eigenvalues[i].imag())
					> std::numeric_limits<scalar_real_type>::epsilon()) {
					result.returned_complex_count++;
				}
			}
			result.returned_count = result.returned_real_count + result.returned_complex_count;
			if (result.converged_count == 0 && result.converged) result.converged_count = result.returned_count;
			if (!result.residuals_absolute.empty()) {
				result.residual_norm_absolute = *std::max_element(result.residuals_absolute.begin(), result.residuals_absolute.end());
			}
			if (!result.residuals_relative.empty()) {
				result.residual_norm_relative = *std::max_element(result.residuals_relative.begin(), result.residuals_relative.end());
			}
		}

		static std::string orthogonalization_to_string(const orthogonalization_method method) {
			switch (method) {
			case orthogonalization_method::modified_gram_schmidt: return "modified_gram_schmidt";
			case orthogonalization_method::classical_gram_schmidt_twice: return "classical_gram_schmidt_twice";
			}
			return "modified_gram_schmidt";
		}

		static void sort_eigenvalues(std::vector<_T>& values) {
			std::sort(values.begin(), values.end(), [](const _T& a, const _T& b) {
				return vcp::tsparse_scalar::real_part(a) < vcp::tsparse_scalar::real_part(b);
			});
		}

		static void sort_eigenpairs(eig_result<_T>& result) {
			if (result.eigenvalues.empty()) return;  // no real eigenvalues to sort; preserve complex_eigenvalues
			std::vector<std::size_t> order(result.eigenvalues.size());
			for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
			std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
				return vcp::tsparse_scalar::real_part(result.eigenvalues[a]) <
				       vcp::tsparse_scalar::real_part(result.eigenvalues[b]);
			});
			std::vector<_T> values(order.size());
			std::vector<typename eig_result<_T>::eigenvalue_type> cvalues;
			if (result.complex_eigenvalues.size() == result.eigenvalues.size()) cvalues.resize(order.size());
			std::vector<scalar_real_type> residuals_abs;
			if (result.residuals_absolute.size() == result.eigenvalues.size()) residuals_abs.resize(order.size());
			std::vector<scalar_real_type> residuals_rel;
			if (result.residuals_relative.size() == result.eigenvalues.size()) residuals_rel.resize(order.size());
			std::vector<std::vector<_T> > vectors;
			if (result.eigenvectors.size() == result.eigenvalues.size()) vectors.resize(order.size());
			for (std::size_t i = 0; i < order.size(); i++) {
				values[i] = result.eigenvalues[order[i]];
				if (!cvalues.empty()) cvalues[i] = result.complex_eigenvalues[order[i]];
				if (!residuals_abs.empty()) residuals_abs[i] = result.residuals_absolute[order[i]];
				if (!residuals_rel.empty()) residuals_rel[i] = result.residuals_relative[order[i]];
				if (!vectors.empty()) vectors[i] = result.eigenvectors[order[i]];
			}
			result.eigenvalues.swap(values);
			if (!cvalues.empty()) result.complex_eigenvalues.swap(cvalues);
			else if (!result.complex_eigenvalues.empty() && result.complex_eigenvalues.size() == result.eigenvalues.size())
				result.complex_eigenvalues.clear();
			if (!residuals_abs.empty()) result.residuals_absolute.swap(residuals_abs);
			else result.residuals_absolute.clear();
			if (!residuals_rel.empty()) result.residuals_relative.swap(residuals_rel);
			else result.residuals_relative.clear();
			if (!vectors.empty()) result.eigenvectors.swap(vectors);
		}

		static void select_eigenpairs(eig_result<_T>& result, const std::size_t k, const eig_target target,
		                              const scalar_real_type& shift = scalar_real_type(0)) {
			const std::vector<std::size_t> order =
				vcp::tsparse_eigen_selection::select_real_eigenpairs(result.eigenvalues, k, target, shift);
			std::vector<_T> values;
			std::vector<std::vector<_T> > vectors;
			std::vector<scalar_real_type> residuals_abs;
			std::vector<scalar_real_type> residuals_rel;
			values.reserve(order.size());
			if (result.eigenvectors.size() == result.eigenvalues.size()) vectors.reserve(order.size());
			if (result.residuals_absolute.size() == result.eigenvalues.size()) residuals_abs.reserve(order.size());
			if (result.residuals_relative.size() == result.eigenvalues.size()) residuals_rel.reserve(order.size());
			for (std::size_t i = 0; i < order.size(); i++) {
				const std::size_t j = order[i];
				values.push_back(result.eigenvalues[j]);
				if (result.eigenvectors.size() == result.eigenvalues.size()) vectors.push_back(result.eigenvectors[j]);
				if (result.residuals_absolute.size() == result.eigenvalues.size()) residuals_abs.push_back(result.residuals_absolute[j]);
				if (result.residuals_relative.size() == result.eigenvalues.size()) residuals_rel.push_back(result.residuals_relative[j]);
			}
			result.eigenvalues.swap(values);
			if (result.eigenvectors.size() == vectors.size()) result.eigenvectors.swap(vectors);
			if (result.residuals_absolute.size() == residuals_abs.size()) result.residuals_absolute.swap(residuals_abs);
			if (result.residuals_relative.size() == residuals_rel.size()) result.residuals_relative.swap(residuals_rel);
			populate_real_complex_eigenvalues(result);
			set_result_counts(result, k);
		}

		static void populate_real_complex_eigenvalues(eig_result<_T>& result) {
			if (!result.complex_eigenvalues.empty() && result.eigenvalues.empty()) return;
			result.complex_eigenvalues.clear();
			result.complex_eigenvalues.reserve(result.eigenvalues.size());
			for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
				result.complex_eigenvalues.push_back(typename eig_result<_T>::eigenvalue_type(result.eigenvalues[i]));
			}
		}

		static scalar_real_type residual_tolerance(const scalar_real_type& tol, const std::size_t scale) {
			const scalar_real_type base(tol);
			const scalar_real_type scaled = base * scalar_real_type(scale) * scalar_real_type(10);
			return scaled > base ? scaled : base;
		}

		static eig_result<_T> dense_eig(std::vector<std::vector<_T> > dense, const eig_options_type& options) {
			eig_result<_T> result = is_dense_symmetric(dense, options.tol * scalar_real_type(10))
				? jacobi_eig_dense(dense, options.max_iter, options.tol)
				: qr_eig_dense(dense, options.max_iter, options.tol);
			result.method = eig_solver_method::dense_fallback_explicit;
			result.used_method = eig_method_to_string(eig_solver_method::dense_fallback_explicit);
			result.used_dense_fallback = true;
			return result;
		}

		static eig_result<_T> jacobi_eig_dense(std::vector<std::vector<_T> > a, const std::size_t max_iter, const scalar_real_type& tol) {
			return convert_dense_result(vcp::tsparse_dense_linalg::jacobi_eig_dense(a, max_iter, tol), eig_solver_method::dense_fallback_explicit);
		}

		static eig_result<_T> qr_eig_dense(std::vector<std::vector<_T> > a, const std::size_t max_iter, const scalar_real_type& tol) {
			return convert_dense_result(vcp::tsparse_dense_linalg::qr_eig_dense(a, max_iter, tol), eig_solver_method::dense_fallback_explicit);
		}

		eig_result<_T> power_eig(const std::size_t max_iter, const scalar_real_type& tol) const {
			validate_eig_input("spmatrix::power_eig");
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());
			std::vector<_T> x(n, _T(1));
			eig_result<_T> result;
			result.converged = false;
			result.iterations = 0;
			result.residual_norm_absolute = std::numeric_limits<scalar_real_type>::infinity();
			result.method = eig_solver_method::arnoldi;
			const scalar_real_type tolerance(tol);
			for (std::size_t iter = 1; iter <= max_iter; iter++) {
				std::vector<_T> y = A.mul_vec(x);
				const scalar_real_type ny = norm_value(y);
				if (ny <= std::numeric_limits<scalar_real_type>::epsilon()) {
					vcp::throw_error<vcp::numerical_error>("spmatrix::power_eig: zero vector breakdown");
				}
				for (std::size_t i = 0; i < n; i++) x[i] = y[i] / _T(ny);
				y = A.mul_vec(x);
				const scalar_real_type lambda = dot_value(x, y);
				const scalar_real_type residual_value = eigenpair_residual_norm_value(A, _T(lambda), x);
				result.residual_norm_absolute = residual_value;
				result.iterations = iter;
				if (residual_value <= tolerance) {
					result.converged = true;
					result.eigenvalues.assign(1, _T(lambda));
					result.eigenvectors.assign(1, x);
					break;
				}
				result.eigenvalues.assign(1, _T(lambda));
				result.eigenvectors.assign(1, x);
			}
			populate_real_complex_eigenvalues(result);
			set_eig_diagnostics(result, eig_solver_method::arnoldi, 1, result.iterations * 2,
				std::string(), std::string());
			return result;
		}

		eig_result<_T> lanczos_eigs(const std::size_t k, const eig_options_type& options) const {
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());
			const std::size_t m_limit = krylov_subspace_dim(n, k, options.subspace_dim);
			struct apply_type {
				const spmatrix* A;
				void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = A->mul_vec(x); }
			} apply = { &A };
			const vcp::tsparse_eigensolvers::lanczos_tridiagonalization<_T> decomp =
				vcp::tsparse_eigensolvers::build_lanczos_tridiagonalization<_T>(n, m_limit, options.max_iter,
					options.tol, true, apply);
			const std::size_t m = decomp.alpha.size();
			if (m == 0) vcp::throw_error<vcp::numerical_error>("spmatrix::lanczos_eigs: no Krylov vector generated");
			const std::size_t projected_iter = vcp::tsparse_eigensolvers::projected_max_iter(options.max_iter, m);
			vcp::tsparse_dense_linalg::dense_eigen_result<_T> small =
				vcp::tsparse_eigensolvers::extract_lanczos_ritz(decomp, projected_iter, options.tol / scalar_real_type(10));
			vcp::tsparse_dense_linalg::select_eigenpairs(small, k, options.target == eig_target::largest_algebraic || options.target == eig_target::largest_magnitude);
			eig_result<_T> result = convert_dense_result(small, eig_solver_method::lanczos);
			result.eigenvectors = lift_ritz_vectors(decomp.basis, result.eigenvectors, n);
			result.iterations = m;
			const scalar_real_type residual_value = max_eigenpair_residual_value(A, result.eigenvalues, result.eigenvectors);
			result.residual_norm_absolute = residual_value;
			result.converged = result.converged && residual_value <= scalar_real_type(options.tol);
			populate_real_complex_eigenvalues(result);
			set_eig_diagnostics(result, eig_solver_method::lanczos, m, decomp.matrix_vector_products,
				decomp.breakdown_reason, std::string());
			return result;
		}

		eig_result<_T> arnoldi_eigs(const std::size_t k, const eig_options_type& options) const {
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());
			const std::size_t m_limit = krylov_subspace_dim(n, k, options.subspace_dim);
			struct apply_type {
				const spmatrix* A;
				void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = A->mul_vec(x); }
			} apply = { &A };
			const vcp::tsparse_eigensolvers::arnoldi_factorization<_T> decomp =
				vcp::tsparse_eigensolvers::build_arnoldi_factorization<_T>(n, m_limit, options.max_iter,
					options.tol, true, options.orthogonalization, apply);
			const std::size_t m = decomp.basis_size;
			if (m == 0) vcp::throw_error<vcp::numerical_error>("spmatrix::arnoldi_eigs: no Krylov vector generated");
			const std::vector<std::vector<_T> > Hm = vcp::tsparse_eigensolvers::square_hessenberg(decomp.hessenberg, m);
			const std::size_t projected_iter = vcp::tsparse_eigensolvers::projected_max_iter(options.max_iter, m);
			vcp::tsparse_dense_linalg::dense_eigen_result<_T> small =
				vcp::tsparse_eigensolvers::extract_arnoldi_ritz(decomp, projected_iter, options.tol / scalar_real_type(10));
			vcp::tsparse_dense_linalg::select_eigenpairs(small, k, options.target == eig_target::largest_algebraic);
			eig_result<_T> result = convert_dense_result(small, eig_solver_method::arnoldi);
			if (result.eigenvectors.size() != result.eigenvalues.size()) {
				result.eigenvectors.clear();
				for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
					result.eigenvectors.push_back(dense_eigenvector_inverse_iteration(Hm, result.eigenvalues[i]));
				}
			}
			result.eigenvectors = lift_ritz_vectors(decomp.basis, result.eigenvectors, n);
			result.iterations = m;
			const scalar_real_type residual_value = max_eigenpair_residual_value(A, result.eigenvalues, result.eigenvectors);
			result.residual_norm_absolute = residual_value;
			result.converged = result.converged && residual_value <= scalar_real_type(options.tol);
			populate_real_complex_eigenvalues(result);
			set_eig_diagnostics(result, eig_solver_method::arnoldi, m, decomp.matrix_vector_products,
				decomp.breakdown_reason, std::string());
			return result;
		}

		// ------------------------------------------------------------------
		// Standard Lanczos implementation using tsparse_lanczos
		// ------------------------------------------------------------------
		eig_result<_T> lanczos_eigs_new(const std::size_t k, const eig_options_type& options) const {
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());
			const std::size_t sdim = (options.subspace_dim == 0) ? std::max(k + 5, std::min(n, std::size_t(30))) : options.subspace_dim;
			const scalar_real_type shift_val = options.shift;

			// Build apply functor for shift-invert if needed
			if (options.method == eig_solver_method::shift_invert_lanczos) {
				return shift_invert_lanczos_eigs(k, options, shift_val);
			}

			struct apply_fn {
				const spmatrix* mat;
				void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
			} apply = { &A };

			auto pkg = vcp::tsparse_lanczos::lanczos_eigs_standard<_T, apply_fn>(
				n, k, sdim, options.max_iter, options.tol,
				options.random_seed, options.random_start,
				options.target, shift_val, options.compute_residual_history, apply);

			return lanczos_package_to_result(pkg, k, eig_solver_method::lanczos);
		}

		template <class ApplyA>
		eig_result<_T> lanczos_package_to_result(
			const vcp::tsparse_lanczos::lanczos_result_package<_T, ApplyA>& pkg,
			const std::size_t k,
			const eig_solver_method method) const
		{
			eig_result<_T> result;
			result.eigenvalues = pkg.eigenvalues;
			result.eigenvectors = pkg.eigenvectors;
			result.converged = pkg.converged;
			result.iterations = pkg.iterations;
			result.converged_count = pkg.converged_count;
			result.returned_count = pkg.returned_count;
			result.matrix_vector_products = pkg.mv_count;
			result.residuals_absolute = pkg.residuals_abs;
			result.residuals_relative = pkg.residuals_rel;
			result.residual_history_absolute = pkg.history_abs;
			result.residual_history_relative = pkg.history_rel;
				result.breakdown_reason = pkg.breakdown_reason;
				result.failure_reason = pkg.failure_reason;
			if (!pkg.residuals_abs.empty()) {
				result.residual_norm_absolute = *std::max_element(pkg.residuals_abs.begin(), pkg.residuals_abs.end());
			}
			populate_real_complex_eigenvalues(result);
			if (!result.eigenvectors.empty()) {
				spmatrix A = this->as_csr();
				result.residuals_absolute = eigenpair_residuals(A, result.eigenvalues, result.eigenvectors);
				result.residuals_relative = eigenpair_relative_residuals(A, result.eigenvalues, result.eigenvectors);
			}
			set_eig_diagnostics(result, method, result.returned_count, result.matrix_vector_products,
				result.breakdown_reason, result.failure_reason);
			result.used_shift_invert = (method == eig_solver_method::shift_invert_lanczos
			                         || method == eig_solver_method::shift_invert_arnoldi);
			set_result_counts(result, k);
			return result;
		}

		eig_result<_T> shift_invert_lanczos_eigs(const std::size_t k, const eig_options_type& options,
		                                          const scalar_real_type& sigma) const {
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());

			// Build (A - sigma I) as sparse matrix
			spmatrix shifted = A;
			for (std::size_t i = 0; i < n; i++) {
				const _T cur = shifted.get(static_cast<index_type>(i), static_cast<index_type>(i));
				shifted.set(static_cast<index_type>(i), static_cast<index_type>(i), cur - _T(sigma));
			}
			shifted.finalize();
			spmatrix shiftedCSR = shifted.as_csr();

			// ILU(0) factorization of (A - sigma I)
			typedef vcp::tsparse_factorization::ilu0_data<_T, index_type> ILU;
			ILU ilu = vcp::tsparse_factorization::ilu0_factorize<_T, index_type>(
				shiftedCSR.outer_index(), shiftedCSR.inner_index(), shiftedCSR.values(), n,
				vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(14));
			if (ilu.singular_or_unstable) {
				eig_result<_T> result;
				result.requested_count = k;
				result.method = eig_solver_method::shift_invert_lanczos;
				result.used_method = eig_method_to_string(eig_solver_method::shift_invert_lanczos);
				result.used_shift_invert = true;
				result.used_dense_fallback = false;
				result.status = "factorization_failed";
				result.failure_reason = "ILU zero or near-zero pivot detected";
				result.message = result.failure_reason;
				result.factorization_diagnostics = ilu.diagnostics;
				result.factorization_zero_pivots = ilu.zero_pivots;
				set_result_counts(result, k);
				return result;
			}

			// Shift-invert apply: y = (A - sigma I)^{-1} x (via GMRES with ILU prec)
			const std::size_t inner_max = std::min(n, std::size_t(100));
			const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
			std::size_t linear_solve_count = 0;
			std::size_t inner_failure_count = 0;
			std::size_t inner_iteration_count = 0;
			scalar_real_type inner_residual_norm = scalar_real_type(0);

			struct apply_fn {
				const spmatrix* mat;
				const ILU* ilu_ptr;
				std::size_t inner_max;
				scalar_real_type inner_tol;
				scalar_real_type sigma;
				std::size_t* solve_count;
				std::size_t* failure_count;
				std::size_t* iteration_count;
				scalar_real_type* max_inner_residual;
				void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
					const std::size_t nn = x.size();
					struct av { const spmatrix* m; scalar_real_type s; void operator()(const std::vector<_T>& u, std::vector<_T>& v) const {
						v = m->mul_vec(u);
						for(std::size_t i=0;i<v.size();i++) v[i]-=_T(s)*u[i]; } } av_fn = {mat, sigma};
					struct pv { const ILU* p; void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
						z = vcp::tsparse_factorization::ilu0_solve(*p, r); } } pv_fn = {ilu_ptr};
					auto gr = vcp::tsparse_factorization::gmres_solve<_T, av, pv>(
						av_fn, pv_fn, x, inner_max / 10 + 1, inner_tol, std::min(nn, inner_max));
					(*solve_count)++;
					(*iteration_count) += gr.iterations;
					if (gr.residual_norm > *max_inner_residual) *max_inner_residual = gr.residual_norm;
					if (!gr.converged) (*failure_count)++;
					y = gr.x;
				}
			} apply_si = { &A, &ilu, inner_max, inner_tol, sigma, &linear_solve_count,
				&inner_failure_count, &inner_iteration_count, &inner_residual_norm };

			const std::size_t sdim = (options.subspace_dim == 0) ? std::max(k + 5, std::min(n, std::size_t(30))) : options.subspace_dim;
			auto pkg = vcp::tsparse_lanczos::lanczos_eigs_standard<_T, apply_fn>(
				n, k, sdim, options.max_iter, options.tol,
				options.random_seed, options.random_start,
				eig_target::largest_magnitude, scalar_real_type(0), options.compute_residual_history, apply_si);

			// For shift-invert, eigenvalues of B^{-1}A are 1/(lambda - sigma);
			// actual eigenvalues = 1/mu + sigma
			for (std::size_t i = 0; i < pkg.eigenvalues.size(); i++) {
				const scalar_real_type mu = tsparse_scalar::real_part(pkg.eigenvalues[i]);
				if (tsparse_scalar::abs_value(mu) > scalar_real_type(0)) {
					pkg.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
				}
			}
			eig_result<_T> result = lanczos_package_to_result(pkg, k, eig_solver_method::shift_invert_lanczos);
			result.linear_solves = linear_solve_count;
			result.inner_iterations = inner_iteration_count;
			result.inner_failure_count = inner_failure_count;
			result.inner_residual_norm = inner_residual_norm;
			result.factorization_diagnostics = ilu.diagnostics;
			result.factorization_zero_pivots = ilu.zero_pivots;
			result.used_shift_invert = true;
			if (inner_failure_count != 0) {
				result.converged = false;
				result.status = "inner_solve_failed";
				result.failure_reason = "inner GMRES solve failed";
				result.inner_failure_reason = result.failure_reason;
				result.message = result.failure_reason;
			}
			set_result_counts(result, k);
			return result;
		}

		// ------------------------------------------------------------------
		// Standard Arnoldi eigensolver.
		// ------------------------------------------------------------------
		eig_result<_T> arnoldi_eigs_new(const std::size_t k, const eig_options_type& options) const {
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());
			if (options.method == eig_solver_method::shift_invert_arnoldi) {
				return shift_invert_arnoldi_eigs(k, options, options.shift);
			}
			const std::size_t sdim = (options.subspace_dim == 0)
				? std::max(k + 5, std::min(n, std::size_t(30)))
				: options.subspace_dim;
			const std::size_t max_restarts = options.max_iter + static_cast<std::size_t>(options.max_iter * k);

			struct apply_fn {
				const spmatrix* mat;
				void operator()(const std::vector<_T>& x, std::vector<_T>& y) const { y = mat->mul_vec(x); }
			} apply_op = { &A };

			const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
				vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
					n, k, sdim, max_restarts,
					options.tol,
					options.orthogonalization,
					true,
					true,
					options.random_seed,
					options.random_start,
					options.target,
					tsparse_scalar::real_part(options.shift),
					options.compute_residual_history,
					apply_op);

			eig_result<_T> result;
			result.method = eig_solver_method::arnoldi;
			result.used_method = "arnoldi";
			result.iterations = pkg.iterations;
			result.matrix_vector_products = pkg.mv_count;
			result.used_subspace_dim = sdim;
			result.breakdown_reason = pkg.breakdown_reason;
			result.failure_reason = pkg.failure_reason;
			result.converged_count = pkg.converged_count;
			result.returned_count = pkg.returned_count;
			result.used_orthogonalization = orthogonalization_to_string(options.orthogonalization);
				result.residual_history_absolute = pkg.history_abs;
				result.residual_history_relative = pkg.history_rel;

			// Real eigenvalues and eigenvectors
			result.eigenvalues = pkg.eigenvalues;
			result.eigenvectors = pkg.eigenvectors;

			// Complex eigenvalues
			for (std::size_t i = 0; i < pkg.complex_eigenvalues.size(); i++) {
				result.complex_eigenvalues.push_back(
					typename eig_result<_T>::eigenvalue_type(
						pkg.complex_eigenvalues[i].first,
						pkg.complex_eigenvalues[i].second));
			}

			// Residuals
			if (!result.eigenvectors.empty()) {
				result.residuals_absolute = eigenpair_residuals(A, result.eigenvalues, result.eigenvectors);
				result.residuals_relative = eigenpair_relative_residuals(A, result.eigenvalues, result.eigenvectors);
				const scalar_real_type res_val = max_eigenpair_residual_value(A, result.eigenvalues, result.eigenvectors);
				result.residual_norm_absolute = res_val;
			} else if (!pkg.residuals_abs.empty()) {
				result.residual_norm_absolute = pkg.residuals_abs[0];
			}

			const bool has_complex = pkg.has_complex;
			result.converged = pkg.converged && !has_complex && (result.eigenvalues.size() >= k);
			if (has_complex) {
				result.status = "complex_ritz_values";
				result.message = "complex Ritz values detected";
				if (result.failure_reason.empty())
					result.failure_reason = "complex Ritz values in requested subset";
			}
			else {
				set_eig_diagnostics(result, eig_solver_method::arnoldi, sdim, result.matrix_vector_products,
					result.breakdown_reason, result.failure_reason);
			}
			set_result_counts(result, k);

			return result;
		}

		eig_result<_T> shift_invert_arnoldi_eigs(const std::size_t k, const eig_options_type& options,
		                                         const scalar_real_type& sigma) const {
			spmatrix A = this->as_csr();
			const std::size_t n = static_cast<std::size_t>(rowsize());
			spmatrix shifted = A;
			for (std::size_t i = 0; i < n; i++) {
				const _T cur = shifted.get(static_cast<index_type>(i), static_cast<index_type>(i));
				shifted.set(static_cast<index_type>(i), static_cast<index_type>(i), cur - _T(sigma));
			}
			shifted.finalize();
			spmatrix shiftedCSR = shifted.as_csr();
			typedef vcp::tsparse_factorization::ilu0_data<_T, index_type> ILU;
			ILU ilu = vcp::tsparse_factorization::ilu0_factorize<_T, index_type>(
				shiftedCSR.outer_index(), shiftedCSR.inner_index(), shiftedCSR.values(), n,
				vcp::tsparse_scalar::decimal_power_negative<scalar_real_type>(14));
			if (ilu.singular_or_unstable) {
				eig_result<_T> result;
				result.requested_count = k;
				result.method = eig_solver_method::shift_invert_arnoldi;
				result.used_method = eig_method_to_string(eig_solver_method::shift_invert_arnoldi);
				result.used_orthogonalization = orthogonalization_to_string(options.orthogonalization);
				result.used_shift_invert = true;
				result.used_dense_fallback = false;
				result.status = "factorization_failed";
				result.failure_reason = "ILU zero or near-zero pivot detected";
				result.message = result.failure_reason;
				result.factorization_diagnostics = ilu.diagnostics;
				result.factorization_zero_pivots = ilu.zero_pivots;
				set_result_counts(result, k);
				return result;
			}
			const std::size_t inner_max = std::min(n, std::size_t(100));
			const scalar_real_type inner_tol = options.tol / scalar_real_type(1000);
			std::size_t linear_solve_count = 0;
			std::size_t inner_failure_count = 0;
			std::size_t inner_iteration_count = 0;
			scalar_real_type inner_residual_norm = scalar_real_type(0);
			struct apply_fn {
				const spmatrix* mat;
				const ILU* ilu_ptr;
				std::size_t inner_max;
				scalar_real_type inner_tol;
				scalar_real_type sigma;
				std::size_t* solve_count;
				std::size_t* failure_count;
				std::size_t* iteration_count;
				scalar_real_type* max_inner_residual;
				void operator()(const std::vector<_T>& x, std::vector<_T>& y) const {
					const std::size_t nn = x.size();
					struct av { const spmatrix* m; scalar_real_type s; void operator()(const std::vector<_T>& u, std::vector<_T>& v) const {
						v = m->mul_vec(u);
						for (std::size_t i = 0; i < v.size(); i++) v[i] -= _T(s) * u[i]; } } av_fn = {mat, sigma};
					struct pv { const ILU* p; void operator()(const std::vector<_T>& r, std::vector<_T>& z) const {
						z = vcp::tsparse_factorization::ilu0_solve(*p, r); } } pv_fn = {ilu_ptr};
					auto gr = vcp::tsparse_factorization::gmres_solve<_T, av, pv>(
						av_fn, pv_fn, x, inner_max / 10 + 1, inner_tol, std::min(nn, inner_max));
					(*solve_count)++;
					(*iteration_count) += gr.iterations;
					if (gr.residual_norm > *max_inner_residual) *max_inner_residual = gr.residual_norm;
					if (!gr.converged) (*failure_count)++;
					y = gr.x;
				}
			} apply_si = { &A, &ilu, inner_max, inner_tol, sigma, &linear_solve_count,
				&inner_failure_count, &inner_iteration_count, &inner_residual_norm };
			const std::size_t sdim = (options.subspace_dim == 0)
				? std::max(k + 5, std::min(n, std::size_t(30)))
				: options.subspace_dim;
			const vcp::tsparse_arnoldi::arnoldi_result_package<_T> pkg =
				vcp::tsparse_arnoldi::arnoldi_eigs_standard<_T>(
					n, k, sdim, options.max_iter + options.max_iter * k,
					options.tol, options.orthogonalization, true, true,
					options.random_seed, options.random_start, eig_target::largest_magnitude, scalar_real_type(0),
						options.compute_residual_history, apply_si);
			eig_result<_T> result;
			result.requested_count = k;
			result.method = eig_solver_method::shift_invert_arnoldi;
			result.used_method = eig_method_to_string(eig_solver_method::shift_invert_arnoldi);
			result.used_orthogonalization = orthogonalization_to_string(options.orthogonalization);
			result.iterations = pkg.iterations;
			result.matrix_vector_products = pkg.mv_count;
			result.linear_solves = linear_solve_count;
			result.inner_iterations = inner_iteration_count;
			result.inner_failure_count = inner_failure_count;
			result.inner_residual_norm = inner_residual_norm;
			result.factorization_diagnostics = ilu.diagnostics;
			result.factorization_zero_pivots = ilu.zero_pivots;
			result.used_subspace_dim = sdim;
			result.breakdown_reason = pkg.breakdown_reason;
			result.failure_reason = pkg.failure_reason;
			result.converged_count = pkg.converged_count;
				result.residual_history_absolute = pkg.history_abs;
				result.residual_history_relative = pkg.history_rel;
			result.eigenvalues = pkg.eigenvalues;
			for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
				const scalar_real_type mu = tsparse_scalar::real_part(result.eigenvalues[i]);
				if (tsparse_scalar::abs_value(mu) > scalar_real_type(0)) result.eigenvalues[i] = _T(scalar_real_type(1) / mu + sigma);
			}
			result.eigenvectors = pkg.eigenvectors;
			result.used_shift_invert = true;
			result.used_dense_fallback = false;
			if (!result.eigenvectors.empty()) {
				result.residuals_absolute = eigenpair_residuals(A, result.eigenvalues, result.eigenvectors);
				result.residuals_relative = eigenpair_relative_residuals(A, result.eigenvalues, result.eigenvectors);
			}
			result.converged = pkg.converged && result.eigenvalues.size() >= k && inner_failure_count == 0;
			if (inner_failure_count != 0) {
				result.status = "inner_solve_failed";
				result.failure_reason = "inner GMRES solve failed";
				result.inner_failure_reason = result.failure_reason;
				result.message = result.failure_reason;
			}
			else {
				set_eig_diagnostics(result, eig_solver_method::shift_invert_arnoldi, sdim,
					result.matrix_vector_products, result.breakdown_reason, result.failure_reason);
			}
			set_result_counts(result, k);
			return result;
		}

		static std::vector<std::vector<_T> > lift_ritz_vectors(const std::vector<std::vector<_T> >& V,
		                                                       const std::vector<std::vector<_T> >& small_vectors,
		                                                       const std::size_t n) {
			return vcp::tsparse_dense_linalg::lift_ritz_vectors(V, small_vectors, n);
		}

		static scalar_real_type eigenpair_residual_norm_value(const spmatrix& A, const _T& lambda, const std::vector<_T>& v) {
			std::vector<_T> r = A.mul_vec(v);
			for (std::size_t i = 0; i < r.size(); i++) r[i] -= lambda * v[i];
			return norm_value(r);
		}

		static scalar_real_type frobenius_norm_value(const spmatrix& A) {
			spmatrix C = A.as_csr();
			scalar_real_type s(0);
			const std::vector<_T>& val = C.values();
			for (std::size_t i = 0; i < val.size(); i++) {
				const scalar_real_type a = tsparse_scalar::abs_value(val[i]);
				s += a * a;
			}
			return tsparse_scalar::sqrt_value(s);
		}

		static scalar_real_type eigenpair_relative_residual_norm_value(const spmatrix& A, const _T& lambda,
		                                                              const std::vector<_T>& v) {
			const scalar_real_type abs_res = eigenpair_residual_norm_value(A, lambda, v);
			const scalar_real_type vn = norm_value(v);
			const scalar_real_type denom = frobenius_norm_value(A) * vn
				+ tsparse_scalar::abs_value(lambda) * vn
				+ std::numeric_limits<scalar_real_type>::epsilon();
			return abs_res / denom;
		}

		static scalar_real_type eigenpair_residual_norm(const spmatrix& A, const _T& lambda, const std::vector<_T>& v) {
			return eigenpair_residual_norm_value(A, lambda, v);
		}

		static scalar_real_type max_eigenpair_residual(const spmatrix& A, const std::vector<_T>& eigenvalues,
		                                              const std::vector<std::vector<_T> >& eigenvectors) {
			return max_eigenpair_residual_value(A, eigenvalues, eigenvectors);
		}

		static scalar_real_type max_eigenpair_residual_value(const spmatrix& A, const std::vector<_T>& eigenvalues,
		                                                    const std::vector<std::vector<_T> >& eigenvectors) {
			scalar_real_type maximum(0);
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
				return std::numeric_limits<scalar_real_type>::infinity();
			}
			for (std::size_t i = 0; i < eigenvalues.size(); i++) {
				const scalar_real_type residual = eigenpair_residual_norm_value(A, eigenvalues[i], eigenvectors[i]);
				if (residual > maximum) maximum = residual;
			}
			return maximum;
		}

		static std::vector<scalar_real_type> eigenpair_residuals(const spmatrix& A, const std::vector<_T>& eigenvalues,
		                                                        const std::vector<std::vector<_T> >& eigenvectors) {
			std::vector<scalar_real_type> residuals;
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
			residuals.reserve(eigenvalues.size());
			for (std::size_t i = 0; i < eigenvalues.size(); i++) {
				residuals.push_back(eigenpair_residual_norm(A, eigenvalues[i], eigenvectors[i]));
			}
			return residuals;
		}

		static std::vector<scalar_real_type> eigenpair_relative_residuals(const spmatrix& A,
		                                                                 const std::vector<_T>& eigenvalues,
		                                                                 const std::vector<std::vector<_T> >& eigenvectors) {
			std::vector<scalar_real_type> residuals;
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
			residuals.reserve(eigenvalues.size());
			for (std::size_t i = 0; i < eigenvalues.size(); i++) {
				residuals.push_back(eigenpair_relative_residual_norm_value(A, eigenvalues[i], eigenvectors[i]));
			}
			return residuals;
		}

		static scalar_real_type max_generalized_eigenpair_residual(const spmatrix& A, const spmatrix& B,
		                                                          const std::vector<_T>& eigenvalues,
		                                                          const std::vector<std::vector<_T> >& eigenvectors) {
			return max_generalized_eigenpair_residual_value(A, B, eigenvalues, eigenvectors);
		}

		static scalar_real_type generalized_eigenpair_residual_norm_value(const spmatrix& A, const spmatrix& B,
		                                                                 const _T& eigenvalue,
		                                                                 const std::vector<_T>& eigenvector) {
			std::vector<_T> r = A.mul_vec(eigenvector);
			std::vector<_T> bv = B.mul_vec(eigenvector);
			for (std::size_t i = 0; i < r.size(); i++) r[i] -= eigenvalue * bv[i];
			return norm_value(r);
		}

		static scalar_real_type generalized_eigenpair_relative_residual_norm_value(const spmatrix& A, const spmatrix& B,
		                                                                          const _T& eigenvalue,
		                                                                          const std::vector<_T>& eigenvector) {
			const scalar_real_type abs_res = generalized_eigenpair_residual_norm_value(A, B, eigenvalue, eigenvector);
			const scalar_real_type vn = norm_value(eigenvector);
			const scalar_real_type denom = frobenius_norm_value(A) * vn
				+ tsparse_scalar::abs_value(eigenvalue) * frobenius_norm_value(B) * vn
				+ std::numeric_limits<scalar_real_type>::epsilon();
			return abs_res / denom;
		}

		static scalar_real_type max_generalized_eigenpair_residual_value(const spmatrix& A, const spmatrix& B,
		                                                                const std::vector<_T>& eigenvalues,
		                                                                const std::vector<std::vector<_T> >& eigenvectors) {
			scalar_real_type maximum(0);
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
				return std::numeric_limits<scalar_real_type>::infinity();
			}
			for (std::size_t i = 0; i < eigenvalues.size(); i++) {
				const scalar_real_type residual = generalized_eigenpair_residual_norm_value(A, B, eigenvalues[i], eigenvectors[i]);
				if (residual > maximum) maximum = residual;
			}
			return maximum;
		}

		static std::vector<scalar_real_type> generalized_eigenpair_residuals(const spmatrix& A, const spmatrix& B,
		                                                                    const std::vector<_T>& eigenvalues,
		                                                                    const std::vector<std::vector<_T> >& eigenvectors) {
			std::vector<scalar_real_type> residuals;
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
			residuals.reserve(eigenvalues.size());
			for (std::size_t p = 0; p < eigenvalues.size(); p++) {
				residuals.push_back(generalized_eigenpair_residual_norm_value(A, B, eigenvalues[p], eigenvectors[p]));
			}
			return residuals;
		}

		static std::vector<scalar_real_type> generalized_eigenpair_relative_residuals(const spmatrix& A, const spmatrix& B,
		                                                                             const std::vector<_T>& eigenvalues,
		                                                                             const std::vector<std::vector<_T> >& eigenvectors) {
			std::vector<scalar_real_type> residuals;
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
			residuals.reserve(eigenvalues.size());
			for (std::size_t p = 0; p < eigenvalues.size(); p++) {
				residuals.push_back(generalized_eigenpair_relative_residual_norm_value(A, B, eigenvalues[p], eigenvectors[p]));
			}
			return residuals;
		}

		static scalar_real_type b_inner_product_value(const std::vector<_T>& x, const spmatrix& B, const std::vector<_T>& y) {
			const std::vector<_T> by = B.mul_vec(y);
			return dot_value(x, by);
		}

		static scalar_real_type b_inner_product(const std::vector<_T>& x, const spmatrix& B, const std::vector<_T>& y) {
			return b_inner_product_value(x, B, y);
		}

		static scalar_real_type b_norm_value(const std::vector<_T>& x, const spmatrix& B) {
			const scalar_real_type v = b_inner_product_value(x, B, x);
			if (v < scalar_real_type(0)) {
				vcp::throw_error<vcp::numerical_error>("spmatrix::b_norm: B-inner product is negative");
			}
			return vcp::tsparse_scalar::sqrt_value(v > scalar_real_type(0) ? v : scalar_real_type(0));
		}

		static scalar_real_type b_norm(const std::vector<_T>& x, const spmatrix& B) {
			return b_norm_value(x, B);
		}

		static void validate_b_inner_lanczos_basis(const spmatrix& B, const scalar_real_type& tol) {
			if (!B.is_symmetric(tol)) {
				vcp::throw_error<vcp::domain_error>("spmatrix::B-inner Lanczos: B must be symmetric");
			}
		}

		static scalar_real_type max_dense_eigen_residual(const std::vector<std::vector<_T> >& A,
		                                                const std::vector<_T>& eigenvalues,
		                                                const std::vector<std::vector<_T> >& eigenvectors) {
			return max_dense_eigen_residual_value(A, eigenvalues, eigenvectors);
		}

		static scalar_real_type dense_eigenpair_residual_norm_value(const std::vector<std::vector<_T> >& A,
		                                                           const _T& eigenvalue,
		                                                           const std::vector<_T>& eigenvector) {
			std::vector<_T> r(eigenvector.size(), _T(0));
			for (std::size_t i = 0; i < A.size(); i++) {
				for (std::size_t j = 0; j < A[i].size(); j++) r[i] += A[i][j] * eigenvector[j];
				r[i] -= eigenvalue * eigenvector[i];
			}
			return norm_value(r);
		}

		static scalar_real_type max_dense_eigen_residual_value(const std::vector<std::vector<_T> >& A,
		                                                      const std::vector<_T>& eigenvalues,
		                                                      const std::vector<std::vector<_T> >& eigenvectors) {
			scalar_real_type maximum(0);
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) {
				return std::numeric_limits<scalar_real_type>::infinity();
			}
			for (std::size_t i = 0; i < eigenvalues.size(); i++) {
				const scalar_real_type residual = dense_eigenpair_residual_norm_value(A, eigenvalues[i], eigenvectors[i]);
				if (residual > maximum) maximum = residual;
			}
			return maximum;
		}

		static std::vector<scalar_real_type> dense_eigenpair_residuals(const std::vector<std::vector<_T> >& A,
		                                                              const std::vector<_T>& eigenvalues,
		                                                              const std::vector<std::vector<_T> >& eigenvectors) {
			std::vector<scalar_real_type> residuals;
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
			residuals.reserve(eigenvalues.size());
			for (std::size_t i = 0; i < eigenvalues.size(); i++) {
				residuals.push_back(dense_eigenpair_residual_norm_value(A, eigenvalues[i], eigenvectors[i]));
			}
			return residuals;
		}

		static bool is_dense_symmetric(const std::vector<std::vector<_T> >& A, const scalar_real_type& tol) {
			return vcp::tsparse_dense_linalg::is_dense_symmetric(A, tol);
		}

		static std::vector<_T> dense_eigenvector_inverse_iteration(const std::vector<std::vector<_T> >& A, const _T& lambda) {
			return vcp::tsparse_dense_linalg::dense_eigenvector_inverse_iteration(A, lambda);
		}

		static std::size_t krylov_subspace_dim(const std::size_t n, const std::size_t k, const std::size_t requested) {
			return vcp::tsparse_eigensolvers::krylov_subspace_dim(n, k, requested);
		}

		static scalar_real_type lower_offdiag_norm(const std::vector<std::vector<_T> >& A) {
			return vcp::tsparse_dense_linalg::lower_offdiag_norm(A);
		}

		static std::vector<std::vector<_T> > matmul_dense(const std::vector<std::vector<_T> >& A, const std::vector<std::vector<_T> >& B) {
			return vcp::tsparse_dense_linalg::matmul_dense(A, B);
		}

		static std::vector<_T> solve_upper_triangular(const std::vector<std::vector<_T> >& R,
		                                              const std::vector<_T>& rhs,
		                                              const std::size_t n) {
			return vcp::tsparse_dense_linalg::solve_upper_triangular(R, rhs, n);
		}

		static std::vector<_T> solve_dense_gaussian(std::vector<std::vector<_T> > A, std::vector<_T> b) {
			return vcp::tsparse_dense_linalg::solve_dense_gaussian(A, b);
		}

	};

	template <typename _T, typename _Index = int> class spmatrix_builder {
	public:
		typedef _T value_type;
		typedef _Index index_type;
		typedef vcp::spmatrix<value_type, vcp::spmats<value_type, index_type> > matrix_type;

		spmatrix_builder(const index_type rows, const index_type cols)
			: matrix_(rows, cols) {}

		void reserve(const std::size_t n) {
			if (n > static_cast<std::size_t>(std::numeric_limits<index_type>::max())) {
				vcp::throw_error<vcp::invalid_argument>("spmatrix_builder::reserve: size exceeds index range");
			}
			matrix_.reserve(static_cast<index_type>(n));
		}

		void add(const index_type i, const index_type j, const value_type& value) {
			matrix_.add(i, j, value);
		}

		void set(const index_type i, const index_type j, const value_type& value) {
			matrix_.set(i, j, value);
		}

		matrix_type finalize_csr() const {
			matrix_type result = matrix_;
			result.finalize();
			return result;
		}

		matrix_type finalize_csc() const {
			matrix_type result = matrix_;
			result.to_csc();
			return result;
		}

	private:
		matrix_type matrix_;
	};
}

#endif
