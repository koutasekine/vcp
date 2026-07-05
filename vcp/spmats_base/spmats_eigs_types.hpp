// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATS_EIGS_TYPES_HPP
#define VCP_SPMATS_EIGS_TYPES_HPP

#include <complex>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_solvers.hpp>
#include <vcp/tsparse/tsparse_sparse_lu.hpp>

namespace vcp {

	// -----------------------------------------------------------------------
	// Type traits for scalar classification
	// -----------------------------------------------------------------------

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

	// -----------------------------------------------------------------------
	// Linear solver method / preconditioner type
	// -----------------------------------------------------------------------

	enum class linear_solver_method {
		jacobi,
		gauss_seidel,
		conjugate_gradient,
		bicgstab,
		gmres,
		sparse_lu
	};

	enum class preconditioner_type {
		none,
		jacobi
	};

	// eig_method, eigs_target, generalized_eig_method, matrix_structure_hint,
	// orthogonalization_method, eig_solver_method, eig_target
	// are all defined in <vcp/tsparse/tsparse_eigs.hpp> (included above).

	// -----------------------------------------------------------------------
	// linear_solve_options<T>
	// -----------------------------------------------------------------------

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
		vcp::sparse_lu_options<T> sparse_lu;

		linear_solve_options()
			: method(linear_solver_method::conjugate_gradient), max_iter(1000),
			  tol(vcp::tsparse_scalar::decimal_power_negative<real_type>(12)),
			  check_symmetric(true), use_relative_residual(true), restart(30),
			  preconditioner(preconditioner_type::none) {}
	};

	// -----------------------------------------------------------------------
	// eig_shift_invert_solver (E-A1)
	//
	// Inner linear solver used by the shift-invert eigensolver paths
	// (standard shift_invert_lanczos / shift_invert_arnoldi and the
	// generalized shift-invert operator).
	// -----------------------------------------------------------------------

	enum class eig_shift_invert_solver {
		sparse_lu,    // default: direct solve (factorize (A - sigma*B) once,
		              // solve repeatedly; solve-time IR follows the
		              // eig_options::shift_invert_lu options)
		ilu0_gmres    // compatibility opt-in: legacy ILU(0)-preconditioned
		              // GMRES inner solve (byte-identical legacy path)
	};

	// -----------------------------------------------------------------------
	// eig_options<T>
	// -----------------------------------------------------------------------

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
		// E-A1: inner solver selection for the shift-invert paths
		// (default: sparse_lu direct solve; D-1).
		eig_shift_invert_solver shift_invert_solver;
		// E-A1: sparse LU options passed through to the shift-invert LU
		// operator (read only when shift_invert_solver == sparse_lu).
		// Default construction = auto_select / threshold_partial /
		// iterative_refinement = true (D-3).
		vcp::sparse_lu_options<T> shift_invert_lu;

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
			  compute_residual_history(false),
			  shift_invert_solver(eig_shift_invert_solver::sparse_lu) {}
	};

	// -----------------------------------------------------------------------
	// linear_solve_result<T>
	// -----------------------------------------------------------------------

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

		// SLU-GT1 D5: residual fields are initialized to real_type(0), NOT an
		// infinity sentinel.  `converged` is the ONLY validity witness: when
		// converged == false the residual_norm / absolute_ / relative_ fields
		// are undefined values and must not be read.
		linear_solve_result()
			: converged(false), iterations(0), residual_norm(real_type(0)),
			  absolute_residual_norm(real_type(0)),
			  relative_residual_norm(real_type(0)),
			  initial_residual_norm(real_type(0)), method(linear_solver_method::conjugate_gradient) {}
	};

	// -----------------------------------------------------------------------
	// eig_result<T>
	//
	// Meaning of converged == true (EIG-0 contract, C-1/C-2/C-4):
	//   (C-1) every returned eigenpair satisfies the residual acceptance test
	//         re-evaluated with the EXACT operator at termination, and
	//   (C-2) at termination no unconverged Ritz candidate was certainly
	//         visible inside (more target-preferred than) the returned set.
	//   (C-4) converged == true is NOT a completeness guarantee: Krylov
	//         subspace methods cannot see eigenspaces orthogonal to the start
	//         vector, so "the k target-side eigenvalues were all found" cannot
	//         be certified by any such solver.  Rigorous enclosure /
	//         completeness belongs to a future verification layer
	//         (Lehmann-Goerisch line), not to this flag.
	// -----------------------------------------------------------------------

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
		// E-A1 E4: the inner_* / factorization_* diagnostics below carry
		// per-solver meanings on the shift-invert paths, depending on
		// eig_options::shift_invert_solver:
		//
		//   field                      | ilu0_gmres (legacy)     | sparse_lu (default)
		//   ---------------------------+-------------------------+--------------------------------
		//   inner_iterations           | total GMRES iterations  | total solve-time IR iterations
		//                              |                         | (0 if IR off)
		//   inner_failure_count        | GMRES non-convergences  | solves whose IR did not reach
		//                              |                         | the IR tolerance
		//   inner_residual_norm        | max final GMRES residual| IR final_residual of the LAST
		//                              |                         | solve (0 if IR off)
		//   inner_failure_reason       | GMRES failure reason    | normally empty (factorization
		//                              |                         | failure is reported via
		//                              |                         | status/failure_reason instead)
		//   factorization_diagnostics  | ILU(0) diagnostics      | sparse_lu_status string +
		//                              |                         | n/nnz summary
		//   factorization_zero_pivots  | ILU(0) zero pivot count | sparse_lu_info::
		//                              |                         | within_panel_zero_pivot_count
		//                              |                         | (populated on the supernodal
		//                              |                         | path; 0 on baseline GP where
		//                              |                         | zero pivots abort via status)
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
			  // SLU-GT1 D5: initialized to real_type(0); valid only when
			  // `converged` (or an explicit residual computation) sets them.
			  residual_norm_absolute(real_type(0)),
			  residual_norm_relative(real_type(0)),
			  residual_history_absolute(), residual_history_relative(),
			  method(eig_solver_method::lanczos), status(), message(), failure_reason(),
			  breakdown_reason(), used_method(), used_orthogonalization(),
			  used_dense_fallback(false), used_shift_invert(false), used_generalized_operator(false),
			  used_subspace_dim(0), inner_iterations(0), inner_failure_count(0),
			  inner_residual_norm(real_type(0)), inner_failure_reason(),
			  factorization_diagnostics(), factorization_zero_pivots(0) {}
	};

} // namespace vcp

#endif // VCP_SPMATS_EIGS_TYPES_HPP
