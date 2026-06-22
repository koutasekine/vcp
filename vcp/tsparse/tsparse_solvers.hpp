// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_SOLVERS_HPP
#define VCP_TSPARSE_SOLVERS_HPP

#include <algorithm>
#include <cstddef>
#include <vector>

#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
	namespace tsparse_solvers {
		template <typename T>
		struct residual_control {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			real_type rhs_norm;
			real_type denominator;
			real_type threshold;
			bool relative;

			residual_control()
				: rhs_norm(real_type(0)), denominator(real_type(1)), threshold(real_type(0)), relative(true) {}
		};

		template <typename T>
		residual_control<T> make_residual_control(const std::vector<T>& b,
		                                          const typename tsparse_scalar::real_type<T>::type& tol,
		                                          const bool use_relative_residual) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			residual_control<T> control;
			control.rhs_norm = tsparse_scalar::real_norm_value(b);
			control.denominator = control.rhs_norm > real_type(1) ? control.rhs_norm : real_type(1);
			control.relative = use_relative_residual;
			control.threshold = use_relative_residual ? tol * control.denominator : tol;
			return control;
		}

		template <class SparseMatrix, typename T>
		typename tsparse_scalar::real_type<T>::type residual_norm_value(const SparseMatrix& A, const std::vector<T>& x, const std::vector<T>& b) {
			std::vector<T> r = A.mul_vec(x);
			for (std::size_t i = 0; i < r.size(); i++) r[i] -= b[i];
			return tsparse_scalar::real_norm_value(r);
		}

		template <class Result, class SparseMatrix, typename T>
		void set_linear_residual_fields(Result& result, const SparseMatrix& A, const std::vector<T>& b) {
			const typename tsparse_scalar::real_type<T>::type absolute = residual_norm_value(A, result.x, b);
			const typename tsparse_scalar::real_type<T>::type initial = tsparse_scalar::real_norm_value(b);
			const typename tsparse_scalar::real_type<T>::type denom = initial > typename tsparse_scalar::real_type<T>::type(1)
				? initial : typename tsparse_scalar::real_type<T>::type(1);
			result.absolute_residual_norm = absolute;
			result.initial_residual_norm = initial;
			result.relative_residual_norm = result.absolute_residual_norm / denom;
			result.residual_norm = result.absolute_residual_norm;
			result.solution = result.x;
		}
	}
}

#endif
