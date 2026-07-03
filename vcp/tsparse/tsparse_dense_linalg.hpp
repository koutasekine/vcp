// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_DENSE_LINALG_HPP
#define VCP_TSPARSE_DENSE_LINALG_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
	namespace tsparse_dense_linalg {
		template <typename T>
		struct dense_lu_factorization {
			std::vector<std::vector<T> > lu;
			std::vector<std::size_t> pivots;
			bool singular;

			dense_lu_factorization() : singular(false) {}
		};

		template <typename T>
		struct dense_eigen_result {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			typedef std::complex<typename tsparse_scalar::real_type<T>::type> eigenvalue_type;
			std::vector<T> eigenvalues;
			std::vector<eigenvalue_type> complex_eigenvalues;
			std::vector<std::vector<T> > eigenvectors;
			std::vector<real_type> residuals;
			bool converged;
			std::size_t iterations;
			real_type residual_norm;

			// SLU-GT1 D5: residual_norm initialized to real_type(0); valid
			// only when `converged` is true (undefined otherwise -- do not read).
			dense_eigen_result()
				: converged(false), iterations(0), residual_norm(real_type(0)) {}
		};

		// SLU-GT1 D5: an empty input has no maximum; T(0) is returned as a
		// defensive placeholder, NOT a sentinel.  Callers must reject empty
		// inputs beforehand (currently this helper has no callers).
		template <typename T>
		T max_value(const std::vector<T>& values) {
			if (values.empty()) return T(0);
			T v(0);
			for (std::size_t i = 0; i < values.size(); i++) v = std::max(v, values[i]);
			return v;
		}

		template <typename T>
		void populate_real_complex_eigenvalues(dense_eigen_result<T>& result) {
			result.complex_eigenvalues.clear();
			result.complex_eigenvalues.reserve(result.eigenvalues.size());
			for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
				result.complex_eigenvalues.push_back(typename dense_eigen_result<T>::eigenvalue_type(result.eigenvalues[i]));
			}
		}

		template <typename T>
		void sort_eigenpairs(dense_eigen_result<T>& result) {
			std::vector<std::size_t> order(result.eigenvalues.size());
			for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
			std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
				return tsparse_scalar::real_part(result.eigenvalues[a]) <
				       tsparse_scalar::real_part(result.eigenvalues[b]);
			});
			std::vector<T> values(order.size());
			std::vector<typename dense_eigen_result<T>::eigenvalue_type> cvalues;
			if (result.complex_eigenvalues.size() == result.eigenvalues.size()) cvalues.resize(order.size());
			std::vector<typename dense_eigen_result<T>::real_type> residuals;
			if (result.residuals.size() == result.eigenvalues.size()) residuals.resize(order.size());
			std::vector<std::vector<T> > vectors;
			if (result.eigenvectors.size() == result.eigenvalues.size()) vectors.resize(order.size());
			for (std::size_t i = 0; i < order.size(); i++) {
				values[i] = result.eigenvalues[order[i]];
				if (!cvalues.empty()) cvalues[i] = result.complex_eigenvalues[order[i]];
				if (!residuals.empty()) residuals[i] = result.residuals[order[i]];
				if (!vectors.empty()) vectors[i] = result.eigenvectors[order[i]];
			}
			result.eigenvalues.swap(values);
			if (!cvalues.empty()) result.complex_eigenvalues.swap(cvalues);
			else result.complex_eigenvalues.clear();
			if (!residuals.empty()) result.residuals.swap(residuals);
			else result.residuals.clear();
			if (!vectors.empty()) result.eigenvectors.swap(vectors);
		}

		template <typename T>
		void select_eigenpairs(dense_eigen_result<T>& result, const std::size_t k, const bool largest) {
			if (k > result.eigenvalues.size()) {
				vcp::throw_error<vcp::invalid_argument>("tsparse_dense_linalg::select_eigenpairs: invalid k");
			}
			const eig_target target = largest ? eig_target::largest_algebraic : eig_target::smallest_algebraic;
			const std::vector<std::size_t> order =
				vcp::tsparse_eigen_selection::select_real_eigenpairs(
					result.eigenvalues, k, target, typename dense_eigen_result<T>::real_type(0));
			std::vector<T> values;
			std::vector<typename dense_eigen_result<T>::eigenvalue_type> cvalues;
			std::vector<typename dense_eigen_result<T>::real_type> residuals;
			std::vector<std::vector<T> > vectors;
			for (std::size_t i = 0; i < order.size(); i++) {
				const std::size_t j = order[i];
				values.push_back(result.eigenvalues[j]);
				if (result.complex_eigenvalues.size() == result.eigenvalues.size()) cvalues.push_back(result.complex_eigenvalues[j]);
				if (result.residuals.size() == result.eigenvalues.size()) residuals.push_back(result.residuals[j]);
				if (result.eigenvectors.size() == result.eigenvalues.size()) vectors.push_back(result.eigenvectors[j]);
			}
			result.eigenvalues.swap(values);
			if (result.complex_eigenvalues.size() == cvalues.size()) result.complex_eigenvalues.swap(cvalues);
			if (result.residuals.size() == residuals.size()) result.residuals.swap(residuals);
			if (result.eigenvectors.size() == vectors.size()) result.eigenvectors.swap(vectors);
		}

		template <typename T>
		bool is_dense_symmetric(const std::vector<std::vector<T> >& A, const typename tsparse_scalar::real_type<T>::type& tol) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const real_type tolerance = tol;
			for (std::size_t i = 0; i < A.size(); i++) {
				for (std::size_t j = i + 1; j < A.size(); j++) {
					if (tsparse_scalar::abs_value(A[i][j] - A[j][i]) > tolerance) return false;
				}
			}
			return true;
		}

		template <typename T>
		typename tsparse_scalar::real_type<T>::type lower_offdiag_norm(const std::vector<std::vector<T> >& A) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			real_type s(0);
			for (std::size_t i = 0; i < A.size(); i++) {
				for (std::size_t j = 0; j < i; j++) {
					const real_type v = tsparse_scalar::abs_value(A[i][j]);
					s += v * v;
				}
			}
			return tsparse_scalar::sqrt_value(s);
		}

		template <typename T>
		std::vector<std::vector<T> > matmul_dense(const std::vector<std::vector<T> >& A,
		                                          const std::vector<std::vector<T> >& B) {
			const std::size_t n = A.size();
			std::vector<std::vector<T> > C(n, std::vector<T>(n, T(0)));
			for (std::size_t i = 0; i < n; i++) {
				for (std::size_t k = 0; k < n; k++) {
					for (std::size_t j = 0; j < n; j++) C[i][j] += A[i][k] * B[k][j];
				}
			}
			return C;
		}

		template <typename T>
		std::vector<T> solve_upper_triangular(const std::vector<std::vector<T> >& R,
		                                      const std::vector<T>& rhs,
		                                      const std::size_t n) {
			std::vector<T> y(n, T(0));
			for (std::size_t kk = 0; kk < n; kk++) {
				const std::size_t i = n - 1 - kk;
				T s = rhs[i];
				for (std::size_t j = i + 1; j < n; j++) s -= R[i][j] * y[j];
				typedef typename tsparse_scalar::real_type<T>::type real_type;
				if (!(tsparse_scalar::abs_value(R[i][i]) > vcp::tsparse_scalar::epsilon<real_type>())) {
					vcp::throw_error<vcp::numerical_error>("tsparse_dense_linalg::solve_upper_triangular: singular matrix");
				}
				y[i] = s / R[i][i];
			}
			return y;
		}

		template <typename T, typename Tol>
		dense_lu_factorization<T> dense_lu_factor(std::vector<std::vector<T> > A, const Tol& tol) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const std::size_t n = A.size();
			const real_type tolerance(tol);
			dense_lu_factorization<T> factor;
			factor.lu = A;
			factor.pivots.resize(n);
			factor.singular = false;
			for (std::size_t k = 0; k < n; k++) {
				std::size_t pivot = k;
				real_type pivot_abs = tsparse_scalar::abs_value(factor.lu[k][k]);
				for (std::size_t i = k + 1; i < n; i++) {
					const real_type v = tsparse_scalar::abs_value(factor.lu[i][k]);
					if (v > pivot_abs) {
						pivot_abs = v;
						pivot = i;
					}
				}
				factor.pivots[k] = pivot;
				if (pivot_abs <= tolerance) {
					factor.singular = true;
					return factor;
				}
				if (pivot != k) std::swap(factor.lu[pivot], factor.lu[k]);
				for (std::size_t i = k + 1; i < n; i++) {
					factor.lu[i][k] /= factor.lu[k][k];
					for (std::size_t j = k + 1; j < n; j++) factor.lu[i][j] -= factor.lu[i][k] * factor.lu[k][j];
				}
			}
			return factor;
		}

		template <typename T>
		std::vector<T> dense_lu_solve(const dense_lu_factorization<T>& factor, std::vector<T> rhs) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			if (factor.singular) {
				vcp::throw_error<vcp::numerical_error>("tsparse_dense_linalg::dense_lu_solve: singular matrix");
			}
			const std::size_t n = factor.lu.size();
			for (std::size_t k = 0; k < n; k++) {
				const std::size_t pivot = factor.pivots[k];
				if (pivot != k) std::swap(rhs[pivot], rhs[k]);
			}
			for (std::size_t i = 0; i < n; i++) {
				for (std::size_t j = 0; j < i; j++) rhs[i] -= factor.lu[i][j] * rhs[j];
			}
			for (std::size_t kk = 0; kk < n; kk++) {
				const std::size_t i = n - 1 - kk;
				for (std::size_t j = i + 1; j < n; j++) rhs[i] -= factor.lu[i][j] * rhs[j];
				if (!(tsparse_scalar::abs_value(factor.lu[i][i]) > vcp::tsparse_scalar::epsilon<real_type>())) {
					vcp::throw_error<vcp::numerical_error>("tsparse_dense_linalg::dense_lu_solve: singular matrix");
				}
				rhs[i] /= factor.lu[i][i];
			}
			return rhs;
		}

		template <typename T>
		std::vector<T> solve_dense_gaussian(std::vector<std::vector<T> > A, std::vector<T> b) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const std::size_t n = A.size();
			for (std::size_t k = 0; k < n; k++) {
				std::size_t pivot = k;
				real_type pivot_abs = tsparse_scalar::abs_value(A[k][k]);
				for (std::size_t i = k + 1; i < n; i++) {
					const real_type v = tsparse_scalar::abs_value(A[i][k]);
					if (v > pivot_abs) {
						pivot_abs = v;
						pivot = i;
					}
				}
				if (!(pivot_abs > vcp::tsparse_scalar::epsilon<real_type>())) {
					vcp::throw_error<vcp::numerical_error>("tsparse_dense_linalg::solve_dense_gaussian: singular matrix");
				}
				if (pivot != k) {
					std::swap(A[pivot], A[k]);
					std::swap(b[pivot], b[k]);
				}
				for (std::size_t i = k + 1; i < n; i++) {
					const T factor = A[i][k] / A[k][k];
					A[i][k] = T(0);
					for (std::size_t j = k + 1; j < n; j++) A[i][j] -= factor * A[k][j];
					b[i] -= factor * b[k];
				}
			}
			std::vector<T> x(n, T(0));
			for (std::size_t ii = 0; ii < n; ii++) {
				const std::size_t i = n - 1 - ii;
				T s = b[i];
				for (std::size_t j = i + 1; j < n; j++) s -= A[i][j] * x[j];
				x[i] = s / A[i][i];
			}
			return x;
		}

		template <typename T>
		typename tsparse_scalar::real_type<T>::type dense_eigenpair_residual_norm_value(const std::vector<std::vector<T> >& A,
		                                                                               const T& eigenvalue,
		                                                                               const std::vector<T>& eigenvector) {
			std::vector<T> r(eigenvector.size(), T(0));
			for (std::size_t i = 0; i < A.size(); i++) {
				for (std::size_t j = 0; j < A[i].size(); j++) r[i] += A[i][j] * eigenvector[j];
				r[i] -= eigenvalue * eigenvector[i];
			}
			return tsparse_scalar::real_norm_value(r);
		}

		template <typename T>
		typename tsparse_scalar::real_type<T>::type max_dense_eigenpair_residual_value(const std::vector<std::vector<T> >& A,
		                                                                              const std::vector<T>& eigenvalues,
		                                                                              const std::vector<std::vector<T> >& eigenvectors) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			// SLU-GT1 D5: defensive branch only.  Both callers (qr_eig_dense)
			// build eigenvectors with eigenvalues.size() entries before calling,
			// so this branch is structurally unreachable; real_type(0) is a
			// placeholder, not a sentinel.  New callers must guarantee
			// non-empty, size-matched inputs.
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return real_type(0);
			real_type maximum(0);
			for (std::size_t p = 0; p < eigenvalues.size(); p++) {
				const real_type residual = dense_eigenpair_residual_norm_value(A, eigenvalues[p], eigenvectors[p]);
				if (residual > maximum) maximum = residual;
			}
			return maximum;
		}

		template <typename T>
		std::vector<typename tsparse_scalar::real_type<T>::type> dense_eigenpair_residuals(const std::vector<std::vector<T> >& A,
		                                                                                  const std::vector<T>& eigenvalues,
		                                                                                  const std::vector<std::vector<T> >& eigenvectors) {
			std::vector<typename tsparse_scalar::real_type<T>::type> residuals;
			if (eigenvalues.empty() || eigenvectors.size() != eigenvalues.size()) return residuals;
			residuals.reserve(eigenvalues.size());
			for (std::size_t p = 0; p < eigenvalues.size(); p++) {
				residuals.push_back(dense_eigenpair_residual_norm_value(A, eigenvalues[p], eigenvectors[p]));
			}
			return residuals;
		}

		template <typename T>
		std::vector<T> dense_eigenvector_inverse_iteration(const std::vector<std::vector<T> >& A, const T& lambda) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const std::size_t n = A.size();
			std::vector<T> y(n, T(1));
			const T shift = lambda + T(tsparse_scalar::decimal_power_negative<real_type>(10));
			for (std::size_t iter = 0; iter < 8; iter++) {
				std::vector<std::vector<T> > M = A;
				for (std::size_t i = 0; i < n; i++) M[i][i] -= shift;
				try {
					y = solve_dense_gaussian(M, y);
				}
				catch (const std::exception&) {
					break;
				}
				const real_type ny = tsparse_scalar::real_norm_value(y);
				if (!(ny > vcp::tsparse_scalar::epsilon<real_type>())) break;
				for (std::size_t i = 0; i < n; i++) y[i] /= T(ny);
			}
			return y;
		}

		template <typename T>
		std::vector<T> dense_eigenvector_inverse_iteration_from(
			const std::vector<std::vector<T> >& A,
			const T& lambda,
			const std::vector<T>& y_init)
		{
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const std::size_t n = A.size();
			std::vector<T> y = y_init;
			if (y.size() != n) y.assign(n, T(1));
			const T shift = lambda + T(tsparse_scalar::decimal_power_negative<real_type>(10));
			for (std::size_t iter = 0; iter < 8; iter++) {
				std::vector<std::vector<T> > M = A;
				for (std::size_t i = 0; i < n; i++) M[i][i] -= shift;
				try {
					y = solve_dense_gaussian(M, y);
				}
				catch (const std::exception&) {
					break;
				}
				const real_type ny = tsparse_scalar::real_norm_value(y);
				if (!(ny > vcp::tsparse_scalar::epsilon<real_type>())) break;
				for (std::size_t i = 0; i < n; i++) y[i] /= T(ny);
			}
			return y;
		}

		template <typename T>
		std::vector<T> small_real_eigenvalues(const std::vector<std::vector<T> >& A) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const std::size_t n = A.size();
			std::vector<T> values;
			if (n == 0) return values;
			if (n == 1) {
				values.push_back(A[0][0]);
				return values;
			}
			if (n == 2) {
				const real_type tr = tsparse_scalar::real_part(A[0][0] + A[1][1]);
				const real_type det = tsparse_scalar::real_part(A[0][0] * A[1][1] - A[0][1] * A[1][0]);
				const real_type disc = tr * tr - real_type(4) * det;
				if (disc < -tsparse_scalar::decimal_power_negative<real_type>(14)) {
					vcp::throw_error<vcp::domain_error>("tsparse_dense_linalg::dense_qr: real API cannot represent complex eigenvalues");
				}
				if (disc < real_type(0)) {
					values.push_back(T(tr / real_type(2)));
					values.push_back(T(tr / real_type(2)));
				}
				else {
					const real_type s = tsparse_scalar::sqrt_value(disc);
					values.push_back(T((tr - s) / real_type(2)));
					values.push_back(T((tr + s) / real_type(2)));
				}
				return values;
			}
			const real_type a00 = tsparse_scalar::real_part(A[0][0]);
			const real_type a01 = tsparse_scalar::real_part(A[0][1]);
			const real_type a02 = tsparse_scalar::real_part(A[0][2]);
			const real_type a10 = tsparse_scalar::real_part(A[1][0]);
			const real_type a11 = tsparse_scalar::real_part(A[1][1]);
			const real_type a12 = tsparse_scalar::real_part(A[1][2]);
			const real_type a20 = tsparse_scalar::real_part(A[2][0]);
			const real_type a21 = tsparse_scalar::real_part(A[2][1]);
			const real_type a22 = tsparse_scalar::real_part(A[2][2]);
			const real_type tr = a00 + a11 + a22;
			const real_type c2 = a00 * a11 + a00 * a22 + a11 * a22 - a01 * a10 - a02 * a20 - a12 * a21;
			const real_type det = a00 * (a11 * a22 - a12 * a21) - a01 * (a10 * a22 - a12 * a20) + a02 * (a10 * a21 - a11 * a20);
			const real_type p = c2 - tr * tr / real_type(3);
			const real_type q = -real_type(2) * tr * tr * tr / real_type(27) + tr * c2 / real_type(3) - det;
			using std::acos;
			using std::cos;
			const real_type pi = acos(real_type(-1));
			if (tsparse_scalar::abs_value(p) <= tsparse_scalar::decimal_power_negative<real_type>(30)) {
				const real_type root = tr / real_type(3);
				values.push_back(T(root));
				values.push_back(T(root));
				values.push_back(T(root));
				return values;
			}
			const real_type discr = q * q / real_type(4) + p * p * p / real_type(27);
			if (discr <= real_type(0)) {
				real_type arg = (real_type(3) * q / (real_type(2) * p)) * tsparse_scalar::sqrt_value(-real_type(3) / p);
				if (arg < real_type(-1)) arg = real_type(-1);
				if (arg > real_type(1)) arg = real_type(1);
				const real_type phi = acos(arg) / real_type(3);
				const real_type scale = real_type(2) * tsparse_scalar::sqrt_value(-p / real_type(3));
				for (std::size_t k = 0; k < 3; k++) {
					values.push_back(T(scale * cos(phi - real_type(2) * pi * real_type(k) / real_type(3)) + tr / real_type(3)));
				}
			}
			else {
				vcp::throw_error<vcp::domain_error>("tsparse_dense_linalg::dense_qr: real API cannot represent complex eigenvalues");
			}
			return values;
		}

		template <typename T>
		dense_eigen_result<T> jacobi_eig_dense(std::vector<std::vector<T> > A,
		                                       const std::size_t max_iter,
		                                       const typename tsparse_scalar::real_type<T>::type& tol) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const real_type tolerance = tol;
			const std::size_t n = A.size();
			dense_eigen_result<T> result;
			result.eigenvectors.assign(n, std::vector<T>(n, T(0)));
			for (std::size_t i = 0; i < n; i++) result.eigenvectors[i][i] = T(1);
			result.converged = false;
			result.iterations = 0;
			result.residual_norm = real_type(0);
			for (std::size_t iter = 1; iter <= max_iter; iter++) {
				std::size_t p = 0, q = 0;
				real_type max_off(0);
				// SLU-GT1.1: convergence must be CERTIFIED -- every off-diagonal
				// certainly <= tolerance.  The old `max_off <= tolerance` check is
				// a fake-convergence hole for intervals: when every off-diagonal
				// straddles 0, no certainly-> comparison fires, max_off stays
				// [0,0], and the certainly-<= test wrongly declares convergence.
				// For double the two are decision-equivalent
				// (max <= tol <=> all aij <= tol).
				bool all_small = true;
				for (std::size_t i = 0; i < n; i++) {
					for (std::size_t j = i + 1; j < n; j++) {
						const real_type aij = tsparse_scalar::abs_value(A[i][j]);
						if (aij > max_off) {
							max_off = aij;
							p = i;
							q = j;
						}
						if (!(aij <= tolerance)) all_small = false;
					}
				}
				result.residual_norm = max_off;
				result.iterations = iter - 1;
				if (all_small) {
					result.converged = true;
					break;
				}
				// SLU-GT1.1 F-7: certified-only -- if the pivot magnitude cannot
				// be certified positive, the rotation's division by apq and the
				// sqrt arguments below cannot be certified either.  Do not fake
				// progress: stop with converged=false (honest non-convergence).
				// For double, !all_small implies some aij > tolerance >= 0, so
				// max_off > 0 and this branch is unreachable.
				if (!(max_off > real_type(0))) break;
				const real_type app = tsparse_scalar::real_part(A[p][p]);
				const real_type aqq = tsparse_scalar::real_part(A[q][q]);
				const real_type apq = tsparse_scalar::real_part(A[p][q]);
				const real_type tau = (aqq - app) / (real_type(2) * apq);
				const real_type sign = tau >= real_type(0) ? real_type(1) : real_type(-1);
				// SLU-GT1.1 F-7: square via abs so the sqrt arguments are
				// certified >= 1 for interval real_type; |x|*|x| is bit-identical
				// to x*x for double.
				const real_type atau = tsparse_scalar::abs_value(tau);
				const real_type t = sign / (atau + tsparse_scalar::sqrt_value(real_type(1) + atau * atau));
				const real_type at = tsparse_scalar::abs_value(t);
				const real_type c = real_type(1) / tsparse_scalar::sqrt_value(real_type(1) + at * at);
				const real_type s = t * c;
				for (std::size_t k = 0; k < n; k++) {
					if (k != p && k != q) {
						const T akp = A[k][p];
						const T akq = A[k][q];
						A[k][p] = T(c) * akp - T(s) * akq;
						A[p][k] = A[k][p];
						A[k][q] = T(s) * akp + T(c) * akq;
						A[q][k] = A[k][q];
					}
				}
				const T new_app = T(c * c) * A[p][p] - T(real_type(2) * c * s) * A[p][q] + T(s * s) * A[q][q];
				const T new_aqq = T(s * s) * A[p][p] + T(real_type(2) * c * s) * A[p][q] + T(c * c) * A[q][q];
				A[p][p] = new_app;
				A[q][q] = new_aqq;
				A[p][q] = T(0);
				A[q][p] = T(0);
				for (std::size_t k = 0; k < n; k++) {
					const T vkp = result.eigenvectors[k][p];
					const T vkq = result.eigenvectors[k][q];
					result.eigenvectors[k][p] = T(c) * vkp - T(s) * vkq;
					result.eigenvectors[k][q] = T(s) * vkp + T(c) * vkq;
				}
				result.iterations = iter;
			}
			result.eigenvalues.resize(n);
			for (std::size_t i = 0; i < n; i++) result.eigenvalues[i] = A[i][i];
			std::vector<std::vector<T> > vectors(n, std::vector<T>(n, T(0)));
			for (std::size_t j = 0; j < n; j++) {
				for (std::size_t i = 0; i < n; i++) vectors[j][i] = result.eigenvectors[i][j];
			}
			result.eigenvectors.swap(vectors);
			populate_real_complex_eigenvalues(result);
			sort_eigenpairs(result);
			return result;
		}

		template <typename T>
		dense_eigen_result<T> qr_eig_dense(std::vector<std::vector<T> > A,
		                                   const std::size_t max_iter,
		                                   const typename tsparse_scalar::real_type<T>::type& tol) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const real_type tolerance = tol;
			const std::size_t n = A.size();
			const std::vector<std::vector<T> > original = A;
			dense_eigen_result<T> result;
			result.converged = false;
			result.iterations = 0;
			// SLU-GT1 D5: real_type(0), not an infinity sentinel; `converged`
			// is the validity witness for residual_norm.
			result.residual_norm = real_type(0);
			result.eigenvectors.clear();
			if (n <= 3) {
				result.eigenvalues = small_real_eigenvalues(A);
				for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
					result.eigenvectors.push_back(dense_eigenvector_inverse_iteration(original, result.eigenvalues[i]));
				}
				result.converged = true;
				result.iterations = 0;
				result.residuals = dense_eigenpair_residuals(original, result.eigenvalues, result.eigenvectors);
				const real_type residual_value = max_dense_eigenpair_residual_value(original, result.eigenvalues, result.eigenvectors);
				result.residual_norm = residual_value;
				populate_real_complex_eigenvalues(result);
				sort_eigenpairs(result);
				return result;
			}
			for (std::size_t iter = 1; iter <= max_iter; iter++) {
				const T shift = A[n - 1][n - 1];
				for (std::size_t i = 0; i < n; i++) A[i][i] -= shift;
				std::vector<std::vector<T> > q(n, std::vector<T>(n, T(0)));
				std::vector<std::vector<T> > r(n, std::vector<T>(n, T(0)));
				for (std::size_t j = 0; j < n; j++) {
					std::vector<T> v(n);
					for (std::size_t i = 0; i < n; i++) v[i] = A[i][j];
					for (std::size_t k = 0; k < j; k++) {
						T rij = T(0);
						for (std::size_t i = 0; i < n; i++) rij += q[i][k] * v[i];
						r[k][j] = rij;
						for (std::size_t i = 0; i < n; i++) v[i] -= rij * q[i][k];
					}
					real_type nv = tsparse_scalar::real_norm_value(v);
					if (!(nv > vcp::tsparse_scalar::epsilon<real_type>())) {
						std::fill(v.begin(), v.end(), T(0));
						v[j] = T(1);
						for (std::size_t k = 0; k < j; k++) {
							T rij = T(0);
							for (std::size_t i = 0; i < n; i++) rij += q[i][k] * v[i];
							for (std::size_t i = 0; i < n; i++) v[i] -= rij * q[i][k];
						}
						nv = tsparse_scalar::real_norm_value(v);
						if (!(nv > vcp::tsparse_scalar::epsilon<real_type>())) continue;
					}
					r[j][j] = T(nv);
					for (std::size_t i = 0; i < n; i++) q[i][j] = v[i] / T(nv);
				}
				A = matmul_dense(r, q);
				for (std::size_t i = 0; i < n; i++) A[i][i] += shift;
				const real_type offdiag = lower_offdiag_norm(A);
				result.residual_norm = offdiag;
				result.iterations = iter;
				if (offdiag <= tolerance) {
					result.converged = true;
					break;
				}
			}
			result.eigenvalues.resize(n);
			for (std::size_t i = 0; i < n; i++) result.eigenvalues[i] = A[i][i];
			if (result.converged) {
				for (std::size_t i = 0; i < result.eigenvalues.size(); i++) {
					result.eigenvectors.push_back(dense_eigenvector_inverse_iteration(original, result.eigenvalues[i]));
				}
				result.residuals = dense_eigenpair_residuals(original, result.eigenvalues, result.eigenvectors);
				const real_type residual_value = max_dense_eigenpair_residual_value(original, result.eigenvalues, result.eigenvectors);
				result.residual_norm = residual_value;
				const real_type residual_limit = tolerance > tolerance * real_type(n) * real_type(10)
					? tolerance : tolerance * real_type(n) * real_type(10);
				result.converged = residual_value <= residual_limit;
			}
			populate_real_complex_eigenvalues(result);
			sort_eigenpairs(result);
			return result;
		}

		template <typename T>
		dense_eigen_result<T> dense_eig(std::vector<std::vector<T> > dense, const bool symmetric,
		                                const std::size_t max_iter,
		                                const typename tsparse_scalar::real_type<T>::type& tol) {
			if (symmetric) return jacobi_eig_dense(dense, max_iter, tol);
			return qr_eig_dense(dense, max_iter, tol);
		}

		template <typename T>
		std::vector<std::vector<T> > lift_ritz_vectors(const std::vector<std::vector<T> >& V,
		                                               const std::vector<std::vector<T> >& small_vectors,
		                                               const std::size_t n) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			std::vector<std::vector<T> > vectors;
			vectors.reserve(small_vectors.size());
			for (std::size_t p = 0; p < small_vectors.size(); p++) {
				std::vector<T> v(n, T(0));
				for (std::size_t j = 0; j < small_vectors[p].size() && j < V.size(); j++) {
					for (std::size_t i = 0; i < n; i++) v[i] += V[j][i] * small_vectors[p][j];
				}
				const real_type nv = tsparse_scalar::real_norm_value(v);
				if (nv > vcp::tsparse_scalar::epsilon<real_type>()) {
					for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
				}
				vectors.push_back(v);
			}
			return vectors;
		}

		// lift_ritz_vectors_checked: dimension-safe variant.
		// Throws vcp::dimension_error if any basis vector has wrong size or if
		// coefficient count != V.size().
		// Throws vcp::numerical_error if the resulting Ritz vector has near-zero norm.
		// Used by Phase 8 projected eigensolver.
		template <typename T>
		std::vector<T> lift_ritz_vector_checked(const std::vector<std::vector<T> >& V,
		                                        const std::vector<T>& y,
		                                        const std::size_t n,
		                                        const typename tsparse_scalar::real_type<T>::type& zero_tol) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			if (y.size() != V.size()) {
				vcp::throw_error<vcp::dimension_error>(
					"tsparse_dense_linalg::lift_ritz_vector_checked: "
					"coefficient size does not match basis size");
			}
			std::vector<T> v(n, T(0));
			for (std::size_t j = 0; j < V.size(); j++) {
				if (V[j].size() != n) {
					vcp::throw_error<vcp::dimension_error>(
						"tsparse_dense_linalg::lift_ritz_vector_checked: "
						"basis vector dimension does not match n");
				}
				for (std::size_t i = 0; i < n; i++) v[i] += V[j][i] * y[j];
			}
			const real_type nv = tsparse_scalar::real_norm_value(v);
			if (nv <= zero_tol) {
				vcp::throw_error<vcp::numerical_error>(
					"tsparse_dense_linalg::lift_ritz_vector_checked: "
					"Ritz vector has near-zero norm");
			}
			for (std::size_t i = 0; i < n; i++) v[i] /= T(nv);
			return v;
		}
	}
}

#endif
