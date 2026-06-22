// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_SCALAR_HPP
#define VCP_TSPARSE_SCALAR_HPP

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <vector>

#include <vcp/error.hpp>

namespace vcp {
	namespace tsparse_scalar {
		template <typename T> struct is_complex : std::false_type {};
		template <typename T> struct is_complex<std::complex<T> > : std::true_type {};

		template <typename T> struct real_type {
			typedef T type;
		};
		template <typename T> struct real_type<std::complex<T> > {
			typedef T type;
		};

		template <typename T>
		inline typename std::enable_if<!is_complex<T>::value, typename real_type<T>::type>::type real_part(const T& x) {
			return x;
		}

		template <typename T>
		inline typename std::enable_if<is_complex<T>::value, typename real_type<T>::type>::type real_part(const T& x) {
			return x.real();
		}

		template <typename T>
		inline typename real_type<T>::type abs_value(const T& x) {
			using std::abs;
			return abs(x);
		}

		template <typename T>
		inline T sqrt_value(const T& x) {
			using std::sqrt;
			return sqrt(x);
		}

		template <typename T>
		inline T hypot_value(const T& x, const T& y) {
			return sqrt_value(x * x + y * y);
		}

		template <typename T>
		inline T epsilon() {
			return std::numeric_limits<T>::epsilon();
		}

		template <typename T>
		inline T infinity() {
			return std::numeric_limits<T>::infinity();
		}

		template <typename T>
		inline typename std::enable_if<!std::is_floating_point<T>::value, bool>::type is_finite(const T&) {
			return true;
		}

		template <typename T>
		inline typename std::enable_if<std::is_floating_point<T>::value, bool>::type is_finite(const T& x) {
			return std::isfinite(x);
		}

		template <typename T>
		inline T decimal_power_negative(const unsigned int exponent) {
			T value(1);
			for (unsigned int i = 0; i < exponent; i++) value /= T(10);
			return value;
		}

		template <typename T>
		inline typename std::enable_if<!is_complex<T>::value, T>::type conjugate_if_needed(const T& x) {
			return x;
		}

		template <typename T>
		inline typename std::enable_if<is_complex<T>::value, T>::type conjugate_if_needed(const T& x) {
			using std::conj;
			return conj(x);
		}

		template <typename T>
		T real_dot(const std::vector<T>& x, const std::vector<T>& y) {
			if (x.size() != y.size()) {
				vcp::throw_error<vcp::dimension_error>("tsparse_scalar::real_dot: dimension mismatch");
			}
			T sum = T(0);
			for (std::size_t i = 0; i < x.size(); i++) sum += x[i] * y[i];
			return sum;
		}

		template <typename T>
		T hermitian_dot(const std::vector<T>& x, const std::vector<T>& y) {
			if (x.size() != y.size()) {
				vcp::throw_error<vcp::dimension_error>("tsparse_scalar::hermitian_dot: dimension mismatch");
			}
			T sum = T(0);
			for (std::size_t i = 0; i < x.size(); i++) sum += conjugate_if_needed(x[i]) * y[i];
			return sum;
		}

		template <typename T>
		typename real_type<T>::type real_dot_value(const std::vector<T>& x, const std::vector<T>& y) {
			return real_part(real_dot(x, y));
		}

		template <typename T>
		typename real_type<T>::type hermitian_norm_value(const std::vector<T>& x) {
			const T v = hermitian_dot(x, x);
			const typename real_type<T>::type r = real_part(v);
			return sqrt_value(r > typename real_type<T>::type(0) ? r : typename real_type<T>::type(0));
		}

		template <typename T>
		typename real_type<T>::type real_norm_value(const std::vector<T>& x) {
			const typename real_type<T>::type v = real_dot_value(x, x);
			return sqrt_value(v > typename real_type<T>::type(0) ? v : typename real_type<T>::type(0));
		}
	}
}

#endif
