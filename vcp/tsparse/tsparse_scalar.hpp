// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_SCALAR_HPP
#define VCP_TSPARSE_SCALAR_HPP

// ---------------------------------------------------------------------------
// Module scalar contract (SLU-GT1 D8)
//
// The tsparse / spmats_base generic code requires of T exactly the mats<T>
// operator set: arithmetic (+,-,*,/), comparisons, and ADL-resolved
// abs/log/sqrt/exp.  T/R-dependent values are never converted to double and
// no branch depends on what T is.
//
// Precision-guaranteed scalars: float / double / long double / kv::dd /
// kv::mpfr<N>.  kv::interval<TT> compiles and runs, but numerical success of
// the LU factorization is NOT guaranteed: gates are certified-only
// (!(x > tol) = "cannot certify x > tol"), so a pivot containing 0 is
// rejected BEFORE division and reported as numerical_singularity /
// zero_pivot.  On classes where interval Gaussian elimination succeeds
// (e.g. M-matrices) an enclosure is returned.
// ---------------------------------------------------------------------------

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

		// SLU-GT1 D4: std::numeric_limits<T> がプライマリテンプレート(未特殊化)
		// かどうかの検出。プライマリテンプレートは必ず is_specialized (= false) を
		// 定義する。一方、is_specialized を省略するユーザ特殊化(kv::dd,
		// kv::mpfr<N> は epsilon()/infinity() 等の関数のみ定義する)は定義上
		// 「特殊化済み」であり、メンバ不在は SFINAE で test_(...) 側に落ちる。
		// is_interval のような型族判定ではなく、numeric_limits の形状のみを見る。
		template <typename T>
		class numeric_limits_is_unspecialized {
			template <typename U>
			static typename std::enable_if<
				!std::numeric_limits<U>::is_specialized, char>::type test_(int);
			template <typename U>
			static long test_(...);
		public:
			static const bool value = (sizeof(test_<T>(0)) == sizeof(char));
		};

		template <typename T>
		inline T epsilon() {
			// numeric_limits の特殊化がない型(例: kv::interval)では T(0) を返す。
			// ゲート側は certified-only 意味論(!(x > eps) =「x が正と保証できない」)
			// になる。これは設計判断であり、プライマリテンプレートへの暗黙依存ではない。
			return numeric_limits_is_unspecialized<T>::value
				? T(0)
				: std::numeric_limits<T>::epsilon();
		}

		// SLU-GT1 D5: infinity<T>() was removed.  +infinity is not expressible
		// in the module scalar requirement set (arithmetic, certainly
		// comparisons, ADL abs/log/sqrt/exp), and for numeric_limits-
		// unspecialized types (e.g. kv::interval) the primary template would
		// silently return [0,0].  Field validity is expressed by bool flags
		// (converged / *_checked), never by numeric sentinels.

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
