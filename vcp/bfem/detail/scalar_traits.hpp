// vcp/bfem/detail/scalar_traits.hpp
// L4 seam alignment (design L4_concept_design_v0.1 section 10, item 3):
// mechanization of the C-1 scalar contract as a static_assert battery.
//
// C-1 requires of a space scalar T:
//  - ring operations +, -, * (and unary -), plus / inside the declared
//    division budget (C-3.1);
//  - construction of small exact constants (T(int), default T());
//  - convert_traits<T>: the enclose-once rational -> T conversion entry
//    point exists with the from_rational(bigint, bigint) -> T signature;
//  - operation counting must remain possible: the battery below checks
//    EXACTLY the countable operation set (the test counting_scalar
//    implements precisely these operators), so any T passing this battery
//    keeps the counting acceptance path (C-3.1) available.
//
// bfem_scalar_traits<T>::require() is called at every space construction
// point (fe_space / rt_space / broken_space / vfe_space / c1_space). It is
// an ADDITIVE compile-time check only: no code is generated, no behavior
// changes, and every scalar type accepted before L4 passes the battery.

#ifndef VCP_BFEM_DETAIL_SCALAR_TRAITS_HPP
#define VCP_BFEM_DETAIL_SCALAR_TRAITS_HPP

#include <type_traits>
#include <utility>

#include <vcp/bfem/rational.hpp>
#include <vcp/bfem/convert_traits.hpp>

namespace vcp {
namespace bfem {

namespace detail {

template <typename...> struct sc_void { typedef void type; };

template <typename T, typename = void>
struct sc_has_add : std::false_type {};
template <typename T>
struct sc_has_add<T, typename sc_void<decltype(
    std::declval<const T&>() + std::declval<const T&>())>::type>
    : std::is_convertible<decltype(std::declval<const T&>()
                                   + std::declval<const T&>()), T> {};

template <typename T, typename = void>
struct sc_has_sub : std::false_type {};
template <typename T>
struct sc_has_sub<T, typename sc_void<decltype(
    std::declval<const T&>() - std::declval<const T&>())>::type>
    : std::is_convertible<decltype(std::declval<const T&>()
                                   - std::declval<const T&>()), T> {};

template <typename T, typename = void>
struct sc_has_mul : std::false_type {};
template <typename T>
struct sc_has_mul<T, typename sc_void<decltype(
    std::declval<const T&>() * std::declval<const T&>())>::type>
    : std::is_convertible<decltype(std::declval<const T&>()
                                   * std::declval<const T&>()), T> {};

template <typename T, typename = void>
struct sc_has_div : std::false_type {};
template <typename T>
struct sc_has_div<T, typename sc_void<decltype(
    std::declval<const T&>() / std::declval<const T&>())>::type>
    : std::is_convertible<decltype(std::declval<const T&>()
                                   / std::declval<const T&>()), T> {};

template <typename T, typename = void>
struct sc_has_neg : std::false_type {};
template <typename T>
struct sc_has_neg<T, typename sc_void<decltype(
    -std::declval<const T&>())>::type>
    : std::is_convertible<decltype(-std::declval<const T&>()), T> {};

// convert_traits<T>::from_rational(bigint, bigint) is declared with a
// T-valued result (the enclose-once entry point; C-1). Unevaluated-context
// check: specializations replacing the primary template must keep the
// signature.
template <typename T, typename = void>
struct sc_has_convert : std::false_type {};
template <typename T>
struct sc_has_convert<T, typename sc_void<decltype(
    convert_traits<T>::from_rational(std::declval<const bigint&>(),
                                     std::declval<const bigint&>()))>::type>
    : std::is_convertible<decltype(
          convert_traits<T>::from_rational(std::declval<const bigint&>(),
                                           std::declval<const bigint&>())),
          T> {};

} // namespace detail

// ---------------------------------------------------------------------------
// bfem_scalar_traits<T>: instantiate (via require()) to check the C-1
// contract; each violated clause fails with its own named message.
// ---------------------------------------------------------------------------
template <typename T>
struct bfem_scalar_traits {
    static_assert(std::is_default_constructible<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires a default constructor T()");
    static_assert(std::is_constructible<T, int>::value,
                  "bfem_scalar_traits<T>: C-1 requires construction from int (T(int))");
    static_assert(detail::sc_has_add<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires the ring operator + (const T& + const T& -> T)");
    static_assert(detail::sc_has_sub<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires the ring operator - (const T& - const T& -> T)");
    static_assert(detail::sc_has_mul<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires the ring operator * (const T& * const T& -> T)");
    static_assert(detail::sc_has_neg<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires the unary operator - (-const T& -> T)");
    static_assert(detail::sc_has_div<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires operator / (division inside the declared budget, C-3-1)");
    static_assert(detail::sc_has_convert<T>::value,
                  "bfem_scalar_traits<T>: C-1 requires convert_traits<T>::from_rational(bigint, bigint) -> T (enclose-once)");

    // C-1 counting clause: the battery above is exactly the countable
    // operation set (counting_scalar implements precisely these operators),
    // so a passing T keeps operation counting possible by construction.
    static const bool counting_declared = true;

    // call at a space construction point to force the checks
    static void require() {}
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_DETAIL_SCALAR_TRAITS_HPP
