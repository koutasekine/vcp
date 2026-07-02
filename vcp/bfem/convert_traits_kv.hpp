// vcp/bfem/convert_traits_kv.hpp
// Layer 0: opt-in helper header for kv types (external design v0.3, 6.3).
//
// Layer 0 itself never includes kv; users who want kv-specific behaviour of
// convert_traits include this header explicitly.
//
// NOTE (recorded design decision, L0 gate report):
//   The primary template convert_traits<T> already fulfils the mandatory
//   enclosure contract of external design 5.1 for kv::interval<F>:
//   horner_limbs performs integer construction steps that are exact in any
//   binary interval format wide enough for 16 bit chunks, and the single
//   final division is outward rounded by kv. For numerators/denominators
//   within the mantissa of F the result is already the tightest 1 ulp
//   enclosure, so no tighter specialization exists to write today.
//   This header is the designated location for future specializations
//   (e.g. correctly rounded kv::mpfr conversion); it currently only pulls in
//   the kv interval types so that a user TU has everything needed to
//   instantiate typed_registry<D, kv::interval<F>>.

#ifndef VCP_BFEM_CONVERT_TRAITS_KV_HPP
#define VCP_BFEM_CONVERT_TRAITS_KV_HPP

#include <kv/interval.hpp>
#include <kv/rdouble.hpp>

#include <vcp/bfem/convert_traits.hpp>

#endif // VCP_BFEM_CONVERT_TRAITS_KV_HPP
