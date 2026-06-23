// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_SPMATS_POLICY_TRAITS_HPP
#define VCP_SPMATS_POLICY_TRAITS_HPP

#include <cstddef>
#include <type_traits>
#include <vector>

namespace vcp {

	// -----------------------------------------------------------------------
	// spmats_policy_traits
	//
	// Documents the minimum requirements that a storage policy P must
	// satisfy in order to be used as a base class for spmatrix<T, P>.
	//
	// Required public members of P:
	//   typedef index_type;
	//   typedef value_type;
	//   typedef format_type;
	//
	//   index_type rowsize() const;
	//   index_type columnsize() const;
	//   index_type nnz() const;
	//   index_type stored_nnz() const;
	//   bool is_finalized() const;
	//   bool is_sorted() const;
	//   bool is_unique() const;
	//   format_type format() const;
	//
	//   void resize(index_type, index_type);
	//   void clear();
	//   void reserve(index_type);
	//   void add(index_type, index_type, const value_type&);
	//   void set(index_type, index_type, const value_type&);
	//   value_type get(index_type, index_type) const;
	//   void finalize();
	//   void sort_coo();
	//   void normalize_coo();
	//   void to_csr();
	//   void to_csc();
	//   P as_csr() const;
	//   P as_csc() const;
	//
	//   std::vector<value_type> mul_vec(const std::vector<value_type>&) const;
	//   void mul_vec(const value_type*, value_type*) const;
	//   std::vector<value_type> trans_mul_vec(const std::vector<value_type>&) const;
	//   void trans_mul_vec(const value_type*, value_type*) const;
	//
	//   P transpose() const;
	//
	//   const std::vector<index_type>& outer_index() const;
	//   const std::vector<index_type>& inner_index() const;
	//   const std::vector<value_type>& values() const;
	//   const std::vector<index_type>& coo_rows() const;
	//   const std::vector<index_type>& coo_columns() const;
	//   const std::vector<value_type>& coo_values() const;
	//
	//   void assign_csr(index_type, index_type,
	//                   const std::vector<index_type>&,
	//                   const std::vector<index_type>&,
	//                   const std::vector<value_type>&);
	//   void assign_csc(index_type, index_type,
	//                   const std::vector<index_type>&,
	//                   const std::vector<index_type>&,
	//                   const std::vector<value_type>&);
	//
	// Required policy computation methods (implemented in spmats_eigs.hpp,
	// spmats_product.hpp, spmats_lss.hpp):
	//
	//   eig_result<T> policy_eigs_with_info(const P&, std::size_t, const eig_options<T>&) const;
	//
	//   template <class Prec>
	//   eig_result<T> policy_eigs_with_info(const P&, std::size_t, const eig_options<T>&,
	//                                        const Prec&) const;
	//
	//   eig_result<T> policy_generalized_eigs_with_info(const P&, const P&,
	//                                                    std::size_t, const eig_options<T>&) const;
	//
	//   template <class Prec>
	//   eig_result<T> policy_generalized_eigs_with_info(const P&, const P&,
	//                                                    std::size_t, const eig_options<T>&,
	//                                                    const Prec&) const;
	//
	//   P policy_add(const P&, const P&) const;
	//   P policy_sub(const P&, const P&) const;
	//   P policy_mul(const P&, const P&) const;
	//   std::vector<T> policy_mul_vec(const P&, const std::vector<T>&) const;
	//   std::vector<T> policy_left_mul_vec(const std::vector<T>&, const P&) const;
	//   P policy_scalar_mul(const P&, const T&) const;
	//   P policy_neg(const P&) const;
	//
	//   linear_solve_result<T> policy_lss_with_info(const P&, const std::vector<T>&,
	//                                                const linear_solve_options<T>&) const;
	//   std::vector<T> policy_lss(const P&, const std::vector<T>&,
	//                              const linear_solve_options<T>&) const;
	//
	//   bool is_symmetric(const real_type& tol) const;
	//   bool is_symmetric() const;
	// -----------------------------------------------------------------------

	// C++11 static_assert check: verifies that P has the required typedefs.
	// Call this in a context where P is a complete type.
	template <class P>
	struct spmats_policy_check {
		typedef typename P::index_type index_type;
		typedef typename P::value_type value_type;
		typedef typename P::format_type format_type;

		// index_type must be a signed integer type
		static_assert(std::is_signed<index_type>::value,
		              "spmats policy: index_type must be a signed integer type");
	};

} // namespace vcp

#endif // VCP_SPMATS_POLICY_TRAITS_HPP
