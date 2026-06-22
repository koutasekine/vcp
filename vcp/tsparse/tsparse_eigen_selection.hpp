// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_EIGEN_SELECTION_HPP
#define VCP_TSPARSE_EIGEN_SELECTION_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <vector>

#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_eigen_selection {

template <typename R>
R target_distance(const std::complex<R>& z, const eig_target target, const R& shift) {
	switch (target) {
	case eig_target::largest_magnitude:
	case eig_target::smallest_magnitude:
		return tsparse_scalar::abs_value(z);
	case eig_target::largest_algebraic:
	case eig_target::smallest_algebraic:
		return z.real();
	case eig_target::target_magnitude:
		return tsparse_scalar::abs_value(tsparse_scalar::abs_value(z) - tsparse_scalar::abs_value(shift));
	case eig_target::target_real:
		return tsparse_scalar::abs_value(z.real() - shift);
	default:
		return z.real();
	}
}

template <typename R>
std::vector<std::size_t> select_eigen_indices(const std::vector<std::complex<R> >& values,
                                              const std::size_t k,
                                              const eig_target target,
                                              const R& shift) {
	std::vector<std::size_t> order(values.size());
	for (std::size_t i = 0; i < order.size(); i++) order[i] = i;
	std::sort(order.begin(), order.end(), [&](const std::size_t a, const std::size_t b) {
		const R ka = target_distance(values[a], target, shift);
		const R kb = target_distance(values[b], target, shift);
		if (ka != kb) {
			return (target == eig_target::largest_magnitude || target == eig_target::largest_algebraic)
				? ka > kb : ka < kb;
		}
		if (values[a].real() != values[b].real()) return values[a].real() < values[b].real();
		if (values[a].imag() != values[b].imag()) return values[a].imag() < values[b].imag();
		return a < b;
	});
	if (order.size() > k) order.resize(k);
	return order;
}

template <typename T>
std::vector<std::size_t> select_eigen_indices_from_real(const std::vector<T>& values,
                                                        const std::size_t k,
                                                        const eig_target target,
                                                        const typename tsparse_scalar::real_type<T>::type& shift) {
	typedef typename tsparse_scalar::real_type<T>::type R;
	std::vector<std::complex<R> > complex_values;
	complex_values.reserve(values.size());
	for (std::size_t i = 0; i < values.size(); i++) {
		complex_values.push_back(std::complex<R>(tsparse_scalar::real_part(values[i]), R(0)));
	}
	return select_eigen_indices(complex_values, k, target, shift);
}

template <typename T>
void sort_eigenpairs_by_target(std::vector<T>& values,
                               std::vector<std::vector<T> >& vectors,
                               std::vector<typename tsparse_scalar::real_type<T>::type>& residuals_abs,
                               std::vector<typename tsparse_scalar::real_type<T>::type>& residuals_rel,
                               const eig_target target,
                               const typename tsparse_scalar::real_type<T>::type& shift) {
	const std::vector<std::size_t> order = select_eigen_indices_from_real(values, values.size(), target, shift);
	std::vector<T> sorted_values;
	std::vector<std::vector<T> > sorted_vectors;
	std::vector<typename tsparse_scalar::real_type<T>::type> sorted_abs;
	std::vector<typename tsparse_scalar::real_type<T>::type> sorted_rel;
	sorted_values.reserve(order.size());
	if (vectors.size() == values.size()) sorted_vectors.reserve(order.size());
	if (residuals_abs.size() == values.size()) sorted_abs.reserve(order.size());
	if (residuals_rel.size() == values.size()) sorted_rel.reserve(order.size());
	for (std::size_t i = 0; i < order.size(); i++) {
		const std::size_t j = order[i];
		sorted_values.push_back(values[j]);
		if (vectors.size() == values.size()) sorted_vectors.push_back(vectors[j]);
		if (residuals_abs.size() == values.size()) sorted_abs.push_back(residuals_abs[j]);
		if (residuals_rel.size() == values.size()) sorted_rel.push_back(residuals_rel[j]);
	}
	values.swap(sorted_values);
	if (vectors.size() == order.size()) vectors.swap(sorted_vectors);
	if (residuals_abs.size() == order.size()) residuals_abs.swap(sorted_abs);
	if (residuals_rel.size() == order.size()) residuals_rel.swap(sorted_rel);
}

template <typename T>
std::vector<std::size_t> select_real_eigenpairs(const std::vector<T>& values,
                                                const std::size_t k,
                                                const eig_target target,
                                                const typename tsparse_scalar::real_type<T>::type& shift) {
	return select_eigen_indices_from_real(values, k, target, shift);
}

template <typename R>
std::vector<std::size_t> select_complex_ritz_values(const std::vector<std::complex<R> >& values,
                                                    const std::size_t k,
                                                    const eig_target target,
                                                    const R& shift) {
	return select_eigen_indices(values, k, target, shift);
}

} // namespace tsparse_eigen_selection
} // namespace vcp

#endif
