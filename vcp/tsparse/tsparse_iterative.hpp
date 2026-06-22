// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_ITERATIVE_HPP
#define VCP_TSPARSE_ITERATIVE_HPP

#include <vector>

namespace vcp {
	namespace tsparse_iterative {
		template <typename T>
		double relative_residual_denominator(const std::vector<T>& b_norm_holder, const double bnorm) {
			(void)b_norm_holder;
			return bnorm > 1.0 ? bnorm : 1.0;
		}
	}
}

#endif
