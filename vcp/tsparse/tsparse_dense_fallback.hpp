// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_DENSE_FALLBACK_HPP
#define VCP_TSPARSE_DENSE_FALLBACK_HPP

#include <cstddef>
#include <limits>

#include <vcp/error.hpp>

namespace vcp {
	namespace tsparse_dense_fallback {
		inline void check_dense_size(const std::size_t rows, const std::size_t cols,
		                             const std::size_t max_entries, const bool allow) {
			if (!allow) {
				vcp::throw_error<vcp::state_error>("dense fallback is disabled");
			}
			if (rows != 0 && cols > (std::numeric_limits<std::size_t>::max)() / rows) {
				vcp::throw_error<vcp::invalid_argument>("dense fallback size overflow");
			}
			if (rows * cols > max_entries) {
				vcp::throw_error<vcp::state_error>("dense fallback size exceeds max_dense_size");
			}
		}
	}
}

#endif
