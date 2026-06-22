// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_EXPERIMENTAL_HPP
#define VCP_TSPARSE_EXPERIMENTAL_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace vcp {
namespace tsparse_experimental {

struct block_lanczos_workspace {
	std::size_t block_size;
	std::vector<std::size_t> active_columns;
	block_lanczos_workspace() : block_size(0), active_columns() {}
};

struct block_arnoldi_workspace {
	std::size_t block_size;
	std::vector<std::size_t> active_columns;
	block_arnoldi_workspace() : block_size(0), active_columns() {}
};

struct thick_restart_workspace {
	std::size_t retained_dimension;
	std::vector<std::size_t> retained_indices;
	thick_restart_workspace() : retained_dimension(0), retained_indices() {}
};

struct locking_workspace {
	std::vector<std::size_t> locked_indices;
	std::string policy_name;
	locking_workspace() : locked_indices(), policy_name() {}
};

struct selective_reorthogonalization_workspace {
	std::vector<std::size_t> monitored_indices;
	std::vector<long double> estimated_loss;
	selective_reorthogonalization_workspace() : monitored_indices(), estimated_loss() {}
};

} // namespace tsparse_experimental
} // namespace vcp

#endif
