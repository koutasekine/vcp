// VCP Library
// http ://verified.computation.jp

#pragma once

#ifndef VCP_HWROUND_GUARD_HPP
#define VCP_HWROUND_GUARD_HPP

#include <kv/hwround.hpp>

namespace vcp {

	class hwround_guard {
	public:
		hwround_guard() = default;
		~hwround_guard() noexcept {
			kv::hwround::roundnear();
		}

		hwround_guard(const hwround_guard&) = delete;
		hwround_guard& operator=(const hwround_guard&) = delete;
	};

}

#endif
