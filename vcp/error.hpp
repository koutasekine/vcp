// VCP Library
// http ://verified.computation.jp
//   
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met :
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and / or other materials provided with the distribution.
// * Neither the name of the Kouta Sekine nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED.IN NO EVENT SHALL KOUTA SEKINE BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#pragma once

#ifndef VCP_ERROR_HPP
#define VCP_ERROR_HPP

#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace vcp {

	class error : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	class invalid_argument : public error {
	public:
		using error::error;
	};

	class dimension_error : public error {
	public:
		using error::error;
	};

	class index_error : public error {
	public:
		using error::error;
	};

	class state_error : public error {
	public:
		using error::error;
	};

	class domain_error : public error {
	public:
		using error::error;
	};

	class numerical_error : public error {
	public:
		using error::error;
	};

	class verification_error : public numerical_error {
	public:
		using numerical_error::numerical_error;
	};

	class backend_error : public numerical_error {
	public:
		backend_error(const std::string& backend, const std::string& routine, int code, const std::string& message)
			: numerical_error(message), backend_(backend), routine_(routine), code_(code) {}

		const std::string& backend() const noexcept { return backend_; }
		const std::string& routine() const noexcept { return routine_; }
		int code() const noexcept { return code_; }

	private:
		std::string backend_;
		std::string routine_;
		int code_;
	};

	class blas_error : public backend_error {
	public:
		blas_error(const std::string& routine, int code, const std::string& message)
			: backend_error("blas", routine, code, message) {}
	};

	class lapack_error : public backend_error {
	public:
		lapack_error(const std::string& routine, int info, const std::string& message)
			: backend_error("lapack", routine, info, message) {}

		int info() const noexcept { return code(); }
	};

	class io_error : public error {
	public:
		using error::error;
	};

	namespace detail {
		inline void append_message(std::ostringstream&) {}

		template <class Head, class... Tail>
		void append_message(std::ostringstream& oss, Head&& head, Tail&&... tail) {
			oss << std::forward<Head>(head);
			append_message(oss, std::forward<Tail>(tail)...);
		}
	}

	template <class Exception, class... Args>
	[[noreturn]] void throw_error(Args&&... args) {
		std::ostringstream oss;
		detail::append_message(oss, std::forward<Args>(args)...);
		throw Exception(oss.str());
	}
}

#endif
