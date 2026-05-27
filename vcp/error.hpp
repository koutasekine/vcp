// VCP Library
// http ://verified.computation.jp

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
