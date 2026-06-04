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

#ifndef VCP_PORTABLE_FIO_HPP
#define VCP_PORTABLE_FIO_HPP

#include <cstdio>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <limits>
#include <string>

#include <vcp/error.hpp>
#include <vcp/matrix.hpp>

namespace kv {
	class dd;
	template <class T> class interval;
}

namespace vcp {
	namespace portable_fio_detail {
		static const unsigned char magic[8] = {
			'V', 'C', 'P', 'M', 'A', 'T', '2', '\0'
		};
		static const std::uint32_t version = 2;
		static const std::uint32_t layout_column_major = 1;

		enum type_id {
			type_int32 = 1,
			type_double = 2,
			type_interval_double = 3,
			type_kv_dd = 4,
			type_interval_kv_dd = 5
		};

		inline const char* type_name(const std::uint32_t type) {
			switch (type) {
			case type_int32:
				return "int32";
			case type_double:
				return "double";
			case type_interval_double:
				return "interval<double>";
			case type_kv_dd:
				return "kv::dd";
			case type_interval_kv_dd:
				return "interval<kv::dd>";
			default:
				return "unknown";
			}
		}

		inline std::string filename_with_suffix(const char* name) {
			if (name == nullptr) {
				vcp::throw_error<vcp::io_error>("portable matrix I/O: null file name");
			}
			return std::string(name) + ".vcpmat";
		}

		inline std::string temporary_filename(const std::string& filename) {
			return filename + ".tmp";
		}

		inline void crc64_update(
			std::uint64_t& crc,
			const unsigned char* data,
			const std::size_t size
		) {
			static const std::uint64_t polynomial = 0x42f0e1eba9ea3693ULL;
			for (std::size_t i = 0; i < size; i++) {
				crc ^= static_cast<std::uint64_t>(data[i]) << 56;
				for (int j = 0; j < 8; j++) {
					if ((crc & 0x8000000000000000ULL) != 0) {
						crc = (crc << 1) ^ polynomial;
					}
					else {
						crc <<= 1;
					}
				}
			}
		}

		inline void write_bytes(
			std::ostream& os,
			const unsigned char* data,
			const std::streamsize size,
			const std::string& filename,
			const char* what,
			std::uint64_t* crc = nullptr
		) {
			os.write(reinterpret_cast<const char*>(data), size);
			if (!os) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: failed to write ", what, " to ", filename);
			}
			if (crc != nullptr) {
				crc64_update(*crc, data, static_cast<std::size_t>(size));
			}
		}

		inline void read_bytes(
			std::istream& is,
			unsigned char* data,
			const std::streamsize size,
			const std::string& filename,
			const char* what,
			std::uint64_t* crc = nullptr
		) {
			is.read(reinterpret_cast<char*>(data), size);
			if (!is) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: failed to read ", what, " from ", filename);
			}
			if (crc != nullptr) {
				crc64_update(*crc, data, static_cast<std::size_t>(size));
			}
		}

		inline void write_le32(
			std::ostream& os,
			const std::uint32_t value,
			const std::string& filename,
			const char* what,
			std::uint64_t* crc = nullptr
		) {
			unsigned char bytes[4];
			for (int i = 0; i < 4; i++) {
				bytes[i] = static_cast<unsigned char>((value >> (8 * i)) & 0xffU);
			}
			write_bytes(os, bytes, 4, filename, what, crc);
		}

		inline std::uint32_t read_le32(
			std::istream& is,
			const std::string& filename,
			const char* what,
			std::uint64_t* crc = nullptr
		) {
			unsigned char bytes[4];
			read_bytes(is, bytes, 4, filename, what, crc);
			std::uint32_t value = 0;
			for (int i = 0; i < 4; i++) {
				value |= static_cast<std::uint32_t>(bytes[i]) << (8 * i);
			}
			return value;
		}

		inline void write_le64(
			std::ostream& os,
			const std::uint64_t value,
			const std::string& filename,
			const char* what,
			std::uint64_t* crc = nullptr
		) {
			unsigned char bytes[8];
			for (int i = 0; i < 8; i++) {
				bytes[i] = static_cast<unsigned char>((value >> (8 * i)) & 0xffU);
			}
			write_bytes(os, bytes, 8, filename, what, crc);
		}

		inline std::uint64_t read_le64(
			std::istream& is,
			const std::string& filename,
			const char* what,
			std::uint64_t* crc = nullptr
		) {
			unsigned char bytes[8];
			read_bytes(is, bytes, 8, filename, what, crc);
			std::uint64_t value = 0;
			for (int i = 0; i < 8; i++) {
				value |= static_cast<std::uint64_t>(bytes[i]) << (8 * i);
			}
			return value;
		}

		inline void ensure_binary64(const std::string& filename) {
			if (sizeof(double) != 8 || !std::numeric_limits<double>::is_iec559) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: IEEE754 binary64 double is required: ", filename);
			}
		}

		inline std::uint64_t double_to_bits(const double value) {
			std::uint64_t bits;
			std::memcpy(&bits, &value, sizeof(bits));
			return bits;
		}

		inline double bits_to_double(const std::uint64_t bits) {
			double value;
			std::memcpy(&value, &bits, sizeof(value));
			return value;
		}

		inline void write_double(
			std::ostream& os,
			const double value,
			const std::string& filename,
			const char* what,
			std::uint64_t& crc
		) {
			write_le64(os, double_to_bits(value), filename, what, &crc);
		}

		inline double read_double(
			std::istream& is,
			const std::string& filename,
			const char* what,
			std::uint64_t& crc
		) {
			return bits_to_double(read_le64(is, filename, what, &crc));
		}

		inline void write_int32(
			std::ostream& os,
			const int value,
			const std::string& filename,
			const char* what,
			std::uint64_t& crc
		) {
			if (value < (std::numeric_limits<std::int32_t>::min)()
				|| value > (std::numeric_limits<std::int32_t>::max)()) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: int value is outside int32 range in ", filename);
			}
			std::uint32_t bits;
			if (value < 0) {
				const std::uint32_t magnitude = static_cast<std::uint32_t>(
					-static_cast<std::int64_t>(value));
				bits = 0U - magnitude;
			}
			else {
				bits = static_cast<std::uint32_t>(value);
			}
			write_le32(os, bits, filename, what, &crc);
		}

		inline int read_int32(
			std::istream& is,
			const std::string& filename,
			const char* what,
			std::uint64_t& crc
		) {
			std::uint32_t bits = read_le32(is, filename, what, &crc);
			std::int64_t value;
			if (bits <= 0x7fffffffU) {
				value = static_cast<std::int64_t>(bits);
			}
			else {
				value = -static_cast<std::int64_t>((~bits) + 1U);
			}
			if (value < static_cast<std::int64_t>((std::numeric_limits<int>::min)())
				|| value > static_cast<std::int64_t>((std::numeric_limits<int>::max)())) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: int32 value cannot be represented as int in ", filename);
			}
			return static_cast<int>(value);
		}

		inline std::uint64_t checked_element_count(
			const std::uint64_t rows,
			const std::uint64_t columns,
			const std::string& filename
		) {
			if ((rows == 0) != (columns == 0)) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: invalid zero matrix size in ", filename,
					": ", rows, " x ", columns);
			}
			if (rows != 0 && columns > (std::numeric_limits<std::uint64_t>::max)() / rows) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: matrix size overflow in ", filename);
			}
			return rows * columns;
		}

		inline std::uint64_t checked_payload_bytes(
			const std::uint64_t elements,
			const std::uint64_t bytes_per_element,
			const std::string& filename
		) {
			if (bytes_per_element != 0
				&& elements > (std::numeric_limits<std::uint64_t>::max)() / bytes_per_element) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: payload size overflow in ", filename);
			}
			return elements * bytes_per_element;
		}

		inline std::uint64_t bytes_per_element(const std::uint32_t type) {
			switch (type) {
			case type_int32:
				return 4;
			case type_double:
				return 8;
			case type_interval_double:
				return 16;
			case type_kv_dd:
				return 16;
			case type_interval_kv_dd:
				return 32;
			default:
				return 0;
			}
		}

		inline std::uint64_t expected_payload_bytes(
			const std::uint32_t type,
			const std::uint64_t rows,
			const std::uint64_t columns,
			const std::string& filename
		) {
			const std::uint64_t unit = bytes_per_element(type);
			if (unit == 0) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: unsupported type id in ", filename, ": ", type);
			}
			return checked_payload_bytes(checked_element_count(rows, columns, filename), unit, filename);
		}

		inline void validate_save_size(
			const int rows,
			const int columns,
			const std::string& filename
		) {
			if (rows < 0 || columns < 0) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: invalid matrix size in ", filename,
					": ", rows, " x ", columns);
			}
			checked_element_count(
				static_cast<std::uint64_t>(rows),
				static_cast<std::uint64_t>(columns),
				filename);
		}

		inline void validate_load_size(
			const std::uint64_t rows,
			const std::uint64_t columns,
			const std::string& filename
		) {
			const std::uint64_t elements = checked_element_count(rows, columns, filename);
			if (rows > static_cast<std::uint64_t>((std::numeric_limits<int>::max)())
				|| columns > static_cast<std::uint64_t>((std::numeric_limits<int>::max)())
				|| elements > static_cast<std::uint64_t>((std::numeric_limits<int>::max)())) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: matrix is too large for vcp::matrix in ", filename,
					": ", rows, " x ", columns);
			}
		}

		inline void write_header(
			std::ostream& os,
			const std::uint32_t type,
			const std::uint64_t rows,
			const std::uint64_t columns,
			const std::uint64_t payload_bytes,
			const std::string& filename,
			std::uint64_t& crc
		) {
			write_bytes(os, magic, 8, filename, "magic", &crc);
			write_le32(os, version, filename, "version", &crc);
			write_le32(os, type, filename, "type id", &crc);
			write_le32(os, layout_column_major, filename, "layout", &crc);
			write_le32(os, 0, filename, "reserved field", &crc);
			write_le64(os, rows, filename, "row size", &crc);
			write_le64(os, columns, filename, "column size", &crc);
			write_le64(os, payload_bytes, filename, "payload byte size", &crc);
		}

		struct header {
			std::uint32_t type;
			std::uint64_t rows;
			std::uint64_t columns;
			std::uint64_t payload_bytes;
		};

		inline header read_header(
			std::istream& is,
			const std::string& filename,
			std::uint64_t& crc
		) {
			unsigned char read_magic[8];
			read_bytes(is, read_magic, 8, filename, "magic", &crc);
			if (std::memcmp(read_magic, magic, 8) != 0) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: invalid magic in ", filename);
			}
			const std::uint32_t file_version = read_le32(is, filename, "version", &crc);
			if (file_version != version) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: unsupported version in ", filename,
					": ", file_version);
			}
			header h;
			h.type = read_le32(is, filename, "type id", &crc);
			const std::uint32_t layout = read_le32(is, filename, "layout", &crc);
			if (layout != layout_column_major) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: unsupported layout in ", filename,
					": ", layout);
			}
			const std::uint32_t reserved = read_le32(is, filename, "reserved field", &crc);
			if (reserved != 0) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: invalid reserved field in ", filename);
			}
			h.rows = read_le64(is, filename, "row size", &crc);
			h.columns = read_le64(is, filename, "column size", &crc);
			h.payload_bytes = read_le64(is, filename, "payload byte size", &crc);
			validate_load_size(h.rows, h.columns, filename);
			return h;
		}

		inline void check_type(
			const std::uint32_t file_type,
			const std::uint32_t expected_type,
			const std::string& filename
		) {
			if (file_type != expected_type) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: data type mismatch in ", filename,
					": file type is ", type_name(file_type),
					", requested type is ", type_name(expected_type));
			}
		}

		inline void check_payload_size(
			const header& h,
			const std::string& filename
		) {
			const std::uint64_t expected = expected_payload_bytes(
				h.type, h.rows, h.columns, filename);
			if (h.payload_bytes != expected) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: payload byte size mismatch in ", filename,
					": ", h.payload_bytes, " != ", expected);
			}
		}

		inline void close_output(std::ofstream& os, const std::string& filename) {
			os.close();
			if (!os) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: failed to close output file: ", filename);
			}
		}

		inline void replace_output_file(const std::string& temporary, const std::string& filename) {
#ifdef _WIN32
			std::remove(filename.c_str());
#endif
			if (std::rename(temporary.c_str(), filename.c_str()) != 0) {
				std::remove(temporary.c_str());
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: failed to replace output file: ", filename);
			}
		}

		template <class Writer>
		void save_matrix(
			const int rows,
			const int columns,
			const std::uint32_t type,
			const char* name,
			Writer writer
		) {
			std::string filename = filename_with_suffix(name);
			std::string temporary = temporary_filename(filename);
			validate_save_size(rows, columns, filename);
			if (type != type_int32) {
				ensure_binary64(filename);
			}

			const std::uint64_t urows = static_cast<std::uint64_t>(rows);
			const std::uint64_t ucolumns = static_cast<std::uint64_t>(columns);
			const std::uint64_t payload_bytes = expected_payload_bytes(
				type, urows, ucolumns, filename);

			std::ofstream savefile;
			savefile.open(temporary.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
			if (!savefile.is_open()) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: cannot open output file: ", temporary);
			}

			std::uint64_t crc = 0;
			write_header(savefile, type, urows, ucolumns, payload_bytes, temporary, crc);
			writer(savefile, temporary, crc);
			write_le64(savefile, crc, temporary, "CRC64");
			close_output(savefile, temporary);
			replace_output_file(temporary, filename);
		}

		template <class Reader>
		void load_matrix(
			const std::uint32_t expected_type,
			const char* name,
			Reader reader
		) {
			std::string filename = filename_with_suffix(name);
			if (expected_type != type_int32) {
				ensure_binary64(filename);
			}

			std::ifstream loadfile;
			loadfile.open(filename.c_str(), std::ios::in | std::ios::binary);
			if (!loadfile.is_open()) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: cannot open input file: ", filename);
			}

			std::uint64_t crc = 0;
			header h = read_header(loadfile, filename, crc);
			check_type(h.type, expected_type, filename);
			check_payload_size(h, filename);
			reader(loadfile, filename, h, crc);
			const std::uint64_t stored_crc = read_le64(loadfile, filename, "CRC64");
			if (stored_crc != crc) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: CRC64 mismatch in ", filename);
			}
			if (loadfile.peek() != std::char_traits<char>::eof()) {
				vcp::throw_error<vcp::io_error>(
					"portable matrix I/O: trailing data after CRC64 in ", filename);
			}
		}
	}

	template< class _P > void save_portable(const vcp::matrix< int, _P >& A, const char* name) {
		portable_fio_detail::save_matrix(
			A.rowsize(),
			A.columnsize(),
			portable_fio_detail::type_int32,
			name,
			[&A](std::ostream& os, const std::string& filename, std::uint64_t& crc) {
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						portable_fio_detail::write_int32(
							os, A(i, j), filename, "matrix<int> element", crc);
					}
				}
			});
	}

	template< class _P > void load_portable(vcp::matrix< int, _P >& A, const char* name) {
		portable_fio_detail::load_matrix(
			portable_fio_detail::type_int32,
			name,
			[&A](std::istream& is, const std::string& filename,
				const portable_fio_detail::header& h, std::uint64_t& crc) {
				A.zeros(static_cast<int>(h.rows), static_cast<int>(h.columns));
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						A(i, j) = portable_fio_detail::read_int32(
							is, filename, "matrix<int> element", crc);
					}
				}
			});
	}

	template< class _P > void save_portable(const vcp::matrix< double, _P >& A, const char* name) {
		portable_fio_detail::save_matrix(
			A.rowsize(),
			A.columnsize(),
			portable_fio_detail::type_double,
			name,
			[&A](std::ostream& os, const std::string& filename, std::uint64_t& crc) {
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						portable_fio_detail::write_double(
							os, A(i, j), filename, "matrix<double> element", crc);
					}
				}
			});
	}

	template< class _P > void load_portable(vcp::matrix< double, _P >& A, const char* name) {
		portable_fio_detail::load_matrix(
			portable_fio_detail::type_double,
			name,
			[&A](std::istream& is, const std::string& filename,
				const portable_fio_detail::header& h, std::uint64_t& crc) {
				A.zeros(static_cast<int>(h.rows), static_cast<int>(h.columns));
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						A(i, j) = portable_fio_detail::read_double(
							is, filename, "matrix<double> element", crc);
					}
				}
			});
	}

	template< class _P > void save_portable(
		const vcp::matrix< kv::interval< double >, _P >& A,
		const char* name
	) {
		portable_fio_detail::save_matrix(
			A.rowsize(),
			A.columnsize(),
			portable_fio_detail::type_interval_double,
			name,
			[&A](std::ostream& os, const std::string& filename, std::uint64_t& crc) {
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						portable_fio_detail::write_double(
							os, A(i, j).lower(), filename, "interval<double> lower bound", crc);
						portable_fio_detail::write_double(
							os, A(i, j).upper(), filename, "interval<double> upper bound", crc);
					}
				}
			});
	}

	template< class _P > void load_portable(
		vcp::matrix< kv::interval< double >, _P >& A,
		const char* name
	) {
		portable_fio_detail::load_matrix(
			portable_fio_detail::type_interval_double,
			name,
			[&A](std::istream& is, const std::string& filename,
				const portable_fio_detail::header& h, std::uint64_t& crc) {
				A.zeros(static_cast<int>(h.rows), static_cast<int>(h.columns));
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						const double lower = portable_fio_detail::read_double(
							is, filename, "interval<double> lower bound", crc);
						const double upper = portable_fio_detail::read_double(
							is, filename, "interval<double> upper bound", crc);
						A(i, j).assign(lower, upper);
					}
				}
			});
	}

	template< class _P > void save_portable(const vcp::matrix< kv::dd, _P >& A, const char* name) {
		portable_fio_detail::save_matrix(
			A.rowsize(),
			A.columnsize(),
			portable_fio_detail::type_kv_dd,
			name,
			[&A](std::ostream& os, const std::string& filename, std::uint64_t& crc) {
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						portable_fio_detail::write_double(
							os, A(i, j).a1, filename, "kv::dd high word", crc);
						portable_fio_detail::write_double(
							os, A(i, j).a2, filename, "kv::dd low word", crc);
					}
				}
			});
	}

	template< class _P > void load_portable(vcp::matrix< kv::dd, _P >& A, const char* name) {
		portable_fio_detail::load_matrix(
			portable_fio_detail::type_kv_dd,
			name,
			[&A](std::istream& is, const std::string& filename,
				const portable_fio_detail::header& h, std::uint64_t& crc) {
				A.zeros(static_cast<int>(h.rows), static_cast<int>(h.columns));
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						A(i, j).a1 = portable_fio_detail::read_double(
							is, filename, "kv::dd high word", crc);
						A(i, j).a2 = portable_fio_detail::read_double(
							is, filename, "kv::dd low word", crc);
					}
				}
			});
	}

	template< class _P > void save_portable(
		const vcp::matrix< kv::interval< kv::dd >, _P >& A,
		const char* name
	) {
		portable_fio_detail::save_matrix(
			A.rowsize(),
			A.columnsize(),
			portable_fio_detail::type_interval_kv_dd,
			name,
			[&A](std::ostream& os, const std::string& filename, std::uint64_t& crc) {
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						portable_fio_detail::write_double(
							os, A(i, j).lower().a1, filename,
							"interval<kv::dd> lower high word", crc);
						portable_fio_detail::write_double(
							os, A(i, j).lower().a2, filename,
							"interval<kv::dd> lower low word", crc);
						portable_fio_detail::write_double(
							os, A(i, j).upper().a1, filename,
							"interval<kv::dd> upper high word", crc);
						portable_fio_detail::write_double(
							os, A(i, j).upper().a2, filename,
							"interval<kv::dd> upper low word", crc);
					}
				}
			});
	}

	template< class _P > void load_portable(
		vcp::matrix< kv::interval< kv::dd >, _P >& A,
		const char* name
	) {
		portable_fio_detail::load_matrix(
			portable_fio_detail::type_interval_kv_dd,
			name,
			[&A](std::istream& is, const std::string& filename,
				const portable_fio_detail::header& h, std::uint64_t& crc) {
				A.zeros(static_cast<int>(h.rows), static_cast<int>(h.columns));
				for (int j = 0; j < A.columnsize(); j++) {
					for (int i = 0; i < A.rowsize(); i++) {
						A(i, j).lower().a1 = portable_fio_detail::read_double(
							is, filename, "interval<kv::dd> lower high word", crc);
						A(i, j).lower().a2 = portable_fio_detail::read_double(
							is, filename, "interval<kv::dd> lower low word", crc);
						A(i, j).upper().a1 = portable_fio_detail::read_double(
							is, filename, "interval<kv::dd> upper high word", crc);
						A(i, j).upper().a2 = portable_fio_detail::read_double(
							is, filename, "interval<kv::dd> upper low word", crc);
					}
				}
			});
	}

	template< class _P >
	void save_portable(const vcp::matrix< int, _P >& A, const std::string& name) {
		vcp::save_portable(A, name.c_str());
	}

	template< class _P >
	void load_portable(vcp::matrix< int, _P >& A, const std::string& name) {
		vcp::load_portable(A, name.c_str());
	}

	template< class _P >
	void save_portable(const vcp::matrix< double, _P >& A, const std::string& name) {
		vcp::save_portable(A, name.c_str());
	}

	template< class _P >
	void load_portable(vcp::matrix< double, _P >& A, const std::string& name) {
		vcp::load_portable(A, name.c_str());
	}

	template< class _P >
	void save_portable(
		const vcp::matrix< kv::interval< double >, _P >& A,
		const std::string& name
	) {
		vcp::save_portable(A, name.c_str());
	}

	template< class _P >
	void load_portable(
		vcp::matrix< kv::interval< double >, _P >& A,
		const std::string& name
	) {
		vcp::load_portable(A, name.c_str());
	}

	template< class _P >
	void save_portable(const vcp::matrix< kv::dd, _P >& A, const std::string& name) {
		vcp::save_portable(A, name.c_str());
	}

	template< class _P >
	void load_portable(vcp::matrix< kv::dd, _P >& A, const std::string& name) {
		vcp::load_portable(A, name.c_str());
	}

	template< class _P >
	void save_portable(
		const vcp::matrix< kv::interval< kv::dd >, _P >& A,
		const std::string& name
	) {
		vcp::save_portable(A, name.c_str());
	}

	template< class _P >
	void load_portable(
		vcp::matrix< kv::interval< kv::dd >, _P >& A,
		const std::string& name
	) {
		vcp::load_portable(A, name.c_str());
	}
}

#endif // VCP_PORTABLE_FIO_HPP
