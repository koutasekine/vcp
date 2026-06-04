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

#ifndef VCP_FIO_HPP
#define VCP_FIO_HPP

#include <string>
#include <fstream>
#include <limits>
#include <cstdio>

#include <vcp/error.hpp>

namespace vcp {
	namespace fio_detail {
		inline std::string filename_with_suffix(const char* name, const char* suffix) {
			if (name == nullptr) {
				vcp::throw_error<vcp::io_error>("matrix file I/O: null file name");
			}
			return std::string(name) + suffix;
		}

		inline std::string temporary_filename(const std::string& filename) {
			return filename + ".tmp";
		}

		inline void validate_matrix_size(const int row, const int column, const std::string& filename) {
			if (row < 0 || column < 0) {
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: invalid matrix size in ", filename,
					": ", row, " x ", column);
			}
			if ((row == 0) != (column == 0)) {
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: invalid zero matrix size in ", filename,
					": ", row, " x ", column);
			}
			if (row != 0 && column > (std::numeric_limits<int>::max)() / row) {
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: matrix size overflow in ", filename,
					": ", row, " x ", column);
			}
		}

		inline void write_bytes(
			std::ostream& os,
			const char* data,
			const std::streamsize size,
			const std::string& filename,
			const char* what
		) {
			os.write(data, size);
			if (!os) {
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: failed to write ", what, " to ", filename);
			}
		}

		inline void read_bytes(
			std::istream& is,
			char* data,
			const std::streamsize size,
			const std::string& filename,
			const char* what
		) {
			is.read(data, size);
			if (!is) {
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: failed to read ", what, " from ", filename);
			}
		}

		template <class T>
		void write_value(
			std::ostream& os,
			const T& value,
			const std::string& filename,
			const char* what
		) {
			write_bytes(os, reinterpret_cast<const char*>(&value), sizeof(T), filename, what);
		}

		template <class T>
		void read_value(
			std::istream& is,
			T& value,
			const std::string& filename,
			const char* what
		) {
			read_bytes(is, reinterpret_cast<char*>(&value), sizeof(T), filename, what);
		}

		inline void close_output(std::ofstream& os, const std::string& filename) {
			os.close();
			if (!os) {
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: failed to close output file: ", filename);
			}
		}

		inline void replace_output_file(const std::string& temporary, const std::string& filename) {
#ifdef _WIN32
			std::remove(filename.c_str());
#endif
			if (std::rename(temporary.c_str(), filename.c_str()) != 0) {
				std::remove(temporary.c_str());
				vcp::throw_error<vcp::io_error>(
					"matrix file I/O: failed to replace output file: ", filename);
			}
		}
	}

	template< class _P > void save(vcp::matrix< int, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_int");
		std::string temporary = fio_detail::temporary_filename(filename);

		std::ofstream savefile;
		savefile.open(temporary, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			vcp::throw_error<vcp::io_error>("save matrix<int>: cannot open file: ", temporary);
		}
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		fio_detail::validate_matrix_size(row, column, filename);
		fio_detail::write_bytes(savefile, "integer", 7, temporary, "type tag");
		fio_detail::write_value(savefile, row, temporary, "row size");
		fio_detail::write_value(savefile, column, temporary, "column size");

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				fio_detail::write_value(savefile, A(i, j), temporary, "matrix<int> element");
			}
		}
		fio_detail::close_output(savefile, temporary);
		fio_detail::replace_output_file(temporary, filename);
	}
	template< class _P > void load(vcp::matrix< int, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_int");

		std::ifstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			vcp::throw_error<vcp::io_error>("load matrix<int>: cannot open file: ", filename);
		}
		char d[7];
		fio_detail::read_bytes(loadfile, d, 7, filename, "type tag");
		std::string check_type(d, 7);
		if (check_type != "integer") {
			vcp::throw_error<vcp::io_error>(
				"load matrix<int>: data type mismatch: ", check_type);
		}
		int row, column;
		fio_detail::read_value(loadfile, row, filename, "row size");
		fio_detail::read_value(loadfile, column, filename, "column size");
		fio_detail::validate_matrix_size(row, column, filename);

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				fio_detail::read_value(loadfile, A(i, j), filename, "matrix<int> element");
			}
		}
		loadfile.close();
	}

	template< class _P > void save(vcp::matrix< double, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_d");
		std::string temporary = fio_detail::temporary_filename(filename);
		
		std::ofstream savefile;
		savefile.open(temporary, std::ios::out | std::ios::binary| std::ios::trunc);
		if (!savefile.is_open()) {
			vcp::throw_error<vcp::io_error>("save matrix<double>: cannot open file: ", temporary);
		}
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		fio_detail::validate_matrix_size(row, column, filename);
		fio_detail::write_bytes(savefile, "double", 6, temporary, "type tag");
		fio_detail::write_value(savefile, row, temporary, "row size");
		fio_detail::write_value(savefile, column, temporary, "column size");

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				fio_detail::write_value(savefile, A(i, j), temporary, "matrix<double> element");
			}
		}
		fio_detail::close_output(savefile, temporary);
		fio_detail::replace_output_file(temporary, filename);
	}
	template< class _P > void load(vcp::matrix< double, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_d");

		std::ifstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			vcp::throw_error<vcp::io_error>("load matrix<double>: cannot open file: ", filename);
		}
		char d[6];
		fio_detail::read_bytes(loadfile, d, 6, filename, "type tag");
		std::string check_type(d, 6);
		if (check_type != "double") {
			vcp::throw_error<vcp::io_error>(
				"load matrix<double>: data type mismatch: ", check_type);
		}
		int row, column;
		fio_detail::read_value(loadfile, row, filename, "row size");
		fio_detail::read_value(loadfile, column, filename, "column size");
		fio_detail::validate_matrix_size(row, column, filename);

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				fio_detail::read_value(loadfile, A(i, j), filename, "matrix<double> element");
			}
		}
		loadfile.close();
	}

#ifdef INTERVAL_HPP
	template< class _P > void save(vcp::matrix< kv::interval< double >, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_id");
		std::string temporary = fio_detail::temporary_filename(filename);

		std::ofstream savefile;
		savefile.open(temporary, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			vcp::throw_error<vcp::io_error>("save matrix<interval<double>>: cannot open file: ", temporary);
		}
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		fio_detail::validate_matrix_size(row, column, filename);
		fio_detail::write_bytes(savefile, "interval_double", 15, temporary, "type tag");
		fio_detail::write_value(savefile, row, temporary, "row size");
		fio_detail::write_value(savefile, column, temporary, "column size");

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				double lower = A(i, j).lower();
				double upper = A(i, j).upper();
				fio_detail::write_value(savefile, lower, temporary, "interval<double> lower bound");
				fio_detail::write_value(savefile, upper, temporary, "interval<double> upper bound");
			}
		}
		fio_detail::close_output(savefile, temporary);
		fio_detail::replace_output_file(temporary, filename);
	}
	template< class _P > void load(vcp::matrix< kv::interval< double >, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_id");

		std::ifstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			vcp::throw_error<vcp::io_error>("load matrix<interval<double>>: cannot open file: ", filename);
		}
		char d[15];
		fio_detail::read_bytes(loadfile, d, 15, filename, "type tag");
		std::string check_type(d, 15);
		if (check_type != "interval_double") {
			vcp::throw_error<vcp::io_error>(
				"load matrix<interval<double>>: data type mismatch: ", check_type);
		}
		int row, column;
		fio_detail::read_value(loadfile, row, filename, "row size");
		fio_detail::read_value(loadfile, column, filename, "column size");
		fio_detail::validate_matrix_size(row, column, filename);

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				double lower, upper;
				fio_detail::read_value(loadfile, lower, filename, "interval<double> lower bound");
				fio_detail::read_value(loadfile, upper, filename, "interval<double> upper bound");
				A(i, j).assign(lower, upper);
			}
		}
		loadfile.close();
	}
#endif

#ifdef DD_HPP
	template< class _P > void save(vcp::matrix< kv::dd, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_kvdd");
		std::string temporary = fio_detail::temporary_filename(filename);

		std::ofstream savefile;
		savefile.open(temporary, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			vcp::throw_error<vcp::io_error>("save matrix<kv::dd>: cannot open file: ", temporary);
		}
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		fio_detail::validate_matrix_size(row, column, filename);
		fio_detail::write_bytes(savefile, "kv_doubledouble", 15, temporary, "type tag");
		fio_detail::write_value(savefile, row, temporary, "row size");
		fio_detail::write_value(savefile, column, temporary, "column size");

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				fio_detail::write_value(savefile, A(i, j).a1, temporary, "kv::dd high word");
				fio_detail::write_value(savefile, A(i, j).a2, temporary, "kv::dd low word");
			}
		}
		fio_detail::close_output(savefile, temporary);
		fio_detail::replace_output_file(temporary, filename);
	}
	template< class _P > void load(vcp::matrix< kv::dd, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_kvdd");

		std::ifstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			vcp::throw_error<vcp::io_error>("load matrix<kv::dd>: cannot open file: ", filename);
		}
		char d[15];
		fio_detail::read_bytes(loadfile, d, 15, filename, "type tag");
		std::string check_type(d, 15);
		if (check_type != "kv_doubledouble") {
			vcp::throw_error<vcp::io_error>(
				"load matrix<kv::dd>: data type mismatch: ", check_type);
		}
		int row, column;
		fio_detail::read_value(loadfile, row, filename, "row size");
		fio_detail::read_value(loadfile, column, filename, "column size");
		fio_detail::validate_matrix_size(row, column, filename);

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				fio_detail::read_value(loadfile, A(i, j).a1, filename, "kv::dd high word");
				fio_detail::read_value(loadfile, A(i, j).a2, filename, "kv::dd low word");
			}
		}
		loadfile.close();
	}
#endif

#if defined(INTERVAL_HPP) && defined(DD_HPP)
	template< class _P > void save(vcp::matrix< kv::interval< kv::dd >, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_ikvdd");
		std::string temporary = fio_detail::temporary_filename(filename);

		std::ofstream savefile;
		savefile.open(temporary, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			vcp::throw_error<vcp::io_error>("save matrix<interval<kv::dd>>: cannot open file: ", temporary);
		}
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		fio_detail::validate_matrix_size(row, column, filename);
		fio_detail::write_bytes(savefile, "interval_kv_dd", 14, temporary, "type tag");
		fio_detail::write_value(savefile, row, temporary, "row size");
		fio_detail::write_value(savefile, column, temporary, "column size");

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				kv::dd lower = A(i, j).lower();
				kv::dd upper = A(i, j).upper();
				fio_detail::write_value(savefile, lower.a1, temporary, "interval<kv::dd> lower high word");
				fio_detail::write_value(savefile, lower.a2, temporary, "interval<kv::dd> lower low word");
				fio_detail::write_value(savefile, upper.a1, temporary, "interval<kv::dd> upper high word");
				fio_detail::write_value(savefile, upper.a2, temporary, "interval<kv::dd> upper low word");
			}
		}
		fio_detail::close_output(savefile, temporary);
		fio_detail::replace_output_file(temporary, filename);
	}
	template< class _P > void load(vcp::matrix< kv::interval< kv::dd >, _P >& A, const char* name) {
		std::string filename = fio_detail::filename_with_suffix(name, ".matrix_ikvdd");

		std::ifstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			vcp::throw_error<vcp::io_error>("load matrix<interval<kv::dd>>: cannot open file: ", filename);
		}
		char d[14];
		fio_detail::read_bytes(loadfile, d, 14, filename, "type tag");
		std::string check_type(d, 14);
		if (check_type != "interval_kv_dd") {
			vcp::throw_error<vcp::io_error>(
				"load matrix<interval<kv::dd>>: data type mismatch: ", check_type);
		}
		int row, column;
		fio_detail::read_value(loadfile, row, filename, "row size");
		fio_detail::read_value(loadfile, column, filename, "column size");
		fio_detail::validate_matrix_size(row, column, filename);

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				kv::dd lower, upper;
				fio_detail::read_value(loadfile, lower.a1, filename, "interval<kv::dd> lower high word");
				fio_detail::read_value(loadfile, lower.a2, filename, "interval<kv::dd> lower low word");
				fio_detail::read_value(loadfile, upper.a1, filename, "interval<kv::dd> upper high word");
				fio_detail::read_value(loadfile, upper.a2, filename, "interval<kv::dd> upper low word");
				A(i, j).assign(lower, upper);
			}
		}
		loadfile.close();
	}
#endif

	template< class _P > void save(vcp::matrix< int, _P >& A, const std::string& name) {
		vcp::save(A, name.c_str());
	}
	template< class _P > void load(vcp::matrix< int, _P >& A, const std::string& name) {
		vcp::load(A, name.c_str());
	}
	template< class _P > void save(vcp::matrix< double, _P >& A, const std::string& name) {
		vcp::save(A, name.c_str());
	}
	template< class _P > void load(vcp::matrix< double, _P >& A, const std::string& name) {
		vcp::load(A, name.c_str());
	}

#ifdef DD_HPP
	template< class _P > void save(vcp::matrix< kv::dd, _P >& A, const std::string& name) {
		vcp::save(A, name.c_str());
	}
	template< class _P > void load(vcp::matrix< kv::dd, _P >& A, const std::string& name) {
		vcp::load(A, name.c_str());
	}
#endif

#ifdef INTERVAL_HPP
	template< class _P > void save(vcp::matrix< kv::interval< double >, _P >& A, const std::string& name) {
		vcp::save(A, name.c_str());
	}
	template< class _P > void load(vcp::matrix< kv::interval< double >, _P >& A, const std::string& name) {
		vcp::load(A, name.c_str());
	}
#endif

#if defined(INTERVAL_HPP) && defined(DD_HPP)
	template< class _P > void save(vcp::matrix< kv::interval< kv::dd >, _P >& A, const std::string& name) {
		vcp::save(A, name.c_str());
	}
	template< class _P > void load(vcp::matrix< kv::interval< kv::dd >, _P >& A, const std::string& name) {
		vcp::load(A, name.c_str());
	}
#endif

}
#endif //VCP_FIO_HPP
