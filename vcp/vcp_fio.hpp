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

namespace vcp {
	template< class _P > void save(vcp::matrix< int, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_int");
		filename = name + filename;

		std::fstream savefile;
		savefile.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			std::cout << "ERROR : save : matrix int" << std::endl;
			exit(1);
		}
		savefile.write("integer", 7);
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		savefile.write((char*)&row, sizeof(int));
		savefile.write((char*)&column, sizeof(int));

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				savefile.write((char*)&A(i, j), sizeof(int));
			}
		}
		savefile.close();
	}
	template< class _P > void load(vcp::matrix< int, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_int");
		filename = name + filename;

		std::fstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			std::cout << "ERROR : load : matrix int : Cannot Open the file : " << name << std::endl;
			exit(1);
		}
		char d[7];
		loadfile.read(d, 7);
		std::string check_type(d);
		if (check_type != "integer") {
			std::cout << "ERROR : load : No match the data type: check_type : " << check_type << std::endl;
			exit(1);
		}
		int row, column;
		loadfile.read((char*)& row, sizeof(int));
		loadfile.read((char*)& column, sizeof(int));

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				loadfile.read((char*)&A(i, j), sizeof(int));
			}
		}
		loadfile.close();
	}

	template< class _P > void save(vcp::matrix< double, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_d");
		filename = name + filename;
		
		std::fstream savefile;
		savefile.open(filename, std::ios::out | std::ios::binary| std::ios::trunc);
		if (!savefile.is_open()) {
			std::cout << "ERROR : save : matrix double" << std::endl;
			exit(1);
		}
		savefile.write("double", 6);
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		savefile.write((char*) &row, sizeof(int));
		savefile.write((char*) &column, sizeof(int));

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				savefile.write((char*)&A(i, j), sizeof(double));
			}
		}
		savefile.close();
	}
	template< class _P > void load(vcp::matrix< double, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_d");
		filename = name + filename;

		std::fstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			std::cout << "ERROR : load : matrix double : Cannot Open the file : " << name << std::endl;
			exit(1);
		}
		char d[6];
		loadfile.read(d, 6);
		std::string check_type(d);
		if (check_type != "double") {
			std::cout << "ERROR : load : No match the data type: check_type : " << check_type << std::endl;
			exit(1);
		}
		int row, column;
		loadfile.read((char*)& row, sizeof(int));
		loadfile.read((char*)& column, sizeof(int));

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				loadfile.read((char*)&A(i, j), sizeof(double));
			}
		}
		loadfile.close();
	}

#ifdef INTERVAL_HPP
	template< class _P > void save(vcp::matrix< kv::interval< double >, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_id");
		filename = name + filename;

		std::fstream savefile;
		savefile.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			std::cout << "ERROR : save : matrix interval double" << std::endl;
			exit(1);
		}
		savefile.write("interval_double", 15);
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		savefile.write((char*)&row, sizeof(int));
		savefile.write((char*)&column, sizeof(int));

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				savefile.write((char*)&A(i, j).lower(), sizeof(double));
				savefile.write((char*)&A(i, j).upper(), sizeof(double));
			}
		}
		savefile.close();
	}
	template< class _P > void load(vcp::matrix< kv::interval< double >, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_id");
		filename = name + filename;

		std::fstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			std::cout << "ERROR : load : matrix interval double : Cannot Open the file : " << name << std::endl;
			exit(1);
		}
		char d[15];
		loadfile.read(d, 15);
		std::string check_type(d);
		if (check_type.find("interval_double") == std::string::npos ) {
			std::cout << "ERROR : load : No match the data type: check_type : " << check_type << std::endl;
			exit(1);
		}
		int row, column;
		loadfile.read((char*)& row, sizeof(int));
		loadfile.read((char*)& column, sizeof(int));

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				loadfile.read((char*)&A(i, j).lower(), sizeof(double));
				loadfile.read((char*)&A(i, j).upper(), sizeof(double));
			}
		}
		loadfile.close();
	}
#endif

#ifdef DD_HPP
	template< class _P > void save(vcp::matrix< kv::dd, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_kvdd");
		filename = name + filename;

		std::fstream savefile;
		savefile.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			std::cout << "ERROR : save : matrix kv::dd" << std::endl;
			exit(1);
		}
		savefile.write("kv_doubledouble", 15);
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		savefile.write((char*)&row, sizeof(int));
		savefile.write((char*)&column, sizeof(int));

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				savefile.write((char*)&A(i, j).a1, sizeof(double));
				savefile.write((char*)&A(i, j).a2, sizeof(double));
			}
		}
		savefile.close();
	}
	template< class _P > void load(vcp::matrix< kv::dd, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_kvdd");
		filename = name + filename;

		std::fstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			std::cout << "ERROR : load : matrix kv::dd : Cannot Open the file : " << name << std::endl;
			exit(1);
		}
		char d[15];
		loadfile.read(d, 15);
		std::string check_type(d);
		if (check_type != "kv_doubledouble") {
			std::cout << "ERROR : load : No match the data type: check_type : " << check_type << std::endl;
			exit(1);
		}
		int row, column;
		loadfile.read((char*)& row, sizeof(int));
		loadfile.read((char*)& column, sizeof(int));

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				loadfile.read((char*)&A(i, j).a1, sizeof(double));
				loadfile.read((char*)&A(i, j).a2, sizeof(double));
			}
		}
		loadfile.close();
	}
#endif

#if defined(INTERVAL_HPP) && defined(DD_HPP)
	template< class _P > void save(vcp::matrix< kv::interval< kv::dd >, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_ikvdd");
		filename = name + filename;

		std::fstream savefile;
		savefile.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!savefile.is_open()) {
			std::cout << "ERROR : save : matrix interval kv::dd" << std::endl;
			exit(1);
		}
		savefile.write("interval_kv_dd", 14);
		int row, column;
		row = A.rowsize();
		column = A.columnsize();
		savefile.write((char*)&row, sizeof(int));
		savefile.write((char*)&column, sizeof(int));

		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				savefile.write((char*)&A(i, j).lower().a1, sizeof(double));
				savefile.write((char*)&A(i, j).lower().a2, sizeof(double));
				savefile.write((char*)&A(i, j).upper().a1, sizeof(double));
				savefile.write((char*)&A(i, j).upper().a2, sizeof(double));
			}
		}
		savefile.close();
	}
	template< class _P > void load(vcp::matrix< kv::interval< kv::dd >, _P >& A, const char* name) {
		std::string str(name);
		std::string filename(".matrix_ikvdd");
		filename = name + filename;

		std::fstream loadfile;
		loadfile.open(filename, std::ios::in | std::ios::binary);
		if (!loadfile.is_open()) {
			std::cout << "ERROR : load : matrix interval kv::dd : Cannot Open the file : " << name << std::endl;
			exit(1);
		}
		char d[14];
		loadfile.read(d, 14);
		std::string check_type(d);
		if (check_type != "interval_kv_dd") {
			std::cout << "ERROR : load : No match the data type: check_type : " << check_type << std::endl;
			exit(1);
		}
		int row, column;
		loadfile.read((char*)& row, sizeof(int));
		loadfile.read((char*)& column, sizeof(int));

		A.zeros(row, column);
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < column; j++) {
				loadfile.read((char*)&A(i, j).lower().a1, sizeof(double));
				loadfile.read((char*)&A(i, j).lower().a2, sizeof(double));
				loadfile.read((char*)&A(i, j).upper().a1, sizeof(double));
				loadfile.read((char*)&A(i, j).upper().a2, sizeof(double));
			}
		}
		loadfile.close();
	}
#endif

#ifdef INTERVAL_HPP
	template< class _T, class _P > void save(vcp::matrix< kv::interval< _T >, _P >& A, const std::string name) {
		const char* name2 = name.c_str();
		vcp::save(A, name2);
	}
#endif

}
#endif //VCP_FIO_HPP