#pragma once

#ifndef VCP_IMATS_ASSIST_HPP
#define VCP_IMATS_ASSIST_HPP

namespace vcp {
#ifdef INTERVAL_HPP
	template<typename _T, class _P1, class _P2> void mid(const vcp::matrix< kv::interval< _T >, _P1 >& A, vcp::matrix< _T, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				B(i, j) = mid(A(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> void rad(const vcp::matrix< kv::interval< _T >, _P1 >& A, vcp::matrix< _T, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				B(i, j) = rad(A(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> void midrad(const vcp::matrix< kv::interval< _T >, _P1 >& A, vcp::matrix< _T, _P2 >& B, vcp::matrix< _T, _P2 >& C) {
		B.zeros(A.rowsize(), A.columnsize());
		C.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				midrad(A(i, j), B(i, j), C(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> void interval(const vcp::matrix< _T, _P1 >& A, vcp::matrix< kv::interval< _T >, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				B(i, j) = kv::interval< _T >(A(i, j));
			}
		}
	}

	template<typename _T, class _P1, class _P2> vcp::matrix< kv::interval< _T >, _P1 > operator*(const vcp::matrix< kv::interval< _T >, _P1 >& A, const  vcp::matrix< _T, _P2 >& B) {
		if (A.columnsize() != B.rowsize()) {
			std::cout << "vmatmul:error " << A.columnsize() << " != " << B.rowsize() << std::endl;
			exit(1);
		}

		vcp::matrix< kv::interval< _T >, _P1 > C;
		C.zeros(A.rowsize(), B.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				for (int k = 0; k < B.columnsize(); k++) {
					C(i, k) = A(i, j) * B(j, k);
				}
			}
		}

		return std::move(C);
	}
	template<typename _T, class _P1, class _P2> vcp::matrix< kv::interval< _T >, _P2 > operator*(const vcp::matrix< _T, _P1 >& A, const vcp::matrix< kv::interval< _T >, _P2 >& B) {
		if (A.columnsize() != B.rowsize()) {
			std::cout << "vmatmul:error " << A.columnsize() << " != " << B.rowsize() << std::endl;
			exit(1);
		}
		vcp::matrix< kv::interval< _T >, _P2 > C;
		C.zeros(A.rowsize(), B.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				for (int k = 0; k < B.columnsize(); k++) {
					C(i, k) = A(i, j) * B(j, k);
				}
			}
		}

		return std::move(C);
	}
#endif
}
#endif // VCP_IMATS_ASSIST_HPP