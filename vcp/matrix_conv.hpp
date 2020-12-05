#pragma once

#ifndef VCP_MATRIX_CONV_HPP
#define VCP_MATRIX_CONV_HPP

#include <vcp/vcp_converter.hpp>

namespace vcp {
	template<typename _T1, typename _T2, class _P1, class _P2> void convert(const vcp::matrix< _T1, _P1 >& A, vcp::matrix< _T2, _P2 >& B) {
		B.zeros(A.rowsize(), A.columnsize());
		for (int i = 0; i < A.rowsize(); i++) {
			for (int j = 0; j < A.columnsize(); j++) {
				convert(A(i, j), B(i, j));
			}
		}
	}
}
#endif