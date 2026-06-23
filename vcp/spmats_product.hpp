// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spmats_product.hpp
// Policy method implementations for arithmetic operations on spmats<_T,_Index>.
// This file is included inside the namespace vcp {} block, AFTER the closing
// brace of spmats<_T,_Index>, via spmats.hpp.

#ifndef VCP_SPMATS_PRODUCT_HPP
#define VCP_SPMATS_PRODUCT_HPP

#include <vcp/tsparse/tsparse_spgemm.hpp>

namespace vcp {

// ---------------------------------------------------------------------------
// policy_add: C = A + B  (element-wise sparse addition)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_add(
	const spmats<_T, _Index>& A,
	const spmats<_T, _Index>& B) const
{
	spmats<_T, _Index> Ac = A.as_csr();
	spmats<_T, _Index> Bc = B.as_csr();
	spmats<_T, _Index> C;
	C.resize(A.rowsize(), A.columnsize());
	C.reserve(Ac.stored_nnz() + Bc.stored_nnz());
	vcp::tsparse_spgemm::csr_csr_linear_combination(
		Ac.rowsize(), Ac.columnsize(),
		Ac.outer_index(), Ac.inner_index(), Ac.values(),
		Bc.outer_index(), Bc.inner_index(), Bc.values(),
		_T(1), _T(1),
		[&](const _Index i, const _Index j, const _T& val) {
			C.add(i, j, val);
		});
	C.finalize();
	return C;
}

// ---------------------------------------------------------------------------
// policy_sub: C = A - B  (element-wise sparse subtraction)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_sub(
	const spmats<_T, _Index>& A,
	const spmats<_T, _Index>& B) const
{
	spmats<_T, _Index> Ac = A.as_csr();
	spmats<_T, _Index> Bc = B.as_csr();
	spmats<_T, _Index> C;
	C.resize(A.rowsize(), A.columnsize());
	C.reserve(Ac.stored_nnz() + Bc.stored_nnz());
	vcp::tsparse_spgemm::csr_csr_linear_combination(
		Ac.rowsize(), Ac.columnsize(),
		Ac.outer_index(), Ac.inner_index(), Ac.values(),
		Bc.outer_index(), Bc.inner_index(), Bc.values(),
		_T(1), _T(-1),
		[&](const _Index i, const _Index j, const _T& val) {
			C.add(i, j, val);
		});
	C.finalize();
	return C;
}

// ---------------------------------------------------------------------------
// policy_mul: C = A * B  (sparse matrix-matrix multiply)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_mul(
	const spmats<_T, _Index>& A,
	const spmats<_T, _Index>& B) const
{
	spmats<_T, _Index> Ac = A.as_csr();
	spmats<_T, _Index> Bc = B.as_csr();
	spmats<_T, _Index> C;
	C.resize(A.rowsize(), B.columnsize());
	vcp::tsparse_spgemm::csr_csr_multiply(
		Ac.rowsize(), Ac.columnsize(), Bc.columnsize(),
		Ac.outer_index(), Ac.inner_index(), Ac.values(),
		Bc.outer_index(), Bc.inner_index(), Bc.values(),
		[&](const _Index i, const _Index j, const _T& val) {
			C.add(i, j, val);
		});
	C.finalize();
	return C;
}

// ---------------------------------------------------------------------------
// policy_mul_vec: y = A * x
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
std::vector<_T> spmats<_T, _Index>::policy_mul_vec(
	const spmats<_T, _Index>& A,
	const std::vector<_T>& x) const
{
	return A.mul_vec(x);
}

// ---------------------------------------------------------------------------
// policy_left_mul_vec: y = x^T * A  (= A^T * x)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
std::vector<_T> spmats<_T, _Index>::policy_left_mul_vec(
	const std::vector<_T>& x,
	const spmats<_T, _Index>& A) const
{
	return A.trans_mul_vec(x);
}

// ---------------------------------------------------------------------------
// policy_scalar_mul: B = alpha * A
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_scalar_mul(
	const spmats<_T, _Index>& A,
	const _T& alpha) const
{
	spmats<_T, _Index> Ac = A.as_csr();
	spmats<_T, _Index> C;
	C.resize(Ac.rowsize(), Ac.columnsize());
	const std::vector<_Index>& outer = Ac.outer_index();
	const std::vector<_Index>& inner = Ac.inner_index();
	const std::vector<_T>& val = Ac.values();
	C.reserve(Ac.stored_nnz());
	for (_Index i = 0; i < Ac.rowsize(); i++) {
		for (_Index p = outer[static_cast<std::size_t>(i)];
		     p < outer[static_cast<std::size_t>(i + 1)]; p++) {
			const _T v = alpha * val[static_cast<std::size_t>(p)];
			if (!(v == _T(0))) {
				C.add(i, inner[static_cast<std::size_t>(p)], v);
			}
		}
	}
	C.finalize();
	return C;
}

// ---------------------------------------------------------------------------
// policy_neg: B = -A
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_neg(
	const spmats<_T, _Index>& A) const
{
	return policy_scalar_mul(A, _T(-1));
}

} // namespace vcp

#endif // VCP_SPMATS_PRODUCT_HPP
