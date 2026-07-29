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

#include <atomic>
#include <cstddef>
#include <exception>

#include <vcp/tsparse/tsparse_spgemm.hpp>

namespace vcp {

// ---------------------------------------------------------------------------
// SPC-P1 shared shape (policy_add / policy_sub / policy_mul_impl):
// inputs are read directly when already finalized CSR (zero-copy); only a
// finalized CSC input is converted via as_csr() (a necessary cost).  The
// output is built by receiving emit directly into CSR arrays and installing
// them with assign_csr() -- no COO staging, no finalize(), no re-sort.
// This DEPENDS on the emit ordering contract of the tsparse_spgemm kernels:
// emit is called in ascending row i, ascending column j within each row,
// duplicate-merged and exact-zero-free (std::sort(touched) + workspace sweep,
// tsparse_spgemm.hpp L49-57 / L87-93).  Rows with no emitted entry are fixed
// up by the monotone fill pass after the kernel call.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// policy_add: C = A + B  (element-wise sparse addition)
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_add(
	const spmats<_T, _Index>& A,
	const spmats<_T, _Index>& B) const
{
	if (!A.is_finalized()) A.finalize();
	if (!B.is_finalized()) B.finalize();
	const spmats<_T, _Index>* Ap; spmats<_T, _Index> Ac_storage;
	if (A.format() == vcp::sparse_csr) { Ap = &A; }
	else { Ac_storage = A.as_csr(); Ap = &Ac_storage; }
	const spmats<_T, _Index>* Bp; spmats<_T, _Index> Bc_storage;
	if (B.format() == vcp::sparse_csr) { Bp = &B; }
	else { Bc_storage = B.as_csr(); Bp = &Bc_storage; }
	const std::size_t nrows = index_to_size(A.rowsize(), "spmats::policy_add");
	std::vector<_Index> new_outer(nrows + 1, _Index(0));
	std::vector<_Index> new_inner;
	std::vector<_T> new_value;
	const std::size_t cap = index_to_size(Ap->stored_nnz(), "spmats::policy_add")
	                      + index_to_size(Bp->stored_nnz(), "spmats::policy_add");
#if VCP_SPARSE_USE_OPENMP
	// SPOMP-1 A1 (design §3.1/§4.1): two-pass parallel builder above the
	// work threshold (work = nnz(A)+nnz(B) = cap); bit-identical to the
	// emit path below (identical per-row algorithm and order).  The emit
	// path below stays untouched and serves the below-threshold case.
	if (cap >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD)) {
		vcp::tsparse_spgemm::csr_csr_linear_combination_par(
			Ap->rowsize(), Ap->columnsize(),
			Ap->outer_index(), Ap->inner_index(), Ap->values(),
			Bp->outer_index(), Bp->inner_index(), Bp->values(),
			_T(1), _T(1),
			new_outer, new_inner, new_value);
		spmats<_T, _Index> C;
		C.clear();
		C.assign_csr(A.rowsize(), A.columnsize(), new_outer, new_inner, new_value);
		return C;
	}
#endif
	new_inner.reserve(cap);
	new_value.reserve(cap);
	vcp::tsparse_spgemm::csr_csr_linear_combination(
		Ap->rowsize(), Ap->columnsize(),
		Ap->outer_index(), Ap->inner_index(), Ap->values(),
		Bp->outer_index(), Bp->inner_index(), Bp->values(),
		_T(1), _T(1),
		[&](const _Index i, const _Index j, const _T& val) {
			new_inner.push_back(j);
			new_value.push_back(val);
			new_outer[static_cast<std::size_t>(i) + 1] =
				size_to_index(new_inner.size(), "spmats::policy_add");
		});
	for (std::size_t r = 1; r <= nrows; r++) {
		if (new_outer[r] < new_outer[r - 1]) new_outer[r] = new_outer[r - 1];
	}
	spmats<_T, _Index> C;
	C.clear();
	C.assign_csr(A.rowsize(), A.columnsize(), new_outer, new_inner, new_value);
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
	if (!A.is_finalized()) A.finalize();
	if (!B.is_finalized()) B.finalize();
	const spmats<_T, _Index>* Ap; spmats<_T, _Index> Ac_storage;
	if (A.format() == vcp::sparse_csr) { Ap = &A; }
	else { Ac_storage = A.as_csr(); Ap = &Ac_storage; }
	const spmats<_T, _Index>* Bp; spmats<_T, _Index> Bc_storage;
	if (B.format() == vcp::sparse_csr) { Bp = &B; }
	else { Bc_storage = B.as_csr(); Bp = &Bc_storage; }
	const std::size_t nrows = index_to_size(A.rowsize(), "spmats::policy_sub");
	std::vector<_Index> new_outer(nrows + 1, _Index(0));
	std::vector<_Index> new_inner;
	std::vector<_T> new_value;
	const std::size_t cap = index_to_size(Ap->stored_nnz(), "spmats::policy_sub")
	                      + index_to_size(Bp->stored_nnz(), "spmats::policy_sub");
#if VCP_SPARSE_USE_OPENMP
	// SPOMP-1 A1 (design §3.1/§4.1): same parallel-builder switch as
	// policy_add, with beta = -1.
	if (cap >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD)) {
		vcp::tsparse_spgemm::csr_csr_linear_combination_par(
			Ap->rowsize(), Ap->columnsize(),
			Ap->outer_index(), Ap->inner_index(), Ap->values(),
			Bp->outer_index(), Bp->inner_index(), Bp->values(),
			_T(1), _T(-1),
			new_outer, new_inner, new_value);
		spmats<_T, _Index> C;
		C.clear();
		C.assign_csr(A.rowsize(), A.columnsize(), new_outer, new_inner, new_value);
		return C;
	}
#endif
	new_inner.reserve(cap);
	new_value.reserve(cap);
	vcp::tsparse_spgemm::csr_csr_linear_combination(
		Ap->rowsize(), Ap->columnsize(),
		Ap->outer_index(), Ap->inner_index(), Ap->values(),
		Bp->outer_index(), Bp->inner_index(), Bp->values(),
		_T(1), _T(-1),
		[&](const _Index i, const _Index j, const _T& val) {
			new_inner.push_back(j);
			new_value.push_back(val);
			new_outer[static_cast<std::size_t>(i) + 1] =
				size_to_index(new_inner.size(), "spmats::policy_sub");
		});
	for (std::size_t r = 1; r <= nrows; r++) {
		if (new_outer[r] < new_outer[r - 1]) new_outer[r] = new_outer[r - 1];
	}
	spmats<_T, _Index> C;
	C.clear();
	C.assign_csr(A.rowsize(), A.columnsize(), new_outer, new_inner, new_value);
	return C;
}

// ---------------------------------------------------------------------------
// policy_mul: C = A * B  (sparse matrix-matrix multiply)
// NVI outer: non-virtual, finalizes A and B, then delegates to the virtual
// policy_mul_impl. Must never be overridden — override policy_mul_impl.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_mul(
	const spmats<_T, _Index>& A,
	const spmats<_T, _Index>& B) const
{
	if (!A.is_finalized()) A.finalize();
	if (!B.is_finalized()) B.finalize();
	return policy_mul_impl(A, B);
}

// ---------------------------------------------------------------------------
// policy_mul_impl: virtual algorithm body (spgemm). Custom policies override
// this, not policy_mul.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_mul_impl(
	const spmats<_T, _Index>& A,
	const spmats<_T, _Index>& B) const
{
	if (!A.is_finalized()) A.finalize();
	if (!B.is_finalized()) B.finalize();
	const spmats<_T, _Index>* Ap; spmats<_T, _Index> Ac_storage;
	if (A.format() == vcp::sparse_csr) { Ap = &A; }
	else { Ac_storage = A.as_csr(); Ap = &Ac_storage; }
	const spmats<_T, _Index>* Bp; spmats<_T, _Index> Bc_storage;
	if (B.format() == vcp::sparse_csr) { Bp = &B; }
	else { Bc_storage = B.as_csr(); Bp = &Bc_storage; }
	const std::size_t nrows = index_to_size(A.rowsize(), "spmats::policy_mul_impl");
	std::vector<_Index> new_outer(nrows + 1, _Index(0));
	std::vector<_Index> new_inner;
	std::vector<_T> new_value;
#if VCP_SPARSE_USE_OPENMP
	// SPOMP-1 A2 (design §3.1/§4.1): two-pass parallel builder above the
	// work threshold (work = nnz(A)+nnz(B)); bit-identical to the emit
	// path below (identical per-row algorithm and order).
	{
		const std::size_t par_work = index_to_size(Ap->stored_nnz(), "spmats::policy_mul_impl")
		                           + index_to_size(Bp->stored_nnz(), "spmats::policy_mul_impl");
		if (par_work >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD)) {
			vcp::tsparse_spgemm::csr_csr_multiply_par(
				Ap->rowsize(), Ap->columnsize(), Bp->columnsize(),
				Ap->outer_index(), Ap->inner_index(), Ap->values(),
				Bp->outer_index(), Bp->inner_index(), Bp->values(),
				new_outer, new_inner, new_value);
			spmats<_T, _Index> C;
			C.clear();
			C.assign_csr(A.rowsize(), B.columnsize(), new_outer, new_inner, new_value);
			return C;
		}
	}
#endif
	// spgemm output size has no cheap a-priori bound: no reserve, amortized
	// push_back (the previous COO path was push_back-based too).
	vcp::tsparse_spgemm::csr_csr_multiply(
		Ap->rowsize(), Ap->columnsize(), Bp->columnsize(),
		Ap->outer_index(), Ap->inner_index(), Ap->values(),
		Bp->outer_index(), Bp->inner_index(), Bp->values(),
		[&](const _Index i, const _Index j, const _T& val) {
			new_inner.push_back(j);
			new_value.push_back(val);
			new_outer[static_cast<std::size_t>(i) + 1] =
				size_to_index(new_inner.size(), "spmats::policy_mul_impl");
		});
	for (std::size_t r = 1; r <= nrows; r++) {
		if (new_outer[r] < new_outer[r - 1]) new_outer[r] = new_outer[r - 1];
	}
	spmats<_T, _Index> C;
	C.clear();
	C.assign_csr(A.rowsize(), B.columnsize(), new_outer, new_inner, new_value);
	return C;
}

// ---------------------------------------------------------------------------
// policy_mul_vec: y = A * x
// Auto-finalizes A at the policy layer entry. The low-level A.mul_vec()
// itself keeps its require_finalized() throwing behavior unchanged (see
// spmats_finalize_policy.md §5); only this policy-layer entry point gains
// auto-finalize.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
std::vector<_T> spmats<_T, _Index>::policy_mul_vec(
	const spmats<_T, _Index>& A,
	const std::vector<_T>& x) const
{
	if (!A.is_finalized()) A.finalize();
	return A.mul_vec(x);
}

// ---------------------------------------------------------------------------
// policy_left_mul_vec: y = x^T * A  (= A^T * x)
// Same auto-finalize policy as policy_mul_vec; A.trans_mul_vec() itself
// keeps its require_finalized() throwing behavior unchanged.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
std::vector<_T> spmats<_T, _Index>::policy_left_mul_vec(
	const std::vector<_T>& x,
	const spmats<_T, _Index>& A) const
{
	if (!A.is_finalized()) A.finalize();
	return A.trans_mul_vec(x);
}

// ---------------------------------------------------------------------------
// SPC-P1 destructive pattern-invariant policies: policy_mulsm / policy_divms
// / policy_minusm.  In-place, zero extra allocation: the stored value array
// of *this is updated directly.  Entries that become exact zero (alpha==0,
// underflow, zeroing division) are compacted away on the fly so the
// finalized-storage invariant (no stored zeros) is preserved.  The zero test
// is the established kernel idiom `!(v == _T(0))` (same as tsparse_spgemm /
// tcoo_remove_zeros / the previous functional policies).
// The compaction loop is format-neutral: outer/inner are relative names, so
// the identical code is correct for both finalized CSR and finalized CSC
// (the major-axis count differs only).  A not-yet-finalized (COO) matrix is
// finalized first.  Note the read cursor p0 is carried from the previous
// row's original end: outer[i+1] is overwritten with the compacted end, so
// it must not be re-read as the next row's start once compaction occurred.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
void spmats<_T, _Index>::policy_mulsm(const _T& alpha)
{
#if VCP_SPARSE_USE_OPENMP
	// SPOMP-1 B1 (design §4.2, pass-2 modified -- see
	// sandbox/docs/reports/SPOMP-1_stop_report.md): phase 1 multiplies every
	// stored value IN PLACE (one multiplication per element, exactly the
	// operation of the sequential path below) and counts survivors, in
	// parallel.  If nothing became exact zero (the common case) the pattern
	// is unchanged and outer/inner/value are already the final state.  Only
	// when zeros appeared does the rare, memory-bound front-packing run --
	// SEQUENTIALLY: the design's parallel pass 2 has a read/write race once
	// any entry is dropped (a later line's write range [cnt_i, cnt_{i+1})
	// can reach into an earlier line's not-yet-read source range), so the
	// sequential compaction of the A5 precedent is used instead.  No
	// O(rows) count array is needed in this form (the SPC-P1 zero-extra-
	// allocation invariant holds unmodified).  Bit-identical either way.
	if (!is_finalized()) finalize();
	const std::size_t nouter = (fmt == vcp::sparse_csr)
		? index_to_size(row, "spmats::policy_mulsm")
		: index_to_size(column, "spmats::policy_mulsm");
	std::atomic<bool> caught(false);
	std::exception_ptr eptr;
	const std::ptrdiff_t nnz_n = static_cast<std::ptrdiff_t>(value.size());
	std::ptrdiff_t survivors = 0;
#pragma omp parallel for schedule(static) reduction(+:survivors) if (static_cast<std::size_t>(nnz_n) >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD))
	for (std::ptrdiff_t k = 0; k < nnz_n; k++) {
		try {
			const std::size_t sp = static_cast<std::size_t>(k);
			value[sp] = alpha * value[sp];
			if (!(value[sp] == _T(0))) survivors++;
		}
		catch (...) {
			bool expected = false;
			if (caught.compare_exchange_strong(expected, true)) {
				eptr = std::current_exception();
			}
		}
	}
	if (caught.load()) std::rethrow_exception(eptr);
	if (survivors == nnz_n) return;
	std::size_t out = 0;
	index_type p0 = outer.empty() ? index_type(0) : outer[0];
	for (std::size_t i = 0; i < nouter; i++) {
		const index_type p1 = outer[i + 1];
		for (index_type p = p0; p < p1; p++) {
			const std::size_t sp = static_cast<std::size_t>(p);
			if (!(value[sp] == _T(0))) {
				value[out] = value[sp];
				inner[out] = inner[sp];
				out++;
			}
		}
		outer[i + 1] = size_to_index(out, "spmats::policy_mulsm");
		p0 = p1;
	}
	value.resize(out);
	inner.resize(out);
#else
	if (!is_finalized()) finalize();
	const std::size_t nouter = (fmt == vcp::sparse_csr)
		? index_to_size(row, "spmats::policy_mulsm")
		: index_to_size(column, "spmats::policy_mulsm");
	std::size_t out = 0;
	index_type p0 = outer.empty() ? index_type(0) : outer[0];
	for (std::size_t i = 0; i < nouter; i++) {
		const index_type p1 = outer[i + 1];
		for (index_type p = p0; p < p1; p++) {
			const std::size_t sp = static_cast<std::size_t>(p);
			const _T v = alpha * value[sp];
			if (!(v == _T(0))) {
				value[out] = v;
				inner[out] = inner[sp];
				out++;
			}
		}
		outer[i + 1] = size_to_index(out, "spmats::policy_mulsm");
		p0 = p1;
	}
	value.resize(out);
	inner.resize(out);
#endif
}

template <typename _T, typename _Index>
void spmats<_T, _Index>::policy_divms(const _T& alpha)
{
#if VCP_SPARSE_USE_OPENMP
	// SPOMP-1 B1: same two-phase scheme as policy_mulsm above (in-place
	// parallel divide + survivor count; sequential compaction only when
	// exact zeros appeared -- see the policy_mulsm comment / stop report).
	if (!is_finalized()) finalize();
	const std::size_t nouter = (fmt == vcp::sparse_csr)
		? index_to_size(row, "spmats::policy_divms")
		: index_to_size(column, "spmats::policy_divms");
	std::atomic<bool> caught(false);
	std::exception_ptr eptr;
	const std::ptrdiff_t nnz_n = static_cast<std::ptrdiff_t>(value.size());
	std::ptrdiff_t survivors = 0;
#pragma omp parallel for schedule(static) reduction(+:survivors) if (static_cast<std::size_t>(nnz_n) >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD))
	for (std::ptrdiff_t k = 0; k < nnz_n; k++) {
		try {
			const std::size_t sp = static_cast<std::size_t>(k);
			value[sp] = value[sp] / alpha;
			if (!(value[sp] == _T(0))) survivors++;
		}
		catch (...) {
			bool expected = false;
			if (caught.compare_exchange_strong(expected, true)) {
				eptr = std::current_exception();
			}
		}
	}
	if (caught.load()) std::rethrow_exception(eptr);
	if (survivors == nnz_n) return;
	std::size_t out = 0;
	index_type p0 = outer.empty() ? index_type(0) : outer[0];
	for (std::size_t i = 0; i < nouter; i++) {
		const index_type p1 = outer[i + 1];
		for (index_type p = p0; p < p1; p++) {
			const std::size_t sp = static_cast<std::size_t>(p);
			if (!(value[sp] == _T(0))) {
				value[out] = value[sp];
				inner[out] = inner[sp];
				out++;
			}
		}
		outer[i + 1] = size_to_index(out, "spmats::policy_divms");
		p0 = p1;
	}
	value.resize(out);
	inner.resize(out);
	return;
#else
	if (!is_finalized()) finalize();
	const std::size_t nouter = (fmt == vcp::sparse_csr)
		? index_to_size(row, "spmats::policy_divms")
		: index_to_size(column, "spmats::policy_divms");
	std::size_t out = 0;
	index_type p0 = outer.empty() ? index_type(0) : outer[0];
	for (std::size_t i = 0; i < nouter; i++) {
		const index_type p1 = outer[i + 1];
		for (index_type p = p0; p < p1; p++) {
			const std::size_t sp = static_cast<std::size_t>(p);
			const _T v = value[sp] / alpha;
			if (!(v == _T(0))) {
				value[out] = v;
				inner[out] = inner[sp];
				out++;
			}
		}
		outer[i + 1] = size_to_index(out, "spmats::policy_divms");
		p0 = p1;
	}
	value.resize(out);
	inner.resize(out);
#endif
}

template <typename _T, typename _Index>
void spmats<_T, _Index>::policy_minusm()
{
#if VCP_SPARSE_USE_OPENMP
	// SPOMP-1 A4 (design §3.1): negation cannot create a zero from a
	// stored non-zero (-x == 0 iff x == 0), so the pattern is invariant
	// and no compaction pass is needed -- plain elementwise parallel map
	// over the value array (work = nnz).  Bit-identical to the sequential
	// path below: one negation per element, stored at the same position.
	if (!is_finalized()) finalize();
	std::atomic<bool> caught(false);
	std::exception_ptr eptr;
	const std::ptrdiff_t nnz_n = static_cast<std::ptrdiff_t>(value.size());
#pragma omp parallel for schedule(static) if (static_cast<std::size_t>(nnz_n) >= static_cast<std::size_t>(VCP_SPMATS_OMP_THRESHOLD))
	for (std::ptrdiff_t k = 0; k < nnz_n; k++) {
		try {
			const std::size_t sp = static_cast<std::size_t>(k);
			value[sp] = -value[sp];
		}
		catch (...) {
			bool expected = false;
			if (caught.compare_exchange_strong(expected, true)) {
				eptr = std::current_exception();
			}
		}
	}
	if (caught.load()) std::rethrow_exception(eptr);
#else
	// Negation cannot create a zero from a stored non-zero (-x == 0 iff
	// x == 0), but the unified compaction shape is kept for code sharing
	// (the branch cost is negligible).
	if (!is_finalized()) finalize();
	const std::size_t nouter = (fmt == vcp::sparse_csr)
		? index_to_size(row, "spmats::policy_minusm")
		: index_to_size(column, "spmats::policy_minusm");
	std::size_t out = 0;
	index_type p0 = outer.empty() ? index_type(0) : outer[0];
	for (std::size_t i = 0; i < nouter; i++) {
		const index_type p1 = outer[i + 1];
		for (index_type p = p0; p < p1; p++) {
			const std::size_t sp = static_cast<std::size_t>(p);
			const _T v = -value[sp];
			if (!(v == _T(0))) {
				value[out] = v;
				inner[out] = inner[sp];
				out++;
			}
		}
		outer[i + 1] = size_to_index(out, "spmats::policy_minusm");
		p0 = p1;
	}
	value.resize(out);
	inner.resize(out);
#endif
}

// ---------------------------------------------------------------------------
// policy_scalar_mul: B = alpha * A
// SPC-P1: one as_csr() copy (the returned matrix itself) + delegation to the
// destructive policy_mulsm.  Still category 4 (self-contained): as_csr() is
// correct regardless of A's finalize state and A is never modified.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_scalar_mul(
	const spmats<_T, _Index>& A,
	const _T& alpha) const
{
	spmats<_T, _Index> C = A.as_csr();
	C.policy_mulsm(alpha);
	return C;
}

// ---------------------------------------------------------------------------
// policy_scalar_div: B = A / alpha  (divide each non-zero by alpha)
// SPC-P1: one as_csr() copy + delegation to the destructive policy_divms.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_scalar_div(
	const spmats<_T, _Index>& A,
	const _T& alpha) const
{
	spmats<_T, _Index> C = A.as_csr();
	C.policy_divms(alpha);
	return C;
}

// ---------------------------------------------------------------------------
// policy_neg: B = -A
// SPC-P1: one as_csr() copy + delegation to the destructive policy_minusm.
// ---------------------------------------------------------------------------
template <typename _T, typename _Index>
spmats<_T, _Index> spmats<_T, _Index>::policy_neg(
	const spmats<_T, _Index>& A) const
{
	spmats<_T, _Index> C = A.as_csr();
	C.policy_minusm();
	return C;
}

} // namespace vcp

#endif // VCP_SPMATS_PRODUCT_HPP
