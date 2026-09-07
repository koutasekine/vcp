// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
// Copyright(c) 2017, Kouta Sekine <k.sekine@computation.jp>
// All rights reserved.

#pragma once

#ifndef VCP_MATS2_HPP
#define VCP_MATS2_HPP

#include <algorithm>
#include <type_traits>
#include <vector>

#include <vcp/mats.hpp>
#include <vcp/tblas/tblas.hpp>
#include <vcp/tlapack/tlapack.hpp>

namespace vcp {

	namespace mats2_detail {
		inline void check_lapack_info(const char* name, const int info) {
			if (info < 0) {
				vcp::throw_error<vcp::domain_error>(name, ": invalid argument: ", -info);
			}
			if (info > 0) {
				vcp::throw_error<vcp::domain_error>(name, ": failed with info=", info);
			}
		}
	}

	template <typename _T> class mats2 : public mats< _T > {
		static_assert(!std::is_same<_T, bool>::value, "mats2<bool> is not supported. Use mbool for comparison results.");

		typedef mats< _T > base_type;

		void copy_upper_to_lower() {
			for (int j = 0; j < this->column; j++) {
				for (int i = j + 1; i < this->row; i++) {
					this->v[i + this->row * j] = this->v[j + this->row * i];
				}
			}
		}

		void zero_lower_triangle() {
			for (int j = 0; j < this->column; j++) {
				for (int i = j + 1; i < this->row; i++) {
					this->v[i + this->row * j] = _T(0);
				}
			}
		}

		void set_diagonal_from_vector(const std::vector< _T >& w) {
			const int nn = static_cast<int>(w.size());
			this->zeros(nn);
			for (int i = 0; i < nn; i++) {
				this->v[i + this->row * i] = w[i];
			}
		}

		// MATS-N64-B (ii): scale n contiguous elements (incx == 1) by alpha.
		// n <= INT_MAX: exactly the former single tscal call (bit-identical path).
		// n >  INT_MAX: 2^30-element slices; scal is element-wise (no reduction), so
		// every element is still computed as the same single product.
		static void scal_all(const vcp::index_t n, const _T& alpha, _T* x) {
			if (n <= 2147483647) {
				vcp::tscal(static_cast<int>(n), alpha, x, 1);
				return;
			}
			const vcp::index_t chunk = 1073741824;
			for (vcp::index_t off = 0; off < n; off += chunk) {
				const vcp::index_t len = (n - off < chunk) ? (n - off) : chunk;
				vcp::tscal(static_cast<int>(len), alpha, x + off, 1);
			}
		}

	public:
		using base_type::base_type;

		mats2() : base_type() {}
		~mats2() = default;
		mats2(const mats2&) = default;
		mats2(mats2&&) = default;
		mats2& operator=(const mats2&) = default;
		mats2& operator=(mats2&&) = default;

		// C = A*B
		virtual void mulmm(const mats2< _T >& B, mats2< _T >& c) const {
			if (this->type == 'S' && (B.type == 'C' || B.type == 'R' || B.type == 'M')) {
				c = B;
				scal_all(B.n, this->v[0], c.v.data());
				return;
			}
			if ((this->type == 'C' || this->type == 'R' || this->type == 'M') && B.type == 'S') {
				c = *this;
				scal_all(this->n, B.v[0], c.v.data());
				return;
			}
			if (this->type == 'S' && B.type == 'S') {
				c.zeros(1, 1);
				c.type = 'S';
				c.v[0] = this->v[0] * B.v[0];
				return;
			}
			if (this->column != B.row) {
				vcp::throw_error<vcp::dimension_error>(
					"mats2::mulmm: dimension mismatch: (", this->row, ", ", this->column,
					") * (", B.row, ", ", B.column, ")");
			}

			if (this->type == 'R' && B.type == 'C') {
				c.zeros(1, 1);
				c.type = 'S';
				c.v[0] = vcp::tdot(static_cast<int>(this->column), this->v.data(), 1, B.v.data(), 1);
				return;
			}
			if (this->type == 'C' && B.type == 'R') {
				c.zeros(static_cast<int>(this->row), static_cast<int>(B.column));
				vcp::tgemm('N', 'N', static_cast<int>(this->row), static_cast<int>(B.column), 1, _T(1), this->v.data(), static_cast<int>(this->row),
					B.v.data(), static_cast<int>(B.row), _T(0), c.v.data(), static_cast<int>(c.row));
				return;
			}
			if (this->type == 'M' && B.type == 'C') {
				c.zeros(static_cast<int>(this->row), 1);
				vcp::tgemm('N', 'N', static_cast<int>(this->row), 1, static_cast<int>(this->column), _T(1), this->v.data(), static_cast<int>(this->row),
					B.v.data(), static_cast<int>(B.row), _T(0), c.v.data(), static_cast<int>(c.row));
				return;
			}
			if (this->type == 'R' && B.type == 'M') {
				c.zeros(1, static_cast<int>(B.column));
				vcp::tgemm('N', 'N', 1, static_cast<int>(B.column), static_cast<int>(this->column), _T(1), this->v.data(), 1,
					B.v.data(), static_cast<int>(B.row), _T(0), c.v.data(), 1);
				return;
			}
			if (this->type == 'M' && B.type == 'M') {
				c.zeros(static_cast<int>(this->row), static_cast<int>(B.column));
				vcp::tgemm('N', 'N', static_cast<int>(this->row), static_cast<int>(B.column), static_cast<int>(this->column), _T(1), this->v.data(), static_cast<int>(this->row),
					B.v.data(), static_cast<int>(B.row), _T(0), c.v.data(), static_cast<int>(c.row));
				return;
			}
			vcp::throw_error<vcp::dimension_error>(
				"mats2::mulmm: unsupported matrix product: type ", this->type, " * ", B.type);
		}

		// C = transpose(A)*A
		virtual void mulltmm(mats2< _T >& c) const {
			if (this->type == 'S' || this->type == 'C' || this->type == 'R') {
				base_type tmp = *this;
				base_type out;
				tmp.mulltmm(out);
				c.row = out.row;
				c.column = out.column;
				c.n = out.n;
				c.type = out.type;
				c.v = out.v;
				return;
			}
			if (this->type != 'M') {
				vcp::throw_error<vcp::state_error>("mats2::mulltmm: unsupported matrix type: ", this->type);
			}
			c.zeros(static_cast<int>(this->column), static_cast<int>(this->column));
			vcp::tsyrk('U', 'T', static_cast<int>(this->column), static_cast<int>(this->row), _T(1), this->v.data(), static_cast<int>(this->row),
				_T(0), c.v.data(), static_cast<int>(c.row));
			c.copy_upper_to_lower();
		}

		void linearsolve(const mats2< _T >& b, mats2< _T >& x) {
			if (this->row != this->column || b.row != this->row) {
				vcp::throw_error<vcp::dimension_error>(
					"mats2::linearsolve: invalid dimensions: A=(", this->row, ", ", this->column,
					"), b=(", b.row, ", ", b.column, ")");
			}
			x = b;
			std::vector<int> ipiv(static_cast<std::size_t>(this->row), 0);
			const int info = vcp::tgesv(static_cast<int>(this->row), static_cast<int>(b.column), this->v.data(), static_cast<int>(this->row),
				ipiv.data(), x.v.data(), static_cast<int>(x.row));
			mats2_detail::check_lapack_info("mats2::linearsolve/tgesv", info);
		}

		virtual void inv() {
			if (this->row != this->column) {
				vcp::throw_error<vcp::dimension_error>(
					"mats2::inv: matrix must be square: ", this->row, " != ", this->column);
			}
			std::vector<int> ipiv(static_cast<std::size_t>(this->row), 0);
			int info = vcp::tgetrf(static_cast<int>(this->row), static_cast<int>(this->column), this->v.data(), static_cast<int>(this->row), ipiv.data());
			mats2_detail::check_lapack_info("mats2::inv/tgetrf", info);
			info = vcp::tgetri(static_cast<int>(this->row), this->v.data(), static_cast<int>(this->row), ipiv.data());
			mats2_detail::check_lapack_info("mats2::inv/tgetri", info);
		}

		virtual void Cholesky() {
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::Cholesky: matrix must be symmetric");
			}
			const int info = vcp::tpotrf('U', static_cast<int>(this->row), this->v.data(), static_cast<int>(this->row));
			mats2_detail::check_lapack_info("mats2::Cholesky/tpotrf", info);
			this->zero_lower_triangle();
		}

		virtual void eigsym(int itep = 1) {
			(void)itep;
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::eigsym: matrix must be symmetric");
			}
			std::vector< _T > w(static_cast<std::size_t>(this->row), _T(0));
			const int info = vcp::tsyev('N', 'U', static_cast<int>(this->row), this->v.data(), static_cast<int>(this->row), w.data());
			mats2_detail::check_lapack_info("mats2::eigsym/tsyev", info);
			this->set_diagonal_from_vector(w);
		}

		void eigsym(mats2< _T >& V, int itep = 1) {
			(void)itep;
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::eigsym: matrix must be symmetric");
			}
			std::vector< _T > w(static_cast<std::size_t>(this->row), _T(0));
			const int info = vcp::tsyev('V', 'U', static_cast<int>(this->row), this->v.data(), static_cast<int>(this->row), w.data());
			mats2_detail::check_lapack_info("mats2::eigsym/tsyev", info);
			V = *this;
			this->set_diagonal_from_vector(w);
		}

		void eigsymge(mats2< _T >& B, int itep = 1) {
			(void)itep;
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::eigsymge: matrix A must be symmetric");
			}
			if (!B.is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::eigsymge: matrix B must be symmetric");
			}
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>(
					"mats2::eigsymge: dimension mismatch: A=(", this->row, ", ", this->column,
					"), B=(", B.row, ", ", B.column, ")");
			}
			std::vector< _T > w(static_cast<std::size_t>(this->row), _T(0));
			const int info = vcp::tsygv(1, 'N', 'U', static_cast<int>(this->row), this->v.data(), static_cast<int>(this->row),
				B.v.data(), static_cast<int>(B.row), w.data());
			mats2_detail::check_lapack_info("mats2::eigsymge/tsygv", info);
			this->set_diagonal_from_vector(w);
		}

		void eigsymge(mats2< _T >& B, mats2< _T >& V, int itep = 1) {
			(void)itep;
			if (!this->is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::eigsymge: matrix A must be symmetric");
			}
			if (!B.is_symmetric()) {
				vcp::throw_error<vcp::domain_error>("mats2::eigsymge: matrix B must be symmetric");
			}
			if (this->row != B.row || this->column != B.column) {
				vcp::throw_error<vcp::dimension_error>(
					"mats2::eigsymge: dimension mismatch: A=(", this->row, ", ", this->column,
					"), B=(", B.row, ", ", B.column, ")");
			}
			std::vector< _T > w(static_cast<std::size_t>(this->row), _T(0));
			const int info = vcp::tsygv(1, 'V', 'U', static_cast<int>(this->row), this->v.data(), static_cast<int>(this->row),
				B.v.data(), static_cast<int>(B.row), w.data());
			mats2_detail::check_lapack_info("mats2::eigsymge/tsygv", info);
			V = *this;
			this->set_diagonal_from_vector(w);
		}
	};
}

#endif // VCP_MATS2_HPP
