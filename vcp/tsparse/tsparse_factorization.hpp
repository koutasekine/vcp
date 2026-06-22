// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_FACTORIZATION_HPP
#define VCP_TSPARSE_FACTORIZATION_HPP

#include <algorithm>
#include <cstddef>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
namespace tsparse_factorization {

// ---------------------------------------------------------------------------
// ILU(0) factorization stored in CSR format
// L has unit diagonal (implicit), U has explicit diagonal
// Combined storage: lu_val[k] is L below diagonal, U on/above diagonal
// ---------------------------------------------------------------------------
template <typename T, typename Index = int>
struct ilu0_data {
	typedef T value_type;
	typedef Index index_type;
	typedef typename tsparse_scalar::real_type<T>::type real_type;

	std::size_t n;
	std::vector<Index> row_ptr;    // size n+1
	std::vector<Index> col_ind;    // size nnz
	std::vector<T> lu_val;         // combined L+U values (L lower, U upper)
	std::vector<Index> diag_pos;   // position of diagonal in each row (size n)
	bool factorized;
	bool singular_or_unstable;
	std::size_t zero_pivots;
	std::size_t first_zero_pivot;
	real_type pivot_tol;
	std::string diagnostics;

	ilu0_data() : n(0), factorized(false), singular_or_unstable(false),
		zero_pivots(0), first_zero_pivot(0),
		pivot_tol(tsparse_scalar::decimal_power_negative<real_type>(14)),
		diagnostics() {}
};

template <typename T, typename Index>
void ilu0_record_zero_pivot(ilu0_data<T, Index>& result, const std::size_t pivot_index) {
	if (result.zero_pivots == 0) result.first_zero_pivot = pivot_index;
	result.zero_pivots++;
	result.singular_or_unstable = true;
}

template <typename T, typename Index>
std::string ilu0_make_diagnostics(const ilu0_data<T, Index>& data) {
	std::ostringstream os;
	if (data.zero_pivots == 0) {
		os << "factorized; zero_pivot_count=0; pivot_tolerance=" << data.pivot_tol;
	}
	else {
		os << "singular_or_unstable; zero_pivot_count=" << data.zero_pivots
		   << "; first_zero_pivot_index=" << data.first_zero_pivot
		   << "; pivot_tolerance=" << data.pivot_tol;
	}
	return os.str();
}

// ---------------------------------------------------------------------------
// Build ILU(0) from CSR data in-place (modifies lu_val)
// Input: sorted CSR matrix (row_ptr, col_ind, values)
// Output: fills result.lu_val with L (below diag) and U (on/above diag) entries
// ---------------------------------------------------------------------------
template <typename T, typename Index>
ilu0_data<T, Index> ilu0_factorize(
	const std::vector<Index>& row_ptr,
	const std::vector<Index>& col_ind,
	const std::vector<T>& values,
	const std::size_t n,
	const typename tsparse_scalar::real_type<T>::type& pivot_tol =
		tsparse_scalar::decimal_power_negative<typename tsparse_scalar::real_type<T>::type>(14))
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	ilu0_data<T, Index> result;
	result.n = n;
	result.row_ptr = row_ptr;
	result.col_ind = col_ind;
	result.lu_val = values;
	result.pivot_tol = pivot_tol;
	result.zero_pivots = 0;
	result.singular_or_unstable = false;
	result.first_zero_pivot = n;

	// Find diagonal positions
	result.diag_pos.assign(n, Index(-1));
	for (std::size_t i = 0; i < n; i++) {
		for (Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
			if (col_ind[static_cast<std::size_t>(p)] == static_cast<Index>(i)) {
				result.diag_pos[i] = p;
				break;
			}
		}
		if (result.diag_pos[i] == Index(-1)) ilu0_record_zero_pivot(result, i);
	}

	// col -> pos in row i (working array for fill lookup)
	// Since ILU(0) has no fill-in, we only update existing entries.
	std::vector<Index> col_to_pos(n, Index(-1));

	for (std::size_t i = 0; i < n; i++) {
		// Build col -> position map for row i
		for (Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
			col_to_pos[static_cast<std::size_t>(col_ind[static_cast<std::size_t>(p)])] = p;
		}

		// For each k < i where A[i][k] != 0 (lower triangle)
		for (Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
			const std::size_t k = static_cast<std::size_t>(col_ind[static_cast<std::size_t>(p)]);
			if (k >= i) break; // sorted, so we can break at diagonal

			// Check if diagonal of row k exists
			if (result.diag_pos[k] == Index(-1)) continue;
			const T u_kk = result.lu_val[static_cast<std::size_t>(result.diag_pos[k])];
			const R abs_ukk = tsparse_scalar::abs_value(u_kk);
			if (abs_ukk <= pivot_tol) {
				ilu0_record_zero_pivot(result, k);
				continue;
			}

			// l_ik = a_ik / u_kk
			result.lu_val[static_cast<std::size_t>(p)] /= u_kk;
			const T l_ik = result.lu_val[static_cast<std::size_t>(p)];

			// For each j > k where a_kj != 0 AND a_ij exists (ILU(0) pattern)
			for (Index q = result.diag_pos[k]; q < row_ptr[k + 1]; q++) {
				const std::size_t j = static_cast<std::size_t>(col_ind[static_cast<std::size_t>(q)]);
				const Index pos_ij = col_to_pos[j];
				if (pos_ij == Index(-1)) continue; // no fill in ILU(0)
				result.lu_val[static_cast<std::size_t>(pos_ij)] -= l_ik * result.lu_val[static_cast<std::size_t>(q)];
			}
		}

		// Reset col_to_pos for row i
		for (Index p = row_ptr[i]; p < row_ptr[i + 1]; p++) {
			col_to_pos[static_cast<std::size_t>(col_ind[static_cast<std::size_t>(p)])] = Index(-1);
		}
		if (result.diag_pos[i] != Index(-1)
		 && tsparse_scalar::abs_value(result.lu_val[static_cast<std::size_t>(result.diag_pos[i])]) <= pivot_tol) {
			ilu0_record_zero_pivot(result, i);
		}
	}

	result.factorized = true;
	result.diagnostics = ilu0_make_diagnostics(result);
	return result;
}

// ---------------------------------------------------------------------------
// Forward substitution: solve L z = b (L has implicit unit diagonal)
// ---------------------------------------------------------------------------
template <typename T, typename Index>
void ilu0_forward_solve(const ilu0_data<T, Index>& f, std::vector<T>& z) {
	const std::size_t n = f.n;
	for (std::size_t i = 0; i < n; i++) {
		T s = z[i];
		for (Index p = f.row_ptr[i]; p < f.diag_pos[i]; p++) {
			s -= f.lu_val[static_cast<std::size_t>(p)] * z[static_cast<std::size_t>(f.col_ind[static_cast<std::size_t>(p)])];
		}
		z[i] = s; // L has unit diagonal
	}
}

// ---------------------------------------------------------------------------
// Backward substitution: solve U x = z
// ---------------------------------------------------------------------------
template <typename T, typename Index>
void ilu0_backward_solve(const ilu0_data<T, Index>& f, std::vector<T>& x) {
	typedef typename tsparse_scalar::real_type<T>::type R;
	const std::size_t n = f.n;
	for (std::size_t ii = 0; ii < n; ii++) {
		const std::size_t i = n - 1 - ii;
		T s = x[i];
		const Index d = f.diag_pos[i];
		for (Index p = d + 1; p < f.row_ptr[i + 1]; p++) {
			s -= f.lu_val[static_cast<std::size_t>(p)] * x[static_cast<std::size_t>(f.col_ind[static_cast<std::size_t>(p)])];
		}
		const T u_ii = (d >= 0) ? f.lu_val[static_cast<std::size_t>(d)] : T(0);
		const R abs_d = tsparse_scalar::abs_value(u_ii);
		if (abs_d <= f.pivot_tol) {
			vcp::throw_error<vcp::numerical_error>(
				"tsparse_factorization::ilu0_backward_solve: zero or near-zero pivot");
		}
		x[i] = s / u_ii;
	}
}

// ---------------------------------------------------------------------------
// Apply ILU(0) preconditioner: solve (LU) x = b
// ---------------------------------------------------------------------------
template <typename T, typename Index>
std::vector<T> ilu0_solve(const ilu0_data<T, Index>& f, const std::vector<T>& b) {
	if (!f.factorized) {
		vcp::throw_error<vcp::state_error>("tsparse_factorization::ilu0_solve: factorization is not ready");
	}
	if (f.singular_or_unstable) {
		vcp::throw_error<vcp::numerical_error>("tsparse_factorization::ilu0_solve: ", f.diagnostics);
	}
	std::vector<T> x = b;
	ilu0_forward_solve(f, x);
	ilu0_backward_solve(f, x);
	return x;
}

template <class SparseMatrix>
class ilu0_factorization {
public:
	typedef typename SparseMatrix::value_type value_type;
	typedef typename SparseMatrix::index_type index_type;
	typedef typename tsparse_scalar::real_type<value_type>::type real_type;

	explicit ilu0_factorization(const SparseMatrix& A)
		: matrix_(A.as_csr()), data_(), pivot_tolerance_(tsparse_scalar::decimal_power_negative<real_type>(14)),
		  diagnostics_() {}

	void factorize() {
		data_ = ilu0_factorize<value_type, index_type>(
			matrix_.outer_index(), matrix_.inner_index(), matrix_.values(),
			static_cast<std::size_t>(matrix_.rowsize()), pivot_tolerance_);
		diagnostics_ = data_.diagnostics;
	}

	void solve(const std::vector<value_type>& b, std::vector<value_type>& x) const {
		if (!data_.factorized) {
			vcp::throw_error<vcp::state_error>("tsparse_factorization::ilu0_factorization::solve: factorization is not ready");
		}
		if (b.size() != static_cast<std::size_t>(matrix_.rowsize())) {
			vcp::throw_error<vcp::dimension_error>("tsparse_factorization::ilu0_factorization::solve: dimension mismatch");
		}
		if (data_.singular_or_unstable) {
			vcp::throw_error<vcp::numerical_error>(
				"tsparse_factorization::ilu0_factorization::solve: ", diagnostics_);
		}
		x = ilu0_solve(data_, b);
	}

	void apply(const std::vector<value_type>& b, std::vector<value_type>& x) const {
		solve(b, x);
	}

	bool is_factorized() const {
		return data_.factorized;
	}

	std::size_t zero_pivot_count() const {
		return data_.zero_pivots;
	}

	bool singular_or_unstable() const {
		return data_.singular_or_unstable;
	}

	std::size_t first_zero_pivot_index() const {
		return data_.first_zero_pivot;
	}

	real_type pivot_tolerance() const {
		return pivot_tolerance_;
	}

	std::string diagnostics() const {
		return diagnostics_;
	}

private:
	SparseMatrix matrix_;
	ilu0_data<value_type, index_type> data_;
	real_type pivot_tolerance_;
	std::string diagnostics_;
};

// ---------------------------------------------------------------------------
// Standalone GMRES (no dependency on spmatrix)
// Used internally for shift-invert inner solves
// ---------------------------------------------------------------------------
template <typename T, typename ApplyA, typename ApplyPrec>
struct gmres_result {
	std::vector<T> x;
	bool converged;
	std::size_t iterations;
	typename tsparse_scalar::real_type<T>::type residual_norm;
};

template <typename T, typename ApplyA, typename ApplyPrec>
gmres_result<T, ApplyA, ApplyPrec> gmres_solve(
	ApplyA apply_A,            // apply_A(x, y): y = A x
	ApplyPrec apply_prec,      // apply_prec(r, z): z = M^{-1} r
	const std::vector<T>& b,
	const std::size_t max_iter,
	const typename tsparse_scalar::real_type<T>::type& tol,
	const std::size_t restart)
{
	typedef typename tsparse_scalar::real_type<T>::type R;
	const std::size_t n = b.size();
	const std::size_t m = std::min(restart, n);
	gmres_result<T, ApplyA, ApplyPrec> result;
	result.x.assign(n, T(0));
	result.converged = false;
	result.iterations = 0;
	result.residual_norm = tsparse_scalar::real_norm_value(b);

	if (result.residual_norm <= tol) {
		result.converged = true;
		return result;
	}

	for (std::size_t outer = 0; outer < max_iter && !result.converged; outer++) {
		// Compute residual r = b - A x
		std::vector<T> Ax(n, T(0));
		apply_A(result.x, Ax);
		std::vector<T> r(n);
		for (std::size_t i = 0; i < n; i++) r[i] = b[i] - Ax[i];
		std::vector<T> z(n);
		apply_prec(r, z);

		const R beta = tsparse_scalar::real_norm_value(z);
		if (beta <= tol) { result.converged = true; break; }

		std::vector<std::vector<T> > V(m + 1, std::vector<T>(n, T(0)));
		std::vector<std::vector<T> > H(m + 1, std::vector<T>(m, T(0)));
		std::vector<R> cs(m, R(0)), sn(m, R(0));
		std::vector<T> g(m + 1, T(0));
		g[0] = T(beta);
		for (std::size_t i = 0; i < n; i++) V[0][i] = z[i] / T(beta);

		std::size_t k = 0;
		for (; k < m; k++) {
			std::vector<T> w(n), wz(n);
			apply_A(V[k], w);
			apply_prec(w, wz);
			for (std::size_t j = 0; j <= k; j++) {
				R hkj = tsparse_scalar::real_dot_value(wz, V[j]);
				H[j][k] = T(hkj);
				for (std::size_t i = 0; i < n; i++) wz[i] -= T(hkj) * V[j][i];
			}
			const R hnext = tsparse_scalar::real_norm_value(wz);
			H[k + 1][k] = T(hnext);
			if (hnext > std::numeric_limits<R>::epsilon()) {
				for (std::size_t i = 0; i < n; i++) V[k + 1][i] = wz[i] / T(hnext);
			}
			for (std::size_t j = 0; j < k; j++) {
				T hj = H[j][k], hj1 = H[j+1][k];
				H[j][k]   =  T(cs[j]) * hj + T(sn[j]) * hj1;
				H[j+1][k] = -T(sn[j]) * hj + T(cs[j]) * hj1;
			}
			const R h0 = tsparse_scalar::abs_value(H[k][k]);
			const R h1 = tsparse_scalar::abs_value(H[k+1][k]);
			const R rho = tsparse_scalar::hypot_value(h0, h1);
			if (rho <= std::numeric_limits<R>::epsilon()) { cs[k] = R(1); sn[k] = R(0); }
			else { cs[k] = tsparse_scalar::real_part(H[k][k]) / rho; sn[k] = tsparse_scalar::real_part(H[k+1][k]) / rho; }
			H[k][k] = T(rho);
			H[k+1][k] = T(0);
			const T gk = g[k];
			g[k]   = T(cs[k]) * gk;
			g[k+1] = T(-sn[k]) * gk;
			result.iterations++;
			result.residual_norm = tsparse_scalar::abs_value(g[k+1]);
			if (result.residual_norm <= tol) { k++; result.converged = true; break; }
		}
		// Back substitution for H(1:k,1:k) y = g(1:k)
		std::vector<T> y(k, T(0));
		for (std::size_t ii = 0; ii < k; ii++) {
			const std::size_t i = k - 1 - ii;
			T s = g[i];
			for (std::size_t j = i + 1; j < k; j++) s -= H[i][j] * y[j];
			if (tsparse_scalar::abs_value(H[i][i]) > std::numeric_limits<R>::epsilon())
				y[i] = s / H[i][i];
		}
		for (std::size_t j = 0; j < k; j++)
			for (std::size_t i = 0; i < n; i++) result.x[i] += V[j][i] * y[j];
	}
	return result;
}

} // namespace tsparse_factorization
} // namespace vcp

#endif
