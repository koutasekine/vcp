// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License

#pragma once

#ifndef VCP_TSPARSE_EIGENSOLVERS_HPP
#define VCP_TSPARSE_EIGENSOLVERS_HPP

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <string>
#include <vector>

#include <vcp/error.hpp>
#include <vcp/tsparse/tsparse_dense_linalg.hpp>
#include <vcp/tsparse/tsparse_eigen_selection.hpp>
#include <vcp/tsparse/tsparse_eigs.hpp>
#include <vcp/tsparse/tsparse_scalar.hpp>

namespace vcp {
	namespace tsparse_eigensolvers {
		template <typename T>
		struct lanczos_tridiagonalization {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			std::vector<std::vector<T> > basis;
			std::vector<T> alpha;
			std::vector<T> beta;
			real_type last_beta;
			std::size_t matrix_vector_products;
			std::string breakdown_reason;

			lanczos_tridiagonalization()
				: last_beta(real_type(0)), matrix_vector_products(0), breakdown_reason() {}
		};

		template <typename T>
		struct arnoldi_factorization {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			std::vector<std::vector<T> > basis;
			std::vector<std::vector<T> > hessenberg;
			std::size_t basis_size;
			real_type last_h;
			std::size_t matrix_vector_products;
			std::string breakdown_reason;

			arnoldi_factorization()
				: basis_size(0), last_h(real_type(0)), matrix_vector_products(0), breakdown_reason() {}
		};

		inline std::size_t krylov_subspace_dim(const std::size_t n, const std::size_t k, const std::size_t requested) {
			if (requested != 0) {
				if (requested <= k) {
					vcp::throw_error<vcp::invalid_argument>("tsparse_eigensolvers::krylov_subspace_dim: subspace_dim must be larger than k");
				}
				return std::min(n, requested);
			}
			return std::min(n, std::max<std::size_t>(2 * k + 10, std::max<std::size_t>(20, k + 1)));
		}

		inline std::size_t projected_max_iter(const std::size_t max_iter, const std::size_t dimension) {
			return std::max<std::size_t>(max_iter, 100 * dimension * dimension);
		}

		template <typename T>
		std::vector<std::vector<T> > lanczos_tridiagonal_matrix(const std::vector<T>& alpha, const std::vector<T>& beta) {
			const std::size_t m = alpha.size();
			std::vector<std::vector<T> > Tm(m, std::vector<T>(m, T(0)));
			for (std::size_t i = 0; i < m; i++) {
				Tm[i][i] = alpha[i];
				if (i + 1 < m) {
					Tm[i][i + 1] = beta[i];
					Tm[i + 1][i] = beta[i];
				}
			}
			return Tm;
		}

		template <typename T>
		std::vector<std::vector<T> > square_hessenberg(const std::vector<std::vector<T> >& H, const std::size_t m) {
			std::vector<std::vector<T> > Hm(m, std::vector<T>(m, T(0)));
			for (std::size_t col = 0; col < m; col++) {
				for (std::size_t row = 0; row < m; row++) Hm[row][col] = H[row][col];
			}
			return Hm;
		}

		template <typename T, class Apply>
		lanczos_tridiagonalization<T> build_lanczos_tridiagonalization(const std::size_t n,
		                                                               const std::size_t m_limit,
		                                                               const std::size_t max_iter,
		                                                               const typename tsparse_scalar::real_type<T>::type& tol,
		                                                               const bool reorthogonalize,
		                                                               Apply apply) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const real_type tolerance = tol;
			lanczos_tridiagonalization<T> decomp;
			decomp.basis.reserve(m_limit + 1);
			std::vector<T> q(n, T(0));
			for (std::size_t i = 0; i < n; i++) q[i] = T(i + 1);
			const real_type nq = tsparse_scalar::real_norm_value(q);
			if (nq <= (std::numeric_limits<real_type>::epsilon)()) {
				decomp.breakdown_reason = "initial vector has zero norm";
				return decomp;
			}
			for (std::size_t i = 0; i < n; i++) q[i] /= T(nq);
			std::vector<T> q_prev(n, T(0));
			real_type beta_prev(0);
			for (std::size_t iter = 0; iter < m_limit && iter < max_iter; iter++) {
				decomp.basis.push_back(q);
				std::vector<T> z;
				apply(q, z);
				decomp.matrix_vector_products++;
				if (iter > 0) {
					for (std::size_t i = 0; i < n; i++) z[i] -= T(beta_prev) * q_prev[i];
				}
				const real_type a = tsparse_scalar::real_dot_value(q, z);
				decomp.alpha.push_back(T(a));
				for (std::size_t i = 0; i < n; i++) z[i] -= T(a) * q[i];
				if (reorthogonalize) {
					for (std::size_t j = 0; j < decomp.basis.size(); j++) {
						const real_type h = tsparse_scalar::real_dot_value(z, decomp.basis[j]);
						for (std::size_t i = 0; i < n; i++) z[i] -= T(h) * decomp.basis[j][i];
					}
				}
				decomp.last_beta = tsparse_scalar::real_norm_value(z);
				if (decomp.last_beta <= tolerance) {
					decomp.breakdown_reason = "happy breakdown";
					break;
				}
				decomp.beta.push_back(T(decomp.last_beta));
				q_prev = q;
				for (std::size_t i = 0; i < n; i++) q[i] = z[i] / T(decomp.last_beta);
				beta_prev = decomp.last_beta;
			}
			return decomp;
		}

		template <typename T, class Apply>
		arnoldi_factorization<T> build_arnoldi_factorization(const std::size_t n,
		                                                     const std::size_t m_limit,
		                                                     const std::size_t max_iter,
		                                                     const typename tsparse_scalar::real_type<T>::type& tol,
		                                                     const bool reorthogonalize,
		                                                     const orthogonalization_method orthogonalization,
		                                                     Apply apply) {
			typedef typename tsparse_scalar::real_type<T>::type real_type;
			const real_type tolerance = tol;
			arnoldi_factorization<T> decomp;
			(void)orthogonalization;
			decomp.basis.assign(m_limit + 1, std::vector<T>(n, T(0)));
			decomp.hessenberg.assign(m_limit + 1, std::vector<T>(m_limit, T(0)));
			for (std::size_t i = 0; i < n; i++) decomp.basis[0][i] = T(i + 1);
			const real_type n0 = tsparse_scalar::real_norm_value(decomp.basis[0]);
			if (n0 <= (std::numeric_limits<real_type>::epsilon)()) {
				decomp.breakdown_reason = "initial vector has zero norm";
				return decomp;
			}
			for (std::size_t i = 0; i < n; i++) decomp.basis[0][i] /= T(n0);
			std::size_t m = 0;
			while (m < m_limit && m < max_iter) {
				std::vector<T> w;
				apply(decomp.basis[m], w);
				decomp.matrix_vector_products++;
				for (std::size_t j = 0; j <= m; j++) {
					decomp.hessenberg[j][m] = T(tsparse_scalar::real_dot_value(w, decomp.basis[j]));
					for (std::size_t i = 0; i < n; i++) w[i] -= decomp.hessenberg[j][m] * decomp.basis[j][i];
				}
				if (reorthogonalize) {
					for (std::size_t j = 0; j <= m; j++) {
						const real_type h2 = tsparse_scalar::real_dot_value(w, decomp.basis[j]);
						decomp.hessenberg[j][m] += T(h2);
						for (std::size_t i = 0; i < n; i++) w[i] -= T(h2) * decomp.basis[j][i];
					}
				}
				decomp.last_h = tsparse_scalar::real_norm_value(w);
				decomp.hessenberg[m + 1][m] = T(decomp.last_h);
				m++;
				if (decomp.last_h <= tolerance || m >= m_limit || m >= n) {
					if (decomp.last_h <= tolerance) decomp.breakdown_reason = "happy breakdown";
					break;
				}
				for (std::size_t i = 0; i < n; i++) decomp.basis[m][i] = w[i] / T(decomp.last_h);
			}
			decomp.basis_size = m;
			if (decomp.basis.size() > m + 1) decomp.basis.resize(m + 1);
			return decomp;
		}

		template <typename T>
		tsparse_dense_linalg::dense_eigen_result<T> extract_lanczos_ritz(const lanczos_tridiagonalization<T>& decomp,
		                                                                const std::size_t projected_iter,
		                                                                const typename tsparse_scalar::real_type<T>::type& tol) {
			const std::vector<std::vector<T> > Tm = lanczos_tridiagonal_matrix(decomp.alpha, decomp.beta);
			return tsparse_dense_linalg::jacobi_eig_dense(Tm, projected_iter, tol);
		}

		template <typename T>
		tsparse_dense_linalg::dense_eigen_result<T> extract_arnoldi_ritz(const arnoldi_factorization<T>& decomp,
		                                                                const std::size_t projected_iter,
		                                                                const typename tsparse_scalar::real_type<T>::type& tol) {
			const std::vector<std::vector<T> > Hm = square_hessenberg(decomp.hessenberg, decomp.basis_size);
			return tsparse_dense_linalg::is_dense_symmetric(Hm, tsparse_scalar::decimal_power_negative<typename tsparse_scalar::real_type<T>::type>(10))
				? tsparse_dense_linalg::jacobi_eig_dense(Hm, projected_iter, tol)
				: tsparse_dense_linalg::qr_eig_dense(Hm, projected_iter, tol);
		}

		// -----------------------------------------------------------------------
		// Francis double-shift QR step on upper Hessenberg H[0..n-1][0..n-1]
		// -----------------------------------------------------------------------
		template <typename T>
		void francis_qr_step(std::vector<std::vector<T> >& H, const std::size_t n) {
		typedef typename tsparse_scalar::real_type<T>::type R;
		if (n < 2) return;
		const R s1 = tsparse_scalar::real_part(H[n-1][n-1] + H[n-2][n-2]);
		const R s2 = tsparse_scalar::real_part(H[n-1][n-1] * H[n-2][n-2]
		                                      - H[n-1][n-2] * H[n-2][n-1]);
		R x = tsparse_scalar::real_part(H[0][0]*H[0][0] + H[0][1]*H[1][0]
		                                - s1*H[0][0] + s2);
		R y = tsparse_scalar::real_part(H[1][0] * (H[0][0] + H[1][1] - s1));
		R z = (n > 2) ? tsparse_scalar::real_part(H[1][0] * H[2][1]) : R(0);
		for (std::size_t k = 0; k + 1 < n; k++) {
			const std::size_t len = (k + 2 < n) ? 3 : 2;
			if (len < 3) z = R(0);
			const R nrm = (len == 3) ? tsparse_scalar::sqrt_value(x*x + y*y + z*z)
			                         : tsparse_scalar::sqrt_value(x*x + y*y);
			if (nrm <= std::numeric_limits<R>::epsilon()) break;
			const R sign_x = (x >= R(0)) ? R(1) : R(-1);
			R u0 = x + sign_x * nrm, u1 = y, u2 = z;
			const R un = (len == 3) ? tsparse_scalar::sqrt_value(u0*u0 + u1*u1 + u2*u2)
			                        : tsparse_scalar::sqrt_value(u0*u0 + u1*u1);
			if (un <= std::numeric_limits<R>::epsilon()) break;
			u0 /= un; u1 /= un; if (len == 3) u2 /= un; else u2 = R(0);
			// Apply from left
			const std::size_t jstart = (k > 0) ? k - 1 : 0;
			for (std::size_t j = jstart; j < n; j++) {
				const R d = u0*tsparse_scalar::real_part(H[k][j])
				          + u1*tsparse_scalar::real_part(H[k+1][j])
				          + ((len==3) ? u2*tsparse_scalar::real_part(H[k+2][j]) : R(0));
				H[k][j]   -= T(R(2)*d*u0);
				H[k+1][j] -= T(R(2)*d*u1);
				if (len == 3) H[k+2][j] -= T(R(2)*d*u2);
			}
			// Apply from right
			const std::size_t iend = std::min(k + len + 1, n);
			for (std::size_t i = 0; i < iend; i++) {
				const R d = u0*tsparse_scalar::real_part(H[i][k])
				          + u1*tsparse_scalar::real_part(H[i][k+1])
				          + ((len==3) ? u2*tsparse_scalar::real_part(H[i][k+2]) : R(0));
				H[i][k]   -= T(R(2)*d*u0);
				H[i][k+1] -= T(R(2)*d*u1);
				if (len == 3) H[i][k+2] -= T(R(2)*d*u2);
			}
			x = tsparse_scalar::real_part(H[k+1][k]);
			y = (k + 2 < n) ? tsparse_scalar::real_part(H[k+2][k]) : R(0);
			z = (k + 3 < n) ? tsparse_scalar::real_part(H[k+3][k]) : R(0);
		}
	}

	// -----------------------------------------------------------------------
	// Extract all eigenvalues (complex) from an upper Hessenberg matrix.
	// Never throws; returns complex<R> values.
	// -----------------------------------------------------------------------
	template <typename T>
	std::vector<std::complex<typename tsparse_scalar::real_type<T>::type> >
	hessenberg_complex_eigenvalues(std::vector<std::vector<T> > H,
	                               const std::size_t max_iter_per_dim,
	                               const typename tsparse_scalar::real_type<T>::type& tol)
	{
		typedef typename tsparse_scalar::real_type<T>::type R;
		typedef std::complex<R> C;
		const std::size_t n_orig = H.size();
		std::vector<C> result;
		result.reserve(n_orig);
		std::size_t active = n_orig;
		const std::size_t max_total = max_iter_per_dim * (n_orig + 1);
		std::size_t iters = 0;

		while (active > 0 && iters < max_total) {
			if (active == 1) {
				result.push_back(C(tsparse_scalar::real_part(H[0][0]), R(0)));
				break;
			}
			if (active == 2) {
				const R a = tsparse_scalar::real_part(H[0][0]);
				const R b = tsparse_scalar::real_part(H[0][1]);
				const R c = tsparse_scalar::real_part(H[1][0]);
				const R d = tsparse_scalar::real_part(H[1][1]);
				const R tr = a + d;
				const R det2 = a*d - b*c;
				const R disc = tr*tr - R(4)*det2;
				if (disc >= R(0)) {
					const R s = tsparse_scalar::sqrt_value(disc);
					result.push_back(C((tr + s) / R(2), R(0)));
					result.push_back(C((tr - s) / R(2), R(0)));
				} else {
					const R s = tsparse_scalar::sqrt_value(-disc);
					result.push_back(C(tr / R(2),  s / R(2)));
					result.push_back(C(tr / R(2), -s / R(2)));
				}
				break;
			}

			// Check bottom subdiagonal for deflation
			const R eps_b = tol * (tsparse_scalar::abs_value(H[active-1][active-1])
			                      + tsparse_scalar::abs_value(H[active-2][active-2]));
			const R sub_b = tsparse_scalar::abs_value(H[active-1][active-2]);
			if (sub_b <= eps_b || sub_b <= tol * tol) {
				result.push_back(C(tsparse_scalar::real_part(H[active-1][active-1]), R(0)));
				active--;
				H.resize(active);
				for (std::size_t i = 0; i < active; i++) H[i].resize(active);
				continue;
			}

			// Check second subdiagonal for 2x2 deflation
			if (active >= 3) {
				const R eps_c = tol * (tsparse_scalar::abs_value(H[active-2][active-2])
				                      + tsparse_scalar::abs_value(H[active-3][active-3]));
				const R sub_c = tsparse_scalar::abs_value(H[active-2][active-3]);
				if (sub_c <= eps_c || sub_c <= tol * tol) {
					const R a = tsparse_scalar::real_part(H[active-2][active-2]);
					const R b = tsparse_scalar::real_part(H[active-2][active-1]);
					const R c = tsparse_scalar::real_part(H[active-1][active-2]);
					const R d = tsparse_scalar::real_part(H[active-1][active-1]);
					const R tr = a + d;
					const R det2 = a*d - b*c;
					const R disc = tr*tr - R(4)*det2;
					if (disc >= R(0)) {
						const R s = tsparse_scalar::sqrt_value(disc);
						result.push_back(C((tr+s)/R(2), R(0)));
						result.push_back(C((tr-s)/R(2), R(0)));
					} else {
						const R s = tsparse_scalar::sqrt_value(-disc);
						result.push_back(C(tr/R(2),  s/R(2)));
						result.push_back(C(tr/R(2), -s/R(2)));
					}
					active -= 2;
					H.resize(active);
					for (std::size_t i = 0; i < active; i++) H[i].resize(active);
					continue;
				}
			}

			// Apply Francis double-shift QR step
			francis_qr_step(H, active);
			iters++;
		}

		// Fallback for remaining unconverged
		if (result.size() < n_orig) {
			for (std::size_t i = 0; i < active && result.size() < n_orig; i++) {
				result.push_back(C(tsparse_scalar::real_part(H[i][i]), R(0)));
			}
		}
		return result;
	}

	// -----------------------------------------------------------------------
	// Orthogonalize vector w against a set of orthonormal vectors Q.
	// Returns the orthogonalized w; coefficients (h) not returned here.
	// -----------------------------------------------------------------------
	template <typename T>
	void orthogonalize_cgs(const std::vector<std::vector<T> >& Q,
	                        const std::size_t m,
	                        std::vector<T>& w,
	                        std::vector<T>& h)
	{
		typedef typename tsparse_scalar::real_type<T>::type R;
		const std::size_t n = w.size();
		h.assign(m, T(0));
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] = T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
	}

	template <typename T>
	void orthogonalize_mgs(const std::vector<std::vector<T> >& Q,
	                        const std::size_t m,
	                        std::vector<T>& w,
	                        std::vector<T>& h)
	{
		typedef typename tsparse_scalar::real_type<T>::type R;
		const std::size_t n = w.size();
		h.assign(m, T(0));
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] += T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
	}

	template <typename T>
	void orthogonalize_cgs2(const std::vector<std::vector<T> >& Q,
	                         const std::size_t m,
	                         std::vector<T>& w,
	                         std::vector<T>& h)
	{
		typedef typename tsparse_scalar::real_type<T>::type R;
		const std::size_t n = w.size();
		h.assign(m, T(0));
		// First pass
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] = T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
		// Second pass
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] += T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
	}

	template <typename T>
	void orthogonalize_mgs2(const std::vector<std::vector<T> >& Q,
	                         const std::size_t m,
	                         std::vector<T>& w,
	                         std::vector<T>& h)
	{
		typedef typename tsparse_scalar::real_type<T>::type R;
		const std::size_t n = w.size();
		h.assign(m, T(0));
		// First pass (MGS)
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] = T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
		// Second pass (MGS again)
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] += T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
	}

	// DGKS: re-orthogonalize if ||w|| drops more than 1/sqrt(2)
	template <typename T>
	void orthogonalize_dgks(const std::vector<std::vector<T> >& Q,
	                         const std::size_t m,
	                         std::vector<T>& w,
	                         std::vector<T>& h)
	{
		typedef typename tsparse_scalar::real_type<T>::type R;
		const std::size_t n = w.size();
		h.assign(m, T(0));
		const R w_norm_before = tsparse_scalar::real_norm_value(w);
		// First pass (MGS)
		for (std::size_t j = 0; j < m; j++) {
			const R c = tsparse_scalar::real_dot_value(Q[j], w);
			h[j] = T(c);
			for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
		}
		const R w_norm_after = tsparse_scalar::real_norm_value(w);
		// Reorthogonalize if significant loss of orthogonality
		const R threshold = R(1) / tsparse_scalar::sqrt_value(R(2));
		if (w_norm_after < threshold * w_norm_before) {
			for (std::size_t j = 0; j < m; j++) {
				const R c = tsparse_scalar::real_dot_value(Q[j], w);
				h[j] += T(c);
				for (std::size_t i = 0; i < n; i++) w[i] -= T(c) * Q[j][i];
			}
		}
	}

	// Dispatch to the selected orthogonalization method
	template <typename T>
	void orthogonalize(const std::vector<std::vector<T> >& Q,
	                   const std::size_t m,
	                   std::vector<T>& w,
	                   std::vector<T>& h,
	                   const orthogonalization_method method)
	{
		switch (method) {
		case orthogonalization_method::classical_gram_schmidt_twice:
			orthogonalize_cgs2(Q, m, w, h); return;
		case orthogonalization_method::modified_gram_schmidt:
			orthogonalize_mgs(Q, m, w, h); return;
		}
		orthogonalize_mgs(Q, m, w, h);
	}

	// Select Ritz indices by target (works on complex eigenvalues)
	// Returns sorted indices (best first)
	template <typename T>
	std::vector<std::size_t> select_ritz_indices(
		const std::vector<std::complex<typename tsparse_scalar::real_type<T>::type> >& eigs,
		const std::size_t k,
		const eig_target target,
		const typename tsparse_scalar::real_type<T>::type& shift_val)
	{
		return vcp::tsparse_eigen_selection::select_complex_ritz_values(
			eigs, k, target, shift_val);
		}

	} // namespace tsparse_eigensolvers
} // namespace vcp

#endif
