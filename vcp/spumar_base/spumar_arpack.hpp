// VCP Library
// http ://verified.computation.jp
//
// VCP Library is licensed under the BSD 3 - clause "New" or "Revised" License
//
// spumar_base/spumar_arpack.hpp
// Single reverse-communication driver for ARPACK-NG (dsaupd/dseupd,
// dnaupd/dneupd) with injected OP / B callbacks (design §3: one loop
// implementation; the standard and generalized delegations differ only in
// bmat / mode / the injected operators).
//
// Determinism: the start vector is caller-supplied (info = 1) with a fixed
// recipe (EIG-5 harness precedent), so repeated runs are reproducible.
//
// This header is part of the spumar delegation layer: it declares the
// external ARPACK Fortran symbols.  No vcp core header may include it
// (dependency isolation, spumar design v1.1 G4 / B-2).

#pragma once

#ifndef VCP_SPUMAR_ARPACK_HPP
#define VCP_SPUMAR_ARPACK_HPP

#include <cstddef>
#include <functional>
#include <string>
#include <vector>

extern "C" {
void dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
             double* resid, int* ncv, double* v, int* ldv, int* iparam,
             int* ipntr, double* workd, double* workl, int* lworkl,
             int* info);
void dseupd_(int* rvec, char* howmny, int* select, double* d, double* z,
             int* ldz, double* sigma, char* bmat, int* n, char* which,
             int* nev, double* tol, double* resid, int* ncv, double* v,
             int* ldv, int* iparam, int* ipntr, double* workd,
             double* workl, int* lworkl, int* ierr);
void dnaupd_(int* ido, char* bmat, int* n, char* which, int* nev, double* tol,
             double* resid, int* ncv, double* v, int* ldv, int* iparam,
             int* ipntr, double* workd, double* workl, int* lworkl,
             int* info);
void dneupd_(int* rvec, char* howmny, int* select, double* dr, double* di,
             double* z, int* ldz, double* sigmar, double* sigmai,
             double* workev, char* bmat, int* n, char* which, int* nev,
             double* tol, double* resid, int* ncv, double* v, int* ldv,
             int* iparam, int* ipntr, double* workd, double* workl,
             int* lworkl, int* ierr);
}

namespace vcp {
namespace spumar_detail {

	// OP callback: y = OP * x.  bx_hint is non-null ONLY in generalized
	// mode 3 with ido == 1, where ARPACK has already placed B*x at
	// workd(ipntr(3)) — the operator may consume it instead of re-applying B.
	typedef std::function<void(const double* x, double* y, const double* bx_hint)> arpack_op_fn;
	// B callback: y = B * x (bmat='G' only; never called for bmat='I').
	typedef std::function<void(const double* x, double* y)> arpack_b_fn;

	struct arpack_outcome {
		int aupd_info;   // dsaupd/dnaupd final info
		int eupd_ierr;   // dseupd/dneupd ierr (0 if eupd not reached)
		int nconv;       // iparam[4]: converged Ritz values
		int niter;       // iparam[2]: actual Arnoldi/Lanczos update iterations
		int nopx;        // iparam[8]: number of OP*x operations
		int ncv_used;
		std::vector<int> iparam;              // full snapshot (S-5)
		// eigenvalues: symmetric -> dr only (ascending, dseupd order);
		// nonsymmetric -> dr/di adjacent conjugate pairs (dneupd order).
		std::vector<double> dr, di;
		// eigenvectors as columns; nonsymmetric complex pairs follow the
		// dneupd convention: columns (j, j+1) hold Re / Im of the (+im)
		// eigenvector.
		std::vector<std::vector<double> > z;

		arpack_outcome()
			: aupd_info(0), eupd_ierr(0), nconv(0), niter(0), nopx(0), ncv_used(0) {}

		bool aupd_ok() const { return aupd_info == 0; }
		bool honest_nonconv() const { return aupd_info == 1; }
	};

	// -------------------------------------------------------------------
	// The single reverse-communication loop (OP / B injection).
	// aupd = dsaupd_ or dnaupd_ (identical C signatures).
	// -------------------------------------------------------------------
	typedef void (*aupd_fn)(int*, char*, int*, char*, int*, double*, double*,
		int*, double*, int*, int*, int*, double*, double*, int*, int*);

	inline int arpack_aupd_loop_(aupd_fn aupd, char bmat, int n,
		const char* which2, int nev, double tol, int ncv, int mode,
		int max_iter,
		std::vector<double>& resid, std::vector<double>& v,
		std::vector<int>& iparam, std::vector<int>& ipntr,
		std::vector<double>& workd, std::vector<double>& workl,
		const arpack_op_fn& opx, const arpack_b_fn& bx)
	{
		int ido = 0;
		int ldv = n;
		int lworkl = static_cast<int>(workl.size());
		int info = 1; // caller-supplied deterministic start vector
		char which[3] = {which2[0], which2[1], '\0'};

		resid.assign(static_cast<std::size_t>(n), 0.0);
		for (int i = 0; i < n; i++) {
			resid[static_cast<std::size_t>(i)] =
				1.0 + 0.01 * static_cast<double>((i * 17) % 101);
		}

		iparam.assign(11, 0);
		iparam[0] = 1;        // exact shifts
		iparam[2] = max_iter; // max update iterations
		iparam[6] = mode;     // 1 = regular, 3 = shift-invert
		ipntr.assign(14, 0);

		while (true) {
			aupd(&ido, &bmat, &n, which, &nev, &tol, resid.data(), &ncv,
			     v.data(), &ldv, iparam.data(), ipntr.data(), workd.data(),
			     workl.data(), &lworkl, &info);
			if (ido == -1 || ido == 1) {
				const double* x = &workd[static_cast<std::size_t>(ipntr[0] - 1)];
				double* y = &workd[static_cast<std::size_t>(ipntr[1] - 1)];
				const double* bx_hint =
					(bmat == 'G' && mode == 3 && ido == 1)
						? &workd[static_cast<std::size_t>(ipntr[2] - 1)]
						: static_cast<const double*>(0);
				opx(x, y, bx_hint);
			} else if (ido == 2) {
				const double* x = &workd[static_cast<std::size_t>(ipntr[0] - 1)];
				double* y = &workd[static_cast<std::size_t>(ipntr[1] - 1)];
				bx(x, y);
			} else {
				break;
			}
		}
		return info;
	}

	// -------------------------------------------------------------------
	// Symmetric drive: dsaupd loop + dseupd extraction.
	// -------------------------------------------------------------------
	inline arpack_outcome arpack_drive_symmetric(int n, int nev, int ncv,
		const char* which, char bmat, int mode, double sigma, double tol,
		int max_iter, const arpack_op_fn& opx, const arpack_b_fn& bx)
	{
		arpack_outcome out;
		out.ncv_used = ncv;
		std::vector<double> resid, workd(static_cast<std::size_t>(3) * n, 0.0);
		std::vector<double> v(static_cast<std::size_t>(n) * ncv, 0.0);
		std::vector<double> workl(static_cast<std::size_t>(ncv) * (ncv + 8), 0.0);
		std::vector<int> iparam, ipntr;

		out.aupd_info = arpack_aupd_loop_(&dsaupd_, bmat, n, which, nev, tol,
			ncv, mode, max_iter, resid, v, iparam, ipntr, workd, workl, opx, bx);
		out.iparam = iparam;
		out.nconv = iparam[4];
		out.niter = iparam[2];
		out.nopx = iparam[8];
		// info = 0: success; info = 1: honest non-convergence (nconv values
		// are still extractable).  Other codes: no extraction.
		if (out.aupd_info != 0 && out.aupd_info != 1) return out;
		if (out.nconv <= 0) return out;

		int rvec = 1;
		char howmny = 'A';
		char bmat_c = bmat;
		char which_c[3] = {which[0], which[1], '\0'};
		int ldv = n;
		int lworkl = static_cast<int>(workl.size());
		std::vector<int> select(static_cast<std::size_t>(ncv), 0);
		std::vector<double> d(static_cast<std::size_t>(nev), 0.0);
		std::vector<double> z(static_cast<std::size_t>(n) * nev, 0.0);
		double sig = sigma;
		int nev_c = nev;
		double tol_c = tol;
		dseupd_(&rvec, &howmny, select.data(), d.data(), z.data(), &ldv, &sig,
		        &bmat_c, &n, which_c, &nev_c, &tol_c, resid.data(), &ncv,
		        v.data(), &ldv, iparam.data(), ipntr.data(), workd.data(),
		        workl.data(), &lworkl, &out.eupd_ierr);
		if (out.eupd_ierr != 0) return out;

		const int m = out.nconv < nev ? out.nconv : nev;
		out.dr.assign(d.begin(), d.begin() + m);
		out.di.assign(static_cast<std::size_t>(m), 0.0);
		out.z.resize(static_cast<std::size_t>(m));
		for (int j = 0; j < m; j++) {
			out.z[static_cast<std::size_t>(j)].assign(
				z.begin() + static_cast<std::ptrdiff_t>(j) * n,
				z.begin() + static_cast<std::ptrdiff_t>(j + 1) * n);
		}
		return out;
	}

	// -------------------------------------------------------------------
	// Nonsymmetric drive: dnaupd loop + dneupd extraction.
	// dneupd may return nev+1 values when a conjugate pair straddles the
	// nev boundary; dr/di/z are sized accordingly and the caller decides.
	// -------------------------------------------------------------------
	inline arpack_outcome arpack_drive_nonsymmetric(int n, int nev, int ncv,
		const char* which, char bmat, int mode, double sigma_r, double sigma_i,
		double tol, int max_iter, const arpack_op_fn& opx, const arpack_b_fn& bx)
	{
		arpack_outcome out;
		out.ncv_used = ncv;
		std::vector<double> resid, workd(static_cast<std::size_t>(3) * n, 0.0);
		std::vector<double> v(static_cast<std::size_t>(n) * ncv, 0.0);
		std::vector<double> workl(static_cast<std::size_t>(3) * ncv * ncv + 6 * ncv, 0.0);
		std::vector<int> iparam, ipntr;

		out.aupd_info = arpack_aupd_loop_(&dnaupd_, bmat, n, which, nev, tol,
			ncv, mode, max_iter, resid, v, iparam, ipntr, workd, workl, opx, bx);
		out.iparam = iparam;
		out.nconv = iparam[4];
		out.niter = iparam[2];
		out.nopx = iparam[8];
		if (out.aupd_info != 0 && out.aupd_info != 1) return out;
		if (out.nconv <= 0) return out;

		int rvec = 1;
		char howmny = 'A';
		char bmat_c = bmat;
		char which_c[3] = {which[0], which[1], '\0'};
		int ldv = n;
		int lworkl = static_cast<int>(workl.size());
		std::vector<int> select(static_cast<std::size_t>(ncv), 0);
		std::vector<double> dr(static_cast<std::size_t>(nev) + 1, 0.0);
		std::vector<double> di(static_cast<std::size_t>(nev) + 1, 0.0);
		std::vector<double> z(static_cast<std::size_t>(n) * (nev + 1), 0.0);
		std::vector<double> workev(static_cast<std::size_t>(3) * ncv, 0.0);
		double sigr = sigma_r, sigi = sigma_i;
		int nev_c = nev;
		double tol_c = tol;
		dneupd_(&rvec, &howmny, select.data(), dr.data(), di.data(), z.data(),
		        &ldv, &sigr, &sigi, workev.data(), &bmat_c, &n, which_c,
		        &nev_c, &tol_c, resid.data(), &ncv, v.data(), &ldv,
		        iparam.data(), ipntr.data(), workd.data(), workl.data(),
		        &lworkl, &out.eupd_ierr);
		if (out.eupd_ierr != 0) return out;

		// dneupd reports the possibly-extended count in iparam[4].
		int m = iparam[4];
		if (m > nev + 1) m = nev + 1;
		if (m < 0) m = 0;
		out.nconv = m;
		out.dr.assign(dr.begin(), dr.begin() + m);
		out.di.assign(di.begin(), di.begin() + m);
		out.z.resize(static_cast<std::size_t>(m));
		for (int j = 0; j < m; j++) {
			out.z[static_cast<std::size_t>(j)].assign(
				z.begin() + static_cast<std::ptrdiff_t>(j) * n,
				z.begin() + static_cast<std::ptrdiff_t>(j + 1) * n);
		}
		return out;
	}

} // namespace spumar_detail
} // namespace vcp

#endif // VCP_SPUMAR_ARPACK_HPP
