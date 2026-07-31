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

// ---------------------------------------------------------------------------
// vcp/spnewton.hpp
//
// Newton's method for problems whose Jacobian is SPARSE.
//
// Relation to vcp/newton.hpp
// --------------------------
// vcp::Newton<_T,_PM> keeps the Jacobian in a dense vcp::matrix<_T,_PM> and
// solves the correction equation with the dense lss().  vcp::SpNewton keeps
// the UNKNOWN in the same dense vcp::matrix<_T,_PM> (an n x 1 column) but the
// JACOBIAN in vcp::spmatrix<_T,_SP>.  Everything else -- the iteration body,
// the correction-based stopping rule, the virtual f() / Df() /
// setting_newton() protocol -- is identical to vcp::Newton on purpose, so
// that porting a dense problem to the sparse solver is a mechanical change of
// the Df() return type.
//
// This is a SEPARATE header rather than an extension of newton.hpp because
//  (i)  the virtual signature of Df() differs, so unifying the two would be a
//       source-breaking change for every existing derived class, and
//  (ii) newton.hpp today depends only on <limits> and vcp_metafunction.hpp;
//       pulling spmatrix.hpp / spmats.hpp into it would make every dense-only
//       user pay for the sparse stack.
//
// Full Newton, not chord
// ----------------------
// The Jacobian is rebuilt AND refactorized at every step.  No LU handle is
// cached.  Freezing the Jacobian (the chord / modified Newton method) would
// tie the region of convergence to the initial point and would destroy the
// quadratic convergence that doubles as a self-check on the user's Df()
// implementation.  The correction-based stopping rule inherited from
// vcp::Newton is only meaningful under quadratic convergence, so this choice
// also keeps that rule valid.
//
// Linear solver
// -------------
// The default method is linear_solver_method::sparse_lu (a DIRECT solve),
// deliberately not auto_select: auto_select resolves a symmetric matrix to
// conjugate_gradient, and a Newton Jacobian is frequently symmetric but
// indefinite, and in any case ill-conditioned enough that CG stalls.  The
// options object is exposed through setting_linear_solve_options() for
// callers who know better.
//
// The dense <-> std::vector bridge
// --------------------------------
// vcp::spmatrix's solve() is std::vector<_T> in / std::vector<_T> out.  The
// right-hand side needs no copy at all (matrix::vecpointer() returns the
// underlying storage by const reference); only the correction has to be
// written back, through matrix::data(), with std::copy.  Both accessors are
// public members of vcp::matrix and work for every policy derived from
// vcp::mats<_T>.
// ---------------------------------------------------------------------------

#pragma once

#ifndef VCP_SPNEWTON_HPP
#define VCP_SPNEWTON_HPP

#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>

#include <vcp/vcp_metafunction.hpp>
#include <vcp/error.hpp>
#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>
#include <vcp/spmats.hpp>

namespace vcp {

    template < typename _T, typename _PM, typename _SP = vcp::spmats< _T > >
    class SpNewton {
    protected:
        bool flag_Convergence;
        int newton_max_iteration;
        int iteration_newton;
        _T newton_tol;
        _T Correction_term;
        vcp::linear_solve_options< _T > linear_solve_option;

    public:
        typedef vcp::matrix< _T, _PM >   vector_type;
        typedef vcp::spmatrix< _T, _SP > spmatrix_type;

        SpNewton() {
            this->newton_max_iteration = 100;
            this->flag_Convergence = false;
            this->newton_tol = _T(4) * std::numeric_limits< _T >::epsilon();
            // direct solve by default (see header note)
            this->linear_solve_option.method = vcp::linear_solver_method::sparse_lu;
        }

        // SpNewton is used through virtual f() / Df(); unlike vcp::Newton it
        // declares a virtual destructor so that deletion through a base
        // pointer is well defined.
        virtual ~SpNewton() {}

        void setting_newton_tol(const int n) {
            this->newton_tol = _T(n) * std::numeric_limits< _T >::epsilon();
        }

        void setting_newton_max_iteration(const int n) {
            this->newton_max_iteration = n;
        }

        // ---- sparse-specific: linear solver control -------------------
        void setting_linear_solve_options(const vcp::linear_solve_options< _T >& options) {
            this->linear_solve_option = options;
        }

        const vcp::linear_solve_options< _T >& get_linear_solve_options() const {
            return this->linear_solve_option;
        }

        bool is_convergence() {
            return this->flag_Convergence;
        }

        bool is_convergence(void) const {
            return this->flag_Convergence;
        }

        // ---- user hooks (same protocol as vcp::Newton) ----------------
        // setting_newton is called with the current iterate before every
        // f() / Df() pair; store whatever the two need.
        virtual void setting_newton(vector_type& xx) {
        }

        // residual F(x), returned as an n x 1 column
        virtual vector_type f() {
            vector_type x;
            x.zeros(1);
            return x;
        }

        // Jacobian F'(x), returned as an n x n sparse matrix.  It may be left
        // in the COO stage; solve_nls finalizes its own copy.
        virtual spmatrix_type Df() {
            spmatrix_type A(1, 1);
            return A;
        }

        // ---- the iteration -------------------------------------------
        vector_type solve_nls(const vector_type& uh) {
            vector_type x_old, x_new, fdi_fx;
            x_new = uh;
            this->flag_Convergence = false;
            using std::abs;
            for (iteration_newton = 0; iteration_newton < this->newton_max_iteration; iteration_newton++) {
                x_old = x_new;

                this->setting_newton(x_new);
                vector_type    fuh = this->f();
                spmatrix_type  J   = this->Df();

                fdi_fx = this->solve_correction(J, fuh);
                x_new = x_old - fdi_fx;
                vector_type s = max(abs(fdi_fx));
                Correction_term = s(0);
                s = max(abs(x_old));
                Correction_term /= s(0);
                if (Correction_term <= this->newton_tol) {
                    this->flag_Convergence = true;
                    this->disp_convergence();
                    return x_new;
                }
                this->disp_continue();
            }
            this->flag_Convergence = false;
            std::cout << "Not Convergence : i = " << this->newton_max_iteration << ", " << Correction_term << " <= " << this->newton_tol << std::endl;
            return x_new;
        }

        virtual void disp_convergence() {
            std::cout << "Convergence : i = " << iteration_newton << ", " << Correction_term << " <= " << this->newton_tol << std::endl;
        }

        virtual void disp_continue() {
            std::cout << "i = " << iteration_newton << ", " << Correction_term << " <= " << this->newton_tol << std::endl;
        }

    protected:
        // The single point where SpNewton differs from vcp::Newton:
        // J \ b with a sparse J.  J is taken by (non-const) reference to a
        // local copy owned by solve_nls, so finalize() here never touches an
        // object the caller still holds.
        vector_type solve_correction(spmatrix_type& J, const vector_type& b) {
            const int n = b.rowsize();
            if (b.columnsize() != 1) {
                vcp::throw_error< vcp::invalid_argument >(
                    "vcp::SpNewton::solve_nls", ": f() must return an n x 1 column vector");
            }
            if (J.rowsize() != n || J.columnsize() != n) {
                vcp::throw_error< vcp::invalid_argument >(
                    "vcp::SpNewton::solve_nls", ": Df() must return an n x n matrix matching f()");
            }
            // spmatrix::solve accepts a COO-stage matrix, but only mul_vec /
            // trans_mul_vec enforce the finalized state, so finalizing the
            // local copy here makes the state deterministic for any policy _SP.
            J.finalize();

            // right-hand side: no copy (vecpointer returns the storage itself)
            const std::vector< _T > d = J.solve(b.vecpointer(), this->linear_solve_option);

            vector_type x;
            x.zeros(n, 1);
            std::copy(d.begin(), d.end(), x.data());
            return x;
        }
    };

}

#endif // VCP_SPNEWTON_HPP
