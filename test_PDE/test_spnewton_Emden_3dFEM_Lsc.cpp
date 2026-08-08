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
// test_spnewton_Emden_3dFEM_Lsc.cpp
//
// 3D counterpart of test_spnewton_Emden_2dFEM.cpp.
//
//   - Delta u = u^2   in Omega,      u = 0 on  d Omega
//
// Omega is an L-shaped COLUMN ("Lsc" = L-shaped section column): the 2D
// L-shaped cross-section of the 2D sample, extruded along z.  Its re-entrant
// feature is an EDGE running parallel to z, not a corner.  The Fichera corner
// (a re-entrant vertex) cannot be produced by vcp::bfem::meshgen3, which only
// meshes extruded domains -- that is the separate MG-4 track.
//
// Everything above the mesh is dimension-generic, so the differences from the
// 2D sample are exactly three:
//     meshgen.hpp      -> meshgen3.hpp
//     polygon_domain   -> extruded_domain   (cross-section + [z0, z1])
//     fe_space<2,...>  -> fe_space<3,...>
//
// Weak form, with V_h = P^p functions vanishing on the boundary:
//     F(u)_i = ( grad u_h, grad psi_i ) - ( u_h^2, psi_i ) = 0
//     F'(u)  = ( grad psi_j, grad psi_i ) - ( 2 u_h psi_j, psi_i )
// The Newton iteration is vcp::SpNewton::solve_nls (sparse Jacobian).
//
// Build (headers only):
//     g++ -O2 -std=c++11 -I<path to vcp root> test_spnewton_Emden_3dFEM_Lsc.cpp
// add -fopenmp to use the parallel assembly paths of fe_space.
//
// SIZE WARNING.  At the delivered defaults (p = 3, h = 2^-3) this problem has
// 852,815 unknowns.  3D sparse LU fill-in is far heavier than in 2D: an
// extrapolation from measured runs at 1,265 and 11,891 unknowns puts the
// factorization at roughly 40 GB and several hours per Newton step, single
// threaded.  That extrapolation is NOT a measurement.  Lower mesh_level for a
// machine that cannot hold it; see the notes at the parameters below.
// ---------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <algorithm>

#include <vcp/mats.hpp>
#include <vcp/matrix.hpp>
#include <vcp/matrix_assist.hpp>
#include <vcp/spmats.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/spnewton.hpp>

#include <vcp/bfem/meshgen3.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/graphics.hpp>

#include <vcp/vcp_timer.hpp>

template < typename _T, typename _PM, typename _SP >
struct EMDEN3DFEM : public vcp::SpNewton< _T, _PM, _SP > {

    typedef vcp::bfem::fe_space< 3, _T, _PM, _SP >          space_type;
    typedef vcp::bfem::dirichlet_reduction< _T, _PM, _SP >  reduction_type;
    typedef typename space_type::function_type              function_type;

    int p;                                   // FEM degree: the k of P^k
    _T  h;                                   // mesh size
    int Number_of_elements;
    int Number_of_dof;

    std::unique_ptr< space_type >     Vh;    // fe_space (no default ctor)
    std::unique_ptr< reduction_type > Dirichlet;

    vcp::spmatrix< _T, _SP > DL;             // reduced stiffness ( grad psi_j, grad psi_i )
    vcp::matrix< _T, _PM >   u;              // current iterate, reduced
    vcp::matrix< _T, _PM >   u_full;         // same, expanded with zeros on the boundary

    vcp::bfem::poly1< _T > f_u;              // f(s)  = s^2
    vcp::bfem::poly1< _T > df_u;             // f'(s) = 2 s

    EMDEN3DFEM()
        : p(0), h(_T(0)), Number_of_elements(0), Number_of_dof(0),
          f_u(make_square()), df_u(make_double()) {}

    // ---- vcp::SpNewton hooks -------------------------------------------
    void setting_newton(vcp::matrix< _T, _PM >& uh) override {
        u = uh;
        u_full = Dirichlet->expand(uh);
    }

    // F(u) = DL * u - ( u_h^2, psi_i )
    vcp::matrix< _T, _PM > f() override {
        vcp::matrix< _T, _PM > DLu;
        DLu.zeros(u.rowsize(), 1);
        const std::vector< _T > t = DL.mul_vec(u.vecpointer());
        std::copy(t.begin(), t.end(), DLu.data());

        vcp::matrix< _T, _PM > uh2phi =
            Dirichlet->reduce(Vh->load(f_u, current_function(), p));
        return DLu - uh2phi;
    }

    // F'(u) = DL - ( 2 u_h psi_j, psi_i )
    // (the factor 2 is already carried by df_u, hence no "2 *" here)
    vcp::spmatrix< _T, _SP > Df() override {
        vcp::spmatrix< _T, _SP > uhphiphi =
            Dirichlet->reduce(Vh->weighted_mass(df_u, current_function(), p));
        uhphiphi.finalize();
        return DL - uhphiphi;
    }

    // ---- set-up ---------------------------------------------------------
    // degree : the k of P^k
    // mesh_h : every edge of the tetrahedral mesh satisfies len^2 <= mesh_h^2
    // z0, z1 : the extrusion interval
    void first_execute(const int degree, const _T& mesh_h,
                       const _T& z0, const _T& z1) {
        p = degree;
        h = mesh_h;

        // 1. L-shaped column: the cross-section is given by its vertices only
        //    (counter-clockwise), extruded along z
        vcp::bfem::extruded_domain< _T > Omega;
        Omega.base.outer.push_back(vertex(-1, -1));
        Omega.base.outer.push_back(vertex( 1, -1));
        Omega.base.outer.push_back(vertex( 1,  0));
        Omega.base.outer.push_back(vertex( 0,  0));
        Omega.base.outer.push_back(vertex( 0,  1));
        Omega.base.outer.push_back(vertex(-1,  1));
        Omega.z0 = z0;
        Omega.z1 = z1;

        // 2. tetrahedral mesh (bfem)
        vcp::bfem::meshgen_status< 3 > status;
        vcp::bfem::mesh< 3, _T > Th = vcp::bfem::generate_mesh(Omega, h, status);

        // 3. P^p space and homogeneous Dirichlet reduction
        Vh.reset(new space_type(Th, p));
        Number_of_elements = Vh->num_elements();
        Number_of_dof      = Vh->ndof(p);
        Dirichlet.reset(new reduction_type(Number_of_dof,
                                           Vh->dofs(p).boundary_dofs()));

        // 4. stiffness (assembled once: it does not depend on u)
        DL = Dirichlet->reduce(Vh->stiffness(p));
        DL.finalize();
    }

    int reduced_size() const { return Dirichlet->reduced_size(); }

    // u_h expanded to the full space, as an fe_function
    function_type current_function() {
        return Vh->function_from_coeffs(p, u_full);
    }

    // Initial guess shaped like the torsion function w  ( - Lap w = 1,
    // w = 0 on the boundary ), rescaled so that max w = peak.  A CONSTANT
    // initial value does not work in 3D: every constant tried either collapses
    // onto the trivial solution u = 0 or diverges (see the note in main).
    // Costs one extra sparse factorization of DL.
    vcp::matrix< _T, _PM > torsion_initial_value(const _T& peak) {
        std::vector< _T > one(1, _T(1));
        vcp::bfem::poly1< _T > constant_one =
            vcp::bfem::poly1< _T >::from_coeffs(one);
        vcp::matrix< _T, _PM > F =
            Dirichlet->reduce(Vh->load(constant_one, Vh->zero_function(p), p));

        vcp::linear_solve_options< _T > options;
        options.method = vcp::linear_solver_method::sparse_lu;
        const std::vector< _T > wv = DL.solve(F.vecpointer(), options);

        vcp::matrix< _T, _PM > w;
        w.zeros(reduced_size(), 1);
        std::copy(wv.begin(), wv.end(), w.data());
        return (peak / max(abs(w))(0)) * w;
    }

private:
    static std::array< _T, 2 > vertex(const int x, const int y) {
        std::array< _T, 2 > v;
        v[0] = _T(x);
        v[1] = _T(y);
        return v;
    }
    static vcp::bfem::poly1< _T > make_square() {          // s^2
        std::vector< _T > a(3, _T(0));
        a[2] = _T(1);
        return vcp::bfem::poly1< _T >::from_coeffs(a);
    }
    static vcp::bfem::poly1< _T > make_double() {          // 2 s
        std::vector< _T > a(2, _T(0));
        a[1] = _T(2);
        return vcp::bfem::poly1< _T >::from_coeffs(a);
    }
};

int main(void) {
    typedef double                      TYPE;
    typedef vcp::mats< TYPE >           POLICY;      // dense policy  (vectors)
    typedef vcp::spmats< TYPE >         SPPOLICY;    // sparse policy (Jacobian)

    // ---- discretization parameters (all variable) -----------------------
    const int  fem_degree = 3;                       // P^k : k = 3
    const TYPE mesh_size  = TYPE(1) / TYPE(8);       // h   = 2^-3
    const TYPE z_bottom   = TYPE(0);                 // extrusion interval
    const TYPE z_top      = TYPE(1);

    // ---- initial value ---------------------------------------------------
    // Peak of the torsion-shaped initial guess.  The basin of attraction of
    // the positive solution is NARROW in 3D: on the coarsest mesh, peaks of
    // 12 and 16 and 25 all fail (12 falls onto u = 0, the others diverge)
    // while 20 converges.  20 also converges one refinement finer.  It has
    // NOT been checked at the defaults above; retune if Newton fails.
    const TYPE initial_peak = TYPE(20);

    // ---- convergence tolerance -------------------------------------------
    // SpNewton's default 4 * eps is a LOWER BOUND on what any problem can
    // reach; the attainable tolerance is a property of the problem and belongs
    // here.  Measured floors: 64 suffices at ~1.2e4 unknowns, and the 2D
    // sample needed 512 at ~2.9e5.  The value below is an EXTRAPOLATION to
    // ~8.5e5 unknowns and is NOT verified; raise it if Newton stalls just
    // above the tolerance, lower it for coarser meshes.
    const int  newton_tol_factor = 1024;

    EMDEN3DFEM< TYPE, POLICY, SPPOLICY > N;

    vcp::time.tic();
    N.first_execute(fem_degree, mesh_size, z_bottom, z_top);
    std::cout << "L-shaped column,  P^" << fem_degree << ",  h = " << mesh_size
              << ",  z in [" << z_bottom << ", " << z_top << "]" << std::endl;
    std::cout << "Number of elements   : " << N.Number_of_elements << std::endl;
    std::cout << "Number of dof        : " << N.Number_of_dof << std::endl;
    std::cout << "Number of unknowns   : " << N.reduced_size() << std::endl;
    vcp::time.toc();

    vcp::time.tic();
    vcp::matrix< TYPE, POLICY > uh = N.torsion_initial_value(initial_peak);
    vcp::time.toc();

    N.setting_newton_tol(newton_tol_factor);

    vcp::time.tic();
    uh = N.solve_nls(uh);
    vcp::time.toc();

    vcp::matrix< TYPE, POLICY > uh_full = N.Dirichlet->expand(uh);
    std::cout << "Convergence          : " << (N.is_convergence() ? "true" : "false") << std::endl;
    std::cout << "max | uh |           : " << max(abs(uh_full))(0) << std::endl;

    // full coefficient vector (large: ~9e5 entries at the default settings)
    // std::cout << uh_full << std::endl;

    // graphics output (GRF-2): (x, y, z, uh) samples + tetrahedron
    // connectivity (0-based; MATLAB tetramesh needs cells + 1).  The full
    // mesh is large (~8e5 point rows at the default h); pass an element list
    // (output_uh_for_graphics(*N.Vh, u, elems)) to sample a subregion.
    {
        vcp::bfem::graphics_output< 3, TYPE, POLICY > g =
            vcp::bfem::output_uh_for_graphics(*N.Vh, N.current_function());
        std::ofstream fp("emden_3dfem_lsc_points.dat");
        fp.precision(17);
        fp << g.points;
        std::ofstream fc("emden_3dfem_lsc_cells.dat");
        fc << g.cells;
        std::cout << "Graphics             : " << g.num_points()
                  << " points / " << g.num_cells()
                  << " cells -> emden_3dfem_lsc_points.dat, emden_3dfem_lsc_cells.dat"
                  << std::endl;
    }

    return 0;
}
