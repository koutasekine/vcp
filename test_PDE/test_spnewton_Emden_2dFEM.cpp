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
// test_newton3_Emden_2dFEM.cpp
//
// FEM counterpart of test_newton2_Emden.cpp.
//
//   - Delta u = u^2   in Omega,      u = 0 on  d Omega
//
// Omega is the L-shaped domain given ONLY by its vertices below.  The
// spectral (Legendre) discretization of test_newton2_Emden.cpp is replaced by
// a P^p Lagrange finite element space on a triangulation produced by
// vcp::bfem::generate_mesh.  The Newton iteration is vcp::SpNewton::solve_nls
// (sparse Jacobian; see vcp/spnewton.hpp).
//
// Weak form, with V_h = P^p functions vanishing on the boundary:
//     F(u)_i = ( grad u_h, grad psi_i ) - ( u_h^2, psi_i ) = 0
//     F'(u)  = ( grad psi_j, grad psi_i ) - ( 2 u_h psi_j, psi_i )
// so the stiffness matrix plays the role of DL in test_newton2_Emden.cpp,
// fe_space::load supplies (u_h^2, psi_i) and fe_space::weighted_mass supplies
// (2 u_h psi_j, psi_i).  The homogeneous Dirichlet condition is imposed by
// dirichlet_reduction (constrained rows AND columns removed).
//
// Build (this file only depends on headers):
//     g++ -O2 -std=c++11 -I<path to vcp root> test_newton3_Emden_2dFEM.cpp
// add -fopenmp to use the parallel assembly paths of fe_space.
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

#include <vcp/bfem/meshgen.hpp>
#include <vcp/bfem/fe_space.hpp>
#include <vcp/bfem/dirichlet.hpp>
#include <vcp/bfem/poly1.hpp>
#include <vcp/bfem/graphics.hpp>

#include <vcp/vcp_timer.hpp>

template < typename _T, typename _PM, typename _SP >
struct EMDEN2DFEM : public vcp::SpNewton< _T, _PM, _SP > {

    typedef vcp::bfem::fe_space< 2, _T, _PM, _SP >          space_type;
    typedef vcp::bfem::dirichlet_reduction< _T, _PM, _SP >  reduction_type;
    typedef typename space_type::function_type              function_type;

    int p;                                   // FEM degree: the k of P^k
    _T  h;                                   // mesh size
    int Number_of_elements;
    int Number_of_dof;

    std::unique_ptr< space_type >     Vh;    // fe_space (no default ctor)
    std::unique_ptr< reduction_type > Dirichlet;

    vcp::spmatrix< _T, _SP > DL;             // reduced stiffness  ( grad psi_j, grad psi_i )
    vcp::matrix< _T, _PM >   u;              // current iterate, reduced
    vcp::matrix< _T, _PM >   u_full;         // same, expanded with zeros on the boundary

    vcp::bfem::poly1< _T > f_u;              // f(s)  = s^2
    vcp::bfem::poly1< _T > df_u;             // f'(s) = 2 s

    EMDEN2DFEM()
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
    // mesh_h : every edge of the triangulation satisfies len^2 <= mesh_h^2
    void first_execute(const int degree, const _T& mesh_h) {
        p = degree;
        h = mesh_h;

        // 1. L-shaped domain: vertices only (counter-clockwise)
        vcp::bfem::polygon_domain< _T > Omega;
        Omega.outer.push_back(vertex(-1, -1));
        Omega.outer.push_back(vertex( 1, -1));
        Omega.outer.push_back(vertex( 1,  0));
        Omega.outer.push_back(vertex( 0,  0));
        Omega.outer.push_back(vertex( 0,  1));
        Omega.outer.push_back(vertex(-1,  1));

        // 2. triangulation (bfem)
        vcp::bfem::meshgen_status< 2 > status;
        vcp::bfem::mesh< 2, _T > Th = vcp::bfem::generate_mesh(Omega, h, status);

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

    // ---- discretization parameters (both variable) ---------------------
    const int  fem_degree = 3;                       // P^k : k = 3
    const TYPE mesh_size  = TYPE(1) / TYPE(32);      // h   = 2^-5

    // ---- initial value ---------------------------------------------------
    const TYPE initial_value = TYPE(8);

    EMDEN2DFEM< TYPE, POLICY, SPPOLICY > N;

    vcp::time.tic();
    N.first_execute(fem_degree, mesh_size);
    std::cout << "P^" << fem_degree << ",  h = " << mesh_size << std::endl;
    std::cout << "Number of elements   : " << N.Number_of_elements << std::endl;
    std::cout << "Number of dof        : " << N.Number_of_dof << std::endl;
    std::cout << "Number of unknowns   : " << N.reduced_size() << std::endl;
    vcp::time.toc();

    vcp::matrix< TYPE, POLICY > uh;
    uh.zeros(N.reduced_size(), 1);

    // The one hand-written loop in this file: the initial value.  Replace the
    // right-hand side to try a spatially varying initial guess.
    for (int i = 0; i < uh.rowsize(); i++) {
        uh(i, 0) = initial_value;
    }

    // SpNewton's default 4 * eps is a LOWER BOUND on what any problem can
    // reach (f(), Df(), the subtraction and the division each contribute
    // rounding); the tolerance that is actually attainable is a property of
    // the problem and belongs here, not in the library.  With ~3e5 unknowns
    // the correction term stalls at about 5e-14 relative, so 512 is given.
    // Coarser meshes need a smaller number (64 suffices at ~4.5e3 unknowns).
    N.setting_newton_tol(512);

    vcp::time.tic();
    uh = N.solve_nls(uh);
    vcp::time.toc();

    vcp::matrix< TYPE, POLICY > uh_full = N.Dirichlet->expand(uh);
    std::cout << "Convergence          : " << (N.is_convergence() ? "true" : "false") << std::endl;
    std::cout << "max | uh |           : " << max(abs(uh_full))(0) << std::endl;

    // full coefficient vector (large: ~3e5 entries at the default settings)
    // std::cout << uh_full << std::endl;

    // graphics output (GRF-2): (x, y, uh) samples + triangle connectivity,
    // ready for matplotlib Triangulation / MATLAB trisurf (indices are
    // 0-based; MATLAB needs cells + 1).  div = 1 samples the vertices only;
    // raise it to see the P^k shape inside the elements.
    {
        vcp::bfem::graphics_output< 2, TYPE, POLICY > g =
            vcp::bfem::output_uh_for_graphics(*N.Vh, N.current_function());
        std::ofstream fp("emden_2dfem_points.dat");
        fp.precision(17);
        fp << g.points;
        std::ofstream fc("emden_2dfem_cells.dat");
        fc << g.cells;
        std::cout << "Graphics             : " << g.num_points()
                  << " points / " << g.num_cells()
                  << " cells -> emden_2dfem_points.dat, emden_2dfem_cells.dat"
                  << std::endl;
    }

    return 0;
}
