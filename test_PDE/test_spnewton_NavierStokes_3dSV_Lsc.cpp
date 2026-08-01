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
// test_spnewton_NavierStokes_3dSV_Lsc.cpp
//
// Steady incompressible Navier-Stokes on the L-shaped column, discretized
// with Scott-Vogelius elements on an Alfeld (barycentric) split, solved by
// vcp::SpNewton.
//
//   - nu Lap u + (u . grad) u + grad p = f     in Omega
//                             div u    = 0     in Omega
//                                  u   = 0     on d Omega
//
// Omega is the L-shaped COLUMN: the 2D L-shaped cross-section extruded along
// z (re-entrant EDGE, not a corner; the Fichera corner is the MG-4 track).
//
// Discretization:
//   mesh      generate_mesh(extruded_domain, h)  then  alfeld_refine
//   velocity  vfe_space, vector P^m, homogeneous Dirichlet
//   pressure  broken_space, discontinuous P^(m-1)
//   pressure constraints  sv_pressure_constraints -> linear_reduction,
//                         plus one pinned dof (the constant null space)
//
// Block system with X = (U, P):
//   F(X) = [ nu A U + N(U) - B^t P - F_ext ]      A   = vector stiffness
//          [ B U                           ]      B   = divergence
//                                                 N(U)= (u.grad u, v)
//   F'(X) = [ nu A + C(U) + D(U)   -B^t ]         C(w)= (w.grad u, v)
//           [ B                      0  ]         D(w)= (u.grad w, v)
//
// Block assembly uses the MATLAB-style horzcat / vercat / transpose free
// functions of vcp::spmatrix and vcp::matrix, so this file contains no
// hand-written loop at all.
//
// The Scott-Vogelius property to look for in the output is ||div u||^2: it
// should be at the level of rounding, not at the level of a discretization
// error.
//
// Build (headers only):
//     g++ -O2 -std=c++11 -I<path to vcp root> test_spnewton_NavierStokes_3dSV_Lsc.cpp
// add -fopenmp for the parallel assembly paths of fe_space.
//
// SIZE WARNING.  The saddle-point Jacobian is far heavier to factorize than
// the scalar problems of the Emden samples: at the delivered defaults the
// system is large enough that only a machine of kemeko class can hold it.
// Lower velocity_degree / raise mesh_size to run it elsewhere.
//
// SV STABILITY.  Following the library's policy, no precondition on the pair
// (mesh, degree) is checked here; sv_pressure_constraints reports a stability
// HINT only, which this file prints and does not act on.  If an unstable pair
// is chosen the saddle-point matrix can be singular and the sparse LU will
// report a failure -- that is the expected behaviour, not a defect.
// ---------------------------------------------------------------------------

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
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/sv/alfeld.hpp>
#include <vcp/bfem/sv/vfe_space.hpp>
#include <vcp/bfem/sv/sv_assemble.hpp>
#include <vcp/bfem/sv/sv_rows.hpp>
#include <vcp/bfem/sv/linear_reduction.hpp>

#include <vcp/vcp_timer.hpp>

template < typename _T, typename _PM, typename _SP >
struct NAVIERSTOKES3DSV : public vcp::SpNewton< _T, _PM, _SP > {

    typedef vcp::bfem::fe_space< 3, _T, _PM, _SP >          scalar_space;
    typedef vcp::bfem::vfe_space< 3, _T, _PM, _SP >         velocity_space;
    typedef vcp::bfem::broken_space< 3, _T, _PM, _SP >      pressure_space;
    typedef vcp::bfem::dirichlet_reduction< _T, _PM, _SP >  dirichlet_type;
    typedef vcp::bfem::linear_reduction< _T, _PM, _SP >     constraint_type;
    typedef vcp::bfem::vfe_function< 3, _T, _PM >           velocity_function;
    typedef vcp::spmatrix< _T, _SP >                        spmatrix_type;
    typedef vcp::matrix< _T, _PM >                          vector_type;

    int m;                                   // velocity degree: the k of P^k
    _T  h;                                   // mesh size (before Alfeld)
    _T  nu;                                  // viscosity
    int Number_of_base_elements;
    int Number_of_elements;                  // after Alfeld
    int nU, nP;                              // reduced velocity / pressure sizes

    std::unique_ptr< scalar_space >    Vs;   // scalar space (wrapped)
    std::unique_ptr< velocity_space >  Vh;
    std::unique_ptr< pressure_space >  Qh;
    std::unique_ptr< dirichlet_type >  NoSlip;      // velocity Dirichlet
    std::unique_ptr< constraint_type > SvRows;      // SV pressure constraints
    std::unique_ptr< dirichlet_type >  PinPressure; // remove the constant mode

    spmatrix_type A;                         // nu * vector stiffness (reduced)
    spmatrix_type B;                         // divergence (reduced)
    vector_type   Fext;                      // body force (reduced)
    vector_type   U, Pr, U_full;             // current iterate

    NAVIERSTOKES3DSV()
        : m(0), h(_T(0)), nu(_T(1)), Number_of_base_elements(0),
          Number_of_elements(0), nU(0), nP(0) {}

    // ---- vcp::SpNewton hooks -------------------------------------------
    void setting_newton(vector_type& x) override {
        U.zeros(nU, 1);
        Pr.zeros(nP, 1);
        std::copy(x.vecpointer().begin(), x.vecpointer().begin() + nU, U.data());
        std::copy(x.vecpointer().begin() + nU, x.vecpointer().end(), Pr.data());
        U_full = NoSlip->expand(U);
    }

    velocity_function current_velocity() {
        return Vh->function_from_coeffs(m, U_full);
    }

    vector_type f() override {
        velocity_function u = current_velocity();
        const std::vector< _T > AU  = A.mul_vec(U.vecpointer());
        const std::vector< _T > BU  = B.mul_vec(U.vecpointer());
        const std::vector< _T > BtP = Pr.vecpointer() * B;      // B^t P
        vector_type N = NoSlip->reduce(
            vcp::bfem::advection_vector(*Vh, u, u, m));

        // velocity rows: nu A U + N(U) - B^t P - F_ext ; pressure rows: B U
        return vercat(as_column(AU, nU) + N - as_column(BtP, nU) - Fext,
                      as_column(BU, nP));
    }

    spmatrix_type Df() override {
        velocity_function u = current_velocity();
        spmatrix_type C = vcp::bfem::assemble_advection(*Vh, u, m);
        C.finalize();
        spmatrix_type Dw = vcp::bfem::assemble_advection_derivative(*Vh, u, m);
        Dw.finalize();
        spmatrix_type CD = C + Dw;
        CD.finalize();
        spmatrix_type Juu = A + NoSlip->reduce(CD);
        Juu.finalize();

        // [ J_uu  -B^t ]
        // [ B       0   ]     (MATLAB-style horzcat / vercat)
        spmatrix_type Bt = transpose(B);          Bt.finalize();
        Bt = _T(-1) * Bt;                         Bt.finalize();
        spmatrix_type Zero(nP, nP);               Zero.finalize();
        spmatrix_type top = horzcat(Juu, Bt);     top.finalize();
        spmatrix_type bottom = horzcat(B, Zero);  bottom.finalize();
        spmatrix_type J = vercat(top, bottom);    J.finalize();
        return J;
    }

    // ---- set-up ---------------------------------------------------------
    void first_execute(const int degree, const _T& mesh_h, const _T& viscosity,
                       const _T& z0, const _T& z1) {
        m  = degree;
        h  = mesh_h;
        nu = viscosity;

        // 1. L-shaped column: cross-section by vertices only, extruded in z
        vcp::bfem::extruded_domain< _T > Omega;
        Omega.base.outer.push_back(vertex(-1, -1));
        Omega.base.outer.push_back(vertex( 1, -1));
        Omega.base.outer.push_back(vertex( 1,  0));
        Omega.base.outer.push_back(vertex( 0,  0));
        Omega.base.outer.push_back(vertex( 0,  1));
        Omega.base.outer.push_back(vertex(-1,  1));
        Omega.z0 = z0;
        Omega.z1 = z1;

        vcp::bfem::mesh< 3, _T > base = vcp::bfem::generate_mesh(Omega, h);
        Number_of_base_elements = base.num_elements();

        // 2. Alfeld (barycentric) split: each tetrahedron into 4 children
        vcp::bfem::mesh< 3, _T > Th = vcp::bfem::alfeld_refine(base);
        Number_of_elements = Th.num_elements();

        // 3. velocity P^m (vector) and pressure broken P^(m-1)
        Vs.reset(new scalar_space(Th, m));
        Vh.reset(new velocity_space(Th, *Vs));
        Qh.reset(new pressure_space(Th, m - 1));

        // 4. velocity: homogeneous Dirichlet on the whole boundary
        NoSlip.reset(new dirichlet_type(Vh->ndof(m), Vh->boundary_dofs(m)));

        // 5. pressure: Scott-Vogelius constraint rows, then pin one dof to
        //    remove the constant null space (enclosed flow)
        vcp::bfem::sv_pressure_constraints< 3, _T > svc(Th, m);
        vcp::bfem::sv_stability_hint hint =
            svc.hint(vcp::bfem::sv_provenance_alfeld);
        std::cout << "SV singular edges     : " << svc.num_singular_edges() << std::endl;
        std::cout << "SV constraint rows    : " << (int)svc.rows().size() << std::endl;
        std::cout << "SV stability hint     : " << hint.reason << std::endl;

        SvRows.reset(new constraint_type(Qh->ndof(), svc.rows()));
        std::vector< int > pinned(1, 0);
        PinPressure.reset(new dirichlet_type(SvRows->reduced_size(), pinned));

        nU = NoSlip->reduced_size();
        nP = PinPressure->reduced_size();

        // 6. matrices (independent of u, assembled once)
        spmatrix_type Af = vcp::bfem::assemble_vector_stiffness(*Vh, m);
        Af.finalize();
        Af = nu * Af;
        A  = NoSlip->reduce(Af);
        A.finalize();

        spmatrix_type Bf = vcp::bfem::assemble_div_velocity(*Qh, *Vh, m);
        Bf.finalize();
        spmatrix_type B1 = SvRows->reduce_rows(Bf);       B1.finalize();
        spmatrix_type B2 = PinPressure->reduce_rows(B1);  B2.finalize();
        B  = NoSlip->reduce_cols(B2);
        B.finalize();

        // 7. body force f = (w, 0, 0) with  -Lap w = 1,  w = 0 on the boundary.
        //    A CONSTANT force would be a gradient and would be absorbed
        //    entirely by the pressure, giving u = 0 exactly -- useless as a
        //    test.  The torsion function is not a gradient field.
        Fext = body_force_from_torsion();
    }

private:
    // std::vector -> n x 1 vcp::matrix (std::copy; no hand-written loop)
    static vector_type as_column(const std::vector< _T >& v, int n) {
        vector_type c;
        c.zeros(n, 1);
        std::copy(v.begin(), v.end(), c.data());
        return c;
    }

    static std::array< _T, 2 > vertex(const int x, const int y) {
        std::array< _T, 2 > v;
        v[0] = _T(x);
        v[1] = _T(y);
        return v;
    }

    vector_type body_force_from_torsion() {
        std::vector< _T > one(1, _T(1));
        vcp::bfem::poly1< _T > constant_one =
            vcp::bfem::poly1< _T >::from_coeffs(one);

        dirichlet_type sd(Vs->ndof(m), Vs->dofs(m).boundary_dofs());
        spmatrix_type K = sd.reduce(Vs->stiffness(m));
        K.finalize();
        vector_type Fw =
            sd.reduce(Vs->load(constant_one, Vs->zero_function(m), m));

        vcp::linear_solve_options< _T > options;
        options.method = vcp::linear_solver_method::sparse_lu;
        const std::vector< _T > wv = K.solve(Fw.vecpointer(), options);

        vector_type wr;
        wr.zeros(sd.reduced_size(), 1);
        std::copy(wv.begin(), wv.end(), wr.data());
        vector_type w_full = sd.expand(wr);

        std::vector< _T > id(2, _T(0));
        id[1] = _T(1);                                   // f(s) = s
        vcp::bfem::poly1< _T > identity =
            vcp::bfem::poly1< _T >::from_coeffs(id);
        vector_type sload =
            Vs->load(identity, Vs->function_from_coeffs(m, w_full), m);

        // component-major: component 0 occupies [0, N)
        vector_type F;
        F.zeros(Vh->ndof(m), 1);
        std::copy(sload.vecpointer().begin(), sload.vecpointer().end(), F.data());
        return NoSlip->reduce(F);
    }
};

int main(void) {
    typedef double                      TYPE;
    typedef vcp::mats< TYPE >           POLICY;
    typedef vcp::spmats< TYPE >         SPPOLICY;

    // ---- discretization parameters (all variable) -----------------------
    const int  velocity_degree = 3;                  // P^k velocity, k = 3
                                                     // (3D Scott-Vogelius on an
                                                     //  Alfeld split wants k >= 3)
    const TYPE mesh_size       = TYPE(1) / TYPE(8);  // h = 2^-3, before Alfeld
    const TYPE viscosity       = TYPE(1);            // nu
    const TYPE z_bottom        = TYPE(0);
    const TYPE z_top           = TYPE(1);

    // ---- convergence tolerance -------------------------------------------
    // SpNewton's default 4 * eps is a lower bound on what any problem can
    // reach; the attainable value is a property of the problem.  Raise it if
    // Newton stalls just above the tolerance.
    const int  newton_tol_factor = 1024;

    NAVIERSTOKES3DSV< TYPE, POLICY, SPPOLICY > N;

    vcp::time.tic();
    N.first_execute(velocity_degree, mesh_size, viscosity, z_bottom, z_top);
    std::cout << "L-shaped column,  velocity P^" << velocity_degree
              << ",  pressure broken P^" << (velocity_degree - 1)
              << ",  h = " << mesh_size << ",  nu = " << viscosity << std::endl;
    std::cout << "Base tetrahedra       : " << N.Number_of_base_elements << std::endl;
    std::cout << "Alfeld tetrahedra     : " << N.Number_of_elements << std::endl;
    std::cout << "Velocity unknowns     : " << N.nU << std::endl;
    std::cout << "Pressure unknowns     : " << N.nP << std::endl;
    vcp::time.toc();

    N.setting_newton_tol(newton_tol_factor);

    // Stokes solution as the initial guess: X = 0 makes the first Newton step
    // solve exactly the Stokes problem, since N(0) = C(0) = D(0) = 0.
    vcp::matrix< TYPE, POLICY > x;
    x.zeros(N.nU + N.nP, 1);

    vcp::time.tic();
    x = N.solve_nls(x);
    vcp::time.toc();

    N.setting_newton(x);
    std::cout << "Convergence           : " << (N.is_convergence() ? "true" : "false") << std::endl;
    std::cout << "max | uh |            : " << max(abs(N.U_full))(0) << std::endl;
    std::cout << "|| div uh ||^2        : "
              << vcp::bfem::div_norm_sq(*N.Vh, N.current_velocity()) << std::endl;

    return 0;
}
