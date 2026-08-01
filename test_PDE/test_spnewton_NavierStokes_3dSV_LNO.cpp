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
// test_spnewton_NavierStokes_3dSV_LNO.cpp
//
// Steady incompressible Navier-Stokes with Scott-Vogelius elements on the
// domain of the numerical example in
//
//   X. Liu, M.T. Nakao, S. Oishi, "Computer-assisted proof for the stationary
//   solution existence of the Navier-Stokes equation over 3D domains",
//   Commun. Nonlinear Sci. Numer. Simul. 108 (2022) 106223, Section 5.3.
//
//     Omega = ( (0,1)^2 \ [0.25, 0.5]^2 ) x (0, 0.5)
//     f     = ( 15 (1-y)^2, 0, 10 z^2 ),   epsilon = 0.25
//     velocity degree k = 3, pressure degree d = k - 1 = 2 (both as in [LNO])
//
// The cross-section is a square with a square hole -- a POLYGON WITH A HOLE
// (polygon_domain::holes), extruded along z.  Its re-entrant features are the
// four edges of the hole, parallel to z, so no Fichera-type mesh (MG-4) is
// needed.
//
// What matches [LNO] and what does not:
//   MATCHES:  domain, viscosity, body force, element degrees (k=3, d=2),
//             Alfeld (barycentric) refinement for the Scott-Vogelius pair.
//   DIFFERS:  mesh generation.  [LNO] subdivides the domain into cubes, cuts
//             each cube into 5 tetrahedra, then applies Zhang's barycentric
//             refinement (600 tetrahedra at h = 0.25).  bfem meshes the
//             polygonal cross-section (ear clipping + red refinement),
//             extrudes prisms, cuts each prism into 3 tetrahedra, and applies
//             alfeld_refine.  Same element family, different base mesh, so
//             dof counts do NOT match [LNO] Table 1.  Also, [LNO] removes the
//             pressure constant mode via the zero-mean space X_{h,0}; here it
//             is removed by pinning one dof.  The approximate solutions of
//             both discretizations converge to the same exact solution as
//             h -> 0, but their finite-h values differ.
//
// Reference values of the approximate solution reported in [LNO] Sec. 5.3
// (on their 600-element mesh):
//     ||u|| = 0.0356,  ||grad u|| = 0.4543,
//     ||u||_inf = 0.1466,  ||grad u||_inf = 11.8763
// This file prints the L2 quantities for comparison.
//
// The body force f = (15(1-y)^2, 0, 10 z^2) depends on the COORDINATES, not
// on a finite element function, so fe_space::load(poly1, fe_function, m)
// cannot be used directly.  Instead the load vector is assembled from
// coordinate functions built by interpolation: for the P^m LAGRANGE basis the
// coefficient vector of any polynomial of degree <= m is simply its values at
// the Lagrange nodes, and fe_space exposes the node coordinates through
// dof_points(m).  The composition 15(1-y)^2 is then evaluated through the
// poly1 machinery on the coordinate function y, which is exact because
// deg( (1-y)^2 ) = 2 <= m.
//
// Build (headers only):
//     g++ -O2 -std=c++11 -I<path to vcp root> test_spnewton_NavierStokes_3dSV_LNO.cpp
// add -fopenmp for the parallel assembly paths of fe_space.
//
// SIZE NOTE.  At mesh_level = 2 (h = 0.25 before Alfeld) the mesh has 6144
// base tetrahedra -- 41x the 600 of [LNO], because the ear-clipping route
// refines the whole cross-section to resolve the hole.  Verified to run in a
// 4 GB single-thread container at mesh_level = 0 only; larger levels need a
// kemeko-class machine.
// ---------------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <array>
#include <memory>
#include <algorithm>
#include <numeric>

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
#include <vcp/bfem/multi_index.hpp>
#include <vcp/bfem/rt/broken_space.hpp>
#include <vcp/bfem/sv/alfeld.hpp>
#include <vcp/bfem/sv/vfe_space.hpp>
#include <vcp/bfem/sv/sv_assemble.hpp>
#include <vcp/bfem/sv/sv_rows.hpp>
#include <vcp/bfem/sv/linear_reduction.hpp>

#include <vcp/vcp_timer.hpp>

template < typename _T, typename _PM, typename _SP >
struct NAVIERSTOKES3DSVLNO : public vcp::SpNewton< _T, _PM, _SP > {

    typedef vcp::bfem::fe_space< 3, _T, _PM, _SP >          scalar_space;
    typedef vcp::bfem::vfe_space< 3, _T, _PM, _SP >         velocity_space;
    typedef vcp::bfem::broken_space< 3, _T, _PM, _SP >      pressure_space;
    typedef vcp::bfem::dirichlet_reduction< _T, _PM, _SP >  dirichlet_type;
    typedef vcp::bfem::linear_reduction< _T, _PM, _SP >     constraint_type;
    typedef vcp::bfem::vfe_function< 3, _T, _PM >           velocity_function;
    typedef vcp::spmatrix< _T, _SP >                        spmatrix_type;
    typedef vcp::matrix< _T, _PM >                          vector_type;

    int m;                                   // velocity degree (k of [LNO])
    _T  nu;                                  // epsilon of [LNO]
    int Number_of_base_elements;
    int Number_of_elements;                  // after Alfeld
    int nU, nP;

    std::unique_ptr< scalar_space >    Vs;
    std::unique_ptr< velocity_space >  Vh;
    std::unique_ptr< pressure_space >  Qh;
    std::unique_ptr< dirichlet_type >  NoSlip;
    std::unique_ptr< constraint_type > SvRows;
    std::unique_ptr< dirichlet_type >  PinPressure;

    std::unique_ptr< vcp::bfem::mesh< 3, _T > > Th;  // kept for coordinate functions
    spmatrix_type Kfull;                     // vector stiffness, UNSCALED, full
    spmatrix_type Mscal;                     // scalar mass, full
    spmatrix_type A;                         // nu * vector stiffness (reduced)
    spmatrix_type B;                         // divergence (reduced)
    vector_type   Fext;                      // body force (reduced)
    vector_type   U, Pr, U_full;

    NAVIERSTOKES3DSVLNO()
        : m(0), nu(_T(1)), Number_of_base_elements(0),
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
        const std::vector< _T > BtP = Pr.vecpointer() * B;
        vector_type N = NoSlip->reduce(
            vcp::bfem::advection_vector(*Vh, u, u, m));

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

        spmatrix_type Bt = transpose(B);          Bt.finalize();
        Bt = _T(-1) * Bt;                         Bt.finalize();
        spmatrix_type Zero(nP, nP);               Zero.finalize();
        spmatrix_type top    = horzcat(Juu, Bt);  top.finalize();
        spmatrix_type bottom = horzcat(B, Zero);  bottom.finalize();
        spmatrix_type J = vercat(top, bottom);    J.finalize();
        return J;
    }

    // ---- set-up ---------------------------------------------------------
    // degree : velocity degree k;  mesh_h : edge-length bound of the
    // triangulation of the cross-section, before the Alfeld refinement
    void first_execute(const int degree, const _T& mesh_h, const _T& viscosity) {
        m  = degree;
        nu = viscosity;

        // 1. cross-section of [LNO]: unit square with the hole [1/4,1/2]^2.
        //    Outer boundary counter-clockwise, hole clockwise; vertices only.
        const _T c0 = _T(0);
        const _T c1 = _T(1);
        const _T q  = _T(1) / _T(4);
        const _T hf = _T(1) / _T(2);

        vcp::bfem::polygon_domain< _T > base;
        base.outer.push_back(pt(c0, c0));
        base.outer.push_back(pt(c1, c0));
        base.outer.push_back(pt(c1, c1));
        base.outer.push_back(pt(c0, c1));

        std::vector< std::array< _T, 2 > > hole;
        hole.push_back(pt(q,  q));
        hole.push_back(pt(q,  hf));
        hole.push_back(pt(hf, hf));
        hole.push_back(pt(hf, q));
        base.holes.push_back(hole);

        vcp::bfem::extruded_domain< _T > Omega;
        Omega.base = base;
        Omega.z0   = c0;
        Omega.z1   = hf;                                  // z in (0, 1/2)

        vcp::bfem::mesh< 3, _T > basemesh = vcp::bfem::generate_mesh(Omega, mesh_h);
        Number_of_base_elements = basemesh.num_elements();

        // 2. Alfeld (barycentric) split for the Scott-Vogelius pair
        Th.reset(new vcp::bfem::mesh< 3, _T >(vcp::bfem::alfeld_refine(basemesh)));
        Number_of_elements = Th->num_elements();

        // 3. velocity P^m (vector), pressure broken P^(m-1)
        Vs.reset(new scalar_space(*Th, m));
        Vh.reset(new velocity_space(*Th, *Vs));
        Qh.reset(new pressure_space(*Th, m - 1));

        NoSlip.reset(new dirichlet_type(Vh->ndof(m), Vh->boundary_dofs(m)));

        vcp::bfem::sv_pressure_constraints< 3, _T > svc(*Th, m);
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

        // 4. matrices independent of u
        Kfull = vcp::bfem::assemble_vector_stiffness(*Vh, m);
        Kfull.finalize();
        Mscal = Vs->mixed_mass(m, m);
        Mscal.finalize();
        spmatrix_type Af = nu * Kfull;
        Af.finalize();
        A  = NoSlip->reduce(Af);
        A.finalize();

        spmatrix_type Bf = vcp::bfem::assemble_div_velocity(*Qh, *Vh, m);
        Bf.finalize();
        spmatrix_type B1 = SvRows->reduce_rows(Bf);       B1.finalize();
        spmatrix_type B2 = PinPressure->reduce_rows(B1);  B2.finalize();
        B  = NoSlip->reduce_cols(B2);
        B.finalize();

        // 5. body force of [LNO]:  f = ( 15(1-y)^2, 0, 10 z^2 )
        Fext = body_force_LNO();
    }

    // ---- norms of the velocity part (for comparison with [LNO]) ---------
    // ||u||^2 = sum_c u_c^T M u_c (scalar mass per component);
    // ||grad u||^2 = U^T K U (unscaled vector stiffness).  Both matrices are
    // already assembled; only mul_vec + inner products are needed.
    _T velocity_l2_norm() {
        using std::sqrt;
        const int N = Vs->ndof(m);
        const std::vector< _T >& uf = U_full.vecpointer();
        _T s = _T(0);
        for (int c = 0; c < 3; c++) {
            std::vector< _T > uc(uf.begin() + c * N, uf.begin() + (c + 1) * N);
            const std::vector< _T > Mu = Mscal.mul_vec(uc);
            s += std::inner_product(Mu.begin(), Mu.end(), uc.begin(), _T(0));
        }
        return sqrt(s);
    }
    _T velocity_h1_seminorm() {
        using std::sqrt;
        const std::vector< _T > Ku = Kfull.mul_vec(U_full.vecpointer());
        return sqrt(std::inner_product(Ku.begin(), Ku.end(),
                                       U_full.vecpointer().begin(), _T(0)));
    }

private:
    static std::array< _T, 2 > pt(const _T& x, const _T& y) {
        std::array< _T, 2 > v;
        v[0] = x;
        v[1] = y;
        return v;
    }

    // Coordinate function x_d as an fe_function.  fe_function stores
    // BERNSTEIN coefficients (fe_function.hpp header note); by the linear
    // precision of the Bernstein basis, the degree-m Bernstein coefficients
    // of a LINEAR polynomial are its values at the barycentric lattice
    // points alpha/m.  So the coefficient of global dof g is the coordinate
    // of its lattice point, reconstructed from the element/local-rank pairs
    // of the dofmap.  Exact for linear functions -- y and z are linear.
    vector_type coordinate_function(const int d) {
        const vcp::bfem::dofmap< 3 >& dm = Vs->dofs(m);
        const vcp::bfem::index_map< 3 > im(m);
        vector_type c;
        c.zeros(Vs->ndof(m), 1);
        for (int e = 0; e < Th->num_elements(); e++) {
            const std::array< int, 4 >& tet = Th->element(e);
            for (int r = 0; r < dm.local_size(); r++) {
                const vcp::bfem::multi_index< 3 > al = im.unrank(r);
                _T x = _T(0);
                for (int i = 0; i < 4; i++) {
                    x += (_T(al.a[i]) / _T(m)) * Th->vertex(tet[i])[d];
                }
                c(dm.global_dof(e, r), 0) = x;
            }
        }
        return c;
    }

    // std::vector -> n x 1 vcp::matrix (std::copy; no hand-written loop)
    static vector_type as_column(const std::vector< _T >& v, int n) {
        vector_type col;
        col.zeros(n, 1);
        std::copy(v.begin(), v.end(), col.data());
        return col;
    }

    vector_type body_force_LNO() {
        // f1 = 15 (1-y)^2 = 15 - 30 y + 15 y^2  (degree 2 <= m: exact)
        std::vector< _T > a1(3, _T(0));
        a1[0] = _T(15);
        a1[1] = _T(-30);
        a1[2] = _T(15);
        vcp::bfem::poly1< _T > f1 = vcp::bfem::poly1< _T >::from_coeffs(a1);

        // f3 = 10 z^2  (degree 2 <= m: exact)
        std::vector< _T > a3(3, _T(0));
        a3[2] = _T(10);
        vcp::bfem::poly1< _T > f3 = vcp::bfem::poly1< _T >::from_coeffs(a3);

        typename scalar_space::function_type ycoord =
            Vs->function_from_coeffs(m, coordinate_function(1));
        typename scalar_space::function_type zcoord =
            Vs->function_from_coeffs(m, coordinate_function(2));

        vector_type l1 = Vs->load(f1, ycoord, m);   // ( 15(1-y)^2, psi_i )
        vector_type l3 = Vs->load(f3, zcoord, m);   // ( 10 z^2,    psi_i )

        // component-major layout: component c occupies [c*N, (c+1)*N)
        const int N = Vs->ndof(m);
        vector_type F;
        F.zeros(Vh->ndof(m), 1);
        std::copy(l1.vecpointer().begin(), l1.vecpointer().end(), F.data());
        std::copy(l3.vecpointer().begin(), l3.vecpointer().end(), F.data() + 2 * N);
        return NoSlip->reduce(F);
    }
};

int main(void) {
    typedef double                      TYPE;
    typedef vcp::mats< TYPE >           POLICY;
    typedef vcp::spmats< TYPE >         SPPOLICY;

    // ---- parameters, all as in [LNO] Sec. 5.3 unless noted ---------------
    const int  velocity_degree = 3;                  // k = 3 (pressure d = 2)
    const TYPE viscosity       = TYPE(1) / TYPE(4);  // epsilon = 0.25
    // Cross-section mesh size BEFORE Alfeld refinement.  [LNO] uses cubes of
    // h = 0.25; bfem's edge-length contract makes 2^-2 the closest analogue,
    // but the resulting element count is larger (see SIZE NOTE above).
    const TYPE mesh_size       = TYPE(1) / TYPE(4);

    // ---- convergence tolerance -------------------------------------------
    // SpNewton's default 4 * eps is a lower bound on what any problem can
    // reach; the attainable value is a property of the problem.
    const int  newton_tol_factor = 512;

    NAVIERSTOKES3DSVLNO< TYPE, POLICY, SPPOLICY > N;

    vcp::time.tic();
    N.first_execute(velocity_degree, mesh_size, viscosity);
    std::cout << "LNO domain ((0,1)^2 \\ [1/4,1/2]^2) x (0,1/2)" << std::endl;
    std::cout << "velocity P^" << velocity_degree
              << ",  pressure broken P^" << (velocity_degree - 1)
              << ",  epsilon = " << viscosity << std::endl;
    std::cout << "Base tetrahedra       : " << N.Number_of_base_elements << std::endl;
    std::cout << "Alfeld tetrahedra     : " << N.Number_of_elements << std::endl;
    std::cout << "Velocity unknowns     : " << N.nU << std::endl;
    std::cout << "Pressure unknowns     : " << N.nP << std::endl;
    vcp::time.toc();

    N.setting_newton_tol(newton_tol_factor);

    // Stokes solution as the initial guess (X = 0: the first Newton step is
    // exactly the Stokes problem since the advection terms vanish at 0)
    vcp::matrix< TYPE, POLICY > x;
    x.zeros(N.nU + N.nP, 1);

    vcp::time.tic();
    x = N.solve_nls(x);
    vcp::time.toc();

    N.setting_newton(x);
    std::cout << "Convergence           : " << (N.is_convergence() ? "true" : "false") << std::endl;
    std::cout << "|| uh ||   (LNO: 0.0356)  : " << N.velocity_l2_norm() << std::endl;
    std::cout << "|| grad uh || (LNO: 0.4543) : " << N.velocity_h1_seminorm() << std::endl;
    std::cout << "max | uh coeff |          : " << max(abs(N.U_full))(0) << std::endl;
    std::cout << "|| div uh ||^2            : "
              << vcp::bfem::div_norm_sq(*N.Vh, N.current_velocity()) << std::endl;

    return 0;
}
