// vcp/bfem/cr1/vcr1_space.hpp
// CONST-G (P1): the VECTOR Crouzeix-Raviart space (CR)^D built on top of the
// scalar cr1_space -- component-major numbering, block-diagonal broken H^1
// stiffness a_h(., .) and mass r(., .), the bundled Dirichlet facet list and
// the mesh size h^2 forwarded from the scalar space.
//
// Conforms to: CONST-G design v1.1 (sections 1, 3, 7) and the P0 ruling of
//              2026-09-07 (file placement: cr1 layer, NO constant here).
//
// ---------------------------------------------------------------------------
// (i) What this header is and is not
// ---------------------------------------------------------------------------
// It is the space V_h^(1) = (V_h^1(Omega))^D of Xie-Xie-Liu 2018, equation
// (7): D copies of the scalar CR space, one per velocity component, with the
// bilinear forms of section 2.2 of that paper,
//
//     a_h(u_h, v_h) = sum_K int_K grad u_h : grad v_h,
//     r(u_h, v_h)   = int_Omega u_h . v_h,
//
// which, in the componentwise l^2 bundling of stokes_constants.hpp (CONST-E
// Lemma 1), are EXACTLY the block-diagonal repetition of the scalar broken
// stiffness and of the scalar mass matrix: a_h(u, v) = sum_d M_h(u_d, v_d)
// and r(u, v) = sum_d N(u_d, v_d).  Nothing is recomputed: the blocks are
// the scalar matrices of cr1_space, scattered D times.
//
// It is NOT the divergence-free subspace V_{0,h}^(1) of equation (9): the
// element-mean divergence constraint, its exact null space and the reduced
// pencil live in vcr1_divfree.hpp.  And, like cr1_space (its ruling R5),
// this header owns NO constant and NO lower bound formula: those are
// consumed in the constants layer (stokes_constants.hpp), which is the only
// place where alpha_1 and eigenvalue_lower_bound are reached.
//
// ---------------------------------------------------------------------------
// (ii) Numbering (normative)
// ---------------------------------------------------------------------------
// Scalar facet dof i in [0, n) of component d in [0, D) carries the global
// vector dof
//
//     g = d * n + i           (component-major; n = scalar ndof).
//
// The Dirichlet list is the scalar boundary list repeated per component and
// is ascending.  Reducing a block-diagonal matrix of this space with
// dirichlet_reduction<T, P, SP>(ndof(), boundary_dofs()) numbers the
// unconstrained dofs ascending, hence the REDUCED vector dof of component d
// and reduced scalar dof r is
//
//     g_red = d * n_r + r      (n_r = the scalar reduced size),
//
// i.e. the reduced numbering is component-major over the scalar reduced
// numbering.  vcr1_divfree.hpp relies on exactly this identity when it
// places the constraint columns; it is stated here once and verified by the
// CONST-G gates.
//
// ---------------------------------------------------------------------------
// (iii) Sources
// ---------------------------------------------------------------------------
//  - V_h^(1) = (V_h^1)^2 and the forms a_h, r: Xie-Xie-Liu, JJIAM 35 (2018)
//    335-354, section 2.2, equations (7) and the display after (8)
//    [verbatim-verified, CONST-G P0];
//  - the block-diagonal identity is the componentwise summation of
//    stokes_constants.hpp Lemma 1 applied to the two forms [self-contained,
//    elementary].
//
// Lexical policy (poisson_constants.hpp design 6.1, inherited): no decimal
// literal, no `double` / `float` token in the code of this file.
// No OpenMP pragma and no hand-written sparse construction: the block
// scatter goes through detail::coo_buffer / detail::spm_adapter like every
// assembly of bfem.

#ifndef VCP_BFEM_CR1_VCR1_SPACE_HPP
#define VCP_BFEM_CR1_VCR1_SPACE_HPP

#include <vector>
#include <cstddef>
#include <stdexcept>

#include <vcp/matrix.hpp>
#include <vcp/spmatrix.hpp>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/fe_space.hpp>        // detail::coo_buffer / detail::spm_adapter
#include <vcp/bfem/cr1/cr1_space.hpp>

namespace vcp {
namespace bfem {

template <int D, typename T, typename P = vcp::mats<T>, class SP = vcp::spmats<T> >
class vcr1_space {
    static_assert(D == 2 || D == 3, "bfem::vcr1_space: only D == 2 or D == 3");
public:
    typedef cr1_space<D, T, P, SP> scalar_space_type;
    typedef typename scalar_space_type::spmatrix_t spmatrix_t;

    explicit vcr1_space(const mesh<D, T>& msh) : scalar_(msh) {}

    // ---- observers ----
    const scalar_space_type& scalar() const { return scalar_; }
    int num_components() const { return D; }
    int scalar_ndof() const { return scalar_.ndof(); }
    int ndof() const { return D * scalar_.ndof(); }
    int num_elements() const { return scalar_.num_elements(); }

    // numbering (ii)
    int vector_dof(int d, int i) const {
        if (d < 0 || d >= D || i < 0 || i >= scalar_.ndof())
            throw std::invalid_argument("bfem::vcr1_space::vector_dof: out of range");
        return d * scalar_.ndof() + i;
    }
    int component_of(int g) const {
        if (g < 0 || g >= ndof())
            throw std::invalid_argument("bfem::vcr1_space::component_of: out of range");
        return g / scalar_.ndof();
    }
    int scalar_dof_of(int g) const {
        if (g < 0 || g >= ndof())
            throw std::invalid_argument("bfem::vcr1_space::scalar_dof_of: out of range");
        return g % scalar_.ndof();
    }

    // the scalar Dirichlet facet list repeated per component; ascending and
    // duplicate free because the scalar list is and the blocks are disjoint
    std::vector<int> boundary_dofs() const {
        const std::vector<int> sb = scalar_.boundary_dofs();
        std::vector<int> out;
        out.reserve(static_cast<std::size_t>(D) * sb.size());
        for (int d = 0; d < D; ++d)
            for (std::size_t k = 0; k < sb.size(); ++k)
                out.push_back(d * scalar_.ndof() + sb[k]);
        return out;
    }

    // ---- bilinear forms: block-diagonal repetition of the scalar ones ----
    spmatrix_t stiffness() const { return block_diagonal_(scalar_.stiffness()); }
    spmatrix_t mass() const { return block_diagonal_(scalar_.mass()); }

    // ---- mesh size: the scalar master API (cr1_space ruling R4) ----
    T max_edge_length_sq() const { return scalar_.max_edge_length_sq(); }

private:
    struct block_scatter {
        detail::coo_buffer<T>* buf;
        int n;
        template <typename I>
        void operator()(I i, I j, const T& v) const {
            for (int d = 0; d < D; ++d)
                buf->push(d * n + static_cast<int>(i), d * n + static_cast<int>(j), v);
        }
    };

    spmatrix_t block_diagonal_(const spmatrix_t& K) const {
        const int n = scalar_.ndof();
        detail::coo_buffer<T> buf;
        block_scatter sc = { &buf, n };
        detail::spm_adapter<T, SP>::for_each_entry(K, sc);
        buf.combine();
        return detail::spm_adapter<T, SP>::build(D * n, D * n, buf);
    }

    scalar_space_type scalar_;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_CR1_VCR1_SPACE_HPP
