// vcp/bfem/rt/rt_element_op.hpp
// RT Layer 2: element operations rt_element_op<D,T,P> (RG1-RG6).
//
// Conforms to: RT-L2 external design v0.2 (sections 1, 2, 3) and
//              RT-L2 internal design v0.2 (sections 2-7).
//
// Piola decomposition (external design 1.1, normative): with
// sigma o F = (1/det) B sigma_hat (signed det, W-RT2), every local quantity
// factors into "geometry factor x RT-L0 reference table":
//   RT mass   : inv_absdet * c_D * sum_{dd'} G_{dd'} R^{(dd')}      (T-R4)
//   div mass  : orient * c_D * DivM                                  (T-R5)
//   cross     : orient * c_D * X'                                    (T-R6')
//   div sigma : coefficients = inv_det_signed * (DivCoef c)          (T-R2)
// c_D = 1/D! (static constant via rational_to), G = B^T B.
//
// D generalization (phase 5c, D5C-5): every formula above is dimension
// uniform (the cross identity (B sigma_hat) . (B^{-T} grad_hat v_hat) =
// sigma_hat . grad_hat v_hat is metric free; the mass carries (1/|det|) G;
// div carries orient c_D). D = 3 only changes the block counts: G = B^T B
// has 6 independent components, c_D = 1/3! = 1/6, and the local mass
// counting norm becomes 6 dim^2 + O(dim) (upper-triangle fold of 6 blocks;
// 2D stays 3 dim^2 + 3, S-RT2-3 unchanged).
//
// Division contract (W-RT1, RC-1): THIS LAYER ADDS NO DIVISION. All inverse
// factors are products of the single L2 geometry division (inv_det). The
// counting test bfem_rt_l2_count_tests.cpp enforces this mechanically.
//
// Sign contract (W-RT2): kernels are sign-free with respect to the global
// edge orientation; only the orient factor of the element appears here.
// RT-L3's signed scatter/gather applies the inter-element signs.
//
// Thread model: one instance is NOT thread safe; clone per thread (the
// shared registries are safe). Buffers are reused; the element loop is
// allocation free from the second element on (RC-7).

#ifndef VCP_BFEM_RT_RT_ELEMENT_OP_HPP
#define VCP_BFEM_RT_RT_ELEMENT_OP_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/matrix.hpp>

#include <vcp/bfem/geometry.hpp>
#include <vcp/bfem/ref_stiffness.hpp>
#include <vcp/bfem/bpoly.hpp>
#include <vcp/bfem/rt/rt_tables.hpp>
#include <vcp/bfem/rt/rt_typed_tables.hpp>

namespace vcp {
namespace bfem {

// -------------------------------------------------------------------------
// rt_local_coeffs<T> (external design 2.1): local RT field coefficients in
// the reference DOF basis, canonical order (RT-L0 section 1.1). No signs
// applied (RT-L3 absorbs them).
// -------------------------------------------------------------------------
template <typename T>
struct rt_local_coeffs {
    int k;                     // RT order
    std::vector<T> c;          // length rt_registry<D>::dim(k)
};

template <int D, typename T, typename P = vcp::mats<T> >
class rt_element_op {
    static_assert(D == 2 || D == 3, "bfem::rt_element_op: only D == 2 or D == 3");
public:
    rt_element_op()
        : geom_(detail::geometry_access::make_empty<D, T>()),
          G_(), f_mass_(), f_div_(), inv_det_signed_(), gp_(),
          dbuf_(), sbuf_(), geometry_set_(false) {}

    // rt_geom_factors (internal design 2): multiplications and sign flips
    // only -- the single division of the element already happened inside
    // element_geometry::from_vertices.
    void set_geometry(const element_geometry<D, T>& g) {
        geom_ = g;
        // G = B^T B (D (D+1) / 2 independent components: 3 for D = 2,
        // 6 for D = 3)
        for (int d = 0; d < D; ++d) {
            for (int dp = d; dp < D; ++dp) {
                T dot = g.edge_matrix(0, d) * g.edge_matrix(0, dp);
                for (int r = 1; r < D; ++r)
                    dot += g.edge_matrix(r, d) * g.edge_matrix(r, dp);
                G_[static_cast<std::size_t>(gidx(d, dp))] = dot;
            }
        }
        f_mass_ = g.inv_absdet() * c_D();
        f_div_ = (g.orientation() > 0) ? c_D() : -c_D();
        inv_det_signed_ = (g.orientation() > 0) ? g.inv_absdet() : -g.inv_absdet();
        // gp_{ij} = |T| grad lambda_i . grad lambda_j (RG5's S_uu factors)
        for (int i = 0; i <= D; ++i) {
            for (int j = i; j <= D; ++j) {
                T dot = g.grad_lambda(i, 0) * g.grad_lambda(j, 0);
                for (int r = 1; r < D; ++r)
                    dot += g.grad_lambda(i, r) * g.grad_lambda(j, r);
                gp_[static_cast<std::size_t>(detail::ref_stiff_block_index(D, i, j))] =
                    g.measure() * dot;
            }
        }
        geometry_set_ = true;
    }

    // ---- RG1: local RT mass, dim(k) x dim(k) ----
    void local_rt_mass(int k, vcp::matrix<T, P>& out) {
        require_geometry();
        check_k(k, "local_rt_mass");
        const typename typed_rt_registry<D, T>::block_table& R =
            typed_rt_registry<D, T>::comp_mass(k);
        const int dim = R.block_rows();
        out.zeros(dim, dim);
        for (int d = 0; d < D; ++d) {
            for (int dp = d; dp < D; ++dp) {
                const T w = f_mass_ * G_[static_cast<std::size_t>(gidx(d, dp))];
                rt_block_view<T> B = R.block(d, dp);
                if (d == dp) {
                    for (int i = 0; i < dim; ++i)
                        for (int j = 0; j < dim; ++j)
                            out(i, j) += w * B.at(i, j);
                } else {
                    // R^{(d'd)} = R^{(dd')^T}: fold the symmetric pair
                    for (int i = 0; i < dim; ++i)
                        for (int j = 0; j < dim; ++j)
                            out(i, j) += w * (B.at(i, j) + B.at(j, i));
                }
            }
        }
    }

    // ---- RG2: local div mass, N(D,l) x dim(k); rows = broken P_l test side
    //      (barycentric definition, W-RT9) ----
    void local_div_mass(int k, int l, vcp::matrix<T, P>& out) {
        require_geometry();
        check_k(k, "local_div_mass");
        if (l < 0)
            throw std::invalid_argument("bfem::rt_element_op::local_div_mass: l < 0");
        const typename typed_rt_registry<D, T>::block_table& DM =
            typed_rt_registry<D, T>::div_mass(k, l);
        out.zeros(DM.rows(), DM.cols());
        for (int i = 0; i < DM.rows(); ++i)
            for (int j = 0; j < DM.cols(); ++j)
                out(i, j) = f_div_ * DM.at(i, j);
    }

    // ---- RG3: local cross (sigma_j, grad psi_alpha), dim(k) x N(D,n);
    //      rows = RT side (W-RT6); computed via T-R6' (W-RT3) ----
    void local_cross_grad(int k, int n, vcp::matrix<T, P>& out) {
        require_geometry();
        check_k(k, "local_cross_grad");
        if (n < 1)
            throw std::invalid_argument("bfem::rt_element_op::local_cross_grad: n < 1");
        const typename typed_rt_registry<D, T>::mat_table& X =
            typed_rt_registry<D, T>::cross_grad_contracted(k, n);
        out.zeros(X.rows(), X.cols());
        for (int j = 0; j < X.rows(); ++j)
            for (int a = 0; a < X.cols(); ++a)
                out(j, a) = f_div_ * X.at(j, a);
    }

    // ---- RG4: physical div sigma_h as a degree k bpoly ----
    // physical div sigma = (1/det) div_hat sigma_hat: coefficients are
    // (DivCoef c_sigma) scaled by the SIGNED inverse determinant.
    void div_as_bpoly(const rt_local_coeffs<T>& sig, bpoly<D, T>& out) {
        require_geometry();
        validate_sig(sig, "div_as_bpoly");
        const typename typed_rt_registry<D, T>::div_table& Dv =
            typed_rt_registry<D, T>::divergence(sig.k);
        detail::bpoly_access::prepare(out, sig.k, false);
        std::vector<T>& o = detail::bpoly_access::vec(out);
        for (int r = 0; r < Dv.rows(); ++r) {
            T acc = Dv.at(r, 0) * sig.c[0];
            for (int j = 1; j < Dv.cols(); ++j)
                acc += Dv.at(r, j) * sig.c[static_cast<std::size_t>(j)];
            o[static_cast<std::size_t>(r)] = acc * inv_det_signed_;
        }
    }

    // ---- RG5: || sigma_h - grad u_h ||^2_{L2(T)} (W-RT4) ----
    // = M_ss - 2 C_su + S_uu; n = 0 (constant u) returns M_ss only (B-1).
    T local_flux_error_sq(const rt_local_coeffs<T>& sig, const bpoly<D, T>& u) {
        require_geometry();
        validate_sig(sig, "local_flux_error_sq");
        const int k = sig.k;
        const int n = u.degree();
        const typename typed_rt_registry<D, T>::block_table& R =
            typed_rt_registry<D, T>::comp_mass(k);
        const int dim = R.block_rows();
        // M_ss = inv_absdet c_D sum_{dd'} G_{dd'} (c^T R^{(dd')} c)
        T Mss(0);
        for (int d = 0; d < D; ++d) {
            for (int dp = d; dp < D; ++dp) {
                rt_block_view<T> B = R.block(d, dp);
                T q(0);
                for (int i = 0; i < dim; ++i) {
                    T row = B.at(i, 0) * sig.c[0];
                    for (int j = 1; j < dim; ++j)
                        row += B.at(i, j) * sig.c[static_cast<std::size_t>(j)];
                    q += sig.c[static_cast<std::size_t>(i)] * row;
                }
                if (d != dp) q += q;             // symmetric pair, no extra mult
                Mss += (f_mass_ * G_[static_cast<std::size_t>(gidx(d, dp))]) * q;
            }
        }
        if (n == 0) return Mss;                  // grad u == 0 (B-1 early path)
        // C_su = orient c_D (c_sigma^T X'(k, n) c_u)
        const typename typed_rt_registry<D, T>::mat_table& X =
            typed_rt_registry<D, T>::cross_grad_contracted(k, n);
        const std::vector<T>& cu = detail::bpoly_access::vec(u);
        T cf(0);
        for (int j = 0; j < X.rows(); ++j) {
            T row = X.at(j, 0) * cu[0];
            for (int a = 1; a < X.cols(); ++a)
                row += X.at(j, a) * cu[static_cast<std::size_t>(a)];
            cf += sig.c[static_cast<std::size_t>(j)] * row;
        }
        const T Csu = f_div_ * cf;
        // S_uu = sum_{i<=j} gp_{ij} (c_u^T R_stiff^{(ij)} c_u)
        const typename detail::ref_stiffness_cache<D, T>::tensor& RS =
            detail::ref_stiffness_cache<D, T>::get(n);
        const int N = RS.N;
        T Suu(0);
        for (int i = 0; i <= D; ++i) {
            for (int j = i; j <= D; ++j) {
                T q(0);
                for (int a = 0; a < N; ++a) {
                    T row = RS.at(i, j, a, 0) * cu[0];
                    for (int b = 1; b < N; ++b)
                        row += RS.at(i, j, a, b) * cu[static_cast<std::size_t>(b)];
                    q += cu[static_cast<std::size_t>(a)] * row;
                }
                if (i != j) q += q;              // symmetric pair
                Suu += gp_[static_cast<std::size_t>(
                    detail::ref_stiff_block_index(D, i, j))] * q;
            }
        }
        return Mss - Csu - Csu + Suu;            // -2C by two subtractions
    }

    // ---- RG6: || div sigma_h + w ||^2_{L2(T)}, w any-degree bpoly ----
    T local_div_residual_sq(const rt_local_coeffs<T>& sig, const bpoly<D, T>& w) {
        require_geometry();
        validate_sig(sig, "local_div_residual_sq");
        div_as_bpoly(sig, dbuf_);
        add_into(sbuf_, dbuf_, w);               // L1 (auto degree elevation)
        return geom_.measure() * ::vcp::bfem::inner(sbuf_, sbuf_);
    }

private:
    element_geometry<D, T> geom_;                // copy (L2 V1 convention)
    std::array<T, static_cast<std::size_t>(D * (D + 1) / 2)> G_;   // upper triangle of B^T B
    T f_mass_;                                   // inv_absdet * c_D
    T f_div_;                                    // orient * c_D
    T inv_det_signed_;                           // orient > 0 ? inv|det| : -inv|det|
    std::array<T, static_cast<std::size_t>((D + 1) * (D + 2) / 2)> gp_;
    bpoly<D, T> dbuf_, sbuf_;                    // RG6 buffers (reused)
    bool geometry_set_;

    // upper-triangle row-major rank of (d, dp), d <= dp:
    // D = 2: (0,0)->0, (0,1)->1, (1,1)->2 (the historical d + dp)
    static int gidx(int d, int dp) {
        return d * D - d * (d - 1) / 2 + (dp - d);
    }
    static const T& c_D() {
        static const T c = rational_to<T>(1, detail::factorial_of<D>::value);
        return c;
    }
    void require_geometry() const {
        if (!geometry_set_)
            throw std::logic_error("bfem::rt_element_op: set_geometry not called");
    }
    static void check_k(int k, const char* where) {
        if (k < 0) {
            std::string msg("bfem::rt_element_op::");
            msg += where;
            msg += ": k < 0";
            throw std::invalid_argument(msg);
        }
    }
    void validate_sig(const rt_local_coeffs<T>& sig, const char* where) const {
        check_k(sig.k, where);
        if (static_cast<int>(sig.c.size()) != rt_registry<D>::dim(sig.k)) {
            std::string msg("bfem::rt_element_op::");
            msg += where;
            msg += ": coefficient size != dim(k)";
            throw std::invalid_argument(msg);
        }
    }
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_RT_RT_ELEMENT_OP_HPP
