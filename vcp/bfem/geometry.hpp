// vcp/bfem/geometry.hpp
// Layer 2: element geometry (G1) -- element_geometry<D,T>, geometry_traits<T>,
// degenerate_element.
//
// Conforms to: L2 external design v0.3 (section 3) and
//              L2 internal design v0.3 (section 2).
//
// Division-once contract (V6): the single T division of the whole layer is
// step 3 of the construction below (inv_det = 1/det). The sign decision
// (step 2) precedes it, so a degenerate or sign-indefinite det never reaches
// the division (no exploded intervals are ever created).

#ifndef VCP_BFEM_GEOMETRY_HPP
#define VCP_BFEM_GEOMETRY_HPP

#include <array>
#include <string>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/convert_traits.hpp>

namespace vcp {
namespace bfem {

class degenerate_element : public std::runtime_error {
public:
    explicit degenerate_element(const std::string& msg)
        : std::runtime_error(msg) {}
};

// ---------------------------------------------------------------------------
// geometry_traits<T> (V2): decide the sign of det or throw.
// Default: works for point types and also for kv::interval (whose operator<
// is a certain comparison); the kv specialization in geometry_traits_kv.hpp
// makes the interval intent explicit (opt-in).
// ---------------------------------------------------------------------------
template <typename T>
struct geometry_traits {
    static int sign(const T& det) {
        if (det < T(0)) return -1;
        if (T(0) < det) return +1;
        throw degenerate_element(
            "bfem::element_geometry: degenerate or sign-indefinite element");
    }
};

namespace detail {

struct geometry_access;   // element_op plumbing (empty construction, cof)

// det and adjugate, closed forms for D = 2, 3 (internal design 2.1).
// NOTE: D == 3 is implemented per internal design 10.4 (recommended option:
// write it now); its tests are deferred to the 3D enablement phase and it is
// NOT covered by the current gate (documented in the L2 gate report).
template <int D, typename T>
struct geom_kernel;

template <typename T>
struct geom_kernel<2, T> {
    static T det(const std::array<std::array<T, 2>, 2>& B) {
        return B[0][0] * B[1][1] - B[0][1] * B[1][0];
    }
    static void adjugate(const std::array<std::array<T, 2>, 2>& B,
                         std::array<std::array<T, 2>, 2>& adj) {
        adj[0][0] = B[1][1];
        adj[0][1] = -B[0][1];
        adj[1][0] = -B[1][0];
        adj[1][1] = B[0][0];
    }
};

template <typename T>
struct geom_kernel<3, T> {
    static T det(const std::array<std::array<T, 3>, 3>& B) {
        return B[0][0] * (B[1][1] * B[2][2] - B[1][2] * B[2][1])
             - B[0][1] * (B[1][0] * B[2][2] - B[1][2] * B[2][0])
             + B[0][2] * (B[1][0] * B[2][1] - B[1][1] * B[2][0]);
    }
    static void adjugate(const std::array<std::array<T, 3>, 3>& B,
                         std::array<std::array<T, 3>, 3>& adj) {
        adj[0][0] = B[1][1] * B[2][2] - B[1][2] * B[2][1];
        adj[0][1] = -(B[0][1] * B[2][2] - B[0][2] * B[2][1]);
        adj[0][2] = B[0][1] * B[1][2] - B[0][2] * B[1][1];
        adj[1][0] = -(B[1][0] * B[2][2] - B[1][2] * B[2][0]);
        adj[1][1] = B[0][0] * B[2][2] - B[0][2] * B[2][0];
        adj[1][2] = -(B[0][0] * B[1][2] - B[0][2] * B[1][0]);
        adj[2][0] = B[1][0] * B[2][1] - B[1][1] * B[2][0];
        adj[2][1] = -(B[0][0] * B[2][1] - B[0][1] * B[2][0]);
        adj[2][2] = B[0][0] * B[1][1] - B[0][1] * B[1][0];
    }
};

template <int D>
struct factorial_of {  };
template <> struct factorial_of<2> { static const long long value = 2; };
template <> struct factorial_of<3> { static const long long value = 6; };

} // namespace detail

// ---------------------------------------------------------------------------
// element_geometry<D, T> (G1). Value type; coordinates are received in T
// (conversion/enclosure of raw input is the caller's responsibility, V1).
// ---------------------------------------------------------------------------
template <int D, typename T>
class element_geometry {
    static_assert(D == 2 || D == 3,
                  "bfem::element_geometry: only D == 2 or D == 3");
public:
    // Construction order (internal design 2.2; the numbered steps):
    //  0. keep vertices and the edge matrix B (Y6)
    //  1. det (multiplications/additions only)
    //  2. sign via geometry_traits (throws on degeneracy -- BEFORE dividing)
    //  3. inv_det = 1/det                       <- the only division (V6)
    //  4-6. inv_absdet, |T| = |det|/D!
    //  7-8. grad_lambda (adjugate x inv_det), grad_lambda_0 = -sum
    static element_geometry from_vertices(
        const std::array<std::array<T, D>, D + 1>& verts) {
        element_geometry g;
        g.verts_ = verts;                                        // step 0
        for (int c = 0; c < D; ++c)
            for (int r = 0; r < D; ++r)
                g.B_[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] =
                    verts[static_cast<std::size_t>(c + 1)][static_cast<std::size_t>(r)]
                    - verts[0][static_cast<std::size_t>(r)];
        T det = detail::geom_kernel<D, T>::det(g.B_);            // step 1
        g.orient_ = geometry_traits<T>::sign(det);               // step 2
        T inv_det = T(1) / det;                                  // step 3 (division)
        g.inv_absdet_ = (g.orient_ > 0) ? inv_det : -inv_det;    // step 4
        T absdet = (g.orient_ > 0) ? det : -det;                 // step 5
        g.measure_ = absdet * dfact_inv();                       // step 6
        std::array<std::array<T, D>, D> adj;
        detail::geom_kernel<D, T>::adjugate(g.B_, adj);
        for (int d = 0; d < D; ++d) {                            // steps 7-8
            T s0(0);
            T c0(0);
            for (int i = 1; i <= D; ++i) {
                // cof_i = det * grad_lambda_i (kept for the stiffness path,
                // internal design 3); grad_lambda_i = adj row (i-1) * inv_det
                const T& a = adj[static_cast<std::size_t>(i - 1)][static_cast<std::size_t>(d)];
                g.cof_[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] = a;
                g.grad_[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] =
                    a * inv_det;
                s0 -= g.grad_[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
                c0 -= a;
            }
            g.grad_[0][static_cast<std::size_t>(d)] = s0;
            g.cof_[0][static_cast<std::size_t>(d)] = c0;
        }
        return g;
    }

    const T& measure() const { return measure_; }        // |T| (always positive)
    const T& inv_absdet() const { return inv_absdet_; }
    int orientation() const { return orient_; }          // sign of det
    // gradient of lambda_i, component d (cofactor / det; multiplications only)
    const T& grad_lambda(int i, int d) const {
        assert(i >= 0 && i <= D && d >= 0 && d < D);
        return grad_[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
    }

    // ---- v0.3 additions (Y6): receptacle for the future Piola transform ----
    const T& edge_matrix(int r, int c) const {           // B = [v_1-v_0, ..., v_D-v_0]
        assert(r >= 0 && r < D && c >= 0 && c < D);
        return B_[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)];
    }
    const std::array<std::array<T, D>, D + 1>& vertices() const { return verts_; }

private:
    element_geometry()
        : verts_(), B_(), inv_absdet_(), measure_(), grad_(), cof_(), orient_(0) {}

    // 1/D! through convert_traits (enclose-once). Magic static per (D, T):
    // its internal division belongs to the L0 conversion budget, not to the
    // per-element division count (internal design 2.2, step 6 note).
    static const T& dfact_inv() {
        static const T c = rational_to<T>(1, detail::factorial_of<D>::value);
        return c;
    }

    std::array<std::array<T, D>, D + 1> verts_;
    std::array<std::array<T, D>, D> B_;
    T inv_absdet_;
    T measure_;
    std::array<std::array<T, D>, D + 1> grad_;   // grad_lambda[i][d]
    std::array<std::array<T, D>, D + 1> cof_;    // det * grad_lambda[i][d]
    int orient_;

    friend struct detail::geometry_access;
};

namespace detail {
// plumbing for element_op: empty-state construction and the cofactor rows
struct geometry_access {
    template <int D, typename T>
    static element_geometry<D, T> make_empty() { return element_geometry<D, T>(); }
    template <int D, typename T>
    static const T& cof(const element_geometry<D, T>& g, int i, int d) {
        return g.cof_[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
    }
};
} // namespace detail

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_GEOMETRY_HPP
