// vcp/bfem/meshgen.hpp
// MG-1 Phase A: mesh generation layer for bfem (2D polygonal domains,
// holes allowed, non-convex allowed) with uniform red refinement for the
// mesh size control "every edge satisfies len^2 <= h^2".
//
// Conforms to: MG-1 design v1.2.
//
// Header rules (design section 4.1, checked mechanically by G-MG-1/2):
//  - the only numeric constructions are T(0), T(1), T(2), ... (small ints)
//  - allowed operators on T: binary + - * /, unary -, <, ==, copy/assign
//  - distances are always squared quantities (no square roots)
//  - vertices are identified by combinatorial keys only, never by
//    coordinate comparison (design section 4.3)
//  - every sign decision goes through meshgen_traits<T>::sign
//    (design section 4.4(a)); size decisions use the exception-free
//    conservative rule of section 4.4(b)

#ifndef VCP_BFEM_MESHGEN_HPP
#define VCP_BFEM_MESHGEN_HPP

#include <vector>
#include <array>
#include <map>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <vcp/bfem/mesh.hpp>

namespace vcp {
namespace bfem {

// ---------------------------------------------------------------------------
// exceptions (design section 5.5)
// ---------------------------------------------------------------------------
class meshgen_error : public std::invalid_argument {
public:
    explicit meshgen_error(const char* msg) : std::invalid_argument(msg) {}
};

class meshgen_degenerate : public std::runtime_error {
public:
    explicit meshgen_degenerate(const char* msg) : std::runtime_error(msg) {}
};

class meshgen_limit : public std::runtime_error {
public:
    explicit meshgen_limit(const char* msg) : std::runtime_error(msg) {}
};

// ---------------------------------------------------------------------------
// domain description (design section 3.1)
// ---------------------------------------------------------------------------
template <typename T>
struct polygon_domain {                                    // D = 2
    std::vector<std::array<T, 2> >               outer;    // outer closed polyline
    std::vector<std::vector<std::array<T, 2> > > holes;    // hole polylines
};

// ---------------------------------------------------------------------------
// options (design section 6)
// ---------------------------------------------------------------------------
struct meshgen_options {
    bool validate_input;   // default true (self-intersection check)
    int  max_refine;       // default 30
    meshgen_options() : validate_input(true), max_refine(30) {}
};

// ---------------------------------------------------------------------------
// boundary facet origin (design section 6, C5)
//   loop: 2D: 0 = outer, 1.. = holes[loop-1]
//         3D encoding is defined in meshgen3.hpp
// ---------------------------------------------------------------------------
struct meshgen_facet_source {
    int loop;
    int segment;
    meshgen_facet_source() : loop(0), segment(0) {}
    meshgen_facet_source(int l, int s) : loop(l), segment(s) {}
};

template <int D>
struct meshgen_status {
    std::vector<int> boundary_vertices;                 // ascending, no dups
    std::vector<std::array<int, D> > boundary_facets;   // D=2: edges, D=3: faces
    std::vector<meshgen_facet_source> facet_source;     // same order as facets
    int  refine_steps;
    bool delaunay_complete;                             // true on non-Delaunay path
    meshgen_status()
        : boundary_vertices(), boundary_facets(), facet_source(),
          refine_steps(0), delaunay_complete(true) {}
};

// ---------------------------------------------------------------------------
// meshgen_traits<T>::sign (design section 4.4(a); same regime as the
// geometry_traits of vcp/bfem/geometry.hpp: for kv::interval both < are
// certain comparisons, so an indefinite sign falls through to the throw)
// ---------------------------------------------------------------------------
template <typename T>
struct meshgen_traits {
    static int sign(const T& x) {
        if (x < T(0)) return -1;
        if (T(0) < x) return +1;
        throw meshgen_degenerate(
            "vcp::bfem::meshgen: zero or sign-indefinite value");
    }
};

namespace meshgen_detail {

// orient2d determinant: (b-a) x (c-a); positive iff a,b,c counterclockwise
template <typename T>
T orient2d(const std::array<T, 2>& a,
           const std::array<T, 2>& b,
           const std::array<T, 2>& c) {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

template <typename T>
int orient2d_sign(const std::array<T, 2>& a,
                  const std::array<T, 2>& b,
                  const std::array<T, 2>& c) {
    return meshgen_traits<T>::sign(orient2d(a, b, c));
}

// twice the signed area of a closed loop (shoelace, additions and
// multiplications only; positive iff counterclockwise)
template <typename T>
T twice_signed_area(const std::vector<std::array<T, 2> >& loop) {
    T s(0);
    const int n = static_cast<int>(loop.size());
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1 == n) ? 0 : i + 1;
        s = s + (loop[static_cast<std::size_t>(i)][0] * loop[static_cast<std::size_t>(j)][1]
               - loop[static_cast<std::size_t>(j)][0] * loop[static_cast<std::size_t>(i)][1]);
    }
    return s;
}

// squared distance between two points (the only length measure of the layer)
template <typename T>
T dist2(const std::array<T, 2>& a, const std::array<T, 2>& b) {
    const T dx = b[0] - a[0];
    const T dy = b[1] - a[1];
    return dx * dx + dy * dy;
}

// ---------------------------------------------------------------------------
// Phase 2 (MG-1b): ear clipping + hole bridging.
//
// Adopted degeneracy convention (directive Phase 2 latitude clause; recorded
// in the completion report): decisions that must tolerate an exact-zero or
// uncertain outcome without catching exceptions use one-sided *certain*
// comparisons below; "not certain" always falls to the conservative side
// (candidate ear rejected, vertex treated as blocking). Decisions where a
// silent wrong branch could corrupt the mesh (nearest-crossing selection,
// visible-vertex selection) throw meshgen_degenerate when undecidable.
// Full three-way sign decisions still go through meshgen_traits<T>::sign.
// ---------------------------------------------------------------------------

template <typename T>
bool certainly_less(const T& a, const T& b) { return a < b; }

template <typename T>
bool certainly_equal(const T& a, const T& b) { return a == b; }

template <typename T>
bool certainly_pos(const T& x) { return T(0) < x; }

template <typename T>
bool certainly_neg(const T& x) { return x < T(0); }

// unordered vertex-index pair (int comparison, not a T sign decision)
inline std::pair<int, int> edge_key(int a, int b) {
    return (a < b) ? std::pair<int, int>(a, b) : std::pair<int, int>(b, a);
}

// side of the horizontal +x ray at height my:
//  -1 = certainly below, +1 = certainly above, 0 = exactly on the ray.
// Perturbation rule (directive Phase 2): on-ray vertices count as the upper
// side wherever a two-way classification is needed.
template <typename T>
int ray_side(const T& y, const T& my) {
    if (certainly_less(y, my)) return -1;
    if (certainly_less(my, y)) return +1;
    if (certainly_equal(y, my)) return 0;
    throw meshgen_degenerate("vcp::bfem::meshgen: undecidable ray side");
}

// p certainly outside the closed CCW triangle (a,b,c)?
template <typename T>
bool certainly_outside(const std::array<T, 2>& a, const std::array<T, 2>& b,
                       const std::array<T, 2>& c, const std::array<T, 2>& p) {
    return certainly_neg(orient2d(a, b, p))
        || certainly_neg(orient2d(b, c, p))
        || certainly_neg(orient2d(c, a, p));
}

// coarse mesh: the shared intermediate representation handed to Phase 3
// (refinement) and Phase 4 (Delaunay flips). Triangles are CCW. segment_of
// maps an unordered input-boundary vertex pair to its (loop, segment) origin
// in ORIGINAL input order (loop 0 = outer, 1.. = holes[loop-1]).
template <typename T>
struct coarse_mesh {
    std::vector<std::array<T, 2> > vertices;
    std::vector<std::array<int, 3> > triangles;
    std::map<std::pair<int, int>, meshgen_facet_source> segment_of;
};

// ear clipping of a weakly simple CCW polygon given as a sequence of vertex
// ids (bridge vertices appear twice). O(N^2) scans; conservative blocking.
template <typename T>
void ear_clip(const std::vector<std::array<T, 2> >& V,
              std::vector<int> seq,
              std::vector<std::array<int, 3> >& out) {
    while (static_cast<int>(seq.size()) > 3) {
        const int n = static_cast<int>(seq.size());
        bool clipped = false;
        for (int i = 0; i < n && !clipped; ++i) {
            const int ia = seq[static_cast<std::size_t>((i + n - 1) % n)];
            const int ib = seq[static_cast<std::size_t>(i)];
            const int ic = seq[static_cast<std::size_t>((i + 1) % n)];
            const std::array<T, 2>& A = V[static_cast<std::size_t>(ia)];
            const std::array<T, 2>& B = V[static_cast<std::size_t>(ib)];
            const std::array<T, 2>& C = V[static_cast<std::size_t>(ic)];
            if (!certainly_pos(orient2d(A, B, C))) continue;   // reflex or degenerate
            bool blocked = false;
            for (int j = 0; j < n; ++j) {
                const int id = seq[static_cast<std::size_t>(j)];
                if (id == ia || id == ib || id == ic) continue;
                if (!certainly_outside(A, B, C, V[static_cast<std::size_t>(id)])) {
                    blocked = true;
                    break;
                }
            }
            if (blocked) continue;
            std::array<int, 3> tri = {{ia, ib, ic}};
            out.push_back(tri);
            seq.erase(seq.begin() + i);
            clipped = true;
        }
        if (!clipped)
            throw meshgen_degenerate("vcp::bfem::meshgen: no clippable ear");
    }
    const std::array<T, 2>& A = V[static_cast<std::size_t>(seq[0])];
    const std::array<T, 2>& B = V[static_cast<std::size_t>(seq[1])];
    const std::array<T, 2>& C = V[static_cast<std::size_t>(seq[2])];
    if (!certainly_pos(orient2d(A, B, C)))
        throw meshgen_degenerate("vcp::bfem::meshgen: degenerate final triangle");
    std::array<int, 3> tri = {{seq[0], seq[1], seq[2]}};
    out.push_back(tri);
}

// bridge one hole into the augmented polygon walk `seq` (CCW).
// hw is the hole walk in CW order, rotated so hw[0] is its largest-x vertex.
// Method: +x ray from hw[0], nearest certainly-crossing edge (division-free
// cross-form comparisons), then the Eberly visible-vertex selection with the
// reflex-in-triangle / smallest-angle / smallest-distance fallback.
template <typename T>
void bridge_hole(const std::vector<std::array<T, 2> >& V,
                 std::vector<int>& seq,
                 const std::vector<int>& hw) {
    const int m = hw[0];
    const T& mx = V[static_cast<std::size_t>(m)][0];
    const T& my = V[static_cast<std::size_t>(m)][1];
    const int n = static_cast<int>(seq.size());

    // --- nearest certainly-crossing edge -----------------------------------
    int best_i = -1;              // seq position of the crossed edge start
    int best_up_pos = -1;         // seq position of its upper endpoint
    int best_up_side = 0;         // ray_side of that endpoint (0 = on ray)
    T bestN(0);
    T bestD(1);
    for (int i = 0; i < n; ++i) {
        const int u = seq[static_cast<std::size_t>(i)];
        const int v = seq[static_cast<std::size_t>((i + 1) % n)];
        const int su = ray_side(V[static_cast<std::size_t>(u)][1], my);
        const int sv = ray_side(V[static_cast<std::size_t>(v)][1], my);
        if ((su < 0) == (sv < 0)) continue;   // needs exactly one strictly-below end
        const int low = (su < 0) ? u : v;
        const int up = (su < 0) ? v : u;
        const int up_pos = (su < 0) ? (i + 1) % n : i;
        const int up_side = (su < 0) ? sv : su;
        const T& lx = V[static_cast<std::size_t>(low)][0];
        const T& ly = V[static_cast<std::size_t>(low)][1];
        const T& ux = V[static_cast<std::size_t>(up)][0];
        const T& uy = V[static_cast<std::size_t>(up)][1];
        const T d = uy - ly;                       // certainly positive
        const T N = lx * (uy - my) + ux * (my - ly);   // = xI * d
        const T mxd = mx * d;
        if (certainly_less(N, mxd)) continue;      // crossing behind the ray start
        if (!certainly_less(mxd, N))
            throw meshgen_degenerate(
                "vcp::bfem::meshgen: outer boundary touches a hole vertex");
        bool take = false;
        if (best_i < 0) {
            take = true;
        } else {
            const T lhs = N * bestD;
            const T rhs = bestN * d;
            if (certainly_less(lhs, rhs)) take = true;
            else if (certainly_less(rhs, lhs)) take = false;
            else if (certainly_equal(lhs, rhs)) take = false;   // keep first
            else throw meshgen_degenerate(
                "vcp::bfem::meshgen: undecidable nearest ray crossing");
        }
        if (take) {
            best_i = i;
            best_up_pos = up_pos;
            best_up_side = up_side;
            bestN = N;
            bestD = d;
        }
    }
    if (best_i < 0)
        throw meshgen_error(
            "vcp::bfem::meshgen: hole is not enclosed by the outer boundary");

    // --- visible vertex ----------------------------------------------------
    int w_pos = -1;
    if (best_up_side == 0) {
        // the ray hits the upper endpoint itself: it is the nearest boundary
        // point on the ray, hence visible (perturbation rule)
        w_pos = best_up_pos;
    } else {
        const int eu = seq[static_cast<std::size_t>(best_i)];
        const int ev = seq[static_cast<std::size_t>((best_i + 1) % n)];
        // P = the greater-x endpoint of the crossed edge (ties: the upper one)
        int p_pos;
        const T& eux = V[static_cast<std::size_t>(eu)][0];
        const T& evx = V[static_cast<std::size_t>(ev)][0];
        if (certainly_less(eux, evx)) p_pos = (best_i + 1) % n;
        else if (certainly_less(evx, eux)) p_pos = best_i;
        else if (certainly_equal(eux, evx)) p_pos = best_up_pos;
        else throw meshgen_degenerate(
            "vcp::bfem::meshgen: undecidable crossed-edge endpoint order");
        const int p = seq[static_cast<std::size_t>(p_pos)];
        const int p_side = ray_side(V[static_cast<std::size_t>(p)][1], my);
        // I = ray/edge intersection; bestD is certainly positive
        std::array<T, 2> M = {{mx, my}};
        std::array<T, 2> I = {{bestN / bestD, my}};
        const std::array<T, 2>& P = V[static_cast<std::size_t>(p)];
        // CCW triangle (M,I,P) resp. (M,P,I) depending on the side of P
        const std::array<T, 2>& t0 = M;
        const std::array<T, 2>& t1 = (p_side > 0) ? I : P;
        const std::array<T, 2>& t2 = (p_side > 0) ? P : I;
        // blockers: not-certainly-convex polygon corners inside the triangle
        int best_w = -1;
        T bw_num(0), bw_den(1), bw_d2(0);
        for (int j = 0; j < n; ++j) {
            const int id = seq[static_cast<std::size_t>(j)];
            if (id == eu || id == ev) continue;
            const int jp = seq[static_cast<std::size_t>((j + n - 1) % n)];
            const int jn = seq[static_cast<std::size_t>((j + 1) % n)];
            const std::array<T, 2>& X = V[static_cast<std::size_t>(id)];
            if (certainly_pos(orient2d(V[static_cast<std::size_t>(jp)], X,
                                       V[static_cast<std::size_t>(jn)])))
                continue;                          // strictly convex corner
            if (certainly_outside(t0, t1, t2, X)) continue;
            // candidate: compare by angle to the ray (tan cross-form), then
            // by squared distance to M
            const T ya = (p_side > 0) ? (X[1] - my) : (my - X[1]);
            const T bx = X[0] - mx;
            const T d2 = dist2(M, X);
            bool take = false;
            if (best_w < 0) {
                take = true;
            } else {
                const T lhs = ya * bw_den;
                const T rhs = bw_num * bx;
                if (certainly_less(lhs, rhs)) take = true;
                else if (certainly_less(rhs, lhs)) take = false;
                else if (certainly_equal(lhs, rhs)) {
                    if (certainly_less(d2, bw_d2)) take = true;
                    else if (certainly_less(bw_d2, d2)) take = false;
                    else if (certainly_equal(d2, bw_d2)) take = false;
                    else throw meshgen_degenerate(
                        "vcp::bfem::meshgen: undecidable blocker distance");
                } else throw meshgen_degenerate(
                    "vcp::bfem::meshgen: undecidable blocker angle");
            }
            if (take) {
                best_w = j;
                bw_num = ya;
                bw_den = bx;
                bw_d2 = d2;
            }
        }
        w_pos = (best_w >= 0) ? best_w : p_pos;
    }

    // --- splice: ... W | m h1 .. h_{k-1} m W | next ... ---------------------
    std::vector<int> ins(hw);
    ins.push_back(m);
    ins.push_back(seq[static_cast<std::size_t>(w_pos)]);
    seq.insert(seq.begin() + (w_pos + 1), ins.begin(), ins.end());
}

// build the coarse triangulation of a polygon_domain (design section 5.2,
// steps 1-4): validation, orientation normalization, hole bridging, ear
// clipping. Output triangles are CCW; no vertices are created.
// walk_out receives the bridged polygon walk (bridge edges are the vertex
// pairs traversed twice) for callers that need the constraint-edge set.
template <typename T>
coarse_mesh<T> build_coarse_mesh(const polygon_domain<T>& dom,
                                 const meshgen_options& opt,
                                 std::vector<int>& walk_out) {
    coarse_mesh<T> cm;
    const int nholes = static_cast<int>(dom.holes.size());
    std::vector<std::vector<int> > loop_ids;   // global ids in INPUT order

    for (int l = 0; l <= nholes; ++l) {
        const std::vector<std::array<T, 2> >& loop =
            (l == 0) ? dom.outer : dom.holes[static_cast<std::size_t>(l - 1)];
        const int nl = static_cast<int>(loop.size());
        if (nl < 3)
            throw meshgen_error("vcp::bfem::meshgen: loop with fewer than 3 vertices");
        std::vector<int> ids;
        for (int s = 0; s < nl; ++s) {
            ids.push_back(static_cast<int>(cm.vertices.size()));
            cm.vertices.push_back(loop[static_cast<std::size_t>(s)]);
        }
        // adjacent duplicates, cyclically (certainly equal -> input error)
        for (int s = 0; s < nl; ++s) {
            const std::array<T, 2>& a = loop[static_cast<std::size_t>(s)];
            const std::array<T, 2>& b = loop[static_cast<std::size_t>((s + 1) % nl)];
            if (certainly_equal(a[0], b[0]) && certainly_equal(a[1], b[1]))
                throw meshgen_error(
                    "vcp::bfem::meshgen: adjacent duplicate vertices in a loop");
        }
        loop_ids.push_back(ids);
        // input segment origins
        for (int s = 0; s < nl; ++s) {
            cm.segment_of[edge_key(ids[static_cast<std::size_t>(s)],
                                   ids[static_cast<std::size_t>((s + 1) % nl)])] =
                meshgen_facet_source(l, s);
        }
    }

    // self-intersection check (certain crossings only; design 5.2-1)
    if (opt.validate_input) {
        std::vector<std::pair<int, int> > edges;
        for (std::size_t l = 0; l < loop_ids.size(); ++l) {
            const std::vector<int>& ids = loop_ids[l];
            const int nl = static_cast<int>(ids.size());
            for (int s = 0; s < nl; ++s)
                edges.push_back(std::pair<int, int>(
                    ids[static_cast<std::size_t>(s)],
                    ids[static_cast<std::size_t>((s + 1) % nl)]));
        }
        const int ne = static_cast<int>(edges.size());
        for (int e = 0; e < ne; ++e) {
            for (int f = e + 1; f < ne; ++f) {
                const int a = edges[static_cast<std::size_t>(e)].first;
                const int b = edges[static_cast<std::size_t>(e)].second;
                const int c = edges[static_cast<std::size_t>(f)].first;
                const int d = edges[static_cast<std::size_t>(f)].second;
                if (a == c || a == d || b == c || b == d) continue;
                const std::array<T, 2>& A = cm.vertices[static_cast<std::size_t>(a)];
                const std::array<T, 2>& B = cm.vertices[static_cast<std::size_t>(b)];
                const std::array<T, 2>& C = cm.vertices[static_cast<std::size_t>(c)];
                const std::array<T, 2>& D = cm.vertices[static_cast<std::size_t>(d)];
                const T d1 = orient2d(A, B, C);
                const T d2 = orient2d(A, B, D);
                const T d3 = orient2d(C, D, A);
                const T d4 = orient2d(C, D, B);
                const bool s12 = (certainly_neg(d1) && certainly_pos(d2))
                              || (certainly_pos(d1) && certainly_neg(d2));
                const bool s34 = (certainly_neg(d3) && certainly_pos(d4))
                              || (certainly_pos(d3) && certainly_neg(d4));
                if (s12 && s34)
                    throw meshgen_error(
                        "vcp::bfem::meshgen: self-intersecting input boundary");
            }
        }
    }

    // orientation-normalized walks: outer CCW, holes CW; coordinates of a
    // global id never change, only the walking order
    std::vector<std::vector<int> > walks;
    for (int l = 0; l <= nholes; ++l) {
        const std::vector<int>& ids = loop_ids[static_cast<std::size_t>(l)];
        const int nl = static_cast<int>(ids.size());
        std::vector<std::array<T, 2> > coords;
        for (int s = 0; s < nl; ++s)
            coords.push_back(cm.vertices[static_cast<std::size_t>(
                ids[static_cast<std::size_t>(s)])]);
        const int area_sign = meshgen_traits<T>::sign(twice_signed_area(coords));
        const int want = (l == 0) ? +1 : -1;
        std::vector<int> w(ids);
        if (area_sign != want) {
            for (int s = 0; s < nl; ++s)
                w[static_cast<std::size_t>(s)] =
                    ids[static_cast<std::size_t>(nl - 1 - s)];
        }
        walks.push_back(w);
    }

    // rotate each hole walk to start at its largest-x vertex
    for (int l = 1; l <= nholes; ++l) {
        std::vector<int>& w = walks[static_cast<std::size_t>(l)];
        const int nl = static_cast<int>(w.size());
        int best = 0;
        for (int s = 1; s < nl; ++s) {
            const T& xb = cm.vertices[static_cast<std::size_t>(
                w[static_cast<std::size_t>(best)])][0];
            const T& xs = cm.vertices[static_cast<std::size_t>(
                w[static_cast<std::size_t>(s)])][0];
            if (certainly_less(xb, xs)) best = s;
            else if (certainly_less(xs, xb)) continue;
            else if (certainly_equal(xs, xb)) continue;   // keep first
            else throw meshgen_degenerate(
                "vcp::bfem::meshgen: undecidable hole extreme vertex");
        }
        std::vector<int> r;
        for (int s = 0; s < nl; ++s)
            r.push_back(w[static_cast<std::size_t>((best + s) % nl)]);
        w = r;
    }

    // bridge holes in order of decreasing largest-x (selection order;
    // ties keep input order; undecidable order -> degenerate)
    std::vector<int> order;
    {
        std::vector<bool> used(static_cast<std::size_t>(nholes), false);
        for (int k = 0; k < nholes; ++k) {
            int pick = -1;
            for (int l = 0; l < nholes; ++l) {
                if (used[static_cast<std::size_t>(l)]) continue;
                if (pick < 0) { pick = l; continue; }
                const T& xp = cm.vertices[static_cast<std::size_t>(
                    walks[static_cast<std::size_t>(pick + 1)][0])][0];
                const T& xl = cm.vertices[static_cast<std::size_t>(
                    walks[static_cast<std::size_t>(l + 1)][0])][0];
                if (certainly_less(xp, xl)) pick = l;
                else if (certainly_less(xl, xp)) continue;
                else if (certainly_equal(xl, xp)) continue;   // keep earlier
                else throw meshgen_degenerate(
                    "vcp::bfem::meshgen: undecidable hole ordering");
            }
            used[static_cast<std::size_t>(pick)] = true;
            order.push_back(pick);
        }
    }

    std::vector<int> seq(walks[0]);
    for (std::size_t k = 0; k < order.size(); ++k)
        bridge_hole(cm.vertices, seq,
                    walks[static_cast<std::size_t>(order[k] + 1)]);

    walk_out = seq;
    ear_clip(cm.vertices, seq, cm.triangles);
    return cm;
}

template <typename T>
coarse_mesh<T> build_coarse_mesh(const polygon_domain<T>& dom,
                                 const meshgen_options& opt) {
    std::vector<int> walk;
    return build_coarse_mesh(dom, opt, walk);
}

// ---------------------------------------------------------------------------
// Phase 3 (MG-1c): uniform red refinement + size control + boundary output.
//
// Midpoints are identified by the unordered parent-index pair only (design
// section 4.3); no coordinate comparison occurs. Stopping rule (design
// section 4.4(b)): keep refining unless "len2 < h2" is certainly true for
// every edge; more than opt.max_refine sweeps raises meshgen_limit.
// Boundary lineage (C5): a midpoint of a boundary edge splits it into two
// children that inherit the parent's (loop, segment) origin; no geometric
// decision is involved.
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// MG-2a: shared boundary-lineage primitives. The bodies below are the
// verbatim lineage / boundary-output blocks that used to live inside
// refine_and_finalize; they are shared functions so that the local
// refinement layer (vcp/bfem/meshgen_adapt.hpp) reuses the SAME
// implementation instead of duplicating it (design MG-2 section 4.4).
// No behavioral change (guarded by G-MG2-0 byte-identity).
// ---------------------------------------------------------------------------

// lineage update after a FULL sweep in which every boundary edge was split:
// the two children of each boundary edge inherit the parent's origin. mid
// maps every split edge to its midpoint id; a boundary edge missing from
// mid is an internal error.
inline void inherit_boundary_lineage_after_sweep(
    const std::map<std::pair<int, int>, int>& mid,
    std::map<std::pair<int, int>, meshgen_facet_source>& bnd) {
    std::map<std::pair<int, int>, meshgen_facet_source> nb;
    for (std::map<std::pair<int, int>, meshgen_facet_source>::const_iterator
             it = bnd.begin(); it != bnd.end(); ++it) {
        std::map<std::pair<int, int>, int>::const_iterator mit =
            mid.find(it->first);
        if (mit == mid.end())
            throw meshgen_degenerate(
                "vcp::bfem::meshgen: internal: unsplit boundary edge");
        nb[edge_key(it->first.first, mit->second)] = it->second;
        nb[edge_key(mit->second, it->first.second)] = it->second;
    }
    bnd.swap(nb);
}

// lineage update for ONE split edge (u,v) with midpoint m (the local
// refinement primitive): if (u,v) is a recorded boundary-origin edge it is
// replaced by its two children, which inherit the parent's origin;
// otherwise bnd is untouched. Pure combinatorics, no T arithmetic.
inline void split_boundary_lineage(
    std::map<std::pair<int, int>, meshgen_facet_source>& bnd,
    int u, int v, int m) {
    std::map<std::pair<int, int>, meshgen_facet_source>::iterator it =
        bnd.find(edge_key(u, v));
    if (it == bnd.end()) return;
    const meshgen_facet_source src = it->second;
    bnd.erase(it);
    bnd[edge_key(u, m)] = src;
    bnd[edge_key(m, v)] = src;
}

// boundary output: edges with adjacency count == 1 (same counting as F2),
// each matched with its recorded origin. Fills boundary_facets /
// facet_source / boundary_vertices; refine_steps is the caller's business.
inline void finalize_boundary_output(
    const std::vector<std::array<int, 3> >& tris,
    const std::map<std::pair<int, int>, meshgen_facet_source>& bnd,
    meshgen_status<2>& status) {
    std::map<std::pair<int, int>, int> cnt;
    for (std::size_t t = 0; t < tris.size(); ++t) {
        const std::array<int, 3>& e = tris[t];
        for (int k = 0; k < 3; ++k)
            cnt[edge_key(e[static_cast<std::size_t>(k)],
                         e[static_cast<std::size_t>((k + 1) % 3)])] += 1;
    }
    status.boundary_facets.clear();
    status.facet_source.clear();
    status.boundary_vertices.clear();
    std::map<int, bool> bverts;
    for (std::map<std::pair<int, int>, int>::const_iterator it = cnt.begin();
         it != cnt.end(); ++it) {
        if (it->second != 1) continue;
        std::map<std::pair<int, int>, meshgen_facet_source>::const_iterator
            src = bnd.find(it->first);
        if (src == bnd.end())
            throw meshgen_degenerate(
                "vcp::bfem::meshgen: internal: boundary edge without origin");
        std::array<int, 2> f = {{it->first.first, it->first.second}};
        status.boundary_facets.push_back(f);
        status.facet_source.push_back(src->second);
        bverts[it->first.first] = true;
        bverts[it->first.second] = true;
    }
    for (std::map<int, bool>::const_iterator it = bverts.begin();
         it != bverts.end(); ++it)
        status.boundary_vertices.push_back(it->first);
}

// one uniform red sweep: every edge is split at its midpoint (identified by
// the unordered parent-index pair only), every triangle is replaced by its
// four CCW children, and the boundary lineage is updated through
// inherit_boundary_lineage_after_sweep. Extracted verbatim from
// refine_and_finalize (MG-2a); the T-arithmetic order is unchanged.
template <typename T>
void red_refine_sweep(std::vector<std::array<T, 2> >& vertices,
                      std::vector<std::array<int, 3> >& tris,
                      std::map<std::pair<int, int>, meshgen_facet_source>& bnd) {
    std::map<std::pair<int, int>, int> mid;
    std::vector<std::array<int, 3> > nt;
    for (std::size_t t = 0; t < tris.size(); ++t) {
        const std::array<int, 3>& e = tris[t];
        int m[3];
        for (int k = 0; k < 3; ++k) {
            const int u = e[static_cast<std::size_t>(k)];
            const int v = e[static_cast<std::size_t>((k + 1) % 3)];
            const std::pair<int, int> key = edge_key(u, v);
            std::map<std::pair<int, int>, int>::iterator it = mid.find(key);
            if (it != mid.end()) {
                m[k] = it->second;
            } else {
                std::array<T, 2> c;
                c[0] = (vertices[static_cast<std::size_t>(u)][0]
                      + vertices[static_cast<std::size_t>(v)][0]) / T(2);
                c[1] = (vertices[static_cast<std::size_t>(u)][1]
                      + vertices[static_cast<std::size_t>(v)][1]) / T(2);
                m[k] = static_cast<int>(vertices.size());
                vertices.push_back(c);
                mid[key] = m[k];
            }
        }
        // red split: 3 corner children + the central child, all CCW
        std::array<int, 3> t0 = {{e[0], m[0], m[2]}};
        std::array<int, 3> t1 = {{m[0], e[1], m[1]}};
        std::array<int, 3> t2 = {{m[2], m[1], e[2]}};
        std::array<int, 3> t3 = {{m[0], m[1], m[2]}};
        nt.push_back(t0);
        nt.push_back(t1);
        nt.push_back(t2);
        nt.push_back(t3);
    }
    // boundary lineage: children inherit the parent's origin
    inherit_boundary_lineage_after_sweep(mid, bnd);
    tris.swap(nt);
}

// h2 is the SQUARED edge-length threshold (callers pass h*h; the 3D layer
// passes (h*h)/2 so that prism diagonals also satisfy the h contract)
template <typename T>
void refine_and_finalize(coarse_mesh<T>& cm, const T& h2,
                         const meshgen_options& opt,
                         meshgen_status<2>& status,
                         std::vector<std::array<T, 2> >& vertices,
                         std::vector<std::array<int, 3> >& elements) {
    std::vector<std::array<int, 3> > tris(cm.triangles);
    std::map<std::pair<int, int>, meshgen_facet_source> bnd(cm.segment_of);
    int steps = 0;

    while (true) {
        bool all_small = true;
        for (std::size_t t = 0; t < tris.size() && all_small; ++t) {
            const std::array<int, 3>& e = tris[t];
            for (int k = 0; k < 3 && all_small; ++k) {
                const std::array<T, 2>& a =
                    cm.vertices[static_cast<std::size_t>(e[static_cast<std::size_t>(k)])];
                const std::array<T, 2>& b =
                    cm.vertices[static_cast<std::size_t>(e[static_cast<std::size_t>((k + 1) % 3)])];
                if (!certainly_less(dist2(a, b), h2)) all_small = false;
            }
        }
        if (all_small) break;
        if (steps == opt.max_refine)
            throw meshgen_limit("vcp::bfem::meshgen: max_refine exceeded");

        red_refine_sweep(cm.vertices, tris, bnd);
        ++steps;
    }

    finalize_boundary_output(tris, bnd, status);
    status.refine_steps = steps;

    vertices = cm.vertices;
    elements = tris;
}

} // namespace meshgen_detail

// ---------------------------------------------------------------------------
// public entry points, Phase A (design section 6)
// ---------------------------------------------------------------------------
template <typename T>
void generate_mesh_lists(const polygon_domain<T>& dom, const T& h,
                         std::vector<std::array<T, 2> >& vertices,
                         std::vector<std::array<int, 3> >& elements,
                         meshgen_status<2>& status,
                         const meshgen_options& opt = meshgen_options()) {
    if (!meshgen_detail::certainly_pos(h))
        throw meshgen_error("vcp::bfem::meshgen: h must be certainly positive");
    status = meshgen_status<2>();   // refine_steps 0, delaunay_complete true
    meshgen_detail::coarse_mesh<T> cm =
        meshgen_detail::build_coarse_mesh(dom, opt);
    meshgen_detail::refine_and_finalize(cm, h * h, opt, status,
                                        vertices, elements);
}

template <typename T>
mesh<2, T> generate_mesh(const polygon_domain<T>& dom, const T& h,
                         meshgen_status<2>& status,
                         const meshgen_options& opt = meshgen_options()) {
    std::vector<std::array<T, 2> > vertices;
    std::vector<std::array<int, 3> > elements;
    generate_mesh_lists(dom, h, vertices, elements, status, opt);
    return mesh<2, T>::from_lists(vertices, elements);
}

template <typename T>
mesh<2, T> generate_mesh(const polygon_domain<T>& dom, const T& h,
                         const meshgen_options& opt = meshgen_options()) {
    meshgen_status<2> status;
    return generate_mesh(dom, h, status, opt);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_HPP
