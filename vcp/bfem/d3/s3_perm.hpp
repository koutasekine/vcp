// vcp/bfem/d3/s3_perm.hpp
// Phase 5a (3D common infrastructure): S_3 permutation encoding (normative).
//
// Conforms to: 3D common internal design v0.1 (section 3).
//
// The 6 permutations of {0,1,2} are encoded 0..5 in the lexicographic order
// of the one-line notation (sigma(0), sigma(1), sigma(2)):
//
//   code | one-line  | parity
//   -----+-----------+-------
//     0  | (0,1,2)   |  +1
//     1  | (0,2,1)   |  -1
//     2  | (1,0,2)   |  -1
//     3  | (1,2,0)   |  +1
//     4  | (2,0,1)   |  +1
//     5  | (2,1,0)   |  -1
//
// Everything is a static table lookup; no allocation, no arithmetic beyond
// int indexing. This header is the single authority for the face-orientation
// permutation used by topology3 / dof_build3 / trace_index3.

#ifndef VCP_BFEM_D3_S3_PERM_HPP
#define VCP_BFEM_D3_S3_PERM_HPP

#include <array>
#include <cassert>

namespace vcp {
namespace bfem {
namespace detail {

struct s3_perm {
    // sigma(i) for permutation `code` (one-line notation, table above)
    static int image(int code, int i) {
        assert(code >= 0 && code < 6 && i >= 0 && i < 3);
        static const int tab[6][3] = {
            {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
        };
        return tab[code][i];
    }

    // parity: +1 for even, -1 for odd
    static int parity(int code) {
        assert(code >= 0 && code < 6);
        static const int par[6] = { 1, -1, -1, 1, 1, -1 };
        return par[code];
    }

    // code from the one-line images (s0, s1, s2) = (sigma(0), sigma(1), sigma(2))
    static int from_images(int s0, int s1, int s2) {
        assert(s0 >= 0 && s0 < 3 && s1 >= 0 && s1 < 3 && s2 >= 0 && s2 < 3);
        assert(s0 != s1 && s1 != s2 && s0 != s2);
        // lexicographic rank of the one-line word
        int r = s0 * 2;                 // 2 words share each leading image
        if (s1 > s2) r += 1;            // second image larger => later word
        (void)s2;
        return r;
    }

    // apply to multi-index components: out[sigma(i)] = in[i]
    // (the component seen at local position i moves to canonical position
    // sigma(i); this is the normative direction used by the face convention)
    static void apply(int code, const std::array<int, 3>& in,
                      std::array<int, 3>& out) {
        for (int i = 0; i < 3; ++i)
            out[static_cast<std::size_t>(image(code, i))] =
                in[static_cast<std::size_t>(i)];
    }

    // composition: compose(a, b) is the code of "apply b first, then a",
    // i.e. (a o b)(i) = a(b(i)); parity(compose(a,b)) = parity(a)*parity(b)
    static int compose(int a, int b) {
        assert(a >= 0 && a < 6 && b >= 0 && b < 6);
        return from_images(image(a, image(b, 0)),
                           image(a, image(b, 1)),
                           image(a, image(b, 2)));
    }

    // inverse: image(inverse(c), image(c, i)) == i
    static int inverse(int code) {
        assert(code >= 0 && code < 6);
        static const int inv[6] = { 0, 1, 2, 4, 3, 5 };
        return inv[code];
    }
};

} // namespace detail
} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_D3_S3_PERM_HPP
