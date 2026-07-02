// vcp/bfem/mesh.hpp
// Layer 3: mesh<D,T> (H1) -- thin value type: vertex coordinates in T plus
// the element-node connectivity.
//
// Conforms to: L3 external design v0.2 (section 3.1).
// Checks are index bounds only (X1); orientation is NOT normalized (L2
// accepts both orientations); overlap/coverage is an input contract.

#ifndef VCP_BFEM_MESH_HPP
#define VCP_BFEM_MESH_HPP

#include <vector>
#include <array>
#include <utility>
#include <stdexcept>
#include <cassert>

namespace vcp {
namespace bfem {

template <int D, typename T>
class mesh {
public:
    static mesh from_lists(std::vector<std::array<T, D> > vertices,
                           std::vector<std::array<int, D + 1> > elements) {
        int nv = static_cast<int>(vertices.size());
        for (std::size_t e = 0; e < elements.size(); ++e) {
            for (int k = 0; k <= D; ++k) {
                int v = elements[e][static_cast<std::size_t>(k)];
                if (v < 0 || v >= nv)
                    throw std::invalid_argument("bfem::mesh: element vertex index out of range");
            }
        }
        mesh m;
        m.verts_ = std::move(vertices);
        m.elems_ = std::move(elements);
        return m;
    }

    int num_vertices() const { return static_cast<int>(verts_.size()); }
    int num_elements() const { return static_cast<int>(elems_.size()); }
    const std::array<T, D>& vertex(int v) const {
        assert(v >= 0 && v < num_vertices());
        return verts_[static_cast<std::size_t>(v)];
    }
    const std::array<int, D + 1>& element(int e) const {
        assert(e >= 0 && e < num_elements());
        return elems_[static_cast<std::size_t>(e)];
    }

private:
    mesh() : verts_(), elems_() {}
    std::vector<std::array<T, D> > verts_;
    std::vector<std::array<int, D + 1> > elems_;
};

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESH_HPP
