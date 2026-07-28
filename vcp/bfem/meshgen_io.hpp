// vcp/bfem/meshgen_io.hpp
// MG-3 (B): exact, machine-independent binary mesh IO (design v1.4, sec. 5).
//
// File format (5.1): little-endian u32/u64 built by arithmetic shifts (no
// memory-image dumps), FNV-1a 64 checksum over everything before the
// checksum field (P2), optional meshgen_status section (P3).
//
// Value records dispatch through meshgen_io_traits<T>; the primary template
// below is the type-generic base-2 path ("b2fp"): sign i8 + exp i64 +
// nbits u32 + mantissa bits in u32 chunks (bit i of the mantissa, i = 0 at
// the leading 1, sits at bit (i mod 32) of chunk (i div 32)). The
// decomposition is the F8 detail of meshgen_convert.hpp; the
// reconstruction is its verbatim inverse walk -- every intermediate equals
// a value that existed in T during decomposition, and the exponent is
// applied ONE halving/doubling at a time (rule 4: never build 2^|e|).
// rational / kv::interval specializations live in meshgen_io_tags.hpp
// (Phase 4), which also carries the type-name tags.
//
// meshgen_status<D> is only forward-declared here: the two-argument
// overloads never touch it, so this header works mesh-alone without any
// MG-1 header; the status-taking overloads are instantiated only from TUs
// that obtained a status (i.e. that include vcp/bfem/meshgen.hpp).
//
// Lexical policy (design 5.4): no decimal literals, no floating point type
// tokens, no sqrt/abs/min/max, integer literals only.

#ifndef VCP_BFEM_MESHGEN_IO_HPP
#define VCP_BFEM_MESHGEN_IO_HPP

#include <array>
#include <cstdint>
#include <fstream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <vcp/bfem/mesh.hpp>
#include <vcp/bfem/meshgen_convert.hpp>

namespace vcp {
namespace bfem {

class meshgen_io_error : public std::runtime_error {
public:
    explicit meshgen_io_error(const char* msg) : std::runtime_error(msg) {}
};

template <int D>
struct meshgen_status;   // full definition: vcp/bfem/meshgen.hpp (MG-1)

namespace meshgen_io_detail {

// ---- FNV-1a 64 running hash -----------------------------------------------
struct fnv1a {
    std::uint64_t h;
    fnv1a() : h(14695981039346656037ull) {}
    void feed(unsigned char b) {
        h = h ^ static_cast<std::uint64_t>(b);
        h = h * 1099511628211ull;
    }
};

// ---- checksummed little-endian writer -------------------------------------
struct owriter {
    std::ostream& os;
    fnv1a sum;
    explicit owriter(std::ostream& s) : os(s), sum() {}

    void byte(unsigned char b) {
        os.put(static_cast<char>(b));
        if (!os)
            throw meshgen_io_error("vcp::bfem::meshgen_io: write failed");
        sum.feed(b);
    }
    void u32(std::uint32_t v) {
        for (int i = 0; i < 4; ++i)
            byte(static_cast<unsigned char>((v >> (8 * i)) & 0xFFu));
    }
    void u64(std::uint64_t v) {
        for (int i = 0; i < 8; ++i)
            byte(static_cast<unsigned char>((v >> (8 * i)) & 0xFFu));
    }
    void i8(int v) { byte(static_cast<unsigned char>(v & 0xFF)); }
    void i32(int v) {
        u32(static_cast<std::uint32_t>(static_cast<std::int64_t>(v) &
                                       0xFFFFFFFFll));
    }
    void i64(long long v) { u64(static_cast<std::uint64_t>(v)); }
    void raw_u64_unsummed(std::uint64_t v) {   // the checksum field itself
        for (int i = 0; i < 8; ++i) {
            os.put(static_cast<char>((v >> (8 * i)) & 0xFFu));
            if (!os)
                throw meshgen_io_error("vcp::bfem::meshgen_io: write failed");
        }
    }
    void text(const std::string& s) {
        for (std::size_t i = 0; i < s.size(); ++i)
            byte(static_cast<unsigned char>(s[i]));
    }
};

// ---- checksummed little-endian reader -------------------------------------
struct ireader {
    std::istream& is;
    fnv1a sum;
    explicit ireader(std::istream& s) : is(s), sum() {}

    unsigned char byte() {
        int c = is.get();
        if (c < 0)
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: truncated or unreadable file");
        unsigned char b = static_cast<unsigned char>(c & 0xFF);
        sum.feed(b);
        return b;
    }
    std::uint32_t u32() {
        std::uint32_t v = 0;
        for (int i = 0; i < 4; ++i)
            v = v | (static_cast<std::uint32_t>(byte()) << (8 * i));
        return v;
    }
    std::uint64_t u64() {
        std::uint64_t v = 0;
        for (int i = 0; i < 8; ++i)
            v = v | (static_cast<std::uint64_t>(byte()) << (8 * i));
        return v;
    }
    int i8() {
        unsigned char b = byte();
        return b < 128u ? static_cast<int>(b) : static_cast<int>(b) - 256;
    }
    int i32() {
        std::uint32_t u = u32();
        std::int64_t v = static_cast<std::int64_t>(u);
        if (u >= 2147483648u) v = v - 4294967296ll;
        return static_cast<int>(v);
    }
    long long i64() {
        std::uint64_t u = u64();
        if ((u & 9223372036854775808ull) != 0ull)
            return -static_cast<long long>(~u + 1ull);
        return static_cast<long long>(u);
    }
    std::uint64_t raw_u64_unsummed() {   // the stored checksum field
        std::uint64_t v = 0;
        for (int i = 0; i < 8; ++i) {
            int c = is.get();
            if (c < 0)
                throw meshgen_io_error(
                    "vcp::bfem::meshgen_io: truncated or unreadable file");
            v = v | (static_cast<std::uint64_t>(
                         static_cast<unsigned char>(c & 0xFF)) << (8 * i));
        }
        return v;
    }
};

} // namespace meshgen_io_detail

// ---------------------------------------------------------------------------
// meshgen_io_traits<T>: value-record serializer. Primary template = the
// type-generic base-2 path "b2fp". Specializations (rational, interval)
// are added by meshgen_io_tags.hpp.
// ---------------------------------------------------------------------------
template <typename T>
struct meshgen_io_traits {
    static std::string tag() { return std::string("b2fp"); }

    static void write_value(meshgen_io_detail::owriter& w, const T& x) {
        // NaN / non-finiteness rejection (design 5.3)
        if (!(x == x))
            throw meshgen_io_error("vcp::bfem::meshgen_io: NaN rejected");
        if (!(x - x == T(0)))
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: non-finite value rejected");
        meshgen_convert_detail::fp_decomp d;
        try {
            d = meshgen_convert_detail::decompose_base2<T>(x);
        } catch (const meshgen_convert_error&) {
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: mantissa decomposition did not "
                "terminate (non-base-2 scalar for the b2fp record)");
        }
        w.i8(d.sign);
        w.i64(d.exp);
        const std::uint32_t nbits = static_cast<std::uint32_t>(d.bits.size());
        w.u32(nbits);
        std::uint32_t chunk = 0;
        for (std::uint32_t i = 0; i < nbits; ++i) {
            if (d.bits[i] != 0)
                chunk = chunk | (static_cast<std::uint32_t>(1) << (i % 32));
            if ((i % 32) == 31 || i + 1 == nbits) {
                w.u32(chunk);
                chunk = 0;
            }
        }
    }

    static T read_value(meshgen_io_detail::ireader& r, bool verify_roundtrip) {
        meshgen_convert_detail::fp_decomp d;
        d.sign = r.i8();
        d.exp = r.i64();
        const std::uint32_t nbits = r.u32();
        if (!(d.sign == 0 || d.sign == 1 || d.sign == -1))
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt value record (sign)");
        if ((d.sign == 0) != (nbits == 0))
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt value record (zero form)");
        d.bits.reserve(nbits);
        std::uint32_t chunk = 0;
        for (std::uint32_t i = 0; i < nbits; ++i) {
            if ((i % 32) == 0) chunk = r.u32();
            d.bits.push_back(
                static_cast<int>((chunk >> (i % 32)) & 1u));
        }
        if (nbits != 0 && (d.bits[0] != 1 || d.bits[nbits - 1] != 1))
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: corrupt value record (mantissa form)");

        if (d.sign == 0) return T(0);

        // verbatim inverse of the decomposition walk (F8): fold the bits
        // from the last one back to the leading one, then apply the
        // exponent one halving/doubling at a time (rule 4). Every
        // intermediate below equals a value that existed in T while the
        // record was being produced (for a matching T).
        const T one = T(1);
        const T two = T(2);
        T m = T(0);
        for (std::uint32_t i = nbits; i > 1; --i)
            m = (m + T(d.bits[i - 1])) / two;
        T y = one + m;
        long long e = d.exp;
        while (e > 0) { y = y * two; e = e - 1; }
        while (e < 0) { y = y / two; e = e + 1; }
        if (d.sign < 0) y = T(0) - y;

        if (verify_roundtrip) {
            // re-decomposition self-check (design 5.2): catches a T whose
            // precision cannot carry the record (e.g. mpfr<N> mismatch)
            meshgen_convert_detail::fp_decomp v;
            try {
                v = meshgen_convert_detail::decompose_base2<T>(y);
            } catch (const meshgen_convert_error&) {
                throw meshgen_io_error(
                    "vcp::bfem::meshgen_io: verify_roundtrip re-decomposition "
                    "failed");
            }
            bool same = (v.sign == d.sign) && (v.exp == d.exp) &&
                        (v.bits.size() == d.bits.size());
            for (std::size_t i = 0; same && i < v.bits.size(); ++i)
                if (v.bits[i] != d.bits[i]) same = false;
            if (!same)
                throw meshgen_io_error(
                    "vcp::bfem::meshgen_io: verify_roundtrip mismatch (value "
                    "not exactly representable in the loading type)");
        }
        return y;
    }
};

// ---------------------------------------------------------------------------
// save_mesh / load_mesh
// ---------------------------------------------------------------------------
namespace meshgen_io_detail {

template <int D, typename T>
void save_mesh_head(owriter& w, const mesh<D, T>& m) {
    w.text(std::string("VCPMESH0"));                       // magic
    w.u32(1u);                                             // format
    w.u32(static_cast<std::uint32_t>(D));
    const std::string tg = meshgen_io_traits<T>::tag();
    w.u32(static_cast<std::uint32_t>(tg.size()));
    w.text(tg);
    w.u64(static_cast<std::uint64_t>(m.num_vertices()));
    w.u64(static_cast<std::uint64_t>(m.num_elements()));
    for (int v = 0; v < m.num_vertices(); ++v)
        for (int d = 0; d < D; ++d)
            meshgen_io_traits<T>::write_value(
                w, m.vertex(v)[static_cast<std::size_t>(d)]);
    for (int e = 0; e < m.num_elements(); ++e)
        for (int k = 0; k <= D; ++k)
            w.i32(m.element(e)[static_cast<std::size_t>(k)]);
}

template <int D, typename T>
mesh<D, T> load_mesh_head(ireader& r, bool verify_roundtrip) {
    static const char* magic = "VCPMESH0";
    for (int i = 0; i < 8; ++i)
        if (r.byte() != static_cast<unsigned char>(magic[i]))
            throw meshgen_io_error("vcp::bfem::meshgen_io: bad magic");
    if (r.u32() != 1u)
        throw meshgen_io_error("vcp::bfem::meshgen_io: unknown format version");
    if (r.u32() != static_cast<std::uint32_t>(D))
        throw meshgen_io_error("vcp::bfem::meshgen_io: dimension mismatch");
    const std::uint32_t taglen = r.u32();
    if (taglen > 4096u)
        throw meshgen_io_error("vcp::bfem::meshgen_io: corrupt type tag");
    std::string tg;
    for (std::uint32_t i = 0; i < taglen; ++i)
        tg.push_back(static_cast<char>(r.byte()));
    if (tg != meshgen_io_traits<T>::tag())
        throw meshgen_io_error("vcp::bfem::meshgen_io: type tag mismatch");
    const std::uint64_t nv = r.u64();
    const std::uint64_t ne = r.u64();
    std::vector<std::array<T, D> > verts;
    for (std::uint64_t v = 0; v < nv; ++v) {
        std::array<T, D> p;
        for (int d = 0; d < D; ++d)
            p[static_cast<std::size_t>(d)] =
                meshgen_io_traits<T>::read_value(r, verify_roundtrip);
        verts.push_back(p);
    }
    std::vector<std::array<int, D + 1> > elems;
    for (std::uint64_t e = 0; e < ne; ++e) {
        std::array<int, D + 1> t;
        for (int k = 0; k <= D; ++k)
            t[static_cast<std::size_t>(k)] = r.i32();
        elems.push_back(t);
    }
    return mesh<D, T>::from_lists(verts, elems);
}

inline void load_skip_status_section(ireader& r, int D) {
    const std::uint64_t nbv = r.u64();
    for (std::uint64_t i = 0; i < nbv; ++i) (void)r.i32();
    const std::uint64_t nbf = r.u64();
    for (std::uint64_t i = 0; i < nbf; ++i)
        for (int d = 0; d < D; ++d) (void)r.i32();
    for (std::uint64_t i = 0; i < nbf; ++i) {
        (void)r.i32();   // facet_source.loop
        (void)r.i32();   // facet_source.segment
    }
    (void)r.i32();       // refine_steps
    (void)r.u32();       // delaunay_complete
}

inline void verify_checksum(ireader& r) {
    const std::uint64_t expect = r.sum.h;
    if (r.raw_u64_unsummed() != expect)
        throw meshgen_io_error("vcp::bfem::meshgen_io: checksum mismatch");
}

} // namespace meshgen_io_detail

// ---- save: without status (works mesh-alone, no MG-1 header needed) -------
template <int D, typename T>
void save_mesh(std::ostream& os, const mesh<D, T>& m) {
    meshgen_io_detail::owriter w(os);
    meshgen_io_detail::save_mesh_head<D, T>(w, m);
    w.u32(0u);                                             // no status section
    w.raw_u64_unsummed(w.sum.h);
}

// ---- save: with optional status (P3). Instantiating this overload needs
// the full meshgen_status<D> from vcp/bfem/meshgen.hpp in the calling TU.
template <int D, typename T>
void save_mesh(std::ostream& os, const mesh<D, T>& m,
               const meshgen_status<D>* status) {
    meshgen_io_detail::owriter w(os);
    meshgen_io_detail::save_mesh_head<D, T>(w, m);
    if (status == 0) {
        w.u32(0u);
    } else {
        w.u32(1u);
        w.u64(static_cast<std::uint64_t>(status->boundary_vertices.size()));
        for (std::size_t i = 0; i < status->boundary_vertices.size(); ++i)
            w.i32(status->boundary_vertices[i]);
        if (status->facet_source.size() != status->boundary_facets.size())
            throw meshgen_io_error(
                "vcp::bfem::meshgen_io: facet_source/boundary_facets size "
                "mismatch in the status to be saved");
        w.u64(static_cast<std::uint64_t>(status->boundary_facets.size()));
        for (std::size_t i = 0; i < status->boundary_facets.size(); ++i)
            for (int d = 0; d < D; ++d)
                w.i32(status->boundary_facets[i][static_cast<std::size_t>(d)]);
        for (std::size_t i = 0; i < status->facet_source.size(); ++i) {
            w.i32(status->facet_source[i].loop);
            w.i32(status->facet_source[i].segment);
        }
        w.i32(status->refine_steps);
        w.u32(status->delaunay_complete ? 1u : 0u);
    }
    w.raw_u64_unsummed(w.sum.h);
}

// ---- load: without status (section, if present, is parsed and dropped) ----
template <int D, typename T>
mesh<D, T> load_mesh(std::istream& is, bool verify_roundtrip = true) {
    meshgen_io_detail::ireader r(is);
    mesh<D, T> m =
        meshgen_io_detail::load_mesh_head<D, T>(r, verify_roundtrip);
    const std::uint32_t flag = r.u32();
    if (flag > 1u)
        throw meshgen_io_error("vcp::bfem::meshgen_io: corrupt status flag");
    if (flag == 1u) meshgen_io_detail::load_skip_status_section(r, D);
    meshgen_io_detail::verify_checksum(r);
    return m;
}

// ---- load: with status out-parameter (untouched when the file carries no
// status section). Needs the full meshgen_status<D> in the calling TU.
template <int D, typename T>
mesh<D, T> load_mesh(std::istream& is, meshgen_status<D>* status,
                     bool verify_roundtrip = true) {
    meshgen_io_detail::ireader r(is);
    mesh<D, T> m =
        meshgen_io_detail::load_mesh_head<D, T>(r, verify_roundtrip);
    const std::uint32_t flag = r.u32();
    if (flag > 1u)
        throw meshgen_io_error("vcp::bfem::meshgen_io: corrupt status flag");
    if (flag == 1u) {
        if (status == 0) {
            meshgen_io_detail::load_skip_status_section(r, D);
        } else {
            status->boundary_vertices.clear();
            status->boundary_facets.clear();
            status->facet_source.clear();
            const std::uint64_t nbv = r.u64();
            for (std::uint64_t i = 0; i < nbv; ++i)
                status->boundary_vertices.push_back(r.i32());
            const std::uint64_t nbf = r.u64();
            for (std::uint64_t i = 0; i < nbf; ++i) {
                std::array<int, D> f;
                for (int d = 0; d < D; ++d)
                    f[static_cast<std::size_t>(d)] = r.i32();
                status->boundary_facets.push_back(f);
            }
            for (std::uint64_t i = 0; i < nbf; ++i) {
                status->facet_source.resize(status->facet_source.size() + 1);
                status->facet_source.back().loop = r.i32();
                status->facet_source.back().segment = r.i32();
            }
            status->refine_steps = r.i32();
            status->delaunay_complete = (r.u32() != 0u);
        }
    }
    meshgen_io_detail::verify_checksum(r);
    return m;
}

// ---- file-path convenience (design section 3) -----------------------------
template <int D, typename T>
void save_mesh(const std::string& path, const mesh<D, T>& m) {
    std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);
    if (!ofs)
        throw meshgen_io_error("vcp::bfem::meshgen_io: cannot open for write");
    save_mesh<D, T>(ofs, m);
}

template <int D, typename T>
void save_mesh(const std::string& path, const mesh<D, T>& m,
               const meshgen_status<D>* status) {
    std::ofstream ofs(path.c_str(), std::ios::out | std::ios::binary);
    if (!ofs)
        throw meshgen_io_error("vcp::bfem::meshgen_io: cannot open for write");
    save_mesh<D, T>(ofs, m, status);
}

template <int D, typename T>
mesh<D, T> load_mesh_file(const std::string& path,
                          bool verify_roundtrip = true) {
    std::ifstream ifs(path.c_str(), std::ios::in | std::ios::binary);
    if (!ifs)
        throw meshgen_io_error("vcp::bfem::meshgen_io: cannot open for read");
    return load_mesh<D, T>(ifs, verify_roundtrip);
}

template <int D, typename T>
mesh<D, T> load_mesh_file(const std::string& path, meshgen_status<D>* status,
                          bool verify_roundtrip = true) {
    std::ifstream ifs(path.c_str(), std::ios::in | std::ios::binary);
    if (!ifs)
        throw meshgen_io_error("vcp::bfem::meshgen_io: cannot open for read");
    return load_mesh<D, T>(ifs, status, verify_roundtrip);
}

} // namespace bfem
} // namespace vcp

#endif // VCP_BFEM_MESHGEN_IO_HPP
