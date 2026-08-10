#!/usr/bin/env python3
# test_PDE/vcp_bfem_dat_to_vtu.py
#
# GRF-3: convert the (points, cells) .dat pair written by
# vcp::bfem::output_uh_for_graphics (GRF-1) into a ParaView-readable .vtu,
# plus a small ParaView Python script that sets up a sensible default view.
#
# Dependencies: numpy only (zlib / base64 / struct / re / fractions are in the
# standard library).  ParaView itself is needed only to VIEW the result.
#
# ============================================================================
# THE PICTURE IS NOT VERIFIED.
#
# For interval or rational input this script is a convenience for LOOKING AT
# verified output; it is not a verified visualisation.  What is guaranteed:
#   * the .dat values themselves enclose the exact values (GRF-1 contract),
#   * parsing is lossless: interval bounds are printed with 17 significant
#     digits and round-trip through double exactly; rational values are
#     converted with OUTWARD rounding so that lower <= exact <= upper,
#   * the .vtu contains exactly those numbers (verifiable by read-back).
# What is NOT guaranteed, and cannot be:
#   * picking one bound (or the midpoint) destroys the enclosure statement --
#     the pair (midpoint coordinate, lower value) encloses nothing,
#   * ParaView interpolates linearly inside each cell; interpolated values are
#     not enclosures of the finite element function,
#   * colour mapping quantises (typically 256 levels) and the screen shows no
#     interval width at all.
# For interval input the radius is therefore ALSO written, as a separate array
# (<name>_rad).  Looking at where the radius grows is the informative view;
# the value picture is not evidence of anything.
# ============================================================================
#
# Accepted input: ONLY the GRF-1 pair, i.e. a points file with (D + ncomp)
# columns and a cells file with (D + 1) columns of 0-based row indices.
#   points 3 cols + cells 3 cols -> 2D scalar (P^k, C^1)
#   points 4 cols + cells 4 cols -> 3D scalar (P^k, broken pressure)
#   points 6 cols + cells 4 cols -> 3D vector (SV velocity, ncomp = 3)
# NOT accepted: the ldbase (spectral) output of
# Legendre_Bases_Generator::output_uh_for_graphics(Div) -- that is a single
# matrix on a tensor-product grid with no cell list -- nor data.log,
# Data_*/ files, or vcp_fio binary saves.
#
# The file names are arbitrary: nothing is inferred from them, and the cells
# file is never guessed.  Guessing it would silently pick up another mesh's
# connectivity, which stays inside the index range and produces a meaningless
# picture without any error.
#
# Usage (keep this script in the same folder as the .dat files; then no paths
# are needed anywhere)
#   python3 vcp_bfem_dat_to_vtu.py cells.dat points.dat
#   python3 vcp_bfem_dat_to_vtu.py ns_3dsv_lsc_cells.dat \
#           ns_3dsv_lsc_velocity_points.dat -n velocity
#   python3 vcp_bfem_dat_to_vtu.py cells.dat p_points.dat --value lower
#
# Then:  pvpython vcp_bfem_view_vtu.py <name>.vtu
# or in one step:  python3 vcp_bfem_view_dat.py cells.dat points.dat

import argparse
import base64
import os
import re
import struct
import sys
import zlib
from fractions import Fraction

import numpy as np

VTK_TRIANGLE, VTK_TETRA = 5, 10
_BLOCK = 1 << 15
_IV = re.compile(r"\[\s*([^,\]]+?)\s*,\s*([^\]]+?)\s*\]")


# ---------------------------------------------------------------- reading ----

def _sniff(path):
    """Return 'double' | 'interval' | 'rational' from the first data line."""
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if "[" in line:
                return "interval"
            if "/" in line:
                return "rational"
            return "double"
    sys.exit("%s: empty file" % path)


def _read_double(path):
    a = np.loadtxt(path, ndmin=2)
    return a, None


def _read_interval(path):
    lo, hi = [], []
    with open(path) as fh:
        for line in fh:
            m = _IV.findall(line)
            if not m:
                continue
            lo.append([float(a) for a, _ in m])
            hi.append([float(b) for _, b in m])
    if not lo:
        sys.exit("%s: looked like interval data but no [a,b] found" % path)
    return np.array(lo), np.array(hi)


def _read_rational(path):
    """Exact parse, then OUTWARD rounding to double."""
    rows = []
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if t:
                rows.append([Fraction(x) for x in t])
    if not rows:
        sys.exit("%s: no rational rows" % path)
    lo = np.empty((len(rows), len(rows[0])))
    hi = np.empty_like(lo)
    for i, r in enumerate(rows):
        if len(r) != lo.shape[1]:
            sys.exit("%s: ragged rows" % path)
        for j, v in enumerate(r):
            f = float(v)                      # nearest
            fv = Fraction(f)
            lo[i, j] = np.nextafter(f, -np.inf) if fv > v else f
            hi[i, j] = np.nextafter(f, np.inf) if fv < v else f
    return lo, hi


def read_points(path):
    kind = _sniff(path)
    lo, hi = {"double": _read_double,
              "interval": _read_interval,
              "rational": _read_rational}[kind](path)
    return kind, lo, hi


def read_cells(path):
    """Cell table: non-negative integers, D+1 columns.  Anything else is
    almost certainly the points file passed in the wrong position -- reject,
    never silently reorder (integer-valued point data would make an automatic
    decision wrong, and a wrong guess produces a plausible but meaningless
    picture)."""
    try:
        raw = np.loadtxt(path, ndmin=2)
    except ValueError:
        sys.exit("%s: not a numeric table.  Argument order is "
                 "<cells> <points>." % path)
    if raw.shape[1] not in (3, 4):
        sys.exit("%s: cells must have 3 or 4 columns, found %d.  Argument "
                 "order is <cells> <points>." % (path, raw.shape[1]))
    if not np.all(raw == np.floor(raw)) or raw.min() < 0:
        sys.exit("%s: cells must be non-negative integers.  Argument order "
                 "is <cells> <points>." % path)
    return raw.astype(np.int64)


def degenerate_check(X, C, D):
    """Every cell must have non-zero volume.  A cells table from another mesh
    typically produces flat or inverted cells here.  Same caveat as
    consistency_check: this catches the common mistake, not every one."""
    P = np.asarray(X, float)[C]
    E = P[:, 1:, :] - P[:, :1, :]
    det = np.linalg.det(E) if D == 3 else (E[:, 0, 0] * E[:, 1, 1]
                                          - E[:, 0, 1] * E[:, 1, 0])
    bad = int(np.count_nonzero(det == 0.0))
    if bad:
        sys.exit("%d of %d cells have zero volume: the cells file does not "
                 "match these points." % (bad, len(C)))


def consistency_check(C, npoint, D):
    """npoint = nsel * C(div+D, D) and ncell = nsel * div^D must hold for one
    integer pair (nsel, div).  A cells file from a DIFFERENT mesh usually
    fails this.  It is not a proof: a coincidence in both counts would pass,
    so this rejects the common mistake, not every mistake."""
    from math import comb
    ncell = len(C)
    for div in range(1, 33):
        npt = comb(div + D, D)
        cpe = div ** D
        if npoint % npt or ncell % cpe:
            continue
        if npoint // npt == ncell // cpe:
            return npoint // npt, div
    sys.exit("point rows (%d) and cell rows (%d) are not consistent with any "
             "(element count, div) for D = %d: the two files probably come "
             "from different meshes." % (npoint, ncell, D))


# ------------------------------------------------------------- vtu writer ----

def _enc(a):
    raw = np.ascontiguousarray(a).tobytes()
    blocks = [raw[i:i + _BLOCK] for i in range(0, len(raw), _BLOCK)] or [b""]
    comp = [zlib.compress(b) for b in blocks]
    head = struct.pack("<3Q", len(blocks), _BLOCK, len(blocks[-1]))
    head += struct.pack("<%dQ" % len(comp), *(len(c) for c in comp))
    return base64.b64encode(head).decode() + base64.b64encode(b"".join(comp)).decode()


def write_vtu(path, X, cells, cell_type, arrays):
    n, m, k = len(X), len(cells), cells.shape[1]
    X = np.asarray(X, np.float64)
    if X.shape[1] == 2:
        X = np.column_stack([X, np.zeros(n)])
    body = ['<?xml version="1.0"?>',
            '<VTKFile type="UnstructuredGrid" version="1.0" '
            'byte_order="LittleEndian" header_type="UInt64" '
            'compressor="vtkZLibDataCompressor">',
            '<UnstructuredGrid>',
            '<Piece NumberOfPoints="%d" NumberOfCells="%d">' % (n, m),
            '<Points>',
            '<DataArray type="Float64" NumberOfComponents="3" format="binary">',
            _enc(X), '</DataArray>', '</Points>', '<Cells>']
    for name, arr, ty in (("connectivity", np.asarray(cells, np.int64).ravel(), "Int64"),
                          ("offsets", np.arange(k, k * m + 1, k, dtype=np.int64), "Int64"),
                          ("types", np.full(m, cell_type, np.uint8), "UInt8")):
        body += ['<DataArray type="%s" Name="%s" format="binary">' % (ty, name),
                 _enc(arr), '</DataArray>']
    body.append('</Cells>')
    if arrays:
        body.append('<PointData Scalars="%s">' % next(iter(arrays)))
        for name, v in arrays.items():
            v = np.asarray(v, np.float64)
            nc = 1 if v.ndim == 1 else v.shape[1]
            # NumberOfComponents is omitted for scalars: VTK defaults to 1 and
            # readers then hand back a 1-D array instead of (n, 1).
            comp = '' if nc == 1 else ' NumberOfComponents="%d"' % nc
            body += ['<DataArray type="Float64" Name="%s"%s format="binary">'
                     % (name, comp),
                     _enc(v), '</DataArray>']
        body.append('</PointData>')
    body += ['</Piece>', '</UnstructuredGrid>', '</VTKFile>']
    with open(path, "w") as fh:
        fh.write("\n".join(body))


# ---------------------------------------------------------------- merging ----
#
# GRF-1 emits shared points once per incident element on purpose.  Merging is
# safe for a CONTINUOUS field and destroys the information for a
# DISCONTINUOUS one (the SV pressure jump was measured at 77 % of the field
# magnitude).  Rather than trusting file names, decide by measurement: group
# exactly coincident coordinates and look at the spread of the values.

def merge_decision(X, V, rtol):
    uniq, first, inv = np.unique(X, axis=0, return_index=True, return_inverse=True)
    if len(uniq) == len(X):
        return False, 0.0, uniq, first, inv           # nothing to merge
    order = np.argsort(inv, kind="stable")
    g, Vs = inv[order], V[order]
    cut = np.flatnonzero(np.diff(g)) + 1
    spread = 0.0
    for seg in np.split(Vs, cut):
        if len(seg) > 1:
            spread = max(spread, float(np.abs(seg - seg[0]).max()))
    scale = float(np.abs(V).max()) or 1.0
    return spread <= rtol * scale, spread / scale, uniq, first, inv


# -------------------------------------------------------------------- main ---

def main():
    ap = argparse.ArgumentParser(
        description="GRF-1 (points, cells) .dat -> .vtu for ParaView. "
                    "The rendered picture is NOT verified; see the script header.")
    # cells first: the topology is fixed while the field varies (the SV
    # velocity and pressure share one cells file), and it mirrors "geometry
    # before values" -- the same order as inside a points row.
    ap.add_argument("cells", help="cell connectivity .dat (D+1 integer columns)")
    ap.add_argument("points", help="sample points .dat (D coordinates + ncomp values)")
    ap.add_argument("-n", "--name", default=None,
                    help="data array name (default: from the points file name)")
    ap.add_argument("-o", "--out", default=None, help="output .vtu path")
    ap.add_argument("--value", choices=["mid", "lower", "upper"], default="mid",
                    help="interval/rational: which bound becomes the value "
                         "array (default mid). Coordinates always use the "
                         "midpoint; picking a bound voids any enclosure claim.")
    ap.add_argument("--merge", choices=["auto", "always", "never"], default="auto",
                    help="node merging (default auto: merge only if coincident "
                         "points carry equal values, i.e. a continuous field)")
    ap.add_argument("--merge-rtol", type=float, default=0.0,
                    help="relative tolerance for the auto decision (default 0: "
                         "merge only on exact agreement)")
    a = ap.parse_args()

    C = read_cells(a.cells)
    kind, lo, hi = read_points(a.points)
    D = C.shape[1] - 1
    if D not in (2, 3):
        sys.exit("cells must have 3 or 4 columns (got %d)" % C.shape[1])
    ncomp = lo.shape[1] - D
    if ncomp < 1:
        sys.exit("points has %d columns but cells implies D = %d: this does not "
                 "look like GRF-1 output (ldbase graphics output is not "
                 "supported)" % (lo.shape[1], D))
    if C.max() >= len(lo):
        sys.exit("cell indices reach %d but there are only %d point rows: "
                 "0-based GRF-1 output expected, and the two files must come "
                 "from the same run" % (C.max(), len(lo)))
    nsel, div = consistency_check(C, len(lo), D)

    name = a.name or re.sub(r"_points$", "", os.path.splitext(
        os.path.basename(a.points))[0]) or "u"
    out = a.out or (name + ".vtu")

    if kind == "double":
        X, V, rad = lo[:, :D], lo[:, D:], None
    else:
        X = 0.5 * (lo[:, :D] + hi[:, :D])
        pick = {"mid": 0.5 * (lo + hi), "lower": lo, "upper": hi}[a.value]
        V = pick[:, D:]
        rad = 0.5 * (hi[:, D:] - lo[:, D:])
    print("input: %s, D = %d, ncomp = %d, %d point rows, %d cells "
          "(%d elements, div = %d)"
          % (kind, D, ncomp, len(X), len(C), nsel, div))
    if kind != "double":
        print("       value array = %s; radius also written as %s_rad"
              % (a.value, name))
        print("       NOTE: the rendered picture is not verified "
              "(see the script header)")

    degenerate_check(X, C, D)
    ok, rel, uniq, first, inv = merge_decision(X, V, a.merge_rtol)
    do = {"always": True, "never": False, "auto": ok}[a.merge]
    if len(uniq) != len(X):
        print("       coincident points: %d -> %d, relative value spread %.3g "
              "-> merge %s" % (len(X), len(uniq), rel, "yes" if do else "no"))
        if not do and a.merge == "auto":
            print("       (kept unmerged: the field is discontinuous there)")
    if do:
        X, V, C = uniq, V[first], inv[C]
        if rad is not None:
            rad = rad[first]

    arrays = {name: (V[:, 0] if ncomp == 1 else V)}
    if ncomp > 1:
        # magnitude as a plain scalar (robust ColorBy target) and the unit
        # direction, so the glyph can be uniform in length
        mag = np.linalg.norm(V, axis=1)
        nz = mag > 0.0
        d = np.zeros_like(V)
        d[nz] = V[nz] / mag[nz, None]
        arrays[name + "_mag"] = mag
        arrays[name + "_dir"] = d
    if rad is not None:
        arrays[name + "_rad"] = (rad[:, 0] if ncomp == 1 else rad)
    write_vtu(out, X, C, VTK_TRIANGLE if D == 2 else VTK_TETRA, arrays)
    print("wrote %s (%.2f MB)" % (out, os.path.getsize(out) / 1e6))
    print("view it with:  pvpython vcp_bfem_view_vtu.py %s" % out)


if __name__ == "__main__":
    main()
