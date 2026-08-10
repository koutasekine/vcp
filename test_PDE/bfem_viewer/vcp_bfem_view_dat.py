#!/usr/bin/env python3
# test_PDE/vcp_bfem_view_dat.py
#
# GRF-3: convert a GRF-1 (cells, points) .dat pair and display it, in one step.
# A thin wrapper -- it only calls the other two scripts:
#
#     vcp_bfem_dat_to_vtu.py   (python3, needs numpy)   .dat -> .vtu
#     vcp_bfem_view_vtu.py     (pvpython)               .vtu -> screen
#
# Run with python3, NOT pvpython: the converter needs numpy and it is not
# verified that pvpython sees the system numpy.  The viewer is then launched
# as a separate pvpython process, so the numpy question never arises.
#
#     python3 vcp_bfem_view_dat.py cells.dat points.dat
#     python3 vcp_bfem_view_dat.py ns_3dsv_lsc_cells.dat \
#             ns_3dsv_lsc_velocity_points.dat -n velocity
#
# Keep all three scripts in the same folder as the .dat files; then no paths
# are needed anywhere.  The wrapper finds its two siblings next to ITSELF (not
# in the current directory), so it also works when called by a path.
#
# File names are arbitrary and the cells file is never guessed: see the header
# of vcp_bfem_dat_to_vtu.py.
#
# THE PICTURE IS NOT VERIFIED -- see the header of vcp_bfem_dat_to_vtu.py.

import argparse
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONVERT = os.path.join(HERE, "vcp_bfem_dat_to_vtu.py")
VIEW = os.path.join(HERE, "vcp_bfem_view_vtu.py")


def main():
    ap = argparse.ArgumentParser(
        description="Convert a GRF-1 .dat pair and show it in ParaView. "
                    "The rendered picture is NOT verified.")
    ap.add_argument("cells", help="cell connectivity .dat (D+1 integer columns)")
    ap.add_argument("points", help="sample points .dat (D coordinates + values)")
    ap.add_argument("-n", "--name", default=None, help="data array name")
    ap.add_argument("--value", choices=["mid", "lower", "upper"], default=None,
                    help="interval/rational input: which bound to display")
    ap.add_argument("--merge", choices=["auto", "always", "never"], default=None,
                    help="node merging (default auto)")
    ap.add_argument("--glyph", choices=["uniform", "magnitude"], default=None,
                    help="vector arrow length (default uniform)")
    ap.add_argument("--screenshot", default=None,
                    help="write a PNG instead of opening a window")
    ap.add_argument("--keep-vtu", action="store_true",
                    help="say where the .vtu was left (it is kept either way)")
    a = ap.parse_args()

    for p in (CONVERT, VIEW):
        if not os.path.exists(p):
            sys.exit("%s is missing: all three scripts must sit in the same "
                     "folder (expected next to %s)"
                     % (os.path.basename(p), os.path.basename(__file__)))

    cmd = [sys.executable, CONVERT, a.cells, a.points]
    if a.name:
        cmd += ["-n", a.name]
    if a.value:
        cmd += ["--value", a.value]
    if a.merge:
        cmd += ["--merge", a.merge]
    r = subprocess.run(cmd)
    if r.returncode != 0:
        sys.exit(r.returncode)

    stem = a.name or os.path.splitext(os.path.basename(a.points))[0]
    if stem.endswith("_points"):
        stem = stem[:-len("_points")]
    vtu = stem + ".vtu"
    if not os.path.exists(vtu):
        sys.exit("expected %s but it was not written; run the converter alone "
                 "to see why" % vtu)
    if a.keep_vtu:
        print("vtu kept at", os.path.abspath(vtu))

    pv = shutil.which("pvpython") or shutil.which("pvpython.exe")
    if not pv:
        print("\npvpython was not found on PATH, so the viewer cannot be "
              "started.\n%s is ready: open it in ParaView, or install "
              "pvpython\n(on Debian/Ubuntu: the python3-paraview package; "
              "'paraview' alone does not contain pvpython)." % vtu)
        return
    cmd = [pv, VIEW, vtu]
    if a.glyph:
        cmd += ["--glyph", a.glyph]
    if a.screenshot:
        cmd += ["--screenshot", a.screenshot]
    sys.exit(subprocess.run(cmd).returncode)


if __name__ == "__main__":
    main()
