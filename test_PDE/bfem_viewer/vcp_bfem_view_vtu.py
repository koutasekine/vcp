#!/usr/bin/env pvpython
# test_PDE/vcp_bfem_view_vtu.py
#
# GRF-3: display a .vtu written by vcp_bfem_dat_to_vtu.py in ParaView.
#
# MUST be run with pvpython (or pasted into ParaView's Python Shell), because
# it needs paraview.simple:
#
#     pvpython vcp_bfem_view_vtu.py velocity.vtu
#     VCP_VIEW_SCREENSHOT=out.png pvpython vcp_bfem_view_vtu.py velocity.vtu
#
# Keep this script in the same folder as the .vtu and no paths are needed.
# Everything it needs -- array names, component counts, value ranges, the
# domain size -- is read from the file at run time, so there is nothing to
# regenerate when the tooling is updated.  numpy is NOT used: it is not
# verified that pvpython sees the system numpy, and none is needed.
#
# Defaults, and why:
#   scalar  -> Surface coloured by the value (the coloured surface already
#              shows the domain shape)
#   vector  -> Arrows (Glyph) of UNIFORM length, coloured by magnitude, over
#              the boundary feature edges.  Length is driven by the unit
#              direction array "<name>_dir" that the converter writes; scaling
#              length by magnitude leaves most arrows invisible (measured on
#              the SV velocity: 91 % shorter than 2 % of the domain).  The
#              outline is FeatureEdges, not an Outline representation, because
#              the latter draws only the bounding box and hides features such
#              as an L-shaped notch.
#
# THE PICTURE IS NOT VERIFIED.  For interval or rational input the enclosure
# is guaranteed only up to the .vtu; ParaView interpolates inside cells and
# quantises colour.  See the header of vcp_bfem_dat_to_vtu.py.

import argparse
import os
import sys

try:
    from paraview.simple import *          # noqa: F401,F403
except ImportError:
    sys.exit("paraview.simple is not importable: run this with pvpython, "
             "e.g.  pvpython " + os.path.basename(__file__) + " out.vtu\n"
             "(on Debian/Ubuntu pvpython comes from the python3-paraview "
             "package; 'paraview' alone does not contain it)")

# arrow length as a fraction of the domain size, and the fallback factor when
# length is scaled by magnitude instead
UNIFORM_FRACTION = 0.06
MAGNITUDE_FRACTION = 0.12


def pick_arrays(src):
    """Return (name, ncomp, has_dir, has_mag) for the field to display.

    The converter writes "<name>", and for vectors also "<name>_mag" and
    "<name>_dir".  Prefer a genuine vector; otherwise take the first array
    that is not one of the helpers.
    """
    pd = src.PointData
    names = [pd.GetArray(i).Name for i in range(pd.GetNumberOfArrays())]
    helpers = {"_mag", "_dir", "_rad"}
    base = [n for n in names
            if not any(n.endswith(h) for h in helpers)]
    if not base:
        sys.exit("no displayable point array in this file (found: %s)"
                 % ", ".join(names) if names else "no point arrays at all")
    vec = [n for n in base if pd.GetArray(n).GetNumberOfComponents() > 1]
    name = vec[0] if vec else base[0]
    nc = pd.GetArray(name).GetNumberOfComponents()
    return name, nc, (name + "_dir") in names, (name + "_mag") in names


def domain_size(src):
    b = src.GetDataInformation().GetBounds()
    span = max(b[1] - b[0], b[3] - b[2], b[5] - b[4])
    centre = [0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])]
    return (span or 1.0), centre


def show_scalar(src, view, name):
    d = Show(src, view)
    d.SetRepresentationType("Surface")
    ColorBy(d, ("POINTS", name))
    d.RescaleTransferFunctionToDataRange(True, False)
    d.SetScalarBarVisibility(view, True)


def show_vector(src, view, name, has_dir, has_mag, span, uniform):
    # domain outline: the actual boundary edges, not the bounding box
    try:
        ed = Show(FeatureEdges(Input=ExtractSurface(Input=src)), view)
        ed.AmbientColor = ed.DiffuseColor = [0.6, 0.6, 0.6]
    except Exception as e:
        print("FeatureEdges unavailable (%s); falling back to Outline" % e)
        ed = Show(src, view)
        ed.SetRepresentationType("Outline")

    g = Glyph(Input=src, GlyphType="Arrow")
    if uniform and has_dir:
        g.OrientationArray = ["POINTS", name + "_dir"]
        g.ScaleArray = ["POINTS", name + "_dir"]
        g.ScaleFactor = UNIFORM_FRACTION * span
    else:
        if uniform and not has_dir:
            print("no '%s_dir' array in this file: falling back to "
                  "length-by-magnitude, so most arrows will look tiny"
                  % name)
        vmax = max(abs(v) for v in _magnitude_range(src, name)) or 1.0
        g.OrientationArray = ["POINTS", name]
        g.ScaleArray = ["POINTS", name]
        g.ScaleFactor = MAGNITUDE_FRACTION * span / vmax
    g.GlyphMode = "All Points"
    gd = Show(g, view)
    # 2-tuple for a plain scalar; the 3-tuple form with "Magnitude" is only
    # needed when colouring by the vector itself (ColorBy in 5.11.2 maps
    # value[2] == "Magnitude" to component -1)
    if has_mag:
        ColorBy(gd, ("POINTS", name + "_mag"))
    else:
        ColorBy(gd, ("POINTS", name, "Magnitude"))
    gd.RescaleTransferFunctionToDataRange(True, False)
    gd.SetScalarBarVisibility(view, True)


def _magnitude_range(src, name):
    a = src.PointData.GetArray(name)
    nc = a.GetNumberOfComponents()
    # component -1 is the magnitude range in ParaView's array information
    try:
        return a.GetRange(-1)
    except Exception:
        lo = min(a.GetRange(c)[0] for c in range(nc))
        hi = max(a.GetRange(c)[1] for c in range(nc))
        return (lo, hi)


def main():
    ap = argparse.ArgumentParser(
        description="Display a vcp_bfem_dat_to_vtu.py output in ParaView. "
                    "The rendered picture is NOT verified.")
    ap.add_argument("vtu", help="the .vtu to display")
    ap.add_argument("-a", "--array", default=None,
                    help="point array to display (default: the vector array, "
                         "else the first non-helper array)")
    ap.add_argument("--glyph", choices=["uniform", "magnitude"],
                    default="uniform",
                    help="vector arrow length: uniform (default) or "
                         "proportional to magnitude")
    ap.add_argument("--screenshot", default=None,
                    help="write a PNG instead of opening a window "
                         "(same as VCP_VIEW_SCREENSHOT)")
    a = ap.parse_args()

    path = os.path.abspath(a.vtu)
    if not os.path.exists(path):
        sys.exit("%s not found (run it in the folder holding the .vtu)" % a.vtu)

    src = XMLUnstructuredGridReader(FileName=[path])
    src.UpdatePipeline()
    view = GetActiveViewOrCreate("RenderView")

    name, nc, has_dir, has_mag = pick_arrays(src)
    if a.array:
        pd = src.PointData
        if pd.GetArray(a.array) is None:
            sys.exit("no point array named %r in %s" % (a.array, a.vtu))
        name = a.array
        nc = pd.GetArray(name).GetNumberOfComponents()
        has_dir = pd.GetArray(name + "_dir") is not None
        has_mag = pd.GetArray(name + "_mag") is not None
    span, centre = domain_size(src)
    print("%s: array %r, %d component(s), domain size %.6g"
          % (os.path.basename(path), name, nc, span))

    if nc == 1:
        show_scalar(src, view, name)
    else:
        show_vector(src, view, name, has_dir, has_mag, span,
                    a.glyph == "uniform")

    # ResetCamera() alone leaves an axis-aligned view with no sense of depth
    try:
        ResetCameraToDirection(centre, [-1.0, -1.0, -0.6], [0.0, 0.0, 1.0],
                               view)
    except Exception:
        ResetCamera()
    Render()

    shot = a.screenshot or os.environ.get("VCP_VIEW_SCREENSHOT")
    if shot:
        SaveScreenshot(shot, view)
        print("wrote", shot)
        return
    # pvpython exits when the script ends, which would close the window at
    # once; Interact() blocks until it is closed ('q' also quits)
    try:
        Interact(view)
    except Exception as e:
        print("Interact() unavailable (%s); pass --screenshot <path.png> "
              "to save an image instead." % e)


if __name__ == "__main__":
    main()
