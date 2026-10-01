"""Render the pictures a run's [postprocess] sections asked for.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" postprocess.py [output-dir]

Reads `<output-dir>/postprocess.json`, which the solver writes, and produces one
PNG per request beside it. The solver renders nothing and this script parses no
input file: the manifest is the whole interface, and it carries the resolved
body tags, the metre-converted geometry and the array names so that nothing here
has to know how the solver spells things.

WHAT THE DISPLAYS MEAN. A vector phasor Ehat = P + jQ has
E(t) = P cos(wt) - Q sin(wt), whose tip traces an ELLIPSE with semi-axes a >= b.
Three different numbers are all reasonably called "the magnitude", and this is
the whole reason the manifest names one explicitly:

    complex_magnitude   N = sqrt(|P|^2+|Q|^2) = sqrt(a^2+b^2) = sqrt(2)*RMS
    magnitude_at_phase  |P cos(th) - Q sin(th)|            depends on th
    peak                a                                  what it reaches
    axial_ratio         b/a in [0,1]; 0 linear, 1 circular

with b <= magnitude_at_phase <= peak <= complex_magnitude <= sqrt(2)*peak.
docs/ComplexVectorPhasorConcept.md derives all of it. Every figure this script
writes is ANNOTATED with which one it is and at what phase, because the
difference is invisible in the picture and a screenshot outlives its caption.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import json
import math
import os
import sys

import numpy as np   # imported AFTER paraview.simple on purpose: that star
                     # import shadows several builtins, and np.* is unambiguous

paraview.simple._DisableFirstRenderCameraReset()

OUT = sys.argv[1] if len(sys.argv) > 1 else "output"
MANIFEST = os.path.join(OUT, "postprocess.json")

# This project's colour-map names -> the ParaView preset that implements each.
#
# The indirection is not decoration: preset names are ParaView's and they move
# between versions. "Viridis (matplotlib)" was the spelling in one release and
# is simply absent in 6.1, where it is "Viridis" -- which this script found out
# by failing on it. Keeping our own names means a version bump is one edit in
# this table rather than a hunt through every case.
PRESET = {
    "rainbow": "Rainbow Uniform",
    "jet": "Jet",
    "turbo": "Turbo",
    "cool_to_warm": "Cool to Warm",
    "viridis": "Viridis",
    "blue_to_red": "Blue to Red Rainbow",
    "black_body": "Black-Body Radiation",
    "grayscale": "Grayscale",
    "x_ray": "X Ray",
}

# Labels for the annotation stamped on every frame.
WHAT = {
    "complex_magnitude": "complex magnitude  sqrt(|P|^2+|Q|^2) = sqrt(a^2+b^2) = sqrt(2) x RMS",
    "magnitude_at_phase": "instantaneous magnitude |P cos(th) - Q sin(th)|",
    "peak": "peak over the cycle: the ellipse's semi-major axis a",
    "axial_ratio": "axial ratio b/a    0 = linear, 1 = circular",
    "phase": "phase, degrees",
    "vector": "instantaneous vector  P cos(th) - Q sin(th)",
    "real": "real part P",
    "imag": "imaginary part Q",
}


# ---------------------------------------------------------------------------
# the derived quantity
#
# Computed in a ProgrammableFilter rather than a Calculator expression. `peak`
# and `axial_ratio` need a cross product and a discriminant, which a Calculator
# string can technically express and nobody can read; and `b` has to be computed
# the stable way (below), which a Calculator cannot do at all.
# ---------------------------------------------------------------------------
SCRIPT = r'''
import numpy as np
from paraview.numpy_support import vtk_to_numpy, numpy_to_vtk

inp = self.GetInputDataObject(0, 0)
out = self.GetOutputDataObject(0)
out.ShallowCopy(inp)

src = out.GetCellData() if USE_CELL else out.GetPointData()
P = vtk_to_numpy(src.GetArray(PREFIX + "_real")).astype(float)
Q = vtk_to_numpy(src.GetArray(PREFIX + "_imag")).astype(float)
th = math.radians(PHASE_DEG)

if P.ndim == 1:                       # a scalar field: phi
    P = P.reshape(-1, 1)
    Q = Q.reshape(-1, 1)

if DISPLAY == "vector":
    res = P * np.cos(th) - Q * np.sin(th)
elif DISPLAY == "real":
    res = P
elif DISPLAY == "imag":
    res = Q
elif DISPLAY == "magnitude_at_phase":
    res = np.linalg.norm(P * np.cos(th) - Q * np.sin(th), axis=1)
elif DISPLAY == "complex_magnitude":
    res = np.sqrt((P ** 2).sum(1) + (Q ** 2).sum(1))
elif DISPLAY == "phase":
    k = {"x": 0, "y": 1, "z": 2}.get(COMPONENT, 0)
    res = np.degrees(np.arctan2(Q[:, k], P[:, k]))
    # BLANK WHERE THERE IS NO FIELD TO HAVE A PHASE. Without this, an insulator
    # -- where J is exactly zero -- does not come out as "no data": it comes out
    # as a confident +180 degrees, because the stored real part is NEGATIVE ZERO
    # and atan2(+0.0, -0.0) is pi, not 0. Half of case 06's PP8 was that: a
    # uniform magenta field covering every air cell, indistinguishable from a
    # real measurement. NaN renders in the LUT's NaN colour instead, which reads
    # as absent rather than as 180.
    cmag = np.hypot(P[:, k], Q[:, k])
    big = cmag.max()
    if big > 0:
        res = np.where(cmag > 1e-9 * big, res, np.nan)
else:
    # peak and axial_ratio both need the ellipse's semi-axes.
    #
    #   S = (|P|^2+|Q|^2)/2,  C = (|P|^2-|Q|^2)/2,  D = P.Q,  R = hypot(C, D)
    #   a = sqrt(S + R)                       well conditioned
    #   b = |P x Q| / a                       STABLE
    #
    # b = sqrt(S - R) is the obvious formula and it is wrong to use: S ~ R
    # exactly when the field is nearly linear, which is most of any ordinary
    # domain, and the subtraction loses half the digits. The product invariant
    # a*b = |P x Q| avoids the cancellation entirely.
    N2 = (P ** 2).sum(1) + (Q ** 2).sum(1)
    S = 0.5 * N2
    C = 0.5 * ((P ** 2).sum(1) - (Q ** 2).sum(1))
    D = (P * Q).sum(1)
    a = np.sqrt(S + np.hypot(C, D))
    if DISPLAY == "peak":
        res = a
    else:
        # b = |P x Q| / a, so the RATIO b/a needs a twice. Dividing once gives
        # the semi-minor axis in the field's own units, which is a different
        # quantity and looks plausible on a colour bar -- case 06 is what caught
        # it: the analytic two-wire superposition says b/a must reach exactly
        # 1.0 at (0, +/-d/2), and this read 2.7e-05 there instead.
        cr = np.sqrt((np.cross(P, Q) ** 2).sum(1)) if P.shape[1] == 3 else np.zeros(len(a))
        safe = np.where(a > 0, a, 1.0)
        res = np.where(a > 0, cr / (safe * safe), 0.0)

arr = numpy_to_vtk(np.ascontiguousarray(res), deep=1)
arr.SetName("result")
(out.GetCellData() if USE_CELL else out.GetPointData()).AddArray(arr)
'''


def derived(source, req, use_cell):
    """A filter carrying the requested quantity as an array called `result`."""
    prefix = req["cell_prefix"] if use_cell else req["point_prefix"]
    pf = ProgrammableFilter(Input=source)
    pf.Script = ("import math\n"
                 "USE_CELL = %s\nPREFIX = %r\nDISPLAY = %r\n"
                 "PHASE_DEG = %r\nCOMPONENT = %r\n" %
                 (use_cell, str(prefix), str(req["display"]),
                  float(req.get("phase_deg", 0.0)), str(req.get("component", "x")))) + SCRIPT
    pf.RequestInformationScript = ""
    pf.RequestUpdateExtentScript = ""
    UpdatePipeline(proxy=pf)
    return pf


# ---------------------------------------------------------------------------
# scene
# ---------------------------------------------------------------------------
view = CreateView("RenderView")
view.ViewSize = [1100, 850]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.UseColorPaletteForBackground = 0


def clear():
    for s in GetSources().values():
        Hide(s, view)


def stamp(req, manifest, extra=""):
    """Say on the picture which quantity it is. A figure outlives its caption."""
    lines = ["%s   field %s   %s" % (req["name"], req["field"], req["geometry"]),
             WHAT.get(req["display"], req["display"])]
    if req["display"] in ("magnitude_at_phase", "vector"):
        lines.append("at phase wt = %g deg,  f = %g Hz" % (req.get("phase_deg", 0.0),
                                                           manifest["frequency_hz"]))
    else:
        lines.append("phase-independent,  f = %g Hz" % manifest["frequency_hz"])
    if req["field"] == "phi":
        lines.append("Phi is GAUGE-DEPENDENT -- only its jump across a port is physical")
    if extra:
        lines.append(extra)
    t = Text(Text="\n".join(lines))
    d = Show(t, view)
    d.WindowLocation = "Upper Left Corner"
    d.FontSize = 14
    d.Color = [0.1, 0.1, 0.1]
    return t


def colour(display_proxy, name, assoc, log_ok=True, cmap="rainbow"):
    """Colour by `result`, auto-ranged, log only when the span earns it."""
    info = display_proxy.Input.GetDataInformation()
    arr = (info.GetPointDataInformation() if assoc == "POINTS"
           else info.GetCellDataInformation()).GetArrayInformation("result")
    if arr is not None and arr.GetNumberOfComponents() > 1:
        ColorBy(display_proxy, (assoc, "result", "Magnitude"))
    else:
        ColorBy(display_proxy, (assoc, "result"))
    display_proxy.RescaleTransferFunctionToDataRange(True, False)
    lut = GetColorTransferFunction("result")
    # Put the control points back to linear BEFORE applying a preset. The LUT is
    # shared across every figure this script writes, and MapControlPointsToLogSpace
    # permanently respaces it -- a linear figure drawn after a log one otherwise
    # inherits log spacing with UseLogScale off and comes out quietly miscoloured.
    lut.UseLogScale = 0
    lut.MapControlPointsToLinearSpace()
    lut.ApplyPreset(PRESET.get(cmap, PRESET["rainbow"]), True)
    # WHITE, and deliberately a colour the scale cannot produce. The masked
    # region has no value, so it must not be drawn in anything the data could
    # take -- painting it at the bottom of the map (dark blue) would read as a
    # measurement at the low end of the range, which is the same mistake as the
    # false 180 degrees this masking exists to remove. White matches the page,
    # so the empty region simply is not there.
    lut.NanColor = [1.0, 1.0, 1.0]
    lut.NanOpacity = 1.0
    display_proxy.RescaleTransferFunctionToDataRange(True, False)
    lo, hi = lut.RGBPoints[0], lut.RGBPoints[-4]
    if log_ok and lo > 0 and hi / lo > 100.0:
        lut.MapControlPointsToLogSpace()
        lut.UseLogScale = 1
    bar = GetScalarBar(lut, view)
    bar.Title = name
    bar.ComponentTitle = ""
    bar.TitleColor = [0, 0, 0]
    bar.LabelColor = [0, 0, 0]
    bar.TitleFontSize = 12
    bar.LabelFontSize = 11
    display_proxy.SetScalarBarVisibility(view, True)


def shoot(path):
    Render()
    SaveScreenshot(path, view, ImageResolution=view.ViewSize)
    print("  wrote %s" % path)


# ---------------------------------------------------------------------------
def render_mesh_request(req, manifest):
    """geometry = body or plane -> an image."""
    vtk = os.path.join(OUT, req["vtk"])
    if not os.path.isfile(vtk):
        print("  SKIP %s: %s is not there" % (req["name"], req["vtk"]))
        return
    clear()
    src = LegacyVTKReader(FileNames=[vtk])

    # NODAL by default, per-tet on request. Point data is Gouraud-shaded so the
    # colour ramps across each triangle; cell data is one flat fill per element
    # and every tet reads as a facet. The manifest carries the choice; the
    # fallback is only for a manifest written before `data` existed.
    use_cell = req.get("data", "nodal") == "per_tet" and bool(req["cell_prefix"])

    if req["geometry"] == "body":
        base = Threshold(Input=src)
        base.Scalars = ["CELLS", "body_tag"]
        base.LowerThreshold = base.UpperThreshold = float(req["body_tag"])
        base.ThresholdMethod = "Between"
    else:
        base = src

    pf = derived(base, req, use_cell)

    if req["geometry"] == "plane":
        normal = {"xy": [0, 0, 1], "yz": [1, 0, 0], "zx": [0, 1, 0]}[req["plane"]]
        sl = Slice(Input=pf)
        sl.SliceType = "Plane"
        sl.SliceType.Normal = normal
        sl.SliceType.Origin = [req["offset_m"] * n for n in normal]
        UpdatePipeline(proxy=sl)
        shown = sl
        view.CameraPosition = [-n * 1.0 for n in normal] if req["plane"] == "yz" \
            else [n * 1.0 for n in normal]
        view.CameraViewUp = [0, 1, 0] if req["plane"] != "xy" else [0, 1, 0]
        if req["plane"] == "zx":
            view.CameraViewUp = [0, 0, 1]
    else:
        shown = pf
        view.CameraPosition = [1.0, -1.6, 1.0]
        view.CameraViewUp = [0, 0, 1]
    view.CameraFocalPoint = [0, 0, 0]

    assoc = "CELLS" if use_cell else "POINTS"
    if req["display"] in ("vector", "real", "imag") and req["is_vector"]:
        # A vector: arrows, coloured and scaled by their own length.
        #
        # GLYPHS LIVE ON POINTS. Handing Glyph a CELL array silently produces
        # nothing, so cell data is moved to the tet centroids first -- which is
        # also where the value actually belongs, B and E being constant per tet.
        d0 = Show(shown, view)
        d0.Representation = "Outline" if req["geometry"] == "body" else "Surface"
        d0.Opacity = 0.12
        seeds = shown
        if use_cell:
            seeds = CellCenters(Input=shown)
            seeds.VertexCells = 1
            UpdatePipeline(proxy=seeds)
        fetched = sm.Fetch(seeds)
        n = fetched.GetNumberOfPoints()
        vals = vtk_to_numpy(fetched.GetPointData().GetArray("result"))
        biggest = float(np.abs(vals).max()) if n else 0.0
        b = seeds.GetDataInformation().GetBounds()
        # Spelled out rather than sum(... for i in range(3)): `from
        # paraview.simple import *` REPLACES the builtins sum, min, max and abs
        # with VTK's array versions, which take a VTK array and reject a
        # generator. The failure is a baffling "'float' object has no attribute
        # 'astype'" from deep inside numpy_interface, nowhere near the call.
        diag = math.sqrt((b[1] - b[0]) ** 2 + (b[3] - b[2]) ** 2 + (b[5] - b[4]) ** 2)

        g = Glyph(Input=seeds, GlyphType="Arrow")
        g.OrientationArray = ["POINTS", "result"]
        g.ScaleArray = ["POINTS", "result"]
        # ScaleFactor multiplies the VALUE, so it has to carry the units away:
        # a field of 1e5 A/m^2 in a 2 cm box needs 1e-7, not 1. Left at its
        # default of 0 every arrow has zero length and the frame renders empty,
        # which is exactly how this first came out.
        g.ScaleFactor = (0.04 * diag / biggest) if biggest > 0 else 1.0
        g.GlyphMode = "Every Nth Point"
        # `max` is one of the builtins the star import above replaces, so this
        # is spelled with a conditional: max(1, k) reads 1 as an array and k as
        # an axis index, and fails with "axis 15 is out of bounds".
        g.Stride = n // 1500 if n > 1500 else 1     # ~1500 arrows whatever the mesh
        UpdatePipeline(proxy=g)
        gd = Show(g, view)
        colour(gd, "%s %s" % (req["field"], req["display"]), "POINTS", log_ok=False,
               cmap=req.get("colormap", "rainbow"))
    else:
        d0 = Show(shown, view)
        d0.Representation = "Surface"
        colour(d0, "%s %s" % (req["field"], req["display"]), assoc,
               log_ok=req["display"] not in ("phase", "real", "imag"),
               cmap=req.get("colormap", "rainbow"))

    stamp(req, manifest)
    ResetCamera(view)
    shoot(os.path.join(OUT, "%s_%s_%s.png" % (req["name"], req["field"], req["display"])))


def render_points_request(req, manifest):
    """geometry = points -> a line plot, not a render."""
    import vtk as _vtk
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    path = os.path.join(OUT, req["vtk"])
    if not os.path.isfile(path):
        print("  SKIP %s: %s is not there" % (req["name"], req["vtk"]))
        return
    src = LegacyVTKReader(FileNames=[path])
    # Probing interpolates POINT data, so the nodal arrays are what a probe can
    # answer with -- the cell arrays would give a piecewise-constant staircase.
    pf = derived(src, req, use_cell=False)
    UpdatePipeline(proxy=pf)
    data = sm.Fetch(pf)

    pts = _vtk.vtkPoints()
    for q in req["points_m"]:
        pts.InsertNextPoint(q[0], q[1], q[2])
    poly = _vtk.vtkPolyData()
    poly.SetPoints(pts)
    pr = _vtk.vtkProbeFilter()
    pr.SetInputData(poly)
    pr.SetSourceData(data)
    pr.Update()
    res = vtk_to_numpy(pr.GetOutput().GetPointData().GetArray("result"))
    valid = vtk_to_numpy(pr.GetOutput().GetPointData().GetArray("vtkValidPointMask"))
    if res.ndim > 1:
        res = np.linalg.norm(res, axis=1)

    q = np.array(req["points_m"], dtype=float)
    arc = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))])

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(arc * 1e3, res, "o-", lw=1.4, ms=5, color="#1f5fa9")
    out_of_mesh = (valid == 0)
    if out_of_mesh.any():
        ax.plot(arc[out_of_mesh] * 1e3, res[out_of_mesh], "x", ms=10, color="#c0522a",
                label="outside the mesh -- probe returned nothing")
        ax.legend(fontsize=8, frameon=False)
    ax.set_xlabel("distance along the probe path (mm)")
    ax.set_ylabel("%s  %s" % (req["field"], req["display"]))
    title = "%s   %s   %s" % (req["name"], req["field"], WHAT.get(req["display"], req["display"]))
    if req["display"] in ("magnitude_at_phase", "vector"):
        title += "\nat phase wt = %g deg,  f = %g Hz" % (req.get("phase_deg", 0.0),
                                                         manifest["frequency_hz"])
    else:
        title += "\nphase-independent,  f = %g Hz" % manifest["frequency_hz"]
    if req["field"] == "phi":
        title += "   (Phi is GAUGE-DEPENDENT)"
    ax.set_title(title, fontsize=9.5)
    ax.grid(alpha=0.25, lw=0.6)
    fig.tight_layout()
    name = os.path.join(OUT, "%s_%s_%s.png" % (req["name"], req["field"], req["display"]))
    fig.savefig(name, dpi=150, facecolor="white")
    plt.close(fig)
    print("  wrote %s" % name)


def main():
    if not os.path.isfile(MANIFEST):
        print("no %s -- the input file asked for no [postprocess] sections" % MANIFEST)
        return 0
    manifest = json.load(open(MANIFEST))
    print("%s: %d request(s) at %g Hz"
          % (MANIFEST, len(manifest["requests"]), manifest["frequency_hz"]))
    failed = 0
    for req in manifest["requests"]:
        try:
            if req["geometry"] == "points":
                render_points_request(req, manifest)
            else:
                render_mesh_request(req, manifest)
        except Exception as exc:                       # one bad request must not
            failed += 1                                # lose the others
            print("  FAILED %s: %s" % (req["name"], exc))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
