"""Figures for 01_OneCylinder: mesh, |Phi|, |B|, |E| on the mid-length plane,
and J as coloured vectors.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_plots.py
    pvbatch make_plots.py  output  output/plots        # explicit in/out

Run it from this folder after the case has been solved. Images land in
output/plots/. The script lives HERE rather than in output/ because output/ is
regenerated (and git-ignored) on every run.

Reading notes that the pictures cannot carry themselves:

  * |B| uses a LINEAR colour scale. A log ramp compresses the 1/r decay from
    the conductor surface to the box wall into the top few colours and hides
    the ring entirely -- the feature the plot exists to show.

  * |E| over the whole domain uses the PER-CELL array. It is exact and has no
    material interface to straddle; the price is that it renders faceted.
    Inside the wire the nodal array is used, where every tet is copper and the
    question does not arise.

  * J is drawn on a LONGITUDINAL slice, not the cross-section. The current is
    axial, so on a mid-length cut every arrow points at the viewer and nothing
    is visible.
"""

from paraview.simple import *
import os
import sys

HERE = sys.argv[1] if len(sys.argv) > 1 else "output"
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "plots")

ZMID = 20.0e-3     # mid-length of the 40 mm cylinder
A = 10.0e-3        # conductor radius
W = 200.0e-3       # box side
L = 40.0e-3        # cylinder length

if not os.path.isdir(OUT):
    os.makedirs(OUT)

# Plot-level smoothing, off by default. Each pass pushes the nodal field onto
# cells and back, replacing a node's value by an average over its element patch
# -- what a viewer's "smooth" option does. MEASURED on 02 for |B|:
#
#   azimuthal sd/mean        raw    1 pass  2 passes  4 passes
#     at the surface       0.080     0.054     0.037     0.020
#     in the far field     0.067     0.064     0.067     0.069
#   mean |B| at the surface  1.078 T   1.029     1.010     0.988
#
# It works near the conductor, where the scatter is high-frequency noise from a
# piecewise-constant B, and buys that by FLATTENING THE PEAK -- 8.3 % after four
# passes. In the far field it does nothing, because that scatter is the mesh
# being coarse against 1/r, which averaging cannot recover.
#
# So it is cosmetic where it helps and useless where it does not. Left at 0 so
# a value read off a figure is the value the solver produced.
# STALE. Measured on case 02, with bands 0.75 mm wide inside which |B| genuinely
# doubles, so radial variation was counted as scatter; and on case 02's mesh
# before its lc_core fix. See docs/FIELD_POSTPROCESSING.md, "The interior mesh
# constraint". The conclusion holds, the numbers should not be quoted. This case
# has the same core-sizing defect in milder form (about 4 elements per radius).
SMOOTH_PASSES = 0

# Colour-map resolution. None (or 0) means a genuinely continuous lookup table;
# an integer N gives N banded contours, which is what Ansys plots use by default.
# Both are honest -- this is a presentation choice and changes no number.
#
# Continuous is the standing choice: banding hides the per-element scatter by
# quantising it, and reads as a contour map rather than a physical field. The
# tradeoff is that B = curl A is CONSTANT per tetrahedron with first-order edge
# elements -- the lowest-order quantity in the formulation -- so it carries 2-6 %
# azimuthal scatter where an axisymmetric problem permits none, and a continuous
# ramp renders every bit of that as visible texture. See
# docs/FIELD_POSTPROCESSING.md. Pass legend(..., bands=11) for one banded figure.
BANDS = None

paraview.simple._DisableFirstRenderCameraReset()


def read(name):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        raise SystemExit("no such file: %s\nsolve the case first:\n"
                         "    ..\\run_case.bat cylinder_50hz.aphi" % p)
    return LegacyVTKReader(registrationName=name, FileNames=[p])


phi_src = read("potential.vtk")
b_src = read("B_field.vtk")
e_src = read("E_field.vtk")
j_src = read("J_field.vtk")

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1000, 900]
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.Background2 = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1

_shown = []


def clear():
    # Hiding a representation does NOT hide its colour legend: a stale one
    # stays drawn and the next lands on top of it.
    for d in _shown:
        try:
            d.SetScalarBarVisibility(view, False)
            d.Visibility = 0
        except Exception:
            pass
    del _shown[:]


def smoothed(src, passes=None):
    """Apply SMOOTH_PASSES point<->cell round trips. Identity when 0."""
    n = SMOOTH_PASSES if passes is None else passes
    cur = src
    for _ in range(n):
        p2c = PointDatatoCellData(Input=cur)
        p2c.ProcessAllArrays = 1
        cur = CellDatatoPointData(Input=p2c)
        cur.ProcessAllArrays = 1
    return cur

def cut(src, normal, origin, only_wire=False, passes=None):
    s = src
    if only_wire:
        s = Threshold(Input=src)
        s.Scalars = ["CELLS", "body_tag"]
        s.LowerThreshold = 1.0
        s.UpperThreshold = 1.0
        s.ThresholdMethod = "Between"
    sl = Slice(Input=smoothed(s, passes))
    sl.SliceType = "Plane"
    sl.SliceType.Origin = origin
    sl.SliceType.Normal = normal
    return sl


def legend(array, title, fmt="{:.3g}", bands=None):
    """Colour map and scalar bar. `bands` defaults to the BANDS constant above;
    pass 11 for banded contours on a single figure."""
    lut = GetColorTransferFunction(array)
    # Rainbow Uniform, not Jet. They look alike, but Jet is not perceptually
    # uniform: it compresses the greens and stretches the cyans, so it draws
    # contour-like bands that are not in the data. Same blue-low-to-red-high
    # reading, without inventing structure. Cases 04 and 05 already use it,
    # and so does the shared tools/postprocess.py.
    lut.ApplyPreset("Rainbow Uniform", True)   # blue low -> red high
    n = BANDS if bands is None else bands
    if n:
        lut.Discretize = 1
        lut.NumberOfTableValues = n
    else:
        lut.Discretize = 0        # a true continuous LUT, not 256 steps
    lut.UseLogScale = 0
    bar = GetScalarBar(lut, view)
    bar.Title = title
    bar.ComponentTitle = ""
    bar.WindowLocation = "Any Location"
    bar.Position = [0.855, 0.22]
    bar.ScalarBarLength = 0.56
    bar.ScalarBarThickness = 20
    bar.TitleColor = [0, 0, 0]
    bar.TitleFontSize = 20
    bar.LabelFontSize = 16
    bar.LabelColor = [0, 0, 0]
    bar.AutomaticLabelFormat = 1
    bar.AddRangeLabels = 1
    bar.RangeLabelFormat = fmt
    return lut


def look(half, along="z"):
    view.CameraParallelScale = half
    if along == "z":
        view.CameraFocalPoint = [0.0, 0.0, ZMID]
        view.CameraPosition = [0.0, 0.0, ZMID + 1.0]
        view.CameraViewUp = [0.0, 1.0, 0.0]
    else:                                   # looking along -y, z upward
        view.CameraFocalPoint = [0.0, 0.0, L / 2.0]
        view.CameraPosition = [0.0, -1.0, L / 2.0]
        view.CameraViewUp = [0.0, 0.0, 1.0]


def shot(name):
    Render(view)
    p = os.path.join(OUT, name)
    SaveScreenshot(p, view, ImageResolution=view.ViewSize)
    print("  wrote %-22s %8d bytes" % (name, os.path.getsize(p)))


def surface(src, array, assoc, title, comp=None, edges=False, flat=True):
    d = Show(src, view, "GeometryRepresentation")
    d.Representation = "Surface With Edges" if edges else "Surface"
    if array:
        ColorBy(d, (assoc, array, comp) if comp else (assoc, array))
        d.RescaleTransferFunctionToDataRange(True, False)
        d.SetScalarBarVisibility(view, True)
        legend(array, title)
    else:
        ColorBy(d, None)
        # Surface With Edges takes its solid colour from AmbientColor when
        # lighting is off, and its lines from EdgeColor. Light fill, dark
        # lines -- the reverse prints as a black block.
        d.AmbientColor = [0.93, 0.94, 0.96]
        d.DiffuseColor = [0.93, 0.94, 0.96]
        d.EdgeColor = [0.10, 0.12, 0.16]
        d.LineWidth = 1.0
    d.NonlinearSubdivisionLevel = 1 if edges else 3
    if flat:
        d.Ambient = 1.0
        d.Diffuse = 0.0
        d.Specular = 0.0
    _shown.append(d)
    return d


# --- 1. the mesh, whole domain ----------------------------------------------
clear()
look(0.5 * W * 1.05)
surface(cut(phi_src, [0, 0, 1], [0, 0, ZMID]), None, None, None, edges=True)
shot("01_mesh_domain.png")

# --- 2. the mesh, the wire alone --------------------------------------------
clear()
look(A * 1.25)
surface(cut(phi_src, [0, 0, 1], [0, 0, ZMID], only_wire=True), None, None, None, edges=True)
shot("02_mesh_wire.png")

# --- 3. Phi over the cylinder SURFACE --------------------------------------
# Not a cross-section: the potential is driven along the axis, so the face that
# shows it is the lateral surface, seen side on.
clear()
look(0.62 * L, along="y")
wire3d = Threshold(Input=phi_src)
wire3d.Scalars = ["CELLS", "body_tag"]
wire3d.LowerThreshold = 1.0
wire3d.UpperThreshold = 1.0
wire3d.ThresholdMethod = "Between"
surface(wire3d, "phi_real", "POINTS", "Re(Phi)  (V)")
# Diverging, centred on 0.5 V. Jet spends its last third on oranges and reds that
# read as one colour, so a smooth ramp renders as two blocks; a diverging map puts
# its resolution where Phi is actually changing. Range pinned to the two applied
# terminal voltages so the midpoint IS 0.5 rather than wherever the data lands.
_p = GetColorTransferFunction("phi_real")
_p.ApplyPreset("Cool to Warm", True)
_p.RescaleTransferFunction(0.0, 1.0)
shot("03_phi_on_surface.png")

clear()
look(0.62 * L, along="y")
surface(wire3d, "phi_magnitude", "POINTS", "Phi magnitude (V)")
_q = GetColorTransferFunction("phi_magnitude")
_q.ApplyPreset("Cool to Warm", True)
_q.RescaleTransferFunction(0.0, 1.0)
shot("03b_phi_magnitude_surface.png")

# --- 4. |B|, LINEAR scale ----------------------------------------------------
clear()
look(0.5 * W * 1.05)
d = surface(cut(b_src, [0, 0, 1], [0, 0, ZMID]), "B_magnitude", "POINTS", "B magnitude (T)")
lut = GetColorTransferFunction("B_magnitude")
lut.MapControlPointsToLinearSpace()
shot("04_B_magnitude.png")

# --- 4b. the same, zoomed on the conductor ---------------------------------
# The box is 20x the wire radius, so at full extent the ring is a few pixels.
clear()
look(3.0 * A)
surface(cut(b_src, [0, 0, 1], [0, 0, ZMID]), "B_magnitude", "POINTS", "B magnitude (T)")
GetColorTransferFunction("B_magnitude").MapControlPointsToLinearSpace()
shot("04b_B_magnitude_zoom.png")

# --- 5. |E|, per cell over the whole domain ---------------------------------
clear()
look(0.5 * W * 1.05)
surface(cut(e_src, [0, 0, 1], [0, 0, ZMID]), "E_cell_real", "CELLS", "E magnitude (V/m)",
        comp="Magnitude")
lut = GetColorTransferFunction("E_cell_real")
lut.MapControlPointsToLogSpace()
lut.UseLogScale = 1          # E spans 1.4 decades of smooth decay -- log suits it
shot("05_E_magnitude.png")

# --- 5b. the same, zoomed --------------------------------------------------
clear()
look(3.0 * A)
surface(cut(e_src, [0, 0, 1], [0, 0, ZMID]), "E_cell_real", "CELLS", "E magnitude (V/m)",
        comp="Magnitude")
_l = GetColorTransferFunction("E_cell_real")
_l.MapControlPointsToLogSpace(); _l.UseLogScale = 1
shot("05b_E_magnitude_zoom.png")

# --- 6. |E| in the wire alone, nodal ----------------------------------------
clear()
look(A * 1.25)
surface(cut(e_src, [0, 0, 1], [0, 0, ZMID], only_wire=True), "E_magnitude", "POINTS",
        "E magnitude (V/m)")
shot("06_E_in_wire.png")

# --- 7. J as vectors, coloured by magnitude ---------------------------------
# Longitudinal cut: the current is axial, so a cross-section shows arrows end on.
clear()
look(0.62 * L, along="y")
wire = Threshold(Input=j_src)
wire.Scalars = ["CELLS", "body_tag"]
wire.LowerThreshold = 1.0
wire.UpperThreshold = 1.0
wire.ThresholdMethod = "Between"
plane = Slice(Input=wire)
plane.SliceType = "Plane"
plane.SliceType.Origin = [0.0, 0.0, L / 2.0]
plane.SliceType.Normal = [0.0, 1.0, 0.0]

back = Show(plane, view, "GeometryRepresentation")
back.Representation = "Surface"
ColorBy(back, None)
back.AmbientColor = [0.9, 0.9, 0.9]
back.DiffuseColor = [0.93, 0.93, 0.93]
back.Ambient = 1.0
back.Diffuse = 0.0
_shown.append(back)

g = Glyph(registrationName="J", Input=plane, GlyphType="Arrow")
g.OrientationArray = ["POINTS", "J_real"]
g.ScaleArray = ["POINTS", "No scale array"]
g.ScaleFactor = 0.055 * L
g.GlyphMode = "All Points"
gd = Show(g, view, "GeometryRepresentation")
gd.Representation = "Surface"
ColorBy(gd, ("POINTS", "J_magnitude"))
gd.RescaleTransferFunctionToDataRange(True, False)
gd.SetScalarBarVisibility(view, True)
gd.Ambient = 1.0
gd.Diffuse = 0.0
legend("J_magnitude", "J magnitude (A/m2)")
_shown.append(gd)
shot("07_J_vectors.png")

print("\n  all figures in %s" % OUT)
