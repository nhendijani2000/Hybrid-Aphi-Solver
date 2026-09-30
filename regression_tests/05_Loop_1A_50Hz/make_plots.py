"""Figures for 05_Loop_1A_50Hz.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_plots.py [output-dir]

Writes into fig/. Four groups:

  * the POTENTIAL over the whole torus, seen from +z. This is the picture that
    shows what an internal cut actually does: Phi ramps smoothly all the way
    round the ring and falls off a cliff at the cut, which is the terminal
    voltage. Nothing in cases 01-04 looks like this.
  * the CURRENT as vectors on the ring, to show it circulating.
  * E and B on the z = 0 plane (the xy plane through the tube centreline),
    zoomed to the ring and again over the whole domain.
  * the y = 0 plane through the PORT, which cuts the tube twice and shows the
    cross-section the cut is made on.

A note on colour ranges, which bit hard in case 04: |E| inside a conductor with
no skin effect is nearly uniform, so an auto range there magnifies noise in the
fifth digit into a full rainbow. Where that applies the range is pinned and the
caption says so.
"""
from paraview.simple import *
from paraview import servermanager as sm
import os
import sys
import math

paraview.simple._DisableFirstRenderCameraReset()

OUT = sys.argv[1] if len(sys.argv) > 1 else "output"
FIG = "fig"
if not os.path.isdir(FIG):
    os.makedirs(FIG)

R, a, Rd, Hd = 6.5e-3, 0.8e-3, 30e-3, 30e-3

# A slice exactly at z = 0 is DEGENERATE here. The tube cross-section is a
# 20-gon with vertices at z = +/- a sin(theta), two of which land exactly on
# z = 0, so that plane grazes a whole ring of coplanar faces and the slice
# comes out ragged -- visibly so in the lower half, which looked at first
# like an asymmetric mesh. Offsetting by an eighth of a tube radius keeps the
# cut essentially through the centreline while missing every vertex.
ZCUT = 0.125 * a

view = CreateView("RenderView")
view.ViewSize = [1100, 1000]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.Background = [1, 1, 1]
view.UseColorPaletteForBackground = 0


def clear():
    for s in GetSources().values():
        try:
            Hide(s, view)
        except Exception:
            pass


def look(scale, along="z"):
    """Park the camera on an axis, looking at the origin."""
    d = 0.2
    if along == "z":
        view.CameraPosition = [0, 0, d]
        view.CameraViewUp = [0, 1, 0]
    elif along == "y":
        view.CameraPosition = [0, d, 0]
        view.CameraViewUp = [0, 0, 1]
    else:
        view.CameraPosition = [d, 0, 0]
        view.CameraViewUp = [0, 0, 1]
    view.CameraFocalPoint = [0, 0, 0]
    view.CameraParallelScale = scale


def src(fn):
    return LegacyVTKReader(FileNames=[os.path.join(OUT, fn)])


def ring_of(fn):
    t = Threshold(Input=src(fn))
    t.Scalars = ["CELLS", "body_tag"]
    t.LowerThreshold = 1.0
    t.UpperThreshold = 1.0
    t.ThresholdMethod = "Between"
    return t


def bar(lut, title, pos=(0.86, 0.30), length=0.42):
    b = GetScalarBar(lut, view)
    b.Title = title
    b.ComponentTitle = ""
    b.TitleColor = [0, 0, 0]
    b.LabelColor = [0, 0, 0]
    b.WindowLocation = "Any Location"
    b.Position = list(pos)
    b.ScalarBarLength = length
    return b


def shot(name, w=1100, h=1000):
    Render(view)
    SaveScreenshot(os.path.join(FIG, name), view, ImageResolution=[w, h])
    print("wrote %s/%s" % (FIG, name))


# ---------------------------------------------------------------------------
# 1. POTENTIAL over the torus, from +z -- the picture of what a cut does
# ---------------------------------------------------------------------------
clear()
look(1.30 * (R + a))
w = ring_of("potential.vtk")
UpdatePipeline(proxy=w)
d = Show(w, view)
ColorBy(d, ("POINTS", "phi_real"))
d.Representation = "Surface"
# Flat, unlit. With normal 3D lighting one side of the tube darkens by more
# than Phi changes across a colour band, which reads as physics and is not.
d.Ambient = 1.0
d.Diffuse = 0.0
d.Specular = 0.0
d.RescaleTransferFunctionToDataRange(True, False)
lut = GetColorTransferFunction("phi_real")
lut.ApplyPreset("Rainbow Uniform", True)
bar(lut, "Re(Phi)  V")
d.SetScalarBarVisibility(view, True)
shot("phi_ring_top.png")

# the same from an angle, so the tube reads as a tube
clear()
d = Show(w, view)
ColorBy(d, ("POINTS", "phi_real"))
d.Representation = "Surface"
d.Ambient = 1.0
d.Diffuse = 0.0
d.Specular = 0.0
d.SetScalarBarVisibility(view, True)
view.CameraPosition = [0.9 * (R + a), -1.5 * (R + a), 1.1 * (R + a)]
view.CameraFocalPoint = [0, 0, 0]
view.CameraViewUp = [0, 0, 1]
view.CameraParallelScale = 1.25 * (R + a)
shot("phi_ring_iso.png")

# ---------------------------------------------------------------------------
# 2. CURRENT vectors on the ring
# ---------------------------------------------------------------------------
clear()
look(1.30 * (R + a))
wj = ring_of("J_field.vtk")
cc = CellCenters(Input=wj)
cc.VertexCells = 1
UpdatePipeline(proxy=cc)
g = Glyph(Input=cc, GlyphType="Arrow")
g.OrientationArray = ["POINTS", "J_cell_real"]
g.ScaleArray = ["POINTS", "No scale array"]
g.ScaleFactor = 0.9e-3
g.GlyphMode = "Every Nth Point"
g.Stride = 11            # 16177 cells would be a solid smear
UpdatePipeline(proxy=g)
dj = Show(g, view)
ColorBy(dj, ("POINTS", "J_cell_real", "Magnitude"))
dj.RescaleTransferFunctionToDataRange(True, False)
lutj = GetColorTransferFunction("J_cell_real")
lutj.ApplyPreset("Rainbow Uniform", True)
bar(lutj, "|Re(J)|  A/m2")
dj.SetScalarBarVisibility(view, True)
# a faint ring underneath, so the arrows have something to sit on
dw = Show(ring_of("J_field.vtk"), view)
ColorBy(dw, None)
dw.Representation = "Surface"
dw.DiffuseColor = [0.88, 0.89, 0.92]
dw.Opacity = 0.35
dw.SetScalarBarVisibility(view, False)
shot("j_vectors.png")

# ---------------------------------------------------------------------------
# 3. E and B on the z = 0 plane (the xy plane through the tube centreline)
# ---------------------------------------------------------------------------
def plane_shot(fn, array, title, name, scale, preset="Rainbow Uniform",
               logscale=False, only_ring=False, pin=None):
    clear()
    look(scale)
    base = ring_of(fn) if only_ring else src(fn)
    calc = Calculator(Input=base)
    calc.AttributeType = "Cell Data"
    calc.ResultArrayName = "mag"
    calc.Function = ("sqrt(%s_real_X^2+%s_real_Y^2+%s_real_Z^2"
                     "+%s_imag_X^2+%s_imag_Y^2+%s_imag_Z^2)"
                     % (array, array, array, array, array, array))
    sl = Slice(Input=calc)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, ZCUT]
    sl.SliceType.Normal = [0.0, 0.0, 1.0]
    UpdatePipeline(proxy=sl)
    dd = Show(sl, view)
    ColorBy(dd, ("CELLS", "mag"))
    dd.RescaleTransferFunctionToDataRange(True, False)
    l = GetColorTransferFunction("mag")
    # EVERY "mag" FIGURE SHARES THIS ONE LOOKUP TABLE, and
    # MapControlPointsToLogSpace PERMANENTLY respaces its control points. A
    # linear figure drawn after a log one therefore inherits log-spaced control
    # points with UseLogScale off, and its colours are simply wrong -- washed
    # out mid-range, which is exactly how the port-plane figure first came out.
    # Put it back to linear FIRST, then choose.
    l.UseLogScale = 0
    l.MapControlPointsToLinearSpace()
    l.ApplyPreset(preset, True)
    dd.RescaleTransferFunctionToDataRange(True, False)   # on the DISPLAY, not the LUT
    if pin is not None:
        l.RescaleTransferFunction(pin[0], pin[1])
    if logscale:
        l.MapControlPointsToLogSpace()
        l.UseLogScale = 1
    bar(l, title)
    dd.SetScalarBarVisibility(view, True)
    shot(name)


plane_shot("B_field.vtk", "B_cell", "|B|  T", "b_xy_ring.png", 1.6 * (R + a))
plane_shot("B_field.vtk", "B_cell", "|B|  T", "b_xy_domain.png", 1.05 * Rd,
           logscale=True)
# PIN |E| IN THE RING VIEW. Phi jumps across the cut, which is a zero-thickness
# interface, so grad(Phi) -- and therefore |E| -- SPIKES there: 1.5 V/m against
# the conductor's own J/sigma = 8.6e-3 V/m, a factor of 175. On an auto range
# that one spike takes the whole colour map and the conductor reads as flat
# blue. Pinned to the conductor's own scale, the cut simply saturates and the
# picture shows what it should. Pinned to 3.0e-2, which spans BOTH the
# conductor's J/sigma = 8.6e-3 and the surrounding air's 2.0e-2 -- the air is
# HIGHER, and not because of -j*omega*A (that is only 2.6e-4): it is the
# electrostatic field of the ring's own varying surface potential, 0.354 mV
# from one side of the cut to the other.
plane_shot("E_field.vtk", "E_cell", "|E|  V/m", "e_xy_ring.png", 1.6 * (R + a),
           pin=(0.0, 3.0e-2))
plane_shot("E_field.vtk", "E_cell", "|E|  V/m", "e_xy_domain.png", 1.05 * Rd,
           logscale=True)

# ---------------------------------------------------------------------------
# 4. the y = 0 plane THROUGH THE PORT -- it cuts the tube twice
# ---------------------------------------------------------------------------
clear()
look(1.6 * (R + a), along="y")
calc = Calculator(Input=src("B_field.vtk"))
calc.AttributeType = "Cell Data"
calc.ResultArrayName = "mag"
calc.Function = ("sqrt(B_cell_real_X^2+B_cell_real_Y^2+B_cell_real_Z^2"
                 "+B_cell_imag_X^2+B_cell_imag_Y^2+B_cell_imag_Z^2)")
sl = Slice(Input=calc)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]
UpdatePipeline(proxy=sl)
dd = Show(sl, view)
ColorBy(dd, ("CELLS", "mag"))
dd.RescaleTransferFunctionToDataRange(True, False)
l = GetColorTransferFunction("mag")
l.UseLogScale = 0
l.MapControlPointsToLinearSpace()   # see the note in plane_shot
l.ApplyPreset("Rainbow Uniform", True)
# PIN THE RANGE. The slice spans the whole 30 mm domain but the view is zoomed
# to 12 mm, so an auto range is set by a far-field minimum three decades below
# anything visible and the entire picture lands in the top decade -- uniformly
# red, which is how this first came out. 0 .. 3.2e-4 is the near field's own
# range and spends the colours where the eye is looking.
l.RescaleTransferFunction(0.0, 3.2e-4)
l.UseLogScale = 0
bar(l, "|B|  T")
dd.SetScalarBarVisibility(view, True)
shot("b_port_plane.png")

# ---------------------------------------------------------------------------
# 5. the mesh
# ---------------------------------------------------------------------------
clear()
look(1.05 * Rd)
ms = Slice(Input=src("J_field.vtk"))
ms.SliceType = "Plane"
ms.SliceType.Origin = [0.0, 0.0, ZCUT]
ms.SliceType.Normal = [0.0, 0.0, 1.0]
UpdatePipeline(proxy=ms)
dm = Show(ms, view)
ColorBy(dm, None)
dm.Representation = "Surface With Edges"
dm.DiffuseColor = [0.86, 0.88, 0.92]
dm.EdgeColor = [0.15, 0.17, 0.22]
dm.SetScalarBarVisibility(view, False)
shot("mesh_domain.png")

clear()
look(1.35 * (R + a))
mr = Slice(Input=ring_of("J_field.vtk"))
mr.SliceType = "Plane"
mr.SliceType.Origin = [0.0, 0.0, ZCUT]
mr.SliceType.Normal = [0.0, 0.0, 1.0]
UpdatePipeline(proxy=mr)
dmr = Show(mr, view)
ColorBy(dmr, None)
dmr.Representation = "Surface With Edges"
dmr.DiffuseColor = [0.86, 0.88, 0.92]
dmr.EdgeColor = [0.15, 0.17, 0.22]
dmr.SetScalarBarVisibility(view, False)
shot("mesh_ring.png")

print("done")
