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

R, a, Rd, Hd = 6.5e-3, 0.8e-3, 90e-3, 90e-3

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
        # From -y, NOT +y. Looking from +y puts +x on the LEFT of the screen,
        # so the port at theta = 0 lands on the left and the plain joint at
        # theta = pi on the right, which reads backwards. From -y, +x is on the
        # right and the port is the right-hand crossing.
        view.CameraPosition = [0, -d, 0]
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
               logscale=False, only_ring=False, pin=None, nodal=True):
    """`nodal` picks which array the colour comes from, and it is a real choice.

    NODAL (point data) is Gouraud-shaded: the colour ramps across each triangle
    between its vertex values, so the picture is smooth. PER-TET (cell data) is
    one flat fill per element, so every tet reads as a facet.

    Smooth is better everywhere EXCEPT at a material interface, and there it is
    actively wrong. A node sitting on the conductor surface is touched by tets on
    both sides, and the nodal value is their volume-weighted average -- a blend of
    the conductor's 8.7e-03 V/m and the air's 4.2e-02, two numbers that are
    different because the field genuinely JUMPS there. Which tets happen to touch
    each node varies around the ring, so the blend varies, and the boundary comes
    out as a ragged fringe instead of the clean discontinuity it is. The solver
    already flags those nodes: 14,736 of them, 5 % of the mesh, carry
    `material_interface = 1`.

    So: B is drawn nodal (mu_r = 1 in both bodies, so there is no interface for it
    to fall over), and E is drawn BOTH ways -- per-tet for the report's own
    figures, and nodal beside it so the artefact can be seen rather than described.
    Every validated number comes from the per-cell arrays via
    analytic_comparison.py and is untouched by any of this.
    """
    clear()
    look(scale)
    base = ring_of(fn) if only_ring else src(fn)
    pre = array if nodal else array + "_cell"
    calc = Calculator(Input=base)
    calc.AttributeType = "Point Data" if nodal else "Cell Data"
    calc.ResultArrayName = "mag"
    calc.Function = ("sqrt(%s_real_X^2+%s_real_Y^2+%s_real_Z^2"
                     "+%s_imag_X^2+%s_imag_Y^2+%s_imag_Z^2)"
                     % (pre, pre, pre, pre, pre, pre))
    sl = Slice(Input=calc)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, ZCUT]
    sl.SliceType.Normal = [0.0, 0.0, 1.0]
    UpdatePipeline(proxy=sl)
    dd = Show(sl, view)
    ColorBy(dd, ("POINTS" if nodal else "CELLS", "mag"))
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
    l.Discretize = 0   # continuous, not 256 steps
    dd.RescaleTransferFunctionToDataRange(True, False)   # on the DISPLAY, not the LUT
    if pin is not None:
        l.RescaleTransferFunction(pin[0], pin[1])
    if logscale:
        l.MapControlPointsToLogSpace()
        l.UseLogScale = 1
    bar(l, title)
    dd.SetScalarBarVisibility(view, True)
    shot(name)


plane_shot("B_field.vtk", "B", "|B|  T", "b_xy_ring.png", 1.6 * (R + a))
# PINNED, NOT AUTO-RANGED, and for a reason that only shows up now the domain
# is 90 mm: an auto range is set by the single faintest cell in the frame, which
# sits in a corner and moves whenever the domain or the grading changes. The
# 90 mm and 30 mm figures would then be drawn on different scales and could not
# honestly be compared. These bounds span both meshes: |B| runs 1.3e-08 to
# 3.2e-04 T on the 90 mm cut.
plane_shot("B_field.vtk", "B", "|B|  T", "b_xy_domain.png", 1.05 * Rd,
           logscale=True, pin=(1.0e-8, 4.0e-4))
# PIN |E| IN THE RING VIEW. Phi jumps across the cut, which is a zero-thickness
# interface, so grad(Phi) -- and therefore |E| -- SPIKES there: 1.5 V/m against
# the conductor's own J/sigma = 8.6e-3 V/m, a factor of 175. On an auto range
# that one spike takes the whole colour map and the conductor reads as flat
# blue. Pinned to the conductor's own scale, the cut simply saturates and the
# picture shows what it should. 3.0e-2 spans the conductor's J/sigma = 8.6e-3
# and the air around the port.
#
# THE AIR IS NOT GENERALLY BRIGHTER THAN THE CONDUCTOR, which an earlier version
# of this comment claimed. On the joint side the air just outside the tube is
# 0.0081 V/m against the conductor's 0.0087 -- LOWER. It is only near the port
# that it runs high, and that is the cut, not the ring. What the air's field is
# NOT is -j*omega*A, which is 2.6e-4 here; it is electrostatic, from the
# 0.354 mV the loop carries across the cut. See section 4.3.
# Drawn BOTH ways. The per-tet one is the figure the report uses, because the
# conductor surface is a real discontinuity and per-tet renders it as one; the
# nodal one sits beside it in section 4.5 so the fringe can be seen.
plane_shot("E_field.vtk", "E", "|E|  V/m", "e_xy_ring.png", 1.6 * (R + a),
           pin=(0.0, 3.0e-2), nodal=False)
plane_shot("E_field.vtk", "E", "|E|  V/m", "e_xy_ring_nodal.png", 1.6 * (R + a),
           pin=(0.0, 3.0e-2), nodal=True)
plane_shot("E_field.vtk", "E", "|E|  V/m", "e_xy_domain.png", 1.05 * Rd,
           logscale=True, pin=(1.0e-5, 2.0), nodal=False)

# ---------------------------------------------------------------------------
# 4. the y = 0 plane THROUGH THE PORT -- it cuts the tube twice
# ---------------------------------------------------------------------------
clear()
look(1.6 * (R + a), along="y")
calc = Calculator(Input=src("B_field.vtk"))
calc.AttributeType = "Point Data"
calc.ResultArrayName = "mag"
calc.Function = ("sqrt(B_real_X^2+B_real_Y^2+B_real_Z^2"
                 "+B_imag_X^2+B_imag_Y^2+B_imag_Z^2)")
sl = Slice(Input=calc)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]
UpdatePipeline(proxy=sl)
dd = Show(sl, view)
ColorBy(dd, ("POINTS", "mag"))
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


# ---------------------------------------------------------------------------
# 6. THE GEOMETRY ITSELF -- the ring inside the air domain
# ---------------------------------------------------------------------------
# What was actually modelled, which every other figure only shows a slice of.
# The domain boundary is drawn transparent so the ring inside is visible;
# ExtractSurface on the FULL dataset returns only the outer faces, because the
# ring is interior to the mesh and its faces are shared.
def geometry(name, clip_half, elev_scale):
    clear()
    dom = ExtractSurface(Input=src("J_field.vtk"))
    if clip_half:
        c = Clip(Input=dom)
        c.ClipType = "Plane"
        c.ClipType.Origin = [0.0, 0.0, 0.0]
        c.ClipType.Normal = [0.0, 1.0, 0.0]
        c.Invert = 0
        dom = c
    UpdatePipeline(proxy=dom)
    dd = Show(dom, view)
    ColorBy(dd, None)
    dd.Representation = "Surface"
    dd.DiffuseColor = [0.62, 0.66, 0.82]
    dd.Opacity = 0.16
    dd.SetScalarBarVisibility(view, False)

    rg = ring_of("J_field.vtk")
    UpdatePipeline(proxy=rg)
    dr = Show(rg, view)
    ColorBy(dr, None)
    dr.Representation = "Surface"
    dr.DiffuseColor = [0.72, 0.42, 0.20]      # copper
    dr.Opacity = 1.0
    dr.SetScalarBarVisibility(view, False)

    view.CameraPosition = [1.15 * Rd, -1.55 * Rd, elev_scale * Rd]
    view.CameraFocalPoint = [0, 0, 0]
    view.CameraViewUp = [0, 0, 1]
    view.CameraParallelScale = 1.02 * Rd
    shot(name)


geometry("geometry.png", False, 0.85)
geometry("geometry_cut.png", True, 0.85)


# ---------------------------------------------------------------------------
# 7. E ON THE PLANE NORMAL TO THE TORUS THROUGH THE PORT
# ---------------------------------------------------------------------------
# "Normal to the torus" at theta = 0 means the plane whose normal is the
# centreline tangent there, which is +y -- so it is the y = 0 plane, and that
# is the CUT PLANE itself. It slices the tube twice: at the port (theta = 0,
# +x) and at theta = pi (-x).
#
# THE SLICE SITS EXACTLY ON A DISCONTINUITY. Phi jumps across the cut, so E
# does too, and a plane lying exactly in it samples cells from both sides. B is
# continuous there and does not care (section 4.2); E very much does. Both are
# drawn: the plane itself, and a slice one tube radius to the +y side of it
# where the field is single valued, so the two can be compared.
def e_on_plane(name, yoff, focal, scale, pin, title, logscale=False):
    clear()
    # per-tet: these cut the conductor surface too, and the fringe above is why
    calc = Calculator(Input=src("E_field.vtk"))
    calc.AttributeType = "Cell Data"
    calc.ResultArrayName = "emag"
    calc.Function = ("sqrt(E_cell_real_X^2+E_cell_real_Y^2+E_cell_real_Z^2"
                     "+E_cell_imag_X^2+E_cell_imag_Y^2+E_cell_imag_Z^2)")
    sl = Slice(Input=calc)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, yoff, 0.0]
    sl.SliceType.Normal = [0.0, 1.0, 0.0]
    UpdatePipeline(proxy=sl)
    dd = Show(sl, view)
    ColorBy(dd, ("CELLS", "emag"))
    l = GetColorTransferFunction("emag")
    l.UseLogScale = 0
    l.MapControlPointsToLinearSpace()
    l.ApplyPreset("Rainbow Uniform", True)
    l.RescaleTransferFunction(pin[0], pin[1])
    if logscale:
        l.MapControlPointsToLogSpace()
        l.UseLogScale = 1
    bar(l, title)
    dd.SetScalarBarVisibility(view, True)
    view.CameraPosition = [focal[0], focal[1] - 0.2, focal[2]]   # from -y: +x on the right
    view.CameraFocalPoint = list(focal)
    view.CameraViewUp = [0, 0, 1]
    view.CameraParallelScale = scale
    shot(name)


# the whole plane: both tube cross-sections
e_on_plane("e_port_plane.png", 0.0, (0, 0, 0), 1.6 * (R + a),
           (0.0, 3.0e-2), "|E|  V/m")
# THE SAME CUT, FAR OUT. Log scaled and it has to be: |E| runs from 1.95 V/m at
# the cut face to 1.8e-5 V/m at the 90 mm wall, five decades.
#
# FRAMED AT 30 mm, NOT AT Rd, AND THAT IS DELIBERATE. The y = 0 cut is 2Rd wide
# by Hd tall -- 180 x 90 mm now -- and the render is square, so framing it at
# the full domain puts the whole ring inside about 8 % of the frame height and
# the picture becomes a dot on a blue field. The square xy views above do not
# have this problem and are still framed at the full Rd, which is where to look
# for "the far field is dead at the wall". Here 30 mm is the window where
# anything happens, and it keeps this figure directly comparable with the one
# the 30 mm control produces. The caption says which window it is.
FAR = 30e-3
e_on_plane("e_port_domain.png", 0.0, (0, 0, 0), 1.03 * FAR,
           (1.0e-5, 2.0), "|E|  V/m", logscale=True)
# B on the same cut over the whole domain, for comparison -- it is continuous
# across the port, so it has none of E's trouble there.
def b_on_plane(name, scale, pin, logscale=True):
    clear()
    calc = Calculator(Input=src("B_field.vtk"))
    calc.AttributeType = "Point Data"
    calc.ResultArrayName = "bmagp"
    calc.Function = ("sqrt(B_real_X^2+B_real_Y^2+B_real_Z^2"
                     "+B_imag_X^2+B_imag_Y^2+B_imag_Z^2)")
    sl = Slice(Input=calc)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, 0.0]
    sl.SliceType.Normal = [0.0, 1.0, 0.0]
    UpdatePipeline(proxy=sl)
    dd = Show(sl, view)
    ColorBy(dd, ("POINTS", "bmagp"))
    l = GetColorTransferFunction("bmagp")
    l.UseLogScale = 0
    l.MapControlPointsToLinearSpace()
    l.ApplyPreset("Rainbow Uniform", True)
    l.RescaleTransferFunction(pin[0], pin[1])
    if logscale:
        l.MapControlPointsToLogSpace()
        l.UseLogScale = 1
    bar(l, "|B|  T")
    dd.SetScalarBarVisibility(view, True)
    look(scale, along="y")
    shot(name)


b_on_plane("b_port_domain.png", 1.03 * FAR, (1.0e-8, 4.0e-4))
# Zoomed onto the PORT cross-section. LOG SCALE, because this view spans more
# than two decades: the conductor sits at J/sigma = 8.6e-3 while the air right
# at the cut reaches 1.7 V/m. A linear range that shows the conductor saturates
# every air cell around it into one flat colour, which is how this first came
# out.
e_on_plane("e_port_zoom.png", 0.0, (R, 0, 0), 3.0 * a,
           (5.0e-3, 2.0), "|E|  V/m", logscale=True)
# The same zoom one tube radius to the +y side, OFF the cut, where Phi is
# single valued. The difference between the two is the discontinuity.
e_on_plane("e_port_zoom_off.png", a, (R, a, 0), 3.0 * a,
           (5.0e-3, 2.0), "|E|  V/m", logscale=True)

print("done")
