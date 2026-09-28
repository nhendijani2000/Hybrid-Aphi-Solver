"""Per-cell |B| against nodal |B|, same framing, for case 02.

    "C:\Program Files\ParaView 6.1.1\bin\pvbatch.exe" pv_bcell.py

B = curl A is CONSTANT over each tetrahedron with first-order Whitney edge
elements, so B_cell_* is the field the solver actually produced -- nothing
combined, nothing smoothed. The nodal B_* arrays are the lumped-mass L2
projection of that onto the nodes (a volume-weighted average over each node's
element patch). These figures put the two side by side.

Every pair is rendered twice:
  * banded   (11 colours, what the report uses)
  * continuous (256 colours, which hides nothing)
and the zoomed pair is PINNED to a common range so the colours mean the same
thing in both images.
"""

from paraview.simple import *
import os
import sys

HERE = sys.argv[1] if len(sys.argv) > 1 else "output"
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "plots")

ZMID = 20.0e-3
A = 1.5e-3
W = 40.0e-3

if not os.path.isdir(OUT):
    os.makedirs(OUT)

paraview.simple._DisableFirstRenderCameraReset()

p = os.path.join(HERE, "B_field.vtk")
if not os.path.exists(p):
    raise SystemExit("no such file: %s\nsolve the case first:\n"
                     "    ..\run_case.bat cylinder_50hz.aphi" % p)
b_src = LegacyVTKReader(registrationName="B", FileNames=[p])

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1000, 900]
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.Background2 = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1

_shown = []


def clear():
    for d in _shown:
        try:
            d.SetScalarBarVisibility(view, False)
            d.Visibility = 0
        except Exception:
            pass
    del _shown[:]


def look(half):
    view.CameraParallelScale = half
    view.CameraFocalPoint = [0.0, 0.0, ZMID]
    view.CameraPosition = [0.0, 0.0, ZMID + 1.0]
    view.CameraViewUp = [0.0, 1.0, 0.0]


def shot(name):
    Render(view)
    q = os.path.join(OUT, name)
    SaveScreenshot(q, view, ImageResolution=view.ViewSize)
    print("  wrote %-34s %8d bytes" % (name, os.path.getsize(q)))


# The slice, once. |B_cell| is the phasor amplitude sqrt(|Re|^2 + |Im|^2),
# matching how the solver builds the nodal B_magnitude scalar -- otherwise the
# comparison would be |Re(B)| against |B| and the cell picture would look
# wrongly small.
sl = Slice(Input=b_src)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0, 0, ZMID]
sl.SliceType.Normal = [0, 0, 1]

calc = Calculator(Input=sl)
calc.AttributeType = "Cell Data"
calc.ResultArrayName = "B_cell_magnitude"
calc.Function = "sqrt(mag(B_cell_real)^2 + mag(B_cell_imag)^2)"

# Zoom region, so a range pinned for the close-up reflects what is on screen
# rather than the whole 40 mm box.
box = Clip(Input=calc)
box.ClipType = "Box"
box.ClipType.Position = [-3 * A, -3 * A, ZMID - 1e-6]
box.ClipType.Length = [6 * A, 6 * A, 2e-6]
box.Invert = 1


def legend(array, title, bands):
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset("Jet", True)
    lut.Discretize = 1
    lut.NumberOfTableValues = bands
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
    bar.RangeLabelFormat = "{:.3g}"
    return lut


def draw(src, array, assoc, title, bands, rng=None):
    d = Show(src, view, "GeometryRepresentation")
    d.Representation = "Surface"
    ColorBy(d, (assoc, array))
    lut = legend(array, title, bands)
    if rng is None:
        d.RescaleTransferFunctionToDataRange(True, False)
    else:
        lut.RescaleTransferFunction(rng[0], rng[1])
    lut.MapControlPointsToLinearSpace()
    d.SetScalarBarVisibility(view, True)
    d.NonlinearSubdivisionLevel = 3
    d.Ambient = 1.0
    d.Diffuse = 0.0
    d.Specular = 0.0
    _shown.append(d)
    return d


def rng_of(src, array, assoc):
    UpdatePipeline(proxy=src)
    info = src.GetCellDataInformation() if assoc == "CELLS" \
        else src.GetPointDataInformation()
    return info[array].GetComponentRange(0)


# Common range for the zoomed pair: the union of the two, so neither is
# clipped and the colours are directly comparable.
r_cell = rng_of(box, "B_cell_magnitude", "CELLS")
r_node = rng_of(box, "B_magnitude", "POINTS")
pin = (min(r_cell[0], r_node[0]), max(r_cell[1], r_node[1]))
print("  zoom range  cell %.4g .. %.4g" % (r_cell[0], r_cell[1]))
print("  zoom range  node %.4g .. %.4g" % (r_node[0], r_node[1]))
print("  pinned      %.4g .. %.4g" % (pin[0], pin[1]))

for bands, tag in ((11, "banded"), (256, "continuous")):
    # --- zoom, pinned range, the pair that answers the question -------------
    clear()
    look(3.0 * A)
    draw(box, "B_cell_magnitude", "CELLS", "|B| per cell (T)", bands, pin)
    shot("08_B_cell_zoom_%s.png" % tag)

    clear()
    look(3.0 * A)
    draw(box, "B_magnitude", "POINTS", "|B| nodal (T)", bands, pin)
    shot("08b_B_nodal_zoom_%s.png" % tag)

    # --- whole domain, each to its own range --------------------------------
    clear()
    look(0.5 * W * 1.05)
    draw(calc, "B_cell_magnitude", "CELLS", "|B| per cell (T)", bands)
    shot("08c_B_cell_domain_%s.png" % tag)

    clear()
    look(0.5 * W * 1.05)
    draw(calc, "B_magnitude", "POINTS", "|B| nodal (T)", bands)
    shot("08d_B_nodal_domain_%s.png" % tag)

print("done -> %s" % OUT)
