"""Render the conductor's potential from potential.vtk, and save a ParaView state.

Run it with ParaView's own interpreter, not the system Python:

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" tools/pv_wire.py
    pvbatch tools/pv_wire.py  potential.vtk  out/            # explicit paths

It writes wire_3d.png, wire_side.png and wire_cut.png, plus
potential_wire.pvsm, which `File -> Load State` reopens with everything already
set up. `docs/POSTPROCESSING_PLAN.md` section 9 is what the pictures show.

Three things here are not cosmetic, and each was a bug first:

  * NonlinearSubdivisionLevel. The cells are VTK_QUADRATIC_TETRA. At the
    default of 1 ParaView draws each as if it were linear, discarding the P2
    mid-edge values the solve actually determined -- the picture looks fine and
    is the linear interpolant.

  * Flat shading on the field views (Ambient 1, Diffuse 0). With lighting on,
    the edge darkening reads as a lower potential. It is not.

  * Re-asserting the scalar bar after the final Show(). ParaView's
    "automatically hide unused scalar bars" turns the legend off when the clip
    is hidden, and Show() alone does not bring it back, so the state saves with
    Visibility = 0 and opens with no legend.

The colour range is pinned to [0, 1] rather than rescaled to the data, because
0 V and 1 V are the prescribed terminal voltages: the colours then mean the
same thing across runs. A case whose potential leaves that range will saturate
at the ends -- rescale to the data range to see it.
"""

import os
import sys

from paraview.simple import *  # noqa: F403

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VTK = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "potential.vtk")
OUT = sys.argv[2] if len(sys.argv) > 2 else HERE
STATE = os.path.join(HERE, "potential_wire.pvsm")

WIRE_TAG = 1  # body_tag: 1 is the copper wire, 2 is the air box

if not os.path.exists(VTK):
    raise SystemExit("no such file: %s\nrun solve_mesh on a case first" % VTK)
if not os.path.isdir(OUT):
    os.makedirs(OUT)

paraview.simple._DisableFirstRenderCameraReset()  # noqa: F405
src = LegacyVTKReader(registrationName="potential.vtk", FileNames=[VTK])  # noqa: F405

wire = Threshold(registrationName="wire only", Input=src)  # noqa: F405
wire.Scalars = ["CELLS", "body_tag"]
wire.LowerThreshold = float(WIRE_TAG)
wire.UpperThreshold = float(WIRE_TAG)
wire.ThresholdMethod = "Between"
wire.UpdatePipeline()

info = wire.GetDataInformation()
print("wire cells %d   points %d" % (info.GetNumberOfCells(), info.GetNumberOfPoints()))
bounds = info.GetBounds()
print("bounds x[%.4g %.4g] y[%.4g %.4g] z[%.4g %.4g]" % bounds)
if info.GetNumberOfCells() == 0:
    raise SystemExit("threshold selected nothing -- is body_tag %d the conductor?" % WIRE_TAG)

# Centre of the conductor, so the cameras below do not depend on the mesh.
FOCAL = [
    0.5 * (bounds[0] + bounds[1]),
    0.5 * (bounds[2] + bounds[3]),
    0.5 * (bounds[4] + bounds[5]),
]
span = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])

view = GetActiveViewOrCreate("RenderView")  # noqa: F405
view.ViewSize = [900, 1000]
view.OrientationAxesVisibility = 1
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]
view.Background2 = [1.0, 1.0, 1.0]

disp = Show(wire, view, "UnstructuredGridRepresentation")  # noqa: F405
disp.Representation = "Surface"
ColorBy(disp, ("POINTS", "phi_real"))  # noqa: F405
disp.NonlinearSubdivisionLevel = 3
disp.SetScalarBarVisibility(view, True)

lut = GetColorTransferFunction("phi_real")  # noqa: F405
lut.ApplyPreset("Viridis", True)
lut.RescaleTransferFunction(0.0, 1.0)

bar = GetScalarBar(lut, view)  # noqa: F405
bar.Title = "Re(Phi)"
bar.ComponentTitle = "[V]"
bar.WindowLocation = "Any Location"
bar.Position = [0.855, 0.25]
bar.ScalarBarLength = 0.5
bar.ScalarBarThickness = 18
bar.TitleColor = [0.0, 0.0, 0.0]
bar.LabelColor = [0.0, 0.0, 0.0]
bar.AutomaticLabelFormat = 1
bar.AddRangeLabels = 1
# std::format style: a printf format works but is deprecated in ParaView 6.1.
bar.RangeLabelFormat = '{:.2f}'


def shot(name):
    Render(view)  # noqa: F405
    path = os.path.join(OUT, name)
    SaveScreenshot(path, view, ImageResolution=view.ViewSize, TransparentBackground=0)  # noqa: F405
    print("wrote %-16s %7d bytes" % (name, os.path.getsize(path)))


# 1. From outside, off-axis so the view-up is never parallel to the view normal
#    (ParaView silently resets the view-up when it is, and the result is tipped).
view.CameraParallelProjection = 0
view.CameraFocalPoint = FOCAL
view.CameraPosition = [FOCAL[0] + 2.6 * span, FOCAL[1] - 2.4 * span, FOCAL[2] + 1.2 * span]
view.CameraViewUp = [0.0, 0.0, 1.0]
shot("wire_3d.png")

# 2. Straight on, parallel projection, no lighting: the one view in which the
#    colour is the value and a linear gradient looks linear.
disp.Ambient = 1.0
disp.Diffuse = 0.0
disp.Specular = 0.0
view.CameraParallelProjection = 1
view.CameraParallelScale = 0.62 * span
view.CameraFocalPoint = FOCAL
view.CameraPosition = [FOCAL[0], FOCAL[1] - 4.0 * span, FOCAL[2]]
view.CameraViewUp = [0.0, 0.0, 1.0]
shot("wire_side.png")

# 3. Cut open on y = 0. Viewed straight on this is pixel-for-pixel identical to
#    the side view, because with no skin effect Phi has no radial variation at
#    all -- which is the result, but makes an unreadable picture. Oblique, so
#    the cut face is visible and the interior is demonstrably what is shown.
clip = Clip(registrationName="half", Input=wire)  # noqa: F405
clip.ClipType = "Plane"
clip.ClipType.Origin = FOCAL
clip.ClipType.Normal = [0.0, 1.0, 0.0]
Hide(wire, view)  # noqa: F405
dclip = Show(clip, view, "UnstructuredGridRepresentation")  # noqa: F405
dclip.Representation = "Surface"
ColorBy(dclip, ("POINTS", "phi_real"))  # noqa: F405
dclip.NonlinearSubdivisionLevel = 3
dclip.SetScalarBarVisibility(view, True)
dclip.Ambient = 1.0
dclip.Diffuse = 0.0
dclip.Specular = 0.0
view.CameraParallelProjection = 0
view.CameraFocalPoint = FOCAL
view.CameraPosition = [FOCAL[0] + 2.4 * span, FOCAL[1] + 2.6 * span, FOCAL[2] + 1.1 * span]
view.CameraViewUp = [0.0, 0.0, 1.0]
shot("wire_cut.png")

# Leave the state showing the solid conductor.
Hide(clip, view)  # noqa: F405
Show(wire, view)  # noqa: F405
# Hiding the clip let ParaView auto-hide the legend; Show() does not bring it
# back, and without this the state opens with no colour bar.
GetDisplayProperties(wire, view).SetScalarBarVisibility(view, True)  # noqa: F405
view.CameraParallelProjection = 0
view.CameraPosition = [FOCAL[0] + 2.6 * span, FOCAL[1] - 2.4 * span, FOCAL[2] + 1.2 * span]
view.CameraViewUp = [0.0, 0.0, 1.0]
Render(view)  # noqa: F405
SaveState(STATE)  # noqa: F405
print("wrote state", STATE)
