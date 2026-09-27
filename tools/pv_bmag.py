from paraview.simple import *
from paraview import servermanager as sm
import math, os

HERE = r"C:\Research\APhi_Solver_Project_LowFrequency_EDA\APhi_Solver"
OUT  = r"C:\Users\nasta\AppData\Local\Temp\claude\C--Research-APhi-Solver-Project-LowFrequency-EDA\a5512f13-3d9d-427a-9ffb-314df45d376f\scratchpad"
ZMID = 5.0e-4

paraview.simple._DisableFirstRenderCameraReset()
src = LegacyVTKReader(registrationName='potential.vtk',
                      FileNames=[os.path.join(HERE, 'potential.vtk')])
sl = Slice(registrationName='mid', Input=src)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.0, 0.0, ZMID]
sl.SliceType.Normal = [0.0, 0.0, 1.0]

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [980, 900]
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]; view.Background2 = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.CameraParallelScale = 1.05e-3
view.CameraFocalPoint = [0.0, 0.0, ZMID]
view.CameraPosition = [0.0, 0.0, ZMID + 0.01]
view.CameraViewUp = [0.0, 1.0, 0.0]

shown = []
def clear():
    for d in shown:
        try:
            d.SetScalarBarVisibility(view, False); d.Visibility = 0
        except Exception: pass
    del shown[:]

def draw(lo, hi, log, name):
    clear()
    d = Show(sl, view, 'GeometryRepresentation')
    d.Representation = 'Surface'
    ColorBy(d, ('POINTS', 'B_magnitude'))
    d.NonlinearSubdivisionLevel = 3
    d.Ambient = 1.0; d.Diffuse = 0.0; d.Specular = 0.0
    d.SetScalarBarVisibility(view, True)
    shown.append(d)
    lut = GetColorTransferFunction('B_magnitude')
    lut.ApplyPreset('Jet', True)
    lut.UseLogScale = 0
    if log:
        lut.RescaleTransferFunction(lo, hi); lut.MapControlPointsToLogSpace(); lut.UseLogScale = 1
    else:
        lut.MapControlPointsToLinearSpace(); lut.RescaleTransferFunction(lo, hi)
    bar = GetScalarBar(lut, view)
    bar.Title = '|B|  [T]'; bar.ComponentTitle = ''
    bar.WindowLocation = 'Any Location'; bar.Position = [0.855, 0.22]
    bar.ScalarBarLength = 0.56; bar.ScalarBarThickness = 20
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.AutomaticLabelFormat = 1; bar.AddRangeLabels = 1
    Render(view)
    p = os.path.join(OUT, name)
    SaveScreenshot(p, view, ImageResolution=view.ViewSize)
    print("wrote %-26s %7d bytes" % (name, os.path.getsize(p)))

d = Show(sl, view); sl.UpdatePipeline()
rng = sl.PointData['B_magnitude'].GetRange()
print("  |B| on the mid plane: %.4f .. %.4f T" % rng)
Hide(sl, view)

draw(0.0, rng[1], False, 'bmag_linear.png')
draw(rng[0], rng[1], True, 'bmag_log.png')
