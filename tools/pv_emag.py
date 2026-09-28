from paraview.simple import *
import os

HERE = r"C:\Research\APhi_Solver_Project_LowFrequency_EDA\APhi_Solver"
OUT  = r"C:\Users\nasta\AppData\Local\Temp\claude\C--Research-APhi-Solver-Project-LowFrequency-EDA\a5512f13-3d9d-427a-9ffb-314df45d376f\scratchpad"
ZMID = 5.0e-4

paraview.simple._DisableFirstRenderCameraReset()
src = LegacyVTKReader(registrationName='E_field.vtk',
                      FileNames=[os.path.join(HERE, 'E_field.vtk')])
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

def draw(array, assoc, comp, lo, hi, log, title, name):
    clear()
    d = Show(sl, view, 'GeometryRepresentation')
    d.Representation = 'Surface'
    ColorBy(d, (assoc, array, comp) if comp else (assoc, array))
    d.NonlinearSubdivisionLevel = 3
    d.Ambient = 1.0; d.Diffuse = 0.0; d.Specular = 0.0
    d.SetScalarBarVisibility(view, True)
    shown.append(d)
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset('Jet', True)
    lut.UseLogScale = 0
    if log:
        lut.RescaleTransferFunction(lo, hi); lut.MapControlPointsToLogSpace(); lut.UseLogScale = 1
    else:
        lut.MapControlPointsToLinearSpace(); lut.RescaleTransferFunction(lo, hi)
    bar = GetScalarBar(lut, view)
    bar.Title = title; bar.ComponentTitle = ''
    bar.WindowLocation = 'Any Location'; bar.Position = [0.855, 0.22]
    bar.ScalarBarLength = 0.56; bar.ScalarBarThickness = 20
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.AutomaticLabelFormat = 1; bar.AddRangeLabels = 1
    Render(view)
    p = os.path.join(OUT, name)
    SaveScreenshot(p, view, ImageResolution=view.ViewSize)
    print("wrote %-26s %7d bytes" % (name, os.path.getsize(p)))

sl.UpdatePipeline()
rc = sl.CellData['E_cell_real'].GetRange(-1)     # -1 = magnitude
rn = sl.PointData['E_magnitude'].GetRange()
print("  |Re(E)| per cell  : %10.3f .. %10.3f V/m" % rc)
print("  |E| nodal         : %10.3f .. %10.3f V/m   <- inflated at the interface" % rn)

draw('E_cell_real', 'CELLS', 'Magnitude', 0.0, rc[1], False, '|E|  [V/m]', 'emag_linear.png')
draw('E_cell_real', 'CELLS', 'Magnitude', rc[0], rc[1], True, '|E|  [V/m]', 'emag_log.png')
# The nodal version, for the contrast: same plane, same scale type.
draw('E_magnitude', 'POINTS', None, 0.0, rn[1], False, '|E| nodal  [V/m]', 'emag_nodal.png')
