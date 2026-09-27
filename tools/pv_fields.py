from paraview.simple import *
import os, sys

HERE = r"C:\Research\APhi_Solver_Project_LowFrequency_EDA\APhi_Solver"
OUT  = r"C:\Users\nasta\AppData\Local\Temp\claude\C--Research-APhi-Solver-Project-LowFrequency-EDA\a5512f13-3d9d-427a-9ffb-314df45d376f\scratchpad"

paraview.simple._DisableFirstRenderCameraReset()
src = LegacyVTKReader(registrationName='potential.vtk',
                      FileNames=[os.path.join(HERE, 'potential.vtk')])

wire = Threshold(registrationName='wire', Input=src)
wire.Scalars = ['CELLS', 'body_tag']
wire.LowerThreshold = 1.0; wire.UpperThreshold = 1.0; wire.ThresholdMethod = 'Between'

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [900, 1000]
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]; view.Background2 = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 1

FOCAL = [0.0, 0.0, 0.0005]

def show(source, array, assoc, title, fmt='{:.4g}'):
    for r in GetRepresentations().values():
        try: r.Visibility = 0
        except Exception: pass
    d = Show(source, view, 'UnstructuredGridRepresentation')
    d.Representation = 'Surface'
    ColorBy(d, (assoc, array))
    d.NonlinearSubdivisionLevel = 3
    d.Ambient = 1.0; d.Diffuse = 0.0; d.Specular = 0.0
    d.RescaleTransferFunctionToDataRange(True, False)
    d.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset('Viridis', True)
    bar = GetScalarBar(lut, view)
    bar.Title = title; bar.ComponentTitle = ''
    bar.WindowLocation = 'Any Location'; bar.Position = [0.80, 0.25]
    bar.ScalarBarLength = 0.5; bar.ScalarBarThickness = 18
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.AutomaticLabelFormat = 1; bar.AddRangeLabels = 1
    bar.RangeLabelFormat = fmt
    return d

view.CameraParallelProjection = 1
view.CameraParallelScale = 0.00062
view.CameraFocalPoint = FOCAL
view.CameraPosition = [0.0, -0.0040, 0.0005]
view.CameraViewUp = [0.0, 0.0, 1.0]

def shot(name):
    Render(view)
    p = os.path.join(OUT, name)
    SaveScreenshot(p, view, ImageResolution=view.ViewSize)
    print("wrote %-22s %7d bytes" % (name, os.path.getsize(p)))

# 1. The nodal E magnitude, unfiltered -- what a naive look shows.
show(wire, 'E_magnitude', 'POINTS', '|E|  [V/m]')
shot('field_E_raw.png')

# 2. The same, with the interface nodes removed.
clean = Threshold(registrationName='interior', Input=wire)
clean.Scalars = ['POINTS', 'material_interface']
clean.LowerThreshold = 0.0; clean.UpperThreshold = 0.0; clean.ThresholdMethod = 'Between'
show(clean, 'E_magnitude', 'POINTS', '|E|  [V/m]')
shot('field_E_clean.png')

# 3. J, per cell, which never straddles an interface.
show(wire, 'J_real', 'CELLS', 'Re(J)  [A/m2]')
shot('field_J.png')

SaveState(os.path.join(HERE, 'potential_fields.pvsm'))
print("wrote state")
