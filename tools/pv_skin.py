from paraview.simple import *
from paraview import servermanager as sm
import math, os, sys

VTK  = sys.argv[1]
TAG  = sys.argv[2]
FREQ = float(sys.argv[3])
OUT  = r"C:\Users\nasta\AppData\Local\Temp\claude\C--Research-APhi-Solver-Project-LowFrequency-EDA\a5512f13-3d9d-427a-9ffb-314df45d376f\scratchpad"
ZMID = 5.0e-4

paraview.simple._DisableFirstRenderCameraReset()
src = LegacyVTKReader(registrationName='v', FileNames=[VTK])
wire3d = Threshold(registrationName='wire', Input=src)
wire3d.Scalars = ['CELLS', 'body_tag']
wire3d.LowerThreshold = 1.0; wire3d.UpperThreshold = 1.0; wire3d.ThresholdMethod = 'Between'
sl = Slice(registrationName='mid', Input=wire3d)
sl.SliceType = 'Plane'
sl.SliceType.Origin = [0.0, 0.0, ZMID]
sl.SliceType.Normal = [0.0, 0.0, 1.0]
sl.UpdatePipeline()

# |E| against radius, from the per-cell values of the whole conductor near mid height.
cc = CellCenters(registrationName='cc', Input=wire3d); cc.VertexCells = 1
d = sm.Fetch(cc)
pts = d.GetPoints(); ER = d.GetPointData().GetArray('E_cell_real')
EI = d.GetPointData().GetArray('E_cell_imag')
bins = {}
for k in range(d.GetNumberOfPoints()):
    x, y, z = pts.GetPoint(k)
    if z < 4.0e-4 or z > 6.0e-4: continue
    r = math.hypot(x, y)
    a3 = ER.GetTuple(k); b3 = EI.GetTuple(k)
    mag = math.sqrt(sum(a3[i]**2 + b3[i]**2 for i in range(3)))
    bins.setdefault(int(r / 2.5e-5), []).append(mag)

mu, sig = 4e-7 * math.pi, 5.8e7
delta = math.sqrt(2.0 / (2 * math.pi * FREQ * mu * sig))
print("  f = %g Hz   skin depth %.4f mm   a/delta = %.2f" % (FREQ, delta * 1e3, 2.0e-4 / delta))
ks = sorted(bins)
base = sum(bins[ks[0]]) / len(bins[ks[0]])
print("    r (um)     n      |E| V/m    relative to axis")
for kb in ks:
    v = bins[kb]; m = sum(v) / len(v)
    print("   %3d - %3d  %5d  %11.4g        %6.3f" % (kb*25, kb*25+25, len(v), m, m / base))

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [900, 860]
view.UseColorPaletteForBackground = 0
view.Background = [1,1,1]; view.Background2 = [1,1,1]
view.OrientationAxesVisibility = 0
view.CameraParallelProjection = 1
view.CameraParallelScale = 2.45e-4
view.CameraFocalPoint = [0,0,ZMID]; view.CameraPosition = [0,0,ZMID+0.01]
view.CameraViewUp = [0,1,0]

rep = Show(sl, view, 'GeometryRepresentation')
rep.Representation = 'Surface'
ColorBy(rep, ('CELLS', 'E_cell_real', 'Magnitude'))
rep.NonlinearSubdivisionLevel = 3
rep.Ambient = 1.0; rep.Diffuse = 0.0; rep.Specular = 0.0
rep.RescaleTransferFunctionToDataRange(True, False)
rep.SetScalarBarVisibility(view, True)
lut = GetColorTransferFunction('E_cell_real')
lut.ApplyPreset('Jet', True); lut.UseLogScale = 0
bar = GetScalarBar(lut, view)
bar.Title = '|E|  [V/m]'; bar.ComponentTitle = ''
bar.WindowLocation = 'Any Location'; bar.Position = [0.84, 0.22]
bar.ScalarBarLength = 0.56; bar.ScalarBarThickness = 20
bar.TitleColor = [0,0,0]; bar.LabelColor = [0,0,0]
bar.AutomaticLabelFormat = 1; bar.AddRangeLabels = 1
Render(view)
p = os.path.join(OUT, 'emag_wire_%s.png' % TAG)
SaveScreenshot(p, view, ImageResolution=view.ViewSize)
print("  wrote", os.path.basename(p))
