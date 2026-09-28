from paraview.simple import *
import os

HERE = r"C:\Research\APhi_Solver_Project_LowFrequency_EDA\APhi_Solver"
OUT  = r"C:\Users\nasta\AppData\Local\Temp\claude\C--Research-APhi-Solver-Project-LowFrequency-EDA\a5512f13-3d9d-427a-9ffb-314df45d376f\scratchpad"
ZMID = 5.0e-4

paraview.simple._DisableFirstRenderCameraReset()

# One file per field, so B and E come from two readers. They share a mesh, so
# the two slices are geometrically identical and the views line up.
b_src = LegacyVTKReader(registrationName='B_field.vtk',
                        FileNames=[os.path.join(HERE, 'B_field.vtk')])
e_src = LegacyVTKReader(registrationName='E_field.vtk',
                        FileNames=[os.path.join(HERE, 'E_field.vtk')])


def cut(src, tag):
    """Horizontal cut at mid height, and the same cut restricted to the copper."""
    whole = Slice(registrationName='mid plane ' + tag, Input=src)
    whole.SliceType = 'Plane'
    whole.SliceType.Origin = [0.0, 0.0, ZMID]
    whole.SliceType.Normal = [0.0, 0.0, 1.0]

    wire3d = Threshold(registrationName='wire ' + tag, Input=src)
    wire3d.Scalars = ['CELLS', 'body_tag']
    wire3d.LowerThreshold = 1.0
    wire3d.UpperThreshold = 1.0
    wire3d.ThresholdMethod = 'Between'
    inner = Slice(registrationName='wire mid plane ' + tag, Input=wire3d)
    inner.SliceType = 'Plane'
    inner.SliceType.Origin = [0.0, 0.0, ZMID]
    inner.SliceType.Normal = [0.0, 0.0, 1.0]
    return whole, inner


b_whole, b_wire = cut(b_src, 'B')
e_whole, e_wire = cut(e_src, 'E')

# The B plots below use the B reader, the E plots the E reader.
whole, wire = b_whole, b_wire

view = GetActiveViewOrCreate('RenderView')
view.ViewSize = [980, 900]
view.UseColorPaletteForBackground = 0
view.Background = [1.0, 1.0, 1.0]; view.Background2 = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 1
view.CameraParallelProjection = 1

_shown = []

def clear():
    # Hiding a representation does NOT hide its scalar bar: a stale legend
    # stays drawn and the next one lands on top of it.
    for d in _shown:
        try:
            d.SetScalarBarVisibility(view, False)
            d.Visibility = 0
        except Exception:
            pass
    del _shown[:]

def draw(source, array, assoc, title, comp=None, log=False):
    clear()
    d = Show(source, view, 'GeometryRepresentation')
    d.Representation = 'Surface'
    if comp: ColorBy(d, (assoc, array, comp))
    else:    ColorBy(d, (assoc, array))
    d.NonlinearSubdivisionLevel = 3
    d.Ambient = 1.0; d.Diffuse = 0.0; d.Specular = 0.0
    d.RescaleTransferFunctionToDataRange(True, False)
    d.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset('Jet', True)          # blue low -> red high
    if log:
        lut.MapControlPointsToLogSpace(); lut.UseLogScale = 1
    else:
        lut.UseLogScale = 0
    bar = GetScalarBar(lut, view)
    bar.Title = title; bar.ComponentTitle = ''
    bar.WindowLocation = 'Any Location'; bar.Position = [0.84, 0.24]
    bar.ScalarBarLength = 0.52; bar.ScalarBarThickness = 18
    bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]
    bar.AutomaticLabelFormat = 1; bar.AddRangeLabels = 1
    _shown.append(d)
    return d

def look(half):
    view.CameraParallelScale = half
    view.CameraFocalPoint = [0.0, 0.0, ZMID]
    view.CameraPosition = [0.0, 0.0, ZMID + 0.01]
    view.CameraViewUp = [0.0, 1.0, 0.0]

def shot(name):
    Render(view)
    p = os.path.join(OUT, name)
    SaveScreenshot(p, view, ImageResolution=view.ViewSize)
    print("wrote %-24s %7d bytes" % (name, os.path.getsize(p)))

# ---- whole domain, 2 x 2 mm -------------------------------------------------
look(1.15e-3)
draw(whole, 'B_magnitude', 'POINTS', '|B|  [T]', log=True)
shot('slice_B_domain.png')

draw(e_whole, 'E_cell_real', 'CELLS', 'Re(E)  [V/m]', comp='Magnitude', log=True)
shot('slice_E_domain.png')

# ---- the wire alone, 0.2 mm radius ------------------------------------------
look(2.6e-4)
draw(wire, 'B_magnitude', 'POINTS', '|B|  [T]')
shot('slice_B_wire.png')

draw(e_wire, 'E_cell_real', 'CELLS', 'Re(E)  [V/m]', comp='Magnitude')
shot('slice_E_wire.png')

# ---- B as vectors, to show the circulation ----------------------------------
look(1.15e-3)
clear()
d = Show(whole, view, 'GeometryRepresentation')
d.Representation = 'Surface'
ColorBy(d, ('POINTS', 'B_magnitude'))
d.Ambient = 1.0; d.Diffuse = 0.0; d.Specular = 0.0
d.RescaleTransferFunctionToDataRange(True, False)
lut = GetColorTransferFunction('B_magnitude')
lut.ApplyPreset('Jet', True); lut.MapControlPointsToLogSpace(); lut.UseLogScale = 1
d.SetScalarBarVisibility(view, True)
bar = GetScalarBar(lut, view)
bar.Title = '|B|  [T]'; bar.ComponentTitle = ''
bar.WindowLocation = 'Any Location'; bar.Position = [0.84, 0.24]
bar.ScalarBarLength = 0.52; bar.ScalarBarThickness = 18
bar.TitleColor = [0, 0, 0]; bar.LabelColor = [0, 0, 0]

g = Glyph(registrationName='Bvec', Input=whole, GlyphType='Arrow')
g.OrientationArray = ['POINTS', 'B_real']
g.ScaleArray = ['POINTS', 'No scale array']
g.ScaleFactor = 1.1e-4
g.GlyphMode = 'Every Nth Point'
g.Stride = 7
_shown.append(d)
gd = Show(g, view, 'GeometryRepresentation')
_shown.append(gd)
gd.AmbientColor = [0, 0, 0]; gd.DiffuseColor = [0, 0, 0]
ColorBy(gd, None)
shot('slice_B_vectors.png')

SaveState(os.path.join(HERE, 'potential_slices.pvsm'))
print("wrote state")
