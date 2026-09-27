from paraview.simple import *
from paraview import servermanager as sm
import math, os

HERE = r"C:\Research\APhi_Solver_Project_LowFrequency_EDA\APhi_Solver"
src = LegacyVTKReader(registrationName='p', FileNames=[os.path.join(HERE, 'potential.vtk')])
cc = CellCenters(registrationName='cc', Input=src)
cc.VertexCells = 1
d = sm.Fetch(cc)
pts = d.GetPoints()
B = d.GetPointData().GetArray('B_cell_real')
n = d.GetNumberOfPoints()

I, a, MU0 = 7205.5168, 2.0e-4, 4e-7 * math.pi
# The meshed conductor is a 24-gon, so the enclosed-current fraction inside
# radius r is (area of the 24-gon clipped to r) / (its full area). Near the axis
# that is just pi r^2 / A_poly until r reaches the apothem.
apothem = a * math.cos(math.pi / 24.0)
A_poly = 0.5 * 24 * a * a * math.sin(2.0 * math.pi / 24.0)

rows = {}
for k in range(n):
    x, y, z = pts.GetPoint(k)
    if z < 4.0e-4 or z > 6.0e-4:      # near mid height, away from the ends
        continue
    r = math.hypot(x, y)
    bx, by, bz = B.GetTuple(k)
    mag = math.sqrt(bx * bx + by * by + bz * bz)
    # B should be purely azimuthal: e_phi = (-y, x, 0)/r
    if r > 1e-9:
        b_phi = (-y * bx + x * by) / r
        b_rad = (x * bx + y * by) / r
    else:
        b_phi, b_rad = mag, 0.0
    rows.setdefault(int(r / 2.5e-5), []).append((r, mag, b_phi, b_rad, bz))

print("   r range (um)    n     |B| meas      |B| Ampere      ratio    |B_r|/|B|   |B_z|/|B|")
for kb in sorted(rows):
    v = rows[kb]
    rm = sum(t[0] for t in v) / len(v)
    mm = sum(t[1] for t in v) / len(v)
    if rm <= apothem:
        pred = MU0 * I * (math.pi * rm * rm / A_poly) / (2 * math.pi * rm)
    else:
        pred = MU0 * I / (2 * math.pi * rm)
    rad = sum(abs(t[3]) for t in v) / len(v)
    axi = sum(abs(t[4]) for t in v) / len(v)
    print("  %5.0f - %5.0f  %5d  %11.4f  %13.4f  %9.4f  %10.2e  %10.2e"
          % (kb * 25, kb * 25 + 25, len(v), mm, pred, mm / pred if pred else 0,
             rad / mm if mm else 0, axi / mm if mm else 0))
