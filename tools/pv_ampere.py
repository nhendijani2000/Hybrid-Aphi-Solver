from paraview.simple import *
from paraview import servermanager as sm
import math, os, sys

HERE = sys.argv[1] if len(sys.argv) > 1 else "."
A_MM = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
I_AMP = float(sys.argv[3]) if len(sys.argv) > 3 else 144505.189
NGON = int(sys.argv[4]) if len(sys.argv) > 4 else 36
ZLO = float(sys.argv[5]) if len(sys.argv) > 5 else 16.0e-3
ZHI = float(sys.argv[6]) if len(sys.argv) > 6 else 24.0e-3
src = LegacyVTKReader(registrationName="B", FileNames=[os.path.join(HERE, "B_field.vtk")])
cc = CellCenters(registrationName='cc', Input=src)
cc.VertexCells = 1
d = sm.Fetch(cc)
pts = d.GetPoints()
BR = d.GetPointData().GetArray('B_cell_real')
BI = d.GetPointData().GetArray('B_cell_imag')
n = d.GetNumberOfPoints()

# The PHASOR magnitude, sqrt(|Bx|^2 + |By|^2 + |Bz|^2), not the length of the
# real part. At 50 Hz on a 10 mm conductor the current is 47270 - 136555j A:
# Re(I) is 0.327 of |I|, so comparing |Re(B)| against mu0*|I|/(2 pi r) reports a
# ratio of 0.33 everywhere and looks like the field is three times too small.
# It is not; the two sides were simply different quantities.

I, a, MU0 = I_AMP, A_MM * 1e-3, 4e-7 * math.pi
# The meshed conductor is a 24-gon, so the enclosed-current fraction inside
# radius r is (area of the 24-gon clipped to r) / (its full area). Near the axis
# that is just pi r^2 / A_poly until r reaches the apothem.
apothem = a * math.cos(math.pi / NGON)
A_poly = 0.5 * NGON * a * a * math.sin(2.0 * math.pi / NGON)

rows = {}
for k in range(n):
    x, y, z = pts.GetPoint(k)
    if z < ZLO or z > ZHI:            # near mid height, away from the ends
        continue
    r = math.hypot(x, y)
    rx, ry, rz = BR.GetTuple(k)
    ix, iy, iz = BI.GetTuple(k)
    mag = math.sqrt(rx * rx + ix * ix + ry * ry + iy * iy + rz * rz + iz * iz)
    # B should be purely azimuthal: e_phi = (-y, x, 0)/r. Both parts projected,
    # then recombined, so these are phasor magnitudes like `mag` above.
    if r > 1e-9:
        b_phi = math.hypot((-y * rx + x * ry) / r, (-y * ix + x * iy) / r)
        b_rad = math.hypot((x * rx + y * ry) / r, (x * ix + y * iy) / r)
    else:
        b_phi, b_rad = mag, 0.0
    b_ax = math.hypot(rz, iz)
    rows.setdefault(int(r / (a / 8.0)), []).append((r, mag, b_phi, b_rad, b_ax))

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
