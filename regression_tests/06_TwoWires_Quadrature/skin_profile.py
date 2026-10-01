"""The radial |J| profile against Kelvin, for case 06.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" skin_profile.py [output-dir]

Prints the table README.md quotes. Named skin_profile.py and not verify_*.py on
purpose: ../check.bat globs verify_*.py and keeps the last match, so a second
file under that name would silently displace the real verifier.

BOTH SIDES MUST BE ANCHORED AT THE SAME RADIUS, and getting that wrong is what
this file exists to prevent repeating.

Kelvin gives Jz(r)/Jz(a) = J0(kr)/J0(ka), a profile normalised at the surface
r = a. A measurement cannot be normalised there: a cell band "at the surface" is
r > 0.92a, whose MEAN radius is 0.9567a. Dividing the measurement by J(0.9567a)
while dividing the theory by J(a) compares two profiles pinned at different
places, and since |J0(k*0.9567a)/J0(k*a)| = 0.9660, it inflates every measured
ratio by 1/0.9660 = +3.5 %.

That is an offset with no h in it, so refining the mesh cannot remove it. It
produced a convincing false story: errors of 7.2 % that "improved" to 3.9 % under
refinement and then stalled, with an apparent convergence order collapsing to
0.37 -- all of it the constant 3.5 % offset dominating a real error ten times
smaller. Two controls finally placed it: a SINGLE isolated wire showed the same
error, ruling out the neighbouring wire, and narrowing the axial sample band
changed nothing, ruling out end effects. An error that survives removing the
neighbour, ignores position and ignores mesh size is in the instrument.

Normalised consistently, the real numbers at 4.2 elements per delta are +0.39 %
and +0.06 degrees.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "output")

A, D, LZ, SIG, F = 2.0e-3, 10.0e-3, 20.0e-3, 5.8e7, 2500.0
MU0 = 4e-7 * math.pi
OMEGA = 2.0 * math.pi * F
DELTA = math.sqrt(2.0 / (OMEGA * MU0 * SIG))
K = (1.0 - 1.0j) / DELTA


def bessel_j0(z):
    """J0 of a complex argument, by its defining series."""
    term, total = 1.0 + 0j, 1.0 + 0j
    for m in range(1, 60):
        term *= -(z * z) / (4.0 * m * m)
        total += term
        if abs(term) < 1e-18 * abs(total):
            break
    return total


def profile(out_dir, wire_tag=1, x0=-D / 2):
    s = LegacyVTKReader(FileNames=[os.path.join(out_dir, "J_field.vtk")])
    sz = CellSize(Input=s)
    sz.ComputeVolume = 1
    cc = CellCenters(Input=sz)
    cc.VertexCells = 1
    UpdatePipeline(proxy=cc)
    d = sm.Fetch(cc)
    q = vtk_to_numpy(d.GetPoints().GetData())
    vol = np.abs(vtk_to_numpy(d.GetPointData().GetArray("Volume")))
    bt = vtk_to_numpy(d.GetPointData().GetArray("body_tag"))
    J = (vtk_to_numpy(d.GetPointData().GetArray("J_cell_real"))[:, 2]
         + 1j * vtk_to_numpy(d.GetPointData().GetArray("J_cell_imag"))[:, 2])
    h = (12.0 * vol / math.sqrt(2.0)) ** (1.0 / 3.0)
    r = np.hypot(q[:, 0] - x0, q[:, 1])
    mid = (q[:, 2] > 0.35 * LZ) & (q[:, 2] < 0.65 * LZ)
    sel = mid & (bt == wire_tag)
    surf = sel & (r > 0.92 * A)
    return q, r, J, h, sel, surf


q, r, J, h, sel, surf = profile(OUT)
rbar = float(r[surf].mean())
Ja = J[surf].mean()

print("  delta = %.4f mm,  a/delta = %.3f,  median element %.3f mm  (%.1f per delta)"
      % (DELTA * 1e3, A / DELTA, np.median(h[sel]) * 1e3, DELTA / np.median(h[sel])))
print("  reference band r > 0.92a has mean radius %.4f a -- BOTH sides are"
      % (rbar / A))
print("  normalised there, not at a. See the note at the top of this file.")
print()
print("    r/a         measured   Kelvin     err        lag       Kelvin    err")
for lo, hi in ((0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.92)):
    m = sel & (r >= lo * A) & (r < hi * A)
    if m.sum() < 20:
        continue
    Jm = J[m].mean()
    ref = bessel_j0(K * r[m].mean()) / bessel_j0(K * rbar)
    lag = math.degrees(np.angle(Jm / Ja))
    rlag = math.degrees(np.angle(ref))
    print("    %.2f-%.2f    %8.4f   %8.4f  %+6.2f %%  %+7.2f  %+7.2f  %+5.2f deg"
          % (lo, hi, abs(Jm) / abs(Ja), abs(ref),
             100 * (abs(Jm) / abs(Ja) / abs(ref) - 1), lag, rlag, lag - rlag))

# proximity: the neighbour redistributes current around the circumference. This
# is real physics and Kelvin does not contain it, which is why the comparison
# above uses the AZIMUTHAL MEAN at each radius -- that removes the dipole term.
face = sel & (r > 0.9 * A) & ((q[:, 0] - (-D / 2)) > 0)
away = sel & (r > 0.9 * A) & ((q[:, 0] - (-D / 2)) < 0)
print()
print("  proximity effect, wire 1: |J| facing the neighbour / away = %.3f"
      % (abs(J[face].mean()) / abs(J[away].mean())))
