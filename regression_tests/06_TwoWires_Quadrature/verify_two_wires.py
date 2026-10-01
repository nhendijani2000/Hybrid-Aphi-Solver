"""Physics checks for 06_TwoWires_Quadrature.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" verify_two_wires.py [output-dir]

Run by ../check.bat, which prefers a case-local verifier over the shared
../verify.py -- the rod verifier assumes one conductor and one drive, and this
case has two of each.

WHAT IS WORTH ASSERTING HERE. The case exists to demonstrate the [postprocess]
options, so the assertions are the ones that would break a demonstration:

  1. each wire really carries 1 A, and
  2. they really are a quarter cycle apart -- this is the only case in the suite
     that uses a port `phase_deg` at all, so nothing else would catch it
     silently doing nothing, and
  3. the field really does become CIRCULARLY polarised, because that is the
     whole reason for the geometry. The analytic two-wire superposition is
     unambiguous about this: with I1 = I and I2 = jI,

         b/a = 0      exactly, everywhere on the line of centres
         b/a = 1      exactly, at (0, +/- d/2)

     Both are asserted. The second also pins N/a = sqrt(2) there, which is the
     equality case of b <= |E(th)| <= a <= N <= sqrt(2) a and ties the three
     magnitudes of docs/ComplexVectorPhasorConcept.md to a measurement.

The first version of this check read b/a = 2.7e-05 where theory says 1.0, and
the fault was in tools/postprocess.py: it divided |PxQ| by a once, which gives
the semi-minor axis b, not the ratio. That is what this case is for.
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

A_W, D, LZ = 2.0e-3, 10.0e-3, 20.0e-3
MU0 = 4e-7 * math.pi

g_checks = 0
g_failures = 0


def check(name, value, target, tol, unit=""):
    global g_checks, g_failures
    g_checks += 1
    ok = abs(value - target) <= tol
    if not ok:
        g_failures += 1
    print("  %-34s %-15.7g vs %-15.7g tol %-10.4g %s%s"
          % (name, value, target, tol, "ok  " if ok else "FAIL", unit))


def at_least(name, value, floor, unit=""):
    global g_checks, g_failures
    g_checks += 1
    ok = value >= floor
    if not ok:
        g_failures += 1
    print("  %-34s %-15.7g >= %-15.7g %12s%s"
          % (name, value, floor, "ok  " if ok else "FAIL", unit))


def cells(fn):
    s = LegacyVTKReader(FileNames=[os.path.join(OUT, fn)])
    sz = CellSize(Input=s)
    sz.ComputeVolume = 1
    cc = CellCenters(Input=sz)
    cc.VertexCells = 1
    UpdatePipeline(proxy=cc)
    return sm.Fetch(cc)


print("verifying %s" % HERE)
print("  output   %s" % OUT)
print("  two wires, 1 A each, a quarter cycle apart\n")

# --- the currents, and the quarter cycle ------------------------------------
d = cells("J_field.vtk")
q = vtk_to_numpy(d.GetPoints().GetData())
vol = np.abs(vtk_to_numpy(d.GetPointData().GetArray("Volume")))
bt = vtk_to_numpy(d.GetPointData().GetArray("body_tag"))
J = (vtk_to_numpy(d.GetPointData().GetArray("J_cell_real"))[:, 2]
     + 1j * vtk_to_numpy(d.GetPointData().GetArray("J_cell_imag"))[:, 2])

# A slab at mid-length, away from the end caps. Integrating Jz over a slab and
# dividing by its thickness gives the current through it: sum(Jz dV)/dz.
lo, hi = 0.4 * LZ, 0.6 * LZ
slab = (q[:, 2] > lo) & (q[:, 2] < hi)
print("current and phase")
args = {}
for tag, name in ((1, "wire1"), (2, "wire2")):
    m = slab & (bt == tag)
    I = (J[m] * vol[m]).sum() / (hi - lo)
    args[tag] = math.degrees(np.angle(I))
    check("|I| in %s, A" % name, abs(I), 1.0, 2.0e-2)
dphi = (args[2] - args[1] + 540.0) % 360.0 - 180.0
check("phase(wire2) - phase(wire1), deg", dphi, 90.0, 1.0)

# --- the polarization ellipse, which is the point of the case ---------------
db = cells("B_field.vtk")
qb = vtk_to_numpy(db.GetPoints().GetData())
btb = vtk_to_numpy(db.GetPointData().GetArray("body_tag"))
P = vtk_to_numpy(db.GetPointData().GetArray("B_cell_real"))
Q = vtk_to_numpy(db.GetPointData().GetArray("B_cell_imag"))

N2 = (P ** 2).sum(1) + (Q ** 2).sum(1)
S = 0.5 * N2
C = 0.5 * ((P ** 2).sum(1) - (Q ** 2).sum(1))
Dd = (P * Q).sum(1)
a = np.sqrt(S + np.hypot(C, Dd))
cr = np.sqrt((np.cross(P, Q) ** 2).sum(1))
safe = np.where(a > 0, a, 1.0)
ar = np.where(a > 0, cr / (safe * safe), 0.0)       # b/a -- NOT cr/a, see above

# air is Physical Volume 3 here, not 2 -- this case has THREE bodies
# (wire1 = 1, wire2 = 2, air = 3), unlike every other case in the suite.
air = btb == 3
mid = (qb[:, 2] > lo) & (qb[:, 2] < hi)
near = air & mid & (np.hypot(qb[:, 0], qb[:, 1]) < 0.9 * D)

print("\npolarization ellipse")
# The circular points sit at (0, +/- d/2). Take the cells nearest them.
spot = near & (np.abs(qb[:, 0]) < 0.15 * D) & (np.abs(np.abs(qb[:, 1]) - 0.5 * D) < 0.15 * D)
at_least("max b/a near (0, +/-d/2)", float(ar[spot].max()), 0.90)
# and N/a there must be its maximum, sqrt(2)
k = int(np.argmax(np.where(spot, ar, -1.0)))
check("N/a at the most circular cell", math.sqrt(N2[k]) / a[k], math.sqrt(2.0), 0.05)

# On the line joining the wires both contributions point the same way, so the
# sum is linear however they are phased: b/a = 0 EXACTLY, at y = 0. Away from
# the axis it grows, so a BAND does not have a zero median -- a first draft of
# this check asserted 0 and was wrong about that.
#
# A SECOND DRAFT WAS WRONG TOO, AND LESS OBVIOUSLY. It compared the cells'
# MEDIAN against an analytic value obtained by averaging over the window's AREA.
# Those are the same number only if the cells are spread uniformly over the
# window, and they are not: refine the mesh and the cell population shifts, so
# the statistic moved from 0.0658 to 0.0509 on a change that altered no physics.
# A test that moves when the mesh does, while the field does not, is measuring
# the mesh.
#
# So the analytic is evaluated AT EACH CELL'S OWN CENTROID and compared cell by
# cell. Whatever the mesh does, both sides are sampled identically.
axis = near & (np.abs(qb[:, 1]) < 0.08 * D) & (np.abs(qb[:, 0]) < 0.35 * D)


def analytic_axial_ratio(x, y):
    """b/a of two infinite wires at x = -+d/2 carrying I and jI."""
    r1 = np.stack([x + 0.5 * D, y], -1)
    r2 = np.stack([x - 0.5 * D, y], -1)
    s1 = (r1 ** 2).sum(-1)
    s2 = (r2 ** 2).sum(-1)
    k = MU0 / (2.0 * math.pi)
    Pa = k * np.stack([-r1[..., 1], r1[..., 0]], -1) / s1[..., None]
    Qa = k * np.stack([-r2[..., 1], r2[..., 0]], -1) / s2[..., None]
    n2 = (Pa ** 2).sum(-1) + (Qa ** 2).sum(-1)
    c = 0.5 * ((Pa ** 2).sum(-1) - (Qa ** 2).sum(-1))
    dd = (Pa * Qa).sum(-1)
    aa = np.sqrt(0.5 * n2 + np.hypot(c, dd))
    cra = np.abs(Pa[..., 0] * Qa[..., 1] - Pa[..., 1] * Qa[..., 0])
    return np.where(aa > 0, cra / (aa * aa), 0.0)


ref = analytic_axial_ratio(qb[axis, 0], qb[axis, 1])
resid = ar[axis] - ref
print("    (%d cells in the window, analytic median %.4f over the same cells)"
      % (axis.sum(), float(np.median(ref))))
check("median b/a - analytic, per cell", float(np.median(resid)), 0.0, 0.03)

print("\n%s" % ("all %d physics checks passed" % g_checks if g_failures == 0
                else "FAILED %d of %d checks" % (g_failures, g_checks)))
sys.exit(1 if g_failures else 0)
