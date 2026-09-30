"""Assert the PHYSICS of the loop case.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" verify_loop.py [output-dir]

WHY THIS IS NOT ../verify.py. That one is built around a straight rod: it finds
the drive by differencing Phi between two END CAPS, and bins J by distance from
the z axis. A closed ring has no caps and no axis in that sense -- the current
circulates, and the terminal voltage lives in a jump ACROSS THE CUT. Contorting
the rod verifier to cover both would make neither clear, so check.bat prefers a
case-local verify_*.py where one exists.

WHAT IT CHECKS, and why each one is here:

 1. Current through the CUT equals the 1 A injected. The port's own promise.
 2. Current through theta = pi/2, where NOTHING is prescribed, equals it too.
    This is the check no previous case could make: a closed loop is multiply
    connected, the current has to circulate the whole way round, and the
    tree-cotree gauge has to permit that. A gauge that wrongly killed the loop
    would show up here and nowhere else.
 3. R against the EXACT dc resistance of this revolved polygon -- not the naive
    2 pi R/(sigma A). Current flows azimuthally so the path is 2 pi rho and the
    conductance is sigma * integral(dA / 2 pi rho): the current crowds to the
    inner radius. Ignoring that is a 0.44 % error, larger than the solver's own.
 4. L, PINNED to the measured value rather than asserted against theory --
    see the note in expected.txt. It is 4 % below the closed form and that is
    an open item, so this check catches CHANGE, it does not claim correctness.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import numpy as np
import math
import os
import sys

paraview.simple._DisableFirstRenderCameraReset()

CASE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(CASE, "output")

failures = []
checks = 0


def check(name, value, target, tol):
    global checks
    checks += 1
    ok = abs(value - target) <= tol
    print("  %-30s %-15.7g vs %-15.7g tol %-9.3g  %s (%.0f%% of tol)"
          % (name, value, target, tol, "ok  " if ok else "FAIL",
             100 * abs(value - target) / tol if tol > 0 else 0))
    if not ok:
        failures.append(name)


exp = {}
for line in open(os.path.join(CASE, "expected.txt")):
    line = line.split("#")[0].strip()
    if not line:
        continue
    p = line.split()
    vals = []
    for v in p[1:]:
        try:
            vals.append(float(v))
        except ValueError:
            vals.append(v)
    exp[p[0]] = vals

R = exp["loop_radius_mm"][0] * 1e-3
a = exp["wire_radius_mm"][0] * 1e-3
sigma = exp["sigma"][0]
freq = exp["frequency"][0]
omega = 2.0 * math.pi * freq

print("verifying %s" % CASE)
print("  output   %s" % OUT)
print("  driven by %g A through the cut" % exp["drive_current"][0])
print()

for f in ("J_field.vtk", "potential.vtk"):
    if not os.path.exists(os.path.join(OUT, f)):
        raise SystemExit("missing %s -- solve the case first" % f)


def ring(fn):
    s = LegacyVTKReader(FileNames=[os.path.join(OUT, fn)])
    t = Threshold(Input=s)
    t.Scalars = ["CELLS", "body_tag"]
    t.LowerThreshold = 1.0
    t.UpperThreshold = 1.0
    t.ThresholdMethod = "Between"
    return t


def current_through(normal, keep_normal):
    """Integrate J.n over the half-plane cut of the ring."""
    sl = Slice(Input=ring("J_field.vtk"))
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0.0, 0.0, 0.0]
    sl.SliceType.Normal = list(normal)
    UpdatePipeline(proxy=sl)
    cl = Clip(Input=sl)
    cl.ClipType = "Plane"
    cl.ClipType.Origin = [0.0, 0.0, 0.0]
    cl.ClipType.Normal = list(keep_normal)
    cl.Invert = 0
    UpdatePipeline(proxy=cl)
    iv = IntegrateVariables(Input=cl)
    iv.DivideCellDataByVolume = 0
    UpdatePipeline(proxy=iv)
    cd = sm.Fetch(iv).GetCellData()
    jr = cd.GetArray("J_cell_real").GetTuple(0)
    ji = cd.GetArray("J_cell_imag").GetTuple(0)
    k = int(np.argmax(np.abs(normal)))
    return complex(jr[k], ji[k])


print("port and circulation")
I_cut = current_through((0, 1, 0), (1, 0, 0))      # the cut, at theta = 0
I_far = -current_through((1, 0, 0), (0, 1, 0))     # theta = pi/2, nothing set
check("|I| through the cut", abs(I_cut),
      exp["drive_current"][0], exp["drive_current"][1])
check("|I| at theta=pi/2", abs(I_far),
      exp["drive_current"][0], exp["circulation_tol"][0])

# --- the terminal voltage: the Phi jump across the cut ----------------------
# One side of the cut is the port's 0 V reference; the other floats to V. Phi is
# single valued on each side, so max-min over the ring IS the terminal voltage.
pw = ring("potential.vtk")
UpdatePipeline(proxy=pw)
d = sm.Fetch(pw)
pr = vtk_to_numpy(d.GetPointData().GetArray("phi_real"))
pi_ = vtk_to_numpy(d.GetPointData().GetArray("phi_imag"))
# Average over each plateau rather than picking single extreme nodes, so one
# stray node cannot set the answer, and take BOTH components from the same set.
hi = pr > pr.max() - 1e-12 * max(abs(pr.max()), 1.0) - 1e-15
lo = pr < pr.min() + 1e-12 * max(abs(pr.max()), 1.0) + 1e-15
V = complex(pr[hi].mean() - pr[lo].mean(), pi_[hi].mean() - pi_[lo].mean())

# DIVIDE BY THE IMPOSED CURRENT, NOT THE MEASURED ONE, and here is why.
# omega*L/R is 0.017 on this loop: the reactance is 1.7 % of the resistance. So
# Im(Z) = Im(V/I) is a small difference sitting on a large product, and the
# measured current's own 0.02 % imaginary part leaks into it -- it moves L by
# 1.3 %, fifty times its own size. The port IMPOSES exactly 1 A; that the
# measured current comes back as 0.99995 is the check above, not the
# definition. Using the imposed value makes L agree with the completely
# independent energy route, integral(|B|^2/mu)/|I|^2, to four digits.
I_drive = complex(exp["drive_current"][0], 0.0)
Z = V / I_drive

print("\nimpedance")
print("  V = %.6e %+.6e j V   (the Phi jump across the cut)" % (V.real, V.imag))
print("  Z = %.6e %+.6e j ohm" % (Z.real, Z.imag))

# Exact dc resistance of the revolved polygon, straight off the mesh:
# integral(dA/rho) = (1/2pi) integral(dV/rho^2), which assumes no shape.
szr = CellSize(Input=ring("J_field.vtk"))
szr.ComputeVolume = 1
ccr = CellCenters(Input=szr)
ccr.VertexCells = 1
UpdatePipeline(proxy=ccr)
dr = sm.Fetch(ccr)
qr = vtk_to_numpy(dr.GetPoints().GetData())
vol = np.abs(vtk_to_numpy(dr.GetPointData().GetArray("Volume")))
rho = np.hypot(qr[:, 0], qr[:, 1])
R_exact = 2 * math.pi / (sigma * (vol / rho ** 2).sum() / (2 * math.pi))
check("R vs exact for this shape", Z.real / R_exact, 1.0, exp["R_rel_tol"][0])
check("L (nH)", Z.imag / omega * 1e9, exp["L_nH"][0], exp["L_nH"][1])

print()
if failures:
    print("FAILED %d of %d checks: %s" % (len(failures), checks, ", ".join(failures)))
    sys.exit(1)
print("all %d physics checks passed" % checks)
