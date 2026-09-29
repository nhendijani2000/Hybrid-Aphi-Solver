"""Assert the PHYSICS of a solved regression case, not just that it ran.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" verify.py <case-dir> [output-dir]

Reads <case-dir>/expected.txt, measures the same quantities from the case's
output, and exits non-zero on the first that drifts. Prints every check either
way, with the margin, so a check that is close to failing is visible before it
fails.

WHY THIS EXISTS. R, L, Phi and J(0)/J(a) were re-checked BY HAND, through
ParaView, a dozen times in one day while the solver was changed underneath them
-- supervariables, the analysis reuse, two MUMPS backends. Nothing in
run-tests.bat would have caught it if one of those had moved the physics.

BOTH V AND I ARE MEASURED, never assumed, whichever way the case is driven:

  * I from the volume integral of J over the conductor, I = -(1/L) int J.z dV.
    A volume form rather than a face integral because J = sigma E is exact per
    tetrahedron and has no material interface to straddle.
  * V as the COMPLEX potential difference between the two end caps. Taking only
    the real part silently throws the reactance away -- on a current-driven
    case that reports L = 0 while everything else looks right.

Measuring both means a voltage-driven case also checks that the prescribed V
came back as prescribed, and a current-driven case checks that the injected I
is the I that actually flows. That second one is port current
self-consistency, SOLVER_PLAN.md Sec. 8 item 7.

WHAT IT DOES NOT DO. It does not run the solve -- use check.bat for that. The
separation is deliberate: the solve is minutes, this is seconds, and
re-verifying existing output without re-solving is what makes it usable while
iterating on post-processing.
"""

from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import numpy as np
import math
import os
import sys

paraview.simple._DisableFirstRenderCameraReset()

CASE = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(CASE, "output")

failures = []
checks = 0


def check(name, value, target, tol):
    global checks
    checks += 1
    ok = abs(value - target) <= tol
    margin = abs(value - target) / tol if tol > 0 else 0.0
    print("  %-30s %-15.7g vs %-15.7g tol %-9.3g  %s (%.0f%% of tol)"
          % (name, value, target, tol, "ok  " if ok else "FAIL", 100 * margin))
    if not ok:
        failures.append(name)


def check_max(name, value, limit):
    global checks
    checks += 1
    ok = value <= limit
    print("  %-30s %-15.7g <= %-15.7g           %s (%.0f%% of limit)"
          % (name, value, limit, "ok  " if ok else "FAIL", 100 * value / limit))
    if not ok:
        failures.append(name)


# --- expected.txt ----------------------------------------------------------
exp = {}
exp_path = os.path.join(CASE, "expected.txt")
if not os.path.exists(exp_path):
    raise SystemExit("no expected.txt in %s" % CASE)
for line in open(exp_path):
    line = line.split("#")[0].strip()
    if not line:
        continue
    parts = line.split()
    exp[parts[0]] = [float(v) for v in parts[1:]]

a = exp["a_mm"][0] * 1e-3
L = exp["length_mm"][0] * 1e-3
ngon = int(exp["ngon"][0])
sigma = exp["sigma"][0]
freq = exp["frequency"][0]
omega = 2.0 * math.pi * freq

drive_v = exp.get("drive_voltage")
drive_i = exp.get("drive_current")
if (drive_v is None) == (drive_i is None):
    raise SystemExit("expected.txt needs exactly one of drive_voltage / drive_current")

print("verifying %s" % os.path.abspath(CASE))
print("  output   %s" % os.path.abspath(OUT))
print("  driven by %s" % ("%g V" % drive_v[0] if drive_v else "%g A" % drive_i[0]))
print()

for f in ("J_field.vtk", "potential.vtk"):
    if not os.path.exists(os.path.join(OUT, f)):
        raise SystemExit("missing %s -- solve the case first" % os.path.join(OUT, f))


def wire_of(path):
    src = LegacyVTKReader(FileNames=[os.path.join(OUT, path)])
    t = Threshold(Input=src)
    t.Scalars = ["CELLS", "body_tag"]
    t.LowerThreshold = 1.0
    t.UpperThreshold = 1.0
    t.ThresholdMethod = "Between"
    return t


# --- I, from the volume integral of J --------------------------------------
wire = wire_of("J_field.vtk")
iv = IntegrateVariables(Input=wire)
iv.DivideCellDataByVolume = 0
UpdatePipeline(proxy=iv)
cd = sm.Fetch(iv).GetCellData()
jr = cd.GetArray("J_cell_real").GetTuple(0)
ji = cd.GetArray("J_cell_imag").GetTuple(0)
I = complex(-jr[2] / L, -ji[2] / L)

# --- V, complex, cap to cap ------------------------------------------------
pw = wire_of("potential.vtk")
UpdatePipeline(proxy=pw)
d = sm.Fetch(pw)
pts = vtk_to_numpy(d.GetPoints().GetData())
pr = vtk_to_numpy(d.GetPointData().GetArray("phi_real"))
pi = vtk_to_numpy(d.GetPointData().GetArray("phi_imag"))
z = pts[:, 2]
top = z > L - 1e-9
bot = z < 1e-9
if top.sum() < 10 or bot.sum() < 10:
    raise SystemExit("end caps not found: %d top, %d bottom nodes" % (top.sum(), bot.sum()))
V = complex(pr[top].mean() - pr[bot].mean(), pi[top].mean() - pi[bot].mean())
Z = V / I

print("measured")
print("  I = %.6f %+.6f j A        |I| = %.6f A" % (I.real, I.imag, abs(I)))
print("  V = %.6e %+.6e j V   |V| = %.6e V" % (V.real, V.imag, abs(V)))
print("  Z = %.6e %+.6e j ohm" % (Z.real, Z.imag))
print()

# --- the drive came back as driven -----------------------------------------
print("drive")
if drive_i is not None:
    # Port current self-consistency: inject 1 A at one cap, and exactly that
    # must flow. Charge does not accumulate in a conductor, so anything else
    # means current is being created or destroyed in the formulation.
    check("|I| collected vs injected", abs(I), drive_i[0], drive_i[1])
else:
    check("|V| measured vs prescribed", abs(V), drive_v[0], drive_v[1])

# --- R and L ---------------------------------------------------------------
print("\nR and L")
area = 0.5 * ngon * a * a * math.sin(2.0 * math.pi / ngon)
R_exact = L / (sigma * area)
check_max("R rel error vs exact", abs(Z.real - R_exact) / R_exact, exp["R_rel_error_max"][0])
check("L (nH)", Z.imag / omega * 1e9, exp["L_nH"][0], exp["L_nH"][1])

# --- Phi against the exact z/L --------------------------------------------
# Normalised by the real terminal voltage, so the same limit applies whether
# the case is driven with 1 V or with a current that produces 98 microvolts.
print("\nPhi")
dev = np.abs(pr / V.real - z / L)
check_max("Phi/V vs z/L, worst", float(dev.max()), exp["phi_vs_zL_worst_max"][0])

# A second, STRICTER form: pointwise relative error against the exact linear
# profile, rather than error normalised by the terminal voltage. The two differ
# because the normalised form divides a worst-case deviation by the LARGEST
# potential in the problem, so a deviation sitting near the low-potential end
# is judged against a value it never approaches. The pointwise form asks what
# it should: how wrong is Phi HERE, relative to what it should be HERE.
#
# Guarded to z > 0.1 L. The exact profile passes through zero at the reference
# cap, and a relative error against zero is not a number.
mrel = z > 0.1 * L
exact = V.real * z / L
point_rel = np.abs(pr[mrel] / exact[mrel] - 1.0)
check_max("Phi pointwise rel, z>0.1L", float(point_rel.max()),
          exp["phi_pointwise_rel_max"][0])
# --- J(0)/J(a) against Bessel ---------------------------------------------
print("\nJ profile")
cc = CellCenters(Input=wire)
cc.VertexCells = 1
calc = Calculator(Input=cc)
calc.AttributeType = "Point Data"
calc.ResultArrayName = "Jz"
calc.Function = "sqrt(J_cell_real_Z^2 + J_cell_imag_Z^2)"
UpdatePipeline(proxy=calc)
dj = sm.Fetch(calc)
q = vtk_to_numpy(dj.GetPoints().GetData())
jz = vtk_to_numpy(dj.GetPointData().GetArray("Jz"))
r = np.hypot(q[:, 0], q[:, 1])
zc = q[:, 2]
mid = (zc > 0.3 * L) & (zc < 0.7 * L)
core = mid & (r < 0.25 * a)
surf = mid & (r > 0.90 * a)
if core.sum() < 5 or surf.sum() < 5:
    raise SystemExit("too few J samples: core %d, surface %d" % (core.sum(), surf.sum()))
check("J(0)/J(a)", float(jz[core].mean() / jz[surf].mean()),
      exp["J_ratio_exact"][0], exp["J_ratio_tol"][0])

# --- verdict ---------------------------------------------------------------
print()
if failures:
    print("FAILED %d of %d checks: %s" % (len(failures), checks, ", ".join(failures)))
    sys.exit(1)
print("all %d physics checks passed" % checks)
