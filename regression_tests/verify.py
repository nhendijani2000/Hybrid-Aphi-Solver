"""Assert the PHYSICS of a solved regression case, not just that it ran.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" verify.py <case-dir> [output-dir]

Reads <case-dir>/expected.txt, measures the same quantities from the case's
output, and exits non-zero on the first that drifts. Prints every check either
way, so a pass is as informative as a failure.

WHY THIS EXISTS. R, L, Phi and J(0)/J(a) were re-checked BY HAND, through
ParaView, a dozen times in one day while the solver was being changed
underneath them -- supervariables, the analysis reuse, two MUMPS backends.
Nothing in run-tests.bat would have caught it if one of those had moved the
physics. Every check here is one that was previously done by eye.

WHAT IT DOES NOT DO. It does not run the solve. Run the case first:

    cd 02_Ansys_Cylinder_50Hz && ..\\run_case.bat cylinder_50hz.aphi

That separation is deliberate: the solve is minutes, the verification is
seconds, and being able to re-verify an existing output without re-solving is
what makes it usable while iterating on the post-processing.
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


def check(name, value, target, tol, unit=""):
    """Absolute-difference check. Prints the margin, not just pass/fail --
    a check passing by 1 % of its tolerance is worth seeing before it fails."""
    global checks
    checks += 1
    ok = abs(value - target) <= tol
    margin = abs(value - target) / tol if tol > 0 else 0.0
    print("  %-28s %-14.6g vs %-14.6g tol %-9.3g  %s  (%.0f%% of tol)"
          % (name, value, target, tol, "ok  " if ok else "FAIL", 100 * margin))
    if not ok:
        failures.append(name)


def check_max(name, value, limit):
    global checks
    checks += 1
    ok = value <= limit
    print("  %-28s %-14.6g <= %-14.6g            %s  (%.0f%% of limit)"
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
volt = exp["voltage"][0]
freq = exp["frequency"][0]
omega = 2.0 * math.pi * freq

print("verifying %s" % os.path.abspath(CASE))
print("  output   %s" % os.path.abspath(OUT))
print()

for f in ("J_field.vtk", "potential.vtk"):
    if not os.path.exists(os.path.join(OUT, f)):
        raise SystemExit("missing %s -- solve the case first" % os.path.join(OUT, f))

# --- R and L from the volume integral of J ---------------------------------
# I = -(1/L) * integral J_z dV over the conductor. A volume form, not a face
# integral: J = sigma E is exact per tet and has no material interface to
# straddle, where a face integral would need J on the boundary itself.
src = LegacyVTKReader(registrationName="J", FileNames=[os.path.join(OUT, "J_field.vtk")])
wire = Threshold(Input=src)
wire.Scalars = ["CELLS", "body_tag"]
wire.LowerThreshold = 1.0
wire.UpperThreshold = 1.0
wire.ThresholdMethod = "Between"
iv = IntegrateVariables(Input=wire)
iv.DivideCellDataByVolume = 0
UpdatePipeline(proxy=iv)
cd = sm.Fetch(iv).GetCellData()
jr = cd.GetArray("J_cell_real").GetTuple(0)
ji = cd.GetArray("J_cell_imag").GetTuple(0)

I_re, I_im = -jr[2] / L, -ji[2] / L
den = I_re * I_re + I_im * I_im
Z_re, Z_im = volt * I_re / den, -volt * I_im / den

area = 0.5 * ngon * a * a * math.sin(2.0 * math.pi / ngon)
R_exact = L / (sigma * area)

print("R and L")
check_max("R rel error vs exact", abs(Z_re - R_exact) / R_exact, exp["R_rel_error_max"][0])
check("L (nH)", Z_im / omega * 1e9, exp["L_nH"][0], exp["L_nH"][1])

# --- Phi against the exact z/L --------------------------------------------
print("\nPhi")
psrc = LegacyVTKReader(registrationName="P", FileNames=[os.path.join(OUT, "potential.vtk")])
pw = Threshold(Input=psrc)
pw.Scalars = ["CELLS", "body_tag"]
pw.LowerThreshold = 1.0
pw.UpperThreshold = 1.0
pw.ThresholdMethod = "Between"
UpdatePipeline(proxy=pw)
d = sm.Fetch(pw)
p = vtk_to_numpy(d.GetPoints().GetData())
phi = vtk_to_numpy(d.GetPointData().GetArray("phi_real"))
dev = np.abs(phi - p[:, 2] / L)
check_max("Phi vs z/L, worst", float(dev.max()), exp["phi_vs_zL_worst_max"][0])

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
z = q[:, 2]
mid = (z > 0.3 * L) & (z < 0.7 * L)
core = mid & (r < 0.25 * a)
surf = mid & (r > 0.90 * a)
if core.sum() < 5 or surf.sum() < 5:
    raise SystemExit("too few samples for the J profile: core %d, surface %d"
                     % (core.sum(), surf.sum()))
ratio = float(jz[core].mean() / jz[surf].mean())
check("J(0)/J(a)", ratio, exp["J_ratio_exact"][0], exp["J_ratio_tol"][0])

# --- verdict ---------------------------------------------------------------
print()
if failures:
    print("FAILED %d of %d checks: %s" % (len(failures), checks, ", ".join(failures)))
    sys.exit(1)
print("all %d physics checks passed" % checks)
