"""Terminal current, R and L, by integrating J = sigma E over the conductor.

    pvbatch tools/pv_extract_rl.py <output_dir> [a_mm] [L_mm] [N] [V] [f_Hz]

Defaults are regression_tests/01_OneCylinder: a = 10 mm, L = 40 mm, a 36-gon,
driven by 1 V at 50 Hz.

Reads J_field.vtk, which carries the per-cell J. Per-cell is the only place J
is defined: sigma is single-valued inside a tet and is not at a node on a
conductor/insulator boundary.

The current is taken as -(1/L) * integral(J_z dV). The sign matters: Phi is
higher at the top, so the conduction current runs in -z, and the current out of
the high-potential terminal is its negative. Getting that backwards reports a
negative resistance, which is how it was first caught.

R is compared against the POLYGON's cross-section, not pi*a^2 -- the meshed
conductor is an N-gon and is 0.51 % smaller at N = 36.
"""

from paraview.simple import *
from paraview import servermanager as sm
import math
import os
import sys

HERE = sys.argv[1] if len(sys.argv) > 1 else "."
A_MM = float(sys.argv[2]) if len(sys.argv) > 2 else 10.0
L_MM = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
NGON = int(sys.argv[4]) if len(sys.argv) > 4 else 36
VOLT = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
FREQ = float(sys.argv[6]) if len(sys.argv) > 6 else 50.0

vtk = os.path.join(HERE, "J_field.vtk")
if not os.path.exists(vtk):
    raise SystemExit("no such file: %s\nrun the case first" % vtk)

src = LegacyVTKReader(registrationName="J", FileNames=[vtk])
wire = Threshold(registrationName="wire", Input=src)
wire.Scalars = ["CELLS", "body_tag"]
wire.LowerThreshold = 1.0
wire.UpperThreshold = 1.0
wire.ThresholdMethod = "Between"

iv = IntegrateVariables(registrationName="iv", Input=wire)
iv.DivideCellDataByVolume = 0
data = sm.Fetch(iv)
cd = data.GetCellData()
vol = cd.GetArray("Volume").GetValue(0)
jr = cd.GetArray("J_cell_real").GetTuple(0)
ji = cd.GetArray("J_cell_imag").GetTuple(0)

a = A_MM * 1e-3
L = L_MM * 1e-3
sigma = 5.8e7

# Current OUT of the high-potential terminal: see the note above on the sign.
I_re, I_im = -jr[2] / L, -ji[2] / L
den = I_re * I_re + I_im * I_im
Z_re, Z_im = VOLT * I_re / den, -VOLT * I_im / den

area = 0.5 * NGON * a * a * math.sin(2.0 * math.pi / NGON)
R_exact = L / (sigma * area)

print("  conductor volume     %.6e m3   (exact %.6e)" % (vol, area * L))
print("  integral J_z dV      %.6e A.m" % jr[2])
print()
print("  terminal current I   %.4f %+.4f j A   |I| = %.4f A"
      % (I_re, I_im, math.hypot(I_re, I_im)))
print("  impedance Z = V/I    %.6e %+.6e j ohm" % (Z_re, Z_im))
print()
print("  R  measured          %.6e ohm" % Z_re)
print("  R  exact L/(sigma A) %.6e ohm" % R_exact)
print("  relative error       %.3e" % (abs(Z_re - R_exact) / R_exact))
print()
print("  L  from Im(Z)/omega  %.6e H  (%.4f nH)"
      % (Z_im / (2 * math.pi * FREQ), Z_im / (2 * math.pi * FREQ) * 1e9))
