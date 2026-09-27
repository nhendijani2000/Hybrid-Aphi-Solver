from paraview.simple import *
from paraview import servermanager as sm
import math, os

HERE = r"C:\Research\APhi_Solver_Project_LowFrequency_EDA\APhi_Solver"
src = LegacyVTKReader(registrationName='p', FileNames=[os.path.join(HERE, 'potential.vtk')])
wire = Threshold(registrationName='w', Input=src)
wire.Scalars = ['CELLS', 'body_tag']
wire.LowerThreshold = 1.0; wire.UpperThreshold = 1.0; wire.ThresholdMethod = 'Between'

iv = IntegrateVariables(registrationName='iv', Input=wire)
iv.DivideCellDataByVolume = 0
d = sm.Fetch(iv)
cd = d.GetCellData()
vol = cd.GetArray('Volume').GetValue(0)
jr = cd.GetArray('J_real').GetTuple(0)
ji = cd.GetArray('J_imag').GetTuple(0)

L = 1.0e-3     # wire length, m
V = 1.0        # prescribed terminal voltage

# I = (1/L) * integral of J_z over the wire volume.
# Current OUT of the high-potential terminal: Phi is 1 V at the top, so the
# conduction current runs in -z, and the terminal current is its negative.
I_re, I_im = -jr[2] / L, -ji[2] / L
I_mag = math.hypot(I_re, I_im)

# R from the real power, which is the definition that survives a reactive part:
# the series impedance seen by the terminals is Z = V / I.
den = I_re * I_re + I_im * I_im
Z_re, Z_im = V * I_re / den, -V * I_im / den

# Expected: the mesh wire is a 24-gon inscribed in r = 0.2 mm, not a circle.
N, a, sigma = 24, 2.0e-4, 5.8e7
area = 0.5 * N * a * a * math.sin(2.0 * math.pi / N)
R_exact = L / (sigma * area)

print("  wire volume          %.6e m3   (exact %.6e)" % (vol, area * L))
print("  integral J_z dV      %.6e A.m" % jr[2])
print()
print("  terminal current I   %.4f %+.4f j A   |I| = %.4f A" % (I_re, I_im, I_mag))
print("  impedance Z = V/I    %.6e %+.6e j ohm" % (Z_re, Z_im))
print()
print("  R  measured          %.6e ohm" % Z_re)
print("  R  exact  L/(sigma*A)%.6e ohm" % R_exact)
print("  relative error       %.3e   -> %.6f %%" % (abs(Z_re - R_exact) / R_exact,
                                                    100 * abs(Z_re - R_exact) / R_exact))
print()
print("  L  from Im(Z)/omega  %.6e H   (omega = %.6f)" % (Z_im / (2 * math.pi * 50),
                                                          2 * math.pi * 50))

L_ext_target = (4e-7 * math.pi * L / (2 * math.pi)) * math.log(1.0787 * 2.0e-3 / (2 * a))
L_int_target = 4e-7 * math.pi * L / (8 * math.pi)
L_target = L_ext_target + L_int_target
L_meas = Z_im / (2 * math.pi * 50)
print("  L  target ext+int    %.6e H  (%.4f nH = %.4f ext + %.4f int)"
      % (L_target, L_target * 1e9, L_ext_target * 1e9, L_int_target * 1e9))
print("  L  measured          %.6e H  (%.4f nH)" % (L_meas, L_meas * 1e9))
print("  difference           %.3f %%" % (100 * abs(L_meas - L_target) / L_target))
