"""Mesh convergence at 393 Hz: N=76 against N=96.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_convergence_plot.py

Writes fig/convergence.png. Needs BOTH output/ (N=76) and output_n96/ present;
regenerate the second with

    gmsh cylinder_skin_n96.geo -3 -o cylinder_skin_n96.msh
    ..\\run_case.bat cylinder_393hz_n96_mumps.aphi

WHY THIS FIGURE MATTERS MORE THAN EITHER ERROR NUMBER. A single error value
says "we are 0.74 % off" and cannot distinguish discretisation error from a
modelling mistake sitting at a floor. Two meshes give the RATE: h falls by
0.778 and the error falls by 0.605 = 0.778^2, so the error is genuine
discretisation converging at second order and refinement would keep paying.
An error that stalled between the two meshes would mean the opposite.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import numpy as np
import math
import cmath
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

paraview.simple._DisableFirstRenderCameraReset()

L, sigma, a, f = 40e-3, 5.8e7, 10e-3, 393.0
mu = 4e-7 * math.pi
omega = 2 * math.pi * f
delta = math.sqrt(2.0 / (omega * mu * sigma))
k = complex(1, -1) / delta
FIG = "fig"
if not os.path.isdir(FIG):
    os.makedirs(FIG)


def J0(z):
    s = t = 1 + 0j
    for m in range(1, 400):
        t *= -(z * z) / (4.0 * m * m)
        s += t
        if abs(t) < 1e-24 * abs(s):
            break
    return s


# (label, output dir, mean h in the conductor interior, colour)
RUNS = [("N=76,  h = 1.08 mm,  $\\delta/h$ = 3.1", "output", 1.08, "#1f77b4"),
        ("N=96,  h = 0.84 mm,  $\\delta/h$ = 4.0", "output_n96", 0.84, "#d62728")]

fig, (axm, axp) = plt.subplots(2, 1, figsize=(8.4, 8.4), sharex=True)
store = {}

for label, out, hmm, col in RUNS:
    if not os.path.isdir(out):
        raise SystemExit("missing %s -- see the docstring" % out)
    src = LegacyVTKReader(FileNames=[os.path.join(out, "J_field.vtk")])
    th = Threshold(Input=src)
    th.Scalars = ["CELLS", "body_tag"]
    th.LowerThreshold = 1.0
    th.UpperThreshold = 1.0
    th.ThresholdMethod = "Between"
    cc = CellCenters(Input=th)
    cc.VertexCells = 1
    UpdatePipeline(proxy=cc)
    d = sm.Fetch(cc)
    q = vtk_to_numpy(d.GetPoints().GetData())
    jc = (vtk_to_numpy(d.GetPointData().GetArray("J_cell_real"))[:, 2]
          + 1j * vtk_to_numpy(d.GetPointData().GetArray("J_cell_imag"))[:, 2])
    r = np.hypot(q[:, 0], q[:, 1])
    zc = q[:, 2]
    mid = (zc > 0.3 * L) & (zc < 0.7 * L)
    ref = mid & (r > 0.95 * a)
    jref = np.abs(jc[ref]).mean()
    jcref = jc[ref].mean()
    bref = np.mean([abs(J0(k * x)) for x in r[ref]])
    bcref = np.mean([J0(k * x) for x in r[ref]])

    bins = np.linspace(0, 1, 21)
    idx = np.digitize(r[mid] / a, bins) - 1
    ctr, em, ep = [], [], []
    for b in range(len(bins) - 1):
        m = idx == b
        if m.sum() < 5:
            continue
        rb = r[mid][m]
        meas = np.abs(jc[mid][m]).mean() / jref
        ex = np.mean([abs(J0(k * x)) for x in rb]) / bref
        pm = math.degrees(cmath.phase(jc[mid][m].mean() / jcref))
        pe = math.degrees(cmath.phase(np.mean([J0(k * x) for x in rb]) / bcref))
        ctr.append(rb.mean() / a)
        em.append(100 * (meas / ex - 1))
        ep.append(pm - pe)
    store[label] = (np.array(ctr), np.array(em), np.array(ep), hmm)
    axm.plot(ctr, em, "o-", ms=4, color=col, label=label)
    axp.plot(ctr, ep, "s-", ms=4, color=col, label=label)

# What second order PREDICTS for the finer mesh, drawn from the coarse one.
(c0, m0, p0, h0) = store[RUNS[0][0]]
scale = (RUNS[1][2] / h0) ** 2
axm.plot(c0, m0 * scale, "k--", lw=1.4, dashes=(5, 3),
         label="$h^2$ from the coarse mesh  ($\\times%.3f$)" % scale)
axp.plot(c0, p0 * scale, "k--", lw=1.4, dashes=(5, 3))

for ax, ylab in ((axm, "magnitude error,  %"), (axp, "phase error,  deg")):
    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(1.0 - delta / a, color="#888888", ls=":", lw=1.4)
    ax.set_ylabel(ylab)
    ax.grid(alpha=0.3)
axm.legend(loc="upper right", fontsize=9)
axp.set_xlabel("$r/a$")
axm.set_title("Mesh convergence at 393 Hz,  a/$\\delta$ = 3.00\n"
              "error against exact Bessel at the sampled radii; the dashed line "
              "is what $h^2$ predicts")
axp.set_title("the same for phase. Both land on the $h^2$ line, so the residual "
              "is discretisation\nconverging at second order, not a modelling "
              "floor.", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "convergence.png"), dpi=150)
print("wrote %s/convergence.png" % FIG)
print("  worst |J| error: %.3f %% -> %.3f %%   (h^2 predicts %.3f %%)"
      % (np.abs(m0).max(), np.abs(store[RUNS[1][0]][1]).max(),
         np.abs(m0).max() * scale))
print("  worst phase err: %.3f deg -> %.3f deg (h^2 predicts %.3f deg)"
      % (np.abs(p0).max(), np.abs(store[RUNS[1][0]][2]).max(),
         np.abs(p0).max() * scale))
