"""Mesh, solve and measure case 06 at three refinements. Writes convergence.json.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_convergence_data.py

Takes a few minutes: three meshes and three solves. The result feeds
make_report_figs.py's fig/convergence.png, and the report's "which mesh" section
quotes it.

WHY THREE AND NOT TWO. Two points fit any power law you care to assume -- they
cannot test one. An earlier version of this case claimed first-order convergence
in h from two meshes; the third point showed the apparent order collapsing, which
turned out to be a constant offset in the COMPARISON rather than anything in the
solve (see the note at the top of tools/skin_profile.py). Three points do not
establish an order either, but they can expose a claim that is not there.
"""
import json
import math
import os
import re
import subprocess
import sys

import numpy as np
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

GMSH = os.path.expandvars(
    r"%LOCALAPPDATA%\Microsoft\WinGet\Packages"
    r"\gmsh.gmsh_Microsoft.Winget.Source_8wekyb3d8bbwe"
    r"\gmsh-4.13.1-Windows64\gmsh.exe")
LEVELS = ("0.60", "0.40", "0.25")

A, D, LZ, SIG, F = 2.0e-3, 10.0e-3, 20.0e-3, 5.8e7, 2500.0
MU0 = 4e-7 * math.pi
DELTA = math.sqrt(2.0 / (2 * math.pi * F * MU0 * SIG))
K = (1 - 1j) / DELTA


def j0(z):
    t, s = 1.0 + 0j, 1.0 + 0j
    for m in range(1, 80):
        t *= -(z * z) / (4.0 * m * m)
        s += t
        if abs(t) < 1e-18 * abs(s):
            break
    return s


def build(lc):
    tag = lc.replace(".", "")
    geo, msh = "conv_%s.geo" % tag, "conv_%s.msh" % tag
    aphi, out = "conv_%s.aphi" % tag, "conv_out_%s" % tag
    src = open("two_wires.geo", encoding="utf-8").read()
    open(geo, "w", encoding="utf-8").write(
        re.sub(r"^lc_wire = [0-9.]+;.*$", "lc_wire = %s;" % lc, src, flags=re.M))
    if not os.path.isfile(GMSH):
        raise SystemExit("gmsh not found at %s" % GMSH)
    subprocess.run([GMSH, geo, "-3", "-o", msh], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    src = open("two_wires.aphi", encoding="utf-8").read()
    src = src.replace("file        = two_wires.msh", "file        = " + msh)
    src = src.replace("directory   = output", "directory   = " + out)
    open(aphi, "w", encoding="utf-8").write(src)
    subprocess.run(["cmd", "/c", os.path.join("..", "run_case.bat"), aphi],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return out


def measure(out):
    s = LegacyVTKReader(FileNames=[os.path.join(out, "J_field.vtk")])
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
    r = np.hypot(q[:, 0] + 0.5 * D, q[:, 1])
    mid = (q[:, 2] > 0.35 * LZ) & (q[:, 2] < 0.65 * LZ)
    w1 = mid & (bt == 1)
    a = math.sqrt((vol[w1].sum() / (0.30 * LZ)) / math.pi)
    surf = w1 & (r > 0.92 * a)
    core = w1 & (r < 0.2 * a)
    Ja = J[surf].mean()
    # the reference is a MEAN over cells, so Bessel is meaned over the same radii
    bref = np.mean([j0(K * x) for x in r[surf]])
    ref = np.mean([j0(K * x) for x in r[core]]) / bref
    Jc = J[core].mean()
    return (float(np.median(h[w1])) * 1e3,
            100.0 * (abs(Jc) / abs(Ja) / abs(ref) - 1.0),
            math.degrees(np.angle(Jc / Ja)) - math.degrees(np.angle(ref)),
            int(w1.sum()))


rows = []
for lc in LEVELS:
    print("  lc_wire = %s ..." % lc)
    out = build(lc)
    rows.append(measure(out))
    print("     h = %.3f mm (%.1f per delta), %d tets, |J| %+.3f %%, lag %+.3f deg"
          % (rows[-1][0], DELTA * 1e3 / rows[-1][0], rows[-1][3], rows[-1][1], rows[-1][2]))

json.dump({"lc_mm": [float(x) for x in LEVELS],
           "h_mm": [r[0] for r in rows],
           "mag_err_pct": [r[1] for r in rows],
           "lag_err_deg": [r[2] for r in rows],
           "tets": [r[3] for r in rows],
           "delta_mm": DELTA * 1e3},
          open("convergence.json", "w"), indent=2)
print("  wrote convergence.json")
print("  (the conv_* meshes and outputs are left in place; delete them when done)")
