"""The radial |J| profile of a round conductor against Kelvin, for any case.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" skin_profile.py <output-dir>

Everything it needs is in the output: the frequency from the `.out` header, the
conductors from the per-cell `sigma` array, and each conductor's axis and radius
from the geometry of its own cells. Nothing is passed in and nothing is assumed
about which case produced it, so the same table can be read across cases.

Assumes each conductor is a straight round rod. It reports the radius it
inferred and how round the cross-section actually is, so a conductor that is
neither will say so rather than quietly produce a plausible table.

BOTH SIDES MUST BE ANCHORED AT THE SAME RADIUS, which is the trap this file
exists to stop repeating.

Kelvin gives Jz(r)/Jz(a) = J0(kr)/J0(ka), normalised at the surface r = a. A
measurement cannot be: a cell band "at the surface" is r > 0.92a, whose MEAN
radius is about 0.957a. Dividing the measurement by J(0.957a) while dividing the
theory by J(a) compares two profiles pinned at different places, and inflates
every measured ratio by roughly 1/0.966 = +3.5 %.

That offset contains no h, so refining the mesh cannot remove it -- and it once
produced a thoroughly convincing false result here: errors of 7.2 % that
"improved" to 3.9 % under refinement and then stalled, with an apparent
convergence order collapsing to 0.37, and a tidy physical story to explain the
floor. All of it was the constant offset sitting on a real error ten times
smaller. Two controls placed it: a single isolated wire showed the identical
error, and narrowing the axial sample band changed nothing.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import math
import os
import re
import sys

import numpy as np

MU0 = 4e-7 * math.pi
OUT = sys.argv[1] if len(sys.argv) > 1 else "output"


def bessel_j0(z):
    """J0 of a complex argument, by its defining series."""
    term, total = 1.0 + 0j, 1.0 + 0j
    for m in range(1, 80):
        term *= -(z * z) / (4.0 * m * m)
        total += term
        if abs(term) < 1e-18 * abs(total):
            break
    return total


def frequency_from_header(out_dir):
    """The solver writes `frequency <f> Hz` into every .out header."""
    for name in ("J_field.out", "E_field.out", "potential.out"):
        p = os.path.join(out_dir, name)
        if not os.path.isfile(p):
            continue
        with open(p, "r") as f:
            for _ in range(40):
                line = f.readline()
                if not line:
                    break
                m = re.search(r"frequency\s+([0-9.eE+-]+)\s*Hz", line)
                if m:
                    return float(m.group(1))
    raise SystemExit("could not find the frequency in %s/*.out" % out_dir)


f_hz = frequency_from_header(OUT)
omega = 2.0 * math.pi * f_hz

s = LegacyVTKReader(FileNames=[os.path.join(OUT, "J_field.vtk")])
sz = CellSize(Input=s)
sz.ComputeVolume = 1
cc = CellCenters(Input=sz)
cc.VertexCells = 1
UpdatePipeline(proxy=cc)
d = sm.Fetch(cc)

q = vtk_to_numpy(d.GetPoints().GetData())
vol = np.abs(vtk_to_numpy(d.GetPointData().GetArray("Volume")))
bt = vtk_to_numpy(d.GetPointData().GetArray("body_tag"))
sig = vtk_to_numpy(d.GetPointData().GetArray("sigma"))
Jv = (vtk_to_numpy(d.GetPointData().GetArray("J_cell_real"))
      + 1j * vtk_to_numpy(d.GetPointData().GetArray("J_cell_imag")))
h = (12.0 * vol / math.sqrt(2.0)) ** (1.0 / 3.0)

print("  %s   f = %g Hz" % (OUT, f_hz))

for tag in sorted(set(bt[sig > 0.0].tolist())):
    body = (bt == tag) & (sig > 0.0)
    if body.sum() < 200:
        continue
    sigma = float(np.median(sig[body]))
    delta = math.sqrt(2.0 / (omega * MU0 * sigma))
    k = (1.0 - 1.0j) / delta

    # The axis is the direction the conductor is longest in.
    span = [q[body, i].max() - q[body, i].min() for i in range(3)]
    ax = int(np.argmax(span))
    per = [i for i in range(3) if i != ax]
    c0, c1 = q[body, per[0]].mean(), q[body, per[1]].mean()
    rr = np.hypot(q[:, per[0]] - c0, q[:, per[1]] - c1)
    L = span[ax]
    zlo = q[body, ax].min() + 0.35 * L
    zhi = q[body, ax].min() + 0.65 * L
    sel = body & (q[:, ax] > zlo) & (q[:, ax] < zhi)

    # THE RADIUS COMES FROM THE AREA, not from the outermost cell. Cell CENTRES
    # never reach the surface -- the outermost sits half an element inside it --
    # so a percentile of their radii under-reads a by a couple of per cent, and
    # under-reading a inflates the reported a/delta. The cross-section area is
    # unbiased: a = sqrt(area/pi). The ratio of the two is printed below as a
    # roundness check, since it is only meaningful if the section really is
    # round.
    area = vol[sel].sum() / (zhi - zlo)
    a = math.sqrt(area / math.pi)
    a_cells = float(np.percentile(rr[body], 99.5))
    print()
    print("  body_tag %d   sigma %.3g S/m   a = %.3f mm   L = %.1f mm   axis %s"
          % (tag, sigma, a * 1e3, L * 1e3, "xyz"[ax]))
    print("    delta = %.4f mm,  a/delta = %.3f,  median element %.3f mm (%.1f per delta)"
          % (delta * 1e3, a / delta, float(np.median(h[sel])) * 1e3,
             delta / float(np.median(h[sel]))))
    print("    outermost cell centre at %.4f a, as expected just inside the surface"
          % (a_cells / a))

    J = Jv[:, ax]
    surf = sel & (rr > 0.92 * a)
    if surf.sum() < 20:
        print("    too few cells near the surface to anchor on")
        continue
    rbar = float(rr[surf].mean())
    Ja = J[surf].mean()
    print("    anchored at r = %.4f a, BOTH sides -- not at a. See the note above."
          % (rbar / a))
    print("    r/a         measured   Kelvin     err        lag      Kelvin    err")
    for lo, hi in ((0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.92)):
        m = sel & (rr >= lo * a) & (rr < hi * a)
        if m.sum() < 20:
            continue
        Jm = J[m].mean()
        ref = bessel_j0(k * rr[m].mean()) / bessel_j0(k * rbar)
        lag = math.degrees(np.angle(Jm / Ja))
        rlag = math.degrees(np.angle(ref))
        print("    %.2f-%.2f    %8.4f   %8.4f  %+6.2f %%  %+7.2f  %+7.2f  %+5.2f deg"
              % (lo, hi, abs(Jm) / abs(Ja), abs(ref),
                 100 * (abs(Jm) / abs(Ja) / abs(ref) - 1), lag, rlag, lag - rlag))
