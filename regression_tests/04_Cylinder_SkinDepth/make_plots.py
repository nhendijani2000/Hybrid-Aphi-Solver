"""Figures for 04_Cylinder_SkinDepth.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_plots.py [output-dir]

Writes into fig/. Two kinds:

  * ParaView renderings of |J| over the conductor cross-section -- the picture
    of the current crowding into the surface layer, which is the whole point of
    running this case at 393 Hz instead of 50 Hz.
  * A matplotlib radial profile of |J(r)| against the EXACT Bessel curve, which
    is the quantitative version of the same thing. The scatter is plotted per
    cell rather than binned, so element noise is visible instead of averaged
    away -- binning into radial bands is what once hid a real mesh defect
    behind a smooth-looking curve.

Read the colour-range warning in the report before drawing conclusions from the
renderings: at this frequency |J| spans nearly 4x, so the auto range is
meaningful here, unlike in the 50 Hz cases where it magnified noise.
"""

from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import numpy as np
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

paraview.simple._DisableFirstRenderCameraReset()

OUT = sys.argv[1] if len(sys.argv) > 1 else "output_mumps"
FIG = "fig"
if not os.path.isdir(FIG):
    os.makedirs(FIG)

a, L, sigma, f, ngon = 10e-3, 40e-3, 5.8e7, 393.0, 76
mu = 4e-7 * math.pi
omega = 2 * math.pi * f
delta = math.sqrt(2.0 / (omega * mu * sigma))
k = complex(1, -1) / delta


def J0(z):
    s = t = 1 + 0j
    for m in range(1, 300):
        t *= -(z * z) / (4.0 * m * m)
        s += t
        if abs(t) < 1e-22 * abs(s):
            break
    return s


def wire(fn):
    src = LegacyVTKReader(FileNames=[os.path.join(OUT, fn)])
    t = Threshold(Input=src)
    t.Scalars = ["CELLS", "body_tag"]
    t.LowerThreshold = 1.0
    t.UpperThreshold = 1.0
    t.ThresholdMethod = "Between"
    return t


# --- the renderings --------------------------------------------------------
w = wire("J_field.vtk")
calc = Calculator(Input=w)
calc.AttributeType = "Cell Data"
calc.ResultArrayName = "Jmag"
calc.Function = "sqrt(J_cell_real_Z^2 + J_cell_imag_Z^2)"

sl = Slice(Input=calc)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, 0.0, L / 2]
sl.SliceType.Normal = [0.0, 0.0, 1.0]
UpdatePipeline(proxy=sl)

view = CreateView("RenderView")
view.ViewSize = [1100, 1000]
view.OrientationAxesVisibility = 0
view.CameraPosition = [0, 0, L / 2 + 0.05]
view.CameraFocalPoint = [0, 0, L / 2]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelProjection = 1
view.CameraParallelScale = 1.08 * a
view.Background = [1, 1, 1]
view.UseColorPaletteForBackground = 0

d = Show(sl, view)
ColorBy(d, ("CELLS", "Jmag"))
d.RescaleTransferFunctionToDataRange(True, False)
lut = GetColorTransferFunction("Jmag")
lut.ApplyPreset("Rainbow Uniform", True)
bar = GetScalarBar(lut, view)
bar.Title = "|J|  A/m2"
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]
d.SetScalarBarVisibility(view, True)

Render(view)
SaveScreenshot(os.path.join(FIG, "j_cross_section.png"), view,
               ImageResolution=[1100, 1000])
print("wrote %s/j_cross_section.png  -- |J| at mid height, auto range" % FIG)

# The same slice with the range pinned to 0..max, so the eye reads the ratio
# rather than a stretched band. At 393 Hz the two look similar, which is itself
# the point: the dynamic range is real here.
lo, hi = lut.RGBPoints[0], lut.RGBPoints[-4]
lut.RescaleTransferFunction(0.0, hi)
Render(view)
SaveScreenshot(os.path.join(FIG, "j_cross_section_zeroed.png"), view,
               ImageResolution=[1100, 1000])
print("wrote %s/j_cross_section_zeroed.png  -- same data, range pinned to 0" % FIG)

# A longitudinal slice, to show there is no z dependence to worry about.
sl2 = Slice(Input=calc)
sl2.SliceType = "Plane"
sl2.SliceType.Origin = [0.0, 0.0, 0.0]
sl2.SliceType.Normal = [0.0, 1.0, 0.0]
UpdatePipeline(proxy=sl2)
Hide(sl, view)
d2 = Show(sl2, view)
ColorBy(d2, ("CELLS", "Jmag"))
d2.RescaleTransferFunctionToDataRange(True, False)
d2.SetScalarBarVisibility(view, True)
view.CameraPosition = [0.05, 0, L / 2]
view.CameraFocalPoint = [0, 0, L / 2]
view.CameraViewUp = [0, 0, 1]
view.CameraParallelScale = 0.60 * L
Render(view)
SaveScreenshot(os.path.join(FIG, "j_longitudinal.png"), view,
               ImageResolution=[700, 1100])
print("wrote %s/j_longitudinal.png  -- no z dependence" % FIG)

# --- the radial profile, per cell, against exact Bessel --------------------
cc = CellCenters(Input=wire("J_field.vtk"))
cc.VertexCells = 1
cal = Calculator(Input=cc)
cal.AttributeType = "Point Data"
cal.ResultArrayName = "Jz"
cal.Function = "sqrt(J_cell_real_Z^2 + J_cell_imag_Z^2)"
UpdatePipeline(proxy=cal)
dj = sm.Fetch(cal)
q = vtk_to_numpy(dj.GetPoints().GetData())
jz = vtk_to_numpy(dj.GetPointData().GetArray("Jz"))
r = np.hypot(q[:, 0], q[:, 1])
zc = q[:, 2]
mid = (zc > 0.3 * L) & (zc < 0.7 * L)

ref = mid & (r > 0.95 * a)
jref = jz[ref].mean()
bref = np.array([abs(J0(k * x)) for x in r[ref]]).mean()

rr = np.linspace(0, a, 400)
exact = np.array([abs(J0(k * x)) for x in rr]) / bref

fig, (ax, axe) = plt.subplots(2, 1, figsize=(8.2, 8.0), sharex=True,
                              gridspec_kw={"height_ratios": [3, 1]})
# The exact curve is DASHED and drawn on top. A solid line here hides the
# scatter completely -- the agreement is inside the marker size -- and a figure
# where the data is invisible under the reference is not showing agreement, it
# is just showing the reference.
ax.scatter(r[mid] / a, jz[mid] / jref, s=4.0, alpha=0.45, color="#4c9be8",
           label="solver, one point per cell (%d cells)" % mid.sum(),
           rasterized=True, zorder=2)
ax.plot(rr / a, exact, "k--", lw=1.6, dashes=(6, 4), zorder=3,
        label=r"exact  $|J_0(kr)/J_0(ka)|$,  $k=(1-j)/\delta$")
ax.axvline(1.0 - delta / a, color="#d62728", ls="--", lw=1.3,
           label=r"one skin depth in,  $r = a-\delta$")
ax.set_ylabel(r"$|J(r)|$  normalised to the $r>0.95a$ band")
ax.set_title("Copper rod, a = 10 mm, 393 Hz:  a/$\\delta$ = 3.00, "
             "$\\delta$ = 3.33 mm\ncurrent crowds into the surface layer, "
             "|J| falls 3.7x to the axis")
ax.legend(loc="upper left", fontsize=9, framealpha=0.95)
ax.grid(alpha=0.3)
ax.set_ylim(0, 1.15)

# error panel, binned only here -- the scatter above is the honest view
bins = np.linspace(0, 1, 21)
idx = np.digitize(r[mid] / a, bins) - 1
me, ex, ctr = [], [], []
for b in range(len(bins) - 1):
    m = idx == b
    if m.sum() < 5:
        continue
    rb = r[mid][m]
    me.append(jz[mid][m].mean() / jref)
    ex.append(np.array([abs(J0(k * x)) for x in rb]).mean() / bref)
    ctr.append(rb.mean() / a)
err = 100.0 * (np.array(me) / np.array(ex) - 1.0)
axe.plot(ctr, err, "o-", ms=4, color="#2ca02c")
axe.axhline(0, color="k", lw=0.8)
axe.axvline(1.0 - delta / a, color="#d62728", ls="--", lw=1.3)
axe.set_xlabel("$r/a$")
axe.set_ylabel("error  %")
axe.grid(alpha=0.3)
axe.set_title("solver minus exact, Bessel evaluated at the radii actually "
              "sampled", fontsize=9)

fig.tight_layout()
fig.savefig(os.path.join(FIG, "j_radial_profile.png"), dpi=150)
print("wrote %s/j_radial_profile.png  -- worst binned error %.2f %%"
      % (FIG, np.abs(err).max()))

# --- R_ac/R_dc against Kelvin across frequency -----------------------------
# Drawn from the exact formula, with our two measured points on it, so the
# figure shows WHERE this case sits on the curve rather than just its value.
def ber_bei(x):
    b = x / 2.0
    rs = i = 0.0
    tr = 1.0
    for kk in range(0, 300):
        if kk:
            tr *= -(b ** 4) / ((2 * kk - 1) * (2 * kk)) ** 2
        rs += tr
        if abs(tr) < 1e-20 * max(abs(rs), 1e-300):
            break
    ti = b * b
    for kk in range(0, 300):
        if kk:
            ti *= -(b ** 4) / ((2 * kk) * (2 * kk + 1)) ** 2
        i += ti
        if abs(ti) < 1e-20 * max(abs(i), 1e-300):
            break
    return rs, i


def ber_bei_p(x):
    b = x / 2.0
    rp = 0.0
    t = -4.0 * b ** 3 * 0.5 / 4.0
    for kk in range(1, 300):
        if kk > 1:
            t *= -(b ** 4) * (4.0 * kk) / ((4.0 * kk - 4) * ((2 * kk - 1) * (2 * kk)) ** 2)
        rp += t
        if abs(t) < 1e-20 * max(abs(rp), 1e-300):
            break
    ip = 0.0
    t = 2 * b * 0.5
    for kk in range(0, 300):
        if kk:
            t *= -(b ** 4) * (4.0 * kk + 2) / ((4.0 * kk - 2) * ((2 * kk) * (2 * kk + 1)) ** 2)
        ip += t
        if abs(t) < 1e-20 * max(abs(ip), 1e-300):
            break
    return rp, ip


def kelvin(u):
    br, bi = ber_bei(u)
    bp, ip = ber_bei_p(u)
    return (u / 2) * (br * ip - bi * bp) / (bp ** 2 + ip ** 2)


ads = np.linspace(0.05, 4.0, 300)
kv = [kelvin(math.sqrt(2) * x) for x in ads]
fig2, ax2 = plt.subplots(figsize=(8.2, 5.2))
ax2.plot(ads, kv, "k-", lw=2, label="exact Kelvin  $R_{ac}/R_{dc}$")
pts = [(1.070, 1.02633, "01, 50 Hz"),
       (2.002, 1.26476, "01's mesh, 175 Hz"),
       (2.621, 1.57436, "01's mesh, 300 Hz"),
       (3.000, 1.77191, "04, 393 Hz")]
for x, y, lab in pts:
    ax2.plot([x], [y], "o", ms=9, mfc="none", mew=2,
             color="#d62728" if "04" in lab else "#1f77b4")
    ax2.annotate(lab, (x, y), textcoords="offset points", xytext=(9, -11),
                 fontsize=9)
ax2.set_xlabel(r"$a/\delta$")
ax2.set_ylabel(r"$R_{ac}/R_{dc}$")
ax2.set_title("Where each run sits on the exact curve.\n"
              "At a/$\\delta$=1.07 the effect is 2.6 % of R; at 3.00 it is 77 %.")
ax2.grid(alpha=0.3)
ax2.legend(loc="upper left")
fig2.tight_layout()
fig2.savefig(os.path.join(FIG, "kelvin_curve.png"), dpi=150)
print("wrote %s/kelvin_curve.png" % FIG)
