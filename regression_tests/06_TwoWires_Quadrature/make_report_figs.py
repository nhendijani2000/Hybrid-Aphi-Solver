"""Geometry, mesh and analysis figures for the case 06 report.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_report_figs.py [output-dir]

Writes into fig/. The post-processing pictures are NOT made here -- they come
from the solver's own manifest through tools/postprocess.py, which is the point
of the case; `collect_postprocess_figs.py` copies them in afterwards.

Figures that are only ever a picture (geometry, mesh) are rendered here, and the
two analysis plots (the skin profile and the mesh convergence) are matplotlib,
because they are graphs of numbers rather than views of a field.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

paraview.simple._DisableFirstRenderCameraReset()

OUT = sys.argv[1] if len(sys.argv) > 1 else "output"
FIG = "fig"
if not os.path.isdir(FIG):
    os.makedirs(FIG)

A, D, LZ, W, H = 2.0e-3, 10.0e-3, 20.0e-3, 30.0e-3, 30.0e-3
SIG, F = 5.8e7, 2500.0
MU0 = 4e-7 * math.pi
DELTA = math.sqrt(2.0 / (2 * math.pi * F * MU0 * SIG))
K = (1 - 1j) / DELTA
INK = "#1b1b1b"

view = CreateView("RenderView")
view.ViewSize = [1100, 900]
view.OrientationAxesVisibility = 0
view.Background = [1, 1, 1]
view.UseColorPaletteForBackground = 0


def clear():
    for s in GetSources().values():
        Hide(s, view)


def shot(name):
    Render()
    SaveScreenshot(os.path.join(FIG, name), view, ImageResolution=view.ViewSize)
    print("  wrote fig/%s" % name)


def src():
    return LegacyVTKReader(FileNames=[os.path.join(OUT, "J_field.vtk")])


def body(tag):
    t = Threshold(Input=src())
    t.Scalars = ["CELLS", "body_tag"]
    t.LowerThreshold = t.UpperThreshold = float(tag)
    t.ThresholdMethod = "Between"
    return t


# ---------------------------------------------------------------- geometry
def geometry(name, cut):
    """The two wires inside the air box, the box drawn transparent."""
    clear()
    view.CameraParallelProjection = 1
    air = body(3)
    if cut:
        c = Clip(Input=air)
        c.ClipType = "Plane"
        c.ClipType.Origin = [0, 0, 0]
        c.ClipType.Normal = [0, 1, 0]
        air = c
    da = Show(air, view)
    ColorBy(da, None)          # or it colours by whatever array is first
    da.Representation = "Surface"
    da.Opacity = 0.10
    da.AmbientColor = [0.2, 0.4, 0.8]
    da.DiffuseColor = [0.2, 0.4, 0.8]
    for tag, col in ((1, [0.85, 0.45, 0.15]), (2, [0.15, 0.55, 0.75])):
        w = body(tag)
        if cut:
            c = Clip(Input=w)
            c.ClipType = "Plane"
            c.ClipType.Origin = [0, 0, 0]
            c.ClipType.Normal = [0, 1, 0]
            w = c
        d = Show(w, view)
        ColorBy(d, None)       # solid colour per body, so the two are telling apart
        d.Representation = "Surface"
        d.DiffuseColor = col
        d.AmbientColor = col
    view.CameraPosition = [1.6 * W, -2.2 * W, 1.5 * LZ]
    view.CameraFocalPoint = [0, 0, 0.5 * LZ]
    view.CameraViewUp = [0, 0, 1]
    ResetCamera(view)
    shot(name)


geometry("geometry.png", False)
geometry("geometry_cut.png", True)


# -------------------------------------------------------------------- mesh
def mesh_shot(name, scale, only_wire=False):
    clear()
    view.CameraParallelProjection = 1
    base = body(1) if only_wire else src()
    sl = Slice(Input=base)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = [0, 0, 0.5 * LZ]
    sl.SliceType.Normal = [0, 0, 1]
    UpdatePipeline(proxy=sl)
    d = Show(sl, view)
    d.Representation = "Surface With Edges"
    ColorBy(d, None)
    d.AmbientColor = [0.80, 0.82, 0.86]
    d.DiffuseColor = [0.80, 0.82, 0.86]
    d.EdgeColor = [0.15, 0.17, 0.22]
    d.LineWidth = 1.0
    # centre on the wire when that is what is being shown -- wire 1 is at
    # x = -d/2, not at the origin, and a camera on the origin cuts it in half
    cx = -0.5 * D if only_wire else 0.0
    view.CameraPosition = [cx, 0, 0.5 * LZ + 1.0]
    view.CameraFocalPoint = [cx, 0, 0.5 * LZ]
    view.CameraViewUp = [0, 1, 0]
    view.CameraParallelScale = scale
    shot(name)


mesh_shot("mesh_domain.png", 0.52 * W)
mesh_shot("mesh_wire.png", 1.35 * A, only_wire=True)


# ------------------------------------------------- the skin profile, plotted
def bessel_j0(z):
    t, s = 1.0 + 0j, 1.0 + 0j
    for m in range(1, 80):
        t *= -(z * z) / (4.0 * m * m)
        s += t
        if abs(t) < 1e-18 * abs(s):
            break
    return s


s = LegacyVTKReader(FileNames=[os.path.join(OUT, "J_field.vtk")])
sz = CellSize(Input=s)
sz.ComputeVolume = 1
cc = CellCenters(Input=sz)
cc.VertexCells = 1
UpdatePipeline(proxy=cc)
dd = sm.Fetch(cc)
q = vtk_to_numpy(dd.GetPoints().GetData())
vol = np.abs(vtk_to_numpy(dd.GetPointData().GetArray("Volume")))
bt = vtk_to_numpy(dd.GetPointData().GetArray("body_tag"))
Jz = (vtk_to_numpy(dd.GetPointData().GetArray("J_cell_real"))[:, 2]
      + 1j * vtk_to_numpy(dd.GetPointData().GetArray("J_cell_imag"))[:, 2])
r = np.hypot(q[:, 0] + 0.5 * D, q[:, 1])
mid = (q[:, 2] > 0.35 * LZ) & (q[:, 2] < 0.65 * LZ)
w1 = mid & (bt == 1)
area = vol[w1].sum() / (0.30 * LZ)
a_eff = math.sqrt(area / math.pi)
surf = w1 & (r > 0.92 * a_eff)
Ja = Jz[surf].mean()
bref = np.mean([bessel_j0(K * x) for x in r[surf]])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.6, 7.0), sharex=True)
rr = np.linspace(0.02, 1.0, 300) * a_eff
exact = np.array([bessel_j0(K * x) for x in rr]) / bref
ax1.plot(rr / a_eff, np.abs(exact), "-", lw=2, color="k", zorder=5,
         label=r"exact  $J_0(kr)/J_0(k\bar r)$")
ax2.plot(rr / a_eff, np.degrees(np.angle(exact)), "-", lw=2, color="k", zorder=5)
step = max(1, int(w1.sum() // 4000))
sub = np.where(w1)[0][::step]
ax1.plot(r[sub] / a_eff, np.abs(Jz[sub] / Ja), ".", ms=2, alpha=0.30,
         color="#1f5fa9", label="per tet (sub-sampled)")
ax2.plot(r[sub] / a_eff, np.degrees(np.angle(Jz[sub] / Ja)), ".", ms=2, alpha=0.30,
         color="#1f5fa9")
edges = np.linspace(0.0, 0.96, 17)
for ax, fn in ((ax1, lambda v: np.abs(v)), (ax2, lambda v: np.degrees(np.angle(v)))):
    xs, ys = [], []
    for i in range(len(edges) - 1):
        m = w1 & (r >= edges[i] * a_eff) & (r < edges[i + 1] * a_eff)
        if m.sum() < 20:
            continue
        xs.append(0.5 * (edges[i] + edges[i + 1]))
        ys.append(fn(Jz[m].mean() / Ja))
    ax.plot(xs, ys, "o", ms=6, color="#c0522a", zorder=6, label="binned mean")
ax1.set_ylabel(r"$|J_z| \;/\; |J_z|$ on the $r>0.92a$ band")
ax2.set_ylabel("phase relative to that band  (deg)")
ax2.set_xlabel(r"$r/a$")
ax1.legend(fontsize=9, frameon=False, loc="upper left")
ax1.set_title("Case 06: the radial current profile against Kelvin"
              "\n$a/\\delta$ = %.3f, %d tets in the mid-length band of wire 1"
              % (a_eff / DELTA, int(w1.sum())), fontsize=11, color=INK)
for ax in (ax1, ax2):
    ax.grid(alpha=0.25, lw=0.6)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "skin_profile.png"), dpi=150, facecolor="white")
plt.close(fig)
print("  wrote fig/skin_profile.png")


# ------------------------------------------- convergence, if the data is there
CONV = "convergence.json"
if os.path.isfile(CONV):
    c = json.load(open(CONV))
    h = np.array(c["h_mm"])
    em = np.abs(c["mag_err_pct"])
    el = np.abs(c["lag_err_deg"])
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.loglog(DELTA * 1e3 / h, em, "o-", lw=1.6, ms=7, color="#1f5fa9",
              label="core $|J|$ error, %")
    ax.loglog(DELTA * 1e3 / h, el, "s-", lw=1.6, ms=7, color="#c0522a",
              label="core phase error, deg")
    for x, y, t in zip(DELTA * 1e3 / h, em, c["tets"]):
        ax.annotate("%d tets" % t, (x, y), fontsize=8, color="#666",
                    xytext=(4, 6), textcoords="offset points")
    ax.axvline(DELTA * 1e3 / h[-1], color="#2b7a4b", lw=1.0, ls="--")
    ax.annotate("shipped", (DELTA * 1e3 / h[-1], em[0]), fontsize=9,
                color="#2b7a4b", xytext=(-52, 0), textcoords="offset points")
    ax.set_xlabel(r"elements per skin depth, $\delta/h$")
    ax.set_ylabel("error against Kelvin, at the core")
    ax.set_title("Case 06: what refinement actually buys", fontsize=11, color=INK)
    ax.grid(alpha=0.25, lw=0.6, which="both")
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "convergence.png"), dpi=150, facecolor="white")
    plt.close(fig)
    print("  wrote fig/convergence.png")
else:
    print("  (no %s -- run make_convergence_data.py for fig/convergence.png)" % CONV)
