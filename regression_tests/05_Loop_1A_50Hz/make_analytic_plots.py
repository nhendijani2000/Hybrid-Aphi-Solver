"""The analytic comparison figures for section 5 of LoopValidation.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" make_analytic_plots.py

Writes fig/analytic_E.png and fig/analytic_B.png. Needs output/ from
loop_50hz.aphi; uses output_big/ from loop_50hz_big.aphi as well when it is
there, which is what shows the outer boundary is the far-field error.

The closed forms and where they come from are derived in analytic_comparison.py's
docstring; this only plots them. Every number the report quotes comes from that
script, not from here.
"""
from paraview.simple import *
from paraview import servermanager as sm
from paraview.numpy_support import vtk_to_numpy
import numpy as np
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R, a, sigma, I = 6.5e-3, 0.8e-3, 5.8e7, 1.0
mu0 = 4e-7 * math.pi
FIG = "fig"
INK, GRID = "#1b1b1b", "#d8d8d8"
C_EX, C_S, C_B = "#000000", "#1f5fa9", "#c0522a"


def KE(m):
    """K(k), E(k) with m = k^2, by the AGM. See analytic_comparison.py."""
    aa, bb, cc, acc, p = 1.0, math.sqrt(1.0 - m), math.sqrt(m), 0.5 * m, 0.5
    for _ in range(60):
        if abs(cc) < 1e-17:
            break
        aa, bb, cc = 0.5 * (aa + bb), math.sqrt(aa * bb), 0.5 * (aa - bb)
        p *= 2.0
        acc += p * cc * cc
    K = math.pi / (2.0 * aa)
    return K, K * (1.0 - acc)


def B_filament(rho, z):
    """|B| of a filamentary loop of radius R carrying I, at (rho, z)."""
    Q = (R + rho) ** 2 + z * z
    D = (R - rho) ** 2 + z * z
    m = 4.0 * R * rho / Q
    out = np.empty(np.size(rho))
    for i in range(np.size(rho)):
        K, E = KE(min(m[i], 1.0 - 1e-15))
        pref = mu0 * I / (2.0 * math.pi) / math.sqrt(Q[i])
        bz = pref * (K + E * (R * R - rho[i] ** 2 - z[i] ** 2) / D[i])
        br = 0.0
        if rho[i] > 1e-12:
            br = pref * (z[i] / rho[i]) * (
                -K + E * (R * R + rho[i] ** 2 + z[i] ** 2) / D[i])
        out[i] = math.hypot(bz, br)
    return out


def load(out, fn, arrs):
    """Cell-centre positions, volumes, body tags and |field| for one file."""
    s = LegacyVTKReader(FileNames=[os.path.join(out, fn)])
    sz = CellSize(Input=s)
    sz.ComputeVolume = 1
    cc = CellCenters(Input=sz)
    cc.VertexCells = 1
    UpdatePipeline(proxy=cc)
    d = sm.Fetch(cc)
    q = vtk_to_numpy(d.GetPoints().GetData())
    vol = np.abs(vtk_to_numpy(d.GetPointData().GetArray("Volume")))
    bt = vtk_to_numpy(d.GetPointData().GetArray("body_tag"))
    re = vtk_to_numpy(d.GetPointData().GetArray(arrs[0]))
    im = vtk_to_numpy(d.GetPointData().GetArray(arrs[1]))
    return q, vol, bt, np.sqrt((re ** 2).sum(1) + (im ** 2).sum(1))


def axis_probe(out, zs):
    """|B| on the loop axis, and the exact superposition over the mesh."""
    import vtk
    s = LegacyVTKReader(FileNames=[os.path.join(out, "B_field.vtk")])
    UpdatePipeline(proxy=s)
    data = sm.Fetch(s)
    pts = vtk.vtkPoints()
    for z0 in zs:
        pts.InsertNextPoint(0.0, 0.0, z0)
    poly = vtk.vtkPolyData()
    poly.SetPoints(pts)
    pr = vtk.vtkProbeFilter()
    pr.SetInputData(poly)
    pr.SetSourceData(data)
    pr.Update()
    o = pr.GetOutput()
    br = vtk_to_numpy(o.GetPointData().GetArray("B_real"))
    bi = vtk_to_numpy(o.GetPointData().GetArray("B_imag"))
    return np.sqrt((br ** 2).sum(1) + (bi ** 2).sum(1))


HAVE_BIG = os.path.isdir("output_big") and os.path.isfile("output_big/B_field.vtk")

# ===========================================================================
# figure 1 -- E inside the conductor against the exact 1/rho law
# ===========================================================================
q, vol, bt, E = load("output", "E_field.vtk", ("E_cell_real", "E_cell_imag"))
x, y, z = q[:, 0], q[:, 1], q[:, 2]
rho = np.hypot(x, y)
cond = bt == 1
K = (vol[cond] / rho[cond] ** 2).sum() / (2 * math.pi)
th = np.arctan2(y, x) % (2 * math.pi)

fig, (A, B) = plt.subplots(1, 2, figsize=(11.0, 4.3))

rr = np.linspace((R - a) * 0.995, (R + a) * 1.005, 400)
A.plot(rr * 1e3, I / (sigma * rr * K) * 1e3, color=C_EX, lw=2.0, zorder=5,
       label=r"exact  $E=I/(\sigma\rho K)$")
edges = np.linspace(R - a, R + a, 22)
mid = 0.5 * (edges[:-1] + edges[1:])
med = np.array([np.median(E[cond & (rho >= edges[i]) & (rho < edges[i + 1])])
                if (cond & (rho >= edges[i]) & (rho < edges[i + 1])).sum() > 20
                else np.nan for i in range(len(mid))])
A.plot(mid * 1e3, med * 1e3, "o", ms=5.5, color=C_S, zorder=6,
       label="solver, median per band")
A.set_ylim(7.4, 10.35)
A.axvline((R - a) * 1e3, color=GRID, lw=1.0)
A.axvline((R + a) * 1e3, color=GRID, lw=1.0)
# labels along the BOTTOM: the curve and the legend both live up top
A.annotate("inner wall $\\rho=R-a$", ((R - a) * 1e3, 7.47), fontsize=8,
           color="#666", ha="left", va="bottom", xytext=(4, 0),
           textcoords="offset points")
A.annotate("outer wall $\\rho=R+a$", ((R + a) * 1e3, 7.47), fontsize=8,
           color="#666", ha="right", va="bottom", xytext=(-4, 0),
           textcoords="offset points")
A.annotate("ratio across the tube  $(R+a)/(R-a) = 1.2807$",
           (0.46, 0.155), xycoords="axes fraction", fontsize=8.5,
           color="#666", ha="center")
A.set_xlabel(r"cylindrical radius  $\rho$  (mm)")
A.set_ylabel(r"$|E|$  (mV/m)")
A.set_title(r"(a)  $E$ across the tube: it goes as $1/\rho$, not uniform",
            fontsize=10.5, color=INK)
A.legend(fontsize=8.5, frameon=False, loc="upper right")
A.grid(alpha=0.25, lw=0.6)

# panel b -- the residual, by azimuth, which is where the port shows up
for lab, sel, col, mk in (
        ("generic azimuth, $\\theta=\\pi/2$",
         np.abs(((th - math.pi / 2 + math.pi) % (2 * math.pi)) - math.pi) < 0.25,
         "#2b7a4b", "s"),
        ("the joint, $\\theta=\\pi$", (np.abs(y) < 0.35e-3) & (x < 0), C_S, "o"),
        ("THE PORT, $\\theta=0$", (np.abs(y) < 0.35e-3) & (x > 0), C_B, "^")):
    pm, pe = [], []
    for i in range(len(mid)):
        m = sel & cond & (rho >= edges[i]) & (rho < edges[i + 1])
        if m.sum() < 8:
            continue
        pm.append(mid[i] * 1e3)
        pe.append(100.0 * np.median(E[m] / (I / (sigma * rho[m] * K)) - 1.0))
    B.plot(pm, pe, mk + "-", ms=4.5, lw=1.2, color=col, label=lab)
B.axhline(0.0, color=C_EX, lw=1.2)
B.set_xlabel(r"cylindrical radius  $\rho$  (mm)")
B.set_ylabel("error against the exact law  (%)")
B.set_title("(b)  the same residual, split by azimuth", fontsize=10.5, color=INK)
B.legend(fontsize=8.5, frameon=False)
B.grid(alpha=0.25, lw=0.6)
fig.suptitle("Case 05: $E$ inside the conductor against the exact interior solution"
             "  ($K=\\int dA/\\rho$ taken from the mesh)", fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(FIG, "analytic_E.png"), dpi=150,
            facecolor="white", bbox_inches="tight")
plt.close(fig)
print("  wrote fig/analytic_E.png")

# ===========================================================================
# figure 2 -- B, on the axis and through the air
# ===========================================================================
fig, (A, B) = plt.subplots(1, 2, figsize=(11.0, 4.3))

zs = np.array([0.0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8.5, 10]) * 1e-3
qb, volb, btb, _ = load("output", "B_field.vtk", ("B_cell_real", "B_cell_imag"))
cb = btb == 1
rhob = np.hypot(qb[:, 0], qb[:, 1])
exact = np.array([(mu0 * I / (4 * math.pi * K)) *
                  (volb[cb] / (rhob[cb] ** 2 + (z0 - qb[cb, 2]) ** 2) ** 1.5).sum()
                  for z0 in zs])
A.plot(zs * 1e3, exact * 1e6, color=C_EX, lw=2.0, zorder=5,
       label="exact, superposed over the real\ncurrent distribution")
A.plot(zs * 1e3, axis_probe("output", zs) * 1e6, "o--", ms=5, lw=1.1,
       color=C_B, label="solver, $R_d=H_d=30$ mm")
if HAVE_BIG:
    A.plot(zs * 1e3, axis_probe("output_big", zs) * 1e6, "s-", ms=5, lw=1.1,
           color=C_S, label="solver, $R_d=H_d=90$ mm")
A.set_xlabel("height on the loop axis  $z$  (mm)")
A.set_ylabel(r"$|B|$  ($\mu$T)")
A.set_title("(a)  $B$ on the loop axis", fontsize=10.5, color=INK)
A.legend(fontsize=8.5, frameon=False)
A.grid(alpha=0.25, lw=0.6)

# the two curves above differ by a couple of per cent, which is invisible at
# this scale and is the whole point -- so plot the residual as an inset
ins = A.inset_axes([0.08, 0.09, 0.38, 0.28])
ins.set_facecolor("white")          # opaque, or the main curve shows through
ins.set_zorder(10)
ins.patch.set_alpha(1.0)
for sp in ins.spines.values():
    sp.set_edgecolor("#999")
ins.axhline(0.0, color=C_EX, lw=1.0)
ins.plot(zs * 1e3, 100 * (axis_probe("output", zs) / exact - 1), "o-", ms=3,
         lw=1.0, color=C_B)
if HAVE_BIG:
    ins.plot(zs * 1e3, 100 * (axis_probe("output_big", zs) / exact - 1), "s-",
             ms=3, lw=1.0, color=C_S)
ins.set_ylim(-24, 8)
ins.set_title("error, %", fontsize=8, color="#555", pad=2)
ins.tick_params(labelsize=7, length=2)
ins.grid(alpha=0.25, lw=0.5)

bands = ((2, 3), (3, 5), (5, 8), (8, 12), (12, 20), (20, 35))
for lab, out, col, mk in (("$R_d=H_d=30$ mm", "output", C_B, "o"),
                          ("$R_d=H_d=90$ mm", "output_big", C_S, "s")):
    if out == "output_big" and not HAVE_BIG:
        continue
    qq, vv, tt, BB = load(out, "B_field.vtk", ("B_cell_real", "B_cell_imag"))
    rr2 = np.hypot(qq[:, 0], qq[:, 1])
    dd = np.hypot(rr2 - R, qq[:, 2])
    airm = tt == 2
    px, py = [], []
    for lo, hi in bands:
        m = airm & (dd >= lo * a) & (dd < hi * a)
        if m.sum() < 30:
            continue
        px.append(math.sqrt(lo * hi))
        py.append(100.0 * np.median(BB[m] / B_filament(rr2[m], qq[m, 2]) - 1.0))
    B.plot(px, py, mk + "-", ms=5.5, lw=1.3, color=col, label=lab)
B.axhline(0.0, color=C_EX, lw=1.2)
B.axhspan(-1, 1, color="#2b7a4b", alpha=0.09, lw=0)
B.annotate("$\\pm1$ %", (2.1, 1.0), fontsize=8, color="#2b7a4b", va="bottom")
B.set_xscale("log")
B.set_xlabel("distance from the tube axis  $d/a$")
B.set_ylabel("error against the filament closed form  (%)")
B.set_title("(b)  $B$ through the air, and the outer wall's reach",
            fontsize=10.5, color=INK)
B.legend(fontsize=8.5, frameon=False, loc="upper left")
B.grid(alpha=0.25, lw=0.6, which="both")
fig.suptitle("Case 05: $B$ against its two closed forms — the far-field error is the"
             " OUTER BOUNDARY, not the mesh", fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(FIG, "analytic_B.png"), dpi=150,
            facecolor="white", bbox_inches="tight")
plt.close(fig)
print("  wrote fig/analytic_B.png")
