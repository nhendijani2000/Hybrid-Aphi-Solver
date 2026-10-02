"""Assert that the tree-cotree gauge moves the potentials and nothing else.

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvbatch.exe" verify_gauge.py

Run by ../check.bat, which finds it as this case's own verify_*.py. It performs
its own solves -- five of them -- so it works in `verify-only` mode too and does
not depend on which .aphi check.bat happened to run first.

WHAT IS BEING TESTED. The A-Phi formulation is invariant under

    A -> A + grad(psi)        Phi -> Phi - j*omega*psi

so B = curl A is unchanged (curl grad = 0) and E = -j*omega*A - grad(Phi) is
unchanged (the two shifts cancel). The tree-cotree gauge removes the freedom by
forcing A = 0 on a spanning tree's edges, which picks ONE representative -- so a
different tree must give the same fields and different potentials.

This is exact in the discrete setting, not only the continuum: the gradient of a
nodal function is exactly representable in the Whitney edge space (on edge (i,j)
its DOF is psi_j - psi_i), so the discrete gauge freedom is exactly the discrete
gradients and discrete curl-grad is exactly zero. The tolerances below are
round-off, not discretisation.

BOTH DIRECTIONS ARE ASSERTED, and that is deliberate. A test that only checked
invariance would still pass if the gauge freedom were accidentally REMOVED --
by over-constraining Phi, say -- which would be a real regression sitting behind
a green check. So the potentials are required to move by at least a floor, and
the fields to move by at most a ceiling, with eleven orders of magnitude of
daylight between the two thresholds.

HOW DIFFERENT TREES ARE OBTAINED: by permuting the order the nodes are listed in
the .msh. See permute_nodes.py for why that changes the tree while leaving the
problem identical.
"""
import math
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
sys.path.insert(0, HERE)
from permute_nodes import permute          # noqa: E402

MU0 = 4e-7 * math.pi
RUN_CASE = os.path.join(os.path.dirname(HERE), "run_case.bat")

# The five solves: (label, source .aphi, mesh, output dir).
VOLT = [("v_base", "gauge.aphi", "gauge.msh", "g_v_base"),
        ("v_permA", "gauge.aphi", "g_permA.msh", "g_v_permA"),
        ("v_permB", "gauge.aphi", "g_permB.msh", "g_v_permB")]
CURR = [("i_base", "gauge_1a.aphi", "gauge.msh", "g_i_base"),
        ("i_permA", "gauge_1a.aphi", "g_permA.msh", "g_i_permA")]
SEEDS = {"g_permA.msh": 11111, "g_permB.msh": 99999}

GENERATED = (["g_permA.msh", "g_permB.msh"]
             + ["g_%s.aphi" % t[0] for t in VOLT + CURR]
             + [t[3] for t in VOLT + CURR])

failures = []
notes = []


def expected():
    """Thresholds from expected.txt -- the same format the other cases use."""
    exp = {}
    with open("expected.txt") as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    exp[parts[0]] = float(parts[1])
                except ValueError:
                    exp[parts[0]] = parts[1]
    return exp


def check(name, got, lo=None, hi=None, fmt="%.6g"):
    ok = True
    bound = ""
    if lo is not None:
        ok = ok and got >= lo
        bound += " >= " + fmt % lo
    if hi is not None:
        ok = ok and got <= hi
        bound += " <= " + fmt % hi
    print("  %-46s %-14s %-24s %s"
          % (name, fmt % got, bound.strip(), "ok" if ok else "FAIL"))
    if not ok:
        failures.append(name)


def cleanup():
    for p in GENERATED:
        if os.path.isdir(p):
            shutil.rmtree(p, ignore_errors=True)
        elif os.path.isfile(p):
            try:
                os.remove(p)
            except OSError:
                pass


def solve_all():
    """Mesh if needed, make the permuted meshes, run the five solves."""
    if not os.path.isfile("gauge.msh"):
        raise SystemExit("gauge.msh missing -- run ../run_case.bat gauge.aphi once "
                         "to mesh it, or mesh gauge.geo by hand")
    for dst, seed in SEEDS.items():
        n, moved = permute("gauge.msh", dst, seed)
        notes.append("%s: %d of %d nodes moved (%.1f %%)"
                     % (dst, moved, n, 100.0 * moved / n))

    for label, src, msh, out in VOLT + CURR:
        aphi = "g_%s.aphi" % label
        text = open(src).read()
        text = re.sub(r"^file\s*=.*$", "file        = " + msh, text,
                      count=1, flags=re.M)
        text = re.sub(r"^directory\s*=.*$", "directory   = " + out, text,
                      count=1, flags=re.M)
        with open(aphi, "w", newline="\n") as f:
            f.write(text)
        r = subprocess.run(["cmd", "/c", RUN_CASE, aphi],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(r.stdout.decode("utf-8", "replace")[-2000:])
            raise SystemExit("solve failed for %s" % aphi)
        nz = re.search(rb"unknowns\s+(\d+)\s+(\d+) stored nonzeros", r.stdout)
        if nz:
            notes.append("%-8s %s unknowns, %s stored nonzeros"
                         % (label, nz.group(1).decode(), nz.group(2).decode()))


def load(out, fname, ncomp):
    pos, val = [], []
    with open(os.path.join(out, fname)) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            p = line.split()
            pos.append([float(p[1]), float(p[2]), float(p[3])])
            c = [float(x) for x in p[4:4 + 2 * ncomp]]
            val.append([complex(c[2 * k], c[2 * k + 1]) for k in range(ncomp)])
    return np.array(pos), np.array(val)


def pair(out_a, out_b, fname, ncomp):
    """Match two solves node by node on coordinates.

    The node ORDER differs between runs -- that is the whole point -- so they
    cannot be compared index by index. Coordinates are written from the same
    doubles in every run, so rounding at 1e-12 m (a picometre against a mesh in
    millimetres) is unambiguous.
    """
    pa, va = load(out_a, fname, ncomp)
    pb, vb = load(out_b, fname, ncomp)
    ib = {tuple(np.round(p, 12)): i for i, p in enumerate(pb)}
    ia, jb = [], []
    for i, p in enumerate(pa):
        j = ib.get(tuple(np.round(p, 12)))
        if j is not None:
            ia.append(i)
            jb.append(j)
    return pa[ia], va[ia], vb[jb]


# ---------------------------------------------------------------------------
exp = expected()
print("verifying %s" % HERE)
print("  one problem, three spanning trees: Phi and A must move, fields must not")
print()
solve_all()
for n in notes:
    print("  %s" % n)
print()

# --- the regime. Phi only moves appreciably away from the resistive limit, so
# --- a pass here means nothing unless omega*L/R is actually large. Measured
# --- from the current-driven solve, where Phi on the driven face IS V.
omega = 2.0 * math.pi * 50.0
pos, phi = load("g_i_base", "potential.out", 1)
r = np.hypot(pos[:, 0], pos[:, 1])
z = pos[:, 2]
A_M, L_M = 10.0e-3, 40.0e-3
drv = (r <= A_M * 1.0001) & (np.abs(z - L_M) < 1e-9)
ref = (r <= A_M * 1.0001) & (np.abs(z) < 1e-9)
Z = (phi[drv, 0].mean() - phi[ref, 0].mean()) / 1.0
wLR = abs(Z.imag) / abs(Z.real)
print("regime")
check("omega*L/R (must be gauge dominated)", wLR, lo=exp["wLR_min"])
print()

# --- the potentials must MOVE ----------------------------------------------
print("the potentials are gauge dependent -- they must MOVE")
for nm, fname, nc in (("Phi", "potential.out", 1), ("A", "A_field.out", 3)):
    for other in ("g_v_permA", "g_v_permB"):
        _, a, b = pair("g_v_base", other, fname, nc)
        rel = np.abs(b - a).max() / np.abs(a).max()
        check("%s, base vs %s, relative change" % (nm, other[4:]), rel,
              lo=exp["potential_moves_min"])
print()

# --- the fields must NOT ----------------------------------------------------
print("the fields are gauge invariant -- they must NOT move")
for nm, fname in (("E", "E_field.out"), ("B", "B_field.out"),
                  ("H", "H_field.out"), ("J", "J_field.out")):
    for other in ("g_v_permA", "g_v_permB"):
        _, a, b = pair("g_v_base", other, fname, 3)
        rel = np.abs(b - a).max() / np.abs(a).max()
        check("%s, base vs %s, relative change" % (nm, other[4:]), rel,
              hi=exp["field_invariance_max"], fmt="%.3e")
print()

# --- psi vanishes on the outer boundary ------------------------------------
# Recovered from the shift itself: delta_Phi = -j*omega*delta_psi.
print("the gauge function psi lives in the INTERIOR and vanishes on the boundary")
pos, a, b = pair("g_v_base", "g_v_permA", "potential.out", 1)
dpsi = np.abs(1j * (b[:, 0] - a[:, 0]) / omega)
h, tol = 100.0e-3, 1e-9
x, y, zc = pos[:, 0], pos[:, 1], pos[:, 2]
outer = ((np.abs(np.abs(x) - h) < tol) | (np.abs(np.abs(y) - h) < tol)
         | (np.abs(zc) < tol) | (np.abs(zc - L_M) < tol))
check("max |psi| on the outer boundary (Wb)", dpsi[outer].max(),
      hi=exp["psi_boundary_max"], fmt="%.3e")
check("max |psi| in the interior (Wb)", dpsi[~outer].max(),
      lo=exp["psi_interior_min"], fmt="%.3e")
print()

# --- the terminals, under both drives ---------------------------------------
print("the terminals are exact in every gauge, under EITHER drive")
for out in ("g_v_base", "g_v_permA", "g_v_permB"):
    pos, v = load(out, "potential.out", 1)
    r = np.hypot(pos[:, 0], pos[:, 1])
    w = r <= A_M * 1.0001
    bot = w & (np.abs(pos[:, 2]) < 1e-9)
    top = w & (np.abs(pos[:, 2] - L_M) < 1e-9)
    check("voltage %s: max |Phi - 0| on reference" % out[4:],
          np.abs(v[bot, 0]).max(), hi=exp["terminal_exact_max"], fmt="%.3e")
    check("voltage %s: max |Phi - 1| on driven" % out[4:],
          np.abs(v[top, 0] - 1.0).max(), hi=exp["terminal_exact_max"], fmt="%.3e")

vals = []
for out in ("g_i_base", "g_i_permA"):
    pos, v = load(out, "potential.out", 1)
    r = np.hypot(pos[:, 0], pos[:, 1])
    w = r <= A_M * 1.0001
    top = w & (np.abs(pos[:, 2] - L_M) < 1e-9)
    bot = w & (np.abs(pos[:, 2]) < 1e-9)
    vals.append(v[top, 0].mean() - v[bot, 0].mean())
    check("current %s: equipotential spread, driven face" % out[4:],
          np.abs(v[top, 0] - v[top, 0].mean()).max(),
          hi=exp["terminal_exact_max"], fmt="%.3e")
dV = abs(vals[1] - vals[0]) / abs(vals[0])
print("    V = %.6e V (base), %.6e V (permA)" % (abs(vals[0]), abs(vals[1])))
check("current: relative change in terminal V", dV,
      hi=exp["field_invariance_max"], fmt="%.3e")

print()
cleanup()
if failures:
    print("GAUGE CHECK FAILED: %s" % ", ".join(failures))
    sys.exit(1)
print("all gauge checks passed")
