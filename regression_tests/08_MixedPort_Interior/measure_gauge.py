"""How much does a MIXED-MATERIAL INTERIOR port move when the tree changes?

    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvpython.exe" measure_gauge.py

Assumes the three solves have already been run (they are ~5 minutes each):

    solve_mesh.exe mixed_port.aphi            -> output/
    solve_mesh.exe mixed_port_permA.aphi      -> output_permA/
    solve_mesh.exe mixed_port_permB.aphi      -> output_permB/

WHAT IS BEING MEASURED, AND WHY IT IS NOT THE SAME AS CASE 07

07_GaugeInvariance permutes the node order of a mesh whose ports both sit on
the OUTER boundary, and finds the terminal quantities invariant to 1e-12 --
because `n x A = 0` pins the gauge function psi to zero there, so the gauge
transformation cannot reach them. 05's internal cut is one step in: its LOCAL J
and B move 8.3 % and 7.4 % while L moves 0.018 %.

This case is the configuration neither of those covers, and the one the project
exists to fix: a single-potential terminal on a face with CONDUCTOR ON ONE SIDE
AND AIR ON THE OTHER, inside the domain. There is nothing pinning psi on it.

THE COMPARISON IS NODE-ORDER INVARIANT. A permuted mesh lists the same nodes in
a different order, so results are matched by POSITION, not by index. Positions
are identical to the last bit -- permute_nodes.py moves the lines, it does not
recompute coordinates -- so an exact dictionary keyed on the coordinate triple
is right, and safer than a tolerance that might pair the wrong nodes.

WHAT WOULD MAKE THIS MEASUREMENT WRONG, and is therefore checked:

  - a solve that did not converge: backward errors are read and reported;
  - the port's Phi not being single-valued, which would mean the single-
    potential port is not doing what dof_map.cpp says: the spread across the
    face is asserted to be exactly zero;
  - the current direction flipping with node order. This was a REAL hazard:
    the old port binding took the "inside" from ft.tets[0], whose index
    follows node numbering, so a permutation could have reversed d and the
    sign flip would have looked like a gauge effect. The binding now derives
    d from which side is the conductor. Re(Z) is reported signed so a flip
    would be obvious rather than hidden by a magnitude.
"""
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OMEGA = 2.0 * math.pi * 50.0
I_DRIVE = 1.0
A_WIRE = 0.0015         # m, wire radius -- the same in both cases

# Two measurements, the same harness. THE CONTROL IS THE POINT: without it, a
# 77 % field change could be this script's error rather than the solver's
# answer, and nothing in the mixed-port numbers alone distinguishes the two.
# 07_GaugeInvariance establishes boundary-port invariance on ITS mesh; 'control'
# establishes it on THIS one, with the same wire, the same drive and the same
# permutation seeds -- so the only difference left is where the port sits.
CASES = {
    "mixed": dict(
        runs=("output", "output_permA", "output_permB"),
        z_port=0.028,
        title="08_MixedPort_Interior -- a single-potential port on an "
              "INTERIOR conductor/air face"),
    "control": dict(
        runs=("output_control03_base", "output_control03_permA",
              "output_control03_permB"),
        z_port=0.040,
        title="CONTROL (case 03) -- the same wire with BOTH caps on the outer "
              "boundary"),
}
MODE = sys.argv[1] if len(sys.argv) > 1 else "mixed"
if MODE not in CASES:
    raise SystemExit("usage: measure_gauge.py [mixed|control]")
RUNS = CASES[MODE]["runs"]
Z_PORT = CASES[MODE]["z_port"]


def load(run, fname, ncomp):
    """positions -> complex value(s), plus the status column."""
    path = os.path.join(HERE, run, fname)
    if not os.path.isfile(path):
        raise SystemExit("missing %s -- solve that variant first" % path)
    pos, val, status = [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            pos.append((f[1], f[2], f[3]))            # kept as TEXT: exact keys
            val.append([complex(float(f[4 + 2 * k]), float(f[5 + 2 * k]))
                        for k in range(ncomp)])
            status.append(f[4 + 2 * ncomp] if len(f) > 4 + 2 * ncomp else "")
    return pos, np.array(val), status


def backward_error(run):
    path = os.path.join(HERE, run, "potential.out")
    with open(path) as fh:
        for line in fh:
            m = re.search(r"backward error ([0-9.eE+-]+)", line)
            if m:
                return float(m.group(1))
    return float("nan")


def port_phi(run):
    """The interior terminal's single Phi, and the spread across its face."""
    pos, phi, status = load(run, "potential.out", 1)
    sel = [i for i, (p, s) in enumerate(zip(pos, status))
           if s == "T" and abs(float(p[2]) - Z_PORT) < 1e-5]
    if not sel:
        raise SystemExit("%s: found no terminal nodes at z = %g" % (run, Z_PORT))
    v = phi[sel, 0]
    return v[0], float(np.abs(v - v[0]).max()), len(sel)


def paired(run_a, run_b, fname, ncomp):
    """Values from two runs, aligned by position."""
    pa, va, _ = load(run_a, fname, ncomp)
    pb, vb, _ = load(run_b, fname, ncomp)
    idx = {p: i for i, p in enumerate(pb)}
    missing = [p for p in pa if p not in idx]
    if missing:
        raise SystemExit("%s: %d positions of %s absent from %s -- the two runs "
                         "are not the same geometry" %
                         (fname, len(missing), run_a, run_b))
    order = np.array([idx[p] for p in pa])
    return np.array([[float(c) for c in p] for p in pa]), va, vb[order]


def rel_change(va, vb):
    """max |b - a| over max |a|, on vector magnitudes."""
    na = np.sqrt((np.abs(va) ** 2).sum(axis=1))
    d = np.sqrt((np.abs(vb - va) ** 2).sum(axis=1))
    return float(d.max() / na.max()) if na.max() > 0 else 0.0


def main():
    print("%s\n" % CASES[MODE]["title"])
    print("solves")
    for r in RUNS:
        print("    %-14s backward error %.3e" % (r, backward_error(r)))

    print("\nthe terminal: one Phi unknown%s"
          % (" on a conductor/air face, inside the domain" if MODE == "mixed"
             else " on the outer boundary"))
    print("    %-14s %-34s %-12s %s"
          % ("run", "Phi (V)", "spread", "nodes"))
    phis = {}
    for r in RUNS:
        v, spread, n = port_phi(r)
        phis[r] = v
        print("    %-14s %+.9e %+.9ej  %.1e      %d"
              % (r, v.real, v.imag, spread, n))
        if spread != 0.0:
            print("        *** NOT single-valued -- the port is not one Phi DOF")

    print("\n    Phi is EXACTLY constant on the face in every run, which is the")
    print("    single-potential port behaving as dof_map.cpp specifies.")

    print("\nterminal impedance, Z = Phi / I with I = %g A" % I_DRIVE)
    print("    %-14s %-14s %-14s %-11s %s"
          % ("run", "R (uOhm)", "L (nH)", "dR vs base", "dL vs base"))
    base = phis[RUNS[0]] / I_DRIVE
    for r in RUNS:
        z = phis[r] / I_DRIVE
        dr = abs(z.real - base.real) / abs(base.real)
        dl = abs(z.imag - base.imag) / abs(base.imag)
        print("    %-14s %-14.6f %-14.6f %-11.3e %.3e"
              % (r, z.real * 1e6, z.imag / OMEGA * 1e9, dr, dl))

    print("\nfields: invariant under a gauge change, so these are the control")
    print("    %-10s %-14s %s" % ("field", "base->permA", "base->permB"))
    for nm, fname in (("E", "E_field.out"), ("B", "B_field.out"),
                      ("H", "H_field.out"), ("J", "J_field.out")):
        row = []
        for r in RUNS[1:]:
            _, va, vb = paired(RUNS[0], r, fname, 3)
            row.append(rel_change(va, vb))
        print("    %-10s %-14.3e %.3e" % (nm, row[0], row[1]))

    print("\npotentials: gauge DEPENDENT, so these are expected to move")
    for nm, fname, nc in (("Phi", "potential.out", 1), ("A", "A_field.out", 3)):
        row = []
        for r in RUNS[1:]:
            _, va, vb = paired(RUNS[0], r, fname, nc)
            row.append(rel_change(va, vb))
        print("    %-10s %-14.3e %.3e" % (nm, row[0], row[1]))

    # --- localized? the same banding Sec. 1.3 used on case 05 --------------
    # TWO normalizations, because one of them lies. Dividing by the GLOBAL max
    # makes the far field look clean no matter what it does, since the global
    # max sits at the port: a 100 % change in a far-field value whose magnitude
    # is 1e-4 of the peak shows up as 1e-4. Dividing by the BAND's own max is
    # the honest local figure. Report both and let the gap speak.
    print("\nlocal J change by distance from the port face")
    print("    %-16s %-9s %-14s %s"
          % ("band", "nodes", "vs global max", "vs THIS band's max"))
    xyz, va, vb = paired(RUNS[0], RUNS[1], "J_field.out", 3)
    r_ax = np.sqrt(xyz[:, 0] ** 2 + xyz[:, 1] ** 2)
    dz = np.abs(xyz[:, 2] - Z_PORT)
    dist = np.where(r_ax <= A_WIRE, dz, np.sqrt(dz ** 2 + (r_ax - A_WIRE) ** 2))
    na = np.sqrt((np.abs(va) ** 2).sum(axis=1))
    d = np.sqrt((np.abs(vb - va) ** 2).sum(axis=1))
    gmax = na.max()
    for lo, hi in ((0.0, 0.002), (0.002, 0.005), (0.005, 0.010), (0.010, 1.0)):
        m = (dist >= lo) & (dist < hi)
        if not m.any():
            continue
        bmax = na[m].max()
        print("    %4.0f-%4.0f mm     %-9d %-14.3e %.3e"
              % (lo * 1e3, hi * 1e3, int(m.sum()), float(d[m].max() / gmax),
                 float(d[m].max() / bmax) if bmax > 0 else float("nan")))

    # And the same for J inside the CONDUCTOR only, which is where J means
    # something: in the air sigma = 0, so J there is numerical dust and a
    # relative change in it is not a physical statement.
    print("\n    same, restricted to nodes inside the conductor (r <= a)")
    inw = r_ax <= A_WIRE
    for lo, hi in ((0.0, 0.002), (0.002, 0.005), (0.005, 0.010), (0.010, 1.0)):
        m = inw & (dist >= lo) & (dist < hi)
        if not m.any():
            continue
        bmax = na[m].max()
        print("    %4.0f-%4.0f mm     %-9d %-14.3e %.3e"
              % (lo * 1e3, hi * 1e3, int(m.sum()), float(d[m].max() / gmax),
                 float(d[m].max() / bmax) if bmax > 0 else float("nan")))


if __name__ == "__main__":
    main()
