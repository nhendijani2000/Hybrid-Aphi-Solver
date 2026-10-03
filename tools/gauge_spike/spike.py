"""Does the generalized-Lorenz gauge condition well enough to be worth building?

    gmsh tiny.geo -3 -o tiny.msh
    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvpython.exe" spike.py

A throwaway investigation kept because of what it found. It assembles the
matrices of `docs/GENERALIZED_LORENZ_GAUGE.md` §5 in Python on a few hundred
tets -- small enough for dense SVDs, so condition numbers are exact rather than
estimated -- and answers three questions without touching a line of C++:

  1. Is the §7 block system exactly complex-symmetric?
  2. Does the gauge condition as well as Chew reports?
  3. Does any of that survive finite conductivity?

**The answers are yes, yes, and no.** See README.md for the findings and
`GENERALIZED_LORENZ_GAUGE.md` §11 for what they mean.

CONVENTION: e^{+jwt}, so E = -jwA - grad(Phi). Chew writes e^{-iwt}; every i in
his paper is -j here. docs/GENERALIZED_LORENZ_GAUGE.md §0.

WHY THE CONTROLS MATTER. The first two runs of this spike produced condition
numbers of 1e38 and concluded the formulation was unusable. Both were wrong: the
first was my own block scaling, the second was finite conductivity, which Chew's
paper does not contain. Only running his material setup unchanged -- dielectric
everywhere, no sigma -- showed the assembly was correct all along. The controls
are the reason any of the numbers below can be believed.
"""
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MU0 = 4e-7 * math.pi
EPS0 = 8.8541878128e-12
LOC_EDGES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
FREQS = (1e0, 1e2, 1e4, 1e6, 1e8, 1e10)


# ---------------------------------------------------------------- mesh ----
def read_msh(path):
    txt = open(path).read().split("\n")
    i = txt.index("$Nodes")
    nblocks = int(txt[i + 1].split()[0])
    tag2xyz, p = {}, i + 2
    for _ in range(nblocks):
        cnt = int(txt[p].split()[3])
        p += 1
        tags = [int(t) for t in txt[p:p + cnt]]
        p += cnt
        for t, line in zip(tags, txt[p:p + cnt]):
            tag2xyz[t] = [float(v) * 1e-3 for v in line.split()]      # mm -> m
        p += cnt

    e0 = txt.index("$Entities")
    npts, ncur, nsur, nvol = (int(x) for x in txt[e0 + 1].split())
    q = e0 + 2 + npts + ncur + nsur
    ent2phys = {}
    for _ in range(nvol):
        f = txt[q].split()
        ent2phys[int(f[0])] = int(f[8]) if int(f[7]) else 0
        q += 1

    i = txt.index("$Elements")
    nblocks = int(txt[i + 1].split()[0])
    p, tets, mats = i + 2, [], []
    for _ in range(nblocks):
        b = txt[p].split()
        etag, etype, cnt = int(b[1]), int(b[2]), int(b[3])
        p += 1
        for line in txt[p:p + cnt]:
            v = [int(x) for x in line.split()]
            if etype == 4:
                tets.append(v[1:5])
                mats.append(ent2phys.get(etag, 0))
        p += cnt

    tags = sorted(tag2xyz)
    remap = {t: k for k, t in enumerate(tags)}
    return (np.array([tag2xyz[t] for t in tags]),
            np.array([[remap[v] for v in t] for t in tets]),
            np.array(mats))


def build_edges(tets):
    s = set()
    for t in tets:
        for a, b in LOC_EDGES:
            s.add((min(t[a], t[b]), max(t[a], t[b])))
    edges = sorted(s)
    return edges, {e: k for k, e in enumerate(edges)}


def tet_geometry(nodes, t):
    M = np.ones((4, 4))
    M[:, 1:] = nodes[t]
    V = abs(np.linalg.det(M)) / 6.0
    return np.linalg.inv(M)[1:, :].T, V


def boundary_edges(nodes, tets, edges):
    from collections import defaultdict
    cnt = defaultdict(int)
    for t in tets:
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            cnt[tuple(sorted((t[f[0]], t[f[1]], t[f[2]])))] += 1
    on = set()
    for f, c in cnt.items():
        if c == 1:
            for a, b in ((0, 1), (0, 2), (1, 2)):
                on.add((min(f[a], f[b]), max(f[a], f[b])))
    eset = {e: k for k, e in enumerate(edges)}
    return sorted(eset[e] for e in on if e in eset)


# ------------------------------------------------------------ assembly ----
def assemble(nodes, tets, mats, eps_of, mu_of, chi_of):
    """K1 (30), K2 (31), K_EN (32), K_NE (24) volume part, G (23).

    chi is evaluated PER ELEMENT, which §11 shows is necessary once conductors
    exist. The surface term of (24) is omitted -- see README.
    """
    edges, eidx = build_edges(tets)
    Ne, Nn = len(edges), len(nodes)
    K1 = np.zeros((Ne, Ne), dtype=complex)
    K2 = np.zeros((Ne, Ne), dtype=complex)
    KEN = np.zeros((Ne, Nn), dtype=complex)
    KNE = np.zeros((Nn, Ne), dtype=complex)
    G = np.zeros((Nn, Nn), dtype=complex)

    for t, mat in zip(tets, mats):
        g, V = tet_geometry(nodes, t)
        eps, nu, chi = eps_of(mat), 1.0 / mu_of(mat), chi_of(mat)
        gi, sg = [], []
        for a, b in LOC_EDGES:
            na, nb = t[a], t[b]
            gi.append(eidx[(min(na, nb), max(na, nb))])
            sg.append(1.0 if na < nb else -1.0)
        curls = [2.0 * np.cross(g[a], g[b]) for a, b in LOC_EDGES]

        def mass(i, j):
            return V / 10.0 if i == j else V / 20.0

        for m, (am, bm) in enumerate(LOC_EDGES):
            for n, (an, bn) in enumerate(LOC_EDGES):
                K1[gi[m], gi[n]] += sg[m] * sg[n] * nu * V * np.dot(curls[m], curls[n])
                K2[gi[m], gi[n]] += sg[m] * sg[n] * eps * (
                    mass(am, an) * np.dot(g[bm], g[bn])
                    - mass(am, bn) * np.dot(g[bm], g[an])
                    - mass(bm, an) * np.dot(g[am], g[bn])
                    + mass(bm, bn) * np.dot(g[am], g[an]))
            for n in range(4):
                v = (V / 4.0) * (np.dot(g[bm], g[n]) - np.dot(g[am], g[n]))
                KEN[gi[m], t[n]] += sg[m] * eps * v
                KNE[t[n], gi[m]] += -sg[m] * (eps / chi) * v
        for i in range(4):
            for j in range(4):
                G[t[i], t[j]] += mass(i, j)
    return edges, K1, K2, KEN, KNE, G


def cond_of(M):
    s = np.linalg.svd(M, compute_uv=False)
    return np.inf if s[-1] == 0 else s[0] / s[-1]


def sweep(title, nodes, tets, mats, eps_of, chi_of, note=""):
    print("\n=== %s ===" % title)
    if note:
        print("    %s" % note)
    print("    %-9s | %-12s | %-12s | %-12s" %
          ("freq", "cond gauged", "cond ungauged", "eps contrast"))
    print("    " + "-" * 56)
    for f in FREQS:
        w = 2.0 * math.pi * f
        ef, cf = (lambda m: eps_of(m, w)), (lambda m: chi_of(m, w))
        edges, K1, K2, KEN, KNE, G = assemble(
            nodes, tets, mats, ef, lambda m: MU0, cf)
        ke = np.array([i for i in range(len(edges))
                       if i not in set(boundary_edges(nodes, tets, edges))])
        kn = np.arange(len(nodes))
        Kuu = (K1 - w ** 2 * K2)[np.ix_(ke, ke)]
        S = Kuu - KEN[np.ix_(ke, kn)] @ np.linalg.solve(
            G[np.ix_(kn, kn)], KNE[np.ix_(kn, ke)])
        e1, e2 = abs(ef(1)), abs(ef(2))
        print("    %-9.0e | %-12.3e | %-12.3e | %-12.3e"
              % (f, cond_of(S), cond_of(Kuu), max(e1, e2) / min(e1, e2)))


def main():
    msh = os.path.join(HERE, "tiny.msh")
    if not os.path.isfile(msh):
        raise SystemExit("run:  gmsh tiny.geo -3 -o tiny.msh")
    nodes, tets, mats = read_msh(msh)
    print("mesh: %d nodes, %d tets  (%d conductor, %d dielectric)"
          % (len(nodes), len(tets), (mats == 1).sum(), (mats == 2).sum()))
    SIGMA, EPS_R = 5.8e7, 4.5

    sweep("CONTROL 1 -- Chew's materials: dielectric only, no sigma",
          nodes, tets, mats,
          lambda m, w: EPS0 * EPS_R,
          lambda m, w: MU0 * (EPS0 * EPS_R) ** 2,
          "the run that validates the assembly: expect FLAT gauged, "
          "blowing-up ungauged")

    sweep("CONTROL 2 -- mild contrast, eps_r 1 against 4.5, still no sigma",
          nodes, tets, mats,
          lambda m, w: EPS0 if m == 1 else EPS0 * EPS_R,
          lambda m, w: MU0 * (EPS0 if m == 1 else EPS0 * EPS_R) ** 2)

    sweep("OURS -- copper against dielectric, sigma folded into eps_eff",
          nodes, tets, mats,
          lambda m, w: (EPS0 - 1j * SIGMA / w) if m == 1 else EPS0 * EPS_R,
          lambda m, w: MU0 * ((EPS0 - 1j * SIGMA / w) if m == 1
                              else EPS0 * EPS_R) ** 2,
          "eps_eff = eps - j sigma/omega is NOT in Chew's paper")

    # --- chi constant vs per element, and the symmetry it costs -------------
    print("\n=== chi constant (his §II.C) against chi = mu eps^2 (his §II.A) ===")
    print("    %-9s | %-14s | %-14s | %-10s"
          % ("freq", "chi constant", "chi = mu eps^2", "K_NE vs K_EN^T"))
    print("    " + "-" * 60)
    for f in FREQS:
        w = 2.0 * math.pi * f

        def ef(m):
            return (EPS0 - 1j * SIGMA / w) if m == 1 else EPS0 * EPS_R

        out = {}
        for key, cf in (("const", lambda m: MU0 * EPS0 ** 2),
                        ("var", lambda m: MU0 * ef(m) ** 2)):
            edges, K1, K2, KEN, KNE, G = assemble(
                nodes, tets, mats, ef, lambda m: MU0, cf)
            ke = np.array([i for i in range(len(edges))
                           if i not in set(boundary_edges(nodes, tets, edges))])
            kn = np.arange(len(nodes))
            Kuu = (K1 - w ** 2 * K2)[np.ix_(ke, ke)]
            A, B = KEN[np.ix_(ke, kn)], KNE[np.ix_(kn, ke)]
            out[key] = cond_of(Kuu - A @ np.linalg.solve(G[np.ix_(kn, kn)], B))
            if key == "var":
                sc = np.linalg.norm(A.T) / max(np.linalg.norm(B), 1e-300)
                out["sym"] = np.abs(B * sc - A.T).max() / max(np.abs(A).max(), 1e-300)
        print("    %-9.0e | %-14.3e | %-14.3e | %-10.2e"
              % (f, out["const"], out["var"], out["sym"]))
    print("\n    K_NE vs K_EN^T: ~1e-16 means §7's symmetrization still applies;")
    print("    O(1) means a per-element chi has broken it.")


if __name__ == "__main__":
    main()
