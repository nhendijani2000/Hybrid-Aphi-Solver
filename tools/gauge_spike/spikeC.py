"""Spike C: does a materially-weighted gauge protect an INTERIOR TERMINAL?

    gmsh stub.geo -3 -o stub.msh
    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvpython.exe" spikeC.py

THE QUESTION. 08_MixedPort_Interior measures our solver losing a
mixed-material interior port: terminal L moves 3.1 % and J 77 % with the choice
of spanning tree, against 1e-12 and 1e-11 for the same wire with both caps on
the boundary. GAUGE_CHOICE.md 15.9 then derives that NO tree can fix it -- every
route to pinning psi on an interior face either leaves it path dependent, hence
unprotected, or imposes a spurious flux constraint, hence wrong.

So the gauge has to change. 13 found that Jochum, Farle & Dyczij-Edlinger carry
an explicit GAUGE INTERFACE CONDITION at the conductor/insulator boundary, their
(32), weighted by the materials of both sides -- and that their 5.5 puts the
gauge unknown psi ON that interface, which is exactly where our port sits.
compare.py then measured their system flat across ten decades, complex-symmetric
to 1e-16, and insensitive to the sigma/omega-eps contrast that wrecks Chew's.

None of that involved a PORT. This does.

MILESTONE 1, WHICH IS WHAT THIS FILE DOES FIRST AND ONLY.
Reproduce the FAILURE in miniature, with plain tree-cotree. If a few-hundred-tet
spike cannot show an interior terminal's potential moving with the tree, then it
cannot show a gauge fixing it either, and every later number here would be
unfalsifiable. Two of this project's five spike runs have already produced
confident wrong conclusions for want of exactly this discipline
(RESULTS_COMPARE.md, "two errors in this run").

Only once that is in hand does Jochum's system get a terminal bolted onto it.

CONVENTION. e^{+jwt}, so E = -jwA - grad(Phi) and B = curl A.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from spike import (read_msh, build_edges, tet_geometry, boundary_edges,      # noqa
                   LOC_EDGES, MU0, EPS0)

CONDUCTOR, DIELECTRIC = 1, 2
SIGMA, EPS_R = 5.8e7, 4.5
FREQ = 50.0
I_DRIVE = 1.0


# ------------------------------------------------------------- terminal ----
def interior_terminal(nodes, tets, mats):
    """The conductor's +z face: interior, conductor below, dielectric above.

    This is 08_MixedPort_Interior's port in miniature. Found from the mesh
    rather than from a tag, so embedded.geo needs no edit: take every face
    whose two tets straddle the material interface, keep those whose outward
    normal (conductor -> dielectric) points +z, and that is the top face of the
    embedded cube.
    """
    from collections import defaultdict
    face_tets = defaultdict(list)
    for t, tet in enumerate(tets):
        for f in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
            key = tuple(sorted((tet[f[0]], tet[f[1]], tet[f[2]])))
            face_tets[key].append(t)

    faces, zs = [], []
    for key, ts in face_tets.items():
        if len(ts) != 2:
            continue
        m = [mats[t] for t in ts]
        if sorted(m) != [CONDUCTOR, DIELECTRIC]:
            continue
        # Conductor-side tet's apex tells us which way is "into the conductor".
        c_tet = ts[0] if m[0] == CONDUCTOR else ts[1]
        apex = [v for v in tets[c_tet] if v not in key][0]
        centroid = nodes[list(key)].mean(axis=0)
        if nodes[apex][2] < centroid[2] - 1e-12:     # conductor below => +z face
            faces.append(key)
            zs.append(centroid[2])
    if not faces:
        raise SystemExit("found no conductor/dielectric +z face")
    # Keep the single planar sheet at the highest z (the cube's top).
    zmax = max(zs)
    sheet = [f for f, z in zip(faces, zs) if z > zmax - 1e-9]
    tnodes = sorted({v for f in sheet for v in f})
    return sheet, tnodes


# ------------------------------------------------------------------ tree ----
def spanning_tree(nodes, edges, bedges, seed):
    """Boundary-first spanning tree, with the outer boundary contracted.

    `seed` reorders the edge scan, which is what makes two different trees over
    the same geometry -- the same lever permute_nodes.py pulls on the real
    solver, where reordering the .msh changes internal edge numbering and hence
    the traversal.
    """
    on_boundary = set()
    for k in bedges:
        on_boundary.update(edges[k])

    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    SUPER = -1
    parent[SUPER] = SUPER
    for n in range(len(nodes)):
        parent[n] = SUPER if n in on_boundary else n

    order = list(range(len(edges)))
    if seed is not None:
        np.random.RandomState(seed).shuffle(order)

    tree = []
    for k in order:
        u, v = edges[k]
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
            tree.append(k)
    return sorted(tree)


# -------------------------------------------------------------- assembly ----
def assemble_aphi(nodes, tets, mats, omega):
    """OUR formulation: curl(nu curl A) + alpha A + beta grad(Phi) = J, and its
    divergence. Edge A, P1 nodal Phi (the real solver uses P2; P1 is enough to
    show a gauge effect and keeps the dense algebra small).

        alpha = j w sigma - w^2 eps      beta = sigma + j w eps
    """
    edges, eidx = build_edges(tets)
    Ne, Nn = len(edges), len(nodes)
    K = np.zeros((Ne, Ne), dtype=complex)      # curl-curl + alpha mass
    C = np.zeros((Ne, Nn), dtype=complex)      # beta * edge x grad(node)
    G = np.zeros((Nn, Nn), dtype=complex)      # beta * grad.grad

    for t, mat in zip(tets, mats):
        g, V = tet_geometry(nodes, t)
        sig = SIGMA if mat == CONDUCTOR else 0.0
        eps = EPS0 if mat == CONDUCTOR else EPS0 * EPS_R
        nu = 1.0 / MU0
        alpha = 1j * omega * sig - omega ** 2 * eps
        beta = sig + 1j * omega * eps
        gi, sg = [], []
        for a, b in LOC_EDGES:
            na, nb = t[a], t[b]
            gi.append(eidx[(min(na, nb), max(na, nb))])
            sg.append(1.0 if na < nb else -1.0)
        curls = [2.0 * np.cross(g[a], g[b]) for a, b in LOC_EDGES]
        mass = lambda i, j: V / 10.0 if i == j else V / 20.0

        for m, (am, bm) in enumerate(LOC_EDGES):
            for n, (an, bn) in enumerate(LOC_EDGES):
                s = sg[m] * sg[n]
                K[gi[m], gi[n]] += s * (
                    nu * V * np.dot(curls[m], curls[n])
                    + alpha * (mass(am, an) * np.dot(g[bm], g[bn])
                               - mass(am, bn) * np.dot(g[bm], g[an])
                               - mass(bm, an) * np.dot(g[am], g[bn])
                               + mass(bm, bn) * np.dot(g[am], g[an])))
            for n in range(4):
                C[gi[m], t[n]] += sg[m] * beta * (V / 4.0) * (
                    np.dot(g[bm], g[n]) - np.dot(g[am], g[n]))
        for i in range(4):
            for j in range(4):
                G[t[i], t[j]] += beta * V * np.dot(g[i], g[j])
    return edges, K, C, G


# ---------------------------------------------------------------- solve ----
def solve_tree_cotree(nodes, tets, mats, omega, seed):
    """One solve under the tree-cotree gauge. Returns the terminal's Phi.

    The terminal is ONE aggregated Phi unknown over its face -- Stysch 6.6.2's
    A = P A_hat P^T, which is also how the real solver's PhiDof::Port behaves
    (every node of the surface reads one unknown). Phi = 0 on the outer
    boundary is the reference; 1 A is injected at the terminal.
    """
    edges, K, C, G = assemble_aphi(nodes, tets, mats, omega)
    Ne, Nn = len(edges), len(nodes)
    bedges = boundary_edges(nodes, tets, edges)
    bnodes = sorted({v for k in bedges for v in edges[k]})
    _, tnodes = interior_terminal(nodes, tets, mats)

    tree = set(spanning_tree(nodes, edges, bedges, seed))
    ke = np.array([k for k in range(Ne)
                   if k not in tree and k not in set(bedges)])

    # Phi unknowns: free interior nodes, then ONE for the terminal.
    bset, tset = set(bnodes), set(tnodes)
    free_n = [i for i in range(Nn) if i not in bset and i not in tset]
    col = {n: j for j, n in enumerate(free_n)}
    nf, term_col = len(free_n), len(free_n)

    P = np.zeros((Nn, nf + 1))
    for n in free_n:
        P[n, col[n]] = 1.0
    for n in tnodes:
        P[n, term_col] = 1.0          # the aggregation: one shared unknown

    Kuu = K[np.ix_(ke, ke)]
    Cu = C[np.ix_(ke, np.arange(Nn))] @ P
    Gp = P.T @ G @ P

    # THE SECOND ROW IS CONTINUITY, AND IT CARRIES A jw THAT IS EASY TO DROP.
    #
    #     div[ beta (jw A + grad Phi) ] = 0   =>   jw C^T a + G phi = 0
    #
    # The first version of this file assembled C^T a + G phi and got a terminal
    # potential of -4.3e8 - 1.5e9j V for 1 A into a copper stub whose
    # resistance is about 6 microohms -- fourteen orders out. The magnitude
    # check is what caught it; the tree-dependence it "reproduced" was
    # meaningless.
    #
    # Dividing the row by jw restores symmetry, since then (2,1) = C^T against
    # (1,2) = C, and rescales the drive to I/(jw). That is the same device
    # Zhao & Fu apply as -j/omega (GAUGE_CHOICE.md Sec. 12.1) and Balian et al.
    # as a congruence (Sec. 11.1): scale the equation, scale the unknown's
    # right-hand side to match, and symmetry survives.
    jw = 1j * omega
    M = np.block([[Kuu, Cu], [Cu.T, Gp / jw]])
    b = np.zeros(M.shape[0], dtype=complex)
    b[len(ke) + term_col] = I_DRIVE / jw

    asym = np.abs(M - M.T).max() / np.abs(M).max()
    x = np.linalg.solve(M, b)
    return x[len(ke) + term_col], len(ke), nf + 1, asym


C0 = 1.0 / math.sqrt(MU0 * EPS0)
ETA0 = math.sqrt(MU0 / EPS0)


def assemble_jochum(nodes, tets, mats, k0):
    """Jochum (18)-(19) in his non-dimensional variables, as compare.py.

        kappa = sigma_0 + j k0 eps_r,  sigma_0 = eta0 sigma,  nu_r = 1

    Non-dimensionalizing is not cosmetic: compare.py's first run transcribed
    his structure into SI and reported the condition number RISING with
    frequency, because his A is a scaled potential and that scaling is part of
    the stabilization (RESULTS_COMPARE.md).
    """
    edges, eidx = build_edges(tets)
    Ne, Nn = len(edges), len(nodes)
    z = lambda r, c: np.zeros((r, c), dtype=complex)
    S, M_kap = z(Ne, Ne), z(Ne, Ne)
    C_kap, G_kap, G_epsN = z(Ne, Nn), z(Nn, Nn), z(Nn, Nn)

    for t, mat in zip(tets, mats):
        g, V = tet_geometry(nodes, t)
        cond = mat == CONDUCTOR
        er = 1.0 if cond else EPS_R
        kap = (ETA0 * SIGMA + 1j * k0 * er) if cond else 1j * k0 * er
        gi, sg = [], []
        for a, b in LOC_EDGES:
            na, nb = t[a], t[b]
            gi.append(eidx[(min(na, nb), max(na, nb))])
            sg.append(1.0 if na < nb else -1.0)
        curls = [2.0 * np.cross(g[a], g[b]) for a, b in LOC_EDGES]
        mass = lambda i, j: V / 10.0 if i == j else V / 20.0

        for m, (am, bm) in enumerate(LOC_EDGES):
            for n, (an, bn) in enumerate(LOC_EDGES):
                s = sg[m] * sg[n]
                S[gi[m], gi[n]] += s * V * np.dot(curls[m], curls[n])
                M_kap[gi[m], gi[n]] += s * kap * (
                    mass(am, an) * np.dot(g[bm], g[bn])
                    - mass(am, bn) * np.dot(g[bm], g[an])
                    - mass(bm, an) * np.dot(g[am], g[bn])
                    + mass(bm, bn) * np.dot(g[am], g[an]))
            for n in range(4):
                C_kap[gi[m], t[n]] += sg[m] * kap * (V / 4.0) * (
                    np.dot(g[bm], g[n]) - np.dot(g[am], g[n]))
        for i in range(4):
            for j in range(4):
                gg = V * np.dot(g[i], g[j])
                G_kap[t[i], t[j]] += kap * gg
                if not cond:
                    G_epsN[t[i], t[j]] += er * gg
    return edges, S, M_kap, C_kap, G_kap, G_epsN


def solve_jochum(nodes, tets, mats, omega, seed):
    """Jochum's region-split gauge, with the SAME interior terminal.

    THE HYPOTHESIS THIS TESTS. 15.9 derives that no tree can protect an
    interior terminal, because psi there is fixed by a tree path and so is
    path dependent. Jochum does not let the tree decide: psi is a genuine
    UNKNOWN, carried on the conductor/insulator interface by his Sec. 5.5 --
    and our terminal lies ON that interface. If psi at the terminal is
    determined by the equations rather than by the traversal, the terminal
    potential should stop moving.

    The tree still appears, because A = A_c + grad(psi) needs A_c in the
    cotree space (his (14), realized by tree-cotree). So the tree is still
    permuted; what changes is whether the ANSWER depends on it.
    """
    k0 = omega / C0
    edges, S, M_kap, C_kap, G_kap, G_epsN = assemble_jochum(
        nodes, tets, mats, k0)
    Ne, Nn = len(edges), len(nodes)
    bedges = boundary_edges(nodes, tets, edges)
    bnodes = set(v for k in bedges for v in edges[k])
    _, tnodes = interior_terminal(nodes, tets, mats)

    tree = set(spanning_tree(nodes, edges, bedges, seed))
    ke = np.array([k for k in range(Ne)
                   if k not in tree and k not in set(bedges)])

    # psi per Jochum 5.5: on the interface and inside the conductor, dropped
    # strictly inside the insulator -- compare.py measured that restriction at
    # nine orders, and the block reading says an insulator-interior psi row
    # vanishes as k0 -> 0.
    node_mats = [set() for _ in range(Nn)]
    for t, m in zip(tets, mats):
        for n in t:
            node_mats[n].add(m)
    psi_n = np.array([n for n in range(Nn) if n not in bnodes
                      and (len(node_mats[n]) > 1
                           or node_mats[n] == {CONDUCTOR})])

    # V: free interior nodes plus ONE aggregated terminal unknown.
    tset = set(tnodes)
    free_n = [i for i in range(Nn) if i not in bnodes and i not in tset]
    nf = len(free_n)
    P = np.zeros((Nn, nf + 1))
    for j, n in enumerate(free_n):
        P[n, j] = 1.0
    for n in tnodes:
        P[n, nf] = 1.0

    jw = 1j * k0
    ix = np.ix_
    A11 = (S + jw * M_kap)[ix(ke, ke)]
    A12 = jw * C_kap[ix(ke, psi_n)]
    A13 = C_kap[ix(ke, np.arange(Nn))] @ P
    A22 = jw * G_kap[ix(psi_n, psi_n)]
    A23 = G_kap[ix(psi_n, np.arange(Nn))] @ P
    A33 = P.T @ G_epsN @ P

    M = np.block([[A11, A12, A13],
                  [A12.T, A22, A23],
                  [A13.T, A23.T, A33]])
    b = np.zeros(M.shape[0], dtype=complex)
    # The terminal's row is its current balance. eta0 carries the same
    # non-dimensionalization the operator got.
    b[len(ke) + len(psi_n) + nf] = I_DRIVE * ETA0 / jw

    asym = np.abs(M - M.T).max() / np.abs(M).max()
    x = np.linalg.solve(M, b)
    # Back to volts: Phi_tilde = Phi / (j mu0 c0) is not used here -- the
    # terminal row was scaled by eta0/jw, so x is already in volts.
    return x[len(ke) + len(psi_n) + nf], len(ke), len(psi_n), nf + 1, asym


def main():
    msh = os.path.join(HERE, "stub.msh")
    if not os.path.isfile(msh):
        raise SystemExit("run:  gmsh stub.geo -3 -o stub.msh")
    nodes, tets, mats = read_msh(msh)
    omega = 2.0 * math.pi * FREQ

    sheet, tnodes = interior_terminal(nodes, tets, mats)
    print("MILESTONE 1 -- can a spike reproduce the interior-terminal failure?\n")
    print("mesh      %d nodes, %d tets (%d conductor, %d dielectric)"
          % (len(nodes), len(tets), (mats == CONDUCTOR).sum(),
             (mats == DIELECTRIC).sum()))
    print("terminal  the conductor's +z face: %d faces, %d nodes, one shared "
          "Phi unknown" % (len(sheet), len(tnodes)))
    print("          conductor below, dielectric above -- interior and mixed\n")

    print("    %-10s %-34s %-11s %s"
          % ("tree", "terminal Phi (V)", "vs base", "asym"))
    print("    " + "-" * 70)
    base = None
    for label, seed in (("base", None), ("permA", 12345), ("permB", 777)):
        phi, nA, nP, asym = solve_tree_cotree(nodes, tets, mats, omega, seed)
        if base is None:
            base = phi
            rel = 0.0
        else:
            rel = abs(phi - base) / abs(base)
        print("    %-10s %+.9e %+.9ej  %-11.3e %.1e"
              % (label, phi.real, phi.imag, rel, asym))

    print("\n    unknowns: %d cotree A, %d Phi (last one is the terminal)"
          % (nA, nP))
    r_exact = 0.006 / (SIGMA * 0.003 * 0.003)
    print("    R = %.4f uOhm against the exact L/(sigma A) = %.4f uOhm  "
          "(%.2f %%)"
          % (base.real * 1e6, r_exact * 1e6,
             100.0 * abs(base.real - r_exact) / r_exact))

    print("\n    A large 'vs base' means the spike reproduces what")
    print("    08_MixedPort_Interior measures, and Jochum can then be")
    print("    bolted onto the same setup and compared. A tiny one means")
    print("    this spike cannot see the effect and must be fixed first.")

    milestone2(nodes, tets, mats, omega, r_exact)
    milestone3(nodes, tets, mats, omega, r_exact)
    milestone4(nodes, tets, mats, omega, r_exact)
    milestone5(nodes, tets, mats, r_exact)




def milestone2(nodes, tets, mats, omega, r_exact):
    """Jochum on the same stub, same terminal, same permuted trees."""
    print("\n" + "=" * 72)
    print("MILESTONE 2 -- Jochum's region-split gauge, same terminal")
    print("=" * 72)
    print("    %-10s %-34s %-11s %s"
          % ("tree", "terminal Phi (V)", "vs base", "asym"))
    print("    " + "-" * 70)
    base, rows = None, []
    for label, seed in (("base", None), ("permA", 12345), ("permB", 777)):
        phi, nA, nps, nV, asym = solve_jochum(nodes, tets, mats, omega, seed)
        rel = 0.0 if base is None else abs(phi - base) / abs(base)
        if base is None:
            base = phi
        rows.append((label, phi, rel))
        print("    %-10s %+.9e %+.9ej  %-11.3e %.1e"
              % (label, phi.real, phi.imag, rel, asym))
    print("\n    unknowns: %d cotree A, %d psi (interface + conductor), "
          "%d V" % (nA, nps, nV))
    print("\n    MAGNITUDE CHECK, which decides whether the invariance above")
    print("    means anything: R must come out near the exact %.3f uOhm."
          % (r_exact * 1e6))
    print("    Jochum R = %.4f uOhm  (%.2f %% from exact)"
          % (base.real * 1e6, 100.0 * abs(base.real - r_exact) / r_exact))
    return rows




# ----------------------------------------------------------------- Chew ----
def solve_chew(nodes, tets, mats, omega, seed, alpha=1.0):
    """Chew's DECOUPLED Phi equation, with the same interior terminal.

        div(eps grad Phi) + chi w^2 Phi = -rho,     chi = alpha mu eps^2
                                                   (GENERALIZED_LORENZ_GAUGE (2),(4))

    WHY THIS IS STRUCTURALLY DIFFERENT FROM JOCHUM. Equation (2) contains NO A.
    Phi is determined by its own Helmholtz problem, so a terminal is an ordinary
    boundary condition ON THAT EQUATION rather than a perturbation of a gauge
    row -- which is what milestone 2 found Jochum has no answer for. And since
    no tree enters (2) at all, Phi cannot depend on the tree: the invariance is
    exact by construction, not by good conditioning.

    `seed` is accepted and IGNORED, deliberately. If (2) is solved as written,
    passing a different tree must change nothing, and showing that is the point.

    SIGMA. Chew's paper has no finite conductivity -- his conductors are PEC
    boundary conditions and his eps varies by 4.5 (RESULTS_COMPARE.md). Ours is
    copper, so eps -> eps_eff = eps - j sigma/omega, which is section 1's
    addition and which spike.py already measured wrecking his conditioning. The
    question HERE is not conditioning but whether the resulting Phi is right.

    Weak form, testing with lam and integrating the divergence by parts:

        -int eps_eff grad(Phi).grad(lam)  +  w^2 int chi Phi lam  =  -int rho lam

    and the boundary term int eps_eff dPhi/dn lam dS is where the terminal
    current enters: with eps_eff ~ -j sigma/w, eps_eff dPhi/dn ~ (j/w) J.n, so
    injecting I puts (j/w) I on the terminal's row.
    """
    Nn = len(nodes)
    K = np.zeros((Nn, Nn), dtype=complex)      # the (2) operator
    for t, mat in zip(tets, mats):
        g, V = tet_geometry(nodes, t)
        if mat == CONDUCTOR:
            eps = EPS0 - 1j * SIGMA / omega
        else:
            eps = EPS0 * EPS_R
        chi = alpha * MU0 * eps ** 2
        mass = lambda i, j: V / 10.0 if i == j else V / 20.0
        for i in range(4):
            for j in range(4):
                K[t[i], t[j]] += (-eps * V * np.dot(g[i], g[j])
                                  + omega ** 2 * chi * mass(i, j))

    edges, _ = build_edges(tets)
    bedges = boundary_edges(nodes, tets, edges)
    bnodes = set(v for k in bedges for v in edges[k])
    _, tnodes = interior_terminal(nodes, tets, mats)

    tset = set(tnodes)
    free_n = [i for i in range(Nn) if i not in bnodes and i not in tset]
    nf = len(free_n)
    P = np.zeros((Nn, nf + 1))
    for j, n in enumerate(free_n):
        P[n, j] = 1.0
    for n in tnodes:
        P[n, nf] = 1.0

    M = P.T @ K @ P
    b = np.zeros(nf + 1, dtype=complex)
    b[nf] = 1j * I_DRIVE / omega

    asym = np.abs(M - M.T).max() / np.abs(M).max()
    x = np.linalg.solve(M, b)
    return x[nf], nf + 1, asym


def milestone3(nodes, tets, mats, omega, r_exact):
    print("\n" + "=" * 72)
    print("MILESTONE 3 -- Chew's decoupled Phi equation, same terminal")
    print("=" * 72)
    print("    %-10s %-34s %-11s %s"
          % ("tree", "terminal Phi (V)", "vs base", "asym"))
    print("    " + "-" * 70)
    base = None
    for label, seed in (("base", None), ("permA", 12345), ("permB", 777)):
        phi, nP, asym = solve_chew(nodes, tets, mats, omega, seed)
        rel = 0.0 if base is None else abs(phi - base) / abs(base)
        if base is None:
            base = phi
        print("    %-10s %+.9e %+.9ej  %-11.3e %.1e"
              % (label, phi.real, phi.imag, rel, asym))
    print("\n    %d Phi unknowns, no A and no tree in this equation at all,"
          % nP)
    print("    so the tree-invariance above is exact by construction.")
    print("\n    MAGNITUDE CHECK -- the only thing that decides whether it means")
    print("    anything.  exact R = %.4f uOhm" % (r_exact * 1e6))
    print("    Chew |Phi|/I = %.6e Ohm,  Re = %.4e,  Im = %.4e"
          % (abs(base), base.real, base.imag))
    print("    ratio to exact R: %.4e" % (abs(base) / r_exact))

    # alpha is Chew's free parameter; he recommends 1. Does the answer depend
    # on it? It must not, if Phi is a physical potential.
    print("\n    dependence on Chew's free parameter alpha (he recommends 1):")
    print("    %-12s %-20s %s" % ("alpha", "|Phi| (V)", "vs alpha=1"))
    ref = None
    for al in (1.0, 1e-3, 1e3):
        phi, _, _ = solve_chew(nodes, tets, mats, omega, None, alpha=al)
        if ref is None:
            ref = abs(phi)
        print("    %-12.0e %-20.9e %.3e"
              % (al, abs(phi), abs(abs(phi) - ref) / ref))




# ------------------------------------------- Phi as a post-process of E ----
def solve_phi_from_E(nodes, tets, mats, omega, seed):
    """Stysch 3.1 / 6.3: solve the FIELD first, then Phi from its own BVP.

        int eps grad(Phi).grad(lam)  =  - int eps E . grad(lam)

    i.e. Phi's source is div(eps E), their (6.25) with the g term folded into E.

    WHY THIS SHOULD BE TREE-INDEPENDENT, AND IS THE POINT OF THE TEST.
    E is gauge invariant -- 08_MixedPort_Interior's control measures it at
    1e-11 for a boundary port, and theory says so for any port. If Phi is
    determined by its OWN boundary value problem whose only input is E, then
    Phi inherits E's invariance no matter how the tree was chosen. Phi stops
    being a tree artefact and becomes a derived field.

    This is the third route to a unique Phi in GAUGE_CHOICE.md: Chew via chi,
    Ostrowski & Hiptmair via the EQS gauge, Stysch via the Lorenz PDE. It is
    also what milestone 3's alpha -> 0 limit degenerated into, except that
    milestone 3 had no source and therefore no reactance.

    Step 1 is still the ordinary tree-cotree solve, so the tree is still
    permuted; what is tested is whether the ANSWER survives it.
    """
    # --- step 1: the tree-cotree solve, same as milestone 1 -----------------
    edges, K, C, G = assemble_aphi(nodes, tets, mats, omega)
    Ne, Nn = len(edges), len(nodes)
    bedges = boundary_edges(nodes, tets, edges)
    bnodes = set(v for k in bedges for v in edges[k])
    _, tnodes = interior_terminal(nodes, tets, mats)
    tree = set(spanning_tree(nodes, edges, bedges, seed))
    ke = np.array([k for k in range(Ne)
                   if k not in tree and k not in set(bedges)])

    tset = set(tnodes)
    free_n = [i for i in range(Nn) if i not in bnodes and i not in tset]
    nf = len(free_n)
    P = np.zeros((Nn, nf + 1))
    for j, n in enumerate(free_n):
        P[n, j] = 1.0
    for n in tnodes:
        P[n, nf] = 1.0

    jw = 1j * omega
    Kuu = K[np.ix_(ke, ke)]
    Cu = C[np.ix_(ke, np.arange(Nn))] @ P
    M1 = np.block([[Kuu, Cu], [Cu.T, (P.T @ G @ P) / jw]])
    b1 = np.zeros(M1.shape[0], dtype=complex)
    b1[len(ke) + nf] = I_DRIVE / jw
    x1 = np.linalg.solve(M1, b1)

    a_full = np.zeros(Ne, dtype=complex)
    a_full[ke] = x1[:len(ke)]
    phi1 = P @ x1[len(ke):]

    # --- step 2: E per tet, which is gauge invariant ------------------------
    eidx = {e: k for k, e in enumerate(edges)}
    Kp = np.zeros((Nn, Nn), dtype=complex)
    rhs = np.zeros(Nn, dtype=complex)
    for t, mat in zip(tets, mats):
        g, V = tet_geometry(nodes, t)
        eps = (EPS0 - 1j * SIGMA / omega) if mat == CONDUCTOR else EPS0 * EPS_R
        # A and grad(Phi) at the tet centroid, so E = -jw A - grad(Phi).
        Avec = np.zeros(3, dtype=complex)
        for m, (am, bm) in enumerate(LOC_EDGES):
            na, nb = t[am], t[bm]
            k = eidx[(min(na, nb), max(na, nb))]
            s = 1.0 if na < nb else -1.0
            # Whitney-1 at the centroid: lam = 1/4 each.
            Avec += s * a_full[k] * 0.25 * (g[bm] - g[am])
        gradphi = sum(phi1[t[i]] * g[i] for i in range(4))
        E = -jw * Avec - gradphi
        for i in range(4):
            for j in range(4):
                Kp[t[i], t[j]] += eps * V * np.dot(g[i], g[j])
            rhs[t[i]] += -eps * V * np.dot(E, g[i])

    Mp = P.T @ Kp @ P
    bp = P.T @ rhs
    xp = np.linalg.solve(Mp, bp)
    return xp[nf], x1[len(ke) + nf]


def milestone4(nodes, tets, mats, omega, r_exact):
    print("\n" + "=" * 72)
    print("MILESTONE 4 -- Phi as its OWN BVP driven by E (Stysch 3.1)")
    print("=" * 72)
    print("    %-10s %-30s %-11s %-30s %s"
          % ("tree", "Phi from its own BVP", "vs base", "Phi from the A-Phi solve",
             "vs base"))
    print("    " + "-" * 104)
    b2 = b1 = None
    for label, seed in (("base", None), ("permA", 12345), ("permB", 777)):
        p2, p1 = solve_phi_from_E(nodes, tets, mats, omega, seed)
        r2 = 0.0 if b2 is None else abs(p2 - b2) / abs(b2)
        r1 = 0.0 if b1 is None else abs(p1 - b1) / abs(b1)
        if b2 is None:
            b2, b1 = p2, p1
        print("    %-10s %+.6e%+.6ej  %-11.3e %+.6e%+.6ej  %.3e"
              % (label, p2.real, p2.imag, r2, p1.real, p1.imag, r1))
    print("\n    left  = Phi solved from div(eps grad Phi) = div(eps E)")
    print("    right = Phi read straight out of the gauged A-Phi system")
    print("            (milestone 1, which moves 3.8 %)")
    print("\n    R = %.6f uOhm against exact %.6f  (%.3f %%)"
          % (b2.real * 1e6, r_exact * 1e6,
             100.0 * abs(b2.real - r_exact) / r_exact))
    print("    X = %.6e Ohm   ->   L = %.4f nH"
          % (b2.imag, b2.imag / omega * 1e9))




# ------------------------- Jochum WITH a terminal, per JOCHUM_PORT_DERIVATION 6 ----
def solve_jochum_port(nodes, tets, mats, omega, seed):
    """Jochum with psi aggregated on the terminal too, and the drive on PSI's row.

    TESTING docs/JOCHUM_PORT_DERIVATION.md Sec. 6.

    Milestone 2 put the current on V's row and got a terminal potential that was
    identical across three trees, had zero reactance, and did not move from
    50 Hz to 50 GHz. Sec. 5 of that report explains it: his Sec. 5.2 tests
    Ampere by w_c AND by grad(psi-bar), and the gauge by V-bar, so

        psi's row is CONTINUITY  and  V's row is the GAUGE

    which inverts the usual A-Phi correspondence. The V_T row carries the gauge
    flux n.(kappa A); the current is n.(kappa E) and appears in psi's row. The
    drive was perturbing the gauge instead of injecting current.

    Sec. 6's proposal: a terminal needs TWO conditions, so give it TWO unknowns.
    Aggregate psi on Gamma_T as well as V --

        V_T   one unknown,  row = aggregated gauge
        psi_T one unknown,  row = aggregated continuity = I

    -- which balances the count (2 unknowns, 2 rows) and keeps symmetry, since
    both aggregations are congruences P(.)P^T. Aggregating psi is a GAUGE
    choice, not a physical constraint: it is what contracting a surface in the
    tree does, and what n x A = 0 achieves on the outer boundary.

    THE DRIVE'S SCALE. In the non-dimensional system the Ampere equation is
    nu_r curl curl A~ + kappa~(j k0 A~ + grad Phi) = eta0 J. Testing by
    grad(psi-bar) and integrating by parts puts eta0 * I on the terminal's row.
    If that factor is wrong the recovered R is wrong by exactly it, which the
    magnitude check below will show rather than hide.
    """
    k0 = omega / C0
    edges, S, M_kap, C_kap, G_kap, G_epsN = assemble_jochum(
        nodes, tets, mats, k0)
    Ne, Nn = len(edges), len(nodes)
    bedges = boundary_edges(nodes, tets, edges)
    bnodes = set(v for k in bedges for v in edges[k])
    _, tnodes = interior_terminal(nodes, tets, mats)
    tset = set(tnodes)

    tree = set(spanning_tree(nodes, edges, bedges, seed))
    ke = np.array([k for k in range(Ne)
                   if k not in tree and k not in set(bedges)])

    node_mats = [set() for _ in range(Nn)]
    for t, m in zip(tets, mats):
        for n in t:
            node_mats[n].add(m)

    # psi per Sec. 5.5: interface and conductor, never strictly inside the
    # insulator. The terminal IS on the interface, so psi_T survives the DC
    # limit -- the expectation Sec. 6.1 flags as untested.
    psi_all = [n for n in range(Nn) if n not in bnodes
               and (len(node_mats[n]) > 1 or node_mats[n] == {CONDUCTOR})]
    psi_free = [n for n in psi_all if n not in tset]
    P_psi = np.zeros((Nn, len(psi_free) + 1))
    for j, n in enumerate(psi_free):
        P_psi[n, j] = 1.0
    for n in tnodes:
        P_psi[n, len(psi_free)] = 1.0          # psi_T
    psi_T = len(psi_free)

    v_free = [n for n in range(Nn) if n not in bnodes and n not in tset]
    P_v = np.zeros((Nn, len(v_free) + 1))
    for j, n in enumerate(v_free):
        P_v[n, j] = 1.0
    for n in tnodes:
        P_v[n, len(v_free)] = 1.0              # V_T
    v_T = len(v_free)

    jw = 1j * k0
    ix = np.ix_
    A11 = (S + jw * M_kap)[ix(ke, ke)]
    A12 = jw * (C_kap[ix(ke, np.arange(Nn))] @ P_psi)
    A13 = C_kap[ix(ke, np.arange(Nn))] @ P_v
    A22 = jw * (P_psi.T @ G_kap @ P_psi)
    A23 = P_psi.T @ G_kap @ P_v
    A33 = P_v.T @ G_epsN @ P_v

    M = np.block([[A11, A12, A13],
                  [A12.T, A22, A23],
                  [A13.T, A23.T, A33]])
    b = np.zeros(M.shape[0], dtype=complex)
    b[len(ke) + psi_T] = I_DRIVE * ETA0        # the current, on PSI's row

    asym = np.abs(M - M.T).max() / np.abs(M).max()
    x = np.linalg.solve(M, b)
    V_T = x[len(ke) + len(psi_free) + 1 + v_T]
    return V_T, len(ke), len(psi_free) + 1, len(v_free) + 1, asym


def milestone5(nodes, tets, mats, r_exact):
    print("\n" + "=" * 72)
    print("MILESTONE 5 -- JOCHUM_PORT_DERIVATION Sec. 6: psi aggregated too,")
    print("               drive moved from V's row to psi's")
    print("=" * 72)
    print("    %-9s %-11s %-30s %-11s %s"
          % ("freq", "k0", "terminal V (V)", "vs base", "asym"))
    print("    " + "-" * 76)
    for f in (5e1, 5e4, 5e7):
        base = None
        for lab, sd in (("base", None), ("permA", 12345), ("permB", 777)):
            w = 2.0 * math.pi * f
            V, nA, nps, nV, asym = solve_jochum_port(nodes, tets, mats, w, sd)
            rel = 0.0 if base is None else abs(V - base) / abs(base)
            if base is None:
                base = V
            tag = "%.0e %s" % (f, lab) if lab != "base" else "%.0e base" % f
            print("    %-21s %+.6e%+.6ej  %-11.3e %.1e"
                  % (tag, V.real, V.imag, rel, asym))
        print()
    print("    ACCEPTANCE (RESULTS_SPIKEC.md): R must be %.6f uOhm, the"
          % (r_exact * 1e6))
    print("    reactance must be present and FREQUENCY DEPENDENT, and the")
    print("    terminal potential must not move with the tree.")
    print("\n    Milestone 2's failure signature -- zero reactance, unchanged")
    print("    across decades -- means the drive is still in the gauge row and")
    print("    Sec. 6 is wrong.")


if __name__ == "__main__":
    main()
