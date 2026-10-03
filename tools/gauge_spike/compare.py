"""Chew's generalized Lorenz against Jochum's region-split gauge, same mesh.

    gmsh tiny.geo -3 -o tiny.msh
    "C:\\Program Files\\ParaView 6.1.1\\bin\\pvpython.exe" compare.py

WHY THIS EXISTS. `spike.py` measured Chew's formulation alone and reported three
problems: finite conductivity wrecks the conditioning, a position-dependent chi
breaks the symmetry of the section 7 block form, and the surface term of (24)
was omitted with no account of what it does.

docs/GAUGE_CHOICE.md section 13 then found Jochum, Farle & Dyczij-Edlinger
(SCEE 2014, pp. 63-71), which addresses all three AT ONCE, in our own
discretization -- edge A, nodal scalars, tree-cotree -- and claims complex-
symmetric matrices with losses and a unique solution. It is therefore a
known-good reference to measure against, which spike.py had nothing of.

THE FALSIFIABLE PREDICTION THIS TESTS. Their section 5.5 gives a one-line
recipe: "the computationally cheapest choice of gauge in Omega_N is to set all
FE coefficients x_psi associated with psi basis functions in the interior of
Omega_N to zero. In this case (26) still contributes to unknowns on Gamma."

Reading the blocks says why that is not an optimization but a necessity. A psi
degree of freedom strictly inside the insulator has, as omega -> 0:

    (2,1) = jw C_kappa^T  ->  jw * jw C_eps  ->  0
    (2,2) = jw G_kappa     ->  -w^2 G_eps     ->  0
    (2,3) = G_kappa        ->  jw G_eps       ->  0

-- its entire row vanishes, so the matrix is singular at DC. A psi degree of
freedom ON the interface touches conductor elements where kappa -> sigma, which
is finite, so its row survives. Hence:

    psi everywhere                       -> must degrade as omega -> 0
    psi zeroed inside Omega_N, kept on   -> must stay flat
    the interface Gamma

If the measurement contradicts that, the reading of section 13 is wrong and
GAUGE_CHOICE.md section 13 needs revisiting before anything is built.

CONVENTION. e^{+jwt}, so E = -jwA - grad(Phi) and B = curl A. Chew writes
e^{-iwt}; every i of his is -j here. Jochum writes e^{+ik0..} with a scaled
potential; his structure is reproduced in OUR variables rather than his, so the
equation numbers below refer to his paper but the symbols are ours.
"""
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from spike import (read_msh, build_edges, tet_geometry, boundary_edges,   # noqa
                   cond_of, LOC_EDGES, MU0, EPS0)

FREQS = (1e0, 1e2, 1e4, 1e6, 1e8, 1e10)
SIGMA, EPS_R = 5.8e7, 4.5
CONDUCTOR, DIELECTRIC = 1, 2

C0 = 1.0 / math.sqrt(MU0 * EPS0)          # speed of light
ETA0 = math.sqrt(MU0 / EPS0)              # characteristic impedance, ~377 ohm

# JOCHUM'S NON-DIMENSIONALIZATION, AND WHY IT IS NOT OPTIONAL.
#
# The first run of this script transcribed his block structure into our SI
# variables and reported a condition number RISING from 8.8e6 at 1 Hz to 3.2e23
# at 10 GHz -- backwards, and nothing to do with his formulation. The cause was
# the transcription: his A is a SCALED potential and his material coefficients
# are dimensionless, and that scaling is part of the stabilization, exactly the
# trap Balian et al. warn about ("must be applied on the analytical level").
#
# Substituting A = A_tilde / c0 into our Ampere equation
#
#     curl(nu curl A) + kappa (jw A + grad Phi) = J,     kappa = sigma + jw eps
#
# and multiplying through by mu0 c0 = eta0 gives
#
#     nu_r curl curl A_tilde + (sigma_0 + j k0 eps_r)(j k0 A_tilde + grad Phi)
#         = eta0 J
#
# with nu_r = mu0/mu, eps_r = eps/eps0, sigma_0 = eta0 sigma and k0 = w/c0 --
# which IS his (18a). Over this sweep k0 spans 2e-8 to 2e2 where w spans 6 to
# 6e10, so the operator's own dynamic range falls by eight orders before any
# matrix is built.
def k0_of(freq):
    return 2.0 * math.pi * freq / C0


# ----------------------------------------------------------------- tree ----
def spanning_tree(nodes, tets, edges, bedges):
    """Tree edge indices, with every boundary node contracted to one vertex.

    Stysch section 6.6.1, citing Klis 2015: "for each disconnected part Gamma_i
    of Gamma_el union Gamma_a, all mesh nodes in Gamma_i must be mapped onto one
    single vertex v_i in the graph". Here the whole outer boundary is one
    component (a box), so one super-vertex.

    Also from there, the ordering rule: span the INSULATOR first, then add trees
    inside each conductor, "to avoid numerical artifacts" at high frequency when
    the skin effect empties the conductor interiors.
    """
    eset = {e: k for k, e in enumerate(edges)}
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

    # Which material(s) touch each edge, so the insulator can be spanned first.
    touches = {}
    for t, mat in zip(tets, MATS_GLOBAL):
        for a, b in LOC_EDGES:
            key = (min(t[a], t[b]), max(t[a], t[b]))
            touches.setdefault(key, set()).add(mat)

    insulator_first = sorted(
        range(len(edges)),
        key=lambda k: 0 if touches.get(edges[k], set()) == {DIELECTRIC} else 1)

    tree = []
    for k in insulator_first:
        u, v = edges[k]
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
            tree.append(k)
    return sorted(tree)


# ------------------------------------------------------------- assembly ----
def assemble(nodes, tets, mats, kappa_of, eps_of, nu_of):
    """Every element matrix both formulations need, in one pass.

        S      = int nu (curl w_m).(curl w_n)          edge x edge
        M_kap  = int kappa w_m.w_n                     edge x edge
        C_kap  = int kappa w_m.grad(lam_n)             edge x node
        G_kap  = int kappa grad(lam_m).grad(lam_n)     node x node
        G_epsN = int_{Omega_N} eps grad.grad           node x node
        M_nod  = int lam_m lam_n                       node x node   (Chew's G)
        K2_eps = int eps w_m.w_n                       edge x edge   (Chew)
        C_eps  = int eps w_m.grad(lam_n)               edge x node   (Chew)
    """
    edges, eidx = build_edges(tets)
    Ne, Nn = len(edges), len(nodes)
    z = lambda r, c: np.zeros((r, c), dtype=complex)
    S, M_kap, K2_eps = z(Ne, Ne), z(Ne, Ne), z(Ne, Ne)
    C_kap, C_eps = z(Ne, Nn), z(Ne, Nn)
    G_kap, G_epsN, M_nod = z(Nn, Nn), z(Nn, Nn), z(Nn, Nn)

    for t, mat in zip(tets, mats):
        g, V = tet_geometry(nodes, t)
        kap, eps, nu = kappa_of(mat), eps_of(mat), nu_of(mat)
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
                S[gi[m], gi[n]] += s * nu * V * np.dot(curls[m], curls[n])
                w_w = (mass(am, an) * np.dot(g[bm], g[bn])
                       - mass(am, bn) * np.dot(g[bm], g[an])
                       - mass(bm, an) * np.dot(g[am], g[bn])
                       + mass(bm, bn) * np.dot(g[am], g[an]))
                M_kap[gi[m], gi[n]] += s * kap * w_w
                K2_eps[gi[m], gi[n]] += s * eps * w_w
            for n in range(4):
                v = (V / 4.0) * (np.dot(g[bm], g[n]) - np.dot(g[am], g[n]))
                C_kap[gi[m], t[n]] += sg[m] * kap * v
                C_eps[gi[m], t[n]] += sg[m] * eps * v
        for i in range(4):
            for j in range(4):
                gg = V * np.dot(g[i], g[j])
                G_kap[t[i], t[j]] += kap * gg
                M_nod[t[i], t[j]] += mass(i, j)
                if mat != CONDUCTOR:
                    G_epsN[t[i], t[j]] += eps * gg
    return dict(edges=edges, S=S, M_kap=M_kap, C_kap=C_kap, G_kap=G_kap,
                G_epsN=G_epsN, M_nod=M_nod, K2_eps=K2_eps, C_eps=C_eps)


# -------------------------------------------------------------- Jochum ----
def jochum_system(A, k0, keep_a, keep_p, keep_v):
    """Jochum (18)-(19) as one symmetric block system.

        Ampere tested by w_c          -> the A_c row            (21)
        Ampere tested by grad(psi)    -> the psi row            (22)
        gauge / Gauss tested by V     -> the V row           (23),(27)

        | S + jk0 M_kap    jk0 C_kap    C_kap   | | a |
        | jk0 C_kap^T      jk0 G_kap    G_kap   | | p |  =  0
        | C_kap^T          G_kap        G_epsN  | | v |

    in the non-dimensional variables above, with kappa = sigma_0 + j k0 eps_r.

    The (3,3) block is eps-weighted over the INSULATOR ONLY: in Omega_C the
    third equation is the gauge div(kappa A) = 0, which carries no V, while in
    Omega_N it is Gauss's law, which does. One expression covers both because
    kappa = j k0 eps_r there.

    Symmetry is a consequence of where the test functions go -- Jochum section
    5.2 -- not of any scaling applied afterwards.
    """
    jw = 1j * k0
    a, p, v = keep_a, keep_p, keep_v
    ix = np.ix_
    return np.block([
        [(A["S"] + jw * A["M_kap"])[ix(a, a)], jw * A["C_kap"][ix(a, p)],
         A["C_kap"][ix(a, v)]],
        [jw * A["C_kap"].T[ix(p, a)], jw * A["G_kap"][ix(p, p)],
         A["G_kap"][ix(p, v)]],
        [A["C_kap"].T[ix(v, a)], A["G_kap"][ix(v, p)], A["G_epsN"][ix(v, v)]],
    ])


# ---------------------------------------------------------------- Chew ----
def chew_schur(A, omega, chi_of, mats, keep_e, keep_n):
    """Chew's section 7 block system, reduced by the Phi block.

    Reproduces spike.py's measurement so the two formulations are compared on
    identical element matrices rather than across two scripts. K_NE carries
    (eps/chi) per element, which spike.py section 11 found is what forces chi to
    be position-dependent.
    """
    Kuu = (A["S"] - omega ** 2 * A["K2_eps"])[np.ix_(keep_e, keep_e)]
    KNE = -A["C_eps_over_chi"].T
    S = Kuu - A["C_eps"][np.ix_(keep_e, keep_n)] @ np.linalg.solve(
        A["M_nod"][np.ix_(keep_n, keep_n)], KNE[np.ix_(keep_n, keep_e)])
    return S


def sym_residual(M):
    """Relative asymmetry. ~1e-16 means complex-symmetric to round-off."""
    d = np.abs(M - M.T).max()
    s = np.abs(M).max()
    return d / s if s > 0 else 0.0


def equilibrate(M, sizes):
    """Balance the diagonal blocks by a CONGRUENCE: D M D, D block-constant.

    WHY THIS IS NECESSARY AND WHY IT IS LEGITIMATE.
    The raw block system is hopelessly scaled, and it is our doing, not the
    formulation's. With h ~ 1.5e-3 m the curl-curl block scales as 1/h ~ 670
    while the (2,2) block scales as k0 * kappa * h, which is ~1e-8 at 1 Hz --
    eleven orders between two diagonal blocks of the same matrix. That is the
    identical trap spike.py fell into on its first run, where chi ~ 1e-28 put
    the (2,2) block 1e40 below the (1,1) and the run reported 1e46.

    Balian et al. 2023 (GAUGE_CHOICE.md section 11.1) prescribe exactly this
    remedy and the form it must take: scale the equation AND the unknown by the
    same factor, a_k = b_k, so the operation is a congruence D M D and complex
    symmetry is preserved exactly. Here d_k = 1 / sqrt(||M_kk||), which gives
    every diagonal block unit norm. sym_residual is reported after the scaling
    precisely so that claim is checked rather than assumed.

    This is a statement about conditioning the linear algebra, not about the
    physics: a congruence leaves the generalized eigenvalue problem, and hence
    the solution, unchanged.

    Frobenius rather than spectral norms: equilibration only needs the right
    order of magnitude, and a spectral norm per block costs an extra SVD per
    block per frequency, which dominated the runtime of this script.
    """
    d = np.ones(M.shape[0])
    off = 0
    for n in sizes:
        if n == 0:
            continue
        blk = M[off:off + n, off:off + n]
        nrm = np.linalg.norm(blk) / math.sqrt(n)      # RMS entry scale
        d[off:off + n] = 1.0 / math.sqrt(nrm) if nrm > 0 else 1.0
        off += n
    return (M * d) * d[:, None]


# ---------------------------------------------------------------- main ----
MATS_GLOBAL = None


def main():
    global MATS_GLOBAL
    msh = os.path.join(HERE, "embedded.msh")
    if not os.path.isfile(msh):
        raise SystemExit("run:  gmsh embedded.geo -3 -o embedded.msh")
    nodes, tets, mats = read_msh(msh)
    MATS_GLOBAL = mats
    edges, _ = build_edges(tets)
    bedges = set(boundary_edges(nodes, tets, edges))
    bnodes = set()
    for k in bedges:
        bnodes.update(edges[k])

    tree = set(spanning_tree(nodes, tets, edges, sorted(bedges)))
    cotree_int = np.array(sorted(k for k in range(len(edges))
                                 if k not in tree and k not in bedges))
    int_nodes = np.array(sorted(n for n in range(len(nodes))
                                if n not in bnodes))
    print("mesh: %d nodes, %d tets (%d conductor, %d dielectric), %d edges"
          % (len(nodes), len(tets), (mats == CONDUCTOR).sum(),
             (mats == DIELECTRIC).sum(), len(edges)))
    print("tree-cotree: %d tree, %d boundary, %d interior cotree edges; "
          "%d interior nodes" % (len(tree), len(bedges), len(cotree_int),
                                 len(int_nodes)))

    # Which interior nodes sit ON the conductor/dielectric interface, and which
    # are strictly inside the insulator. Jochum section 5.5 keeps psi only on
    # the former.
    node_mats = [set() for _ in range(len(nodes))]
    for t, m in zip(tets, mats):
        for n in t:
            node_mats[n].add(m)
    on_gamma = np.array([n for n in int_nodes if len(node_mats[n]) > 1])
    in_cond = np.array([n for n in int_nodes
                        if node_mats[n] == {CONDUCTOR}])
    psi_jochum = np.array(sorted(set(on_gamma) | set(in_cond)))
    print("psi degrees of freedom: %d on the interface, %d inside the "
          "conductor, %d strictly inside the insulator (dropped)\n"
          % (len(on_gamma), len(in_cond),
             len(int_nodes) - len(psi_jochum)))

    # Dimensionless, per the note at the top of this file: sigma_0 = eta0 sigma,
    # eps_r, nu_r = 1, kappa = sigma_0 + j k0 eps_r.
    kap = lambda m, k0: (ETA0 * SIGMA + 1j * k0) if m == CONDUCTOR \
        else 1j * k0 * EPS_R
    eps_r = lambda m: 1.0 if m == CONDUCTOR else EPS_R
    nu_r = lambda m: 1.0
    eps_eff = lambda m, w: (EPS0 - 1j * SIGMA / w) if m == CONDUCTOR \
        else EPS0 * EPS_R

    def jochum_sweep(label, sigma, er_c, er_n, note=""):
        """One Jochum sweep. ALWAYS run with his own materials first.

        spike.py's first two runs declared Chew's formulation unusable, and both
        were wrong -- one was a scaling, one was a material he never used. The
        only thing that caught it was re-running his published setup unchanged.
        Same rule here.
        """
        print("\n--- %s" % label)
        if note:
            print("    %s" % note)
        print("    %-9s | %-10s | %-12s | %-12s | %-10s | %s"
              % ("freq", "k0", "cond (5.5)", "cond psi-all", "kap ratio",
                 "sym res"))
        print("    " + "-" * 76)
        kf = lambda m, k0: (ETA0 * sigma + 1j * k0 * er_c) if m == CONDUCTOR \
            else 1j * k0 * er_n
        ef = lambda m: er_c if m == CONDUCTOR else er_n
        for f in FREQS:
            k0 = k0_of(f)
            A = assemble(nodes, tets, mats, lambda m: kf(m, k0), ef, nu_r)
            M55 = equilibrate(
                jochum_system(A, k0, cotree_int, psi_jochum, int_nodes),
                (len(cotree_int), len(psi_jochum), len(int_nodes)))
            Mall = equilibrate(
                jochum_system(A, k0, cotree_int, int_nodes, int_nodes),
                (len(cotree_int), len(int_nodes), len(int_nodes)))
            r = abs(kf(CONDUCTOR, k0)) / abs(kf(DIELECTRIC, k0))
            print("    %-9.0e | %-10.2e | %-12.3e | %-12.3e | %-10.2e | %.1e"
                  % (f, k0, cond_of(M55), cond_of(Mall), r,
                     sym_residual(M55)))

    print("=" * 78)
    print("JOCHUM -- non-dimensional, psi restricted per his section 5.5")
    print("=" * 78)

    jochum_sweep("CONTROL: Jochum's own section 6.1 cavity materials",
                 sigma=1.0, er_c=1.0, er_n=1.0,
                 note="sigma_C = 1 S/m, eps_C = eps_N = 1 -- the run that "
                      "validates the assembly")

    jochum_sweep("OURS: copper against eps_r 4.5",
                 sigma=SIGMA, er_c=1.0, er_n=EPS_R,
                 note="sigma_C = 5.8e7, seven orders above his")

    print("\n    cond (5.5): psi on the interface and inside the conductor only.")
    print("    cond psi-all: psi on every interior node. The block reading says")
    print("    an insulator-interior psi row vanishes as k0 -> 0, so psi-all")
    print("    must be worse, and by a margin that does not close.")
    print("    kap ratio: |kappa_C| / |kappa_N|, the material contrast.")

    print("\n" + "=" * 78)
    print("CHEW -- same mesh, same element matrices, eps_eff for sigma")
    print("=" * 78)
    print("    %-9s | %-13s | %-13s | %s"
          % ("freq", "cond gauged", "cond ungauged", "eps contrast"))
    print("    " + "-" * 62)
    for f in FREQS:
        w = 2.0 * math.pi * f
        A = assemble(nodes, tets, mats, lambda m: kap(m, w),
                     lambda m: eps_eff(m, w), lambda m: 1.0 / MU0)
        # K_NE = -(eps/chi) C^T with chi = mu eps^2 per element, so eps/chi = nu/eps.
        A2 = assemble(nodes, tets, mats, lambda m: kap(m, w),
                      lambda m: (1.0 / MU0) / eps_eff(m, w),
                      lambda m: 1.0 / MU0)
        A["C_eps_over_chi"] = A2["C_eps"]
        Kuu = (A["S"] - w ** 2 * A["K2_eps"])[np.ix_(cotree_int, cotree_int)]
        Sc = chew_schur(A, w, None, mats, cotree_int,
                        np.arange(len(nodes)))
        e1, e2 = abs(eps_eff(1, w)), abs(eps_eff(2, w))
        print("    %-9.0e | %-13.3e | %-13.3e | %.3e"
              % (f, cond_of(Sc), cond_of(Kuu), max(e1, e2) / min(e1, e2)))

    print("\nNOTE: Chew is measured on the tree-cotree-reduced edge space here,")
    print("not on the full space spike.py used, so his numbers are not directly")
    print("comparable to section 11 of GENERALIZED_LORENZ_GAUGE.md -- only to")
    print("Jochum's above, which is the point of this run.")


if __name__ == "__main__":
    main()
