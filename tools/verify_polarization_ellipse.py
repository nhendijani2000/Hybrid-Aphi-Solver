"""Check the closed forms in docs/ComplexVectorPhasorConcept.md.

    pvpython verify_polarization_ellipse.py        (needs only numpy)

A vector phasor Ehat = P + jQ has E(u) = P cos u - Q sin u, whose tip traces an
ellipse. The document derives its semi-axes in closed form:

    S = (|P|^2+|Q|^2)/2 ,  C = (|P|^2-|Q|^2)/2 ,  D = P.Q ,  R = sqrt(C^2+D^2)
    a = sqrt(S+R)   semi-major (the peak of |E(t)|)
    b = sqrt(S-R)   semi-minor          <- but see below
    N = sqrt(|P|^2+|Q|^2) = sqrt(a^2+b^2)    the stored "complex magnitude"

Closed forms are derived; a brute-force sweep over u is not. If they disagree,
the derivation is wrong. That is what makes section 15 of the document a
measurement rather than an assertion.

TWO WAYS TO GET b, AND ONLY ONE OF THEM WORKS. sqrt(S-R) cancels catastrophically
when S ~ R, which is precisely the nearly-linear case -- and that is most of the
domain in every case in this suite. The product invariant a*b = |PxQ| gives

    b = |PxQ| / a

with a well conditioned, so this one is stable. Test 4 below exists to exhibit
the cancellation; every other test uses the stable form, because that is what
the document tells an implementer to write.
"""
import numpy as np

TH = np.linspace(0.0, 2.0 * np.pi, 200001)
rng = np.random.default_rng(20261001)


def closed(P, Q):
    """(a, b, N) with b from the STABLE product invariant."""
    N2 = P @ P + Q @ Q
    S, C, D = 0.5 * N2, 0.5 * (P @ P - Q @ Q), P @ Q
    a = np.sqrt(S + np.hypot(C, D))
    b = np.linalg.norm(np.cross(P, Q)) / a if a > 0 else 0.0
    return a, b, np.sqrt(N2)


def b_cancelling(P, Q):
    """The obvious formula, kept only to measure how badly it loses digits."""
    S = 0.5 * (P @ P + Q @ Q)
    R = np.hypot(0.5 * (P @ P - Q @ Q), P @ Q)
    return np.sqrt(max(S - R, 0.0))


def brute(P, Q):
    """(max, min) of |E(u)| by sweeping u. No closed form used."""
    m = np.linalg.norm(np.outer(np.cos(TH), P) - np.outer(np.sin(TH), Q), axis=1)
    return m.max(), m.min()


def draws(n, force_linear=0.25, scales=(1.0, 1e-3, 1e3)):
    for _ in range(n):
        P = rng.normal(size=3) * rng.choice(scales)
        Q = rng.normal(size=3) * rng.choice(scales)
        if rng.random() < force_linear:
            Q = P * rng.normal()          # exactly linear: |PxQ| = 0, b = 0
        yield P, Q


print("verifying docs/ComplexVectorPhasorConcept.md\n")

ea = eb = 0.0
for P, Q in draws(4000):
    a, b, _ = closed(P, Q)
    am, bm = brute(P, Q)
    ea = max(ea, abs(a - am) / am)
    eb = max(eb, abs(b - bm) / am)
print("1. a vs sweep maximum                        %.1e" % ea)
print("2. b vs sweep minimum                        %.1e  <- grid-limited. |E| is"
      % eb)
print("                                                     quadratic near its min,")
print("                                                     so a discrete sweep")
print("                                                     resolves it to ~one step.")

e1 = e2 = ec = 0.0
for P, Q in draws(20000):
    a, b, N = closed(P, Q)
    e1 = max(e1, abs(a * a + b * b - N * N) / (N * N))
    # normalise by N^2, NOT by |PxQ|: that is exactly 0 for the linear draws,
    # and dividing by it reports a meaningless 1e8 from pure roundoff.
    # Test 4 must use the INDEPENDENTLY derived b = sqrt(S-R). Using the
    # stable b would make it a tautology -- b is defined as |PxQ|/a there, so
    # a*b = |PxQ| would hold by construction and the test would check nothing.
    e2 = max(e2, abs(a * b_cancelling(P, Q) - np.linalg.norm(np.cross(P, Q))) / (N * N))
    ec = max(ec, abs(b_cancelling(P, Q) - b) / N)
print("3. identity  a^2 + b^2 = N^2  (stable b)     %.1e  <- genuine: a from S+R," % e1)
print("                                                     b from |PxQ|/a, and")
print("                                                     a^4 - N^2 a^2 + |PxQ|^2 = 0")
print("4. identity  a*b = |PxQ|  (b from sqrt(S-R)) %.1e  <- sqrt(eps) = 1.5e-8," % e2)
print("5. the two b formulas against each other     %.1e     the cancellation that" % ec)
print("                                                     section 11 warns about.")

e = 0.0
for P, Q in draws(2000):
    _, _, N = closed(P, Q)
    E = np.outer(np.cos(TH[:-1]), P) - np.outer(np.sin(TH[:-1]), Q)
    e = max(e, abs((np.linalg.norm(E, axis=1) ** 2).mean() - N * N / 2) / (N * N / 2))
print("6. <|E|^2> over a cycle = N^2 / 2            %.1e" % e)

e = 0.0
for P, Q in draws(2000, force_linear=0.0, scales=(1.0,)):
    C, D = 0.5 * (P @ P - Q @ Q), P @ Q
    dl = -0.5 * np.arctan2(D, C)
    Pp = P * np.cos(dl) - Q * np.sin(dl)
    Qp = P * np.sin(dl) + Q * np.cos(dl)
    a, b, N = closed(P, Q)
    hi, lo = sorted((np.linalg.norm(Pp), np.linalg.norm(Qp)))[::-1]
    e = max(e, abs(Pp @ Qp) / (N * N), abs(hi - a) / N, abs(lo - b) / N)
print("7. delta = -psi/2 orthogonalises, and")
print("   |P'|, |Q'| are the semi-axes              %.1e" % e)

bad = 0
for P, Q in draws(4000):
    a, b, N = closed(P, Q)
    am, bm = brute(P, Q)
    if not (a <= N * (1 + 1e-9) and N <= np.sqrt(2) * a * (1 + 1e-9)):
        bad += 1
    if not (bm >= b * (1 - 1e-6) - 1e-9 * N and am <= a * (1 + 1e-9)):
        bad += 1
print("8. a <= N <= sqrt(2) a  and  b <= |E| <= a   %d violations in 4000 draws" % bad)

print()
print("   the three worked examples of section 8")
print("   %-22s %8s %8s %8s %8s %8s" % ("", "a", "b", "N", "b/a", "N/a"))
for name, E in (("linear      (3,0,0)", np.array([3.0, 0, 0], dtype=complex)),
                ("elliptical  (3,2j,0)", np.array([3.0, 2.0j, 0])),
                ("circular    3(x+jy)", np.array([3.0, 3.0j, 0]))):
    a, b, N = closed(E.real.copy(), E.imag.copy())
    print("   %-22s %8.4f %8.4f %8.4f %8.4f %8.4f" % (name, a, b, N, b / a, N / a))
