# 07_GaugeInvariance

**One problem, five solves, three spanning trees.** The only case in the suite
that compares solves against *each other* rather than against an analytic
answer — so there is no discretisation error in the comparison at all, and the
tolerances are round-off.

```bat
..\check.bat 07_GaugeInvariance
```

About a minute: five solves of a 19,686-unknown problem on the internal solver,
no MUMPS and no oneAPI needed.

## What it asserts

The A‑Φ formulation is invariant under

```
A  →  A + ∇ψ            Φ  →  Φ − jωψ
```

`B = ∇×A` is unchanged because the curl of a gradient is identically zero, and
`E = −jωA − ∇Φ` is unchanged because the two shifts cancel against each other.
The tree‑cotree gauge removes the freedom by forcing `A = 0` on a spanning
tree's edges, which picks **one** representative — so a different tree must give
the same fields and different potentials.

| | measured | asserted |
|---|---|---|
| Φ moves | 0.144, 0.153 | ≥ 0.02 |
| A moves | 0.199, 0.218 | ≥ 0.02 |
| E, B, H, J move | 4.7e‑13 … 3.2e‑12 | ≤ 1e‑9 |
| ψ on the outer boundary | 3.6e‑16 Wb | ≤ 1e‑12 |
| ψ in the interior | 4.6e‑04 Wb | ≥ 1e‑6 |
| terminals, voltage drive | 0.0 exactly | ≤ 1e‑15 |
| terminals, current drive | spread 4e‑22, ΔV 2.7e‑13 | ≤ 1e‑15, 1e‑9 |

**Both directions are asserted, and that is the point.** A test that only
checked invariance would still pass if the gauge freedom were accidentally
*removed* — by over‑constraining Φ, say — which is a real regression that would
sit happily behind a green check. So the potentials are required to move by at
least a floor and the fields by at most a ceiling, with eleven orders of
magnitude of daylight between the two.

## How different trees are obtained

By permuting the order the nodes are listed in the `.msh`. Nothing about that is
obvious from the file, so, in the solver's own terms:

- `Accumulator::add_node` assigns each node its internal index from
  `mesh.nodes.size()` — **in file order**. The gmsh tag is only a lookup key.
- `Mesh::build` sorts edges by their canonical vertex pair, so edge numbering
  follows node numbering.
- `build_tree_cotree` roots at the lowest‑numbered node of Dirichlet component 0
  and walks neighbours in edge‑index order.

So reordering the file moves the root and the traversal, giving a genuinely
different tree over a problem that is unchanged in every physical respect. The
permutation is applied to the (tag, coordinate) *pairs*, so each node keeps its
own tag and coordinates and `$Elements` needs no edit.

It is confirmed rather than assumed: the solver reports the same 19,686 unknowns
with **different stored nonzeros** — 406,808 / 406,866 / 406,828 — because
different edges were gauged away.

> Renumbering the *tags* in place would not work. The reader ignores tag values
> for ordering, so the file would look shuffled and yield the identical tree — a
> test that silently proves nothing.

## Why this geometry, and why it is coarse

Φ only moves appreciably away from the resistive limit: the gauge‑sensitive
fraction is about `ωL/|Z| = sin(arg Z)`, derived and measured in
[06_TwoWires_Quadrature's report §7](../06_TwoWires_Quadrature/TwoWiresValidation.md).
These proportions at 50 Hz give **`ωL/R = 2.84`**, where changing the tree moves
Φ by 15 % of the applied volt. A long thin wire at the same frequency would give
`ωL/R ≈ 0.07`, Φ would barely move, and every assertion would pass for the wrong
reason — so `verify_gauge.py` **measures the ratio and asserts the regime**
rather than trusting the geometry.

The mesh is deliberately loose. Nothing here is compared against an analytic
field, so discretisation error cancels exactly between the two solves being
differenced and a fine mesh would buy nothing but five slower solves. The skin
effect is `01_OneCylinder`'s job and the Bessel profile is `04`'s.

## Why the terminals survive, including under a current source

This is the part worth reading, because the obvious guess is wrong.

A voltage‑driven case pins Φ by Dirichlet at **both** terminals, forcing `ψ = 0`
at both, so the terminal drop cannot move — a weaker test than it looks. A
**current**‑driven case pins Φ at the reference terminal only; the driven
terminal's potential is a free unknown. If terminal voltage were going to move
with the gauge anywhere, it would be there.

It does not, and the reason is geometric rather than anything about port types:
`n × A = 0` forces `ψ` constant on the outer boundary, the reference port pins
it to zero, the box surface is connected, and the wire spans its full height — so
**`ψ = 0` on the whole outer boundary, both terminals included.** The case
measures `ψ` directly, recovered from the shift as `ψ = jΔΦ/ω`, and finds it
twelve orders of magnitude smaller on the boundary than in the interior. The
gauge distorts the *shape* of Φ and pins its endpoints.

> **The protection is geometric, so it has a limit.** An electrode floating in
> the *interior* of the domain would have no such pinning, and its potential
> would be gauge‑dependent. Every case in this suite puts its terminals on the
> boundary; a future one that does not would need this re‑examined.

## Files

| | |
|---|---|
| `gauge.geo`, `gauge.msh` | a short fat copper rod in an air box, 8,365 tets |
| `gauge.aphi` | the voltage drive, 1 V at 50 Hz |
| `gauge_1a.aphi` | the current drive, 1 A — same geometry |
| `permute_nodes.py` | reorders a 4.1 `.msh`'s node list; geometrically a no‑op |
| `verify_gauge.py` | runs all five solves, compares, asserts, cleans up |
| `expected.txt` | the thresholds, and what sets each one |

`verify_gauge.py` performs its own solves rather than relying on the one
`check.bat` runs, so it behaves identically under `check.bat 07_GaugeInvariance
verify-only`. It writes into `g_*` directories and removes them on the way out.
