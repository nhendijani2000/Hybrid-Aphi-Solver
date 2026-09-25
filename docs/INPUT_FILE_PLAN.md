# Input file preparation — plan

> **Status, 25 Sept 2026.** Steps 1–6 of §6 are **built and tested**: the
> `Problem` struct, the parser, the cylinder mesh, binding, conduction paths
> / gauge / cut side labels, and the DOF map. Roughly 1280 checks across
> eleven test executables. What remains is element matrices, global
> assembly, and the solve — which are Phase 04 of `docs/ROADMAP.md`, not
> this document.
>
> Read the dated revision notes below as a record of how the design was
> arrived at, including two mistakes that were caught and corrected (the
> `port_sign` that turned out to be identically +1, and the "exactly one
> potential reference" rule that rejected a voltage-driven resistor). The
> §7 decisions that remain open are noted there.
>
> *Moved into the repository 25 Sept, from `Claude outputs/` beside it,
> which is not under version control — the code comments referencing it
> were pointing outside the repo.*

*Opus 5.5, 23 Sept 2026. Plan for review; nothing here is implemented yet.
Follows the decisions of 23 Sept: internal ports are grounded on one side
(no floating ports for now); Φ is left free on the domain box; the input is a
plain text file organised in sections; the cylinder is the first test case.*

*Revised (Opus 5, 23 Sept 2026) after review — six changes, all recorded in
place: the internal-cut rim check is now regime-dependent and an **error** at
DC (§3.2); `L` has a stated tolerance and the faceting correction is applied
to it as well as to `R` (§5.1); an **ω > 0 skin-effect stage** is added to
the cylinder test (§5.2), which is what actually exercises the "Φ free on the
box" decision; the deliberate omissions are written down rather than left to
be rediscovered (§1.3); the `type`/`formulation` interaction is documented
(§1.2); and the Gmsh version is pinned (§5.3). Then revised again the same
day: materials became **`[Body …]` sections, one per mesh volume** (§1.2),
and §3.3 now names the **containers binding must build**, since `Mesh`
itself has no bodies and no surfaces. The complex-amplitude grammar was then
settled as **magnitude plus `phase_deg`** (§1.2), which closes the format:
see the status note at the end of §7. Finally, 24 Sept: **`direction` and
`ground` were removed from the input file** and are derived at bind time.
That makes the port sign identically +1 and retires `port_sign` altogether
(§2.1), at the cost of an arbitrary — though deterministic — sign on an
internal port's reported I and V. Later the same day that cost was bought
back: **`current_direction` returns as a *hint*** — any vector within 90° of
the intended flow, since only its sign is used — so a user can say "roughly
+y" without ever computing a normal, and the dry run prints the resolved
vector back for checking (§2.1, §3.3). §2.4 records what went and why.*

The pipeline this plan covers ends at a **validated problem bound to the
mesh**. DOF assignment is the next plan, and consumes what this one produces.

```
input file ──parse──▶ Problem ──bind to mesh──▶ bound problem ──▶ (next: DOF map)
                         ▲                            ▲
                    syntax + semantic          names, faces, sides,
                      validation               topology checks
```

---

## 1. Format

A plain text file, one section per topic, `key = value` inside each.

```ini
# Cylinder in a square box -- DC resistance and inductance.

[mesh]
file        = ../meshes/cylinder_box.msh   # relative to this input file
length_unit = mm                           # m | mm | um | nm

[analysis]
type        = dc                           # dc | frequency
# frequencies = 1e3 1e6 1e9                # Hz; required when type = frequency
# formulation = full_wave                  # full_wave | reduced

[boundary]
outer       = flux_tangential              # n x A = 0 on the box, Phi free

[Body B1]
volume      = wire                         # Physical Volume name (or tag)
sigma       = 5.8e7                        # S/m
eps_r       = 1
mu_r        = 1

[Body B2]
volume      = air
sigma       = 0

[port P1]
type        = boundary_current
surface     = wire_bottom                  # Physical Surface name (or tag)
current     = 1.0                          # A, peak
# phase_deg = 0

[port P2]
type        = boundary_voltage
surface     = wire_top
voltage     = 0.0                          # V, peak; this port is the 0 V reference
```

The same problem swept over frequency — stage 2 of §5.2. Only `[analysis]`
differs from the file above; same mesh, same bodies, same ports:

```ini
# Cylinder in a square box -- skin effect vs frequency.

[mesh]
file        = ../meshes/cylinder_box.msh
length_unit = mm

[analysis]
type        = frequency
frequencies = 1e3 1e4 1e5 1e6 1e7      # Hz -- delta/a = 10.4 ... 0.105
formulation = full_wave                # or: reduced, which must agree here

[boundary]
outer       = flux_tangential

[Body B1]
volume      = wire
sigma       = 5.8e7
eps_r       = 1
mu_r        = 1

[Body B2]
volume      = air
sigma       = 0

[port P1]
type        = boundary_current
surface     = wire_bottom
current     = 1.0
phase_deg   = 0

[port P2]
type        = boundary_voltage
surface     = wire_top
voltage     = 0.0
```

One solve per frequency, with the same amplitudes and phases at each.

An internal port, e.g. across a loop:

```ini
[port P1]
type              = internal_current
surface           = loop_cut
current_direction = +y        # roughly which way positive current flows
current           = 1.0
```

`current_direction` is **a hint, not the normal** — see §2.1. The cut's
normal comes from its own faces; the hint only resolves the remaining ±
sign, so any vector within 90° of what you mean gives the same answer. For a
ring in the x–y plane cut at the half-plane `y = 0, x > 0`, the cut's normal
is ±ŷ, so `+y` says "counterclockwise seen from +z" — a statement you can
make by looking at your own `.geo`, with nothing to compute.

Omitting it is legal for a single internal port: the direction is then
derived and the port's reported `I` and `V` carry an arbitrary sign (`Z`,
`R`, `L`, `P` do not). With two or more internal ports it is required.

### 1.1 Grammar

| rule | detail |
|---|---|
| comments | `#` to end of line, anywhere |
| section header | `[kind]` or `[kind name]` |
| entries | `key = value`; whitespace around `=` ignored |
| lists | whitespace-separated: `surface = cut1 cut2`, `frequencies = 1e3 1e4` |
| numbers | plain or scientific (`5.8e7`) |
| case | section kinds, keys and keyword values are case-insensitive; **names are case-sensitive** (they must match Gmsh physical names) |
| repetition | `[mesh]`, `[analysis]`, `[boundary]` once each; `[Body …]` and `[port …]` repeatable, names unique |
| unknown section or key | **error**, naming file and line, listing the valid keys — a typo must never be silently ignored |
| reserved sections | `[solver]`, `[output]`: recognised, rejected with "not supported yet" until they do something |

### 1.2 Sections and keys

**`[mesh]`**

| key | required | values |
|---|---|---|
| `file` | yes | path, relative to the input file |
| `length_unit` | yes | `m`, `mm`, `um`, `nm` — Gmsh files carry no units |

**`[analysis]`**

| key | required | values |
|---|---|---|
| `type` | yes | `dc` or `frequency` |
| `frequencies` | see below | an explicit list, Hz, all > 0; forbidden for `dc` |
| `sweep` | see below | `log` (default) or `linear` |
| `f_start`, `f_stop` | see below | Hz, both > 0, `f_stop` > `f_start` |
| `f_points` | see below | integer ≥ 2 — **total** points, both endpoints included |
| `formulation` | no | `full_wave` (default) or `reduced`; `frequency` only |

**Two ways to give the frequencies, decided 23 Sept.** `type = frequency`
needs exactly one of them:

```ini
frequencies = 1e3 1e4 1e5 1e6 1e7     # explicit list
```
```ini
sweep    = log                         # log | linear
f_start  = 1e3
f_stop   = 1e7
f_points = 41                          # total, endpoints included
```

Giving both is an error. Both forms are kept because they do different
jobs: the explicit list pins the exact points a regression test asserts (and
lets one awkward frequency be reproduced without recomputing a sweep),
while the sweep form is what any real run wants — writing 41 values by hand
is intolerable.

**The parser expands a sweep into the explicit list immediately**, so
`Problem` only ever carries `frequencies` and no code downstream of the
parser learns a second way to say the same thing.

`type` and `formulation` are not independent — they jointly decide where Φ
lives, which in turn decides the rim check in §3.2. The valid combinations,
stated here so the parser can reject the rest rather than silently pick one:

| `type` | `formulation` | Φ lives on | note |
|---|---|---|---|
| `dc` | *(forbidden)* | conductors (σ > 0) | at ω = 0 the Φ equation `−∇·(σ∇Φ) = 0` is identically empty where σ = 0, so Φ DOFs in an insulator would be zero rows. Not a choice — forced |
| `frequency` | `full_wave` (default) | whole domain | ε retained |
| `frequency` | `reduced` | conductors (σ > 0) | ε → 0, the eddy-current variant of `FORMULATION.md` §3 |

Giving `formulation` alongside `type = dc` is an error, not a silent
override: it almost always means the user believed they had a choice there.

**`[boundary]`**

| key | required | values |
|---|---|---|
| `outer` | no | `flux_tangential` (default): `n×A = 0`, Φ free. `pec`, `abc` recognised and rejected as not yet supported |

`outer` is a single global setting: every face of the domain boundary gets
the same treatment, found from the mesh topology (`boundary_edge_mask`, which
needs no tags). Per-surface conditions — a PEC ground plane, an ABC on three
faces and a PEC floor on the fourth — will need `[boundary <surface-name>]`
sections. The grammar should leave that spelling free now (a `[boundary …]`
section *with* a name is reserved, not repurposed), so adding it later does
not invalidate files written today.

**`[Body name]`** — **one section per mesh volume.** The section name (`B1`,
`wire_body`, whatever) is the user's own handle for that body; `volume` is
what binds it to the mesh, and must match a Physical Volume exactly.

| key | required | values |
|---|---|---|
| `volume` | yes | exactly one Physical Volume name or tag |
| `sigma` | yes | S/m, ≥ 0 — required, so a conductor can never default to an insulator |
| `eps_r` | no | > 0, default 1 |
| `mu_r` | no | > 0, default 1 |

Decided 23 Sept. Two alternatives were considered and rejected:

- *`[material copper]` with a `volumes = wire via1 via2` list* — one material
  definition shared by several volumes. Rejected: the indirection buys
  nothing at this scale, and the body list is more useful read directly as
  "what the model is made of". Several volumes of one material simply get
  several `[Body …]` sections.
- *`material = 5.8e7 1` as a positional σ/ε pair* — rejected in favour of
  named keys: `mu_r` has nowhere to go in a pair, later additions (loss
  tangent, frequency-dependent σ, anisotropy) extend cleanly as new keys but
  break a positional list, and a σ/ε swap in a positional pair is exactly
  the kind of error that solves quietly and wrongly.

**`[port name]`**

| key | boundary ports | internal ports | values |
|---|---|---|---|
| `type` | yes | yes | `boundary_current`, `boundary_voltage`, `internal_current`, `internal_voltage` |
| `surface` | yes | yes | one or more Physical Surface names or tags |
| `current` | current ports | current ports | magnitude, A (peak) |
| `voltage` | voltage ports | voltage ports | magnitude, V (peak) |
| `phase_deg` | no | no | degrees, default 0 — applies to `current` or `voltage`, whichever the port has |

| `current_direction` | forbidden | optional (required if ≥ 2 internal ports) | axis shorthand (`+x`, `-y`, …) or a 3-vector; a **hint**, see below |

A current port given `voltage`, or a voltage port given `current`, is an
error — the wrong key for the type almost always means the wrong type.

**`current_direction` is a hint, not a normal.** The cut's normal `n̂` is
computed from its faces; the hint `v` only resolves the sign:

```
d = sign(v · n̂) · n̂
```

so `+y`, `0 1 0` and `0.1 0.9 -0.2` are all the same input. A user never has
to produce a mesh normal — which is the point, since that is not something
anyone can write down by hand. Rejected if `|v̂·n̂| < 0.1`: the hint then lies
nearly *in* the cut plane and picks no side.

Absent, the direction is derived (§3.3) and the sign of that port's reported
`I` and `V` becomes arbitrary. That is a **warning** for one internal port
and an **error** for two or more, where the *relative* signs between ports
would otherwise make any coupling term meaningless.

**There is still no ground key** (removed 24 Sept; see §2.4) — the grounded
side follows from the direction. And boundary ports take no direction at
all: inward is positive, and the mesh already knows which way that is.

**Complex amplitude — decided 23 Sept: magnitude plus `phase_deg`, for
current and voltage alike.**

`phase_deg` is **one key that applies to whichever amplitude the port
carries** — there is no `current_phase_deg`/`voltage_phase_deg` pair, since
a port has exactly one excitation:

```ini
[port P1]                          [port P2]
type      = boundary_current       type      = boundary_voltage
surface   = wire_bottom            surface   = feed
current   = 1.0    # A, peak       voltage   = 3.3   # V, peak
phase_deg = 0      # degrees       phase_deg = 90    # degrees
```

In both cases the driven complex amplitude is
`A · exp(j · phase_deg · π/180)`, under this project's `e^{+jωt}` convention
(`FORMULATION.md` §1.1) — so a positive `phase_deg` *leads*. The port sign
`s` of §2 multiplies that value; it is not folded into the phase.

Rules that follow, identically for both:

- One amplitude and phase per port, applied **at every frequency** of a
  sweep. A frequency-dependent excitation is not expressible and is not
  needed by anything planned.
- A negative magnitude is allowed, not an error: `−1.0` at 0° and `1.0` at
  180° are the same excitation, and forbidding one of them would only
  surprise people.
- At `type = dc`, `phase_deg` must be absent or `0`. A nonzero phase at DC
  is an error — there is no phase to have, so the file says something its
  author cannot have meant. This is why the cylinder example omits the key
  entirely.
- `phase_deg` without a `current` or `voltage` to attach to is an error, as
  is giving it on a port whose amplitude key is missing.

*Rejected: a re/im pair (`current = 1.0 0.0`).* It is marginally less
ambiguous to machine-read, but magnitude-and-phase is how the quantity is
actually thought about at a port, and keeping every value in the file a
single scalar keeps the parser's value grammar to "number, keyword, name, or
3-vector" -- in fact, with `direction` since removed, to "number, keyword
or name" alone.

### 1.3 Deliberately not in this format yet

Written down so these read as decisions rather than oversights, and so they
are not rediscovered later as gaps:

| omitted | why | when |
|---|---|---|
| **gauge choice** and **frequency scaling** — the `(gauge, frequency_scaling)` pair `ROADMAP.md` Phase 04 step 6 wants exposed | both are solver-conditioning knobs, not problem description. Method A and one scaling are hard-wired until there is a validated solve to compare against; exposing a knob with one tested setting invites turning it | `[solver]`, once Phase 05 has something to compare |
| **volumetric `J_imp`** — `ROADMAP.md` Phase 04 step 7's "simple uniform current density" | superseded. That step predates the port abstraction and exists only to have *some* right-hand side; ports give a better-posed excitation with an analytic check attached. Step 7 should be re-pointed at ports rather than implemented as written | not planned |
| **output destination / field export** | the first milestones are scalar checks (`R`, `L`, `Z(ω)`) printed by the CLI. Field export needs a format decision (VTU) that nothing yet depends on | `[output]`, when there are fields worth looking at |
| **per-surface boundary conditions** | see `[boundary]` above — the spelling is reserved | with PEC/ABC |

---

## 2. Orientation and sign conventions

**Revised 24 Sept 2026 — there is now no sign factor at all.** The earlier
design let the input file ground either side of a cut, which needed a
per-port `s = ±1` multiplying both the prescribed potential and the port
row. Grounding the minus side *always* makes `s ≡ +1`, so `ground`,
`GroundSide` and `port_sign` were removed rather than left as a branch no
input could reach. The only thing that choice ever controlled was Φ's
arbitrary offset — never `V`, `I`, `Z` or `P`. The history is kept in §2.4.

### 2.1 Internal ports

Every cut has a direction **d**, and it fixes everything:

- the **minus** side is the one **d** emits from; it is **held at 0 V**;
- the **plus** side is the one **d** points into; it carries the port's
  single unknown potential;
- `I` = current crossing the cut along **d**;
- `V` = Φ(plus) − Φ(minus) = Φ(plus).

So `P = ½ Re(V I*)` is the power the source delivers: current runs − → +
inside the gap, then round the conductor from + back to −, which is the
battery convention. A current port puts `I` straight into its row's
right-hand side; a voltage port prescribes `V` on the plus side and reads
`I` back as the row residual. Nothing is multiplied by anything.

**How d is fixed (revised 24 Sept — the user gives a hint).** The cut's
normal `n̂` comes from its own faces, which are planar by §3.2's check. Only
the ± sign is left, and that is what `current_direction` supplies:

```
d = sign(v · n̂) · n̂          (resolve_cut_direction)
```

The essential property: **the hint is used only for its sign**, so it need
not be the normal — anything within 90° of the intended direction gives an
identical `d`. `+y`, `0 1 0` and `0.1 0.9 -0.2` are the same input. That is
what makes this answerable by a user, who can see from their own geometry
which way current should go but cannot compute a mesh normal. Rejected when
`|v̂·n̂| < 0.1`, where the hint lies nearly in the cut plane and picks no
side.

**With no hint**, `canonical_orientation` picks a deterministic sign — the
normal oriented so its largest-magnitude component is positive, ties broken
in x, y, z order. Arbitrary with respect to intent, but *reproducible*,
which is the only requirement on it: a regression test's reported sign must
not drift between builds on one unchanged mesh. The consequence must then be
reported with the result:

> The **sign** of that port's `I` and `V` is arbitrary — mesh-determined,
> though stable for a given mesh. `Z = V/I`, `R`, `L`, `|I|` and `P` are
> unaffected, because `I` and `V` flip together.

Hence the escalation: **warn** when one internal port has no hint, **error**
when two or more do. Relative signs between ports are not arbitrary in any
useful sense — a mutual term computed from them would be meaningless — so
that case is worth refusing rather than warning about.

**The grounded side follows from d**, so it too is arbitrary without a hint.
Harmless for a floating loop, which has no other reference; but a conductor
carrying both a grounded cut side and a 0 V boundary port has two references
and is over-constrained (§3.2, error 9).

*Deferred, not solved:* `current_toward = <point or named surface>` reads
better for curved geometry ("current flows toward the output pad"), but
needs a landmark unambiguously on one side, which fails for a loop where
both sides of the cut are the same body. Worth adding later as a second
spelling. A bare `reverse = true` was rejected outright: it records no
intent, so a file read a month later cannot be understood without running
it.

### 2.2 Boundary ports

- `I` = current **entering** the domain through the terminal, i.e. along the
  inward normal. That normal is derivable from the mesh (a boundary face has
  exactly one adjacent tet, which is the inside), so **nothing is written in
  the file** — supplying a direction would create a second, possibly
  contradicting, source of truth.
- `V` = Φ(terminal) − 0.
- `P = ½ Re(V I*)` is again the power delivered.

Unlike an internal port, a boundary port's sign is **not** arbitrary: inward
is positive, and the mesh knows which way that is.

A `boundary_voltage` port with `voltage = 0` is the ground terminal. In the
cylinder test, P1 injects +1 A at the bottom and P2 holds the top at 0 V; the
current read back from P2 must be **−1 A** (1 A leaving). That is the first
self-consistency check the solver will have.

### 2.3 Recorded limitation

Grounding a whole side of an internal cut makes both sides equipotential,
which enforces the *total* current across the cut but not the pointwise
current continuity the tree gauge relies on (review of 22 Sept, §3c). It is
exact when the cut is a true equipotential — a straight conductor, or a
ring cut along a radial plane, which is the planned loop test. Kept as a
known limitation, to be measured once a bent conductor is tested.

### 2.4 What was removed, and why it is not coming back

For the record, since a `±1` is exactly the kind of thing that grows back
once someone notices a convention could have gone the other way:

| removed | was |
|---|---|
| `ground = minus \| plus` | an input key choosing which side of a cut sat at 0 V |
| `GroundSide` | its enum |
| `port_sign(const Port&)` | returned `−1` when the plus side was grounded, `+1` otherwise |
| `port_signed_amplitude` | `port_sign × port_amplitude` |
| `direction = 0 1 0` | an input key giving **d** *exactly* -- replaced 24 Sept by `current_direction`, a hint that need only be within 90 degrees (§2.1) |

The `s = −1` branch was reachable only through `ground = plus`, which
described the same physics as `ground = minus` with **d** reversed. Removing
the key made the branch unreachable, and an unreachable branch cannot be
tested against reality — it can only be "fixed" later in the wrong
direction.

**Floating internal ports, when they arrive, are not this.** They have *both*
sides unknown and no ground at all, which is a different mechanism (an extra
unknown plus a net-current-zero row), not a resurrection of `s = −1`.

The test that replaced the sign table asserts the property rather than the
value: an excitation must depend on amplitude and phase alone, never on the
port's kind. If internal and boundary ports ever diverge there, a sign
factor has crept back in.

---

## 3. Validation

All errors carry file and line. Two stages, because some checks need only
the text and some need the mesh.

### 3.1 Without the mesh (parse stage)

**Errors**
1. Unknown section, unknown key, duplicate key, duplicate section.
2. A required key missing; a value of the wrong kind (not a number, not a
   known keyword, a vector without 3 components).
3. Out-of-range values: `sigma < 0`, `eps_r <= 0`, `mu_r <= 0`, a frequency <= 0,
   a zero `current_direction`, or a `current_direction` on a boundary port.
4. `frequencies` or any `sweep`/`f_*` key given for `dc`; neither form given
   for `frequency`; or both forms given at once. Also `f_stop` <= `f_start`,
   or `f_points` < 2.
5. A key that does not belong to the port type (§1.2).
6. No `[Body …]` section, or no port.
7. The same surface used by two ports.
8. **No potential reference:** no voltage port and no internal port (every
   internal port is grounded). Φ would be determined only up to a constant
   and the matrix would be singular. The message suggests adding a
   `boundary_voltage` port with `voltage = 0`.

### 3.2 Against the mesh (binding stage)

**Errors**
1. A `volume` or `surface` entry names no physical group of the right
   dimension (3 for volumes, 2 for surfaces), or that group has no elements.
2. A physical volume present in the mesh that no `[Body …]` claims — or that
   two bodies claim.
3. A `boundary_*` port with any face that is not on the domain boundary; an
   `internal_*` port with any face that is not interior.
4. **An internal cut that is not planar.** Since **d** is now derived rather
   than supplied, this is a check that the cut's faces agree on a normal:
   taking `n̂₀` from the lowest-indexed cut face, every face must satisfy

   ```
   | |n̂_f · n̂₀| − 1 | < 1e-8
   ```

   The magnitude is what is tested: a face normal built from ascending node
   indices has an arbitrary sign, and the sign is not needed anyway — the
   plus side comes from the fourth-vertex test, not from `n̂`. The tolerance
   admits an angular deviation of `√(2e-8) ≈ 1.4e-4 rad ≈ 0.008°`, so it
   passes coordinate round-off on a planar cut and rejects any real
   curvature. **Recorded limitation:** a cut through a bent or twisted
   conductor is therefore rejected. Both near-term cuts — the ring's radial
   plane — are planar.
5. Terminals of two different ports sharing a node — a short circuit.
6. **DC only:** a port surface touching no conductor tet (σ > 0). It can
   carry no current.
7. **DC only:** a conductor piece (connected set of σ > 0 tets) that carries
   a port but has no potential reference.
8. **DC / reduced only — incomplete cut.** See §3.2.1.
9. **A conducting path with no potential reference.** Each connected
   conducting region is its own Φ problem at DC, so each needs at least one
   Dirichlet condition or Φ there is fixed only up to a constant. A voltage
   port supplies one; so does an internal port, which grounds its minus
   side; a boundary current port does not. **At least one, never exactly
   one** — two voltage ports on one conductor is a voltage-driven resistor,
   a well-posed Dirichlet problem whose current falls out as a result.

   *Corrected 25 Sept.* The first implementation required exactly one and
   rejected the voltage-driven case, and the test asserted that rejection,
   so the suite defended the bug. The cylinder passes either way, which is
   why it went unnoticed.

   A path with **no port at all** is floating: warned about rather than
   refused (an unconnected shield is legitimate), with its lowest-numbered
   node pinned to 0 V so the matrix stays non-singular. A floating conductor
   carries no DC current, so the pinned value is arbitrary. Full wave needs
   none of this: `(σ + jωε)` couples everything through the dielectric and
   one global reference suffices.

**Warnings**
- DC: faces of an internal cut lying outside conductors (they carry nothing
  at DC).
- Full wave: an internal cut is a delta gap. See §3.2.1.

#### 3.2.1 The internal-cut rim check

**The rule.** For an internal cut port, every **rim edge** — an edge
belonging to exactly one face of the cut — must lie on the boundary of Φ's
support. Because Φ's support depends on the regime (§1.2), so does the
verdict:

| regime | Φ's support | rim on its boundary | rim inside it |
|---|---|---|---|
| `dc`, or `frequency` + `reduced` | conductors (σ > 0) | clean — **the normal case** | **error**: the cut does not fully span its conductor |
| `frequency` + `full_wave` | whole domain | only if the cut reaches the outer box | **warning**: delta gap |

**The normal case needs no special handling.** A cut that is a complete
conductor cross-section has its rim on the conductor's lateral surface —
i.e. on the *edge* of where Φ lives, not inside it — so the jump closes
cleanly. No singularity, nothing to fix. This is what an internal port
normally is.

**Why the DC case is an error and not a warning.** If the tagged surface
stops partway through the conductor, or is drawn slightly off and misses the
far wall, part of the rim sits *inside* the conductor. Current then flows
around the uncut part and the port is partly shorted — but the matrix stays
non-singular, the solve succeeds, and it returns a plausible resistance that
is simply wrong. That is the silent-wrong-answer class this project keeps
finding, so it is worth failing loudly.

**Why the full-wave case cannot be an error.** Once Φ lives in the
dielectric too, a conductor cross-section *always* has its rim interior to
Φ's support: a path around the rim through the dielectric does not pick up
the jump, so Φ is discontinuous along the rim, the field there is singular,
and the computed gap capacitance drifts with mesh refinement. That is the
standard delta-gap idealization — tolerable and deliberate, hence a warning.
It would only be avoidable by extending the cut to the outer box, which
models something different (a source driving the whole box cross-section).

**Implementation.** The rim is found from the cut's face list alone (count
each edge's incident cut faces; rim = count 1). An edge is *on* the boundary
of Φ's support if any tet touching it lies outside that support, or if it
lies on the domain boundary. Boundary terminal ports need no such check:
their rim is on the domain boundary by construction.

### 3.3 What binding produces — the containers

**Why this step exists at all.** `Mesh` gives a complete edge/face/tet
topology but **no bodies and no surfaces**. Taking stock of what the reader
actually leaves behind (`include/aphi_solver/mesh.hpp`):

| present | missing |
|---|---|
| `nodes`, `tets` | **no body grouping** — only `tet_tags[t]`, one integer per tet; "the tets of body *k*" has to be built by scanning |
| `edges` (unique, `i<j`), `faces` (unique, sorted) | **no surface grouping** — `tagged_boundary_faces` is a *flat* list of node triples in file order, not grouped by tag and **not resolved to global face indices** |
| `tet_edges`, `tet_edge_signs`, `tet_faces` | **no name → tag lookup** — `physical_name(dim,tag)` only goes tag → name |
| `face_tets` (1 = boundary, 2 = interior, >2 detected) | **no node or edge set per surface** — needed for Φ slaving and the `n×A = 0` mask |
| `tet_tags`, `tagged_boundary_faces`, `physical_names` | **no material per tet**, no conductor connectivity for the DC checks |

So binding is exactly the step that turns the flat mesh into the containers
the rest of the solver wants:

```
bodies[]          name, tag, sigma, eps_r, mu_r, tets[]
surfaces[]        name, tag, faces[] (global indices), nodes[], edges[]
body_of_tet[]     per-tet index into bodies[]   (size = num_tets)
ports[]           faces[], derived d, plus_side_tet per face, node set
phi_support       which tets Phi lives on (regime-dependent, per §1.2)
warnings[]
```

plus `length_unit` applied — mesh coordinates scaled to metres, once, at
bind time.

**Two notes on the mesh reader, both settled 23 Sept:**

1. **`tagged_boundary_faces` is misnamed** and should become `tagged_faces`.
   The reader stores *any* tagged triangle without checking whether it lies
   on the boundary — which is precisely what makes internal ports possible.
   The name suggests a restriction the code does not have.
2. **Discarding triangle winding is fine; an earlier review was wrong to
   flag it.** The reader sorts each triangle's nodes, so the original
   orientation is lost — but it is not needed. **d** is built from the cut's
   own geometry and then the plus side follows from it: for each cut face,
   take its two adjacent tets from `face_tets` and pick the one whose fourth
   vertex lies on the `+d` side. No winding, no geometric tolerance.

**Deriving d (24 Sept).** Since a cut must be planar (§3.2, error 4), every
face shares one normal `n̂` up to sign, so only the sign is in question.

| case | rule |
|---|---|
| internal port **with** `current_direction` | `d = sign(v·n̂)·n̂` — `resolve_cut_direction`, rejecting `\|v̂·n̂\| < 0.1` |
| internal port **without** | `d = canonical_orientation(n̂)` — sign-fixed so the largest-magnitude component is positive, ties in x, y, z order. Arbitrary but reproducible |
| boundary port | `d` = the inward normal. Nothing to choose: a boundary face has exactly one adjacent tet, which is the inside |

**The dry run prints the resolved direction**, which is what closes the loop
for a user who cannot compute a normal — the tool proposes, they check:

```
ports
  P1  internal_current   loop_cut   412 faces
      d = (0, 1, 0)   [from current_direction = +y]
      plus side: 1196 tets    minus side (0 V): 1183 tets
```

and, with no hint given, the same line reads
`d = (0, 1, 0)   [derived -- sign of I and V is arbitrary; set current_direction to fix]`.
So the workflow never requires calculating anything: run it, read the
vector, and add a hint if it points the wrong way.

Assigning a side to tets that touch the cut only at an edge or a node
(flood-fill) belongs to the DOF map, not here.

---

## 4. Code layout

No new dependency: the grammar is small enough that a hand-written parser
gives better error messages than a general format library would.

| file | contents |
|---|---|
| `include/aphi_solver/problem.hpp` | `Problem` struct: mesh settings, analysis, boundary, bodies, ports; enums; `port_amplitude`, `phi_on_conductors_only`, `length_scale` |
| `include/aphi_solver/input_file.hpp`, `src/input_file.cpp` | `parse_input_file(path)` → `Problem`; `InputError` carrying file and line |
| `include/aphi_solver/problem_binding.hpp`, `src/problem_binding.cpp` | `bind_to_mesh(Problem, Mesh&)` → the §3.3 containers + warnings |
| `src/main.cpp` | `aphi_solver input.aphi`: read, bind, print a summary (bodies, ports, face counts, warnings). A dry run — no solve yet |
| `docs/INPUT_FILE.md` | user reference: grammar, every key, sign conventions, examples |
| `examples/cylinder_box.aphi` | the §1 example |
| `tools/cylinder_box.geo` → `meshes/cylinder_box.msh` | the first test geometry (§5) |
| `tests/test_input_file.cpp` | parser tests |
| `tests/test_problem_binding.cpp` | binding tests on the cylinder mesh and small hand-built meshes |

---

## 5. The cylinder test geometry

A wire of radius `a` along z, spanning the full height `ℓ` of a square box of
side `W`. Physical groups: volumes `wire`, `air`; surfaces `wire_bottom`,
`wire_top` (the wire's end disks, on the box). The rest of the box needs no
tag — it is the outer boundary, found from the mesh.

**Extruded**, not free-meshed: a 2-D mesh of the cross-section extruded in
layers, so the conductor is an exact prism.

**The closed form is exact here, not asymptotic.** With `n×A = 0` on all six
faces the exact 2-D square-coax field satisfies `B·n = 0` on the four walls
(that is what the conformal map arranges) *and* `B_z = 0` on top and bottom,
so those are the field's own boundary conditions and the solution is
genuinely z-invariant. The box is the return conductor. This is a stronger
test than the coax it replaces: there is no truncation error to converge
away, because nothing is being truncated.

### 5.1 Stage 1 — DC

```
R = ℓ / (σ · A_mesh)          A_mesh = the meshed cross-section area
L = (μ₀ℓ / 2π) [ ln(1.0787 · W / 2a) + 1/4 ]
```

`1.0787 = √2 / ∫₀¹(1−t⁴)^(−½) dt` is the conformal radius of a square's
interior at its centre, in units of the half-side: the Schwarz-Christoffel
map sends `t = 1` to a corner, which sits at the half-diagonal `(W/2)√2`, so
`r_c = 1.0787·(W/2) = 0.5394 W`. Then `L_ext = (μ₀/2π) ln(r_c/a)`. Sanity
check: `r_c` falls between the inscribed radius `0.5 W` and the circumscribed
`0.707 W`, nearer the inscribed one. The integral is half the lemniscate
constant, `1.3110288`.

> **Not to be confused with** the *logarithmic capacity* of a square,
> `0.5902 W`. That is the exterior problem — a charged square in free space —
> and it is the wrong number here. Easy mistake; it was made once during
> review.

**Accuracy is not the same for the two quantities, and the test should not
pretend otherwise.**

| | expected | why |
|---|---|---|
| `R` | **exact to round-off** | the meshed conductor is an exact prism, so the true potential `Φ = V·z/ℓ` is linear, lies in P1 ⊂ P2, and is reproduced exactly. Compare against `A_mesh`, the *faceted* cross-section area, not `πa²` |
| `L` | **within ~0.5 %**, converging under refinement | two systematic errors, both one-sided and both ignored by an `a`-based formula: the meshed wire is an n-gon whose conformal radius is slightly below `a` (O(1/n²)), and the closed form itself drops O((a/r_c)⁴) |

For the sizes below, `(a/r_c)⁴ ≈ 1.2e−3`, which entering through the
logarithm is ≈ 0.06 % in `L`; the faceting adds a comparable amount. Assert
0.5 % and require the value to move toward the closed form as the mesh is
refined — an unconverged pass is not a pass.

Sizes: `W = 2 mm`, `a = 0.2 mm`, `ℓ = 1 mm`, copper (σ = 5.8e7).

**Built and measured, 24 Sept** — `tools/cylinder_box.geo` → 
`meshes/cylinder_box.msh` (Gmsh 4.13.1, 614 KB, tracked):

| | measured |
|---|---|
| nodes / edges / faces / tets | 2912 / 19587 / 32716 / 16040 |
| Euler χ = V−E+F−T | **1** ✓ |
| total volume | **4.0000000000 mm³** = 2×2×1 exactly |
| wire (tag 1) | 7535 tets, 0.1242331416 mm³ |
| air (tag 2) | 8505 tets, 3.8757668584 mm³ |
| `A_mesh` | **0.1242331416 mm²** |
| `A_mesh / A_24-gon` | **1.0000000000** — the meshed wire *is* the exact polygon prism |
| `A_mesh / πa²` | 0.988616 (1.14 % low) |
| `wire_bottom` / `wire_top` | 122 triangles each, all on the domain boundary; 74 vertices each |
| outer boundary | 1272 faces |

**Targets for the DC solve:**

```
R = 0.1388 mOhm      exact to round-off, against A_mesh
                     (0.1372 mOhm is the ideal-circle value and is NOT
                      what this mesh should give -- a 1.14 % trap)
L = 0.3870 nH        within ~0.5 %
```

**The wire is a regular 24-gon, not a circle**, and that is what earns the
exactness. The linear potential `Φ = V·z/ℓ` is the exact solution only if it
also satisfies `∂Φ/∂n = 0` on the conductor's lateral surface — which needs
that surface vertical. Mesh a true cylinder and a boundary triangle
generally has its three vertices at three different heights *and* angles, so
its normal picks up a z-component and the linear field stops being the
discrete solution. A polygon's lateral faces are planar and vertical **by
construction**, so any triangulation of them has a horizontal normal.

That `A_mesh / A_polygon = 1` to full double precision is the load-bearing
fact, and `tests/test_gmsh_reader.cpp` pins it at `1e-12` relative rather
than with a loose tolerance. *Negative control:* regenerating the mesh as a
20-gon fails that check and the `πa²` ratio check, and nothing else.

**What stage 1 does *not* test.** At DC, Φ lives only in the wire — the air
has σ = 0, so it carries no Φ DOFs at all and the box walls carry none
either. The decision to leave **Φ free on the box therefore has no effect
here**; it is vacuous until ω > 0. A green stage 1 must not be read as
having validated it. That is what stage 2 is for.

### 5.2 Stage 2 — ω > 0, skin effect

The same geometry, same mesh, no new tags. The exact solution stays
z-invariant, and the closed form is the classical internal impedance of a
round wire:

```
Z(ω)/ℓ = Z_int(ω) + jω·L_ext
Z_int   = (γ / (2πaσ)) · I₀(γa) / I₁(γa),     γ = sqrt(jωμσ)
L_ext   = (μ₀/2π) ln(1.0787·W / 2a)
```

`L_ext` is frequency-independent: the `n×A = 0` walls are a lossless flux
barrier, so there is no wall skin effect. As ω → 0, `Z_int → R_DC +
jω μ₀/8π`, so **stage 1 is the ω → 0 limit of this same curve** — the `+1/4`
of §5.1 is exactly that `μ₀/8π`.

**Frequency range.** With the sizes above, δ = a at **109 kHz**
(`ω = 2/(a²μσ)`). A sweep of **1 kHz → 10 MHz** takes δ/a from 10.4 down to
0.105, covering the whole transition. λ at 10 MHz is 30 m against a 1 mm
structure, so no wave effect enters anywhere in that range: the sweep
isolates the `jωσ` term.

**What this buys:**
- It exercises **Φ free on the box** — the gap left by stage 1.
- `full_wave` and `reduced` must agree closely across this whole sweep,
  which cross-validates the "Φ everywhere" path against "Φ in conductors
  only" on a problem with a known answer.
- Skin effect is demanding: it needs correct edge orientation, a correct
  mass term, and enough radial resolution. A mesh convergence study in
  radius is part of the test, not optional.

**Reference values, practically.** `I₀`/`I₁` of complex argument are awkward
— at 10 MHz, `|γa| ≈ 13.5`, where the power series loses accuracy to
cancellation. Compute the references offline to high precision and hard-code
them in the test beside the formula and the method used, the way the tree's
27,732 count is pinned. Do not put a complex Bessel implementation in the
solver for a test's sake.

### 5.3 Stage 3 — displacement current (later)

Above the range of stage 2 the box becomes a transmission line, needing
`C' = 2πε / ln(1.0787·W / 2a)` and the box's own cutoff checked before any
result is trusted. Worthwhile, but a separate piece of work and not a
dependency of the input-file or DOF work.

### 5.4 Mesh generation

Use **Gmsh 4.13.1** — the version that produced `validation_box_v22.msh`,
`validation_box_v41.msh` and `coax_via.msh` — and record the version in a
header comment in the `.geo`. A different version can renumber entities and
size elements differently, which would silently change every count the tests
assert. (4.15.2 is also present on this machine; it is not the one to use.)

---

## 6. Implementation sequence

Each step builds and has its own tests.

1. **`Problem` struct and its conventions. Done 24 Sept** -- 63 checks in
   `tests/test_problem.cpp`: port classification; the §2 convention, phrased
   as a property (an excitation depends on amplitude and phase alone, never
   on the port's kind, so a reintroduced sign factor fails the test);
   complex amplitudes against hand values; the §1.2 Φ-support table; length
   scales; `resolve_cut_direction` and `canonical_orientation` (§3.3);
   and the cylinder problem assembled by hand, as a check that the
   struct can express the first real problem before a parser exists.
   *Negative control:* forcing the then-existing `port_sign` to a constant
   broke 7 checks across 4 test functions — after which the decision to
   derive `direction` retired that function entirely.
2. **Parser.** Tests: the §1 example parses to the expected struct; comments,
   case-insensitivity and defaults; each §3.1 error triggers with the right
   line number; a negative control per check (a valid file stays valid).
3. **Cylinder geometry. Done 24 Sept.** `.geo` (Gmsh 4.13.1, version in a header comment),
   mesh, track the `.msh` as a fixture; confirm names, volumes and that
   `wire_bottom`/`wire_top` are boundary faces. Record `A_mesh` and the
   facet count — §5.1's `R` is asserted against `A_mesh`, not `πa²`.
4. **Binding. Done 24 Sept.** Tests on the cylinder: names resolve, both volumes assigned,
   face counts, scaling to metres. Each §3.2 error on a small hand-built
   mesh: unassigned volume, internal face on a boundary port, shared terminal
   node, direction parallel to a cut, and — as its own case — the §3.2.1 rim
   check in all four regime/rim combinations, including a deliberately
   incomplete cut that must be rejected at DC and warned about at full wave.
5. **CLI dry run** and **`docs/INPUT_FILE.md`**.

The stage-1 and stage-2 solves themselves are not part of this plan — they
belong to Phase 04, after the DOF map. They are specified in §5 now so the
geometry is built once, with both stages in mind, rather than twice.

Then the DOF-map plan, starting from the bound cylinder problem.

---

## 7. Decisions for you

1. **File extension.** Suggest `.aphi` (still plain text).
2. **`sigma` required** in every `[Body …]` (suggested), or default 0?
3. **Amplitudes as peak** (suggested; matches `½ Re(V I*)`) or RMS?
4. **Units.** SI numbers only, with `length_unit` for the mesh (suggested for
   now), or also accept suffixes like `1 GHz` / `5.8e7 S/m`?
5. ~~Materials by `[material name]` + `volumes = …`, or one section per
   volume?~~ **Decided 23 Sept: `[Body name]`, one section per mesh volume,
   with named material keys.** See §1.2.
6. **Cylinder sizes** in §5.1.

*Status (23 Sept, Opus 5): all six review items applied. **The format is
settled** — decision 5 (`[Body …]`, one section per mesh volume, named
material keys) and the complex-amplitude grammar (magnitude plus
`phase_deg`) are both decided and written into §1.2. Decisions 1–4 and 6 —
extension, required `sigma`, peak amplitudes, units, cylinder sizes — still
carry my recommendations and are unanswered, but none of them blocks
starting: each changes a literal or a default, not the structure. Ready to
implement from §6 step 1.*
