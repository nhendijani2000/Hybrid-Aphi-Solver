# Plan: DOF setup + a Phase-04-minimal input file

> **SUPERSEDED (Sept 2026) — see
> `docs/PHASE_04_REVIEW.md`, and `docs/ROADMAP.md`
> Phase 03.5 / Phase 04, which are the source of truth.**
>
> Kept for the record. What carried forward: Section 2's mesh-tag gap and
> Section 5's gauge-reduction index-space gap were both real and are now
> Phase 03.5 steps 3 and 4. What changed: this draft missed the missing
> Whitney edge-orientation sign (a live bug, since fixed — commit
> `919cce5`) and understated the sparse-matrix problem as "template the
> gauge functions" when the actual blocker is that `APhiBlockSystem` is
> dense. Its open questions 1 and 2 are answered by one templated
> `Sparse<T>`; question 4 is answered by a DC straight-wire milestone
> before the coaxial via. Its Section 6 sequence is replaced.

*Draft for review -- Sept 2026. Not yet implemented, not yet committed to
`docs/`. Grounded in `docs/ROADMAP.md` Phase 04, `docs/FORMULATION.md`
Section 5, `docs/CONDITIONING.md`, and the actual current code (checked
just now, not from memory).*

## 1. Scope

Per your choice: **Phase-04-minimal**. The input file and DOF setup need to
carry exactly what `docs/ROADMAP.md` Phase 04 and `docs/FORMULATION.md`
Section 5 already say Phase 04 needs -- not the full ports/sources/
conductor-type machinery `docs/ROADMAP.md` Phase 06 explicitly owns:

- mesh reference
- material regions: (sigma, eps_r, mu_r) per region, tagged to the mesh
- gauge choice: Albanese-Rubinacci vs. Munteanu unsymmetric
  (`docs/TREE_COTREE_GAUGE.md`)
- frequency-scaling choice: Formulation 1/2/3
  (`docs/CONDITIONING.md`)
- operating frequency (a single value -- Phase 05 owns frequency sweeps)
- a simple uniform/lumped current-density source (Phase 04 step 4's own
  words: "a simple uniform current density is enough here")
- PEC boundary tagging (already needed by Phase 03's tree-cotree, just not
  yet wired to a real mesh's own tags)

Explicitly **not** in this pass (reserved, not designed in detail):
stranded vs. massive conductor source models, a real port abstraction,
S-parameter/impedance extraction, lossy-dielectric frequency models beyond
a flat sigma. Those are Phase 06. The schema below reserves empty slots for
them so the file format doesn't change shape later, but their content
isn't decided now.

## 2. A gap found while grounding this plan: mesh ingestion discards all tags today

Checked `src/gmsh_reader.cpp` directly (not from memory) before writing
this. Right now, for every element in the `$Elements` section:

```cpp
elem_line >> elm_number >> elm_type >> num_tags;
for (long long tg = 0; tg < num_tags; ++tg) {
    long long discard_tag;
    elem_line >> discard_tag;
}
```

Every tag is read and thrown away -- including a **tet's own** physical-group
tag, not just boundary-surface tags. And non-tet elements (the triangles
Gmsh uses to tag boundary surfaces -- PEC walls, ports, the outer
truncation) are skipped entirely; `Mesh` has no field to hold them at all.

This means region/boundary tagging isn't "not fully designed yet" -- it's
**completely absent** from mesh ingestion today. Since DOF setup, material
lookup, and PEC handling all depend on knowing which mesh entities belong to
which named region, this has to be step 0, before the input file can mean
anything.

**Fix, scoped small:**
- Add `std::vector<int> tet_tags;` to `Mesh` (one physical-group tag per
  tet, parallel to `Mesh::tets`).
- Add `std::vector<std::pair<std::array<int,3>, int>> tagged_boundary_faces;`
  (or similar) capturing Gmsh's 3-node triangle elements (elm-type 2) and
  their first tag, remapped to this codebase's 0-based node indices --
  these are the PEC/port/outer-boundary surface tags.
- `read_gmsh_msh` keeps the tet's first tag (Gmsh's convention: tag 1 is the
  physical-group id) instead of discarding it, and now also parses
  elm-type 2 (3-node triangle) elements into `tagged_boundary_faces` instead
  of skipping every non-tet element.
- No interpretation happens here -- `Mesh` just carries integers. What tag
  `7` *means* (copper, PEC, port-1) is entirely the input file's job. This
  keeps `gmsh_reader.cpp` dependency-free and format-agnostic, matching its
  existing scope.

This is a small, independently testable change (a new
`tests/test_gmsh_reader.cpp` case: a hand-written `.msh` with tagged tets
and boundary triangles, check `tet_tags`/`tagged_boundary_faces` come back
right) -- not a big detour, but a real prerequisite.

## 3. Input file schema (JSON)

Worked example using `docs/FORMULATION.md` Section 7.1's own EDA target
problem (coaxial via/power-plane analog) -- not a synthetic placeholder,
the actual problem Phase 04-05 need to reproduce `Z_0`, `L'`, `C'`, `R_DC`
for:

```json
{
  "mesh": "meshes/coax_via.msh",

  "regime": "full",
  "// regime": "'full' keeps eps everywhere (Section 1); 'low_frequency_reduced' forces eps=0 and confines Phi to regions with sigma>0 (Section 3). Not a separate code path -- a material-property and domain-restriction input, per FORMULATION.md Section 3.",

  "regions": {
    "1": { "name": "via_conductor",   "sigma": 5.8e7,  "eps_r": 1.0, "mu_r": 1.0 },
    "2": { "name": "return_conductor","sigma": 5.8e7,  "eps_r": 1.0, "mu_r": 1.0 },
    "3": { "name": "fr4_dielectric",  "sigma": 0.0,    "eps_r": 4.3, "mu_r": 1.0,
           "loss_tangent": 0.02 }
  },
  "// regions": "keys are the mesh's own Gmsh physical-group tags (Section 2 above). loss_tangent, when present, derives a frequency-dependent sigma = omega*eps_r*eps_0*tan(delta) at assembly time (Section 7.1) instead of a flat sigma.",

  "boundary_conditions": {
    "10": { "type": "pec" },
    "11": { "type": "pec" }
  },
  "// boundary_conditions": "tag -> BC type. Only 'pec' is meaningful in this Phase-04-minimal pass (needed by tree-cotree, Phase 03). 'abc' (Phase 07) and 'port' (Phase 06) are valid future values but not implemented against yet.",

  "gauge": "albanese_rubinacci",
  "// gauge": "'albanese_rubinacci' | 'munteanu_unsymmetric', per docs/TREE_COTREE_GAUGE.md",

  "frequency_scaling": "natural",
  "// frequency_scaling": "'natural' (Formulation 3) | 'symmetric_row_scaled' (Formulation 1) | 'scaled_scalar_potential' (Formulation 2), per docs/CONDITIONING.md. Matches ConditioningStrategy in conditioning.hpp.",

  "frequency_hz": 1.0e9,

  "source": {
    "type": "uniform_current_density",
    "region": 1,
    "J": [0.0, 0.0, 1.0e6]
  },
  "// source": "Phase 04 step 4's own scope: 'a simple uniform current density is enough here'. Real stranded/massive conductor and port excitation models are Phase 06.",

  "ports": [],
  "excitations": [],
  "// reserved": "empty on purpose -- Phase 06 defines these. Keeping the keys present now means Phase 06 adds entries, not a format migration."
}
```

Notes:
- `//`-prefixed keys are comments-as-data (plain JSON has no comment
  syntax); the actual parser ignores any key starting with `//`. Cheap and
  keeps the schema self-documenting without a second spec document to keep
  in sync. Easy to drop later if you'd rather use a JSON5/JSONC library
  instead -- flagging as a style choice, not a hard commitment.
- `regions`/`boundary_conditions` keyed by tag-as-string because JSON
  object keys must be strings; parsed back to `int` on load.
- This format needs a JSON library. Recommend a small header-only one
  (e.g. `nlohmann::json`, vendored via CMake `FetchContent` into
  `third_party/`, which `.gitignore` already has a slot for -- "Third-party
  fetched sources... not committed"). This is ordinary software tooling,
  not EM-specific IP, so it doesn't touch the project's literature-only
  constraint on the *solver* content.

## 4. DOF numbering plan

Grounded directly in `docs/FORMULATION.md` Section 5.2's own formulas:

```
N_A   = N_edges                                    (regime: full)
      = N_edges(Omega_c union boundary edges)       (regime: low_frequency_reduced)
N_Phi = N_vertices + N_edges                        (full -- P2 has one DOF
                                                      per vertex AND per edge
                                                      midpoint)
      = N_vertices(Omega_c) + N_edges(Omega_c)       (low_frequency_reduced)
```

`Omega_c` (conductive region) is **derived**, not separately tagged: a
region counts as conductive iff its `sigma > 0` in the input file. This
avoids a second, redundant classification the input file would otherwise
need to keep consistent with the material table by hand.

Proposed concrete numbering, contiguous blocks (matches how
`conditioning.hpp`'s `APhiBlockSystem` already keeps `K_AA`/`K_APhi`/
`K_PhiA`/`K_PhiPhi` as **separate** matrices rather than one flattened
global matrix -- this plan follows that existing convention rather than
introducing a second one):

- A-DOFs: one per mesh edge in `Omega_A` (all edges, or only
  `Omega_c`-touching edges under the reduced regime), in `Mesh::edges`
  order, filtered.
- Phi-DOFs: one per mesh vertex in `Omega_Phi`, followed by one per mesh
  edge in `Omega_Phi` (matching the P2 element's vertex-then-edge-midpoint
  node convention already fixed in `basis_functions.hpp`'s `local_node`
  0..9 ordering), both filtered to `Omega_c` under the reduced regime.
- Each DOF set keeps a `global_index -> local_dof_index` map (`-1` for "not
  a DOF", e.g. a vertex outside `Omega_c` under the reduced regime, or a
  vertex on a PEC-only boundary) -- the same pattern `GaugeVariant` and
  `EssentialIncidenceMatrix` already use (`cotree_local_index`,
  `free_group_index`), so this is a natural extension of an existing idiom,
  not a new one.

**PEC handling**: `FORMULATION.md` Section 5.4 says PEC-boundary tangential-
edge DOFs are removed from the unknown set. This is *already* how
`build_tree_cotree`'s `is_pec` argument works (it grounds the spanning tree
through PEC-tagged nodes). Once Section 2's boundary-face tags exist, `is_pec`
for `build_tree_cotree` becomes derivable directly: a node is PEC iff it's a
vertex of a boundary face tagged `"type": "pec"` in the input file. No new
PEC machinery needed -- just wiring the existing parameter to the new tag
data instead of a hand-built `std::vector<bool>` (as `compare_gauges.cpp`'s
`--pec-nodes` currently does for test purposes).

## 5. A real architectural gap this plan surfaces: gauge reduction assumes a bare edge-indexed matrix

Checked `include/aphi_solver/gauge_variants.hpp` directly. Both
`build_albanese_rubinacci_gauge(const SparseMatrix& M, const
TreeCotreeResult& tc)` and `build_munteanu_unsymmetric_gauge(...)` assume
`M`'s row/column index space **is** the mesh's global edge index space,
1:1 -- there's no concept of "a matrix where only some rows/columns
correspond to gauge-reducible A-DOFs and the rest (Phi rows/columns) pass
through untouched." That assumption is correct for `test_gauge_variants.cpp`
and `compare_gauges` (both only ever hand it `M = C^T*C`, size
`N_edges x N_edges`), but it doesn't fit the real Phase 04 system, where
gauge reduction has to act on just the A-DOF rows/columns of a larger
coupled system while leaving Phi rows/columns alone.

Two ways to close this gap:

**(a) Generalize the existing functions** to take an explicit "which global
indices are A-DOFs" set, so they can be handed either a bare `N_edges x
N_edges` matrix (today's callers, unchanged) or a slice of a larger block
system.

**(b) Keep `gauge_variants.hpp` exactly as-is**, and add new, separate
assembly-side code (a `dof_map.hpp`/`assembly.hpp`, wherever Phase 04's
weak-form code lands) that builds the full `APhiBlockSystem`, then applies
the *same reduction rule* (drop tree-edge columns from `K_AA`/`K_PhiA`;
project/select `K_AA`/`K_APhi`'s rows the same way `build_*_gauge` already
does) directly against the block structure, calling into
`compute_essential_incidence_matrix`/`TreeCotreeResult` for the shared
tree/cotree bookkeeping rather than re-deriving it.

I lean toward **(a)**: the actual reduction math (test space W, trial space
V, per `docs/TREE_COTREE_GAUGE.md`) is identical either way, and (b) risks
re-implementing that same math a second time against a different matrix
shape -- exactly the kind of duplication `docs/ENGINEERING_STANDARDS.md`
already has us avoiding. But this needs your sign-off since it means
changing an existing, already-tested (19/19), already-shipped-in-a-CLI-tool
function signature, not just adding new code.

**A second, separate gap in the same area**: `gauge_variants.hpp`/`.cpp`
operate on the real-valued `SparseMatrix` (`incidence.hpp`), but Phase 04's
actual assembled blocks are `ComplexMatrix` (`conditioning.hpp`) -- the
frequency-domain coefficients are complex. Gauge reduction is pure index
selection/combination, so it's valid over `ComplexMatrix` exactly the same
way, but the current code physically can't be called on one without either
a templated rewrite or a second, complex-typed copy of the same three
functions. Also needs your call before implementation starts.

## 6. Proposed implementation sequence

Each step independently buildable and testable, in this order:

1. **Mesh tag ingestion** (Section 2 above): `Mesh::tet_tags`,
   `Mesh::tagged_boundary_faces`, `read_gmsh_msh` keeps tet tags and parses
   triangle elements. New `tests/test_gmsh_reader.cpp` cases.
2. **Input file schema + parser**: a `Problem` struct (regions, boundary
   conditions, gauge/scaling/frequency/source choices) and a
   `read_problem_json(path, mesh)` that also cross-checks every tag the
   input file mentions actually exists in the loaded mesh (a real, cheap
   correctness check worth having from day one -- a typo'd region tag
   should fail loudly, not silently assemble a region with default/zero
   material properties). Vendor the JSON library first.
3. **Resolve Section 5's architectural question** (gauge reduction over a
   block system, real vs. complex) -- a short spike/decision, not full
   implementation, before DOF setup code depends on a signature that might
   change under it.
4. **DOF setup**: given `Mesh` + `Problem`, compute `Omega_c` (derived from
   `sigma > 0`), build the A-DOF and Phi-DOF index maps from Section 4,
   derive `is_pec` from `boundary_conditions`, call `build_tree_cotree`.
   Testable independently of assembly: given a known small mesh + a known
   `Problem`, assert the exact DOF counts and index maps by hand
   (mirrors how `test_tree_cotree.cpp` already works).
5. **Only then**: element-level weak-form matrices (Phase 04 step 1) and
   global assembly (step 2), now with real DOF maps and real material
   coefficients to assemble against, plus the uniform-current-density RHS
   (step 4).

Steps 1-2 are the direct answer to "we need a plan for a practical input
file before DOF setup"; step 4 is DOF setup itself; step 3 is the one piece
of real design risk this plan surfaced that's worth resolving deliberately
rather than discovering mid-implementation.

## 7. Open questions for you

1. Section 5: generalize `gauge_variants.hpp`'s existing functions (option
   a), or add separate block-aware code that reuses their derivation
   without changing their signature (option b)?
2. Real vs. complex gauge reduction: template the existing functions, or
   write a second complex-typed copy?
3. JSON library choice: `nlohmann::json` (most common, single header) is
   the default assumption above -- any preference, or is that fine?
4. Confirm the coaxial-via problem (`docs/FORMULATION.md` Section 7.1) as
   the actual first mesh+input-file pair to build, since a real target
   (rather than another synthetic cube) is what finally exercises this
   end-to-end.
