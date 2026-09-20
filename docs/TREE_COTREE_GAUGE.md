# Tree-cotree gauge variants: derivation

*Written Sept 2026, in response to a direct question about where the
Method D (Munteanu "unsymmetric") reduced-matrix formula actually comes
from. `docs/REFERENCES.md` records what was verified directly against I.
Munteanu, "Tree-cotree condensation properties," and discloses that her own
equations (19)-(20), as extracted from the PDF, showed what looks like a
symbol collision ("F" used for both "number of mesh faces" and "the
essential incidence matrix" in different sections) -- too risky to copy
verbatim. This document is the derivation that was done instead, from her
stated general framework and explicit basis choices, worked through fully
and independently rather than summarized. It matches `src/gauge_variants.cpp`
exactly; if the two ever disagree, the code's behavior is the bug.*

## 1. Setup and notation

Let `a ∈ R^E` be the vector of edge-based (Whitney/Nédélec, first-order)
degrees of freedom for the magnetic vector potential **A**, one per mesh
edge (`E` = `mesh.num_edges()`). The curl-curl system this project will
eventually assemble in Phase 04 has the form

```
M a = j,      M = C^T ν C
```

where `C` is the discrete curl matrix (`build_curl_matrix`, edges → faces)
and `ν` is the (symmetric, positive semi-definite) material reluctivity
matrix. `M` is symmetric and positive *semi*-definite, and singular: its
null space is exactly the set of discrete gradient fields (Bossavit 1998,
already cited in this project for the `C·G = 0` identity that Phase 02's
tests verify) -- physically, adding `grad(φ)` to **A** for any scalar `φ`
doesn't change **B** = curl(**A**), so it doesn't change what `M a` measures
either. Gauging means adding enough independent linear constraints to
remove exactly that null space, without removing any physically meaningful
solution.

Phase 03 step 1 (`build_tree_cotree`) already does the topological half of
this: it builds a spanning tree/forest over the mesh's node graph (with
each PEC body collapsed to a single reference/ground group, per
`docs/ROADMAP.md` Phase 03 step 4), and partitions the `E` edges into
**tree edges** (the ones used by that spanning tree/forest) and **cotree
edges** (everything else). Relabel/reorder `a`'s components accordingly:

```
a = [a_c; a_t],   a_c ∈ R^(num_cotree),   a_t ∈ R^(num_tree)
```

This is purely a re-indexing (no change of basis), matching
`GaugeVariant::cotree_local_index` in the code.

Define the **reduced gradient matrix** `G` (size `E × num_free`, `num_free`
= number of non-reference DOF groups = `num_tree` exactly, per the
tree-edge-count invariant `tree_edge_count == num_groups -
num_reference_groups` already checked in `tests/test_tree_cotree.cpp`): row
`e = (i,j)` (canonical `i<j`) has `-1` at the column for `group(i)` *if that
group is free* (0 if it's a reference/PEC body, since its potential is
fixed and not a free unknown), and `+1` at the column for `group(j)` under
the same rule. This is exactly `build_gradient_matrix`'s own convention
(`incidence.hpp`), just with reference-body columns dropped. Partition its
rows the same cotree/tree way:

```
G = [G_c; G_t],   G_c: num_cotree x num_free,   G_t: num_tree x num_free (square, num_tree = num_free)
```

**Claim: `G_t` is invertible.** This is the standard fact that the
(ground-collapsed) incidence matrix of a spanning tree is nonsingular: every
free group is reached by *exactly one* tree path back to some reference
root (that's what a tree/forest means), so `G_t`'s rows/columns can be
ordered (by BFS discovery order -- exactly `TreeCotreeResult::parent_group`
/ `discovering_edge`, already computed) into a lower-triangular structure
with `±1` on the diagonal after a suitable permutation. This is *why*
`compute_essential_incidence_matrix` can solve `G_t x = b` by a single
top-down tree walk (`solve_tree_system` in `gauge_variants.cpp`) instead of
a general matrix inversion -- it isn't an optimization on top of a generic
solve, it's a direct reading of this triangular structure.

## 2. The "essential incidence matrix" F, from first principles

Define `F = G_c * G_t^{-1}` (size `num_cotree x num_free`) -- Munteanu's own
name and definition, Sec. III.A, independently reconstructed below.

**Claim: the condition `G^T a = 0` is equivalent to `a_t = -F^T a_c`.**

Write `G^T = [G_c^T | G_t^T]` (columns of `G^T` are edges, so this is the
same cotree/tree split applied to columns instead of rows). Then

```
G^T a = G_c^T a_c + G_t^T a_t
```

Setting this to zero and solving for `a_t` (using that `G_t`, hence `G_t^T`,
is invertible):

```
G_t^T a_t = -G_c^T a_c
a_t = -(G_t^T)^{-1} G_c^T a_c = -(G_c G_t^{-1})^T a_c = -F^T a_c
```

This derives the paper's stated relation (her eq. 14, "`G_0^T a = 0 ⇒ a_t =
-F^T a_c`", which *was* extracted cleanly and unambiguously, unlike eqs.
19-20) directly from the definition of `F`, rather than needing to trust
that extraction at all.

## 3. Why `{a : G^T a = 0}` is a valid gauge choice

Call `S_s = {a ∈ R^E : G^T a = 0}` ("solenoidal", Munteanu's term, Sec. V).
Two things need checking for this to be a legitimate gauge-fixing: `S_s`
must have the right dimension, and it must meet every physically-equivalent
"coset" `a + Range(G)` (the set of vector potentials that differ only by a
gradient, hence give the same **B** field) in *exactly one point*.

**Dimension.** `G` has full column rank `num_free` (immediate from `G_t`
being an invertible `num_free x num_free` block of it), so `G^T` has rank
`num_free` too, and `S_s = ker(G^T)` has dimension `E - num_free =
num_cotree`. The parametrization `a = [a_c; -F^T a_c]` from Section 2, as
`a_c` ranges over all of `R^num_cotree`, traces out exactly this
`num_cotree`-dimensional space.

**Trivial intersection with `Range(G)`.** Suppose `Gφ ∈ S_s` for some free-
group potential vector `φ`, i.e. `G^T(Gφ) = 0`. Dot both sides with `φ`:

```
0 = φ^T G^T G φ = (Gφ) · (Gφ) = |Gφ|^2   ⟹   Gφ = 0
```

So `S_s ∩ Range(G) = {0}`. Since `dim(S_s) + dim(Range(G)) = num_cotree +
num_free = E`, this gives `R^E = S_s ⊕ Range(G)` (a direct sum) -- every
coset `a + Range(G)` meets `S_s` in exactly one point. That point is what
Method D solves for.

(The same argument, run on `{a : a_t = 0}` instead of `S_s`, shows *why*
Method A's cruder gauge choice is equally valid: if `Gφ` has `(Gφ)_t = 0`,
i.e. `G_t φ = 0`, then `G_t` invertible forces `φ = 0`, hence `Gφ = 0` --
the same trivial-intersection argument, just against a different
`num_free`-codimension-`num_cotree` subspace. Both methods are legitimate
gauge choices for the same underlying reason: `G_t`'s invertibility.)

## 4. The reduced system

Both variants share the same **test space**: Munteanu's stated choice
(Sec. V) is the plain cotree-row selector `W = [I; 0]` (an `E x num_cotree`
matrix whose columns are the standard basis vectors for the cotree
positions) -- i.e., "test" the full system `M a = j` by keeping only its
cotree rows, discarding the tree rows entirely. This is a genuine modeling
choice attributable to her (not re-derivable from anything more
fundamental) and is what both `build_albanese_rubinacci_gauge` and
`build_munteanu_unsymmetric_gauge` do identically via
`select_cotree_entries`.

Where the two methods differ is the **trial space** -- how `a_t` is
expressed in terms of `a_c` before substituting into those cotree rows.

Write `M` in cotree/tree block form (symmetric, so `M_tc = M_ct^T`):

```
M = [ M_cc   M_ct ]
    [ M_tc   M_tt ]
```

**Method A (Albanese-Rubinacci):** trial space is literally `a_t = 0`.
Substituting into the cotree rows of `M a = j`:

```
(M a)_cotree = M_cc a_c + M_ct·0 = M_cc a_c  =  j_c
```

so the reduced system is simply `M_cc a_c = j_c` -- exactly `M`'s own
cotree-cotree principal submatrix, which is why it preserves `M`'s sparsity
exactly (it's a genuine submatrix, not a product of anything).

**Method D (Munteanu unsymmetric):** trial space is `a = L^T a_c` where
`L = [I | -F]`, i.e. `a_t = -F^T a_c` (Section 2's solenoidal condition).
Substituting into the cotree rows of `M a = j`:

```
(M a)_cotree = M_cc a_c + M_ct a_t = M_cc a_c + M_ct(-F^T a_c) = (M_cc - M_ct F^T) a_c  =  j_c
```

so the reduced system is `(M_cc - M_ct F^T) a_c = j_c` -- this is exactly
`build_munteanu_unsymmetric_gauge`'s `product` matrix (computed there as
`(select cotree rows of M) * L^T`, the same expression written out in dense
matrix form rather than block-algebra form).

## 5. Why Method D is "unsymmetric" (and Method E would be "symmetric")

`M_cc` is symmetric (a principal submatrix of a symmetric matrix always
is). But `M_ct F^T` has no reason to be: `(M_ct F^T)^T = F M_ct^T = F M_tc`,
and there's nothing forcing `F M_tc = M_ct F^T` in general, since `F` comes
purely from the mesh's tree topology (Section 2) and has no relationship to
`M`'s specific entries. So `M_cc - M_ct F^T` is generically **not**
symmetric, even though the original `M` is -- exactly why Munteanu calls
this variant "unsymmetric", and a useful independent check that this
derivation lines up with her own naming rather than having gone wrong
somewhere.

For contrast (not implemented in this project, since it isn't one of the
two variants asked for, but worth noting for why the name split makes
sense): Munteanu's "symmetric" variant (her Method E) instead uses `W = V =
L^T` for *both* test and trial spaces, giving a reduced matrix `L M L^T`.
That form is automatically symmetric for *any* `L`, since `(L M L^T)^T = L
M^T L^T = L M L^T` whenever `M = M^T` -- symmetric by construction,
regardless of what `L` is, which is exactly the distinguishing feature the
name refers to.

## 6. Measured behaviour on real meshes (Sept 2026)

The condition-number ordering in Section 5's contrast is Munteanu's, checked
here only as a direction of comparison. Once `tools/compare_gauges` could
run on real meshes, three sizes were measured with `M = C^T C` (vacuum,
`nu = 1`):

| mesh | cotree DOFs | kappa_A | kappa_D | kappa_D/kappa_A | nnz_D/nnz_A |
|---|---|---|---|---|---|
| cube_4 | 480 | 629.8 | 152.8 | 0.243 | 2.97 |
| cube_6 | 1512 | 1876.7 | 525.6 | 0.280 | 4.37 |
| cube_9 | 4860 | 5731.0 | 2041.7 | 0.356 | 6.49 |

**The ordering reproduces Munteanu's (`kappa_D < kappa_A`), and the numbers
are converged**, not artifacts of a truncated estimate: re-running
`estimate_condition_number` with its outer iteration cap raised from 50 to
5000 moves `kappa_A`/`kappa_D` by less than one part in 10^5 at both cube_4
and cube_6, and the estimator itself is independently validated against
diagonal matrices with analytically known condition numbers
(`tests/test_gauge_variants.cpp`).

**Why Method D is better conditioned**, in terms of Section 4's algebra:
the two methods project onto different subspaces. Method A restricts to
`{a : a_t = 0}` -- an arbitrary *coordinate* subspace, since forcing the
potential to vanish on a spanning tree bears no relation to the operator.
Method D's trial space `range(L^T) = {a : a_t = -F^T a_c}` is exactly
`{a : G^T a = 0}` (Section 3), the discrete gauge-consistent subspace. D
eliminates along a structurally meaningful direction; A slices arbitrarily.

**But a lower condition number does not make D the better choice here, and
both trends run against it.** Across those three meshes D's conditioning
advantage *shrinks* (0.243 -> 0.280 -> 0.356) while its fill-in *grows*
(2.97x -> 4.37x -> 6.49x). On top of that, D's reduced matrix is
non-symmetric by construction (Section 5), which forecloses COCG/COCR and
symmetric-indefinite factorization (`docs/CONDITIONING.md`'s decision
matrix). Paying 6.5x the nonzeros and a more expensive solver class for a
2.8x conditioning edge that is eroding with mesh size is a poor trade, and
it is trending worse -- which is part of why Phase 04 targets
Albanese-Rubinacci first (`docs/ROADMAP.md`).

**The caveat that outweighs all of the above:** `M = C^T C` is a vacuum
curl-curl stand-in, *not* the A-Phi system. The assembled Phase 04 matrix
carries the `(j*omega*sigma - omega^2*eps)` mass term and the A-Phi
coupling blocks, and that mass term regularizes precisely the curl-curl
null space this whole comparison is about. So every number in the table is
a statement about `C^T C`, not about the matrix this solver will actually
factor. Re-run the comparison in Phase 05 against real assembled systems
before treating any of it as a gauge recommendation.

## 7. What this derivation does and doesn't rely on the paper for

Attributed to Munteanu's paper, and taken as a design choice rather than
re-derived: the specific *names* Albanese-Rubinacci/unsymmetric for these
two constructions, the choice of test space `W = [I;0]` for both, and the
reported condition-number ordering used only as a direction-of-comparison
sanity check (not as ground truth for these tiny test meshes).

Derived independently in this document, and not dependent on trusting any
specific equation's OCR/extraction fidelity: `G_t`'s invertibility, the
equivalence of `G^T a = 0` and `a_t = -F^T a_c`, the transversality argument
for both gauge choices, the block-algebra reduction to `M_cc a_c = j_c`
(Method A) and `(M_cc - M_ct F^T) a_c = j_c` (Method D), and the
symmetric-vs-unsymmetric contrast in Section 5.

`tests/test_gauge_variants.cpp` checks the *result* of this derivation two
independent ways: Method A's reduced matrix is checked directly against
`M`'s dense cotree submatrix, and both methods are checked to recover a
vector potential whose curl matches `curl(z)` for a manufactured, self-
consistent source `j = M*z` -- a check that would fail if the block algebra
above had a sign or transpose error, regardless of whether the general
approach was sound.
