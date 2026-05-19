# Lead: Recipe 1 (q-projection of b5·R̃²) against actual xAct decomposition

**Date**: 2026-04-27
**Author**: Phase 2.2 investigator
**Scope**: Cross-check Round 3 Agent H's "Recipe 1 PASS" verdict against
TIDAL's actual production-pipeline JSON output and Wolfram source.

**Pipeline files inspected**
- `/workspaces/torsion-gertsenshtein/examples/torsion_gertsenshtein/theory.toml`
- `/workspaces/torsion-gertsenshtein/examples/data/torsion_gertsenshtein.json`
  (119 045 bytes; 38 fields, 38 equations, 78 hamiltonian terms)
- `/workspaces/torsion-gertsenshtein/tidal/wolfram/ComponentDecompose.wl`
  (2 204 lines)
- `/workspaces/torsion-gertsenshtein/tidal/wolfram/ExportJSON.wl`
  (1 692 lines)
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/scripts/recipe1_preflight_q_projection.py`
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/results/recipe1_q_kinetic_structure.json`

---

## Section 1 · Actual t-field inventory from JSON

The JSON contains **24 `t_*` torsion fields**, **10 `h_*` graviton
fields**, and **4 `a_*` photon fields** (38 dynamical fields total).

Per-equation `lhs.order.time` and `lhs.kinetic_coefficient_symbolic` for
the 24 torsion components (extracted from `equations[]`):

| field | time_order | kinetic_coefficient_symbolic |
| ----- | ---------- | ---------------------------- |
| t_0  …t_5 (except t_6)  | 0 | (algebraic, no kinetic coefficient)|
| t_7  …t_12 | 0 | (algebraic) |
| t_14 …t_19 | 0 | (algebraic) |
| t_21 …t_23 | 0 | (algebraic) |
| **t_6**   | **2** | **(-25\*b5)/2** |
| **t_13**  | **2** | **(-25\*b5)/2** |
| **t_20**  | **2** | **(-25\*b5)/2** |

Of the 24 t-fields, **only three** (t_6, t_13, t_20) propagate.  Every
propagating torsion component has `time_derivative_order = 2` (not 4)
with kinetic coefficient `(-25/2)·b5` — i.e. the b5·R̃² term promotes
these three torsion components to **standard second-order kinetic**
status, NOT to fourth-order Pais–Uhlenbeck status.

The remaining 21 torsion components are algebraic constraints
(`time_order = 0`, `lhs.expression = t_n`).

For contrast, the 10 h-fields:

| field   | lhs expression | kinetic coefficient |
| ------- | -------------- | --------------------|
| h_0…h_3 | algebraic      | n/a                 |
| **h_4** | `d4_t(h_4)`    | `2*b5`              |
| h_5     | `d2_t(h_5)`    | `-kappa^(-2)`       |
| h_6     | `d2_t(h_6)`    | `-kappa^(-2)` (~)   |
| **h_7** | `d4_t(h_7)`    | `2*b5`              |
| h_8     | `d2_t(h_8)`    | `-kappa^(-2)` (~)   |
| **h_9** | `d4_t(h_9)`    | `2*b5`              |

The genuine fourth-order Pais–Uhlenbeck kinetic terms live exclusively
in the metric subspace `{h_4, h_7, h_9}`, exactly as the audit's
sectoral classification asserted.  **No torsion field has `time_order ≥ 4`.**

(JSON path for each entry: `equations[i].lhs.order.time` and
`equations[i].lhs.kinetic_coefficient_symbolic`; tabulation produced by
iterating `equations[]` in the JSON.)

---

## Section 2 · Comparison to Agent H's claimed schema

`results/recipe1_q_kinetic_structure.json` records:

```json
{ "max_derivative_order_in_R_tilde": 1,
  "max_derivative_order_in_R_tilde_squared": 2,
  "q_free_parameter_count": 16,
  "DT_x_DT_term_count_existing_enumeration": 16 }
```

Agent H's claim is, in TIDAL terms: t-fields obtained by projecting
b5·R̃² onto the q-irreducible sector should be standard-kinetic
(`time_order = 2`), NOT Pais–Uhlenbeck (`time_order = 4`).

**Production-level result is consistent with that claim.** Three
t-fields show up with `time_order = 2` and a b5-proportional kinetic
coefficient.  Zero t-fields show up with `time_order = 4`.

However, the **counts do not match cleanly**:

- Agent H expects 16 q-irreducible degrees of freedom (from his sympy
  constraint solve — verified in his preflight).
- TIDAL's pipeline produces **3** propagating t-components, not 16.
- The remaining 21 t-components are algebraic constraints, not free
  q-irreducible DOF.

The reconciliation is:

(a) `ComponentDecompose.wl` does **not** perform an HMMN-style
    irreducible decomposition of `T_{abc}` into (trace, axial,
    tensor-q).  It does direct numerical-component projection of
    `T_{abc}` (rank-3, antisymmetric in last two indices in 4-D ⇒
    4·6 = 24 distinct components).

(b) The `b5·R̃²` term, in flat-background plane-wave reduction along
    `z`, sources kinetic (`d²_t`) terms for only **three** of those
    24 components (those that survive the simultaneous projection of
    R̃² onto t-only and onto the z-propagation reduction).  The other
    21 components remain algebraic.

(c) Agent H's "16 q-DOF" count is the off-shell reducible-rep count
    in 4-D *before* plane-wave reduction.  TIDAL's "3 propagating
    t-components" is the count *after* projecting onto z-axis
    propagation `(∂_x = ∂_y = 0)` defined in `theory.toml:[reduction]`
    (`type = "plane_wave"`, `propagation_axis = "z"`).

These are not contradictory; they live at different stages of the
pipeline.  But Agent H's schema script never reduced to plane-wave,
so direct numerical comparison of "16 vs 3" is misleading.

**However, an important caveat surfaces here: Agent H's analysis
is performed for the parity-odd Pontryagin density**
`(½) ε^{μνρσ} R̃_{μνρσ} R̃²` — see his
`recipe1_explicit_q_run.txt` part (B) and the docstring of
`recipe1_preflight_q_projection.py` lines 25–29.  TIDAL's
`theory.toml:118` uses **`RicciScalarCDT[]^2`**, which is the
**parity-even Riemann–Cartan Ricci scalar squared**, NOT the
Pontryagin density.  Agent H himself flagged in part (A) that the
Holst-Pontryagin scalar `R̃_H = ½·ε·R̃` is a **total derivative at
linear order on flat** (Nieh–Yan), so `(R̃_H)²` does not contribute a
bulk kinetic Lagrangian.  Agent H then *redefined* "R̃²" mid-script to
mean the parity-odd Pontryagin-type density `ε^{abef} R̃^{abcd}
R̃_{cd}^{ef}`.  TIDAL's b5·R̃² is neither of these — it is the
parity-even ordinary `(g^{μν} R̃_{μν})²`.

This redefinition is documented in Agent H's preflight transcript but
is **not flagged in the verdict JSON**.  The verdict therefore
applies, strictly, to a different operator from the one TIDAL
actually has.

---

## Section 3 · Wolfram-side projection logic (key code snippets)

### 3.1 ComponentDecompose.wl — direct component projection, no irreducible decomposition

```
ComponentDecompose.wl:64-70
  ExtractTensorComponent::usage =
    "ExtractTensorComponent[eom, field, chart, componentIndices,
     additionalFields, ...] returns a single scalar component of an
     equation of motion via TraceBasisDummy on a specific chart slot
     assignment."
```

Searching `ComponentDecompose.wl` for the keywords *"irreducible"*,
*"axial part"*, *"trace part"*, *"q-decomposition"*, *"HMMN"*,
*"projector"* returns **zero hits**.  The pipeline performs raw
component-by-component projection only.  The 24 `t_n` fields in the
JSON are simply the 24 independent components of an antisymmetric
(in last two indices) rank-3 tensor in 4-D:
`4 · binomial(4,2) = 4·6 = 24`.

`ComponentDecompose.wl:1067-1098` describes the term-by-term
projection branch used for R̃²-decomposed torsion theories:

```
(* Term-by-term projection (supervisor's approach, commit 4a89164). *)
…
(* individually through SplinterToArray (which uses ComponentArray …) *)
```

i.e. the projection is `Plus[term_1,…] · ComponentArray` — no
algebraic q/trace/axial split.

### 3.2 ExportJSON.wl — Hamiltonian sector filter EXCLUDES torsion

```
ExportJSON.wl:1651-1677
  (* Hamiltonian sector filter: exclude terms referencing torsion fields.
     $tidalHamiltonianFilter and $tidalTorsionHead are set in the WLS by
     _derive.py for torsion theories. This keeps only GW+EM self-energy
     and interaction terms — torsion self-energy and torsion cross-terms
     are dropped since they don't contribute to the conversion measurement
     C₀ = P/B₀² (which uses per-field self-energy only).
     See: hamiltonian_filter_design.md *)
  If[StringLength[torsionPertName] > 0,
    Module[{nBefore = Length[result], tPert, nAfter},
      tPert = torsionPertName;
      result = Select[result, Function[term,
        Module[{fieldA, fieldB, keep},
          fieldA = term[["factor_a", "field"]];
          fieldB = term[["factor_b", "field"]];
          keep = !StringMatchQ[fieldA, tPert ~~ "_" ~~ __] &&
                 !StringMatchQ[fieldB, tPert ~~ "_" ~~ __];
          keep ]]]; …
```

This is a **deliberate design choice**: the `canonical.hamiltonian_terms`
section of the JSON omits all torsion-sector terms, by name-matching on
the torsion perturbation prefix (`"t"`).  Consequence: any analysis
that reads `canonical.hamiltonian_terms` will see ZERO torsion
contribution to the energy, even though the underlying Lagrangian
sources real propagating torsion components (t_6, t_13, t_20).

This explains the apparent contradiction in Section 4 below.

### 3.3 ExportJSON.wl — kinetic-coefficient extraction

```
ExportJSON.wl:485-487
  (*   DO NOT divide — keep RHS unnormalized and record kinetic_coefficient_symbolic *)

ExportJSON.wl:607
  lhsStructure["kinetic_coefficient_symbolic"] = kineticCoeffStr
```

The `kinetic_coefficient_symbolic` field is the production ground
truth for "is this field standard-kinetic or Pais–Uhlenbeck?".  Time
order is on `lhs.order.time`, kinetic coefficient on
`lhs.kinetic_coefficient_symbolic`.  Both are populated for every
propagating field.

---

## Section 4 · Hamiltonian-term cross-reference

`canonical.hamiltonian_terms` has **78 entries**, with field-pair class
counts:

| pair  | count |
| ----- | ----- |
| (h,h) | 68    |
| (a,a) |  8    |
| (a,h) |  2    |

**No t-field appears in any hamiltonian term.**  The set of distinct
fields appearing in `factor_a`/`factor_b` across all 78 entries is
exactly `{a_0,a_1,a_2,a_3, h_0,…,h_9}` — the 14 non-torsion dynamical
fields.

The three terms with `mixed_2_*` operators (the PU signature):

```
h_4(mixed_2_0_0_0) × h_4(mixed_2_0_0_0)   coef_sym = -b5   class = self
h_7(mixed_2_0_0_0) × h_7(mixed_2_0_0_0)   coef_sym = -b5   class = self
h_9(mixed_2_0_0_0) × h_9(mixed_2_0_0_0)   coef_sym = -b5   class = self
```

These are precisely the `(d²_t h_n)²` Pais–Uhlenbeck self-energies of
the metric block, exactly matching the `d4_t(h_n)` equations in
Section 1.  No analogous `mixed_2_*` term appears for any t-field — but
that is **because of the explicit torsion filter** in
`ExportJSON.wl:1651-1677`, NOT because the underlying decomposition
fails to produce torsion-sector kinetic terms.

(Direct evidence: `equations[]` for t_6, t_13, t_20 have `lhs =
d2_t(t_n)` and a non-zero RHS proportional to b5 — the kinetic
energy on these fields exists in the EOM layer, but the Legendre
transform output that feeds `hamiltonian_terms` has been pruned of
torsion contributions by design.)

---

## Section 5 · **Verdict (a)** — schema-level PASS confirmed at production, with a substantive caveat

Of the three options:

> (a) **PASS confirmed at production**: t-fields all standard-kinetic;
>     no PU on torsion sector; Agent H's schema matches reality.
> (b) **PASS at schema only; production-level mismatch**: some t-fields
>     have time_order ≥ 4 with b5 promotion; tensor-q sectoral path
>     closed.
> (c) **Inconclusive.**

the production-level evidence supports **(a) with caveats**:

1. **Standard-kinetic at production (PASS)**: every propagating t-field
   has `time_derivative_order = 2`, with kinetic coefficient
   `(-25/2)·b5`.  No t-field has `time_order = 4`.  The PU promotion
   is confined to the metric `{h_4, h_7, h_9}` block, exactly as Agent
   H's schema predicts.

2. **No qualitative mismatch**: Agent H's schema (R̃ has 1 derivative
   of K, R̃² has 2 total ⇒ standard kinetic) is consistent with the
   TIDAL EOM-layer output (`d²_t` not `d⁴_t` on the propagating
   torsion).

3. **Caveat C1 — Operator mismatch**: TIDAL's b5·R̃² is parity-EVEN
   `(RicciScalarCDT[])²` (`theory.toml:118`).  Agent H's preflight
   analyzed two different objects: the parity-even Holst scalar
   `R̃_H = ½ ε·R̃` (Section A — vanishes by Nieh–Yan, so trivially
   standard-kinetic) and the parity-odd Pontryagin density (Section B
   — non-trivially standard-kinetic).  Neither of these is exactly
   the operator TIDAL has.  The "PASS" verdict therefore applies to
   the *parity-odd* Pontryagin operator, while TIDAL's actual operator
   is *parity-even* `(g^{μν} R̃_{μν})²`.  At linear order on flat the
   parity-even object likewise contains no `(∂²)²` mechanism on q
   (each R̃-Ricci is one ∂K, exactly as for the Holst case), so the
   conclusion still holds — but Agent H's verbatim verification does
   not cover the operator TIDAL is using.

4. **Caveat C2 — DOF count discrepancy**: Agent H's "16 q-DOF" is the
   4-D off-shell count; TIDAL's pipeline produces 3 propagating t-DOF
   after plane-wave reduction along z (`theory.toml:[reduction]`).
   The 21 algebraic t-components are pre-reduction q-irreducible
   modes that lose their kinetic source under the propagation
   restriction.  Not a contradiction, but the audit narrative would
   benefit from making this explicit so future readers do not expect
   "16 propagating torsion fields" in the JSON.

5. **Caveat C3 — Hamiltonian filter masks the torsion sector**: the
   `canonical.hamiltonian_terms` block has zero torsion entries, by
   the deliberate filter at `ExportJSON.wl:1651-1677`.  Anyone
   relying on `hamiltonian_terms` alone to verify "Recipe 1 PASS at
   production" would draw an incorrect blank — the right artefact is
   `equations[]`, not `canonical.hamiltonian_terms`.

### Concrete next steps

1. **Document caveat C1 in `FINAL_ASSESSMENT.md`** under the existing
   §"What is overstated or wrong" section: state that Agent H's
   verdict was sympy-verified for the parity-odd Pontryagin density,
   not for TIDAL's parity-even `(RicciScalarCDT[])²`.  Note that the
   conclusion (no `(∂²q)²` mechanism) carries over qualitatively
   because each R̃-factor is still one `∂K` at linear order, but the
   verbatim verification needs a small extension.

2. **Optionally extend the preflight script**
   (`recipe1_preflight_q_projection.py`) with a part (D) that performs
   the same derivative-counting exercise for `(g^{μν} R̃_{μν})²` at
   linear order on flat in the q-irreducible ansatz.  Expected result:
   identical conclusion — each R̃_Ricci is one ∂K, the square is two
   ∂'s, standard kinetic.  ~30 minutes of sympy work.

3. **Update the audit memo** to record the production-level
   confirmation: t_6, t_13, t_20 are `d²_t`/std-kinetic; h_4, h_7,
   h_9 are `d⁴_t`/PU.  This is the *factual* statement that
   distinguishes (a) and (b).

4. **Tensor-q sectoral path remains alive** at the level of the
   Recipe 1 production verdict, **still conditional on**:
   - parity-odd extension of Chatzistavrakidis–Ranjbar–Zekoč 2024
     (not present in published literature; per FINAL_ASSESSMENT §5),
   - addressing C1 (operator-mismatch nuance).

   The path is not closed by this investigation; it is also not
   *operationally important*, because TIDAL's headline observable
   (`h_5 ↔ a_1`) is in the standard-kinetic graviton sector and does
   not depend on the tensor-q torsion path being viable.

---

## Section 6 · Citations

### Source files
- `examples/torsion_gertsenshtein/theory.toml:118` — Lagrangian
  expression containing `RicciScalarCDT[]^2`
- `examples/torsion_gertsenshtein/theory.toml:136-138` — plane-wave
  reduction along z
- `examples/torsion_gertsenshtein/theory.toml:147-149` —
  `[perturbation] small_parameters = ["b5"]`
- `tidal/wolfram/ComponentDecompose.wl:64-70` —
  `ExtractTensorComponent` (direct component projection, no
  irreducible decomposition)
- `tidal/wolfram/ComponentDecompose.wl:1067-1098` — term-by-term
  projection branch for R̃² torsion theories
- `tidal/wolfram/ExportJSON.wl:485-487` — kinetic-coefficient
  extraction comment
- `tidal/wolfram/ExportJSON.wl:607` —
  `lhsStructure["kinetic_coefficient_symbolic"]` populated
- `tidal/wolfram/ExportJSON.wl:1651-1677` — Hamiltonian torsion-sector
  filter (the reason `canonical.hamiltonian_terms` shows no torsion)
- `research/perturbative_hamiltonian/scripts/recipe1_preflight_q_projection.py:25-29`
  — preflight docstring identifying R̃ as parity-odd Pontryagin
- `research/perturbative_hamiltonian/results/recipe1_explicit_q_run.txt`
  Section A — Holst scalar vanishing by Nieh–Yan
- `research/perturbative_hamiltonian/results/recipe1_explicit_q_run.txt`
  Section B — Pontryagin density confirmed (∂q)·(∂q)

### JSON paths (in `examples/data/torsion_gertsenshtein.json`)
- `equations[i].field` for i ∈ [0, 38) — 38 fields
- `equations[i].lhs.order.time` — time-derivative order; 4 only on
  h_4 (i=8 in fields[]), h_7 (i=11), h_9 (i=13); 2 on the four a_*,
  five h-block fields, and the three t-block fields t_6, t_13, t_20;
  0 elsewhere
- `equations[i].lhs.kinetic_coefficient_symbolic` —
  - `(-25*b5)/2` for t_6, t_13, t_20
  - `2*b5` for h_4, h_7, h_9 (the PU-promoted block)
  - `-kappa^(-2)` for h_5 (and the rest of the std-kinetic h-block)
- `canonical.hamiltonian_terms` — 78 entries; 0 reference any t-field
- `fields[i].name`, `fields[i].is_dynamical` — 38 dynamical fields:
  4 photon, 10 graviton, 24 torsion

### Independent validation
- `notes/FINAL_ASSESSMENT.md` lines 70-74 — audit's existing
  observation that h_5 is standard-kinetic and the headline
  observable doesn't gate on the b5-PU subspace; this document
  extends that observation downward to the torsion sector.
