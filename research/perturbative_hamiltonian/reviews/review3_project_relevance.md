# Review 3 — Project Relevance & TIDAL Integration Audit

**Reviewer**: critical project-relevance review of Rounds 1-3 (10 agents, 18 sympy
scripts) on the constraint-promotion barrier in PGT b5·R̃².
**Date**: 2026-04-26
**Posture**: project's advocate, not the investigation's apologist. Verify
whether toy-model results translate to TIDAL's *actual* physics goal —
graviton-photon conversion probability in PGT torsion theories.

---

## TL;DR for the user

**Path A (Vainberg-Tonti) and Path B (sectoral Stückelberg) as currently
formulated DO NOT solve TIDAL's physics problem.** They are
academically interesting results, but neither delivers the operational
deliverable TIDAL needs: a `hamiltonian_terms` block on
`examples/data/torsion_gertsenshtein.json` that the existing
`tidal/measurement/_energy.py` evaluator can consume to compute
graviton-photon conversion probabilities. The investigation has
produced **three publishable papers' worth of theoretical content**
(structural classifications, no-go theorems, sectoral recipes) but
**zero TIDAL-pipeline-usable artefacts**. The metric h₄/h₇/h₉
Pais-Uhlenbeck blocker — which is the only one of the four sectors
that matters for the Gertsenshtein effect — is precisely the one all
three rounds confirmed remains genuinely unsolved.

If the project goal is "publish theoretical-physics papers about the
constraint-promotion barrier", the investigation succeeded.
If the project goal is "produce a working Hamiltonian for
torsion_gertsenshtein.json so we can measure conversion probability",
the investigation **did not advance us**.

---

## Q1 audit — Does Path A (Vainberg-Tonti) give us something USABLE?

**Verdict: NO.**

What VT actually produces (Round 3 Agent I, Phase 4-5):

- A 339-term polynomial Lagrangian L_VT for a **3+3 toy** modelling the
  metric h₄/h₇/h₉ subspace.
- Krupka-Voicu Theorem 1 verified: EL(L_VT) = ε exactly.
- **L_VT inherits the Pais-Uhlenbeck higher-derivative structure of the
  parent theory.** Highest jet order is 5 (canonical IBP form) at b5≠0.
- **Phase-space dimension jumps from 6 to 30** as b5 turns on — the
  Ostrogradsky rank-jump that motivated this whole investigation.

What TIDAL needs (`tidal/measurement/_energy.py` lines 1080-1395, plus
`tidal/symbolic/json_loader.py` HamiltonianTerm schema):

- `canonical.hamiltonian_terms` is a flat list of bilinear products
  `coeff · op_a(field_a) · op_b(field_b)`, with operators drawn from
  `{identity, gradient_*, laplacian, laplacian_*, mixed_T_S1_..., ...}`.
- The energy evaluator treats these as **classical phase-space
  observables**: factor_a and factor_b are evaluated on stored
  field/velocity snapshots.
- Operators with time_order ≥ 2 (`mixed_2_0_0_0`, ...) are resolved
  via `_compute_acceleration_from_eom` — i.e., by **re-evaluating the
  EOM at snapshot time**. This is the precise issue the TeX writeup
  flags as "not legitimate canonical observables".

The mismatch:

1. **L_VT is 339 terms with jet-order-5 dependence.** TIDAL's JSON
   schema cannot represent this. The schema's deepest operator is
   `mixed_T_S1_S2_...` which encodes *one* time derivative count and
   *one* spatial multi-index. Encoding `Y_0·Y_5` (zeroth derivative
   times fifth time derivative) in `hamiltonian_terms` would require
   either (a) extending the schema to support `mixed_5_0_0_0`-style
   operators (which the validator at `_validate.py:189-207` would then
   *correctly* flag as constraint-promotion residue, defeating the
   purpose), or (b) rewriting the energy evaluator to handle Routhian
   denominators of the form `2·M_1²·M_2²·M_3²` and to
   evaluate up-to-fifth time derivatives. Neither has been attempted.
2. **The Routhian projector denominator 2·M_1²·M_2²·M_3² introduces
   field-rational structure** that TIDAL's polynomial-coefficient
   pipeline does not support. This is not a parser limitation; it is
   a fundamental schema assumption — the entire xAct/Wolfram pipeline
   produces polynomial coefficients.
3. **Even if we lifted L_VT into the schema**, evaluating it for energy
   would require sampling 5th-order time derivatives of the dynamical
   fields. The simulator's state stores `(q, v, A_0)` — there is no
   well-defined `∂_t^5 q` from snapshot data without re-evaluating the
   EOM five times in succession, which is the "EOM-defines-Hamiltonian"
   antipattern flagged in §pr-cb-conclusion.

**Bottom line on Path A**: it solves a Lagrangian-side existence
question (variational completion exists) but does not produce a
TIDAL-pipeline-usable Hamiltonian. Agent I's own Phase 6 verdict
states this explicitly: "Path A's L_eff inherits the rank-jump." The
Hamiltonian-side rank-jump is the entire reason the project needs a
working Hamiltonian; VT does not remove it.

**What VT *would* give us if we plumbed it in**: a "Lagrangian for the
already-PS-reduced y-only EOMs." But TIDAL already runs the
PS-reduced EOMs successfully on the simulation side
(`tidal/solver/perturbative_driver.py`). The thing TIDAL does not have
is **a Hamiltonian observable on the (q, π) phase space of the
already-reduced system**. VT gives us a Lagrangian for the same
EOMs, not a Hamiltonian. The Legendre transform of L_VT puts us right
back at the rank-jump.

---

## Q2 audit — Does Path B (sectoral Stückelberg) give us something USABLE?

**Verdict: NO for the actual `torsion_gertsenshtein` theory.**

What Path B delivers (Round 3 synthesis):

| Sector | Recipe | Status |
|--------|--------|--------|
| Axial | Bopp-Podolsky scalar-aux | ✅ verified (Agent F) |
| Trace | Barker conformal embedding | ✅ literature |
| Tensor-q | Curtright Stückelberg | ✅ verified with caveats (Agent J) |
| **Metric h₄/h₇/h₉** | **None** | **❌ blocked (Agent A+D+G no-gos)** |

What `torsion_gertsenshtein.json` actually contains (verified via
direct JSON inspection):

- 38 fields total: `a` (4 photon components), `h` (10 metric-perturbation
  components), `t` (24 torsion components).
- **The constraint-promoted fields are h_4, h_7, h_9 — metric
  perturbation components, not torsion.** The JSON LHS has
  `kinetic_coefficient_symbolic = "2*b5"` for these three equations
  with `time_order = 4`. The Hamiltonian terms 75/76/77 have
  `factor_a.field = h_4/h_7/h_9` with `mixed_2_0_0_0` operators. This
  is exactly the structure the TeX writeup describes (Eq.
  pr-cb-Mb5).
- **Three of 24 torsion components have `time_order = 2` with
  `kinetic_coefficient_symbolic` mentioning b5** — i.e. b5-dependent
  but standard-second-order kinetics. These are *not* constraint
  promotion in Path B's tensor-q sense; they are second-order
  dynamical fields whose mass scales with b5.

The Path B sectoral coverage maps onto the JSON as:

- Axial sector: covers AT MOST the t-fields' axial projection (a
  single scalar Goldstone plus the axial vector). The JSON has 24
  t-fields representing the *full* torsion tensor; the axial
  projection is ~4 of these.
- Trace sector: covers the trace projection — another ~4 of the
  t-fields, in principle.
- Tensor-q sector: covers the 16-component q (2,1) Curtright tensor —
  another ~16 of the t-fields. The Curtright recipe assumes
  standard-kinetic q (Recipe 1 PASS), so this is *applicable* only to
  the b5·R̃²-induced (∂q)² standard-kinetic terms. The JSON's actual
  3 b5-dependent torsion equations may or may not align with this —
  not verified.
- **Metric h₄/h₇/h₉**: NOT covered by any Path B sector. This is the
  showstopper.

The Gertsenshtein effect requires the **trace and non-TT modes of the
metric perturbation** (per the TOML comment "NO TT gauge on h —
preserves trace for R̃² graviton-torsion mixing (TT kills this
channel)"). Specifically, h_0/h_4/h_7/h_9 are some of the diagonal
metric components in the chosen index conventions; h_4/h_7/h_9 are
exactly the ones promoted by b5.

**These are precisely the modes Path B explicitly cannot handle.** If
Path B "solves" axial + trace + tensor-q (covers 24 torsion + 4
photon = 28 fields) but does *not* solve the metric h₄/h₇/h₉ subspace
(3 of 10 metric fields), then the photon-graviton conversion
amplitude — which is mediated by a quadratic h-coupling to the F²
photon term and an h-h coupling through R̃² — **cannot be cleanly
expressed as a Hamiltonian observable in the Path B recipe**.

**Bottom line on Path B**: covers the parts of the theory that don't
matter for the headline observable, leaves the parts that do matter
("metric Pais-Uhlenbeck") explicitly unsolved.

---

## Q3 audit — Toy-model faithfulness to PGT

**Verdict: the toys capture the *barrier* but not the *physics*. This
matters because Path A/B publishability is about the barrier; TIDAL's
goal is about the physics.**

Concrete mismatches:

### Round 1/2 toys (Agent A, B, E, F)

The toy Lagrangians are (paraphrasing):

```
T2: L = ½(∂_tφ)² + ½(∂_tχ)² − ½m²φ² − λφh − λ_χ χ ∂_t h
       − ½h² + ½b5(∂_t²h)²
T3: T2 plus a second constraint-promoted h field
T4: 3+3 (3 dynamical y_a, 3 constraint-promoted h_c)
```

What's faithful to PGT:

- The `½b5(∂_t²h)²` term IS the linearized projection of `b5·R̃²` onto
  the metric subspace (after component decomposition + IBP).
- The constraint-promotion mechanism (M_0 = 0, M_1 ∝ b5) IS what TIDAL
  hits at b5·R̃² = 0.
- The 3-field cardinality of T4 matches h₄/h₇/h₉.

What's NOT faithful:

1. **Gauge structure absent.** PGT has both diffeomorphism gauge (4
   parameters) and local Lorentz gauge (6 parameters). The toys are
   fully gauge-fixed. The Gertsenshtein effect's amplification scaling
   `P = sin²(κB₀D/2)` depends on the gauge-invariant transverse
   modes; whether Path A's L_VT respects diffeomorphism invariance
   has not been investigated. (Agent I's joint VT consistency check
   verified Helmholtz residue but not diffeomorphism invariance of
   L_VT.)
2. **Cross-sector couplings absent.** In real PGT, the metric h
   couples to torsion T through R̃ (linearly) and through R̃² (quadratically),
   AND to the photon a through F². The toys represent each sector in
   isolation. The actual amplification physics depends on the
   metric-torsion-photon *coupling chain* — which the toys don't
   model. So even if Path B works for isolated sectors, the
   *cross-sector* Hamiltonian remains unaddressed.
3. **q-sector "standard kinetic" claim is parity-odd-incomplete.**
   Agent J's caveat C1 flags 38 ε·DT·DT cross-terms not covered by
   Chatzistavrakidis-Ranjbar-Zekoč. The actual JSON uses
   `RicciScalarCDT[]^2` which is parity-odd-mixed by virtue of
   torsion-Riemann's parity-odd Holst component. Recipe 1 PASS holds
   only for the parity-even projection of the q-sector; the
   parity-odd part is **explicitly open research** per Agent J.
4. **Recipe 1 PASS may not hold for the actual metric coupling.**
   Agent H verified that q-irreducible projection of b5·R̃² is
   standard-kinetic for q. But the actual JSON shows **b5·R̃² promotes
   *metric* components h₄/h₇/h₉ to time_order=4, not torsion
   components**. This is the metric-PU subcase, not the q-Curtright
   subcase. The investigation correctly identifies this (Round 2
   Agent G's "Subcase A vs Subcase B" reframing) — but the
   *conclusion* is that Path B does not address Subcase A. The toys
   correctly identified the barrier; they do not bridge it.

The toys are **faithful as falsifiers of proposed solutions**. Agent
A/D's no-go theorems on toys reliably predict the failure of those
methods on the full theory. The toys are **not faithful as templates
for constructive solutions to the full theory** because they exclude
gauge structure, cross-sector couplings, and parity-odd content.

---

## Q4 audit — Are h₄/h₇/h₉ really constraint-promoted in PGT?

**Verdict: YES, confirmed by direct JSON inspection.**

From `examples/data/torsion_gertsenshtein.json`:

```
field h_4: time_order=4, kinetic_coefficient_symbolic="2*b5"
field h_7: time_order=4, kinetic_coefficient_symbolic="2*b5"
field h_9: time_order=4, kinetic_coefficient_symbolic="2*b5"
```

The Hamiltonian terms (75/76/77) are:
```
coeff=1.0, factor_a=h_X with mixed_2_0_0_0, factor_b=h_X with mixed_2_0_0_0
```

i.e. `(∂_t² h_X)²` — exactly the Pais-Uhlenbeck signature.

These are **components of the metric perturbation h_{μν}** (rank-2
symmetric tensor, 10 independent components in 4D after
component decomposition: h_0..h_9). h_4/h_7/h_9 are three specific
spatial-spatial components in the chosen TT-violating gauge. The
TOML comment confirms: "NO TT gauge on h — preserves trace for R̃²
graviton-torsion mixing".

**Not torsion components.** The torsion has its own t_0..t_23 fields
(24 components per TOML's PGT count: 4×3 + ... irreducible
decomposition). Only 3 of the 24 t-fields are b5-dependent at
order=2 (standard-kinetic with b5-mass).

The constraint-promotion barrier IS the metric subspace. The Path B
sectoral classification (axial / trace / tensor-q / metric) labels
the metric piece as the unsolved one. **The investigation correctly
identifies the showstopper sector and correctly claims it cannot be
bridged.**

---

## Q5 audit — Parity-odd Curtright extension caveat

**Verdict: relevant to torsion_gertsenshtein but a SECONDARY blocker.**

The JSON's parent Lagrangian uses `RicciScalarCDT[]^2`, which after
component decomposition includes parity-odd contributions through the
Holst-like `R̃ ⊃ ε·DT` projection. Per Agent J's caveat C1, 38
ε·DT·DT cross-terms in the q-sector are NOT covered by Chatzistavrakidis-
Ranjbar-Zekoč 2024. So Path B-tensor-q would require:

1. Either a parity-odd extension of the Curtright Stückelberg
   construction (Agent J flagged as "open research, 1-2 weeks").
2. Or a proof that these terms vanish on-shell (not yet attempted).

For the torsion_gertsenshtein theory specifically, the parity-odd
content matters because R̃² inherently contains parity-odd
contractions (via the dual εR ε R structure). So this is **not** a
caveat that goes away for our theory; it is a real gap.

However, this is *secondary* to the metric h₄/h₇/h₉ blocker. Even if
the q-sector parity-odd extension is published next month, the
metric subspace remains unsolved. So this caveat doesn't change the
overall verdict — Path B as a whole still doesn't deliver the
Hamiltonian for the full theory.

---

## Q6 audit — What does the TeX writeup say is needed?

**Verdict: the TeX writeup's conclusions are essentially unchanged by
the Round 1-3 results. The writeup correctly identifies the barrier
and the writeup's "what solving looks like" criteria are NOT met by
Path A or Path B.**

The writeup's §pr-cb-future enumerates four required properties for a
"successful" reduction:

> (i) H depends only on (q, π), with no Ostrogradsky higher momenta.
> (ii) Hamilton's equations from H reproduce the Parker-Simon-reduced EOM.
> (iii) The reduction is well-defined for the multi-field gauge-theory
> case relevant to PGT torsion.
> (iv) The construction is theory-agnostic, not requiring per-theory
> fine-tuning.

Path A scoring:

- (i): **FAIL.** Path A's L_VT inherits PU structure, so its Legendre
  transform has Ostrogradsky higher momenta (phase space jumps from
  6 to 30). Agent I's Phase 6 confirms this.
- (ii): **PASS** at the toy level (Krupka-Voicu Theorem 1 verified).
  Generalisation to full PGT not verified (38 fields, gauge structure).
- (iii): **PARTIAL.** T4 is 3+3 with no gauge structure; full PGT has
  diffeomorphism + Lorentz gauge. Agent I's "Does NOT generalise"
  caveats #1/#3 explicitly flag this gap.
- (iv): **PASS.** VT is theory-agnostic.

Path B scoring:

- (i): **PASS** for axial, trace, tensor-q. **FAIL** for metric.
- (ii): **PASS** at the toy level for the three covered sectors;
  cross-sector Hamiltonian (which is what the conversion observable
  actually needs) not addressed.
- (iii): **PARTIAL** — toys are sector-isolated; cross-sector recipe
  not built.
- (iv): **PASS** for the three covered sectors.

The writeup's existing conclusion (§pr-cb-conclusion) is:

> Simulation alone is not measurement-ready: a canonical Hamiltonian
> observable must be a function of phase-space data (q, π) only, but
> the pre-LPS b5·R̃² JSONs carry mixed_2_* operators in their
> hamiltonian_terms, which encode ∂_t² q factors that are defined by
> the EOM, not by the canonical state. These are not legitimate
> canonical observables.

**The Round 1-3 investigation does not change this conclusion for
the theories that matter.** What it adds:

- A clearer sectoral classification (axial/trace/tensor-q/metric).
- Three convergent no-go theorems for the metric subspace
  (publishable as a stand-alone result, doesn't help us).
- Three constructive recipes for sectors that aren't the bottleneck.

The writeup's three "research directions" — (a) algebraic-constraint
chain-rule expansion, (b) Parker-Simon Lagrangian analogue, (c)
antiderivative reconstruction — are **not addressed** by the
investigation. Path A is direction (d) (Krupka-Voicu + Stückelberg);
Path B is direction (e) (Lyakhovich-style sectoral Stückelberg).
Neither (d) nor (e) was claimed to be guaranteed; the writeup said
both were promising leads. The investigation has determined that
both leads, when pursued, **bottom out at the metric h₄/h₇/h₉
subspace**, which the no-go theorems already covered.

The writeup does need a small update to record the no-go theorem
strengthening (Round 2 Agent G's dual no-go) and the explicit
phase-space-dimension count from Agent I (6 → 30, factor 5). This
upgrades "metric subspace blocked by Agent A+D" to "metric subspace
blocked by Agent A+D+G + quantitative Ostrogradsky count". A useful
clarification but not a change in conclusion.

---

## Q7 audit — Project goal alignment

**Verdict: misaligned. The investigation pursued the academic problem
(can we build a perturbative Hamiltonian recipe?) rather than the
project problem (how do we measure conversion probability for the
torsion_gertsenshtein theory?).**

The user's stated requirements:

1. **"Without [a working Hamiltonian], we cannot do anything with the
   simulations."** Path A and Path B do not deliver a working
   Hamiltonian for the full theory. The constraint-promoted metric
   sector remains unbridged.

2. **"The torsion-Gertsenshtein conversion is the headline physics
   goal."** Conversion measurement requires energy of the photon
   sector vs. graviton sector. The graviton-photon coupling goes
   through h_{0i} and h_{ij} components — including h_4/h_7/h_9 (the
   constraint-promoted ones). So the conversion observable
   *requires* a clean Hamiltonian for the constraint-promoted
   subspace.

3. **"TIDAL is theory-agnostic and shouldn't restrict to specific
   subclasses."** Path B is explicitly sectoral; it would require
   per-sector Wolfram-pipeline integration. Even if Path B were
   fully implemented (3 sectors × 1-2 weeks each + integration), it
   would only handle theories where the constraint-promotion sector
   is axial/trace/tensor-q. The torsion_gertsenshtein theory's
   *metric* constraint-promotion is exactly the case Path B doesn't
   handle.

What an operationally useful investigation would look like:

1. **Compute the conversion probability *despite* the missing
   Hamiltonian.** The user has noted (CLAUDE.md §"t_end independence
   test") that the EOM-based simulation does produce P(t)
   timeseries, just without a true canonical-observable
   interpretation. A pragmatic answer: identify which conversion
   observables can be computed without a canonical Hamiltonian.
2. **Identify a workaround for the energy measurement of the
   metric subspace.** For example: project onto TT modes only
   (deliberately *gauge-discard* the constraint-promoted modes for
   energy measurement), and use the gauge-violating h_4/h_7/h_9
   modes only for tracking conversion *amplitude*. This is what the
   pre-LPS JSONs effectively do; it works but is "EOM-defines-energy"
   rather than legitimate.
3. **Numerically benchmark vs analytical Boccaletti formula** to
   determine whether the EOM-based P measurement is trustworthy
   despite the Hamiltonian gap. The user's notes report
   `P = sin²(κB₀D/2)` confirmed to 0.04% in the gertsenshtein
   (no-torsion) example — which suggests the conversion observable
   is robust against the Hamiltonian issue when the conversion is
   small. The campaign findings note `propagating_model_finding.md`
   (A=1.0, zero amplification at boundary) and
   `dark_photon_amplification_campaign_v0.31.md` (NULL amplification)
   indicate the *physics* answer for these theories is "no
   amplification", which makes the Hamiltonian question less
   pressing for *answering the science question* — but the user
   wanted measurement-grade rigor.

The investigation's three publishable papers (per Round 3 synthesis)
are:

- A: VT applied to PGT b5·R̃² → "we have a Lagrangian for the EOM, but
  it inherits PU structure for the metric subspace".
- B: Sectoral Stückelberg → "we cover 3 of 4 sectors; the 4th is
  blocked".
- C: Three convergent no-go theorems → "the metric subspace cannot be
  perturbatively reduced".

**These are publishable physics-research-frontier results. They are
not TIDAL operational deliverables.**

---

## Specific showstoppers for project goals

### Showstopper #1: L_VT cannot be plumbed into TIDAL's JSON schema

The JSON `hamiltonian_terms` schema supports operator names like
`mixed_T_S1_S2_...` with T = time-derivative count and S_i = spatial
order along axis i. The schema does not support:

- Routhian denominators (rational-coefficient terms) — TIDAL is a
  polynomial-coefficient pipeline.
- Per-term time-derivative orders > 2 — `_resolve_time_derivative`
  in `_energy.py:430-464` returns `None` for time_order ≥ 3.
- Schema violations of the LPS invariant are *flagged as errors* by
  `_check_perturbative_consistency` (lines 189-207).

To use Path A's L_VT, one would have to either (a) extend the schema
to allow up-to-jet-5 operators (which then defeats the purpose of LPS
— we'd be writing what LPS is meant to eliminate), or (b)
post-process L_VT into a sequence of auxiliary fields plus
second-order operators. (b) is essentially Path B but for the metric
sector — and Path B has been shown to fail there.

### Showstopper #2: Path B doesn't cover the metric sector

The Gertsenshtein effect's physics is mediated by metric perturbations
coupling to F² (the (R̃)² → metric → photon channel). The Path B
sectoral coverage explicitly excludes the metric h₄/h₇/h₉ subspace,
which is where the constraint promotion happens.

### Showstopper #3: Cross-sector couplings are unmodelled

Even if Path B were extended to cover the metric subspace (somehow,
contradicting the no-go theorems), the recipe is sectoral. The actual
torsion_gertsenshtein Hamiltonian has cross-sector couplings:

- h ↔ photon a (through F²)
- h ↔ torsion t (through R̃)
- h ↔ h (through R̃²)

The cross-couplings are not isolated to one sector. A sector-by-sector
recipe needs a glue that handles cross terms. None of the agents
addressed this.

### Showstopper #4: Operationalisation cost

Even setting aside the theoretical blockers, *implementing* Path A or
Path B in TIDAL would require:

- Wolfram-side: extending `PerturbativeReduction.wl` to emit the new
  operator forms.
- Schema-side: extending `json_loader.py` HamiltonianTerm to handle
  new term structures.
- Energy-side: extending `_energy.py` to evaluate the new operators
  (with consistent gauge-invariance handling).
- Test-side: re-deriving the affected example JSONs and verifying
  conservation tests.

Per the per-feature work cost in CLAUDE.md, this is multiple weeks
of careful pipeline work — and would only be worth it if the
underlying theoretical recipe gave a Hamiltonian we could trust. The
investigation has not produced such a recipe for the metric sector.

---

## Toy-to-PGT inferential gap quantification

| Aspect | Toys | PGT b5·R̃² | Gap |
|--------|------|------------|-----|
| Field count | 2-6 | 38 | 6× to 19× scale-up |
| Gauge structure | None | Diff + Lorentz | Untested |
| Cross-sector coupling | None | h-t-A all coupled | Untested |
| Parity-odd content | None | ε·R·R contractions | Untested for q-sector |
| Spatial dimension | 0 (time-only) or 1+1D | 1+1D after plane-wave reduction | Probably OK |
| Christoffel structure | None | Flat at 0th order, perturbative | Verified flat-Minkowski |
| Constraint-promoted sector | h_c (PU) | h_4/h_7/h_9 (PU) | **MATCHED** — and confirms the no-go |

The barrier identification is robust (toys correctly identify the
PGT barrier). The constructive recipes are NOT robust; they are
sector-isolated and the actual theory has cross-sector physics.

---

## What additional work would make the investigation USEFUL for TIDAL?

**Tier 1: short-term, makes existing code more honest (1-2 weeks)**

1. Update `docs/tex/perturbative_reduction_constraint_barrier.tex`
   to reflect Round 1-3 findings:
   - Strengthen the no-go from 2 (Agent A + D) to 3 convergent
     arguments (+ Agent G dual no-go).
   - Add the explicit Ostrogradsky phase-space count (6 → 30).
   - Cite the new arXiv references (Chatzistavrakidis-Ranjbar-Zekoč
     2024; the FJ Schur-bordering 2026 paper).
   - State explicitly that VT and sectoral Stückelberg are not
     bridges across the metric subspace barrier.
2. Extend `_check_perturbative_consistency` to *additionally* warn
   that for the affected theories, the conversion measurement should
   be cross-validated against a non-Hamiltonian observable (e.g.
   raw amplitude ratios, t_end-independence checks).

**Tier 2: pragmatic engineering for the science (2-4 weeks)**

3. **Build a non-Hamiltonian conversion measurement pathway.** The
   "P_target / P_source" amplitude ratio is gauge-invariant and does
   not require a Hamiltonian. Document this as the canonical
   path for constraint-promotion theories until a true Hamiltonian
   becomes available. The existing `tidal/measurement/_conversion.py`
   already does this; we just need to flag it as the
   *primary* observable for affected theories rather than treating
   energy as primary.
4. **Project-onto-TT-modes energy variant**: an alternative energy
   measurement that explicitly excludes h_4/h_7/h_9 from the energy
   sum, with a clear documented caveat that it loses
   "constraint-promoted-mode" content. Strictly an underestimate but
   gauge-honest.

**Tier 3: research with real project payoff (months)**

5. **Operationalise Path B for the actually-covered sectors.** If the
   axial/trace/tensor-q recipes can be implemented in
   `PerturbativeReduction.wl` (3-4 weeks of Wolfram work + tests), we'd
   gain *partial* Hamiltonian coverage for theories where the
   constraint promotion is in those sectors. This doesn't help
   torsion_gertsenshtein (metric promotion) but might help future
   theories. **Caveat:** this is research investment for a
   not-yet-justified target — the user has not requested this.
6. **Tackle the metric subspace blocker directly via
   non-Stückelberg approaches.** The investigation explicitly closed
   off Stückelberg (irreducible + reducible). Approaches not yet
   tried: (a) Asymptotic-series resummation of the constraint chain,
   (b) Constraint-mode integration-out at the path-integral level
   (Wilsonian), (c) Direct Dirac-Bergmann at b5≠0 with the b5=0
   limit handled as a separate Hamiltonian. None of these is
   "engineering"; all are research with high failure risk.

**Tier 4: things NOT to do**

7. Do *not* try to make L_VT consumable by extending TIDAL's schema
   to allow jet-5 operators. This contradicts the LPS invariant and
   the validator design.
8. Do *not* publish Path A as "we have a working Hamiltonian for PGT
   b5·R̃²" — that's a misrepresentation. It's a Lagrangian that
   inherits PU structure; the Legendre transform has rank-jumps.
9. Do *not* redirect the project budget to writing the three Round 3
   papers without an explicit decision that publication is more
   valuable than continued physics work. The investigation generated
   genuinely novel physics-research-frontier content but the user's
   stated priority is the conversion observable.

---

## Final recommendation

The Round 1-3 investigation **definitively answered** an academic
question:

> "Can the constraint-promotion barrier in PGT b5·R̃² be bridged by
> known perturbative-Hamiltonian techniques?"

The answer is: NO for the metric subspace, YES for sectoral cases that
aren't bottlenecks for TIDAL.

This is a publishable physics-frontier result and should be written up
**as such** — not as "we solved TIDAL's barrier". The TeX writeup
already correctly frames the situation; minor updates suffice.

For the project's *actual* goal (torsion-Gertsenshtein conversion
measurement), the investigation has **not advanced us**. The
recommended next steps are operational/pragmatic (Tier 1-2 above),
not the further theoretical work the Round 3 synthesis recommends.

The user should be told: **the academic problem is closed; the
project problem remains open and is best addressed by accepting the
non-Hamiltonian conversion observable as the primary measurement
pathway for the affected theories**. This is what the existing
codebase already does (`measurement/_conversion.py` is amplitude-
based; energy is secondary). Document this explicitly, ship the
publication-grade no-go theorem paper if the user wants, and move on
to the science.

---

## Cross-references

- TeX writeup: `docs/tex/perturbative_reduction_constraint_barrier.tex`
- Pipeline integration points:
  - `tidal/wolfram/PerturbativeReduction.wl` — LPS Throw mechanism
  - `tidal/wolfram/ExportJSON.wl` — JSON schema (mixed_T_S1_... operators)
  - `tidal/symbolic/json_loader.py` — HamiltonianTerm class,
    canonicalize_kinetic_for_perturbation
  - `tidal/measurement/_energy.py` — energy evaluator, _resolve_time_derivative,
    _evaluate_mixed_factor
  - `tidal/measurement/_conversion.py` — amplitude-based conversion
    measurement (the operationally usable pathway)
  - `tidal/cli/_validate.py` — _check_perturbative_consistency
- Round summaries:
  - `notes/round1_synthesis.md`
  - `notes/round2_synthesis.md`
  - `notes/round3_synthesis.md`
- Theory under review:
  - `examples/torsion_gertsenshtein/theory.toml`
  - `examples/data/torsion_gertsenshtein.json`
- Related project memory:
  - `~/.claude/projects/-workspaces-torsion-gertsenshtein/memory/torsion_gertsenshtein_findings.md`
  - `dark_photon_amplification_campaign_v0.31.md` (NULL amplification result)
  - `propagating_model_finding.md` (A=1.0 boundary)
