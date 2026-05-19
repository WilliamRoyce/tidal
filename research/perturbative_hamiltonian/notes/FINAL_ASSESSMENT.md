# Final Consolidated Assessment

**Date**: 2026-04-26 (initial); **amended 2026-04-27** with Phase 2 lead
investigations (5 single-agent leads + synthesis).
**Scope**: 10-agent investigation (Rounds 1-3) + 3 critical reviews
(1, 2, 3) + 4 meta-reviews (K, L, M, N) + 5 Phase 2 leads
+ 2 correction agents (W1, W2) + 1 cleanup agent (W3) + a personal
BN-citation correction.
**Output**: All artefacts in
`research/perturbative_hamiltonian/{scripts,notes,results,reviews,meta_reviews}/`.
**Phase 2 verdict (added 2026-04-27)**: **(c) Definitive gap** — no
constructive published bridge for TIDAL's b5=0 critical surface in
PGT b5·R̃², but the phenomenon IS published in adjacent frameworks.
**Status**: Theoretical investigation closed.  Phase 3 (operational
closure) deferred per user direction.  Documentation corrections
(A-G) applied 2026-04-27.

> **Phase 2 amendments to read alongside this document**:
> `notes/phase2_synthesis.md` — convergent verdict from 5 lead investigations.
> `notes/lead_glavan_zlosnik_lin.md` — GZL verdict (c) NOT TRANSFERABLE; surfaced 3 published analogues.
> `notes/lead_recipe1_xact.md` — Recipe 1 (a) PASS at production with C1 operator-mismatch caveat.
> `notes/lead_lyakhovich_forward_citations.md` — (b) NO NEW LEAD; Aashish-Saif misattributed.
> `notes/lead_bn1983_direct_read.md` — BN 1983 PRD covers TIDAL class but explicitly defers; "40-year acknowledged-as-open".
> `notes/lead_bjhk_metric_affine.md` — (b) NO TRANSFER; 1911.08846 misattributed.

---

## Executive verdict

The constraint-promotion barrier in PGT b5·R̃² **is real** as a theoretical
problem, but it **does not block** TIDAL's flagship Gertsenshtein-conversion
measurement. The original investigation's framing significantly overstated
the strength of its positive results and misidentified the critical
sector. The audits then over-corrected in some places. The verified
picture below distinguishes what survives from what doesn't.

**Bottom-line recommendation**: ~3 days of operational work to ship the
project in a defensible state. Don't write the three Round 3 papers as
framed. Publication C (no-go theorems) is salvageable with revisions.

---

## What is genuinely true (survives all audits)

### Theoretical content

1. **Sectoral classification of the constraint-promotion phenomenon** —
   the 4-way decomposition (axial / trace / tensor-q / metric h₄,₇,₉)
   is a genuine intellectual contribution and survives all reviews.

2. **Three convergent no-go arguments** for the metric h₄/h₇/h₉
   Pais-Uhlenbeck subspace — Round 1 Agent A (Lyakhovich rank-jump
   `det(M) ∝ b5`), Agent D (reducible-Stückelberg structural argument),
   Round 2 Agent G (dual no-go for 2-form auxiliaries). These cover
   *local first-order auxiliary lifts*. Each individual argument has a
   minor caveat, but they converge on the same conclusion. Review 1
   confirmed via re-execution.

3. **VT integral converges for polynomial degree-1 source forms** at
   toy level — Round 2 Agent E + Round 3 Agent I sympy-verified.
   Voicu 2020 4D-GB pathology is genuinely avoided in the linearised
   regime.

4. **Helmholtz residue δE = 0 for variational source forms** — Round 1
   Agent B confirmed; this is a corollary of "the EOMs came from a
   Lagrangian", not a new theorem.

5. **The 25-year history of named-but-unsolved** — Blagojević-Nikolić
   (1983-84), Yo-Nester-Ni (1999, 2002), Blagojević-Cvetković (2018) —
   is real as historical framing. Meta-N caveat below: BC 2018's
   framework doesn't apply to TIDAL's higher-derivative case, so the
   "25-year unsolved for TIDAL's case" is more accurate than "25-year
   unsolved" simpliciter.

6. **The Bopp-Podolsky scalar-auxiliary lift for the axial sector** at
   linear-flat order is internally self-consistent (Agent C + Agent F
   sympy-verified, with Agent F's correction from vector to scalar
   auxiliary).

7. **Curtright Stückelberg construction (arXiv:2411.16928) at 1+1D toy
   level**: gauge invariance `δF̊ = 0` and rank uniformity
   `det(H_kin) = 1 - λ_a²` are sympy-verified (Agent J + Review 1 C3).

### Project-relevance content

1. **TIDAL's headline observable uses `--source h_5 --target a_1`**.
   `h_5` is `time_order=2, kinetic_coefficient_symbolic="-kappa^(-2)"`
   — a STANDARD-KINETIC graviton, **not constraint-promoted**. This
   is the single most important factual finding of the entire audit
   process (Meta-L). The `d²_t(h_5)` and `d²_t(a_1)` equations form a
   closed 2-equation sub-system; the b5-promoted block (h₄/h₇/h₉) does
   not enter the headline measurement.

2. **`tidal/measurement/_conversion.py` IS Hamiltonian-based**, not
   amplitude-based as Review 3 incorrectly claimed. It evaluates
   `canonical.hamiltonian_terms` and computes
   `P(t) = E_target(t) / E_source(0)` (Meta-L). The existing pipeline
   produces measurement-grade `P(t)` to ~1e-6 precision.

3. **Existing campaigns already answer the science question.** NULL
   amplification across multiple regimes:
   - 276-run dark-photon sweep (`dark_photon_amplification_campaign_v0.31.md`)
     matches Boccaletti baseline to 6.7×10⁻⁶
   - Propagating-PGT model (`propagating_model_finding.md`): A=1.0 within
     stability window
   - Stage A+B PGT campaign (`project_pgt_campaign_stageA_B_results.md`)
     closed with NULL physics

   The science finding is "no amplification at minimal / dark-photon /
   propagating couplings", not "we don't know yet".

---

## What is overstated or wrong (must be corrected)

### Original investigation overstatements

1. **"Krupka-Voicu Theorem 1 verified symbolically"** — there is no
   Theorem 1 in arXiv:1406.6646. The paper has Definition 1
   (canonical variational completion). What was "verified" is the
   tautology that an already-variational source form admits a
   canonical completion. Reviews 2 + Meta-K confirmed.

2. **"Path A resolves the barrier for full PGT b5·R̃²"** — Round 3
   Agent I's own Phase 6 concedes that L_VT inherits the parent
   Pais-Uhlenbeck structure for the metric subspace; the Legendre
   transform retains the rank-jump. Path A produces a Lagrangian for
   the EOM, not a Hamiltonian for the phase space. Review 1 confirmed.

3. **"Path B-trace established literature (Barker 2024)"** — Barker
   identifies trace torsion as the Weyl *gauge field*, not a
   Goldstone, and explicitly excludes parity-odd terms (the actual
   TIDAL b5·R̃² use case). The construction has not been done for
   parity-odd PGT. Review 2 + Meta-K confirmed.

4. **"Path B-tensor-q applies (Chatzistavrakidis-Ranjbar-Zekoč 2024)"**
   — CRZ handles only parity-even free fields and the m → 0 Goldstone
   limit. TIDAL b5·R̃² is parity-odd by construction; Agent J's
   translation requires a parity-odd extension that does not exist in
   published literature. Reviews 2 + Meta-K confirmed.

5. **"AM 2020 cross-validation strong"** — Aoki-Mukohyama's φ has
   non-canonical kinetic `(∂φ)²/(1+φ²)` with sinh-canonicalisation,
   and lives inside an infinite-tower dRGT-bigravity. The bare `𝒳²`
   term *breaks* AM's ghost-freedom. Agent F's
   `φ²/(2b)` Lagrange-multiplier construction is structurally different.
   "Cross-validation" should be downgraded to "qualitative consistency
   at linearised order". Review 2 confirmed.

6. **"Three publishable papers' worth"** — Publication A's headline
   contradicts its own Phase 6 result; Publication B has overstated
   literature claims; only Publication C (no-go theorems) is
   defensible, and even that needs (i) Agent D's structural argument
   tightened, (ii) Agent G's dual no-go generalised beyond 2-forms.

### NEW pathology found by Review 1 (not flagged by original investigation)

7. **L_VT diverges as M_c² → 0** (Review 1 C5 — sympy-verified).
   Path A's L_VT has explicit `1/M_c²` and `1/M_c⁴` poles. PGT
   critical-mass surfaces (Karananas 2014, Blagojević 2018) are
   exactly where Path A breaks. Round 3 Agent I did not flag this.

8. **F-J cross-validation OVERSTATED** (Review 1 C8). Agent F's
   `det(M_aux) = 1/b²` *diverges* at b=0 — the rank-jump is
   relocated, not removed. Agent J's `det(H_kin) = 1 - λ_a²` is
   genuinely b5-independent. They are NOT parallel mechanisms.

9. **Phase-space jump factor 3 vs 5 — convention drift** (Review 1
   C6). Round 3 Agent I's main script reports 6 → 18 (factor 3);
   Phase 6 v2 reports 6 → 30 (factor 5). The synthesis cited factor
   5 without flagging the convention.

### Review 2 over-corrections (Meta-K confirmed)

10. **"Cabo Bizet-Bartocci 2026 attribution fabricated"** — *Review 2's
    own invention*. The original investigation cites
    arXiv:2601.22007 (Aashish-Saif 2026), not 2602.12114, and not
    "Cabo Bizet-Bartocci". Meta-K grep'd the entire corpus and
    confirmed.

11. **"BC 2018 Appendix D verbatim quote hallucinated"** — *Review 2
    was wrong*. Meta-K retrieved the BC PDF via `pdftotext` and found
    the verbatim quote at lines 2685-2688. Quote is real. **HOWEVER**
    (per Meta-N), the original investigation truncated the next
    sentence which inverts the framing: BC's parametric diagonalisation
    is singular at b̄=0, but the physics is recovered by reading F
    directly. Appendix D is constructive, not a no-go.

### Review 3 over-corrections (Meta-L confirmed)

12. **"`_conversion.py` is amplitude-based"** — Review 3 wrong; it's
    Hamiltonian-based via `compute_energy_timeseries`. Meta-L
    confirmed.

13. **"TIDAL is polynomial-coefficient"** — Review 3 wrong; TIDAL's
    `_resolve_symbolic_coeff` uses `eval` and accepts rational
    parameter expressions like `1/(2*kappa^2)`. Meta-L confirmed.

14. **"The metric h_4/h_7/h_9 subspace is THE Gertsenshtein
    bottleneck"** — Review 3 wrong. The actual headline channel is
    `h_5 ↔ a_1`, which is standard-kinetic. h_4/h_7/h_9 promotion
    affects total system energy conservation but does not gate the
    conversion observable. **This is the single most important
    correction.**

### Meta-N: BC Appendix D structure is different from TIDAL's case

15. **BC Appendix D's framework does not apply to TIDAL** (Meta-N).
    Three structural reasons:
    - Wrong derivative order: BC is 2nd-order Dirac-ADM; TIDAL is
      4th-order Ostrogradsky from b5·R̃²
    - Wrong constraint topology: BC's `cn → 0` MAKES a constraint
      primary; TIDAL's `b5 → 0` BREAKS algebraic constraints into
      propagating modes — inverse limits
    - BC's `b5` is NOT TIDAL's `b5`; naming collision is coincidental

---

## Concrete documentation corrections

### `docs/tex/perturbative_reduction_constraint_barrier.tex`

1. **Replace "Krupka-Voicu Theorem 1"** with "Krupka-Voicu Definition 1
   (canonical variational completion)" throughout.
2. **Repair the BC 2018 Appendix D quote**: keep the verbatim quote
   but ADD the constructive follow-up sentence:
   > "However, since the matrix F for b̄ = 0 is already diagonal, the
   > critical parameters cn can be obtained directly from F. The same
   > conclusion also holds for the form of H⊥ᶠ."
   And reframe: BC's Appendix D is constructive for *its* 2nd-order
   case but does not apply to TIDAL's higher-derivative Ostrogradsky case.
3. **Verify or remove the Cabo Bizet-Bartocci citation** if it
   appears anywhere — the original investigation appears to have
   cited Aashish-Saif arXiv:2601.22007 instead, but double-check the
   TeX writeup. (Review 2 invented this audit target; Meta-K confirmed
   it doesn't exist in original notes.)
4. **Downgrade Path B-trace** from "established literature" to "open
   research conditional on parity-odd Barker extension".
5. **Downgrade Path B-tensor-q** from "applies with caveats" to
   "conditional on parity-odd extension of Chatzistavrakidis-
   Ranjbar-Zekoč 2024 that does not exist in published literature".
6. **Add Voicu linearity-in-highest-derivatives gate** as a separate
   preflight, distinct from VT integral convergence. The two are
   independent failure modes per Voicu 2020 §4 + Appendix A.
7. **Add the M_c² → 0 caveat** to Path A's verdict — L_VT diverges at
   constraint-mass critical surfaces (Review 1 C5).
8. **Tone down the F-J cross-validation** to "qualitative consistency
   at linearised order; structurally different mechanisms" (Review 1
   C8).
9. **Add the operational primary section**: document
   `tidal/measurement/_conversion.py` as the Hamiltonian-based
   conversion observable for affected theories. Cite existing
   campaign findings as evidence the observable produces
   measurement-grade `P(t)`.

### Memory file `perturbative_reduction_hamiltonian_lit.md`

Append a "verified picture" entry summarising the audit results
(sectoral classification real; metric subspace genuinely blocked for
local first-order auxiliaries; amplitude/Hamiltonian conversion
observable is the operational primary; existing campaigns already
delivered the science answer).

### GitHub issue #321

Update body to reflect: the recommended path forward is *documentation
of the limitation*, not pursuit of a Stückelberg recipe. Cite the
review findings.

---

## Recommended next steps (Meta-M's plan, supported by all audits)

### Step 1 (1 day, lowest risk) — Documentation corrections + ship

(a) Apply the 9 TeX corrections above.
(b) Add the operational-primary section for `_conversion.py`.
(c) Update `_check_perturbative_consistency` CLI hint to point users
    to `--what conversion` rather than `--what energy` for theories
    with `mixed_2_*` operators.
(d) Update issue #321 with the verified picture.
(e) Append memory-file entry.

**Risk**: trivial. **Payoff**: very high — converts a "stuck research
project" into a "science-completed project with documented limitation".

### Step 2 (½ day + 1h HPC) — Validate existing campaign findings

(a) Run `t_end`-independence sweep on `torsion_gertsenshtein` (P(t) vs
    P(2t) at each parameter point — the `A(2t)/A(t) ≈ 1` diagnostic
    per CLAUDE.md).
(b) Grid-resolution convergence sweep at fixed parameters (128, 256,
    512, 1024).
(c) Cross-check Boccaletti analytical formula at b5=0 sub-limit.

**Risk**: moderate (could uncover tachyonic-instability artefact,
which is itself useful information). **Payoff**: high — produces
measurement-grade evidence that NULL amplification findings are
robust.

### Step 3 (1 day + 30 min HPC) — Light-mediator regime sweep

50-point logarithmically-spaced 1D sweep on the surviving `h_5 ↔ a_1`
channel: chosen mediator mass from `m → 10⁻³` to `m ~ 1`. Measure
P_max, plot A = P_max/P_GR.

**Risk**: moderate. **Payoff**: high if signal exists; medium if
null.

### Sequencing

Step 1 first → Step 2 → Step 3. Total ~3 days agent time + ~2 hours
HPC.

After this the flagship deliverable is in a defensible, documented,
measurement-grade state.

---

## What NOT to do

1. **Do not write the three Round 3 papers as currently framed.**
   - Publication A: headline contradicts its own Phase 6
   - Publication B: overstated literature claims
   - Publication C: salvageable, but consumes 4-6 weeks for a result
     that doesn't unblock TIDAL

2. **Do not extend TIDAL JSON schema for jet-5 operators** (would
   contradict the LPS invariant and validator design; pointless given
   the headline observable doesn't need it).

3. **Do not pursue Tier-3 research lines** (asymptotic resummation,
   Wilsonian integration, BV-BFV homological reduction) — high
   failure risk, no near-term payoff for project goals.

4. **Do not chase parity-odd extensions of Chatzistavrakidis-
   Ranjbar-Zekoč or Barker** unless the project later needs the
   axial/trace/tensor-q sectors specifically. They are not on the
   critical path for the headline observable.

---

## Optional Tier-2 work (only if supervisor asks)

If the supervisor specifically requests publication-grade theoretical
output:

A. **Rescue Publication C**: tighten Agent D's structural argument to
   handle the starting-point ambiguity Review 2 §4 flags; generalise
   Agent G's dual no-go beyond 2-form auxiliaries. Then write up as a
   no-go letter. ~4-6 weeks.

B. **Forward-citation lead from Meta-N**: investigate Glavan-Zlosnik-
   Lin 2024 (arXiv:2311.17459, "Hamiltonian analysis of metric-affine
   R²"). Methodologically closest published paper to TIDAL's problem.
   Different framework but worth a focused read if the project later
   needs Path A or Path B work to advance.

Both are research-value-only; neither advances TIDAL's flagship physics
question.

---

## Manual-retrieval bookkeeping

Per `MANUAL_RETRIEVAL_NEEDED.md`:

- BC 2018 Appendix D **resolved** (Meta-K + Meta-N retrieved via `pdftotext`)
- Tier 1 remaining: Blagojević-Nikolić 1983 (PRD 28:2455, pre-arXiv),
  Blagojević-Nikolić 1983 (Nuovo Cim B **73**:258, doi:10.1007/BF02721794) — earlier mis-cited as Nuovo Cim B 84:25 (1984).  User-retrieved 2026-04-27 to `literature/BF02721794/`.  The actual paper is on $R + T^2$~PGT (linear in $R$, quadratic in torsion), NOT $R+R^2+T^2$ as previously claimed.
  These are needed to verify the "25-year unsolved" history claim
  rigorously. Manual retrieval recommended before any Publication C
  citation.
- Yo-Nester-Ni 1999/2002 should be arXiv-accessible.
- Hehl-McCrea-Mielke-Ne'eman 1995 (gr-qc/9402012) arXiv-accessible.
- Tier 2 items (pre-arXiv classics) are not load-bearing for the
  surviving publishable content.

---

## Summary in one paragraph

**The 10-agent investigation produced real mathematical content but
its headline framing significantly overstated what was achieved. The
3-review audit then over-corrected in a few places, prompting 4
meta-reviews that pinned down the verified picture. The metric
h₄/h₇/h₉ Pais-Uhlenbeck subspace is genuinely blocked for local
first-order auxiliary lifts (3 convergent no-gos), but it does not
gate TIDAL's headline observable (`h_5 ↔ a_1` is standard-kinetic).
The existing `_conversion.py` pipeline is Hamiltonian-based and
produces measurement-grade `P(t)`. The science question has been
answered by existing campaigns (NULL amplification across multiple
regimes). The right next move is ~3 days of operational work
(documentation corrections + validation sweeps + light-mediator
sweep), not the academic publication trajectory the original
synthesis recommended. Don't write the three papers as framed; only
Publication C is defensible after substantial revision and consumes
4-6 weeks for a result that doesn't unblock anything.**

---

## Audit trail

All artefacts in `research/perturbative_hamiltonian/`:

- `notes/round1_synthesis.md`, `round2_synthesis.md`, `round3_synthesis.md`
  — the original investigation's claims
- `notes/round2_agentE/F/G_*.md`, `round3_agentH/I/J_*.md` — per-agent
  writeups
- `scripts/*.py` — 18+ sympy scripts (all re-execute cleanly per
  Review 1)
- `reviews/review1_mathematical_verification.md` — Review 1 + 8
  C-checks
- `reviews/review2_literature_interpretation.md` — Review 2
- `reviews/review3_project_relevance.md` — Review 3
- `reviews/scripts_review/*.py` — 8 new C-check scripts
- `meta_reviews/meta_review_K_literature_claims.md` — verifies Review 2
- `meta_reviews/meta_review_L_pipeline_claims.md` — verifies Review 3
- `meta_reviews/meta_review_M_next_steps.md` — synthesis + plan
- `meta_reviews/meta_review_N_BC_AppendixD_content.md` — BC Appendix D
  deep dive (the user's lead, vindicated with twist)
- `MANUAL_RETRIEVAL_NEEDED.md` — papers that need physical/library
  access
- `notes/FINAL_ASSESSMENT.md` — this document

Total agent-hours: ~10 hours wall-clock across all 14 agents.
Cost-effective for the depth of audit produced.
