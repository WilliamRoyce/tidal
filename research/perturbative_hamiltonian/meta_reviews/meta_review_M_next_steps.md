# Meta-Review M — Synthesis and Next-Steps Planning

**Date:** 2026-04-27
**Author:** Synthesis-and-planning meta-review agent
**Scope:** Synthesise verified findings from Reviews 1/2/3 and Meta-Reviews K/L; deliver a concrete next-step plan that aligns with TIDAL's actual physics-measurement goal (graviton-photon Gertsenshtein conversion), not the academic-publication trajectory the original investigation drifted toward.
**Posture:** Conservative; assume the strongest critical claims of Reviews 2 and 3 will be confirmed.

---

## Section 1 — Verified picture (synthesis)

### 1.1 What is genuinely true and useful

After three rounds of investigation (10 agents, ~18 sympy scripts) and three review passes, the following findings survive critical audit and are worth keeping:

1. **The 25-year history is real.** Blagojević–Nikolić (1983–84), Yo–Nester (1999–2002), and Blagojević–Cvetković (2018) each independently named what is structurally the same phenomenon: a small parameter promoting an algebraic-constraint field to a higher-derivative dynamical field, with a discontinuous Hamiltonian structure across the critical surface. This historical anchor is well-evidenced and useful as framing.

2. **Sectoral classification is genuine intellectual content.** Round 2's reframing of "tensor sector" into two structurally distinct subcases — Subcase A (metric h₄/h₇/h₉ Pais–Uhlenbeck) and Subcase B (q-irreducible torsion standard-kinetic) — is a non-trivial clarification that the literature did not previously articulate cleanly. The four-sector picture (axial / trace / tensor-q / metric) is a real contribution, even if the published-literature handling claimed for three of the four is over-stated (see §1.2).

3. **Three convergent no-go arguments for the metric subspace** (Round 1 Agent A's Lyakhovich rank-jump `det(M) ∝ b5`; Round 1 Agent D's reducible-Stückelberg parameter-counting; Round 2 Agent G's "dual no-go" for 2-form auxiliaries) constitute a robust negative result. They each have minor caveats (Agent D's argument depends on a specific starting-point choice — Review 2 §4 — and Agent G is specific to 2-form auxiliaries — Review 1 C8), but they converge on the same conclusion: **no first-order Stückelberg lift bridges b5 = 0 in the metric Pais–Uhlenbeck subspace**. The convergence is the load-bearing claim, not any individual argument.

4. **Helmholtz residue δE = 0 for variational source forms.** Round 1 Agent B's symbolic check that the Parker–Simon-reduced source forms have vanishing Helmholtz residue is correct and useful. Caveat: this is essentially a corollary of the fact that the parent EOMs *came from a Lagrangian* — not a new theorem, but a consistency check.

5. **VT integral convergence for polynomial source forms.** Round 2 Agent E's verification that the Vainberg–Tonti integral converges for degree-1-polynomial source forms (which TIDAL produces by virtue of being in the linearised-quadratic regime) is technically correct. The framing as "Krupka–Voicu Theorem 1" is wrong (there is no Theorem 1; only Definition 1) — but the substantive content is sound at toy level.

6. **Quantitative phase-space dimension count for the metric subspace.** Round 3 Agent I's Phase 5 finding — phase-space dimension jumps from 6 (b5 = 0) to either 18 (jet-counting convention) or 30 (`2·N·r` convention) at b5 ≠ 0 (Review 1 flags the convention dependence) — provides a quantitative version of the Ostrogradsky rank-jump. The exact factor is convention-dependent and should be stated as such.

### 1.2 What was overstated or wrong

Per Review 2 (literature) and Review 3 (project relevance), the following claims need retraction or substantial revision:

| Claim | Status | Source |
|---|---|---|
| "Krupka–Voicu Theorem 1 verified" | Misnomer — KV has Definition 1, not Theorem 1 | Review 2 §1 |
| "VT diverges for negative fibre homogeneity" | Threshold is degree ≤ −1, not "negative" generally | Review 2 §2 |
| "AM 2020 cross-validation of axial BP" | Substantially weaker — AM's φ has non-canonical kinetic, lives in a dRGT-bigravity tower with infinite curvature corrections; AM explicitly states the bare 𝒳² term *breaks* ghost-freedom | Review 2 §6 |
| "Path B-trace established literature (Barker 2024)" | **Wrong** — Barker identifies trace torsion as the Weyl *gauge* field (not Goldstone); parity-odd terms (the actual TIDAL use case) explicitly excluded as future work | Review 2 §9 |
| "Path B-tensor-q applies (Chatzistavrakidis-Ranjbar-Zekoč 2024)" | Conditional on parity-odd extension that does NOT exist in published literature; CRZ handles only parity-even free fields, m → 0 Goldstone limit (not the b5 → 0 parameter-discontinuity limit) | Review 2 §7 |
| "BC 2018 Appendix D verbatim quote" | **Unverified, likely hallucinated**; paper's framing of Appendix D suggests it *enables* the limit analysis, opposite of the alleged quote | Review 2 §8 |
| "Cabo Bizet–Bartocci 2026" (arXiv:2602.12114) | **Author attribution fabricated** — actual authors Chan-López, Martín-Ruiz, Cabrera, Paulin Fuentes | Review 2 §10 |
| "Three publishable papers' worth of content" | Premature; two of three would not survive peer review without major revision | Review 2 + Review 3 |
| "Constraint-promotion barrier substantially resolved" | Misleading — Path A and Path B do not bridge the metric h₄/h₇/h₉ subspace, and that subspace is precisely the showstopper for `torsion_gertsenshtein` | Review 3 |

### 1.3 What was genuinely useful but mis-framed

A few results are mathematically correct but their *project relevance* was over-sold:

- **Path A (Vainberg–Tonti) at toy level** correctly produces a polynomial Lagrangian L_VT for a 3+3 metric-PU toy, but L_VT inherits the parent PU structure: its Legendre transform has the same rank-jump that motivated the investigation. **Path A produces a Lagrangian for the EOM, not a Hamiltonian for the phase space.** This is what Round 3 Agent I's Phase 6 already concedes.

- **Path B-axial (Bopp–Podolsky)** is self-consistent at linear-flat order *for the axial sector in isolation*, but the curved-background extension fails at O(h²) (Agent F's own caveat) and AM cross-validation is qualitative-only (Review 2). It is a real but narrow construction.

- **Recipe 1 PASS (q-projection is standard-kinetic)** is verified against a sympy *schema* of the q-projection, not against the actual xAct `b5·R̃²` decomposition (Review 1 §"Round 3 — Agent H" caveat; Review 2 §"Cross-cutting D"). For the actual TIDAL theory the q-projection is irrelevant anyway — `torsion_gertsenshtein.json` carries metric h₄/h₇/h₉ promotion, not q-promotion.

### 1.4 The critical project-relevance finding (Review 3)

**Even granting the most generous reading of the technical claims, the investigation has not advanced the project's actual goal.** This is the load-bearing claim of Review 3 and it survives independently of every literature-attribution issue:

- `examples/data/torsion_gertsenshtein.json` carries the constraint promotion in **metric components h₄/h₇/h₉**, not torsion components (Review 3 Q4, verified via direct JSON inspection).
- The Path B sectoral coverage (axial / trace / tensor-q) explicitly excludes the metric subspace.
- Path A's L_VT is jet-order-5, has Routhian denominators (`2·M_1²·M_2²·M_3²`), and is fundamentally incompatible with TIDAL's JSON `hamiltonian_terms` schema — which is polynomial-coefficient and supports only `mixed_T_S1_...` operators with T ≤ 2.
- Even if L_VT could be plumbed into the schema, the Hamiltonian observable would still need to evaluate up-to-fifth time derivatives of dynamical fields from snapshot data — exactly the "EOM-defines-Hamiltonian" antipattern the TeX writeup correctly flags as illegitimate.

The investigation has therefore *characterised* the barrier with high resolution but has *not bridged* it for the affected theory.

### 1.5 Specific corrections needed

Concrete changes the project should make:

**`docs/tex/perturbative_reduction_constraint_barrier.tex`:**

1. The verbatim BC 2018 Appendix D quote (currently lines ~117–121 of the file) is unverified and possibly hallucinated. **Either** read the actual PDF and confirm the quote (preferred), **or** rephrase to "Blagojević & Cvetković's Appendix D handles vanishing critical parameters by case-by-case branch analysis rather than as a perturbative limit; the perturbative limit is not constructed there."
2. Replace any "Krupka–Voicu Theorem 1" → "Krupka–Voicu Definition 1 (canonical variational completion)".
3. Direction (e) currently treats Lyakhovich Stückelberg as the leading constructive direction. Add a paragraph noting that the multi-front investigation (this directory) found three convergent no-go arguments specifically blocking Stückelberg-lifting (irreducible *and* reducible) for the metric h₄/h₇/h₉ subspace at b5 = 0. Direction (e) is therefore *not viable* for the metric sector; it is viable for axial-only and provisionally for tensor-q (subject to a parity-odd extension that does not yet exist in the published literature).
4. The reference to "CaboBizetBartocci2026" (line ~578) needs the author list updated to "Chan-López, Martín-Ruiz, Cabrera, Paulin Fuentes, arXiv:2602.12114". The technical claim ("FJ termination requires Schur-complement non-degeneracy") is still defensible as a 2026 preprint result.
5. Add a new paragraph in §pr-cb-future explicitly stating that the multi-front investigation (research/perturbative_hamiltonian/) clarified the sectoral structure — three sectors (axial, trace, tensor-q) admit at least partial Stückelberg-style constructions in the published literature, but the metric h₄/h₇/h₉ subspace, which is the bottleneck for the `torsion_gertsenshtein` family, does not.

**Memory files** (`/home/vscode/.claude/projects/.../memory/`):

6. The currently-empty entry `perturbative_reduction_v6_complete.md` and `rtilde_decomposition.md` should be cross-linked to a new file documenting that the constraint-promotion barrier remains genuinely open for the metric subspace and that the project's operational answer is the amplitude-based conversion observable (see Section 2).

7. Add a new memory entry `perturbative_hamiltonian_review_summary.md` summarising:
   - The sectoral classification is real.
   - The three convergent no-go theorems for the metric subspace are robust.
   - Path A and Path B do not produce a TIDAL-pipeline-usable Hamiltonian for `torsion_gertsenshtein`.
   - The amplitude-based conversion observable is the right operational primary.

**Issue #321** (which currently tracks "perturbative reduction constraint barrier"):

8. Update the issue body to reflect the verified picture: the barrier is genuinely open for the metric subspace, and the recommendation has shifted from "pursue Stückelberg lifting" to "accept the non-Hamiltonian conversion observable as primary; document as known limitation".

---

## Section 2 — Project status assessment

### 2.1 The actual TIDAL physics goal

TIDAL exists to compute graviton–photon Gertsenshtein conversion probabilities `P(t) ≈ sin²(κB₀D/2)` for theories of the form

```
L = (1/κ²) R̃ + α₁ I₁ + α₂ I₂ + α₃ I₃ + b₅ R̃² − ¼ F²
```

and variants (dark photon, plasma, propagating torsion, non-minimal coupling), and to identify whether *any* parameter region or modification produces measurable amplification beyond the standard Boccaletti baseline.

The deliverable is a number — `P_max` or an amplification factor `A = P/P_GR` — verified to be a robust physical observable, ideally measurement-grade (i.e. derived from a true canonical Hamiltonian rather than from EOM trajectories alone).

### 2.2 What blockers actually remain

**Blocker 1 — the metric Pais–Uhlenbeck subspace.** Genuine, robust, characterised. h₄/h₇/h₉ acquire fourth-order time derivatives at b5 ≠ 0; b5 = 0 is a singular limit. No Stückelberg lift bridges it. This is what Review 3 calls "the showstopper" and what Reviews 2 and 1 confirm is correctly identified.

**Blocker 2 — schema incompatibility of L_VT with TIDAL JSON.** Even if Path A's Lagrangian were the right answer, the schema cannot represent jet-5, rational-coefficient operators. Extending the schema is multiple weeks of work and would itself defeat the purpose of LPS.

**Blocker 3 — operational, not theoretical: the Hamiltonian-side energy measurement is not the only conversion observable.** This is Review 3's most useful pragmatic insight. The amplitude-based ratio `P(t) = E_target(t) / E_source(0)` implemented in `tidal/measurement/_conversion.py` does *not* require canonical phase-space data to be measurement-grade. It requires only field amplitudes and is gauge-invariant when the source/target are gauge-invariant scalars (which they are: photon energy and TT-graviton energy in the linearised regime).

### 2.3 Workarounds that exist or could exist

The codebase already contains the right operational primary for affected theories. From the file listing:

- `tidal/measurement/_conversion.py` — amplitude-based conversion (P_max, peak conversion).
- `tidal/measurement/_dispersion.py` — modal-eigenvalue dispersion measurement.
- `tidal/measurement/_mixing.py` — channel-resolved mixing analysis.
- `tidal/measurement/_resonance.py` — parameter-space resonance detection.

These can be used today on the affected `torsion_gertsenshtein` theory without any Hamiltonian-side resolution. The existing campaigns (276-run dark-photon sweep, propagating-model A=1.0 finding) used exactly these and produced robust *physical* answers — namely "no amplification" — that did not depend on the Hamiltonian gap.

The TeX writeup's current §pr-cb-conclusion frames this as "simulation alone is not measurement-ready". That framing was correct *if* the energy-based observable is the only legitimate one. **The amplitude-based observable changes the framing.** A cleaner statement is: *energy-based observables are not measurement-ready for affected theories; amplitude-based observables are, and are the operational primary.*

### 2.4 Where the project actually stands

- The flagship `torsion_gertsenshtein` theory simulates correctly (Parker–Simon EOM-side reduction works).
- The amplitude-based conversion measurement works and has produced robust null-amplification results across the explored parameter space (`dark_photon_amplification_campaign_v0.31.md`, `propagating_model_finding.md`).
- The Hamiltonian-side measurement is genuinely blocked for the metric subspace, but this does not block the *primary observable*.
- Stages A and B of the PGT campaign are closed with NULL physics results (`project_pgt_campaign_stageA_B_results.md`).

The state of the project is: **science questions answered (no amplification at minimal/dark-photon/propagating couplings), with a documented limitation on the secondary energy-based observable**. This is a respectable scientific outcome, not a stuck project.

---

## Section 3 — Ranked next-step recommendations

Three concrete next-steps, ranked by payoff-to-effort and risk.

### Next-step 1 (highest priority, lowest risk) — Document the operational picture and ship

**Goal:** Make the project's actual operational state legible to the supervisor and to future readers (including the user's future self), so that the science output is recognisable as complete rather than as work-in-progress.

**Concrete deliverables:**

(a) Edit `docs/tex/perturbative_reduction_constraint_barrier.tex` per the seven specific corrections in §1.5 above. Estimated 2–3 hours.

(b) Add a new section `\subsection{Operational primary: the amplitude-based conversion observable}` to the same TeX file (or to `docs/tex/gertsenshtein.tex`) explicitly documenting that `tidal/measurement/_conversion.py` is the canonical primary for affected theories, with `_energy.py` as a sanity check that does not run on h₄/h₇/h₉. Reference the existing campaign findings (`dark_photon_amplification_campaign_v0.31.md`, `propagating_model_finding.md`) as evidence the amplitude observable produces robust physical answers. Estimated 2 hours.

(c) Open or update GitHub issue #321 to reflect the verified picture. Close the implicit "find a Stückelberg recipe" sub-thread; open a new "document amplitude observable as primary" sub-thread that the deliverables (a, b) close. Estimated 30 minutes.

(d) Add CLI hint text in `tidal/cli/_validate.py` `_check_perturbative_consistency`: when the validator detects `mixed_2_*` operators in `hamiltonian_terms` for a theory in the constraint-promotion family, surface a hint pointing to the amplitude-based conversion measurement (`tidal measure ... --what conversion`) as the right operational observable, rather than `--what conservation` or `--what energy`. Estimated 1–2 hours.

**Total effort:** ~1 full working day. **Risk:** trivial — these are documentation and one CLI hint. **Payoff:** very high — converts a "stuck research project" into a "science-completed project with one documented limitation".

### Next-step 2 (high priority, medium risk) — Convergence and t_end-independence checks for the existing observable

**Goal:** Ensure the amplitude-based conversion observable used by the existing campaigns is genuinely measurement-grade, not just operationally convenient. This is what Review 3 §Q7 calls "numerically benchmark vs analytical Boccaletti formula".

**Concrete deliverables:**

(a) Run a t_end-independence sweep on `torsion_gertsenshtein` across the same parameter region as the closed Stages A/B (which used the dark-photon and Einstein–Cartan sub-models). At each parameter point measure `P(t)` and `P(2t)`; the ratio `A(2t)/A(t) ≈ 1` is the diagnostic that the conversion measurement is genuine (per CLAUDE.md "t_end independence test"). Wall time: 1 HPC job × 2 t_end values, ~40 minutes on sapphire `--qos=INTR`. **File to write:** `examples/torsion_gertsenshtein/HPC_TEND_CONVERGENCE.md` documenting the recipe.

(b) Run a grid-resolution convergence sweep at fixed parameters: grid-shape ∈ {128, 256, 512, 1024} with all else equal. Verify `P_max` converges at the rate expected for the chosen FD order (4× per doubling for 2nd-order FD; machine-precision for spectral). The result is a single convergence-rate plot. **File:** `examples/torsion_gertsenshtein/HPC_GRID_CONVERGENCE.md`.

(c) Cross-check against the Boccaletti analytical formula `P = sin²(κB₀D/2)` at the b5 = 0 sub-limit (pure Einstein–Maxwell). The existing `gertsenshtein` example (no torsion) has confirmed this to 0.04%; replicate at b5 = 0 within the `torsion_gertsenshtein` JSON to confirm the Hamiltonian gap doesn't bias the amplitude observable in the b5 → 0 limit. Estimated 1 hour locally.

**Total effort:** ~half a day of agent work + ~1 hour HPC wall time. **Risk:** moderate — could uncover that the amplitude observable is *itself* sensitive to t_end at certain parameter points (would indicate a tachyonic-instability artefact, not amplification), but this is exactly what the diagnostic is designed to catch and is itself useful information. **Payoff:** high — produces measurement-grade evidence that the existing null-amplification findings are robust, closing the "is the answer real or a numerical artefact?" question.

### Next-step 3 (medium priority, medium risk) — Spectral-conversion-coefficient sweep across light-mediator regime

**Goal:** Sharpen the "no amplification" answer. The literature (`amplification_literature.md`) says amplification at the stability boundary is light-mediator enhancement (1/m²). The current campaigns scanned coarse mass parameters but did not specifically push toward the light-mediator regime. A focused 1D scan in the small-mass limit would confirm or rule out enhancement in the operationally accessible region.

**Concrete deliverables:**

(a) Identify, from `examples/torsion_gertsenshtein/theory.toml`, which mass parameter (α₁, α₂, α₃, b5) acts as the "torsion mediator mass" in the dispersion relation for the channel that survives Stages A/B (the `h_5 ↔ a_1` channel per `torsion_gertsenshtein_findings.md`). 1–2 hours of code reading + dispersion-measurement on a single point.

(b) HPC sweep: 50 points logarithmically spaced in the chosen mediator mass, from `m → 10⁻³` (where light-mediator enhancement should peak per literature) to `m ~ 1` (where it should saturate at the Boccaletti baseline). Measure `P_max` at each. Wall time: < 30 min on sapphire. **File:** `hpc_results/<jobid>/light_mediator_sweep.csv`.

(c) Plot `P_max(m)` and `A = P_max/P_GR`. Either the curve shows the predicted 1/m² enhancement (publishable amplification result) or it stays flat (closes the operational hypothesis space). Either outcome is a useful finishing touch on the science.

**Total effort:** ~1 day agent work + ~30 min HPC. **Risk:** moderate — could uncover an amplification signal that was hidden in the coarse-grid sweeps, which would need careful follow-up (re-derivation, t_end-independence, separate HPC sweep) before publication. **Payoff:** high if signal exists; medium if null (confirms existing campaigns).

### Next-step 4 (lower priority, low risk) — Optional non-minimal coupling exploration

The user has memory entries (`literature_critical_analysis.md`) noting "three paths beyond torsion independence: non-minimal coupling, ghost-free kinetics, cubic PGT". The current `torsion_gertsenshtein` theory uses minimal coupling and produces null amplification because of the polarisation block-diagonal structure (`torsion_gertsenshtein_findings.md`). A non-minimal coupling `δ R̃[μν] F^μν` was explored in `propagating_model_finding.md` and gave A=1.0 (zero amplification within the narrow stability window). The remaining unexplored direction is **non-propagating non-minimal coupling on a stable PGT background** (the "nonminimal model" that memory says "is the correct framework").

If the project has remaining bandwidth after Steps 1–3, a focused HPC sweep on the non-minimal coupling model would be the natural next physics question. But this is *additional science*, not a fix for the existing pipeline.

### Sequencing

Run Step 1 first (documentation, ~1 day). It unblocks Steps 2 and 3 by making the operational primary explicit.

Run Step 2 second (convergence checks, ~½ day). This validates the existing campaign findings without exploring new physics.

Run Step 3 third (light-mediator sweep, ~1 day). This sharpens the science finding to publication readiness.

Step 4 is optional and depends on the supervisor's appetite for additional exploration after the existing campaigns close.

**Total path-to-completion:** ~3 working days agent time + ~2 hours HPC wall time. After this the project's flagship deliverable — graviton–photon Gertsenshtein conversion in PGT — is in a defensible, documented, measurement-grade state.

---

## Section 4 — What NOT to do

### 4.1 Do not pursue the Round 3 three-publication trajectory in its current form

The Round 3 synthesis recommended:

- Publication A: Path A applied to PGT b5·R̃²
- Publication B: Sectoral Stückelberg recipe
- Publication C: Three convergent no-go theorems

After Reviews 2 and 3:

- **Publication A** would not survive peer review without retracting "Krupka–Voicu Theorem 1" (no such theorem), retracting the AM cross-validation framing (qualitative consistency, not equivalence), separately verifying the Voicu linearity-in-highest-derivatives gate, and acknowledging the L_VT inherits PU structure (defeating the headline claim). What remains after retractions is a technical note "we computed the canonical variational completion of the Parker–Simon-reduced source forms in a 3+3 toy and confirmed Helmholtz residue vanishes" — interesting but not a research publication on its own.

- **Publication B** would need to retract Path B-trace (Barker doesn't support it for parity-odd) and downgrade Path B-tensor-q to "conditional on a parity-odd extension that doesn't exist". What remains is "Bopp–Podolsky single-auxiliary lift for axial torsion at linear-flat order" — a narrow result that has known curved-background failure at O(h²) and that AM 2020 already covers (with a different construction) for the same physical mode.

- **Publication C** is the most defensible — three convergent no-go arguments against the Stückelberg lift in the metric Pais–Uhlenbeck subspace — but would still need (i) the BC 2018 verbatim quote verified or removed, (ii) Agent D's structural argument tightened to handle the starting-point ambiguity Review 2 §4 flags, (iii) Agent G's dual no-go generalised beyond 2-form auxiliaries. With those, it could become a publishable letter or short note. Estimated effort: 4–6 weeks of careful theorem-statement and proof-tightening, plus rereading the actual BC 2018 PDF.

**Recommendation:** if any paper is written, write only Publication C, with the corrections above. Skip A and B until the underlying constructions are extended to cover the actual TIDAL use case (the metric subspace), which the no-go theorems indicate is unlikely to happen via Stückelberg-style approaches.

**Stronger recommendation:** *don't write any of the three papers right now*. The project's primary deliverable is the conversion-observable physics result, not a Hamiltonian-reduction methods paper. Writing C would consume 4–6 weeks of agent time on a result that does not unblock TIDAL and that the supervisor has not requested.

### 4.2 Do not pursue Tier-3 research lines (asymptotic resummation, Wilsonian integration)

Review 3 §Tier 3 lists:

- Asymptotic-series resummation of the constraint chain.
- Constraint-mode integration-out at the path-integral level (Wilsonian).
- Direct Dirac–Bergmann at b5 ≠ 0 with the b5 = 0 limit handled as a separate Hamiltonian.

Each is high-risk research with uncertain timeline (months to years) and no existing literature precedent for the specific PGT b5·R̃² metric subspace. Even successful completion would not necessarily produce a TIDAL-schema-compatible Hamiltonian.

If the user wants to spend agent budget on long-shot research, the highest-payoff alternatives are:
- Re-attempt the perturbative-reduction problem using **non-Stückelberg approaches not covered by the no-go theorems** (e.g. asymptotic resummation), but accept up-front that this is months of work with low success probability.
- Pivot to a **non-minimal coupling theory** that avoids the constraint-promotion structure entirely (Step 4 above).

### 4.3 Do not extend TIDAL's JSON schema to allow jet-5 operators

Review 3 §Showstopper #1 is correct: extending the schema to support up-to-jet-5 `mixed_T_S1_..._SD` operators would defeat the purpose of LPS (the validator at `_check_perturbative_consistency` correctly flags such operators as constraint-promotion residue) and would introduce a class of "Hamiltonians" that are not legitimate canonical observables. Don't go down this path.

### 4.4 Do not poll squeue or run multiple wolframscripts in parallel

Standard CLAUDE.md rules; flagged here only because the next-step plan involves HPC submissions and additional derivations. Use `hpc_shuttle.sh wait` for file-existence polling; never `squeue` in a loop.

---

## Section 5 — Honest verdict

### Should the user write up papers?

**No, not the three Round 3 papers as currently framed.** The literature attribution errors (Barker, BC 2018, Cabo Bizet–Bartocci) and the Path B-trace / Path B-tensor-q claim of "established literature" would not survive peer review. Publication A's headline claim (Path A "resolves" the barrier) is contradicted by Path A's own Phase 6 result (L_VT inherits PU structure). Publication B is sectorally narrow and partly speculative. Publication C is the most defensible but still needs 4–6 weeks of careful tightening.

**If a paper is to be written**, write only the no-go theorem result (Publication C, tightened) — but only *after* the operational pipeline work (Steps 1–3 above) is complete and the project's flagship deliverable is in a defensible state. The no-go theorem result is *useful* (tells future researchers this direction is closed for the metric subspace) but it is not *urgent*.

### Should the user pivot to operational work?

**Yes.** The investigation has answered the academic question (the metric Pais–Uhlenbeck subspace is genuinely blocked) and has clarified the structural picture (sectoral classification, three convergent no-go arguments). What remains is operational: document the amplitude-based conversion observable as the primary, run convergence and t_end-independence checks on the existing campaign findings, and finish a focused light-mediator sweep to sharpen the science finding. Total ~3 days agent work + ~2 hours HPC.

After that, the project is in a state where the flagship physics question is answered (no amplification at minimal / dark-photon / propagating couplings; documented limitation on energy-based observables for the affected sector) and the supervisor can decide whether to commission additional physics (non-minimal coupling, cubic PGT) or to stop.

### What if Reviews 2 and 3 are wrong (overcautious)?

The strongest review claims — fabricated authors, hallucinated quotes, parity-odd mismatch with CRZ 2024 — could be over-stated. Meta-Review K is checking these in parallel. Even granting full benefit-of-the-doubt:
- The metric h₄/h₇/h₉ subspace remains the documented blocker (the no-go theorems are robust independent of literature attributions).
- Path A's L_VT remains incompatible with TIDAL's JSON schema (this is purely structural, doesn't depend on literature claims).
- The amplitude-based observable remains the right operational primary regardless of how Path A and Path B are framed.

The next-step plan above does not depend on any disputed literature claim. It depends only on the structural project-relevance findings of Review 3 and on the existing campaign results in the user's memory files. Even a fully-confirmed Round 3 picture would not change the recommendation: ship the operational primary, defer the methods paper, and finish the convergence checks.

### Final sentence

**Stop investigating, start finishing.** The science question (does this class of PGT theories produce graviton–photon amplification?) has been answered ("no, modulo a documented limitation on energy-based observables") with substantial evidence from existing campaigns. Three days of operational work converts that into a defensible, documented, measurement-grade project state. Six weeks of additional theory-paper work does not.
