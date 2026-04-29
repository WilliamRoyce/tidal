# Phase 2 Synthesis — Lead Investigation Verdict

**Date**: 2026-04-27
**Scope**: 5 focused single-agent investigations of audit-identified leads
**Verdict**: **(c) Definitive gap** — no constructive published bridge for TIDAL's
b5=0 critical surface in PGT b5·R̃², but phenomenon IS documented in adjacent
frameworks.

---

## Outcome by lead

| Lead | Local evidence | Verdict |
|------|---------------|---------|
| **2.1** Glavan-Zlosnik-Lin 2024 (arXiv:2311.17459) | `notes/lead_glavan_zlosnik_lin.md` | **(c) NOT TRANSFERABLE** — three independent obstructions: symmetry mismatch (GZL uses projective + Weyl, TIDAL doesn't); higher-derivative-trick singularity (GZL requires `f''(φ) ≠ 0`, vanishes at b5=0); geometry mismatch (metric-affine vs vielbein PGT, different gauge algebras) |
| **2.2** Recipe 1 vs xAct production pipeline | `notes/lead_recipe1_xact.md` | **(a) PASS at production**, with three caveats: C1 operator-mismatch (Agent H verified parity-odd Pontryagin density; TIDAL uses parity-even RicciScalarCDT²); C2 DOF-count discrepancy (16 off-shell vs 3 post-plane-wave); C3 hamiltonian_terms filter (must read `equations[]`, not the filtered-of-torsion `hamiltonian_terms`) |
| **2.3** Lyakhovich-Sharapov forward citations 2022-2026 | `notes/lead_lyakhovich_forward_citations.md` | **(b) NO NEW LEAD** — Agent D's structural no-go is robust under 2022-2026 literature; closest follow-up (arXiv:2303.02616, multi-layer reducible Stückelberg) strengthens rather than overturns the no-go; **Aashish-Saif arXiv:2601.22007 was misattributed by Round 2 Agent G — does not cite Lyakhovich at all** |
| **2.4** BN 1983 direct read | `notes/lead_bn1983_direct_read.md` | **(b)-with-amendment** — BN 1983 PRD 28:2455 (companion to the Nuovo Cim paper) DOES cover R+T²+R² PGT including parity-odd `(εR)²`; BUT explicitly defers the massless-tordion case (TIDAL's b5→0 case); sharpens "25-year unsolved" → **"40-year acknowledged-as-open"** (BN 1983 PRD → YNN 1999/2002 → BC 2018) |
| **2.5** Beltrán Jiménez/Heisenberg/Koivisto metric-affine corpus | `notes/lead_bjhk_metric_affine.md` | **(b) NO TRANSFER** — three papers all scope-restricted (parity-even IDG; FLRW-only teleparallel; f(Q) "wrong direction"); **misattribution: arXiv:1911.08846 is NOT a BJHK paper** (Cruz-Dombriz/Maldonado Torralba/Mazumdar) |

---

## Two unexpected positive findings (audit corrections)

These tighten the framing rather than open a new path:

### F1 — Recent published analogues in metric-only quadratic gravity (Phase 2.1)

The constraint-promotion phenomenon **is documented in 2024-2025 peer-reviewed
literature** for metric-only quadratic gravity, by credible authors:

- **Barker-Glavan 2025 (arXiv:2510.08201)** "Spectrum of pure R² gravity":
  > *"the mechanism for this phenomenon is a change in the nature of the
  > constraints upon linearisation: ten second-class constraints of the full
  > theory become first-class. ... Such backgrounds are surfaces of strong
  > coupling in field space, where the dynamics of perturbations becomes
  > nonperturbative."*
- **Bellorin 2025 (arXiv:2506.07305)** "Hamiltonian equations of motion of
  quadratic gravity":
  > *"One cannot obtain the ADM Hamiltonian of general relativity as a smooth
  > limit of the case we have considered."*
- **Karananas 2024 (arXiv:2407.09598, 2408.16818)** "Particle content of
  (scalar curvature)² gravities": pure R² is "infinitely strongly-coupled at
  Minkowski"; 2408.16818 explicitly covers Einstein-Cartan quadratic
  parity-even sector.

**Implication**: TIDAL's contribution is identifying the same barrier in the
parity-odd PGT-torsion variant where it has not been published.  The "barrier
is real" narrative is now anchored in current peer-reviewed literature, not a
tenuous historical claim.

### F2 — BN 1983 PRD is the right historical anchor (Phase 2.4)

The audit's revised framing ("BN 1983 covers R+T² only — not TIDAL's class")
was wrong on scope.  BN 1983 PRD 28:2455 (the companion to the Nuovo Cim
paper) covers **R+T²+R² PGT** including a parity-odd `(εR)²` invariant, and
explicitly defers the case TIDAL hits.

The single most important verbatim quote (`/tmp/BN1983_PRD.txt:1052-1062`):

> *"In the case corresponding to the massless tordions (with a new kind of
> gauge symmetry in the theory) more detailed considerations are necessary
> for a complete understanding of the situation."*

This is a 1983 published acknowledgement that the specific problem TIDAL
faces (massless-rank-jump in higher-derivative PGT) was an open research
question 40 years ago.  Stronger framing than the 25-year-from-BC2018 chain.

---

## Audit corrections to apply to canonical documentation

### Correction A — replace BN 1983 historical anchor

**Where**: `docs/tex/perturbative_reduction_constraint_barrier.tex`,
§Historical context.

**Change**: cite BN 1983 PRD (`BlagojevicNikolic1983`) as the historical anchor
that covers TIDAL's Lagrangian class; cite BN 1983 Nuovo Cim
(`BlagojevicNikolic1983Nuovo`) only as the R+T² stepping-stone precursor.
Add the verbatim 1983 deferral quote.

**Effect**: sharpens "25-year unsolved" → "40-year acknowledged-as-open".
Stronger evidence that the specific problem is recognized as open.

### Correction B — add direction-inversion framing

**Where**: same TeX section as Correction A.

**Change**: add an explicit paragraph noting that BN 1983 / BC 2018 use the
direction `c_n → 0 ⇒ if-constraint becomes a primary constraint`, while
TIDAL faces the inverse direction `b5 ≠ 0 ⇒ Lagrange-multiplier becomes a
4th-order propagating field`.  Both are rank-jump phenomena; opposite sides
of the critical surface.

**Effect**: explains why BN 1983's case-by-case dispatch (Table I) and BC
2018's F-matrix construction don't directly transfer — TIDAL needs the
inverse direction not handled there.

### Correction C — cite recent published analogues

**Where**: §Future work or §Methods that fail (decide based on TeX
flow).

**Change**: add a paragraph noting that the constraint-promotion phenomenon is
recognised in 2024-2025 peer-reviewed literature for metric-only quadratic
gravity (Barker-Glavan 2510.08201; Bellorin 2506.07305; Karananas
2408.16818).  Each surfaces the same barrier ("strong coupling in field
space"; "no smooth GR limit"; "infinitely strongly-coupled at Minkowski").
TIDAL identifies the parity-odd PGT-torsion variant.

**Effect**: anchors the "barrier is real" claim in current literature.
Supports Publication C if pursued.

### Correction D — remove/restrict Aashish-Saif citation

**Where**: any doc citing arXiv:2601.22007 in the constraint-promotion
context.

**Change**: Phase 2.3 verified arXiv:2601.22007 does NOT cite Lyakhovich and
addresses a structurally different singularity category (field-configuration
vacuum, not parameter critical surface).  Remove from Stückelberg-lead
sections; restrict any citation to its actual scope (Lorentz-violating
antisymmetric-tensor Stückelberg).

### Correction E — fix BJHK metric-affine misattribution

**Where**: any doc citing arXiv:1911.08846 as a BJHK paper.

**Change**: Phase 2.5 verified the actual authors are Cruz-Dombriz/Maldonado
Torralba/Mazumdar, framework is infinite-derivative gravity with non-symmetric
connection, scope is parity-even only.  Update author attribution and scope
notes.  Beltrán Jiménez-Maldonado Torralba 2019 (arXiv:1910.07506) is the
genuine BJHK-adjacent reference if cited.

### Correction F — Recipe 1 operator-mismatch caveat

**Where**: §Open research directions or wherever Recipe 1 PASS is cited.

**Change**: note that Agent H's PASS verdict was verified for the parity-odd
Pontryagin density `ε^{abef} R̃^{abcd} R̃_{cd}^{ef}`, while TIDAL's actual
operator is the parity-even `RicciScalarCDT[]^2` per
`examples/torsion_gertsenshtein/theory.toml:118`.  The conclusion (no
`(∂²q)²` at linear order) carries qualitatively because each Ricci-type R̃
factor is still one `∂K` at linear order, but the verbatim sympy verification
does not cover TIDAL's operator.  Recommended: ~30-min extension of
`recipe1_preflight_q_projection.py` adding a part (D) for
`(g^{μν} R̃_{μν})²`.

### Correction G — note GZL's f''(φ)=0 obstruction

**Where**: §Future work or §Open directions where Path A / Krupka-Voicu is
discussed.

**Change**: note that GZL 2024's standard scalar-auxiliary lift requires
`f''(φ) ≠ 0`; for `f = R + b5·R²` this is `f'' = 2b5`, which vanishes at b5=0.
GZL's own construction is non-invertible at the constraint-promotion limit.
This reinforces (rather than weakens) the audit's verdict that no published
technique bridges the b5=0 surface.

---

## Phase 2 deliverables (for verification)

5 lead writeups in `notes/`:
- ✅ `lead_glavan_zlosnik_lin.md` (Phase 2.1)
- ✅ `lead_recipe1_xact.md` (Phase 2.2)
- ✅ `lead_lyakhovich_forward_citations.md` (Phase 2.3)
- ✅ `lead_bn1983_direct_read.md` (Phase 2.4)
- ✅ `lead_bjhk_metric_affine.md` (Phase 2.5)
- ✅ `phase2_synthesis.md` (this file)

4 new arXiv papers downloaded to `literature/` (from Phase 2.1 forward citations):
- ✅ `2510.08201` Barker-Glavan 2025 (R² spectrum)
- ✅ `2506.07305` Bellorin 2025 (quadratic gravity Hamiltonian)
- ✅ `2407.09598` Karananas 2024 (R² particle content)
- ✅ `2408.16818` Karananas 2024 (metric-affine R²)

`MANUAL_RETRIEVAL_NEEDED.md` updated by Phase 2.3 + Phase 2.5 with
non-load-bearing forward-citation papers.

---

## Phase 2 decision-point outcome

Per the plan's three-outcome decision tree:

(a) **Lead found** — N/A.
(b) **No lead found, close** — N/A.
(c) **Definitive gap** — **YES, this is the verdict**.

> The published literature does NOT contain a bridge for TIDAL's specific
> case (PGT torsion + b5·R̃² + parity-odd + 4D), but the phenomenon IS
> published in adjacent frameworks (metric-only quadratic gravity per Phase
> 2.1; PGT 2nd-order per BN 1983 PRD with explicit deferral of TIDAL's
> case per Phase 2.4).
>
> The metric h₄/h₇/h₉ Pais-Uhlenbeck subspace remains genuinely blocked for
> local first-order auxiliary lifts (3 convergent no-gos from Round 1+2; not
> overturned by Round 3 Phase 2 leads).
>
> TIDAL's headline observable (h_5 ↔ a_1) is unaffected (standard-kinetic);
> `_conversion.py` is Hamiltonian-based and produces measurement-grade P(t).

---

## Recommended next-step (per Phase 2 decision-point (c))

The plan's instruction for verdict (c):

> *"Document the strengthened verdict, then proceed to Phase 3 (operational
> closure) without further theoretical work."*

**However, the user explicitly directed earlier (2026-04-27): "do not do hpc
sweeps and investigations like this, we havent resolved the blocking
issues".**  Phase 3 (operational closure) is therefore deferred until user
explicitly authorises it.

The **non-deferred** next step is to apply Corrections A-G above to the
canonical TeX writeup, the references.bib, and `FINAL_ASSESSMENT.md`.  This is
documentation work, not operational work — it can be done at user discretion.

If at any future point the user does authorise Phase 3, the existing
sub-steps (3.1 t_end-independence + grid convergence + Boccaletti
cross-check; 3.2 light-mediator regime sweep on h_5↔a_1; 3.3 doc closure +
version bump + issue close) remain ready to run.

---

## Cross-references

- `notes/FINAL_ASSESSMENT.md` — pre-Phase-2 audit picture
- `reviews/review{1,2,3}_*.md` — original audit
- `meta_reviews/meta_review_{K,L,M,N}_*.md` — audit-of-audit
- `MANUAL_RETRIEVAL_NEEDED.md` — local-cache status
- `docs/tex/perturbative_reduction_constraint_barrier.tex` — canonical writeup
  (corrections A-G needed)
- `docs/tex/references.bib` — bib (entries needed for Barker-Glavan, Bellorin,
  Karananas papers + 1910.07506 + 2308.02250)
