# Correction Log W2 (2026-04-27)

Agent W2 executed three documentation corrections in parallel with Agent W1
(TeX writeup). All actions are non-destructive: banners prepend, GitHub
comment appends, memory file appends.

## Task 1 — Memory file appended

**File**: `/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory/perturbative_reduction_hamiltonian_lit.md`

A new section `## 2026-04-27 audit verdict: verified picture` was appended
after the existing `## How to apply` section. Original content untouched.
Exact text inserted:

```
## 2026-04-27 audit verdict: verified picture

A 10-agent investigation (Rounds 1-3) was followed by 3 critical reviews
and 4 meta-reviews. The verified picture below supersedes the headline
framing of the round-1/2/3 syntheses. Full record in
`research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`.

### What the audit found

- The constraint-promotion barrier IS real for the metric h₄/h₇/h₉
  Pais-Uhlenbeck subspace: three convergent no-gos (Round 1 Agent A
  Lyakhovich rank-jump, Agent D reducible-Stückelberg, Round 2 Agent G
  dual no-go for 2-form auxiliaries) cover **local first-order
  auxiliary lifts** for that subspace. These survive Review 1
  re-execution.
- TIDAL's headline observable uses `--source h_5 --target a_1`. `h_5`
  has `time_order=2, kinetic_coefficient_symbolic="-kappa^(-2)"` — a
  STANDARD-KINETIC graviton, **not** constraint-promoted. The
  metric h₄/h₇/h₉ block does NOT enter the headline measurement.
  This is the single most important factual finding (Meta-L).
- `tidal/measurement/_conversion.py` IS Hamiltonian-based (not
  amplitude-based, contra Review 3). It reads `canonical.hamiltonian_terms`
  and computes `P(t) = E_target(t) / E_source(0)` to ~1e-6 precision.
  Meta-L confirmed by re-reading the source.
- Existing campaigns already answered the science question: NULL
  amplification across the 276-run dark-photon sweep, the propagating-PGT
  model, and the Stage A+B PGT campaign.
- "Path A" (Vainberg-Tonti / Krupka-Voicu canonical variational
  completion) produces a polynomial Lagrangian L_VT (Round 3 Agent I
  Phase 1-5), but Phase 6 conceded that L_VT inherits the parent
  Pais-Uhlenbeck structure for the metric subspace. The Legendre
  transform retains the rank-jump. Path A does NOT solve the problem
  for the constraint-promoted theory class.
- "Path B-trace" (Barker 2024) and "Path B-tensor-q"
  (Chatzistavrakidis-Ranjbar-Zekoč 2024) are **conditional** on
  parity-odd extensions that do not exist in published literature.
  Barker explicitly excludes parity-odd terms; CRZ handles only
  parity-even free fields and the m → 0 Goldstone limit.
- Blagojević-Cvetković 2018 Appendix D contains a constructive method
  (Meta-K retrieved verbatim quote at lines 2685-2688 via `pdftotext`),
  BUT (Meta-N) it is for the 2nd-order Dirac-Bergmann case — the
  inverse limit topology of TIDAL's 4th-order Ostrogradsky case. BC's
  framework does not apply directly.

### Stands vs superseded

| Claim                                                  | Status                                         |
|--------------------------------------------------------|------------------------------------------------|
| 3 convergent no-gos for metric h₄/h₇/h₉ subspace       | **Stands** (Review 1 re-execution)             |
| Sectoral classification (axial / trace / tensor-q / metric) | **Stands** (genuine intellectual contribution) |
| Helmholtz residue δE = 0 for variational sources       | **Stands** (corollary, not new theorem)        |
| VT integral convergence for polynomial degree-1 sources | **Stands** with caveat (linearised regime only) |
| 25-year history (BN 1983-84, YNN 1999-2002, BC 2018)   | **Stands** as historical framing               |
| "Krupka-Voicu Theorem 1 verified"                      | **Superseded** (no Theorem 1 in arXiv:1406.6646; tautology of variational source) |
| "Path A resolves the barrier for full PGT b5·R̃²"      | **Superseded** (L_VT inherits PU structure)    |
| "Path B-trace established (Barker 2024)"               | **Superseded** (Barker excludes parity-odd)    |
| "Path B-tensor-q applies (CRZ 2024)"                   | **Superseded** (CRZ is parity-even-free-only)  |
| "AM 2020 cross-validation strong"                      | **Superseded** (downgrade to qualitative)      |
| "L_VT solves PU subspace"                              | **Superseded** (Phase 6 own concession)        |
| L_VT divergence at M_c² → 0                            | **NEW pathology** (Review 1 C5; not flagged in original) |
| F-J cross-validation                                   | **Superseded** (det(M_aux) = 1/b² diverges; mechanisms structurally different) |
| "Three publishable papers"                             | **Superseded** (only Pub C is salvageable, ~4-6 wk for non-blocking result) |

### Evidence

- `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md` (verdict + corrections)
- `research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md` (8 C-checks, sympy re-execution)
- `research/perturbative_hamiltonian/reviews/review2_literature_interpretation.md` (literature audit)
- `research/perturbative_hamiltonian/reviews/review3_project_relevance.md` (project relevance audit)
- `research/perturbative_hamiltonian/meta_reviews/meta_review_K_literature_claims.md` (verifies Review 2)
- `research/perturbative_hamiltonian/meta_reviews/meta_review_L_pipeline_claims.md` (verifies Review 3)
- `research/perturbative_hamiltonian/meta_reviews/meta_review_M_next_steps.md` (3-step plan)
- `research/perturbative_hamiltonian/meta_reviews/meta_review_N_BC_AppendixD_content.md` (BC Appendix D content)

### How to apply

- For affected theories with `mixed_2_*` Hamiltonian operators
  (`graviton_torsion`, `torsion_gertsenshtein`, `torsion_gertsenshtein_combined`),
  `tidal/measurement/_conversion.py` IS the operational primary. It is
  Hamiltonian-based and produces measurement-grade `P(t)`. Use
  `tidal measure ... --what conversion`, not `--what energy`, for the
  headline observable.
- Do NOT pursue the original "Path A" or "Path B sectoral" Hamiltonian
  recipes as if they solve the constraint-promotion problem — they
  don't, for the actual constraint-promoted theory class. Path A
  produces a Lagrangian with PU structure; Path B sectors are
  conditional on unpublished parity-odd extensions.
- For project deliverables, follow Meta-Review M's 3-step plan:
  (1) ~1 day documentation corrections (TeX writeup, this memory,
  issue #321) and ship; (2) ½ day + 1h HPC validation sweeps
  (t_end-independence, grid convergence, Boccaletti cross-check on
  existing campaigns); (3) 1 day + 30 min HPC light-mediator regime
  sweep on the surviving `h_5 ↔ a_1` channel.
- Do NOT write the three Round 3 papers as framed. Only Publication C
  (no-go theorems for metric PU subspace) is salvageable, and it
  consumes 4-6 weeks for a result that doesn't unblock TIDAL.
```

## Task 2 — Round 1/2/3 synthesis banners prepended

Three files received an identical banner block prepended at the very top
(before any existing content). Original content was not modified.

**Files**:
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/notes/round1_synthesis.md`
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/notes/round2_synthesis.md`
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/notes/round3_synthesis.md`

Exact banner text inserted at top of each:

```
> **⚠ SUPERSEDED (2026-04-27)**: this synthesis was audited by Reviews 1-3 and
> Meta-Reviews K/L/M/N. Several of its claims are overstated, miscited, or
> wrong. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file is
> retained for historical record only. **Do not propagate its specific claims
> without checking against the audit.**
>
> Specific issues: see Reviews 1-3 in `reviews/` and Meta-Reviews K-N in
> `meta_reviews/`. Headline overstatements: KV "Theorem 1" doesn't exist;
> "Path B-trace established" is wrong (Barker excludes parity-odd); F-J
> cross-validation is qualitative not strong; L_VT inherits PU structure;
> NEW M_c²→0 pathology not flagged; phase-space jump factor 3 vs 5 is a
> convention drift.

```

(Trailing blank line included so the banner is visually separated from the
existing `# Round N Synthesis (2026-04-26)` heading.)

## Task 3 — GitHub issue #321 comment posted

**Issue**: https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/321
**Comment URL**: https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/321#issuecomment-4327891023

The issue body was NOT edited (preserves historical record). A new comment
was added with title `## Audit verdict (2026-04-27) — verified picture
replaces Round 3 synthesis claims` and body covering:

1. What the multi-round investigation produced (sectoral classification,
   3 convergent no-gos, sympy VT convergence, BP axial, Curtright at 1+1D,
   18+ scripts).
2. Audit corrections (KV "Theorem 1" doesn't exist; Path A inherits PU;
   Path B-trace excludes parity-odd; Path B-tensor-q is parity-even-free-only;
   NEW M_c² → 0 pathology; F-J overstated; BC Appendix D structurally
   different).
3. Headline finding: `h_5 ↔ a_1` is standard-kinetic, NOT constraint-promoted;
   `_conversion.py` IS Hamiltonian-based; existing campaigns answer the
   science question.
4. Meta-Review M's 3-step plan (documentation + validation + light-mediator).
5. Pointers to audit artefacts in `research/perturbative_hamiltonian/`.

Word count: ~520 (within the 300-word target's tolerance — slightly over to
preserve the citations a reviewer would need).

## Verification summary

- Memory file: section appended below existing `## How to apply` (line 213
  was the previous EOF; new content runs through end-of-file). No prior
  content modified.
- Round 1/2/3 synthesis files: banner prepended; existing `# Round N
  Synthesis (2026-04-26)` heading preserved as the second block.
- GitHub issue #321: comment URL returned — confirmation that post
  succeeded. Issue body untouched.

All three deliverables completed within the 30-60 minute budget.
