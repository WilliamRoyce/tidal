# Perturbative Hamiltonian Reduction for Constraint-Promotion PGT

## Audit status (2026-04-27)

A 10-agent investigation (Rounds 1-3) was followed by 3 critical reviews
(Reviews 1-3) and 4 meta-reviews (K, L, M, N). The audit corrected
several headline overstatements. **Use the table below to navigate
which artefacts are still authoritative.**

| Artefact group | Audit status | Where to look |
| --- | --- | --- |
| `notes/FINAL_ASSESSMENT.md` | **authoritative** | start here |
| `reviews/review{1,2,3}_*.md` | **authoritative** | mathematical / literature / project-relevance audits |
| `meta_reviews/meta_review_{K,L,M,N}_*.md` | **authoritative** | corrects review over-corrections, plans next steps |
| `notes/round{1,2,3}_synthesis.md` | superseded (banner present) | historical only |
| `notes/round2_agent{E,F,G}_*.md`, `notes/round3_agent{I,J}_*.md` | superseded (banners present) | per-agent writeups; see file-level audit corrections |
| `scripts/*.py` (21 files) | sympy execution verified clean (Review 1) | math correct, original framing overstated |
| `reviews/scripts_review/*.py` (8 files, C1–C8) | audit's own re-verification scripts | independent checks, not original investigation |
| `results/*.json`, `results/*.md` | annotated (`_audit` key or banner) | per-file corrections noted |
| `docs/tex/perturbative_reduction_constraint_barrier.tex` | corrected by W1/W2 (see `reviews/correction_log_W{1,2}.md`) | canonical writeup |
| GitHub issue #321 | updated to reflect verified picture | tracking |

**One-paragraph verdict.** The metric h₄/h₇/h₉ Pais-Uhlenbeck subspace
is genuinely blocked for local first-order auxiliary lifts (3
convergent no-gos), but it does NOT gate TIDAL's headline observable
(`h_5 ↔ a_1` is standard-kinetic per Meta-L). The
`tidal/measurement/_conversion.py` pipeline is Hamiltonian-based
(NOT amplitude-based as Review 3 claimed) and produces measurement-grade
P(t). The science question has been answered NULL by existing campaigns
across multiple regimes. The recommended next move is operational
(documentation corrections + validation sweeps + light-mediator sweep),
NOT the academic publication trajectory the original synthesis
recommended. See `notes/FINAL_ASSESSMENT.md` for full detail.

---

**Research investigation**, parallel to `research/general_quadratic_lagrangian/`.
This directory holds the in-progress analysis of the
**constraint-promotion barrier** in TIDAL's perturbative-Hamiltonian
reduction, specifically for theories of the form

```
L = (1/κ²)R̃ + α₁I₁ + α₂I₂ + α₃I₃ + b5·R̃² - ¼F²
```

where R̃ is the Riemann-Cartan Ricci scalar and the small parameter b5
promotes algebraic-constraint fields (h_4, h_7, h_9) to 4th-order
dynamical at b5 ≠ 0.

## Status (2026-04-26)

A four-agent deep theoretical investigation (see
`docs/tex/perturbative_reduction_constraint_barrier.tex` and
`/home/vscode/.claude/plans/review-the-work-done-shiny-mitten.md`)
has produced:

### Closed (genuinely blocked, do not pursue)

- ✗ JLM-style algebraic substitution
- ✗ Hamiltonian reconstruction from EOM via antiderivative (Helmholtz of `q̈`-dependent force forbidden)
- ✗ FJ rank-jump perturbation theory (Schur termination requires non-degeneracy)
- ✗ **Stückelberg (irreducible AND reducible) for generic PGT** — three convergent arguments confirm `det(M) ∝ b5^N` regardless of reducibility structure (Agent A + Agent D)

### Open (positively viable)

- ✅ **Path A — Krupka-Voicu / Vainberg-Tonti for the FULL PGT b5·R̃²**
  Helmholtz residue δE = 0 generically when M and N matrices come from the same Lagrangian (Agent B). Remaining preflight: VT integral convergence (Voicu 2020 4D-GB cautionary).

- ✅ **Path B — Sector-by-sector partial extensions** (Agent C)
  - Axial torsion sector: explicit Bopp-Podolsky single-auxiliary Lagrangian, clean b5→0 limit, no rank-jump
  - Trace torsion sector: works via Barker et al. 2024 conformal embedding (parity-even part)
  - Tensor torsion sector: blocked

## Layout

- `scripts/` — sympy / Python investigation code (executable, no TIDAL pipeline modification)
- `notes/` — derivations, analysis writeups (Markdown / TeX)
- `results/` — JSON or text outputs from script runs

## Constraints

- **READ-ONLY against the TIDAL pipeline.** Do not modify
  `tidal/wolfram/PerturbativeReduction.wl`, `tidal/solver/perturbative_driver.py`,
  `tidal/wolfram/ExportJSON.wl`, or any test files.  This is purely
  investigative work.
- **All outputs land in this directory.**  Sympy scripts are standalone;
  notes are Markdown unless they explicitly need TeX features.
- Cross-reference the docs (`docs/tex/perturbative_reduction_constraint_barrier.tex`)
  but do not modify them yet — landing investigation results into the main
  doc is a separate step done after results stabilize.

## Investigation rounds

### Round 1 (complete) — see `notes/round1_synthesis.md`

Four parallel agents covered:
- Agent A: Lyakhovich Stückelberg recipe applied to PGT toy
- Agent B: Helmholtz residue computation (Vainberg-Tonti preflight)
- Agent C: Hinterbichler-Saravani extension to PGT torsion
- Agent D: Lyakhovich-Sharapov reducible-generator analysis

### Round 2 (active 2026-04-26)

- Agent E: VT integral convergence test on minimal PGT toy
- Agent F: Independent verification of axial-torsion Bopp-Podolsky lift
- Agent G: Novel alternative directions (BV-BFV homological, Kontsevich, alternative auxiliary fields for tensor sector)

## Key references

| arXiv ID | Authors | Why |
|----------|---------|-----|
| 1804.05556 | Blagojević-Cvetković 2018 | Modern "critical parameters" classification; Appendix D no-limit quote |
| 1406.6646 | Krupka-Voicu 2015 | Vainberg-Tonti canonical variational completion |
| 2009.05459 | Voicu 2020 | 4D-Gauss-Bonnet VT divergence cautionary |
| 1508.02401 | Hinterbichler-Saravani 2015 | Stückelberg quadratic curvature (metric-only) |
| 2009.11739 | Aoki-Mukohyama 2020 | PGT-bigravity equivalence (3D works, 4D needs infinite tower) |
| 2102.10579 | Lyakhovich 2021 | General Stückelberg recipe |
| 2106.09355 | Abakumova-Lyakhovich 2021 | Reducible Stückelberg |
| 2406.12826 | Barker et al. 2024 "Every PGT is conformal" | Trace torsion as Goldstone of broken Weyl |
| 2501.00133 | (Bopp-Podolsky reduction) | Single-auxiliary lift template |
| 2206.00658 | Barker 2022 | HiGGS validation tool |
| 2210.15980 | Chiou-Geiller-Wang 2022 | Tetrad-gravity DOF discontinuity (twin) |
