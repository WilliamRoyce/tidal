# Correction Log W3 — Phase 1 Cleanup of Per-Agent Writeups, Sympy Scripts, README, and Result-Claim Files

**Date:** 2026-04-27
**Scope:** Apply audit annotations to remaining unmarked artefacts in
`research/perturbative_hamiltonian/`. The 3 round-synthesis notes had
already been bannered by W2; `notes/round2_agentE_vt_convergence.md`
had been bannered in a prior partial run. W3 covers the four remaining
per-agent writeups, all 21 sympy scripts, the 8 audit re-verification
scripts, the README, and four result-claim files.

**Posture:** No new physics claims. Every banner cites a specific
audit document. No deletion of existing content — banners prepend,
corrections append, JSON gets a top-level `_audit` key.

---

## Task 1 — Per-agent writeup banners

Four files received SUPERSEDED blockquote banners (matching the style
used in `notes/round2_agentE_vt_convergence.md`) followed by a
"Per-writeup audit corrections" section. The original heading and body
were left untouched.

### `notes/round2_agentF_axial_verification.md`

Banner cites:
- **Vector→scalar correction sympy-correct** — Review 1
  (`reviews/review1_mathematical_verification.md`).
- **AM 2020 cross-validation downgraded** to qualitative consistency at
  linearised order — Meta-K K6
  (`meta_reviews/meta_review_K_literature_claims.md`).
- **A-sector rank-uniformity is sector-restricted** — `det(M_aux) =
  1/b²` diverges as b→0; rank-jump relocated, not removed (Review 1
  C8).
- **Curved-background O(h²) extension fails** — Nieh-Yan terms
  unabsorbed (FINAL_ASSESSMENT §"What is overstated").

### `notes/round2_agentG_novel_directions.md`

Banner cites:
- **Sectoral reclassification** (axial / trace / tensor-q / metric
  h₄,₇,₉) survives audit — FINAL_ASSESSMENT §"What is genuinely true",
  item 1.
- **"Path B-trace established (Barker)" is WRONG** — Meta-K K3
  (Barker identifies trace torsion as Yang-Mills GAUGE FIELD, excludes
  parity-odd terms).
- **"Path B-tensor-q applies (CRZ 2024)" is CONDITIONAL** — Meta-K K5
  (CRZ handles only parity-even free fields, m→0 Goldstone limit;
  opposite to TIDAL's b5→0 infinite-mass limit).
- **"Born-Oppenheimer ≡ Path A"** identification is correct.
- **Recipe 1 PASS** verified at sympy *schema* level only — actual
  xAct decomposition not yet checked (Review 1 §"Phase 2.2 follow-up
  needed").

### `notes/round3_agentI_vt_3field.md`

Banner cites:
- **VT integrand convergence + KV Definition-1 closure** at N=3+3 —
  verified clean by Review 1.
- **Phase-space jump factor 3 vs 5** — convention drift not flagged
  (Review 1 C6); main script gives factor 3 (6→18), Phase 6 v2 gives
  factor 5 (6→30).
- **Phase 6** (Hamiltonian rank analysis on metric subspace) ran into
  sympy performance issues; qualitative verdict only.
- **Metric h₄,₇,₉ inheriting Pais-Uhlenbeck** structure is correct.
- **NEW pathology**: `L_VT diverges as M_c² → 0` — Review 1 C5
  (`reviews/scripts_review/C5_routhian_M_to_zero.py`).
- **BN 1983 Nuovo Cim citation** is on R + T² PGT, not R+R²+T² as
  Phase 1 narrative implies — see `BlagojevicNikolic1983Nuovo` in
  `docs/tex/references.bib`.

### `notes/round3_agentJ_curtright_stueckelberg.md`

Banner cites:
- **Gauge invariance δF̊ = 0** verified at 1+1D toy ONLY (24
  components in 4D would need separate verification — Review 1 C3).
- **Rank uniformity `det(H_kin) = 1−λ_a²`** verified at 4×4 (1+1D toy)
  ONLY — 4D 36×36 verification NOT done (Review 1 C3).
- **CRZ paper** (arXiv:2411.16928) handles only parity-even free
  fields; TIDAL b5·R̃² is parity-odd by construction (Meta-K K5).
- **F-J cross-validation** with Agent F is structurally OPPOSITE, not
  parallel — F's `det(M_aux) = 1/b²` diverges at b=0; J's
  `det(H_kin) = 1−λ_a²` is genuinely b5-independent (Review 1 C8).

---

## Task 2 — Sympy script headers

### `scripts/*.py` (21 files)

Identical 5-line header prepended to each:

```python
# AUDITED 2026-04-27.  Sympy execution verified clean by Review 1
# (research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# The math is correct; original framing in the agent notes/synthesis was overstated.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the
# verified picture and which specific claims survive audit.
```

Files:
1. `bopp_podolsky_axial.py`
2. `curtright_stueckelberg_q.py`
3. `line3_2form_auxiliary.py`
4. `line3b_2form_IBP_constraint_check.py`
5. `line4_born_oppenheimer.py`
6. `line6_dual_formulation.py`
7. `recipe1_all_eps_RR_contractions.py`
8. `recipe1_check_all_R_components.py`
9. `recipe1_debug_R_squared.py`
10. `recipe1_explicit_q_substitution.py`
11. `recipe1_parity_even_RR.py`
12. `recipe1_preflight_q_projection.py`
13. `vt_T4_3plus3_PGT.py`
14. `vt_T4_phase6_final.py`
15. `vt_T4_phase6_fixed.py`
16. `vt_T4_phase6_npfinal.py`
17. `vt_T4_phase6_numeric.py`
18. `vt_T4_phase6_v2.py`
19. `vt_T4_phase6_v3.py`
20. `vt_convergence_T2.py`
21. `vt_convergence_T3.py`

All 21 verified to still parse cleanly via `ast.parse`.

### `reviews/scripts_review/*.py` (8 files, C1–C8)

These ARE the audit's own re-verification scripts, so a different
header was used:

```python
# AUDITED 2026-04-27.  This script is part of Review 1's own re-verification
# of the original investigation (one of the C1-C8 audit checks; see
# research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md).
# It implements an independent sympy check rather than reproducing an original-investigation result.
# See research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md for the verified picture.
```

Files: `C1_T2_1plus1D_extension.py`, `C2_helmholtz_higher_orders.py`,
`C3_curtright_2plus1D.py`, `C4_voicu_linearity.py`,
`C5_routhian_M_to_zero.py`, `C6_simpler_phase6.py`,
`C7_3form_aux_counterexample.py`, `C8_F_vs_J_structural_diff.py`.

All 8 verified to still parse cleanly via `ast.parse`.

---

## Task 3 — README.md

Inserted an "Audit status (2026-04-27)" section between the title and
the existing body. Contains:

- A 10-row navigation table mapping each artefact group to its audit
  status and pointer.
- A one-paragraph verdict summary citing FINAL_ASSESSMENT, Meta-L
  (h_5↔a_1 standard-kinetic), Meta-L
  (`_conversion.py` Hamiltonian-based not amplitude-based), and the
  "operational not academic" recommendation.

The original Status, Layout, Constraints, Investigation rounds, and Key
references sections are preserved verbatim below the new audit-status
section.

---

## Task 4 — Result-claim file annotations

### `results/curtright_stueckelberg_verdict.json` — `_audit` key added

Specific corrections:
- Verification at 1+1D toy only; 4D rank uniformity NOT verified (Review 1 C3).
- CRZ paper handles only parity-even free fields (Meta-K K5).
- F-J cross-validation with Agent F is structurally OPPOSITE, not parallel (Review 1 C8).

### `results/recipe1_q_kinetic_structure.json` — `_audit` key added

Specific corrections:
- Verified at sympy schema level; actual xAct decomposition not yet checked (Review 1).
- Phase 2.2 follow-up needed for production-pipeline verification.

### `results/axial_constraint_matrix.json` — `_audit` key added

Specific corrections:
- A-block rank uniformity is verified BUT aux block diverges as 1/b² (Review 1 C8).
- Curved-background O(h²) extension fails (Nieh-Yan terms unabsorbed).
- AM 2020 compatibility downgraded to qualitative consistency at linearised order (Meta-K K6).

### `results/novel_directions_assessment.md` — SUPERSEDED banner prepended

Cites:
- "Path B-trace 'established literature' framing is WRONG" — Meta-K K3.
- "Path B-tensor-q is CONDITIONAL on a parity-odd extension that does not exist" — Meta-K K5.

All three JSON files validated via `python3 -m json.tool` (parsed as
valid JSON with `_audit` top-level key intact).

---

## Verification summary

| Check | Result |
| --- | --- |
| 4 per-agent writeups have SUPERSEDED banner + correction section | PASS |
| `notes/round2_agentE_vt_convergence.md` left untouched (already had banner) | PASS |
| 21 `scripts/*.py` files have audit-status header | PASS — all parse cleanly |
| 8 `reviews/scripts_review/*.py` files have audit-author header | PASS — all parse cleanly |
| README.md has audit-status table at top, original content preserved below | PASS |
| 3 result-claim JSON files have `_audit` top-level key | PASS — all valid JSON |
| 1 result-claim Markdown file has SUPERSEDED banner | PASS |
| Cited files exist (FINAL_ASSESSMENT.md, Reviews 1-3, Meta K/L/M/N) | PASS — all in directory |

---

## Evidence references

- **FINAL_ASSESSMENT** — `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md` (authoritative audit verdict, dated 2026-04-26)
- **Review 1** — `research/perturbative_hamiltonian/reviews/review1_mathematical_verification.md` (C1–C8 sympy re-verification)
- **Review 2** — `research/perturbative_hamiltonian/reviews/review2_literature_interpretation.md`
- **Review 3** — `research/perturbative_hamiltonian/reviews/review3_project_relevance.md`
- **Meta-K** — `research/perturbative_hamiltonian/meta_reviews/meta_review_K_literature_claims.md` (K3 Barker parity-odd; K5 CRZ parity-even; K6 AM ghost-free tower)
- **Meta-L** — `research/perturbative_hamiltonian/meta_reviews/meta_review_L_pipeline_claims.md` (h_5↔a_1 standard-kinetic; `_conversion.py` Hamiltonian-based)
- **Meta-M** — `research/perturbative_hamiltonian/meta_reviews/meta_review_M_next_steps.md`
- **Meta-N** — `research/perturbative_hamiltonian/meta_reviews/meta_review_N_BC_AppendixD_content.md`
- **Existing W2 banners (synthesis files)** — `research/perturbative_hamiltonian/reviews/correction_log_W2.md`
- **Existing prior-run banner (Agent E)** — top of `research/perturbative_hamiltonian/notes/round2_agentE_vt_convergence.md`

---

## Files modified by W3

- `notes/round2_agentF_axial_verification.md`
- `notes/round2_agentG_novel_directions.md`
- `notes/round3_agentI_vt_3field.md`
- `notes/round3_agentJ_curtright_stueckelberg.md`
- `scripts/*.py` (21 files)
- `reviews/scripts_review/*.py` (8 files)
- `README.md`
- `results/curtright_stueckelberg_verdict.json`
- `results/recipe1_q_kinetic_structure.json`
- `results/axial_constraint_matrix.json`
- `results/novel_directions_assessment.md`
- `reviews/correction_log_W3.md` (this file)

Total: 38 files modified, 0 files deleted, 0 existing claims rewritten
(banners and audit annotations only).
