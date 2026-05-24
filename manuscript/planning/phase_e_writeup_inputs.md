# Phase E — Manuscript writeup inputs (inventory)

Direct feed to the Phase 2 drafting session (see `report_plan.md`). One row per
artefact this campaign will land into the MSci report. Filled in as Phase E
results arrive.

**Companion:** [docs/PHASE_E_TRACKER.md](../../docs/PHASE_E_TRACKER.md) is the
campaign tracker (status, jobids, verdicts). This file is the writeup-side index.

## §Introduction

| Status | Asset | File path | Words | Notes |
|--------|-------|-----------|-------|-------|
| none | Localised-geometry motivation paragraph | — | ~50–100 | Why we pivot from uniform-field; cite the t_end-conventional headline problem from Phase B. |

## §Methodology — geometry pivot

| Status | Asset | File path | Words | Notes |
|--------|-------|-----------|-------|-------|
| none | Dual-Gaussian schematic + wavepacket overlay figure | `manuscript/figures/phase_e_geometry_schematic.pdf` (TBD) | ~100 caption | Source TOML: `examples/gertsenshtein/theory_ungauged_e_dual_gaussian.toml` |
| none | Frequency-regime table (k·σ_B≈15 etc.) | (inline) | ~80 | Copied from `docs/V3_PHASE_E_DESIGN.md` |
| none | TT-IC + theory unconstrained pattern paragraph | — | ~80 | Cite issue #167; cite `theory_ungauged.toml` |
| none | Stability-diagnostics description | — | ~150 | 4 metrics from `tidal/measurement/_phase_e_transit.py` |

## §Results — per-theory headline numbers

| Status | Asset | File path | Words | Notes |
|--------|-------|-----------|-------|-------|
| none | Headline amplification/suppression table (one row per theory: log(Z)_amp, log(Z)_sup, BF, D_KL, A_max, verdict, stability) | `manuscript/tables/phase_e_headline.tex` (TBD) | ~80 caption + ~250 prose | Built from `hpc_results/<jobid>/phase_e/<theory>/inference.json` + `stability.json` |
| none | E.cal corner figure | `manuscript/figures/phase_e_cal_corner.pdf` (TBD) | ~50 caption | Positive control — verifies P > 1e-4 · Boccaletti |
| none | E.EH corner figure | `manuscript/figures/phase_e_eh_corner.pdf` (TBD) | ~50 caption | QED F⁴ cross-check |
| none | E.T1 corner figure | `manuscript/figures/phase_e_t1_corner.pdf` (TBD) | ~50 caption | DP-plasma |
| none | E.T2 corner figure | `manuscript/figures/phase_e_t2_corner.pdf` (TBD) | ~50 caption | Einstein-Cartan baseline |
| none | E.T4 corner figure | `manuscript/figures/phase_e_t4_corner.pdf` (TBD) | ~50 caption | Ricci-EM nonminimal δ₁ |
| none | E.T5 corner figure | `manuscript/figures/phase_e_t5_corner.pdf` (TBD) | ~50 caption | Bahamonde 5D |
| none | E.NP corner figure | `manuscript/figures/phase_e_np_corner.pdf` (TBD) | ~50 caption | Non-propagating-torsion variant (ξ=0) |
| none | E.T7s / E.T8s corner figures | `manuscript/figures/phase_e_t{7,8}s_corner.pdf` (TBD) | ~80 captions | R²-stripped complete-{even,odd} |

## §Results — atlas projections (only theories ≤5D)

| Status | Asset | File path | Words | Notes |
|--------|-------|-----------|-------|-------|
| none | E.EH atlas (2D) | `manuscript/figures/phase_e_eh_atlas.pdf` (TBD) | ~80 caption | Methodology shake-down |
| none | E.T2 atlas | `manuscript/figures/phase_e_t2_atlas.pdf` (TBD) | ~80 caption | |
| none | E.T1 atlas | `manuscript/figures/phase_e_t1_atlas.pdf` (TBD) | ~80 caption | |
| none | E.T5 atlas (5D, headline) | `manuscript/figures/phase_e_t5_atlas.pdf` (TBD) | ~100 caption + ~150 prose | Headline atlas figure |

## §Discussion — Phase B vs Phase E comparison

| Status | Asset | File path | Words | Notes |
|--------|-------|-----------|-------|-------|
| none | Phase B vs Phase E shift table (per-coupling D_KL, A_max shift, verdict-change rows) | `docs/PHASE_E_GEOMETRY.md` (TBD; built at Stage 7.1) | ~100 caption + ~300 prose | Decision rule: < 0.5 nats per coupling → Phase B headline OK; > 0.5 → Phase E becomes canonical |
| none | Finite-interaction-time stability narrative (catastrophic vs finite-growth distinction) | — | ~200 prose | Cite stability diagnostics in `_phase_e_transit.py`; reference per-theory `stability.json` verdicts |
| none | "Why localised geometry is canonical going forward" closing paragraph | — | ~100 prose | Headline number choice |

## §Appendix C — TIDAL software

| Status | Asset | File path | Words | Notes |
|--------|-------|-----------|-------|-------|
| none | Modal solver #367 fix subsection (Phase E enabler) | — | ~250 prose | Cite commits c2b4e00 (v0.42.0) + 7dceb0c (v0.42.1); cite the regression test |
| none | Cubed-sphere atlas methodology subsection | — | ~300 prose | Source: `tidal/inference/_sphere.py` + `_atlas.py`; cite `R_ATLAS=0.4` rationale from `_geometry.env` |
| none | Phase E stability diagnostics subsection | — | ~200 prose | 4 metrics + their use of the t_check_1/t_check_2 checkpoints; cite `_phase_e_transit.py` |

## Status legend

- **none**: not yet produced
- **draft**: produced; awaiting analyst sign-off
- **final**: signed off, in manuscript

## How to maintain this file

Toggle the `Status` column in place as artefacts land. Add notes inline if an
asset gets re-rendered or replaced. Coupled to the SAME commit that produces
the artefact — never a separate "writeup inputs update" commit.
