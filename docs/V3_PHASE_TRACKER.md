# v3 Phase Tracker

**Created:** 2026-05-10
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md)
**Purpose:** running checklist of v3 phases — updated as each task converges (mirrors how Phase 6.J/K/L/M were tracked in CAMPAIGN.md).

## Phase A — Soft-penalty refactor + compactified-prior + de-pruning (in progress)

| Step | Description | Status | Issue / Notes |
|---|---|---|---|
| A.0a | `docs/V3_ARCHITECTURE.md` written | ✅ done | 2026-05-10 |
| A.0b | `docs/V3_PHASE_TRACKER.md` written | ✅ done | (this file) |
| A.0c | `docs/V3_PHASE_E_DESIGN.md` written | ✅ done | Deferred-but-documented |
| A.0d | `CAMPAIGN.md` updated with v3 architecture pointer | ✅ done | |
| A.0e | GitHub issues created for v3-A/A-γ/B/D/E/C tracking | ✅ done | #345–#355 |
| A.0f | Obsolete HPC jobs (28982006, 28985879) cancelled | ⏳ in progress | |
| A.0g | Phase A.0 persistence committed | ⏳ pending | |
| A.1 | Soft-penalty refactor in `tidal/inference/_likelihood.py` | ⏳ pending | [#345](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/345) — lines 499–511, 581–590, 554, 557, 573 |
| A.2 | Lagrangian de-pruning audit + per-param `arctan_uniform` scripts | ⏳ pending | [#346](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/346) — 12 paired scripts |
| A.3 | Corner-plot upper-triangle removal in `tidal/inference/_visualize.py` | ⏳ pending | [#347](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/347) |
| A.4 | Soft-penalty tests (`test_likelihood_*.py`) | ⏳ pending | [#348](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/348) — 3 new test modules |
| A.5 | D1 v1 chain replay sanity check on `hpc_results/28520217/` | ⏳ pending | [#349](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/349) |
| A.6 | Commit + version bump | ⏳ pending | |

## Phase A-γ — γ_conversion (deferred to its own session)

| Step | Description | Status |
|---|---|---|
| A-γ.1 | Refactored γ_conversion with multi-t_test sampling + log-zero clamp in `_stability.py` | ⏳ deferred |
| A-γ.2 | `gamma_conversion` metric in `parse_likelihood()` (opt-in) | ⏳ deferred |
| A-γ.3 | 10-point D1 P_max correlation validation | ⏳ deferred |
| A-γ.4 | Decision: ship γ_conversion for amp campaigns or keep as metadata | ⏳ deferred |

GH issues `[v3-A-γ] γ_conversion: refactored definition + multi-t_test + log-zero clamp` and `[v3-A-γ] γ_conversion: P_max correlation validation` track the work.

## Phase B — Campaign re-runs

| Step | Description | Status |
|---|---|---|
| B.0 | Decision: per-param `arctan_uniform` vs cubed-sphere joint prior for publication numbers (small-scale comparison) | ⏳ pending (gated on parallel session landing joint prior) |
| B.1 | D1 amp v3 smoke (mid-res INTR) | ⏳ pending |
| B.2 | D1 amp v3 publication run | ⏳ pending |
| B.3 | D1 sup v3 paired with B.2 | ⏳ pending |
| B.4 | Stage A v3 paired re-run | ⏳ pending |
| B.5 | D2.0–D2.3 v3 paired re-runs | ⏳ pending |
| B.6 | Bounds-dependence cross-check at L=75 under v3 | ⏳ pending |

## Phase C — Cubed-sphere coupling-space chart (HANDLED BY PARALLEL SESSION)

See `docs/V3_PHASE_C_REFERENCE.md` once written by the parallel session. This phase is shipped externally to the v3 plan owned here.

Open follow-up question for supervisor: coupling grouping (monolithic / per-Lagrangian-symmetry-class / per-SPO-sector / per-parity / other). GH issue `[v3-C] Coupling-grouping question for supervisor follow-up` tracks.

## Phase D — Manuscript implications

| Step | Description | Status |
|---|---|---|
| D.1 | `manuscript/sections/computational_approach.tex` — add §"Tachyon-permissive inference architecture" | ⏳ pending (after Phase B) |
| D.2 | Update results sections with v3 numbers; comparison table v2 vs v3 | ⏳ pending (after Phase B) |
| D.3 | `docs/PHASE_6_COMPARISON.md` — append "v3 architecture" section | ⏳ pending (after Phase B) |

## Phase E — Localised geometry pivot (DEFERRED, gated on Phase B convergence)

See [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md). Replaces plane-wave + uniform B₀ with Gaussian wavepacket + Gaussian B-field profile; finite interaction time bounds P_max physically.

## Sign-off log

(Append entries as phases converge — date, step ID, summary, commit hash.)

- 2026-05-10 — A.0a — `docs/V3_ARCHITECTURE.md` written; canonical architecture reference live.
- 2026-05-10 — A.0b — `docs/V3_PHASE_TRACKER.md` written (this file).
