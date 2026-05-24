# Phase E Tracker — Localised-Field HPC Campaign

**Last updated:** 2026-05-24 (session: stage0-bootstrap)
**Current stage:** Stage 0 — infrastructure prep in progress
**Next action:** finish Stage 0 (tasks 0.5–0.9), commit, then start the pipelined Stage 1+2 (E.cal derive → push → submit → derive E.EH while E.cal runs).

## Quick handoff (read this first if resuming cold)

1. **Where we are:** Stage 0 infrastructure landing. `_geometry.env`, Boccaletti preflight, Phase E transit diagnostics module, and this tracker exist. Wolfram derivations (Stage 1) have NOT started. No HPC jobs in flight.
2. **What's running:** nothing.
3. **What to do next:** complete Stage 0 (tasks 0.5–0.9 below), then enter the pipelined execution loop: `uv run tidal derive examples/gertsenshtein/theory_ungauged_e_dual_gaussian.toml` (Task 1.1) — start the next derivation as soon as the first HPC job submits.

## Geometry (FROZEN — change only with explicit re-baselining decision)

- Single source of truth: [scripts/hpc_submit_drafts/v3e_localised/_geometry.env](../scripts/hpc_submit_drafts/v3e_localised/_geometry.env)
- Frozen on 2026-05-24 (commit hash will be added when Stage 0 lands).
- Pre-flight check: `python scripts/v3e_boccaletti_preflight.py` (must exit 0 before any HPC submit).
- Any change here invalidates all completed Phase E runs and requires a re-baseline.

## Decision log (append-only, dated)

| Date | Decision | Rationale | Owner |
|------|----------|-----------|-------|
| 2026-05-24 | `BPEAK=0.01` | κ·Bpeak·σ_B·√(2π)/2 ≈ 0.063 — perturbative, far from any Boccaletti node at nπ (n≥1) | stage0-bootstrap |
| 2026-05-24 | `R_ATLAS=0.4` | brackets ~35% of T4 v3 MAP magnitude (\|MAP\|≈1.16); near-perturbative, small enough not to swamp tachyonic boundaries | stage0-bootstrap |
| 2026-05-24 | `K_CARRIER=3` | k·σ_B = 15 → high-frequency Gertsenshtein regime; 5.4 cells/λ at Δz=0.39 | stage0-bootstrap |
| 2026-05-24 | `T_END=60`, `T_CHECK_1=40`, `T_CHECK_2=60` | wavepacket front clears B-field at t≈50; t_check_1 is mid-transit (A-plateau diagnostic), t_check_2 is post-transit | stage0-bootstrap |
| 2026-05-24 | `BC=periodic` + dual-Gaussian B-field | required for modal solver (post-#367 fix in v0.42.0/0.42.1); periodicity preserved by `zc1 + zc2 = L` | stage0-bootstrap |
| 2026-05-24 | TT-IC + theory unconstrained (no TT gauge) | per issue #167; validated pattern in `examples/gertsenshtein/theory_ungauged.toml` | stage0-bootstrap |
| 2026-05-24 | R²/R̃² family explicitly stripped from theory roster | Ostrogradsky ghosts via constraint promotion; not unblocked by geometry change. Stripped variants of T7/T8 retained (parity-odd safe couplings) | stage0-bootstrap |

## Stage progress (toggle ☐→☑ in place; add "(jobid 12345)" inline as HPC jobs land)

### Stage 0 — Infrastructure

- [x] 0.1 `scripts/hpc_submit_drafts/v3e_localised/_geometry.env` created
- [x] 0.2 `scripts/v3e_boccaletti_preflight.py` created + 6 unit tests pass
- [x] 0.3 `tidal/measurement/_phase_e_transit.py` created (4 diagnostics) + 5 unit tests pass
- [x] 0.4 `docs/PHASE_E_TRACKER.md` (this file) created
- [ ] 0.5 Phase E memory files + MEMORY.md index updated
- [ ] 0.6 `docs/V3_PHASE_E_DESIGN.md` updated (E.1–E.3 in-progress, stale #367 cvode note removed, frequency-regime table added)
- [ ] 0.7 `manuscript/planning/phase_e_writeup_inputs.md` created
- [ ] 0.8 Modal solver verified on existing `examples/data/gertsenshtein_e0_dual_gaussian.json`
- [ ] 0.9 Stage 0 commit lands (hash → decision log)

### Stage 1 — Wolfram derivations (interleaved with Stage 2+)

- [ ] 1.1 E.cal `gertsenshtein/theory_ungauged_e_dual_gaussian.toml` (~5 min)
- [ ] 1.1b E.EH `euler_heisenberg/theory_e_dual_gaussian.toml` (~10 min)
- [ ] 1.2 E.T1 `dark_photon_plasma/theory_e_dual_gaussian.toml` (~10 min)
- [ ] 1.3 E.T2 `torsion_gertsenshtein/theory_einstein_cartan_e_dual_gaussian.toml` (~10 min)
- [ ] 1.4 E.T4 `torsion_gertsenshtein/theory_nonminimal_e_dual_gaussian.toml` (~30 min)
- [ ] 1.5 E.T5.3 `torsion_gertsenshtein/theory_general_nonminimal_e_dual_gaussian.toml` (~30–60 min; reused by T5/T5.1/T5.2/NP)
- [ ] 1.6 E.T7s `theory_complete_even_minus_R2_e_dual_gaussian.toml` (~30–60 min)
- [ ] 1.7 E.T8s `theory_complete_odd_minus_R2_e_dual_gaussian.toml` (~30–60 min)

### Stage 2 — Calibration (HARD GATE — must PASS before Stage 3)

- [ ] 2.1 Push E.cal + E.T2 JSONs to HPC
- [ ] 2.2 E.cal corner pair on INTR (amp+sup)
- [ ] 2.3 Verify E.cal stability.json: `P_gertsenshtein > 1e-4 · sin²(arg)` (denominator well-resolved)
- [ ] 2.4 E.T2 corner pair on INTR
- [ ] 2.5 Verify E.T2 stability.json: pre-arrival vacuum reads `< 1e-6 · h₀²` (noise floor)
- [ ] 2.6 Tracker + per-jobid READMEs updated

### Stage 3 — Wave 1 corner plots

- [ ] 3.1 E.EH corner pair (σ, ρ)
- [ ] 3.2 E.T1 corner pair
- [ ] 3.3 E.T4 corner pair
- [ ] 3.4 Tracker + writeup-inputs updated with figure paths

### Stage 4 — Wave 1 atlas (≤5D theories only)

- [ ] 4.1 E.EH atlas (2D → 4 faces × 1 tile) — methodology shake-down
- [ ] 4.2 E.T2 atlas (2–3D → 6 faces × 1 tile)
- [ ] 4.3 E.T1 atlas (4D → 8 faces × 1 tile)
- [ ] 4.4 Render atlas pdfs, READMEs, tracker update

### Stage 5 — Wave 2 (Bahamonde family + NP)

- [ ] 5.1 E.T5 (Bahamonde 5D) corner + atlas (10 faces × 1 tile)
- [ ] 5.2 E.T5.1 (Barker 6D) corner-only (reuses E.T5.3 JSON, prior restricted)
- [ ] 5.3 E.T5.2 (Shapiro 8D) corner-only (reuses E.T5.3 JSON)
- [ ] 5.4 E.T5.3 (Complete 9D) corner-only
- [ ] 5.5 E.NP (algebraic 8D) corner-only — same JSON as E.T5.3, prior pins ξ=0
- [ ] 5.6 Render + tracker update

### Stage 6 — Wave 3 (R²-stripped variants)

- [ ] 6.1 E.T7s corner (atlas if ≤5D after stripping)
- [ ] 6.2 E.T8s corner (parity-odd; new physics)
- [ ] 6.3 Render + tracker update

### Stage 7 — Synthesis + manuscript inputs

- [ ] 7.1 Phase B vs Phase E comparison table → `docs/PHASE_E_GEOMETRY.md`
- [ ] 7.2 Headline results table → `phase_e_writeup_inputs.md`
- [ ] 7.3 Canonical figures copied to `manuscript/figures/`
- [ ] 7.4 `docs/V3_PHASE_TRACKER.md` + `CAMPAIGN.md` Stage E section
- [ ] 7.5 Memory final summary + `sync-claude-memory.sh backup`
- [ ] 7.6 `manuscript/planning/report_plan.md` + `report_checklist.md` updated

## Theory roster — one row per theory

Verdict column: `null` / `amplification` / `suppression` / `anomalous` / `—`.
Stability column: `PASS` / `SOFT-PENALIZED` / `CATASTROPHIC`.

| Stage | Theory | Derivation | Corner amp | Corner sup | Atlas | Stability | Verdict |
|-------|--------|------------|------------|------------|-------|-----------|---------|
| E.cal | gertsenshtein_ungauged | ☐ | ☐ jobid | ☐ jobid | n/a (positive control) | — | — |
| E.EH  | Euler-Heisenberg + EM (σ, ρ) | ☐ | ☐ jobid | ☐ jobid | ☐ jobid (2D atlas) | — | — |
| E.T1  | DP-plasma | ☐ | ☐ jobid | ☐ jobid | ☐ jobid | — | — |
| E.T2  | Einstein-Cartan minimal | ☐ | ☐ jobid | ☐ jobid | ☐ jobid | — | — |
| E.T4  | Ricci-EM nonminimal (δ₁) | ☐ | ☐ jobid | ☐ jobid | n/a (1D) | — | — |
| E.T5  | Bahamonde 5D | ☐ (uses E.T5.3) | ☐ jobid | ☐ jobid | ☐ jobid (5D) | — | — |
| E.T5.1 | Barker 6D | reuses E.T5.3 | ☐ jobid | ☐ jobid | n/a | — | — |
| E.T5.2 | Shapiro 8D | reuses E.T5.3 | ☐ jobid | ☐ jobid | n/a | — | — |
| E.T5.3 | Complete 9D | ☐ | ☐ jobid | ☐ jobid | n/a | — | — |
| E.NP  | Non-prop (ξ=0 prior on E.T5.3) | reuses E.T5.3 | ☐ jobid | ☐ jobid | n/a | — | — |
| E.T7s | Complete-even minus R²/R̃² | ☐ | ☐ jobid | ☐ jobid | (≤5D only) | — | — |
| E.T8s | Complete-odd minus R²/R̃² | ☐ | ☐ jobid | ☐ jobid | (≤5D only) | — | — |

## Per-theory artefact index (append per HPC pull)

_None yet — fill as `hpc_results/<jobid>/phase_e/<theory>/README.md` files land._

## Open blockers (dated, drop when resolved)

- [ ] (none yet)

## Resolved blockers (archive)

- _none yet_

## Update protocol (codified — every session follows this)

1. Toggle ☐→☑ in place as tasks complete.
2. Add jobids inline in roster cells when HPC submits land.
3. Append to decision-log only — never rewrite history.
4. Append to per-session footer at session end (see below).
5. Update "Last updated" + "Current stage" + "Next action" lines at the top with every commit.
6. NEVER delete a row — strike-through with `~~` and a reason if a decision changes.

Tracker updates are coupled to the SAME commit that did the work (derivation, code,
pull). Never a separate "tracker update" commit.

## Per-session footer (append at session end; never edit prior entries)

### Session 2026-05-24 (stage0-bootstrap)

- Completed: 0.1, 0.2, 0.3, 0.4
- In flight at session end: none
- Decisions made: BPEAK=0.01, R_ATLAS=0.4, K_CARRIER=3, T_END=60, BC=periodic, TT-IC unconstrained, R²/R̃² stripped from roster (see decision log above)
- Issues raised: _none yet_
