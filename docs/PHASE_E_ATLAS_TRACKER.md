# Phase E Atlas Tracker — Uniform-Field Cubed-Sphere Campaign

**Last updated:** 2026-05-29 (session: stage-a0-bootstrap)
**Current stage:** Stage A0 — infrastructure
**Next action:** push v3_atlas/ to HPC, then submit Phase α S1 (T4 + T1).

## Quick handoff (read first if resuming cold)

1. **Where we are**: Stage A0 infrastructure (config + runner + wrapper) committed locally; not yet pushed to HPC.
2. **What's running**: nothing.
3. **What to do next**: `bash scripts/hpc_shuttle.sh push`, then `bash scripts/hpc_submit_drafts/v3_atlas/submit_atlas_slot.sh t4 t1`.

## Context

- Phase E **localised** is paused (see `docs/PHASE_E_TRACKER.md` PAUSED banner). Localised infrastructure preserved unchanged in `scripts/hpc_submit_drafts/v3e_localised/`.
- Pivot direction (user, 2026-05-29): demonstrate **cubed-sphere atlas methodology** on **v3 uniform-field** posteriors. Hard deadline 2026-06-01.
- Approved plan: `/home/vscode/.claude/plans/currently-we-are-working-goofy-oasis.md`.

## Decision log (append-only)

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-05-29 | Pivot localised → uniform-field atlas | Localised overnight jobs 29714807/29715239 hit walltime; pre-#384 per-eval cost makes landscape PolyChord infeasible in remaining 3 days |
| 2026-05-29 | `R_LO=0.35, R_HI=0.45` (narrow shell) | Sphere-volume bias r^(N-1) bounded by 1.286^8 ≈ 6 even at ndim=9 (cf. 10^54 for the smoke tile 29204991's [1e-3, 1e3]) |
| 2026-05-29 | amp likelihood only (no sup) | Methodology demo; halves cost; angular structure visible from amp alone |
| 2026-05-29 | M=1 single tile per face | First-pass coverage; sub-tile subdivision deferred |
| 2026-05-29 | Order T4 → T1 → T5.0 → T5.1 (Phase α INTR); T5.2 → T5.3 (Phase β overnight) | Minimum-viable T4 first; multi-theory packing into INTR slots maximises throughput |
| 2026-05-29 | INTR validate before standard-queue scale (HARD GATE) | Avoid wasting 12+ h overnight runs on broken pipelines (lesson from 29714807) |

## Stage progress

### Stage A0 — Infrastructure

- [x] A0.1 `scripts/hpc_submit_drafts/v3_atlas/_atlas_config.env`
- [x] A0.2 `scripts/hpc_submit_drafts/v3_atlas/run_atlas_slot.sh` (multi-theory parallel launcher)
- [x] A0.3 `scripts/hpc_submit_drafts/v3_atlas/submit_atlas_slot.sh` (local wrapper: book INTR + ssh-launch)
- [x] A0.4 `docs/PHASE_E_ATLAS_TRACKER.md` skeleton (this file)
- [ ] A0.5 PAUSED banner on `docs/PHASE_E_TRACKER.md`
- [ ] A0.6 push to HPC + verify version aligned
- [x] A0.7 commit Stage A0 infrastructure

### Stage A0b — fixes uncovered during first launch

- [x] CLI: add "atlas" to `--type` choices in `tidal/cli/__init__.py` (dispatch existed in `_plot_command.py:322` but choices list rejected it)
- [x] Runner: source venv + load modules + use absolute path `${TIDAL_ROOT}/.venv/bin/tidal` so srun's PMI handoff can exec
- [x] Runner: `--output` path is the SURVEY directory not face-specific; `tidal sample --joint-prior` writes to `<survey>/<face_label>_tile<sub>/` automatically (psalter convention)

### Stage A1 — Phase α (INTR; v3 order; HARD GATE between slots)

- [ ] **S1**: book INTR + launch `submit_atlas_slot.sh t4 t1` (jobid: ☐)
- [ ] **S1 verify**: per-tile inference.json ≥75% faces complete; atlas.pdf renders for both theories
- [ ] **HARD GATE**: at least one valid atlas.pdf from S1 before proceeding
- [ ] **S2**: `submit_atlas_slot.sh t5_0 t5_1` (jobid: ☐)
- [ ] **S2 verify**: render atlas.pdfs; check `.resume` if any chains incomplete

### Stage A1b — Phase β (overnight standard-queue, ONLY after Phase α succeeds)

- [ ] Decision: proceed to β based on S1+S2 results
- [ ] **B1**: submit_atlas_long.sh t5_2 t5_3 (jobid: ☐) — standard-queue sapphire, ~3-hour walltime, cpu=448
- [ ] B1 verify: render + document
- [ ] Optional: re-run T4 at PUB quality (nlive=200) in B1 for sharper headline fig

### Stage A2 — Synthesis (after last completed theory)

- [ ] A2.1 cross-theory comparison table (`docs/PHASE_E_ATLAS_GEOMETRY.md`)
- [ ] A2.2 canonical figure for manuscript Results section (T4 atlas → `manuscript/figures/atlas_t4.pdf`)
- [ ] A2.3 `manuscript/planning/phase_e_writeup_inputs.md` final pass
- [ ] A2.4 memory backup (`bash .devcontainer/scripts/sync-claude-memory.sh backup`)

## Theory roster — fill jobids + verdicts as each lands

| Theory | ndim | faces | jobid | r_lo | r_hi | n_completed | atlas.pdf | verdict |
|--------|------|-------|-------|------|------|-------------|-----------|---------|
| T4 | 4 | 8 | 29855408 (S1) | 0.35 | 0.45 | 8/8 | ☑ atlas_t4.pdf | informative-amp; δ₁ dir suppressed (logZ=−1.82 face 4) |
| T1 | 4 | 8 | 29855408 (S1) | 0.35 | 0.45 | 8/8 | ☑ atlas_t1.pdf | strong amp; mA² dir face 02m suppressed (logZ=−0.01 vs +9) |
| T5.0 | 5 | 10 | ☐ | 0.35 | 0.45 | ☐ | ☐ | ☐ |
| T5.1 | 6 | 12 | ☐ | 0.35 | 0.45 | ☐ | ☐ | ☐ |
| T5.2 (bonus) | 8 | 16 | ☐ | 0.35 | 0.45 | ☐ | ☐ | ☐ |
| T5.3 (bonus) | 9 | 18 | ☐ | 0.35 | 0.45 | ☐ | ☐ | ☐ |

Verdict column: "informative-amp" / "flat" / "multimodal" / "—".

## Open blockers (dated, drop when resolved)

- (none yet)

## Per-session footer (append at end of each session, never edit prior entries)

### Session 2026-05-29 (stage-a0-bootstrap)

- Completed: A0.1–A0.4 (config, runner, wrapper, tracker skeleton)
- In flight: none
- Decisions: pivot to atlas-on-uniform; narrow shell r=[0.35,0.45]; multi-theory per INTR slot
