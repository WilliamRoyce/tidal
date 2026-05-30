# Phase E Atlas Tracker — Uniform-Field Cubed-Sphere Campaign

**Last updated:** 2026-05-30 (session: stage-a1-s3-fresh-solo)
**Current stage:** Stage A1 S3 — FRESH SOLO per theory (T5.0 then T5.1)
**Next action:** book INTR slot, launch T5.0 fresh-solo (10 chains × 11 ranks = 110 cores), monitor mid-slot

## Quick handoff (read first if resuming cold)

1. **Where we are**: T4 + T1 atlases delivered (S1, jobid 29855408). T5.0 + T5.1 had 3 failed resume attempts due to PolyChord re-verification contention with 22 simultaneous srun steps. Pivoting to FRESH SOLO per theory.
2. **What's running**: nothing (last slot 29878165 expired 11:21 BST).
3. **What to do next**: edit runner for FRESH mode + solo ranks, push, book INTR slot, launch T5.0 alone, **actively monitor mid-slot every 15 min**.

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

- [x] **S1**: `submit_atlas_slot.sh t4 t1` — jobid **29855408**, both atlases rendered, T4 + T1 informative
- [x] **S1 verify**: 8/8 + 8/8 complete; atlas_t4.pdf + atlas_t1.pdf in `hpc_results/29855408/figures/`; per-face corners rendered
- [x] **HARD GATE**: passed (T4 + T1 atlases delivered)
- [x] **S2**: `submit_atlas_slot.sh t5_0 t5_1` — jobid **29855408** (S2 phase, packed 4 ranks/face) — hit walltime, all 22 chains have .resume
- [x] **S2 resume #1**: jobid **29858312** (1h INTR, 4 ranks/face uniform) — T5.0 +260 ndead, T5.1 +220 ndead; none completed
- [x] **S2 resume #2**: jobid **29878165** (1h INTR, smart-pack T5.1=6 ranks, T5.0=4 ranks) — only ~12 of 22 chains progressed; ~10 chains stuck in resume-verification for the whole slot. **Diagnosis: 22 simultaneous `srun --exclusive --read-resume` steps starve each other during PolyChord's live-point re-verification phase.**
- [ ] **S3 strategy B (FRESH SOLO per theory)**: T5.0 alone (10×11=110 cores) → T5.1 alone (12×9=108 cores)

#### Learnings from S2 (CRITICAL — apply to all future multi-chain INTR slots)

1. **Do not launch >16 simultaneous srun-resume steps on a single node**. PolyChord re-verifies all live points after `--read-resume`, and 22 chains × 4-6 ranks competing for cores+memory means most chains never finish re-verification within 1h walltime.
2. **For FRESH chains**, the bottleneck is per-chain compute (likelihood evals), so packing is fine — but stick to ≤12 chains per slot at higher ranks.
3. **Active mid-slot monitoring is mandatory**: every 15 min check that all chains' `_chains/tidal.stats` timestamps are advancing. If any chain's stats hasn't updated in 30+ min while others have, that chain is stuck — react (cancel + rebook with different strategy).
4. **Smart-pack the rank distribution by per-theory ndim and convergence target, not by total core count alone**. Lesson: 22 cores fit on 112-core node mathematically but contention is real.

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
| T5.0 | 5 | 10 | 29855408 (S2) + 29858312 + 29878165 (3 INTR slots, all hit walltime) | 0.35 | 0.45 | 0/10 (~4 chains at/past target ndead 625) | ☐ | pending S3 fresh-solo |
| T5.1 | 6 | 12 | 29855408 (S2) + 29858312 + 29878165 (3 INTR slots, all hit walltime) | 0.35 | 0.45 | 0/12 (~5 chains at/past target ndead 750) | ☐ | pending S3 fresh-solo |
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
