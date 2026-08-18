# Phase E Atlas Tracker — Uniform-Field Cubed-Sphere Campaign

**Last updated:** 2026-05-30 (session: stage-a1-s4-standard-queue)
**Current stage:** Stage A1 — T4/T1/T5.0 DONE; T5.1/T5.2/T5.3 running in standard-queue
**Next action:** monitor 29890992-29890994 (sapphire standard); pull + render as each completes

## Quick handoff (read first if resuming cold)

1. **Where we are**: T4, T1, **T5.0** atlases all DELIVERED. T5.1/T5.2/T5.3 submitted to standard-queue sapphire (single-node each, validated runner). 4 INTR slots taught us: 22-simultaneous resume contention is real, ~10·nlive convergence (not 5·nlive), one INTR slot can't finish ndim ≥ 5 atlases — hence standard-queue.
2. **What's running**: 29890992 (t5_1, 6h), 29890993 (t5_2, 9h), 29890994 (t5_3, 12h) — all RUNNING since 16:55 BST 2026-05-30.
3. **What to do next**: poll periodically (file-existence on slurm log + inference.json count); pull each as it completes; render atlas + per-face corners; commit per-theory.

## Context

- Phase E **localized** is paused (see `docs/PHASE_E_TRACKER.md` PAUSED banner). Localized infrastructure preserved unchanged in `scripts/hpc_submit_drafts/v3e_localised/`.
- Pivot direction (user, 2026-05-29): demonstrate **cubed-sphere atlas methodology** on **v3 uniform-field** posteriors. Hard deadline 2026-06-01.
- Approved plan: `/home/vscode/.claude/plans/currently-we-are-working-goofy-oasis.md`.

## Decision log (append-only)

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-05-29 | Pivot localized → uniform-field atlas | Localized overnight jobs 29714807/29715239 hit walltime; pre-#384 per-eval cost makes landscape PolyChord infeasible in remaining 3 days |
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
- [x] **S3a T5.0 FRESH SOLO**: jobid **29879757** (1h INTR, 10 × 11 ranks = 110 cores, cpu-r-12). All 10 chains advanced cleanly (20-30 ndead/min/chain, no stalls). Final ndead 686-1466. **0/10 inference.json** — slot hit walltime before precision_criterion=0.01 triggered. Rich .resume checkpoints from uncontended fresh run.
- [x] **S3b T5.0 RESUME + T5.1 FRESH sequential**: jobid **29882283** (1h INTR, cpu-r-111, TIDAL_SEQUENTIAL=1, runner added auto-detect resume). **5/10 T5.0 chains COMPLETED** (03p, 04m, 04p, 05m, 05p have inference.json). 5 incomplete (01m=1108, 01p=1146, 02m=1068, 02p=968, 03m=1104). T5.1 never started — T5.0 dominated the slot.
- [x] **S3c T5.0 RESUME-FINISH (5 incomplete only)**: slot **29884704** with `TIDAL_FACES_T5_0="1 2 3 4 6"` override — **5/5 chains converged** (ndead 1956-2116). t5_0_DONE sentinel triggered; T5.1 fresh started at 15:43 BST with ~5 min progress (ndead 150-310) before slot expiry.

### Stage A1b — Phase β (standard-queue sapphire, FRESH per theory, single-node)

- [x] **β.1 T5.1**: jobid **29890992** — `submit_atlas_standard.sh t5_1` (108 cores, 6h walltime, sapphire). RUNNING since 2026-05-30 16:55 BST.
- [x] **β.2 T5.2**: jobid **29890993** — same wrapper, t5_2 (112 cores, 9h walltime). RUNNING since 2026-05-30 16:55 BST.
- [x] **β.3 T5.3**: jobid **29890994** — same wrapper, t5_3 (108 cores, 12h walltime). RUNNING since 2026-05-30 16:55 BST.
- [ ] β.4: poll periodically (file-existence + sacct one-shot); pull + render each as it completes.
- [ ] β.5: render atlas + per-face corners per theory (mirror T5.0 figure pattern).

#### S4 standard-queue strategy notes

- **One sbatch per theory**: independent scheduling, isolated failure modes. All 3 submitted in <30s and all started RUNNING immediately (no queue wait under HANDLEY-SL3-CPU).
- **Single sapphire node per theory**: avoids multi-node MPI complexity (which had caused PMI issues in past hpc work).
- **TIDAL_FRESH=1**: defensive clean start; no resume-contention risk (the issue that killed slot 29878165).
- **Parallel mode within sbatch** (not sequential): each theory's faces all run concurrently; no cross-theory contention since they're in different sbatches.
- **Walltime cushion**: 33-50% beyond T5.0-extrapolated estimate. T5.0 took ~3h total at 11 ranks/face; T5.3 at 6 ranks/face with 1.5× more work per dead point is the conservative 12h ceiling.

#### S3 additional learnings (slot 29882283)

- **PolyChord precision_criterion=0.01 convergence ndead ≈ 1700+ for T5.0 (ndim=5, nlive=125)** — that is ~14·nlive, NOT 5·nlive. Heuristic was way too optimistic. Use 10·nlive minimum for budget planning; expect up to 15·nlive.
- **Sequential mode validated**: T5.0 ran cleanly through its waitloop; T5.1 was correctly queued and never launched (no contention).
- **Cross-chain ndead variance is large**: 04p completed first (resume started at 1466, converged ~1675). 02p was slowest (969 final). Same nlive/num_repeats, same likelihood, but different starting random state → different number of iterations needed. Budget for the slowest chain.
- **Auto-detect resume works**: chains with .resume staged in atlas_t5_0/ correctly resumed; T5.1 (no staged dirs) would have started fresh if the slot had reached it.

#### S3 LEARNINGS so far (apply to S3b/c)

- **FRESH SOLO works**: 10 chains × 11 ranks/face with no other-theory contention → all chains advanced cleanly at 20-30 ndead/min/chain. No stuck chains. This validates the strategy.
- **PolyChord precision_criterion=0.01 needs more dead points than 5·nlive heuristic suggests**: with nlive=125, the heuristic said convergence at ndead~625, but at ndead=1466 some chains still hadn't converged. Plan for ~10·nlive worst-case wall time.
- **Mid-slot monitoring caught the success**: iter 1→6 trajectory showed all chains advancing in lockstep (~+140 ndead per 5 min for working chains). Compare to S2 resume slots where some chains showed flat ndead — the diff is obvious in 5-min cadence and lets us react.

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
| T5.0 | 5 | 10 | 29879757 (fresh-solo) + 29882283 (resume-seq, 5/10) + 29884704 (resume-finish, 5/10) | 0.35 | 0.45 | **10/10** | ☑ atlas_t5_0.pdf | DELIVERED — see hpc_results/29884704/figures/ |
| T5.1 | 6 | 12 | 29890992 (standard-queue sapphire, 6h, 108 cores) | 0.35 | 0.45 | RUNNING | ☐ | RUNNING since 2026-05-30 16:55 BST |
| T5.2 (bonus) | 8 | 16 | 29890993 (standard-queue sapphire, 9h, 112 cores) | 0.35 | 0.45 | RUNNING | ☐ | RUNNING since 2026-05-30 16:55 BST |
| T5.3 (bonus) | 9 | 18 | 29890994 (standard-queue sapphire, 12h, 108 cores) | 0.35 | 0.45 | RUNNING | ☐ | RUNNING since 2026-05-30 16:55 BST |
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
