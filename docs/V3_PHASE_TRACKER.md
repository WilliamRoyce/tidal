# v3 Phase Tracker

**Created:** 2026-05-10
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md)
**Purpose:** running checklist of v3 phases — updated as each task converges (mirrors how Phase 6.J/K/L/M were tracked in CAMPAIGN.md).

## Phase A — Soft-penalty refactor + compactified-prior + de-pruning (in progress)

| Step | Description | Status | Issue / Notes |
| --- | --- | --- | --- |
| A.0a | `docs/V3_ARCHITECTURE.md` written | ✅ done | 2026-05-10 |
| A.0b | `docs/V3_PHASE_TRACKER.md` written | ✅ done | (this file) |
| A.0c | `docs/V3_PHASE_E_DESIGN.md` written | ✅ done | Deferred-but-documented |
| A.0d | `CAMPAIGN.md` updated with v3 architecture pointer | ✅ done | |
| A.0e | GitHub issues created for v3-A/A-γ/B/D/E/C tracking | ✅ done | #345–#355 |
| A.0f | Obsolete HPC jobs (28982006, 28985879) cancelled | ✅ done | 2026-05-10 |
| A.0g | Phase A.0 persistence committed | ✅ done | 62c7ac9 |
| A.1 | Soft-penalty refactor in `tidal/inference/_likelihood.py` | ✅ done | [#345](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/345) — Hwang-Noh + P_max>2 cap removed; soft floor `-100 + Normal(0, σ)` for sim/NaN/exception; distinct `run_status` tags. New `--gated` and `--soft-floor-noise SIGMA` flags. |
| A.2 | Lagrangian de-pruning audit + per-param `arctan_uniform` scripts | ✅ done | [#346](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/346) — `docs/lagrangian_depruning_audit.md` + 12 v3_permissive scripts |
| A.3 | Corner-plot upper-triangle removal in `tidal/inference/_visualize.py` | ✅ done | [#347](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/347) — `_hide_upper_triangle()` helper |
| A.4 | Soft-penalty tests (`test_likelihood_*.py`) | ✅ done | [#348](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/348) — 3 new modules, 36 new test cases, all green |
| A.5 | D1 v1 chain replay sanity check | ✅ done | [#349](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/349) — `examples/data/v3_d1_replay/`; v3 admits 14.3% more samples than v2 |
| A.6 | Commit + version bump | ✅ done | 400a455 (v0.39.0) |

## Phase A-γ — γ_conversion (deferred to its own session)

| Step | Description | Status |
| --- | --- | --- |
| A-γ.1 | Refactored γ_conversion with multi-t_test sampling + log-zero clamp in `_stability.py` | ⏳ deferred |
| A-γ.2 | `gamma_conversion` metric in `parse_likelihood()` (opt-in) | ⏳ deferred |
| A-γ.3 | 10-point D1 P_max correlation validation | ⏳ deferred |
| A-γ.4 | Decision: ship γ_conversion for amp campaigns or keep as metadata | ⏳ deferred |

GH issues `[v3-A-γ] γ_conversion: refactored definition + multi-t_test + log-zero clamp` and `[v3-A-γ] γ_conversion: P_max correlation validation` track the work.

## PolyChord sampling budget reference

Previous campaigns used `nlive=1500–2400` to compensate for divergence-guard rejections
(many prior samples were rejected as tachyonic, requiring a large live set to explore the
posterior). v3 permissive sampling removes that need — PolyChord's own recommendations apply
directly.

| Theory | ndim | Landscape nlive | Publication nlive | num_repeats |
| --- | --- | --- | --- | --- |
| D2.0 Bahamonde | 5 | **125** (25×ndim) | 250 (50×ndim) | 10 (2×ndim) |
| D2.1 Barker | 6 | **150** | 300 | 12 |
| D2.2 Shapiro | 8 | **200** | 400 | 16 |
| D2.3 Full | 9 | **225** | 450 | 18 |

**Landscape pass** = current B.5 goal — sufficient for Bayes-factor comparison and MAP
identification. **Publication pass** = future rerun for tighter logZ (±0.1 vs ±0.3) and
credible-interval reporting in the manuscript. All `v3_permissive/*.sh` scripts use landscape
settings; future scripts go in `v3_publication/`.

## Phase B — Campaign re-runs

| Step | Description | Status | Notes |
| --- | --- | --- | --- |
| B.0a | Joint-prior smoke comparison (D1 amp, 1 tile, INTR) | ⚠️ done — **does not pass adoption criterion** | 29204991 (8 min INTR). log Z=+13.76±0.09 (matches B.1 +13.29), ESS=1398 (matches B.1 1410), but **MAP α₁=+582, α₂=−575, δ₁=+575** vs B.1's ±0.5 — joint prior explores a different region than per-param. Per-coupling D_KL much lower (0.27–0.34 vs 2.96–3.20). Joint prior's r∈[1e-3,1e3] × sphere volume scaling concentrates at r_hi=1000; needs r_hi tuning. **Decision: keep per-param arctan for B.5.** See `docs/comparison/d1_amp_joint_v2_v3.md` |
| B.1 | D1 amp v3 smoke (mid-res INTR) | ✅ done | 29149987 (23 min INTR, grid=128/nlive=300). log Z=+13.29±0.13 (vs v2 +0.72; +12.5 nats), ESS=1410/3651, MAP δ₁=4.80 (outside v2 prior). Per-coupling D_KL: α₁=3.20, α₂=2.96, α₃=0.02, δ₁=1.44 nats. 100% success run_status. Post-hoc probe sweep flags 45% of prior as v2-tachyonic. See `docs/comparison/d1_amp_v2_v3.md` |
| B.2 | D1 amp v3 publication run | 🟡 pending | 29189748 (standard 8h) — submitted 2026-05-10 |
| B.3 | D1 sup v3 paired with B.2 | 🟡 pending | 29189761 (standard 12h) — submitted 2026-05-10 |
| B.4a | Stage A amp INTR smoke | ✅ done | 29189966 (8 min INTR). log Z=+9.31±0.13 (vs v2 −0.07; +9.4 nats), ESS=2649/6464. **v2 "null verdict" was an artefact of narrow priors** — v3 reveals joint D_KL=6.59 nats. MAP shifts ~1.1–2.0σ on all 4 params. See `docs/comparison/stage_a_amp_v2_v3.md` |
| B.4b | Stage A sup INTR smoke | ✅ done | 29199129 (22 min INTR). logZ=+4.08±0.12 (vs v2 +0.66; +3.4 nats). ESS=4961. Joint D_KL=7.30 nats. MAP mA2=605 (v2: 0.97 — massive prior shift). Per-param D_KL: deltam=3.02, xi=1.05, mA2=0.86, alpha3=0.63 nats. |
| B.4-full-amp | Stage A amp v3 publication | 🟡 pending | 29205968 (standard 6h) — submitted 2026-05-11 |
| B.4-full-sup | Stage A sup v3 publication | 🟡 pending | 29205982 (standard 12h) — submitted 2026-05-11 |
| B.5 | D2.0–D2.3 v3 paired re-runs (8 chains) | ✅ done (landscape pass, 2026-05-17) | **All 8 chains captured at landscape nlive=25×ndim** — replaces broken nlive=1500–2400 runs that suffered from divergence-guard-era over-resolution. **Single-theory INTR**: D2.0 amp 29229768 (logZ=+9.29±0.06), D2.1 amp 29256858 (logZ=+9.58±0.05), D2.0 sup **29468539** (9 min, logZ=+5.382±0.200, ESS=1302). **Interactive INTR batch 29468763** (5/5 sequential @ 76 ranks, ~58 min): D2.1 sup (8m, logZ=+2.415±0.170), D2.2 amp (10m, logZ=+8.889±0.135), D2.2 sup (22m, logZ=+6.212±0.176), D2.3 amp (10m, logZ=+10.034±0.142), D2.3 sup SIGTERM at walltime. **Interactive INTR batch 29471255** (Phase 1: D2.3 sup at 76 ranks 26m, logZ=+2.662±0.156, ESS=2870; Phase 2: parallel mpirun test — D2.0 sup pub + D2.1 sup pub at 38 ranks each completed in 19m, **1.7× speedup over sequential**). |
| B.6 | Bounds-dependence cross-check at L=75 under v3 | ✅ done | 29205638 (25 min INTR). log Z=+13.18±0.12 (nearly identical to B.1's +13.29). Per-coupling D_KL essentially unchanged. **v3 washes out v2's bounds-dependence non-monotonicity** (v2 L=75 = +0.84, v2 L=100 = +0.68 — non-mono; v3 L=75 ≈ L=100). |

## Phase C — Cubed-sphere coupling-space chart (HANDLED BY PARALLEL SESSION)

See `docs/V3_PHASE_C_REFERENCE.md` once written by the parallel session. This phase is shipped externally to the v3 plan owned here.

Open follow-up question for supervisor: coupling grouping (monolithic / per-Lagrangian-symmetry-class / per-SPO-sector / per-parity / other). GH issue `[v3-C] Coupling-grouping question for supervisor follow-up` tracks.

## Phase D — Manuscript implications

| Step | Description | Status |
| --- | --- | --- |
| D.1 | `manuscript/sections/computational_approach.tex` — add §"Tachyon-permissive inference architecture" | ⏳ pending (after Phase B) |
| D.2 | Update results sections with v3 numbers; comparison table v2 vs v3 | ⏳ pending (after Phase B) |
| D.3 | `docs/PHASE_6_COMPARISON.md` — append "v3 architecture" section | ⏳ pending (after Phase B) |

## Phase E — Localised geometry pivot (DEFERRED, gated on Phase B convergence)

See [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md). Replaces plane-wave + uniform B₀ with Gaussian wavepacket + dual-Gaussian B-field profile; finite interaction time bounds P_max physically.

| Step | Description | Status | Notes |
| --- | --- | --- | --- |
| E.0 | Dual-Gaussian localised B-field empirical validation | ✅ done (2026-05-16) | `examples/gertsenshtein/theory_e0_dual_gaussian.toml` derived. sigB=5, zc1=25, zc2=75, L=100. A(20)/A(10)=1.00 < 1.05 ✓. P_max=0.00348 at t=18 (CVODE). Modal auto-selects but diverges → GH #367; use `--scheme cvode`. |
| E.1 | Wavepacket IC switch for all v3 campaign scripts | ⏳ deferred | After Phase B convergence |
| E.2 | Localised B-field theory TOMLs for all D1/D2/StageA models | ⏳ deferred | After Phase B convergence |
| E.3 | σ_w × σ_B tuning scan | ⏳ deferred | After E.2 |
| E.4 | Re-baseline campaign scripts in `v3e_localised/` | ⏳ deferred | After E.3 |
| E.5 | Re-run highest-value chains under E geometry | ⏳ deferred | D1 amp first |
| E.6 | Phase B vs Phase E comparison + publication decision | ⏳ deferred | After E.5 |

## Sign-off log

(Append entries as phases converge — date, step ID, summary, commit hash.)

- 2026-05-10 — A.0a — `docs/V3_ARCHITECTURE.md` written; canonical architecture reference live.
- 2026-05-10 — A.0b — `docs/V3_PHASE_TRACKER.md` written (this file).
- 2026-05-10 — A.0c — `docs/V3_PHASE_E_DESIGN.md` written; Phase E deferred-but-documented.
- 2026-05-10 — A.0d — `CAMPAIGN.md` updated with v3 architecture pointer and HPC-job cancellations.
- 2026-05-10 — A.0e — GitHub issues #345–#355 created; `v3-architecture` label added.
- 2026-05-10 — A.0f — HPC jobs 28982006 + 28985879 cancelled.
- 2026-05-10 — A.0g — Phase A.0 persistence committed in 62c7ac9.
- 2026-05-10 — A.1 — Soft-penalty refactor in `_likelihood.py`: Hwang-Noh and `P_max>2` cap removed; soft floor `−100 + Normal(0, σ)` for sim/NaN/exception; distinct run_status tags. New `--gated` and `--soft-floor-noise SIGMA` flags wired through `_sample.py` + `cli/__init__.py`.
- 2026-05-10 — A.2 — Lagrangian de-pruning audit (`docs/lagrangian_depruning_audit.md`) committed: no de-pruning needed for Phase A's 6 campaigns; T6/EH applications deferred. Twelve `scripts/hpc_submit_drafts/v3_permissive/` campaign scripts written for D1, Stage A, D2.0–D2.3 paired chains.
- 2026-05-10 — A.3 — Corner-plot upper-triangle scatter rendering disabled via new `_hide_upper_triangle()` helper in `tidal/inference/_visualize.py`.
- 2026-05-10 — A.4 — Three new test modules `tests/test_likelihood_{soft_floor_noise,no_hwang_noh,permissive}.py` (36 cases); 107/107 inference tests pass; ruff/format clean.
- 2026-05-10 — A.5 — D1 v1/v2 chain replay landed at `examples/data/v3_d1_replay/`. Headline: v3 architecture admits 14.3% more samples than v2 (1428/10000 tachyonic samples now contribute via the v2→v3 admission shift on `28982029`). 28520217 v1 chain pre-dates the rejected-prior sidecar so the comparison uses 28982029 instead.
- 2026-05-12 — B.4b — Stage A sup v3 INTR smoke (29199129, 22 min): logZ=+4.08±0.12 (vs v2 +0.66; +3.4 nats), ESS=4961, joint D_KL=7.30 nats. MAP mA2=605 vs v2 MAP 0.97 — massive prior-shift. deltam D_KL=3.02 nats strongest. v3 unambiguously reveals structure masked by v2's narrow priors.
- 2026-05-12 — B.5 (D2.1 amp) — D2.1 Barker amp v3 INTR (29256858, 27 min): logZ=+9.58±0.05 (vs v2 +0.62; +8.96 nats), ESS=9517, joint D_KL=4.13 nats. β₂ MAP=−11.7 (outside v2's narrow [-3..-0.3] prior, 2.5σ shift). All 5 sign-symmetric params D_KL ≈ 2.5-2.8 nats. 27 clusters (multi-modal). See `docs/comparison/d21_barker_amp_v2_v3.md`.
- 2026-05-16 — B.5 (policy) — D2.0 sup truncated at session 7 (supervisor-approved): 58.5K dead pts, 174/401 active clusters, logZ=+134.47±0.453. `*_resume.sh` scripts deleted. 1-INTR-session-per-chain policy baked into `v3_permissive/README.md`. See `docs/comparison/d20_bahamonde_sup_v2_v3.md`.
- 2026-05-16 — E.0 — Dual-Gaussian localised B-field validated: `theory_e0_dual_gaussian.toml` derived (6 fields, 22 H terms). CVODE smoke tests pass; A(20)/A(10)=1.00 < 1.05 ✓; P_max=0.00348 at t=18. Modal solver auto-selects but diverges (GH #367 filed); workaround: `--scheme cvode`.
- 2026-05-17 — B.5 done — All D2.0–D2.3 amp+sup landscape pass captured at nlive=25×ndim (commit d786851). Validation run D2.0 sup 29468539 (9 min) confirmed default settings work. Interactive INTR batch 29468763 ran 4/5 sequential in ~58 min; batch 29471255 closed out D2.3 sup + pioneered 2-way parallel mpirun pattern (1.7× speedup over sequential, both rc=0). Bonus publication-quality reruns of D2.0 sup (nlive=250: logZ=+5.837±0.143, ESS=2404) and D2.1 sup (nlive=300: logZ=+2.524±0.158, ESS=2895) included in 29471255. Parallel mpirun pattern documented in `scripts/hpc_submit_drafts/v3_permissive/run_interactive_batch_2.sh` for future v3_publication and Phase E batches.
