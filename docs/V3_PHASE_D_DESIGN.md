# Phase D — Manuscript implications of the v3 architecture

> **Correction (2026-08-19, additive — the tables below are preserved as written):**
> the `arctan_uniform` prior columns describe the range as degrees mapped through
> `tan`, giving ±57.3. The sampler never did that: `low`/`high` are ignored
> entirely and the support is fixed at ±tan(π/2 − 0.05) ≈ **±19.98** (GH #425).
> Read every ±57.3 in this file as ±19.98. Marginal D_KL values quoted here also
> pre-date the GH #420 estimator fix — see `docs/RESULTS_AMENDMENTS.md`.

**Created:** 2026-05-11
**Status:** PENDING — gated on Phase B convergence
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md), [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md)

## Why this document exists

The v3 architecture pivot (post-2026-05-08 supervisor meeting) materially changes the manuscript: methodology, headline numbers, and the interpretive frame all need updates. This document captures the manuscript update sequence so future sessions can execute it without re-deriving the plan from scratch.

## Update sequence

Order matters — methodology section must land before results can reference it, results section must land before discussion can interpret it.

### D.1 — `manuscript/sections/computational_approach.tex` (methodology)

Add §"Tachyon-permissive inference architecture" describing:

- The three principles (probe-as-measurement-not-gate, soft penalties, compactified priors)
- The soft-floor logL = −100 + Normal(0, σ_explore) mechanism with rationale (sampler-gradient information in failure regions)
- The `run_status` taxonomy (success, simulation_diverged, simulation_failed, kinetic_error, metric_nan, metric_missing, logl_minus_inf, below_noise_floor, exception) — defined in `tidal.measurement._run_stages.RunStatus`. `tachyonic_gated` was dropped with the `--gated` flag in v0.49.6 and survives only in chains recorded before then; `simulation_diverged` was only aspirational until v0.49.6 (GH #480).
- Compactified priors: `arctan_uniform:-89:89` for sign-symmetric dimensionless couplings; `log_uniform:1e-3:1e3` for positive-definite kinetic coefficients
- The cubed-sphere joint-prior chart for higher-dimensional coupling spaces (per parallel session's psalter-aligned implementation in `tidal/inference/_sphere.py`)
- Cite `docs/V3_ARCHITECTURE.md` (commit reference) and `docs/meetings/2026-05-08_supervisor.md` (supervisor pivot)

Text can be drafted *now* (architecture is settled); placeholder numbers `[B.2 log Z]` get filled in when chains land.

### D.2 — `manuscript/sections/results.tex` (numbers + comparison)

Two changes:

1. **Replace v2 numbers** with v3 chain results. For each model (D1, Stage A, D2.0–D2.3), report: log Z (v3), joint D_KL (v3), per-coupling marginal D_KL table, MAP location, posterior 95% credible interval.
2. **Add §"v2 vs v3 comparison" table** showing per-model what changed. Use the auto-generated tables from `scripts/v3_v2_comparison.py` (e.g., `docs/comparison/d1_amp_v2_v3.md`).

The headline framing shifts from "log Z = X (compelling evidence)" to **"per-coupling marginal D_KL = Y nats (this coupling carries Y nats of signal-shape information)"**. The supervisor's framing: results become "this operator gives structure / this operator doesn't do much" rather than "A_max = 38".

### D.3 — `docs/PHASE_6_COMPARISON.md` (campaign comparison)

Append §"Phase 7 / v3 architecture" with:

- Architecture pivot summary (one paragraph, pointer to V3_ARCHITECTURE.md)
- Per-model comparison table (rows: D1 amp/sup, Stage A amp/sup, D2.0-D2.3 amp/sup; columns: v2 log Z, v3 log Z, v2 joint D_KL, v3 joint D_KL)
- Per-model per-coupling shift table
- One-paragraph interpretation per model — which couplings became measurable that weren't, which posterior modes appeared that v2 couldn't reach

Generated using `scripts/v3_v2_comparison.py` per chain pair, then aggregated.

### D.4 — `manuscript/sections/discussion.tex` (interpretation)

Updates needed:

- Acknowledge the methodology shift in framing (v2 chains conditioned on stability; v3 integrates over wider support)
- Address the **log Z incomparability across architectures** (footnote: v3 includes prior support v2 excluded, so log Z differences mix posterior shape with prior compactness)
- The supervisor's framing question: when the v3 architecture admits previously-rejected regions, what does the resulting posterior *mean* physically? Discussion section is where the linearized-theory-validity caveats appear
- If Phase E results are in: discuss geometry-dependence and which numbers headline
- If Phase A-γ results are in: discuss γ_conversion as a fast proxy and its validation

### D.5 — `manuscript/sections/abstract.tex` + `manuscript/sections/introduction.tex` (framing)

Lightest touch — update headline numbers in abstract; introduction paragraph on methodology gets a "tachyon-permissive" qualifier. Cross-references to the methodology section.

## v2-vs-v3 comparison table template

Per model, generated by `scripts/v3_v2_comparison.py`:

```text
# <Model> — v2 vs v3 comparison

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z          | +0.720 ± 0.024 | +13.288 ± 0.125 | +12.57 nats |
| ESS            | 5997           | 1410            | -4588       |
| Joint D_KL     | 0.76 nats      | 4.77 nats       | +4.01 nats  |
| n_samples      | 12524          | 3651            | —           |

| Param  | v2 prior        | v3 prior                   | v2 MAP  | v3 MAP  | shift (σ) | v2 D_KL | v3 D_KL |
| ---    | ---             | ---                        | ---     | ---     | ---       | ---     | ---     |
| alpha1 | uniform[-1..1]  | arctan[-89°..89°] → ±57.3  | -0.422  | -0.539  | -0.18     | 0.02    | 3.20    |
| ...
```

Live exemplar: [docs/comparison/d1_amp_v2_v3.md](comparison/d1_amp_v2_v3.md).

## Geometry decision rule (gates D.2 numbers)

**If Phase E lands before manuscript submission:**

- Compute per-coupling D_KL shift `|D_KL_E - D_KL_B|` per coupling
- If max shift < 0.5 nats: Phase B (plane-wave) numbers headline; Phase E appears as appendix sanity-check
- If max shift ≥ 0.5 nats: Phase E (localized wavepacket) numbers headline; Phase B appears as plane-wave-conditioned reference

**If Phase E does not land before submission:**

- Phase B numbers headline with explicit "plane-wave geometry" qualifier in §results
- Discussion section explicitly flags the unbounded-growth artifact for tachyonic samples (γ_eff·t_end exponential) as a known limitation to be resolved by Phase E in follow-up work

## Chain-pull workflow (post-HPC-job)

Every HPC chain that lands must be processed through the same workflow to ensure consistency:

1. **Pull artifacts**: `bash scripts/hpc_shuttle.sh pull <jobid> --src /rds/user/wr286/hpc-work/tidal/hpc_results/<jobid>/<chain_name>` → drops files in `hpc_results/<jobid>/`
2. **Generate corner plot** (mandatory): `uv run tidal plot hpc_results/<jobid> --type corner --show-rejected-prior --output hpc_results/<jobid>/corner_v3.png`
3. **Verify chain health**: inspect `inference.json` (log Z, ESS, MAP, per-coupling D_KL) and `results.csv` for `run_status` fraction breakdown — flag any > 5% non-success fraction (excluding `below_noise_floor`)
4. **Generate v2-vs-v3 comparison** (when v2 reference exists): `uv run python scripts/v3_v2_comparison.py --v2 hpc_results/<v2_jobid> --v3 hpc_results/<v3_jobid> --label "<model>" --output docs/comparison/<model>_v2_v3.md`
5. **Update `docs/V3_PHASE_TRACKER.md`** with status + log Z + ESS for that step
6. **Post results to GitHub issue #355** (Phase B umbrella)
7. **Commit** the comparison doc + tracker update

The plotting and comparison are *mandatory* — a chain without corner plot is half-pulled. Future sessions should never have to ask "did the corner plot get made?"; the answer is always yes, in `hpc_results/<jobid>/corner_*.png`.

## Pointer to Phase A-γ deferred session

The γ_conversion coupling-aware probe (Phase A-γ, GH issues #350, #351) is deferred to its own session. Phase D does **not** block on A-γ:

- If A-γ ships *before* manuscript submission: D.2 can mention γ_conversion as the fast metric used for some chains (with P_max as the verified ranking)
- If A-γ ships *after* submission: Phase D ships under P_max only; γ_conversion appears in follow-up work

The architecture (Phase A) and campaign re-runs (Phase B) are sufficient for a publishable manuscript without γ_conversion.

## Done condition

- `manuscript/sections/computational_approach.tex` has §"Tachyon-permissive inference architecture"
- `manuscript/sections/results.tex` references v3 numbers and contains §"v2 vs v3 comparison"
- `docs/PHASE_6_COMPARISON.md` has §"Phase 7 / v3 architecture" with per-model tables
- `manuscript/sections/discussion.tex` addresses log Z incomparability and methodology shift
- All `docs/comparison/<model>_v2_v3.md` files committed for the campaign models (D1 amp/sup, Stage A amp/sup, D2.0–D2.3 amp/sup)
- Supervisor reviewing draft

## References

- [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md) — canonical architecture reference
- [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md) — Phase E design (gates geometry decision)
- [V3_PHASE_TRACKER.md](V3_PHASE_TRACKER.md) — running checklist
- [scripts/v3_v2_comparison.py](../scripts/v3_v2_comparison.py) — comparison-table generator
- [docs/comparison/d1_amp_v2_v3.md](comparison/d1_amp_v2_v3.md) — live exemplar (B.1 vs v2 28896653)
