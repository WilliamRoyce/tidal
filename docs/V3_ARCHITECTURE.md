# v3 Inference Architecture (post-2026-05-08 supervisor pivot)

**Created:** 2026-05-10
**Branch:** `hpc/pgt-survey`
**Status:** Phase A in progress; supersedes the canonical-probe ("v2") architecture
**Reference plan:** `/home/vscode/.claude/plans/binary-snacking-chipmunk.md`
**Supervisor meeting:** [docs/meetings/2026-05-08_supervisor.md](meetings/2026-05-08_supervisor.md)

This document is the canonical architecture reference for the v3 era. It survives context loss / session boundaries so any future agent or person picking up the campaign has a complete picture without needing the plan file.

## Why v3?

The 8-May-2026 supervisor meeting (notes [§3](meetings/2026-05-08_supervisor.md#3-stability-filtering-in-the-linearised-regime)) explicitly raised whether the v2 hard tachyonic rejection (τ=0.15 probe gate) is over-conservative. In conversation the supervisor went further and asked for three coupled architectural changes that supersede the canonical-probe direction:

1. **Compactify the parameter space** — replace narrow uniform priors (e.g. α₁ ∈ [-1,1]) with priors covering the full real line, density concentrated near zero. Whiteboard sketch: extract a scale factor from a group of couplings, parameterise the remainder as angles on a sphere, project a cube wireframe onto the sphere for a quasi-rectilinear angular grid. The supervisor sent his unpublished `psalter` Python package (in `psalter.tar.gz`, gitignored) as the canonical reference for this scheme; a parallel session is implementing it (see "Parallel session" below).
2. **Soft penalty instead of `-inf`** — graded penalties for sim divergence/NaN/exception so PolyChord sees a stability gradient. Optional Gaussian noise on the floor (`logL = floor + Normal(0, σ_explore)`) prevents the sampler from seeing a flat plateau in failure regions.
3. **Tachyon-permissive sampling** — don't gate on growing modes. Map the structure of the amplification landscape including unstable regions; downstream interpretation handles whether each growth mode is real physics or a linearisation artefact.

Two further refinements emerged during planning:

4. **Hwang-Noh validity is not a sim-correctness criterion.** P_max < 1 is a statement about whether the linearised theory is *physically meaningful*, not numerical correctness. The Hwang-Noh inline gate (`P_max > 0.5 → -inf`) is removed entirely. P_max is recorded faithfully across all samples; downstream analysis interprets.
5. **Don't prune the Lagrangian with analytical inertness arguments.** Operators currently omitted as boundary terms (e.g. parity-odd FF̃) should be included with free coefficients; the chain's null verdict on a known-inert operator is a *positive* methodology demonstration. Increases dimensionality of D2 parity-odd (was 22 free params) and possibly EH; compactified priors handle the inflation gracefully.

The campaign output shifts from headline log Z and A_max numbers to per-coupling marginal D_KL on the chain posterior — "this operator gives structure / this operator doesn't do much" rather than "A_max = 38".

## Soft-penalty table

| Failure mode | v2 | v3 | run_status tag |
|---|---|---|---|
| Tachyonic (γ_eff > τ=0.15) | `-inf` | **REMOVED**. Sim runs; γ_eff recorded as metadata only. | `success` |
| Non-perturbative (P_max > 0.5) | `-inf` | **REMOVED**. P_max recorded; no upper cap at any value. | `success` |
| Sim divergence (`SimulationDivergedError`) | `-inf` | `logL = -100 + Normal(0, σ_explore)` | `simulation_diverged` |
| Sim exit_code ≠ 0 | `-inf` | same floor + noise | `simulation_failed` |
| NaN/Inf metric | `-inf` | same floor + noise; **distinct tag** for post-chain inspection | `metric_nan` |
| Metric missing from sim output | `-inf` | `-inf` (genuine bug, not parameter-space signal) | `metric_missing` |
| Sim returns finite logL **below** `SOFT_FLOOR_LOGL` | n/a (no analogue in v2) | logL kept verbatim (no clamp); distinct tag so post-chain analysis can filter sub-noise-floor samples from "physical" min/max summaries. Issue #356. | `below_noise_floor` |

Default `σ_explore = 1.0` nat; tunable via `--soft-floor-noise SIGMA`. The `--gated` flag preserves v2 hard-rejection behaviour for reproducibility.

`below_noise_floor` is observational metadata only — sample weights at logL ≈ −101 are ~1e−50 (vs ~1e−2 at MAP), so posterior inference is unaffected. The tag protects diagnostic summaries (corner-plot A-range, headline tables) from quoting numerical-noise floors as physical bounds. The simulation's effective noise floor scales with IC amplitude and solver precision; for typical IC=1e−2 and double-precision modal solver, P_max ≈ 1e−34 is the natural threshold and SOFT_FLOOR_LOGL = −100 (P_max ≈ 1e−44) provides a comfortable margin.

## Compactified prior table (per-param `arctan_uniform` safety-net path)

| Coupling type | v2 | v3 |
|---|---|---|
| Sign-symmetric dimensionless (α, β, δ, χ, ζ) | `uniform:-1:1` or `uniform:-2:2` | `arctan_uniform:-89:89` (degrees → tan(±89°) ≈ ±57; full real line, Cauchy density at 0) |
| Positive-definite kinetic (ξ, α₃, mA²) | `log_uniform:0.05:2` | `log_uniform:1e-3:1e3` (already compactified — 6 decades is "effectively ±∞") |

The cubed-sphere joint-prior path (parallel session) is the alternative; both ship.

## Likelihood metric

- **Amp campaigns**: `--likelihood "P_max:maximize"`. Sim always runs. The simulation faithfully resolves whether tachyonic eigenvalues couple source → target (the eigenvalue alone is coupling-blind: a torsion-sector tachyon decoupled from h_5 / a_1 produces no physical conversion despite huge γ).
- **Sup campaigns**: `--likelihood "P_max:minimize"`. Sim always runs. Suppression integrates IC overlap, phase mismatch, destructive interference — all of which require the actual sim.

A coupling-aware probe metric (`gamma_conversion = log(|target_amp|)/t_test` from the Padé probe output) is **deferred to its own future session** (Phase A-γ) — failure-mode investigation identified showstoppers (log-zero, source-decay-masking, P_max-correlation-unverified) that need refactoring + empirical validation before it can be load-bearing.

## Parallel session

`/home/vscode/.claude/plans/read-this-plan-file-scalable-church.md` is implementing the cubed-sphere joint prior (Phase C) using the supervisor's unpublished `psalter` package as the reference design. That session ships:

- `tidal/inference/_sphere.py` — cubed-sphere geometry (face indexing, gnomonic cube→sphere, sub-tiles, random rotation).
- `RadialAngularPrior` joint-prior class in `tidal/inference/_prior.py`.
- `tidal sample --joint-prior` per-tile NS mode in `tidal/cli/_sample.py`.
- `tidal plot --type atlas` rendering 2N face panels in `tidal/cli/_plot.py`.
- `docs/V3_PHASE_C_REFERENCE.md` — handoff doc (canonical reference for the cubed-sphere conventions).

These files are **owned by that session** — do not edit from this session.

## Phase tracker (updated as work converges)

See [V3_PHASE_TRACKER.md](V3_PHASE_TRACKER.md).

## Phase E (localised geometry, deferred)

See [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md). Gated on Phase B convergence — provides physical regulator for unbounded P_max in tachyonic regions.

## Phase A-γ (γ_conversion, deferred)

GitHub issues `[v3-A-γ] γ_conversion: refactored definition + multi-t_test + log-zero clamp` and `[v3-A-γ] γ_conversion: P_max correlation validation`. Refactored definition: `γ_conversion = log(|target_amp(t_test)|)/t_test` with multi-t_test minimum-sampling and `log(amp_floor)/t_test` clamp for decoupled targets. Validation: 10-point D1 sample, Spearman ρ vs `log(P_max/P_GR)/t_end`. If ρ > 0.7, ship as opt-in amp metric for ~1000× speedup.

## What's reused from v2

- **The probe** (`tidal/measurement/_stability.py`, profile `unit-ic-all-k-0.15`) stays as the metadata measurement engine — no longer a gate, but γ_eff / k_dominant / n_tachyonic_modes / borderline_stability all flow into chain CSVs for diagnostic and post-hoc filtering.
- **`arctan_uniform` prior** (`tidal/inference/_prior.py:80–84,107–115,138–142`) is the existing compactification primitive — Phase A.2 just rewrites campaign scripts to use it.
- **Stage C truth-table** (`stage_c_truth_table.csv`) and **D1 v2 results** (`hpc_results/28520217/`, `28883112/`, `28896653/`, etc.) are preserved as historical evidence; v3 rerun comparisons frame v2 as "posterior conditioned on stability and perturbativity" and v3 as "posterior over the full coupling space".

## Known limitations (v3.1, to be addressed in v3.2)

The architecture as deployed in Phase B chains carries two architectural compromises that are intentionally retained pending supervisor review. They are documented here so chain interpretation is honest.

1. **Asymmetric per-parameter priors**: v2-inherited `log_uniform:1e-3:1e3` (positive-only) is still used for couplings named α₃ in D1 and mA², ξ, α₃ in Stage A. The v3 architecture should let the chain *learn* whether a parameter must be positive (ghost constraint) from D_KL evidence — but currently we impose positivity by prior fiat. Whether each of these parameters is truly kinetic-coefficient (positivity required to avoid ghosts) or mass-like (sign-symmetric in v3 tachyon-permissive policy) needs theory-by-theory classification. Observed effect in B.4b Stage A sup: chain MAP at α₃ = 0.0012 hits the prior lower bound, indicating posterior support beyond the prior. See [V3_2_DESIGN_INVESTIGATION.md](V3_2_DESIGN_INVESTIGATION.md) for the planned resolution (per-coupling abs-fold in the cubed-sphere joint prior, gated on supervisor classification, with PSALTer integration as the ultimate goal).

2. ~~**Visualization: 95%-credible contours show posterior cliffs as sharp edges**~~. **Resolved (v0.39.1, d784bf1)**: levels extended to `[0.99, 0.95, 0.68]`; the 99% outer ring now shows tail falloff. GH issue #361 closed.

## What's discarded from v2

- The premise that the probe rejects samples (now: the probe records γ_eff for them).
- The Hwang-Noh perturbativity gate.
- Pre-flight numerical-overflow filtering (per supervisor critique: γ·t_end > N is just a precision limit, not a physical signal — running at smaller t_end resolves it).
- Plans to ship a "v3 = consistent-IC-solve" probe path (parallel investigation that died on the v2-blind-spot finding).

## References

- [docs/meetings/2026-05-08_supervisor.md](meetings/2026-05-08_supervisor.md) — origin meeting notes
- [docs/PHASE_6_COMPARISON.md](PHASE_6_COMPARISON.md) — last v2 phase (6.J L=75 bounds-stability)
- [CAMPAIGN.md](../CAMPAIGN.md) — overall campaign log
- [docs/V3_PHASE_TRACKER.md](V3_PHASE_TRACKER.md) — per-phase progress
- [docs/V3_PHASE_D_DESIGN.md](V3_PHASE_D_DESIGN.md) — manuscript update sequence
- [docs/V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md) — localised geometry design (deferred)
- [docs/V3_2_DESIGN_INVESTIGATION.md](V3_2_DESIGN_INVESTIGATION.md) — v3.2 design investigation (symmetric priors + PSALTer)
- `docs/V3_PHASE_C_REFERENCE.md` (parallel session) — cubed-sphere reference
