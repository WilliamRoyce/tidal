# D1 v1/v2 chain replay — Phase A.5 sanity check

**Created:** 2026-05-10
**Phase:** v3 architecture A.5 (chain replay; precedes Phase B HPC re-runs)
**Companion to:** [docs/V3_ARCHITECTURE.md](V3_ARCHITECTURE.md)

## Headline finding

Replaying the D1 amp v2 chain `28982029` (INTR-reduced, grid=128/nlive=600/τ=0.15) under the v3 architecture admits **14.3% more samples** that v2 rejected via the pre-flight tachyonic probe.

| Category | Count | Fraction | v2 verdict | v3 verdict |
|---|---|---|---|---|
| In-chain success | 5000 | 50.0% | admit | admit (unchanged) |
| Rejected as `success` | 3572 | 35.7% | rejected (reason: not specified in this slice) | unchanged |
| Rejected as `tachyonic` | 1395 (finite γ) + 33 (γ=∞ overflow) = 1428 | 14.3% | reject (-inf via probe gate) | **admit** (probe is metadata only) |

The 1428 tachyonic samples have `γ_eff ∈ [0.15, ∞]` per the probe — by construction they're above the v2 τ=0.15 gate. Under v3 they enter the chain with whatever P_max the simulation produces; their actual contribution to the posterior is determined by the sim, not the probe.

## What this validates

1. **The architecture admits the predicted samples.** v3 doesn't just shuffle existing samples — it widens the chain coverage measurably (50% → 64%).
2. **The probe metadata path works.** All 5000 in-chain success samples carry `tachyonic_excess` metadata in the v2 sidecar; v3 propagates the same data through the chain CSV without gating.
3. **The Hwang-Noh case is unobserved here.** Neither the in-chain nor the prior-rejected samples carry `non_perturbative` status — they were either successful or probe-rejected. The Hwang-Noh removal is verified at the unit-test level (`tests/test_likelihood_no_hwang_noh.py`).

## What this does NOT validate

- The new logL values for the admitted tachyonic samples — those need a real Phase B sim run because the v2 chain CSV doesn't store P_max for prior-rejected samples (the sim was never run for them under v2).
- The soft-floor noise behaviour — verified at unit-test level only.
- The end-to-end PolyChord behaviour with the new architecture — Phase B's smoke run (B.1) is the first integration test.

## Notes on the original 28520217 chain

The earlier D1 v1 chain `hpc_results/28520217/` (run before the rejected-sample sidecar landed in #258) does not have the rejected-prior CSV. All 3514 of its samples are tagged `success`, so the v2-vs-v3 admission comparison can't be made against that chain — it pre-dates the diagnostic infrastructure.

The replay therefore uses 28982029 (v2 INTR-reduced, post-#258 metadata) as the operational reference. The headline numbers from 28520217 (`log Z = -2.26 ± 0.07`, `D_KL = 1.79 nats`) remain the published v1 anchors; v3 will produce numerically distinct headline numbers because it integrates over a different effective posterior support.

## Files

- `v2_admission_shift.png` — two-panel figure: v2 in-chain logL distribution + v2 rejected-tachyonic γ_eff histogram with the τ=0.15 gate marked
- `v2_admission_shift_summary.csv` — admission counts per category
- `v1_logl_distribution.png`, `run_status_summary.csv` — original 28520217 replay artefacts (kept for reference; show that 28520217 pre-dates rejected-sample tracking)

## References

- [docs/V3_ARCHITECTURE.md](V3_ARCHITECTURE.md) — soft-penalty table
- [docs/V3_PHASE_TRACKER.md](V3_PHASE_TRACKER.md) — A.5 status
- GH issue [#349](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/349)
- `hpc_results/28982029/d1_amp_v2_intr_reduced/_rejected_prior.csv` — source data
