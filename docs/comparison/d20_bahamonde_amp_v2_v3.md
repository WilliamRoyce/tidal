# d20_bahamonde_amp — v2 vs v3 comparison

> **Correction (2026-08-19, additive — the tables below are preserved as written):**
> the `arctan_uniform` prior columns describe the range as degrees mapped through
> `tan`, giving ±57.3. The sampler never did that: `low`/`high` are ignored
> entirely and the support is fixed at ±tan(π/2 − 0.05) ≈ **±19.98** (GH #425).
> Read every ±57.3 in this file as ±19.98. Marginal D_KL values quoted here also
> pre-date the GH #420 estimator fix — see `docs/RESULTS_AMENDMENTS.md`.

- **v2 reference**: `hpc_results/28598736/`
- **v3 chain**: `hpc_results/29229768/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | +0.616 ± 0.001 | +9.288 ± 0.058 | +8.67 nats |
| ESS | 877 | 8338 | +7461 |
| Joint D_KL | 0.00 nats | 4.31 nats | +4.30 nats |
| n_samples | 3450 | 18124 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| beta1 | uniform[-1.5..1] | arctan[-89°..89°] → ±57.3 | -0.0993 | +2.82 | +0.81 | 0.24 | 2.57 |
| beta2 | uniform[-3..-0.3] | arctan[-89°..89°] → ±57.3 | -1.55 | -0.23 | +0.27 | 0.09 | 2.47 |
| beta3 | uniform[-1..1] | arctan[-89°..89°] → ±57.3 | +0.959 | -0.37 | -0.27 | 0.02 | 2.37 |
| delta1 | uniform[-0.025..0.025] | arctan[-89°..89°] → ±57.3 | -0.0204 | +0.567 | +0.48 | 0.02 | 2.96 |
| xi | log[1e-02..1e+01] | log[1e-03..1e+03] | +1.85 | +0.548 | -0.01 | 0.13 | 0.22 |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
