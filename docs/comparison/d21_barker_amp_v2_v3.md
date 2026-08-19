# D2.1 Barker amp — v2 vs v3 comparison

> **Correction (2026-08-19, additive — the tables below are preserved as written):**
> the `arctan_uniform` prior columns describe the range as degrees mapped through
> `tan`, giving ±57.3. The sampler never did that: `low`/`high` are ignored
> entirely and the support is fixed at ±tan(π/2 − 0.05) ≈ **±19.98** (GH #425).
> Read every ±57.3 in this file as ±19.98. Marginal D_KL values quoted here also
> pre-date the GH #420 estimator fix — see `docs/RESULTS_AMENDMENTS.md`.

- **v2 reference**: `hpc_results/d21_barker_amp/`
- **v3 chain**: `hpc_results/d21_barker_amp_v3/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | +0.618 ± 0.001 | +9.581 ± 0.049 | +8.96 nats |
| ESS | 875 | 9517 | +8642 |
| Joint D_KL | 0.00 nats | 4.13 nats | +4.13 nats |
| n_samples | 3458 | 21212 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| beta1 | uniform[-1.5..1] | arctan[-89°..89°] → ±57.3 | +0.734 | -0.796 | -0.44 | 0.30 | 2.60 |
| beta2 | uniform[-3..-0.3] | arctan[-89°..89°] → ±57.3 | -1.91 | -11.7 | -2.54 | 0.07 | 2.54 |
| beta3 | uniform[-1..1] | arctan[-89°..89°] → ±57.3 | -0.671 | +0.585 | +0.45 | 0.02 | 2.68 |
| chi | uniform[-0.009..0.009] | arctan[-89°..89°] → ±57.3 | -0.00207 | +1.75 | +0.46 | 0.03 | 2.60 |
| delta1 | uniform[-0.025..0.025] | arctan[-89°..89°] → ±57.3 | -0.0235 | -1.83 | -0.92 | 0.03 | 2.84 |
| xi | log[1e-02..1e+01] | log[1e-03..1e+03] | +0.0105 | +108 | +0.58 | 0.16 | 0.30 |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
