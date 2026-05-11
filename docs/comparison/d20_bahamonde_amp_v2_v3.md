# d20_bahamonde_amp — v2 vs v3 comparison

- **v2 reference**: `hpc_results/28598736/`
- **v3 chain**: `hpc_results/29207374/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | +0.616 ± 0.001 | +7.709 ± 0.098 | +7.09 nats |
| ESS | 877 | 2761 | +1884 |
| Joint D_KL | 0.00 nats | 4.35 nats | +4.35 nats |
| n_samples | 3450 | 5675 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| beta1 | uniform[-1.5..1] | arctan[-89°..89°] → ±57.3 | -0.0993 | -0.716 | -0.16 | 0.24 | 2.51 |
| beta2 | uniform[-3..-0.3] | arctan[-89°..89°] → ±57.3 | -1.55 | -1.19 | +0.09 | 0.09 | 2.55 |
| beta3 | uniform[-1..1] | arctan[-89°..89°] → ±57.3 | +0.959 | +1.62 | +0.25 | 0.02 | 2.91 |
| delta1 | uniform[-0.025..0.025] | arctan[-89°..89°] → ±57.3 | -0.0204 | -0.3 | -0.25 | 0.02 | 3.00 |
| xi | log[1e-02..1e+01] | (v2-only) | +1.85 | +nan | +nan | 0.13 | nan |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
