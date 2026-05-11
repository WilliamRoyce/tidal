# D1 amp — v2 vs v3 comparison

- **v2 reference**: `hpc_results/28896653/`
- **v3 chain**: `hpc_results/29149987/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | +0.720 ± 0.024 | +13.288 ± 0.125 | +12.57 nats |
| ESS | 5997 | 1410 | -4588 |
| Joint D_KL | 0.76 nats | 4.77 nats | +4.01 nats |
| n_samples | 12524 | 3651 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| alpha1 | uniform[-1..1] | arctan[-89°..89°] → ±57.3 | -0.422 | -0.539 | -0.18 | 0.02 | 3.20 |
| alpha2 | uniform[-2..2] | arctan[-89°..89°] → ±57.3 | -0.594 | -0.65 | -0.05 | 0.02 | 2.96 |
| alpha3 | log[5e-02..2e+00] | log[1e-03..1e+03] | +0.204 | +0.0018 | -0.00 | 0.00 | 0.02 |
| delta1 | uniform[-2..2] | arctan[-89°..89°] → ±57.3 | -0.847 | +4.8 | +0.58 | 0.10 | 1.44 |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
