# d1_amp_joint_smoke — v2 vs v3 comparison

- **v2 reference**: `hpc_results/28896653/`
- **v3 chain**: `hpc_results/29204991/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | +0.720 ± 0.024 | +13.762 ± 0.090 | +13.04 nats |
| ESS | 5997 | 1398 | -4599 |
| Joint D_KL | 0.76 nats | 3.18 nats | +2.41 nats |
| n_samples | 12524 | 3229 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| alpha1 | uniform[-1..1] | joint sphere S^2 × log[1e-03..1e+03] | -0.422 | +582 | +3.80 | 0.02 | 0.34 |
| alpha2 | uniform[-2..2] | joint sphere S^2 × log[1e-03..1e+03] | -0.594 | -575 | -2.43 | 0.02 | 0.28 |
| alpha3 | log[5e-02..2e+00] | log[1e-03..1e+03] | +0.204 | +0.183 | -0.00 | 0.00 | 0.01 |
| delta1 | uniform[-2..2] | joint sphere S^2 × log[1e-03..1e+03] | -0.847 | +575 | +1.86 | 0.10 | 0.27 |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
