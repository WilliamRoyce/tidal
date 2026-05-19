# stage_a_amp — v2 vs v3 comparison

- **v2 reference**: `hpc_results/28474676/`
- **v3 chain**: `hpc_results/29189966/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | -0.073 ± 0.007 | +9.307 ± 0.126 | +9.38 nats |
| ESS | 937 | 2649 | +1712 |
| Joint D_KL | 0.02 nats | 6.59 nats | +6.57 nats |
| n_samples | 2457 | 6464 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| alpha3 | log[1e-03..5e-01] | log[1e-03..1e+03] | +0.0015 | +4.46 | +1.10 | 0.03 | 0.95 |
| deltam | uniform[-0.5..0.5] | arctan[-89°..89°] → ±57.3 | +0.477 | +6.91 | +1.95 | 0.05 | 2.46 |
| mA2 | log[1e-03..1e+00] | log[1e-03..1e+03] | +0.646 | +231 | +1.60 | 0.05 | 0.55 |
| xi | log[5e-02..2e+01] | log[1e-03..1e+03] | +0.0571 | +2.47 | +1.76 | 0.26 | 0.59 |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
