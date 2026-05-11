# stage_a_sup — v2 vs v3 comparison

- **v2 reference**: `hpc_results/28477675/`
- **v3 chain**: `hpc_results/29199129/`

## Headline

| Metric | v2 | v3 | Δ |
| --- | --- | --- | --- |
| log Z | +0.661 ± 0.049 | +4.079 ± 0.117 | +3.42 nats |
| ESS | 2724 | 4961 | +2237 |
| Joint D_KL | 1.98 nats | 7.30 nats | +5.32 nats |
| n_samples | 6633 | 9484 | — |

## Per-coupling

| Param | v2 prior | v3 prior | v2 MAP | v3 MAP | MAP shift (v3 σ) | v2 D_KL | v3 D_KL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| alpha3 | log[1e-03..5e-01] | log[1e-03..1e+03] | +0.00109 | +0.00118 | +0.00 | 0.10 | 0.63 |
| deltam | uniform[-0.5..0.5] | arctan[-89°..89°] → ±57.3 | -0.446 | +0.719 | +1.76 | 0.03 | 3.02 |
| mA2 | log[1e-03..1e+00] | log[1e-03..1e+03] | +0.969 | +605 | +2.80 | 0.37 | 0.86 |
| xi | log[5e-02..2e+01] | log[1e-03..1e+03] | +0.796 | +2.07 | +0.67 | 0.41 | 1.05 |

## Notes

* log Z values are not directly comparable across architectures — v2 chains conditioned on the stability gate; v3 integrates over wider compactified support. Use the comparison as a magnitude/direction guide, not a Bayes factor.
* MAP shift in v3-σ units quantifies how far v2's MAP sits inside v3's posterior. |shift| > 2 indicates v3 found a region v2's prior couldn't sample.
* Per-coupling marginal D_KL is the headline metric for v3 — measures how much each coupling carries signal-shape information.
