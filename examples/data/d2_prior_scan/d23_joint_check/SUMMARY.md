# D2.3 Joint Pre-Flight Check — Summary

Scan: random 9-D draws from the composed D2.3 prior, with D2.0 soft constraints applied, then probed via `check_conversion_stability` at the canonical settings (κ=1, B₀=0.01, h_5→a_1, ic_wavevector=2π/L=0.0628…).

Source: `examples/data/d2_prior_scan/d23_joint_check/d23_joint_survival.csv`  
Recipe: `docs/tex/pgt_stability_priors.tex` §D2.3  
Template: `examples/data/d2_prior_scan/zeta_scan/_analyze.py`

## Survival Rates

| Stage | Count | Rate |
|-------|------:|-----:|
| Total draws | 2000 | 100.00% |
| Post β-constraints | 1450 | 72.50% |
| Post stability probe (given β) | 500 | 34.48% |
| **Effective survival** | **500** | **25.00%** |

Reference rates: D2.0/D2.1 prior stability sweep ≈ 60% rejected (40% effective survival); D2.2 sweep 59.1% rejected (40.9% effective).

## Top-5 Surviving Samples by γ_eff

These are the samples that passed both the β-constraints and the stability probe but with the largest max_excess values (closest to the γ_eff > 0.3 rejection threshold). If the top-5 cluster on a particular (param₁, param₂) corner, that signals a near-instability worth a tightening or constraint.

| sample | β₁ | β₂ | β₃ | ξ | δ₁ | χ | ζ₁ | ζ₂ | ζ₃ | γ_eff |
|--------|----|----|----|---|----|---|----|----|----|-------|
| 1972 | -0.332109 | -2.71584 | 0.0168648 | 0.0100433 | -0.0118488 | 0.00222964 | -0.0345195 | 0.00589322 | 0.0429269 | 17.6856 |
| 1132 | -1.09113 | -2.71504 | 0.838479 | 0.0194955 | -0.00305399 | 0.00172623 | 0.0155613 | 0.000698442 | 0.00275594 | 17.5072 |
| 356 | -1.37085 | -1.81128 | -0.713227 | 0.018081 | -0.0155674 | -0.00434243 | 0.0147957 | -0.0075401 | -1.32527e-05 | 17.1629 |
| 1947 | -0.672078 | -2.85968 | 0.12094 | 0.0158268 | 0.0212307 | -0.00373339 | -0.00978328 | -0.00351173 | 0.0453874 | 17.1023 |
| 1167 | -0.372401 | -2.87408 | 0.282858 | 0.0121421 | 0.0192944 | -0.007192 | -0.014928 | 0.00627542 | 0.0230396 | 17.0358 |

## Decision Verdict

**Verdict: PROCEED**

effective survival ≥ 5% → D2.3-B unchanged

Gate thresholds: PROCEED ≥ 5%, PROCEED-WITH-WARN ∈ [1%, 5%), STOP < 1% effective survival.

