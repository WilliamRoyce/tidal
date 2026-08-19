# Marginal D_KL Recompute Report (GH #420)

**Status: for review. No document under `manuscript/` or `CAMPAIGN.md` has been
modified — this report is the comparison record; propagation is a separate,
explicitly user-initiated step.**

Generated 2026-08-18 on branch `chore/420-dkl-recompute` after the #420 estimator
fix (PR #426), by re-running `scripts/analysis/recompute_parameter_kl.py` over
`scripts/analysis/d_kl_chains.toml` (chains local in `hpc_results/`, no HPC time).
Pre-fix values are preserved verbatim from the on-disk
`parameter_importance.json` snapshots taken before the re-run.

## Headline

**Per-coupling rankings do NOT survive the fix.** Of the 26 chain-directions
compared (13 chains × amp/sup), **19 change their dominant coupling**, and 220
per-parameter verdicts (STRONG/moderate/weak at 1.0/0.1 nats) flip. The working
hypothesis that the bug inflated absolute values while preserving relative
comparisons is empirically false: the old estimator saturated arctan-prior
marginals near its ceiling log(40) ≈ 3.69 nats, compressing exactly the
differences the rankings relied on, while chains without recorded priors were
scored against fabricated or self-referential references. Any conclusion of the
form "coupling X dominates" or "all couplings inert" that rests on a pre-fix
marginal D_KL needs to be re-checked against the tables below before being
quoted.

Numbers that survive unchanged: joint D_KL, log Z, d_G (computed by anesthetic,
never affected), all `uniform`/`log_uniform` marginals (bit-identical by
construction of the fix, verified: ricci_em beta3 0.0049 → 0.0049), and the
cross amp-vs-sup KL (a different estimator, not prior-referenced).

## The defects (what was wrong, in brief)

1. **arctan_uniform** marginals were histogrammed in linear space against a
   uniform reference over the recorded bounds (±89) — bounds the sampler never
   used (support is fixed at ±tan(π/2 − 0.05) ≈ ±20, see GH #425). Null
   posteriors reported ~2.5 nats; constrained ones saturated near log(40) ≈ 3.69.
2. **normal** priors used `range=(mean, std)` as a histogram range (garbage or
   NaN). No campaign chain used a normal prior; listed for completeness.
3. **radial_angular** records were silently skipped by the prior parser (no
   `name`/`low`/`high` keys), so every cubed-sphere coupling was scored against
   a uniform reference over its own sample range (affects the atlas runs, which
   are outside this inventory).
4. An `or`-chain dropped object-path priors with `low == 0.0` into the same
   fallback (falsy zero).
5. **Recompute-script rename mismatch** (fixed on this branch):
   `d_kl_chains.toml` relabels chain columns for some chains (ricci_em
   `alpha*` → `beta*`), and the name-based prior lookup missed every renamed
   column. Masked pre-fix because the saturated estimator produced
   near-identical wrong numbers either way.
6. **Fabricated priors for pre-schema chains** (mitigated on this branch): six
   chain dirs have no `inference.json`, so the script fabricated
   all-arctan(±89) priors. The TOML header promised a `[chain.priors]`
   override for this case, but the script never read it. The override is now
   implemented (`prior_overrides`), and the true priors were recovered from
   the campaign sbatch templates: T7 (chi_closure) and T6 (parity_odd_full)
   actually used `xi = log_uniform:1e-3:1e3` — under the fabricated arctan
   prior, xi showed a spurious dominant 1.5–1.7 nats in T7; with the correct
   prior it drops to noise level.

The fixed estimator also self-checks from now on: importance.json carries a
`consistency` block — superadditivity (for product priors
`sum(marginals) <= joint D_KL` exactly; the pre-fix outputs violated it) with an
N_eff-scaled bias allowance, histogram-ceiling saturation flags, and the Kish
effective sample size `n_eff` — surfaced as WARNING lines in the CLI table.

## How to read the tables

- **Verdicts** use `format_importance_table` bands: STRONG > 1.0, moderate >
  0.1, weak otherwise. Campaign thresholds (`docs/campaign_plan.md:140-149`)
  additionally key on 0.05 / 0.005 / 0.1 nats for stage pass/fail.
- **N_eff and the noise floor**: the histogram marginal estimator is biased up
  by ≈ (n_bins − 1)/(2·N_eff) per parameter. For healthy chains
  (N_eff ≳ 2000) that is < 0.01 nats; for the timeout/floor-contaminated
  rescue chains (T9 amp: N_eff ≈ 21 → floor ≈ 0.9 nats) **every marginal is
  floor, not signal** — this fully accounts for the residual raw
  sum-of-marginals > joint in those rows and matches their already-recorded
  "overflow-contaminated" status in `docs/V3_PHASE_TRACKER.md`.

## Dominant-coupling summary

| chain | dir | dominant: old → new | raw sum>joint: old → new | N_eff |
|---|---|---|---|---|
| dark_photon_plasma | amp | deltam → xi ⚠ | no → no | 7443 |
| dark_photon_plasma | sup | deltam → deltam | no → no | 10675 |
| ricci_em | amp | beta1 → delta1 ⚠ | yes → no | 5694 |
| ricci_em | sup | delta1 → delta1 | no → no | 5224 |
| ym_pgt_bahamonde | amp | delta1 → xi ⚠ | yes → no | 1063 |
| ym_pgt_bahamonde | sup | delta1 → xi ⚠ | yes → no | 2404 |
| ym_pgt_barker | amp | delta1 → chi ⚠ | yes → no | 1340 |
| ym_pgt_barker | sup | chi → chi | yes → no | 2895 |
| ym_pgt_shapiro | amp | zeta2 → zeta2 | yes → no | 1341 |
| ym_pgt_shapiro | sup | beta1 → zeta1 ⚠ | yes → no | 3302 |
| ym_pgt_full_nonminimal | amp | delta1 → xi ⚠ | yes → no | 1297 |
| ym_pgt_full_nonminimal | sup | chi → chi | yes → no | 2870 |
| np_t5 | amp | zeta2 → beta1 ⚠ | yes → no | 1146 |
| np_t5 | sup | chi → zeta2 ⚠ | yes → no | 2867 |
| chi_closure | amp | chi5 → chi1 ⚠ | yes → yes | 118 |
| chi_closure | sup | chi1 → chi1 | yes → no | 248 |
| np_chi_closure | amp | chi10 → chi1 ⚠ | yes → yes | 80 |
| np_chi_closure | sup | chi1 → chi1 | yes → no | 211 |
| xi_kinetic_closure | amp | chi8 → chi1 ⚠ | yes → yes | 21 |
| xi_kinetic_closure | sup | chi7 → chi1 ⚠ | yes → yes | 102 |
| parity_odd_minimal | amp | beta2 → d21 ⚠ | yes → no | 944 |
| parity_odd_minimal | sup | beta2 → beta1 ⚠ | no → no | 642 |
| parity_odd_full | amp | xi → beta3 ⚠ | yes → yes | 6 |
| parity_odd_full | sup | d21 → beta2 ⚠ | yes → yes | 14 |
| eh_gertsenshtein | amp | rho → sigma ⚠ | yes → yes | 104 |
| eh_gertsenshtein | sup | rho → sigma ⚠ | yes → yes | 103 |

## Per-chain comparison

### dark_photon_plasma (29205968+29205982, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 7.015; N_eff ≈ 7443 (per-param noise floor ≈ 0.00 nats); sum of marginals 4.776 → 2.179 (raw superadditivity ok → ok)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| mA2 | log_uniform | 0.505 | 0.505 | +0.000 | moderate → moderate |
| deltam | arctan_uniform | 2.767 | 0.169 | -2.597 | STRONG → moderate **FLIP** |
| xi | log_uniform | 0.848 | 0.848 | +0.000 | moderate → moderate |
| alpha3 | log_uniform | 0.657 | 0.657 | +0.000 | moderate → moderate |

Ranking (most→least constrained): old `deltam > xi > alpha3 > mA2` → new `xi > alpha3 > mA2 > deltam` — **dominant coupling changed**

**sup** — joint D_KL 7.226; N_eff ≈ 10675 (per-param noise floor ≈ 0.00 nats); sum of marginals 6.144 → 4.617 (raw superadditivity ok → ok)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| mA2 | log_uniform | 1.250 | 1.250 | +0.000 | STRONG → STRONG |
| deltam | arctan_uniform | 2.983 | 1.455 | -1.527 | STRONG → STRONG |
| xi | log_uniform | 1.238 | 1.238 | +0.000 | STRONG → STRONG |
| alpha3 | log_uniform | 0.673 | 0.673 | +0.000 | moderate → moderate |

Ranking (most→least constrained): old `deltam > mA2 > xi > alpha3` → new `deltam > mA2 > xi > alpha3` — order preserved at the top

### ricci_em (29189748+29189761, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 4.645; N_eff ≈ 5694 (per-param noise floor ≈ 0.00 nats); sum of marginals 7.557 → 2.533 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 3.179 | 0.522 | -2.657 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.944 | 0.145 | -2.799 | STRONG → moderate **FLIP** |
| beta3 | log_uniform | 0.005 | 0.005 | +0.000 | weak → weak |
| delta1 | arctan_uniform | 1.430 | 1.861 | +0.431 | STRONG → STRONG |

Ranking (most→least constrained): old `beta1 > beta2 > delta1 > beta3` → new `delta1 > beta1 > beta2 > beta3` — **dominant coupling changed**

**sup** — joint D_KL 13.159; N_eff ≈ 5224 (per-param noise floor ≈ 0.00 nats); sum of marginals 8.791 → 1.450 (raw superadditivity ok → ok)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.870 | 0.270 | -2.600 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.770 | 0.258 | -2.511 | STRONG → moderate **FLIP** |
| beta3 | log_uniform | 0.147 | 0.147 | +0.000 | moderate → moderate |
| delta1 | arctan_uniform | 3.004 | 0.774 | -2.230 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `delta1 > beta1 > beta2 > beta3` → new `delta1 > beta1 > beta2 > beta3` — order preserved at the top

### ym_pgt_bahamonde (29507332, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 6.015; N_eff ≈ 1063 (per-param noise floor ≈ 0.02 nats); sum of marginals 12.636 → 3.050 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.934 | 0.595 | -2.339 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.667 | 0.313 | -2.353 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.817 | 0.185 | -2.632 | STRONG → moderate **FLIP** |
| xi | log_uniform | 1.201 | 1.201 | +0.000 | STRONG → STRONG |
| delta1 | arctan_uniform | 3.017 | 0.755 | -2.262 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `delta1 > beta1 > beta3 > beta2 > xi` → new `xi > delta1 > beta1 > beta2 > beta3` — **dominant coupling changed**

**sup** — joint D_KL 5.832; N_eff ≈ 2404 (per-param noise floor ≈ 0.01 nats); sum of marginals 13.133 → 2.808 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.930 | 0.275 | -2.655 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.878 | 0.116 | -2.762 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.500 | 0.030 | -2.470 | STRONG → weak **FLIP** |
| xi | log_uniform | 1.767 | 1.767 | +0.000 | STRONG → STRONG |
| delta1 | arctan_uniform | 3.057 | 0.619 | -2.437 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `delta1 > beta1 > beta2 > beta3 > xi` → new `xi > delta1 > beta1 > beta2 > beta3` — **dominant coupling changed**

### ym_pgt_barker (29507332, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 8.144; N_eff ≈ 1340 (per-param noise floor ≈ 0.01 nats); sum of marginals 16.739 → 3.479 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.364 | 0.419 | -1.945 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.664 | 0.108 | -2.556 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.973 | 0.246 | -2.727 | STRONG → moderate **FLIP** |
| xi | arctan_uniform | 2.647 | 0.920 | -1.727 | STRONG → moderate **FLIP** |
| delta1 | arctan_uniform | 3.091 | 0.764 | -2.327 | STRONG → moderate **FLIP** |
| chi | arctan_uniform | 3.000 | 1.022 | -1.978 | STRONG → STRONG |

Ranking (most→least constrained): old `delta1 > chi > beta3 > beta2 > xi > beta1` → new `chi > xi > delta1 > beta1 > beta3 > beta2` — **dominant coupling changed**

**sup** — joint D_KL 6.553; N_eff ≈ 2895 (per-param noise floor ≈ 0.01 nats); sum of marginals 15.717 → 4.362 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.832 | 0.160 | -2.671 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.653 | 0.076 | -2.577 | STRONG → weak **FLIP** |
| beta3 | arctan_uniform | 2.546 | 0.008 | -2.538 | STRONG → weak **FLIP** |
| xi | arctan_uniform | 1.592 | 1.592 | +0.000 | STRONG → STRONG |
| delta1 | arctan_uniform | 2.982 | 0.606 | -2.376 | STRONG → moderate **FLIP** |
| chi | arctan_uniform | 3.112 | 1.919 | -1.193 | STRONG → STRONG |

Ranking (most→least constrained): old `chi > delta1 > beta1 > beta2 > beta3 > xi` → new `chi > xi > delta1 > beta1 > beta2 > beta3` — order preserved at the top

### ym_pgt_shapiro (29468763, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 4.541; N_eff ≈ 1341 (per-param noise floor ≈ 0.01 nats); sum of marginals 18.751 → 0.829 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.370 | 0.141 | -2.229 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.484 | 0.028 | -2.456 | STRONG → weak **FLIP** |
| beta3 | arctan_uniform | 2.616 | 0.097 | -2.519 | STRONG → weak **FLIP** |
| xi | log_uniform | 0.152 | 0.152 | +0.000 | moderate → moderate |
| delta1 | arctan_uniform | 2.867 | 0.119 | -2.748 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.674 | 0.072 | -2.601 | STRONG → weak **FLIP** |
| zeta2 | arctan_uniform | 2.962 | 0.154 | -2.808 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.626 | 0.065 | -2.561 | STRONG → weak **FLIP** |

Ranking (most→least constrained): old `zeta2 > delta1 > zeta1 > zeta3 > beta3 > beta2 > beta1 > xi` → new `zeta2 > xi > beta1 > delta1 > beta3 > zeta1 > zeta3 > beta2` — order preserved at the top

**sup** — joint D_KL 6.720; N_eff ≈ 3302 (per-param noise floor ≈ 0.01 nats); sum of marginals 19.419 → 2.239 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 3.035 | 0.403 | -2.632 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.829 | 0.139 | -2.690 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.743 | 0.108 | -2.635 | STRONG → moderate **FLIP** |
| xi | log_uniform | 0.440 | 0.440 | +0.000 | moderate → moderate |
| delta1 | arctan_uniform | 2.940 | 0.209 | -2.732 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 1.657 | 0.583 | -1.074 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.969 | 0.265 | -2.705 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.805 | 0.092 | -2.713 | STRONG → weak **FLIP** |

Ranking (most→least constrained): old `beta1 > zeta2 > delta1 > beta2 > zeta3 > beta3 > zeta1 > xi` → new `zeta1 > xi > beta1 > zeta2 > delta1 > beta2 > beta3 > zeta3` — **dominant coupling changed**

### ym_pgt_full_nonminimal (29468763, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 3.916; N_eff ≈ 1297 (per-param noise floor ≈ 0.02 nats); sum of marginals 21.219 → 0.853 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.558 | 0.052 | -2.506 | STRONG → weak **FLIP** |
| beta2 | arctan_uniform | 2.466 | 0.036 | -2.430 | STRONG → weak **FLIP** |
| beta3 | arctan_uniform | 2.630 | 0.046 | -2.584 | STRONG → weak **FLIP** |
| xi | log_uniform | 0.180 | 0.180 | +0.000 | moderate → moderate |
| delta1 | arctan_uniform | 2.926 | 0.111 | -2.815 | STRONG → moderate **FLIP** |
| chi | arctan_uniform | 2.307 | 0.153 | -2.154 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.630 | 0.096 | -2.534 | STRONG → weak **FLIP** |
| zeta2 | arctan_uniform | 2.810 | 0.105 | -2.705 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.711 | 0.074 | -2.637 | STRONG → weak **FLIP** |

Ranking (most→least constrained): old `delta1 > zeta2 > zeta3 > zeta1 > beta3 > beta1 > beta2 > chi > xi` → new `xi > chi > delta1 > zeta2 > zeta1 > zeta3 > beta1 > beta3 > beta2` — **dominant coupling changed**

**sup** — joint D_KL 6.486; N_eff ≈ 2870 (per-param noise floor ≈ 0.01 nats); sum of marginals 21.619 → 2.806 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.795 | 0.411 | -2.384 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.638 | 0.182 | -2.456 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.612 | 0.032 | -2.580 | STRONG → weak **FLIP** |
| xi | log_uniform | 0.311 | 0.311 | +0.000 | moderate → moderate |
| delta1 | arctan_uniform | 2.711 | 0.060 | -2.651 | STRONG → weak **FLIP** |
| chi | arctan_uniform | 3.011 | 1.434 | -1.576 | STRONG → STRONG |
| zeta1 | arctan_uniform | 1.963 | 0.238 | -1.725 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.858 | 0.086 | -2.772 | STRONG → weak **FLIP** |
| zeta3 | arctan_uniform | 2.720 | 0.051 | -2.670 | STRONG → weak **FLIP** |

Ranking (most→least constrained): old `chi > zeta2 > beta1 > zeta3 > delta1 > beta2 > beta3 > zeta1 > xi` → new `chi > beta1 > xi > zeta1 > beta2 > zeta2 > delta1 > zeta3 > beta3` — order preserved at the top

### np_t5 (29700462, plane_wave)

Priors: **no inference.json** — reconstructed from the campaign sbatch template (all arctan ±89). Verify against scripts/hpc_templates before quoting.

**amp** — joint D_KL 6.204; N_eff ≈ 1146 (per-param noise floor ≈ 0.02 nats); sum of marginals 21.154 → 2.548 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.189 | 0.640 | -1.549 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.352 | 0.250 | -2.102 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.467 | 0.163 | -2.304 | STRONG → moderate **FLIP** |
| delta1 | arctan_uniform | 2.763 | 0.189 | -2.574 | STRONG → moderate **FLIP** |
| chi | arctan_uniform | 2.812 | 0.350 | -2.462 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.820 | 0.310 | -2.510 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.989 | 0.383 | -2.606 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.763 | 0.264 | -2.499 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `zeta2 > zeta1 > chi > delta1 > zeta3 > beta3 > beta2 > beta1` → new `beta1 > zeta2 > chi > zeta1 > zeta3 > beta2 > delta1 > beta3` — **dominant coupling changed**

**sup** — joint D_KL 4.351; N_eff ≈ 2867 (per-param noise floor ≈ 0.01 nats); sum of marginals 20.461 → 0.809 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.676 | 0.086 | -2.589 | STRONG → weak **FLIP** |
| beta2 | arctan_uniform | 2.539 | 0.036 | -2.503 | STRONG → weak **FLIP** |
| beta3 | arctan_uniform | 2.750 | 0.084 | -2.666 | STRONG → weak **FLIP** |
| delta1 | arctan_uniform | 2.734 | 0.066 | -2.668 | STRONG → weak **FLIP** |
| chi | arctan_uniform | 2.893 | 0.166 | -2.727 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.556 | 0.044 | -2.511 | STRONG → weak **FLIP** |
| zeta2 | arctan_uniform | 2.055 | 0.219 | -1.835 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.258 | 0.106 | -2.152 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `chi > beta3 > delta1 > beta1 > zeta1 > beta2 > zeta3 > zeta2` → new `zeta2 > chi > zeta3 > beta1 > beta3 > delta1 > zeta1 > beta2` — **dominant coupling changed**

### chi_closure (29682868, plane_wave)

Priors: **no inference.json** — reconstructed from the campaign sbatch template (all arctan ±89; overrides: {'xi': 'log_uniform:1e-3:1e3'}). Verify against scripts/hpc_templates before quoting.

**amp** — joint D_KL 3.749; N_eff ≈ 118 (per-param noise floor ≈ 0.17 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 44.527 → 3.778 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.594 | 0.222 | -2.372 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.666 | 0.176 | -2.491 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.677 | 0.162 | -2.515 | STRONG → moderate **FLIP** |
| xi | arctan_uniform | 2.447 | 0.180 | -2.267 | STRONG → moderate **FLIP** |
| delta1 | arctan_uniform | 2.550 | 0.203 | -2.347 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.720 | 0.202 | -2.517 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.505 | 0.119 | -2.386 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.372 | 0.191 | -2.181 | STRONG → moderate **FLIP** |
| chi1 | arctan_uniform | 1.677 | 0.670 | -1.007 | STRONG → moderate **FLIP** |
| chi2 | arctan_uniform | 2.374 | 0.217 | -2.157 | STRONG → moderate **FLIP** |
| chi3 | arctan_uniform | 2.543 | 0.183 | -2.360 | STRONG → moderate **FLIP** |
| chi4 | arctan_uniform | 2.398 | 0.157 | -2.241 | STRONG → moderate **FLIP** |
| chi5 | arctan_uniform | 2.741 | 0.148 | -2.593 | STRONG → moderate **FLIP** |
| chi6 | arctan_uniform | 2.505 | 0.177 | -2.327 | STRONG → moderate **FLIP** |
| chi7 | arctan_uniform | 2.375 | 0.220 | -2.155 | STRONG → moderate **FLIP** |
| chi8 | arctan_uniform | 2.386 | 0.223 | -2.164 | STRONG → moderate **FLIP** |
| chi9 | arctan_uniform | 2.439 | 0.110 | -2.329 | STRONG → moderate **FLIP** |
| chi10 | arctan_uniform | 2.557 | 0.216 | -2.341 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `chi5 > zeta1 > beta3 > beta2 > beta1 > chi10 > delta1 > chi3 > zeta2 > chi6 > xi > chi9 > chi4 > chi8 > chi7 > chi2 > zeta3 > chi1` → new `chi1 > chi8 > beta1 > chi7 > chi2 > chi10 > delta1 > zeta1 > zeta3 > chi3 > xi > chi6 > beta2 > beta3 > chi4 > chi5 > zeta2 > chi9` — **dominant coupling changed**

**sup** — joint D_KL 2.579; N_eff ≈ 248 (per-param noise floor ≈ 0.08 nats); sum of marginals 46.701 → 1.458 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.733 | 0.098 | -2.634 | STRONG → weak **FLIP** |
| beta2 | arctan_uniform | 2.725 | 0.086 | -2.639 | STRONG → weak **FLIP** |
| beta3 | arctan_uniform | 2.608 | 0.053 | -2.555 | STRONG → weak **FLIP** |
| xi | arctan_uniform | 2.382 | 0.089 | -2.294 | STRONG → weak **FLIP** |
| delta1 | arctan_uniform | 2.528 | 0.066 | -2.461 | STRONG → weak **FLIP** |
| zeta1 | arctan_uniform | 2.598 | 0.075 | -2.523 | STRONG → weak **FLIP** |
| zeta2 | arctan_uniform | 2.379 | 0.081 | -2.298 | STRONG → weak **FLIP** |
| zeta3 | arctan_uniform | 2.428 | 0.079 | -2.349 | STRONG → weak **FLIP** |
| chi1 | arctan_uniform | 2.822 | 0.123 | -2.700 | STRONG → moderate **FLIP** |
| chi2 | arctan_uniform | 2.659 | 0.060 | -2.598 | STRONG → weak **FLIP** |
| chi3 | arctan_uniform | 2.527 | 0.077 | -2.451 | STRONG → weak **FLIP** |
| chi4 | arctan_uniform | 2.573 | 0.074 | -2.498 | STRONG → weak **FLIP** |
| chi5 | arctan_uniform | 2.552 | 0.062 | -2.491 | STRONG → weak **FLIP** |
| chi6 | arctan_uniform | 2.684 | 0.077 | -2.607 | STRONG → weak **FLIP** |
| chi7 | arctan_uniform | 2.810 | 0.098 | -2.712 | STRONG → weak **FLIP** |
| chi8 | arctan_uniform | 2.559 | 0.090 | -2.469 | STRONG → weak **FLIP** |
| chi9 | arctan_uniform | 2.703 | 0.085 | -2.618 | STRONG → weak **FLIP** |
| chi10 | arctan_uniform | 2.431 | 0.085 | -2.346 | STRONG → weak **FLIP** |

Ranking (most→least constrained): old `chi1 > chi7 > beta1 > beta2 > chi9 > chi6 > chi2 > beta3 > zeta1 > chi4 > chi8 > chi5 > delta1 > chi3 > chi10 > zeta3 > xi > zeta2` → new `chi1 > beta1 > chi7 > chi8 > xi > beta2 > chi10 > chi9 > zeta2 > zeta3 > chi6 > chi3 > zeta1 > chi4 > delta1 > chi5 > chi2 > beta3` — order preserved at the top

### np_chi_closure (29705560, plane_wave)

Priors: **no inference.json** — reconstructed from the campaign sbatch template (all arctan ±89). Verify against scripts/hpc_templates before quoting.

**amp** — joint D_KL 3.918; N_eff ≈ 80 (per-param noise floor ≈ 0.24 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 42.123 → 4.504 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.600 | 0.220 | -2.380 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.503 | 0.291 | -2.211 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.467 | 0.294 | -2.173 | STRONG → moderate **FLIP** |
| delta1 | arctan_uniform | 2.477 | 0.294 | -2.183 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.681 | 0.204 | -2.477 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.621 | 0.223 | -2.398 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.301 | 0.261 | -2.040 | STRONG → moderate **FLIP** |
| chi1 | arctan_uniform | 1.720 | 0.715 | -1.005 | STRONG → moderate **FLIP** |
| chi2 | arctan_uniform | 2.394 | 0.287 | -2.107 | STRONG → moderate **FLIP** |
| chi3 | arctan_uniform | 2.578 | 0.200 | -2.378 | STRONG → moderate **FLIP** |
| chi4 | arctan_uniform | 2.503 | 0.203 | -2.300 | STRONG → moderate **FLIP** |
| chi5 | arctan_uniform | 2.680 | 0.186 | -2.494 | STRONG → moderate **FLIP** |
| chi6 | arctan_uniform | 2.538 | 0.167 | -2.371 | STRONG → moderate **FLIP** |
| chi7 | arctan_uniform | 2.366 | 0.245 | -2.121 | STRONG → moderate **FLIP** |
| chi8 | arctan_uniform | 2.420 | 0.258 | -2.162 | STRONG → moderate **FLIP** |
| chi9 | arctan_uniform | 2.553 | 0.206 | -2.348 | STRONG → moderate **FLIP** |
| chi10 | arctan_uniform | 2.722 | 0.252 | -2.470 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `chi10 > zeta1 > chi5 > zeta2 > beta1 > chi3 > chi9 > chi6 > chi4 > beta2 > delta1 > beta3 > chi8 > chi2 > chi7 > zeta3 > chi1` → new `chi1 > delta1 > beta3 > beta2 > chi2 > zeta3 > chi8 > chi10 > chi7 > zeta2 > beta1 > chi9 > zeta1 > chi4 > chi3 > chi5 > chi6` — **dominant coupling changed**

**sup** — joint D_KL 2.653; N_eff ≈ 211 (per-param noise floor ≈ 0.09 nats); sum of marginals 44.293 → 1.583 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.763 | 0.102 | -2.661 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.681 | 0.095 | -2.586 | STRONG → weak **FLIP** |
| beta3 | arctan_uniform | 2.542 | 0.076 | -2.465 | STRONG → weak **FLIP** |
| delta1 | arctan_uniform | 2.460 | 0.086 | -2.374 | STRONG → weak **FLIP** |
| zeta1 | arctan_uniform | 2.508 | 0.095 | -2.412 | STRONG → weak **FLIP** |
| zeta2 | arctan_uniform | 2.514 | 0.089 | -2.424 | STRONG → weak **FLIP** |
| zeta3 | arctan_uniform | 2.473 | 0.076 | -2.396 | STRONG → weak **FLIP** |
| chi1 | arctan_uniform | 2.820 | 0.145 | -2.675 | STRONG → moderate **FLIP** |
| chi2 | arctan_uniform | 2.777 | 0.098 | -2.679 | STRONG → weak **FLIP** |
| chi3 | arctan_uniform | 2.448 | 0.122 | -2.325 | STRONG → moderate **FLIP** |
| chi4 | arctan_uniform | 2.611 | 0.062 | -2.548 | STRONG → weak **FLIP** |
| chi5 | arctan_uniform | 2.627 | 0.098 | -2.529 | STRONG → weak **FLIP** |
| chi6 | arctan_uniform | 2.610 | 0.075 | -2.534 | STRONG → weak **FLIP** |
| chi7 | arctan_uniform | 2.796 | 0.111 | -2.684 | STRONG → moderate **FLIP** |
| chi8 | arctan_uniform | 2.425 | 0.086 | -2.340 | STRONG → weak **FLIP** |
| chi9 | arctan_uniform | 2.698 | 0.091 | -2.606 | STRONG → weak **FLIP** |
| chi10 | arctan_uniform | 2.543 | 0.073 | -2.470 | STRONG → weak **FLIP** |

Ranking (most→least constrained): old `chi1 > chi7 > chi2 > beta1 > chi9 > beta2 > chi5 > chi4 > chi6 > chi10 > beta3 > zeta2 > zeta1 > zeta3 > delta1 > chi3 > chi8` → new `chi1 > chi3 > chi7 > beta1 > chi5 > chi2 > zeta1 > beta2 > chi9 > zeta2 > delta1 > chi8 > beta3 > zeta3 > chi6 > chi10 > chi4` — order preserved at the top

### xi_kinetic_closure (29694142, plane_wave)

Priors: **no inference.json** — reconstructed from the campaign sbatch template (all arctan ±89). Verify against scripts/hpc_templates before quoting.

**amp** — joint D_KL 4.129; N_eff ≈ 21 (per-param noise floor ≈ 0.94 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 83.767 → 24.253 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.887 | 0.770 | -2.117 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.753 | 0.600 | -2.153 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.592 | 0.906 | -1.686 | STRONG → moderate **FLIP** |
| xi1 | arctan_uniform | 2.221 | 0.794 | -1.426 | STRONG → moderate **FLIP** |
| xi2 | arctan_uniform | 2.631 | 0.597 | -2.033 | STRONG → moderate **FLIP** |
| xi3 | arctan_uniform | 2.650 | 0.743 | -1.907 | STRONG → moderate **FLIP** |
| xi4 | arctan_uniform | 2.551 | 0.626 | -1.925 | STRONG → moderate **FLIP** |
| xi5 | arctan_uniform | 2.721 | 0.557 | -2.164 | STRONG → moderate **FLIP** |
| xi6 | arctan_uniform | 2.810 | 0.772 | -2.038 | STRONG → moderate **FLIP** |
| xi7 | arctan_uniform | 2.814 | 0.791 | -2.023 | STRONG → moderate **FLIP** |
| xi8 | arctan_uniform | 2.868 | 0.818 | -2.051 | STRONG → moderate **FLIP** |
| xi9 | arctan_uniform | 2.861 | 0.728 | -2.132 | STRONG → moderate **FLIP** |
| xi10 | arctan_uniform | 2.428 | 0.755 | -1.673 | STRONG → moderate **FLIP** |
| xi12 | arctan_uniform | 2.478 | 0.684 | -1.793 | STRONG → moderate **FLIP** |
| xi13 | arctan_uniform | 2.507 | 0.703 | -1.804 | STRONG → moderate **FLIP** |
| xi14 | arctan_uniform | 2.865 | 0.612 | -2.252 | STRONG → moderate **FLIP** |
| xi15 | arctan_uniform | 2.474 | 0.869 | -1.605 | STRONG → moderate **FLIP** |
| xi16 | arctan_uniform | 2.464 | 0.973 | -1.490 | STRONG → moderate **FLIP** |
| delta1 | arctan_uniform | 2.658 | 0.748 | -1.910 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.792 | 0.691 | -2.101 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.658 | 0.772 | -1.886 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.643 | 0.765 | -1.877 | STRONG → moderate **FLIP** |
| chi1 | arctan_uniform | 2.019 | 1.186 | -0.833 | STRONG → STRONG |
| chi2 | arctan_uniform | 2.586 | 0.766 | -1.820 | STRONG → moderate **FLIP** |
| chi3 | arctan_uniform | 2.809 | 0.704 | -2.105 | STRONG → moderate **FLIP** |
| chi4 | arctan_uniform | 2.638 | 0.856 | -1.781 | STRONG → moderate **FLIP** |
| chi5 | arctan_uniform | 2.330 | 0.768 | -1.562 | STRONG → moderate **FLIP** |
| chi6 | arctan_uniform | 2.576 | 0.763 | -1.813 | STRONG → moderate **FLIP** |
| chi7 | arctan_uniform | 2.393 | 0.663 | -1.730 | STRONG → moderate **FLIP** |
| chi8 | arctan_uniform | 2.924 | 0.712 | -2.211 | STRONG → moderate **FLIP** |
| chi9 | arctan_uniform | 2.574 | 0.966 | -1.608 | STRONG → moderate **FLIP** |
| chi10 | arctan_uniform | 2.596 | 0.593 | -2.003 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `chi8 > beta1 > xi8 > xi14 > xi9 > xi7 > xi6 > chi3 > zeta1 > beta2 > xi5 > delta1 > zeta2 > xi3 > zeta3 > chi4 > xi2 > chi10 > beta3 > chi2 > chi6 > chi9 > xi4 > xi13 > xi12 > xi15 > xi16 > xi10 > chi7 > chi5 > xi1 > chi1` → new `chi1 > xi16 > chi9 > beta3 > xi15 > chi4 > xi8 > xi1 > xi7 > xi6 > zeta2 > beta1 > chi5 > chi2 > zeta3 > chi6 > xi10 > delta1 > xi3 > xi9 > chi8 > chi3 > xi13 > zeta1 > xi12 > chi7 > xi4 > xi14 > beta2 > xi2 > chi10 > xi5` — **dominant coupling changed**

**sup** — joint D_KL 2.240; N_eff ≈ 102 (per-param noise floor ≈ 0.19 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 82.799 → 5.421 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.627 | 0.200 | -2.428 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.669 | 0.223 | -2.446 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.694 | 0.158 | -2.536 | STRONG → moderate **FLIP** |
| xi1 | arctan_uniform | 2.486 | 0.189 | -2.298 | STRONG → moderate **FLIP** |
| xi2 | arctan_uniform | 2.575 | 0.181 | -2.394 | STRONG → moderate **FLIP** |
| xi3 | arctan_uniform | 2.600 | 0.159 | -2.442 | STRONG → moderate **FLIP** |
| xi4 | arctan_uniform | 2.574 | 0.125 | -2.449 | STRONG → moderate **FLIP** |
| xi5 | arctan_uniform | 2.603 | 0.196 | -2.407 | STRONG → moderate **FLIP** |
| xi6 | arctan_uniform | 2.537 | 0.175 | -2.363 | STRONG → moderate **FLIP** |
| xi7 | arctan_uniform | 2.634 | 0.147 | -2.487 | STRONG → moderate **FLIP** |
| xi8 | arctan_uniform | 2.295 | 0.202 | -2.093 | STRONG → moderate **FLIP** |
| xi9 | arctan_uniform | 2.467 | 0.144 | -2.324 | STRONG → moderate **FLIP** |
| xi10 | arctan_uniform | 2.681 | 0.181 | -2.500 | STRONG → moderate **FLIP** |
| xi12 | arctan_uniform | 2.452 | 0.120 | -2.331 | STRONG → moderate **FLIP** |
| xi13 | arctan_uniform | 2.496 | 0.152 | -2.344 | STRONG → moderate **FLIP** |
| xi14 | arctan_uniform | 2.448 | 0.154 | -2.295 | STRONG → moderate **FLIP** |
| xi15 | arctan_uniform | 2.588 | 0.193 | -2.395 | STRONG → moderate **FLIP** |
| xi16 | arctan_uniform | 2.557 | 0.149 | -2.409 | STRONG → moderate **FLIP** |
| delta1 | arctan_uniform | 2.602 | 0.170 | -2.432 | STRONG → moderate **FLIP** |
| zeta1 | arctan_uniform | 2.594 | 0.162 | -2.432 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.412 | 0.168 | -2.244 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.485 | 0.154 | -2.331 | STRONG → moderate **FLIP** |
| chi1 | arctan_uniform | 2.726 | 0.231 | -2.495 | STRONG → moderate **FLIP** |
| chi2 | arctan_uniform | 2.781 | 0.177 | -2.604 | STRONG → moderate **FLIP** |
| chi3 | arctan_uniform | 2.540 | 0.173 | -2.367 | STRONG → moderate **FLIP** |
| chi4 | arctan_uniform | 2.649 | 0.132 | -2.518 | STRONG → moderate **FLIP** |
| chi5 | arctan_uniform | 2.760 | 0.186 | -2.573 | STRONG → moderate **FLIP** |
| chi6 | arctan_uniform | 2.665 | 0.175 | -2.490 | STRONG → moderate **FLIP** |
| chi7 | arctan_uniform | 2.914 | 0.171 | -2.744 | STRONG → moderate **FLIP** |
| chi8 | arctan_uniform | 2.583 | 0.158 | -2.425 | STRONG → moderate **FLIP** |
| chi9 | arctan_uniform | 2.537 | 0.146 | -2.391 | STRONG → moderate **FLIP** |
| chi10 | arctan_uniform | 2.567 | 0.174 | -2.393 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `chi7 > chi2 > chi5 > chi1 > beta3 > xi10 > beta2 > chi6 > chi4 > xi7 > beta1 > xi5 > delta1 > xi3 > zeta1 > xi15 > chi8 > xi2 > xi4 > chi10 > xi16 > chi3 > xi6 > chi9 > xi13 > xi1 > zeta3 > xi9 > xi12 > xi14 > zeta2 > xi8` → new `chi1 > beta2 > xi8 > beta1 > xi5 > xi15 > xi1 > chi5 > xi2 > xi10 > chi2 > chi6 > xi6 > chi10 > chi3 > chi7 > delta1 > zeta2 > zeta1 > xi3 > beta3 > chi8 > zeta3 > xi14 > xi13 > xi16 > xi7 > chi9 > xi9 > chi4 > xi4 > xi12` — **dominant coupling changed**

### parity_odd_minimal (29515407, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 8.441; N_eff ≈ 944 (per-param noise floor ≈ 0.02 nats); sum of marginals 10.409 → 2.880 (raw superadditivity violated → restored)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.997 | 0.587 | -2.410 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 3.095 | 0.657 | -2.438 | STRONG → moderate **FLIP** |
| beta3 | arctan_uniform | 2.530 | 0.075 | -2.455 | STRONG → weak **FLIP** |
| d21 | arctan_uniform | 1.786 | 1.560 | -0.226 | STRONG → STRONG |

Ranking (most→least constrained): old `beta2 > beta1 > beta3 > d21` → new `d21 > beta2 > beta1 > beta3` — **dominant coupling changed**

**sup** — joint D_KL 13.267; N_eff ≈ 642 (per-param noise floor ≈ 0.03 nats); sum of marginals 12.112 → 5.092 (raw superadditivity ok → ok)

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 3.102 | 1.759 | -1.343 | STRONG → STRONG |
| beta2 | arctan_uniform | 3.398 | 1.629 | -1.769 | STRONG → STRONG |
| beta3 | arctan_uniform | 2.571 | 0.528 | -2.042 | STRONG → moderate **FLIP** |
| d21 | arctan_uniform | 3.042 | 1.176 | -1.866 | STRONG → STRONG |

Ranking (most→least constrained): old `beta2 > beta1 > d21 > beta3` → new `beta1 > beta2 > d21 > beta3` — **dominant coupling changed**

### parity_odd_full (29687506, plane_wave)

Priors: **no inference.json** — reconstructed from the campaign sbatch template (all arctan ±89; overrides: {'xi': 'log_uniform:1e-3:1e3'}). Verify against scripts/hpc_templates before quoting.

**amp** — joint D_KL 5.933; N_eff ≈ 6 (per-param noise floor ≈ 3.25 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 56.854 → 32.260 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.871 | 1.623 | -1.248 | STRONG → STRONG |
| beta2 | arctan_uniform | 2.634 | 1.543 | -1.091 | STRONG → STRONG |
| beta3 | arctan_uniform | 2.679 | 1.821 | -0.858 | STRONG → STRONG |
| xi | arctan_uniform | 3.256 | 1.608 | -1.648 | STRONG → STRONG |
| delta1 | arctan_uniform | 3.127 | 1.497 | -1.630 | STRONG → STRONG |
| chi | arctan_uniform | 3.040 | 1.812 | -1.229 | STRONG → STRONG |
| zeta1 | arctan_uniform | 2.422 | 1.720 | -0.702 | STRONG → STRONG |
| zeta2 | arctan_uniform | 2.979 | 1.621 | -1.358 | STRONG → STRONG |
| zeta3 | arctan_uniform | 3.028 | 1.595 | -1.434 | STRONG → STRONG |
| d14 | arctan_uniform | 2.827 | 1.532 | -1.295 | STRONG → STRONG |
| d15 | arctan_uniform | 2.876 | 1.692 | -1.183 | STRONG → STRONG |
| d17 | arctan_uniform | 2.757 | 1.498 | -1.259 | STRONG → STRONG |
| d19 | arctan_uniform | 3.005 | 1.521 | -1.484 | STRONG → STRONG |
| d20 | arctan_uniform | 2.766 | 1.531 | -1.235 | STRONG → STRONG |
| d21 | arctan_uniform | 2.717 | 1.786 | -0.931 | STRONG → STRONG |
| tildezeta1 | arctan_uniform | 2.848 | 1.515 | -1.333 | STRONG → STRONG |
| tildezeta2 | arctan_uniform | 2.670 | 1.503 | -1.168 | STRONG → STRONG |
| tildezeta3 | arctan_uniform | 2.582 | 1.541 | -1.041 | STRONG → STRONG |
| tildezeta5 | arctan_uniform | 2.704 | 1.554 | -1.150 | STRONG → STRONG |
| tildezeta6 | arctan_uniform | 3.068 | 1.749 | -1.319 | STRONG → STRONG |

Ranking (most→least constrained): old `xi > delta1 > tildezeta6 > chi > zeta3 > d19 > zeta2 > d15 > beta1 > tildezeta1 > d14 > d20 > d17 > d21 > tildezeta5 > beta3 > tildezeta2 > beta2 > tildezeta3 > zeta1` → new `beta3 > chi > d21 > tildezeta6 > zeta1 > d15 > beta1 > zeta2 > xi > zeta3 > tildezeta5 > beta2 > tildezeta3 > d14 > d20 > d19 > tildezeta1 > tildezeta2 > d17 > delta1` — **dominant coupling changed**

**sup** — joint D_KL 5.103; N_eff ≈ 14 (per-param noise floor ≈ 1.42 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 54.751 → 20.267 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| beta1 | arctan_uniform | 2.863 | 0.889 | -1.974 | STRONG → moderate **FLIP** |
| beta2 | arctan_uniform | 2.025 | 1.217 | -0.808 | STRONG → STRONG |
| beta3 | arctan_uniform | 2.608 | 0.954 | -1.653 | STRONG → moderate **FLIP** |
| xi | arctan_uniform | 2.919 | 1.029 | -1.890 | STRONG → STRONG |
| delta1 | arctan_uniform | 2.680 | 1.046 | -1.633 | STRONG → STRONG |
| chi | arctan_uniform | 2.414 | 1.076 | -1.337 | STRONG → STRONG |
| zeta1 | arctan_uniform | 2.916 | 0.987 | -1.929 | STRONG → moderate **FLIP** |
| zeta2 | arctan_uniform | 2.509 | 0.995 | -1.514 | STRONG → moderate **FLIP** |
| zeta3 | arctan_uniform | 2.211 | 1.130 | -1.081 | STRONG → STRONG |
| d14 | arctan_uniform | 3.014 | 1.162 | -1.852 | STRONG → STRONG |
| d15 | arctan_uniform | 2.990 | 1.046 | -1.944 | STRONG → STRONG |
| d17 | arctan_uniform | 2.819 | 0.925 | -1.893 | STRONG → moderate **FLIP** |
| d19 | arctan_uniform | 2.908 | 1.017 | -1.891 | STRONG → STRONG |
| d20 | arctan_uniform | 2.688 | 0.988 | -1.700 | STRONG → moderate **FLIP** |
| d21 | arctan_uniform | 3.091 | 1.017 | -2.074 | STRONG → STRONG |
| tildezeta1 | arctan_uniform | 2.812 | 0.937 | -1.874 | STRONG → moderate **FLIP** |
| tildezeta2 | arctan_uniform | 2.956 | 0.948 | -2.008 | STRONG → moderate **FLIP** |
| tildezeta3 | arctan_uniform | 2.616 | 0.989 | -1.627 | STRONG → moderate **FLIP** |
| tildezeta5 | arctan_uniform | 2.924 | 0.945 | -1.978 | STRONG → moderate **FLIP** |
| tildezeta6 | arctan_uniform | 2.792 | 0.970 | -1.822 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `d21 > d14 > d15 > tildezeta2 > tildezeta5 > xi > zeta1 > d19 > beta1 > d17 > tildezeta1 > tildezeta6 > d20 > delta1 > tildezeta3 > beta3 > zeta2 > chi > zeta3 > beta2` → new `beta2 > d14 > zeta3 > chi > delta1 > d15 > xi > d19 > d21 > zeta2 > tildezeta3 > d20 > zeta1 > tildezeta6 > beta3 > tildezeta2 > tildezeta5 > tildezeta1 > d17 > beta1` — **dominant coupling changed**

### complete_parity_odd

Skipped — no recomputed file (chain never converged).

### eh_gertsenshtein (29700083, plane_wave)

Priors: recorded in inference.json.

**amp** — joint D_KL 0.010; N_eff ≈ 104 (per-param noise floor ≈ 0.19 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 5.055 → 0.344 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| rho | arctan_uniform | 2.643 | 0.153 | -2.490 | STRONG → moderate **FLIP** |
| sigma | arctan_uniform | 2.411 | 0.191 | -2.221 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `rho > sigma` → new `sigma > rho` — **dominant coupling changed**

**sup** — joint D_KL 0.010; N_eff ≈ 103 (per-param noise floor ≈ 0.19 nats — **floor-dominated, marginals below the floor are noise**); sum of marginals 5.332 → 0.574 (raw superadditivity raw-violated (within N_eff bias allowance))

| param | prior | old | new | Δ | verdict old → new |
|---|---|---|---|---|---|
| rho | arctan_uniform | 2.696 | 0.236 | -2.460 | STRONG → moderate **FLIP** |
| sigma | arctan_uniform | 2.635 | 0.338 | -2.298 | STRONG → moderate **FLIP** |

Ranking (most→least constrained): old `rho > sigma` → new `sigma > rho` — **dominant coupling changed**


## Where pre-fix numbers are recorded (inventory only — nothing edited)

- `manuscript/sections/results.tex` — the per-coupling D tables
  `KLTableEMRCParityPair` (:254), `KLTableYMUnified` (:441),
  `KLTableChiClosureComparison` (~:502), `KLTableDarkPhoton`, cited from :129
  and used for physics statements at :205, :280-284, :337, :460-469, :539-556.
  Verified example: the amp column of `KLTableEMRCParityPair` (3.18/2.94/1.43)
  is `hpc_results/29189748/d1_amp_v3` pre-fix output, which this recompute
  changes to 0.52/0.15/1.86 with the dominant coupling moving from β₁ to δ₁.
- `manuscript/sections/listings/d_kl_*.tex` — generated by
  `scripts/analysis/extract_d_kl.py` from the (now refreshed)
  `parameter_importance.json` files; NOT regenerated by this branch.
- `CAMPAIGN.md` — marginal D_KL quoted in verdict lines at :82, :108, :147-148,
  :183, :196-202, :211, :217-219, :229-233 and the run table :294-344.
- `docs/campaign_plan.md:140-149` — the formal stage thresholds; several stage
  verdicts key on `max D_KL(params)`, which the fix changes in both directions
  (saturated arctan values drop; some previously-fallback values rise).
- Figure scripts that *select* couplings by ranking marginal D_KL
  (`scripts/figures/kl_carrier_corner.py`, `overlay_corner_pair.py`,
  `manuscript/poster/scripts/result_landscape.py`) will pick different
  couplings on regeneration.
- The wider `hpc_results/` set (115 runs with affected prior kinds) keeps its
  pre-fix `importance.json`; recompute on demand when any of them is quoted.
  **Staleness tell**: a pre-v0.48.8 `importance.json` has no `consistency`
  block — if that key is absent, the marginals were computed with the broken
  estimator and must not be quoted. (`tidal analyze <dir> --inference
  --importance` recomputes live with the fixed estimator; the joint D_KL,
  log Z and d_G in the same file were always correct.) Audited consumers:
  the corner-plot headline quotes only the joint `d_kl` (safe), the
  `plot_importance` bar chart always receives a freshly-computed result
  (safe), and `tidal/inference/_atlas.py` does not read marginals at all.

## Caveats

- `complete_parity_odd` was skipped (chain never converged — pre-existing).
- `eh_gertsenshtein` has near-zero information in both directions
  (joint ≈ 0.01 nats); its marginals are noise-level before and after.
- For chains without `inference.json` the priors are reconstructed from sbatch
  templates; each such section is labeled. If any of those runs was launched
  with modified priors (not via the committed template), the reconstruction is
  wrong for the modified parameters — worth a one-time check against the HPC
  submit logs before quoting.
- The raw sum>joint column ignores the N_eff bias allowance on purpose (it is
  the reader-facing red flag); the JSON `consistency.superadditivity_ok` field
  applies the allowance and is the machine-facing verdict.

## Suggested next steps (user decisions)

1. Review the dominant-coupling flips against the specific claims in
   results.tex §4 — especially every "X carries the information" sentence, and
   treat the floor-dominated rescue chains (T7/T9/NP-ceven, low N_eff) as
   "insufficient effective samples for per-coupling claims" rather than as
   ranked results.
2. If accepted: regenerate `manuscript/sections/listings/d_kl_*.tex` via
   `scripts/analysis/extract_d_kl.py`, re-render the ranked-selection figures,
   and revisit CAMPAIGN.md stage verdicts against `docs/campaign_plan.md`
   thresholds.
3. Decide whether the floor-dominated rescue chains warrant clean re-runs on
   HPC (their logZ conclusions are recorded as contaminated in
   `docs/V3_PHASE_TRACKER.md` independently of this bug).
