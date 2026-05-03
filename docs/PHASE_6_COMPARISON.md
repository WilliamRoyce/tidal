# Phase 6 — Comparison Ledger: Old Inference Chains vs Canonical-Probe Reruns

**Created:** 2026-05-03
**Branch:** `hpc/pgt-survey`, HEAD = `e9ed623` (Stage D2 closed)
**Purpose:** Audit old inference chains (28520217, 28519675, 28477675) before
re-running them under the canonical D2.x probe + canonical numerical settings.

---

## Spec drift checks

| Spec file | Last commit (any) | Mtime | Drift since original run? |
|-----------|-------------------|-------|---------------------------|
| `examples/data/torsion_gertsenshtein_nonminimal.json` | `d249031` (chore: un-ignore theory JSONs in examples/data/) | 2026-04-11 12:51 | **No drift.** Only one commit ever touched the path; file content stable since the original D1 derivation. |
| `examples/data/dark_photon_plasma.json` | `9bb67e9` (sign-flip 2026-04-24, pre-dates 28477675 by 2 days) | 2026-04-24 20:47 | **No drift.** Last edit predates 28477675 (submitted 2026-04-26). No edits since 2026-04-26. |
| `examples/dark_photon_plasma/theory.toml` | `9bb67e9` (same sign-flip) | 2026-04-24 20:46 | **No drift.** Same audit. |

**Verdict:** Both specs are bit-identical to what the original three target jobs
used. Comparison ledger is not contaminated by spec drift.

---

## Original sbatch + cmd records (verbatim from `scripts/.hpc_jobs`)

### 28477675 — `dp_suppress_v5` (Stage A v5 sup)
- Submitted: `2026-04-26T19:13:56+00:00`
- Template: `scripts/hpc_templates/polychord_intr.sbatch`
- Spec: `examples/data/dark_photon_plasma.json`
- Priors: `mA2=log_uniform:0.001:1.0`, `deltam=uniform:-0.5:0.5`, `xi=log_uniform:0.05:20.0`, `alpha3=log_uniform:0.001:0.5`
- Constraints: none
- Likelihood: `P_max:minimize`, baseline `sin(kappa*B0*t_end/2)**2`
- Numerical: `--grid-shape 64 --bounds 0:50 --periodic`, `--t-end 10.0`, `--snapshots 2`
- IC: plane-wave on `h_5`, `--ic-wavevector 2.0`, `--ic-amplitude 1e-2`
- Conversion: `--source h_5 --target a_1`
- Sampler: `polychord`, `--nlive 400 --num-repeats 5 --precision-criterion 0.01`

### 28519675 — `t4_suppress_v5` (D1 sup)
- Submitted: `2026-04-27T19:42:04+00:00`
- Template: `scripts/hpc_templates/polychord_standard.sbatch`
- Spec: `examples/data/torsion_gertsenshtein_nonminimal.json`
- Priors: `alpha1=uniform:-1:1`, `alpha2=uniform:-2:2`, `alpha3=log_uniform:0.05:2.0`, `delta1=uniform:-2:2`
- Constraints: none
- Likelihood: `P_max:minimize`
- Numerical, IC, sampler: identical to 28477675 (grid 64, bounds 0:50, k_IC=2.0, nlive=400, num_repeats=5, prec=0.01, snapshots=2)

### 28520217 — `t4_amplify_v5b` (D1 amp)
- Submitted: `2026-04-27T19:59:50+00:00`
- Template: `scripts/hpc_templates/polychord_intr.sbatch`
- Spec: `examples/data/torsion_gertsenshtein_nonminimal.json`
- Priors, numerical, IC, sampler: identical to 28519675
- Likelihood: `P_max:maximize` (the only difference vs 28519675)

**Cross-cutting note:** All three chains used `--grid-shape 64 --bounds 0:50 --ic-wavevector 2.0`.
The canonical D2.x probe uses `--grid-shape 256 --bounds 0:100 --ic-wavevector 0.0628 (≈2π/100)`.
The probe itself (`check_conversion_stability`) scans **all** k modes with a
unit IC at the source slot, so the grid difference primarily increases k-bin
resolution (canonical sees ~129 modes; original saw ~33). The probe call below
uses the canonical (256, 0:100) grid.

---

## Canonical-probe MAP-survival check

**Probe settings:** `t_test=20.0`, `threshold=0.3` (γ_eff > 0.3 ⇒ tachyonic),
`source=h_5`, `target=a_1`, `grid=256`, `bounds=0:100`, `periodic=True`.
Implementation: `tidal.measurement._stability.check_conversion_stability` at HEAD `e9ed623`.

| Chain | Original MAP | γ_eff (canonical) | Verdict |
|-------|--------------|-------------------|---------|
| **28520217 D1 amp** | α₁=−0.703, α₂=−0.795, α₃=+0.130, δ₁=+0.074 | **∞ (overflow)** | **REJECTED** |
| **28519675 D1 sup** | α₁=+0.333, α₂=+0.916, α₃=+0.916, δ₁=−1.291 | 0.101 | SURVIVES |
| **28477675 Stage A v5 sup** | mA₂=+0.969, δm=−0.446, ξ=+0.796, α₃=+0.001 | 0.101 | SURVIVES |

### Finding to flag

**The published D1 amp MAP (28520217) is in the tachyonic-overflow region of
the canonical probe.** A canonical-probe-rejecting chain rerun cannot land on
this MAP — the maximum-likelihood point under the constrained likelihood will
shift to whichever stable region of (α₁, α₂, α₃, δ₁) attains the highest P_max
without triggering the probe. The current published amplification ceiling
(A_max = 1.26 at this MAP, log Z = −2.26) is therefore *contingent on the
absence of probe filtering at chain time*. The new chain will have:

- **Smaller effective prior volume** (the tachyonic strip carved out).
- **Different MAP location** (likely on the outer boundary of the surviving region).
- **Different log Z** (the integral shrinks; could shift in either direction depending on whether the rejected region had high or low likelihood).

The D1 sup and Stage A v5 sup MAPs both lie comfortably inside the stable
region (γ_eff ≈ 0.10, three times below the threshold of 0.3) — those chains
should re-converge to nearly the same MAPs and log Z.

---

## Comparison ledger

| Chain | Original Job | Original log Z | Original D_KL (joint / max marginal) | Original MAP | Canonical-Probe MAP-survival | Replay survival (wt / unwt) | Δlog Z (predicted) | Predicted new log Z | Planned baseline-rerun job | Planned hi-res-rerun job | New log Z | New D_KL | Verdict (predicted) |
|-------|--------------|----------------|--------------------------------------|--------------|-------------------------------|------------------------------|--------------------|---------------------|----------------------------|--------------------------|-----------|----------|---------------------|
| **D1 amp** (T4 Ricci-EM) | 28520217 | −2.261 ± 0.066 | 1.789 / δ₁=1.679 | α₁=−0.703, α₂=−0.795, α₃=+0.130, δ₁=+0.074 | **REJECTED** (γ_eff = ∞) | **0.965 / 0.703** | **−0.036** | **−2.297** | **28789437** (icelake-standard 6h) | TBD | — | — | **MAP shift, Δlog Z small (~−0.04). Amplification claim survives but ceiling will move.** |
| **D1 sup** (T4 Ricci-EM) | 28519675 | +15.911 ± 0.152 | 8.911 / δ₁=0.870 | α₁=+0.333, α₂=+0.916, α₃=+0.916, δ₁=−1.291 | SURVIVES (γ_eff = 0.101) | **1.000 / 0.995** | **−0.000** | **+15.911** | **28789439** (icelake-standard 6h) | TBD | — | — | Reproduce original — no change expected. |
| **Stage A v5 sup** (dark-photon plasma CDT) | 28477675 | +0.654 ± 0.056 | 1.978 / ξ=0.414 | mA₂=+0.969, δm=−0.446, ξ=+0.796, α₃=+0.001 | SURVIVES (γ_eff = 0.101) | **0.858 / 0.942** | **−0.154** | **+0.500** | TBD | TBD | — | — | Reproduce original within 0.15 nats; α₃≳0.1 region shrinks. |

### Original marginal D_KL (full)

```
28520217 (D1 amp):
  alpha1: 0.0530 nats
  alpha2: 0.0678 nats
  alpha3: 0.0115 nats
  delta1: 1.6794 nats   <- dominates; chain is mostly informative on δ₁

28519675 (D1 sup):
  alpha1: 0.2957 nats
  alpha2: 0.2429 nats
  alpha3: 0.2787 nats
  delta1: 0.8702 nats

28477675 (Stage A v5 sup):
  mA2:    0.3656 nats
  deltam: 0.0336 nats
  xi:     0.4143 nats
  alpha3: 0.1040 nats
```

### Joint D_KL and Gaussianised D_G

```
28520217: D_KL = 1.789 ± 0.055,  D_G = 1.462 ± 0.075
28519675: D_KL = 8.911 ± 0.132,  D_G = 6.858 ± 0.410
28477675: D_KL = 1.978 ± 0.155,  D_G = 20.689 ± 1.040  (highly non-Gaussian posterior)
```

---

## Phase 6.B — Per-sample probe replay on the full chain

Replayed each weighted posterior sample in `tidal.txt` through the canonical
probe (settings exactly as above) to measure the survival rate, not just the
MAP. The hand-wave Δlog Z bound is `log(weighted_survival)`; this is also the
*tight* bound when the original chain explored the surviving region adequately,
because the survivors' likelihoods are unchanged and only the prior mass
shrinks.

Replay code: `/tmp/replay_probe.py` (saves per-sample verdicts + γ_eff to
`/tmp/probe_replay.npz`); analysis: `/tmp/replay_analysis.py`.

### Survival rates

| Chain | N | Unweighted pass | Weighted pass | Δlog Z | Predicted new log Z | Tachyonic-overflow count | γ_eff(pass).max | γ_eff(fail).min |
|-------|---|-----------------|---------------|--------|---------------------|--------------------------|-----------------|-----------------|
| 28520217 D1 amp | 3248 | 2284/3248 = **0.703** | 749.81/777.41 = **0.965** | **−0.036** | **−2.297** | 911 (28%) | 0.2997 | 0.3025 |
| 28519675 D1 sup | 6527 | 6496/6527 = 0.995 | 66.190/66.191 = **1.000** | −0.0001 | +15.911 | 3 (0.05%) | 0.2979 | 0.3005 |
| 28477675 Stage A v5 sup | 6128 | 5770/6128 = 0.942 | 774.88/903.45 = **0.858** | **−0.154** | **+0.500** | 351 (5.7%) | 0.1018 | 0.6269 |

**Critical observation:** for D1 amp, 28% of *raw* samples have γ_eff = ∞ and
another 1.6% have finite γ_eff > 0.3 (964 fails total = 30% of samples), but
the *weighted* fail mass is only 3.5%. PolyChord assigned very low posterior
weight to the unstable region — the chain was already mostly avoiding it
through the likelihood, just not strictly. The MAP, however, sits at a sharp
peak that lies inside the rejected strip (γ_eff = ∞).

The `γ_eff(pass).max` ≈ 0.30 and `γ_eff(fail).min` ≈ 0.30 separation is sharp:
the threshold cleanly bisects pass/fail with no near-miss ambiguity for D1.
For Stage A v5 the gap is wider (0.10 vs 0.63) — fails are decisively
tachyonic, no borderline cases.

### Rejection-region map — D1 amp (28520217)

#### 1D stratification (weighted rejection rate per bin)

```
δ₁ bin           | reject | n     | weight
[-2.0, -1.5)     | 0.223  |  76   |   0.58
[-1.5, -1.0)     | 0.243  |  87   |   0.54
[-1.0, -0.5)     | 0.002  | 117   |   4.32
[-0.5,  0.0)     | 0.035  | 1335  | 383.13   <-- MAP region (δ₁=+0.074)
[ 0.0, +0.5)     | 0.036  | 1322  | 380.17   <-- MAP region
[+0.5, +1.0)     | 0.004  | 150   |   7.59
[+1.0, +1.5)     | 0.154  |  76   |   0.42
[+1.5, +2.0)     | 0.164  |  85   |   0.66

α₂ bin           | reject | n     | weight
[-2.0, -1.5)     | 0.206  | 664   |  57.62   <-- α₂ near MAP (=−0.795)
[-1.5, -1.0)     | 0.131  | 502   |  61.90
[-1.0, -0.5)     | 0.063  | 420   |  83.16   <-- MAP α₂ band
[-0.5,  0.0)     | 0.022  | 324   |  77.86
[ 0.0, +0.5)     | 0.006  | 311   | 110.53
[+0.5, +1.0)     | 0.000  | 324   | 113.11
[+1.0, +1.5)     | 0.000  | 325   | 116.84
[+1.5, +2.0)     | 0.000  | 378   | 156.40

α₁ bin           | reject | n     | weight
[-1.0, -0.5)     | 0.155  | 1187  | 129.81   <-- MAP α₁ band (=−0.703)
[-0.5,  0.0)     | 0.045  |  760  | 162.01
[ 0.0, +0.5)     | 0.001  |  570  | 214.26
[+0.5, +1.0)     | 0.000  |  731  | 271.33
```

#### 2D rejection map (δ₁ × α₂, weighted)

```
              α₂∈[-2,-1)    α₂∈[-1,+0)    α₂∈[+0,+1)    α₂∈[+1,+2)
δ₁∈[-2,-1)    reject=0.31   reject=0.29   reject=0.29   reject=0.00
δ₁∈[-1,+0)    reject=0.15   reject=0.04   reject=0.00   reject=0.00   <-- MAP cell
δ₁∈[+0,+1)    reject=0.19   reject=0.04   reject=0.00   reject=0.00
δ₁∈[+1,+2)    reject=0.22   reject=0.20   reject=0.14   reject=0.00
```

The MAP at (α₂=−0.795, δ₁=+0.074) sits in the (α₂∈[−1,0), δ₁∈[0,+1)) cell —
average rejection 4%, but the MAP itself is exactly on the *boundary* of the
unstable region.  The rerun MAP will move toward the α₂≳0 region (rejection
≤1%) where the same likelihood structure is available with full stability
margin.

### Rejection-region map — D1 sup (28519675)

```
δ₁ bin           | reject | n
[-2.0, -1.5)     | 0.000  | 442
[-1.5, -1.0)     | 0.000  | 2227
[-1.0, -0.5)     | 0.000  | 1234
[-0.5,  0.0)     | 0.000  | 208
[ 0.0, +0.5)     | 0.006  |  81
[+0.5, +1.0)     | 0.000  | 1027
[+1.0, +1.5)     | 0.000  | 961
[+1.5, +2.0)     | 0.000  | 347

α₂ bin           | reject | n
all bins         | 0.000  | every sample
```

Suppression chain is stable everywhere it explored — confirms expectation.
The 31 rejections are isolated artefacts (γ_eff just above 0.3) with negligible
posterior weight.

### Rejection-region map — Stage A v5 sup (28477675)

```
δm bin           | reject | n
[-0.50, -0.30)   | 0.161  | 2072
[-0.30, -0.10)   | 0.144  |  832
[-0.10, +0.10)   | 0.165  |  324
[+0.10, +0.30)   | 0.135  |  877
[+0.30, +0.50)   | 0.116  | 2023

α₃ bin           | reject | n
[ 0.000, 0.100)  | 0.115  | 5625   <-- MAP α₃ band (=+0.001)
[ 0.100, 0.200)  | 0.277  |  230
[ 0.200, 0.300)  | 0.374  |  141
[ 0.300, 0.400)  | 0.324  |   70
[ 0.400, 0.500)  | 0.145  |   62
```

α₃ rejection scales with α₃ magnitude — the Proca mass term destabilises the
torsion mode at finite α₃ (the CDT trace channel `m² = +2·α₃`). Rejection rate
peaks at α₃≈0.25 (37%), where the conversion-stability gap closes. The MAP
α₃=0.001 is in the safest bin (12% rejection); rerun should preserve the MAP
position.

### Quantitative bound on Δlog Z

For each chain the formal bound is

    Δlog Z = log(Σ w_i · 𝟙[stable_i] / Σ w_i)

This is a *tight* bound under the assumption that the original chain explored
the surviving region adequately (no missing modes). The chain explored the
unstable region too — but assigned it low weight because the likelihood was
either bounded by `P_max ≤ 1` or driven by destructive interference. For all
three chains the unweighted-vs-weighted survival ratio confirms the chain
preferentially weighted the stable region.

Expected accuracy of these predictions:
- **D1 sup**: bound is essentially exact (Δ < 0.01 nats); rerun is a true reproducibility check.
- **D1 amp**: bound assumes surviving samples are sufficient to converge; with only 70% unweighted survival there is some risk PolyChord cannot re-converge from this support alone, but the weighted pass-rate is high enough (96%) that the original chain's stable-region statistics should suffice.
- **Stage A v5 sup**: bound's reliability is intermediate; the α₃-dependent rejection introduces mild non-uniformity in the surviving region.

### Implications for Phase 6.C launches

- **No soft constraints needed.** With weighted survival ≥ 86% on all three
  chains, the canonical probe acts as a hard prior cut without significantly
  shrinking the explorable region. PolyChord can run with the probe gating
  the likelihood directly (return −∞ on probe-rejection); no soft β-style
  constraints required.
- **Budget unchanged.** With Δlog Z bounded by 0.16 nats, log Z reproducibility
  is well within typical PolyChord error bars at nlive=400 (±0.07–0.15) and
  nlive=800 (±0.05). The reruns will be of similar wall-time to the originals.
- **D1 amp MAP location is the science finding.** The expected MAP shift —
  toward α₂≳0, α₁≳+0.5, |δ₁| small — needs to be *measured*, not predicted from
  the old chain's posterior alone. The old chain's posterior near the new MAP
  is sparse (the chain spent its weight near the old MAP), so the old chain
  cannot predict the new MAP location with confidence even though it can
  predict log Z to ~0.04 nats.

---

## Pre-launch sanity checklist (compute not yet authorised)

- [x] Original sbatch + cmd lines retrieved (3/3)
- [x] Spec drift checks performed (no drift)
- [x] Canonical-probe MAP-survival run on all three MAPs
- [x] Finding flagged: **D1 amp published MAP REJECTED** by canonical probe
- [x] Per-sample replay through canonical probe (3/3 chains, 15903 samples total)
- [x] Survival rates + Δlog Z estimates committed
- [x] Rejection-region maps for D1 amp, D1 sup, Stage A v5 sup
- [x] Predicted verdicts populated
- [x] Decide rerun grid/bounds for each chain — chose canonical 256/0:100 (probe and chain matched)
- [x] Plan baseline-rerun submissions (same priors + canonical probe gate, **no** soft constraints)
- [ ] Plan hi-res-rerun submissions (canonical probe + grid 256, t_end 10, snapshots 2)
- [ ] Estimate budget — three chains × {baseline, hi-res} = 6 submissions

---

## Phase 6.C — Baseline rerun submissions (2026-05-03)

**Decision (recorded):** No soft β-style constraints. Both chains use the same
priors as the originals so the log Z comparison isolates the architectural-
correction effect (canonical probe + Hwang–Noh gate + canonical grid) from
any prior-narrowing effect.

**Submission tally:**

| Chain | New job | Template | Partition / QOS | Wall budget | nlive | Probe | Constraints | State |
|-------|---------|----------|-----------------|-------------|-------|-------|-------------|-------|
| D1 amp clean (std) | **28789437** | `polychord_standard.sbatch` | icelake / cpu2 | 06:00:00 | 400 | canonical (t_test=20, threshold=0.3, all-k unit-IC) + inline P_max>0.5 Hwang–Noh gate | none | PENDING |
| D1 sup clean (std) | **28789439** | `polychord_standard.sbatch` | icelake / cpu2 | 06:00:00 | 400 | same | none | PENDING |
| D1 amp clean (INTR cross-check) | **28789579** | `polychord_intr.sbatch` | icelake / intr | 01:00:00 | 400 | same | none | **COMPLETED** (21 min) |

**Wall-time anticipation note:** the original 28520217 took 3:30 wall; with the
canonical probe rejecting ~30% of raw samples (per Phase 6.B replay), the new
chain may need 1.4–1.6× more likelihood evaluations to converge to the same
nlive. INTR's 1h limit was therefore deemed too tight for a single round;
submitted to icelake standard 6h as primary. Sup chain has ≈100% probe survival
— expected wall close to original 16:30.

**INTR cross-check (28789579):** submitted in parallel with the standard D1
amp job, using `--campaign d1_amp_intr` so PolyChord's `tidal.resume` is
written to the stable path
`hpc_results/campaigns/d1_amp_intr/d1_amp_intr/_chains/`. On INTR's 1 h
timeout, resubmit with `--read-resume` and the same campaign dir to continue
from the checkpoint. Whichever chain (standard or INTR-resumed) reaches
PolyChord's stop criterion first becomes the primary; the other is the
cross-check.

**Pre-flight check (~30 s post-submit):** both jobs PENDING with reason
`None` (clean queue), no fast-fail. SLURM logs not yet written (jobs not
started).

**Note on the modal A-template cache WIP:** Reverted on 2026-05-03 (commits
`62352c4` and `a068453`); revert commits `b4e70c6` and `f994329`. Full
post-mortem in
[MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md](MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md).
The probe runs above were performed with that WIP stashed at `e9ed623`, so they
are unaffected by the revert. The `grid.bounds.ravel()` reconciliation flagged
earlier is no longer outstanding (the offending code is gone). Bottom line for
this audit: probe results in this ledger remain valid; nothing to redo.

---

## Phase 6.C — INTR result: 28789579 D1 amp clean (2026-05-03)

**Completed in 21 minutes** (nlive=0; fully converged, not timed out).

| Quantity | Value |
|----------|-------|
| log Z (new INTR) | **+2.135 ± 0.059** |
| log Z (original 28520217) | −2.261 ± 0.066 |
| Δlog Z (new − old) | **+4.396 nats** |
| New MAP | α₁=+0.838, α₂=−1.830, α₃=+0.841, δ₁=+1.943 |
| Old MAP (28520217) | α₁=−0.703, α₂=−0.795, α₃=+0.130, δ₁=+0.074 |
| New MAP P_max | 0.105 (A ≈ 42×) |
| Old MAP A | 1.25× (barely above GR) |
| logL range (new) | [−10.3, +5.3] → A ∈ [3×10⁻⁵, 200] |
| logL range (old) | [−19.1, +0.23] → A ∈ [5×10⁻⁹, 1.25] |
| Prior stability sweep | 957/5000 rejected (19.1% tachyonic) |
| n_posterior / n_like | 4444 / 246371 |

**P_max perturbativity breakdown of the new posterior:**

| P_max range | Weighted fraction |
|-------------|-------------------|
| [0, 0.01) | 9.4% |
| [0.01, 0.05) | 14.3% |
| [0.05, 0.10) | 12.4% |
| [0.10, 0.20) | 20.6% |
| [0.20, 0.30) | 17.4% |
| [0.30, 0.50) | 26.0% |
| **P_max > 0.1 (borderline)** | **64%** |
| **P_max > 0.3 (near-saturated)** | **26%** |

**Physical interpretation of the log Z shift:**

The +4.4 nat shift is NOT purely the probe correction predicted in Phase 6.B (−0.036 nats).
Two compounding architectural changes explain the radical difference:

1. **k_IC change**: Original D1 runs used k_IC = 2.0 (grid 64, L=50). The canonical
   reruns use k_IC = 2π/100 ≈ 0.0628 (grid 256, L=100). This is a factor of ~32×
   change in the probed wavenumber. In the torsion-nonminimal model, the Ricci-EM
   identity coupling is k-independent but torsion gradient couplings scale with k;
   the mode structure of the conversion landscape differs substantially.

2. **Stability probe effect**: Old run had no probe gate → explored tachyonic parameter
   regions where the graviton instability caused the source field to grow faster than
   the target, leading to *suppression* (A < 1). With the probe, tachyonic regions
   return −∞ and the chain samples only the stable sector where genuine photon
   amplification accumulates.

**Combined effect**: The old run was suppression-dominated (tachyonic graviton growth
outpaced photon build-up). The new canonical run reveals a large stable-sector region
with genuine amplification up to A ≈ 200 (P_max ≈ 0.5, just below the Hwang-Noh gate).

**Perturbativity concern:** 64% of the new posterior has P_max > 0.1, with 26% above 0.3.
The linearized conversion measurement is valid for P_max ≪ 1; at P_max = 0.4 the
perturbative approximation is questionable. The log Z = +2.135 is real (within the
accepted Hwang-Noh criterion of P_max < 0.5) but a cross-check at smaller B₀ is needed
to confirm the A values hold under perturbativity. The standard icelake jobs
(28789437/28789439) will provide a second draw; if they agree, they corroborate the
canonical result.

**Corner plot:** `hpc_results/28789579/corner_28789579_amp_clean.png`
