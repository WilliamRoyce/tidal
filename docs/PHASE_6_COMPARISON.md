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

| Chain | Original Job | Original log Z | Original D_KL (joint / max marginal) | Original MAP | Canonical-Probe MAP-survival | Planned baseline-rerun job | Planned hi-res-rerun job | New log Z | New D_KL | Δ | Verdict |
|-------|--------------|----------------|--------------------------------------|--------------|-------------------------------|----------------------------|--------------------------|-----------|----------|---|---------|
| **D1 amp** (T4 Ricci-EM) | 28520217 | −2.261 ± 0.066 | 1.789 / δ₁=1.679 | α₁=−0.703, α₂=−0.795, α₃=+0.130, δ₁=+0.074 | **REJECTED** (γ_eff = ∞) | TBD | TBD | — | — | — | — |
| **D1 sup** (T4 Ricci-EM) | 28519675 | +15.911 ± 0.152 | 8.911 / δ₁=0.870 | α₁=+0.333, α₂=+0.916, α₃=+0.916, δ₁=−1.291 | SURVIVES (γ_eff = 0.101) | TBD | TBD | — | — | — | — |
| **Stage A v5 sup** (dark-photon plasma CDT) | 28477675 | +0.654 ± 0.056 | 1.978 / ξ=0.414 | mA₂=+0.969, δm=−0.446, ξ=+0.796, α₃=+0.001 | SURVIVES (γ_eff = 0.101) | TBD | TBD | — | — | — | — |

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

## Pre-launch sanity checklist (compute not yet authorised)

- [x] Original sbatch + cmd lines retrieved (3/3)
- [x] Spec drift checks performed (no drift)
- [x] Canonical-probe MAP-survival run on all three MAPs
- [x] Finding flagged: **D1 amp published MAP REJECTED** by canonical probe
- [ ] Decide rerun grid/bounds for each chain (canonical 256/0:100 vs original 64/0:50 vs both)
- [ ] Plan baseline-rerun submissions (same priors + canonical probe gate)
- [ ] Plan hi-res-rerun submissions (canonical probe + grid 256, t_end 10, snapshots 2)
- [ ] Estimate budget — three chains × {baseline, hi-res} = 6 submissions

**Note on the modal A-template cache WIP:** Reverted on 2026-05-03 (commits
`62352c4` and `a068453`); revert commits `b4e70c6` and `f994329`. Full
post-mortem in
[MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md](MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md).
The probe runs above were performed with that WIP stashed at `e9ed623`, so they
are unaffected by the revert. The `grid.bounds.ravel()` reconciliation flagged
earlier is no longer outstanding (the offending code is gone). Bottom line for
this audit: probe results in this ledger remain valid; nothing to redo.
