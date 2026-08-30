# Results Amendments

**Purpose.** The thesis manuscript (`manuscript/`) is a frozen archive: it records
what was claimed, on what evidence, at submission time, and is never edited. This
document is the living record of amendments — what has since been superseded, what
stands, and why — so the current state of results is always readable without
rewriting history. Each amendment cites its evidence and the GitHub trail.
`CAMPAIGN.md` mirrors each amendment as an append-only dated entry.

Conventions: `stands` = the archived claim is supported by the corrected analysis;
`superseded` = the archived number/claim is replaced by the corrected value here;
`floor` = no claim is supportable either way at the chain's effective sample size;
`strengthened` = the corrected analysis supports the archived claim more strongly.

---

## Amendment 1 (2026-08-18) — Marginal D_KL estimator correction (GH #420)

### What happened

`compute_parameter_importance` transformed only `log_uniform` priors into the space
where the prior is uniform before histogramming. `arctan_uniform` columns (most
campaign couplings) were scored in linear space against a uniform reference over
recorded bounds the sampler never used, inflating and **saturating** marginals near
the estimator ceiling log(40) ≈ 3.69 nats; `radial_angular` columns fell to a
self-referential fallback. Fixed in v0.48.8 (PR #426); all 13 viable publication
chains recomputed (PR #431); full evidence in `docs/dkl_recompute_report.md`.
Two further defects were found and fixed during remediation: fabricated all-arctan
priors for pre-schema chains (#434 — T7/T6 actually used log-uniform ξ per the
sbatch templates) and a column-rename mismatch in the recompute script. The
estimator now self-checks (superadditivity, saturation, noise floor, N_eff —
GH #433).

### What is affected, and what is not

The archived per-coupling tables report three quantities per coupling. They are
NOT equally affected:

| quantity | symbol in the archive | status |
|---|---|---|
| amplification-only marginal D_KL | D+ | **superseded** — recomputed below |
| suppression-only marginal D_KL | D− | **superseded** — recomputed below |
| amplification-vs-suppression cross divergence | D₊‖₋ | **unchanged** (prior-free estimator, #420 never touched it) |

Also unchanged and quotable directly from the archive: joint D_KL, log Z, Bayes
factors, d_G (all anesthetic-computed), and every `uniform`/`log_uniform` marginal
(bit-identical by construction; verified e.g. β₃ 0.0049 → 0.0049, ξ 1.20 → 1.20).

**Load-bearing consequence:** most of the thesis's operator-identification prose
keys on D₊‖₋, not on D+/D−, so those claims stand. Better: the corrected D+/D−
marginals now largely *agree* with the D₊‖₋ rankings (e.g. EM-RC amp: corrected
top coupling is δ₁, matching its D₊‖₋ = 20.8; Bahamonde: corrected top is ξ,
matching D₊‖₋ = 5.55) — a coherence the broken estimator obscured. What is
genuinely superseded is every quoted D+/D− number, plus the small set of claims
that rested on marginal-to-marginal comparisons (see the claim table).

### Claim-by-claim status (results.tex)

| archive site | claim | basis | status |
|---|---|---|---|
| :203-207 | δ₁ an order of magnitude above the other couplings (D₊‖₋ = 20.8) | cross | stands |
| :238-242 (`KLTableEMRCParityPair`) | D+/D− cells | marginal | superseded — corrected table below |
| :278-284 | parity-even concentrates on δ₁ (20.8); parity-odd vertex drops to 8.8 and β₁/β₂ become comparable (5.44/4.35) | cross | stands |
| :336-342 | each YM truncation concentrates D on a different operator: ξ (Bahamonde, Barker), β₁ (Shapiro), χ (union) | cross | stands (the claim keys on D₊‖₋, which is unchanged); corrected amp marginals agree for Bahamonde (ξ tops, 1.20) but NOT for Barker (χ tops at 1.02, ξ second at 0.83 under its true log-uniform prior — see the provenance note below) or Shapiro (ζ₂ tops) — marginal and cross rankings are different questions and only the cross one is claimed |
| :422-430 (`KLTableYMUnified`) | D+/D− cells | marginal | superseded — corrected table below |
| :457-465 | χ₁ dominant axis, largest D₊‖₋ = 4.2 | cross | stands (D+/D− cells floor-limited, see below) |
| :536-541 | NP control matches YM union "within the bootstrap uncertainty of the per-coupling marginal D in every operator direction" | **marginal** | **superseded** — corrected marginals differ by up to 1.27 nats (sup χ: 1.43 propagating vs 0.17 NP; amp β₁: 0.05 vs 0.64), far beyond floor at N_eff ≈ 1100–2900. The qualitative control conclusion (amplification without propagating torsion) is separately supported by the corner overlays and log Z and is NOT withdrawn here, but this quantitative support line no longer holds; a future campaign should re-examine it |
| :543-552 | same control claim at the parity-even closure (17D vs 18D) | marginal | floor — both chains have N_eff 80–248 (floors 0.08–0.24 nats); corrected differences ≤ 0.13 nats are within floor, so the comparison has no resolving power either way |
| :554-559 | dominant operator: β₁ in the YM control, χ₁ in the closure control | cross | stands (corrected NP amp marginal also puts β₁ top, 0.64) |
| :572-583 (`KLTableChiClosureComparison`) | D+/D− cells | marginal | superseded and floor-limited — corrected table below |
| :633-643 | T9 kinetic closure "per-coupling D figures not quantitatively meaningful" (overflow) | — | **strengthened & quantified**: N_eff = 21 (amp), noise floor ≈ 0.93 nats/param — no per-coupling value in this sector is signal (GH #433) |
| :756-765 | dark photon: Δm at 3.4 nats, m_A² at 1.5; "marginal D values … not directly cross-comparable" | cross (+ an archived caveat) | stands — the archived caveat anticipated exactly the per-prior-shape problem #420 fixed |
| :786-789 (`KLTableDarkPhoton`) | D+/D− cells | marginal | superseded — corrected table below |
| :836-846 | Euler-Heisenberg null: marginals "statistically indistinguishable from the prior" | marginal | **strengthened** — corrected marginals sit at the chain's noise floor (0.15–0.34 vs floor ≈ 0.19 at N_eff ≈ 104), i.e. no detectable information gain, and the joint D_KL ≈ 0.01 nats; the commented-out KLTableEH values 2.4–2.7 were pure estimator artifact |

### Corrected per-coupling tables

Values in nats; "archived" = as printed in the thesis tables (equal to the pre-fix
`parameter_importance.json`, to rounding); "corrected" = post-#420 recompute
(PR #431). `(≤ floor)` marks values at or below the chain's per-parameter noise
floor `(n_bins−1)/(2·N_eff)` — estimator noise, not constraint (GH #433).

### `KLTableEMRCParityPair`

**Parity-even EM-RC (delta1)** (`hpc_results/29189748/d1_amp_v3`; N_eff amp/sup ≈ 5694/5224)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 3.18 | 0.52 | 2.87 | 0.27 | 0.71 |
| beta2 | 2.94 | 0.15 | 2.77 | 0.26 | 0.54 |
| beta3 | 0.00 | 0.00 | 0.15 | 0.15 | 0.20 |
| delta1 | 1.43 | 1.86 | 3.00 | 0.77 | 20.78 |

Top-3 by D+ (amp marginal): archived `beta1 > beta2 > delta1` → corrected `delta1 > beta1 > beta2`.

**Parity-odd (d21)** (`hpc_results/29515407/t6_minimal_amp_v3`; N_eff amp/sup ≈ 944/642)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 3.00 | 0.59 | 3.10 | 1.76 | 5.44 |
| beta2 | 3.10 | 0.66 | 3.40 | 1.63 | 4.35 |
| beta3 | 2.53 | 0.08 | 2.57 | 0.53 | 0.81 |
| d21 | 1.79 | 1.56 | 3.04 | 1.18 | 8.83 |

Top-3 by D+ (amp marginal): archived `beta2 > beta1 > beta3` → corrected `d21 > beta2 > beta1`.

### `KLTableYMUnified`

**Bahamonde** (`hpc_results/29507332/d20_bahamonde_amp_v3_pub`; N_eff amp/sup ≈ 1063/2404)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.93 | 0.59 | 2.93 | 0.28 | 0.94 |
| beta2 | 2.67 | 0.31 | 2.88 | 0.12 | 0.52 |
| beta3 | 2.82 | 0.18 | 2.50 | 0.03 | 0.24 |
| xi | 1.20 | 1.20 | 1.77 | 1.77 | 5.55 |
| delta1 | 3.02 | 0.75 | 3.06 | 0.62 | 0.08 |

Top-3 by D+ (amp marginal): archived `delta1 > beta1 > beta3` → corrected `xi > delta1 > beta1`.

**Barker** (`hpc_results/29507332/d21_barker_amp_v3_pub`; N_eff amp/sup ≈ 1340/2895)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.36 | 0.42 | 2.83 | 0.16 | 1.30 |
| beta2 | 2.66 | 0.11 | 2.65 | 0.08 | 0.21 |
| beta3 | 2.97 | 0.25 | 2.55 | 0.01 | 0.22 |
| xi | 2.65 | 0.83 | 1.59 | 1.59 | 4.87 |
| delta1 | 3.09 | 0.76 | 2.98 | 0.61 | 0.52 |
| chi | 3.00 | 1.02 | 3.11 | 1.92 | 0.64 |

Top-3 by D+ (amp marginal): archived `delta1 > chi > beta3` → corrected `chi > xi > delta1`.

**Shapiro** (`hpc_results/29468763/d22_shapiro_amp_v3`; N_eff amp/sup ≈ 1341/3302)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.37 | 0.14 | 3.04 | 0.40 | 4.77 |
| beta2 | 2.48 | 0.03 | 2.83 | 0.14 | 1.14 |
| beta3 | 2.62 | 0.10 | 2.74 | 0.11 | 0.32 |
| xi | 0.15 | 0.15 | 0.44 | 0.44 | 0.53 |
| delta1 | 2.87 | 0.12 | 2.94 | 0.21 | 0.20 |
| zeta1 | 2.67 | 0.07 | 1.66 | 0.58 | 0.54 |
| zeta2 | 2.96 | 0.15 | 2.97 | 0.26 | 0.16 |
| zeta3 | 2.63 | 0.07 | 2.80 | 0.09 | 0.25 |

Top-3 by D+ (amp marginal): archived `zeta2 > delta1 > zeta1` → corrected `zeta2 > xi > beta1`.

**Union** (`hpc_results/29468763/d23_full_amp_v3`; N_eff amp/sup ≈ 1297/2870)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.56 | 0.05 | 2.80 | 0.41 | 0.77 |
| beta2 | 2.47 | 0.04 | 2.64 | 0.18 | 0.44 |
| beta3 | 2.63 | 0.05 | 2.61 | 0.03 | 0.13 |
| xi | 0.18 | 0.18 | 0.31 | 0.31 | 0.40 |
| delta1 | 2.93 | 0.11 | 2.71 | 0.06 | 0.17 |
| chi | 2.31 | 0.15 | 3.01 | 1.43 | 3.86 |
| zeta1 | 2.63 | 0.10 | 1.96 | 0.24 | 0.38 |
| zeta2 | 2.81 | 0.11 | 2.86 | 0.09 | 0.25 |
| zeta3 | 2.71 | 0.07 | 2.72 | 0.05 | 0.27 |

Top-3 by D+ (amp marginal): archived `delta1 > zeta2 > zeta3` → corrected `xi > chi > delta1`.

**NP control (xi=0)** (`hpc_results/29700462/np_amp_v1`; N_eff amp/sup ≈ 1146/2867)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.19 | 0.64 | 2.68 | 0.09 | 3.48 |
| beta2 | 2.35 | 0.25 | 2.54 | 0.04 | 0.66 |
| beta3 | 2.47 | 0.16 | 2.75 | 0.08 | 0.43 |
| delta1 | 2.76 | 0.19 | 2.73 | 0.07 | 0.26 |
| chi | 2.81 | 0.35 | 2.89 | 0.17 | 0.57 |
| zeta1 | 2.82 | 0.31 | 2.56 | 0.04 | 0.23 |
| zeta2 | 2.99 | 0.38 | 2.05 | 0.22 | 0.74 |
| zeta3 | 2.76 | 0.26 | 2.26 | 0.11 | 0.43 |

Top-3 by D+ (amp marginal): archived `zeta2 > zeta1 > chi` → corrected `beta1 > zeta2 > chi`.

### `KLTableChiClosureComparison`

**NP control (17D)** (`hpc_results/29705560/np_ceven_amp_v1`; N_eff amp/sup ≈ 80/211)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.60 | 0.22 (≤ floor) | 2.76 | 0.10 | 0.57 |
| beta2 | 2.50 | 0.29 | 2.68 | 0.09 | 0.77 |
| beta3 | 2.47 | 0.29 | 2.54 | 0.08 (≤ floor) | 0.49 |
| delta1 | 2.48 | 0.29 | 2.46 | 0.09 (≤ floor) | 0.85 |
| zeta1 | 2.68 | 0.20 (≤ floor) | 2.51 | 0.10 | 0.31 |
| zeta2 | 2.62 | 0.22 (≤ floor) | 2.51 | 0.09 (≤ floor) | 0.38 |
| zeta3 | 2.30 | 0.26 | 2.47 | 0.08 (≤ floor) | 0.89 |
| chi1 | 1.72 | 0.71 | 2.82 | 0.15 | 5.78 |
| chi2 | 2.39 | 0.29 | 2.78 | 0.10 | 0.89 |
| chi3 | 2.58 | 0.20 (≤ floor) | 2.45 | 0.12 | 0.76 |
| chi4 | 2.50 | 0.20 (≤ floor) | 2.61 | 0.06 (≤ floor) | 0.63 |
| chi5 | 2.68 | 0.19 (≤ floor) | 2.63 | 0.10 | 0.27 |
| chi6 | 2.54 | 0.17 (≤ floor) | 2.61 | 0.08 (≤ floor) | 0.57 |
| chi7 | 2.37 | 0.25 | 2.80 | 0.11 | 1.44 |
| chi8 | 2.42 | 0.26 | 2.43 | 0.09 (≤ floor) | 0.96 |
| chi9 | 2.55 | 0.21 (≤ floor) | 2.70 | 0.09 (≤ floor) | 0.49 |
| chi10 | 2.72 | 0.25 | 2.54 | 0.07 (≤ floor) | 0.29 |

Top-3 by D+ (amp marginal): archived `chi10 > zeta1 > chi5` → corrected `chi1 > delta1 > beta3`.

**Propagating (18D)** (`hpc_results/29682868/t7_amp_v2`; N_eff amp/sup ≈ 118/248)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| beta1 | 2.59 | 0.22 | 2.73 | 0.10 | 0.52 |
| beta2 | 2.67 | 0.18 | 2.73 | 0.09 | 0.46 |
| beta3 | 2.68 | 0.16 (≤ floor) | 2.61 | 0.05 (≤ floor) | 0.40 |
| xi | 2.45 | 0.18 | 2.38 | 0.09 | 0.44 |
| delta1 | 2.55 | 0.20 | 2.53 | 0.07 (≤ floor) | 0.36 |
| zeta1 | 2.72 | 0.20 | 2.60 | 0.08 (≤ floor) | 0.35 |
| zeta2 | 2.51 | 0.12 (≤ floor) | 2.38 | 0.08 | 0.45 |
| zeta3 | 2.37 | 0.19 | 2.43 | 0.08 | 0.42 |
| chi1 | 1.68 | 0.67 | 2.82 | 0.12 | 4.18 |
| chi2 | 2.37 | 0.22 | 2.66 | 0.06 (≤ floor) | 0.85 |
| chi3 | 2.54 | 0.18 | 2.53 | 0.08 (≤ floor) | 0.47 |
| chi4 | 2.40 | 0.16 (≤ floor) | 2.57 | 0.07 (≤ floor) | 0.85 |
| chi5 | 2.74 | 0.15 (≤ floor) | 2.55 | 0.06 (≤ floor) | 0.22 |
| chi6 | 2.50 | 0.18 | 2.68 | 0.08 (≤ floor) | 0.30 |
| chi7 | 2.38 | 0.22 | 2.81 | 0.10 | 0.66 |
| chi8 | 2.39 | 0.22 | 2.56 | 0.09 | 0.41 |
| chi9 | 2.44 | 0.11 (≤ floor) | 2.70 | 0.09 | 0.78 |
| chi10 | 2.56 | 0.22 | 2.43 | 0.09 | 0.37 |

Top-3 by D+ (amp marginal): archived `chi5 > zeta1 > beta3` → corrected `chi1 > chi8 > beta1`.

### `KLTableDarkPhoton`

**Dark-photon plasma** (`hpc_results/29205968/stage_a_amp_v3_pub`; N_eff amp/sup ≈ 7443/10675)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| mA2 | 0.50 | 0.50 | 1.25 | 1.25 | 1.45 |
| deltam | 2.77 | 0.17 | 2.98 | 1.46 | 3.38 |
| xi | 0.85 | 0.85 | 1.24 | 1.24 | 0.00 |
| alpha3 | 0.66 | 0.66 | 0.67 | 0.67 | 0.49 |

Top-3 by D+ (amp marginal): archived `deltam > xi > alpha3` → corrected `xi > alpha3 > mA2`.

### `KLTableEH (commented out in the archive)`

**Euler-Heisenberg** (`hpc_results/29700083/eh_gert_amp_v1`; N_eff amp/sup ≈ 104/103)

| coupling | D+ archived | D+ corrected | D− archived | D− corrected | D₊‖₋ (unchanged) |
|---|---|---|---|---|---|
| rho | 2.64 | 0.15 (≤ floor) | 2.70 | 0.24 | 0.45 |
| sigma | 2.41 | 0.19 | 2.64 | 0.34 | 1.66 |

Top-3 by D+ (amp marginal): archived `rho > sigma` → corrected `sigma > rho`.

### Rescue-chain status (GH #433)

T7 (chi_closure), T9 (xi_kinetic_closure), the 17D NP control, **and T6-full
(parity_odd_full — the worst case, omitted from the first revision of this
record: N_eff = 6.0 amp / 13.8 sup, floors 3.25/1.42 nats, all 20 of 20
parameters floor-dominated in both directions)** have effective sample sizes of
6–248, i.e. per-parameter noise floors of ~0.08–3.25 nats: **no per-coupling
marginal from these chains is quotable as signal**, before or after the #420
fix. (The Euler-Heisenberg chain is also floor-limited at N_eff ≈ 104 — floors
≈ 0.19 — which is consistent with, and part of, its null verdict.) This sharpens the archive's own hedges (T9 was already marked "not
quantitatively meaningful"). Re-runs with the overflow mitigations are explicitly
deferred to a future campaign; the fixed pipeline now flags the condition
automatically (`floor_dominated_params` in the consistency block, `(≤ floor)`
annotations in the CLI table).

### Provenance notes

- **Barker amp** (`29507332/d21_barker_amp_v3_pub`): its `inference.json` has NO
  `priors` key (`reconstructed: true`), so the first revision of this record
  scored ξ against a fabricated arctan prior (0.92); under the true log-uniform
  prior (sup chain's recorded priors + the campaign sbatch template) ξ = 0.83.
  Rankings unaffected (χ stays top). `parameter_importance.json` now records
  `priors_provenance` per direction, and the recompute warns on every
  fabrication (post-merge review, #434 follow-up).

- Chains without recorded priors (pre-schema) are scored against priors
  reconstructed from the committed sbatch templates (T7/T6: log-uniform ξ); a
  one-time check against HPC submit logs would make this airtight
  (`docs/dkl_recompute_report.md` caveats).
- The ricci_em chain used a one-sided broad β₃ prior (GH #387): the β₃ row of
  `KLTableEMRCParityPair` is not comparable across the parity pair, independent
  of this amendment.
- Staleness tell for anything outside the 13 recomputed chains: a
  pre-v0.48.8 `importance.json` has no `consistency` block; recompute before
  quoting (`tidal analyze <dir> --inference --importance`).

### For a future campaign

The estimator is fixed and self-checking in the pipeline (v0.48.8+); the
amendments above are analysis-level only — no simulation was re-run. A future
campaign should: re-run the floor-limited chains with overflow mitigations before
making per-coupling claims there; re-examine the YM NP-control marginal-matching
claim with healthy chains; use matched β priors for the parity pair (GH #387);
and quote D+/D− only alongside the `consistency` block (N_eff, floors).

### Trail

GH #420 (estimator bug) · #425 (arctan bounds) · #426/#431/#435 (fix, recompute,
hardening PRs) · #432 (this propagation) · #433 (noise floor) · #434 (fabricated
priors) · `docs/dkl_recompute_report.md` (full evidence, per-chain tables,
old-vs-new for both amp and sup directions).

---

## Amendment 2 (2026-08-27) — Localized E.cal calibration on the corrected operator (GH #449)

### What happened

The localized (position-dependent background) solver path carried a chain of
defects when the Phase E calibration was recorded: the rfft convolution basis
(#445), the missing velocity-row mass scale on deferred constraint terms (#444),
and the inter-constraint time-derivative mishandling of the ungauged-gravity
class (#457, with #458/#459/#460 on the same path). All are fixed (v0.49.x,
branch `feat/ws2-localized-path-audit`). The corrected operator exposed that the
FULL ungauged localized pencil cannot be evolved faithfully in double precision —
the far-field background restores the linearized-diffeomorphism freedom below
machine precision, and the probe-based gauge quotient discarded independent
equations there (#474). The observable-sector closure route (#468) therefore
evolves the exactly closed {h_5, a_1} sector on its own (aa8a3b38), which is the
sector the Phase E conversion observable reads.

### The corrected number

Frozen Phase E geometry (`scripts/hpc_submit_drafts/v3e_localised/_geometry.env`:
L = 100, N = 128, periodic; dual Gaussians at 25/75 with σ_B = 5, Bpeak = 0.01;
wavepacket on h_5 at x_c = 8, σ_w = 3, k = 2, h0 = 0.01; κ = 1):

| quantity | archived (PHASE_E_TRACKER 1.1, 2026-05-24) | corrected (2026-08-27) | status |
|---|---|---|---|
| P_peak(h_5 → a_1) | 0.0036 | **0.003912** (t_end = 40) | superseded |
| Boccaletti sin²(κ·Bpeak·σ_B·√(2π)/2) | 0.0039 | 0.003922 | — |
| ratio P / sin² | ~0.92 | **0.9975** | strengthened |
| t_end independence | not checked | P_peak = 0.003914 at t_end = 80 (A(80)/A(40) = 1.0006) | new |
| final P after the sign-flipped second Gaussian (t = 80) | — | 7.8e-6 (phase cancels, ∫B dz = 0 over the box) | new |

The archived calibration claim ("matches Boccaletti within ~10%") **stands** in
substance and is **strengthened**: the corrected operator reproduces the
path-integrated Boccaletti prediction to a quarter of a percent, with no
resolution- or t_end-dependent growth (the #455 spurious-operator signature is
absent on the closure route).

### What is affected, and what is not

- The thesis manuscript's results rest on the uniform-background plane-wave
  benchmark (Amendment 1 context; `results.tex`), which never ran the localized
  path — nothing in the archive is superseded beyond the internal Phase E
  bookkeeping number above.
- The E.T2 amplification run (HPC 29640051) was produced on the defective
  operator and is re-judged only when the localized programme resumes.
- Torsion-roster localized specs need the same re-measurement; the closure fact
  `closure({h_5, a_1}) == {h_5, a_1}` is pinned for E.cal only
  (`tests/test_spec_restriction.py`) and must be verified per spec.
- Total energy on closure-restricted runs is deliberately unavailable (81
  Hamiltonian terms touch omitted fields); conversion is the exact observable.

### Trail

GH #445 · #444 · #457/#458/#459/#460 (defects) · #468 (closure route, aa8a3b38) ·
#474 (quotient root cause) · #473 (general staircase reduction, future) · #449
(this quantification) · `docs/PHASE_E_TRACKER.md` 1.1 (superseded line annotated).
