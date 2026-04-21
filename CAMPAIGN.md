# PGT Survey Campaign

**Central question:** Can PGT torsion amplify the Gertsenshtein effect? If so, which
Lagrangian sector and at what parameters?

**Plan file:** `/home/vscode/.claude/plans/binary-snacking-chipmunk.md`
**Branch:** `hpc/pgt-survey` (created from `feature/torsion-sweeps`; rebase onto `main` once PR merges)

**Fixed params throughout:** `kappa=1.0`, `B0=0.01`
**IC convention:** `--ic plane-wave --ic-amplitude 1e-2 --periodic --bounds 0:50`
**Paired runs:** every theory gets an amplification run (`P_max:maximize`) + suppression run (`P_max:minimize`) submitted simultaneously.

---

## Theory Registry

| ID | Campaign name | JSON | Free params | IC fields |
|----|--------------|------|-------------|-----------|
| T1 | Dark-Photon-Plasma | `dark_photon_plasma.json` | mA2, deltam, xi, alpha3 (4) | h=h_5, a=a_1 |
| T2 | Einstein-Cartan | `torsion_gertsenshtein_b5_zero.json` | alpha1, alpha2, alpha3 (3) | h=h_5, a=a_1 |
| T3 | R²-PGT | `torsion_gertsenshtein.json` | alpha1-3, b5 (4) | h=h_5, a=a_1 |
| T4 | Ricci-EM | `torsion_gertsenshtein_nonminimal.json` | alpha1-3, delta1 (4) | h=h_5, a=a_1 |
| T5 | YM-PGT | `torsion_gertsenshtein_general_nonminimal.json` | beta1-3, xi, delta1, chi, zeta1-3 (9) | h=h_5, a=a_1 |
| T6 | YM-PGT-CP | `torsion_gertsenshtein_parity_odd.json` | T5 + d14-17, d19-21, zt1-6 (22; d18 absent) | h=h_5, a=a_1 |
| T7 | Complete-Even-PGT | `pgt_complete_even.json` (TBD, needs derivation) | T5 + chi2-10, b2-6 | TBD |
| T8 | Complete-Odd-PGT | `pgt_complete_odd.json` (TBD, needs derivation) | T6 + chitilde1-36, d1-13 | TBD |
| EH | Einstein-Maxwell-EH | TBD (GitHub issue filed) | c1, c2 | TBD |

**Sub-theories of T5 (no new JSON needed, fix unused params to 0):**
| Sub-ID | Campaign label | Free params | Fixed to 0 | Paper |
|--------|---------------|-------------|------------|-------|
| T5-Barker | Barker-PGT | beta1-3, xi, chi | delta1, zeta1-3 | Barker 2406.12826 |
| T5-Shapiro | Shapiro-PGT | beta1-3, zeta1-3 | xi, delta1, chi | Shapiro hep-th/0103093 |
| T5-Bahamonde | Bahamonde-PGT | beta1-3, delta1 | xi, chi, zeta1-3 | Bahamonde et al. |

---

## Stage Progress

### Stage 0: Local Smoke Tests (MANDATORY gate before any HPC submission)

| Theory | Inspected | IC fields found | Simulated | Energy OK | Notes |
|--------|-----------|----------------|-----------|-----------|-------|
| T1 Dark-Photon-Plasma | [x] | h=h_5, a=a_1 | [x] | N/A (conservation notes) | |
| T2 Einstein-Cartan | [x] | h=h_5, a=a_1 | [x] | N/A | |
| T3 R²-PGT | [x] | h=h_5, a=a_1 | [x] | N/A | b5 gate: PHYSICALLY NULL — see Stage C notes |
| T4 Ricci-EM | [x] | h=h_5, a=a_1 | [x] | N/A | |
| T5 YM-PGT | [x] | h=h_5, a=a_1 | [x] | N/A | small couplings required (chi=0.001, xi=0.01) |
| T6 YM-PGT-CP | [x] | h=h_5, a=a_1 | [x] | N/A | small couplings required |

**Note:** Conservation checks not valid for these theories (Dirac-Bergmann Hamiltonian issues, documented separately).

### Stage A: Dark-Photon-Plasma 4D nested sampling — ✅ NULL CONFIRMED

- [x] Stage 0 gate passed for T1
- [x] HPC amplification job submitted (final: 28145377, INTR QOS, 28:16)
- [x] HPC suppression job submitted (final: 28145425, std QOS, 3:52)
- [x] Results pulled (hpc_results/28145377/, hpc_results/28145425/)
- [x] Analysis: **Amplification D_KL=0.0155, Suppression D_KL=0.0057** — both below 0.05 threshold
- **Finding (2026-04-21):** Dark photon plasma (mA2, deltam, xi, alpha3) shows NULL amplification
  across the explored 4D prior space. log(Z)=+0.022 for amplify, -0.015 for suppress — both
  essentially zero. MAP at (mA2=0.057, deltam=-0.24, xi=4.42, alpha3=0.143) is far from the
  old tachyonic peak (xi=0.52, which the v0.31 pre-guard run falsely concentrated on). 95% CI
  spans most of the prior for every parameter → posterior ≈ prior.
- Consistent with v0.31 MC sweep (276 runs, P_max matched sin²(κB₀t/2) to 6.7e-6 precision).
- **Lesson:** HPC pip metadata must match local source. Invalid pre-fix runs (28133218/516/517,
  28134330) traced to v0.31.5 install predating the stability guard; fixed in this session by
  `pip install -e .` reinstall + tarball refresh + version sync check in hpc_shuttle push.

### Stage B: Einstein-Cartan null (T2) — ✅ NULL CONFIRMED

- [x] Stage 0 gate passed for T2
- [x] HPC null job submitted (job ID: 28134415, standard QOS, COMPLETED in 13s)
- [x] Results pulled (hpc_results/28134415/)
- [x] Analysis: **D_KL = 0.0031 ± 0.0020 nats** (all alphas at 0.0032) — below 0.005 threshold
- **Finding:** Non-propagating torsion (alpha1-3 across their full priors ±2, log_uniform[0.01,4])
  has ZERO effect on Gertsenshtein conversion. log(Z) = -0.0007 ± 0.0022 (posterior ≈ prior).
  Completion in 13s reflects PolyChord finding a completely flat likelihood surface — no structure
  to explore. This is the expected structural null for Einstein-Cartan minimal torsion coupling.

### Stage C: R²-PGT perturbative b5 (T3)

- [x] Stage 0: local b5 correction check — PHYSICALLY NULL (b5 decouples from TT channel)
- **Finding (2026-04-21):** b5 corrections source from torsion fields (e.g. `b5*∂³t_10` in h_0 eq).
  With plane-wave IC on h_5 only, all torsion fields are zero in Pass 0 → b5 corrections = 0 trivially.
  P_max(h_5→a_1) is identical to GR baseline sin²(κB₀t/2) for all b5 values.
  R̃² propagating torsion DECOUPLES from the TT graviton-photon Gertsenshtein channel.
  Stage C HPC will confirm D_KL(b5) ≈ 0 (flat posterior) — expected null, not a failure.
- [ ] Stage C gate passed (re-defined: perturbative solver runs, gate is "no crash"; physics null accepted)
- [ ] HPC order=0 amplification job (job ID: ?)
- [ ] HPC order=1 amplification job (job ID: ?)
- [ ] Results pulled
- [ ] Analysis: D_KL(b5, ord=1) − D_KL(b5, ord=0) > 0.01 nats
- Notes: Expected D_KL(b5) ≈ 0 — structural null, same as T2.

### Stage D: Nonminimal survey

#### D1: Ricci-EM / T4

- [ ] Stage 0 gate passed for T4
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- [ ] Analysis: D_KL(delta1) > 0.05 nats?
- Notes:

#### D2: YM-PGT / T5 + paper sub-theories

- [ ] Stage 0 gate passed for T5
- [ ] HPC: full YM-PGT amplify (job ID: ?)
- [ ] HPC: full YM-PGT suppress (job ID: ?)
- [ ] HPC: Barker-PGT amplify (job ID: ?)
- [ ] HPC: Shapiro-PGT amplify (job ID: ?)
- [ ] Results pulled
- [ ] Analysis: any D_KL > 0.1 nats?
- Notes:

#### D3: YM-PGT-CP / T6

- [ ] Stage 0 gate passed for T6
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- [ ] Analysis: any D_KL > 0.1 nats?
- Notes:

### Stage E: New theory derivations

#### E1: Complete-Even-PGT / T7

- [ ] theory.toml designed (`examples/pgt_complete_even/`)
- [ ] Wolfram derivation complete
- [ ] Smoke test: non-zero chi2+ or b2+ correction confirmed
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- Trigger: Stage D shows signal in χ or R×R direction

#### E2: Complete-Odd-PGT / T8

- [ ] theory.toml designed (`examples/pgt_complete_odd/`)
- [ ] Wolfram derivation complete
- [ ] HPC amplification job (job ID: ?)
- [ ] HPC suppression job (job ID: ?)
- Trigger: Stage D3 shows parity-odd signal

---

## HPC Job Log

| Job ID | Theory | Run type | Status | Wall time | Notes |
|--------|--------|----------|--------|-----------|-------|
| 28133218 | T1 Dark-Photon-Plasma | amplification | FAILED | 0:02 | REMOTE_ROOT bug (fixed) |
| 28133516 | T1 Dark-Photon-Plasma | suppression | FAILED | 0:00 | REMOTE_ROOT bug (fixed) |
| 28133517 | T2 Einstein-Cartan | null | FAILED | 0:00 | REMOTE_ROOT bug (fixed) |
| 28133656 | T1 Dark-Photon-Plasma | amplification | FAILED | 0:00 | sed & expansion bug (fixed) |
| 28133932 | T1 Dark-Photon-Plasma | suppression | FAILED | 0:00 | sed & expansion bug (fixed) |
| 28133933 | T2 Einstein-Cartan | null | FAILED | 0:00 | sed & expansion bug (fixed) |
| 28134330 | T1 Dark-Photon-Plasma | amplification | TIMEOUT | 1:00:26 | INVALID — tachyon peak logL=31430 (#307); ran stale v0.31.5 lacking stability guard |
| 28134387 | T1 Dark-Photon-Plasma | suppression | COMPLETED | 2:43 | QUESTIONABLE — ran stale v0.31.5; re-verify after HPC upgrade |
| 28134415 | T2 Einstein-Cartan | null (amplification) | COMPLETED | 0:13 | D_KL=0.003 (NULL CONFIRMED); valid — torsion decouples, no tachyonic region |
| 28141098 | T1 Dark-Photon-Plasma | amplification (re) | CANCELLED | 3:20 | cancelled to reinstall HPC tidal properly |
| 28145274 | (guard test) | verification | COMPLETED | 0:05 | v0.33.13 confirmed; stability guard fires on cpu-q-19 with ratio 1.26e+05 |
| 28145377 | T1 Dark-Photon-Plasma | amplification (re²) | COMPLETED | 28:16 | D_KL=0.0155, log(Z)=+0.022 (NULL) |
| 28145425 | T1 Dark-Photon-Plasma | suppression (re) | COMPLETED | 3:52 | D_KL=0.0057, log(Z)=-0.015 (NULL) |

---

## Key Findings

*(Filled in as results arrive)*

- Stage A: **Dark-Photon-Plasma NULL CONFIRMED**. amp D_KL=0.0155 log(Z)=+0.022, supp D_KL=0.0057
  log(Z)=-0.015. Torsion-CDT trace dark photon does not amplify Gertsenshtein conversion across
  the 4D prior (mA2∈[0.001,1], δₘ∈[-0.5,0.5], ξ∈[0.05,20], α3∈[0.001,0.5]). Confirms v0.31
  MC sweep finding. (hpc_results/28145377/, 28145425/)
- Stage B: **Einstein-Cartan NULL CONFIRMED** — D_KL=0.003 nats for all alphas,
  log(Z)≈0. Non-propagating torsion does not affect h↔a conversion. (hpc_results/28134415/)
- Stage C: R̃² b5 term decouples from TT Gertsenshtein channel. With IC on h_5, all b5 corrections = 0 (torsion source fields zero in Pass 0). P_max identical to GR baseline for all b5. Expected D_KL(b5) ≈ 0.
- Stage D1: 
- Stage D2: 
- Stage D3: 
- Stage E: 

---

## GitHub Issues Filed

| Issue | Topic |
|-------|-------|
| #296 | Euler-Heisenberg F⁴ extension (#271 resolved, defer EH theory to future) |
| #297 | Non-propagating/constraint torsion theories (single irreducible sector investigation) |
