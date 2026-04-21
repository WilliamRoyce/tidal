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

### Stage A: Dark-Photon-Plasma 4D nested sampling

- [x] Stage 0 gate passed for T1
- [x] HPC amplification job submitted (job ID: 28134330, INTR, RUNNING cpu-q-580)
- [x] HPC suppression job submitted (job ID: 28134387, standard QOS, RUNNING cpu-q-182)
- [ ] Results pulled
- [ ] Analysis: D_KL(params) computed — threshold > 0.05 nats
- Notes:

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

---

## Key Findings

*(Filled in as results arrive)*

- Stage A: T1 suppression D_KL=0.017 (weak signal, MAP at alpha3→0 boundary).
  Amplification run still in progress (28134330, 52min).
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
