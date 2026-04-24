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

### Stage A: Dark-Photon-Plasma 4D nested sampling — ⚠️ ALL 3 OLD RUNS WRONG-REGIME; REDOING UNDER NEW SIGN CONVENTION

**2026-04-24 convention flip:** The CDT Lagrangian was changed from `+alpha3·I3` →
`-alpha3·I3` so that `alpha3 > 0` now corresponds to the **stable Proca dark photon**
regime (physical mass² = +2·alpha3 > 0). Previously `alpha3 > 0` was the tachyonic-
spatial-trace regime, which doesn't match the "torsion as a dark photon" physics intent.
All three historical Stage A runs used `alpha3 = log_uniform(0.001, 0.5) > 0` under the
OLD convention, so they sampled the tachyonic regime. They are **archived as valid
simulations of the wrong physics regime** and are being replaced by new-convention reruns.

- [x] Stage 0 gate passed for T1
- [x] HPC amplification 28226826 (amplify, old convention) — archived; was a NULL in tachyonic regime
- [x] HPC suppression 28216072 (suppress, old convention; pre-IC-snap) — IC-leakage invalidated, archived
- [x] HPC suppression 28365129 (suppress retry 1, old convention; post-IC-snap) — FAILED on /tmp-full on compute node cpu-q-553, archived
- [x] HPC suppression 28366464 (suppress retry 2, old convention; post-IC-snap) — completed 3:33, pulled, archived as wrong-regime. log(Z)=-0.069, D_KL(xi)=0.23
- [ ] HPC amplify rerun under NEW convention (pending)
- [ ] HPC suppress rerun under NEW convention (pending)

- **Amplify finding (2026-04-22, 28226826):** log(Z)=+0.118±0.006, joint D_KL=0.043,
  D_KL(xi)=0.226, D_KL(mA2)=0.032. MAP at (mA2=0.34, deltam=-0.21, xi=1.08, alpha3=0.054).
  Posterior concentrates at small mA2 and xi≲1 (stability boundary); within the stable region,
  no parameter combination produces P_max > P_GR. **Amplify null is valid.**

- **Suppress finding superseded (2026-04-22, 28216072):** the 0/2319 count was NOT physical —
  investigation traced it to plane-wave IC spectral leakage. `--ic-wavevector 2.0` on a
  `bounds 0:50 --grid-shape 32` grid leaks amplitude onto every discrete Fourier mode
  (the continuous `cos(2.0·x)` is not periodic on [0, 50]). With `deltam ≠ 0` the Gertsenshtein
  block merges with torsion-trace sub-blocks containing redundant TorsionCDT ghost modes;
  those ghosts are tachyonic at low k (`Re(λ) ≈ 0.95` at k=0) and the leakage × exp growth
  trips the modal divergence guard for most of the interesting 4D prior. PolyChord then
  samples only the narrow stable fraction (alpha3 → 0 decoupling limit) and reports the
  trivial P_max = P_GR null. This is a discretization artefact, not physics.

- **Fix applied (this session, 2026-04-22):** `tidal/cli/_simulate.py::_plane_wave_slots`
  now auto-snaps `--ic-wavevector` to the nearest discrete Fourier mode on periodic axes
  (clamped below Nyquist). Verified empirically: at (mA2=0.955, deltam=0.01, xi=0.274,
  alpha3=0.123), the `--ic-wavevector 2.0` case now succeeds (snapped to 1.885) and gives
  P_max ≈ 0.00176 — the genuine ~97% plasma suppression. Cross-check: FV model (10-field
  formulation, no ghost sector) **agrees bit-exactly with TorsionCDT** under the correct
  parameter mapping `mT2 = -2·alpha3` (relative Δ = 2.8e-14, verified 2026-04-24 and
  documented in Issue #318). The initial "1.2% tension" reported earlier this session
  was a sign error in the theory.toml equivalence-comment; the modal solver is
  machine-precision as designed. Docs: `docs/tex/plane_wave_ic.tex`,
  `examples/dark_photon_plasma/theory.toml` (line 18-28 for the corrected equivalence map).

- **Lesson:** HPC pip metadata must match local source. Invalid pre-fix runs (28133218/516/517,
  28134330) traced to v0.31.5 install predating the stability guard; fixed in an earlier session.
- **Retired t_end=10 runs:** 28145377 (amp, D_KL=0.0155), 28145425 (sup, D_KL=0.0057) —
  hpc_results/ dirs deleted as part of this session's cleanup.

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
| 28145377 | T1 Dark-Photon-Plasma | amplification (re², t_end=10) | COMPLETED | 28:16 | D_KL=0.0155, log(Z)=+0.022 — pre-oscillation regime, plasma-mass suppression invisible |
| 28145425 | T1 Dark-Photon-Plasma | suppression (re, t_end=10) | COMPLETED | 3:52 | D_KL=0.0057, log(Z)=-0.015 — same regime caveat |
| 28215825 | T1 Dark-Photon-Plasma | suppress (tend50 attempt 1) | CANCELLED | 1:02 | INTR; cancelled by accident while targeting amplify |
| 28215827 | T1 Dark-Photon-Plasma | amplify (tend50 attempt 1) | CANCELLED | 1:43 | std QOS; cancelled — should have been INTR |
| 28216041 | T1 Dark-Photon-Plasma | amplification (tend50) | TIMEOUT | 1:00:27 | INTR; PolyChord global log(Z)=+0.106±0.005 converged but 12/15 clusters still active at 1h wall |
| 28216072 | T1 Dark-Photon-Plasma | suppression (tend50) | INVALID | 0:05:18 | IC spectral-leakage bug; most samples rejected by divergence guard. Superseded by post-IC-snap rerun (pending). |
| 28226826 | T1 Dark-Photon-Plasma | amplification (tend50, resume) | COMPLETED | ~2h | std QOS; D_KL=0.043, log(Z)=+0.118±0.006, D_KL(xi)=0.226, MAP(mA2=0.34,xi=1.08,a3=0.054). Amplify null still valid after IC-snap fix. |
| 28365129 | T1 Dark-Photon-Plasma | suppression rerun with IC-snap fix | FAILED | 0:00:16 | Venv tarball extraction to compute-node /tmp failed with "no space left on device" on node cpu-q-553. Node-specific; resubmitted as 28366464. |
| 28366464 | T1 Dark-Photon-Plasma | suppression rerun with IC-snap fix (retry) | SUBMITTED | — | INTR QOS, 1h wall; same command as 28365129, hoping for a different node. |

---

## Key Findings

*(Filled in as results arrive)*

- Stage A: **Amplify NULL confirmed; Suppress rerun pending IC-snap fix.**
  - **Amplify (28226826)** — joint D_KL=0.043, log(Z)=+0.118±0.006, D_KL(xi)=0.226, D_KL(mA2)=0.032.
    MAP at (mA2=0.34, δₘ=-0.21, ξ=1.08, α3=0.054). No P_max > P_GR in the stable region. Valid.
  - **Suppress (28216072) INVALIDATED (2026-04-22).** The 0/2319 count with `P_max < P_GR` was
    NOT physics — root cause is plane-wave IC spectral leakage on the periodic grid. Off-grid
    `--ic-wavevector 2.0` leaks amplitude onto every discrete Fourier mode; at mA2·deltam
    combinations where the merged Gertsenshtein-torsion block has low-k tachyons (redundant
    TorsionCDT directions, Re(λ)≈0.95 at k=0), the leaked amplitude is exponentially amplified
    and trips the divergence guard, leaving PolyChord with only the trivial alpha3→0
    decoupling slice of the prior. **Fix applied**: `tidal/cli/_simulate.py::_plane_wave_slots`
    now auto-snaps `--ic-wavevector` to the nearest discrete Fourier mode on periodic axes
    (pass `--ic-no-snap` for legacy). Direct `tidal simulate` at
    (mA2=0.955, δₘ=0.01, ξ=0.274, α3=0.123) with the fix gives P_max≈0.00176 — the expected
    97% plasma suppression below P_GR=0.0612. **FV ↔ TorsionCDT equivalence verified
    bit-exact** (Δ = 2.8e-14) after re-deriving FV from a fresh theory.toml and applying
    the correct parameter mapping `mT2 = -2·α₃` (Issue #318 closed 2026-04-24; the
    initial 1.2% tension was a sign error in `dark_photon_plasma/theory.toml:19`,
    now corrected to read `mT2 = -2·alpha3`).
  - **Suppress rerun** submitted with full original 4D prior; awaiting results.
- Stage B: Einstein-Cartan (T2). Joint D_KL=0.003, log(Z)≈0 — but **corrected marginals show
  structure**: α1=0.11, α2=0.06, α3=0.07 (all informative). Joint ≪ sum because the posterior
  is broad on a 2D stability ridge in 3D prior. P_max flat on ridge (null amplification) but
  the ridge itself is smaller than the prior — stability guard excludes unstable α combinations.
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
