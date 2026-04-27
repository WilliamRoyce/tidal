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
- [x] HPC amplify rerun under NEW convention: **job 28367920** submitted 2026-04-24, std QOS, 3h wall (**COMPLETED 20:01**)
- [x] HPC suppress rerun under NEW convention: **job 28367934** submitted 2026-04-24, INTR QOS, 1h wall (**COMPLETED 2:19**)
- [x] Results pulled for both
- [x] Two bugs found during analysis (now fixed); results re-interpreted

### Bugs identified and fixed during 2026-04-24 post-run analysis

1. **Issue #319**: `parse_likelihood()` silently dropped `--baseline-formula` for
   `maximize`/`minimize` types (only `extremize` propagated it). Result: HPC
   `log_likelihood` field was raw `P_max` (or `−P_max` for minimize), not
   `log(P_max/P_GR)` as documented. Fixed in commit `d34a204`. v3 results need
   re-interpretation under "raw P_max" semantics.
2. **Issue #320**: Modal solver eigendecomposition is **ill-conditioned at high-xi
   resonant-instability regions** of the dark photon plasma model. At "the same"
   parameter point (precision-equivalent inputs), max `Re(λ)` can vary by 4× (e.g.
   168 → 316 → 559) and the divergence pre-check fires inconsistently. v3 amplify's
   "high-amplification" top samples lie precisely in these ill-conditioned regions,
   so their reported `P_max` values are not trustworthy as physics.

### Re-interpreted results (under raw `P_max` semantics, post-#319 understanding)

- **Suppress 28367934 (real signal, robust)**:
  - Median `P_max ≈ 0.0036 = 6%` of P_GR (94% suppression).
  - **2283 / 2294 (99.5%) samples with `P_max < P_GR`** — strong genuine suppression.
  - Posterior concentrates at small `alpha3` (decoupling limit) — well-conditioned region.
  - D_KL(xi)=0.22 dominates (stability cutoff); other D_KLs ~0.02.
  - **This is the first non-null Stage A result**. The claim "0/2294" reported
    pre-fix was wrong — it came from inverting the wrong likelihood formula.

- **Amplify 28367920 (high-xi signal NOT trustworthy — Issue #320)**:
  - 41% of samples have `P_max > P_GR`, 18% have `P_max > 1` (unphysical).
  - Top samples cluster at `xi ∈ [5, 8]` (high-end of `log_uniform(0.05, 20)`),
    `|deltam| ≈ 0.5` (prior boundary), `alpha3` and `mA2` variable.
  - Local reproduction of the top-1 sample (`mA2=0.00295, deltam=0.054, xi=5.48,
    alpha3=0.0098`) showed `tidal simulate` (disk path) DIVERGES at the modal
    solver's pre-check, while `tidal sample` (in-memory path) returns finite
    `P_max ≈ 1.7-2.0` — a path-dependent contradiction that Issue #320 traces to
    eigenvalue ill-conditioning at resonant-instability windows.
  - Conclusion: the "amplification signal" at high xi is numerical-noise
    eigenvalue amplification, NOT real physics. v3 amplify cannot be used as-is.

- [ ] Submit Stage A v4 on Tier 1 + 1.5 modal solver (path-D + augmented-exp Pass 1, v0.31+)
  - Tier 1 (committed v0.31+): `_evolve_per_mode` retired the eigendecomposition default
    in favour of unconditional `scipy.linalg.expm(M·dt)` precompute + matvec (Higham 2009
    Padé scaling-and-squaring). Robust for arbitrary `cond(V)`; benchmark wall time
    0.49-1.06× of eigendecomposition on the campaign workload envelope. Closes the
    eigenvector ill-conditioning failure mode (#320 root cause).
  - Tier 1.5 (also v0.31+): `_evolve_duhamel_per_mode` rewritten to use the augmented
    matrix exponential `exp(t·[[A, S], [0, A]])·[0; y₀]` (Al-Mohy & Higham 2011 §5.2).
    Pass 0 + Pass 1 share a single canonical robust backend; no eigendecomposition
    anywhere in the modal solver's per-mode path.
  - Verification: 2044 tests pass, including 6 new robustness tests in
    `tests/test_modal_robust_evolution.py` (synthetic ill-conditioned 4×4, Padé
    machine-precision against mpmath reference, augmented-exp Pass 1 vs Duhamel
    kernel, Pass 1 IC-zero check, CDT path D smoke).
  - v4 plan: `B₀=0.01, t_end=10, snapshots=2`, `--nlive 1000 --num-repeats 10
    --precision-criterion 0.005`, corrected `--baseline-formula 'sin(kappa*B0*t_end/2)**2'`
    (post-#319 parser fix). Same 4D prior as v3 — no constraint-based tightening
    needed because path D + augmented-exp make CDT numerically robust everywhere.
  - Cross-validation: top-amplification samples should agree with FV scan to
    `Δ/P_max < 1e-9`. Real expected range per FV scan: A ∈ [0.02, 1.6],
    P_max ∈ [0.001, 0.10]. v3 amplify (28367920) "P_max up to 2.0" was CDT
    ill-conditioning contamination — superseded.

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

- [x] Stage 0 gate passed for T4
- [x] HPC amplification job: 28520217 (INTR; 28519213 first attempt — files vanished post-completion, mystery; resubmitted as v5b)
- [ ] HPC suppression job: 28519675 (standard QOS, PD in queue)
- [x] Analysis: **D_KL(delta1) = 1.679 nats — overwhelmingly the dominant constraint (96% of joint D_KL)**
- Notes:
  - **Stability-guard refactor (post #322, commit c10aa8a)** was a precondition; old eigenvalue+pinv path 100% rejected T4.
  - Priors: α₁ ∈ uniform[-1, 1], α₂ ∈ uniform[-2, 2], α₃ ∈ log_uniform[0.05, 2], δ₁ ∈ uniform[-2, 2].
  - Padé probe 15.2% rejection rate across the prior; mostly at α₂ < -1 (51% at α₂=-2) — consistent with the lower analytic instability boundary at α₂ = -7/(4κ²) = -1.75.
  - **Amplify (28520217) verdict**: log Z = -2.26 ± 0.07 (model strongly disfavoured vs null), joint D_KL = 1.79 ± 0.06 nats, ESS = 1356.
  - **Marginal D_KL**: δ₁ = 1.68 nats (dominant), α₁ = 0.05, α₂ = 0.07, α₃ = 0.01. The δ₁ coupling carries virtually all the information; the PGT torsion-mass parameters are barely constrained.
  - **A range in chain**: [4.95×10⁻⁹, 1.255]. Maximum amplification 1.26 — **does NOT match `docs/AMPLIFICATION_INVESTIGATION.md`'s A_coupling ≈ 1.9 Schur-complement prediction**. Top samples cluster at small |δ₁| ∈ [0.07, 0.15], α₂ ∈ [-1.9, -0.8]; the predicted optimum at δ₁ ≈ 1, α₂ ≈ -1 was sampled (15% rejection rate at δ₁=1) but the inference did not select it.
  - **A range minimum (4.95×10⁻⁹) FAR EXCEEDS the predicted suppression valley (10⁻³)** at δ₁ ≈ ±1.3, varying α₂ — suggests the destructive-interference mechanism is even stronger than the linearised Schur complement estimate.
  - **Verdict**: T4 is a **strong suppressor** of Gertsenshtein conversion (avg A ≈ exp(-2.26) ≈ 0.10), with a **modest amplification ceiling** (max A = 1.26) — physical mechanism dominated by the δ₁ R̃ₘᵤᵥ Fᵘᵥ coupling.  Suppress (28519675) pending.

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
| 28366464 | T1 Dark-Photon-Plasma | suppression rerun with IC-snap fix (retry) | COMPLETED (WRONG REGIME) | 0:03:33 | INTR QOS, ran fine on a different node. log(Z)=-0.069, D_KL(xi)=0.235, D_KL(alpha3)=0.037. Archived — OLD Lagrangian sign convention, sampled tachyonic regime not stable Proca. |
| 28367920 | T1 Dark-Photon-Plasma | amplification, NEW Lagrangian sign convention (stable Proca) | RUNNING | ≥5:12 | std QOS, 3h wall. Replaces 28226826 semantically: alpha3 = log_uniform(0.001, 0.5) now sweeps stable-Proca regime (m² = +2·alpha3 > 0). |
| 28367934 | T1 Dark-Photon-Plasma | suppression, NEW Lagrangian sign convention (stable Proca) | COMPLETED | 0:02:19 | log(Z)=-0.056±0.004, 2294 samples, 100% success. 0/2294 with P_max<P_GR (posterior concentrates at alpha3→0 decoupling limit). But top-P_max samples reach 5.88× P_GR — amplification signal. D_KL(xi)=0.22 drives stability. |
| 28418115 | T1 Dark-Photon-Plasma | amplify v4 — Tier 1 + 1.5 fixed solver, B0=0.01 t_end=10 snapshots=2 nlive=1000 num_repeats=10 prec=0.005 ntasks=76 | PENDING | — | Submitted 2026-04-25 20:24Z, std QOS, 6h wall, polychord_standard.sbatch. Replaces v3 28367920 whose "P_max up to 2.0" was CDT eigenvector ill-conditioning artefact (#320). v0.36.0 modal solver uses unconditional `scipy.linalg.expm` Padé + augmented-exp Pass 1 — robust for arbitrary cond(V). |
| 28418421 | T1 Dark-Photon-Plasma | suppress v4 — Tier 1 + 1.5 fixed solver, same params as amplify but minimize | TIMEOUT | 1:00:19 | INTR wall hit before convergence. **partial chains valid**: 2030 dead samples, anesthetic gives log Z (partial) = +0.523, D_KL (partial) = 1.12 nats — same magnitude as amplify, very different from v3 suppress's near-trivial posterior. Confirms the post-fix landscape is **structurally informative** (not concentrated at α₃→0 decoupling). v4 nlive=1000 num_repeats=10 was too aggressive for INTR. Resubmitted as 28461922 at v3 resolution. |
| 28461922 | T1 Dark-Photon-Plasma | suppress v4 RESUBMIT — v3 resolution (nlive=400, num_repeats=5, prec=0.01) on Tier 1+1.5 solver | COMPLETED | 0:38 | Submitted 2026-04-26 09:47Z, INTR. log Z = +0.66 ± 0.05, D_KL = 1.95 nats. Posterior concentrates at decoupling corner (small alpha3). Ghost contamination still present without stability guard — superseded by v5. |
| 28474676 | T1 Dark-Photon-Plasma | **amplify v5** — stability-guard fix (commits e361113, 3855fc1, 7ef182d): pre-flight tachyonic eigenvalue check in `_evaluate_likelihood` rejects Re(λ) > 0.3 samples before simulation. Includes per-sample run_status metadata and post-hoc prior-only stability sweep for visualisation. | COMPLETED | 0:24 | Submitted 2026-04-26, INTR (v3 resolution: nlive=400, num_repeats=5, prec=0.01, ntasks=76). **log Z = -0.073 ± 0.007**, **D_KL = 0.024 nats** (vs v4's 1.025 — 40× collapse confirms the v4 amplify posterior was ghost contamination). Max A in chain = 1.079 — modest decoupling-corner amplification, well below the A ≈ 1.58 finding at t=5. Corner plot: `hpc_results/28474676/corner_amplify_v5.png`. |
| 28477675 | T1 Dark-Photon-Plasma | **suppress v5** — same fix, P_max:minimize | COMPLETED | 0:38 | Submitted 2026-04-26, INTR. **log Z = +0.66 ± 0.05**, **D_KL = 1.98 nats** — structurally informative (similar to v4's 1.95 nats; the suppress region IS robust to the stability guard because the decoupling corner gives both small P AND no growing modes). Max suppression P_GR/P ≈ 1e5 at MAP (mA2≈0.97, deltam≈-0.45, xi≈0.80, alpha3≈0.001). Corner plot: `hpc_results/28477675/corner_suppress_v5.png`. |
| 28519213 | T4 Ricci-EM | **amplify v5** (Padé-probe stability guard, post-#322 commit c10aa8a) | COMPLETED | 0:03:23 | INTR; ran cleanly but inference.json + results.csv vanished post-completion (cause unknown, only corner_amplify.png remained on /rds — possibly cleanup race). Log Z = −2.26, D_KL = 1.72 from the corner image. Resubmitted as 28520217. |
| 28520217 | T4 Ricci-EM | **amplify v5b** — same priors, files captured | COMPLETED | 0:03:30 | INTR. **log Z = −2.26 ± 0.07** (model strongly disfavoured vs null, Bayes factor ≈ 0.10 against), **joint D_KL = 1.79 ± 0.06 nats**, ESS = 1356. Marginal D_KL: δ₁ = **1.68 nats** (dominant); α₁₋₃ < 0.07 nats. A range [4.95×10⁻⁹, 1.255]. Top samples cluster at small \|δ₁\| ∈ [0.07, 0.15]; predicted A_coupling ≈ 1.9 (`docs/AMPLIFICATION_INVESTIGATION.md`) NOT confirmed — the inference finds max A ≈ 1.26.  Bottom of chain shows extreme suppression (~10⁻⁹ at \|δ₁\| ≈ 1.3), much stronger than the predicted 10⁻³ valley. Corner: `hpc_results/28520217/corner_amplify_v5.png`. |
| 28519675 | T4 Ricci-EM | **suppress v5** | PENDING | — | Standard QOS (parallel with INTR amplify), 3h wall, submitted 2026-04-27. |

---

## Key Findings

*(Filled in as results arrive)*

- Stage D1 v5 (T4 Ricci-EM, 2026-04-27): **δ₁ R̃ₘᵤᵥFᵘᵥ coupling drives strong destructive interference, modest amplification ceiling.**
  - **Stability-guard refactor required first** (#322, commit c10aa8a). Old eigenvalue+pinv path 100% rejected T4 due to high cond(V); new Padé matrix-exponential probe correctly handles IC-decoupled growing modes — same numerical machinery as the modal solver's Pass 0 path-D evolution.
  - **Amplify (28520217) verdict**: log Z = −2.26 ± 0.07 (Bayes factor ≈ 0.10 vs null, model strongly disfavoured for amplification), joint D_KL = 1.79 ± 0.06 nats. **Marginal D_KL is concentrated almost entirely on δ₁ (1.68 nats; α₁₋₃ each < 0.07 nats)** — confirming `R̃ₘᵤᵥFᵘᵥ` is the dominant coupling mechanism, with the PGT torsion-mass parameters (α₁₋₃) playing a negligible role within the explored range.
  - **A range [4.95×10⁻⁹, 1.255]** spans 8 decades. Top samples (A ≈ 1.26) cluster at small \|δ₁\| ∈ [0.07, 0.15] with α₂ ∈ [-1.9, -0.8].  Bottom samples (A ≈ 10⁻⁹) cluster at \|δ₁\| ∈ [1.2, 1.3], any α₂.
  - **Comparison to `docs/AMPLIFICATION_INVESTIGATION.md` predictions**:
    - Schur-complement A_coupling ≈ 1.9: **NOT confirmed**. Max A in chain is 1.26.  The predicted optimum at δ₁ ≈ 1, α₂ ≈ -1 was sampled (15% Padé-probe rejection rate at δ₁=1) but the inference did not select it — the actual simulation gives lower amplification than the linearised Schur estimate.
    - Suppression valley A ≈ 10⁻³: **MASSIVELY EXCEEDED**.  Inference finds A ≈ 10⁻⁹ in the deep valley — the destructive-interference mechanism (δ₁ R̃·F coupling shifting the photon's effective dispersion away from resonance) is much more effective than the analytic estimate.
    - Tachyonic boundary at α₂ ≈ -0.91: Padé probe rejects 18-51% of α₂ < -1 region (consistent), but 0% rejection at α₂ > 0 (where the original investigation expected divergences via Gaussian-IC contamination).  The plane-wave IC is much cleaner than the Gaussian IC used in `AMPLIFICATION_INVESTIGATION.md`'s heatmap.
  - **Physical interpretation**: T4's δ₁ coupling shifts the photon effective mass through the Schur complement.  In one direction (small δ₁) you get modest amplification through resonance enhancement (A up to 1.26).  In the other (\|δ₁\| ≈ 1) you get destructive interference suppressing conversion by 9 orders of magnitude.  The model is **not a Gertsenshtein amplifier** — its Bayes factor against the null is 1:10.  But it is a **strong nonminimal modifier** of the conversion (D_KL = 1.79 nats — substantial structure).
  - **Suppress (28519675) PENDING** on standard QOS — will confirm/refine the suppression-valley structure.
  - Corner plot: `hpc_results/28520217/corner_amplify_v5.png`.

- Stage A v5 (stability-guard fix, 2026-04-26): **Definitive null on amplify; informative suppress at decoupling corner.**
  - **Pre-flight stability guard wired into inference** (commits e361113, 3855fc1, 7ef182d, 9e6776e). `_evaluate_likelihood` now calls `check_conversion_stability(conservative=True)` before any simulation; samples with Re(λ) > 0.3 in the source-containing block return logL=-inf with `run_status='tachyonic'` and the maximum growth rate. Conservative path skips the IC-coupling filter when cond(V) > 1e12 (typical for CDT models) to prevent false-negatives. Fixed coupling_floor bug where cond ≥ 1e14 forced floor=1.0 → all growing modes silently skipped.
  - **Per-sample rejection metadata** flows through MC and nested-sampling pipelines as InferenceResult.metrics. PolyChord rejects -inf samples upfront (they never enter the chain), so for visualisation we run a post-hoc prior-only stability sweep (5000 prior draws × ~14 ms each ≈ 70 s, < 0.5% of inference wall time) saved as `_rejected_prior.csv`. Corner-plot tool overlays these rejected prior samples in red on the upper triangle, alongside posterior contours on the lower triangle — cosmology-paper convention (Planck/DES/getdist) of solid two-tone fills (light = 95% CR, dark = 68% CR), identical across all panels.
  - **Amplify (28474676) — null result confirmed.** log Z = -0.073 ± 0.007, joint D_KL = 0.024 nats (vs v4's 1.025 — 40× collapse), Bayes factor 0.93 vs GR baseline. Max A in chain = 1.079 (vs v4's 1454). 86.6% of accepted samples have logL > 0 (modest amplification above GR), but 0 samples reach the A ≈ 1.58 threshold reported at t=5 in v4. Top-10 cluster: large mA2, small xi, small α₃, deltam ≈ 0 — the decoupling corner. **Interpretation:** v4's "structurally informative" amplify posterior was entirely ghost contamination; once positive eigenvalues are excluded, no genuine kinematic amplification exists in the prior range. The previous t=5 A ≈ 1.58 result was either (a) ghost growth not yet visible at t=5, or (b) genuine modest amplification rejected by the Re(λ) > 0.3 threshold. Distinguished by a future t_end=5 INTR run.
  - **Suppress (28477675) — informative null on the decoupling axis.** log Z = +0.66 ± 0.05, joint D_KL = 1.98 nats — comparable to v4 suppress's 1.95 nats. Max suppression factor P_GR/P ≈ 1e5 at MAP (mA2≈0.97, deltam≈-0.45, xi≈0.80, α3≈0.001). The suppress posterior **is structurally informative** because the small-α₃ decoupling corner gives both (i) small P_max (favoured by minimize) and (ii) zero growing modes (passes the stability guard). The two-fold consistency of the v4 and v5 suppress D_KL values confirms decoupling is the genuine suppression mechanism, not a ghost artefact.
  - **Visualisation improvements** for these and future inference runs: append logL as final corner row/column (anesthetic auto-handles weight-aware KDE with TeX label r'$\ln\mathcal{L}$'); legend with 68%/95% CL contour patch + MAP+CI + prior + rejected-region entries; suptitle with log Z ± err, max A (mode-aware), D_KL nats; MAP solid line + 68% CI grey band on each diagonal; solid two-tone fills (Planck-blue palette) post-hoc applied via `_force_solid_credible_fills` to override anesthetic's per-panel-normalised gradient default. Documented in `tidal/inference/_visualize.py`. See `hpc_results/28474676/corner_amplify_v5.png` and `hpc_results/28477675/corner_suppress_v5.png`.

- Stage A v4 (Tier 1 + 1.5 path-D solver, 2026-04-26): **Wall-time vindication + tachyonic-instability discovery.** [SUPERSEDED BY v5: amplify posterior was ghost contamination.]
  - **Amplify (28418115) COMPLETED in 41 min** despite 5× higher resolution than v3 (which took ~3h). Confirms user's hypothesis: v3's runtime was dominated by numerical-instability churn that path-D + augmented-exp eliminated.
  - **log(Z) = +0.141 ± 0.021** (v3: +0.118 ± 0.006), **joint D_KL = 1.025 nats** (v3: 0.043 — 24× more posterior structure now visible). MAP at light-mediator corner (mA2=0.001, deltam=+0.28, xi=0.32, alpha3=0.003).
  - **Top samples reach P_max ≈ 3.6 at t=10, A ≈ 1454.** P_max > 1 means the linearised theory is being evolved past its validity window. Verified independently via `tidal simulate`: same value. Verified bit-exact across CDT ↔ FV formulations: identical P_max(t) and identical divergence at t≥15. Tier 1 + 1.5 path-D fix is correct; the apparent "amplification" is **a real prediction of the linearised model**, not a CDT-specific numerical artefact.
  - **Diagnostic verdicts (top sample):**
    - B₀ scaling test: A constant at 1454 ± 1% across B₀ ∈ [0.001, 0.03] → linear regime, equations are linear.
    - t_end test: A(5)=1.58, A(10)=1454, A(14)=3×10⁶, A≥15 → SimulationDivergedError. Per CLAUDE.md t-independence test, A(2t)/A(t) >> 1 = exponential growth.
    - Eigenvalue analysis: every k-mode has Re(λ) > 0 (γ scales linearly in k from 2.2 at k=0 to 94 at k=4). The h₅ component of unstable eigenvectors is ~10⁻²¹ (machine-noise direct projection), but the B₀ source term in the EOM couples h₅ → unstable modes at order B₀, giving observed effective growth rate γ_eff ≈ 1.84.
    - **Within strict linear regime (t≤5, P_max < 0.1):** valid amplification A ≈ 1.58 — modest, real.
  - **Interpretation:** the dark-photon-plasma model has positive linearised eigenvalues across the prior in the light-mediator corner. Within linearised theory this is exponential conversion; physically, unitarity bounds P ≤ 1 so back-reaction must saturate the conversion at some scale we haven't computed. The amplification IS real in linearised theory, but extracting an observable "amplification factor" requires either (a) restricting the prior to the perturbatively stable region (Re(λ) ≤ 0 with strong physical coupling), or (b) implementing the non-linear back-reaction.
  - **Suppress (28418421) TIMEOUT at 1h INTR.** Partial chains (2030 dead samples) analysed via anesthetic: log Z (partial) = +0.52, D_KL (partial) = 1.12 nats — same scale as amplify, NOT the v3 trivial decoupling concentration. Resubmitted as 28461922 at v3 resolution to validate workflow.
  - **Plots:** `hpc_results/28418115/diagnostics/Pt_top_sample.png` shows the P(t) trajectory across linear → transition → tachyonic regimes. `hpc_results/28418115/corner_amplify.png` is the 4D corner.

- Stage A (v3 — superseded): **Amplify NULL confirmed; Suppress rerun pending IC-snap fix.**
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
