# v3 Phase Tracker

**Created:** 2026-05-10
**Companion to:** [V3_ARCHITECTURE.md](V3_ARCHITECTURE.md)
**Purpose:** running checklist of v3 phases — updated as each task converges (mirrors how Phase 6.J/K/L/M were tracked in CAMPAIGN.md).

## Phase A — Soft-penalty refactor + compactified-prior + de-pruning (in progress)

| Step | Description | Status | Issue / Notes |
| --- | --- | --- | --- |
| A.0a | `docs/V3_ARCHITECTURE.md` written | ✅ done | 2026-05-10 |
| A.0b | `docs/V3_PHASE_TRACKER.md` written | ✅ done | (this file) |
| A.0c | `docs/V3_PHASE_E_DESIGN.md` written | ✅ done | Deferred-but-documented |
| A.0d | `CAMPAIGN.md` updated with v3 architecture pointer | ✅ done | |
| A.0e | GitHub issues created for v3-A/A-γ/B/D/E/C tracking | ✅ done | #345–#355 |
| A.0f | Obsolete HPC jobs (28982006, 28985879) cancelled | ✅ done | 2026-05-10 |
| A.0g | Phase A.0 persistence committed | ✅ done | 62c7ac9 |
| A.1 | Soft-penalty refactor in `tidal/inference/_likelihood.py` | ✅ done | [#345](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/345) — Hwang-Noh + P_max>2 cap removed; soft floor `-100 + Normal(0, σ)` for sim/NaN/exception; distinct `run_status` tags. New `--gated` and `--soft-floor-noise SIGMA` flags. |
| A.2 | Lagrangian de-pruning audit + per-param `arctan_uniform` scripts | ✅ done | [#346](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/346) — `docs/lagrangian_depruning_audit.md` + 12 v3_permissive scripts |
| A.3 | Corner-plot upper-triangle removal in `tidal/inference/_visualize.py` | ✅ done | [#347](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/347) — `_hide_upper_triangle()` helper |
| A.4 | Soft-penalty tests (`test_likelihood_*.py`) | ✅ done | [#348](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/348) — 3 new modules, 36 new test cases, all green |
| A.5 | D1 v1 chain replay sanity check | ✅ done | [#349](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/349) — `examples/data/v3_d1_replay/`; v3 admits 14.3% more samples than v2 |
| A.6 | Commit + version bump | ✅ done | 400a455 (v0.39.0) |

## Phase A-γ — γ_conversion (deferred to its own session)

| Step | Description | Status |
| --- | --- | --- |
| A-γ.1 | Refactored γ_conversion with multi-t_test sampling + log-zero clamp in `_stability.py` | ⏳ deferred |
| A-γ.2 | `gamma_conversion` metric in `parse_likelihood()` (opt-in) | ⏳ deferred |
| A-γ.3 | 10-point D1 P_max correlation validation | ⏳ deferred |
| A-γ.4 | Decision: ship γ_conversion for amp campaigns or keep as metadata | ⏳ deferred |

GH issues `[v3-A-γ] γ_conversion: refactored definition + multi-t_test + log-zero clamp` and `[v3-A-γ] γ_conversion: P_max correlation validation` track the work.

## Complete theory inventory and campaign order

> **Breadth-first policy (2026-05-22):** ALL theories with existing JSONs must have landscape surveys before ANY pub quality runs. Pub quality requires explicit user approval after all landscapes complete.

| Priority | Theory / label | ndim | Landscape nlive | Pub nlive | num_repeats | grid | JSON | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ✅ 1 | D1 Ricci-EM | 4 | 100 | 200 | 8 | 512 | torsion_gertsenshtein_nonminimal.json | ✅ B.2/B.3 (nlive=1200) |
| ✅ 2 | Stage A DPP | 4 | 100 | 200 | 8 | 64 | dark_photon_plasma.json | ✅ B.4-full (nlive=1200) |
| ✅ 3 | D2.0 Bahamonde | 5 | 125 | 250 | 10 | 64 | torsion_gertsenshtein_general_nonminimal.json | ✅ B.5 landscape |
| ✅ 4 | D2.1 Barker | 6 | 150 | 300 | 12 | 64 | same | ✅ B.5 landscape |
| ✅ 5 | D2.2 Shapiro | 8 | 200 | 400 | 16 | 64 | same | ✅ B.5 landscape |
| ✅ 6 | D2.3 Full T5 | 9 | 225 | 450 | 18 | 64 | same | ✅ B.5 landscape |
| ✅ 7 | T6-minimal parity-odd | 4 | 100 | 200 | 8 | 64 | torsion_gertsenshtein_parity_odd_minimal.json | ✅ B.5-t6 landscape (29515407) |
| ✅ 8 | T6 parity-odd YM-PGT | 20 | 500 | 1000 | 40 | 64 | torsion_gertsenshtein_parity_odd.json | v1 (29567416) overflow-contaminated. ❌ v2 rescue (29687506) also floor-dominated: **amp logZ=−15.01±0.023, sup logZ=−14.99±0.025**. Root cause: parity-odd sector (d14-21, zt1-6) generates Re(λ)~1406 at t_end=1, still above γ>709 overflow threshold. Corner plots in `hpc_results/29687506/`. **Phase E required.** |
| ✅ 9 | T7 Complete-Even PGT | 18 | 450 | 900 | 36 | 64 | torsion_gertsenshtein_complete_even.json | v1 (29588680) overflow-contaminated. ✅ **v2 rescue (29682868): amp logZ=+10.30±0.077, sup logZ=−3.18±0.067** (still-active at timeout). Genuine amplification signal. Corner plots: `hpc_results/29682868/corner_t7_v2_{amp,sup}.png`. |
| ❌ 10 | T8 Complete-Odd PGT | 53 | 1325 | 2650 | 106 | 64 | torsion_gertsenshtein_complete_odd.json | ❌ B.5-rescue failed (29684804, 1hr, **0 dead points**). xt sector eigenvalues γ≈5.5×10⁵ → expm overflows at ALL viable t_end. **Plane-wave geometry fundamentally incompatible. Phase E required.** |
| ✅ | T9 Complete-Even + full ξ kinetic | 32 | 800 | 1600 | 64 | 64 | torsion_gertsenshtein_complete_even_full_xi.json | v1 (29596287) overflow-contaminated. 🟡 **v2 rescue (29694142): amp logZ=−0.18±0.072, sup logZ=−9.61±0.085** (still-active at timeout). Both negative/null — partially floor-contaminated (ξ sector overflow persists at t_end=1 for fraction of prior); T7 is a sub-model so amp should be positive in a clean run. Corner plots: `hpc_results/29694142/corner_t9_v2_{amp,sup}.png`. |
| ✅ | EH+Gert (F⁴+Gertsenshtein) | 2 | 50 | 100 | 4 | 64 | gertsenshtein_eh.json | ✅ **B.5-rescue (29700083): amp logZ=+0.184±0.0008, sup logZ=−0.185±0.0004**. Near-null and symmetric — EH F⁴ corrections (ρ,σ) produce negligible effect on Gertsenshtein conversion at B₀=0.01. QED point (σ/ρ=7/16) not distinguishable from GR baseline. Corner plots: `hpc_results/29700083/corner_eh_gert_v1_{amp,sup}.png`. **Topological control (EH+Top):** `lambdatop * F∧F̃` (Pontryagin density) added as 3rd coupling in `gertsenshtein_eh_top.json` (derived 2026-05-25). Result: `lambdatop` **absent from EOM and Hamiltonian** — symbolic pipeline correctly identifies total derivative and drops it. No HPC run needed; derivation is the proof. |
| ✅ | NP (non-propagating torsion, ξ=0) | 8 | 200 | 400 | 16 | 64 | torsion_gertsenshtein_general_nonminimal.json | ✅ **B.5-rescue (29700462): amp logZ=+10.98±0.176, sup logZ=+2.67±0.139** (TIMEOUT 1hr, still-active; 2655/3994 dead points). Both positive — **torsion-EM coupling alone (without propagating torsion kinetic term) is sufficient for Gertsenshtein enhancement above GR.** Comparable amp logZ to T7 (+10.30). Corner plots: `hpc_results/29700462/corner_np_v1_{amp,sup}.png`. |
| 🔄 | NP.T7 (non-propagating complete-even, ξ=0) | 17 | 425 | 850 | 34 | 64 | torsion_gertsenshtein_complete_even.json | 🔄 **B.5-post: scheduled** (2026-05-25). Extends NP.T5 (29700462, ndim=8) by adding chi1-10 (R̃×∂T basis) — the 10 interaction terms present in T7 but absent from general_nonminimal. xi=0 pinned via `--param`. ndim=17. Direct isolation: does T7v2's amp logZ=+10.30 require torsion kinetics (xi), or does the chi coupling sector alone drive the enhancement? Script: `polychord_intr_np_ceven_landscape.sbatch`. |

**T6 parity-odd params (20):** beta1-3, xi (log_uniform:1e-3:1e3), delta1, chi, zeta1-3 + d14, d15, d17, d19-21, zt1-3, zt5-6. (d16, zt4 vanish in plane-wave 1+1D reduction.)
**T6-minimal params (4):** beta1, beta2, beta3, d21.
**T7 design:** T5 base + chi1–chi10 (full R̃×∂T sector, 10 terms). No R×R — curvature-squared blocked by pipeline. ndim=18 (beta1-3, xi, delta1, zeta1-3, chi1-10).
**T8 design:** T6 base + xt1–xt36 (ε·R̃×∂T sector, 36 terms). No ε·R×R — blocked. `complete_odd.json` derived 2026-05-23 (38 fields, 48 H-terms). Deferred to Phase E.
**T9 design (derived 2026-05-23):** T7 base + full ξ₁–ξ₁₆ (∇T)² kinetic sector (16 terms, replaces single Barker ξ). Completes all ghost-free parity-even sectors. ndim=32 (xi11 vanishes in plane-wave reduction; 15/16 ξ terms survive). Source: `research/lagrangian_enumeration/explicit_terms_raw.txt` DT×DT block. Barker's ξ = LC{ξ₆,ξ₇,ξ₁₂,ξ₁₆} — all 16 independent kinetic invariants required for complete survey. Landscape sbatch: `scripts/hpc_templates/polychord_intr_t9_landscape.sbatch`.
**Kinetic gap note (2026-05-23):** T5/T6/T7/T8 all include only Barker's single ξ (1 of 16 independent (∇T)² invariants). T9 addresses this gap.
**EH+Gert design (2026-05-25):** GR+Maxwell Gertsenshtein base (graviton h_ij + photon A_μ, B₀ background coupling) + σ(F·F)² + ρ(F·F̃)² (Euler-Heisenberg F⁴ corrections). Both h_ij (graviton) and A_μ (photon) fields present — B₀ coupling drives graviton-photon conversion. ndim=2 (ρ, σ). QED point: σ/ρ = 7/16 (Dunne/Adler). Pure-EM `euler_heisenberg.json` has NO graviton fields and cannot test Gertsenshtein conversion — new derivation from `examples/gertsenshtein/theory_eh.toml` required.
**NP.T5 design (2026-05-25):** D2.3 (Full T5, general_nonminimal.json) with xi=0 pinned. Free: beta1-3, delta1, chi, zeta1-3 (8 params). Algebraic/Einstein-Cartan-like torsion; no torsion propagation. Result: amp logZ=+10.98 — proves torsion-EM coupling alone drives Gertsenshtein enhancement.
**NP.T7 design (2026-05-25):** T7 (complete_even.json) with xi=0 pinned. Free: beta1-3, delta1, zeta1-3, chi1-10 (17 params). Adds the chi1-10 (R̃×∂T basis) sector absent from NP.T5. Direct test: if NP.T7 amp ≈ T7v2 amp (+10.30), the chi coupling structure drives T7's enhancement without needing torsion propagation. If NP.T7 ≪ T7v2, xi (kinetic) is the necessary ingredient. Script: `polychord_intr_np_ceven_landscape.sbatch`.

## Phase B — Campaign re-runs

> **Landscape-first policy (2026-05-22):** All 6 theories must have landscape (25×ndim) chains before any publication-quality (50×ndim) work. D2.x landscape ✅ done. D1 and Stage A landscape ⏳ in progress (INTR 29514476). Pub quality requires explicit user approval.

| Step | Description | Status | Notes |
| --- | --- | --- | --- |
| B.0a | Joint-prior smoke comparison (D1 amp, 1 tile, INTR) | ⚠️ done — **does not pass adoption criterion** | 29204991 (8 min INTR). log Z=+13.76±0.09 (matches B.1 +13.29), ESS=1398 (matches B.1 1410), but **MAP α₁=+582, α₂=−575, δ₁=+575** vs B.1's ±0.5 — joint prior explores a different region than per-param. Per-coupling D_KL much lower (0.27–0.34 vs 2.96–3.20). Joint prior's r∈[1e-3,1e3] × sphere volume scaling concentrates at r_hi=1000; needs r_hi tuning. **Decision: keep per-param arctan for B.5.** See `docs/comparison/d1_amp_joint_v2_v3.md` |
| B.1 | D1 amp v3 smoke (mid-res INTR) | ✅ done | 29149987 (23 min INTR, grid=128/nlive=300). log Z=+13.29±0.13 (vs v2 +0.72; +12.5 nats), ESS=1410/3651, MAP δ₁=4.80 (outside v2 prior). Per-coupling D_KL: α₁=3.20, α₂=2.96, α₃=0.02, δ₁=1.44 nats. 100% success run_status. Post-hoc probe sweep flags 45% of prior as v2-tachyonic. See `docs/comparison/d1_amp_v2_v3.md` |
| B.2 | D1 amp v3 publication run | ✅ done | 29189748 (standard 8h, completed 2026-05-17, wall 1h43min). logZ=+13.374±0.060, ESS=5694, nlive=1200, grid=512. Consistent with B.1 smoke (+13.29); tight posterior confirms large prior shift vs v2 (+0.72). |
| B.3 | D1 sup v3 paired with B.2 | ✅ done | 29189761 (standard 12h, completed 2026-05-17, wall 6h37min). logZ=+11.395±0.097, ESS=5224, nlive=1200. High Bayes factor on suppression — suppression valley as deep as v2 but spread over wider prior volume. |
| B.4a | Stage A amp INTR smoke | ✅ done | 29189966 (8 min INTR). log Z=+9.31±0.13 (vs v2 −0.07; +9.4 nats), ESS=2649/6464. **v2 "null verdict" was an artefact of narrow priors** — v3 reveals joint D_KL=6.59 nats. MAP shifts ~1.1–2.0σ on all 4 params. See `docs/comparison/stage_a_amp_v2_v3.md` |
| B.4b | Stage A sup INTR smoke | ✅ done | 29199129 (22 min INTR). logZ=+4.08±0.12 (vs v2 +0.66; +3.4 nats). ESS=4961. Joint D_KL=7.30 nats. MAP mA2=605 (v2: 0.97 — massive prior shift). Per-param D_KL: deltam=3.02, xi=1.05, mA2=0.86, alpha3=0.63 nats. |
| B.4-full-amp | Stage A amp v3 publication | ✅ done | 29205968 (standard 6h, completed 2026-05-17, wall 1h37min). logZ=+9.034±0.068, ESS=7443, nlive=1200, grid=256. Consistent with B.4a smoke (+9.31). |
| B.4-full-sup | Stage A sup v3 publication | ✅ done | 29205982 (standard 12h, completed 2026-05-17, wall 2h47min). logZ=+2.924±0.083, ESS=10675, nlive=1200. Consistent with B.4b smoke (+4.08 — small shift from enlarged prior volume at pub nlive). |
| B.5-stageA | Stage A (DPP) landscape pass | ✅ redundant — superseded by B.4-full (2026-05-22) | **Note:** B.4-full (29205968/29205982) at nlive=1200 (300×ndim) already provides landscape and pub quality for Stage A. The 29514476 re-run at nlive=100 was submitted before recognising B.4-full was complete. Results: amp logZ=+9.641±0.223, ESS=647; sup logZ=+4.109±0.213, ESS=1369. Consistent with B.4a/B.4b smokes. **Gold standard: use B.4-full-amp/sup (nlive=1200) for publication.** |
| B.5-d1 | D1 (Ricci-EM) landscape pass | ✅ redundant — superseded by B.2/B.3 (2026-05-22) | **Note:** B.2 (29189748) and B.3 (29189761) at nlive=1200 (300×ndim) already provide landscape and pub quality for D1. The 29514476 Phase 2 D1 re-run was CANCELLED before completing. **Gold standard: use B.2/B.3 (nlive=1200) for publication.** |
| B.5 | D2.0–D2.3 v3 paired re-runs (8 chains) | ✅ done (landscape pass, 2026-05-17) | **All 8 chains captured at landscape nlive=25×ndim** — replaces broken nlive=1500–2400 runs that suffered from divergence-guard-era over-resolution. **Single-theory INTR**: D2.0 amp 29229768 (logZ=+9.29±0.06), D2.1 amp 29256858 (logZ=+9.58±0.05), D2.0 sup **29468539** (9 min, logZ=+5.382±0.200, ESS=1302). **Interactive INTR batch 29468763** (5/5 sequential @ 76 ranks, ~58 min): D2.1 sup (8m, logZ=+2.415±0.170), D2.2 amp (10m, logZ=+8.889±0.135), D2.2 sup (22m, logZ=+6.212±0.176), D2.3 amp (10m, logZ=+10.034±0.142), D2.3 sup SIGTERM at walltime. **Interactive INTR batch 29471255** (Phase 1: D2.3 sup at 76 ranks 26m, logZ=+2.662±0.156, ESS=2870; Phase 2: parallel mpirun test — D2.0 sup pub + D2.1 sup pub at 38 ranks each completed in 19m, **1.7× speedup over sequential**). |
| B.5-t6 | T6-minimal + T6 parity-odd landscape pass | ✅ done (2026-05-23) | **T6 full fresh: 29567416** (1hr INTR, TIMEOUT). logZ=-15.07±0.022 (amp), -15.17±0.023 (sup). **100% floor hits** — all 1040/1033 dead points in floor noise band [-17.8,-13.9]. Configurations with d19/d20/d21≠0 cause exponential field growth → SimulationDivergedError → soft floor (this is NOT a prior filter — Phase A removed the tachyon filter; the divergence is numerical). **Interpretation (revised 2026-05-23): floor-dominated = large amplification candidates. Phase E localised geometry required** to evaluate these with bounded transit-time growth. T6-minimal: 29515407/t6_minimal_amp/sup (logZ see stats; corner plots generated). |
| B.5-t7 | T7 Complete-Even PGT landscape | ✅ done (2026-05-23) | **29588680** (1hr INTR, TIMEOUT). logZ=-14.93±0.028 (amp), -14.92±0.028 (sup). **74% floor hits** (26% non-floor, max logL=-13.3). chi1-10 sector is less universally divergent than T6's ε·Riem×F. **Interpretation (revised 2026-05-23): floor-hit configurations are amplification candidates** — exponential growth means large potential signal, not null physics. Phase E required to properly evaluate. ndim=18, nlive=450. Full 18-param corner plots: `hpc_results/29588680/*/corner_full.png`. |
| B.5-t8 | T8 Complete-Odd PGT landscape | ❌ not viable in plane-wave geometry | T8 JSON derived 2026-05-23 (T6+xt1-36, 38 fields, 48 H-terms). Attempted B.5-rescue (29684804, t_end=1, 1hr INTR): **0 dead points**. xt sector (ε·R̃×∂T) generates eigenvalues γ~5.5×10⁵ → expm overflows float64 at t_end=1 (threshold γ>709 insufficient). t_end<0.001 would be needed but yields P≈10⁻¹⁷ (below soft floor). Plane-wave geometry fundamentally incompatible with T8. **Phase E (localised dual-Gaussian) required.** |
| B.5-rescue | Uniform-field rescue at t_end=1 (T6v2, T7v2, T8v1, T9v2, EH+Gert, NP) | 🟡 mostly complete (NP pending) | t_end=10 v1 runs identified as overflow-contaminated: modal solver's expm overflows float64 when max_eigenvalue × t_end > 708. At t_end=10 this triggers at γ>71 → 74–100% floor rates for T6/T7/T9. Modal solver expm cost is O(1) in t_end (empirically confirmed: ~0.1ms across t_end ∈ [0.1,10]). Rescue at t_end=1 safe to γ>709 with same machine-precision A=P/P_GR quality. **Also adding NP (non-propagating, xi=0) as new uniform-field run** — reuses `torsion_gertsenshtein_general_nonminimal.json` with `--param xi=0.0` pinned; free params β₁₋₃, δ₁, χ, ζ₁₋₃ (ndim=8). Script: `polychord_intr_np_landscape.sbatch`. |
| B.5-pub | D2.0–D2.3 v3 publication pass (50×ndim) | ⏸ PAUSED — awaiting T6/T7/T8 landscapes + explicit user approval | **Policy change (2026-05-22):** ALL landscape (25×ndim) surveys must complete before ANY pub quality runs. Pub runs require explicit user approval to resume. Partial results captured: D2.0 sup/amp pub (✅ done — 29471255, 29507332); D2.1 sup/amp pub (✅/🟡 — 29471255, partial 63% at 29507332); D2.2/D2.3 amp pub (🟡 40% at 29511699, tidal.resume exists — awaiting approval). Remaining: D2.1 amp pub resume, D2.2/D2.3 amp pub complete, D2.2/D2.3 sup pub. **DO NOT restart until user approves.** |
| B.6 | Bounds-dependence cross-check at L=75 under v3 | ✅ done | 29205638 (25 min INTR). log Z=+13.18±0.12 (nearly identical to B.1's +13.29). Per-coupling D_KL essentially unchanged. **v3 washes out v2's bounds-dependence non-monotonicity** (v2 L=75 = +0.84, v2 L=100 = +0.68 — non-mono; v3 L=75 ≈ L=100). |

## Phase B.5-post — Forward plan (after rescue campaign)

| Idea | Description | Gate |
| --- | --- | --- |
| EH+Gert atlas pilot | EH+Gert (ndim=2) is an ideal Phase C atlas proof-of-concept: 6 faces × 4 ranks = 24 ranks, fits one INTR slot; all faces in parallel (~10 min). Maps QED line σ/ρ=7/16 across atlas faces. | EH+Gert landscape gives valid (non-floor) logZ |
| NP.T7 landscape | Non-propagating T7: complete_even.json with xi=0 pinned, chi1-10 free (ndim=17). Isolates whether T7's amp logZ=+10.30 requires kinetic ξ. Direct NP.T5 (8D, +10.98) vs NP.T7 (17D) vs T7v2 (18D, +10.30) comparison. Script: `polychord_intr_np_ceven_landscape.sbatch`. Est. ~25-35 min INTR. | NP result (29700462) positive — follow-on motivated |

## Phase C — Cubed-sphere coupling-space chart (HANDLED BY PARALLEL SESSION)

See `docs/V3_PHASE_C_REFERENCE.md` once written by the parallel session. This phase is shipped externally to the v3 plan owned here.

Open follow-up question for supervisor: coupling grouping (monolithic / per-Lagrangian-symmetry-class / per-SPO-sector / per-parity / other). GH issue `[v3-C] Coupling-grouping question for supervisor follow-up` tracks.

## Phase D — Manuscript implications

| Step | Description | Status |
| --- | --- | --- |
| D.1 | `manuscript/sections/computational_approach.tex` — add §"Tachyon-permissive inference architecture" | ⏳ pending (after Phase B) |
| D.2 | Update results sections with v3 numbers; comparison table v2 vs v3 | ⏳ pending (after Phase B) |
| D.3 | `docs/PHASE_6_COMPARISON.md` — append "v3 architecture" section | ⏳ pending (after Phase B) |

## Phase E — Localised geometry pivot (DEFERRED, gated on Phase B convergence)

See [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md). Replaces plane-wave + uniform B₀ with Gaussian wavepacket + dual-Gaussian B-field profile; finite interaction time bounds P_max physically.

| Step | Description | Status | Notes |
| --- | --- | --- | --- |
| E.0 | Dual-Gaussian localised B-field empirical validation | ✅ done (2026-05-16) | `examples/gertsenshtein/theory_e0_dual_gaussian.toml` derived. sigB=5, zc1=25, zc2=75, L=100. A(20)/A(10)=1.00 < 1.05 ✓. P_max=0.00348 at t=18 (CVODE). Modal auto-selects but diverges → GH #367; use `--scheme cvode`. |
| E.1 | Wavepacket IC switch for all v3 campaign scripts | ⏳ deferred | After Phase B convergence |
| E.2 | Localised B-field theory TOMLs for all D1/D2/StageA models | ⏳ deferred | After Phase B convergence |
| E.3 | σ_w × σ_B tuning scan | ⏳ deferred | After E.2 |
| E.4 | Re-baseline campaign scripts in `v3e_localised/` | ⏳ deferred | After E.3 |
| E.5 | Re-run highest-value chains under E geometry | ⏳ deferred | D1 amp first |
| E.6 | Phase B vs Phase E comparison + publication decision | ⏳ deferred | After E.5 |

### Phase E theory roster

Planned localised-geometry runs (see `docs/PHASE_E_TRACKER.md` for full staging). All reuse the derived E-geometry JSON `torsion_gertsenshtein_general_nonminimal_e_dual_gaussian.json` (E.T5.3, 38 fields, 82 H-terms) unless noted.

| Theory | ndim | Description | Status | Physics motivation |
| --- | --- | --- | --- | --- |
| E.D1 | 4 | D1 Ricci-EM (non-minimal δ₁) | ⏳ Stage 1 | Highest Phase B logZ (+13.4 amp, +11.4 sup); Phase E priority 1 |
| E.StageA | 4 | Stage A DPP (dark-photon-plasma) | ⏳ Stage 1 | Phase B logZ +9.0/+2.9; control cross-check |
| E.T5 (E.D2.0) | 5 | Bahamonde sub-theory (β₁₋₃, ξ, δ₁) | ⏳ Stage 5 | Phase B landscape + pub chains complete; localised comparison |
| E.T5.1 (E.D2.1) | 6 | Barker sub-theory (+ χ) | ⏳ Stage 5 | Sub-theory with one extra coupling |
| E.T5.2 (E.D2.2) | 8 | Shapiro sub-theory (β₁₋₃, ξ, δ₁, ζ₁₋₃) | ⏳ Stage 5 | All torsion-EM ∇T×F couplings, with kinetic ξ |
| E.T5.3 (E.D2.3) | 9 | Full T5 (β₁₋₃, ξ, δ₁, ζ₁₋₃, χ) | ⏳ Stage 5 | Complete T5; primary localised survey target |
| **E.NP** | **8** | **Non-propagating torsion (ξ=0, all other D2.3 couplings free)** | **⏳ Stage 5** | **Key comparison: does torsion-EM coupling alone (without torsion propagation) produce Gertsenshtein enhancement? Reuses E.T5.3 JSON with `--param xi=0.0` pinned. Answers whether ξ is necessary for any amplification signal.** |
| **E.NP.T7** | **17** | **Non-propagating complete-even (T7 + ξ=0, chi1-10 free)** | **⏳ Stage 5** | **Richer NP theory: adds chi1-10 (R̃×∂T sector) absent from E.NP. Directly isolates whether T7's enhancement (+10.30) requires kinetic ξ or is driven by chi couplings alone. Uniform-field counterpart: NP.T7 (B.5-post, polychord_intr_np_ceven_landscape.sbatch). Reuses complete_even.json with `--param xi=0.0`.** |
| E.T6 | 20 | T6 parity-odd (B.5-rescue needed first) | ⏳ Stage 3 | Phase B floor-dominated; Phase E required for meaningful evaluation |
| E.T7 | 18 | T7 Complete-Even (B.5-rescue needed first) | ⏳ Stage 3 | T7v2 amp logZ=+10.30 — top candidate after T7 rescue |
| E.T8 | 53 | T8 Complete-Odd (plane-wave incompatible) | ⏳ Stage 6 | Only viable in Phase E geometry; xt sector eigenvalues overflow plane-wave |
| E.T9 | 32 | T9 Complete-Even + full ξ kinetic | ⏳ Stage 3 | B.5-rescue pending; full 16-invariant kinetic sector |
| E.EH+Gert | 2 | EH+Gert (F⁴ corrections, ρ/σ) | ⏳ post-B.5 | Atlas pilot candidate (2D → all 10 cubed-sphere faces feasible in 1 INTR slot) |

**E.NP physics rationale (2026-05-25):** D2.3 (Full T5) has ξ as its only propagating torsion DOF. Setting ξ=0 gives non-propagating / algebraic torsion (Einstein-Cartan-like): torsion is fully determined by local matter content, with no independent wave equation. All torsion-EM couplings (δ₁, ζ₁₋₃, χ) and torsion mass terms (β₁₋₃) remain free. If E.NP shows amplification, it proves the enhancement is driven by the torsion-EM coupling structure alone, not by torsion propagation. If E.NP is null but E.T5.3 is not, it isolates torsion propagation as the necessary ingredient. Script: `scripts/hpc_submit_drafts/v3e_localised/e_np_algebraic_amp.sh`.

## Sign-off log

(Append entries as phases converge — date, step ID, summary, commit hash.)

- 2026-05-10 — A.0a — `docs/V3_ARCHITECTURE.md` written; canonical architecture reference live.
- 2026-05-10 — A.0b — `docs/V3_PHASE_TRACKER.md` written (this file).
- 2026-05-10 — A.0c — `docs/V3_PHASE_E_DESIGN.md` written; Phase E deferred-but-documented.
- 2026-05-10 — A.0d — `CAMPAIGN.md` updated with v3 architecture pointer and HPC-job cancellations.
- 2026-05-10 — A.0e — GitHub issues #345–#355 created; `v3-architecture` label added.
- 2026-05-10 — A.0f — HPC jobs 28982006 + 28985879 cancelled.
- 2026-05-10 — A.0g — Phase A.0 persistence committed in 62c7ac9.
- 2026-05-10 — A.1 — Soft-penalty refactor in `_likelihood.py`: Hwang-Noh and `P_max>2` cap removed; soft floor `−100 + Normal(0, σ)` for sim/NaN/exception; distinct run_status tags. New `--gated` and `--soft-floor-noise SIGMA` flags wired through `_sample.py` + `cli/__init__.py`.
- 2026-05-10 — A.2 — Lagrangian de-pruning audit (`docs/lagrangian_depruning_audit.md`) committed: no de-pruning needed for Phase A's 6 campaigns; T6/EH applications deferred. Twelve `scripts/hpc_submit_drafts/v3_permissive/` campaign scripts written for D1, Stage A, D2.0–D2.3 paired chains.
- 2026-05-10 — A.3 — Corner-plot upper-triangle scatter rendering disabled via new `_hide_upper_triangle()` helper in `tidal/inference/_visualize.py`.
- 2026-05-10 — A.4 — Three new test modules `tests/test_likelihood_{soft_floor_noise,no_hwang_noh,permissive}.py` (36 cases); 107/107 inference tests pass; ruff/format clean.
- 2026-05-10 — A.5 — D1 v1/v2 chain replay landed at `examples/data/v3_d1_replay/`. Headline: v3 architecture admits 14.3% more samples than v2 (1428/10000 tachyonic samples now contribute via the v2→v3 admission shift on `28982029`). 28520217 v1 chain pre-dates the rejected-prior sidecar so the comparison uses 28982029 instead.
- 2026-05-12 — B.4b — Stage A sup v3 INTR smoke (29199129, 22 min): logZ=+4.08±0.12 (vs v2 +0.66; +3.4 nats), ESS=4961, joint D_KL=7.30 nats. MAP mA2=605 vs v2 MAP 0.97 — massive prior-shift. deltam D_KL=3.02 nats strongest. v3 unambiguously reveals structure masked by v2's narrow priors.
- 2026-05-12 — B.5 (D2.1 amp) — D2.1 Barker amp v3 INTR (29256858, 27 min): logZ=+9.58±0.05 (vs v2 +0.62; +8.96 nats), ESS=9517, joint D_KL=4.13 nats. β₂ MAP=−11.7 (outside v2's narrow [-3..-0.3] prior, 2.5σ shift). All 5 sign-symmetric params D_KL ≈ 2.5-2.8 nats. 27 clusters (multi-modal). See `docs/comparison/d21_barker_amp_v2_v3.md`.
- 2026-05-16 — B.5 (policy) — D2.0 sup truncated at session 7 (supervisor-approved): 58.5K dead pts, 174/401 active clusters, logZ=+134.47±0.453. `*_resume.sh` scripts deleted. 1-INTR-session-per-chain policy baked into `v3_permissive/README.md`. See `docs/comparison/d20_bahamonde_sup_v2_v3.md`.
- 2026-05-16 — E.0 — Dual-Gaussian localised B-field validated: `theory_e0_dual_gaussian.toml` derived (6 fields, 22 H terms). CVODE smoke tests pass; A(20)/A(10)=1.00 < 1.05 ✓; P_max=0.00348 at t=18. Modal solver auto-selects but diverges (GH #367 filed); workaround: `--scheme cvode`.
- 2026-05-17 — B.5 done — All D2.0–D2.3 amp+sup landscape pass captured at nlive=25×ndim (commit d786851). Validation run D2.0 sup 29468539 (9 min) confirmed default settings work. Interactive INTR batch 29468763 ran 4/5 sequential in ~58 min; batch 29471255 closed out D2.3 sup + pioneered 2-way parallel mpirun pattern (1.7× speedup over sequential, both rc=0). Bonus publication-quality reruns of D2.0 sup (nlive=250: logZ=+5.837±0.143, ESS=2404) and D2.1 sup (nlive=300: logZ=+2.524±0.158, ESS=2895) included in 29471255. Parallel mpirun pattern documented in `scripts/hpc_submit_drafts/v3_permissive/run_interactive_batch_2.sh` for future v3_publication and Phase E batches.
- 2026-05-22 — B.5-pub (D2.0 amp pub) — D2.0 Bahamonde amp pub complete: 29507332 (logZ=+7.679, ESS=1063). D2.1 amp pub SIGTERM'd at 63% (3451 dead pts, reconstructed logZ=+5.32±0.17 — resume required). D2.2+D2.3 amp pub launched: 29511699 (38+38 parallel, expect SIGTERM at ~75%; resume via run_pub_batch_d2amp_p2_resume.sh). Corner plots: hpc_results/29507332/{d20,d21}*/corner*.png.
- 2026-05-23 — B.5-t6 — T6 fresh landscape: 29567416 TIMEOUT 1hr, logZ=-15.07±0.022 (amp)/-15.17±0.023 (sup). 100% of 1040/1033 dead points at soft floor. Floor hits = SimulationDivergedError from exponential growth in infinite plane-wave (not prior filter — Phase A removed tachyon filter). **Interpretation: floor-dominated = large amplification candidates; Phase E localised geometry required to resolve.** T6-minimal corner plots generated (29515407). T7 re-derivation completed (~5 min, chi1-10).
- 2026-05-23 — B.5-t7 — T7 landscape: 29588680 TIMEOUT 1hr, logZ=-14.93±0.028 (amp)/-14.92±0.028 (sup). 74% floor hits (same mechanism as T6: numerical divergence in plane-wave). **Interpretation: floor-hit configurations are amplification candidates; chi1-10 sector has more stable configurations than T6 but still requires Phase E for meaningful survey.** Full 18-param corner plots generated. All B.5 landscape passes complete.
- 2026-05-23 — B.5-t9 — T9 Complete-Even + full ξ kinetic: TOML written, derived (~7m, ndim=32, xi11 vanishes), GR smoke ✓, landscape 29596287 TIMEOUT 1hr. logZ=-15.70±0.025 (amp), -15.68±0.026 (sup). ndead=818/800 nlive (nlike~90k). Floor-dominated — same mechanism as T6/T7: exponential growth in plane-wave geometry → Phase E required. All B.5 landscape passes now complete including T9. Closes kinetic sector gap.
- 2026-05-25 — B.5-rescue (diagnosis) — t_end=10 identified as root cause of T6/T7/T9 floor-dominated results: modal solver expm overflows float64 when max eigenvalue × t_end > 708. At t_end=10, theories with γ>71 overflow (74–100% of T6/T7/T9 parameter space). Modal expm cost is O(1) in t_end (empirically confirmed: ~0.09–0.20ms across t_end ∈ [0.1,10], noise-dominated). Rescue campaign planned at t_end=1 (safe to γ>709; A=P/P_GR t_end-independent at machine precision). T8 first attempt and EH+Gert (new derivation: GR+Maxwell+F⁴, ndim=2) included. V1 runs (29567416, 29588680, 29596287) remain archived as evidence of parameter-space instability structure.
- 2026-05-25 — B.5-rescue (T7v2, T8 failure, T6 running) — T7 v2 (29682868) confirmed fixed: amp logZ=+10.30±0.077 vs v1 −14.93 floor — rescue works for T7. T8 attempt (29684804) failed: **0 dead points in 1hr INTR**. Root cause: xt1-36 sector (ε·R̃×∂T) generates eigenvalues γ~5.5×10⁵ → expm still overflows at t_end=1 (threshold γ>709 insufficient by 3 orders of magnitude). Entire 75MB log = scipy RuntimeWarning overflow in matmul. t_end<0.001 required but P≈10⁻¹⁷ at that scale (below soft floor) — plane-wave geometry fundamentally incompatible with T8. T8 deferred to Phase E. T6 rescue (29687506) now running at t_end=1.
- 2026-05-25 — B.5-rescue (results, all theories) — T6 v2 (29687506): amp logZ=−15.01±0.023, sup=−14.99±0.025 — still floor-dominated; parity-odd sector (d14-21, zt1-6) generates Re(λ)~1406 → overflow persists at t_end=1. Parity-odd PGT incompatible with plane-wave geometry. T9 v2 (29694142): amp=−0.18±0.072, sup=−9.61±0.085 — partially rescued (moved from −15.7 floor) but still contaminated; both negative indicates ξ sector still causing partial overflow; T7 is a sub-model so amp should be positive in clean run. EH+Gert (29700083): amp=+0.184±0.0008, sup=−0.185±0.0004 — near-null symmetric result; EH F⁴ corrections negligible for Gertsenshtein at B₀=0.01. NP (29700462, ξ=0 pinned, 8D, general_nonminimal): amp=+10.98±0.176, sup=+2.67±0.139 — **strong positive signal; torsion-EM coupling alone (without propagating torsion kinetic term) is sufficient for Gertsenshtein enhancement.** Corner plots generated for all completed runs.
- 2026-05-25 — B.5-post (NP.T7 planned) — NP result (+10.98) motivates richer non-propagating investigation: NP.T7 uses complete_even.json (T7 basis) with xi=0 pinned. Adds chi1-10 (R̃×∂T sector, 10 terms) absent from NP.T5. ndim=17, nlive=425, nrep=34. Key question: does T7v2's amp logZ=+10.30 require kinetic ξ, or does the chi coupling structure alone drive the enhancement? Three-way comparison: NP.T5 (8D, +10.98) vs NP.T7 (17D, TBD) vs T7v2 (18D, +10.30). Script: `polychord_intr_np_ceven_landscape.sbatch`.
- 2026-05-22 — B.2/B.3/B.4-full — Pre-maintenance standard queue jobs confirmed COMPLETED (sacct check post-maintenance). D1 amp pub 29189748 (logZ=+13.374±0.060, ESS=5694, wall 1h43min); D1 sup 29189761 (logZ=+11.395±0.097, ESS=5224, wall 6h37min); Stage A amp pub 29205968 (logZ=+9.034±0.068, ESS=7443, wall 1h37min); Stage A sup pub 29205982 (logZ=+2.924±0.083, ESS=10675, wall 2h47min). All nlive=1200 (300×ndim — publication quality and above). D2.x landscape corner plots + D1/Stage A pub corner plots generated for all chains. Phase B.5-pub (D2.x pub pass) added as forward-plan step.
