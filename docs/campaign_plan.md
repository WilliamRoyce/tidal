# PGT Survey Campaign — Strategic Plan (stages C onward)

**Central question:** Can PGT torsion amplify the Gertsenshtein effect? If so, which
Lagrangian sector and at what parameters?

**Stage progress tracker:** [CAMPAIGN.md](../CAMPAIGN.md) (theory registry, HPC job log,
key findings, GitHub issues).

**Submit scripts (persistent):**
- [scripts/campaign/submit_campaign_AB.sh](../scripts/campaign/submit_campaign_AB.sh) — Stage A (T1 paired) + Stage B (T2 null)
- [scripts/campaign/submit_campaign_C.sh](../scripts/campaign/submit_campaign_C.sh) — Stage C (T3 b5 perturbative, both orders)
- [scripts/campaign/submit_campaign_D.sh](../scripts/campaign/submit_campaign_D.sh) — Stage D (T4 Ricci-EM / T5 YM-PGT / T6 YM-PGT-CP, paired)

**Branch:** `hpc/pgt-survey`. Fixed params throughout: `kappa=1.0`, `B0=0.01`.

---

## Theory Inventory

All sectors involving the EM field F (relevant for h↔a conversion), quadratic in fields,
from `research/general_quadratic_lagrangian.tex`:

| Research sector | Count (after Bianchi) | Params (LaTeX) | In existing theory? |
|----------------|----------------------|----------------|---------------------|
| T² (torsion mass) | 3 | α₁,α₂,α₃ → beta1-3 in T5 | T2–T6 ✓ |
| F² (Maxwell kinetic) | 1 | — | T1–T6 ✓ |
| Ft² (torsion-trace kinetic) | 1 | ξ | T1,T5,T6 ✓ |
| F·Ft (kinetic mixing) | 1 | δₘ | T1 ✓ |
| R̃×F (nonmin Ricci-EM) | 1 | δ₁ | T4–T6 ✓ |
| R̃×Ft (Barker χ₁) | 1 of ≤10 | χ | T5,T6 ✓ (χ₁ only) |
| ∂T×F | 3 | ζ₁,ζ₂,ζ₃ | T5,T6 ✓ |
| R̃×∂T beyond χ₁ | ≤9 more | χ₂…χ₁₀ | **MISSING** from all |
| R×R (curv-squared) | 6 | α₁…α₆ (b5 = α₁) | b5 only in T3 |
| ε·R̃×F | 3 | d₁₉,d₂₀,d₂₁ | T6 ✓ |
| ε·∂T×F | 6 | zt₁–zt₆ | T6 ✓ |
| ε·R̃×∂T | ≤36 | χ̃₁–χ̃₃₆ | **MISSING** from all |
| ε·R×R | 13 | d₁–d₁₃ | **MISSING** from all |

All R×F×F and similar cubic-in-fields terms are deliberately excluded (campaign stays
quadratic in all field strengths). **T5 (general_nonminimal)** includes χ₁ but χ₂–χ₁₀ are
absent. **T6 (parity_odd)** adds ε·R̃×F and ε·∂T×F but lacks ε·R̃×∂T (up to 36 terms). The
R×R sector has 6 terms but only b5·R̃² = α₁ is derived.

**Production-ready theories** — details in [CAMPAIGN.md](../CAMPAIGN.md) Theory Registry
(T1-T6 with JSON paths, free params, IC fields h=h_5, a=a_1 for all). Stage E theories
(T7 Complete-Even, T8 Complete-Odd) need new Wolfram derivations: base on T5/T6
theory.toml, add the missing χ, α, χ̃ sectors listed above.

**Sub-theories of T5** — see [CAMPAIGN.md](../CAMPAIGN.md) Sub-theories table.
Barker-PGT (β₁₋₃, ξ, χ), Shapiro-PGT (β₁₋₃, ζ₁₋₃), Bahamonde-PGT (β₁₋₃, δ₁) all run as
Stage D2 sub-runs with unused params fixed to 0 — no new JSON needed.

---

## Campaign Stages (C / D / E)

**Stages 0, A, B already complete** — see [CAMPAIGN.md](../CAMPAIGN.md). Stage A at
t_end=10 is being rerun at t_end=50 to expose plasma-mass suppression of mA² (separate
plan file).

### Stage C — R²-PGT perturbative b5 (T3)

- Gate: Stage 0 local check at b5 ∈ {0, 1e-3, 5e-3} already found PHYSICALLY NULL — b5
  decouples from TT channel under h_5 IC, see issue #299; Stage C HPC runs anyway as
  consistency check — expected flat posterior.
- Submit: `bash scripts/campaign/submit_campaign_C.sh` (runs both perturbative orders 0
  and 1). Priors: `b5=log_uniform:1e-4:1e-2, alpha1,2=uniform:-1:1, alpha3=log_uniform:0.05:2`,
  nlive=400, apply **`--t-end 50.0`** (lesson from Stage A — edit the script before
  submit).
- Success criterion: D_KL(b5, ord=1) − D_KL(b5, ord=0) > 0.01 nats would indicate a real
  perturbative correction; expected value near 0 (confirms #299).

### Stage D — Nonminimal survey (T4 → T5 → T6)

- Submit: `bash scripts/campaign/submit_campaign_D.sh d1` (Ricci-EM, 4 params, paired
  amp+supp), then `d2` (YM-PGT 9 params + Barker/Shapiro sub-theories), then `d3`
  (YM-PGT-CP 22 params — d18 genuinely absent, it's the EM θ-term total derivative).
- Apply `--t-end 50.0` to all (edit scripts before submit — they currently have
  `--t-end 10.0`).
- D2 sub-theories submitted in parallel with full D2: full YM-PGT (amp+supp),
  Barker-PGT amp, Shapiro-PGT amp. All use same JSON/ICs.
- D3 is the 22-param run; nlive ≥ 500; standard QOS 3h wall likely needed.
- Submit paired (amplify + suppress) for D1 and D3. D2 includes sub-theory variants.
- Success: any D_KL > 0.1 nats → real signal in that sector. Run T8 derivation ONLY if
  D3 shows parity-odd signal.

### Stage E — Wolfram derivations for T7 / T8

- **T7 Complete-Even-PGT (high priority)**: base on T5 theory.toml + add χ₂–χ₁₀
  (8 terms from `research/general_quadratic_lagrangian.tex` §R̃×∂T) + α₂–α₆ from §R×R
  (4 independent in 4D due to Gauss-Bonnet). Design `examples/pgt_complete_even/theory.toml`;
  put all α, χ in `[perturbation]`. Derivation: `tidal derive ... --timeout 0`
  (~35–45 min). Before HPC: smoke test, confirm ≥1 new χ or α correction non-zero.
- **T8 Complete-Odd-PGT (conditional)**: base on T6 + add ε·R̃×∂T (≤36 χ̃ terms) +
  ε·R×R (13 d terms). Run T8 derivation ONLY IF Stage D3 shows parity-odd signal.
  Literature treats all ε sectors as one family — all-at-once, not piecemeal.
  Derivation may exceed 1h; profile on smoke test first.

---

## Simulation Conditions (apply uniformly)

- `--ic plane-wave` (no Gaussian, no noise)
- `--param B0=0.01` (small-P regime)
- `--ic-amplitude 1e-2` (≪ B0, valid linearization)
- `--periodic` `--bounds 0:50`
- `--grid-shape 64` (HPC) / 32 (smoke)
- **`--t-end 50.0`** for any theory with a mass scale (mA², ξ, δ₁, ζ*, χ*); `--t-end 10.0`
  only acceptable for pure-torsion nulls like T2 where no mass scale exists. Modal
  solver cost is independent of t_end (analytic `exp(λ·t)` evaluation at snapshot times)
  — rule of thumb: t_end ≳ 5× the longest relevant oscillation period. At t_end=50 the
  stability guard cutoff tightens: modes with tiny positive Re(λ) that pass at t_end=10
  may be rejected at t_end=50 — physically correct, those ARE unstable on relevant
  timescales.

---

## Likelihood Strategy — paired amplify + suppress

- **Amplify run:** `--likelihood 'P_max:maximize'` — posterior ∝ P_max, finds amplification
  peaks (P > P_GR).
- **Suppress run:** `--likelihood 'P_max:minimize'` — `logL = -P_max`, posterior
  concentrated where P_max → 0. Handley's principled choice; no σ tuning.
- Baseline for both: `--baseline-formula 'sin(kappa*B0*t_end/2)**2'`.
- **Exception**: Stage B (null) ran amp only.
- **Clustering**: remove `--no-clustering` for Stage D/E runs (multi-modal posteriors
  possible). Keep `--no-clustering` on Stage C and Stage B.

---

## QOS heuristic

- Fast/known runtime (suppression, low-D nulls) → `polychord_intr.sbatch`, 1h wall.
  INTR has `MaxSubmitPU=1` — only one INTR job in queue at a time.
- Slow/unknown (amplification, high-D, Stage D2/D3) → `polychord_standard.sbatch`, 3h wall.

---

## Success Criteria

| Stage | Metric | Threshold | Significance |
|-------|--------|-----------|-------------|
| A | max D_KL(params) | > 0.05 nats | Dark photon resonance found |
| B | max D_KL(params) | < 0.005 nats | Minimal PGT null confirmed |
| C | D_KL(b5, ord=1) − D_KL(b5, ord=0) | > 0.01 nats | b5 correction measurable |
| D1 | D_KL(delta1) | > 0.05 nats | Nonminimal Ricci-EM effect measured |
| D2 | any D_KL | > 0.1 nats | Parity-even survey: first signal |
| D3 | any D_KL | > 0.1 nats | Parity-odd survey complete |
| E1 | D_KL(any χ₂+ or α₂+) | > 0.1 nats | Combined χ+R×R measurable |

---

## Non-Propagating / Constraint Torsion Theories

Deferred (issue #297). Non-propagating sectors easy to add: fix kinetic params to 0 in
existing JSONs. Trigger: Stage B posterior concentrates near single-sector corner (didn't
happen at t_end=10; re-evaluate after Stage D).

---

## Euler-Heisenberg Extension

Deferred (issue #296). #271 xAct F⁴ blocker resolved. 2D PolyChord run:
L = R̃/κ² − ¼F² + c₁(FF)² + c₂(FF̃)². Tests QED-vacuum-polarization enhancement of
Gertsenshtein independent of torsion; relevant for IAXO/CAST comparisons.

---

## Out of Scope (for entire campaign)

- R×F×F, T×T×F, R×R×F etc. (cubic in field strengths)
- Multi-node sweeps (`tidal sweep` is single-node multiprocessing.Pool only)
- 3+1D theories (plane-wave reduction sufficient)

---

## Critical Files

- JSON theories: `examples/data/{dark_photon_plasma, torsion_gertsenshtein_b5_zero,
  torsion_gertsenshtein, torsion_gertsenshtein_nonminimal,
  torsion_gertsenshtein_general_nonminimal, torsion_gertsenshtein_parity_odd}.json`
- Stage E bases: `examples/torsion_gertsenshtein_general_nonminimal/theory.toml` (T7 base),
  `examples/torsion_gertsenshtein_parity_odd/theory.toml` (T8 base)
- Sector definitions: `research/general_quadratic_lagrangian.tex`,
  `research/enumeration_physical.json`
- Submit scripts: `scripts/campaign/submit_campaign_{AB,C,D}.sh`
- HPC templates: `scripts/hpc_templates/polychord_{intr,standard}.sbatch`
- HPC shuttle: `scripts/hpc_shuttle.sh`

---

## Memory / key references

- [project_pgt_campaign_stageA_B_results.md](https://) (Claude memory) — Stage A+B
  results summary
- Issue #299 — b5 decouples from TT channel (Stage C expected null)
- Issue #300 — TT plane-wave IC is not WLOG (future: run with torsion IC too)
- Issue #307 — P_max > 2.0 likelihood cap (defense against stability-guard misses)
- Issue #308 — marginal D_KL bug (fixed)
