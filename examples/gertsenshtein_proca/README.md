# Plasma-Gertsenshtein Effect (Massive Photon Extension)

Graviton-photon conversion in a background magnetic field with an effective
photon mass — the astrophysically-relevant regime where the photon acquires
a plasma mass from the surrounding medium.

**Physics summary.** In vacuum, the Gertsenshtein effect (1962) converts
gravitational waves into photons in a magnetic field with probability
`P = sin²(κB₀D/2)`. In a realistic plasma environment — magnetar magnetospheres,
primordial epochs, intracluster media — the photon acquires an effective mass
`mA` from the plasma, and the conversion transitions into the
Raffelt-Stodolsky (1988) two-state mixing regime with a Lorentzian resonance
in the detuning `Δm² = mA² - κ²B₀²/2`:

```
P_max = (2·κ·B₀·k)² / [(2·κ·B₀·k)² + Δm²²]
```

## Quickstart

```bash
bash reproduce_figures.sh              # default: reuse cached sweep data
bash reproduce_figures.sh --fresh      # force rerun all sweeps (~25 min)
```

Produces three figures in `examples/data/gertsenshtein_proca_figures/`:

- **fig1_lorentzian_1d.png** — 1D Raffelt-Stodolsky Lorentzian at B₀=0.10
  with analytical overlay (40-point scan across mA² ∈ [0, 1.0])
- **fig2_resonance_map_2d.png** — 2D heatmap P_max(B₀, mA²) with predicted
  resonance line (300-point grid across B₀ × mA²)
- **fig3_mA2_family.png** — P_max comparison at 6 discrete mA² values
  spanning on- to off-resonance at fixed B₀=0.10

The script is **HPC-friendly**: sweeps may be outsourced to a cluster and
the output directories copied back to a local checkout. `reproduce_figures.sh`
without flags detects existing data and goes straight to plotting.

## Files

| File | Purpose | Formula / Reference |
|------|---------|---------------------|
| `theory.toml` | Einstein-Maxwell + perturbation-level Proca Lagrangian | `L = (1/κ²)R - (1/4)F² - (mA²/2)a·a` |
| `run.sh` | Single baseline run (smoke test) | — |
| `sweep_B0_family.sh` | B₀ sweep at multiple mA² values | Raffelt-Stodolsky suppressed oscillation |
| `sweep_detuning.sh` | Detuning scan | Raffelt-Stodolsky Lorentzian cross-section |
| `sweep_resonance_1d.sh` | **1D resonance scan at B₀=0.10** (produces fig1) | `P_max = coupling²/(coupling²+Δm²²)` |
| `sweep_resonance.sh` | **2D resonance map** (300 sims, produces fig2) | Predicted line `mA² = κ²B₀²/2` |
| `sweep_mA2_family.sh` | **mA² family traces at B₀=0.10** (produces fig3) | On/off-resonance comparison |
| `reproduce_figures.sh` | **Single-command figure reproduction** | Combines all three sweeps |

## Physics motivation

### Why perturbation-level Proca?

The Proca mass term in `theory.toml` is applied to the perturbation field
`a_μ` (O(ε)) rather than to the full Maxwell field `A_μ`. This is a
deliberate physical choice:

- **Effective photon mass is dispersion, not Lagrangian mass.** An
  astrophysical plasma mass modifies the photon dispersion relation
  `ω² = k² + ωₚ²` but does not add a Lorentz-invariant mass to the
  background field carrying the static `B₀`.
- **Full-field mass breaks perturbation expansion.** Writing `-mA²/2·A·A`
  on `A = Ā + εa` expands to a term `εmA²·Ā·a` that couples the background
  to the perturbation — for a uniform `B₀` this generates position-dependent
  coefficients in the graviton EOM (root cause of issue #142 in earlier
  TIDAL derivations).
- **xPert handles the perturbation-only form cleanly.** Declaring `a` as
  O(ε) via `DefTensorPerturbation` makes `-mA²/2·a·a` an O(ε²) Lagrangian
  term that cleanly enters L^(2) without contaminating the background.

This is equivalent to the EOM-level treatment in Domcke, Garcia-Cely & Lee
2025 (`Literature/2507.16609/`) via variational calculus.

## Validation (v0.30.1)

All physics is locked in by regression tests in
[tests/test_gertsenshtein_h5_regression.py](../../tests/test_gertsenshtein_h5_regression.py)
(`TestGertsenshteinProcaRegression` class, 6 tests). Numerical validation
from the 1D sweep at B₀=0.10:

| Quantity | Observed | Theory | Agreement |
|---|---|---|---|
| Peak `P_max` (on resonance) | 1.00 | 1.00 | exact |
| HWHM (Δm²) | 0.405 | 0.402 | 0.7% |
| `P_max` at mA²=0.41 | 0.521 | 0.500 | 4% |
| `P_max` at mA²=1.00 | 0.158 | 0.141 | 12% |

2D sweep HWHM tracks theory within 9% across B₀ ∈ [0.09, 0.25].
Small-B₀ saturation (P_max < 1 for B₀ < 0.05) is expected physics:
`κB₀D/2 < π/4` means insufficient propagation for full Rabi oscillation.

## Further reading

- **Full derivation**: [docs/tex/gertsenshtein_proca.tex](../../docs/tex/gertsenshtein_proca.tex)
- **Base Gertsenshtein**: [docs/tex/gertsenshtein.tex](../../docs/tex/gertsenshtein.tex)
- **Architecture (canonical pipeline)**: [docs/tex/pipeline.tex](../../docs/tex/pipeline.tex)
  §"Canonicalization of Deferred Derived Fields"
- **Primary references** (transcribed equations, no PDF): [docs/references.md](../../docs/references.md)
  "Primary Classical References (pre-arXiv)"
- **Related arXiv papers** locally stored:
  - [Literature/2310.04150/](../../Literature/2310.04150/) — Hwang & Noh 2023
  - [Literature/2406.17853/](../../Literature/2406.17853/) — Dandoy et al. 2024
  - [Literature/2507.16609/](../../Literature/2507.16609/) — Domcke et al. 2025
  - [Literature/2004.02714/](../../Literature/2004.02714/) — Ejlli 2020
  - [Literature/2405.11786/](../../Literature/2405.11786/) — Hwang & Noh 2024
