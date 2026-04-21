# Dark Photon Torsion — Local Smoke Test

Recorded on: 2026-04-14 (TIDAL v0.31.0, post-#256 solver unification).

The dark photon torsion model uses the TorsionCDT trace Lagrangian

```
L = (1/κ²) R + α I₃ − ¼ ξ Ftorsion² + δₘ F·Ftorsion − ¼ F²
```

with constants `[kappa, B0, alpha, xi, deltam]`. The modal solver's unified
evolution-matrix builder (`_build_evolution_matrices`, see `#256`) handles the
rank-deficient mass matrix from the trace projection automatically via per-mode
mass eigendecomposition and Schur elimination.

## B.1 — Local xi sweep (5 points)

```bash
uv run tidal sweep examples/data/torsion_dark_photon.json \
  --sweep "xi=0.01:0.5:5" \
  --measure peak_conversion --source h_5 --target a_1 \
  --grid-shape 256 --bounds 0:100 --periodic \
  --ic plane-wave --ic-wavevector 2.0 --ic-amplitude 0.1 --ic-component h_5 \
  --t-end 50 \
  --param kappa=1.0 --param B0=0.01 --param alpha=0.5 --param deltam=0.1 \
  --output /tmp/tidal_tests/dp_smoke_sweep
```

**Result:** 5/5 runs `success`, total wall time ~5 s.

| xi     | P_max                   | status  |
|--------|-------------------------|---------|
| 0.01   | 0.061208312187801976    | success |
| 0.1325 | 0.061208312187784115    | success |
| 0.255  | 0.061208312187800054    | success |
| 0.3775 | 0.06120831218780377     | success |
| 0.5    | 0.06120831218780689     | success |

All `P_max ≈ 0.061208` (analytical `sin²(κB₀t_end/2) = sin²(0.25) = 0.061209`
to 4 digits). The ξ parameter affects the torsion sector but does not alter
the h_5 → a_1 Gertsenshtein channel at this nominal point — consistent with
the decoupled sector structure of the trace-vector Lagrangian.

## B.2 — t_end independence check at nominal point

```bash
# (α=0.5, ξ=1.0, δₘ=0.1, k=2)
# t_end=25: P_max = 0.015544
# t_end=50: P_max = 0.061208
# ratio = 3.938 ≈ 4 (sin² growth in the linearized regime)
```

Analytical prediction for pure Gertsenshtein channel:
`P(t) = sin²(κ B₀ t / 2)`.
- `sin²(0.01·25/2) = sin²(0.125) = 0.015552` vs measured 0.015544 → 0.05% off
- `sin²(0.01·50/2) = sin²(0.25) = 0.061209` vs measured 0.061208 → 0.002% off

**Conclusion:** the TorsionCDT trace dark photon runs cleanly through the
unified modal solver, reproduces the standard Gertsenshtein formula to
4-digit precision at the nominal `(α=0.5, ξ=1, δₘ=0.1)` point, and shows
textbook sin² growth (ratio 4) over a 2× increase in `t_end`. Ready for
parameter sweeps.

## Known cosmetic warnings (not regressions)

- `UserWarning: mass_matrix/coupling_matrix inconsistent with identity
  operator terms` from `tidal/symbolic/json_loader.py` — unrelated to
  #256, affects diagnostic output only, not the simulation.
- `Modal solver (per-mode): eigenvalues with positive real parts (max
  Re(λ)≈1.02)` — spurious gauge/null-space modes that are correctly
  frozen downstream by the gauge filter `|λ| > 1e12` and the
  `_suppress_tachyonic_noise` uncoupled-tachyon suppressor. The
  simulation runs cleanly despite the warning.
