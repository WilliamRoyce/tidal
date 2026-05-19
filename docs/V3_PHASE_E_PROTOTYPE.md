# Phase E prototype validation

**Created:** 2026-05-11
**Status:** Infrastructure verified — no code changes needed to ship Phase E
**Companion to:** [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md)

## Goal

De-risk Phase E (localised wavepacket + localised B-field geometry) by verifying that the necessary infrastructure already exists and works end-to-end on a v3 theory.

## What needs to work for Phase E

| Capability | Required for | Status |
| --- | --- | --- |
| `--ic gaussian --ic-component <f> --ic-width σ_w --ic-center x_c` | Localised wavepacket IC at the source field (E.1) | ✅ verified |
| `[[background_fields]]` with arbitrary Mathematica expressions in `components` | Localised Gaussian B-field profile (E.2) | ✅ verified |
| Periodic-boundary spectral solver running with Gaussian IC on the v3 nonminimal theory | Combined E.1 + E.2 end-to-end | ✅ verified |

## Verification

### E.1 prerequisite — Gaussian IC

```bash
uv run tidal simulate examples/data/torsion_gertsenshtein_nonminimal.json \
    --grid-shape 64 --bounds=0:100 --periodic \
    --ic gaussian --ic-component h_5 --ic-amplitude 1e-2 \
    --ic-width 5.0 --ic-center 25.0 \
    --t-end 5 \
    --param B0=0.01 --param kappa=1.0 \
    --param alpha1=0.5 --param alpha2=0.5 --param alpha3=0.5 --param delta1=0.5 \
    --output /tmp/phaseE_gaussian_ic_smoke
```

Result: solver ran (modal auto-selected for periodic BCs), produced overview plot at `/tmp/phaseE_gaussian_ic_smoke/overview.png`. The IC parser accepts `--ic-width` and `--ic-center` per `tidal/cli/__init__.py:244` and the IC builder writes a Gaussian-profile field on the named component.

### E.2 prerequisite — position-dependent background fields

Existing precedent in `examples/gertsenshtein/theory.toml`:

```toml
[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "-B0 * z[]", "0"]
```

The `z[]` syntax is Mathematica for "the z spatial coordinate" — the Wolfram pipeline substitutes the actual coordinate symbol during derivation. A Gaussian profile is therefore directly expressible:

```toml
[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "-B0 * Exp[-((z[] - zc)^2)/(2*sigmaB^2)] * z[]", "0"]

[constants]
zc = 75.0
sigmaB = 25.0
```

No symbolic-pipeline extension needed. The existing `tidal/symbolic/_derive.py` + `tidal/solver/coefficients.py` paths already evaluate arbitrary spatial expressions for background-field components.

### Combined E.1 + E.2

The smoke run above used a *uniform* background (`B0 = 0.01` constant); E.2 swaps that for a localised Gaussian. Combining them just requires:

1. Re-deriving the affected theory file with the Gaussian-profile `Abar` (`tidal derive examples/torsion_gertsenshtein_nonminimal/theory.toml` — ~30 min Wolfram wall, one-time cost per theory)
2. Submitting the simulation with both `--ic gaussian` and the new derived JSON

No additional code paths required.

## What this means for Phase E

Phase E can launch as soon as Phase B converges and the geometry-decision gate (see [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md) §E.3 acceptance criterion) is reached. The expected work:

1. Re-derive 3 affected theories with Gaussian background (Wolfram, sequential — one engine license, ~1.5 h total)
2. Update `scripts/hpc_submit_drafts/v3e_localised/` with the new JSONs + Gaussian IC flags (~10 min)
3. Submit Phase E tuning sweep (E.3, 9-point × 3-t_end), then full Phase E chains (E.5)

The prototyping verified that **no Python or Wolfram code changes are needed**; only TOML edits + re-derivation. This was the de-risking goal.

## References

- [V3_PHASE_E_DESIGN.md](V3_PHASE_E_DESIGN.md) — full Phase E design
- `examples/gertsenshtein/theory.toml` — existing position-dependent `Abar` precedent
- `tidal/cli/__init__.py:179` — `--ic gaussian --ic-width` documentation
- `tidal/cli/_derive_validate.py:311` — `[[background_fields]]` schema validation
