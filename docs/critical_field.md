# Critical Field Analysis

**Purpose**: Determine the minimum magnetic field strength required for full
graviton-photon conversion in a given physical setup, and compare conversion
efficiency between theories.

## Physics Motivation

The Gertsenshtein effect converts gravitons to photons (and vice versa) in a
background magnetic field. The key experimental question is: **how strong
must the magnetic field be to achieve full conversion?**

For the standard Einstein-Maxwell (E-M) theory, the conversion probability
after a wavepacket transits a localized B-field region is given by the
Boccaletti formula:

```
P_final = sin^2(kappa/2 * integral(B(z) dz))
```

Full conversion (P_final = 1) occurs at the first peak of the sin^2
oscillation, when:

```
kappa/2 * integral(B(z) dz) = pi/2
```

For a Gaussian profile B(z) = Bpeak * exp(-z^2 / 2R^2), this gives:

```
Bpeak_min = pi / (kappa * R * sqrt(pi/2))
```

For theories with additional couplings (torsion, Chern-Simons, axion-photon),
the conversion formula differs. The critical field analysis finds B_min
numerically from sweep data, enabling direct comparison between theories.

## Key Concepts

### Why P_final, not P_max

During transit through a localized B-field, energy oscillates (Rabi
oscillation) between graviton and photon. A momentary P = 1 inside the
field region is **not** full conversion -- the energy may oscillate back
before exit.

The physically meaningful quantity is **P_final**: the permanent net
conversion after the wavepacket exits the B-field region. This is what a
downstream detector would measure, and what the Boccaletti formula predicts.

For uniform B in a periodic domain, P_final = P(t_end) = sin^2(kappa * B0 *
t_end / 2), which is the conversion at the observation time.

### Amplification Factor

To compare theories, define the amplification factor:

```
amplification = B_min(E-M) / B_min(new theory)
```

at the same threshold (P_final >= 0.99). If:
- amplification > 1: new theory enhances conversion (needs less B)
- amplification < 1: new theory weakens conversion
- amplification = 1: identical to E-M

### Threshold

The default threshold is 0.99 (full conversion within numerical tolerance).
This is a physics-defined target, not an arbitrary number. Regions of
parameter space where no B field achieves P_final >= 0.99 produce B_min =
NaN, shown as blank spots on heatmaps -- which is itself informative.

## Algorithm

Given sweep results with a field-strength parameter (e.g. B0) as one of the
swept variables:

1. **Group** rows by unique outer parameter combinations (all swept params
   except the field-strength param)
2. **Sort** each group by the field-strength parameter (ascending)
3. **Scan** upward for the first crossing where `metric >= threshold`
4. **Interpolate** linearly between the two bracketing grid points for
   sub-grid accuracy (optional, enabled by default)
5. **Estimate errors** from three sources (see below)
6. **Output** a reduced `SweepResults` with the field-strength parameter
   collapsed, containing `B_min`, `inv_B_min`, and error columns

The output plugs directly into the existing `tidal plot --type sweep`
infrastructure -- no plot code changes needed. 1D line plots, 2D heatmaps,
and 3+ dimension parallel coordinates all work automatically.

## Error Estimation

Three independent error sources, combined in quadrature:

### Grid spacing error (B_min_err_grid)

Always available. Conservative half-interval:

```
err_grid = (B_hi - B_lo) / 2
```

This is the fundamental resolution limit of the sweep grid and dominates for
coarse sweeps.

### Metric propagation error (B_min_err_metric)

Available when ensemble replicates exist. Propagates replicate standard
deviation through the linear interpolation formula using the Jacobian:

```
B_min = B_lo + (T - m_lo)/(m_hi - m_lo) * dB

dB_min/dm_lo = -dB * (m_hi - T) / (m_hi - m_lo)^2
dB_min/dm_hi = -dB * (T - m_lo) / (m_hi - m_lo)^2

err_metric^2 = (dB/dm_lo)^2 * std_lo^2 + (dB/dm_hi)^2 * std_hi^2
```

For n_replicates >= 5, a more robust alternative: compute B_min independently
per replicate, take std of the distribution.

### Interpolation model error (B_min_err_interp)

Available when 3+ points exist near the crossing. Compares linear and
quadratic interpolants:

```
err_interp = |B_min_linear - B_min_quadratic|
```

The quadratic interpolant uses three points (bracketing pair + one neighbor)
to fit a parabola through the metric values.

### Combined error

```
err_B_min = sqrt(err_grid^2 + err_metric^2 + err_interp^2)
err_{1/B} = err_B / B_min^2
```

### Quality flags

Each output point carries a `crossing_quality` flag:

| Flag | Meaning |
|------|---------|
| `good` | Clean crossing with interpolation, small errors |
| `coarse` | Large grid gap at crossing (err_grid > 0.1 * B_min) |
| `edge` | Threshold exceeded at smallest parameter value |
| `none` | No crossing found (B_min = NaN) |

## Analytical Reference Formulas

Built-in E-M reference formulas compute B_min analytically for known
geometries, useful for comparison:

| Formula | Geometry | P_EM formula | Required params |
|---------|----------|-------------|-----------------|
| `boccaletti` | Gaussian B(z) | sin^2(kappa/2 * B * R * sqrt(pi/2)) | kappa, R |
| `uniform` | Uniform B, periodic | sin^2(kappa * B * t_end / 2) | kappa, t_end |

For arbitrary geometries (dipolar, 2D/3D), use `--threshold` with a value
computed from a reference simulation or literature.

## Usage

### Single theory characterization

```bash
# 1. Run a parameter sweep including the field-strength parameter
uv run tidal sweep examples/data/gertsenshtein_localized.json \
  --sweep "Bpeak=0.01:0.5:30" --sweep "R=1:20:10" \
  --measure peak_conversion --output /tmp/sweep_torsion

# 2. Extract critical field -- collapses Bpeak dimension
#    Uses Boccaletti formula at B_ref=0.1 to compute threshold
uv run tidal analyze /tmp/sweep_torsion \
  --critical-field Bpeak \
  --reference-formula boccaletti --reference-B 0.1 \
  --output /tmp/critical_torsion

# 3. Plot 1/B_min heatmap (outer params determine plot type)
#    1 outer param -> line plot
#    2 outer params -> 2D heatmap
#    3+ outer params -> parallel coordinates
uv run tidal plot /tmp/critical_torsion --type sweep --metric inv_B_min
```

### Theory comparison (amplification factor)

```bash
# Run both E-M and new theory sweeps with identical parameter grids
uv run tidal sweep examples/data/gertsenshtein_localized.json \
  --sweep "Bpeak=0.01:0.5:30" --sweep "R=1:20:10" \
  --measure peak_conversion --output /tmp/sweep_em

uv run tidal sweep examples/data/torsion_gertsenshtein.json \
  --sweep "Bpeak=0.01:0.5:30" --sweep "R=1:20:10" \
  --measure peak_conversion --output /tmp/sweep_torsion

# Extract B_min for both with the same threshold
uv run tidal analyze /tmp/sweep_em \
  --critical-field Bpeak --threshold 0.99 --output /tmp/critical_em
uv run tidal analyze /tmp/sweep_torsion \
  --critical-field Bpeak --threshold 0.99 --output /tmp/critical_torsion

# Compare the results.csv files:
# amplification(R) = B_min_em(R) / B_min_torsion(R)
```

### Simple numeric threshold

```bash
# Directly specify a threshold (e.g. from literature)
uv run tidal analyze /tmp/sweep_torsion \
  --critical-field Bpeak --threshold 0.5 --output /tmp/critical_simple
```

## CLI Reference

```
tidal analyze SWEEP_DIR --critical-field PARAM [options]

Required:
  SWEEP_DIR               Path to sweep output directory
  --critical-field PARAM  Field-strength parameter to threshold on
  --output DIR            Output directory for reduced results

Options:
  --threshold FLOAT       Metric threshold (default: 0.99)
  --metric STR            Metric to threshold on (default: P_final)
  --no-interpolate        Disable linear interpolation
  --reference-formula     {boccaletti,uniform}
                          Compute threshold from E-M analytical formula
  --reference-B FLOAT     Reference B value for analytical formula
```

## Output Files

The output directory contains standard sweep result files:

| File | Contents |
|------|----------|
| `results.csv` | One row per outer-param combination, with B_min, inv_B_min, errors |
| `results.json` | Same data in JSON format |
| `sweep.json` | Provenance metadata (original sweep path, threshold, metric) |

### Output columns

| Column | Description |
|--------|-------------|
| (outer params) | Remaining swept parameter values |
| `B_min` | Minimum field strength for conversion >= threshold |
| `inv_B_min` | 1/B_min (default heatmap metric -- brighter = easier conversion) |
| `B_min_err` | Combined error (quadrature of grid + interp + metric) |
| `B_min_err_grid` | Error from sweep grid spacing |
| `B_min_err_metric` | Error from replicate variance (NaN if no replicates) |
| `B_min_err_interp` | Error from interpolation model |
| `inv_B_min_err` | Propagated error on 1/B_min |
| `crossing_quality` | Quality flag: good, coarse, edge, none |

## Implementation

### Key files

| File | Role |
|------|------|
| `tidal/measurement/_critical_field.py` | Core algorithm: crossing detection, error estimation, analytical formulas |
| `tidal/cli/_analyze.py` | CLI dispatch for `tidal analyze --critical-field` |
| `tidal/cli/__init__.py` | CLI argument definitions |
| `tests/test_critical_field.py` | 32 unit tests covering all edge cases |

### Architecture

The critical field analysis is a **post-hoc cross-run reduction**, not a
per-simulation measurement. It operates on completed sweep data:

```
SweepResults (N swept params)
    -> compute_critical_field()
        -> CriticalFieldResult
            -> critical_field_to_sweep_results()
                -> SweepResults (N-1 swept params)
                    -> existing tidal plot infrastructure
```

The key design insight: by producing a standard `SweepResults`, the entire
existing plot infrastructure (1D line, 2D heatmap, parallel coordinates)
works automatically without any plot code changes.

## Future Work

- **Arbitrary geometry formulas**: User-provided analytical expressions for
  reference thresholds (tracked in GitHub issue)
- **Paired comparison mode**: `--reference DIR` to load a reference sweep and
  compute amplification factors directly
- **Transit completeness check**: Detect whether the wavepacket has fully
  exited the B-field region before trusting P_final
- **Replicate-aware B_min**: Compute B_min per replicate for robust
  uncertainty estimation when ensemble data is available

## References

1. Boccaletti, D. et al. (1970) Nuovo Cim. 70B, 129-146
2. Dandoy, V. et al. (2024) arXiv:2406.17853
3. Domcke, V. & Garcia-Cely, C. (2023) arXiv:2301.02072
