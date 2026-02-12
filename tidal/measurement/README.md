# Measurement Module

Measurement and analysis tools for coupled-field simulations produced by the
TIDAL pipeline.

## Overview

This module sits at the end of the TIDAL pipeline:

```
Lagrangian (xAct) --> JSON spec --> PDE simulation (py-pde) --> measurement
```

It extracts quantitative physics from simulation output: per-field energy
decomposition, wave conversion probability, Fourier spectral analysis, and
energy conservation diagnostics. The primary use case is measuring the
Gertsenshtein effect — the conversion of electromagnetic waves to
gravitational waves (and vice versa) in coupled field systems.

All measurements operate on `SimulationData`, a uniform abstraction over
the full time history of field and momentum arrays. Data can come from a
live `MemoryStorage` (straight from the solver) or from a saved `.npz` file.

## Architecture

```
tidal/measurement/
    __init__.py           Public API (16 exports)
    _io.py                SimulationData: uniform data abstraction
    _energy.py            Per-field canonical energy, interaction energy
    _conversion.py        Conversion probability P(t) = E_target(t) / E_source(0)
    _spectral.py          FFT spectral decomposition, mode tracking
    _diagnostics.py       Energy conservation checks, summary statistics
```

Data flows through the module as:

```
MemoryStorage / NPZ file
        |
        v
  SimulationData          <-- _io.py: extracts full time history
        |
   +---------+-----------+-----------+
   |         |           |           |
   v         v           v           v
 Energy   Conversion   Spectral   Diagnostics
 (_energy)  (_conversion) (_spectral) (_diagnostics)
```

## Quick Start

```python
from tidal.measurement import (
    SimulationData,
    compute_conversion_probability,
    check_energy_conservation,
    compute_spectrum,
)

# From a live simulation
data = SimulationData.from_storage(storage, spec, grid, parameters)

# Or from a saved NPZ file (reads grid metadata automatically)
data = SimulationData.from_npz_auto("output.npz", spec)

# Measure wave conversion (primary Gertsenshtein observable)
result = compute_conversion_probability(data, "phi_0", "chi_0")
print(f"Peak conversion: {result.probability.max():.6f}")

# Check energy conservation
diag = check_energy_conservation(data, threshold=1e-3)
print(f"Conserved: {diag.is_conserved}, max drift: {diag.max_relative_error:.2e}")

# Spectral decomposition at a single snapshot
snap = compute_spectrum(data.fields["phi_0"][0], data.grid_spacing, data.periodic)
```

## Core Data Structures

### SimulationData

Frozen dataclass storing the **full time history** of a simulation. Both
`fields` and `momenta` have shape `(n_snapshots, *grid_shape)` — one
spatial array per recorded time. This captures transient dynamics like
Rabi oscillations where peak conversion may occur mid-simulation.

| Field | Type | Description |
|-------|------|-------------|
| `times` | `ndarray (n,)` | Snapshot times |
| `fields` | `dict[str, ndarray]` | `name -> (n, *grid)` field arrays |
| `momenta` | `dict[str, ndarray]` | `name -> (n, *grid)` conjugate momenta |
| `grid_spacing` | `tuple[float, ...]` | Cell size per axis `(dx, dy, ...)` |
| `grid_bounds` | `tuple[tuple[float, float], ...]` | Domain bounds per axis |
| `periodic` | `tuple[bool, ...]` | Per-axis periodicity flags |
| `spec` | `EquationSystem` | JSON-derived equation specification |
| `parameters` | `dict[str, float]` | Resolved parameter values |

**Constructors:**

- `from_storage(storage, spec, grid, parameters)` — from live `MemoryStorage`
- `from_npz(path, spec, grid_spacing, grid_bounds, periodic, parameters)` — from `.npz` with explicit grid parameters
- `from_npz_auto(path, spec)` — from `.npz` with grid metadata and parameters read from the file

**Properties:** `n_snapshots`, `volume_element`, `dynamical_fields`

### FieldEnergy

Energy decomposition for a single field at one snapshot.

| Field | Description |
|-------|-------------|
| `kinetic` | `0.5 * integral(pi^2) dV` |
| `gradient` | `0.5 * integral(|grad(phi)|^2) dV` |
| `mass` | `0.5 * m^2 * integral(phi^2) dV` |
| `total` | Sum of kinetic + gradient + mass |

### SystemEnergy

Energy for the full coupled system at one snapshot.

| Field | Description |
|-------|-------------|
| `per_field` | `dict[str, FieldEnergy]` per dynamical field |
| `interaction` | `0.5 * sum_{i!=j} C_ij * integral(phi_i * phi_j) dV` |
| `total` | Sum of all per-field energies + interaction |

### ConversionResult

Result of a conversion probability measurement.

| Field | Description |
|-------|-------------|
| `times` | Snapshot times |
| `probability` | `P(t) = E_target(t) / E_source(0)` |
| `source_energy` | Source field energy over time |
| `target_energy` | Target field energy over time |
| `total_energy` | Source + target energy over time |
| `relative_energy_error` | `(E_total(t) - E_total(0)) / E_total(0)` |
| `source_field` | Name of the source field |
| `target_field` | Name of the target field |

### SpectralSnapshot

Fourier decomposition of a field at one time.

| Field | Description |
|-------|-------------|
| `wavenumbers` | Radially binned `|k|` values |
| `power_spectrum` | `|phi_hat(k)|^2` averaged over shells of constant `|k|` |

### EnergyDiagnostics

Energy conservation diagnostic result.

| Field | Description |
|-------|-------------|
| `times` | Snapshot times |
| `total_energy` | Total system energy at each snapshot |
| `relative_error` | `(E(t) - E(0)) / E(0)` |
| `max_relative_error` | Peak relative energy drift |
| `is_conserved` | `max_relative_error < threshold` |

## Physics

### Canonical Hamiltonian Energy

Per-field energy uses the canonical Hamiltonian for a scalar field:

```
E_i = integral dV [ 1/2 pi_i^2 + 1/2 |grad(phi_i)|^2 + 1/2 m_ii^2 phi_i^2 ]
```

where `pi_i` is the conjugate momentum (time derivative of `phi_i`), and
`m_ii` is the diagonal mass matrix entry. Constraint fields
(`time_derivative_order == 0`) are excluded — they are not dynamical
degrees of freedom.

**Gradient computation:** For periodic axes, the gradient uses spectral
(FFT) differentiation (`ik * FFT(phi)`) for exact accuracy. For
non-periodic axes, 2nd-order central finite differences are used instead.

### Interaction Energy

Coupling between fields contributes interaction energy:

```
E_int = 1/2 sum_{i != j} C_ij integral phi_i phi_j dV
```

where `C_ij = coupling_matrix[i][j]`. The matrix convention is:
`matrix[i][j] = -(coefficient of identity(field_j) in equation_i)`.
The factor of `1/2` prevents double-counting since the sum iterates
over all ordered pairs.

### Conversion Probability

The primary Gertsenshtein observable:

```
P(t) = E_target(t) / E_source(0)
```

The source field starts with nonzero energy (e.g., a Gaussian or
plane-wave initial condition). The target field starts at zero. Coupling
terms transfer energy between them over time. `P(t)` measures how much
of the original source energy has converted to the target field.

For two coupled oscillators with equal masses, this reduces to a Rabi
oscillation with beat frequency `delta_omega = omega_+ - omega_-` where
`omega_+/-` are the normal mode frequencies of the coupled system.

### Spectral Decomposition

Fields are decomposed via `numpy.fft.rfftn` (real-to-complex FFT). For
multi-dimensional grids, the power spectrum is radially averaged by
binning in shells of constant `|k| = sqrt(kx^2 + ky^2 + ...)`.

Per-mode spectral energy:

```
E(k) = 1/2 [ |pi_hat_k|^2 + (k^2 + m^2) |phi_hat_k|^2 ]
```

For non-periodic axes, a Hann window is applied before the FFT (with a
`UserWarning` noting that amplitudes are approximate).

## Design Decisions

### Energy Floor (`_ENERGY_FLOOR = 1e-12`)

Floating-point energy values from numerical integration can be tiny but
nonzero (e.g., `1e-30`). All zero-energy comparisons use a threshold
instead of `== 0.0` to avoid division by near-zero values or incorrect
branching. This affects the conversion probability denominator and the
relative energy error computation.

### Position-Dependent Coefficient Rejection

The energy formulas assume `m^2` and `C_ij` are scalar constants. If the
JSON spec has position-dependent mass or coupling terms (identity
operators with non-empty `coordinate_dependent`), the module raises
`ValueError` with a clear message rather than computing wrong results
silently. Future work will extend the energy computation to handle
spatially varying coefficients.

### Constraint Field Handling

Fields with `time_derivative_order == 0` (e.g., `A_0` in Coulomb gauge)
are constraint fields — they are determined by elliptic equations, not
evolved in time. They have no conjugate momentum and contribute no
dynamical energy. All energy functions skip these fields automatically.

### Symbolic Coefficient Resolution

Mass and coupling matrix entries may be symbolic expressions (e.g.,
`"-m2"`, `"g/2"`). The module resolves these using the `parameters` dict
passed to `SimulationData`, falling back to the pre-computed numeric
matrix if symbolic resolution fails. This allows the same JSON spec to
be used with different parameter values.

### Gradient Method Selection

| Axis Type | Method | Accuracy |
|-----------|--------|----------|
| Periodic | FFT spectral gradient (`ik * FFT`) | Exact (spectral) |
| Non-periodic | `numpy.gradient` (2nd-order central) | O(dx^2) |

## Integration with the CLI

### From a Live Simulation

After `tidal simulate` runs, the `PlotContext` object holds all
simulation data. Convert it to `SimulationData` for measurement:

```python
ctx: PlotContext = ...  # from simulation
data = ctx.to_simulation_data()
result = compute_conversion_probability(data, "phi_0", "chi_0")
```

### From a Saved NPZ File

`tidal simulate --output result.npz` saves field snapshots, momentum
snapshots, grid metadata, and resolved parameters. Load it back:

```python
from tidal.measurement import SimulationData
from tidal.symbolic.json_loader import load_equation_system

spec = load_equation_system("examples/data/coupled_scalars.json")
data = SimulationData.from_npz_auto("result.npz", spec)
```

The NPZ file contains:

| Key Pattern | Content |
|-------------|---------|
| `times` | Snapshot time values |
| `{name}_t{idx}` | Field snapshot (e.g., `phi_0_t0`) |
| `pi_{name}_t{idx}` | Momentum snapshot (e.g., `pi_phi_0_t0`) |
| `grid_spacing` | Cell sizes per axis |
| `grid_bounds` | Domain bounds per axis |
| `grid_periodic` | Per-axis periodicity flags |
| `_param_names` | Parameter name array |
| `_param_values` | Parameter value array |

## API Reference

### Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `compute_field_energy(field, momentum, m2, spacing, periodic)` | `FieldEnergy` | Energy of one field at one snapshot |
| `compute_system_energy(data, t_idx)` | `SystemEnergy` | Full system energy at snapshot `t_idx` |
| `compute_energy_timeseries(data)` | `(times, per_field, interaction, total)` | Energy at every snapshot |
| `compute_conversion_probability(data, source, target)` | `ConversionResult` | Conversion probability `P(t)` |
| `compute_spectrum(field, spacing, periodic)` | `SpectralSnapshot` | Radially-averaged power spectrum |
| `compute_spectral_energy(field, momentum, m2, spacing, periodic)` | `(wavenumbers, energy)` | Per-mode energy `E(k)` |
| `compute_mode_amplitudes(data, field_name)` | `(times, k, amplitudes)` | Track `|phi_hat(k)|` over time |
| `check_energy_conservation(data, threshold)` | `EnergyDiagnostics` | Conservation check with threshold |
| `summarize(data)` | `dict` | Full measurement summary |

### Error Handling

All functions follow the project's fail-fast convention:

- `ValueError` for invalid field names, zero source energy, position-dependent coefficients, NaN/Inf data, out-of-range indices, non-positive thresholds
- `FileNotFoundError` for missing NPZ files

## Limitations and Future Work

- **Position-dependent coefficients:** Energy computation requires constant `m^2` and `C_ij`. Spatially varying coefficients raise `ValueError`. Extending this requires integrating `m^2(x) phi^2(x)` on the grid (Phase 4).
- **CLI integration:** No `tidal measure` subcommand yet. The `PlotContext.to_simulation_data()` bridge is in place for Phase 2.
- **Per-mode spectral conversion:** `P(k, t)` — tracking which Fourier modes participate in conversion — is planned for Phase 3.
- **Plotting utilities:** Dedicated conversion curve and spectral waterfall plots are planned for Phase 5.

## Tests

31 tests in `tests/test_measurement.py` covering:

- `SimulationData` construction from `MemoryStorage` and NPZ (including `from_npz_auto`)
- Per-field and system energy computation with analytical validation
- Rabi oscillation: pointwise comparison against exact analytical `P(t)` curve
- Energy conservation for exact coupled-oscillator data
- Decoupled fields (zero coupling gives zero conversion)
- Spectral peak detection in 1D and 2D
- Position-dependent coefficient rejection
- Near-zero energy threshold (`_ENERGY_FLOOR`)
- Edge cases: empty storage, missing files, invalid field names, NaN/Inf data
