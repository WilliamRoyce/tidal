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
gravitational waves (and vice versa) in coupled field systems
(Gertsenshtein 1962; Domcke & Garcia-Cely 2023, [arXiv:2301.02072](https://arxiv.org/abs/2301.02072)).

All measurements operate on `SimulationData`, a uniform abstraction over
the full time history of field and momentum arrays. Data can come from a
live `MemoryStorage` (straight from the solver) or from a saved snapshot directory.

## Architecture

```
tidal/measurement/
    __init__.py                Public API (25 exports)
    _io.py                     SimulationData: uniform data abstraction
    _energy.py                 Per-field canonical energy, interaction energy
    _conversion.py             Conversion probability P(t) = E_target(t) / E_source(0)
    _mixing.py                 Mixing length extraction and mixing spectrum
    _spectral.py               FFT spectral decomposition, mode tracking
    _spectral_conversion.py    Per-mode spectral conversion P(k,t)
    _dispersion.py             Dispersion relation omega(k) extraction
    _diagnostics.py            Energy conservation checks, summary statistics
    _utils.py                  Shared utilities (_normalize_group)
```

Data flows through the module as:

```
MemoryStorage / Snapshot directory
        |
        v
  SimulationData          <-- _io.py: extracts full time history
        |
   +------+-------+-------+---------+-------------+----------+-----------+
   |      |       |       |         |             |          |           |
   v      v       v       v         v             v          v           v
 Energy  Conv   Mixing  Spectral  SpectralConv  Dispersion  Diagnostics
 (_energy)      (_mix)  (_spec)   (_spec_conv)  (_disp)     (_diagnostics)
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

# Or from a saved snapshot directory (memory-mapped, O(1) RAM)
data = SimulationData.from_directory("output_dir", spec)

# Measure wave conversion (primary Gertsenshtein observable)
result = compute_conversion_probability(data, "phi_0", "chi_0")
print(f"Peak conversion: {result.probability.max():.6f}")

# Extract mixing length (model-independent)
from tidal.measurement import compute_mixing_length, compute_mixing_spectrum
mixing = compute_mixing_length(result)
print(f"L_mix = {mixing.mixing_length:.4f} +/- {mixing.mixing_length_uncertainty:.4f}")
print(f"Dominant frequency: {mixing.dominant_frequency:.4f}")

# Mixing spectrum — which oscillation frequencies participate?
spectrum = compute_mixing_spectrum(result)
print(f"Dominant freq: {spectrum.dominant_frequency:.4f}")

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
- `from_directory(path, spec)` — from snapshot directory (memory-mapped, O(1) RAM)
- `load(path, spec)` — universal entry point (delegates to `from_directory`)

**Methods:**

- `save(path)` — save to snapshot directory (inverse of `load()`)

**Properties:** `n_snapshots`, `volume_element`, `dynamical_fields`

### FieldEnergy

Energy density decomposition for a single field at one snapshot.
All values are spatially-averaged energy densities (intensive quantities).

| Field | Description |
|-------|-------------|
| `kinetic` | `0.5 * mean(pi^2)` |
| `gradient` | `0.5 * mean(|grad(phi)|^2)` |
| `mass` | `0.5 * mean(m^2 * phi^2)` |
| `total` | Sum of kinetic + gradient + mass |

### SystemEnergy

Energy density for the full coupled system at one snapshot.
All values are spatially-averaged energy densities.

| Field | Description |
|-------|-------------|
| `per_field` | `dict[str, FieldEnergy]` per dynamical field |
| `interaction` | `v_total - sum(per_field self-potentials)` — captures all coupling types |
| `total` | Complete Hamiltonian density: `kinetic + v_virial + v_constraint_self` |

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

### MixingResult

Spectral mixing length extracted from the dominant peak of the temporal FFT of
P(t).  This correctly identifies the physically meaningful oscillation timescale
even in multi-scale systems where rapid noise oscillations sit on top of a
slower mixing envelope.

| Field | Description |
|-------|-------------|
| `mixing_length` | `pi / omega_dom` — half-period of the dominant oscillation |
| `mixing_length_uncertainty` | Propagated from HWHM: `(pi / omega^2) * (FWHM/2)` |
| `dominant_frequency` | `omega_dom` — angular frequency of strongest spectral peak |
| `frequency_fwhm` | FWHM of the dominant peak (rad/time) |
| `max_conversion` | `max(P(t))` — peak conversion probability |
| `peaks` | All detected `SpectralPeak` entries, sorted by power descending |

### SpectralPeak

A detected peak in the mixing power spectrum.  Each peak represents a frequency
at which P(t) oscillates.

| Field | Description |
|-------|-------------|
| `frequency` | Angular frequency `omega` (rad/time) |
| `power` | `|P_hat(omega)|^2` — spectral power |
| `mixing_length` | `pi / omega` — half-period at this frequency |
| `fwhm` | Full width at half maximum (rad/time) |
| `mixing_length_uncertainty` | `(pi / omega^2) * (FWHM/2)` |

### MixingSpectrum

Temporal frequency decomposition of P(t).

| Field | Description |
|-------|-------------|
| `frequencies` | Angular frequencies `omega` of P(t) oscillation (excl. DC) |
| `power` | `|P_hat(omega)|^2` at each frequency |
| `dominant_frequency` | `omega` of strongest oscillation peak |
| `dominant_mixing_length` | `pi / dominant_frequency` — half-period at dominant freq |
| `rayleigh_resolution` | `2*pi/T` — fundamental frequency resolution from observation window |

### SpectralSnapshot

Fourier decomposition of a field at one time.

| Field | Description |
|-------|-------------|
| `wavenumbers` | Radially binned `|k|` values |
| `power_spectrum` | `|phi_hat(k)|^2` averaged over shells of constant `|k|` |

### SpectralConversion

Per-mode spectral conversion probability P(k,t).  Tracks which Fourier
modes participate in energy conversion between source and target fields.

| Field | Type | Description |
|-------|------|-------------|
| `times` | `ndarray (n,)` | Snapshot times |
| `wavenumbers` | `ndarray (n_modes,)` | Radially binned `|k|` values |
| `probability` | `ndarray (n, n_modes)` | `P(k,t) = E_target(k,t) / E_source(k,0)` |
| `source_spectral_energy` | `ndarray (n, n_modes)` | Source spectral energy per mode over time |
| `target_spectral_energy` | `ndarray (n, n_modes)` | Target spectral energy per mode over time |
| `active_modes` | `ndarray (n_modes,)` bool | `True` for k-bins with nonzero initial source energy |
| `source_field` | `str` | Source field name (or comma-joined group) |
| `target_field` | `str` | Target field name (or comma-joined group) |

### DispersionResult

Dispersion relation extracted from simulation output via spacetime FFT.

| Field | Type | Description |
|-------|------|-------------|
| `wavenumbers` | `ndarray (n_modes,)` | Radially binned `|k|` values |
| `frequencies` | `ndarray (n_freq,)` | Angular frequencies `omega` (rad/time), excluding DC |
| `power` | `ndarray (n_modes, n_freq)` | `S(k, omega) = |A_hat(k, omega)|^2` |
| `peak_frequencies` | `ndarray (n_modes,)` | Dominant `omega` at each k-bin (0 for inactive modes) |
| `peak_powers` | `ndarray (n_modes,)` | Spectral power at dominant frequency per k-bin |
| `field_name` | `str` | Which field this dispersion was computed for |
| `rayleigh_resolution` | `float` | `2*pi/T` — minimum resolvable frequency difference |

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

> For the full mathematical derivation, see [HAMILTONIAN.md](HAMILTONIAN.md).

### Hamiltonian Energy Density

The system Hamiltonian density is reconstructed automatically from the
Euler-Lagrange equations in the JSON spec — no manual per-system
formulas.  All values are **spatially-averaged energy densities**
(intensive quantities, independent of domain size).  The complete
formula is:

```
<epsilon> = 1/2 sum_{dynamical} <pi_sim^2>       (kinetic density)
          + <v_virial>                             (potential from EOM)
          + <v_constraint_self>                    (temporal components)
```

**Per-field energy density** uses the standard scalar decomposition:

```
<epsilon_i> = 1/2 <pi_i^2> + 1/2 <|grad(phi_i)|^2> + 1/2 m_ii^2 <phi_i^2>
```

where `pi_i` is the simulation momentum (time derivative of `phi_i`),
`m_ii` is the diagonal mass matrix entry, and `<...>` denotes the
spatial average (mean over grid points). Constraint fields
(`time_derivative_order == 0`) are excluded from per-field energy.

**Gradient computation:** For periodic axes, the gradient uses spectral
(FFT) differentiation (`ik * FFT(phi)`) for exact accuracy (following
the approach used in spectral PDE solvers such as Dedalus; Burns et al.
2020). For non-periodic axes, 2nd-order central finite differences are
used instead.

### Virial Potential Density

The virial potential density captures ALL cross-field coupling (identity,
derivative, constraint-mediated) using Euler's homogeneous function
theorem for degree-2 functionals:

```
<v_virial> = -1/2 sum_{i: dynamical} <phi_i * RHS_i^{spatial}>
```

where `RHS_i^{spatial}` is the right-hand side of the i-th dynamical
equation with `first_derivative_t` (gyroscopic) and `pi_N`
(velocity-dependent) terms excluded.

### Constraint Field Self-Energy Density

Temporal gauge components (A_0 in electrodynamics) have **negative**
self-energy due to the Minkowski metric g^{00} = -1:

```
<v_constraint> = sum_{j: constraint} [-1/2 <|grad(C_j)|^2>
                                      -1/2 m_j^2 <C_j^2>]
```

### Interaction Energy

Interaction energy is defined as the total potential minus each
dynamical field's self-potential:

```
E_int = V_virial + V_constraint_self - sum_i (gradient_i + mass_i)
```

This automatically captures identity coupling, derivative coupling
(e.g., `gSV * grad(phi) . A`), Chern-Simons terms, and
constraint-mediated interactions.

### Canonical vs. Simulation Momenta (Important!)

The simulation stores `pi_sim = d_t phi` for all fields.  For scalar
fields, this IS the canonical momentum.  For vector gauge fields,
`pi_canonical = F^{0i} = d_t A_i - d_i A_0 != pi_sim`.

However, when computing H via the Legendre transform, the cross-terms
`(d_i A_0)(d_t A_i)` cancel exactly between the `pi * phi_dot` and
`-L` parts.  The result uses `pi_sim` (not `pi_canonical`) with the
gauge correction absorbed into `V_constraint_self`.  See
[HAMILTONIAN.md](HAMILTONIAN.md) Section 2 for the full derivation.

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

### Mixing Spectrum and Mixing Length

The mixing spectrum is the temporal FFT of the conversion probability
`P(t)`, revealing which oscillation frequencies participate in energy
exchange between coupled fields. For a two-field Rabi system, the
dominant spectral peak sits at the beat frequency `delta_omega =
omega_+ - omega_-`; multi-field or multi-scale systems show richer
spectra with multiple peaks.

The **mixing length** `L_mix = pi / omega_dom` is the half-period of
the dominant oscillation — the characteristic distance (or time) over
which energy transfers between fields.

**Uncertainty estimation:** The uncertainty in `L_mix` is propagated
from the half-width at half-maximum (HWHM = FWHM/2) of the dominant
spectral peak via `dL = (pi / omega^2) * HWHM`. HWHM is the standard
spectroscopic convention for peak position uncertainty — it represents
the distance from the peak center to the half-power point. A sharp
peak (narrow FWHM) gives a precise mixing length; a broad peak means
the oscillation frequency is less well-defined.

The FWHM is floored at the Rayleigh resolution (`2*pi/T`) since no
measurement can resolve features narrower than the fundamental
resolution limit.  To improve frequency resolution, increase the
simulation duration `T` or decrease the snapshot interval `dt`.

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

### Spectral Conversion P(k,t)

The per-mode spectral conversion probability decomposes the scalar
conversion `P(t)` into a 2D array `P(k,t) = E_target(k,t) / E_source(k, t=0)`,
revealing which Fourier modes participate in energy conversion.

For systems with derivative coupling (gradient/cross-derivative cross-field
terms), the mixing angle `theta(k)` is k-dependent, so different modes have
different oscillation amplitudes — not just different frequencies. The
spectral conversion makes this k-dependence visible.

CLI: `tidal measure result_dir/ --what spectral_conversion --source phi_0 --target chi_0`

### Dispersion Relation

The dispersion relation `omega(k)` maps wavenumber to oscillation frequency,
encoding the dynamics of wave propagation. For a free Klein-Gordon field,
`omega(k) = sqrt(k^2 + m^2)`.

**Algorithm:**
1. Compute spatial FFT (`rfftn`) at each snapshot, producing complex
   coefficients `phi_hat(k, t)`
2. Apply temporal FFT to `phi_hat(k, t)` along the time axis, yielding
   the 2D power spectrum `S(k, omega) = |FFT_t[phi_hat(k,t)]|^2`
3. Radially bin `S(k, omega)` by `|k|`
4. Detect dominant frequency per k-bin via `argmax`

**Key insight:** The temporal FFT operates on *complex* spatial Fourier
coefficients, not on amplitudes `|phi_hat|`. Using amplitudes (always
positive) would produce frequency doubling: the temporal FFT of
`|cos(omega*t)|` has a peak at `2*omega`, not `omega`.

**Resolution limit:** The minimum resolvable frequency difference is the
Rayleigh resolution `2*pi/T` where `T` is the simulation duration. Longer
simulations improve frequency resolution.

CLI: `tidal measure result_dir/ --what dispersion --source phi_0`

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
evolved in time. They have no conjugate momentum and are excluded from
per-field energy. However, they DO contribute to the total system energy
through `V_constraint_self` (negative gradient and mass energy due to
the Minkowski metric g^{00} = -1). Cross-field terms involving
constraint fields in dynamical equations are captured automatically by
the virial formula.

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
| Non-periodic | Central differences with Dirichlet ghost cells | O(dx^2) |

For non-periodic axes, the derivative stencils use anti-symmetric ghost
cell padding (`f_ghost = -f_interior`) to match py-pde's Dirichlet BC
convention. This gives `(f[i+1] - f[i-1]) / (2dx)` for first
derivatives and `(f[i+1] - 2f[i] + f[i-1]) / dx^2` for second
derivatives, including at boundary cells where the ghost cell is implied.

## Where Simulation Data Is Saved

| Command | Persisted to disk? | Location |
|---------|-------------------|----------|
| `tidal simulate spec.json` | PNG plot only | `{spec_dir}/{stem}_output.png` |
| `tidal simulate spec.json --no-plot` | Nothing | Text summary to stdout |
| `tidal simulate spec.json --output mydir` | Snapshot directory | `mydir/` |
| `tidal simulate spec.json --output plot.png` | PNG plot only | `plot.png` |

Data is only saved to a snapshot directory when `--output` is a path without
an image extension. The example `measure_conversion.py` scripts save to
`examples/data/{name}_output/` using `create_snapshot_callback`.

## Integration with the CLI

### From a Live Simulation

After `tidal simulate` runs, the `PlotContext` object holds all
simulation data. Convert it to `SimulationData` for measurement:

```python
ctx: PlotContext = ...  # from simulation
data = ctx.to_simulation_data()
result = compute_conversion_probability(data, "phi_0", "chi_0")
```

### From a Snapshot Directory

`tidal simulate --output result_dir` streams field snapshots to disk
via `SnapshotWriter` (O(1) memory). Load it back with memory-mapped
arrays:

```python
from tidal.measurement import SimulationData
from tidal.symbolic.json_loader import load_equation_system

spec = load_equation_system("examples/data/coupled_scalars.json")
data = SimulationData.from_directory("result_dir", spec)
# Or equivalently:
data = SimulationData.load("result_dir", spec)
```

The directory contains:

| File | Content |
|------|---------|
| `metadata.json` | Grid shape, spacing, bounds, parameters, field list |
| `times.npy` | Snapshot time values, shape `(n_snapshots,)` |
| `{name}.npy` | Field snapshots, shape `(n_snapshots, *grid_shape)` |
| `pi_{name}.npy` | Momentum snapshots, shape `(n_snapshots, *grid_shape)` |

## API Reference

### Functions

| Function | Returns | Description |
|----------|---------|-------------|
| `compute_field_energy(field, momentum, m2, spacing, periodic, gradient_axes=None)` | `FieldEnergy` | Energy of one field; `gradient_axes` selects operator-aware axes |
| `compute_system_energy(data, t_idx)` | `SystemEnergy` | Full system energy at snapshot (operator-aware gradient per field) |
| `compute_energy_timeseries(data)` | `(times, per_field, interaction, total)` | Energy at every snapshot |
| `compute_conversion_probability(data, source, target)` | `ConversionResult` | Conversion probability `P(t)` |
| `compute_group_conversion(data, source_fields, target_fields)` | `ConversionResult` | Multi-field group conversion |
| `compute_mixing_length(conversion, min_prominence=0.01)` | `MixingResult` | Mixing length from dominant spectral peak of P(t) |
| `compute_mixing_spectrum(conversion)` | `MixingSpectrum` | Temporal FFT of P(t) — all mixing frequencies |
| `compute_spectrum(field, spacing, periodic)` | `SpectralSnapshot` | Radially-averaged power spectrum |
| `compute_spectral_energy(field, momentum, m2, spacing, periodic)` | `(wavenumbers, energy)` | Per-mode energy `E(k)` |
| `compute_spectral_conversion(data, source, target)` | `SpectralConversion` | Per-mode conversion `P(k,t) = E_target(k,t) / E_source(k,0)` |
| `compute_group_spectral_conversion(data, source_fields, target_fields)` | `SpectralConversion` | Multi-field group spectral conversion |
| `compute_mode_amplitudes(data, field_name)` | `(times, k, amplitudes)` | Track `|phi_hat(k)|` over time |
| `compute_dispersion(data, field_name, min_amplitude=1e-12)` | `DispersionResult` | Dispersion relation `omega(k)` via spacetime FFT |
| `check_energy_conservation(data, threshold)` | `EnergyDiagnostics` | Conservation check with threshold |
| `create_snapshot_callback(output_dir, spec, grid, ...)` | `(SnapshotWriter, callback)` | Create disk-backed snapshot callback for simulation |
| `compute_snapshot_count(t_end, snapshot_interval)` | `int` | Exact snapshot count (includes t=0) |
| `summarize(data)` | `dict` | Full measurement summary |

### Error Handling

All functions follow the project's fail-fast convention:

- `ValueError` for invalid field names, zero source energy, position-dependent coefficients, NaN/Inf data, out-of-range indices, non-positive thresholds
- `FileNotFoundError` for missing snapshot directories

## Limitations

- **Dirichlet BCs + cross_derivative operators:** The discrete ghost-cell stencil breaks self-adjointness at boundary cells, causing ~30% energy drift. Periodic BCs conserve to ~1e-10. See [HAMILTONIAN.md](HAMILTONIAN.md) Section 7, item 5 for the full analysis and SBP as a future remedy.
- **Position-dependent coefficients:** Energy computation requires constant `m²` and coupling coefficients. Spatially varying coefficients raise `ValueError`. Future: evaluate at each grid point during virial integration.
- **Quadratic Lagrangians only:** The virial formula is exact for degree-2 potentials. Higher-order Lagrangians would need explicit potential density integration.

## References

- Gertsenshtein (1962), "Wave resonance of light and gravitational waves", JETP 14, 84
- Domcke & Garcia-Cely (2023), "A simple derivation of the Gertsenshtein effect", [arXiv:2301.02072](https://arxiv.org/abs/2301.02072)
- Burns et al. (2020), "Dedalus: A Flexible Framework for Numerical Simulations with Spectral Methods", Phys. Rev. Research 2, 023068, [arXiv:1905.10388](https://arxiv.org/abs/1905.10388)
- Zwicker (2020), "py-pde: A Python package for solving partial differential equations", JOSS 5(48), 2158

See [`docs/references.md`](../../docs/references.md) for the full citation list.

## Tests

150+ tests in `tests/test_measurement.py` and `tests/test_cli.py` covering:

- `SimulationData` construction from `MemoryStorage` and snapshot directories
- Per-field and system energy computation with analytical validation
- Rabi oscillation: pointwise comparison against exact analytical `P(t)` curve
- Energy conservation for exact coupled-oscillator data
- Decoupled fields (zero coupling gives zero conversion)
- Spectral peak detection in 1D and 2D
- Position-dependent coefficient rejection
- Near-zero energy threshold (`_ENERGY_FLOOR`)
- Edge cases: empty storage, missing files, invalid field names, NaN/Inf data
