# Next Features: Parameter Sweeps & Measurements

**Created:** March 2026
**Branch:** `feature/parameter-sweep` (commit `6150ec6`, 1242 tests passing)
**Status:** All phases (A–D) complete

This document tracks proposed features for TIDAL's parameter sweep framework and measurement modules. Each feature includes scope, generality analysis, implementation details, and acceptance criteria.

---

## Table of Contents

1. [Current State Summary](#current-state-summary)
2. [Generality Constraints](#generality-constraints)
3. [Feature F1: Sweep TOML Configuration](#f1-sweep-toml-configuration)
4. [Feature F2: Adaptive Parameter Sampling](#f2-adaptive-parameter-sampling)
5. [Feature F3: Group/Phase Velocity & Resonance Analysis](#f3-groupphase-velocity--resonance-analysis)
6. [Feature F4: Sobol Sensitivity Analysis](#f4-sobol-sensitivity-analysis)
7. [Feature F5: SweepResults Query Methods](#f5-sweepresults-query-methods)
8. [Feature F6: Spectrum Scalar Aggregation in Sweeps](#f6-spectrum-scalar-aggregation-in-sweeps)
9. [Feature F7: Run Status Tracking & Failure Classification](#f7-run-status-tracking--failure-classification)
10. [Feature F8: Advanced Sweep Visualization](#f8-advanced-sweep-visualization)
11. [Dependency Graph](#dependency-graph)
12. [Implementation Phases](#implementation-phases)

---

## Current State Summary

The `tidal sweep` command is fully operational with:

- **9 measurement types**: summary, energy, conversion, mixing, dispersion, conservation, effective_mass, asymptotic, peak_conversion
- **Sampling**: linspace, logspace, explicit values, Cartesian product of multiple parameters
- **Execution**: sequential, parallel (`multiprocessing.Pool`), resume from interruption
- **Output**: `results.csv`, `results.json`, `sweep.json` (provenance)
- **Safety**: dry-run mode, sweep size limits, parameter/measurement validation, ETA tracking
- **Convergence**: `--converge` mode with Richardson extrapolation order estimation

### Key Files

| File                                        | Purpose                                        |
| ------------------------------------------- | ---------------------------------------------- |
| `tidal/cli/_sweep.py`                       | Core sweep orchestration (~1350 lines)         |
| `tidal/cli/_sweep_panels.py`                | Sweep plot renderers (~550 lines)              |
| `tidal/measurement/_sweep_results.py`       | SweepResults dataclass (~570 lines)            |
| `tidal/cli/_measure.py`                     | Measurement CLI dispatcher                     |
| `tidal/measurement/_energy.py`              | Energy computation (~1700 lines, most general) |
| `tidal/measurement/_conversion.py`          | Conversion probability P(t)                    |
| `tidal/measurement/_mixing.py`              | Mixing length L_mix                            |
| `tidal/measurement/_asymptotic.py`          | Asymptotic scattering observables              |
| `tidal/measurement/_dispersion.py`          | Dispersion relation omega(k)                   |
| `tidal/measurement/_effective_mass.py`      | Effective mass m2_eff                          |
| `tidal/measurement/_spectral.py`            | Spatial power spectrum                         |
| `tidal/measurement/_spectral_conversion.py` | Per-mode spectral conversion P(k,t)            |
| `tidal/cli/_sweep_config.py`                | TOML sweep config parser (~475 lines)          |
| `tests/test_sweep.py`                       | Sweep tests (~1200 lines)                      |
| `tests/test_sweep_config.py`                | TOML config tests (~300 lines)                 |
| `tidal/measurement/_sensitivity.py`         | Sobol/Morris sensitivity analysis              |
| `tidal/cli/_analyze.py`                     | `tidal analyze` CLI handler                    |
| `tests/test_new_measurements.py`            | Measurement tests (~391 lines)                 |
| `tests/test_sensitivity_viz.py`             | Sensitivity + visualization tests (23 tests)   |

---

## Generality Constraints

TIDAL supports theories from flat 1+1D scalars to curved 3+1D tensor fields. Measurement modules have different generality levels:

### Tier 1: General (any spacetime, any theory)

These modules use the volume-element-aware energy computation and operate on scalar timeseries. They work for all 22 examples including curved spacetimes (polar_kg, sphere_kg, spherical_kg), background fields, and constraint systems.

| Module            | Why General                                                                                              |
| ----------------- | -------------------------------------------------------------------------------------------------------- |
| `_energy.py`      | Evaluates Hamiltonian terms with volume_element, position-dependent coefficients, mixed operators        |
| `_conversion.py`  | Computes P(t) = E_target(t) / E_source(0) using general energy                                           |
| `_mixing.py`      | Temporal FFT of P(t) timeseries — no spatial assumptions                                                 |
| `_diagnostics.py` | Energy drift from general energy computation                                                             |
| `_asymptotic.py`  | Energy-based totals are general; directional split (P_transmitted/P_reflected) uses FFT and is flat-only |

### Tier 2: Flat + Spatially Homogeneous Only

These modules use spatial FFT, which requires translation invariance (constant coefficients, periodic BCs, Cartesian coordinates). They enforce this via `_check_no_position_dependent_terms()` and raise `ValueError` for incompatible systems.

| Module                    | Why Restricted                                                 |
| ------------------------- | -------------------------------------------------------------- | ---------- | ----------------------- |
| `_dispersion.py`          | omega(k) via space-time FFT requires global Fourier eigenmodes |
| `_effective_mass.py`      | Wraps dispersion: m2_eff = omega^2 - k^2                       |
| `_spectral.py`            | Power spectrum                                                 | phi_hat(k) | ^2 requires spatial FFT |
| `_spectral_conversion.py` | Per-mode P(k,t) requires spectral energy decomposition         |

### Impact on New Features

Every proposed feature below is annotated with its generality tier:

- **[GENERAL]** — works for any theory/spacetime
- **[FLAT+HOMOGENEOUS]** — requires flat spacetime with constant coefficients
- **[GENERAL*]** — mostly general, with noted limitations for specific sub-features

When a feature is FLAT+HOMOGENEOUS, the specification notes what a general alternative would look like (e.g., wave packet tracking instead of FFT-based group velocity).

---

## F1: Sweep TOML Configuration

**Priority:** 1 (immediate)
**Generality:** [GENERAL] — infrastructure feature, theory-agnostic
**Status:** Complete (Phase A, commit `5b4019c`)
**Depends on:** Nothing

### Motivation

Current sweep invocations require long CLI commands with many flags. A TOML configuration file would:

- Make sweeps reproducible and version-controllable (commit sweep.toml alongside theory.toml)
- Allow complex multi-parameter setups without shell quoting headaches
- Follow the existing TIDAL pattern (theory.toml for derivation, sweep.toml for parameter studies)
- Enable adaptive sampling configuration (F2) which has too many options for CLI flags alone

### Specification

**Usage:** `tidal sweep --config examples/coupled_scattering/sweep_coupling.toml`

CLI flags override TOML values, so `--config sweep.toml --param g0=0.3` uses TOML for everything except g0.

**TOML schema:**

```toml
# Required: path to JSON equation spec (resolved relative to TOML file)
spec = "../data/coupled_scattering.json"

# Fixed parameters (not swept)
[parameters]
mPhi2 = 1.0
mChi2 = 1.0
R = 5.0

# Swept parameters — each [sweep.NAME] section defines one swept axis
[sweep.g0]
# Range mode (generates values via linspace/logspace):
start = 0.01
stop = 1.0
count = 10
scale = "log"          # "linear" (default) | "log" | "adaptive" (see F2)

# Alternative: explicit values mode:
# values = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0]

# Multiple swept params -> Cartesian product
# [sweep.mChi2]
# start = 0.5
# stop = 4.0
# count = 5

# Simulation settings (same flags as `tidal simulate`)
[simulation]
grid_shape = 128
bounds = [0, 100]
periodic = true
ic = "gaussian"
ic_component = "phi_0"
ic_center = 25.0
ic_width = 5.0
t_end = 20.0
# snapshots = 0.5   # optional snapshot interval

# Measurement configuration
[measurement]
types = ["conversion", "mixing", "conservation", "asymptotic"]
source = "phi_0"
target = "chi_0"
energy_threshold = 1e-3

# Output directory (resolved relative to TOML file)
[output]
path = "../data/sweep_coupling"

# Execution settings
[execution]
parallel = 4
resume = true

# Convergence mode (mutually exclusive with [sweep.*] sections)
# [convergence]
# grid_sizes = [32, 64, 128, 256]
```

### Implementation

**New file: `tidal/cli/_sweep_config.py`**

```python
@dataclass
class SweepConfig:
    spec_path: Path
    swept_params: dict[str, list[float]]
    fixed_params: dict[str, float]
    sim_settings: dict[str, Any]
    measurements: list[str]
    source: str | None
    target: str | None
    output: Path
    parallel: int | None
    resume: bool
    energy_threshold: float
    converge_sizes: list[int] | None
    # Adaptive config (for F2)
    adaptive_config: dict[str, dict[str, Any]] | None

def load_sweep_config(path: Path) -> SweepConfig:
    """Load and validate a sweep TOML configuration file."""
```

**Key design decisions:**

- Uses stdlib `tomllib` (Python 3.11+, already in project requirements)
- `spec` path resolved relative to TOML file location (same as theory.toml `[output].path`)
- Unknown TOML keys produce warnings (not errors) for forward compatibility
- Missing `spec` key is an error
- `[sweep.PARAM]` sections parsed into the same `(name, values)` tuples as `parse_sweep_spec()`
- The `scale = "adaptive"` value is reserved for F2

**Modified files:**

| File                                              | Change                                                           |
| ------------------------------------------------- | ---------------------------------------------------------------- |
| `tidal/cli/_sweep_config.py`                      | New: TOML parser + SweepConfig dataclass                         |
| `tidal/cli/_sweep.py`                             | Merge SweepConfig into existing flow at top of `sweep_command()` |
| `tidal/cli/__init__.py`                           | Add `--config` flag to sweep subparser                           |
| `tests/test_sweep_config.py`                      | New: parsing, validation, CLI override, path resolution          |
| `examples/coupled_scattering/sweep_coupling.toml` | New: example TOML for 1D coupling sweep                          |
| `examples/coupled_scalars/sweep_coupling.toml`    | New: example TOML for coupled scalars                            |

### Acceptance Criteria

- [x] `tidal sweep --config sweep.toml` runs sweep from TOML
- [x] `tidal sweep --config sweep.toml --param g0=0.3` overrides TOML parameter
- [x] `tidal sweep --config sweep.toml --dry-run` shows correct plan
- [x] Path resolution works relative to TOML file location
- [x] Unknown TOML keys warn but don't error
- [x] Missing `spec` key produces clear error
- [x] Both range mode and explicit values mode work in TOML
- [x] Convergence mode (`[convergence]`) works in TOML
- [x] Tests pass: parsing, validation, override, path resolution
- [x] Lint clean (`ruff check`)

---

## F2: Adaptive Parameter Sampling

**Priority:** 2 (immediate, after F1)
**Generality:** [GENERAL] — infrastructure feature, theory-agnostic
**Status:** Complete (F2a in Phase A commit `5b4019c`, F2b in Phase B commit `3133f49`)
**Depends on:** F1 (TOML config provides the natural interface for adaptive settings)

### Motivation

Uniform grids (linspace/logspace) waste samples in flat regions and miss sharp features. For Gertsenshtein physics, P(g0) has a sharp knee near the instability threshold g0 = sqrt(mPhi2 \* mChi2), where a few extra points around the transition are worth more than a dense uniform grid everywhere.

### F2a: Interval Refinement (1D Sweeps)

**How it works:**

1. Run initial coarse grid (e.g. 5 points, linspace between bounds)
2. Compute interest score for each interval based on metric curvature
3. Bisect the highest-interest interval, run new simulation at midpoint
4. Repeat until budget exhausted or max interval score < threshold

**Interest score for interval [a, b] with values f(a), f(b), and midpoint m with f(m):**

```
score = |f(a) - 2*f(m) + f(b)| / (|f(b) - f(a)| + epsilon)
```

This is the discrete second derivative (curvature) normalized by the first derivative. High score = sharp feature worth resolving. The midpoint f(m) is computed by running a new simulation at parameter value (a+b)/2.

**Constraint:** Adaptive mode works for 1D sweeps only (single swept parameter). Multi-parameter adaptive would require Bayesian optimization or similar, deferred to future work.

**CLI syntax:**

```bash
tidal sweep spec.json --sweep "g0=0.01:1.0:adaptive" \
  --adaptive-metric P_max \
  --adaptive-budget 30 \
  --adaptive-initial 5 \
  --adaptive-threshold 0.01
```

**TOML syntax (preferred):**

```toml
[sweep.g0]
start = 0.01
stop = 1.0
scale = "adaptive"
initial_count = 5       # coarse grid points (default: 5)
max_count = 30          # total budget including initial (default: 20)
metric = "P_max"        # metric driving refinement
threshold = 0.01        # stop when max interval score < this
```

**Implementation in `_sweep.py`:**

```python
def _execute_adaptive(
    args: Namespace,
    spec_path: Path,
    param_name: str,
    bounds: tuple[float, float],
    initial_n: int,
    max_n: int,
    metric_key: str,
    threshold: float,
    # ... remaining sweep config ...
) -> list[dict[str, Any]]:
    """Adaptive refinement sweep for a single parameter."""
    # 1. Run initial grid
    points = np.linspace(bounds[0], bounds[1], initial_n).tolist()
    rows = [_run_single_at_value(p, ...) for p in points]

    # 2. Iterative refinement
    while len(points) < max_n:
        scores = _interval_scores(points, rows, metric_key)
        if max(scores) < threshold:
            break
        idx = int(np.argmax(scores))
        midpoint = (points[idx] + points[idx + 1]) / 2
        new_row = _run_single_at_value(midpoint, ...)
        points.insert(idx + 1, midpoint)
        rows.insert(idx + 1, new_row)
        _save_incremental(...)  # crash recovery

    return rows  # sorted by parameter value
```

### F2b: Latin Hypercube / Sobol Sequences (Multi-D)

For 3+ parameter sweeps where Cartesian grids are infeasible. Uses `scipy.stats.qmc` (already a transitive dependency).

**TOML syntax:**

```toml
[sweep]
strategy = "latin_hypercube"   # or "sobol" or "grid" (default)
n_samples = 100

[sweep.g0]
start = 0.01
stop = 1.0

[sweep.mChi2]
start = 0.5
stop = 4.0

[sweep.R]
start = 1.0
stop = 20.0
```

**Implementation:**

- `_generate_samples()` function using `scipy.stats.qmc.LatinHypercube` or `Sobol`
- Returns list of parameter dicts (same format as Cartesian product output)
- `SweepResults.sampling_strategy` field for provenance
- Falls back to Cartesian product if strategy not specified

### Modified Files

| File                                  | Change                                                                                                         |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `tidal/cli/_sweep.py`                 | `_execute_adaptive()`, `_interval_scores()`, `_generate_samples()`, update `parse_sweep_spec()` for "adaptive" |
| `tidal/cli/__init__.py`               | `--adaptive-metric`, `--adaptive-budget`, `--adaptive-initial`, `--adaptive-threshold` flags                   |
| `tidal/cli/_sweep_config.py`          | Parse adaptive and strategy settings from TOML                                                                 |
| `tidal/measurement/_sweep_results.py` | `sampling_strategy` field                                                                                      |
| `tests/test_sweep.py`                 | Adaptive parsing, refinement logic, budget enforcement, LHS/Sobol generation                                   |

### Acceptance Criteria

- [x] `--sweep "g0=0.01:1.0:adaptive"` runs adaptive refinement
- [x] Adaptive concentrates points around sharp features (verifiable on coupled_scattering P(g0) curve)
- [x] Budget limit enforced (never exceeds `--adaptive-budget`)
- [x] Threshold stopping works (stops early if all intervals smooth)
- [x] Results sorted by parameter value in output CSV
- [x] Incremental save works during adaptive refinement
- [x] `strategy = "latin_hypercube"` generates space-filling samples
- [x] `strategy = "sobol"` generates low-discrepancy sequences
- [x] `SweepResults.sampling_strategy` records which strategy was used
- [x] Tests pass for parsing, refinement, budget, LHS, Sobol

---

## F3: Group/Phase Velocity & Resonance Analysis

**Priority:** 3 (after F1-F2)
**Generality:** [FLAT+HOMOGENEOUS] — FFT-based, inherits dispersion constraints
**Status:** Complete (Phase C, commit `2a6f04e`)
**Depends on:** Nothing (uses existing `_dispersion.py`)

### Motivation

For Gertsenshtein-type scattering, the key physics questions are:

1. **Which modes convert?** — Modes near resonance (omega_source(k) ~ omega_target(k)) convert efficiently
2. **How fast do wave packets travel?** — Group velocity mismatch determines coherence length
3. **What's the conversion bandwidth?** — How many modes lie within the resonance window?

Currently TIDAL extracts omega(k) from dispersion but never computes velocities or identifies resonant modes.

### Generality Note

These features inherit the flat+homogeneous constraint from `_dispersion.py`. For curved spacetimes:

- **Group velocity alternative:** Track wave packet centroid position over time (real-space). Compute v_g = dx_centroid/dt. Works for any geometry but requires well-separated wave packets.
- **Resonance alternative:** Compare energy transfer rates between fields at different spatial locations. Identify regions of maximum coupling.
- **Implementation:** These curved-spacetime alternatives are deferred to future work. The FFT-based approach covers the primary use cases (coupled_scalars, coupled_scattering, coupled_proca, proca) which are all flat.

### F3a: Group & Phase Velocity

**New module: `tidal/measurement/_velocity.py`**

```python
@dataclass
class VelocityResult:
    """Velocity analysis from dispersion relation."""
    wavenumbers: np.ndarray         # k values
    group_velocity: np.ndarray      # dw/dk per mode
    phase_velocity: np.ndarray      # w/k per mode
    group_velocity_mean: float      # weighted mean v_g
    phase_velocity_mean: float      # weighted mean v_p

def compute_velocities(
    dispersion: DispersionResult,
    *,
    smoothing_window: int = 3,
) -> VelocityResult:
    """Extract group velocity dw/dk and phase velocity w/k.

    Uses Savitzky-Golay filter on omega(k) to compute smooth
    numerical derivative for group velocity.
    """
```

**For two-field systems, add velocity mismatch:**

```python
@dataclass
class VelocityMismatchResult:
    """Velocity mismatch between two fields."""
    source_velocity: VelocityResult
    target_velocity: VelocityResult
    mismatch: np.ndarray            # |v_g_source(k) - v_g_target(k)|
    max_mismatch: float
    coherence_length: float         # L_coh ~ 1 / (k * delta_v_g)

def compute_velocity_mismatch(
    data: SimulationData,
    source_field: str,
    target_field: str,
) -> VelocityMismatchResult:
```

**Sweep scalars:** `v_group_mean`, `v_phase_mean`, `v_mismatch_max`, `coherence_length`

### F3b: Resonance Analysis

**New module: `tidal/measurement/_resonance.py`**

```python
@dataclass
class ResonanceResult:
    """Resonance analysis for coupled fields."""
    wavenumbers: np.ndarray                  # k values
    omega_source: np.ndarray                 # omega_source(k)
    omega_target: np.ndarray                 # omega_target(k)
    resonance_mismatch: np.ndarray           # |omega_s - omega_t| per mode
    resonant_modes: np.ndarray               # bool mask: which modes are near-resonant
    n_resonant_modes: int                    # count of resonant modes
    conversion_bandwidth: float              # FWHM of resonance window in k-space
    peak_conversion_k: float                 # k with max conversion efficiency
    adiabaticity_parameter: float | None     # L_coupling * delta_k (if coupling_length provided)

def compute_resonance_analysis(
    data: SimulationData,
    source_field: str,
    target_field: str,
    *,
    resonance_threshold: float = 0.1,       # |delta_omega/omega| < threshold
    coupling_length: float | None = None,    # for adiabaticity computation
) -> ResonanceResult:
    """Identify resonant modes and compute conversion bandwidth.

    A mode k is resonant when |omega_source(k) - omega_target(k)| / omega_avg(k)
    < resonance_threshold. The conversion bandwidth is the FWHM of the
    resonance mismatch curve.

    The adiabaticity parameter Psi = coupling_length * delta_k_resonance
    determines whether conversion is adiabatic (Psi >> 1) or sudden (Psi << 1).
    """
```

**Sweep scalars:** `n_resonant_modes`, `conversion_bandwidth`, `k_peak_conversion`, `adiabaticity_parameter`

### Modified Files

| File                              | Change                                                                                             |
| --------------------------------- | -------------------------------------------------------------------------------------------------- |
| `tidal/measurement/_velocity.py`  | New: velocity computation                                                                          |
| `tidal/measurement/_resonance.py` | New: resonance analysis                                                                            |
| `tidal/measurement/__init__.py`   | Export new functions                                                                               |
| `tidal/cli/_measure.py`           | Wire `velocity` and `resonance` measurement types                                                  |
| `tidal/cli/_sweep.py`             | Add `velocity` and `resonance` to `_SWEEP_MEASUREMENTS`, add scalar extraction to `_measure_run()` |
| `tests/test_velocity.py`          | New: velocity extraction, mismatch                                                                 |
| `tests/test_resonance.py`         | New: resonance detection, bandwidth, adiabaticity                                                  |

### Acceptance Criteria

- [x] `tidal measure output/ --what velocity --source phi_0` extracts group and phase velocity
- [x] `tidal measure output/ --what resonance --source phi_0 --target chi_0` identifies resonant modes
- [x] Group velocity matches analytic dw/dk for free Klein-Gordon (v_g = k/sqrt(k^2 + m^2))
- [x] Resonance correctly identifies k-modes where omega_source ~ omega_target
- [x] Conversion bandwidth narrows as mass difference increases (physical expectation)
- [ ] Adiabaticity parameter computed when coupling_length provided (deferred)
- [x] Both measurements wire into sweep with scalar summaries
- [x] Clear error message for curved-spacetime systems (position-dependent terms)

---

## F4: Sobol Sensitivity Analysis

**Priority:** 4
**Generality:** [GENERAL] — operates on sweep result scalars, no field data assumptions
**Status:** Complete (Phase D, commit `6150ec6`)
**Depends on:** F2b (requires specific sampling designs like Saltelli)

### Motivation

For multi-parameter sweeps, Sobol sensitivity indices tell you which parameters actually matter. A physicist sweeping 5 parameters (g0, mPhi2, mChi2, R, grid_shape) wants to know that coupling strength and mass ratio control 90% of the variance in conversion probability, while grid resolution and coupling width are less important.

### Specification

**Optional dependency:** [SALib](https://github.com/SALib/SALib) for Sobol/Morris analysis. Import guarded with helpful error message if not installed.

**New CLI subcommand:** `tidal analyze`

```bash
# Sobol sensitivity analysis
tidal analyze sweep_output/ --sensitivity sobol --metric P_max --bootstrap 100

# Morris screening (faster, good for many parameters)
tidal analyze sweep_output/ --sensitivity morris --metric P_max
```

**Output:**

```
Sobol Sensitivity Analysis: P_max
=================================
Parameter    S1 (main)    ST (total)    Interaction
---------    ---------    ----------    -----------
g0           0.45 +/- 0.02   0.52 +/- 0.02   moderate
mChi2        0.38 +/- 0.02   0.41 +/- 0.02   weak
mPhi2        0.08 +/- 0.01   0.12 +/- 0.01   weak
R            0.02 +/- 0.01   0.05 +/- 0.01   none
```

**New module: `tidal/measurement/_sensitivity.py`**

```python
@dataclass
class SensitivityResult:
    """Sobol or Morris sensitivity analysis result."""
    method: str                              # "sobol" or "morris"
    metric: str                              # which metric was analyzed
    param_names: list[str]
    # Sobol indices (None for Morris)
    s1: np.ndarray | None                    # first-order indices
    st: np.ndarray | None                    # total-order indices
    s1_conf: np.ndarray | None               # confidence intervals
    st_conf: np.ndarray | None
    # Morris indices (None for Sobol)
    mu_star: np.ndarray | None               # abs mean elementary effect
    sigma: np.ndarray | None                 # std of elementary effects

def compute_sobol_indices(
    results: SweepResults,
    metric: str,
    *,
    n_bootstrap: int = 100,
) -> SensitivityResult:

def compute_morris_screening(
    results: SweepResults,
    metric: str,
) -> SensitivityResult:
```

**Sampling requirement:** Sobol analysis requires specific sampling designs (Saltelli scheme). Could add `--sweep-strategy saltelli` to F2b, or post-process existing Latin hypercube/Sobol data.

### Modified Files

| File                                | Change                                      |
| ----------------------------------- | ------------------------------------------- |
| `tidal/measurement/_sensitivity.py` | New: SALib integration                      |
| `tidal/cli/__init__.py`             | New `analyze` subparser                     |
| `tidal/cli/_analyze.py`             | New: analyze command handler                |
| `tests/test_sensitivity.py`         | New: Sobol/Morris tests with synthetic data |

### Acceptance Criteria

- [x] `tidal analyze sweep_output/ --sensitivity sobol --metric P_max` produces index table
- [x] Clear error if SALib not installed
- [x] Sobol S1 + interactions = ST (within confidence intervals)
- [x] Morris screening ranks parameters by importance
- [x] Works with any SweepResults (any theory, any metric)

---

## F5: SweepResults Query Methods

**Priority:** 5
**Generality:** [GENERAL] — data structure feature
**Status:** Complete (Phase B, commit `3133f49`)
**Depends on:** Nothing

### Motivation

Currently `SweepResults` only has `.column(name)` for data extraction. Any filtering, grouping, or statistical analysis requires manual numpy/pandas work. Adding a few query methods makes programmatic analysis much easier.

### Specification

Add to `tidal/measurement/_sweep_results.py`:

```python
def filter(self, **kwargs: float) -> SweepResults:
    """Filter rows by parameter values.

    Example: results.filter(mChi2=1.0) returns only rows where mChi2 == 1.0.
    For float comparison, uses np.isclose with default tolerances.
    Returns a new SweepResults with filtered rows (preserving all other metadata).
    """

def group_by(self, param: str) -> dict[float, SweepResults]:
    """Group results by a parameter value.

    Returns dict mapping parameter values to SweepResults subsets.
    Useful for slicing a 2D sweep into 1D slices.
    """

def best(self, metric: str, *, maximize: bool = True) -> dict[str, Any]:
    """Return the row with the best (max or min) metric value.

    Example: results.best("P_max") returns the row with highest P_max.
    Example: results.best("max_energy_error", maximize=False) returns
    the row with lowest energy error.
    """

def summary(self) -> str:
    """Return a human-readable summary table.

    Shows: n_runs, parameter ranges, metric statistics (mean, std, min, max)
    for all available metrics.
    """
```

### Modified Files

| File                                  | Change                                      |
| ------------------------------------- | ------------------------------------------- |
| `tidal/measurement/_sweep_results.py` | Add filter, group_by, best, summary methods |
| `tests/test_sweep.py`                 | Tests for each new method                   |

### Acceptance Criteria

- [x] `results.filter(mChi2=1.0)` returns correct subset
- [x] `results.group_by("g0")` returns dict with correct grouping
- [x] `results.best("P_max")` returns row with highest P_max
- [x] `results.best("max_energy_error", maximize=False)` returns row with lowest error
- [x] `results.summary()` produces readable output
- [x] Filter preserves metadata (swept_params, fixed_params, etc.)
- [x] Empty filter result doesn't crash

---

## F6: Spectrum Scalar Aggregation in Sweeps

**Priority:** 6
**Generality:** [FLAT+HOMOGENEOUS] — inherits FFT constraints from `_spectral.py` and `_spectral_conversion.py`
**Status:** Complete (Phase C, commit `2a6f04e`)
**Depends on:** Nothing

### Motivation

`spectrum` and `spectral_conversion` are the two measurement types available in `tidal measure` but not in `tidal sweep` (9 out of 11 types are supported). They produce 2D arrays, which don't fit sweep's scalar-row format. But useful scalar summaries can be extracted.

### Specification

Add to `_measure_run()` in `_sweep.py`:

**For `spectrum`:**

```python
# Extract scalar summaries from spectral snapshot
spec_result = compute_spectrum(data, source_field)
metrics["n_active_modes"] = int(np.sum(spec_result.power > threshold * spec_result.power.max()))
metrics["peak_k"] = float(spec_result.wavenumbers[np.argmax(spec_result.power)])
metrics["peak_power"] = float(spec_result.power.max())
```

**For `spectral_conversion`:**

```python
# Extract scalar summaries from per-mode conversion
sc_result = compute_spectral_conversion(data, source_field, target_field)
final_P_k = sc_result.conversion[-1]  # P(k) at final time
metrics["P_k_max"] = float(final_P_k.max())
metrics["k_max_conversion"] = float(sc_result.wavenumbers[np.argmax(final_P_k)])
# Conversion bandwidth: FWHM of P(k) at final time
half_max = final_P_k.max() / 2
above_half = final_P_k > half_max
if np.any(above_half):
    k_range = sc_result.wavenumbers[above_half]
    metrics["conversion_bandwidth"] = float(k_range[-1] - k_range[0])
else:
    metrics["conversion_bandwidth"] = 0.0
```

### Modified Files

| File                    | Change                                                                                                          |
| ----------------------- | --------------------------------------------------------------------------------------------------------------- |
| `tidal/cli/_sweep.py`   | Add `"spectrum"` and `"spectral_conversion"` to `_SWEEP_MEASUREMENTS`, add extraction logic to `_measure_run()` |
| `tidal/cli/__init__.py` | Update help text to list 11 measurement types                                                                   |
| `tests/test_sweep.py`   | Test spectrum/spectral_conversion scalar extraction                                                             |

### Acceptance Criteria

- [x] `--measure spectrum` in sweep produces `n_active_modes`, `peak_k`, `peak_power` columns
- [x] `--measure spectral_conversion` produces `P_k_max`, `k_max_conversion`, `spectral_conversion_bandwidth`
- [x] Clear error for curved-spacetime systems (inherits from `_check_no_position_dependent_terms()`)
- [x] 13 measurement types now supported in sweeps (was 9, now includes velocity, resonance, spectrum, spectral_conversion)

---

## F7: Run Status Tracking & Failure Classification

**Priority:** 7
**Generality:** [GENERAL] — infrastructure feature
**Status:** Complete (Phase B, commit `3133f49`)
**Depends on:** Nothing

### Motivation

Not all parameter combinations converge. In coupled field sweeps near instability thresholds, 10-30% of runs may fail (solver timeout, energy divergence, CFL violation). Currently, failed runs produce empty metric rows with no explanation. Users need to know what failed and why.

### Specification

**Add to each sweep row:**

| Column             | Type        | Description                                                           |
| ------------------ | ----------- | --------------------------------------------------------------------- |
| `run_status`       | str         | "success", "timeout", "diverged", "solver_error", "measurement_error" |
| `error_message`    | str or None | Brief error description if status != "success"                        |
| `solver_exit_code` | int         | 0 for success, nonzero for failure                                    |

**Add to `SweepResults`:**

```python
@property
def failure_rate(self) -> float:
    """Fraction of runs that failed (status != 'success')."""

def successful_rows(self) -> SweepResults:
    """Return only successful runs."""

def failed_rows(self) -> SweepResults:
    """Return only failed runs (for debugging)."""
```

**Failure classification logic in `_run_single()`:**

```python
try:
    exit_code, wall_time = _simulate_run(...)
    if exit_code != 0:
        return {"run_status": "solver_error", "error_message": f"exit code {exit_code}", ...}
    metrics = _measure_run(...)
    return {"run_status": "success", **metrics}
except TimeoutError:
    return {"run_status": "timeout", "error_message": "simulation exceeded time limit", ...}
except Exception as exc:
    return {"run_status": "diverged", "error_message": str(exc)[:200], ...}
```

### Modified Files

| File                                  | Change                                                   |
| ------------------------------------- | -------------------------------------------------------- |
| `tidal/cli/_sweep.py`                 | Add status tracking to `_run_single()`, `_build_row()`   |
| `tidal/measurement/_sweep_results.py` | Add `failure_rate`, `successful_rows()`, `failed_rows()` |
| `tests/test_sweep.py`                 | Test failure tracking, classification                    |

### Acceptance Criteria

- [x] Successful runs have `run_status = "success"` in CSV
- [x] Failed runs have appropriate status and error message
- [x] `results.failure_rate` returns correct fraction
- [x] `results.successful_rows()` filters correctly
- [x] `results.failed_rows()` returns failed runs for debugging
- [x] Partial sweep results preserved (failed runs don't block saving)

---

## F8: Advanced Sweep Visualization

**Priority:** 8 (stretch)
**Generality:** [GENERAL] — operates on sweep result data
**Status:** Complete (Phase D, commit `6150ec6`)
**Depends on:** F5 (SweepResults query methods useful for data preparation)

### Motivation

For 3+ parameter sweeps, 1D/2D plots are insufficient. Standard scientific visualization tools offer parallel coordinates, tornado charts, and scatter matrices that help identify patterns in high-dimensional parameter spaces.

### Specification

**New plot types in `_sweep_panels.py`:**

| Plot Type        | Description                                                                                            | Use Case                                                                      |
| ---------------- | ------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------- |
| `sweep-parallel` | Parallel coordinates: each vertical axis = parameter, each polyline = run, color = metric              | Multi-parameter exploration: which parameter combinations produce high P_max? |
| `sweep-tornado`  | Horizontal bars: for each parameter, bar shows min-to-max metric range when only that parameter varies | Quick importance ranking: which parameters have the largest effect?           |
| `sweep-scatter`  | Pairwise scatter matrix: each parameter pair as 2D scatter, diagonal = 1D histogram, color = metric    | Interaction detection: do parameters combine nonlinearly?                     |

**CLI usage:**

```bash
tidal plot sweep_output/ --type sweep-parallel --metric P_max
tidal plot sweep_output/ --type sweep-tornado --metric P_max
tidal plot sweep_output/ --type sweep-scatter --metric P_max
```

### Modified Files

| File                         | Change                                    |
| ---------------------------- | ----------------------------------------- |
| `tidal/cli/_sweep_panels.py` | New render functions for each plot type   |
| `tidal/cli/_plot_command.py` | Register new plot types in `_SWEEP_TYPES` |
| `tests/test_sweep.py`        | Smoke tests for new plot types            |

### Acceptance Criteria

- [x] `--type sweep-parallel` produces readable parallel coordinates plot
- [x] `--type sweep-tornado` shows parameter importance ranking
- [x] `--type sweep-scatter` shows pairwise relationships
- [x] All three work with 2+ swept parameters
- [x] Color mapping is clear and includes colorbar

---

## Dependency Graph

```
F1 (TOML Config)
 |
 v
F2 (Adaptive Sampling)  -------> F4 (Sensitivity Analysis)
                                     |
                                     v (requires specific sampling)
F3 (Velocity/Resonance) - standalone
F5 (Query Methods)      - standalone --> F8 (Advanced Viz)
F6 (Spectrum Scalars)   - standalone
F7 (Failure Tracking)   - standalone
```

Most features are independent. The main dependency chain is F1 -> F2 -> F4 (TOML enables adaptive config, and Sobol analysis requires specific sampling designs).

---

## Implementation Phases

| Phase | Features      | Est. Scope | Branch                    | Status   |
| ----- | ------------- | ---------- | ------------------------- | -------- |
| A     | F1 + F2a      | Medium     | `feature/parameter-sweep` | Complete |
| B     | F2b + F5 + F7 | Medium     | `feature/parameter-sweep` | Complete |
| C     | F3 + F6       | Medium     | `feature/parameter-sweep` | Complete |
| D     | F4 + F8       | Medium     | `feature/parameter-sweep` | Complete |

All 8 features are implemented. 1242 tests passing, lint clean.
