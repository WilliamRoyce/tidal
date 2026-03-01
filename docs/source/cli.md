# CLI Reference

TIDAL provides the `tidal` command-line tool with seven subcommands for the full derive-to-measure workflow.
**Zero additional dependencies** — uses only `argparse` and `tomllib` from the Python standard library.

## Installation

```bash
# Install with uv (recommended)
uv sync --all-extras

# Or via pip
pip install -e ".[dev]"

# Verify
tidal --version
```

## Subcommands

### `tidal derive` — Derive equations from a Lagrangian

Generates a Wolfram Language script (`.wls`) from a TOML configuration file and optionally runs it with `wolframscript`.

```bash
# Generate .wls only
tidal derive examples/scalar_field/theory.toml

# Generate .wls and run wolframscript via the generated run.sh
bash examples/scalar_field/run.sh
```

### `tidal simulate` — Run a PDE simulation

Loads a JSON equation specification and runs a numerical simulation with configurable parameters.

```bash
# Basic simulation with defaults
tidal simulate examples/data/klein_gordon_1d.json

# With parameter overrides and custom settings
tidal simulate examples/data/coupled_scalars.json \
    --param m_phi2=1.0 --param g=0.5 \
    --t-end 20.0 --dt 0.005 \
    --grid-shape 256 \
    --ic gaussian --ic-width 5.0 \
    --bc periodic \
    --scheme scipy

# Constraint-only solve (no time evolution)
tidal simulate examples/data/chern_simons.json --mode constraint
```

**Options:**

| Flag | Description | Default |
|------|-------------|---------|
| `--param KEY=VALUE` | Override symbolic coefficient | (from JSON) |
| `--t-end FLOAT` | Simulation end time | 10.0 |
| `--dt FLOAT` | Time step | 0.01 |
| `--grid-shape INT` | Grid points per axis | 128 |
| `--ic {gaussian,plane-wave,zero,formula,file,noise}` | Initial condition preset | gaussian |
| `--ic-width FLOAT` | Gaussian pulse width | 5.0 |
| `--ic-component FIELD` | Field component for IC (multi-field systems) | (first field) |
| `--ic-center X[,Y,Z]` | Gaussian center position | (domain center) |
| `--ic-amplitude FLOAT` | Gaussian amplitude | 1.0 |
| `--ic-wavevector K[,K,K]` | Wavevector for travelling wave packets | — |
| `--ic-formula EXPR` | Custom IC formula (Python expression with x, y, z, np) | — |
| `--ic-formula-velocity EXPR` | Velocity formula for `--ic formula` (enables travelling waves) | — |
| `--ic-field FIELD:EXPR` | Per-field IC override formula (repeatable; `FIELD:velocity:EXPR` for velocity) | — |
| `--ic-file PATH` | Load IC from `.npy` file or simulation output directory | — |
| `--ic-noise-seed N` | Random seed for reproducible `--ic noise` | — |
| `--bc TYPES` | Boundary conditions (comma-separated per axis) | periodic |
| `--scheme {ida,leapfrog,cvode,scipy,auto}` | Solver scheme | auto |
| `--mode {evolve,constraint}` | Simulation mode | evolve |
| `--plot` / `--no-plot` | Enable/disable plotting | --plot |
| `--output PATH` | Output path (directory for snapshot data; image extension for plot-only) | — |

### `tidal inspect` — Display equation system info

Shows a summary of a JSON equation specification without running a simulation.

```bash
tidal inspect examples/data/chern_simons.json
```

Output includes: spacetime dimension, field names, equation count, operator types, mass/coupling matrices.

### `tidal list` — Discover available JSON specs

Searches for JSON specification files in `examples/data/`.

```bash
tidal list
```

### `tidal validate` — Validate a JSON specification

Checks that a JSON file conforms to the expected equation system schema.

```bash
tidal validate examples/data/klein_gordon_1d.json

# JSON output for scripting
tidal validate examples/data/klein_gordon_1d.json --json
```

### `tidal measure` — Extract physics measurements

Loads a snapshot directory produced by `tidal simulate --output <dir>` and runs measurement analyses:
energy decomposition, conversion probability, mixing length, spectral analysis, and conservation diagnostics.

```bash
# Full summary (energy + conservation + auto-detect conversion + mixing)
tidal measure result_dir/ --spec spec.json

# Specific measurements
tidal measure result_dir/ --what conversion,mixing \
    --source phi_0 --target chi_0

# JSON output for scripting
tidal measure result_dir/ --what energy,conservation --json

# Save 2x3 measurement plot
tidal measure result_dir/ --output measurements.png
```

The JSON spec can be auto-discovered from `metadata.json` in the snapshot directory (stored by `tidal simulate`) or provided explicitly via `--spec`.

**Options:**

| Flag | Description | Default |
|------|-------------|---------|
| `--spec PATH` | JSON equation specification | (auto-discovered from metadata.json) |
| `--what TYPES` | Measurement types (comma-separated) | summary |
| `--source FIELDS` | Source field(s) for conversion (comma-separated) | (auto-detect) |
| `--target FIELDS` | Target field(s) for conversion (comma-separated) | (auto-detect) |
| `--param KEY=VALUE` | Override parameter from metadata.json | (from metadata.json) |
| `--energy-threshold T` | Conservation threshold | 1e-3 |
| `--output PATH` | Save plot (.png/.pdf) | — |
| `--json` | JSON output instead of text | — |
| `--quiet` | Suppress progress messages | — |

**Available measurement types:** `summary`, `energy`, `conversion`, `mixing`, `spectrum`, `spectral_conversion`, `dispersion`, `conservation`, `effective_mass`, `asymptotic`, `peak_conversion`.

### `tidal sweep` — Run parameter sweeps and convergence studies

Automates running `simulate` + `measure` across a parameter grid, collecting scalar metrics into a portable CSV/JSON results table.

```bash
# 1D sweep: coupling strength
tidal sweep examples/data/coupled_scattering.json \
  --sweep "g0=0.1:0.9:10" \
  --measure conversion,mixing --source phi_0 --target chi_0 \
  --grid-shape 256 --bounds=-40:40 --periodic \
  --ic gaussian --ic-component phi_0 --ic-center=-20 --ic-width 3 --ic-wavevector 3 \
  --t-end 60 --output sweeps/coupling_scan

# 2D sweep: coupling x mass (cartesian product)
tidal sweep examples/data/coupled_scattering.json \
  --sweep "g0=0.1:0.9:5" --sweep "mChi2=0.5:4.0:5" \
  --measure conversion --source phi_0 --target chi_0 \
  --output sweeps/mass_coupling_2d

# Grid convergence study
tidal sweep examples/data/coupled_scattering.json \
  --converge "32,64,128,256,512" \
  --measure conservation \
  --output sweeps/convergence
```

**Sweep specification syntax:**

| Format | Description | Example |
| ------ | ----------- | ------- |
| `PARAM=START:STOP:N` | N linearly-spaced values | `g0=0.1:1.0:10` |
| `PARAM=START:STOP:N:log` | N log-spaced values | `m2=0.01:100:20:log` |
| `PARAM=V1,V2,...` | Explicit values | `g0=0.1,0.5,0.9` |

Multiple `--sweep` flags produce a cartesian product grid.

**Output structure:**

```text
sweeps/coupling_scan/
├── sweep.json          # Sweep definition + provenance
├── results.csv         # Self-contained metrics table
├── results.json        # Same data in JSON
├── g0_0.100/           # Individual run output
│   ├── metadata.json
│   ├── times.npy
│   ├── phi_0.npy, v_phi_0.npy, ...
├── g0_0.200/
│   └── ...
```

The CSV is **self-contained and portable**: every row includes all parameters (swept + fixed), simulation settings, and measured metrics so the data can be analyzed outside TIDAL.

**Options:**

| Flag | Description | Default |
|------|-------------|---------|
| `--sweep SPEC` | Parameter sweep specification (repeatable) | — |
| `--converge SIZES` | Grid sizes for convergence study | — |
| `--measure TYPES` | Measurement types (comma-separated) | summary |
| `--source FIELDS` | Source field(s) for conversion measurements | — |
| `--target FIELDS` | Target field(s) for conversion measurements | — |
| `--resume` | Skip completed runs (checks for metadata.json) | — |
| `--parallel N` | Number of parallel workers | 1 |

All `tidal simulate` flags (`--param`, `--grid-shape`, `--bounds`, `--ic`, etc.) are passed through unchanged.

**Visualizing sweep results** — use `tidal plot`:

```bash
# 1D sweep → line plot
tidal plot sweeps/coupling_scan/ --type sweep --metric P_max

# 2D sweep → heatmap
tidal plot sweeps/mass_coupling_2d/ --type sweep --metric P_max

# Overlay timeseries from each run
tidal plot sweeps/coupling_scan/ --type sweep-compare --metric conversion

# Convergence → log-log with fitted order
tidal plot sweeps/convergence/ --type convergence --metric max_energy_error
```

## TOML Configuration (`theory.toml`)

The `tidal derive` command reads a TOML file that defines the field theory:

```toml
[spacetime]
dimension = 3          # 2 = 1+1D, 3 = 2+1D, 4 = 3+1D
signature = [-1, 1, 1]
coordinates = ["t", "x", "y"]

[fields]
names = ["phi"]
ranks = [0]            # 0 = scalar, 1 = vector, 2 = tensor

[constants]
names = ["m2"]
values = [1.0]

[lagrangian]
expression = "1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - 1/2 m2 phi[]^2"

[parameters]
output_json = "examples/data/my_theory.json"

# Optional: diagonal metric for curvilinear coordinates
[metric]
diagonal = [-1, 1, "x[]^2"]
```

### Linearization (`[linearization]`)

For theories defined by linearizing a tensor expression (e.g., linearized gravity via xPert), use `[linearization]` **instead of** `[lagrangian]`. The two sections are mutually exclusive.

```toml
[linearization]
# Tensor expression to linearize (uses CD placeholder, auto-prefixed)
expression = "Einstein[CD][-a, -b]"
# Which [[fields]] entry is the metric perturbation
perturbation_field = "h"
```

This generates xPert setup (`SetupMetricPerturbation`), linearization (`LinearizeTensorExpression`), and automatic notation conversion before passing to the standard decomposition pipeline.

**Example** (`examples/gravitational_waves/theory.toml`):

```toml
[theory]
name = "Linearized Gravity"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[linearization]
expression = "Einstein[CD][-a, -b]"
perturbation_field = "h"

[output]
path = "../data/linearized_gravity.json"
```

### Derived Fields

For theories that need intermediate tensor definitions (e.g., the electromagnetic field strength tensor), use `[[derived_fields]]`:

```toml
[[derived_fields]]
name = "F"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"
```

This generates `DefTensor` + `MakeRule` in the `.wls` output, with automatic substitution before decomposition.

## Common Workflows

### Derive and simulate a scalar field theory

```bash
# 1. Write theory.toml (see examples/scalar_field/theory.toml)

# 2. Derive equations
tidal derive examples/scalar_field/theory.toml

# 3. Inspect the result
tidal inspect examples/data/klein_gordon_1d.json

# 4. Simulate
tidal simulate examples/data/klein_gordon_1d.json --t-end 20 --ic gaussian
```

### Parameter sweep

```bash
# Automated sweep (recommended)
tidal sweep examples/data/coupled_scalars.json \
    --sweep "mPhi2=0.5:4.0:8" \
    --measure conversion --source phi_0 --target chi_0 \
    --output sweeps/mass_scan

# Or manually with a shell loop
for m2 in 0.5 1.0 2.0 4.0; do
    tidal simulate examples/data/coupled_scalars.json \
        --param mPhi2=$m2 --output "sweep_m2_${m2}.png"
done
```

### Solver schemes

- `--scheme auto` (default): auto-selects best solver based on equation structure
- `--scheme cvode`: SUNDIALS/CVODE — adaptive BDF/Adams for wave systems
- `--scheme ida`: SUNDIALS/IDA — handles all equation types including algebraic constraints
- `--scheme leapfrog`: symplectic Störmer-Verlet — fixed dt, zero energy drift
- `--scheme scipy`: scipy.integrate.solve_ivp — DOP853/RK45/Radau/BDF

### Boundary conditions

Specify per-axis, comma-separated:

```bash
# Periodic in x, Neumann in y
tidal simulate spec.json --bc periodic,neumann
```

Available types: `periodic`, `neumann`, `dirichlet`.

### Full derive-to-measure pipeline

```bash
# 1. Derive equations from a Lagrangian
tidal derive examples/coupled_scalars/theory.toml

# 2. Simulate and save to snapshot directory
tidal simulate examples/data/coupled_scalars.json \
    --param mPhi2=1.0 --param mChi2=4.0 --param gCpl=0.5 \
    --t-end 20.0 --output coupled_scalars_output

# 3. Extract measurements: conversion probability and mixing length
tidal measure coupled_scalars_output/ \
    --what conversion,mixing \
    --source phi_0 --target chi_0

# 4. Save measurement plot
tidal measure coupled_scalars_output/ --output measurements.png
```
