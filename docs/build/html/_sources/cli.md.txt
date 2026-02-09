# CLI Reference

TIDAL provides the `tidal` command-line tool with five subcommands for the full derive-to-simulate workflow.
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

# Generate and run wolframscript
tidal derive examples/scalar_field/theory.toml --run
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
    --grid-size 256 \
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
| `--grid-size INT` | Grid points per axis | 128 |
| `--ic {gaussian,plane-wave,zero,formula}` | Initial condition preset | gaussian |
| `--ic-width FLOAT` | Gaussian pulse width | 5.0 |
| `--ic-formula EXPR` | Custom IC formula (Python expression) | — |
| `--bc TYPES` | Boundary conditions (comma-separated per axis) | periodic |
| `--scheme {scipy,runge-kutta}` | Solver scheme | runge-kutta |
| `--mode {evolution,constraint}` | Simulation mode | evolution |
| `--plot` / `--no-plot` | Enable/disable plotting | --plot |
| `--output PATH` | Output file path | — |

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

### Derived Fields

For theories that need intermediate tensor definitions (e.g., the electromagnetic field strength tensor), use `[[derived_fields]]`:

```toml
[[derived_fields]]
name = "F"
rank = 2
symmetry = "Antisymmetric"
rule = "F[-a, -b] -> CD[-a][A[-b]] - CD[-b][A[-a]]"
```

This generates `DefTensor` + `MakeRule` in the `.wls` output, with automatic substitution before decomposition.

## Common Workflows

### Derive and simulate a scalar field theory

```bash
# 1. Write theory.toml (see examples/scalar_field/theory.toml)

# 2. Derive equations
tidal derive examples/scalar_field/theory.toml --run

# 3. Inspect the result
tidal inspect examples/data/klein_gordon_1d.json

# 4. Simulate
tidal simulate examples/data/klein_gordon_1d.json --t-end 20 --ic gaussian
```

### Parameter sweep without re-deriving

```bash
# Same JSON, different parameters
for m2 in 0.5 1.0 2.0 4.0; do
    tidal simulate examples/data/coupled_scalars.json \
        --param m_phi2=$m2 --output "sweep_m2_${m2}.png"
done
```

### Solver schemes

- `--scheme runge-kutta`: py-pde `ExplicitSolver` — fixed-step, predictable timing
- `--scheme scipy`: py-pde `ScipySolver` — adaptive step, better for stiff systems

### Boundary conditions

Specify per-axis, comma-separated:

```bash
# Periodic in x, Neumann in y
tidal simulate spec.json --bc periodic,neumann
```

Available types: `periodic`, `neumann`, `dirichlet`.
