# Lagrangian-to-PDE Simulation Pipeline

This document describes the implementation of the symbolic Lagrangian-to-PDE pipeline for GitHub Issue #33.

## Overview

The pipeline transforms symbolic Lagrangians (expressed in Mathematica/xAct) into numerical PDE simulations using py-pde. The key principle is that **NO physics is hardcoded** - all equation structure comes from symbolic derivation.

## Architecture

```
Mathematica/xAct          JSON              Python/py-pde
┌─────────────┐      ┌──────────┐      ┌───────────────┐
│ Lagrangian  │      │          │      │               │
│ L = ...     │ ───> │ Equation │ ───> │ PDEFromSpec   │
│             │      │ System   │      │               │
│ xTensor     │      │ (JSON)   │      │ FieldCollection
└─────────────┘      └──────────┘      └───────────────┘
     │                     │                    │
     │ Euler-Lagrange     │ Parse/Validate     │ Solve
     │ Linearize          │ Build PDE          │ py-pde
     │ Decompose          │                    │
     └────────────────────┴────────────────────┘
```

## Directory Structure

```
tidal/
├── symbolic/                # Python-side symbolic processing
│   ├── json_loader.py       # Load equations from JSON → EquationSystem
│   ├── pde_builder.py       # Build PDEBase from spec → PDEFromSpec
│   └── __init__.py
├── cli/                     # CLI (`tidal` command)
│   ├── __main__.py          # Entry point (derive/simulate/inspect/list/validate)
│   ├── _derive.py           # TOML → .wls → wolframscript → JSON
│   ├── _simulate.py         # JSON → PDE → solve → plot
│   ├── _inspect.py          # Display equation system info
│   ├── _list.py             # Discover available JSON specs
│   ├── _validate.py         # JSON spec validation
│   └── _plot.py             # Plotting utilities
├── vectorfield/             # Multi-component field utilities
│   ├── config.py            # ComponentFieldParams
│   ├── initial_conditions.py # Gaussian pulses, plane waves
│   └── __init__.py
└── wolfram/                 # Mathematica/xAct packages
    ├── EulerLagrange.wl     # Euler-Lagrange derivation
    ├── CommonUtilities.wl   # Shared utilities, epsilon evaluation
    ├── ComponentDecompose.wl # Decompose to components
    └── ExportJSON.wl        # JSON export, mass/coupling extraction
```

## Usage

### 1. Mathematica Side: Derive Equations from Lagrangian

The recommended workflow uses `tidal derive` with a TOML configuration file:

```bash
# Derive equations from a theory TOML file (generates .wls, runs wolframscript)
tidal derive examples/electromagnetic/theory.toml
```

Or use the Wolfram pipeline directly:

```mathematica
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];

(* Setup 1+1D Minkowski spacetime *)
DefManifold[M2, 2, IndexRange[a, h]];
DefChart[cart, M2, {0, 1}, {t[], x[]}];
DefMetric[-1, eta[-a, -b], CD, PrintAs -> "η"];

(* Define EM vector potential and derive Euler-Lagrange equations *)
DefTensor[A[-a], M2];
L = -1/2 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]];
eomA = EulerLagrangeEquation[L, A[a], CD];

(* Decompose to components and export *)
comps = DecomposeToComponents[eomA, A[-a], cart, {}];
jsonStructure = BuildMultiFieldJSONStructure[comps, metadata];
```

### 2. Python Side: Load and Simulate

```python
from pde import CartesianGrid
from tidal.symbolic import build_pde_from_json
from tidal.vectorfield import ComponentGaussianPulse

# Load equations derived from Lagrangian
pde = build_pde_from_json("examples/data/em_1d.json")

# Setup simulation
grid = CartesianGrid([(0, 100)], 256, periodic=True)

# Create initial conditions
pulse = ComponentGaussianPulse(
    center=(50.0,),
    width=5.0,
    active_components={"A_1": 1.0}  # Excite spatial component
)
initial = pulse.create(grid, pde.spec)

# Run simulation
result = pde.solve(initial, t_range=25.0, dt=0.01)
```

## JSON Schema

The JSON file is the contract between symbolic and numerical layers:

```json
{
  "metadata": {
    "lagrangian_expr": "-1/4 F[-a,-b] F[a,b]",
    "gauge": "lorenz",
    "linearized": true
  },
  "spacetime": {
    "dimension": 2,
    "signature": [-1, 1],
    "coordinates": ["t", "x"]
  },
  "fields": [
    {"name": "A_0", "index": 0, "is_dynamical": true},
    {"name": "A_1", "index": 1, "is_dynamical": true}
  ],
  "equations": [
    {
      "field": "A_0",
      "lhs": "d2_t(A_0)",
      "rhs": {
        "type": "linear_combination",
        "terms": [
          {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"}
        ]
      }
    }
  ],
  "coupling": {
    "mass_matrix": [[0, 0], [0, 0]],
    "coupling_matrix": [[0, 0], [0, 0]]
  }
}
```

## Phase 12: Auto-Computed Mass/Coupling Matrices

As of Phase 12 (February 2026), mass and coupling matrices are **automatically computed** from equation terms, eliminating the need for manual matrix specification.

### Matrix Convention

```
mass_matrix[i][j] = -(coefficient of identity(field_j) in equation_i)
coupling_matrix[i][j] = -(coefficient of identity(field_j) in equation_i, where i≠j)
```

### Symbolic Coefficient Preservation

The pipeline preserves both numeric and symbolic forms of coefficients:

```json
{
  "coefficient": -1.0,
  "coefficient_symbolic": "-m2",
  "operator": "identity",
  "field": "phi_0"
}
```

This enables **runtime parameter sweeps** without regenerating JSON files — symbolic coefficients are evaluated dynamically using the `_mathematica_to_python` expression evaluator.

### Auto-Computation Workflow

**Wolfram Side (`ExportJSON.wl`):**
- `ExtractMassCouplingFromEquations[fieldEquations]` parses RHS terms
- Identifies identity operator terms for each field
- Extracts both numeric and symbolic coefficients
- Exports `mass_matrix`, `coupling_matrix`, `mass_matrix_symbolic`, `coupling_matrix_symbolic`

**Python Side (`json_loader.py`):**
- `EquationSystem.from_dict()` calls `_compute_matrices_from_terms()`
- Returns 4-tuple: `(mass_numeric, coupling_numeric, mass_symbolic, coupling_symbolic)`
- Defense-in-depth: Both Wolfram and Python compute matrices to catch inconsistencies
- `__post_init__` guard: UserWarning if constructor-provided matrices differ from computed values

### Example: Coupled Scalars with Symbolic Coefficients

```mathematica
(* Wolfram: Define Lagrangian with symbolic parameters *)
L = 1/2 D[-a][phi[]] eta[a,b] D[-b][phi[]] - 1/2 m_phi^2 phi[]^2 +
    1/2 D[-a][chi[]] eta[a,b] D[-b][chi[]] - 1/2 m_chi^2 chi[]^2 -
    g phi[] chi[]
```

**Generated JSON:**
```json
{
  "equations": [
    {
      "field": "phi_0",
      "rhs": {
        "terms": [
          {"coefficient": -1.0, "coefficient_symbolic": "-m_phi2", "operator": "identity", "field": "phi_0"},
          {"coefficient": -0.5, "coefficient_symbolic": "-g", "operator": "identity", "field": "chi_0"},
          {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"}
        ]
      }
    }
  ],
  "coupling": {
    "mass_matrix": [[1.0, 0.0], [0.0, 4.0]],
    "mass_matrix_symbolic": [["m_phi2", "0"], ["0", "m_chi2"]],
    "coupling_matrix": [[0.0, 0.5], [0.5, 0.0]],
    "coupling_matrix_symbolic": [["0", "g"], ["g", "0"]]
  }
}
```

**Python Simulation:**
```python
# Sweep coupling parameter without regenerating JSON
for g_value in [0.1, 0.5, 1.0, 2.0]:
    pde = build_pde_from_json("coupled_scalars.json", parameters={"g": g_value, "m_phi2": 1.0, "m_chi2": 4.0})
    result = pde.solve(initial_state, t_range=10.0)
```

## Examples

### Klein-Gordon (Validation)

The Klein-Gordon field validates the full pipeline:

1. Mathematica derives: `Box[phi] - m^2 phi = 0`
2. Export to JSON
3. Python builds `PDEFromSpec`
4. Simulate and verify dispersion relation matches theory
5. Energy conservation holds ✓

### Electromagnetic Field (1+1D)

From `L = -1/4 F_μν F^μν`:

1. Derive Maxwell equations in Lorenz gauge
2. Get wave equations: `∂²A_μ/∂t² = ∂²A_μ/∂x²`
3. Simulate massless wave propagation at c=1

## Key Classes

### `PDEFromSpec`

Dynamically constructs PDE from equation specification:

```python
class PDEFromSpec(PDEBase):
    def __init__(self, spec: EquationSystem, parameters: dict | None = None):
        self.spec = spec  # All physics comes from here
        # parameters override symbolic coefficients at runtime

    def evolution_rate(self, state, t=0.0):
        # Build RHS from spec.equations
        for eq in self.spec.equations:
            for term in eq.rhs_terms:
                # Apply operator specified in JSON
                ...
```

### `EquationSystem`

Parsed equation specification:

```python
@dataclass(frozen=True)
class EquationSystem:
    n_components: int
    component_names: tuple[str, ...]
    equations: tuple[ComponentEquation, ...]
    mass_matrix: tuple[tuple[float, ...], ...]
    metadata: dict
```

### `ComponentGaussianPulse`

Initial conditions for multi-component fields:

```python
pulse = ComponentGaussianPulse(
    center=(50.0,),
    width=5.0,
    active_components={"A_1": 1.0}
)
state = pulse.create(grid, spec)
```

## Testing

### Unit Tests

- `tests/test_json_loader.py` - JSON parsing and validation
- `tests/test_pde_builder.py` - PDE construction from spec (includes EM pipeline tests)
- `tests/test_vectorfield.py` - Multi-component utilities
- `tests/test_measurement.py` - Energy, conversion, spectral analysis
- `tests/test_cli.py` - CLI subcommands and TOML-to-WLS generation

### Validation

```bash
# Run full test suite (Python + Wolfram)
./scripts/full_test.sh

# Validate end-to-end pipeline
./scripts/validate_pipeline.sh
```

Checks:
- ✓ All modules import correctly
- ✓ JSON files load and validate
- ✓ PDEs build from specifications
- ✓ Simulations run and produce finite results
- ✓ Physics is correct (wave speed, energy conservation)

## Design Principles

1. **No Hardcoded Physics**: All equation structure from JSON
2. **Symbolic First**: Equations derived, not guessed
3. **Layered Architecture**: Clean separation Mathematica ↔ Python
4. **Type Safety**: Frozen dataclasses, type annotations
5. **Testable**: Each layer validates independently
6. **Extensible**: Same pipeline for scalar, vector, tensor fields

## Completed Extensions

- ✅ **Coupled systems**: Multi-field Lagrangians (coupled_scalars, scalar_vector_coupling)
- ✅ **Higher dimensions**: 2+1D and 3+1D spacetimes (20 examples)
- ✅ **Gauge constraints**: Constraint equations with `time_derivative_order=0`
- ✅ **Rank 3+ tensors**: Antisymmetric rank-3 with symmetry reduction (massive_3form)
- ✅ **Curvilinear coordinates**: Polar, spherical, cylindrical with Christoffel auto-detection
- ✅ **CLI**: `tidal` command with derive/simulate/measure/inspect/list/validate subcommands
- ✅ **Background fields**: Non-dynamical `[[background_fields]]` for position-dependent coupling
- ✅ **Measurement module**: Energy, conversion P(t), spectral P(k,t), dispersion, mixing length
- ✅ **Auto-computed matrices**: Mass/coupling matrices from identity terms (Phase 12)

## Future Extensions

- **Nonlinear theories**: Yang-Mills, Einstein-Hilbert (beyond linear perturbation)
- **Automatic gauge fixing**: Lorenz/Coulomb in Wolfram layer
- **Code generation**: Generate optimized Numba kernels
- **CI**: GitHub Actions for Wolfram tests on PRs

## References

- **xAct**: https://xact.es/
- **py-pde**: https://py-pde.readthedocs.io/
- **GitHub Issue #33**: https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/33
