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
torsion_gertsenshtein/
├── symbolic/                # Python-side symbolic processing
│   ├── json_loader.py       # Load equations from JSON
│   ├── pde_builder.py       # Build PDEBase from spec
│   └── __init__.py
├── vectorfield/             # Multi-component field utilities
│   ├── config.py            # ComponentFieldParams
│   ├── initial_conditions.py # Gaussian pulses, plane waves
│   └── __init__.py
└── wolfram/                 # Mathematica/xAct packages
    ├── EulerLagrange.wl     # Euler-Lagrange derivation
    ├── Linearize.wl         # Linearization
    ├── ComponentDecompose.wl # Decompose to components
    ├── ExportJSON.wl        # JSON export
    ├── LagrangianPipeline.wl # Main entry point
    └── examples/
        ├── klein_gordon.wls # KG validation
        └── em_lagrangian_1d.wls # EM example
```

## Usage

### 1. Mathematica Side: Derive Equations from Lagrangian

```mathematica
<< TorsionGertsenshtein`LagrangianPipeline`;

(* Setup 1+1D Minkowski spacetime *)
{M, eta, CD, cart} = SetupMinkowski1D[];

(* Define EM vector potential *)
DefTensor[A[-a], M];
DefTensor[F[-a, -b], M, Antisymmetric[{-a, -b}]];
FieldStrengthRule = F[-a, -b] -> CD[-a][A[-b]] - CD[-b][A[-a]];

(* EM Lagrangian: L = -1/4 F_ab F^ab *)
EMLagrangian = -1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d];
EMLagrangianExpanded = EMLagrangian /. FieldStrengthRule;

(* Process and export *)
json = ProcessLagrangian[
    EMLagrangianExpanded, A[-a], CD, cart,
    "OutputPath" -> "examples/data/em_1d.json"
];
```

### 2. Python Side: Load and Simulate

```python
from pde import CartesianGrid
from torsion_gertsenshtein.symbolic import build_pde_from_json
from torsion_gertsenshtein.vectorfield import ComponentGaussianPulse

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

## Examples

### Klein-Gordon (Validation)

The Klein-Gordon field provides cross-validation:

1. Mathematica derives: `Box[phi] - m^2 phi = 0`
2. Export to JSON
3. Python builds `PDEFromSpec`
4. Compare to existing `KleinGordonPDE`
5. Results must match ✓

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
    def __init__(self, spec: EquationSystem):
        self.spec = spec  # All physics comes from here

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
- `tests/test_pde_builder.py` - PDE construction from spec
- `tests/test_vectorfield.py` - Multi-component utilities
- `tests/test_em_pipeline.py` - End-to-end integration

### Validation

```bash
python validate_implementation.py
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

## Future Extensions

- **Coupled systems**: Multi-field Lagrangians
- **Nonlinear theories**: Yang-Mills, Einstein-Hilbert
- **Higher dimensions**: 2+1D, 3+1D spacetimes
- **Gauge constraints**: Enforce constraints dynamically
- **Code generation**: Generate optimized Numba kernels
- **Automated pipeline**: `wolframclient` bridge for full automation

## References

- **xAct**: https://xact.es/
- **py-pde**: https://py-pde.readthedocs.io/
- **GitHub Issue #33**: https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/33
