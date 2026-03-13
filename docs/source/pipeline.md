# Symbolic Pipeline

TIDAL implements a two-stage pipeline that transforms symbolic Lagrangians into numerical PDE simulations.
**No physics is hardcoded** — all equation structure comes from symbolic derivation.

## Data Flow

```
Mathematica/xAct              JSON                Python/SUNDIALS
┌──────────────────┐    ┌──────────────┐    ┌──────────────────┐
│ Lagrangian L(φ)  │    │   Equation   │    │ EquationSystem   │
│                  │ -> │   System     │ -> │                  │
│ Euler-Lagrange   │    │   (JSON)     │    │ FieldSet + RHS   │
│ Decompose        │    │              │    │ IDA/CVODE/LF     │
└──────────────────┘    └──────────────┘    └──────────────────┘
```

**Stage 1 (Mathematica/xAct)** derives field equations symbolically from a Lagrangian.
**Stage 2 (Python/SUNDIALS + numpy)** loads the specification and runs a numerical simulation.

## Wolfram Modules

| Module | Purpose | Key Function |
|--------|---------|-------------|
| `CommonUtilities.wl` | Shared helpers: CD-to-Derivative conversion, Christoffel/epsilon evaluation | `ConvertCDToDerivatives`, `EvaluateChristoffelComponents` |
| `EulerLagrange.wl` | Derive equations of motion from a Lagrangian | `EulerLagrangeEquation[L, field, cd]` |
| `ComponentDecompose.wl` | Convert tensor EOM to scalar component equations | `DecomposeToComponents[eom, field, chart, additionalFields]` |
| `ExportJSON.wl` | Serialize component equations to JSON | `BuildMultiFieldJSONStructure[fieldEquations, metadata]` |
| `Linearize.wl` | Perturbation theory (xPert) for linearized gravity | `LinearizeTensorExpression[expr]` |

## JSON Specification

The JSON file is the contract between symbolic and numerical layers:

```json
{
  "spacetime": { "dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"] },
  "fields": [
    { "name": "phi_0", "index": 0, "is_dynamical": true }
  ],
  "equations": [
    {
      "field": "phi_0",
      "lhs": { "expression": "d2_t(phi_0)", "order": { "time": 2 } },
      "rhs": {
        "type": "linear_combination",
        "terms": [
          { "coefficient": 1.0, "operator": "laplacian", "field": "phi_0" },
          { "coefficient": -1.0, "coefficient_symbolic": "-m2", "operator": "identity", "field": "phi_0" }
        ]
      }
    }
  ],
  "coupling": {
    "mass_matrix": [[1.0]],
    "coupling_matrix": [[0.0]]
  },
  "metadata": { "lagrangian_expr": "1/2 (d phi)^2 - 1/2 m^2 phi^2" }
}
```

For the full schema reference, see the [JSON Schema Guide](https://github.com/WilliamRoyce/torsion-gertsenshtein/blob/main/docs/JSON_SCHEMA_GUIDE.md).

## Python-Side Classes

### `EquationSystem` (json_loader.py)

Parsed equation specification loaded from JSON:

```python
from tidal.symbolic.json_loader import load_equation_system

# Load equation specification
spec = load_equation_system("examples/data/klein_gordon_1d.json")

# Runtime parameter overrides are passed to the solver via --param or programmatically
```

Key features:
- Frozen dataclass for immutability
- Auto-computes mass/coupling matrices from equation terms
- Preserves symbolic coefficient expressions for runtime parameter sweeps

### Solver Classes (tidal/solver/)

The solver layer builds the numerical system from the equation specification:

- `FieldSet.from_spec(spec)` — field storage and boundary conditions
- `CoefficientEvaluator(spec)` — 4-level cached coefficient resolution
- `RHSEvaluator(spec)` — evaluates RHS from JSON terms (no hardcoded physics)
- `StateLayout.from_spec(spec)` — maps fields to state vector slots

Supports second-order (wave), first-order, and constraint (order-0) equations in the same system.

## Supported Operators

Operators map JSON terms to numpy spatial operations (finite-difference stencils). All support cross-field references.

| Operator | Min Dim | Description |
|----------|---------|-------------|
| `identity` | 1D | Field value (mass/coupling terms) |
| `laplacian` | 1D | Full spatial Laplacian |
| `laplacian_x`, `_y`, `_z` | 1D/2D/3D | Directional second derivative |
| `gradient_x`, `_y`, `_z` | 1D/2D/3D | Directional first derivative |
| `cross_derivative_xy`, `_xz`, `_yz` | 2D/3D | Mixed partial derivative |
| `first_derivative_t` | 1D | Time derivative (friction/damping terms) |
| `biharmonic` | 1D | Fourth-order Laplacian |

## Christoffel Auto-Detection

For curvilinear coordinates, the pipeline automatically detects whether Christoffel symbols are needed:

- **Constant metric** (flat Minkowski, static conformal): all Christoffels = 0
- **Non-constant metric** (polar, spherical, expanding universe): computed from standard formula

This is controlled by `DecomposeToComponents` with the `"MetricMatrix"` option. Override with `"ComputeChristoffels" -> True/False` for explicit control.

## Mass/Coupling Matrix Auto-Computation

As of Phase 12, mass and coupling matrices are **automatically computed** from equation terms on both the Wolfram and Python sides:

- **Convention**: `matrix[i][j] = -(coefficient of identity(field_j) in equation_i)`
- **Symbolic preservation**: `mass_matrix_symbolic` / `coupling_matrix_symbolic` preserve exact Mathematica expressions
- **Runtime sweeps**: Override symbolic coefficients via `parameters={"m2": 1.0, "g": 0.5}` without regenerating JSON
