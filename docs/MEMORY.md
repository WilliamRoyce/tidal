# Project Memory: TIDAL Lagrangian-to-PDE Pipeline

## Document Maintenance

**IMPORTANT:** This file should be updated throughout the project lifecycle to capture learnings and patterns.

**When to update this file:**

- After implementing a new physics example (scalar, vector, tensor fields)
- When discovering new xAct/Wolfram patterns or pitfalls
- After generalizing infrastructure (dimension support, operator handling)
- When solving tricky bugs that reveal architectural insights
- After making critical design decisions that affect future work

**Companion files:**

- `troubleshooting.md` - Error patterns and debugging techniques
- `chern-simons-notes.md` - Example-specific implementation details

---

## Overview

This project implements a symbolic physics pipeline: Lagrangian (xAct/Mathematica) → JSON → PDE simulation (py-pde). Users should NEVER manually hardcode equations - all PDEs must derive from Lagrangian via symbolic computation.

## Critical Architecture Decisions

### Multi-Field Coupling (WORKING)

- **Cross-field terms:** Use gradient operators in JSON (e.g., `gradient_x(A_2)` in A_0 equation)
- **Field transformation:** `DecomposeToComponents` accepts `additionalFields` parameter to transform all coupled fields to coordinate form
- **JSON detection:** `IdentifyMultiFieldTerm` uses case-insensitive symbol matching and function head extraction

### Dimension Handling (GENERALIZED)

- **Dynamic dimension:** `ComponentDecompose.wl` now uses `dim = Length[ScalarsOfChart[chart]]` instead of hardcoded `dim = 2`
- **Derivative conversion:** `CommonUtilities.wl` handles both 2-arg (1+1D) and 3-arg (2+1D) Derivative forms
- **Supports:** 1+1D (t,x), 2+1D (t,x,y), 3+1D (t,x,y,z), extensible to higher

### Mass/Coupling Matrix Auto-Computation (Phase 12)

- **Convention:** `matrix[i][j] = -(coefficient of identity(field_j) in equation_i)`
- **Wolfram:** `ExtractMassCouplingFromEquations` in ExportJSON.wl auto-computes at export
- **Python:** `_compute_matrices_from_terms` in json_loader.py auto-computes at load (defense-in-depth)
- **Symbolic preservation:** `mass_matrix_symbolic` / `coupling_matrix_symbolic` preserve exact Mathematica expressions; only evaluated at runtime

### Christoffel Symbol Computation (Auto-Detection)

- **Auto-detection (default):** `DecomposeToComponents[eom, field, chart, {}, "MetricMatrix" -> matrix]`
- **Constant metric** (flat or static conformal) → All Christoffels = 0
- **Non-constant metric** (time/position-dependent) → Compute from standard formula
- **Override available:** `"ComputeChristoffels" -> True/False` for explicit control

## xAct/Wolfram Patterns

### Symbol Management

- **Always check existence:** `If[!xTensorQ[M2], DefManifold[M2, ...]]`
- **Use shared symbols:** M2, M3, eta, CD for examples to avoid conflicts

### Epsilon Tensors (AUTOMATED)

- `DefMetric` automatically creates `epsiloneta3[-a,-b,-c]`
- `EvaluateEpsilonComponents` in CommonUtilities.wl → numeric ±1
- Minkowski (-,+,+): covariant ε_012 = -1, contravariant ε^012 = +1

### xAct Cycles Context (CRITICAL)

- xAct's `Cycles` is `xAct`xPerm`Cycles`, NOT `System`Cycles`
- All pattern matching must use `_xAct`xPerm`Cycles` and `Head[x] === xAct`xPerm`Cycles`

### Lagrangian Construction Rules

1. Use **explicit metric tensors** (`eta[a,b]`) for index raising
2. Use **`DefConstantSymbol[m2]`** for mass/coupling constants (not bare `Symbol`)
3. **Parenthesize multi-line expressions** in .wls: `L = (term1 + term2);`
4. **Expand field strength before decompose:** `L = -1/2 CD[-a][A[-b]] ...` (not via F substitution)

## CLI (`tidal` Command)

The `tidal` CLI provides 5 subcommands with zero new dependencies:

| Command                    | Description                                |
| -------------------------- | ------------------------------------------ |
| `tidal derive theory.toml` | Generate .wls from TOML, run wolframscript |
| `tidal simulate spec.json` | Full simulation with smart defaults        |
| `tidal inspect spec.json`  | Display equation system info               |
| `tidal list`               | Discover available JSON specs              |
| `tidal validate spec.json` | Validate JSON spec structure               |

**TOML Configuration:**

- `theory.toml` with spacetime, fields, constants, Lagrangian, parameters
- `[[derived_fields]]` for intermediate tensors (e.g., field strength F_ab)
- IC presets: `gaussian`, `plane-wave`, `zero`, `formula`
- Per-axis boundary conditions via `--bc neumann,periodic`

## Pipeline Module Integration

### Loading Pattern

```mathematica
Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];
Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];
Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];
Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];
```

### Cross-Field Decomposition

```mathematica
phiComponents = DecomposeToComponents[eomPhi, phi[], cart, {chi[]}];
chiComponents = DecomposeToComponents[eomChi, chi[], cart, {phi[]}];
```

### JSON Export

```mathematica
fieldEquations = {{"phi_0", eqPhi}, {"chi_0", eqChi}};
jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];
```

## Python Operators

`identity`, `laplacian`, `laplacian_{x,y,z}`, `gradient_{x,y,z}`, `cross_derivative_{xy,xz,yz}`, `first_derivative_t`, `biharmonic`

- All support cross-field application and momentum (`pi_i`) references

## Example Implementations

| Example                   | Dim  | Key Features                                                               |
| ------------------------- | ---- | -------------------------------------------------------------------------- |
| `scalar_field/`           | 1+1D | KG, mass term, dispersion                                                  |
| `electromagnetic/`        | 1+1D | Maxwell, Lorenz gauge                                                      |
| `proca/`                  | 1+1D | Massive vector (Proca mass)                                                |
| `coupled_scalars/`        | 1+1D | Cross-field coupling, mass matrix                                          |
| `chern_simons/`           | 2+1D | Epsilon tensor, A_0 constraint                                             |
| `elasticity/`             | 2+1D | Anisotropic laplacian, cross_derivative_xy                                 |
| `curved_spacetime/`       | 2+1D | Hubble friction, time-dependent coefficients                               |
| `sphere_kg/`              | 2+1D | Position-dependent coefficients, S²                                        |
| `polar_kg/`               | 2+1D | Polar coordinates, Christoffel auto-detection                              |
| `scalar_vector_coupling/` | 2+1D | Mixed-rank cross-field (scalar+vector), 4 constants                        |
| `scalar_field_3d/`        | 3+1D | Full 4D KG                                                                 |
| `spherical_kg/`           | 3+1D | Spherical coordinates, trig coefficients                                   |
| `cylindrical_kg/`         | 3+1D | Cylindrical, mixed curved/flat                                             |
| `gravitational_waves/`    | 3+1D | xPert linearization, TT gauge, constraints                                 |
| `massive_3form/`          | 3+1D | Rank-3 antisymmetric, symmetry reduction                                   |
| `electrostatics/`         | 2+1D | Poisson constraint, no time evolution                                      |
| `massive_gravity/`        | 2+1D | Linearized massive gravity, xPert, coupled constraints, SVD regularization |
| `coupled_proca/`          | 2+1D | Two massive vectors, coupled Helmholtz, Gauss-Seidel, Dirichlet BCs        |

## Testing Guidelines

### Test Counts

- **850 Python tests** + **~108 Wolfram tests** passing
- Run: `uv run pytest tests/` and `./scripts/full_test.sh`

### Verification Pattern

1. Wolfram derivation produces N component equations
2. JSON has correct dimension, signature, coordinates
3. JSON terms show proper cross-field references and coefficients
4. Python simulation runs without errors
5. Physical behavior matches expectations (e.g., energy transfer)

## Known Issues & Solutions

| Issue                     | Cause                                | Solution                                        |
| ------------------------- | ------------------------------------ | ----------------------------------------------- |
| Component equations = 0   | Field strength not expanded          | Replace F with derivative expr before decompose |
| Cross-field not detected  | Fields not in coordinate form        | Pass `additionalFields={...}`                   |
| "Symbol X already used"   | xAct kernel caching                  | Use `xTensorQ[X]` check before `DefManifold`    |
| Mixed 2/3-arg Derivatives | `ConvertCDToDerivatives` forms       | Parser handles both                             |
| xAct Cycles context       | `System`Cycles`vs`xAct`xPerm`Cycles` | Use `Head[x] === xAct`xPerm`Cycles`             |

## Quick Reference

### File Locations

- Wolfram modules: `tidal/wolfram/`
- Python pipeline: `tidal/symbolic/`
- CLI: `tidal/cli/`
- Examples: `examples/{name}/`, JSON: `examples/data/*.json`

### Key Functions (Wolfram)

- `DecomposeToComponents[eom, field, chart, additionalFields, "MetricMatrix" -> matrix]`
- `BuildMultiFieldJSONStructure[fieldEquations, metadata]`
- `ExtractMassCouplingFromEquations[fieldEquations]`
- `IsNonConstantMetric[metricMatrix, coords]`
- `EvaluateEpsilonComponents[epsilon, chart]`

### Key Functions (Python)

- `build_pde_from_json(json_path, parameters={...})`
- `EquationSystem.from_dict(data)` — auto-computes mass/coupling matrices
- `PDEFromSpec(spec, parameters=...)` — dynamic PDE construction
- `LHSStructure.from_data(data)`, `parse_field_name(name)`

---

See detailed notes in: `/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory/`
