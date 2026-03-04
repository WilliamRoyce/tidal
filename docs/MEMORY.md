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
- `gauge_fixing.md` - Gauge fixing tutorial and preset reference
- `background_fields.md` - Background field architecture documentation
- `adaptive_timestepping.md` - Solver architecture and tolerance control
- `constraint_fields.md` - Constraint field physics and solver treatment
- `solver_migration.md` - Migration from py-pde to native SUNDIALS/numpy architecture

---

## Overview

This project implements a symbolic physics pipeline: Lagrangian (xAct/Mathematica) → JSON → native PDE solver (SUNDIALS IDA/CVODE, leapfrog, scipy) with numpy spatial operators. Users should NEVER manually hardcode equations — all PDEs must derive from the Lagrangian via symbolic computation.

The project operates exclusively in the **linearised regime**: all Lagrangians are quadratic, producing linear PDEs. Higher-order time derivatives (> 2) are out of scope (Ostrogradsky instabilities).

## Critical Architecture Decisions

### Solver Architecture (SUNDIALS + numpy)

The PDE time-stepping layer uses four backends, automatically selected based on equation structure:

| Backend      | Library                  | Use Case                            | Key Feature                                    |
| ------------ | ------------------------ | ----------------------------------- | ---------------------------------------------- |
| **IDA**      | SUNDIALS (scikit-sundae) | DAE systems (algebraic constraints) | Implicit Newton iteration, sparse Jacobian     |
| **CVODE**    | SUNDIALS (scikit-sundae) | Adaptive ODE (wave equations)       | BDF, tolerance control (`--rtol`/`--atol`)     |
| **Leapfrog** | Native (numpy)           | Symplectic integration              | Exact energy conservation (shadow Hamiltonian) |
| **scipy**    | scipy.integrate          | General-purpose adaptive ODE        | DOP853, Radau, BDF via `solve_ivp`             |

**Spatial operators** are pure numpy (`tidal/solver/operators.py`): 2nd-order FD stencils for laplacian, gradient, cross_derivative, with periodic / Dirichlet / Neumann BC support.

**Key solver classes** (all in `tidal/solver/`):

- `FieldSet` (`fields.py`): typed container, contiguous flat array, zero-copy views
- `CoefficientEvaluator` (`coefficients.py`): 4-level cache (L0 preresolved → L1 expression → L2 spatial grid → L3 per-call)
- `RHSEvaluator` (`rhs.py`): unified operator + coefficient application
- `SpecValidator` (`validation.py`): operator dimensions, field references, CFL, mass sign
- `StateLayout` (`state.py`): maps field names to flat-array slices
- `GridInfo` (`grid.py`): minimal frozen dataclass replacing py-pde's `CartesianGrid`

**Constraint pre-solve** (`constraint_solve.py`): three-tier solver (FFT → sparse matrix probe → automatic selection) runs before IDA to ensure consistent initial conditions. Gauge regularisation for singular Poisson with periodic BCs.

**References:** Hindmarsh et al. (2005) for SUNDIALS; Hairer et al. (2006) for symplectic integrators. See `docs/solver_migration.md` for the full migration architecture, `docs/adaptive_timestepping.md` for tolerance control.

### Multi-Field Coupling

- **Cross-field terms:** gradient operators in JSON (e.g., `gradient_x(A_2)` in A_0 equation)
- **Field transformation:** `DecomposeToComponents` accepts `additionalFields` parameter
- **JSON detection:** `IdentifyMultiFieldTerm` uses case-insensitive matching + function head extraction

### Dimension Handling (Generalised)

- **Dynamic dimension:** `ComponentDecompose.wl` uses `dim = Length[ScalarsOfChart[chart]]`
- **Derivative conversion:** `CommonUtilities.wl` handles both 2-arg (1+1D) and 3-arg (2+1D) Derivative forms
- **Supports:** 1+1D (t,x), 2+1D (t,x,y), 3+1D (t,x,y,z), extensible to higher

### E-L Velocity Form (Current)

State vector stores velocities `v = dq/dt` directly (not canonical momenta). The JSON `equations[]` array contains original Euler-Lagrange equations as derived by the Wolfram pipeline. No Legendre transform, no `field_rates`, no K-inversion at runtime.

- **State:** `(q, v, A₀)` where `v_X = dX/dt`
- **Kinematic:** `dq/dt = v` (trivial identity)
- **Dynamic:** `dv/dt = RHS` (E-L equations from JSON)
- **JSON canonical structure:** `hamiltonian_terms` (bilinear H for energy) + optional `volume_element` (sqrt|g_spatial| for curved coordinates)
- **Velocity naming:** `v_{field_name}` (e.g., `v_phi_0`, `v_A_1`), slot kind `"velocity"`
- **SimulationData:** `.velocities` attribute stores `v = dq/dt` (renamed from `.momenta`)

**Volume-weighted energy:** For curved coordinates, `E = ∫ H(x) √|g_spatial| d^n x`. The Wolfram pipeline computes `sqrt|det(g_spatial)|` from `MetricMatrix[[2;;, 2;;]]` and injects it as `canonical.volume_element` in the JSON. Omitted for flat spacetimes (fast path: scalar 1.0 multiply).

### Mass/Coupling Matrix Auto-Computation (Phase 12)

- **Convention:** `matrix[i][j] = -(coefficient of identity(field_j) in equation_i)`
- **Wolfram:** `ExtractMassCouplingFromEquations` in ExportJSON.wl auto-computes at export
- **Python:** `_compute_matrices_from_terms` in json_loader.py auto-computes at load (defence-in-depth)
- **Symbolic preservation:** `mass_matrix_symbolic` / `coupling_matrix_symbolic` preserve exact Mathematica expressions

### Christoffel Symbol Computation (Auto-Detection)

- **Auto-detection (default):** `DecomposeToComponents[eom, field, chart, {}, "MetricMatrix" -> matrix]`
- **Constant metric** (flat or static conformal) → All Christoffels = 0
- **Non-constant metric** (time/position-dependent) → Compute from standard formula
- **Override available:** `"ComputeChristoffels" -> True/False` for explicit control

### Background Fields (Position-Dependent Coupling)

Non-dynamical tensors declared via `[[background_fields]]` TOML. Survive as position-dependent coefficients in the EOM. 4-level caching in `CoefficientEvaluator`. Supports scalar, vector, and tensor backgrounds.

See `docs/background_fields.md` for the full pipeline trace.

### Gauge Fixing (Optional)

Per-field `[[gauge]]` TOML section. Named presets: Lorenz, de Donder, Coulomb, temporal, axial. Custom expressions supported. Expression-based extensibility — adding a new gauge requires one function in `GaugeFix.wl` + one registry entry in `_derive.py`.

See `docs/gauge_fixing.md` for tutorial and developer guide.

## xAct/Wolfram Patterns

### Symbol Management

- **Always check existence:** `If[!xTensorQ[M2], DefManifold[M2, ...]]`
- **Use shared symbols:** M2, M3, eta, CD for examples to avoid conflicts

### Epsilon Tensors (Automated)

- `DefMetric` automatically creates `epsiloneta3[-a,-b,-c]`
- `EvaluateEpsilonComponents` in CommonUtilities.wl → numeric ±1
- Minkowski (-,+,+): covariant ε_012 = −1, contravariant ε^012 = +1

### xAct Cycles Context (Critical)

- xAct's `Cycles` is `xAct`xPerm`Cycles`, NOT `System`Cycles`
- All pattern matching must use `_xAct`xPerm`Cycles` and `Head[x] === xAct`xPerm`Cycles`

### Lagrangian Construction Rules

1. Use **explicit metric tensors** (`eta[a,b]`) for index raising
2. Use **`DefConstantSymbol[m2]`** for mass/coupling constants (not bare `Symbol`)
3. **Parenthesize multi-line expressions** in .wls: `L = (term1 + term2);`
4. **Expand field strength before decompose:** `L = -1/2 CD[-a][A[-b]] ...` (not via F substitution)

## CLI (`tidal` Command)

The `tidal` CLI provides 9 subcommands with zero new dependencies (stdlib argparse + tomllib):

| Command                     | Description                                                           |
| --------------------------- | --------------------------------------------------------------------- |
| `tidal derive theory.toml`  | Generate .wls from TOML, run wolframscript to produce JSON            |
| `tidal simulate spec.json`  | Full simulation with plotting (`--param`, `--ic`, `--bc`, `--scheme`) |
| `tidal measure result_dir/` | Extract physics measurements (energy, conversion, mixing, spectra)    |
| `tidal inspect spec.json`   | Display equation system info (fields, operators, parameters)          |
| `tidal list`                | Discover available JSON specs in `examples/data/`                     |
| `tidal validate spec.json`  | Validate JSON equation specification structure                        |
| `tidal plot result_dir/`    | Standalone plotting from simulation output directories                |
| `tidal sweep spec.json`     | Automated parameter sweep with parallel execution and measurement     |
| `tidal analyze sweep_dir/`  | Sobol/Morris sensitivity analysis on completed sweep results          |

**Simulation resume:** `tidal simulate spec.json --resume output_dir/ [--snapshot N] [--t-additional T] --t-end T`
Continues from any saved snapshot. Inherits grid/params/BC from metadata. Output goes to new directory.

**TOML Configuration:**

- `theory.toml` with spacetime, fields, constants, Lagrangian, parameters
- `[[derived_fields]]` for intermediate tensors (e.g., field strength F_ab)
- `[[background_fields]]` for non-dynamical tensors (e.g., external B-field)
- `[[gauge]]` for optional per-field gauge fixing (Lorenz, de Donder, etc.)
- IC presets: `gaussian`, `plane-wave`, `zero`, `formula`, `file`, `noise`
- Per-axis boundary conditions via `--bc neumann,periodic`
- Solver selection via `--scheme cvode|ida|scipy|leapfrog`
- Tolerance control via `--rtol`, `--atol`

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

`identity`, `laplacian`, `laplacian_{x,y,z}`, `gradient_{x,y,z}`, `cross_derivative_{xy,xz,yz}`, `first_derivative_t`, `biharmonic`, `mixed_T_S1_S2_...`

- All support cross-field application and velocity (`v_i`) references. `first_derivative_t` resolves to velocity slot at runtime. `mixed_*` operators decompose into time derivatives (velocity/EOM RHS) + spatial gradients.

## Example Implementations

| Example                   | Dim  | Key Features                                                       |
| ------------------------- | ---- | ------------------------------------------------------------------ |
| `scalar_field/`           | 1+1D | KG, mass term, dispersion                                          |
| `electromagnetic/`        | 1+1D | Maxwell, Lorenz gauge                                              |
| `proca/`                  | 1+1D | Massive vector (Proca mass)                                        |
| `coupled_scalars/`        | 1+1D | Cross-field coupling, mass matrix                                  |
| `scalar_potential_well/`  | 1+1D | Background potential well, `[[background_fields]]`                 |
| `chern_simons/`           | 2+1D | Epsilon tensor, A_0 constraint                                     |
| `elasticity/`             | 2+1D | Anisotropic laplacian, cross_derivative_xy                         |
| `curved_spacetime/`       | 2+1D | Hubble friction, time-dependent coefficients                       |
| `sphere_kg/`              | 2+1D | Position-dependent coefficients, S²                                |
| `polar_kg/`               | 2+1D | Polar coordinates, Christoffel auto-detection                      |
| `electrostatics/`         | 2+1D | Poisson constraint, no time evolution                              |
| `scalar_vector_coupling/` | 2+1D | Mixed-rank cross-field (scalar+vector), 4 constants                |
| `massive_gravity/`        | 2+1D | Linearized massive gravity, xPert, coupled constraints             |
| `coupled_proca/`          | 2+1D | Two massive vectors, coupled Helmholtz, periodic BCs               |
| `coupled_scattering/`     | 2+1D | Position-dependent Gaussian coupling, background fields            |
| `proca_background/`       | 2+1D | Lorentzian scalar BG, two Proca vectors, constraint+BG integration |
| `vector_background/`      | 2+1D | Tanh domain wall vector BG, ComponentValue mechanism               |
| `scalar_field_3d/`        | 3+1D | Full 4D KG                                                         |
| `spherical_kg/`           | 3+1D | Spherical coordinates, trig coefficients                           |
| `cylindrical_kg/`         | 3+1D | Cylindrical, mixed curved/flat                                     |
| `gravitational_waves/`    | 3+1D | xPert linearization, TT gauge, constraints                         |
| `massive_3form/`          | 3+1D | Rank-3 antisymmetric, symmetry reduction 64→4                      |

## Testing Guidelines

### Test Counts

- **~1242 Python tests** + **~115 Wolfram tests** passing
- Run: `uv run pytest tests/` and `./scripts/full_test.sh`

### Verification Pattern

1. Wolfram derivation produces N component equations
2. JSON has correct dimension, signature, coordinates
3. JSON terms show proper cross-field references and coefficients
4. Python simulation runs without errors
5. Physical behaviour matches expectations (e.g., energy transfer)

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
- Solver: `tidal/solver/` (ida.py, leapfrog.py, cvode.py, fields.py, operators.py, grid.py)
- CLI: `tidal/cli/`
- Measurement: `tidal/measurement/`
- Examples: `examples/{name}/`, JSON: `examples/data/*.json`

### Key Functions (Wolfram)

- `DecomposeToComponents[eom, field, chart, additionalFields, "MetricMatrix" -> matrix]`
- `BuildMultiFieldJSONStructure[fieldEquations, metadata]`
- `ExtractMassCouplingFromEquations[fieldEquations]`
- `IsNonConstantMetric[metricMatrix, coords]`
- `EvaluateEpsilonComponents[epsilon, chart]`

### Key Functions (Python)

- `EquationSystem.from_dict(data)` — loads JSON, auto-computes mass/coupling matrices
- `FieldSet.from_spec(spec, grid)` — creates typed field container
- `CoefficientEvaluator(spec, grid)` — evaluates position-dependent coefficients
- `RHSEvaluator(spec, coeff_eval, grid)` — applies operators with coefficients
- `StateLayout.from_spec(spec)` — maps field names to flat-array slices
- `GridInfo(bounds, shape, periodic)` — minimal grid descriptor

---

See detailed developer notes in: `/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/memory/`
