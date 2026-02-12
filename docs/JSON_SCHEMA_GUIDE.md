# JSON Schema Reference

Complete reference for the JSON specification format used as the interface between the Wolfram/xAct symbolic layer and the Python/py-pde simulation layer.

**Source of truth:** [`tidal/symbolic/json_loader.py`](../tidal/symbolic/json_loader.py)

---

## Top-Level Structure

```json
{
  "metadata":  { ... },
  "spacetime": { ... },
  "fields":    [ ... ],
  "equations": [ ... ],
  "coupling":  { ... }
}
```

| Key | Required | Description |
|-----|----------|-------------|
| `metadata` | No | Source info, default parameters, gauge choice |
| `spacetime` | **Yes** | Dimension, signature, coordinates |
| `fields` | **Yes** | List of field components |
| `equations` | **Yes** | Equations of motion (one per field) |
| `coupling` | **Yes** | Mass/coupling matrices (auto-computed from terms) |

---

## 1. `metadata`

Informational block documenting the Lagrangian origin. Not parsed by the PDE solver except for `parameters`.

| Key | Type | Description |
|-----|------|-------------|
| `source` | string | `"xAct"` (derived symbolically) or `"analytical"` (hand-constructed) |
| `lagrangian_expr` | string | Human-readable Lagrangian expression |
| `derived_from` | string | Always `"Euler-Lagrange"` for xAct-derived equations |
| `gauge` | string | Gauge choice: `"none"`, `"lorenz"`, `"coulomb"`, `"de_donder"` |
| `linearized` | boolean | Whether equations are linearized perturbations |
| `construction` | string | `"analytical"` for hand-built JSON |
| `parameters` | object | Default parameter values, e.g. `{"H": 0.1, "m2": 1.0}` |

**Example:**
```json
{
  "metadata": {
    "source": "xAct",
    "lagrangian_expr": "-1/2 g^ab nabla_a phi nabla_b phi - 1/2 m^2 phi^2",
    "derived_from": "Euler-Lagrange",
    "gauge": "none",
    "parameters": {"dSH": 0.1, "dSm2": 1.0}
  }
}
```

---

## 2. `spacetime`

Defines the geometric background.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `dimension` | integer | **Yes** | Total spacetime dimension: 2 (1+1D), 3 (2+1D), 4 (3+1D) |
| `signature` | array[int] | No | Metric signature, e.g. `[-1, 1, 1]` for Minkowski |
| `coordinates` | array[string] | No | Coordinate names in order. First is always time. Defaults to `["t","x","y","z"]` truncated to `dimension`. |
| `metric_type` | string | No | `"minkowski"`, `"conformal_flat"`, `"product_conformal"`, `"polar"`, etc. |
| `conformal_factor` | string | No | Symbolic conformal factor (Mathematica form) |
| `metric_components` | string | No | Metric description, e.g. `"diag(-1, 1, x^2)"` |

**Examples:**
```json
{"spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]}}
```
```json
{"spacetime": {"dimension": 3, "coordinates": ["t", "x", "y"],
               "metric_type": "conformal_flat", "conformal_factor": "exp(H*t)"}}
```
```json
{"spacetime": {"dimension": 3, "coordinates": ["t", "x", "y"],
               "metric_type": "polar", "metric_components": "diag(-1, 1, x^2)"}}
```

---

## 3. `fields`

Array of field components. Order determines the state layout.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `name` | string | **Yes** | Unique field name (see naming conventions) |
| `index` | integer | **Yes** | Numeric index (0, 1, 2, ...) |
| `is_dynamical` | boolean | **Yes** | Whether the field evolves in time (always `true`) |

### Field Naming Conventions

| Format | Pattern | Examples |
|--------|---------|----------|
| Standard | `Base_Index` | `A_0`, `phi_1`, `u_2` |
| Tensor | `Base_Component_Index` | `stress_xy_0`, `u_x_1` |
| Compact | `BaseIndex` | `phi0`, `A1`, `psi2` |
| Simple | `Base` | `phi`, `psi`, `rho` (index defaults to 0) |

**Validation:** All names must be unique. Index values must be sequential starting from 0.

**Example:**
```json
{
  "fields": [
    {"name": "A_0", "index": 0, "is_dynamical": true},
    {"name": "A_1", "index": 1, "is_dynamical": true},
    {"name": "A_2", "index": 2, "is_dynamical": true}
  ]
}
```

---

## 4. `equations`

Array of equations of motion, one per field. Length must equal the number of fields.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `field` | string | **Yes** | Field name this equation governs (must match a `fields[].name`) |
| `lhs` | object | **Yes** | Left-hand side structure |
| `rhs` | object | **Yes** | Right-hand side (operator terms) |
| `constraint_solver` | object | No | Elliptic constraint config (only valid when `lhs.order.time == 0`) |

### 4a. LHS Structure

Specifies what time derivative appears on the left-hand side.

| Key | Type | Required | Default | Description |
|-----|------|----------|---------|-------------|
| `expression` | string | **Yes** | — | Human-readable LHS, e.g. `"d2_t(phi_0)"`. Informational only. |
| `order.time` | integer | **Yes** | 2 | Time derivative order |
| `order.space` | integer | No | 0 | Space derivative order on LHS (typically 0) |

**PDE type is determined by `order.time`:**

| `order.time` | PDE Type | Example |
|--------------|----------|---------|
| 0 | Elliptic (constraint) | Poisson: `nabla^2 phi = rho` |
| 1 | Parabolic | Heat: `d_t phi = nabla^2 phi` |
| 2 | Hyperbolic (wave) | Wave: `d2_t phi = nabla^2 phi - m^2 phi` |
| 3+ | Higher-order | Supported but uncommon |

**Examples:**
```json
{"lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}}}
{"lhs": {"expression": "A_0",         "order": {"time": 0, "space": 0}}}
```

### 4b. RHS Structure

The right-hand side is always a linear combination of operator terms.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `type` | string | **Yes** | Must be `"linear_combination"` |
| `terms` | array | **Yes** | List of operator terms (can be empty for trivial equations) |

### 4c. Operator Terms

Each term represents `coefficient * operator(field)`.

| Key | Type | Required | Default | Description |
|-----|------|----------|---------|-------------|
| `coefficient` | float | **Yes** | — | Numeric coefficient |
| `operator` | string | **Yes** | — | Differential operator name (see [Operator Reference](#5-operator-reference)) |
| `field` | string | **Yes** | — | Target field or momentum reference |
| `coefficient_symbolic` | string | No | `null` | Symbolic expression (Mathematica InputForm) |
| `time_dependent` | boolean | No | `false` | Whether coefficient depends on time |
| `coordinate_dependent` | array[string] | No | `[]` | Coordinate names the coefficient depends on |

**Example (constant coefficient):**
```json
{"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"}
```

**Example (parametric — resolved from user-supplied parameters):**
```json
{"coefficient": -1.0, "operator": "identity", "field": "phi_0",
 "coefficient_symbolic": "-m2"}
```

**Example (time-dependent — de Sitter mass term):**
```json
{"coefficient": 1.0, "operator": "identity", "field": "dSphi_0",
 "coefficient_symbolic": "-(dSm2*E^(2*dSH*t[]))",
 "time_dependent": true, "coordinate_dependent": ["t"]}
```

**Example (position-dependent — 1/r in polar coordinates):**
```json
{"coefficient": 1.0, "operator": "gradient_x", "field": "polphi_0",
 "coefficient_symbolic": "x[]^(-1)",
 "coordinate_dependent": ["x"]}
```

**Example (cross-field reference):**
```json
{"coefficient": -0.5, "operator": "identity", "field": "chi_0"}
```
This term appears in the `phi_0` equation but acts on `chi_0`.

**Example (momentum reference — mixed time-space derivative):**
```json
{"coefficient": -1.0, "operator": "gradient_x", "field": "pi_1"}
```
`pi_1` references the time derivative (momentum) of field at index 1. This encodes `d_x(d_t(A_1))`.

---

## 5. Operator Reference

### Static Operators

| Operator | Min Grid Dim | Description | Math |
|----------|-------------|-------------|------|
| `identity` | 1 | Field itself (mass/coupling) | `f` |
| `laplacian` | 1 | Full spatial Laplacian | `nabla^2 f` |
| `laplacian_x` | 1 | Second derivative in x | `d^2f/dx^2` |
| `laplacian_y` | 2 | Second derivative in y | `d^2f/dy^2` |
| `laplacian_z` | 3 | Second derivative in z | `d^2f/dz^2` |
| `gradient_x` | 1 | First derivative in x | `df/dx` |
| `gradient_y` | 2 | First derivative in y | `df/dy` |
| `gradient_z` | 3 | First derivative in z | `df/dz` |
| `cross_derivative_xy` | 2 | Mixed partial in x and y | `d^2f/dxdy` |
| `cross_derivative_xz` | 3 | Mixed partial in x and z | `d^2f/dxdz` |
| `cross_derivative_yz` | 3 | Mixed partial in y and z | `d^2f/dydz` |
| `first_derivative_t` | 1 | First time derivative | `df/dt` |
| `biharmonic` | 1 | Fourth-order Laplacian | `nabla^4 f` |

### Generic Operators (Pattern-Matched)

| Pattern | Example | Description |
|---------|---------|-------------|
| `derivative_N_x` | `derivative_3_x` | Nth derivative along x: `d^3f/dx^3` |
| `derivative_N_y` | `derivative_5_y` | Nth derivative along y: `d^5f/dy^5` |
| `derivative_N_z` | `derivative_4_z` | Nth derivative along z: `d^4f/dz^4` |
| `derivative_Nx_My` | `derivative_2x_1y` | Mixed: `d^3f/dx^2 dy` |

### Special: `first_derivative_t`

Used for Hubble friction in curved spacetime. Represents the first time derivative of a field. Unlike spatial operators, it accesses the momentum variable directly from the state vector.

### Field References

A term's `field` value can be:
- **A field name** from the `fields` list (e.g., `"phi_0"`, `"A_1"`, `"chi_0"`)
- **A momentum reference** `"pi_N"` where N is the field index (0-based). This accesses `d_t(field_N)`. Used for mixed time-space derivatives, e.g., `gradient_x(pi_1)` = `d_x(d_t(A_1))`.

**Validation:** Momentum index N must satisfy `0 <= N < n_components`.

---

## 6. Coefficient Types

### Constant

No symbolic expression. The numeric `coefficient` value is used directly.

```json
{"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"}
```

### Parametric

A `coefficient_symbolic` that is a simple parameter name (or its negation). Resolved at runtime from the `parameters` dict passed to `build_pde_from_json()`.

```json
{"coefficient": -1.0, "operator": "identity", "field": "phi_0",
 "coefficient_symbolic": "-polm2"}
```

At runtime: `build_pde_from_json(path, parameters={"polm2": 0.5})` resolves this to `+0.5`.

Convention: a leading `-` in the symbolic name means the parameter value is negated. So `"-polm2"` with `polm2=0.5` gives `coefficient = +0.5`.

### Time-Dependent

Coefficient varies with time. Requires both `time_dependent: true` and `coordinate_dependent: ["t"]`.

```json
{"coefficient": 1.0, "operator": "identity", "field": "dSphi_0",
 "coefficient_symbolic": "-(dSm2*E^(2*dSH*t[]))",
 "time_dependent": true, "coordinate_dependent": ["t"]}
```

The expression is evaluated at each timestep using Python's `eval()` with math functions and parameter values injected.

### Position-Dependent

Coefficient varies with spatial coordinates. Listed in `coordinate_dependent`.

```json
{"coefficient": 1.0, "operator": "laplacian_y", "field": "polphi_0",
 "coefficient_symbolic": "x[]^(-2)",
 "coordinate_dependent": ["x"]}
```

Produces an ndarray matching the grid shape. The `x[]` is Mathematica's xCoba notation for the coordinate function — converted to plain `x` at runtime.

### Mixed (Time + Position)

Both time and spatial dependence. Set `time_dependent: true` and list all coordinates in `coordinate_dependent`.

```json
{"coefficient_symbolic": "exp(2*H*t[])*x[]^(-2)",
 "time_dependent": true, "coordinate_dependent": ["t", "x"]}
```

### Mathematica-to-Python Conversion

Symbolic expressions are stored in Mathematica InputForm. Automatic conversion handles:

| Mathematica | Python | Notes |
|-------------|--------|-------|
| `E^(...)` | `exp(...)` | Euler's number |
| `Power[x,y]` | `(x)**(y)` | Function form |
| `^` | `**` | Infix form |
| `Sin[x]`, `Cos[x]`, `Tan[x]` | `sin(x)`, `cos(x)`, `tan(x)` | |
| `Cot[x]`, `Sec[x]`, `Csc[x]` | `cot(x)`, `sec(x)`, `csc(x)` | scipy.special |
| `ArcSin[x]`, `ArcCos[x]` | `arcsin(x)`, `arccos(x)` | |
| `ArcTan[x, y]` | `arctan2(y, x)` | Argument order swapped! |
| `Sinh[x]`, `Cosh[x]`, `Tanh[x]` | `sinh(x)`, `cosh(x)`, `tanh(x)` | |
| `Log[x]`, `Sqrt[x]`, `Abs[x]` | `log(x)`, `sqrt(x)`, `abs(x)` | |
| `Erf[x]` | `erf(x)` | scipy.special |
| `BesselJ[n,x]`, `BesselY[n,x]` | `jv(n,x)`, `yv(n,x)` | scipy.special |
| `t[]`, `x[]`, `y[]` | `t`, `x`, `y` | xCoba coordinate symbols |
| `[`, `]` | `(`, `)` | Bracket conversion |

---

## 7. `coupling`

Mass and coupling matrices extracted from `identity` operator terms. Only symbolic matrices are stored in JSON; numeric values are auto-computed at load time from symbolic expressions and metadata parameters.

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `mass_matrix_symbolic` | array[array[string\|null]] | **Yes** | n x n symbolic mass matrix entries |
| `coupling_matrix_symbolic` | array[array[string\|null]] | **Yes** | n x n symbolic coupling matrix entries |

**Note:** Prior to the Numeric Matrix Cleanup, JSON also contained numeric `mass_matrix` and `coupling_matrix` arrays. These have been removed — numeric values are now computed by `EquationSystem.from_dict()` using symbolic expressions resolved against `metadata.parameters`.

### Convention

```
matrix[i][j] = -(coefficient of identity(field_j) in equation_i)
```

This makes mass-squared positive for the standard Lagrangian sign convention: `d2_t phi = ... - m^2 phi` produces `mass_matrix[i][i] = m^2`.

- Diagonal entries (i == j) go in `mass_matrix_symbolic`
- Off-diagonal entries (i != j) go in `coupling_matrix_symbolic`
- Null entries indicate zero coupling

### Auto-Computation

**Important:** `EquationSystem.from_dict()` always recomputes numeric matrices from the equation terms via `_compute_matrices_from_terms`. The `_resolve_symbolic_coeff` helper resolves symbolic expressions (e.g., `"-m2"`) against the `metadata.parameters` dict to produce correct numeric values.

**Example (coupled scalars):**

Given `phi_0` equation has term `{"coefficient": -1.0, "operator": "identity", "field": "phi_0", "coefficient_symbolic": "-m2"}` and `{"coefficient": -0.5, "operator": "identity", "field": "chi_0", "coefficient_symbolic": "-g"}`:

```json
{
  "coupling": {
    "mass_matrix_symbolic": [["-m2", null], [null, "-mchi2"]],
    "coupling_matrix_symbolic": [[null, "-g"], ["-g", null]]
  }
}
```

---

## 8. `constraint_solver`

Per-equation configuration for elliptic constraint solving. Supports Poisson, Helmholtz, algebraic, anisotropic, and coupled multi-field constraints.

**Valid only when `lhs.order.time == 0`.** Raises `ValueError` otherwise.

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `enabled` | boolean | `false` | Whether to solve at each timestep |
| `method` | string | `"auto"` | Solver method: `"auto"` (selects best), `"poisson"`, `"fft"`, `"sparse"` |
| `max_iterations` | integer | `100` | Max iterations for iterative solvers (Gauss-Seidel) |
| `tolerance` | float | `1e-10` | Convergence tolerance for iterative solvers |
| `boundary_conditions` | object | `{}` | Per-axis BCs keyed by coordinate name |

### Method Selection

| Method | When Used | Description |
|--------|-----------|-------------|
| `"auto"` | Default | FFT for periodic grids, sparse matrix for non-periodic |
| `"fft"` | Periodic BCs | Fast spectral solve; supports coupled block solve via SVD |
| `"sparse"` | Non-periodic BCs | Sparse matrix assembly + direct solve |
| `"poisson"` | Legacy alias | Equivalent to `"auto"` |

**Coupled constraints:** When multiple constraint equations reference each other's fields (e.g., coupled Proca A_0/B_0), the solver automatically detects coupling and solves them simultaneously via block matrix (FFT) or Gauss-Seidel iteration (sparse).

### Boundary Condition Values

| Key | Type | Description |
|-----|------|-------------|
| `type` | string | `"periodic"`, `"dirichlet"`, or `"neumann"` |
| `value` | float | Fixed value (Dirichlet only) |
| `derivative` | float | Fixed normal derivative (Neumann only) |

**Example (Helmholtz with Dirichlet BCs):**
```json
{
  "constraint_solver": {
    "enabled": true,
    "method": "auto",
    "max_iterations": 30,
    "tolerance": 1e-10,
    "boundary_conditions": {
      "x": {"type": "dirichlet", "value": 0.0},
      "y": {"type": "dirichlet", "value": 0.0}
    }
  }
}
```

**Example (Poisson with periodic BCs):**
```json
{
  "constraint_solver": {
    "enabled": true,
    "method": "auto",
    "boundary_conditions": {
      "x": {"type": "periodic"},
      "y": {"type": "periodic"}
    }
  }
}
```

---

## 9. Validation Rules

These checks run automatically when loading JSON via `load_equation_system()`.

### Schema Validation (`validate_json_schema`)

1. Top-level keys `spacetime`, `fields`, `equations` must be present
2. `spacetime.dimension` must be an integer
3. `fields` must be non-empty; each entry must have a `name` key
4. `equations` must be non-empty; each entry must have `field` and `rhs` keys

### Construction Validation (`EquationSystem.from_dict`)

5. All field names must be unique
6. Each equation's `field` must match a name in the `fields` list
7. `rhs.type` must be `"linear_combination"`

### Operator Validation (`OperatorTerm.from_dict`)

8. Each `operator` must be a static operator OR match a generic pattern (`derivative_N_x`, `derivative_Nx_My`)
9. Each term must have `coefficient`, `operator`, and `field` keys

### Field Reference Validation (`_validate_field_references`)

10. Regular field references must match a name in the `fields` list
11. Momentum references (`pi_N`) must have a valid numeric index: `0 <= N < n_components`

### Constraint Validation (`ComponentEquation.__post_init__`)

12. `constraint_solver.enabled: true` is only valid when `lhs.order.time == 0`

### Matrix Validation (`EquationSystem.__post_init__`)

13. `mass_matrix` and `coupling_matrix` must be n x n (matching `n_components`)
14. A warning (not error) is emitted if matrices don't match auto-computed values

---

## 10. Complete Examples

### Klein-Gordon (1+1D, Minimal)

Source: [`examples/data/klein_gordon_1d.json`](../examples/data/klein_gordon_1d.json)

```json
{
  "metadata": {
    "source": "xAct",
    "lagrangian_expr": "-1/2 CD[-a][phi[]] CD[a][phi[]] - 1/2 m2 phi[]^2",
    "derived_from": "Euler-Lagrange",
    "gauge": "none"
  },
  "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
  "fields": [{"name": "phi_0", "index": 0, "is_dynamical": true}],
  "equations": [{
    "field": "phi_0",
    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
    "rhs": {
      "type": "linear_combination",
      "terms": [
        {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
        {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"}
      ]
    }
  }],
  "coupling": {"mass_matrix_symbolic": [[null]], "coupling_matrix_symbolic": [[null]]}
}
```

### Coupled Scalars (Cross-Field Coupling)

Source: [`examples/data/coupled_scalars.json`](../examples/data/coupled_scalars.json)

Two scalar fields phi and chi with mutual coupling. Note how `phi_0`'s equation includes `identity` acting on `chi_0`:

```json
{
  "fields": [
    {"name": "phi_0", "index": 0, "is_dynamical": true},
    {"name": "chi_0", "index": 1, "is_dynamical": true}
  ],
  "equations": [
    {
      "field": "phi_0",
      "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
      "rhs": {"type": "linear_combination", "terms": [
        {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
        {"coefficient": -0.5, "operator": "identity", "field": "chi_0"},
        {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"}
      ]}
    },
    {
      "field": "chi_0",
      "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
      "rhs": {"type": "linear_combination", "terms": [
        {"coefficient": -4.0, "operator": "identity", "field": "chi_0"},
        {"coefficient": -0.5, "operator": "identity", "field": "phi_0"},
        {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"}
      ]}
    }
  ],
  "coupling": {
    "mass_matrix_symbolic": [["-m2", null], [null, "-mchi2"]],
    "coupling_matrix_symbolic": [[null, "-g"], ["-g", null]]
  }
}
```

### Polar KG (Position-Dependent Coefficients)

Source: [`examples/data/polar_kg.json`](../examples/data/polar_kg.json)

Christoffel corrections produce 1/r gradient and 1/r^2 angular Laplacian:

```json
{
  "spacetime": {"dimension": 3, "coordinates": ["t", "x", "y"],
                "metric_type": "polar", "metric_components": "diag(-1, 1, x^2)"},
  "equations": [{
    "field": "polphi_0",
    "lhs": {"expression": "d2_t(polphi_0)", "order": {"time": 2}},
    "rhs": {"type": "linear_combination", "terms": [
      {"coefficient": -1.0, "operator": "identity", "field": "polphi_0",
       "coefficient_symbolic": "-polm2"},
      {"coefficient": 1.0, "operator": "gradient_x", "field": "polphi_0",
       "coefficient_symbolic": "x[]^(-1)", "coordinate_dependent": ["x"]},
      {"coefficient": 1.0, "operator": "laplacian_x", "field": "polphi_0"},
      {"coefficient": 1.0, "operator": "laplacian_y", "field": "polphi_0",
       "coefficient_symbolic": "x[]^(-2)", "coordinate_dependent": ["x"]}
    ]}
  }]
}
```

### De Sitter KG (Time-Dependent Coefficients)

Source: [`examples/data/de_sitter_kg.json`](../examples/data/de_sitter_kg.json)

Hubble friction and exponentially growing mass from conformal factor:

```json
{
  "spacetime": {"dimension": 3, "coordinates": ["t", "x", "y"],
                "metric_type": "conformal_flat", "conformal_factor": "exp(H*t)"},
  "equations": [{
    "field": "dSphi_0",
    "lhs": {"expression": "d2_t(dSphi_0)", "order": {"time": 2}},
    "rhs": {"type": "linear_combination", "terms": [
      {"coefficient": 1.0, "operator": "identity", "field": "dSphi_0",
       "coefficient_symbolic": "-(dSm2*E^(2*dSH*t[]))",
       "time_dependent": true, "coordinate_dependent": ["t"]},
      {"coefficient": -1.0, "operator": "first_derivative_t", "field": "dSphi_0",
       "coefficient_symbolic": "-dSH"},
      {"coefficient": 1.0, "operator": "laplacian_x", "field": "dSphi_0"},
      {"coefficient": 1.0, "operator": "laplacian_y", "field": "dSphi_0"}
    ]}
  }]
}
```

### Chern-Simons (Constraint + Momentum References)

Source: [`examples/data/chern_simons_3d.json`](../examples/data/chern_simons_3d.json)

A_0 is a constraint (time_order=0), A_1 and A_2 are second-order. Uses `pi_0`, `pi_1`, `pi_2` momentum references and `first_derivative_t`:

```json
{
  "fields": [
    {"name": "A_0", "index": 0, "is_dynamical": true},
    {"name": "A_1", "index": 1, "is_dynamical": true},
    {"name": "A_2", "index": 2, "is_dynamical": true}
  ],
  "equations": [
    {
      "field": "A_0",
      "lhs": {"expression": "A_0", "order": {"time": 0, "space": 0}},
      "rhs": {"type": "linear_combination", "terms": [
        {"coefficient": 0.5, "operator": "gradient_x", "field": "A_2"},
        {"coefficient": -0.5, "operator": "gradient_y", "field": "A_1"},
        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_0"},
        {"coefficient": 1.0, "operator": "laplacian_y", "field": "A_0"},
        {"coefficient": -1.0, "operator": "gradient_x", "field": "pi_1"},
        {"coefficient": -1.0, "operator": "gradient_y", "field": "pi_2"}
      ]}
    },
    {
      "field": "A_1",
      "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2, "space": 0}},
      "rhs": {"type": "linear_combination", "terms": [
        {"coefficient": 0.5, "operator": "first_derivative_t", "field": "A_2"},
        {"coefficient": -1.0, "operator": "gradient_x", "field": "pi_0"},
        {"coefficient": -0.5, "operator": "gradient_y", "field": "A_0"},
        {"coefficient": 1.0, "operator": "laplacian_y", "field": "A_1"},
        {"coefficient": -1.0, "operator": "cross_derivative_xy", "field": "A_2"}
      ]}
    }
  ]
}
```

### Electrostatics (Constraint Solver)

Source: [`examples/data/electrostatics_2d.json`](../examples/data/electrostatics_2d.json)

Poisson equation solved as a constraint at each timestep:

```json
{
  "equations": [{
    "field": "phi",
    "lhs": {"expression": "phi", "order": {"time": 0, "space": 0}},
    "rhs": {"type": "linear_combination", "terms": [
      {"coefficient": 1.0, "operator": "laplacian", "field": "phi"},
      {"coefficient": 1.0, "operator": "identity", "field": "rho"}
    ]},
    "constraint_solver": {
      "enabled": true,
      "method": "poisson",
      "boundary_conditions": {
        "x": {"type": "dirichlet", "value": 0.0},
        "y": {"type": "dirichlet", "value": 0.0}
      }
    }
  }]
}
```

---

## 11. Python API

### Loading and Simulating

```python
from tidal.symbolic import load_equation_system, build_pde_from_json

# Load spec (validates JSON)
spec = load_equation_system("examples/data/polar_kg.json")

# Build PDE with runtime parameters
pde = build_pde_from_json("examples/data/polar_kg.json", parameters={"polm2": 0.5})
```

### Key `EquationSystem` Properties

| Property | Type | Description |
|----------|------|-------------|
| `n_components` | `int` | Number of field components |
| `dimension` | `int` | Spacetime dimension |
| `spatial_dimension` | `int` | `dimension - 1` |
| `component_names` | `tuple[str, ...]` | Field names in order |
| `time_orders` | `tuple[int, ...]` | Per-component time derivative order |
| `state_size` | `int` | Total state variables (2nd-order: 2, others: 1) |
| `state_layout` | `tuple[tuple[str,str], ...]` | `(name, "field"/"momentum")` pairs |
| `effective_coordinates` | `tuple[str, ...]` | Coordinate names (inferred if not set) |
| `spatial_coordinates` | `tuple[str, ...]` | Non-time coordinates |
| `mass_matrix` | `tuple[tuple[float,...],...]` | Auto-computed mass matrix |
| `coupling_matrix` | `tuple[tuple[float,...],...]` | Auto-computed coupling matrix |
| `mass_matrix_symbolic` | `tuple[tuple[str\|None,...],...]` | Symbolic mass entries |
