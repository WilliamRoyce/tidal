# Solver Migration: py-pde → SUNDIALS/IDA

## Overview

Replace TIDAL's py-pde numerical backend with SUNDIALS/IDA (via scikit-sundae)
for native DAE support, eliminating the need for bolted-on constraint solvers,
symbolic K^{-1} inversion, and manual Hamiltonian splitting.

**Status**: Implementation complete (February 2026). IDA, CVODE, leapfrog, and scipy backends all operational. State management (`FieldSet`, `StateLayout`), spatial operators, and grid infrastructure (`GridInfo`) fully implemented.

---

## 1. Motivation

### 1.1 The Problem

TIDAL derives equations of motion from Lagrangians via xAct/Mathematica. These
equations are naturally a **mixed system**:

- **2nd-order wave equations** (Klein-Gordon, Proca, linearized gravity):
  `∂²_t φ = ∇²φ - m²φ`
- **1st-order evolution** (Chern-Simons, curved spacetime friction):
  `∂_t A_i = RHS`
- **Algebraic constraints** (Gauss's law, Hamiltonian constraint):
  `0 = ∇²A_0 - ∂_x π_1 - ∂_y π_2`

This is a **differential-algebraic equation (DAE)**, classified as index-1 when
the constraint Jacobian is nonsingular.

py-pde solves only explicit ODEs of the form `du/dt = F(u)`. To force our DAE
into this framework, we accumulated workarounds:

| Workaround                                            | Location                   | Problem                                                       |
| ----------------------------------------------------- | -------------------------- | ------------------------------------------------------------- |
| 3-pass `evolution_rate()`                             | `pde_builder.py:2919-2993` | Fragile ordering; constraints must be solved before evolution |
| Custom constraint solvers (FFT, matrix, Gauss-Seidel) | `pde_builder.py:1846-2750` | ~900 lines of custom elliptic solver code                     |
| K^{-1} symbolic inversion                             | `_derive.py:2110-2166`     | Complicated expressions; fails for large kinetic blocks       |
| Virtual momenta                                       | `pde_builder.py:2759-2793` | Ad-hoc mechanism for 1st-order equations                      |
| Manual (q, π) splitting                               | `pde_builder.py:842-983`   | State layout management complexity                            |

### 1.2 The Solution

**SUNDIALS/IDA** (Hindmarsh et al., 2005) solves the general implicit DAE:

```
F(t, y, y') = 0
```

where `y` is a flat state vector containing **all** variables (fields, momenta,
constraints). The `algebraic_idx` parameter marks which components are algebraic
(constraints). IDA uses variable-order BDF with Newton iteration, automatically
handling:

- Mixed differential + algebraic variables
- Non-diagonal kinetic matrices (mass matrices)
- Stiff coupling
- Consistent initial conditions

**scikit-sundae** (NREL, BSD-3) provides a clean Python interface to SUNDIALS/IDA.

### 1.3 Research Basis

Three parallel research surveys informed this design:

1. **PDE Framework Landscape** — Assessed Dedalus, FEniCS, PETSc, SUNDIALS,
   diffrax, Devito, scipy_dae, diffeqpy against TIDAL's requirements.
   Conclusion: SUNDIALS/IDA is the best fit for DAE time integration with
   clean Python API. (Dedalus is the best full-framework alternative but
   requires spectral spatial discretization — future Phase E.)

2. **Constrained Hamiltonian Evolution** — Surveyed Einstein Toolkit, SpECTRE,
   NRPy+, Kranc, Z4c/CCZ4, MEEP/Yee, lattice gauge theory. Consensus:
   **constraint damping beats constraint solving** (free evolution + damping
   source terms, not elliptic solves at every timestep).

3. **TIDAL Architecture Analysis** — Mapped exact py-pde API surface. Found
   clean interface boundary: py-pde provides 5 classes (~15 methods); TIDAL
   provides all physics. 80% of spatial operators already implemented as
   custom numpy. Position-dependent coefficient system is pure numpy/eval.

### 1.4 Key References

- Hindmarsh, Brown, Grant, Lee, Serban, Shumaker, Woodward (2005), "SUNDIALS:
  Suite of Nonlinear and Differential/Algebraic Equation Solvers", ACM TOMS
  31(3):363-396.
- Dedner, Kemm, Kröner, Munz, Schnitzer (2002), "Hyperbolic Divergence
  Cleaning for MHD", J. Comput. Phys. 175(2):645-673.
- Gundlach, Calabrese, Hinder, Martin-Garcia (2005), "Constraint damping in
  the Z4 formulation and harmonic gauge", Class. Quantum Grav. 22:3767.
- Bernuzzi, Hilditch (2010), "Constraint violation in free evolution schemes",
  Phys. Rev. D 81:084003.
- Hairer, Lubich, Wanner (2006), "Geometric Numerical Integration", Springer.
- Zwicker (2020), "py-pde: A Python package for solving partial differential
  equations", JOSS 5(48):2158.
- Brown (2009), "Covariant formulations of BSSN and the standard gauge",
  Phys. Rev. D 79:104029.
- Alic, Bona-Casas, Bona, Rezzolla, Palenzuela (2013), "Conformal and
  covariant formulation of the Z4 system", Phys. Rev. D 85:064040.
- Yoshida (1990), "Construction of higher order symplectic integrators",
  Phys. Lett. A 150:262-268.

---

## 2. Architecture

### 2.1 Component Replacement Map

| Component         | Old (py-pde)                        | New                                 | Rationale                                               |
| ----------------- | ----------------------------------- | ----------------------------------- | ------------------------------------------------------- |
| Time integration  | `PDEBase.solve()` (RK4/scipy)       | **scikit-sundae IDA**               | Native DAE; handles constraints + mass matrices         |
| Field containers  | `FieldCollection` / `ScalarField`   | **Plain numpy arrays** + `GridInfo` | Remove unnecessary wrapper layer                        |
| Grid              | `CartesianGrid`                     | **`GridInfo` dataclass**            | Tailored to TIDAL (bounds, shape, periodic, dx, coords) |
| Spatial operators | `field.laplace(bc)` / custom numpy  | **Custom numpy module**             | Already 80% custom; consolidate                         |
| Storage           | `MemoryStorage` / `CallbackTracker` | **Existing `SnapshotWriter`**       | Already numpy-based disk-backed storage                 |
| Kinetic matrix    | K^{-1} symbolic inversion           | **K passed directly** to IDA        | Simpler expressions; IDA handles implicit solve         |
| Symplectic        | N/A                                 | **Störmer-Verlet** (custom)         | Energy-conserving option for Hamiltonian systems        |

### 2.2 Solver Selection

Users select **one scheme** per simulation via `--scheme`:

| Scheme          | Solver            | Use Case            | Constraints        | Energy            |
| --------------- | ----------------- | ------------------- | ------------------ | ----------------- |
| `ida` (default) | SUNDIALS/IDA      | All systems         | Algebraic (native) | Dissipative OK    |
| `leapfrog`      | Störmer-Verlet    | Hamiltonian systems | Not supported      | Conserved exactly |
| `scipy`         | scipy `solve_ivp` | Simple ODEs         | Not supported      | Dissipative OK    |

### 2.3 Data Flow

```
theory.toml
    │
    ▼ (tidal derive)
wolframscript (xAct)
    │ Lagrangian → EOM → K, S, RHS terms, H
    ▼
spec.json (kinetic_matrix, spatial_momenta, equations, hamiltonian)
    │
    ▼ (tidal simulate --scheme ida)
json_loader.py → EquationSystem (kinetic_matrix, spatial_momenta, equations)
    │
    ▼
ida.py → build_residual_fn(spec, grid)
    │         F_i = K_{ij} · y'_j - (π_i - S_i)     [field slots]
    │         F_i = y'_i - RHS_i                      [momentum slots]
    │         F_i = RHS_i                              [constraint slots]
    ▼
sksundae.ida.IDA(residual_fn, algebraic_idx, ...)
    │
    ▼
result.y → SnapshotWriter → .npy files → tidal plot / tidal measure
```

---

## 3. Key Design Decision: K Instead of K^{-1}

### 3.1 Current Pipeline (K^{-1})

The Wolfram pipeline currently (`_derive.py:2054-2220`):

1. Computes canonical momenta: `π_i = ∂L/∂(∂_t q_i)`
2. Extracts kinetic matrix: `K_{ij} = ∂π_i/∂(∂_t q_j)`
3. **Inverts K symbolically**: `K^{-1} = Inverse[K]`
4. Computes field*rates: `∂_t q_i = K^{-1}*{ij} · (π_j - S_j)`
5. Exports K^{-1} coefficients embedded in `field_rates` JSON

This produces complicated expressions for non-diagonal K:

- Chern-Simons (2×2): `dA_1/dt = π_1 - (κ/2)·A_2 + ∇_x A_0`
- Linearized gravity (7×7): coefficients like `-4/3`, `-2/3` from rational K^{-1}

### 3.2 New Pipeline (K directly)

With IDA, Hamilton's 1st equation in residual form:

```
K_{ij} · ∂_t q_j - π_i + S_i = 0
```

No K^{-1} needed. IDA's Newton iteration solves the implicit system.

The Wolfram pipeline simplifies:

1. Compute π_i — same as before
2. Extract K — same as before
3. ~~Invert K~~ — removed
4. ~~Compute field_rates~~ — removed
5. Export K entries + spatial momenta S directly

### 3.3 JSON Format Change

**Old** (`field_rates` with K^{-1} embedded):

```json
"field_rates": {
  "A_1": [
    {"coefficient": 1.0, "operator": "identity", "field": "pi_1"},
    {"coefficient": -0.5, "operator": "identity", "field": "A_2",
     "coefficient_symbolic": "-kappa/2"},
    {"coefficient": 1.0, "operator": "gradient_x", "field": "A_0"}
  ]
}
```

**New** (`kinetic_matrix` + `spatial_momenta`):

```json
"kinetic_matrix": {
  "entries": [
    {"i": 0, "j": 0, "value": 1.0},
    {"i": 0, "j": 1, "value": 0.5, "symbolic": "kappa/2"},
    {"i": 1, "j": 0, "value": -0.5, "symbolic": "-kappa/2"},
    {"i": 1, "j": 1, "value": 1.0}
  ],
  "dimension": 2
},
"spatial_momenta": {
  "A_1": [{"coefficient": 1.0, "operator": "gradient_x", "field": "A_0"}],
  "A_2": [{"coefficient": 1.0, "operator": "gradient_y", "field": "A_0"}]
}
```

### 3.4 Impact on Störmer-Verlet

Leapfrog needs explicit velocities `∂_t q = f(π)`. For non-diagonal K, solve
the small N_dyn × N_dyn system `K · v = π - S` per timestep using
`np.linalg.solve()`. This is 2×2 to 10×10 — microseconds per call.

---

## 4. Position-Dependent Coefficients

The coefficient evaluation system is **pure numpy/eval** with no py-pde
dependency. It transfers to the new architecture unchanged.

### 4.1 Pipeline

1. **JSON**: `coefficient_symbolic` string + `coordinate_dependent` metadata
2. **Conversion**: `mathematica_to_python()` in `_eval_utils.py` (regex + string
   transforms: `E^x → exp(x)`, `Power[x,y] → x**y`, etc.)
3. **Evaluation**: `eval(py_expr, namespace)` where namespace includes numpy
   functions and grid coordinate arrays
4. **Caching**: 4-level hierarchy:
   - L0: preresolved constants (at `__init__`)
   - L1: Mathematica→Python string cache (`_expr_cache`)
   - L2: spatial grid arrays (`_spatial_cache`) — position-dependent,
     time-independent terms evaluated once
   - L3: per-call evaluation for time+position-dependent terms

### 4.2 Grid Coordinate Injection

Currently: `grid.cell_coords` from py-pde's `CartesianGrid`
New: `grid.cell_coords` from TIDAL's `GridInfo` — same numpy array, different
source.

### 4.3 Verified Examples

| Example            | Coefficient Type           | Expression                        |
| ------------------ | -------------------------- | --------------------------------- |
| sphere_kg          | x-dependent Laplacian      | `x[]^2/(2*sphR^2)`                |
| coupled_scattering | Gaussian spatial coupling  | `-(g0/E^((x[]^2+y[]^2)/(2*R^2)))` |
| curved_spacetime   | Hubble friction (time-dep) | `-2*H(t)`                         |
| vector_background  | tanh domain wall           | `Tanh[(x[]-x0)/w]`                |

All evaluated by `_resolve_coefficient_at_point()` which returns `float` (scalar)
or `np.ndarray` (grid-shaped) depending on coordinate dependence.

### 4.4 Constraint Solver Interaction

- **FFT solver**: Rejects position-dependent self-terms (raises `ValueError`).
  Only handles constant-coefficient Poisson-type constraints.
- **Matrix solver**: Handles position-dependent source terms. Operator matrices
  assembled with constant coefficients.
- **IDA**: Position-dependent coefficients appear in the residual function
  naturally. No special handling needed — `compute_rhs()` evaluates them per
  call, using the spatial cache for time-independent terms.

---

## 5. Implementation Checklist

### Step 0: Documentation

- [x] Create this documentation file (`docs/solver_migration.md`)

### Phase 1: Foundation (GridInfo + Operators)

#### 1a. GridInfo dataclass

- [ ] Create `tidal/solver/__init__.py`
- [ ] Create `tidal/solver/grid.py` with `GridInfo` dataclass
  - Properties: `bounds`, `shape`, `periodic`, `dx`, `cell_coords`, `ndim`, `num_points`
  - Test: `tests/test_solver_grid.py`
- [ ] Verify `GridInfo.cell_coords` matches `CartesianGrid.cell_coords` output
      for all grid configurations (1D, 2D, 3D, periodic/non-periodic)

#### 1b. Spatial operators

- [ ] Create `tidal/solver/operators.py`
  - Functions: `laplacian`, `gradient`, `directional_laplacian`,
    `cross_derivative`, `identity`, `biharmonic`
  - Source: Refactor from `pde_builder.py:71-304`
  - BC handling: ghost-cell padding for Neumann/Dirichlet
  - Test: `tests/test_solver_operators.py`
- [ ] Verify numerical equivalence with py-pde operators on test grids
- [ ] Operator registry mapping string names → functions

#### 1c. Position-dependent coefficients

- [ ] Verify `_resolve_coefficient_at_point()` works with `GridInfo`
  - Change type signature: `GridBase → GridInfo`
  - `cell_coords` property must have compatible shape
- [ ] Run sphere_kg, coupled_scattering, curved_spacetime, vector_background
      coefficient evaluation tests
- [ ] Verify `_spatial_cache` produces identical arrays

### Phase 2: Wolfram Pipeline Changes

#### 2a. Modify `_derive.py`

- [ ] In `_wls_canonical_pipeline()` (lines 2054-2220):
  - Keep K computation (lines 2089-2096)
  - Keep K validation / det(K) check (lines 2101-2107)
  - Keep spatial momentum S computation (lines 2120-2126)
  - Remove K^{-1} inversion (lines 2110-2112)
  - Remove K^{-1} field_rates expansion (lines 2128-2166)
  - Add: emit K entries to WLS metadata for JSON export
  - Add: emit S terms to WLS metadata for JSON export

#### 2b. Modify `ExportJSON.wl`

- [ ] Add `"kinetic_matrix"` section to canonical JSON block
  - Sparse entries: `{"i": i, "j": j, "value": numeric, "symbolic": string}`
  - Dimension field
- [ ] Add `"spatial_momenta"` section
  - Per-field list of OperatorTerms (same format as existing RHS terms)
- [ ] Remove or make optional the `"field_rates"` section

#### 2c. Modify `json_loader.py`

- [ ] Parse `kinetic_matrix` entries into numpy array on `CanonicalStructure`
- [ ] Parse `spatial_momenta` into dict of OperatorTerm tuples
- [ ] Update `CanonicalStructure` dataclass
- [ ] Backward compatibility: if `field_rates` present (old JSON), still parse

#### 2d. Verification

- [ ] Re-derive chern_simons with new pipeline — verify K matches:
      K = [[1, κ/2], [-κ/2, 1]]
- [ ] Re-derive coupled_proca — verify K is diagonal (identity)
- [ ] Re-derive linearized_gravity — verify K is 7×7 with known structure
- [ ] All existing Wolfram tests pass

### Phase 3: SUNDIALS/IDA Integration

#### 3a. State management

- [ ] Create `tidal/solver/state.py`
  - `state_to_flat()`: dict of field arrays → flat numpy array
  - `flat_to_fields()`: flat numpy array → dict of field arrays
  - `StateLayout` class: maps slot indices to field names and types
  - Test: `tests/test_solver_state.py`

#### 3b. IDA residual builder

- [ ] Create `tidal/solver/ida.py`
  - `build_residual_fn(spec, grid)`: returns IDA-compatible residual function
  - Residual logic:
    - Constraint slots: `F_i = RHS_i` (algebraic)
    - Momentum slots: `F_i = y'_i - RHS_i` (Hamilton's 2nd)
    - Field slots: `F_i = K_{ij} · y'_j - (π_i - S_i)` (Hamilton's 1st)
    - 1st-order slots: `F_i = y'_i - RHS_i`
  - `build_algebraic_indices(spec, grid)`: returns list of algebraic indices
  - `build_sparsity_pattern(spec, grid)`: returns sparse matrix for Jacobian
  - Test: `tests/test_solver_ida.py`

#### 3c. Solver setup and integration

- [ ] IDA solver creation with appropriate options:
  - `linsolver='sparse'` with sparsity pattern
  - `calc_initcond='yp0'` for consistent initial conditions
  - Tolerances: `rtol=1e-8, atol=1e-10`
  - Thread control: `nthreads=-1`
- [ ] Snapshot callback: extract fields from flat state, write via SnapshotWriter
- [ ] Error handling: IDA convergence failures → clear error messages

#### 3d. CLI integration

- [ ] Modify `_simulate.py`:
  - `--scheme ida` (default): dispatch to SUNDIALS/IDA
  - `--scheme leapfrog`: dispatch to Störmer-Verlet
  - `--scheme scipy`: dispatch to scipy solve_ivp (no constraints)
  - Remove all py-pde solver dispatch code
  - Build `GridInfo` instead of `CartesianGrid`
- [ ] Modify `pyproject.toml`:
  - Remove `py-pde` from dependencies
  - Add `scikit-sundae` as required dependency
- [ ] Update `--help` text for new scheme options

#### 3e. Verification

- [ ] EM plane-wave (`--scheme ida`): amplitude stable (< 1.2× initial)
- [ ] EM plane-wave: A_0 constraint satisfied automatically (no separate solver)
- [ ] Coupled Proca (`--scheme ida`): coupled Helmholtz constraints solved
- [ ] Coupled scalars: energy conservation matches py-pde baseline
- [ ] Chern-Simons with κ: correct dispersion relation
- [ ] sphere_kg: position-dependent coefficients handled correctly
- [ ] curved_spacetime: time-dependent Hubble friction works

### Phase 4: Störmer-Verlet Integrator

- [ ] Create `tidal/solver/leapfrog.py`
  - `stormer_verlet()` function (~20 lines)
  - `solve_kinetic(K, rhs)`: identity shortcut + `np.linalg.solve`
  - Snapshot callback support
  - Test: `tests/test_solver_leapfrog.py`
- [ ] CLI `--scheme leapfrog` dispatch
- [ ] Verify: KG wave — energy drift < 1e-10 over 1000 oscillation periods
- [ ] Verify: coupled_scalars — energy conserved with leapfrog vs decaying with ida

### Phase 5: Test Migration

- [ ] Replace `from pde import ...` in all 17 test files
  - `CartesianGrid → GridInfo`
  - `ScalarField → np.ndarray`
  - `FieldCollection → dict[str, np.ndarray]`
  - `MemoryStorage → SnapshotWriter` or direct npy
- [ ] Update `conftest.py` fixtures
- [ ] Verify: `uv run pytest tests/ -q --tb=no` passes (all ~1067 tests)

### Phase 6: Cleanup

- [ ] Remove py-pde from pyproject.toml
- [ ] `uv run ruff check tidal/ tests/` — 0 violations
- [ ] `uv run pyright tidal/` — 0 errors (strict mode)
- [ ] Update `docs/NEXT_PHASES.md` to reflect completed migration
- [ ] Update MEMORY.md with new architecture notes

### Future: Phase 7 (FieldSet — Typed Field Container)

- [ ] Create `tidal/solver/fieldset.py` with `FieldSet` class
- [ ] Own field data: `dict[str, np.ndarray]` with typed access and grid shape enforcement
- [ ] Consolidate field name validation: `parse_field_name()`, momentum naming conventions
      (`pi_phi_0`), index validation — single source of truth
- [ ] Handle serialization: flat vector ↔ named fields (from `state.py`), NPZ I/O
      (from `_io.py`), snapshot packing
- [ ] Bridge py-pde boundary: `to_field_collection()` / `from_field_collection()` only
      at CLI output edge
- [ ] Replace raw `dict[str, np.ndarray]` throughout solver and measurement modules
- [ ] Enforce invariants: grid shape consistency, field/momentum pairing, immutable metadata

### Future: Phase 8 (Constraint Damping / Dedner GLM)

- [ ] TOML: `[gauge] method = "dedner"`
- [ ] Auto-compute c_h, c_p from grid/domain
- [ ] Wolfram pipeline: add auxiliary ψ field, modify EOM
- [ ] Constraint monitoring: L2/L∞ norms at each snapshot
- [ ] `constraints.npy` output alongside field data

---

## 6. IDA Technical Details

### 6.1 Residual Function Signature

```python
def resfn(t: float, y: np.ndarray, yp: np.ndarray, res: np.ndarray) -> None:
    """
    t:   Current time
    y:   State vector (N,) — all fields, momenta, constraints flattened
    yp:  Time derivatives (N,) — y' values
    res: Output residual (N,) — write F(t,y,y') in-place
    """
```

### 6.2 Solver Setup

```python
from sksundae.ida import IDA

solver = IDA(
    resfn,
    algebraic_idx=algebraic_indices,  # constraint grid points
    calc_initcond='yp0',              # auto-consistent initial yp
    linsolver='sparse',               # SuperLU_MT (N ≤ 10^5)
    sparsity=sparsity_matrix,         # Jacobian sparsity pattern
    nthreads=-1,                      # all cores for SuperLU_MT
    rtol=1e-8,
    atol=1e-10,
)

result = solver.solve(t_span, y0, yp0)
# result.t: (n_times,) time points
# result.y: (n_times, N) solution
# result.success: bool
# result.nfev: residual evaluation count
```

### 6.3 Sparsity Pattern

The Jacobian `J_{ij} = ∂F_i/∂y_j + α · ∂F_i/∂y'_j` has sparsity determined by
the equation coupling structure:

- Spatial operators (laplacian, gradient) couple nearest neighbors on the grid
- Cross-field terms couple different field slots at the same grid point
- The kinetic matrix K couples field slots at the same grid point

The sparsity pattern can be pre-computed from the equation specification without
evaluating any numerical values.

### 6.4 Performance Expectations

| Grid Size | N (state) | Linear Solver          | Expected Step Time |
| --------- | --------- | ---------------------- | ------------------ |
| 64        | ~640      | dense                  | < 1 ms             |
| 64×64     | ~40,000   | sparse (SuperLU_MT)    | ~10 ms             |
| 128×128   | ~160,000  | sparse or GMRES        | ~100 ms            |
| 256×256   | ~650,000  | GMRES + preconditioner | ~1 s               |

### 6.5 Index-1 DAE Requirement

IDA requires index-1 DAEs. TIDAL's constraint equations are index-1 when:

- Gauss's law `∇²A_0 = source`: invertible if BCs are periodic or Dirichlet
  (Neumann on all boundaries has a kernel — requires gauge fixing)
- Coupled Helmholtz: `(∇² - m²)A_0 = source`: always invertible for m² > 0

If a constraint is structurally singular (index > 1), IDA will fail with a
convergence error. This is the correct behavior — fail fast, not silently wrong.

---

## 7. Störmer-Verlet Technical Details

### 7.1 Algorithm

For separable Hamiltonian `H = T(π) + V(q)`:

```
Hamilton's equations:
  dq_i/dt = ∂T/∂π_i = K^{-1}_{ij} (π_j - S_j)    [velocity from momentum]
  dπ_i/dt = -∂V/∂q_i = RHS_i(q, ∇q)               [force from configuration]

Störmer-Verlet (one step):
  π_{n+1/2} = π_n + (dt/2) · force(q_n)             [half-kick]
  v_{n+1/2} = solve(K, π_{n+1/2} - S(q_n))          [velocity from momentum]
  q_{n+1}   = q_n + dt · v_{n+1/2}                   [drift]
  π_{n+1}   = π_{n+1/2} + (dt/2) · force(q_{n+1})   [half-kick]
```

### 7.2 Properties

- **Symplectic**: preserves the symplectic 2-form exactly
- **Time-reversible**: integrator is its own adjoint
- **Shadow Hamiltonian**: preserves H̃ = H + O(dt²) exactly → no energy drift
- **2nd-order accurate**: global error O(dt²)
- **No energy drift**: unlike RK4 which drifts O(dt⁴) per step, accumulating

### 7.3 Limitations

- **Not for dissipative systems**: absorbing BCs, friction terms, Dedner damping
  all break the Hamiltonian structure
- **No algebraic constraints**: cannot handle Gauss's law directly
- **Fixed timestep**: no adaptive stepping (CFL condition must be satisfied)
- **Separable H required**: H = T(π) + V(q), not H(q, π) with mixed terms

### 7.4 Higher-Order Methods

4th-order Yoshida (composition of leapfrog steps with different dt):

```python
# Yoshida (1990) coefficients
c1 = c4 = 1 / (2 * (2 - 2**(1/3)))
c2 = c3 = (1 - 2**(1/3)) / (2 * (2 - 2**(1/3)))
d1 = d3 = 1 / (2 - 2**(1/3))
d2 = -2**(1/3) / (2 - 2**(1/3))
# Apply 3 leapfrog substeps with these coefficients
```

---

## 8. Constraint Handling Strategy

### 8.1 Short Term: IDA Algebraic Equations

IDA solves constraints as algebraic equations within the DAE. This is equivalent
to constraint projection at every timestep — the Newton iteration ensures the
algebraic equations are satisfied to within tolerance.

### 8.2 Medium Term: Constraint Damping

For systems where hyperbolic propagation of constraint errors is preferred
(e.g., absorbing boundaries), add damping source terms:

```
Modified EOM: ∂_t π_i = RHS_i - κ · C_i
```

where `C_i = 0` is the constraint and κ is the damping parameter.

### 8.3 Long Term: Dedner GLM

For vector fields, Dedner divergence cleaning adds an auxiliary scalar ψ:

```
∂_t A_i += ∂_i ψ
∂_t ψ = -c_h² (∇·A) - (c_p²/c_h²) ψ
```

Parameters c_h (propagation speed) and c_p (damping rate) are auto-computed
from the grid spacing and domain size.

---

## 9. Files Reference

### New Files

| File                             | Purpose                       |
| -------------------------------- | ----------------------------- |
| `tidal/solver/__init__.py`       | Solver package                |
| `tidal/solver/grid.py`           | `GridInfo` dataclass          |
| `tidal/solver/operators.py`      | Spatial operators (numpy)     |
| `tidal/solver/state.py`          | State flattening/unflattening |
| `tidal/solver/ida.py`            | SUNDIALS/IDA integration      |
| `tidal/solver/leapfrog.py`       | Störmer-Verlet integrator     |
| `tests/test_solver_grid.py`      | GridInfo tests                |
| `tests/test_solver_operators.py` | Operator tests                |
| `tests/test_solver_state.py`     | State management tests        |
| `tests/test_solver_ida.py`       | IDA integration tests         |
| `tests/test_solver_leapfrog.py`  | Leapfrog tests                |

### Modified Files

| File                            | Changes                                   |
| ------------------------------- | ----------------------------------------- |
| `tidal/symbolic/pde_builder.py` | Decouple from py-pde types                |
| `tidal/symbolic/json_loader.py` | Parse `kinetic_matrix`, `spatial_momenta` |
| `tidal/cli/_simulate.py`        | New solver dispatch                       |
| `tidal/cli/_derive.py`          | Export K (not K^{-1})                     |
| `tidal/wolfram/ExportJSON.wl`   | Emit new canonical sections               |
| `tidal/measurement/_writer.py`  | Remove py-pde, constraint norms           |
| `tidal/measurement/_io.py`      | Remove py-pde, use GridInfo               |
| `pyproject.toml`                | py-pde → scikit-sundae                    |
| `tests/conftest.py`             | Replace fixtures                          |
| 17 test files                   | Replace py-pde imports                    |

### Unchanged Files

| File                                  | Why                     |
| ------------------------------------- | ----------------------- |
| `tidal/symbolic/_eval_utils.py`       | Already pure numpy/eval |
| `tidal/wolfram/CommonUtilities.wl`    | No solver dependency    |
| `tidal/wolfram/ComponentDecompose.wl` | No solver dependency    |
| `tidal/wolfram/EulerLagrange.wl`      | No solver dependency    |
