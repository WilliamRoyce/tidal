# Adaptive Time-Stepping

**Status:** Implemented (CVODE + scipy + IDA + leapfrog)
**Date:** 2026-02-23

Tolerance-controlled adaptive time-stepping replaces manual `dt` selection
with error-controlled solvers that automatically choose step sizes. Users
specify accuracy targets (via `--rtol`/`--atol`) instead of guessing a
stable `dt`.

## Problem Statement

Fixed-step integrators (leapfrog) require the user to manually satisfy the
CFL condition: `dt < dx / c_max`. This has three problems:

1. **Fragility.** A `dt` that works for one parameter set can be unstable
   for another. Getting it wrong produces garbage output that can look
   plausible.

2. **Inefficiency.** A `dt` chosen for the fastest mode wastes compute on
   slow phases. Resonant mixing (the Gertsenshtein effect) has beat
   frequencies much slower than the carrier — fixed steps resolve every
   carrier cycle even when nothing interesting is happening. Berlin et al.
   (2024, arXiv:2405.08865) document these numerical challenges.

3. **No accuracy guarantee.** Leapfrog energy error is O(dt²) — the
   shadow Hamiltonian drifts by an amount proportional to dt². There is no
   built-in error control to bound this.

## Four Solver Paths

### `--scheme auto` (default): Intelligent Selection

Auto-selection **always** picks a tolerance-controlled adaptive solver.
Leapfrog is opt-in only via `--scheme leapfrog`.

Detection algorithm:

| Priority | Condition | Solver | Reason |
|----------|-----------|--------|--------|
| 1 | Constraint equations (time_order=0) | IDA | Algebraic constraints need DAE residual form |
| 2 | Modal eligible (flat metric + periodic BCs + time-independent + supported operators) | **Modal** | Machine-precision eigendecomposition, O(1) in t_end. See [modal_solver.md](modal_solver.md) |
| 3 | First-order equations (time_order=1) | IDA | Diffusion/transport need implicit handling |
| 4 | Dissipation (`first_derivative_t` in RHS) | IDA | Friction breaks symplecticity, BDF handles well |
| 5 | No canonical Hamiltonian | CVODE (warning) | Missing canonical section — likely hand-crafted JSON |
| 6 | Pure wave (2nd-order, Hamiltonian) | CVODE | Adaptive BDF, no DAE overhead |

### `--scheme cvode`: SUNDIALS CVODE (adaptive ODE)

Pure ODE solver using SUNDIALS CVODE via scikit-sundae. Supports BDF
(stiff, order 1-5) and Adams (non-stiff, order 1-12).

```bash
tidal simulate spec.json --scheme cvode                 # BDF default
tidal simulate spec.json --scheme cvode --method Adams   # Non-stiff
tidal simulate spec.json --scheme cvode --rtol 1e-6      # Custom tolerance
```

**Best for:** Wave systems (the auto-selected default for pure wave
equations). BDF handles the mild stiffness from spatial discretization
(eigenvalues scale with 1/dx²). In-place RHS evaluation avoids
per-call allocation.

Reference: Hindmarsh et al. (2005), "SUNDIALS: Suite of Nonlinear and
Differential/Algebraic Equation Solvers", ACM TOMS 31(3).

### `--scheme modal`: Fourier Modal Solver (exact eigendecomposition)

Transforms spatial grid to Fourier space, builds per-mode evolution matrix,
and eigendecomposes for exact y(t) = exp(A·t)·y₀ solutions. No time-stepping,
no CFL condition, no tolerance parameter needed.

```bash
tidal simulate spec.json --scheme modal                  # Auto for eligible
```

**Best for:** Periodic-domain wave systems with time-independent coefficients
(the auto-selected default for eligible systems). Machine-precision accuracy
(~10⁻¹⁶ error). Cost is O(1) in simulation time — eigendecomposition is done
once, evaluation at any t is just `exp(λ·Δt)`. Speedup over CVODE: 57-1451×
on coupled_scalars (N=64-256), 3.4-7.7× on Gertsenshtein (N=256, t=50-500).

**Requires:** Flat metric, all-periodic BCs, time-independent coefficients,
supported Fourier operators. Position-dependent coefficients are handled via
k-space convolution. See [modal_solver.md](modal_solver.md) for full details.

### `--scheme ida`: SUNDIALS IDA (adaptive DAE)

DAE solver using SUNDIALS IDA. Handles algebraic constraints via residual
form F(t, y, y') = 0.

```bash
tidal simulate spec.json --scheme ida                    # Auto for DAE
tidal simulate spec.json --scheme ida --rtol 1e-6        # Custom tolerance
```

**Best for:** Constraint systems (Chern-Simons, coupled Proca), dissipative
systems (Hubble friction), and any system with mixed time-derivative orders.

### `--scheme scipy`: scipy solve_ivp (adaptive ODE)

Wraps `scipy.integrate.solve_ivp`. Default method is DOP853 (8th-order
Dormand-Prince), an explicit one-step method with excellent accuracy per
function evaluation.

```bash
tidal simulate spec.json --scheme scipy                  # DOP853 default
tidal simulate spec.json --scheme scipy --method RK45    # Lower-order explicit
tidal simulate spec.json --scheme scipy --method Radau   # Implicit (stiff)
tidal simulate spec.json --scheme scipy --method BDF     # Implicit multi-step
```

**Best for:** Smooth non-stiff wave problems where explicit high-order
methods outperform implicit BDF (fewer function evaluations per accuracy
unit). Also useful when scipy's Radau is preferred over SUNDIALS.

Reference: Dormand & Prince (1980), "A family of embedded Runge-Kutta
formulae", J. Comp. Appl. Math. 6, 19-26.

### `--scheme leapfrog`: Symplectic (opt-in, fixed dt)

Stormer-Verlet symplectic integrator. Preserves a shadow Hamiltonian to
machine precision — zero secular energy drift.

```bash
tidal simulate spec.json --scheme leapfrog               # CFL-auto dt
tidal simulate spec.json --scheme leapfrog --dt 0.01     # Manual dt
```

**Best for:** Long-time Hamiltonian evolution where zero energy drift
matters more than per-step accuracy. Not suitable for dissipative systems,
constraints, or when tolerance control is needed.

Reference: Hairer, Lubich & Wanner (2006), "Geometric Numerical
Integration", Springer, Chapter VI.

### Choosing a Solver

| Scenario | Recommended |
|----------|-------------|
| Default (any system) | `--scheme auto` (always adaptive) |
| Wave equation, production accuracy | `--scheme cvode --rtol 1e-8` |
| Constraint system (A_0, gauge fields) | `--scheme ida` (auto-selected) |
| Dissipative system (friction, damping) | `--scheme ida` (auto-selected) |
| Smooth non-stiff, high accuracy | `--scheme scipy --rtol 1e-10` |
| Long-time energy conservation | `--scheme leapfrog` |
| Stiff system (m²dx² >> 1) | `--scheme scipy --method Radau` |
| Convergence study | `--scheme cvode --rtol 1e-12 --atol 1e-14` |

## CLI Reference

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--scheme` | `{auto,ida,leapfrog,cvode,scipy}` | `auto` | Solver selection |
| `--rtol` | float | 1e-8 | Relative tolerance (adaptive solvers) |
| `--atol` | float | 1e-10 | Absolute tolerance (adaptive solvers) |
| `--method` | str | auto | Integration method (cvode: BDF/Adams; scipy: DOP853/RK45/Radau/BDF) |
| `--max-step` | float | auto | Maximum step size |
| `--dt` | float | CFL | Fixed time step (leapfrog only) |

### Default Tolerances

All adaptive solvers use `rtol=1e-8`, `atol=1e-10` by default. These are
production-quality defaults that bound the temporal integration error well
below spatial discretization error for typical grid resolutions (32-256
points).

### Validation Rules

- `--rtol` must be positive
- `--atol` must be positive
- `--max-step` must be positive
- `--method` is only meaningful with `--scheme cvode` or `--scheme scipy`

## Architecture

### Shared Infrastructure

All four solvers share core infrastructure, minimizing code duplication:

- **`compute_force()`** (`leapfrog.py`): Evaluates dpi/dt = F(q) for all
  momentum slots. Used by leapfrog, CVODE, and scipy.
- **`compute_velocity()`** (`leapfrog.py`): Evaluates dq/dt = K^{-1}(pi-S)
  for all field slots. Used by leapfrog, CVODE, and scipy.
- **`RHSEvaluator`** (`rhs.py`): Unified operator + coefficient application.
  Used by all four solvers.
- **`CoefficientEvaluator`** (`coefficients.py`): 4-level cache (L0-L3).
  Used by all adaptive solvers.
- **`build_jacobian_sparsity()`** (`sparsity.py`): Jacobian sparsity pattern
  for implicit solvers. Used by IDA, CVODE (sparse linsolver), and scipy
  (Radau/BDF).
- **`StateLayout`**, **`FieldSet`**, **`GridInfo`**: State management shared
  by all solvers.

### Result Format

All solvers return the same dictionary:

```python
{
    "t": np.ndarray,     # (num_snapshots,) time points
    "y": np.ndarray,     # (num_snapshots, total_size) state snapshots
    "success": bool,     # True if integration completed
    "message": str,      # Human-readable status
}
```

### Linear Solver Selection (SUNDIALS)

Both CVODE and IDA use the same thresholds for automatic linear solver
selection based on system size:

| System Size | Linear Solver | Rationale |
|-------------|---------------|-----------|
| N ≤ 2,000 | Dense | Direct solve, minimal overhead |
| 2,000 < N ≤ 200,000 | Sparse (SuperLU) | Exploits sparsity pattern |
| N > 200,000 | GMRES (iterative) | Avoids O(N²) memory |

### Sparsity Reuse

The IDA sparsity builder (`build_jacobian_sparsity()`) produces a pattern
for `J = dF/dy + cj * dF/dy'`. For CVODE (ODE form), the needed pattern
is `J = df/dy`. The IDA pattern is a **superset** of CVODE's true pattern:
extra diagonal entries from `cj * I` are always nonzero anyway. This means
the same sparsity builder works for both solvers — a key simplification.

For scipy's implicit methods (Radau, BDF), the same pattern is passed as
`jac_sparsity` to `solve_ivp`.

## Design Decisions

### Why CVODE as default for wave systems (not leapfrog)?

- **Tolerance-controlled:** Energy error bounded by rtol, not O(dt²)
- **Adaptive:** Automatically picks step sizes, no CFL tuning
- **Safe default:** Works correctly for all wave systems out of the box
- **Stiffness handling:** BDF handles mild stiffness from 1/dx² eigenvalues
- Leapfrog remains available for users who need symplectic integration

### Why CVODE over scipy as default?

- Method-of-lines PDE discretization is mildly stiff (eigenvalues near
  imaginary axis scale with 1/dx²). BDF handles this robustly.
- DOP853 would need a CFL-like max_step guard to prevent instability
  during quiescent phases.
- Same SUNDIALS ecosystem as IDA — consistent tolerance semantics.
- In-place RHS (`rhsfn(t, y, yp)`) avoids allocation per call.

### Why also provide scipy?

- DOP853 (8th-order explicit) excels for smooth non-stiff waves — fewer
  function evaluations than BDF for the same accuracy.
- One-step method = no startup cost (relevant for short simulations).
- Good option when implicit BDF is overkill.

### Why four solvers?

Each has a unique, irreplaceable strength:

| Solver | Unique Capability |
|--------|-------------------|
| Leapfrog | Symplectic — zero secular energy drift |
| CVODE | Implicit adaptive ODE — mildly stiff, tolerance-controlled |
| IDA | DAE — algebraic constraints (irreplaceable for Chern-Simons) |
| scipy | Explicit adaptive — DOP853 8th-order for smooth problems |

## Usage Examples

```bash
# Default — auto-selects best adaptive solver
tidal simulate spec.json

# Custom tolerances
tidal simulate spec.json --rtol 1e-6 --atol 1e-8

# CVODE with Adams (non-stiff, high-order)
tidal simulate spec.json --scheme cvode --method Adams

# High-accuracy explicit (DOP853)
tidal simulate spec.json --scheme scipy --rtol 1e-10

# IDA with custom tolerances
tidal simulate spec.json --scheme ida --rtol 1e-6

# Symplectic (opt-in, no tolerance control)
tidal simulate spec.json --scheme leapfrog
```

## References

- **Hindmarsh et al. (2005)**, "SUNDIALS: Suite of Nonlinear and
  Differential/Algebraic Equation Solvers", ACM TOMS 31(3). — CVODE and
  IDA algorithms, BDF/Adams methods, SUNDIALS architecture.

- **Hairer, Norsett & Wanner (1993)**, _Solving Ordinary Differential
  Equations I: Nonstiff Problems_, Springer. — DOP853, embedded Runge-Kutta
  theory, error estimation.

- **Hairer & Wanner (1996)**, _Solving Ordinary Differential Equations II:
  Stiff and Differential-Algebraic Problems_, Springer. — BDF, Radau,
  stiffness theory, DAE index theory.

- **Dormand & Prince (1980)**, "A family of embedded Runge-Kutta formulae",
  J. Comp. Appl. Math. 6, 19-26. — The RK method family underlying DOP853.

- **Hairer, Lubich & Wanner (2006)**, _Geometric Numerical Integration_,
  Springer. — Symplectic integrators (Stormer-Verlet/leapfrog).

- **Berlin et al. (2024)**, "Numerical analysis of resonant axion-photon
  mixing", arXiv:2405.08865. — Numerical challenges in resonant wave mixing
  analogous to TIDAL's Gertsenshtein simulation.

- **Courant, Friedrichs & Lewy (1928)**, "Über die partiellen
  Differenzengleichungen der mathematischen Physik", Math. Ann. 100, 32-74.
  — The CFL stability condition.

- **scikit-sundae** (NREL, BSD-3 license) — Python wrapper for SUNDIALS
  CVODE and IDA solvers.
