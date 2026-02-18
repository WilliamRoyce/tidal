# Phase F: Adaptive Time-Stepping

**Status:** Implemented
**Date:** 2026-02-18
**Implementation:** ~300 new lines across 4 files, 27 new tests

Adaptive time-stepping replaces hardcoded fixed time steps with
error-controlled solvers that automatically select step sizes. This gives
users accuracy guarantees (via tolerances) instead of hoping a manually
chosen `dt` is "small enough".

## Problem Statement

Before Phase F, every TIDAL simulation used a fixed time step:

```python
DT = 0.01  # manually chosen, hope it's stable
pde.solve(state, t_range=T_END, dt=DT, scheme="runge-kutta")
```

This has three fundamental problems:

1. **Fragility.** A `dt` that works for one parameter set can be unstable
   for another. The user has to know the CFL condition and recompute `dt`
   whenever they change grid resolution or wave speed. Getting it wrong
   produces garbage output that can look plausible.

2. **Inefficiency.** A `dt` chosen for the fastest mode wastes compute on
   slow phases of the simulation. Resonant mixing (the Gertsenshtein
   effect) has beat frequencies much slower than the carrier — fixed steps
   resolve every carrier cycle even when nothing interesting is happening.
   Berlin et al. (2024, arXiv:2405.08865) document these numerical
   challenges in resonant axion-photon mixing.

3. **Stiffness.** When mass terms satisfy m^2 dx^2 >> 1, explicit methods
   require impossibly small steps. Implicit methods (Radau, BDF) solve
   this, but need Jacobian information to be efficient.

## Design Goals

Phase F adds three solver paths while maintaining the project's
zero-new-dependency principle (scipy is already a transitive dependency
of py-pde):

- **Fixed-step explicit**: unchanged default for backwards compatibility
- **Adaptive explicit**: py-pde's native RKF error estimator
- **Adaptive scipy**: full `solve_ivp` arsenal (DOP853, Radau, BDF, etc.)

Supporting infrastructure:

- CFL auto-computation as safety guard and initial step hint
- Jacobian sparsity pattern for implicit solver performance
- Energy conservation monitor (opt-in) and blow-up detection (always on)
- Automatic stiffness advisory

## Three Solver Paths

### Path A: Fixed-Step Explicit (default)

```bash
tidal simulate spec.json --dt 0.01
tidal simulate spec.json                  # dt auto-computed from CFL
```

The original behaviour. Uses py-pde's `ExplicitSolver` with constant `dt`.
CFL stability is checked at startup; violations produce warnings on stderr.

When `--dt` is omitted, the CLI computes `dt = 0.5 * min(dx)` as a
conservative default. The `0.5` safety factor (`_CFL_FACTOR`) provides
margin below the CFL limit.

**Best for:** Quick exploratory runs, debugging, known-stable configurations.

### Path B: Adaptive Explicit

```bash
tidal simulate spec.json --scheme runge-kutta --adaptive
tidal simulate spec.json --scheme runge-kutta --adaptive --tolerance 1e-6
```

Uses py-pde's built-in Runge-Kutta-Fehlberg error estimator. The solver
compares a 4th-order and 5th-order solution at each step and adjusts
`dt` to keep the local truncation error below `tolerance`.

The `--dt` value (or CFL auto-dt) serves as the initial step guess.
The `--tolerance` flag controls accuracy (default: py-pde's 1e-4).

**Best for:** Systems with moderate scale separation where explicit
methods are stable but fixed steps waste compute.

### Path C: Adaptive Scipy (recommended for production)

```bash
tidal simulate spec.json --scheme scipy
tidal simulate spec.json --scheme scipy --method DOP853
tidal simulate spec.json --scheme scipy --method DOP853 --rtol 1e-8 --atol 1e-10
tidal simulate spec.json --scheme scipy --method Radau   # implicit, for stiff systems
```

Forwards all keyword arguments to `scipy.integrate.solve_ivp` via
py-pde's `ScipySolver(**kwargs)` passthrough. Default method is DOP853
(8th-order Dormand-Prince), an explicit one-step method with dense output.

The CFL limit is automatically used as `max_step` to prevent the adaptive
solver from overshooting the spatial stability boundary. For implicit
methods (Radau, BDF), the Jacobian sparsity pattern is computed and
passed to `solve_ivp` for efficient Newton iterations.

**Best for:** Production runs, stiff systems, long simulations, and any
case where accuracy matters more than speed.

### Choosing a Path

| Scenario | Recommended Path |
|----------|-----------------|
| Quick test, debugging | Path A (fixed) |
| Standard wave equation, moderate accuracy | Path C (`--scheme scipy`) |
| High-precision convergence study | Path C with `--rtol 1e-10 --atol 1e-12` |
| Stiff system (m^2 dx^2 >> 1) | Path C with `--method Radau` |
| Multi-step implicit (stiff + history) | Path C with `--method BDF` |
| Backward compatibility with old scripts | Path A with explicit `--dt` |

## CLI Reference

All new flags are on the `tidal simulate` subcommand. Defined in
`tidal/cli/__init__.py:206-267`.

| Flag | Type | Default | Requires | Purpose |
|------|------|---------|----------|---------|
| `--scheme` | `{runge-kutta, scipy}` | `runge-kutta` | — | Solver backend selection |
| `--adaptive` | bool flag | `False` | `--scheme runge-kutta` | Enable adaptive step-size control |
| `--method` | `{RK45, RK23, DOP853, Radau, BDF, LSODA}` | `DOP853` | `--scheme scipy` | scipy `solve_ivp` integration method |
| `--rtol` | float | scipy default (1e-3) | `--scheme scipy` | Relative error tolerance |
| `--atol` | float | scipy default (1e-6) | `--scheme scipy` | Absolute error tolerance |
| `--tolerance` | float | py-pde default (1e-4) | `--adaptive` | Error tolerance for adaptive explicit RK |
| `--max-step` | float | auto (CFL) | any | Maximum allowed step size |
| `--energy-monitor` | float | disabled | any | Halt if `\|dE/E0\|` exceeds threshold |

### Validation Rules

Validated in `_validate_solver_params()` at `tidal/cli/_simulate.py:756-804`:

- `--method`, `--rtol`, `--atol` require `--scheme scipy`
- `--tolerance` requires `--scheme runge-kutta --adaptive`
- All tolerance/step values must be positive
- `--energy-monitor` threshold must be positive

## Implementation Details

### CFL Limit Computation

**Location:** `PDEFromSpec.cfl_limit(grid)` at `tidal/symbolic/pde_builder.py:2968-3011`

Scans all second-order-in-time equations for Laplacian and biharmonic
operators. Computes:

- **Wave CFL:** `dt_max = dx_min / sqrt(max_laplacian_coeff)` — the
  classical Courant condition where `sqrt(coeff)` is the wave speed.
- **Biharmonic stability:** `dt_max = dx_min^4 / (2 * max_biharmonic_coeff)`

Returns the most restrictive limit, or `None` if no wave/biharmonic terms
exist (e.g., constraint-only or first-order systems).

Cross-field Laplacian terms (e.g., `laplacian(chi)` in the `phi` equation)
are included — omitting them would underestimate the maximum wave speed.
Directional Laplacians (`laplacian_x`, `laplacian_y`) are also included.

**Note on `time_derivative_order`:** The JSON stores the **original
Lagrangian derivative order**, not the internal first-order form. Wave
equations from the Lagrangian get `time_derivative_order = 2` in JSON
(set by `DetectLHSTimeOrder` in `ExportJSON.wl`), even though
`pde_builder.py` internally decomposes them into first-order (field,
momentum) pairs for time evolution. Constraints get order 0, true
parabolic equations get order 1. The `< 2` check correctly identifies
wave equations that need CFL checking — it does not skip them despite the
internal first-order representation.

### Stability Checking

**Location:** `PDEFromSpec.check_stability(dt, grid)` at `tidal/symbolic/pde_builder.py:2909-2966`

Same algorithm as `cfl_limit()`, but compares a proposed `dt` against
each limit and returns a list of human-readable warning strings. Used by
Path A (fixed-step) to print warnings at simulation startup.

### Jacobian Sparsity Pattern

**Location:** `PDEFromSpec.jacobian_sparsity(grid)` at `tidal/symbolic/pde_builder.py:3013-3121`

For implicit solvers (Radau, BDF), providing the Jacobian sparsity pattern
avoids O(N) finite-difference evaluations per Newton step. Without it,
`solve_ivp` estimates the full N x N Jacobian by perturbing each state
variable — prohibitively expensive for large grids.

The sparsity pattern is a binary matrix indicating which state variables
influence which rates. TIDAL computes it from the equation structure:

**State vector layout** (for a system with fields phi, chi, each second-order):
```
[ phi_0, phi_1, ..., phi_{N-1},    # field values
  pi_phi_0, ..., pi_phi_{N-1},     # momenta (d/dt phi)
  chi_0, ..., chi_{N-1},           # field values
  pi_chi_0, ..., pi_chi_{N-1} ]   # momenta
```

**Sparsity rules:**
- `d/dt field = momentum` — diagonal block (field ← momentum)
- `d/dt momentum = laplacian(field)` — tridiagonal block (3-point stencil)
- `d/dt momentum = identity(other_field)` — diagonal block (cross-field)
- `d/dt momentum = gradient(other_field)` — tridiagonal block (cross-field)
- Periodic BCs add corner entries (first ↔ last grid point wrapping)

**Limitations:** Only 1D grids are supported. For multi-dimensional grids
(2D, 3D), the Laplacian stencil involves neighbours at offsets +/-1,
+/-N_x, +/-N_x*N_y, etc. — requiring a multi-diagonal band structure
that the current implementation does not construct. Multi-D grids return
`None`, causing `solve_ivp` to fall back to dense Jacobian estimation.

### Stiffness Advisory

**Location:** `_stiffness_advisory()` at `tidal/cli/_simulate.py:807-853`

Automatically warns when the system appears stiff and an explicit method
is selected. Computes a dimensionless stiffness ratio:

```
stiffness_ratio = max(|m^2|) * dx_min^2 / max(|c^2|)
```

where `m^2` is the largest mass matrix entry and `c^2` is the largest
Laplacian coefficient (wave speed squared). This ratio measures how many
mass oscillation periods fit within one spatial CFL step. When it
exceeds 100, explicit methods waste many steps resolving mass oscillations.

The advisory prints:
```
Note: system may be stiff (m^2*dx^2/c^2=450). Consider --scheme scipy --method Radau.
```

The system does NOT automatically switch methods — the user retains full
control over solver selection.

### Energy Conservation Monitor

**Location:** `_build_energy_monitor()` at `tidal/cli/_simulate.py:710-744`

An opt-in callback (`--energy-monitor THRESHOLD`) that tracks the L2 norm
of the full state vector as an energy proxy:

```python
E(t) = sum(state.data^2)   # L2 norm squared
```

At each snapshot, it computes the relative drift `|E(t) - E(0)| / E(0)`.
If this exceeds the threshold, the simulation halts with `RuntimeError`.

The L2 norm is conserved (up to discretization error) for Hamiltonian
systems with symplectic-like time stepping. It is NOT the true physical
Hamiltonian energy — but for detecting numerical instability (the primary
use case), any conserved-ish quantity suffices.

### Blow-Up Detection

**Location:** `simulate_command()` at `tidal/cli/_simulate.py:974-993`

Always-on callback that checks `max(|state|)` at each snapshot. If any
field component exceeds `max(initial_max * 1e6, 1e6)`, the simulation
halts. This catches tachyonic instabilities (m^2 < 0) and numerical
explosions early, before they produce multi-GB output files full of `inf`.

The threshold is 6 orders of magnitude above the initial state or 10^6,
whichever is larger. The absolute floor of 10^6 handles the case where
the initial state is zero (e.g., chi in a phi-only IC).

## py-pde ScipySolver Architecture

Understanding how py-pde wraps scipy is important for choosing methods.

### "Many Calls" Architecture (current)

py-pde's `Controller` calls `solve_ivp` **once per tracker interval**
(snapshot interval). Each call is independent:

```python
# Pseudocode for py-pde's Controller loop (pde/solvers/scipy.py:87-109)
for t_start, t_end in snapshot_intervals:
    solve_ivp(rhs, t_span=(t_start, t_end), y0=state, **kwargs)
    # state updated, trackers called, next interval
```

### Impact on Method Choice

| Method | Type | Restart Penalty | Recommendation |
|--------|------|----------------|----------------|
| DOP853 | Explicit, one-step, order 8 | None | Default choice |
| RK45 | Explicit, one-step, order 5 | None | Faster per step, lower order |
| RK23 | Explicit, one-step, order 3 | None | Rough estimates only |
| Radau | Implicit, one-step, order 5 | Minimal | Best for stiff systems |
| BDF | Implicit, multi-step, order 1-5 | **Significant** — restarts at order 1, loses history | Penalized by many-calls architecture |
| LSODA | Auto-switch Adams/BDF | **Significant** — loses stiff/non-stiff detection | Penalized by many-calls architecture |

**Why DOP853 is the default:** It is an 8th-order one-step method with
excellent accuracy per function evaluation. Being one-step, it suffers
no penalty from py-pde's per-interval restarts. The 8th-order accuracy
means it can take much larger steps than RK45 for smooth wave problems.

**For stiff systems:** Use Radau (implicit one-step). It handles stiffness
well and has minimal restart penalty. BDF is theoretically better for
sustained stiff integration, but the restart penalty from py-pde's
architecture degrades it significantly.

A "single call" bypass (calling `solve_ivp` once for the entire simulation
with `t_eval` for snapshots) would eliminate the BDF restart penalty. This
is deferred — see Future Extensions.

## Example Migration

Phase F also updated all example scripts to use adaptive time-stepping
as the recommended default.

### Before (fixed-step)

**run.sh:**
```bash
tidal simulate ../data/coupled_scalars.json \
  --dt 0.01 \
  --grid-shape 256 --bounds 0:100 --periodic \
  --t-end 20.0
```

**Python:**
```python
DT = 0.01
pde.solve(initial_state, t_range=T_END, dt=DT, scheme="runge-kutta", tracker=tracker)
```

### After (adaptive scipy)

**run.sh:**
```bash
tidal simulate ../data/coupled_scalars.json \
  --scheme scipy \
  --grid-shape 256 --bounds 0:100 --periodic \
  --t-end 20.0
```

**Python:**
```python
pde.solve(initial_state, t_range=T_END, solver="scipy", method="DOP853", tracker=tracker)
```

**Files updated:** 21 run.sh scripts (removed `--dt`, added `--scheme scipy`)
and 24 Python scripts (removed `DT` constant, replaced solver call).

## Test Coverage

27 new tests across three categories.

### Argument Parsing & Validation (17 tests)

Located in `tests/test_cli.py`. Cover all flag combinations, type
checking, and mutual-exclusion rules:

- `--method` without `--scheme scipy` → `ValueError`
- `--rtol` without `--scheme scipy` → `ValueError`
- `--tolerance` without `--adaptive` → `ValueError`
- Negative tolerance values → `ValueError`
- Valid combinations parse correctly

### PDE Builder Methods (7 tests)

Located in `tests/test_pde_builder.py`:

| Test | What it verifies |
|------|-----------------|
| `test_wave_equation_cfl` | CFL = dx/c for unit wave speed |
| `test_massive_wave_cfl` | CFL ignores mass term (identity operator) |
| `test_no_wave_returns_none` | First-order equations → `None` |
| `test_two_field_picks_smallest` | Most restrictive limit across fields |
| `test_directional_laplacian_cfl` | `laplacian_x` contributes to CFL |
| `test_cross_field_laplacian_cfl` | Cross-field Laplacian included |
| `test_multi_d_returns_none` | Jacobian sparsity returns `None` for 2D+ |

### Integration Tests (3 tests)

Located in `tests/test_cli.py`. Run actual simulations:

| Test | What it verifies |
|------|-----------------|
| `test_simulate_scipy_adaptive_dop853` | Scipy DOP853 with custom rtol/atol |
| `test_simulate_adaptive_rk` | py-pde adaptive explicit RK |
| `test_simulate_scipy_with_max_step` | Scipy with max-step limit |

### Monitoring Tests (2 tests)

| Test | What it verifies |
|------|-----------------|
| `test_energy_monitor_fires` | Energy monitor triggers on impossible threshold |
| `test_blowup_detection_fires` | Blow-up detection catches tachyonic instability |

## Design Decisions & Rationale

### Why DOP853 as default (not RK45)?

DOP853 (Dormand-Prince order 8) delivers ~3 orders of magnitude better
accuracy per step than RK45 (order 5) for smooth problems. Wave equations
with polynomial ICs are smooth, so the higher order pays off. The cost
per step is ~13 function evaluations (vs ~6 for RK45), but the larger
permissible step sizes more than compensate.

For the Gertsenshtein problem specifically, the resonant beat frequency
creates smooth envelopes that DOP853 captures with far fewer steps than
lower-order methods.

### Why CFL as `max_step` guard?

Even adaptive methods can overshoot. If the solver's error estimate
happens to be small (e.g., during a quiescent phase), it may try a step
larger than the spatial stability limit. Setting `max_step = CFL_limit`
prevents this. The adaptive solver can take smaller steps than `max_step`
but never larger.

### Why Jacobian sparsity for 1D only?

The 3-point Laplacian stencil on a 1D grid produces a tridiagonal band
in the Jacobian. The pattern is simple: neighbours are at indices i-1,
i, i+1 (plus periodic wrapping). For a 2D grid of shape (Nx, Ny), the
stencil involves neighbours at offsets +/-1 and +/-Nx in the flattened
vector — a penta-diagonal band. For 3D, it's hepta-diagonal. Each
dimension adds two more diagonals.

Constructing the correct multi-diagonal band structure requires knowledge
of the grid shape and periodic boundary conditions in each dimension.
This is straightforward but was deferred to keep Phase F focused. For
multi-D, returning `None` causes `solve_ivp` to use dense Jacobian
estimation, which is slower but correct.

### Why L2 norm as energy proxy (not true Hamiltonian)?

The true Hamiltonian for a coupled wave system is:

```
H = sum_i integral[ (1/2)(pi_i)^2 + (1/2)(grad phi_i)^2 + V(phi) ] dx
```

Computing this requires knowledge of the potential V(phi), gradient
operators, and integration weights — all equation-specific. The L2 norm
`sum(state^2)` is equation-agnostic and still detects numerical instability
(the primary use case: if the solver diverges, L2 norm grows rapidly).

For true energy conservation checks, the measurement module's
`compute_system_energy()` provides the physical Hamiltonian — but it
operates post-hoc on saved data, not during simulation.

### Why blow-up threshold = max(initial * 1e6, 1e6)?

The threshold needs to be:
- **High enough** to not trigger on legitimate physics (large but finite amplification)
- **Low enough** to catch actual divergence before it produces inf/nan
- **Nonzero** even when the initial state is zero (e.g., chi starts at 0)

Six orders of magnitude covers the range between "largest physically
plausible amplification" (resonant amplification rarely exceeds 10^3)
and "definitely diverging" (10^6 growth in finite time is not physical).
The absolute floor of 10^6 handles zero-initialized fields.

### Why stiffness threshold = 100?

The ratio `m^2 dx^2 / c^2` counts how many mass oscillation periods per
CFL step. When this exceeds ~10, explicit methods start wasting steps.
At 100, the waste factor is severe enough to warrant a user-visible
advisory. The threshold is conservative — some systems are manageable at
ratios up to ~50 — but false positives (advisory when not needed) are
harmless, while false negatives (no advisory for a stiff system) waste
the user's time.

### Why no automatic method selection?

Automatic method selection (e.g., "if stiff, use Radau") would require
reliable stiffness detection, which is an open research problem. The
stiffness ratio heuristic is a rough indicator — it can miss systems
that are stiff due to boundary conditions or nonlinear interactions.
Giving the user an advisory + control is safer than silent auto-switching.

## Known Limitations

1. **Jacobian sparsity: 1D only.** Multi-dimensional grids fall back to
   dense Jacobian estimation, which is O(N) per Newton step. For large 2D/3D
   grids with implicit methods, this can be prohibitively slow.

2. **BDF restart penalty.** py-pde's "many calls" architecture restarts
   BDF at order 1 each snapshot interval. For stiff systems that benefit
   from high-order BDF, this negates much of the advantage. Use Radau
   (implicit one-step) instead until the single-call bypass is implemented.

3. **Energy monitor is approximate.** The L2 norm proxy does not equal the
   true physical Hamiltonian. It can drift for non-symplectic integration
   even when the system is stable (e.g., BDF introduces numerical
   dissipation). Use the measurement module for post-hoc energy analysis.

4. **No adaptive mesh refinement.** Time-step adaptation only. Spatial
   resolution is fixed throughout the simulation. AMR would require
   significant py-pde extensions or a solver migration.

5. **Position-dependent coefficients.** The CFL and stiffness estimates
   use the constant coefficient from the JSON spec. For position-dependent
   coefficients (background fields), the actual wave speed varies across
   the grid. The estimates use the symbolic constant, which may over- or
   under-estimate the true stability limit.

## Future Extensions

1. **Single-call `solve_ivp` bypass.** Call `solve_ivp` once for the
   entire simulation with `t_eval=np.linspace(0, T_END, n_snapshots)`.
   Requires ~100 lines of custom snapshot collection, bypassing py-pde's
   tracker system. Eliminates BDF restart penalty.

2. **Multi-D Jacobian sparsity.** Construct the correct banded pattern
   for 2D/3D grids (penta-diagonal for 2D, hepta-diagonal for 3D, plus
   periodic wrapping). Requires knowledge of grid shape per dimension.

3. **True Hamiltonian energy monitor.** Compute the physical Hamiltonian
   during simulation using the equation structure. Requires integrating
   `compute_system_energy()` into the callback — currently it operates
   on `SimulationData`, not raw state vectors.

4. **Adaptive mesh refinement.** Refine spatial resolution near wavefronts
   or coupling regions. Would require migrating from py-pde's uniform
   Cartesian grids to an AMR-capable backend (e.g., PETSc, AMReX).

## References

- **Hairer, Norsett & Wanner (1993)**, *Solving Ordinary Differential
  Equations I: Nonstiff Problems*, Springer. — DOP853 method, embedded
  Runge-Kutta theory, error estimation.

- **Hairer & Wanner (1996)**, *Solving Ordinary Differential Equations II:
  Stiff and Differential-Algebraic Problems*, Springer. — BDF, Radau,
  stiffness theory, Jacobian sparsity exploitation.

- **Dormand & Prince (1980)**, "A family of embedded Runge-Kutta formulae",
  J. Comp. Appl. Math. 6, 19-26. — The RK method family underlying
  RK45 and DOP853.

- **Berlin et al. (2024)**, "Numerical analysis of resonant axion-photon
  mixing", arXiv:2405.08865. — Documents numerical challenges in resonant
  wave mixing directly analogous to TIDAL's Gertsenshtein simulation.

- **Zwicker (2020)**, "py-pde: A Python package for solving partial
  differential equations", JOSS 5(48), 2158. — The PDE backend used by
  TIDAL; ScipySolver architecture.

- **Courant, Friedrichs & Lewy (1928)**, "Uber die partiellen
  Differenzengleichungen der mathematischen Physik", Math. Ann. 100,
  32-74. — The CFL stability condition.
