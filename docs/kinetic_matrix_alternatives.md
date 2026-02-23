# Non-Diagonal Kinetic Matrix: Alternative Approaches

When the kinetic matrix K\_{ij} = D[pi_i, vel_j] is non-diagonal (e.g., linearized
gravity, massive gravity), the standard Hamiltonian split dq/dt = pi fails. This
document records all investigated approaches and their tradeoffs.

## Current Implementation: Symbolic K^{-1} Inversion (Wolfram)

**Status:** Implemented and working (commit 906b250, hardened in cb7c12f).

Compute K symbolically in Wolfram, invert via `Inverse[K]`, emit K^{-1}\_{ij}
coefficients to JSON. Field rates become:

```
dq_i/dt = Sum_j K^{-1}_{ij} (pi_j - S_j)
```

where S_j = spatial part of pi_j (momenta evaluated with all velocities = 0).

**Pros:**

- Gauge-independent — works for any theory without gauge fixing
- Exact symbolic inversion — no numerical errors
- Pre-computed at derive-time — zero runtime cost
- Handles constant and parameter-dependent K (via coefficient_symbolic)

**Cons:**

- K^{-1} entries can be large rational expressions for big systems
- Requires det(K) != 0 (throws if singular, indicating constraint structure)
- _Position-dependent K (from curved metrics) would need grid evaluation at runtime_

**Works for:** Proca (2×2 diagonal), massive gravity (4×4), gravitational waves
(7×7 dynamical block), elasticity (2×2 diagonal with symbolic rho).

---

## Alternative 1: Mass Matrix ODE Solvers (Julia / PETSc / SUNDIALS)

**Status:** Investigated, not implemented. Viable future option.

Instead of inverting K, reformulate as a mass-matrix ODE:

```
M · dy/dt = f(y, t)
```

where M encodes the kinetic coupling. Specialised solvers handle M directly.

### Julia DifferentialEquations.jl

- `ODEProblem(M, f, u0, tspan)` with Rosenbrock or SDIRK methods
- Supports mass matrix M (constant or time-dependent)
- **Ref:** Rackauckas & Nie (2017), JORS 5(1), 15
- **How it would work:** Export K from Wolfram to JSON, construct M in Julia,
  use `Rodas5()` or `Rosenbrock23()` for time integration
- **Advantage:** Handles stiff systems, implicit methods, arbitrary M
- **Disadvantage:** Requires Julia runtime, not Python; would need new backend

### PETSc TS

- `DMTSCreateRHSMassMatrix()` for non-identity mass matrix
- Implicit Runge-Kutta, BDF, ARKIMEX schemes
- **Ref:** Balay et al. (2024), PETSc/TAO Users Manual
- **Advantage:** Industrial-strength, scalable, MPI-parallel
- **Disadvantage:** Heavy dependency, C API, complex setup

### SUNDIALS IDA

- General DAE solver: F(t, y, y') = 0 — avoids K^{-1} entirely
- Variable-order BDF with Newton iteration
- **Ref:** Hindmarsh et al. (2005), ACM TOMS 31(3)
- **Advantage:** Most general formulation, handles constraints naturally
- **Disadvantage:** Implicit (requires Jacobian), overkill for linear systems

### When to consider:

- Position-dependent K (curved spacetime with non-constant metric)
- Very large K (>50×50) where symbolic inversion produces huge expressions
- Stiff systems where explicit methods are CFL-limited
- If py-pde is replaced with a Julia or C++ backend

---

## Alternative 2: Gauge Fixing (Reduce K to Diagonal)

**Status:** TT gauge implemented (Phase B, uncommitted). De Donder planned.

Instead of inverting K, apply gauge conditions that eliminate off-diagonal
kinetic coupling, making K diagonal (or even identity).

### TT Gauge (Transverse-Traceless)

- Eliminates non-physical components entirely
- Reduces 10 GR components → 2 physical polarizations (h\_+, h_x)
- Type B constraint: temporal zeroing + traceless + transverse
- **Works for:** Gravitational waves in vacuum (massless)
- **Doesn't work for:** Massive gravity (need all components for mass term)

### De Donder (Harmonic) Gauge

- Adds L_gf = -(1/2ξ)(∂^a h_ab - ½ ∂_b h)² to Lagrangian
- Diagonalizes K and exposes constraint structure (h_0, h_0i become constraints)
- Type A (Lagrangian term): pipeline infrastructure exists, needs testing
- **Works for:** Both massless and massive gravity
- **Disadvantage:** Gauge-dependent; results must be checked for gauge artifacts

### When to consider:

- When physical interpretation of individual components is needed
- When reducing computational cost (fewer components to evolve)
- When connecting to observables in specific gauges (e.g., TT for GW detection)

---

## Alternative 3: Symplectic Integrators with Pre-Computed K^{-1}

**Status:** Investigated, not implemented.

For constant K (flat spacetime, constant parameters), pre-compute K^{-1} once
and use explicit symplectic integrator (leapfrog/Verlet) modified for K^{-1}:

```
q_{n+1} = q_n + dt * K^{-1} · p_{n+1/2}
p_{n+1} = p_n + dt * F(q_{n+1})
```

- **Ref:** Hairer & Lubich (2003), ZAMM 83(1)
- **Advantage:** Energy-conserving (symplectic), explicit, no Jacobian
- **Disadvantage:** Only works for constant K; our current explicit RK already works
- **When to consider:** Long-time simulations where energy conservation is critical

---

## Alternative 4: Runtime Numeric Inversion (NumPy)

**Status:** Not implemented. Simplest fallback.

Build K numerically from JSON coefficients at runtime, invert with `numpy.linalg.inv`.

- **Advantage:** Simple, no Wolfram dependency, handles parameter sweeps
- **Disadvantage:** Numerical (not exact), repeated inversion per timestep for
  position-dependent K, doesn't leverage symbolic structure
- **When to consider:** If symbolic K^{-1} expressions become too large

---

## Alternative 5: Mass Lumping (FEM Approximation)

**Status:** Not applicable to our spectral/FD discretization.

In FEM, the mass matrix M arises from basis function integration. "Mass lumping"
replaces M with a diagonal approximation (row-sum or higher-order lumping).

- **Ref:** Rathgeber et al. (2016), ACM TOMS 43(3)
- **Not applicable** because our K comes from the Lagrangian structure, not from
  spatial discretization. K is typically small (10×10 for GR) and exact inversion
  is feasible.

---

## Alternative 6: MultiModeCode-Style Field-Space Metric

**Status:** Investigated for reference.

In multi-field inflation, the kinetic matrix is the "field-space metric" G*{IJ}.
MultiModeCode evolves transport equations with non-trivial G*{IJ} directly.

- **Ref:** Dias et al. (2015), JCAP 2015(01):030
- **Relevance:** Similar mathematical structure to our K\_{ij}
- **Difference:** They evolve ODEs (homogeneous fields), we evolve PDEs (spatial grid)
- **When to consider:** If we ever add multi-field inflation examples

---

## Decision Matrix

| Approach                       | Implemented | Exact | Runtime Cost | Handles Position-Dep K | Handles DAEs |
| ------------------------------ | :---------: | :---: | :----------: | :--------------------: | :----------: |
| Symbolic K^{-1} (current)      |     ✅      |  ✅   |     Zero     |       Partial\*        |      ❌      |
| Julia DifferentialEquations.jl |     ❌      |  ❌   |    Medium    |           ✅           |      ✅      |
| PETSc TS                       |     ❌      |  ❌   |     Low      |           ✅           |      ✅      |
| SUNDIALS IDA                   |     ❌      |  ❌   |    Medium    |           ✅           |      ✅      |
| Gauge fixing                   |   Partial   |  ✅   |     Zero     |          N/A           |     N/A      |
| Symplectic K^{-1}              |     ❌      |  ✅   |     Low      |           ❌           |      ❌      |
| NumPy inversion                |     ❌      |  ❌   |     High     |           ✅           |      ❌      |

\*Position-dependent K: coefficient_symbolic + evaluate_coefficient() supports it,
but not yet tested with non-constant K^{-1} entries.
