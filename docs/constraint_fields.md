# Constraint Fields in TIDAL

**Status:** Active
**Date:** 2026-02-19
**Applies to:** Any Lagrangian that produces mixed time-derivative orders

---

## Overview

When the Wolfram pipeline derives equations of motion via `VarD` (Euler-Lagrange
variation), the resulting component equations may have different time-derivative
orders on the LHS:

| `time_derivative_order` | Type                 | Mathematical character             | Solver treatment                         |
| ----------------------- | -------------------- | ---------------------------------- | ---------------------------------------- |
| 2                       | **Dynamical** (wave) | Hyperbolic: `∂²_t φ = RHS`         | ODE-integrated (Runge-Kutta / DOP853)    |
| 1                       | **First-order**      | Parabolic/advective: `∂_t φ = RHS` | ODE-integrated                           |
| 0                       | **Constraint**       | Elliptic: `0 = RHS`                | Solved at each instant (FFT / iterative) |

Constraint fields are **not** integrated forward in time. They are determined
by an elliptic equation at each instant, given the current state of the
dynamical fields. This is a fundamental feature of the Lagrangian, not a
pipeline limitation.

---

## Physics: Why Constraint Equations Arise

### The Proca Example

Consider the Proca Lagrangian for a massive vector field `A_μ` in 2+1D with
Minkowski signature `(-,+,+)`:

```
L = -1/2 (∂_a A_b) η^{ac} η^{bd} (∂_c A_d)
  + 1/2 (∂_a A_b) η^{ac} η^{bd} (∂_d A_c)
  - m²/2 A_a η^{ab} A_b
```

The Euler-Lagrange equation `∂_b F^{ba} + m² A^a = 0` yields three component
equations. For `a = 0` (the temporal component):

```
∂_x F^{x0} + ∂_y F^{y0} + m² A^0 = 0
```

Expanding the field strength `F^{i0} = η^{ii} η^{00} (∂_i A_0 - ∂_t A_i)`:

```
-∂_x(∂_x A_0 - ∂_t A_1) - ∂_y(∂_y A_0 - ∂_t A_2) - m² A_0 = 0
```

Rearranging:

```
∇² A_0 - m² A_0 = -∂_x(∂_t A_1) - ∂_y(∂_t A_2)
```

**There is no `∂²_t A_0` term.** The metric signature `η^{00} = -1` causes
the temporal kinetic term to vanish from the `a = 0` equation. The result is a
**Helmholtz equation** — an elliptic PDE that determines `A_0` at each instant
from the momenta of the spatial components `A_1`, `A_2`.

### Analogy: Coulomb's Equation in Electrostatics

This is analogous to the electric potential `V` in electrostatics. Poisson's
equation `∇²V = -ρ/ε₀` determines `V` instantaneously from the charge density
`ρ`. If `ρ` changes (charges move), `V` changes — but `V` has no independent
dynamics. It is a "slave" to the dynamical degrees of freedom.

The same pattern appears in:

- **Numerical relativity:** The Hamiltonian and momentum constraints determine
  the lapse and shift (Baumgarte & Shapiro, _Numerical Relativity_, Cambridge,
  2010, Chapter 3).
- **Incompressible fluid dynamics:** The pressure `p` satisfies a Poisson
  equation `∇²p = -ρ ∂_i u_j ∂_j u_i` and is determined at each instant from
  the velocity field (Chorin, 1968; Temam, 1969).
- **Gauge field theories:** The temporal component `A_0` of any gauge field
  satisfies a constraint (Gauss's law), not a wave equation.

### General Rule

For any Lagrangian with Lorentzian signature, the EOM for the **temporal
component** of a vector or tensor field will typically have
`time_derivative_order = 0`. This is not an error — it is a direct consequence
of the signature. The pipeline correctly identifies and classifies these
equations.

---

## Solver Architecture

### State Vector Layout

Constraint fields are included in the ODE state vector so that dynamical
equations can reference them. The state layout for coupled Proca
(`A_0, A_1, π_1, A_2, π_2, B_0, B_1, π_4, B_2, π_5`) includes both
constraint fields (`A_0`, `B_0`) and dynamical fields with their momenta.

### Evolution Rate

At each `evolution_rate(state, t)` call (`pde_builder.py`):

1. **Pass 1a — Constraint solve:** `_evolve_constraints()` solves the elliptic
   equations for all constraint fields, updating their values in-place.
   It then estimates `∂_t(constraint_field)` via finite difference (see below).
2. **Pass 1b — First-order momenta:** `_evolve_first_order()` computes virtual
   momenta for `time_derivative_order = 1` fields.
3. **Pass 2 — Rates:** `_evolve_second_order()` builds the full rate vector.
   For constraint fields, the rate is zero (the integrator should not modify
   them — the constraint solver handles that).

### The ODE rate for constraint fields is zero

```python
# pde_builder.py, _evolve_second_order():
if eq.time_derivative_order >= 2:
    rates[field_slot] = momentum.copy()
    rates[momentum_slot] = compute_rhs(...)
elif eq.time_derivative_order == 1:
    rates[field_slot] = virtual_momenta[field_name]
else:
    rates[field_slot] = ScalarField(grid, data=0.0)  # constraint
```

The zero rate tells the ODE integrator "don't change this field." The
constraint solver overrides the field value at the start of the next
`evolution_rate()` call.

---

## Constraint Momentum Estimation

### The Problem

Although constraint fields have no evolution equation, they **do change in
time** because their source terms (dynamical field momenta) evolve. If a
dynamical equation references `gradient_x(pi_0)` — i.e., `∂_x(∂_t A_0)` —
this term should be nonzero.

Prior to this fix, the code set `virtual_momenta["A_0"] = 0` unconditionally,
dropping the `∂_x(∂_t A_0)` coupling from the dynamical equations.

### The Fix: Backward Finite Difference

After solving the constraint at time `t`, we compare the new field value with
the stored value from the previous `evolution_rate()` call:

```
π_N(t) ≈ [field_N(t) − field_N(t_prev)] / (t − t_prev)
```

This is a **backward (causal) finite difference**, first-order accurate in
`Δt`. On the first call (`t = 0` or no previous data), we return `π_N = 0`.

### Why This Works Generally

The finite-difference approach is **theory-agnostic**:

- No knowledge of the constraint equation form is needed
- No Wolfram pipeline changes required
- No gauge choices — it is purely kinematic (observe the field changed,
  estimate the rate)
- Works for scalar, vector, or tensor constraint fields
- Works for coupled constraints (e.g., `A_0 ↔ B_0` in coupled Proca)
- Works with any adaptive ODE solver (DOP853, Radau, BDF)

### Implementation

Persistent state on `PDEFromSpec`:

```python
self._prev_constraint_fields: dict[str, NumericArray] = {}
self._prev_constraint_t: float | None = None
```

Method `_estimate_constraint_momenta()` runs after the constraint solve,
updating `virtual_momenta` in-place.

### Accuracy and Limitations

- **First-order in Δt:** The error in `π_N` is `O(Δt)`, where `Δt` is the
  time between consecutive `evolution_rate()` calls. For DOP853 (13 stages
  per step), these calls are closely spaced within each step, so the
  approximation is generally accurate.
- **First call:** At `t = 0`, `π_N = 0`. This is correct when initial
  conditions have zero momenta (constraint sources are zero).
- **Adaptive step rejection:** If the solver retries a step, the stored
  previous state may come from a rejected stage. The finite-difference
  estimate remains valid because the state and time are self-consistent.

### Future Enhancement

For higher accuracy, the constraint equation could be differentiated in time
analytically, yielding another elliptic equation for `π_N` with exact source
terms. This would eliminate the `O(Δt)` error but requires additional
elliptic solves. See Phase B (gauge fixing) for related work.

---

## References

1. Baumgarte, T. W. & Shapiro, S. L. (2010). _Numerical Relativity: Solving
   Einstein's Equations on the Computer_. Cambridge University Press.
   — Chapters 2-3: constraint equations in 3+1 decomposition.

2. Chorin, A. J. (1968). "Numerical solution of the Navier-Stokes equations."
   _Mathematics of Computation_, 22(104), 745-762.
   — Projection method: pressure as an elliptic constraint.

3. Jackson, J. D. (1999). _Classical Electrodynamics_ (3rd ed.). Wiley.
   — Section 6.3: Coulomb gauge and the scalar potential as a constraint.

4. Ruehl, W. & Yurov, V. P. (2004). "Proca equations derived from first
   principles." _arXiv:hep-th/0412059_.
   — Structure of Proca EOM and constraint analysis.
