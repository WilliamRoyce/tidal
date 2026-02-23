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
| 2                       | **Dynamical** (wave) | Hyperbolic: `∂²_t φ = RHS`         | ODE-integrated (CVODE / IDA / scipy)     |
| 1                       | **First-order**      | Parabolic/advective: `∂_t φ = RHS` | ODE-integrated                           |
| 0                       | **Constraint**       | Elliptic: `0 = RHS`                | Algebraic in IDA; pre-solved (FFT / sparse) |

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

### Solver Treatment

**With IDA (DAE solver — default for systems with constraints):**

Constraint equations are treated as algebraic equations within the DAE
system. IDA's Newton iteration simultaneously solves for all field values
(dynamical + constraint) at each timestep. The constraint pre-solve module
(`tidal/solver/constraint_solve.py`) ensures consistent initial conditions
before IDA starts. See `tidal/solver/ida.py` for the residual construction.

**Residual structure:** For constraint fields, the IDA residual is:
```
F_i = RHS_i(state)   (algebraic: no time derivative on LHS)
```
For dynamical fields:
```
F_i = yp_i - RHS_i(state)   (differential: yp = dy/dt)
```

**Constraint pre-solve** (`tidal/solver/constraint_solve.py`): Three-tier
solver runs before IDA to find constraint field values consistent with the
initial dynamical field state:
1. **Tier 1 (FFT)** — periodic BCs, constant coefficients (O(N log N))
2. **Tier 2 (Sparse probe)** — non-periodic BCs or variable coefficients (O(N²) build, O(N) solve)
3. **Automatic selection** — `_select_method()` chooses the fastest applicable tier

**Gauge regularisation** for singular Poisson with periodic BCs: the zero
Fourier mode is pinned (`u_hat[0,...,0] = 0`) to fix the gauge freedom.
This is numerical regularisation, not physics gauge fixing — observables
depend on derivatives of the constraint field, not its absolute value.

---

## Constraint Momentum Estimation

### The Problem

Although constraint fields have no evolution equation, they **do change in
time** because their source terms (dynamical field momenta) evolve. If a
dynamical equation references `gradient_x(pi_0)` — i.e., `∂_x(∂_t A_0)` —
this term should be nonzero.

Prior to this fix, the code set `virtual_momenta["A_0"] = 0` unconditionally,
dropping the `∂_x(∂_t A_0)` coupling from the dynamical equations.

### The Fix: Analytical Constraint Momentum via Time Differentiation

Rather than estimating `π_N` via finite differences (which introduces a time
lag that causes k-dependent instability through the constraint-dynamical
feedback loop), we compute `π_N` **exactly** by differentiating the constraint
equation in time.

#### Mathematical Derivation

Given a constraint equation of the form:

```
L[A_0] = -S(π_1, π_2, ..., B_0, ...)
```

where `L` is the self-operator (e.g., `∇² − m²`) and `S` is the source
(terms referencing other fields and their momenta), differentiating both sides
with respect to time yields:

```
L[∂_t A_0] = L[π_0] = -∂_t S
```

This is another elliptic equation for `π_0`, with the **same self-operator**
`L` as the original constraint, and a new source `-∂_t S` computed from the
time derivatives of the dynamical fields.

#### Computing ∂_t S

The time derivative of each source term is obtained by replacing fields with
their time derivatives:

| Source term type | Replacement | Source of replacement |
| --- | --- | --- |
| `c · op(π_j)` where j is 2nd-order | `c · op(∂_t π_j) = c · op(RHS_j)` | Dynamical equation j |
| `c · op(h_j)` where h_j is 2nd-order field | `c · op(π_j)` | Momentum from state |
| `c · op(constraint_field_j)` | `c · op(π_j)` | Off-diagonal coupling |

#### Handling Back-Coupling

The dynamical equation `RHS_j` may itself contain `gradient(π_0)` terms
(the very momentum we are solving for). These terms are moved to the LHS,
creating an **augmented operator**:

```
L_aug[π_0] = L[π_0] + back_coupling[π_0] = -∂_t S_base
```

where `∂_t S_base` is computed with all constraint momenta set to zero.

The back-coupling contribution is a product of Fourier multipliers. For
example, in coupled Proca:

- Source term: `-1 · gradient_x(π_1)` in A_0 constraint
- A_1 dynamical eq contains: `-1 · gradient_x(π_0)`
- Back-coupling multiplier: `(-1) · ik_x · (-1) · ik_x = -k_x²`
- Similarly for y: `(-1) · ik_y · (-1) · ik_y = -k_y²`
- Total: augmented diagonal = `(-k² - m²) + (-k²) = -2k² - m²`

#### Block FFT Solve for Coupled Constraints

When multiple constraint fields are coupled (e.g., `A_0 ↔ B_0` in coupled
Proca), the system is solved simultaneously via batched SVD:

```
M_aug(k) · [π_A0, π_B0, ...]^T = -[∂_t S_A0, ∂_t S_B0, ...]^T
```

where `M_aug(k)` has:
- **Diagonal entries:** self-operator multiplier + back-coupling
- **Off-diagonal entries:** from cross-constraint coupling (e.g., `gcoup`)

This uses the same SVD infrastructure and Tikhonov regularization as the
existing block FFT constraint solver (`_solve_coupled_constraints_fft`).

### Why This Is Better Than Finite Differences

The previous approach estimated `π_N` via backward FD:
`π_N ≈ (field(t) - field(t_prev)) / Δt`. This introduced a **time lag** that
turned the constraint-dynamical feedback loop into a delay differential
equation (DDE). Perturbative analysis of the DDE characteristic equation
showed a growth rate `σ ∝ k⁴ Δt / m²` — exponential blow-up of high-k modes
that no time step reduction could cure.

The analytical approach eliminates the time lag entirely:

- **No persistent state:** The method is a pure function of `(state, t)`,
  inherently compatible with adaptive ODE solvers (no anchor/caching needed).
- **Exact at every instant:** No `O(Δt)` error — the momentum is the exact
  solution of the differentiated constraint equation.
- **k-independent stability:** The augmented operator `L_aug` has strictly
  negative eigenvalues for all k, preventing instability at any wavelength.

### Implementation

Method `_compute_analytical_constraint_momenta()` runs after the constraint
solve, computing virtual momenta via `_solve_constraint_momenta_fft()`.

The method is theory-agnostic: it reads the equation structure from the JSON
spec and automatically identifies source terms, back-coupling chains, and
cross-constraint coupling without any theory-specific knowledge.

### Limitations

- **Periodic BCs only:** Requires FFT. For non-periodic grids, a future
  matrix-based analytical solve extension is needed.
- **Position-dependent self-term coefficients:** Same limitation as the FFT
  constraint solver — raises `ValueError` suggesting `method='matrix'`.
  Falls back to `π = 0` with a debug-level log message.
- **Time-dependent source coefficients:** `∂_t(c(t) · op(field))` requires
  `∂_t c`, which is not yet implemented. Raises `ValueError` if detected.
  (No current example has this case.)

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
