---
paths:
  - "tidal/solver/**"
  - "tidal/measurement/**"
---
> **SCOPE (2026-09-04): this file describes the LEGACY solver, `tidal/solver/`.**
>
> It is accurate for that code, which still exists and still runs — but the auto-selection
> ladder below is **not** the architecture going forward. WS3 (#492) replaces it with a
> different design: two engines over one shared core (an oscillation-resolving mode-equation
> solver for O2, an eikonal amplitude engine with coherence-patch averaging for O3), chosen
> by bake-off rather than by a fixed priority list. See `docs/cosmology/solver_design.md`.
>
> **Live wrong-answer hazard (#517).** Modal eligibility is checked **only in
> `can_use_modal`**, not inside the solver. A *library* call that constructs the modal solver
> directly on a time-dependent spec does not hit that check and will silently freeze the
> coefficients at `t = 0`. That has been harmless while FRW was out of scope; it stops being
> harmless the moment time-dependent specs become normal. Do not assume the refusals at
> `modal.py:371-372` and `:409-413` protect a non-CLI caller.

# Solver Architecture Rules

## Solver Selection (auto, by priority)
1. **Modal** — flat metric + periodic BCs + time-independent + 15 supported operators → Padé matrix-exp (path D, `scipy.linalg.expm`). Position-dependent kinetics supported via real-space M⁻¹(x) folded into the convolution paths, which run in the FULL fft basis (GH #427/#445); the per-mode genEig engine, the stability probe, Pass-1 Duhamel, and modal-jax still refuse them. Eigendecomposition retired v0.31+. Perturbative eligibility is judged on the canonicalized base spec.
2. **IDA** — constraints present (time_order=0 fields) → implicit SUNDIALS DAE
3. **IDA** — first-order equations (time_order=1, diffusion) → implicit BDF
4. **IDA** — dissipation (`first_derivative_t` in RHS) → implicit BDF
5. **CVODE** — pure wave equations → adaptive BDF (default fallback)

**Leapfrog and scipy are manual-only** (`--scheme leapfrog` / `--scheme scipy`); never auto-selected.

## Key Classes
- `FieldSet` — field metadata, slot indices, rebind(). Pre-computed offsets for hot paths.
- `CoefficientEvaluator` — 4-level cache (raw → param-resolved → position-dep → field-combined). `resolve()` evaluates with actual params.
- `RHSEvaluator` — evaluates dq/dt, dv/dt from equation terms. Pre-resolved operators/BCs.
- `StateLayout` — slot groups (dynamic, velocity, constraint). `from_spec()` factory.
- `GridInfo` — grid shape, bounds, dx (cached), ghost cells. `min(shape) >= fd_order + 1`.

## Operator System
- Identity, laplacian, gradient_{x,y,z}, cross_derivative, first_derivative_t, biharmonic
- Mixed: `mixed_T_S1_S2_...` via `_evaluate_mixed_factor` in `_energy.py`
- FD stencils: `--fd-order 2|4|6` (Fornberg 1988). Spectral: `--spectral` (auto for periodic).

## Energy Measurement
- `hamiltonian_terms` in JSON + optional `volume_element`
- Raw coefficients assume params=1.0 — always use `CoefficientEvaluator.resolve()`
- Position-dependent mass breaks energy measurement — use conversion measurement instead

## Constraint IC
- `ensure_consistent_ic()` in `constraint_solve.py` — runs after `_build_initial_y0()`
- Solves algebraic constraints for time_order=0 fields
