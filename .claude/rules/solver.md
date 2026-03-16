---
paths:
  - "tidal/solver/**"
  - "tidal/measurement/**"
---

# Solver Architecture Rules

## Solver Selection (auto, by priority)
1. **Modal** — flat metric + periodic BCs + time-independent + supported operators → eigendecomposition
2. **Leapfrog** — no constraints + no first_derivative_t + no dissipation → symplectic
3. **CVODE** — no constraints → adaptive BDF
4. **IDA** — has constraints (time_order=0 fields) → implicit DAE
5. **scipy** — fallback

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
