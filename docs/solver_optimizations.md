# Solver Performance Optimizations

This document describes the solver-level performance optimizations in TIDAL,
covering both Python micro-optimizations (Rounds 1-6) and algorithmic
improvements (higher-order FD stencils, Yoshida integrator).

## Summary of optimizations

| Phase | Feature | Key Benefit |
|-------|---------|-------------|
| Rounds 1-6 | Python hot-path micro-optimizations | Leapfrog 1.66x, scipy 1.29x |
| Phase 1 | `--fd-order 4\|6` (higher-order FD stencils) | 500x spatial error reduction per order at same N |
| Phase 2 | `--leapfrog-order 4` (Yoshida integrator) | 1.5x wall-clock at equal accuracy |

---

## Phase 1: Higher-order finite difference stencils

**CLI flag**: `--fd-order 2|4|6` (default: 4)

4th and 6th-order central-difference stencils for all spatial operators
(gradient, Laplacian, cross-derivative, biharmonic). Higher orders use
wider stencils and more ghost cells, but dramatically improve spatial
accuracy -- typically enabling halving N for the same accuracy.

### Stencil coefficients

**1st derivative (central):**

| Order | Points | Coefficients (x 1/dx) |
|-------|--------|----------------------|
| 2 | 3 | -1/2, 0, +1/2 |
| 4 | 5 | 1/12, -2/3, 0, 2/3, -1/12 |
| 6 | 7 | -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60 |

**2nd derivative (central):**

| Order | Points | Coefficients (x 1/dx^2) |
|-------|--------|------------------------|
| 2 | 3 | 1, -2, 1 |
| 4 | 5 | -1/12, 4/3, -5/2, 4/3, -1/12 |
| 6 | 7 | 1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90 |

Ref: Fornberg (1988), Mathematics of Computation 51(184), pp. 699-706.

### Spatial accuracy (N=64, gradient of sin(x))

| FD Order | Max Error | Improvement vs order 2 |
|----------|-----------|------------------------|
| 2 | 1.6e-3 | baseline |
| 4 | 3.1e-6 | 500x |
| 6 | 6.4e-9 | 500x |

### Architecture

- Module-level state: `_fd_order`, `_n_ghost`, `set_fd_order()`, `get_fd_order()`
- Variable-width ghost cells: 1/2/3 per side for orders 2/4/6
- All ghost-cell slice tuples pre-cached in `_PadEntry.__init__` for zero
  per-call allocation overhead
- Cross-derivative via successive `gradient()` calls -- automatically inherits FD order
- Biharmonic = `laplacian(laplacian(data))` -- inherits order automatically
- Recursive ghost fill for non-periodic BCs (layer-by-layer from interior outward)

### Energy measurement consistency

The energy measurement module (`_energy.py`) delegates all spatial operators
to `operators.py` instead of maintaining independent duplicate stencils. This
ensures energy measurement tracks the discrete (shadow) Hamiltonian exactly
when FD order changes -- critical for conservation diagnostics.

Ref: Hairer, Lubich & Wanner (2006), "Geometric Numerical Integration",
Ch. VI, Springer -- stencil-measurement consistency.

### Constraint solver FFT multipliers

The FFT constraint solver (`constraint_solve.py`) uses FD-order-aware modified
wavenumber Fourier symbols:

- Order 2: `k_eff^2 = -(2/h^2)(1 - cos(kh))`
- Order 4: `k_eff^2 = (1/h^2)(-cos(2kh)/6 + 8cos(kh)/3 - 5/2)`
- Order 6: `k_eff^2 = (1/h^2)(cos(3kh)/45 - 3cos(2kh)/10 + 3cos(kh) - 49/18)`

This ensures the constraint IC solution is spectrally consistent with the
FD stencils used during time evolution.

### Sparsity pattern

`sparsity.py` returns wider offset sets based on `get_n_ghost()`:
- Order 2: [-1, 0, +1]
- Order 4: [-2, -1, 0, +1, +2]
- Order 6: [-3, -2, -1, 0, +1, +2, +3]

### Persistence

FD order is saved in `metadata.json` via `_writer.py` and restored by `_io.py`
on load, so that measurement tools use the correct stencils.

### CLI default

The CLI default is `--fd-order 4` (5-point Fornberg stencil), which provides
O(dx^4) convergence at ~10% per-step overhead vs order 2. The module-level
default stays at 2 for backward compatibility with library/test code.

### Files modified

`operators.py`, `sparsity.py`, `_energy.py`, `_writer.py`, `_io.py`,
`_simulate.py`, `cli/__init__.py`, `constraint_solve.py`,
`tests/test_fd_order.py` (18 tests), `tests/test_measurement.py` (3 tests updated)

### References

- Fornberg, B. (1988). "Generation of Finite Difference Formulas on Arbitrarily
  Spaced Grids". Mathematics of Computation 51(184), pp. 699-706.
- LeVeque, R.J. (2007). "Finite Difference Methods for Ordinary and Partial
  Differential Equations". SIAM. Ch. 1-3.
- Hairer, E., Lubich, C. & Wanner, G. (2006). "Geometric Numerical Integration".
  Springer. Ch. VI.

---

## Phase 2: Yoshida 4th-order symplectic leapfrog

**CLI flag**: `--leapfrog-order 2|4` (default: 2)

Yoshida triple-composition of Stormer-Verlet sub-steps:

    S4(dt) = S2(w1*dt) o S2(w2*dt) o S2(w3*dt)

where w1 = w3 = 1/(2 - 2^(1/3)) ~ 1.3512, w2 = -2^(1/3)/(2 - 2^(1/3)) ~ -1.7024.

The negative middle weight w2 causes a backward sub-step that cancels the
leading O(dt^2) error term, achieving O(dt^4) accuracy.

### Key insight: speedup comes from accuracy, not per-step cost

Each individual Yoshida step is approximately 3x more expensive than a
Stormer-Verlet step (6 force evaluations vs 2). However, the O(dt^4)
accuracy means far fewer total steps are needed to reach a given error
target. The speedup is realized when comparing wall-clock time at
**equal accuracy**, not equal step count:

| Integrator | dt | Steps | Force evals | Wall-clock | Error |
|------------|--------|-------|-------------|------------|---------|
| Verlet | 0.020 | 1000 | 2001 | 0.071s | 7.5e-5 |
| Verlet | 0.010 | 2000 | 4001 | 0.118s | 1.9e-5 |
| Verlet | 0.005 | 4000 | 8001 | 0.224s | 4.5e-6 |
| Yoshida | 0.050 | 400 | 2400 | 0.117s | 7.8e-6 |
| Yoshida | 0.040 | 500 | 3000 | 0.147s | 3.3e-6 |
| Yoshida | 0.020 | 1000 | 6000 | 0.280s | 3.8e-7 |

**Equal-accuracy comparison (~4e-6 error):**
- Verlet dt=0.005: 4000 steps, 8001 force evals, 0.21s
- Yoshida dt=0.040: 500 steps, 3000 force evals, 0.14s
- **Speedup: 1.5x wall-clock, 2.7x fewer force evaluations**

Benchmark: coupled scalars, N=256, fd_order=4, t_end=20.

### Convergence order verification

Temporal convergence measured by comparing to a fine-dt reference solution
on the same spatial grid (eliminating spatial error):

- **Verlet slopes: ~2.0** (O(dt^2)) -- confirmed
- **Yoshida slopes: ~4.0** (O(dt^4)) -- confirmed

### Energy conservation

Both integrators preserve a shadow Hamiltonian to machine precision (zero
secular drift). The shadow Hamiltonian error scales with the integrator order:
- Verlet: |dE/E| < 1e-4 at dt=0.01 over t=10
- Yoshida: |dE/E| < 1e-6 at dt=0.01 over t=10 (100x tighter)

### Architecture

- Module constants: `_W1`, `_W2`, `YOSHIDA_WEIGHTS` (exported)
- `solve_leapfrog()` accepts `order=2|4` keyword argument (default 2)
- Order 2 path unchanged (KDK with force caching, 1 force eval/step)
- Order 4: 3 KDK sub-steps per outer step, 6 force evals total
  (no inter-sub-step force caching possible because position changes)
- Sub-step time tracked via `t_sub` variable for correct handling of
  time-dependent coefficients (background fields, time-varying sources).
  This makes the implementation general -- not limited to time-independent
  systems like Klein-Gordon.

### CFL auto-adjustment

The effective CFL limit for Yoshida is reduced by max(|wi|) ~ 1.70 because
the middle sub-step has |w2| > 1. The CLI (`_simulate.py`) auto-adjusts dt
by this factor when `--leapfrog-order 4` is specified, before snapshot
configuration to ensure correct writer pre-allocation.

### Files modified

`leapfrog.py`, `cli/__init__.py`, `cli/_simulate.py`,
`tests/test_solver_leapfrog.py` (17 tests total, 9 new Yoshida tests)

### References

- Yoshida, H. (1990). "Construction of higher order symplectic integrators".
  Physics Letters A, 150(5-7), pp. 262-268.
- Hairer, E., Lubich, C. & Wanner, G. (2006). "Geometric Numerical Integration".
  Springer. Ch. VI: Symplectic Integration of Hamiltonian Systems.

---

## Python micro-optimizations (Rounds 1-6)

These optimizations target the Python hot path in the solver loop. IDA and
CVODE spend ~78% of wall time in C code (SUNDIALS), so Python-side gains
are limited for those backends. Leapfrog (pure Python) benefits the most.

### Round 1: Operator caching in CoefficientEvaluator

Precompute L0-resolved coefficients at spec load time.

### Round 2: Grid-shape buffer reuse

Pre-allocated output buffers for spatial operators in `operators.py`.

### Round 3: Slot groups, FieldSet reuse, RHS buffers, regex

- `StateLayout` slot groups: 4 `@cached_property` methods returning tuples
  of `(int, slice, str)` -- eliminates per-timestep branching in all 4 solvers
- `FieldSet.rebind()`: zero-allocation flat array swap, replaces
  `FieldSet.from_flat()` per timestep
- RHS term buffer: `np.multiply(coeff, operated, out=temp)` eliminates
  N-1 temporary arrays per `evaluate()`
- Pre-compiled regex patterns in `_eval_utils.py`

### Round 4: Hot-path operator optimizations

- `np.roll` replaced by ghost-cell padding via `_pad_axis()` + slice views
  (1.7-2.4x per operator)
- `grid.dx` as `@cached_property` (13.4x faster access)
- Pre-resolved operator dispatch in `RHSEvaluator.__init__`
- Fused `laplacian()` normalizes BCs once for all axes

IDA speedup: 1.13-1.40x.

### Round 5: Per-call allocation elimination

- `_str_to_axis_bc` singleton cache (6.1x faster)
- `_PadEntry` + `_PadBufferCache`: pre-allocated padded buffers + cached
  slice tuples per (shape, axis)
- In-place stencil arithmetic eliminates 2 temporaries per Laplacian

Scipy end-to-end: 1.29x (22.4% faster).

### Round 6: Inner solver hot-path

- KDK force caching in leapfrog: halves force evaluations (2N -> N+1)
- Zero-copy drift via `StateLayout.drift_slot_pairs`
- IDA `begin_timestep` dedup
- FieldSet pre-computed offsets
- `CoefficientEvaluator.resolve()` fast path

Leapfrog speedup: 1.66x.

---

## Planned: Phase 3 -- FFT spectral operators

**CLI flag**: `--spectral` (not yet implemented)

For all-periodic BCs, replace FD operators with FFT-based spectral operators
for exponential convergence. Expected N=64-128 instead of N=512-1024 for
the same (or better) accuracy on smooth problems.

Ref: Burns et al. (2020). "Dedalus: A flexible framework for numerical
simulations with spectral methods". Physical Review Research 2:023068.
