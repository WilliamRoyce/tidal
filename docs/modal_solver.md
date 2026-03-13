# Fourier Modal Solver

## Overview

The Fourier modal solver (`tidal/solver/modal.py`) provides exact, machine-precision solutions for linear PDE systems with periodic boundary conditions and time-independent coefficients. It transforms the spatial grid to Fourier space, builds a per-mode evolution matrix, and eigendecomposes to obtain the exact solution y(t) = exp(A·t)·y₀ — eliminating spatial discretization error entirely.

**Key advantages over PDE solvers (CVODE/IDA):**
- Machine-precision accuracy (~1e-14 error, limited only by eigendecomposition)
- No CFL condition or adaptive timestepping needed
- O(N²) per output time after O(N³) precomputation
- Exact for any time point — no time-stepping error accumulation

**Limitations:**
- Requires periodic boundary conditions (FFT periodicity)
- Requires time-independent coefficients (evolution matrix must be constant)
- Requires flat (Minkowski) metric — curved metrics couple modes via position-dependent coefficients that decay too slowly in k-space
- Currently pure NumPy/SciPy — no GPU acceleration

## How It Works

### 1. First-Order Reformulation

TIDAL stores state in Euler-Lagrange velocity form (q, v) where v = dq/dt. A second-order PDE ∂²t(φ) = L[φ] becomes:

```
dφ/dt = v
dv/dt = L[φ]
```

The `StateLayout` organizes field+velocity slot pairs. The modal solver reads this directly.

### 2. FFT to Fourier Space

Each field and velocity slot is transformed: u(x) → û(k) via `numpy.fft.rfftn`. Spatial operators become multiplications:

| Operator | Fourier multiplier |
|----------|-------------------|
| identity | 1 |
| laplacian | -k² |
| laplacian_x | -kₓ² |
| gradient_x | i·kₓ |
| cross_derivative_xy | -kₓ·kᵧ |
| biharmonic | k⁴ |

**Fourier convention:** k = 2π·`numpy.fft.rfftfreq(N, d=dx)` — consistent with `operators.py` spectral mode. NOT the modified-wavenumber convention from `constraint_solve.py`.

### 3. Build Evolution Matrix

**Constant coefficients (all terms position-independent):** Each k-mode decouples. For a 2-field second-order system, each mode yields a 4×4 matrix:

```
d/dt [φ̂_k, v̂φ_k, χ̂_k, v̂χ_k]ᵀ = M(k) · [φ̂_k, v̂φ_k, χ̂_k, v̂χ_k]ᵀ
```

Result: N_modes independent small matrices. Total: O(N_modes × N_fields³).

**Position-dependent coefficients (convolution coupling):** When a coefficient c(x) multiplies a field, the product in k-space becomes a convolution:

```
FFT[c(x) · u(x)] = ĉ ∗ û  (k-space convolution)
```

This couples different k-modes, creating a full (N_modes × N_slots)² matrix. Built via probe vectors (unit impulse per mode): reconstruct to physical space, multiply by coefficient, FFT back. For Gaussian c(x), the kernel ĉ(q) decays exponentially → banded matrix.

### 4. Eigendecomposition and Exact Evolution

The evolution matrix A is constant → exact solution via eigendecomposition:

```python
eigenvalues, V = np.linalg.eig(A)
V_inv = np.linalg.inv(V)
y0_eigen = V_inv @ y0_k

for t in t_eval:
    y_k(t) = V @ (exp(eigenvalues * (t - t0)) * y0_eigen)
    y(t) = IFFT(y_k(t))
```

This is **exact** — no time-stepping error, no CFL condition, no tolerance parameter needed.

### 5. Tiered Matrix Algorithms

| Matrix size | Algorithm | Reference |
|-------------|-----------|-----------|
| Per-mode blocks ≤ 12×12 | Batch dense `np.linalg.eig` | Golub & Van Loan (1996), §4.8 |
| Full matrix ≤ 2000×2000 | Dense eigendecomposition | `np.linalg.eig` |
| Full matrix > 2000×2000 | Schur decomposition | Moler & Van Loan (2003), SIAM Review 45(1):3-49 |

Schur decomposition (T, Z = schur(A)) provides better numerical stability for large matrices: eigenvalues from diag(T), evolution via Z·diag(exp(λ·t))·Zᴴ (unitary Z gives exact inverse).

## Auto-Selection

The modal solver is auto-selected (step 2 in `_resolve_scheme()`) when ALL conditions are met:

1. **Flat metric** — `spec.canonical.volume_element is None` (checked first; curved metrics like spherical/cylindrical have non-None volume elements)
2. **No constraint equations** — no time_order=0 equations
3. **All BCs periodic** — `grid.periodic` all True
4. **All operators supported** — all RHS operators in `_EXACT_MULTIPLIERS`
5. **No time-dependent coefficients** — `term.time_dependent` all False

Position-dependent coefficients are OK (handled via convolution). The full auto-selection hierarchy:

```
1. Constraints (time_order=0)         → IDA
2. Modal eligible (flat+periodic+...) → modal  ← NEW
3. First-order (time_order=1)         → IDA
4. Dissipation (first_derivative_t)   → IDA
5. Missing canonical (wave eq)        → warning + CVODE
6. Pure wave, Hamiltonian             → CVODE
```

Users can explicitly request `--scheme modal` (validated for eligibility) or override with `--scheme cvode`.

## Usage

```bash
# Auto-selected for eligible systems:
tidal simulate examples/data/coupled_scalars.json --periodic --t-end 10

# Explicit selection:
tidal simulate examples/data/coupled_scalars.json --scheme modal --periodic --t-end 10

# In sweeps:
tidal sweep examples/data/coupled_scalars.json --sweep "gCpl=0.1:1.0:20" \
  --measure conversion --scheme modal
```

## First-Order Systems (Diffusion)

The modal solver handles first-order equations naturally:

```
∂t(u) = D·∂²x(u)  →  per-mode: dû_k/dt = -D·k²·û_k
```

Solution: û_k(t) = û_k(0)·exp(-D·k²·t) — exponential decay per mode, diagonal matrix.

## Error Model

The **only** error source is eigendecomposition precision (~1e-14 for well-conditioned matrices). For parameter sweeps, numerical uncertainty is effectively zero — the dominant error becomes the physics model itself.

Comparison with other solvers on single-field KG (k=1, m²=1, N=64, t=1):

| Solver | Error vs analytical | Note |
|--------|---------------------|------|
| Modal | 5×10⁻¹⁶ | Machine precision |
| CVODE (rtol=1e-10) | 2.8×10⁻⁴ | Truncation error accumulates |

### Performance Benchmarks

Coupled scalars (2 fields, 31 snapshots, t_end=3, periodic BCs):

| Grid | Modal | CVODE (rtol=1e-10) | Speedup | Max diff |
|------|-------|--------------------|---------|----------|
| N=64 | 0.003 s | 0.16 s | **57×** | 4.0×10⁻⁴ |
| N=128 | 0.003 s | 0.40 s | **157×** | 5.2×10⁻⁴ |
| N=256 | 0.003 s | 4.27 s | **1451×** | 6.3×10⁻⁴ |

Key observations:
- **Modal time is O(1) in grid size** — eigendecomposition cost is per-mode (4×4 matrices), independent of N. The IFFT reconstruction is O(N log N) but negligible.
- **CVODE scales as O(N²+)** — spatial operator evaluation is O(N) per timestep, and more grid points require more timesteps to resolve the same physics.
- **Max diff is CVODE error, not modal error** — modal solutions are exact to machine precision (~10⁻¹⁶); the 10⁻⁴ difference is entirely CVODE truncation error.
- **Speedup grows with N** — at N=256 the modal solver is >1000× faster, making it transformative for parameter sweeps where hundreds of simulations are needed.

## Implementation Details

### Files

| File | Purpose |
|------|---------|
| `tidal/solver/modal.py` | Core solver (~760 lines) |
| `tidal/cli/_simulate.py` | Auto-selection + dispatch |
| `tidal/cli/__init__.py` | `"modal"` in `--scheme` choices |
| `tests/test_solver_modal.py` | 25 tests (11 eligibility + 14 correctness) |

### Key Functions

| Function | Purpose |
|----------|---------|
| `can_use_modal(spec, grid, bc)` | Eligibility check (shared with auto-selection) |
| `solve_modal(spec, grid, y0, ...)` | Main entry point |
| `_fft_slots / _ifft_slots` | State ↔ Fourier transforms |
| `_build_per_mode_matrices` | Per-mode evolution matrices (constant coefficients) |
| `_build_convolution_matrix` | Full convolution matrix (position-dependent) |
| `_evolve_per_mode` | Batch eigendecomposition solver |
| `_evolve_full_matrix` | Full matrix solver (dense or Schur) |

### Reused Infrastructure

| Component | Location | How Used |
|-----------|----------|----------|
| `_get_wavenumbers()` | `operators.py` | k-grid construction |
| `is_periodic_bc()` | `operators.py` | Eligibility gate |
| `StateLayout` | `state.py` | Field↔state mapping |
| `CoefficientEvaluator` | `coefficients.py` | Coefficient resolution |
| `SimulationProgress` | `progress.py` | tqdm progress bar |

**NOT reused:** `constraint_solve.py:_MULTIPLIERS` — uses modified wavenumbers matching FD stencils. Modal solver defines its own exact multiplier registry.

## Constraints and Limitations

1. **Constraint equations** (time_order=0): Incompatible — algebraic constraints are not ODEs. Systems with constraints auto-select IDA.

2. **Curved metrics** (spherical, cylindrical, polar): Not applicable — non-periodic domains and position-dependent operators with slowly-decaying Fourier transforms create dense all-to-all coupling.

3. **Time-dependent coefficients**: Would require ODE integration in modal space (`scipy.integrate.solve_ivp`). Currently rejected by eligibility check (no known TIDAL examples need this).

4. **Spectral operators**: Redundant when using modal solver (both operate in k-space). The `--spectral` flag is silently disabled when modal is selected.

## References

- Moler, C. & Van Loan, C. (2003). "Nineteen dubious ways to compute the exponential of a matrix, twenty-five years later." *SIAM Review*, 45(1):3-49. DOI: [10.1137/S00361445024180](https://doi.org/10.1137/S00361445024180)
- Golub, G.H. & Van Loan, C.F. (1996). *Matrix Computations*, 3rd ed. Johns Hopkins University Press. ISBN: 978-0801854149
- Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*, Springer, §4. DOI: [10.1007/3-540-30666-8](https://doi.org/10.1007/3-540-30666-8)
- Burns, K.J. et al. (2020). "Dedalus: A flexible framework for numerical simulations with spectral methods." *Physical Review Research*, 2:023068. DOI: [10.1103/PhysRevResearch.2.023068](https://doi.org/10.1103/PhysRevResearch.2.023068)
- Raffelt, G. & Stodolsky, L. (1988). "Mixing of photons with low-mass particles in magnetic fields." *Physical Review D*, 37:1237. DOI: [10.1103/PhysRevD.37.1237](https://doi.org/10.1103/PhysRevD.37.1237)
- Dormand, J.R. & Prince, P.J. (1980). "A family of embedded Runge-Kutta formulae." *Journal of Computational and Applied Mathematics*, 6. DOI: [10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X(80)90013-3)
- Hindmarsh, A.C. et al. (2005). "SUNDIALS: Suite of nonlinear and differential/algebraic equation solvers." *ACM TOMS*, 31(3):363-396. DOI: [10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020)
