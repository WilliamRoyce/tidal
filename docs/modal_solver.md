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

### 4. Time Evolution Algorithms

Two paths depending on coefficient structure:

**Per-mode eigendecomposition (constant coefficients):**

When all coefficients are constant, the evolution matrix is block-diagonal in mode space. Each mode's small block (≤ 12×12) is eigendecomposed independently:

```python
eigenvalues, V = np.linalg.eig(A_k)  # per-mode block
y_k(t) = V @ (exp(eigenvalues * t) * V_inv @ y0_k)
```

This has no time-stepping error and no CFL condition. O(N³_fields) per mode, O(N_modes × N³_fields) total.

### Nyquist Mode Zeroing (Critical for Energy Conservation)

The per-mode eigendecomposition achieves **machine-precision energy conservation** (|dE/E| ≈ 2e-14) after zeroing the Nyquist mode in the initial conditions. Without this, conservation degrades to ~1.5e-5.

**Root cause**: The rfft Nyquist bin (k = N/2, the last mode in even-N grids) must have a purely real Fourier coefficient for real-valued fields. However, the per-mode evolution matrix A_k has complex entries from the gradient coupling (ik term). When exp(A_k·t) evolves the Nyquist mode, it creates imaginary components in the Fourier coefficient. The `irfft` reconstruction silently drops these imaginary parts (by construction — the Nyquist bin of a real signal IS real). This truncation discards energy from the Nyquist mode at every snapshot.

**Investigation**: Per-mode analysis showed that ALL modes except the Nyquist are exactly conserved to machine precision (dH_k = 0 for k < N/2). The Nyquist mode's Hamiltonian grew by 500× (from 0.029 to 16.8) in the complex eigendecomp, but the stored (irfft-reconstructed) field showed shrinking amplitude — the discarded imaginary part contained the "missing" energy.

**Fix**: Zero the Nyquist mode's IC before eigendecomposition. This is standard practice in pseudospectral methods — the Nyquist frequency aliases with its negative-frequency conjugate and cannot faithfully represent directional (travelling wave) content. It is the highest resolvable frequency and carries negligible physical content for well-resolved simulations.

**Results** (coupled_scalars, Gertsenshtein gradient coupling):

| Grid | Before fix | After fix |
|------|-----------|-----------|
| N=128 | |dE/E| = 2.8e-5 | |dE/E| = 2.2e-14 |
| N=256 | |dE/E| = 1.8e-5 | |dE/E| = 2.2e-14 |
| N=512 | |dE/E| = 1.7e-5 | |dE/E| = 2.2e-14 |
| N=1024 | |dE/E| = 1.6e-5 | |dE/E| = 2.2e-14 |

Full Gertsenshtein (6 fields): |dE/E| = 1.3e-14.

**References**:
- Boyd, J.P. (2001). *Chebyshev and Fourier Spectral Methods*, 2nd ed. Dover. §11.5 (aliasing and the Nyquist mode).
- Canuto, C. et al. (2006). *Spectral Methods: Fundamentals in Single Domains*. Springer. §3.2 (dealiasing via mode truncation).
- Burns et al. (2020). Dedalus framework zeroes the top 1/3 of modes ("2/3 dealiasing rule") to prevent aliasing in nonlinear terms. Our Nyquist zeroing is the minimal version for linear systems.

**Krylov matrix exponential (position-dependent coefficients):**

Position-dependent convolution produces a full (N_slots × N_modes)² matrix that is generally **non-normal**: gradient operators (ik) combined with real convolution kernels create eigenvalues with significant positive real parts despite conservative physics. Eigendecomposition overflows because individual exp(λ·t) diverge even though exp(A·t)·y₀ is bounded (pseudospectral phenomenon; Trefethen & Embree 2005, Ch. 14).

The fix: `scipy.sparse.linalg.expm_multiply` computes exp(A·t)·y₀ directly via scaling + truncated Taylor series in matrix-vector products, which is backward-stable for non-normal matrices:

```python
from scipy.sparse.linalg import expm_multiply
y_all = expm_multiply(A_full, y0_flat, start=0, stop=t_end, num=n_snapshots)
```

Ref: Al-Mohy & Higham (2011), "Computing the Action of the Matrix Exponential", SIAM J. Sci. Comput. 33(2):488-511.

### 5. Algorithm Selection

| Coefficient type | Matrix structure | Algorithm | Time (Gertsenshtein, N=512) |
|------------------|-----------------|-----------|---------------------------|
| Constant | Per-mode blocks (12×12) | Eigendecomposition | ~1.5s |
| Position-dependent | Full convolution (3000×3000) | `expm_multiply` | ~21s |

The routing is automatic — `_has_position_dependent_terms()` determines which path is used. The per-mode path is ~14x faster because it avoids building the full matrix.

## Auto-Selection

The modal solver is auto-selected (step 2 in `_resolve_scheme()`) when ALL conditions are met:

1. **Flat metric** — `spec.canonical.volume_element is None` (checked first; curved metrics like spherical/cylindrical have non-None volume elements)
2. **No constraint equations** — no time_order=0 equations
3. **All BCs periodic** — `grid.periodic` all True
4. **All operators supported** — all RHS operators in `_EXACT_MULTIPLIERS`
5. **No time-dependent coefficients** — `term.time_dependent` all False

Position-dependent coefficients are OK (handled via convolution). The full auto-selection hierarchy:

```
1. Modal eligible (flat+periodic+time-independent+supported operators)
   → modal (handles constraints via Schur complement if present)
2. Constraints (time_order=0) not modal-eligible → IDA
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
- **Modal time is O(1) in simulation time** — eigendecomposition is done once; evaluating at any t is just `exp(λ·dt)`. CVODE must take more timesteps for longer runs. See scaling table below.

### Gertsenshtein Effect (6 fields, gradient coupling)

B₀ sweep (N=256, L=100, t_end=50, k=2.01, κ=1, 10 points), P_final compared to analytical P = sin²(κB₀D/2) × k²/(k²+κ²B₀²):

| B₀ | P(modal) | P(CVODE) | P(analytical) | Error |
|----|----------|----------|---------------|-------|
| 0.010 | 0.0612 | 0.0612 | 0.0612 | 2×10⁻⁶ |
| 0.031 | 0.4924 | 0.4924 | 0.4923 | 1×10⁻⁴ |
| 0.052 | 0.9313 | 0.9312 | 0.9307 | 6×10⁻⁴ |
| 0.073 | 0.9325 | 0.9325 | 0.9314 | 1×10⁻³ |
| 0.094 | 0.4949 | 0.4948 | 0.4940 | 9×10⁻⁴ |
| 0.116 | 0.0624 | 0.0624 | 0.0623 | 5×10⁻⁵ |
| 0.137 | 0.0740 | 0.0740 | 0.0734 | 6×10⁻⁴ |
| 0.158 | 0.5180 | 0.5180 | 0.5143 | 4×10⁻³ |
| 0.179 | 0.9435 | 0.9435 | 0.9360 | 7×10⁻³ |
| 0.200 | 0.9181 | 0.9181 | 0.9105 | 8×10⁻³ |

- **Both solvers agree to ~10⁻⁵** with each other — identical physics, different numerics
- **RMS error vs analytical: 0.0036** (0.36%) — dominated by the effective-mass correction approximation, not solver error
- **Modal: 0.49s/point, CVODE: 1.66s/point → 3.4× speedup** at N=256
- Speedup is lower than coupled_scalars because the 6-field Gertsenshtein system has 12×12 per-mode matrices (vs 4×4 for 2 fields) and 6× more IFFT reconstructions per snapshot

### Scaling with Simulation Time (t_end)

Gertsenshtein (6 fields, N=256, 101 snapshots, B₀=0.1):

| t_end | Modal | CVODE | Speedup |
|-------|-------|-------|---------|
| 10 | 1.4 s | 1.7 s | **1.2×** |
| 50 | 1.5 s | 3.1 s | **2.0×** |
| 100 | 1.6 s | 4.3 s | **2.8×** |
| 200 | 2.1 s | 7.2 s | **3.4×** |
| 500 | 1.6 s | 12.2 s | **7.7×** |

- **Modal cost is constant in t_end** (~1.5s) — eigendecomposition is done once, and evaluating exp(λ·Δt) at each output time is O(1) regardless of Δt. The ~1.5s is Python startup + FFT + eigendecomposition overhead.
- **CVODE cost scales linearly with t_end** — each additional unit of simulation time requires proportionally more adaptive timesteps.
- **Speedup grows linearly with t_end** — for a 40-point parameter sweep at t_end=500, modal saves ~7 minutes vs CVODE. For t_end=5000 (astrophysical timescales), the advantage would be ~50×.

## Implementation Details

### Files

| File | Purpose |
|------|---------|
| `tidal/solver/modal.py` | Core solver (~760 lines) |
| `tidal/cli/_simulate.py` | Auto-selection + dispatch |
| `tidal/cli/__init__.py` | `"modal"` in `--scheme` choices |
| `tests/test_solver_modal.py` | 35 tests (12 eligibility + 14 correctness + 4 stability + 5 constraint) |

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

## Block-Aware Eigendecomposition

Multi-field systems often have block-diagonal per-mode matrices (e.g. Gertsenshtein: h₅↔a₁ and h₇↔a₂ as independent 4×4 blocks). The modal solver detects these blocks via union-find on the coupling graph and eigendecomposes each independently. This:

1. **Prevents degenerate-eigenvalue mixing** — `np.linalg.eig` on a full matrix with repeated eigenvalues across independent blocks can mix eigenvectors, projecting nonzero-IC components onto zero-IC blocks
2. **Enables zero-IC block skipping** — blocks where all initial conditions are zero produce exact zeros without eigendecomposition
3. **Reduces computation** — many small eigendecompositions instead of one large one

For systems with gradient coupling (e.g. Gertsenshtein h↔a), low-k modes may have genuinely positive real eigenvalues (physical parametric conversion). A warning is issued when `max(Re(λ)) × Δt > 30` (approaching overflow).

## Constraint Elimination (Fourier Schur Complement)

Systems with constraint equations (time_order=0, e.g., Proca A₀) are now handled via algebraic elimination in Fourier space. The constraint self-operator L (e.g., m²-∇²) has exact Fourier multipliers, so L⁻¹ is trivially 1/(m²+k²) per mode.

**How it works**: For a mixed system with dynamical (d) and constraint (c) fields:

```
d/dt[d] = A_dd·d + A_dc·c     (dynamical equations)
    0   = S_cd·d + S_cc·c      (constraint defines c algebraically)
```

Substituting c = -S_cc⁻¹·S_cd·d gives the reduced system:

```
d/dt[d] = (A_dd - A_dc·S_cc⁻¹·S_cd)·d
```

The constraint velocity coupling (when dynamical equations reference v_A₀ = ∂ₜA₀) creates an implicit equation resolved by matrix inversion of (I - A_dc_vel·S_cc⁻¹·S_cd) — a small per-mode operation.

At each output time, constraint fields are reconstructed: c(k) = -S_cc⁻¹·S_cd·d(k).

**Performance**: On coupled_proca_3d (16×16, t=5): modal 2.0s vs IDA 19.0s — **9.5× speedup**. IDA fails entirely at 32×32 (convergence issues), while modal runs at any grid size.

**References**: Hairer & Wanner (1996), *Solving ODEs II*, Ch. VII; Ascher & Petzold (1998), *Computer Methods for ODEs/DAEs*, §10.2.

## Constraints and Limitations

1. **Constraint equations** (time_order=0): Handled via Fourier Schur complement (see above) when all constraint operators have exact Fourier multipliers. Otherwise falls back to IDA.

2. **Curved metrics** (spherical, cylindrical, polar): Not applicable — non-periodic domains and position-dependent operators with slowly-decaying Fourier transforms create dense all-to-all coupling.

3. **Time-dependent coefficients**: Would require ODE integration in modal space (`scipy.integrate.solve_ivp`). Currently rejected by eligibility check (no known TIDAL examples need this).

4. **Spectral operators and energy measurement**: The modal solver works in k-space, but the energy measurement evaluates the Hamiltonian in physical space. The `--spectral` flag is preserved (not disabled) for modal so that the energy measurement uses FFT-based gradient operators matching the modal solver's conserved Hamiltonian. Without this, FD gradients produce conservation errors that increase with N (the FD Hamiltonian differs from the Fourier Hamiltonian). For gradient×gradient terms in the Hamiltonian, the measurement uses direct gradient product ⟨∂f, ∂g⟩ rather than IBP ⟨-f, ∂²g⟩ when spectral is active, because rfft Nyquist-mode handling creates O(1/N) discrepancy between the two.

5. **Per-mode eigendecomposition conditioning**: For coupled multi-field systems where fields have similar dispersion (e.g., massless graviton-photon), the 4×4 per-mode block A = [[0, I], [K, 0]] has near-degenerate eigenvalue pairs, causing cond(V) up to 10^291. This limits energy conservation to ~1.5e-5 (eigenvector reconstruction error). The solution (future work) is to eigendecompose the 2×2 K matrix and use Hamiltonian cos/sin structure, reducing conditioning to O(1).

## References

- Moler, C. & Van Loan, C. (2003). "Nineteen dubious ways to compute the exponential of a matrix, twenty-five years later." *SIAM Review*, 45(1):3-49. DOI: [10.1137/S00361445024180](https://doi.org/10.1137/S00361445024180)
- Golub, G.H. & Van Loan, C.F. (1996). *Matrix Computations*, 3rd ed. Johns Hopkins University Press. ISBN: 978-0801854149
- Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*, Springer, §4. DOI: [10.1007/3-540-30666-8](https://doi.org/10.1007/3-540-30666-8)
- Burns, K.J. et al. (2020). "Dedalus: A flexible framework for numerical simulations with spectral methods." *Physical Review Research*, 2:023068. DOI: [10.1103/PhysRevResearch.2.023068](https://doi.org/10.1103/PhysRevResearch.2.023068)
- Raffelt, G. & Stodolsky, L. (1988). "Mixing of photons with low-mass particles in magnetic fields." *Physical Review D*, 37:1237. DOI: [10.1103/PhysRevD.37.1237](https://doi.org/10.1103/PhysRevD.37.1237)
- Dormand, J.R. & Prince, P.J. (1980). "A family of embedded Runge-Kutta formulae." *Journal of Computational and Applied Mathematics*, 6. DOI: [10.1016/0771-050X(80)90013-3](https://doi.org/10.1016/0771-050X(80)90013-3)
- Hindmarsh, A.C. et al. (2005). "SUNDIALS: Suite of nonlinear and differential/algebraic equation solvers." *ACM TOMS*, 31(3):363-396. DOI: [10.1145/1089014.1089020](https://doi.org/10.1145/1089014.1089020)
- Al-Mohy, A.H. & Higham, N.J. (2011). "Computing the action of the matrix exponential, with an application to exponential integrators." *SIAM J. Sci. Comput.*, 33(2):488-511. DOI: [10.1137/100788860](https://doi.org/10.1137/100788860)
- Trefethen, L.N. & Embree, M. (2005). *Spectra and Pseudospectra: The Behavior of Nonnormal Matrices and Operators*. Princeton University Press.
- Van Loan, C. (1978). "Computing integrals involving the matrix exponential." *IEEE Trans. Automatic Control*, 23(3):395-404.
- Higham, N.J. (2008). *Functions of Matrices: Theory and Computation*. SIAM. Ch. 12 (matrix cosine/sine).
