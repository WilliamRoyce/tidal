# Next Major Implementation Phases for TIDAL

**Created:** February 2026
**Status:** Planning complete, Phase A in progress

## Context

TIDAL (Tensor Integration and Derivation for Any Lagrangian) has completed its core pipeline: Lagrangian (xAct/Mathematica) -> JSON spec -> PDE simulation (py-pde) -> measurement/analysis. With 1010 tests, 18 working examples spanning 1+1D to 3+1D, a full CLI, and a comprehensive measurement module (energy, conversion P(t), dispersion, mixing, spectral), the project is mature and ready for its next major advances.

The project's core research motivation is the **Gertsenshtein effect** (electromagnetic <-> gravitational wave conversion in external magnetic fields). The phases below are ordered by their impact toward enabling realistic Gertsenshtein simulations, while also broadening TIDAL's general utility as a field theory simulation framework.

---

## Phase A: Background Fields (External Source Terms)

**Priority: HIGHEST -- unlocks the fundamental Gertsenshtein physics**
**Status:** In progress

### What and Why

The Gertsenshtein effect requires a static external magnetic field B_0(x) that acts as a catalyst for photon-graviton conversion. Currently TIDAL only supports dynamical fields -- there is no mechanism for specifying a non-dynamical background field that survives as a spatially-varying coefficient in the EOM. This is the single most critical missing piece for the project's core research goal.

### What It Enables

- Lagrangians like `L = -1/4 F_ab F^ab + alpha * B_0(x) * F_ab * h^ab` where B_0 is prescribed
- Any field theory with external sources, fixed potentials, or background configurations
- Probe-field approximations (test field on fixed background)

### Implementation Details

1. **TOML**: New `[[background_fields]]` section with name, type, and profile expression
2. **Wolfram**: `DefTensor[B0[], M]` as non-dynamical (not varied in EL derivation -- `VarD` already only varies specified fields). Background field functions naturally survive as coefficients in the EOM
3. **ExportJSON**: Background field terms export as `coordinate_dependent` symbolic coefficients, reusing existing infrastructure (`"sin(2*Pi*x/L)"`, `"B0_val"`, etc.)
4. **Python**: Zero changes needed -- `_resolve_coefficient_at_point` already handles arbitrary position-dependent symbolic expressions

### Key Files

- `tidal/cli/_derive.py` -- TOML parsing and WLS code generation
- `tidal/wolfram/ExportJSON.wl` -- Term classification, coefficient extraction
- `tidal/wolfram/ComponentDecompose.wl` -- Component expansion (background fields as coefficients)
- `tidal/wolfram/EulerLagrange.wl` -- `VarD` only varies dynamical fields (background fields survive)

### Scope: Medium (~3-5 days)
### Dependencies: None (builds on existing position-dependent coefficient infrastructure)

---

## Phase B: Automatic Gauge Fixing

**Priority: HIGH -- required for any vector field theory to produce clean wave equations**
**Status:** Planned

### What and Why

Currently gauge fixing is done manually in WLS scripts: the electromagnetic example manually constructs the Lorenz-gauge wave equation instead of letting the pipeline derive and gauge-fix automatically. This doesn't scale to coupled systems like EM-gravity. Every Gertsenshtein simulation needs a gauge-fixed EM sector.

### What It Enables

- One-line `[gauge]` TOML config instead of manual Wolfram manipulation
- Lorenz gauge for EM: adds `L_gf = -(1/2xi)(div A)^2` to Lagrangian, reducing Maxwell to wave equations
- Coulomb gauge: A_0 = 0 + div(A) = 0 constraint (maps to existing constraint solver)
- De Donder (harmonic) gauge for linearized gravity: reduces 10-component h_ab to clean wave equations
- Scalable to any new vector/tensor field theory

### Implementation Details

1. **TOML**: `[gauge]` section with `type = "lorenz"|"coulomb"|"de_donder"`, `fields = ["A"]`
2. **New Wolfram module** `GaugeFix.wl`: `ApplyLorenzGaugeFix[L, field, covd, xi]` injects gauge-fixing term before EL derivation. `ApplyDeDonderGauge[L, h, metric, covd]` for gravitational perturbations
3. **`_derive.py`**: WLS generator reads `[gauge]` and emits the appropriate gauge-fixing call
4. **Validation**: Verify gauge-fixed Maxwell equations match the manually constructed EM example

### Scope: Medium-Large (~5-8 days, Lorenz + de Donder)
### Dependencies: None; de Donder is a natural extension of Lorenz

---

## Phase C: Parameter Sweep and Convergence Analysis

**Priority: HIGH -- essential for systematic physics studies and publication-quality results**
**Status:** Planned

### What and Why

Currently each simulation is a single run with one parameter set and one grid resolution. There is no way to sweep coupling constants to map out `P(g)`, or to verify numerical convergence by running at multiple resolutions. Reviewers expect convergence analysis in any computational physics paper.

### What It Enables

- Conversion probability vs. coupling strength curves: `P_max(g)`, `L_mix(g)`
- Convergence order verification via Richardson extrapolation
- Resolution-independence validation for all measurement quantities
- Automated multi-run workflows with collected results

### Implementation Details

1. **CLI**: `tidal sweep spec.json --param "g=0.1:1.0:10"` and `tidal converge spec.json --grids "32,64,128,256"`
2. **New module** `tidal/cli/_sweep.py`: Orchestrates multiple simulation runs
3. **Post-processing**: Automatically calls measurement module on each run, collects scalar results into summary JSON
4. **Convergence**: Computes `||u_h - u_{h/2}||` norms, estimates convergence order
5. **Output**: JSON summary + optional multi-panel plots

### Scope: Medium (~4-6 days)
### Dependencies: None (can be developed in parallel with other phases)

---

## Phase D: Nonlinear (Self-Interaction) Terms

**Priority: MEDIUM -- extends TIDAL from linear to general field theory**
**Status:** Planned

### What and Why

TIDAL currently handles only quadratic Lagrangians (linear PDEs). The JSON spec hardcodes `rhs.type = "linear_combination"`, and `IdentifyMultiFieldTerm` in ExportJSON.wl discards nonlinear terms. This prevents studying phi^4 theory, gravitational self-interaction, nonlinear wave steepening, and the strong-field Gertsenshtein regime.

### What It Enables

- phi^4 theory: `L = -1/2 (d phi)^2 - 1/2 m^2 phi^2 - lambda/4! phi^4`
- Kerr nonlinearity in optics, nonlinear wave phenomena
- Strong-field photon-graviton conversion (beyond linear approximation)

### Implementation Details

1. **JSON**: New operator types `"power_N"` for field^N terms
2. **Wolfram**: Modify `IdentifyMultiFieldTerm` to classify nonlinear terms by field occurrence count
3. **Python**: Pointwise nonlinear operators in `_OPERATOR_REGISTRY`
4. **Energy**: Explicit potential density integration (virial theorem doesn't hold for non-quadratic)
5. **CFL**: Nonlinear CFL estimation based on max field amplitude

### Scope: Large (~6-10 days)
### Dependencies: None; Phase A provides an alternative for weak-field Gertsenshtein

---

## Phase E: Coupled EM-Gravity Gertsenshtein Example

**Priority: HIGH -- the culmination of the project's research goal**
**Status:** Planned

### What and Why

This is the integration example that combines Phases A and B into the first fully automated, Lagrangian-derived simulation of the Gertsenshtein effect. The repository is literally named "torsion-gertsenshtein" -- this example is the raison d'etre.

### What It Enables

- End-to-end Gertsenshtein simulation from a single `theory.toml`
- Validation against the analytical formula: P ~ (alpha * B0 * L)^2
- Full measurement suite: P(t), P(k,t), mixing length, dispersion

### Implementation Details

1. **theory.toml**: EM field A_mu, linearized h_ab, background B0, gauge fixing
2. **Simulation**: Plane wave in A, zero in h, measure energy transfer
3. **Validation**: Compare P(t) against analytical Gertsenshtein formula
4. **Measurement**: Existing `compute_conversion_probability`, `compute_mixing_length`, `compute_spectral_conversion`

### Scope: Medium (~3-5 days)
### Dependencies: Phases A and B

---

## Phase F: Spectral (Fourier) Spatial Discretization

**Priority: MEDIUM -- significant accuracy and performance improvement**
**Status:** Planned

### What and Why

TIDAL currently uses py-pde's 2nd-order finite-difference spatial discretization. For wave propagation (the core use case), spectral methods (FFT-based) offer exponential convergence for smooth solutions on periodic domains. Comparable frameworks like Dedalus use spectral methods as their foundation.

### What It Enables

- Exponential spatial convergence for periodic wave problems (vs. O(h^2) for FD)
- Much smaller grids for the same accuracy (e.g., 64 points vs. 512)
- Exact Laplacian in Fourier space (no numerical dispersion)
- Resolves Dirichlet + cross_derivative energy drift issue

### Implementation Details

1. **New spatial backend**: `SpectralPDE` class computing derivatives via FFT
2. **CLI flag**: `--spatial spectral`
3. **All operators**: Spectral implementations (laplacian, gradient, cross_derivative)
4. **Dealiasing**: 2/3 rule for nonlinear terms (Phase D compatibility)

### Scope: Large (~7-10 days)
### Dependencies: None

---

## Phase G: Adaptive Time-Stepping and Efficiency

**Priority: MEDIUM -- required for production-quality long-duration runs**
**Status:** Planned

### What and Why

The Gertsenshtein conversion timescale can be vastly different from the wave oscillation timescale. Fixed-step integration either wastes time on the fast scale or under-resolves the slow scale. The existing `--scheme scipy` wraps `solve_ivp` but doesn't expose tolerance controls or method choices.

### What It Enables

- Efficient long-duration simulations with scale separation
- Stiff system support (large m^2)
- Event-based stopping
- Publication-quality accuracy control

### Implementation Details

1. **CLI**: `--rtol 1e-8 --atol 1e-10 --method DOP853`
2. **Event detection**: Stop on P > threshold or energy conservation failure
3. **Snapshot interpolation**: Uniform-time resampling for FFT-based measurements
4. **Implicit methods**: BDF/Radau for stiff systems

### Scope: Small-Medium (~2-4 days)
### Dependencies: None

---

## Implementation Order

```
Phase A (Background Fields)     ─┐
Phase B (Gauge Fixing)           ├─→ Phase E (Gertsenshtein Example)
                                 ─┘
Phase C (Parameter Sweep)        ─── Independent, can be done in parallel
Phase G (Adaptive Time-Stepping) ─── Independent, quick win
Phase D (Nonlinear Terms)        ─── Independent, lower priority for Gertsenshtein
Phase F (Spectral Methods)       ─── Independent, large scope
```

**Critical path to Gertsenshtein:** A -> B -> E (10-18 days total)

## Verification Plan

After each phase:
1. **Unit tests**: 15-30 new tests per phase, maintaining 0 ruff/pyright errors
2. **Integration tests**: End-to-end TOML -> JSON -> simulation -> measurement
3. **Physics validation**: Compare against analytical solutions where available
4. **Energy conservation**: Verify dE/E < 1e-6 for new features with periodic BCs
5. **Example parity**: New example with `theory.toml` + `run.sh`
6. **Documentation**: CHANGELOG entry, README update
