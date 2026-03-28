# Next Major Implementation Phases for TIDAL

**Created:** February 2026
**Last Updated:** March 2026
**Status:** Phases A, B, C, D, E (FFT), F, J complete; Torsion (PGT) complete (v0.18.0); Torsion-Gertsenshtein investigated (v0.22.8); Phases G–I planned
**Version:** 0.22.9 | **Tests:** 1,700 collected | **Examples:** 20 working (1+1D to 3+1D)

## Context

TIDAL (Tensor Integration and Derivation for Any Lagrangian) has completed its core pipeline: Lagrangian (xAct/Mathematica) → JSON spec → native PDE solver (SUNDIALS IDA/CVODE, leapfrog, scipy) → measurement/analysis. With 20 working examples spanning 1+1D to 3+1D, a full CLI with 9 subcommands (`tidal derive|simulate|measure|inspect|list|validate|plot|sweep|analyze`), a comprehensive measurement module (12 types: energy, conversion, mixing, spectrum, dispersion, conservation, effective_mass, asymptotic, peak_conversion, velocity, resonance, summary), and a complete parameter sweep framework with sensitivity analysis, the project is mature and ready for its next major advances.

The project's core research motivation is the **Gertsenshtein effect** (electromagnetic ↔ gravitational wave conversion in external magnetic fields). The project operates exclusively in the **linearised regime** — all Lagrangians are quadratic, producing linear PDEs. The phases below are ordered by their impact toward enabling realistic Gertsenshtein simulations, while also broadening TIDAL's general utility as a linearised field theory simulation framework.

Design decisions and feature choices below are informed by established scientific codebases — principally **Dedalus** (Burns et al. 2020), **MEEP** (Oskooi et al. 2010), and **FEniCS** (Baratta et al. 2023) — and by the Gertsenshtein-effect literature. Full citations are maintained in [`docs/references.md`](references.md).

---

## Phase A: Background Fields (External Source Terms) ✓

**Priority: HIGHEST — unlocks the fundamental Gertsenshtein physics**
**Status:** Complete

Background fields allow non-dynamical tensors (e.g., an external magnetic field B₀(x)) to appear in the Lagrangian without being varied in the Euler–Lagrange derivation. They survive as spatially-varying coefficients in the equations of motion.

**Key deliverables (all complete):**

- `[[background_fields]]` TOML section with scalar, vector, and tensor support
- Wolfram: `DefTensor` + `ReplaceAll` (scalar) / `ComponentValue` + `ToBasis` (vector/tensor)
- Python: 3-level caching (L0 preresolved → L1 expression → L2 spatial grid → L3 per-call)
- 3 working examples: `coupled_scattering/`, `proca_background/`, `vector_background/`

See `docs/tex/background_fields.tex` for the full architecture documentation.

---

## Phase B: Optional Gauge Fixing Toolkit

**Priority: MEDIUM — useful convenience for simplifying vector/tensor equations, not required**
**Status:** Complete

### What and Why

Gauge fixing simplifies equation structure for theories with gauge symmetry (massless vectors, linearized gravity). It is never required — TIDAL's existing pipeline handles gauge-invariant theories correctly, and all measurement quantities (energy, conversion, mixing) are gauge-invariant. However, explicit gauge fixing can be desirable to:

- Reduce coupled Maxwell equations to uncoupled wave equations (Lorenz gauge)
- Reduce 10-component linearized Einstein equations to clean wave equations (de Donder gauge)
- Eliminate unphysical degrees of freedom for cleaner simulations

Gauge fixing is **always opt-in and per-field**: a multi-field theory (A, B, h) can have different gauge choices for each field, or no gauge fixing at all.

### Architecture: Expression-Based Extensibility

The gauge system follows the same pattern as `[[derived_fields]]`: users can write **arbitrary Wolfram expressions** as gauge-fixing terms. Named gauges (Lorenz, de Donder, etc.) are built-in presets — convenience sugar over the same expression mechanism. Adding a new gauge preset requires only one new function in `GaugeFix.wl` and one registry entry in `_derive.py`.

**Two mechanisms:**

- **Type A (Lagrangian term):** An expression added to L before Euler-Lagrange derivation — changes EOM structure
- **Type B (Constraint):** A constraint imposed on the EOM after derivation — eliminates degrees of freedom

### TOML Configuration

```toml
# Named preset (convenience)
[[gauge]]
field = "A"
type = "lorenz"
xi = 1.0              # optional gauge parameter (default 1.0 = Feynman gauge)

# Custom Lagrangian term (full flexibility)
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
expression = "-(1/(2*xi)) * eta[a,b] CD[-a][A[-b]] eta[c,d] CD[-c][A[-d]]"

# Custom constraint
[[gauge]]
field = "A"
type = "custom"
mechanism = "constraint"
expression = "eta[a,b] CD[-a][A[-b]]"   # set to zero
```

### Built-In Presets

| Preset      | Mechanism       | Fields      | Expression                     | Effect                             |
| ----------- | --------------- | ----------- | ------------------------------ | ---------------------------------- |
| `lorenz`    | lagrangian_term | vector      | `-(1/2ξ)(∂_μ A^μ)²`            | Maxwell → uncoupled wave equations |
| `de_donder` | lagrangian_term | sym. rank-2 | `-(1/2ξ)(∂_a h^a_b - ½∂_b h)²` | Lin. Einstein → uncoupled waves    |
| `temporal`  | constraint      | vector      | `A_0 = 0`                      | Eliminates temporal component      |
| `coulomb`   | constraint      | vector      | `∇·A = 0`                      | Transversality constraint          |
| `axial`     | constraint      | vector      | `A_n = 0`                      | Eliminates one spatial component   |

New presets are trivially added: write a `Build*GaugeTerm` function in `GaugeFix.wl` + add one entry to `_GAUGE_PRESETS` in `_derive.py`. See `docs/tex/gauge_fixing.tex` for a full tutorial and developer guide.

### Implementation Sub-Phases

**B1: Core framework + Lorenz proof-of-concept** (~4–5 days)

- Expression-based `[[gauge]]` TOML parsing + `_validate_gauge()` in `_derive.py`
- `_GAUGE_PRESETS` registry (extensible dict mapping names → builder functions)
- `_WlsContext.gauge` field
- `GaugeFix.wl`: `AddGaugeFixingTerm` (core primitive) + `BuildLorenzGaugeTerm` (first preset)
- WLS generation: `_wls_gauge_fixing_type_a()` handles both custom expressions and named presets
- Custom expression path reuses `_substitute_field_names()` (existing infrastructure)
- Dynamic gauge metadata in `_wls_metadata_and_export()`
- Tests: validation, WLS generation, custom expression handling, Wolfram unit tests
- Tutorial: `docs/tex/gauge_fixing.tex` — quick start, preset reference, custom expression walkthrough, "adding new presets" developer guide (includes inline TOML examples for Lorenz preset and custom expressions)

**B2: Additional presets + constraint mechanism** (~3–5 days)

- `GaugeFix.wl`: `BuildDeDonderGaugeTerm`, `BuildTemporalGaugeConstraint`, `BuildCoulombGaugeConstraint`, `BuildAxialGaugeConstraint`
- Type B WLS generation: `_wls_gauge_fixing_type_b()` (post-EOM constraint application)
- Constraint mechanism reuses existing `constraint_solver` infrastructure
- Update `examples/gravitational_waves/` with optional de Donder config
- Additional examples as appropriate

### Key Files

- **NEW** `tidal/wolfram/GaugeFix.wl` — Core primitive + preset builder functions
- **NEW** `docs/tex/gauge_fixing.tex` — Tutorial, preset reference, custom expression guide, developer recipe
- `tidal/cli/_derive.py` — `_GAUGE_PRESETS` registry, TOML validation, WLS generation
- `tidal/wolfram/ExportJSON.wl` — gauge metadata passthrough (already works)

### Scope: Medium (~7–10 days total across B1–B2)

### Dependencies: None

---

## Phase C: Parameter Sweep & Convergence Analysis ✓

**Priority: HIGH — essential for systematic physics studies and publication-quality results**
**Status:** Complete

### What Was Delivered

A comprehensive parameter sweep framework implemented as two CLI commands (`tidal sweep` and `tidal analyze`) with 8 feature areas (F1–F8). See [`docs/next-features.md`](next-features.md) for the full feature specification.

**CLI Commands:**

- `tidal sweep spec.json --sweep "g=0.1:1.0:10" --measure conversion` — unified command for parameter sweeps AND convergence studies (via `--converge "32,64,128,256"` flag)
- `tidal analyze sweep_dir/ --sensitivity sobol --metric P_max` — post-hoc Sobol/Morris sensitivity analysis

**Features delivered:**

- **F1: TOML sweep configuration** — `--config sweep.toml` for reproducible, version-controlled sweep definitions
- **F2a: Adaptive sampling** — `--adaptive-metric`, `--adaptive-budget`, `--adaptive-threshold` for automatic refinement in interesting parameter regions
- **F2b: Latin Hypercube / Sobol sampling** — `--sweep-strategy latin_hypercube|sobol`, `--n-samples N` for multi-dimensional parameter spaces
- **F3: Velocity + resonance analysis** — group/phase velocity mismatch and resonance condition detection
- **F4: Sobol/Morris sensitivity analysis** — first-order, total-order, and interaction indices via SALib
- **F5: SweepResults query methods** — programmatic access to sweep data with filtering and aggregation
- **F6: Spectrum scalar aggregation** — all 12 measurement types supported in sweeps (including spectrum via scalar summaries)
- **F7: Run status tracking + resume** — `--resume` for interrupted sweeps, failure classification
- **F8: Advanced visualization** — 6 plot types: sweep (auto 1D/2D/multi), sweep-compare, convergence, sweep-parallel, sweep-tornado, sweep-scatter
- **Convergence mode** — Richardson extrapolation for convergence order estimation

**12 measurement types in sweeps:** summary, energy, conversion, mixing, spectrum, dispersion, conservation, effective_mass, asymptotic, peak_conversion, velocity, resonance

**7 working example scripts** in `coupled_scattering/`, `coupled_scalars/`, `scalar_field/`

**Key files:** `tidal/cli/_sweep.py` (~1350 lines), `tidal/cli/_analyze.py`, `tidal/cli/_sweep_panels.py`, `tidal/cli/_sweep_config.py`, `tidal/measurement/_velocity.py`, `tidal/measurement/_resonance.py`, `tidal/measurement/_sensitivity.py`

**Not implemented (deferred):** Grid Convergence Index (GCI) per Roache's formulation, Method of Manufactured Solutions (`tidal/verification/mms.py` — planned but not built), explicit ‖u_h − u_{h/2}‖ error norms.

### References

- Roache (1998), _Verification and Validation in Computational Science and Engineering_
- SALib (Herman & Usher), Sensitivity Analysis Library for Sobol/Morris methods

### Dependencies: None

---

## Phase D: Coupled EM-Gravity Gertsenshtein Example ✓

**Priority: HIGH — the culmination of the project's research goal**
**Status:** Complete (v0.22.8). Torsion-independence proven for minimal PGT.

### What and Why

This is the integration example that combines Phase A (and optionally Phase B) into the first fully automated, Lagrangian-derived simulation of the Gertsenshtein effect. The repository is named "torsion-gertsenshtein" — this example is the raison d'être.

### Delivered

- **Pipeline extension**: `[[linearization.matter_perturbations]]` in TOML — uses xPert's `DefTensorPerturbation` for matter fields alongside `SetupMetricPerturbation`
- **End-to-end simulation** from `examples/gertsenshtein/theory.toml` — graviton-photon conversion via Einstein-Maxwell Lagrangian
- **Uniform B₀ validation (Phase E)**: P = sin²(κB₀t/2) confirmed via 40-point B₀ sweep (N=1024, RMS < 0.012). Corrected P&R (2023) error: missing √(4π) in coupling. Confirmed by Dandoy/Lella (arXiv:2406.17853).
- **Localized B-field validation (Phase F2)**: Gaussian B_x(z) via `theory_localized.toml`. Boccaletti formula P = sin²(κ/2 × ∫B dz) validated: P_numerical = 0.3436 vs P_Boccaletti = 0.3432 (0.04% agreement). 48-point sweep max error < 0.003.
- **Gauge-independence validation**: Ungauged EM-only theory (`theory_ungauged.toml`, 14 components) confirms P_peak = 0.997, identical to TT-gauged result. Active channel: h_5 ↔ a_1 (h× ↔ a_x).

### Torsion-Gertsenshtein Investigation (v0.22.8)

Combined PGT + Einstein-Maxwell theory derived (`examples/torsion_gertsenshtein/theory.toml`, 23 components after Ostrogradsky reduction). Key findings:

- **Polarisation block-diagonal structure (#200):** The plane-wave reduced equations decompose into two completely decoupled channels with zero cross-talk:
  - **h× ↔ a_x** (h_5 ↔ a_1): torsion-independent, stable, P = sin²(κB₀D/2) unchanged
  - **trace ↔ a_y** (h_4/h_7/h_9 ↔ a_2 ↔ torsion): torsion-dependent, but ghost-unstable (Ostrogradsky from R̃²)
- **Torsion-independence of standard Gertsenshtein (#199):** The h× ↔ a_x equations contain NO torsion parameters (α₁, α₂, α₃, b₅). Algebraically identical to EM-only.
- **Double inaccessibility:** The torsion-sensitive trace channel is both (1) algebraically unreachable from TT initial conditions and (2) ghost-unstable from Ostrogradsky.
- **Coupling-instability tension:** The b₅R̃² term simultaneously enables propagating torsion AND introduces Ostrogradsky ghosts — a structural consequence of Ostrogradsky's theorem.
- **Field-filtered energy measurement:** Added `fields` parameter to energy computation chain, allowing conversion measurement on Ostrogradsky theories by evaluating only source/target field Hamiltonian terms.
- **PGT formulation verified correct:** TIDAL varies L independently w.r.t. h, a, and t via VarD — the standard metric-affine/Palatini approach for the torsion sector.

### Remaining / Blocked

- **Radial dipolar (Phase F3a)**: `theory_radial.toml` derived (spherical coords). BLOCKED on compute.
- **Plasma detuning (Phase F1)**: BLOCKED — xPert spurious z²-terms from background 4-potential.
- **Non-minimal torsion-EM coupling:** Required to break the polarisation block-diagonal structure. T·F, T·(*F) terms would create direct photon-torsion interaction. Future work.
- **Ghost-free parameter conditions:** Literature (Sezgin & van Nieuwenhuizen 1980, Nikiforova et al. 2009, Barker 2024) provides sector-specific conditions. No universal closed-form for general (α₁, α₂, α₃, b₅).

### References

- Gertsenshtein (1962), "Wave resonance of light and gravitational waves", JETP 14, 84
- Domcke & Garcia-Cely (2023), "A simple derivation of the Gertsenshtein effect", [arXiv:2301.02072](https://arxiv.org/abs/2301.02072) — thin-magnet formula
- Domcke & Garcia-Cely (2023), "On graviton-photon conversions in magnetic environments", [arXiv:2310.04150](https://arxiv.org/abs/2310.04150)
- Dandoy & Lella (2024), "Graviton-photon oscillations", [arXiv:2406.17853](https://arxiv.org/abs/2406.17853) — confirms correct coupling normalization
- Berlin et al. (2024), "Numerical analysis of resonant axion-photon mixing", [arXiv:2405.08865](https://arxiv.org/abs/2405.08865)
- Sezgin & van Nieuwenhuizen (1980), "New ghost-free gravity Lagrangians", Phys. Rev. D 21:3269 — linearised PGT ghost conditions
- Nikiforova et al. (2009), "Stability of the massive torsion modes", [arXiv:0905.4007](https://arxiv.org/abs/0905.4007)
- Barker (2024), "Every Poincaré gauge theory is conformal", [arXiv:2406.12826](https://arxiv.org/abs/2406.12826) — no-ghost β₃ > 0

### Dependencies: Phase A (complete) required; Phase B (complete) simplifies EM and gravity equations

---

## Phase E: Spectral (Fourier) Spatial Discretisation ✓

**Priority: MEDIUM — significant accuracy and performance improvement**
**Status:** Complete (FFT operators); Chebyshev remaining

### What and Why

TIDAL supports FFT spectral operators (`--spectral`) alongside finite-difference stencils (2nd/4th/6th order via `--fd-order`). For wave propagation on periodic domains, spectral methods offer exponential convergence for smooth solutions. Auto-enabled when all boundary conditions are periodic.

### Delivered

- FFT spectral operators in `tidal/solver/operators.py` — machine-precision accuracy for smooth fields
- CLI flag: `--spectral` (auto-enabled for all-periodic BCs)
- All operators: spectral implementations (laplacian, gradient, cross_derivative)
- Higher-order FD stencils (`--fd-order 4|6`) via Fornberg (1988) coefficients — 500x error reduction per order doubling

### Remaining

- **Chebyshev basis** for non-periodic directions (mixed Fourier-Chebyshev, following Dedalus architecture)

### References

- Burns et al. (2020), "Dedalus: A Flexible Framework for Numerical Simulations with Spectral Methods", Phys. Rev. Research 2, 023068
- Fornberg (1988), "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids", Math. Comp. 51(184), 699–706

---

## Phase F: Adaptive Time-Stepping and Efficiency ✓

**Priority: MEDIUM — required for production-quality long-duration runs**
**Status:** Complete

Tolerance-controlled adaptive time-stepping via SUNDIALS CVODE (BDF) and IDA (DAE) replaces manual `dt` selection with error-controlled solvers (Hindmarsh et al. 2005). Users specify accuracy targets via `--rtol`/`--atol` instead of guessing a stable `dt`. This eliminates the fundamental fragility of fixed-step integration: a `dt` that satisfies CFL for one parameter set can be unstable for another, producing plausible but incorrect output. Berlin et al. (2024, arXiv:2405.08865) document these numerical challenges for resonant mixing simulations.

**Key deliverables (all complete):**

- Four solver paths: **CVODE** (adaptive BDF ODE), **IDA** (implicit DAE), **scipy** (`solve_ivp` with DOP853/Radau/BDF), **leapfrog** (Störmer-Verlet symplectic)
- CLI: `--scheme cvode|ida|scipy|leapfrog`, `--rtol`, `--atol`
- Automatic solver selection: systems with algebraic constraints → IDA; pure wave equations → CVODE or leapfrog
- Snapshot interpolation for non-uniform timesteps from adaptive solvers
- Smart CFL estimation for leapfrog; tolerance control for CVODE/IDA
- Sparse Jacobian computation for IDA to prevent out-of-memory on large systems

See `docs/tex/adaptive_timestepping.tex` for the full architecture documentation.

### References

- Hindmarsh et al. (2005), "SUNDIALS: Suite of Nonlinear and Differential/Algebraic Equation Solvers", ACM TOMS 31(3)
- Hairer, Lubich & Wanner (2006), _Geometric Numerical Integration_, Springer, 2nd ed.

---

## Phase G: Absorbing Boundaries (PML / Sponge Layers)

**Priority: MEDIUM-HIGH — needed for finite interaction region Gertsenshtein simulations**
**Status:** Planned

### What and Why

Periodic boundary conditions cannot model a finite interaction region (e.g., a magnet of length L in an otherwise infinite domain). Waves reaching the boundary re-enter the domain, contaminating the signal. Absorbing boundary layers damp outgoing waves without reflection, enabling open-domain simulations.

The Perfectly Matched Layer (PML) technique (Bérenger 1994) is the gold standard for absorbing boundaries in wave simulations, used extensively in MEEP (Oskooi et al. 2010) and other FDTD codes. A simpler alternative — the **sponge layer** — can be implemented as an exponentially ramped dissipation term, which maps naturally onto TIDAL's existing background-field coefficient infrastructure.

### What It Enables

- Open-domain wave propagation (outgoing waves absorbed, not reflected)
- Finite-magnet Gertsenshtein simulations: B₀(x) non-zero only in [x₁, x₂], with absorbers outside
- Scattering problems: incoming plane wave, measure transmitted/reflected amplitudes
- Direct comparison with analytic Gertsenshtein predictions for finite interaction lengths

### Implementation Details

1. **Sponge layer (Phase G.1)**: Add a dissipation term `−σ(x) ∂_t u` where σ(x) ramps from 0 inside the domain to σ_max in the absorbing region. Implementable as a `coordinate_dependent` coefficient in the JSON spec — no new Python operators needed
2. **TOML config**: `[absorbing_boundary]` section with `type = "sponge"`, `width`, `strength`, `profile = "quadratic"|"cubic"`
3. **Wolfram**: `_derive.py` injects the dissipation term into the EOM before export
4. **Full PML (Phase G.2)**: Complex coordinate stretching in the frequency domain; requires split-field formulation. Larger scope, deferred to G.2
5. **Validation**: Measure reflection coefficient R(ω) and verify R < 10⁻⁴ for well-resolved frequencies

### References

- Bérenger (1994), "A perfectly matched layer for the absorption of electromagnetic waves", J. Computational Physics 114, 185–200
- Johnson (2007), "Notes on Perfectly Matched Layers (PMLs)", MIT
- Oskooi et al. (2010), "MEEP: A flexible free-software package for electromagnetic simulations by the FDTD method"

### Scope: Medium (~4–6 days for sponge layer; +5–8 days for full PML)

### Dependencies: Phase A (complete) provides the `coordinate_dependent` coefficient infrastructure

---

## Phase H: HDF5/XDMF Output

**Priority: MEDIUM — publication-standard data export and interoperability**
**Status:** Planned

### What and Why

TIDAL currently uses raw numpy memory-mapped arrays for disk-backed storage. While efficient, this format is not interoperable with standard scientific visualisation tools (ParaView, VisIt, yt). The HDF5 + XDMF combination is the de facto standard for PDE simulation output — used by Dedalus, FEniCS, and most production codes. Adopting it enables direct post-processing in established toolchains without format conversion.

### What It Enables

- Direct loading in ParaView, VisIt, and yt for 3D volume rendering and slicing
- Self-describing data with metadata (grid info, time stamps, field names, parameters)
- Efficient parallel I/O for future MPI-parallel extensions
- On-the-fly analysis tasks (following Dedalus's analysis framework model)

### Implementation Details

1. **New storage backend**: `HDF5Storage` class in `tidal/measurement/_io.py` alongside existing memmap
2. **XDMF descriptor**: Auto-generated `.xdmf` file describing the HDF5 data layout for ParaView/VisIt
3. **CLI flag**: `--output-format hdf5` (default remains memmap for backwards compatibility)
4. **Metadata**: Store JSON spec, parameters, git hash, creation time in HDF5 attributes
5. **Migration**: `tidal convert output_dir/ --to hdf5` for existing simulation data

### References

- XDMF, "XDMF Model and Format", [xdmf.org](https://www.xdmf.org/index.php/XDMF_Model_and_Format)
- Burns et al. (2020), Dedalus HDF5 analysis output framework

### Scope: Medium (~4–6 days)

### Dependencies: None

---

## Phase I: Eigenvalue / Dispersion Solver

**Priority: MEDIUM — identifies propagating modes without full time-domain simulation**
**Status:** Planned

### What and Why

For the torsion PGT (Poincaré Gauge Theory) research goal, it is essential to identify which parameter windows support propagating (non-tachyonic, non-ghost) modes before committing to expensive time-domain simulations. An eigenvalue solver computes the dispersion relation ω(k) directly from the linearised equation system, revealing mode speeds, stability boundaries, and resonance conditions.

Dedalus (Burns et al. 2020) provides a native eigenvalue problem (EVP) capability that has proven invaluable for hydrodynamic stability analysis. TIDAL's linearised structure — mass matrices, coupling matrices, and spatial operators already extracted — is well-suited for a similar capability.

### What It Enables

- **Dispersion relations** ω(k) for all field components from the JSON spec alone
- **Stability analysis**: Detect tachyonic modes (ω² < 0) and ghosts (wrong-sign kinetic term) at the linear algebra level
- **Parameter window scanning**: Sweep coupling constants and identify viable regions before simulation
- **Group/phase velocity computation**: dω/dk for wave packet propagation predictions
- **Resonance identification**: Find k-values where mode speeds match (critical for Gertsenshtein conversion)

### Implementation Details

1. **New module** `tidal/analysis/dispersion.py`: Constructs the generalised eigenvalue problem (M − ω² I) u = 0 from the JSON spec's mass and coupling matrices plus spatial operator eigenvalues
2. **CLI**: `tidal dispersion spec.json --k-range "0:10:100" --param "m2=1.0"`
3. **Output**: ω(k) curves as JSON + optional matplotlib plots
4. **Validation**: Compare ω(k) against known analytic results (e.g., ω² = k² + m² for Klein-Gordon)

### References

- Burns et al. (2020), "Dedalus: A Flexible Framework for Numerical Simulations with Spectral Methods" — native EVP capability

### Scope: Medium (~4–6 days)

### Dependencies: None (uses existing mass/coupling matrix infrastructure from Phase 12)

---

## Phase J: Constraint Pre-Solve ✓

**Priority: LOW → Completed**
**Status:** Complete

### What and Why

When a DAE system has algebraic constraints (time_order=0 fields like A_0 in electromagnetism or Chern-Simons), the constraint field values must be consistent with the initial conditions of the dynamical fields. Previously IDA failed when constraints had nontrivial sources (e.g., Chern-Simons with Gaussian IC on A_1).

### Implemented Solution: Three-Tier Solver Architecture

**`tidal/solver/constraint_solve.py`** — automatic constraint pre-solve before IDA starts:

1. **Tier 1: FFT (O(N log N))** — for periodic BCs with constant-coefficient self-operators. Uses modified wavenumbers `k_mod = (2/dx)·sin(k·dx/2)` for exact FD-consistency. Handles coupled constraints (block system at each wavenumber).
2. **Tier 2: Operator Probing → Sparse Matrix (O(N²) build, O(N) solve)** — universal fallback for non-periodic BCs, position-dependent self-coefficients, or unknown operators. Probes `apply_operator()` with unit vectors to build the exact operator matrix.
3. **Automatic selection**: `_select_method()` checks periodicity, coefficient constancy, and multiplier availability to choose the fastest applicable tier.

**Gauge regularisation for singular Poisson (periodic BCs):**
- Pure-Laplacian constraints with periodic BCs have a null space (constant functions)
- FFT: sets zero-mode `u_hat[0,...,0] = 0` (zero-mean gauge)
- IDA: auto-detects via `_is_pure_laplacian()` and pins one DOF (`A_0[0] = 0`) in the residual
- This is numerical regularisation, not physics gauge fixing — observables (E, B) depend on derivatives of A_0, not A_0 itself (standard FEniCS/Firedrake/PETSc practice)

**Key files:**
- **NEW** `tidal/solver/constraint_solve.py` (~450 lines) — three-tier solver
- **MODIFIED** `tidal/solver/ida.py` — pre-solve integration + gauge regularisation
- **NEW** `tests/test_solver_constraint_solve.py` — 25 tests (unit + integration + IDA)

### References

- Standard FFT-based Poisson solvers; see e.g. Numerical Recipes (Press et al. 2007), Sec. 19.4
- Dedalus (Burns et al. 2020) uses spectral methods for constraint equations natively
- FEniCS/Firedrake null-space handling for Poisson with Neumann/periodic BCs

---

## Known Limitations

1. ~~**Chern-Simons IDA failure**~~: **Resolved by Phase J.** Constraint pre-solve + gauge regularisation handles all DAE systems, including those with nontrivially violated algebraic constraints and singular Laplacian Jacobians.

2. **Non-periodic BCs for constraint mode**: The `--mode constraint` path works with periodic BCs but may fail with Dirichlet/Neumann BCs for certain systems. Phase J's Tier 2 (operator probing) supports non-periodic BCs for the pre-solve step.

---

## Implementation Order

```
Phase A (Background Fields)      ─── COMPLETE
Phase B (Gauge Fixing, optional) ─── COMPLETE
Phase C (Sweep & Convergence)    ─── COMPLETE
Phase F (Adaptive Time-Stepping) ─── COMPLETE
Phase J (Constraint Pre-Solve)   ─── COMPLETE
Phase D (Gertsenshtein Example)  ─── Requires A; B optional for cleaner equations
Phase G (Absorbing Boundaries)   ─── Independent, uses Phase A infrastructure
Phase H (HDF5/XDMF Output)      ─── Independent, interoperability
Phase I (Eigenvalue/Dispersion)  ─── Independent, analysis capability
Phase E (Spectral Methods)       ─── Independent, large scope
```

**Critical path to Gertsenshtein:** A (done), B (done), C (done), F (done) → D (~3–5 days)

**Recommended order for maximum impact:**

1. **D** (Gertsenshtein Example) — the project's raison d'être, unblocked by A+B+C+F
2. **G** (Absorbing Boundaries) — extends D to realistic finite-magnet geometries
3. **I** (Eigenvalue/Dispersion) — analysis tool for parameter exploration
4. **H** (HDF5/XDMF Output) — interoperability with standard tools
5. **E** (Spectral Methods) — large scope, significant accuracy payoff

---

## Verification Plan

After each phase:

1. **Unit tests**: 15–30 new tests per phase, maintaining 0 ruff/pyright errors
2. **Integration tests**: End-to-end TOML → JSON → simulation → measurement
3. **Physics validation**: Compare against analytical solutions where available
4. **Energy conservation**: Verify dE/E < 10⁻⁶ for new features with periodic BCs
5. **Analytic benchmarks**: For Phase D specifically, automated comparison against the thin-magnet Gertsenshtein formula (Domcke & Garcia-Cely 2023)
6. **Convergence verification**: For Phase C, demonstrate expected convergence order (2nd for FD, exponential for spectral) using GCI methodology (Roache 1998)
7. **Example parity**: New example with `theory.toml` + `run.sh`
8. **Documentation**: CHANGELOG entry, README update
9. **Citations**: Document which external codebases/papers informed design decisions
