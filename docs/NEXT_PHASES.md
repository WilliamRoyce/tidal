# Next Major Implementation Phases for TIDAL

**Created:** February 2026
**Last Updated:** February 2026
**Status:** Phase A complete; Phases B–I planned
**Version:** 0.4.0 | **Tests:** 957+ collected | **Examples:** 22 working (1+1D to 3+1D)

## Context

TIDAL (Tensor Integration and Derivation for Any Lagrangian) has completed its core pipeline: Lagrangian (xAct/Mathematica) → JSON spec → PDE simulation (py-pde) → measurement/analysis. With 22 working examples spanning 1+1D to 3+1D, a full CLI (`tidal derive|simulate|measure|inspect|list|validate`), and a comprehensive measurement module (energy, conversion P(t), dispersion, mixing, spectral), the project is mature and ready for its next major advances.

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

See `docs/background_fields.md` for the full architecture documentation.

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

New presets are trivially added: write a `Build*GaugeTerm` function in `GaugeFix.wl` + add one entry to `_GAUGE_PRESETS` in `_derive.py`. See `docs/gauge_fixing.md` for a full tutorial and developer guide.

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
- Tutorial: `docs/gauge_fixing.md` — quick start, preset reference, custom expression walkthrough, "adding new presets" developer guide (includes inline TOML examples for Lorenz preset and custom expressions)

**B2: Additional presets + constraint mechanism** (~3–5 days)

- `GaugeFix.wl`: `BuildDeDonderGaugeTerm`, `BuildTemporalGaugeConstraint`, `BuildCoulombGaugeConstraint`, `BuildAxialGaugeConstraint`
- Type B WLS generation: `_wls_gauge_fixing_type_b()` (post-EOM constraint application)
- Constraint mechanism reuses existing `constraint_solver` infrastructure
- Update `examples/gravitational_waves/` with optional de Donder config
- Additional examples as appropriate

### Key Files

- **NEW** `tidal/wolfram/GaugeFix.wl` — Core primitive + preset builder functions
- **NEW** `docs/gauge_fixing.md` — Tutorial, preset reference, custom expression guide, developer recipe
- `tidal/cli/_derive.py` — `_GAUGE_PRESETS` registry, TOML validation, WLS generation
- `tidal/wolfram/ExportJSON.wl` — gauge metadata passthrough (already works)

### Scope: Medium (~7–10 days total across B1–B2)

### Dependencies: None

---

## Phase C: Parameter Sweep & Convergence Analysis

**Priority: HIGH — essential for systematic physics studies and publication-quality results**
**Status:** Planned

### What and Why

Currently each simulation is a single run with one parameter set and one grid resolution. There is no way to sweep coupling constants to map out P(g), or to verify numerical convergence by running at multiple resolutions. Reviewers expect convergence analysis in any computational physics paper.

Standard V&V methodology (Roache 1998; NASA GRC grid convergence tutorial; AIAA G-077-1998) requires demonstrating that numerical solutions converge to the continuum limit at the expected order. This phase implements those requirements as automated CLI commands.

### What It Enables

- **Parameter sweeps**: Conversion probability vs. coupling strength curves P_max(g), L_mix(g)
- **Richardson extrapolation**: Estimate the grid-independent solution from runs at 2–3 resolutions
- **Grid Convergence Index (GCI)**: Quantify numerical uncertainty per Roache's formulation
- **Method of Manufactured Solutions (MMS)**: Inject known analytic solutions to verify operator discretisation accuracy
- **Resolution-independence validation** for all measurement quantities
- **Automated multi-run workflows** with collected results

### Implementation Details

1. **CLI**: `tidal sweep spec.json --param "g=0.1:1.0:10"` and `tidal converge spec.json --grids "32,64,128,256"`
2. **New module** `tidal/cli/_sweep.py`: Orchestrates multiple simulation runs
3. **Post-processing**: Automatically calls measurement module on each run, collects scalar results into summary JSON
4. **Convergence**: Computes ‖u*h − u*{h/2}‖ norms, estimates convergence order p, computes GCI
5. **MMS module**: `tidal/verification/mms.py` — generates source terms from a prescribed analytic solution
6. **Output**: JSON summary + optional multi-panel plots (P vs. parameter, error vs. resolution)

### References

- Roache (1998), _Verification and Validation in Computational Science and Engineering_
- NASA GRC, "Examining Spatial (Grid) Convergence"
- AIAA G-077-1998, _Guide for V&V of CFD Simulations_

### Scope: Medium (~5–7 days)

### Dependencies: None (can be developed in parallel with other phases)

---

## Phase D: Coupled EM-Gravity Gertsenshtein Example

**Priority: HIGH — the culmination of the project's research goal**
**Status:** Planned

### What and Why

This is the integration example that combines Phase A (and optionally Phase B) into the first fully automated, Lagrangian-derived simulation of the Gertsenshtein effect. The repository is named "torsion-gertsenshtein" — this example is the raison d'être.

### What It Enables

- End-to-end Gertsenshtein simulation from a single `theory.toml`
- Validation against the analytical thin-magnet formula: P ≈ (κ B₀ L / 2)² (Domcke & Garcia-Cely 2023)
- Full measurement suite: P(t), P(k,t), mixing length, dispersion
- **Automated analytic benchmark test**: a pytest parametrised test comparing numerical P against the thin-magnet prediction at multiple coupling strengths, ensuring agreement to within 5% for weak coupling

### Implementation Details

1. **theory.toml**: EM field A_μ, linearised metric perturbation h_ab, background B₀, optional gauge fixing (Lorenz for EM, de Donder for gravity — simplifies equations but not required)
2. **Simulation**: Plane wave in A, zero in h, measure energy transfer A → h
3. **Validation**: Compare P(t) against the analytic Gertsenshtein formula at steady state
4. **Measurement**: Existing `compute_conversion_probability`, `compute_mixing_length`, `compute_spectral_conversion`
5. **Benchmark test**: `tests/test_gertsenshtein_benchmark.py` with parametrised coupling values

### References

- Gertsenshtein (1962), "Wave resonance of light and gravitational waves", JETP 14, 84
- Domcke & Garcia-Cely (2023), "A simple derivation of the Gertsenshtein effect", [arXiv:2301.02072](https://arxiv.org/abs/2301.02072) — thin-magnet formula
- Domcke & Garcia-Cely (2023), "On graviton-photon conversions in magnetic environments", [arXiv:2310.04150](https://arxiv.org/abs/2310.04150)
- Berlin et al. (2024), "Numerical analysis of resonant axion-photon mixing", [arXiv:2405.08865](https://arxiv.org/abs/2405.08865) — numerical methods for resonant mixing

### Scope: Medium (~3–5 days)

### Dependencies: Phase A (complete) required; Phase B optional (simplifies EM and gravity equations)

---

## Phase E: Spectral (Fourier) Spatial Discretisation

**Priority: MEDIUM — significant accuracy and performance improvement**
**Status:** Planned

### What and Why

TIDAL currently uses py-pde's 2nd-order finite-difference spatial discretisation. For wave propagation (the core use case), spectral methods (FFT-based) offer exponential convergence for smooth solutions on periodic domains. The Dedalus framework (Burns et al. 2020) demonstrates that spectral methods are the natural foundation for PDE solvers targeting wave physics.

### What It Enables

- Exponential spatial convergence for periodic wave problems (vs. O(h²) for FD)
- Much smaller grids for the same accuracy (e.g., 64 points vs. 512)
- Exact Laplacian in Fourier space (no numerical dispersion)
- Resolves the Dirichlet + cross_derivative energy drift issue

### Implementation Details

1. **New spatial backend**: `SpectralPDE` class computing derivatives via FFT
2. **CLI flag**: `--spatial spectral`
3. **All operators**: Spectral implementations (laplacian, gradient, cross_derivative)
4. **Chebyshev basis**: For non-periodic directions (mixed Fourier-Chebyshev, following Dedalus architecture)

### References

- Burns et al. (2020), "Dedalus: A Flexible Framework for Numerical Simulations with Spectral Methods", Phys. Rev. Research 2, 023068

### Scope: Large (~7–10 days)

### Dependencies: None

---

## Phase F: Adaptive Time-Stepping and Efficiency

**Priority: MEDIUM — required for production-quality long-duration runs**
**Status:** Planned

### What and Why

The Gertsenshtein conversion timescale can be vastly different from the wave oscillation timescale. Fixed-step integration either wastes time on the fast scale or under-resolves the slow scale. The existing `--scheme scipy` wraps `solve_ivp` but doesn't expose tolerance controls or method choices.

### What It Enables

- Efficient long-duration simulations with scale separation
- Stiff system support (large m²)
- Event-based stopping (e.g., stop when P exceeds a threshold)
- Publication-quality accuracy control

### Implementation Details

1. **CLI**: `--rtol 1e-8 --atol 1e-10 --method DOP853`
2. **Event detection**: Stop on P > threshold or energy conservation failure
3. **Snapshot interpolation**: Uniform-time resampling for FFT-based measurements
4. **Implicit methods**: BDF/Radau for stiff systems

### Scope: Small-Medium (~2–4 days)

### Dependencies: None

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

## Phase J: Constraint Pre-Solve (FFT-Based Initial Conditions)

**Priority: LOW — only needed for constraint systems with nontrivially violated ICs**
**Status:** Planned

### What and Why

When a DAE system has algebraic constraints (time_order=0 fields like A_0 in electromagnetism or Chern-Simons), the constraint field values must be consistent with the initial conditions of the dynamical fields. Currently, IDA handles this via `calc_initcond="yp0"`, which works when the constraint is trivially satisfied at t=0 (e.g., EM with zero initial momenta → Gauss's law gives A_0=0). However, it fails when the constraint has a nontrivial source:

- **Chern-Simons**: A_0 satisfies `laplacian(A_0) = kappa * d_y(A_1) - kappa * d_x(A_2)`. With a Gaussian in A_1, the RHS is nonzero, so A_0=0 is inconsistent. IDA's Newton solver fails because the Laplacian with periodic BCs has a zero eigenmode (singular Jacobian).
- **Any DAE** where the algebraic equation has a nontrivial source term at t=0.

### Proposed Solution

1. **FFT-based Poisson/Helmholtz solver**: Before passing ICs to IDA, solve each constraint equation spectrally. For periodic BCs, FFT naturally handles the zero mode by setting k=0 component to zero (unique up to a constant, which is unphysical for gauge fields).
2. **Integration with `solve_ida`**: Add an optional `pre_solve_constraints=True` parameter that runs the spectral solve before `calc_initcond`.
3. **CLI flag**: `--pre-solve-constraints` to enable this for nontrivial systems.

### References

- Standard FFT-based Poisson solvers; see e.g. Numerical Recipes (Press et al. 2007), Sec. 19.4
- Dedalus (Burns et al. 2020) uses spectral methods for constraint equations natively

### Scope: Small (~2–3 days)

### Dependencies: None (uses existing numpy/scipy FFT infrastructure)

---

## Known Limitations

1. **Chern-Simons IDA failure**: Systems where algebraic constraints are nontrivially violated at t=0 fail with `IDACalcIC - The line search failed`. Root cause: singular Laplacian with periodic BCs in the Jacobian. **Fix: Phase J** (constraint pre-solve). Workaround: choose initial conditions that trivially satisfy the constraint.

2. **Non-periodic BCs for constraint mode**: The `--mode constraint` path works with periodic BCs but may fail with Dirichlet/Neumann BCs for certain systems. **Future improvement**.

---

## Implementation Order

```
Phase A (Background Fields)      ─── COMPLETE
Phase D (Gertsenshtein Example)  ─── Requires A; B optional for cleaner equations
Phase B (Gauge Fixing, optional) ─── COMPLETE
Phase C (Sweep & Convergence)    ─── Independent, high priority
Phase F (Adaptive Time-Stepping) ─── Independent, quick win
Phase G (Absorbing Boundaries)   ─── Independent, uses Phase A infrastructure
Phase J (Constraint Pre-Solve)   ─── Independent, fixes Chern-Simons IDA failure
Phase H (HDF5/XDMF Output)      ─── Independent, interoperability
Phase I (Eigenvalue/Dispersion)  ─── Independent, analysis capability
Phase E (Spectral Methods)       ─── Independent, large scope
```

**Critical path to Gertsenshtein:** A (done), B (done) → D (~3–5 days)

**Recommended order for maximum impact:**

1. **D** (Gertsenshtein Example) — the project's raison d'être, unblocked by A+B
2. **C** (Sweep & Convergence) — required for publication-quality validation of D
3. **G** (Absorbing Boundaries) — extends D to realistic finite-magnet geometries
4. **F** (Adaptive Time-Stepping) — quick win for production runs
5. **I** (Eigenvalue/Dispersion) — analysis tool for parameter exploration
6. **H** (HDF5/XDMF Output) — interoperability with standard tools
7. **E** (Spectral Methods) — large scope, significant accuracy payoff

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
