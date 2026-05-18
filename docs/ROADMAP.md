# TIDAL Pipeline Roadmap

This document outlines the planned improvements and features for the TIDAL symbolic physics pipeline project.

**Last Updated:** April 2026 (v0.33.9)
**Project Status:** Phase 13+ Complete; perturbative reduction v6 complete (Stage 7 closed 2026-04-20, issue #271 resolved). Euler–Heisenberg quartic-EM path fully supported; matter-only derivative-only theories now unblocked. 1,991 Python tests + ~115 Wolfram tests.

## Overview

The project is in a mature state with a robust symbolic pipeline from Lagrangian to PDE simulation. This roadmap tracks refinements toward production-grade reliability and extensibility.

**Total Identified Improvements:** 25 issues
**Estimated Effort:** 100-150 hours (6-8 sprint cycles)

---

## Implementation Phases

### Phase 1: Critical Fixes ✅ COMPLETE

**Goal:** Address security concerns and high-impact bugs

| Issue                                                       | Priority    | Type     | Status  |
| ----------------------------------------------------------- | ----------- | -------- | ------- |
| [#67] Replace Assertions with Explicit Error Handling       | 🔴 Critical | Bug      | ✅ Done |
| [#68] Strengthen Mathematica Expression Evaluation Security | 🔴 Critical | Security | ✅ Done |
| [#75] Validate Grid Dimensions During PDE Construction      | 🟠 High     | Bug      | ✅ Done |

**Delivered:**

- All `assert isinstance()` replaced with explicit `TypeError` raises (Issue #67)
- `_validate_eval_result` checks for NaN/Inf/complex (Issue #68)
- `_validate_operator_dimensions()` at construction time (Issue #75)

---

### Phase 2: Testing & CI Infrastructure — 🔄 PRIMARY REMAINING FOCUS

**Goal:** Improve test coverage and CI reliability
**Status:** Partially complete — coverage runs in CI, but key gaps remain

| Issue                                                      | Priority    | Type    | Status       |
| ---------------------------------------------------------- | ----------- | ------- | ------------ |
| [#TBD] Add Animation Module Test Coverage                  | 🟠 High     | Testing | 🔄 Remaining |
| [#74] Add Code Coverage Reporting to CI                    | 🟠 High     | CI/CD   | ✅ Done      |
| [#TBD] Add Wolfram Tests to GitHub Actions CI              | 🔴 Critical | CI/CD   | 🔄 Remaining |
| [#TBD] Add Tests for Observers, Profiling, Runners Modules | 🟡 Medium   | Testing | 🔄 Remaining |

**Delivered so far:**

- Coverage runs in CI (`--cov=tidal --cov-report=term-missing --cov-report=xml`)
- Coverage uploaded to Codecov via `codecov/codecov-action@v5` (PR #74)
- Codecov badge in README (links to Codecov dashboard)
- Root Makefile with `make test-coverage` target

**Still needed:**

- `test_animation_builder.py` with 10-15 test cases
- codecov.io integration for dynamic badge and PR reports
- Optional Wolfram test workflow (label-triggered or weekly)
- ~20-30 tests for observers, profiling, runners modules

---

### Phase 3: Features & Documentation ✅ COMPLETE

**Goal:** Validate 3+1D support and improve user documentation

| Issue                                             | Priority    | Type          | Status   |
| ------------------------------------------------- | ----------- | ------------- | -------- |
| [#71] Add 3+1D Spacetime Examples                 | 🔴 Critical | Documentation | ✅ Done  |
| [#TBD] Document JSON Schema with Detailed Guide   | 🟠 High     | Documentation | ✅ Done  |
| [#TBD] Improve Coefficient Resolution Performance | 🟠 High     | Performance   | ✅ Done  |
| [#TBD] Create Architecture Diagrams               | 🟡 Medium   | Documentation | Deferred |

**Delivered:**

- 3+1D Klein-Gordon example (Issue #71), plus spherical_kg, cylindrical_kg, gravitational_waves, massive_3form
- `docs/tex/json_schema.tex` with complete field reference
- Unified `_resolve_coefficient_at_point` evaluator (Phase 10b)

---

### Phase 4: Advanced Features ✅ MOSTLY COMPLETE

**Goal:** Implement major feature requests

| Issue                                               | Priority    | Type          | Status             |
| --------------------------------------------------- | ----------- | ------------- | ------------------ |
| [#70] Support Rank-3+ Tensor Decomposition          | 🔴 Critical | Feature       | ✅ Done (Phase 13) |
| [#TBD] Implement Automatic Gauge Fixing             | 🔴 Critical | Feature       | ✅ Done (Phase B)  |
| [#79] Handle Mixed Time-Space Cross-Derivatives     | 🟡 Medium   | Bug           | ✅ Done            |
| [#TBD] Expand \_mathematica_to_python Function Set  | 🟡 Medium   | Feature       | ✅ Done            |
| [#TBD] Add Non-Cartesian Coordinate System Examples | 🟡 Medium   | Documentation | ✅ Done            |
| [#TBD] Add Convergence and Stability Stress Tests   | 🟡 Medium   | Testing       | ✅ Done            |
| [#TBD] Add Full Pipeline Validation to CI           | 🟡 Medium   | CI/CD         | Remaining          |

**Delivered:**

- Rank-3+ tensor support: `ReplaceHigherRankFieldComponents`, massive_3form example (Phase 13)
- Mixed time-space derivatives: `ClassifySpatialProfile`, `ExtractSpatialOperatorFromMixed` (Issue #79)
- Curvilinear examples: polar_kg, spherical_kg, cylindrical_kg with Christoffel auto-detection
- CFL `check_stability()`, physics validation (energy conservation + analytical solutions)

**Remaining:**

- Wolfram tests in GitHub Actions CI

---

### Phase 5: Polish & Optimization (Partially Complete)

**Goal:** Code quality improvements and nice-to-have features

| Issue                                                      | Priority | Type          | Status    |
| ---------------------------------------------------------- | -------- | ------------- | --------- |
| [#85] Refactor derivative classification                   | 🟢 Low   | Refactoring   | ✅ Done   |
| [#TBD] Add Parameter Sweep Examples                        | 🟢 Low   | Documentation | ✅ Done   |
| [#TBD] Add Python 3.12+ Testing to CI Matrix               | 🟢 Low   | CI/CD         | Remaining |
| [#TBD] Support Elliptic PDE Solving (Constraint Equations) | 🟢 Low   | Feature       | ✅ Done   |

**Delivered:**

- Unified `ExtractDerivativeProfile` function (Issue #85, ~108 lines saved)
- Elliptic PDE solving: constraint equations with `time_derivative_order=0`, `--mode constraint` CLI flag

**Remaining:**

- Python 3.12 CI matrix

---

### Additional Completed Work (not in original roadmap)

| Feature                                                                                     | Status      |
| ------------------------------------------------------------------------------------------- | ----------- |
| CLI (`tidal` command) — 9 subcommands (derive, inspect, simulate, measure, list, validate, plot, sweep, analyze) | ✅ Complete |
| Measurement module — 12 types: energy, conversion, mixing, spectrum, dispersion, conservation, effective_mass, asymptotic, peak_conversion, velocity, resonance, summary | ✅ Complete |
| Parameter Sweep Framework (Phase C) — `tidal sweep` + `tidal analyze`, TOML config, adaptive/LHS/Sobol sampling, sensitivity analysis, convergence mode, 6 plot types, parallel execution | ✅ Complete |
| `theory.toml` configuration with `[[derived_fields]]`                                       | ✅ Complete |
| Scalar-vector coupling stress test (mixed-rank cross-field)                                 | ✅ Complete |
| Massive 3-form example (rank-3 antisymmetric tensor)                                        | ✅ Complete |
| Critical Review Pass 1 & 2 (fail-fast, CFL, operator plugin API)                            | ✅ Complete |
| Auto-computed mass/coupling matrices (Phase 12)                                             | ✅ Complete |
| Project rename to TIDAL (package, CLI, imports, docs)                                       | ✅ Complete |
| TIDAL logo (SVG) in README and Sphinx docs                                                  | ✅ Complete |
| Root Makefile with 13 convenience targets                                                   | ✅ Complete |
| `docs/COMMUNITY.md` — community guidelines and support channels                             | ✅ Complete |
| Sphinx extensions: autosummary, doctest                                                     | ✅ Complete |
| README Community & Support section                                                          | ✅ Complete |
| Performance optimizations (6 rounds + Phases 1-3): hot-path elimination, ghost-cell padding, fused operators, higher-order FD stencils (`--fd-order 4\|6`), Yoshida 4th-order leapfrog (`--leapfrog-order 4`), FFT spectral operators (`--spectral`) | ✅ Complete |
| Analytical Jacobian — precomputed dF/dy + dF/dyp for time-independent systems (5.3x IDA speedup) | ✅ Complete |
| Multi-field perturbation pipeline — `[[linearization.matter_perturbations]]` TOML for xPert DefTensorPerturbation | ✅ Complete |
| Curved-metric pipeline — non-constant metrics (spherical, cylindrical) in derive → JSON workflow | ✅ Complete |
| TT gauge constraint elimination — transverse constraints replace constrained EOM in-place (Wolfram-side) | ✅ Complete |
| Simulation progress bar — tqdm-based `SimulationProgress` class, auto-suppressed in sweep inner runs | ✅ Complete |
| Simulation resume — `--resume DIR [--snapshot N] [--t-additional T]` checkpoint loading | ✅ Complete |
| Gertsenshtein effect (Phase D) — graviton-photon conversion validated: uniform B₀ (sin²(κB₀t/2)), localized Gaussian (Boccaletti formula). **Torsion-independence proven**: minimal PGT coupling cannot amplify Gertsenshtein for any vacuum GW (#199, #200) | ✅ Complete |
| Torsion-Gertsenshtein investigation — combined PGT+EM theory (23 components), polarisation block-diagonal structure, field-filtered energy measurement for Ostrogradsky theories | ✅ Complete |

---

## Priority Legend

| Symbol | Priority     | Criteria                                    |
| ------ | ------------ | ------------------------------------------- |
| 🔴     | **CRITICAL** | Security, correctness, or major user impact |
| 🟠     | **HIGH**     | Important features or significant gaps      |
| 🟡     | **MEDIUM**   | Improvements that enhance robustness        |
| 🟢     | **LOW**      | Nice-to-have enhancements                   |

---

## Labels for GitHub Issues

Issues should be tagged with appropriate labels:

### Priority Labels

- `priority: critical` - 🔴 Security/correctness issues
- `priority: high` - 🟠 Important improvements
- `priority: medium` - 🟡 Robustness enhancements
- `priority: low` - 🟢 Nice-to-have features

### Type Labels

- `bug` - Something isn't working
- `enhancement` - New feature or request
- `documentation` - Improvements or additions to documentation
- `testing` - Test additions or improvements
- `ci-cd` - Continuous integration/deployment
- `performance` - Performance optimization
- `security` - Security vulnerability
- `refactoring` - Code quality improvement

### Component Labels

- `wolfram` - Mathematica/xAct pipeline
- `python` - Python simulation code
- `examples` - Example additions/improvements
- `animation` - Visualization features
- `validation` - Input/output validation

---

## Current Focus

**As of v0.6.0:**

- ✅ Phase 13+ completed: All core pipeline features implemented
- ✅ Solver migration: SUNDIALS IDA/CVODE + leapfrog + scipy replaces py-pde
- ✅ Adaptive time-stepping: tolerance-controlled integration (Phase F)
- ✅ Background fields: position-dependent coefficients (Phase A)
- ✅ Gauge fixing: optional per-field Lorenz/de Donder/Coulomb/temporal/axial (Phase B)
- ✅ Constraint pre-solve: FFT/sparse three-tier solver (Phase J)
- ✅ Parameter sweep framework (Phase C): `tidal sweep` + `tidal analyze`, 12 measurement types, adaptive/LHS/Sobol sampling, Sobol/Morris sensitivity, convergence mode, 6 plot types
- ✅ CLI (`tidal` command) implemented: 9 subcommands, zero new dependencies
- ✅ Measurement module: 12 types (energy, conversion, mixing, spectrum, dispersion, conservation, effective_mass, asymptotic, peak_conversion, velocity, resonance, summary)
- ✅ 20 working examples spanning 1+1D through 3+1D
- ✅ 1,700 Python tests + ~115 Wolfram tests passing, 0 ruff violations, 0 pyright errors
- ✅ 21 of 25 original issues resolved (84%)
- ✅ **Phase D (Gertsenshtein):** Complete — all physical regimes validated, torsion-independence proven for minimal PGT
- 🔄 **Primary remaining focus:** Non-minimal torsion-EM coupling investigation (T·F terms), ghost-free parameter conditions
- 🔄 **Secondary:** Phase 2 (Wolfram CI), Phase G (Absorbing Boundaries)

---

## Future Vision (Beyond Current Roadmap)

### Long-Term Goals

1. **Gertsenshtein Effect (Phase D)**
   - Coupled EM-gravity simulation from a single `theory.toml`
   - Validation against analytical thin-magnet formula (Domcke & Garcia-Cely 2023)
   - Automated analytic benchmark tests

2. **Poincaré Gauge Theory (Torsion)** — ✅ COMPLETE (v0.18.0)
   - ✅ PGT Lagrangian support (T² invariants + R̃²) — COMPLETE
   - ✅ Component-level E-L derivation (5s vs 77min) — COMPLETE (now default for ALL theories)
   - ✅ Ostrogradsky reduction (4th→2nd order) — COMPLETE (automatic on JSON load)
   - ✅ Simulation via generalized mass-matrix modal solver (v0.16.0, #165 resolved)
   - ✅ QZ generalized eigenvalue for velocity coupling singularity (v0.16.1, #166 resolved)
   - ✅ 4D R̃² derivation pipeline with canonical Hamiltonian (v0.17.0, #170)
   - ✅ Full general quadratic Lagrangian (α₁I₁ + α₂I₂ + α₃I₃ + b₅R̃²) — COMPLETE (v0.18.0)
   - ✅ Machine-precision energy conservation (max|dE/E| = 0) for all parameter values tested
   - 🔄 Dispersion relation validation against Nikiforova/Barker predictions (#169)
   - 🔄 Ghost propagator analysis for Ostrogradsky stability (#164)
   - Parameter window scanning for viable mode configurations

3. **Non-Abelian Gauge Theories**
   - Yang-Mills equations (linearized sector)
   - SU(2) and SU(3) gauge groups

4. **Spectral Spatial Discretisation (Phase E)** — ✅ FFT COMPLETE
   - ~~FFT-based operators for exponential convergence on periodic domains~~ **Done** (`--spectral`, auto-enabled for all-periodic BCs)
   - Chebyshev for non-periodic directions (following Dedalus architecture) — remaining

5. **Performance Scaling**
   - GPU acceleration for large grids
   - Adaptive mesh refinement

---

## Contributing to the Roadmap

See tracked issues with the [`roadmap`](https://github.com/WilliamRoyce/torsion-gertsenshtein/labels/roadmap) label.

**To propose new features:**

1. Check existing issues to avoid duplicates
2. Open a feature request with the `enhancement` label
3. Describe use case, motivation, and potential implementation
4. Maintainers will triage and add to roadmap if appropriate

**To claim an issue:**

1. Comment on the issue expressing interest
2. Wait for maintainer assignment/approval
3. Follow [CONTRIBUTING.md](../CONTRIBUTING.md) guidelines

---

## Versioning Strategy

This project follows [Semantic Versioning](https://semver.org/):

- **MAJOR** (1.0.0): Rank-3+ tensor support, automatic gauge fixing
- **MINOR** (0.X.0): New features (3+1D examples, JSON schema extensions)
- **PATCH** (0.0.X): Bug fixes, documentation improvements

**Current Version:** 0.41.4
**Previous Milestones:** 0.3.0 delivered Phase 3 + CLI + rename to TIDAL; 0.4.0 delivered solver migration + gauge fixing + background fields + adaptive timestepping + constraint pre-solve; 0.5.0 delivered parameter sweep framework (Phase C) with 12 measurements, sensitivity analysis, and advanced visualization
**Next Major Release (1.0.0):** Phase D (Gertsenshtein example) + Wolfram CI

---

## Questions or Feedback?

- **Open an issue:** For roadmap suggestions
- **Start a discussion:** For feature brainstorming
- **Join development:** See [CONTRIBUTING.md](../CONTRIBUTING.md)

---

_This roadmap is based on comprehensive codebase analysis identifying 25 improvement areas. Priorities and timelines may adjust based on community feedback and contributions._
