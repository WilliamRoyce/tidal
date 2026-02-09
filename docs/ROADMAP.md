# TIDAL Pipeline Roadmap

This document outlines the planned improvements and features for the torsion-gertsenshtein symbolic physics pipeline project.

**Last Updated:** February 2026
**Project Status:** Phase 13+ Complete, CLI done (743 Python tests + ~108 Wolfram tests)

## Overview

The project is in a mature state with a robust symbolic pipeline from Lagrangian to PDE simulation. This roadmap tracks refinements toward production-grade reliability and extensibility.

**Total Identified Improvements:** 25 issues
**Estimated Effort:** 100-150 hours (6-8 sprint cycles)

---

## Implementation Phases

### Phase 1: Critical Fixes ✅ COMPLETE
**Goal:** Address security concerns and high-impact bugs

| Issue | Priority | Type | Status |
|-------|----------|------|--------|
| [#67] Replace Assertions with Explicit Error Handling | 🔴 Critical | Bug | ✅ Done |
| [#68] Strengthen Mathematica Expression Evaluation Security | 🔴 Critical | Security | ✅ Done |
| [#75] Validate Grid Dimensions During PDE Construction | 🟠 High | Bug | ✅ Done |

**Delivered:**
- All `assert isinstance()` replaced with explicit `TypeError` raises (Issue #67)
- `_validate_eval_result` checks for NaN/Inf/complex (Issue #68)
- `_validate_operator_dimensions()` at construction time (Issue #75)

---

### Phase 2: Testing & CI Infrastructure (1-2 Weeks)
**Goal:** Improve test coverage and CI reliability
**Estimated Duration:** 2-3 weeks

| Issue | Priority | Type | Effort |
|-------|----------|------|--------|
| [#TBD] Add Animation Module Test Coverage | 🟠 High | Testing | Medium |
| [#TBD] Add Code Coverage Reporting to CI | 🟠 High | CI/CD | Low |
| [#TBD] Add Wolfram Tests to GitHub Actions CI | 🔴 Critical | CI/CD | High |
| [#TBD] Add Tests for Observers, Profiling, Runners Modules | 🟡 Medium | Testing | Medium |

**Key Deliverables:**
- `test_animation_builder.py` with 10-15 test cases
- Codecov integration in GitHub Actions
- Optional Wolfram test workflow (label-triggered or weekly)
- ~20-30 tests for observers, profiling, runners modules
- Coverage increase from ~70% to 80%+

---

### Phase 3: Features & Documentation ✅ COMPLETE
**Goal:** Validate 3+1D support and improve user documentation

| Issue | Priority | Type | Status |
|-------|----------|------|--------|
| [#71] Add 3+1D Spacetime Examples | 🔴 Critical | Documentation | ✅ Done |
| [#TBD] Document JSON Schema with Detailed Guide | 🟠 High | Documentation | ✅ Done |
| [#TBD] Improve Coefficient Resolution Performance | 🟠 High | Performance | ✅ Done |
| [#TBD] Create Architecture Diagrams | 🟡 Medium | Documentation | Deferred |

**Delivered:**
- 3+1D Klein-Gordon example (Issue #71), plus spherical_kg, cylindrical_kg, gravitational_waves, massive_3form
- `docs/JSON_SCHEMA_GUIDE.md` with complete field reference
- Unified `_resolve_coefficient_at_point` evaluator (Phase 10b)

---

### Phase 4: Advanced Features ✅ MOSTLY COMPLETE
**Goal:** Implement major feature requests

| Issue | Priority | Type | Status |
|-------|----------|------|--------|
| [#70] Support Rank-3+ Tensor Decomposition | 🔴 Critical | Feature | ✅ Done (Phase 13) |
| [#TBD] Implement Automatic Gauge Fixing | 🔴 Critical | Feature | Remaining |
| [#79] Handle Mixed Time-Space Cross-Derivatives | 🟡 Medium | Bug | ✅ Done |
| [#TBD] Expand _mathematica_to_python Function Set | 🟡 Medium | Feature | ✅ Done |
| [#TBD] Add Non-Cartesian Coordinate System Examples | 🟡 Medium | Documentation | ✅ Done |
| [#TBD] Add Convergence and Stability Stress Tests | 🟡 Medium | Testing | ✅ Done |
| [#TBD] Add Full Pipeline Validation to CI | 🟡 Medium | CI/CD | Remaining |

**Delivered:**
- Rank-3+ tensor support: `ReplaceHigherRankFieldComponents`, massive_3form example (Phase 13)
- Mixed time-space derivatives: `ClassifySpatialProfile`, `ExtractSpatialOperatorFromMixed` (Issue #79)
- Curvilinear examples: polar_kg, spherical_kg, cylindrical_kg with Christoffel auto-detection
- CFL `check_stability()`, physics validation (energy conservation + analytical solutions)

**Remaining:**
- Automatic gauge fixing (Lorenz/Coulomb in Wolfram layer)
- Wolfram tests in GitHub Actions CI

---

### Phase 5: Polish & Optimization (Partially Complete)
**Goal:** Code quality improvements and nice-to-have features

| Issue | Priority | Type | Status |
|-------|----------|------|--------|
| [#85] Refactor derivative classification | 🟢 Low | Refactoring | ✅ Done |
| [#TBD] Add Parameter Sweep Examples | 🟢 Low | Documentation | Remaining |
| [#TBD] Add Python 3.12+ Testing to CI Matrix | 🟢 Low | CI/CD | Remaining |
| [#TBD] Support Elliptic PDE Solving (Constraint Equations) | 🟢 Low | Feature | ✅ Done |

**Delivered:**
- Unified `ExtractDerivativeProfile` function (Issue #85, ~108 lines saved)
- Elliptic PDE solving: constraint equations with `time_derivative_order=0`, `--mode constraint` CLI flag

**Remaining:**
- Parameter sweep examples
- Python 3.12 CI matrix

---

### Additional Completed Work (not in original roadmap)

| Feature | Status |
|---------|--------|
| CLI (`tidal` command) — 5 subcommands, 147 tests | ✅ Complete |
| `theory.toml` configuration with `[[derived_fields]]` | ✅ Complete |
| Scalar-vector coupling stress test (mixed-rank cross-field) | ✅ Complete |
| Massive 3-form example (rank-3 antisymmetric tensor) | ✅ Complete |
| Critical Review Pass 1 & 2 (fail-fast, CFL, operator plugin API) | ✅ Complete |
| Auto-computed mass/coupling matrices (Phase 12) | ✅ Complete |

---

## Priority Legend

| Symbol | Priority | Criteria |
|--------|----------|----------|
| 🔴 | **CRITICAL** | Security, correctness, or major user impact |
| 🟠 | **HIGH** | Important features or significant gaps |
| 🟡 | **MEDIUM** | Improvements that enhance robustness |
| 🟢 | **LOW** | Nice-to-have enhancements |

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

**As of February 2026:**
- ✅ Phase 13+ completed: All core pipeline features implemented
- ✅ CLI (`tidal` command) implemented: 5 subcommands, 147 tests, zero new dependencies
- ✅ 18 working examples spanning 1+1D through 3+1D
- ✅ 743 Python tests + ~108 Wolfram tests passing, 0 ruff violations, 0 pyright errors
- 🔄 Remaining: Automatic gauge fixing, Wolfram CI, parameter sweep examples

---

## Future Vision (Beyond Current Roadmap)

### Long-Term Goals
1. **Full General Relativity Support**
   - Riemann curvature tensor (rank 4)
   - Einstein field equations
   - Black hole spacetime examples

2. **Non-Abelian Gauge Theories**
   - Yang-Mills equations
   - SU(2) and SU(3) gauge groups
   - Automatic gauge fixing for non-Abelian theories

3. **Quantum Field Theory Features**
   - Renormalization group equations
   - Loop corrections automation
   - Feynman diagram generation

4. **Interactive Visualization**
   - Real-time parameter adjustment
   - Web-based simulation viewer
   - Interactive 3D visualization

5. **Performance Scaling**
   - GPU acceleration for large grids
   - Distributed computing support
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

**Current Version:** 0.2.4
**Next Minor Release (0.3.0):** Phase 3 completion (3+1D examples, JSON schema guide)
**Next Major Release (1.0.0):** Phase 4 completion (rank-3+ tensors, gauge automation)

---

## Questions or Feedback?

- **Open an issue:** For roadmap suggestions
- **Start a discussion:** For feature brainstorming
- **Join development:** See [CONTRIBUTING.md](../CONTRIBUTING.md)

---

*This roadmap is based on comprehensive codebase analysis identifying 25 improvement areas. Priorities and timelines may adjust based on community feedback and contributions.*
