# Torsion-Gertsenshtein Pipeline Roadmap

This document outlines the planned improvements and features for the torsion-gertsenshtein symbolic physics pipeline project.

**Last Updated:** February 2026
**Project Status:** Phase 10b Complete (323 tests passing)

## Overview

The project is in a mature state with a robust symbolic pipeline from Lagrangian to PDE simulation. This roadmap tracks refinements toward production-grade reliability and extensibility.

**Total Identified Improvements:** 25 issues
**Estimated Effort:** 100-150 hours (6-8 sprint cycles)

---

## Implementation Phases

### Phase 1: Critical Fixes (Next Sprint)
**Goal:** Address security concerns and high-impact bugs
**Estimated Duration:** 1-2 weeks

| Issue | Priority | Type | Effort |
|-------|----------|------|--------|
| [#TBD] Replace Assertions with Explicit Error Handling | 🔴 Critical | Bug | Medium |
| [#TBD] Strengthen Mathematica Expression Evaluation Security | 🔴 Critical | Security | Medium |
| [#TBD] Validate Grid Dimensions During PDE Construction | 🟠 High | Bug | Low |

**Key Deliverables:**
- All `assert` statements replaced with explicit error checks
- NaN/Inf validation for coefficient evaluation
- Grid dimension validation in `PDEFromSpec.__init__()`

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

### Phase 3: Features & Documentation (1 Month)
**Goal:** Validate 3+1D support and improve user documentation
**Estimated Duration:** 3-4 weeks

| Issue | Priority | Type | Effort |
|-------|----------|------|--------|
| [#TBD] Add 3+1D Spacetime Examples | 🔴 Critical | Documentation | High |
| [#TBD] Document JSON Schema with Detailed Guide | 🟠 High | Documentation | Medium |
| [#TBD] Improve Coefficient Resolution Performance | 🟠 High | Performance | Medium |
| [#TBD] Create Architecture Diagrams | 🟡 Medium | Documentation | Low |

**Key Deliverables:**
- 3+1D Klein-Gordon example with visualization
- `docs/JSON_SCHEMA_GUIDE.md` with complete field reference
- Coefficient caching layer (10-20% speedup for `sphere_kg/`)
- Architecture diagrams (data flow, module dependencies, state structure)

---

### Phase 4: Advanced Features (2-3 Months)
**Goal:** Implement major feature requests
**Estimated Duration:** 8-12 weeks

| Issue | Priority | Type | Effort |
|-------|----------|------|--------|
| [#TBD] Support Rank-3+ Tensor Decomposition | 🔴 Critical | Feature | Very High |
| [#TBD] Implement Automatic Gauge Fixing | 🔴 Critical | Feature | High |
| [#TBD] Handle Mixed Time-Space Cross-Derivatives Properly | 🟡 Medium | Bug | Medium |
| [#TBD] Expand _mathematica_to_python Function Set | 🟡 Medium | Feature | Low |
| [#TBD] Add Non-Cartesian Coordinate System Examples | 🟡 Medium | Documentation | Medium |
| [#TBD] Add Convergence and Stability Stress Tests | 🟡 Medium | Testing | Medium |
| [#TBD] Add Full Pipeline Validation to CI | 🟡 Medium | CI/CD | Medium |

**Key Deliverables:**
- Rank-3+ tensor support in `ComponentDecompose.wl`
- Lorenz and Coulomb gauge automation in `Linearize.wl`
- Nested operator support for mixed derivatives
- Extended Mathematica function set (inverse trig, hyperbolic, special functions)
- Spherical/cylindrical coordinate examples
- Convergence test suite with O(N^-2) validation
- Weekly JSON regeneration workflow

---

### Phase 5: Polish & Optimization (Ongoing)
**Goal:** Code quality improvements and nice-to-have features
**Estimated Duration:** Ongoing

| Issue | Priority | Type | Effort |
|-------|----------|------|--------|
| [#TBD] Refactor ContainsTimeDerivative and IsMixedTimeSpaceDerivative | 🟢 Low | Refactoring | Low |
| [#TBD] Add Parameter Sweep Examples | 🟢 Low | Documentation | Low |
| [#TBD] Add Python 3.12+ Testing to CI Matrix | 🟢 Low | CI/CD | Very Low |
| [#TBD] Add PGF Export Module Tests | 🟢 Low | Testing | Low |
| [#TBD] Cache Coordinate Symbols and Chart Dimension | 🟢 Low | Performance | Low |
| [#TBD] Add Debugging and Performance Tuning Guides | 🟢 Low | Documentation | Medium |
| [#TBD] Support Elliptic PDE Solving (Constraint Equations) | 🟢 Low | Feature | High |

**Key Deliverables:**
- Unified `ExtractDerivativeProfile` function
- Parameter sweep example with joblib parallelization
- Python 3.12 CI matrix (allow_failure: true)
- PGF export test suite
- Memoized coordinate symbol lookups (10-20% speedup)
- `docs/DEBUGGING.md` and `docs/PERFORMANCE_TUNING.md`
- Elliptic PDE solver integration

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
- ✅ Phase 10b completed: Unified coefficient resolution, 323 tests passing
- 🔄 Phase 1 planning: Critical bug fixes and security improvements
- 📋 Issue creation in progress: All 25 issues identified and ready to track

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
