<p align="center">
  <img src="docs/TIDAL_Logo_TikZ_Figure.svg" alt="TIDAL: Tensor Integration and Derivation for Any Lagrangian" width="500">
</p>
<p align="center"><em>Tensor Integration and Derivation for Any Lagrangian</em></p>

<p align="center">

[![CI Tests](https://github.com/WilliamRoyce/torsion-gertsenshtein/workflows/test/badge.svg)](https://github.com/WilliamRoyce/torsion-gertsenshtein/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![Code style: ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![Type checked: pyright](https://img.shields.io/badge/type%20checked-pyright-informational)](https://github.com/microsoft/pyright)
[![Coverage](https://img.shields.io/badge/coverage-see%20CI-brightgreen)](https://github.com/WilliamRoyce/torsion-gertsenshtein/actions/workflows/test.yml)

</p>

View the `tidal` package documentation [here](https://williamroyce.github.io/torsion-gertsenshtein/).

A research codebase for exploring **electromagnetic ↔ gravitational wave conversion** ([Gertsenshtein effect](https://arxiv.org/abs/2301.02072); Gertsenshtein 1962, Domcke & Garcia-Cely 2023) and potential **amplification mechanisms** in gravity theories with **torsion** (Poincaré gauge theory; parity-even quadratic invariants). The repository includes:

- **A native PDE solver framework** (SUNDIALS IDA/CVODE + leapfrog + scipy, with numpy spatial operators) for time-domain simulations with **1,343 Python tests + ~115 Wolfram tests**.
- A symbolic pipeline (Mathematica + xAct) for **deriving linearized field equations** and exporting them to Python-friendly JSON specifications.
- Documentation and experiments for **mixing mechanisms** and **hyperbolicity/causality checks** relevant to the effect.

> TL;DR: start with Klein–Gordon toy systems in 1+1D → grow to coupled EM/metric/torsion perturbations → test conversion and stability in controlled scenarios.

---

## Community & Support

- **Questions & Ideas**: [GitHub Discussions](https://github.com/WilliamRoyce/torsion-gertsenshtein/discussions) — Ask questions, share use cases, discuss physics
- **Bug Reports**: [Issue Tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues) — Report bugs with the `bug` label
- **Feature Requests**: [Issue Tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues) — Propose features with the `enhancement` label
- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and workflow

For more, visit the [Documentation](https://williamroyce.github.io/torsion-gertsenshtein/).

---

## Current Status (usable today)

- **Lagrangian-to-PDE pipeline (`tidal.symbolic`, `tidal.wolfram`)**: complete symbolic-to-numerical pipeline for deriving field equations from Lagrangian densities. Uses Mathematica/xAct for symbolic derivation (Euler-Lagrange equations, linearization via xPert, component decomposition) → JSON export → native Python solvers for PDE construction and time integration. **Zero hardcoded physics** in the numerical layer — all equations derived symbolically. Includes **25 working examples** spanning 1+1D through 3+1D: scalars, vectors, tensors (rank 3+), coupled multi-field systems, curvilinear coordinates, curved spacetimes, and background-field scattering. See [examples/README.md](examples/README.md) for complete documentation.
- **Solver architecture (`tidal.solver`)**: four time-integration backends — **IDA** (SUNDIALS DAE solver for systems with algebraic constraints; Hindmarsh et al. 2005), **CVODE** (SUNDIALS BDF adaptive ODE with tolerance control), **leapfrog** (Störmer-Verlet symplectic integrator; Hairer et al. 2006), and **scipy** (`solve_ivp` with DOP853/Radau/BDF). Automatic solver selection based on equation structure: systems with constraints route to IDA, pure wave equations to CVODE or leapfrog. Pure numpy spatial operators (`tidal/solver/operators.py`) with 2nd-order finite-difference stencils. Three-tier constraint pre-solve (FFT → sparse matrix → automatic selection) with gauge regularization for singular Poisson problems.
- **CLI (`tidal` command)**: unified command-line interface with 9 subcommands — `tidal derive` (Lagrangian → JSON via TOML config), `tidal simulate` (JSON → PDE simulation with plotting), `tidal measure` (post-hoc measurement extraction from snapshot directories), `tidal inspect` (equation system info), `tidal list` (discover available specs), `tidal validate` (JSON spec validation), `tidal plot` (standalone plotting from simulation output), `tidal sweep` (automated parameter sweeps with adaptive sampling, convergence analysis, and sensitivity analysis), `tidal analyze` (post-hoc Sobol/Morris sensitivity analysis of sweep results). Supports `theory.toml` configs with `[[derived_fields]]`, `[[background_fields]]`, and optional `[[gauge]]` sections. Zero new dependencies (stdlib argparse + tomllib).
- **Measurement module (`tidal.measurement`)**: 13 post-hoc analysis types — summary, energy, conversion P(t), mixing length, spectrum, spectral conversion P(k,t), dispersion omega(k), conservation diagnostics, effective mass, asymptotic scattering, peak conversion, group/phase velocity, and resonance analysis. Disk-backed snapshot storage for long simulations via `SnapshotWriter`.
- **Dev environment**: container-first, [`uv`] for Python (3.11 pinned), Wolfram Engine 14.3 with xAct tensor framework, optional ffmpeg; Sphinx docs skeleton; type-checked codebase with pytest test suite.
- **Professional development infrastructure**: 1,343 Python tests + ~115 Wolfram tests, 5 utility scripts for streamlined workflows (`run_wolfram_tests.sh`, `run_examples.sh`, `full_test.sh`, `validate_pipeline.sh`, `lint_wolfram.sh`), comprehensive documentation with module headers and usage strings, robust kernel caching handling for reliable test execution. 0 ruff violations, 0 pyright errors (strict mode).

This README describes the current capabilities, how to run the examples, and planned improvements.

---

# Project Scope and Milestones

This section clarifies the project scope and milestones.

## Symbolic Computing Infrastructure

The project now includes a complete symbolic tensor algebra pipeline for deriving linearized field equations:

- **Wolfram Engine 14.3**: Headless installation with free license activation
- **xAct Tensor Algebra Suite**: State-of-the-art packages for General Relativity computations
  - **xCore**: Generic programming tools and core functionality
  - **xPerm**: Large group permutation manipulation (GLIBC-compatible binary)
  - **xTensor**: Abstract tensor computations (flagship package)
  - **xCoba**: Coordinate-based tensor computations for component calculations
- **Automated Setup**: Container-friendly installation scripts with verification
- **Compatibility Fixes**: Handles GLIBC version mismatches by recompiling xPerm binary

### Usage Example

```wolfram
Needs["xAct`xCoba`"];
DefManifold[M, 4, IndexRange[a, z]];
DefChart[cart, M, {0, 1, 2, 3}, {t[], x[], y[], z[]}];
DefMetric[-1, g[-a, -b], CD, {";", "∇"}, PrintAs -> "g"];
```

See [`scripts/README.md`](scripts/README.md) for complete setup instructions and [`scripts/verify-wolfram-setup.sh`](scripts/verify-wolfram-setup.sh) for verification.

**Development Tools:**

- Comprehensive test suite with ~115 Wolfram unit tests covering all pipeline modules
- 5 utility scripts for workflow automation (test execution, example derivations, pipeline validation)
- Robust kernel caching support ensuring reliable test execution across multiple runs

## Objectives

- Baseline re-derivation of the standard Gertsenshtein effect (Einstein–Maxwell) and its tiny conversion amplitude.
- Extend the gravitational sector to parity-even quadratic PGT with torsion; identify propagating modes and viable parameter windows.
- Linearized PDE system in a flat metric background with constant external magnetic field and (if allowed) homogeneous torsion background. Extract mixing terms.
- Well-posedness: characteristic analysis, hyperbolicity, and causality (characteristic speeds).
- Numerical experiments: 1+1D toy models mapping EM/GR/torsion modes to coupled scalars; verify conversion scaling and stability; then scale up in fidelity.

## Recent Improvements

- **Solver migration to SUNDIALS (February 2026)**: Replaced py-pde with a native solver architecture: SUNDIALS IDA for DAE systems (algebraic constraints), CVODE for adaptive BDF time-stepping with tolerance control, Störmer-Verlet leapfrog for symplectic integration, and scipy `solve_ivp` for general-purpose adaptive ODE. Pure numpy spatial operators. Automatic solver selection based on equation structure. Three-tier constraint pre-solve (FFT, sparse matrix, automatic). See [docs/solver_migration.md](docs/solver_migration.md) and [docs/adaptive_timestepping.md](docs/adaptive_timestepping.md).
- **Parameter Sweep Framework (March 2026)**: Complete `tidal sweep` and `tidal analyze` commands with 13 measurement types, TOML sweep configuration, adaptive + Latin Hypercube + Sobol sampling, Sobol/Morris sensitivity analysis, run status tracking, and advanced visualization (parallel coordinates, tornado, scatter plots). See [docs/next-features.md](docs/next-features.md).
- **Phase 4-13+ Pipeline Evolution (February 2026)**: All critical implementation complete. **Phase 12**: Auto-computed mass/coupling matrices with symbolic preservation. **Phase 13**: Rank 3+ tensor support. **CLI**: Full `tidal` command with 9 subcommands. **Background fields**: `[[background_fields]]` TOML for non-dynamical tensors with position-dependent coefficients and 4-level caching. **Gauge fixing**: Optional per-field `[[gauge]]` TOML (Lorenz, de Donder, Coulomb, temporal, axial). **1,343 Python tests + ~115 Wolfram tests passing**. See [CHANGELOG.md](CHANGELOG.md) for complete history.
- **Lagrangian-to-PDE pipeline (February 2026)**: Complete symbolic derivation pipeline: Mathematica/xAct → JSON → native Python solvers. Canonical momentum pipeline with symbolic K^{-1} inversion for non-diagonal kinetic matrices. Lagrangian-first linearization via xPert (`Perturbation[L, 2]`). See [CHANGELOG.md](CHANGELOG.md) for details.

## Future Development

- **Gertsenshtein example (Phase D)**: Coupled EM-gravity simulation — the project's primary research target, now unblocked by Phases A (background fields), B (gauge fixing), and C (parameter sweeps).
- **Continuous Integration**: GitHub Actions workflow for automated Wolfram test execution on pull requests.
- **Spectral spatial discretization (Phase E)**: FFT-based operators for exponential convergence on periodic domains (following Dedalus architecture).
- **Absorbing boundaries (Phase G)**: Sponge layers and PML (Bérenger 1994) for finite-magnet interaction regions.
- **Extended physics examples**: Coupled EM/torsion systems for Poincaré gauge theory research.

---

# Development Environment

## Quickstart

This project uses **uv** for Python version/venv/dependencies and is container-first.

```bash
# Ensure Python 3.11 is used (numba/llvmlite friendly)
uv python pin 3.11

# Install runtime + (optional) dev dependencies from pyproject.toml
uv sync --all-extras

# Smoke test: can we import the package?
uv run python -c "import tidal; print('OK')"
```

## Dev Container (VS Code / Codespaces)

This repo includes a Debian-based **VS Code Dev Container**. It ensures a consistent toolchain and avoids host-machine drift.

- Open the folder in VS Code
- Command Palette → Dev Containers: Reopen in Container
- Once inside the container:

```bash
uv python pin 3.11
uv sync --all-extras
```

Common CLI tools pre-installed in the container: `git`, `node`, `npm`, `eslint`, `apt`, `dpkg`, `curl`, `wget`, `ssh`, `rsync`, `gpg`, `tree`, `find`, `grep`, `zip`, `tar`, `gzip`, etc.

## Running the Examples

### Lagrangian-to-PDE Pipeline Examples

The repository includes a complete symbolic-to-numerical pipeline for deriving field equations from Lagrangians and simulating them numerically. **25 examples** cover scalars, vectors, rank-3 tensors, coupled multi-field systems, curvilinear coordinates, curved spacetimes, and background-field scattering.

```bash
# Each example has a run.sh showing the full derive → inspect → simulate workflow:
cd examples/scalar_field && bash run.sh

# Or use the CLI directly:
tidal derive examples/scalar_field/theory.toml        # derive equations from Lagrangian
tidal inspect examples/data/klein_gordon_1d.json       # inspect equation structure
tidal simulate examples/data/klein_gordon_1d.json \    # simulate
  --param m2=1.0 --ic gaussian --t-end 20
tidal list                                             # discover all available JSON specs
tidal validate examples/data/klein_gordon_1d.json      # validate JSON spec structure
```

**CLI Subcommands:**

| Command                     | Description                                                                    |
| --------------------------- | ------------------------------------------------------------------------------ |
| `tidal derive theory.toml`  | Generate .wls from TOML, run wolframscript to produce JSON                     |
| `tidal simulate spec.json`  | Full simulation with plotting (supports `--param`, `--ic`, `--bc`, `--scheme`) |
| `tidal measure result_dir/` | Extract physics measurements (energy, conversion, mixing length, spectra)      |
| `tidal inspect spec.json`   | Display equation system info (fields, operators, parameters)                   |
| `tidal list`                | Discover all available JSON specs in `examples/data/`                          |
| `tidal validate spec.json`  | Validate JSON equation specification structure                                 |
| `tidal plot result_dir/`    | Standalone plotting from simulation output directories                         |
| `tidal sweep spec.json`     | Parameter sweeps, convergence studies, and adaptive sampling                   |
| `tidal analyze sweep_dir/`  | Post-hoc sensitivity analysis (Sobol/Morris) of sweep results                  |

**TOML Configuration** (`theory.toml`):

- Define spacetime dimension, metric, fields, constants, and Lagrangian expression
- `[[derived_fields]]` section for intermediate tensors (e.g., field strength `F_ab = CD[-a][A[-b]] - CD[-b][A[-a]]`)
- Runtime parameters with default values in `[parameters]` section

**Pipeline Examples:**

| Example                   | Dim  | Key Features                                                                         |
| ------------------------- | ---- | ------------------------------------------------------------------------------------ |
| `scalar_field/`           | 1+1D | Klein-Gordon, mass term, dispersion                                                  |
| `electromagnetic/`        | 1+1D | Maxwell, Lorenz gauge, massless waves                                                |
| `proca/`                  | 1+1D | Massive vector field (Proca mass)                                                    |
| `coupled_scalars/`        | 1+1D | Cross-field coupling, mass matrix, energy transfer                                   |
| `chern_simons/`           | 2+1D | Epsilon tensor, topological mass, A_0 constraint                                     |
| `elasticity/`             | 2+1D | Anisotropic laplacian, cross_derivative_xy                                           |
| `curved_spacetime/`       | 2+1D | De Sitter, Hubble friction, time-dependent coefficients                              |
| `sphere_kg/`              | 2+1D | KG on S², position-dependent coefficients                                            |
| `polar_kg/`               | 2+1D | Polar coordinates, Christoffel auto-detection                                        |
| `electrostatics/`         | 2+1D | Poisson equation, constraint solver                                                  |
| `scalar_vector_coupling/` | 2+1D | Mixed-rank cross-field (scalar+vector), 4 constants, CS+coupling                     |
| `scalar_field_3d/`        | 3+1D | Full 4D KG, 32^3 grid                                                                |
| `spherical_kg/`           | 3+1D | Spherical coordinates, trig coefficients                                             |
| `cylindrical_kg/`         | 3+1D | Cylindrical coordinates, mixed curved/flat                                           |
| `gravitational_waves/`    | 3+1D | xPert linearization, TT gauge, constraints                                           |
| `massive_3form/`          | 3+1D | Rank-3 antisymmetric tensor, symmetry reduction                                      |
| `massive_gravity/`        | 2+1D | Linearized massive gravity, Fierz-Pauli mass, xPert, coupled constraints             |
| `coupled_proca/`          | 2+1D | Two massive vectors, coupled Helmholtz constraints, periodic BCs                     |
| `coupled_scattering/`     | 2+1D | Position-dependent Gaussian coupling, background fields, wave scattering             |
| `scalar_potential_well/`  | 1+1D | Background potential well, `[[background_fields]]`, bound states                     |
| `cylindrical_kg_1d/`      | 1+1D | Cylindrical coordinates, plane-wave dimensional reduction                            |
| `gravitational_waves_1d/` | 1+1D | Linearized gravity, plane-wave 1D reduction                                          |
| `spherical_kg_1d/`        | 1+1D | Spherical coordinates, plane-wave dimensional reduction                              |
| `proca_background/`       | 2+1D | Lorentzian scalar background, two Proca vectors, constraint+BG integration           |
| `vector_background/`      | 2+1D | Tanh domain wall vector background, ComponentValue mechanism, sign-changing coupling |

See [examples/README.md](examples/README.md) for complete documentation and verification that the Python layer contains zero hardcoded physics.

## (Optional) Video Support

For MP4 via Matplotlib’s FFMpegWriter:

```bash
# inside the dev container
sudo apt-get update && sudo apt-get install -y ffmpeg
```

If `ffmpeg` is unavailable, the example falls back to a GIF via Pillow.

## Tests

The project includes a comprehensive test suite with **1,343 Python tests + ~115 Wolfram tests**.

### Python Tests (1,343 tests)

```bash
# Run all Python tests with pytest
uv run pytest -v

# Run a specific test module
uv run pytest tests/test_json_loader.py -v

# Run with coverage report (HTML)
uv run pytest --cov=tidal --cov-report=html
open htmlcov/index.html  # View detailed HTML report

# Run with coverage report (terminal)
uv run pytest --cov=tidal --cov-report=term-missing

# Run with coverage report (XML for CI)
uv run pytest --cov=tidal --cov-report=xml
```

### Wolfram Tests (~115 tests)

```bash
# Run all Wolfram unit tests
./scripts/run_wolfram_tests.sh

# Run individual test files
wolframscript -file tests/wolfram/test_euler_lagrange.wls
wolframscript -file tests/wolfram/test_common_utilities.wls
wolframscript -file tests/wolfram/test_export_json.wls
```

### Complete Test Suite

```bash
# Run both Python and Wolfram tests
./scripts/full_test.sh

# Validate end-to-end pipeline (Lagrangian → JSON → simulation)
./scripts/validate_pipeline.sh

# Check Wolfram module syntax
./scripts/lint_wolfram.sh
```

**Test Coverage:**

- Symbolic derivation (Euler-Lagrange, component decomposition, JSON export)
- PDE construction and operator identification
- Initial conditions and boundary conditions
- Multi-field coupling and energy transfer
- Parameter sweep and convergence analysis
- 13 measurement types (energy, conversion, spectrum, sensitivity analysis, etc.)
- Edge cases (empty grids, invalid bounds, division by zero)
- Path traversal protection and validation

See [`scripts/README.md`](scripts/README.md) for detailed documentation of utility scripts.

---

## Documentation

The repo builds Sphinx docs and deploys to GitHub Pages via Actions.

### Local Build

```bash
# auto-generate API docs
uv run sphinx-apidoc --force --module-first -o docs/source/ tidal/

# build HTML
( cd docs && make html )

# open locally (container users: use $BROWSER to open in host browser)
python -m http.server -d docs/build/html 8000
# then navigate to http://localhost:8000 or run:
# $BROWSER http://localhost:8000
```

On push to `main`, CI builds the docs and deploys to:

```bash
https://williamroyce.github.io/torsion-gertsenshtein/
```

Use `$BROWSER <url>` (from within the devcontainer) to open the project documentation link in the host browser.

---

## Logging and Profiling

The library uses Python's `logging` module (no print statements). To see info-level logs (solver progress, profiling summaries):

```python
import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
```

Run examples with `profile=True` in the `run` call to see timing breakdowns (initialization delay, solver overhead, etc.).

---

## Symbolic Computing Setup

For symbolic tensor algebra computations (deriving linearized field equations):

```bash
# 1. Download Wolfram Engine installer
#    Visit https://www.wolfram.com/engine/
#    Place installer in third_party/ directory

# 2. Install and activate Wolfram Engine
sudo ./scripts/install-wolfram-engine.sh
./scripts/activate-wolfram.sh

# 3. Install xAct/xCoba tensor algebra packages
./scripts/install-xact-xcoba.sh

# 4. Verify complete setup
./scripts/verify-wolfram-setup.sh
```

The verification script checks:

- Wolfram Engine installation and activation
- xAct package installation (xCore, xPerm, xTensor, xCoba)
- xPerm binary GLIBC compatibility
- Full smoke test with tensor operations

See [`scripts/README.md`](scripts/README.md) for detailed setup instructions.

---

## Troubleshooting

- **Import errors in VS Code** (e.g., numpy not found): ensure the interpreter is the repo's venv (`.venv/bin/python3`), then reload the window.
- **`llvmlite/numba` build failures**: stick to **Python 3.11** (`uv python pin 3.11`).
- **FileNotFoundError: ffmpeg** — install ffmpeg (see apt command above) or let the example produce a GIF.
- **Type-checker warnings about third-party stubs** — run examples anyway; code uses TYPE_CHECKING guards and runtime-safe casts where necessary.
- **Pages 404 or deploy errors**: ensure Settings → Pages → Source = GitHub Actions and Actions → Workflow permissions = Read/Write.
- **Animation has low frame count**: ensure `snapshot_interval` in the `run` call matches your desired temporal resolution (e.g., set to `dt` for every integrator step). Increase `fps` in `choose_writer_and_out` for smoother playback.
- **Logging messages not visible**: call `logging.basicConfig(level=logging.INFO, ...)` at the start of your script or in `main()`.
- **Wolfram Engine not activated**: run `./scripts/activate-wolfram.sh` and enter your Wolfram ID credentials (free account at wolfram.com).
- **xPerm GLIBC errors** (`GLIBC_2.38 not found`): run `./scripts/install-xact-xcoba.sh` to recompile the binary for your system.
- **xAct packages not loading**: ensure xAct is installed in `~/.WolframEngine/Applications/xAct/` — run verification script for diagnosis.

---

## Contributing

- Open an issue or submit a PR.
- **Test requirements**: All changes must maintain 100% test pass rate (1,343 Python + ~115 Wolfram tests). New features require corresponding unit tests in both Python and Wolfram layers where applicable.
- Run `./scripts/full_test.sh` before submitting PRs to verify all tests pass.
- Follow the project's type-checking and linting conventions (keyword-only booleans, explicit type annotations, no print in library code).

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

This project builds on:

- [SUNDIALS](https://computing.llnl.gov/projects/sundials) — IDA (DAE) and CVODE (BDF) solvers via [scikit-sundae](https://github.com/NREL/scikit-sundae) (Hindmarsh et al. 2005).
- The [xAct/xTensor ecosystem](http://www.xact.es/) — symbolic tensor algebra (Martín-García et al.) powering the Lagrangian-to-PDE derivation pipeline.
- [xPert](https://www.researchgate.net/publication/1740524) — metric perturbation theory (Brizuela et al. 2009) for linearization.
- [`uv`](https://github.com/astral-sh/uv) — fast Python environment management.
- Originally built on [py-pde](https://py-pde.readthedocs.io/) (Zwicker, JOSS 2020); finite-difference stencil conventions retained in TIDAL's native operators.

Design decisions are informed by [Dedalus](https://arxiv.org/abs/1905.10388) (Burns et al. 2020), [MEEP](https://meep.readthedocs.io/) (Oskooi et al. 2010), and [FEniCS](https://fenicsproject.org/) (Baratta et al. 2023). The core physics targets the Gertsenshtein effect (Gertsenshtein 1962; [Domcke & Garcia-Cely 2023](https://arxiv.org/abs/2301.02072)). See [`docs/references.md`](docs/references.md) for the full citation list.

[`uv`]: https://github.com/astral-sh/uv
