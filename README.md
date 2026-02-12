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

A research codebase for exploring **electromagnetic ↔ gravitational wave conversion** (Gertsenshtein effect) and potential **amplification mechanisms** in gravity theories with **torsion** (Poincaré gauge theory; parity-even quadratic invariants). The repository includes:

- **A lightweight PDE sandbox** (built on [`py-pde`]) for rapid prototyping and numerics with **850 passing Python tests + ~108 Wolfram tests**.
- A symbolic pipeline (Mathematica + xAct) for **deriving linearized field equations** and exporting them to Python-friendly forms.
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

- **Lagrangian-to-PDE pipeline (`tidal.symbolic`, `tidal.wolfram`)**: complete symbolic-to-numerical pipeline for deriving field equations from Lagrangian densities. Uses Mathematica/xAct for symbolic derivation (Euler-Lagrange equations, linearization, component decomposition) → JSON export → Python/py-pde for dynamic PDE construction and simulation. **Zero hardcoded physics** in numerical layer - all equations derived symbolically. Includes **18 working examples** spanning 1+1D through 3+1D: scalars, vectors, tensors (rank 3+), coupled multi-field systems, curvilinear coordinates, and curved spacetimes. See [examples/README.md](examples/README.md) for complete documentation.
- **CLI (`tidal` command)**: unified command-line interface with 5 subcommands — `tidal derive` (Lagrangian → JSON via TOML config), `tidal simulate` (JSON → PDE simulation with plotting), `tidal inspect` (equation system info), `tidal list` (discover available specs), `tidal validate` (JSON spec validation). Supports `theory.toml` configs with `[[derived_fields]]` for intermediate tensor definitions (e.g., field strength F_ab). Zero new dependencies (stdlib argparse + tomllib).
- **PDE sandbox (`tidal.kgsim`)**: lightweight PDE examples and utilities for experimenting with first-order Klein–Gordon systems built on top of the py-pde library. This repository collects a small simulation toolkit (kgsim) with helpers for grids, initial data, observers, profiling, and runs; 1D & 2D examples with snapshot/video export and animation support.
- **Multi-field coupled systems**: support for N-field coupled Klein–Gordon PDEs with arbitrary mass matrices and coupling terms in both 1D and 2D; includes `multi_gaussian_2d` initializer for spatially separated or overlapping Gaussian pulses; test suite validates symmetry preservation, energy transfer, and decoupled limits.
- **2D coupled simulations**: full support for coupled field evolution in 2D with visualization showing energy transfer between fields; dimension-agnostic PDE implementations work seamlessly with any grid dimensionality.
- **Animation and visualization**: side-by-side spacetime plots (φ(x,t) heatmaps), 1D evolution animations (φ(x) vs time), 2D field evolution with dual-panel animations showing both coupled fields, and high-fps video export (MP4 via ffmpeg or GIF fallback).
- **Dev environment**: container-first, [`uv`] for Python (3.11 pinned), Wolfram Engine 14.3 with xAct tensor framework, optional ffmpeg; Sphinx docs skeleton; type-checked codebase with pytest test suite.
- **Professional development infrastructure**: 850 Python tests + ~108 Wolfram tests, 5 utility scripts for streamlined workflows (`run_wolfram_tests.sh`, `run_examples.sh`, `full_test.sh`, `validate_pipeline.sh`, `lint_wolfram.sh`), comprehensive documentation with module headers and usage strings, robust kernel caching handling for reliable test execution. 0 ruff violations, 0 pyright errors (strict mode).

This README describes the current capabilities, how to run the examples, and planned improvements.

---

# Project Scope and Milestones

This section clarifies where we’re going beyond KG demos.

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

- Comprehensive test suite with ~108 Wolfram unit tests covering all pipeline modules
- 5 utility scripts for workflow automation (test execution, example derivations, pipeline validation)
- Robust kernel caching support ensuring reliable test execution across multiple runs

## Objectives

- Baseline re-derivation of the standard Gertsenshtein effect (Einstein–Maxwell) and its tiny conversion amplitude.
- Extend the gravitational sector to parity-even quadratic PGT with torsion; identify propagating modes and viable parameter windows.
- Linearized PDE system in a flat metric background with constant external magnetic field and (if allowed) homogeneous torsion background. Extract mixing terms.
- Well-posedness: characteristic analysis, hyperbolicity, and causality (characteristic speeds).
- Numerical experiments: 1+1D toy models mapping EM/GR/torsion modes to coupled scalars; verify conversion scaling and stability; then scale up in fidelity.

## Recent Improvements

- **Phase 4-13+ Pipeline Evolution (February 2026)**: ✅ **ALL CRITICAL IMPLEMENTATION COMPLETE**. Fixed Wolfram test symbol conflicts, created 5 development utility scripts, completed module header documentation. **Phase 12**: Auto-computed mass/coupling matrices with symbolic preservation. **Phase 13**: Rank 3+ tensor support. **CLI**: Full `tidal` command with 5 subcommands (derive, inspect, simulate, list, validate) and `theory.toml` support including `[[derived_fields]]`. **Stress tests**: Scalar-vector coupling (mixed-rank cross-field, 4 constants, 4x4 matrices), massive 3-form (rank-3 antisymmetric tensor). **850 Python tests + ~108 Wolfram tests passing**. See [CHANGELOG.md](CHANGELOG.md) for complete history.
- **Lagrangian-to-PDE pipeline (February 2026)**: Complete symbolic derivation pipeline implemented with Mathematica/xAct → JSON → Python/py-pde. Fixed component extraction (free index detection with `IndicesOf[]`), operator identification (laplacian vs identity/mass terms), and RHS extraction. Created EM and Klein-Gordon examples demonstrating massless vs massive field dynamics. Added comprehensive documentation proving zero hardcoded physics in numerical layer. See [CHANGELOG.md](CHANGELOG.md) for detailed Phase 11 implementation notes.
- **2D coupled-field support**: `multi_gaussian_2d` initializer for N-field systems on 2D grids; dimension-agnostic coordinate handling for both 1D and 2D; comprehensive test suite for 2D coupled dynamics (symmetry preservation, energy transfer, decoupled limits).
- **Enhanced coordinate handling**: automatic detection of coordinate array format (flattened vs grid-shaped) in `gaussian_pulse` and `multi_gaussian_2d`; works seamlessly with different py-pde grid configurations.
- **Coupled multi-field PDEs**: `make_coupled_kg_pde` builder for N-field systems with mass/coupling matrices; validation for matrix dimensions, finite masses, and non-negativity.
- **Enhanced plotting**: `choose_writer_and_out` accepts custom output paths and FPS; side-by-side evolution plots with shared colorbars; clean plots with optional axis/tick removal.
- **Profiling and logging**: timer utilities with detailed profiling summaries; logging-based output (no print statements in library code); configurable verbosity.
- **Type safety**: explicit type annotations for observers, trackers, and callbacks; keyword-only boolean arguments to avoid positional ambiguity.
- **Test coverage**: edge-case tests for initial conditions (non-1D/2D grids, non-positive widths, empty amplitudes, mismatched parameter lengths); symmetry tests for coupled-field evolution in both 1D and 2D.

## Future Development

- **Continuous Integration**: GitHub Actions workflow for automated Wolfram test execution on pull requests.
- **Extended physics examples**: Coupled EM/torsion systems, non-Abelian field theories, Yang-Mills gauge theories.
- **Gallery of example runs** with parameter sweeps, convergence studies, and performance benchmarks.
- **Advanced gauge automation**: Lorenz/Coulomb gauge fixing in Wolfram layer.
- **Nonlinear extensions**: Beyond linear perturbation theory.

---

# Development Environment

## Quickstart

This project uses **uv** for Python version/venv/dependencies and is container-first.

```bash
# Ensure Python 3.11 is used (numba/llvmlite friendly)
uv python pin 3.11

# Install runtime + (optional) dev dependencies from pyproject.toml
uv sync --all-extras

# Smoke test: can we import the package and py-pde?
uv run python -c "import tidal, pde; print('OK')"
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

This repository includes a development container configured for Debian. When you open the workspace in a devcontainer (VS Code Remote - Containers / Codespaces) common tools (git, node/npm, eslint, apt, etc.) are already available on PATH which simplifies setup.

All commands from the repo root:

```bash
# 1D KG example: Gaussian pulse spacetime evolution (heatmap)
uv run python examples/klein_gordon/1d_gaussian_pulse.py

# 1D KG animation: time-evolving line plot of φ(x) (MP4 or GIF)
uv run python examples/klein_gordon/1d_gaussian_pulse_anim.py

# 1D KG with potential step (mass variation)
uv run python examples/klein_gordon/1d_step_mass.py

# 1D potential barrier example
uv run python examples/klein_gordon/1d_barrier.py

# 2D KG example: radial ring pulse (writes snapshots and a video/gif)
uv run python examples/klein_gordon/2d_ring_pulse.py

# Coupled two-field symmetric KG system in 1D (side-by-side spacetime plots)
uv run python examples/klein_gordon/2field_coupled.py

# **NEW** Coupled two-field KG system in 2D (animated dual-panel evolution)
uv run python examples/klein_gordon/2d_2field_coupled.py

# **NEW** 3+1D Klein-Gordon example (full 4D spacetime)
uv run python examples/scalar_field_3d/kg_3d_simulation.py
```

Outputs are written to `outputs/` (created automatically if missing).

### Lagrangian-to-PDE Pipeline Examples

The repository includes a complete symbolic-to-numerical pipeline for deriving field equations from Lagrangians and simulating them numerically. **18 examples** cover scalars, vectors, rank-3 tensors, coupled multi-field systems, curvilinear coordinates, and curved spacetimes.

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

| Command | Description |
|---------|-------------|
| `tidal derive theory.toml` | Generate .wls from TOML, run wolframscript to produce JSON |
| `tidal simulate spec.json` | Full simulation with plotting (supports `--param`, `--ic`, `--bc`, `--scheme`) |
| `tidal inspect spec.json` | Display equation system info (fields, operators, parameters) |
| `tidal list` | Discover all available JSON specs in `examples/data/` |
| `tidal validate spec.json` | Validate JSON equation specification structure |

**TOML Configuration** (`theory.toml`):
- Define spacetime dimension, metric, fields, constants, and Lagrangian expression
- `[[derived_fields]]` section for intermediate tensors (e.g., field strength `F_ab = CD[-a][A[-b]] - CD[-b][A[-a]]`)
- Runtime parameters with default values in `[parameters]` section

**Pipeline Examples:**

| Example | Dim | Key Features |
|---------|-----|--------------|
| `scalar_field/` | 1+1D | Klein-Gordon, mass term, dispersion |
| `electromagnetic/` | 1+1D | Maxwell, Lorenz gauge, massless waves |
| `proca/` | 1+1D | Massive vector field (Proca mass) |
| `coupled_scalars/` | 1+1D | Cross-field coupling, mass matrix, energy transfer |
| `chern_simons/` | 2+1D | Epsilon tensor, topological mass, A_0 constraint |
| `elasticity/` | 2+1D | Anisotropic laplacian, cross_derivative_xy |
| `curved_spacetime/` | 2+1D | De Sitter, Hubble friction, time-dependent coefficients |
| `sphere_kg/` | 2+1D | KG on S², position-dependent coefficients |
| `polar_kg/` | 2+1D | Polar coordinates, Christoffel auto-detection |
| `electrostatics/` | 2+1D | Poisson equation, constraint solver |
| `scalar_vector_coupling/` | 2+1D | Mixed-rank cross-field (scalar+vector), 4 constants, CS+coupling |
| `scalar_field_3d/` | 3+1D | Full 4D KG, 32^3 grid |
| `spherical_kg/` | 3+1D | Spherical coordinates, trig coefficients |
| `cylindrical_kg/` | 3+1D | Cylindrical coordinates, mixed curved/flat |
| `gravitational_waves/` | 3+1D | xPert linearization, TT gauge, constraints |
| `massive_3form/` | 3+1D | Rank-3 antisymmetric tensor, symmetry reduction |
| `massive_gravity/` | 2+1D | Linearized massive gravity, Fierz-Pauli mass, xPert, coupled constraints |
| `coupled_proca/` | 2+1D | Two massive vectors, coupled Helmholtz constraints, Dirichlet BCs |

See [examples/README.md](examples/README.md) for complete documentation and verification that the Python layer contains zero hardcoded physics.

### Animation Features

- **1D animations** (`1d_gaussian_pulse_anim.py`): show φ(x) vs x at each time step; configurable FPS and snapshot interval; supports both MP4 (ffmpeg) and GIF (Pillow).
- **2D single-field animations** (`2d_ring_pulse.py`): radial collapse/expansion of a ring pulse; imshow-based heatmaps written frame-by-frame with tqdm progress bars.
- **2D coupled-field animations** (`2d_2field_coupled.py`): dual-panel visualization showing both fields evolving simultaneously with shared color scale; demonstrates energy transfer between coupled fields with different masses; spatially separated initial Gaussians show propagation and interaction.
- **1D coupled-field plots** (`2field_coupled.py`): two side-by-side spacetime heatmaps (one per field) with a shared colorbar; clean axis formatting; test for symmetry preservation.

### Customization

- **FPS control**: pass `fps=<int>` to `choose_writer_and_out` (default: auto-calculated from snapshot count and t_end).
- **Snapshot resolution**: set `snapshot_interval` in the `run` call to match your desired temporal resolution (e.g., `dt` for every integrator step).
- **Output format**: MP4 if ffmpeg is available; GIF fallback otherwise. Specify output path (without extension) to `choose_writer_and_out`.

## (Optional) Video Support

For MP4 via Matplotlib’s FFMpegWriter:

```bash
# inside the dev container
sudo apt-get update && sudo apt-get install -y ffmpeg
```

If `ffmpeg` is unavailable, the example falls back to a GIF via Pillow.

## Tests

The project includes a comprehensive test suite with **850 Python tests + ~108 Wolfram tests**.

### Python Tests (850 tests)

```bash
# Run all Python tests with pytest
uv run pytest -v

# Run a specific test module
uv run pytest tests/test_py_pde_smoke.py -v

# Run with coverage report (HTML)
uv run pytest --cov=tidal --cov-report=html
open htmlcov/index.html  # View detailed HTML report

# Run with coverage report (terminal)
uv run pytest --cov=tidal --cov-report=term-missing

# Run with coverage report (XML for CI)
uv run pytest --cov=tidal --cov-report=xml
```

### Wolfram Tests (~108 tests)

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
- Edge cases (empty grids, invalid bounds, division by zero)
- Path traversal protection and validation

See [`scripts/README.md`](scripts/README.md) for detailed documentation of utility scripts.

---

## Documentation

The repo builds Sphinx docs and deploys to GitHub Pages via Actions.

### Local Build:

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
- **Test requirements**: All changes must maintain 100% test pass rate (850 Python + ~108 Wolfram tests). New features require corresponding unit tests in both Python and Wolfram layers where applicable.
- Run `./scripts/full_test.sh` before submitting PRs to verify all tests pass.
- Follow the project's type-checking and linting conventions (keyword-only booleans, explicit type annotations, no print in library code).

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

This project builds on:

- [`py-pde`](https://py-pde.readthedocs.io/) for the PDE solver framework.
- [`uv`](https://github.com/astral-sh/uv) for fast Python environment management.
- The [xAct/xTensor ecosystem](http://www.xact.es/) for symbolic tensor algebra (xCore, xPerm, xTensor, xCoba) powering the complete Lagrangian-to-PDE derivation pipeline.

[`py-pde`]: https://py-pde.readthedocs.io/
[`uv`]: https://github.com/astral-sh/uv
