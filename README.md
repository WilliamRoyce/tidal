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

A symbolic-to-numerical framework for **linearized field theory** — define a Lagrangian, derive the PDEs automatically, and simulate. Built for the **Gertsenshtein effect** (electromagnetic ↔ gravitational wave conversion; [Domcke & Garcia-Cely 2023](https://arxiv.org/abs/2301.02072)) and **torsion wave physics** (Poincare gauge theory). The graviton-photon conversion probability P = sin²(κB₀D/2) has been validated numerically to 0.04% against the Boccaletti (1970) formula. The repository includes:

- A **symbolic derivation pipeline** (Mathematica/xAct) that derives linearized field equations from any Lagrangian and exports them as JSON specifications — zero hardcoded physics.
- **Five solver backends** (SUNDIALS IDA/CVODE, Fourier modal, leapfrog, scipy) with analytical Jacobians, FFT spectral operators, and 2nd/4th/6th-order FD stencils.
- **20 working examples** spanning 1+1D to 3+1D: scalars, vectors, rank-3 tensors, coupled multi-field systems, curvilinear coordinates, curved spacetimes, background-field scattering, graviton-photon conversion, and graviton-torsion mixing.
- **1,701 Python tests + ~121 Wolfram tests**, 0 ruff violations, 0 pyright errors (strict mode).

> Define a Lagrangian in TOML → derive linearized PDEs symbolically → simulate with adaptive solvers → measure conversion, spectra, and scattering.

---

## Community & Support

- **Questions & Ideas**: [GitHub Discussions](https://github.com/WilliamRoyce/torsion-gertsenshtein/discussions) — Ask questions, share use cases, discuss physics
- **Bug Reports**: [Issue Tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues) — Report bugs with the `bug` label
- **Feature Requests**: [Issue Tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues) — Propose features with the `enhancement` label
- **Contributing**: See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup and workflow

For more, visit the [Documentation](https://williamroyce.github.io/torsion-gertsenshtein/).

---

## Capabilities

### Symbolic Pipeline (`tidal.symbolic`, `tidal.wolfram`)

Complete Lagrangian-to-PDE derivation: TOML config → Mathematica/xAct (Euler-Lagrange, linearization via xPert, component decomposition) → JSON specification → Python solvers. Supports `[[derived_fields]]` (e.g., field strength tensors), `[[background_fields]]` (external magnetic fields, potentials), `[[gauge]]` (Lorenz, de Donder, Coulomb, temporal, axial), and `[torsion]` (Poincare gauge theory with propagating torsion). All equations derived symbolically — the numerical layer contains zero physics.

### Solver Architecture (`tidal.solver`)

Five time-integration backends with automatic selection based on equation structure:

| Backend | Library | Use Case | Key Feature |
| ------- | ------- | -------- | ----------- |
| **IDA** | SUNDIALS | DAE (algebraic constraints) | Implicit Newton, 3-tier analytical Jacobian |
| **CVODE** | SUNDIALS | Adaptive ODE (waves) | BDF, tolerance control, sparse Jacobian |
| **Modal** | numpy/scipy | Exact spectral (periodic, time-independent) | Machine-precision via eigendecomposition |
| **Leapfrog** | numpy | Symplectic | Exact energy conservation (Yoshida 4th-order) |
| **scipy** | scipy.integrate | General-purpose | DOP853, Radau, BDF via `solve_ivp` |

Spatial operators: 2nd/4th/6th-order finite-difference stencils + FFT spectral operators (auto-enabled for all-periodic BCs). Three-tier constraint IC pre-solve (FFT → sparse matrix → automatic). Analytical Jacobian with three active tiers: dense (N ≤ 2K), sparse CSC with SuperLU_MT (N ≤ 200K), and GMRES with JVP (N > 200K).

### CLI (`tidal` command)

Unified command-line interface with 10 subcommands:

| Command | Description |
| ------- | ----------- |
| `tidal derive theory.toml` | Generate .wls from TOML, run wolframscript to produce JSON |
| `tidal simulate spec.json` | Full simulation with plotting (`--param`, `--ic`, `--bc`, `--scheme`) |
| `tidal measure result_dir/` | Extract physics measurements (energy, conversion, mixing, spectra) |
| `tidal inspect spec.json` | Display equation system info (fields, operators, parameters) |
| `tidal list` | Discover all available JSON specs in `examples/data/` |
| `tidal validate spec.json` | Validate JSON equation specification structure |
| `tidal plot result_dir/` | Standalone plotting from simulation output directories |
| `tidal sweep spec.json` | Parameter sweeps, convergence studies, adaptive sampling |
| `tidal analyze sweep_dir/` | Post-hoc sensitivity analysis (Sobol/Morris) of sweep results |
| `tidal doctor` | Environment diagnostics (Wolfram, dependencies, xAct) |

Supports `theory.toml` configs with `[[derived_fields]]`, `[[background_fields]]`, optional `[[gauge]]`, and `[torsion]` sections. Zero new dependencies (stdlib argparse + tomllib).

### Measurement Module (`tidal.measurement`)

13 post-hoc analysis types: summary, energy, conversion P(t), mixing length, spectrum, spectral conversion P(k,t), dispersion ω(k), conservation diagnostics, effective mass, asymptotic scattering, peak conversion, group/phase velocity, and resonance analysis. Critical field analysis for theory comparison (amplification factors). Disk-backed snapshot storage for long simulations via `SnapshotWriter`.

---

## Research Results

### Gertsenshtein Effect (Validated)

The standard graviton-photon conversion has been derived from the Einstein-Maxwell Lagrangian L = (1/κ²)R − (1/4)F² via the TIDAL pipeline and validated:

- **Uniform B-field**: P = sin²(κB₀D/2) confirmed to 0.36% RMS across a 40-point B₀ sweep (N=1024, κ=1)
- **Localized Gaussian B-field**: Boccaletti (1970) formula P = sin²(κ/2 × ∫B₀ dz) confirmed to 0.04% across a 48-point 2D sweep
- **Literature comparison**: Identifies and documents the √(4π) normalization error in Palessandro & Rothman (2023); independently confirmed by Dandoy, Lella et al. (2024)

See `docs/tex/gertsenshtein.tex`, `docs/tex/gertsenshtein_formula.tex`, and `docs/tex/gertsenshtein_localized.tex` for the full physics, derivation, and validation.

### Torsion (In Progress)

Poincare gauge theory with propagating torsion via `[torsion]` TOML section. The `graviton_torsion/` example derives the full 3+1D PGT Lagrangian (R + α₁T² + α₂T² + α₃T²) with torsion perturbations alongside metric perturbations. See `docs/tex/torsion.tex`.

### Remaining Research Targets

- **Plasma detuning** (Phase F1): Requires gauge-invariant photon mass mechanism — blocked by gauge-potential coupling artifact
- **Magnetar/FRB scattering** (Phase F3): Dipolar B(r) ∝ 1/r³ in radial coordinates
- **Absorbing boundaries** (Phase G): Sponge layers and PML for finite-magnet interaction regions
- **Torsion-EM mixing**: Full graviton-torsion-photon conversion in background magnetic field — the project's ultimate goal

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

The repository includes **20 working examples** covering scalars, vectors, rank-3 tensors, coupled multi-field systems, curvilinear coordinates, curved spacetimes, background-field scattering, and graviton-photon/torsion conversion.

```bash
# Each example has a run.sh showing the full derive → inspect → simulate workflow:
cd examples/coupled_scalars && bash run.sh

# Or use the CLI directly:
tidal derive examples/coupled_scalars/theory.toml    # derive equations from Lagrangian
tidal inspect examples/data/coupled_scalars.json      # inspect equation structure
tidal simulate examples/data/coupled_scalars.json \   # simulate
  --param m2=1.0 --ic gaussian --t-end 20
tidal list                                             # discover all available JSON specs
tidal validate examples/data/coupled_scalars.json      # validate JSON spec structure
```

**TOML Configuration** (`theory.toml`):

- Define spacetime dimension, metric, fields, constants, and Lagrangian expression
- `[[derived_fields]]` section for intermediate tensors (e.g., field strength `F_ab = CD[-a][A[-b]] - CD[-b][A[-a]]`)
- `[[background_fields]]` for external fields (magnetic field, potentials) with position-dependent coefficients
- `[[gauge]]` for per-field gauge fixing (Lorenz, de Donder, Coulomb, temporal, axial)
- `[torsion]` for Poincare gauge theory with propagating torsion
- Runtime parameters with default values in `[parameters]` section

**Pipeline Examples:**

| Example | Dim | Key Features |
| ------- | --- | ------------ |
| `chern_simons/` | 2+1D | Epsilon tensor, topological mass, A_0 constraint |
| `coupled_proca/` | 2+1D | Two massive vectors, coupled Helmholtz constraints, periodic BCs |
| `coupled_scalars/` | 1+1D | Cross-field coupling, mass matrix, energy transfer |
| `coupled_scattering/` | 2+1D | Position-dependent Gaussian coupling, background fields, wave scattering |
| `curved_spacetime/` | 2+1D | De Sitter, Hubble friction, time-dependent coefficients |
| `cylindrical_kg/` | 3+1D | Cylindrical coordinates, mixed curved/flat |
| `cylindrical_kg_1d/` | 1+1D | Cylindrical coordinates, plane-wave dimensional reduction |
| `elasticity/` | 2+1D | Anisotropic laplacian, cross_derivative_xy |
| `gertsenshtein/` | 1+1D | Einstein-Maxwell graviton-photon conversion, multi-field perturbation |
| `gravitational_waves/` | 3+1D | xPert linearization, TT gauge, constraints |
| `gravitational_waves_1d/` | 1+1D | Linearized gravity, plane-wave 1D reduction |
| `graviton_torsion/` | 3+1D | PGT Lagrangian, torsion perturbations, graviton-torsion mixing |
| `massive_3form/` | 3+1D | Rank-3 antisymmetric tensor, symmetry reduction |
| `massive_gravity/` | 2+1D | Linearized massive gravity, Fierz-Pauli mass, xPert, coupled constraints |
| `polar_kg/` | 2+1D | Polar coordinates, Christoffel auto-detection |
| `proca_background/` | 2+1D | Lorentzian scalar background, two Proca vectors, constraint+BG integration |
| `scalar_potential_well/` | 1+1D | Background potential well, `[[background_fields]]`, bound states |
| `scalar_vector_coupling/` | 2+1D | Mixed-rank cross-field (scalar+vector), 4 constants, CS+coupling |
| `sphere_kg/` | 2+1D | KG on S², position-dependent coefficients |
| `spherical_kg_1d/` | 1+1D | Spherical coordinates, plane-wave dimensional reduction |

## (Optional) Video Support

For MP4 via Matplotlib's FFMpegWriter:

```bash
# inside the dev container
sudo apt-get update && sudo apt-get install -y ffmpeg
```

If `ffmpeg` is unavailable, the example falls back to a GIF via Pillow.

## Tests

The project includes a comprehensive test suite with **1,701 Python tests + ~121 Wolfram tests**.

### Python Tests (1,701 tests)

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

### Wolfram Tests (~121 tests)

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

### LaTeX Technical Documentation (`docs/tex/`)

25 self-contained LaTeX fragments covering physics, architecture, features, and operational guides. Each file uses shared macros from `preamble.tex` and can be included in an Overleaf report via `\input{fragment_name}`. See [`docs/README.md`](docs/README.md) for the complete index.

### Sphinx API Documentation

The repo builds Sphinx docs and deploys to GitHub Pages via Actions.

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

The verification script checks Wolfram Engine activation, xAct package installation (xCore, xPerm, xTensor, xCoba, xPert), xPerm binary GLIBC compatibility, and runs a full smoke test with tensor operations.

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
- **Environment diagnostics**: run `tidal doctor` for a comprehensive check of Wolfram, Python dependencies, and xAct installation.

See `docs/tex/troubleshooting.tex` for a comprehensive error encyclopedia covering Wolfram/xAct and Python solver issues.

---

## Contributing

- Open an issue or submit a PR.
- **Test requirements**: All changes must maintain 100% test pass rate (1,701 Python + ~121 Wolfram tests). New features require corresponding unit tests in both Python and Wolfram layers where applicable.
- Run `./scripts/full_test.sh` before submitting PRs to verify all tests pass.
- Follow the project's type-checking and linting conventions (keyword-only booleans, explicit type annotations, no print in library code).

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

This project builds on:

- [SUNDIALS](https://computing.llnl.gov/projects/sundials) — IDA (DAE) and CVODE (BDF) solvers via [scikit-sundae](https://github.com/NREL/scikit-sundae) (Hindmarsh et al. 2005).
- The [xAct/xTensor ecosystem](http://www.xact.es/) — symbolic tensor algebra (Martin-Garcia et al.) powering the Lagrangian-to-PDE derivation pipeline.
- [xPert](https://www.researchgate.net/publication/1740524) — metric perturbation theory (Brizuela et al. 2009) for linearization.
- [`uv`](https://github.com/astral-sh/uv) — fast Python environment management.
- Originally built on [py-pde](https://py-pde.readthedocs.io/) (Zwicker, JOSS 2020); finite-difference stencil conventions retained in TIDAL's native operators.

Design decisions are informed by [Dedalus](https://arxiv.org/abs/1905.10388) (Burns et al. 2020), [MEEP](https://meep.readthedocs.io/) (Oskooi et al. 2010), and [FEniCS](https://fenicsproject.org/) (Baratta et al. 2023). The core physics targets the Gertsenshtein effect (Gertsenshtein 1962; [Domcke & Garcia-Cely 2023](https://arxiv.org/abs/2301.02072)). See [`docs/references.md`](docs/references.md) and `docs/tex/references.bib` for the full citation list.

[`uv`]: https://github.com/astral-sh/uv
