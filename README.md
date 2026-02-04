# Torsion Gertsenshtein

View the `torsion_gertsenshtein` package documentation [here](https://williamroyce.github.io/torsion-gertsenshtein/).

A research codebase for exploring **electromagnetic ↔ gravitational wave conversion** (Gertsenshtein effect) and potential **amplification mechanisms** in gravity theories with **torsion** (Poincaré gauge theory; parity-even quadratic invariants). The repository will host:

- A lightweight PDE sandbox (built on [`py-pde`]) for rapid prototyping and numerics.
- A symbolic pipeline (Mathematica + xAct) for **deriving linearized field equations** and exporting them to Python-friendly forms.
- Documentation and experiments for **mixing mechanisms** and **hyperbolicity/causality checks** relevant to the effect.

> TL;DR: start with Klein–Gordon toy systems in 1+1D → grow to coupled EM/metric/torsion perturbations → test conversion and stability in controlled scenarios.

---

## Current Status (usable today)

- **Lagrangian-to-PDE pipeline (`torsion_gertsenshtein.symbolic`, `torsion_gertsenshtein.wolfram`)**: complete symbolic-to-numerical pipeline for deriving field equations from Lagrangian densities. Uses Mathematica/xAct for symbolic derivation (Euler-Lagrange equations, linearization, component decomposition) → JSON export → Python/py-pde for dynamic PDE construction and simulation. **Zero hardcoded physics** in numerical layer - all equations derived symbolically. Includes working examples: EM field (massless waves) and Klein-Gordon (massive scalar field with dispersion). See [examples/README.md](examples/README.md) for complete documentation.
- **PDE sandbox (`torsion_gertsenshtein.kgsim`)**: lightweight PDE examples and utilities for experimenting with first-order Klein–Gordon systems built on top of the py-pde library. This repository collects a small simulation toolkit (kgsim) with helpers for grids, initial data, observers, profiling, and runs; 1D & 2D examples with snapshot/video export and animation support.
- **Multi-field coupled systems**: support for N-field coupled Klein–Gordon PDEs with arbitrary mass matrices and coupling terms in both 1D and 2D; includes `multi_gaussian_2d` initializer for spatially separated or overlapping Gaussian pulses; test suite validates symmetry preservation, energy transfer, and decoupled limits.
- **2D coupled simulations**: full support for coupled field evolution in 2D with visualization showing energy transfer between fields; dimension-agnostic PDE implementations work seamlessly with any grid dimensionality.
- **Animation and visualization**: side-by-side spacetime plots (φ(x,t) heatmaps), 1D evolution animations (φ(x) vs time), 2D field evolution with dual-panel animations showing both coupled fields, and high-fps video export (MP4 via ffmpeg or GIF fallback).
- **Dev environment**: container-first, [`uv`] for Python (3.11 pinned), Wolfram Engine 14.3 with xAct tensor framework, optional ffmpeg; Sphinx docs skeleton; type-checked codebase with pytest test suite.

This README describes the current capabilities, how to run the examples, and planned improvements.

---

# Project Scope and Milestones

This section clarifies where we’re going beyond KG demos.

## Objectives

- Baseline re-derivation of the standard Gertsenshtein effect (Einstein–Maxwell) and its tiny conversion amplitude.
- Extend the gravitational sector to parity-even quadratic PGT with torsion; identify propagating modes and viable parameter windows.
- Linearized PDE system in a flat metric background with constant external magnetic field and (if allowed) homogeneous torsion background. Extract mixing terms.
- Well-posedness: characteristic analysis, hyperbolicity, and causality (characteristic speeds).
- Numerical experiments: 1+1D toy models mapping EM/GR/torsion modes to coupled scalars; verify conversion scaling and stability; then scale up in fidelity.

## Recent Improvements

- **Lagrangian-to-PDE pipeline (February 2026)**: Complete symbolic derivation pipeline implemented with Mathematica/xAct → JSON → Python/py-pde. Fixed component extraction (free index detection with `IndicesOf[]`), operator identification (laplacian vs identity/mass terms), and RHS extraction. Created EM and Klein-Gordon examples demonstrating massless vs massive field dynamics. Added comprehensive documentation proving zero hardcoded physics in numerical layer. See [CHANGELOG.md](CHANGELOG.md) for detailed Phase 11 implementation notes.
- **2D coupled-field support**: `multi_gaussian_2d` initializer for N-field systems on 2D grids; dimension-agnostic coordinate handling for both 1D and 2D; comprehensive test suite for 2D coupled dynamics (symmetry preservation, energy transfer, decoupled limits).
- **Enhanced coordinate handling**: automatic detection of coordinate array format (flattened vs grid-shaped) in `gaussian_pulse` and `multi_gaussian_2d`; works seamlessly with different py-pde grid configurations.
- **Coupled multi-field PDEs**: `make_coupled_kg_pde` builder for N-field systems with mass/coupling matrices; validation for matrix dimensions, finite masses, and non-negativity.
- **Enhanced plotting**: `choose_writer_and_out` accepts custom output paths and FPS; side-by-side evolution plots with shared colorbars; clean plots with optional axis/tick removal.
- **Profiling and logging**: timer utilities with detailed profiling summaries; logging-based output (no print statements in library code); configurable verbosity.
- **Type safety**: explicit type annotations for observers, trackers, and callbacks; keyword-only boolean arguments to avoid positional ambiguity.
- **Test coverage**: edge-case tests for initial conditions (non-1D/2D grids, non-positive widths, empty amplitudes, mismatched parameter lengths); symmetry tests for coupled-field evolution in both 1D and 2D.

## Future Aims / TODOs

- Implement a Numba RHS for `KleinGordonPDE` and `make_coupled_kg_pde` to enable the `numba` backend for coupled systems.
- Expand test suite: observers, profiling utilities, boundary condition handling.
- Improve and publish type stubs for py-pde usages or vendor a narrow Protocol for solver/field interfaces to reduce casts.
- Add CI (GitHub Actions) to run linting, type-checking (pyright/mypy), tests, and build docs on every push.
- Expand initial condition library: plane waves, solitons, custom profiles.
- Gallery of example runs with parameter sweeps and convergence studies.

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
uv run python -c "import torsion_gertsenshtein, pde; print('OK')"
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
```

Outputs are written to `outputs/` (created automatically if missing).

### Lagrangian-to-PDE Pipeline Examples

The repository includes a complete symbolic-to-numerical pipeline for deriving field equations from Lagrangians and simulating them numerically:

```bash
# Full pipeline demonstration (both EM and Klein-Gordon)
bash examples/demo_full_pipeline.sh

# EM field: Stage 1 (derive from Lagrangian L = -1/4 F_μν F^μν)
cd examples/electromagnetic
wolframscript -file em_lagrangian_1d.wls

# EM field: Stage 2 (simulate massless wave propagation)
uv run python examples/electromagnetic/em_from_lagrangian.py

# Klein-Gordon: Stage 1 (derive from Lagrangian L = -1/2 (∂φ)² - 1/2 m²φ²)
cd examples/scalar_field
wolframscript -file klein_gordon.wls

# Klein-Gordon: Stage 2 (simulate massive field with dispersion)
uv run python examples/scalar_field/kg_from_lagrangian.py
```

**Key Features:**
- **No hardcoded physics**: All equations dynamically loaded from JSON
- **Symbolic derivation**: Mathematica/xAct computes Euler-Lagrange equations
- **JSON interface**: Well-defined schema for equation specification
- **Dynamic PDE construction**: Python builds solver from specification
- **Verified examples**: EM (massless) vs Klein-Gordon (massive) demonstrate different physics from different Lagrangians

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

A minimal smoke test suite is included to ensure `py-pde` and the solver path work, and to validate edge cases in initial conditions and multi-field systems:

```bash
# run all tests with pytest
uv run pytest -v

# run a specific test module
uv run pytest tests/test_py_pde_smoke.py -v

# run with coverage report
uv run pytest --cov=torsion_gertsenshtein --cov-report=html
```

Add more tests under `tests/` (ICs, utils, observers).

---

## Documentation

The repo builds Sphinx docs and deploys to GitHub Pages via Actions.

### Local Build:

```bash
# auto-generate API docs
uv run sphinx-apidoc --force --module-first -o docs/source/ torsion_gertsenshtein/

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

## Troubleshooting

- **Import errors in VS Code** (e.g., numpy not found): ensure the interpreter is the repo's venv (`.venv/bin/python3`), then reload the window.
- **`llvmlite/numba` build failures**: stick to **Python 3.11** (`uv python pin 3.11`).
- **FileNotFoundError: ffmpeg** — install ffmpeg (see apt command above) or let the example produce a GIF.
- **Type-checker warnings about third-party stubs** — run examples anyway; code uses TYPE_CHECKING guards and runtime-safe casts where necessary.
- **Pages 404 or deploy errors**: ensure Settings → Pages → Source = GitHub Actions and Actions → Workflow permissions = Read/Write.
- **Animation has low frame count**: ensure `snapshot_interval` in the `run` call matches your desired temporal resolution (e.g., set to `dt` for every integrator step). Increase `fps` in `choose_writer_and_out` for smoother playback.
- **Logging messages not visible**: call `logging.basicConfig(level=logging.INFO, ...)` at the start of your script or in `main()`.

---

## Contributing

- Open an issue or submit a PR.
- If adding features that touch numerical kernels, include unit/regression tests.
- Follow the project's type-checking and linting conventions (keyword-only booleans, explicit type annotations, no print in library code).

---

## License

No LICENSE file included in this repository. Add a LICENSE file if you intend to open-source this work, before distributing artifacts or accepting external contributions.

---

## Acknowledgements

This project builds on:

- [`py-pde`](https://py-pde.readthedocs.io/) for the PDE solver framework.
- [`uv`](https://github.com/astral-sh/uv) for fast Python environment management.
- The xAct/xTensor ecosystem (for symbolic tensor algebra, not yet integrated).

[`py-pde`]: https://py-pde.readthedocs.io/
[`uv`]: https://github.com/astral-sh/uv
