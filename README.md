# Project Template

View the `torsion_gertsenshtein` package documentation [here](https://williamroyce.github.io/torsion-gertsenshtein/).

Lightweight PDE examples and utilities for experimenting with a Klein–Gordon system built on top of the py-pde library. This repository collects a small simulation toolkit (kgsim) plus runnable examples (1D and 2D) used during development.

This README describes the current capabilities, how to run the examples, and planned improvements.

## Current functionality

Package: `torsion_gertsenshtein.kgsim`

- `config.py`
  - `KGParameters` — physics parameters (mass).
  - `GridConfig`, `SimulationConfig` — concise dataclasses for grid and runtime settings.
- `grids.py`
  - `make_grid(cfg)` — construct a `py-pde` CartesianGrid from `GridConfig`.
- `initial_conditions.py`
  - `gaussian_pulse(...)` — 1D Gaussian initial condition (φ, π).
  - `ring_pulse_2d(...)` — 2D ring Gaussian initial condition (φ, π).
- `equations.py`
  - `KleinGordonPDE` — Klein–Gordon PDE in first-order form (φ, π).
- `utils.py`
  - small helpers for BC inference, coordinate helpers, and typed arithmetic helpers (e.g. `mul_scalar_field`, `sub_scalar_fields`).
- `observers.py`
  - `total_energy_observer(mass)` — computes total field energy; intended for trackers.
- `runners.py`
  - `run(...)` — convenience entry to run a simulation, attach trackers/observers, and normalize solver returns.
  - Supports `extra_observer` and `snapshot_interval` (record snapshots at a fixed cadence).

Examples (under `examples/`):

- `klein_gordon/1d_gaussian_pulse.py` — 1D evolution and static plots.
- `klein_gordon/2d_ring_pulse.py` — 2D ring pulse example with snapshot collection and video export.
- `py_pde_laplacian_demo.py` — simple demo adapted from py-pde.

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

## Running the examples

This repository includes a development container configured for Debian. When you open the workspace in a devcontainer (VS Code Remote - Containers / Codespaces) common tools (git, node/npm, eslint, apt, etc.) are already available on PATH which simplifies setup.

All commands from the repo root:

```bash
# 1D example
uv run python examples/klein_gordon/1d_gaussian_pulse.py

# 2D example (writes snapshots and a video/gif)
uv run python examples/klein_gordon/2d_ring_pulse.py
```

Outputs are written to `outputs/`.

## (Optional) video support

For MP4 via Matplotlib’s FFMpegWriter:

```bash
# inside the dev container
sudo apt-get update && sudo apt-get install -y ffmpeg
```

If `ffmpeg` is unavailable, the example falls back to a GIF via Pillow.

## Tests

A minimal smoke test is included to ensure `py-pde` and the solver path work:

```bash
uv run pytest -q
```

Add more tests under `tests/` (ICs, utils, observers).

## Documentation

The repo builds Sphinx docs and deploys to GitHub Pages via Actions.

### Local build:

```bash
# auto-generate API docs
uv run sphinx-apidoc --force --module-first -o docs/source/ torsion_gertsenshtein/

# build HTML
( cd docs && make html )

# open locally
python -m http.server -d docs/build/html 8000
```

On push to `main`, CI builds the docs and deploys to:

```bash
$BROWSER https://williamroyce.github.io/torsion-gertsenshtein/
```

Use `$BROWSER <url>` (from within the devcontainer) to open the project documentation link in the host browser.

## Troubleshooting

- Import errors in VS Code (e.g., numpy not found): ensure the interpreter is the repo’s venv (`.venv/bin/python3`), then reload the window.
- `llvmlite/numba` build failures: stick to **Python 3.11** (`uv python pin 3.11`).
- FileNotFoundError: ffmpeg — install ffmpeg (see apt command above) or let the example produce a GIF.
- Type-checker warnings about third-party stubs — run examples anyway; code uses TYPE_CHECKING guards and runtime-safe casts where necessary.
- Pages 404 or deploy errors: ensure Settings → Pages → Source = GitHub Actions and Actions → Workflow permissions = Read/Write.

## Future aims / TODOs

- Implement a Numba RHS for `KleinGordonPDE` to enable the `numba` backend.
- Add a small test suite (initial_conditions, utils, observers).
- Improve and publish type stubs for py-pde usages or vendor a narrow Protocol for solver/field interfaces to reduce casts.
- Add CI (GitHub Actions) to run linting, tests, and build docs.
- Expand initial condition library and example gallery.

## Contributing

- Open an issue or submit a PR.
- If adding features that touch numerical kernels, include unit/regression tests.

## License

No license file included in this repository. Add a LICENSE file if you intend to open-source this work.
