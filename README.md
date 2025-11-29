# Torsion Gertsenshtein

View the `torsion_gertsenshtein` package documentation [here](https://williamroyce.github.io/torsion-gertsenshtein/).

A research codebase for exploring **electromagnetic ↔ gravitational wave conversion** (Gertsenshtein effect) and potential **amplification mechanisms** in gravity theories with **torsion** (Poincaré gauge theory; parity-even quadratic invariants). The repository will host:

- A lightweight PDE sandbox (built on [`py-pde`]) for rapid prototyping and numerics.
- A symbolic pipeline (Mathematica + xAct) for **deriving linearized field equations** and exporting them to Python-friendly forms.
- Documentation and experiments for **mixing mechanisms** and **hyperbolicity/causality checks** relevant to the effect.

> TL;DR: start with Klein–Gordon toy systems in 1+1D → grow to coupled EM/metric/torsion perturbations → test conversion and stability in controlled scenarios.

---

## Current Status (usable today)

- **PDE sandbox (`torsion_gertsenshtein.kgsim`)**: lightweight PDE examples and utilities for experimenting with a first-order Klein–Gordon system built on top of the py-pde library. This repository collects a small simulation toolkit (kgsim) with helpers for grids, initial data, observers, and runs; 1D & 2D examples and snapshot/video export used during development.
- **Dev environment**: container-first, [`uv`] for Python (3.11 pinned), optional ffmpeg; Sphinx docs skeleton; basic tests.

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

## Future Aims / TODOs

- Implement a Numba RHS for `KleinGordonPDE` to enable the `numba` backend.
- Add a small test suite (initial_conditions, utils, observers).
- Improve and publish type stubs for py-pde usages or vendor a narrow Protocol for solver/field interfaces to reduce casts.
- Add CI (GitHub Actions) to run linting, tests, and build docs.
- Expand initial condition library and example gallery.

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

## Running the Examples

This repository includes a development container configured for Debian. When you open the workspace in a devcontainer (VS Code Remote - Containers / Codespaces) common tools (git, node/npm, eslint, apt, etc.) are already available on PATH which simplifies setup.

All commands from the repo root:

```bash
# 1D example
uv run python examples/klein_gordon/1d_gaussian_pulse.py

# 2D example (writes snapshots and a video/gif)
uv run python examples/klein_gordon/2d_ring_pulse.py
```

Outputs are written to `outputs/`.

## (Optional) Video Support

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

---

## Documentation

The repo builds Sphinx docs and deploys to GitHub Pages via Actions.

### Local Build:

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
https://williamroyce.github.io/torsion-gertsenshtein/
```

Use `$BROWSER <url>` (from within the devcontainer) to open the project documentation link in the host browser.

---

## Troubleshooting

- Import errors in VS Code (e.g., numpy not found): ensure the interpreter is the repo’s venv (`.venv/bin/python3`), then reload the window.
- `llvmlite/numba` build failures: stick to **Python 3.11** (`uv python pin 3.11`).
- FileNotFoundError: ffmpeg — install ffmpeg (see apt command above) or let the example produce a GIF.
- Type-checker warnings about third-party stubs — run examples anyway; code uses TYPE_CHECKING guards and runtime-safe casts where necessary.
- Pages 404 or deploy errors: ensure Settings → Pages → Source = GitHub Actions and Actions → Workflow permissions = Read/Write.

## Contributing

- Open an issue or submit a PR.
- If adding features that touch numerical kernels, include unit/regression tests.

## License

No LICENSE file included in this repository. Add a LICENSE file if you intend to open-source this work, before distributing artifacts or accepting external contributions.
