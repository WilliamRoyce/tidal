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

## Running the examples

This repository includes a development container configured for Debian. When you open the workspace in a devcontainer (VS Code Remote - Containers / Codespaces) common tools (git, node/npm, eslint, apt, etc.) are already available on PATH which simplifies setup.

Quick steps (recommended)

- Open the project in the devcontainer (recommended). The container already provides many CLI tools; you only need to install Python runtime dependencies.

- Install Python dependencies (inside the devcontainer). You can use the container Python directly or create a virtualenv:

```bash
# optional: create a venv if you prefer isolation
python -m venv .venv
source .venv/bin/activate

# install minimal runtime deps
pip install --upgrade pip
pip install py-pde numpy matplotlib
# optional for speed / video output:
pip install numba imageio[ffmpeg] pillow
```

- (Optional) If you want MP4 output from Matplotlib's FFMpegWriter, ensure ffmpeg is installed in the container:

```bash
# inside devcontainer
sudo apt update && sudo apt install -y ffmpeg
```

- Run an example (from repository root):

```bash
python examples/klein_gordon/1d_gaussian_pulse.py
python examples/klein_gordon/2d_ring_pulse.py
```

Outputs (images / video) are written to `outputs/` by the examples.

Notes

- If ffmpeg is missing, the 2D example falls back to GIF output via Pillow.
- Use `$BROWSER <url>` (from within the devcontainer) to open the project documentation link in the host browser, e.g.:

```bash
$BROWSER https://williamroyce.github.io/torsion-gertsenshtein/
```

Troubleshooting

- FileNotFoundError: ffmpeg — install ffmpeg (see apt command above) or let the example produce a GIF.
- Type-checker warnings about third-party stubs — run examples anyway; code uses TYPE_CHECKING guards and runtime-safe casts where necessary.

## Known issues & design tradeoffs

- py-pde typing stubs are incomplete in places. The code uses guarded imports under `TYPE_CHECKING` and casts to keep type checkers happy while remaining safe at runtime.
- The `numba` backend requires a specialized RHS implementation. If the PDE does not provide `_make_pde_rhs_numba`, `runners.run` will fall back to `numpy`.
- Boundary conditions: py-pde expects structured BC descriptors (dicts); helpers attempt to infer a suitable BC, but explicit BCs are more robust.
- Observers and trackers: the runner exposes `snapshot_interval` to let callers request less-frequent or more-frequent snapshotting without modifying observer code.

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

No license file included in this repository. Add a LICENSE file if you intend to make this open source.
