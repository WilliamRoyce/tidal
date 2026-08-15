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
[![codecov](https://codecov.io/gh/WilliamRoyce/torsion-gertsenshtein/branch/main/graph/badge.svg)](https://codecov.io/gh/WilliamRoyce/torsion-gertsenshtein)

</p>

## About

TIDAL derives and integrates the linearized equations of motion of **any tensorial Lagrangian** on a given background. You write the Lagrangian in a TOML file; TIDAL performs the variation symbolically in Mathematica/xAct, decomposes the result into components, exports a JSON specification, and evolves it numerically. The numerical layer contains no physics — every equation it solves was derived from your Lagrangian.

Fields may be scalars, vectors or tensors up to rank 3, in 1+1D through 3+1D, on flat or curved spacetimes, with background fields, gauge fixing, and Poincaré gauge theory (propagating torsion) supported directly in the theory specification.

> Define a Lagrangian in TOML → derive linearized PDEs symbolically → simulate with adaptive solvers → measure conversion, spectra and scattering → sweep or sample over parameter space.

Written by William Royce (`wr286@cantab.ac.uk`), developed in the Astrophysics Group, Cavendish Laboratory, University of Cambridge.

API documentation: <https://williamroyce.github.io/torsion-gertsenshtein/>

## Citing TIDAL

A software paper describing TIDAL is in preparation. Until it appears, please cite this repository:

```text
W. Royce, "TIDAL: Tensor Integration and Derivation for Any Lagrangian",
https://github.com/WilliamRoyce/torsion-gertsenshtein
```

`tidal --cite` prints this along with the citations for SUNDIALS, SciPy, NumPy and xAct, which you should cite alongside TIDAL when you use the corresponding components. Physics work using TIDAL is also in preparation.

## Example

A complete run, from Lagrangian to measured physics. The theory is an effective 1+1D graviton–photon system — two scalars coupled by a gradient term `B0 * h * ∂ₓa` — for which the conversion probability has a known closed form.

**1. Derive the equations of motion** from `examples/coupled_scalars/theory.toml`:

```bash
tidal derive examples/coupled_scalars/theory.toml
```

**2. Inspect what came out of the symbolic pipeline:**

```console
$ tidal inspect examples/data/coupled_scalars.json

Spacetime:
  Dimension: 2 (1+1D)
  Coordinates: ('t', 'x')

Fields (2 components):
  a_0          dynamical    time_order=2
  h_0          dynamical    time_order=2

Equations:
  d2_t(a_0) = [-omegaP2] identity(a_0) [B0] gradient_x(h_0) +1 laplacian_x(a_0)
  d2_t(h_0) = [B0^2] identity(h_0) [mg2/kappa^2] identity(h_0) [B0] gradient_x(a_0) [-kappa^(-2)] laplacian_x(h_0)

Required parameters:
  B0  (in: a_0, h_0)
  kappa  (in: h_0)
  mg2  (in: h_0)
  omegaP2  (in: a_0)
```

**3. Simulate.** The solver is chosen automatically from the structure of the equations:

```console
$ tidal simulate examples/data/coupled_scalars.json \
    --param kappa=1.0 --param B0=0.1 --param omegaP2=0.0 --param mg2=0.0 \
    --grid-shape 256 --bounds 0:100 --periodic \
    --ic plane-wave --ic-component h_0 --ic-wavevector 2.0 --ic-amplitude 0.1 \
    --t-end 50.0 --output run/

  Auto-selected solver: modal
  21 snapshots stored

  a_0: peak 0.0000 → 0.0560
  h_0: peak 0.0924 → 0.0750 (ratio: 0.8115)
```

**4. Measure.** Energy transfer from one field to the other, and a conservation diagnostic:

```console
$ tidal measure run/ --what conversion,conservation --source h_0 --target a_0 \
    --param kappa=1.0 --param B0=0.1 --param omegaP2=0.0 --param mg2=0.0

Energy Conservation: PASS
  max |dE/E|           2.57e-15
  threshold            1e-03

Conversion (h_0 -> a_0):
  Peak P(t)            0.995633
  at t                 32.50
```

The analytic result for this system is `P(t) = sin²(κB₀t/2)`, peaking at `t = π/(κB₀) ≈ 31.4`. Nothing about that formula appears anywhere in the code.

## General use

### Defining a theory

A `theory.toml` specifies the spacetime, the fields, the constants and the Lagrangian:

```toml
[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["m2"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - m2/2 phi[]^2"
```

Optional sections extend this:

- `[[derived_fields]]` — intermediate tensors, e.g. a field strength `F_ab = CD[-a][A[-b]] - CD[-b][A[-a]]`
- `[[background_fields]]` — external fields (magnetic fields, potentials, vectors) giving position-dependent coefficients
- `[[gauge]]` — per-field gauge fixing: Lorenz, de Donder, Coulomb, temporal, axial
- `[torsion]` — Poincaré gauge theory with propagating torsion
- `[parameters]` — runtime defaults, overridable with `--param`

### Solver backends

Five time-integration backends, selected automatically from the structure of the equations:

| Backend      | Library         | Use Case                                    | Key Feature                                   |
| ------------ | --------------- | ------------------------------------------- | --------------------------------------------- |
| **IDA**      | SUNDIALS        | DAE (algebraic constraints)                 | Implicit Newton, 3-tier analytical Jacobian   |
| **CVODE**    | SUNDIALS        | Adaptive ODE (waves)                        | BDF, tolerance control, sparse Jacobian       |
| **Modal**    | numpy/scipy     | Exact spectral (periodic, time-independent) | Machine-precision via eigendecomposition      |
| **Leapfrog** | numpy           | Symplectic                                  | Exact energy conservation (Yoshida 4th-order) |
| **scipy**    | scipy.integrate | General-purpose                             | DOP853, Radau, BDF via `solve_ivp`            |

Spatial operators: 2nd/4th/6th-order finite-difference stencils and FFT spectral operators (auto-enabled for all-periodic BCs). Analytical Jacobians in three tiers — dense, sparse CSC with SuperLU_MT, and GMRES with a Jacobian-vector product for the largest systems. Constraint initial conditions are pre-solved automatically.

### Commands

| Command                     | Description                                                           |
| --------------------------- | --------------------------------------------------------------------- |
| `tidal derive theory.toml`  | Generate .wls from TOML, run wolframscript to produce JSON            |
| `tidal simulate spec.json`  | Run a simulation (`--param`, `--ic`, `--bc`, `--scheme`)              |
| `tidal measure result_dir/` | Extract physics measurements from simulation output                   |
| `tidal sweep spec.json`     | Parameter sweeps, convergence studies, adaptive sampling              |
| `tidal sample spec.json`    | Bayesian inference — Monte Carlo or nested sampling                   |
| `tidal analyze sweep_dir/`  | Post-hoc sensitivity analysis (Sobol/Morris) of sweep results         |
| `tidal plot result_dir/`    | Standalone plotting from simulation output directories                |
| `tidal inspect spec.json`   | Display equation system info (fields, operators, parameters)          |
| `tidal validate spec.json`  | Validate a JSON equation specification                                |
| `tidal list`                | Discover available JSON specs in `examples/data/`                     |
| `tidal doctor`              | Environment diagnostics (Wolfram, dependencies, xAct)                 |

Measurements available through `tidal measure`: energy, conversion `P(t)`, mixing length, spectrum, spectral conversion `P(k,t)`, dispersion `ω(k)`, conservation diagnostics, effective mass, asymptotic scattering, peak conversion, group/phase velocity, and resonance analysis.

`tidal sample` performs Bayesian parameter estimation over a theory's coupling space, with priors, inequality constraints, and either Monte Carlo or nested sampling (dynesty or PolyChord). Install with `pip install tidal[inference]`.

### Worked examples

Each directory under `examples/` contains a `theory.toml`, the generated `.wls`, and a `run.sh` showing the full derive → inspect → simulate → measure workflow.

| Example                    | Dim  | Key Features                                                         |
| -------------------------- | ---- | -------------------------------------------------------------------- |
| `chern_simons/`            | 2+1D | Epsilon tensor, topological mass, `A_0` constraint                   |
| `coupled_scalars/`         | 1+1D | Cross-field gradient coupling, mass matrix, energy transfer          |
| `coupled_scattering/`      | 3+1D | Position-dependent Gaussian coupling, background fields, scattering  |
| `curved_spacetime/`        | 1+1D | de Sitter and conformal static metrics, Hubble friction              |
| `dark_photon_plasma/`      | 3+1D | Torsion dark photon with kinetic mixing and photon plasma mass       |
| `elasticity/`              | 2+1D | Navier–Cauchy, anisotropic laplacian, cross-derivative operators     |
| `euler_heisenberg/`        | 3+1D | Vacuum birefringence, quartic `(F·F)²` QED correction                |
| `gertsenshtein/`           | 3+1D | Einstein–Maxwell graviton–photon conversion, multi-field xPert       |
| `gertsenshtein_proca/`     | 3+1D | Graviton–photon conversion with a photon effective mass              |
| `gravitational_waves/`     | 3+1D | xPert linearization, TT gauge, constraints                           |
| `graviton_torsion/`        | 3+1D | General quadratic PGT, torsion perturbation, graviton–torsion mixing |
| `massive_3form/`           | 3+1D | Rank-3 antisymmetric tensor, symmetry reduction                      |
| `massive_gravity/`         | 2+1D | Fierz–Pauli mass, xPert, coupled constraints                         |
| `proca_background/`        | 2+1D | Two massive vectors, Lorentzian scalar background                    |
| `scalar_vector_coupling/`  | 2+1D | Mixed-rank cross-field coupling, Chern–Simons and divergence terms   |
| `spherical_kg_1d/`         | 3+1D | Spherical coordinates, plane-wave dimensional reduction              |
| `torsion_dark_photon/`     | 3+1D | Propagating torsion, kinetic mixing, non-minimal coupling            |
| `torsion_dark_photon_fv/`  | 3+1D | Massive Proca dark photon with kinetic mixing                        |
| `torsion_gertsenshtein/`   | 3+1D | Graviton–photon conversion with PGT torsion                          |

Run one directly:

```bash
cd examples/coupled_scalars && bash run.sh
```

## Quickstart

TIDAL uses **uv** for Python version, environment and dependency management, and is container-first. Python 3.11 is required (numba/llvmlite compatibility).

```bash
uv python pin 3.11
uv sync --all-extras
uv run python -c "import tidal; print('OK')"
```

Simulation, measurement, plotting, sweeps and inference all work with this alone. Deriving new equations from a Lagrangian additionally requires Wolfram Engine and xAct — see below.

### Dev container (VS Code / Codespaces)

The repository includes a Debian-based dev container, which is the supported path for a consistent toolchain:

- Open the folder in VS Code
- Command Palette → Dev Containers: Reopen in Container
- Inside the container: `uv python pin 3.11 && uv sync --all-extras`

### Optional: video output

For MP4 output via Matplotlib's `FFMpegWriter`:

```bash
sudo apt-get update && sudo apt-get install -y ffmpeg
```

Without `ffmpeg`, animations fall back to GIF via Pillow.

## Symbolic computing setup

Required only for `tidal derive` — deriving linearized field equations from a Lagrangian. Everything downstream of a JSON specification runs without it.

```bash
# 1. Download the Wolfram Engine installer from https://www.wolfram.com/engine/
#    and place it in third_party/

# 2. Install and activate Wolfram Engine
sudo ./scripts/install-wolfram-engine.sh
./scripts/activate-wolfram.sh

# 3. Install the xAct/xCoba tensor algebra packages
./scripts/install-xact-xcoba.sh

# 4. Verify the complete setup
./scripts/verify-wolfram-setup.sh
```

The verification script checks Wolfram Engine activation, xAct package installation (xCore, xPerm, xTensor, xCoba, xPert), xPerm binary GLIBC compatibility, and runs a smoke test with tensor operations. `tidal doctor` performs the same diagnosis at any time.

Note that a Wolfram Engine license permits **one** `wolframscript` session at a time; do not run `tidal derive` in parallel.

## Logging and profiling

The library uses Python's `logging` module and never prints from library code. To see solver progress and profiling summaries:

```python
import logging
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
```

Pass `profile=True` to a `run` call for a timing breakdown.

## Troubleshooting

- **Import errors in VS Code** (e.g. numpy not found): ensure the interpreter is the repository venv (`.venv/bin/python3`), then reload the window.
- **`llvmlite`/`numba` build failures**: use Python 3.11 (`uv python pin 3.11`).
- **`FileNotFoundError: ffmpeg`**: install ffmpeg as above, or let the example produce a GIF.
- **Type-checker warnings about third-party stubs**: run anyway; the code uses `TYPE_CHECKING` guards and runtime-safe casts where necessary.
- **GitHub Pages 404 or deploy errors**: Settings → Pages → Source must be GitHub Actions, and Actions → Workflow permissions must be Read/Write.
- **Animation has too few frames**: match `snapshot_interval` in the `run` call to your desired temporal resolution, and increase `fps` in `choose_writer_and_out`.
- **Logging messages not visible**: call `logging.basicConfig(level=logging.INFO, ...)` at the start of your script.
- **Wolfram Engine not activated**: run `./scripts/activate-wolfram.sh` and enter your Wolfram ID (a free account suffices).
- **xPerm GLIBC errors** (`GLIBC_2.38 not found`): run `./scripts/install-xact-xcoba.sh` to recompile the binary for your system.
- **xAct packages not loading**: xAct must be installed in `~/.WolframEngine/Applications/xAct/`; run the verification script to diagnose.
- **Anything else**: `tidal doctor` checks Wolfram, Python dependencies and xAct in one pass.

## Origins

TIDAL was written for a Cambridge Part III / MSci project in the Astrophysics Group, Cavendish Laboratory, investigating graviton–photon conversion in Poincaré gauge theory. The project was assessed at 97.03% and nominated by the examiners for both the Theory and Computing prizes.

> "A truly remarkable project — hugely ambitious and brilliantly executed."
>
> "Both the supervisor and the assessor have been PtIII project assessors for many years and have never before encountered such an outstanding project."
>
> — Cambridge Part III examiners' joint report, June 2026

Development continues, towards a distributable package and the associated publications.

## Getting help

- **Questions and ideas**: [GitHub Discussions](https://github.com/WilliamRoyce/torsion-gertsenshtein/discussions)
- **Bugs**: [Issue tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues), with the `bug` label
- **Feature requests**: [Issue tracker](https://github.com/WilliamRoyce/torsion-gertsenshtein/issues), with the `enhancement` label
- **API reference**: <https://williamroyce.github.io/torsion-gertsenshtein/>

## Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, coding standards and the PR checklist. In short: run `./scripts/full_test.sh` before submitting, keep the test suite green, and add tests in both the Python and Wolfram layers where applicable.

## License

MIT. See [LICENSE](LICENSE).

## Acknowledgements

This project builds on:

- [SUNDIALS](https://computing.llnl.gov/projects/sundials) — IDA (DAE) and CVODE (BDF) solvers via [scikit-sundae](https://github.com/NREL/scikit-sundae) (Hindmarsh et al. 2005).
- The [xAct/xTensor ecosystem](http://www.xact.es/) — symbolic tensor algebra (Martín-García et al.) powering the Lagrangian-to-PDE derivation pipeline.
- [xPert](https://www.researchgate.net/publication/1740524) — metric perturbation theory (Brizuela et al. 2009) for linearization.
- [PolyChord](https://github.com/PolyChord/PolyChordLite) (Handley, Hobson & Lasenby 2015) and [dynesty](https://github.com/joshspeagle/dynesty) (Speagle 2020) — nested sampling; [anesthetic](https://github.com/handley-lab/anesthetic) (Handley 2019) for posterior analysis.
- [`uv`](https://github.com/astral-sh/uv) — fast Python environment management.
- Originally built on [py-pde](https://py-pde.readthedocs.io/) (Zwicker, JOSS 2020); finite-difference stencil conventions are retained in TIDAL's native operators.

Design decisions are informed by [Dedalus](https://arxiv.org/abs/1905.10388) (Burns et al. 2020), [MEEP](https://meep.readthedocs.io/) (Oskooi et al. 2010) and [FEniCS](https://fenicsproject.org/) (Baratta et al. 2023). The physics targets the Gertsenshtein effect (Gertsenshtein 1962; [Domcke & Garcia-Cely 2023](https://arxiv.org/abs/2301.02072)).
