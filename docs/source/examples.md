# Examples

TIDAL includes 27 pipeline examples spanning 1+1D through 3+1D spacetimes with scalar, vector, and tensor fields.

## Pipeline Examples

All pipeline examples follow the same pattern: `tidal derive theory.toml` derives equations from a Lagrangian (via wolframscript), exports to JSON, and a Python script or CLI command runs the simulation.

| Example | Dim | Fields | Key Features |
|---------|-----|--------|--------------|
| `scalar_field/` | 1+1D | phi_0 | Basic KG equation, mass term, dispersion |
| `electromagnetic/` | 1+1D | A_0, A_1 | Maxwell equations, Lorenz gauge |
| `proca/` | 1+1D | A_0, A_1 | Massive vector field (Proca mass term) |
| `coupled_scalars/` | 1+1D | phi_0, chi_0 | Cross-field coupling, mass matrix, parameter sweeps |
| `chern_simons/` | 2+1D | A_0, A_1, A_2 | Topological term, epsilon tensor, A_0 constraint |
| `elasticity/` | 2+1D | u_0, u_1 | Anisotropic Laplacian, cross_derivative_xy |
| `curved_spacetime/` | 2+1D | phi_0 | De Sitter Hubble friction, time-dependent coefficients |
| `sphere_kg/` | 2+1D | phi_0 | Stereographic S2, position-dependent coefficients |
| `polar_kg/` | 2+1D | phi_0 | Polar coordinates, 1/r Christoffel correction |
| `electrostatics/` | 2+1D | phi_0 | Poisson equation, constraint solve |
| `scalar_vector_coupling/` | 2+1D | phi_0, A_0, A_1, A_2 | Mixed-rank cross-field (scalar+vector), 4 constants |
| `scalar_field_3d/` | 3+1D | phi_0 | Full 4D Klein-Gordon |
| `spherical_kg/` | 3+1D | phi_0 | Spherical coordinates, trig coefficients (Cot, Csc) |
| `cylindrical_kg/` | 3+1D | phi_0 | Cylindrical coordinates, mixed curved/flat |
| `massive_gravity/` | 2+1D | h_ij | Linearized massive gravity, m² mass term, xPert linearization |
| `gravitational_waves/` | 3+1D | h_ij | xPert linearization, TT gauge, rank-2 tensor, constraints |
| `massive_3form/` | 3+1D | C_ijk | Rank-3 antisymmetric tensor, symmetry reduction 64 to 4, KG per component |
| `coupled_proca/` | 2+1D | A_i, B_i | Two massive vectors, coupled Helmholtz constraints, periodic BCs |
| `coupled_scattering/` | 2+1D | phi_0, chi_0 | Position-dependent Gaussian coupling, background fields, wave scattering |
| `scalar_potential_well/` | 1+1D | phi_0 | Background potential well, `[[background_fields]]`, bound states |
| `cylindrical_kg_1d/` | 1+1D | phi_0 | Cylindrical coordinates, plane-wave 1D reduction |
| `gravitational_waves_1d/` | 1+1D | h_ij | Linearized gravity, plane-wave 1D reduction |
| `spherical_kg_1d/` | 1+1D | phi_0 | Spherical coordinates, plane-wave 1D reduction |
| `proca_background/` | 2+1D | A_i, B_i | Lorentzian scalar background, two Proca vectors, constraint+BG integration |
| `vector_background/` | 2+1D | phi_0, A_i | Tanh domain wall vector background, ComponentValue mechanism, sign-changing coupling |

## Running Examples

### With the CLI

```bash
# List available JSON specifications
tidal list

# Run a simulation
tidal simulate examples/data/klein_gordon_1d.json --t-end 20 --ic gaussian

# Derive from TOML (requires wolframscript)
tidal derive examples/scalar_field/theory.toml --run
```

### With run.sh scripts

Most examples include a `run.sh` script that runs the full pipeline:

```bash
cd examples/scalar_field
bash run.sh
```

### Running examples

Each example directory has a `run.sh` that runs the full pipeline:

```bash
cd examples/scalar_field && bash run.sh
cd examples/scalar_vector_coupling && bash run.sh
```

## Example Completeness

| Example | theory.toml | run.sh |
|---------|:-----------:|:------:|
| scalar_field | Yes | Yes |
| electromagnetic | Yes | Yes |
| proca | Yes | Yes |
| coupled_scalars | Yes | Yes |
| chern_simons | Yes | Yes |
| elasticity | Yes | Yes |
| curved_spacetime | Yes | Yes |
| sphere_kg | Yes | Yes |
| polar_kg | Yes | Yes |
| electrostatics | Yes | Yes |
| scalar_vector_coupling | Yes | Yes |
| scalar_field_3d | Yes | Yes |
| spherical_kg | Yes | Yes |
| cylindrical_kg | Yes | Yes |
| massive_gravity | Yes | Yes |
| gravitational_waves | Yes | Yes |
| massive_3form | Yes | Yes |
| coupled_proca | Yes | Yes |
| coupled_scattering | Yes | Yes |
| scalar_potential_well | Yes | Yes |
| proca_background | Yes | Yes |
| vector_background | Yes | Yes |
| cylindrical_kg_1d | Yes | Yes |
| gravitational_waves_1d | Yes | Yes |
| spherical_kg_1d | Yes | Yes |

## Parameter Sweeps

TIDAL includes example scripts demonstrating the `tidal sweep` framework for automated parameter studies, convergence analysis, and sensitivity analysis.

### Sweep Examples

| Script | Example | Features Demonstrated |
|--------|---------|----------------------|
| `coupled_scattering/sweep.sh` | 1D coupling sweep | Basic CLI sweep (F1) |
| `coupled_scattering/sweep_2d.sh` | 2D coupling × mass | Cartesian product sweep (F1) |
| `coupled_scattering/convergence.sh` | Grid convergence | `--converge` flag |
| `coupled_scattering/sweep_adaptive.sh` | Adaptive sampling | `--adaptive-metric` (F2a) |
| `coupled_scattering/sweep_advanced.toml` | Full TOML config | Adaptive + resume + parallel (F1+F2a+F7) |
| `coupled_scattering/sweep_scattering.sh` | Asymptotic observables | asymptotic + peak_conversion measurements |
| `coupled_scalars/sweep_coupling.sh` | 1D coupling sweep | Basic CLI sweep (F1) |
| `coupled_scalars/sweep_mass_ratio.sh` | Mass ratio resonance | velocity + resonance measurements (F3) |
| `coupled_scalars/sweep_2d_sensitivity.toml` | Latin Hypercube 2D | LHS sampling (F2b) |
| `coupled_scalars/sweep_2d_sensitivity.sh` | Sensitivity analysis | Sobol/Morris + advanced viz (F4+F8) |
| `scalar_field/sweep_convergence.sh` | Single-field convergence | Conservation convergence study |

### Running Sweep Examples

```bash
# Run from the example directory
cd examples/coupled_scattering && bash sweep_adaptive.sh

# Or use a TOML config directly
tidal sweep --config examples/coupled_scattering/sweep_advanced.toml

# Visualize results
tidal plot sweeps/output/ --type sweep --metric P_max
```

### Measurement Suitability by Example Type

| Category | Examples | Suitable Measurements |
|----------|---------|----------------------|
| Coupled (flat, homogeneous) | `coupled_scalars`, `coupled_proca`, `scalar_vector_coupling` | All 13 types |
| Scattering (flat, background) | `coupled_scattering`, `proca_background`, `vector_background` | energy, conversion, mixing, asymptotic, peak_conversion, conservation |
| Single-field / curved | `scalar_field`, `proca`, `polar_kg`, etc. | energy, conservation, spectrum, dispersion |

### Feature Coverage

| Feature | Description | Example |
|---------|-------------|---------|
| F1 | CLI + TOML sweep configs | `sweep_coupling.toml`, `sweep_advanced.toml` |
| F2a | Adaptive sampling | `sweep_adaptive.sh` |
| F2b | Latin Hypercube sampling | `sweep_2d_sensitivity.toml` |
| F3 | Velocity + resonance analysis | `sweep_mass_ratio.sh` |
| F4 | Sobol/Morris sensitivity | `sweep_2d_sensitivity.sh` |
| F7 | Run status tracking + resume | `sweep_advanced.toml` |
| F8 | Advanced visualization | `sweep_2d_sensitivity.sh` |

## JSON Specification Files

All generated JSON files live in `examples/data/`:

```bash
ls examples/data/*.json
```

Each JSON file is a self-contained equation system specification that can be loaded with:

```python
from tidal.symbolic import build_pde_from_json

pde = build_pde_from_json("examples/data/klein_gordon_1d.json")
```

See the [Pipeline](pipeline.md) page for details on the JSON schema.
