# Examples

TIDAL includes 20 pipeline examples spanning 1+1D through 3+1D spacetimes with scalar, vector, and tensor fields.

## Pipeline Examples

All pipeline examples follow the same pattern: a Wolfram script (`.wls`) derives equations from a Lagrangian, exports to JSON, and a Python script or CLI command runs the simulation.

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

### With Python scripts

Examples with a `simulation.py` or `kg_from_lagrangian.py` can be run directly:

```bash
uv run python examples/scalar_field/kg_from_lagrangian.py
uv run python examples/scalar_vector_coupling/simulation.py
```

## Example Completeness

| Example | theory.toml | run.sh | simulation.py |
|---------|:-----------:|:------:|:-------------:|
| scalar_field | Yes | Yes | Yes |
| electromagnetic | Yes | Yes | Yes |
| proca | Yes | Yes | Yes |
| coupled_scalars | Yes | Yes | Yes |
| chern_simons | Yes | Yes | Yes |
| elasticity | Yes | Yes | Yes |
| curved_spacetime | Yes | Yes | Yes |
| sphere_kg | Yes | Yes | — |
| polar_kg | Yes | Yes | — |
| electrostatics | Yes | Yes | — |
| scalar_vector_coupling | Yes | Yes | Yes |
| scalar_field_3d | Yes | Yes | — |
| spherical_kg | Yes | Yes | — |
| cylindrical_kg | Yes | Yes | — |
| massive_gravity | Yes | Yes | Yes |
| gravitational_waves | Yes | Yes | — |
| massive_3form | Yes | Yes | Yes |
| coupled_proca | Yes | Yes | Yes |
| coupled_scattering | Yes | Yes | Yes |
| scalar_potential_well | Yes | Yes | Yes |

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
