# Lagrangian-to-PDE Pipeline Examples

This directory contains examples demonstrating the **complete symbolic-to-numerical pipeline** for deriving and simulating field equations from Lagrangians.

## Pipeline Overview

The pipeline has two stages:

### Stage 1: Symbolic Derivation (Mathematica/xAct)

Lagrangian → Euler-Lagrange equations → Component PDEs → JSON export

### Stage 2: Numerical Simulation (Python/py-pde)

JSON → Dynamic PDE construction → Simulation → Visualization

**Key Point**: No physics is hardcoded in Python. All equation structure comes from the JSON files that were symbolically derived from Lagrangians.

### CLI Workflow (Recommended)

The `tidal` command provides a unified interface for both stages:

```bash
# Derive equations from a TOML configuration
tidal derive examples/scalar_field/theory.toml

# Inspect the resulting equation system
tidal inspect examples/data/klein_gordon_1d.json

# Run simulation with parameter override
tidal simulate examples/data/klein_gordon_1d.json --param m2=1.0 --ic gaussian --t-end 20

# Discover all available JSON specs
tidal list

# Validate a JSON spec
tidal validate examples/data/klein_gordon_1d.json
```

Each example also has a `run.sh` script demonstrating the equivalent CLI workflow.

---

## Available Examples

### 1. Electromagnetic Field (Vector Field)

**Lagrangian**: `L = -1/4 F_μν F^μν`

**Stage 1 - Derivation**:

```bash
tidal derive examples/electromagnetic/theory.toml
```

- Derives Maxwell equations from Lagrangian
- Applies Lorenz gauge: `∂_μ A^μ = 0`
- Decomposes to wave equations: `∂²A_μ/∂t² = ∂²A_μ/∂x²`
- Exports to: `examples/data/em_1d.json`

**Stage 2 - Simulation**:

```bash
cd examples/electromagnetic && bash run.sh
```

- Loads equation specification from JSON
- Builds PDE dynamically (no hardcoded wave equation)
- Simulates Gaussian pulse in A₁ component
- Outputs: heatmap and amplitude plots in `examples/data/`

**Physics**: Massless EM waves propagating at c = 1, clean wave splitting

---

### 2. Klein-Gordon Field (Scalar Field)

**Lagrangian**: `L = -1/2 (∂φ)² - 1/2 m²φ²`

**Stage 1 - Derivation**:

```bash
tidal derive examples/scalar_field/theory.toml
```

- Derives Klein-Gordon equation from Lagrangian
- Decomposes to: `∂²φ/∂t² = ∂²φ/∂x² - m²φ` (with m² = 1)
- Exports to: `examples/data/klein_gordon_1d.json`

**Stage 2 - Simulation**:

```bash
cd examples/scalar_field && bash run.sh
```

- Loads equation specification from JSON
- Builds PDE with both Laplacian and mass terms
- Simulates Gaussian pulse with dispersion
- Outputs: heatmap, amplitude, snapshot, and profile plots in `examples/data/`

**Physics**: Massive scalar field with dispersion (amplitude decreases due to spreading)

---

## Verifying No Hardcoded Physics

The Python simulation scripts use `build_pde_from_json()` which:

1. **Loads** the JSON specification
2. **Parses** operator types (laplacian, identity, etc.)
3. **Builds** PDE class dynamically from specification
4. **Applies** operators based on JSON, not hardcoded logic

### Key Code Structure

```python
# Load equation specification from JSON
spec = load_equation_system("examples/data/em_1d.json")

# Build PDE dynamically (NO hardcoded equations)
pde = build_pde_from_json("examples/data/em_1d.json")

# The PDE's evolution_rate method uses:
#   for term in eq.rhs_terms:  # From JSON spec
#       operator = self.operators[term.operator]  # Dynamic dispatch
#       result += term.coefficient * operator(field)
```

### Proof of Dynamic Construction

Compare the two examples:

- **EM**: JSON has `{"operator": "laplacian", "coefficient": 1.0}` only
- **KG**: JSON has `{"operator": "laplacian", "coefficient": 1.0}` AND `{"operator": "identity", "coefficient": -1.0}`

The Python code doesn't "know" about mass terms - it reads them from JSON and applies them.

---

## Running Full Pipeline End-to-End

### EM Example

```bash
# Stage 1: Derive from Lagrangian
tidal derive examples/electromagnetic/theory.toml

# Verify JSON was created
cat examples/data/em_1d.json | jq '.equations[0].rhs.terms'

# Stage 2: Simulate and plot
cd examples/electromagnetic && bash run.sh
```

### Klein-Gordon Example

```bash
# Stage 1: Derive from Lagrangian
tidal derive examples/scalar_field/theory.toml

# Verify JSON was created
cat examples/data/klein_gordon_1d.json | jq '.equations[0].rhs.terms'

# Stage 2: Simulate and plot
cd examples/scalar_field && bash run.sh
```

---

## Comparison: EM vs Klein-Gordon

| Aspect                | EM (Vector)      | Klein-Gordon (Scalar)       |
| --------------------- | ---------------- | --------------------------- |
| **Lagrangian**        | `-1/4 F_μν F^μν` | `-1/2 (∂φ)² - 1/2 m²φ²`     |
| **Field Components**  | 2 (A₀, A₁)       | 1 (φ)                       |
| **Operators in JSON** | Laplacian only   | Laplacian + Identity (mass) |
| **Wave Speed**        | c = 1 (massless) | Dispersive (massive)        |
| **Amplitude**         | Conserved        | Decreases (spreading)       |
| **Gauge**             | Lorenz gauge     | N/A                         |

Both examples demonstrate that **different Lagrangians produce different field equations**, and the Python simulation correctly implements the symbolically-derived physics **without manual intervention**.

---

### 3. Curvilinear Coordinate Examples (Polar, Spherical, Cylindrical)

These examples demonstrate that the pipeline handles non-Cartesian coordinate systems purely through the metric definition. Coordinate names remain generic (`x`, `y`, `z`) while the metric encodes the geometry. Christoffel corrections are auto-detected from metric type: constant metrics have Γ=0, non-constant metrics (position/time-dependent) trigger explicit Christoffel computation via the standard formula.

#### Polar Coordinates (2+1D)

**Metric**: `ds² = -dt² + dx² + x² dy²` (x=r, y=θ)

```bash
# Stage 1: Derive from Lagrangian
tidal derive examples/polar_kg/theory.toml

# Stage 2: Simulate
cd examples/polar_kg && bash run.sh
```

**Key features**: `gradient_x` with 1/r coefficient, `laplacian_y` with 1/r², mixed periodic BCs.

#### Spherical Coordinates (3+1D)

**Metric**: `ds² = -dt² + dx² + x² dy² + x² sin²(y) dz²` (x=r, y=θ, z=φ)

```bash
tidal derive examples/spherical_kg/theory.toml
cd examples/spherical_kg && bash run.sh
```

**Key features**: 6 RHS terms, trigonometric coefficients (`Cot`, `Csc`), `gradient_y` with cot(θ)/r².

#### Cylindrical Coordinates (3+1D)

**Metric**: `ds² = -dt² + dx² + x² dy² + dz²` (x=r, y=θ, z=z)

```bash
tidal derive examples/cylindrical_kg/theory.toml
cd examples/cylindrical_kg && bash run.sh
```

**Key features**: Mixed curved (r, θ) and flat (z) spatial directions, `laplacian_z` with constant coefficient.

---

### 4. Klein-Gordon 3+1D (Full 4D Spacetime)

**Lagrangian**: `L = -1/2 η^{ab} ∂_a φ ∂_b φ - 1/2 m² φ²`

**Spacetime**: 3+1D flat Minkowski, signature (-,+,+,+), coordinates (t,x,y,z)

**Stage 1 - Derivation**:

```bash
tidal derive examples/scalar_field_3d/theory.toml
```

- Derives Klein-Gordon equation in full 4D spacetime
- Decomposes to: `∂²φ/∂t² = ∂²φ/∂x² + ∂²φ/∂y² + ∂²φ/∂z² - m²φ`
- Exports to: `examples/data/klein_gordon_3d.json`

**Stage 2 - Simulation**:

```bash
cd examples/scalar_field_3d && bash run.sh
```

- Loads equation specification from JSON
- 32³ = 32,768 cell 3D grid with periodic boundary conditions
- 3D Gaussian pulse at rest (momentum = 0)
- Runge-Kutta (RK4) time integration
- Output: 4-panel visualization (z-profile, xy-slices, amplitude)

**Visualization**:

- Panel 1: φ(z) profile at x=y=center (initial vs final)
- Panel 2: x-y slice of initial field at z=center
- Panel 3: x-y slice of final field at z=center
- Panel 4: max|φ| over time (amplitude decay curve)

**Physics**: Massive scalar field in 3D space exhibits spherical wave propagation with amplitude decay ~ 1/r. Initial peak 0.964 decays to 0.362 (ratio ~0.375) as pulse spreads outward.

**Key features**: Full 3+1D pipeline demonstration, directional laplacians (`laplacian_x`, `laplacian_y`, `laplacian_z`), 3D spatial visualization.

---

### 5. Proca Field (Massive Vector, 1+1D)

**Lagrangian**: `L = -1/4 F_ab F^ab - 1/2 m² A_a A^a`

```bash
tidal derive examples/proca/theory.toml
cd examples/proca && bash run.sh
```

**Key features**: Massive vector field, Proca mass term, uses `[[derived_fields]]` for field strength tensor F_ab.

---

### 6. Coupled Scalar Fields (1+1D)

**Lagrangian**: `L = 1/2(∂φ)² - 1/2 m_φ² φ² + 1/2(∂χ)² - 1/2 m_χ² χ² - g φ χ`

```bash
tidal derive examples/coupled_scalars/theory.toml
cd examples/coupled_scalars && bash run.sh
```

**Key features**: Cross-field coupling via `identity` operator on other field, mass matrix, mode-mixing, energy transfer between fields.

---

### 7. Chern-Simons Gauge Theory (2+1D)

**Lagrangian**: `L = -1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a ∂_b A_c`

```bash
tidal derive examples/chern_simons/theory.toml
cd examples/chern_simons && bash run.sh
```

**Key features**: Epsilon tensor (automated), topological mass, A_0 constraint equation (time_derivative_order=0), cross-field gradient coupling.

---

### 8. Elasticity (Navier-Cauchy, 2+1D)

**Lagrangian**: Anisotropic elastic medium

```bash
cd examples/elasticity && bash run.sh
```

**Key features**: Anisotropic laplacian (`laplacian_x`, `laplacian_y` with different coefficients), `cross_derivative_xy` operator.

---

### 9. Curved Spacetime (De Sitter, 2+1D)

**Lagrangian**: KG scalar on de Sitter background with conformal time

```bash
cd examples/curved_spacetime && bash run.sh
```

**Key features**: Time-dependent coefficients (Hubble friction `exp(2Ht)`), Christoffel auto-detection, `first_derivative_t` operator.

---

### 10. Klein-Gordon on 2-Sphere (2+1D)

**Lagrangian**: KG scalar on S² (stereographic projection)

```bash
tidal derive examples/sphere_kg/theory.toml
cd examples/sphere_kg && bash run.sh
```

**Key features**: Position-dependent coefficients (stereographic metric), `_resolve_coefficient_at_point` evaluator.

---

### 11. Electrostatics (Poisson Equation, 2+1D)

**Lagrangian**: Electrostatic potential energy functional

```bash
cd examples/electrostatics && bash run.sh
```

**Key features**: Constraint solver (time_derivative_order=0), `--mode constraint` CLI flag.

---

### 12. Gravitational Waves (Linearized Gravity, 3+1D)

**Lagrangian**: Einstein-Hilbert linearized around Minkowski

```bash
cd examples/gravitational_waves && bash run.sh
```

**Key features**: xPert linearization, TT gauge, constraint equations, field-aware LHS detection, 6 independent metric perturbation components.

---

### 13. Massive 3-Form (Rank-3 Tensor, 3+1D)

**Lagrangian**: Rank-3 antisymmetric tensor with mass term

```bash
tidal derive examples/massive_3form/theory.toml
```

**Key features**: Rank-3 antisymmetric tensor, symmetry reduction (64 → 4 independent components), KG equation per component, `DefConstantSymbol` for mass.

---

### 14. Scalar-Vector Coupling (Stress Test, 2+1D)

**Lagrangian**: `L = KG(φ) + Maxwell(A) + Proca(A) + CS(A) + gSV φ ∂_a A^a`

```bash
tidal derive examples/scalar_vector_coupling/theory.toml
cd examples/scalar_vector_coupling && bash run.sh
```

**Key features**: Mixed-rank cross-field coupling (scalar + vector), 4 symbolic constants (phim2, Am2, kCS, gSV), 4x4 mass/coupling matrices, cross-field `first_derivative_t` and `gradient` operators, `[[derived_fields]]` for F_ab, A_0 constraint.

---

### 15. Massive Gravity (Linearized, 2+1D)

**Field Equation**: `G^(1)_ab[h] - m^2 h_ab = 0`

```bash
tidal derive examples/massive_gravity/theory.toml
cd examples/massive_gravity && bash run.sh
```

**Key features**: xPert linearization, Fierz-Pauli mass term, 6 symmetric tensor components (h_tt constraint, h_tx/h_ty first-order, h_xx/h_xy/h_yy evolution), coupled constraint solver with FFT+SVD regularization, dispersion `ω²=k²+m²`.

---

### 16. Coupled Proca (Two Massive Vectors, 2+1D)

**Lagrangian**: `L = -1/4 F^A F^A - mA²/2 A² - 1/4 F^B F^B - mB²/2 B² + g A·B`

```bash
tidal derive examples/coupled_proca/theory.toml
cd examples/coupled_proca && bash run.sh
```

**Key features**: Two massive vector fields with cross-coupling, coupled Helmholtz constraints (A_0, B_0), periodic boundary conditions, FFT constraint solver, 6×6 mass/coupling matrices.

---

### 17. Coupled Scattering (Background Fields, 2+1D)

**Lagrangian**: `L = KG(φ) + KG(χ) + G(x,y) φ χ` where G is a Gaussian background coupling

```bash
tidal derive examples/coupled_scattering/theory.toml
cd examples/coupled_scattering && bash run.sh
```

**Key features**: Position-dependent Gaussian coupling via `[[background_fields]]`, wave packet scattering, conversion probability measurement, `_eval_utils.py` coefficient evaluation.

---

### 18. Scalar Potential Well (Background Fields, 1+1D)

**Lagrangian**: `L = -1/2 (∂φ)² - 1/2 V(x) φ²` where V(x) is a background potential

```bash
tidal derive examples/scalar_potential_well/theory.toml
cd examples/scalar_potential_well && bash run.sh
```

**Key features**: Non-dynamical background scalar field via `[[background_fields]]`, position-dependent mass term, Gaussian potential well profile, bound state dynamics.

---

## Example Completeness

All examples use `theory.toml` + `run.sh` for the full derive → simulate → plot pipeline. No standalone Python scripts are needed.

| Example                | `theory.toml` | `run.sh` |
| ---------------------- | :-----------: | :------: |
| scalar_field           |       Y       |    Y     |
| electromagnetic        |       Y       |    Y     |
| proca                  |       Y       |    Y     |
| coupled_scalars        |       Y       |    Y     |
| chern_simons           |       Y       |    Y     |
| elasticity             |       Y       |    Y     |
| curved_spacetime       |      \*       |    Y     |
| sphere_kg              |       Y       |    Y     |
| polar_kg               |       Y       |    Y     |
| electrostatics         |       Y       |    Y     |
| scalar_vector_coupling |       Y       |    Y     |
| scalar_field_3d        |       Y       |    Y     |
| spherical_kg           |       Y       |    Y     |
| cylindrical_kg         |       Y       |    Y     |
| gravitational_waves    |       Y       |    Y     |
| massive_3form          |       Y       |    Y     |
| massive_gravity        |       Y       |    Y     |
| coupled_proca          |       Y       |    Y     |
| coupled_scattering     |       Y       |    Y     |
| scalar_potential_well  |       Y       |    Y     |
| proca_background       |       Y       |    Y     |
| vector_background      |       Y       |    Y     |

_Note: `curved_spacetime` uses two separate TOML files (`de_sitter.toml`, `conformal_static.toml`) instead of a single `theory.toml`._

---

## Validation

Run the full validation suite:

```bash
cd /workspaces/torsion-gertsenshtein

# Validate all JSON specs with the CLI
tidal validate examples/data/klein_gordon_1d.json

```

This verifies:

- JSON loading and parsing
- Dynamic PDE construction from specifications
- Simulation produces physically correct results
- Operators (laplacian, identity, gradient, cross_derivative, etc.) correctly identified

---

## Creating New Examples

### Option A: TOML Configuration (Recommended)

Create a `theory.toml` in `examples/new_field/`:

```toml
[theory]
name = "My New Theory"

[spacetime]
dimension = 2  # 1+1D
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["m2"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - m2/2 phi[]^2"

[parameters]
m2 = 1.0

[output]
path = "../data/my_field.json"
```

Then derive and simulate:

```bash
tidal derive examples/new_field/theory.toml
tidal simulate examples/data/my_field.json --param m2=1.0 --ic gaussian
```

For intermediate tensors (e.g., field strength), use `[[derived_fields]]`:

```toml
[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"
```

### Option B: Manual Wolfram Script

For cases requiring custom Wolfram logic (gauge fixing, xPert linearization):

1. **Create Mathematica script** (e.g., `examples/new_field/derive_equations.wls`):
   - Define manifold and metric
   - Construct Lagrangian using xAct
   - Call `EulerLagrangeEquation`
   - Call `DecomposeToComponents`
   - Call `BuildMultiFieldJSONStructure` and export

2. **Create `run.sh`** to simulate and plot:

   ```bash
   tidal simulate ../data/new_field.json --param m2=1.0 --ic gaussian --output ../data/new_field_output
   tidal plot ../data/new_field_output --type heatmap --output ../data/new_field_heatmap.png --quiet
   tidal plot ../data/new_field_output --type amplitude --output ../data/new_field_amplitude.png --quiet
   ```

The pipeline handles the rest automatically!

---

## Directory Structure

```
examples/
├── README.md                      # This file
├── data/                          # Generated JSON specifications (20 files)
│   ├── em_1d.json                # EM equations
│   ├── klein_gordon_1d.json      # KG equations
│   ├── klein_gordon_3d.json      # KG equations (3+1D)
│   ├── proca_1d.json             # Proca (massive vector)
│   ├── coupled_scalars.json      # Coupled scalar fields
│   ├── chern_simons_3d.json      # Chern-Simons (2+1D)
│   ├── navier_cauchy_2d.json     # Elasticity (2+1D)
│   ├── de_sitter_kg.json         # Curved spacetime KG
│   ├── sphere_kg.json            # KG on 2-sphere
│   ├── polar_kg.json             # KG in polar coordinates
│   ├── electrostatics_laplace.json # Electrostatics wave equation
│   ├── scalar_vector_coupling.json # Scalar-vector stress test
│   ├── spherical_kg.json         # KG in spherical coordinates
│   ├── cylindrical_kg.json       # KG in cylindrical coordinates
│   ├── linearized_gravity.json   # Gravitational waves
│   ├── massive_3form.json        # Rank-3 antisymmetric tensor
│   ├── conformal_kg_static.json  # Conformal scalar field
│   ├── massive_gravity_3d.json   # Linearized massive gravity (2+1D)
│   └── coupled_proca_3d.json     # Coupled Proca (2+1D)
├── electromagnetic/               # EM field (1+1D)
├── scalar_field/                  # Scalar field (1+1D)
├── proca/                         # Massive vector field (1+1D)
├── coupled_scalars/               # Coupled scalar fields (1+1D)
├── chern_simons/                  # Chern-Simons gauge theory (2+1D)
├── elasticity/                    # Anisotropic elasticity (2+1D)
├── curved_spacetime/              # De Sitter spacetime (2+1D)
├── sphere_kg/                     # KG on 2-sphere, stereographic (2+1D)
├── polar_kg/                      # KG in polar coordinates (2+1D)
├── electrostatics/                # Poisson equation, constraint solver
├── scalar_vector_coupling/        # Scalar+vector cross-field stress test (2+1D)
├── scalar_field_3d/               # Scalar field (3+1D)
├── spherical_kg/                  # KG in spherical coordinates (3+1D)
├── cylindrical_kg/                # KG in cylindrical coordinates (3+1D)
├── gravitational_waves/           # Linearized gravity (3+1D)
├── massive_3form/                 # Rank-3 antisymmetric tensor (3+1D)
├── massive_gravity/               # Linearized massive gravity (2+1D)
├── coupled_proca/                 # Two coupled massive vectors (2+1D)
├── coupled_scattering/            # Position-dependent coupling, background fields (2+1D)
├── proca_background/              # Scalar background, two Proca vectors (2+1D)
├── scalar_potential_well/         # Background potential well, bound states (1+1D)
└── vector_background/             # Vector domain wall background (2+1D)
```

---

## Dependencies

**Mathematica/xAct** (Stage 1):

- Wolfram Engine 14.3+
- xAct 1.3.0 (xTensor, xCoba, xPert)

**Python** (Stage 2):

- py-pde ≥ 0.38
- numpy, matplotlib
- tidal package

---

## Key Insight

The separation into two stages is **intentional and beneficial**:

1. **Stage 1 (symbolic)**: Mathematica/xAct excels at tensor algebra and symbolic manipulation
2. **Stage 2 (numerical)**: Python/py-pde excels at efficient numerical simulation

The JSON format serves as a **well-defined interface** between the symbolic and numerical worlds, allowing each tool to do what it does best.

This design means:

- ✅ No hardcoded physics in Python
- ✅ Equations verified symbolically before simulation
- ✅ Easy to add new field theories (just write new Lagrangian)
- ✅ Separation of concerns (symbolic vs numerical expertise)
