# Lagrangian-to-PDE Pipeline Examples

This directory contains examples demonstrating the **complete symbolic-to-numerical pipeline** for deriving and simulating field equations from Lagrangians.

## Pipeline Overview

The pipeline has two stages:

### Stage 1: Symbolic Derivation (Mathematica/xAct)
Lagrangian → Euler-Lagrange equations → Component PDEs → JSON export

### Stage 2: Numerical Simulation (Python/py-pde)
JSON → Dynamic PDE construction → Simulation → Visualization

**Key Point**: No physics is hardcoded in Python. All equation structure comes from the JSON files that were symbolically derived from Lagrangians.

---

## Available Examples

### 1. Electromagnetic Field (Vector Field)

**Lagrangian**: `L = -1/4 F_μν F^μν`

**Stage 1 - Derivation**:
```bash
cd examples/electromagnetic
wolframscript -file em_lagrangian_1d.wls
```
- Derives Maxwell equations from Lagrangian
- Applies Lorenz gauge: `∂_μ A^μ = 0`
- Decomposes to wave equations: `∂²A_μ/∂t² = ∂²A_μ/∂x²`
- Exports to: `examples/data/em_1d.json`

**Stage 2 - Simulation**:
```bash
python examples/electromagnetic/em_from_lagrangian.py
```
- Loads equation specification from JSON
- Builds PDE dynamically (no hardcoded wave equation)
- Simulates Gaussian pulse in A₁ component
- Output: `outputs/em_from_lagrangian_output.png`

**Physics**: Massless EM waves propagating at c = 1, clean wave splitting

---

### 2. Klein-Gordon Field (Scalar Field)

**Lagrangian**: `L = -1/2 (∂φ)² - 1/2 m²φ²`

**Stage 1 - Derivation**:
```bash
cd examples/scalar_field
wolframscript -file klein_gordon.wls
```
- Derives Klein-Gordon equation from Lagrangian
- Decomposes to: `∂²φ/∂t² = ∂²φ/∂x² - m²φ` (with m² = 1)
- Exports to: `examples/data/klein_gordon_1d.json`

**Stage 2 - Simulation**:
```bash
python examples/scalar_field/kg_from_lagrangian.py
```
- Loads equation specification from JSON
- Builds PDE with both Laplacian and mass terms
- Simulates Gaussian pulse with dispersion
- Output: `outputs/kg_from_lagrangian_output.png`

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
cd /workspaces/torsion-gertsenshtein/examples/electromagnetic
wolframscript -file em_lagrangian_1d.wls

# Verify JSON was created
cat ../data/em_1d.json | jq '.equations[0].rhs.terms'

# Stage 2: Simulate
python em_from_lagrangian.py

# Check output
ls -lh ../../outputs/em_from_lagrangian_output.png
```

### Klein-Gordon Example
```bash
# Stage 1: Derive from Lagrangian
cd /workspaces/torsion-gertsenshtein/examples/scalar_field
wolframscript -file klein_gordon.wls

# Verify JSON was created
cat ../data/klein_gordon_1d.json | jq '.equations[0].rhs.terms'

# Stage 2: Simulate
python kg_from_lagrangian.py

# Check output
ls -lh ../../outputs/kg_from_lagrangian_output.png
```

---

## Comparison: EM vs Klein-Gordon

| Aspect | EM (Vector) | Klein-Gordon (Scalar) |
|--------|-------------|----------------------|
| **Lagrangian** | `-1/4 F_μν F^μν` | `-1/2 (∂φ)² - 1/2 m²φ²` |
| **Field Components** | 2 (A₀, A₁) | 1 (φ) |
| **Operators in JSON** | Laplacian only | Laplacian + Identity (mass) |
| **Wave Speed** | c = 1 (massless) | Dispersive (massive) |
| **Amplitude** | Conserved | Decreases (spreading) |
| **Gauge** | Lorenz gauge | N/A |

Both examples demonstrate that **different Lagrangians produce different field equations**, and the Python simulation correctly implements the symbolically-derived physics **without manual intervention**.

---

## Validation

Run the full validation suite:
```bash
cd /workspaces/torsion-gertsenshtein
python validate_implementation.py
```

This verifies:
- JSON loading and parsing
- Dynamic PDE construction from specifications
- Simulation produces physically correct results
- Operators (laplacian, identity) correctly identified

---

## Creating New Examples

To add a new field theory:

1. **Create Mathematica script** (e.g., `examples/new_field/derive_equations.wls`):
   - Define manifold and metric
   - Construct Lagrangian using xAct
   - Call `EulerLagrangeEquation`
   - Call `DecomposeToComponents`
   - Call `BuildJSONStructure` and export

2. **Create Python simulation** (e.g., `examples/new_field/simulate.py`):
   - Load JSON: `spec = load_equation_system("../data/new_field.json")`
   - Build PDE: `pde = build_pde_from_json("../data/new_field.json")`
   - Create initial conditions
   - Run simulation
   - Visualize results

The pipeline handles the rest automatically!

---

## Directory Structure

```
examples/
├── README.md                      # This file
├── data/                          # Generated JSON specifications
│   ├── em_1d.json                # EM equations (from Stage 1)
│   └── klein_gordon_1d.json      # KG equations (from Stage 1)
├── electromagnetic/               # EM field examples
│   ├── em_lagrangian_1d.wls     # Stage 1: Derive from Lagrangian
│   └── em_from_lagrangian.py     # Stage 2: Simulate
├── scalar_field/                  # Scalar field examples
│   ├── klein_gordon.wls          # Stage 1: Derive from Lagrangian
│   └── kg_from_lagrangian.py     # Stage 2: Simulate
└── outputs/                       # Visualization outputs
    ├── em_from_lagrangian_output.png
    └── kg_from_lagrangian_output.png
```

---

## Dependencies

**Mathematica/xAct** (Stage 1):
- Wolfram Engine 14.3+
- xAct 1.3.0 (xTensor, xCoba, xPert)

**Python** (Stage 2):
- py-pde ≥ 0.38
- numpy, matplotlib
- torsion_gertsenshtein package

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
