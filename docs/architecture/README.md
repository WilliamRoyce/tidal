# Architecture

Visual guide to the TIDAL Lagrangian-to-PDE pipeline.
All PDEs derive from a Lagrangian via symbolic computation — no physics is hardcoded in Python.

---

## 1. Data Flow Pipeline

The pipeline has two stages separated by a JSON specification file.

**Stage 1 (Mathematica/xAct)** derives field equations symbolically from a Lagrangian.
**Stage 2 (Python/py-pde)** loads the specification and runs a numerical simulation.

```mermaid
flowchart LR
    subgraph stage1["Stage 1 — Symbolic (Mathematica / xAct)"]
        L["Lagrangian<br/>L(φ, ∂φ, g)"] --> EL["EulerLagrange.wl<br/>δL/δφ = 0"]
        EL --> CD["ComponentDecompose.wl<br/>Tensor → scalar components"]
        CD --> EJ["ExportJSON.wl<br/>BuildMultiFieldJSONStructure"]
    end

    EJ --> JSON[("JSON Specification<br/>examples/data/*.json")]

    subgraph stage2["Stage 2 — Numerical (Python / py-pde)"]
        JSON --> JL["json_loader.py<br/>load_equation_system()"]
        JL --> ES["EquationSystem<br/>(frozen dataclass)"]
        ES --> PB["pde_builder.py<br/>PDEFromSpec(PDEBase)"]
        PB --> SOLVE["py-pde solver<br/>(scipy RK45)"]
        SOLVE --> OUT["Simulation results<br/>& plots"]
    end
```

### Wolfram Module Roles

| Module | Purpose | Key Export |
|--------|---------|-----------|
| `CommonUtilities.wl` | Shared helpers: CD→Derivative conversion, Christoffel/epsilon evaluation, metric components | `ConvertCDToDerivatives`, `EvaluateChristoffelComponents` |
| `EulerLagrange.wl` | Derive equations of motion from a Lagrangian | `EulerLagrangeEquation[L, field, cd]` |
| `ComponentDecompose.wl` | Convert tensor EOM to scalar component equations | `DecomposeToComponents[eom, field, chart, additionalFields]` |
| `ExportJSON.wl` | Serialize component equations to JSON | `BuildMultiFieldJSONStructure[fieldEquations, metadata]` |
| `Linearize.wl` | Perturbation theory (xPert) for linearized gravity | `LinearizeTensorExpression[expr]` |

### JSON Specification Structure

```
{
  "spacetime": { dimension, signature, coordinates },
  "fields":    [{ name, index, is_dynamical }],
  "equations": [{
    field_name,
    lhs: { expression, order: { time: N } },
    rhs: [{ coefficient, coefficient_symbolic, operator, field,
             time_dependent, coordinate_dependent }]
  }],
  "coupling":  { mass_matrix, coupling_matrix },
  "metadata":  { lagrangian, gauge, linearized, ... }
}
```

---

## 2. Python Module Dependencies

```mermaid
flowchart TD
    subgraph pkg["tidal"]
        init["__init__.py"]

        subgraph sym["symbolic/"]
            sinit["__init__.py"]
            jl["json_loader.py<br/><i>EquationSystem, OperatorTerm,<br/>ComponentEquation, LHSStructure</i>"]
            pb["pde_builder.py<br/><i>PDEFromSpec, build_pde_from_json,<br/>create_initial_state</i>"]
        end

        subgraph cli["cli/"]
            cmain["__main__.py<br/><i>Entry point</i>"]
            cderiv["_derive.py<br/><i>TOML → .wls → JSON</i>"]
            csim["_simulate.py<br/><i>JSON → PDE → solve</i>"]
        end

        subgraph meas["measurement/"]
            menergy["_energy.py<br/><i>Hamiltonian energy</i>"]
            mconv["_conversion.py<br/><i>Field conversion P(t)</i>"]
        end

        vf["vectorfield/<br/><i>ComponentGaussianPulse,<br/>ComponentPlaneWave</i>"]
        utils["utils.py<br/><i>normalize_solve_result</i>"]
    end

    init --> sinit
    init --> vf
    sinit --> jl
    sinit --> pb
    pb -- "imports dataclasses" --> jl
    cmain --> cderiv
    cmain --> csim
    csim --> pb

    subgraph ext["External Libraries"]
        pypde["py-pde<br/><i>PDEBase, ScalarField,<br/>FieldCollection, CartesianGrid</i>"]
        numpy["numpy / scipy"]
    end

    pb --> pypde
    pb --> numpy
```

The `symbolic/` subpackage is the core pipeline — it loads JSON specs and builds PDE solvers.
The `cli/` subpackage provides the `tidal` command-line interface (derive, simulate, inspect, measure, list, validate).
The `measurement/` subpackage provides post-hoc analysis tools (energy, conversion, spectra, mixing).

---

## 3. State Structure (FieldCollection Layout)

The simulation state is a `FieldCollection` — a flat list of `ScalarField` objects.
The layout depends on each component's **time derivative order** in the JSON specification:

| Time Order | Slots Used | Rate Computation |
|-----------|------------|-----------------|
| 2 (wave-like) | field + momentum | d/dt field = momentum; d/dt momentum = RHS |
| 1 (first-order) | field only | d/dt field = RHS (via virtual momentum) |
| 0 (constraint) | field only | d/dt field = 0 (solved algebraically) |

### Example: Chern-Simons (A_0 constraint, A_1 and A_2 first-order)

```mermaid
flowchart TB
    subgraph state["FieldCollection — state vector"]
        direction LR
        s0["Slot 0<br/><b>A_0</b><br/>(field)"]
        s1["Slot 1<br/><b>A_1</b><br/>(field)"]
        s2["Slot 2<br/><b>A_2</b><br/>(field)"]
    end

    subgraph layout["state_layout (from EquationSystem)"]
        direction TB
        l0["A_0: order=0 → constraint<br/>1 slot, rate = 0"]
        l1["A_1: order=1 → first-order<br/>1 slot, rate = RHS"]
        l2["A_2: order=1 → first-order<br/>1 slot, rate = RHS"]
    end

    subgraph rates["evolution_rate() output"]
        direction LR
        r0["d/dt A_0 = 0"]
        r1["d/dt A_1 = RHS₁"]
        r2["d/dt A_2 = RHS₂"]
    end

    layout -.-> state
    state --> rates
```

### Example: Coupled Scalars (two second-order fields)

```mermaid
flowchart TB
    subgraph state2["FieldCollection — state vector"]
        direction LR
        a0["Slot 0<br/><b>phi_0</b><br/>(field)"]
        a1["Slot 1<br/><b>pi_phi_0</b><br/>(momentum)"]
        a2["Slot 2<br/><b>chi_0</b><br/>(field)"]
        a3["Slot 3<br/><b>pi_chi_0</b><br/>(momentum)"]
    end

    subgraph rates2["evolution_rate() output"]
        direction LR
        b0["d/dt phi_0<br/>= pi_phi_0"]
        b1["d/dt pi_phi_0<br/>= Σ coeff × op(phi, chi)"]
        b2["d/dt chi_0<br/>= pi_chi_0"]
        b3["d/dt pi_chi_0<br/>= Σ coeff × op(phi, chi)"]
    end

    state2 --> rates2
```

Cross-field coupling: a term in the phi equation can reference `chi_0` via `{"operator": "laplacian", "field": "chi_0"}`.

---

## 4. Execution Flow: `evolution_rate()`

Called by the ODE solver at every timestep. Uses a two-pass algorithm to handle mixed time orders.

```mermaid
flowchart TB
    entry(["evolution_rate(state, t)"])
    entry --> bc["Resolve boundary conditions<br/>(cached after first call — B5)"]
    bc --> pass1

    subgraph pass1["Pass 1: Compute Virtual Momenta"]
        direction TB
        p1a{"Any 1st-order<br/>components?"}
        p1a -- yes --> p1b["Evaluate RHS →<br/>virtual_momenta#91;name#93;"]
        p1a -- no --> p1c
        p1b --> p1c{"Any constraint<br/>components?"}
        p1c -- yes --> p1d["Solve constraint<br/>or set zero field"]
        p1c -- no --> p1end[ ]
        p1d --> p1end
    end

    pass1 --> pass2

    subgraph pass2["Pass 2: Build Rate Vector"]
        direction TB
        p2a["<b>2nd-order</b><br/>rates#91;field_slot#93; = momentum<br/>rates#91;momentum_slot#93; = _compute_rhs()"]
        p2b["<b>1st-order</b><br/>rates#91;field_slot#93; = virtual_momenta#91;name#93;"]
        p2c["<b>Constraint</b><br/>rates#91;field_slot#93; = 0"]
        p2a --> p2b --> p2c
    end

    pass2 --> ret(["Return FieldCollection(rates)"])

    p2a -. "calls" .-> rhs

    subgraph rhs["_compute_rhs_for_component()"]
        direction TB
        r1["Pre-extract coord_arrays<br/>for position-dependent terms"]
        r2["Initialize per-timestep<br/>coefficient cache"]
        r3["For each term in equation spec:"]
        r4["1. Resolve target field from state"]
        r5["2. Apply operator → ScalarField"]
        r6["3. Resolve coefficient<br/>(pre-resolved → cache → eval)"]
        r7["4. result += coefficient × operated"]
        r1 --> r2 --> r3 --> r4 --> r5 --> r6 --> r7
    end
```

### Coefficient Resolution Hierarchy

Each RHS term has a coefficient that may be constant, time-dependent, or position-dependent.
Resolution follows a fast-path hierarchy:

```
1. Pre-resolved (B4)     — constant coefficients computed once at __init__
2. Per-timestep cache (C1) — deduplicates eval() for shared symbolic expressions
3. Full evaluation         — _resolve_coefficient_at_point(term, t, grid, coord_arrays)
     ├─ Constant: simple parameter lookup → float
     ├─ Time-only: eval(expr, {t=t, params}) → float
     └─ Position-dependent: eval(expr, {t, x, y, z}) → ndarray (grid-shaped)
```

---

## Working Examples

| Example | Dim | Fields | Key Features |
|---------|-----|--------|--------------|
| `scalar_field/` | 1+1D | phi_0 | Basic scalar wave, mass term |
| `electromagnetic/` | 1+1D | A_0, A_1 | Vector field, Lorenz gauge |
| `coupled_scalars/` | 1+1D | phi_0, chi_0 | Cross-field coupling, mass matrix |
| `chern_simons/` | 2+1D | A_0, A_1, A_2 | Epsilon tensor, constraint equation |
| `elasticity/` | 2+1D | u_0, u_1 | Anisotropic Laplacian, cross derivatives |
| `curved_spacetime/` | 2+1D | phi_0 | De Sitter Hubble friction, time-dependent coefficients |
| `sphere_kg/` | 2+1D | phi_0 | Stereographic S2, position-dependent coefficients |
| `polar_kg/` | 2+1D | phi_0 | Polar coordinates, 1/r Christoffel correction |
| `spherical_kg/` | 3+1D | phi_0 | Spherical coordinates, trig coefficients (Cot, Csc) |
| `cylindrical_kg/` | 3+1D | phi_0 | Mixed curved (r,theta) and flat (z) directions |
| `gravitational_waves/` | 3+1D | h_ij | xPert linearization, TT gauge, rank-2 tensor |

---

## Supported Operators

Operators map JSON terms to py-pde field operations. All support cross-field references.

| Operator | Min Dim | Description |
|----------|---------|-------------|
| `identity` | 1D | Field itself (mass terms) |
| `laplacian` | 1D | Full spatial Laplacian |
| `laplacian_x`, `_y`, `_z` | 1D/2D/3D | Directional second derivative |
| `gradient_x`, `_y`, `_z` | 1D/2D/3D | Directional first derivative |
| `cross_derivative_xy`, `_xz`, `_yz` | 2D/3D | Mixed partial derivative |
| `first_derivative_t` | 1D | Time derivative (friction terms) |
| `biharmonic` | 1D | Fourth-order Laplacian |
| `derivative_N_x` | 1D | Nth-order derivative along axis |
