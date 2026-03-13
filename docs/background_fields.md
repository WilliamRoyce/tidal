# Background Fields: Position-Dependent Coupling

Background fields are non-dynamical tensors that appear in the Lagrangian
but are NOT varied in the Euler-Lagrange derivation. They survive as
(possibly position-dependent) coefficients in the equations of motion.

This feature is essential for simulating the Gertsenshtein effect
(Gertsenshtein 1962; Domcke & Garcia-Cely 2023), where a static external
magnetic field B₀(x) catalyses photon-graviton conversion. It also enables
general probe-field approximations and theories with fixed external sources.

## Pipeline Trace

### 1. TOML Declaration

```toml
[[background_fields]]
name = "G"
type = "scalar"
components = ["g0 * Exp[-(x[]^2 + y[]^2) / (2 * R^2)]"]
```

Component values can be:

- **Numbers**: `0`, `1.0` (constant)
- **Strings**: Wolfram expressions evaluated symbolically (may depend on coordinates)

#### ⚠ Index Convention: `components` always specifies **covariant** (lowered-index) values

For a vector background field `A`, `components[μ]` defines **`A_μ`** (covariant, lower index).

```
components = [A_0, A_1, ..., A_{d-1}]   ← COVARIANT (A with subscript)
```

This is the convention throughout TIDAL's Wolfram pipeline:
- `SetBackgroundFieldDownValues` sets `field[{μ, -chart}] = components[μ]` (covariant `-chart`)
- `EvaluatePDBackgroundField` processes covariant `{μ, -chart}` forms
- xAct's `ComponentValue` stores covariant components

**Contravariant forms are derived automatically:**
- For **Minkowski** metrics (`metric = "minkowski"`): `A^μ ≈ A_μ` (spatial components identical; temporal `A^0 = -A_0`, but typically zero)
- For **curved diagonal** metrics (`metric = "diagonal"`): `A^μ = g^{μμ} A_μ = A_μ / g_{μμ}`, computed via `Simplify[(A_μ) / (g_{μμ})]`

**Example: Dipolar magnetic field in spherical coordinates** (`diag[-1, 1, r², r²sin²θ]`):
```
A_θ = Bpeak*z0³/(2r²)    ← what you put in components[2]
A^θ = g^{θθ} A_θ = A_θ/r² = Bpeak*z0³/(2r⁴)   ← computed automatically
```

**Gauge potentials**: Use the covariant gauge potential `A_μ` in `components`.
The Faraday tensor `F_{μν} = ∂_μ A_ν - ∂_ν A_μ` is metric-independent (exterior derivative),
so `F` is determined correctly from covariant `A_μ` alone.

### 2. WLS Generation (`_derive.py`)

- `_wls_fields()` emits `DefTensor[G[], M2]` + `ComponentValue[G[{0,-Cart}], expr]`
- **Scalar backgrounds**: Explicit `ReplaceAll` via `_wls_scalar_background_substitution()` before decomposition (ToBasis doesn't trigger ComponentValue for rank-0 tensors)
- **Vector/tensor backgrounds**: Handled by `ComponentValue` + `ToBasis` in `DecomposeToComponents`

### 3. Wolfram Pipeline

- `VarD` is called **only** on dynamical fields — background fields are correctly NOT varied
- The background field's symbolic expression survives as a coefficient in the EOM

### 4. JSON Export (`ExportJSON.wl`)

- `IsCoordinateDependentCoefficient` detects xCoba symbols in coefficients
- Sets `coordinate_dependent: ["x", "y"]` and `coefficient_symbolic: "g0*Exp[...]"`
- `time_dependent: true` if the coefficient contains `t[]`

### 5. Python PDE Evaluation (`pde_builder.py`)

- `_mathematica_to_python()` converts Mathematica InputForm to Python (~30 functions)
- `_resolve_coefficient_at_point()` evaluates on numpy grid arrays via `eval()` with restricted namespace (`__builtins__: {}`)

### 6. Measurement Module (`_energy.py`)

- `evaluate_coefficient()` in `_eval_utils.py` provides standalone evaluation
- `_resolve_mass_squared()` returns ndarray for position-dependent mass terms
- `_compute_virial_potential()` evaluates position-dependent coefficients on the grid

## Caching Architecture

Four levels of caching prevent redundant evaluation:

| Level | Cache            | What                                             | Lifetime                                   |
| ----- | ---------------- | ------------------------------------------------ | ------------------------------------------ |
| L0    | `_preresolved`   | Fully constant coefficients (no `eval()` needed) | Instance (construction)                    |
| L1    | `_expr_cache`    | Mathematica→Python string conversion             | Instance (singleton)                       |
| L2    | `_spatial_cache` | Grid arrays for spatial-only coefficients        | Instance (populated on first RHS call)     |
| L3    | `coeff_cache`    | Per-RHS-call results                             | Single `_compute_rhs_for_component()` call |

**Key behavior**: Spatial-only coefficients (`position_dependent=True`, `time_dependent=False`) are evaluated once on the grid and cached persistently across all timesteps. Time-dependent coefficients are re-evaluated each substep. Fully constant coefficients (neither position- nor time-dependent) are resolved once at construction time (L0) with no `eval()` at runtime.

## Supported Patterns

| Pattern                   | Example                                | Status                           |
| ------------------------- | -------------------------------------- | -------------------------------- |
| Constant background       | `components = ["B0"]`                  | Fully supported                  |
| Position-dependent scalar | `components = ["g0 * Exp[-r^2]"]`      | Fully supported                  |
| Position-dependent vector | `components = [0, 0, "B0 * Sin[x[]]"]` | Fully supported                  |
| Time-dependent background | `components = ["g0 * Cos[t[]]"]`       | Works (Python-level tested)      |
| Background field gradient | `CD[-a][G[]]` in Lagrangian            | **NOT supported** — raises error |

## Background Field Validity: Externally Imposed Fields

**TIDAL does NOT validate that background fields satisfy the background equations of
motion.** Validation is purely structural (type, rank, component count, name collisions
— see `_derive_validate.py`). The user is responsible for ensuring physical consistency.

For linearized perturbation theory, background fields are treated as **externally imposed**
— maintained by some external source (e.g., a magnet, stellar interior currents, or an
applied potential). The perturbation equations are valid regardless of whether the
background configuration is self-consistent, because:

1. The background is NOT varied in the Euler-Lagrange derivation (`VarD` acts only on
   dynamical perturbation fields)
2. The perturbation equations depend on the background only through its values and
   derivatives at each grid point
3. Whether those values satisfy the background EOM is irrelevant to the linearized dynamics

**Example**: A spatially-varying magnetic field B₀(z) implies ∇×B = J ≠ 0, violating
vacuum Maxwell equations. This is not a problem — B₀(z) is interpreted as the field
produced by external currents (e.g., a solenoid or stellar magnetosphere). This treatment
is standard in the graviton-photon mixing literature:
- Raffelt & Stodolsky (1988, PRD 37:1237) — treat B(z) as external
- Boccaletti et al. (1970, Nuovo Cimento 70B:129) — arbitrary B(z) profiles
- Domcke & Garcia-Cely (2024, JCAP 05:051) — inhomogeneous B(z)

**When self-consistency matters**: For strong-field or nonlinear regimes (not supported
by TIDAL's linearized framework), the background must satisfy its own EOM. In TIDAL's
linearized regime, this is never required.

## Energy Measurement Limitation with Position-Dependent Coefficients

**Known issue**: Energy measurements (`tidal measure --what energy`) do not properly
account for spatially-varying mass or potential terms. The Hamiltonian density computation
assumes translation-invariant coefficients. This affects:

- Position-dependent mass terms (e.g., `m²(z) φ²`)
- Localized coupling regions (e.g., Gaussian B₀(z) in Gertsenshtein scattering)
- Potential wells (`examples/scalar_potential_well/`)

**Workaround**: Use conversion probability measurement (`--what conversion`) as the
primary validation metric — it measures amplitude ratios, which are unaffected by this
issue. Energy conservation checks should use the `--what conservation` measurement,
which tracks relative change rather than absolute values.

## Known Limitations

1. **No background field gradients**: `CD[-a][G[]]` in the Lagrangian will raise
   `ValueError` at TOML validation time. The scalar `ReplaceAll` substitution
   doesn't handle derivatives of background fields.

2. **Time-dependent backgrounds less tested**: The Python evaluation pipeline
   handles time-dependent expressions correctly, but the Wolfram-side behavior
   (xAct handling time derivatives of time-dependent background expressions)
   is not end-to-end tested with wolframscript.

3. **Sign change diagnostic**: If a position-dependent mass term changes sign
   across the grid, a `UserWarning` is emitted (possible tachyonic instability).
   This is informational, not a hard error.

## Adding New Mathematica Functions

To support a new Mathematica function (e.g., `Gamma[x]`), update
`_eval_utils.py` only — `pde_builder.py` delegates to it automatically:

1. **String conversion**: Add to `_FUNCTION_MAP` in `_eval_utils.py`:

   ```python
   ("Gamma", "gamma"),
   ```

2. **Evaluation namespace**: Add to `build_eval_namespace()` in `_eval_utils.py`:

   ```python
   ns["gamma"] = special.gamma
   ```

3. **Test**: Add a conversion test in `test_background_fields.py`.

`PDEFromSpec._mathematica_to_python()` and `PDEFromSpec._build_base_namespace()`
both delegate to `_eval_utils.py`, so no changes are needed in `pde_builder.py`.

## References

- Gertsenshtein (1962), "Wave resonance of light and gravitational waves", JETP 14, 84 — original prediction requiring external B-field as background
- Domcke & Garcia-Cely (2023), "A simple derivation of the Gertsenshtein effect", [arXiv:2301.02072](https://arxiv.org/abs/2301.02072) — modern derivation with inhomogeneous background profiles
- Martín-García et al., "xAct: Efficient tensor computer algebra for Mathematica", [xact.es](http://www.xact.es/) — symbolic tensor algebra (`VarD`, `ComponentValue`, `ToBasis`)
- Zwicker (2020), "py-pde: A Python package for solving partial differential equations", JOSS 5(48), 2158 — original PDE backend; FD stencil conventions retained in TIDAL's native operators

See [`docs/references.md`](references.md) for the full citation list.
