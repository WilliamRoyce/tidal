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
- Zwicker (2020), "py-pde: A Python package for solving partial differential equations", JOSS 5(48), 2158 — PDE solver backend

See [`docs/references.md`](references.md) for the full citation list.
