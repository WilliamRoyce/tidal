# Design Document: xPert Linearization Support in `tidal derive`

## 1. Motivation and Context

### The Gap

TIDAL's `tidal derive` CLI generates Wolfram Language scripts (`.wls`) from TOML configuration files. The standard path uses Euler-Lagrange variation (`VarD`) to derive equations of motion from a Lagrangian. However, one important class of field theories -- **linearized gravity** -- cannot be expressed this way because the equations are obtained by **perturbation theory** rather than direct variation.

The gravitational waves example (`examples/gravitational_waves/`) was the only pipeline example that required a handwritten `.wls` script. Its derivation uses xAct's xPert package for systematic metric perturbation, a fundamentally different approach from VarD.

### The Insight

After thorough analysis of the handwritten `linearized_gravity.wls`, the key realization was that xPert linearization is a **clean, isolated 3-call extension** to the existing pipeline:

```
Standard VarD:    [lagrangian] -> VarD[field, CD][L] -> DecomposeToComponents -> JSON
xPert:            [linearization] -> SetupMetricPerturbation -> LinearizeTensorExpression -> convert notation -> DecomposeToComponents -> JSON
```

The difference is **only in the equation derivation step**. Everything before (spacetime, fields, constants) and after (component decomposition, JSON export) is identical. Once xPert produces the linearized tensor equation and we convert from `h[LI[1], -a, -b]` notation to plain `H[-a, -b]`, the existing pipeline handles everything unchanged.

This meant the feature could be implemented as a targeted extension with minimal risk to existing functionality.

## 2. Architecture Decisions

### 2.1 TOML Section Design: `[linearization]` vs `[lagrangian]`

**Decision:** Make `[linearization]` and `[lagrangian]` mutually exclusive sections.

**Rationale:**
- A theory is EITHER defined by a Lagrangian (VarD path) OR by linearizing a tensor expression (xPert path). These are fundamentally different derivation strategies.
- Mutual exclusivity prevents ambiguous configurations and simplifies validation.
- The TOML schema remains clean: config has exactly one of `[lagrangian]` or `[linearization]`.

**Alternative considered:** A single `[equations]` section with a `method = "euler-lagrange" | "linearization"` field. Rejected because the required parameters differ significantly (expression vs perturbation_field, etc.), making a unified section awkward.

### 2.2 Expression Substitution

**Decision:** Reuse `_substitute_field_names` for the linearization expression, with the same substitution rules as Lagrangian expressions.

**Rationale:**
- The user writes `Einstein[CD][-a, -b]` and the `CD` gets prefixed to `{prefix}CD`, just like in Lagrangian expressions.
- Field names, metric names (`eta`), and chart placeholders are all substituted uniformly.
- This means linearization expressions follow the exact same conventions as Lagrangian expressions -- no new syntax to learn.

**Key finding:** The `eta[` -> `{prefix}Eta[` substitution already existed in `_substitute_field_names`, added as part of the general symbol prefixing. This meant the massive gravity expression `Einstein[CD][-a, -b] + m2 eta[-a, -b]` would work with zero changes.

### 2.3 xPert Symbol Naming

**Decision:** Auto-generate internal xPert symbols as `{prefix}hpert` (perturbation) and `{prefix}Epsilon` (expansion parameter).

**Rationale:**
- xPert requires distinct symbols for the perturbation (`h[LI[1], ...]`) and the pipeline's plain tensor (`H[...]`).
- Using `{prefix}hpert` (lowercase) prevents collision with the pipeline field `{prefix}H` (uppercase).
- The expansion parameter `{prefix}Epsilon` is internal to xPert and never appears in the JSON output.
- No user-facing configuration needed -- these are implementation details.

### 2.4 Notation Conversion Bridge

**Decision:** Explicit pattern rule to convert xPert notation to pipeline notation.

The generated WLS includes:
```mathematica
linExprPlain = linExpr /. {prefix}hpert[LI[1], idx__] :> {prefix}H[idx];
```

**Rationale:**
- xPert produces expressions with `h[LI[1], -a, -b]` (order-indexed perturbation notation).
- The pipeline's `DecomposeToComponents` expects plain tensors `H[-a, -b]`.
- An explicit substitution rule is transparent, debuggable, and reliable.
- `Simplify[]` is applied after conversion to clean up the expression.

### 2.5 Constants with Linearization

**Decision:** Constants (`[constants].names`) work identically with both `[lagrangian]` and `[linearization]`.

**Rationale:**
- Constants are defined via `DefConstantSymbol` before any equation derivation.
- They flow through `LinearizeTensorExpression` as symbolic scalars (xPert treats them correctly).
- They appear in the JSON as symbolic coefficient strings, evaluable at runtime via `[parameters]`.
- No special handling needed -- the existing constant machinery is derivation-agnostic.

### 2.6 Package Loading

**Decision:** Conditionally load xPert and Linearize.wl only when `[linearization]` is present.

```python
def _wls_packages(pipeline_path: str, *, load_xpert: bool = False) -> list[str]:
```

**Rationale:**
- xPert is a large xAct package -- loading it when not needed would waste time and risk symbol conflicts.
- `Linearize.wl` (the pipeline's linearization module) is only relevant for xPert workflows.
- Standard VarD derivations continue to load exactly the same packages as before.

## 3. Implementation Summary

### Files Modified (for `[linearization]` feature)

**`tidal/cli/_derive.py`** -- All changes:
1. `_WlsContext` dataclass: added `linearization: dict[str, Any] | None` field
2. `_validate_linearization()`: new validation function for `[linearization]` section
3. `_validate_config()`: accepts either `[lagrangian]` or `[linearization]` (mutually exclusive)
4. `_wls_packages()`: `load_xpert: bool` parameter to conditionally load xPert + Linearize.wl
5. `_wls_linearization()`: new function generating xPert setup + linearize + convert + decompose
6. `generate_wls()`: dispatch to linearization or VarD path based on config
7. `_wls_metadata_and_export()`: `"linearized" -> True` metadata for linearization mode
8. `_substitute_field_names()`: `CD]` -> `{prefix}CD]` substitution (for `Einstein[CD]` pattern)

### Files Created (for massive gravity example)

- `examples/massive_gravity/theory.toml` -- TOML config with `[linearization]` + `[constants]` + `[parameters]`
- `examples/massive_gravity/massive_gravity.wls` -- Handwritten Wolfram derivation
- `examples/massive_gravity/simulation.py` -- Python simulation with dispersive propagation
- `examples/massive_gravity/run.sh` -- CLI equivalent commands

### Tests Added

- 5 validation tests: missing expression, missing perturbation_field, unknown field, mutual exclusivity, neither present
- 2 dry-run tests: basic linearization (GW-like), massive gravity (constants + metric reference)
- Total: 752 Python tests, 0 ruff violations, 0 pyright errors

## 4. Risk Analysis

### What Could Go Wrong

| Risk | Mitigation | Status |
|------|-----------|--------|
| Existing TOMLs break | `[linearization]` is opt-in; no changes to VarD path | Verified via backward-compat dry-runs |
| Expression substitution misses something | All substitution rules tested; `eta[`, `CD[`, `CD]` all verified | Dry-run shows correct output |
| xPert notation not fully converted | Explicit pattern rule catches all `LI[1]` occurrences | Proven in handwritten .wls |
| Constants not recognized by xPert | `DefConstantSymbol` called before linearization | Verified in massive_gravity dry-run |
| Python JSON loader breaks | JSON structure is identical to VarD output; no Python changes needed | N/A (zero Python pipeline changes) |

### Backward Compatibility

- All 16 existing TOMLs continue to work unchanged (verified via dry-runs of `scalar_field` and `gravitational_waves`)
- VarD code path has zero modifications
- JSON schema is unchanged
- Python `build_pde_from_json` needs zero changes

## 5. The Massive Gravity Example

### Why This Specific Example

The massive gravity example was chosen to exercise features that `gravitational_waves/` does not:

| Feature | gravitational_waves | massive_gravity |
|---------|--------------------|-----------------|
| Dimension | 3+1D (10 components) | 2+1D (6 components) |
| Constants | None | m2 (DefConstantSymbol) |
| Parameters | None | m2 = 1.0 (runtime) |
| Mass matrix | Zero | Non-zero diagonal |
| Metric in expression | No | Yes (eta[-a,-b]) |

### Physics Background

In 2+1D, pure GR has no local propagating degrees of freedom (the Weyl tensor vanishes identically). Adding a mass term `m^2 h_ab` creates a propagating massive mode with dispersion relation omega^2 = k^2 + m^2.

The expression `G_ab + m^2 g_ab = 0` is linearized by xPert as:
- `LinearizeTensorExpression[Einstein[CD][-a,-b]]` -> `G^(1)_ab[h]` (linearized Einstein)
- `LinearizeTensorExpression[m2 * eta[-a,-b]]` -> `m2 * h_ab` (first-order metric perturbation)

The background inconsistency (m^2 * eta != 0 at zeroth order) is irrelevant at the linear level -- xPert extracts the O(epsilon) part regardless. This is standard practice in massive gravity.

### Alternative Examples Considered

See `examples/massive_gravity/DESIGN_RATIONALE.md` for the full evaluation of 7 candidate examples and why massive gravity in 2+1D was selected.

## 6. Future Extensions

The `[linearization]` feature currently supports single-field metric perturbation. Natural extensions include:

1. **Multi-field linearization** -- Perturb metric AND matter fields simultaneously (e.g., Einstein-Maxwell)
2. **Non-flat backgrounds** -- Linearize around de Sitter, FLRW, or other curved backgrounds
3. **Higher-order perturbation** -- Second-order xPert expansion (gravitational self-interaction)
4. **Gauge automation** -- Automatic de Donder or Lorenz gauge fixing in the WLS
5. **Ricci/Weyl linearization** -- Alternative geometric tensors (already supported by `LinearizeTensorExpression` but untested via TOML)
