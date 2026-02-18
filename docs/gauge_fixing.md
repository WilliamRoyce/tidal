# Gauge Fixing: Optional Per-Field Toolkit

Gauge fixing simplifies equation structure for theories with gauge symmetry
(massless vectors, linearized gravity). It is **never required** — TIDAL's
existing pipeline handles gauge-invariant theories correctly, and all
measurement quantities (energy, conversion, mixing) are gauge-invariant.

## Quick Start

Add a `[[gauge]]` section to your `theory.toml`:

```toml
[[fields]]
name = "A"
type = "vector"

[[gauge]]
field = "A"
type = "lorenz"
xi = 1.0              # optional, default 1.0 (Feynman gauge)
```

This adds the Lorenz gauge-fixing term `-(1/2xi)(d_mu A^mu)^2` to the
Lagrangian before Euler-Lagrange derivation. For Maxwell equations, this
reduces coupled equations to uncoupled wave equations for each component.

No `[[gauge]]` section means no gauge fixing (the default). Existing
examples continue to work unchanged.

## Named Presets

### Type A — Lagrangian term (added to L before EL derivation)

| Preset | Fields | Expression | Effect |
|--------|--------|------------|--------|
| `lorenz` | vector | `-(1/2xi)(d_mu A^mu)^2` | Maxwell -> uncoupled wave equations |
| `de_donder` | sym. rank-2 | `-(1/2xi)(d_a h^a_b - 1/2 d_b h)^2` | Lin. Einstein -> uncoupled waves |

### Type B — Constraint (applied after EOM derivation)

| Preset | Fields | Constraint | Effect |
|--------|--------|------------|--------|
| `temporal` | vector | `A_0 = 0` | Zeroes temporal component in all equations |
| `coulomb` | vector | `div A_spatial = 0` | Adds transversality constraint equation |
| `axial` | vector | `A_n = 0` | Zeroes last spatial component in all equations |

Type A presets use Wolfram builder functions in `GaugeFix.wl`.
Type B presets operate at the component level after `DecomposeToComponents`,
modifying or augmenting the equation system directly.

## Custom Gauge Terms

For full flexibility, use `type = "custom"` with an arbitrary Wolfram
expression and `mechanism = "lagrangian_term"`:

```toml
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
expression = "-(1/(2*xi)) * eta[a,b] CD[-a][A[-b]] eta[c,d] CD[-c][A[-d]]"
```

The `expression` field uses the same syntax as `[lagrangian].expression` and
`[[derived_fields]].definition` — field names, `eta`, `CD`, and chart
placeholders are all substituted automatically.

> **Note:** Custom constraint gauges (`mechanism = "constraint"`) are not yet
> supported. Use a named preset (`temporal`, `coulomb`, `axial`) for
> constraint-based gauge fixing.

### Worked Example: Lorenz Gauge Two Ways

**Using the preset:**
```toml
[[gauge]]
field = "A"
type = "lorenz"
xi = 1.0
```

**Equivalent custom expression:**
```toml
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
expression = "-(1/2) * eta[a,b] CD[-a][A[-b]] eta[c,d] CD[-c][A[-d]]"
```

Both produce the same gauge-fixed Lagrangian. The preset is more readable;
the custom form shows the underlying mechanism and serves as a template for
non-standard gauge choices.

## Two Mechanisms

**Type A — Lagrangian term** (`mechanism = "lagrangian_term"`):
- An expression is added to L *before* Euler-Lagrange derivation
- Changes the structure of the equations of motion
- The gauge-fixed Lagrangian is canonicalized (`ToCanonical` + `ContractMetric`)
  before EL derivation
- Example: Lorenz gauge adds `-(1/2xi)(div A)^2`, turning Maxwell into wave equations
- Presets: `lorenz`, `de_donder`

**Type B — Constraint** (`mechanism = "constraint"`):
- Applied *after* EOM derivation and component decomposition
- For `temporal` and `axial`: substitutes the constrained component (and all its
  derivatives) with zero via `ReplaceAll` on the component equations
- For `coulomb`: appends a spatial divergence constraint equation
  (`div A_spatial = 0`) to the equation system as a `time_order=0` equation
- The original EOM structure is preserved; constraints modify or augment it
- Presets: `temporal`, `coulomb`, `axial`

For named presets, the mechanism is inferred automatically (e.g., `lorenz`
is always `lagrangian_term`, `temporal` is always `constraint`).

## Multi-Field Example

Different fields can have different gauge choices, or no gauge at all:

```toml
[[fields]]
name = "A"
type = "vector"

[[fields]]
name = "B"
type = "vector"

# Lorenz gauge on A only
[[gauge]]
field = "A"
type = "lorenz"

# B has no gauge fixing (default)
```

## Adding New Presets (Developer Guide)

The process depends on the gauge mechanism.

### Type A (Lagrangian term) — 3 steps

**Step 1:** Write a builder function in `tidal/wolfram/GaugeFix.wl`:

```mathematica
BuildMyGaugeTerm[field_, metric_, covd_, xi_:1] := Module[
  {gaugeTerm},

  If[!NumericQ[xi] || xi <= 0,
    Throw["BuildMyGaugeTerm: xi must be positive, got " <> ToString[xi]]
  ];

  (* Your gauge-fixing expression here *)
  gaugeTerm = (* ... *);

  gaugeTerm
];
```

The function must return an *expression* (a Wolfram symbolic expression),
not perform an action. The pipeline calls `AddGaugeFixingTerm[L, yourTerm]`
to inject it into the Lagrangian.

**Step 2:** Add a registry entry in `tidal/cli/_derive.py`:

```python
_GAUGE_PRESETS["my_gauge"] = {
    "mechanism": "lagrangian_term",
    "builder": "BuildMyGaugeTerm",    # Wolfram function name
    "requires": "vector",             # or "tensor"
}
```

**Step 3:** (Optional) Add validation rules in `_validate_gauge()`.

### Type B (Constraint) — 2 steps

Type B presets operate at the component level in Python — no Wolfram builder
function needed.

**Step 1:** Add a registry entry (no `builder` key):

```python
_GAUGE_PRESETS["my_constraint"] = {
    "mechanism": "constraint",
    "requires": "vector",
}
```

**Step 2:** Add a handler in `_wls_gauge_fixing_type_b()` in `_derive.py`.
Use `_type_b_zero_component()` for substitution-type constraints or
`_type_b_coulomb_constraint()` as a template for differential constraints.

## FAQ

**Do I need gauge fixing?**
No. All existing TIDAL examples work without it. Gauge fixing is a
computational convenience — it simplifies equation structure but doesn't
change the physics. Massless vectors already get implicit constraint
structure (A_0 becomes time_order=0). Massive vectors have no gauge
freedom at all.

**Are measurements affected by gauge choice?**
No. Energy density, conversion probability, mixing length, and spectral
quantities are all gauge-invariant physical observables.

**Can I mix presets and custom gauges?**
Yes. Each `[[gauge]]` entry is independent. You can use `type = "lorenz"`
on one field and `type = "custom"` on another.

**Can I mix Type A and Type B gauges?**
Yes. For example, `lorenz` on field A (Type A) and `temporal` on field B
(Type B) is valid. Type A modifies the Lagrangian before EL; Type B
modifies the component equations after decomposition.

**What happens if I gauge-fix a field that has no gauge symmetry?**
For massive fields (Proca), the mass term already breaks gauge symmetry.
Adding a gauge-fixing term is mathematically valid but physically
unnecessary and may produce unexpected equations. Validation prevents
gauge-fixing scalar fields (which never have gauge symmetry).
